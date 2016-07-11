/*Copyright (C) 2015 Olivier Delaneau, Halit Ongen, Emmanouil T. Dermitzakis
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

#include "ase_data.h"

typedef struct {     			// auxiliary data structure
	samFile * fp;				// the file handle
	bam_hdr_t * hdr;			// the file header
	hts_itr_t * iter;			// NULL if a region not specified
} aux_t;

static int read_bam(void *data, bam1_t *b) {
	aux_t * aux = (aux_t*) data;
	return (aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b));
}

void ase_data::readSequences(string fbam, string fvcf) {
	aux_t * data = (aux_t *) malloc (sizeof(aux_t));

	vrb.title("Opening BAM file [" + fbam + "] and writing VCF file [" + fvcf + "]");

	//BAM OPEN
	data->fp = sam_open(fbam.c_str(), "r");
    if (data->fp == 0) vrb.error("Cannot open BAM file!");
    data->hdr = sam_hdr_read(data->fp);
    if (data->hdr == 0) vrb.error("Cannot parse BAM header!");
    hts_idx_t *idx = sam_index_load(data->fp, fbam.c_str());
    if (idx == NULL) vrb.error("Cannot load BAM index!");

    //VCF OPEN
	htsFile * bcf_fd = NULL;
	bool compressed_vcf = fvcf.substr(fvcf.find_last_of(".") + 1) == "gz";
	if (compressed_vcf) bcf_fd = bcf_open(fvcf.c_str(), "wz");
	else bcf_fd = bcf_open(fvcf.c_str(), "wu");
	if (bcf_fd == NULL) vrb.error("Impossible to create VCF file");
	else if (compressed_vcf) vrb.bullet("BGZIP for VCF compression is ON");
	else vrb.bullet("BGZIP compression for VCF is OFF (add .gz to filename to activate it)");
	bcf_hdr_t * bcf_hdr = bcf_hdr_init("w");
	kstring_t str = {0,0,0};
	ksprintf(&str, "##QTLtools ase Version=%s\n",QTLTOOLS_VERSION);
	bcf_hdr_append(bcf_hdr, str.s);
    free(str.s);

	//VCF INFO
	vrb.bullet("Writing VCF header [INFO, CONTIG, FORMAT, SAMPLES]");
	bcf_hdr_append(bcf_hdr,"##INFO=<ID=M_FAI,Number=1,Type=Integer,Description=\"Number of reads failing mapping filters\">");
	if (!param_dup_rd) bcf_hdr_append(bcf_hdr,"##INFO=<ID=M_DUP,Number=1,Type=Integer,Description=\"Number of duplicate reads removed\">");
	bcf_hdr_append(bcf_hdr,"##INFO=<ID=M_SUC,Number=1,Type=Integer,Description=\"Number of reads passing mapping filters\">");
	bcf_hdr_append(bcf_hdr,"##INFO=<ID=B_DEL,Number=1,Type=Integer,Description=\"Number of reads with deletion\">");
	bcf_hdr_append(bcf_hdr,"##INFO=<ID=B_DIS,Number=1,Type=Integer,Description=\"Number of reads with base discordance\">");
	bcf_hdr_append(bcf_hdr,"##INFO=<ID=B_QUA,Number=1,Type=Integer,Description=\"Number of reads with low base quality\">");
	bcf_hdr_append(bcf_hdr,"##INFO=<ID=B_REF,Number=1,Type=Integer,Description=\"Number of reads matching REF allele\">");
	bcf_hdr_append(bcf_hdr,"##INFO=<ID=B_ALT,Number=1,Type=Integer,Description=\"Number of reads matching ALT allele\">");

	//VCF CONTIG
	set < string > contig_set (regions.begin(), regions.end());
	vrb.bullet("Writing " + stb.str(contig_set.size()) + " CONTIG fields");
	for (set < string > :: iterator it_c = contig_set.begin() ; it_c != contig_set.end() ; ++ it_c) {
		string tmp_str = "##contig=<ID=" + *it_c + ">";
		bcf_hdr_append(bcf_hdr, tmp_str.c_str());
	}

    //VCF FORMAT
    bcf_hdr_append(bcf_hdr,"##FORMAT=<ID=AS,Number=1,Type=Float,Description=\"Binomial test for ASE (i.e. ref/(ref+alt) != 0.5)\">");

    //VCF SAMPLES
    bcf_hdr_add_sample(bcf_hdr, sample_id[0].c_str());
    bcf_hdr_add_sample(bcf_hdr, NULL);

    //VCF HEADER
    bcf_hdr_write(bcf_fd, bcf_hdr);

	//Loop across regions
    bcf1_t * bcf_rec = bcf_init1();
	for (int reg = 0; reg < regions.size() ; reg ++) {

		//Jump to region
		data->iter = sam_itr_querys(idx, data->hdr, regions[reg].c_str()); // set the iterator
		if (data->iter == NULL) vrb.error("Problem jumping to region [" + regions[reg] + "]");
		else vrb.bullet("Scanning region [" + regions[reg] + "]");

		int beg = data->iter->beg;
		int end = data->iter->end;

		//Pile up reads
		const bam_pileup1_t * v_plp;
		int n_plp = 0, tid, pos, i_site = 0;
		bam_plp_t s_plp = bam_plp_init(read_bam, (void*)data);
		while (((v_plp = bam_plp_auto(s_plp, &tid, &pos, &n_plp)) != 0) && i_site < variants[reg].size()) {
		    int chr = bam_name2id(data->hdr, variants[reg][i_site].chr.c_str());
			if (pos < beg || pos >= end) continue;
			while (i_site < variants[reg].size() && (chr != tid || pos > variants[reg][i_site].pos)) { i_site ++; }
			if (tid == chr && pos == variants[reg][i_site].pos) {
				int m_fai = 0, m_dup = 0, m_suc = 0;
				int b_del = 0, b_ref = 0, b_alt = 0, b_dis = 0, b_qua = 0;

				//STEP1: Parse sequencing reads
				for (int iread = 0 ; iread < n_plp ; iread ++) {
					bool failed_qc = false;
					const bam_pileup1_t * p = v_plp + iread;

					if (p->b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL)) failed_qc = true;
					else if ((int)p->b->core.qual < param_min_mapQ) failed_qc = true;
					else if (p->b->core.flag & BAM_FPAIRED) {
						if (!(p->b->core.flag & BAM_FPROPER_PAIR)) failed_qc = true;
						else if (p->b->core.flag & BAM_FMUNMAP) failed_qc = true;
						else if ((p->b->core.flag & BAM_FREVERSE) == (p->b->core.flag & BAM_FMREVERSE)) failed_qc = true;
					}

					if (failed_qc) m_fai ++;
					else {
						if (p->b->core.flag & BAM_FDUP) m_dup ++;
						if (!param_dup_rd || !(p->b->core.flag & BAM_FDUP)) {
							if (p->is_del || p->is_refskip || p->indel == 1) b_del++;
							else if (bam_get_qual(p->b)[p->qpos] < param_min_baseQ) b_qua++;
							else {
								char base = ase_getBase(bam_seqi(bam_get_seq(p->b), p->qpos));
								bool isRef = (base == variants[reg][i_site].ref);
								bool isAlt = (base == variants[reg][i_site].alt);
								if (isRef) b_ref++;
								if (isAlt) b_alt++;
								if (!isRef && !isAlt) b_dis++;
							}
							m_suc ++;
						}
					}
				}

				//STEP2: Write VCF record
				if ((b_ref + b_alt) >= param_min_cov) {
					bcf_clear1(bcf_rec);
					bcf_rec->rid = bcf_hdr_name2id(bcf_hdr, variants[reg][i_site].chr.c_str());
					bcf_rec->pos = variants[reg][i_site].pos;
					bcf_update_id(bcf_hdr, bcf_rec, variants[reg][i_site].sid.c_str());
					string str_alleles = "N,N";
					str_alleles[0] = variants[reg][i_site].ref;
					str_alleles[2] = variants[reg][i_site].alt;
					bcf_update_alleles_str(bcf_hdr, bcf_rec, str_alleles.c_str());
					bcf_rec->qual = 100;
					int tmpi = bcf_hdr_id2int(bcf_hdr, BCF_DT_ID, "PASS");
					bcf_update_filter(bcf_hdr, bcf_rec, &tmpi, 1);
					bcf_update_info_int32(bcf_hdr, bcf_rec, "M_FAI", &m_fai, 1);
					if (!param_dup_rd) bcf_update_info_int32(bcf_hdr, bcf_rec, "M_DUP", &m_dup, 1);
					bcf_update_info_int32(bcf_hdr, bcf_rec, "M_SUC", &m_suc, 1);
					bcf_update_info_int32(bcf_hdr, bcf_rec, "B_DEL", &b_del, 1);
					bcf_update_info_int32(bcf_hdr, bcf_rec, "B_DIS", &b_dis, 1);
					bcf_update_info_int32(bcf_hdr, bcf_rec, "B_QUA", &b_qua, 1);
					bcf_update_info_int32(bcf_hdr, bcf_rec, "B_REF", &b_ref, 1);
					bcf_update_info_int32(bcf_hdr, bcf_rec, "B_ALT", &b_alt, 1);
					float tmpf = ase_binomialTest(b_ref, b_ref+b_alt, 0.5);
					bcf_update_format_float(bcf_hdr, bcf_rec, "AS", &tmpf, 1);
					bcf_write1(bcf_fd, bcf_hdr, bcf_rec);
				}
			}
		}
		bam_plp_reset(s_plp);
		bam_plp_destroy(s_plp);
    }

    bam_hdr_destroy(data->hdr);
    hts_idx_destroy(idx);
    if (data->fp) sam_close(data->fp);
    hts_itr_destroy(data->iter);
    free(data);
    bcf_destroy1(bcf_rec);
    hts_close(bcf_fd);
    bcf_hdr_destroy(bcf_hdr);
}
