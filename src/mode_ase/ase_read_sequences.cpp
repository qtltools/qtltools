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
	int mq;
} aux_t;

static int ase_read_bam(void *data, bam1_t *b) {
	aux_t * aux = (aux_t*) data;
	return (aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b));
}

static int mplp_func(void *data, bam1_t *b){
	aux_t * aux = (aux_t*) data;
	int ret, skip = 0;
	do{
		ret =  (aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b));
		if (ret < 0) break;
		if (b->core.tid < 0 || (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP))) {/*cerr << "FAIL" << endl*/; skip = 1; continue;}
		skip = 0;
		if (b->core.qual < aux->mq) {/*cerr << "MAPQ" << endl*/; skip = 1;}
		else if ((b->core.flag&BAM_FPAIRED) && !(b->core.flag&BAM_FPROPER_PAIR)) {/*cerr << "PAIR" << endl*/; skip = 1;}

	}while(skip);
	return ret;
}

static inline void pileup_seq(const bam_pileup1_t * p){
	if(p->is_head){
		cerr << "^";
		cerr << (char) (p->b->core.qual > 93? 126 : p->b->core.qual + 33);
	}
	if(!p->is_del){
		char c = p->qpos < p->b->core.l_qseq ? seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos)] : 'N';
		if (c == '=') c = bam_is_rev(p->b)? ',' : '.';
        else c = bam_is_rev(p->b)? tolower(c) : toupper(c);
		cerr << c;
	}else cerr << (char) (p->is_refskip ? (bam_is_rev(p->b)? '<' : '>') : '*');
	if (p->indel > 0) {
        cerr << '+' << p->indel; 
        for (int j = 1; j <= p->indel; ++j) {
            int c = seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos + j)];
            cerr << (char) (bam_is_rev(p->b)? tolower(c) : toupper(c));
        }
    } else if (p->indel < 0) {
        cerr << p->indel;
        for (int j = 1; j <= -p->indel; ++j) {
            int c = 'N';
            cerr << (char) (bam_is_rev(p->b)? tolower(c) : toupper(c));
        }
    }
    if (p->is_tail) cerr << '$';
}

void ase_data::parseBamMpileup(void * d){
	int n = 1;
	aux_t **data;
	data = (aux_t**) calloc(n, sizeof(aux_t*));
	data[0] = (aux_t *) d;
	const bam_pileup1_t **plp;
	plp = (const bam_pileup1_t **)calloc(n, sizeof(bam_pileup1_t*));
	int tid, pos, *n_plp;
	n_plp = (int*) calloc(n, sizeof(int));
	bam_mplp_t iter;
	iter = bam_mplp_init(n, ase_read_bam, (void**)data);
	bam_mplp_set_maxcnt(iter, max_depth);

	int beg,end,ret;
	if (data[0]->iter == NULL){
		beg = INT_MIN;
		end = INT_MAX;
	}else{
		beg = data[0]->iter->beg;
		end = data[0]->iter->end;
	}
	unsigned long linecount = 0;
	while ( (ret=bam_mplp_auto(iter, &tid, &pos, n_plp, plp)) > 0) {
		if (pos < beg || pos > end) continue;
		linecount++;
		if(linecount % 50000000 == 0) vrb.bullet("Processed " + stb.str(linecount) + " positions");
		string chr = data[0]->hdr->target_name[tid];
		ase_site temp(chr,pos);
		auto av_it = all_variants.find(temp);
		if (av_it != all_variants.end()){
			//unsigned int m_fai = 0, m_dup = 0, m_suc = 0, m_skip = 0;
			//unsigned int b_del = 0, b_ref = 0, b_alt = 0, b_dis = 0, b_qua = 0;
			//unsigned int m_sec = 0, m_umap = 0, m_fqc = 0, m_npp = 0, m_mapq = 0,m_mum = 0 , m_mio = 0;
			unsigned int b_ref = 0, b_alt = 0, b_dis = 0;
			mapping_stats ms;
			//STEP1: Parse sequencing reads
			if (n_plp[0] >= max_depth) vrb.warning(av_it->sid + " depth " + stb.str(n_plp[0]) + " is >= max-depth potential data loss!");
			for (int iread = 0 ; iread < n_plp[0] ; iread ++) {
				bool failed_qc = false;
				const bam_pileup1_t * p = plp[0] + iread;
				
				if (p->b->core.tid < 0 || (p->b->core.flag & BAM_FUNMAP)) {ms.unmapped++; failed_qc = true;}
				else if (p->b->core.flag & BAM_FSECONDARY) {ms.secondary++; failed_qc = true;}
				else if (p->b->core.flag & BAM_FQCFAIL) {ms.fail_qc++; if(!keep_failqc) failed_qc = true;}
				else if ((int)p->b->core.qual < param_min_mapQ) {ms.fail_mapping++; failed_qc = true;}
				else if (p->b->core.flag & BAM_FPAIRED) {
					if (p->b->core.flag & BAM_FMUNMAP) {ms.mate_unmapped++; if(!keep_orphan) failed_qc = true;}
					else if (!(p->b->core.flag & BAM_FPROPER_PAIR)) {ms.not_pp++; if(check_proper_pair) failed_qc = true;}
					//else if ((p->b->core.flag & BAM_FREVERSE) == (p->b->core.flag & BAM_FMREVERSE)) {ms.orientation++; failed_qc = true;}
				}



				/*if (p->b->core.tid < 0 || (p->b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL))) failed_qc = true;
				else if ((int)p->b->core.qual < param_min_mapQ) failed_qc = true;
				else if (p->b->core.flag & BAM_FPAIRED) {
					if (!(p->b->core.flag & BAM_FPROPER_PAIR)) failed_qc = true;
					//else if (p->b->core.flag & BAM_FMUNMAP) failed_qc = true;
					//else if ((p->b->core.flag & BAM_FREVERSE) == (p->b->core.flag & BAM_FMREVERSE)) failed_qc = true;
				}*/

				if (!failed_qc){
					if (p->b->core.flag & BAM_FDUP){
						ms.duplicate++;
						if(param_dup_rd) continue;
					}
					if (p->indel != 0 && !(p->is_del || p->is_refskip)){
						ms.indel++;
						if(param_rm_indel) continue;
					}
					if (p->is_del || p->is_refskip) {
						ms.skipped++;
						continue;
					}
					if ((p->qpos < p->b->core.l_qseq ? bam_get_qual(p->b)[p->qpos] : 0) < param_min_baseQ) {
						ms.fail_baseq++;
						continue;
					}
					char base = ase_getBase(bam_seqi(bam_get_seq(p->b), p->qpos));
					bool isRef = (base == av_it->ref);
					bool isAlt = (base == av_it->alt);
					if (isRef) b_ref++;
					if (isAlt) b_alt++;
					if (!isRef && !isAlt) b_dis++;
				}


				/*if (!failed_qc){
					if (p->b->core.flag & BAM_FDUP) ms.duplicate ++;
					else if (p->indel != 0 && !(p->is_del || p->is_refskip)) ms.indel++;
					if (!param_dup_rd || !(p->b->core.flag & BAM_FDUP)) {
						if (p->is_del || p->is_refskip) ms.skipped++;
						else if (!param_rm_indel || p->indel == 0){
							if ((p->qpos < p->b->core.l_qseq ? bam_get_qual(p->b)[p->qpos] : 0) < param_min_baseQ) ms.fail_baseq++;
							else {
								char base = ase_getBase(bam_seqi(bam_get_seq(p->b), p->qpos));
								bool isRef = (base == av_it->ref);
								bool isAlt = (base == av_it->alt);
								if (isRef) b_ref++;
								if (isAlt) b_alt++;
								if (!isRef && !isAlt) b_dis++;
							}
						}
					}
				}*/
			}
			//cerr << m_fai << " " << m_dup << " " << b_del << " " << b_qua << " " << m_suc << endl;
			ase_site current = *av_it;
			current.setCounts(b_ref,b_alt,b_dis, ms);
			passing_variants.push_back(current);
		}
	}
	bam_mplp_destroy(iter);
}

void ase_data::parseBam(void * d){
	aux_t * data = (aux_t *) d;
	//Pile up reads
	const bam_pileup1_t * v_plp;
	int n_plp = 0, tid, pos;
	bam_plp_t s_plp = bam_plp_init(ase_read_bam, (void*)data);
	bam_plp_set_maxcnt(s_plp, max_depth);

	int beg,end;
	if (data->iter == NULL){
		beg = INT_MIN;
		end = INT_MAX;
	}else{
		beg = data->iter->beg;
		end = data->iter->end;
	}

	unsigned long linecount = 0;
	while((v_plp = bam_plp_auto(s_plp, &tid, &pos, &n_plp)) != 0){
		if (pos < beg || pos > end) continue;
		linecount++;
		if(linecount % 50000000 == 0) vrb.bullet("Processed " + stb.str(linecount) + " positions");
		string chr = data->hdr->target_name[tid];
		ase_site temp(chr,pos);
		auto av_it = all_variants.find(temp);
		if (av_it != all_variants.end()){
			unsigned int b_ref = 0, b_alt = 0, b_dis = 0;
			mapping_stats ms;

			//STEP1: Parse sequencing reads
			for (int iread = 0 ; iread < n_plp ; iread ++) {
				bool failed_qc = false;
				const bam_pileup1_t * p = v_plp + iread;

				if (p->b->core.tid < 0 || (p->b->core.flag & BAM_FUNMAP)) {ms.unmapped++; failed_qc = true;}
				else if (p->b->core.flag & BAM_FSECONDARY) {ms.secondary++; failed_qc = true;}
				else if (p->b->core.flag & BAM_FQCFAIL) {ms.fail_qc++; if(!keep_failqc) failed_qc = true;}
				else if ((int)p->b->core.qual < param_min_mapQ) {ms.fail_mapping++; failed_qc = true;}
				else if (p->b->core.flag & BAM_FPAIRED) {
					if (p->b->core.flag & BAM_FMUNMAP) {ms.mate_unmapped++; if(!keep_orphan) failed_qc = true;}
					else if (!(p->b->core.flag & BAM_FPROPER_PAIR)) {ms.not_pp++; if(check_proper_pair) failed_qc = true;}
					//else if ((p->b->core.flag & BAM_FREVERSE) == (p->b->core.flag & BAM_FMREVERSE)) {ms.orientation++; failed_qc = true;}
				}

				if (!failed_qc){
					if (p->b->core.flag & BAM_FDUP){
						ms.duplicate++;
						if(param_dup_rd) continue;
					}
					if (p->indel != 0 && !(p->is_del || p->is_refskip)){
						ms.indel++;
						if(param_rm_indel) continue;
					}
					if (p->is_del || p->is_refskip) {
						ms.skipped++;
						continue;
					}
					if ((p->qpos < p->b->core.l_qseq ? bam_get_qual(p->b)[p->qpos] : 0) < param_min_baseQ) {
						ms.fail_baseq++;
						continue;
					}
					char base = ase_getBase(bam_seqi(bam_get_seq(p->b), p->qpos));
					bool isRef = (base == av_it->ref);
					bool isAlt = (base == av_it->alt);
					if (isRef) b_ref++;
					if (isAlt) b_alt++;
					if (!isRef && !isAlt) b_dis++;
				}
			}
			ase_site current = *av_it;
			current.setCounts(b_ref,b_alt,b_dis,ms);
			passing_variants.push_back(current);
		}
	}
	bam_plp_reset(s_plp);
	bam_plp_destroy(s_plp);
}

void ase_data::readSequences(string fbam) {
	timer current_timer;
	aux_t * data = (aux_t *) malloc (sizeof(aux_t));
	data->iter = NULL;

	vrb.title("Opening BAM file [" + fbam + "]");

	//BAM OPEN
	data->fp = sam_open(fbam.c_str(), "r");
    if (data->fp == 0) vrb.error("Cannot open BAM file!");
    data->hdr = sam_hdr_read(data->fp);
    if (data->hdr == 0) vrb.error("Cannot parse BAM header!");
    hts_idx_t *idx = sam_index_load(data->fp, fbam.c_str());
    if (idx == NULL) vrb.error("Cannot load BAM index!");
	data->mq = param_min_mapQ;

    if (my_regions.size()){
    	for (int reg = 0; reg < my_regions.size(); reg++){
    		data->iter = sam_itr_querys(idx, data->hdr, my_regions[reg].get().c_str()); // set the iterator
    		if (data->iter == NULL) vrb.error("Problem jumping to region [" + my_regions[reg].get() + "]");
    		else vrb.bullet("Scanning region [" + my_regions[reg].get() + "] " + stb.str(reg+1) + " / " + stb.str(my_regions.size()));
    		parseBamMpileup((void *) data);
    	}
    }else parseBamMpileup((void *) data);

	vrb.bullet(stb.str(passing_variants.size()) + " variant found");

    bam_hdr_destroy(data->hdr);
    hts_idx_destroy(idx);
    if (data->fp) sam_close(data->fp);
    hts_itr_destroy(data->iter);
    free(data);
	vrb.bullet("Time taken: " + stb.str(current_timer.abs_time()) + " seconds");
}

void ase_data::calculateRefToAltBias(string olog){
	output_file fdo;
	vrb.title("Calculating reference allele mapping bias");
	if (olog != ""){
		vrb.bullet("Writing failed variants to [" + olog + "]");
		fdo.append(olog);
		if (fdo.fail()) vrb.error("Cannot open file [" + olog +"]");
	}
	vector < unsigned int > all_total_counts;
	map < string , vector <unsigned int> > alleles_total_counts;
	vector < ase_site * > filtered_variants;
	unsigned int filtered_cov = 0 , filtered_bas = 0;
	for (int i = 0 ; i < passing_variants.size(); i++){
		if (passing_variants[i].total_count < param_min_cov_for_ref_alt) { filtered_cov++; if(olog!="") fdo<<"BIAS_COV " << passing_variants[i].sid << endl; continue;}
		if (param_both_alleles_seen_bias && (passing_variants[i].alt_count == 0 || passing_variants[i].ref_count == 0 )) {filtered_bas++; if(olog!="") fdo<<"BIAS_BOTH_ALLELES " << passing_variants[i].sid << endl; continue;}
		filtered_variants.push_back(&passing_variants[i]);
		if (alleles_total_counts.count(passing_variants[i].alleles)) alleles_total_counts[passing_variants[i].alleles].push_back(passing_variants[i].total_count);
		else alleles_total_counts[passing_variants[i].alleles] = vector <unsigned int> (1, passing_variants[i].total_count);
		all_total_counts.push_back(passing_variants[i].total_count);
	}

	if (filtered_cov) vrb.bullet(stb.str(filtered_cov) + " sites with coverage less than " + stb.str(param_min_cov_for_ref_alt));
	if (filtered_bas) vrb.bullet(stb.str(filtered_bas) + " sites where only one allele was observed");

	bool calculate_all = false;
	vector <string> add;
	for (auto it = alleles_total_counts.begin(); it != alleles_total_counts.end(); it++){
		if (it->second.size() < param_min_sites_for_ref_alt) {
			calculate_all = true;
			add.push_back(it->first);
			vrb.bullet("Bias for " + it->first + " will be calculated from all sites since " + stb.str(it->second.size()) + " is too low");
		}
		else{
			sort(alleles_total_counts[it->first].begin(), alleles_total_counts[it->first].end());
			unsigned int max_index = alleles_total_counts[it->first].size() * param_sample;
			unsigned int max = alleles_total_counts[it->first][max_index];
			unsigned int refc = 0, altc = 0;
			for (int i = 0; i < filtered_variants.size(); i++){
				if (filtered_variants[i]->alleles != it->first) continue;
				if (filtered_variants[i]->total_count <= max){
					refc += filtered_variants[i]->ref_count;
					altc += filtered_variants[i]->alt_count;
				}else{
					double ratio = (double) filtered_variants[i]->ref_count / (double) filtered_variants[i]->total_count;
					for (int r = 0; r < max; r++){
						if (rng.getDouble() <= ratio) refc++;
						else altc++;
					}
				}
			}
			ref_to_alt_bias[it->first] = (double) refc / ((double) altc + (double) refc);
			vrb.bullet("Bias for " + it->first + " from " + stb.str(it->second.size()) + " sites is " + stb.str(ref_to_alt_bias[it->first]));
		}
	}
	if (calculate_all){
		sort(all_total_counts.begin(),all_total_counts.end());
		unsigned int max_index = all_total_counts.size() * param_sample;
		unsigned int max = all_total_counts[max_index];
		unsigned int refc = 0, altc = 0;
		for (int i = 0; i < filtered_variants.size(); i++){
			if (filtered_variants[i]->total_count <= max){
				refc += filtered_variants[i]->ref_count;
				altc += filtered_variants[i]->alt_count;
			}else{
				double ratio = (double) filtered_variants[i]->ref_count / (double) filtered_variants[i]->total_count;
				for (int r = 0; r < max; r++){
					if (rng.getDouble() <= ratio) refc++;
					else altc++;
				}
			}
		}
		double all_bias = (double) refc / ((double) altc + (double) refc);
		vrb.bullet("Bias from all sites is " + stb.str(all_bias));
		for (int e = 0; e < add.size(); e++) ref_to_alt_bias[add[e]] = all_bias;
	}
}

void ase_data::calculateASE(string fout , string olog){
	vrb.title("Calculating ASE");
	output_file fdoo(fout);
	if (fdoo.fail()) vrb.error("Cannot open file [" + fout +"]");
	output_file fdo;
	if (olog != ""){
		vrb.bullet("Writing failed variants to [" + olog + "]");
		fdo.append(olog);
		if (fdo.fail()) vrb.error("Cannot open file [" + olog +"]");
	}
	fdoo <<"CHR\tPOS\tID\tREF\tALT\tALLELES\tREF_COUNT\tNONREF_COUNT\tTOTAL_COUNT\tMAR\tOTHER_COUNT\tREF_BIAS\tWEIGHTED_REF_COUNT\tWEIGHTED_NONREF_COUNT\tPVAL\tUNMAPPED\tSECONDARY\tFAILQC\tFAIL_MAPQ\tNOT_PROPER_PAIR\tMATE_UNMAPPED\tWRONG_ORIENTATION\tDUPLICATE\tINDEL\tSKIPPED\tFAIL_BASEQ"<<endl;
	for (int i = 0; i < passing_variants.size(); i++){
		if (passing_variants[i].total_count >= param_min_cov){
			if (param_both_alleles_seen && (passing_variants[i].alt_count == 0 || passing_variants[i].ref_count == 0 )){if(olog!="") fdo<<"ASE_BOTH_ALLELES " << passing_variants[i].sid << endl; continue;};
			if (ref_to_alt_bias[passing_variants[i].alleles] >= 0 && ref_to_alt_bias[passing_variants[i].alleles] <= 1) passing_variants[i].calculatePval(ref_to_alt_bias[passing_variants[i].alleles]);
			else vrb.error("Reference allele mapping bias is incorrect [" + stb.str(ref_to_alt_bias[passing_variants[i].alleles]) + "]");
			fdoo << passing_variants[i] << endl;
		}else if(olog!="") fdo<<"ASE_COV " << passing_variants[i].sid << endl;
	}

}


/*void ase_data::readSequences(string fbam, string fvcf) {
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
		bam_plp_t s_plp = bam_plp_init(ase_read_bam, (void*)data);
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
}*/
