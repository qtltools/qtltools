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
		if (b->core.tid < 0 || (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP))) { skip = 1; continue;}
		skip = 0;
		if (b->core.qual < aux->mq)  {skip = 1;}
		else if ((b->core.flag&BAM_FPAIRED) && !(b->core.flag&BAM_FPROPER_PAIR)) { skip = 1;}

	}while(skip);
	return ret;
}

/*static inline void pileup_seq(const bam_pileup1_t * p){
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
}*/



void ase_data::parseBam(void * d){
	aux_t * data = (aux_t *) d;
	//Pile up reads
	const bam_pileup1_t * v_plp;
	int n_plp = 0, tid, pos;
	bam_plp_t s_plp = legacy_options ? bam_plp_init(mplp_func, (void*)data) : bam_plp_init(ase_read_bam, (void*)data);
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
			set <string> as;
			//STEP1: Parse sequencing reads
			if (n_plp >= max_depth * 0.95) vrb.warning(av_it->sid + " depth " + stb.str(n_plp) + " is >= 95% max-depth, potential data loss!");
			for (int iread = 0 ; iread < n_plp ; iread ++) {
				bool failed_qc = false;
				const bam_pileup1_t * p = v_plp + iread;

				int baseq = p->qpos < p->b->core.l_qseq ? bam_get_qual(p->b)[p->qpos] : 0;
				if (illumina13) baseq = baseq > 31 ? baseq - 31 : 0;

				if (p->b->core.tid < 0 || (p->b->core.flag & BAM_FUNMAP)) {ms.unmapped++; failed_qc = true;}
				else if (p->b->core.flag & BAM_FSECONDARY) {ms.secondary++; failed_qc = true;}
				else if ((int)p->b->core.qual < param_min_mapQ) {ms.fail_mapping++; failed_qc = true;}
				else if (p->is_del || p->is_refskip) {ms.skipped++; failed_qc = true;}
				else if (baseq < param_min_baseQ) {ms.fail_baseq++; failed_qc = true;}
				else{
					if (p->b->core.flag & BAM_FQCFAIL) {ms.fail_qc++; if(!keep_failqc) failed_qc = true;}
					if (!failed_qc && (p->b->core.flag & BAM_FDUP)) {ms.duplicate++; if(param_dup_rd) failed_qc = true;}
					if (!failed_qc && p->indel != 0) {ms.indel++; if(param_rm_indel) failed_qc = true;}
					if (!failed_qc && (p->b->core.flag & BAM_FPAIRED)) {
						if (p->b->core.flag & BAM_FMUNMAP) {ms.mate_unmapped++; if(!keep_orphan) failed_qc = true;}
						else if (
								(p->b->core.tid != p->b->core.mtid) ||  //must be on the same chr
								((p->b->core.flag & BAM_FREVERSE) && (p->b->core.flag & BAM_FMREVERSE)) || // read 1 and 2 cannot be on the same strand
								((p->b->core.pos < p->b->core.mpos) && (p->b->core.flag & BAM_FREVERSE)) || //read 1 must be on the +ve strand
								((p->b->core.pos > p->b->core.mpos) && !(p->b->core.flag & BAM_FREVERSE)) //read 2 must be on the -ve strand
								) { ms.orientation++; if(check_orientation) failed_qc = true;}
						else if (!(p->b->core.flag & BAM_FPROPER_PAIR)) {ms.not_pp++; if(check_proper_pair) failed_qc = true;}
					}
				}

				if (!failed_qc){
					char base = ase_getBase(bam_seqi(bam_get_seq(p->b), p->qpos));
					as.insert(string(1,base));
					bool isRef = (base == av_it->ref);
					bool isAlt = (base == av_it->alt);
					if (isRef) b_ref++;
					if (isAlt) b_alt++;
					if (!isRef && !isAlt) b_dis++;
				}
			}
			ase_site current = *av_it;
			current.setCounts(b_ref,b_alt,b_dis,as,ms);
			if (n_plp >= max_depth * 0.95) current.concern += "PD,";
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
    		data->iter = sam_itr_querys(idx, data->hdr, my_regions[reg].get_string().c_str()); // set the iterator
    		if (data->iter == NULL) vrb.error("Problem jumping to region [" + my_regions[reg].get_string() + "]");
    		else vrb.bullet("Scanning region [" + my_regions[reg].get_string() + "] " + stb.str(reg+1) + " / " + stb.str(my_regions.size()));
    		parseBam((void *) data);
    	}
    }else parseBam((void *) data);

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
	for (auto it = all_allele_combinations.begin(); it != all_allele_combinations.end(); it++){
		if (alleles_total_counts.count(*it) == 0 || alleles_total_counts[*it].size() < param_min_sites_for_ref_alt) {
			calculate_all = true;
			add.push_back(*it);
			vrb.bullet("Bias for " + *it + " will be calculated from all sites since " + stb.str(alleles_total_counts.count(*it) ? alleles_total_counts[*it].size() : 0) + " is too low");
		}
		else{
			sort(alleles_total_counts[*it].begin(), alleles_total_counts[*it].end());
			unsigned int max_index = alleles_total_counts[*it].size() * param_sample - 1;
			unsigned int max = alleles_total_counts[*it][max_index];
			unsigned int refc = 0, altc = 0, sc = 0;
			for (int i = 0; i < filtered_variants.size(); i++){
				if (filtered_variants[i]->alleles != *it) continue;
				if (filtered_variants[i]->total_count <= max){
					refc += filtered_variants[i]->ref_count;
					altc += filtered_variants[i]->alt_count;
				}else{
					sc++;
					double ratio = (double) filtered_variants[i]->ref_count / (double) filtered_variants[i]->total_count;
					for (int r = 0; r < max; r++){
						if (rng.getDouble() <= ratio) refc++;
						else altc++;
					}
				}
			}
			ref_to_alt_bias[*it] = (double) refc / ((double) altc + (double) refc);
			vrb.bullet("Bias for " + *it + " from " + stb.str(alleles_total_counts[*it].size()) + " sites is " + stb.str(ref_to_alt_bias[*it]));
			if (sc) vrb.bullet("Subsampled " + stb.str(sc) + " sites to " + stb.str(max));
		}
	}
	if (calculate_all){
		sort(all_total_counts.begin(),all_total_counts.end());
		unsigned int max_index = all_total_counts.size() * param_sample - 1;
		unsigned int max = all_total_counts[max_index];
		unsigned int refc = 0, altc = 0, sc = 0;
		for (int i = 0; i < filtered_variants.size(); i++){
			if (filtered_variants[i]->total_count <= max){
				refc += filtered_variants[i]->ref_count;
				altc += filtered_variants[i]->alt_count;
			}else{
				sc++;
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
		if (sc) vrb.bullet("Subsampled " + stb.str(sc) + " sites to " + stb.str(max));
	}
}

void ase_data::calculateASE(string fout , string olog){
	vrb.title("Calculating ASE and writing to [" + fout + "]");
	output_file fdoo(fout);
	if (fdoo.fail()) vrb.error("Cannot open file [" + fout +"]");
	output_file fdo;
	if (olog != ""){
		vrb.bullet("Writing failed variants to [" + olog + "]");
		fdo.append(olog);
		if (fdo.fail()) vrb.error("Cannot open file [" + olog +"]");
	}
	fdoo <<"INDIVIDUAL\tRSID\tCHR\tPOS\tALLELES\tBOTH_ALLELES_SEEN\tMIN_ALLELE_RATIO\tREF_COUNT\tNONREF_COUNT\tTOTAL_COUNT\tWEIGHTED_REF_COUNT\tWEIGHTED_NONREF_COUNT\tWRC_MINUS_WNC\tALLELES_SEEN\tREF_ALLELE\tALT_ALLELE\tOTHER_COUNT\tREF_RATIO\tPVALUE\tCONCERN\tEXON_INFO";
	if (print_stats) fdoo <<"\tUNMAPPED\tSECONDARY\tFAIL_MAPQ\tSKIPPED\tFAIL_BASEQ\tFAILQC\tDUPLICATE\tINDEL\tMATE_UNMAPPED\tWRONG_ORIENTATION\tNOT_PROPER_PAIR";
	fdoo <<endl;
	for (int i = 0; i < passing_variants.size(); i++){
		if (passing_variants[i].total_count >= param_min_cov){
			if (param_both_alleles_seen && (passing_variants[i].alt_count == 0 || passing_variants[i].ref_count == 0 )){if(olog!="") fdo<<"ASE_BOTH_ALLELES " << passing_variants[i].sid << endl; continue;};
			if (ref_to_alt_bias[passing_variants[i].alleles] >= 0 && ref_to_alt_bias[passing_variants[i].alleles] <= 1) passing_variants[i].calculatePval(ref_to_alt_bias[passing_variants[i].alleles]);
			else vrb.error("Reference allele mapping bias is incorrect [" + stb.str(ref_to_alt_bias[passing_variants[i].alleles]) + "]");
			assignGenesToAseSite(passing_variants[i]);
			fdoo << sample_id[0] << "\t" << passing_variants[i];
			if (print_stats) fdoo << "\t" << passing_variants[i].stats;
			fdoo << endl;
		}else if(olog!="") fdo<<"ASE_COV " << passing_variants[i].sid << endl;
	}

}


