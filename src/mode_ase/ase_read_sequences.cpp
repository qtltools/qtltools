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
	bool keep_orphan,check_orientation,check_proper_pair;
	int fflag;
} aux_t;

/*static int ase_read_bam(void *data, bam1_t *b) {
	aux_t * aux = (aux_t*) data;
	int ret, skip = 0;
	do{
		ret =  (aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b));
		if (ret < 0) break;
		if (b->core.tid < 0 || (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY ))) { skip = 1; continue;}
		skip = 0;
		if (b->core.qual < aux->mq)  {skip = 1;} //having the mapq filter here improves performance drastically
	}while(skip);
	return ret;
}*/


static int ase_read_bam(void *data, bam1_t *b) {
	aux_t * aux = (aux_t*) data;
	int ret, skip = 0;
	do{
		ret =  (aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b));
		if (ret < 0) break;
		else if (b->core.tid < 0 || (b->core.flag & aux->fflag) || b->core.qual < aux->mq) {skip = 1; continue;}
		else if(b->core.flag & BAM_FPAIRED){
			if (b->core.flag & BAM_FMUNMAP){
				if (!aux->keep_orphan) { skip = 1;continue;}
			}else{
				if (aux->check_orientation && ((b->core.tid !=b->core.mtid) ||  //must be on the same chr
					((b->core.flag & BAM_FREVERSE) && (b->core.flag & BAM_FMREVERSE)) || // read 1 and 2 cannot be on the same strand
					((b->core.pos < b->core.mpos) && (b->core.flag & BAM_FREVERSE)) || //read 1 must be on the +ve strand
					((b->core.pos > b->core.mpos) && !(b->core.flag & BAM_FREVERSE)) //read 2 must be on the -ve strand
				   )) { skip = 1; continue;}
				if (aux->check_proper_pair && !(b->core.flag & BAM_FPROPER_PAIR)) { skip = 1; continue;}
			}
		}
		skip = 0;
	}while(skip);
	return ret;
}


/*static int mplp_func(void *data, bam1_t *b){
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
}*/

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
	//bam_plp_t s_plp = legacy_options ? bam_plp_init(mplp_func, (void*)data) : bam_plp_init(ase_read_bam, (void*)data) ;
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
	unsigned long linecount = 0, prevstart = 0, prevend = 0;
	string prevchr = "";
	timer local;
	unsigned int max_depth_limit = max_depth * 0.90;
	while((v_plp = bam_plp_auto(s_plp, &tid, &pos, &n_plp)) != 0){
		linecount++;
		if (linecount % 50000000 == 0) vrb.bullet(stb.str(linecount) + " positions read, " + stb.str((double) linecount / (double) local.abs_time()) + " per second.");
		if (pos < beg || pos > end) continue;
		bool depth_prb = false;
		string chr = data->hdr->target_name[tid];

		//Check for position exceeding max depth
		if (chr != prevchr){
			if (region_length == 0 && !bam_region.isSet( ) && ase_chrs == bam_chrs){
				vrb.bullet("Processing chromosome [" + chr + "]");
			}
			if (prevstart){
				depth_exceeded.push_back(ase_basic_block(prevchr,prevstart > depth_flag_length ? prevstart - depth_flag_length : 1 , prevend + depth_flag_length <= ULONG_MAX ? prevend + depth_flag_length : ULONG_MAX));
			}
			if (n_plp >= max_depth_limit){
				depth_prb = true;
				prevstart = prevend = pos+1;
			}else prevstart = prevend = 0;
			prevchr = chr;
		}else if (n_plp >= max_depth_limit){
			depth_prb = true;
			if(prevstart){
				if (pos == prevend){
					prevend = pos+1;
				}else{
					depth_exceeded.push_back(ase_basic_block(prevchr,prevstart > depth_flag_length ? prevstart - depth_flag_length : 1 , prevend + depth_flag_length <= ULONG_MAX ? prevend + depth_flag_length : ULONG_MAX));
					prevstart = prevend = pos+1;
				}
			}else prevstart = prevend = pos+1;
		}else if (prevstart){
			depth_exceeded.push_back(ase_basic_block(prevchr,prevstart > depth_flag_length ? prevstart - depth_flag_length : 1 , prevend + depth_flag_length <= ULONG_MAX ? prevend + depth_flag_length : ULONG_MAX));
			prevstart = prevend = 0;
		}


		ase_site temp(chr,pos);
		auto av_it = all_variants.find(temp);
		if (av_it != all_variants.end()){
			unsigned int b_ref = 0, b_alt = 0, b_dis = 0;
			mapping_stats ms;
			set <string> as;
			//STEP1: Parse sequencing reads
			if (print_warnings && depth_prb) vrb.warning(av_it->sid + " depth " + stb.str(n_plp) + " is >= 90% max-depth, potential data loss!");
			ms.depth = n_plp;
			for (int iread = 0 ; iread < n_plp ; iread ++) {
				const bam_pileup1_t * p = v_plp + iread;

				int baseq = p->qpos < p->b->core.l_qseq ? bam_get_qual(p->b)[p->qpos] : 0;
				if (illumina13) baseq = baseq > 31 ? baseq - 31 : 0;

				if (p->is_del || p->is_refskip) ms.skipped++; //skipped read
				else if (baseq < param_min_baseQ) ms.fail_baseq++; //fail base q
				else if (param_rm_indel && p->indel != 0) ms.indel++; //read containing an indel
				else{
					char base = getBase(bam_seqi(bam_get_seq(p->b), p->qpos));
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
			//check ASE site
			if (current.other_count && current.ref_count == 0 && current.alt_count == 0) {
				if (print_warnings) vrb.warning("No ref or alt allele for " + current.getName() + " in " + stb.str(current.other_count) + " sites");
				current.concern += "NRA,";
			}else if (current.other_count){
				current.concern += "DA,";
				if (current.other_count > current.alt_count) {
					if (print_warnings) vrb.warning("More discordant alleles than alt alleles for " + current.getName() + " " + stb.str(current.other_count) + " > " + stb.str(current.alt_count));
					current.concern += "MDTA,";
				}
				if (current.other_count > current.ref_count) {
					if (print_warnings) vrb.warning("More discordant alleles than ref alleles for " + current.getName() + " " + stb.str(current.other_count) + " > " + stb.str(current.ref_count));
					current.concern += "MDTR,";
				}
			}
			if (!current.ref_count || !current.alt_count) {current.concern += "BANS,";}
			else if (current.mar < 0.02) {current.concern += "LMAR,";}
			passing_variants.push_back(current);
		}
	}//while end

	if (prevstart){
		depth_exceeded.push_back(ase_basic_block(prevchr,prevstart > depth_flag_length ? prevstart - depth_flag_length : 1 , prevend + depth_flag_length <= ULONG_MAX ? prevend + depth_flag_length : ULONG_MAX));
	}
	bam_plp_reset(s_plp);
	bam_plp_destroy(s_plp);
	vrb.bullet(stb.str(linecount) + " positions read, " + stb.str((double) linecount / (double) local.abs_time()) + " per second.");
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
    data->fflag = (BAM_FUNMAP | BAM_FSECONDARY);
    if (!keep_failqc) data->fflag |= BAM_FQCFAIL;
    if (param_dup_rd) data->fflag |= BAM_FDUP;
	data->mq = param_min_mapQ;
	data->keep_orphan = keep_orphan;
	data->check_orientation = check_orientation;
	data->check_proper_pair = check_proper_pair;



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
	if (passing_variants.size() == 0) vrb.leave("Cannot find usable variants in target region!");
	if (depth_exceeded.size()){
		vector <ase_basic_block> temp;
		mergeContiguousBlocks(depth_exceeded, temp);
		depth_exceeded = temp;
		//for (int i = 0 ; i < depth_exceeded.size() ; i++) cerr << depth_exceeded[i].get_string() << endl;
		unsigned int n = 0;
		for (int i = 0 ; i < passing_variants.size() ; i++){
			if (ase_basic_block(passing_variants[i]).find_this_in_bool(depth_exceeded)){
				n++;
				passing_variants[i].concern += "PD,";
			}
		}
		if (n) vrb.warning("Pileup depth is > 90% of max depth in " + stb.str(depth_exceeded.size()) + " regions. Due to intricacies of HTSlib, and to be on the safe side, " + stb.str(n) + " variants within " + stb.str(depth_flag_length) + " bp around these regions will be labeled with the PD concern. HIGHLY RECOMMENDED to increase --max-depth and rerun the analysis.");
	}
}



