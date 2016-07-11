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

#include "match_data.h"

typedef struct {     			// auxiliary data structure
	samFile * fp;				// the file handle
	bam_hdr_t * hdr;			// the file header
	hts_itr_t * iter;			// NULL if a region not specified
	int min_mapQ;				// mapQ filter
	bool dup_rd;			// Do we consider duplicate read or not?
} aux_t;

static int read_bam(void *data, bam1_t *b) {
	aux_t * aux = (aux_t*) data;
	int ret;
	while (1) {
		ret = aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
		if (ret < 0) break;
		if (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL)) continue;
		if (!aux->dup_rd && (b->core.flag & BAM_FDUP)) continue;
		if (b->core.flag & BAM_FPAIRED) {
			if (! (b->core.flag & BAM_FPROPER_PAIR)) continue;
			if (b->core.flag & BAM_FMUNMAP) continue;
			if ((b->core.flag & BAM_FREVERSE) == (b->core.flag & BAM_FMREVERSE)) continue;
		}
		if ((int)b->core.qual < aux->min_mapQ) continue;
		break;
    }
    return ret;
}

void match_data::readSequences(string fbam) {
	aux_t * data = (aux_t *) malloc (sizeof(aux_t));

	vrb.title("Processing BAM file [" + fbam + "]");
	data->fp = sam_open(fbam.c_str(), "r");
    if (data->fp == 0) vrb.error("Cannot open file!");
    data->min_mapQ = param_min_mapQ;
    data->dup_rd = param_dup_rd;
    data->hdr = sam_hdr_read(data->fp);
    if (data->hdr == 0) vrb.error("Cannot parse header!");
    hts_idx_t *idx = sam_index_load(data->fp, fbam.c_str());
    if (idx == NULL) vrb.error("Cannot load index!");

	//Loop across regions
	for (int reg = 0; reg < regions.size() ; reg ++) {

		//Jump to region
		data->iter = sam_itr_querys(idx, data->hdr, regions[reg].c_str()); // set the iterator
		if (data->iter == NULL) vrb.error("Problem jumping to region [" + regions[reg] + "]");
		else vrb.bullet("scanning region [" + regions[reg] + "]");

		int beg = data->iter->beg;
		int end = data->iter->end;

		//Pile up reads
		const bam_pileup1_t * v_plp;
		int n_plp = 0, tid, pos, i_site = 0;
		bam_plp_t s_plp = bam_plp_init(read_bam, (void*)data);
		while (((v_plp = bam_plp_auto(s_plp, &tid, &pos, &n_plp)) != 0) && i_site < sites[reg].size()) {
		    int chr = bam_name2id(data->hdr, sites[reg][i_site].chr.c_str());
			if (pos < beg || pos >= end) continue;
			while (i_site < sites[reg].size() && (chr != tid || pos > sites[reg][i_site].pos)) { i_site ++; }
			if (tid == chr && pos == sites[reg][i_site].pos) {
				for (int red = 0 ; red < n_plp ; red ++) {
					const bam_pileup1_t * p = v_plp + red;
					sites[reg][i_site].ctot ++;
					if (p->is_del || p->is_refskip || p->indel == 1) sites[reg][i_site].cdel ++;
					else if (bam_get_qual(p->b)[p->qpos] < param_min_baseQ) sites[reg][i_site].cqual ++;
					else {
						char base = match_getBase(bam_seqi(bam_get_seq(p->b), p->qpos));
						bool isRef = (base == sites[reg][i_site].ref);
						bool isAlt = (base == sites[reg][i_site].alt);
						if (isRef) sites[reg][i_site].cref ++;
						if (isAlt) sites[reg][i_site].calt ++;
						if (!isRef && !isAlt) sites[reg][i_site].cdisc ++;
					}
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
}
