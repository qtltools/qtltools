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

#include "bamstat_data.h"

int bamstat_data::keepRead(bam1_t * b) {
        if (b->core.flag & BAM_FUNMAP) return 1;
        if (b->core.flag & BAM_FSECONDARY) return 2;
        if (b->core.flag & BAM_FQCFAIL) return 3;
        if (b->core.flag & BAM_FDUP) return 4;

        if (b->core.flag & BAM_FPAIRED) {
                if (! (b->core.flag & BAM_FPROPER_PAIR)) return 5;
                if (b->core.flag & BAM_FMUNMAP) return 6;
                if ((b->core.flag & BAM_FREVERSE) == (b->core.flag & BAM_FMREVERSE)) return 7;
        }
        if ((int)b->core.qual < param_min_mapQ) return 8;
        return 0;
}

void bamstat_data::readSequences(string fbam) {
	vrb.title("Opening BAM file ["  + fbam + "]");
	samFile * fd = sam_open(fbam.c_str(), "r");
	if (fd == 0) vrb.error("Failed to open file");
    vrb.bullet("reading header");
    bam_hdr_t * header = sam_hdr_read(fd);
    if (header == 0) vrb.error("Failed to read header");
    vrb.bullet("reading index file");
    hts_idx_t *idx = sam_index_load(fd, fbam.c_str());
    if (idx == NULL) vrb.error("Failed to load index");

    vrb.bullet("STEP1: Iterate across reads");
    bam1_t * b = bam_init1();
	bool pair_end_seq = false;
	unsigned int n_unmap = 0;
	unsigned int n_second = 0;
	unsigned int n_qcfail = 0;
	unsigned int n_dup = 0;
	unsigned int n_proper = 0;
	unsigned int n_munmapped = 0;
	unsigned int n_strand = 0;
	unsigned int n_mappingqual = 0;
	unsigned int n_bothmappingqual = 0;
	unordered_map < string, pair < unsigned char, unsigned char > > M;
	unordered_map < string, pair < unsigned char, unsigned char > > :: iterator itM;
    while (sam_read1(fd, header, b) >= 0) {
        int error_code = keepRead(b);
        pair_end_seq = pair_end_seq || (b->core.flag & BAM_FPAIRED);
        if (error_code == 0) {
        	string uid = string(bam_get_qname(b));
        	itM = M.find(uid);
        	if (itM == M.end()) M.insert(make_pair(string(bam_get_qname(b)), make_pair(1,0)));
            else itM->second.first ++;
            n_keep++;
        } else if (error_code == 1) n_unmap ++;
        else if (error_code == 2) n_second ++;
        else if (error_code == 3) n_qcfail ++;
        else if (error_code == 4) {
        	if (!param_dup_rd) {
            	string uid = string(bam_get_qname(b));
            	itM = M.find(uid);
            	if (itM == M.end()) M.insert(make_pair(string(bam_get_qname(b)), make_pair(1,0)));
                else itM->second.first ++;
                n_keep++;
        	}
        	n_dup ++;
        } else if (error_code == 5) n_proper ++;
		else if (error_code == 6) n_munmapped ++;
		else if (error_code == 7) n_strand ++;
		else if (error_code == 8) n_mappingqual ++;
		n_total_reads ++;
    }
    bam_destroy1(b);

    vrb.bullet("  + #n_total_reads = " + stb.str(n_total_reads));
    vrb.bullet("  - #n_unmap_reads = " + stb.str(n_unmap));
    vrb.bullet("  - #n_secondary_reads = " + stb.str(n_second));
    vrb.bullet("  - #n_qcfail_reads = " + stb.str(n_qcfail));
    if (param_dup_rd) vrb.bullet("  - #n_dup_reads = " + stb.str(n_dup));
    else vrb.bullet("    #n_dup_reads = " + stb.str(n_dup));
    if (pair_end_seq) vrb.bullet("  - #n_improperpairs_reads = " + stb.str(n_proper));
    if (pair_end_seq) vrb.bullet("  - #n_mate_unmapped_reads = " + stb.str(n_munmapped));
    if (pair_end_seq) vrb.bullet("  - #n_strand_pb_reads = " + stb.str(n_strand));
    vrb.bullet("  - #n_mapping_qual_reads = " + stb.str(n_mappingqual));
    vrb.bullet("  = #n_good_reads = " + stb.str(n_keep));

    if (pair_end_seq) {
    	vrb.bullet("  = #n_total_pairs = "  + stb.str(M.size()));
    	for (itM = M.begin() ; itM != M.end() ;) {
    		if (itM->second.first != 2) {
    			itM = M.erase(itM);
    			n_bothmappingqual ++;
    		} else {
    			itM->second.first = 0;
    			++itM;
    		}
    	}
    	vrb.bullet("  - #n_both_not_mapped_pairs = " + stb.str(n_bothmappingqual));
    	vrb.bullet("  = #n_good_pairs = " + stb.str(M.size()));
    } else for (itM = M.begin() ; itM != M.end() ; ++itM) itM->second.first = 0;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    vrb.bullet("STEP2: Iterate across regions");
    b = bam_init1();
    for (int r = 0 ; r < R.size() ; r ++) {
    	string r_str = R[r].toString();
        hts_itr_t *iter = sam_itr_querys(idx, header, r_str.c_str());
        if (iter == NULL) vrb.error("Failed to parse region [" + r_str + "]");
        else {
        	unsigned int n_tmp = 0;
        	while (sam_itr_next(fd, iter, b) >= 0) {
        		string uid = string(bam_get_qname(b));
        		itM = M.find(uid);
        		if (itM != M.end()) {
        			if (b->core.flag & BAM_FPAIRED) {
        				if ((b->core.flag & BAM_FREAD1) && itM->second.first < 255) itM->second.first ++;
        				if ((b->core.flag & BAM_FREAD2) && itM->second.second < 255) itM->second.second ++;
        			} else if (itM->second.first < 255) itM->second.first ++;
        			R[r].n_covering_reads ++;
        		}
        		n_tmp ++;
        	}
        }
        hts_itr_destroy(iter);

        vrb.progress((r+1) * 1.0 / R.size());
    }
    bam_destroy1(b);

    //Counting 'exonic' reads
    vrb.title("Counting pairs overlapping annotations");
    for (itM = M.begin() ; itM != M.end() ; ++ itM) {
        if (itM->second.first > 0) n_overlap_reads ++;
        if (itM->second.second > 0) n_overlap_reads ++;
    }
    vrb.bullet("#n_overlap_reads = " + stb.str(n_overlap_reads) +  " / " + stb.str(n_keep));

    //Counting covered regions
    vrb.title("Counting annotations overlapped by reads");
    for (int r = 0 ; r < R.size() ; r ++) {
    	if (R[r].n_covering_reads > 0) n_overlap_annotations ++;
    }
    vrb.bullet("#n_overlap_annotations = " + stb.str(n_overlap_annotations) +  " / " + stb.str(R.size()));

    hts_idx_destroy(idx);
    bam_hdr_destroy(header);
    sam_close(fd);

}
