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

#include "quan_data.h"


inline vector < my_gene > quan2_data::getOverlappingGenesQ(my_block & bl){
	if (!genome.count(bl.chr)) return vector < my_gene >(0);
	set <int> list;
	set <int>::iterator it;
	for (int sb = bl.start/binsize; sb <= bl.end/binsize; sb++){
		if (genome[bl.chr].count(sb)){
			for (int i = 0 ; i < genome[bl.chr][sb].size(); i++){
				if (genes[genome[bl.chr][sb][i]].overlap(bl.start,bl.end)) list.insert(genome[bl.chr][sb][i]);
			}
		}
	}
	vector < my_gene > res(0);
	for(it=list.begin(); it != list.end();it++) res.push_back(genes[(*it)]);
	return res;
}

my_cont_blocks quan2_data::keep_read(bam1_t *b) {
	my_cont_blocks block; //Filter is set to pass by default so no need to set it to pass explicitly here
	block.name = bam_get_qname(b);
	block.core = b->core;

	//Check if secondary first so that we can ignore this alingment
	if (b->core.flag & BAM_FSECONDARY) {
		block.filter = SECD;
#ifdef DEBUG
		cerr << "SECD\t" << block.name << endl;
#endif
		return block;
	}
	//Check if mapped
	if (b->core.flag & BAM_FUNMAP) {
		block.filter = UNMAP;
#ifdef DEBUG
		cerr << "UNMAP\t" << block.name << endl;
#endif
		return block;
	}
	//Check if fails qc
	if(filter.fail_qc && (b->core.flag & BAM_FQCFAIL)){
		block.filter = FAILQC;
#ifdef DEBUG
		cerr << "FAILQC\t" << block.name << endl;
#endif
		return block;
	}
	//Check if duplicate
	if (filter.dup_remove && (b->core.flag & BAM_FDUP)) {
		block.filter = DUP;
#ifdef DEBUG
		cerr << "DUP\t" << block.name << endl;
#endif
		return block;
	}
	//Check the mapping quality
	if ((int)b->core.qual < filter.min_mapQ) {
		block.filter = MAPQ;
#ifdef DEBUG
		cerr << "MAPQ\t" << block.name << endl;
#endif
		return block;
	}
	//If paired check basic orientation and if the mate is mapped
	if (b->core.flag & BAM_FPAIRED) {
		if ((b->core.flag & BAM_FREVERSE) == (b->core.flag & BAM_FMREVERSE) || (b->core.flag & BAM_FMUNMAP) ){
			block.filter = NPP;
#ifdef DEBUG
			cerr << "NPP\t" << block.name << endl;
#endif
			return block;
		}
		if (filter.proper_pair && !(b->core.flag & BAM_FPROPER_PAIR)){
			block.filter = NPP;
#ifdef DEBUG
			cerr << "NPP\t" << block.name << endl;
#endif
			return block;
		}
	}
	//This is adapted from sam.c in htslib original author Heng Li
	if (filter.max_mismatch_count >= 0 || filter.max_mismatch_count_total >= 0){
		unsigned int mc = 0;
		unsigned int count = 0;
		bool NMfound = false;
		bool nonInt = false;
		stats.checkedNM++;
		uint8_t *s = bam_get_aux(b);
		while (s+4 <= b->data + b->l_data) {
			uint8_t type, key[3];
			key[0] = s[0]; key[1] = s[1]; key[2] = '\0';
			string keys((const char*)key);
			s += 2; type = *s++;
			//cerr << block.name << " "<< type << " " << key << endl;
			if (type == 'A') {
				if (keys=="NM"){
					//mc = *s;
					//NMfound = true;
					count++;
					nonInt = true;
				}
				++s;
			} else if (type == 'C') {
				if (keys=="NM"){
					mc = *s;
					NMfound = true;
					count++;
					//nonInt = true;
				}
				++s;
			} else if (type == 'c') {
				if (keys=="NM"){
					mc = *(int8_t*)s;
					NMfound = true;
					count++;
					//nonInt = true;
				}
				++s;
			} else if (type == 'S') {
				if (s+2 <= b->data + b->l_data) {
					if (keys=="NM"){
						mc = *(uint16_t*)s;
						NMfound = true;
						count++;
						//nonInt = true;
					}
					s += 2;
				} else break;
			} else if (type == 's') {
				if (s+2 <= b->data + b->l_data) {
					if (keys=="NM"){
						mc = *(int16_t*)s;
						NMfound = true;
						count++;
						//nonInt = true;
					}
					s += 2;
				} else break;
			} else if (type == 'I') {
				if (s+4 <= b->data + b->l_data) {
					if (keys=="NM"){
						mc = *(uint32_t*)s;
						NMfound = true;
						count++;
					}
					s += 4;
				} else break;
			} else if (type == 'i') {
				if (s+4 <= b->data + b->l_data) {
					if (keys=="NM"){
						mc = *(int32_t*)s;
						NMfound = true;
						count++;
					}
					s += 4;
				} else break;
			} else if (type == 'f') {
				if (s+4 <= b->data + b->l_data) {
					if (keys=="NM"){
						//mc = *(float*)s;
						//NMfound = true;
						count++;
						nonInt = true;
					}
					s += 4;
				} else break;

			} else if (type == 'd') {
				if (s+8 <= b->data + b->l_data) {
					if (keys=="NM") {
						//mc = *(double*)s;
						//NMfound = true;
						count++;
						nonInt = true;
					}
					s += 8;
				}else break;
			} else if (type == 'Z' || type == 'H') {
				if (keys == "NM"){
					count++;
					nonInt = true;
				}
				while (s < b->data + b->l_data && *s) s++;
				if (s >= b->data + b->l_data)
					break;
				++s;
			} else if (type == 'B') {
				if (keys == "NM"){
					count++;
					nonInt = true;
				}
				uint8_t sub_type = *(s++);
				int32_t n;
				memcpy(&n, s, 4);
				s += 4; // no point to the start of the array
				if (s + n >= b->data + b->l_data)
					break;
				for (int i = 0; i < n; ++i) { // FIXME: for better performance, put the loop after "if"
					if ('c' == sub_type)      {  ++s; }
					else if ('C' == sub_type) {  ++s; }
					else if ('s' == sub_type) {  s += 2; }
					else if ('S' == sub_type) {  s += 2; }
					else if ('i' == sub_type) {  s += 4; }
					else if ('I' == sub_type) {  s += 4; }
					else if ('f' == sub_type) {  s += 4; }
				}
			}
		}
		block.mmc = mc;
		block.NMfound = NMfound;
		if (!NMfound) stats.noNM++;
		if (nonInt) stats.nonIntNM++;
		if (count > 1) stats.multiNM++;
		//Check the number of mismatches
		if ((!filter.fraction_mm && mc > filter.max_mismatch_count) || (filter.fraction_mm && (double) mc / (double) (b->core.l_qseq) > filter.max_mismatch_count)){
			block.filter = MISM;
#ifdef DEBUG
			cerr << "MISM\t" << block.name << endl;
#endif
			return block;
		}
	}

	//If we made it this far then we are OK thus break the alignment into continuous blocks based on the cigar string
    const bam1_core_t *c = &b->core;
    if (c->n_cigar) { // cigar
        unsigned int bS = b->core.pos+1;
        unsigned int bL = 0;
        uint32_t *cigar = bam_get_cigar(b);
#ifdef DEBUG
        cerr<< block.name << "\t" << b->core.flag << "\t" << c->n_cigar;
#endif
        for (int i = 0; i < c->n_cigar; ++i) {
            int l = bam_cigar_oplen(cigar[i]);
            char c = bam_cigar_opchr(cigar[i]);
#ifdef DEBUG
            cerr << " " << l << c;
#endif
            if(c=='S' || c =='H' || c =='I' || c=='P') continue;
            else if ((c=='N' || c=='D') && l){
                if (bL){
					block.starts.push_back(bS);
					block.ends.push_back(bS+bL-1);
					block.lengths.push_back(bL);
					block.base_count.push_back(vector <my_basic_overlap> (0));
					block.read_length+=bL;
                }
                bS+=bL+l;
                bL = 0;
            }else bL += l;
        }
#ifdef DEBUG
        cerr << endl;
#endif
        if(bL){
			block.starts.push_back(bS);
			block.ends.push_back(bS+bL-1);
			block.lengths.push_back(bL);
			block.base_count.push_back(vector <my_basic_overlap> (0));
			block.read_length+=bL;
        }
    }
#ifdef DEBUG
        cerr << block << endl;
#endif
    return(block);
}


void quan2_data::readBam(string fbam){
	vrb.title("Opening BAM file ["  + fbam + "]");
	samFile * fd = sam_open(fbam.c_str(), "r");
	if (fd == 0) vrb.error("Failed to open file");
    vrb.bullet("reading header");
    bam_hdr_t * header = sam_hdr_read(fd);
    if (header == 0) vrb.error("Failed to read header");
    vrb.bullet("reading index file");
    hts_idx_t *idx = sam_index_load(fd, fbam.c_str());
    if (idx == NULL) vrb.error("Failed to load index");
    vrb.bullet("reading BAM file");

    map < string, my_cont_blocks> read_sink;
    bam1_t * b = bam_init1();
    unsigned int line_count = 0;
    hts_itr_t * iter = NULL;
    if(region.isSet()){
    	iter = sam_itr_querys(idx, header, region.get().c_str());
    	if(iter==NULL){
    	    bam_destroy1(b);
    	    hts_itr_destroy(iter);
    	    hts_idx_destroy(idx);
    	    bam_hdr_destroy(header);
    	    sam_close(fd);
    		vrb.error("Problem jumping to region [" + region.get() + "]");
    	}
    }
    //GET THE NUMBER OF READS FROM INDEX
    uint64_t totr = 0;
    for (int i = 0; i < header->n_targets; ++i) {
        uint64_t u, v;
        hts_idx_get_stat(idx, i, &u, &v);
        totr+= u + v;
    }
    totr += hts_idx_get_n_no_coor(idx);
    if (!region.isSet()) vrb.bullet("Expecting " + stb.str(totr) + " lines");
    //while (sam_read1(fd, header, b) >= 0) {
    int ret;
    while (1){
		ret = iter ? sam_itr_next(fd, iter, b) : sam_read1(fd, header, b);
		if (ret < 0) break;
    	line_count++;
    	if (!region.isSet()) vrb.progress( ((float) line_count) / (float) totr);
    	else if (line_count % 1000000 == 0) vrb.bullet(stb.str(line_count) + " lines read");

    	my_cont_blocks B = keep_read(b);

    	//Skip if it is a secondary alignment
    	if (B.filter == SECD){
    		stats.secondary++;
    		continue;
    	}
        if (b->core.tid >= 0){ //if this is less than 0 it means it is unmapped
        	if(b->core.tid >= header->n_targets) vrb.error("Chromosome for [" + B.name + "] is not found in the BAM header");
        	B.chr = string(header->target_name[b->core.tid]);
        }
    	stats.total++;
    	//Check if paired end
    	if (b->core.flag & BAM_FPAIRED){
    		//Check if we encountered the read before
            if(read_sink.count(B.name)){
                my_cont_blocks A = read_sink[B.name]; //Get the first mate
                read_sink.erase(B.name); //We are done with this read
                //Check if any of the reads fail filters
                if (A.filter == UNMAP || B.filter == UNMAP){
                	stats.unmapped +=2;
                	continue;
                }
                if (A.filter == FAILQC || B.filter == FAILQC){
                	stats.failqc +=2;
                	continue;
                }
                if (A.filter == DUP || B.filter == DUP){
                	stats.dup +=2;
                	continue;
                }
                if (A.filter == MAPQ || B.filter == MAPQ){
                	stats.mapQ +=2;
                	continue;
                }
                if (A.filter == NPP || B.filter == NPP){
                	stats.unpaired +=2;
                	continue;
                }
                if (A.filter == MISM || B.filter == MISM){
                	stats.mismatch +=2;
                	continue;
                }
                //Check if the total mismatches in the pair are acceptable
                if (filter.max_mismatch_count_total >= 0 && ((!filter.fraction_mmt && A.mmc + B.mmc > filter.max_mismatch_count_total) || (filter.fraction_mmt && (double) (A.mmc + B.mmc) / (double) (A.core.l_qseq + B.core.l_qseq) > filter.max_mismatch_count_total) )){
#ifdef DEBUG
                	cerr << "MISMT\t" << B.name<<endl;
#endif
                	stats.mismatch+=2;
                	continue;
                }
                //We are good so far, let's check if we are truly correctly orientated
                if(A.core.tid != B.core.tid || (A.core.flag & BAM_FREVERSE) || !(B.core.flag & BAM_FREVERSE) || ( filter.old_wrong_split && B.core.pos <= A.core.pos) || ( !filter.old_wrong_split && B.core.pos < A.core.pos)){
#ifdef DEBUG
                	cerr << "NPPP\t" << B.name << "\tc:" << (A.core.tid != B.core.tid) << "\t1s:" << (A.core.flag & BAM_FREVERSE) << "\t2s:" << !(B.core.flag & BAM_FREVERSE) << "\to:" << ( filter.old_wrong_split && B.core.pos <= A.core.pos) << "\tn:" << ( !filter.old_wrong_split && B.core.pos < A.core.pos) << "\t" << A.core.flag << "\t" << B.core.flag << endl;
#endif
                	stats.unpaired+=2;
                    continue;
                }
                //Passed all filters, good on the reads
                stats.good+=2;
                if (filter.merge) A.merge(B); //if we are merging overlapping mate pair do that
                if (A.merged || B.merged ) stats.merged +=2;
                //Check we overlap with a gene
                my_block temp = my_block(A.chr, A.starts[0], B.ends.back());
                vector < my_gene > overlapping_genes = getOverlappingGenesQ(temp);
                if (overlapping_genes.size()){
                	//overlapping genes found
                	bool both_found = false;
                	int total_genes_found = 0;
                	double total_read1_count = 0.0 , total_read2_count = 0.0;
                    for (int g = 0 ; g < overlapping_genes.size(); g++){
#ifdef DEBUG
                    	cerr << "CHECKING " << A << " " << B << " in " << temp << " for " << overlapping_genes[g].gene_id <<" " << g+1 << "/" << overlapping_genes.size()<< endl;
#endif
                        vector < unsigned int > exon_overlap1,exon_overlap2;
                        vector <double> exon_count1,exon_count2;
                        bool all_found1 = true, all_found2 = true , any_found1 = false, any_found2 = false;
                        double exon_overlap1_length_total = 0.0 , exon_overlap2_length_total = 0.0;
                        //check if the first mate is overlapping any exons
                        //iterate through all blocks without gaps
                        for (int i = 0 ; i < A.starts.size() ; i++){
                            int idx = -1;
                            if (overlapping_genes[g].overlap(A.starts[i],A.ends[i])){ //if this block is overlapping with a gene
                                for (int e = 0 ; e < overlapping_genes[g].exons.size();e++){ // iterate trough all the exons of the gene
                                    if(overlapping_genes[g].exons[e].overlap(A.starts[i],A.ends[i])){ //if we overlap with this exon
                                        idx = e;
                                        any_found1 = true; //at least one block of this read overlaps with an exon
                                        exon_overlap1.push_back(idx); //record which exon we overlap
                                        unsigned int rs = max(overlapping_genes[g].exons[idx].start, A.starts[i]); //restrictive start position of the overlap
                                        unsigned int re = min(overlapping_genes[g].exons[idx].end, A.ends[i]); //restrictive end position of the overlap
                                        if ( filter.old_wrong_split){ //This replicates a bug in the old quantification script in the Dermitzakis lab
                                        	//Check if the start and end positions fall into different bins, as defined in the old script and if so multi count them
                                        	int mul = ( (int) A.starts[i] / 100000 != (int) A.ends[i] / 100000 && (int) overlapping_genes[g].exons[idx].start / 100000  != (int) overlapping_genes[g].exons[idx].end / 100000 ) ?  2 : 1;
                                        	//The buggy way of calculating overlap
                                        	unsigned int ovrlp = overlapping_genes[g].exons[idx].end - A.starts[i] < A.ends[i] - A.starts[i] ? (overlapping_genes[g].exons[idx].end - A.starts[i] + 1) * mul : A.lengths[i] * mul;
                                        	exon_count1.push_back((double) ovrlp);
                                        	exon_overlap1_length_total += ovrlp;
                                        }else {
                                        	unsigned int ovrlp = re - rs + 1; //overlap length
                                        	double total_count = 0.0;
                                        	for (int mbo = 0 ; mbo < A.base_count[i].size(); mbo++) total_count += (double) A.base_count[i][mbo].overlap(rs,re); //Check if any of the bases of the overlap include bases that are covered by both of the mate pairs
                                        	exon_count1.push_back((total_count * 0.5) + ((double) ovrlp - total_count)); //If there are bases covered by both mate pairs then these contribute half of a base that is covered by only one of the mate pairs
                                        	exon_overlap1_length_total += (double) ovrlp;
                                        	//The denominator will be the total length of the read that overlaps with exons, which may be less than the aligned read length, so that they add up to 1
                                        }
                                    }
                                }
                            }
                            if (idx == -1) all_found1 = false; //If one of the block
                        }//mate1 loop

                        //check if the second mate is overlapping any exons
                        for (int i = 0 ; i < B.starts.size(); i++){
                            int idx = -1;
                            if (overlapping_genes[g].overlap(B.starts[i],B.ends[i])){
                                for (int e = 0 ; e < overlapping_genes[g].exons.size(); e++){
                                    if(overlapping_genes[g].exons[e].overlap(B.starts[i],B.ends[i])){
                                        idx = e;
                                        any_found2 = true;
                                        exon_overlap2.push_back(idx);
                                        unsigned int rs = max(overlapping_genes[g].exons[idx].start, B.starts[i]);
                                        unsigned int re = min(overlapping_genes[g].exons[idx].end, B.ends[i]);
                                        if ( filter.old_wrong_split) {
                                        	int mul = ( (int) B.starts[i] / 100000 != (int) B.ends[i] / 100000 && (int) overlapping_genes[g].exons[idx].start / 100000  != (int) overlapping_genes[g].exons[idx].end / 100000 ) ?  2 : 1;
                                        	unsigned int ovrlp = overlapping_genes[g].exons[idx].end - B.starts[i] < B.ends[i] - B.starts[i] ? (overlapping_genes[g].exons[idx].end - B.starts[i] + 1) * mul: B.lengths[i] * mul;
                                        	exon_count2.push_back((double) ovrlp);
                                        	exon_overlap2_length_total += (double) ovrlp;
                                        }else{
                                        	unsigned int ovrlp = re - rs + 1;
                                        	double total_count = 0.0;
                                        	for (int mbo = 0 ; mbo < B.base_count[i].size(); mbo++) total_count += (double) B.base_count[i][mbo].overlap(rs,re);
                                        	exon_count2.push_back((total_count * 0.5) + ((double) ovrlp - total_count));
                                        	exon_overlap2_length_total += (double) ovrlp;
                                        }
                                    }
                                }
                            }
                            if (idx == -1) all_found2 = false;
                        }//mate2 loop
                        if (!all_found1 && filter.check_consistency){ //If we require all blocks of a read to be aligned to an exon
#ifdef DEBUG
                        	cerr << "NCONSALL1\t" << A.name <<endl;
                        	cerr << overlapping_genes[g];
                        	cerr << A << endl;
                        	cerr << B << endl;
#endif
                        	continue;
                        }

                        if (!any_found1){ //If any of the blocks overlap with an exon
#ifdef DEBUG
                        	cerr << "NCONSANY1\t" << A.name<<endl;
                        	cerr << overlapping_genes[g];
                        	cerr << A << endl;
                        	cerr << B << endl;
#endif
                        	continue;
                        }


                        if (!all_found2 && filter.check_consistency){
#ifdef DEBUG
                        	cerr << "NCONSALL2\t" << B.name<<endl;
                        	cerr << overlapping_genes[g];
                        	cerr << A << endl;
                        	cerr << B << endl;
#endif
                        	continue;
                        }
                        if (!any_found2){
#ifdef DEBUG
                        	cerr << "NCONSANY2\t" << B.name<<endl;
                        	cerr << overlapping_genes[g];
                        	cerr << A << endl;
                        	cerr << B << endl;
#endif
                        	continue;
                        }

                        both_found= true;
                        total_genes_found++;
                        float tot_re1 = 0.0, tot_re2 = 0.0;
                        for (int i = 0 ; i < exon_overlap1.size(); i++) {
                        	double ttlc = exon_count1[i] / exon_overlap1_length_total;
                        	genes[genes_map[overlapping_genes[g].gene_id]].exons[exon_overlap1[i]].read_count += ttlc;
                        	tot_re1 += ttlc;
                        	total_read1_count += ttlc;
#ifdef DEBUG
                        	cerr << genes[genes_map[overlapping_genes[g].gene_id]].exons[exon_overlap1[i]].name << "\t" << A.name << "\t" << ttlc << "\t" << exon_count1[i] << "\t" << exon_overlap1_length_total << "\tA" << endl;
#endif
                        }
                        for (int i = 0 ; i < exon_overlap2.size(); i++) {
                        	double ttlc = exon_count2[i] / exon_overlap2_length_total;
                        	genes[genes_map[overlapping_genes[g].gene_id]].exons[exon_overlap2[i]].read_count += ttlc ;
                        	tot_re2 += ttlc;
                        	total_read2_count += ttlc;
#ifdef DEBUG
                        	cerr << genes[genes_map[overlapping_genes[g].gene_id]].exons[exon_overlap2[i]].name << "\t" << B.name << "\t" << ttlc << "\t" << exon_count2[i] << "\t" << exon_overlap2_length_total << "\tB" << endl;
#endif
                        }
                        genes[genes_map[overlapping_genes[g].gene_id]].read_count += tot_re1 + tot_re2;
#ifdef DEBUG
                        cerr << overlapping_genes[g].gene_id << "\t" << A.name << "\t" << tot_re1 <<endl;
                        cerr << overlapping_genes[g].gene_id << "\t" << B.name << "\t" << tot_re2 <<endl;
#endif
                    }//gene loop
                    if(both_found) {
                    	stats.exonicint += 2;
                    	stats.exonic_multi += total_read1_count + total_read2_count;
                    	stats.exonicint_multi += 2 * total_genes_found;
                    }
                    else stats.notexon+=2;
                }else{
                	//non-exonic
#ifdef DEBUG
                	cerr << "NONEXONIC " << A.name << " " << temp << endl;
#endif
                	stats.notexon+=2;
                }
            }else{
            	//First of the mate pairs
            	read_sink[B.name] = B;
            }
    	}else{
    		//Single end

            //Check if the read fails filters
            if (B.filter == UNMAP){
            	stats.unmapped +=1;
            	continue;
            }
            if (B.filter == FAILQC){
            	stats.failqc +=1;
            	continue;
            }
            if (B.filter == DUP){
            	stats.dup +=1;
            	continue;
            }
            if (B.filter == MAPQ){
            	stats.mapQ +=1;
            	continue;
            }
            if (B.filter == NPP){
            	stats.unpaired +=1;
            	continue;
            }
            if (B.filter == MISM){
            	stats.mismatch +=1;
            	continue;
            }

    		if (filter.max_mismatch_count_total >= 0 && filter.fraction_mmt && (double) (B.mmc) / (double) (B.core.l_qseq) > filter.max_mismatch_count_total) {
    			stats.mismatch+=1;
    			continue;
    		}
    		stats.good +=1;

    		//Check we overlap with a gene
    		my_block temp = my_block(B.chr, B.starts[0], B.ends.back());
            vector < my_gene > overlapping_genes = getOverlappingGenesQ(temp);
            float total_read2_count = 0.0;
            int total_genes_found = 0;
    		if (overlapping_genes.size()){
    			//overlapping genes found
    			bool both_found = false;
    			for (int g = 0 ; g < overlapping_genes.size(); g++){
    				vector < int >exon_overlap2;
					vector <float> exon_count2;
    				float exon_overlap2_length_total = 0.0;
    				bool all_found2 = true , any_found2 = false;
                    //check if the second mate is overlapping any exons
                    for (int i = 0 ; i < B.starts.size(); i++){
                        int idx = -1;
                        if (overlapping_genes[g].overlap(B.starts[i],B.ends[i])){
                            for (int e = 0 ; e < overlapping_genes[g].exons.size(); e++){
                                if(overlapping_genes[g].exons[e].overlap(B.starts[i],B.ends[i])){
                                    idx = e;
                                    any_found2 = true;
                                    exon_overlap2.push_back(idx);
                                    unsigned int rs = max(overlapping_genes[g].exons[idx].start, B.starts[i]);
                                    unsigned int re = min(overlapping_genes[g].exons[idx].end, B.ends[i]);
                                    if ( filter.old_wrong_split) {
                                    	int mul = ( (int) B.starts[i] / 100000 != (int) B.ends[i] / 100000 && (int) overlapping_genes[g].exons[idx].start / 100000  != (int) overlapping_genes[g].exons[idx].end / 100000 ) ?  2 : 1;
                                    	unsigned int ovrlp = overlapping_genes[g].exons[idx].end - B.starts[i] < B.ends[i] - B.starts[i] ? (overlapping_genes[g].exons[idx].end - B.starts[i] + 1) * mul: B.lengths[i] * mul;
                                    	exon_count2.push_back((float) ovrlp);
                                    	exon_overlap2_length_total += (float) ovrlp;
                                    }else{
                                    	unsigned int ovrlp = re - rs + 1;
                                    	float total_count = 0.0;
                                    	for (int mbo = 0 ; mbo < B.base_count[i].size(); mbo++) total_count += (float) B.base_count[i][mbo].overlap(rs,re);
                                    	exon_count2.push_back((total_count * 0.5) + ((float) ovrlp - total_count));
                                    	exon_overlap2_length_total += (float) ovrlp;
                                    }
                                }
                            }
                        }
                        if (idx == -1) all_found2 = false;
                    }//mate2 loop
                    if (!all_found2 && filter.check_consistency) continue;
                    if (!any_found2) continue;

                    both_found = true;
                    total_genes_found++;
                    float tot_re2 = 0.0;
                    for (int i = 0 ; i < exon_overlap2.size(); i++) {
                    	float ttlc = exon_count2[i] / exon_overlap2_length_total;
                    	genes[genes_map[overlapping_genes[g].gene_id]].exons[exon_overlap2[i]].read_count += ttlc ;
                    	tot_re2 += ttlc;
                    	total_read2_count += ttlc;
                    }
                    genes[genes_map[overlapping_genes[g].gene_id]].read_count += tot_re2;
    			}
                if(both_found) {
                	stats.exonicint += 1;
                	stats.exonic_multi += total_read2_count;
                	stats.exonicint_multi += total_genes_found;
                }
                else stats.notexon+=1;
    		}else stats.notexon++;
    	}//single end
    }//BAM loop
    vrb.bullet("DONE: " + stb.str(line_count) + " lines read");
    if (read_sink.size() && !region.isSet()) vrb.warning(stb.str(read_sink.size()) + " unmatched mate pairs found");
    if (stats.checkedNM && stats.noNM) vrb.warning("Checked " + stb.str(stats.checkedNM) + " reads for number of mismatches (NM tag) and " + stb.str(stats.noNM) + " reads did NOT have a proper NM tag. These reads are treated as if they had no mismatches and hence included in the quantifications.");
    if (stats.checkedNM && stats.nonIntNM) vrb.warning("Checked " + stb.str(stats.checkedNM) + " reads for number of mismatches (NM tag) and " + stb.str(stats.nonIntNM) + " reads did NOT have an integer NM tag.");
    if (stats.checkedNM && stats.multiNM) vrb.warning("Checked " + stb.str(stats.checkedNM) + " reads for number of mismatches (NM tag) and " + stb.str(stats.multiNM) + " reads had multiple NM tags");
    bam_destroy1(b);
    hts_idx_destroy(idx);
    bam_hdr_destroy(header);
    sam_close(fd);
}

