/*
 * quan2_readBAM.cpp
 *
 *  Created on: Aug 21, 2017
 *      Author: halit
 */

#include "quan2_data.h"


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
		uint8_t *s = bam_get_aux(b);
		while (s+4 <= b->data + b->l_data) {
			uint8_t type, key[3];
			key[0] = s[0]; key[1] = s[1]; key[2] = '\0';
			string keys((const char*)key);
			s += 2; type = *s++;
			if (type == 'A') {
				if (keys=="NM") mc = *s;
				++s;
			} else if (type == 'C') {
				if (keys=="NM") mc = *s;
				++s;
			} else if (type == 'c') {
				if (keys=="NM") mc = *(int8_t*)s;
				++s;
			} else if (type == 'S') {
				if (s+2 <= b->data + b->l_data) {
					if (keys=="NM") mc = *(uint16_t*)s;
					s += 2;
				} else break;
			} else if (type == 's') {
				if (s+2 <= b->data + b->l_data) {
					if (keys=="NM") mc = *(int16_t*)s;
					s += 2;
				} else break;
			} else if (type == 'I') {
				if (s+4 <= b->data + b->l_data) {
					if (keys=="NM") mc = *(uint32_t*)s;
					s += 4;
				} else break;
			} else if (type == 'i') {
				if (s+4 <= b->data + b->l_data) {
					if (keys=="NM") mc = *(int32_t*)s;
					s += 4;
				} else break;
			} else if (type == 'f') {
				if (s+4 <= b->data + b->l_data) {
					if (keys=="NM") mc = *(float*)s;
					s += 4;
				} else break;

			} else if (type == 'd') {
				if (s+8 <= b->data + b->l_data) {
					if (keys=="NM") mc = *(double*)s;
					s += 8;
				}else break;
			} else if (type == 'Z' || type == 'H') {
				while (s < b->data + b->l_data && *s) s++;
				if (s >= b->data + b->l_data)
					break;
				++s;
			} else if (type == 'B') {
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
        cerr<< block.name << "\t" << b->core.flag << "\t";
#endif
        for (int i = 0; i < c->n_cigar; ++i) {
            int l = bam_cigar_oplen(cigar[i]);
            char c = bam_cigar_opchr(cigar[i]);
#ifdef DEBUG
            cerr << l << c;
#endif
            if(c=='S' || c =='H' || c =='I' || c=='P') continue;
            else if ((c=='N' || c=='D') && l){
                if (bL){
					block.starts.push_back(bS);
					block.ends.push_back(bS+bL-1);
					block.lengths.push_back(bL);
					block.block_overlap.push_back(1.0);
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
			block.block_overlap.push_back(1.0);
			block.read_length+=bL;
        }
    }
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
                if(A.core.mtid != B.core.tid || (A.core.flag & BAM_FREVERSE) || !(B.core.flag & BAM_FREVERSE) || ( filter.old_wrong_split && B.core.pos <= A.core.pos) || ( !filter.old_wrong_split && B.core.pos < A.core.pos)){
#ifdef DEBUG
                	cerr << "NPPP\t" << B.name << "\tc:" << (A.core.mtid != B.core.mtid) << "\t1s:" << (A.core.flag & BAM_FREVERSE) << "\t2s:" << !(B.core.flag & BAM_FREVERSE) << "\to:" << ( filter.old_wrong_split && B.core.pos <= A.core.pos) << "\tn:" << ( !filter.old_wrong_split && B.core.pos < A.core.pos) << "\t" << A.core.flag << "\t" << B.core.flag << endl;
#endif
                	stats.unpaired++;
                    continue;
                }
                //Passed all filters, good on the reads
                stats.good+=2;
                if (filter.merge) A.merge(B); //if we are merging overlapping mate pair do that
                //Check we overlap with a gene
                my_block temp = my_block(A.chr, A.starts[0], B.ends.back());
                vector < my_gene > overlapping_genes = getOverlappingGenesQ(temp);
                if (overlapping_genes.size()){
                	//overlapping genes found
                	bool both_found = false;
                	int total_genes_found = 0;
                    for (int g = 0 ; g < overlapping_genes.size(); g++){
#ifdef DEBUG
                    	cerr << "CHECKING " << A << " " << B << " in " << temp << " for " << overlapping_genes[g].gene_id <<" " << g+1 << "/" << overlapping_genes.size()<< endl;
#endif
                        vector < int > exon_overlap1,exon_overlap2,exon_overlap1_length,exon_overlap2_length,exon_map1,exon_map2;
                        bool all_found1 = true, all_found2 = true , any_found1 = false, any_found2 = false;
                        unsigned long int exon_overlap1_length_total = 0 , exon_overlap2_length_total = 0;
                        //check if the first mate is overlapping any exons
                        for (int i = 0 ; i < A.starts.size() ; i++){
                            int idx = -1;
                            if (overlapping_genes[g].overlap(A.starts[i],A.ends[i])){
                                for (int e = 0 ; e < overlapping_genes[g].exons.size();e++){
                                    if(overlapping_genes[g].exons[e].overlap(A.starts[i],A.ends[i])){
                                        idx = e;
                                        any_found1 = true;
                                        exon_overlap1.push_back(idx);
                                        if ( filter.old_wrong_split){
                                        	int mul = ( (int) A.starts[i] / 100000 != (int) A.ends[i] / 100000 && (int) overlapping_genes[g].exons[idx].start / 100000  != (int) overlapping_genes[g].exons[idx].end / 100000 ) ?  2 : 1;
                                        	exon_overlap1_length.push_back(overlapping_genes[g].exons[idx].end - A.starts[i] < A.ends[i] - A.starts[i] ? (overlapping_genes[g].exons[idx].end - A.starts[i] + 1) * mul : A.lengths[i] * mul);
                                        }else exon_overlap1_length.push_back(min(overlapping_genes[g].exons[idx].end, A.ends[i]) - max(overlapping_genes[g].exons[idx].start, A.starts[i]) + 1);
                                        exon_overlap1_length_total += exon_overlap1_length.back();
                                        exon_map1.push_back(i);
                                    }
                                }
                            }
                            if (idx == -1) all_found1 = false;
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
                                        if ( filter.old_wrong_split) {
                                        	int mul = ( (int) B.starts[i] / 100000 != (int) B.ends[i] / 100000 && (int) overlapping_genes[g].exons[idx].start / 100000  != (int) overlapping_genes[g].exons[idx].end / 100000 ) ?  2 : 1;
                                        	exon_overlap2_length.push_back(overlapping_genes[g].exons[idx].end - B.starts[i] < B.ends[i] - B.starts[i] ? (overlapping_genes[g].exons[idx].end - B.starts[i] + 1) * mul: B.lengths[i] * mul);
                                        }else exon_overlap2_length.push_back(min(overlapping_genes[g].exons[idx].end, B.ends[i]) - max(overlapping_genes[g].exons[idx].start, B.starts[i]) + 1);
                                        exon_overlap2_length_total += exon_overlap2_length.back();
                                        exon_map2.push_back(i);
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

                        if (!any_found1){
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
                        for (int i = 0 ; i < exon_overlap1.size(); i++) {
                        	genes[genes_map[overlapping_genes[g].gene_id]].exons[exon_overlap1[i]].read_count += (double)exon_overlap1_length[i] / (double)exon_overlap1_length_total * A.block_overlap[exon_map1[i]];
#ifdef DEBUG
                        	cerr << genes[genes_map[overlapping_genes[g].gene_id]].exons[exon_overlap1[i]].name << "\t" << A.name << "\t" <<(double)exon_overlap1_length[i] / (double)exon_overlap1_length_total * A.block_overlap[exon_map1[i]] << "\t" << exon_overlap1_length[i] << "\t" << exon_overlap1_length_total << "\t" << A.block_overlap[exon_map1[i]] << "\t" << A.starts[exon_map1[i]] << "\t" << A.ends[exon_map1[i]] << "\tA" << endl;
#endif
                        }
                        for (int i = 0 ; i < exon_overlap2.size(); i++) {
                        	genes[genes_map[overlapping_genes[g].gene_id]].exons[exon_overlap2[i]].read_count += (double)exon_overlap2_length[i] / (double)exon_overlap2_length_total * B.block_overlap[exon_map2[i]];
#ifdef DEBUG
                        	cerr << genes[genes_map[overlapping_genes[g].gene_id]].exons[exon_overlap2[i]].name << "\t" << B.name << "\t" <<(double)exon_overlap2_length[i] / (double)exon_overlap2_length_total * B.block_overlap[exon_map2[i]] << "\t" << exon_overlap2_length[i] << "\t" << exon_overlap2_length_total << "\t" << B.block_overlap[exon_map2[i]] << "\t" << B.starts[exon_map2[i]] << "\t" << B.ends[exon_map2[i]] << "\tB" << endl;
#endif
                        }
                        genes[genes_map[overlapping_genes[g].gene_id]].read_count += A.total_contribution + B.total_contribution;
#ifdef DEBUG
                        cerr << overlapping_genes[g].gene_id << "\t" << A.name << "\t" << A.total_contribution <<endl;
                        cerr << overlapping_genes[g].gene_id << "\t" << B.name << "\t" << B.total_contribution <<endl;
#endif
                    }//gene loop
                    if(both_found) {
                    	stats.exonic += A.total_contribution + B.total_contribution;
                    	stats.exonicint += 2;
                    	stats.exonic_multi += (A.total_contribution + B.total_contribution) * total_genes_found;
                    	stats.exonicint_multi += 2 * total_genes_found;
                    }
                    else stats.notexon+=2;
                }else{
                	//non-exonic
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
    		if (overlapping_genes.size()){
    			//overlapping genes found
    			bool both_found = false;
    			for (int g = 0 ; g < overlapping_genes.size(); g++){
    				vector < int >exon_overlap2,exon_overlap2_length,exon_map2;
    				int exon_overlap2_length_total = 0;
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
    								if ( filter.old_wrong_split) {
    									int mul = ( (int) B.starts[i] / 100000 != (int) B.ends[i] / 100000 && (int) overlapping_genes[g].exons[idx].start / 100000  != (int) overlapping_genes[g].exons[idx].end / 100000 ) ?  2 : 1;
    									exon_overlap2_length.push_back(overlapping_genes[g].exons[idx].end - B.starts[i] < B.ends[i] - B.starts[i] ? (overlapping_genes[g].exons[idx].end - B.starts[i] + 1) * mul: B.lengths[i] * mul);
    								}else exon_overlap2_length.push_back(min(overlapping_genes[g].exons[idx].end, B.ends[i]) - max(overlapping_genes[g].exons[idx].start, B.starts[i]) + 1);
    								exon_overlap2_length_total += exon_overlap2_length.back();
    								exon_map2.push_back(i);
    							}
    						}
    					}
    					if (idx == -1) all_found2 = false;
    				}//mate2 loop
                    if (!all_found2 && filter.check_consistency) continue;
                    if (!any_found2) continue;

                    both_found = true;

                    for (int i = 0 ; i < exon_overlap2.size(); i++) {
                    	genes[genes_map[overlapping_genes[g].gene_id]].exons[exon_overlap2[i]].read_count += (double)exon_overlap2_length[i] / (double)exon_overlap2_length_total * B.block_overlap[exon_map2[i]];
                    }
                    genes[genes_map[overlapping_genes[g].gene_id]].read_count += B.total_contribution;
    			}
                if(both_found) {stats.exonic += B.total_contribution; stats.exonicint++;}
                else stats.notexon+=1;
    		}else stats.notexon++;
    	}//single end
    }//BAM loop
    vrb.bullet("DONE: " + stb.str(line_count) + " lines read");
    if (read_sink.size() && !region.isSet()) vrb.warning(stb.str(read_sink.size()) + " unmatched mate pairs found");
    bam_destroy1(b);
    hts_idx_destroy(idx);
    bam_hdr_destroy(header);
    sam_close(fd);
}

