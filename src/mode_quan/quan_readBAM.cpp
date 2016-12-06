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

int quan_data::read_bam(void *data, bam1_t *b, quan_stats &f, unsigned int &mmc) {
	aux_tq * aux = (aux_tq*) data;
	int ret;
    
	while (1) {
		ret = aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
		if (ret < 0) break;
        f.total++;
        string name = bam_get_qname(b);
		if (f.failed.count(name)){
            if (debug) cerr << f.failed[name] << "\t" << name << endl;
            switch(f.failed[name]){
                case 'u': f.unmapped++; break;
                case 'd': f.dup++; break;
                case 'm': f.mapQ++; break;
                case 'p': f.unpaired++; break;
                case 'M': f.mismatch++; break;
            }
			continue;
		}
        
		if (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY)) {
			f.failed[name] = 'u';
            if (debug) cerr << f.failed[name] << "\t" << name << endl;
			f.unmapped++;
			continue;
		}

		if(aux->fail_qc && b->core.flag & BAM_FQCFAIL){
			f.failed[name] = 'u';
            if (debug) cerr << f.failed[name] << "\t" << name << endl;
			f.unmapped++;
			continue;
		}
        
		if (aux->dup_remove && (b->core.flag & BAM_FDUP)) {
			f.failed[name] = 'd';
            if (debug) cerr << f.failed[name] << "\t" << name << endl;
			f.dup++;
			continue;
		}
        
		if ((int)b->core.qual < aux->min_mapQ) {
			f.failed[name] = 'm';
            if (debug) cerr << f.failed[name] << "\t" << name << endl;
			f.mapQ++;
			continue;
		}
        
		if (b->core.flag & BAM_FPAIRED) {
			if ((b->core.flag & BAM_FREVERSE) == (b->core.flag & BAM_FMREVERSE) || (b->core.flag & BAM_FMUNMAP) ){
				f.failed[name] = 'p';
                if (debug) cerr << f.failed[name] << "\t" << name << endl;
				f.unpaired++;
				continue;
			}
            if (proper_pair && !(b->core.flag & BAM_FPROPER_PAIR)){
                f.failed[name] = 'p';
                if (debug) cerr << f.failed[name] << "\t" << name << endl;
                f.unpaired++;
                continue;
            }
		}
        
		if (aux->max_mismatch_count >= 0 || aux->max_mismatch_count_total >= 0){
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
			if ((!fraction_mm && mc > aux->max_mismatch_count) || (fraction_mm && (double) mc / (double) (b->core.l_qseq) > aux->max_mismatch_count)){
				f.failed[name] = 'M';
                if (debug) cerr << f.failed[name] << "\t" << name << endl;
				f.mismatch++;
				continue;
			}
			mmc = mc;
		}
		break;
    }
    return ret;
}

void quan_data::readBams(){
    vrb.title("Reading BAMs");
    
    //ALLOCATE
    stats = vector < quan_stats >(bams.size());
    for (int gr = 0 ; gr < gene_grps.size(); gr++) gene_grps[gr].allocate(bams.size());
    
    for (int bm = 0 ; bm < bams.size(); bm++){
        string fbam = bams[bm];
        vrb.bullet(fbam + " [" + stb.str(bm+1) + " / " + stb.str(bams.size()) +"]");
        aux_tq * data = (aux_tq *) malloc (sizeof(aux_tq));
        data->min_mapQ = min_mapQ;
        data->dup_remove = dup_remove;
        data->max_mismatch_count = max_mismatch_count;
        data->max_mismatch_count_total = max_mismatch_count_total;
        data->max_intron_length = max_intron_length;
        data->fail_qc = fail_qc;
        data->fp = sam_open(fbam.c_str(), "r");
        if (data->fp == 0) vrb.error("Cannot open file! [" + fbam + "]");
        data->hdr = sam_hdr_read(data->fp);
        if (data->hdr == 0) vrb.error("Cannot parse header![" + fbam + "]");
        hts_idx_t *idx = sam_index_load(data->fp, fbam.c_str());
        if (idx == NULL) vrb.error("Cannot load index![" + fbam + ".bai]");
        data->idx = idx;
        quan_stats stat;
        for (int gr = 0 ; gr < gene_grps.size(); gr++){
            vrb.progress( ((float) gr + 1.0) / (float) gene_grps.size());
            data->iter = sam_itr_querys(data->idx, data->hdr, gene_grps[gr].region.c_str()); // set the iterator
            if (data->iter == NULL) {
            	vrb.warning("Problem jumping to region [" + gene_grps[gr].region + "]");
                hts_itr_destroy(data->iter);
            	continue;
            }
            map < string , quan_block > read_sink;
            bam1_t *b = bam_init1();
            int r;
            unsigned int mmc = 0;
            while((r=read_bam((void*)data,b,stat,mmc)) >= 0){
                string name = bam_get_qname(b);
                const bam1_core_t *c = &b->core;
                if (c->n_cigar) { // cigar
                    quan_block B;
                    B.mmc = mmc;
                    unsigned int bS = b->core.pos+1;
                    unsigned int bL = 0;
                    uint32_t *cigar = bam_get_cigar(b);
                    for (int i = 0; i < c->n_cigar; ++i) {
                        int l = bam_cigar_oplen(cigar[i]);
                        char c = bam_cigar_opchr(cigar[i]);
                        if(c=='S' || c =='H' || c=='P') continue;
                        else if (c=='N' && l){
                            B.starts.push_back(bS);
                            B.ends.push_back(bS+bL-1);
                            B.lengths.push_back(bL);
                            B.block_overlap.push_back(1.0);
                            B.read_length+=bL;
                            bS+=bL+l;
                            bL = 0;
                        }else bL += l;
                    }
                    B.starts.push_back(bS);
                    B.ends.push_back(bS+bL-1);
                    B.lengths.push_back(bL);
                    B.block_overlap.push_back(1.0);
                    B.read_length+=bL;
                    B.core = b->core;
                    if (b->core.flag & BAM_FPAIRED){
                        if(read_sink.count(name)){
                            quan_block A = read_sink[name];
                            read_sink.erase(name);
                            if (max_mismatch_count_total >= 0 && ((!fraction_mmt && A.mmc + B.mmc > max_mismatch_count_total) || (fraction_mmt && (double) (A.mmc + B.mmc) / (double) (A.core.l_qseq + B.core.l_qseq) > max_mismatch_count_total) )){
                                if (debug){
                                    cerr << "M\t" << name<<endl;
                                    cerr << "M\t" << name<<endl;
                                }
                            	stat.mismatch+=2;
                            	continue;
                            }
                            if(A.core.mtid != B.core.tid || (A.core.flag & BAM_FREVERSE) || !(B.core.flag & BAM_FREVERSE) || ( old_wrong_split && B.core.pos <= A.core.pos) || ( !old_wrong_split && B.core.pos < A.core.pos)){
                                stat.failed[name] = 'p';
                                if (debug) cerr << stat.failed[name] << "\t" << name << endl;
                                stat.unpaired++;
                                continue;
                            }

                            if (merge) A.merge(B);
                            stat.good +=2;
                            bool both_found = false;
                            for (int g = 0 ; g < gene_grps[gr].genes.size(); g++){
                                vector < int > exon_overlap1,exon_overlap2,exon_overlap1_length,exon_overlap2_length,exon_map1,exon_map2;
                                bool all_found1 = true, all_found2 = true , any_found1 = false, any_found2 = false;
                                unsigned long int exon_overlap1_length_total = 0 , exon_overlap2_length_total = 0;
                                for (int i = 0 ; i < A.starts.size() ; i++){
                                    int idx = -1;
                                    if (gene_grps[gr].genes[g].overlap(A.starts[i],A.ends[i])){
                                        for (int e = 0 ; e < gene_grps[gr].genes[g].exons.size();e++){
                                            if(gene_grps[gr].genes[g].exons[e].overlap(A.starts[i],A.ends[i])){
                                                idx = e;
                                                any_found1 = true;
                                                exon_overlap1.push_back(idx);
                                                if ( old_wrong_split) exon_overlap1_length.push_back(gene_grps[gr].genes[g].exons[idx].end - A.starts[i] < A.ends[i] - A.starts[i] ? gene_grps[gr].genes[g].exons[idx].end - A.starts[i] + 1 : A.lengths[i]);
                                                else exon_overlap1_length.push_back(min(gene_grps[gr].genes[g].exons[idx].end, A.ends[i]) - max(gene_grps[gr].genes[g].exons[idx].start, A.starts[i]) + 1);
                                                exon_overlap1_length_total += exon_overlap1_length.back();
                                                exon_map1.push_back(i);
                                            }
                                        }
                                    }
                                    if (idx == -1) all_found1 = false;
                                }


                                for (int i = 0 ; i < B.starts.size(); i++){
                                    int idx = -1;
                                    if (gene_grps[gr].genes[g].overlap(B.starts[i],B.ends[i])){
                                        for (int e = 0 ; e < gene_grps[gr].genes[g].exons.size(); e++){
                                            if(gene_grps[gr].genes[g].exons[e].overlap(B.starts[i],B.ends[i])){
                                                idx = e;
                                                any_found2 = true;
                                                exon_overlap2.push_back(idx);
                                                if ( old_wrong_split) exon_overlap2_length.push_back(gene_grps[gr].genes[g].exons[idx].end - B.starts[i] < B.ends[i] - B.starts[i] ? gene_grps[gr].genes[g].exons[idx].end - B.starts[i] + 1 : B.lengths[i]);
                                                else exon_overlap2_length.push_back(min(gene_grps[gr].genes[g].exons[idx].end, B.ends[i]) - max(gene_grps[gr].genes[g].exons[idx].start, B.starts[i]) + 1);
                                                exon_overlap2_length_total += exon_overlap2_length.back();
                                                exon_map2.push_back(i);
                                            }
                                        }
                                    }
                                    if (idx == -1) all_found2 = false;
                                }

                                if (!all_found1 && check_consistency){
                                    if (debug){
                                        cerr << "NCONSALL1\t" << name<<endl;
                                        cerr << gene_grps[gr].genes[g];
                                        cerr << A;
                                        cerr << B;
                                    }
                                    continue;
                                }
                                //for debugging
                                if (!any_found1){
                                    if (debug){
                                        cerr << "NCONSANY1\t" << name<<endl;
                                        cerr << gene_grps[gr].genes[g];
                                        cerr << A;
                                        cerr << B;
                                    }
                                    continue;
                                }
                                
                                if (!all_found2 && check_consistency){
                                    if (debug){
                                        cerr << "NCONSALL2\t" << name<<endl;
                                        cerr << gene_grps[gr].genes[g];
                                        cerr << A;
                                        cerr << B;
                                    }
                                    continue;
                                }
                                if (!any_found2){
                                    if(debug){
                                        cerr << "NCONSANY2\t" << name<<endl;
                                        cerr << gene_grps[gr].genes[g];
                                        cerr << A;
                                        cerr << B;
                                    }
                                    continue;
                                }
          
                                
                                both_found= true;
                                for (int i = 0 ; i < exon_overlap1.size(); i++) {
                                    gene_grps[gr].genes[g].exons[exon_overlap1[i]].read_count[bm] += (double)exon_overlap1_length[i] / (double)exon_overlap1_length_total * A.block_overlap[exon_map1[i]];
                                    if (debug) cerr << gene_grps[gr].genes[g].exons[exon_overlap1[i]].name << "\t" << name << "\t" <<(double)exon_overlap1_length[i] / (double)exon_overlap1_length_total * A.block_overlap[exon_map1[i]] << endl;
                                }
                                for (int i = 0 ; i < exon_overlap2.size(); i++) {
                                    gene_grps[gr].genes[g].exons[exon_overlap2[i]].read_count[bm] += (double)exon_overlap2_length[i] / (double)exon_overlap2_length_total * B.block_overlap[exon_map2[i]];
                                    if (debug) cerr << gene_grps[gr].genes[g].exons[exon_overlap2[i]].name << "\t" << name << "\t" <<(double)exon_overlap2_length[i] / (double)exon_overlap2_length_total * B.block_overlap[exon_map2[i]] << endl;
                                }
                                gene_grps[gr].genes[g].read_count[bm]+= A.total_contribution + B.total_contribution;
                                if (debug){
                                    cerr << gene_grps[gr].genes[g].ID << "\t" << name << "\t" << A.total_contribution <<endl;
                                    cerr << gene_grps[gr].genes[g].ID << "\t" << name << "\t" << B.total_contribution <<endl;
                                }

                            }//gene loop
                            if(both_found) stat.exonic+= A.total_contribution + B.total_contribution;
                            else stat.notexon+=2;
                        }else read_sink[name]= B;
                    }else{
                    	if (max_mismatch_count_total >= 0 && fraction_mmt && (double) (B.mmc) / (double) (B.core.l_qseq) > max_mismatch_count_total) {
                    		stat.mismatch+=1;
                    		continue;
                    	}
                        stat.good +=1;
                        bool both_found = false;
                        for (int g = 0 ; g < gene_grps[gr].genes.size(); g++){
                            vector < int >exon_overlap2,exon_overlap2_length,exon_map2;
                            int exon_overlap2_length_total = 0;
                            bool all_found2 = true , any_found2 = false;
                            for (int i = 0 ; i < B.starts.size(); i++){
                                int idx = -1;
                                if (gene_grps[gr].genes[g].overlap(B.starts[i],B.ends[i])){
                                    for (int e = 0 ; e < gene_grps[gr].genes[g].exons.size(); e++){
                                        if(gene_grps[gr].genes[g].exons[e].overlap(B.starts[i],B.ends[i])){
                                            idx = e;
                                            any_found2 = true;
                                            exon_overlap2.push_back(idx);
                                            exon_overlap2_length.push_back(min(gene_grps[gr].genes[g].exons[idx].end, B.ends[i]) - max(gene_grps[gr].genes[g].exons[idx].start, B.starts[i]) + 1);
                                            exon_overlap2_length_total += exon_overlap2_length[i];
                                            exon_map2.push_back(i);
                                        }
                                    }
                                }
                                if (idx == -1) all_found2 = false;
                            }
                            if (!all_found2 && check_consistency) continue;
                            if (!any_found2) continue;
                            
                            both_found = true;

                            for (int i = 0 ; i < exon_overlap2.size(); i++) {
                                gene_grps[gr].genes[g].exons[exon_overlap2[i]].read_count[bm] += (double)exon_overlap2_length[i] / (double)exon_overlap2_length_total * B.block_overlap[exon_map2[i]];
                            }
                            gene_grps[gr].genes[g].read_count[bm] += B.total_contribution;
                        }//gene loop
                        if(both_found) stat.exonic += B.total_contribution;
                        else stat.notexon+=1;
                    }
                }//cigar
            }//while loop
            bam_destroy1(b);
            stat.notexon += read_sink.size(); //Orphan mate pairs will not map to the same gene
            if (debug) for (map < string , quan_block >::iterator it = read_sink.begin() ; it != read_sink.end(); it++) cerr << "O\t" << it->first << endl;
            read_sink.clear();
            hts_itr_destroy(data->iter);
        } //gene_grp floop
        stat.failed.clear();
        stats[bm] = stat;
        bam_hdr_destroy(data->hdr);
        hts_idx_destroy(data->idx);
        if (data->fp) sam_close(data->fp);
        free(data);
        
    }//bam loop
    
}

