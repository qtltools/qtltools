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

#ifndef _QUAN_DATA_H
#define _QUAN_DATA_H

#include "../common/data.h"


class quan_stats{
public:
	map < string , char > failed;
	unsigned long int mapQ,mismatch,unmapped,unpaired,dup,notexon,good,total,exonicint;
    double exonic;
    quan_stats(){mapQ=0;mismatch=0;unmapped=0;dup=0;unpaired=0;notexon=0;exonic=0.0;total=0;good=0;exonicint=0;}
};

class quan_exon{
public:
    unsigned int start;
    unsigned int end;
    unsigned int length;
    string chr;
    string name;
    string gene_id;
    string gene_name;
    string gene_type;
    short int strand;
    vector < double > read_count;
    bool merged;
    quan_exon(){chr="NA";name="NA";start=0;end=0;gene_id="NA";gene_name="NA";gene_type="NA";strand = 0;;merged=false;length=0;}
    ~quan_exon(){read_count.clear();}
    quan_exon(string c , unsigned int s, unsigned int e, string gi, string gn, string gt, string st, unsigned int nof){chr=c;start=s;end=e;gene_id=gi;name = gi + "_" + stb.str(start) + "_" + stb.str(end);gene_name = gn; gene_type =gt; strand = st == "-" ? -1 : 1; length = end - start + 1; merged=false; read_count=vector < double >(nof,0.0);}
    quan_exon(string c , unsigned int s, unsigned int e, string gi, string gn, string gt, string st){chr=c;start=s;end=e;gene_id=gi;name = gi + "_" + stb.str(start) + "_" + stb.str(end);gene_name = gn; gene_type =gt; strand = st == "-" ? -1 : 1; length = end - start + 1; merged=false;}
    void allocate(unsigned int nof) { read_count=vector < double >(nof,0.0); }
    
    void resize(unsigned int s, unsigned int e){start = s; end = e;length = end - start + 1; name = gene_id + "_" + stb.str(start) + "_" + stb.str(end);merged=true;}
    void merge(vector < quan_exon > &exs){
        unsigned int sr = start;
        unsigned int er = end;
        for (int i =0 ; i < exs.size(); i++){
            if (exs[i].start < sr) sr = exs[i].start;
            if (exs[i].end > er) er = exs[i].end;
        }
        resize(sr,er);
    }
    
    bool overlap(unsigned int s, unsigned int e){
    	if (s <= end && e >= start) return true;
    	else return false;
    }

    bool overlap(string c, unsigned int s, unsigned int e){
    	if (c == chr && s <= end && e >= start) return true;
    	else return false;
    }

    bool operator < (quan_exon const & p) const {
        if (chr < p.chr) return true;
        if (chr > p.chr) return false;
        if (start < p.start) return true;
        if (start >= p.start) return false;
        return false;
    }
    
    friend ostream& operator<<(ostream& out, const quan_exon& p){
        out << p.chr << "\t" << p.start << "\t" << p.end << "\t" << p.name << "\t" << p.length << "\t" << (p.strand == 1 ? "+" : "-") << "\t"  << p.merged << endl;
        return out;
    }
};

class quan_gene{
public:
    string ID;
    string name;
    string chr;
    string region;
    unsigned int start;
    unsigned int end;
    unsigned int tss;
    short int strand;
    unsigned int length;
    vector < double > read_count;
    vector < quan_exon > exons;
    quan_gene(){chr="NA";name="NA";start=0;end=0;ID="NA";strand = 0;region="NA"; length = 0 ;tss=0;};
    quan_gene(unsigned int nof){chr="NA";name="NA";start=0;end=0;ID="NA";strand = 0;region="NA"; length = 0 ;tss=0;read_count=vector < double >(nof,0.0);};
    ~quan_gene(){exons.clear();read_count.clear();}
    void allocate(unsigned int nof) {
        read_count=vector < double >(nof,0.0);
        for (int i = 0 ; i < exons.size(); i++) exons[i].allocate(nof);
    }
    void assign(quan_exon &e){
        if (strand && e.strand != strand) vrb.error("Strand mismatch");
        if (chr != "NA" && e.chr != chr) vrb.error("Chr mismatch");
        if (name != "NA" && e.gene_name != name) vrb.error("Name mismatch");
        if (ID != "NA" && e.gene_id != ID) vrb.error("ID mismatch");
        strand = e.strand;
        name = e.gene_name;
        chr = e.chr;
        ID = e.gene_id;
        if (start == 0 || e.start < start ) start = e.start;
        if (end == 0 || e.end > end) end = e.end;
        tss = strand == -1 ? end : start;
        vector < quan_exon > overlaping;
        vector < int > overlaping_idx;
        unsigned int tot = 0;
        sort(exons.begin(), exons.end());
        for (int i =0 ; i < exons.size() ; i++){
            if (exons[i].end >= e.start && exons[i].start <= e.end ){
                overlaping.push_back(exons[i]);
                overlaping_idx.push_back(i);
                tot += exons[i].length;
            }
        }
        if (overlaping.size()){
            e.merge(overlaping);
            exons.erase(exons.begin() + overlaping_idx[0], exons.begin() + (overlaping_idx.back() + 1));
            length -= tot;
        }
        exons.push_back(e);
        length += e.length;
        region = chr + ":" + stb.str(start) + "-" + stb.str(end);
    }
    
    bool operator < (quan_gene const & p) const {
        if (chr < p.chr) return true;
        if (chr > p.chr) return false;
        if (tss < p.tss) return true;
        if (tss >= p.tss) return false;
        return false;
    }
    
    bool overlap(unsigned int s, unsigned int e){
    	if (s <= end && e >= start) return true;
    	else return false;
    }

    bool overlap(string c, unsigned int s, unsigned int e){
    	if (c == chr && s <= end && e >= start) return true;
    	else return false;
    }
    
    bool overlap(genomic_region &gr){
        if (gr.chr == chr && gr.start <= end && gr.end >= start) return true;
        else return false;
    }

    friend ostream& operator<<(ostream& out, quan_gene& g){
        out << g.chr << "\t" << g.start << "\t" << g.end << "\t" << g.ID << "\t" << g.length << "\t" << (g.strand == 1 ? "+" : "-") << "\t" << 2<<endl;
        sort(g.exons.begin(), g.exons.end());
        for (int i =0; i < g.exons.size(); i++ ) out << g.exons[i];
        return out;
    }
    
};

class quan_gene_grp{
public:
    vector < quan_gene > genes;
    string chr;
    unsigned int start;
    unsigned int end;
    string region;
    quan_gene_grp(){start=0;end=0;region="";chr="";}
    ~quan_gene_grp(){genes.clear();}
    void allocate(unsigned int nof){ for(int i = 0 ; i < genes.size(); i++) genes[i].allocate(nof);}
    bool overlap(unsigned int s, unsigned int e){
        if (s <= end && e >= start) return true;
        else return false;
    }
    
    bool overlap(string c, unsigned int s, unsigned int e){
        if (c == chr && s <= end && e >= start) return true;
        else return false;
    }
    
    bool overlap(genomic_region &gr){
        if (gr.chr == chr && gr.start <= end && gr.end >= start) return true;
        else return false;
    }
};

typedef struct {     			// auxiliary data structure
    samFile * fp;				// the file handle
    bam_hdr_t * hdr;			// the file header
    hts_itr_t * iter;			// NULL if a region not specified
    hts_idx_t * idx;
    int min_mapQ;				// mapQ filter
    bool dup_remove,fail_qc;	// remove duplicates, failed QC reads
    int max_intron_length;
    double max_mismatch_count;
    double max_mismatch_count_total;
} aux_tq;

class quan_block{
public:
    vector < unsigned int > starts;
    vector < unsigned int > ends;
    vector < unsigned int > lengths;
    vector < double > block_overlap;
    double total_contribution;
    unsigned int read_length;
    unsigned int mmc;
    bam1_core_t core;
    quan_block(){read_length=0;mmc=0; total_contribution = 1.0;}
    void merge(quan_block &B){
        for (int i =0 ; i < starts.size(); i++){
            for (int j =i; j < B.starts.size(); j++){
                if(starts[i] <= B.ends[j] && ends[i] >= B.starts[j]){
                    unsigned int overlap = min(ends[i], B.ends[j]) - max(starts[i], B.starts[j]) + 1;
                    block_overlap[i] = 0.5 + ((lengths[i] - overlap) / lengths[i] * 0.5);
                    B.block_overlap[j] = 0.5 + ((B.lengths[j] - overlap) / B.lengths[j] * 0.5);
                    break;
                }
            }
        }
        double total = 0.0;
        for (int i = 0 ; i < block_overlap.size(); i++) total += (double) lengths[i] * block_overlap[i];
        total_contribution = total / (double) read_length;
        total = 0.0;
        for (int i = 0 ; i < B.block_overlap.size(); i++) total += (double) B.lengths[i] * B.block_overlap[i];
        B.total_contribution = total / (double) B.read_length;
    }
    friend ostream& operator<<(ostream& out, quan_block& g){
        out << g.read_length << "\t" << g.mmc << "\t" << g.core.pos << "\t" << g.core.mpos << "\t" << g.total_contribution;
        for (int i =0; i < g.starts.size(); i++ ) out << "\t" << g.starts[i] << "," << g.ends[i] << " " << g.block_overlap[i];
        out << endl;
        return out;
    }
};

class quan_data : public data {
public:
    map < string, unsigned int > genes_map;
    vector < quan_gene > genes;
    vector < quan_gene_grp > gene_grps;
    vector < string > bams;
    vector < quan_stats > stats;
    vector < string > samples;
    set < string > gene_types;
    
    //FILTERS
    int max_intron_length;
    double max_mismatch_count,max_mismatch_count_total;
    unsigned int min_mapQ,max_read_length,min_exon;
    bool dup_remove,proper_pair,check_consistency,debug,merge,fraction_mm,fraction_mmt,fail_qc,old_wrong_split;
    quan_data(){max_intron_length=0;max_mismatch_count=0.0;max_mismatch_count_total=0.0;min_mapQ=0;max_read_length=0;dup_remove=false;proper_pair=false;check_consistency=false;debug=false;merge=true;fraction_mm=false;fraction_mmt=false;fail_qc=false;old_wrong_split=false;}
    
    void setChunk(int,int);
    void setRegion(string);
    genomic_region region;
    
    
    void read_Sample_Names(vector < string > &);
    void readGTF(string,unsigned int);
    void groupGenes();
    void readBams();
    int read_bam(void *, bam1_t *, quan_stats &, unsigned int &mmc);
    void printBEDcount(string);
    void printBEDrpkm(string);
    void printStats(string);
};

void quan_main(vector < string > & );

#endif /* quan_data_h */
