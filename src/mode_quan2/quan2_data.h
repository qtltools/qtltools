/*
 * quan2_data.h
 *
 *  Created on: Aug 21, 2017
 *      Author: halit
 */

#ifndef SRC_MODE_QUAN2_QUAN2_DATA_H_
#define SRC_MODE_QUAN2_QUAN2_DATA_H_

#include "../common/data.h"

enum FILTER { PASS, UNMAP, SECD, FAILQC , DUP, NPP, MAPQ, STRAND, MISM};

class my_stats{
public:
	unsigned long int mapQ,mismatch,unmapped,unpaired,dup,notexon,good,total,exonicint,secondary,failqc;
    double exonic;
    my_stats(){mapQ=0;mismatch=0;unmapped=0;dup=0;unpaired=0;notexon=0;exonic=0.0;total=0;good=0;exonicint=0;secondary=0;failqc=0;}
};


class my_cont_blocks{
public:
	string chr;
    vector < unsigned int > starts;
    vector < unsigned int > ends;
    vector < unsigned int > lengths;
    vector < double > block_overlap;
    double total_contribution;
    unsigned int read_length;
    unsigned int mmc;
    FILTER filter;
    bam1_core_t core;
    string name;
    my_cont_blocks(){read_length=0;mmc=0; total_contribution = 1.0;filter=PASS;name="";chr="NA";}
    void merge(my_cont_blocks &B){
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
    friend ostream& operator<<(ostream& out, my_cont_blocks& g){
        out << g.read_length << "\t" << g.mmc << "\t" << g.core.pos << "\t" << g.core.mpos << "\t" << g.total_contribution << g.chr;
        for (int i =0; i < g.starts.size(); i++ ) out << "\t" << g.starts[i] << "," << g.ends[i] << " " << g.block_overlap[i];
        out << endl;
        return out;
    }
};

class filters{
public:
	int max_intron_length;
    double max_mismatch_count,max_mismatch_count_total;
    unsigned int min_mapQ,max_read_length,min_exon;
    bool dup_remove,proper_pair,check_consistency,debug,merge,fraction_mm,fraction_mmt,fail_qc,old_wrong_split;
    filters(){
    	max_intron_length=0;
    	max_mismatch_count=0.0;
    	max_mismatch_count_total=0.0;
    	min_mapQ=0;max_read_length=0;
    	dup_remove=false;
    	proper_pair=false;
    	check_consistency=false;
    	debug=false;
    	merge=true;
    	fraction_mm=false;
    	fraction_mmt=false;
    	fail_qc=false;
    	old_wrong_split=false;
    	min_exon=0;
    }
};

class my_block{
public:
	string chr,gene_name,gene_id,name;
    unsigned int start;
    unsigned int end;
    unsigned int length;
    short int strand;
    double read_count;
    my_block(){start=0,end=0,length=0,strand=0;read_count=0.0;chr="NA",gene_name="NA";gene_id="NA",name="NA";}
    my_block(string c , unsigned int s, unsigned int e){start=s,end=e,length=0,strand=0;read_count=0.0;chr=c,gene_name="NA";gene_id="NA",name="NA";}

    bool operator < (my_block const & p) const {
        if (chr < p.chr) return true;
        if (chr > p.chr) return false;
        if (start < p.start) return true;
        if (start > p.start) return false;
        if (end < p.end) return true;
        return false;
    }

};

class my_exon: public my_block{
public:
    string gene_type;
    bool merged;
    my_exon(){chr="NA";name="NA";start=0;end=0;gene_id="NA";gene_name="NA";gene_type="NA";strand = 0;;merged=false;length=0;read_count=0.0;}
    my_exon(string c , unsigned int s, unsigned int e, string gi, string gn, string gt, string st){chr=c;start=s;end=e;gene_id=gi;name = gi + "_" + stb.str(start) + "_" + stb.str(end);gene_name = gn; gene_type =gt; strand = st == "-" ? -1 : 1; length = end - start + 1; merged=false;read_count=0.0;}

    void resize(unsigned int s, unsigned int e){start = s; end = e;length = end - start + 1; name = gene_id + "_" + stb.str(start) + "_" + stb.str(end);merged=true;}
    void merge(vector < my_exon > &exs){
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

    bool operator < (my_exon const & p) const {
        if (chr < p.chr) return true;
        if (chr > p.chr) return false;
        if (start < p.start) return true;
        if (start > p.start) return false;
        if (end < p.end) return true;
        return false;
    }

    friend ostream& operator<<(ostream& out, const my_exon& p){
        out << p.chr << "\t" << p.start << "\t" << p.end << "\t" << p.name << "\t" << p.length << "\t" << (p.strand == 1 ? "+" : "-") << "\t"  << p.merged << endl;
        return out;
    }
};

class my_gene: public my_block{
public:
    string region;
    unsigned int tss;
    vector < my_exon > exons;
    my_gene(){chr="NA";name="NA";start=0;end=0;gene_id="NA";strand = 0;region="NA"; length = 0 ;tss=0;read_count=0.0;gene_name = name;};
    ~my_gene(){exons.clear();}

    void assign(my_exon &e){
        if (strand && e.strand != strand) vrb.error("Strand mismatch");
        if (chr != "NA" && e.chr != chr) vrb.error("Chr mismatch");
        if (gene_name != "NA" && e.gene_name != gene_name) vrb.error("Name mismatch");
        if ( gene_id != "NA" && e.gene_id != gene_id) vrb.error("ID mismatch");
        strand = e.strand;
        gene_name = e.gene_name;
        chr = e.chr;
        gene_id = e.gene_id;
        if (start == 0 || e.start < start ) start = e.start;
        if (end == 0 || e.end > end) end = e.end;
        tss = strand == -1 ? end : start;
        vector < my_exon > overlaping;
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

    bool operator < (my_gene const & p) const {
        if (chr < p.chr) return true;
        if (chr > p.chr) return false;
        if (start < p.start) return true;
        if (start > p.start) return false;
        if (end < p.end) return true;
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

    friend ostream& operator<<(ostream& out, my_gene& g){
        out << g.chr << "\t" << g.start << "\t" << g.end << "\t" << g.gene_id << "\t" << g.length << "\t" << (g.strand == 1 ? "+" : "-") << "\t" << 2<<endl;
        sort(g.exons.begin(), g.exons.end());
        for (int i =0; i < g.exons.size(); i++ ) out << g.exons[i];
        return out;
    }

};

class quan2_data : public data {
public:
    map < string, unsigned int > genes_map;
    vector < my_gene > genes;
    string bam;
    my_stats  stats;
    string  sample;
    set < string > gene_types;
    filters filter;

    //void setChunk(int,int);
    //void setRegion(string);
    //genomic_region region;


    //void read_Sample_Names(vector < string > &);
    inline vector <my_gene> binary_find(vector < my_gene > &, my_block &);
    void readGTF(string);
    //void groupGenes();
    void readBam(string);
    my_cont_blocks keep_read(bam1_t *b);
    void printBEDcount(string);
    void printBEDrpkm(string);
    void printStats(string);
};

void quan2_main(vector < string > & );

inline static bool cmp_blocks (const my_block &a, const my_block &b){
	if (a.chr == b.chr && a.end >= b.start && a.start <= b.end) return false;
	else return (a < b);
}

#endif /* SRC_MODE_QUAN2_QUAN2_DATA_H_ */
