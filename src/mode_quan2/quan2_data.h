/*
 * quan2_data.h
 *
 *  Created on: Aug 21, 2017
 *      Author: halit
 */

#ifndef SRC_MODE_QUAN2_QUAN2_DATA_H_
#define SRC_MODE_QUAN2_QUAN2_DATA_H_

#include "../common/data.h"

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

#endif /* SRC_MODE_QUAN2_QUAN2_DATA_H_ */
