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

#ifndef _ASE_DATA_H
#define _ASE_DATA_H

//INCLUDES
#include "../common/data.h"

class mapping_stats{
public:
	unsigned int fail_mapping,fail_baseq,indel, duplicate, unmapped, secondary, fail_qc,skipped, not_pp, mate_unmapped, orientation;
	mapping_stats(){
		fail_mapping = fail_baseq = indel = duplicate = 0;
		unmapped = secondary = fail_qc = skipped = not_pp = mate_unmapped = orientation = 0;
	}
	friend ostream& operator<<(ostream& out, mapping_stats& g){
		out << g.unmapped << "\t" << g.secondary << "\t" << g.fail_qc << "\t" << g.fail_mapping << "\t" << g.not_pp << "\t"
		<< g.mate_unmapped << "\t" << g.orientation << "\t" << g.duplicate << "\t" << g.indel << "\t"
		<< g.skipped << "\t" << g.fail_baseq;       
		return out;
	}
};


class ase_site {
public:
	unsigned int pos,ref_count,alt_count,total_count,other_count;
	double weighted_ref_count,weighted_alt_count,pval,ref_allele_mapping_bias,mar;
	string alleles;
	string chr, sid;
	char ref, alt;
	mapping_stats stats;

	ase_site(string _chr,unsigned int _pos){
		chr = _chr;
		sid = "NA";
		pos = _pos;
		ref = 'N';
		alt = 'N';
		alleles = ref;
		alleles += alt;
		other_count =ref_count= alt_count = total_count = 0;
		weighted_ref_count = weighted_alt_count = 0.0;
		pval = 1.0;
		ref_allele_mapping_bias = 0.5;
		mar = 0.0 / 0.0;
	}

	ase_site (string _chr, string _sid, unsigned int _pos, string _ref, string _alt) {
		chr = _chr;
		sid = _sid;
		pos = _pos;
		ref = _ref[0];
		alt = _alt[0];
		alleles = ref;
		alleles += alt;
		other_count =ref_count= alt_count = total_count = 0;
		weighted_ref_count = weighted_alt_count = 0.0;
		pval = 1.0;
		ref_allele_mapping_bias = 0.5;
		mar = 0.0 / 0.0;
	}


	void setCounts(const unsigned int r, const unsigned a, const unsigned o, const mapping_stats _stats = mapping_stats()){
		ref_count = r;
		alt_count = a;
		other_count = o;
		total_count = r + a;
		mar = (double) (r < a ? r : a) / (double) total_count;
		stats = _stats;
		if (other_count && ref_count == 0 && alt_count == 0) vrb.warning("No ref or alt allele for " + this->getName());
	}

	void calculatePval(double rab = 0.5){
		ref_allele_mapping_bias = rab;
		weighted_ref_count = ref_count * (1 - ref_allele_mapping_bias);
		weighted_alt_count = alt_count * ref_allele_mapping_bias;
		int y = 0;
		if (ref_allele_mapping_bias == 0) {pval = (ref_count == 0); return;}
		if (ref_allele_mapping_bias == 1) {pval = (ref_count == total_count); return;}
		double relErr = 1 + 1e-07;
		double d = dbinom(ref_count, total_count, ref_allele_mapping_bias, 0);
		double m = total_count * ref_allele_mapping_bias;
		if (ref_count == m) {pval = 1.0;return;}
		if (ref_count < m) {
			for (int i = (int)ceil (m); i <= total_count ; i++) y += (dbinom(i, total_count, ref_allele_mapping_bias, 0) <= d * relErr);
			pval = pbinom(ref_count, total_count, ref_allele_mapping_bias, 1, 0) + pbinom(total_count - y, total_count, ref_allele_mapping_bias, 0, 0);
		} else {
			for (int i = 0 ; i <= (int)floor(m) ; i++) y += (dbinom(i, total_count, ref_allele_mapping_bias, 0) <= d * relErr);
			pval = pbinom(y - 1, total_count, ref_allele_mapping_bias, 1, 0) + pbinom(ref_count - 1, total_count, ref_allele_mapping_bias, 0, 0);
		}
	}

	string getName(){
		return chr + ":" + stb.str(pos+1) + ":" + alleles + ":" + sid;
	}

	bool operator < (ase_site const & a) const {
		if (chr < a.chr) return true;
		if (chr > a.chr) return false;
		if (pos < a.pos) return true;
		//if (pos > a.pos) return false;
		//if (ref < a.ref) return true;
		//if (ref > a.ref) return false;
		//if (alt < a.alt) return true;
		return false;
	}

    friend ostream& operator<<(ostream& out, ase_site& g){
        out << g.chr << "\t" << g.pos + 1 << "\t" << g.sid << "\t" << g.ref << "\t" << g.alt << "\t" << g.alleles << "\t" << g.ref_count
            << "\t" << g.alt_count << "\t" << g.total_count << "\t" << g.mar << "\t" << g.other_count << "\t" << g.ref_allele_mapping_bias
			<< "\t" << g.weighted_ref_count << "\t" << g.weighted_alt_count << "\t" << g.pval << "\t" << g.stats;
        return out;
    }

};

class ase_basic_block{
public:
	unsigned long  start, end, length;
	bool merged;
	string chr;
	ase_basic_block(){start = 0; end = 0 ; chr="NA"; length = 0; merged=false;}
	ase_basic_block(string c, unsigned long int s , unsigned long int e, bool m = false){
		assert(e >= s);
		start = s;
		end = e;
		chr = c;
		length = end - start + 1;
		merged = m;
	}
	ase_basic_block(ase_site &in){
		ase_basic_block(in.chr,in.pos,in.pos);
	}

	bool contains(unsigned long p){
		if (p<=end && p>=start) return true;
		else return false;
	}
	bool contains(string c ,unsigned long p){
		if (c == chr && p<=end && p>=start) return true;
		else return false;
	}

	bool overlap(ase_basic_block &b){
		return (chr == b.chr && b.end >= start && b.start <= end);
	}

	bool overlap(string c , unsigned long int s , unsigned long int e){
		return (chr == c && e >= start && s <= end);
	}

	vector < ase_basic_block > subtract(ase_basic_block &b){
		if (overlap(b)){
			if (start >= b.start && end <= b.end) return vector <ase_basic_block>(0);
			if (start <  b.start && end <= b.end) return {ase_basic_block(chr,start,b.start-1)};
			if (start >= b.start && end >  b.end) return {ase_basic_block(chr,b.end+1,end)};
			return { ase_basic_block(chr,start,b.start-1),ase_basic_block(chr,b.end+1,end)};
		}else return {*this};
	}

	vector < ase_basic_block > subtract(vector < ase_basic_block > &b ){
		vector < ase_basic_block> s {*this};
		for (int i = 0 ; i < b.size(); i++){
			if (i > 0 && (b[i-1].overlap(b[i]) || b[i] < b[i-1])) {
				cerr << "Cannot subtract overlapping or unsorted blocks" << endl;
				exit(10);
			}
			vector < ase_basic_block > temp = s.back().subtract(b[i]);
			s.pop_back();
			if (!temp.size()) return s;
			s.insert(s.end(),temp.begin(),temp.end());
		}
		return s;
	}

	bool operator < (const ase_basic_block &b) const {
		if (chr < b.chr) return true;
		if (chr > b.chr) return false;
		if (start < b.start) return true;
		if (start > b.start) return false;
		if (end < b.end) return true;
		else return false;
	}

	static bool cmp_blocks (const ase_basic_block &a, const ase_basic_block &b){
		if (a.chr == b.chr && a.end >= b.start && a.start <= b.end) return false;
		else return (a < b);
	}

	ase_basic_block merge(ase_basic_block &b){
		if (overlap(b)) {
			return ase_basic_block(b.chr, start < b.start ? start : b.start, end > b.end ? end : b.end , true);
		}else return ase_basic_block();
	}

	ase_basic_block merge_nocheck(ase_basic_block &b){
			return ase_basic_block(b.chr, start < b.start ? start : b.start, end > b.end ? end : b.end, true);
	}

	string get_string(){
		return chr + ":" + stb.str(start)+ "-" + stb.str(end);
	}

	string get_string_merged(){
		return chr + ":" + stb.str(start)+ "-" + stb.str(end) + "-" + stb.str(merged);
	}

	bool null(){
		return (start == 0 && end == 0 && chr == "NA");
	}

	template <class B> vector <B> find_this_in(vector <B> &input){
		auto lower = lower_bound(input.begin(),input.end(),*this,this->cmp_blocks);
		auto upper = upper_bound(input.begin(),input.end(),*this,this->cmp_blocks);
		if (lower >= upper) return vector <B>(0);
		else return vector <B>(lower,upper);
	}

	template <class B> bool find_this_in_bool(vector <B> &input){
		auto lower = lower_bound(input.begin(),input.end(),*this,this->cmp_blocks);
		auto upper = upper_bound(input.begin(),input.end(),*this,this->cmp_blocks);
		if (lower >= upper) return false;
		else return true;
	}
};


class my_region{
public:
	string chr;
	unsigned int start,end,count;
	my_region(){
		chr = "";
		start = end = count = 0;
	}
	my_region(string c, unsigned int s){
		chr = c;
		start = s;
		end = s;
		count = 1;
	}
	my_region(string c, unsigned int s, unsigned int e){
		chr = c;
		start = s;
		end = e;
		count = 1;
	}

	my_region(genomic_region &in){
		chr = in.chr;
		start = in.start;
		end = in.end;
		count = 1;
	}

	string get(){
		return chr + ":" + stb.str(start) + "-" + stb.str(end);
	}
};

class ase_data : public data {
public :
	//PARAMETERS
	unsigned int param_min_mapQ;
	unsigned int param_min_baseQ;
	unsigned int param_min_cov;
	unsigned int param_min_cov_for_ref_alt;
	unsigned int param_min_sites_for_ref_alt;
	unsigned int region_length;
	int max_depth;
	float param_min_gp;
	float param_min_iq;
	float param_min_pval;
	float param_sample;
	bool param_dup_rd,param_both_alleles_seen,param_both_alleles_seen_bias,fix_chr,param_rm_indel;
	bool keep_orphan,check_proper_pair,keep_failqc,legacy_options,auto_flip,check_orientation;
	string param_imputation_score_label,param_genotype_likelihood_label;

	//DATA
	vector < string > regions;
	vector < my_region > my_regions;
	vector < vector < ase_site > > variants;
	set < ase_site >  all_variants;
	vector < ase_site >  passing_variants;
	map < string , double> ref_to_alt_bias;
	vector <string> all_allele_combinations;
	vector < ase_basic_block> blacklisted_regions;
	set <string> add_chr,remove_chr;
	set <string> bam_chrs;
	map <string, string> genome;

	//CONSTRUCTOR/DESTRUCTOR
	ase_data() {
		//these defaults get overwritten in ase_main.cpp
		param_min_mapQ = 10;
		param_min_baseQ = 13;
		param_min_pval = 1.0;
		param_min_gp = 0.99;
		param_min_iq = 0.90;
		param_min_cov = 16;
		param_min_cov_for_ref_alt = 10;
		param_dup_rd = false;
		param_both_alleles_seen = param_both_alleles_seen_bias = false;
		fix_chr = false;
		param_rm_indel = true;
		param_min_sites_for_ref_alt = 100;
		param_sample = 1.0;
		all_allele_combinations = {"AC", "AG" , "AT" , "CA" , "CG" , "CT" , "GA" , "GC" , "GT" , "TA" , "TC" , "TG"};
		for (int i = 0 ; i < all_allele_combinations.size(); i++ ) ref_to_alt_bias[all_allele_combinations[i]] = -1.0;
		param_imputation_score_label = "INFO";
		param_genotype_likelihood_label = "GL";
		region_length = 1000000;
		max_depth = 16000;
		keep_orphan = false;
		check_proper_pair = false;
		keep_failqc = false;
		legacy_options = false;
		auto_flip = false;
		check_orientation = false;
	}

	~ase_data() {
		regions.clear();
		variants.clear();
	}

	//
	void readGenotypes(string, string);
	void readGenotypes2(string v, string r, string l = "");
	void readSequences(string);
	void calculateRefToAltBias(string l = "");
	void calculateASE(string o , string l = "");
	void getRegions();
	void parseBam(void *);
	void parseBamMpileup(void *);
	void readBlacklist(string);
	void compareChrs(string, string);
	void readGenome(string);
	char complement(string &);


};

void ase_main(vector < string > & );

inline char ase_getBase (int code) {
	switch (code) {
	case 1: return 'A';
	case 2: return 'C';
	case 4: return 'G';
	case 8: return 'T';
	case 15: return 'N';
	}
	return -1;
}

inline double ase_binomialTest(int x, int n, float p) {
	int y = 0;
	if (p == 0) return (x == 0);
	if (p == 1) return (x == n);
	double relErr = 1 + 1e-07;
	double d = dbinom(x, n, p, 0);
	double m = n * p;
	if (x == m) return 1.0;
	if (x < m) {
		for (int i = (int)ceil (m); i <= n ; i++) y += (dbinom(i, n, p, 0) <= d * relErr);
		return pbinom(x, n, p, 1, 0) + pbinom(n - y, n, p, 0, 0);
	} else {
		for (int i = 0 ; i <= (int)floor(m) ; i++) y += (dbinom(i, n, p, 0) <= d * relErr);
		return pbinom(y - 1, n, p, 1, 0) + pbinom(x - 1, n, p, 0, 0);
	}
}

#endif
