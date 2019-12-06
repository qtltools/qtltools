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

class mapping_stats_full{
public:
	unsigned int fail_baseq,indel, duplicate, fail_qc,skipped, not_pp, mate_unmapped, orientation, depth, mapq , secondary;
	mapping_stats_full(){
		fail_baseq = indel = duplicate = fail_qc = skipped = not_pp = mate_unmapped = orientation = depth = mapq = secondary = 0;
	}

	void clear(){
		fail_baseq = indel = duplicate = fail_qc = skipped = not_pp = mate_unmapped = orientation = depth = mapq = secondary =0;
	}

	friend ostream& operator<<(ostream& out, mapping_stats_full& g){
		out << g.secondary << "\t" << g.mapq << "\t"<< g.fail_qc << "\t" << g.duplicate  << "\t"
		<< g.mate_unmapped  << "\t" << g.orientation << "\t" << g.not_pp  << "\t" << g.skipped  << "\t"
		<< g.fail_baseq   << "\t" << g.indel << "\t" << g.depth;
		return out;
	}
};

class mapping_stats{
public:
	unsigned int fail_baseq,indel,skipped,depth;
	mapping_stats(){
		fail_baseq = indel = skipped= depth = 0;
	}

	void clear(){
		fail_baseq = indel = skipped= depth = 0;
	}

	friend ostream& operator<<(ostream& out, mapping_stats& g){
		out << g.skipped  << "\t" << g.fail_baseq   << "\t" << g.indel << "\t" << g.depth;
		return out;
	}

};

class bias_stats{
public:
	string source;
	unsigned int ref_count,alt_count,total_sites,subsambled_sites,subsampled_to;
	double bias;
	bias_stats(){source="NA",bias=0.5,ref_count=alt_count=total_sites=subsambled_sites=subsampled_to=0;}
	bias_stats(unsigned int _r,unsigned int _a,unsigned int _s,unsigned int _t,unsigned int _m, double _b, string _g = "ALL_SITES"){
		ref_count = _r;
		alt_count = _a;
		total_sites = _s;
		subsambled_sites = _t;
		subsampled_to = _m;
		bias = _b;
		source = _g;
	}
	friend ostream& operator<<(ostream& out, bias_stats& g){
		out << g.ref_count << "\t" << g.alt_count << "\t" << g.total_sites  << "\t" << g.subsambled_sites << "\t" << g.subsampled_to  << "\t" << g.bias  << "\t" << g.source;
		return out;
	}
};


class ase_site {
public:
	unsigned int pos; //should be 1-based
	string chr;
	mutable char ref, alt;
	mutable unsigned int ref_count,alt_count,total_count,other_count;
	mutable double weighted_ref_count,weighted_alt_count,pval,ref_allele_mapping_bias,mar,discordant_pval,expected_error;
	mutable string alleles, sid;
	mutable mapping_stats stats;
	mutable set <string> alleles_seen;
	mutable string concern,genes;

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
		discordant_pval = pval = 1.0;
		ref_allele_mapping_bias = 0.5;
		mar = 0.0 / 0.0;
		concern = "";
		genes = "NA";
		expected_error = 0.0;
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
		discordant_pval = pval = 1.0;
		ref_allele_mapping_bias = 0.5;
		mar = 0.0 / 0.0;
		concern = "";
		genes = "NA";
		expected_error = 0.0;
	}

	~ase_site(){
		alleles_seen.clear();
	}


	void setCounts(const unsigned int r, const unsigned a, const unsigned o, const set <string> &as ,  const mapping_stats _stats = mapping_stats() , const double _em = 0.0){
		ref_count = r;
		alt_count = a;
		other_count = o;
		total_count = r + a;
		mar = (double) (r < a ? r : a) / (double) total_count;
		stats = _stats;
		alleles_seen = as;
		expected_error = _em;
	}

	void calculatePval(double rab = 0.5){
		ref_allele_mapping_bias = rab;
		weighted_ref_count = ref_count * (1 - ref_allele_mapping_bias);
		weighted_alt_count = alt_count * ref_allele_mapping_bias;
		if (other_count) discordant_pval = ppois(other_count-1, expected_error , 0, 0);
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

	string getName() const {
		return chr + ":" + stb.str(pos) + ":" + alleles + ":" + sid;
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

    friend ostream& operator<<(ostream& out, const ase_site& g) {
        out << g.sid  << "\t" << g.chr << "\t" << g.pos << "\t" << g.alleles << "\t" << (g.ref_count && g.alt_count) << "\t" << g.mar  << "\t" << g.ref_count
            << "\t" << g.alt_count << "\t" << g.total_count << "\t" << g.weighted_ref_count << "\t" << g.weighted_alt_count << "\t" << g.weighted_ref_count - g.weighted_alt_count << "\t";
        if (g.alleles_seen.size()){
        	ostringstream stream;
        	copy(g.alleles_seen.begin(), g.alleles_seen.end(), ostream_iterator<string>(stream, ":"));
        	out << stream.str().substr(0,stream.str().size() - 1);
        }else out << "NA";
		out << "\t" << g.ref << "\t" << g.alt << "\t" << g.other_count << "\t" << g.expected_error << "\t" << (g.discordant_pval <= std::numeric_limits<double>::min() ? std::numeric_limits<double>::min() : g.discordant_pval) << "\t" << g.ref_allele_mapping_bias << "\t" << (g.pval<= std::numeric_limits<double>::min() ? std::numeric_limits<double>::min() : g.pval) << "\t" << (g.concern == "" ? "PASS" : g.concern.substr(0,g.concern.size()-1)) << "\t" << g.genes;
        return out;
    }

};

class basic_block{
public:
	unsigned long  start, end; //should be 1-based
	mutable unsigned long length;
	mutable bool merged;
	string chr;
	basic_block(){start = 0; end = 0 ; chr="NA"; length = 0; merged=false;}
	basic_block(string c, unsigned long int s , unsigned long int e, bool m = false){
		assert(e >= s);
		start = s;
		end = e;
		chr = c;
		length = end - start + 1;
		merged = m;
	}
	basic_block(ase_site &in, bool m = false){
		start = in.pos;
		end = in.pos;
		chr = in.chr;
		length = end - start + 1;
		merged = m;
	}

	bool contains(unsigned long p) const{
		if (p<=end && p>=start) return true;
		else return false;
	}
	bool contains(string c ,unsigned long p) const{
		if (c == chr && p<=end && p>=start) return true;
		else return false;
	}

	bool overlap(basic_block &b) const{
		return (chr == b.chr && b.end >= start && b.start <= end);
	}

	bool overlap(genomic_region &b) const{
		return (chr == b.chr && b.end >= start && b.start <= end);
	}

	bool overlap(string c , unsigned long int s , unsigned long int e) const{
		return (chr == c && e >= start && s <= end);
	}

	bool contiguous(basic_block &b) const{
		return (overlap(b) || end + 1 == b.start || b.end + 1 == start);
	}

	vector < basic_block > subtract(basic_block &b) const{
		if (overlap(b)){
			if (start >= b.start && end <= b.end) return vector <basic_block>(0);
			if (start <  b.start && end <= b.end) return {basic_block(chr,start,b.start-1)};
			if (start >= b.start && end >  b.end) return {basic_block(chr,b.end+1,end)};
			return { basic_block(chr,start,b.start-1),basic_block(chr,b.end+1,end)};
		}else return {*this};
	}

	vector < basic_block > subtract(vector < basic_block > &b ) const{
		vector < basic_block> s {*this};
		for (int i = 0 ; i < b.size(); i++){
			if (i > 0 && (b[i-1].overlap(b[i]) || b[i] < b[i-1])) {
				cerr << "Cannot subtract overlapping or unsorted blocks" << endl;
				exit(10);
			}
			vector < basic_block > temp = s.back().subtract(b[i]);
			s.pop_back();
			if (!temp.size()) return s;
			s.insert(s.end(),temp.begin(),temp.end());
		}
		return s;
	}

	bool operator < (const basic_block &b) const {
		if (chr < b.chr) return true;
		if (chr > b.chr) return false;
		if (start < b.start) return true;
		if (start > b.start) return false;
		if (end < b.end) return true;
		else return false;
	}

	static bool cmp_blocks (const basic_block &a, const basic_block &b) {
		if (a.chr == b.chr && a.end >= b.start && a.start <= b.end) return false;
		else return (a < b);
	}

	basic_block merge(basic_block &b) const{
		if (overlap(b)) {
			return basic_block(b.chr, start < b.start ? start : b.start, end > b.end ? end : b.end , true);
		}else return basic_block();
	}

	basic_block merge_nocheck(basic_block &b) const{
			return basic_block(b.chr, start < b.start ? start : b.start, end > b.end ? end : b.end, true);
	}

	string get_string() const {
		return chr + ":" + stb.str(start)+ "-" + stb.str(end);
	}

	string get_string_merged() const {
		return chr + ":" + stb.str(start)+ "-" + stb.str(end) + "-" + stb.str(merged);
	}

	bool null() const {
		return (start == 0 && end == 0 && chr == "NA");
	}

	template <class B> vector <B> find_this_in(vector <B> &input) const {
		auto lower = lower_bound(input.begin(),input.end(),*this,this->cmp_blocks);
		auto upper = upper_bound(input.begin(),input.end(),*this,this->cmp_blocks);
		if (lower >= upper) return vector <B>(0);
		else return vector <B>(lower,upper);
	}

	template <class B> bool find_this_in_bool(vector <B> &input) const{
		auto lower = lower_bound(input.begin(),input.end(),*this,this->cmp_blocks);
		auto upper = upper_bound(input.begin(),input.end(),*this,this->cmp_blocks);
		if (lower >= upper) return false;
		else return true;
	}
};


class ase_region: public basic_block{
public:
	mutable unsigned int count;
	ase_region() : basic_block(){
		count = 0;
	}
	ase_region(string c, unsigned int s) : basic_block(c,s,s){
		count = 1;
	}
	ase_region(string c, unsigned int s, unsigned int e) : basic_block(c,s,e){
		count = 1;
	}
	ase_region(genomic_region &in) : basic_block(in.chr,in.start,in.end){
		count = 1;
	}
};

class ase_exon: public basic_block{
public:
	mutable string id;
	ase_exon() : basic_block(){id="NA";}
	ase_exon(unsigned int s, unsigned int e) : basic_block("NA",s,e){
		id = "NA";
	}

	ase_exon(string gid, string tid, string gn, unsigned int s, unsigned int e) : basic_block("NA",s,e){
		id = gid + ":" + tid + ":" + stb.str(start) + "_" + stb.str(end) + ":" + gn;
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
	bool param_dup_rd,param_both_alleles_seen,param_both_alleles_seen_bias,fix_chr,param_rm_indel,fix_id;
	bool keep_orphan,check_proper_pair,keep_failqc,legacy_options,auto_flip,check_orientation,print_stats,illumina13,on_the_fly,keep_discordant,print_warnings;
	string param_imputation_score_label,param_genotype_likelihood_label;
	static const int binsize = 10000, depth_flag_length = 10000;
	static constexpr float depth_fraction = 0.9f;

	genomic_region vcf_region,bam_region;

	//DATA
	vector < ase_region > my_regions;
	vector < vector < ase_site > > variants;
	set < ase_site >  all_variants;
	vector < ase_site >  passing_variants;
	map < string , double> ref_to_alt_bias;
	vector <string> all_allele_combinations;
	vector < basic_block> blacklisted_regions;
	set <string> add_chr,remove_chr;
	set <string> bam_chrs,vcf_chrs,ase_chrs;
	map <string, string> genome;
	map < string , map < unsigned int , vector < ase_exon > > > annotation;
	vector < basic_block> depth_exceeded;

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
		print_stats = false;
		illumina13 = false;
		fix_id = false;
		on_the_fly = false;
		keep_discordant = false;
		print_warnings = false;
	}

	~ase_data() {
		my_regions.clear();
		variants.clear();
		all_variants.clear();
		passing_variants.clear();
		ref_to_alt_bias.clear();
		all_allele_combinations.clear();
		blacklisted_regions.clear();
		add_chr.clear();
		remove_chr.clear();
		bam_chrs.clear();
		vcf_chrs.clear();
		ase_chrs.clear();
		genome.clear();
		annotation.clear();
	}

	//
	void readGenotypes(string v, string l = "");
	void readSequences(string);
	void readGTF(string );
	void readBlacklist(string);
	void readGenome(string);
	void calculateRefToAltBias(string o , string l = "");
	void calculateASE(string o , string l = "");
	void getRegions();
	void collapseRegions();
	void parseBam(void *);
	void compareChrs(string, string,string);
	void assignGenesToAseSite(ase_site &);
	inline char complement(string &in){
		if (in == "A") return 'T';
		if (in == "T") return 'A';
		if (in == "G") return 'C';
		if (in == "C") return 'G';
		return 'N';
	}
	inline char getBase(int code) {
		switch (code) {
		case 1: return 'A';
		case 2: return 'C';
		case 4: return 'G';
		case 8: return 'T';
		case 15: return 'N';
		}
		return -1;
	}
	inline void mergeContiguousBlocks(vector <basic_block> &input, vector <basic_block> &output){
		output.clear();
		if (input.size() > 1){
			sort(input.begin(),input.end());
			basic_block prev = input[0];
			for (int i = 1 ; i < input.size(); i++){
				if (prev.contiguous(input[i])) {
					prev = prev.merge_nocheck(input[i]);
				}else {
					output.push_back(prev);
					prev = input[i];
				}
			}
			output.push_back(prev);
		}else output = input;
	}
};

void ase_main(vector < string > & );




#endif
