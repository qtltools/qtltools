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

#ifndef _RTC_DATA_H
#define _RTC_DATA_H

#define __INTERVAL_CENTRIC_RTC //do interval centric analysis

//ANALYSIS MODES
#define RTC_MODE1	1
#define RTC_MODE2	2
#define RTC_MODE3	3
#define RTC_MODE4	4

//AGGREGATION MODES
#define GRP_NONE	0
#define GRP_BEST	1
#define GRP_PCA1	2
#define GRP_MEAN	3

#define __RTC_NA__ (0.0/0.0)

//INCLUDES
#include "../common/data.h"


class coldspot{
public:
    string chr;
    int start;
    int end;
    int idx;
    string type;
    vector <int> coldspot_variant_idx;
    coldspot(){chr="";start=-1;end=-1;idx=-1;type="NA";}
    coldspot(string c, int s, int e, int i,string t){chr=c;start=s;end=e;idx=i; type = t ;}
    //~coldspot(){coldspot_variant_idx.clear();}
    friend ostream& operator<<(ostream& out, const coldspot& p){
        out << "ID: " << p.idx << " CHR: " << p.chr << " START: " << p.start << " END: " << p.end << " TYPE: " << p.type << " VARIANTS:";
        for (int i =0 ; i < p.coldspot_variant_idx.size(); i++) out << " " << p.coldspot_variant_idx[i];
        out << endl;
        return out;
    }
    long long unsigned int getMemoryUsage(){
    	return 3*sizeof(int) + type.capacity()*sizeof(char) + chr.capacity()*sizeof(char) + coldspot_variant_idx.capacity()*sizeof(int);
    }
};



class pairsToTestForRTC{
public:
    int test_snp_idx;
    int eqtl_snp_idx;
    int test_snp_rank;
    int eqtl_snp_rank;
    vector < int > other_conditional_signal_idx;
    int test_snp_coldspot_idx;
    int eqtl_snp_coldspot_idx;
    double Dprime;
    double R2;
    int pheno_idx;
    pairsToTestForRTC(){test_snp_idx=-1;eqtl_snp_idx=-1;test_snp_coldspot_idx=-1;eqtl_snp_coldspot_idx=-1;Dprime=-1.0;R2=0.0;test_snp_rank=0;eqtl_snp_rank=0;pheno_idx=-1;}
    pairsToTestForRTC(int tsi , int esi , vector <int> & ocsi, int tsci, int esci,double dp, double r2,int tr, int er , int pidx){eqtl_snp_rank=er;test_snp_rank=tr;test_snp_idx=tsi;eqtl_snp_idx=esi;test_snp_coldspot_idx=tsci;eqtl_snp_coldspot_idx=esci;other_conditional_signal_idx = ocsi;Dprime=dp;R2=r2;pheno_idx=pidx;}
    //~pairsToTestForRTC(){other_conditional_signal_idx.clear();}
    friend ostream& operator<<(ostream& out, const pairsToTestForRTC& p){
        out << "TSI: " << p.test_snp_idx << " ESI: " << p.eqtl_snp_idx << " TSCI: " << p.test_snp_coldspot_idx << " ESCI: " << p.eqtl_snp_coldspot_idx << " D: " << p.Dprime << " R: " <<p.R2 << " O:";
        for (int i =0 ; i < p.other_conditional_signal_idx.size(); i++) out << " " << p.other_conditional_signal_idx[i];
        out << endl;
        return out;
    }
    long long unsigned int getMemoryUsage(){
    	return 7*sizeof(int) + 2*sizeof(double) + other_conditional_signal_idx.capacity()*sizeof(int);
    }
};

class rtc_sample_results{
public:
	double gtoe_h0,gt_h0,gtoe_h1,gt_h1,count_h0,count_h1,pval;
	unsigned long int unique_h0,unique_h1,variants;
	float medianR2, median_h0, median_h1;
	double rtc_bin_start,rtc_bin_end,rtc_bin_h0_proportion, rtc_bin_h1_proportion;
	string h0,h1;
	rtc_sample_results(){pval=gtoe_h0=__RTC_NA__;gt_h0=__RTC_NA__;gtoe_h1=__RTC_NA__;gt_h1=__RTC_NA__;count_h0=0.0;count_h1=0.0;unique_h0=0,unique_h1=0; median_h0 = median_h1 = medianR2=__RTC_NA__; rtc_bin_start = 0.0,rtc_bin_end=1.0,rtc_bin_h0_proportion=1.0, rtc_bin_h1_proportion=1.0 ;h1=h0="";variants=0;}
    long long unsigned int getMemoryUsage(){
    	return sizeof(long int) + 11*sizeof(double) + 3*sizeof(float) + h0.capacity()*sizeof(char) + h1.capacity()*sizeof(char);
    }
};

class rtc_data : public data {
public:
	//PARAMETERS
	unsigned int mode;
	unsigned int grp_mode;
	unsigned int cis_window;
	int pvalue_column,variant_column,phenotype_column,rank_column,best_column,coldspot_count,group_column;
	static const int bin_size = 1000000;
    double Dprime_cutoff;
    double R2_cutoff;
    unsigned long int sample_iterations,max_sample_iterations;
    set <string> unphased_regions,unfound_regions,no_variance_regions,unmatched_alleles,unfound_ids,unfound_phenotypes;
    unsigned int DprimeR2inMem;
    static const int DprimePrintFreq = 10;
    static const int normal_output_columns = 18;
    bool calculate_Dprime_R2;

    //ADDITIONAL VCF
    filter stats_vcf_sample_filter;
    string stats_vcf_file;
    vector <int> stats_mappingS;
    int stats_n_includedS;

	//REGIONS
	genomic_region regionPhenotype;
	genomic_region regionGenotype;

	//GENOTYPES
	int genotype_count;									//variant site number
	vector < vector < float > > genotype_val;			//variant site genotype dosages
	vector < string > genotype_chr;						//variant site chromosome
	vector < string > genotype_id;						//variant site IDs
	vector < int > genotype_start;						//variant site start positions
	vector < int > genotype_end;						//variant site end positions
	unordered_map < string, int > genotype_id_to_idx;
    vector < string > genotype_alleles;

	//PHENOTYPES
	int phenotype_count;								//phenotype number
	vector < vector < float > > phenotype_val;			//phenotype values
	vector < string > phenotype_id;						//phenotype ids
	vector < string > phenotype_grp;					//phenotype groups
	vector < string > phenotype_chr;					//phenotype chromosomes
	vector < int > phenotype_start;						//phenotype start positions
	vector < int > phenotype_end;						//phenotype end positions
    unordered_map < string, int > phenotype_id_to_idx;

	//PHENOTYPE GROUPS
	vector < vector < unsigned int > > group_idx;		//group index to phenotype indexes
	vector < double > group_var;						//group variance explained by PC1
	vector < int > group_size;							//number of phenotypes in group
	map < string, unsigned int > group_id;

	//COVARIATES & INTERACTION
	int covariate_count;								//covariate number
	vector < vector < string > > covariate_val;			//covariate values
	vector < string > covariate_id;						//covariate ids

	//RTC
    map < string, map < int,  vector <coldspot *> > >  coldspot_bins_p;
    vector < coldspot > all_coldspots;
    //map < string, vector < int > > coldspot_end_idx;
    map < string ,vector < pairsToTestForRTC > > pheno_eqtls;
    unordered_map < string , unordered_map< string, vector <double> > > DprimeRsquareSink;
    //unordered_map < string , vector < float > > phenotypeSink; //UNUSED
    //unordered_map < string , vector < float > > phenotypeSinkRaw; //UNUSED
    unordered_map < int , vector < float > > genotypeSink;
    //unordered_map < int , unordered_map < int, float > > RsquareSink;
    vector < float> RsquareSink;


	//CONSTRUCTOR / DESTRUCTOR
	rtc_data();
	~rtc_data();
	void clear();

	//DATA REGION
	bool setPhenotypeRegion(string);
	void setPhenotypeRegion(int, int);
	void deduceGenotypeRegion(int);

	//READ DATA
	void readGenotypes(string);
	void readGenotypesVCF(string);
	void readGenotypesBED(string);
    void scanGenotypes(string);
	void scanGenotypesVCF(string);
	void scanGenotypesBED(string);
    int readGenotypesVCFStats(string, string &, vector <string> &);
	void readPhenotypes(string);
	void scanPhenotypes(string);
	void readCovariates(string);
    void readSampleInclusionStats(string);
    void readSampleExclusionStats(string);
    void copyIncludeExclude();
    void setStatsVCF(string);
    void checkStatsVCF();
    void createTransLists();
    void collapsePhenotypes();
    bool readRTCline(const string &buufer, set < int > &their);
    bool readRTCline(const string &buffer, string &pheno, string &snp, string &best, string &group, int &rank);
    bool readRTCline(const string &buffer, string &pheno, string &snp, string &group,unsigned int &line_count);

	//GENOTYPE & PHENOTYPE MANAGEMENT
	void clusterizePhenotypes(int);
	void imputeGenotypes();
    void imputePhenotypes();
    void residualizePhenotypes();
    void normalTransformPhenotypes();
    void normalTransform(vector < float > &);
	void normalize(vector < float > &);
	void normalize(vector < vector < float > > &);
	long long unsigned int getMemoryUsage();

	//COMPUTATION METHODS [ALL INLINES FOR SPEED]
	double getCorrelation(vector < float > &, vector < float > &);
    double getCorrelation(vector < float > &, vector < float > &, int);
	double getPvalue(double, double);
	double getPvalue(double, vector < double > &);
	double getSlope(double, double, double);
	void regression(vector < float > & X, vector < float > & Y, double & slope);
	float median(vector < float > &);

	//ANALYSIS
    void readHotspots(string);
    int getColdspot(string, int);
    void mapVariantsToColdspots();
    void calculateRTC(string);
    rtc_sample_results sampleRTC(vector <int> &, vector <float> &, int, double,int pI = -1);
    void generatePhenotype (double,int, vector <float> &);
    void generatePhenotype(vector<float> &, linReg &, vector <float>&);
    void mergeqtl_cis_conditional(string, string);
    void mergeqtl_cis(string, string);
    void mergeqtl_trans_conditional(string, string);
    void mergeqtl_trans(string, string);
    void gwas_cis_conditional(string, string);
    void gwas_cis(string, string);
    void gwas_trans_conditional(string, string);
    void gwas_trans(string, string);
    vector <double> getDprimeRsquare(string,int,string,int, string al1 = "" , string al2 = "", int idx1 = -1 , int idx2 = -1);
    double getRsquare(int,int,int,int);
    double getRsquare(int,int);
    void printPTTFR();
    vector < float > correct( vector < float> , vector <float> );
    int getBestVariant(vector <int> &, int, double &);
    int bcf_all_phased(const bcf_hdr_t *, bcf1_t *);
    void probability(vector < float > &, vector < float > & , double , rtc_sample_results &);

};

//***************************************************************//
//******************** DECLARE FUNCTIONS ************************//
//***************************************************************//
void rtc_main(vector < string > &);

//***************************************************************//
//******************** INLINE FUNCTIONS *************************//
//***************************************************************//

inline double rtc_data::getCorrelation(vector < float > & vec1, vector < float > & vec2) {
	int i = 0;
	int repeat = (sample_count / 4);
	int left = (sample_count % 4);
	double sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0;

	while (repeat --) {
		sum0 += vec1[i] * vec2[i];
		sum1 += vec1[i+1] * vec2[i+1];
		sum2 += vec1[i+2] * vec2[i+2];
		sum3 += vec1[i+3] * vec2[i+3];
		i += 4;
	}

	switch (left) {
	case 3:	sum0 += vec1[i+2] * vec2[i+2];
	case 2:	sum0 += vec1[i+1] * vec2[i+1];
	case 1:	sum0 += vec1[i+0] * vec2[i+0];
	case 0: ;
	}

	return sum0 + sum1 + sum2 + sum3;
}

inline double rtc_data::getCorrelation(vector < float > & vec1, vector < float > & vec2, int sc) {
    int i = 0;
    int repeat = (sc / 4);
    int left = (sc % 4);
    double sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0;
    
    while (repeat --) {
        sum0 += vec1[i] * vec2[i];
        sum1 += vec1[i+1] * vec2[i+1];
        sum2 += vec1[i+2] * vec2[i+2];
        sum3 += vec1[i+3] * vec2[i+3];
        i += 4;
    }
    
    switch (left) {
        case 3:	sum0 += vec1[i+2] * vec2[i+2];
        case 2:	sum0 += vec1[i+1] * vec2[i+1];
        case 1:	sum0 += vec1[i+0] * vec2[i+0];
        case 0: ;
    }
    
    return sum0 + sum1 + sum2 + sum3;
}

inline double rtc_data::getPvalue(double corr, double df) {
	return pf(df * corr * corr / (1 - corr * corr), 1, df, 0, 0);
}

inline double rtc_data::getPvalue(double ncorr, vector < double > & pcorr) {
	unsigned int n_hits = 0;
	for (int p = 0 ; p < pcorr.size() ; p++) if (abs(pcorr[p]) >= abs(ncorr)) n_hits++;
	return ((n_hits + 1) * 1.0 / (pcorr.size() + 1.0));
}

inline double rtc_data::getSlope(double nominal_correlation, double gsd, double psd) {
	if (gsd < 1e-16 || psd < 1e-16) return 0;
	else return nominal_correlation * psd / gsd;
}

inline void rtc_data::regression(vector < float > & X, vector < float > & Y, double & slope) {
	vector < float > Xtmp = X;
	vector < float > Ytmp = Y;
	double sdXtmp = basic_stats(Xtmp).sd();
	double sdYtmp = basic_stats(Ytmp).sd();
	normalize(Xtmp);
	normalize(Ytmp);
	double correlation = getCorrelation(Xtmp, Ytmp);
	//pvalue = getPvalue(correlation, sample_count - 2);
	slope = getSlope(correlation, sdXtmp, sdYtmp);
}

inline double rtc_data::getRsquare(int i, int j){
	vector < float > v1;
    if (genotypeSink.count(i)){
    	v1 = genotypeSink[i];
    }else{
		v1 = genotype_val[i];
		normalize(v1);
		if (DprimeR2inMem >=2) genotypeSink[i] = v1;
    }
    vector < float > v2;
    if (genotypeSink.count(j)){
    	v2 = genotypeSink[j];
    }else{
		v2 = genotype_val[j];
		normalize(v2);
		if (DprimeR2inMem >=2) genotypeSink[j] = v2;
    }
    double r = getCorrelation(v1 , v2 ,v1.size());
    return r * r;
}

inline double rtc_data::getRsquare(int i, int j,int start, int size){
	if (i==j) return 1.0;
	int first = -1 , second = -1;
	if (i > j){
		first = j - start + 1;
		second = i - start + 1;
	}else{
		first = i - start + 1;
		second = j - start + 1;
	}
	int total = size * (size - 1) / 2;
	if (!RsquareSink.size() && DprimeR2inMem >=2) RsquareSink = vector <float> (total, __RTC_NA__ );
	int diff = (size - first) * (size - first + 1) / 2;
	int index = total - diff + second - first - 1;
	if (index < RsquareSink.size() && !isnan(RsquareSink[index])) return RsquareSink[index];
	vector < float > v1;
    if (genotypeSink.count(i)){
    	v1 = genotypeSink[i];
    }else{
		v1 = genotype_val[i];
		normalize(v1);
		if (DprimeR2inMem >=2) genotypeSink[i] = v1;
    }
    vector < float > v2;
    if (genotypeSink.count(j)){
    	v2 = genotypeSink[j];
    }else{
		v2 = genotype_val[j];
		normalize(v2);
		if (DprimeR2inMem >=2) genotypeSink[j] = v2;
    }
    double r = getCorrelation(v1 , v2 ,v1.size());
    if (DprimeR2inMem >=2) RsquareSink[index] = r * r;
    return r * r;
}

inline int rtc_data::getBestVariant(vector <int> &genotype_idx, int phenotype_idx, double &pval){
	vector < float > y = phenotype_val[phenotype_idx];
	if (options.count("normal")) normalTransform(y);
	normalize(y);
	double bestR = 0.0;
	int bestV = -1;
	for (int g = 0; g< genotype_idx.size(); g++){
		vector < float > x = genotype_val[genotype_idx[g]];
		normalize(x);
		double R = abs(getCorrelation(x,y));
		if (R > bestR){
			bestR = R;
			bestV = genotype_idx[g];
		}
	}
	pval = getPvalue(bestR,sample_count-2);
	return bestV;
}

inline float rtc_data::median(vector <float> &scores){
  float median;
  size_t size = scores.size();
  if (!size) return __RTC_NA__;
  sort(scores.begin(), scores.end());
  if (size  % 2 == 0){
      median = (scores[size / 2 - 1] + scores[size / 2]) / 2;
  }else {
      median = scores[size / 2];
  }
  return median;
}

inline bool rtc_data::readRTCline(const string &buffer, set < int > &their){
	vector < string > str;
    stb.split(buffer, str);
    string snp = str[0];
    if (!genotype_id_to_idx.count(snp) ) {
    	unfound_ids.insert(snp);
    	return false;
    }
    their.insert(genotype_id_to_idx[snp]);
	return true;
}

inline bool rtc_data::readRTCline(const string &buffer, string &pheno, string &snp, string &best, string &group, int &rank){
	vector < string > str;
    stb.split(buffer, str);
    if (str.size() < 4) vrb.error("Wrong QTLtools output file format");
	if (str.size() < 4) vrb.error("Wrong QTLtools output file format");
	if (str[0] == "__UNION__"){
		if (str[2].substr(0,15) == "__UNION_FILLER_") return false;
		pheno = str[1];
		snp = str[2];
		best = "1";
		rank = atoi(str[3].c_str());
		group = pheno;
	}else{
		if (rank_column >= str.size()) vrb.error("rank column = " + stb.str(pvalue_column+1) + " but found " + stb.str(str.size()) + " columns in the following line:\n" + buffer);
		if (variant_column >= str.size()) vrb.error("variant column = " + stb.str(variant_column+1) + " but found " + stb.str(str.size()) + " columns in the following line:\n" + buffer);
		if (phenotype_column >= str.size()) vrb.error("phenotype column = " + stb.str(phenotype_column+1) + " but found " + stb.str(str.size()) + " columns in the following line:\n" + buffer);
		if (best_column >= str.size()) vrb.error("best column = " + stb.str(best_column+1) + " but found " + stb.str(str.size()) + " columns in the following line:\n" + buffer);
		if (group_column >= str.size()) vrb.error("group column = " + stb.str(group_column+1) + " but found " + stb.str(str.size()) + " columns in the following line:\n" + buffer);
		pheno = str[phenotype_column];
		group = str[group_column];
		snp = str[variant_column];
		best = str[best_column];
		rank = atoi(str[rank_column].c_str());
	}
	string test = grp_mode > 0 ? group : pheno;
	//cerr << pheno << " " << snp << " " << group << " " << test << " " << group_id.count(test) << endl;
	if(best != "1"){
		return false;
	}else if (!genotype_id_to_idx.count(snp)){
		unfound_ids.insert(snp);
		if (!phenotype_id_to_idx.count(test) && !group_id.count(test)) unfound_phenotypes.insert(test);
		return false;
	}else if (!phenotype_id_to_idx.count(test) && !group_id.count(test)){
		unfound_phenotypes.insert(test);
		return false;
	}
	return true;
}

inline bool rtc_data::readRTCline(const string &buffer, string &pheno, string &snp, string &group,unsigned int &line_count){
	vector < string > str;
    stb.split(buffer, str);
    if (str[0] == "__UNION__"){
    	if (str[2].substr(0,15) == "__UNION_FILLER_") return false;
    	pheno = str[1];
        snp = str[2];
        group = pheno;
    }else{
        if (!line_count && str.size() > normal_output_columns && !options.count("conditional")) vrb.warning("Looks like a conditional QTLtools output yet no --conditional provided, is this desired?");
		if (variant_column >= str.size()) vrb.error("variant column = " + stb.str(variant_column+1) + " but found " + stb.str(str.size()) + " columns in the following line:\n" + buffer);
		if (phenotype_column >= str.size()) vrb.error("phenotype column = " + stb.str(phenotype_column+1) + " but found " + stb.str(str.size()) + " columns in the following line:\n" + buffer);
		if (group_column >= str.size()) vrb.error("group column = " + stb.str(group_column+1) + " but found " + stb.str(str.size()) + " columns in the following line:\n" + buffer);
		pheno = str[phenotype_column];
		snp = str[variant_column];
		group = str[group_column];
    }
    line_count++;
    string test = grp_mode > 0 ? group : pheno;
    if (!genotype_id_to_idx.count(snp)){
    	unfound_ids.insert(snp);
    	if (!phenotype_id_to_idx.count(test) && !group_id.count(test)) unfound_phenotypes.insert(test);
    	return false;
    }else if (!phenotype_id_to_idx.count(test) && !group_id.count(test)){
    	unfound_phenotypes.insert(test);
    	return false;
    }
	return true;
}

#endif
