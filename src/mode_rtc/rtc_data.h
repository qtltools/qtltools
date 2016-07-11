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

//ANALYSIS MODES
#define RTC_MODE1	1
#define RTC_MODE2	2
#define RTC_MODE3	3
#define RTC_MODE4	4

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
    friend ostream& operator<<(ostream& out, const coldspot& p){
        out << "ID: " << p.idx << " CHR: " << p.chr << " START: " << p.start << " END: " << p.end << " TYPE: " << p.type << " VARIANTS:";
        for (int i =0 ; i < p.coldspot_variant_idx.size(); i++) out << " " << p.coldspot_variant_idx[i];
        out << endl;
        return out;
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
    pairsToTestForRTC(){test_snp_idx=-1;eqtl_snp_idx=-1;test_snp_coldspot_idx=-1;eqtl_snp_coldspot_idx=-1;Dprime=-1.0;R2=0.0;test_snp_rank=0;eqtl_snp_rank=0;}
    pairsToTestForRTC(int tsi , int esi , vector <int> & ocsi, int tsci, int esci,double dp, double r2,int tr, int er){eqtl_snp_rank=er;test_snp_rank=tr;test_snp_idx=tsi;eqtl_snp_idx=esi;test_snp_coldspot_idx=tsci;eqtl_snp_coldspot_idx=esci;other_conditional_signal_idx = ocsi;Dprime=dp;R2=r2;}
    friend ostream& operator<<(ostream& out, const pairsToTestForRTC& p){
        out << "TSI: " << p.test_snp_idx << " ESI: " << p.eqtl_snp_idx << " TSCI: " << p.test_snp_coldspot_idx << " ESCI: " << p.eqtl_snp_coldspot_idx << " D: " << p.Dprime << " R: " <<p.R2 << " O:";
        for (int i =0 ; i < p.other_conditional_signal_idx.size(); i++) out << " " << p.other_conditional_signal_idx[i];
        out << endl;
        return out;
    }
};

class rtc_data : public data {
public:
	//PARAMETERS
	unsigned int mode;
	unsigned int cis_window;
	int pvalue_column,variant_column,phenotype_column,rank_column,best_column,coldspot_count;
	static const int bin_size = 1000000;
    double Dprime_cutoff;
    double R2_cutoff;
    unsigned long int sample_iterations;
    set <string> unphased_regions,unfound_regions,no_variance_regions,unmatched_alleles,unfound_ids,unfound_phenotypes;
    bool DprimeR2inMem;
    const int DprimePrintFreq = 10;
    const int normal_output_columns = 18;
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

	//COVARIATES & INTERACTION
	int covariate_count;								//covariate number
	vector < vector < string > > covariate_val;			//covariate values
	vector < string > covariate_id;						//covariate ids

	//RTC
    map < string, map < int,  vector <coldspot *> > >  coldspot_bins_p;
    vector < coldspot *> all_coldspots_p;
    map < string, vector < int > > coldspot_end_idx;
    map < int ,vector < pairsToTestForRTC > > pheno_eqtls;
    unordered_map < string , unordered_map< string, vector <double> > > DprimeRsquareSink;
    //unordered_map < string , vector < float > > phenotypeSink; //UNUSED
    //unordered_map < string , vector < float > > phenotypeSinkRaw; //UNUSED
    unordered_map < int , vector < float > > genotypeSink;
    unordered_map < int , unordered_map < int , double > > RsquareSink;

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

	//GENOTYPE & PHENOTYPE MANAGEMENT
	void clusterizePhenotypes(int);
	void imputeGenotypes();
    void imputePhenotypes();
    void residualizePhenotypes();
    void normalTransformPhenotypes();
    void normalTransform(vector < float > &);
	void normalize(vector < float > &);
	void normalize(vector < vector < float > > &);

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
    vector <double> sampleRTC(vector <int> &, vector <float> &, int, double,int pI = -1);
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
    double getRsquare(int,int);
    void printPTTFR();
    vector < float > correct( vector < float> , vector <float> );
    int getBestVariant(vector <int> &, int, double &);
    int bcf_all_phased(const bcf_hdr_t *, bcf1_t *);

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
    if (RsquareSink.count(i) && RsquareSink[i].count(j)) return RsquareSink[i][j];
    if (RsquareSink.count(j) && RsquareSink[j].count(i)) return RsquareSink[j][i];
	vector < float > v1;
    if (genotypeSink.count(i)){
    	v1 = genotypeSink[i];
    }else{
		v1 = genotype_val[i];
		normalize(v1);
		if (DprimeR2inMem) genotypeSink[i] = v1;
    }
    vector < float > v2;
    if (genotypeSink.count(j)){
    	v2 = genotypeSink[j];
    }else{
		v2 = genotype_val[j];
		normalize(v2);
		if (DprimeR2inMem) genotypeSink[j] = v2;
    }
    double r = getCorrelation(v1 , v2 ,v1.size());
    if (DprimeR2inMem) RsquareSink[i][j] = r * r;
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
  sort(scores.begin(), scores.end());
  if (size  % 2 == 0){
      median = (scores[size / 2 - 1] + scores[size / 2]) / 2;
  }else {
      median = scores[size / 2];
  }
  return median;
}

#endif
