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

#ifndef _rep_data_H
#define _rep_data_H

//ANALYSIS MODES
#define CIS_PERM	1
#define CIS_NOMI	2
#define CIS_COND	3

//AGGREGATION MODES
#define GRP_NONE	0
#define GRP_BEST	1
#define GRP_PCA1	2
#define GRP_MEAN	3

//INCLUDES
#include "../common/data.h"

class rep_data : public data {
public:
	//GENOTYPES
	int genotype_count;									//variant site number
	vector < vector < float > > genotype_val;			//variant site genotype dosages
	vector < string > genotype_chr;						//variant site chromosome
	vector < string > genotype_id;						//variant site IDs
	vector < int > genotype_start;						//variant site start positions
	vector < int > genotype_end;						//variant site end positions

	//PHENOTYPES
	int phenotype_count;								//phenotype number
	vector < vector < float > > phenotype_val;			//phenotype values
	vector < string > phenotype_id;						//phenotype ids
	vector < string > phenotype_chr;					//phenotype chromosomes
	vector < int > phenotype_start;						//phenotype start positions
	vector < int > phenotype_end;						//phenotype end positions
	vector < bool > phenotype_neg;						//phenotype is on the negative strand

	//COVARIATES
	int covariate_count;									//covariate number
	vector < vector < string > > covariate_val;				//covariate values
	vector < string > covariate_id;							//covariate ids

	//QTLS
	int qtl_count;
	pair < vector < string >, vector < string > > qtl_ids;
	pair < vector < int >, vector < int > > qtl_idx;

	//CONSTRUCTOR / DESTRUCTOR
	rep_data();
	~rep_data();
	void clear();

	//READ DATA
	void readGenotypes(string);
	void readGenotypesVCF(string);
	void readGenotypesBED(string);
	void readPhenotypes(string);
	void readCovariates(string);
	void readQTLs(string);

	//GENOTYPE & PHENOTYPE MANAGEMENT
	void imputeGenotypes();
    void imputePhenotypes();
    void residualizePhenotypes();
    void normalTransformPhenotypes();
    void normalTransform(vector < float > &);
	void normalize(vector < float > &);
	void normalize(vector < vector < float > > &);
	void mapping();

	//COMPUTATION METHODS [ALL INLINES FOR SPEED]
	double getCorrelation(vector < float > &, vector < float > &);
	double getPvalue(double, double);
	double getPvalue(double, vector < double > &);
	double getSlope(double, double, double);
	void regression(vector < float > & X, vector < float > & Y, double & pvalue, double & slope);

	//ANALYSIS
	void writeHeader(string);
	void runNominalPass(string);
};

//***************************************************************//
//******************** DECLARE FUNCTIONS *************************//
//***************************************************************//
void rep_main(vector < string > &);

//***************************************************************//
//******************** INLINE FUNCTIONS *************************//
//***************************************************************//

inline double rep_data::getCorrelation(vector < float > & vec1, vector < float > & vec2) {
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

inline double rep_data::getPvalue(double corr, double df) {
	double pval = pf(df * corr * corr / (1 - corr * corr), 1, df, 0, 0);
	if (pval <= std::numeric_limits<double>::min()) pval =std::numeric_limits<double>::min();
	return pval;
}

inline double rep_data::getPvalue(double ncorr, vector < double > & pcorr) {
	unsigned int n_hits = 0;
	for (int p = 0 ; p < pcorr.size() ; p++) if (abs(pcorr[p]) >= abs(ncorr)) n_hits++;
	return ((n_hits + 1) * 1.0 / (pcorr.size() + 1.0));
}

inline double rep_data::getSlope(double nominal_correlation, double gsd, double psd) {
	if (gsd < 1e-16 || psd < 1e-16) return 0;
	else return nominal_correlation * psd / gsd;
}

inline void rep_data::regression(vector < float > & X, vector < float > & Y, double & pvalue, double & slope) {
	vector < float > Xtmp = X;
	vector < float > Ytmp = Y;
	double sdXtmp = basic_stats(Xtmp).sd();
	double sdYtmp = basic_stats(Ytmp).sd();
	normalize(Xtmp);
	normalize(Ytmp);
	double correlation = getCorrelation(Xtmp, Ytmp);
	pvalue = getPvalue(correlation, sample_count - 2);
	slope = getSlope(correlation, sdXtmp, sdYtmp);
}

#endif
