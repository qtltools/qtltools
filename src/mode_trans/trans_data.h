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

#ifndef _TRANS_DATA_H
#define _TRANS_DATA_H

//ANALYSIS MODES
#define TRANS_MODE1		1
#define TRANS_MODE2		2
#define TRANS_MODE3		3
#define TRANS_MODE4		4

//INCLUDES
#include "../common/data.h"

class trans_data : public data {
public:
	//PARAMETERS
	unsigned int mode;
	double cis_window;
	double correlation_threshold;
	double beta_ml1, beta_ml2;
	unsigned int n_bins;
    unsigned long int start_line,end_line;
    
    //REGIONS
    genomic_region regionPhenotype;

	//PHENOTYPES
	int phenotype_count;								//phenotype number
	vector < vector < float > > phenotype_val;			//phenotype values
	vector < string > phenotype_id;						//phenotype ids
	vector < string > phenotype_chr;					//phenotype chromosomes
	vector < int > phenotype_start;						//phenotype start positions
	vector < int > phenotype_end;						//phenotype end positions
	vector < bool > phenotype_neg;						//phenotype is on the negative strand

	//COVARIATES
	int covariate_count;								//covariate number
	vector < vector < string > > covariate_val;			//covariate values
	vector < string > covariate_id;						//covariate ids

	//CONSTRUCTOR / DESTRUCTOR
    trans_data();
    ~trans_data();
    void clear();
    

	//READ OR GENERATE DATA
	void readPhenotypes(string);
    void scanPhenotypes(string);
    void setPhenotypeLines(int,int);
	void readSampleFromVCF(string);
	void checkSampleInBED(string);
	void checkSampleInCOV(string);
	void readCovariates(string);

	//GENOTYPE & PHENOTYPE MANAGEMENT
	void computeDosages(int *, float *);
	void imputeGenotypes(float *);
	void imputePhenotypes();
	void normalizePhenotypes();
	void normalize(float *);
	void normalize(vector < float > &);
	void residualizePhenotypes();
	void shufflePhenotypes();
	void samplePhenotypes(unsigned int);
	void normalTranformPhenotypes();

	//COMPUTATION METHODS [ALL INLINES FOR SPEED]
	double fastCorrelation(vector < float > &, float *);
	double slowCorrelation(vector < float > &, float *);
	double getNominalPvalue(double, double);
	double getAdjustedPvalue(double);
	void getCorrelationThreshold (double);
	int learnBetaParameters(vector < double > & , double &, double &);
	void buildNullDistribution(string);

	//ANALYSIS
	void runTransPass(string, string);
};

//***************************************************************//
//******************** DECLARE FUNCTIONS ************************//
//***************************************************************//
void trans_main(vector < string > &);

//***************************************************************//
//******************** INLINE FUNCTIONS *************************//
//***************************************************************//

/*
 * Fast implementation of inner_product optimized for 64 bytes cache lines.
 */
inline double trans_data::fastCorrelation(vector < float > & P, float * G) {
	int i = 0;
	int repeat = (sample_count / 4);
	int left = (sample_count % 4);
	double sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0;

	while (repeat --) {
		sum0 += P[i] * G[i];
		sum1 += P[i+1] * G[i+1];
		sum2 += P[i+2] * G[i+2];
		sum3 += P[i+3] * G[i+3];
		i += 4;
	}

	switch (left) {
	case 3:	sum0 += P[i+2] * G[i+2];
	case 2:	sum0 += P[i+1] * G[i+1];
	case 1:	sum0 += P[i+0] * G[i+0];
	case 0: ;
	}

	return sum0 + sum1 + sum2 + sum3;
}

inline double trans_data::slowCorrelation(vector < float > & P, float * G) {
	double sum = 0.0;
	for (int e = 0 ; e < sample_count ; e ++) sum += P[e] * G[e];
	return sum;
}

inline double trans_data::getNominalPvalue(double corr, double df) {
	double pval = pf(df * corr * corr / (1 - corr * corr), 1, df, 0, 0);
	if (pval <= std::numeric_limits<double>::min()) pval =std::numeric_limits<double>::min();
	return pval;
}

inline double trans_data::getAdjustedPvalue(double pv) {
	if (beta_ml1 == 0.0 || beta_ml2 == 0.0) return -1;
	else return pbeta(pv, beta_ml1, beta_ml2, 1, 0);
}

#endif
