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

#ifndef _GWAS_DATA_H
#define _GWAS_DATA_H

//INCLUDES
#include "../common/data.h"

class gwas_data : public data {
public:
	//PHENOTYPES
	int phenotype_count;								//phenotype number
	vector < vector < float > > phenotype_val;			//phenotype values
	vector < string > phenotype_id;						//phenotype ids

	//COVARIATES
	int covariate_count;								//covariate number
	vector < vector < string > > covariate_val;			//covariate values
	vector < string > covariate_id;						//covariate ids

	//CONSTRUCTOR / DESTRUCTOR
    gwas_data();
    ~gwas_data();
    void clear();
    

	//READ OR GENERATE DATA
	void readPhenotypes(string);
	void readCovariates(string);

	//GENOTYPE & PHENOTYPE MANAGEMENT
	void computeDosages(int *, float *);
	void imputeGenotypes(float *);
	void imputePhenotypes();
	void imputeValues(vector < float > &);
	void normalizePhenotypes();
	void normalize(float *);
	void normalize(vector < float > &);
	void residualizePhenotypes();
	void normalTranformPhenotypes();

	//COMPUTATION METHODS [ALL INLINES FOR SPEED]
	double fastCorrelation(vector < float > &, float *);
	double fastCorrelation(vector < float > &, vector < float > &);
	double slowCorrelation(vector < float > &, float *);
	double getNominalPvalue(double, double);

	//ANALYSIS
	void runGwasPassOnVCF(string, string);
	void runGwasPassOnBED(string, string);
};

//***************************************************************//
//******************** DECLARE FUNCTIONS ************************//
//***************************************************************//
void gwas_main(vector < string > &);

//***************************************************************//
//******************** INLINE FUNCTIONS *************************//
//***************************************************************//

/*
 * Fast implementation of inner_product optimized for 64 bytes cache lines.
 */
inline double gwas_data::fastCorrelation(vector < float > & P, float * G) {
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

inline double gwas_data::fastCorrelation(vector < float > & P, vector < float > & G) {
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

inline double gwas_data::slowCorrelation(vector < float > & P, float * G) {
	double sum = 0.0;
	for (int e = 0 ; e < sample_count ; e ++) sum += P[e] * G[e];
	return sum;
}

inline double gwas_data::getNominalPvalue(double corr, double df) {
	double pval = pf(df * corr * corr / (1 - corr * corr), 1, df, 0, 0);
	if (pval <= std::numeric_limits<double>::min()) pval =std::numeric_limits<double>::min();
	return pval;
}

#endif
