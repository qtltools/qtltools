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

#ifndef _CORRECT_DATA_H
#define _CORRECT_DATA_H

//ANALYSIS MODES
#define CORRECT_VCF	1
#define CORRECT_BED	2

//INCLUDES
#include "../common/data.h"

class correct_data : public data {
public:
	//PARAMETERS
	unsigned int mode;
	bool normalize;
	bool residualize;

	//COVARIATES
	int covariate_count;								//covariate number
	vector < vector < string > > covariate_val;			//covariate values
	vector < string > covariate_id;						//covariate ids
	vector < vector < string > > covariate_target;
	residualizer * covariate_engine;					//covariate engine machinery

	//QTL COVARIATES


	//CONSTRUCTOR / DESTRUCTOR
	correct_data();
	~correct_data();
	void clear();

	//READ DATA
	void readCovariates(string);
	void readQTLCovariates(string, string);

	//GENOTYPE & PHENOTYPE MANAGEMENT
	void imputeMissing(vector < float > &);
    void normalTransform(vector < float > &);
	void imputeMissing(float *);
    void normalTransform(float *);

    //PROCESS DATA
    void processBED(string, string);
    void processVCF(string, string);
};

//***************************************************************//
//******************** DECLARE FUNCTIONS *************************//
//***************************************************************//
void correct_main(vector < string > &);

#endif
