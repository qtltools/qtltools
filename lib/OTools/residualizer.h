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

#ifndef _RESIDUALIZER_H
#define _RESIDUALIZER_H

#define R_QR_TOLERANCE 1e-7

//ERROR CODE FOR COVARIATE PARSING
#define COV_OKAY	0
#define COV_MIXD	1
#define COV_DROP	2
#define COV_NVAR	3
#define COV_NCOV	4
#define COV_CORR	5

//STL INCLUDES
#include <vector>
#include <set>
#include <sstream>
#include <string>

//EIGEN INCLUDES
#include <Eigen/Dense>
#include <Eigen/LU>

//EIGEN NAMESPACE
using namespace Eigen;

class residualizer {
public:
	unsigned int n_samples;
	unsigned int n_covariates;
	MatrixXd covarM;
    MatrixXd PQR_Q;
    MatrixXd PQR_Q_A;
    ColPivHouseholderQR < MatrixXd > PQR;

	residualizer(int _n_samples) : n_samples (_n_samples) {
    	covarM.resize(n_samples,1);
    	covarM.col(0) = VectorXd::Ones(n_samples);
    	n_covariates = 0;
    }

    ~residualizer() {
    	n_samples = 0;
    	n_covariates = 0;
    	covarM.resize(0,0);
        PQR_Q.resize(0,0);
        PQR_Q_A.resize(0,0);
    }

    unsigned int push(vector < string > & covariate) {
    	set < string > factors;
    	set < unsigned int > i_yesmissing, i_nonmissing;

    	//MAP MISSING
    	for (int i = 0 ; i < covariate.size() ; i++) {
    		if (covariate[i] == "NA") i_yesmissing.insert(i);
    		else i_nonmissing.insert(i);
    	}

    	//TEST FOR NUMERIC
    	bool isNumeric = false, isAlphabetic = false;
		for (set < unsigned int >::iterator itNM = i_nonmissing.begin(); itNM != i_nonmissing.end() ; ++itNM) {
			float value;
			std::istringstream in(covariate[*itNM]);
			if (!(in >> value)) {
				factors.insert(covariate[*itNM]);
				isAlphabetic = true;
			} else isNumeric = true;
		}
		if (isNumeric && isAlphabetic) return COV_MIXD;

		//FILL IN VALUES
		vector < vector < float > > additional_hcov;
		if (factors.size() == 0) {
			additional_hcov = vector < vector < float > > (1, vector < float > (n_samples, 0.0));
			for (int i = 0 ; i < covariate.size() ; i++) additional_hcov[0][i] = std::stof(covariate[i]);
		} else if (factors.size() > 1) {
			factors.erase(factors.begin());
			for (set < string > ::iterator itF = factors.begin(); itF != factors.end() ; itF++) {
				additional_hcov.push_back(vector < float > (n_samples, 0.0));
				for (int i = 0 ; i < n_samples ; i++) additional_hcov.back()[i] = (covariate[i] == (*itF));
			}
		} else return COV_DROP;

		//IMPUTE MISSING
		if (i_yesmissing.size() > 0) {
			for (unsigned int c = 0 ; c < additional_hcov.size() ; c ++) {
				double sum_row = 0;
				for (set < unsigned int >::iterator itNM = i_nonmissing.begin(); itNM != i_nonmissing.end() ; ++itNM) sum_row += additional_hcov[c][*itNM];
				for (set < unsigned int >::iterator itYM = i_yesmissing.begin(); itYM != i_yesmissing.end() ; ++itYM) additional_hcov[c][*itYM] = sum_row / i_nonmissing.size();
			}
		}

		//ADD RESULTING COVARIATES
		for (int c = 0 ; c < additional_hcov.size() ; c ++) if (!push(additional_hcov[c])) return COV_NVAR;

		return COV_OKAY;
    }

    bool push(vector < float > & covariate) {
    	bool isVariable = false;
    	for (unsigned int e = 1 ; e < covariate.size() ; e ++) if (covariate[e] != covariate[e-1]) isVariable = true;
    	if (!isVariable) return false;
    	n_covariates ++;
    	covarM.conservativeResize(n_samples, n_covariates+1);
    	for(int i = 0 ; i < n_samples ; i ++) covarM(i, n_covariates) = covariate[i];
    	return true;
    }

    unsigned int build() {
    	if (n_covariates == 0) return COV_NCOV;
    	PQR = ColPivHouseholderQR<MatrixXd>(covarM);
    	PQR.setThreshold(R_QR_TOLERANCE);
    	if (PQR.rank() != n_covariates + 1) {
    		PQR_Q = PQR.householderQ();
    		PQR_Q_A = PQR.householderQ().adjoint();
    		return COV_CORR;
    	}
    	return COV_OKAY;
    }

    unsigned int residualize(vector < float > & data) {
    	if (n_covariates == 0) return COV_NCOV;

    	bool isVariable = false;
    	for (unsigned int e = 1 ; e < data.size() ; e ++) if (data[e] != data[e-1]) isVariable = true;
    	if (!isVariable) return COV_NVAR;

    	//FILL IN DATA
    	VectorXd counts(n_samples);
    	for(int i = 0; i < n_samples ; i ++) counts(i) = data[i];

    	//CORRECTION
    	if (PQR.rank() == n_covariates + 1) {
    		VectorXd m_coef = PQR.solve(counts);
    		VectorXd fitted = covarM * m_coef;
    		VectorXd e = counts - fitted;
    		for (int i = 0; i < e.size(); i ++) data[i] = e(i);
    	} else {
    		VectorXd effects(PQR_Q_A * counts);
    		effects.tail(n_samples - PQR.rank()).setZero();
    		VectorXd fitted = PQR_Q * effects;
    		VectorXd e = counts - fitted;
    		for (int i = 0; i < e.size(); i ++) data[i] = (float)e(i);
    	}

    	return COV_OKAY;
    }


    unsigned int residualize(float * data) {
    	if (n_covariates == 0) return COV_NCOV;

    	bool isVariable = false;
    	for (unsigned int e = 1 ; e < n_samples ; e ++) if (data[e] != data[e-1]) isVariable = true;
    	if (!isVariable) return COV_NVAR;

    	//FILL IN DATA
    	VectorXd counts(n_samples);
    	for(int i = 0; i < n_samples ; i ++) counts(i) = data[i];

    	//CORRECTION
    	if (PQR.rank() == n_covariates + 1) {
    		VectorXd m_coef = PQR.solve(counts);
    		VectorXd fitted = covarM * m_coef;
    		VectorXd e = counts - fitted;
    		for (int i = 0; i < e.size(); i ++) data[i] = e(i);
    	} else {
    		VectorXd effects(PQR_Q_A * counts);
    		effects.tail(n_samples - PQR.rank()).setZero();
    		VectorXd fitted = PQR_Q * effects;
    		VectorXd e = counts - fitted;
    		for (int i = 0; i < e.size(); i ++) data[i] = (float)e(i);
    	}

    	return COV_OKAY;
    }
};

#endif
