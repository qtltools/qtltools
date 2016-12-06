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

#include "correct_data.h"

correct_data::correct_data() {
	normalize = false;
	residualize = false;
	mode = 0;
	sample_count = 0;
	covariate_engine = NULL;
	covariate_count = 0;
}

correct_data::~correct_data() {
	normalize = false;
	residualize = false;
	sample_count = 0;
	sample_id.clear();
	covariate_count = 0;
	covariate_val.clear();
	covariate_id.clear();
	if (covariate_engine != NULL) delete covariate_engine;
	covariate_engine = NULL;
}

void correct_data::imputeMissing(vector < float > & V) {
	double mean = 0.0;
	int c_mean = 0;
	for (int s = 0; s < V.size() ; s ++) if (V[s] != bcf_float_missing) { mean += V[s]; c_mean ++; }
	mean /= c_mean;
	for (int s = 0; s < V.size() ; s ++) if (V[s] == bcf_float_missing) V[s] = mean;
}

void correct_data::normalTransform(vector < float > & V) {
	vector < float > R;
	myranker::rank(V, R);
	double max = 0;
	for (int s = 0 ; s < sample_count ; s ++) {
		R[s] = R[s] - 0.5;
		if (R[s] > max) max = R[s];
	}
	max = max + 0.5;
	for (int s = 0 ; s < sample_count ; s ++) {
		R[s] /= max;
		V[s] = qnorm(R[s], 0.0, 1.0, 1, 0);
	}
}

void correct_data::imputeMissing(float * V) {
	double mean = 0.0;
	int c_mean = 0;
	for (int s = 0; s < sample_count ; s ++) if (V[s] != bcf_float_missing) { mean += V[s]; c_mean ++; }
	mean /= c_mean;
	for (int s = 0; s < sample_count ; s ++) if (V[s] == bcf_float_missing) V[s] = mean;
}

void correct_data::normalTransform(float * V) {
	vector < float > R;
	myranker::rank(V, sample_count, R);
	double max = 0;
	for (int s = 0 ; s < sample_count ; s ++) {
		R[s] = R[s] - 0.5;
		if (R[s] > max) max = R[s];
	}
	max = max + 0.5;
	for (int s = 0 ; s < sample_count ; s ++) {
		R[s] /= max;
		V[s] = qnorm(R[s], 0.0, 1.0, 1, 0);
	}
}
