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

#include "cis_data.h"

void cis_data::imputeGenotypes() {
	for (int g = 0; g < genotype_count ; g ++) {
		double mean = 0.0;
		int c_mean = 0;
		for (int s = 0; s < sample_count ; s ++) {
			if (genotype_val[g][s] != bcf_float_missing) {
				mean += genotype_val[g][s];
				c_mean ++;
			}
		}
		mean /= c_mean;
		for (int s = 0; s < sample_count ; s ++) if (genotype_val[g][s] == bcf_float_missing) genotype_val[g][s] = mean;
	}
}

void cis_data::imputePhenotypes() {
	for (int p = 0; p < phenotype_count ; p ++) {
		double mean = 0.0;
		int c_mean= 0;
		for (int s = 0; s < sample_count; s ++) {
			if (phenotype_val[p][s] != bcf_float_missing) {
				mean += phenotype_val [p][s];
				c_mean ++;
			}
		}
		mean /= c_mean;
		for (int s = 0; s < sample_count ; s ++) if (phenotype_val[p][s] == bcf_float_missing) phenotype_val[p][s] = mean;
	}
}

void cis_data::normalTransform(vector < float > & V) {
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

void cis_data::normalTransformPhenotypes() {
	vrb.title("Match phenotypes to Normal distribution");
	for (int p = 0; p < phenotype_count ; p ++) normalTransform(phenotype_val[p]);
}

void cis_data::normalize(vector < float > & X) {
	double mean = 0.0, sum = 0.0;
	for (int s = 0; s < sample_count ; s ++) mean += X[s];
	mean /= sample_count;
	for (int s = 0; s < sample_count ; s ++) {
		X[s] -= mean;
		sum += X[s] * X[s];
	}
	sum = sqrt(sum);
	if (sum == 0) sum = 1;
	for (int s = 0; s < sample_count ; s ++) X[s] /= sum;
}

void cis_data::normalize(vector < vector < float > > & X) {
	for (int x = 0 ; x < X.size() ; x++) {
		double mean = 0.0, sum = 0.0;
		for (int s = 0; s < sample_count ; s ++) mean += X[x][s];
		mean /= sample_count;
		for (int s = 0; s < sample_count ; s ++) {
			X[x][s] -= mean;
			sum += X[x][s] * X[x][s];
		}
		sum = sqrt(sum);
		if (sum == 0) sum = 1;
		for (int s = 0; s < sample_count ; s ++) X[x][s] /= sum;
	}
}
