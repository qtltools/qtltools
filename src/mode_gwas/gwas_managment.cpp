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

#include "gwas_data.h"

void gwas_data::computeDosages(int * vGT, float * vDS) {
	for (int s = 0; s < sample_count; s ++) {
		if (vGT[2*s+0] == bcf_gt_missing || vGT[2*s+1] == bcf_gt_missing) vDS[s] = bcf_float_missing;
		else vDS[s] = bcf_gt_allele(vGT[2*s+0]) + bcf_gt_allele(vGT[2*s+1]);
	}
}

void gwas_data::imputeGenotypes(float * G) {
	double mean = 0.0;
	int c_mean= 0;
	for (int s = 0; s < sample_count; s ++) {
		if (G[s] != bcf_float_missing) {
			mean += G[s];
			c_mean ++;
		}
	}
	mean /= c_mean;
	for (int s = 0; s < sample_count ; s ++) if (G[s] == bcf_float_missing) G[s] = mean;
}

void gwas_data::imputePhenotypes() {
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

void gwas_data::imputeValues(vector < float > & V) {
	double mean = 0.0;
	int c_mean= 0;
	for (int s = 0; s < sample_count; s ++) {
		if (V[s] != bcf_float_missing) {
			mean += V[s];
			c_mean ++;
		}
	}
	mean /= c_mean;
	for (int s = 0; s < sample_count ; s ++) if (V[s] == bcf_float_missing) V[s] = mean;
}

void gwas_data::normalize(float * G) {
	double mean = 0.0, sum = 0.0;
	for (int s = 0; s < sample_count ; s ++) mean += G[s];
	mean /= sample_count;
	for (int s = 0; s < sample_count ; s ++) {
		G[s] -= mean;
		sum += G[s] * G[s];
	}
	sum = sqrt(sum);
	if (sum == 0) sum = 1;
	for (int s = 0; s < sample_count ; s ++) G[s] /= sum;
}

void gwas_data::normalize(vector < float > & G) {
	double mean = 0.0, sum = 0.0;
	for (int s = 0; s < sample_count ; s ++) mean += G[s];
	mean /= sample_count;
	for (int s = 0; s < sample_count ; s ++) {
		G[s] -= mean;
		sum += G[s] * G[s];
	}
	sum = sqrt(sum);
	if (sum == 0) sum = 1;
	for (int s = 0; s < sample_count ; s ++) G[s] /= sum;
}

void gwas_data::normalizePhenotypes() {
	for (int p = 0 ; p < phenotype_count ; p ++) normalize(phenotype_val[p]);
}

void gwas_data::normalTranformPhenotypes() {
	vrb.title("Normal transform phenotypes");
	for (int p = 0 ; p < phenotype_count ; p ++) {
		vector < float > R;
		myranker::rank(phenotype_val[p], R);
		double max = 0;
		for (int s = 0 ; s < sample_count ; s ++) {
			R[s] = R[s] - 0.5;
			if (R[s] > max) max = R[s];
		}
		max = max + 0.5;
		for (int s = 0 ; s < sample_count ; s ++) {
			R[s] /= max;
			phenotype_val[p][s] = qnorm(R[s], 0.0, 1.0, 1, 0);
		}
	}
}

void gwas_data::residualizePhenotypes() {
	vrb.title("Residualize phenotypes for covariates");
	residualizer covariate_engine (sample_count);
	for (int c = 0 ; c < covariate_count ; c ++) covariate_engine.push(covariate_val[c]);
	covariate_engine.build();
	for (unsigned int p = 0 ; p < phenotype_count ; p ++) covariate_engine.residualize(phenotype_val[p]);
	vrb.bullet("#covariates = " + stb.str(covariate_count));
}

