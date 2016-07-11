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

#include "trans_data.h"

void trans_data::computeDosages(int * vGT, float * vDS) {
	for (int s = 0; s < sample_count; s ++) {
		if (vGT[2*s+0] == bcf_gt_missing || vGT[2*s+1] == bcf_gt_missing) vDS[s] = bcf_float_missing;
		else vDS[s] = bcf_gt_allele(vGT[2*s+0]) + bcf_gt_allele(vGT[2*s+1]);
	}
}

void trans_data::imputeGenotypes(float * G) {
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

void trans_data::imputePhenotypes() {
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

void trans_data::normalize(float * G) {
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

void trans_data::normalize(vector < float > & G) {
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

void trans_data::normalizePhenotypes() {
	for (int p = 0 ; p < phenotype_count ; p ++) normalize(phenotype_val[p]);
}

void trans_data::normalTranformPhenotypes() {
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

void trans_data::shufflePhenotypes() {
	vector < float > phenotype_tmp = vector < float > (sample_count, 0.0);
	vector < unsigned int > O = vector < unsigned int > (sample_count, 0);
	for (int i = 0 ; i < sample_count ; i ++) O[i] = i;
	shuffle(O.begin(), O.end(), rng.getEngine());
	for (unsigned int p = 0 ; p < phenotype_count ; p ++) {
		for (unsigned int i = 0 ; i < sample_count ; i ++) phenotype_tmp[i] = phenotype_val[p][O[i]];
		phenotype_val[p] = phenotype_tmp;
	}
}

void trans_data::samplePhenotypes(unsigned int N) {
	int phenotype_count_tmp = N;
	vector < vector < float > > phenotype_val_tmp = vector < vector < float > > (N);
	vector < string > phenotype_id_tmp = vector < string > (N);
	vector < string > phenotype_chr_tmp = vector < string > (N);
	vector < int > phenotype_start_tmp = vector < int > (N);
	vector < int > phenotype_end_tmp = vector < int > (N);

	for (int n = 0; n < phenotype_count_tmp; n ++) {
		unsigned int ridx = rng.getInt(phenotype_count);
		phenotype_val_tmp[n] = phenotype_val[ridx];
		shuffle(phenotype_val_tmp[n].begin(), phenotype_val_tmp[n].end(), rng.getEngine());
		phenotype_id_tmp[n] = phenotype_id[ridx];
		phenotype_chr_tmp[n] = phenotype_chr[ridx];
		phenotype_start_tmp[n] = phenotype_start[ridx];
		phenotype_end_tmp[n] = phenotype_end[ridx];
	}
	phenotype_count = phenotype_count_tmp;
	phenotype_val = phenotype_val_tmp;
	phenotype_id = phenotype_id_tmp;
	phenotype_chr = phenotype_chr_tmp;
	phenotype_start = phenotype_start_tmp;
	phenotype_end = phenotype_end_tmp;
}

void trans_data::getCorrelationThreshold(double pvalue) {
	vrb.title("Calculate correlation threshold");
	vrb.bullet("thres = " + stb.str(pvalue));
	if (beta_ml1 == 0 && beta_ml2 == 0) {
		double p = qf(pvalue, 1, sample_count - 2, 0, 0);
		correlation_threshold = sqrt(p / (sample_count - 2 + p));
	} else {
		double p0 = qbeta(pvalue, beta_ml1, beta_ml2, 1, 0);
		vrb.bullet("npvt = " + stb.str(p0));
		double p1 = qf(p0, 1, sample_count - 2, 0, 0);
		vrb.bullet("nfst = " + stb.str(p1));
		correlation_threshold = sqrt(p1 / (sample_count - 2 + p1));
	}
	vrb.bullet("corr = " + stb.str(correlation_threshold, 4));
}

void trans_data::residualizePhenotypes() {
	vrb.title("Residualize phenotypes for covariates");
	residualizer covariate_engine (sample_count);
	for (int c = 0 ; c < covariate_count ; c ++) covariate_engine.push(covariate_val[c]);
	covariate_engine.build();
	for (unsigned int p = 0 ; p < phenotype_count ; p ++) covariate_engine.residualize(phenotype_val[p]);
	vrb.bullet("#covariates = " + stb.str(covariate_count));
}

