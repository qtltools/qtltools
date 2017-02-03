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

#include "rep_data.h"

rep_data::rep_data() {
	sample_count = 0;
	genotype_count = 0;
	phenotype_count = 0;
	covariate_count = 0;
	qtl_count = 0;
}

void rep_data::clear() {
	sample_count = 0;
	sample_id.clear();
	genotype_count = 0;
	genotype_val.clear();
	genotype_chr.clear();
	genotype_id.clear();
	genotype_start.clear();
	genotype_end.clear();
	phenotype_count = 0;
	phenotype_val.clear();
	phenotype_id.clear();
	phenotype_chr.clear();
	phenotype_start.clear();
	phenotype_end.clear();
	phenotype_neg.clear();
	covariate_count = 0;
	covariate_val.clear();
	covariate_id.clear();
}

rep_data::~rep_data() {
	clear();
}

void rep_data::residualizePhenotypes() {
	vrb.title("Residualize phenotypes for covariates");
	residualizer covariate_engine (sample_count);
	for (int c = 0 ; c < covariate_count ; c ++) covariate_engine.push(covariate_val[c]);
	covariate_engine.build();
	for (unsigned int p = 0 ; p < phenotype_count ; p ++) covariate_engine.residualize(phenotype_val[p]);
	vrb.bullet("#covariates = " + stb.str(covariate_engine.n_covariates));
}

