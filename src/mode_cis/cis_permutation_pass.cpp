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

void cis_data::runPermutationPass(string fout) {

	//STEP0: INITIALIZATION OF IO
	output_file fdo (fout);
	if (fdo.fail()) vrb.error("Cannot open file [" + fout + "]");

	//STEP1: INITIALIZE A WORKING COPY OF GENOTYPES
	vector < double > genotype_sd = vector < double > (genotype_count, 0.0);
	for (unsigned int v = 0 ; v < genotype_count ; v ++) {
		genotype_sd[v] = basic_stats(genotype_val[v]).sd();
		normalize(genotype_val[v]);
	}

	//STEP2: INITIALIZE A WORKING COPY OF PHENOTYPES
	vector < double > phenotype_sd = vector < double > (phenotype_count, 0.0);
	for (unsigned int p = 0 ; p < phenotype_count ; p ++) {
		phenotype_sd[p] = basic_stats(phenotype_val[p]).sd();
		normalize(phenotype_val[p]);
	}

	//STEP3: MAIN SWEEP THROUGH PHENOTYPES
	for (unsigned int i_group = 0 ; i_group < group_idx.size() ; i_group ++) {

		//STEP4: VERBOSE PROCESSED PHENOTYPES
		if (grp_mode == GRP_NONE) vrb.title("Processing phenotype [" + phenotype_id[group_idx[i_group][0]] + "] [" + stb.str(i_group+1) + "/" + stb.str(group_idx.size()) + "]");
		else {
			vrb.title("Processing group of phenotypes [" + phenotype_grp[group_idx[i_group][0]] + "] [" + stb.str(i_group+1) + "/" + stb.str(group_idx.size()) + "]");
			vrb.bullet("#phenotypes in group = " + stb.str(group_size[i_group]));
			if (grp_mode == GRP_PCA1) vrb.bullet("variance explained by PC1 = " + stb.str(group_var[i_group], 3));
		}

		//STEP5: ENUMERATE ALL VARIANTS IN CIS
		vector < unsigned int > variant_indexes;
		vector < int > variant_distances;
		for (unsigned int v = 0 ; v < genotype_count ; v ++) {
			if (!full_test && phenotype_chr[group_idx[i_group][0]] != genotype_chr[v]) continue;
			int ps = (phenotype_start[group_idx[i_group][0]]>cis_window)?(phenotype_start[group_idx[i_group][0]]-cis_window):0;
			int pe = phenotype_end[group_idx[i_group][0]] + cis_window;

			if (full_test || (genotype_start[v] <= pe && ps <= genotype_end[v])) {
				int cisdistance = 0;
				if (genotype_start[v] <= phenotype_end[group_idx[i_group][0]] && phenotype_start[group_idx[i_group][0]] <= genotype_end[v]) cisdistance = 0;
				else if (genotype_end[v] < phenotype_start[group_idx[i_group][0]]) cisdistance = genotype_end[v] - phenotype_start[group_idx[i_group][0]];
				else cisdistance = genotype_start[v] - phenotype_end[group_idx[i_group][0]];
				if (phenotype_neg[group_idx[i_group][0]]) cisdistance *= -1;
				variant_indexes.push_back(v);
				variant_distances.push_back(cisdistance);
			}
		}
		vrb.bullet("#variants in cis = " + stb.str(variant_indexes.size()));

		//STEP7: NO VARIANTS IN CIS: OUTPUT NAs
		if (variant_indexes.size() == 0) {
			if (grp_mode == GRP_NONE) fdo << phenotype_id[group_idx[i_group][0]];
			else fdo << phenotype_grp[group_idx[i_group][0]];
			fdo << " " << phenotype_chr[group_idx[i_group][0]];
			fdo << " " << phenotype_start[group_idx[i_group][0]];
			fdo << " " << phenotype_end[group_idx[i_group][0]];
			fdo << " " << (phenotype_neg[group_idx[i_group][0]]?"-":"+");
			switch (grp_mode) {
			case GRP_BEST: fdo << " NA " << stb.str(group_size[i_group]); break;
			case GRP_PCA1: fdo << " " << stb.str(group_var[i_group], 3) << " " << stb.str(group_size[i_group]); break;
			case GRP_MEAN: fdo << " " << stb.str(group_size[i_group]); break;
			}
			fdo << " " << variant_indexes.size() << " NA NA NA NA NA NA NA NA NA NA NA NA NA" << endl;
		}

		//STEP8: VARIANTS IN CIS FOUND: FULL COMPUTATIONS
		else {
			double best_nominal_correlation = 0.0;
			int best_nominal_variant_abs = -1;
			int best_nominal_phenotype_abs = -1;
			int best_nominal_distance = cis_window;
			vector < double > best_permuted_correlations = vector < double > (n_permutations, 0.0);

			//STEP9: NOMINAL PASS
			for (unsigned int p = 0 ; p < group_idx[i_group].size() ; p ++) {
				for (unsigned int v = 0 ; v < variant_indexes.size() ; v ++) {
					double curr_correlation = getCorrelation(genotype_val[variant_indexes[v]], phenotype_val[group_idx[i_group][p]]);
					if (abs(curr_correlation) > abs(best_nominal_correlation) || (abs(curr_correlation) == abs(best_nominal_correlation) && abs(variant_distances[v]) < abs(best_nominal_distance))) {
						best_nominal_correlation = curr_correlation;
						best_nominal_variant_abs = variant_indexes[v];
						best_nominal_phenotype_abs = group_idx[i_group][p];
						best_nominal_distance = variant_distances[v];
					}
				}
			}

			//STEP10: PERMUTATION PASS
			vector < float > phenotype_curr = vector < float > (sample_count, 0.0);
			vector < unsigned int > permuted_indexes = vector < unsigned int > (sample_count, 0);
			for (unsigned int i = 0 ; i < sample_count ; i ++) permuted_indexes[i] = i;
			for (int perm = 0 ; perm < n_permutations ; perm ++) {
				shuffle(permuted_indexes.begin(), permuted_indexes.end(), rng.getEngine());
				for (unsigned int p = 0 ; p < group_idx[i_group].size() ; p ++) {
					for (unsigned int i = 0 ; i < sample_count ; i ++) phenotype_curr[i] = phenotype_val[group_idx[i_group][p]][permuted_indexes[i]];
					for (unsigned int v = 0 ; v < variant_indexes.size() ; v ++) {
						double curr_correlation = getCorrelation(genotype_val[variant_indexes[v]], phenotype_curr);
						if (abs(curr_correlation) > abs(best_permuted_correlations[perm])) best_permuted_correlations[perm] = curr_correlation;
					}
				}
			}

			//STEP11: COMPUTE BASIC STATS FOR BEST HIT
			double dof_true = sample_count - 2;
			double pval_nom = getPvalue(best_nominal_correlation, dof_true);
			double pval_slope = getSlope(best_nominal_correlation, genotype_sd[best_nominal_variant_abs], phenotype_sd[best_nominal_phenotype_abs]);
			double pval_r2 = best_nominal_correlation * best_nominal_correlation;
			double pval_se = 0.0;
			if(std_err) pval_se = getSE(pval_r2, genotype_sd[best_nominal_variant_abs], phenotype_sd[best_nominal_phenotype_abs]);

			//STEP12: VERBOSE BEST HIT
			if (grp_mode == GRP_BEST)
				vrb.bullet("Best hit: [id=" + genotype_id[best_nominal_variant_abs] + ", d=" + stb.str(best_nominal_distance) + ", p=" + phenotype_id[best_nominal_phenotype_abs] + ", pv=" + stb.str(pval_nom) + ", s=" + stb.str(pval_slope, 4) + "]");
			else vrb.bullet("Best hit: [id=" + genotype_id[best_nominal_variant_abs] + ", d=" + stb.str(best_nominal_distance) + ", pv=" + stb.str(pval_nom) + ", s=" + stb.str(pval_slope, 4) + "]");

			//STEP13: PROCESS DEGREES OF FREEDOM
			double dof_esti = dof_true;
			double variance_best_permuted_correlations = basic_stats(best_permuted_correlations).variance();
			if (variance_best_permuted_correlations != 0) learnDegreeOfFreedom(best_permuted_correlations, dof_esti);
			//vrb.bullet("DOF: [t=" + stb.str(dof_true, 1) + ", e=" + stb.str(dof_esti, 1) + "]");

			//STEP14: COMPUTE BEST PERMUTATION HIT
			vector < double > best_permuted_pvalues = vector < double > (n_permutations, 0.0);
			for (int perm = 0 ; perm < n_permutations ; perm ++) best_permuted_pvalues[perm] = getPvalue(best_permuted_correlations[perm], dof_esti);
			double mean_best_permuted_pvalues, variance_best_permuted_pvalues;
			mean_best_permuted_pvalues = basic_stats(best_permuted_pvalues).mean();
			variance_best_permuted_pvalues = basic_stats(best_permuted_pvalues).variance();

			//STEP15: LEARN BETA PARAMETERS
			double beta_mm1 = mean_best_permuted_pvalues * (mean_best_permuted_pvalues * (1 - mean_best_permuted_pvalues ) / variance_best_permuted_pvalues - 1);
			double beta_mm2 = beta_mm1 * (1 / mean_best_permuted_pvalues - 1);
			double beta_ml1 = beta_mm1;
			double beta_ml2 = beta_mm2;
			try {
				learnBetaParameters(best_permuted_pvalues, beta_ml1, beta_ml2);
			} catch (const std::exception & e) {
				vrb.bullet("Maximum Likelihood estimation failed, use Moment Matching instead!");
				beta_ml1 = beta_mm1;
				beta_ml2 = beta_mm2;
			}
			vrb.bullet("Beta parameters: [s1=" + stb.str(beta_ml1) + ", s2=" + stb.str(beta_ml2) +"]");

			//STEP16: COMPUTE ADJUSTED PVALUES
			double pval_emp = getPvalue(best_nominal_correlation, best_permuted_correlations);
			double pval_bml = pbeta(getPvalue(best_nominal_correlation, dof_esti), beta_ml1, beta_ml2, 1, 0);

			//STEP17: VERBOSE ADJUSTED PVALUES
			vrb.bullet("Adjusted p-values: [emp=" + stb.str(pval_emp) + ", beta=" + stb.str(pval_bml) + "]");

			//STEP18: PRINT RESULTS IN FILE
			if (grp_mode == GRP_NONE) fdo << phenotype_id[group_idx[i_group][0]];
			else fdo << phenotype_grp[group_idx[i_group][0]];
			fdo << " " << phenotype_chr[group_idx[i_group][0]];
			fdo << " " << phenotype_start[group_idx[i_group][0]];
			fdo << " " << phenotype_end[group_idx[i_group][0]];
			fdo << " " << (phenotype_neg[group_idx[i_group][0]]?"-":"+");
			switch (grp_mode) {
			case GRP_BEST: fdo << " " << phenotype_id[best_nominal_phenotype_abs] << " " << stb.str(group_size[i_group]); break;
			case GRP_PCA1: fdo << " " << stb.str(group_var[i_group], 3) << " " << stb.str(group_size[i_group]); break;
			case GRP_MEAN: fdo << " " << stb.str(group_size[i_group]); break;
			}
			fdo << " " << variant_indexes.size();
			fdo << " " << best_nominal_distance;
			fdo << " " << genotype_id[best_nominal_variant_abs];
			fdo << " " << genotype_chr[best_nominal_variant_abs];
			fdo << " " << genotype_start[best_nominal_variant_abs];
			fdo << " " << genotype_end[best_nominal_variant_abs];
			fdo << " " << dof_true;
			fdo << " " << dof_esti;
			fdo << " " << beta_ml1;
			fdo << " " << beta_ml2;
			fdo << " " << pval_nom;
			fdo << " " << pval_r2;
			fdo << " " << pval_slope;
			if(std_err) fdo << " " << pval_se;
			fdo << " " << pval_emp;
			fdo << " " << pval_bml;
			fdo << endl;
		}
	}
	fdo.close();
}
