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

void cis_data::runNominalPass(string fout) {

	//STEP0: INITIALIZATION OF IO
	output_file fdo (fout);
	if (fdo.fail()) vrb.error("Cannot open file [" + fout + "]");

	//STEP2: INITIALIZE A WORKING COPY OF GENOTYPES
	vector < double > genotype_sd = vector < double > (genotype_count, 0.0);
	for (unsigned int v = 0 ; v < genotype_count ; v ++) {
		genotype_sd[v] = basic_stats(genotype_val[v]).sd();
		normalize(genotype_val[v]);
	}

	//STEP3: INITIALIZE A WORKING COPY OF PHENOTYPES
	vector < double > phenotype_sd = vector < double > (phenotype_count, 0.0);
	for (unsigned int p = 0 ; p < phenotype_count ; p ++) {
		phenotype_sd[p] = basic_stats(phenotype_val[p]).sd();
		normalize(phenotype_val[p]);
	}

	//STEP4: MAIN SWEEP THROUGH PHENOTYPES
	for (unsigned int i_group = 0 ; i_group < group_idx.size() ; i_group ++) {

		//STEP4: VERBOSE PROCESSED PHENOTYPES
		if (grp_mode == GRP_NONE) vrb.title("Processing phenotype [" + phenotype_id[group_idx[i_group][0]] + "] [" + stb.str(i_group+1) + "/" + stb.str(group_idx.size()) + "]");
		else {
			vrb.title("Processing group of phenotypes [" + phenotype_grp[group_idx[i_group][0]] + "] [" + stb.str(i_group+1) + "/" + stb.str(group_idx.size()) + "]");
			vrb.bullet("#phenotypes in group = " + stb.str(group_size[i_group]));
			if (grp_mode == GRP_PCA1) vrb.bullet("variance explained by PC1 = " + stb.str(group_var[i_group], 3));
		}

		//STEP6: ENUMERATE ALL VARIANTS IN CIS
		vector < unsigned int > variant_indexes;
		vector < int > variant_distances;
		for (unsigned int v = 0 ; v < genotype_count ; v ++) {
            if (phenotype_chr[group_idx[i_group][0]] != genotype_chr[v]) continue;
			int ps = (phenotype_start[group_idx[i_group][0]]>cis_window)?(phenotype_start[group_idx[i_group][0]]-cis_window):0;
			int pe = phenotype_end[group_idx[i_group][0]] + cis_window;

			if (genotype_start[v] <= pe && ps <= genotype_end[v]) {
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

		//STEP7: VARIANTS IN CIS FOUND: FULL COMPUTATIONS
		double dof_true;
		if (variant_indexes.size() > 0) {
			double best_nominal_correlation = 0.0;
			int best_nominal_variant_abs = -1;
			int best_nominal_variant_rel = -1;
			int best_nominal_phenotype_abs = -1;
			int best_nominal_phenotype_rel = -1;
			int best_nominal_distance = cis_window;
			vector < double > pval_nom = vector < double > (variant_indexes.size() * group_idx[i_group].size(), 0.0);
			vector < double > pval_slope = vector < double > (variant_indexes.size() * group_idx[i_group].size(), 0.0);
			vector < double > pval_r2 = vector < double > (variant_indexes.size() * group_idx[i_group].size(), 0.0);
			vector < double > pval_se ;
			if(std_err) pval_se = vector < double > (variant_indexes.size() * group_idx[i_group].size(), 0.0);

			//STEP8: ASSOCIATION TESTING
			dof_true = sample_count - 2;
			for (unsigned int p = 0 ; p < group_idx[i_group].size() ; p ++) {
				for (unsigned int v = 0 ; v < variant_indexes.size() ; v ++) {
					double curr_correlation = getCorrelation(genotype_val[variant_indexes[v]], phenotype_val[group_idx[i_group][p]]);
					double r2 = curr_correlation * curr_correlation;
					pval_r2[v*group_idx[i_group].size()+p] = r2;
					pval_nom[v*group_idx[i_group].size()+p] = getPvalue(curr_correlation, dof_true);
					pval_slope[v*group_idx[i_group].size()+p] = getSlope(curr_correlation, genotype_sd[variant_indexes[v]], phenotype_sd[group_idx[i_group][p]]);
					//cerr << genotype_id[variant_indexes[v]] << " " << phenotype_id[group_idx[i_group][p]] << " ";
					if (std_err) pval_se[v*group_idx[i_group].size()+p] = getSE(r2, genotype_sd[variant_indexes[v]], phenotype_sd[group_idx[i_group][p]]);
					if (abs(curr_correlation) > abs(best_nominal_correlation) || (abs(curr_correlation) == abs(best_nominal_correlation) && abs(variant_distances[v]) < abs(best_nominal_distance))) {
						best_nominal_correlation = curr_correlation;
						best_nominal_variant_rel = v;
						best_nominal_variant_abs = variant_indexes[v];
						best_nominal_phenotype_rel = p;
						best_nominal_phenotype_abs = group_idx[i_group][p];
						best_nominal_distance = variant_distances[v];
					}
				}
			}

			//STEP9: VERBOSE BEST HIT
			if (grp_mode == GRP_BEST)
				vrb.bullet("Best hit: [id=" + genotype_id[best_nominal_variant_abs] + ", d=" + stb.str(best_nominal_distance) + ", p=" + phenotype_id[best_nominal_phenotype_abs] + ", pv=" + stb.str(pval_nom[best_nominal_variant_rel]) + ", s=" + stb.str(pval_slope[best_nominal_variant_rel], 4) + "]");
			else vrb.bullet("Best hit: [id=" + genotype_id[best_nominal_variant_abs] + ", d=" + stb.str(best_nominal_distance) + ", pv=" + stb.str(pval_nom[best_nominal_variant_rel]) + ", s=" + stb.str(pval_slope[best_nominal_variant_rel], 4) + "]");

			//STEP10: PRINT RESULTS IN FILE
			for (unsigned int p = 0 ; p < group_idx[i_group].size() ; p ++) {
				for (unsigned int v = 0 ; v < variant_indexes.size() ; v ++) {
					bool toBeWritten = false;
					if ((phenotype_threshold.size() == 0) && (pval_nom[v * group_idx[i_group].size() + p] <= threshold)) toBeWritten = true;
					if ((phenotype_threshold.size() != 0) && (pval_nom[v * group_idx[i_group].size() + p] <= phenotype_threshold[group_idx[i_group][p]])) toBeWritten = true;
					if (toBeWritten) {
						if (grp_mode == GRP_NONE) fdo << phenotype_id[group_idx[i_group][p]];
						else fdo << phenotype_grp[group_idx[i_group][p]];
						fdo << " " << phenotype_chr[group_idx[i_group][p]];
						fdo << " " << phenotype_start[group_idx[i_group][p]];
						fdo << " " << phenotype_end[group_idx[i_group][p]];
						fdo << " " << (phenotype_neg[group_idx[i_group][p]]?"-":"+");
						switch (grp_mode) {
						case GRP_BEST: fdo << " " << phenotype_id[group_idx[i_group][p]] << " " << stb.str(group_size[i_group]); break;
						case GRP_PCA1: fdo << " " << stb.str(group_var[i_group], 3) << " " << stb.str(group_size[i_group]); break;
						case GRP_MEAN: fdo << " " << stb.str(group_size[i_group]); break;
						}
						fdo << " " << variant_indexes.size();
						fdo << " " << variant_distances[v];
						fdo << " " << genotype_id[variant_indexes[v]];
						fdo << " " << genotype_chr[variant_indexes[v]];
						fdo << " " << genotype_start[variant_indexes[v]];
						fdo << " " << genotype_end[variant_indexes[v]];
						fdo << " " << pval_nom[v * group_idx[i_group].size() + p];
						fdo << " " << pval_r2[v * group_idx[i_group].size() + p];
						fdo << " " << pval_slope[v * group_idx[i_group].size() + p];
						if (std_err) fdo << " " << pval_se[v * group_idx[i_group].size() + p];
						fdo << " " << ((v == best_nominal_variant_rel && p == best_nominal_phenotype_rel)?"1":"0");
						fdo << endl;
					}
				}
			}
		}
	}
	fdo.close();
}
