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

void cis_data::runConditionalPass(string fout) {

	//STEP0: INITIALIZATION OF IO
	output_file fdo (fout);
	if (fdo.fail()) vrb.error("Cannot open file [" + fout + "]");

	//STEP2: MAIN SWEEP THROUGH PHENOTYPES
	for (unsigned int i_group = 0 ; i_group < group_idx.size() ; i_group ++) {

		//STEP4: VERBOSE PROCESSED PHENOTYPES
		if (grp_mode == GRP_NONE) vrb.title("Processing phenotype [" + phenotype_id[group_idx[i_group][0]] + "] [" + stb.str(i_group+1) + "/" + stb.str(group_idx.size()) + "]");
		else {
			vrb.title("Processing group of phenotypes [" + phenotype_grp[group_idx[i_group][0]] + "] [" + stb.str(i_group+1) + "/" + stb.str(group_idx.size()) + "]");
			vrb.bullet("#phenotypes in group = " + stb.str(group_size[i_group]));
			if (grp_mode == GRP_PCA1) vrb.bullet("variance explained by PC1 = " + stb.str(group_var[i_group], 3));
		}

		//STEP4: ENUMERATE ALL VARIANTS IN CIS
		vector < unsigned int > variant_indexes;
		vector < int > variant_distances;
		for (unsigned int v = 0 ; v < genotype_count ; v ++) {
			if (phenotype_chr[group_idx[i_group][0]] != genotype_chr[v]) continue;
			int ps = (phenotype_start[group_idx[i_group][0]]>cis_window)?(phenotype_start[group_idx[i_group][0]]-cis_window):0;
			int pe = phenotype_end[group_idx[i_group][0]] + cis_window;

			if (genotype_start[v] <= pe && ps <= genotype_end[v]) {
				int cisdistance = 0;
				if (genotype_start[v] <= phenotype_end[group_idx[i_group][0]] && phenotype_start[group_idx[i_group][0]] <= genotype_end[v]) cisdistance = 0;
				else if (genotype_end[v] < phenotype_start[group_idx[i_group][0]]) cisdistance = (genotype_end[v] - phenotype_start[group_idx[i_group][0]]);
				else cisdistance = genotype_start[v] - phenotype_end[group_idx[i_group][0]];
				if (phenotype_neg[group_idx[i_group][0]]) cisdistance *= -1;
				variant_indexes.push_back(v);
				variant_distances.push_back(cisdistance);
			}
		}
		vrb.bullet("#variants in cis = " + stb.str(variant_indexes.size()));
		vrb.bullet("Nominal significance threshold = " + stb.str(phenotype_threshold[group_idx[i_group][0]]));

		//STEP5: VARIANTS IN CIS FOUND: PERFORM COMPUTATIONS
		if (variant_indexes.size() > 0) {

			//STEP6: FORWARD PASS
			bool fdone = true;
			unsigned int fhits = 0, fsignals = 0;
			vector < int > fbest_idx;
			vector < double > fbest_pvalue;
			vector < vector < double > > fpvalue, fslope;
			vector < vector < float > > phenotype_curr = vector < vector < float > > (group_idx[i_group].size());
			for (unsigned int p = 0 ; p < group_idx[i_group].size() ; p ++) phenotype_curr[p] = phenotype_val[group_idx[i_group][p]];

			do {
				if (fsignals > 0) {
					residualizer conditional_engine(sample_count);
					conditional_engine.push(genotype_val[fbest_idx.back()]);
					conditional_engine.build();
					for (unsigned int p = 0 ; p < group_idx[i_group].size() ; p ++) {
						conditional_engine.residualize(phenotype_curr[p]);
						normalTransform(phenotype_curr[p]);
					}
				}

				fdone = true;
				fbest_idx.push_back(-1);
				fbest_pvalue.push_back(1.0);
				fpvalue.push_back(vector < double > (group_idx[i_group].size() * variant_indexes.size(), 1.0));
				fslope.push_back(vector < double > (group_idx[i_group].size() * variant_indexes.size(), 0.0));

				for (unsigned int p = 0 ; p < group_idx[i_group].size() ; p ++) {
					for (unsigned int v = 0 ; v < variant_indexes.size() ; v ++) {
						unsigned int rel_idx = v * group_idx[i_group].size() + p;
						regression(genotype_val[variant_indexes[v]], phenotype_curr[p], fpvalue.back()[rel_idx], fslope.back()[rel_idx]);
						if (fpvalue.back()[rel_idx] <= phenotype_threshold[group_idx[i_group][p]]) {
							if (fdone) fsignals ++;
							fdone = false;
							fhits ++;
						}
						if (fpvalue.back()[rel_idx] <= fbest_pvalue.back()) {
							fbest_pvalue.back() = fpvalue.back()[rel_idx];
							fbest_idx.back() = variant_indexes[v];
						}
					}
				}
			} while (!fdone);

			fpvalue.pop_back();
			fslope.pop_back();
			fbest_idx.pop_back();
			fbest_pvalue.pop_back();

			vrb.bullet("Forward pass: [ni=" + stb.str(fsignals) + ", nh=" + stb.str(fhits) + "]");

			//STEP7: IF THERE IS SIGNIFICANT QTLs
			if (fsignals == 0) vrb.bullet("No backward pass");
			else {

				//STEP8: BACKWARD PASS
				bool bdone = true;
				unsigned int bhits = 0, bsignals = 0;
				vector < double > bbest_pvalue = vector < double > (fsignals, 1.0);
				vector < vector < double > > bpvalue = vector < vector < double > > (fsignals, vector < double > (group_idx[i_group].size() * variant_indexes.size(), 1.0));
				vector < vector < double > > bslope = vector < vector < double > > (fsignals, vector < double > (group_idx[i_group].size() * variant_indexes.size(), 0.0));

				for (unsigned int i_sig = 0 ; i_sig < fsignals ; i_sig ++) {

					//Dump phenotypes
					vector < vector < float > > phenotype_curr = vector < vector < float > > (group_idx[i_group].size());
					for (unsigned int p = 0 ; p < group_idx[i_group].size() ; p ++) phenotype_curr[p] = phenotype_val[group_idx[i_group][p]];

					//Iterative correction
					for (unsigned int s = 0 ; s < fsignals ; s ++) {
						if (s != i_sig) {
							residualizer conditional_engine (sample_count);
							conditional_engine.push(genotype_val[fbest_idx[s]]);
							conditional_engine.build();
							for (unsigned int p = 0 ; p < group_idx[i_group].size() ; p ++) {
								conditional_engine.residualize(phenotype_curr[p]);
								normalTransform(phenotype_curr[p]);
							}
						}
					}

					bdone = true;
					for (unsigned int p = 0 ; p < group_idx[i_group].size() ; p ++) {
						for (unsigned int v = 0 ; v < variant_indexes.size() ; v ++) {
							unsigned int rel_idx = v*group_idx[i_group].size() + p;
							regression(genotype_val[variant_indexes[v]], phenotype_curr[p], bpvalue[i_sig][rel_idx], bslope[i_sig][rel_idx]);
							if (bpvalue[i_sig][rel_idx] <= phenotype_threshold[group_idx[i_group][p]]) {
								if (bdone) bsignals ++;
								bdone = false;
								bhits ++;
							}
							if (bpvalue[i_sig][rel_idx] <= bbest_pvalue[i_sig]) bbest_pvalue[i_sig] = bpvalue[i_sig][rel_idx];
						}
					}
				}

				vrb.bullet("Backward pass: [ni=" + stb.str(bsignals) + ", nh=" + stb.str(bhits) + "]");

				//STEP9: WRITE OUTPUT
				for (unsigned int i_sig = 0 ; i_sig < fsignals ; i_sig ++) {
					vector < unsigned int > fbest_region, bbest_region;
					for (unsigned int p = 0 ; p < group_idx[i_group].size() ; p ++) {
						for (unsigned int v = 0 ; v < variant_indexes.size() ; v ++) {
							unsigned int rel_idx = v * group_idx[i_group].size() + p;
							if (fpvalue[i_sig][rel_idx] == fbest_pvalue[i_sig]) fbest_region.push_back(rel_idx);
							if (bpvalue[i_sig][rel_idx] == bbest_pvalue[i_sig]) bbest_region.push_back(rel_idx);
						}
					}
					random_shuffle(fbest_region.begin(), fbest_region.begin());
					random_shuffle(bbest_region.begin(), bbest_region.begin());
					for (unsigned int p = 0 ; p < group_idx[i_group].size() ; p ++) {
						for (unsigned int v = 0 ; v < variant_indexes.size() ; v ++) {
							unsigned int rel_idx = v * group_idx[i_group].size() + p;
							if (fpvalue[i_sig][rel_idx] <= phenotype_threshold[group_idx[i_group][p]] || bpvalue[i_sig][rel_idx] <= phenotype_threshold[group_idx[i_group][p]]) {
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
								fdo << " " << i_sig;
								fdo << " " << fpvalue[i_sig][rel_idx];
								fdo << " " << fslope[i_sig][rel_idx];
								fdo << " " << ((rel_idx == fbest_region[0])?"1":"0");
								fdo << " " << ((fpvalue[i_sig][rel_idx] <= phenotype_threshold[group_idx[i_group][p]])?"1":"0");
								fdo << " " << bpvalue[i_sig][rel_idx];
								fdo << " " << bslope[i_sig][rel_idx];
								fdo << " " << ((rel_idx == bbest_region[0])?"1":"0");
								fdo << " " << ((bpvalue[i_sig][rel_idx] <= phenotype_threshold[group_idx[i_group][p]])?"1":"0");
								fdo << endl;
							}
						}
					}
				}
			}
		}
	}
	fdo.close();
}
