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

#include "rtc_data.h"


void rtc_data::collapsePhenotypes() {
	group_idx.clear();

	//PASS0: check that groups are specified for aggragation methods



	//PASS1: regroup phenotypes by group ID
	//map < string, unsigned int > group_id;
	map < string, unsigned int >::iterator group_it;
	for (int p = 0 ; p < phenotype_count ; p ++) {
		if (grp_mode != GRP_NONE) {
			group_it = group_id.find(phenotype_grp[p]);
			if (group_it == group_id.end()) {
				group_idx.push_back(vector < unsigned int > (1, p));
				group_var.push_back(1.0);
				group_size.push_back(1);
				group_id.insert(pair < string, unsigned int > (phenotype_grp[p], group_id.size()));
			} else {
				group_idx[group_it->second].push_back(p);
				group_size[group_it->second]++;
			}
		} else {
			group_idx.push_back(vector < unsigned int > (1, p));
			group_var.push_back(1.0);
			group_size.push_back(1);
		}
	}

	//PASS2: sort & stats
	basic_stats bspg;
	for (int g = 0 ; g < group_idx.size() ; g ++) {
		sort(group_idx[g].begin(), group_idx[g].end());
		bspg.push(group_idx[g].size());
	}
	if (grp_mode != GRP_NONE) {
		vrb.title("Regrouping phenotypes within groups");
		vrb.bullet("#phenotypes = " + stb.str(phenotype_count));
		vrb.bullet("#groups = " + stb.str(group_idx.size()));
		vrb.bullet("#phenotypes per group = " + stb.str(bspg.mean(), 2) + " +/-" + stb.str(bspg.sd(), 2));
	}

	//PASS3: pca1 and mean
	basic_stats bsvg;
	for (int g = 0 ; g < group_idx.size() ; g ++) {
		if (group_idx[g].size() > 1) {
			if (grp_mode == GRP_MEAN) {
				for (int s = 0 ; s < sample_count ; s ++) {
					for (int p = 1 ; p < group_idx[g].size() ; p ++) phenotype_val[group_idx[g][0]][s] += phenotype_val[group_idx[g][p]][s];
					phenotype_val[group_idx[g][0]][s] /= group_idx[g].size();
				}
				group_idx[g].erase(group_idx[g].begin() + 1, group_idx[g].end());
			} else if (grp_mode == GRP_PCA1) {
				pca P (sample_count, group_idx[g].size());
				P.fill(phenotype_val, group_idx[g]);
				P.run(false, true, true);
				P.get(0, phenotype_val[group_idx[g][0]]);
				group_var[g] = P.getVariance(0);
				bsvg.push(group_var[g]);
				group_idx[g].erase(group_idx[g].begin() + 1, group_idx[g].end());
			}
		} else if (grp_mode == GRP_PCA1) group_var[g] = 1.0;
	}
	if (grp_mode == GRP_PCA1) vrb.bullet("variance explained by PC1 per group = " + stb.str(bsvg.mean(), 3) + " +/-" + stb.str(bsvg.sd(), 3));
}
