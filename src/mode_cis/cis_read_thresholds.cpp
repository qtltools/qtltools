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

void cis_data::readThresholds(string fres) {
	string buffer; vector < string > tokens;

	//1.0 Allocation
	phenotype_threshold = vector < double > (phenotype_count, -1);
	vector < bool > phenotype_mask = vector < bool > (phenotype_count, false);

	//2.0 Read results
	vrb.title("Reading nominal thresholds in [" + fres + "]");
	input_file fdr(fres);
	if (fdr.fail()) vrb.error("Cannot not open file!");

	while (getline(fdr, buffer)) {
		stb.split(buffer, tokens);
		if (tokens.size() < 2) vrb.error("Incorrect number of columns!");

		vector < int> phenotype_idx;

		if (grp_mode != GRP_NONE) {
			for (int p = 0 ; p < phenotype_count ; p ++) if (phenotype_grp[p] == tokens[0]) phenotype_idx.push_back(p);
		} else {
			for (int p = 0 ; p < phenotype_count; p ++) if (phenotype_id[p] == tokens[0]) phenotype_idx.push_back(p);
		}

		for (int i = 0 ; i < phenotype_idx.size() ; i++) {
			if (tokens[1] != "NA") phenotype_threshold[phenotype_idx[i]] = atof(tokens[1].c_str());
			phenotype_mask[phenotype_idx[i]] = true;
		}
	}
	fdr.close();

	//3.0 Make sure that each MP has a qvalue
	int n_set= 0, n_unset = 0;
	for (int p = 0 ; p < phenotype_count ; p ++) {
		if (phenotype_mask[p]) n_set ++;
		else n_unset ++;
	}
	vrb.bullet("#phenotypes set = " + stb.str(n_set));
	if (n_unset > 0) vrb.error("Cannot find thresholds for " + stb.str(n_unset) + " phenotypes!");
}
