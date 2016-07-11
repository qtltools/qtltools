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

#include "extract_data.h"

void extract_data::imputeMissing() {
	unsigned int n_missing = 0, n_nmissing = 0;
	vrb.title("Impute missing data with mean");
	for (int v = 0; v < variable_val.size() ; v ++) {
		double mean = 0.0; int c_mean= 0;
		for (int s = 0; s < sample_count; s ++) {
			if (variable_val[v][s] != "NA") {
				mean += atof(variable_val[v][s].c_str());
				c_mean ++;
			}
		}
		mean /= c_mean;
		for (int s = 0; s < sample_count ; s ++) {
			if (variable_val[v][s] == "NA") {
				variable_val[v][s] = stb.str(mean);
				n_missing ++;
			} else n_nmissing ++;
		}
	}
	vrb.bullet("#non_missing_data_points = " + stb.str(n_nmissing));
	vrb.bullet("#imputed_data_points = " + stb.str(n_missing) + " (=" + stb.str(n_missing * 100.0 / (n_nmissing + n_missing)) + "%)");
}
