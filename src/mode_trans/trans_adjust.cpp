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

void trans_data::buildNullDistribution(string fnull) {
	string buffer;
	vector < string > tokens;
	vector < double > null_pvalues;

	//Open BED file
	vrb.title("Reading null p-values in [" + fnull + "]");
	input_file fd(fnull);
	if (fd.fail()) vrb.error("Cannot open file!");
	while (getline(fd, buffer)) {
		stb.split(buffer, tokens);
		null_pvalues.push_back(stof(tokens.back()));
	}
	fd.close();
	vrb.bullet("#null p-values = " + stb.str(null_pvalues.size()));

	double mean_null_pvalues = basic_stats(null_pvalues).mean();
	double variance_null_pvalues = basic_stats(null_pvalues).variance();
	vrb.bullet("Mean=" + stb.str(mean_null_pvalues) + " Var=" + stb.str(variance_null_pvalues));

	double beta_mm1 = mean_null_pvalues * (mean_null_pvalues * (1 - mean_null_pvalues ) / variance_null_pvalues - 1);
	double beta_mm2 = beta_mm1 * (1 / mean_null_pvalues - 1);
	beta_ml1 = beta_mm1;
	beta_ml2 = beta_mm2;
	vrb.bullet("Beta parameters (MM) : [s1=" + stb.str(beta_ml1) + ", s2=" + stb.str(beta_ml2) +"]");

	try {
		learnBetaParameters(null_pvalues, beta_ml1, beta_ml2);
	} catch (const std::exception & e) {
		vrb.bullet("Maximum Likelihood estimation failed, use Moment Matching instead!");
		beta_ml1 = beta_mm1;
		beta_ml2 = beta_mm2;
	}
	vrb.bullet("Beta parameters (ML): [s1=" + stb.str(beta_ml1) + ", s2=" + stb.str(beta_ml2) +"]");
}
