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

void trans_data::readCovariates(string fcov) {
	string buffer;vector < string > tokens;
	int n_includedS = 0;
	int n_excludedS = 0;
	int n_missingS = 0;
	int n_includedC = 0;
	int n_excludedC = 0;
	vector < int > mappingS;

	vrb.title("Reading covariates in [" + fcov + "]");
	input_file fd (fcov);
	if (fd.fail()) vrb.error("Cannot open file!");

	//Read samples
	getline(fd, buffer);
	stb.split(buffer, tokens);
	for (int t = 1 ; t < tokens.size() ; t ++) {
		if (filter_sample.check(tokens[t])) {
			mappingS.push_back(findSample(tokens[t]));
			if (mappingS.back() >= 0) n_includedS ++;
			else n_missingS ++;
		} else {
			mappingS.push_back(-1);
			n_excludedS++;
		}
	}
	vrb.bullet(stb.str(n_includedS) + " samples included");
	if (n_excludedS > 0) vrb.bullet(stb.str(n_excludedS) + " samples excluded by user");
	if (n_missingS > 0) vrb.bullet(stb.str(n_missingS) + " samples without phenotype data");
	if (n_includedS != sample_count) vrb.error("Cannot find covariates for " + stb.str(sample_count - n_includedS) + " samples!");

	//Read covariates
	while(getline(fd, buffer)) {
        if (buffer.size() == 0) continue;
		stb.split(buffer, tokens);
		if (tokens.size() < 2) vrb.error("Incorrect number of columns!");
		if (filter_covariate.check(tokens[0])) {
			covariate_val.push_back(vector < string > (sample_count));
			for (int t = 1 ; t < tokens.size() ; t ++) if (mappingS[t-1] >= 0) covariate_val.back()[mappingS[t-1]] = tokens[t];
			n_includedC ++;
		} else n_excludedC ++;
	}

	//Finalise
	covariate_count = n_includedC;
	vrb.bullet(stb.str(n_includedC) + " covariate(s) included");
	if (n_excludedC > 0) vrb.bullet(stb.str(n_excludedC) + " covariate(s) excluded");
	fd.close();
}

