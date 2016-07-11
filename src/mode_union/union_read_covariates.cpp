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

#include "union_data.h"

void union_data::readCovariates(string fcov) {
    covariate_val.clear();
    covariate_count =0;
    string buffer;
	vector < string > str;
	int n_includedS = 0;
	int n_includedC = 0;
	int n_excludedC = 0;
	vector < int > mappingS;

	//vrb.title("Reading covariates in [" + fcov + "]");
	input_file fd (fcov);
	if (fd.fail()) vrb.error("Cannot open file!");

	//Read samples
	getline(fd, buffer);
	if (buffer.size() == 0) vrb.error("No header line detected!");
	stb.split(buffer, str	);
	for (int t = 1 ; t < str.size() ; t ++) {
		mappingS.push_back(findSample(str[t]));
		if (mappingS.back() >= 0) n_includedS++;
	}

	//Read covariates
	while(getline(fd, buffer)) {
		stb.split(buffer, str);
		if (str.size() < 2) vrb.error("Incorrect number of columns!");
		if (filter_covariate.check(str[0])) {
			covariate_val.push_back(vector < string > (sample_count));
			for (int t = 1 ; t < str.size() ; t ++) if (mappingS[t-1] >= 0) covariate_val.back()[mappingS[t-1]] = str[t];
            n_includedC ++;
		} else n_excludedC ++;
	}

	//Finalise
	covariate_count = n_includedC;
	//vrb.bullet(stb.str(n_includedC) + " covariates included");
	//if (n_excludedC > 0) vrb.bullet(stb.str(n_excludedC) + " covariates excluded");
	fd.close();
}
