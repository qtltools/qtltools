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

#include "gwas_data.h"

void gwas_data::readPhenotypes(string ftxt) {
	int n_includedP = 0;
	int n_excludedP = 0;
	string buffer;
	vector < int > mappingS;
	vector < string > tokens;

	//Open BED file
	vrb.title("Reading phenotype data in [" + ftxt + "]");
	input_file fd(ftxt);

	getline(fd, buffer);
	stb.split(buffer, tokens);
	if (tokens.size() < 2) vrb.error("Incorrect number of columns!");
	for (int t = 1 ; t < tokens.size() ; t ++) mappingS.push_back(findSample(tokens[t]));

    //Read phenotypes
	while (getline(fd, buffer)) {
		stb.split(buffer, tokens);
		if (tokens.size() < 2) vrb.error("Incorrect number of columns!");
		if (filter_phenotype.check(tokens[0])) {
			phenotype_id.push_back(tokens[0]);
			phenotype_val.push_back(vector < float > (sample_count, 0.0));
			for (int t = 1 ; t < tokens.size() ; t ++) {
				if (mappingS[t-1] >= 0) {
					if (tokens[t] == "NA") phenotype_val.back()[mappingS[t-1]] = bcf_float_missing;
					else phenotype_val.back()[mappingS[t-1]] = stof(tokens[t]);
				}
			}
			n_includedP++;
		} else n_excludedP++;
	}

	fd.close();
	phenotype_count = phenotype_id.size();
	vrb.bullet(stb.str(n_includedP) + " phenotypes included");
	if (n_excludedP > 0) vrb.bullet(stb.str(n_excludedP) + " phenotypes excluded by user");
    if (phenotype_count == 0) vrb.leave("Cannot find phenotypes!");
}
