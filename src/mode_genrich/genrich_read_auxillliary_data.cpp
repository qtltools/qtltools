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

#include "genrich_data.h"

void genrich_data::readQTL(string fqtl) {
	string buffer;
	vector < string > tokens;
	unsigned int n_var_found = 0, n_var_unfound = 0;

	vrb.title("Reading QTLs in [" + fqtl + "]");
	input_file fd (fqtl);
	while (getline(fd, buffer)) {
		if (buffer[0] != '#') {
			stb.split(buffer, tokens);
			string uid = tokens[0] + "_" + tokens[1];
			unordered_map < string, unsigned int > :: iterator it_genotype_id = genotype_uuid.find(uid);

			if (it_genotype_id == genotype_uuid.end()) n_var_unfound ++;
			else {
				genotype_qtl[it_genotype_id->second] = true;
				n_var_found ++;
			}
		}
	}
	fd.close();
	vrb.bullet("Number of QTL found in reference = " + stb.str(n_var_found));
	vrb.bullet("Number of QTL not found in reference = " + stb.str(n_var_unfound));
}

void genrich_data::readGWAS(string fgwas) {
	string buffer;
	vector < string > tokens;
	unsigned int n_var_found = 0, n_var_unfound = 0;

	vrb.title("Reading GWAS hits in [" + fgwas + "]");
	input_file fd (fgwas);
	while (getline(fd, buffer)) {
		if (buffer[0] != '#') {
			stb.split(buffer, tokens);
			string uid = tokens[0] + "_" + tokens[1];
			unordered_map < string, unsigned int > :: iterator it_genotype_id = genotype_uuid.find(uid);

			if (it_genotype_id == genotype_uuid.end()) n_var_unfound ++;
			else {
				genotype_gwas[it_genotype_id->second] = true;
				n_var_found ++;
			}
		}
	}
	fd.close();
	vrb.bullet("Number of GWAS hits found in reference = " + stb.str(n_var_found));
	vrb.bullet("Number of GWAS hits not found in reference = " + stb.str(n_var_unfound));
}
