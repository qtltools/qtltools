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

#include "rep_data.h"

void rep_data::readPhenotypes(string fbed) {
	int n_includedS = 0;
	int n_includedP = 0;
	int n_excludedP = 0;
	int n_negativeStrd = 0;
	vector < int > mappingS;

	//Open BED file
	vrb.title("Reading phenotype data in [" + fbed + "]");
	htsFile *fp = hts_open(fbed.c_str(),"r");
	if (!fp) vrb.error("Cannot open file");
	tbx_t *tbx = tbx_index_load(fbed.c_str());
	if (!tbx) vrb.error("Cannot open index file");
	kstring_t str = {0,0,0};
	if (hts_getline(fp, KS_SEP_LINE, &str) <= 0 || !str.l || str.s[0] != tbx->conf.meta_char ) vrb.error("Cannot read header line!");

	//Process sample names
	vector < string > tokens;
	stb.split(string(str.s), tokens);
	if (tokens.size() < 7) vrb.error("Incorrect number of columns!");
	for (int t = 6 ; t < tokens.size() ; t ++) {
		mappingS.push_back(findSample(tokens[t]));
		if (mappingS.back() >= 0) n_includedS++;
	}

	//Read phenotypes
    unsigned int linecount =0;
    while (hts_getline(fp, KS_SEP_LINE, &str) >= 0) {
    	linecount ++;
		if (linecount % 100000 == 0) vrb.bullet("Read " + stb.str(linecount) + " lines");
		stb.split(string(str.s), tokens);
		if (str.l && str.s[0] != tbx->conf.meta_char) {
			if (tokens.size() < 7) vrb.error("Incorrect number of columns!");
			if (filter_phenotype.check(tokens[3])) {
				phenotype_id.push_back(tokens[3]);
				phenotype_chr.push_back(tokens[0]);
				phenotype_start.push_back(atoi(tokens[1].c_str()) + 1);
				phenotype_end.push_back(atoi(tokens[2].c_str()));
				phenotype_neg.push_back(tokens[5] == "-");
				if (phenotype_neg.back()) n_negativeStrd ++;
				phenotype_val.push_back(vector < float > (sample_count, 0.0));
				for (int t = 6 ; t < tokens.size() ; t ++) {
					if (mappingS[t-6] >= 0) {
						if (tokens[t] == "NA") phenotype_val.back()[mappingS[t-6]] = bcf_float_missing;
						else phenotype_val.back()[mappingS[t-6]] = stof(tokens[t]);
					}
				}
				n_includedP++;
			} else n_excludedP ++;
		}
	}

	//Finalize & verbose
	tbx_destroy(tbx);
	if (hts_close(fp)) vrb.error("Cannot properly close file");
	phenotype_count = phenotype_id.size();
	vrb.bullet(stb.str(n_includedP) + " phenotypes included");
	if (n_excludedP > 0) vrb.bullet(stb.str(n_excludedP) + " phenotypes excluded by user");
	if (n_negativeStrd > 0 ) vrb.bullet(stb.str(n_negativeStrd) + " phenotypes are on the negative strand");
    if (phenotype_count == 0) vrb.leave("Cannot find phenotypes in target region!");
}
