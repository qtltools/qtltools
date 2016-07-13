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

#include "fdensity_data.h"

void fdensity_data::readQTL(string fbed) {
	//Open BED file
	vrb.title("Reading QTL in [" + fbed + "]");
	htsFile *fp = hts_open(fbed.c_str(),"r");
	if (!fp) vrb.error("Cannot open file!");

	//Read annotations
	kstring_t str = {0,0,0};
	vector < string > tokens;
	while (hts_getline(fp, KS_SEP_LINE, &str) >= 0) {
		if (str.s[0] != '#') {
			stb.split(string(str.s), tokens);
			if (tokens.size() < 6) vrb.error("Incorrect number of columns!");
			tss_id.push_back(tokens[3]);
			tss_chr.push_back(tokens[0]);
			tss_pos.push_back(atoi(tokens[1].c_str()) + 1);
			tss_neg.push_back(tokens[5] == "-");
		}
	}

	//Finalize & verbose
	hts_close(fp);
	tss_count = tss_chr.size();
	vrb.bullet("#QTLs = " + stb.str(tss_chr.size()));
}
