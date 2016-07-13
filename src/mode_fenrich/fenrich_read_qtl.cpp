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

#include "fenrich_data.h"

void fenrich_data::readQTL(string fqtl) {
	string buffer; vector < string > str;

	//Read QTL
	vrb.title("Reading QTL in [" + fqtl + "]");
	input_file fdq(fqtl);
	if (fdq.fail()) vrb.error("Cannot open file!");
	while (getline(fdq, buffer)) {
		if (buffer[0] != '#') {
			stb.split(buffer, str);
			if (str.size() < 6) vrb.error("Incorrect number of columns, observed = " + stb.str(str.size())  + " expected = 5");
			int idx_tss = findTSS(str[4]);
			if (idx_tss < 0) vrb.error("Unknown phenotype id!");
			qtl_pos.push_back(atoi(str[1].c_str()) - tss_pos[idx_tss]);
			qtl_order.push_back(idx_tss);
		}
	}
	fdq.close();
	qtl_count = qtl_pos.size();
	vrb.bullet("#qtl = " + stb.str(qtl_count));
}
