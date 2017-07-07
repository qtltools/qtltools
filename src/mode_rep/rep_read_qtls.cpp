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

void rep_data::readQTLs(string fqtl) {
	string buffer; vector < string > tokens;
	qtl_count = 0;

	vrb.title("Reading QTL data in [" + fqtl + "]");
	input_file fdr(fqtl);
	if (fdr.fail()) vrb.error("Cannot not open file!");
	while (getline(fdr, buffer)) {
		stb.split(buffer, tokens);
		if (tokens.size() < 2) vrb.error("Incorrect number of columns, needs to be 2!");
		qtl_ids.first.push_back(tokens[0]);
		qtl_ids.second.push_back(tokens[1]);
		qtl_count++;
	}
	fdr.close();
	vrb.bullet("#qtls = " + stb.str(qtl_count));
	if (qtl_ids.first.size() == 0) vrb.error("No QTLs found!");
}
