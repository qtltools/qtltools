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

void extract_data::writeOUT(string fout) {

	//
	string filename_header = fout + ".header.txt";
	vrb.title("Writing header in [" + filename_header + "]");
	output_file fdh(filename_header);
	fdh << "id chr start end" << endl;
	for (int v = 0 ; v < variable_id.size() ; v ++) {
		fdh << variable_id[v] << " " << variable_chr[v] << " " << variable_start[v] << " " << variable_end[v] << endl;
	}
	fdh.close();
	vrb.bullet("#variables = " + stb.str(variable_id.size()));

	string filename_content = fout + ".content.txt.gz";
	vrb.title("Writing content in [" + filename_content + "]");
	output_file fdc(filename_content);
	fdc << "sample";
	for (int v = 0 ; v < variable_id.size() ; v ++) fdc << " " << variable_id[v];
	fdc << endl;

	for (int s = 0 ; s < sample_count ; s ++) {
		fdc << sample_id[s];
		for (int v = 0 ; v < variable_id.size() ; v ++) fdc << " " << variable_val[v][s];
		fdc << endl;
	}
	fdc.close();
	vrb.bullet("#samples = " + stb.str(sample_count));
}
