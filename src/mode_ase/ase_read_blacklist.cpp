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

#include "ase_data.h"

void ase_data::readBlacklist(string fgtf) {
	string buffer;
	vector < string > str;
	vector < ase_basic_block> input;
	unsigned long long mem = sizeof(vector < ase_basic_block>);
	unsigned long long max_step = 500000000;
	vrb.title("Reading blacklist in [" + fgtf + "]");
	input_file fd (fgtf);
	if (fd.fail()) vrb.error("Cannot open file!");
	int linecount = 0;
	set <string> found_c,missed_c;
	while(getline(fd, buffer)) {
		linecount++;
		if (linecount % 2000000 == 0) vrb.bullet(stb.str(linecount) + " lines read");
		if (buffer.size() == 0 || buffer[0] == '#') continue;
		stb.split(buffer, str);
		if (str.size() < 3) vrb.error("Incorrect number of columns: " + stb.str(str.size()));
        string chr = str[0];
		if(bam_chrs.count(chr) == 0){
			if(fix_chr && chr.size() > 3 && chr.substr(0,3) == "chr" && bam_chrs.count(chr.substr(3))){
				chr = chr.substr(3);
				found_c.insert(chr);
			}else if (fix_chr && bam_chrs.count("chr" + chr)){
				chr = "chr" + chr;
				found_c.insert(chr);
			}else{
				missed_c.insert(chr);
			}
		}else found_c.insert(chr);
        unsigned int start = atoi(str[1].c_str());
        unsigned int end = atoi(str[2].c_str());
        ase_basic_block abb(chr,start+1,end);
        if (bam_region.isSet() && !abb.overlap(bam_region)) continue;
        input.push_back(ase_basic_block(abb));
        mem += sizeof(ase_basic_block) + chr.capacity();
        if (on_the_fly && mem > max_step && input.size() > 1){
        	vector <ase_basic_block> temp;
        	sort(input.begin(),input.end());
        	ase_basic_block prev = input[0];
        	for (int i = 1 ; i < input.size(); i++){
        		if (prev.contiguous(input[i])) {
        			prev = prev.merge_nocheck(input[i]);
        		}else {
        			temp.push_back(prev);
        			prev = input[i];
        		}
        	}
        	temp.push_back(prev);
        	vrb.bullet("Reduced from " + stb.str(input.size()) + " to " + stb.str(temp.size()));
        	input = temp;
        	max_step += max_step;
        }

	}
	if(found_c.size() == 0) vrb.error("No chromosomes match between BED and BAM. Try --fix-chr!");
	if(missed_c.size()) vrb.warning(stb.str(missed_c.size()) + " BED chromosomes are missing from the BAM file. Found " + stb.str(found_c.size()) + " chromosomes.");
	blacklisted_regions.clear();
	if (input.size() > 1){
		sort(input.begin(),input.end());
		ase_basic_block prev = input[0];
		for (int i = 1 ; i < input.size(); i++){
			if (prev.contiguous(input[i])) {
				prev = prev.merge_nocheck(input[i]);
			}else {
				blacklisted_regions.push_back(prev);
				prev = input[i];
			}
		}
		blacklisted_regions.push_back(prev);
	}else blacklisted_regions = input;
	vrb.bullet(stb.str(blacklisted_regions.size()) + " non-overlapping regions read.");
}
