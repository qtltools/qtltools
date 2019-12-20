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
	timer current_timer;
	string buffer;
	vector < string > str;
	vector < basic_block> input;
	bool already_sorted = true, need_merging = false;
	unsigned long long max_step = (500 * 1024 * 1024) / (sizeof(basic_block)+16);
	unsigned long long next_step = max_step;
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
        basic_block abb(chr,start+1,end);
        if (bam_region.isSet() && !abb.overlap(bam_region)) continue;
        if (input.size()){
        	if(abb < input.back()){ already_sorted = false; need_merging = true;}
        	else if (abb.contiguous(input.back())) need_merging = true;
        }
        input.push_back(basic_block(abb));
        if (on_the_fly && input.size() >= next_step ){
        	if (need_merging){
				size_t ps = blacklisted_regions.size();
				mergeContiguousBlocks(input, !already_sorted);
				//mergeContiguousBlocks(input, blacklisted_regions);
				vrb.bullet("Reduced from " + stb.str(input.size() + ps) + " to " + stb.str(blacklisted_regions.size()));
				//input = blacklisted_regions;
				input = vector <basic_block>(0);
				already_sorted = true;
				need_merging = false;
        	}else{
        		next_step += max_step;
        	}

        }

	}
	if(found_c.size() == 0) vrb.error("No chromosomes match between BED and BAM. Try --fix-chr!");
	if(missed_c.size()) vrb.warning(stb.str(missed_c.size()) + " BED chromosomes are missing from the BAM file. Found " + stb.str(found_c.size()) + " chromosomes.");
	//if (need_merging) mergeContiguousBlocks(input,blacklisted_regions);
	if (need_merging) mergeContiguousBlocks(input,!already_sorted);
	else blacklisted_regions.insert(blacklisted_regions.end(),input.begin(), input.end());
	vrb.bullet(stb.str(blacklisted_regions.size()) + " non-overlapping regions read.");
	vrb.bullet("Time taken: " + stb.str(current_timer.high_res_abs_time()) + " seconds");
}
