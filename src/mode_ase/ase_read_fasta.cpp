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
const string MISSING_CHR = "___NA___";
const unsigned long BUFFER_SIZE = 250000000;

void ase_data::readGenome(string fin) {
	string buffer;
	vrb.title("Reading genome sequence in [" + fin + "]");
	input_file fd (fin);
	if (fd.fail()) vrb.error("Cannot open file!");
	string chr = MISSING_CHR;
	int found_c = 0 , missed_c = 0;
	bool read = false;
	while (getline(fd, buffer)) {
		if (buffer=="") continue;
		if (buffer[0] == '>'){
			buffer.erase(0,1);
			if(buffer == MISSING_CHR) vrb.error("Chromosome name cannot be " + MISSING_CHR);
			if (chr != MISSING_CHR) {
				if ((!bam_region.isSet() || bam_region.chr == chr) && read ) {
					vrb.print("    - " + chr + " " + stb.str(genome[chr].size()) + " bp read.");
					genome[chr].shrink_to_fit();
				}
			}
			chr = buffer;
			read = false;
			if (bam_region.isSet() && genome.size()) break;
			if(bam_chrs.count(chr) == 0){
				if(fix_chr && chr.size() > 3 && chr.substr(0,3) == "chr" && bam_chrs.count(chr.substr(3))){
					found_c++;
					chr = chr.substr(3);
					read = true;
				}else if (fix_chr && bam_chrs.count("chr" + chr)){
					found_c++;
					chr = "chr" + chr;
					read = true;
				}else{
					missed_c++;
				}
			}else{
				found_c++;
				read = true;
			}
			if ((!bam_region.isSet() || bam_region.chr == chr) && read ) genome[chr].reserve(BUFFER_SIZE);
		}else{
			if(chr == MISSING_CHR) vrb.error("Chromosome name missing for a sequence");
			if (!read || (bam_region.isSet()  && bam_region.chr != chr) ) continue;
			transform(buffer.begin(), buffer.end(),buffer.begin(), ::toupper);
			//copy(buffer.begin(), buffer.end(), back_inserter(genome[chr]));
			genome[chr] += buffer;
		}
	}

	if ((!bam_region.isSet() || bam_region.chr == chr) && read ) {
		vrb.print("    - " + chr + " "+ stb.str(genome[chr].size()) + " bp read.");
		genome[chr].shrink_to_fit();
	}

	if(found_c == 0) vrb.error("No chromosomes match between FASTA and BAM. Try --fix-chr!");
	if(missed_c) vrb.warning(stb.str(missed_c) + " FASTA chromosomes are missing from the BAM file. Found " + stb.str(found_c) + " chromosomes.");
}



