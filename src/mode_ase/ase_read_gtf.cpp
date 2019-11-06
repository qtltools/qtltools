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

void ase_data::readGTF(string fgtf){
    string buffer;
    vector < string > str;

    vrb.title("Reading exons in [" + fgtf + "]");
    input_file fd (fgtf);
    if (fd.fail()) vrb.error("Cannot open file!");
    int linecount = 0;
    int found_c = 0 , missed_c = 0;
    while(getline(fd, buffer)) {
        linecount++;
        if (linecount % 500000 == 0) vrb.bullet(stb.str(linecount) + " lines read");
        if (buffer[0] == '#') continue;
        stb.split(buffer, str, " \t;");
        if (str.size() < 10) vrb.error("Incorrect number of columns: " + stb.str(str.size()));
        if (str.size() % 2 ) vrb.error("Unmatched attributes: " + buffer);
        if (str[2] != "exon") continue;
        string chr = str[0];
		if(bam_chrs.count(chr) == 0){
			if(fix_chr && chr.size() > 3 && chr.substr(0,3) == "chr" && bam_chrs.count(chr.substr(3))){
				found_c++;
				chr = chr.substr(3);
			}else if (fix_chr && bam_chrs.count("chr" + chr)){
				found_c++;
				chr = "chr" + chr;
			}else{
				missed_c++;
			}
		}else{
			found_c++;
		}
        unsigned int start = atoi(str[3].c_str());
        unsigned int end = atoi(str[4].c_str());
        string strand = str[6];
        string gene_id="",gene_name="",trans_id="";
        for (int i = 8 ; i < str.size(); i+=2 ){
            str[i+1].erase(remove(str[i+1].begin(), str[i+1].end(), '"'), str[i+1].end());
            str[i+1].erase(remove(str[i+1].begin(), str[i+1].end(), ';'), str[i+1].end());
            if (str[i] == "gene_name") gene_name = str[i+1];
            if (str[i] == "gene_id") gene_id = str[i+1];
            if (str[i] == "transcript_id") trans_id = str[i+1];
        }
        if (gene_id=="") vrb.error("gene_id attribute is required: " + buffer);
        if (trans_id=="") vrb.error("transcript_id attribute is required: " + buffer);
        //cerr << gene_id << " " << gene_name << " " << chr << " "  << start << " " << end << " " << strand << " " << type << endl;
    	unsigned int sb = start / binsize;
    	unsigned int eb = end / binsize;
    	while (sb <= eb){
    		annotation[chr][sb].push_back(ase_exon(gene_id,trans_id,gene_name,start,end));
    		sb++;
    	}

    }

	if(found_c == 0) vrb.error("No chromosomes match between GTF and BAM. Try --fix-chr!");
	if(missed_c) vrb.warning(stb.str(missed_c) + " chromosomes are missing from the BAM file");
}

void ase_data::assignGenesToAseSite(ase_site &in){
	if (annotation.size() && annotation.count(in.chr)){
		unsigned int pos1based = in.pos + 1;
		unsigned int b = pos1based / binsize;
		if (annotation[in.chr].count(b)){
			string anno = "";
			for (int i =0 ; i < annotation[in.chr][b].size(); i++){
				if (annotation[in.chr][b][i].contains(pos1based)){
					if (anno != "") anno += ";";
					anno += annotation[in.chr][b][i].id;
				}
			}
			if (anno != "") in.genes = anno;
		}
	}
}
