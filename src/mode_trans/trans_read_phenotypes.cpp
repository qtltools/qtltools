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

#include "trans_data.h"

void trans_data::readPhenotypes(string fbed) {
    int n_includedS = 0, n_excludedS = 0, n_excludedU = 0, n_excludedP = 0, n_negativeStrd = 0;
    vector < int > mappingS;
    
    //Open BED file
    vrb.title("Reading phenotype data in [" + fbed + "] from line " + stb.str(start_line) + " to " + stb.str(end_line));
    htsFile *fp = hts_open(fbed.c_str(),"r");
    if (!fp) vrb.error("Cannot open file!");
    tbx_t * tbx = tbx_index_load(fbed.c_str());
    if (!tbx) vrb.error("Cannot open index file");
    kstring_t str = {0,0,0};
    if (hts_getline(fp, KS_SEP_LINE, &str) <= 0 || !str.l || str.s[0] != '#' ) vrb.error("Cannot read header line");
    
    //Process sample names
    vector < string > tokens;
    stb.split(string(str.s), tokens);
    if (tokens.size() < 7) vrb.error("Incorrect number of columns!");
    for (int t = 6 ; t < tokens.size() ; t ++) {
        if (filter_sample.check(tokens[t])) {
            mappingS.push_back(findSample(tokens[t]));
            if (mappingS.back() < 0) n_excludedS ++;
            else n_includedS ++;
        } else n_excludedU ++;
    }
    vrb.bullet(stb.str(n_includedS) + " samples included");
    if (n_excludedU > 0) vrb.bullet(stb.str(n_excludedU) + " samples excluded by user");
    if (n_excludedS > 0) vrb.bullet(stb.str(n_excludedS) + " samples without genotype data");
    if (n_includedS != sample_count) vrb.error("Cannot find phenotype data for " + stb.str(sample_count - n_includedS) + " samples!");
    
    unsigned long int linecount = 1;
    //Read phenotypes
    while (hts_getline(fp, KS_SEP_LINE, &str) >= 0) {
        if (str.l && str.s[0] != tbx->conf.meta_char) {
            stb.split(string(str.s), tokens);
            if (filter_phenotype.check(tokens[3])) {
                if (linecount >= start_line){
                    phenotype_id.push_back(tokens[3]);
                    phenotype_chr.push_back(tokens[0]);
                    phenotype_start.push_back(atoi(tokens[1].c_str()) + 1);
                    phenotype_end.push_back(atoi(tokens[2].c_str()));
                    phenotype_neg.push_back(tokens[5] == "-");
                    if (phenotype_neg.back()) n_negativeStrd ++;
                    phenotype_val.push_back(vector < float > (sample_count, 0.0));
                    for (int t = 6 ; t < tokens.size() ; t ++) if (mappingS[t-6] >= 0) {
                        if (tokens[t] == "NA") phenotype_val.back()[mappingS[t-6]] = bcf_float_missing;
                        else phenotype_val.back()[mappingS[t-6]] = stof(tokens[t]);
                    }
                }
                if (end_line != 0 && linecount >= end_line) break;
                linecount++;
            } else n_excludedP ++;
        }
    }
    
    //Finalize & verbose
    tbx_destroy(tbx);
    phenotype_count = phenotype_id.size();
    vrb.bullet(stb.str(phenotype_count) + " phenotypes included");
    if (n_excludedP > 0) vrb.bullet(stb.str(n_excludedP) + " phenotypes excluded by user");
    if (n_negativeStrd > 0) vrb.bullet(stb.str(n_negativeStrd) + " phenotypes is on the negative strand");
    hts_close(fp);
}


void trans_data::scanPhenotypes(string fbed) {
    int n_includedP = 0;
    int n_excludedP = 0;
    
    //Open BED file
    vrb.title("Reading phenotype data in [" + fbed + "]");
    htsFile *fp = hts_open(fbed.c_str(),"r");
    if (!fp) vrb.error("Cannot open file");
    tbx_t * tbx = tbx_index_load(fbed.c_str());
    if (!tbx) vrb.error("Cannot open index file");
    
    //Read header
    kstring_t str = {0,0,0};
    if (!hts_getline(fp, KS_SEP_LINE, &str) || !str.l || str.s[0] != tbx->conf.meta_char ) vrb.error("Cannot read header line");
    
    //Scan file
    vector < string > tokens;
    while (hts_getline(fp, KS_SEP_LINE, &str) >= 0) {
        if (str.l && str.s[0] != tbx->conf.meta_char) {
            stb.split(string(str.s), tokens," 	",4);
            if (tokens.size() < 4) vrb.error("Incorrect number of columns!");
            if (filter_phenotype.check(tokens[3])) {
                n_includedP++;
            } else n_excludedP ++;
        }
    }
    
    //Finalize & verbose
    tbx_destroy(tbx);
    if (hts_close(fp)) vrb.error("Cannot properly close file");
    phenotype_count = n_includedP;
    vrb.bullet(stb.str(n_includedP) + " phenotypes included");
    if (n_excludedP > 0) vrb.bullet(stb.str(n_excludedP) + " phenotypes excluded by user");
    if (phenotype_count == 0) vrb.leave("Cannot find phenotypes in region!");
}
