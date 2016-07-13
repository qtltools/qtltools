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

#include "rtc_data.h"

void rtc_data::mergeqtl_trans_conditional(string frtc1, string frtc2){
    string buffer;
    vector < string > str;
    map < string , vector < int > > file1,pheno1;
    map < string , map < int , int > > rank1;
    map < string , vector < int >  >::iterator it;
    unsigned int count = 0;
    vrb.title("Reading conditional trans QTLtools output in [" + frtc1 + "]");
    input_file fd (frtc1);
    if (fd.fail()) vrb.error("Cannot open file!");
    string pheno,snp,best,group;
    int rank;
    while(getline(fd, buffer)) {
    	if(readRTCline(buffer,pheno,snp,best,group,rank)){
			file1[group].push_back(genotype_id_to_idx[snp]);
			rank1[group][genotype_id_to_idx[snp]] = rank;
			pheno1[group].push_back(phenotype_id_to_idx[pheno]);
    	}
    }
    fd.close();

    vrb.title("Reading conditional trans QTLtools output in [" + frtc2 + "]");
    input_file fd2 (frtc2);
    if (fd2.fail()) vrb.error("Cannot open file!");
    map < string , vector < int > > file2,pheno2;
    map < string , map < int , int > > rank2;
    while(getline(fd2, buffer)) {
    	if(readRTCline(buffer,pheno,snp,best,group,rank)){
    		if (file1.count(group)==0) continue;
			file2[group].push_back(genotype_id_to_idx[snp]);
			rank2[group][genotype_id_to_idx[snp]] = rank;
			pheno2[group].push_back(phenotype_id_to_idx[pheno]);
    	}
    }
    fd2.close();
    vrb.bullet(stb.str(file2.size()) + " common phenotypes found");
    vrb.title("Merging and calculating D' and R2");
    unsigned int event_count=0,pairs_tested =0;
    for (it = file2.begin(); it!= file2.end(); it++){
        event_count++;
        if(event_count % DprimePrintFreq == 0) {
        	vrb.bullet(stb.str(event_count) + " common phenotypes processed [" + stb.str(pairs_tested) + " pairs]");
        	pairs_tested =0;
        }
        vector <int > our = it->second;
        string pheno = it->first;
        vector < int > ourpi = pheno2[pheno];
        vector <int> their = file1[pheno];
        for (int o = 0 ; o < our.size(); o++){
            vector <int> others;
            for (int oo = 0 ; oo < our.size(); oo++) if (o != oo) others.push_back(our[oo]);
            int eqtl_snp_csi = getColdspot(genotype_chr[our[o]], genotype_start[our[o]]);
            for (int t = 0 ;t < their.size(); t++){
                if (genotype_chr[their[t]] == genotype_chr[our[o]] && abs(genotype_start[their[t]] - genotype_start[our[o]]) <= cis_window ){
                    int test_snp_csi = getColdspot(genotype_chr[their[t]], genotype_start[their[t]]);
                    vector <double> info= genotype_alleles.size() ? getDprimeRsquare(genotype_chr[their[t]], genotype_start[their[t]] , genotype_chr[our[o]], genotype_start[our[o]],genotype_alleles[their[t]], genotype_alleles[our[o]],their[t],our[o]) : getDprimeRsquare(genotype_chr[their[t]], genotype_start[their[t]] , genotype_chr[our[o]], genotype_start[our[o]],"","",their[t],our[o]);
                    pheno_eqtls[pheno].push_back(pairsToTestForRTC(their[t], our[o], others, test_snp_csi, eqtl_snp_csi, info[0] , info[1],rank1[pheno][their[t]],rank2[pheno][our[o]],ourpi[o]));
                    count++;
                    pairs_tested++;
                }
            }
        }
    }
    vrb.bullet(stb.str(count) + " potential merge events found for " + stb.str(pheno_eqtls.size()) + " phenotypes.");
    if (pheno_eqtls.size() == 0 ) vrb.leave("No merge events found!");
}

void rtc_data::mergeqtl_trans(string frtc1, string frtc2){
    string buffer;
    vector < string > str;
    map < string , vector < int > > file1,pheno1;
    map < string , vector < int >  >::iterator it;
    unsigned int count = 0;
    vrb.title("Reading trans QTLtools output in [" + frtc1 + "]");
    input_file fd (frtc1);
    if (fd.fail()) vrb.error("Cannot open file!");
    string pheno,snp,group;
    unsigned int line_count = 0;
    while(getline(fd, buffer)) {
    	if(readRTCline(buffer,pheno,snp,group,line_count)){
			file1[group].push_back(genotype_id_to_idx[snp]);
			pheno1[group].push_back(phenotype_id_to_idx[pheno]);
    	}
    }
    fd.close();
    line_count = 0;
    vrb.title("Reading trans QTLtools output in [" + frtc2 + "]");
    input_file fd2 (frtc2);
    if (fd2.fail()) vrb.error("Cannot open file!");
    map < string , vector < int > > file2,pheno2;
    while(getline(fd2, buffer)) {
    	if(readRTCline(buffer,pheno,snp,group,line_count)){
    		if (file1.count(group)==0) continue;
			file2[group].push_back(genotype_id_to_idx[snp]);
			pheno2[group].push_back(phenotype_id_to_idx[pheno]);
    	}
    }
    fd2.close();
    vrb.bullet(stb.str(file2.size()) + " common phenotypes found");
    vrb.title("Merging and calculating D' and R2");
    unsigned int event_count=0,pairs_tested =0;
    for (it = file2.begin(); it!= file2.end(); it++){
        event_count++;
        if(event_count % DprimePrintFreq == 0) {
        	vrb.bullet(stb.str(event_count) + " common phenotypes processed [" + stb.str(pairs_tested) + " pairs]");
        	pairs_tested =0;
        }
        vector <int > our = it->second;
        string pheno = it->first;
        vector < int > ourpi = pheno2[pheno];
        vector <int> their = file1[pheno];
        for (int o = 0 ; o < our.size(); o++){
            vector <int> others;
            for (int oo = 0 ; oo < our.size(); oo++) if (o != oo) others.push_back(our[oo]);
            int eqtl_snp_csi = getColdspot(genotype_chr[our[o]], genotype_start[our[o]]);
            for (int t = 0 ;t < their.size(); t++){
                if (genotype_chr[their[t]] == genotype_chr[our[o]] && abs(genotype_start[their[t]] - genotype_start[our[o]]) <= cis_window ){
                    int test_snp_csi = getColdspot(genotype_chr[their[t]], genotype_start[their[t]]);
                    vector <double> info= genotype_alleles.size() ? getDprimeRsquare(genotype_chr[their[t]], genotype_start[their[t]] , genotype_chr[our[o]], genotype_start[our[o]],genotype_alleles[their[t]], genotype_alleles[our[o]],their[t],our[o]) : getDprimeRsquare(genotype_chr[their[t]], genotype_start[their[t]] , genotype_chr[our[o]], genotype_start[our[o]],"","",their[t],our[o]);
                    pheno_eqtls[pheno].push_back(pairsToTestForRTC(their[t], our[o], others, test_snp_csi, eqtl_snp_csi, info[0] , info[1],0,0,ourpi[o]));
                    count++;
                    pairs_tested++;
                }
            }
        }
    }
    vrb.bullet(stb.str(count) + " potential merge events found for " + stb.str(pheno_eqtls.size()) + " phenotypes.");
    if (pheno_eqtls.size() == 0 ) vrb.leave("No merge events found!");

}

