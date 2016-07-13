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

#include "union_data.h"

void union_data::unions(string frtc1,int i){
	string buffer;
    vector < string > str;
    vrb.title("Reading QTLtools output in [" + frtc1 + "]");
    input_file fd (frtc1);
    if (fd.fail()) vrb.error("Cannot open file!");
    while(getline(fd, buffer)) {
        stb.split(buffer, str);
        if (str.size() < 2) vrb.error("Wrong fastQTL output file format");
        if (variant_column >= str.size()) vrb.error("variant column = " + stb.str(variant_column+1) + " but found " + stb.str(str.size()) + " columns in the following line:\n" + buffer);
        if (phenotype_column >= str.size()) vrb.error("phenotype column = " + stb.str(phenotype_column+1) + " but found " + stb.str(str.size()) + " columns in the following line:\n" + buffer);
        string pheno = str[phenotype_column];
        string snp = str[variant_column];
        if (!genotype_id_to_idx.count(snp) || !phenotype_id_to_idx.count(pheno) ) continue;
        int cs = getColdspot(genotype_chr[genotype_id_to_idx[snp]],genotype_start[genotype_id_to_idx[snp]]);
        if( cs < 0) continue;
        if (!toUnite.count(cs) || !toUnite[cs].count(pheno) ) {
        	//myPhenotype P(phenotype_chr[phenotype_id_to_idx[pheno]],phenotype_start[phenotype_id_to_idx[pheno]],no_of_files);
        	toUnite[cs][pheno] = myPhenotype(phenotype_chr[phenotype_id_to_idx[pheno]],phenotype_start[phenotype_id_to_idx[pheno]],no_of_files) ;
        }
        toUnite[cs][pheno].assign(snp,0,i);
    }
    fd.close();
}

void union_data::unions_conditional(string frtc1,int i){
	string buffer;
    vector < string > str;
    vrb.title("Reading QTLtools output in [" + frtc1 + "]");
    input_file fd (frtc1);
    if (fd.fail()) vrb.error("Cannot open file!");
    while(getline(fd, buffer)) {
        stb.split(buffer, str);
        if (str.size() < 4) vrb.error("Wrong fastQTL output file format");
        if (rank_column >= str.size()) vrb.error("rank column = " + stb.str(pvalue_column+1) + " but found " + stb.str(str.size()) + " columns in the following line:\n" + buffer);
		if (variant_column >= str.size()) vrb.error("variant column = " + stb.str(variant_column+1) + " but found " + stb.str(str.size()) + " columns in the following line:\n" + buffer);
		if (phenotype_column >= str.size()) vrb.error("phenotype column = " + stb.str(phenotype_column+1) + " but found " + stb.str(str.size()) + " columns in the following line:\n" + buffer);
		if (best_column >= str.size()) vrb.error("best column = " + stb.str(best_column+1) + " but found " + stb.str(str.size()) + " columns in the following line:\n" + buffer);
		string pheno = str[phenotype_column];
		string snp = str[variant_column];
		string best = str[best_column];
		int rank = atoi(str[rank_column].c_str());
		if (!genotype_id_to_idx.count(snp) || !phenotype_id_to_idx.count(pheno) || best != "1" ) continue;
        int cs = getColdspot(genotype_chr[genotype_id_to_idx[snp]],genotype_start[genotype_id_to_idx[snp]]);
        //cerr << cs << " " << pheno << " " << snp << " "<< frtc1 << " " << i << " "<< genotype_chr[genotype_id_to_idx[snp]] << " " << genotype_start[genotype_id_to_idx[snp]] << " " << genotype_id_to_idx[snp]<<endl;
        if( cs < 0) continue;
        if (!toUnite.count(cs) || !toUnite[cs].count(pheno) ) {
        	toUnite[cs][pheno] = myPhenotype(phenotype_chr[phenotype_id_to_idx[pheno]],phenotype_start[phenotype_id_to_idx[pheno]],no_of_files) ;
        }
        toUnite[cs][pheno].assign(snp,rank,i);
    }
    fd.close();
}

void union_data::create_unions(vector <string> &hitFiles , vector <string> &vcfFiles , vector <string> &bedFiles , vector <string> &covFiles){
	for (int i = 0 ; i < no_of_files; i++){
		readSampleFromBED(bedFiles[i]);
		if (vcfFiles.size() > 1) readSampleFromVCF(vcfFiles[i]);
		else readSampleFromVCF(vcfFiles[0]);
		if (options.count("cov")) readSampleFromCOV(covFiles[i]);
		mergeSampleLists();
		scanPhenotypes(bedFiles[i]);
		if(vcfFiles.size()==1){
			if(i==0) scanGenotypes(vcfFiles[i]);
		}else scanGenotypes(vcfFiles[i]);
		if(options.count("conditional")) unions_conditional(hitFiles[i],i);
		else unions(hitFiles[i],i);
		clearSamples();
	}
}

void union_data::find_unions(vector <string> &hitFiles , vector <string> &vcfFiles , vector <string> &bedFiles , vector <string> &covFiles){
	vrb.title("Find the best variant in all significant regions [" + stb.str(toUnite.size()) +"]");
	map < int ,map <string, myPhenotype > >::iterator outer_it;
	map <string, myPhenotype >::iterator inner_it;
	vector < results > results_vector(no_of_files);
	int outer_count = 1;
	for (outer_it = toUnite.begin(); outer_it != toUnite.end(); outer_it++){
		if (outer_count % 10 ==0) vrb.bullet("Region [" + stb.str(outer_count) + "/" + stb.str(toUnite.size()) +"]");
		outer_count++;
		int cs = outer_it->first;
		if (cs < 0) continue;
		vector < genotypes_holder > geno_sink(no_of_files);
		string regionG = all_coldspots_p[cs]->chr + ":" + stb.str(all_coldspots_p[cs]->start) + "-" + stb.str(all_coldspots_p[cs]->end);
		for (inner_it = outer_it->second.begin(); inner_it != outer_it->second.end(); inner_it++){
			string pheno = inner_it->first;
			string chr = inner_it->second.pheno_chr;
			int pos = inner_it->second.pheno_pos;
			string regionP = chr + ":" + stb.str(pos) + "-" + stb.str(pos+1);
			myPhenotype *P = &(inner_it->second);
			for (int i =0; i < no_of_files;i++){
				//cerr << bedFiles[i] << " " << pheno << " " << regionP << " " << regionG << endl;
				if (P->found[i]) {
					for (int u = 0 ; u < P->genotypes[i].size(); u++) results_vector[i].assign(pheno,P->genotypes[i][u],P->ranks[i][u],0.0,cs,regionG);
					for (int u = P->genotypes[i].size() ; u < P->max_independent_signal ; u++) results_vector[i].assign(pheno,"__UNION_FILLER_MAX_INDEP__",-1,0.0,cs,regionG);
				}else{
					readSampleFromBED(bedFiles[i],true);
					if (vcfFiles.size() > 1) readSampleFromVCF(vcfFiles[i],true);
					else readSampleFromVCF(vcfFiles[0],true);
					if (options.count("cov")) readSampleFromCOV(covFiles[i],true);
					mergeSampleLists(true);
					readPhenotypes(bedFiles[i],regionP);
                    if(!phenotype_id_to_idx.count(pheno)){
                    	for (int u = 0 ; u < P->max_independent_signal ; u++) results_vector[i].assign(pheno,"__UNION_FILLER_MISS_PHENO__",-1,0.0,cs,regionG);
                        clearSamples();
                        continue;
                    }
                    if (options.count("cov")){
                        readCovariates(covFiles[i]);
                        residualizePhenotypes();
                    }
					int idx = 0;
					if(vcfFiles.size()==1){
						if (!geno_sink[0].genotypes.size()){
							readGenotypes(vcfFiles[i],regionG);
							geno_sink[0].genotypes = genotype_val;
							geno_sink[0].ids = genotype_id;
						}
					}else{
						if (!geno_sink[i].genotypes.size()){
							readGenotypes(vcfFiles[i],regionG);
							geno_sink[i].genotypes = genotype_val;
							geno_sink[i].ids = genotype_id;
						}
						idx = i;
					}
					int geno_count = geno_sink[idx].genotypes.size();
					if (geno_count){
						double pval;
						string bestSNP = getBestVariant(geno_sink[idx],phenotype_id_to_idx[pheno], pval );
                        //if (bestSNP == "NA") cerr << pheno << " " << phenotype_id[phenotype_id_to_idx[pheno]] << " " << idx << " " << vcfFiles[idx] << " " << geno_sink[idx].genotypes.size() << endl;
						results_vector[i].assign(pheno,bestSNP, -1 , pval,cs,regionG);
						for (int u = 1 ; u < P->max_independent_signal ; u++) results_vector[i].assign(pheno,"__UNION_FILLER_MAX_INDEP__",-1,0.0,cs,regionG);
					}else{
						for (int u = 0 ; u < P->max_independent_signal ; u++) results_vector[i].assign(pheno,"__UNION_FILLER_MISS_GENO__",-1,0.0,cs,regionG);
					}
					clearSamples();
				}
			}
		}
        geno_sink.clear();
	}
	for (int i = 0 ; i < no_of_files; i++){
		string name = hitFiles[i].substr(hitFiles[i].find_last_of("/") + 1);
		//add prefix
		if (options.count("out-suffix")) name += options["out-suffix"].as <string> ();
		output_file fout(name + ".union");
		if (fout.fail()) vrb.error("Cannot open [" + name + ".union]");
		fout << results_vector[i];
	}

}

