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

void rtc_data::createTransLists(){
    vrb.title("Creating a list of variants and phenotypes for trans analyses");
    unsigned int count = 0;
    map < int ,vector < pairsToTestForRTC > >::iterator it;
    for (it = pheno_eqtls.begin(); it != pheno_eqtls.end(); it++){
        filter_phenotype.addInclusion(phenotype_id[it->first]);
        for (int p = 0 ; p < it->second.size(); p++){
            vector < int > genotype_idx_to_test;
            if (it->second[p].eqtl_snp_coldspot_idx >= 0 && it->second[p].test_snp_coldspot_idx >= 0 && genotype_chr[it->second[p].test_snp_idx] == genotype_chr[it->second[p].eqtl_snp_idx] && (it->second[p].eqtl_snp_coldspot_idx == it->second[p].test_snp_coldspot_idx || it->second[p].Dprime >= Dprime_cutoff)){
                int si,ei;
                if (it->second[p].eqtl_snp_coldspot_idx < it->second[p].test_snp_coldspot_idx){
                    si = it->second[p].eqtl_snp_coldspot_idx;
                    ei = it->second[p].test_snp_coldspot_idx;
                }else{
                    ei = it->second[p].eqtl_snp_coldspot_idx;
                    si = it->second[p].test_snp_coldspot_idx;
                }
                for (int csi = si ; csi <= ei; csi++){
                    genotype_idx_to_test.insert(genotype_idx_to_test.end(),all_coldspots_p[csi]->coldspot_variant_idx.begin(),all_coldspots_p[csi]->coldspot_variant_idx.end());
                }
            } else if (it->second[p].eqtl_snp_coldspot_idx == -1 && it->second[p].test_snp_coldspot_idx == -1 && genotype_chr[it->second[p].test_snp_idx] == genotype_chr[it->second[p].eqtl_snp_idx]){
                genotype_idx_to_test=coldspot_end_idx[genotype_chr[it->second[p].test_snp_idx]];
            } else if (it->second[p].eqtl_snp_coldspot_idx >= -1 && it->second[p].test_snp_coldspot_idx == -1 && genotype_chr[it->second[p].test_snp_idx] == genotype_chr[it->second[p].eqtl_snp_idx] && it->second[p].Dprime >= Dprime_cutoff){
                string chr = genotype_chr[it->second[p].eqtl_snp_idx];
                int start_idx = it->second[p].eqtl_snp_coldspot_idx;
                for (int csi = start_idx ; csi < all_coldspots_p.size() && all_coldspots_p[csi]->chr == chr; csi++){
                    genotype_idx_to_test.insert(genotype_idx_to_test.end(),all_coldspots_p[csi]->coldspot_variant_idx.begin(),all_coldspots_p[csi]->coldspot_variant_idx.end());
                }
                genotype_idx_to_test.insert(genotype_idx_to_test.end(),coldspot_end_idx[chr].begin(),coldspot_end_idx[chr].end());
            } else if (it->second[p].eqtl_snp_coldspot_idx == -1 && it->second[p].test_snp_coldspot_idx >= -1 && genotype_chr[it->second[p].test_snp_idx] == genotype_chr[it->second[p].eqtl_snp_idx] && it->second[p].Dprime >= Dprime_cutoff){
                string chr = genotype_chr[it->second[p].test_snp_idx];
                int start_idx = it->second[p].test_snp_coldspot_idx;
                for (int csi = start_idx ; csi < all_coldspots_p.size() && all_coldspots_p[csi]->chr == chr; csi++){
                    genotype_idx_to_test.insert(genotype_idx_to_test.end(),all_coldspots_p[csi]->coldspot_variant_idx.begin(),all_coldspots_p[csi]->coldspot_variant_idx.end());
                }
                genotype_idx_to_test.insert(genotype_idx_to_test.end(),coldspot_end_idx[chr].begin(),coldspot_end_idx[chr].end());
            }else{
                continue;
            }
            count++;
            for (int i = 0 ; i < genotype_idx_to_test.size(); i ++) filter_genotype.addInclusion(genotype_id[genotype_idx_to_test[i]]);
        }
    }
    vrb.bullet(stb.str(count) + " actual RTC calculations found.");
}


vector <double> rtc_data::getDprimeRsquare(string chr1, int pos1, string chr2 , int pos2, string al1, string al2, int idx1 , int idx2){
    vector <double > result(2,-9.0);

    if(!calculate_Dprime_R2){
    	if ( idx1 >= 0 && idx2 >= 0) result[1] = getRsquare(idx1,idx2);
        return result;
    }

    if (chr1 == chr2 && pos1 == pos2){
        result = {1.0,1.0};
        return result;
    }
    vector < string > geno1,geno2;
    string alleles1,alleles2;
    string region1 = chr1 + ":" + stb.str(pos1) + "-" + stb.str(pos1);
    string region2 = chr2 + ":" + stb.str(pos2) + "-" + stb.str(pos2);
    if (DprimeR2inMem){
        if (DprimeRsquareSink.count(region1) && DprimeRsquareSink[region1].count(region2) ) return DprimeRsquareSink[region1][region2];
        if (DprimeRsquareSink.count(region2) && DprimeRsquareSink[region1].count(region1) ) return DprimeRsquareSink[region2][region1];
    }
    int s1 = readGenotypesVCFStats(region1, alleles1, geno1);
    int s2 = readGenotypesVCFStats(region2, alleles2, geno2);
    //if (!s1) vrb.warning(region1 + " does not exist in [" +stats_vcf_file +"]");
    //if (!s2) vrb.warning(region2 + " does not exist in [" +stats_vcf_file +"]");
    if (s1 == 0) unfound_regions.insert(region1);
    if (s2 == 0) unfound_regions.insert(region2);
    if (s1 == -1) unphased_regions.insert(region1);
    if (s2 == -1) unphased_regions.insert(region2);
    if ((s1 == -1 || s2 == -1 )){
    	if(!options.count("individual-Dprime")){
			calculate_Dprime_R2 = false;
			vrb.warning("Unphased genotypes encounted stopping D' calculations");
    	}
    	if ( idx1 >= 0 && idx2 >= 0) result[1] = getRsquare(idx1,idx2);
    }
    if (s1 <= 0 || s2 <= 0){
    	if ( idx1 >= 0 && idx2 >= 0) result[1] = getRsquare(idx1,idx2);
        if (DprimeR2inMem) DprimeRsquareSink[region1][region2] = result;
        return result;
    }
    //CALCULATE D' AND R2 TAKEN FROM VCFTOOLS variant_file_output.cpp
    double x11=0, x12=0, x21=0, x22=0;
    double X=0, X2=0, Y=0, Y2=0, XY=0;
    double sx, sy;
    double rel_x11, p1, p2, q1, q2, Dmax;
    double var1, var2, cov12;
    double r2,D, Dprime;
    int chr_count = 0;
    int allele1, allele2;
    for ( int i = 0 ; i < geno1.size(); i++){
        if (geno1[i] == "NA" || geno2[i] == "NA") continue;
        int g1a1 = geno1[i][0] - '0';
        int g1a2 = geno1[i][1] - '0';
        int g2a1 = geno2[i][0] - '0';
        int g2a2 = geno2[i][1] - '0';
        for (unsigned int c=0; c<2; c++){
            if (c==0){
                allele1 = g1a1;
                allele2 = g2a1;
            }else{
                allele1 = g1a2;
                allele2 = g2a2;
            }

            if ((allele1 < 0) || (allele2 < 0))
                continue;

            if (allele1 == 0 && allele2 == 0){
                x11++;
            } else if (allele1 == 0 && allele2 != 0){
                x12++;
            } else if (allele1 != 0 && allele2 == 0){
                x21++;
            } else { // (allele1 !=0 && allele2 != 0)
                x22++;
            }

            sx=0, sy=0;
            if (allele1 == 0)
                sx += 1;

            if (allele2 == 0)
                sy += 1;

            X += sx; Y += sy;
            XY += sx*sy;
            sx *= sx; sy *= sy;
            X2 += sx;
            Y2 += sy;

            chr_count++;
        }

    }
    rel_x11 = x11/double(chr_count);
    p1 = (x11 + x12)/double(chr_count);
    p2 = (x21 + x22)/double(chr_count);
    q1 = (x11 + x21)/double(chr_count);
    q2 = (x12 + x22)/double(chr_count);
    D = rel_x11 - p1*q1;
    if (D < 0)
        Dmax = min(p1*q1,p2*q2);
    else
        Dmax = min(p1*q2,p2*q1);
    Dprime = D/Dmax;

    X /= chr_count; X2 /= chr_count;
    Y /= chr_count; Y2 /= chr_count;
    XY /= chr_count;
    var1 = X2 - X*X;
    var2 = Y2 - Y*Y;
    cov12 = XY - X*Y;
    //if (var1 == 0) vrb.warning(region1 + " has zero variance in [" +stats_vcf_file +"]");
    //if (var2 == 0) vrb.warning(region2 + " has zero variance in [" +stats_vcf_file +"]");
    if (var1 == 0) no_variance_regions.insert(region1);
    if (var2 == 0) no_variance_regions.insert(region1);
    if (var1 == 0 || var2 == 0){
        if (DprimeR2inMem) DprimeRsquareSink[region1][region2] = result;
        return result;
    }
    r2 = cov12 * cov12 / (var1 * var2);
    result = {abs(Dprime),r2};
    if (DprimeR2inMem) DprimeRsquareSink[region1][region2] = result;
    if (al1 != "" && alleles1 != al1) unmatched_alleles.insert(region1 + " " + al1 +" "+ alleles1);
    if (al2 != "" && alleles2 != al2) unmatched_alleles.insert(region2 + " " + al2 +" "+ alleles2);
    return result;

}


void rtc_data::printPTTFR(){
    map < int ,vector < pairsToTestForRTC > >::iterator it;
    for (it = pheno_eqtls.begin() ; it != pheno_eqtls.end(); it++)
        for(int i = 0 ; i < (it->second).size() ; i++)
            cout << phenotype_id[it->first] << " " << (it->second)[i];
}

void rtc_data::mapVariantsToColdspots(){
    vrb.title("Mapping variants to coldspots");
    for (int g = 0 ; g < genotype_count; g++ ){
        if (g+1 % 100000 == 0 ) vrb.bullet(stb.str(g+1) + " genotypes mapped");
        if (coldspot_bins_p.find(genotype_chr[g]) != coldspot_bins_p.end()){
            int max = (coldspot_bins_p[genotype_chr[g]].rbegin()->second).back()->end;
            if (genotype_start[g] > max){
                coldspot_end_idx[genotype_chr[g]].push_back(g);
                continue;
            }
            int bin = genotype_start[g] / bin_size;
            if (coldspot_bins_p[genotype_chr[g]].find(bin) != coldspot_bins_p[genotype_chr[g]].end()){
                for (int c = 0 ; c < coldspot_bins_p[genotype_chr[g]][bin].size(); c++){
                    if (genotype_start[g] >= coldspot_bins_p[genotype_chr[g]][bin][c]->start && genotype_start[g] <= coldspot_bins_p[genotype_chr[g]][bin][c]->end ){
                        coldspot_bins_p[genotype_chr[g]][bin][c]->coldspot_variant_idx.push_back(g);
                        break;
                    }
                }
            }
        }
    }
    vrb.bullet(stb.str(genotype_count) + " genotypes mapped to coldspots");
}

void rtc_data::calculateRTC(string fout){
    vrb.title("Calculating RTC");
    output_file fdo (fout);
    if (fdo.fail()) vrb.error("Cannot open file [" + fout + "]");
    map < int ,vector < pairsToTestForRTC > >::iterator it;
    unsigned int count_pheno = 1;
    unsigned int done = 0;
    for (it = pheno_eqtls.begin(); it != pheno_eqtls.end(); it++){
        vrb.title("Processing phenotype [" + phenotype_id[it->first] + "] " + stb.str(count_pheno) + " / " + stb.str(pheno_eqtls.size()) );
        vrb.bullet(stb.str(it->second.size()) + " pairwise tests");
        count_pheno++;
        for (int p = 0 ; p < it->second.size(); p++){
            fdo << genotype_id[it->second[p].test_snp_idx] << " " ;
            fdo << genotype_id[it->second[p].eqtl_snp_idx] << " " ;
            fdo << phenotype_id[it->first] << " ";
            fdo << genotype_chr[it->second[p].test_snp_idx] << " " ;
            fdo << genotype_start[it->second[p].test_snp_idx] << " " ;
            fdo << it->second[p].test_snp_rank << " ";
            fdo << genotype_chr[it->second[p].eqtl_snp_idx] << " " ;
            fdo << genotype_start[it->second[p].eqtl_snp_idx] << " " ;
            fdo << it->second[p].eqtl_snp_rank << " ";
            fdo << phenotype_chr[it->first] << " ";
            fdo << phenotype_start[it->first] << " ";
            fdo << abs(genotype_start[it->second[p].test_snp_idx] - genotype_start[it->second[p].eqtl_snp_idx]) << " ";
            if (genotype_chr[it->second[p].test_snp_idx] == phenotype_chr[it->first]) fdo << abs(genotype_start[it->second[p].test_snp_idx] - phenotype_start[it->first]) << " ";
            else fdo << "NA ";
            fdo << it->second[p].test_snp_coldspot_idx << " " ;
            fdo << it->second[p].eqtl_snp_coldspot_idx << " " ;
            vector < int > genotype_idx_to_test;
            int pI = options.count("debug") ? it->first : -1 ;
            if (it->second[p].eqtl_snp_coldspot_idx >= 0 && it->second[p].test_snp_coldspot_idx >= 0 && genotype_chr[it->second[p].test_snp_idx] == genotype_chr[it->second[p].eqtl_snp_idx] && (it->second[p].eqtl_snp_coldspot_idx == it->second[p].test_snp_coldspot_idx || it->second[p].Dprime >= Dprime_cutoff)){
                if (it->second[p].test_snp_idx == it->second[p].eqtl_snp_idx){
                    fdo << all_coldspots_p[it->second[p].test_snp_coldspot_idx]->start << " ";
                    fdo << all_coldspots_p[it->second[p].test_snp_coldspot_idx]->end << " ";
                    if (sample_iterations > 0) fdo << "1 1 1 0 1" << endl;
                    else fdo << "1 1 1" << endl;
                    done++;
                    continue;
                }
                int si,ei;
                if (it->second[p].eqtl_snp_coldspot_idx < it->second[p].test_snp_coldspot_idx){
                    si = it->second[p].eqtl_snp_coldspot_idx;
                    ei = it->second[p].test_snp_coldspot_idx;
                }else{
                    ei = it->second[p].eqtl_snp_coldspot_idx;
                    si = it->second[p].test_snp_coldspot_idx;
                }
                fdo << all_coldspots_p[si]->start << " ";
                fdo << all_coldspots_p[ei]->end << " ";
                for (int csi = si ; csi <= ei; csi++){
                    genotype_idx_to_test.insert(genotype_idx_to_test.end(),all_coldspots_p[csi]->coldspot_variant_idx.begin(),all_coldspots_p[csi]->coldspot_variant_idx.end());
                }
            } else if (it->second[p].eqtl_snp_coldspot_idx == -1 && it->second[p].test_snp_coldspot_idx == -1 && genotype_chr[it->second[p].test_snp_idx] == genotype_chr[it->second[p].eqtl_snp_idx]){
                genotype_idx_to_test=coldspot_end_idx[genotype_chr[it->second[p].test_snp_idx]];
                fdo << genotype_start[genotype_idx_to_test[0]] << " ";
                fdo << genotype_start[genotype_idx_to_test.back()] << " ";
            } else if (it->second[p].eqtl_snp_coldspot_idx >= -1 && it->second[p].test_snp_coldspot_idx == -1 && genotype_chr[it->second[p].test_snp_idx] == genotype_chr[it->second[p].eqtl_snp_idx] && it->second[p].Dprime >= Dprime_cutoff){
                string chr = genotype_chr[it->second[p].eqtl_snp_idx];
                int start_idx = it->second[p].eqtl_snp_coldspot_idx;
                for (int csi = start_idx ; csi < all_coldspots_p.size() && all_coldspots_p[csi]->chr == chr; csi++){
                    genotype_idx_to_test.insert(genotype_idx_to_test.end(),all_coldspots_p[csi]->coldspot_variant_idx.begin(),all_coldspots_p[csi]->coldspot_variant_idx.end());
                }
                genotype_idx_to_test.insert(genotype_idx_to_test.end(),coldspot_end_idx[chr].begin(),coldspot_end_idx[chr].end());
                fdo << all_coldspots_p[start_idx]->start << " ";
                fdo << genotype_start[genotype_idx_to_test.back()] << " ";
            } else if (it->second[p].eqtl_snp_coldspot_idx == -1 && it->second[p].test_snp_coldspot_idx >= -1 && genotype_chr[it->second[p].test_snp_idx] == genotype_chr[it->second[p].eqtl_snp_idx] && it->second[p].Dprime >= Dprime_cutoff){
                string chr = genotype_chr[it->second[p].test_snp_idx];
                int start_idx = it->second[p].test_snp_coldspot_idx;
                for (int csi = start_idx ; csi < all_coldspots_p.size() && all_coldspots_p[csi]->chr == chr; csi++){
                    genotype_idx_to_test.insert(genotype_idx_to_test.end(),all_coldspots_p[csi]->coldspot_variant_idx.begin(),all_coldspots_p[csi]->coldspot_variant_idx.end());
                }
                genotype_idx_to_test.insert(genotype_idx_to_test.end(),coldspot_end_idx[chr].begin(),coldspot_end_idx[chr].end());
                fdo << all_coldspots_p[start_idx]->start << " ";
                fdo << genotype_start[genotype_idx_to_test.back()] << " ";
            }else{
                if (sample_iterations > 0){
                    if (it->second[p].Dprime == -9) fdo << "NA NA NA NA NA NA NA"<<endl;
                    else fdo << "NA NA NA " << it->second[p].Dprime << " " << it->second[p].R2 << " NA NA" << endl;
                }else{
                    if (it->second[p].Dprime == -9) fdo << "NA NA NA NA NA"<<endl;
                    else fdo << "NA NA NA " << it->second[p].Dprime << " " << it->second[p].R2 << endl;
                }
                continue;
            }
            done++;
            vector < double > corrs(genotype_idx_to_test.size());
            double test_snp_corr = 0.0 ;
            vector < float > genotype_eqtl;
            if (genotypeSink.count(it->second[p].eqtl_snp_idx)){
            	genotype_eqtl = genotypeSink[it->second[p].eqtl_snp_idx];
            }else{
				genotype_eqtl = genotype_val[it->second[p].eqtl_snp_idx];
				normalize(genotype_eqtl);
				if(DprimeR2inMem) genotypeSink[it->second[p].eqtl_snp_idx] = genotype_eqtl;
            }
            vector < float > phenotype_eqtl;
            vector < float > raw_phenotype_eqtl;

            phenotype_eqtl = phenotype_val[it->first];
            if (it->second[p].other_conditional_signal_idx.size()){
                residualizer covariate_engine (sample_count);
                for (int o = 0 ; o < it->second[p].other_conditional_signal_idx.size(); o++){
                    covariate_engine.push(genotype_val[it->second[p].other_conditional_signal_idx[o]]);
                }
                covariate_engine.build();
                covariate_engine.residualize(phenotype_eqtl);
            }
            if (sample_iterations > 0 ) {
                raw_phenotype_eqtl = phenotype_eqtl;

            }
            if (options.count("normal")) normalTransform(phenotype_eqtl);
            normalize(phenotype_eqtl);


            for (int s = 0 ; s < genotype_idx_to_test.size() ; s++){
                vector < float > test = genotype_val[genotype_idx_to_test[s]];
                normalize(test);
                vector <float> new_pheno = correct(test,phenotype_eqtl);
                if (options.count("normal")) normalTransform(new_pheno);
                normalize(new_pheno);
                corrs[s] = abs(getCorrelation(genotype_eqtl, new_pheno));
                if( genotype_idx_to_test[s] == it->second[p].test_snp_idx) test_snp_corr = corrs[s];
            }
            sort(corrs.begin(),corrs.end());
            int rank = -1;
            for (int i = 0 ; i<corrs.size() && corrs[i] <= test_snp_corr; i++) if(corrs[i] == test_snp_corr) rank = i;
            double RTC = ((double) corrs.size() - (double) rank) / (double) corrs.size();
            string dprime = it->second[p].Dprime != -9 ? stb.str(it->second[p].Dprime) : "NA";
            string rsquared = it->second[p].R2 != -9 ? stb.str(it->second[p].R2) : "NA";
            if (sample_iterations >0 ) {
                vector <double > res = sampleRTC(genotype_idx_to_test, raw_phenotype_eqtl, it->second[p].eqtl_snp_idx, RTC,pI);
            	fdo << RTC <<" " << dprime << " " << rsquared << " " << res.at(0) << " " << res.at(1)  << endl;
            }else fdo << RTC <<" " << dprime << " " << rsquared << endl;
        }
    }
    vrb.print("\n\n  * " + stb.str(done) + " actual RTCs calculated.");
    fdo.close();
}
