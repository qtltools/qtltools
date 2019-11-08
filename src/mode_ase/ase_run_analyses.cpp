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


void ase_data::calculateRefToAltBias(string olog){
	output_file fdo;
	vrb.title("Calculating reference allele mapping bias");
	if (olog != ""){
		vrb.bullet("Writing failed variants to [" + olog + "]");
		fdo.append(olog);
		if (fdo.fail()) vrb.error("Cannot open file [" + olog +"]");
	}
	vector < unsigned int > all_total_counts;
	map < string , vector <unsigned int> > alleles_total_counts;
	vector < ase_site * > filtered_variants;
	unsigned int filtered_cov = 0 , filtered_bas = 0;
	for (int i = 0 ; i < passing_variants.size(); i++){
		if (passing_variants[i].total_count < param_min_cov_for_ref_alt) { filtered_cov++; if(olog!="") fdo<<"BIAS_COV " << passing_variants[i].sid << endl; continue;}
		if (param_both_alleles_seen_bias && (passing_variants[i].alt_count == 0 || passing_variants[i].ref_count == 0 )) {filtered_bas++; if(olog!="") fdo<<"BIAS_BOTH_ALLELES " << passing_variants[i].sid << endl; continue;}
		filtered_variants.push_back(&passing_variants[i]);
		if (alleles_total_counts.count(passing_variants[i].alleles)) alleles_total_counts[passing_variants[i].alleles].push_back(passing_variants[i].total_count);
		else alleles_total_counts[passing_variants[i].alleles] = vector <unsigned int> (1, passing_variants[i].total_count);
		all_total_counts.push_back(passing_variants[i].total_count);
	}

	if (filtered_cov) vrb.bullet(stb.str(filtered_cov) + " sites with coverage less than " + stb.str(param_min_cov_for_ref_alt));
	if (filtered_bas) vrb.bullet(stb.str(filtered_bas) + " sites where only one allele was observed");

	bool calculate_all = false;
	vector <string> add;
	for (auto it = all_allele_combinations.begin(); it != all_allele_combinations.end(); it++){
		if (alleles_total_counts.count(*it) == 0 || alleles_total_counts[*it].size() < param_min_sites_for_ref_alt) {
			calculate_all = true;
			add.push_back(*it);
			vrb.bullet("Bias for " + *it + " will be calculated from all sites since " + stb.str(alleles_total_counts.count(*it) ? alleles_total_counts[*it].size() : 0) + " is too low");
		}
		else{
			sort(alleles_total_counts[*it].begin(), alleles_total_counts[*it].end());
			unsigned int max_index = alleles_total_counts[*it].size() * param_sample - 1;
			unsigned int max = alleles_total_counts[*it][max_index];
			unsigned int refc = 0, altc = 0, sc = 0;
			for (int i = 0; i < filtered_variants.size(); i++){
				if (filtered_variants[i]->alleles != *it) continue;
				if (filtered_variants[i]->total_count <= max){
					refc += filtered_variants[i]->ref_count;
					altc += filtered_variants[i]->alt_count;
				}else{
					sc++;
					double ratio = (double) filtered_variants[i]->ref_count / (double) filtered_variants[i]->total_count;
					for (int r = 0; r < max; r++){
						if (rng.getDouble() <= ratio) refc++;
						else altc++;
					}
				}
			}
			ref_to_alt_bias[*it] = (double) refc / ((double) altc + (double) refc);
			vrb.bullet("Bias for " + *it + " from " + stb.str(alleles_total_counts[*it].size()) + " sites is " + stb.str(ref_to_alt_bias[*it]));
			if (sc) vrb.bullet("Subsampled " + stb.str(sc) + " sites to " + stb.str(max));
		}
	}
	if (calculate_all){
		sort(all_total_counts.begin(),all_total_counts.end());
		unsigned int max_index = all_total_counts.size() * param_sample - 1;
		unsigned int max = all_total_counts[max_index];
		unsigned int refc = 0, altc = 0, sc = 0;
		for (int i = 0; i < filtered_variants.size(); i++){
			if (filtered_variants[i]->total_count <= max){
				refc += filtered_variants[i]->ref_count;
				altc += filtered_variants[i]->alt_count;
			}else{
				sc++;
				double ratio = (double) filtered_variants[i]->ref_count / (double) filtered_variants[i]->total_count;
				for (int r = 0; r < max; r++){
					if (rng.getDouble() <= ratio) refc++;
					else altc++;
				}
			}
		}
		double all_bias = (double) refc / ((double) altc + (double) refc);
		vrb.bullet("Bias from all sites is " + stb.str(all_bias));
		for (int e = 0; e < add.size(); e++) ref_to_alt_bias[add[e]] = all_bias;
		if (sc) vrb.bullet("Subsampled " + stb.str(sc) + " sites to " + stb.str(max));
	}
}

void ase_data::calculateASE(string fout , string olog){
	vrb.title("Calculating ASE and writing to [" + fout + "]");
	output_file fdoo(fout);
	if (fdoo.fail()) vrb.error("Cannot open file [" + fout +"]");
	output_file fdo;
	if (olog != ""){
		vrb.bullet("Writing failed variants to [" + olog + "]");
		fdo.append(olog);
		if (fdo.fail()) vrb.error("Cannot open file [" + olog +"]");
	}
	fdoo <<"INDIVIDUAL\tRSID\tCHR\tPOS\tALLELES\tBOTH_ALLELES_SEEN\tMIN_ALLELE_RATIO\tREF_COUNT\tNONREF_COUNT\tTOTAL_COUNT\tWEIGHTED_REF_COUNT\tWEIGHTED_NONREF_COUNT\tWRC_MINUS_WNC\tALLELES_SEEN\tREF_ALLELE\tALT_ALLELE\tOTHER_COUNT\tREF_RATIO\tPVALUE\tCONCERN\tEXON_INFO";
	if (print_stats) fdoo <<"\tUNMAPPED\tSECONDARY\tFAIL_MAPQ\tSKIPPED\tFAIL_BASEQ\tFAILQC\tDUPLICATE\tINDEL\tMATE_UNMAPPED\tWRONG_ORIENTATION\tNOT_PROPER_PAIR";
	fdoo <<endl;
	for (int i = 0; i < passing_variants.size(); i++){
		if (passing_variants[i].total_count >= param_min_cov){
			if (param_both_alleles_seen && (passing_variants[i].alt_count == 0 || passing_variants[i].ref_count == 0 )){if(olog!="") fdo<<"ASE_BOTH_ALLELES " << passing_variants[i].sid << endl; continue;};
			if (ref_to_alt_bias[passing_variants[i].alleles] >= 0 && ref_to_alt_bias[passing_variants[i].alleles] <= 1) passing_variants[i].calculatePval(ref_to_alt_bias[passing_variants[i].alleles]);
			else vrb.error("Reference allele mapping bias is incorrect [" + stb.str(ref_to_alt_bias[passing_variants[i].alleles]) + "]");
			assignGenesToAseSite(passing_variants[i]);
			fdoo << sample_id[0] << "\t" << passing_variants[i];
			if (print_stats) fdoo << "\t" << passing_variants[i].stats;
			fdoo << endl;
		}else if(olog!="") fdo<<"ASE_COV " << passing_variants[i].sid << endl;
	}

}
