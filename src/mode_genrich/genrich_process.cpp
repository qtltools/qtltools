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

#include "genrich_data.h"

void genrich_data::overlapGWASandQTL(string fout) {

	vrb.title("Calculating observed overlap between QTLs and GWAS hits");
	vector < unsigned int > qtl_idx, gwas_idx;
	for (int v = 0 ; v < genotype_pos.size() ; v ++) {
		if (genotype_qtl[v]) qtl_idx.push_back(v);
		if (genotype_gwas[v]) gwas_idx.push_back(v);
	}
	vrb.bullet("#qtls=" + stb.str(qtl_idx.size()));
	vrb.bullet("#gwas=" + stb.str(gwas_idx.size()));

	vrb.title("Classifying null variants");
	vector < vector < unsigned int > > null_sets = vector < vector < unsigned int > > (bin_min_maf.size());
	for (int v = 0 ; v < genotype_pos.size() ; v ++) if (!genotype_qtl[v] && genotype_bin[v] >= 0) null_sets[genotype_bin[v]].push_back(v);
	basic_stats bs_null_count;
	for (int b = 0 ; b < null_sets.size() ; b ++) {
		bs_null_count.push(null_sets[b].size());
		//cerr << b << " " << null_sets[b].size() << " " << 	bin_min_maf[b] << " " << bin_max_maf[b] << " " << bin_min_dist[b] << " " << bin_max_dist[b] << endl;
	}
	vrb.bullet("#null variants per bin = " + stb.str(bs_null_count.mean(),3) + " +/-" + stb.str(bs_null_count.sd(), 3));

	vrb.title("Discarding poorly populated bins");
	int n_discarded = 0, n_removed_qtl = 0;
	for (int b = 0 ; b < null_sets.size() ; b ++) if (null_sets[b].size() < 3) {
		for (int q = 0 ; q < qtl_idx.size() ; q ++) if (genotype_bin[qtl_idx[q]] == b) {
			genotype_bin[qtl_idx[q]] = -1;
			n_removed_qtl ++;
		}
		n_discarded ++;
	}
	vrb.bullet("#bin discarded=" + stb.str(n_discarded) + " out of " + stb.str(null_sets.size()));
	vrb.bullet("#qtl discarded=" + stb.str(n_removed_qtl) + " out of " + stb.str(qtl_idx.size()));

	vrb.title("Calculating overlap between OBSERVED set of variants and GWAS hits");
	vector < bool > overlap_qtl = vector < bool > (qtl_idx.size() , false);
	for (int q = 0 ; q < qtl_idx.size() ; q ++) {
		if (genotype_bin[qtl_idx[q]] >= 0) {
			for (int g = 0 ; g < gwas_idx.size() && !overlap_qtl[q] ; g ++) overlap_qtl[q] = isSameSignal(gwas_idx[g], qtl_idx[q]);
		}
	}
	unsigned int n_obs_overlap = 0;
	for (int q = 0 ; q < qtl_idx.size() ; q ++) if (overlap_qtl[q]) n_obs_overlap ++;
	vrb.bullet("#observed overlap=" + stb.str(n_obs_overlap) + " (=" + stb.str(n_obs_overlap * 100.0 / qtl_idx.size(), 2) + " %)");

	vrb.title("Calculating overlap between NULL sets of variants and GWAS hits");
	vrb.bullet("#permutations=" + stb.str(n_permutations));
	basic_stats bs_null_overlap;
	vector < unsigned int > n_exp_overlap = vector < unsigned int > (n_permutations, 0);
	for (int p = 0 ; p < n_permutations ; p ++) {

		//step1: sample sequence of null variants
		//cerr << "1. Permutation " << p << endl;
		vector < int > seq_null_qtl;
		for (int q = 0 ; q < qtl_idx.size() ; q ++) {
			if (genotype_bin[qtl_idx[q]] >= 0) {
				unsigned int idx_bin = genotype_bin[qtl_idx[q]];
				unsigned int idx_rnd = rng.getInt(null_sets[idx_bin].size());
				//cerr << q << " " << idx_bin << " " << idx_rnd << " " <<  null_sets[idx_bin].size() << endl;
				seq_null_qtl.push_back(null_sets[idx_bin][idx_rnd]);
			} else seq_null_qtl.push_back(-1);
		}
		//cerr << "1. Size = " << seq_null_qtl.size() << endl;

		//step2: work out overlap
		//cerr << "2. Overlap " << p << endl;
		overlap_qtl = vector < bool > (qtl_idx.size() , false);
		for (int q = 0 ; q < qtl_idx.size() ; q ++) {
			if (seq_null_qtl[q] >= 0) {
				for (int g = 0 ; g < gwas_idx.size() && !overlap_qtl[q] ; g ++) overlap_qtl[q] = isSameSignal(gwas_idx[g], seq_null_qtl[q]);
			}
		}

		//step3: count overlaps
		//cerr << "3. Count overlaps " << p << endl;
		for (int q = 0 ; q < qtl_idx.size() ; q ++) if (overlap_qtl[q]) n_exp_overlap[p] ++;
		bs_null_overlap.push(n_exp_overlap[p]);
		//cerr << "3. #overlaps = " << n_exp_overlap[p] << endl;
		vrb.bullet("permutation = " + stb.str(p) + " overlaps = " + stb.str(n_exp_overlap[p]));
	}
	sort(n_exp_overlap.begin(), n_exp_overlap.end());
	vrb.bullet("#null overlap = " + stb.str(bs_null_overlap.mean(),3) + " +/-" + stb.str(bs_null_overlap.sd(), 3));

	unsigned int n_bigger = 0, n_smaller = 0;
	for (int p = 0 ; p < n_permutations ; p++) {
		if (n_obs_overlap <= n_exp_overlap[p]) n_bigger ++;
		if (n_obs_overlap >= n_exp_overlap[p]) n_smaller ++;
	}
	double epval = min((min(n_smaller, n_bigger) * 2.0 + 1.0) / (n_permutations + 1.0), 1.0);
	vrb.bullet("Empirical p-value = " + stb.str(epval));

	double obs_freq = n_obs_overlap * 1.0 / qtl_idx.size();
	double exp_freq_med = n_exp_overlap[(int)round(n_exp_overlap.size() * 0.500)] * 1.0 / qtl_idx.size();
	double exp_freq_upv = n_exp_overlap[(int)round(n_exp_overlap.size() * 0.975)] * 1.0 / qtl_idx.size();
	double exp_freq_dnv = n_exp_overlap[(int)round(n_exp_overlap.size() * 0.025)] * 1.0 / qtl_idx.size();
	double odd_ratio_med = (obs_freq * (1 - exp_freq_med)) / (exp_freq_med * (1 - obs_freq));
	double odd_ratio_upv = (obs_freq * (1 - exp_freq_upv)) / (exp_freq_upv * (1 - obs_freq));
	double odd_ratio_dnv = (obs_freq * (1 - exp_freq_dnv)) / (exp_freq_dnv * (1 - obs_freq));
	vrb.bullet("Odd ratio = " + stb.str(odd_ratio_med, 4) + " [" + stb.str(odd_ratio_dnv, 4) + "," + stb.str(odd_ratio_upv, 4) + "]");

	string fout_sum = fout + ".summary.txt";
	vrb.title ("Writing summary of enrichment analysis in [" + fout + "]");
	output_file fdo_sum (fout_sum.c_str());
	if (fdo_sum.fail()) vrb.error("Cannot open output file!");
	fdo_sum << n_obs_overlap << " " << qtl_idx.size() << " " << bs_null_overlap.mean() << " " << bs_null_overlap.sd() << " " << epval << " " << odd_ratio_dnv << " " << odd_ratio_med << " " << odd_ratio_upv << endl;
	fdo_sum.close();

	string fout_full = fout + ".full.txt.gz";
	vrb.title ("Writing full enrichment analysis outcome in [" + fout + "]");
	output_file fdo_full (fout_full.c_str());
	if (fdo_full.fail()) vrb.error("Cannot open output file!");
	fdo_full << n_obs_overlap << " " << qtl_idx.size() << endl;
	for (int p = 0 ; p < n_permutations ; p ++) fdo_full << n_exp_overlap[p] << " " << qtl_idx.size() << endl;
	fdo_full.close();
}

