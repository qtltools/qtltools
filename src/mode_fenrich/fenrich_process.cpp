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

#include "fenrich_data.h"

void fenrich_data::runEnrichmentPass(string fout) {

	//1. Nominal pass
	vrb.title("Enrichment analysis");
	unsigned int obs_overlaps = countOverlaps();
	vrb.bullet("#observed overlaps = " + stb.str(obs_overlaps) + " / " + stb.str(qtl_count) + " (" + stb.str(obs_overlaps * 100.0 / qtl_count, 2) + "%)");

	//2. Permutation pass
	basic_stats null_stat;
	vector < int > rorder = vector < int > (tss_count, -1);
	for ( int t = 0 ; t < tss_count ; t ++) rorder[t] = t;
	vector < unsigned int > null_overlaps;
	for (int p = 0 ; p < n_permutation ; p++) {
		random_shuffle(rorder.begin(), rorder.end());
		qtl_order = vector < int > (rorder.begin() , rorder.begin() + qtl_pos.size());
		unsigned int no = countOverlaps();
		null_overlaps.push_back(no);
		null_stat.push(no * 1.0);
	}
	sort(null_overlaps.begin(), null_overlaps.end());
	vrb.bullet("#null overlaps = " + stb.str(null_stat.mean(), 2) + " +/- " + stb.str(null_stat.sd(), 2) + " (" + stb.str(null_stat.mean() * 100.0 / qtl_count, 2) + "% +/- " + stb.str(null_stat.sd() * 100.0 / qtl_count, 2) + "%)");

	//3. calculate empirical p-value
	unsigned int n_smaller = 0, n_bigger = 0;
	for (int p = 0 ; p < n_permutation ; p++) {
		if (obs_overlaps >= null_overlaps[p]) n_smaller ++;
		if (obs_overlaps <= null_overlaps[p]) n_bigger ++;
	}
	double epval = min((min(n_smaller, n_bigger) * 2.0 + 1) / (n_permutation + 1), 1.0);
	vrb.bullet("empirical p-value = " + stb.str(epval));

	//4. Compute odd ratios
	double obs_freq = obs_overlaps * 1.0 / qtl_count;
	double exp_freq_med = null_overlaps[(int)round(null_overlaps.size() * 0.500)] * 1.0 / qtl_count;
	double exp_freq_upv = null_overlaps[(int)round(null_overlaps.size() * 0.975)] * 1.0 / qtl_count;
	double exp_freq_dnv = null_overlaps[(int)round(null_overlaps.size() * 0.025)] * 1.0 / qtl_count;
	double odd_ratio_med = (obs_freq * (1 - exp_freq_med)) / (exp_freq_med * (1 - obs_freq));
	double odd_ratio_upv = (obs_freq * (1 - exp_freq_upv)) / (exp_freq_upv * (1 - obs_freq));
	double odd_ratio_dnv = (obs_freq * (1 - exp_freq_dnv)) / (exp_freq_dnv * (1 - obs_freq));
	vrb.bullet("Odd ratio = " + stb.str(odd_ratio_med, 4) + " [" + stb.str(odd_ratio_dnv, 4) + "," + stb.str(odd_ratio_upv, 4) + "]");

	//5. Write output
	output_file fdo (fout.c_str());
	if (fdo.fail()) vrb.error("Cannot open output file!");
	fdo << obs_overlaps << " " << qtl_count << " " << null_stat.mean() << " " << null_stat.sd() << " " << epval << " " << odd_ratio_dnv << " " << odd_ratio_med << " " << odd_ratio_upv << endl;
	fdo.close();
}
