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

void genrich_data::binningAllVariants() {
	vrb.title("Bin all variants");
	vrb.bullet("Distance +/- " + stb.str(bin_distance) + "bp and MAF +/- " + stb.str(bin_maf, 4));
	for (int v = 0 ; v < genotype_chr.size() ; v ++) {
		if (genotype_qtl[v]) {
			int idx_bin = -1;
			unsigned int n_bin = bin_min_maf.size();
			for (int b = 0 ; b < n_bin && idx_bin < 0 ; b++) {
				bool in_freq = (genotype_maf[v] >= bin_min_maf[b] && genotype_maf[v] < bin_max_maf[b]);
				bool in_dist = (genotype_dist[v] >= bin_min_dist[b] && genotype_dist[v] < bin_max_dist[b]);
				if (in_freq && in_dist) idx_bin = b;
			}
			if (idx_bin < 0) {
				idx_bin = bin_min_maf.size();
				bin_min_maf.push_back(genotype_maf[v] - bin_maf);
				if (bin_min_maf.back() < 0) bin_min_maf.back() = 0;
				bin_max_maf.push_back(genotype_maf[v] + bin_maf);
				if (bin_max_maf.back() >= 0.5) bin_max_maf.back() = 0.5;
				bin_min_dist.push_back(genotype_dist[v] - bin_distance);
				bin_max_dist.push_back(genotype_dist[v] + bin_distance);
			}
			genotype_bin[v] = idx_bin;
		}
	}
	vrb.bullet("Number of bins made from QTL data = " + stb.str(bin_min_maf.size()));

	unsigned int n_binned = 0, n_nbinned = 0;
	for (int v = 0 ; v < genotype_chr.size() ; v ++) {
		if (!genotype_qtl[v]) {
			int idx_bin = -1;
			unsigned int n_bin = bin_min_maf.size();
			for (int b = 0 ; b < n_bin && idx_bin < 0 ; b++) {
				bool in_freq = (genotype_maf[v] >= bin_min_maf[b] && genotype_maf[v] < bin_max_maf[b]);
				bool in_dist = (genotype_dist[v] >= bin_min_dist[b] && genotype_dist[v] < bin_max_dist[b]);
				if (in_freq && in_dist) idx_bin = b;
			}
			genotype_bin[v] = idx_bin;
			if (idx_bin >= 0) n_binned ++;
			else n_nbinned ++;
		}
	}
	vrb.bullet("Number of reference variants falling within bins = " + stb.str(n_binned));
	vrb.bullet("Number of reference variants outside of any bins = " + stb.str(n_nbinned));
}
