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

#include "match_data.h"

void match_data::writeOutput(string filename) {

	vrb.title("Write summary report in ["  + filename + "]");
	output_file fd (filename);

	fd << "SampleID n_geno_missing n_het_total n_hom_total n_het_covered n_hom_covered n_het_consistent n_hom_consistent perc_het_consistent perc_hom_consistent n_het_in_ase" << endl;

	for (int i = 0 ; i < sample_count ; i ++) {
		unsigned int n_mis_tot = 0, n_het_tot = 0, n_hom_tot = 0, n_het_cov = 0, n_hom_cov = 0, n_het_fit = 0, n_hom_fit = 0, n_het_ase = 0;
		for (int r = 0; r < regions.size() ; r ++) {
			for (int s = 0 ; s < sites[r].size() ; s ++) {
				unsigned int coverage = sites[r][s].cref + sites[r][s].calt;

				//A. Missing genotype
				if (!gen_ref[r][s][i] && !gen_alt[r][s][i]) n_mis_tot ++;

				//B. Homozygous
				if (gen_ref[r][s][i] != gen_alt[r][s][i]) {
					n_hom_tot ++;
					if (coverage >= param_min_cov) {
						n_hom_cov ++;
						if (sites[r][s].cref == 0 || sites[r][s].calt == 0) n_hom_fit ++;
					}
				}

				//C. Heterozygous
				if (gen_ref[r][s][i] && gen_alt[r][s][i]) {
					n_het_tot ++;
					if (coverage >= param_min_cov) {
						n_het_cov ++;
						if (sites[r][s].cref > 0 && sites[r][s].calt > 0) {
							n_het_fit ++;
							double bvalue = match_binomialTest(sites[r][s].cref, sites[r][s].cref + sites[r][s].calt, 0.5);
							if (bvalue < param_min_pval) n_het_ase ++;
						}
					}
				}
			}
		}

		fd << sample_id[i];
		fd << " " << n_mis_tot;
		fd << " " << n_het_tot;
		fd << " " << n_hom_tot;
		fd << " " << n_het_cov;
		fd << " " << n_hom_cov;
		fd << " " << n_het_fit;
		fd << " " << n_hom_fit;
		fd << " " << n_het_fit*1.0/n_het_cov;
		fd << " " << n_hom_fit*1.0/n_hom_cov;
		fd << " " << n_het_ase;
		fd << endl;
	}
	fd.close();
}
