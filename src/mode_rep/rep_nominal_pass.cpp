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

#include "rep_data.h"

void rep_data::runNominalPass(string fout) {

	//STEP0: INITIALIZATION OF IO
	output_file fdo (fout);
	if (fdo.fail()) vrb.error("Cannot open file [" + fout + "]");

	//STEP2: INITIALIZE A WORKING COPY OF GENOTYPES
	vector < double > genotype_sd = vector < double > (genotype_count, 0.0);
	for (unsigned int v = 0 ; v < genotype_count ; v ++) {
		genotype_sd[v] = basic_stats(genotype_val[v]).sd();
		normalize(genotype_val[v]);
	}

	//STEP3: INITIALIZE A WORKING COPY OF PHENOTYPES
	vector < double > phenotype_sd = vector < double > (phenotype_count, 0.0);
	for (unsigned int p = 0 ; p < phenotype_count ; p ++) {
		phenotype_sd[p] = basic_stats(phenotype_val[p]).sd();
		normalize(phenotype_val[p]);
	}

	//STEP4: MAIN SWEEP THROUGH PHENOTYPES
	for (unsigned int p = 0 ; p < qtl_count ; p ++) {

		int i_phenotype = qtl_idx.first[p];
		int i_genotype = qtl_idx.second[p];

		if (i_phenotype >= 0 && i_genotype >= 0) {

			//STEP4: VERBOSE PROCESSED PHENOTYPES
			vrb.title("Testing [" + phenotype_id[i_phenotype] + "] against [" + genotype_id[i_genotype] + "] / Progress [" + stb.str(p+1) + "/" + stb.str(qtl_count) + "]");
			double curr_correlation = getCorrelation(genotype_val[i_genotype], phenotype_val[i_phenotype]);
			double pvalue = getPvalue(curr_correlation, sample_count - 2);
			double slope = getSlope(curr_correlation, genotype_sd[i_genotype], phenotype_sd[i_phenotype]);

			//STEP9: VERBOSE BEST HIT
			vrb.bullet("p=" + stb.str(pvalue) + ", s=" + stb.str(slope, 4));

			//STEP10: PRINT RESULTS IN FILE
			fdo << phenotype_id[i_phenotype];
			fdo << " " << phenotype_chr[i_phenotype];
			fdo << " " << phenotype_start[i_phenotype];
			fdo << " " << phenotype_end[i_phenotype];
			fdo << " " << (phenotype_neg[i_phenotype]?"-":"+");
			fdo << " " << genotype_id[i_genotype];
			fdo << " " << genotype_chr[i_genotype];
			fdo << " " << genotype_start[i_genotype];
			fdo << " " << genotype_end[i_genotype];
			fdo << " " << pvalue;
			fdo << " " << slope;
			fdo << endl;
		}
	}
	fdo.close();
}
