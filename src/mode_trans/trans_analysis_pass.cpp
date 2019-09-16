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

void trans_data::runTransPass(string fvcf, string fout) {
	string fhits = fout + ".hits.txt.gz";
	string fbins = fout + ".bins.txt.gz";
	string fbest = fout + ".best.txt.gz";

	vrb.title("Sweep through genotype data in [" + fvcf + "]");
	bcf_sweep_t * sw = bcf_sweep_init(fvcf.c_str());
	if (!sw) vrb.error("Cannot open file for reading [" + fvcf + "]");
	bcf_hdr_t * hdr  = bcf_sweep_hdr(sw);
	if (!hdr) vrb.error("Cannot read header!");

	unsigned long n_pos_tests = 0, n_variants = 0, n_curr_tests = 0, n_prev_tests = 0;
	vector < double > best_hits = vector < double > (phenotype_count, 0);
	vector < string > best_var = vector < string > (phenotype_count, "");
	vector < unsigned long > bins_hits = vector < unsigned long > (n_bins, 0);

    int mDS = 0, mGT = 0;
    float * vDS = NULL;
    int * vGT = NULL;
    bcf1_t * rec;
    double step_bins = n_bins * 1.0 / correlation_threshold;

	timer testing_speed_timer;
	testing_speed_timer.clock();

	int mode_GT_DS = -1;	//-1 means unset / 0 means dosages / 1 means genotypes
    output_file fdhits (fhits);
    if (fdhits.fail()) vrb.error("Cannot open file for writing [" + fhits + "]");
    while ( (rec = bcf_sweep_fwd(sw)) ) {
		bcf_unpack(rec, BCF_UN_STR);
		string vid = string(rec->d.id);
		string chr = string(bcf_hdr_id2name(hdr, rec->rid));
		int pos = rec->pos + 1;
		if (rec->n_allele == 2) {
			bool variant_is_to_be_processed = false;
			if (mode_GT_DS != 0 && (bcf_get_genotypes(hdr, rec, &vGT, &mGT) == (2 * sample_count))) {
				if (vDS == NULL) vDS = (float*)malloc(sample_count*sizeof(float));
				computeDosages(vGT, vDS);
				variant_is_to_be_processed = true;
				mode_GT_DS = 1;
			} else if (mode_GT_DS != 1 && (bcf_get_format_float(hdr, rec, "DS", &vDS, &mDS) == sample_count)) {
				variant_is_to_be_processed = true;
				mode_GT_DS = 0;
			}

			if (variant_is_to_be_processed) {
				imputeGenotypes(vDS);
				normalize(vDS);
				for (int p = 0 ; p < phenotype_count ; p ++) {
					if (chr != phenotype_chr[p] || abs(phenotype_start[p] - pos) > cis_window) {
						double rcorr = fastCorrelation(phenotype_val[p], vDS);
						double acorr = abs (rcorr);
						if (acorr >= correlation_threshold) {
							double npval = getNominalPvalue(acorr, sample_count - 2);
							double apval = getAdjustedPvalue(npval);
							fdhits << phenotype_id[p] << " " << phenotype_chr[p] << " " << phenotype_start[p] << " " << vid << " " << chr << " " << pos << " " << npval << " " << apval << " " << rcorr << endl;
							n_pos_tests ++;
						} else {
							unsigned int idx_bin = (unsigned int) floor(acorr * step_bins);
							bins_hits[idx_bin] ++;
						}
						if (acorr > best_hits[p]) {
							best_hits[p] = acorr;
							best_var[p] = vid;
						}

						n_curr_tests ++;

						if (n_curr_tests % 1000000 == 0) {
							unsigned int elapsed_msecs = testing_speed_timer.rel_time();
							vrb.bullet("#variants=" + stb.str(n_variants) + "\t#hits=" + stb.str(n_pos_tests) + "/" + stb.str(n_curr_tests) + "\tspeed=" + stb.str((n_curr_tests - n_prev_tests) * 1.0 / (elapsed_msecs * 1000), 2) + "MT/s");
							testing_speed_timer.clock();
							n_prev_tests = n_curr_tests;
						}
					}
				}
			}
		}
		n_variants ++;
	}
	unsigned int elapsed_msecs = testing_speed_timer.rel_time();
	vrb.bullet("#variants=" + stb.str(n_variants) + "\t#hits=" + stb.str(n_pos_tests) + "/" + stb.str(n_curr_tests) + "\tspeed=" + stb.str((n_curr_tests - n_prev_tests) * 1.0 / (elapsed_msecs * 1000), 2) + "MT/s");
	n_prev_tests = n_curr_tests;
	bcf_sweep_destroy(sw);
	fdhits.close();

	output_file fdbins(fbins);
    if (fdbins.fail()) vrb.error("Cannot open file for writing [" + fbins + "]");
	for (int b = 0 ; b < n_bins ; b ++) {
		double start_corr = (b * 1.0 / n_bins) * correlation_threshold;
		double end_corr = ((b+1) * 1.0 / n_bins) * correlation_threshold;
		fdbins << b << " " << start_corr << " " << end_corr << " " << getNominalPvalue(start_corr, sample_count - 2) << " " << getNominalPvalue(end_corr, sample_count - 2) << " " << bins_hits[b] << endl;
	}
	fdbins.close();

	output_file fdbest(fbest);
	if (fdbest.fail()) vrb.error("Cannot open file for writing [" + fbest + "]");
	for (int p = 0 ; p < phenotype_count ; p ++) {
		double npval = getNominalPvalue(best_hits[p], sample_count - 2);
		double apval = getAdjustedPvalue(npval);
		fdbest << phenotype_id[p] << " " << apval << " " << npval << " " << best_var[p] << endl;
	}
	fdbest.close();
}
