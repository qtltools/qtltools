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

#include "gwas_data.h"

void gwas_data::runGwasPass(string fvcf, string fout) {
	vrb.title("Sweep through genotype data in [" + fvcf + "]");
	bcf_sweep_t * sw = bcf_sweep_init(fvcf.c_str());
	if (!sw) vrb.error("Cannot open file for reading [" + fvcf + "]");
	bcf_hdr_t * hdr  = bcf_sweep_hdr(sw);
	if (!hdr) vrb.error("Cannot read header!");

	//Specify samples to be read
	string concat_id=sample_id[0];
	for (int i = 1; i < sample_id.size() ; i ++) concat_id+=","+sample_id[i];
	bcf_hdr_set_samples(hdr, concat_id.c_str(), 0);
	vrb.bullet("#samples=" + stb.str(sample_count));

	unsigned long n_pos_tests = 0, n_variants = 0, n_curr_tests = 0, n_prev_tests = 0;
    int mDS = 0, mGT = 0;
    float * vDS = NULL;
    int * vGT = NULL;
    bcf1_t * rec;

	timer testing_speed_timer;
	testing_speed_timer.clock();

	output_file fdhits (fout);
    if (fdhits.fail()) vrb.error("Cannot open file for writing [" + fout + "]");
    while ( (rec = bcf_sweep_fwd(sw)) ) {
		bcf_unpack(rec, BCF_UN_STR);
		string vid = string(rec->d.id);
		string chr = string(bcf_hdr_id2name(hdr, rec->rid));
		int pos = rec->pos + 1;
		if (rec->n_allele == 2) {
			bool variant_is_to_be_processed = false;
			int n_captured = bcf_get_format_float(hdr, rec, "DS", &vDS, &mDS);
			if (n_captured == sample_count) {
				imputeGenotypes(vDS);
				normalize(vDS);
				for (int p = 0 ; p < phenotype_count ; p ++) {
					double rcorr = fastCorrelation(phenotype_val[p], vDS);
					double npval = getNominalPvalue(abs(rcorr), sample_count - 2);
					fdhits << phenotype_id[p] << " " << chr << " " << pos << " " << vid << " " << npval << " " << rcorr << endl;
					if (n_curr_tests % 1000 == 0) {
						unsigned int elapsed_msecs = testing_speed_timer.rel_time();
						vrb.bullet("#variants=" + stb.str(n_variants) + "\t#tests=" + stb.str(n_curr_tests) + "\t" + stb.str(n_curr_tests) + "\tspeed=" + stb.str((n_curr_tests - n_prev_tests) * 1.0 / (elapsed_msecs * 1000), 2) + "MT/s");
						testing_speed_timer.clock();
						n_prev_tests = n_curr_tests;
					}
				}
				n_curr_tests ++;
			}
		}
		n_variants ++;
	}
	unsigned int elapsed_msecs = testing_speed_timer.rel_time();
	vrb.bullet("#variants=" + stb.str(n_variants) + "\t#tests=" + stb.str(n_curr_tests) + "\t" + stb.str(n_curr_tests) + "\tspeed=" + stb.str((n_curr_tests - n_prev_tests) * 1.0 / (elapsed_msecs * 1000), 2) + "MT/s");
	n_prev_tests = n_curr_tests;
	bcf_sweep_destroy(sw);
	fdhits.close();
}
