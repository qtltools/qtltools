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

void gwas_data::runGwasPassOnVCF(string fvcf, string fout) {
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

	unsigned long n_variants = 0, n_curr_tests = 0, n_prev_tests = 0;
    int mDS = 0;
    float * vDS = NULL;
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
			int n_captured = bcf_get_format_float(hdr, rec, "DS", &vDS, &mDS);
			if (n_captured == sample_count) {
				imputeGenotypes(vDS);
				normalize(vDS);
				for (int p = 0 ; p < phenotype_count ; p ++) {
					double rcorr = fastCorrelation(phenotype_val[p], vDS);
					double npval = getNominalPvalue(abs(rcorr), sample_count - 2);
					fdhits << phenotype_id[p] << " " << chr << " " << pos << " " << vid << " " << npval << " " << rcorr << endl;
				}
				if ((n_curr_tests/phenotype_count) % 1000 == 0) {
					unsigned int elapsed_msecs = testing_speed_timer.rel_time();
					vrb.bullet("#variants=" + stb.str(n_variants) + "\t#tests=" + stb.str(n_curr_tests) + "\t" + stb.str(n_curr_tests) + "\tspeed=" + stb.str((n_curr_tests - n_prev_tests) * 1.0 / (elapsed_msecs * 1000), 2) + "MT/s");
					testing_speed_timer.clock();
					n_prev_tests = n_curr_tests;
				}
				n_curr_tests += phenotype_count;
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


void gwas_data::runGwasPassOnBED(string fbed, string fout) {
	vrb.title("Sweep through BED file [" + fbed + "]");
	htsFile *fp = hts_open(fbed.c_str(),"r");
	if (!fp) vrb.error("Cannot open file");
	tbx_t *tbx = tbx_index_load(fbed.c_str());
	if (!tbx) vrb.error("Cannot open index file");
	kstring_t str = {0,0,0};
	if (hts_getline(fp, KS_SEP_LINE, &str) <= 0 || !str.l || str.s[0] != tbx->conf.meta_char ) vrb.error("Cannot read header line!");

	//Process sample names
	vector < string > tokens;
	vector < int > mappingS;
	stb.split(string(str.s), tokens);
	if (tokens.size() < 7) vrb.error("Incorrect number of columns!");
	for (int t = 6 ; t < tokens.size() ; t ++) mappingS.push_back(findSample(tokens[t]));
	vrb.bullet("#samples=" + stb.str(sample_count));

	//Read phenotypes
	int n_tested_variants = 0, n_ntested_variants = 0;
	vector < float > values = vector < float > (sample_count, 0.0);
	output_file fdhits (fout);
    if (fdhits.fail()) vrb.error("Cannot open file for writing [" + fout + "]");
	while (hts_getline(fp, KS_SEP_LINE, &str) >= 0) {
		stb.split(string(str.s), tokens);
		if (str.l && str.s[0] != tbx->conf.meta_char) {
			if (tokens.size() < 7) vrb.error("Incorrect number of columns!");

			string vid = tokens[3];
			string chr = tokens[0];
			string sta = tokens[1];
			string sto = tokens[2];

			if (filter_genotype.check(vid)) {
				for (int t = 6 ; t < tokens.size() ; t ++) {
					if (mappingS[t-6] >= 0) {
						if (tokens[t] == "NA") values[mappingS[t-6]] = bcf_float_missing;
						else values[mappingS[t-6]] = stof(tokens[t]);
					}
				}
				imputeValues(values);
				normalize(values);
				for (int p = 0 ; p < phenotype_count ; p ++) {
					double rcorr = fastCorrelation(phenotype_val[p], values);
					double npval = getNominalPvalue(abs(rcorr), sample_count - 2);
					fdhits << phenotype_id[p] << " " << chr << " " << sta << " " << sto << " " << vid << " " << npval << " " << rcorr << endl;
				}
				n_tested_variants ++;
			} else n_ntested_variants ++;
		}
	}

	//Finalize & verbose
	tbx_destroy(tbx);
	if (hts_close(fp)) vrb.error("Cannot properly close file");
	vrb.bullet(stb.str(n_tested_variants) + " variables tested");
	if (n_ntested_variants > 0) vrb.bullet(stb.str(n_ntested_variants) + " variables not tested");
	fdhits.close();
}







