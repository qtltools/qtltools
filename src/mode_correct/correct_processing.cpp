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

#include "correct_data.h"

void correct_data::processBED(string fin, string fout) {
	int n_includedP = 0, n_excludedP_user = 0;

	vrb.title("Initialize residualizer");
	if (covariate_target.size() == 0) {
		covariate_engine = new residualizer (sample_count);
		for (int c = 0 ; c < covariate_count ; c ++) covariate_engine->push(covariate_val[c]);
		covariate_engine->build();
		vrb.bullet("#covariates = " + stb.str(covariate_count));
	}

	//Open BED file
	vrb.title("Reading phenotype data in [" + fin + "] and writing [" + fout + "]");
	output_file fdo (fout.c_str());
	if (fdo.fail()) vrb.error("CAnnot open file for writing!");

	htsFile *fp = hts_open(fin.c_str(),"r");
	if (!fp) vrb.error("Cannot open file for reading!");
	kstring_t str = {0,0,0};
	if (hts_getline(fp, KS_SEP_LINE, &str) <= 0 || !str.l || str.s[0] != '#' ) vrb.error("Cannot read header!");
	fdo << "#chr\tstart\tend\tid\tgid\tstrd";

	//Read and map sample names
	vector < int > mappingS;
	vector < string > tokens;
	stb.split(string(str.s), tokens);
	if (tokens.size() < 7) vrb.error("Incorrect number of columns!");
	for (int t = 6 ; t < tokens.size() ; t ++) {
		mappingS.push_back(findSample(tokens[t]));
		//if (mappingS.back() >= 0) fdo << "\t" << tokens[t];
	}
	for (int i = 0 ; i < sample_count ; i ++) fdo << "\t" << sample_id[i];
	fdo << endl;

	//Read phenotypes
	while (hts_getline(fp, KS_SEP_LINE, &str) >= 0) {
		stb.split(string(str.s), tokens);
		if (filter_phenotype.check(tokens[3])) {
			vector < float > values = vector < float > (sample_count, 0.0);
			fdo << tokens[0] << "\t" << tokens[1] << "\t" << tokens[2] << "\t" << tokens[3] << "\t" << tokens[4] << "\t" << tokens[5];
			for (int t = 6 ; t < tokens.size() ; t ++) if (mappingS[t-6] >= 0) values[mappingS[t-6]] = ((tokens[t] != "NA")?stof(tokens[t]):bcf_float_missing);
			imputeMissing(values);
			if (residualize) {
				if (covariate_target.size() > 0) {
					covariate_engine = new residualizer (sample_count);
					int all_covariates = 0;
					int qtl_covariates = 0;
					for (int c = 0 ; c < covariate_count ; c ++) {
						for (int g = 0 ; g < covariate_target[c].size() ; g ++) {
							if (covariate_target[c][g] == "ALL") {
								covariate_engine->push(covariate_val[c]);
								all_covariates ++;
							}
							if (covariate_target[c][g] == tokens[3]) {
								covariate_engine->push(covariate_val[c]);
								qtl_covariates ++;
							}
						}
					}
					covariate_engine->build();
					vrb.bullet (tokens[3] + " residualized for #cov_common = " + stb.str(all_covariates) + "\t#cov_qtl=" + stb.str(qtl_covariates));
				}
				covariate_engine->residualize(values);
				if (covariate_target.size() > 0) {
					delete covariate_engine;
					covariate_engine = NULL;
				}
			}
			if (normalize) normalTransform(values);
			for (int s = 0 ; s < sample_count ; s++) fdo << "\t" << values[s];
			n_includedP ++ ;
			fdo << endl;
		} else n_excludedP_user ++;
	}

	//Finalize & verbose
	vrb.bullet(stb.str(n_includedP) + " phenotypes");
	if (n_excludedP_user > 0) vrb.bullet(stb.str(n_excludedP_user) + " phenotypes excluded by user");
	hts_close(fp);
	fdo.close();
}

void correct_data::processVCF(string fin, string fout) {
	int n_includedG = 0, n_excludedG_user = 0, n_excludedG_miss = 0, n_excludedG_mult = 0;

	if (covariate_target.size() == 0) {
		covariate_engine = new residualizer (sample_count);
		for (int c = 0 ; c < covariate_count ; c ++) covariate_engine->push(covariate_val[c]);
		covariate_engine->build();
		vrb.bullet("#covariates = " + stb.str(covariate_count));
	}

	vrb.title("Reading genotype data in [" + fin + "] and writing [" + fout + "]");
	bcf_sweep_t * sw = bcf_sweep_init(fin.c_str());
	if (!sw) vrb.error("Cannot open file for reading [" + fin + "]");
	bcf_hdr_t * hdr_old  = bcf_sweep_hdr(sw);
	if (!hdr_old) vrb.error("Cannot read header!");

	htsFile * fp = hts_open(fout.c_str(),"wz");
	if (!fp) vrb.error("Cannot open file for writing [" + fout + "]");

	//Update sample ids in hdr
	int * imap = (int *) malloc(sample_count * sizeof(int));
	char ** samples = (char **) malloc(sample_count * sizeof(char *));
	for (int i = 0 ; i < sample_count ; i ++) samples[i] = hdr_old->samples[bcf_hdr_id2int(hdr_old, BCF_DT_SAMPLE, sample_id[i].c_str())];
	bcf_hdr_t * hdr_new = bcf_hdr_subset(hdr_old, sample_count, samples, imap);
	bcf_hdr_write(fp, hdr_new);

	//Read the VCF
	bcf1_t * rec = bcf_init1();
	int mDS = 0; float * vDS = NULL;
    while ( (rec = bcf_sweep_fwd(sw)) ) {
    	bcf_subset(hdr_old, rec, sample_count, imap);
    	bcf_unpack(rec, BCF_UN_STR);
		if (rec->n_allele != 2) n_excludedG_mult++;
		else if (bcf_get_format_float(hdr_new, rec, "DS", &vDS, &mDS) != sample_count) n_excludedG_miss++;
		else if (!filter_genotype.check(string(rec->d.id))) n_excludedG_user++;
		else {
			imputeMissing(vDS);
			if (residualize) covariate_engine->residualize(vDS);
			if (normalize) normalTransform(vDS);
			bcf_update_format_float(hdr_new, rec, "DS", vDS, sample_count);
			bcf_write1(fp, hdr_new, rec);
			n_includedG ++;
		}
	}
    free(vDS);
    bcf_sweep_destroy(sw);
	bcf_destroy1(rec);
	bcf_hdr_destroy(hdr_old);
	bcf_hdr_destroy(hdr_new);
	if (hts_close(fp)) vrb.error("Cannot close properly file!");
	vrb.bullet(stb.str(n_includedG) + " variants corrected");
	if (n_excludedG_user) vrb.bullet(stb.str(n_excludedG_user) + " variants excluded by user");
	if (n_excludedG_miss) vrb.bullet(stb.str(n_excludedG_miss) + " missing DS variants excluded");
	if (n_excludedG_mult) vrb.bullet(stb.str(n_excludedG_mult) + " multi-allelic variants excluded");
}
