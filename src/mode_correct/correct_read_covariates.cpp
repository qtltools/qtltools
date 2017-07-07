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

void correct_data::readCovariates(string fcov) {
	string buffer; vector < string > tokens;
	vector < int > mappingS;
	int n_includedS = 0, n_includedC = 0, n_excludedC = 0;

	vrb.title("Reading covariates in [" + fcov + "]");
	input_file fd (fcov);
	if (fd.fail()) vrb.error("Cannot open file!");

	//Read samples
	getline(fd, buffer);
	if (buffer.size() == 0) vrb.error("No header line detected!");
	stb.split(buffer, tokens);
	for (int t = 1 ; t < tokens.size() ; t ++) {
		mappingS.push_back(findSample(tokens[t]));
		if (mappingS.back() >= 0) n_includedS ++;
	}

	//Read covariates
	while (getline(fd, buffer)) {
		stb.split(buffer, tokens);
		if (tokens.size() < 2) vrb.error("Wrong Incorrect number of columns!");
		if (filter_covariate.check(tokens[0])) {
			covariate_val.push_back(vector < string > (sample_count, "0"));
			for (int t = 1 ; t < tokens.size() ; t ++) if (mappingS[t-1] >= 0) covariate_val.back()[mappingS[t-1]] = tokens[t];
            n_includedC ++;
		} else n_excludedC ++;
	}

	//Finalise
	covariate_count = n_includedC;
	vrb.bullet(stb.str(n_includedC) + " covariate(s) included");
	if (n_excludedC > 0) vrb.bullet(stb.str(n_excludedC) + " covariate(s) excluded by user");
	fd.close();
}


void correct_data::readQTLCovariates(string fqtl, string fvcf) {
	string buffer; vector < string > tokens;

	//File qtl
	map < string, vector < string > > map_qtl;
	map < string, vector < string > > :: iterator it_map_qtl;
	vrb.title("Reading QTLs in [" + fqtl + "]");
	input_file fd (fqtl);
	if (fd.fail()) vrb.error("Cannot open file!");
	int n_variant_covariates = 0;
	int n_variant_duplicates = 0;
	while (getline(fd, buffer)) {
		stb.split(buffer, tokens);
		if (tokens.size() != 2) vrb.error("Wrong Incorrect number of columns, expected 2!");
		string variant = tokens[1];
		string gene = tokens[0];

		it_map_qtl = map_qtl.find(variant);
		if (it_map_qtl == map_qtl.end()) {
			vector < string > value = vector < string >(1, gene);
			map_qtl.insert(pair < string, vector < string > > (tokens[1], value));
			n_variant_covariates ++;
		} else {
			it_map_qtl->second.push_back(gene);
			n_variant_duplicates ++;
		}
	}
	vrb.bullet(stb.str(n_variant_covariates) + " unique variants used as covariates");
	vrb.bullet(stb.str(n_variant_duplicates) + " unique variants act on multiple genes");
	fd.close();

	//
	for (int c = 0 ; c < covariate_count ; c ++) covariate_target.push_back(vector < string > (1, "ALL"));

	//File VCF
	vector < int > mappingS;
	int n_includedG = 0, n_excludedG = 0, n_includedS = 0;
	bcf_srs_t * sr =  bcf_sr_init();
	if(!(bcf_sr_add_reader (sr, fvcf.c_str()))) {
		switch (sr->errnum) {
		case not_bgzf: vrb.error("File not compressed with bgzip!"); break;
		case idx_load_failed: vrb.error("Impossible to load index file!"); break;
		case file_type_error: vrb.error("File format not detected by htslib!"); break;
		default : vrb.error("Unknown error!");
		}
	}
	int n_samples = bcf_hdr_nsamples(sr->readers[0].header);
	for (int i0 = 0 ; i0 < n_samples ; i0 ++) {
		mappingS.push_back(findSample(string(sr->readers[0].header->samples[i0])));
		if (mappingS.back() >= 0) n_includedS++;
	}
    unsigned int linecount=0;
	int ngt, ngt_arr = 0, nds, nds_arr = 0, * gt_arr = NULL, nsl, nsl_arr = 0, * sl_arr = NULL;
	float * ds_arr = NULL;
	bcf1_t * line;
	while(bcf_sr_next_line (sr)) {
        linecount ++;
        if (linecount % 100000 == 0) vrb.bullet("Read " + stb.str(linecount) + " lines");
		line =  bcf_sr_get_line(sr, 0);
		if (line->n_allele == 2) {
			ngt = bcf_get_genotypes(sr->readers[0].header, line, &gt_arr, &ngt_arr);
			nds = bcf_get_format_float(sr->readers[0].header, line,"DS", &ds_arr, &nds_arr);
			if (nds == n_samples || ngt == 2*n_samples) {
				bcf_unpack(line, BCF_UN_STR);
				string sid = string(line->d.id);

				it_map_qtl = map_qtl.find(sid);

				if (it_map_qtl != map_qtl.end() && filter_covariate.check(sid)) {
					covariate_val.push_back(vector < string > (sample_count, "0"));
					for(int i = 0 ; i < n_samples ; i ++) {
						if (mappingS[i] >= 0) {
							if (nds > 0) covariate_val.back()[mappingS[i]] = stb.str(ds_arr[i]);
							else {
								if (gt_arr[2*i+0] == bcf_gt_missing || gt_arr[2*i+1] == bcf_gt_missing) covariate_val.back()[mappingS[i]] ="NA";
								else covariate_val.back()[mappingS[i]] = stb.str(bcf_gt_allele(gt_arr[2*i+0]) + bcf_gt_allele(gt_arr[2*i+1]));
							}
						}
					}
					covariate_target.push_back(vector < string >());
					for (int g = 0 ; g  < it_map_qtl->second.size() ; g ++) covariate_target.back().push_back(it_map_qtl->second[g]);
					covariate_count++;
					n_includedG++;
				} else n_excludedG ++;
			}
		}
	}

	//Finalize
	free(gt_arr);
	free(ds_arr);
	bcf_sr_destroy(sr);
	vrb.bullet(stb.str(n_includedG) + " variants included");
	if (n_excludedG > 0) vrb.bullet(stb.str(n_excludedG) + " variants excluded");
}
