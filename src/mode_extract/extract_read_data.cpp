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

#include "extract_data.h"

void extract_data::readBED(string fbed) {
	int n_includedP = 0;
	int n_excludedP = 0;
	vector < int > mappingS;

	//Open BED file
	vrb.title("Reading phenotype data in [" + fbed + "]");
	htsFile *fp = hts_open(fbed.c_str(),"r");
	if (!fp) vrb.error("Cannot open file");
	tbx_t *tbx = tbx_index_load(fbed.c_str());
	if (!tbx) vrb.error("Cannot open index file");
	kstring_t str = {0,0,0};
	if (hts_getline(fp, KS_SEP_LINE, &str) <= 0 || !str.l || str.s[0] != tbx->conf.meta_char ) vrb.error("Cannot read header line!");

	//Process sample names
	vector < string > tokens;
	stb.split(string(str.s), tokens);
	if (tokens.size() < 7) vrb.error("Incorrect number of columns!");
	for (int t = 6 ; t < tokens.size() ; t ++) mappingS.push_back(findSample(tokens[t]));

    //Read phenotypes
    if (regionData.chr != "NA"){
        hts_itr_t *itr = tbx_itr_querys(tbx, regionData.get().c_str());
        vrb.bullet("target region [" + regionData.get() + "]");
        if (!itr) vrb.error("Cannot jump to region!");

        //Read data
        while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
            stb.split(string(str.s), tokens);
            if (tokens.size() < 7) vrb.error("Incorrect number of columns!");
            if (filter_phenotype.check(tokens[3])) {
                variable_id.push_back(tokens[3]);
                variable_chr.push_back(tokens[0]);
                variable_start.push_back(atoi(tokens[1].c_str()) + 1);
                variable_end.push_back(atoi(tokens[2].c_str()));
                variable_val.push_back(vector < string > (sample_count));
                for (int t = 6 ; t < tokens.size() ; t ++) if (mappingS[t-6] >= 0) variable_val.back()[mappingS[t-6]] = tokens[t];
                n_includedP++;
            } else n_excludedP ++;
        }
        tbx_itr_destroy(itr);
    }else{
        while (hts_getline(fp, KS_SEP_LINE, &str) >= 0) {
            stb.split(string(str.s), tokens);
            if (str.l && str.s[0] != tbx->conf.meta_char) {
                if (tokens.size() < 5) vrb.error("Incorrect number of columns!");
                if (filter_phenotype.check(tokens[3])) {
                    variable_id.push_back(tokens[3]);
                    variable_chr.push_back(tokens[0]);
                    variable_start.push_back(atoi(tokens[1].c_str()) + 1);
                    variable_end.push_back(atoi(tokens[2].c_str()));
                    variable_val.push_back(vector < string > (sample_count));
                    for (int t = 6 ; t < tokens.size() ; t ++) if (mappingS[t-6] >= 0) variable_val.back()[mappingS[t-6]] = tokens[t];
                    n_includedP++;
                } else n_excludedP ++;
            }
        }
    }

	//Finalize & verbose
	tbx_destroy(tbx);
	if (hts_close(fp)) vrb.error("Cannot properly close file");
	vrb.bullet(stb.str(n_includedP) + " phenotypes included");
	if (n_excludedP > 0) vrb.bullet(stb.str(n_excludedP) + " phenotypes excluded by user");
    if (n_includedP == 0) vrb.warning("Cannot find phenotypes to extract!");
}

void extract_data::readVCF(string fvcf) {
	int n_includedG = 0;
	int n_excludedG_mult = 0;
	int n_excludedG_void = 0;
	int n_excludedG_user = 0;
	vector < int > mappingS;

	//Opening files
	vrb.title("Reading genotype data in [" + fvcf + "]");
	bcf_srs_t * sr =  bcf_sr_init();
    if (regionData.chr != "NA"){
        vrb.bullet("target region [" + regionData.get() + "]");
        if (bcf_sr_set_regions(sr, regionData.get().c_str(), 0) == -1) vrb.error("Cannot jump to region!");
    }
	if(!(bcf_sr_add_reader (sr, fvcf.c_str()))) {
		switch (sr->errnum) {
		case not_bgzf: vrb.error("File not compressed with bgzip!");
		case idx_load_failed: vrb.error("Impossible to load index file!");
		case file_type_error: vrb.error("File format not detected by htslib!");
		default : vrb.error("Unknown error!");
		}
	}

	//Sample processing
	int n_samples = bcf_hdr_nsamples(sr->readers[0].header);
	for (int i0 = 0 ; i0 < n_samples ; i0 ++) mappingS.push_back(findSample(string(sr->readers[0].header->samples[i0])));

	//Read genotype data
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
				if (filter_genotype.check(sid)) {
					variable_id.push_back(sid);
					variable_chr.push_back(string(bcf_hdr_id2name(sr->readers[0].header, line->rid)));
					string genotype_ref = string(line->d.allele[0]);
					variable_start.push_back(line->pos + 1);
					nsl = bcf_get_info_int32(sr->readers[0].header, line, "END", &sl_arr, &nsl_arr);
					if (nsl >= 0 && nsl_arr == 1) variable_end.push_back(sl_arr[0]);
					else variable_end.push_back(variable_start.back() + genotype_ref.size() - 1);
					variable_val.push_back(vector < string > (sample_count));
					for(int i = 0 ; i < n_samples ; i ++) {
						if (mappingS[i] >= 0) {
							if (nds > 0) variable_val.back()[mappingS[i]] = stb.str(ds_arr[i]);
							else {
								if (gt_arr[2*i+0] == bcf_gt_missing || gt_arr[2*i+1] == bcf_gt_missing) variable_val.back()[mappingS[i]] = "NA";
								else variable_val.back()[mappingS[i]] = stb.str(bcf_gt_allele(gt_arr[2*i+0]) + bcf_gt_allele(gt_arr[2*i+1]));
							}
						}
					}
					n_includedG++;
				} else n_excludedG_user ++;
			} else n_excludedG_void ++;
		} else n_excludedG_mult ++;
	}

	//Finalize
	free(gt_arr);
	free(ds_arr);
	bcf_sr_destroy(sr);
	vrb.bullet(stb.str(n_includedG) + " variants included");
	if (n_excludedG_user > 0) vrb.bullet(stb.str(n_excludedG_user) + " variants excluded by user");
	if (n_excludedG_mult > 0) vrb.bullet(stb.str(n_excludedG_mult) + " multi-allelic variants excluded");
	if (n_excludedG_void > 0) vrb.bullet(stb.str(n_excludedG_void) + " uninformative variants excluded [no GT/DS]");
    if (n_includedG == 0) vrb.leave("Cannot find genotypes to extract!");
}

void extract_data::readCOV(string fcov) {
	string buffer;
	vector < string > str;
	int n_includedC = 0;
	int n_excludedC = 0;
	vector < int > mappingS;

	vrb.title("Reading covariates in [" + fcov + "]");
	input_file fd (fcov);
	if (fd.fail()) vrb.error("Cannot open file!");

	//Read samples
	getline(fd, buffer);
	if (buffer.size() == 0) vrb.error("No header line detected!");
	stb.split(buffer, str	);
	for (int t = 1 ; t < str.size() ; t ++) mappingS.push_back(findSample(str[t]));

	//Read covariates
	while(getline(fd, buffer)) {
		stb.split(buffer, str);
		if (str.size() < 2) vrb.error("Incorrect number of columns!");
		if (filter_covariate.check(str[0])) {
            variable_id.push_back(str[0]);
            variable_chr.push_back(string("NA"));
            variable_start.push_back(-1);
            variable_end.push_back(-1);
            variable_val.push_back(vector < string > (sample_count));
			for (int t = 1 ; t < str.size() ; t ++) if (mappingS[t-1] >= 0) variable_val.back()[mappingS[t-1]] = str[t];
			n_includedC ++;
		} else n_excludedC ++;
	}

	//Finalise
	vrb.bullet(stb.str(n_includedC) + " covariates included");
	if (n_excludedC > 0) vrb.bullet(stb.str(n_excludedC) + " covariates excluded");
	fd.close();
}
