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

#include "pca_data.h"

void pca_data::readData(string filename) {
    vrb.title("Reading genotype data in [" + filename + "]");
    htsFile * fp = hts_open(filename.c_str(),"r");
    enum htsExactFormat fileformat = fp->format.format;
    hts_close(fp);
    if (fileformat == bcf) {
        vrb.bullet("File format detected: BCF");
        readDataVCF(filename);
    } else if (fileformat == vcf) {
        vrb.bullet("File format detected: VCF");
        readDataVCF(filename);
    } else if (fileformat == sam) {
        vrb.bullet("File format detected: BED");
        readDataBED(filename);
    } else vrb.error("File format not supported!");
}

void pca_data::readDataVCF(string fvcf) {
    int n_includedG = 0;
    int n_excludedG_mult = 0;
    int n_excludedG_void = 0;
    int n_excludedG_user = 0;
    int n_includedS = 0;
    vector < int > mappingS;
    
    //Opening files
    bcf_srs_t * sr =  bcf_sr_init();
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
    for (int i0 = 0 ; i0 < n_samples ; i0 ++) {
        mappingS.push_back(findSample(string(sr->readers[0].header->samples[i0])));
        if (mappingS.back() >= 0) n_includedS++;
    }
    
    PCA._xXf.resize(sample_count, __RESIZE_CHUNK__);
    
    //Read genotype data
    //int ngt, ngt_arr = 0, nds, nds_arr = 0, * gt_arr = NULL, nsl, nsl_arr = 0, * sl_arr = NULL;
    int ngt, ngt_arr = 0, nds, nds_arr = 0, * gt_arr = NULL;
    float * ds_arr = NULL;
    bcf1_t * line;
    unsigned int linecount = 0;
    string pChr = "";
    int pPos = 0;
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
                string chr  = string(bcf_hdr_id2name(sr->readers[0].header, line->rid));
                int pos = line->pos + 1;
                if (filter_genotype.check(sid) && ((pChr == chr && abs(pos - pPos) >= distance_separator) || pChr != chr)) {
                    vector < float > temp(sample_count, 0.0);
                    int total = 0 ;
                    int count = 0;
                    for(int i = 0 ; i < n_samples ; i ++) {
                        if (mappingS[i] >= 0) {
                            if (nds > 0) {
                            	temp[mappingS[i]] = ds_arr[i];
        						count+=2;
        						total+=temp[mappingS[i]];
                            } else {
                                if (gt_arr[2*i+0] == bcf_gt_missing || gt_arr[2*i+1] == bcf_gt_missing) temp[mappingS[i]] = bcf_float_missing;
                                else {
                                	temp[mappingS[i]] = bcf_gt_allele(gt_arr[2*i+0]) + bcf_gt_allele(gt_arr[2*i+1]);
            						count+=2;
            						total+=temp[mappingS[i]];
                                }
                            }
                        }
                    }
        			double af = (double) total / (double) count;
        			if (maf_cutoff > af || 1.0-maf_cutoff < af){
        				n_excludedG_user ++;
        				continue;
        			}
                    pChr = chr;
                    pPos = pos;
                    //data_id.push_back(sid);
                    //data_chr.push_back(chr);
                    //string genotype_ref = string(line->d.allele[0]);
                    //data_start.push_back(pos);
                    //nsl = bcf_get_info_int32(sr->readers[0].header, line, "END", &sl_arr, &nsl_arr);
                    //if (nsl >= 0 && nsl_arr == 1) data_end.push_back(sl_arr[0]);
                    //else data_end.push_back(data_start.back() + genotype_ref.size() - 1);
                    if (n_includedG >= PCA._xXf.cols()) resizeData();
                    for (int i = 0 ; i < temp.size(); i++) PCA._xXf(i,n_includedG) = temp[i];
                    n_includedG++;
                } else n_excludedG_user ++;
            } else n_excludedG_void ++;
        } else n_excludedG_mult ++;
    }
    
    //Finalize
    free(gt_arr);
    free(ds_arr);
    bcf_sr_destroy(sr);
    data_count = n_includedG;
    vrb.bullet(stb.str(n_includedG) + " variants included");
    if (n_excludedG_user > 0) vrb.bullet(stb.str(n_excludedG_user) + " variants excluded by user");
    if (n_excludedG_mult > 0) vrb.bullet(stb.str(n_excludedG_mult) + " multi-allelic variants excluded");
    if (n_excludedG_void > 0) vrb.bullet(stb.str(n_excludedG_void) + " uninformative variants excluded [no GT/DS]");
    if (data_count == 0) vrb.leave("Cannot find genotypes in target region!");
    finalizeData(n_includedG);
}

void pca_data::readDataBED(string fbed) {
	string buffer;
	int n_includedG = 0;
	int n_excludedG_user = 0;
	int n_includedS = 0;
	int n_excludedS = 0;
	int n_missingS = 0;
	vector < int > mappingS;

	//Opening files
	htsFile *fp = hts_open(fbed.c_str(),"r");
	if (!fp) vrb.error("Cannot open file!");
	tbx_t * tbx = tbx_index_load(fbed.c_str());
	if (!tbx) vrb.error("Cannot load index file!");
	kstring_t str = {0,0,0};
	if (hts_getline(fp, KS_SEP_LINE, &str) <= 0 || !str.l || str.s[0] != tbx->conf.meta_char ) vrb.error("Cannot read header line!");

	//Process sample names
	vector < string > tokens;
	stb.split(string(str.s), tokens);
	if (tokens.size() < 7) vrb.error("Incorrect number of columns!");
	for (int i0 = 6 ; i0 < tokens.size() ; i0 ++) {
		string sid = tokens[i0];
		if (filter_sample.check(sid)) {
			mappingS.push_back(findSample(sid));
			if (mappingS.back() >= 0) n_includedS ++;
			else n_missingS ++;
		} else {
			mappingS.push_back(-1);
			n_excludedS ++;
		}
	}
	vrb.bullet(stb.str(n_includedS) + " samples included");
	if (n_excludedS > 0) vrb.bullet(stb.str(n_excludedS) + " samples excluded by user");
	if (n_missingS > 0) vrb.bullet(stb.str(n_missingS) + " samples without phenotype data");
	if (n_includedS != sample_count) vrb.error("Cannot find genotype for " + stb.str(sample_count - n_includedS) + " samples!");

    PCA._xXf.resize(sample_count,__RESIZE_CHUNK__);

    unsigned int linecount = 0;
    string pChr = "";
    int pPos = 0;
	while (hts_getline(fp, KS_SEP_LINE, &str) >= 0) {
		linecount ++;
		if (linecount % 100000 == 0) vrb.bullet("Read " + stb.str(linecount) + " lines");
		stb.split(string(str.s), tokens);
		if (tokens.size() < 7) vrb.error("Incorrect number of columns!");
        string chr  = tokens[0];
        int pos = atoi(tokens[1].c_str()) + 1;
		if (filter_genotype.check(tokens[3]) && ((pChr == chr && abs(pos - pPos) >= distance_separator) || pChr != chr)) {
			vector < float > temp(sample_count, 0.0);
            int total = 0 ;
            int count = 0;
			for (int t = 6 ; t < tokens.size() ; t ++) {
				if (mappingS[t-6] >= 0) {
					if (tokens[t] == "NA") temp[mappingS[t-6]] = bcf_float_missing;
					else {
						temp[mappingS[t-6]] = stof(tokens[t]);
						count+=2;
						total+=temp[mappingS[t-6]];
					}
				}
			}
			double af = (double) total / (double) count;
			if (maf_cutoff > af || 1.0-maf_cutoff < af){
				n_excludedG_user ++;
				continue;
			}
			//data_id.push_back(tokens[3]);
			//data_chr.push_back(chr);
			//data_start.push_back(pos);
			//data_end.push_back(atoi(tokens[2].c_str()));
            if (n_includedG >= PCA._xXf.cols()) resizeData();
            for (int i = 0 ; i < temp.size(); i++) PCA._xXf(i,n_includedG) = temp[i];
			n_includedG++;
		} else n_excludedG_user ++;
	}


	//Finalize & verbose
	tbx_destroy(tbx);
	if (hts_close(fp)) vrb.error("Cannot properly close file!");
	data_count = n_includedG;
	vrb.bullet(stb.str(n_includedG) + " variants included");
	if (n_excludedG_user > 0) vrb.bullet(stb.str(n_excludedG_user) + " variants excluded by user");
    if (data_count == 0) vrb.leave("Cannot find variants in target region!");
    finalizeData(n_includedG);
}

void pca_data::readDataPhenoBED(string fbed) {
    int n_includedS = 0, n_excludedS = 0, n_excludedU = 0, n_excludedP = 0, n_includedP = 0;
    vector < int > mappingS;

    //Open BED file
    vrb.title("Reading phenotype data in [" + fbed + "]");
    htsFile *fp = hts_open(fbed.c_str(),"r");
    if (!fp) vrb.error("Cannot open file!");
    tbx_t * tbx = tbx_index_load(fbed.c_str());
    if (!tbx) vrb.error("Cannot open index file");
    kstring_t str = {0,0,0};
    if (hts_getline(fp, KS_SEP_LINE, &str) <= 0 || !str.l || str.s[0] != '#' ) vrb.error("Cannot read header line");

    //Process sample names
    vector < string > tokens;
    stb.split(string(str.s), tokens);
    if (tokens.size() < 7) vrb.error("Incorrect number of columns!");
    for (int t = 6 ; t < tokens.size() ; t ++) {
        if (filter_sample.check(tokens[t])) {
            mappingS.push_back(findSample(tokens[t]));
            if (mappingS.back() < 0) n_excludedS ++;
            else n_includedS ++;
        } else n_excludedU ++;
    }
    vrb.bullet(stb.str(n_includedS) + " samples included");
    if (n_excludedU > 0) vrb.bullet(stb.str(n_excludedU) + " samples excluded by user");
    if (n_excludedS > 0) vrb.bullet(stb.str(n_excludedS) + " samples without genotype data");
    if (n_includedS != sample_count) vrb.error("Cannot find phenotype data for " + stb.str(sample_count - n_includedS) + " samples!");
    PCA._xXf.resize(sample_count, __RESIZE_CHUNK__);
    unsigned long int linecount = 1;
    //Read phenotypes
    while (hts_getline(fp, KS_SEP_LINE, &str) >= 0) {
        if (str.l && str.s[0] != tbx->conf.meta_char) {
            stb.split(string(str.s), tokens);
            if (filter_phenotype.check(tokens[3])) {
				//data_id.push_back(tokens[3]);
				//data_chr.push_back(tokens[0]);
				//data_start.push_back(atoi(tokens[1].c_str()) + 1);
				//data_end.push_back(atoi(tokens[2].c_str()));
                if (n_includedP >= PCA._xXf.cols()) resizeData();
				for (int t = 6 ; t < tokens.size() ; t ++) if (mappingS[t-6] >= 0) {
					if (tokens[t] == "NA") PCA._xXf(mappingS[t-6],n_includedP) = bcf_float_missing;
					else PCA._xXf(mappingS[t-6],n_includedP) = stof(tokens[t]);
				}
                linecount++;
                n_includedP++;
            } else n_excludedP ++;
        }
    }

    //Finalize & verbose
    tbx_destroy(tbx);
    data_count = n_includedP;
    vrb.bullet(stb.str(data_count) + " phenotypes included");
    if (n_excludedP > 0) vrb.bullet(stb.str(n_excludedP) + " phenotypes excluded by user");
    if (data_count == 0) vrb.leave("Cannot find phenotypes in target region!");
    hts_close(fp);
    finalizeData(n_includedP);
}


