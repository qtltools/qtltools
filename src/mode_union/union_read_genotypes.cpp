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

#include "union_data.h"

void union_data::readGenotypes(string filename, string region) {
	//vrb.title("Reading genotype data in [" + filename + "]");
    htsFile * fp = hts_open(filename.c_str(),"r");
    enum htsExactFormat fileformat = fp->format.format;
    hts_close(fp);
    if (fileformat == bcf) {
    	//vrb.bullet("File format detected: BCF");
    	readGenotypesVCF(filename,region);
    } else if (fileformat == vcf) {
    	//vrb.bullet("File format detected: VCF");
    	readGenotypesVCF(filename,region);
    } else if (fileformat == sam) {
    	//vrb.bullet("File format detected: BED");
    	readGenotypesBED(filename,region);
    } else vrb.error("File format not supported!");
}

void union_data::scanGenotypes(string filename) {
	vrb.title("Scanning genotype data in [" + filename + "]");
    htsFile * fp = hts_open(filename.c_str(),"r");
    enum htsExactFormat fileformat = fp->format.format;
    hts_close(fp);
    if (fileformat == bcf) {
    	vrb.bullet("File format detected: BCF");
    	scanGenotypesVCF(filename);
    } else if (fileformat == vcf) {
    	vrb.bullet("File format detected: VCF");
    	scanGenotypesVCF(filename);
    } else if (fileformat == sam) {
    	vrb.bullet("File format detected: BED");
    	scanGenotypesBED(filename);
    } else vrb.error("File format not supported!");
}


void union_data::readGenotypesVCF(string fvcf,string region) {
	int n_includedG = 0;
	int n_excludedG_mult = 0;
	int n_excludedG_void = 0;
	int n_excludedG_user = 0;
	int n_includedS = 0;
	vector < int > mappingS;
	genotype_id.clear();
	genotype_chr.clear();
	genotype_start.clear();
	genotype_end.clear();
	genotype_val.clear();
	genotype_count=0;
	genotype_id_to_idx.clear();

	//Opening files
	bcf_srs_t * sr =  bcf_sr_init();

    //vrb.bullet("target region [" + regionGenotype.get() + "]");
    //if (bcf_sr_set_regions(sr, regionGenotype.get().c_str(), 0) == -1) vrb.error("Cannot jump to region!");
	bcf_sr_set_regions(sr, region.c_str(), 0);
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


	//Read genotype data
	int ngt, ngt_arr = 0, nds, nds_arr = 0, * gt_arr = NULL, nsl, nsl_arr = 0, * sl_arr = NULL;
	float * ds_arr = NULL;
	bcf1_t * line;
    unsigned int linecount = 0;
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
					genotype_id.push_back(sid);
					genotype_chr.push_back(string(bcf_hdr_id2name(sr->readers[0].header, line->rid)));
					string genotype_ref = string(line->d.allele[0]);
					genotype_start.push_back(line->pos + 1);
					nsl = bcf_get_info_int32(sr->readers[0].header, line, "END", &sl_arr, &nsl_arr);
					if (nsl >= 0 && nsl_arr == 1) genotype_end.push_back(sl_arr[0]);
					else genotype_end.push_back(genotype_start.back() + genotype_ref.size() - 1);
					genotype_val.push_back(vector < float > (sample_count, 0.0));

					for(int i = 0 ; i < n_samples ; i ++) {
						if (mappingS[i] >= 0) {
							if (nds > 0) genotype_val.back()[mappingS[i]] = ds_arr[i];
							else {
								if (gt_arr[2*i+0] == bcf_gt_missing || gt_arr[2*i+1] == bcf_gt_missing) genotype_val.back()[mappingS[i]] = bcf_float_missing;
								else genotype_val.back()[mappingS[i]] = bcf_gt_allele(gt_arr[2*i+0]) + bcf_gt_allele(gt_arr[2*i+1]);
							}
						}
					}
                    pair < string, int > temp (sid,n_includedG);
                    genotype_id_to_idx.insert(temp);
					n_includedG++;
				} else n_excludedG_user ++;
			} else n_excludedG_void ++;
		} else n_excludedG_mult ++;
	}

	//Finalize
	free(gt_arr);
	free(ds_arr);
	bcf_sr_destroy(sr);
	genotype_count = n_includedG;
	//vrb.bullet(stb.str(n_includedG) + " variants included");
	//if (n_excludedG_user > 0) vrb.bullet(stb.str(n_excludedG_user) + " variants excluded by user");
	//if (n_excludedG_mult > 0) vrb.bullet(stb.str(n_excludedG_mult) + " multi-allelic variants excluded");
	//if (n_excludedG_void > 0) vrb.bullet(stb.str(n_excludedG_void) + " uninformative variants excluded [no GT/DS]");
    //if (genotype_count == 0) vrb.leave("Cannot find genotypes in target region!");
}

void union_data::readGenotypesBED(string fbed,string region) {
	string buffer;
	int n_includedG = 0;
	int n_excludedG_user = 0;
	int n_includedS = 0;
	int n_excludedS = 0;
	int n_missingS = 0;
	vector < int > mappingS;
	genotype_id.clear();
	genotype_chr.clear();
	genotype_start.clear();
	genotype_end.clear();
	genotype_val.clear();
	genotype_count=0;
	genotype_id_to_idx.clear();
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
	//vrb.bullet(stb.str(n_includedS) + " samples included");
	//if (n_excludedS > 0) vrb.bullet(stb.str(n_excludedS) + " samples excluded by user");
	//if (n_missingS > 0) vrb.bullet(stb.str(n_missingS) + " samples without phenotype data");
	//if (n_includedS != sample_count) vrb.error("Cannot find genotype for " + stb.str(sample_count - n_includedS) + " samples!");

    unsigned int linecount = 0;

	//Jump to interesting region

	hts_itr_t *itr = tbx_itr_querys(tbx, region.c_str());
	//vrb.bullet("target region [" + regionGenotype.get() + "]");
	//if (!itr) vrb.error("Cannot jump to region!");
	while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
		linecount ++;
		if (linecount % 100000 == 0) vrb.bullet("Read " + stb.str(linecount) + " lines");
		stb.split(string(str.s), tokens);
		if (tokens.size() < 7) vrb.error("Incorrect number of columns!");
		if (filter_genotype.check(tokens[3])) {
			genotype_id.push_back(tokens[3]);
			genotype_chr.push_back(tokens[0]);
			genotype_start.push_back(atoi(tokens[1].c_str()) + 1);
			genotype_end.push_back(atoi(tokens[2].c_str()));
			genotype_val.push_back(vector < float > (sample_count, 0.0));
			for (int t = 6 ; t < tokens.size() ; t ++) {
				if (mappingS[t-6] >= 0) {
					if (tokens[t] == "NA") genotype_val.back()[mappingS[t-6]] = bcf_float_missing;
					else genotype_val.back()[mappingS[t-6]] = stof(tokens[t]);
				}
			}
			pair < string, int > temp (tokens[3],n_includedG);
			genotype_id_to_idx.insert(temp);
			n_includedG++;
		} else n_excludedG_user ++;
	}
	tbx_itr_destroy(itr);


	//Finalize & verbose
	tbx_destroy(tbx);
	if (hts_close(fp)) vrb.error("Cannot properly close file!");
	genotype_count = n_includedG;
	//vrb.bullet(stb.str(n_includedG) + " variants included");
	//if (n_excludedG_user > 0) vrb.bullet(stb.str(n_excludedG_user) + " variants excluded by user");
    //if (genotype_count == 0) vrb.leave("Cannot find variants in target region!");
}


void union_data::scanGenotypesVCF(string fvcf) {
    int n_includedG = 0;
    int n_excludedG_mult = 0;
    int n_excludedG_void = 0;
    int n_excludedG_user = 0;
    int n_includedS = 0;
    vector < int > mappingS;

    //Opening files
    bcf_srs_t * sr =  bcf_sr_init();
    
    if ( regionGenotype.chr != "NA"){
        vrb.bullet("target region [" + regionGenotype.get() + "]");
        if (bcf_sr_set_regions(sr, regionGenotype.get().c_str(), 0) == -1) vrb.error("Cannot jump to region!");
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
    for (int i0 = 0 ; i0 < n_samples ; i0 ++) {
        mappingS.push_back(findSample(string(sr->readers[0].header->samples[i0])));
        if (mappingS.back() >= 0) n_includedS++;
    }

    //Read genotype data
    int ngt, ngt_arr = 0, nds, nds_arr = 0, * gt_arr = NULL, nsl, nsl_arr = 0, * sl_arr = NULL;
    float * ds_arr = NULL;
    bcf1_t * line;
    unsigned int linecount = 0;
    while(bcf_sr_next_line (sr)) {
        linecount ++;
        if (linecount % 1000000 == 0) vrb.bullet("Read " + stb.str(linecount) + " lines");
        line =  bcf_sr_get_line(sr, 0);
        if (line->n_allele == 2) {
            ngt = bcf_get_genotypes(sr->readers[0].header, line, &gt_arr, &ngt_arr);
            nds = bcf_get_format_float(sr->readers[0].header, line,"DS", &ds_arr, &nds_arr);
            if (nds == n_samples || ngt == 2*n_samples) {
                bcf_unpack(line, BCF_UN_STR);
                string sid = string(line->d.id);
                if (genotype_id_to_idx.count(sid)) continue;
                if (filter_genotype.check(sid)) {
                    genotype_id.push_back(sid);
                    genotype_chr.push_back(string(bcf_hdr_id2name(sr->readers[0].header, line->rid)));
                    string genotype_ref = string(line->d.allele[0]);
                    genotype_start.push_back(line->pos + 1);
                    nsl = bcf_get_info_int32(sr->readers[0].header, line, "END", &sl_arr, &nsl_arr);
                    if (nsl >= 0 && nsl_arr == 1) genotype_end.push_back(sl_arr[0]);
                    else genotype_end.push_back(genotype_start.back() + genotype_ref.size() - 1);
                    pair < string, int > temp (sid,genotype_id_to_idx.size());
                    genotype_id_to_idx.insert(temp);
                    n_includedG++;
                } else n_excludedG_user ++;
            } else n_excludedG_void ++;
        } else n_excludedG_mult ++;
    }

    //Finalize
    free(gt_arr);
    free(ds_arr);
    bcf_sr_destroy(sr);
    genotype_count += n_includedG;
    vrb.bullet(stb.str(n_includedG) + " new variants included");
    if (n_excludedG_user > 0) vrb.bullet(stb.str(n_excludedG_user) + " variants excluded by user");
    if (n_excludedG_mult > 0) vrb.bullet(stb.str(n_excludedG_mult) + " multi-allelic variants excluded");
    if (n_excludedG_void > 0) vrb.bullet(stb.str(n_excludedG_void) + " uninformative variants excluded [no GT/DS]");
    if (genotype_count == 0) vrb.leave("Cannot find genotypes in target region!");
}

void union_data::scanGenotypesBED(string fbed) {
	string buffer;
	int n_includedG = 0;
	int n_excludedG_user = 0;

	//Opening files
	htsFile *fp = hts_open(fbed.c_str(),"r");
	if (!fp) vrb.error("Cannot open file!");
	tbx_t * tbx = tbx_index_load(fbed.c_str());
	if (!tbx) vrb.error("Cannot load index file!");
	kstring_t str = {0,0,0};
	if (hts_getline(fp, KS_SEP_LINE, &str) <= 0 || !str.l || str.s[0] != tbx->conf.meta_char ) vrb.error("Cannot read header line!");

	//Read genotype data
	vector < string > tokens;
    unsigned int linecount = 0;
    //Jump to interesting region
    if (regionGenotype.chr != "NA"){
        hts_itr_t *itr = tbx_itr_querys(tbx, regionGenotype.get().c_str());
        vrb.bullet("target region [" + regionGenotype.get() + "]");
        if (!itr) vrb.error("Cannot jump to region!");
        while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
            linecount ++;
            if (linecount % 1000000 == 0) vrb.bullet("Read " + stb.str(linecount) + " lines");
            stb.split(string(str.s), tokens);
            if (tokens.size() < 7) vrb.error("Incorrect number of columns!");
            if (genotype_id_to_idx.count(tokens[3])) continue;
            if (filter_genotype.check(tokens[3])) {
                genotype_id.push_back(tokens[3]);
                genotype_chr.push_back(tokens[0]);
                genotype_start.push_back(atoi(tokens[1].c_str()) + 1);
                genotype_end.push_back(atoi(tokens[2].c_str()));
                pair < string, int > temp (tokens[3],genotype_id_to_idx.size());
                genotype_id_to_idx.insert(temp);
                n_includedG++;
            } else n_excludedG_user ++;
        }
        tbx_itr_destroy(itr);
    }else{
        while (hts_getline(fp, KS_SEP_LINE, &str) >= 0) {
            linecount ++;
            if (linecount % 1000000 == 0) vrb.bullet("Read " + stb.str(linecount) + " lines");
            stb.split(string(str.s), tokens);
            if (str.l && str.s[0] != tbx->conf.meta_char) {
                if (tokens.size() < 7) vrb.error("Incorrect number of columns!");
                if (genotype_id_to_idx.count(tokens[3])) continue;
                if (filter_genotype.check(tokens[3])) {
                    genotype_id.push_back(tokens[3]);
                    genotype_chr.push_back(tokens[0]);
                    genotype_start.push_back(atoi(tokens[1].c_str()) + 1);
                    genotype_end.push_back(atoi(tokens[2].c_str()));
                    pair < string, int > temp (tokens[3],genotype_id_to_idx.size());
                    genotype_id_to_idx.insert(temp);
                    n_includedG++;
                } else n_excludedG_user ++;
            }
        }
    }

	//Finalize & verbose
	tbx_destroy(tbx);
    genotype_count += n_includedG;
	if (hts_close(fp)) vrb.error("Cannot properly close file!");
	vrb.bullet(stb.str(n_includedG) + " new variants included");
	if (n_excludedG_user > 0) vrb.bullet(stb.str(n_excludedG_user) + " variants excluded by user");
    if (n_includedG  == 0) vrb.leave("Cannot find variants in target region!");
}


