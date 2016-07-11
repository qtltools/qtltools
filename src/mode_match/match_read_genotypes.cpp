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

void match_data::readGenotypes(string filename, string str_regions) {
	int n_includedG = 0;
	int n_excludedG_mult = 0;
	int n_excludedG_snpv = 0;
	int n_excludedG_void = 0;
	int n_excludedG_user = 0;
	int n_excludedG_impq = 0;
	vector < int > mappingS;

	vrb.title("Reading VCF [" + filename + "]");
	bcf_srs_t * sr =  bcf_sr_init();
	sr->collapse = COLLAPSE_NONE;

	//Jump to regions if necessary
	if (str_regions.size() > 0) {
		if (bcf_sr_set_regions(sr, str_regions.c_str(), 0) == -1) vrb.error("Failed to jump to region [" + str_regions + "]");
		else vrb.bullet("scanning region(s) [" + str_regions + "]");
	} else vrb.bullet("scanning full VCF file");

	//Add readers
	if(!(bcf_sr_add_reader (sr, filename.c_str()))) {
		switch (sr->errnum) {
		case not_bgzf: vrb.error("Not compressed with bgzip");
		case idx_load_failed: vrb.error("Impossible to load index file");
		case file_type_error: vrb.error("Unrecognized file format");
		default: vrb.error("Unknown error when opening");
		}
	}

	//Sample processing
	unsigned int n_samples_in_file = bcf_hdr_nsamples(sr->readers[0].header);
	for (int i = 0 ; i < n_samples_in_file ; i ++) mappingS.push_back(findSample(string(sr->readers[0].header->samples[i])));

	//Init needed data
	int ngp = 0, ngt = 0, niq = 0, ngt_arr = 0, ngp_arr = 0, niq_arr = 0;
	int * gt_arr = NULL;
	float * gp_arr = NULL, * iq_arr = NULL;
	bcf1_t * line;

	//Parse VCF
	map < string, unsigned int > region_map;
	map < string, unsigned int > :: iterator region_map_it;
	while(bcf_sr_next_line (sr)) {
		line =  bcf_sr_get_line(sr, 0);
		if (line->n_allele > 2) n_excludedG_mult ++;
		else {
			bcf_unpack(line, BCF_UN_STR);
			string sid = string(line->d.id);
			if (!filter_genotype.check(sid)) n_excludedG_user ++;
			else {
				string curr_chr = bcf_hdr_id2name(sr->readers[0].header, line->rid);				//chr
				unsigned int pos = line->pos;														//pos
				string ref = string(line->d.allele[0]);												//ref
				string alt = string(line->d.allele[1]);												//alt
				niq = bcf_get_info_float(sr->readers[0].header, line, "IQ", &iq_arr, &niq_arr);		//imp score
				unsigned int region_idx;
				region_map_it = region_map.find(curr_chr);
				if (region_map_it == region_map.end()) {
					vrb.bullet("new chromosome discovered [" + curr_chr + "]");
					region_map.insert(pair < string, unsigned int > (curr_chr, regions.size()));
					region_idx = regions.size();
					regions.push_back(curr_chr);
					sites.push_back(vector < site > ());
					gen_ref.push_back(vector < vector < bool > > ());
					gen_alt.push_back(vector < vector < bool > > ());
				} else region_idx = region_map_it->second;

				if (ref.size() > 1 || alt.size() > 1) n_excludedG_snpv ++;
				else if (niq > 0 && iq_arr[0] < param_min_iq) n_excludedG_impq ++;
				else {
					ngt = bcf_get_genotypes(sr->readers[0].header, line, &gt_arr, &ngt_arr);
					ngp = bcf_get_format_float(sr->readers[0].header, line,"GP", &gp_arr, &ngp_arr);
					if (ngt != n_samples_in_file * 2) n_excludedG_void ++;
					else {
						sites[region_idx].push_back(site (curr_chr, pos, ref, alt));
						gen_ref[region_idx].push_back(vector < bool > (sample_count, false));
						gen_alt[region_idx].push_back(vector < bool > (sample_count, false));
						for (int i = 0 ; i < n_samples_in_file ; i ++) {
							if (mappingS[i] >= 0) {
								bool miss = false;
								if (gt_arr[2*i+0] == bcf_gt_missing || gt_arr[2*i+1] == bcf_gt_missing) miss = true;
								else if (ngp == 3 * n_samples_in_file && gp_arr[3*i+0] != bcf_float_missing && gp_arr[3*i+1] != bcf_float_missing && gp_arr[3*i+2] != bcf_float_missing && gp_arr[3*i+0] < param_min_gp && gp_arr[3*i+1] < param_min_gp && gp_arr[3*i+2] < param_min_gp) miss = true;
								if (!miss) {
									gen_ref[region_idx].back()[mappingS[i]] = !(bcf_gt_allele(gt_arr[2*i+0]) && bcf_gt_allele(gt_arr[2*i+1]));
									gen_alt[region_idx].back()[mappingS[i]] = (bcf_gt_allele(gt_arr[2*i+0]) || bcf_gt_allele(gt_arr[2*i+1]));
								}
							}
						}
						n_includedG ++;
					}
				}
			}
		}
	}
	vrb.bullet(stb.str(n_includedG) + " variants included");
	if (n_excludedG_user > 0) vrb.bullet(stb.str(n_excludedG_user) + " variants excluded by user");
	if (n_excludedG_mult > 0) vrb.bullet(stb.str(n_excludedG_mult) + " multi-allelic variants excluded");
	if (n_excludedG_snpv > 0) vrb.bullet(stb.str(n_excludedG_snpv) + " multi-nucleotidic variants excluded");
	if (n_excludedG_impq > 0) vrb.bullet(stb.str(n_excludedG_user) + " badly imputed variants excluded");
	if (n_excludedG_void > 0) vrb.bullet(stb.str(n_excludedG_void) + " variants without GT field excluded");
	if (sites.size() == 0) vrb.leave("Cannot find usable variants in target region!");
	free(gt_arr);
	bcf_sr_destroy(sr);
}
