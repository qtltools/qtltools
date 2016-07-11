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

#include "ase_data.h"

void ase_data::readGenotypes(string filename, string str_regions) {
	int n_includedG = 0;
	int n_excludedG_mult = 0;
	int n_excludedG_snpv = 0;
	int n_excludedG_void = 0;
	int n_excludedG_user = 0;
	int n_excludedG_impq = 0;
	int n_excludedG_impp = 0;
	int n_excludedG_homo = 0;
	int n_excludedG_miss = 0;

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
	int index_sample = -1;
	unsigned int n_samples_in_file = bcf_hdr_nsamples(sr->readers[0].header);
	for (int i = 0 ; i < n_samples_in_file ; i ++) if (strcmp(sr->readers[0].header->samples[i], sample_id[0].c_str()) == 0) index_sample = i;
	if (index_sample < 0) vrb.error("Unexpected error: sample unfound!");
	else vrb.bullet("index of [" + sample_id[0] + "] = " + stb.str(index_sample));

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
				string sid = string(line->d.id);													//sid
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
					variants.push_back(vector < ase_site > ());
				} else region_idx = region_map_it->second;

				if (ref.size() > 1 || alt.size() > 1) n_excludedG_snpv ++;
				else if (niq > 0 && iq_arr[0] < param_min_iq) n_excludedG_impq ++;
				else {
					ngt = bcf_get_genotypes(sr->readers[0].header, line, &gt_arr, &ngt_arr);
					ngp = bcf_get_format_float(sr->readers[0].header, line,"GP", &gp_arr, &ngp_arr);
					if (ngt != n_samples_in_file * 2) n_excludedG_void ++;
					else if (gt_arr[2*index_sample+0] == bcf_gt_missing || gt_arr[2*index_sample+1] == bcf_gt_missing) n_excludedG_miss ++;
					else if (ngp == 3 * n_samples_in_file && gp_arr[3*index_sample+0] != bcf_float_missing && gp_arr[3*index_sample+1] != bcf_float_missing && gp_arr[3*index_sample+2] != bcf_float_missing && gp_arr[3*index_sample+0] < param_min_gp && gp_arr[3*index_sample+1] < param_min_gp && gp_arr[3*index_sample+2] < param_min_gp) n_excludedG_impp ++;
					else if (bcf_gt_allele(gt_arr[2*index_sample+0]) == bcf_gt_allele(gt_arr[2*index_sample+1])) n_excludedG_homo ++;
					else {
						variants[region_idx].push_back(ase_site (curr_chr, sid, pos, ref, alt));
						n_includedG ++;
					}
				}
			}
		}
	}
	vrb.bullet(stb.str(n_includedG) + " heterozygous genotypes included");
	if (n_excludedG_user > 0) vrb.bullet(stb.str(n_excludedG_user) + " variants excluded by user");
	if (n_excludedG_mult > 0) vrb.bullet(stb.str(n_excludedG_mult) + " multi-allelic variants excluded");
	if (n_excludedG_snpv > 0) vrb.bullet(stb.str(n_excludedG_snpv) + " multi-nucleotidic variants excluded");
	if (n_excludedG_impq > 0) vrb.bullet(stb.str(n_excludedG_user) + " badly imputed variants excluded");
	if (n_excludedG_impp > 0) vrb.bullet(stb.str(n_excludedG_impp) + " badly imputed genotypes excluded");
	if (n_excludedG_void > 0) vrb.bullet(stb.str(n_excludedG_void) + " variants without GT field excluded");
	if (n_excludedG_miss > 0) vrb.bullet(stb.str(n_excludedG_miss) + " missing genotypes excluded");
	if (n_excludedG_homo > 0) vrb.bullet(stb.str(n_excludedG_homo) + " homozygous genotypes excluded");
	if (variants.size() == 0) vrb.leave("Cannot find usable variants in target region!");
	free(gt_arr);
	bcf_sr_destroy(sr);
}
