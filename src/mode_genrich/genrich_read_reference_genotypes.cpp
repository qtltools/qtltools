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

#include "genrich_data.h"

void genrich_data::readReferenceGenotypes(string fvcf) {
	vector < int > mappingS;

	//Opening files
	vrb.title("Reading variant list in [" + fvcf + "] MAF=" + stb.str(threshold_maf));
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
	int included_sample = 0;
	int n_samples = bcf_hdr_nsamples(sr->readers[0].header);
	for (int i = 0 ; i < n_samples ; i ++) {
		mappingS.push_back(findSample(string(sr->readers[0].header->samples[i])));
		if (mappingS.back() >= 0) included_sample ++;
	}
	vrb.bullet("#samples = " + stb.str(included_sample));

	//Variant processing
	unsigned int n_excludedV_mult = 0, n_excludedV_void = 0, n_excludedV_rare = 0, n_excludedV_uchr = 0, n_line = 0, n_excludedV_toofar = 0;
	int ngt, ngt_arr = 0, *gt_arr = NULL;
	bcf1_t * line;
	while(bcf_sr_next_line (sr)) {
		line =  bcf_sr_get_line(sr, 0);
		if (line->n_allele == 2) {
			bcf_unpack(line, BCF_UN_STR);
			string sid = string(line->d.id);
			string chr = string(bcf_hdr_id2name(sr->readers[0].header, line->rid));
			int chr_idx = findCHR(chr);
			if (chr_idx >= 0) {
				unsigned int pos = line->pos + 1;
				ngt = bcf_get_genotypes(sr->readers[0].header, line, &gt_arr, &ngt_arr);
				if (ngt == 2*n_samples) {
					double freq = 0.0, tot = 0.0;
					for(int i = 0 ; i < n_samples ; i ++) {
						assert(gt_arr[2*i+0] != bcf_gt_missing && gt_arr[2*i+1] != bcf_gt_missing);
						if (mappingS[i] >= 0) {
							freq += bcf_gt_allele(gt_arr[2*i+0]) + bcf_gt_allele(gt_arr[2*i+1]);
							tot += 2.0;
						}
					}
					double maf = freq / tot;
					if (maf > 0.5) maf = 1.0 - maf;
					if (maf >= threshold_maf) {
						int dist_tss = getDistance(chr_idx, pos);
						if (dist_tss < 1e6) {
							string tmp_id = chr + "_" + stb.str(pos);
							genotype_uuid.insert(pair < string, unsigned int > (tmp_id, genotype_pos.size()));
							genotype_chr.push_back(chr_idx);
							genotype_pos.push_back(pos);
							genotype_maf.push_back(maf);
							genotype_dist.push_back(dist_tss);
							genotype_haps.push_back(vector < bool > (2 * included_sample, false));
							for(int i = 0 ; i < n_samples ; i ++) {
								if (mappingS[i] >= 0) {
									genotype_haps.back()[2 * mappingS[i] + 0] = bcf_gt_allele(gt_arr[2 * i + 0]);
									genotype_haps.back()[2 * mappingS[i] + 1] = bcf_gt_allele(gt_arr[2 * i + 1]);
								}
							}
						} else n_excludedV_toofar++;
					} else n_excludedV_rare ++;
				} else n_excludedV_void ++;
			} else n_excludedV_uchr ++;
		} else n_excludedV_mult ++;

		if (n_line % 100000 == 0) vrb.bullet("#lines = " + stb.str(n_line));

		n_line ++;
 	}
	genotype_qtl = vector < bool > (genotype_pos.size(), false);
	genotype_gwas = vector < bool > (genotype_pos.size(), false);
	genotype_bin = vector < int > (genotype_pos.size(), -1);

	//Finalize
	bcf_sr_destroy(sr);
	vrb.bullet(stb.str(genotype_pos.size()) + " variants included");
	if (n_excludedV_mult > 0) vrb.bullet(stb.str(n_excludedV_mult) + " multi-allelic variants excluded");
	if (n_excludedV_uchr > 0) vrb.bullet(stb.str(n_excludedV_uchr) + " variants with unreferenced chromosome in --tss");
	if (n_excludedV_rare > 0) vrb.bullet(stb.str(n_excludedV_rare) + " maf filtered variants");
	if (n_excludedV_toofar > 0) vrb.bullet(stb.str(n_excludedV_toofar) + " too far variants");
}
