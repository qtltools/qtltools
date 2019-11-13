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

void ase_data::compareChrs(string vcf, string bam, string str_regions){

	vrb.title("Getting chromosomes from BAM [" + bam + "]");
	samFile * fd = sam_open(bam.c_str(), "r");
	if (fd == 0) vrb.error("Failed to open file");
    bam_hdr_t * header = sam_hdr_read(fd);
    if (header == 0) vrb.error("Failed to read header");
    hts_idx_t *idx = sam_index_load(fd, bam.c_str());
    if (idx == NULL) vrb.error("Failed to load index");

    for (int i = 0; i < header->n_targets; i++){
    	bam_chrs.insert(string(header->target_name[i]));
    }
    if (header->n_targets == 0) vrb.error("No chromosomes in BAM header");

    hts_idx_destroy(idx);
    bam_hdr_destroy(header);

	vrb.title("Getting chromosomes from VCF [" + vcf + "]");
	bcf_srs_t * sr =  bcf_sr_init();
	if(!(bcf_sr_add_reader (sr, vcf.c_str()))) {
		switch (sr->errnum) {
		case not_bgzf: vrb.error("Not compressed with bgzip");
		case idx_load_failed: vrb.error("Impossible to load index file");
		case file_type_error: vrb.error("Unrecognized file format");
		default: vrb.error("Unknown error when opening");
		}
	}

	int found_c = 0 , missed_c = 0;
	if (sr->readers[0].header->n[BCF_DT_CTG] == 0) vrb.error("No chromosomes in VCF header");
	for (int  i = 0; i < sr->readers[0].header->n[BCF_DT_CTG]; i++ ){
		string chr(sr->readers[0].header->id[BCF_DT_CTG][i].key);
		vcf_chrs.insert(chr);
		if (bam_chrs.count(chr) == 0){
			if(fix_chr && chr.size() > 3 && chr.substr(0,3) == "chr" && bam_chrs.count(chr.substr(3))){
				found_c++;
				remove_chr.insert(chr);
				vrb.bullet("Removing chr from VCF chromosome [" + chr + "] to match the BAM file");
			}else if (fix_chr && bam_chrs.count("chr" + chr)){
				found_c++;
				add_chr.insert(chr);
				vrb.bullet("Adding chr to VCF chromosome [" + chr + "] to match the BAM file");
			}else{
				missed_c++;
			}
		}else found_c++;
	}
	bcf_sr_destroy(sr);

	if(found_c == 0) vrb.error("No chromosomes match between VCF and BAM. Try --fix-chr!");
	if(missed_c) vrb.warning(stb.str(missed_c) + " VCF chromosomes are missing from the BAM file. Found " + stb.str(found_c) + " chromosomes.");

	if (str_regions.size() > 0) {
		if (!vcf_region.parse(str_regions)) vrb.error("Unable to parse region: " + str_regions);
		if (vcf_chrs.count(vcf_region.chr) == 0){
			if(fix_chr && vcf_region.chr.size() > 3 && vcf_region.chr.substr(0,3) == "chr" && vcf_chrs.count(vcf_region.chr.substr(3))){
				vcf_region.chr = vcf_region.chr.substr(3);
				vrb.warning("Changing region [" + str_regions + "] to [" + vcf_region.get() + "] to match the VCF chromosomes!");
			}else if (fix_chr && vcf_chrs.count("chr" + vcf_region.chr)){
				vcf_region.chr = "chr" + vcf_region.chr;
				vrb.warning("Changing region [" + str_regions + "] to [" + vcf_region.get() + "] to match the VCF chromosomes!");
			}else{
				vrb.error("Chromosome " + vcf_region.chr + " is not in the VCF!");
			}
		}

		bam_region.parse(str_regions);
		if (add_chr.count(bam_region.chr)) bam_region.chr = "chr" + bam_region.chr;
		if (remove_chr.count(bam_region.chr)) bam_region.chr = bam_region.chr.substr(3);

		if (region_length == 0){
			my_regions.push_back(ase_region(bam_region));
			vrb.bullet("Setting BAM region to [" + my_regions.back().get_string() + "]");
		}
	}
}

void ase_data::getRegions(){

	vrb.title("Calculating regions");

	unsigned int pp = 0;
	for (auto it = all_variants.begin(); it != all_variants.end(); it++){
		unsigned int one_based = it->pos+1;
		if (my_regions.size()==0 || my_regions.back().chr != it->chr || one_based - my_regions.back().start + 1 > region_length ){
			ase_region newr(it->chr, one_based);
			if (my_regions.size()) my_regions.back().end = pp;
			my_regions.push_back(newr);
		}else my_regions.back().count++;
		pp = one_based;
	}
	my_regions.back().end = pp;
	vrb.bullet(stb.str(my_regions.size()) + " regions found");
}

void ase_data::collapseRegions(){
	vrb.title("Collapsing regions");
	if (bam_chrs != ase_chrs){
		vrb.bullet("Not all BAM chromosomes have an ASE site, will query individual chromosomes");
		for (auto it = ase_chrs.begin(); it != ase_chrs.end(); it++){
			my_regions.push_back(ase_region(*it,POS_MIN,POS_MAX));
		}
		vrb.bullet(stb.str(my_regions.size()) + " chromosomes will be processed");
	}
}

void ase_data::assignGenesToAseSite(ase_site &in){
	if (annotation.size() && annotation.count(in.chr)){
		unsigned int pos1based = in.pos + 1;
		unsigned int b = pos1based / binsize;
		if (annotation[in.chr].count(b)){
			string anno = "";
			for (int i =0 ; i < annotation[in.chr][b].size(); i++){
				if (annotation[in.chr][b][i].contains(pos1based)){
					if (anno != "") anno += ";";
					anno += annotation[in.chr][b][i].id;
				}
			}
			if (anno != "") in.genes = anno;
		}
	}
}





