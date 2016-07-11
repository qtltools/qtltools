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

#include "trans_data.h"

void trans_data::readSampleFromVCF(string fname) {
	vrb.title("Read sample list from [" + fname + "]");
	bcf_sweep_t * sw = bcf_sweep_init(fname.c_str());
	if (!sw) vrb.error("Cannot open file!");
	bcf_hdr_t * hdr  = bcf_sweep_hdr(sw);
	if (!hdr) vrb.error("Cannot read vcf header!");
	unsigned int n_sample = bcf_hdr_nsamples(hdr);
	for (int i = 0 ; i < n_sample ; i ++) {
		sample_id.push_back(string(hdr->samples[i]));
		sample_count ++;
	}
	vrb.bullet("#samples = " + stb.str(sample_count));
	bcf_sweep_destroy(sw);
}

void trans_data::checkSampleInBED(string fname) {
	unsigned int n_included = 0, n_excluded = 0, n_missing = 0;
	vrb.title("Checking sample list in [" + fname + "]");
	htsFile *fp = hts_open(fname.c_str(),"r");
	if (!fp) vrb.error("Cannot open file");
	kstring_t str = {0,0,0};
	vector < string > tokens;
	if (hts_getline(fp, KS_SEP_LINE, &str) <= 0 || !str.l || str.s[0] != '#' ) vrb.error("Cannot read BED header!");
	stb.split(string(str.s), tokens);
	if (tokens.size() < 7) vrb.error("Incorrect number of columns!");
	vector < bool > sample_found = vector < bool > (sample_count, false);
	for (int t = 6 ; t < tokens.size() ; t ++) {
		int sidx = findSample(tokens[t]);
		if (sidx >= 0) {
			sample_found[sidx] = true;
			n_included ++;
		} else n_excluded ++;
	}
	for (int i = 0 ; i < sample_count ; i ++) if (!sample_found[i]) n_missing ++;
	if (n_missing > 0) {
		vrb.bullet("#samples included = " + stb.str(n_included));
		vrb.bullet("#samples excluded = " + stb.str(n_excluded));
		vrb.bullet("#samples missing = " + stb.str(n_missing));
		vrb.error("Some sample have genotype data and no phenotype data. Trans analysis does allow this for speed purposes.");
	} else if (n_excluded > 0) {
		vrb.bullet("#samples included = " + stb.str(n_included));
		vrb.bullet("#samples excluded = " + stb.str(n_excluded));
	} else vrb.bullet("#samples = " + stb.str(n_included));
	hts_close(fp);
}

void trans_data::checkSampleInCOV(string fname) {
	string buffer; vector < string > tokens;
	unsigned int n_included = 0, n_excluded = 0, n_missing = 0;
	vrb.title("Checking sample list in [" + fname + "]");
	input_file fd (fname);
	if (fd.fail()) vrb.error("Cannot open file!");
	getline(fd, buffer);
	if (buffer.size() == 0) vrb.error("No header line detected!");
	stb.split(buffer, tokens);
	if (tokens.size() < 2) vrb.error("Incorrect number of columns!");
	vector < bool > sample_found = vector < bool > (sample_count, false);
	for (int t = 1 ; t < tokens.size() ; t ++) {
		int sidx = findSample(tokens[t]);
		if (sidx >= 0) {
			sample_found[sidx] = true;
			n_included ++;
		} else n_excluded ++;
	}
	for (int i = 0 ; i < sample_count ; i ++) if (!sample_found[i]) n_missing ++;
	if (n_missing > 0) {
		vrb.bullet("#samples included = " + stb.str(n_included));
		vrb.bullet("#samples excluded = " + stb.str(n_excluded));
		vrb.bullet("#samples missing = " + stb.str(n_missing));
		vrb.error("Some sample have genotype data and no phenotype data!");
	} else if (n_excluded > 0) {
		vrb.bullet("#samples included = " + stb.str(n_included));
		vrb.bullet("#samples excluded = " + stb.str(n_excluded));
	} else vrb.bullet("#samples = " + stb.str(n_included));
	fd.close();
}
