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

#include "data.h"

//SAMPLE NAMES
void data::readSampleFromVCF(string fname, bool silent) {
	unsigned int n_included = 0, n_excluded = 0;
	if (!silent) vrb.title("Read sample list from [" + fname + "]");
	bcf_sweep_t * sw = bcf_sweep_init(fname.c_str());
	if (!sw) vrb.error("Cannot open file!");
	bcf_hdr_t * hdr  = bcf_sweep_hdr(sw);
	if (!hdr) vrb.error("Cannot read vcf header!");
	unsigned int n_sample = bcf_hdr_nsamples(hdr);
	for (int i = 0 ; i < n_sample ; i ++) {
		string sid = string(hdr->samples[i]);
		if (filter_sample.check(sid)) {
			map < string, unsigned int > :: iterator it_SO = sample_occurrence.find(sid);
			if (it_SO == sample_occurrence.end()) sample_occurrence.insert(make_pair(sid, 1));
			else it_SO->second ++;
			n_included ++;
		} else n_excluded ++;
	}
	if (!silent){
		if (n_excluded == 0) vrb.bullet("#samples = " + stb.str(n_included));
		else {
			vrb.bullet("#samples included by user = " + stb.str(n_included));
			vrb.bullet("#samples excluded by user = " + stb.str(n_excluded));
		}
	}
	file_count ++;
	bcf_sweep_destroy(sw);
}

void data::readSampleFromBED(string fname, bool silent) {
	unsigned int n_included = 0, n_excluded = 0;
	if (!silent) vrb.title("Reading sample list from [" + fname + "]");
	htsFile *fp = hts_open(fname.c_str(),"r");
	if (!fp) vrb.error("Cannot open file");
	kstring_t str = {0,0,0};
	vector < string > tokens;
	if (hts_getline(fp, KS_SEP_LINE, &str) <= 0 || !str.l || str.s[0] != '#' ) vrb.error("Cannot read BED header!");
	stb.split(string(str.s), tokens);
	if (tokens.size() < 7) vrb.error("Incorrect number of columns!");
	for (int t = 6 ; t < tokens.size() ; t ++) {
		string sid = tokens[t];
		if (filter_sample.check(sid)) {
			map < string, unsigned int > :: iterator it_SO = sample_occurrence.find(sid);
			if (it_SO == sample_occurrence.end()) sample_occurrence.insert(make_pair(sid, 1));
			else it_SO->second ++;
			n_included ++;
		} else n_excluded ++;
	}
	if (!silent){
		if (n_excluded == 0) vrb.bullet("#samples = " + stb.str(n_included));
		else {
			vrb.bullet("#samples included by user = " + stb.str(n_included));
			vrb.bullet("#samples excluded by user = " + stb.str(n_excluded));
		}
	}
	file_count ++;
	hts_close(fp);
}

void data::readSampleFromCOV(string fname, bool silent) {
	string buffer; vector < string > tokens;
	unsigned int n_included = 0, n_excluded = 0;
	if (!silent) vrb.title("Reading sample list from [" + fname + "]");
	input_file fd (fname);
	if (fd.fail()) vrb.error("Cannot open file!");
	getline(fd, buffer);
	if (buffer.size() == 0) vrb.error("No header line detected!");
	stb.split(buffer, tokens);
	if (tokens.size() < 2) vrb.error("Incorrect number of columns!");
	for (int t = 1 ; t < tokens.size() ; t ++) {
		string sid = tokens[t];
		if (filter_sample.check(sid)) {
			map < string, unsigned int > :: iterator it_SO = sample_occurrence.find(sid);
			if (it_SO == sample_occurrence.end()) sample_occurrence.insert(make_pair(sid, 1));
			else it_SO->second ++;
			n_included ++;
		} else n_excluded ++;
	}
	if (!silent){
		if (n_excluded == 0) vrb.bullet("#samples = " + stb.str(n_included));
		else {
			vrb.bullet("#samples included by user = " + stb.str(n_included));
			vrb.bullet("#samples excluded by user = " + stb.str(n_excluded));
		}
	}
	file_count ++;
	fd.close();
}

void data::readSampleFromTXT(string fname) {
	string buffer; vector < string > tokens;
	unsigned int n_included = 0, n_excluded = 0;
	vrb.title("Reading sample list from [" + fname + "]");
	input_file fd (fname);
	if (fd.fail()) vrb.error("Cannot open file!");
	while(getline(fd, buffer)) {
		stb.split(buffer, tokens);
		string sid = tokens[0];
		if (filter_sample.check(sid)) {
			map < string, unsigned int > :: iterator it_SO = sample_occurrence.find(sid);
			if (it_SO == sample_occurrence.end()) sample_occurrence.insert(make_pair(sid, 1));
			else it_SO->second ++;
			n_included ++;
		} else n_excluded ++;
	}
	if (n_excluded == 0) vrb.bullet("#samples = " + stb.str(n_included));
	else {
		vrb.bullet("#samples included by user = " + stb.str(n_included));
		vrb.bullet("#samples excluded by user = " + stb.str(n_excluded));
	}
	file_count ++;
	fd.close();
}

void data::readSampleFromSTR(string fname) {
	string buffer; vector < string > tokens;
	unsigned int n_included = 0, n_excluded = 0;
	vrb.title("Parse sample list [" + fname + "]");
	stb.split(fname, tokens, ",");
	for (int t = 0 ; t < tokens.size() ; t++) if (filter_sample.check(tokens[t])) {
		map < string, unsigned int > :: iterator it_SO = sample_occurrence.find(tokens[t]);
		if (it_SO == sample_occurrence.end()) sample_occurrence.insert(make_pair(tokens[t], 1));
		else it_SO->second ++;
		n_included ++;
	} else n_excluded ++;
	if (n_excluded == 0) vrb.bullet("#samples = " + stb.str(n_included));
	else {
		vrb.bullet("#samples included by user = " + stb.str(n_included));
		vrb.bullet("#samples excluded by user = " + stb.str(n_excluded));
	}
	file_count ++;
}


void data::mergeSampleLists(bool silent) {
	unsigned int n_incomplete = 0;
	if (!silent){
		vrb.title("Merging sample sets from all input files");
		vrb.bullet("#files = " + stb.str(file_count));
	}
	for (map < string, unsigned int > :: iterator it_SO = sample_occurrence.begin() ; it_SO != sample_occurrence.end() ; ++it_SO) {
		if (it_SO->second == file_count) {
			sample_id.push_back(it_SO->first);
			sample_count ++;
		} else n_incomplete ++;
	}
	if (!silent){
		vrb.bullet("#samples in all files = " + stb.str(sample_count));
		if (n_incomplete > 0) vrb.bullet("#samples NOT in all files (i.e. ignored) = " + stb.str(n_incomplete));
	}
}

int data::findSample(string sid) {
	for (unsigned int i = 0 ; i < sample_count ; i ++) if (sample_id[i] == sid) return i;
	return -1;
}

void data::declareBasicOptions() {
	boost::program_options::options_description opt_basic ("\x1B[32mBasics\33[0m");
	opt_basic.add_options()
		("help", "Produces option description")
		("seed", boost::program_options::value< unsigned int >()->default_value(15112011), "Random number seed. Useful to replicate runs.")
		("log", boost::program_options::value< string >(), "Output on screen goes to this file.")
		("silent", "Disable screen output");

	boost::program_options::options_description opt_exclusion ("\x1B[32mData Exclusion/Inclusion\33[0m");
	opt_exclusion.add_options()
		("exclude-samples", boost::program_options::value< string >(), "List of samples to exclude.")
		("include-samples", boost::program_options::value< string >(), "List of samples to include.")
		("exclude-sites", boost::program_options::value< string >(), "List of sites to exclude.")
		("include-sites", boost::program_options::value< string >(), "List of sites to include.")
		("exclude-positions", boost::program_options::value< string >(), "List of positions to exclude.")
		("include-positions", boost::program_options::value< string >(), "List of positions to include.")
		("exclude-phenotypes", boost::program_options::value< string >(), "List of phenotypes to exclude.")
		("include-phenotypes", boost::program_options::value< string >(), "List of phenotypes to include.")
		("exclude-covariates", boost::program_options::value< string >(), "List of covariates to exclude.")
		("include-covariates", boost::program_options::value< string >(), "List of covariates to include.");

	option_descriptions.add(opt_basic).add(opt_exclusion);
}

void data::processBasicOptions() {
	if (options.count("exclude-samples")) readSampleExclusion(options["exclude-samples"].as < string > ());
	if (options.count("include-samples")) readSampleInclusion(options["include-samples"].as < string > ());
	if (options.count("exclude-sites")) readGenotypeExclusion(options["exclude-sites"].as < string > ());
	if (options.count("include-sites")) readGenotypeInclusion(options["include-sites"].as < string > ());
	if (options.count("exclude-positions")) readPositionExclusion(options["exclude-positions"].as < string > ());
	if (options.count("include-positions")) readPositionInclusion(options["include-positions"].as < string > ());
	if (options.count("exclude-phenotypes")) readPhenotypeExclusion(options["exclude-phenotypes"].as < string > ());
	if (options.count("include-phenotypes")) readPhenotypeInclusion(options["include-phenotypes"].as < string > ());
	if (options.count("exclude-covariates")) readCovariateExclusion(options["exclude-covariates"].as < string > ());
	if (options.count("include-covariates")) readCovariateInclusion(options["include-covariates"].as < string > ());

	vrb.title("Initialize random number generator");
	if (!options["seed"].defaulted()) vrb.bullet("User specified seed is " + stb.str(options["seed"].as < unsigned int > ()));
	else vrb.bullet("Built-in seed is 15112011");
	rng.setSeed(options["seed"].as < unsigned int > ());
	vrb.bullet("First Integer = " + stb.str(rng.getInt(32768)));
	vrb.bullet("First Double = " + stb.str(rng.getDouble()));
}

void data::readSampleExclusion(string fname){
	vrb.title("Read sample exclusion list [" + fname + "]");
	int ret = filter_sample.readExclusion(fname);
	if (ret < 0) vrb.error("Cannot open file!");
	else vrb.bullet(stb.str(ret) + " samples");
}

void data::readSampleInclusion(string fname){
	vrb.title("Read sample inclusion list [" + fname + "]");
	int ret = filter_sample.readInclusion(fname);
	if (ret < 0) vrb.error("Cannot open file!");
	else vrb.bullet(stb.str(ret) + " samples");
}

void data::readGenotypeExclusion(string fname){
	vrb.title("Read variant exclusion list [" + fname + "]");
	int ret = filter_genotype.readExclusion(fname);
	if (ret < 0) vrb.error("Cannot open file!");
	else vrb.bullet(stb.str(ret) + " variants");
}

void data::readGenotypeInclusion(string fname){
	vrb.title("Read variant inclusion list [" + fname + "]");
	int ret = filter_genotype.readInclusion(fname);
	if (ret < 0) vrb.error("Cannot open file!");
	else vrb.bullet(stb.str(ret) + " variants");
}

void data::readPositionExclusion(string fname){
	vrb.title("Read position exclusion list [" + fname + "]");
	int ret = filter_position.readExclusion(fname, true);
	if (ret < 0) vrb.error("Cannot open file!");
	else vrb.bullet(stb.str(ret) + " positions");
}

void data::readPositionInclusion(string fname){
	vrb.title("Read position inclusion list [" + fname + "]");
	int ret = filter_position.readInclusion(fname, true);
	if (ret < 0) vrb.error("Cannot open file!");
	else vrb.bullet(stb.str(ret) + " positions");
}

void data::readPhenotypeExclusion(string fname){
	vrb.title("Read phenotype exclusion list [" + fname + "]");
	int ret = filter_phenotype.readExclusion(fname);
	if (ret < 0) vrb.error("Cannot open file!");
	else vrb.bullet(stb.str(ret) + " phenotypes");
}

void data::readPhenotypeInclusion(string fname){
	vrb.title("Read phenotype inclusion list [" + fname + "]");
	int ret = filter_phenotype.readInclusion(fname);
	if (ret < 0) vrb.error("Cannot open file!");
	else vrb.bullet(stb.str(ret) + " phenotypes");
}

void data::readCovariateExclusion(string fname){
	vrb.title("Read covariate exclusion list [" + fname + "]");
	int ret = filter_covariate.readExclusion(fname);
	if (ret < 0) vrb.error("Cannot open file!");
	else vrb.bullet(stb.str(ret) + " covariates");
}

void data::readCovariateInclusion(string fname){
	vrb.title("Read covariate inclusion list [" + fname + "]");
	int ret = filter_covariate.readInclusion(fname);
	if (ret < 0) vrb.error("Cannot open file!");
	else vrb.bullet(stb.str(ret) + " covariates");
}
