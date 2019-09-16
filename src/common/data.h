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

#ifndef _DATA_H
#define _DATA_H

#include "otools.h"
#include "filter.h"

#define QTLTOOLS_VERSION "1.2"

class data {
public:
	//SAMPLES
	unsigned sample_count;
	vector < string > sample_id;
	unsigned file_count;
	map < string, unsigned int > sample_occurrence;

	//FILTERS
	filter filter_sample;
	filter filter_phenotype;
	filter filter_genotype;
	filter filter_position;
	filter filter_covariate;

	//OPTIONS
	boost::program_options::options_description option_descriptions;
	boost::program_options::variables_map options;

	//CONSTRUCTOR
	data() { file_count = 0; sample_count = 0; };
	~data() { sample_id.clear(); sample_occurrence.clear(); };

	//SAMPLE NAMES
	void readSampleFromTXT(string);
	void readSampleFromSTR(string);
	void readSampleFromVCF(string, bool silent = false);
	void readSampleFromBED(string, bool silent = false);
	void readSampleFromCOV(string, bool silent = false);
	void mergeSampleLists(bool silent = false);
	int findSample(string);

	//COMMON OPTIONS
	void declareBasicOptions();
	void processBasicOptions();

	//INCLUSION/EXCLUSION LISTS
	void readSampleExclusion(string);
	void readSampleInclusion(string);
	void readGenotypeExclusion(string);
	void readGenotypeInclusion(string);
	void readPositionExclusion(string);
	void readPositionInclusion(string);
	void readPhenotypeExclusion(string);
	void readPhenotypeInclusion(string);
	void readCovariateExclusion(string);
	void readCovariateInclusion(string);
};

#endif
