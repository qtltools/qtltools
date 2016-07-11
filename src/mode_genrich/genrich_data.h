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

#ifndef _GENRICH_DATA_H
#define _GENRICH_DATA_H

//ANALYSIS MODES

//INCLUDES
#include "../common/data.h"

class genrich_data : public data {
public:

	//PARAMETERS
	unsigned int n_permutations;
	float threshold_ld;
	float threshold_maf;
	float bin_maf;
	unsigned int bin_distance;

	//DATA FOR CHROMOSOME ID
	vector < string > chromosome_id;
	unordered_map < string, unsigned int > chromosome_idx;

	//DATA FOR PHENOTYPES
	vector <  IntervalTree < pair < bool, bool > > > phenotype_pos;

	//DATA FOR VARIANTS
	vector < unsigned int > genotype_chr;
	vector < int > genotype_pos;
	vector < float > genotype_maf;
	vector < bool > genotype_qtl;
	vector < bool > genotype_gwas;
	vector < int > genotype_dist;
	vector < int > genotype_bin;
	unordered_map < string, unsigned int > genotype_uuid;
	vector < vector < bool > > genotype_haps;

	//QTL BINS
	vector < float > bin_min_maf;
	vector < float > bin_max_maf;
	vector < int > bin_min_dist;
	vector < int > bin_max_dist;

	//CONSTRUCTOR / DESTRUCTOR
	genrich_data() {};
	~genrich_data() {};

	//READ DATA
	void readReferenceGenotypes(string);
	void readQTL(string fqtl);
	void readGWAS(string fgwas);
	void readPhenotypes(string fgwas);

	//PROCESSES
	void binningAllVariants();
	void overlapGWASandQTL(string);

	//ROUTINES
	bool isSameSignal(unsigned int, unsigned int);
	int getDistance(unsigned int, int);
	int findCHR (string &);
};

//***************************************************************//
//******************** DECLARE FUNCTIONS *************************//
//***************************************************************//
void genrich_main(vector < string > &);

#endif
