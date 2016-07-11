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

#ifndef _FENRICH_DATA_H
#define _FENRICH_DATA_H

//INCLUDES
#include "../common/data.h"

class fenrich_data : public data {
public:
	//PARAMETERS
	unsigned int n_permutation;

	//QTL
	int qtl_count;							//QTL number
	vector < int > qtl_pos;					//QTL variant start position
	vector < int > qtl_order;

	//Annotations
	int ann_count;									//Annotation number
	vector < int > ann_start;						//Annotation start position
	vector < int > ann_end;							//Annotation end position
	vector < string > ann_chr;						//Annotation chromosome

	//Tss
	int tss_count;
	vector < string > tss_id;
	vector < string > tss_chr;
	vector < int > tss_pos;
	vector < bool > tss_neg;

	//QTL functional neighborhood
	vector <  IntervalTree < bool > > R;

	//CONSTRUCTOR / DESTRUCTOR
	fenrich_data() {}
	~fenrich_data() {}

	//READ DATA
	void readQTL(string);
	void readAnnotation(string);
	void readTSS(string);

	//DATA MANAGEMENT
	int findTSS(string &);
	void mapAnnotation2QTL();

	//COMPUTATIONNAL ROUTINES
	unsigned int countOverlaps();

	//ANALYSIS
	void runEnrichmentPass(string);
};

//***************************************************************//
//******************** DECLARE FUNCTIONS ************************//
//***************************************************************//
void fenrich_main(vector < string > &);

#endif
