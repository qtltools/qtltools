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

#ifndef _FDENSITY_DATA_H
#define _FDENSITY_DATA_H

//INCLUDES
#include "../common/data.h"

class fdensity_data : public data {
public:
	//PARAMETERS
	int window, bin;

	//Annotations
	int ann_count;									//Annotation number
	vector < int > ann_start;						//Annotation start position
	vector < int > ann_end;							//Annotation end position
	vector < string > ann_chr;						//Annotation chromosome

	//TSS
	int tss_count;
	vector < string > tss_id;
	vector < string > tss_chr;
	vector < int > tss_pos;
	vector < bool > tss_neg;

	//TSS functional neighborhood
	vector <  IntervalTree < bool > > Itree;

	//CONSTRUCTOR / DESTRUCTOR
	fdensity_data() {}
	~fdensity_data() {}

	//READ DATA
	void readAnnotation(string);
	void readQTL(string);

	//DATA MANAGEMENT
	void buildIntervalTrees();

	//ANALYSIS
	void runDensityCalculation(string);
};

//***************************************************************//
//******************** DECLARE FUNCTIONS ************************//
//***************************************************************//
void fdensity_main(vector < string > &);

#endif
