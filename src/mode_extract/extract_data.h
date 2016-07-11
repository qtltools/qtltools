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

#ifndef _EXTRACT_DATA_H
#define _EXTRACT_DATA_H

//INCLUDES
#include "../common/data.h"

class extract_data : public data {
public:

	//REGIONS
	genomic_region regionData;

	//DATA
	vector < string > variable_id;
	vector < string > variable_chr;
	vector < int > variable_start;
	vector < int > variable_end;
	vector < vector < string > > variable_val;

	//CONSTRUCTOR / DESTRUCTOR
	extract_data() {}
	~extract_data() {}
	void clear() { variable_val.clear(); }

	//READ & WRITE DATA
	void readBED(string);
	void readVCF(string);
	void readCOV(string);
	void writeOUT(string);

	//DATA MANAGMENT
	void imputeMissing();

};

//***************************************************************//
//******************** DECLARE FUNCTIONS *************************//
//***************************************************************//
void extract_main(vector < string > &);

#endif
