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

#ifndef pca_data_h
#define pca_data_h

#define __RESIZE_CHUNK__ 50000

//INCLUDES
#include "../common/data.h"
#include "pca_pca.h"

//dataS
class pca_data : public data {
public:
    double maf_cutoff;
    int distance_separator;
    Pca PCA;
    
    int data_count;								//data number
    MatrixXf data_val;			//data values
    /*vector < string > data_id;						//data ids
    vector < string > data_chr;					//data chromosomes
    vector < int > data_start;						//data start positions
    vector < int > data_end;						//data end positions*/
    
    pca_data(){maf_cutoff = 0.0 ; distance_separator = 0;}
    
    void resizeData();
    void finalizeData(int);
    void imputeData();
    void readData(string);
    void readDataVCF(string);
    void readDataBED(string);
    void readDataPhenoBED(string);
    void printPCA(string);
    
};

//***************************************************************//
//******************** DECLARE FUNCTIONS *************************//
//***************************************************************//
void pca_main(vector < string > &);


#endif /* pca_data_h */
