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

#ifndef _BAMSTAT_DATA_H
#define _BAMSTAT_DATA_H

//INCLUDES
#include "../common/data.h"

class bamstat_region {
public:
        string chr;
        unsigned int start;
        unsigned int end;
        unsigned int n_covering_reads;

        bamstat_region (string _chr, unsigned int _start, unsigned int _end) {
                chr = _chr;
                start = _start;
                end = _end;
                n_covering_reads = 0;
        }

        string toString() {
                return chr + ":" + stb.str(start) + "-" + stb.str(end);
        }
};

class bamstat_data : public data {
public :

        //DATA
        vector < bamstat_region > R;
        bool param_dup_rd;
        unsigned int param_min_mapQ;
        unsigned int n_total_reads;
        unsigned int n_mapped_reads;
        unsigned int n_overlap_reads;
        unsigned int n_overlap_annotations;
        unsigned int n_keep;

        //CONSTRUCTOR/DESTRUCTOR
        bamstat_data () {
        		param_dup_rd = false;
                param_min_mapQ = 0;
                n_total_reads = 0;
                n_mapped_reads = 0;
                n_overlap_reads = 0;
                n_overlap_annotations = 0;
                n_keep = 0;
        }

        ~bamstat_data () { }

        //
        int keepRead(bam1_t *);
        void readAnnotationsBED(string);
        void readSequences(string);
        void writeOutput(string);
};

//***************************************************************//
//******************** DECLARE FUNCTIONS *************************//
//***************************************************************//
void bamstat_main(vector < string > &);


#endif
