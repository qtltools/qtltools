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

void trans_data::setPhenotypeLines(int k, int K) {
    //STEP0: check input values
    if (K < 1) vrb.error("Number of chunks needs to be > 0");
    if (K > phenotype_count) vrb.error("Number of chunks (" + stb.str(K) + ") is greater than the number of phenotypes (" + stb.str(phenotype_count) + ")");
    if (k < 0) vrb.error("Chunk index needs to be > 0");
    if (k > K) vrb.error("Chunk index needs to be smaller than or equal to the total number of chunks [=" + stb.str(K) + "]");
    
    unsigned long int max_length =0 ;
    if (phenotype_count % K == 0) max_length = phenotype_count / K;
    else for ( unsigned long int l = 1 ; l * (K-1) < phenotype_count; l++ ) max_length = l;
    if (K * max_length < phenotype_count){
    	int diff = phenotype_count - (K * max_length);
    	if (k <= diff){
        	start_line = (k-1) * (max_length + 1) + 1;
        	end_line = k * (max_length + 1);
    	}else{
    		int prev = diff * (max_length + 1);
        	start_line = (k-diff-1) * max_length + 1 + prev;
        	end_line = (k-diff) * max_length + prev;
    	}
    }else{
    	start_line = (k-1) * max_length + 1;
    	end_line = k * max_length;
    }
    if (end_line > phenotype_count) end_line = phenotype_count;
}
