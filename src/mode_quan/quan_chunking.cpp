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

#include "quan_data.h"


void quan_data::setRegion(string r){
    if (!region.parse(r)) vrb.error("Unable to parse [" + r +"]");
    vector < quan_gene_grp > new_genes_grp;
    for (int g = 0 ;g < gene_grps.size(); g++) if(gene_grps[g].overlap(region)) new_genes_grp.push_back(gene_grps[g]);
    gene_grps = new_genes_grp;
    vrb.bullet("Number of gene groups in [" + region.get() +"] = " + stb.str(gene_grps.size()));
}

void quan_data::setChunk(int k, int K){
    //STEP0: check input values
    if (K < 1) vrb.error("Number of chunks needs to be > 0");
    if (K > gene_grps.size()) vrb.error("Number of chunks (" + stb.str(K) + ") is greater than the number of gene groups (" + stb.str(gene_grps.size()) + ")");
    if (k < 0) vrb.error("Chunk index needs to be > 0");
    if (k > K) vrb.error("Chunk index needs to be smaller than or equal to the total number of chunks [=" + stb.str(K) + "]");
    
    unsigned long int max_length =0 ;
    if (gene_grps.size() % K == 0) max_length = gene_grps.size() / K;
    else for ( unsigned long int l = 1 ; l * (K-1) < gene_grps.size(); l++ ) max_length = l;
    unsigned long int start_idx = (k-1) * max_length;
    unsigned long int end_idx = k * max_length;
    if (end_idx > gene_grps.size()) end_idx = gene_grps.size();
    gene_grps = vector < quan_gene_grp > (gene_grps.begin()+start_idx, gene_grps.begin()+end_idx);
    vrb.bullet("Number of gene groups in chunk [" + stb.str(k) + " / " + stb.str(K) +"] = " + stb.str(gene_grps.size()));
}