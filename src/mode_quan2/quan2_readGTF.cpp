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

void quan_data::readGTF(string fgtf, unsigned int nof){
    string buffer;
    vector < string > str;

    vrb.title("Reading exons in [" + fgtf + "]");
    input_file fd (fgtf);
    if (fd.fail()) vrb.error("Cannot open file!");
    int linecount = 0;
    while(getline(fd, buffer)) {
        linecount++;
        if (linecount % 500000 == 0) vrb.bullet(stb.str(linecount) + " lines read");
        if (buffer[0] == '#') continue;
        stb.split(buffer, str, " \t;");
        if (str.size() < 10) vrb.error("Incorrect number of columns: " + stb.str(str.size()));
        if (str.size() % 2 ) vrb.error("Unmatched attributes: " + buffer);
        if (str[2] != "exon") continue;
        string chr = str[0];
        unsigned int start = atoi(str[3].c_str());
        unsigned int end = atoi(str[4].c_str());
        string strand = str[6];
        string gene_type="",gene_id="",gene_name="",trans_type="";
        for (int i = 8 ; i < str.size(); i+=2 ){
            str[i+1].erase(remove(str[i+1].begin(), str[i+1].end(), '"'), str[i+1].end());
            str[i+1].erase(remove(str[i+1].begin(), str[i+1].end(), ';'), str[i+1].end());
            if (str[i] == "gene_name") gene_name = str[i+1];
            if (str[i] == "gene_id") gene_id = str[i+1];
            if (str[i] == "gene_type") gene_type = str[i+1];
            if (str[i] == "transcript_type") trans_type = str[i+1];
        }
        if (gene_id=="") vrb.error("gene_id attribute is required: " + buffer);
        if (gene_type != "" && gene_types.size() && !gene_types.count(gene_type)) continue;
        if (trans_type != "" && gene_types.size() && !gene_types.count(trans_type)) continue;
        //cerr << gene_id << " " << gene_name << " " << chr << " "  << start << " " << end << " " << strand << " " << type << endl;
        quan_exon E(chr, start, end, gene_id, gene_name, gene_type, strand);
        if (!genes_map.count(gene_id)){
            genes_map[gene_id] = genes.size();
            genes.push_back(quan_gene());
            genes.back().assign(E);
        }else genes[genes_map[gene_id]].assign(E);

    }
    groupGenes();
    //for (int i = 0 ; i < genes.size(); i++) cout << genes[i];
}


void quan_data::groupGenes(){
    vrb.title("Sorting and grouping genes");
    sort(genes.begin(),genes.end());
    string pC ="";
    unsigned int pE = 0;
    unsigned int pS = UINT_MAX;
    vector < quan_gene > temp;
    for (int i = 0 ; i < genes.size(); i++){
        if ( pC != "" && (genes[i].chr != pC || genes[i].start - pE > max_read_length)){
            gene_grps.push_back(quan_gene_grp());
            gene_grps.back().genes = temp;
            gene_grps.back().chr = pC;
            gene_grps.back().start = pS;
            gene_grps.back().end = pE;
            gene_grps.back().region = pC + ":" + stb.str(pS) + "-" + stb.str(pE);
            temp.clear();
            temp.push_back(genes[i]);
            pC = genes[i].chr;
            pE = genes[i].end;
            pS = genes[i].start;
        }else{
            pC = genes[i].chr;
            pE = pE < genes[i].end ? genes[i].end : pE;
            pS = pS < genes[i].start ? pS : genes[i].start;
            temp.push_back(genes[i]);
        }
    }
    gene_grps.push_back(quan_gene_grp());
    gene_grps.back().genes = temp;
    gene_grps.back().chr = pC;
    gene_grps.back().start = pS;
    gene_grps.back().end = pE;
    gene_grps.back().region = pC + ":" + stb.str(pS) + "-" + stb.str(pE);
    temp.clear();
    vrb.bullet("Number of genes = " + stb.str(genes.size()));
    vrb.bullet("Number of genes groups = " + stb.str(gene_grps.size()));
    genes.clear();
    genes_map.clear();
}
