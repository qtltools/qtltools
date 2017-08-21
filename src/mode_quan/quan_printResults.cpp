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


void quan_data::printBEDcount(string fout){
    vrb.title("Printing counts");
    string ext ="";
    string prefix = fout;
    if (fout.substr(fout.find_last_of(".") + 1) == "gz") {
    	ext = ".gz";
    	prefix = fout.substr(0,fout.find_last_of("."));
    }
    if (fout.substr(fout.find_last_of(".") + 1) == "bz2"){
    	ext = ".bz2";
    	prefix = fout.substr(0,fout.find_last_of("."));
    }
	output_file fdo(prefix+".gene.count.bed" + ext);
    output_file fdoe(prefix+".exon.count.bed" + ext);
    if (fdo.fail()) vrb.error("Cannot open file [" + prefix + ".gene.count.bed" + ext + "]");
    if (fdoe.fail()) vrb.error("Cannot open file [" + prefix + ".exon.count.bed" + ext + "]");
    fdo.precision(10);
    fdoe.precision(10);
    fdo << "#chr\tstart\tend\tgene\tinfo\tstrand";
    fdoe << "#chr\tstart\tend\texon\tgeneID\tstrand";
    for (int i = 0 ; i < samples.size(); i++) {fdo << "\t" << samples[i]; fdoe << "\t" << samples[i];}
    fdo<<endl;
    fdoe<<endl;
    for (int gr= 0; gr < gene_grps.size(); gr++){
        for (int g = 0 ; g < gene_grps[gr].genes.size(); g++){
        	string chr = gene_grps[gr].genes[g].chr;
            if (chr.substr(0,3) == "chr") chr.erase(0,3);
            fdo << chr;
            fdo << "\t" << gene_grps[gr].genes[g].tss-1;
            fdo << "\t" << gene_grps[gr].genes[g].tss;
            fdo << "\t" << gene_grps[gr].genes[g].ID;
            fdo << "\tL=" << gene_grps[gr].genes[g].length << ";T=" << gene_grps[gr].genes[g].exons[0].gene_type << ";R=" << gene_grps[gr].genes[g].region << ";N=" << gene_grps[gr].genes[g].exons[0].gene_name;
            fdo << "\t" << (gene_grps[gr].genes[g].strand == -1 ? "-" : "+");
            for (int i = 0 ; i < bams.size(); i++) fdo << "\t" << gene_grps[gr].genes[g].read_count[i];
            fdo << endl;
            for (int e = 0 ; e < gene_grps[gr].genes[g].exons.size(); e++){
            	if (gene_grps[gr].genes[g].exons[e].length < min_exon) continue;
                fdoe << chr;
                fdoe << "\t" << gene_grps[gr].genes[g].tss-1;
                fdoe << "\t" << gene_grps[gr].genes[g].tss;
                fdoe << "\t" << gene_grps[gr].genes[g].exons[e].name;
                fdoe << "\t" << gene_grps[gr].genes[g].ID;
                fdoe << "\t" << (gene_grps[gr].genes[g].strand == -1 ? "-" : "+");
                for (int i = 0 ; i < bams.size(); i++) fdoe << "\t" << gene_grps[gr].genes[g].exons[e].read_count[i];
                fdoe << endl;
            }
        }
    }
}


void quan_data::printBEDrpkm(string fout){
	vrb.title("Printing RPKM");
    string ext ="";
    string prefix = fout;
    if (fout.substr(fout.find_last_of(".") + 1) == "gz") {
    	ext = ".gz";
    	prefix = fout.substr(0,fout.find_last_of("."));
    }
    if (fout.substr(fout.find_last_of(".") + 1) == "bz2"){
    	ext = ".bz2";
    	prefix = fout.substr(0,fout.find_last_of("."));
    }
	output_file fdo(prefix+".gene.rpkm.bed" + ext);
    output_file fdoe(prefix+".exon.rpkm.bed" + ext);
    if (fdo.fail()) vrb.error("Cannot open file [" + prefix + ".gene.rpkm.bed" + ext + "]");
    if (fdoe.fail()) vrb.error("Cannot open file [" + prefix + ".exon.rpkm.bed" + ext + "]");


    fdo.precision(10);
    fdoe.precision(10);
    fdo << "#chr\tstart\tend\tgene\tinfo\tstrand";
    fdoe << "#chr\tstart\tend\texon\tgeneID\tstrand";
    for (int i = 0 ; i < samples.size(); i++) {fdo << "\t" << samples[i]; fdoe << "\t" << samples[i];}
    fdo<<endl;
    fdoe<<endl;
    for (int gr= 0; gr < gene_grps.size(); gr++){
        for (int g = 0 ; g < gene_grps[gr].genes.size(); g++){
            string chr = gene_grps[gr].genes[g].chr;
            if (chr.substr(0,3) == "chr") chr.erase(0,3);
            fdo << chr;
            fdo << "\t" << gene_grps[gr].genes[g].tss-1;
            fdo << "\t" << gene_grps[gr].genes[g].tss;
            fdo << "\t" << gene_grps[gr].genes[g].ID;
            fdo << "\tL=" << gene_grps[gr].genes[g].length << ";T=" << gene_grps[gr].genes[g].exons[0].gene_type << ";R=" << gene_grps[gr].genes[g].region << ";N=" << gene_grps[gr].genes[g].exons[0].gene_name;
            fdo << "\t" << (gene_grps[gr].genes[g].strand == -1 ? "-" : "+");
            for (int i = 0 ; i < bams.size(); i++) fdo << "\t" << ((gene_grps[gr].genes[g].read_count[i] * 1000.0) / (double) gene_grps[gr].genes[g].length) * (1000000.0 / (double)stats[i].exonic);
            fdo << endl;
            for (int e = 0 ; e < gene_grps[gr].genes[g].exons.size(); e++){
            	if (gene_grps[gr].genes[g].exons[e].length < min_exon) continue;
                fdoe << chr;
                fdoe << "\t" << gene_grps[gr].genes[g].tss-1;
                fdoe << "\t" << gene_grps[gr].genes[g].tss;
                fdoe << "\t" << gene_grps[gr].genes[g].exons[e].name;
                fdoe << "\t" << gene_grps[gr].genes[g].ID;
                fdoe << "\t" << (gene_grps[gr].genes[g].strand == -1 ? "-" : "+");
                for (int i = 0 ; i < bams.size(); i++) fdoe << "\t" << ((gene_grps[gr].genes[g].exons[e].read_count[i] * 1000.0) / (double) gene_grps[gr].genes[g].exons[e].length) * (1000000.0 / (double)stats[i].exonic);
                fdoe << endl;
            }
        }
    }
}

void quan_data::printStats(string fout){
	vrb.title("Printing stats");
    string ext ="";
    string prefix = fout;
    if (fout.substr(fout.find_last_of(".") + 1) == "gz") {
    	ext = ".gz";
    	prefix = fout.substr(0,fout.find_last_of("."));
    }
    if (fout.substr(fout.find_last_of(".") + 1) == "bz2"){
    	ext = ".bz2";
    	prefix = fout.substr(0,fout.find_last_of("."));
    }
	output_file fdo(prefix+".stats" + ext);
	if (fdo.fail()) vrb.error("Cannot open file [" + prefix + ".stats" + ext + "]");
	fdo << "sample\tunmmaped_in_genic_regions\tduplicate_reads_in_genic_regions\tfails_mapQ_in_genic_regions:" << min_mapQ <<"\tunpaired_in_genic_regions\tfails_mismatch_in_genic_regions:" << max_mismatch_count << ":" << max_mismatch_count_total << "\tgood_reads_in_genic_regions\tnot_exonic_in_genic_regions\texonic_in_genic_regions\ttotal_reads_in_genic_regions" << endl;
	for (int i = 0 ; i < samples.size(); i++) fdo << samples[i] <<"\t" << stats[i].unmapped << "\t" << stats[i].dup << "\t" << stats[i].mapQ << "\t"<< stats[i].unpaired <<"\t" << stats[i].mismatch << "\t" << stats[i].good << "\t" << stats[i].notexon << "\t" << stats[i].exonicint << "\t" << stats[i].total << endl;
}
