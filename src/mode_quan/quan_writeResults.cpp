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

void quan2_data::printBEDcount(string fout){
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
    string hash_to_add = !options.count("no-hash") ? "." + hash : "";

	output_file fdo(prefix+hash_to_add+".gene.count.bed" + ext);
    output_file fdoe(prefix+hash_to_add+".exon.count.bed" + ext);
    if (fdo.fail()) vrb.error("Cannot open file [" + prefix + ".gene.count.bed" + ext + "]");
    if (fdoe.fail()) vrb.error("Cannot open file [" + prefix + ".exon.count.bed" + ext + "]");
    fdo.precision(10);
    fdoe.precision(10);
    fdo << "#chr\tstart\tend\tgene\tinfo\tstrand\t" << sample << endl;
    fdoe << "#chr\tstart\tend\texon\tgeneID\tstrand\t" << sample << endl;
    for (int g = 0 ; g < genes.size(); g++){
    	if(region.isSet() && !genes[g].overlap(region)) continue;
    	string chr = genes[g].chr;
    	if (chr.substr(0,3) == "chr") chr.erase(0,3);
    	fdo << chr;
    	fdo << "\t" << genes[g].tss-1;
    	fdo << "\t" << genes[g].tss;
    	fdo << "\t" << genes[g].gene_id;
    	fdo << "\tL=" << genes[g].length << ";T=" << genes[g].exons[0].gene_type << ";R=" << genes[g].region << ";N=" << genes[g].exons[0].gene_name;
    	fdo << "\t" << (genes[g].strand == -1 ? "-" : "+");
    	fdo << "\t" << genes[g].read_count;
    	fdo << endl;
    	for (int e = 0 ; e < genes[g].exons.size(); e++){
    		if (genes[g].exons[e].length < filter.min_exon) continue;
    		fdoe << chr;
    		fdoe << "\t" << genes[g].tss-1;
    		fdoe << "\t" << genes[g].tss;
    		fdoe << "\t" << genes[g].exons[e].name;
    		fdoe << "\t" << genes[g].gene_id;
    		fdoe << "\t" << (genes[g].strand == -1 ? "-" : "+");
    		fdoe << "\t" << genes[g].exons[e].read_count;
    		fdoe << endl;
    	}
    }
}


void quan2_data::printBEDrpkm(string fout){
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
    string hash_to_add = !options.count("no-hash") ? "." + hash : "";
	output_file fdo(prefix+hash_to_add+".gene.rpkm.bed" + ext);
    output_file fdoe(prefix+hash_to_add+".exon.rpkm.bed" + ext);
    if (fdo.fail()) vrb.error("Cannot open file [" + prefix + ".gene.rpkm.bed" + ext + "]");
    if (fdoe.fail()) vrb.error("Cannot open file [" + prefix + ".exon.rpkm.bed" + ext + "]");


    fdo.precision(10);
    fdoe.precision(10);
    fdo << "#chr\tstart\tend\tgene\tinfo\tstrand\t" << sample << endl;
    fdoe << "#chr\tstart\tend\texon\tgeneID\tstrand\t" << sample << endl;

    for (int g = 0 ; g < genes.size(); g++){
    	if(region.isSet() && !genes[g].overlap(region)) continue;
    	string chr = genes[g].chr;
    	if (chr.substr(0,3) == "chr") chr.erase(0,3);
    	fdo << chr;
    	fdo << "\t" << genes[g].tss-1;
    	fdo << "\t" << genes[g].tss;
    	fdo << "\t" << genes[g].gene_id;
    	fdo << "\tL=" << genes[g].length << ";T=" << genes[g].exons[0].gene_type << ";R=" << genes[g].region << ";N=" << genes[g].exons[0].gene_name;
    	fdo << "\t" << (genes[g].strand == -1 ? "-" : "+");
    	fdo << "\t" << ((genes[g].read_count * 1000.0) / (double) genes[g].length) * (1000000.0 / (double) stats.exonic_multi);
    	fdo << endl;
    	for (int e = 0 ; e < genes[g].exons.size(); e++){
    		if (genes[g].exons[e].length < filter.min_exon) continue;
    		fdoe << chr;
    		fdoe << "\t" << genes[g].tss-1;
    		fdoe << "\t" << genes[g].tss;
    		fdoe << "\t" << genes[g].exons[e].name;
    		fdoe << "\t" << genes[g].gene_id;
    		fdoe << "\t" << (genes[g].strand == -1 ? "-" : "+");
    		fdoe << "\t" << ((genes[g].exons[e].read_count * 1000.0) / (double) genes[g].exons[e].length) * (1000000.0 / (double)stats.exonic_multi);
    		fdoe << endl;
    	}
    }
}

void quan2_data::printBEDtpm(string fout){
	vrb.title("Printing TPM");
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
    string hash_to_add = !options.count("no-hash") ? "." + hash : "";
	output_file fdo(prefix+hash_to_add+".gene.tpm.bed" + ext);
    output_file fdoe(prefix+hash_to_add+".exon.tpm.bed" + ext);
    if (fdo.fail()) vrb.error("Cannot open file [" + prefix + ".gene.tpm.bed" + ext + "]");
    if (fdoe.fail()) vrb.error("Cannot open file [" + prefix + ".exon.tpm.bed" + ext + "]");


    fdo.precision(10);
    fdoe.precision(10);
    fdo << "#chr\tstart\tend\tgene\tinfo\tstrand\t" << sample << endl;
    fdoe << "#chr\tstart\tend\texon\tgeneID\tstrand\t" << sample << endl;
    double rpk_gene = 0.0, rpk_exon = 0.0;
    for (int g = 0 ; g < genes.size(); g++){
    	rpk_gene += (genes[g].read_count * 1000.0) / (double) genes[g].length;
    	for (int e = 0 ; e < genes[g].exons.size(); e++){
    		rpk_exon += (genes[g].exons[e].read_count * 1000.0) / (double) genes[g].exons[e].length;
    	}
    }
    for (int g = 0 ; g < genes.size(); g++){
    	if(region.isSet() && !genes[g].overlap(region)) continue;
    	string chr = genes[g].chr;
    	if (chr.substr(0,3) == "chr") chr.erase(0,3);
    	fdo << chr;
    	fdo << "\t" << genes[g].tss-1;
    	fdo << "\t" << genes[g].tss;
    	fdo << "\t" << genes[g].gene_id;
    	fdo << "\tL=" << genes[g].length << ";T=" << genes[g].exons[0].gene_type << ";R=" << genes[g].region << ";N=" << genes[g].exons[0].gene_name;
    	fdo << "\t" << (genes[g].strand == -1 ? "-" : "+");
    	fdo << "\t" << ((genes[g].read_count * 1000.0) / (double) genes[g].length) * (1000000.0 / rpk_gene);
    	fdo << endl;
    	for (int e = 0 ; e < genes[g].exons.size(); e++){
    		if (genes[g].exons[e].length < filter.min_exon) continue;
    		fdoe << chr;
    		fdoe << "\t" << genes[g].tss-1;
    		fdoe << "\t" << genes[g].tss;
    		fdoe << "\t" << genes[g].exons[e].name;
    		fdoe << "\t" << genes[g].gene_id;
    		fdoe << "\t" << (genes[g].strand == -1 ? "-" : "+");
    		fdoe << "\t" << ((genes[g].exons[e].read_count * 1000.0) / (double) genes[g].exons[e].length) * (1000000.0 / rpk_exon);
    		fdoe << endl;
    	}
    }
}

void quan2_data::printStats(string fout){
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
    string hash_to_add = !options.count("no-hash") ? "." + hash : "";
	output_file fdo(prefix+hash_to_add+".stats" + ext);
	if (fdo.fail()) vrb.error("Cannot open file [" + prefix + ".stats" + ext + "]");
	fdo.precision(16);
	fdo << "sample\t" << sample << endl;
	fdo << "total_reads\t" << stats.total << endl;
	fdo << "total_secondary_alingments\t" << stats.secondary << endl;
	fdo << "total_unmapped\t" << stats.unmapped << endl;
	fdo << "total_failqc\t" << stats.failqc << endl;
	fdo << "total_duplicate\t" << stats.dup << endl;
	fdo << "total_mapQ_less_than_" << filter.min_mapQ << "\t" << stats.mapQ << endl;
	fdo << "total_notpaired\t" << stats.unpaired << endl;
	fdo << "total_mismatches_greater_than_" << filter.max_mismatch_count << "_" << filter.max_mismatch_count_total << "\t" << stats.mismatch << endl;
	fdo << "total_merged_reads\t" << stats.merged << endl;
	fdo << "total_good\t" << stats.good << endl;
	fdo << "total_exonic\t" << stats.exonicint << endl;
	//fdo << "total_exonic_after_merge\t" << stats.exonic << endl;
	fdo << "total_exonic_multi_counting\t" << stats.exonicint_multi << endl;
	fdo << "total_exonic_multi_counting_after_merge_(used_for_rpkm)\t" << stats.exonic_multi << endl;
	fdo << "good_over_total\t" << (double)stats.good / (double)stats.total << endl;
	fdo << "exonic_over_total\t" << (double)stats.exonicint / (double)stats.total << endl;
	fdo << "exonic_over_good\t" << (double)stats.exonicint / (double)stats.good << endl;
}


string quan2_data::convertToBase(unsigned long long num , unsigned base ){
    char digits[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_-+=;~";
    int i;
    char buf[514]; //enough to hold a 512 bit in base 2
    
    if (base < 2 || base > 68) return "";
    if (!num) return "0";
    
    buf[513] = '\0';
    i = 513;
    
    while (num) {
        buf[--i] = digits[num % base];
        num /= base;
    }
    return string(strdup(buf + i));
}
