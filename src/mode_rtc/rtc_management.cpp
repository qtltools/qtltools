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

#include "rtc_data.h"

void rtc_data::imputeGenotypes() {
	vrb.title("Imputing missing genotypes");
	for (int g = 0; g < genotype_count ; g ++) {
		double mean = 0.0;
		int c_mean = 0;
		for (int s = 0; s < sample_count ; s ++) {
			if (genotype_val[g][s] != bcf_float_missing && !std::isnan(genotype_val[g][s])) {
				mean += genotype_val[g][s];
				c_mean ++;
			}
		}
		mean /= c_mean;
		for (int s = 0; s < sample_count ; s ++) if (genotype_val[g][s] == bcf_float_missing || std::isnan(genotype_val[g][s])) genotype_val[g][s] = mean;
	}
}

void rtc_data::imputePhenotypes() {
	vrb.title("Imputing missing phenotypes");
	for (int p = 0; p < phenotype_count ; p ++) {
		double mean = 0.0;
		int c_mean= 0;
		for (int s = 0; s < sample_count; s ++) {
			if (phenotype_val[p][s] != bcf_float_missing && !std::isnan(phenotype_val[p][s])) {
				mean += phenotype_val [p][s];
				c_mean ++;
			}
		}
		mean /= c_mean;
		for (int s = 0; s < sample_count ; s ++) if (phenotype_val[p][s] == bcf_float_missing || std::isnan(phenotype_val[p][s])) phenotype_val[p][s] = mean;
	}
}

void rtc_data::normalTransform(vector < float > & V) {
	vector < float > R;
	myranker::rank(V, R);
	double max = 0;
	for (int s = 0 ; s < sample_count ; s ++) {
		R[s] = R[s] - 0.5;
		if (R[s] > max) max = R[s];
	}
	max = max + 0.5;
	for (int s = 0 ; s < sample_count ; s ++) {
		R[s] /= max;
		V[s] = qnorm(R[s], 0.0, 1.0, 1, 0);
	}
}

void rtc_data::normalTransformPhenotypes() {
	vrb.title("Match phenotypes to Normal distribution");
	for (int p = 0; p < phenotype_count ; p ++) normalTransform(phenotype_val[p]);
}

void rtc_data::normalize(vector < float > & X) {
	double mean = 0.0, sum = 0.0;
	for (int s = 0; s < X.size() ; s ++) mean += X[s];
	mean /= X.size();
	for (int s = 0; s < X.size() ; s ++) {
		X[s] -= mean;
		sum += X[s] * X[s];
	}
	sum = sqrt(sum);
	if (sum == 0) sum = 1;
	for (int s = 0; s < X.size() ; s ++) X[s] /= sum;
}


void rtc_data::normalize(vector < vector < float > > & X) {
	for (int x = 0 ; x < X.size() ; x++) {
		double mean = 0.0, sum = 0.0;
		for (int s = 0; s < sample_count ; s ++) mean += X[x][s];
		mean /= sample_count;
		for (int s = 0; s < sample_count ; s ++) {
			X[x][s] -= mean;
			sum += X[x][s] * X[x][s];
		}
		sum = sqrt(sum);
		if (sum == 0) sum = 1;
		for (int s = 0; s < sample_count ; s ++) X[x][s] /= sum;
	}
}

vector <float> rtc_data::correct(vector < float >  X, vector < float >  Y) {
    vector < float > R(sample_count);
    double corr = getCorrelation(X, Y);
    for (int s = 0 ; s < sample_count ; s ++) R[s] = Y[s] - corr * X[s];
    return R;
}

// true if all samples are phased.
// haploid genotypes are considered phased
// ./. => not phased, .|. => phased
// From vcfview.c
int rtc_data::bcf_all_phased(const bcf_hdr_t *header, bcf1_t *line){
    bcf_unpack(line, BCF_UN_FMT);
    bcf_fmt_t *fmt_ptr = bcf_get_fmt(header, line, "GT");
    int all_phased = 1;
    if ( fmt_ptr ){
        int i, isample;
        for (isample=0; isample<line->n_sample; isample++){
            if (stats_mappingS[isample] != 1) continue;
            int sample_phased = 0;
            #define BRANCH_INT(type_t,vector_end) { \
                type_t *p = (type_t*) (fmt_ptr->p + isample*fmt_ptr->size); \
                for (i=0; i<fmt_ptr->n; i++) \
                { \
                    if (fmt_ptr->n == 1 || (p[i] == vector_end && i == 1)) { sample_phased = 1; break; } /* haploid phased by definition */ \
                    if ( p[i] == vector_end ) { break; }; /* smaller ploidy */ \
                    if ( bcf_gt_is_missing(p[i]) ) continue; /* missing allele */ \
                    if ((p[i])&1) { \
                        sample_phased = 1; \
                        break; \
                    } \
                } \
            }
            switch (fmt_ptr->type) {
                case BCF_BT_INT8:  BRANCH_INT(int8_t,  bcf_int8_vector_end); break;
                case BCF_BT_INT16: BRANCH_INT(int16_t, bcf_int16_vector_end); break;
                case BCF_BT_INT32: BRANCH_INT(int32_t, bcf_int32_vector_end); break;
                default: fprintf(stderr, "[E::%s] todo: fmt_type %d\n", __func__, fmt_ptr->type); exit(1); break;
            }
            #undef BRANCH_INT
            if (!sample_phased) {
                all_phased = 0;
                break;
            }
        }
    }else all_phased = 0;
    return all_phased;
}


void rtc_data::readSampleExclusionStats(string fname){
    vrb.title("Read sample exclusion list [" + fname + "]");
    int ret = stats_vcf_sample_filter.readExclusion(fname);
    if (ret < 0) vrb.error("Cannot open file!");
    else vrb.bullet(stb.str(ret) + " samples");
}

void rtc_data::readSampleInclusionStats(string fname){
    vrb.title("Read sample inclusion list [" + fname + "]");
    int ret = stats_vcf_sample_filter.readInclusion(fname);
    if (ret < 0) vrb.error("Cannot open file!");
    else vrb.bullet(stb.str(ret) + " samples");
}

void rtc_data::copyIncludeExclude(){
    stats_vcf_sample_filter = filter_sample;
}

void rtc_data::setStatsVCF(string file){
    stats_vcf_file = file;
}

long long unsigned int rtc_data::getMemoryUsage(){
	long long unsigned int mem = 0;
	ifstream file("/proc/self/status");
	if (file.is_open()){
	    string buffer;
	    vector < string > str;
	    while(getline(file, buffer)) {
	        stb.split(buffer, str);
	        if (str[0] == "VmRSS:"){
	        	mem += atoi(str[1].c_str()) * 1024;
	        	break;
	        }
	    }
	    file.close();
	}else mem += (genotype_val.capacity() * genotype_val[0].capacity() + phenotype_val.capacity() * phenotype_val[0].capacity()) * sizeof(float);

	long long unsigned int max_variants = 0;
	for (int c = 0 ; c < all_coldspots.size(); c++){
		if (all_coldspots[c].coldspot_variant_idx.size() > max_variants) max_variants = all_coldspots[c].coldspot_variant_idx.size();
		mem += all_coldspots[c].getMemoryUsage();
	}
	cerr << max_variants << endl;
	if(DprimeR2inMem) mem += ((sample_iterations * 2 * max_variants * phenotype_val[0].capacity() * sizeof(float))  + (max_variants * genotype_val[0].capacity() * sizeof(float)) );
	mem += max_variants * (max_variants - 1) / 2 * sizeof(float);
	return mem;
}
