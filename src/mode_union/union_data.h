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

#ifndef _UNION_DATA_H
#define _UNION_DATA_H

//INCLUDES
#include "../common/data.h"

class coldspot_u{
public:
    string chr;
    int start;
    int end;
    int idx;
    string type;
    coldspot_u(){chr="";start=-1;end=-1;idx=-1;type="NA";}
    coldspot_u(string c, int s, int e, int i,string t){chr=c;start=s;end=e;idx=i; type = t ;}
    friend ostream& operator<<(ostream& out, const coldspot_u& p){
        out << "ID: " << p.idx << " CHR: " << p.chr << " START: " << p.start << " END: " << p.end << " TYPE: " << p.type << endl;
        return out;
    }
};

class results{
public:
    vector < string > pheno;
    vector < string > geno;
    vector < int > rank;
    vector < double > pval;
    vector < int > csi;
    vector < string > regions;
    void assign(string p, string s, int r, double pv,int cs , string re){pheno.push_back(p);geno.push_back(s);rank.push_back(r);pval.push_back(pv);csi.push_back(cs); regions.push_back(re);}
    friend ostream& operator<<(ostream& out, const results& r){
        for (int i =0 ; i < r.pheno.size(); i++) out << "__UNION__ " << r.pheno[i] << " " << r.geno[i] << " " << r.rank[i] << " 1 " << r.pval[i] << " "<< r.csi[i] << " " << r.regions[i] << "\n";
        return out;
    }
};


struct genotypes_holder{
	vector < vector < float> > genotypes;
	vector < string > ids;
};

class myPhenotype{
public:
	vector < vector < int > > ranks;
	vector < vector < string > > genotypes;
	vector < bool > found;
	int max_independent_signal;
	string pheno_chr;
	int pheno_pos;
	myPhenotype(){pheno_chr = "";pheno_pos = -1; max_independent_signal=0;}
	~myPhenotype(){found.clear();ranks.clear();genotypes.clear();}
	myPhenotype(string c , int p , int f){ranks = vector < vector < int > >(f); genotypes = vector < vector < string > >(f); pheno_chr = c ; pheno_pos = p; found = vector < bool >(f,false);max_independent_signal=1;}
	void assign(string g , int r , int i){genotypes[i].push_back(g); ranks[i].push_back(r); found[i] = true; if(ranks[i].size() > max_independent_signal) max_independent_signal = ranks[i].size();}
};


class union_data : public data {
public:
	//PARAMETERS
	int pvalue_column,variant_column,phenotype_column,rank_column,best_column,coldspot_count;
	static const int bin_size = 1000000;
    int no_of_files;


	//REGIONS
	genomic_region regionPhenotype;
	genomic_region regionGenotype;

	//GENOTYPES
	int genotype_count;									//variant site number
	vector < vector < float > > genotype_val;			//variant site genotype dosages
	vector < string > genotype_chr;						//variant site chromosome
	vector < string > genotype_id;						//variant site IDs
	vector < int > genotype_start;						//variant site start positions
	vector < int > genotype_end;						//variant site end positions
	unordered_map < string, int > genotype_id_to_idx;

	//PHENOTYPES
	int phenotype_count;								//phenotype number
	vector < vector < float > > phenotype_val;			//phenotype values
	vector < string > phenotype_id;						//phenotype ids
	vector < string > phenotype_grp;					//phenotype groups
	vector < string > phenotype_chr;					//phenotype chromosomes
	vector < int > phenotype_start;						//phenotype start positions
	vector < int > phenotype_end;						//phenotype end positions
    unordered_map < string, int > phenotype_id_to_idx;

	//COVARIATES & INTERACTION
	int covariate_count;								//covariate number
	vector < vector < string > > covariate_val;			//covariate values
	vector < string > covariate_id;						//covariate ids

	//RTC
    map < string, map < int,  vector <coldspot_u *> > >  coldspot_bins_p;
    vector < coldspot_u *> all_coldspots_p;
    map < int , map < string , myPhenotype> > toUnite;

	//CONSTRUCTOR / DESTRUCTOR
	union_data();
	~union_data();
	void clear();
    void clearNotHotspot();
	void clearSamples();

	//DATA REGION
	bool setPhenotypeRegion(string);
	bool setGenotypeRegion(string);
    void deduceGenotypeRegion(int);
	void setPhenotypeRegion(int, int);

	//READ DATA
	void readGenotypes(string,string);
	void readGenotypesVCF(string,string);
	void readGenotypesBED(string,string);
    void scanGenotypes(string);
	void scanGenotypesVCF(string);
	void scanGenotypesBED(string);
	void readPhenotypes(string,string);
	void scanPhenotypes(string);
	void readCovariates(string);

	//GENOTYPE & PHENOTYPE MANAGEMENT
	void clusterizePhenotypes(int);
	void imputeGenotypes();
    void imputePhenotypes();
    void residualizePhenotypes();
    void normalTransformPhenotypes();
    void normalTransform(vector < float > &);
	void normalize(vector < float > &);
	void normalize(vector < vector < float > > &);

	//COMPUTATION METHODS [ALL INLINES FOR SPEED]
	double getCorrelation(vector < float > &, vector < float > &);
    double getCorrelation(vector < float > &, vector < float > &, int);
	double getPvalue(double, double);

	//ANALYSIS
    void readHotspots(string);
    int getColdspot(string, int);
    void mapVariantsToColdspots();
    string getBestVariant(genotypes_holder&, int, double &);
    void unions(string,int);
    void unions_conditional(string,int);
    void create_unions(vector <string> & , vector <string> & , vector <string> & , vector <string> &);
    void find_unions(vector <string> & , vector <string> & , vector <string> & , vector <string> &);

};

//***************************************************************//
//******************** DECLARE FUNCTIONS ************************//
//***************************************************************//
void union_main(vector < string > &);

//***************************************************************//
//******************** INLINE FUNCTIONS *************************//
//***************************************************************//

inline double union_data::getCorrelation(vector < float > & vec1, vector < float > & vec2) {
	int i = 0;
	int repeat = (sample_count / 4);
	int left = (sample_count % 4);
	double sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0;

	while (repeat --) {
		sum0 += vec1[i] * vec2[i];
		sum1 += vec1[i+1] * vec2[i+1];
		sum2 += vec1[i+2] * vec2[i+2];
		sum3 += vec1[i+3] * vec2[i+3];
		i += 4;
	}

	switch (left) {
	case 3:	sum0 += vec1[i+2] * vec2[i+2];
	case 2:	sum0 += vec1[i+1] * vec2[i+1];
	case 1:	sum0 += vec1[i+0] * vec2[i+0];
	case 0: ;
	}

	return sum0 + sum1 + sum2 + sum3;
}

inline double union_data::getCorrelation(vector < float > & vec1, vector < float > & vec2, int sc) {
    int i = 0;
    int repeat = (sc / 4);
    int left = (sc % 4);
    double sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0;

    while (repeat --) {
        sum0 += vec1[i] * vec2[i];
        sum1 += vec1[i+1] * vec2[i+1];
        sum2 += vec1[i+2] * vec2[i+2];
        sum3 += vec1[i+3] * vec2[i+3];
        i += 4;
    }

    switch (left) {
        case 3:	sum0 += vec1[i+2] * vec2[i+2];
        case 2:	sum0 += vec1[i+1] * vec2[i+1];
        case 1:	sum0 += vec1[i+0] * vec2[i+0];
        case 0: ;
    }

    return sum0 + sum1 + sum2 + sum3;
}

inline double union_data::getPvalue(double corr, double df) {
	return pf(df * corr * corr / (1 - corr * corr), 1, df, 0, 0);
}

inline string union_data::getBestVariant(genotypes_holder &genotype, int phenotype_idx, double &pval){
	vector < float > y = phenotype_val[phenotype_idx];
	if (options.count("normal")) normalTransform(y);
	normalize(y);
	double bestR = 0.0;
	string bestV = "NA";
	int size = 3;
	for (int g = 0; g< genotype.genotypes.size(); g++){
		vector < float > x = genotype.genotypes[g];
		size = x.size();
		normalize(x);
		double R = abs(getCorrelation(x,y));
		if (R > bestR){
			bestR = R;
			bestV = genotype.ids[g];
		}
	}
	pval = getPvalue(bestR,size-2);
	return bestV;
}

#endif
