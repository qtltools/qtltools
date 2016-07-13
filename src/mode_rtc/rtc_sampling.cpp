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

void rtc_data::generatePhenotype (double slope ,int geno_idx, vector <float> & new_pheno){
    double fake_intercept = 100.0;
    new_pheno = vector <float>(sample_count,0.0);
    for (int i = 0 ; i < sample_count; i++) new_pheno[i] = slope * genotype_val[geno_idx][i] + fake_intercept + rnorm(0.0, 1.0);
}

void rtc_data::generatePhenotype(vector<float> &X, linReg &linreg, vector <float> &np){
	int N = X.size();
	np = vector < float > (N,0.0);
    vector < float > residuals = linreg.residuals;
    random_shuffle(residuals.begin(),residuals.end());
    for (int i = 0 ; i<N; i++) np[i] = X[i] * linreg.beta + linreg.yIntercept + residuals[i];
}

rtc_sample_results rtc_data::sampleRTC(vector <int> & genotype_idx, vector <float> & phenotype, int eqtl_idx, double RTC, int pI){
	rtc_sample_results result;
    if (genotype_idx.size() < 4) return result;
    //double slope = 0.0;
    //regression(genotype_val[eqtl_idx], phenotype, slope);
    linReg linreg(genotype_val[eqtl_idx], phenotype);
    double count = 0.0;
    double count2 = 0.0;
    double gt_h0 = 0.0, gtoe_h0 = 0.0 , gt_h1 = 0.0, gtoe_h1 = 0.0;
    unordered_map < int , vector < vector < float > > > pseudo_phenos;
    unordered_map < int , map < int , map < int, double > > > h0s;
    unordered_map < int , map < int , map < int, double > > > h1s;
    set < string > h0ss,h1ss;
    int better_hit = 0, pseudo_hit = 0;
    float medianR2;
    vector < float > h0,h1;
    long int trials = 0;
    //calculate median r2
    //if (pI >= 0){
		vector < float > r2s;
		for (int i =0 ; i < genotype_idx.size(); i++){
			for (int j =i+1; j < genotype_idx.size(); j++){
				r2s.push_back(getRsquare(genotype_idx[i], genotype_idx[j],genotype_idx[0],genotype_idx.size()));
			}
		}
		medianR2 = median(r2s);
        vector < float > ().swap(r2s);
    //}
    ////////////////////
    while ( count < sample_iterations && trials < max_sample_iterations){
        trials++;
        //select random causal eQTL
        int r_eqtl_causal = genotype_idx[rng.getInt(genotype_idx.size())];
        
        //Find variants linked to this random eQTL
        vector <int> possible_selections1;
        for (int s = 0 ; s < genotype_idx.size(); s++){
            if (genotype_idx[s] == r_eqtl_causal) continue;
            if (getRsquare(r_eqtl_causal, genotype_idx[s],genotype_idx[0],genotype_idx.size()) >= R2_cutoff) possible_selections1.push_back(genotype_idx[s]);
        }
        
        //If there are no linked ones continue
        if(!possible_selections1.size()) {
        	if (pI >= 0) cerr << phenotype_id[pI] << "-" << genotype_id[eqtl_idx] << " H0 " << genotype_id[r_eqtl_causal] << " " << count << " " << sample_iterations << " 0 NA 0 NA NA " << gtoe_h0 << " " << gt_h0 << " " << count << " " << RTC << " " << medianR2<< " NA NA NA NA" <<endl;
        	continue;
        }
        //Select a random linked one
        int r_eqtl = possible_selections1[rng.getInt(possible_selections1.size())];
        
        //select second random causal
        int r_other_causal = r_eqtl_causal;
        while(r_other_causal == r_eqtl_causal || r_other_causal == r_eqtl) r_other_causal = genotype_idx[rng.getInt(genotype_idx.size())];
        
        //Find variants linked to second random causal
        vector <int> possible_selections2;
        for (int s = 0 ; s < genotype_idx.size(); s++){
            if (genotype_idx[s] == r_other_causal) continue;
            if (genotype_idx[s] != r_eqtl_causal && genotype_idx[s] != r_eqtl && getRsquare(r_other_causal, genotype_idx[s],genotype_idx[0],genotype_idx.size()) >= R2_cutoff) possible_selections2.push_back(genotype_idx[s]);
        }
        
        //If there are no linked ones continue
        if(!possible_selections2.size()) {
        	if (pI >= 0) cerr << phenotype_id[pI] << "-" << genotype_id[eqtl_idx] << " H0 " << genotype_id[r_eqtl_causal] << " " << count << " " << sample_iterations << " " << possible_selections1.size() << " " << genotype_id[r_eqtl] << " 0 NA NA " <<  gtoe_h0 << " " << gt_h0 << " " << count << " " << RTC << " " << medianR2 << " NA NA NA NA"<< endl;
        	continue;
        }
        //Select a random linked one
        int r_other = possible_selections2[rng.getInt(possible_selections2.size())];
        count++;
        h0ss.insert(stb.str(r_eqtl_causal) + stb.str(r_eqtl) + stb.str(r_other));
        if(h0s.count(r_eqtl_causal) && h0s[r_eqtl_causal].count(r_eqtl) && h0s[r_eqtl_causal][r_eqtl].count(r_other)){
        	better_hit++;
        	h0.push_back(h0s[r_eqtl_causal][r_eqtl][r_other]);
        	if (h0s[r_eqtl_causal][r_eqtl][r_other] >= RTC) {
        		gtoe_h0++;
        		if (h0s[r_eqtl_causal][r_eqtl][r_other] > RTC) gt_h0++;
        	}
        	if (pI >= 0) cerr << phenotype_id[pI] << "-" << genotype_id[eqtl_idx] << " H0 " << genotype_id[r_eqtl_causal] << " " << count << " " << sample_iterations << " " << possible_selections1.size() << " " << genotype_id[r_eqtl] << " " << possible_selections2.size() << " " << genotype_id[r_other] << " " <<  h0s[r_eqtl_causal][r_eqtl][r_other] << " " << gtoe_h0 << " " << gt_h0 << " " << count << " " << RTC << " " << medianR2 << " " << getRsquare(r_eqtl_causal, r_other_causal,genotype_idx[0],genotype_idx.size()) << " " << getRsquare(r_eqtl_causal, r_eqtl,genotype_idx[0],genotype_idx.size()) << " " << getRsquare(r_other_causal, r_other,genotype_idx[0],genotype_idx.size()) << " " <<  getRsquare(r_other, r_eqtl,genotype_idx[0],genotype_idx.size()) << endl;
        }else{
			vector < double > corrs(genotype_idx.size());
			double test_snp_corr = 0.0 ;
			vector < float > genotype_eqtl;
			if (genotypeSink.count(r_eqtl)){
				genotype_eqtl = genotypeSink[r_eqtl];
			}else{
				genotype_eqtl = genotype_val[r_eqtl];
				normalize(genotype_eqtl);
				if (DprimeR2inMem >= 2) genotypeSink[r_eqtl] = genotype_eqtl;
			}
			if (pseudo_phenos.count(r_eqtl_causal)){
				pseudo_hit++;
				for (int s = 0 ; s < genotype_idx.size() ; s++){
					corrs[s] = abs(getCorrelation(genotype_eqtl, pseudo_phenos[r_eqtl_causal][s]));
					if( genotype_idx[s] == r_other) test_snp_corr = corrs[s];
				}
			}else{
				vector < float > pseudo_pheno;
				//generatePhenotype(slope, r_eqtl_causal, pseudo_pheno);
				generatePhenotype(genotype_val[r_eqtl_causal], linreg, pseudo_pheno);
				if (options.count("normal")) normalTransform(pseudo_pheno);
				normalize(pseudo_pheno);
				for (int s = 0 ; s < genotype_idx.size() ; s++){
					vector < float > test;
					if(genotypeSink.count(genotype_idx[s])){
						test = genotypeSink[genotype_idx[s]];
					}else{
						test = genotype_val[genotype_idx[s]];
						normalize(test);
						if (DprimeR2inMem >= 2) genotypeSink[genotype_idx[s]] = test;
					}
					vector <float> new_pheno = correct(test,pseudo_pheno);
					if (options.count("normal")) normalTransform(new_pheno);
					normalize(new_pheno);
					if (DprimeR2inMem) pseudo_phenos[r_eqtl_causal].push_back(new_pheno);
					corrs[s] = abs(getCorrelation(genotype_eqtl, new_pheno));
					if( genotype_idx[s] == r_other) test_snp_corr = corrs[s];
				}
			}

			sort(corrs.begin(),corrs.end());
			int rank = -1;
			for (int i = 0 ; i<corrs.size() && corrs[i] <= test_snp_corr; i++) if(corrs[i] == test_snp_corr) rank = i;
			double rtc = ((double) corrs.size() - (double) rank) / (double) corrs.size();
			h0.push_back(rtc);
			if (rtc >= RTC) {
				gtoe_h0++;
				if (rtc > RTC) gt_h0++;
			}
			if (DprimeR2inMem) h0s[r_eqtl_causal][r_eqtl][r_other]= rtc;
        	if (pI >= 0) cerr << phenotype_id[pI] << "-" << genotype_id[eqtl_idx] << " H0 " << genotype_id[r_eqtl_causal] << " " << count << " " << sample_iterations << " " << possible_selections1.size() << " " << genotype_id[r_eqtl] << " " << possible_selections2.size() << " " << genotype_id[r_other] << " " << rtc << " " << gtoe_h0 << " " << gt_h0 << " " << count << " " << RTC << " " << medianR2 << " " << getRsquare(r_eqtl_causal, r_other_causal,genotype_idx[0],genotype_idx.size()) << " " << getRsquare(r_eqtl_causal, r_eqtl,genotype_idx[0],genotype_idx.size()) << " " << getRsquare(r_other_causal, r_other,genotype_idx[0],genotype_idx.size()) << " " <<  getRsquare(r_other, r_eqtl,genotype_idx[0],genotype_idx.size()) << endl;
        }
        
    }
    trials = 0;
    while ( count2 < sample_iterations && trials < max_sample_iterations ){
        trials++;
        //select random causal eQTL
        int r_eqtl_causal = genotype_idx[rng.getInt(genotype_idx.size())];
        
        //Find variants linked to this true eQTL
        vector <int> possible_selections1;
        for (int s = 0 ; s < genotype_idx.size(); s++){
            if (genotype_idx[s] == r_eqtl_causal) continue;
            if (getRsquare(r_eqtl_causal, genotype_idx[s],genotype_idx[0],genotype_idx.size()) >= R2_cutoff) possible_selections1.push_back(genotype_idx[s]);
        }
        
        //If there are less than 2 linked ones continue
        if(possible_selections1.size() < 2) {
            if (pI >= 0) cerr << phenotype_id[pI] << "-" << genotype_id[eqtl_idx] << " H1 " << genotype_id[r_eqtl_causal] << " " << count2 << " " << sample_iterations << " " << possible_selections1.size() << " NA NA NA NA " << gtoe_h1 << " " << gt_h1 << " " << count2 << " " << RTC << " " << medianR2<< " NA NA NA NA"<<endl;
            continue;
        }
        unsigned int rngi = rng.getInt(possible_selections1.size());
        int r_eqtl = possible_selections1[rngi];
        possible_selections1.erase(possible_selections1.begin() + rngi);
        int r_other = possible_selections1[rng.getInt(possible_selections1.size())];
        count2++;
        h1ss.insert(stb.str(r_eqtl_causal) + stb.str(r_eqtl) + stb.str(r_other));
        if(h1s.count(r_eqtl_causal) && h1s[r_eqtl_causal].count(r_eqtl) && h1s[r_eqtl_causal][r_eqtl].count(r_other)){
            better_hit++;
            h1.push_back(h1s[r_eqtl_causal][r_eqtl][r_other]);
            if (h1s[r_eqtl_causal][r_eqtl][r_other] >= RTC) {
            	gtoe_h1++;
				if(h1s[r_eqtl_causal][r_eqtl][r_other] > RTC) gt_h1++;
            }
            if (pI >= 0) cerr << phenotype_id[pI] << "-" << genotype_id[eqtl_idx] << " H1 " << genotype_id[r_eqtl_causal] << " " << count2 << " " << sample_iterations << " " << possible_selections1.size() << " " << genotype_id[r_eqtl] << " " << "NA" << " " << genotype_id[r_other] << " " <<  h1s[r_eqtl_causal][r_eqtl][r_other] << " " << gtoe_h1 << " " << gt_h1  << " " << count2 << " " << RTC << " " << medianR2 << " NA " << getRsquare(r_eqtl_causal, r_eqtl,genotype_idx[0],genotype_idx.size()) << " " << getRsquare(r_eqtl_causal, r_other,genotype_idx[0],genotype_idx.size()) << " " <<  getRsquare(r_other, r_eqtl,genotype_idx[0],genotype_idx.size()) << endl;
        }else{
        
            vector < double > corrs(genotype_idx.size());
            double test_snp_corr = 0.0 ;
            vector < float > genotype_eqtl;
            if (genotypeSink.count(r_eqtl)){
                genotype_eqtl = genotypeSink[r_eqtl];
            }else{
                genotype_eqtl = genotype_val[r_eqtl];
                normalize(genotype_eqtl);
                if (DprimeR2inMem >= 2) genotypeSink[r_eqtl] = genotype_eqtl;
            }
            if (pseudo_phenos.count(r_eqtl_causal)){
                pseudo_hit++;
                for (int s = 0 ; s < genotype_idx.size() ; s++){
                    corrs[s] = abs(getCorrelation(genotype_eqtl, pseudo_phenos[r_eqtl_causal][s]));
                    if( genotype_idx[s] == r_other) test_snp_corr = corrs[s];
                }
            }else{
                vector < float > pseudo_pheno;
                //generatePhenotype(slope, r_eqtl_causal, pseudo_pheno);
                generatePhenotype(genotype_val[r_eqtl_causal], linreg, pseudo_pheno);
                if (options.count("normal")) normalTransform(pseudo_pheno);
                normalize(pseudo_pheno);
                for (int s = 0 ; s < genotype_idx.size() ; s++){
                    vector < float > test;
                    if(genotypeSink.count(genotype_idx[s])){
                        test = genotypeSink[genotype_idx[s]];
                    }else{
                        test = genotype_val[genotype_idx[s]];
                        normalize(test);
                        if (DprimeR2inMem >= 2) genotypeSink[genotype_idx[s]] = test;
                    }
                    vector <float> new_pheno = correct(test,pseudo_pheno);
                    if (options.count("normal")) normalTransform(new_pheno);
                    normalize(new_pheno);
                    if (DprimeR2inMem) pseudo_phenos[r_eqtl_causal].push_back(new_pheno);
                    corrs[s] = abs(getCorrelation(genotype_eqtl, new_pheno));
                    if( genotype_idx[s] == r_other) test_snp_corr = corrs[s];
                }
            }

            sort(corrs.begin(),corrs.end());
            int rank = -1;
            for (int i = 0 ; i<corrs.size() && corrs[i] <= test_snp_corr; i++) if(corrs[i] == test_snp_corr) rank = i;
            double rtc = ((double) corrs.size() - (double) rank) / (double) corrs.size();
            h1.push_back(rtc);
            if (rtc >= RTC) {
            	gtoe_h1++;
            	if(rtc > RTC) gt_h1++;
            }
            if (DprimeR2inMem) h1s[r_eqtl_causal][r_eqtl][r_other]= rtc;
            if (pI >= 0) cerr << phenotype_id[pI] << "-" << genotype_id[eqtl_idx] << " H1 " << genotype_id[r_eqtl_causal] << " " << count2 << " " << sample_iterations << " " << possible_selections1.size() << " " << genotype_id[r_eqtl] << " " << "NA" << " " << genotype_id[r_other] << " " << rtc << " " << gtoe_h1 << " " << gt_h1  << " " << count2 << " " << RTC << " " << medianR2 << " NA " << getRsquare(r_eqtl_causal, r_eqtl) << " " << getRsquare(r_eqtl_causal, r_other) << " " <<  getRsquare(r_other, r_eqtl) << endl;
        }
        
    }
    //if (DprimeR2inMem) vrb.bullet(stb.str(better_hit) + " " + stb.str(pseudo_hit));
    result.gtoe_h0 =  gtoe_h0;
    result.gt_h0 =  gt_h0;
    result.gtoe_h1 =  gtoe_h1;
    result.gt_h1 =  gt_h1;
    result.count_h0 =  count;
    result.count_h1 =  count2;
    result.unique_h0 = h0ss.size();
    result.unique_h1 = h1ss.size();
    result.medianR2 = medianR2;
    result.median_h0 = median(h0);
    result.median_h1 = median(h1);
	result.pval = (gtoe_h0 + 1.0) / (count + 1.0);
	const char * sep = ",";
	stringstream h0sss,h1sss;
	if(h0.size()){
		copy(h0.begin(),h0.end(),ostream_iterator<float>(h0sss,sep));
		result.h0 = h0sss.str();
		//result.h0.pop_back();
		result.h0.resize (result.h0.size () - 1);
	} else result.h0 = "NA";
	if(h1.size()){
		copy(h1.begin(),h1.end(),ostream_iterator<float>(h1sss,sep));
		result.h1 = h1sss.str();
		result.h1.resize (result.h1.size () - 1);
	} else result.h1 = "NA";
	if (count && count2){
		probability(h0,h1,RTC, result);
    }
    unordered_map < int , vector < vector < float > > >().swap(pseudo_phenos);
    unordered_map < int , map < int , map < int, double > > >().swap(h0s);
    unordered_map < int , map < int , map < int, double > > >().swap(h1s);
    vector < float > ().swap(h1);
    vector < float > ().swap(h0);
    return result;
}

void rtc_data::probability(vector < float > &h0, vector < float > &h1, double RTC, rtc_sample_results & res){
	vector < float > all = h0;
	all.insert(all.end(),h1.begin(),h1.end());
	sort(all.begin(),all.end());
	int step = all.size() / 10 ;
	double diff = 2.0;
	int index = 0;
	for (unsigned long int i = 0 ; i < all.size(); i++){
		double d = abs(RTC - all[i]);
		if (d < diff){
			diff = d;
			index = i;
		}
	}
	int s = index - step >= 0 ? index - step : 0;
	int e = index + step >= all.size() ? all.size() - 1 : index+step;
	res.rtc_bin_start = all[s];
	res.rtc_bin_end = all[e];
	double c_h0 = 0.0 , c_h1 = 0.0;
	for (int i = 0 ; i < h0.size() ; i++) if (h0[i] >= all[s] && h0[i] <= all[e]) c_h0++;
	for (int i = 0 ; i < h1.size() ; i++) if (h1[i] >= all[s] && h1[i] <= all[e]) c_h1++;
	res.rtc_bin_h0_proportion = c_h0 / h0.size();
	res.rtc_bin_h1_proportion = c_h1 / h1.size();
}
