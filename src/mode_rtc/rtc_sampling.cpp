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

vector <double> rtc_data::sampleRTC(vector <int> & genotype_idx, vector <float> & phenotype, int eqtl_idx, double RTC, int pI){
    vector < double > result = {__RTC_NA__,__RTC_NA__};
    if (genotype_idx.size() < 4) return result;
    //double slope = 0.0;
    //regression(genotype_val[eqtl_idx], phenotype, slope);
    linReg linreg(genotype_val[eqtl_idx], phenotype);
    double count = 0.0;
    double better = 0.0;
    double count2 = 0.0;
    double worse = 0.0;
    unordered_map < int , vector < vector < float > > > pseudo_phenos;
    unordered_map < int , map < int , map < int, double > > > betters;
    unordered_map < int , map < int , map < int, double > > > worses;
    int better_hit = 0, pseudo_hit = 0;
    float medianR2;
    //calculate median r2
    if (pI >= 0){
		vector < float > r2s;
		for (int i =0 ; i < genotype_idx.size(); i++){
			for (int j =i+1; j < genotype_idx.size(); j++){
				r2s.push_back(getRsquare(genotype_idx[i], genotype_idx[j]));
			}
		}
		medianR2 = median(r2s);
        r2s.clear();
    }
    ////////////////////
    for (int si = 0 ; si < sample_iterations; si++){
        
        //select random causal eQTL
        int r_eqtl_causal = genotype_idx[rng.getInt(genotype_idx.size())];
        
        //Find variants linked to this random eQTL
        vector <int> possible_selections1;
        for (int s = 0 ; s < genotype_idx.size(); s++){
            if (genotype_idx[s] == r_eqtl_causal) continue;
            if (getRsquare(r_eqtl_causal, genotype_idx[s]) >= R2_cutoff) possible_selections1.push_back(genotype_idx[s]);
        }
        
        //If there are no linked ones continue
        if(!possible_selections1.size()) {
        	if (pI >= 0) cerr << phenotype_id[pI] << "-" << genotype_id[eqtl_idx] << " H0 " << genotype_id[r_eqtl_causal] << " " << si << " " << sample_iterations << " 0 NA 0 NA NA " << better << " " << count << " " << RTC << " " << medianR2<< " NA NA NA NA" <<endl;
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
            if (genotype_idx[s] != r_eqtl_causal && genotype_idx[s] != r_eqtl && getRsquare(r_other_causal, genotype_idx[s]) >= R2_cutoff) possible_selections2.push_back(genotype_idx[s]);
        }
        
        //If there are no linked ones continue
        if(!possible_selections2.size()) {
        	if (pI >= 0) cerr << phenotype_id[pI] << "-" << genotype_id[eqtl_idx] << " H0 " << genotype_id[r_eqtl_causal] << " " << si << " " << sample_iterations << " " << possible_selections1.size() << " " << genotype_id[r_eqtl] << " 0 NA NA " << better << " " << count << " " << RTC << " " << medianR2 << " NA NA NA NA"<< endl;
        	continue;
        }
        //Select a random linked one
        int r_other = possible_selections2[rng.getInt(possible_selections2.size())];
        count++;

        if(betters.count(r_eqtl_causal) && betters[r_eqtl_causal].count(r_eqtl) && betters[r_eqtl_causal][r_eqtl].count(r_other)){
        	better_hit++;
        	if (betters[r_eqtl_causal][r_eqtl][r_other] >= RTC) better++;
        	if (pI >= 0) cerr << phenotype_id[pI] << "-" << genotype_id[eqtl_idx] << " H0 " << genotype_id[r_eqtl_causal] << " " << si << " " << sample_iterations << " " << possible_selections1.size() << " " << genotype_id[r_eqtl] << " " << possible_selections2.size() << " " << genotype_id[r_other] << " " <<  betters[r_eqtl_causal][r_eqtl][r_other] << " " << better << " " << count << " " << RTC << " " << medianR2 << " " << getRsquare(r_eqtl_causal, r_other_causal) << " " << getRsquare(r_eqtl_causal, r_eqtl) << " " << getRsquare(r_other_causal, r_other) << " " <<  getRsquare(r_other, r_eqtl) << endl;
        }else{
			vector < double > corrs(genotype_idx.size());
			double test_snp_corr = 0.0 ;
			vector < float > genotype_eqtl;
			if (genotypeSink.count(r_eqtl)){
				genotype_eqtl = genotypeSink[r_eqtl];
			}else{
				genotype_eqtl = genotype_val[r_eqtl];
				normalize(genotype_eqtl);
				if (DprimeR2inMem) genotypeSink[r_eqtl] = genotype_eqtl;
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
						if (DprimeR2inMem) genotypeSink[genotype_idx[s]] = test;
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
			if (rtc >= RTC) better++;
			if (DprimeR2inMem) betters[r_eqtl_causal][r_eqtl][r_other]= rtc;
        	if (pI >= 0) cerr << phenotype_id[pI] << "-" << genotype_id[eqtl_idx] << " H0 " << genotype_id[r_eqtl_causal] << " " << si << " " << sample_iterations << " " << possible_selections1.size() << " " << genotype_id[r_eqtl] << " " << possible_selections2.size() << " " << genotype_id[r_other] << " " << rtc << " " << better << " " << count << " " << RTC << " " << medianR2 << " " << getRsquare(r_eqtl_causal, r_other_causal) << " " << getRsquare(r_eqtl_causal, r_eqtl) << " " << getRsquare(r_other_causal, r_other) << " " <<  getRsquare(r_other, r_eqtl) << endl;
        }
        
    }
    
    for (int si = 0 ; si < sample_iterations; si++){
        
        //select random causal eQTL
        int r_eqtl_causal = genotype_idx[rng.getInt(genotype_idx.size())];
        
        //Find variants linked to this true eQTL
        vector <int> possible_selections1;
        for (int s = 0 ; s < genotype_idx.size(); s++){
            if (genotype_idx[s] == r_eqtl_causal) continue;
            if (getRsquare(r_eqtl_causal, genotype_idx[s]) >= R2_cutoff) possible_selections1.push_back(genotype_idx[s]);
        }
        
        //If there are less than 2 linked ones continue
        if(possible_selections1.size() < 2) {
            if (pI >= 0) cerr << phenotype_id[pI] << "-" << genotype_id[eqtl_idx] << " H1 " << genotype_id[r_eqtl_causal] << " " << si << " " << sample_iterations << " " << possible_selections1.size() << " NA NA NA NA " << worse << " " << count2 << " " << RTC << " " << medianR2<< " NA NA NA NA"<<endl;
            continue;
        }
        unsigned int rngi = rng.getInt(possible_selections1.size());
        int r_eqtl = possible_selections1[rngi];
        possible_selections1.erase(possible_selections1.begin() + rngi);
        int r_other = possible_selections1[rng.getInt(possible_selections1.size())];
        count2++;
        if(worses.count(r_eqtl_causal) && worses[r_eqtl_causal].count(r_eqtl) && worses[r_eqtl_causal][r_eqtl].count(r_other)){
            better_hit++;
            if (worses[r_eqtl_causal][r_eqtl][r_other] < RTC) worse++;
            if (pI >= 0) cerr << phenotype_id[pI] << "-" << genotype_id[eqtl_idx] << " H1 " << genotype_id[r_eqtl_causal] << " " << si << " " << sample_iterations << " " << possible_selections1.size() << " " << genotype_id[r_eqtl] << " " << "NA" << " " << genotype_id[r_other] << " " <<  worses[r_eqtl_causal][r_eqtl][r_other] << " " << worse << " " << count2 << " " << RTC << " " << medianR2 << " NA " << getRsquare(r_eqtl_causal, r_eqtl) << " " << getRsquare(r_eqtl_causal, r_other) << " " <<  getRsquare(r_other, r_eqtl) << endl;
        }else{
        
            vector < double > corrs(genotype_idx.size());
            double test_snp_corr = 0.0 ;
            vector < float > genotype_eqtl;
            if (genotypeSink.count(r_eqtl)){
                genotype_eqtl = genotypeSink[r_eqtl];
            }else{
                genotype_eqtl = genotype_val[r_eqtl];
                normalize(genotype_eqtl);
                if (DprimeR2inMem) genotypeSink[r_eqtl] = genotype_eqtl;
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
                        if (DprimeR2inMem) genotypeSink[genotype_idx[s]] = test;
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
            if (rtc < RTC) worse++;
            if (DprimeR2inMem) worses[r_eqtl_causal][r_eqtl][r_other]= rtc;
            if (pI >= 0) cerr << phenotype_id[pI] << "-" << genotype_id[eqtl_idx] << " H1 " << genotype_id[r_eqtl_causal] << " " << si << " " << sample_iterations << " " << possible_selections1.size() << " " << genotype_id[r_eqtl] << " " << "NA" << " " << genotype_id[r_other] << " " << rtc << " " << worse << " " << count2 << " " << RTC << " " << medianR2 << " NA " << getRsquare(r_eqtl_causal, r_eqtl) << " " << getRsquare(r_eqtl_causal, r_other) << " " <<  getRsquare(r_other, r_eqtl) << endl;
        }
        
    }
    //if (DprimeR2inMem) vrb.bullet(stb.str(better_hit) + " " + stb.str(pseudo_hit));
    result = { (better + 1.0) / (count + 1.0) , (worse + 1.0) / (count2 + 1.0) };
    return result;
}


