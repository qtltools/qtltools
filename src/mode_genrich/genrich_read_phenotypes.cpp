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

#include "genrich_data.h"

void genrich_data::readPhenotypes(string fbed) {
	int n_includedP = 0;
	int n_excludedP = 0;
	int n_negativeStrd = 0;
	vector < string > tokens;

	vector < int > vchr, vpos;
	vector < bool > vneg;

	//Open BED file
	vrb.title("Reading phenotype coordinates in [" + fbed + "]");
	htsFile *fp = hts_open(fbed.c_str(),"r");
	if (!fp) vrb.error("Cannot open file");
	kstring_t str = {0,0,0};
	if (hts_getline(fp, KS_SEP_LINE, &str) <= 0 || !str.l || str.s[0] != '#') vrb.error("Cannot read header line!");

	while (hts_getline(fp, KS_SEP_LINE, &str) >= 0) {
		stb.split(string(str.s), tokens);
		if (str.l && str.s[0] != '#') {
			if (tokens.size() < 6) vrb.error("Incorrect number of columns!");
			if (filter_phenotype.check(tokens[3])) {
				int chr_idx = findCHR(tokens[0]);
				if (chr_idx < 0) {
					chr_idx = chromosome_id.size();
					chromosome_id.push_back(tokens[0]);
					chromosome_idx.insert(pair < string , unsigned int > (tokens[0], chr_idx)) ;
				}

				vchr.push_back(chr_idx);
				vpos.push_back((atoi(tokens[1].c_str()) + atoi(tokens[2].c_str()))/2);
				vneg.push_back(tokens[5] == "-");
				if (vneg.back()) n_negativeStrd ++;
				n_includedP ++;
			} else n_excludedP ++;
		}
	}
	if (hts_close(fp)) vrb.error("Cannot properly close file");
	vrb.bullet(stb.str(n_includedP) + " phenotypes included");
	if (n_excludedP > 0) vrb.bullet(stb.str(n_excludedP) + " phenotypes excluded by user");
	if (n_negativeStrd > 0 ) vrb.bullet(stb.str(n_negativeStrd) + " phenotypes are on the negative strand");
	if (n_includedP == 0) vrb.leave("Cannot find phenotypes!");
	vrb.bullet("Detected number of distinct chromosomes = " + stb.str(chromosome_id.size()));
	vrb.bullet("Number of phenotypes on the negative strand = " + stb.str(n_negativeStrd));

	vrb.title("Convert phenotypes coordinates into a set of interval trees");
	vector < vector < Interval < pair < bool, bool > > > > Tvec = vector < vector < Interval < pair < bool, bool > > > > (chromosome_id.size(), vector < Interval < pair < bool, bool > > > ());
	for (int t = 0 ; t < vpos.size() ; t ++) {
		if (t == 0 || vchr[t-1] != vchr[t]) {
			Tvec[vchr[t]].push_back(Interval < pair < bool, bool > > (-1000000000, vpos[t] - 1, pair < bool, bool >(false, vneg[t])));
		} else if ((t == vpos.size() - 1) || vchr[t] != vchr[t+1]) {
			Tvec[vchr[t]].push_back(Interval < pair < bool, bool > > (vpos[t-1], vpos[t] - 1, pair < bool, bool >(vneg[t-1], vneg[t])));
			Tvec[vchr[t]].push_back(Interval < pair < bool, bool > > (vpos[t], 1000000000, pair < bool, bool >(vneg[t], false)));
		} else {
			Tvec[vchr[t]].push_back(Interval < pair < bool, bool > > (vpos[t-1], vpos[t] - 1, pair < bool, bool >(vneg[t-1], vneg[t])));
		}
	}
	phenotype_pos = vector <  IntervalTree < pair < bool, bool > > > (chromosome_id.size(), IntervalTree < pair < bool, bool > > ());
	for (int c = 0 ; c < chromosome_id.size() ; c ++) phenotype_pos[c] = IntervalTree < pair < bool, bool > > (Tvec[c]);
}
