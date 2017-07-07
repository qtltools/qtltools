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


#include "fdensity_data.h"

void fdensity_data::buildIntervalTrees() {

	//0. Initilization
	vrb.title("Mapping annotations to QTLs ");

	//1. Enumerate chr
	map < string, int > chr2idx;
	for (int t = 0 ; t < tss_count ; t ++) {
		map < string, int >::iterator itC = chr2idx.find(tss_chr[t]);
		if (itC == chr2idx.end()) chr2idx.insert(pair < string, int > (tss_chr[t], chr2idx.size()));
	}
	unsigned int n_chr = chr2idx.size();
	vrb.bullet("#detected chromosomes in QTL data = " + stb.str(n_chr));

	//2. Build chromosomal interval trees
	unsigned int chr_unfound = 0;
	vector < vector < Interval < bool > > > Ivec = vector < vector < Interval < bool > > > (n_chr, vector < Interval < bool > > ());
	for (int a = 0 ; a < ann_count ; a ++) {
		map < string, int >::iterator itC = chr2idx.find(ann_chr[a]);
		if (itC == chr2idx.end()) chr_unfound ++;
		else Ivec[itC->second].push_back(Interval < bool > (ann_start[a], ann_end[a], true));
	}
	vrb.bullet("#annotations NOT mapped to QTL chromosomes = " + stb.str(chr_unfound));
	if ((ann_count - chr_unfound) > 0) vrb.bullet("#annotations mapped to QTL chromosomes = " + stb.str(ann_count - chr_unfound));
	else vrb.error("None of the annotations have been found to be located on the same chromosomes than QTL!");

	vector <  IntervalTree < bool > > Ctree = vector <  IntervalTree < bool > > (n_chr, IntervalTree < bool > ());
	for (int c = 0 ; c < n_chr ; c ++) Ctree[c] = IntervalTree < bool > (Ivec[c]);

	//4. Build functional neighborhoods

	int n_neg = 0, n_pos = 0, n_mid = 0;
	int s_neg = 0, s_pos = 0, s_mid = 0;

	basic_stats Rstat;
	Itree = vector <  IntervalTree < bool > > (tss_count);
	for (int t = 0 ; t < tss_count ; t ++) {
		map < string, int >::iterator itC = chr2idx.find(tss_chr[t]);
		assert(itC != chr2idx.end());
		vector < Interval < bool > > ann_in_cis;
		Ctree[itC->second].findOverlapping(tss_pos[t] - window, tss_pos[t] + window, ann_in_cis);

		IntervalStartSorter < bool > intervalStartSorter;
		sort(ann_in_cis.begin(), ann_in_cis.end(), intervalStartSorter);

        Rstat.push(ann_in_cis.size() * 1.0);
		vector < Interval < bool > > Rvec;

		for (int a = 0 ; a < ann_in_cis.size() ; a ++) {

			if (!tss_neg[t]) Rvec.push_back(Interval < bool > (ann_in_cis[a].start - tss_pos[t], ann_in_cis[a].stop - tss_pos[t], true));
			else Rvec.push_back(Interval < bool > (-1 * (ann_in_cis[a].stop - tss_pos[t]), -1 * (ann_in_cis[a].start - tss_pos[t]), true));
			assert(Rvec.back().start <= Rvec.back().stop);

			if (Rvec.back().start < 0 && Rvec.back().stop < 0) {
				s_neg += Rvec.back().stop - Rvec.back().start;
				n_neg ++;
			} else if (Rvec.back().start > 0 && Rvec.back().stop > 0) {
				s_pos += Rvec.back().stop - Rvec.back().start;
				n_pos ++;
			} else {
				s_mid += Rvec.back().stop - Rvec.back().start;
				n_mid ++;
			}
		}
		sort(Rvec.begin(), Rvec.end(), intervalStartSorter);
		Itree[t] = IntervalTree < bool > (Rvec);
	}

	vrb.bullet("#annotated cis-windows = " + stb.str(Rstat.size()));
	vrb.bullet("#annotations per cis-window = " + stb.str(Rstat.mean(), 2) + " +/- " + stb.str(Rstat.sd(), 2));
}

void fdensity_data::runDensityCalculation(string fout) {

	vrb.title("Density analysis");
	output_file fdo (fout.c_str());
	for (int w = -1 * window ; w < window ; w += bin) {
		int wfrom = w;
		int wto = w + bin - 1;
		int n_annotation = 0;


		for (int t = 0 ; t < tss_count ; t ++) {
			vector < Interval < bool > > ann_in_bin;
			Itree[t].findOverlapping(wfrom, wto, ann_in_bin);
			n_annotation += (ann_in_bin.size() > 0);
		}

		fdo << wfrom << " " << wto << " " << n_annotation << endl;
	}

	fdo.close();
}
