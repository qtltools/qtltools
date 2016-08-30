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

#include "fenrich_data.h"

int fenrich_data::findTSS(string & tss_str) {
	for (int t = 0 ; t < tss_count ; t ++) if (tss_id[t] == tss_str) return t;
	return -1;
}

void fenrich_data::mapAnnotation2QTL() {

	//0. Initilization
	vrb.title("Mapping annotations to TSS ");

	//1. Enumerate chr
	map < string, int > chr2idx;
	for (int t = 0 ; t < tss_count ; t ++) {
		map < string, int >::iterator itC = chr2idx.find(tss_chr[t]);
		if (itC == chr2idx.end()) chr2idx.insert(pair < string, int > (tss_chr[t], chr2idx.size()));
	}
	unsigned int n_chr = chr2idx.size();
	vrb.bullet("#detected chromosomes in TSS data = " + stb.str(n_chr));

	//2. Build chromosomal interval trees
	unsigned int chr_unfound = 0;
	vector < vector < Interval < bool > > > Tvec = vector < vector < Interval < bool > > > (n_chr, vector < Interval < bool > > ());
	for (int a = 0 ; a < ann_count ; a ++) {
		map < string, int >::iterator itC = chr2idx.find(ann_chr[a]);
		if (itC == chr2idx.end()) chr_unfound ++;
		else Tvec[itC->second].push_back(Interval < bool > (ann_start[a], ann_end[a], true));
	}
	vrb.bullet("#annotations NOT mapped to TSS chromosomes = " + stb.str(chr_unfound));
	if ((ann_count - chr_unfound) > 0) vrb.bullet("#annotations mapped to TSS chromosomes = " + stb.str(ann_count - chr_unfound));
	else vrb.error("None of the annotations have been found to be located on the same chromosomes than TSS!");
	vector <  IntervalTree < bool > > Ttree = vector <  IntervalTree < bool > > (n_chr, IntervalTree < bool > ());
	for (int c = 0 ; c < n_chr ; c ++) Ttree[c] = IntervalTree < bool > (Tvec[c]);

	//3. Work out cis-window size in QTL data
	unsigned int max_distance = 0;
	for (int q = 0 ; q < qtl_count ; q ++) if (abs(qtl_pos[q]) > max_distance) max_distance = abs(qtl_pos[q]);
	vrb.bullet("estimated cis-window size from the data = " + stb.str(max_distance));

	//4. Build functional neighborhoods
	basic_stats Rstat;
	R = vector <  IntervalTree < bool > > (tss_count);
	for (int t = 0 ; t < tss_count ; t ++) {
		map < string, int >::iterator itC = chr2idx.find(tss_chr[t]);
		assert(itC != chr2idx.end());
		vector < Interval < bool > > ann_in_cis;
		Ttree[itC->second].findOverlapping(tss_pos[t] - max_distance - 10000, tss_pos[t] + max_distance + 10000, ann_in_cis);
		IntervalStartSorter < bool > intervalStartSorter;
		sort(ann_in_cis.begin(), ann_in_cis.end(), intervalStartSorter);
		Rstat.push(ann_in_cis.size() * 1.0);
		vector < Interval < bool > > Rvec;

		for (int a = 0 ; a < ann_in_cis.size() ; a ++) {
			if (!tss_neg[t]) Rvec.push_back(Interval < bool > (ann_in_cis[a].start - tss_pos[t], ann_in_cis[a].stop - tss_pos[t], true));
			else Rvec.push_back(Interval < bool > (-1 * (ann_in_cis[a].stop - tss_pos[t]), -1 * (ann_in_cis[a].start - tss_pos[t]), true));
			assert(Rvec.back().start <= Rvec.back().stop);
		}
		sort(Rvec.begin(), Rvec.end(), intervalStartSorter);
		R[t] = IntervalTree < bool > (Rvec);
	}
	vrb.bullet("#annotated cis-windows = " + stb.str(Rstat.size()));
	vrb.bullet("#annotations per cis-window = " + stb.str(Rstat.mean(), 2) + " +/- " + stb.str(Rstat.sd(), 2));
}

unsigned int fenrich_data::countOverlaps() {
	unsigned int n_overlaps = 0;
	for (int q = 0 ; q < qtl_count ; q ++) {
		if (R[qtl_order[q]].checkOverlapping(qtl_pos[q])) n_overlaps ++;
	}
	return n_overlaps;
}
