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

bool genrich_data::isSameSignal(unsigned int idx1, unsigned int idx2) {
	if (idx1 == idx2) return true;
	if (genotype_chr[idx1] != genotype_chr[idx2]) return false;
	if (abs(genotype_pos[idx1] - genotype_pos[idx2]) > 1000000) return false;

	double h11 = 0.0, hX1 = 0.0, h1X = 0.0;
	for (int i = 0 ; i < 2 * sample_count ; i ++) {
		if (genotype_haps[idx1][i]) h1X ++;
		if (genotype_haps[idx2][i]) hX1 ++;
		if (genotype_haps[idx1][i] && genotype_haps[idx2][i]) h11 ++;
	}
	h11 /= 2 * sample_count;
	hX1 /= 2 * sample_count;
	h1X /= 2 * sample_count;
	double r2 = (h11 - h1X * hX1) * (h11 - h1X * hX1) / (h1X * (1 - h1X) * hX1 * (1 - hX1));
	//if (!(r2 >= 0 && r2 <= 1.0)) vrb.warning("LD r2 of " + stb.str(r2));
	if (r2 >= threshold_ld) return true;
	else return false;
}

int genrich_data::getDistance(unsigned int chr, int pos) {
	vector < Interval < pair < bool, bool > > > phenotype_pair;
	phenotype_pos[chr].findOverlapping(pos, phenotype_pair);
	assert(phenotype_pair.size() == 1);
	int x1 = abs((int)pos - (int)phenotype_pair[0].start);
	int x2 = abs((int)pos - (int)phenotype_pair[0].stop);
	if (x1 < x2) return ((phenotype_pair[0].value.first)?(-1*x1):x1);
	else return ((phenotype_pair[0].value.second)?(-1*x2):x2);
}

int genrich_data::findCHR (string & chr) {
	unordered_map < string, unsigned int > :: iterator it = chromosome_idx.find(chr);
	if (it == chromosome_idx.end()) return -1;
	else return it->second;
}
