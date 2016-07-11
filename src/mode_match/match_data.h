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

#ifndef _MATCH_DATA_H
#define _MATCH_DATA_H

//INCLUDES
#include "../common/data.h"

class site {
public:
	string chr;
	unsigned int pos;
	unsigned int idx;
	char ref;
	char alt;
	unsigned int cref;
	unsigned int calt;
	unsigned int cdisc;
	unsigned int cqual;
	unsigned int cdel;
	unsigned int ctot;

	site (string _chr, unsigned int _pos, string _ref, string _alt) {
		chr = _chr;
		pos = _pos;
		ref = _ref[0];
		alt = _alt[0];
		cref = 0;
		calt = 0;
		cdisc = 0;
		cqual = 0;
		cdel = 0;
		ctot = 0;
	}
};

class match_data : public data {
public :
	//PARAMETERS
	unsigned int param_min_mapQ;
	unsigned int param_min_baseQ;
	unsigned int param_min_cov;
	float param_min_pval;
	float param_min_gp;
	float param_min_iq;
	bool param_dup_rd;

	//DATA
	vector < string > regions;
	vector < vector < site > > sites;
	vector < vector < vector < bool > > > gen_ref, gen_alt;

	//CONSTRUCTOR/DESTRUCTOR
	match_data();
	~match_data();

	//
	void readGenotypes(string, string);
	void readSequences(string);
	void writeOutput(string);
};

void match_main(vector < string > & );

inline char match_getBase (int code) {
	switch (code) {
	case 1: return 'A';
	case 2: return 'C';
	case 4: return 'G';
	case 8: return 'T';
	case 15: return 'N';
	}
	return -1;
}

inline double match_binomialTest(int x, int n, float p) {
	int y = 0;
	if (p == 0) return (x == 0);
	if (p == 1) return (x == n);
	double relErr = 1 + 1e-07;
	double d = dbinom(x, n, p, 0);
	double m = n * p;
	if (x == m) return 1.0;
	if (x < m) {
		for (int i = (int)ceil (m); i <= n ; i++) y += (dbinom(i, n, p, 0) <= d * relErr);
		return pbinom(x, n, p, 1, 0) + pbinom(n - y, n, p, 0, 0);
	} else {
		for (int i = 0 ; i <= (int)floor(m) ; i++) y += (dbinom(i, n, p, 0) <= d * relErr);
		return pbinom(y - 1, n, p, 1, 0) + pbinom(x - 1, n, p, 0, 0);
	}
}

#endif
