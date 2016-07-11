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

#ifndef _GENOMIC_REGION_H
#define _GENOMIC_REGION_H

#define POS_MIN 0000000000
#define POS_MAX 1000000000

#include <string>
#include <exception>

using namespace std;

class genomic_region {
public:
	string chr;
	unsigned int start;
	unsigned int end;

	genomic_region() {
		chr = "NA";
		start = POS_MIN;
		end = POS_MAX;
	}

	genomic_region(string _chr, unsigned int _start, unsigned int _end) {
		chr = _chr;
		start = _start;
		end = _end;
		assert(start <= end);
	}

	genomic_region(genomic_region & gr, unsigned int window) {
		chr = gr.chr;
		if (window > gr.start) start = 0;
		else start = gr.start - window;
		end = gr.end + window;
		assert(start <= end);
	}

	bool isSet () {
		return (chr != "NA");
	}

	string get() {
		ostringstream s2( stringstream::out );
		s2 << chr;
		if (start != POS_MIN || end != POS_MAX) s2 << ":" << start << "-" << end;
		return s2.str();
	}

	bool set (string _chr, unsigned int _start, unsigned int _end) {
		chr = _chr;
		start = _start;
		end = _end;
		if (start > end) return false;
		return true;
	}

	bool parse(string str) {
		size_t split_point_chr = str.find_first_of(":");
		size_t split_point_pos = str.find_first_of("-");

		//
		if (split_point_chr == string::npos && split_point_pos == string::npos) {
			chr = str;
		} else if (split_point_chr == string::npos && split_point_pos != string::npos) {
			return false;
		} else if (split_point_chr != string::npos && split_point_pos == string::npos) {
			chr = str.substr(0, split_point_chr);
			try {
				start = std::stoi(str.substr(split_point_chr+1, string::npos));
				end = std::stoi(str.substr(split_point_chr+1, string::npos));
			} catch (const std::exception & e) {
				return false;
			}
		} else {
			chr = str.substr(0, split_point_chr);
			try {
				start = std::stoi(str.substr(split_point_chr+1, split_point_pos));
				end = std::stoi(str.substr(split_point_pos+1, string::npos));
			} catch (const std::exception & e) {
				return false;
			}
		}
		return true;
	}
};

#endif
