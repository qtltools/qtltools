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

#ifndef _STRING_UTILS_H
#define _STRING_UTILS_H

#include <sstream>
#include <iomanip>
#include <string>
#include <vector>

//using namespace std;

class string_utils {
public:
	string_utils () {};
	~string_utils () {};

	int split(const string & str, vector < string > & tokens, string sep = " 	", unsigned int n_max_tokens = 1000000) {
		tokens.clear();
		if (str == ""){
			tokens.push_back("");
			return tokens.size();
		}
		string::size_type p_last = str.find_first_not_of(sep, 0);
		string::size_type p_curr = str.find_first_of(sep, p_last);
		while ((string::npos != p_curr || string::npos != p_last) && tokens.size() < n_max_tokens) {
			tokens.push_back(str.substr(p_last, p_curr - p_last));
			p_last = str.find_first_not_of(sep, p_curr);
			p_curr = str.find_first_of(sep, p_last);
		}
		if (tokens.back()[tokens.back().size()-1] == '\r') tokens.back() = tokens.back().substr(0, tokens.back().size()-1);
		return tokens.size();
	}

	bool numeric(string & str) {
		float n;
		std::istringstream in(str);
		if (!(in >> n)) return false;
		return in.rdbuf()->in_avail() == 0;
	}

	template < class T >
	string str(T n, int prec = -1) {
		ostringstream ss( stringstream::out );
		if (prec >= 0) { ss << setiosflags( ios::fixed ); ss.precision(prec); }
		ss << n;
		return ss.str();
	}

	template < class T >
	string str(vector < T > & v, int prec = -1) {
		ostringstream ss( stringstream::out );
		if (prec >= 0) { ss << setiosflags( ios::fixed ); ss.precision(prec); }
		for (int e = 0 ; e < v.size() ; e ++) ss << (e>0?" ":"") << v[e] ;
		return ss.str();
	}
};

#endif
