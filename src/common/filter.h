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

#ifndef _DATA_FILTER_H
#define _DATA_FILTER_H

#include <unordered_set>
#include <compressed_io.h>

class filter {
protected:
	unordered_set < string > inclusion_map;
	unordered_set < string > exclusion_map;

public:
	filter() {
	}

	~filter() {
		inclusion_map.clear();
		exclusion_map.clear();
	}

	int readInclusion(string file, bool position = false) {
		unsigned int n_ids = 0;
		string buffer;
		vector < string > tokens;
		input_file fd(file);
		if (fd.fail()) return -1;
		while(getline(fd, buffer, '\n')) {
			stb.split(buffer, tokens);
			if (position && tokens.size() != 2) return -2;
			if (position) inclusion_map.insert(tokens[0] + "_" + tokens[1]);
			else inclusion_map.insert(tokens[0]);
			n_ids++;
		}
		fd.close();
		return n_ids;
	}

	int readExclusion(string file, bool position = false) {
		unsigned int n_ids = 0;
		string buffer;
		vector < string > tokens;
		input_file fd(file);
		if (fd.fail()) return -1;
		while(getline(fd, buffer, '\n')) {
			stb.split(buffer, tokens);
			if (position && tokens.size() != 2) return -2;
			if (position) exclusion_map.insert(tokens[0] + "_" + tokens[1]);
			else exclusion_map.insert(tokens[0]);
			exclusion_map.insert(tokens[0]);
			n_ids++;
		}
		fd.close();
		return n_ids;
	}

	bool check(string id) {
		bool included = ((inclusion_map.size() == 0)?true:inclusion_map.count(id));
		bool excluded = ((exclusion_map.size() == 0)?false:exclusion_map.count(id));
		if (!included || excluded) return false;
		return true;
	}
    
    void addInclusion(string value){
        inclusion_map.insert(value);
    }
    
    void addInclusion(vector < string > & values){
    	for (int e = 0 ; e < values.size() ; e ++) inclusion_map.insert(values[e]);
    }

    void addExclusion(string value){
        exclusion_map.insert(value);
    }
};

#endif
