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

#ifndef _COMPRESSED_IO_H
#define _COMPRESSED_IO_H

//STL INCLUDES
#include <iostream>
#include <sstream>
#include <fstream>

//BOOST INCLUDES
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

class input_file : public boost::iostreams::filtering_istream {
protected:
	ifstream file_descriptor;

public:
	input_file(string filename) {
		if (filename.substr(filename.find_last_of(".") + 1) == "gz") {
			file_descriptor.open(filename.c_str(), ios::in | ios::binary);
			push(boost::iostreams::gzip_decompressor());
		} else if (filename.substr(filename.find_last_of(".") + 1) == "bz2") {
			file_descriptor.open(filename.c_str(), ios::in | ios::binary);
			push(boost::iostreams::bzip2_decompressor());
		} else file_descriptor.open(filename.c_str());
		if (!file_descriptor.fail()) push(file_descriptor);
	}

	~input_file() {
		close();
	}

	bool fail() {
		return file_descriptor.fail();
	}

	void close() {
		if (!file_descriptor.fail()) {
			if (!empty()) reset();
			file_descriptor.close();
		}
	}
};

class output_file : public boost::iostreams::filtering_ostream {
protected:
	ofstream file_descriptor;

public:
	output_file(string filename) {
		if (filename.substr(filename.find_last_of(".") + 1) == "gz") {
			file_descriptor.open(filename.c_str(), ios::out | ios::binary);
			push(boost::iostreams::gzip_compressor());
		} else if (filename.substr(filename.find_last_of(".") + 1) == "bz2") {
			file_descriptor.open(filename.c_str(), ios::out | ios::binary);
			push(boost::iostreams::bzip2_compressor());
		} else file_descriptor.open(filename.c_str());
		if (!file_descriptor.fail()) push(file_descriptor);
	}

	~output_file() {
		close();
	}

	bool fail() {
		return file_descriptor.fail();
	}

	void close() {
		if (!file_descriptor.fail()) {
			if (!empty()) reset();
			file_descriptor.close();
		}
	}
};

#endif
