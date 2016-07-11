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

#ifndef _RANDOM_NUMBER_H
#define _RANDOM_NUMBER_H

#include <cfloat>
#include <cstdint>
#include <random>

class random_number_generator {
protected:
	unsigned int seed;
	std::mt19937 randomEngine;
	std::uniform_int_distribution < unsigned int > uniformDistributionInt;
	std::uniform_real_distribution < double > uniformDistributionDouble;

public:

	random_number_generator(unsigned int seed = 15052011) : randomEngine(seed), uniformDistributionInt(0, 32768), uniformDistributionDouble(0, 1.0) {
	}

	~random_number_generator(){
	}

	void setSeed(unsigned int _seed) {
		seed = _seed;
		randomEngine.seed(seed);
	}

	unsigned int getSeed() {
		return seed;
	}

	std::mt19937 & getEngine() {
		return randomEngine;
	}

	unsigned int getInt(unsigned int imin, unsigned int imax) {
		return uniformDistributionInt(randomEngine, std::uniform_int_distribution < unsigned int > {imin, imax}.param());
	}

	unsigned int getInt(unsigned int isize) {
		return getInt(0, isize - 1);
	}

	double getDouble(double fmin, double fmax) {
		return uniformDistributionDouble(randomEngine, std::uniform_real_distribution < double > {fmin, fmax}.param());
	}

	double getDouble() {
		return getDouble(0.0, 1.0);
	}
};

#endif
