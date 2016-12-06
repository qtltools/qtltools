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

#ifndef _OLIVIER_TOOLS_H
#define _OLIVIER_TOOLS_H

//INCLUDE STANDARD TEMPLATE LIBRARY USEFULL STUFFS (STL)
#include <vector>
#include <list>
#include <queue>
#include <stack>
#include <bitset>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <string>
#include <exception>
#include <cassert>
#include <limits>

//INCLUDE BOOST USEFULL STUFFS (BOOST)
#include <boost/program_options.hpp>

//INCLUDE HTS LIBRARY
#include <htslib/hts.h>
#include <htslib/kseq.h>
#include <htslib/sam.h>
extern "C" {
	#include <htslib/vcf_sweep.h>
	#include <htslib/synced_bcf_reader.h>
	#include <htslib/vcf.h>
	#include <htslib/vcfutils.h>
}

//INCLUDE RMATH LIBRARY
#define MATHLIB_STANDALONE
#include <Rmath.h>

//INCLUDES BASE STUFFS
#include "genomic_region.h"
#include "interval_tree.h"
#include "compressed_io.h"
#include "random_number.h"
#include "ranker.h"
#include "residualizer.h"
#include "pca.h"
#include "full_linear_regression.h"
#include <basic_stats.h>
#include <string_utils.h>
#include <timer.h>
#include <verbose.h>

//MAKE SOME TOOL FULLY ACCESSIBLE THROUGHOUT THE SOFTWARE
#ifdef _DECLARE_TOOLBOX_HERE
	random_number_generator rng;	//Random number generator
	string_utils stb;				//Utils for string manipulation
	verbose vrb;					//Verbose
#else
	extern random_number_generator rng;
	extern string_utils stb;
	extern verbose vrb;
#endif

#endif

