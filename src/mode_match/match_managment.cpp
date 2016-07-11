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

#include "match_data.h"

match_data::match_data() {
	param_min_mapQ = 10;
	param_min_baseQ = 5;
	param_min_cov = 10;
	param_min_pval = 0.05;
	param_min_gp = 0.99;
	param_min_iq = 0.90;
}

match_data::~match_data() {
	regions.clear();
	sites.clear();
	gen_ref.clear();
	gen_alt.clear();
}
