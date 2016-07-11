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

#include "trans_data.h"


trans_data::trans_data() {
    sample_count = 0;
    phenotype_count = 0;
    covariate_count = 0;
    beta_ml1 = 0;
    beta_ml2 = 0;
    correlation_threshold = 0;
    start_line=0;
    end_line=0;
}

void trans_data::clear() {
    sample_count = 0;
    sample_id.clear();
    phenotype_count = 0;
    phenotype_val.clear();
    phenotype_id.clear();
    phenotype_chr.clear();
    phenotype_start.clear();
    phenotype_end.clear();
    covariate_count = 0;
    covariate_val.clear();
    covariate_id.clear();
    sample_count = 0;
    phenotype_count = 0;
    covariate_count = 0;
    beta_ml1 = 0;
    beta_ml2 = 0;
    correlation_threshold = 0;
}

trans_data::~trans_data() {
    clear();
}