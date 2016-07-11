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

#include "quan_data.h"


void quan_data::read_Sample_Names(vector <string> &names){
    if (names.size() == 0){
        samples = bams;
    }else if (names.size()==1){
        input_file fd(names[0]);
        if (fd.fail()){
            //Assuming a single name
            fd.close();
            samples.push_back(names[0]);
        }else{
            //Assuming a file with names
            string buffer;
            vector < string > str;
            while(getline(fd, buffer)) {
                stb.split(buffer, str);
                if (str.size()!=1) vrb.error("Expecting a single sample name per line in [" + names[0] + "] but got: " + buffer);
                samples.push_back(str[0]);
            }
            fd.close();
        }
    }else samples = names;
    if (samples.size() != bams.size()) vrb.error("Sample names does not match with BAM files!");
}
