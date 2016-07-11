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

#include "bamstat_data.h"

void bamstat_data::readAnnotationsBED(string fbed) {
        string buffer;
        vector < string > tok;
        vrb.title("Read BED annotations in ["  + fbed  + "]");
        input_file fd (fbed);
        while (getline(fd, buffer)) {
                stb.split(buffer, tok);
                if (tok.size() < 3) vrb.error("Incorrect number of columns in BED file!");
                R.push_back(bamstat_region(tok[0], atoi(tok[1].c_str()), atoi(tok[2].c_str())));
        }
        if (R.size() > 0) vrb.bullet("#annotations = " + stb.str(R.size()));
        else vrb.error("Could not find any annotations in file!");
        fd.close();
}





