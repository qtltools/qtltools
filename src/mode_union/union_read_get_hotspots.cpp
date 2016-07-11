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

#include "union_data.h"

void union_data::readHotspots(string fcov) {
    string buffer;
    vector < string > str;
    int idx = 0;
    vrb.title("Reading hotspots in [" + fcov + "]");
    input_file fd (fcov);
    if(fd.fail()) vrb.error("Cannot open file");
    coldspot_u prev;
    //Read hotspots
    while(getline(fd, buffer)) {
        stb.split(buffer, str);
        if (str.size() < 3) vrb.error("Wrong hotspot file format");
        if (prev.chr != "" && prev.chr != str[0]){
        	coldspot_u *pCs = new coldspot_u(prev.chr,prev.end + 1,INT_MAX,idx,"CS");
            all_coldspots_p.push_back(pCs);
            int sb = (prev.end + 1) / bin_size;
            int eb = INT_MAX / bin_size;
            for (int b = sb; b <= eb; b++) coldspot_bins_p[prev.chr][b].push_back(pCs);
            idx++;
            coldspot_count++;
        }
        int s = prev.chr == str[0] ? prev.end + 1 : 0 ;
        int e = atoi(str[1].c_str());
        int sb = s / bin_size;
        int eb = e / bin_size;
        coldspot_u *pCs = new coldspot_u(str[0],s,e,idx,"CS");
        for (int b = sb; b <= eb; b++) coldspot_bins_p[str[0]][b].push_back(pCs);
        all_coldspots_p.push_back(pCs);
        idx++;
        coldspot_count++;
        prev = coldspot_u(str[0],atoi(str[1].c_str())+1,atoi(str[2].c_str()),-1,"NA");
        s = atoi(str[1].c_str())+1;
        e = atoi(str[2].c_str());
        if (e < s) vrb.error("Hotspot end cannot be smaller than start " + buffer);
        if (prev.chr == str[0] && prev.start > s) vrb.error("Hotspots are not sorted at " + buffer);
        sb = s / bin_size;
        eb = e / bin_size;
        pCs = new coldspot_u(str[0],s,e,idx,"HS");
        for (int b = sb; b <= eb; b++) coldspot_bins_p[str[0]][b].push_back(pCs);
        all_coldspots_p.push_back(pCs);
        idx++;
        coldspot_count++;
    }
	coldspot_u *pCs = new coldspot_u(prev.chr,prev.end + 1,INT_MAX,idx,"CS");
    all_coldspots_p.push_back(pCs);
    int sb = (prev.end + 1) / bin_size;
    int eb = INT_MAX / bin_size;
    for (int b = sb; b <= eb; b++) coldspot_bins_p[prev.chr][b].push_back(pCs);
    idx++;
    coldspot_count++;

    //Finalise
    if (!coldspot_count) vrb.error("No coldspots found");
    vrb.bullet(stb.str(coldspot_count) + " coldspots included");
    fd.close();

    //for (int i = 0 ; i < coldspot_count; i++ ) cerr << (*all_coldspots_p[i]);
}


int union_data::getColdspot(string chr, int pos){
    if (coldspot_bins_p.find(chr) != coldspot_bins_p.end()){
        int max = (coldspot_bins_p[chr].rbegin()->second).back()->end;
        if (pos > max) return -1; //after the last hotspot on this chr
        int bin = pos / bin_size;
        if (coldspot_bins_p[chr].find(bin) != coldspot_bins_p[chr].end()){
            for ( int i = 0 ; i < coldspot_bins_p[chr][bin].size(); i ++ ) if (coldspot_bins_p[chr][bin][i]->start <= pos && coldspot_bins_p[chr][bin][i]->end >= pos) return coldspot_bins_p[chr][bin][i]->idx;
            return -2; //in a hotspot
        }else return -3; //in a hotspot BUT SINCE BIN IS 1MB SHOULD NOT HAPPEN
    }else return -4; //no hospot found for this chr
}
