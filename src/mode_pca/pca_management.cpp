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

#include "pca_data.h"


void pca_data::resizeData(){
    PCA._xXf.conservativeResize(NoChange,PCA._xXf.cols() + __RESIZE_CHUNK__);
}

void pca_data::finalizeData(int size){
    PCA._xXf.conservativeResize(NoChange,size);
}

void pca_data::imputeData() {
    for (int g = 0; g < data_count ; g ++) {
        double mean = 0.0;
        int c_mean = 0;
        for (int s = 0; s < sample_count ; s ++) {
            if (PCA._xXf(s,g) != bcf_float_missing) {
                mean += PCA._xXf(s,g);
                c_mean ++;
            }
        }
        mean /= c_mean;
        for (int s = 0; s < sample_count ; s ++) if (PCA._xXf(s,g) == bcf_float_missing) PCA._xXf(s,g) = mean;
    }
}

void pca_data::printPCA(string prefix){
    output_file pca(prefix + ".pca");
    output_file pcaStats(prefix + ".pca_stats");
    pca << "SampleID";
    for (int i = 0 ; i < sample_count; i++) pca << " " << sample_id[i];
    pca<< endl;
    IOFormat MF(StreamPrecision, DontAlignCols, " " , "\n",  "", "", "", "");
    for (int r = 0 ; r < PCA._pcs.rows(); r++ ) {
        pca << prefix << "_" << PCA.is_center() << "_" << PCA.is_scale() <<  "_" << PCA._method + "_PC" + stb.str(r+1) + " " << PCA._pcs.row(r).format(MF) << endl;
    }
    
    pcaStats << "sd ";
    copy(PCA._sd.begin(), PCA._sd.end(),ostream_iterator<float>(pcaStats," "));
    pcaStats << endl;
    pcaStats << "prop_var ";
    copy(PCA._prop_of_var.begin(), PCA._prop_of_var.end(),ostream_iterator<float>(pcaStats," "));
    pcaStats << endl;
    pcaStats << "cumm_prop ";
    copy(PCA._cum_prop.begin(), PCA._cum_prop.end(),ostream_iterator<float>(pcaStats," "));
    pcaStats << endl;
    pcaStats << "#Kaiser criterion: PC #" << PCA._kaiser << endl;
    pcaStats << "#Thresh995 criterion: PC #" << PCA._thresh995 << endl;
}
