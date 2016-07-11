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

#include "pca_pca.h"
#include <iostream>
#include <iterator>
#include "Eigen/SVD"


using namespace std;
using namespace Eigen;

void Pca::Calculate(const bool is_corr, const bool is_center, const bool is_scale) {
    _ncols = _xXf.cols();
    _nrows = _xXf.rows();
    _is_corr = is_corr;
    _is_center = is_center;
    _is_scale = is_scale;
    
    vrb.title("Calculating PCA with a matrix " + stb.str(_nrows) + " x "  + stb.str(_ncols));
    
    if ((1 == _ncols) || (1 == _nrows)) {
        vrb.error("Row or column count equals 1!");
    }
    // Mean and standard deviation for each column
    VectorXf mean_vector(_ncols);
    mean_vector = _xXf.colwise().mean();
    VectorXf sd_vector(_ncols);
    unsigned int zero_sd_num = 0;
    float denom = static_cast<float>((_nrows > 1)? _nrows - 1: 1);
    for (unsigned int i = 0; i < _ncols; ++i) {
        VectorXf curr_col  = VectorXf::Constant(_nrows, mean_vector(i)); // mean(x) for column x
        curr_col = _xXf.col(i) - curr_col; // x - mean(x)
        curr_col = curr_col.array().square(); // (x-mean(x))^2
        sd_vector(i) = sqrt((curr_col.sum())/denom);
        if (0 == sd_vector(i)) {
            zero_sd_num++;
        }
    }
    vrb.bullet("Calculated mean and sd");
    // If colums with sd == 0 are too many,
    // don't continue calculations
    if (1 > _ncols-zero_sd_num) {
        vrb.error("No variant data found!");
    }
    
    // Delete columns where sd == 0
    MatrixXf tmp(_nrows, _ncols-zero_sd_num);
    VectorXf tmp_mean_vector(_ncols-zero_sd_num);
    unsigned int curr_col_num = 0;
    for (unsigned int i = 0; i < _ncols; ++i) {
        if (0 != sd_vector(i)) {
            tmp.col(curr_col_num) = _xXf.col(i);
            tmp_mean_vector(curr_col_num) = mean_vector(i);
            curr_col_num++;
        } else {
            _eliminated_columns.push_back(i);
        }
    }
    _ncols -= zero_sd_num;
    _xXf = tmp;
    mean_vector = tmp_mean_vector;
    tmp.resize(0, 0); tmp_mean_vector.resize(0);
    
    vrb.bullet( stb.str(zero_sd_num) + " sd==0 columns removed");
    
    // Shift to zero
    if (true == _is_center) {
        for (unsigned int i = 0; i < _ncols; ++i) {
            _xXf.col(i) -= VectorXf::Constant(_nrows, mean_vector(i));
        }
        vrb.bullet("Variables centered");
    }else vrb.warning("Variables are NOT centered");
    
    // Scale to unit variance
    if ( true == _is_scale) {
        for (unsigned int i = 0; i < _ncols; ++i) {
            _xXf.col(i) /= sqrt(_xXf.col(i).array().square().sum()/denom);
        }
        vrb.bullet("Variables scaled to unit variance");
    }else vrb.warning("Variables are NOT scaled");

    // When _nrows < _ncols then svd will be used.
    // If corr is true and _nrows > _ncols then will be used correlation matrix
    // (TODO): What about covariance?
    if ( (_nrows < _ncols) || (false == _is_corr)) { // Singular Value Decomposition is on
        _method = "svd";
        //if ( (_nrows < _ncols) && true == _is_corr) vrb.bullet("Number of rows is less than number of columns defaulting to SVD");
        //if ( (_nrows >= _ncols) && false == _is_corr) vrb.bullet("Number of rows is greater than number of columns, but SVD is forced (--use-cor if you don't want this)");
        vrb.bullet("Running SVD");
        JacobiSVD<MatrixXf> svd(_xXf, ComputeThinV);
        VectorXf eigen_singular_values = svd.singularValues();
        VectorXf tmp_vec = eigen_singular_values.array().square();
        vrb.bullet("SVD complete");
        float tmp_sum = tmp_vec.sum();
        tmp_vec /= tmp_sum;
        // PC's standard deviation and
        // PC's proportion of variance
        _kaiser = 0;
        unsigned int lim = (_nrows < _ncols)? _nrows : _ncols;
        for (unsigned int i = 0; i < lim; ++i) {
            _sd.push_back(eigen_singular_values(i)/sqrt(denom));
            if (_sd[i] >= 1) {
                _kaiser = i + 1;
            }
            _prop_of_var.push_back(tmp_vec(i));
        }
        vrb.bullet("Calculated sd and var for PCs");
        tmp_vec.resize(0);
        // PC's cumulative proportion
        _thresh995 = 1;
        _cum_prop.push_back(_prop_of_var[0]);
        for (unsigned int i = 1; i < _prop_of_var.size(); ++i) {
            _cum_prop.push_back(_cum_prop[i-1]+_prop_of_var[i]);
            if (_cum_prop[i] < 0.995) {
                _thresh995 = i+1;
            }
        }
        vrb.bullet("Calculated cumulative var for PCs");
        // Scores
        MatrixXf eigen_scores = _xXf * svd.matrixV();
        _pcs = eigen_scores.transpose();
        eigen_scores.resize(0, 0);
        vrb.bullet("Done!");
    } else { // COR OR COV MATRICES ARE HERE
        _method = "cor";
        vrb.bullet("PCA with correlation matrix");
        // Calculate covariance matrix
        MatrixXf eigen_cov; // = MatrixXf::Zero(_ncols, _ncols);
        VectorXf sds;
        // (TODO) Should be weighted cov matrix, even if is_center == false
        eigen_cov = (1.0 /(_nrows/*-1*/)) * _xXf.transpose() * _xXf;
        sds = eigen_cov.diagonal().array().sqrt();
        MatrixXf outer_sds = sds * sds.transpose();
        eigen_cov = eigen_cov.array() / outer_sds.array();
        vrb.bullet("Covariance complete");
        outer_sds.resize(0, 0);
        // ?If data matrix is scaled, covariance matrix is equal to correlation matrix

        EigenSolver<MatrixXf> edc(eigen_cov);
        VectorXf eigen_eigenvalues = edc.eigenvalues().real();
        MatrixXf eigen_eigenvectors = edc.eigenvectors().real();


        // The eigenvalues and eigenvectors are not sorted in any particular order.
        // So, we should sort them
        typedef pair<float, int> eigen_pair;
        vector<eigen_pair> ep;
        for (unsigned int i = 0 ; i < _ncols; ++i) {
            ep.push_back(make_pair(eigen_eigenvalues(i), i));
        }
        sort(ep.begin(), ep.end()); // Ascending order by default
        // Sort them all in descending order
        MatrixXf eigen_eigenvectors_sorted = MatrixXf::Zero(eigen_eigenvectors.rows(), eigen_eigenvectors.cols());
        VectorXf eigen_eigenvalues_sorted = VectorXf::Zero(_ncols);
        int colnum = 0;
        int i = ep.size()-1;
        for (; i > -1; i--) {
            eigen_eigenvalues_sorted(colnum) = ep[i].first;
            eigen_eigenvectors_sorted.col(colnum++) += eigen_eigenvectors.col(ep[i].second);
        }

        // We don't need not sorted arrays anymore
        eigen_eigenvalues.resize(0);
        eigen_eigenvectors.resize(0, 0);
        
        _sd.clear(); _prop_of_var.clear(); _kaiser = 0;
        float tmp_sum = eigen_eigenvalues_sorted.sum();
        for (unsigned int i = 0; i < _ncols; ++i) {
            _sd.push_back(sqrt(eigen_eigenvalues_sorted(i)));
            if (_sd[i] >= 1) {
                _kaiser = i + 1;
            }
            _prop_of_var.push_back(eigen_eigenvalues_sorted(i)/tmp_sum);
        }
        vrb.bullet("Calculated sd and var for PCs");

        // PC's cumulative proportion
        _cum_prop.clear(); _thresh995 = 1;
        _cum_prop.push_back(_prop_of_var[0]);
        for (unsigned int i = 1; i < _prop_of_var.size(); ++i) {
            _cum_prop.push_back(_cum_prop[i-1]+_prop_of_var[i]);
            if (_cum_prop[i] < 0.995) {
                _thresh995 = i+1;
            }
        }
        vrb.bullet("Calculated cumulative var for PCs");
        // Scores for PCA with correlation matrix
        // Scale before calculating new values
        for (unsigned int i = 0; i < _ncols; ++i) {
            _xXf.col(i) /= sds(i);
        }
        sds.resize(0);
        MatrixXf eigen_scores = _xXf * eigen_eigenvectors_sorted;
        _pcs = eigen_scores.transpose();
        eigen_scores.resize(0, 0);
        vrb.bullet("Done!");
    }
}
std::vector<float> Pca::sd(void) { return _sd; };
std::vector<float> Pca::prop_of_var(void) {return _prop_of_var; };
std::vector<float> Pca::cum_prop(void) { return _cum_prop; };
std::vector<float> Pca::scores(void) { return _scores; };
std::vector<unsigned int> Pca::eliminated_columns(void) { return _eliminated_columns; }
string Pca::method(void) { return _method; }
unsigned int Pca::kaiser(void) { return _kaiser; };
unsigned int Pca::thresh995(void) { return _thresh995; };
unsigned int Pca::ncols(void) { return _ncols; }
unsigned int Pca::nrows(void) { return _nrows; }
bool Pca::is_scale(void) {  return _is_scale; }
bool Pca::is_center(void) { return _is_center; }
Pca::Pca(void) {
    _nrows      = 0;
    _ncols      = 0;
    // Variables will be scaled by default
    _is_center  = true;
    _is_scale   = true;  
    // By default will be used singular value decomposition
    _method   = "svd";
    _is_corr  = false;
    
    _kaiser   = 0;
    _thresh995 = 1;
}
Pca::~Pca(void) { 
    _xXf.resize(0, 0);
    _x.clear();
}
