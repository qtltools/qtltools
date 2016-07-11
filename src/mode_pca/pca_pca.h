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

#ifndef __fastqtl__pca__
#define __fastqtl__pca__

#include <vector>
#include <string>
#include "Eigen/Dense"
#include "../common/data.h"

class Pca {
public:
    std::vector<float>  _x;   // Initial matrix as vector filled by rows.
    unsigned int  _nrows,     // Number of rows in matrix x.
    _ncols;     // Number of cols in matrix x.
    bool  _is_center,         // Whether the variables should be shifted to be zero centered
    _is_scale,          // Whether the variables should be scaled to have unit variance
    _is_corr;           // PCA with correlation matrix, not covariance
    std::string
    _method;            // svd, cor, cov
    std::vector<unsigned int>
    _eliminated_columns;  // Numbers of eliminated columns
    std::vector<float>  _sd,  // Standard deviation of each component
    _prop_of_var,   // Proportion of variance
    _cum_prop,      // Cumulative proportion
    _scores;        // Rotated values
    unsigned int  _kaiser,    // Number of PC according Kaiser criterion
    _thresh995;  // Number of PC according 95% variance threshol
    Eigen::MatrixXf     _xXf; // Initial matrix as Eigen MatrixXf structure
    Eigen::MatrixXf     _pcs;
    //! Initializing values and performing PCA
    /*!
     The main method for performin Principal Component Analysis
     \param  x     Initial data matrix
     \param  nrows Number of matrix rows
     \param  ncols Number of matrix cols
     \param  is_corr   Correlation matrix will be used instead of covariance matrix
     \param  is_center Whether the variables should be shifted to be zero centered
     \param  is_scale  Whether the variables should be scaled to have unit variance
     \result
     0 if everything is Ok
     -1 if there were some errors
     */
    void Calculate(const bool is_corr = true, const bool is_center = true, const bool is_scale = true);
    //! Return number of rows in initial matrix
    /*!
     \result Number of rows in initial matrix
     */
    unsigned int nrows(void);
    //! Return number of cols in initial matrix
    /*!
     \result Number of cols in initial matrix
     */
    unsigned int ncols(void);
    //! If variables are centered
    /*!
     \result
     true - variables are centered
     false - otherwise
     */
    bool is_center(void);
    //! If variables are scaled
    /*!
     \result
     true - variables are scaled
     false - otherwise
     */
    bool is_scale(void);
    //! Method for calculation of principal components
    /*!
     There are different methods used. The most used is SVD.
     But in some cases it may be correlation or covariance matrices.
     If
     \result
     "svd" - PCA with singular value decomposition
     "cor" - PCA with correlation matrix
     "cov" - PCA with covariance matrix
     */
    std::string method(void);
    //! Returns numbers of eliminated columns
    /*!
     If standard deviation of a column is equal to 0, the column shoud be rejected,
     or PCA will fail.
     \result Numbers of eliminated columns, empty vector otherwise
     */
    std::vector<unsigned int> eliminated_columns(void);
    //! Standard deviation of each principal component
    /*!
     \result Vector of standard deviation for each principal component:
     1st element is sd for 1st PC, 2nd - for 2nd PC and so on.
     */
    std::vector<float> sd(void);
    //! Proportion of variance
    /*!
     \result Vector of variances for each component
     */
    std::vector<float> prop_of_var(void);
    //! Cumulative proportion
    /*!
     \result Vector of cumulative proportions for each components
     */
    std::vector<float> cum_prop(void);
    //! Principal component by the Kaiser criterion
    /*!
     Number of the last component with eigenvalue greater than 1.
     \result Number of the first components we should retain defined by the Kaiser criterion
     */
    unsigned int kaiser(void);
    //! 95% threshold
    /*!
     Retain only PC which cumulative proportion is less than 0.95
     \result Number of PCs should be retain with the 95% threshold criterion
     */
    unsigned int thresh995(void);
    //! Rotated values (scores)
    /*!
     Return calculated scores (coordinates in a new space) as vector. Matrix filled by rows.
     \result Vector of scores
     */
    std::vector<float> scores(void);
    //! Class constructor
    Pca(void);
    //! Class destructor
    ~Pca(void);
};

inline void removeRow(Eigen::MatrixXf& matrix, unsigned int rowToRemove){
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();
    
    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);
    
    matrix.conservativeResize(numRows,numCols);
}


inline void removeColumn(Eigen::MatrixXf& matrix, unsigned int colToRemove){
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;
    
    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);
    
    matrix.conservativeResize(numRows,numCols);
}

#endif /* defined(__fastqtl__pca__) */
