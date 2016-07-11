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

#ifndef full_linear_regression_h
#define full_linear_regression_h

#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#ifndef MATHLIB_STANDALONE
#define MATHLIB_STANDALONE
#endif
#include <Rmath.h>

using namespace std;


class linReg{
public:
    double pval;
    double corr;
    double r;
    double beta;
    double yIntercept;
    double se;
    vector < float > residuals,X,Y;
    linReg(){pval=1;corr=0;beta=0; yIntercept = 0; residuals = vector < float> (0);X = vector < float> (0);Y = vector < float> (0); se = 0;r=0;}
    linReg(double p ,double c, double b, double y , vector <float> &res){pval=p;corr=c;beta=b; yIntercept = y; residuals = res; r = beta < 0.0 ? sqrt(corr) * -1.0 : sqrt(corr);}
    friend ostream& operator<<(ostream& out, const linReg& l){
        out << "X";
        for (int i =0 ; i < l.X.size(); i++) out << "\t" << l.X[i];
        out << "\n";
        out << "Y";
        for (int i =0 ; i < l.Y.size(); i++) out << "\t" << l.Y[i];
        out << "\n";
        out << "Residuals";
        for (int i =0 ; i < l.residuals.size(); i++) out << "\t" << l.residuals[i];
        out << "\n";
        out << "R2\t" << l.corr << "\n";
        out << "R\t" << l.r << "\n";
        out << "Beta\t" << l.beta << "\n";
        out << "SE\t" << l.se << "\n";
        out << "Y-int\t" << l.yIntercept << "\n";
        out << "Pval\t" << l.pval;
        return out;
    }
    
    linReg( vector <float > & x, vector <float > &y){
        assert(x.size() == y.size());
        X = x;
        Y = y;
        double dataSize = (double) x.size();
        residuals = vector < float > ( (unsigned int) dataSize, 0.0);
        double sum_x = 0.0;     //sum of x values
        double sum_y = 0.0;     //sum of y values
        double sum_xy = 0.0;    //sum of x * y
        double sum_xx = 0.0;    //sum of x^2
        double sum_res = 0.0;   //sum of squared residue
        double sum_Yres = 0.0; //sum of squared of the discrepancies
        double AVGy = 0.0;     //mean of y
        double AVGx = 0.0;     //mean of x
        
        for (int i = 0 ; i < dataSize; i++){
            sum_x += X[i];
            sum_y += Y[i];
            sum_xy += X[i] * Y[i];
            sum_xx += X[i] * X[i];
        }
        
        AVGy = sum_y / dataSize;
        AVGx = sum_x / dataSize;
        
        beta = (dataSize * sum_xy - sum_x * sum_y) / (dataSize * sum_xx - sum_x*sum_x);
        yIntercept = AVGy - beta * AVGx;
        
        for (int i = 0 ; i < dataSize; i++){
            residuals[i] = Y[i] - (X[i] * beta + yIntercept);
            sum_Yres += residuals[i] * residuals[i];
            sum_res += pow(Y[i] - AVGy,2);
        }
        corr = (sum_res - sum_Yres) / sum_res;
        r = beta < 0 ? sqrt(corr) * -1.0 : sqrt(corr);
        se = sqrt(sum_Yres / (dataSize-2));
        pval = pf((dataSize-2) * corr / (1 - corr), 1, (dataSize-2), 0, 0);
    }
};

#undef MATHLIB_STANDALONE
#endif /* full_linear_regression_h */
