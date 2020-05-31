/*
 * LegendreApprox1d.h
 *
 *  Created on: Jul 22, 2018
 *      Author: anderson
 */

/*
#############################################################################
#
# Copyright  2018 Chris Anderson
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# For a copy of the GNU General Public License see
# <http://www.gnu.org/licenses/>.
#
#############################################################################
*/

#include "LapackInterface/SCC_LapackMatrix.h"
#include "LapackInterface/SCC_LapackHeaders.h"
#include "LegendrePolyEvaluator.h"

#include <vector>

#ifndef LEGENDRE_APPROX_1D_
#define LEGENDRE_APPROX_1D_

class LegendreApprox1d
{
public:

LegendreApprox1d()
{
	initialize();
}

LegendreApprox1d(const LegendreApprox1d& P)
{
	initialize(P);
}

LegendreApprox1d(int degreeX)
{
	initialize(degreeX);
}

virtual ~LegendreApprox1d(){};

void initialize()
{
	degreeX =  -1;

	Ainv.initialize();
	coeff.clear();
}

void initialize(const LegendreApprox1d& P)
{
    degreeX =  P.degreeX;

    Ainv.initialize(P.Ainv);

    coeff  = P.coeff;

    legendreValues_X   = P.legendreValues_X;
    legendreDvalues_X  = P.legendreDvalues_X;
    legendreEval_X     = P.legendreEval_X;

}

void initialize(int degreeX)
{
    this->degreeX = degreeX;

	int NX = degreeX+1; // Number of equi-spaced data points in X-direction

	coeff.resize(NX,0.0);

	legendreValues_X.resize(NX,0.0);
	legendreDvalues_X.resize(NX,0.0);

    // Create matrix for that maps function values at grid
    // points to coefficients of the interpolation Legendre
    // approximation. This construction uses a mesh size
    // of 1. Scaling to the appropriate mesh size is implemented
    // in the evaluation procedure.

	legendreEval_X.initialize(0,(double)degreeX,degreeX);

    //
    //  Construct matrix coefficients for interpolation at
    //  equispaced data points
    //

    long N = NX;

    SCC::LapackMatrix A(N,N);
    Ainv.initialize(N,N);

    double xI;

	for(int i = 0; i < NX; i++)
	{
	xI = (double)i;
	legendreEval_X.evaluate(xI,legendreValues_X);

		for(int p = 0; p < NX; p++)
		{
		    A(i,p) = legendreValues_X[p];
		}
	}

    if(degreeX == 0) {Ainv(0,0) = 1.0/A(0,0);}
    else
    {
    computeInverse(A, Ainv);
    }
}

double evaluate(double x, double xMin, double xMax,
                std::vector<double>& F)
{
    assert((int)F.size() == (degreeX+1));

    if(degreeX == 0){return F[0];}

    // Create coefficients of the approximation

    coeff = Ainv*F;

	// Evaluate interpolant at relative locations in
	// region [0, degreeX]

    double xI;

	xI = (x-xMin)*(degreeX)/(xMax-xMin);
	legendreEval_X.evaluate(xI,legendreValues_X);

	int NX = degreeX+1;

	double appVal = 0.0;

	for(int p = 0; p < NX; p++)
	{
	appVal +=  coeff[p]*legendreValues_X[p];
	}

	return appVal;
}

double evaluateDerivative(double x, double xMin, double xMax,
                        std::vector<double>& F)
{
    assert((int)F.size() == (degreeX+1));

    double appValX = 0.0;

    if(degreeX == 0){return appValX;}

    // Create coefficients of the approximation

    coeff = Ainv*F;

    double xI;

    xI = (x-xMin)*(degreeX)/(xMax-xMin);

	legendreEval_X.evaluate(xI,legendreValues_X);
	legendreEval_X.evaluateDerivatives(xI,legendreDvalues_X);

	int NX = degreeX+1;

	for(int p = 0; p < NX; p++)
	{
	appValX +=  coeff[p]*legendreDvalues_X[p];
	}

	appValX *= (degreeX/(xMax-xMin));

    return appValX;
}


void evaluate(double x, double xMin, double xMax,
              std::vector<double>& F, double& Fval, double& dFval)
{
    assert((int)F.size() == (degreeX+1));

    Fval  = 0.0;
    dFval = 0.0;

    if(degreeX == 0){Fval = 0.0; dFval= 0.0; return;}

       // Create coefficients of the approximation

    coeff = Ainv*F;

    double xI;

    xI = (x-xMin)*(degreeX)/(xMax-xMin);

	legendreEval_X.evaluate(xI,legendreValues_X);
	legendreEval_X.evaluateDerivatives(xI,legendreDvalues_X);

	int NX = degreeX+1;

    double appVal  = 0.0;
	double appValX = 0.0;

	for(int p = 0; p < NX; p++)
	{
	appVal  +=  coeff[p]*legendreValues_X[p];
	appValX +=  coeff[p]*legendreDvalues_X[p];
	}

	appValX *= (degreeX/(xMax-xMin));

    Fval  =  appVal;
	dFval = appValX;

	return;
}

void computeInverse(const SCC::LapackMatrix& Ainput, SCC::LapackMatrix& Ainv)
{
	assert(Ainput.sizeCheck(Ainput.rows,Ainput.cols));

    char FACT  = 'N'; // Equilibrate, then factor
    char TRANS = 'N'; // No transpose
    long N     = Ainput.rows;

    // Allocate temporaries

    SCC::LapackMatrix A(Ainput);
    SCC::LapackMatrix AF(N,N);
    SCC::LapackMatrix B(N,N);

    double* Aptr  =  A.dataPtr;
    double* AFptr = AF.dataPtr;

    long LDA   = N;
    long LDAF  = N;

    std::vector <long >   IPIV(N);
    long* IPIVptr = &IPIV[0];

    char  EQED;

    std::vector<double>   R(N);
    double* Rptr  = &R[0];

    std::vector<double>    C(N);
    double* Cptr  =  &C[0];

    long NRHS     = N;
    B.setToIdentity();
    double* Bptr  =  B.dataPtr;
    long LDB     =   N;

    // Ainv will be overwritten with the solution
    // so no need to declare X separately

    Ainv.initialize(N,N);
    double* Xptr = Ainv.dataPtr;
    long LDX     = N;

    double         RCOND;
    std::vector<double>  FERR(NRHS);
    std::vector<double>  BERR(NRHS);


    std::vector<double>   WORK(4*N);
    double* WORKptr = &WORK[0];

    std::vector<long>       IWORK(N);
    long* IWORKptr  = &IWORK[0];

    long   INFO = 0;


    dgesvx_(&FACT, &TRANS, &N, &NRHS, Aptr, &LDA, AFptr, &LDAF, IPIVptr,
	&EQED, Rptr, Cptr, Bptr,&LDB, Xptr, &LDX, &RCOND,
    &FERR[0], &BERR[0], WORKptr, IWORKptr, &INFO);


    if(INFO != 0)
    {
        std::cerr << "dgesvx  Failed : INFO = " << INFO  << std::endl;
        exit((int)1);
    }

        //cout << "RCOND " << RCOND << std::endl;
	}

	int degreeX;

    SCC::LapackMatrix     Ainv;
    std::vector<double>       coeff;

	std::vector<double>      legendreValues_X;
	std::vector<double>     legendreDvalues_X;

    LegendrePolyEvaluator legendreEval_X;
};




#endif /* _LegendreApprox1d_ */
