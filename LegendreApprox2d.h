/*
 * LegendreApprox2d.h
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

#ifndef  LEGENDRE_APPROX_2D_
#define  LEGENDRE_APPROX_2D_

class LegendreApprox2d
{
public:

LegendreApprox2d()
{
	initialize();
}

LegendreApprox2d(const LegendreApprox2d& P)
{
	initialize(P);
}

LegendreApprox2d(int degreeX, int degreeY)
{
	initialize(degreeX, degreeY);
}

virtual ~LegendreApprox2d(){};

void initialize()
{
	degreeX =  -1;
	degreeY =  -1;

	Ainv.initialize();
	coeff.clear();
}

void initialize(const LegendreApprox2d& P)
{
    degreeX =  P.degreeX;
    degreeY =  P.degreeY;

    Ainv.initialize(P.Ainv);

    coeff  = P.coeff;

    legendreValues_X = P.legendreValues_X;
    legendreValues_Y = P.legendreValues_Y;

    legendreDvalues_X = P.legendreDvalues_X;
    legendreDvalues_Y = P.legendreDvalues_Y;

    legendreEval_X   = P.legendreEval_X;
    legendreEval_Y   = P.legendreEval_Y;
}

void initialize(int degreeX, int degreeY)
{
    this->degreeX = degreeX;
    this->degreeY = degreeY;

	int NX = degreeX+1; // Number of equi-spaced data points in X-direction
	int NY = degreeY+1; // Number of equi-spaced data points in Y-direction

	coeff.resize(NX*NY,0.0);

	legendreValues_X.resize(NX,0.0);
	legendreValues_Y.resize(NY,0.0);

	legendreDvalues_X.resize(NX,0.0);
	legendreDvalues_Y.resize(NY,0.0);

    // Create matrix for that maps function values at grid
    // points to coefficients of the interpolation Legendre
    // approximation. This construction uses a mesh size
    // of 1. Scaling to the appropriate mesh size is implemented
    // in the evaluation procedure.

	legendreEval_X.initialize(0,(double)degreeX,degreeX);
	legendreEval_Y.initialize(0,(double)degreeY,degreeY);

    //
    //  Construct matrix coefficients for interpolation at
    //  equispaced data points
    //

    long N = NX*NY;

    SCC::LapackMatrix A(N,N);
    Ainv.initialize(N,N);

    double xI; double yJ;

	for(int i = 0; i < NX; i++)
	{
	xI = (double)i;
	legendreEval_X.evaluate(xI,legendreValues_X);
	for(int j = 0; j < NY; j++)
	{
	yJ = (double)j;
	legendreEval_Y.evaluate(yJ,legendreValues_Y);

		for(int p = 0; p < NX; p++)
		{
			for(int q = 0; q < NY; q++)
			{
				A(j + i*NY,q + p*NY) = legendreValues_X[p]*legendreValues_Y[q];
		}}
	}}

    if((degreeX == 0)&&(degreeY == 0)) {Ainv(0,0) = 1.0/A(0,0);}
    else
    {
    computeInverse(A, Ainv);
    }
}

double evaluate(double x, double xMin, double xMax,
                double y, double yMin, double yMax,
                std::vector<double>& F)
{
    assert((int)F.size() == (degreeX+1)*(degreeY+1));

    if((degreeX == 0)&&(degreeY == 0)){return F[0];}

    // Create coefficients of the approximation

    coeff = Ainv*F;

	// Evaluate interpolant at relative locations in
	// region [0, degreeX]x[0,degreeY]

    double xI; double yJ;

	xI = (x-xMin)*(degreeX)/(xMax-xMin);
	legendreEval_X.evaluate(xI,legendreValues_X);

	yJ = (y-yMin)*(degreeY)/(yMax-yMin);
	legendreEval_Y.evaluate(yJ,legendreValues_Y);

	int NX = degreeX+1;
	int NY = degreeY+1;

	double appVal = 0.0;

	for(int p = 0; p < NX; p++)
	{
	for(int q = 0; q < NY; q++)
	{
	appVal +=  coeff[q + p*NY]*legendreValues_X[p]*legendreValues_Y[q];
	}}

	return appVal;
}

void evaluateDerivative(double x, double xMin, double xMax,
                        double y, double yMin, double yMax,
                        std::vector<double>& F, std::vector<double>& dFvalues)
{
    assert((int)F.size() == (degreeX+1)*(degreeY+1));

    dFvalues.resize(2,0.0);

    if((degreeX == 0)&&(degreeY == 0)){dFvalues[0]= 0.0; dFvalues[1] = 0.0; return;}

    // Create coefficients of the approximation

    coeff = Ainv*F;

    double xI; double yJ;

    xI = (x-xMin)*(degreeX)/(xMax-xMin);

	legendreEval_X.evaluate(xI,legendreValues_X);
	legendreEval_X.evaluateDerivatives(xI,legendreDvalues_X);

	yJ = (y-yMin)*(degreeY)/(yMax-yMin);

	legendreEval_Y.evaluate(yJ,legendreValues_Y);
	legendreEval_Y.evaluateDerivatives(yJ,legendreDvalues_Y);

	int NX = degreeX+1;
	int NY = degreeY+1;

	double appValX = 0.0;
	double appValY = 0.0;

	for(int p = 0; p < NX; p++)
	{
	for(int q = 0; q < NY; q++)
	{
	appValX +=  coeff[q + p*NY]*legendreDvalues_X[p]*legendreValues_Y[q];
	appValY +=  coeff[q + p*NY]*legendreValues_X[p]*legendreDvalues_Y[q];
	}}

	appValX *= (degreeX/(xMax-xMin));
	appValY *= (degreeY/(yMax-yMin));

	dFvalues[0] = appValX;
	dFvalues[1] = appValY;
}


void evaluate(double x, double xMin, double xMax,
                        double y, double yMin, double yMax,
                        std::vector<double>& F, double& Fval, std::vector<double>& dFvalues)
{
    assert((int)F.size() == (degreeX+1)*(degreeY+1));

    dFvalues.resize(2,0.0);

    if((degreeX == 0)&&(degreeY == 0)){Fval = 0.0; dFvalues[0]= 0.0; dFvalues[1] = 0.0; return;}

       // Create coefficients of the approximation

    coeff = Ainv*F;

    double xI; double yJ;

    xI = (x-xMin)*(degreeX)/(xMax-xMin);

	legendreEval_X.evaluate(xI,legendreValues_X);
	legendreEval_X.evaluateDerivatives(xI,legendreDvalues_X);

	yJ = (y-yMin)*(degreeY)/(yMax-yMin);

	legendreEval_Y.evaluate(yJ,legendreValues_Y);
	legendreEval_Y.evaluateDerivatives(yJ,legendreDvalues_Y);

	int NX = degreeX+1;
	int NY = degreeY+1;

    double appVal  = 0.0;
	double appValX = 0.0;
	double appValY = 0.0;

	for(int p = 0; p < NX; p++)
	{
	for(int q = 0; q < NY; q++)
	{
	appVal  +=  coeff[q + p*NY]*legendreValues_X[p]*legendreValues_Y[q];

	appValX +=  coeff[q + p*NY]*legendreDvalues_X[p]*legendreValues_Y[q];

	appValY +=  coeff[q + p*NY]*legendreValues_X[p]*legendreDvalues_Y[q];
	}}

	appValX *= (degreeX/(xMax-xMin));
	appValY *= (degreeY/(yMax-yMin));

    Fval        =  appVal;
	dFvalues[0] = appValX;
	dFvalues[1] = appValY;
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

        //cout << "RCOND " << RCOND << endl;
	}


	int degreeX;
	int degreeY;

    SCC::LapackMatrix     Ainv;
    std::vector<double>   coeff;

	std::vector<double>   legendreValues_X;
	std::vector<double>   legendreValues_Y;

	std::vector<double>   legendreDvalues_X;
	std::vector<double>   legendreDvalues_Y;

    LegendrePolyEvaluator legendreEval_X;
    LegendrePolyEvaluator legendreEval_Y;
};




#endif /* _LegendreApprox2d_ */
