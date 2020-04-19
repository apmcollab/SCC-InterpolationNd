/*
 * LegendreApprox3d.h
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
using namespace std;


#ifndef _LegendreApprox3d_
#define _LegendreApprox3d_

class LegendreApprox3d
{
public:

LegendreApprox3d()
{
	initialize();
}

LegendreApprox3d(const LegendreApprox3d& P)
{
	initialize(P);
}

LegendreApprox3d(int degreeX, int degreeY, int degreeZ)
{
	initialize(degreeX, degreeY, degreeZ);
}

virtual ~LegendreApprox3d(){};

void initialize()
{
	degreeX =  -1;
	degreeY =  -1;
	degreeZ =  -1;

	Ainv.initialize();
	coeff.clear();
}

void initialize(const LegendreApprox3d& P)
{
    degreeX =  P.degreeX;
    degreeY =  P.degreeY;
    degreeZ =  P.degreeZ;

    Ainv.initialize(P.Ainv);

    coeff  = P.coeff;

    legendreValues_X = P.legendreValues_X;
    legendreValues_Y = P.legendreValues_Y;
    legendreValues_Z = P.legendreValues_Z;

    legendreDvalues_X = P.legendreDvalues_X;
    legendreDvalues_Y = P.legendreDvalues_Y;
    legendreDvalues_Z = P.legendreDvalues_Z;

    legendreEval_X   = P.legendreEval_X;
    legendreEval_Y   = P.legendreEval_Y;
    legendreEval_Z   = P.legendreEval_Z;
}

void initialize(int degreeX, int degreeY,int degreeZ)
{
    this->degreeX = degreeX;
    this->degreeY = degreeY;
    this->degreeZ = degreeZ;

	int NX = degreeX+1; // Number of equi-spaced data points in X-direction
	int NY = degreeY+1; // Number of equi-spaced data points in Y-direction
	int NZ = degreeZ+1; // Number of equi-spaced data points in Y-direction


	coeff.resize(NX*NY*NZ,0.0);

	legendreValues_X.resize(NX,0.0);
	legendreValues_Y.resize(NY,0.0);

	legendreDvalues_X.resize(NX,0.0);
	legendreDvalues_Y.resize(NY,0.0);

	legendreDvalues_Z.resize(NZ,0.0);
	legendreDvalues_Z.resize(NZ,0.0);

    // Create matrix for that maps function values at grid
    // points to coefficients of the interpolation Legendre
    // approximation. This construction uses a mesh size
    // of 1. Scaling to the appropriate mesh size is implemented
    // in the evaluation procedure.

	legendreEval_X.initialize(0,(double)degreeX,degreeX);
	legendreEval_Y.initialize(0,(double)degreeY,degreeY);
	legendreEval_Z.initialize(0,(double)degreeZ,degreeZ);

    //
    //  Construct matrix coefficients for interpolation at
    //  equispaced data points
    //

    long N = NX*NY*NZ;


    SCC::LapackMatrix A(N,N);
    Ainv.initialize(N,N);


    double xI; double yJ; double zK;

	for(int i = 0; i < NX; i++)
	{
	xI = (double)i;
	legendreEval_X.evaluate(xI,legendreValues_X);
	for(int j = 0; j < NY; j++)
	{
	yJ = (double)j;
	legendreEval_Y.evaluate(yJ,legendreValues_Y);
	for(int k = 0; k < NZ; k++)
	{
	zK = (double)k;
	legendreEval_Z.evaluate(zK,legendreValues_Z);

		for(int p = 0; p < NX; p++)
		{
		for(int q = 0; q < NY; q++)
		{
		for(int r = 0; r < NZ; r++)
		{
			A(k + j*NZ + i*NZ*NY,r + q*NZ + p*NZ*NY) = legendreValues_X[p]*legendreValues_Y[q]*legendreValues_Z[r];
		}}}
	}}}

    if((degreeX == 0)&&(degreeY == 0)&&(degreeZ == 0)) {Ainv(0,0) = 1.0/A(0,0);}
    else
    {
    computeInverse(A, Ainv);
    }
}

double evaluate(double x, double xMin, double xMax,
                double y, double yMin, double yMax,
                double z, double zMin, double zMax,
                vector<double>& F)
{
    assert((int)F.size() == (degreeX+1)*(degreeY+1)*(degreeZ+1));

    if((degreeX == 0)&&(degreeY == 0)&&(degreeZ == 0)){return F[0];}

    // Replace matrix class product with inline invocation of dgemv

    // coeff = Ainv*F;

    // Create coefficients of the approximation

    long rows      = (degreeX+1)*(degreeY+1)*(degreeZ+1);
    long cols      = rows;
    char TRANS     = 'N';
    double ALPHA   = 1.0;
    double BETA    = 0.0;
    long INCX      = 1;
    long INCY      = 1;

    dgemv_(&TRANS,&rows,&cols,&ALPHA,Ainv.getDataPointer(),&rows,&F[0],&INCX,&BETA,&coeff[0],&INCY);

	// Evaluate interpolant at relative locations in
	// region [0, degreeX]x[0,degreeY]

    double xI; double yJ; double zK;

	xI = (x-xMin)*(degreeX)/(xMax-xMin);
	legendreEval_X.evaluate(xI,legendreValues_X);

	yJ = (y-yMin)*(degreeY)/(yMax-yMin);
	legendreEval_Y.evaluate(yJ,legendreValues_Y);

    zK    = (z-zMin)*(degreeZ)/(zMax-zMin);
	legendreEval_Z.evaluate(zK,legendreValues_Z);

	int NX = degreeX+1;
	int NY = degreeY+1;
	int NZ = degreeZ+1;

	double appVal = 0.0;

	for(int p = 0; p < NX; p++)
	{
	for(int q = 0; q < NY; q++)
	{
	for(int r = 0; r < NZ; r++)
	{
	appVal +=  coeff[r + q*NZ + p*NZ*NY]*legendreValues_X[p]*legendreValues_Y[q]*legendreValues_Z[r];
	}}}


	return appVal;
}

void evaluateDerivative(double x, double xMin, double xMax,
                        double y, double yMin, double yMax,
                        double z, double zMin, double zMax,
                        vector<double>& F, vector<double>& dFvalues)
{
    assert((int)F.size() == (degreeX+1)*(degreeY+1)*(degreeZ+1));

    dFvalues.resize(3,0.0);

    if((degreeX == 0)&&(degreeY == 0)&&(degreeZ == 0))
    {dFvalues[0]= 0.0; dFvalues[1] = 0.0; dFvalues[2] = 0.0; return;}

    // Replace matrix class product with inline invocation of dgemv

    // coeff = Ainv*F;

    // Create coefficients of the approximation

    long rows      = (degreeX+1)*(degreeY+1)*(degreeZ+1);
    long cols      = rows;
    char TRANS     = 'N';
    double ALPHA   = 1.0;
    double BETA    = 0.0;
    long INCX      = 1;
    long INCY      = 1;

    dgemv_(&TRANS,&rows,&cols,&ALPHA,Ainv.getDataPointer(),&rows,&F[0],&INCX,&BETA,&coeff[0],&INCY);


    double xI; double yJ; double zK;

    xI = (x-xMin)*(degreeX)/(xMax-xMin);

	legendreEval_X.evaluate(xI,legendreValues_X);
	legendreEval_X.evaluateDerivatives(xI,legendreDvalues_X);

	yJ = (y-yMin)*(degreeY)/(yMax-yMin);

	legendreEval_Y.evaluate(yJ,legendreValues_Y);
	legendreEval_Y.evaluateDerivatives(yJ,legendreDvalues_Y);

	zK    = (z-zMin)*(degreeZ)/(zMax-zMin);
	legendreEval_Z.evaluate(zK,legendreValues_Z);
	legendreEval_Z.evaluateDerivatives(zK,legendreDvalues_Z);


	int NX = degreeX+1;
	int NY = degreeY+1;
	int NZ = degreeZ+1;

	double appValX = 0.0;
	double appValY = 0.0;
	double appValZ = 0.0;

	for(int p = 0; p < NX; p++)
	{
	for(int q = 0; q < NY; q++)
	{
	for(int r = 0; r < NZ; r++)
	{
	appValX +=  coeff[r + q*NZ + p*NZ*NY]*legendreDvalues_X[p]*legendreValues_Y[q]*legendreValues_Z[r];
	appValY +=  coeff[r + q*NZ + p*NZ*NY]*legendreValues_X[p]*legendreDvalues_Y[q]*legendreValues_Z[r];
	appValZ +=  coeff[r + q*NZ + p*NZ*NY]*legendreValues_X[p]*legendreValues_Y[q]*legendreDvalues_Z[r];
	}}}

	appValX *= (degreeX/(xMax-xMin));
	appValY *= (degreeY/(yMax-yMin));
	appValZ *= (degreeZ/(zMax-zMin));

	dFvalues[0] = appValX;
	dFvalues[1] = appValY;
	dFvalues[2] = appValZ;
}


void evaluate(double x, double xMin, double xMax,
                        double y, double yMin, double yMax,
                        double z, double zMin, double zMax,
                        vector<double>& F, double& Fval, vector<double>& dFvalues)
{
    assert((int)F.size() == (degreeX+1)*(degreeY+1)*(degreeZ+1));

    dFvalues.resize(3,0.0);

    if((degreeX == 0)&&(degreeY == 0)&&(degreeZ == 0))
    {Fval = 0.0; dFvalues[0]= 0.0; dFvalues[1] = 0.0; dFvalues[2] = 0.0; return;}

    // Create coefficients of the approximation

    coeff = Ainv*F;

    double xI; double yJ; double zK;

    xI = (x-xMin)*(degreeX)/(xMax-xMin);

	legendreEval_X.evaluate(xI,legendreValues_X);
	legendreEval_X.evaluateDerivatives(xI,legendreDvalues_X);

	yJ = (y-yMin)*(degreeY)/(yMax-yMin);

	legendreEval_Y.evaluate(yJ,legendreValues_Y);
	legendreEval_Y.evaluateDerivatives(yJ,legendreDvalues_Y);

	zK    = (z-zMin)*(degreeZ)/(zMax-zMin);
	legendreEval_Z.evaluate(zK,legendreValues_Z);
	legendreEval_Z.evaluateDerivatives(zK,legendreDvalues_Z);

	int NX = degreeX+1;
	int NY = degreeY+1;
	int NZ = degreeZ+1;

    double appVal  = 0.0;
	double appValX = 0.0;
	double appValY = 0.0;
	double appValZ = 0.0;

	for(int p = 0; p < NX; p++)
	{
	for(int q = 0; q < NY; q++)
	{
	for(int r = 0; r < NZ; r++)
	{
	appVal  +=  coeff[r + q*NZ + p*NZ*NY]*legendreValues_X[p]*legendreValues_Y[q]*legendreValues_Z[r];
	appValX +=  coeff[r + q*NZ + p*NZ*NY]*legendreDvalues_X[p]*legendreValues_Y[q]*legendreValues_Z[r];
	appValY +=  coeff[r + q*NZ + p*NZ*NY]*legendreValues_X[p]*legendreDvalues_Y[q]*legendreValues_Z[r];
	appValZ +=  coeff[r + q*NZ + p*NZ*NY]*legendreValues_X[p]*legendreValues_Y[q]*legendreDvalues_Z[r];
	}}}

	appValX *= (degreeX/(xMax-xMin));
	appValY *= (degreeY/(yMax-yMin));
	appValZ *= (degreeZ/(zMax-zMin));

    Fval        = appVal;

	dFvalues[0] = appValX;
	dFvalues[1] = appValY;
	dFvalues[2] = appValZ;

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

    vector <long >   IPIV(N);
    long* IPIVptr = &IPIV[0];

    char  EQED;

    vector<double>   R(N);
    double* Rptr  = &R[0];

    vector<double>    C(N);
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
    vector<double>  FERR(NRHS);
    vector<double>  BERR(NRHS);


    vector<double>   WORK(4*N);
    double* WORKptr = &WORK[0];

    vector<long>       IWORK(N);
    long* IWORKptr  = &IWORK[0];

    long   INFO = 0;


    dgesvx_(&FACT, &TRANS, &N, &NRHS, Aptr, &LDA, AFptr, &LDAF, IPIVptr,
	&EQED, Rptr, Cptr, Bptr,&LDB, Xptr, &LDX, &RCOND,
    &FERR[0], &BERR[0], WORKptr, IWORKptr, &INFO);


    if(INFO != 0)
    {
        cerr << "dgesvx  Failed : INFO = " << INFO  << endl;
        exit((int)1);
    }

        //cout << "RCOND " << RCOND << endl;
	}


	int degreeX;
	int degreeY;
	int degreeZ;


    SCC::LapackMatrix    Ainv;
    vector<double>       coeff;

	vector<double>      legendreValues_X;
	vector<double>      legendreValues_Y;
	vector<double>      legendreValues_Z;

	vector<double>     legendreDvalues_X;
	vector<double>     legendreDvalues_Y;
	vector<double>     legendreDvalues_Z;

    LegendrePolyEvaluator legendreEval_X;
    LegendrePolyEvaluator legendreEval_Y;
    LegendrePolyEvaluator legendreEval_Z;
};




#endif /* _LegendreApprox3d_ */
