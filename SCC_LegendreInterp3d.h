#include "GridFunctionNd/SCC_GridFunction3d.h"
#include "DoubleVectorNd/SCC_DoubleVector3d.h"

#include "SCC_ProductLegendrePoly3d.h"

#include "DoubleVectorNd/SCC_DoubleVector2d.h"

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <stdexcept>
#include <array>


#ifndef LEGENDRE_INTERP_3D_
#define LEGENDRE_INTERP_3D_


/*
#############################################################################
#
# Copyright  2015 Chris Anderson
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

namespace SCC
{


class LegendreInterp3d
{
public:

LegendreInterp3d()
{
    initialize();
}

void initialize()
{
    legendreInd_X      = 0;
    legendreInd_Y      = 0;
    legendreInd_Z      = 0;
	interpMatrixInvTranspose.initialize();
    evaluationVector.clear();
    weightVector.clear();
}

LegendreInterp3d(const LegendreInterp3d& LegendreInterp3d)
{
	initialize(LegendreInterp3d);
}

LegendreInterp3d(int interpolationOrder)
{
	initialize(interpolationOrder);
}


void initialize(const LegendreInterp3d& LegendreInterp3d)
{
    legendreInd_X            =  LegendreInterp3d.legendreInd_X;
    legendreInd_Y            =  LegendreInterp3d.legendreInd_Y;
    legendreInd_Z            =  LegendreInterp3d.legendreInd_Z;
	interpMatrixInvTranspose =  LegendreInterp3d.interpMatrixInvTranspose;
}

void initialize(int interpolationOrder)
{
    long legendreIndex = interpolationOrder-1;

    legendreInd_X  = legendreIndex;
    legendreInd_Y  = legendreIndex;
    legendreInd_Z  = legendreIndex;

    createInterpolationMatrixInvTranspose();
}

double evaluateInterpolant(double xPos, double yPos, double zPos, const SCC::GridFunction3d& f)
{
	createNodesAndWeightsData(xPos,yPos,zPos,f);

    double interpolationValue  = 0.0;
    long weightIndex           = 0;
    for(long i = xMinIndex; i <= xMaxIndex; i++)
    {
    for(long j = yMinIndex; j <= yMaxIndex; j++)
    {
    for(long k = zMinIndex; k <= zMaxIndex; k++)
    {
    interpolationValue += f(i,j,k)*weightVector[weightIndex];
    weightIndex++;
    }}}

    return interpolationValue;
}


void createNodesAndWeightsData(double xPos, double yPos, double zPos, array<long,3>& coordPanels,
array<double,3>& coordMin, array<double, 3>& coordMax)
{
	// Determine the range of grid indices associated with the interpolant, a range
	// determined so that the interpolated value is centered as much as possible
	// with respect to the interpolation values.


    long xPanel = coordPanels[0];
    long yPanel = coordPanels[1];
    long zPanel = coordPanels[2];

   	double hx = (coordMax[0] - coordMin[0])/xPanel;
    double hy = (coordMax[1] - coordMin[1])/yPanel;
    double hz = (coordMax[2] - coordMin[2])/zPanel;

    string errMessage;

    if((legendreInd_X > xPanel)
     ||(legendreInd_Y > yPanel)
     ||(legendreInd_Z > zPanel))
    {
    errMessage.assign("\nXXX LegendreInterp3d XXX\n");
    errMessage += "In evaluateInterpolant(...) \n";
	errMessage += "SCC::GridFunction3d argument of insufficient size for\n";
    errMessage += "requested interpolation.  \n";
    throw std::runtime_error(errMessage);
    }

    double xMin   = coordMin[0];
    double yMin   = coordMin[1];
    double zMin   = coordMin[2];

    long    xInterpIndex = round((xPos-xMin)/hx);
    long    yInterpIndex = round((yPos-yMin)/hy);
    long    zInterpIndex = round((zPos-zMin)/hz);

    if(xInterpIndex < 0)      xInterpIndex = 0;
    if(yInterpIndex < 0)      yInterpIndex = 0;
    if(zInterpIndex < 0)      zInterpIndex = 0;
    if(xInterpIndex > xPanel) xInterpIndex = xPanel;
    if(yInterpIndex > yPanel) yInterpIndex = yPanel;
    if(zInterpIndex > zPanel) zInterpIndex = zPanel;

    getIndexBounds(xInterpIndex,legendreInd_X,hx,xPanel,xMin,xPos, xMinIndex, xMaxIndex);
    getIndexBounds(yInterpIndex,legendreInd_Y,hy,yPanel,yMin,yPos, yMinIndex, yMaxIndex);
    getIndexBounds(zInterpIndex,legendreInd_Z,hz,zPanel,zMin,zPos, zMinIndex, zMaxIndex);

    ProductLegendrePoly3d productLegendre(-1.0,1.0,-1.0,1.0,-1.0,1.0);
    long interpSystemSize = (legendreInd_X + 1)*(legendreInd_Y + 1)*(legendreInd_Z + 1);

    evaluationVector.resize(interpSystemSize);
    weightVector.resize(interpSystemSize);

    double xEvalPos = xPos - (xMinIndex*hx + xMin);
    double yEvalPos = yPos - (yMinIndex*hy + yMin);
    double zEvalPos = zPos - (zMinIndex*hz + zMin);

    // Scale evaluation point to [-1.0, 1.0]

    if(legendreInd_X == 0) {xEvalPos /= hx;}
    else           {xEvalPos = -1.0 + (2.0*xEvalPos)/(legendreInd_X*hx);}

    if(legendreInd_Y == 0) {yEvalPos /= hy;}
    else           {yEvalPos = -1.0 + (2.0*yEvalPos)/(legendreInd_Y*hy);}

    if(legendreInd_Z == 0) {zEvalPos /= hz;}
    else           {zEvalPos = -1.0 + (2.0*zEvalPos)/(legendreInd_Z*hz);}

    long prodIndex = 0;
    int pxInd; int pyInd; int pzInd;

    for(pxInd = 0; pxInd <= legendreInd_X; pxInd++)
    {
    for(pyInd = 0; pyInd <= legendreInd_Y; pyInd++)
    {
    for(pzInd = 0; pzInd <= legendreInd_Z; pzInd++)
    {
    	evaluationVector[prodIndex] = productLegendre.evaluate(xEvalPos,pxInd,yEvalPos,pyInd,zEvalPos,pzInd);
    	prodIndex++;
    }}}

    // Create weight std::vector

    for(long i = 0; i < interpSystemSize; i++)
    {
    weightVector[i] = 0.0;
    for(long j = 0; j < interpSystemSize; j++)
    {
    weightVector[i] += interpMatrixInvTranspose(i,j)*evaluationVector[j];
    }}
}

void createNodesAndWeightsData(double xPos, double yPos, double zPos, const SCC::GridFunction3d& f)
{
	// Determine the range of grid indices associated with the interpolant, a range
	// determined so that the interpolated value is centered as much as possible
	// with respect to the interpolation values.

	double hx = f.getHx();
    double hy = f.getHy();
    double hz = f.getHz();

    long xPanel = f.getXpanelCount();
    long yPanel = f.getYpanelCount();
    long zPanel = f.getZpanelCount();

    string errMessage;

    if((legendreInd_X > xPanel)
     ||(legendreInd_Y > yPanel)
     ||(legendreInd_Z > zPanel))
    {
    errMessage.assign("\nXXX LegendreInterp3d XXX\n");
    errMessage += "In evaluateInterpolant(...) \n";
	errMessage += "SCC::GridFunction3d argument of insufficient size for\n";
    errMessage += "requested interpolation.  \n";
    throw runtime_error(errMessage);
    }

    double xMin   = f.getXmin();
    double yMin   = f.getYmin();
    double zMin   = f.getZmin();

    long    xInterpIndex = round((xPos-xMin)/hx);
    long    yInterpIndex = round((yPos-yMin)/hy);
    long    zInterpIndex = round((zPos-zMin)/hz);

    if(xInterpIndex < 0)      xInterpIndex = 0;
    if(yInterpIndex < 0)      yInterpIndex = 0;
    if(zInterpIndex < 0)      zInterpIndex = 0;
    if(xInterpIndex > xPanel) xInterpIndex = xPanel;
    if(yInterpIndex > yPanel) yInterpIndex = yPanel;
    if(zInterpIndex > zPanel) zInterpIndex = zPanel;

    getIndexBounds(xInterpIndex,legendreInd_X,hx,xPanel,xMin,xPos, xMinIndex, xMaxIndex);
    getIndexBounds(yInterpIndex,legendreInd_Y,hy,yPanel,yMin,yPos, yMinIndex, yMaxIndex);
    getIndexBounds(zInterpIndex,legendreInd_Z,hz,zPanel,zMin,zPos, zMinIndex, zMaxIndex);

    ProductLegendrePoly3d productLegendre(-1.0,1.0,-1.0,1.0,-1.0,1.0);
    long interpSystemSize = (legendreInd_X + 1)*(legendreInd_Y + 1)*(legendreInd_Z + 1);

    evaluationVector.resize(interpSystemSize);
    weightVector.resize(interpSystemSize);

    double xEvalPos = xPos - (xMinIndex*hx + xMin);
    double yEvalPos = yPos - (yMinIndex*hy + yMin);
    double zEvalPos = zPos - (zMinIndex*hz + zMin);

    // Scale evaluation point to [-1.0, 1.0]

    if(legendreInd_X == 0) {xEvalPos /= hx;}
    else           {xEvalPos = -1.0 + (2.0*xEvalPos)/(legendreInd_X*hx);}

    if(legendreInd_Y == 0) {yEvalPos /= hy;}
    else           {yEvalPos = -1.0 + (2.0*yEvalPos)/(legendreInd_Y*hy);}

    if(legendreInd_Z == 0) {zEvalPos /= hz;}
    else           {zEvalPos = -1.0 + (2.0*zEvalPos)/(legendreInd_Z*hz);}

    long prodIndex = 0;
    int pxInd; int pyInd; int pzInd;

    for(pxInd = 0; pxInd <= legendreInd_X; pxInd++)
    {
    for(pyInd = 0; pyInd <= legendreInd_Y; pyInd++)
    {
    for(pzInd = 0; pzInd <= legendreInd_Z; pzInd++)
    {
    	evaluationVector[prodIndex] = productLegendre.evaluate(xEvalPos,pxInd,yEvalPos,pyInd,zEvalPos,pzInd);
    	prodIndex++;
    }}}

    // Create weight std::vector

    for(long i = 0; i < interpSystemSize; i++)
    {
    weightVector[i] = 0.0;
    for(long j = 0; j < interpSystemSize; j++)
    {
    weightVector[i] += interpMatrixInvTranspose(i,j)*evaluationVector[j];
    }}
}



void getIndexBounds(long interpIndex, int legendreInd, double h, long panel,
double pMin, double pPos, long& minIndex, long& maxIndex)
{
	if(legendreInd  == 0)
	{
    	minIndex = interpIndex; maxIndex = interpIndex;
	}
	else
	{
		if(interpIndex*h + pMin <= pPos)
		{
			maxIndex = interpIndex + (legendreInd + 1)/2;
			minIndex = maxIndex    -  legendreInd;
		}
		else
		{
			maxIndex = interpIndex +  legendreInd/2;
			minIndex = maxIndex    -  legendreInd;
		}
		if(maxIndex > panel){minIndex -= (maxIndex-panel); maxIndex = panel;}
		if(minIndex < 0    ){maxIndex += (-minIndex);          minIndex = 0;}

	}
}

void captureNodesAndWeights (double xPos, double yPos, double zPos,
      array<long,3>& coordPanels, array<double,3>& coordMin, array<double, 3>& coordMax,
      array<long,3>& minIndex, array<long, 3>& maxIndex,  std::vector<double>& weightVector)
{
	  createNodesAndWeightsData(xPos, yPos, zPos, coordPanels,coordMin,coordMax);

	  minIndex[0]  = this->xMinIndex;
	  maxIndex[0]  = this->xMaxIndex;
	  minIndex[1]  = this->yMinIndex;
	  maxIndex[1]  = this->yMaxIndex;
	  minIndex[2]  = this->zMinIndex;
	  maxIndex[2]  = this->zMaxIndex;

	  weightVector = this->weightVector;

	  /*
	  // Pack the weights into a 3d array
	  interpWeights.initialize(legendreInd_X+1,legendreInd_Y+1,legendreInd_Z+1);

	  long weightIndex = 0;
      for(long i = 0; i <= legendreInd_X; i++)
      {
      for(long j = 0; j <= legendreInd_Y; j++)
	  {
	  for(long k = 0; k <= legendreInd_Z; k++)
	  {
	    interpWeights(i,j,k) = weightVector[weightIndex];
	    weightIndex++;
	  }}}
	  */
}


void captureNodesAndWeights
      (double xPos, double yPos, double zPos, const SCC::GridFunction3d& F,
      long&  xMinIndex, long& xMaxIndex,
	  long&  yMinIndex, long& yMaxIndex,
      long&  zMinIndex, long& zMaxIndex,
      SCC::DoubleVector3d& interpWeights)
{
	  createNodesAndWeightsData(xPos, yPos, zPos, F);

	  xMinIndex = this->xMinIndex;
	  xMaxIndex = this->xMaxIndex;
	  yMinIndex = this->yMinIndex;
	  yMaxIndex = this->yMaxIndex;
	  zMinIndex = this->zMinIndex;
	  zMaxIndex = this->zMaxIndex;

	  //
	  // Pack the weights into a 3d array
	  interpWeights.initialize(legendreInd_X+1,legendreInd_Y+1,legendreInd_Z+1);

	  long weightIndex = 0;
      for(long i = 0; i <= legendreInd_X; i++)
      {
      for(long j = 0; j <= legendreInd_Y; j++)
	  {
	  for(long k = 0; k <= legendreInd_Z; k++)
	  {
	    interpWeights(i,j,k) = weightVector[weightIndex];
	    weightIndex++;
	  }}}
}
	// This routine creates the transpose of the
	// linear mapping from function values defined on a
	// grid to coefficients of an  expansion in products
	// of Legendre polynomials.

void createInterpolationMatrixInvTranspose()
{
    // This inverse transpose uses  product Legendre functions
	// defined in the coordinate system
    //
    // [-1.0,1.0]X[-1.0,1.0]X[-1.0,1.0]
    //
    // e.g. the standard coordinate system for the orthonormal
	// set of product Legendre functions.
    //
    ProductLegendrePoly3d productLegendre(-1.0,1.0,-1.0,1.0,-1.0,1.0);

    double hx; double hy; double hz;

    if(legendreInd_X == 0) {hx = 2.0;}
    else {hx = 2.0/(double)legendreInd_X;}

    if(legendreInd_Y == 0) {hy = 2.0;}
    else {hy = 2.0/(double)legendreInd_Y;}

    if(legendreInd_Z == 0) {hz = 2.0;}
    else {hz = 2.0/(double)legendreInd_Z;}

    long interpSystemSize = (legendreInd_X + 1)*(legendreInd_Y + 1)*(legendreInd_Z + 1);

    SCC::DoubleVector2d interpMatrix(interpSystemSize,interpSystemSize);

    long interpMatrixIndex;
    long prodIndex;

    int   pxInd;  int  pyInd;  int  pzInd;
    double xPos; double yPos; double zPos;

    long fxIndex; long fyIndex; long fzIndex;

    interpMatrixIndex = 0;
    for(fxIndex=0; fxIndex <= legendreInd_X; fxIndex++)
    {
    xPos = -1.0 + fxIndex*hx;
    for(fyIndex=0; fyIndex <= legendreInd_Y; fyIndex++)
    {
    yPos = -1.0 + fyIndex*hy;
    for(fzIndex=0; fzIndex <= legendreInd_Z; fzIndex++)
    {
    zPos = -1.0 + fzIndex*hz;

    prodIndex = 0;
    for(pxInd = 0; pxInd <= legendreInd_X; pxInd++)
    {
    for(pyInd = 0; pyInd <= legendreInd_Y; pyInd++)
    {
    for(pzInd = 0; pzInd <= legendreInd_Z; pzInd++)
    {
    interpMatrix(interpMatrixIndex,prodIndex) = productLegendre.evaluate(xPos,pxInd,yPos,pyInd,zPos,pzInd);
    prodIndex++;
    }}}

    interpMatrixIndex++;
    }}}

    // Create inverse transpose of interpolation matrix

    interpMatrixInvTranspose.initialize(interpSystemSize,interpSystemSize);

    interpMatrixInvTranspose.setToValue(0.0);

    for(long i = 0; i < interpSystemSize; i++)
    {
    	interpMatrixInvTranspose(i,i) = 1.0;
    }

    double* MdataPtr    = interpMatrix.getDataPointer();
    double* MinvDataPtr = interpMatrixInvTranspose.getDataPointer();
    SolveLinearSystem(MdataPtr, interpSystemSize, interpSystemSize, MinvDataPtr, interpSystemSize);
	}

	// Data required for constructing interpolation coefficients

    int   legendreInd_X;
    int   legendreInd_Y;
    int   legendreInd_Z;

	SCC::DoubleVector2d interpMatrixInvTranspose;
    std::vector<double>        evaluationVector;
    std::vector<double>        weightVector;

    long xMinIndex; long xMaxIndex;
    long yMinIndex; long yMaxIndex;
    long zMinIndex; long zMaxIndex;

//
//  Utility routine for computing inverse of interpolation matrix
//
//  (*) This routine computes the transpose of the inverse when B
//      is input with the identity matrix since it was created for
//      data stored by columns. For performance improvements this
//      routine should be replaced by a call to a Lapack routine.
//
private:

void SolveLinearSystem(double* MdataPtr, long M, long N,
double* BdataPtr, long P)
{
    long i,j,k;

    double machine_eps = ::pow(10.0,-12.0);

    long Mindex1size = M;
    long Mindex2size = N;
    long Bindex1size = M;
    long Bindex2size = P;
//
// create double temporaries
//
    double* RdataPtr = new double[Mindex1size*Mindex2size];
    double* XdataPtr = new double[Bindex1size*Bindex2size];
//
//  Put input data into temporaries :
//
//  Transpose copy of the matrix, since the original routine was written
//  assuming array data is stored by columns.
//
    for(i = 0; i < Mindex1size; i++)
    {
    for(j = 0; j < Mindex2size; j++)
    {
    *(RdataPtr + i + j*M) = *(MdataPtr + j + i*N);
    }}

//  Standard copy
//
    double* Bdptr;
    double* Xdptr;

    for(Xdptr = XdataPtr, Bdptr = BdataPtr;
    Xdptr < XdataPtr + (M*P); Xdptr++, Bdptr++)
    {*(Xdptr) = *(Bdptr);}

    long m     =   Mindex1size;
    long n     =   Mindex2size;
    long p     =   Bindex2size;

    double tau;

    double* Rptr; double* Sptr;
    double* Cptr; double* Tptr;
    double* Xptr;

    register double* Top;
    register double* piv_elem;

    double *S_data = new double[m];
    double *C_data = new double[m];
    double *T_data = new double[m];

    for( k = 1;  k <= n; k++)
    {
      piv_elem = RdataPtr + (k-1)*m + (k-1);
      for( Rptr  = piv_elem + 1, Sptr = S_data, Cptr = C_data;
           Rptr <= piv_elem + (m-k); Rptr++, Sptr++, Cptr++)
      {
        if( *Rptr == 0.0 ){ *Cptr = 1.0; *Sptr =0.0;}
        else
        { if( std::abs(*Rptr) > std::abs(*piv_elem) )
          { tau   = -(*piv_elem)/(*Rptr);
            *Sptr = 1.0/sqrt(1.0 + tau*tau);
            *Cptr = (*Sptr)*tau;
          }
          else
          { tau  = -(*Rptr)/(*piv_elem);
            *Cptr = 1.0/sqrt(1.0 + tau*tau);
            *Sptr = (*Cptr)*tau;
          }
        }
       *piv_elem = ((*Cptr) * (*piv_elem)) - ((*Sptr) * (*Rptr));
      }

      for( j = k+1; j <= n; j++)
      {
       Top = RdataPtr + (j-1)*m + (k-1);
       for(Rptr = Top + 1, Sptr = S_data, Cptr = C_data, Tptr = T_data;
           Rptr <= Top + (m - k); Rptr++, Sptr++, Cptr++, Tptr++)
            {  *Tptr = (*Sptr) * (*Top);
               *Top  = ((*Cptr) * (*Top)) - ((*Sptr) * (*Rptr));
               *Rptr *= (*Cptr);
            }

       for(Rptr = Top + 1, Tptr = T_data; Rptr <= Top + (m - k);  Rptr++, Tptr++)
           { *Rptr += *Tptr;}
      }
//
//    Transform the right hand side
//
      for( j = 1; j <= p ; j++)
      {
         Top = XdataPtr + (j-1)*m + (k-1);
         for( Xptr = Top + 1, Sptr = S_data, Cptr = C_data, Tptr = T_data;
              Xptr <= Top + (m - k); Xptr++, Sptr++, Cptr++, Tptr++)
              { *Tptr = (*Sptr) * (*Top);
                *Top  = ((*Cptr) * (*Top)) - ((*Sptr) * (*Xptr));
                *Xptr *= (*Cptr);
              }

         for( Xptr  = Top + 1, Tptr = T_data; Xptr <= Top + (m - k);
              Xptr++, Tptr++)
              { *Xptr += *Tptr; }
      }

}
//
//  Estimate the condition number
//
    double R_norm = 0.0;
    double R_col_norm;
    for( k = 1; k <= n; k++)
    {  R_col_norm = 0.0;
       for( Rptr = RdataPtr + (k-1)*m; Rptr < RdataPtr + (k-1)*m + k; Rptr++ )
       {R_col_norm += std::abs(*Rptr);}
       if(R_norm < R_col_norm ) R_norm = R_col_norm;
     }

    long singular_flag = 0;
    for( j=1; j <= n ; j++)
    {if(std::abs(*(RdataPtr + (j-1) + (j-1)*m)) <= machine_eps*R_norm ) singular_flag = 1;}

    if( singular_flag == 1)
    {
      printf(" Matrix system singular or close to singular \n");
      printf("     Computed results may be inaccurate      \n");
     }
//
//  Back Substitute
//
    register double XJ;
    for( k = 1; k <= p; k++)
    {
      for( j=  n; j >= 2; j--)
      {
       XJ = (*(XdataPtr +(k-1)*m + (j-1))) / (*(RdataPtr + (j-1)*m + (j-1)));
       *(XdataPtr +(k-1)*m + (j-1)) = XJ;
        for( Xptr = XdataPtr + (k-1)*m, Rptr = RdataPtr + (j-1)*m;
             Xptr < XdataPtr + (k-1)*m + (j-1); Xptr++, Rptr++)
            *Xptr -= XJ*(*Rptr);
      }
    *(XdataPtr + (k-1)*m) = (*(XdataPtr + (k-1)*m))/(*(RdataPtr));
    }
 //
 // Copy the solution to B
 //

    for(Xdptr = XdataPtr, Bdptr = BdataPtr;
    Xdptr < XdataPtr + (M*P); Xdptr++, Bdptr++)
    {*(Bdptr) = *(Xdptr);}


    delete [] XdataPtr;
    delete [] RdataPtr;
    delete [] S_data;
    delete [] C_data;
    delete [] T_data;
}

};

} // namespace SCC

#endif
