/*
 * LegendreGridFunApprox3d.h
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

#include <iostream>
#include <vector>
#include <functional>
#include <cmath>


#include "SCC_LegendreApprox3d.h"

#ifndef LEGENDRE_GRID_FUN_APPROX_3D_
#define LEGENDRE_GRID_FUN_APPROX_3D_

//
//
// This is a class whose purpose is to evaluate a local Legendre interpolant
// (and possibly the interpolant's derivatives) of discrete equispaced data.
//
// An instance of LegendreApprox3d is used to evaluate the Legendre
// approximation based upon local data values that is extracted by this class.
//
// The equispaced data upon which the approximation is based is access through a pointer,
// e.g. the data values are not captured.
//
//
// The requested interpolation degree must be > 0.
//
// If the number of grid panels is less than the specified interpolation degree
// in any coordinate direction then the interpolation degree in that coordinate
// direction is reduced.
//
namespace SCC
{

class LegendreGridFunApprox3d
{
public :

LegendreGridFunApprox3d()
{
	initialize();
}

LegendreGridFunApprox3d(long degreeX, long xDataPanels, double xDataMin, double xDataMax,
						long degreeY, long yDataPanels, double yDataMin, double yDataMax,
						long degreeZ, long zDataPanels, double zDataMin, double zDataMax,
                        const double* FDataPtr)
{
	initialize(degreeX,xDataPanels, xDataMin, xDataMax,
	           degreeY,yDataPanels, yDataMin, yDataMax,
	           degreeZ,zDataPanels, zDataMin, zDataMax, FDataPtr);
}

LegendreGridFunApprox3d(const LegendreGridFunApprox3d& Lapprox)
{
	initialize(Lapprox);
}


void setFortranFlag(bool fortranFlag = true)
{
	this->fortranFlag = fortranFlag;
}

void clearFortranFlag()
{
	fortranFlag = false;
}

void setDataPtr(const double* FDataPtr = nullptr)
{
	this->FDataPtr = FDataPtr;
}

void clearDataPtr()
{
	this->FDataPtr = nullptr;
}


void initialize()
{
	this->degreeX     = -1;
	this->xDataPanels = 0;
	this->xDataMin    = 0.0;
	this->xDataMax    = 0.0;
	this->hxData      = 0.0;

	this->degreeY     = -1;
	this->yDataPanels = 0;
	this->yDataMin    = 0.0;
	this->yDataMax    = 0.0;
	this->hyData      = 0.0;

	this->degreeZ     = -1;
	this->zDataPanels = 0;
	this->zDataMin    = 0.0;
	this->zDataMax    = 0.0;
	this->hzData      = 0.0;

	this->FDataPtr    = nullptr;

	Fvalues.clear();
	fortranFlag = false;
}

void initialize(long degreeX, long xDataPanels, double xDataMin, double xDataMax,
			    long degreeY, long yDataPanels, double yDataMin, double yDataMax,
			    long degreeZ, long zDataPanels, double zDataMin, double zDataMax,
                const double* FDataPtr)
{
    if(xDataPanels < degreeX){degreeX = xDataPanels;}
	if(yDataPanels < degreeY){degreeY = yDataPanels;}
	if(zDataPanels < degreeZ){degreeZ = zDataPanels;}

    long NX = (degreeX+1);
    long NY = (degreeY+1);
    long NZ = (degreeZ+1);

    if((this->degreeX != degreeX)||
       (this->degreeY != degreeY)||
       (this->degreeZ != degreeZ)
       )
    {
	this->Fvalues.resize(NX*NY*NZ);
	legendreApprox3d.initialize(degreeX,degreeY,degreeZ);
    }



	this->degreeX     = degreeX;
	this->xDataPanels = xDataPanels;
	this->xDataMin    = xDataMin;
	this->xDataMax    = xDataMax;
	this->hxData      = (xDataMax-xDataMin)/(double)xDataPanels;

	this->degreeY     = degreeY;
	this->yDataPanels = yDataPanels;
	this->yDataMin    = yDataMin;
	this->yDataMax    = yDataMax;
	this->hyData      = (yDataMax-yDataMin)/(double)yDataPanels;

	this->degreeZ     = degreeZ;
	this->zDataPanels = zDataPanels;
	this->zDataMin    = zDataMin;
	this->zDataMax    = zDataMax;
	this->hzData      = (zDataMax-zDataMin)/(double)zDataPanels;

    this->FDataPtr    = FDataPtr;      // Shallow capture


	fortranFlag = false;
}

void initialize(const LegendreGridFunApprox3d& Lapprox)
{
	this->degreeX     = Lapprox.degreeX;
	this->xDataPanels = Lapprox.xDataPanels;
	this->xDataMin    = Lapprox.xDataMin;
	this->xDataMax    = Lapprox.xDataMax;
	this->hxData      = (xDataMax-xDataMin)/(double)xDataPanels;

	this->degreeY     = Lapprox.degreeY;
	this->yDataPanels = Lapprox.yDataPanels;
	this->yDataMin    = Lapprox.yDataMin;
	this->yDataMax    = Lapprox.yDataMax;
	this->hyData      = (yDataMax-yDataMin)/(double)yDataPanels;

	this->degreeZ     = Lapprox.degreeZ;
	this->zDataPanels = Lapprox.zDataPanels;
	this->zDataMin    = Lapprox.zDataMin;
	this->zDataMax    = Lapprox.zDataMax;
	this->hzData      = (zDataMax-zDataMin)/(double)zDataPanels;

    this->FDataPtr    = Lapprox.FDataPtr;  // Shallow copy

    long NX = (degreeX+1);
    long NY = (degreeY+1);
    long NZ = (degreeZ+1);

	this->Fvalues.resize(NX*NY*NZ);
	legendreApprox3d  = Lapprox.legendreApprox3d;

	fortranFlag = false;
}

double evaluate(double x,double y,double z)
{
//
// If evaluation point is coincident with a grid point of data being
// interpolated then skip interpolation and return data value directly
//
   long interpIndexX = (long)std::round((x-xDataMin)/hxData);
   if(interpIndexX < 0)            interpIndexX = 0;
   if(interpIndexX >= xDataPanels) interpIndexX = xDataPanels-1;

   long interpIndexY = (long)std::round((y-yDataMin)/hyData);
   if(interpIndexY < 0)            interpIndexY = 0;
   if(interpIndexY >= yDataPanels) interpIndexY = yDataPanels-1;

   long interpIndexZ = (long)std::round((z-zDataMin)/hzData);
   if(interpIndexZ < 0)            interpIndexZ = 0;
   if(interpIndexZ >= zDataPanels) interpIndexZ = zDataPanels-1;

   if( ( std::abs(x - (interpIndexX*hxData + xDataMin)) < (std::abs(x) + 1.0e-06)*1.0e-12 )
       &&
       ( std::abs(y - (interpIndexY*hyData + yDataMin)) < (std::abs(y) + 1.0e-06)*1.0e-12 )
       &&
       ( std::abs(x - (interpIndexZ*hzData + zDataMin)) < (std::abs(z) + 1.0e-06)*1.0e-12 )
       )
   {
	   if(not fortranFlag)
	   {
	       return
	       FDataPtr[interpIndexZ + (interpIndexY*(zDataPanels+1))
	                             + (interpIndexX*(zDataPanels+1)*(yDataPanels+1))];
	   }
	   else
	   {
		   return
		   FDataPtr[interpIndexX + (interpIndexY*(xDataPanels+1))
		                         + (interpIndexZ*(xDataPanels+1)*(yDataPanels+1))];
	   }
   }

   setLocalData(x,y,z);
   return legendreApprox3d.evaluate(x,interpXmin,interpXmax,
                                    y,interpYmin,interpYmax,
		   	   	   	   	   	   	    z,interpZmin,interpZmax,Fvalues);
}

void evaluateDerivative(double x,double y, double z, std::vector<double>& dFvalues)
{
   setLocalData(x,y,z);
   return legendreApprox3d.evaluateDerivative(x,interpXmin,interpXmax,
                                              y,interpYmin,interpYmax,
		   	   	   	   	   	   	   	   	   	  z,interpZmin,interpZmax,
                                                Fvalues,dFvalues);
}

void evaluate(double x, double y, double z, double& Fval, std::vector<double>& dFvalues)
{
	setLocalData(x,y,z);
    legendreApprox3d.evaluate(x,interpXmin,interpXmax,
                              y,interpYmin,interpYmax,
		   	   	   	   	   	  z,interpZmin,interpZmax,
                                Fvalues,Fval,dFvalues);
}


void setLocalData(double x,double y,double z)
{
   long    interpIndexBaseX;
   long    interpIndexStartX;
   long      interpIndexEndX;

   long    interpIndexBaseY;
   long    interpIndexStartY;
   long      interpIndexEndY;

   long    interpIndexBaseZ;
   long    interpIndexStartZ;
   long      interpIndexEndZ;

   // X-coordinate

   interpIndexBaseX = floor((x-xDataMin)/hxData);
   if(interpIndexBaseX < 0)            interpIndexBaseX = 0;
   if(interpIndexBaseX >= xDataPanels) interpIndexBaseX = xDataPanels-1;

   interpIndexStartX = interpIndexBaseX;
   interpIndexEndX   = interpIndexBaseX;


   interpIndexStartX = interpIndexBaseX;
   interpIndexEndX   = interpIndexBaseX + 1;

   while(interpIndexEndX - interpIndexStartX < degreeX)
   {
		   if(interpIndexStartX > 0          ) interpIndexStartX -= 1;
		   if(interpIndexEndX   < xDataPanels) interpIndexEndX   += 1;
   }

   if( (interpIndexEndX - interpIndexStartX) ==  degreeX+1 )
   {
   if( (x - (xDataMin + interpIndexBaseX*hxData))/hxData > 0.5 ){interpIndexStartX +=1;}
   else                                                         {interpIndexEndX   -=1;}
   }

   interpXmin = xDataMin + interpIndexStartX*hxData;
   interpXmax = xDataMin + interpIndexEndX*hxData;

   // Y-coordinate

   interpIndexBaseY = floor((y-yDataMin)/hyData);
   if(interpIndexBaseY < 0)            interpIndexBaseY = 0;
   if(interpIndexBaseY >= yDataPanels) interpIndexBaseY = yDataPanels-1;

   interpIndexStartY = interpIndexBaseY;
   interpIndexEndY   = interpIndexBaseY + 1;

   while(interpIndexEndY - interpIndexStartY < degreeY)
   {
		   if(interpIndexStartY > 0          ) interpIndexStartY -= 1;
		   if(interpIndexEndY   < yDataPanels) interpIndexEndY   += 1;
   }

   if( (interpIndexEndY- interpIndexStartY) ==  degreeY+1 )
   {
   if( (y - (yDataMin + interpIndexBaseY*hyData))/hyData > 0.5 ){interpIndexStartY +=1;}
   else                                                         {interpIndexEndY   -=1;}
   }

   interpYmin = yDataMin + interpIndexStartY*hyData;
   interpYmax = yDataMin + interpIndexEndY*hyData;

   // Z-coordinate

   interpIndexBaseZ = floor((z-zDataMin)/hzData);
   if(interpIndexBaseZ < 0)            interpIndexBaseZ = 0;
   if(interpIndexBaseZ >= zDataPanels) interpIndexBaseZ = zDataPanels-1;

   interpIndexStartZ = interpIndexBaseZ;
   interpIndexEndZ   = interpIndexBaseZ + 1;

   while(interpIndexEndZ - interpIndexStartZ < degreeZ)
   {
		   if(interpIndexStartZ > 0          ) interpIndexStartZ -= 1;
		   if(interpIndexEndZ   < zDataPanels) interpIndexEndZ   += 1;
   }

   if( (interpIndexEndZ- interpIndexStartZ) ==  degreeZ+1 )
   {
   if( (z - (zDataMin + interpIndexBaseZ*hzData))/hzData > 0.5 ){interpIndexStartZ +=1;}
   else                                                         {interpIndexEndZ   -=1;}
   }

   interpZmin = zDataMin + interpIndexStartZ*hzData;
   interpZmax = zDataMin + interpIndexEndZ*hzData;


   long NX = (degreeX+1);
   long NY = (degreeY+1);
   long NZ = (degreeZ+1);

   long NXdata = xDataPanels+1;
   long NYdata = yDataPanels+1;
   long NZdata = zDataPanels+1;

   if(not fortranFlag)
   {

   for(long i = 0; i < NX; i++)
   {
   for(long j = 0; j < NY; j++)
   {
   for(long k = 0; k < NZ; k++)
   {
   Fvalues[k + j*NZ + i*NZ*NY]
   = FDataPtr[(interpIndexStartZ + k) + (interpIndexStartY + j)*NZdata + (interpIndexStartX + i)*NZdata*NYdata];
   }}}

   }

   if( fortranFlag )
   {
   for(long i = 0; i < NX; i++)
   {
   for(long j = 0; j < NY; j++)
   {
   for(long k = 0; k < NZ; k++)
   {
   Fvalues[k + j*NZ + i*NZ*NY]
   = FDataPtr[ (interpIndexStartX + i) + (interpIndexStartY + j)*NXdata + (interpIndexStartZ + k)*NXdata*NYdata];
   }}}

   }
}

    long    degreeX;
    long    xDataPanels;
	double  xDataMin;
	double  xDataMax;
	double  hxData;

    long    degreeY;
    long    yDataPanels;
	double  yDataMin;
	double  yDataMax;
	double  hyData;

    long    degreeZ;
    long    zDataPanels;
	double  zDataMin;
	double  zDataMax;
	double  hzData;


	const double* FDataPtr;

	LegendreApprox3d legendreApprox3d;

	std::vector<double>   Fvalues;

	double           interpXmin;
	double           interpXmax;

	double           interpYmin;
	double           interpYmax;

	double           interpZmin;
	double           interpZmax;

	bool             fortranFlag;
};



} // namespace SCC

#endif /* _LegendreGridFunApprox3d_ */
