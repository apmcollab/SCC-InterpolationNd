/*
 * LegendreGridFunApprox2d.h
 *
 *  Created on: Jul 22, 2018
 *      Author: anderson
 */

/*
#############################################################################
#
# Copyright  2018- Chris Anderson
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



#include "SCC_LegendreApprox2d.h"

#ifndef LEGENDRE_GRID_FUN_APPROX_2D_
#define LEGENDRE_GRID_FUN_APPROX_2D_

//
//
// This is a class whose purpose is to evaluate a local Legendre interpolant
// (and possibly the interpolant's derivatives) of discrete equispaced data.
//
// An instance of LegendreApprox2d is used to evaluate the Legendre
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

class LegendreGridFunApprox2d
{
public :

LegendreGridFunApprox2d()
{
	initialize();
}

LegendreGridFunApprox2d(long degreeX, long xDataPanels, double xDataMin, double xDataMax,
						long degreeY, long yDataPanels, double yDataMin, double yDataMax,
                        const double* FDataPtr)
{
	initialize(degreeX,xDataPanels, xDataMin, xDataMax,
	           degreeY,yDataPanels, yDataMin, yDataMax, FDataPtr);
}

LegendreGridFunApprox2d(const LegendreGridFunApprox2d& Lapprox)
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

void setDataPtr(double* FDataPtr = nullptr)
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

	this->FDataPtr    = nullptr;

	Fvalues.clear();

	fortranFlag = false;
}

void initialize(long degreeX, long xDataPanels, double xDataMin, double xDataMax,
			    long degreeY, long yDataPanels, double yDataMin, double yDataMax,
                const double* FDataPtr)
{
	if(xDataPanels < degreeX){degreeX = xDataPanels;}
	if(yDataPanels < degreeY){degreeY = yDataPanels;}

    long NX = (degreeX+1);
    long NY = (degreeY+1);

    if((this->degreeX != degreeX)||(this->degreeY != degreeY))
    {
	this->Fvalues.resize(NX*NY);
	legendreApprox2d.initialize(degreeX,degreeY);
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

    this->FDataPtr    = FDataPtr;      // Shallow capture

    fortranFlag = false;
}

void initialize(const LegendreGridFunApprox2d& Lapprox)
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

    this->FDataPtr    = Lapprox.FDataPtr;  // Shallow copy

    long NX = (degreeX+1);
    long NY = (degreeY+1);

	this->Fvalues.resize(NX*NY);
	legendreApprox2d  = Lapprox.legendreApprox2d;

    fortranFlag = false;
}

double evaluate(double x,double y)
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


   if( ( std::abs(x - (interpIndexX*hxData + xDataMin)) < (std::abs(x) + 1.0e-06)*1.0e-12 )
       &&
       ( std::abs(y - (interpIndexY*hyData + yDataMin)) < (std::abs(y) + 1.0e-06)*1.0e-12 ) )
   {
	   if(not fortranFlag)
	   {
		   return FDataPtr[interpIndexY  + (interpIndexX)*(yDataPanels+1)];
	   }
	   else
	   {
		   return FDataPtr[interpIndexX  + (interpIndexY)*(xDataPanels+1)];
	   }
   }

   setLocalData(x,y);
   return legendreApprox2d.evaluate(x,interpXmin,interpXmax,y,interpYmin,interpYmax,Fvalues);
}

void evaluateDerivative(double x,double y, std::vector<double>& dFvalues)
{
   setLocalData(x,y);
   return legendreApprox2d.evaluateDerivative(x,interpXmin,interpXmax,y,interpYmin,interpYmax,
                                             Fvalues,dFvalues);
}

void evaluate(double x, double y, double& Fval, std::vector<double>& dFvalues)
{
	setLocalData(x,y);
    legendreApprox2d.evaluate(x,interpXmin,interpXmax,y,interpYmin,interpYmax,
                              Fvalues,Fval,dFvalues);
}


void setLocalData(double x,double y)
{
   long    interpIndexBaseX;
   long    interpIndexStartX;
   long      interpIndexEndX;

   long    interpIndexBaseY;
   long    interpIndexStartY;
   long      interpIndexEndY;

   interpIndexBaseX = (long)floor((x-xDataMin)/hxData);
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

   interpIndexBaseY = (long)floor((y-yDataMin)/hyData);
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

   long NX = (degreeX+1);
   long NY = (degreeY+1);

   long NXdata = xDataPanels+1;
   long NYdata = yDataPanels+1;

   // Capture local data assuming C/C++ data storage convention

   if(not fortranFlag)
   {
	   for(long i = 0; i < NX; i++)
	   {
		   for(long j = 0; j < NY; j++)
		   {
			   Fvalues[j + i*NY] = FDataPtr[(interpIndexStartY + j) + (interpIndexStartX + i)*NYdata];
	   }}
   }

   if(fortranFlag)
   {
   	   for(long i = 0; i < NX; i++)
	   {
		   for(long j = 0; j < NY; j++)
		   {
			   Fvalues[j + i*NY] = FDataPtr[(interpIndexStartX + i) + (interpIndexStartY + j)*NXdata];
	   }}
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

	const double* FDataPtr;

	LegendreApprox2d legendreApprox2d;

	std::vector<double>   Fvalues;

	double           interpXmin;
	double           interpXmax;

	double           interpYmin;
	double           interpYmax;

	bool             fortranFlag;
};


} // namespace SCC


#endif /* _LegendreGridFunApprox2d_ */
