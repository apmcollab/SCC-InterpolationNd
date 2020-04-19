/*
 * LegendreGridFunApprox1d.h
 *
 *  Created on: Jul 22, 2018
 *      Author: anderson
 */

#include <iostream>
#include <vector>
#include <functional>
#include <cmath>



#include "LegendreApprox1d.h"

#ifndef _LegendreGridFunApprox1d_
#define _LegendreGridFunApprox1d_

//
//
// This is a class whose purpose is to evaluate a local Legendre interpolant
// (and possibly the interpolant's derivatives) of discrete equispaced data.
//
// An instance of LegendreApprox1d is used to evaluate the Legendre
// approximation based upon local data values that is extracted by this class.
//
// The equispaced data upon which the approximation is based is access through a pointer,
// e.g. the data values are not captured.
//
//
// The requested interpolation degree must be > 0.
//
//

class LegendreGridFunApprox1d
{
public :

LegendreGridFunApprox1d()
{
	initialize();
}

LegendreGridFunApprox1d(long degreeX, long xDataPanels, double xDataMin, double xDataMax,double* FDataPtr)
{
	initialize(degreeX,xDataPanels, xDataMin, xDataMax,FDataPtr);
}

LegendreGridFunApprox1d(const LegendreGridFunApprox1d& Lapprox)
{
	initialize(Lapprox);
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
	degreeX      = -1;
	xDataPanels = 0;
	xDataMin    = 0.0;
	xDataMax    = 0.0;
	hData       = 0.0;

    FDataPtr    = nullptr;
	Fvalues.clear();
}

void initialize(long degreeX, long xDataPanels, double xDataMin, double xDataMax,double* FDataPtr)
{
    if(this->degreeX != degreeX)
    {
    this->Fvalues.resize(degreeX+1);
	legendreApprox1d.initialize(degreeX);
    }

	this->degreeX      = degreeX;
	this->xDataPanels = xDataPanels;
	this->FDataPtr    = FDataPtr;
	this->xDataMin    = xDataMin;
	this->xDataMax    = xDataMax;
	this->hData       = (xDataMax-xDataMin)/(double)xDataPanels;
}

void initialize(const LegendreGridFunApprox1d& Lapprox)
{
	this->degreeX      = Lapprox.degreeX;
	this->xDataPanels  = Lapprox.xDataPanels;

	this->xDataMin    = Lapprox.xDataMin;
	this->xDataMax    = Lapprox.xDataMax;
	this->hData       = Lapprox.hData;

    this->FDataPtr    = Lapprox.FDataPtr;
	this->Fvalues.resize(degreeX+1);

	legendreApprox1d  = Lapprox.legendreApprox1d;
}

double evaluate(double x)
{
   setLocalData(x);
   return legendreApprox1d.evaluate(x,interpXmin,interpXmax,Fvalues);
}

double evaluateDerivative(double x)
{
   setLocalData(x);
   return legendreApprox1d.evaluateDerivative(x,interpXmin,interpXmax,Fvalues);
}

void evaluate(double x, double& fVal, double &dfVal)
{
	setLocalData(x);
    legendreApprox1d.evaluate(x,interpXmin,interpXmax,Fvalues,fVal,dfVal);

}


void setLocalData(double x)
{
   long    interpIndexBase;
   long    interpIndexStart;
   long      interpIndexEnd;

   interpIndexBase = floor((x-xDataMin)/hData);
   if(interpIndexBase < 0)            interpIndexBase = 0;
   if(interpIndexBase >= xDataPanels) interpIndexBase = xDataPanels-1;

   interpIndexStart = interpIndexBase;
   interpIndexEnd   = interpIndexBase;


   interpIndexStart = interpIndexBase;
   interpIndexEnd   = interpIndexBase + 1;

   while(interpIndexEnd - interpIndexStart < degreeX)
   {
		   if(interpIndexStart > 0          ) interpIndexStart -= 1;
		   if(interpIndexEnd   < xDataPanels) interpIndexEnd   += 1;
   }


   if( (interpIndexEnd - interpIndexStart) ==  degreeX+1 )
   {
   if( (x - (xDataMin + interpIndexBase*hData))/hData > 0.5 ){interpIndexStart +=1;}
   else                                                      {interpIndexEnd   -=1;}
   }

   interpXmin = xDataMin + interpIndexStart*hData;
   interpXmax = xDataMin + interpIndexEnd*hData;

   for(long i = 0; i <= degreeX; i++)
   {
   Fvalues[i] = FDataPtr[interpIndexStart + i];
   }
}

    long    degreeX;
    long    xDataPanels;
	double  xDataMin;
	double  xDataMax;
	double  hData;

    double* FDataPtr;

	LegendreApprox1d legendreApprox1d;
	std::vector<double>   Fvalues;
	double           interpXmin;
	double           interpXmax;
};




#endif /* _LegendreGridFunApprox1d_ */
