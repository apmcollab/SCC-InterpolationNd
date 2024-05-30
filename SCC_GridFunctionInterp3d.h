/*
 * GridFunctionInterp3d.h
 *
 *  Created on: Jul 25, 2022
 *      Author: anderson
 *
 *
 * ====> alpha version <===
 *
 * This class provides member functions that interpolate values from a source
 * SCC::GridFunction3d instance to a SCC::GridFunction3d target instance.
 *
 * The method of interpolation is a tensor product of Legendre functions that
 * interpolate the source data values. The data values used to form the interpolant
 * are chosen so that the evaluation is as centered as much as possible with respect
 * to the source value locations.
 *
 * The XYZ interpolation is a full three dimensional interpolation. If the grid
 * points of the target GridFunction3d lie outside the domain of the source domain
 * then the approximation is an extrapolation. This is not result in an error condition
 * being thrown, but generally results in a very inaccurate approximation.
 *
 * The XY interpolation routine is a full three dimensional interpolation but is formed
 * using only interpolation in the x-y direction for each z-coordinate. This method assumes
 * that the source function and target function possess the same grid structure in the
 * vertical direction. If this assumption isn't satisfied then an error is thrown.
 *
 * The XY interpolation is most useful for functions that are continuous but only piecewise
 * smooth in the vertical direction. In such instances, functions that are continuous but
 * has discontinuous derivatives along specific z-coordinates will will be still be interpolated
 * with high accuracy in the transverse directions.
 *
 */

/*
#############################################################################
#
# Copyright  2022- Chris Anderson
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
#include "GridFunctionNd/SCC_GridFunction2d.h"
#include "GridFunctionNd/SCC_GridFunction3d.h"

#include "InterpolationNd/SCC_LegendreGridFunApprox2d.h"
#include "InterpolationNd/SCC_LegendreGridFunApprox3d.h"

#ifndef GRID_FUNCTION_INTERP_3D_
#define GRID_FUNCTION_INTERP_3D_

namespace SCC
{

class GridFunctionInterp3d
{
public:

	GridFunctionInterp3d()
	{
		initialize();
	}

	GridFunctionInterp3d(GridFunctionInterp3d& V)
	{
		initialize(V);
	}

	GridFunctionInterp3d(int degreeX,int degreeY, int degreeZ, SCC::GridFunction3d& dataFun3d)
	{
		initialize(degreeX,degreeY, degreeZ,dataFun3d);
	}

	void initialize()
	{

    dataFun3d.initialize();
    dataXY.initialize();

	legendreGridFunApprox2d.initialize();
	legendreGridFunApprox3d.initialize();
	}

	void initialize(GridFunctionInterp3d& T)
	{
		dataFun3d.initialize(T.dataFun3d);

		dataXY.initialize(T.dataXY);

		legendreGridFunApprox2d.initialize(T.legendreGridFunApprox2d);
		legendreGridFunApprox3d.initialize(T.legendreGridFunApprox3d);

		// Set pointer to point to *this instance of the data

		legendreGridFunApprox3d.FDataPtr = dataFun3d.getDataPointer();
	}

	void initialize(int degreeX,int degreeY, int degreeZ, SCC::GridFunction3d& dataFun3d)
	{
		 if(dataFun3d.isNull())
		 {
			throw std::runtime_error("\nGridFunctionInterp3d: initialize(...)SCC::GridFunction3d argument not initialized.\n");
		 }

		 this->dataFun3d.initialize(dataFun3d);

		 long xPanels; double xMin; double xMax;
		 long yPanels; double yMin; double yMax;
		 long zPanels; double zMin; double zMax;

         xPanels = dataFun3d.getXpanelCount();
         xMin    = dataFun3d.getXmin();
         xMax    = dataFun3d.getXmax();

         yPanels = dataFun3d.getYpanelCount();
         yMin    = dataFun3d.getYmin();
         yMax    = dataFun3d.getYmax();

         zPanels = dataFun3d.getZpanelCount();
         zMin    = dataFun3d.getZmin();
         zMax    = dataFun3d.getZmax();

         dataXY.initialize(xPanels,xMin,xMax,yPanels,yMin,yMax);

         legendreGridFunApprox3d.initialize(degreeX, xPanels, xMin, xMax,
        		                                  degreeY, yPanels, yMin, yMax,
												  degreeZ, zPanels, zMin, zMax,
												  dataFun3d.getDataPointer());

         legendreGridFunApprox2d.initialize(degreeX, xPanels, xMin, xMax,
        		                               degreeY, yPanels, yMin, yMax,
											   nullptr);
	}

	//
	// Three dimensional interpolation.
	//
	// The domain size of the target function need not be coincident with the
	// source (data) function
	//
	void interpXYZ(SCC::GridFunction3d& outFun)
	{
		 if(outFun.isNull())
		 {
			throw std::runtime_error("\n GridFunctionInterp3d: interpXYZ(...) return argument not initialized.\n");
		 }

		 long xPanels; double xMin; double xMax;
		 long yPanels; double yMin; double yMax;
		 long zPanels; double zMin; double zMax;

		 xPanels = outFun.getXpanelCount();
         xMin    = outFun.getXmin();
         xMax    = outFun.getXmax();

         yPanels = outFun.getYpanelCount();
         yMin    = outFun.getYmin();
         yMax    = outFun.getYmax();

         zPanels = outFun.getZpanelCount();
         zMin    = outFun.getZmin();
         zMax    = outFun.getZmax();

		 double hx; double hy; double hz;

		 hx = (xMax-xMin)/(double)xPanels;
		 hy = (yMax-yMin)/(double)yPanels;
		 hz = (zMax-zMin)/(double)zPanels;

		 double xPos; double yPos; double zPos;

		 for(long p = 0; p <= xPanels; p++)
		 {
	     xPos = xMin + p*hx;

		 for(long q = 0; q <= yPanels;  q++)
		 {
		 yPos = yMin + q*hy;

		 for(long r = 0; r  <= zPanels; r++)
		 {
		 zPos = zMin + r*hz;
		 outFun(p,q,r) = legendreGridFunApprox3d.evaluate(xPos, yPos, zPos);
		 }}}
	}

	//
	// Three dimensional interpolation of data values to a target function that possesses the same vertical structure
	// as the source (data) function. This interpolation just interpolates the 2D slices for each
	// discrete z-coordinate value.
	//
	// This form of interpolation will accurately interpolate functions that are just continuous
	// and piecewise smooth in the vertical direction
	//
	// The transverse domain size of the target function need not be coincident with the
	// source (data) function
	//
	void interpXY(SCC::GridFunction3d& outFun)
	{
		 if(outFun.isNull())
		 {
			throw std::runtime_error("\n GridFunctionInterp3d: interpXYZ(...) return argument not initialized.\n");
		 }

	    // Check for consistency of vertical structure

        SCC::GridFunction1d dataFunZ = dataFun3d.getConstantXYslice(0, 0);
        SCC::GridFunction1d  outFunZ = outFun.getConstantXYslice(0, 0);

        if(not outFunZ.isCoincident(dataFunZ))
        {
        	throw std::runtime_error(std::string("\n GridFunctionInterp3d: interpXY(...) vertical structure \n")
        			                             + " of GridFunction3d argument not consistent with data\n");
        }

        SCC::GridFunction2d valuesXY;

        long xPanels; double xMin; double xMax;
        long yPanels; double yMin; double yMax;
        long zPanels;
		xPanels = outFun.getXpanelCount();
        xMin    = outFun.getXmin();
        xMax    = outFun.getXmax();

        yPanels = outFun.getYpanelCount();
        yMin    = outFun.getYmin();
        yMax    = outFun.getYmax();

        zPanels = outFun.getZpanelCount();

		 double hx; double hy;

		 hx = (xMax-xMin)/(double)xPanels;
		 hy = (yMax-yMin)/(double)yPanels;

		 double xPos; double yPos;

		 valuesXY.initialize(xPanels,xMin,xMax,yPanels,yMin,yMax);

		 for(long r = 0; r <= zPanels; r++)
		 {
             dataFun3d.getConstantZslice(r,dataXY);

		     legendreGridFunApprox2d.FDataPtr = dataXY.getDataPointer();

		     for(long p = 0; p <= xPanels; p++)
		     {
		     xPos = xMin + p*hx;
		     for(long q = 0; q <= yPanels;  q++)
		     {
		     yPos = yMin + q*hy;
		     valuesXY(p,q) = legendreGridFunApprox2d.evaluate(xPos, yPos);
		     }}

			 outFun.setConstantZslice(r,valuesXY);
			 }
	}

	SCC::GridFunction2d getConstantXslice(double xPos)
	{
	double yMin  = legendreGridFunApprox3d.yDataMin;
	double yMax  = legendreGridFunApprox3d.yDataMax;
	double hy    = legendreGridFunApprox3d.hyData;
	long yPanels = legendreGridFunApprox3d.yDataPanels;

    double zMin  = legendreGridFunApprox3d.zDataMin;
	double zMax  = legendreGridFunApprox3d.zDataMax;
	double hz    = legendreGridFunApprox3d.hzData;
	long zPanels = legendreGridFunApprox3d.zDataPanels;

    SCC::GridFunction2d R(yPanels,yMin,yMax,zPanels,zMin,zMax);

	double yPos; double zPos;

    for(long j = 0; j <= yPanels; j++)
    {
    yPos = yMin + j*hy;
    for(long k = 0; k <= zPanels; k++)
    {
    zPos = zMin + k*hz;
    R.Values(j,k) = legendreGridFunApprox3d.evaluate(xPos,yPos,zPos);
    }}
    return R;
	}

	SCC::GridFunction2d getConstantYslice(double yPos)
	{
	double xMin  = legendreGridFunApprox3d.xDataMin;
	double xMax  = legendreGridFunApprox3d.xDataMax;
	double hx    = legendreGridFunApprox3d.hxData;
	long xPanels = legendreGridFunApprox3d.xDataPanels;

    double zMin  = legendreGridFunApprox3d.zDataMin;
	double zMax  = legendreGridFunApprox3d.zDataMax;
	double hz    = legendreGridFunApprox3d.hzData;
	long zPanels = legendreGridFunApprox3d.zDataPanels;

    SCC::GridFunction2d R(xPanels,xMin,xMax,zPanels,zMin,zMax);

	double xPos; double zPos;

    for(long i = 0; i <= xPanels; i++)
    {
    xPos = xMin + i*hx;
    for(long k = 0; k <= zPanels; k++)
    {
    zPos = zMin + k*hz;
    R.Values(i,k) = legendreGridFunApprox3d.evaluate(xPos,yPos,zPos);
    }}
    return R;
	}


	SCC::GridFunction2d getConstantZslice(double zPos)
	{
	double xMin  = legendreGridFunApprox3d.xDataMin;
	double xMax  = legendreGridFunApprox3d.xDataMax;
	double hx    = legendreGridFunApprox3d.hxData;
	long xPanels = legendreGridFunApprox3d.xDataPanels;

    double yMin  = legendreGridFunApprox3d.yDataMin;
	double yMax  = legendreGridFunApprox3d.yDataMax;
	double hy    = legendreGridFunApprox3d.hyData;
	long yPanels = legendreGridFunApprox3d.yDataPanels;

    SCC::GridFunction2d R(xPanels,xMin,xMax,yPanels,yMin,yMax);

	double xPos; double yPos;

    for(long i = 0; i <= xPanels; i++)
    {
    xPos = xMin + i*hx;
    for(long j = 0; j <= yPanels; j++)
    {
    yPos = yMin + j*hy;
    R.Values(i,j) = legendreGridFunApprox3d.evaluate(xPos,yPos,zPos);
    }}
    return R;
	}


	SCC::GridFunction3d   dataFun3d;
    SCC::GridFunction2d   dataXY;

	LegendreGridFunApprox2d legendreGridFunApprox2d;
	LegendreGridFunApprox3d legendreGridFunApprox3d;
};


} // SCC namespace
#endif /* GRID_FUNCTION_INTERP_3D_ */
