
#include "LegendreInterp3d.h"
#include "DoubleVectorNd/UCLAQ_DoubleVector3d.h"
#include "GridFunctionNd/UCLAQ_GridFunction3d.h"

#ifndef _FixedLocationLegendreInterp3d_
#define _FixedLocationLegendreInterp3d_


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


class FixedLocationLegendreInterp3d
{
public:

    FixedLocationLegendreInterp3d():interpWeights()
    {
	xMinIndex=0; xMaxIndex=0;
	yMinIndex=0; yMaxIndex=0;
	zMinIndex=0; zMaxIndex=0;
    }


	FixedLocationLegendreInterp3d(int interpolationOrder,double xPos, double yPos, double zPos,
    const UCLAQ::GridFunction3d& F)
    {
	LegendreInterp3d interp3d(interpolationOrder);
    interp3d.captureNodesAndWeights(xPos,yPos,zPos,F,
    		xMinIndex,xMaxIndex,
			yMinIndex,yMaxIndex,
			zMinIndex,zMaxIndex,
			interpWeights);
    }

	FixedLocationLegendreInterp3d(const FixedLocationLegendreInterp3d& P)
    {
	xMinIndex= P.xMinIndex; xMaxIndex= P.xMaxIndex;
	yMinIndex= P.yMinIndex; yMaxIndex= P.yMaxIndex;
	zMinIndex= P.zMinIndex; zMaxIndex= P.zMaxIndex;
	interpWeights.initialize(P.interpWeights);
    }

	void initialize(int interpolationOrder,double xPos, double yPos, double zPos,
    const UCLAQ::GridFunction3d& F)
    {
	LegendreInterp3d interp3d(interpolationOrder);
    interp3d.captureNodesAndWeights(xPos,yPos,zPos,F,
    		xMinIndex,xMaxIndex,
			yMinIndex,yMaxIndex,
			zMinIndex,zMaxIndex,
			interpWeights);
    }

	double getInterpolatedValue(const UCLAQ::GridFunction3d& F)
	{
    long i; long j; long k;
    long fxIndex; long fyIndex; long fzIndex;

    double interpolationValue  = 0.0;
    for(i = xMinIndex,fxIndex=0; i <= xMaxIndex; i++,fxIndex++)
    {
    for(j = yMinIndex,fyIndex=0; j <= yMaxIndex; j++,fyIndex++)
    {
    for(k = zMinIndex,fzIndex=0; k <= zMaxIndex; k++,fzIndex++)
    {
    interpolationValue += F(i,j,k)*interpWeights(fxIndex,fyIndex,fzIndex);
    }}}

    return interpolationValue;
	}

	long xMinIndex; long xMaxIndex;
	long yMinIndex; long yMaxIndex;
	long zMinIndex; long zMaxIndex;

	UCLAQ::DoubleVector3d interpWeights;

};
#endif
