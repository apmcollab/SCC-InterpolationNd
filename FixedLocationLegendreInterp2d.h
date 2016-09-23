
#include "LegendreInterp2d.h"
#include "DoubleVectorNd/SCC_DoubleVector2d.h"
#include "GridFunctionNd/SCC_GridFunction2d.h"

#ifndef _FixedLocationLegendreInterp2d_
#define _FixedLocationLegendreInterp2d_


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


class FixedLocationLegendreInterp2d
{
public:

	FixedLocationLegendreInterp2d():interpWeights()
    {
	xMinIndex=0; xMaxIndex=0;
	yMinIndex=0; yMaxIndex=0;
    }

	FixedLocationLegendreInterp2d( int interpolationOrder, double xPos, double yPos,
    const SCC::GridFunction2d& F)
    {
	LegendreInterp2d interp2d(interpolationOrder);
    interp2d.captureNodesAndWeights(xPos,yPos,F, xMinIndex,xMaxIndex, yMinIndex,yMaxIndex, interpWeights);
    }


	FixedLocationLegendreInterp2d(const FixedLocationLegendreInterp2d& P)
    {
	xMinIndex= P.xMinIndex; xMaxIndex= P.xMaxIndex;
	yMinIndex= P.yMinIndex; yMaxIndex= P.yMaxIndex;
	interpWeights.initialize(P.interpWeights);
    }

	void initialize(int interpolationOrder, double xPos, double yPos,
    const SCC::GridFunction2d& F)
    {
	LegendreInterp2d interp2d(interpolationOrder);
    interp2d.captureNodesAndWeights(xPos,yPos,F, xMinIndex,xMaxIndex, yMinIndex,yMaxIndex, interpWeights);
    }


	double getInterpolatedValue(const SCC::GridFunction2d& F)
	{
    long i; long j;
    long fxIndex; long fyIndex;

    double interpolationValue  = 0.0;
    for(i = xMinIndex,fxIndex=0; i <= xMaxIndex; i++,fxIndex++)
    {
    for(j = yMinIndex,fyIndex=0; j <= yMaxIndex; j++,fyIndex++)
    {
    interpolationValue += F(i,j)*interpWeights(fxIndex,fyIndex);
    }}

    return interpolationValue;
	}

	long xMinIndex; long xMaxIndex;
	long yMinIndex; long yMaxIndex;

	SCC::DoubleVector2d interpWeights;
};
#endif
