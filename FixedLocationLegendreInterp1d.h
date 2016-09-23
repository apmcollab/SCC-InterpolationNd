
#include "LegendreInterp1d.h"
#include "DoubleVectorNd/UCLAQ_DoubleVector1d.h"
#include "GridFunctionNd/UCLAQ_GridFunction1d.h"

#ifndef _FixedLocationLegendreInterp1d_
#define _FixedLocationLegendreInterp1d_


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


class FixedLocationLegendreInterp1d
{
public:


	FixedLocationLegendreInterp1d():interpWeights()
    {
	xMinIndex=0; xMaxIndex=0;
    }

	FixedLocationLegendreInterp1d(int interpolationOrder, double xPos, const UCLAQ::GridFunction1d& F)
    {
	LegendreInterp1d interp1d(interpolationOrder);
    interp1d.captureNodesAndWeights(xPos,F,xMinIndex,xMaxIndex,interpWeights);
    }


	void initialize(int interpolationOrder, double xPos, const UCLAQ::GridFunction1d& F)
    {
	LegendreInterp1d interp1d(interpolationOrder);
    interp1d.captureNodesAndWeights(xPos,F, xMinIndex,xMaxIndex,interpWeights);
    }


	FixedLocationLegendreInterp1d(const FixedLocationLegendreInterp1d& P)
    {
	xMinIndex= P.xMinIndex; xMaxIndex= P.xMaxIndex;
	interpWeights.initialize(P.interpWeights);
    }

	double getInterpolatedValue(const UCLAQ::GridFunction1d& F)
	{
    long i;
    long fxIndex;

    double interpolationValue  = 0.0;
    for(i = xMinIndex,fxIndex=0; i <= xMaxIndex; i++,fxIndex++)
    {
    interpolationValue += F(i)*interpWeights(fxIndex);
    }

    return interpolationValue;
	}

	long xMinIndex; long xMaxIndex;

	UCLAQ::DoubleVector1d interpWeights;
};
#endif
