/*
 * RefineGrid1d.h
 *
 *  Created on: Jul 31, 2020
 *      Author: anderson
 */
/*
#############################################################################
#
# Copyright  2020- Chris Anderson
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
#ifndef REFINE_GRID_1D_
#define REFINE_GRID_1D_

#include "GridFunctionNd/SCC_GridFunction1d.h"
#include "SCC_LegendreGridFunApprox1d.h"

namespace SCC
{

class RefineGrid1d
{
	public:


    RefineGrid1d()
	{
	initialize();
	}

    RefineGrid1d(const RefineGrid1d& R)
	{
	interp1d.initialize(R.interp1d);
	}

	void initialize()
	{
	interp1d.initialize();
	}

	void refineNX(int degree, int N, SCC::GridFunction1d& F, SCC::GridFunction1d& FN)
	{
	    long xPanels = F.getXpanelCount();
	    double xMin  = F.getXmin();
	    double xMax  = F.getXmax();
	    double   hx  = F.getHx();

		interp1d.initialize(degree,xPanels,xMin,xMax,F.getDataPointer());

		FN.initialize(N*xPanels,xMin,xMax);

		double xPos;
		double hFact = 1.0/(double)N;

	    for(long i = 0; i < xPanels; i++)
	    {
	    	FN(N*i)   = F(i);

	    	for(long k = 1; k < N; k++)
	    	{
	    	xPos = xMin + i*hx + k*hFact*hx;
	    	FN(N*i+k) = interp1d.evaluate(xPos);
	    	}
	    }
	    FN(N*xPanels) = F(xPanels);
	}


	LegendreGridFunApprox1d interp1d;
};

} // namespace SCC
#endif


