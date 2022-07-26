/*
 * RefineGrid1d.h
 *
 *  Created on: Jul 31, 2020
 *      Author: anderson
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


