/*
 * RefineGrid2d.h
 *
 *  Created on: Jul 31, 2020
 *      Author: anderson
 */

#ifdef _OPENMP
#include <omp.h>
#endif


#ifndef REFINE_GRID_2D_
#define REFINE_GRID_2D_

#include "GridFunctionNd/SCC_GridFunction2d.h"
#include "InterpolationNd/LegendreGridFunApprox2d.h"

class RefineGrid2d
{
	public:

	void refineNX(int degreeX, int degreeY, int N, SCC::GridFunction2d& F, SCC::GridFunction2d& FN)
	{
	    long xPanels = F.getXpanelCount();
	    double xMin  = F.getXmin();
	    double xMax  = F.getXmax();
	    double   hx  = F.getHx();


	    long yPanels = F.getYpanelCount();
	    double yMin  = F.getYmin();
	    double yMax  = F.getYmax();
	    double   hy  = F.getHy();

		double xPos;
		double yPos;

#ifndef _OPENMP
		interp2d.initialize(degreeX,xPanels,xMin,xMax,degreeY,yPanels,yMin,yMax,F.getDataPointer());
		FN.initialize(N*xPanels,xMin,xMax,N*yPanels,yMin,yMax);

		double hFact= 1.0/(double)N;

	    for(long i = 0; i < xPanels; i++)
	    {
	   	for(long j = 0; j < yPanels; j++)
	    {
	    	for(long p = 0; p < N; p++)
	    	{
	        xPos = xMin + i*hx + p*hFact*hx;
	        for(long q = 0; q < N; q++)
	    	{
	    	yPos = yMin + j*hy + q*hFact*hy;
	    	if((p == 0)&&(q == 0))
	    	{
	    		FN(N*i,N*j)   = F(i,j);
	    	}
	    	else
	    	{
	    		FN(N*i+p,N*j+q) = interp2d.evaluate(xPos,yPos);
	    	}

	    	}}
	    }}

	    // Fill in top and right side

	    long j = yPanels;
	    yPos   = yMax;

	    for(long i = 0; i < xPanels; i++)
	    {
	    FN(N*i,N*j)   = F(i,j);
	    for(long p = 1; p < N; p++)
	    {
	        xPos = xMin + i*hx + p*hFact*hx;
	        FN(N*i+p,N*j) = interp2d.evaluate(xPos,yPos);
	    }
	    }

	    long i = xPanels;
	    xPos   = xMax;

	    for(long j = 0; j < yPanels; j++)
	    {
	    FN(N*i,N*j)   = F(i,j);
	    for(long q = 1; q < N; q++)
	    {
	        yPos = yMin + j*hy + q*hFact*hy;
	        FN(N*i,N*j+q) = interp2d.evaluate(xPos,yPos);
	    }
	    }

	    // Upper right corner

	    FN(N*xPanels,N*yPanels) = F(xPanels,yPanels);
#endif

#ifdef _OPENMP

	    int curThreadCount = omp_get_max_threads();
	    int tI;

	    interp2dA.resize(curThreadCount);

	    for(long k = 0; k < curThreadCount; k++)
	    {
	    	interp2dA[k].initialize(degreeX,xPanels,xMin,xMax,degreeY,yPanels,yMin,yMax,F.getDataPointer());
	    }

		FN.initialize(N*xPanels,xMin,xMax,N*yPanels,yMin,yMax);

		double hFact= 1.0/(double)N;

	    #pragma omp parallel for  \
        private(xPos,yPos,tI)\
        schedule(static,1)
	    for(long i = 0; i < xPanels; i++)
	    {
        tI = omp_get_thread_num();
	   	for(long j = 0; j < yPanels; j++)
	    {

	    	for(long p = 0; p < N; p++)
	    	{
	        xPos = xMin + i*hx + p*hFact*hx;
	        for(long q = 0; q < N; q++)
	    	{
	    	yPos = yMin + j*hy + q*hFact*hy;
	    	if((p == 0)&&(q == 0))
	    	{
	    		FN(N*i,N*j)   = F(i,j);
	    	}
	    	else
	    	{
	    		FN(N*i+p,N*j+q) = interp2dA[tI].evaluate(xPos,yPos);
	    	}

	    	}}
	    }}

	    // Fill in top and right side
		long i; long j;

	    j = yPanels;
	    yPos   = yMax;

	    for(i = 0; i < xPanels; i++)
	    {
	    FN(N*i,N*j)   = F(i,j);
	    for(long p = 1; p < N; p++)
	    {
	        xPos = xMin + i*hx + p*hFact*hx;
	        FN(N*i+p,N*j) = interp2dA[0].evaluate(xPos,yPos);
	    }
	    }

	    i = xPanels;
	    xPos   = xMax;

	    for(j = 0; j < yPanels; j++)
	    {
	    FN(N*i,N*j)   = F(i,j);
	    for(long q = 1; q < N; q++)
	    {
	        yPos = yMin + j*hy + q*hFact*hy;
	        FN(N*i,N*j+q) = interp2dA[0].evaluate(xPos,yPos);
	    }
	    }

	    // Upper right corner

	    FN(N*xPanels,N*yPanels) = F(xPanels,yPanels);


#endif
	}

#ifndef _OPENMP
	LegendreGridFunApprox2d interp2d;
#endif

#ifdef _OPENMP
	std::vector<LegendreGridFunApprox2d> interp2dA;
#endif
};
#endif


