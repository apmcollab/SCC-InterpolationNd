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

#include <vector>
#include <iostream>
using namespace std;

#include "LegendreGridFunApprox2d.h"

extern "C" void legendre_approx_2d(double x, long xdegree, long xpanels, double xmin, double xmax,
                                   double y, long ydegree, long ypanels, double ymin, double ymax, double* fdataptr,
double* val, double* dval)
{
	// Declared static so the class instance is not instantiated with every call

	static LegendreGridFunApprox2d legendreApprox2d(xdegree,xpanels,xmin,xmax,
	                                                ydegree,ypanels,ymin,ymax,fdataptr);

    // If the structure of the grid data being interpolated changes, or
    // the degree of approximation changes then re-initialize.

    if( ((int)xdegree != legendreApprox2d.degreeX)
    ||  ((int)xpanels != legendreApprox2d.xDataPanels)
    ||  (abs(xmin - legendreApprox2d.xDataMin) > 1.0e-12)
    ||  (abs(xmax - legendreApprox2d.xDataMax) > 1.0e-12)
    ||  ((int)ydegree != legendreApprox2d.degreeY)
    ||  ((int)ypanels != legendreApprox2d.yDataPanels)
    ||  (abs(ymin - legendreApprox2d.yDataMin) > 1.0e-12)
    ||  (abs(ymax - legendreApprox2d.yDataMax) > 1.0e-12)
	)
    {
    legendreApprox2d.initialize(xdegree,xpanels,xmin,xmax,
                                ydegree,ypanels,ymin,ymax,fdataptr);
    }


    legendreApprox2d.setDataPtr(fdataptr);
    legendreApprox2d.setFortranFlag();

    vector<double> dfval(2,0.0);
    legendreApprox2d.evaluate(x, y, val[0], dfval);

    dval[0] = dfval[0];
    dval[1] = dfval[1];

}

