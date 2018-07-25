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

#include "LegendreGridFunApprox1d.h"

extern "C" void legendre_approx_1d(double x, long xdegree, long xpanels, double xmin, double xmax, double* fdataptr,
double* val, double* dval)
{
    // Declared static so the class instance is not instantiated with every call

	static LegendreGridFunApprox1d legendreApprox1d(xdegree,xpanels,xmin,xmax, fdataptr);
 
    // If the degree of approximation has changed induces a re-initialization of the class


    // If the structure of the grid data being interpolated changes, or
    // the degree of approximation changes, then reset values

    if( ((int)xdegree != legendreApprox1d.degreeX)
    ||  ((int)xpanels != legendreApprox1d.xDataPanels)
    ||  (abs(xmin - legendreApprox1d.xDataMin) > 1.0e-12)
    ||  (abs(xmax - legendreApprox1d.xDataMax) > 1.0e-12)
	)
    {
    legendreApprox1d.initialize(xdegree,xpanels,xmin,xmax,fdataptr);
    }

    legendreApprox1d.setDataPtr(fdataptr);
    legendreApprox1d.evaluate(x, val[0], dval[0]);
}

