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


#include "LegendreGridFunApprox3d.h"

extern "C" void legendre_approx_3d(double x, long xdegree, long xpanels, double xmin, double xmax,
                                   double y, long ydegree, long ypanels, double ymin, double ymax,
                                   double z, long zdegree, long zpanels, double zmin, double zmax,
                                   double* fdataptr, double* val, double* dval)
{
	// Declared static so the class instance is not instantiated with every call

	static LegendreGridFunApprox3d legendreApprox3d(xdegree,xpanels,xmin,xmax,
	                                                ydegree,ypanels,ymin,ymax,
	                                                zdegree,zpanels,zmin,zmax,fdataptr);

    // If the structure of the grid data being interpolated changes, or
    // the degree of approximation changes, then re-initialize.

    if( ((int)xdegree != legendreApprox3d.degreeX)
    ||  ((int)xpanels != legendreApprox3d.xDataPanels)
    ||  (abs(xmin - legendreApprox3d.xDataMin) > 1.0e-12)
    ||  (abs(xmax - legendreApprox3d.xDataMax) > 1.0e-12)
    ||  ((int)ydegree != legendreApprox3d.degreeY)
    ||  ((int)ypanels != legendreApprox3d.yDataPanels)
    ||  (abs(ymin - legendreApprox3d.yDataMin) > 1.0e-12)
    ||  (abs(ymax - legendreApprox3d.yDataMax) > 1.0e-12)
    ||  ((int)zdegree != legendreApprox3d.degreeZ)
    ||  ((int)zpanels != legendreApprox3d.zDataPanels)
    ||  (abs(zmin - legendreApprox3d.zDataMin) > 1.0e-12)
    ||  (abs(zmax - legendreApprox3d.zDataMax) > 1.0e-12)
	)
    {
    legendreApprox3d.initialize(xdegree,xpanels,xmin,xmax,
                                ydegree,ypanels,ymin,ymax,
                                zdegree,zpanels,zmin,zmax,fdataptr);
    }


    legendreApprox3d.setDataPtr(fdataptr);
    legendreApprox3d.setFortranFlag();

    vector<double> dfval(3,0.0);

    legendreApprox3d.evaluate(x, y, z, val[0], dfval);

    dval[0] = dfval[0];
    dval[1] = dfval[1];
    dval[2] = dfval[2];
}

