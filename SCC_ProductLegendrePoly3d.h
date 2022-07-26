#include "SCC_LegendrePoly.h"

#ifndef PRODUCT_LEGENDRE_POLY_3D_
#define PRODUCT_LEGENDRE_POLY_3D_

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
namespace SCC
{

class ProductLegendrePoly3d
{
public:
	ProductLegendrePoly3d()
	{
	initialize();
	}

	ProductLegendrePoly3d(const ProductLegendrePoly3d& P)
	{
	a = P.a; b = P.b;
	c = P.c; d = P.d;
	e = P.e; f = P.f;

	polyX = P.polyX;
	polyY = P.polyY;
	polyZ = P.polyZ;
    }

	ProductLegendrePoly3d(double xMin, double xMax, double yMin,
			              double yMax, double zMin, double zMax)
	{
	initialize(xMin,xMax,yMin,yMax,zMin,zMax);
	}

    ProductLegendrePoly3d(double xMin, double xMax, int pxIndex, double yMin, double yMax, int pyIndex,
			              double zMin, double zMax, int pzIndex)
	{
	initialize(xMin,xMax,pxIndex,yMin,yMax,pyIndex,zMin,zMax,pzIndex);
	}

    virtual ~ProductLegendrePoly3d(){};

	void initialize()
	{
	a = -1.0; b =  1.0;
	c = -1.0; d =  1.0;
	e = -1.0; f =  1.0;

	polyX.initialize(a,b);
	polyY.initialize(c,d);
	polyZ.initialize(e,f);
	}

	void initialize(double xMin, double xMax, double yMin,
			        double yMax, double zMin, double zMax)
	{
		this->a = xMin; this->b = xMax;
		this->c = yMin; this->d = yMax;
		this->e = zMin; this->f = zMax;
		polyX.initialize(a,b);
		polyY.initialize(c,d);
		polyZ.initialize(e,f);
	}


	void initialize(double xMin, double xMax, int pxIndex, double yMin, double yMax, int pyIndex,
			        double zMin, double zMax, int pzIndex)
	{
		this->a = xMin; this->b = xMax;
		this->c = yMin; this->d = yMax;
		this->e = zMin; this->f = zMax;
		polyX.initialize(a,b,pxIndex);
		polyY.initialize(c,d,pyIndex);
		polyZ.initialize(e,f,pzIndex);
	}

	void setIndices(int pxIndex, int pyIndex, int pzIndex)
	{
		polyX.setIndex(pxIndex);
		polyY.setIndex(pyIndex);
		polyZ.setIndex(pzIndex);
	}

	double evaluate(double xPos, int pxIndex, double yPos, int pyIndex,  double zPos, int pzIndex) const
	{
    return polyX.evaluate(xPos,pxIndex)*polyY.evaluate(yPos,pyIndex)*polyZ.evaluate(zPos,pzIndex);
	}

	double operator()(double xPos, double yPos, double zPos) const
	{
	return polyX(xPos)*polyY(yPos)*polyZ(zPos);
	}

	std::function<double(double, double , double)> getEvaluationPtr() const
	{
	std::function<double(double, double, double)> F = [this](double x, double y, double z) {return this->operator()(x,y,z);};
	return std::move(F);
	}

	double a; double b;
	double c; double d;
	double e; double f;

	LegendrePoly polyX;
	LegendrePoly polyY;
	LegendrePoly polyZ;
};

} // namespace SCC
#endif
