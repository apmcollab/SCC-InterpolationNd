#ifndef _ProductLegendrePoly2d_
#define _ProductLegendrePoly2d_

#include "LegendrePoly.h"

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


class ProductLegendrePoly2d
{
public:

	ProductLegendrePoly2d()
	{
	initialize();
	}

	ProductLegendrePoly2d(const ProductLegendrePoly2d& P)
	{
	a = P.a; b = P.b;
	c = P.c; d = P.d;

	polyX = P.polyX;
	polyY = P.polyY;
    }

	ProductLegendrePoly2d(double xMin, double xMax, double yMin, double yMax)
	{
	initialize(xMin,xMax,yMin,yMax);
	}

    ProductLegendrePoly2d(double xMin, double xMax, int pxIndex,
    		              double yMin, double yMax, int pyIndex)
	{
	initialize(xMin,xMax,pxIndex,yMin,yMax,pyIndex);
	}

    virtual ~ProductLegendrePoly2d(){};

	void initialize()
	{
	a = -1.0; b =  1.0;
	c = -1.0; d =  1.0;

	polyX.initialize(a,b);
	polyY.initialize(c,d);
	}

	void initialize(double xMin, double xMax, double yMin, double yMax)
	{
		this->a = xMin; this->b = xMax;
		this->c = yMin; this->d = yMax;
		polyX.initialize(a,b);
		polyY.initialize(c,d);
	}


	void initialize(double xMin, double xMax, int pxIndex,
			        double yMin, double yMax, int pyIndex)
	{
		this->a = xMin; this->b = xMax;
		this->c = yMin; this->d = yMax;
		polyX.initialize(a,b,pxIndex);
		polyY.initialize(c,d,pyIndex);
	}

	void setIndices(int pxIndex, int pyIndex)
	{
		polyX.setIndex(pxIndex);
		polyY.setIndex(pyIndex);
	}

	double evaluate(double xPos, int pxIndex, double yPos, int pyIndex) const
	{
    return polyX.evaluate(xPos,pxIndex)*polyY.evaluate(yPos,pyIndex);
	}

	double operator()(double xPos, double yPos) const
	{
	return polyX(xPos)*polyY(yPos);
	}

	std::function<double(double, double)> getEvaluationPtr() const
	{
	std::function<double(double, double)> F = [this](double x, double y) {return this->operator()(x,y);};
	return std::move(F);
	}

	double a; double b;
	double c; double d;

	LegendrePoly polyX;
	LegendrePoly polyY;
};


#endif
