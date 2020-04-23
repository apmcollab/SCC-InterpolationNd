#ifndef _LegendrePoly_
#define _LegendrePoly_
#include <cmath>
#include <functional>
//
// P_[0](x) = 1
// P_[1](x) = x
//
// (n+1)*P_[n+1](x) - (2n+1)*x*P_[n](x) + n*P_[n-1] = 0
//
// (P_[n], P_[n]) = 2/(2*n+1)
//
// When using the mapping from [a,b] to [-1,1]
//
// x = (2*s - (a+b))/(b-a)
//
// (P_[n](s), P_[n](s)) = ((b-a)/2)*(2/(2*n+1)) = (b-a)/(2*n +1)
//

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


class LegendrePoly
{
public:
	LegendrePoly()
	{
	a      = -1.0;
	b      =  1.0;
	index  =  0;
	}

	LegendrePoly(const LegendrePoly& P)
	{
	a      =  P.a;
	b      =  P.b;
	index  =  P.index;
	}

	LegendrePoly(double a, double b)
	{
		this->a      = a;
		this->b      = b;
		this->index  = 0;
	}

	LegendrePoly(double a, double b, int index)
	{
		this->a      = a;
		this->b      = b;
		this->index  = index;
	}

	virtual ~LegendrePoly(){};

	void initialize()
	{
	a      = -1.0;
	b      =  1.0;
	index  =    0;
	}

	void initialize(double a, double b, int index)
	{
		this->a      = a;
		this->b      = b;
		this->index  = index;
	}

	void setIndex(int index)
	{
		this->index = index;
	}

	void initialize(double a, double b)
	{
		this->a = a;
		this->b = b;
	}

	inline double operator()(double s) const
	{
    return evaluate(s,index);
	}


	std::function<double(double)> getEvaluationPtr() const
	{
	std::function<double(double)> F = [this](double x) {return this->operator()(x);};
	return std::move(F);
	}
 

	double evaluate(double s, int index) const
	{
	double x          = (2*s - (a+b))/(b-a);
	double normFactor = 1.0/(std::sqrt((b-a)/(2.0*index +1.0)));
    double Fkm1;
	double Fk;
    double Fkp1;

	Fkm1 = 1.0;
	if(index == 0) return Fkm1*normFactor;
	Fk   = x;
	if(index == 1) return Fk*normFactor;


	long k;
	for(k = 1; k <= index-1; k++)
	{
	Fkp1 =  ((2.0*k+1.0)*x*Fk - k*Fkm1)/(k+1.0);
	Fkm1 = Fk;
	Fk   = Fkp1;
	}

	return Fk*normFactor;
	}

	double   a;
	double   b;
	int  index;
};
#endif

