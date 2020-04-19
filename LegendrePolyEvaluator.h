/*
 * LegendrePolyEvaluator.h
 *
 *  Created on: Jul 22, 2018
 *      Author: anderson
 */
#include <vector>
#include <functional>
#include <cmath>

#ifndef _LegendrePolyEvaluator_
#define _LegendrePolyEvaluator_

//
// Legendre polynomials
//
// P_[0](x) = 1
// P_[1](x) = x
//
// (n+1)*P_[n+1](x) - (2n+1)*x*P_[n](x) + n*P_[n-1] = 0
//
// (P_[n], P_[n]) = 2/(2*n+1)
//
//
// When using the mapping from [a,b] to [-1,1]
//
// x = (2*s - (a+b))/(b-a)
//
// Normalization factor = 1.0/(sqrt((b-a)/(2.0*k + 1.0)));
// ===================================================
//
// Recurrence for un-normalized Legendre polynomials
//
// (n+1)*P_[n+1](x) - (2n+1)*x*P_[n](x) + n*P_[n-1] = 0
//
// Recurrence for un-normalized Legendre polynomials derivatives
//
// P'_[n+1](x) - P'_[n-1](x) = (2*n +1)*P_[n](x)
//
//
// This class provides member functions that return
//
// (a) arrays of normalized Legendre polynomial values
// (b) arrays of the derivatives of the normalized Legendre polynomial values
//
// up to a maximal specified index.
//
// The Legendre polynomials are those scaled to an interval [a,b] with
// a multiplicative normalization factor so they are orthonormal over
// [a,b].
//
// The default interval is the standard interval of definition [-1,1].
//
class LegendrePolyEvaluator
{
    public:

	LegendrePolyEvaluator()
	{
	a      = -1.0;
	b      =  1.0;
	maxIndex  = -1;
	}

	LegendrePolyEvaluator(const LegendrePolyEvaluator& P)
	{
	a         =  P.a;
	b         =  P.b;
	maxIndex  =  P.maxIndex;
	}

	LegendrePolyEvaluator(double a, double b, int maxIndex)
	{
		this->a      = a;
		this->b      = b;
		this->maxIndex  = maxIndex;
	}

	virtual ~LegendrePolyEvaluator(){};

	void initialize()
	{
	a         =  -1.0;
	b         =   1.0;
	maxIndex  =    -1;
	}

	void initialize(double a, double b, int maxIndex)
	{
		this->a      = a;
		this->b      = b;
		this->maxIndex  = maxIndex;
	}

	void setMaxIndex(int maxIndex)
	{
		this->maxIndex = maxIndex;
	}

	void evaluate(double s,std::vector<double>& values)
	{
    double x = (2*s - (a+b))/(b-a);
    double normFactor;

    values.resize(maxIndex+1);

    double Fkm1 = 1;
	double Fk   = x;
    double Fkp1;

    normFactor = 1.0/(sqrt((b-a)));
	if(maxIndex >= 0) values[0] = Fkm1*normFactor;

	normFactor = 1.0/(sqrt((b-a)/(3.0)));
	if(maxIndex >= 1) values[1] = Fk*normFactor;

    for(int k = 1; k <= maxIndex-1; k++)
    {
    	Fkp1 =  ((2.0*k+1.0)*x*Fk - k*Fkm1)/(k+1.0);
    	Fkm1 = Fk;
    	Fk   = Fkp1;
    	normFactor = 1.0/(sqrt((b-a)/(2.0*(k+1) + 1.0)));
    	values[k+1] = normFactor*Fk;
    }
	}

	void evaluateDerivatives(double s,std::vector<double>& derivativeValues)
	{
	//
	// Transform to unit interval
	//
	double x = (2*s - (a+b))/(b-a);

    derivativeValues.resize(maxIndex+1);

    double Fkm1 = 1;
	double Fk   = x;
    double Fkp1;

    // Create un-normalized values

	if(maxIndex >= 0)
	{
	derivativeValues[0] = 0.0;
	}

	if(maxIndex >= 1)
	{
	derivativeValues[1] = 1.0;
	}

    for(int k = 1; k <= maxIndex-1; k++)
    {
    	Fkp1 =  ((2.0*k+1.0)*x*Fk - k*Fkm1)/(k+1.0);
    	Fkm1 = Fk;
    	Fk   = Fkp1;
    	derivativeValues[k+1] = derivativeValues[k-1] + (2.0*(k) + 1.0)*Fkm1;
    }

    // Add normalization and transformation scaling factor

    for(int k = 0; k <= maxIndex; k++)
    {
    derivativeValues[k] *= (1.0/(sqrt((b-a)/(2.0*k + 1.0))))*(2.0/(b-a));
    }
	}

	double              a;
	double              b;
	int          maxIndex;
	std::vector<double> values;
};



#endif /* _LegendrePolyEvaluator_ */
