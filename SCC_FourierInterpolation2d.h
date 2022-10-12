
#include <cmath>


#include "DoubleVectorNd/SCC_DoubleVector2d.h"
#include "GridFunctionNd/SCC_GridFunction2d.h"
#include "FFTW3_InterfaceNd/SCC_fftw3_2d.h"

#ifndef FOURIER_INTERPOLATION_2D_
#define FOURIER_INTERPOLATION_2D_

//
// This class implements Fourier interpolation on double* arrays holding 2D data
// that is stored by ROWS.
//
// The Fourier interpolation procedure consists of transforming the input data,
// and then either padding the output transform with zeros in the case
// that output sampling is greater than the input sampling or truncating the
// transform values in the case the output sampling is less than the input sampling.
//
// Optionally the use of Lanczo's sigma factors can be applied to the output
// transform values.
//
//  Common parameters:
//  nxA, nxB : The number of panels associated with the input data in the x-direction
//  nxA, nxB : The number of panels associated with the input data in the y-direction
//
//
// If periodicDataFlag is set, then the input data is a 2D array in which the
// periodic values are present. The panel counts (nxA, nyA) and (nxB,nyB) are
// the same but the input/output array sizes are assumed to be one larger in
// in each direction. Periodicity of the output values is enforced when this flag
// is set.
//
// Dependencies : SCC::fftw3_2d to provide the requisite transforms.
//
// For efficiency, the values of nxA,nyA, nxB and nyB should be products of small primes
// less than or equal to 7.
//
// (C) Chris Anderson March 12, 2012
//
//

/*
#############################################################################
#
# Copyright  2022- Chris Anderson
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

class FourierInterpolation2d
{
	public :

	FourierInterpolation2d()
	{
	 initialize();
	}

    //
    // nxA, nxB : The number of panels associated with the input data in the x-direction
    // nyA, nyB : The number of panels associated with the input data in the y-direction
    //

	FourierInterpolation2d(long nxA, long nyA, long nxB, long nyB)
	{
	initialize(nxA,nyA,nxB,nyB);
	}

	~FourierInterpolation2d()
	{}

	void initialize()
	{
	imag_A.initialize(); imag_Ag.initialize();
	realtrans_A.initialize(); imagtrans_A.initialize();

	imag_B.initialize(); imag_Bg.initialize();
	realtrans_B.initialize(); imagtrans_B.initialize();

	DFT_A.initialize();
	DFT_B.initialize();

	nxA  = 0; nyA = 0;
	nxB  = 0; nyB = 0;
	}

	void initialize(long nxA, long nyA, long nxB, long nyB)
	{
	 DFT_A.initialize(nxA,nyA);
	 DFT_B.initialize(nxB,nyB);

	 this->nxA  = nxA; this->nyA  = nyA;
	 this->nxB  = nxB; this->nyB  = nyB;

	 realtrans_A.initialize(nxA,nyA); imagtrans_A.initialize(nxA,nyA);
	 realtrans_B.initialize(nxB,nyB); imagtrans_B.initialize(nxB,nyB);
	}


    //
    // Interpolate from nxA X nyA  values to nxB X nyB values
    // If the sigma factor flag is set, then multiply the
    // transform coefficients of the input function by
    // Lanczo's sigma factors before carrying out the inverse
    // transform.
    //
    void interpolate(SCC::GridFunction2d& valA, SCC::GridFunction2d& valB,bool sigmaFactorFlag = false)
    {


    	if((this->nxA != valA.getXpanelCount())||
    	   (this->nxB != valB.getXpanelCount())||
		   (this->nyA != valA.getYpanelCount())||
		   (this->nyB != valB.getYpanelCount()))
    	{
    	initialize(valA.getXpanelCount(),valA.getYpanelCount(),valB.getXpanelCount(),valB.getYpanelCount());
    	}


    	/* Code to check DoubleVector2D version */
        /*
    	SCC::DoubleVector2d vA(nxA,nyA);
    	SCC::DoubleVector2d vB(nxB,nyB);
    	for(long i = 0; i < nxA; i++)
    	{
    	for(long j = 0; j < nyA; j++)
    	{
    		vA(i,j) = valA(i,j);
    	}}

    	interpolate(vA,vB, sigmaFactorFlag);

    	for(long i = 0; i < nxB; i++)
    	{
    	for(long j = 0; j < nyB; j++)
    	{
    		valB(i,j) = vB(i,j);
    	}}
    	valB.enforcePeriodicity();

    	return;
    	*/

     	imag_Ag.initialize(valA.getXpanelCount(),valA.getXmin(),valA.getXmax(),valA.getYpanelCount(),valA.getYmin(),valA.getYmax());
        imag_Ag.setToValue(0.0);

        imag_Bg.initialize(valB.getXpanelCount(),valB.getXmin(),valB.getXmax(),valB.getYpanelCount(),valB.getYmin(),valB.getYmax());
        imag_Bg.setToValue(0.0);


        DFT_A.fftw2d_forward(valA, imag_Ag, realtrans_A, imagtrans_A);


        realtrans_B.setToValue(0.0);
        imagtrans_B.setToValue(0.0);

    	long kAX_Index; long kAY_Index;
    	long kBX_Index; long kBY_Index;

    	long nx = nxA;
    	if(nxA > nxB) nx = nxB;

    	long ny = nyA;
    	if(nyA > nyB) ny = nyB;

        for(long k1 = -(nx/2); k1 <= (nx-1)/2; k1++)
        {
        for(long k2 = -(ny/2); k2 <= (ny-1)/2; k2++)
    	{
        		kAX_Index = k1 + (nxA/2);
        		kAY_Index = k2 + (nyA/2);

        		kBX_Index = k1 + (nxB/2);
        		kBY_Index = k2 + (nyB/2);

        		realtrans_B(kBX_Index,kBY_Index) = realtrans_A(kAX_Index,kAY_Index);
        		imagtrans_B(kBX_Index,kBY_Index) = imagtrans_A(kAX_Index,kAY_Index);
        }}

        // [kY_Index + kX_Index*nyB]
        //
        // When reducing resolution for an even number of panels, the
        // highest cosine mode coefficient must be scaled by 2.0
        //

	    if((nxB%2 == 0)&&(nxA > nxB))
        {
	    kBX_Index = 0;
        for(long k2 = -(ny/2); k2 <= (ny-1)/2; k2++)
    	{
        		kBY_Index = k2 + (nyB/2);
        		realtrans_B(kBX_Index,kBY_Index) *= 2.0; // [kBY_Index]
        		imagtrans_B(kBX_Index,kBY_Index) *= 2.0;
        }}


	    if((nyB%2 == 0)&&(nyA > nyB))
        {
	    kBY_Index = 0;
        for(long k1 = -(nx/2); k1 <= (nx-1)/2; k1++)
    	{
    	        kBX_Index = k1 + (nxB/2);
        		realtrans_B(kBX_Index,kBY_Index) *= 2.0; // [kBX_Index*nyB]
        		imagtrans_B(kBX_Index,kBY_Index) *= 2.0;
        }}

        // Apply sigma factors to the transform coefficients of the output function

    	double pi2N_X = 3.14159265358979323846/(double)(nxB/2); // NB : integer division here
    	double pi2N_Y = 3.14159265358979323846/(double)(nyB/2); // NB : integer division here
    	long kX_Index; long kY_Index;
    	double sigma;
    	if(sigmaFactorFlag)
    	{
    		for(long k1 = -(nxB/2); k1 <= (nxB-1)/2; k1++)
    		{
    	    for(long k2 = -(nyB/2); k2 <= (nyB-1)/2; k2++)
    		{
    			if((k1 != 0)&&(k2 != 0))
    			{
    			kX_Index = k1 + (nxB/2);
    			kY_Index = k2 + (nyB/2);
    			sigma = (std::sin(std::abs(k1*pi2N_X))/std::fabs(k1*pi2N_X))*(std::sin(std::abs(k2*pi2N_Y))/std::abs(k2*pi2N_Y));
    	        realtrans_B(kX_Index,kY_Index)  *= sigma;
    	        imagtrans_B(kX_Index,kY_Index)  *= sigma;
    			}
    		}}
    	}

    	DFT_B.fftw2d_inverse(realtrans_B,imagtrans_B, valB, imag_Bg);
    }



    //
    // Interpolate from nxA X nyA  values to nxB X nyB values
    // If the sigma factor flag is set, then multiply the
    // transform coefficients of the input function by
    // Lanczo's sigma factors before carrying out the inverse
    // transform.
    //
    void interpolate(SCC::DoubleVector2d& valA, SCC::DoubleVector2d& valB,bool sigmaFactorFlag = false)
    {
    	if((this->nxA != valA.getIndex1Size())||
    	   (this->nxB != valB.getIndex1Size())||
		   (this->nyA != valA.getIndex2Size())||
		   (this->nyB != valB.getIndex2Size()))
    	{
    		initialize(valA.getIndex1Size(),valA.getIndex2Size(),valB.getIndex1Size(),valB.getIndex2Size());
    	}

        imag_A.initialize(nxA,nyA);
        imag_A.setToValue(0.0);


   	    imag_B.initialize(nxB,nyB);
   	    imag_B.setToValue(0.0);

        DFT_A.fftw2d_forward(valA, imag_A, realtrans_A, imagtrans_A);

        realtrans_B.setToValue(0.0);
        imagtrans_B.setToValue(0.0);

    	long kAX_Index; long kAY_Index;
    	long kBX_Index; long kBY_Index;

    	long nx = nxA;
    	if(nxA > nxB) nx = nxB;

    	long ny = nyA;
    	if(nyA > nyB) ny = nyB;

        for(long k1 = -(nx/2); k1 <= (nx-1)/2; k1++)
        {
        for(long k2 = -(ny/2); k2 <= (ny-1)/2; k2++)
    	{
        		kAX_Index = k1 + (nxA/2);
        		kAY_Index = k2 + (nyA/2);

        		kBX_Index = k1 + (nxB/2);
        		kBY_Index = k2 + (nyB/2);

        		realtrans_B(kBX_Index,kBY_Index) = realtrans_A(kAX_Index,kAY_Index);
        		imagtrans_B(kBX_Index,kBY_Index) = imagtrans_A(kAX_Index,kAY_Index);
        }}


        // [kY_Index + kX_Index*nyB]
        //
        // When reducing resolution for an even number of panels, the
        // highest cosine mode coefficient must be scaled by 2.0
        //

	    if((nxB%2 == 0)&&(nxA > nxB))
        {
	    kBX_Index = 0;
        for(long k2 = -(ny/2); k2 <= (ny-1)/2; k2++)
    	{
        		kBY_Index = k2 + (nyB/2);
        		realtrans_B(kBX_Index,kBY_Index) *= 2.0; // [kBY_Index]
        		imagtrans_B(kBX_Index,kBY_Index) *= 2.0;
        }}


	    if((nyB%2 == 0)&&(nyA > nyB))
        {
	    kBY_Index = 0;
        for(long k1 = -(nx/2); k1 <= (nx-1)/2; k1++)
    	{
    	        kBX_Index = k1 + (nxB/2);
        		realtrans_B(kBX_Index,kBY_Index) *= 2.0; // [kBX_Index*nyB]
        		imagtrans_B(kBX_Index,kBY_Index) *= 2.0;
        }}

        // Apply sigma factors to the transform coefficients of the output function

    	double pi2N_X = 3.14159265358979323846/(double)(nxB/2); // NB : integer division here
    	double pi2N_Y = 3.14159265358979323846/(double)(nyB/2); // NB : integer division here
    	long kX_Index; long kY_Index;
    	double sigma;
    	if(sigmaFactorFlag)
    	{
    		for(long k1 = -(nxB/2); k1 <= (nxB-1)/2; k1++)
    		{
    	    for(long k2 = -(nyB/2); k2 <= (nyB-1)/2; k2++)
    		{
    			if((k1 != 0)&&(k2 != 0))
    			{
    			kX_Index = k1 + (nxB/2);
    			kY_Index = k2 + (nyB/2);
    			sigma = (std::sin(std::abs(k1*pi2N_X))/std::fabs(k1*pi2N_X))*(std::sin(std::abs(k2*pi2N_Y))/std::abs(k2*pi2N_Y));
    	        realtrans_B(kX_Index,kY_Index)  *= sigma;
    	        imagtrans_B(kX_Index,kY_Index)  *= sigma;
    			}
    		}}
    	}


    	DFT_B.fftw2d_inverse(realtrans_B,imagtrans_B, valB, imag_B);
    }


    SCC::fftw3_2d DFT_A;
    SCC::fftw3_2d DFT_B;

    long    nxA; long nyA;
    long    nxB; long nyB;

    SCC::DoubleVector2d imag_A;
    SCC::GridFunction2d imag_Ag;
    SCC::DoubleVector2d realtrans_A; SCC::DoubleVector2d imagtrans_A;

    SCC::DoubleVector2d imag_B;
    SCC::GridFunction2d imag_Bg;
    SCC::DoubleVector2d realtrans_B; SCC::DoubleVector2d imagtrans_B;
};

} // namespace SCC
#endif
