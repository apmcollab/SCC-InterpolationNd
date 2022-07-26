//
// FourierInterpolation1d.h
//

#include <cmath>


#include "DoubleVectorNd/SCC_DoubleVector1d.h"
#include "GridFunctionNd/SCC_GridFunction1d.h"
#include "FFTW3_InterfaceNd/SCC_fftw3_1d.h"

#ifndef FOURIER_INTERPOLATION_1D_
#define FOURIER_INTERPOLATION_1D_

namespace SCC
{

class FourierInterpolation1d
{
	public :

	FourierInterpolation1d()
	{
	 initialize();
	}

	FourierInterpolation1d(long nxA, long nxB)
	{
	initialize(nxA, nxB);
	}

	~FourierInterpolation1d()
	{}

	void initialize()
	{
	imag_A.initialize(); imag_Ag.initialize();
	realtrans_A.initialize(); imagtrans_A.initialize();

	imag_B.initialize(); imag_Bg.initialize();
	realtrans_B.initialize(); imagtrans_B.initialize();

	 DFT_A.initialize();
	 DFT_B.initialize();

	 nxA  = 0;
	 nxB  = 0;
	}

	void initialize(long nxA, long nxB)
	{
     this->nxA  = nxA;
     this->nxB  = nxB;

	 realtrans_A.initialize(nxA); imagtrans_A.initialize(nxA);
	 realtrans_B.initialize(nxB); imagtrans_B.initialize(nxB);

	 DFT_A.initialize(nxA); DFT_B.initialize(nxB);
	}

    //
    // Interpolate from nxA values to nxB values (nxA < nxB)
    // If the sigma factor flag is set, then multiply the
    // transform coefficients of the input function by
    // Lanczo's sigma factors before carrying out the inverse
    // transform.
    //

	// GridFunciton1d valA : nxA panels
	// GridFunciton1d valB : nxB panels

    void interpolate(SCC::GridFunction1d& valA, SCC::GridFunction1d& valB,bool sigmaFactorFlag = false)
    {
    	if((this->nxA != valA.getXpanelCount())||(this->nxB != valB.getXpanelCount()))
    	{
    	initialize(valA.getXpanelCount(),valB.getXpanelCount());
    	}

    	imag_Ag.initialize(valA.getXpanelCount(),valA.getXmin(),valA.getXmax());
    	imag_Ag.setToValue(0.0);

    	imag_Bg.initialize(valB.getXpanelCount(),valB.getXmin(),valB.getXmax());
    	imag_Bg.setToValue(0.0);


    	DFT_A.fftw1d_forward(valA, imag_Ag, realtrans_A, imagtrans_A);

        long kA_Index;
    	long kB_Index;

    	realtrans_B.setToValue(0.0);
    	imagtrans_B.setToValue(0.0);

    	long nx = nxA;
    	if(nxA > nxB) nx = nxB;

        for(long k1 = -(nx/2); k1 <= (nx-1)/2; k1++)
        {
        		kA_Index = k1 + (nxA/2);
        		kB_Index = k1 + (nxB/2);

        		realtrans_B(kB_Index) = realtrans_A(kA_Index);
        		imagtrans_B(kB_Index) = imagtrans_A(kA_Index);
        }

        //
        // When reducing resolution for an even number of panels, the
        // highest cosine mode coefficient must be scaled by 2.0
        //

        if((nxB%2 == 0)&&(nxA > nxB))
        {
        	realtrans_B(0) *= 2.0;
        	imagtrans_B(0) *= 2.0;
        }

        //
        // To reduce oscillations apply Lanczo's sigma factors.
        //

    	double pi2N = 3.14159265358979323846/(double)(nxB/2); // NB : integer division here
    	double sigma;
    	if(sigmaFactorFlag)
    	{
    		for(long k1 = -(nxB/2); k1 <= (nxB-1)/2; k1++)
    		{
    			if(k1 != 0)
    			{
    			kB_Index = k1 + (nxB/2);
    			sigma = std::sin(std::fabs(k1*pi2N))/std::abs(k1*pi2N);
    	        realtrans_B(kB_Index) *= sigma;
    	        imagtrans_B(kB_Index) *= sigma;
    			}
    		}
    	}


        DFT_B.fftw1d_inverse(realtrans_B,imagtrans_B, valB, imag_Bg);
    }


    //
    // Interpolate from nxA values to nxB values (nxA < nxB)
    // If the sigma factor flag is set, then multiply the
    // transform coefficients of the input function by
    // Lanczo's sigma factors before carrying out the inverse
    // transform.
    //
    void interpolate(SCC::DoubleVector1d& valA, SCC::DoubleVector1d& valB,bool sigmaFactorFlag = false)
    {
    	if((this->nxA != valA.getSize())||(this->nxB != valB.getSize()))
    	{
    	initialize(valA.getSize(),valB.getSize());
    	}

   	   imag_A.initialize(nxA);
       imag_A.setToValue(0.0);


  	   imag_B.initialize(nxB);
  	   imag_B.setToValue(0.0);

    	DFT_A.fftw1d_forward(valA, imag_A, realtrans_A, imagtrans_A);

        long kA_Index;
    	long kB_Index;

    	realtrans_B.setToValue(0.0);
    	imagtrans_B.setToValue(0.0);

    	long nx = nxA;
    	if(nxA > nxB) nx = nxB;

        for(long k1 = -(nx/2); k1 <= (nx-1)/2; k1++)
        {
        		kA_Index = k1 + (nxA/2);
        		kB_Index = k1 + (nxB/2);

        		realtrans_B(kB_Index) = realtrans_A(kA_Index);
        		imagtrans_B(kB_Index) = imagtrans_A(kA_Index);
        }

        //
        // When reducing resolution for an even number of panels, the
        // highest cosine mode coefficient must be scaled by 2.0
        //

        if((nxB%2 == 0)&&(nxA > nxB))
        {
        	realtrans_B(0) *= 2.0;
        	imagtrans_B(0) *= 2.0;
        }

        //
        // To reduce oscillations apply Lanczo's sigma factors.
        //

    	double pi2N = 3.14159265358979323846/(double)(nxB/2); // NB : integer division here
    	double sigma;
    	if(sigmaFactorFlag)
    	{
    		for(long k1 = -(nxB/2); k1 <= (nxB-1)/2; k1++)
    		{
    			if(k1 != 0)
    			{
    			kB_Index = k1 + (nxB/2);
    			sigma = std::sin(std::fabs(k1*pi2N))/std::abs(k1*pi2N);
    	        realtrans_B(kB_Index) *= sigma;
    	        imagtrans_B(kB_Index) *= sigma;
    			}
    		}
    	}


        DFT_B.fftw1d_inverse(realtrans_B,imagtrans_B, valB, imag_B);
    }

    SCC::fftw3_1d DFT_A;
    SCC::fftw3_1d DFT_B;

    long    nxA;
    long    nxB;

    SCC::DoubleVector1d imag_A;
    SCC::GridFunction1d imag_Ag;
    SCC::DoubleVector1d realtrans_A; SCC::DoubleVector1d imagtrans_A;

    SCC::DoubleVector1d imag_B;
    SCC::GridFunction1d imag_Bg;
    SCC::DoubleVector1d realtrans_B; SCC::DoubleVector1d imagtrans_B;
};

} // namespace SCC
#endif
