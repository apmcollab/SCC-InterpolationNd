//
// FourierInterpolation1D.h 
//
// Author: Chris Anderson  
// (C) UCLA 2012 
//
#ifndef _FourierInterpolation1D_
#define _FourierInterpolation1D_
#include "FFTW3_InterfaceNd/SCC_fftw3_1d.h"

class FourierInterpolation1D
{
	public :

	FourierInterpolation1D()
	{
	 real_A = 0; imag_A = 0; real_B = 0; imag_B = 0;
	 sigmaFactorFlag  = false;
	 periodicDataFlag = false;
	 initialize();
	}

	FourierInterpolation1D(long nxA, long nxB)
	{
	real_A = 0; imag_A = 0; real_B = 0; imag_B = 0;
	sigmaFactorFlag  = false;
	periodicDataFlag = false;
	initialize(nxA, nxB);
	}

	~FourierInterpolation1D()
	{
	destroyData();
	}

	void initialize()
	{
	 destroyData();
	 DFT_A.initialize();
	 DFT_B.initialize();
	 nxA  = 0;
	 nxB  = 0;
	 nMax = 0;
	 sigmaFactorFlag  = false;
	 periodicDataFlag = false;
	}

	void initialize(long nxA, long nxB)
	{
	 destroyData();
	 DFT_A.initialize(nxA);
	 DFT_B.initialize(nxB);
	 this->nxA  = nxA;
	 this->nxB  = nxB;
	 this->nMax = (nxA > nxB) ? nxA : nxB;
	 sigmaFactorFlag  = false;
	 periodicDataFlag = false;
	 allocateData();
	}

    //
    // If the periodiceDataFlag is set, then the input/output data
    // is assumed to be one larger than the panel count,
    // and the output of the interpolant has periodicity
    // enforced.
    //
    void setPeriodicDataFlag()
    {
    periodicDataFlag   = true;
    }

    void clearPeriodicDataFlag()
    {
    periodicDataFlag  = false;
    }

    void setSigmaFactorFlag()
    {
    sigmaFactorFlag = true;
    }

    void clearSigmaFactorFlag()
    {
    sigmaFactorFlag = false;
    }

	void allocateData()
	{
	real_A = new double[nMax];
	imag_A = new double[nMax];
	real_B = new double[nMax];
	imag_B = new double[nMax];
	}

    void destroyData()
    {
    if(real_A != 0) delete [] real_A;
    if(imag_A != 0) delete [] imag_A;
    if(real_B != 0) delete [] real_B;
    if(imag_B != 0) delete [] imag_B;
    real_A = 0; imag_A = 0; real_B = 0; imag_B = 0;
    }
    //
    // Interpolate from nxA values to nxB values (nxA < nxB)
    // If the sigma factor flag is set, then multiply the
    // transform coefficients of the input function by
    // Lanczo's sigma factors before carrying out the inverse
    // transform.
    //
    void interpolate(long nxA, double* valA, long nxB, double* valB)
    {

        bool periodicFlagCache = periodicDataFlag;

    	if((this->nxA != nxA)||(this->nxB != nxB))
    	{
    	initialize(nxA,nxB);
    	periodicDataFlag = periodicFlagCache;
    	}

    	for(long i = 0; i < nxA; i++)
    	{
    		imag_A[i] = 0.0;
    	}

    	DFT_A.fftw1d_forward(nxA, valA, imag_A, real_B, imag_B);

        long kA_Index;
    	long kB_Index;
    	double scaleFactor = (double)nxB/(double)nxA;


    	for(long i = 0; i < nxA; i++)
    	{
    	real_B[i] *= scaleFactor;
    	imag_B[i] *= scaleFactor;
    	}

    	for(long i = 0; i < nxB; i++)
    	{
    	real_A[i] = 0.0;
    	imag_A[i] = 0.0;
    	}



    	long nx = nxA;
    	if(nxA > nxB) nx = nxB;

        for(long k1 = -(nx/2); k1 <= (nx-1)/2; k1++)
        {
        		kA_Index = k1 + (nxA/2);
        		kB_Index = k1 + (nxB/2);

        		real_A[kB_Index] = real_B[kA_Index];
        		imag_A[kB_Index] = imag_B[kA_Index];
        }

        //
        // When reducing resolution for an even number of panels, the
        // highest cosine mode coefficient must be scaled by 2.0
        //

        if((nxB%2 == 0)&&(nxA > nxB))
        {
        	real_A[0] *= 2.0;
        	imag_A[0] *= 2.0;
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
    			sigma = sin(fabs(k1*pi2N))/fabs(k1*pi2N);
    	        real_A[kB_Index] *= sigma;
    	        imag_A[kB_Index] *= sigma;
    			}
    		}
    	}


        DFT_B.fftw1d_inverse(nxB, real_A,imag_A, valB, imag_B);

        if(periodicDataFlag) // enforce periodicity
        {
        valB[nxB] = valB[0];
        }
    }

    SCC::fftw3_1d DFT_A;
    SCC::fftw3_1d DFT_B;

    long    nxA;
    long    nxB;
    long   nMax;

    double* real_A;
    double* imag_A;
    double* real_B;
    double* imag_B;

    bool sigmaFactorFlag;
    bool periodicDataFlag;
};

#endif
