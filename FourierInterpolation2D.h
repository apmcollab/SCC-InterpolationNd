#ifndef _FourierInterpolation2D_
#define _FourierInterpolation2D_
#include "FFTW3_InterfaceNd/SCC_fftw3_2d.h"
#include "SCC_fftw3_2d.h"
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
class FourierInterpolation2D
{
	public :

	FourierInterpolation2D()
	{
	 real_A = 0; imag_A = 0; real_B = 0; imag_B = 0;
	 sigmaFactorFlag  = false;
	 periodicDataFlag = false;
	 initialize();
	}

    //
    // nxA, nxB : The number of panels associated with the input data in the x-direction
    // nyA, nyB : The number of panels associated with the input data in the y-direction
    //

	FourierInterpolation2D(long nxA, long nyA, long nxB, long nyB)
	{
	real_A = 0; imag_A = 0; real_B = 0; imag_B = 0;
	sigmaFactorFlag  = false;
	periodicDataFlag = false;
	initialize(nxA,nyA,nxB,nyB);
	}

	~FourierInterpolation2D()
	{
	destroyData();
	}

	void initialize()
	{
	 destroyData();
	 DFT_A.initialize();
	 DFT_B.initialize();
	 nxA  = 0; nyA = 0;
	 nxB  = 0; nyB = 0;
	 nMax = 0;
	 sigmaFactorFlag   = false;
	 periodicDataFlag  = false;
	}

	void initialize(long nxA, long nyA, long nxB, long nyB)
	{
	 destroyData();
	 DFT_A.initialize(nxA,nyA);
	 DFT_B.initialize(nxB,nyB);
	 this->nxA  = nxA;
	 this->nyA  = nyA;
	 this->nxB  = nxB;
	 this->nyB  = nyB;
	 this->nMax = (nxA*nyA > nxB*nyB) ? nxA*nyA : nxB*nyB;
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
    // Interpolate from nxA X nyA  values to nxB X nyB values
    // If the sigma factor flag is set, then multiply the
    // transform coefficients of the input function by
    // Lanczo's sigma factors before carrying out the inverse
    // transform.
    //
    void interpolate(long nxA, long nyA, double* valA, long nxB, long nyB, double* valB)
    {
        bool columnStorageFlag = false;
        bool periodicFlagCache = periodicDataFlag;

    	if((this->nxA != nxA)||(this->nxB != nxB)||(this->nyA != nyA)||(this->nyB != nyB))
    	{
    	initialize(nxA,nyA,nxB,nyB);
    	periodicDataFlag = periodicFlagCache;
    	}

    	for(long i = 0; i < nxA*nyA; i++)
    	{
    		imag_A[i] = 0.0;
    	}

        long nyAp1 = nyA+1;
        long nyBp1 = nyB+1;

        if(periodicDataFlag)
        {
        	for(long i = 0; i < nxA; i++)
        	{
        	for(long j = 0; j < nyA; j++)
        	{
        		real_A[j + i*nyA] = valA[j + i*nyAp1];
        	}}
        	DFT_A.fftw2d_forward(nxA, nyA, real_A, imag_A, real_B, imag_B, columnStorageFlag);
        }
        else
        {
        	DFT_A.fftw2d_forward(nxA, nyA, valA, imag_A, real_B, imag_B, columnStorageFlag);
    	}

    	double scaleFactor = ((double)nxB * (double)nyB)/((double)nxA * (double)nyA);


    	for(long i = 0; i < nxA*nyA; i++)
    	{
    	real_B[i] *= scaleFactor;
    	imag_B[i] *= scaleFactor;
    	}

    	for(long i = 0; i < nxB*nyB; i++)
    	{
    	real_A[i] = 0.0;
    	imag_A[i] = 0.0;
    	}

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

        		real_A[kBY_Index + kBX_Index*nyB] = real_B[kAY_Index + kAX_Index*nyA];
        		imag_A[kBY_Index + kBX_Index*nyB] = imag_B[kAY_Index + kAX_Index*nyA];
        }}

        //
        // When reducing resolution for an even number of panels, the
        // highest cosine mode coefficient must be scaled by 2.0
        //

	    if((nxB%2 == 0)&&(nxA > nxB))
        {
        for(long k2 = -(ny/2); k2 <= (ny-1)/2; k2++)
    	{
        		kBY_Index = k2 + (nyB/2);
        		real_A[kBY_Index] *= 2.0;
        		imag_A[kBY_Index] *= 2.0;
        }}


	    if((nyB%2 == 0)&&(nyA > nyB))
        {
        for(long k1 = -(nx/2); k1 <= (nx-1)/2; k1++)
    	{
    	        kBX_Index = k1 + (nxB/2);
        		real_A[kBX_Index*nyB] *= 2.0;
        		imag_A[kBX_Index*nyB] *= 2.0;
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
    			sigma = (sin(fabs(k1*pi2N_X))/fabs(k1*pi2N_X))*(sin(fabs(k2*pi2N_Y))/fabs(k2*pi2N_Y));
    	        real_A[kY_Index + kX_Index*nyB] *= sigma;
    	        imag_A[kY_Index + kX_Index*nyB] *= sigma;
    			}
    		}}
    	}


        DFT_B.fftw2d_inverse(nxB, nyB, real_A, imag_A, valB, imag_B);

        if(periodicDataFlag)
        {
        	DFT_B.fftw2d_inverse(nxB, nyB, real_A, imag_A, real_B, imag_B);

        	for(long i = 0; i < nxB; i++)
        	{
        	for(long j = 0; j < nyB; j++)
        	{
        		valB[j + i*nyBp1] = real_B[j + i*nyB];
        	}}

            for(long i = 0; i < nxB; i++)
        	{
            valB[nyB + i*nyBp1] = valB[i*nyBp1];
            }

        	for(long j = 0; j < nyB; j++)
        	{
        	valB[j + nxB*nyBp1] = valB[j];
        	}
        	valB[nyB + nxB*nyBp1] = valB[0];
        }
        else
        {
        	DFT_B.fftw2d_inverse(nxB, nyB, real_A, imag_A, valB, imag_B);
    	}
    }


    SCC::fftw3_2d DFT_A;
    SCC::fftw3_2d DFT_B;

    long    nxA; long nyA;
    long    nxB; long nyB;
    long   nMax;

    double* real_A;
    double* imag_A;
    double* real_B;
    double* imag_B;

    bool sigmaFactorFlag;
    bool periodicDataFlag;
};

#endif
