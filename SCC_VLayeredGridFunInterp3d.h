/*
 * VLayeredGridFunInterp3d.h
 *
 *  Created on: Jul 25, 2022
 *      Author: anderson
 *
 *
 *  ====> alpha version <===
 *
 * This class provides member functions that interpolate values from a source
 * SCC::VLayeredGridFun3d instance to a SCC::VLayeredGridFun3d target instance.
 *
 * The method of interpolation is a tensor product of Legendre functions that
 * interpolate the source data values. The data values used to form the interpolant
 * are chosen so that the evaluation is as centered as much as possible with respect
 * to the source value locations. The interpolation formula also never crosses
 * layer interface boundaries.
 *
 * The XYZ interpolation is a full three dimensional interpolation. If the grid
 * points of the target SCC::VLayeredGridFun3d lie outside the domain of the source domain
 * then the approximation is an extrapolation. This is not result in an error condition
 * being thrown, but generally results in a very inaccurate approximation.
 *
 * The XY interpolation routine is a full three dimensional interpolation but is formed
 * using only interpolation in the x-y direction for each z-coordinate. This method assumes
 * that the source function and target function possess the same grid structure in the
 * vertical direction. If this assumption isn't satisfied then an error is thrown.
 *
 * The XY interpolation is most useful for functions that are only piecewise
 * smooth in the vertical direction. In such instances, functions that are continuous but
 * has discontinuous derivatives along specific z-coordinates will will be still be interpolated
 * with high accuracy in the transverse directions.
 *
 * In addition, since the XY interpolation is implemented on a layer by layer basis, if the
 * source function possesses a discontinuity at the layer boundaries, this discontinuity will
 * be propagated to the interpolated function.
 *
 */

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
#include "GridFunctionNd/SCC_GridFunction2d.h"
#include "GridFunctionNd/SCC_GridFunction3d.h"

#include "VLayeredGridFunNd/SCC_VLayeredGridFun3d.h"

#include "SCC_LegendreGridFunApprox2d.h"
#include "SCC_LegendreGridFunApprox3d.h"


#ifdef _OPENMP
#include <omp.h>
#ifdef OPENBLAS_THREADED
#include <cblas.h>
#endif
#endif

#ifndef VLAYERED_GRID_FUN_INTERP_3D_
#define VLAYERED_GRID_FUN_INTERP_3D_

namespace SCC
{
class VLayeredGridFunInterp3d
{
public:

	VLayeredGridFunInterp3d()
	{
		initialize();
	}

	VLayeredGridFunInterp3d(VLayeredGridFunInterp3d& V)
	{
		if(V.layerCount == 0) {initialize(); return;}
		initialize(V);
	}

	VLayeredGridFunInterp3d(int degreeX,int degreeY, int degreeZ, const SCC::VLayeredGridFun3d& dataFun3d)
	{
		initialize(degreeX, degreeY, degreeZ, dataFun3d);
	}

	void initialize()
	{
    layerCount  = 0;
    dataFun3dPtr = nullptr;

    dataXY.initialize();

	legendreGridFunApprox2d.clear();
	legendreGridFunApprox3d.clear();

	#ifdef _OPENMP
	legendreGridFunApprox2dMT.clear();
	dataXYMT.clear();
    dataFun3dPtrMT.clear();
	#endif
	}

	void initialize(VLayeredGridFunInterp3d& T)
	{
		layerCount   = T.layerCount;
		dataFun3dPtr = T.dataFun3dPtr;

		dataXY.initialize(T.dataXY);

		legendreGridFunApprox3d.clear();
		legendreGridFunApprox3d.resize(layerCount);

		legendreGridFunApprox2d.clear();
		legendreGridFunApprox2d.resize(layerCount);

	    for(long k = 0; k < layerCount; k++)
		{
		legendreGridFunApprox2d[k].initialize(T.legendreGridFunApprox2d[k]);
		legendreGridFunApprox3d[k].initialize(T.legendreGridFunApprox3d[k]);

		// Set pointer to point to *this instance of the data

		legendreGridFunApprox3d[k].FDataPtr = dataFun3dPtr->layer[k].getDataPointer();
		}

	    #ifdef _OPENMP
	    int threadCount =  omp_get_max_threads();
	    legendreGridFunApprox2dMT.resize(threadCount);
	    dataXYMT.resize(threadCount);
	    dataFun3dPtrMT.resize(threadCount);

	    for(long k = 0; k < threadCount; k++)
	    {
	    	dataXYMT[k].initialize(T.dataXYMT[k]);
	    	legendreGridFunApprox2dMT[k].initialize(T.legendreGridFunApprox2dMT[k]);
	    	dataFun3dPtrMT[k]  = T.dataFun3dPtrMT[k];
	    }
		#endif
	}

	void initialize(int degreeX,int degreeY, int degreeZ, const SCC::VLayeredGridFun3d& dataFun3d)
	{
		 if(dataFun3d.getLayerCount() == 0)
		 {
			throw std::runtime_error("\n VLayeredGridFunInterp3d: initialize(...) VLayeredGridFun3d argument not initialized.\n");
		 }

		 this->dataFun3dPtr = &dataFun3d;
		 this->layerCount   = dataFun3d.getLayerCount();

		 // For x-y-z interpolation

		 legendreGridFunApprox3d.clear();
		 legendreGridFunApprox3d.resize(layerCount);

		 // For x-y (transverse) interpolation alone assuming coincident z structure

		 legendreGridFunApprox2d.clear();
		 legendreGridFunApprox2d.resize(layerCount);

		 long xPanels; double xMin; double xMax;
		 long yPanels; double yMin; double yMax;
		 long zPanels; double zMin; double zMax;

		 std::vector<long>   zPanelArray;
		 std::vector<double> zBdrys;

		 zPanelArray = dataFun3d.getZpanels();
		 zBdrys      = dataFun3d.getZbdrys();

         xPanels = dataFun3dPtr->layer[0].getXpanelCount();
         xMin    = dataFun3dPtr->layer[0].getXmin();
         xMax    = dataFun3dPtr->layer[0].getXmax();


         yPanels = dataFun3dPtr->layer[0].getYpanelCount();
         yMin    = dataFun3dPtr->layer[0].getYmin();
         yMax    = dataFun3dPtr->layer[0].getYmax();

         dataXY.initialize(xPanels,xMin,xMax,yPanels,yMin,yMax);

		 for(long k = 0; k < layerCount; k++)
		 {
			 zPanels = zPanelArray[k];
			 zMin    = zBdrys[k];
			 zMax    = zBdrys[k+1];

            legendreGridFunApprox3d[k].initialize(degreeX, xPanels, xMin, xMax,
        		                                  degreeY, yPanels, yMin, yMax,
												  degreeZ, zPanels, zMin, zMax,
												  dataFun3dPtr->layer[k].getDataPointer());

            legendreGridFunApprox2d[k].initialize(degreeX, xPanels, xMin, xMax,
        		                                  degreeY, yPanels, yMin, yMax,
												  nullptr);
		 }

		#ifdef _OPENMP
	    int threadCount =  omp_get_max_threads();
	    legendreGridFunApprox2dMT.clear();
	    dataXYMT.clear();
        dataFun3dPtrMT.clear();

	    legendreGridFunApprox2dMT.resize(threadCount);
	    dataXYMT.resize(threadCount);
	    dataFun3dPtrMT.resize(threadCount);

	    for(long k = 0; k < threadCount; k++)
	    {
	    	dataXYMT[k].initialize(xPanels,xMin,xMax,yPanels,yMin,yMax);
	    	legendreGridFunApprox2dMT[k].initialize(degreeX, xPanels, xMin, xMax,
        		                                    degreeY, yPanels, yMin, yMax,
												    nullptr);
	    	dataFun3dPtrMT[k]  = nullptr;
	    }
		#endif

	}

    //
	// Three dimensional interpolation over each layer. The interpolation
	// formula never uses values on both sides of a layer boundary and
	// hence can interpolate to high order source (data) functions that
	// are piecewise smooth over each layer.
	//
	// To accurately interpolate a function that is discontinuous on
	// layer boundaries, the layer boundaries of the target function
	// much be coincident with the layer boundaries of the source (data)
	// function.
	//
	// The transverse domain size of the target function need not be coincident with the
	// source (data) function
	//
	void interpXYZ(SCC::VLayeredGridFun3d& outFun)
	{
		 if(outFun.getLayerCount() == 0)
		 {
			throw std::runtime_error("\n VLayeredGridFunInterp3d: interpXYZ(...) return argument not initialized.\n");
		 }

		 long xPanels; double xMin; double xMax;
		 long yPanels; double yMin; double yMax;
		 long zPanels; double zMin; double zMax;

		 std::vector<long>   zPanelArray;
		 std::vector<double> zBdrys;

		 xPanels = outFun.layer[0].getXpanelCount();
         xMin    = outFun.layer[0].getXmin();
         xMax    = outFun.layer[0].getXmax();

         yPanels = outFun.layer[0].getYpanelCount();
         yMin    = outFun.layer[0].getYmin();
         yMax    = outFun.layer[0].getYmax();

		 zPanelArray = outFun.getZpanels();
		 zBdrys      = outFun.getZbdrys();

		 double hx; double hy; double hz;

		 hx = (xMax-xMin)/(double)xPanels;
		 hy = (yMax-yMin)/(double)yPanels;

		 double xPos; double yPos; double zPos;

		 long outFunLayerCount = outFun.getLayerCount();
		 long dataLayerIndex;

		 std::vector<double> dataZbdrys = dataFun3dPtr->getZbdrys();

		 bool layerCoincidenceFlag = false;

		 if(outFunLayerCount == layerCount)
		 {
			 layerCoincidenceFlag = true;
			 for(long k = 0; k < outFunLayerCount; k++)
		     {
				 if(std::abs(dataZbdrys[k] - zBdrys[k]) > 1.0e-12*std::abs(dataZbdrys[k]))
				 {
					 layerCoincidenceFlag = false;
					 break;
				 }
		     }
		 }

		 for(long k = 0; k < outFunLayerCount; k++)
		 {
			 zPanels = zPanelArray[k];
			 zMin    = zBdrys[k];
			 zMax    = zBdrys[k+1];
		     hz      = (zMax-zMin)/(double)zPanels;

		     for(long p = 0; p <= xPanels; p++)
		     {
		     xPos = xMin + p*hx;

		     for(long q = 0; q <= yPanels;  q++)
		     {
		     yPos = yMin + q*hy;

		     for(long r = 0; r  <= zPanels; r++)
		     {
		     zPos = zMin + r*hz;

		     if(not layerCoincidenceFlag)
		     {
		     // Determine layer of input data to use if the layers are not coincident

		     if(zPos < dataZbdrys[0])           { dataLayerIndex = 0;}
		     if(zPos > dataZbdrys[layerCount])  { dataLayerIndex = layerCount-1;}

		     for(long i = 0; i < layerCount; i++)
		     {

		    	 if((zPos >= dataZbdrys[i])&&(zPos <= zBdrys[i+1]))
		    	 {
		    	 dataLayerIndex = i;
		    	 break;
		    	 }
		     }}
		     else
		     {
		    	 dataLayerIndex = k;
		     }
		     outFun.layer[k](p,q,r) = legendreGridFunApprox3d[dataLayerIndex].evaluate(xPos, yPos, zPos);
		     }}}
		 }
	}

	//
	// 3d interpolation of data values to a target function that possesses the same vertical structure
	// as the source (data) function. This interpolation just interpolates the 2D slices for each
	// discrete z-coordinate value.
	//
	// This form of interpolation will accurately interpolate functions that are just piecewise smooth
	// over each layer; in particular functions that are discontinuous across layer boundaries.
	//
	// The transverse domain size of the target function need not be coincident with the
	// source (data) function
	//
	void interpXY(SCC::VLayeredGridFun3d& outFun)
	{
	    if(outFun.getLayerCount() == 0)
		{
			throw std::runtime_error("\n VgridFunInterp3d : interpXY(...) return argument not initialized.\n");
		}
	    // Check for consistency of vertical structure

        SCC::VLayeredGridFun1d dataFunZ = dataFun3dPtr->getConstantXYslice(0, 0);
        SCC::VLayeredGridFun1d  outFunZ = outFun.getConstantXYslice(0, 0);

        if(not outFunZ.isEqualStructure(dataFunZ))
        {
        	throw std::runtime_error(std::string("\n VLayeredGridFunInterp3d: interpXY(...) vertical structure \n")
        			                             + " of VLayeredGridFun3d argument not consistent with data\n");
        }


        SCC::GridFunction2d valuesXY;

        long xPanels; double xMin; double xMax;
        long yPanels; double yMin; double yMax;
        long zPanels;

        std::vector<long>   zPanelArray;

		xPanels = outFun.layer[0].getXpanelCount();
        xMin    = outFun.layer[0].getXmin();
        xMax    = outFun.layer[0].getXmax();

        yPanels = outFun.layer[0].getYpanelCount();
        yMin    = outFun.layer[0].getYmin();
        yMax    = outFun.layer[0].getYmax();

		 zPanelArray = outFun.getZpanels();

		 double hx; double hy;

		 hx = (xMax-xMin)/(double)xPanels;
		 hy = (yMax-yMin)/(double)yPanels;

		 double xPos; double yPos;

		 #ifndef _OPENMP

		 valuesXY.initialize(xPanels,xMin,xMax,yPanels,yMin,yMax);
		 for(long k = 0; k < layerCount; k++)
		 {
			 zPanels      = zPanelArray[k];
			 for(long r = 0; r <= zPanels; r++)
			 {
             dataFun3dPtr->getConstantZslice(k,r,dataXY);

		     legendreGridFunApprox2d[k].FDataPtr = dataXY.getDataPointer();

		     for(long p = 0; p <= xPanels; p++)
		     {
		     xPos = xMin + p*hx;
		     for(long q = 0; q <= yPanels;  q++)
		     {
		     yPos = yMin + q*hy;
		     valuesXY(p,q) = legendreGridFunApprox2d[k].evaluate(xPos, yPos);
		     }}

			 outFun.setConstantZslice(k,r,valuesXY);
			 }
		 }

       #else

	   #ifdef OPENBLAS_THREADED
	   int openblas_max_threads = openblas_get_num_threads();
	   openblas_set_num_threads(1);
	   #endif

       int tI;

	   for(long k = 0; k < layerCount; k++)
	   {
			 zPanels      = zPanelArray[k];
			 #pragma omp parallel for private(tI,xPos,yPos) schedule(static,1)
			 for(long r = 0; r <= zPanels; r++)
			 {
		     tI = omp_get_thread_num();

		     for(long p = 0; p <= dataFun3dPtr->layer[k].xPanels; p++)
             {
             for(long q = 0; q <= dataFun3dPtr->layer[k].yPanels; q++)
             {
            	 dataXYMT[tI](p,q) = dataFun3dPtr->layer[k](p,q,r);
             }}

		     legendreGridFunApprox2dMT[tI].FDataPtr = dataXYMT[tI].getDataPointer();

		     for(long p = 0; p <= xPanels; p++)
		     {
		     xPos = xMin + p*hx;
		     for(long q = 0; q <= yPanels;  q++)
		     {
		     yPos = yMin + q*hy;
		     outFun.layer[k](p,q,r) = legendreGridFunApprox2dMT[tI].evaluate(xPos, yPos);
		     }}
			 }
		 }


       #endif
	}

   //
	// 3d interpolation of data values to a target function that possesses the same vertical structure
	// as the source (data) function layers from indexBegin to indexEnd.
	//
	// This interpolation just interpolates the 2D slices for each discrete z-coordinate value.
	//
	// This form of interpolation will accurately interpolate functions that are just piecewise smooth
	// over each layer; in particular functions that are discontinuous across layer boundaries.
	//
	// The transverse domain size of the target function need not be coincident with the
	// source (data) function
	//
	void interpXY(long indexBegin, long indexEnd, SCC::VLayeredGridFun3d& outFun)
	{
	    if(outFun.getLayerCount() == 0)
		{
			throw std::runtime_error("\n VgridFunInterp3d : interpXY(...) return argument not initialized.\n");
		}

	    if(outFun.getLayerCount() != (indexEnd-indexBegin) + 1)
	    {
	    	throw std::runtime_error(std::string("\n VLayeredGridFunInterp3d: interpXY(...) vertical structure \n")
        			                             + " of VLayeredGridFun3d argument not consistent with data\n");
	    }
	    // Check for consistency of vertical structure

        SCC::VLayeredGridFun1d dataFunZ = (dataFun3dPtr->getConstantXYslice(0, 0)).extractLayers(indexBegin, indexEnd);
        SCC::VLayeredGridFun1d  outFunZ = outFun.getConstantXYslice(0, 0);

        if(not outFunZ.isEqualStructure(dataFunZ))
        {
        	throw std::runtime_error(std::string("\n VLayeredGridFunInterp3d: interpXY(...) vertical structure \n")
        			                             + " of VLayeredGridFun3d argument not consistent with data\n");
        }


        SCC::GridFunction2d valuesXY;

        long xPanels; double xMin; double xMax;
        long yPanels; double yMin; double yMax;
        long zPanels;

        std::vector<long>   zPanelArray;

		xPanels = outFun.layer[0].getXpanelCount();
        xMin    = outFun.layer[0].getXmin();
        xMax    = outFun.layer[0].getXmax();

        yPanels = outFun.layer[0].getYpanelCount();
        yMin    = outFun.layer[0].getYmin();
        yMax    = outFun.layer[0].getYmax();

		 zPanelArray = outFun.getZpanels();

		 double hx; double hy;

		 hx = (xMax-xMin)/(double)xPanels;
		 hy = (yMax-yMin)/(double)yPanels;

		 double xPos; double yPos;

		 valuesXY.initialize(xPanels,xMin,xMax,yPanels,yMin,yMax);

		#ifndef _OPENMP

		 for(long k = 0; k < outFun.getLayerCount(); k++)
		 {
			 zPanels      = zPanelArray[k];
			 for(long r = 0; r <= zPanels; r++)
			 {

		     dataFun3dPtr->getConstantZslice(k + indexBegin,r,dataXY);

		     legendreGridFunApprox2d[k + indexBegin].FDataPtr = dataXY.getDataPointer();

		     for(long p = 0; p <= xPanels; p++)
		     {
		     xPos = xMin + p*hx;
		     for(long q = 0; q <= yPanels;  q++)
		     {
		     yPos = yMin + q*hy;
		     valuesXY(p,q) = legendreGridFunApprox2d[k+indexBegin].evaluate(xPos, yPos);
		     }}

			 outFun.setConstantZslice(k,r,valuesXY);
			 }
		 }
		#else
	   	#ifdef OPENBLAS_THREADED
		int openblas_max_threads = openblas_get_num_threads();
		openblas_set_num_threads(1);
	   	#endif

		int tI;

		for(long k = 0; k < outFun.getLayerCount(); k++)
		{
			 zPanels      = zPanelArray[k];
			 #pragma omp parallel for private(tI,xPos,yPos) schedule(static,1)
			 for(long r = 0; r <= zPanels; r++)
			 {
		     tI = omp_get_thread_num();

		     for(long p = 0; p <= dataFun3dPtr->layer[k + indexBegin].xPanels; p++)
             {
             for(long q = 0; q <= dataFun3dPtr->layer[k + indexBegin].yPanels; q++)
             {
            	 dataXYMT[tI](p,q) = dataFun3dPtr->layer[k + indexBegin](p,q,r);
             }}


		     legendreGridFunApprox2dMT[tI].FDataPtr = dataXYMT[tI].getDataPointer();

		     for(long p = 0; p <= xPanels; p++)
		     {
		     xPos = xMin + p*hx;
		     for(long q = 0; q <= yPanels;  q++)
		     {
		     yPos = yMin + q*hy;
		     outFun.layer[k](p,q,r) = legendreGridFunApprox2dMT[tI].evaluate(xPos, yPos);
		     }}
			 }
		 }


		#endif



	}

    long layerCount;

    const SCC::VLayeredGridFun3d*    dataFun3dPtr;
    SCC::GridFunction2d              dataXY;

	std::vector<LegendreGridFunApprox2d> legendreGridFunApprox2d;
	std::vector<LegendreGridFunApprox3d> legendreGridFunApprox3d;


	#ifdef _OPENMP
	std::vector<LegendreGridFunApprox2d> legendreGridFunApprox2dMT;
	std::vector<SCC::GridFunction2d>                      dataXYMT;
	std::vector<const SCC::VLayeredGridFun3d*>      dataFun3dPtrMT;
	#endif
};

} // SCC namespace


#endif /* VLAYERED_GRID_FUN_INTERP_3D_ */
