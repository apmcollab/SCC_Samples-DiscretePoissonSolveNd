#include <iostream>
#include <cmath>
#include <functional>
using namespace std;
//
// fftw_DiscretePoissonSolve3d.cpp
//
// A test code that tests the use of FFTW to create a solution the the discrete
// 3 dimensional Laplace equation with homogeneous Dirichlet boundary conditions.
//
// The discrete equations to be solved is
//
// alpha*([D+D-]_x + [D+D-]_y + [D+D-]_z) u = f
//
// where ([D+D-]_x + [D+D-]_y + [D+D-]_z) is the standard 7 point discrete Laplace operator with homogeneous
// Dirichlet boundary conditions.
//
// This program depends on an installation of the FFTW3 libraries. You may need to modify
// the include paths, library paths, and library to link to an specific installation of
// FFTW3. This program was created an tested on an Ubuntu linux system in which the FFTW3
// libraries were installed using the synaptic package manager.
//
// The test program assumes the directory structure
//
//
// [Test directory]
//---------------------------
//      |             |
//   [TestPrograms]      SCC
//                    |
//                -------------------
//              SCC Component directories
//
//
// where the directory SCC contains the SCC components MakeScripts, DoubleVectorNd, GridFunctionNd and
// FFTW3_InterfaceNd source directories.
//
// The command line compilation command for single threaded execution is
//
// g++ fftw_DiscretePoissonSolve3d.cpp -std=c++11 -I./SCC -lfftw3 -o fftw_DiscretePoissonSolve3d.exe
//
// for multi-threaded execution the command is
//
// g++ fftw_DiscretePoissonSolve3d.cpp -D_FFTW_OPENMP -fopenmp -std=c++11 -I./SCC -lfftw3_omp -lfftw3 -o fftw_DiscretePoissonSolve3d.exe
//
// Alternately, if one is using MakeScripts the build command is
//
// make -f fftw_DiscretePoissonSolve3d.mk release
//
// or for multi-threaded execution
//
// make -f fftw_DiscretePoissonSolve3d.mk release FFTW_OPENMP=1
//
// Nov. 27, 2017
//
//
/*
#############################################################################
#
# Copyright 2015-17 Chris Anderson
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
#include "GridFunctionNd/SCC_GridFunction3d.h"
#include "DoubleVectorNd/SCC_DoubleVector3d.h"
#include "FFTW3_InterfaceNd/SCC_fftw3_sin3d.h"
#include "FFTW3_InterfaceNd/SCC_FFT_Nvalues.h"

#include "Timing/ClockIt.h"

#ifdef _FFTW_OPENMP
#include <omp.h>
#endif

//
// This class implements the identity operator at the boundary
// points of the SCC::GridFunction3d and an centered difference approximation
// at the interior points.
//
//
class LaplaceOp3d
{
public:
	LaplaceOp3d()
	{this->alpha = 0.0;}

	LaplaceOp3d(double alpha)
	{this->alpha = alpha;}

	LaplaceOp3d(const LaplaceOp3d& Op)
	{this->alpha = Op.alpha;}

	void apply(const SCC::GridFunction3d& F,SCC::GridFunction3d& Fout)
	{
		Fout.initialize(F);      // Copies values, so Fout has same boundary
		                         // values as F

		double hx = F.getHx();
		double hy = F.getHy();
		double hz = F.getHz();

		long   M  = F.getXpanelCount();
		long   N  = F.getYpanelCount();
		long   P  = F.getZpanelCount();


		// Second difference approximation at interior points

		for(long i = 1; i < M; i++)
		{
        for(long j = 1; j < N; j++)
        {
        for(long k = 1; k < P; k++)
        {
			Fout(i,j,k) = (alpha/(hx*hx))*(F(i-1,j,k) - 2.0*F(i,j,k) + F(i+1,j,k)) +
					    (alpha/(hy*hy))*(F(i,j-1,k) - 2.0*F(i,j,k) + F(i,j+1,k)) +
						(alpha/(hz*hz))*(F(i,j,k-1) - 2.0*F(i,j,k) + F(i,j,k+1));
		}}}
	}

	double alpha;

};


int main()
{
	string threadCountInput = "-1";

	#ifdef _FFTW_OPENMP
    int threadIndex;
    int threadCount;
    if(not threadCountInput.empty())
    {
    threadCount = atoi(threadCountInput.c_str());
    }
    if(threadCount > omp_get_max_threads()){threadCount = omp_get_max_threads();}
    if(threadCount <= 0)                   {threadCount = omp_get_max_threads();}
    omp_set_num_threads(threadCount);

	printf("\n");
    printf("#############\n");
	printf("############# Using OpenMP With %d Threads\n",omp_get_max_threads());
	printf("#############\n");
	printf("\n");
    #endif

	long xPanels  = 43;
	long yPanels  = 34;
    long zPanels  = 23;

    long exponent = 5;  // This is the exponent for the exact solution.
                        // The exponent must be >= 2

    // Reset panel count so efficient FFT's are used; the new
    // panel count is the next larger value that's a product
    // of primes < 13

    SCC::FFT_Nvalues fft_Nvalues;

    printf("Original xPanels = %ld ",xPanels);
    xPanels = fft_Nvalues.getFFT_N(xPanels);
    printf("::: New xPanels = %ld \n",xPanels);

    printf("Original yPanels = %ld ",yPanels);
    yPanels = fft_Nvalues.getFFT_N(yPanels);
    printf("::: New yPanels = %ld \n",yPanels);

    printf("Original zPanels = %ld ",zPanels);
    zPanels = fft_Nvalues.getFFT_N(zPanels);
    printf("::: New zPanels = %ld \n",zPanels);

//   Problem set up.
//
//   Use non-unit domain so that this program demonstrates
//   how to scale the transforms of derivatives appropriately
 
	double xMin   = -3.0;
	double xMax   =  2.0;

	double yMin   = -2.0;
	double yMax   =  2.0;

	double zMin   = -1.0;
	double zMax   =  2.0;

	double LX     = (xMax-xMin);
	double LY     = (yMax-yMin);
	double LZ     = (zMax-zMin);

	double hx    = LX/xPanels;
	double hy    = LY/yPanels;
	double hz    = LZ/zPanels;

	double pi     =  3.141592653589793238;

	double alpha  = -1.0;  // Coefficient of laplace operator
	
	LaplaceOp3d laplaceOp3d(alpha); // Discrete Laplace operator


    SCC::GridFunction3d       f(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
    SCC::GridFunction3d       u(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
    SCC::GridFunction3d  uExact(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
    SCC::GridFunction3d  uError(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
    SCC::GridFunction3d  fExact(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
    SCC::GridFunction3d  fError(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);

    f.setToValue(0.0);
    u.setToValue(0.0);

	// uFunction is the function (1-r^2))^exponent scaled and centered so that
    // the support of u(x,y,z) spans the width of the computational domain.


	double xCent  = (xMin + xMax)/2.0;
	double yCent  = (yMin + yMax)/2.0;
	double zCent  = (zMin + zMax)/2.0;
	double radius      = (xMax-xMin)/2.0;
	radius = ( radius  < (yMax-yMin)/2.0) ? radius: (yMax-yMin)/2.0;
	radius = ( radius  < (zMax-zMin)/2.0) ? radius: (zMax-zMin)/2.0;

    std::function<double(double,double,double)> uFunction = [xCent,yCent,zCent,radius,exponent](double x, double y,double z)
	{
    double r2 = ((x-xCent)*(x-xCent) + (y - yCent)*(y-yCent) + (z - zCent)*(z-zCent))/(radius*radius);
	if(r2 > 1.0) return 0.0;
    return pow(1.0 - r2,exponent);
	};

    // Create a right hand side by applying the discrete Laplace operator to
    // the specified solution.

    uExact.specify(uFunction);

    laplaceOp3d.apply(uExact,fExact);

    // Initialize the FFTW3 interface routine

    SCC::fftw3_sin3d DFT;                    // Discrete sin transform interface

    #ifdef _FFTW_OPENMP
	DFT.initialize(xPanels,yPanels,zPanels,threadCount);
	#else
	DFT.initialize(xPanels,yPanels,zPanels);
	#endif

    // Checking FFT evaluation of the alpha*([D+D-]_x + [D+D-]_y + [D+D-]_z) u

    SCC::DoubleVector3d uTransform(xPanels-1,yPanels-1,zPanels-1);
    uTransform.setToValue(0.0);

    DFT.fftw3d_sin_forward(uExact,uTransform); // Note: We can use a GridFunction3d as an argument since it's extended
                                               // from a DoubleVector3d

    double lambda;

    for(long kx = 1; kx <= xPanels-1; kx++)
    {
    for(long ky = 1; ky <= yPanels-1; ky++)
    {
    for(long kz = 1; kz <= zPanels-1; kz++)
    {

    	lambda = alpha*( (2.0/(hx*hx))*(cos((kx*pi*hx)/LX) - 1.0)
    			       + (2.0/(hy*hy))*(cos((ky*pi*hy)/LY) - 1.0)
					   + (2.0/(hz*hz))*(cos((kz*pi*hz)/LZ) - 1.0));

    	uTransform(kx-1,ky-1,kz-1) *= lambda;
    }}}

    DFT.fftw3d_sin_inverse(uTransform,f);      // Note: We can use a GridFunction3d as an argument since it's extended
                                               // from a DoubleVector3d

    fError = f-  fExact;

    printf("\n");
    printf("Error in FFT evaluation of alpha*(u_xx + u_yy  +u_zz)                     : %10.5e \n",fError.norm2());


    // Solving the equation alpha*([D+D-]_x + [D+D-]_y + [D+D-]_z) u = f  using FFT's

    SCC::DoubleVector3d fTransform(xPanels-1,yPanels-1,zPanels-1);
    fTransform.setToValue(0.0);

    DFT.fftw3d_sin_forward(fExact,fTransform); // Note: We can use a GridFunction3d as an argument since it's extended
                                               // from a DoubleVector3d

    for(long kx = 1; kx <= xPanels-1; kx++)
    {
    for(long ky = 1; ky <= yPanels-1; ky++)
    {
    for(long kz = 1; kz <= zPanels-1; kz++)
    {
    	lambda = alpha*( (2.0/(hx*hx))*(cos((kx*pi*hx)/LX) - 1.0)
    	    		   + (2.0/(hy*hy))*(cos((ky*pi*hy)/LY) - 1.0)
    		           + (2.0/(hz*hz))*(cos((kz*pi*hz)/LZ) - 1.0));

    	fTransform(kx-1,ky-1,kz-1) /= lambda;
    }}}

    DFT.fftw3d_sin_inverse(fTransform,u);      // Note: We can use a GridFunction3d as an argument since it's extended
                                               // from a DoubleVector3d

    uError = u - uExact;

    printf("Error in FFT evaluation of the solution to alpha*([D+D-]_x + [D+D-]_y + [D+D-]_z) u = f : %10.5e \n",uError.norm2());


    //
    // Timing
    //

    long panelIncrement = 50;
    long incrementCount = 3;
    long repetitions    = 5;
    double timeTaken    = 0.0;

    ClockIt clockTimer;


    printf("\n\nTiming Runs (%ld repetitions) \n\n",repetitions);

    for(long k = 0; k < incrementCount; k++)
    {


    xPanels  = panelIncrement*(k+1);
	yPanels  = panelIncrement*(k+1);
    zPanels   =panelIncrement*(k+1);

    f.initialize(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
    u.initialize(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);

    u.specify(uFunction);
    laplaceOp3d.apply(u,f);

    // Solve system

    #ifdef _FFTW_OPENMP
	DFT.initialize(xPanels,yPanels,zPanels,threadCount);
	#else
	DFT.initialize(xPanels,yPanels,zPanels);
	#endif

    fTransform.initialize(xPanels-1,yPanels-1,zPanels-1);
    fTransform.setToValue(0.0);

    clockTimer.start();

    for(long kTimes = 0; kTimes < repetitions; kTimes++)
    {

    DFT.fftw3d_sin_forward(f,fTransform);

    for(long kx = 1; kx <= xPanels-1; kx++)
    {
    for(long ky = 1; ky <= yPanels-1; ky++)
    {
    for(long kz = 1; kz <= zPanels-1; kz++)
    {
    	lambda = alpha*( (2.0/(hx*hx))*(cos((kx*pi*hx)/LX) - 1.0)
    	    		   + (2.0/(hy*hy))*(cos((ky*pi*hy)/LY) - 1.0)
    		           + (2.0/(hz*hz))*(cos((kz*pi*hz)/LZ) - 1.0));

    	fTransform(kx-1,ky-1,kz-1) /= lambda;
    }}}

    DFT.fftw3d_sin_inverse(fTransform,u);      // Note: We can use a GridFunction3d as an argument since it's extended
                                               // from a DoubleVector3d

    clockTimer.stop();
    }

    timeTaken = clockTimer.getSecElapsedTime()/(double)repetitions;
    printf("FFT Solution time (sec) %3ld X %3ld X %3ld :  %10.5e \n",xPanels,yPanels,zPanels,timeTaken);

    }
}

