#include <iostream>
#include <cmath>
#include <functional>
using namespace std;
//
// fftw_DiscretePoissonSolve2d.cpp
//
// A test code that tests the use of FFTW to create a solution the the discrete
// 2 dimensional Laplace equation with homogeneous Dirichlet boundary conditions.
//
// The discrete equations to be solved is
//
// alpha*([D+D-]_x + [D+D-]_y) u  = f
//
// where ([D+D-]_x + [D+D-]_y) is the standard 5 point discrete Laplace operator with homogeneous
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
// The command line compilation command is
//
// g++ fftw_DiscretePoissonSolve2d.cpp -O2  -std=c++11 -I./SCC -lfftw3 -o fftw_DiscretePoissonSolve2d.exe
//
// Alternately, if one is using MakeScripts, the build command executed is
//
// make -f fftw_DiscretePoissonSolve2d.mk release
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
#include "GridFunctionNd/SCC_GridFunction2d.h"
#include "DoubleVectorNd/SCC_DoubleVector2d.h"
#include "FFTW3_InterfaceNd/SCC_fftw3_sin2d.h"
#include "FFTW3_InterfaceNd/SCC_FFT_Nvalues.h"

#include "Timing/ClockIt.h"

//
// This class implements the identity operator at the boundary
// points of the SCC::GridFunction2d and an centered difference approximation
// at the interior points.
//
//
class LaplaceOp2d
{
public:
	LaplaceOp2d()
	{this->alpha = 0.0;}

	LaplaceOp2d(double alpha)
	{this->alpha = alpha;}

	LaplaceOp2d(const LaplaceOp2d& Op)
	{this->alpha = Op.alpha;}

	void apply(const SCC::GridFunction2d& F,SCC::GridFunction2d& Fout)
	{
		Fout.initialize(F); // Copies values, so Fout has same boundary
		                    // values as F

		double hx = F.getHx();
		double hy = F.getHy();

		long   M  = F.getXpanelCount();
		long   N  = F.getYpanelCount();


		// Second difference approximation at interior points

		for(long i = 1; i < M; i++)
		{
       for(long j = 1; j < N; j++)
       {
			Fout(i,j) = (alpha/(hx*hx))*(F(i-1,j) - 2.0*F(i,j) + F(i+1,j)) +
					    (alpha/(hy*hy))*(F(i,j-1) - 2.0*F(i,j) + F(i,j+1));
		}}
	}

	double alpha;

};



int main()
{
	long xPanels  = 43;
	long yPanels  = 34;

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

//   Problem set up.
//
//   Use non-unit domain so that this program demonstrates
//   how to scale the transforms of derivatives appropriately
 
	double xMin   = -3.0;
	double xMax   =  2.0;

	double yMin   = -2.0;
	double yMax   =  2.0;

	double LX    = (xMax-xMin);
	double LY    = (yMax-yMin);

	double hx    = LX/xPanels;
	double hy    = LY/yPanels;

	double pi    =  3.141592653589793238;

	double alpha  = -1.0;  // Coefficient of laplace operator
	
	LaplaceOp2d laplaceOp2d(alpha); // Discrete Laplace operator

    SCC::GridFunction2d      f(xPanels,xMin,xMax,yPanels,yMin,yMax);
    SCC::GridFunction2d      u(xPanels,xMin,xMax,yPanels,yMin,yMax);
    SCC::GridFunction2d  uExact(xPanels,xMin,xMax,yPanels,yMin,yMax);
    SCC::GridFunction2d  uError(xPanels,xMin,xMax,yPanels,yMin,yMax);
    SCC::GridFunction2d  fExact(xPanels,xMin,xMax,yPanels,yMin,yMax);
    SCC::GridFunction2d  fError(xPanels,xMin,xMax,yPanels,yMin,yMax);

    f.setToValue(0.0);
    u.setToValue(0.0);


	// uFunction is the function (1-r^2))^exponent scaled and centered so that
    // the support of u(x,y) spans the width of the computational domain.

	double xCent  = (xMin + xMax)/2.0;
	double yCent  = (yMin + yMax)/2.0;
	double radius = ( (xMax-xMin)  < (yMax-yMin)) ? (xMax-xMin)/2.0 : (yMax-yMin)/2.0;

    std::function<double(double,double)> uFunction = [xCent,yCent,radius,exponent](double x, double y)
	{
    double r2 = ((x-xCent)*(x-xCent) + (y - yCent)*(y-yCent))/(radius*radius);
	if(r2 > 1.0) return 0.0;
    return pow(1.0 - r2,exponent);
	};

    // Create a right hand side by applying the discrete Laplace operator to
    // the specified solution.

    uExact.specify(uFunction);

    laplaceOp2d.apply(uExact,fExact);

    // Initialize the FFTW3 interface routine

    SCC::fftw3_sin2d DFT;            // Discrete sin transform interface
    DFT.initialize(xPanels,yPanels);   // Specify the number of grid panels!

    // Checking spectral evaluation of the alpha*([D+D-]_x + [D+D-]_y)

    SCC::DoubleVector2d uTransform(xPanels-1,yPanels-1);
    uTransform.setToValue(0.0);

    DFT.fftw2d_sin_forward(uExact,uTransform); // Note: We can use a GridFunction2d as an argument since it's extended
                                      // from a DoubleVector2d

    double lambda;

    //
    // Loop over wave numbers. Note indexing of uTransform shifted by 1
    //

    for(long kx = 1; kx <= xPanels-1; kx++)
    {
    for(long ky = 1; ky <= yPanels-1; ky++)
    {
    	lambda = alpha*( (2.0/(hx*hx))*(cos((kx*pi*hx)/LX) - 1.0)
    			      + (2.0/(hy*hy))*(cos((ky*pi*hy)/LY) - 1.0));
    	uTransform(kx-1,ky-1) *= lambda;
    }}

    DFT.fftw2d_sin_inverse(uTransform,f);     // Note: We can use a GridFunction2d as an argument since it's extended
                                      // from a DoubleVector2d

    fError = f-  fExact;

    printf("\n");
    printf("Error in FFT evaluation of alpha*([D+D-]_x + [D+D-]_y) u   : %10.5e \n",fError.norm2());


    // Solving the equation alpha*(u_xx + u_yy) = f using a spectral discretization

    SCC::DoubleVector2d fTransform(xPanels-1,yPanels-1);
    fTransform.setToValue(0.0);

    DFT.fftw2d_sin_forward(fExact,fTransform); // Note: We can use a GridFunction2d as an argument since it's extended
                                               // from a DoubleVector2d


    //
    // Loop over wave numbers. Note indexing of uTransform shifted by 1
    //


    for(long kx = 1; kx <= xPanels-1; kx++)
    {
    for(long ky = 1; ky <= yPanels-1; ky++)
    {
    	lambda = alpha*( (2.0/(hx*hx))*(cos((kx*pi*hx)/LX) - 1.0)
    	    			      + (2.0/(hy*hy))*(cos((ky*pi*hy)/LY) - 1.0));
    	fTransform(kx-1,ky-1) /= lambda;
    }}


    DFT.fftw2d_sin_inverse(fTransform,u);     // Note: We can use a GridFunction2d as an argument since it's extended
                                              // from a DoubleVector2d

    uError = u - uExact;

    printf("Error in FFT valuation of the solution to  alpha*([D+D-]_x + [D+D-]_y) u = f : %10.5e \n",uError.norm2());


   //
    // Timing
    //

    long panelIncrement = 50;
    long incrementCount = 4;
    long repetitions    = 40;
    double timeTaken    = 0.0;

    ClockIt clockTimer;


    printf("\n\nTiming Runs (%ld repetitions) \n\n",repetitions);

    for(long k = 0; k < incrementCount; k++)
    {


    xPanels  = panelIncrement*(k+1);
	yPanels  = panelIncrement*(k+1);

    f.initialize(xPanels,xMin,xMax,yPanels,yMin,yMax);
    u.initialize(xPanels,xMin,xMax,yPanels,yMin,yMax);

    u.specify(uFunction);
    laplaceOp2d.apply(u,f);

    // Solve system

	DFT.initialize(xPanels,yPanels);

    fTransform.initialize(xPanels-1,yPanels-1);
    fTransform.setToValue(0.0);

    clockTimer.start();

    for(long kTimes = 0; kTimes < repetitions; kTimes++)
    {

    DFT.fftw2d_sin_forward(f,fTransform);

    for(long kx = 1; kx <= xPanels-1; kx++)
    {
    for(long ky = 1; ky <= yPanels-1; ky++)
    {
    	lambda = alpha*( (2.0/(hx*hx))*(cos((kx*pi*hx)/LX) - 1.0)
    	    		   + (2.0/(hy*hy))*(cos((ky*pi*hy)/LY) - 1.0));

    	fTransform(kx-1,ky-1) /= lambda;
    }}

    DFT.fftw2d_sin_inverse(fTransform,u);      // Note: We can use a GridFunction3d as an argument since it's extended
                                               // from a DoubleVector3d

    clockTimer.stop();
    }

    timeTaken = clockTimer.getSecElapsedTime()/(double)repetitions;
    printf("FFT Solution time (sec) %3ld X %3ld  :  %10.5e \n",xPanels,yPanels,timeTaken);

    }
}

