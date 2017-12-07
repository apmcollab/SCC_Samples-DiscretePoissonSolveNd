#include <iostream>
#include <cmath>
#include <functional>
using namespace std;
//
// fftw_DiscretePoissonSolve1d.cpp
//
// A test code that tests the use of FFTW to create a solution the the discrete
// 1 dimensional Laplace equation with homogeneous Dirichlet boundary conditions.
//
// The discrete equations to be solved is
//
// alpha*([D+D-]_x) u  = f
//
// where [D+D-]_x is the standard 3 point discrete Laplace operator with homogeneous
// Dirichlet boundary conditions.
//
// Note: This is an inefficient way to solve this particular discrete problem, one would
// typically use  a compact tri-diagonal solver.
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
//       |                |
//   [TestPrograms]       SCC
//                        |
//                   -------------------
//                 SCC Component directories
//
//
// where the directory SCC contains the SCC components MakeScripts, DoubleVectorNd, GridFunctionNd and
// FFTW3_InterfaceNd source directories.
//
// The command line compilation command is
//
// g++ fftw_DiscretePoissonSolve1d.cpp -O2  -std=c++11 -I./SCC -lfftw3 -o fftw_DiscretePoissonSolve1d.exe
//
// Alternately, if one is using MakeScripts, the build command executed is
//
// make -f fftw_DiscretePoissonSolve1d.mk release
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
#include "GridFunctionNd/SCC_GridFunction1d.h"
#include "DoubleVectorNd/SCC_DoubleVector1d.h"
#include "FFTW3_InterfaceNd/SCC_fftw3_sin1d.h"

//
// This class implements the identity operator at the boundary
// points of the SCC::GridFunction1d and an centered difference approximation
// at the interior points.
//
//
class LaplaceOp1d
{
public:
	LaplaceOp1d()
	{this->alpha = 0.0;}

	LaplaceOp1d(double alpha)
	{this->alpha = alpha;}

	LaplaceOp1d(const LaplaceOp1d& Op)
	{this->alpha = Op.alpha;}

	void apply(const SCC::GridFunction1d& F,SCC::GridFunction1d& Fout)
	{
		Fout.initialize(F); // Copies values, so Fout has same boundary
                           // values as F

		double hx = F.getHx();
		long   M  = F.getXpanelCount();


		// Second difference approximation at interior
		// points

		for(long i = 1; i < M; i++)
		{
			Fout(i) = (alpha/(hx*hx))*(F(i-1) - 2.0*F(i) + F(i+1));
		}
	}


	double alpha;

};



int main()
{
	long xPanels  = 20;
    long exponent = 5;  // This is the exponent for the exact solution and must be >= 2

//   Problem set up.
//
//   Use non-unit domain so that this program demonstrates
//   how to scale the transforms of derivatives appropriately
 
	double xMin = -3.0;
	double xMax =  2.0;

	double LX   = (xMax-xMin);
	double hx   = LX/xPanels;

	double pi   =  3.141592653589793238;

	double alpha  = -1.0;           // Coefficient of Laplace operator


	LaplaceOp1d laplaceOp1d(alpha); // Discrete Laplace operator
	
	// Create

    SCC::GridFunction1d       f(xPanels,xMin,xMax);
    SCC::GridFunction1d       u(xPanels,xMin,xMax);
    SCC::GridFunction1d  uExact(xPanels,xMin,xMax);
    SCC::GridFunction1d  uError(xPanels,xMin,xMax);
    SCC::GridFunction1d  fExact(xPanels,xMin,xMax);
    SCC::GridFunction1d  fError(xPanels,xMin,xMax);

    f.setToValue(0.0);
    u.setToValue(0.0);

    //
	// uFunction is the function (1-x^2)^exponent scaled and centered so that
    // the support of u(x) spans the width of the computational domain, e.g.
    // if vanishes on the  boundary of the domain.
    //

    std::function<double(double)> uFunction = [xMin,xMax,exponent](double x)
	{
    double xBar  = (xMin+xMax)/2.0;
    double xStar = (x-xBar)/((xMax-xMin)/2.0);
	if(x <= xMin) return 0.0;
	if(x >= xMax) return 0.0;
    return pow(1.0 - xStar*xStar,exponent);
	};


    // Create a right hand side by applying the discrete Laplace operator to
    // the specified solution.

    uExact.specify(uFunction);              // Evaluate uFunction at nodes of uExact

    laplaceOp1d.apply(uExact,fExact);       // Apply discrete Laplace operator to determine the right hand side

    // Initialize the FFTW3 interface routine

    SCC::fftw3_sin1d DFT;       // Discrete sin transform interface
    DFT.initialize(xPanels);   // Specify the number of grid panels!

    // Checking FFTW evaluation of the alpha*[D+D-]_x u

    SCC::DoubleVector1d uTransform(xPanels-1);
    uTransform.setToValue(0.0);

    DFT.fftw1d_sin_forward(uExact,uTransform); // Note: We can use a GridFunction1d as an argument since it's extended
                                               // from a DoubleVector1d

    double lambda;

    //
    // Loop over wave numbers. Note indexing of uTransform shifted by 1
    //
    for(long kx = 1; kx <= xPanels-1; kx++)
    {
    	lambda            = alpha*((2.0/(hx*hx))*(cos((kx*pi*hx)/LX) - 1.0));
    	uTransform(kx-1) *= lambda;
    }

    DFT.fftw1d_sin_inverse(uTransform,f);      // Note: We can use a GridFunction1d as an argument since it's extended
                                               // from a DoubleVector1d

    fError = f-  fExact;

    printf("Error in FFT evaluation of (alpha*[D+D-]_x) u : %10.5e \n",fError.norm2());


    // Solving the equation alpha*[D+D-]_x  = f using fftw

    SCC::DoubleVector1d fTransform(xPanels-1);
    fTransform.setToValue(0.0);

    DFT.fftw1d_sin_forward(fExact,fTransform); // Note: We can use a GridFunction1d as an argument since it's extended
                                               // from a DoubleVector1d

    for(long kx = 1; kx <= xPanels-1; kx++)
    {
    	lambda            = alpha*((2.0/(hx*hx))*(cos((kx*pi*hx)/LX) - 1.0));
    	fTransform(kx-1) /= lambda;
    }

    DFT.fftw1d_sin_inverse(fTransform,u);      // Note: We can use a GridFunction1d as an argument since it's extended
                                               // from a DoubleVector1d

    uError = u - uExact;

    printf("Error in FFT evaluation of the solution to  alpha*([D+D-]_x) u = f : %10.5e \n",uError.norm2());

}

