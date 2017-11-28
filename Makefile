#
# This is the "meta" makefile that builds the sample executables
# by invoking their respective makescripts.
#
# Usage :
#
# make release                   === builds release versions of all executables
# make debug                     === builds debug versions of all executables
#
# make clean                     === removes object files and executables
# make cleanall                  === removes object files, executables and temporary directories
#
#
# To build FFTW 3d samples that use the multi-threaded capabilities, add the specification
# _FFTW_OPENMP=1 to the make invocation, e.g. make .FFTWsinSolveTest _FFTW_OPENMP=1 
# 
# Oct. 5, 2016
#
SHELL=/bin/sh

# Location of SCC source files and Makescripts files

SCC_Dir=./SCC
MAKESCRIPTS_Dir=$(SCC_Dir)/MakeScripts

# Use BaseCommonConfig.mk if it exists, otherwise use BaseCommonConfig_Default.mk 

ifneq ("$(wildcard $(MAKESCRIPTS_Dir)/BaseCommonConfig.mk)","")
	include $(MAKESCRIPTS_Dir)/BaseCommonConfig.mk
else
	include $(MAKESCRIPTS_Dir)/BaseCommonConfig_Default.mk
endif


ifeq ($(MAKECMDGOALS),release)
BUILD_TYPE=release
endif
ifeq ($(MAKECMDGOALS),debug)
BUILD_TYPE=debug
endif


release : .Starting   .fftw_DiscretePoissonSolve1d   .fftw_DiscretePoissonSolve2d .fftw_DiscretePoissonSolve3d .Finished 

debug   : .Starting   .fftw_DiscretePoissonSolve1d_debug .fftw_DiscretePoissonSolve2d_debug  .fftw_DiscretePoissonSolve3d_debug .Finished  

.Starting :
	###########################################################
	#            Build Started
	########################################################### 
	$(QUIET) echo "Compilation Date : " `date` 
ifeq ($(OpenMP),0)
	$(QUIET) echo "Executable Type  :  Single-threaded" 
else
	$(QUIET) echo "Executable Type  :  OpenMP based multi-threaded" 
endif

 
.fftw_DiscretePoissonSolve1d :
	###########################################################
	#               .fftw_DiscretePoissonSolve1d
	###########################################################
	$(QUIET)$(MAKE) -f fftw_DiscretePoissonSolve1d.mk  release 
.fftw_DiscretePoissonSolve1d_debug :
	###########################################################
	#               .fftw_DiscretePoissonSolve1d_debug 
	########################################################### 
	$(QUIET)$(MAKE) -f fftw_DiscretePoissonSolve1d.mk  debug

.fftw_DiscretePoissonSolve2d :
	###########################################################
	#               .fftw_DiscretePoissonSolve2d
	###########################################################
	$(QUIET)$(MAKE) -f fftw_DiscretePoissonSolve2d.mk  release 
.fftw_DiscretePoissonSolve2d_debug :
	###########################################################
	#               .fftw_DiscretePoissonSolve2d_debug 
	########################################################### 
	$(QUIET)$(MAKE) -f fftw_DiscretePoissonSolve2d.mk  debug


.fftw_DiscretePoissonSolve3d :
	###########################################################
	#               .fftw_DiscretePoissonSolve3d
	###########################################################
	$(QUIET)$(MAKE) -f fftw_DiscretePoissonSolve3d.mk  release 
.fftw_DiscretePoissonSolve3d_debug :
	###########################################################
	#               .fftw_DiscretePoissonSolve3d_debug 
	########################################################### 
	$(QUIET)$(MAKE) -f fftw_DiscretePoissonSolve3d.mk  debug
	

.Finished :
	###########################################################
	#            Build Completed
	########################################################### 
.PHONY : cleanall clean
cleanall  : clean
	rm -rf ./_*
clean: 
	$(QUIET)$(MAKE) -f fftw_DiscretePoissonSolve1d.mk       clean 
	$(QUIET)$(MAKE) -f fftw_DiscretePoissonSolve2d.mk       clean
	$(QUIET)$(MAKE) -f fftw_DiscretePoissonSolve3d.mk       clean
 



