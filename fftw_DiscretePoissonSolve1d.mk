#
# fftw_DiscretePoissonSolve1d program makefile using MakeScripts
#
# This makefile assumes a directory structure
#
# [Test directory]
#---------------------------
#         |                |                
#   [TestProgram]         SCC
#                          |
#                   -------------------
#                   SCC Component directories
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


CPPfiles   +=  fftw_DiscretePoissonSolve1d.cpp
INCLUDES    = -I./
INCLUDES   += -I$(SCC_Dir)

# Specifying FFTW3 libraries and defines when OpenMP version is to be used. 

LIBS     += $(FFTW_LIB)  
LIB_PATH += $(FFTW_PATH)

ifeq ($(_FFTW_OPENMP),1)
CXXDEFINES   += -D_FFTW_OPENMP
endif


RELEASE_EXEC  = fftw_DiscretePoissonSolve1d.exe
DEBUG_EXEC    = fftw_DiscretePoissonSolve1d_debug.exe

RELEASE_DIR  = ./_releasefftw_DiscretePoissonSolve1d
DEBUG_DIR    = ./_debugfftw_DiscretePoissonSolve1d

include $(SCC_Dir)/MakeScripts/ExecutableMake.mk

