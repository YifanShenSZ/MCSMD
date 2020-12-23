###########################################
#                                         #
#        Makefile for SurfGenBound        #
#                                         #
###########################################

compiler = ifort
MyLibDir = /home-4/yshen57@jhu.edu/Library/Fortran-Library
MyLib = $(MyLibDir)/General.f90 $(MyLibDir)/Mathematics.f90 $(MyLibDir)/LinearAlgebra.f90 \
$(MyLibDir)/mkl_rci.f90 $(MyLibDir)/NonlinearOptimization.f90
src = Basic.f90 Wigner.f90 SMD.f90 Main.f90
exe = MCSMD.exe
# -static causes strange bug in ifort 2018
flag = -m64 -xCORE-AVX2 -mtune=core-avx2 -mkl -ipo -O3 -no-prec-div -fp-model fast=2

$(exe): $(MyLib) $(src)
	$(compiler) $(flag) $^ -o $(exe)

clean:
	rm *.mod