# @
# @         @
# @             @       @  @@@@@@@  @ @@  @ @@ @@@    @@@@@   @ @@   @@@@@@@
# @         @   @       @ @       @ @@  @ @@  @   @  @     @  @@  @ @       @
# @         @    @     @  @@@@@@@@@ @     @   @   @ @       @ @     @@@@@@@@@
# @         @  @  @   @   @         @     @   @   @  @     @  @     @
# @@@@@@@@@  @@    @@@     @@@@@@@  @     @   @   @   @@@@@   @      @@@@@@@
#
# @       @                                    @@
# @      @                                      @
# @     @    @@@@@@@  @ @@  @ @@@@@    @@@@@@@  @  @@@@@@@
# @@@@@@    @       @ @@  @ @@     @  @       @ @ @
# @     @   @@@@@@@@@ @     @       @ @@@@@@@@@ @  @@@@@@@
# @      @  @         @     @       @ @         @         @
# @       @  @@@@@@@  @     @       @  @@@@@@@ @@@ @@@@@@@
#
# This makes the Livermore Kernels Benchmark in both Fortran and C
# The executables are as follows:
#       mflops   - Fortran SINGLE precision benchmark
#       dmflops  - Fortran DOUBLE precision benchmark
#       cmflops  - C (kernel) SINGLE precision benchmark
#       dcmflops - C (kernel) DOUBLE precision benchmark
# The C benchmarks use the Fortran support routines and driver found in
# the file "mflops.f".  Most variables are passed though common blocks.
# The C source for the kernels is in "cernel.c" and has structures defined
# to pick up all the variables found in the Fortran common blocks.  We have
# striven to make this interface as portable as possible, but it may need
# some work for any "new" system.
#
# What is your Fortran 77 compiler called...
FC=  xlf
# Use the highest optimization flags your compiler can stand...
FFLAGS=  -O3 -qstrict -qtune=604 -qarch=ppc -qnomaf
FBASE= 
#
# What is your C compiler called....
CC= cc
# Use the highest optimization flags your compiler can stand...
CFLAGS= -O
#
# If you need any special load flags put them here...
LDFLAGS=-lm
#
# Default make rules for Fortran and C
.f.o :
	$(FC) $(FFLAGS) -c $*.f
.c.o :
	$(CC) $(CFLAGS) -c $*.c
 
# This is the rule that make runs by default when you type make.
# It will make the desired version of the benchmark and run it.
# The results are collected in a file whose name starts mf, dmf,
# cmf, dcmf (for the four possible targets) followed by the 
# file name specified in the input file.
basic: mflops
	mflops
# 
all: mflops dmflops cmflops dcmflops
 
#
#       Fortran Version
mflops : mflops.o kernel.o clksec.o
	$(FC) ${FFLAGS} -o mflops mflops.o kernel.o clksec.o
dmflops : dmflops.o dkernel.o dclksec.o
	$(FC) ${FFLAGS} -o dmflops dmflops.o dkernel.o dclksec.o
#
#       C version of Kernel routine.
cmflops : mflops.o cernel.o clksec.o
	$(FC) ${FFLAGS} -o cmflops mflops.o cernel.o clksec.o ${LDFLAGS}
dcmflops : dmflops.o dcernel.o dclksec.o
	$(FC) ${FFLAGS} -o dcmflops dmflops.o dcernel.o dclksec.o ${LDFLAGS}
#
mflops.o : mflops.f
	$(FC) ${FBASE} -c mflops.f
clksec.o : clksec.f
	$(FC) ${FBASE} -c clksec.f
kernel.o : kernel.f
	$(FC) ${FFLAGS} -c kernel.f
cernel.o : cernel.c
	$(CC) ${CFLAGS} -DREAL=FLOAT -c cernel.c
dcernel.o : cernel.c
	cp cernel.c dcernel.c
	$(CC) ${CFLAGS} -DREAL=double -c dcernel.c
dmflops.o : mflops.f
	sed -e 's/cANSI/     /g' -e 's/mflops.out/dmflops.out/' mflops.f > dmflops.f
	$(FC) ${FBASE} -c dmflops.f
dkernel.o : kernel.f
	sed -e 's/cANSI/     /g' kernel.f > dkernel.f
	$(FC) ${FFLAGS} -c dkernel.f
dclksec.o : clksec.f
	sed -e 's/cANSI/     /g' -e 's/mflops.out/dmflops.out/' clksec.f > dclksec.f
	$(FC) ${FBASE} -c dclksec.f
#
# Clean up the junk, but leave the results...
clean:
	rm -f *.o mflops dmflops cmflops dcmflops
	rm -f dmflops.f dkernel.f dcernel.c dclksec.f gmon.out core
#
# Clean up everything
pristine:
	rm -f *.o mflops dmflops cmflops dcmflops input dmf.* mf.* Makefile
	rm -f dmflops.f dkernel.f dcernel.c dclksec.f gmon.out core

