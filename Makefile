FC=gfortran
#FC=ifort
FCFLAGS=-O4 
#FCFLAGS=-g -fbounds-check
#FCFLAGS=-O4 -L/usr/local/pgplot-intel -L/usr/X11R6/lib -lpgplot -lX11
#FCFLAGS=-O4 -L/usr/local/pgplot -L/usr/X11R6/lib -lpgplot -lX11
#FCFLAGS=-g -fbounds-check -L/opt/local/lib -L/usr/X11/lib -lpgplot -lX11
#FCFLAGS=-g -fbounds-check -L/usr/local/lib -lpgplot 
#FCFLAGS=-g -fbounds-check -L/usr/local/lib -lpgplot -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan

OBJS=nrtype.o vars_wvwnd24.o MHD_wvwnd24.o RadCondResBd_wvwnd_semifullT24.o

testrun: $(OBJS) 
	$(FC) $(FCFLAGS) $(OBJS) -o testrun

nrtype.mod: nrtype.o
nrtype.o: nrtype.f90
	$(FC) $(FCFLAGS) -c nrtype.f90

vars.mod: nrtype.mod vars_wvwnd24.f90
vars_wvwnd24.o: nrtype.mod vars_wvwnd24.f90
	$(FC) $(FCFLAGS) -c vars_wvwnd24.f90

MHD_wvwnd24.o: nrtype.mod vars.mod MHD_wvwnd24.f90
	$(FC) $(FCFLAGS) -c MHD_wvwnd24.f90

RadCondResBd_wvwnd_semifullT24.o: nrtype.mod vars.mod RadCondResBd_wvwnd_semifullT24.f90
	$(FC) $(FCFLAGS) -c RadCondResBd_wvwnd_semifullT24.f90

clean:
	@rm -rf testrun *.o
