# MHD-Simulation-for-solar-and-stellar-winds

The codes are written in Fortran90. To compile, please type


#> Make


on terminal. The default compiler is gfortran. If you use a different compiler, change FC in Makefile. After the
codes are complied successfully, the executable file, “testrun” is created. Then, you can run the simulation for solar
winds via


#> ./testrun


If you like to run your job in the background, you can type


#> ./testrun < /dev/null > tmp.log&

The simulation results are printed out in printvars.dat in binary format. The data are printed out every dt=0.01.
See Table for the detail of the variables.
I put some plotting codes in “Codes”. If you can use fortran + pgplot, readdata_tev2.f90 is available.

#gfortran readdata_tev1.f90 -L/usr/local/lib -lpgplot

The library link (-L/usr/local/lib) may need to be modified, depending on your system.
If you are more familiar with python, readplot_sample2.py is available.

#python readplot_sample4.py
<img width="485" alt="Screenshot 2024-03-30 122858" src="https://github.com/mmyas/MHD-Simulation-for-solar-and-stellar-winds/assets/165456825/b9c83d78-c26f-4a20-995e-02390f1501ed">



(The status of the running job is output to tmp.log.)
