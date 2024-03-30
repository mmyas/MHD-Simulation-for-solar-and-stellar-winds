# MHD-Simulation-for-solar-and-stellar-winds

The codes are written in Fortran90. To compile, please type
> Make
on terminal. The default compiler is gfortran. If you use a different compiler, change FC in Makefile. After the
codes are complied successfully, the executable file, “testrun” is created. Then, you can run the simulation for solar
winds via
> ./testrun
If you like to run your job in the background, you can type
> ./testrun < /dev/null > tmp.log&
(The status of the running job is output to tmp.log.)
