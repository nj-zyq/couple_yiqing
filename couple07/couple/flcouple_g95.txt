del couple_g95.exe
rem
rem Compile some recently modified subroutines, update
rem the library of subrotines libcouple.a, and compile
rem and link couple.for to create a new couple.exe
rem
REM g95 -c -O2 -march=pentium4 input.for
REM ar -rcvf libcouple.a input.o
rem 
REM g95 -c -O2 -march=pentium4 bathy.for
REM ar -rcvf libcouple.a bathy.o
rem 
REM g95 -c -O2 -march=pentium4 rwprof.for
REM ar -rcvf libcouple.a rwprof.o
rem 
REM g95 -c -O2 -march=pentium4 cwns.for
REM ar -rcvf libcouple.a cwns.o
rem
REM g95 -c -O2 -march=pentium4 bounds.for
REM ar -rcvf libcouple.a bounds.o
rem
REM g95 -c -O2 -march=pentium4 botdeps.for
REM ar -rcvf libcouple.a botdeps.o
rem 
REM g95 -c -O2 -march=pentium4 amat.for
REM ar -rcvf libcouple.a amat.o
rem 
REM g95 -c -O2 -march=pentium4 galrkn.for
REM ar -rcvf libcouple.a galrkn.o
rem 
g95 -c -O2 -march=pentium4 header.for
ar -rcvf libcouple.a galrkn.o
rem 
g95 -o couple_g95 -O2 -march=pentium4 couple_g95.for -L. libcouple.a
