# couple_yiqing

This is a file to compile COUPLE (an ocean acoustic model) successfully.

## Steps for compiling

- Change the compiler from f77 to gfortran

    ```
    $ gfortran  --version
    GNU Fortran (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0
    Copyright (C) 2019 Free Software Foundation, Inc.
    This is free software; see the source for copying conditions.  There is NO
    warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    ```

- Make sure all the file names are lowercase.for

- Make sure there are no extra spaces in makefile and the name of the files are correct


## Steps for running COUPLE

- Copy the executable files 'test.dat' and  'couple' to another folder '../test'

- There is a really silly thing in COUPLE... The name of the input file must be 'COUPLE.DAT' (can be seen in my comment in 'couple.for'). Thus we should first change the name of the input file from 'test1.dat' to 'COUPLE.DAT', then can we do the calculation. 

    Use these lines to run couple

    ```
    $ cp miao_input.dat COUPLE.DAT
    $ couple
    ```

- The output files are

    ```
    COUPLE.TL COUPLE.CPR COUPLE.PRT
    ```

    Then let's change the name of the files back.... Silly again...
    
    ```
    $ mv COUPLE.TL miao.TL
    $ mv COUPLE.CPR miao.CPR
    $ mv COUPLE.PRT miao.PRT
    ```

- Open you MATLAB, read the results and plot.