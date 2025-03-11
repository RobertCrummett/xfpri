### Source

The forward models are codes writen by Dr Toshio Fukushima

[https://www.researchgate.net/profile/Toshio-Fukushima](https://www.researchgate.net/profile/Toshio-Fukushima)

### Running Fortran Subroutines in C - Intel Compilers

Linking with Intel oneAPI Compiler Toolchain
```console
ifx /c xfpri.f90
cl main.c xfpri.obj
```

In C, the Fortran routines appear to be upper capitalized.
For example, to declare the function `fvpri` in C,
```c
extern double FVPRI(double*,double*,double*,double*,double*,double*);
```
In declarations, do not forget that Fortran is a pass by reference langauge!


It took me a long time to figure out how to link `main.c` to `fxpri.obj`.
This is because I did not know how to look at the symbol table on Windows.
To examine the symbols in `xfpri.obj` (or any `.obj` file for that matter), use the `dumpbin` utility:
```console
dumpbin /SYMBOLS xfpri.obj
```
The declaration in C must match the name in the dump *exactly*.
It is case sensitive.

### Running Fortran Subroutines in C - GNU Compilers

This process went more smoothly.
```console
gfortran -c xfpri.f90
gcc -c main.c
gcc -o app main.o xfpri.o
```

I did have to remove the program from the Fortran code, and change the declaration of the function in C
```c
extern double fvpri_(double*,double*,double*,double*,double*,double*);
```
Not a big issue.
I made sure to keep the Intel versions in this repository, because they were slightly more difficult to figure out.

The equivalent of `dumpbin` on unix platforms is the `nm` utility.
```console
nm xfpri.o
```
