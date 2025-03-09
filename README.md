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
To example the symbols in `xfpri.obj` (or any `.obj` file for that matter), use the `dumpbin` utility:
```console
dumpbin /SYMBOLS xfpri.obj
```
The declaration in C must match the name in the dump *exactly*.
It is case sensitive.
