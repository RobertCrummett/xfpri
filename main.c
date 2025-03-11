/*
 * This is an example of how to use Fortran subroutines
 * in C code, on Intel and GNU compiler toolchains.
 */
#include <math.h>
#include <stdio.h>

/*
 * My Intel Compiler prefers to call the `fvpri` routine
 * `FVPRI`, while gcc prefers to call the `fvpri` routine
 * as `fvpri_`. These idiosyncrasies in compilers are
 * addressed via the header guard and redefinition.
 */
#ifndef _WIN32
#define FVPRI fvpri_
#endif

/*
 * Fortran is a pass by reference language.
 * C is a pass by value language.
 * Fortran subroutines called from C ought to be
 * declared as pass by reference - which explains
 * the pointers everywhere.
 *
 * The value returned by the subroutine FVPRI
 * is passed by value.
 */
extern double FVPRI(double*,double*,double*,double*,double*,double*);

int main() {
	double a = 2.5;
	double b = 3.1;
	double c = 5.0;
	double fx = 3.0 / sqrt(50);
	double fy = 4.0 / sqrt(50);
	double fz = 5.0 / sqrt(50);

	printf("a,b,c=%lf,%lf,%lf\n",a,b,c);
	printf("fx,fy,fz=%lf,%lf,%lf\n",fx,fy,fz);

	for (int i = -1; i < 3; i++) {
		double r = pow(4, i);
		double x = r * fx;
		double y = r * fy;
		double z = r * fz;
		double xm = x - a;
		double xp = x + a;
		double ym = y - b;
		double yp = y + b;
		double zm = z - c;
		double zp = z + c;
		double u = FVPRI(&xp,&xm,&yp,&ym,&zp,&zm);
		printf("%lf %lf\n", r, u);
	}

	return 0;
}
