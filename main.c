#include <math.h>
#include <stdio.h>

#ifndef _WIN32
#define FVPRI fvpri_
#endif

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
