#include <math.h>
#include <stdio.h>

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

/*
program xfpri
!
!   Test driver for 'fpri': the fast evaluator of the prismatic gravitational field
!
implicit integer (i-n)
implicit real*8 (a-h,o-z)
a=2.5d0;b=3.1d0;c=5.0d0
fx=3.d0/sqrt(50.d0)
fy=4.d0/sqrt(50.d0)
fz=5.d0/sqrt(50.d0)
write (*,"(a25,1p3e25.15)") "a,b,c=",a,b,c
write (*,"(a25,1p3e25.15)") "fx,fy,fz=",fx,fy,fz
write (*,"(a25,2a25)") "R","V (fvpri)"
write (*,"(25x,3a25)") "gX (fgpri)","gY (fgpri)","gZ (fgpri)"
write (*,"(25x,3a25)") "GXX (fggpri)","GYY (fggpri)","GZZ (fggpri)"
write (*,"(25x,3a25)") "GYZ (fggpri)","GZX (fggpri)","GXY (fggpri)"
write (*,"(25x,a25)") "V (fvgpri)"
write (*,"(25x,3a25)") "gX (fvgpri)","gY (fvgpri)","gZ (fvgpri)"
write (*,"(25x,a25)") "V (fpri)"
write (*,"(25x,3a25)") "gX (fpri)","gY (fpri)","gZ (fpri)"
write (*,"(25x,3a25)") "GXX (fpri)","GYY (fpri)","GZZ (fpri)"
write (*,"(25x,3a25)") "GYZ (fpri)","GZX (fpri)","GXY (fpri)"
do j=-1,2
    r=4.d0**j
    x=r*fx
    y=r*fy
    z=r*fz
    xm=x-a;xp=x+a
    ym=y-b;yp=y+b
    zm=z-c;zp=z+c
    u=fvpri(xp,xm,yp,ym,zp,zm)
    write (*,"(1p2e25.15)") r,u
    call fgpri(xp,xm,yp,ym,zp,zm,ux,uy,uz)
    write (*,"(25x,1p3e25.15)") ux,uy,uz
    call fggpri(xp,xm,yp,ym,zp,zm,uxx,uxy,uzx,uyy,uyz,uzz)
    write (*,"(25x,1p3e25.15)") uxx,uyy,uzz
    write (*,"(25x,1p3e25.15)") uyz,uzx,uxy
    call fvgpri(xp,xm,yp,ym,zp,zm,u,ux,uy,uz)
    write (*,"(25x,1pe25.15)") u
    write (*,"(25x,1p3e25.15)") ux,uy,uz
    call fpri(xp,xm,yp,ym,zp,zm,u,ux,uy,uz,uxx,uxy,uzx,uyy,uyz,uzz)
    write (*,"(25x,1pe25.15)") u
    write (*,"(25x,1p3e25.15)") ux,uy,uz
    write (*,"(25x,1p3e25.15)") uxx,uyy,uzz
    write (*,"(25x,1p3e25.15)") uyz,uzx,uxy
enddo
stop
end program xfpri
*/
