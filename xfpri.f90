!program xfpri
!!
!!   Test driver for 'fpri': the fast evaluator of the prismatic gravitational field
!!
!implicit integer (i-n)
!implicit real*8 (a-h,o-z)
!a=2.5d0;b=3.1d0;c=5.0d0
!fx=3.d0/sqrt(50.d0)
!fy=4.d0/sqrt(50.d0)
!fz=5.d0/sqrt(50.d0)
!write (*,"(a25,1p3e25.15)") "a,b,c=",a,b,c
!write (*,"(a25,1p3e25.15)") "fx,fy,fz=",fx,fy,fz
!write (*,"(a25,2a25)") "R","V (fvpri)"
!write (*,"(25x,3a25)") "gX (fgpri)","gY (fgpri)","gZ (fgpri)"
!write (*,"(25x,3a25)") "GXX (fggpri)","GYY (fggpri)","GZZ (fggpri)"
!write (*,"(25x,3a25)") "GYZ (fggpri)","GZX (fggpri)","GXY (fggpri)"
!write (*,"(25x,a25)") "V (fvgpri)"
!write (*,"(25x,3a25)") "gX (fvgpri)","gY (fvgpri)","gZ (fvgpri)"
!write (*,"(25x,a25)") "V (fpri)"
!write (*,"(25x,3a25)") "gX (fpri)","gY (fpri)","gZ (fpri)"
!write (*,"(25x,3a25)") "GXX (fpri)","GYY (fpri)","GZZ (fpri)"
!write (*,"(25x,3a25)") "GYZ (fpri)","GZX (fpri)","GXY (fpri)"
!do j=-1,2
!    r=4.d0**j
!    x=r*fx
!    y=r*fy
!    z=r*fz
!    xm=x-a;xp=x+a
!    ym=y-b;yp=y+b
!    zm=z-c;zp=z+c
!    u=fvpri(xp,xm,yp,ym,zp,zm)
!    write (*,"(1p2e25.15)") r,u
!    call fgpri(xp,xm,yp,ym,zp,zm,ux,uy,uz)
!    write (*,"(25x,1p3e25.15)") ux,uy,uz
!    call fggpri(xp,xm,yp,ym,zp,zm,uxx,uxy,uzx,uyy,uyz,uzz)
!    write (*,"(25x,1p3e25.15)") uxx,uyy,uzz
!    write (*,"(25x,1p3e25.15)") uyz,uzx,uxy
!    call fvgpri(xp,xm,yp,ym,zp,zm,u,ux,uy,uz)
!    write (*,"(25x,1pe25.15)") u
!    write (*,"(25x,1p3e25.15)") ux,uy,uz
!    call fpri(xp,xm,yp,ym,zp,zm,u,ux,uy,uz,uxx,uxy,uzx,uyy,uyz,uzz)
!    write (*,"(25x,1pe25.15)") u
!    write (*,"(25x,1p3e25.15)") ux,uy,uz
!    write (*,"(25x,1p3e25.15)") uxx,uyy,uzz
!    write (*,"(25x,1p3e25.15)") uyz,uzx,uxy
!enddo
!stop
!end program xfpri
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8 function fvpri(xp,xm,yp,ym,zp,zm)
!
! Fast computation of gravitational potential of a uniform rectangular prism
! by full utilization of the addition thorems of arctangent and logarithm
!
! Reference: Fukushima, T, (2019) Geophys. J. Int'l, submitted
!   Fast computation of prismatic gravitational field
!
! Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
!
!     Inputs: xp,xm,yp,ym,zp,zm = the relative coordinates of
!                                 the evaluation point referred to
!                                 the vertices of the prism
!     Output: fvpri:       gravitational potential
!
implicit integer (i-n)
implicit real*8 (a-h,o-z)
parameter (PI=3.14159265358979324d0)
xp2=xp*xp;yp2=yp*yp;zp2=zp*zp
xm2=xm*xm;ym2=ym*ym;zm2=zm*zm
xpyp=xp*yp;ypzp=yp*zp;zpxp=zp*xp
xpym=xp*ym;ypzm=yp*zm;zpxm=zp*xm
xmyp=xm*yp;ymzp=ym*zp;zmxp=zm*xp
xmym=xm*ym;ymzm=ym*zm;zmxm=zm*xm
wpp=xp2+yp2
rppp=sqrt(wpp+zp2)
rppm=sqrt(wpp+zm2)
wpm=xp2+ym2
rpmp=sqrt(wpm+zp2)
rpmm=sqrt(wpm+zm2)
wmp=xm2+yp2
rmpp=sqrt(wmp+zp2)
rmpm=sqrt(wmp+zm2)
wmm=xm2+ym2
rmmp=sqrt(wmm+zp2)
rmmm=sqrt(wmm+zm2)
sappp=ypzp;tappp=xp*rppp
sappm=ypzm;tappm=xp*rppm
sapmp=ymzp;tapmp=xp*rpmp
sapmm=ymzm;tapmm=xp*rpmm
sampp=ypzp;tampp=xm*rmpp
sampm=ypzm;tampm=xm*rmpm
sammp=ymzp;tammp=xm*rmmp
sammm=ymzm;tammm=xm*rmmm
s1app=sappp*tappm-tappp*sappm
s1apm=sapmp*tapmm-tapmp*sapmm
s1amp=sampp*tampm-tampp*sampm
s1amm=sammp*tammm-tammp*sammm
tappptappm=tappp*tappm
tapmptapmm=tapmp*tapmm
tampptampm=tampp*tampm
tammptammm=tammp*tammm
t1app=sappp*sappm+tappptappm
t1apm=sapmp*sapmm+tapmptapmm
t1amp=sampp*sampm+tampptampm
t1amm=sammp*sammm+tammptammm
i1app=isign2(t1app*tappptappm,sappp*tappp)
i1apm=isign2(t1apm*tapmptapmm,sapmp*tapmp)
i1amp=isign2(t1amp*tampptampm,sampp*tampp)
i1amm=isign2(t1amm*tammptammm,sammp*tammp)
s2ap=s1app*t1apm-t1app*s1apm
s2am=s1amp*t1amm-t1amp*s1amm
t1appt1apm=t1app*t1apm
t1ampt1amm=t1amp*t1amm
t2ap=s1app*s1apm+t1appt1apm
t2am=s1amp*s1amm+t1ampt1amm
i2ap=i1app-i1apm+isign2(t2ap*t1appt1apm,s1app*t1app)
i2am=i1amp-i1amm+isign2(t2am*t1ampt1amm,s1amp*t1amp)
ap=arctan2(s2ap,t2ap)
am=arctan2(s2am,t2am)
sbppp=zpxp;tbppp=yp*rppp
sbppm=zpxm;tbppm=yp*rmpp
sbpmp=zmxp;tbpmp=yp*rppm
sbpmm=zmxm;tbpmm=yp*rmpm
sbmpp=zpxp;tbmpp=ym*rpmp
sbmpm=zpxm;tbmpm=ym*rmmp
sbmmp=zmxp;tbmmp=ym*rpmm
sbmmm=zmxm;tbmmm=ym*rmmm
s1bpp=sbppp*tbppm-tbppp*sbppm
s1bpm=sbpmp*tbpmm-tbpmp*sbpmm
s1bmp=sbmpp*tbmpm-tbmpp*sbmpm
s1bmm=sbmmp*tbmmm-tbmmp*sbmmm
tbppptbppm=tbppp*tbppm
tbpmptbpmm=tbpmp*tbpmm
tbmpptbmpm=tbmpp*tbmpm
tbmmptbmmm=tbmmp*tbmmm
t1bpp=sbppp*sbppm+tbppptbppm
t1bpm=sbpmp*sbpmm+tbpmptbpmm
t1bmp=sbmpp*sbmpm+tbmpptbmpm
t1bmm=sbmmp*sbmmm+tbmmptbmmm
i1bpp=isign2(t1bpp*tbppptbppm,sbppp*tbppp)
i1bpm=isign2(t1bpm*tbpmptbpmm,sbpmp*tbpmp)
i1bmp=isign2(t1bmp*tbmpptbmpm,sbmpp*tbmpp)
i1bmm=isign2(t1bmm*tbmmptbmmm,sbmmp*tbmmp)
s2bp=s1bpp*t1bpm-t1bpp*s1bpm
s2bm=s1bmp*t1bmm-t1bmp*s1bmm
t1bppt1bpm=t1bpp*t1bpm
t1bmpt1bmm=t1bmp*t1bmm
t2bp=s1bpp*s1bpm+t1bppt1bpm
t2bm=s1bmp*s1bmm+t1bmpt1bmm
i2bp=i1bpp-i1bpm+isign2(t2bp*t1bppt1bpm,s1bpp*t1bpp)
i2bm=i1bmp-i1bmm+isign2(t2bm*t1bmpt1bmm,s1bmp*t1bmp)
bp=arctan2(s2bp,t2bp)
bm=arctan2(s2bm,t2bm)
scppp=xpyp;tcppp=zp*rppp
scppm=xpym;tcppm=zp*rpmp
scpmp=xmyp;tcpmp=zp*rmpp
scpmm=xmym;tcpmm=zp*rmmp
scmpp=xpyp;tcmpp=zm*rppm
scmpm=xpym;tcmpm=zm*rpmm
scmmp=xmyp;tcmmp=zm*rmpm
scmmm=xmym;tcmmm=zm*rmmm
s1cpp=scppp*tcppm-tcppp*scppm
s1cpm=scpmp*tcpmm-tcpmp*scpmm
s1cmp=scmpp*tcmpm-tcmpp*scmpm
s1cmm=scmmp*tcmmm-tcmmp*scmmm
tcppptcppm=tcppp*tcppm
tcpmptcpmm=tcpmp*tcpmm
tcmpptcmpm=tcmpp*tcmpm
tcmmptcmmm=tcmmp*tcmmm
t1cpp=scppp*scppm+tcppptcppm
t1cpm=scpmp*scpmm+tcpmptcpmm
t1cmp=scmpp*scmpm+tcmpptcmpm
t1cmm=scmmp*scmmm+tcmmptcmmm
i1cpp=isign2(t1cpp*tcppptcppm,scppp*tcppp)
i1cpm=isign2(t1cpm*tcpmptcpmm,scpmp*tcpmp)
i1cmp=isign2(t1cmp*tcmpptcmpm,scmpp*tcmpp)
i1cmm=isign2(t1cmm*tcmmptcmmm,scmmp*tcmmp)
s2cp=s1cpp*t1cpm-t1cpp*s1cpm
s2cm=s1cmp*t1cmm-t1cmp*s1cmm
t1cppt1cpm=t1cpp*t1cpm
t1cmpt1cmm=t1cmp*t1cmm
t2cp=s1cpp*s1cpm+t1cppt1cpm
t2cm=s1cmp*s1cmm+t1cmpt1cmm
i2cp=i1cpp-i1cpm+isign2(t2cp*t1cppt1cpm,s1cpp*t1cpp)
i2cm=i1cmp-i1cmm+isign2(t2cm*t1cmpt1cmm,s1cmp*t1cmp)
cp=arctan2(s2cp,t2cp)
cm=arctan2(s2cm,t2cm)
dyzap=ap+dble(i2ap)*PI
dyzam=am+dble(i2am)*PI
dzxbp=bp+dble(i2bp)*PI
dzxbm=bm+dble(i2bm)*PI
dxycp=cp+dble(i2cp)*PI
dxycm=cm+dble(i2cm)*PI
qxpp=argln3(xp,xm,yp2+zp2)
qxpm=argln3(xp,xm,yp2+zm2)
qxmp=argln3(xp,xm,ym2+zp2)
qxmm=argln3(xp,xm,ym2+zm2)
qypp=argln3(yp,ym,zp2+xp2)
qypm=argln3(yp,ym,zp2+xm2)
qymp=argln3(yp,ym,zm2+xp2)
qymm=argln3(yp,ym,zm2+xm2)
qzpp=argln3(zp,zm,xp2+yp2)
qzpm=argln3(zp,zm,xp2+ym2)
qzmp=argln3(zp,zm,xm2+yp2)
qzmm=argln3(zp,zm,xm2+ym2)
dxdpp=log(qxpp)
dxdpm=log(qxpm)
dxdmp=log(qxmp)
dxdmm=log(qxmm)
dyepp=log(qypp)
dyepm=log(qypm)
dyemp=log(qymp)
dyemm=log(qymm)
dzfpp=log(qzpp)
dzfpm=log(qzpm)
dzfmp=log(qzmp)
dzfmm=log(qzmm)
fvpri=-((xp2*dyzap-xm2*dyzam) &
    +(yp2*dzxbp-ym2*dzxbm) &
    +(zp2*dxycp-zm2*dxycm))*0.5d0 &
    +(ypzp*dxdpp-ypzm*dxdpm-ymzp*dxdmp+ymzm*dxdmm) &
    +(zpxp*dyepp-zpxm*dyepm-zmxp*dyemp+zmxm*dyemm) &
    +(xpyp*dzfpp-xpym*dzfpm-xmyp*dzfmp+xmym*dzfmm)
return;end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fgpri(xp,xm,yp,ym,zp,zm,ux,uy,uz)
!
! Fast computation of gravity vector of a uniform rectangular prism
! by full utilization of the addition thorems of arctangent and logarithm
!
! Reference: Fukushima, T, (2019) Geophys. J. Int'l, submitted
!   Fast computation of prismatic gravitational field
!
! Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
!
!     Inputs: xp,xm,yp,ym,zp,zm = the relative coordinates of
!                                 the evaluation point referred to
!                                 the vertices of the prism
!     Output: ux,uy,uz     gravitational acceleration vector, components
!
implicit integer (i-n)
implicit real*8 (a-h,o-z)
parameter (PI=3.14159265358979324d0)
xp2=xp*xp;yp2=yp*yp;zp2=zp*zp
xm2=xm*xm;ym2=ym*ym;zm2=zm*zm
xpyp=xp*yp;ypzp=yp*zp;zpxp=zp*xp
xpym=xp*ym;ypzm=yp*zm;zpxm=zp*xm
xmyp=xm*yp;ymzp=ym*zp;zmxp=zm*xp
xmym=xm*ym;ymzm=ym*zm;zmxm=zm*xm
wpp=xp2+yp2
rppp=sqrt(wpp+zp2)
rppm=sqrt(wpp+zm2)
wpm=xp2+ym2
rpmp=sqrt(wpm+zp2)
rpmm=sqrt(wpm+zm2)
wmp=xm2+yp2
rmpp=sqrt(wmp+zp2)
rmpm=sqrt(wmp+zm2)
wmm=xm2+ym2
rmmp=sqrt(wmm+zp2)
rmmm=sqrt(wmm+zm2)
sappp=ypzp;tappp=xp*rppp
sappm=ypzm;tappm=xp*rppm
sapmp=ymzp;tapmp=xp*rpmp
sapmm=ymzm;tapmm=xp*rpmm
sampp=ypzp;tampp=xm*rmpp
sampm=ypzm;tampm=xm*rmpm
sammp=ymzp;tammp=xm*rmmp
sammm=ymzm;tammm=xm*rmmm
s1app=sappp*tappm-tappp*sappm
s1apm=sapmp*tapmm-tapmp*sapmm
s1amp=sampp*tampm-tampp*sampm
s1amm=sammp*tammm-tammp*sammm
tappptappm=tappp*tappm
tapmptapmm=tapmp*tapmm
tampptampm=tampp*tampm
tammptammm=tammp*tammm
t1app=sappp*sappm+tappptappm
t1apm=sapmp*sapmm+tapmptapmm
t1amp=sampp*sampm+tampptampm
t1amm=sammp*sammm+tammptammm
i1app=isign2(t1app*tappptappm,sappp*tappp)
i1apm=isign2(t1apm*tapmptapmm,sapmp*tapmp)
i1amp=isign2(t1amp*tampptampm,sampp*tampp)
i1amm=isign2(t1amm*tammptammm,sammp*tammp)
s2ap=s1app*t1apm-t1app*s1apm
s2am=s1amp*t1amm-t1amp*s1amm
t1appt1apm=t1app*t1apm
t1ampt1amm=t1amp*t1amm
t2ap=s1app*s1apm+t1appt1apm
t2am=s1amp*s1amm+t1ampt1amm
i2ap=i1app-i1apm+isign2(t2ap*t1appt1apm,s1app*t1app)
i2am=i1amp-i1amm+isign2(t2am*t1ampt1amm,s1amp*t1amp)
ap=arctan2(s2ap,t2ap)
am=arctan2(s2am,t2am)
sbppp=zpxp;tbppp=yp*rppp
sbppm=zpxm;tbppm=yp*rmpp
sbpmp=zmxp;tbpmp=yp*rppm
sbpmm=zmxm;tbpmm=yp*rmpm
sbmpp=zpxp;tbmpp=ym*rpmp
sbmpm=zpxm;tbmpm=ym*rmmp
sbmmp=zmxp;tbmmp=ym*rpmm
sbmmm=zmxm;tbmmm=ym*rmmm
s1bpp=sbppp*tbppm-tbppp*sbppm
s1bpm=sbpmp*tbpmm-tbpmp*sbpmm
s1bmp=sbmpp*tbmpm-tbmpp*sbmpm
s1bmm=sbmmp*tbmmm-tbmmp*sbmmm
tbppptbppm=tbppp*tbppm
tbpmptbpmm=tbpmp*tbpmm
tbmpptbmpm=tbmpp*tbmpm
tbmmptbmmm=tbmmp*tbmmm
t1bpp=sbppp*sbppm+tbppptbppm
t1bpm=sbpmp*sbpmm+tbpmptbpmm
t1bmp=sbmpp*sbmpm+tbmpptbmpm
t1bmm=sbmmp*sbmmm+tbmmptbmmm
i1bpp=isign2(t1bpp*tbppptbppm,sbppp*tbppp)
i1bpm=isign2(t1bpm*tbpmptbpmm,sbpmp*tbpmp)
i1bmp=isign2(t1bmp*tbmpptbmpm,sbmpp*tbmpp)
i1bmm=isign2(t1bmm*tbmmptbmmm,sbmmp*tbmmp)
s2bp=s1bpp*t1bpm-t1bpp*s1bpm
s2bm=s1bmp*t1bmm-t1bmp*s1bmm
t1bppt1bpm=t1bpp*t1bpm
t1bmpt1bmm=t1bmp*t1bmm
t2bp=s1bpp*s1bpm+t1bppt1bpm
t2bm=s1bmp*s1bmm+t1bmpt1bmm
i2bp=i1bpp-i1bpm+isign2(t2bp*t1bppt1bpm,s1bpp*t1bpp)
i2bm=i1bmp-i1bmm+isign2(t2bm*t1bmpt1bmm,s1bmp*t1bmp)
bp=arctan2(s2bp,t2bp)
bm=arctan2(s2bm,t2bm)
scppp=xpyp;tcppp=zp*rppp
scppm=xpym;tcppm=zp*rpmp
scpmp=xmyp;tcpmp=zp*rmpp
scpmm=xmym;tcpmm=zp*rmmp
scmpp=xpyp;tcmpp=zm*rppm
scmpm=xpym;tcmpm=zm*rpmm
scmmp=xmyp;tcmmp=zm*rmpm
scmmm=xmym;tcmmm=zm*rmmm
s1cpp=scppp*tcppm-tcppp*scppm
s1cpm=scpmp*tcpmm-tcpmp*scpmm
s1cmp=scmpp*tcmpm-tcmpp*scmpm
s1cmm=scmmp*tcmmm-tcmmp*scmmm
tcppptcppm=tcppp*tcppm
tcpmptcpmm=tcpmp*tcpmm
tcmpptcmpm=tcmpp*tcmpm
tcmmptcmmm=tcmmp*tcmmm
t1cpp=scppp*scppm+tcppptcppm
t1cpm=scpmp*scpmm+tcpmptcpmm
t1cmp=scmpp*scmpm+tcmpptcmpm
t1cmm=scmmp*scmmm+tcmmptcmmm
i1cpp=isign2(t1cpp*tcppptcppm,scppp*tcppp)
i1cpm=isign2(t1cpm*tcpmptcpmm,scpmp*tcpmp)
i1cmp=isign2(t1cmp*tcmpptcmpm,scmpp*tcmpp)
i1cmm=isign2(t1cmm*tcmmptcmmm,scmmp*tcmmp)
s2cp=s1cpp*t1cpm-t1cpp*s1cpm
s2cm=s1cmp*t1cmm-t1cmp*s1cmm
t1cppt1cpm=t1cpp*t1cpm
t1cmpt1cmm=t1cmp*t1cmm
t2cp=s1cpp*s1cpm+t1cppt1cpm
t2cm=s1cmp*s1cmm+t1cmpt1cmm
i2cp=i1cpp-i1cpm+isign2(t2cp*t1cppt1cpm,s1cpp*t1cpp)
i2cm=i1cmp-i1cmm+isign2(t2cm*t1cmpt1cmm,s1cmp*t1cmp)
cp=arctan2(s2cp,t2cp)
cm=arctan2(s2cm,t2cm)
dyzap=ap+dble(i2ap)*PI
dyzam=am+dble(i2am)*PI
dzxbp=bp+dble(i2bp)*PI
dzxbm=bm+dble(i2bm)*PI
dxycp=cp+dble(i2cp)*PI
dxycm=cm+dble(i2cm)*PI
qxpp=argln3(xp,xm,yp2+zp2)
qxpm=argln3(xp,xm,yp2+zm2)
qxmp=argln3(xp,xm,ym2+zp2)
qxmm=argln3(xp,xm,ym2+zm2)
qypp=argln3(yp,ym,zp2+xp2)
qypm=argln3(yp,ym,zp2+xm2)
qymp=argln3(yp,ym,zm2+xp2)
qymm=argln3(yp,ym,zm2+xm2)
qzpp=argln3(zp,zm,xp2+yp2)
qzpm=argln3(zp,zm,xp2+ym2)
qzmp=argln3(zp,zm,xm2+yp2)
qzmm=argln3(zp,zm,xm2+ym2)
dxydp=log(qxpp/qxmp)
dxydm=log(qxpm/qxmm)
dzxdp=log(qxpp/qxpm)
dzxdm=log(qxmp/qxmm)
dyzep=log(qypp/qymp)
dyzem=log(qypm/qymm)
dxyep=log(qypp/qypm)
dxyem=log(qymp/qymm)
dzxfp=log(qzpp/qzmp)
dzxfm=log(qzpm/qzmm)
dyzfp=log(qzpp/qzpm)
dyzfm=log(qzmp/qzmm)
ux=-(xp*dyzap-xm*dyzam)+(zp*dxyep-zm*dxyem)+(yp*dzxfp-ym*dzxfm)
uy=-(yp*dzxbp-ym*dzxbm)+(xp*dyzfp-xm*dyzfm)+(zp*dxydp-zm*dxydm)
uz=-(zp*dxycp-zm*dxycm)+(yp*dzxdp-ym*dzxdm)+(xp*dyzep-xm*dyzem)
return;end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fggpri(xp,xm,yp,ym,zp,zm,uxx,uxy,uxz,uyy,uyz,uzz)
!
! Fast computation of gravity gradient tensor of a uniform rectangular prism
! by full utilization of the addition thorems of arctangent and logarithm
!
! Reference: Fukushima, T, (2019) Geophys. J. Int'l, submitted
!   Fast computation of prismatic gravitational field
!
! Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
!
!     Inputs: xp,xm,yp,ym,zp,zm = the relative coordinates of
!                                 the evaluation point referred to
!                                 the vertices of the prism
!     Output: uxx,uyy,uzz  gravity gradient tensor, diagonal components
!             uyz,uzx,uxy  gravity gradient tensor, non-diagonal components
!
implicit integer (i-n)
implicit real*8 (a-h,o-z)
parameter (PI=3.14159265358979324d0)
xp2=xp*xp;yp2=yp*yp;zp2=zp*zp
xm2=xm*xm;ym2=ym*ym;zm2=zm*zm
xpyp=xp*yp;ypzp=yp*zp;zpxp=zp*xp
xpym=xp*ym;ypzm=yp*zm;zpxm=zp*xm
xmyp=xm*yp;ymzp=ym*zp;zmxp=zm*xp
xmym=xm*ym;ymzm=ym*zm;zmxm=zm*xm
wpp=xp2+yp2
rppp=sqrt(wpp+zp2)
rppm=sqrt(wpp+zm2)
wpm=xp2+ym2
rpmp=sqrt(wpm+zp2)
rpmm=sqrt(wpm+zm2)
wmp=xm2+yp2
rmpp=sqrt(wmp+zp2)
rmpm=sqrt(wmp+zm2)
wmm=xm2+ym2
rmmp=sqrt(wmm+zp2)
rmmm=sqrt(wmm+zm2)
sappp=ypzp;tappp=xp*rppp
sappm=ypzm;tappm=xp*rppm
sapmp=ymzp;tapmp=xp*rpmp
sapmm=ymzm;tapmm=xp*rpmm
sampp=ypzp;tampp=xm*rmpp
sampm=ypzm;tampm=xm*rmpm
sammp=ymzp;tammp=xm*rmmp
sammm=ymzm;tammm=xm*rmmm
s1app=sappp*tappm-tappp*sappm
s1apm=sapmp*tapmm-tapmp*sapmm
s1amp=sampp*tampm-tampp*sampm
s1amm=sammp*tammm-tammp*sammm
tappptappm=tappp*tappm
tapmptapmm=tapmp*tapmm
tampptampm=tampp*tampm
tammptammm=tammp*tammm
t1app=sappp*sappm+tappptappm
t1apm=sapmp*sapmm+tapmptapmm
t1amp=sampp*sampm+tampptampm
t1amm=sammp*sammm+tammptammm
i1app=isign2(t1app*tappptappm,sappp*tappp)
i1apm=isign2(t1apm*tapmptapmm,sapmp*tapmp)
i1amp=isign2(t1amp*tampptampm,sampp*tampp)
i1amm=isign2(t1amm*tammptammm,sammp*tammp)
s2ap=s1app*t1apm-t1app*s1apm
s2am=s1amp*t1amm-t1amp*s1amm
t1appt1apm=t1app*t1apm
t1ampt1amm=t1amp*t1amm
t2ap=s1app*s1apm+t1appt1apm
t2am=s1amp*s1amm+t1ampt1amm
i2ap=i1app-i1apm+isign2(t2ap*t1appt1apm,s1app*t1app)
i2am=i1amp-i1amm+isign2(t2am*t1ampt1amm,s1amp*t1amp)
sbppp=zpxp;tbppp=yp*rppp
sbppm=zpxm;tbppm=yp*rmpp
sbpmp=zmxp;tbpmp=yp*rppm
sbpmm=zmxm;tbpmm=yp*rmpm
sbmpp=zpxp;tbmpp=ym*rpmp
sbmpm=zpxm;tbmpm=ym*rmmp
sbmmp=zmxp;tbmmp=ym*rpmm
sbmmm=zmxm;tbmmm=ym*rmmm
s1bpp=sbppp*tbppm-tbppp*sbppm
s1bpm=sbpmp*tbpmm-tbpmp*sbpmm
s1bmp=sbmpp*tbmpm-tbmpp*sbmpm
s1bmm=sbmmp*tbmmm-tbmmp*sbmmm
tbppptbppm=tbppp*tbppm
tbpmptbpmm=tbpmp*tbpmm
tbmpptbmpm=tbmpp*tbmpm
tbmmptbmmm=tbmmp*tbmmm
t1bpp=sbppp*sbppm+tbppptbppm
t1bpm=sbpmp*sbpmm+tbpmptbpmm
t1bmp=sbmpp*sbmpm+tbmpptbmpm
t1bmm=sbmmp*sbmmm+tbmmptbmmm
i1bpp=isign2(t1bpp*tbppptbppm,sbppp*tbppp)
i1bpm=isign2(t1bpm*tbpmptbpmm,sbpmp*tbpmp)
i1bmp=isign2(t1bmp*tbmpptbmpm,sbmpp*tbmpp)
i1bmm=isign2(t1bmm*tbmmptbmmm,sbmmp*tbmmp)
s2bp=s1bpp*t1bpm-t1bpp*s1bpm
s2bm=s1bmp*t1bmm-t1bmp*s1bmm
t1bppt1bpm=t1bpp*t1bpm
t1bmpt1bmm=t1bmp*t1bmm
t2bp=s1bpp*s1bpm+t1bppt1bpm
t2bm=s1bmp*s1bmm+t1bmpt1bmm
i2bp=i1bpp-i1bpm+isign2(t2bp*t1bppt1bpm,s1bpp*t1bpp)
i2bm=i1bmp-i1bmm+isign2(t2bm*t1bmpt1bmm,s1bmp*t1bmp)
scppp=xpyp;tcppp=zp*rppp
scppm=xpym;tcppm=zp*rpmp
scpmp=xmyp;tcpmp=zp*rmpp
scpmm=xmym;tcpmm=zp*rmmp
scmpp=xpyp;tcmpp=zm*rppm
scmpm=xpym;tcmpm=zm*rpmm
scmmp=xmyp;tcmmp=zm*rmpm
scmmm=xmym;tcmmm=zm*rmmm
s1cpp=scppp*tcppm-tcppp*scppm
s1cpm=scpmp*tcpmm-tcpmp*scpmm
s1cmp=scmpp*tcmpm-tcmpp*scmpm
s1cmm=scmmp*tcmmm-tcmmp*scmmm
tcppptcppm=tcppp*tcppm
tcpmptcpmm=tcpmp*tcpmm
tcmpptcmpm=tcmpp*tcmpm
tcmmptcmmm=tcmmp*tcmmm
t1cpp=scppp*scppm+tcppptcppm
t1cpm=scpmp*scpmm+tcpmptcpmm
t1cmp=scmpp*scmpm+tcmpptcmpm
t1cmm=scmmp*scmmm+tcmmptcmmm
i1cpp=isign2(t1cpp*tcppptcppm,scppp*tcppp)
i1cpm=isign2(t1cpm*tcpmptcpmm,scpmp*tcpmp)
i1cmp=isign2(t1cmp*tcmpptcmpm,scmpp*tcmpp)
i1cmm=isign2(t1cmm*tcmmptcmmm,scmmp*tcmmp)
s2cp=s1cpp*t1cpm-t1cpp*s1cpm
s2cm=s1cmp*t1cmm-t1cmp*s1cmm
t1cppt1cpm=t1cpp*t1cpm
t1cmpt1cmm=t1cmp*t1cmm
t2cp=s1cpp*s1cpm+t1cppt1cpm
t2cm=s1cmp*s1cmm+t1cmpt1cmm
i2cp=i1cpp-i1cpm+isign2(t2cp*t1cppt1cpm,s1cpp*t1cpp)
i2cm=i1cmp-i1cmm+isign2(t2cm*t1cmpt1cmm,s1cmp*t1cmp)
qxpp=argln3(xp,xm,yp2+zp2)
qxpm=argln3(xp,xm,yp2+zm2)
qxmp=argln3(xp,xm,ym2+zp2)
qxmm=argln3(xp,xm,ym2+zm2)
qypp=argln3(yp,ym,zp2+xp2)
qypm=argln3(yp,ym,zp2+xm2)
qymp=argln3(yp,ym,zm2+xp2)
qymm=argln3(yp,ym,zm2+xm2)
qzpp=argln3(zp,zm,xp2+yp2)
qzpm=argln3(zp,zm,xp2+ym2)
qzmp=argln3(zp,zm,xm2+yp2)
qzmm=argln3(zp,zm,xm2+ym2)
s3a=s2ap*t2am-t2ap*s2am
t3a=s2ap*s2am+t2ap*t2am
i3a=i2ap-i2am+isign2(t2ap*t2am*t3a,s2ap*t2ap)
s3b=s2bp*t2bm-t2bp*s2bm
t3b=s2bp*s2bm+t2bp*t2bm
i3b=i2bp-i2bm+isign2(t2bp*t2bm*t3b,s2bp*t2bp)
s3c=s2cp*t2cm-t2cp*s2cm
t3c=s2cp*s2cm+t2cp*t2cm
i3c=i2cp-i2cm+isign2(t2cp*t2cm*t3c,s2cp*t2cp)
d3a=arctan2(s3a,t3a)+dble(i3a)*PI
d3b=arctan2(s3b,t3b)+dble(i3b)*PI
d3c=arctan2(s3c,t3c)+dble(i3c)*PI
d3d=log((qxpp*qxmm)/(qxpm*qxmp))
d3e=log((qypp*qymm)/(qypm*qymp))
d3f=log((qzpp*qzmm)/(qzpm*qzmp))
uxx=-d3a
uyy=-d3b
uzz=-d3c
uyz=d3d
uxz=d3e
uxy=d3f
return;end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fvgpri(xp,xm,yp,ym,zp,zm,u,ux,uy,uz)
!
! Fast computation of gravitational field of a uniform rectangular prism
! by full utilization of the addition thorems of arctangent and logarithm
!
! Simplified version returning potential and gravity vector
!
! Reference: Fukushima, T, (2019) Geophys. J. Int'l, submitted
!   Fast computation of prismatic gravitational field
!
! Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
!
!     Inputs: xp,xm,yp,ym,zp,zm = the relative coordinates of
!                                 the evaluation point referred to
!                                 the vertices of the prism
!     Output: u:           gravitational potential
!             ux,uy,uz     gravitational acceleration vector, components
!
implicit integer (i-n)
implicit real*8 (a-h,o-z)
parameter (PI=3.14159265358979324d0)
xp2=xp*xp;yp2=yp*yp;zp2=zp*zp
xm2=xm*xm;ym2=ym*ym;zm2=zm*zm
xpyp=xp*yp;ypzp=yp*zp;zpxp=zp*xp
xpym=xp*ym;ypzm=yp*zm;zpxm=zp*xm
xmyp=xm*yp;ymzp=ym*zp;zmxp=zm*xp
xmym=xm*ym;ymzm=ym*zm;zmxm=zm*xm
wpp=xp2+yp2
rppp=sqrt(wpp+zp2)
rppm=sqrt(wpp+zm2)
wpm=xp2+ym2
rpmp=sqrt(wpm+zp2)
rpmm=sqrt(wpm+zm2)
wmp=xm2+yp2
rmpp=sqrt(wmp+zp2)
rmpm=sqrt(wmp+zm2)
wmm=xm2+ym2
rmmp=sqrt(wmm+zp2)
rmmm=sqrt(wmm+zm2)
sappp=ypzp;tappp=xp*rppp
sappm=ypzm;tappm=xp*rppm
sapmp=ymzp;tapmp=xp*rpmp
sapmm=ymzm;tapmm=xp*rpmm
sampp=ypzp;tampp=xm*rmpp
sampm=ypzm;tampm=xm*rmpm
sammp=ymzp;tammp=xm*rmmp
sammm=ymzm;tammm=xm*rmmm
s1app=sappp*tappm-tappp*sappm
s1apm=sapmp*tapmm-tapmp*sapmm
s1amp=sampp*tampm-tampp*sampm
s1amm=sammp*tammm-tammp*sammm
tappptappm=tappp*tappm
tapmptapmm=tapmp*tapmm
tampptampm=tampp*tampm
tammptammm=tammp*tammm
t1app=sappp*sappm+tappptappm
t1apm=sapmp*sapmm+tapmptapmm
t1amp=sampp*sampm+tampptampm
t1amm=sammp*sammm+tammptammm
i1app=isign2(t1app*tappptappm,sappp*tappp)
i1apm=isign2(t1apm*tapmptapmm,sapmp*tapmp)
i1amp=isign2(t1amp*tampptampm,sampp*tampp)
i1amm=isign2(t1amm*tammptammm,sammp*tammp)
s2ap=s1app*t1apm-t1app*s1apm
s2am=s1amp*t1amm-t1amp*s1amm
t1appt1apm=t1app*t1apm
t1ampt1amm=t1amp*t1amm
t2ap=s1app*s1apm+t1appt1apm
t2am=s1amp*s1amm+t1ampt1amm
i2ap=i1app-i1apm+isign2(t2ap*t1appt1apm,s1app*t1app)
i2am=i1amp-i1amm+isign2(t2am*t1ampt1amm,s1amp*t1amp)
ap=arctan2(s2ap,t2ap)
am=arctan2(s2am,t2am)
sbppp=zpxp;tbppp=yp*rppp
sbppm=zpxm;tbppm=yp*rmpp
sbpmp=zmxp;tbpmp=yp*rppm
sbpmm=zmxm;tbpmm=yp*rmpm
sbmpp=zpxp;tbmpp=ym*rpmp
sbmpm=zpxm;tbmpm=ym*rmmp
sbmmp=zmxp;tbmmp=ym*rpmm
sbmmm=zmxm;tbmmm=ym*rmmm
s1bpp=sbppp*tbppm-tbppp*sbppm
s1bpm=sbpmp*tbpmm-tbpmp*sbpmm
s1bmp=sbmpp*tbmpm-tbmpp*sbmpm
s1bmm=sbmmp*tbmmm-tbmmp*sbmmm
tbppptbppm=tbppp*tbppm
tbpmptbpmm=tbpmp*tbpmm
tbmpptbmpm=tbmpp*tbmpm
tbmmptbmmm=tbmmp*tbmmm
t1bpp=sbppp*sbppm+tbppptbppm
t1bpm=sbpmp*sbpmm+tbpmptbpmm
t1bmp=sbmpp*sbmpm+tbmpptbmpm
t1bmm=sbmmp*sbmmm+tbmmptbmmm
i1bpp=isign2(t1bpp*tbppptbppm,sbppp*tbppp)
i1bpm=isign2(t1bpm*tbpmptbpmm,sbpmp*tbpmp)
i1bmp=isign2(t1bmp*tbmpptbmpm,sbmpp*tbmpp)
i1bmm=isign2(t1bmm*tbmmptbmmm,sbmmp*tbmmp)
s2bp=s1bpp*t1bpm-t1bpp*s1bpm
s2bm=s1bmp*t1bmm-t1bmp*s1bmm
t1bppt1bpm=t1bpp*t1bpm
t1bmpt1bmm=t1bmp*t1bmm
t2bp=s1bpp*s1bpm+t1bppt1bpm
t2bm=s1bmp*s1bmm+t1bmpt1bmm
i2bp=i1bpp-i1bpm+isign2(t2bp*t1bppt1bpm,s1bpp*t1bpp)
i2bm=i1bmp-i1bmm+isign2(t2bm*t1bmpt1bmm,s1bmp*t1bmp)
bp=arctan2(s2bp,t2bp)
bm=arctan2(s2bm,t2bm)
scppp=xpyp;tcppp=zp*rppp
scppm=xpym;tcppm=zp*rpmp
scpmp=xmyp;tcpmp=zp*rmpp
scpmm=xmym;tcpmm=zp*rmmp
scmpp=xpyp;tcmpp=zm*rppm
scmpm=xpym;tcmpm=zm*rpmm
scmmp=xmyp;tcmmp=zm*rmpm
scmmm=xmym;tcmmm=zm*rmmm
s1cpp=scppp*tcppm-tcppp*scppm
s1cpm=scpmp*tcpmm-tcpmp*scpmm
s1cmp=scmpp*tcmpm-tcmpp*scmpm
s1cmm=scmmp*tcmmm-tcmmp*scmmm
tcppptcppm=tcppp*tcppm
tcpmptcpmm=tcpmp*tcpmm
tcmpptcmpm=tcmpp*tcmpm
tcmmptcmmm=tcmmp*tcmmm
t1cpp=scppp*scppm+tcppptcppm
t1cpm=scpmp*scpmm+tcpmptcpmm
t1cmp=scmpp*scmpm+tcmpptcmpm
t1cmm=scmmp*scmmm+tcmmptcmmm
i1cpp=isign2(t1cpp*tcppptcppm,scppp*tcppp)
i1cpm=isign2(t1cpm*tcpmptcpmm,scpmp*tcpmp)
i1cmp=isign2(t1cmp*tcmpptcmpm,scmpp*tcmpp)
i1cmm=isign2(t1cmm*tcmmptcmmm,scmmp*tcmmp)
s2cp=s1cpp*t1cpm-t1cpp*s1cpm
s2cm=s1cmp*t1cmm-t1cmp*s1cmm
t1cppt1cpm=t1cpp*t1cpm
t1cmpt1cmm=t1cmp*t1cmm
t2cp=s1cpp*s1cpm+t1cppt1cpm
t2cm=s1cmp*s1cmm+t1cmpt1cmm
i2cp=i1cpp-i1cpm+isign2(t2cp*t1cppt1cpm,s1cpp*t1cpp)
i2cm=i1cmp-i1cmm+isign2(t2cm*t1cmpt1cmm,s1cmp*t1cmp)
cp=arctan2(s2cp,t2cp)
cm=arctan2(s2cm,t2cm)
dyzap=ap+dble(i2ap)*PI
dyzam=am+dble(i2am)*PI
dzxbp=bp+dble(i2bp)*PI
dzxbm=bm+dble(i2bm)*PI
dxycp=cp+dble(i2cp)*PI
dxycm=cm+dble(i2cm)*PI
qxpp=argln3(xp,xm,yp2+zp2)
qxpm=argln3(xp,xm,yp2+zm2)
qxmp=argln3(xp,xm,ym2+zp2)
qxmm=argln3(xp,xm,ym2+zm2)
qypp=argln3(yp,ym,zp2+xp2)
qypm=argln3(yp,ym,zp2+xm2)
qymp=argln3(yp,ym,zm2+xp2)
qymm=argln3(yp,ym,zm2+xm2)
qzpp=argln3(zp,zm,xp2+yp2)
qzpm=argln3(zp,zm,xp2+ym2)
qzmp=argln3(zp,zm,xm2+yp2)
qzmm=argln3(zp,zm,xm2+ym2)
dxdpp=log(qxpp)
dxdpm=log(qxpm)
dxdmp=log(qxmp)
dxdmm=log(qxmm)
dyepp=log(qypp)
dyepm=log(qypm)
dyemp=log(qymp)
dyemm=log(qymm)
dzfpp=log(qzpp)
dzfpm=log(qzpm)
dzfmp=log(qzmp)
dzfmm=log(qzmm)
dxydp=dxdpp-dxdmp
dxydm=dxdpm-dxdmm
dzxdp=dxdpp-dxdpm
dzxdm=dxdmp-dxdmm
dyzep=dyepp-dyemp
dyzem=dyepm-dyemm
dxyep=dyepp-dyepm
dxyem=dyemp-dyemm
dzxfp=dzfpp-dzfmp
dzxfm=dzfpm-dzfmm
dyzfp=dzfpp-dzfpm
dyzfm=dzfmp-dzfmm
u=-((xp2*dyzap-xm2*dyzam) &
    +(yp2*dzxbp-ym2*dzxbm) &
    +(zp2*dxycp-zm2*dxycm))*0.5d0 &
    +(ypzp*dxdpp-ypzm*dxdpm-ymzp*dxdmp+ymzm*dxdmm) &
    +(zpxp*dyepp-zpxm*dyepm-zmxp*dyemp+zmxm*dyemm) &
    +(xpyp*dzfpp-xpym*dzfpm-xmyp*dzfmp+xmym*dzfmm)
ux=-(xp*dyzap-xm*dyzam)+(zp*dxyep-zm*dxyem)+(yp*dzxfp-ym*dzxfm)
uy=-(yp*dzxbp-ym*dzxbm)+(xp*dyzfp-xm*dyzfm)+(zp*dxydp-zm*dxydm)
uz=-(zp*dxycp-zm*dxycm)+(yp*dzxdp-ym*dzxdm)+(xp*dyzep-xm*dyzem)
return;end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fpri(xp,xm,yp,ym,zp,zm,u,ux,uy,uz,uxx,uxy,uxz,uyy,uyz,uzz)
!
! Fast computation of gravitational field of a uniform rectangular prism
! by full utilization of the addition thorems of arctangent and logarithm
!
! Full version returning potential, gravity, and gravity gradient     
!
! Reference: Fukushima, T, (2019) Geophys. J. Int'l, submitted
!   Fast computation of prismatic gravitational field
!
! Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
!
!     Inputs: xp,xm,yp,ym,zp,zm = the relative coordinates of
!                                 the evaluation point referred to
!                                 the vertices of the prism
!     Output: u:           gravitational potential
!             ux,uy,uz     gravitational acceleration vector, components
!             uxx,uyy,uzz  gravity gradient tensor, diagonal components
!             uyz,uzx,uxy  gravity gradient tensor, non-diagonal components
!
implicit integer (i-n)
implicit real*8 (a-h,o-z)
parameter (PI=3.14159265358979324d0)
xp2=xp*xp;yp2=yp*yp;zp2=zp*zp
xm2=xm*xm;ym2=ym*ym;zm2=zm*zm
xpyp=xp*yp;ypzp=yp*zp;zpxp=zp*xp
xpym=xp*ym;ypzm=yp*zm;zpxm=zp*xm
xmyp=xm*yp;ymzp=ym*zp;zmxp=zm*xp
xmym=xm*ym;ymzm=ym*zm;zmxm=zm*xm
wpp=xp2+yp2
rppp=sqrt(wpp+zp2)
rppm=sqrt(wpp+zm2)
wpm=xp2+ym2
rpmp=sqrt(wpm+zp2)
rpmm=sqrt(wpm+zm2)
wmp=xm2+yp2
rmpp=sqrt(wmp+zp2)
rmpm=sqrt(wmp+zm2)
wmm=xm2+ym2
rmmp=sqrt(wmm+zp2)
rmmm=sqrt(wmm+zm2)
sappp=ypzp;tappp=xp*rppp
sappm=ypzm;tappm=xp*rppm
sapmp=ymzp;tapmp=xp*rpmp
sapmm=ymzm;tapmm=xp*rpmm
sampp=ypzp;tampp=xm*rmpp
sampm=ypzm;tampm=xm*rmpm
sammp=ymzp;tammp=xm*rmmp
sammm=ymzm;tammm=xm*rmmm
s1app=sappp*tappm-tappp*sappm
s1apm=sapmp*tapmm-tapmp*sapmm
s1amp=sampp*tampm-tampp*sampm
s1amm=sammp*tammm-tammp*sammm
tappptappm=tappp*tappm
tapmptapmm=tapmp*tapmm
tampptampm=tampp*tampm
tammptammm=tammp*tammm
t1app=sappp*sappm+tappptappm
t1apm=sapmp*sapmm+tapmptapmm
t1amp=sampp*sampm+tampptampm
t1amm=sammp*sammm+tammptammm
i1app=isign2(t1app*tappptappm,sappp*tappp)
i1apm=isign2(t1apm*tapmptapmm,sapmp*tapmp)
i1amp=isign2(t1amp*tampptampm,sampp*tampp)
i1amm=isign2(t1amm*tammptammm,sammp*tammp)
s2ap=s1app*t1apm-t1app*s1apm
s2am=s1amp*t1amm-t1amp*s1amm
t1appt1apm=t1app*t1apm
t1ampt1amm=t1amp*t1amm
t2ap=s1app*s1apm+t1appt1apm
t2am=s1amp*s1amm+t1ampt1amm
i2ap=i1app-i1apm+isign2(t2ap*t1appt1apm,s1app*t1app)
i2am=i1amp-i1amm+isign2(t2am*t1ampt1amm,s1amp*t1amp)
ap=arctan2(s2ap,t2ap)
am=arctan2(s2am,t2am)
sbppp=zpxp;tbppp=yp*rppp
sbppm=zpxm;tbppm=yp*rmpp
sbpmp=zmxp;tbpmp=yp*rppm
sbpmm=zmxm;tbpmm=yp*rmpm
sbmpp=zpxp;tbmpp=ym*rpmp
sbmpm=zpxm;tbmpm=ym*rmmp
sbmmp=zmxp;tbmmp=ym*rpmm
sbmmm=zmxm;tbmmm=ym*rmmm
s1bpp=sbppp*tbppm-tbppp*sbppm
s1bpm=sbpmp*tbpmm-tbpmp*sbpmm
s1bmp=sbmpp*tbmpm-tbmpp*sbmpm
s1bmm=sbmmp*tbmmm-tbmmp*sbmmm
tbppptbppm=tbppp*tbppm
tbpmptbpmm=tbpmp*tbpmm
tbmpptbmpm=tbmpp*tbmpm
tbmmptbmmm=tbmmp*tbmmm
t1bpp=sbppp*sbppm+tbppptbppm
t1bpm=sbpmp*sbpmm+tbpmptbpmm
t1bmp=sbmpp*sbmpm+tbmpptbmpm
t1bmm=sbmmp*sbmmm+tbmmptbmmm
i1bpp=isign2(t1bpp*tbppptbppm,sbppp*tbppp)
i1bpm=isign2(t1bpm*tbpmptbpmm,sbpmp*tbpmp)
i1bmp=isign2(t1bmp*tbmpptbmpm,sbmpp*tbmpp)
i1bmm=isign2(t1bmm*tbmmptbmmm,sbmmp*tbmmp)
s2bp=s1bpp*t1bpm-t1bpp*s1bpm
s2bm=s1bmp*t1bmm-t1bmp*s1bmm
t1bppt1bpm=t1bpp*t1bpm
t1bmpt1bmm=t1bmp*t1bmm
t2bp=s1bpp*s1bpm+t1bppt1bpm
t2bm=s1bmp*s1bmm+t1bmpt1bmm
i2bp=i1bpp-i1bpm+isign2(t2bp*t1bppt1bpm,s1bpp*t1bpp)
i2bm=i1bmp-i1bmm+isign2(t2bm*t1bmpt1bmm,s1bmp*t1bmp)
bp=arctan2(s2bp,t2bp)
bm=arctan2(s2bm,t2bm)
scppp=xpyp;tcppp=zp*rppp
scppm=xpym;tcppm=zp*rpmp
scpmp=xmyp;tcpmp=zp*rmpp
scpmm=xmym;tcpmm=zp*rmmp
scmpp=xpyp;tcmpp=zm*rppm
scmpm=xpym;tcmpm=zm*rpmm
scmmp=xmyp;tcmmp=zm*rmpm
scmmm=xmym;tcmmm=zm*rmmm
s1cpp=scppp*tcppm-tcppp*scppm
s1cpm=scpmp*tcpmm-tcpmp*scpmm
s1cmp=scmpp*tcmpm-tcmpp*scmpm
s1cmm=scmmp*tcmmm-tcmmp*scmmm
tcppptcppm=tcppp*tcppm
tcpmptcpmm=tcpmp*tcpmm
tcmpptcmpm=tcmpp*tcmpm
tcmmptcmmm=tcmmp*tcmmm
t1cpp=scppp*scppm+tcppptcppm
t1cpm=scpmp*scpmm+tcpmptcpmm
t1cmp=scmpp*scmpm+tcmpptcmpm
t1cmm=scmmp*scmmm+tcmmptcmmm
i1cpp=isign2(t1cpp*tcppptcppm,scppp*tcppp)
i1cpm=isign2(t1cpm*tcpmptcpmm,scpmp*tcpmp)
i1cmp=isign2(t1cmp*tcmpptcmpm,scmpp*tcmpp)
i1cmm=isign2(t1cmm*tcmmptcmmm,scmmp*tcmmp)
s2cp=s1cpp*t1cpm-t1cpp*s1cpm
s2cm=s1cmp*t1cmm-t1cmp*s1cmm
t1cppt1cpm=t1cpp*t1cpm
t1cmpt1cmm=t1cmp*t1cmm
t2cp=s1cpp*s1cpm+t1cppt1cpm
t2cm=s1cmp*s1cmm+t1cmpt1cmm
i2cp=i1cpp-i1cpm+isign2(t2cp*t1cppt1cpm,s1cpp*t1cpp)
i2cm=i1cmp-i1cmm+isign2(t2cm*t1cmpt1cmm,s1cmp*t1cmp)
cp=arctan2(s2cp,t2cp)
cm=arctan2(s2cm,t2cm)
dyzap=ap+dble(i2ap)*PI
dyzam=am+dble(i2am)*PI
dzxbp=bp+dble(i2bp)*PI
dzxbm=bm+dble(i2bm)*PI
dxycp=cp+dble(i2cp)*PI
dxycm=cm+dble(i2cm)*PI
qxpp=argln3(xp,xm,yp2+zp2)
qxpm=argln3(xp,xm,yp2+zm2)
qxmp=argln3(xp,xm,ym2+zp2)
qxmm=argln3(xp,xm,ym2+zm2)
qypp=argln3(yp,ym,zp2+xp2)
qypm=argln3(yp,ym,zp2+xm2)
qymp=argln3(yp,ym,zm2+xp2)
qymm=argln3(yp,ym,zm2+xm2)
qzpp=argln3(zp,zm,xp2+yp2)
qzpm=argln3(zp,zm,xp2+ym2)
qzmp=argln3(zp,zm,xm2+yp2)
qzmm=argln3(zp,zm,xm2+ym2)
dxdpp=log(qxpp)
dxdpm=log(qxpm)
dxdmp=log(qxmp)
dxdmm=log(qxmm)
dyepp=log(qypp)
dyepm=log(qypm)
dyemp=log(qymp)
dyemm=log(qymm)
dzfpp=log(qzpp)
dzfpm=log(qzpm)
dzfmp=log(qzmp)
dzfmm=log(qzmm)
dxydp=dxdpp-dxdmp
dxydm=dxdpm-dxdmm
dzxdp=dxdpp-dxdpm
dzxdm=dxdmp-dxdmm
dyzep=dyepp-dyemp
dyzem=dyepm-dyemm
dxyep=dyepp-dyepm
dxyem=dyemp-dyemm
dzxfp=dzfpp-dzfmp
dzxfm=dzfpm-dzfmm
dyzfp=dzfpp-dzfpm
dyzfm=dzfmp-dzfmm
d3a=dyzap-dyzam
d3b=dzxbp-dzxbm
d3c=dxycp-dxycm
d3d=dxydp-dxydm
d3e=dyzep-dyzem
d3f=dzxfp-dzxfm
u=-((xp2*dyzap-xm2*dyzam) &
    +(yp2*dzxbp-ym2*dzxbm) &
    +(zp2*dxycp-zm2*dxycm))*0.5d0 &
    +(ypzp*dxdpp-ypzm*dxdpm-ymzp*dxdmp+ymzm*dxdmm) &
    +(zpxp*dyepp-zpxm*dyepm-zmxp*dyemp+zmxm*dyemm) &
    +(xpyp*dzfpp-xpym*dzfpm-xmyp*dzfmp+xmym*dzfmm)
ux=-(xp*dyzap-xm*dyzam)+(zp*dxyep-zm*dxyem)+(yp*dzxfp-ym*dzxfm)
uy=-(yp*dzxbp-ym*dzxbm)+(xp*dyzfp-xm*dyzfm)+(zp*dxydp-zm*dxydm)
uz=-(zp*dxycp-zm*dxycm)+(yp*dzxdp-ym*dzxdm)+(xp*dyzep-xm*dyzem)
uxx=-d3a
uyy=-d3b
uzz=-d3c
uyz=d3d
uxz=d3e
uxy=d3f
return;end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8 function arctan2(x,y)
!
! 2-argument arctangent function: double precision version
!
! Reference: Fukushima, T, (2019) Geophys. J. Int'l, submitted
!   Fast computation of prismatic gravitational field
!
! Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
!
real*8 PIHALF
parameter (PIHALF=1.5707963267948966d0)
real*8 x,y
if(y.ne.0.d0) then  ! avoid zero division
    arctan2=atan(x/y)  ! compute by definition
elseif(x.gt.0.d0) then
    arctan2=PIHALF
elseif(x.lt.0.d0) then
    arctan2=-PIHALF
else
    arctan2=0.d0
endif
return;end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8 function argln3(x,y,w)
!
! Argument function for the logarithm: double precision version
!
! Reference: Fukushima, T, (2019) Geophys. J. Int'l, submitted
!   Fast computation of prismatic gravitational field
!
! Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
!
real*8 x,y,w
if(y.gt.0.d0) then
    argln3=(x+sqrt(x*x+w))/(y+sqrt(y*y+w))
elseif(x.lt.0.d0) then
    argln3=(-y+sqrt(y*y+w))/(-x+sqrt(x*x+w))
else
    argln3=(x+sqrt(x*x+w))*(-y+sqrt(y*y+w))/(w+1.d-150)    
endif
return;end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function isign2(x,y)
!
! 2-argument sign function: integer version
!
! Reference: Fukushima, T, (2019) Geophys. J. Int'l, submitted
!   Fast computation of prismatic gravitational field
!
! Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
!
real*8 x,y
if(x.gt.0.d0) then
    isign2=0
elseif(y.gt.0.d0) then
    isign2=1
elseif(y.lt.0.d0) then
    isign2=-1
else
    isign2=0
endif
return;end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! a,b,c=    2.500000000000000E+00    3.100000000000000E+00    5.000000000000000E+00
! fx,fy,fz=    4.242640687119285E-01    5.656854249492380E-01    7.071067811865475E-01
! R                V (fvpri)
! gX (fgpri)               gY (fgpri)               gZ (fgpri)
! GXX (fggpri)             GYY (fggpri)             GZZ (fggpri)
! GYZ (fggpri)             GZX (fggpri)             GXY (fggpri)
! V (fvgpri)
! gX (fvgpri)              gY (fvgpri)              gZ (fvgpri)
! V (fpri)
! gX (fpri)                gY (fpri)                gZ (fpri)
! GXX (fpri)               GYY (fpri)               GZZ (fpri)
! GYZ (fpri)               GZX (fpri)               GXY (fpri)
!     2.500000000000000E-01    1.045493502522513E+02
! -6.528861883079162E-01   -6.363485275629387E-01   -3.363222353307181E-01
! -6.157715834856090E+00   -4.503501232741651E+00   -1.905153546761432E+00
! 2.882399480365517E-03    2.575611821341042E-03    6.906167105710924E-03
! 1.045493502522513E+02
! -6.528861883079149E-01   -6.363485275629397E-01   -3.363222353307187E-01
! 1.045493502522513E+02
! -6.528861883079149E-01   -6.363485275629397E-01   -3.363222353307187E-01
! -6.157715834856090E+00   -4.503501232741651E+00   -1.905153546761432E+00
! 2.882399480365505E-03    2.575611821340917E-03    6.906167105710903E-03
!     1.000000000000000E+00    1.029165650223841E+02
! -2.575795929629697E+00   -2.524382543341331E+00   -1.338619757450026E+00
! -6.108157515307118E+00   -4.523852961875530E+00   -1.934360137176524E+00
! 4.629304800468847E-02    4.148802049974923E-02    1.104316094654552E-01
! 1.029165650223841E+02
! -2.575795929629697E+00   -2.524382543341332E+00   -1.338619757450027E+00
! 1.029165650223841E+02
! -2.575795929629697E+00   -2.524382543341332E+00   -1.338619757450027E+00
! -6.108157515307118E+00   -4.523852961875530E+00   -1.934360137176524E+00
! 4.629304800468859E-02    4.148802049974920E-02    1.104316094654552E-01
!     4.000000000000000E+00    7.890362388224828E+01
! -8.038756369499978E+00   -8.720371512433630E+00   -4.742779493682749E+00
! -5.523742173017988E+00   -4.908871520157275E+00   -2.133756921183909E+00
! 7.799700525764759E-01    7.036194800263724E-01    1.883813127813028E+00
! 7.890362388224828E+01
! -8.038756369499978E+00   -8.720371512433630E+00   -4.742779493682749E+00
! 7.890362388224828E+01
! -8.038756369499978E+00   -8.720371512433630E+00   -4.742779493682749E+00
! -5.523742173017988E+00   -4.908871520157276E+00   -2.133756921183909E+00
! 7.799700525764759E-01    7.036194800263725E-01    1.883813127813029E+00
!     1.600000000000000E+01    1.949303202244659E+01
! -5.448977961438974E-01   -7.163264770343847E-01   -8.441665747134663E-01
! -3.225672642878615E-02    2.191649826687918E-03    3.006507660209839E-02
! 9.288204846998578E-02    7.126681559988013E-02    6.262258276887225E-02
! 1.949303202244659E+01
! -5.448977961438991E-01   -7.163264770343840E-01   -8.441665747134663E-01
! 1.949303202244659E+01
! -5.448977961438991E-01   -7.163264770343840E-01   -8.441665747134663E-01
! -3.225672642878614E-02    2.191649826687916E-03    3.006507660209838E-02
! 9.288204846998574E-02    7.126681559988024E-02    6.262258276887223E-02
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
