!program xtnpri
!!
!! Test program package of "tnvpri", "tngpri". & "tnggpri"
!!
!! Computation of graviational potential of a uniform rectangular prism
!! by even order 0, 2, to 16 truncated Taylor expansion
!!
!! Reference: Fukushima, T, (2019) Geophys. J. Int'l, revised
!!   Taylor expansion of prismatic gravitational field
!!
!! Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
!!
!! Used programs: fabset, fcprep, tnvpri, tngpri, tnggpri
!!
!implicit integer (i-n)
!implicit real*8 (a-h,o-z)
!!
!common /cset0/gm
!common /cset2/c200,c020,c002
!common /cset4/c400,c040,c004,c022,c202,c220
!common /cset6/c600,c060,c006,c420,c042,c204,c402,c024,c240,c222
!common /cset8/c800,c080,c008,c620,c062,c206,c602,c026,c260, &
!              c044,c404,c440,c422,c242,c224
!common /csetA/cA00,c0A0,c00A,c820,c082,c208,c802,c028,c280, &
!              c640,c064,c406,c604,c046,c460,c622,c262,c226, &
!              c244,c424,c442
!common /csetC/cC00,c0C0,c00C,cA20,c0A2,c20A,cA02,c02A,c2A0, &
!              c840,c084,c408,c804,c048,c480,c822,c282,c228, &
!              c066,c606,c660,c642,c264,c426,c624,c246,c462,c444
!common /csetE/cE00,c0E0,c00E,cC20,c0C2,c20C,cC02,c02C,c2C0, &
!              cA40,c0A4,c40A,cA04,c04A,c4A0,cA22,c2A2,c22A, &
!              c860,c086,c608,c806,c068,c680,c842,c284,c428, &
!              c824,c248,c482,c266,c626,c662,c644,c464,c446
!common /csetG/cG00,c0G0,c00G,cE20,c0E2,c20E,cE02,c02E,c2E0, &
!              cC40,c0C4,c40C,cc04,c04C,c4c0,cC22,c2C2,c22C, &
!              cA60,c0A6,c60A,cA06,c06A,c6A0,cA42,c2A4,c42A, &
!              cA24,c24A,c4A2,c088,c808,c880,c862,c286,c628, &
!              c826,c268,c682,c844,c484,c448,c466,c646,c664
!!
!dx=5.0d0
!dy=6.2d0
!dz=10.d0
!!
!a=dx*0.5d0;b=dy*0.5d0;c=dz*0.5d0
!!
!write (*,"(a20,1p3e25.15)") "a,b,c=",a,b,c
!!
!n=16
!call fabset(n,a,b)
!call fcprep(n,c)
!!
!write (*,"(a20,1pe25.15)") "GM=",gm
!write (*,"(a20,1p3e25.15)") "c200,c020,c002=",c200,c020,c002
!write (*,"(a20,1p3e25.15)") "c400,c040,c004=",c400,c040,c004
!write (*,"(a20,1p3e25.15)") "c022,c202,c220=",c022,c202,c220
!write (*,"(a20,1p3e25.15)") "c600,c060,c006=",c600,c060,c006
!write (*,"(a20,1p3e25.15)") "c420,c042,c204=",c420,c042,c204
!write (*,"(a20,1p3e25.15)") "c240,c024,c402=",c240,c024,c402
!write (*,"(a20,1pe25.15)") "c222=",c222
!write (*,"(a20,1p3e25.15)") "c800,c080,c008=",c800,c080,c008
!write (*,"(a20,1p3e25.15)") "c620,c062,c206=",c620,c062,c206
!write (*,"(a20,1p3e25.15)") "c260,c026,c602=",c260,c026,c602
!write (*,"(a20,1p3e25.15)") "c044,c404,c440=",c044,c404,c440
!write (*,"(a20,1p3e25.15)") "c422,c242,c224=",c422,c242,c224
!write (*,"(a20,1p3e25.15)") "cA00,c0A0,c00A=",cA00,c0A0,c00A
!write (*,"(a20,1p3e25.15)") "c820,c082,c208=",c820,c082,c208
!write (*,"(a20,1p3e25.15)") "c280,c028,c802=",c280,c028,c802
!write (*,"(a20,1p3e25.15)") "c640,c064,c406=",c640,c064,c406
!write (*,"(a20,1p3e25.15)") "c460,c046,c604=",c460,c046,c604
!write (*,"(a20,1p3e25.15)") "c622,c262,c226=",c622,c262,c226
!write (*,"(a20,1p3e25.15)") "c244,c424,c442=",c244,c424,c442
!write (*,"(a20,1p3e25.15)") "cC00,c0C0,c00C=",cC00,c0C0,c00C
!write (*,"(a20,1p3e25.15)") "cA20,c0A2,c20A=",cA20,c0A2,c20A
!write (*,"(a20,1p3e25.15)") "c2A0,c02A,cA02=",c2A0,c02A,cA02
!write (*,"(a20,1p3e25.15)") "c840,c084,c408=",c840,c084,c408
!write (*,"(a20,1p3e25.15)") "c480,c048,c804=",c480,c048,c804
!write (*,"(a20,1p3e25.15)") "c822,c282,c228=",c822,c282,c228
!write (*,"(a20,1p3e25.15)") "c066,c606,c660=",c066,c606,c660
!write (*,"(a20,1p3e25.15)") "c642,c264,c426=",c642,c264,c426
!write (*,"(a20,1p3e25.15)") "c462,c246,c624=",c462,c246,c624
!write (*,"(a20,1pe25.15)") "c444=",c444
!write (*,"(a20,1p3e25.15)") "cE00,c0E0,c00E=",cE00,c0E0,c00E
!write (*,"(a20,1p3e25.15)") "cC20,c0C2,c20C=",cC20,c0C2,c20C
!write (*,"(a20,1p3e25.15)") "c2C0,c02C,cC02=",c2C0,c02C,cC02
!write (*,"(a20,1p3e25.15)") "cA40,c0A4,c40A=",cA40,c0A4,c40A
!write (*,"(a20,1p3e25.15)") "c4A0,c04A,cA04=",c4A0,c04A,cA04
!write (*,"(a20,1p3e25.15)") "cA22,c2A2,c22A=",cA22,c2A2,c22A
!write (*,"(a20,1p3e25.15)") "c860,c086,c608=",c860,c086,c608
!write (*,"(a20,1p3e25.15)") "c680,c068,c806=",c680,c068,c806
!write (*,"(a20,1p3e25.15)") "c842,c284,c428=",c842,c284,c428
!write (*,"(a20,1p3e25.15)") "c482,c248,c824=",c482,c248,c824
!write (*,"(a20,1p3e25.15)") "c266,c626,c662=",c266,c626,c662
!write (*,"(a20,1p3e25.15)") "c644,c464,c446=",c644,c464,c446
!write (*,"(a20,1p3e25.15)") "cG00,c0G0,c00G=",cG00,c0G0,c00G
!write (*,"(a20,1p3e25.15)") "cE20,c0E2,c20E=",cE20,c0E2,c20E
!write (*,"(a20,1p3e25.15)") "c2E0,c02E,cE02=",c2E0,c02E,cE02
!write (*,"(a20,1p3e25.15)") "cC40,c0C4,c40C=",cC40,c0C4,c40C
!write (*,"(a20,1p3e25.15)") "c4C0,c04C,cC04=",c4C0,c04C,cC04
!write (*,"(a20,1p3e25.15)") "cC22,c2C2,c22C=",cC22,c2C2,c22C
!write (*,"(a20,1p3e25.15)") "cA60,c0A6,c60A=",cA60,c0A6,c60A
!write (*,"(a20,1p3e25.15)") "c6A0,c06A,cA06=",c6A0,c06A,cA06
!write (*,"(a20,1p3e25.15)") "cA42,c2A4,c42A=",cA42,c2A4,c42A
!write (*,"(a20,1p3e25.15)") "c4A2,c24A,cA24=",c4A2,c24A,cA24
!write (*,"(a20,1p3e25.15)") "c088,c808,c880=",c088,c808,c880
!write (*,"(a20,1p3e25.15)") "c862,c286,c628=",c862,c286,c628
!write (*,"(a20,1p3e25.15)") "c682,c268,c826=",c682,c268,c826
!write (*,"(a20,1p3e25.15)") "c844,c484,c448=",c844,c484,c448
!write (*,"(a20,1p3e25.15)") "c466,c646,c664=",c466,c646,c664
!!
!x=30.d0
!y=40.d0
!z=50.d0
!!
!write (*,"(a20,1p3e25.15)") "X,Y,Z=",x,y,z
!!
!do n=0,8
!    v=tnvpri(n,x,y,z)
!    write (*,"(a15,i5,1pe25.15)") "order,V=",2*n,v
!enddo
!do n=0,8
!    call tngpri(n,x,y,z,v,gx,gy,gz)
!    write (*,"(a15,i5,1p3e25.15)") "order,g=",2*n,gx,gy,gz
!enddo
!do n=0,8
!    call tnggpri(n,x,y,z,v,gx,gy,gz,gxx,gxy,gxz,gyy,gyz,gzz)
!    write (*,"(a15,i5,1p6e25.15)") "order,Gamma=",2*n,gxx,gxy,gxz,gyy,gyz,gzz
!enddo
!!
!stop
!end program xtnpri
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fabset(n,a,b)
!
! Preparation of polynomial coefficients, C_{ijkn}, from a and b
! of a uniform rectangular prism
! by even order 0, 2, to 16 truncated Taylor expansion
!
! Reference: Fukushima, T, (2019) Geophys. J. Int'l, revised
!   Taylor expansion of prismatic gravitational field
!
! Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
!
implicit real*8 (a-z)
integer n
!
parameter (F6=1.d0/6.d0,F3=1.d0/3.d0,F120=1.d0/120.d0,F12=1.d0/12.d0)
parameter (F3d40=3.d0/40.d0,F5=1.d0/5.d0,F20=1.d0/20.d0,F4=1.d0/4.d0)
parameter (F3d5=3.d0/5.d0,F3d20=3.d0/20.d0,F112=1.d0/112.d0,F48=1.d0/48.d0)
parameter (F16=1.d0/16.d0,F5d112=5.d0/112.d0,F24=1.d0/24.d0,F2=1.d0/2.d0)
parameter (F7=1.d0/7.d0,F8=1.d0/8.d0,F15d112=15.d0/112.d0,F45d56=45.d0/56.d0)
parameter (F15d14=15.d0/14.d0,F15d8=15.d0/8.d0,F45d8=45.d0/8.d0)
parameter (F45d28=45.d0/28.d0,F5760=1.d0/5760.d0,F96=1.d0/96.d0)
parameter (F7d960=7.d0/960.d0,F5d96=5.d0/96.d0,F35d1152=35.d0/1152.d0)
parameter (F7d60=7.d0/60.d0,F2d3=2.d0/3.d0,F9=1.d0/9.d0,F1440=1.d0/1440.d0)
parameter (F7d480=7.d0/480.d0,F35d288=35.d0/288.d0,F35d36=35.d0/36.d0)
parameter (F14d9=14.d0/9.d0,F192=1.d0/192.d0,F5d32=5.d0/32.d0,F7d64=7.d0/64.d0)
parameter (F5d4=5.d0/4.d0,F35d12=35.d0/12.d0,F35d32=35.d0/32.d0)
parameter (F35d192=35.d0/192.d0,F35d6=35.d0/6.d0,F8448=1.d0/8448.d0)
parameter (F3d128=3.d0/128.d0,F128=1.d0/128.d0,F35d768=35.d0/768.d0)
parameter (F63d2816=63.d0/2816.d0,F384=1.d0/384.d0,F3d8=3.d0/8.d0)
parameter (F5d6=5.d0/6.d0,F11=1.d0/11.d0,F2816=1.d0/2816.d0,F256=1.d0/256.d0)
parameter (F9d128=9.d0/128.d0,F35d256=35.d0/256.d0,F315d2816=315.d0/2816.d0)
parameter (F21d64=21.d0/64.d0,F1575d1408=1575.d0/1408.d0,F5d256=5.d0/256.d0)
parameter (F9d32=9.d0/32.d0,F9d8=9.d0/8.d0,F45d22=45.d0/22.d0)
parameter (F7d4224=7.d0/4224.d0,F35d384=35.d0/384.d0,F21d128=21.d0/128.d0)
parameter (F7d128=7.d0/128.d0,F315d1408=315.d0/1408.d0,F7d384=7.d0/384.d0)
parameter (F175d32=175.d0/32.d0,F525d88=525.d0/88.d0,F7d192=7.d0/192.d0)
parameter (F35d24=35.d0/24.d0,F315d44=315.d0/44.d0,F7d704=7.d0/704.d0)
parameter (F63d64=63.d0/64.d0,F35d64=35.d0/64.d0,F1575d352=1575.d0/352.d0)
parameter (F21d32=21.d0/32.d0,F525d4=525.d0/4.d0,F315d22=315.d0/22.d0)
parameter (F35d1408=35.d0/1408.d0,F35d128=35.d0/128.d0,F315d128=315.d0/128.d0)
parameter (F105d128=105.d0/128.d0,F1575d88=1575.d0/88.d0)
parameter (F175d1408=175.d0/1408.d0,F175d64=175.d0/64.d0)
parameter (F105d32=105.d0/32.d0,F7875d128=7875.d0/128.d0)
parameter (F4725d704=4725.d0/704.d0,F93184=1.d0/93184.d0,F1536=1.d0/1536.d0)
parameter (F11d5120=11.d0/5120.d0,F33d1792=33.d0/1792.d0,F77d3072=77.d0/3072.d0)
parameter (F21d512=21.d0/512.d0,F231d13312=231.d0/13312.d0,F11d640=11.d0/640.d0)
parameter (F33d28=33.d0/28.d0,F11d24=11.d0/24.d0,F13=1.d0/13.d0)
parameter (F46592=1.d0/46592.d0,F11d2560=11.d0/2560.d0,F77d1536=77.d0/1536.d0)
parameter (F693d6656=693.d0/6656.d0,F11d1280=11.d0/1280.d0)
parameter (F2079d1664=2079.d0/1664.d0,F33d224=33.d0/224.d0)
parameter (F11d4=11.d0/4.d0,F33d13=33.d0/13.d0,F5d93184=5.d0/93184.d0)
parameter (F5d1536=5.d0/1536.d0,F11d1024=11.d0/1024.d0,F105d512=105.d0/512.d0)
parameter (F3465d13312=3465.d0/13312.d0,F5d768=5.d0/768.d0,F11d512=11.d0/512.d0)
parameter (F165d1792=165.d0/1792.d0,F385d3072=385.d0/3072.d0)
parameter (F17325d1664=17325.d0/1664.d0,F165d448=165.d0/448.d0)
parameter (F55d192=55.d0/192.d0,F15d4=15.d0/4.d0,F1485d104=1485.d0/104.d0)
parameter (F15d46592=15.d0/46592.d0,F5d512=5.d0/512.d0,F165d256=165.d0/256.d0)
parameter (F495d1792=495.d0/1792.d0,F385d512=385.d0/512.d0)
parameter (F315d512=315.d0/512.d0,F10395d1664=10395.d0/1664.d0)
parameter (F165d512=165.d0/512.d0,F495d112=495.d0/112.d0,F275d32=275.d0/32.d0)
parameter (F1485d4=1485.d0/4.d0,F1485d52=1485.d0/52.d0,F3328=1.d0/3328.d0)
parameter (F7d768=7.d0/768.d0,F77d1280=77.d0/1280.d0,F33d256=33.d0/256.d0)
parameter (F77d192=77.d0/192.d0,F1155d52=1155.d0/52.d0,F385d768=385.d0/768.d0)
parameter (F1155d256=1155.d0/256.d0,F1155d3328=1155.d0/3328.d0)
parameter (F5d3328=5.d0/3328.d0,F385d256=385.d0/256.d0,F105d256=105.d0/256.d0)
parameter (F10395d832=10395.d0/832.d0,F175d768=175.d0/768.d0)
parameter (F17325d416=17325.d0/416.d0,F385d192=385.d0/192.d0)
parameter (F3465d52=3465.d0/52.d0,F275d6656=275.d0/6656.d0)
parameter (F1925d768=1925.d0/768.d0,F1925d512=1925.d0/512.d0)
parameter (F825d256=825.d0/256.d0,F1925d1536=1925.d0/1536.d0)
parameter (F51975d64=51975.d0/64.d0,F51975d832=51975.d0/832.d0)
parameter (F30720=1.d0/30720.d0,F6144=1.d0/6144.d0,F91d30720=91.d0/30720.d0)
parameter (F143d30720=143.d0/30720.d0,F143d6144=143.d0/6144.d0)
parameter (F91d10240=91.d0/10240.d0,F77d2048=77.d0/2048.d0)
parameter (F143d10240=143.d0/10240.d0,F3072=1.d0/3072.d0)
parameter (F91d3840=91.d0/3840.d0,F143d1920=143.d0/1920.d0)
parameter (F143d48=143.d0/48.d0,F91d120=91.d0/120.d0,F7d6=7.d0/6.d0)
parameter (F15=1.d0/15.d0,F1001d10240=1001.d0/10240.d0)
parameter (F143d15360=143.d0/15360.d0,F91d5120=91.d0/5120.d0)
parameter (F7007d5120=7007.d0/5120.d0,F7d6144=7.d0/6144.d0)
parameter (F91d7680=91.d0/7680.d0,F143d3840=143.d0/3840.d0)
parameter (F143d96=143.d0/96.d0,F91d240=91.d0/240.d0,F7d12=7.d0/12.d0)
parameter (F91d30=91.d0/30.d0,F10240=1.d0/10240.d0,F2048=1.d0/2048.d0)
parameter (F273d10240=273.d0/10240.d0,F143d2048=143.d0/2048.d0)
parameter (F819d10240=819.d0/10240.d0,F231d2048=231.d0/2048.d0)
parameter (F3003d10240=3003.d0/10240.d0,F143d5120=143.d0/5120.d0)
parameter (F429d2048=429.d0/2048.d0,F1617d512=1617.d0/512.d0)
parameter (F21021d1280=21021.d0/1280.d0,F1024=1.d0/1024.d0)
parameter (F143d128=143.d0/128.d0,F91d320=91.d0/320.d0,F77d16=77.d0/16.d0)
parameter (F1001d40=1001.d0/40.d0,F5120=1.d0/5120.d0,F143d2560=143.d0/2560.d0)
parameter (F143d1024=143.d0/1024.d0,F273d5120=273.d0/5120.d0)
parameter (F231d1024=231.d0/1024.d0,F21021d2560=21021.d0/2560.d0)
parameter (F143d64=143.d0/64.d0,F91d160=91.d0/160.d0,F7007d8=7007.d0/8.d0)
parameter (F1001d20=1001.d0/20.d0,F35d2048=35.d0/2048.d0,F91d2048=91.d0/2048.d0)
parameter (F715d2048=715.d0/2048.d0,F385d2048=385.d0/2048.d0)
parameter (F1001d2048=1001.d0/2048.d0,F5d1024=5.d0/1024.d0)
parameter (F91d1024=91.d0/1024.d0,F143d512=143.d0/512.d0)
parameter (F637d1024=637.d0/1024.d0,F2695d256=2695.d0/256.d0)
parameter (F7007d128=7007.d0/128.d0,F5d2048=5.d0/2048.d0)
parameter (F715d256=715.d0/256.d0,F91d128=91.d0/128.d0,F385d32=385.d0/32.d0)
parameter (F1001d16=1001.d0/16.d0,F3d2048=3.d0/2048.d0)
parameter (F15d1024=15.d0/1024.d0,F273d2048=273.d0/2048.d0)
parameter (F429d1024=429.d0/1024.d0,F2145d2048=2145.d0/2048.d0)
parameter (F273d1024=273.d0/1024.d0,F1155d2048=1155.d0/2048.d0)
parameter (F21021d1024=21021.d0/1024.d0,F105d2048=105.d0/2048.d0)
parameter (F8085d512=8085.d0/512.d0,F21021d256=21021.d0/256.d0)
parameter (F15d2048=15.d0/2048.d0,F2145d256=2145.d0/256.d0)
parameter (F273d128=273.d0/128.d0,F1155d32=1155.d0/32.d0)
parameter (F3003d16=3003.d0/16.d0,F7d512=7.d0/512.d0,F35d512=35.d0/512.d0)
parameter (F637d512=637.d0/512.d0,F5005d512=5005.d0/512.d0)
parameter (F637d256=637.d0/256.d0,F2695d64=2695.d0/64.d0)
parameter (F7007d32=7007.d0/32.d0,F49d512=49.d0/512.d0)
parameter (F245d512=245.d0/512.d0,F245245d512=245245.d0/512.d0)
parameter (F7007d256=7007.d0/256.d0,F21d1024=21.d0/1024.d0)
parameter (F105d1024=105.d0/1024.d0,F1911d1024=1911.d0/1024.d0)
parameter (F3003d1024=3003.d0/1024.d0,F15015d1024=15015.d0/1024.d0)
parameter (F8085d256=8085.d0/256.d0,F21021d128=21021.d0/128.d0)
parameter (F9009d256=9009.d0/256.d0,F1911d512=1911.d0/512.d0)
parameter (F735735d128=735735.d0/128.d0,F21021d64=21021.d0/64.d0)
parameter (F5013504=1.d0/5013504.d0,F12288=1.d0/12288.d0)
parameter (F8192=1.d0/8192.d0,F13d12288=13.d0/12288.d0)
parameter (F143d147456=143.d0/147456.d0,F39d4096=39.d0/4096.d0)
parameter (F77d8192=77.d0/8192.d0,F143d4096=143.d0/4096.d0)
parameter (F6435d557056=6435.d0/557056.d0,F768=1.d0/768.d0)
parameter (F13d96=13.d0/96.d0,F143d576=143.d0/576.d0,F13d2=13.d0/2.d0)
parameter (F4d3=4.d0/3.d0,F17=1.d0/17.d0,F208896=1.d0/208896.d0)
parameter (F4096=1.d0/4096.d0,F3d4096=3.d0/4096.d0,F13d4096=13.d0/4096.d0)
parameter (F143d12288=143.d0/12288.d0,F117d4096=117.d0/4096.d0)
parameter (F3003d4096=3003.d0/4096.d0,F429d4096=429.d0/4096.d0)
parameter (F6435d69632=6435.d0/69632.d0,F21d4096=21.d0/4096.d0)
parameter (F1001d12288=1001.d0/12288.d0,F9009d4096=9009.d0/4096.d0)
parameter (F231d4096=231.d0/4096.d0,F6435d4352=6435.d0/4352.d0)
parameter (F3d256=3.d0/256.d0,F13d128=13.d0/128.d0,F715d96=715.d0/96.d0)
parameter (F39d16=39.d0/16.d0,F7d2=7.d0/2.d0,F60d17=60.d0/17.d0)
parameter (F417792=1.d0/417792.d0,F3d8192=3.d0/8192.d0)
parameter (F143d24576=143.d0/24576.d0,F1001d4096=1001.d0/4096.d0)
parameter (F45045d139264=45045.d0/139264.d0,F13d2048=13.d0/2048.d0)
parameter (F1573d24576=1573.d0/24576.d0,F1001d256=1001.d0/256.d0)
parameter (F105105d4352=105105.d0/4352.d0,F139264=1.d0/139264.d0)
parameter (F13d256=13.d0/256.d0,F715d768=715.d0/768.d0,F91d2=91.d0/2.d0)
parameter (F1365d34=1365.d0/34.d0,F5d208896=5.d0/208896.d0)
parameter (F5d4096=5.d0/4096.d0,F15d4096=15.d0/4096.d0)
parameter (F65d4096=65.d0/4096.d0,F3575d12288=3575.d0/12288.d0)
parameter (F6435d4096=6435.d0/4096.d0,F3465d4096=3465.d0/4096.d0)
parameter (F15015d4096=15015.d0/4096.d0,F45045d4352=45045.d0/4352.d0)
parameter (F55d208896=55.d0/208896.d0,F275d4096=275.d0/4096.d0)
parameter (F165d4096=165.d0/4096.d0,F715d64=715.d0/64.d0)
parameter (F3575d384=3575.d0/384.d0,F195d16=195.d0/16.d0)
parameter (F175d16=175.d0/16.d0,F1365d17=1365.d0/17.d0)
parameter (F11d626688=11.d0/626688.d0,F11d12288=11.d0/12288.d0)
parameter (F11d4096=11.d0/4096.d0,F1573d36864=1573.d0/36864.d0)
parameter (F77d4096=77.d0/4096.d0,F45045d69632=45045.d0/69632.d0)
parameter (F11d2048=11.d0/2048.d0,F1287d4096=1287.d0/4096.d0)
parameter (F231d256=231.d0/256.d0,F63063d544=63063.d0/544.d0)
parameter (F1573d2304=1573.d0/2304.d0,F77d48=77.d0/48.d0)
parameter (F1001d48=1001.d0/48.d0,F5005d34=5005.d0/34.d0)
parameter (F11d208896=11.d0/208896.d0,F33d4096=33.d0/4096.d0)
parameter (F1573d12288=1573.d0/12288.d0,F135135d4352=135135.d0/4352.d0)
parameter (F3003d64=3003.d0/64.d0,F315315d2176=315315.d0/2176.d0)
parameter (F7865d768=7865.d0/768.d0,F429d64=429.d0/64.d0)
parameter (F15015d34=15015.d0/34.d0,F33d278528=33.d0/278528.d0)
parameter (F99d8192=99.d0/8192.d0,F4719d16384=4719.d0/16384.d0)
parameter (F1287d256=1287.d0/256.d0,F3003d32=3003.d0/32.d0)
parameter (F225225d1088=225225.d0/1088.d0,F4719d4096=4719.d0/4096.d0)
parameter (F1155d8192=1155.d0/8192.d0,F75075d4096=75075.d0/4096.d0)
parameter (F225225d278528=225225.d0/278528.d0,F11d69632=11.d0/69632.d0)
parameter (F99d4096=99.d0/4096.d0,F1573d4096=1573.d0/4096.d0)
parameter (F1155d4096=1155.d0/4096.d0,F225225d4352=225225.d0/4352.d0)
parameter (F99d2048=99.d0/2048.d0,F315315d544=315315.d0/544.d0)
parameter (F1287d32=1287.d0/32.d0,F231d64=231.d0/64.d0)
parameter (F225225d272=225225.d0/272.d0,F55d139264=55.d0/139264.d0)
parameter (F495d8192=495.d0/8192.d0,F7865d8192=7865.d0/8192.d0)
parameter (F15015d256=15015.d0/256.d0,F1576575d4352=1576575.d0/4352.d0)
parameter (F96525d128=96525.d0/128.d0,F17325d128=17325.d0/128.d0)
parameter (F225225d8=225225.d0/8.d0,F675675d544=675675.d0/544.d0)
parameter (F77d104448=77.d0/104448.d0,F11011d6144=11011.d0/6144.d0)
parameter (F9009d2048=9009.d0/2048.d0,F1617d128=1617.d0/128.d0)
parameter (F315315d272=315315.d0/272.d0,F75075d2048=75075.d0/2048.d0)
parameter (F13475d2048=13475.d0/2048.d0,F175175d16=175175.d0/16.d0)
parameter (F525525d1088=525525.d0/1088.d0)
!
common /fcsetgm/gm1
common /fcset200/c2000,c2001
common /fcset020/c0200,c0201
common /fcset002/c0020,c0021
common /fcset400/c4000,c4001,c4002
common /fcset040/c0400,c0401,c0402
common /fcset004/c0040,c0041,c0042
common /fcset022/c0220,c0221,c0222
common /fcset202/c2020,c2021,c2022
common /fcset220/c2200,c2201,c2202
common /fcset600/c6000,c6001,c6002,c6003
common /fcset060/c0600,c0601,c0602,c0603
common /fcset006/c0060,c0061,c0062,c0063
common /fcset420/c4200,c4201,c4202,c4203
common /fcset042/c0420,c0421,c0422,c0423
common /fcset204/c2040,c2041,c2042,c2043
common /fcset402/c4020,c4021,c4022,c4023
common /fcset240/c2400,c2401,c2402,c2403
common /fcset024/c0240,c0241,c0242,c0243
common /fcset222/c2220,c2221,c2222,c2223
common /fcset800/c8000,c8001,c8002,c8003,c8004
common /fcset080/c0800,c0801,c0802,c0803,c0804
common /fcset008/c0080,c0081,c0082,c0083,c0084
common /fcset620/c6200,c6201,c6202,c6203,c6204
common /fcset062/c0620,c0621,c0622,c0623,c0624
common /fcset206/c2060,c2061,c2062,c2063,c2064
common /fcset602/c6020,c6021,c6022,c6023,c6024
common /fcset026/c0260,c0261,c0262,c0263,c0264
common /fcset260/c2600,c2601,c2602,c2603,c2604
common /fcset044/c0440,c0441,c0442,c0443,c0444
common /fcset404/c4040,c4041,c4042,c4043,c4044
common /fcset440/c4400,c4401,c4402,c4403,c4404
common /fcset422/c4220,c4221,c4222,c4223,c4224
common /fcset242/c2420,c2421,c2422,c2423,c2424
common /fcset224/c2240,c2241,c2242,c2243,c2244
common /fcsetA00/cA000,cA001,cA002,cA003,cA004,cA005
common /fcset0A0/c0A00,c0A01,c0A02,c0A03,c0A04,c0A05
common /fcset00A/c00A0,c00A1,c00A2,c00A3,c00A4,c00A5
common /fcset820/c8200,c8201,c8202,c8203,c8204,c8205
common /fcset082/c0820,c0821,c0822,c0823,c0824,c0825
common /fcset208/c2080,c2081,c2082,c2083,c2084,c2085
common /fcset802/c8020,c8021,c8022,c8023,c8024,c8025
common /fcset028/c0280,c0281,c0282,c0283,c0284,c0285
common /fcset280/c2800,c2801,c2802,c2803,c2804,c2805
common /fcset640/c6400,c6401,c6402,c6403,c6404,c6405
common /fcset064/c0640,c0641,c0642,c0643,c0644,c0645
common /fcset406/c4060,c4061,c4062,c4063,c4064,c4065
common /fcset604/c6040,c6041,c6042,c6043,c6044,c6045
common /fcset046/c0460,c0461,c0462,c0463,c0464,c0465
common /fcset460/c4600,c4601,c4602,c4603,c4604,c4605
common /fcset622/c6220,c6221,c6222,c6223,c6224,c6225
common /fcset262/c2620,c2621,c2622,c2623,c2624,c2625
common /fcset226/c2260,c2261,c2262,c2263,c2264,c2265
common /fcset244/c2440,c2441,c2442,c2443,c2444,c2445
common /fcset424/c4240,c4241,c4242,c4243,c4244,c4245
common /fcset442/c4420,c4421,c4422,c4423,c4424,c4425
common /fcsetC00/cC000,cC001,cC002,cC003,cC004,cC005,cC006
common /fcset0C0/c0C00,c0C01,c0C02,c0C03,c0C04,c0C05,c0C06
common /fcset00C/c00C0,c00C1,c00C2,c00C3,c00C4,c00C5,c00C6
common /fcsetA20/cA200,cA201,cA202,cA203,cA204,cA205,cA206
common /fcset0A2/c0A20,c0A21,c0A22,c0A23,c0A24,c0A25,c0A26
common /fcset20A/c20A0,c20A1,c20A2,c20A3,c20A4,c20A5,c20A6
common /fcsetA02/cA020,cA021,cA022,cA023,cA024,cA025,cA026
common /fcset2A0/c2A00,c2A01,c2A02,c2A03,c2A04,c2A05,c2A06
common /fcset02A/c02A0,c02A1,c02A2,c02A3,c02A4,c02A5,c02A6
common /fcset840/c8400,c8401,c8402,c8403,c8404,c8405,c8406
common /fcset084/c0840,c0841,c0842,c0843,c0844,c0845,c0846
common /fcset408/c4080,c4081,c4082,c4083,c4084,c4085,c4086
common /fcset804/c8040,c8041,c8042,c8043,c8044,c8045,c8046
common /fcset048/c0480,c0481,c0482,c0483,c0484,c0485,c0486
common /fcset480/c4800,c4801,c4802,c4803,c4804,c4805,c4806
common /fcset822/c8220,c8221,c8222,c8223,c8224,c8225,c8226
common /fcset282/c2820,c2821,c2822,c2823,c2824,c2825,c2826
common /fcset228/c2280,c2281,c2282,c2283,c2284,c2285,c2286
common /fcset066/c0660,c0661,c0662,c0663,c0664,c0665,c0666
common /fcset606/c6060,c6061,c6062,c6063,c6064,c6065,c6066
common /fcset660/c6600,c6601,c6602,c6603,c6604,c6605,c6606
common /fcset642/c6420,c6421,c6422,c6423,c6424,c6425,c6426
common /fcset264/c2640,c2641,c2642,c2643,c2644,c2645,c2646
common /fcset426/c4260,c4261,c4262,c4263,c4264,c4265,c4266
common /fcset624/c6240,c6241,c6242,c6243,c6244,c6245,c6246
common /fcset462/c4620,c4621,c4622,c4623,c4624,c4625,c4626
common /fcset246/c2460,c2461,c2462,c2463,c2464,c2465,c2466
common /fcset444/c4440,c4441,c4442,c4443,c4444,c4445,c4446
common /fcsetE00/cE000,cE001,cE002,cE003,cE004,cE005,cE006,cE007
common /fcset0E0/c0E00,c0E01,c0E02,c0E03,c0E04,c0E05,c0E06,c0E07
common /fcset00E/c00E0,c00E1,c00E2,c00E3,c00E4,c00E5,c00E6,c00E7
common /fcsetC20/cC200,cC201,cC202,cC203,cC204,cC205,cC206,cC207
common /fcset0C2/c0C20,c0C21,c0C22,c0C23,c0C24,c0C25,c0C26,c0C27
common /fcset20C/c20C0,c20C1,c20C2,c20C3,c20C4,c20C5,c20C6,c20C7
common /fcsetC02/cC020,cC021,cC022,cC023,cC024,cC025,cC026,cC027
common /fcset2C0/c2C00,c2C01,c2C02,c2C03,c2C04,c2C05,c2C06,c2C07
common /fcset02C/c02C0,c02C1,c02C2,c02C3,c02C4,c02C5,c02C6,c02C7
common /fcsetA40/cA400,cA401,cA402,cA403,cA404,cA405,cA406,cA407
common /fcset0A4/c0A40,c0A41,c0A42,c0A43,c0A44,c0A45,c0A46,c0A47
common /fcset40A/c40A0,c40A1,c40A2,c40A3,c40A4,c40A5,c40A6,c40A7
common /fcsetA04/cA040,cA041,cA042,cA043,cA044,cA045,cA046,cA047
common /fcset4A0/c4A00,c4A01,c4A02,c4A03,c4A04,c4A05,c4A06,c4A07
common /fcset04A/c04A0,c04A1,c04A2,c04A3,c04A4,c04A5,c04A6,c04A7
common /fcsetA22/cA220,cA221,cA222,cA223,cA224,cA225,cA226,cA227
common /fcset2A2/c2A20,c2A21,c2A22,c2A23,c2A24,c2A25,c2A26,c2A27
common /fcset22A/c22A0,c22A1,c22A2,c22A3,c22A4,c22A5,c22A6,c22A7
common /fcset860/c8600,c8601,c8602,c8603,c8604,c8605,c8606,c8607
common /fcset086/c0860,c0861,c0862,c0863,c0864,c0865,c0866,c0867
common /fcset608/c6080,c6081,c6082,c6083,c6084,c6085,c6086,c6087
common /fcset806/c8060,c8061,c8062,c8063,c8064,c8065,c8066,c8067
common /fcset680/c6800,c6801,c6802,c6803,c6804,c6805,c6806,c6807
common /fcset068/c0680,c0681,c0682,c0683,c0684,c0685,c0686,c0687
common /fcset842/c8420,c8421,c8422,c8423,c8424,c8425,c8426,c8427
common /fcset284/c2840,c2841,c2842,c2843,c2844,c2845,c2846,c2847
common /fcset428/c4280,c4281,c4282,c4283,c4284,c4285,c4286,c4287
common /fcset824/c8240,c8241,c8242,c8243,c8244,c8245,c8246,c8247
common /fcset482/c4820,c4821,c4822,c4823,c4824,c4825,c4826,c4827
common /fcset248/c2480,c2481,c2482,c2483,c2484,c2485,c2486,c2487
common /fcset266/c2660,c2661,c2662,c2663,c2664,c2665,c2666,c2667
common /fcset626/c6260,c6261,c6262,c6263,c6264,c6265,c6266,c6267
common /fcset662/c6620,c6621,c6622,c6623,c6624,c6625,c6626,c6627
common /fcset644/c6440,c6441,c6442,c6443,c6444,c6445,c6446,c6447
common /fcset464/c4640,c4641,c4642,c4643,c4644,c4645,c4646,c4647
common /fcset446/c4460,c4461,c4462,c4463,c4464,c4465,c4466,c4467
common /fcsetG00/cG000,cG001,cG002,cG003,cG004,cG005,cG006,cG007,cG008
common /fcset0G0/c0G00,c0G01,c0G02,c0G03,c0G04,c0G05,c0G06,c0G07,c0G08
common /fcset00G/c00G0,c00G1,c00G2,c00G3,c00G4,c00G5,c00G6,c00G7,c00G8
common /fcsetE20/cE200,cE201,cE202,cE203,cE204,cE205,cE206,cE207,cE208
common /fcset0E2/c0E20,c0E21,c0E22,c0E23,c0E24,c0E25,c0E26,c0E27,c0E28
common /fcset20E/c20E0,c20E1,c20E2,c20E3,c20E4,c20E5,c20E6,c20E7,c20E8
common /fcsetE02/cE020,cE021,cE022,cE023,cE024,cE025,cE026,cE027,cE028
common /fcset2E0/c2E00,c2E01,c2E02,c2E03,c2E04,c2E05,c2E06,c2E07,c2E08
common /fcset02E/c02E0,c02E1,c02E2,c02E3,c02E4,c02E5,c02E6,c02E7,c02E8
common /fcsetC40/cC400,cC401,cC402,cC403,cC404,cC405,cC406,cC407,cC408
common /fcset0C4/c0C40,c0C41,c0C42,c0C43,c0C44,c0C45,c0C46,c0C47,c0C48
common /fcset40C/c40C0,c40C1,c40C2,c40C3,c40C4,c40C5,c40C6,c40C7,c40C8
common /fcsetC04/cC040,cC041,cC042,cC043,cC044,cC045,cC046,cC047,cC048
common /fcset4C0/c4C00,c4C01,c4C02,c4C03,c4C04,c4C05,c4C06,c4C07,c4C08
common /fcset04C/c04C0,c04C1,c04C2,c04C3,c04C4,c04C5,c04C6,c04C7,c04C8
common /fcsetC22/cC220,cC221,cC222,cC223,cC224,cC225,cC226,cC227,cC228
common /fcset2C2/c2C20,c2C21,c2C22,c2C23,c2C24,c2C25,c2C26,c2C27,c2C28
common /fcset22C/c22C0,c22C1,c22C2,c22C3,c22C4,c22C5,c22C6,c22C7,c22C8
common /fcsetA60/cA600,cA601,cA602,cA603,cA604,cA605,cA606,cA607,cA608
common /fcset0A6/c0A60,c0A61,c0A62,c0A63,c0A64,c0A65,c0A66,c0A67,c0A68
common /fcset60A/c60A0,c60A1,c60A2,c60A3,c60A4,c60A5,c60A6,c60A7,c60A8
common /fcsetA06/cA060,cA061,cA062,cA063,cA064,cA065,cA066,cA067,cA068
common /fcset6A0/c6A00,c6A01,c6A02,c6A03,c6A04,c6A05,c6A06,c6A07,c6A08
common /fcset06A/c06A0,c06A1,c06A2,c06A3,c06A4,c06A5,c06A6,c06A7,c06A8
common /fcsetA42/cA420,cA421,cA422,cA423,cA424,cA425,cA426,cA427,cA428
common /fcset2A4/c2A40,c2A41,c2A42,c2A43,c2A44,c2A45,c2A46,c2A47,c2A48
common /fcset42A/c42A0,c42A1,c42A2,c42A3,c42A4,c42A5,c42A6,c42A7,c42A8
common /fcsetA24/cA240,cA241,cA242,cA243,cA244,cA245,cA246,cA247,cA248
common /fcset4A2/c4A20,c4A21,c4A22,c4A23,c4A24,c4A25,c4A26,c4A27,c4A28
common /fcset24A/c24A0,c24A1,c24A2,c24A3,c24A4,c24A5,c24A6,c24A7,c24A8
common /fcset088/c0880,c0881,c0882,c0883,c0884,c0885,c0886,c0887,c0888
common /fcset808/c8080,c8081,c8082,c8083,c8084,c8085,c8086,c8087,c8088
common /fcset880/c8800,c8801,c8802,c8803,c8804,c8805,c8806,c8807,c8808
common /fcset862/c8620,c8621,c8622,c8623,c8624,c8625,c8626,c8627,c8628
common /fcset286/c2860,c2861,c2862,c2863,c2864,c2865,c2866,c2867,c2868
common /fcset628/c6280,c6281,c6282,c6283,c6284,c6285,c6286,c6287,c6288
common /fcset826/c8260,c8261,c8262,c8263,c8264,c8265,c8266,c8267,c8268
common /fcset682/c6820,c6821,c6822,c6823,c6824,c6825,c6826,c6827,c6828
common /fcset268/c2680,c2681,c2682,c2683,c2684,c2685,c2686,c2687,c2688
common /fcset844/c8440,c8441,c8442,c8443,c8444,c8445,c8446,c8447,c8448
common /fcset484/c4840,c4841,c4842,c4843,c4844,c4845,c4846,c4847,c4848
common /fcset448/c4480,c4481,c4482,c4483,c4484,c4485,c4486,c4487,c4488
common /fcset466/c4660,c4661,c4662,c4663,c4664,c4665,c4666,c4667,c4668
common /fcset646/c6460,c6461,c6462,c6463,c6464,c6465,c6466,c6467,c6468
common /fcset664/c6640,c6641,c6642,c6643,c6644,c6645,c6646,c6647,c6648
!
if(n.ge.0) then
    ab=a*b
!
    gm1=8.d0*ab
endif
!
if(n.ge.2) then
    a2=a*a;b2=b*b;a2pb2=a2+b2
!
    c2000=(2.d0*a2-b2)*F6
    c2001=-F6
    c0200=(2.d0*b2-a2)*F6
    c0201=-F6
    c0020=-a2pb2*F6
    c0021=F3
endif
!
if(n.ge.4) then
    a4=a2*a2;b4=b2*b2;a2b2=a2*b2;a4pb4=a4+b4
!
    c4000=(24.d0*a4-40.d0*a2b2+9.d0*b4)*F120
    c4001=-(4.d0*a2-b2)*F12
    c4002=F3d40
    c0400=(24.d0*b4-40.d0*a2b2+9.d0*a4)*F120
    c0401=-(4.d0*b2-a2)*F12
    c0402=F3d40
    c0040=(9.d0*a4pb4+10.d0*a2b2)*F120
    c0041=-a2pb2*F3
    c0042=F5
    c0220=(3.d0*a4-5.d0*a2b2-12.d0*b4)*F20
    c0221=-(a2-9.d0*b2)*F4
    c0222=-F3d5
    c2020=(3.d0*b4-5.d0*a2b2-12.d0*a4)*F20
    c2021=-(b2-9.d0*a2)*F4
    c2022=-F3d5
    c2200=-(4.d0*a4pb4-15.d0*a2b2)*F3d20
    c2201=-a2pb2*F4
    c2202=F3d20
endif
!
if(n.ge.6) then
    a6=a4*a2;b6=b4*b2;a4b2=a4*b2;a2b4=a2*b4;a6pb6=a6+b6
!
    c6000=(16.d0*a6-56.d0*a4b2+42.d0*a2b4-5.d0*b6)*F112
    c6001=-(24.d0*a4-20.d0*a2b2+3.d0*b4)*F48
    c6002=(6.d0*a2-b2)*F16
    c6003=-F5d112
    c0600=(16.d0*b6-56.d0*a2b4+42.d0*a4b2-5.d0*a6)*F112
    c0601=-(24.d0*b4-20.d0*a2b2+3.d0*a4)*F48
    c0602=(6.d0*b2-a2)*F16
    c0603=-F5d112
    c0060=-(5.d0*a4pb4+2.d0*a2b2)*a2pb2*F112
    c0061=(9.d0*a4pb4+10.d0*a2b2)*F24
    c0062=-a2pb2*F2
    c0063=F7
    c4200=-(120.d0*a6-812.d0*a4b2+707.d0*a2b4-90.d0*b6)*F112
    c4201=(4.d0*a4-50.d0*a2b2+11.d0*b4)*F16
    c4202=(11.d0*a2+4.d0*b2)*F16
    c4203=-F15d112
    c0420=-(15.d0*a6-77.d0*a4b2-28.d0*a2b4+120.d0*b6)*F112
    c0421=(2.d0*a4-25.d0*a2b2+58.d0*b4)*F8
    c0422=(11.d0*a2-101.d0*b2)*F16
    c0423=F45d56
    c2040=(90.d0*a6+77.d0*a4b2+28.d0*a2b4-15.d0*b6)*F112
    c2041=-(101.d0*a4+50.d0*a2b2-11.d0*b4)*F16
    c2042=(29.d0*a2+b2)*F4
    c2043=-F15d14
    c2400=-(120.d0*b6-812.d0*a2b4+707.d0*a4b2-90.d0*a6)*F112
    c2401=(4.d0*b4-50.d0*a2b2+11.d0*a4)*F16
    c2402=(11.d0*b2+4.d0*a2)*F16
    c2403=-F15d112
    c4020=-(15.d0*b6-77.d0*a2b4-28.d0*a4b2+120.d0*a6)*F112
    c4021=(2.d0*b4-25.d0*a2b2+58.d0*a4)*F8
    c4022=(11.d0*b2-101.d0*a2)*F16
    c4023=F45d56
    c0240=(90.d0*b6+77.d0*a2b4+28.d0*a4b2-15.d0*a6)*F112
    c0241=-(101.d0*b4+50.d0*a2b2-11.d0*a4)*F16
    c0242=(29.d0*b2+a2)*F4
    c0243=-F15d14
    c2220=(2.d0*a6pb6-7.d0*a2pb2*a2b2)*F45d56
    c2221=-(3.d0*a4pb4-20.d0*a2b2)*F15d8
    c2222=-a2pb2*F45d8
    c2223=F45d28
endif
!
if(n.ge.8) then
    a8=a6*a2;b8=b6*b2;a6b2=a6*b2;a4b4=a4*b4;a2b6=a2*b6;a8pb8=a8+b8
!
    c8000=(640.d0*a8-3840.d0*a6b2+6048.d0*a4b4-2400.d0*a2b6+175.d0*b8)*F5760
    c8001=-(64.d0*a6-112.d0*a4b2+56.d0*a2b4-5.d0*b6)*F96
    c8002=(144.d0*a4-80.d0*a2b2+9.d0*b4)*F7d960
    c8003=-(8.d0*a2-b2)*F5d96
    c8004=F35d1152
    c0800=(640.d0*b8-3840.d0*a2b6+6048.d0*a4b4-2400.d0*a6b2+175.d0*a8)*F5760
    c0801=-(64.d0*b6-112.d0*a2b4+56.d0*a4b2-5.d0*a6)*F96
    c0802=(144.d0*b4-80.d0*a2b2+9.d0*a4)*F7d960
    c0803=-(8.d0*b2-a2)*F5d96
    c0804=F35d1152
    c0080=(175.d0*a8pb8+300.d0*a4pb4*a2b2+378.d0*a4b4)*F5760
    c0081=-(5.d0*a4pb4+2.d0*a2b2)*a2pb2*F12
    c0082=(9.d0*a4pb4+10.d0*a2b2)*F7d60
    c0083=-a2pb2*F2d3
    c0084=F9
    c6200=-(2240.d0*a8-24240.d0*a6b2+43848.d0*a4b4-18525.d0*a2b6 &
        +1400.d0*b8)*F1440
    c6201=(176.d0*a6-1568.d0*a4b2+1099.d0*a2b4-115.d0*b6)*F96
    c6202=(72.d0*a4+335.d0*a2b2-63.d0*b4)*F7d480
    c6203=-(23.d0*a2+5.d0*b2)*F5d96
    c6204=F35d288
    c0620=(175.d0*a8-1725.d0*a6b2+1512.d0*a4b4+2640.d0*a2b6 &
        -2240.d0*b8)*F1440
    c0621=-(25.d0*a6-469.d0*a4b2+1568.d0*a2b4-1616.d0*b6)*F96
    c0622=-(63.d0*a4-785.d0*a2b2+2088.d0*b4)*F7d480
    c0623=-(23.d0*a2-247.d0*b2)*F5d96
    c0624=-F35d36
    c2060=-(1400.d0*a8+1725.d0*a6b2+1323.d0*a4b4+375.d0*a2b6 &
        -175.d0*b8)*F1440
    c2061=(1235.d0*a6+1099.d0*a4b2+469.d0*a2b4-115.d0*b6)*F96
    c2062=-(261.d0*a4+140.d0*a2b2-9.d0*b4)*F7d60
    c2063=(101.d0*a2+11.d0*b2)*F6
    c2064=-F14d9
    c2600=-(2240.d0*b8-24240.d0*a2b6+43848.d0*a4b4-18525.d0*a6b2 &
        +1400.d0*a8)*F1440
    c2601=(176.d0*b6-1568.d0*a2b4+1099.d0*a4b2-115.d0*a6)*F96
    c2602=(72.d0*b4+335.d0*a2b2-63.d0*a4)*F7d480
    c2603=-(23.d0*b2+5.d0*a2)*F5d96
    c2604=F35d288
    c6020=(175.d0*b8-1725.d0*a2b6+1512.d0*a4b4+2640.d0*a6b2 &
        -2240.d0*a8)*F1440
    c6021=-(25.d0*b6-469.d0*a2b4+1568.d0*a4b2-1616.d0*a6)*F96
    c6022=-(63.d0*b4-785.d0*a2b2+2088.d0*a4)*F7d480
    c6023=-(23.d0*b2-247.d0*a2)*F5d96
    c6024=-F35d36
    c0260=-(1400.d0*b8+1725.d0*a2b6+1323.d0*a4b4+375.d0*a6b2 &
        -175.d0*a8)*F1440
    c0261=(1235.d0*b6+1099.d0*a2b4+469.d0*a4b2-115.d0*a6)*F96
    c0262=-(261.d0*b4+140.d0*a2b2-9.d0*a4)*F7d60
    c0263=(101.d0*b2+11.d0*a2)*F6
    c0264=-F14d9
    c0440=(35.d0*a8-210.d0*a6b2-189.d0*a4b4+240.d0*a2b6+560.d0*b8)*F192
    c0441=-(7.d0*a6-70.d0*a4b2+35.d0*a2b4+232.d0*b6)*F5d32
    c0442=-(9.d0*a4+50.d0*a2b2-711.d0*b4)*F7d64
    c0443=(a2-29.d0*b2)*F5d4
    c0444=F35d12
    c4040=(35.d0*b8-210.d0*a2b6-189.d0*a4b4+240.d0*a6b2+560.d0*a8)*F192
    c4041=-(7.d0*b6-70.d0*a2b4+35.d0*a4b2+232.d0*a6)*F5d32
    c4042=-(9.d0*b4+50.d0*a2b2-711.d0*a4)*F7d64
    c4043=(b2-29.d0*a2)*F5d4
    c4044=F35d12
    c4400=(560.d0*a8pb8-6960.d0*a4pb4*a2b2+14931.d0*a4b4)*F192
    c4401=(8.d0*a4pb4-43.d0*a2b2)*a2pb2*F5d32
    c4402=-(9.d0*a4pb4-100.d0*a2b2)*F7d64
    c4403=-a2pb2*F35d32
    c4404=F35d192
    c4220=(112.d0*a8-672.d0*a6b2-189.d0*a4b4+471.d0*a2b6-56.d0*b8)*F5d96
    c4221=-(224.d0*a6-1778.d0*a4b2+889.d0*a2b4-67.d0*b6)*F5d32
    c4222=-(9.d0*a4+127.d0*a2b2-18.d0*b4)*F35d32
    c4223=(157.d0*a2+67.d0*b2)*F5d32
    c4224=-F35d12
    c2420=(112.d0*b8-672.d0*a2b6-189.d0*a4b4+471.d0*a6b2-56.d0*a8)*F5d96
    c2421=-(224.d0*b6-1778.d0*a2b4+889.d0*a4b2-67.d0*a6)*F5d32
    c2422=-(9.d0*b4+127.d0*a2b2-18.d0*a4)*F35d32
    c2423=(157.d0*b2+67.d0*a2)*F5d32
    c2424=-F35d12
    c2240=-(56.d0*a8pb8-201.d0*a4pb4*a2b2-378.d0*a4b4)*F5d96
    c2241=(157.d0*a4pb4-1046.d0*a2b2)*a2pb2*F5d32
    c2242=-(9.d0*a4pb4-254.d0*a2b2)*F35d32
    c2243=-a2pb2*35.d0
    c2244=F35d6
endif
!
if(n.ge.10) then
    aA=a8*a2;bA=b8*b2;a8b2=a8*b2;a6b4=a6*b4;a4b6=a4*b6;a2b8=a2*b8;aApbA=aA+bA
!
    cA000=(768.d0*aA-7040.d0*a8b2+19008.d0*a6b4-15840.d0*a4b6+3850.d0*a2b8 &
        -189.d0*bA)*F8448
    cA001=-(7040.d0*a8-21120.d0*a6b2+22176.d0*a4b4-6600.d0*a2b6+385.d0*b8) &
        *F8448
    cA002=(96.d0*a6-112.d0*a4b2+42.d0*a2b4-3.d0*b6)*F3d128
    cA003=-(240.d0*a4-100.d0*a2b2+9.d0*b4)*F128
    cA004=(10.d0*a2-b2)*F35d768
    cA005=-F63d2816
    c0A00=(768.d0*bA-7040.d0*a2b8+19008.d0*a4b6-15840.d0*a6b4+3850.d0*a8b2 &
        -189.d0*aA)*F8448
    c0A01=-(7040.d0*b8-21120.d0*a2b6+22176.d0*a4b4-6600.d0*a6b2+385.d0*a8) &
        *F8448
    c0A02=(96.d0*b6-112.d0*a2b4+42.d0*a4b2-3.d0*a6)*F3d128
    c0A03=-(240.d0*b4-100.d0*a2b2+9.d0*a4)*F128
    c0A04=(10.d0*b2-a2)*F35d768
    c0A05=-F63d2816
    c00A0=-(189.d0*aApbA+385.d0*a6pb6*a2b2+594.d0*a2pb2*a4b4)*F8448
    c00A1=(175.d0*a8pb8+300.d0*a4pb4*a2b2+378.d0*a4b4)*F384
    c00A2=-(5.d0*a4pb4+2.d0*a2b2)*a2pb2*F3d8
    c00A3=(9.d0*a4pb4+10.d0*a2b2)*F4
    c00A4=-a2pb2*F5d6
    c00A5=F11
    c8200=-(5760.d0*aA-91520.d0*a8b2+281952.d0*a6b4-249480.d0*a4b6 &
        +62755.d0*a2b8-3150.d0*bA)*F2816
    c8201=(1280.d0*a8-14400.d0*a6b2+20664.d0*a4b4-7140.d0*a2b6+455.d0*b8) &
        *F256
    c8202=-(16.d0*a6-532.d0*a4b2+315.d0*a2b4-28.d0*b6)*F9d128
    c8203=-(180.d0*a4+310.d0*a2b2-51.d0*b4)*F3d128
    c8204=(13.d0*a2+2.d0*b2)*F35d256
    c8205=-F315d2816
    c0820=-(315.d0*aA-5005.d0*a8b2+11880.d0*a6b4+3168.d0*a4b6-14080.d0*a2b8 &
        +5760.d0*bA)*F2816
    c0821=(35.d0*a8-930.d0*a6b2+4788.d0*a4b4-7200.d0*a2b6+4160.d0*b8)*F128
    c0822=(17.d0*a6-315.d0*a4b2+1148.d0*a2b4-1424.d0*b6)*F9d128
    c0823=(6.d0*a4-85.d0*a2b2+270.d0*b4)*F21d64
    c0824=(13.d0*a2-163.d0*b2)*F35d256
    c0825=F1575d1408
    c2080=(3150.d0*aA+5005.d0*a8b2+5544.d0*a6b4+3366.d0*a4b6+770.d0*a2b8 &
        -315.d0*bA)*F2816
    c2081=-(1141.d0*a8+1428.d0*a6b2+1134.d0*a4b4+372.d0*a2b6-91.d0*b8) &
        *F5d256
    c2082=(315.d0*a6+287.d0*a4b2+133.d0*a2b4-15.d0*b6)*F9d32
    c2083=-(89.d0*a4+50.d0*a2b2+b4)*F9d8
    c2084=(65.d0*a2+10.d0*b2)*F2
    c2085=-F45d22
    c2800=-(5760.d0*bA-91520.d0*a2b8+281952.d0*a4b6-249480.d0*a6b4 &
        +62755.d0*a8b2-3150.d0*aA)*F2816
    c2801=(1280.d0*b8-14400.d0*a2b6+20664.d0*a4b4-7140.d0*a6b2+455.d0*a8) &
        *F256
    c2802=-(16.d0*b6-532.d0*a2b4+315.d0*a4b2-28.d0*a6)*F9d128
    c2803=-(180.d0*b4+310.d0*a2b2-51.d0*a4)*F3d128
    c2804=(13.d0*b2+2.d0*a2)*F35d256
    c2805=-F315d2816
    c8020=-(315.d0*bA-5005.d0*a2b8+11880.d0*a4b6+3168.d0*a6b4-14080.d0*a8b2 &
        +5760.d0*aA)*F2816
    c8021=(35.d0*b8-930.d0*a2b6+4788.d0*a4b4-7200.d0*a6b2+4160.d0*a8)*F128
    c8022=(17.d0*b6-315.d0*a2b4+1148.d0*a4b2-1424.d0*a6)*F9d128
    c8023=(6.d0*b4-85.d0*a2b2+270.d0*a4)*F21d64
    c8024=(13.d0*b2-163.d0*a2)*F35d256
    c8025=F1575d1408
    c0280=(3150.d0*bA+5005.d0*a2b8+5544.d0*a4b6+3366.d0*a6b4+770.d0*a8b2 &
        -315.d0*aA)*F2816
    c0281=-(1141.d0*b8+1428.d0*a2b6+1134.d0*a4b4+372.d0*a6b2-91.d0*a8) &
        *F5d256
    c0282=(315.d0*b6+287.d0*a2b4+133.d0*a4b2-15.d0*a6)*F9d32
    c0283=-(89.d0*b4+50.d0*a2b2+a4)*F9d8
    c0284=(65.d0*b2+10.d0*a2)*F2
    c0285=-F45d22
    c6400=(4320.d0*aA-78320.d0*a8b2+274626.d0*a6b4-262251.d0*a4b6 &
        +69300.d0*a2b8-3600.d0*bA)*F7d4224
    c6401=-(16.d0*a8-444.d0*a6b2+1575.d0*a4b4-774.d0*a2b6+60.d0*b8)*F35d384
    c6402=-(42.d0*a6-357.d0*a4b2+28.d0*a2b4+9.d0*b6)*F21d128
    c6403=-(27.d0*a4+470.d0*a2b2-39.d0*b4)*F7d128
    c6404=(28.d0*a2+17.d0*b2)*F35d384
    c6405=-F315d1408
    c0640=-(135.d0*aA-1540.d0*a8b2+891.d0*a6b4+4158.d0*a4b6+880.d0*a2b8 &
        -4320.d0*bA)*F7d4224
    c0641=(85.d0*a8-1410.d0*a6b2+3213.d0*a4b4+2220.d0*a2b6-7120.d0*b8)*F7d384
    c0642=(13.d0*a6-28.d0*a4b2-875.d0*a2b4+2774.d0*b6)*F21d128
    c0643=-(9.d0*a4-430.d0*a2b2+2649.d0*b4)*F21d128
    c0644=-(a2-21.d0*b2)*F175d32
    c0645=-F525d88
    c4060=-(3600.d0*aA+3300.d0*a8b2+891.d0*a6b4-1287.d0*a4b6-935.d0*a2b8 &
        +135.d0*bA)*F7d4224
    c4061=(3150.d0*a8+1935.d0*a6b2-126.d0*a4b4-705.d0*a2b6+70.d0*b8)*F7d192
    c4062=-(2649.d0*a6+875.d0*a4b2-357.d0*a2b4+9.d0*b6)*F21d128
    c4063=(4161.d0*a4+370.d0*a2b2-63.d0*b4)*F7d64
    c4064=-(89.d0*a2+b2)*F35d24
    c4065=F315d44
    c4600=(4320.d0*bA-78320.d0*a2b8+274626.d0*a4b6-262251.d0*a6b4 &
        +69300.d0*a8b2-3600.d0*aA)*F7d4224
    c4601=-(16.d0*b8-444.d0*a2b6+1575.d0*a4b4-774.d0*a6b2+60.d0*a8)*F35d384
    c4602=-(42.d0*b6-357.d0*a2b4+28.d0*a4b2+9.d0*a6)*F21d128
    c4603=-(27.d0*b4+470.d0*a2b2-39.d0*a4)*F7d128
    c4604=(28.d0*b2+17.d0*a2)*F35d384
    c4605=-F315d1408
    c6040=-(135.d0*bA-1540.d0*a2b8+891.d0*a4b6+4158.d0*a6b4+880.d0*a8b2 &
        -4320.d0*aA)*F7d4224
    c6041=(85.d0*b8-1410.d0*a2b6+3213.d0*a4b4+2220.d0*a6b2-7120.d0*a8)*F7d384
    c6042=(13.d0*b6-28.d0*a2b4-875.d0*a4b2+2774.d0*a6)*F21d128
    c6043=-(9.d0*b4-430.d0*a2b2+2649.d0*a4)*F21d128
    c6044=-(b2-21.d0*a2)*F175d32
    c6045=-F525d88
    c0460=-(3600.d0*bA+3300.d0*a2b8+891.d0*a4b6-1287.d0*a6b4-935.d0*a8b2 &
        +135.d0*aA)*F7d4224
    c0461=(3150.d0*b8+1935.d0*a2b6-126.d0*a4b4-705.d0*a6b2+70.d0*a8)*F7d192
    c0462=-(2649.d0*b6+875.d0*a2b4-357.d0*a4b2+9.d0*a6)*F21d128
    c0463=(4161.d0*b4+370.d0*a2b2-63.d0*a4)*F7d64
    c0464=-(89.d0*b2+a2)*F35d24
    c0465=F315d44
    c6220=(1440.d0*aA-13200.d0*a8b2+7326.d0*a6b4+12771.d0*a4b6 &
        -6545.d0*a2b8+450.d0*bA)*F7d704
    c6221=-(1200.d0*a8-12180.d0*a6b2+12789.d0*a4b4-3270.d0*a2b6 &
        +155.d0*b8)*F7d64
    c6222=(74.d0*a6-1421.d0*a4b2+658.d0*a2b4-47.d0*b6)*F63d64
    c6223=(387.d0*a4+1090.d0*a2b2-141.d0*b4)*F21d64
    c6224=-(119.d0*a2+31.d0*b2)*F35d64
    c6225=F1575d352
    c2620=(1440.d0*bA-13200.d0*a2b8+7326.d0*a4b6+12771.d0*a6b4 &
        -6545.d0*a8b2+450.d0*aA)*F7d704
    c2621=-(1200.d0*b8-12180.d0*a2b6+12789.d0*a4b4-3270.d0*a6b2 &
        +155.d0*a8)*F7d64
    c2622=(74.d0*b6-1421.d0*a2b4+658.d0*a4b2-47.d0*a6)*F63d64
    c2623=(387.d0*b4+1090.d0*a2b2-141.d0*a4)*F21d64
    c2624=-(119.d0*b2+31.d0*a2)*F35d64
    c2625=F1575d352
    c2260=(450.d0*aApbA-1705.d0*a6pb6*a2b2-4653.d0*a2pb2*a4b4)*F7d704
    c2261=-(595.d0*a8pb8-3270.d0*a4pb4*a2b2-5922.d0*a4b4)*F7d64
    c2262=(129.d0*a6pb6-1421.d0*a2pb2*a2b2)*F63d64
    c2263=(111.d0*a4pb4+2030.d0*a2b2)*F21d32
    c2264=-a2pb2*F525d4
    c2265=F315d22
    c2440=(270.d0*aA-2475.d0*a8b2-396.d0*a6b4+5049.d0*a4b6+4180.d0*a2b8 &
        -720.d0*bA)*F35d1408
    c2441=-(225.d0*a8-2820.d0*a6b2+2961.d0*a4b4+6090.d0*a2b6-820.d0*b8) &
        *F35d128
    c2442=-(4.d0*a6+329.d0*a4b2-1750.d0*a2b4+125.d0*b6)*F315d128
    c2443=(153.d0*a4-2030.d0*a2b2-375.d0*b4)*F105d128
    c2444=(19.d0*a2+41.d0*b2)*F175d32
    c2445=-F1575d88
    c4240=(270.d0*bA-2475.d0*a2b8-396.d0*a4b6+5049.d0*a6b4+4180.d0*a8b2 &
        -720.d0*aA)*F35d1408
    c4241=-(225.d0*b8-2820.d0*a2b6+2961.d0*a4b4+6090.d0*a6b2-820.d0*a8) &
        *F35d128
    c4242=-(4.d0*b6+329.d0*a2b4-1750.d0*a4b2+125.d0*a6)*F315d128
    c4243=(153.d0*b4-2030.d0*a2b2-375.d0*a4)*F105d128
    c4244=(19.d0*b2+41.d0*a2)*F175d32
    c4245=-F1575d88
    c4420=-(144.d0*a8pb8-1948.d0*a4pb4*a2b2+4423.d0*a4b4)*a2pb2*F175d1408
    c4421=(38.d0*a8pb8-609.d0*a4pb4*a2b2+1575.d0*a4b4)*F175d64
    c4422=(51.d0*a4pb4-380.d0*a2b2)*a2pb2*F315d128
    c4423=-(3.d0*a4pb4-235.d0*a2b2)*F105d32
    c4424=-a2pb2*F7875d128
    c4425=F4725d704
endif
!
if(n.ge.12) then
    aC=aA*a2;bC=bA*b2;aAb2=aA*b2;a8b4=a8*b4;a6b6=a6*b6;a4b8=a4*b8;a2bA=a2*bA
    aCpbC=aC+bC
!
    cC000=(7168.d0*aC-93184.d0*aAb2+384384.d0*a8b4-549120.d0*a6b6 &
        +280280.d0*a4b8-45864.d0*a2bA+1617.d0*bC)*F93184
    cC001=-(1536.d0*aA-7040.d0*a8b2+12672.d0*a6b4-7920.d0*a4b6+1540.d0*a2b8 &
        -63.d0*bA)*F1536
    cC002=(1920.d0*a8-3840.d0*a6b2+3024.d0*a4b4-720.d0*a2b6+35.d0*b8) &
        *F11d5120
    cC003=-(320.d0*a6-280.d0*a4b2+84.d0*a2b4-5.d0*b6)*F33d1792
    cC004=(120.d0*a4-40.d0*a2b2+3.d0*b4)*F77d3072
    cC005=-(12.d0*a2-b2)*F21d512
    cC006=F231d13312
    c0C00=(7168.d0*bC-93184.d0*a2bA+384384.d0*a4b8-549120.d0*a6b6 &
        +280280.d0*a8b4-45864.d0*aAb2+1617.d0*aC)*F93184
    c0C01=-(1536.d0*bA-7040.d0*a2b8+12672.d0*a4b6-7920.d0*a6b4+1540.d0*a8b2 &
        -63.d0*aA)*F1536
    c0C02=(1920.d0*b8-3840.d0*a2b6+3024.d0*a4b4-720.d0*a6b2+35.d0*a8) &
        *F11d5120
    c0C03=-(320.d0*b6-280.d0*a2b4+84.d0*a4b2-5.d0*a6)*F33d1792
    c0C04=(120.d0*b4-40.d0*a2b2+3.d0*a4)*F77d3072
    c0C05=-(12.d0*b2-a2)*F21d512
    c0C06=F231d13312
    c00C0=(1617.d0*aCpbC+3822.d0*a8pb8*a2b2+7007.d0*a4pb4*a4b4 &
        +8580.d0*a6b6)*F93184
    c00C1=-(189.d0*a8pb8+196.d0*a4pb4*a2b2+398.d0*a4b4)*a2pb2*F384
    c00C2=(175.d0*a8pb8+300.d0*a4pb4*a2b2+378.d0*a4b4)*F11d640
    c00C3=-(5.d0*a4pb4+2.d0*a2b2)*a2pb2*F33d28
    c00C4=(9.d0*a4pb4+10.d0*a2b2)*F11d24
    c00C5=-a2pb2
    c00C6=F13
    cA200=-(118272.d0*aC-2597504.d0*aAb2+12172160.d0*a8b4-18429840.d0*a6b6 &
        +9725716.d0*a4b8-1626261.d0*a2bA+58212.d0*bC)*F46592
    cA201=(15744.d0*aA-232320.d0*a8b2+562320.d0*a6b4-405504.d0*a4b6 &
        +85855.d0*a2b8-3717.d0*bA)*F1536
    cA202=-(2560.d0*a8-41520.d0*a6b2+49896.d0*a4b4-14610.d0*a2b6+805.d0*b8) &
        *F11d2560
    cA203=-(360.d0*a6+4144.d0*a4b2-2135.d0*a2b4+165.d0*b6)*F33d1792
    cA204=(204.d0*a4+205.d0*a2b2-30.d0*b4)*F77d1536
    cA205=-(59.d0*a2+7.d0*b2)*F21d512
    cA206=F693d6656
    c0A20=(4851.d0*aC-112749.d0*aAb2+476476.d0*a8b4-308880.d0*a6b6 &
        -512512.d0*a4b8+477568.d0*a2bA-118272.d0*bC)*F46592
    c0A21=-(441.d0*aA-15785.d0*a8b2+117216.d0*a6b4-274032.d0*a4b6 &
        +232320.d0*a2b8-85632.d0*bA)*F1536
    c0A22=-(175.d0*a8-4575.d0*a6b2+24948.d0*a4b4-42600.d0*a2b6+30400.d0*b8) &
        *F11d1280
    c0A23=-(165.d0*a6-3409.d0*a4b2+14336.d0*a2b4-21480.d0*b6)*F33d1792
    c0A24=-(69.d0*a4-1115.d0*a2b2+4164.d0*b4)*F77d1536
    c0A25=-(59.d0*a2-851.d0*b2)*F21d512
    c0A26=-F2079d1664
    c20A0=-(58212.d0*aC+112749.d0*aAb2+161161.d0*a8b4+141570.d0*a6b6 &
        +70070.d0*a4b8+13377.d0*a2bA-4851.d0*bC)*F46592
    c20A1=(53613.d0*aA+85855.d0*a8b2+96426.d0*a6b4+60390.d0*a4b6 &
        +15785.d0*a2b8-3717.d0*bA)*F1536
    c20A2=-(12145.d0*a8+15360.d0*a6b2+12474.d0*a4b4+4440.d0*a2b6-595.d0*b8) &
        *F11d640
    c20A3=(2685.d0*a6+2485.d0*a4b2+1211.d0*a2b4-45.d0*b6)*F33d224
    c20A4=-(95.d0*a4+55.d0*a2b2+4.d0*b4)*F11d4
    c20A5=(223.d0*a2+41.d0*b2)*F4
    c20A6=-F33d13
    c2A00=-(118272.d0*bC-2597504.d0*a2bA+12172160.d0*a4b8-18429840.d0*a6b6 &
        +9725716.d0*a8b4-1626261.d0*aAb2+58212.d0*aC)*F46592
    c2A01=(15744.d0*bA-232320.d0*a2b8+562320.d0*a4b6-405504.d0*a6b4 &
        +85855.d0*a8b2-3717.d0*aA)*F1536
    c2A02=-(2560.d0*b8-41520.d0*a2b6+49896.d0*a4b4-14610.d0*a6b2+805.d0*a8) &
        *F11d2560
    c2A03=-(360.d0*b6+4144.d0*a2b4-2135.d0*a4b2+165.d0*a6)*F33d1792
    c2A04=(204.d0*b4+205.d0*a2b2-30.d0*a4)*F77d1536
    c2A05=-(59.d0*b2+7.d0*a2)*F21d512
    c2A06=F693d6656
    cA020=(4851.d0*bC-112749.d0*a2bA+476476.d0*a4b8-308880.d0*a6b6 &
        -512512.d0*a8b4+477568.d0*aAb2-118272.d0*aC)*F46592
    cA021=-(441.d0*bA-15785.d0*a2b8+117216.d0*a4b6-274032.d0*a6b4 &
        +232320.d0*a8b2-85632.d0*aA)*F1536
    cA022=-(175.d0*b8-4575.d0*a2b6+24948.d0*a4b4-42600.d0*a6b2+30400.d0*a8) &
        *F11d1280
    cA023=-(165.d0*b6-3409.d0*a2b4+14336.d0*a4b2-21480.d0*a6)*F33d1792
    cA024=-(69.d0*b4-1115.d0*a2b2+4164.d0*a4)*F77d1536
    cA025=-(59.d0*b2-851.d0*a2)*F21d512
    cA026=-F2079d1664
    c02A0=-(58212.d0*bC+112749.d0*a2bA+161161.d0*a4b8+141570.d0*a6b6 &
        +70070.d0*a8b4+13377.d0*aAb2-4851.d0*aC)*F46592
    c02A1=(53613.d0*bA+85855.d0*a2b8+96426.d0*a4b6+60390.d0*a6b4 &
        +15785.d0*a8b2-3717.d0*aA)*F1536
    c02A2=-(12145.d0*b8+15360.d0*a2b6+12474.d0*a4b4+4440.d0*a6b2-595.d0*a8) &
        *F11d640
    c02A3=(2685.d0*b6+2485.d0*a2b4+1211.d0*a4b2-45.d0*a6)*F33d224
    c02A4=-(95.d0*b4+55.d0*a2b2+4.d0*a4)*F11d4
    c02A5=(223.d0*b2+41.d0*a2)*F4
    c02A6=-F33d13
    c8400=(266112.d0*aC-6639360.d0*aAb2+34674640.d0*a8b4-55989648.d0*a6b6 &
        +30809779.d0*a4b8-5304936.d0*a2bA+194040.d0*bC)*F5d93184
    c8401=-(4608.d0*aA-101200.d0*a8b2+398376.d0*a6b4-366102.d0*a4b6 &
        +89705.d0*a2b8-4284.d0*bA)*F5d1536
    c8402=-(1520.d0*a8-10320.d0*a6b2-20538.d0*a4b4+13080.d0*a2b6-1015.d0*b8) &
        *F11d1024
    c8403=(132.d0*a6-2345.d0*a4b2+385.d0*a2b4+15.d0*b6)*F165d1792
    c8404=(87.d0*a4+400.d0*a2b2-30.d0*b4)*F385d3072
    c8405=-(23.d0*a2+10.d0*b2)*F105d512
    c8406=F3465d13312
    c0840=(4851.d0*aC-87906.d0*aAb2+203203.d0*a8b4+226512.d0*a6b6 &
        -304304.d0*a4b8-279552.d0*a2bA+266112.d0*bC)*F5d93184
    c0841=-(315.d0*aA-7700.d0*a8b2+33165.d0*a6b4-17028.d0*a4b6-50600.d0*a2b8 &
        +54720.d0*bA)*F5d768
    c0842=-(175.d0*a8-1650.d0*a6b2-10269.d0*a4b4+60360.d0*a2b6-86600.d0*b8) &
        *F11d512
    c0843=(15.d0*a6-1526.d0*a4b2+12943.d0*a2b4-32628.d0*b6)*F165d1792
    c0844=(87.d0*a4-2330.d0*a2b2+13191.d0*b4)*F385d3072
    c0845=(17.d0*a2-347.d0*b2)*F105d128
    c0846=F17325d1664
    c4080=(194040.d0*aC+259896.d0*aAb2+203203.d0*a8b4+25740.d0*a6b6 &
        -70070.d0*a4b8-38220.d0*a2bA+4851.d0*bC)*F5d93184
    c4081=-(87444.d0*aA+89705.d0*a8b2+43164.d0*a6b4-10890.d0*a4b6 &
        -15400.d0*a2b8+1449.d0*bA)*F5d1536
    c4082=(153895.d0*a8+110940.d0*a6b2+20538.d0*a4b4-20100.d0*a2b6+1015.d0*b8) &
        *F11d1024
    c4083=-(8157.d0*a6+3521.d0*a4b2-301.d0*a2b4-33.d0*b6)*F165d448
    c4084=(6495.d0*a4+1150.d0*a2b2-57.d0*b4)*F55d192
    c4085=-(95.d0*a2+4.d0*b2)*F15d4
    c4086=F1485d104
    c4800=(266112.d0*bC-6639360.d0*a2bA+34674640.d0*a4b8-55989648.d0*a6b6 &
        +30809779.d0*a8b4-5304936.d0*aAb2+194040.d0*aC)*F5d93184
    c4801=-(4608.d0*bA-101200.d0*a2b8+398376.d0*a4b6-366102.d0*a6b4 &
        +89705.d0*a8b2-4284.d0*aA)*F5d1536
    c4802=-(1520.d0*b8-10320.d0*a2b6-20538.d0*a4b4+13080.d0*a6b2-1015.d0*a8) &
        *F11d1024
    c4803=(132.d0*b6-2345.d0*a2b4+385.d0*a4b2+15.d0*a6)*F165d1792
    c4804=(87.d0*b4+400.d0*a2b2-30.d0*a4)*F385d3072
    c4805=-(23.d0*b2+10.d0*a2)*F105d512
    c4806=F3465d13312
    c8040=(4851.d0*bC-87906.d0*a2bA+203203.d0*a4b8+226512.d0*a6b6 &
        -304304.d0*a8b4-279552.d0*aAb2+266112.d0*aC)*F5d93184
    c8041=-(315.d0*bA-7700.d0*a2b8+33165.d0*a4b6-17028.d0*a6b4-50600.d0*a8b2 &
        +54720.d0*aA)*F5d768
    c8042=-(175.d0*b8-1650.d0*a2b6-10269.d0*a4b4+60360.d0*a6b2-86600.d0*a8) &
        *F11d512
    c8043=(15.d0*b6-1526.d0*a2b4+12943.d0*a4b2-32628.d0*a6)*F165d1792
    c8044=(87.d0*b4-2330.d0*a2b2+13191.d0*a4)*F385d3072
    c8045=(17.d0*b2-347.d0*a2)*F105d128
    c8046=F17325d1664
    c0480=(194040.d0*bC+259896.d0*a2bA+203203.d0*a4b8+25740.d0*a6b6 &
        -70070.d0*a8b4-38220.d0*aAb2+4851.d0*aC)*F5d93184
    c0481=-(87444.d0*bA+89705.d0*a2b8+43164.d0*a4b6-10890.d0*a6b4 &
        -15400.d0*a8b2+1449.d0*aA)*F5d1536
    c0482=(153895.d0*b8+110940.d0*a2b6+20538.d0*a4b4-20100.d0*a6b2+1015.d0*a8) &
        *F11d1024
    c0483=-(8157.d0*b6+3521.d0*a2b4-301.d0*a4b2-33.d0*a6)*F165d448
    c0484=(6495.d0*b4+1150.d0*a2b2-57.d0*a4)*F55d192
    c0485=-(95.d0*b2+4.d0*a2)*F15d4
    c0486=F1485d104
    c8220=(88704.d0*aC-1153152.d0*aAb2+1841840.d0*a8b4+700128.d0*a6b6 &
        -1632631.d0*a4b8+426153.d0*a2bA-19404.d0*bC)*F15d46592
    c8221=-(38016.d0*aA-494560.d0*a8b2+890208.d0*a6b4-484308.d0*a4b6 &
        +78155.d0*a2b8-2583.d0*bA)*F5d512
    c8222=(920.d0*a8-13488.d0*a6b2+12915.d0*a4b4-3075.d0*a2b6+140.d0*b8) &
        *F165d256
    c8223=(816.d0*a6+17122.d0*a4b2-7175.d0*a2b4+465.d0*b6)*F495d1792
    c8224=-(699.d0*a4+1015.d0*a2b2-120.d0*b4)*F385d512
    c8225=(223.d0*a2+41.d0*b2)*F315d512
    c8226=-F10395d1664
    c2820=(88704.d0*bC-1153152.d0*a2bA+1841840.d0*a4b8+700128.d0*a6b6 &
        -1632631.d0*a8b4+426153.d0*aAb2-19404.d0*aC)*F15d46592
    c2821=-(38016.d0*bA-494560.d0*a2b8+890208.d0*a4b6-484308.d0*a6b4 &
        +78155.d0*a8b2-2583.d0*aA)*F5d512
    c2822=(920.d0*b8-13488.d0*a2b6+12915.d0*a4b4-3075.d0*a6b2+140.d0*a8) &
        *F165d256
    c2823=(816.d0*b6+17122.d0*a2b4-7175.d0*a4b2+465.d0*a6)*F495d1792
    c2824=-(699.d0*b4+1015.d0*a2b2-120.d0*a4)*F385d512
    c2825=(223.d0*b2+41.d0*a2)*F315d512
    c2826=-F10395d1664
    c2280=-(19404.d0*aCpbC-78351.d0*a8pb8*a2b2-280280.d0*a4pb4*a4b4 &
        -398970.d0*a6b6)*F15d46592
    c2281=(14049.d0*a8pb8-92204.d0*a4pb4*a2b2-110746.d0*a4b4)*a2pb2*F5d512
    c2282=-(1631.d0*a8pb8-14676.d0*a4pb4*a2b2-25830.d0*a4b4)*F165d512
    c2283=(51.d0*a4pb4-2018.d0*a2b2)*a2pb2*F495d112
    c2284=(69.d0*a4pb4+562.d0*a2b2)*F275d32
    c2285=-a2pb2*F1485d4
    c2286=F1485d52
    c0660=(1155.d0*aC-15015.d0*aAb2+5005.d0*a8b4+57915.d0*a6b6+44044.d0*a4b8 &
        -32760.d0*a2bA-73920.d0*bC)*F3328
    c0661=-(495.d0*aA-8525.d0*a8b2+15345.d0*a6b4+30789.d0*a4b6-14960.d0*a2b8 &
        -64440.d0*bA)*F7d768
    c0662=(25.d0*a8-2325.d0*a6b2+14679.d0*a4b4-3495.d0*a2b6-54380.d0*b8) &
        *F77d1280
    c0663=(135.d0*a6-2177.d0*a4b2-1631.d0*a2b4+43737.d0*b6)*F33d256
    c0664=(33.d0*a4+340.d0*a2b2-8157.d0*b4)*F77d192
    c0665=-(3.d0*a2-179.d0*b2)*F105d32
    c0666=-F1155d52
    c6060=(1155.d0*bC-15015.d0*a2bA+5005.d0*a4b8+57915.d0*a6b6+44044.d0*a8b4 &
        -32760.d0*aAb2-73920.d0*aC)*F3328
    c6061=-(495.d0*bA-8525.d0*a2b8+15345.d0*a4b6+30789.d0*a6b4-14960.d0*a8b2 &
        -64440.d0*aA)*F7d768
    c6062=(25.d0*b8-2325.d0*a2b6+14679.d0*a4b4-3495.d0*a6b2-54380.d0*a8) &
        *F77d1280
    c6063=(135.d0*b6-2177.d0*a2b4-1631.d0*a4b2+43737.d0*a6)*F33d256
    c6064=(33.d0*b4+340.d0*a2b2-8157.d0*a4)*F77d192
    c6065=-(3.d0*b2-179.d0*a2)*F105d32
    c6066=-F1155d52
    c6600=-(73920.d0*aCpbC-1954680.d0*a8pb8*a2b2+10886876.d0*a4pb4*a4b4 &
        -18763173.d0*a6b6)*F3328
    c6601=-(1080.d0*a8pb8-16040.d0*a4pb4*a2b2+39107.d0*a4b4)*a2pb2*F7d768
    c6602=(220.d0*a8pb8-4665.d0*a4pb4*a2b2+14679.d0*a4b4)*F77d1280
    c6603=(27.d0*a4pb4-244.d0*a2b2)*a2pb2*F165d256
    c6604=(3.d0*a4pb4+155.d0*a2b2)*F385d768
    c6605=-a2pb2*F1155d256
    c6606=F1155d3328
    c6420=-(44352.d0*aC-775320.d0*aAb2+2014012.d0*a8b4+299871.d0*a6b6 &
        -1850849.d0*a4b8+559104.d0*a2bA-27720.d0*bC)*F5d3328
    c6421=(12456.d0*aA-247280.d0*a8b2+865953.d0*a6b4-663003.d0*a4b6 &
        +134530.d0*a2b8-5328.d0*bA)*F35d768
    c6422=(172.d0*a8+735.d0*a6b2-12915.d0*a4b4+5415.d0*a2b6-335.d0*b8) &
        *F385d256
    c6423=-(933.d0*a6-12635.d0*a4b2-1715.d0*a2b4+465.d0*b6)*F165d256
    c6424=-(654.d0*a4+5125.d0*a2b2-165.d0*b4)*F385d768
    c6425=(487.d0*a2+305.d0*b2)*F105d256
    c6426=-F10395d832
    c2640=-(8316.d0*aC-132951.d0*aAb2+218218.d0*a8b4+400257.d0*a6b6 &
        -172172.d0*a4b8-377832.d0*a2bA+44352.d0*bC)*F5d3328
    c2641=(549.d0*aA-11275.d0*a8b2+35739.d0*a6b4+4851.d0*a4b6-49456.d0*a2b8 &
        +5112.d0*bA)*F175d768
    c2642=(55.d0*a8+735.d0*a6b2-12915.d0*a4b4+26241.d0*a2b6-2012.d0*b8) &
        *F385d256
    c2643=-(465.d0*a6-12635.d0*a4b2+46879.d0*a2b4+699.d0*b6)*F165d256
    c2644=-(1005.d0*a4-12230.d0*a2b2-5547.d0*b4)*F385d768
    c2645=-(3885.d0*a2+13440.d0*b2)*F16
    c2646=F17325d416
    c4260=(27720.d0*aC-161616.d0*aAb2-335335.d0*a8b4-199485.d0*a6b6 &
        +55055.d0*a4b8+83265.d0*a2bA-8316.d0*bC)*F5d3328
    c4261=-(18432.d0*aA-134530.d0*a8b2-178695.d0*a6b4-24255.d0*a4b6 &
        +56375.d0*a2b8-4383.d0*bA)*F35d768
    c4262=(1849.d0*a8-20091.d0*a6b2-12915.d0*a4b4+5415.d0*a2b6-218.d0*b8) &
        *F385d256
    c4263=-(699.d0*a6-61229.d0*a4b2-1715.d0*a2b4+933.d0*b6)*F165d256
    c4264=-(1509.d0*a4+5620.d0*a2b2-129.d0*b4)*F385d192
    c4265=(355.d0*a2+173.d0*b2)*F105d32
    c4266=-F3465d52
    c4620=-(44352.d0*bC-775320.d0*a2bA+2014012.d0*a4b8+299871.d0*a6b6 &
        -1850849.d0*a8b4+559104.d0*aAb2-27720.d0*aC)*F5d3328
    c4621=(12456.d0*bA-247280.d0*a2b8+865953.d0*a4b6-663003.d0*a6b4 &
        +134530.d0*a8b2-5328.d0*aA)*F35d768
    c4622=(172.d0*b8+735.d0*a2b6-12915.d0*a4b4+5415.d0*a6b2-335.d0*a8) &
        *F385d256
    c4623=-(933.d0*b6-12635.d0*a2b4-1715.d0*a4b2+465.d0*a6)*F165d256
    c4624=-(654.d0*b4+5125.d0*a2b2-165.d0*a4)*F385d768
    c4625=(487.d0*b2+305.d0*a2)*F105d256
    c4626=-F10395d832
    c6240=-(8316.d0*bC-132951.d0*a2bA+218218.d0*a4b8+400257.d0*a6b6 &
        -172172.d0*a8b4-377832.d0*aAb2+44352.d0*aC)*F5d3328
    c6241=(549.d0*bA-11275.d0*a2b8+35739.d0*a4b6+4851.d0*a6b4-49456.d0*a8b2 &
        +5112.d0*aA)*F175d768
    c6242=(55.d0*b8+735.d0*a2b6-12915.d0*a4b4+26241.d0*a6b2-2012.d0*a8) &
        *F385d256
    c6243=-(465.d0*b6-12635.d0*a2b4+46879.d0*a4b2+699.d0*a6)*F165d256
    c6244=-(1005.d0*b4-12230.d0*a2b2-5547.d0*a4)*F385d768
    c6245=-(3885.d0*b2+13440.d0*a2)*F16
    c6246=F17325d416
    c2460=(27720.d0*bC-161616.d0*a2bA-335335.d0*a4b8-199485.d0*a6b6 &
        +55055.d0*a8b4+83265.d0*aAb2-8316.d0*aC)*F5d3328
    c2461=-(18432.d0*bA-134530.d0*a2b8-178695.d0*a4b6-24255.d0*a6b4 &
        +56375.d0*a8b2-4383.d0*aA)*F35d768
    c2462=(1849.d0*b8-20091.d0*a2b6-12915.d0*a4b4+5415.d0*a6b2-218.d0*a8) &
        *F385d256
    c2463=-(699.d0*b6-61229.d0*a2b4-1715.d0*a4b2+933.d0*a6)*F165d256
    c2464=-(1509.d0*b4+5620.d0*a2b2-129.d0*a4)*F385d192
    c2465=(355.d0*b2+173.d0*a2)*F105d32
    c2466=-F3465d52
    c4440=(1512.d0*aCpbC-19656.d0*a8pb8*a2b2+14833.d0*a4pb4*a4b4 &
        +54522.d0*a6b6)*F275d6656
    c4441=-(324.d0*aApbA-5125.d0*a6pb6*a2b2+9225.d0*a2pb2*a4b4)*F1925d768
    c4442=(163.d0*a8pb8-6150.d0*a4pb4*a2b2+25830.d0*a4b4)*F1925d512
    c4443=(699.d0*a6pb6-7175.d0*a2pb2*a2b2)*F825d256
    c4444=(489.d0*a4pb4+10250.d0*a2b2)*F1925d1536
    c4445=-a2pb2*F51975d64
    c4446=F51975d832
endif
!
if(n.ge.14) then
    aE=aC*a2;bE=bC*b2;aCb2=aC*b2;aAb4=aA*b4;a8b6=a8*b6;a6b8=a6*b8;a4bA=a4*bA
    a2bC=a2*bC;aEpbE=aE+bE
!
    cE000=(2048.d0*aE-35840.d0*aCb2+209664.d0*aAb4-457600.d0*a8b6 &
        +400400.d0*a6b8-137592.d0*a4bA+16170.d0*a2bC-429.d0*bE)*F30720
    cE001=-(7168.d0*aC-46592.d0*aAb2+128128.d0*a8b4-137280.d0*a6b6 &
        +56056.d0*a4b8-7644.d0*a2bA+231.d0*bC)*F6144
    cE002=(2304.d0*aA-7040.d0*a8b2+9504.d0*a6b4-4752.d0*a4b6+770.d0*a2b8 &
        -27.d0*bA)*F91d30720
    cE003=-(3200.d0*a8-4800.d0*a6b2+3024.d0*a4b4-600.d0*a2b6+25.d0*b8) &
        *F143d30720
    cE004=(560.d0*a6-392.d0*a4b2+98.d0*a2b4-5.d0*b6)*F143d6144
    cE005=-(504.d0*a4-140.d0*a2b2+9.d0*b4)*F91d10240
    cE006=(14.d0*a2-b2)*F77d2048
    cE007=-F143d10240
    c0E00=(2048.d0*bE-35840.d0*a2bC+209664.d0*a4bA-457600.d0*a6b8 &
        +400400.d0*a8b6-137592.d0*aAb4+16170.d0*aCb2-429.d0*aE)*F30720
    c0E01=-(7168.d0*bC-46592.d0*a2bA+128128.d0*a4b8-137280.d0*a6b6 &
        +56056.d0*a8b4-7644.d0*aAb2+231.d0*aC)*F6144
    c0E02=(2304.d0*bA-7040.d0*a2b8+9504.d0*a4b6-4752.d0*a6b4+770.d0*a8b2 &
        -27.d0*aA)*F91d30720
    c0E03=-(3200.d0*b8-4800.d0*a2b6+3024.d0*a4b4-600.d0*a6b2+25.d0*a8) &
        *F143d30720
    c0E04=(560.d0*b6-392.d0*a2b4+98.d0*a4b2-5.d0*a6)*F143d6144
    c0E05=-(504.d0*b4-140.d0*a2b2+9.d0*a4)*F91d10240
    c0E06=(14.d0*b2-a2)*F77d2048
    c0E07=-F143d10240
    c00E0=-(429.d0*aCpbC+726.d0*a8pb8*a2b2+1731.d0*a4pb4*a4b4 &
        +1844.d0*a6b6)*a2pb2*F30720
    c00E1=(1617.d0*aCpbC+3822.d0*a8pb8*a2b2+7007.d0*a4pb4*a4b4 &
        +8580.d0*a6b6)*F3072
    c00E2=-(189.d0*a8pb8+196.d0*a4pb4*a2b2+398.d0*a4b4)*a2pb2*F91d3840
    c00E3=(175.d0*a8pb8+300.d0*a4pb4*a2b2+378.d0*a4b4)*F143d1920
    c00E4=-(5.d0*a4pb4+2.d0*a2b2)*a2pb2*F143d48
    c00E5=(9.d0*a4pb4+10.d0*a2b2)*F91d120
    c00E6=-a2pb2*F7d6
    c00E7=F15
    cC200=-(93184.d0*aE-2705920.d0*aCb2+17926272.d0*aAb4 &
        -41412800.d0*a8b6+37437400.d0*a6b8-13140036.d0*a4bA &
        +1567335.d0*a2bC-42042.d0*bE)*F30720
    cC201=(111104.d0*aC-2119936.d0*aAb2+7751744.d0*a8b4-9540960.d0*a6b6 &
        +4232228.d0*a4b8-609882.d0*a2bA+19173.d0*bC)*F6144
    cC202=-(12672.d0*aA-214720.d0*a8b2+432432.d0*a6b4-263736.d0*a4b6 &
        +48235.d0*a2b8-1836.d0*bA)*F91d30720
    cC203=(1600.d0*a8-103200.d0*a6b2+107352.d0*a4b4-27300.d0*a2b6+1325.d0*b8) &
        *F143d30720
    cC204=(1400.d0*a6+6076.d0*a4b2-2779.d0*a2b4+190.d0*b6)*F143d6144
    cC205=-(2268.d0*a4+1570.d0*a2b2-207.d0*b4)*F91d10240
    cC206=(83.d0*a2+8.d0*b2)*F77d2048
    cC207=-F1001d10240
    c0C20=-(3003.d0*aE-95865.d0*aCb2+619164.d0*aAb4-1001000.d0*a8b6 &
        -228800.d0*a6b8+1153152.d0*a4bA-555520.d0*a2bC+93184.d0*bE)*F30720
    c0C21=(924.d0*aC-42861.d0*aAb2+434434.d0*a8b4-1475760.d0*a6b6 &
        +1953952.d0*a4b8-1059968.d0*a2bA+270592.d0*bC)*F3072
    c0C22=(621.d0*aA-21835.d0*a8b2+168696.d0*a6b4-432432.d0*a4b6 &
        +425920.d0*a2b8-196992.d0*bA)*F91d30720
    c0C23=(475.d0*a8-13650.d0*a6b2+83916.d0*a4b4-166800.d0*a2b6+144800.d0*b8) &
        *F143d15360
    c0C24=(265.d0*a6-6139.d0*a4b2+29596.d0*a2b4-52360.d0*b6)*F143d6144
    c0C25=(306.d0*a4-5585.d0*a2b2+24066.d0*b4)*F91d5120
    c0C26=(83.d0*a2-1357.d0*b2)*F77d2048
    c0C27=F7007d5120
    c20C0=(42042.d0*aE+95865.d0*aCb2+167076.d0*aAb4+189475.d0*a8b6 &
        +135850.d0*a6b8+56511.d0*a4bA+9240.d0*a2bC-3003.d0*bE)*F30720
    c20C1=-(44781.d0*aC+87126.d0*aAb2+125411.d0*a8b4+111540.d0*a6b6 &
        +56771.d0*a4b8+12246.d0*a2bA-2739.d0*bC)*F7d6144
    c20C2=(36099.d0*aA+58135.d0*a8b2+65934.d0*a6b4+42174.d0*a4b6 &
        +11935.d0*a2b8-1701.d0*bA)*F91d7680
    c20C3=-(32725.d0*a8+41700.d0*a6b2+34398.d0*a4b4+12900.d0*a2b6-875.d0*b8) &
        *F143d3840
    c20C4=(905.d0*a6+847.d0*a4b2+427.d0*a2b4+5.d0*b6)*F143d96
    c20C5=-(1539.d0*a4+910.d0*a2b2+99.d0*b4)*F91d240
    c20C6=(151.d0*a2+31.d0*b2)*F7d12
    c20C7=-F91d30
    c2C00=-(93184.d0*bE-2705920.d0*a2bC+17926272.d0*a4bA &
        -41412800.d0*a6b8+37437400.d0*a8b6-13140036.d0*aAb4 &
        +1567335.d0*aCb2-42042.d0*aE)*F30720
    c2C01=(111104.d0*bC-2119936.d0*a2bA+7751744.d0*a4b8-9540960.d0*a6b6 &
        +4232228.d0*a8b4-609882.d0*aAb2+19173.d0*aC)*F6144
    c2C02=-(12672.d0*bA-214720.d0*a2b8+432432.d0*a4b6-263736.d0*a6b4 &
        +48235.d0*a8b2-1836.d0*aA)*F91d30720
    c2C03=(1600.d0*b8-103200.d0*a2b6+107352.d0*a4b4-27300.d0*a6b2+1325.d0*a8) &
        *F143d30720
    c2C04=(1400.d0*b6+6076.d0*a2b4-2779.d0*a4b2+190.d0*a6)*F143d6144
    c2C05=-(2268.d0*b4+1570.d0*a2b2-207.d0*a4)*F91d10240
    c2C06=(83.d0*b2+8.d0*a2)*F77d2048
    c2C07=-F1001d10240
    cC020=-(3003.d0*bE-95865.d0*a2bC+619164.d0*a4bA-1001000.d0*a6b8 &
        -228800.d0*a8b6+1153152.d0*aAb4-555520.d0*aCb2+93184.d0*aE)*F30720
    cC021=(924.d0*bC-42861.d0*a2bA+434434.d0*a4b8-1475760.d0*a6b6 &
        +1953952.d0*a8b4-1059968.d0*aAb2+270592.d0*aC)*F3072
    cC022=(621.d0*bA-21835.d0*a2b8+168696.d0*a4b6-432432.d0*a6b4 &
        +425920.d0*a8b2-196992.d0*aA)*F91d30720
    cC023=(475.d0*b8-13650.d0*a2b6+83916.d0*a4b4-166800.d0*a6b2+144800.d0*a8) &
        *F143d15360
    cC024=(265.d0*b6-6139.d0*a2b4+29596.d0*a4b2-52360.d0*a6)*F143d6144
    cC025=(306.d0*b4-5585.d0*a2b2+24066.d0*a4)*F91d5120
    cC026=(83.d0*b2-1357.d0*a2)*F77d2048
    cC027=F7007d5120
    c02C0=(42042.d0*bE+95865.d0*a2bC+167076.d0*a4bA+189475.d0*a6b8 &
        +135850.d0*a8b6+56511.d0*aAb4+9240.d0*aCb2-3003.d0*aE)*F30720
    c02C1=-(44781.d0*bC+87126.d0*a2bA+125411.d0*a4b8+111540.d0*a6b6 &
        +56771.d0*a8b4+12246.d0*aAb2-2739.d0*aC)*F7d6144
    c02C2=(36099.d0*bA+58135.d0*a2b8+65934.d0*a4b6+42174.d0*a6b4 &
        +11935.d0*a8b2-1701.d0*aA)*F91d7680
    c02C3=-(32725.d0*b8+41700.d0*a2b6+34398.d0*a4b4+12900.d0*a6b2-875.d0*a8) &
        *F143d3840
    c02C4=(905.d0*b6+847.d0*a2b4+427.d0*a4b2+5.d0*a6)*F143d96
    c02C5=-(1539.d0*b4+910.d0*a2b2+99.d0*a4)*F91d240
    c02C6=(151.d0*b2+31.d0*a2)*F7d12
    c02C7=-F91d30
    cA400=(256256.d0*aE-8426880.d0*aCb2+61440288.d0*aAb4-150264400.d0*a8b6 &
        +140990850.d0*a6b8-50795199.d0*a4bA+6176940.d0*a2bC-168168.d0*bE) &
        *F10240
    cA401=-(108416.d0*aC-2859584.d0*aAb2+14910896.d0*a8b4-22256520.d0*a6b6 &
        +11138127.d0*a4b8-1743378.d0*a2bA+58212.d0*bC)*F2048
    cA402=-(224.d0*aA+14960.d0*a8b2-127116.d0*a6b4+121638.d0*a4b6 &
        -28380.d0*a2b8+1263.d0*bA)*F273d10240
    cA403=(5200.d0*a8-66600.d0*a6b2-10458.d0*a4b4+17700.d0*a2b6-1425.d0*b8) &
        *F143d10240
    cA404=(70.d0*a6+8351.d0*a4b2-1484.d0*a2b4-10.d0*b6)*F143d2048
    cA405=-(421.d0*a4+1090.d0*a2b2-74.d0*b4)*F819d10240
    cA406=(68.d0*a2+23.d0*b2)*F231d2048
    cA407=-F3003d10240
    c0A40=-(3003.d0*aE-78540.d0*aCb2+344799.d0*aAb4-50050.d0*a8b6 &
        -743600.d0*a6b8+61152.d0*a4bA+542080.d0*a2bC-256256.d0*bE)*F10240
    c0A41=(5313.d0*aC-178542.d0*aAb2+1194193.d0*a8b4-1904760.d0*a6b6 &
        -816816.d0*a4b8+2859584.d0*a2bA-1685376.d0*bC)*F2048
    c0A42=(333.d0*aA-5830.d0*a8b2-8217.d0*a6b4+190674.d0*a4b6-409640.d0*a2b8 &
        +337584.d0*bA)*F91d5120
    c0A43=-(25.d0*a8-8850.d0*a6b2+116109.d0*a4b4-389100.d0*a2b6+525400.d0*b8) &
        *F143d5120
    c0A44=-(95.d0*a6-3612.d0*a4b2+25963.d0*a2b4-65730.d0*b6)*F429d2048
    c0A45=-(1263.d0*a4-31930.d0*a2b2+186063.d0*b4)*F273d10240
    c0A46=-(9.d0*a2-191.d0*b2)*F1617d512
    c0A47=-F21021d1280
    c40A0=-(168168.d0*aE+291060.d0*aCb2+344799.d0*aAb4+203775.d0*a8b6 &
        +7150.d0*a6b8-60606.d0*a4bA-26565.d0*a2bC+3003.d0*bE)*F10240
    c40A1=(617694.d0*aC+871689.d0*aAb2+774774.d0*a8b4+253110.d0*a6b6 &
        -106106.d0*a4b8-89271.d0*a2bA+7854.d0*bC)*F1024
    c40A2=-(558189.d0*aA+611985.d0*a8b2+364914.d0*a6b4+16434.d0*a4b6 &
        -65615.d0*a2b8+3789.d0*bA)*F91d10240
    c40A3=(492975.d0*a8+389100.d0*a6b2+121338.d0*a4b4-33300.d0*a2b6 &
        +175.d0*b8)*F143d5120
    c40A4=-(13135.d0*a6+6517.d0*a4b2+357.d0*a2b4-65.d0*b6)*F143d128
    c40A5=(21099.d0*a4+4910.d0*a2b2-21.d0*b4)*F91d320
    c40A6=-(171.d0*a2+11.d0*b2)*F77d16
    c40A7=F1001d40
    c4A00=(256256.d0*bE-8426880.d0*a2bC+61440288.d0*a4bA-150264400.d0*a6b8 &
        +140990850.d0*a8b6-50795199.d0*aAb4+6176940.d0*aCb2-168168.d0*aE) &
        *F10240
    c4A01=-(108416.d0*bC-2859584.d0*a2bA+14910896.d0*a4b8-22256520.d0*a6b6 &
        +11138127.d0*a8b4-1743378.d0*aAb2+58212.d0*aC)*F2048
    c4A02=-(224.d0*bA+14960.d0*a2b8-127116.d0*a4b6+121638.d0*a6b4 &
        -28380.d0*a8b2+1263.d0*aA)*F273d10240
    c4A03=(5200.d0*b8-66600.d0*a2b6-10458.d0*a4b4+17700.d0*a6b2-1425.d0*a8) &
        *F143d10240
    c4A04=(70.d0*b6+8351.d0*a2b4-1484.d0*a4b2-10.d0*a6)*F143d2048
    c4A05=-(421.d0*b4+1090.d0*a2b2-74.d0*a4)*F819d10240
    c4A06=(68.d0*b2+23.d0*a2)*F231d2048
    c4A07=-F3003d10240
    cA040=-(3003.d0*bE-78540.d0*a2bC+344799.d0*a4bA-50050.d0*a6b8 &
        -743600.d0*a8b6+61152.d0*aAb4+542080.d0*aCb2-256256.d0*aE)*F10240
    cA041=(5313.d0*bC-178542.d0*a2bA+1194193.d0*a4b8-1904760.d0*a6b6 &
        -816816.d0*a8b4+2859584.d0*aAb2-1685376.d0*aC)*F2048
    cA042=(333.d0*bA-5830.d0*a2b8-8217.d0*a4b6+190674.d0*a6b4-409640.d0*a8b2 &
        +337584.d0*aA)*F91d5120
    cA043=-(25.d0*b8-8850.d0*a2b6+116109.d0*a4b4-389100.d0*a6b2+525400.d0*a8) &
        *F143d5120
    cA044=-(95.d0*b6-3612.d0*a2b4+25963.d0*a4b2-65730.d0*a6)*F429d2048
    cA045=-(1263.d0*b4-31930.d0*a2b2+186063.d0*a4)*F273d10240
    cA046=-(9.d0*b2-191.d0*a2)*F1617d512
    cA047=-F21021d1280
    c04A0=-(168168.d0*bE+291060.d0*a2bC+344799.d0*a4bA+203775.d0*a6b8 &
        +7150.d0*a8b6-60606.d0*aAb4-26565.d0*aCb2+3003.d0*aE)*F10240
    c04A1=(617694.d0*bC+871689.d0*a2bA+774774.d0*a4b8+253110.d0*a6b6 &
        -106106.d0*a8b4-89271.d0*aAb2+7854.d0*aC)*F1024
    c04A2=-(558189.d0*bA+611985.d0*a2b8+364914.d0*a4b6+16434.d0*a6b4 &
        -65615.d0*a8b2+3789.d0*aA)*F91d10240
    c04A3=(492975.d0*b8+389100.d0*a2b6+121338.d0*a4b4-33300.d0*a6b2 &
        +175.d0*a8)*F143d5120
    c04A4=-(13135.d0*b6+6517.d0*a2b4+357.d0*a4b2-65.d0*a6)*F143d128
    c04A5=(21099.d0*b4+4910.d0*a2b2-21.d0*a4)*F91d320
    c04A6=-(171.d0*b2+11.d0*a2)*F77d16
    c04A7=F1001d40
    cA220=(256256.d0*aE-4484480.d0*aCb2+12868128.d0*aAb4-4747600.d0*a8b6 &
        -11161150.d0*a6b8+7845201.d0*a4bA-1290135.d0*a2bC+42042.d0*bE) &
        *F5120
    cA221=-(896896.d0*aC-14740544.d0*aAb2+40536496.d0*a8b4 &
        -38181000.d0*a6b6+13140127.d0*a4b8-1478568.d0*a2bA &
        +36267.d0*bC)*F1024
    cA222=(141408.d0*aA-2227280.d0*a8b2+3612708.d0*a6b4-1806354.d0*a4b6 &
        +275165.d0*a2b8-8829.d0*bA)*F91d5120
    cA223=-(16600.d0*a8-667500.d0*a6b2+574749.d0*a4b4-123600.d0*a2b6 &
        +5150.d0*b8)*F143d2560
    cA224=-(15610.d0*a6+91889.d0*a4b2-35021.d0*a2b4+2060.d0*b6)*F143d1024
    cA225=(28737.d0*a4+27080.d0*a2b2-2943.d0*b4)*F273d5120
    cA226=-(1117.d0*a2+157.d0*b2)*F231d1024
    cA227=F21021d2560
    c2A20=(256256.d0*bE-4484480.d0*a2bC+12868128.d0*a4bA-4747600.d0*a6b8 &
        -11161150.d0*a8b6+7845201.d0*aAb4-1290135.d0*aCb2+42042.d0*aE) &
        *F5120
    c2A21=-(896896.d0*bC-14740544.d0*a2bA+40536496.d0*a4b8 &
        -38181000.d0*a6b6+13140127.d0*a8b4-1478568.d0*aAb2 &
        +36267.d0*aC)*F1024
    c2A22=(141408.d0*bA-2227280.d0*a2b8+3612708.d0*a4b6-1806354.d0*a6b4 &
        +275165.d0*a8b2-8829.d0*aA)*F91d5120
    c2A23=-(16600.d0*b8-667500.d0*a2b6+574749.d0*a4b4-123600.d0*a6b2 &
        +5150.d0*a8)*F143d2560
    c2A24=-(15610.d0*b6+91889.d0*a2b4-35021.d0*a4b2+2060.d0*a6)*F143d1024
    c2A25=(28737.d0*b4+27080.d0*a2b2-2943.d0*a4)*F273d5120
    c2A26=-(1117.d0*b2+157.d0*a2)*F231d1024
    c2A27=F21021d2560
    c22A0=(42042.d0*aCpbC-223377.d0*a8pb8*a2b2-580062.d0*a4pb4*a4b4 &
        -892838.d0*a6b6)*a2pb2*F5120
    c22A1=-(258027.d0*aCpbC-1478568.d0*a8pb8*a2b2-5008003.d0*a4pb4*a4b4 &
        -7069920.d0*a6b6)*F1024
    c22A2=(86211.d0*a8pb8-808196.d0*a4pb4*a2b2-998158.d0*a4b4)*a2pb2*F91d5120
    c22A3=-(39025.d0*a8pb8-667500.d0*a4pb4*a2b2-1149498.d0*a4b4)*F143d2560
    c22A4=-(415.d0*a4pb4+17302.d0*a2b2)*a2pb2*F143d64
    c22A5=(4419.d0*a4pb4+25310.d0*a2b2)*F91d160
    c22A6=-a2pb2*F7007d8
    c22A7=F1001d20
    c8600=-(128128.d0*aE-4459840.d0*aCb2+34424208.d0*aAb4-88700040.d0*a8b6 &
        +86763105.d0*a6b8-32299722.d0*a4bA+4031720.d0*a2bC-112112.d0*bE) &
        *F2048
    c8601=(704.d0*aC-34528.d0*aAb2+291720.d0*a8b4-664092.d0*a6b6 &
        +433719.d0*a4b8-81172.d0*a2bA+3080.d0*bC)*F35d2048
    c8602=(1872.d0*aA-39160.d0*a8b2+132462.d0*a6b4-64152.d0*a4b6 &
        +4895.d0*a2b8+126.d0*bA)*F91d2048
    c8603=(88.d0*a8+516.d0*a6b2-8694.d0*a4b4+2820.d0*a2b6-131.d0*b8) &
        *F715d2048
    c8604=-(131.d0*a6-2114.d0*a4b2-770.d0*a2b4+104.d0*b6)*F715d2048
    c8605=-(513.d0*a4+4120.d0*a2b2+18.d0*b4)*F91d2048
    c8606=(53.d0*a2+38.d0*b2)*F385d2048
    c8607=-F1001d2048
    c0860=-(1001.d0*aE-20405.d0*aCb2+46683.d0*aAb4+93665.d0*a8b6 &
        -62920.d0*a6b8-170352.d0*a4bA-24640.d0*a2bC+128128.d0*bE) &
        *F2048
    c0861=(1463.d0*aC-37492.d0*aAb2+151151.d0*a8b4+36894.d0*a6b6 &
        -356356.d0*a4b8-120848.d0*a2bA+445984.d0*bC)*F5d1024
    c0862=-(9.d0*aA-3025.d0*a8b2+34155.d0*a6b4-66231.d0*a4b6 &
        -56100.d0*a2b8+189144.d0*bA)*F91d1024
    c0863=-(130.d0*a8-3525.d0*a6b2+10206.d0*a4b4+40635.d0*a2b6 &
        -155070.d0*b8)*F143d512
    c0864=-(131.d0*a6-623.d0*a4b2-21231.d0*a2b4+121347.d0*b6)*F715d2048
    c0865=(9.d0*a4-2230.d0*a2b2+25353.d0*b4)*F637d1024
    c0866=(5.d0*a2-187.d0*b2)*F2695d256
    c0867=F7007d128
    c6080=(112112.d0*aE+107800.d0*aCb2+11466.d0*aAb4-93665.d0*a8b6 &
        -74360.d0*a6b8-1638.d0*a4bA+14630.d0*a2bC-1001.d0*bE)*F2048
    c6081=-(806344.d0*aC+568204.d0*aAb2-89089.d0*a8b4-403260.d0*a6b6 &
        -110110.d0*a4b8+74984.d0*a2bA-4081.d0*bC)*F5d2048
    c6082=(354942.d0*aA+166815.d0*a8b2-64152.d0*a6b4-68310.d0*a4b6 &
        +16610.d0*a2b8-513.d0*bA)*F91d2048
    c6083=-(606735.d0*a8+162540.d0*a6b2-84294.d0*a4b4-2580.d0*a2b6 &
        +655.d0*b8)*F143d2048
    c6084=(15507.d0*a6+1785.d0*a4b2-623.d0*a2b4+11.d0*b6)*F715d256
    c6085=-(23643.d0*a4+830.d0*a2b2-117.d0*b4)*F91d128
    c6086=(181.d0*a2+b2)*F385d32
    c6087=-F1001d16
    c6800=-(128128.d0*bE-4459840.d0*a2bC+34424208.d0*a4bA-88700040.d0*a6b8 &
        +86763105.d0*a8b6-32299722.d0*aAb4+4031720.d0*aCb2-112112.d0*aE) &
        *F2048
    c6801=(704.d0*bC-34528.d0*a2bA+291720.d0*a4b8-664092.d0*a6b6 &
        +433719.d0*a8b4-81172.d0*aAb2+3080.d0*aC)*F35d2048
    c6802=(1872.d0*bA-39160.d0*a2b8+132462.d0*a4b6-64152.d0*a6b4 &
        +4895.d0*a8b2+126.d0*aA)*F91d2048
    c6803=(88.d0*b8+516.d0*a2b6-8694.d0*a4b4+2820.d0*a6b2-131.d0*a8) &
        *F715d2048
    c6804=-(131.d0*b6-2114.d0*a2b4-770.d0*a4b2+104.d0*a6)*F715d2048
    c6805=-(513.d0*b4+4120.d0*a2b2+18.d0*a4)*F91d2048
    c6806=(53.d0*b2+38.d0*a2)*F385d2048
    c6807=-F1001d2048
    c8060=-(1001.d0*bE-20405.d0*a2bC+46683.d0*a4bA+93665.d0*a6b8 &
        -62920.d0*a8b6-170352.d0*aAb4-24640.d0*aCb2+128128.d0*aE) &
        *F2048
    c8061=(1463.d0*bC-37492.d0*a2bA+151151.d0*a4b8+36894.d0*a6b6 &
        -356356.d0*a8b4-120848.d0*aAb2+445984.d0*aC)*F5d1024
    c8062=-(9.d0*bA-3025.d0*a2b8+34155.d0*a4b6-66231.d0*a6b4 &
        -56100.d0*a8b2+189144.d0*aA)*F91d1024
    c8063=-(130.d0*b8-3525.d0*a2b6+10206.d0*a4b4+40635.d0*a6b2 &
        -155070.d0*a8)*F143d512
    c8064=-(131.d0*b6-623.d0*a2b4-21231.d0*a4b2+121347.d0*a6)*F715d2048
    c8065=(9.d0*b4-2230.d0*a2b2+25353.d0*a4)*F637d1024
    c8066=(5.d0*b2-187.d0*a2)*F2695d256
    c8067=F7007d128
    c0680=(112112.d0*bE+107800.d0*a2bC+11466.d0*a4bA-93665.d0*a6b8 &
        -74360.d0*a8b6-1638.d0*aAb4+14630.d0*aCb2-1001.d0*aE)*F2048
    c0681=-(806344.d0*bC+568204.d0*a2bA-89089.d0*a4b8-403260.d0*a6b6 &
        -110110.d0*a8b4+74984.d0*aAb2-4081.d0*aC)*F5d2048
    c0682=(354942.d0*bA+166815.d0*a2b8-64152.d0*a4b6-68310.d0*a6b4 &
        +16610.d0*a8b2-513.d0*aA)*F91d2048
    c0683=-(606735.d0*b8+162540.d0*a2b6-84294.d0*a4b4-2580.d0*a6b2 &
        +655.d0*a8)*F143d2048
    c0684=(15507.d0*b6+1785.d0*a2b4-623.d0*a4b2+11.d0*a6)*F715d256
    c0685=-(23643.d0*b4+830.d0*a2b2-117.d0*a4)*F91d128
    c0686=(181.d0*b2+a2)*F385d32
    c0687=-F1001d16
    c8420=-(128128.d0*aE-2981440.d0*aCb2+12199824.d0*aAb4-7293000.d0*a8b6 &
        -10842975.d0*a6b8+9113013.d0*a4bA-1627780.d0*a2bC+56056.d0*bE) &
        *F3d2048
    c8421=(150304.d0*aC-3685136.d0*aAb2+17261244.d0*a8b4-21763170.d0*a6b6 &
        +9117108.d0*a4b8-1194557.d0*a2bA+33418.d0*bC)*F15d1024
    c8422=-(7344.d0*aA-330440.d0*a8b2+1806354.d0*a6b4-1415502.d0*a4b6 &
        +279895.d0*a2b8-10737.d0*bA)*F273d2048
    c8423=-(8900.d0*a8-93450.d0*a6b2-124362.d0*a4b4+61800.d0*a2b6-3775.d0*b8) &
        *F429d1024
    c8424=(445.d0*a6-35623.d0*a4b2+602.d0*a2b4+550.d0*b6)*F2145d2048
    c8425=(6966.d0*a4+25015.d0*a2b2-954.d0*b4)*F273d1024
    c8426=-(877.d0*a2+397.d0*b2)*F1155d2048
    c8427=F21021d1024
    c2840=(14014.d0*aE-337645.d0*aCb2+1267812.d0*aAb4+318175.d0*a8b6 &
        -2545400.d0*a6b8-668304.d0*a4bA+1503040.d0*a2bC-128128.d0*bE) &
        *F3d2048
    c2841=-(4367.d0*aC-130078.d0*aAb2+727727.d0*a8b4-763620.d0*a6b6 &
        -859144.d0*a4b8+1052896.d0*a2bA-85184.d0*bC)*F105d2048
    c2842=-(954.d0*aA-2365.d0*a8b2-195426.d0*a6b4+903177.d0*a4b6 &
        -948420.d0*a2b8+67032.d0*bA)*F273d1024
    c2843=(1375.d0*a8-61800.d0*a6b2+450387.d0*a4b4-760950.d0*a2b6 &
        +25500.d0*b8)*F429d1024
    c2844=(1510.d0*a6-35623.d0*a4b2+127512.d0*a2b4+15165.d0*b6)*F2145d2048
    c2845=(10737.d0*a4-131270.d0*a2b2-100143.d0*b4)*F273d2048
    c2846=(31.d0*a2+151.d0*b2)*F8085d512
    c2847=-F21021d256
    c4280=-(56056.d0*aE-334180.d0*aCb2-977067.d0*aAb4-1079650.d0*a8b6 &
        -393250.d0*a6b8+173628.d0*a4bA+152845.d0*a2bC-14014.d0*bE) &
        *F3d2048
    c4281=(325556.d0*aC-2389114.d0*aAb2-5094089.d0*a8b4-3534960.d0*a6b6 &
        +86086.d0*a4b8+910546.d0*a2bA-67529.d0*bC)*F15d2048
    c4282=-(100143.d0*aA-1001880.d0*a8b2-1415502.d0*a6b4-390852.d0*a4b6 &
        +279895.d0*a2b8-13932.d0*bA)*F273d2048
    c4283=(75825.d0*a8-1521900.d0*a6b2-1149498.d0*a4b4+186900.d0*a2b6 &
        +2225.d0*b8)*F429d2048
    c4284=(1275.d0*a6+30177.d0*a4b2+5257.d0*a2b4-445.d0*b6)*F2145d256
    c4285=-(8379.d0*a4+25310.d0*a2b2+459.d0*b4)*F273d128
    c4286=(121.d0*a2+61.d0*b2)*F1155d32
    c4287=-F3003d16
    c4820=-(128128.d0*bE-2981440.d0*a2bC+12199824.d0*a4bA-7293000.d0*a6b8 &
        -10842975.d0*a8b6+9113013.d0*aAb4-1627780.d0*aCb2+56056.d0*aE) &
        *F3d2048
    c4821=(150304.d0*bC-3685136.d0*a2bA+17261244.d0*a4b8-21763170.d0*a6b6 &
        +9117108.d0*a8b4-1194557.d0*aAb2+33418.d0*aC)*F15d1024
    c4822=-(7344.d0*bA-330440.d0*a2b8+1806354.d0*a4b6-1415502.d0*a6b4 &
        +279895.d0*a8b2-10737.d0*aA)*F273d2048
    c4823=-(8900.d0*b8-93450.d0*a2b6-124362.d0*a4b4+61800.d0*a6b2 &
        -3775.d0*a8)*F429d1024
    c4824=(445.d0*b6-35623.d0*a2b4+602.d0*a4b2+550.d0*a6)*F2145d2048
    c4825=(6966.d0*b4+25015.d0*a2b2-954.d0*a4)*F273d1024
    c4826=-(877.d0*b2+397.d0*a2)*F1155d2048
    c4827=F21021d1024
    c8240=(14014.d0*bE-337645.d0*a2bC+1267812.d0*a4bA+318175.d0*a6b8 &
        -2545400.d0*a8b6-668304.d0*aAb4+1503040.d0*aCb2-128128.d0*aE) &
        *F3d2048
    c8241=-(4367.d0*bC-130078.d0*a2bA+727727.d0*a4b8-763620.d0*a6b6 &
        -859144.d0*a8b4+1052896.d0*aAb2-85184.d0*aC)*F105d2048
    c8242=-(954.d0*bA-2365.d0*a2b8-195426.d0*a4b6+903177.d0*a6b4 &
        -948420.d0*a8b2+67032.d0*aA)*F273d1024
    c8243=(1375.d0*b8-61800.d0*a2b6+450387.d0*a4b4-760950.d0*a6b2 &
        +25500.d0*a8)*F429d1024
    c8244=(1510.d0*b6-35623.d0*a2b4+127512.d0*a4b2+15165.d0*a6)*F2145d2048
    c8245=(10737.d0*b4-131270.d0*a2b2-100143.d0*a4)*F273d2048
    c8246=(31.d0*b2+151.d0*a2)*F8085d512
    c8247=-F21021d256
    c2480=-(56056.d0*bE-334180.d0*a2bC-977067.d0*a4bA-1079650.d0*a6b8 &
        -393250.d0*a8b6+173628.d0*aAb4+152845.d0*aCb2-14014.d0*aE) &
        *F3d2048
    c2481=(325556.d0*bC-2389114.d0*a2bA-5094089.d0*a4b8-3534960.d0*a6b6 &
        +86086.d0*a8b4+910546.d0*aAb2-67529.d0*aC)*F15d2048
    c2482=-(100143.d0*bA-1001880.d0*a2b8-1415502.d0*a4b6-390852.d0*a6b4 &
        +279895.d0*a8b2-13932.d0*aA)*F273d2048
    c2483=(75825.d0*b8-1521900.d0*a2b6-1149498.d0*a4b4+186900.d0*a6b2 &
        +2225.d0*a8)*F429d2048
    c2484=(1275.d0*b6+30177.d0*a2b4+5257.d0*a4b2-445.d0*a6)*F2145d256
    c2485=-(8379.d0*b4+25310.d0*a2b2+459.d0*a4)*F273d128
    c2486=(121.d0*b2+61.d0*a2)*F1155d32
    c2487=-F3003d16
    c2660=(2002.d0*aE-35035.d0*aCb2+48321.d0*aAb4+168025.d0*a8b6 &
        +30745.d0*a6b8-181818.d0*a4bA-132440.d0*a2bC+16016.d0*bE) &
        *F7d512
    c2661=-(7007.d0*aC-149968.d0*aAb2+412412.d0*a8b4+477048.d0*a6b6 &
        -623623.d0*a4b8-809900.d0*a2bA+85624.d0*bC)*F35d512
    c2662=(531.d0*aA-22660.d0*a8b2+136620.d0*a6b4-68310.d0*a4b6 &
        -279015.d0*a2b8+23346.d0*bA)*F637d512
    c2663=(235.d0*a8-3336.d0*a6b2-8694.d0*a4b4+65016.d0*a2b6-2709.d0*b8) &
        *F5005d512
    c2664=(43.d0*a6+4361.d0*a4b2-35511.d0*a2b4-2709.d0*b6)*F5005d512
    c2665=-(999.d0*a4-22250.d0*a2b2-11673.d0*b4)*F637d256
    c2666=-(43.d0*a2+139.d0*b2)*F2695d64
    c2667=F7007d32
    c6260=(2002.d0*bE-35035.d0*a2bC+48321.d0*a4bA+168025.d0*a6b8 &
        +30745.d0*a8b6-181818.d0*aAb4-132440.d0*aCb2+16016.d0*aE) &
        *F7d512
    c6261=-(7007.d0*bC-149968.d0*a2bA+412412.d0*a4b8+477048.d0*a6b6 &
        -623623.d0*a8b4-809900.d0*aAb2+85624.d0*aC)*F35d512
    c6262=(531.d0*bA-22660.d0*a2b8+136620.d0*a4b6-68310.d0*a6b4 &
        -279015.d0*a8b2+23346.d0*aA)*F637d512
    c6263=(235.d0*b8-3336.d0*a2b6-8694.d0*a4b4+65016.d0*a6b2-2709.d0*a8) &
        *F5005d512
    c6264=(43.d0*b6+4361.d0*a2b4-35511.d0*a4b2-2709.d0*a6)*F5005d512
    c6265=-(999.d0*b4-22250.d0*a2b2-11673.d0*a4)*F637d256
    c6266=-(43.d0*b2+139.d0*a2)*F2695d64
    c6267=F7007d32
    c6620=(2288.d0*aCpbC-63448.d0*a8pb8*a2b2+366946.d0*a4pb4*a4b4 &
        -643651.d0*a6b6)*a2pb2*F49d512
    c6621=-(3784.d0*aCpbC-115700.d0*a8pb8*a2b2+725439.d0*a4pb4*a4b4 &
        -1328184.d0*a6b6)*F245d512
    c6622=-(1998.d0*a8pb8-36263.d0*a4pb4*a2b2+104573.d0*a4b4)*a2pb2 &
        *F637d512
    c6623=(43.d0*a8pb8-3336.d0*a4pb4*a2b2+17388.d0*a4b4)*F5005d512
    c6624=(235.d0*a4pb4-3119.d0*a2b2)*a2pb2*F5005d512
    c6625=(531.d0*a4pb4+8240.d0*a2b2)*F637d512
    c6626=-a2pb2*F245245d512
    c6627=F7007d256
    c6440=(16016.d0*aE-280280.d0*aCb2+525798.d0*aAb4+797225.d0*a8b6 &
        -386100.d0*a6b8-503139.d0*a4bA+170940.d0*a2bC-8008.d0*bE) &
        *F21d1024
    c6441=-(56056.d0*aC-1106924.d0*aAb2+3044041.d0*a8b4+986700.d0*a6b6 &
        -2385383.d0*a4b8+553462.d0*a2bA-21868.d0*bC)*F105d1024
    c6442=(5778.d0*aA-167255.d0*a8b2+715968.d0*a6b4-357984.d0*a4b6 &
        +36190.d0*a2b8-249.d0*bA)*F1911d1024
    c6443=(5575.d0*a8-34500.d0*a6b2-227808.d0*a4b4+69000.d0*a2b6 &
        -2875.d0*b8)*F3003d1024
    c6444=-(540.d0*a6-16681.d0*a4b2-4606.d0*a2b4+575.d0*b6)*F15015d1024
    c6445=-(5529.d0*a4+30410.d0*a2b2+249.d0*b4)*F1911d1024
    c6446=(111.d0*a2+71.d0*b2)*F8085d256
    c6447=-F21021d128
    c4640=(16016.d0*bE-280280.d0*a2bC+525798.d0*a4bA+797225.d0*a6b8 &
        -386100.d0*a8b6-503139.d0*aAb4+170940.d0*aCb2-8008.d0*aE) &
        *F21d1024
    c4641=-(56056.d0*bC-1106924.d0*a2bA+3044041.d0*a4b8+986700.d0*a6b6 &
        -2385383.d0*a8b4+553462.d0*aAb2-21868.d0*aC)*F105d1024
    c4642=(5778.d0*bA-167255.d0*a2b8+715968.d0*a4b6-357984.d0*a6b4 &
        +36190.d0*a8b2-249.d0*aA)*F1911d1024
    c4643=(5575.d0*b8-34500.d0*a2b6-227808.d0*a4b4+69000.d0*a6b2 &
        -2875.d0*a8)*F3003d1024
    c4644=-(540.d0*b6-16681.d0*a2b4-4606.d0*a4b2+575.d0*a6)*F15015d1024
    c4645=-(5529.d0*b4+30410.d0*a2b2+249.d0*a4)*F1911d1024
    c4646=(111.d0*b2+71.d0*a2)*F8085d256
    c4647=-F21021d128
    c4460=-(8008.d0*aCpbC-117348.d0*a8pb8*a2b2+140007.d0*a4pb4*a4b4 &
        +271118.d0*a6b6)*a2pb2*F21d1024
    c4461=(17094.d0*aCpbC-276731.d0*a8pb8*a2b2+329329.d0*a4pb4*a4b4 &
        +986700.d0*a6b6)*F105d512
    c4462=-(5529.d0*a8pb8-136594.d0*a4pb4*a2b2+494578.d0*a4b4)*a2pb2 &
        *F1911d1024
    c4463=-(225.d0*a8pb8+2875.d0*a4pb4*a2b2-37968.d0*a4b4)*F9009d256
    c4464=(1115.d0*a4pb4-22402.d0*a2b2)*a2pb2*F15015d1024
    c4465=(2889.d0*a4pb4+30410.d0*a2b2)*F1911d512
    c4466=-a2pb2*F735735d128
    c4467=F21021d64
endif
!
if(n.ge.16) then
    aG=aE*a2;bG=bE*b2;aEb2=aE*b2;aCb4=aC*b4;aAb6=aA*b6;a8b8=a8*b8;a6bA=a6*bA
    a4bC=a4*bC;a2bE=a2*bE;aGpbG=aG+bG
!
    cG000=(294912.d0*aG-6684672.d0*aEb2+52641792.d0*aCb4 &
        -162938880.d0*aAb6+217817600.d0*a8b8-128314368.d0*a6bA &
        +31667328.d0*a4bC-2800512.d0*a2bE+57915.d0*bG)*F5013504
    cG001=-(16384.d0*aE-143360.d0*aCb2+559104.d0*aAb4-915200.d0*a8b6 &
        +640640.d0*a6b8-183456.d0*a4bA+18480.d0*a2bC-429.d0*bE)*F12288
    cG002=(86016.d0*aC-372736.d0*aAb2+768768.d0*a8b4-658944.d0*a6b6 &
        +224224.d0*a4b8-26208.d0*a2bA+693.d0*bC)*F8192
    cG003=-(30720.d0*aA-70400.d0*a8b2+76032.d0*a6b4-31680.d0*a4b6 &
        +4400.d0*a2b8-135.d0*bA)*F13d12288
    cG004=(44800.d0*a8-53760.d0*a6b2+28224.d0*a4b4-4800.d0*a2b6 &
        +175.d0*b8)*F143d147456
    cG005=-(2688.d0*a6-1568.d0*a4b2+336.d0*a2b4-15.d0*b6)*F39d4096
    cG006=(672.d0*a4-160.d0*a2b2+9.d0*b4)*F77d8192
    cG007=-(16.d0*a2-b2)*F143d4096
    cG008=F6435d557056
    c0G00=(294912.d0*bG-6684672.d0*a2bE+52641792.d0*a4bC &
        -162938880.d0*a6bA+217817600.d0*a8b8-128314368.d0*aAb6 &
        +31667328.d0*aCb4-2800512.d0*aEb2+57915.d0*aG)*F5013504
    c0G01=-(16384.d0*bE-143360.d0*a2bC+559104.d0*a4bA-915200.d0*a6b8 &
        +640640.d0*a8b6-183456.d0*aAb4+18480.d0*aCb2-429.d0*aE)*F12288
    c0G02=(86016.d0*bC-372736.d0*a2bA+768768.d0*a4b8-658944.d0*a6b6 &
        +224224.d0*a8b4-26208.d0*aAb2+693.d0*aC)*F8192
    c0G03=-(30720.d0*bA-70400.d0*a2b8+76032.d0*a4b6-31680.d0*a6b4 &
        +4400.d0*a8b2-135.d0*aA)*F13d12288
    c0G04=(44800.d0*b8-53760.d0*a2b6+28224.d0*a4b4-4800.d0*a6b2 &
        +175.d0*a8)*F143d147456
    c0G05=-(2688.d0*b6-1568.d0*a2b4+336.d0*a4b2-15.d0*a6)*F39d4096
    c0G06=(672.d0*b4-160.d0*a2b2+9.d0*a4)*F77d8192
    c0G07=-(16.d0*b2-a2)*F143d4096
    c0G08=F6435d557056
    c00G0=(57915.d0*aGpbG+175032.d0*aCpbC*a2b2+424116.d0*a8pb8*a4b4 &
        +716040.d0*a4pb4*a6b6+850850.d0*a8b8)*F5013504
    c00G1=-(429.d0*aEpbE+1155.d0*aApbA*a2b2+2457.d0*a6pb6*a4b4 &
        +3575.d0*a2pb2*a6b6)*F768
    c00G2=(1617.d0*aCpbC+3822.d0*a8pb8*a2b2+7007.d0*a4pb4*a4b4 &
        +8580.d0*a6b6)*F256
    c00G3=-(189.d0*aApbA+385.d0*a6pb6*a2b2+594.d0*a2pb2*a4b4)*F13d96
    c00G4=(175.d0*a8pb8+300.d0*a4pb4*a2b2+378.d0*a4b4)*F143d576
    c00G5=-(5.d0*a6pb6+7.d0*a2pb2*a2b2)*F13d2
    c00G6=(9.d0*a4pb4+10.d0*a2b2)*F7d6
    c00G7=-a2pb2*F4d3
    c00G8=F17
    cE200=-(737280.d0*aG-27365376.d0*aEb2+243468288.d0*aCb4 &
        -796872960.d0*aAb6+1099978880.d0*a8b8-661620960.d0*a6bA &
        +165687984.d0*a4bC-14812083.d0*a2bE+308880.d0*bG)*F208896
    cE201=(118784.d0*aE-2867200.d0*aCb2+14746368.d0*aAb4 &
        -27639040.d0*a8b6+20980960.d0*a6b8-6342336.d0*a4bA &
        +664125.d0*a2bC-15873.d0*bE)*F4096
    cE202=-(129024.d0*aC-2539264.d0*aAb2+7687680.d0*a8b4 &
        -7989696.d0*a6b6+3059056.d0*a4b8-387387.d0*a2bA &
        +10857.d0*bC)*F3d4096
    cE203=(26880.d0*aA-689920.d0*a8b2+1197504.d0*a6b4-633600.d0*a4b6 &
        +102025.d0*a2b8-3465.d0*bA)*F13d4096
    cE204=(4480.d0*a8+97440.d0*a6b2-89712.d0*a4b4+20175.d0*a2b6 &
        -875.d0*b8)*F143d12288
    cE205=-(3360.d0*a6+8512.d0*a4b2-3507.d0*a2b4+215.d0*b6)*F117d4096
    cE206=(48.d0*a4+25.d0*a2b2-3.d0*b4)*F3003d4096
    cE207=-(37.d0*a2+3.d0*b2)*F429d4096
    cE208=F6435d69632
    c0E20=(19305.d0*aG-809523.d0*aEb2+7351344.d0*aCb4-20049120.d0*aAb6 &
        +10890880.d0*a8b8+17821440.d0*a6bA-19740672.d0*a4bC &
        +6057984.d0*a2bE-737280.d0*bG)*F208896
    c0E21=-(1287.d0*aE-75075.d0*aCb2+995904.d0*aAb4-4644640.d0*a8b6 &
        +8968960.d0*a6b8-7617792.d0*a4bA+2867200.d0*a2bC &
        -536576.d0*bE)*F4096
    c0E22=-(429.d0*aC-19539.d0*aAb2+203632.d0*a8b4-741312.d0*a6b6 &
        +1098240.d0*a4b8-702208.d0*a2bA+227328.d0*bC)*F21d4096
    c0E23=-(1935.d0*aA-73975.d0*a8b2+633600.d0*a6b4-1843776.d0*a4b6 &
        +2126080.d0*a2b8-1201920.d0*bA)*F13d4096
    c0E24=-(125.d0*a8-3975.d0*a6b2+27504.d0*a4b4-62880.d0*a2b6 &
        +64640.d0*b8)*F1001d12288
    c0E25=-(5.d0*a6-129.d0*a4b2+704.d0*a2b4-1440.d0*b6)*F9009d4096
    c0E26=-(141.d0*a4-2875.d0*a2b2+14064.d0*b4)*F231d4096
    c0E27=-(37.d0*a2-677.d0*b2)*F429d4096
    c0E28=-F6435d4352
    c20E0=-(308880.d0*aG+809523.d0*aEb2+1661121.d0*aCb4+2297295.d0*aAb6 &
        +2127125.d0*a8b8+1282905.d0*a6bA+459459.d0*a4bC+65637.d0*a2bE &
        -19305.d0*bG)*F208896
    c20E1=(290433.d0*aE+664125.d0*aCb2+1162161.d0*aAb4+1326325.d0*a8b6 &
        +961675.d0*a6b8+410319.d0*a4bA+75075.d0*a2bC-15873.d0*bE)*F4096
    c20E2=-(67683.d0*aC+132132.d0*aAb2+191191.d0*a8b4+171600.d0*a6b6 &
        +89089.d0*a4b8+20748.d0*a2bA-3003.d0*bC)*F3d256
    c20E3=(31185.d0*aA+50435.d0*a8b2+57618.d0*a6b4+37422.d0*a4b6 &
        +11165.d0*a2b8-945.d0*bA)*F13d128
    c20E4=-(707.d0*a8+906.d0*a6b2+756.d0*a4b4+294.d0*a2b6-7.d0*b8)*F715d96
    c20E5=(1565.d0*a6+1477.d0*a4b2+763.d0*a2b4+35.d0*b6)*F39d16
    c20E6=-(333.d0*a4+200.d0*a2b2+27.d0*b4)*F7d2
    c20E7=131.d0*a2+29.d0*b2
    c20E8=-F60d17
    c2E00=-(737280.d0*bG-27365376.d0*a2bE+243468288.d0*a4bC &
        -796872960.d0*a6bA+1099978880.d0*a8b8-661620960.d0*aAb6 &
        +165687984.d0*aCb4-14812083.d0*aEb2+308880.d0*aG)*F208896
    c2E01=(118784.d0*bE-2867200.d0*a2bC+14746368.d0*a4bA &
        -27639040.d0*a6b8+20980960.d0*a8b6-6342336.d0*aAb4 &
        +664125.d0*aCb2-15873.d0*aE)*F4096
    c2E02=-(129024.d0*bC-2539264.d0*a2bA+7687680.d0*a4b8 &
        -7989696.d0*a6b6+3059056.d0*a8b4-387387.d0*aAb2 &
        +10857.d0*aC)*F3d4096
    c2E03=(26880.d0*bA-689920.d0*a2b8+1197504.d0*a4b6-633600.d0*a6b4 &
        +102025.d0*a8b2-3465.d0*aA)*F13d4096
    c2E04=(4480.d0*b8+97440.d0*a2b6-89712.d0*a4b4+20175.d0*a6b2 &
        -875.d0*a8)*F143d12288
    c2E05=-(3360.d0*b6+8512.d0*a2b4-3507.d0*a4b2+215.d0*a6)*F117d4096
    c2E06=(48.d0*b4+25.d0*a2b2-3.d0*a4)*F3003d4096
    c2E07=-(37.d0*b2+3.d0*a2)*F429d4096
    c2E08=F6435d69632
    cE020=(19305.d0*bG-809523.d0*a2bE+7351344.d0*a4bC-20049120.d0*a6bA &
        +10890880.d0*a8b8+17821440.d0*aAb6-19740672.d0*aCb4 &
        +6057984.d0*aEb2-737280.d0*aG)*F208896
    cE021=-(1287.d0*bE-75075.d0*a2bC+995904.d0*a4bA-4644640.d0*a6b8 &
        +8968960.d0*a8b6-7617792.d0*aAb4+2867200.d0*aCb2 &
        -536576.d0*aE)*F4096
    cE022=-(429.d0*bC-19539.d0*a2bA+203632.d0*a4b8-741312.d0*a6b6 &
        +1098240.d0*a8b4-702208.d0*aAb2+227328.d0*aC)*F21d4096
    cE023=-(1935.d0*bA-73975.d0*a2b8+633600.d0*a4b6-1843776.d0*a6b4 &
        +2126080.d0*a8b2-1201920.d0*aA)*F13d4096
    cE024=-(125.d0*b8-3975.d0*a2b6+27504.d0*a4b4-62880.d0*a6b2 &
        +64640.d0*a8)*F1001d12288
    cE025=-(5.d0*b6-129.d0*a2b4+704.d0*a4b2-1440.d0*a6)*F9009d4096
    cE026=-(141.d0*b4-2875.d0*a2b2+14064.d0*a4)*F231d4096
    cE027=-(37.d0*b2-677.d0*a2)*F429d4096
    cE028=-F6435d4352
    c02E0=-(308880.d0*bG+809523.d0*a2bE+1661121.d0*a4bC+2297295.d0*a6bA &
        +2127125.d0*a8b8+1282905.d0*aAb6+459459.d0*aCb4+65637.d0*aEb2 &
        -19305.d0*aG)*F208896
    c02E1=(290433.d0*bE+664125.d0*a2bC+1162161.d0*a4bA+1326325.d0*a6b8 &
        +961675.d0*a8b6+410319.d0*aAb4+75075.d0*aCb2-15873.d0*aE)*F4096
    c02E2=-(67683.d0*bC+132132.d0*a2bA+191191.d0*a4b8+171600.d0*a6b6 &
        +89089.d0*a8b4+20748.d0*aAb2-3003.d0*aC)*F3d256
    c02E3=(31185.d0*bA+50435.d0*a2b8+57618.d0*a4b6+37422.d0*a6b4 &
        +11165.d0*a8b2-945.d0*aA)*F13d128
    c02E4=-(707.d0*b8+906.d0*a2b6+756.d0*a4b4+294.d0*a6b2-7.d0*a8)*F715d96
    c02E5=(1565.d0*b6+1477.d0*a2b4+763.d0*a4b2+35.d0*a6)*F39d16
    c02E6=-(333.d0*b4+200.d0*a2b2+27.d0*a4)*F7d2
    c02E7=131.d0*b2+29.d0*a2
    c02E8=-F60d17
    cC400=(16773120.d0*aG-703352832.d0*aEb2+6829998336.d0*aCb4 &
        -23549760000.d0*aAb6+33633760160.d0*a8b8-20720288160.d0*a6bA &
        +5281693263.d0*a4bC-478654176.d0*a2bE+10090080.d0*bG)*F417792
    cC401=-(559104.d0*aE-18144000.d0*aCb2+124819968.d0*aAb4 &
        -275887040.d0*a8b6+232678160.d0*a6b8-75666591.d0*a4bA &
        +8360275.d0*a2bC-208208.d0*bE)*F4096
    cC402=(284928.d0*aC-11135488.d0*aAb2+71431360.d0*a8b4 &
        -103701312.d0*a6b6+48305257.d0*a4b8-6970964.d0*a2bA &
        +214599.d0*bC)*F3d8192
    cC403=(61440.d0*aA-559680.d0*a8b2-2163744.d0*a6b4+2410650.d0*a4b6 &
        -547525.d0*a2b8+22935.d0*bA)*F13d4096
    cC404=-(25760.d0*a8-633120.d0*a6b2+78561.d0*a4b4+66600.d0*a2b6 &
        -5775.d0*b8)*F143d24576
    cC405=-(10800.d0*a6+138999.d0*a4b2-24269.d0*a2b4+90.d0*b6)*F39d4096
    cC406=(8361.d0*a4+14900.d0*a2b2-921.d0*b4)*F77d8192
    cC407=-(47.d0*a2+13.d0*b2)*F1001d4096
    cC408=F45045d139264
    c0C40=(135135.d0*aG-4798794.d0*aEb2+32833647.d0*aCb4-42962400.d0*aAb6 &
        -62622560.d0*a8b8+81469440.d0*a6bA+43593984.d0*a4bC &
        -57028608.d0*a2bE+16773120.d0*bG)*F417792
    c0C41=-(13013.d0*aE-573650.d0*aCb2+5420961.d0*aAb4-15089360.d0*a8b6 &
        +7275840.d0*a6b8+16703232.d0*a4bA-18144000.d0*a2bC &
        +6895616.d0*bE)*F4096
    c0C42=-(23639.d0*aC-630994.d0*aAb2+1248247.d0*a8b4+18752448.d0*a6b6 &
        -71431360.d0*a4b8+83213312.d0*a2bA-44640512.d0*bC)*F3d8192
    c0C43=-(135.d0*aA+61050.d0*a8b2-1205325.d0*a6b4+5982768.d0*a4b6 &
        -10611040.d0*a2b8+8880000.d0*bA)*F13d2048
    c0C44=(525.d0*a8-27150.d0*a6b2+276381.d0*a4b4-887520.d0*a2b6 &
        +1257760.d0*b8)*F1573d24576
    c0C45=(695.d0*a6-24374.d0*a4b2+176379.d0*a2b4-473520.d0*b6)*F429d4096
    c0C46=(8361.d0*a4-217150.d0*a2b2+1344969.d0*b4)*F77d8192
    c0C47=(13.d0*a2-293.d0*b2)*F1001d256
    c0C48=F105105d4352
    c40C0=(3363360.d0*aG+7079072.d0*aEb2+10944549.d0*aCb4 &
        +10137270.d0*aAb6+4679675.d0*a8b8-119340.d0*a6bA &
        -1205589.d0*a4bC-442442.d0*a2bE+45045.d0*bG)*F139264
    c40C1=-(4692688.d0*aE+8360275.d0*aCb2+10456446.d0*aAb4 &
        +7117825.d0*a8b6+1587300.d0*a6b8-946491.d0*a4bA-573650.d0*a2bC &
        +47047.d0*bE)*F4096
    c40C2=(34520871.d0*aC+50444394.d0*aAb2+48305257.d0*a8b4 &
        +20892300.d0*a6b6-1248247.d0*a4b8-3613974.d0*a2bA &
        +214599.d0*bC)*F3d8192
    c40C3=-(976635.d0*aA+1118645.d0*a8b2+747846.d0*a6b4+135234.d0*a4b6 &
        -72545.d0*a2b8+2025.d0*bA)*F13d256
    c40C4=(86471.d0*a8+72348.d0*a6b2+28098.d0*a4b4-1908.d0*a2b6 &
        -161.d0*b8)*F715d768
    c40C5=-(23125.d0*a6+12502.d0*a4b2+1673.d0*a2b4-80.d0*b6)*F39d16
    c40C6=(24911.d0*a4+6750.d0*a2b2+159.d0*b4)*F21d32
    c40C7=-(37.d0*a2+3.d0*b2)*F91d2
    c40C8=F1365d34
    c4C00=(16773120.d0*bG-703352832.d0*a2bE+6829998336.d0*a4bC &
        -23549760000.d0*a6bA+33633760160.d0*a8b8-20720288160.d0*aAb6 &
        +5281693263.d0*aCb4-478654176.d0*aEb2+10090080.d0*aG)*F417792
    c4C01=-(559104.d0*bE-18144000.d0*a2bC+124819968.d0*a4bA &
        -275887040.d0*a6b8+232678160.d0*a8b6-75666591.d0*aAb4 &
        +8360275.d0*aCb2-208208.d0*aE)*F4096
    c4C02=(284928.d0*bC-11135488.d0*a2bA+71431360.d0*a4b8 &
        -103701312.d0*a6b6+48305257.d0*a8b4-6970964.d0*aAb2 &
        +214599.d0*aC)*F3d8192
    c4C03=(61440.d0*bA-559680.d0*a2b8-2163744.d0*a4b6+2410650.d0*a6b4 &
        -547525.d0*a8b2+22935.d0*aA)*F13d4096
    c4C04=-(25760.d0*b8-633120.d0*a2b6+78561.d0*a4b4+66600.d0*a6b2 &
        -5775.d0*a8)*F143d24576
    c4C05=-(10800.d0*b6+138999.d0*a2b4-24269.d0*a4b2+90.d0*a6)*F39d4096
    c4C06=(8361.d0*b4+14900.d0*a2b2-921.d0*a4)*F77d8192
    c4C07=-(47.d0*b2+13.d0*a2)*F1001d4096
    c4C08=F45045d139264
    cC040=(135135.d0*bG-4798794.d0*a2bE+32833647.d0*a4bC-42962400.d0*a6bA &
        -62622560.d0*a8b8+81469440.d0*aAb6+43593984.d0*aCb4 &
        -57028608.d0*aEb2+16773120.d0*aG)*F417792
    cC041=-(13013.d0*bE-573650.d0*a2bC+5420961.d0*a4bA-15089360.d0*a6b8 &
        +7275840.d0*a8b6+16703232.d0*aAb4-18144000.d0*aCb2 &
        +6895616.d0*aE)*F4096
    cC042=-(23639.d0*bC-630994.d0*a2bA+1248247.d0*a4b8+18752448.d0*a6b6 &
        -71431360.d0*a8b4+83213312.d0*aAb2-44640512.d0*aC)*F3d8192
    cC043=-(135.d0*bA+61050.d0*a2b8-1205325.d0*a4b6+5982768.d0*a6b4 &
        -10611040.d0*a8b2+8880000.d0*aA)*F13d2048
    cC044=(525.d0*b8-27150.d0*a2b6+276381.d0*a4b4-887520.d0*a6b2 &
        +1257760.d0*a8)*F1573d24576
    cC045=(695.d0*b6-24374.d0*a2b4+176379.d0*a4b2-473520.d0*a6)*F429d4096
    cC046=(8361.d0*b4-217150.d0*a2b2+1344969.d0*a4)*F77d8192
    cC047=(13.d0*b2-293.d0*a2)*F1001d256
    cC048=F105105d4352
    c04C0=(3363360.d0*bG+7079072.d0*a2bE+10944549.d0*a4bC &
        +10137270.d0*a6bA+4679675.d0*a8b8-119340.d0*aAb6 &
        -1205589.d0*aCb4-442442.d0*aEb2+45045.d0*aG)*F139264
    c04C1=-(4692688.d0*bE+8360275.d0*a2bC+10456446.d0*a4bA &
        +7117825.d0*a6b8+1587300.d0*a8b6-946491.d0*aAb4-573650.d0*aCb2 &
        +47047.d0*aE)*F4096
    c04C2=(34520871.d0*bC+50444394.d0*a2bA+48305257.d0*a4b8 &
        +20892300.d0*a6b6-1248247.d0*a8b4-3613974.d0*aAb2 &
        +214599.d0*aC)*F3d8192
    c04C3=-(976635.d0*bA+1118645.d0*a2b8+747846.d0*a4b6+135234.d0*a6b4 &
        -72545.d0*a8b2+2025.d0*aA)*F13d256
    c04C4=(86471.d0*b8+72348.d0*a2b6+28098.d0*a4b4-1908.d0*a6b2 &
        -161.d0*a8)*F715d768
    c04C5=-(23125.d0*b6+12502.d0*a2b4+1673.d0*a4b2-80.d0*a6)*F39d16
    c04C6=(24911.d0*b4+6750.d0*a2b2+159.d0*a4)*F21d32
    c04C7=-(37.d0*b2+3.d0*a2)*F91d2
    c04C8=F1365d34
    cC220=(3354624.d0*aG-76038144.d0*aEb2+333123840.d0*aCb4 &
        -373231872.d0*aAb6-160640480.d0*a8b8+390671424.d0*a6bA &
        -153494649.d0*a4bC+17612595.d0*a2bE-432432.d0*bG)*F5d208896
    cC221=-(1490944.d0*aE-30410240.d0*aCb2+118599936.d0*aAb4 &
        -171966080.d0*a8b6+102639680.d0*a6b8-24630606.d0*a4bA &
        +2054745.d0*a2bC-39039.d0*bE)*F5d4096
    cC222=(2177280.d0*aC-39533312.d0*aAb2+97056960.d0*a8b4 &
        -83191680.d0*a6b6+26691665.d0*a4b8-2867865.d0*a2bA &
        +68838.d0*bC)*F15d4096
    cC223=-(562944.d0*aA-13228160.d0*a8b2+19198080.d0*a6b4 &
        -8638740.d0*a4b6+1199825.d0*a2b8-35541.d0*bA)*F65d4096
    cC224=-(13216.d0*a8+430656.d0*a6b2-335979.d0*a4b4+65445.d0*a2b6 &
        -2492.d0*b8)*F3575d12288
    cC225=(5952.d0*a6+19138.d0*a4b2-6685.d0*a2b4+359.d0*b6)*F6435d4096
    cC226=-(4343.d0*a4+2965.d0*a2b2-298.d0*b4)*F3465d4096
    cC227=(115.d0*a2+13.d0*b2)*F15015d4096
    cC228=-F45045d4352
    c2C20=(3354624.d0*bG-76038144.d0*a2bE+333123840.d0*a4bC &
        -373231872.d0*a6bA-160640480.d0*a8b8+390671424.d0*aAb6 &
        -153494649.d0*aCb4+17612595.d0*aEb2-432432.d0*aG)*F5d208896
    c2C21=-(1490944.d0*bE-30410240.d0*a2bC+118599936.d0*a4bA &
        -171966080.d0*a6b8+102639680.d0*a8b6-24630606.d0*aAb4 &
        +2054745.d0*aCb2-39039.d0*aE)*F5d4096
    c2C22=(2177280.d0*bC-39533312.d0*a2bA+97056960.d0*a4b8 &
        -83191680.d0*a6b6+26691665.d0*a8b4-2867865.d0*aAb2 &
        +68838.d0*aC)*F15d4096
    c2C23=-(562944.d0*bA-13228160.d0*a2b8+19198080.d0*a4b6 &
        -8638740.d0*a6b4+1199825.d0*a8b2-35541.d0*aA)*F65d4096
    c2C24=-(13216.d0*b8+430656.d0*a2b6-335979.d0*a4b4+65445.d0*a6b2 &
        -2492.d0*a8)*F3575d12288
    c2C25=(5952.d0*b6+19138.d0*a2b4-6685.d0*a4b2+359.d0*a6)*F6435d4096
    c2C26=-(4343.d0*b4+2965.d0*a2b2-298.d0*a4)*F3465d4096
    c2C27=(115.d0*b2+13.d0*a2)*F15015d4096
    c2C28=-F45045d4352
    c22C0=-(39312.d0*aGpbG-180999.d0*aCpbC*a2b2-957474.d0*a8pb8*a4b4 &
        -2142153.d0*a4pb4*a6b6-2753660.d0*a8b8)*F55d208896
    c22C1=(6279.d0*aCpbC-43638.d0*a8pb8*a2b2-112791.d0*a4pb4*a4b4 &
        -170804.d0*a6b6)*a2pb2*F275d4096
    c22C2=-(91203.d0*aCpbC-746382.d0*a8pb8*a2b2-2426515.d0*a4pb4*a4b4 &
        -3403140.d0*a6b6)*F165d4096
    c22C3=(837.d0*a8pb8-12052.d0*a4pb4*a2b2-15218.d0*a4b4)*a2pb2*F715d64
    c22C4=-(413.d0*a8pb8-22548.d0*a4pb4*a2b2-38178.d0*a4b4)*F3575d384
    c22C5=-(733.d0*a4pb4+11146.d0*a2b2)*a2pb2*F195d16
    c22C6=(729.d0*a4pb4+3394.d0*a2b2)*F175d16
    c22C7=-a2pb2*1820.d0
    c22C8=F1365d17
    cA600=-(8386560.d0*aG-371874048.d0*aEb2+3804192000.d0*aCb4 &
        -13716016704.d0*aAb6+20288348080.d0*a8b8-12846921165.d0*a6bA &
        +3347123472.d0*a4bC-308756448.d0*a2bE+6604416.d0*bG)*F11d626688
    cA601=(163072.d0*aE-6567680.d0*aCb2+56801472.d0*aAb4 &
        -156865280.d0*a8b6+155554685.d0*a6b8-56926233.d0*a4bA &
        +6874560.d0*a2bC-183456.d0*bE)*F11d12288
    cA602=(86016.d0*aC-1857856.d0*aAb2+5669664.d0*a8b4+1572714.d0*a6b6 &
        -4523519.d0*a4b8+1101555.d0*a2bA-45360.d0*bC)*F11d4096
    cA603=-(12480.d0*aA-499840.d0*a8b2+2605878.d0*a6b4-1426590.d0*a4b6 &
        +168685.d0*a2b8-2565.d0*bA)*F143d12288
    cA604=-(6160.d0*a8-61095.d0*a6b2-187614.d0*a4b4+62070.d0*a2b6 &
        -2765.d0*b8)*F1573d36864
    cA605=(285.d0*a6-23401.d0*a4b2-4242.d0*a2b4+642.d0*b6)*F429d4096
    cA606=(4587.d0*a4+19745.d0*a2b2-54.d0*b4)*F77d4096
    cA607=-(77.d0*a2+43.d0*b2)*F1001d4096
    cA608=F45045d69632
    c0A60=(36855.d0*aG-1072071.d0*aEb2+4912677.d0*aCb4+1700595.d0*aAb6 &
        -14974960.d0*a8b8-8274240.d0*a6bA+13160448.d0*a4bC &
        +8316672.d0*a2bE-8386560.d0*bG)*F11d626688
    c0A61=-(11739.d0*aE-414645.d0*aCb2+2737917.d0*aAb4-2912195.d0*a8b6 &
        -6497920.d0*a6b8+5573568.d0*a4bA+6567680.d0*a2bC &
        -7291648.d0*bE)*F11d12288
    c0A62=-(189.d0*aC+82719.d0*aAb2-1490489.d0*a8b4+5646069.d0*a6b6 &
        -2834832.d0*a4b8-9466912.d0*a2bA+12432000.d0*bC)*F11d2048
    c0A63=(2889.d0*aA-113795.d0*a8b2+713295.d0*a6b4+181467.d0*a4b6 &
        -6033280.d0*a2b8+10343904.d0*bA)*F143d6144
    c0A64=(2765.d0*a8-46005.d0*a6b2-284697.d0*a4b4+3263385.d0*a2b6 &
        -8345680.d0*b8)*F1573d36864
    c0A65=(95.d0*a6+9415.d0*a4b2-162183.d0*a2b4+717665.d0*b6)*F1287d4096
    c0A66=-(135.d0*a4-6820.d0*a2b2+65109.d0*b4)*F231d256
    c0A67=-(a2-33.d0*b2)*F21021d128
    c0A68=-F63063d544
    c60A0=-(6604416.d0*aG+9356256.d0*aEb2+6940080.d0*aCb4 &
        -1700595.d0*aAb6-6721715.d0*a8b8-3830814.d0*a6bA &
        +57834.d0*a4bC+598689.d0*a2bE-36855.d0*bG)*F11d626688
    c60A1=(6054048.d0*aE+6874560.d0*aCb2+3304665.d0*aAb4 &
        -2192905.d0*a8b6-2958670.d0*a6b8-496314.d0*a4bA &
        +414645.d0*a2bC-21021.d0*bE)*F11d12288
    c60A2=-(21876624.d0*aC+18975411.d0*aAb2+4523519.d0*a8b4 &
        -6181890.d0*a6b6-2980978.d0*a4b8+912639.d0*a2bA &
        -32109.d0*bC)*F11d4096
    c60A3=(19376955.d0*aA+11965745.d0*a8b2+362934.d0*a6b4 &
        -2605878.d0*a4b6+224015.d0*a2b8+2565.d0*bA)*F143d12288
    c60A4=-(521605.d0*a8+205680.d0*a6b2-22302.d0*a4b4-8520.d0*a2b6 &
        +385.d0*b8)*F1573d2304
    c60A5=(107749.d0*a6+22757.d0*a4b2-2233.d0*a2b4-65.d0*b6)*F143d64
    c60A6=-(41625.d0*a4+3665.d0*a2b2-144.d0*b4)*F77d48
    c60A7=(313.d0*a2+7.d0*b2)*F1001d48
    c60A8=-F5005d34
    c6A00=-(8386560.d0*bG-371874048.d0*a2bE+3804192000.d0*a4bC &
        -13716016704.d0*a6bA+20288348080.d0*a8b8-12846921165.d0*aAb6 &
        +3347123472.d0*aCb4-308756448.d0*aEb2+6604416.d0*aG)*F11d626688
    c6A01=(163072.d0*bE-6567680.d0*a2bC+56801472.d0*a4bA &
        -156865280.d0*a6b8+155554685.d0*a8b6-56926233.d0*aAb4 &
        +6874560.d0*aCb2-183456.d0*aE)*F11d12288
    c6A02=(86016.d0*bC-1857856.d0*a2bA+5669664.d0*a4b8+1572714.d0*a6b6 &
        -4523519.d0*a8b4+1101555.d0*aAb2-45360.d0*aC)*F11d4096
    c6A03=-(12480.d0*bA-499840.d0*a2b8+2605878.d0*a4b6-1426590.d0*a6b4 &
        +168685.d0*a8b2-2565.d0*aA)*F143d12288
    c6A04=-(6160.d0*b8-61095.d0*a2b6-187614.d0*a4b4+62070.d0*a6b2 &
        -2765.d0*a8)*F1573d36864
    c6A05=(285.d0*b6-23401.d0*a2b4-4242.d0*a4b2+642.d0*a6)*F429d4096
    c6A06=(4587.d0*b4+19745.d0*a2b2-54.d0*a4)*F77d4096
    c6A07=-(77.d0*b2+43.d0*a2)*F1001d4096
    c6A08=F45045d69632
    cA060=(36855.d0*bG-1072071.d0*a2bE+4912677.d0*a4bC+1700595.d0*a6bA &
        -14974960.d0*a8b8-8274240.d0*aAb6+13160448.d0*aCb4 &
        +8316672.d0*aEb2-8386560.d0*aG)*F11d626688
    cA061=-(11739.d0*bE-414645.d0*a2bC+2737917.d0*a4bA-2912195.d0*a6b8 &
        -6497920.d0*a8b6+5573568.d0*aAb4+6567680.d0*aCb2 &
        -7291648.d0*aE)*F11d12288
    cA062=-(189.d0*bC+82719.d0*a2bA-1490489.d0*a4b8+5646069.d0*a6b6 &
        -2834832.d0*a8b4-9466912.d0*aAb2+12432000.d0*aC)*F11d2048
    cA063=(2889.d0*bA-113795.d0*a2b8+713295.d0*a4b6+181467.d0*a6b4 &
        -6033280.d0*a8b2+10343904.d0*aA)*F143d6144
    cA064=(2765.d0*b8-46005.d0*a2b6-284697.d0*a4b4+3263385.d0*a6b2 &
        -8345680.d0*a8)*F1573d36864
    cA065=(95.d0*b6+9415.d0*a2b4-162183.d0*a4b2+717665.d0*a6)*F1287d4096
    cA066=-(135.d0*b4-6820.d0*a2b2+65109.d0*a4)*F231d256
    cA067=-(b2-33.d0*a2)*F21021d128
    cA068=-F63063d544
    c06A0=-(6604416.d0*bG+9356256.d0*a2bE+6940080.d0*a4bC &
        -1700595.d0*a6bA-6721715.d0*a8b8-3830814.d0*aAb6 &
        +57834.d0*aCb4+598689.d0*aEb2-36855.d0*aG)*F11d626688
    c06A1=(6054048.d0*bE+6874560.d0*a2bC+3304665.d0*a4bA &
        -2192905.d0*a6b8-2958670.d0*a8b6-496314.d0*aAb4 &
        +414645.d0*aCb2-21021.d0*aE)*F11d12288
    c06A2=-(21876624.d0*bC+18975411.d0*a2bA+4523519.d0*a4b8 &
        -6181890.d0*a6b6-2980978.d0*a8b4+912639.d0*aAb2 &
        -32109.d0*aC)*F11d4096
    c06A3=(19376955.d0*bA+11965745.d0*a2b8+362934.d0*a4b6 &
        -2605878.d0*a6b4+224015.d0*a8b2+2565.d0*aA)*F143d12288
    c06A4=-(521605.d0*b8+205680.d0*a2b6-22302.d0*a4b4-8520.d0*a6b2 &
        +385.d0*a8)*F1573d2304
    c06A5=(107749.d0*b6+22757.d0*a2b4-2233.d0*a4b2-65.d0*a6)*F143d64
    c06A6=-(41625.d0*b4+3665.d0*a2b2-144.d0*a4)*F77d48
    c06A7=(313.d0*b2+7.d0*a2)*F1001d48
    c06A8=-F5005d34
    cA420=-(8386560.d0*aG-250688256.d0*aEb2+1469035008.d0*aCb4 &
        -2069196480.d0*aAb6-540459920.d0*a8b8+2073741345.d0*a6bA &
        -890537571.d0*a4bC+107819712.d0*a2bE-2751840.d0*bG)*F11d208896
    cA421=(2539264.d0*aE-76025600.d0*aCb2+464912448.d0*aAb4 &
        -870995840.d0*a8b6+618295535.d0*a6b8-169368381.d0*a4bA &
        +15788850.d0*a2bC-331968.d0*bE)*F11d4096
    cA422=-(1284864.d0*aC-42695744.d0*aAb2+242642400.d0*a8b4 &
        -303240366.d0*a6b6+122298176.d0*a4b8-15405117.d0*a2bA &
        +416997.d0*bC)*F33d4096
    cA423=-(306240.d0*aA-858880.d0*a8b2-26011854.d0*a6b4 &
        +21596850.d0*a4b6-4128575.d0*a2b8+150435.d0*bA)*F143d4096
    cA424=(108080.d0*a8-2204835.d0*a6b2-702387.d0*a4b4+510150.d0*a2b6 &
        -31150.d0*b8)*F1573d12288
    cA425=(20175.d0*a6+395003.d0*a4b2-27328.d0*a2b4-3030.d0*b6)*F1287d4096
    cA426=-(114906.d0*a4+262625.d0*a2b2-10401.d0*b4)*F231d4096
    cA427=(473.d0*a2+167.d0*b2)*F9009d4096
    cA428=-F135135d4352
    c2A40=-(589680.d0*aG-19756737.d0*aEb2+123064326.d0*aCb4 &
        -120384225.d0*aAb6-262742480.d0*a8b8+203037120.d0*a6bA &
        +196584192.d0*a4bC-129502464.d0*a2bE+8386560.d0*bG)*F11d208896
    c2A41=(136773.d0*aE-5515125.d0*aCb2+46215351.d0*aAb4 &
        -105097135.d0*a8b6+11165440.d0*a6b8+128087232.d0*a4bA &
        -76025600.d0*a2bC+4915456.d0*bE)*F11d4096
    c2A42=(72807.d0*aC-1065792.d0*aAb2-11160149.d0*a8b4 &
        +112718034.d0*a6b6-242642400.d0*a4b8+154970816.d0*a2bA &
        -9601536.d0*bC)*F33d4096
    c2A43=-(13635.d0*aA-935275.d0*a8b2+10798425.d0*a6b4-34989273.d0*a4b6 &
        +33499840.d0*a2b8-1560480.d0*bA)*F143d2048
    c2A44=-(31150.d0*a8-1125975.d0*a6b2+7697088.d0*a4b4-12971235.d0*a2b6 &
        -222320.d0*b8)*F1573d12288
    c2A45=-(16715.d0*a6-395003.d0*a4b2+1447593.d0*a2b4 &
        +347535.d0*b6)*F1287d4096
    c2A46=-(59571.d0*a4-751850.d0*a2b2-831501.d0*b4)*F231d4096
    c2A47=-(19.d0*a2+121.d0*b2)*F3003d64
    c2A48=F315315d2176
    c42A0=(2751840.d0*aG-16930368.d0*aEb2-63800541.d0*aCb4 &
        -99738405.d0*aAb6-75725650.d0*a8b8-18080010.d0*a6bA &
        +11139471.d0*a4bC+6975423.d0*a2bE-589680.d0*bG)*F11d208896
    c42A1=-(2114112.d0*aE-15788850.d0*aCb2-46215351.d0*aAb4 &
        -53671475.d0*a8b6-24317150.d0*a6b8+3197376.d0*a4bA &
        +5515125.d0*a2bC-387387.d0*bE)*F11d4096
    c42A2=(5820507.d0*aC-56456127.d0*aAb2-122298176.d0*a8b4 &
        -93586350.d0*a6b6-11160149.d0*a4b8+15405117.d0*a2bA &
        -804342.d0*bC)*F33d4096
    c42A3=-(3127815.d0*aA-47561195.d0*a8b2-69978546.d0*a6b4 &
        -26011854.d0*a4b6+8084395.d0*a2b8-181575.d0*bA)*F143d4096
    c42A4=(2779.d0*a8-228408.d0*a6b2-190890.d0*a4b4+2928.d0*a2b6 &
        +1351.d0*b8)*F7865d768
    c42A5=(16255.d0*a6+186263.d0*a4b2+51317.d0*a2b4-1595.d0*b6)*F429d64
    c42A6=-(16074.d0*a4+42425.d0*a2b2+2151.d0*b4)*F77d16
    c42A7=(211.d0*a2+109.d0*b2)*F1001d16
    c42A8=-F15015d34
    c4A20=-(8386560.d0*bG-250688256.d0*a2bE+1469035008.d0*a4bC &
        -2069196480.d0*a6bA-540459920.d0*a8b8+2073741345.d0*aAb6 &
        -890537571.d0*aCb4+107819712.d0*aEb2-2751840.d0*aG)*F11d208896
    c4A21=(2539264.d0*bE-76025600.d0*a2bC+464912448.d0*a4bA &
        -870995840.d0*a6b8+618295535.d0*a8b6-169368381.d0*aAb4 &
        +15788850.d0*aCb2-331968.d0*aE)*F11d4096
    c4A22=-(1284864.d0*bC-42695744.d0*a2bA+242642400.d0*a4b8 &
        -303240366.d0*a6b6+122298176.d0*a8b4-15405117.d0*aAb2 &
        +416997.d0*aC)*F33d4096
    c4A23=-(306240.d0*bA-858880.d0*a2b8-26011854.d0*a4b6 &
        +21596850.d0*a6b4-4128575.d0*a8b2+150435.d0*aA)*F143d4096
    c4A24=(108080.d0*b8-2204835.d0*a2b6-702387.d0*a4b4+510150.d0*a6b2 &
        -31150.d0*a8)*F1573d12288
    c4A25=(20175.d0*b6+395003.d0*a2b4-27328.d0*a4b2-3030.d0*a6)*F1287d4096
    c4A26=-(114906.d0*b4+262625.d0*a2b2-10401.d0*a4)*F231d4096
    c4A27=(473.d0*b2+167.d0*a2)*F9009d4096
    c4A28=-F135135d4352
    cA240=-(589680.d0*bG-19756737.d0*a2bE+123064326.d0*a4bC &
        -120384225.d0*a6bA-262742480.d0*a8b8+203037120.d0*aAb6 &
        +196584192.d0*aCb4-129502464.d0*aEb2+8386560.d0*aG)*F11d208896
    cA241=(136773.d0*bE-5515125.d0*a2bC+46215351.d0*a4bA &
        -105097135.d0*a6b8+11165440.d0*a8b6+128087232.d0*aAb4 &
        -76025600.d0*aCb2+4915456.d0*aE)*F11d4096
    cA242=(72807.d0*bC-1065792.d0*a2bA-11160149.d0*a4b8 &
        +112718034.d0*a6b6-242642400.d0*a8b4+154970816.d0*aAb2 &
        -9601536.d0*aC)*F33d4096
    cA243=-(13635.d0*bA-935275.d0*a2b8+10798425.d0*a4b6-34989273.d0*a6b4 &
        +33499840.d0*a8b2-1560480.d0*aA)*F143d2048
    cA244=-(31150.d0*b8-1125975.d0*a2b6+7697088.d0*a4b4-12971235.d0*a6b2 &
        -222320.d0*a8)*F1573d12288
    cA245=-(16715.d0*b6-395003.d0*a2b4+1447593.d0*a4b2 &
        +347535.d0*a6)*F1287d4096
    cA246=-(59571.d0*b4-751850.d0*a2b2-831501.d0*a4)*F231d4096
    cA247=-(19.d0*b2+121.d0*a2)*F3003d64
    cA248=F315315d2176
    c24A0=(2751840.d0*bG-16930368.d0*a2bE-63800541.d0*a4bC &
        -99738405.d0*a6bA-75725650.d0*a8b8-18080010.d0*aAb6 &
        +11139471.d0*aCb4+6975423.d0*aEb2-589680.d0*aG)*F11d208896
    c24A1=-(2114112.d0*bE-15788850.d0*a2bC-46215351.d0*a4bA &
        -53671475.d0*a6b8-24317150.d0*a8b6+3197376.d0*aAb4 &
        +5515125.d0*aCb2-387387.d0*aE)*F11d4096
    c24A2=(5820507.d0*bC-56456127.d0*a2bA-122298176.d0*a4b8 &
        -93586350.d0*a6b6-11160149.d0*a8b4+15405117.d0*aAb2 &
        -804342.d0*aC)*F33d4096
    c24A3=-(3127815.d0*bA-47561195.d0*a2b8-69978546.d0*a4b6 &
        -26011854.d0*a6b4+8084395.d0*a8b2-181575.d0*aA)*F143d4096
    c24A4=(2779.d0*b8-228408.d0*a2b6-190890.d0*a4b4+2928.d0*a6b2 &
        +1351.d0*a8)*F7865d768
    c24A5=(16255.d0*b6+186263.d0*a2b4+51317.d0*a4b2-1595.d0*a6)*F429d64
    c24A6=-(16074.d0*b4+42425.d0*a2b2+2151.d0*a4)*F77d16
    c24A7=(211.d0*b2+109.d0*a2)*F1001d16
    c24A8=-F15015d34
    c0880=(6825.d0*aG-154700.d0*aEb2+353430.d0*aCb4+1047540.d0*aAb6 &
        -255255.d0*a8b8-2333760.d0*a6bA-1576512.d0*a4bC &
        +792064.d0*a2bE+1747200.d0*bG)*F33d278528
    c0881=-(2275.d0*aE-62300.d0*aCb2+242970.d0*aAb4+237380.d0*a8b6 &
        -674245.d0*a6b8-843024.d0*a4bA+330400.d0*a2bC &
        +1176448.d0*bE)*F33d4096
    c0882=(3465.d0*aC-161980.d0*aAb2+1223222.d0*a8b4-1048476.d0*a6b6 &
        -4035031.d0*a4b8+1156064.d0*a2bA+8301216.d0*bC)*F99d8192
    c0883=(1185.d0*aA-18260.d0*a8b2-120978.d0*a6b4+919116.d0*a4b6 &
        -127655.d0*a2b8-3576720.d0*bA)*F429d4096
    c0884=-(105.d0*a8-18860.d0*a6b2+169302.d0*a4b4+46420.d0*a2b6 &
        -1990935.d0*b8)*F4719d16384
    c0885=-(55.d0*a6-1351.d0*a4b2-2779.d0*a2b4+74515.d0*b6)*F1287d256
    c0886=-(207.d0*a4+2950.d0*a2b2-111177.d0*b4)*F231d256
    c0887=(a2-101.d0*b2)*F3003d32
    c0888=F225225d1088
    c8080=(6825.d0*bG-154700.d0*a2bE+353430.d0*a4bC+1047540.d0*a6bA &
        -255255.d0*a8b8-2333760.d0*aAb6-1576512.d0*aCb4 &
        +792064.d0*aEb2+1747200.d0*aG)*F33d278528
    c8081=-(2275.d0*bE-62300.d0*a2bC+242970.d0*a4bA+237380.d0*a6b8 &
        -674245.d0*a8b6-843024.d0*aAb4+330400.d0*aCb2 &
        +1176448.d0*aE)*F33d4096
    c8082=(3465.d0*bC-161980.d0*a2bA+1223222.d0*a4b8-1048476.d0*a6b6 &
        -4035031.d0*a8b4+1156064.d0*aAb2+8301216.d0*aC)*F99d8192
    c8083=(1185.d0*bA-18260.d0*a2b8-120978.d0*a4b6+919116.d0*a6b4 &
        -127655.d0*a8b2-3576720.d0*aA)*F429d4096
    c8084=-(105.d0*b8-18860.d0*a2b6+169302.d0*a4b4+46420.d0*a6b2 &
        -1990935.d0*a8)*F4719d16384
    c8085=-(55.d0*b6-1351.d0*a2b4-2779.d0*a4b2+74515.d0*a6)*F1287d256
    c8086=-(207.d0*b4+2950.d0*a2b2-111177.d0*a4)*F231d256
    c8087=(b2-101.d0*a2)*F3003d32
    c8088=F225225d1088
    c8800=(1747200.d0*aGpbG-79998464.d0*aCpbC*a2b2 &
        +846724032.d0*a8pb8*a4b4-3161820480.d0*a4pb4*a6b6 &
        +4839962985.d0*a8b8)*F33d278528
    c8801=(11648.d0*aCpbC-342048.d0*a8pb8*a2b2+2076144.d0*a4pb4*a4b4 &
        -3735659.d0*a6b6)*a2pb2*F33d4096
    c8802=-(15456.d0*aCpbC-562016.d0*a8pb8*a2b2+4035031.d0*a4pb4*a4b4 &
        -7965672.d0*a6b6)*F99d8192
    c8803=-(240.d0*a8pb8-4955.d0*a4pb4*a2b2+15953.d0*a4b4) &
        *a2pb2*F4719d4096
    c8804=-(105.d0*a8pb8+6640.d0*a4pb4*a2b2-51324.d0*a4b4)*F4719d16384
    c8805=(79.d0*a4pb4-1325.d0*a2b2)*a2pb2*F6435d4096
    c8806=(297.d0*a4pb4+3560.d0*a2b2)*F1155d8192
    c8807=-a2pb2*F75075d4096
    c8808=F225225d278528
    c8620=(5241600.d0*aG-179402496.d0*aEb2+1239755328.d0*aCb4 &
        -2181853440.d0*aAb6-197482285.d0*a8b8+2163624255.d0*a6bA &
        -1045587312.d0*a4bC+136185504.d0*a2bE-3669120.d0*bG)*F11d69632
    c8621=-(1141504.d0*aE-42089600.d0*aCb2+332562048.d0*aAb4 &
        -830792820.d0*a8b6+731307005.d0*a6b8-236076477.d0*a4bA &
        +25121600.d0*a2bC-591136.d0*bE)*F33d4096
    c8622=-(213696.d0*aC-1421056.d0*aAb2-28142114.d0*a8b4 &
        +119382978.d0*a6b6-79108029.d0*a4b8+13375999.d0*a2bA &
        -443184.d0*bC)*F99d4096
    c8623=(136320.d0*aA-3951420.d0*a8b2+16416774.d0*a6b4-3745566.d0*a4b6 &
        -608795.d0*a2b8+61095.d0*bA)*F429d4096
    c8624=(33005.d0*a8-166035.d0*a6b2-2015874.d0*a4b4+449790.d0*a2b6 &
        -11620.d0*b8)*F1573d4096
    c8625=-(3067.d0*a6-105091.d0*a4b2-47614.d0*a2b4+4138.d0*b6)*F6435d4096
    c8626=-(17919.d0*a4+109075.d0*a2b2+3996.d0*b4)*F1155d4096
    c8627=(371.d0*a2+269.d0*b2)*F15015d4096
    c8628=-F225225d4352
    c2860=-(327600.d0*aG-8609055.d0*aEb2+31985415.d0*aCb4 &
        +30501315.d0*aAb6-80235155.d0*a8b8-90380160.d0*a6bA &
        +32695488.d0*a4bC+58216704.d0*a2bE-5241600.d0*bG)*F11d69632
    c2861=(122395.d0*aE-3817625.d0*aCb2+20492745.d0*aAb4 &
        -7914335.d0*a8b6-51368460.d0*a6b8+4263168.d0*a4bA &
        +42089600.d0*a2bC-3517696.d0*bE)*F33d4096
    c2862=-(23310.d0*aC-1547455.d0*aAb2+16014999.d0*a8b4-35569677.d0*a6b6 &
        -14071057.d0*a4b8+55427008.d0*a2bA-4051488.d0*bC)*F99d2048
    c2863=-(31035.d0*aA-824615.d0*a8b2+1872783.d0*a6b4+13774959.d0*a4b6 &
        -31953570.d0*a2b8+1645440.d0*bA)*F429d2048
    c2864=-(11620.d0*a8+166035.d0*a6b2-4978827.d0*a4b4+15342105.d0*a2b6 &
        +81235.d0*b8)*F1573d4096
    c2865=(20365.d0*a6-1028923.d0*a4b2+6053243.d0*a2b4 &
        +1087795.d0*b6)*F1287d4096
    c2866=(11871.d0*a4-224300.d0*a2b2-183051.d0*b4)*F231d256
    c2867=(29.d0*a2+131.d0*b2)*F21021d128
    c2868=-F315315d544
    c6280=-(3669120.d0*aG-30147936.d0*aEb2-67807152.d0*aCb4 &
        -40505985.d0*aAb6+28248220.d0*a8b8+41152410.d0*a6bA &
        +7132860.d0*a4bC-6242145.d0*a2bE+327600.d0*bG)*F11d69632
    c6281=(2670304.d0*aE-25121600.d0*aCb2-40127997.d0*aAb4 &
        -7914335.d0*a8b6+21439990.d0*a6b8+9284730.d0*a4bA &
        -3817625.d0*a2bC+168805.d0*bE)*F33d4096
    c6282=-(6833904.d0*aC-78692159.d0*aAb2-79108029.d0*a8b4 &
        +16230786.d0*a6b6+32029998.d0*a4b8-6830915.d0*a2bA &
        +209055.d0*bC)*F99d4096
    c6283=(3263385.d0*aA-56254385.d0*a8b2-27549918.d0*a6b4 &
        +16416774.d0*a4b6-608795.d0*a2b8-46005.d0*bA)*F429d4096
    c6284=-(81235.d0*a8-17429220.d0*a6b2-1771182.d0*a4b4+1077660.d0*a2b6 &
        -33005.d0*b8)*F1573d4096
    c6285=-(8570.d0*a6+66619.d0*a4b2-854.d0*a2b4-355.d0*b6)*F1287d32
    c6286=(54261.d0*a4+93950.d0*a2b2-1431.d0*b4)*F231d64
    c6287=-(151.d0*a2+49.d0*b2)*F3003d16
    c6288=F225225d272
    c6820=(5241600.d0*bG-179402496.d0*a2bE+1239755328.d0*a4bC &
        -2181853440.d0*a6bA-197482285.d0*a8b8+2163624255.d0*aAb6 &
        -1045587312.d0*aCb4+136185504.d0*aEb2-3669120.d0*aG)*F11d69632
    c6821=-(1141504.d0*bE-42089600.d0*a2bC+332562048.d0*a4bA &
        -830792820.d0*a6b8+731307005.d0*a8b6-236076477.d0*aAb4 &
        +25121600.d0*aCb2-591136.d0*aE)*F33d4096
    c6822=-(213696.d0*bC-1421056.d0*a2bA-28142114.d0*a4b8 &
        +119382978.d0*a6b6-79108029.d0*a8b4+13375999.d0*aAb2 &
        -443184.d0*aC)*F99d4096
    c6823=(136320.d0*bA-3951420.d0*a2b8+16416774.d0*a4b6-3745566.d0*a6b4 &
        -608795.d0*a8b2+61095.d0*aA)*F429d4096
    c6824=(33005.d0*b8-166035.d0*a2b6-2015874.d0*a4b4+449790.d0*a6b2 &
        -11620.d0*a8)*F1573d4096
    c6825=-(3067.d0*b6-105091.d0*a2b4-47614.d0*a4b2+4138.d0*a6)*F6435d4096
    c6826=-(17919.d0*b4+109075.d0*a2b2+3996.d0*a4)*F1155d4096
    c6827=(371.d0*b2+269.d0*a2)*F15015d4096
    c6828=-F225225d4352
    c8260=-(327600.d0*bG-8609055.d0*a2bE+31985415.d0*a4bC &
        +30501315.d0*a6bA-80235155.d0*a8b8-90380160.d0*aAb6 &
        +32695488.d0*aCb4+58216704.d0*aEb2-5241600.d0*aG)*F11d69632
    c8261=(122395.d0*bE-3817625.d0*a2bC+20492745.d0*a4bA &
        -7914335.d0*a6b8-51368460.d0*a8b6+4263168.d0*aAb4 &
        +42089600.d0*aCb2-3517696.d0*aE)*F33d4096
    c8262=-(23310.d0*bC-1547455.d0*a2bA+16014999.d0*a4b8-35569677.d0*a6b6 &
        -14071057.d0*a8b4+55427008.d0*aAb2-4051488.d0*aC)*F99d2048
    c8263=-(31035.d0*bA-824615.d0*a2b8+1872783.d0*a4b6+13774959.d0*a6b4 &
        -31953570.d0*a8b2+1645440.d0*aA)*F429d2048
    c8264=-(11620.d0*b8+166035.d0*a2b6-4978827.d0*a4b4+15342105.d0*a6b2 &
        +81235.d0*a8)*F1573d4096
    c8265=(20365.d0*b6-1028923.d0*a2b4+6053243.d0*a4b2 &
        +1087795.d0*a6)*F1287d4096
    c8266=(11871.d0*b4-224300.d0*a2b2-183051.d0*a4)*F231d256
    c8267=(29.d0*b2+131.d0*a2)*F21021d128
    c8268=-F315315d544
    c2680=-(3669120.d0*bG-30147936.d0*a2bE-67807152.d0*a4bC &
        -40505985.d0*a6bA+28248220.d0*a8b8+41152410.d0*aAb6 &
        +7132860.d0*aCb4-6242145.d0*aEb2+327600.d0*aG)*F11d69632
    c2681=(2670304.d0*bE-25121600.d0*a2bC-40127997.d0*a4bA &
        -7914335.d0*a6b8+21439990.d0*a8b6+9284730.d0*aAb4 &
        -3817625.d0*aCb2+168805.d0*aE)*F33d4096
    c2682=-(6833904.d0*bC-78692159.d0*a2bA-79108029.d0*a4b8 &
        +16230786.d0*a6b6+32029998.d0*a8b4-6830915.d0*aAb2 &
        +209055.d0*aC)*F99d4096
    c2683=(3263385.d0*bA-56254385.d0*a2b8-27549918.d0*a4b6 &
        +16416774.d0*a6b4-608795.d0*a8b2-46005.d0*aA)*F429d4096
    c2684=-(81235.d0*b8-17429220.d0*a2b6-1771182.d0*a4b4+1077660.d0*a6b2 &
        -33005.d0*a8)*F1573d4096
    c2685=-(8570.d0*b6+66619.d0*a2b4-854.d0*a4b2-355.d0*a6)*F1287d32
    c2686=(54261.d0*b4+93950.d0*a2b2-1431.d0*a4)*F231d64
    c2687=-(151.d0*b2+49.d0*a2)*F3003d16
    c2688=F225225d272
    c8440=(3144960.d0*aG-71285760.d0*aEb2+229279680.d0*aCb4 &
        +112656960.d0*aAb6-342977635.d0*a8b8-89882910.d0*a6bA &
        +155049741.d0*a4bC-28365792.d0*a2bE+917280.d0*bG)*F55d139264
    c8441=-(698880.d0*aE-16968000.d0*aCb2+66175200.d0*aAb4 &
        -20101510.d0*a8b6-56505735.d0*a6b8+33354048.d0*a4bA &
        -4666375.d0*a2bC+129584.d0*bE)*F165d4096
    c8442=(1498560.d0*aC-44116800.d0*aAb2+214500286.d0*a8b4 &
        -183857388.d0*a6b6+43190147.d0*a4b8-2029118.d0*a2bA &
        -26187.d0*bC)*F495d8192
    c8443=(42480.d0*aA+773135.d0*a8b2-10607157.d0*a6b4+6335604.d0*a4b6 &
        -879945.d0*a2b8+22335.d0*bA)*F2145d2048
    c8444=-(141085.d0*a8-2370870.d0*a6b2-2718261.d0*a4b4+959940.d0*a2b6 &
        -42770.d0*b8)*F7865d8192
    c8445=-(22595.d0*a6+855232.d0*a4b2+78043.d0*a2b4 &
        -14890.d0*b6)*F6435d4096
    c8446=(434313.d0*a4+1333250.d0*a2b2-11223.d0*b4)*F1155d8192
    c8447=-(191.d0*a2+89.d0*b2)*F15015d256
    c8448=F1576575d4352
    c4840=(3144960.d0*bG-71285760.d0*a2bE+229279680.d0*a4bC &
        +112656960.d0*a6bA-342977635.d0*a8b8-89882910.d0*aAb6 &
        +155049741.d0*aCb4-28365792.d0*aEb2+917280.d0*aG)*F55d139264
    c4841=-(698880.d0*bE-16968000.d0*a2bC+66175200.d0*a4bA &
        -20101510.d0*a6b8-56505735.d0*a8b6+33354048.d0*aAb4 &
        -4666375.d0*aCb2+129584.d0*aE)*F165d4096
    c4842=(1498560.d0*bC-44116800.d0*a2bA+214500286.d0*a4b8 &
        -183857388.d0*a6b6+43190147.d0*a8b4-2029118.d0*aAb2 &
        -26187.d0*aC)*F495d8192
    c4843=(42480.d0*bA+773135.d0*a2b8-10607157.d0*a4b6+6335604.d0*a6b4 &
        -879945.d0*a8b2+22335.d0*aA)*F2145d2048
    c4844=-(141085.d0*b8-2370870.d0*a2b6-2718261.d0*a4b4+959940.d0*a6b2 &
        -42770.d0*a8)*F7865d8192
    c4845=-(22595.d0*b6+855232.d0*a2b4+78043.d0*a4b2 &
        -14890.d0*a6)*F6435d4096
    c4846=(434313.d0*b4+1333250.d0*a2b2-11223.d0*a4)*F1155d8192
    c4847=-(191.d0*b2+89.d0*a2)*F15015d256
    c4848=F1576575d4352
    c4480=(917280.d0*aGpbG-13217568.d0*aCpbC*a2b2-4006611.d0*a8pb8*a4b4 &
        +59232420.d0*a4pb4*a6b6+103973870.d0*a8b8)*F55d139264
    c4481=-(278096.d0*aCpbC-4944471.d0*a8pb8*a2b2+7988148.d0*a4pb4*a4b4 &
        +14890422.d0*a6b6)*a2pb2*F165d4096
    c4482=(1013397.d0*aCpbC-22236032.d0*a8pb8*a2b2+43190147.d0*a4pb4*a4b4 &
        +109817136.d0*a6b6)*F495d8192
    c4483=-(22595.d0*a8pb8-1471460.d0*a4pb4*a2b2+8542898.d0*a4b4) &
        *a2pb2*F6435d4096
    c4484=-(141085.d0*a8pb8-843420.d0*a4pb4*a2b2-13500018.d0*a4b4) &
        *F7865d8192
    c4485=(59.d0*a4pb4-3594.d0*a2b2)*a2pb2*F96525d128
    c4486=(669.d0*a4pb4+5050.d0*a2b2)*F17325d128
    c4487=-a2pb2*F225225d8
    c4488=F675675d544
    c4660=(655200.d0*aG-14851200.d0*aEb2+39118275.d0*aCb4 &
        +71653725.d0*aAb6-51986935.d0*a8b8-130886145.d0*a6bA &
        -35111664.d0*a4bC+28068768.d0*a2bE-1572480.d0*bG)*F77d104448
    c4661=-(291200.d0*aE-7635250.d0*aCb2+29777475.d0*aAb4 &
        +13525655.d0*a8b6-59282795.d0*a6b8-35864829.d0*a4bA &
        +16968000.d0*a2bC-847392.d0*bE)*F77d2048
    c4662=(255675.d0*aC-9925825.d0*aAb2+64059996.d0*a8b4 &
        -54908568.d0*a6b6-107250143.d0*a4b8+32161857.d0*a2bA &
        -1269072.d0*bC)*F231d2048
    c4663=(108075.d0*aA-1040435.d0*a8b2-12671208.d0*a6b4+55099836.d0*a4b6 &
        -7652755.d0*a2b8+27495.d0*bA)*F1001d2048
    c4664=-(21385.d0*a8-1243695.d0*a6b2+6750009.d0*a4b4+2087115.d0*a2b6 &
        -162470.d0*b8)*F11011d6144
    c4665=-(21935.d0*a6-306537.d0*a4b2-824663.d0*a2b4-3055.d0*b6) &
        *F9009d2048
    c4666=-(2049.d0*a4+50500.d0*a2b2+11331.d0*b4)*F1617d128
    c4667=(63.d0*a2+97.d0*b2)*F21021d64
    c4668=-F315315d272
    c6460=(655200.d0*bG-14851200.d0*a2bE+39118275.d0*a4bC &
        +71653725.d0*a6bA-51986935.d0*a8b8-130886145.d0*aAb6 &
        -35111664.d0*aCb4+28068768.d0*aEb2-1572480.d0*aG)*F77d104448
    c6461=-(291200.d0*bE-7635250.d0*a2bC+29777475.d0*a4bA &
        +13525655.d0*a6b8-59282795.d0*a8b6-35864829.d0*aAb4 &
        +16968000.d0*aCb2-847392.d0*aE)*F77d2048
    c6462=(255675.d0*bC-9925825.d0*a2bA+64059996.d0*a4b8 &
        -54908568.d0*a6b6-107250143.d0*a8b4+32161857.d0*aAb2 &
        -1269072.d0*aC)*F231d2048
    c6463=(108075.d0*bA-1040435.d0*a2b8-12671208.d0*a4b6+55099836.d0*a6b4 &
        -7652755.d0*a8b2+27495.d0*aA)*F1001d2048
    c6464=-(21385.d0*b8-1243695.d0*a2b6+6750009.d0*a4b4+2087115.d0*a6b2 &
        -162470.d0*a8)*F11011d6144
    c6465=-(21935.d0*b6-306537.d0*a2b4-824663.d0*a4b2-3055.d0*a6) &
        *F9009d2048
    c6466=-(2049.d0*b4+50500.d0*a2b2+11331.d0*a4)*F1617d128
    c6467=(63.d0*b2+97.d0*a2)*F21021d64
    c6468=-F315315d272
    c6640=-(1572480.d0*aGpbG-43216992.d0*aCpbC*a2b2 &
        +194168016.d0*a8pb8*a4b4-18229185.d0*a4pb4*a6b6 &
        -394964570.d0*a8b8)*F77d104448
    c6641=(550368.d0*aCpbC-17518368.d0*a8pb8*a2b2 &
        +114003939.d0*a4pb4*a4b4-213489754.d0*a6b6)*a2pb2*F77d2048
    c6642=-(229488.d0*aCpbC-11954943.d0*a8pb8*a2b2 &
        +107250143.d0*a4pb4*a4b4-238765956.d0*a6b6)*F231d2048
    c6643=-(197415.d0*a8pb8-4757630.d0*a4pb4*a2b2+17428838.d0*a4b4) &
        *a2pb2*F1001d2048
    c6644=-(21385.d0*a8pb8+283755.d0*a4pb4*a2b2-4031748.d0*a4b4) &
        *F11011d6144
    c6645=(1441.d0*a4pb4-31982.d0*a2b2)*a2pb2*F75075d2048
    c6646=(4383.d0*a4pb4+43630.d0*a2b2)*F13475d2048
    c6647=-a2pb2*F175175d16
    c6648=F525525d1088
endif
!
return;end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fcprep(n,c)
!
! Computation of expansion coefficients, C_{ijk},from C_{ijkn} and c
! of a uniform rectangular prism
! by even order 0, 2, to 16 truncated Taylor expansion
!
! Reference: Fukushima, T, (2019) Geophys. J. Int'l, revised
!   Taylor expansion of prismatic gravitational field
!
! Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
!
! Caution: one must call "fabset" to prepare C_{ijkn} beforehand
!
implicit real*8 (a-z)
integer n
!
common /cset0/gm
common /cset2/c200,c020,c002
common /cset4/c400,c040,c004,c022,c202,c220
common /cset6/c600,c060,c006,c420,c042,c204,c402,c024,c240,c222
common /cset8/c800,c080,c008,c620,c062,c206,c602,c026,c260, &
              c044,c404,c440,c422,c242,c224
common /csetA/cA00,c0A0,c00A,c820,c082,c208,c802,c028,c280, &
              c640,c064,c406,c604,c046,c460,c622,c262,c226, &
              c244,c424,c442
common /csetC/cC00,c0C0,c00C,cA20,c0A2,c20A,cA02,c02A,c2A0, &
              c840,c084,c408,c804,c048,c480,c822,c282,c228, &
              c066,c606,c660,c642,c264,c426,c624,c246,c462,c444
common /csetE/cE00,c0E0,c00E,cC20,c0C2,c20C,cC02,c02C,c2C0, &
              cA40,c0A4,c40A,cA04,c04A,c4A0,cA22,c2A2,c22A, &
              c860,c086,c608,c806,c068,c680,c842,c284,c428, &
              c824,c248,c482,c266,c626,c662,c644,c464,c446
common /csetG/cG00,c0G0,c00G,cE20,c0E2,c20E,cE02,c02E,c2E0, &
              cC40,c0C4,c40C,cc04,c04C,c4c0,cC22,c2C2,c22C, &
              cA60,c0A6,c60A,cA06,c06A,c6A0,cA42,c2A4,c42A, &
              cA24,c24A,c4A2,c088,c808,c880,c862,c286,c628, &
              c826,c268,c682,c844,c484,c448,c466,c646,c664
!
common /fcsetgm/gm1
common /fcset200/c2000,c2001
common /fcset020/c0200,c0201
common /fcset002/c0020,c0021
common /fcset400/c4000,c4001,c4002
common /fcset040/c0400,c0401,c0402
common /fcset004/c0040,c0041,c0042
common /fcset022/c0220,c0221,c0222
common /fcset202/c2020,c2021,c2022
common /fcset220/c2200,c2201,c2202
common /fcset600/c6000,c6001,c6002,c6003
common /fcset060/c0600,c0601,c0602,c0603
common /fcset006/c0060,c0061,c0062,c0063
common /fcset420/c4200,c4201,c4202,c4203
common /fcset042/c0420,c0421,c0422,c0423
common /fcset204/c2040,c2041,c2042,c2043
common /fcset402/c4020,c4021,c4022,c4023
common /fcset240/c2400,c2401,c2402,c2403
common /fcset024/c0240,c0241,c0242,c0243
common /fcset222/c2220,c2221,c2222,c2223
common /fcset800/c8000,c8001,c8002,c8003,c8004
common /fcset080/c0800,c0801,c0802,c0803,c0804
common /fcset008/c0080,c0081,c0082,c0083,c0084
common /fcset620/c6200,c6201,c6202,c6203,c6204
common /fcset062/c0620,c0621,c0622,c0623,c0624
common /fcset206/c2060,c2061,c2062,c2063,c2064
common /fcset602/c6020,c6021,c6022,c6023,c6024
common /fcset026/c0260,c0261,c0262,c0263,c0264
common /fcset260/c2600,c2601,c2602,c2603,c2604
common /fcset044/c0440,c0441,c0442,c0443,c0444
common /fcset404/c4040,c4041,c4042,c4043,c4044
common /fcset440/c4400,c4401,c4402,c4403,c4404
common /fcset422/c4220,c4221,c4222,c4223,c4224
common /fcset242/c2420,c2421,c2422,c2423,c2424
common /fcset224/c2240,c2241,c2242,c2243,c2244
common /fcsetA00/cA000,cA001,cA002,cA003,cA004,cA005
common /fcset0A0/c0A00,c0A01,c0A02,c0A03,c0A04,c0A05
common /fcset00A/c00A0,c00A1,c00A2,c00A3,c00A4,c00A5
common /fcset820/c8200,c8201,c8202,c8203,c8204,c8205
common /fcset082/c0820,c0821,c0822,c0823,c0824,c0825
common /fcset208/c2080,c2081,c2082,c2083,c2084,c2085
common /fcset802/c8020,c8021,c8022,c8023,c8024,c8025
common /fcset028/c0280,c0281,c0282,c0283,c0284,c0285
common /fcset280/c2800,c2801,c2802,c2803,c2804,c2805
common /fcset640/c6400,c6401,c6402,c6403,c6404,c6405
common /fcset064/c0640,c0641,c0642,c0643,c0644,c0645
common /fcset406/c4060,c4061,c4062,c4063,c4064,c4065
common /fcset604/c6040,c6041,c6042,c6043,c6044,c6045
common /fcset046/c0460,c0461,c0462,c0463,c0464,c0465
common /fcset460/c4600,c4601,c4602,c4603,c4604,c4605
common /fcset622/c6220,c6221,c6222,c6223,c6224,c6225
common /fcset262/c2620,c2621,c2622,c2623,c2624,c2625
common /fcset226/c2260,c2261,c2262,c2263,c2264,c2265
common /fcset244/c2440,c2441,c2442,c2443,c2444,c2445
common /fcset424/c4240,c4241,c4242,c4243,c4244,c4245
common /fcset442/c4420,c4421,c4422,c4423,c4424,c4425
common /fcsetC00/cC000,cC001,cC002,cC003,cC004,cC005,cC006
common /fcset0C0/c0C00,c0C01,c0C02,c0C03,c0C04,c0C05,c0C06
common /fcset00C/c00C0,c00C1,c00C2,c00C3,c00C4,c00C5,c00C6
common /fcsetA20/cA200,cA201,cA202,cA203,cA204,cA205,cA206
common /fcset0A2/c0A20,c0A21,c0A22,c0A23,c0A24,c0A25,c0A26
common /fcset20A/c20A0,c20A1,c20A2,c20A3,c20A4,c20A5,c20A6
common /fcsetA02/cA020,cA021,cA022,cA023,cA024,cA025,cA026
common /fcset2A0/c2A00,c2A01,c2A02,c2A03,c2A04,c2A05,c2A06
common /fcset02A/c02A0,c02A1,c02A2,c02A3,c02A4,c02A5,c02A6
common /fcset840/c8400,c8401,c8402,c8403,c8404,c8405,c8406
common /fcset084/c0840,c0841,c0842,c0843,c0844,c0845,c0846
common /fcset408/c4080,c4081,c4082,c4083,c4084,c4085,c4086
common /fcset804/c8040,c8041,c8042,c8043,c8044,c8045,c8046
common /fcset048/c0480,c0481,c0482,c0483,c0484,c0485,c0486
common /fcset480/c4800,c4801,c4802,c4803,c4804,c4805,c4806
common /fcset822/c8220,c8221,c8222,c8223,c8224,c8225,c8226
common /fcset282/c2820,c2821,c2822,c2823,c2824,c2825,c2826
common /fcset228/c2280,c2281,c2282,c2283,c2284,c2285,c2286
common /fcset066/c0660,c0661,c0662,c0663,c0664,c0665,c0666
common /fcset606/c6060,c6061,c6062,c6063,c6064,c6065,c6066
common /fcset660/c6600,c6601,c6602,c6603,c6604,c6605,c6606
common /fcset642/c6420,c6421,c6422,c6423,c6424,c6425,c6426
common /fcset264/c2640,c2641,c2642,c2643,c2644,c2645,c2646
common /fcset426/c4260,c4261,c4262,c4263,c4264,c4265,c4266
common /fcset624/c6240,c6241,c6242,c6243,c6244,c6245,c6246
common /fcset462/c4620,c4621,c4622,c4623,c4624,c4625,c4626
common /fcset246/c2460,c2461,c2462,c2463,c2464,c2465,c2466
common /fcset444/c4440,c4441,c4442,c4443,c4444,c4445,c4446
common /fcsetE00/cE000,cE001,cE002,cE003,cE004,cE005,cE006,cE007
common /fcset0E0/c0E00,c0E01,c0E02,c0E03,c0E04,c0E05,c0E06,c0E07
common /fcset00E/c00E0,c00E1,c00E2,c00E3,c00E4,c00E5,c00E6,c00E7
common /fcsetC20/cC200,cC201,cC202,cC203,cC204,cC205,cC206,cC207
common /fcset0C2/c0C20,c0C21,c0C22,c0C23,c0C24,c0C25,c0C26,c0C27
common /fcset20C/c20C0,c20C1,c20C2,c20C3,c20C4,c20C5,c20C6,c20C7
common /fcsetC02/cC020,cC021,cC022,cC023,cC024,cC025,cC026,cC027
common /fcset2C0/c2C00,c2C01,c2C02,c2C03,c2C04,c2C05,c2C06,c2C07
common /fcset02C/c02C0,c02C1,c02C2,c02C3,c02C4,c02C5,c02C6,c02C7
common /fcsetA40/cA400,cA401,cA402,cA403,cA404,cA405,cA406,cA407
common /fcset0A4/c0A40,c0A41,c0A42,c0A43,c0A44,c0A45,c0A46,c0A47
common /fcset40A/c40A0,c40A1,c40A2,c40A3,c40A4,c40A5,c40A6,c40A7
common /fcsetA04/cA040,cA041,cA042,cA043,cA044,cA045,cA046,cA047
common /fcset4A0/c4A00,c4A01,c4A02,c4A03,c4A04,c4A05,c4A06,c4A07
common /fcset04A/c04A0,c04A1,c04A2,c04A3,c04A4,c04A5,c04A6,c04A7
common /fcsetA22/cA220,cA221,cA222,cA223,cA224,cA225,cA226,cA227
common /fcset2A2/c2A20,c2A21,c2A22,c2A23,c2A24,c2A25,c2A26,c2A27
common /fcset22A/c22A0,c22A1,c22A2,c22A3,c22A4,c22A5,c22A6,c22A7
common /fcset860/c8600,c8601,c8602,c8603,c8604,c8605,c8606,c8607
common /fcset086/c0860,c0861,c0862,c0863,c0864,c0865,c0866,c0867
common /fcset608/c6080,c6081,c6082,c6083,c6084,c6085,c6086,c6087
common /fcset806/c8060,c8061,c8062,c8063,c8064,c8065,c8066,c8067
common /fcset680/c6800,c6801,c6802,c6803,c6804,c6805,c6806,c6807
common /fcset068/c0680,c0681,c0682,c0683,c0684,c0685,c0686,c0687
common /fcset842/c8420,c8421,c8422,c8423,c8424,c8425,c8426,c8427
common /fcset284/c2840,c2841,c2842,c2843,c2844,c2845,c2846,c2847
common /fcset428/c4280,c4281,c4282,c4283,c4284,c4285,c4286,c4287
common /fcset824/c8240,c8241,c8242,c8243,c8244,c8245,c8246,c8247
common /fcset482/c4820,c4821,c4822,c4823,c4824,c4825,c4826,c4827
common /fcset248/c2480,c2481,c2482,c2483,c2484,c2485,c2486,c2487
common /fcset266/c2660,c2661,c2662,c2663,c2664,c2665,c2666,c2667
common /fcset626/c6260,c6261,c6262,c6263,c6264,c6265,c6266,c6267
common /fcset662/c6620,c6621,c6622,c6623,c6624,c6625,c6626,c6627
common /fcset644/c6440,c6441,c6442,c6443,c6444,c6445,c6446,c6447
common /fcset464/c4640,c4641,c4642,c4643,c4644,c4645,c4646,c4647
common /fcset446/c4460,c4461,c4462,c4463,c4464,c4465,c4466,c4467
common /fcsetG00/cG000,cG001,cG002,cG003,cG004,cG005,cG006,cG007,cG008
common /fcset0G0/c0G00,c0G01,c0G02,c0G03,c0G04,c0G05,c0G06,c0G07,c0G08
common /fcset00G/c00G0,c00G1,c00G2,c00G3,c00G4,c00G5,c00G6,c00G7,c00G8
common /fcsetE20/cE200,cE201,cE202,cE203,cE204,cE205,cE206,cE207,cE208
common /fcset0E2/c0E20,c0E21,c0E22,c0E23,c0E24,c0E25,c0E26,c0E27,c0E28
common /fcset20E/c20E0,c20E1,c20E2,c20E3,c20E4,c20E5,c20E6,c20E7,c20E8
common /fcsetE02/cE020,cE021,cE022,cE023,cE024,cE025,cE026,cE027,cE028
common /fcset2E0/c2E00,c2E01,c2E02,c2E03,c2E04,c2E05,c2E06,c2E07,c2E08
common /fcset02E/c02E0,c02E1,c02E2,c02E3,c02E4,c02E5,c02E6,c02E7,c02E8
common /fcsetC40/cC400,cC401,cC402,cC403,cC404,cC405,cC406,cC407,cC408
common /fcset0C4/c0C40,c0C41,c0C42,c0C43,c0C44,c0C45,c0C46,c0C47,c0C48
common /fcset40C/c40C0,c40C1,c40C2,c40C3,c40C4,c40C5,c40C6,c40C7,c40C8
common /fcsetC04/cC040,cC041,cC042,cC043,cC044,cC045,cC046,cC047,cC048
common /fcset4C0/c4C00,c4C01,c4C02,c4C03,c4C04,c4C05,c4C06,c4C07,c4C08
common /fcset04C/c04C0,c04C1,c04C2,c04C3,c04C4,c04C5,c04C6,c04C7,c04C8
common /fcsetC22/cC220,cC221,cC222,cC223,cC224,cC225,cC226,cC227,cC228
common /fcset2C2/c2C20,c2C21,c2C22,c2C23,c2C24,c2C25,c2C26,c2C27,c2C28
common /fcset22C/c22C0,c22C1,c22C2,c22C3,c22C4,c22C5,c22C6,c22C7,c22C8
common /fcsetA60/cA600,cA601,cA602,cA603,cA604,cA605,cA606,cA607,cA608
common /fcset0A6/c0A60,c0A61,c0A62,c0A63,c0A64,c0A65,c0A66,c0A67,c0A68
common /fcset60A/c60A0,c60A1,c60A2,c60A3,c60A4,c60A5,c60A6,c60A7,c60A8
common /fcsetA06/cA060,cA061,cA062,cA063,cA064,cA065,cA066,cA067,cA068
common /fcset6A0/c6A00,c6A01,c6A02,c6A03,c6A04,c6A05,c6A06,c6A07,c6A08
common /fcset06A/c06A0,c06A1,c06A2,c06A3,c06A4,c06A5,c06A6,c06A7,c06A8
common /fcsetA42/cA420,cA421,cA422,cA423,cA424,cA425,cA426,cA427,cA428
common /fcset2A4/c2A40,c2A41,c2A42,c2A43,c2A44,c2A45,c2A46,c2A47,c2A48
common /fcset42A/c42A0,c42A1,c42A2,c42A3,c42A4,c42A5,c42A6,c42A7,c42A8
common /fcsetA24/cA240,cA241,cA242,cA243,cA244,cA245,cA246,cA247,cA248
common /fcset4A2/c4A20,c4A21,c4A22,c4A23,c4A24,c4A25,c4A26,c4A27,c4A28
common /fcset24A/c24A0,c24A1,c24A2,c24A3,c24A4,c24A5,c24A6,c24A7,c24A8
common /fcset088/c0880,c0881,c0882,c0883,c0884,c0885,c0886,c0887,c0888
common /fcset808/c8080,c8081,c8082,c8083,c8084,c8085,c8086,c8087,c8088
common /fcset880/c8800,c8801,c8802,c8803,c8804,c8805,c8806,c8807,c8808
common /fcset862/c8620,c8621,c8622,c8623,c8624,c8625,c8626,c8627,c8628
common /fcset286/c2860,c2861,c2862,c2863,c2864,c2865,c2866,c2867,c2868
common /fcset628/c6280,c6281,c6282,c6283,c6284,c6285,c6286,c6287,c6288
common /fcset826/c8260,c8261,c8262,c8263,c8264,c8265,c8266,c8267,c8268
common /fcset682/c6820,c6821,c6822,c6823,c6824,c6825,c6826,c6827,c6828
common /fcset268/c2680,c2681,c2682,c2683,c2684,c2685,c2686,c2687,c2688
common /fcset844/c8440,c8441,c8442,c8443,c8444,c8445,c8446,c8447,c8448
common /fcset484/c4840,c4841,c4842,c4843,c4844,c4845,c4846,c4847,c4848
common /fcset448/c4480,c4481,c4482,c4483,c4484,c4485,c4486,c4487,c4488
common /fcset466/c4660,c4661,c4662,c4663,c4664,c4665,c4666,c4667,c4668
common /fcset646/c6460,c6461,c6462,c6463,c6464,c6465,c6466,c6467,c6468
common /fcset664/c6640,c6641,c6642,c6643,c6644,c6645,c6646,c6647,c6648
!
if(n.ge.0) then
    gm=gm1*c
endif
if(n.ge.2) then
    c2=c*c
    c200=c2000+c2*c2001
    c020=c0200+c2*c0201
    c002=c0020+c2*c0021
endif
if(n.ge.4) then
    c400=c4000+c2*(c4001+c2*c4002)
    c040=c0400+c2*(c0401+c2*c0402)
    c004=c0040+c2*(c0041+c2*c0042)
    c022=c0220+c2*(c0221+c2*c0222)
    c202=c2020+c2*(c2021+c2*c2022)
    c220=c2200+c2*(c2201+c2*c2202)
endif
if(n.ge.6) then
    c600=c6000+c2*(c6001+c2*(c6002+c2*c6003))
    c060=c0600+c2*(c0601+c2*(c0602+c2*c0603))
    c006=c0060+c2*(c0061+c2*(c0062+c2*c0063))
    c420=c4200+c2*(c4201+c2*(c4202+c2*c4203))
    c042=c0420+c2*(c0421+c2*(c0422+c2*c0423))
    c204=c2040+c2*(c2041+c2*(c2042+c2*c2043))
    c240=c2400+c2*(c2401+c2*(c2402+c2*c2403))
    c024=c0240+c2*(c0241+c2*(c0242+c2*c0243))
    c402=c4020+c2*(c4021+c2*(c4022+c2*c4023))
    c222=c2220+c2*(c2221+c2*(c2222+c2*c2223))
endif
if(n.ge.8) then
    c800=c8000+c2*(c8001+c2*(c8002+c2*(c8003+c2*c8004)))
    c080=c0800+c2*(c0801+c2*(c0802+c2*(c0803+c2*c0804)))
    c008=c0080+c2*(c0081+c2*(c0082+c2*(c0083+c2*c0084)))
    c620=c6200+c2*(c6201+c2*(c6202+c2*(c6203+c2*c6204)))
    c062=c0620+c2*(c0621+c2*(c0622+c2*(c0623+c2*c0624)))
    c206=c2060+c2*(c2061+c2*(c2062+c2*(c2063+c2*c2064)))
    c260=c2600+c2*(c2601+c2*(c2602+c2*(c2603+c2*c2604)))
    c026=c0260+c2*(c0261+c2*(c0262+c2*(c0263+c2*c0264)))
    c602=c6020+c2*(c6021+c2*(c6022+c2*(c6023+c2*c6024)))
    c044=c0440+c2*(c0441+c2*(c0442+c2*(c0443+c2*c0444)))
    c404=c4040+c2*(c4041+c2*(c4042+c2*(c4043+c2*c4044)))
    c440=c4400+c2*(c4401+c2*(c4402+c2*(c4403+c2*c4404)))
    c422=c4220+c2*(c4221+c2*(c4222+c2*(c4223+c2*c4224)))
    c242=c2420+c2*(c2421+c2*(c2422+c2*(c2423+c2*c2424)))
    c224=c2240+c2*(c2241+c2*(c2242+c2*(c2243+c2*c2244)))
endif
if(n.ge.10) then
    cA00=cA000+c2*(cA001+c2*(cA002+c2*(cA003+c2*(cA004+c2*cA005))))
    c0A0=c0A00+c2*(c0A01+c2*(c0A02+c2*(c0A03+c2*(c0A04+c2*c0A05))))
    c00A=c00A0+c2*(c00A1+c2*(c00A2+c2*(c00A3+c2*(c00A4+c2*c00A5))))
    c820=c8200+c2*(c8201+c2*(c8202+c2*(c8203+c2*(c8204+c2*c8205))))
    c082=c0820+c2*(c0821+c2*(c0822+c2*(c0823+c2*(c0824+c2*c0825))))
    c208=c2080+c2*(c2081+c2*(c2082+c2*(c2083+c2*(c2084+c2*c2085))))
    c280=c2800+c2*(c2801+c2*(c2802+c2*(c2803+c2*(c2804+c2*c2805))))
    c802=c8020+c2*(c8021+c2*(c8022+c2*(c8023+c2*(c8024+c2*c8025))))
    c028=c0280+c2*(c0281+c2*(c0282+c2*(c0283+c2*(c0284+c2*c0285))))
    c640=c6400+c2*(c6401+c2*(c6402+c2*(c6403+c2*(c6404+c2*c6405))))
    c064=c0640+c2*(c0641+c2*(c0642+c2*(c0643+c2*(c0644+c2*c0645))))
    c406=c4060+c2*(c4061+c2*(c4062+c2*(c4063+c2*(c4064+c2*c4065))))
    c460=c4600+c2*(c4601+c2*(c4602+c2*(c4603+c2*(c4604+c2*c4605))))
    c604=c6040+c2*(c6041+c2*(c6042+c2*(c6043+c2*(c6044+c2*c6045))))
    c046=c0460+c2*(c0461+c2*(c0462+c2*(c0463+c2*(c0464+c2*c0465))))
    c622=c6220+c2*(c6221+c2*(c6222+c2*(c6223+c2*(c6224+c2*c6225))))
    c262=c2620+c2*(c2621+c2*(c2622+c2*(c2623+c2*(c2624+c2*c2625))))
    c226=c2260+c2*(c2261+c2*(c2262+c2*(c2263+c2*(c2264+c2*c2265))))
    c244=c2440+c2*(c2441+c2*(c2442+c2*(c2443+c2*(c2444+c2*c2445))))
    c424=c4240+c2*(c4241+c2*(c4242+c2*(c4243+c2*(c4244+c2*c4245))))
    c442=c4420+c2*(c4421+c2*(c4422+c2*(c4423+c2*(c4424+c2*c4425))))
endif
if(n.ge.12) then
    cC00=cC000+c2*(cC001+c2*(cC002+c2*(cC003+c2*(cC004+c2*(cC005+c2*cC006)))))
    c0C0=c0C00+c2*(c0C01+c2*(c0C02+c2*(c0C03+c2*(c0C04+c2*(c0C05+c2*c0C06)))))
    c00C=c00C0+c2*(c00C1+c2*(c00C2+c2*(c00C3+c2*(c00C4+c2*(c00C5+c2*c00C6)))))
    cA20=cA200+c2*(cA201+c2*(cA202+c2*(cA203+c2*(cA204+c2*(cA205+c2*cA206)))))
    c0A2=c0A20+c2*(c0A21+c2*(c0A22+c2*(c0A23+c2*(c0A24+c2*(c0A25+c2*c0A26)))))
    c20A=c20A0+c2*(c20A1+c2*(c20A2+c2*(c20A3+c2*(c20A4+c2*(c20A5+c2*c20A6)))))
    c2A0=c2A00+c2*(c2A01+c2*(c2A02+c2*(c2A03+c2*(c2A04+c2*(c2A05+c2*c2A06)))))
    cA02=cA020+c2*(cA021+c2*(cA022+c2*(cA023+c2*(cA024+c2*(cA025+c2*cA026)))))
    c02A=c02A0+c2*(c02A1+c2*(c02A2+c2*(c02A3+c2*(c02A4+c2*(c02A5+c2*c02A6)))))
    c840=c8400+c2*(c8401+c2*(c8402+c2*(c8403+c2*(c8404+c2*(c8405+c2*c8406)))))
    c084=c0840+c2*(c0841+c2*(c0842+c2*(c0843+c2*(c0844+c2*(c0845+c2*c0846)))))
    c408=c4080+c2*(c4081+c2*(c4082+c2*(c4083+c2*(c4084+c2*(c4085+c2*c4086)))))
    c480=c4800+c2*(c4801+c2*(c4802+c2*(c4803+c2*(c4804+c2*(c4805+c2*c4806)))))
    c804=c8040+c2*(c8041+c2*(c8042+c2*(c8043+c2*(c8044+c2*(c8045+c2*c8046)))))
    c048=c0480+c2*(c0481+c2*(c0482+c2*(c0483+c2*(c0484+c2*(c0485+c2*c0486)))))
    c822=c8220+c2*(c8221+c2*(c8222+c2*(c8223+c2*(c8224+c2*(c8225+c2*c8226)))))
    c282=c2820+c2*(c2821+c2*(c2822+c2*(c2823+c2*(c2824+c2*(c2825+c2*c2826)))))
    c228=c2280+c2*(c2281+c2*(c2282+c2*(c2283+c2*(c2284+c2*(c2285+c2*c2286)))))
    c066=c0660+c2*(c0661+c2*(c0662+c2*(c0663+c2*(c0664+c2*(c0665+c2*c0666)))))
    c606=c6060+c2*(c6061+c2*(c6062+c2*(c6063+c2*(c6064+c2*(c6065+c2*c6066)))))
    c660=c6600+c2*(c6601+c2*(c6602+c2*(c6603+c2*(c6604+c2*(c6605+c2*c6606)))))
    c642=c6420+c2*(c6421+c2*(c6422+c2*(c6423+c2*(c6424+c2*(c6425+c2*c6426)))))
    c264=c2640+c2*(c2641+c2*(c2642+c2*(c2643+c2*(c2644+c2*(c2645+c2*c2646)))))
    c426=c4260+c2*(c4261+c2*(c4262+c2*(c4263+c2*(c4264+c2*(c4265+c2*c4266)))))
    c462=c4620+c2*(c4621+c2*(c4622+c2*(c4623+c2*(c4624+c2*(c4625+c2*c4626)))))
    c624=c6240+c2*(c6241+c2*(c6242+c2*(c6243+c2*(c6244+c2*(c6245+c2*c6246)))))
    c246=c2460+c2*(c2461+c2*(c2462+c2*(c2463+c2*(c2464+c2*(c2465+c2*c2466)))))
    c444=c4440+c2*(c4441+c2*(c4442+c2*(c4443+c2*(c4444+c2*(c4445+c2*c4446)))))
endif
if(n.ge.14) then
    cE00=cE000+c2*(cE001+c2*(cE002+c2*(cE003+c2*(cE004+c2*(cE005+c2*(cE006 &
        +c2*cE007))))))
    c0E0=c0E00+c2*(c0E01+c2*(c0E02+c2*(c0E03+c2*(c0E04+c2*(c0E05+c2*(c0E06 &
        +c2*c0E07))))))
    c00E=c00E0+c2*(c00E1+c2*(c00E2+c2*(c00E3+c2*(c00E4+c2*(c00E5+c2*(c00E6 &
        +c2*c00E7))))))
    cC20=cC200+c2*(cC201+c2*(cC202+c2*(cC203+c2*(cC204+c2*(cC205+c2*(cC206 &
        +c2*cC207))))))
    c0C2=c0C20+c2*(c0C21+c2*(c0C22+c2*(c0C23+c2*(c0C24+c2*(c0C25+c2*(c0C26 &
        +c2*c0C27))))))
    c20C=c20C0+c2*(c20C1+c2*(c20C2+c2*(c20C3+c2*(c20C4+c2*(c20C5+c2*(c20C6 &
        +c2*c20C7))))))
    c2C0=c2C00+c2*(c2C01+c2*(c2C02+c2*(c2C03+c2*(c2C04+c2*(c2C05+c2*(c2C06 &
        +c2*c2C07))))))
    cC02=cC020+c2*(cC021+c2*(cC022+c2*(cC023+c2*(cC024+c2*(cC025+c2*(cC026 &
        +c2*cC027))))))
    c02C=c02C0+c2*(c02C1+c2*(c02C2+c2*(c02C3+c2*(c02C4+c2*(c02C5+c2*(c02C6 &
        +c2*c02C7))))))
    cA40=cA400+c2*(cA401+c2*(cA402+c2*(cA403+c2*(cA404+c2*(cA405+c2*(cA406 &
        +c2*cA407))))))
    c0A4=c0A40+c2*(c0A41+c2*(c0A42+c2*(c0A43+c2*(c0A44+c2*(c0A45+c2*(c0A46 &
        +c2*c0A47))))))
    c40A=c40A0+c2*(c40A1+c2*(c40A2+c2*(c40A3+c2*(c40A4+c2*(c40A5+c2*(c40A6 &
        +c2*c40A7))))))
    c4A0=c4A00+c2*(c4A01+c2*(c4A02+c2*(c4A03+c2*(c4A04+c2*(c4A05+c2*(c4A06 &
        +c2*c4A07))))))
    cA04=cA040+c2*(cA041+c2*(cA042+c2*(cA043+c2*(cA044+c2*(cA045+c2*(cA046 &
        +c2*cA047))))))
    c04A=c04A0+c2*(c04A1+c2*(c04A2+c2*(c04A3+c2*(c04A4+c2*(c04A5+c2*(c04A6 &
        +c2*c04A7))))))
    cA22=cA220+c2*(cA221+c2*(cA222+c2*(cA223+c2*(cA224+c2*(cA225+c2*(cA226 &
        +c2*cA227))))))
    c2A2=c2A20+c2*(c2A21+c2*(c2A22+c2*(c2A23+c2*(c2A24+c2*(c2A25+c2*(c2A26 &
        +c2*c2A27))))))
    c22A=c22A0+c2*(c22A1+c2*(c22A2+c2*(c22A3+c2*(c22A4+c2*(c22A5+c2*(c22A6 &
        +c2*c22A7))))))
    c860=c8600+c2*(c8601+c2*(c8602+c2*(c8603+c2*(c8604+c2*(c8605+c2*(c8606 &
        +c2*c8607))))))
    c086=c0860+c2*(c0861+c2*(c0862+c2*(c0863+c2*(c0864+c2*(c0865+c2*(c0866 &
        +c2*c0867))))))
    c608=c6080+c2*(c6081+c2*(c6082+c2*(c6083+c2*(c6084+c2*(c6085+c2*(c6086 &
        +c2*c6087))))))
    c680=c6800+c2*(c6801+c2*(c6802+c2*(c6803+c2*(c6804+c2*(c6805+c2*(c6806 &
        +c2*c6807))))))
    c806=c8060+c2*(c8061+c2*(c8062+c2*(c8063+c2*(c8064+c2*(c8065+c2*(c8066 &
        +c2*c8067))))))
    c068=c0680+c2*(c0681+c2*(c0682+c2*(c0683+c2*(c0684+c2*(c0685+c2*(c0686 &
        +c2*c0687))))))
    c842=c8420+c2*(c8421+c2*(c8422+c2*(c8423+c2*(c8424+c2*(c8425+c2*(c8426 &
        +c2*c8427))))))
    c284=c2840+c2*(c2841+c2*(c2842+c2*(c2843+c2*(c2844+c2*(c2845+c2*(c2846 &
        +c2*c2847))))))
    c428=c4280+c2*(c4281+c2*(c4282+c2*(c4283+c2*(c4284+c2*(c4285+c2*(c4286 &
        +c2*c4287))))))
    c482=c4820+c2*(c4821+c2*(c4822+c2*(c4823+c2*(c4824+c2*(c4825+c2*(c4826 &
        +c2*c4827))))))
    c824=c8240+c2*(c8241+c2*(c8242+c2*(c8243+c2*(c8244+c2*(c8245+c2*(c8246 &
        +c2*c8247))))))
    c248=c2480+c2*(c2481+c2*(c2482+c2*(c2483+c2*(c2484+c2*(c2485+c2*(c2486 &
        +c2*c2487))))))
    c266=c2660+c2*(c2661+c2*(c2662+c2*(c2663+c2*(c2664+c2*(c2665+c2*(c2666 &
        +c2*c2667))))))
    c626=c6260+c2*(c6261+c2*(c6262+c2*(c6263+c2*(c6264+c2*(c6265+c2*(c6266 &
        +c2*c6267))))))
    c662=c6620+c2*(c6621+c2*(c6622+c2*(c6623+c2*(c6624+c2*(c6625+c2*(c6626 &
        +c2*c6627))))))
    c644=c6440+c2*(c6441+c2*(c6442+c2*(c6443+c2*(c6444+c2*(c6445+c2*(c6446 &
        +c2*c6447))))))
    c464=c4640+c2*(c4641+c2*(c4642+c2*(c4643+c2*(c4644+c2*(c4645+c2*(c4646 &
        +c2*c4647))))))
    c446=c4460+c2*(c4461+c2*(c4462+c2*(c4463+c2*(c4464+c2*(c4465+c2*(c4466 &
        +c2*c4467))))))
endif
if(n.ge.16) then
    cG00=cG000+c2*(cG001+c2*(cG002+c2*(cG003+c2*(cG004+c2*(cG005+c2*(cG006 &
        +c2*(cG007+c2*cG008)))))))
    c0G0=c0G00+c2*(c0G01+c2*(c0G02+c2*(c0G03+c2*(c0G04+c2*(c0G05+c2*(c0G06 &
        +c2*(c0G07+c2*c0G08)))))))
    c00G=c00G0+c2*(c00G1+c2*(c00G2+c2*(c00G3+c2*(c00G4+c2*(c00G5+c2*(c00G6 &
        +c2*(c00G7+c2*c00G8)))))))
    cE20=cE200+c2*(cE201+c2*(cE202+c2*(cE203+c2*(cE204+c2*(cE205+c2*(cE206 &
        +c2*(cE207+c2*cE208)))))))
    c0E2=c0E20+c2*(c0E21+c2*(c0E22+c2*(c0E23+c2*(c0E24+c2*(c0E25+c2*(c0E26 &
        +c2*(c0E27+c2*c0E28)))))))
    c20E=c20E0+c2*(c20E1+c2*(c20E2+c2*(c20E3+c2*(c20E4+c2*(c20E5+c2*(c20E6 &
        +c2*(c20E7+c2*c20E8)))))))
    c2E0=c2E00+c2*(c2E01+c2*(c2E02+c2*(c2E03+c2*(c2E04+c2*(c2E05+c2*(c2E06 &
        +c2*(c2E07+c2*c2E08)))))))
    cE02=cE020+c2*(cE021+c2*(cE022+c2*(cE023+c2*(cE024+c2*(cE025+c2*(cE026 &
        +c2*(cE027+c2*cE028)))))))
    c02E=c02E0+c2*(c02E1+c2*(c02E2+c2*(c02E3+c2*(c02E4+c2*(c02E5+c2*(c02E6 &
        +c2*(c02E7+c2*c02E8)))))))
    cC40=cC400+c2*(cC401+c2*(cC402+c2*(cC403+c2*(cC404+c2*(cC405+c2*(cC406 &
        +c2*(cC407+c2*cC408)))))))
    c0C4=c0C40+c2*(c0C41+c2*(c0C42+c2*(c0C43+c2*(c0C44+c2*(c0C45+c2*(c0C46 &
        +c2*(c0C47+c2*c0C48)))))))
    c40C=c40C0+c2*(c40C1+c2*(c40C2+c2*(c40C3+c2*(c40C4+c2*(c40C5+c2*(c40C6 &
        +c2*(c40C7+c2*c40C8)))))))
    c4C0=c4C00+c2*(c4C01+c2*(c4C02+c2*(c4C03+c2*(c4C04+c2*(c4C05+c2*(c4C06 &
        +c2*(c4C07+c2*c4C08)))))))
    cC04=cC040+c2*(cC041+c2*(cC042+c2*(cC043+c2*(cC044+c2*(cC045+c2*(cC046 &
        +c2*(cC047+c2*cC048)))))))
    c04C=c04C0+c2*(c04C1+c2*(c04C2+c2*(c04C3+c2*(c04C4+c2*(c04C5+c2*(c04C6 &
        +c2*(c04C7+c2*c04C8)))))))
    cC22=cC220+c2*(cC221+c2*(cC222+c2*(cC223+c2*(cC224+c2*(cC225+c2*(cC226 &
        +c2*(cC227+c2*cC228)))))))
    c2C2=c2C20+c2*(c2C21+c2*(c2C22+c2*(c2C23+c2*(c2C24+c2*(c2C25+c2*(c2C26 &
        +c2*(c2C27+c2*c2C28)))))))
    c22C=c22C0+c2*(c22C1+c2*(c22C2+c2*(c22C3+c2*(c22C4+c2*(c22C5+c2*(c22C6 &
        +c2*(c22C7+c2*c22C8)))))))
    cA60=cA600+c2*(cA601+c2*(cA602+c2*(cA603+c2*(cA604+c2*(cA605+c2*(cA606 &
        +c2*(cA607+c2*cA608)))))))
    c0A6=c0A60+c2*(c0A61+c2*(c0A62+c2*(c0A63+c2*(c0A64+c2*(c0A65+c2*(c0A66 &
        +c2*(c0A67+c2*c0A68)))))))
    c60A=c60A0+c2*(c60A1+c2*(c60A2+c2*(c60A3+c2*(c60A4+c2*(c60A5+c2*(c60A6 &
        +c2*(c60A7+c2*c60A8)))))))
    c6A0=c6A00+c2*(c6A01+c2*(c6A02+c2*(c6A03+c2*(c6A04+c2*(c6A05+c2*(c6A06 &
        +c2*(c6A07+c2*cA608)))))))
    cA06=cA060+c2*(cA061+c2*(cA062+c2*(cA063+c2*(cA064+c2*(cA065+c2*(cA066 &
        +c2*(cA067+c2*cA068)))))))
    c06A=c06A0+c2*(c06A1+c2*(c06A2+c2*(c06A3+c2*(c06A4+c2*(c06A5+c2*(c06A6 &
        +c2*(c06A7+c2*c06A8)))))))
    cA42=cA420+c2*(cA421+c2*(cA422+c2*(cA423+c2*(cA424+c2*(cA425+c2*(cA426 &
        +c2*(cA427+c2*cA428)))))))
    c2A4=c2A40+c2*(c2A41+c2*(c2A42+c2*(c2A43+c2*(c2A44+c2*(c2A45+c2*(c2A46 &
        +c2*(c2A47+c2*c2A48)))))))
    c42A=c42A0+c2*(c42A1+c2*(c42A2+c2*(c42A3+c2*(c42A4+c2*(c42A5+c2*(c42A6 &
        +c2*(c42A7+c2*c42A8)))))))
    c4A2=c4A20+c2*(c4A21+c2*(c4A22+c2*(c4A23+c2*(c4A24+c2*(c4A25+c2*(c4A26 &
        +c2*(c4A27+c2*c4A28)))))))
    cA24=cA240+c2*(cA241+c2*(cA242+c2*(cA243+c2*(cA244+c2*(cA245+c2*(cA246 &
        +c2*(cA247+c2*cA248)))))))
    c24A=c24A0+c2*(c24A1+c2*(c24A2+c2*(c24A3+c2*(c24A4+c2*(c24A5+c2*(c24A6 &
        +c2*(c24A7+c2*c24A8)))))))
    c088=c0880+c2*(c0881+c2*(c0882+c2*(c0883+c2*(c0884+c2*(c0885+c2*(c0886 &
        +c2*(c0887+c2*c0888)))))))
    c808=c8080+c2*(c8081+c2*(c8082+c2*(c8083+c2*(c8084+c2*(c8085+c2*(c8086 &
        +c2*(c8087+c2*c8088)))))))
    c880=c8800+c2*(c8801+c2*(c8802+c2*(c8803+c2*(c8804+c2*(c8805+c2*(c8806 &
        +c2*(c8807+c2*c8808)))))))
    c862=c8620+c2*(c8621+c2*(c8622+c2*(c8623+c2*(c8624+c2*(c8625+c2*(c8626 &
        +c2*(c8627+c2*c8628)))))))
    c286=c2860+c2*(c2861+c2*(c2862+c2*(c2863+c2*(c2864+c2*(c2865+c2*(c2866 &
        +c2*(c2867+c2*c2868)))))))
    c628=c6280+c2*(c6281+c2*(c6282+c2*(c6283+c2*(c6284+c2*(c6285+c2*(c6286 &
        +c2*(c6287+c2*c6288)))))))
    c682=c6820+c2*(c6821+c2*(c6822+c2*(c6823+c2*(c6824+c2*(c6825+c2*(c6826 &
        +c2*(c6827+c2*c6828)))))))
    c826=c8260+c2*(c8261+c2*(c8262+c2*(c8263+c2*(c8264+c2*(c8265+c2*(c8266 &
        +c2*(c8267+c2*c8268)))))))
    c268=c2680+c2*(c2681+c2*(c2682+c2*(c2683+c2*(c2684+c2*(c2685+c2*(c2686 &
        +c2*(c2687+c2*c2688)))))))
    c844=c8440+c2*(c8441+c2*(c8442+c2*(c8443+c2*(c8444+c2*(c8445+c2*(c8446 &
        +c2*(c8447+c2*c8448)))))))
    c484=c4840+c2*(c4841+c2*(c4842+c2*(c4843+c2*(c4844+c2*(c4845+c2*(c4846 &
        +c2*(c4847+c2*c4848)))))))
    c448=c4480+c2*(c4481+c2*(c4482+c2*(c4483+c2*(c4484+c2*(c4485+c2*(c4486 &
        +c2*(c4487+c2*c4488)))))))
    c466=c4660+c2*(c4661+c2*(c4662+c2*(c4663+c2*(c4664+c2*(c4665+c2*(c4666 &
        +c2*(c4667+c2*c4668)))))))
    c646=c6460+c2*(c6461+c2*(c6462+c2*(c6463+c2*(c6464+c2*(c6465+c2*(c6466 &
        +c2*(c6467+c2*c6468)))))))
    c664=c6640+c2*(c6641+c2*(c6642+c2*(c6643+c2*(c6644+c2*(c6645+c2*(c6646 &
        +c2*(c6647+c2*c6648)))))))
endif
!
return;end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8 function tnvpri(n,x,y,z)
!
! Computation of graviational potential of a uniform rectangular prism
! by an even order truncated Taylor expansion
!
! Reference: Fukushima, T, (2019) Geophys. J. Int'l, revised
!   Taylor expansion of prismatic gravitational field
!
! Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
!
!     Caution: Coefficients, gm and cijk, must be prepared beforehand
!
!     Inputs: n = truncated order                0 <= n <= 16
!             x,y,z = relative coordinates of evaluated point
!
!     Output: tnvpri
!
implicit real*8 (a-h,o-z)
implicit integer (i-n)
!
common /cset0/gm
common /cset2/c200,c020,c002
common /cset4/c400,c040,c004,c022,c202,c220
common /cset6/c600,c060,c006,c420,c042,c204,c402,c024,c240,c222
common /cset8/c800,c080,c008,c620,c062,c206,c602,c026,c260, &
              c044,c404,c440,c422,c242,c224
common /csetA/cA00,c0A0,c00A,c820,c082,c208,c802,c028,c280, &
              c640,c064,c406,c604,c046,c460,c622,c262,c226, &
              c244,c424,c442
common /csetC/cC00,c0C0,c00C,cA20,c0A2,c20A,cA02,c02A,c2A0, &
              c840,c084,c408,c804,c048,c480,c822,c282,c228, &
              c066,c606,c660,c642,c264,c426,c624,c246,c462,c444
common /csetE/cE00,c0E0,c00E,cC20,c0C2,c20C,cC02,c02C,c2C0, &
              cA40,c0A4,c40A,cA04,c04A,c4A0,cA22,c2A2,c22A, &
              c860,c086,c608,c806,c068,c680,c842,c284,c428, &
              c824,c248,c482,c266,c626,c662,c644,c464,c446
common /csetG/cG00,c0G0,c00G,cE20,c0E2,c20E,cE02,c02E,c2E0, &
              cC40,c0C4,c40C,cC04,c04C,c4C0,cC22,c2C2,c22C, &
              cA60,c0A6,c60A,cA06,c06A,c6A0,cA42,c2A4,c42A, &
              cA24,c24A,c4A2,c088,c808,c880,c862,c286,c628, &
              c826,c268,c682,c844,c484,c448,c466,c646,c664
!
x2=x*x;y2=y*y;z2=z*z
r2=x2+y2+z2
ri2=1.d0/r2
ri=sqrt(ri2)
if(n.le.0) then
    tnvpri=gm*ri
    return
endif
!
ri4=ri2*ri2
p2=c200*x2+c020*y2+c002*z2
if(n.le.1) then
    tnvpri=gm*ri*(1.d0+ri4*p2)
    return
endif
!
x4=x2*x2;y4=y2*y2;z4=z2*z2
y2z2=y2*z2;z2x2=z2*x2;x2y2=x2*y2
p4=c400*x4+c040*y4+c004*z4+c022*y2z2+c202*z2x2+c220*x2y2
if(n.le.2) then
    tnvpri=gm*ri*(1.d0+ri4*(p2+ri4*p4))
    return
endif
!
x6=x4*x2;y6=y4*y2;z6=z4*z2
x2y2z2=x2y2*z2
p6=c600*x6+c060*y6+c006*z6 &
    +x4*(c420*y2+c402*z2) &
    +y4*(c042*z2+c240*x2) &
    +z4*(c204*x2+c024*y2) &
    +c222*x2y2z2
if(n.le.3) then
    tnvpri=gm*ri*(1.d0+ri4*(p2+ri4*(p4+ri4*p6)))
    return
endif
!
x8=x6*x2;y8=y6*y2;z8=z6*z2
p8=c800*x8+c080*y8+c008*z8 &
    +x6*(c620*y2+c602*z2) &
    +y6*(c062*z2+c260*x2) &
    +z6*(c206*x2+c026*y2) &
    +x4*(c440*y4+c422*y2z2) &
    +y4*(c044*z4+c242*z2x2) &
    +z4*(c404*x4+c224*x2y2)
if(n.le.4) then
    tnvpri=gm*ri*(1.d0+ri4*(p2+ri4*(p4+ri4*(p6+ri4*p8))))
    return
endif
!
xA=x8*x2;yA=y8*y2;zA=z8*z2
x4y2=x4*y2;y4z2=y4*z2;z4x2=z4*x2
pA=cA00*xA+c0A0*yA+c00A*zA &
    +x8*(c820*y2+c802*z2) &
    +y8*(c082*z2+c280*x2) &
    +z8*(c208*x2+c028*y2) &
    +x6*(c640*y4+c622*y2z2+c604*z4) &
    +y6*(c064*z4+c262*z2x2+c460*x4) &
    +z6*(c406*x4+c226*x2y2+c046*y4) &
    +x4*(c442*y4z2) &
    +y4*(c244*z4x2) &
    +z4*(c424*x4y2)
if(n.le.5) then
    tnvpri=gm*ri*(1.d0+ri4*(p2+ri4*(p4+ri4*(p6+ri4*(p8+ri4*pA)))))
    return
endif
!
xC=xA*x2;yC=yA*y2;zC=zA*z2
x2y4=x2*y4;y2z4=y2*z4;z2x4=z2*x4
x4y4z4=x2y2z2*x2y2z2
pC=cC00*xC+c0C0*yC+c00C*zC &
    +xA*(cA20*y2+cA02*z2) &
    +yA*(c0A2*z2+c2A0*x2) &
    +zA*(c20A*x2+c02A*y2) &
    +x8*(c840*y4+c822*y2z2+c804*z4) &
    +y8*(c084*z4+c282*z2x2+c480*x4) &
    +z8*(c408*x4+c228*x2y2+c048*y4) &
    +x6*(c660*y6+c642*y4z2+c624*y2z4) &
    +y6*(c066*z6+c264*z4x2+c462*z2x4) &
    +z6*(c606*x6+c426*x4y2+c246*x2y4) &
    +c444*x4y4z4
if(n.le.6) then
    tnvpri=gm*ri*(1.d0+ri4*(p2+ri4*(p4+ri4*(p6 &
        +ri4*(p8+ri4*(pA+ri4*pC))))))
    return
endif
!
xE=xC*x2;yE=yC*y2;zE=zC*z2
y4z4=y4*z4;z4x4=z4*x4;x4y4=x4*y4
y6z2=y6*z2;z6x2=z6*x2;x6y2=x6*y2
pE=cE00*xE+c0E0*yE+c00E*zE &
    +xC*(cC20*y2+cC02*z2) &
    +yC*(c0C2*z2+c2C0*x2) &
    +zC*(c20C*x2+c02C*y2) &
    +xA*(cA40*y4+cA22*y2z2+cA04*z4) &
    +yA*(c0A4*z4+c2A2*z2x2+c4A0*x4) &
    +zA*(c40A*x4+c22A*x2y2+c04A*y4) &
    +x8*(c860*y6+c842*y4z2+c824*y2z4+c806*z6) &
    +y8*(c086*z6+c284*z4x2+c482*z2x4+c680*x6) &
    +z8*(c608*x6+c428*x4y2+c248*x2y4+c068*y6) &
    +x6*(c662*y6z2+c644*y4z4) &
    +y6*(c266*z6x2+c464*z4x4) &
    +z6*(c626*x6y2+c446*x4y4)
if(n.le.7) then
    tnvpri=gm*ri*(1.d0+ri4*(p2+ri4*(p4+ri4*(p6 &
        +ri4*(p8+ri4*(pA+ri4*(pC+ri4*pE)))))))
    return
endif
!
xG=xE*x2;yG=yE*y2;zG=zE*z2
y2z6=y2*z6;z2x6=z2*x6;x2y6=x2*y6
pG=cG00*xG+c0G0*yG+c00G*zG &
    +xE*(cE20*y2+cE02*z2) &
    +yE*(c0E2*z2+c2E0*x2) &
    +zE*(c20E*x2+c02E*y2) &
    +xC*(cC40*y4+cC22*y2z2+cC04*z4) &
    +yC*(c0C4*z4+c2C2*z2x2+c4C0*x4) &
    +zC*(c40C*x4+c22C*x2y2+c04C*y4) &
    +xA*(cA60*y6+cA42*y4z2+cA24*y2z4+cA06*z6) &
    +yA*(c0A6*z6+c2A4*z4x2+c4A2*z2x4+c6A0*x6) &
    +zA*(c60A*x6+c42A*x4y2+c24A*x2y4+c06A*y6) &
    +x8*(c880*y8+c862*y6z2+c844*y4z4+c826*y2z6) &
    +y8*(c088*z8+c286*z6x2+c484*z4x4+c682*z2x6) &
    +z8*(c808*x8+c628*x6y2+c448*x4y4+c268*x2y6) &
    +c466*x4*y6*z6+c646*x6*y4*z6+c664*x6*y6*z4

tnvpri=gm*ri*(1.d0+ri4*(p2+ri4*(p4+ri4*(p6 &
    +ri4*(p8+ri4*(pA+ri4*(pC+ri4*(pE+ri4*pG))))))))
!
return;end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tngpri(n,x,y,z,u,ux,uy,uz)
!
! Computation of graviational potential and acceleration vector
! of a uniform rectangular prism
! by an even order truncated Taylor expansion
!
! Reference: Fukushima, T, (2019) Geophys. J. Int'l, revised
!   Taylor expansion of prismatic gravitational field
!
! Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
!
!     Caution: Coefficients, gm and cijk, must be prepared beforehand
!
!     Inputs: n = truncated order                0 <= n <= 16
!             x,y,z = relative coordinates of evaluated point
!
!     Output: u = gravitational potential 
!             ux = x-component gravitational acceleration
!             uy = y-component gravitational acceleration
!             uz = z-component gravitational acceleration
!
implicit real*8 (a-h,o-z)
implicit integer (i-n)
!
common /cset0/gm
common /cset2/c200,c020,c002
common /cset4/c400,c040,c004,c022,c202,c220
common /cset6/c600,c060,c006,c420,c042,c204,c402,c024,c240,c222
common /cset8/c800,c080,c008,c620,c062,c206,c602,c026,c260, &
              c044,c404,c440,c422,c242,c224
common /csetA/cA00,c0A0,c00A,c820,c082,c208,c802,c028,c280, &
              c640,c064,c406,c604,c046,c460,c622,c262,c226, &
              c244,c424,c442
common /csetC/cC00,c0C0,c00C,cA20,c0A2,c20A,cA02,c02A,c2A0, &
              c840,c084,c408,c804,c048,c480,c822,c282,c228, &
              c066,c606,c660,c642,c264,c426,c624,c246,c462,c444
common /csetE/cE00,c0E0,c00E,cC20,c0C2,c20C,cC02,c02C,c2C0, &
              cA40,c0A4,c40A,cA04,c04A,c4A0,cA22,c2A2,c22A, &
              c860,c086,c608,c806,c068,c680,c842,c284,c428, &
              c824,c248,c482,c266,c626,c662,c644,c464,c446
common /csetG/cG00,c0G0,c00G,cE20,c0E2,c20E,cE02,c02E,c2E0, &
              cC40,c0C4,c40C,cc04,c04C,c4c0,cC22,c2C2,c22C, &
              cA60,c0A6,c60A,cA06,c06A,c6A0,cA42,c2A4,c42A, &
              cA24,c24A,c4A2,c088,c808,c880,c862,c286,c628, &
              c826,c268,c682,c844,c484,c448,c466,c646,c664
!
x2=x*x;y2=y*y;z2=z*z
xy=x*y;yz=y*z;zx=z*x
r2=x2+y2+z2
ri2=1.d0/r2
ri=sqrt(ri2)
gmri=gm*ri
ri3=ri*ri2
gmri3=gm*ri3
if(n.le.0) then
    u=gmri
    ur=-gmri3
    ux=ur*x;uy=ur*y;uz=ur*z
    return
endif
!
ri5=ri3*ri2
gmri5=gm*ri5
ri4=ri2*ri2
p2=c200*x2+c020*y2+c002*z2
p2xx=2.d0*c200
p2yy=2.d0*c020
p2zz=2.d0*c002
p2x=p2xx*x
p2y=p2yy*y
p2z=p2zz*z
if(n.le.1) then
    u=gmri*(1.d0+ri4*p2)
    r5i4=ri4*5.d0
    ur=-gmri3*(1.d0+r5i4*p2)
    ux=ur*x+gmri5*p2x
    uy=ur*y+gmri5*p2y
    uz=ur*z+gmri5*p2z
    return
endif
!
x4=x2*x2;y4=y2*y2;z4=z2*z2
y2z2=y2*z2;z2x2=z2*x2;x2y2=x2*y2
p4=c400*x4+c040*y4+c004*z4+c022*y2z2+c202*z2x2+c220*x2y2
a4x2=4.d0*c400*x2
a4x0=2.d0*(c202*z2+c220*y2)
a4y2=4.d0*c040*y2
a4y0=2.d0*(c220*x2+c022*z2)
a4z2=4.d0*c004*z2
a4z0=2.d0*(c022*y2+c202*x2)
p4x=(a4x2+a4x0)*x
p4y=(a4y2+a4y0)*y
p4z=(a4z2+a4z0)*z
if(n.le.2) then
    u=gmri*(1.d0+ri4*(p2+ri4*p4))
    r9i4=ri4*9.d0
    ur=-gmri3*(1.d0+ri4*(5.d0*p2+r9i4*p4))
    ux=ur*x+gmri5*(p2x+ri4*p4x)
    uy=ur*y+gmri5*(p2y+ri4*p4y)
    uz=ur*z+gmri5*(p2z+ri4*p4z)
    return
endif
!
x6=x4*x2;y6=y4*y2;z6=z4*z2
x4y2=x4*y2;y4z2=y4*z2;z4x2=z4*x2
x2y4=x2*y4;y2z4=y2*z4;z2x4=z2*x4
x2y2z2=x2y2*z2
p6=c600*x6+c060*y6+c006*z6+c222*x2y2z2 &
    +c420*x4y2+c402*z2x4+c042*y4z2+c240*x2y4+c204*z4x2+c024*y2z4
a6x4=6.d0*c600*x4
a6x2=4.d0*(c420*x2y2+c402*z2x2)
a6x0=2.d0*(c240*y4+c222*y2z2+c204*z4)
a6y4=6.d0*c060*y4
a6y2=4.d0*(c042*y2z2+c240*x2y2)
a6y0=2.d0*(c024*z4+c222*z2x2+c420*x4)
a6z4=6.d0*c006*z4
a6z2=4.d0*(c204*z2x2+c024*y2z2)
a6z0=2.d0*(c402*x4+c222*x2y2+c042*y4)
p6x=(a6x4+a6x2+a6x0)*x
p6y=(a6y4+a6y2+a6y0)*y
p6z=(a6z4+a6z2+a6z0)*z
if(n.le.3) then
    u=gmri*(1.d0+ri4*(p2+ri4*(p4+ri4*p6)))
    r13i4=ri4*13.d0
    ur=-gmri3*(1.d0+ri4*(5.d0*p2+ri4*(9.d0*p4+r13i4*p6)))
    ux=ur*x+gmri5*(p2x+ri4*(p4x+ri4*p6x))
    uy=ur*y+gmri5*(p2y+ri4*(p4y+ri4*p6y))
    uz=ur*z+gmri5*(p2z+ri4*(p4z+ri4*p6z))
    return
endif
!
x8=x6*x2;y8=y6*y2;z8=z6*z2
x6y2=x6*y2;y6z2=y6*z2;z6x2=z6*x2
x4y4=x4*y4;y4z4=y4*z4;z4x4=z4*x4
x2y6=x2*y6;y2z6=y2*z6;z2x6=z2*x6
x4y2z2=x4*y2z2;y4z2x2=y4*z2x2;z4x2y2=z4*x2y2
p8=c800*x8+c080*y8+c008*z8 &
    +c620*x6y2+c602*z2x6+c062*y6z2+c260*x2y6+c206*z6x2+c026*y2z6 &
    +c440*x4y4+c422*x4y2z2+c044*y4z4+c242*y4z2x2+c404*z4x4+c224*z4x2y2
a8x6=8.d0*c800*x6
a8x4=6.d0*(c620*x4y2+c602*z2x4)
a8x2=4.d0*(c440*x2y4+c422*x2y2z2+c404*z4x2)
a8x0=2.d0*(c260*y6+c242*y4z2+c224*y2z4+c206*z6)
a8y6=8.d0*c080*y6
a8y4=6.d0*(c062*y4z2+c260*x2y4)
a8y2=4.d0*(c044*y2z4+c242*x2y2z2+c440*x4y2)
a8y0=2.d0*(c026*z6+c224*z4x2+c422*z2x4+c620*x6)
a8z6=8.d0*c008*z6
a8z4=6.d0*(c206*z4x2+c026*y2z4)
a8z2=4.d0*(c404*z2x4+c224*x2y2z2+c044*y4z2)
a8z0=2.d0*(c602*x6+c422*x4y2+c242*x2y4+c062*y6)
p8x=(a8x6+a8x4+a8x2+a8x0)*x
p8y=(a8y6+a8y4+a8y2+a8y0)*y
p8z=(a8z6+a8z4+a8z2+a8z0)*z
if(n.le.4) then
    u=gmri*(1.d0+ri4*(p2+ri4*(p4+ri4*(p6+ri4*p8))))
    r17i4=ri4*17.d0
    ur=-gmri3*(1.d0+ri4*(5.d0*p2+ri4*(9.d0*p4+ri4*(13.d0*p6 &
        +r17i4*p8))))
    ux=ur*x+gmri5*(p2x+ri4*(p4x+ri4*(p6x+ri4*p8x)))
    uy=ur*y+gmri5*(p2y+ri4*(p4y+ri4*(p6y+ri4*p8y)))
    uz=ur*z+gmri5*(p2z+ri4*(p4z+ri4*(p6z+ri4*p8z)))
    return
endif
!
xA=x8*x2;yA=y8*y2;zA=z8*z2
x8y2=x8*y2;y8z2=y8*z2;z8x2=z8*x2
x6y4=x6*y4;y6z4=y6*z4;z6x4=z6*x4
x4y6=x4*y6;y4z6=y4*z6;z4x6=z4*x6
x2y8=x2*y8;y2z8=y2*z8;z2x8=z2*x8
x6y2z2=x6*y2z2;y6z2x2=y6*z2x2;z6x2y2=z6*x2y2
x2y4z4=x2*y4z4;y2z4x4=y2*z4x4;z2x4y4=z2*x4y4
pA=cA00*xA+c0A0*yA+c00A*zA &
    +c820*x8y2+c802*z2x8+c082*y8z2+c280*x2y8+c208*z8x2+c028*y2z8 &
    +c640*x6y4+c622*x6y2z2+c604*z4x6 &
    +c064*y6z4+c262*y6z2x2+c460*x4y6 &
    +c406*z6x4+c226*z6x2y2+c046*y4z6 &
    +c244*x2y4z4+c424*y2z4x4+c442*z2x4y4
aAx8=10.d0*cA00*x8
aAx6=8.d0*(c820*x6y2+c802*z2x6)
aAx4=6.d0*(c640*x4y4+c622*x4y2z2+c604*z4x4)
aAx2=4.d0*(c460*x2y6+c442*y4z2x2+c424*z4x2y2+c406*z6x2)
aAx0=2.d0*(c280*y8+c262*y6z2+c244*y4z4+c226*y2z6+c208*z8)
aAy8=10.d0*c0A0*y8
aAy6=8.d0*(c082*y6z2+c280*x2y6)
aAy4=6.d0*(c064*y4z4+c262*y4z2x2+c460*x4y4)
aAy2=4.d0*(c046*y2z6+c244*z4x2y2+c442*x4y2z2+c640*x6y2)
aAy0=2.d0*(c028*z8+c226*z6x2+c424*z4x4+c622*z2x6+c820*x8)
aAz8=10.d0*c00A*z8
aAz6=8.d0*(c208*z6x2+c028*y2z6)
aAz4=6.d0*(c406*z4x4+c226*z4x2y2+c046*y4z4)
aAz2=4.d0*(c604*z2x6+c424*x4y2z2+c244*y4z2x2+c064*y6z2)
aAz0=2.d0*(c802*x8+c622*x6y2+c442*x4y4+c262*x2y6+c082*y8)
pAx=(aAx8+aAx6+aAx4+aAx2+aAx0)*x
pAy=(aAy8+aAy6+aAy4+aAy2+aAy0)*y
pAz=(aAz8+aAz6+aAz4+aAz2+aAz0)*z
if(n.le.5) then
    u=gmri*(1.d0+ri4*(p2+ri4*(p4+ri4*(p6+ri4*(p8+ri4*pA)))))
    r21i4=ri4*21.d0
    ur=-gmri3*(1.d0+ri4*(5.d0*p2+ri4*(9.d0*p4+ri4*(13.d0*p6 &
        +ri4*(17.d0*p8+r21i4*pA)))))
    ux=ur*x+gmri5*(p2x+ri4*(p4x+ri4*(p6x+ri4*(p8x+ri4*pAx))))
    uy=ur*y+gmri5*(p2y+ri4*(p4y+ri4*(p6y+ri4*(p8y+ri4*pAy))))
    uz=ur*z+gmri5*(p2z+ri4*(p4z+ri4*(p6z+ri4*(p8z+ri4*pAz))))
    return
endif
!
xC=xA*x2;yC=yA*y2;zC=zA*z2
xAy2=xA*y2;yAz2=yA*z2;zAx2=zA*x2
x8y4=x8*y4;y8z4=y8*z4;z8x4=z8*x4
x6y6=x6*y6;y6z6=y6*z6;z6x6=z6*x6
x4y8=x4*y8;y4z8=y4*z8;z4x8=z4*x8
x2yA=x2*yA;y2zA=y2*zA;z2xA=z2*xA
x8y2z2=x8*y2z2;y8z2x2=y8*z2x2;z8x2y2=z8*x2y2
x6y4z2=x6*y4z2;y6z4x2=y6*z4x2;z6x4y2=z6*x4y2
x6y2z4=x6*y2z4;y6z2x4=y6*z2x4;z6x2y4=z6*x2y4
x4y4z4=x2y2z2*x2y2z2
pC=cC00*xC+c0C0*yC+c00C*zC &
    +cA20*xAy2+cA02*z2xA+c0A2*yAz2+c2A0*x2yA+c20A*zAx2+c02A*y2zA &
    +c840*x8y4+c822*x8y2z2+c804*z4x8 &
    +c084*y8z4+c282*y8z2x2+c480*x4y8 &
    +c408*z8x4+c228*z8x2y2+c048*y4z8 &
    +c660*x6y6+c642*x6y4z2+c624*x6y2z4 &
    +c066*y6z6+c264*y6z4x2+c462*y6z2x4 &
    +c606*z6x6+c426*z6x4y2+c246*z6x2y4 &
    +c444*x4y4z4
aCxA=12.d0*cC00*xA
aCx8=10.d0*(cA20*x8y2+cA02*z2x8)
aCx6=8.d0*(c840*x6y4+c822*x6y2z2+c804*z4x6)
aCx4=6.d0*(c660*x4y6+c642*z2x4y4+c624*y2z4x4+c606*z6x4)
aCx2=4.d0*(c480*x2y8+c462*y6z2x2+c444*x2y4z4+c426*z6x2y2+c408*z8x2)
aCx0=2.d0*(c2A0*yA+c282*y8z2+c264*y6z4+c246*y4z6+c228*y2z8+c20A*zA)
aCyA=12.d0*c0C0*yA
aCy8=10.d0*(c0A2*y8z2+c2A0*x2y8)
aCy6=8.d0*(c084*y6z4+c282*y6z2x2+c480*x4y6)
aCy4=6.d0*(c066*y4z6+c264*x2y4z4+c462*z2x4y4+c660*x6y4)
aCy2=4.d0*(c048*y2z8+c246*z6x2y2+c444*y2z4x4+c642*x6y2z2+c840*x8y2)
aCy0=2.d0*(c02A*zA+c228*z8x2+c426*z6x4+c624*z4x6+c822*z2x8+cA20*xA)
aCzA=12.d0*c00C*zA
aCz8=10.d0*(c20A*z8x2+c02A*y2z8)
aCz6=8.d0*(c408*z6x4+c228*z6x2y2+c048*y4z6)
aCz4=6.d0*(c606*z4x6+c426*y2z4x4+c246*x2y4z4+c066*y6z4)
aCz2=4.d0*(c804*z2x8+c624*x6y2z2+c444*z2x4y4+c264*y6z2x2+c084*y8z2)
aCz0=2.d0*(cA02*xA+c822*x8y2+c642*x6y4+c462*x4y6+c282*x2y8+c0A2*yA)
pCx=(aCxA+aCx8+aCx6+aCx4+aCx2+aCx0)*x
pCy=(aCyA+aCy8+aCy6+aCy4+aCy2+aCy0)*y
pCz=(aCzA+aCz8+aCz6+aCz4+aCz2+aCz0)*z
if(n.le.6) then
    u=gmri*(1.d0+ri4*(p2+ri4*(p4+ri4*(p6 &
        +ri4*(p8+ri4*(pA+ri4*pC))))))
    r25i4=ri4*25.d0
    ur=-gmri3*(1.d0+ri4*(5.d0*p2+ri4*(9.d0*p4+ri4*(13.d0*p6 &
        +ri4*(17.d0*p8+ri4*(21.d0*pA+r25i4*pC))))))
    ux=ur*x+gmri5*(p2x+ri4*(p4x+ri4*(p6x+ri4*(p8x &
        +ri4*(pAx+ri4*pCx)))))
    uy=ur*y+gmri5*(p2y+ri4*(p4y+ri4*(p6y+ri4*(p8y &
        +ri4*(pAy+ri4*pCy)))))
    uz=ur*z+gmri5*(p2z+ri4*(p4z+ri4*(p6z+ri4*(p8z &
        +ri4*(pAz+ri4*pCz)))))
    return
endif
!
xE=xC*x2;yE=yC*y2;zE=zC*z2
xCy2=xC*y2;yCz2=yC*z2;zCx2=zC*x2
xAy4=xA*y4;yAz4=yA*z4;zAx4=zA*x4
x8y6=x8*y6;y8z6=y8*z6;z8x6=z8*x6
x6y8=x6*y8;y6z8=y6*z8;z6x8=z6*x8
x4yA=x4*yA;y4zA=y4*zA;z4xA=z4*xA
x2yC=x2*yC;y2zC=y2*zC;z2xC=z2*xC
xAy2z2=xA*y2z2;yAz2x2=yA*z2x2;zAx2y2=zA*x2y2
x8y4z2=x8*y4z2;y8z4x2=y8*z4x2;z8x4y2=z8*x4y2
x8y2z4=x8*y2z4;y8z2x4=y8*z2x4;z8x2y4=z8*x2y4
x2y6z6=x2*y6z6;y2z6x6=y2*z6x6;z2x6y6=z2*x6y6
x6y4z4=x6*y4z4;y6z4x4=y6*z4x4;z6x4y4=z6*x4y4
pE=cE00*xE+c0E0*yE+c00E*zE &
    +cC20*xCy2+cC02*z2xC+c0C2*yCz2+c2C0*x2yC+c20C*zCx2+c02C*y2zC &
    +cA40*xAy4+cA22*xAy2z2+cA04*z4xA &
    +c0A4*yAz4+c2A2*yAz2x2+c4A0*x4yA &
    +c40A*zAx4+c22A*zAx2y2+c04A*y4zA &
    +c860*x8y6+c842*x8y4z2+c824*x8y2z4+c806*z6x8 &
    +c086*y8z6+c284*y8z4x2+c482*y8z2x4+c680*x6y8 &
    +c608*z8x6+c428*z8x4y2+c248*z8x2y4+c068*y6z8 &
    +c662*z2x6y6+c644*x6y4z4 &
    +c266*x2y6z6+c464*y6z4x4 &
    +c626*y2z6x6+c446*z6x4y4
aExC=14.d0*cE00*xC
aExA=12.d0*(cC20*xAy2+cC02*z2xA)
aEx8=10.d0*(cA40*x8y4+cA22*x8y2z2+cA04*z4x8)
aEx6=8.d0*(c860*x6y6+c842*x6y4z2+c824*x6y2z4+c806*z6x6)
aEx4=6.d0*(c680*x4y8+c662*y6z2x4+c644*x4y4z4+c626*z6x4y2+c608*z8x4)
aEx2=4.d0*(c4A0*x2yA+c482*y8z2x2+c464*y6z4x2+c446*z6x2y4+c428*z8x2y2 &
        +c40A*zAx2)
aEx0=2.d0*(c2C0*yC+c2A2*yAz2+c284*y8z4+c266*y6z6+c248*y4z8+c22A*y2zA &
        +c20C*zC)
aEyC=14.d0*c0E0*yC
aEyA=12.d0*(c0C2*yAz2+c2C0*x2yA)
aEy8=10.d0*(c0A4*y8z4+c2A2*y8z2x2+c4A0*x4y8)
aEy6=8.d0*(c086*y6z6+c284*y6z4x2+c482*y6z2x4+c680*x6y6)
aEy4=6.d0*(c068*y4z8+c266*z6x2y4+c464*x4y4z4+c662*x6y4z2+c860*x8y4)
aEy2=4.d0*(c04A*y2zA+c248*z8x2y2+c446*z6x4y2+c644*x6y2z4+c842*x8y2z2 &
        +cA40*xAy2)
aEy0=2.d0*(c02C*zC+c22A*zAx2+c428*z8x4+c626*z6x6+c824*z4x8 &
        +cA22*z2xA+cC20*xC)
aEzC=14.d0*c00E*zC
aEzA=12.d0*(c20C*zAx2+c02C*y2zA)
aEz8=10.d0*(c40A*z8x4+c22A*z8x2y2+c04A*y4z8)
aEz6=8.d0*(c608*z6x6+c428*z6x4y2+c248*z6x2y4+c068*y6z6)
aEz4=6.d0*(c806*z4x8+c626*x6y2z4+c446*x4y4z4+c266*y6z4x2+c086*y8z4)
aEz2=4.d0*(cA04*z2xA+c824*x8y2z2+c644*x6y4z2+c464*y6z2x4+c284*y8z2x2 &
        +c0A4*yAz2)
aEz0=2.d0*(cC02*xC+cA22*xAy2+c842*x8y4+c662*x6y6+c482*x4y8 &
        +c2A2*x2yA+c0C2*yC)
pEx=(aExC+aExA+aEx8+aEx6+aEx4+aEx2+aEx0)*x
pEy=(aEyC+aEyA+aEy8+aEy6+aEy4+aEy2+aEy0)*y
pEz=(aEzC+aEzA+aEz8+aEz6+aEz4+aEz2+aEz0)*z
if(n.le.7) then
    u=gmri*(1.d0+ri4*(p2+ri4*(p4+ri4*(p6 &
        +ri4*(p8+ri4*(pA+ri4*(pC+ri4*pE)))))))
    r29i4=ri4*29.d0
    ur=-gmri3*(1.d0+ri4*(5.d0*p2+ri4*(9.d0*p4+ri4*(13.d0*p6 &
        +ri4*(17.d0*p8+ri4*(21.d0*pA+ri4*(25.d0*pC &
        +r29i4*pE)))))))
    ux=ur*x+gmri5*(p2x+ri4*(p4x+ri4*(p6x+ri4*(p8x &
        +ri4*(pAx+ri4*(pCx+ri4*pEx))))))
    uy=ur*y+gmri5*(p2y+ri4*(p4y+ri4*(p6y+ri4*(p8y &
        +ri4*(pAy+ri4*(pCy+ri4*pEy))))))
    uz=ur*z+gmri5*(p2z+ri4*(p4z+ri4*(p6z+ri4*(p8z &
        +ri4*(pAz+ri4*(pCz+ri4*pEz))))))
    return
endif
!
xG=xE*x2;yG=yE*y2;zG=zE*z2
xEy2=xE*y2;yEz2=yE*z2;zEx2=zE*x2
xCy4=xC*y4;yCz4=yC*z4;zCx4=zC*x4
xAy6=xA*y6;yAz6=yA*z6;zAx6=zA*x6
x8y8=x8*y8;y8z8=y8*z8;z8x8=z8*x8
x6yA=x6*yA;y6zA=y6*zA;z6xA=z6*xA
x4yC=x4*yC;y4zC=y4*zC;z4xC=z4*xC
x2yE=x2*yE;y2zE=y2*zE;z2xE=z2*xE
xCy2z2=xC*y2z2;yCz2x2=yC*z2x2;zCx2y2=zC*x2y2
xAy4z2=xA*y4z2;yAz4x2=yA*z4x2;zAx4y2=zA*x4y2
xAy2z4=xA*y2z4;yAz2x4=yA*z2x4;zAx2y4=zA*x2y4
x8y6z2=x8*y6z2;y8z6x2=y8*z6x2;z8x6y2=z8*x6y2
x8y4z4=x8*y4z4;y8z4x4=y8*z4x4;z8x4y4=z8*x4y4
x8y2z6=x8*y2z6;y8z2x6=y8*z2x6;z8x2y6=z8*x2y6
x4y6z6=x4*y6z6;y4z6x6=y4*z6x6;z4x6y6=z4*x6y6
pG=cG00*xG+c0G0*yG+c00G*zG &
    +cE20*xEy2+cE02*z2xE+c0E2*yEz2+c2E0*x2yE+c20E*zEx2+c02E*y2zE &
    +cC40*xCy4+cC22*xCy2z2+cC04*z4xC &
    +c0C4*yCz4+c2C2*yCz2x2+c4C0*x4yC &
    +c40C*zCx4+c22C*zCx2y2+c04C*y4zC &
    +cA60*xAy6+cA42*xAy4z2+cA24*xAy2z4+cA06*z6xA &
    +c0A6*yAz6+c2A4*yAz4x2+c4A2*yAz2x4+c6A0*x6yA &
    +c60A*zAx6+c42A*zAx4y2+c24A*zAx2y4+c06A*y6zA &
    +c880*x8y8+c862*x8y6z2+c844*x8y4z4+c826*x8y2z6 &
    +c088*y8z8+c286*y8z6x2+c484*y8z4x4+c682*y8z2x6 &
    +c808*z8x8+c628*z8x6y2+c448*z8x4y4+c268*z8x2y6 &
    +c466*x4y6z6+c646*y4z6x6+c664*z4x6y6
aGxE=16.d0*cG00*xE
aGxC=14.d0*(cE20*xCy2+cE02*z2xC)
aGxA=12.d0*(cC40*xAy4+cC22*xAy2z2+cC04*z4xA)
aGx8=10.d0*(cA60*x8y6+cA42*x8y4z2+cA24*x8y2z4+cA06*z6x8)
aGx6=8.d0*(c880*x6y8+c862*z2x6y6+c844*x6y4z4+c826*y2z6x6+c808*z8x6)
aGx4=6.d0*(c6A0*x4yA+c682*y8z2x4+c664*y6z4x4+c646*z6x4y4+c628*z8x4y2 &
        +c60A*zAx4)
aGx2=4.d0*(c4C0*x2yC+c4A2*yAz2x2+c484*y8z4x2+c466*x2y6z6+c448*z8x2y4 &
        +c42A*zAx2y2+c40C*zCx2)
aGx0=2.d0*(c2E0*yE+c2C2*yCz2+c2A4*yAz4+c286*y8z6+c268*y6z8 &
        +c24A*y4zA+c22C*y2zC+c20E*zE)
aGyE=16.d0*c0G0*yE
aGyC=14.d0*(c0E2*yCz2+c2E0*x2yC)
aGyA=12.d0*(c0C4*yAz4+c2C2*yAz2x2+c4C0*x4yA)
aGy8=10.d0*(c0A6*y8z6+c2A4*y8z4x2+c4A2*y8z2x4+c6A0*x6y8)
aGy6=8.d0*(c088*y6z8+c286*x2y6z6+c484*y6z4x4+c682*z2x6y6+c880*x8y6)
aGy4=6.d0*(c06A*y4zA+c268*z8x2y4+c466*z6x4y4+c664*x6y4z4+c862*x8y4z2 &
        +cA60*xAy4)
aGy2=4.d0*(c04C*y2zC+c24A*zAx2y2+c448*z8x4y2+c646*y2z6x6+c844*x8y2z4 &
        +cA42*xAy2z2+cC40*xCy2)
aGy0=2.d0*(c02E*zE+c22C*zCx2+c42A*zAx4+c628*z8x6+c826*z6x8 &
        +cA24*z4xA+cC22*z2xC+cE20*xE)
aGzE=16.d0*c00G*zE
aGzC=14.d0*(c20E*zCx2+c02E*y2zC)
aGzA=12.d0*(c40C*zAx4+c22C*zAx2y2+c04C*y4zA)
aGz8=10.d0*(c60A*z8x6+c42A*z8x4y2+c24A*z8x2y4+c06A*y6z8)
aGz6=8.d0*(c808*z6x8+c628*y2z6x6+c448*z6x4y4+c268*x2y6z6+c088*y8z6)
aGz4=6.d0*(cA06*z4xA+c826*x8y2z4+c646*x6y4z4+c466*y6z4x4+c286*y8z4x2 &
        +c0A6*yAz4)
aGz2=4.d0*(cC04*z2xC+cA24*xAy2z2+c844*x8y4z2+c664*z2x6y6+c484*y8z2x4 &
        +c2A4*yAz2x2+c0C4*yCz2)
aGz0=2.d0*(cE02*xE+cC22*xCy2+cA42*xAy4+c862*x8y6+c682*x6y8 &
        +c4A2*x4yA+c2C2*x2yC+c0E2*yE)
pGx=(aGxE+aGxC+aGxA+aGx8+aGx6+aGx4+aGx2+aGx0)*x
pGy=(aGyE+aGyC+aGyA+aGy8+aGy6+aGy4+aGy2+aGy0)*y
pGz=(aGzE+aGzC+aGzA+aGz8+aGz6+aGz4+aGz2+aGz0)*z
u=gmri*(1.d0+ri4*(p2+ri4*(p4+ri4*(p6 &
    +ri4*(p8+ri4*(pA+ri4*(pC+ri4*(pE+ri4*pG))))))))
r33i4=ri4*33.d0
ur=-gmri3*(1.d0+ri4*(5.d0*p2+ri4*(9.d0*p4+ri4*(13.d0*p6 &
    +ri4*(17.d0*p8+ri4*(21.d0*pA+ri4*(25.d0*pC &
    +ri4*(29.d0*pE+r33i4*pG))))))))
ux=ur*x+gmri5*(p2x+ri4*(p4x+ri4*(p6x+ri4*(p8x &
    +ri4*(pAx+ri4*(pCx+ri4*(pEx+ri4*pGx)))))))
uy=ur*y+gmri5*(p2y+ri4*(p4y+ri4*(p6y+ri4*(p8y &
    +ri4*(pAy+ri4*(pCy+ri4*(pEy+ri4*pGy)))))))
uz=ur*z+gmri5*(p2z+ri4*(p4z+ri4*(p6z+ri4*(p8z &
   +ri4*(pAz+ri4*(pCz+ri4*(pEz+ri4*pGz)))))))
!
return;end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tnggpri(n,x,y,z,u,ux,uy,uz,uxx,uxy,uxz,uyy,uyz,uzz)
!
! Computation of graviational potential, acceleration vector, and
! gravity gradient tensor of a uniform rectangular prism
! by an even order truncated Taylor expansion
!
! Reference: Fukushima, T, (2019) Geophys. J. Int'l, revised
!   Taylor expansion of prismatic gravitational field
!
! Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
!
!     Caution: Coefficients, gm and cijk, must be prepared beforehand
!
!     Inputs: n = truncated order                0 <= n <= 16
!             x,y,z = relative coordinates of evaluated point
!
!     Output: u = gravitational potential 
!             ux = x-component gravitational acceleration
!             uy = y-component gravitational acceleration
!             uz = z-component gravitational acceleration
!             uxx = xx-component gravity gradient
!             uxy = xy-component gravity gradient
!             uxz = xz-component gravity gradient
!             uyy = yy-component gravity gradient
!             uyz = yz-component gravity gradient
!             uzz = zz-component gravity gradient
!
implicit real*8 (a-h,o-z)
implicit integer (i-n)
!
common /cset0/gm
common /cset2/c200,c020,c002
common /cset4/c400,c040,c004,c022,c202,c220
common /cset6/c600,c060,c006,c420,c042,c204,c402,c024,c240,c222
common /cset8/c800,c080,c008,c620,c062,c206,c602,c026,c260, &
              c044,c404,c440,c422,c242,c224
common /csetA/cA00,c0A0,c00A,c820,c082,c208,c802,c028,c280, &
              c640,c064,c406,c604,c046,c460,c622,c262,c226, &
              c244,c424,c442
common /csetC/cC00,c0C0,c00C,cA20,c0A2,c20A,cA02,c02A,c2A0, &
              c840,c084,c408,c804,c048,c480,c822,c282,c228, &
              c066,c606,c660,c642,c264,c426,c624,c246,c462,c444
common /csetE/cE00,c0E0,c00E,cC20,c0C2,c20C,cC02,c02C,c2C0, &
              cA40,c0A4,c40A,cA04,c04A,c4A0,cA22,c2A2,c22A, &
              c860,c086,c608,c806,c068,c680,c842,c284,c428, &
              c824,c248,c482,c266,c626,c662,c644,c464,c446
common /csetG/cG00,c0G0,c00G,cE20,c0E2,c20E,cE02,c02E,c2E0, &
              cC40,c0C4,c40C,cc04,c04C,c4c0,cC22,c2C2,c22C, &
              cA60,c0A6,c60A,cA06,c06A,c6A0,cA42,c2A4,c42A, &
              cA24,c24A,c4A2,c088,c808,c880,c862,c286,c628, &
              c826,c268,c682,c844,c484,c448,c466,c646,c664
!
x2=x*x;y2=y*y;z2=z*z
xy=x*y;yz=y*z;zx=z*x
r2=x2+y2+z2
ri2=1.d0/r2
ri=sqrt(ri2)
gmri=gm*ri
ri3=ri*ri2
gmri3=gm*ri3
ri5=ri3*ri2
gmri5=gm*ri5
gm3ri5=gmri5*3.d0
if(n.le.0) then
    u=gmri
    ur=-gmri3
    ux=ur*x;uy=ur*y;uz=ur*z
    uxx=gmri5*(2.d0*x2-y2-z2)
    uyy=gmri5*(2.d0*y2-x2-z2)
    uzz=gmri5*(2.d0*z2-x2-y2)
    uxy=gm3ri5*xy
    uyz=gm3ri5*yz
    uxz=gm3ri5*zx
    return
endif
!
ri4=ri2*ri2
ri7=ri5*ri2
gmri7=gm*ri7
gm2ri7=gmri7*2.d0
p2=c200*x2+c020*y2+c002*z2
p2xx=2.d0*c200
p2yy=2.d0*c020
p2zz=2.d0*c002
p2x=p2xx*x
p2y=p2yy*y
p2z=p2zz*z
if(n.le.1) then
    u=gmri*(1.d0+ri4*p2)
    ur=-gmri3*(1.d0+ri4*5.d0*p2)
    ux=ur*x+gmri5*p2x
    uy=ur*y+gmri5*p2y
    uz=ur*z+gmri5*p2z
    urr=gmri5*(3.d0+ri4*35.d0*p2)
    wx=5.d0*p2x
    wy=5.d0*p2y
    wz=5.d0*p2z
    uxx=-gm2ri7*wx*x
    uyy=-gm2ri7*wy*y
    uzz=-gm2ri7*wz*z
    uxy=-gmri7*(wx*y+wy*x)
    uyz=-gmri7*(wy*z+wz*y)
    uxz=-gmri7*(wz*x+wx*z)
    uxx=uxx+urr*x2+ur+gmri5*p2xx
    uyy=uyy+urr*y2+ur+gmri5*p2yy
    uzz=uzz+urr*z2+ur+gmri5*p2zz
    uxy=uxy+urr*xy
    uyz=uyz+urr*yz
    uxz=uxz+urr*zx
    return
endif
!
ri9=ri7*ri2
gmri9=gm*ri9
x4=x2*x2;y4=y2*y2;z4=z2*z2
y2z2=y2*z2;z2x2=z2*x2;x2y2=x2*y2
p4=c400*x4+c040*y4+c004*z4+c022*y2z2+c202*z2x2+c220*x2y2
a4x2=4.d0*c400*x2
a4x0=2.d0*(c202*z2+c220*y2)
a4y2=4.d0*c040*y2
a4y0=2.d0*(c220*x2+c022*z2)
a4z2=4.d0*c004*z2
a4z0=2.d0*(c022*y2+c202*x2)
p4x=(a4x2+a4x0)*x
p4y=(a4y2+a4y0)*y
p4z=(a4z2+a4z0)*z
p4xx=3.d0*a4x2+a4x0
p4yy=3.d0*a4y2+a4y0
p4zz=3.d0*a4z2+a4z0
p4xy=4.d0*c220*xy
p4yz=4.d0*c022*yz
p4zx=4.d0*c202*zx
if(n.le.2) then
    u=gmri*(1.d0+ri4*(p2+ri4*p4))
    r9i4=ri4*9.d0
    ur=-gmri3*(1.d0+ri4*(5.d0*p2+r9i4*p4))
    ux=ur*x+gmri5*(p2x+ri4*p4x)
    uy=ur*y+gmri5*(p2y+ri4*p4y)
    uz=ur*z+gmri5*(p2z+ri4*p4z)
    urr=gmri5*(3.d0+ri4*(35.d0*p2+ri4*99.d0*p4))
    wx=5.d0*p2x+r9i4*p4x
    wy=5.d0*p2y+r9i4*p4y
    wz=5.d0*p2z+r9i4*p4z
    uxx=-gm2ri7*wx*x
    uyy=-gm2ri7*wy*y
    uzz=-gm2ri7*wz*z
    uxy=-gmri7*(wx*y+wy*x)
    uyz=-gmri7*(wy*z+wz*y)
    uxz=-gmri7*(wz*x+wx*z)
    uxx=uxx+urr*x2+ur+gmri5*(p2xx+ri4*p4xx)
    uyy=uyy+urr*y2+ur+gmri5*(p2yy+ri4*p4yy)
    uzz=uzz+urr*z2+ur+gmri5*(p2zz+ri4*p4zz)
    uxy=uxy+urr*xy+gmri9*p4xy
    uyz=uyz+urr*yz+gmri9*p4yz
    uxz=uxz+urr*zx+gmri9*p4zx
    return
endif
!
x6=x4*x2;y6=y4*y2;z6=z4*z2
x4y2=x4*y2;y4z2=y4*z2;z4x2=z4*x2
x2y4=x2*y4;y2z4=y2*z4;z2x4=z2*x4
x2y2z2=x2y2*z2
p6=c600*x6+c060*y6+c006*z6+c222*x2y2z2 &
    +c420*x4y2+c402*z2x4+c042*y4z2+c240*x2y4+c204*z4x2+c024*y2z4
a6x4=6.d0*c600*x4
a6x2=4.d0*(c420*x2y2+c402*z2x2)
a6x0=2.d0*(c240*y4+c222*y2z2+c204*z4)
a6y4=6.d0*c060*y4
a6y2=4.d0*(c042*y2z2+c240*x2y2)
a6y0=2.d0*(c024*z4+c222*z2x2+c420*x4)
a6z4=6.d0*c006*z4
a6z2=4.d0*(c204*z2x2+c024*y2z2)
a6z0=2.d0*(c402*x4+c222*x2y2+c042*y4)
p6x=(a6x4+a6x2+a6x0)*x
p6y=(a6y4+a6y2+a6y0)*y
p6z=(a6z4+a6z2+a6z0)*z
p6xx=5.d0*a6x4+3.d0*a6x2+a6x0
p6yy=5.d0*a6y4+3.d0*a6y2+a6y0
p6zz=5.d0*a6z4+3.d0*a6z2+a6z0
a4c2=4.d0*c222
p6xy=(8.d0*(c420*x2+c240*y2)+a4c2*z2)*xy
p6yz=(8.d0*(c042*y2+c024*z2)+a4c2*x2)*yz
p6zx=(8.d0*(c204*z2+c402*x2)+a4c2*y2)*zx
if(n.le.3) then
    u=gmri*(1.d0+ri4*(p2+ri4*(p4+ri4*p6)))
    r13i4=ri4*13.d0
    ur=-gmri3*(1.d0+ri4*(5.d0*p2+ri4*(9.d0*p4+r13i4*p6)))
    ux=ur*x+gmri5*(p2x+ri4*(p4x+ri4*p6x))
    uy=ur*y+gmri5*(p2y+ri4*(p4y+ri4*p6y))
    uz=ur*z+gmri5*(p2z+ri4*(p4z+ri4*p6z))
    urr=gmri5*(3.d0+ri4*(35.d0*p2+ri4*(99.d0*p4+ri4*195.d0*p6)))
    wx=5.d0*p2x+ri4*(9.d0*p4x+r13i4*p6x)
    wy=5.d0*p2y+ri4*(9.d0*p4y+r13i4*p6y)
    wz=5.d0*p2z+ri4*(9.d0*p4z+r13i4*p6z)
    uxx=-gm2ri7*wx*x
    uyy=-gm2ri7*wy*y
    uzz=-gm2ri7*wz*z
    uxy=-gmri7*(wx*y+wy*x)
    uyz=-gmri7*(wy*z+wz*y)
    uxz=-gmri7*(wz*x+wx*z)
    uxx=uxx+urr*x2+ur+gmri5*(p2xx+ri4*(p4xx+ri4*p6xx))
    uyy=uyy+urr*y2+ur+gmri5*(p2yy+ri4*(p4yy+ri4*p6yy))
    uzz=uzz+urr*z2+ur+gmri5*(p2zz+ri4*(p4zz+ri4*p6zz))
    uxy=uxy+urr*xy+gmri9*(p4xy+ri4*p6xy)
    uyz=uyz+urr*yz+gmri9*(p4yz+ri4*p6yz)
    uxz=uxz+urr*zx+gmri9*(p4zx+ri4*p6zx)
    return
endif
!
x8=x6*x2;y8=y6*y2;z8=z6*z2
x6y2=x6*y2;y6z2=y6*z2;z6x2=z6*x2
x4y4=x4*y4;y4z4=y4*z4;z4x4=z4*x4
x2y6=x2*y6;y2z6=y2*z6;z2x6=z2*x6
x4y2z2=x4*y2z2;y4z2x2=y4*z2x2;z4x2y2=z4*x2y2
p8=c800*x8+c080*y8+c008*z8 &
    +c620*x6y2+c602*z2x6+c062*y6z2+c260*x2y6+c206*z6x2+c026*y2z6 &
    +c440*x4y4+c422*x4y2z2+c044*y4z4+c242*y4z2x2+c404*z4x4+c224*z4x2y2
a8x6=8.d0*c800*x6
a8x4=6.d0*(c620*x4y2+c602*z2x4)
a8x2=4.d0*(c440*x2y4+c422*x2y2z2+c404*z4x2)
a8x0=2.d0*(c260*y6+c242*y4z2+c224*y2z4+c206*z6)
a8y6=8.d0*c080*y6
a8y4=6.d0*(c062*y4z2+c260*x2y4)
a8y2=4.d0*(c044*y2z4+c242*x2y2z2+c440*x4y2)
a8y0=2.d0*(c026*z6+c224*z4x2+c422*z2x4+c620*x6)
a8z6=8.d0*c008*z6
a8z4=6.d0*(c206*z4x2+c026*y2z4)
a8z2=4.d0*(c404*z2x4+c224*x2y2z2+c044*y4z2)
a8z0=2.d0*(c602*x6+c422*x4y2+c242*x2y4+c062*y6)
p8x=(a8x6+a8x4+a8x2+a8x0)*x
p8y=(a8y6+a8y4+a8y2+a8y0)*y
p8z=(a8z6+a8z4+a8z2+a8z0)*z
p8xx=7.d0*a8x6+5.d0*a8x4+3.d0*a8x2+a8x0
p8yy=7.d0*a8y6+5.d0*a8y4+3.d0*a8y2+a8y0
p8zz=7.d0*a8z6+5.d0*a8z4+3.d0*a8z2+a8z0
p8xy=(12.d0*c620*x4+8.d0*(2.d0*c440*x2y2+c422*z2x2) &
    +4.d0*(3.d0*c260*y4+2.d0*c242*y2z2+c224*z4))*xy
p8yz=(12.d0*c062*y4+8.d0*(2.d0*c044*y2z2+c242*x2y2) &
    +4.d0*(3.d0*c026*z4+2.d0*c224*z2x2+c422*x4))*yz
p8zx=(12.d0*c206*z4+8.d0*(2.d0*c404*z2x2+c224*y2z2) &
    +4.d0*(3.d0*c602*x4+2.d0*c422*x2y2+c242*y4))*zx
if(n.le.4) then
    u=gmri*(1.d0+ri4*(p2+ri4*(p4+ri4*(p6+ri4*p8))))
    r17i4=ri4*17.d0
    ur=-gmri3*(1.d0+ri4*(5.d0*p2+ri4*(9.d0*p4+ri4*(13.d0*p6 &
        +r17i4*p8))))
    ux=ur*x+gmri5*(p2x+ri4*(p4x+ri4*(p6x+ri4*p8x)))
    uy=ur*y+gmri5*(p2y+ri4*(p4y+ri4*(p6y+ri4*p8y)))
    uz=ur*z+gmri5*(p2z+ri4*(p4z+ri4*(p6z+ri4*p8z)))
    urr=gmri5*(3.d0+ri4*(35.d0*p2+ri4*(99.d0*p4+ri4*(195.d0*p6 &
        +ri4*323.d0*p8))))
    wx=5.d0*p2x+ri4*(9.d0*p4x+ri4*(13.d0*p6x+r17i4*p8x))
    wy=5.d0*p2y+ri4*(9.d0*p4y+ri4*(13.d0*p6y+r17i4*p8y))
    wz=5.d0*p2z+ri4*(9.d0*p4z+ri4*(13.d0*p6z+r17i4*p8z))
    uxx=-gm2ri7*wx*x
    uyy=-gm2ri7*wy*y
    uzz=-gm2ri7*wz*z
    uxy=-gmri7*(wx*y+wy*x)
    uyz=-gmri7*(wy*z+wz*y)
    uxz=-gmri7*(wz*x+wx*z)
    uxx=uxx+urr*x2+ur+gmri5*(p2xx+ri4*(p4xx+ri4*(p6xx+ri4*p8xx)))
    uyy=uyy+urr*y2+ur+gmri5*(p2yy+ri4*(p4yy+ri4*(p6yy+ri4*p8yy)))
    uzz=uzz+urr*z2+ur+gmri5*(p2zz+ri4*(p4zz+ri4*(p6zz+ri4*p8zz)))
    uxy=uxy+urr*xy+gmri9*(p4xy+ri4*(p6xy+ri4*p8xy))
    uyz=uyz+urr*yz+gmri9*(p4yz+ri4*(p6yz+ri4*p8yz))
    uxz=uxz+urr*zx+gmri9*(p4zx+ri4*(p6zx+ri4*p8zx))
    return
endif
!
xA=x8*x2;yA=y8*y2;zA=z8*z2
x8y2=x8*y2;y8z2=y8*z2;z8x2=z8*x2
x6y4=x6*y4;y6z4=y6*z4;z6x4=z6*x4
x4y6=x4*y6;y4z6=y4*z6;z4x6=z4*x6
x2y8=x2*y8;y2z8=y2*z8;z2x8=z2*x8
x6y2z2=x6*y2z2;y6z2x2=y6*z2x2;z6x2y2=z6*x2y2
x2y4z4=x2*y4z4;y2z4x4=y2*z4x4;z2x4y4=z2*x4y4
pA=cA00*xA+c0A0*yA+c00A*zA &
    +c820*x8y2+c802*z2x8+c082*y8z2+c280*x2y8+c208*z8x2+c028*y2z8 &
    +c640*x6y4+c622*x6y2z2+c604*z4x6 &
    +c064*y6z4+c262*y6z2x2+c460*x4y6 &
    +c406*z6x4+c226*z6x2y2+c046*y4z6 &
    +c244*x2y4z4+c424*y2z4x4+c442*z2x4y4
aAx8=10.d0*cA00*x8
aAx6=8.d0*(c820*x6y2+c802*z2x6)
aAx4=6.d0*(c640*x4y4+c622*x4y2z2+c604*z4x4)
aAx2=4.d0*(c460*x2y6+c442*y4z2x2+c424*z4x2y2+c406*z6x2)
aAx0=2.d0*(c280*y8+c262*y6z2+c244*y4z4+c226*y2z6+c208*z8)
aAy8=10.d0*c0A0*y8
aAy6=8.d0*(c082*y6z2+c280*x2y6)
aAy4=6.d0*(c064*y4z4+c262*y4z2x2+c460*x4y4)
aAy2=4.d0*(c046*y2z6+c244*z4x2y2+c442*x4y2z2+c640*x6y2)
aAy0=2.d0*(c028*z8+c226*z6x2+c424*z4x4+c622*z2x6+c820*x8)
aAz8=10.d0*c00A*z8
aAz6=8.d0*(c208*z6x2+c028*y2z6)
aAz4=6.d0*(c406*z4x4+c226*z4x2y2+c046*y4z4)
aAz2=4.d0*(c604*z2x6+c424*x4y2z2+c244*y4z2x2+c064*y6z2)
aAz0=2.d0*(c802*x8+c622*x6y2+c442*x4y4+c262*x2y6+c082*y8)
pAx=(aAx8+aAx6+aAx4+aAx2+aAx0)*x
pAy=(aAy8+aAy6+aAy4+aAy2+aAy0)*y
pAz=(aAz8+aAz6+aAz4+aAz2+aAz0)*z
pAxx=9.d0*aAx8+7.d0*aAx6+5.d0*aAx4+3.d0*aAx2+aAx0
pAyy=9.d0*aAy8+7.d0*aAy6+5.d0*aAy4+3.d0*aAy2+aAy0
pAzz=9.d0*aAz8+7.d0*aAz6+5.d0*aAz4+3.d0*aAz2+aAz0
pAxy=(16.d0*c820*x6+12.d0*(2.d0*c640*x4y2+c622*z2x4) &
    +8.d0*(3.d0*c460*x2y4+2.d0*c442*x2y2z2+c424*z4x2) &
    +4.d0*(4.d0*c280*y6+3.d0*c262*y4z2+2.d0*c244*y2z4+c226*z6))*xy
pAyz=(16.d0*c082*y6+12.d0*(2.d0*c064*y4z2+c262*x2y4) &
    +8.d0*(3.d0*c046*y2z4+2.d0*c244*x2y2z2+c442*x4y2) &
    +4.d0*(4.d0*c028*z6+3.d0*c226*z4x2+2.d0*c424*z2x4+c622*x6))*yz
pAzx=(16.d0*c208*z6+12.d0*(2.d0*c406*z4x2+c226*y2z4) &
    +8.d0*(3.d0*c604*z2x4+2.d0*c424*x2y2z2+c244*y4z2) &
    +4.d0*(4.d0*c802*x6+3.d0*c622*x4y2+2.d0*c442*x2y4+c262*y6))*zx
if(n.le.5) then
    u=gmri*(1.d0+ri4*(p2+ri4*(p4+ri4*(p6+ri4*(p8+ri4*pA)))))
    r21i4=ri4*21.d0
    ur=-gmri3*(1.d0+ri4*(5.d0*p2+ri4*(9.d0*p4+ri4*(13.d0*p6 &
        +ri4*(17.d0*p8+r21i4*pA)))))
    ux=ur*x+gmri5*(p2x+ri4*(p4x+ri4*(p6x+ri4*(p8x+ri4*pAx))))
    uy=ur*y+gmri5*(p2y+ri4*(p4y+ri4*(p6y+ri4*(p8y+ri4*pAy))))
    uz=ur*z+gmri5*(p2z+ri4*(p4z+ri4*(p6z+ri4*(p8z+ri4*pAz))))
    urr=gmri5*(3.d0+ri4*(35.d0*p2+ri4*(99.d0*p4+ri4*(195.d0*p6 &
        +ri4*(323.d0*p8+ri4*483.d0*pA)))))
    wx=5.d0*p2x+ri4*(9.d0*p4x+ri4*(13.d0*p6x+ri4*(17.d0*p8x &
        +r21i4*pAx)))
    wy=5.d0*p2y+ri4*(9.d0*p4y+ri4*(13.d0*p6y+ri4*(17.d0*p8y &
        +r21i4*pAy)))
    wz=5.d0*p2z+ri4*(9.d0*p4z+ri4*(13.d0*p6z+ri4*(17.d0*p8z &
        +r21i4*pAz)))
    uxx=-gm2ri7*wx*x
    uyy=-gm2ri7*wy*y
    uzz=-gm2ri7*wz*z
    uxy=-gmri7*(wx*y+wy*x)
    uyz=-gmri7*(wy*z+wz*y)
    uxz=-gmri7*(wz*x+wx*z)
    uxx=uxx+urr*x2+ur+gmri5*(p2xx+ri4*(p4xx+ri4*(p6xx+ri4*(p8xx &
        +ri4*pAxx))))
    uyy=uyy+urr*y2+ur+gmri5*(p2yy+ri4*(p4yy+ri4*(p6yy+ri4*(p8yy &
        +ri4*pAyy))))
    uzz=uzz+urr*z2+ur+gmri5*(p2zz+ri4*(p4zz+ri4*(p6zz+ri4*(p8zz &
        +ri4*pAzz))))
    uxy=uxy+urr*xy+gmri9*(p4xy+ri4*(p6xy+ri4*(p8xy+ri4*pAxy)))
    uyz=uyz+urr*yz+gmri9*(p4yz+ri4*(p6yz+ri4*(p8yz+ri4*pAyz)))
    uxz=uxz+urr*zx+gmri9*(p4zx+ri4*(p6zx+ri4*(p8zx+ri4*pAzx)))
    return
endif
!
xC=xA*x2;yC=yA*y2;zC=zA*z2
xAy2=xA*y2;yAz2=yA*z2;zAx2=zA*x2
x8y4=x8*y4;y8z4=y8*z4;z8x4=z8*x4
x6y6=x6*y6;y6z6=y6*z6;z6x6=z6*x6
x4y8=x4*y8;y4z8=y4*z8;z4x8=z4*x8
x2yA=x2*yA;y2zA=y2*zA;z2xA=z2*xA
x8y2z2=x8*y2z2;y8z2x2=y8*z2x2;z8x2y2=z8*x2y2
x6y4z2=x6*y4z2;y6z4x2=y6*z4x2;z6x4y2=z6*x4y2
x6y2z4=x6*y2z4;y6z2x4=y6*z2x4;z6x2y4=z6*x2y4
x4y4z4=x2y2z2*x2y2z2
pC=cC00*xC+c0C0*yC+c00C*zC &
    +cA20*xAy2+cA02*z2xA+c0A2*yAz2+c2A0*x2yA+c20A*zAx2+c02A*y2zA &
    +c840*x8y4+c822*x8y2z2+c804*z4x8 &
    +c084*y8z4+c282*y8z2x2+c480*x4y8 &
    +c408*z8x4+c228*z8x2y2+c048*y4z8 &
    +c660*x6y6+c642*x6y4z2+c624*x6y2z4 &
    +c066*y6z6+c264*y6z4x2+c462*y6z2x4 &
    +c606*z6x6+c426*z6x4y2+c246*z6x2y4 &
    +c444*x4y4z4
aCxA=12.d0*cC00*xA
aCx8=10.d0*(cA20*x8y2+cA02*z2x8)
aCx6=8.d0*(c840*x6y4+c822*x6y2z2+c804*z4x6)
aCx4=6.d0*(c660*x4y6+c642*z2x4y4+c624*y2z4x4+c606*z6x4)
aCx2=4.d0*(c480*x2y8+c462*y6z2x2+c444*x2y4z4+c426*z6x2y2+c408*z8x2)
aCx0=2.d0*(c2A0*yA+c282*y8z2+c264*y6z4+c246*y4z6+c228*y2z8+c20A*zA)
aCyA=12.d0*c0C0*yA
aCy8=10.d0*(c0A2*y8z2+c2A0*x2y8)
aCy6=8.d0*(c084*y6z4+c282*y6z2x2+c480*x4y6)
aCy4=6.d0*(c066*y4z6+c264*x2y4z4+c462*z2x4y4+c660*x6y4)
aCy2=4.d0*(c048*y2z8+c246*z6x2y2+c444*y2z4x4+c642*x6y2z2+c840*x8y2)
aCy0=2.d0*(c02A*zA+c228*z8x2+c426*z6x4+c624*z4x6+c822*z2x8+cA20*xA)
aCzA=12.d0*c00C*zA
aCz8=10.d0*(c20A*z8x2+c02A*y2z8)
aCz6=8.d0*(c408*z6x4+c228*z6x2y2+c048*y4z6)
aCz4=6.d0*(c606*z4x6+c426*y2z4x4+c246*x2y4z4+c066*y6z4)
aCz2=4.d0*(c804*z2x8+c624*x6y2z2+c444*z2x4y4+c264*y6z2x2+c084*y8z2)
aCz0=2.d0*(cA02*xA+c822*x8y2+c642*x6y4+c462*x4y6+c282*x2y8+c0A2*yA)
pCx=(aCxA+aCx8+aCx6+aCx4+aCx2+aCx0)*x
pCy=(aCyA+aCy8+aCy6+aCy4+aCy2+aCy0)*y
pCz=(aCzA+aCz8+aCz6+aCz4+aCz2+aCz0)*z
pCxx=11.d0*aCxA+9.d0*aCx8+7.d0*aCx6+5.d0*aCx4+3.d0*aCx2+aCx0
pCyy=11.d0*aCyA+9.d0*aCy8+7.d0*aCy6+5.d0*aCy4+3.d0*aCy2+aCy0
pCzz=11.d0*aCzA+9.d0*aCz8+7.d0*aCz6+5.d0*aCz4+3.d0*aCz2+aCz0
pCxy=(20.d0*cA20*x8+16.d0*(2.d0*c840*x6y2+c822*z2x6) &
    +12.d0*(3.d0*c660*x4y4+2.d0*c642*x4y2z2+c624*z4x4) &
    +8.d0*(4.d0*c480*x2y6+3.d0*c462*y4z2x2+2.d0*c444*z4x2y2+c426*z6x2) &
    +4.d0*(5.d0*c2A0*y8+4.d0*c282*y6z2+3.d0*c264*y4z4 &
        +2.d0*c246*y2z6+c228*z8))*xy
pCyz=(20.d0*c0A2*y8+16.d0*(2.d0*c084*y6z2+c282*x2y6) &
    +12.d0*(3.d0*c066*y4z4+2.d0*c264*y4z2x2+c462*x4y4) &
    +8.d0*(4.d0*c048*y2z6+3.d0*c246*z4x2y2+2.d0*c444*x4y2z2+c642*x6y2) &
    +4.d0*(5.d0*c02A*z8+4.d0*c228*z6x2+3.d0*c426*z4x4 &
        +2.d0*c624*z2x6+c822*x8))*yz
pCzx=(20.d0*c20A*z8+16.d0*(2.d0*c408*z6x2+c228*y2z6) &
    +12.d0*(3.d0*c606*z4x4+2.d0*c426*z4x2y2+c246*y4z4) &
    +8.d0*(4.d0*c804*z2x6+3.d0*c624*x4y2z2+2.d0*c444*y4z2x2+c264*y6z2) &
    +4.d0*(5.d0*cA02*x8+4.d0*c822*x6y2+3.d0*c642*x4y4 &
        +2.d0*c462*x2y6+c282*y8))*zx
if(n.le.6) then
    u=gmri*(1.d0+ri4*(p2+ri4*(p4+ri4*(p6 &
        +ri4*(p8+ri4*(pA+ri4*pC))))))
    r25i4=ri4*25.d0
    ur=-gmri3*(1.d0+ri4*(5.d0*p2+ri4*(9.d0*p4+ri4*(13.d0*p6 &
        +ri4*(17.d0*p8+ri4*(21.d0*pA+r25i4*pC))))))
    ux=ur*x+gmri5*(p2x+ri4*(p4x+ri4*(p6x+ri4*(p8x &
        +ri4*(pAx+ri4*pCx)))))
    uy=ur*y+gmri5*(p2y+ri4*(p4y+ri4*(p6y+ri4*(p8y &
        +ri4*(pAy+ri4*pCy)))))
    uz=ur*z+gmri5*(p2z+ri4*(p4z+ri4*(p6z+ri4*(p8z &
        +ri4*(pAz+ri4*pCz)))))
    urr=gmri5*(3.d0+ri4*(35.d0*p2+ri4*(99.d0*p4+ri4*(195.d0*p6 &
        +ri4*(323.d0*p8+ri4*(483.d0*pA+ri4*675.d0*pC))))))
    wx=5.d0*p2x+ri4*(9.d0*p4x+ri4*(13.d0*p6x+ri4*(17.d0*p8x &
        +ri4*(21.d0*pAx+r25i4*pCx))))
    wy=5.d0*p2y+ri4*(9.d0*p4y+ri4*(13.d0*p6y+ri4*(17.d0*p8y &
        +ri4*(21.d0*pAy+r25i4*pCy))))
    wz=5.d0*p2z+ri4*(9.d0*p4z+ri4*(13.d0*p6z+ri4*(17.d0*p8z &
        +ri4*(21.d0*pAz+r25i4*pCz))))
    uxx=-gm2ri7*wx*x
    uyy=-gm2ri7*wy*y
    uzz=-gm2ri7*wz*z
    uxy=-gmri7*(wx*y+wy*x)
    uyz=-gmri7*(wy*z+wz*y)
    uxz=-gmri7*(wz*x+wx*z)
    uxx=uxx+urr*x2+ur+gmri5*(p2xx+ri4*(p4xx+ri4*(p6xx+ri4*(p8xx &
        +ri4*(pAxx+ri4*pCxx)))))
    uyy=uyy+urr*y2+ur+gmri5*(p2yy+ri4*(p4yy+ri4*(p6yy+ri4*(p8yy &
        +ri4*(pAyy+ri4*pCyy)))))
    uzz=uzz+urr*z2+ur+gmri5*(p2zz+ri4*(p4zz+ri4*(p6zz+ri4*(p8zz &
        +ri4*(pAzz+ri4*pCzz)))))
    uxy=uxy+urr*xy+gmri9*(p4xy+ri4*(p6xy+ri4*(p8xy+ri4*(pAxy &
        +ri4*pCxy))))
    uyz=uyz+urr*yz+gmri9*(p4yz+ri4*(p6yz+ri4*(p8yz+ri4*(pAyz &
        +ri4*pCyz))))
    uxz=uxz+urr*zx+gmri9*(p4zx+ri4*(p6zx+ri4*(p8zx+ri4*(pAzx &
        +ri4*pCzx))))
    return
endif
!
xE=xC*x2;yE=yC*y2;zE=zC*z2
xCy2=xC*y2;yCz2=yC*z2;zCx2=zC*x2
xAy4=xA*y4;yAz4=yA*z4;zAx4=zA*x4
x8y6=x8*y6;y8z6=y8*z6;z8x6=z8*x6
x6y8=x6*y8;y6z8=y6*z8;z6x8=z6*x8
x4yA=x4*yA;y4zA=y4*zA;z4xA=z4*xA
x2yC=x2*yC;y2zC=y2*zC;z2xC=z2*xC
xAy2z2=xA*y2z2;yAz2x2=yA*z2x2;zAx2y2=zA*x2y2
x8y4z2=x8*y4z2;y8z4x2=y8*z4x2;z8x4y2=z8*x4y2
x8y2z4=x8*y2z4;y8z2x4=y8*z2x4;z8x2y4=z8*x2y4
x2y6z6=x2*y6z6;y2z6x6=y2*z6x6;z2x6y6=z2*x6y6
x6y4z4=x6*y4z4;y6z4x4=y6*z4x4;z6x4y4=z6*x4y4
pE=cE00*xE+c0E0*yE+c00E*zE &
    +cC20*xCy2+cC02*z2xC+c0C2*yCz2+c2C0*x2yC+c20C*zCx2+c02C*y2zC &
    +cA40*xAy4+cA22*xAy2z2+cA04*z4xA &
    +c0A4*yAz4+c2A2*yAz2x2+c4A0*x4yA &
    +c40A*zAx4+c22A*zAx2y2+c04A*y4zA &
    +c860*x8y6+c842*x8y4z2+c824*x8y2z4+c806*z6x8 &
    +c086*y8z6+c284*y8z4x2+c482*y8z2x4+c680*x6y8 &
    +c608*z8x6+c428*z8x4y2+c248*z8x2y4+c068*y6z8 &
    +c662*z2x6y6+c644*x6y4z4 &
    +c266*x2y6z6+c464*y6z4x4 &
    +c626*y2z6x6+c446*z6x4y4
aExC=14.d0*cE00*xC
aExA=12.d0*(cC20*xAy2+cC02*z2xA)
aEx8=10.d0*(cA40*x8y4+cA22*x8y2z2+cA04*z4x8)
aEx6=8.d0*(c860*x6y6+c842*x6y4z2+c824*x6y2z4+c806*z6x6)
aEx4=6.d0*(c680*x4y8+c662*y6z2x4+c644*x4y4z4+c626*z6x4y2+c608*z8x4)
aEx2=4.d0*(c4A0*x2yA+c482*y8z2x2+c464*y6z4x2+c446*z6x2y4+c428*z8x2y2 &
        +c40A*zAx2)
aEx0=2.d0*(c2C0*yC+c2A2*yAz2+c284*y8z4+c266*y6z6+c248*y4z8+c22A*y2zA &
        +c20C*zC)
aEyC=14.d0*c0E0*yC
aEyA=12.d0*(c0C2*yAz2+c2C0*x2yA)
aEy8=10.d0*(c0A4*y8z4+c2A2*y8z2x2+c4A0*x4y8)
aEy6=8.d0*(c086*y6z6+c284*y6z4x2+c482*y6z2x4+c680*x6y6)
aEy4=6.d0*(c068*y4z8+c266*z6x2y4+c464*x4y4z4+c662*x6y4z2+c860*x8y4)
aEy2=4.d0*(c04A*y2zA+c248*z8x2y2+c446*z6x4y2+c644*x6y2z4+c842*x8y2z2 &
        +cA40*xAy2)
aEy0=2.d0*(c02C*zC+c22A*zAx2+c428*z8x4+c626*z6x6+c824*z4x8 &
        +cA22*z2xA+cC20*xC)
aEzC=14.d0*c00E*zC
aEzA=12.d0*(c20C*zAx2+c02C*y2zA)
aEz8=10.d0*(c40A*z8x4+c22A*z8x2y2+c04A*y4z8)
aEz6=8.d0*(c608*z6x6+c428*z6x4y2+c248*z6x2y4+c068*y6z6)
aEz4=6.d0*(c806*z4x8+c626*x6y2z4+c446*x4y4z4+c266*y6z4x2+c086*y8z4)
aEz2=4.d0*(cA04*z2xA+c824*x8y2z2+c644*x6y4z2+c464*y6z2x4+c284*y8z2x2 &
        +c0A4*yAz2)
aEz0=2.d0*(cC02*xC+cA22*xAy2+c842*x8y4+c662*x6y6+c482*x4y8 &
        +c2A2*x2yA+c0C2*yC)
pEx=(aExC+aExA+aEx8+aEx6+aEx4+aEx2+aEx0)*x
pEy=(aEyC+aEyA+aEy8+aEy6+aEy4+aEy2+aEy0)*y
pEz=(aEzC+aEzA+aEz8+aEz6+aEz4+aEz2+aEz0)*z
pExx=13.d0*aExC+11.d0*aExA+9.d0*aEx8+7.d0*aEx6+5.d0*aEx4+3.d0*aEx2+aEx0
pEyy=13.d0*aEyC+11.d0*aEyA+9.d0*aEy8+7.d0*aEy6+5.d0*aEy4+3.d0*aEy2+aEy0
pEzz=13.d0*aEzC+11.d0*aEzA+9.d0*aEz8+7.d0*aEz6+5.d0*aEz4+3.d0*aEz2+aEz0
pExy=(24.d0*cC20*xA+20.d0*(2.d0*cA40*x8y2+cA22*z2x8) &
    +16.d0*(3.d0*c860*x6y4+2.d0*c842*x6y2z2+c824*z4x6) &
    +12.d0*(4.d0*c680*x4y6+3.d0*c662*z2x4y4+2.d0*c644*y2z4x4 &
        +c626*z6x4) &
    +8.d0*(5.d0*c4A0*x2y8+4.d0*c482*y6z2x2+3.d0*c464*x2y4z4 &
        +2.d0*c446*z6x2y2+c428*z8x2) &
    +4.d0*(6.d0*c2C0*yA+5.d0*c2A2*y8z2+4.d0*c284*y6z4+3.d0*c266*y4z6 &
        +2.d0*c248*y2z8+c22A*zA))*xy
pEyz=(24.d0*c0C2*yA+20.d0*(2.d0*c0A4*y8z2+c2A2*x2y8) &
    +16.d0*(3.d0*c086*y6z4+2.d0*c284*y6z2x2+c482*x4y6) &
    +12.d0*(4.d0*c068*y4z6+3.d0*c266*x2y4z4+2.d0*c464*z2x4y4 &
        +c662*x6y4) &
    +8.d0*(5.d0*c04A*y2z8+4.d0*c248*z6x2y2+3.d0*c446*y2z4x4 &
        +2.d0*c644*x6y2z2+c842*x8y2) &
    +4.d0*(6.d0*c02C*zA+5.d0*c22A*z8x2+4.d0*c428*z6x4+3.d0*c626*z4x6 &
        +2.d0*c824*z2x8+cA22*xA))*yz
pEzx=(24.d0*c20C*zA+20.d0*(2.d0*c40A*z8x2+c22A*y2z8) &
    +16.d0*(3.d0*c608*z6x4+2.d0*c428*z6x2y2+c248*y4z6) &
    +12.d0*(4.d0*c806*z4x6+3.d0*c626*y2z4x4+2.d0*c446*x2y4z4 &
        +c266*y6z4) &
    +8.d0*(5.d0*cA04*z2x8+4.d0*c824*x6y2z2+3.d0*c644*z2x4y4 &
        +2.d0*c464*y6z2x2+c284*y8z2) &
    +4.d0*(6.d0*cC02*xA+5.d0*cA22*x8y2+4.d0*c842*x6y4+3.d0*c662*x4y6 &
        +2.d0*c482*x2y8+c2A2*yA))*zx
if(n.le.7) then
    u=gmri*(1.d0+ri4*(p2+ri4*(p4+ri4*(p6 &
        +ri4*(p8+ri4*(pA+ri4*(pC+ri4*pE)))))))
    r29i4=ri4*29.d0
    ur=-gmri3*(1.d0+ri4*(5.d0*p2+ri4*(9.d0*p4+ri4*(13.d0*p6 &
        +ri4*(17.d0*p8+ri4*(21.d0*pA+ri4*(25.d0*pC &
        +r29i4*pE)))))))
    ux=ur*x+gmri5*(p2x+ri4*(p4x+ri4*(p6x+ri4*(p8x &
        +ri4*(pAx+ri4*(pCx+ri4*pEx))))))
    uy=ur*y+gmri5*(p2y+ri4*(p4y+ri4*(p6y+ri4*(p8y &
        +ri4*(pAy+ri4*(pCy+ri4*pEy))))))
    uz=ur*z+gmri5*(p2z+ri4*(p4z+ri4*(p6z+ri4*(p8z &
        +ri4*(pAz+ri4*(pCz+ri4*pEz))))))
    urr=gmri5*(3.d0+ri4*(35.d0*p2+ri4*(99.d0*p4+ri4*(195.d0*p6 &
        +ri4*(323.d0*p8+ri4*(483.d0*pA+ri4*(675.d0*pC &
        +ri4*899.d0*pE)))))))
    wx=5.d0*p2x+ri4*(9.d0*p4x+ri4*(13.d0*p6x+ri4*(17.d0*p8x &
        +ri4*(21.d0*pAx+ri4*(25.d0*pCx+r29i4*pEx)))))
    wy=5.d0*p2y+ri4*(9.d0*p4y+ri4*(13.d0*p6y+ri4*(17.d0*p8y &
        +ri4*(21.d0*pAy+ri4*(25.d0*pCy+r29i4*pEy)))))
    wz=5.d0*p2z+ri4*(9.d0*p4z+ri4*(13.d0*p6z+ri4*(17.d0*p8z &
        +ri4*(21.d0*pAz+ri4*(25.d0*pCz+r29i4*pEz)))))
    uxx=-gm2ri7*wx*x
    uyy=-gm2ri7*wy*y
    uzz=-gm2ri7*wz*z
    uxy=-gmri7*(wx*y+wy*x)
    uyz=-gmri7*(wy*z+wz*y)
    uxz=-gmri7*(wz*x+wx*z)
    uxx=uxx+urr*x2+ur+gmri5*(p2xx+ri4*(p4xx+ri4*(p6xx+ri4*(p8xx &
        +ri4*(pAxx+ri4*(pCxx+ri4*pExx))))))
    uyy=uyy+urr*y2+ur+gmri5*(p2yy+ri4*(p4yy+ri4*(p6yy+ri4*(p8yy &
        +ri4*(pAyy+ri4*(pCyy+ri4*pEyy))))))
    uzz=uzz+urr*z2+ur+gmri5*(p2zz+ri4*(p4zz+ri4*(p6zz+ri4*(p8zz &
        +ri4*(pAzz+ri4*(pCzz+ri4*pEzz))))))
    uxy=uxy+urr*xy+gmri9*(p4xy+ri4*(p6xy+ri4*(p8xy+ri4*(pAxy &
        +ri4*(pCxy+ri4*pExy)))))
    uyz=uyz+urr*yz+gmri9*(p4yz+ri4*(p6yz+ri4*(p8yz+ri4*(pAyz &
        +ri4*(pCyz+ri4*pEyz)))))
    uxz=uxz+urr*zx+gmri9*(p4zx+ri4*(p6zx+ri4*(p8zx+ri4*(pAzx &
        +ri4*(pCzx+ri4*pEzx)))))
    return
endif
!
xG=xE*x2;yG=yE*y2;zG=zE*z2
xEy2=xE*y2;yEz2=yE*z2;zEx2=zE*x2
xCy4=xC*y4;yCz4=yC*z4;zCx4=zC*x4
xAy6=xA*y6;yAz6=yA*z6;zAx6=zA*x6
x8y8=x8*y8;y8z8=y8*z8;z8x8=z8*x8
x6yA=x6*yA;y6zA=y6*zA;z6xA=z6*xA
x4yC=x4*yC;y4zC=y4*zC;z4xC=z4*xC
x2yE=x2*yE;y2zE=y2*zE;z2xE=z2*xE
xCy2z2=xC*y2z2;yCz2x2=yC*z2x2;zCx2y2=zC*x2y2
xAy4z2=xA*y4z2;yAz4x2=yA*z4x2;zAx4y2=zA*x4y2
xAy2z4=xA*y2z4;yAz2x4=yA*z2x4;zAx2y4=zA*x2y4
x8y6z2=x8*y6z2;y8z6x2=y8*z6x2;z8x6y2=z8*x6y2
x8y4z4=x8*y4z4;y8z4x4=y8*z4x4;z8x4y4=z8*x4y4
x8y2z6=x8*y2z6;y8z2x6=y8*z2x6;z8x2y6=z8*x2y6
x4y6z6=x4*y6z6;y4z6x6=y4*z6x6;z4x6y6=z4*x6y6
pG=cG00*xG+c0G0*yG+c00G*zG &
    +cE20*xEy2+cE02*z2xE+c0E2*yEz2+c2E0*x2yE+c20E*zEx2+c02E*y2zE &
    +cC40*xCy4+cC22*xCy2z2+cC04*z4xC &
    +c0C4*yCz4+c2C2*yCz2x2+c4C0*x4yC &
    +c40C*zCx4+c22C*zCx2y2+c04C*y4zC &
    +cA60*xAy6+cA42*xAy4z2+cA24*xAy2z4+cA06*z6xA &
    +c0A6*yAz6+c2A4*yAz4x2+c4A2*yAz2x4+c6A0*x6yA &
    +c60A*zAx6+c42A*zAx4y2+c24A*zAx2y4+c06A*y6zA &
    +c880*x8y8+c862*x8y6z2+c844*x8y4z4+c826*x8y2z6 &
    +c088*y8z8+c286*y8z6x2+c484*y8z4x4+c682*y8z2x6 &
    +c808*z8x8+c628*z8x6y2+c448*z8x4y4+c268*z8x2y6 &
    +c466*x4y6z6+c646*y4z6x6+c664*z4x6y6
aGxE=16.d0*cG00*xE
aGxC=14.d0*(cE20*xCy2+cE02*z2xC)
aGxA=12.d0*(cC40*xAy4+cC22*xAy2z2+cC04*z4xA)
aGx8=10.d0*(cA60*x8y6+cA42*x8y4z2+cA24*x8y2z4+cA06*z6x8)
aGx6=8.d0*(c880*x6y8+c862*z2x6y6+c844*x6y4z4+c826*y2z6x6+c808*z8x6)
aGx4=6.d0*(c6A0*x4yA+c682*y8z2x4+c664*y6z4x4+c646*z6x4y4+c628*z8x4y2 &
        +c60A*zAx4)
aGx2=4.d0*(c4C0*x2yC+c4A2*yAz2x2+c484*y8z4x2+c466*x2y6z6+c448*z8x2y4 &
        +c42A*zAx2y2+c40C*zCx2)
aGx0=2.d0*(c2E0*yE+c2C2*yCz2+c2A4*yAz4+c286*y8z6+c268*y6z8 &
        +c24A*y4zA+c22C*y2zC+c20E*zE)
aGyE=16.d0*c0G0*yE
aGyC=14.d0*(c0E2*yCz2+c2E0*x2yC)
aGyA=12.d0*(c0C4*yAz4+c2C2*yAz2x2+c4C0*x4yA)
aGy8=10.d0*(c0A6*y8z6+c2A4*y8z4x2+c4A2*y8z2x4+c6A0*x6y8)
aGy6=8.d0*(c088*y6z8+c286*x2y6z6+c484*y6z4x4+c682*z2x6y6+c880*x8y6)
aGy4=6.d0*(c06A*y4zA+c268*z8x2y4+c466*z6x4y4+c664*x6y4z4+c862*x8y4z2 &
        +cA60*xAy4)
aGy2=4.d0*(c04C*y2zC+c24A*zAx2y2+c448*z8x4y2+c646*y2z6x6+c844*x8y2z4 &
        +cA42*xAy2z2+cC40*xCy2)
aGy0=2.d0*(c02E*zE+c22C*zCx2+c42A*zAx4+c628*z8x6+c826*z6x8 &
        +cA24*z4xA+cC22*z2xC+cE20*xE)
aGzE=16.d0*c00G*zE
aGzC=14.d0*(c20E*zCx2+c02E*y2zC)
aGzA=12.d0*(c40C*zAx4+c22C*zAx2y2+c04C*y4zA)
aGz8=10.d0*(c60A*z8x6+c42A*z8x4y2+c24A*z8x2y4+c06A*y6z8)
aGz6=8.d0*(c808*z6x8+c628*y2z6x6+c448*z6x4y4+c268*x2y6z6+c088*y8z6)
aGz4=6.d0*(cA06*z4xA+c826*x8y2z4+c646*x6y4z4+c466*y6z4x4+c286*y8z4x2 &
        +c0A6*yAz4)
aGz2=4.d0*(cC04*z2xC+cA24*xAy2z2+c844*x8y4z2+c664*z2x6y6+c484*y8z2x4 &
        +c2A4*yAz2x2+c0C4*yCz2)
aGz0=2.d0*(cE02*xE+cC22*xCy2+cA42*xAy4+c862*x8y6+c682*x6y8 &
        +c4A2*x4yA+c2C2*x2yC+c0E2*yE)
pGx=(aGxE+aGxC+aGxA+aGx8+aGx6+aGx4+aGx2+aGx0)*x
pGy=(aGyE+aGyC+aGyA+aGy8+aGy6+aGy4+aGy2+aGy0)*y
pGz=(aGzE+aGzC+aGzA+aGz8+aGz6+aGz4+aGz2+aGz0)*z
pGxx=15.d0*aGxE+13.d0*aGxC+11.d0*aGxA+9.d0*aGx8+7.d0*aGx6+5.d0*aGx4 &
    +3.d0*aGx2+aGx0
pGyy=15.d0*aGyE+13.d0*aGyC+11.d0*aGyA+9.d0*aGy8+7.d0*aGy6+5.d0*aGy4 &
    +3.d0*aGy2+aGy0
pGzz=15.d0*aGzE+13.d0*aGzC+11.d0*aGzA+9.d0*aGz8+7.d0*aGz6+5.d0*aGz4 &
    +3.d0*aGz2+aGz0
pGxy=(28.d0*cE20*xC+24.d0*(2.d0*cC40*xAy2+cC22*z2xA) &
    +20.d0*(3.d0*cA60*x8y4+2.d0*cA42*x8y2z2+cA24*z4x8) &
    +16.d0*(4.d0*c880*x6y6+3.d0*c862*x6y4z2+2.d0*c844*x6y2z4 &
        +c826*z6x6) &
    +12.d0*(5.d0*c6A0*x4y8+4.d0*c682*y6z2x4+3.d0*c664*x4y4z4 &
        +2.d0*c646*z6x4y2+c628*z8x4) &
    +8.d0*(6.d0*c4C0*x2yA+5.d0*c4A2*y8z2x2+4.d0*c484*y6z4x2 &
        +3.d0*c466*z6x2y4+2.d0*c448*z8x2y2+c42A*zAx2) &
    +4.d0*(7.d0*c2E0*yC+6.d0*c2C2*yAz2+5.d0*c2A4*y8z4+4.d0*c286*y6z6 &
        +3.d0*c268*y4z8+2.d0*c24A*y2zA+c22C*zC))*xy
pGyz=(28.d0*c0E2*yC+24.d0*(2.d0*c0C4*yAz2+c2C2*x2yA) &
    +20.d0*(3.d0*c0A6*y8z4+2.d0*c2A4*y8z2x2+c4A2*x4y8) &
    +16.d0*(4.d0*c088*y6z6+3.d0*c286*y6z4x2+2.d0*c484*y6z2x4 &
        +c682*x6y6) &
    +12.d0*(5.d0*c06A*y4z8+4.d0*c268*z6x2y4+3.d0*c466*x4y4z4 &
        +2.d0*c664*x6y4z2+c862*x8y4) &
    +8.d0*(6.d0*c04C*y2zA+5.d0*c24A*z8x2y2+4.d0*c448*x6y4z2 &
        +3.d0*c646*x6y2z4+2.d0*c844*x8y2z2+cA42*xAy2) &
    +4.d0*(7.d0*c02E*zC+6.d0*c22C*zAx2+5.d0*c42A*z8x4+4.d0*c628*z6x6 &
        +3.d0*c826*z4x8+2.d0*cA24*z2xA+cC22*xC))*yz
pGzx=(28.d0*c20E*zC+24.d0*(2.d0*c40C*zAx2+c22C*y2zA) &
    +20.d0*(3.d0*c60A*z8x4+2.d0*c42A*z8x2y2+c24A*y4z8) &
    +16.d0*(4.d0*c808*z6x6+3.d0*c628*z6x4y2+2.d0*c448*z6x2y4 &
        +c268*y6z6) &
    +12.d0*(5.d0*cA06*z4x8+4.d0*c826*x6y2z4+3.d0*c646*x4y4z4 &
        +2.d0*c466*y6z4x2+c286*y8z4) &
    +8.d0*(6.d0*cC04*z2xA+5.d0*cA24*x8y2z2+4.d0*c844*y6z4x2 &
        +3.d0*c664*y6z2x4+2.d0*c484*y8z2x2+c2A4*yAz2) &
    +4.d0*(7.d0*cE02*xC+6.d0*cC22*xAy2+5.d0*cA42*x8y4+4.d0*c862*x6y6 &
        +3.d0*c682*x4y8+2.d0*c4A2*x2yA+c2C2*yC))*zx
u=gmri*(1.d0+ri4*(p2+ri4*(p4+ri4*(p6 &
    +ri4*(p8+ri4*(pA+ri4*(pC+ri4*(pE+ri4*pG))))))))
r33i4=ri4*33.d0
ur=-gmri3*(1.d0+ri4*(5.d0*p2+ri4*(9.d0*p4+ri4*(13.d0*p6 &
    +ri4*(17.d0*p8+ri4*(21.d0*pA+ri4*(25.d0*pC &
    +ri4*(29.d0*pE+r33i4*pG))))))))
ux=ur*x+gmri5*(p2x+ri4*(p4x+ri4*(p6x+ri4*(p8x &
    +ri4*(pAx+ri4*(pCx+ri4*(pEx+ri4*pGx)))))))
uy=ur*y+gmri5*(p2y+ri4*(p4y+ri4*(p6y+ri4*(p8y &
    +ri4*(pAy+ri4*(pCy+ri4*(pEy+ri4*pGy)))))))
uz=ur*z+gmri5*(p2z+ri4*(p4z+ri4*(p6z+ri4*(p8z &
    +ri4*(pAz+ri4*(pCz+ri4*(pEz+ri4*pGz)))))))
urr=gmri5*(3.d0+ri4*(35.d0*p2+ri4*(99.d0*p4+ri4*(195.d0*p6 &
    +ri4*(323.d0*p8+ri4*(483.d0*pA+ri4*(675.d0*pC &
    +ri4*(899.d0*pE+ri4*1155.d0*pG))))))))
wx=5.d0*p2x+ri4*(9.d0*p4x+ri4*(13.d0*p6x+ri4*(17.d0*p8x &
    +ri4*(21.d0*pAx+ri4*(25.d0*pCx+ri4*(29.d0*pEx+r33i4*pGx))))))
wy=5.d0*p2y+ri4*(9.d0*p4y+ri4*(13.d0*p6y+ri4*(17.d0*p8y &
    +ri4*(21.d0*pAy+ri4*(25.d0*pCy+ri4*(29.d0*pEy+r33i4*pGy))))))
wz=5.d0*p2z+ri4*(9.d0*p4z+ri4*(13.d0*p6z+ri4*(17.d0*p8z &
    +ri4*(21.d0*pAz+ri4*(25.d0*pCz+ri4*(29.d0*pEz+r33i4*pGz))))))
uxx=-gm2ri7*wx*x
uyy=-gm2ri7*wy*y
uzz=-gm2ri7*wz*z
uxy=-gmri7*(wx*y+wy*x)
uyz=-gmri7*(wy*z+wz*y)
uxz=-gmri7*(wz*x+wx*z)
uxx=uxx+urr*x2+ur+gmri5*(p2xx+ri4*(p4xx+ri4*(p6xx+ri4*(p8xx &
    +ri4*(pAxx+ri4*(pCxx+ri4*(pExx+ri4*pGxx)))))))
uyy=uyy+urr*y2+ur+gmri5*(p2yy+ri4*(p4yy+ri4*(p6yy+ri4*(p8yy &
    +ri4*(pAyy+ri4*(pCyy+ri4*(pEyy+ri4*pGyy)))))))
uzz=uzz+urr*z2+ur+gmri5*(p2zz+ri4*(p4zz+ri4*(p6zz+ri4*(p8zz &
    +ri4*(pAzz+ri4*(pCzz+ri4*(pEzz+ri4*pGzz)))))))
uxy=uxy+urr*xy+gmri9*(p4xy+ri4*(p6xy+ri4*(p8xy+ri4*(pAxy &
    +ri4*(pCxy+ri4*(pExy+ri4*pGxy))))))
uyz=uyz+urr*yz+gmri9*(p4yz+ri4*(p6yz+ri4*(p8yz+ri4*(pAyz &
    +ri4*(pCyz+ri4*(pEyz+ri4*pGyz))))))
uxz=uxz+urr*zx+gmri9*(p4zx+ri4*(p6zx+ri4*(p8zx+ri4*(pAzx &
    +ri4*(pCzx+ri4*(pEzx+ri4*pGzx))))))
!
return;end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              a,b,c=    2.500000000000000E+00    3.100000000000000E+00    5.000000000000000E+00
!                 GM=    3.100000000000000E+02
!     c200,c020,c002=   -3.685000000000000E+00   -2.004999999999999E+00    5.689999999999999E+00
!     c400,c040,c004=    9.530574166666671E+00   -1.880822583333333E+01    7.694636666666657E+00
!     c022,c202,c220=    6.193249000000004E+01   -1.081003100000000E+02    5.091686500000000E+01
!     c600,c060,c006=    4.089895715327382E+02    6.889313945327383E+02   -9.765010144047624E+02
!     c420,c042,c204=   -1.236868576651786E+03   -1.096373973633929E+04    4.290875238035719E+03
!     c240,c024,c402=    6.297688183482146E+02    1.035663997803572E+04   -4.897974996339283E+03
!               c222=    3.642598549821439E+03
!     c800,c080,c008=   -1.351747240312123E+04   -8.307323821056232E+03   -1.649574373881459E+04
!     c620,c062,c206=    3.424164602976201E+04    2.648017758058222E+05    2.374441376463040E+05
!     c260,c026,c602=   -3.219670881624805E+04    2.244366870405043E+05    3.442475812576324E+05
!     c044,c404,c440=   -6.837321612484663E+05   -7.162507877629654E+05    5.876405030671020E+04
!     c422,c242,c224=   -8.662089922866907E+05    1.303663304034592E+05    7.358426618832316E+05
!     cA00,c0A0,c00A=    2.274999576968344E+05   -1.314618333577888E+05    1.505982219689371E+05
!     c820,c082,c208=   -8.316326357409207E+05    5.118548016974139E+06    1.242527605381609E+06
!     c280,c028,c802=    7.972344841263462E+05   -8.019447593983772E+06   -9.405865460616628E+06
!     c640,c064,c406=   -4.888371837589480E+06   -2.037685890106719E+07   -2.125628962208308E+07
!     c460,c046,c604=   -2.107290811108435E+05    2.196626130828875E+07    3.512471467849720E+07
!     c622,c262,c226=    5.261594482628252E+07   -2.105819106887266E+07    9.274696478181395E+07
!     c244,c424,c442=   -2.384103610832330E+07   -2.080263758462113E+08    7.648651378050473E+07
!     cC00,c0C0,c00C=   -2.671376500927814E+05    8.066021011084014E+06    1.023082822600369E+07
!     cA20,c0A2,c20A=   -4.576713838034065E+06   -5.110771571701628E+08   -1.882091141645356E+08
!     c2A0,c02A,cA02=   -2.128022956138252E+07   -4.870255487517104E+08    2.220779874415722E+07
!     c840,c084,c408=    1.754080963587964E+08    3.572348092188994E+09    9.324224991241341E+08
!     c480,c048,c804=   -1.011288648768552E+08    3.173545758527917E+09   -2.547574800763699E+07
!     c822,c282,c228=   -8.464964554412479E+08    1.564383519523342E+09    2.874875142659283E+09
!     c066,c606,c660=   -6.109876220319574E+09   -6.601922606925120E+08    1.169640909664074E+08
!     c642,c264,c426=   -6.665888062542413E+09   -8.377603276498122E+09   -1.620494606508818E+10
!     c462,c246,c624=    1.077146852055824E+09    2.788862066011534E+09    1.061620485460155E+10
!               c444=    1.397185302621633E+10
!     cE00,c0E0,c00E=   -1.290377622552024E+08   -1.267119206264848E+08    6.588658268822457E+07
!     cC20,c0C2,c20C=    1.743648346626192E+09    1.134477060781146E+10   -4.181507079707178E+09
!     c2C0,c02C,cC02=    1.860141691987171E+08   -1.814171944921460E+09    9.998788018597248E+09
!     cA40,c0A4,c40A=   -2.852272763108765E+09   -1.164034737866744E+11    7.558056258307852E+10
!     c4A0,c04A,cA04=    6.342847038065727E+09    4.953987610043535E+10   -9.365880915479024E+10
!     cA22,c2A2,c22A=   -9.796715429867648E+10   -5.033401739550996E+10   -1.775039082377903E+11
!     c860,c086,c608=   -9.352172319455513E+09    2.983912515313831E+11   -2.394647972948075E+11
!     c680,c068,c806=    6.623620016687537E+09   -2.246484728746308E+11    2.499018409238138E+11
!     c842,c284,c428=    2.686348591317279E+11    7.622875474295964E+11    1.908466431836058E+11
!     c482,c248,c824=   -3.847824169632701E+11    1.140432668599846E+12    4.661187981083454E+11
!     c266,c626,c662=   -2.064797802389170E+12   -2.922372216121761E+11    7.639946447750406E+10
!     c644,c464,c446=   -1.444628003808483E+12    1.604652617968178E+12   -1.600246141596890E+11
!     cG00,c0G0,c00G=    4.870406925446485E+09   -2.330274513975320E+09   -5.325118180931003E+09
!     cE20,c0E2,c20E=   -9.291989921145975E+10    2.503775869142189E+11    9.986353660300833E+10
!     c2E0,c02E,cE02=    2.925535476282125E+10    5.391506451086981E+11   -4.915289318421187E+11
!     cC40,c0C4,c40C=   -7.457031422502141E+10   -3.503023922530839E+12    8.755409866248896E+11
!     c4C0,c04C,cC04=   -1.493367349012375E+11   -5.786980159044739E+12    5.971000014006618E+12
!     cC22,c2C2,c22C=    8.903132713592949E+12   -1.766216874009364E+12   -1.434082775062313E+13
!     cA60,c0A6,c60A=    5.470495271311066E+11    1.536849710451325E+13   -9.678356705660154E+12
!     c6A0,c06A,cA06=   -5.933359194191797E+11    2.077208204718378E+13   -1.996237621620191E+13
!     cA42,c2A4,c42A=   -3.284102168115132E+12    6.721223193329869E+11    8.738964546765830E+13
!     c4A2,c24A,cA24=    1.875626329476941E+13    7.035945978919550E+13   -9.465035768140719E+13
!     c088,c808,c880=   -2.756342139232163E+13    2.319944649744181E+13    2.711107208062259E+10
!     c862,c286,c628=   -2.537633873915747E+13    8.019342928191161E+13   -2.140584501736661E+14
!     c682,c268,c826=    2.594100635560549E+13   -1.629678931382746E+14    2.487224278007198E+14
!     c844,c484,c448=    8.807161310875780E+13   -2.055244905997848E+14   -1.202762155732662E+14
!     c466,c646,c664=    3.862808313296832E+14   -1.617652289262518E+14   -2.635115543425295E+12
!              X,Y,Z=    3.000000000000000E+01    4.000000000000000E+01    5.000000000000000E+01
!       order,V=    0    4.384062043356595E+00
!       order,V=    2    4.385412422147190E+00
!       order,V=    4    4.385413021793824E+00
!       order,V=    6    4.385413035385058E+00
!       order,V=    8    4.385413035403303E+00
!       order,V=   10    4.385413035403428E+00
!       order,V=   12    4.385413035403428E+00
!       order,V=   14    4.385413035403428E+00
!       order,V=   16    4.385413035403428E+00
!       order,g=    0   -2.630437226013957E-02   -3.507249634685276E-02   -4.384062043356595E-02
!       order,g=    2   -2.638365626856886E-02   -3.515463964054672E-02   -4.380835812098889E-02
!       order,g=    4   -2.638376088287478E-02   -3.515460398892207E-02   -4.380838383836853E-02
!       order,g=    6   -2.638376163337641E-02   -3.515460559729471E-02   -4.380838400414223E-02
!       order,g=    8   -2.638376163528297E-02   -3.515460560151899E-02   -4.380838400290287E-02
!       order,g=   10   -2.638376163531891E-02   -3.515460560153799E-02   -4.380838400289384E-02
!       order,g=   12   -2.638376163531892E-02   -3.515460560153809E-02   -4.380838400289370E-02
!       order,g=   14   -2.638376163531892E-02   -3.515460560153809E-02   -4.380838400289370E-02
!       order,g=   16   -2.638376163531892E-02   -3.515460560153809E-02   -4.380838400289370E-02
!   order,Gamma=    0   -4.033337079888068E-04    6.313049342433497E-04    7.891311678041872E-04   -3.507249634685276E-05    1.052174890405583E-03    4.384062043356596E-04
!   order,Gamma=    2   -4.019486723109468E-04    6.359683206621120E-04    7.909121579368044E-04   -3.185097883265554E-05    1.053371108038485E-03    4.337996511436029E-04
!   order,Gamma=    4   -4.019417397721261E-04    6.359738287232826E-04    7.909161456812855E-04   -3.185932785888587E-05    1.053370204227808E-03    4.338010676310122E-04
!   order,Gamma=    6   -4.019417441822762E-04    6.359739509184590E-04    7.909161706514961E-04   -3.185924244028891E-05    1.053370319915448E-03    4.338009866225661E-04
!   order,Gamma=    8   -4.019417444595240E-04    6.359739514886640E-04    7.909161707429945E-04   -3.185924184138461E-05    1.053370319939057E-03    4.338009863009094E-04
!   order,Gamma=   10   -4.019417444564358E-04    6.359739514938600E-04    7.909161707456098E-04   -3.185924184028444E-05    1.053370319939620E-03    4.338009862967208E-04
!   order,Gamma=   12   -4.019417444564295E-04    6.359739514938835E-04    7.909161707455901E-04   -3.185924184026720E-05    1.053370319939618E-03    4.338009862966971E-04
!   order,Gamma=   14   -4.019417444564296E-04    6.359739514938836E-04    7.909161707455900E-04   -3.185924184026688E-05    1.053370319939618E-03    4.338009862966968E-04
!   order,Gamma=   16   -4.019417444564296E-04    6.359739514938836E-04    7.909161707455900E-04   -3.185924184026688E-05    1.053370319939618E-03    4.338009862966968E-04
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
