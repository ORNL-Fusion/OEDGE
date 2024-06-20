      program nk2dg
c
c  version : 16.01.96 19:48
c
c=====================================================
c*** Translation of equilibrium data from Naka to dg-compatible format
c***
c*** The input and output files must be pre-connected to the units
c*** fort.1 and fort.2
c=====================================================
      parameter (ngpr=257, ngpz=257)
      real*8 pfm(ngpr,ngpz),rgr(ngpr),zgr(ngpz)
      real*8 rmin,zmin,rmax,zmax,delr,delz,psilim,btorc,rcntc,fg
c=====================================================
c
      call rdeqnk(1,ngpr,ngpz,iret, nr,nz,
     ,           rmin,zmin,rmax,zmax,delr,delz,psilim,
     ,           btorc,rcntc,fg,pfm,rgr,zgr)
      if(iret.ne.0) then
          print *,'==== nk2dg: error in rdeqnk. iret =',iret
          stop
      end if
c
      print *,'psilim = ',psilim
      call wreqdg(2,ngpr,iret, nr,nz, psilim,btorc,rcntc,rgr,zgr,pfm)
      if(iret.ne.0) then
          print *,'==== nk2dg: error in wreqdg. iret = ',iret
      end if
c
      end
