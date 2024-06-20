      program vr2dg
c
c  version : 16.11.95 20:04
c
c=====================================================
c*** Translation of TdeV equilibrium data into dg-compatible format
c***
c*** The input and output files must be pre-connected to the units
c*** fort.1 and fort.2
c=====================================================
      parameter (ngpr=257, ngpz=257)
      real*8 pfm(ngpr,ngpz),rgr(ngpr),zgr(ngpz)
      real*8 rcntc,psilim,btorc
c=====================================================
c
      call rdeqvr(1,ngpr,ngpz,iret,nr,nz,btorc,rcntc,rgr,zgr,pfm)
      if(iret.ne.0) then
          print *,'==== vr2dg: error in rdeqvr. iret =',iret
          stop
      end if
c
      psilim=0.
      print *,'psilim = ',psilim
      call wreqdg(2,ngpr,iret, nr,nz, psilim,btorc,rcntc,rgr,zgr,pfm)
      if(iret.ne.0) then
          print *,'==== vr2dg: error in wreqdg. iret = ',iret
      end if
c
      end
