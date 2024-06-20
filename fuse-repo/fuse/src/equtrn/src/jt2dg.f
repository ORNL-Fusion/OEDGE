      program jt2dg
c
c  version : 09.07.97 17:14
c
c=====================================================
c*** Translation of efit equilibrium data into dg-compatible format
c***
c*** The input and output files must be pre-connected to the units
c*** fort.1 and fort.2, and the field data to the fort.3
c=====================================================
      parameter (ngpr=257, ngpz=257)
      real*8 fg(ngpr),pg(ngpr),ffg(ngpr),ppg(ngpr)
      real*8 pfm(ngpr,ngpz),rgr(ngpr),zgr(ngpz)
      real*8 rdim,zdim,rcntc,redge,zmsmid,rma,zma,psimin,psilim,btorc
      character title*40, date*8
c=====================================================
c
      call rdeqdg(1,ngpr,ngpz,iret, nr,nz,
     ,                                  psilim,btorc,rcntc,rgr,zgr,pfm)
      if(iret.ne.0) then
          print *,'==== jt2dg: error in rdefit. iret =',iret
          stop
      end if
c
      do i=1,nr
        do j=1,nz
          pfm(i,j)=pfm(i,j)-psilim
        end do
      end do
      print *,'psilim = ',psilim
      call wreqdg(2,ngpr,iret, nr,nz, psilim,btorc,rcntc,rgr,zgr,pfm)
      if(iret.ne.0) then
          print *,'==== ef2dg: error in wreqdg. iret = ',iret
      end if
c
      end
