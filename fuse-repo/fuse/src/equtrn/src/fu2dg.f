      program fu2dg
c
c  version : 13.08.2002 15:57
c
c=======================================================================
c*** Translation of equilibrium data from Naka (Fujieda Excel format)
c*** to dg-compatible format
c=======================================================================
      parameter (ngpr=257, ngpz=257)
      real*8 pfm(ngpr,ngpz),rgr(ngpr),zgr(ngpz)
      real*8 rmin,zmin,rmax,zmax,delr,delz,psilim,btorc,rcntc,fg,rbtr
      character*256 filename
      parameter (rbtr=5.3*6.2)
c=======================================================================

      if(iargc().lt.2) then !{
	write(*,*) '1st arg == input dg equilibrium file'
	write(*,*) '2nd arg == output dg equilibrium file'
	write(*,*) '3rd arg == R*Btor (optional, default',rbtr,')'
	stop
      endif !}

      call getarg(1,filename)
      open(1,file=filename)

      btorc=rbtr
      rcntc=1
      call rdeqfu(1,ngpr,ngpz,iret, nr,nz,
     ,           rmin,zmin,rmax,zmax,delr,delz,psilim,
     ,           btorc,rcntc,fg,pfm,rgr,zgr)
      if(iret.ne.0) then !{
          print *,'==== fu2dg: error in rdeqfu. iret =',iret
          stop
      end if !}

      call getarg(2,filename)
      open(2,file=filename)
      if(iargc().gt.2) then !{
        call getarg(3,filename)
        read(filename,*) btorc
        rcntc=1.
      end if !}

      print *,'psilim, btorc, rcntc = ',psilim, btorc, rcntc
      call wreqdg(2,ngpr,iret, nr,nz, psilim,btorc,rcntc,rgr,zgr,pfm)
      if(iret.ne.0) then !{
          print *,'==== fu2dg: error in wreqdg. iret = ',iret
      end if !}
c=======================================================================
      end
