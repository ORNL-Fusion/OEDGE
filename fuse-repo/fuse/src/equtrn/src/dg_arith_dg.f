      program dg_arith_dg
c
c  version : 30.01.96 21:58
c
c=====================================================
c*** Shift of the magnetic flux by a specified increment
c***
c*** The input and output files must be pre-connected to the units
c*** fort.1 and fort.2; the increment is taken from the standard input.
c=====================================================
      parameter (ngpr=257, ngpz=257)
      real*8 pfm(ngpr,ngpz),rgr(ngpr),zgr(ngpz)
      real*8 rcntc,psilim,btorc,shift
      character*256 filename
c=====================================================
c
      if(iargc().ne.2) then
	write(*,*) '1st arg == input dg equilibrium file'
	write(*,*) '2nd arg == output dg equilibrium file'
	stop
      endif

      call getarg(1,filename)
      open(1,file=filename)
      call rdeqdg(1,ngpr,ngpz,iret,nr,nz,btorc,rcntc,rgr,zgr,pfm)
      if(iret.ne.0) then
          print *,'==== dg_arith_dg: error in rdeqdg. iret =',iret
          stop
      end if
c
      write(*,*) 'alpha, beta'
      read(*,*) alpha,beta
      do j=1,nz
          do i=1,nr
              pfm(i,j)=pfm(i,j)+alpha*zgr(j)+beta
          end do
      end do
      psilim=0.
      call getarg(2,filename)
      open(2,file=filename)
      call wreqdg(2,ngpr,iret,nr,nz,psilim,btorc,rcntc,rgr,zgr,pfm)
      if(iret.ne.0) then
          print *,'==== dg_arith_dg: error in wreqdg. iret = ',iret
      end if
c
      end
