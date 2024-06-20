      program ef2dg
c
c  version : 16.11.95 20:04
c
c=====================================================
c*** Translation of efit equilibrium data into dg-compatible format
c***
c*** The input and output files must be pre-connected to the units
c*** fort.1 and fort.2, and the field data to the fort.3
c=====================================================
      parameter (pi=3.14159 26535 89793)
      parameter (ngpr=257, ngpz=257, ngpf=257)
      real*8 gpr(ngpr),gpz(ngpz),
     1 pfl(ngpf),pc(ngpf),pcd(ngpf),pr(ngpf),prd(ngpf),q(ngpf),
     2 apf(3),xpf(3)
      real*8 pfm(ngpr,ngpz)
      real*8 rcntc,btorc
      real*8 Bt,derivative(6)
      character*256 filename
c=====================================================
c
      if(iargc().ne.2) then
	write(*,*) '1st arg == Diva file'
	write(*,*) '2nd arg == output dg equilibrium file'
	stop
      endif
      call getarg(1,filename)
      
      call open_diva(ifail,filename)
      if(ifail.ne.0) then
	write(*,*) 'open_diva: ifail = ',ifail
	stop
      endif

      call dimension_diva(ifail, ngr, ngz, ngf)
      write(*,*) 'dimension_diva: ngr,ngz,ngf = ',ngr,ngz,ngf
      if(ifail.ne.0) then
	write(*,*) 'dimension_diva: ifail = ',ifail
	stop
      endif
      if(ngr.gt.ngpr) write(*,*) 'Increase ngpr to at least ',ngr
      if(ngz.gt.ngpz) write(*,*) 'Increase ngpz to at least ',ngz
      if(ngf.gt.ngpf) write(*,*) 'Increase ngpf to at least ',ngf

      call profile_diva(
     1 ifail,
     2 gpr,ngpr,
     3 gpz,ngpz,
     4 pfl,pc,pcd,pr,prd,q,ngpf,
     5 apf,xpf)
      if(ifail.ne.0) then
	write(*,*) 'profile_diva: ifail = ',ifail
	stop
      endif
      write(*,*) 'profile_diva: apf = ',apf
      write(*,*) 'profile_diva: xpf = ',xpf
      psilim=xpf(3)/(2.0*pi)

      call flux_diva(ifail,apf(1),apf(2),derivative,Bt)
      if(ifail.ne.0) then
	write(*,*) 'flux_diva: ifail = ',ifail
	stop
      endif
      btorc=bt
      rcntc=apf(1)
      write(*,*) 'derivative @ axis = ',derivative

      do j=1,ngz
	do i=1,ngr
	  call flux_diva(ifail,gpr(i),gpz(j),derivative,Bt)
	  if(ifail.ne.0) then
	    write(*,*) 'flux_diva: ifail = ',ifail
	    stop
	  endif
	  pfm(i,j)=derivative(1)/(2.0*pi)
	enddo
      enddo

      call close_diva()
c
      print *,'psilim = ',psilim
      call getarg(2,filename)
      open(2,file=filename)
      psilim=0
      call wreqdg(2,ngpr,iret,ngr,ngz,psilim,
     1 btorc,rcntc,gpr,gpz,pfm)
      if(iret.ne.0) then
          print *,'==== ef2dg: error in wreqdg. iret = ',iret
      end if
c
      end
