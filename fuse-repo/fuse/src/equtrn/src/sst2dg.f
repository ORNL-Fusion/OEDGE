      program sst2dg
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
	write(*,*) '1st arg == sst eq file'
	write(*,*) '2nd arg == output dg equilibrium file'
	stop
      endif
      call getarg(1,filename)
      open(1,file=filename)
      ngr=65
      ngz=65
      zfac=-0.005
      do i=1,ngr
	do j=1,ngz
	  read(1,*) gpr(i),gpz(j),pfm(i,j)
	  pfm(i,j)=pfm(i,j)+gpz(j)*zfac
	enddo
      enddo
      psilim=0.0634-0.355983*zfac
      pfm(1:ngr,1:ngz)=pfm(1:ngr,1:ngz)-psilim
      btorc=3
      rcntc=1.1
c
      print *,'psilim = ',psilim
      call getarg(2,filename)
      open(2,file=filename)
      call wreqdg(2,ngpr,iret,ngr,ngz,psilim,
     1 btorc,rcntc,gpr,gpz,pfm)
      if(iret.ne.0) then
          print *,'==== ef2dg: error in wreqdg. iret = ',iret
      end if
c
      end
