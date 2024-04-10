      program shotfile2dg
c
c  version : 16.11.95 20:04
c
c=====================================================
c*** Translation of shotfile equilibrium data into dg-compatible format
c=====================================================
      parameter (pi=3.14159 26535 89793)
      parameter (ngpr=257, ngpz=257, ngpf=257)
      real*8 gpr(ngpr),gpz(ngpz)
      real*8 pfmd(ngpr,ngpz)
      real*8 rcntc,btorc,psilim
      character*256 arg,filename
      character expnam*4, dianam*3
      integer nshot,nedit,iarg,iargc,lnblnk
      real tshot
      integer mdim,lpfxdim,lin
      parameter(mdim=256,lpfxdim=4,lin=1)
      real ri(0:mdim),zj(0:mdim),pfm(0:mdim,0:mdim),
     1 pfxx(0:lpfxdim),rpfx(0:lpfxdim),zpfx(0:lpfxdim),
     2 rin(lin),zin(lin),br(lin),bz(lin),bt(lin),fpf(lin),fjp(lin)
c=====================================================
c
      dianam='FPP'
      expnam='AUGD'
      m=mdim
      n=mdim
      nedit=0
      iarg=iargc()
      write(*,*) 'iarg',iarg
      if(iarg.ge.2) then
	call getarg(1,arg)
	read(arg,*) nshot
	call getarg(2,arg)
	read(arg,*) tshot
	if(iarg.ge.3) call getarg(3,dianam)
	if(iarg.ge.4) call getarg(4,expnam)
	if(iarg.ge.5) then
	  call getarg(5,arg)
	  read(arg,*) nedit
	endif
      else
	write(*,*) '1st arg == shot number'
	write(*,*) '2nd arg == time'
	write(*,*) '3rd arg == diagnostic (FPP)'
	write(*,*) '4th arg == experiment (AUGD)'
	stop
      endif
      write(*,*) expnam,dianam,nshot,tshot
      write(filename,'(i5,''.'',i4,''.'',a4,''.'',a3,''.'',i2,''.eq'')')
     1 nshot,nint(tshot*1000),expnam,dianam,nedit
      do i=1,10
	if(filename(i:i).eq.' ') filename(i:i)='0'
      enddo
      do i=21,22
	if(filename(i:i).eq.' ') filename(i:i)='0'
      enddo
      do i=11,lnblnk(filename)
	if(filename(i:i).eq.' ') filename(i:i)='_'
      enddo
      
      call kkEQPFM(ierr,expnam,dianam,nshot,nedit,tshot,mdim+1,
     1 m,n,ri,zj,pfm)
      write(*,*) 'kkEQPFM: ierr',ierr
      write(*,*) 'm,n=',m,n
      lpfx=lpfxdim
      call kkEQpfx(iERR,expnam,dianam,nSHOT,nEDIT,tSHOT,
     1 LPFx,PFxx,RPFx,zPFx)
      write(*,*) 'kkEQpfx: ierr',ierr
      rin(1)=rpfx(0)
      zin(1)=zpfx(0)
      write(*,*) 'rpfx',rpfx
      write(*,*) 'zpfx',zpfx
      call kkrzBrzt (iERR,expnam,dianam,nSHOT,nEDIT,tSHOT,
     1 rin,zin,lin,Br,Bz,Bt,fPF,fJp)
      write(*,*) 'kkrzBrzt: ierr',ierr

      rcntc=rin(1)
      btorc=bt(1)
      psilim=pfxx(1)/(2.0*3.14159 26535 89793 23846)
      do i=0,m
	do j=0,n
	  pfmd(i+1,j+1)=pfm(i,j)/(2.0*3.14159 26535 89793 23846)-psilim
	enddo
      enddo
      do i=0,m
	gpr(i+1)=ri(i)
      enddo
      ngr=m+1
      do j=0,n
	gpz(j+1)=zj(j)
      enddo
      ngz=n+1

      print *,'psilim = ',psilim
      print *,'rcntc = ',rcntc
      print *,'btorc = ',btorc
      open(2,file=filename)
      psilim=0
      write(*,*) 'ngr,ngz=',ngr,ngz
      call wreqdg(2,ngpr,iret,ngr,ngz,psilim,
     1 btorc,rcntc,gpr,gpz,pfmd)
      if(iret.ne.0) then
          print *,'==== shotfile2dg: error in wreqdg. iret = ',iret
      end if
c
      end
c
      integer function lnblnk(string)
c
c returns the position of the last non-blank character in "string"
c
      character*(*) string
      integer i
      do i=len(string),1,-1
        if(string(i:i).ne.' ') then
          lnblnk=i
          return
        endif
      enddo
      lnblnk=0
      return
      end
