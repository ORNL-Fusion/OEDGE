      program cropequ
c
c  version : 29.11.2005 16:17
c
      implicit none
      integer ngpr,ngpz
      parameter (ngpr=1025, ngpz=1025)
      real*8 pfm(ngpr,ngpz),rgr(ngpr),zgr(ngpz)
      real*8 rcntc,psimin,psilim,btorc,rmin,rmax,zmin,zmax
      integer i,j,nr,nz,nrn,nzn,iret,iargc
      logical err
      character*256 filename
c      character tab    !###
c----------------------------------------------------------------------

c      tab=char(9)      !###

      if(iargc().ne.2) then
	write(*,*) '1st arg == input dg equilibrium file'
	write(*,*) '2nd arg == output dg equilibrium file'
        write(*,*) 'Input   == nr, Rmin, Rmax, nz, Zmin, Zmax'
	stop
      endif

      call getarg(1,filename)
      open(1,file=filename)
      call rdeqdg(1,ngpr,ngpz,iret, nr,nz,btorc,rcntc,rgr,zgr,pfm)
      if(iret.ne.0) then
	  print *,'==== dg2dg: error in rdeqdg. iret =',iret
	  stop
      end if
      nrn=nr
      nzn=nz
      rmin=rgr(1)
      rmax=rgr(nr)
      zmin=zgr(1)
      zmax=zgr(nz)
      write (*,*) 
      write (*,*) 'Source grid:'
      write (*,*) '   nr,rmin,rmax=', nr,rmin,rmax
      write (*,*) '   nz,zmin,zmax=', nz,zmin,zmax
      write (*,*) 'Input new values (list-directed format) =>'
      read (*,*,end=20) nrn,rmin,rmax,nzn,zmin,zmax
 20   err=.false.
      if(nrn.lt.3 .or. nrn.gt.ngpr) then !{
        write (0,*) 'Wrong input value of nr: ',nrn,'. Must be between',
     ,    ' 3 and ',ngpr
        err=.true.
      end if !}
      if(nzn.lt.3 .or. nzn.gt.ngpz) then !{
        write (0,*) 'Wrong input value of nz: ',nzn,'. Must be between',
     ,    ' 3 and ',ngpz
        err=.true.
      end if !}
      if(err) then !{
        write(*,*) 'Input   == nr, Rmin, Rmax, nz, Zmin, Zmax'
	stop '==> check ngpr and ngpz in cropequ.f'
      end if !}
      rmin=max(rmin,rgr(1))
      rmax=min(rmax,rgr(nr))
      zmin=max(zmin,zgr(1))
      zmax=min(zmax,zgr(nz))
      if(rmax.le.rgr(1)) rmax=rgr(nr)
      if(rmin.ge.rgr(nr)) rmin=rgr(1)
      if(zmax.le.zgr(1)) zmax=zgr(nz)
      if(zmin.ge.zgr(nz)) zmin=zgr(1)

      call crop(pfm,ngpr,ngpz,rgr,nr,zgr,nz,nrn,rmin,rmax,nzn,zmin,zmax)

      psilim=0
      call getarg(2,filename)
      open(2,file=filename)
      call wreqdg(2,ngpr,iret,nr,nz,psilim,btorc,rcntc,rgr,zgr,pfm)
      if(iret.ne.0) then
	  print *,'==== dg2dg: error in wreqdg. iret = ',iret
      end if

      end

c======================================================================
      subroutine crop(pfm,ngpr,ngpz,rgr,nr,zgr,nz,
     ,                                      nrn,rmin,rmax,nzn,zmin,zmax)

      implicit none
      integer ir,ifail,nn,nz2,iz,nr2,nr,ngpz,ngpr,j,i,nz,nrn,nzn
      real*8 pfm(ngpr,ngpz),rgr(ngpr),zgr(ngpz),rmin,rmax,zmin,zmax,h,x
      real   pfms(nr,nz),rgrs(ngpr),zgrs(ngpz)
      real   pfmss(ngpr,ngpz),rgrss(ngpr),zgrss(ngpz)
      real   cx(nr,3,nz),cy(nz,3,nr),cc(nz,3,nr,3),wrkz(nz),wrkr(nr,4)

c----------------------------------------------------------------------
c      character tab
c
c      tab=char(9)
c----------------------------------------------------------------------

      pfms(1:nr,1:nz)=pfm(1:nr,1:nz)
      rgrs(1:nr)=rgr(1:nr)
      zgrs(1:nz)=zgr(1:nz)

c----------------------------------------------------------------------
c      open(91,file='rgrs.dat')
c      write(91,'(1h#)')
c      write(91,'(1p,i4,a1,e15.6)') (i,tab,rgrs(i),i=1,nr)
c      open(91,file='zgrs.dat')
c      write(91,'(1h#)')
c      write(91,'(1p,i4,a1,e15.6)') (i,tab,zgrs(i),i=1,nz)
c      open(91,file='zpfs.dat')
c      write(91,'(2h# ,2a1,257(a1,i3.3,a1))') 'z',tab,('r',i,tab,i=1,nr)
c      do j=1,nz
c	 write(91,'(1p,258(e11.4,a1))') zgrs(j),tab,
c     , 					  (pfms(i,j),tab,i=1,nr)
c      end do
c      open(91,file='rpfs.dat')
c      write(91,'(2h# ,2a1,257(a1,i3.3,a1))') 'r',tab,('z',i,tab,i=1,nz)
c      do j=1,nr
c	 write(91,'(1p,258(e11.4,a1))') rgrs(j),tab,
c     , 					  (pfms(j,i),tab,i=1,nz)
c      end do
c----------------------------------------------------------------------

c      if(ngpr.lt.2*nr-1) stop 'increase ngpr'
c      if(ngpz.lt.2*nz-1) stop 'increase ngpz'

      nn=nr*nz
      call splbcb(rgrs, nr, zgrs, nz, pfms, nr, cx, nn*3, cy, nn*3,
     ,					      cc, nn*9, wrkr, nz, ifail)
      if(ifail.ne.0) then
	write(*,*) 'splbcb: ifail ',ifail
	stop
      endif

      x=rmin
      h=(rmax-rmin)/(nrn-1)
      do i=1,nrn
	rgr(i)=x
        x=x+h
      end do

      x=zmin
      h=(zmax-zmin)/(nzn-1)
      do i=1,nzn
	zgr(i)=x
        x=x+h
      end do

      rgrss(1:nrn)=rgr(1:nrn)
      zgrss(1:nzn)=zgr(1:nzn)

      call evlbcs(rgrss, nrn, zgrss, nzn, pfmss, ngpr,
     ,            rgrs, nr, zgrs, nz, pfms, nr, cx, nn*3, cy, nn*3,
     ,                            cc, nn*9, wrkr, 4*nr, wrkz, nz, ifail)
      if(ifail.ne.0) then
	write(*,*) 'evlbcs: ifail ',ifail
	if(ifail.lt.0) stop
      end if

      nr=nrn
      nz=nzn
      do i=1,nr !{
        do j=1,nz !{
          pfm(i,j)=pfmss(i,j)
        end do !}
      end do !}
c      pfm(1:nr,1:nz)=pfmss(1:nr,1:nz)

c----------------------------------------------------------------------
c      open(91,file='zpss.dat')
c      write(91,'(2h# ,2a1,257(a1,i3.3,a1))') 'z',tab,('r',i,tab,i=1,nr)
c      do j=1,nz
c	 write(91,'(1p,258(e11.4,a1))') zgr(j),tab,
c     , 					 (pfmss(i,j),tab,i=1,nr)
c      end do
c      open(91,file='rpss.dat')
c      write(91,'(2h# ,2a1,257(a1,i3.3,a1))') 'r',tab,('z',i,tab,i=1,nz)
c      do j=1,nr
c	 write(91,'(1p,258(e11.4,a1))') rgr(j),tab,
c     , 					 (pfmss(j,i),tab,i=1,nz)
c      end do
c----------------------------------------------------------------------

      end
