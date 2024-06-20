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
      parameter (ngpr=257, ngpz=257)
      real*8 fg(ngpr),pg(ngpr),ffg(ngpr),ppg(ngpr)
      real*8 pfm(ngpr,ngpz),rgr(ngpr),zgr(ngpz)
      real*8 rdim,zdim,rcntc,redge,zmsmid,rma,zma,psimin,psilim,btorc
      character title*40, date*8
      integer chop
c=====================================================
c
      character*256 filename
c      character tab
c----------------------------------------------------------------------

c      tab=char(9)

      if(iargc().ne.2) then
	write(*,*) '1st arg == input dg equilibrium file'
	write(*,*) '2nd arg == output dg equilibrium file'
	stop
      endif

      call getarg(1,filename)
      open(1,file=filename)
 
      call rdefit(1,ngpr,ngpz,iret, title,date,ipestg,nr,nz,
     ,           rdim,zdim,zmsmid,rcntc,redge,rma,zma,psimin,psilim,
     ,           btorc,fg,pg,ffg,ppg,pfm,rgr,zgr)
      if(iret.ne.0) then
          print *,'==== ef2dg: error in rdefit. iret =',iret
          stop
      end if

c need to add something here to remove the first few columns
      chop=2
      rgr(1:nr-chop)=rgr(1+chop:nr)
      pfm(1:nr-chop,1:nz)=pfm(1+chop:nr,1:nz)
      nr=nr-chop
c of positive flux

      call double(pfm,ngpr,ngpz,rgr,nr,zgr,nz)
      call double(pfm,ngpr,ngpz,rgr,nr,zgr,nz)

c      psilim=0

      call getarg(2,filename)
      open(2,file=filename)

      call wreqdg(2,ngpr,iret,nr,nz,psilim,btorc,rcntc,rgr,zgr,pfm)
      if(iret.ne.0) then
	  print *,'==== dg2dg: error in wreqdg. iret = ',iret
      end if

      end


c======================================================================
      subroutine double(pfm,ngpr,ngpz,rgr,nr,zgr,nz)

      implicit none
      integer ir,ifail,nn,nz2,iz,nr2,nr,ngpz,ngpr,j,i,nz
      real*8 pfm(ngpr,ngpz),rgr(ngpr),zgr(ngpz)
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

      if(ngpr.lt.2*nr-1) stop 'increase ngpr'
      if(ngpz.lt.2*nz-1) stop 'increase ngpz'

      nn=nr*nz
      call splbcb(rgrs, nr, zgrs, nz, pfms, nr, cx, nn*3, cy, nn*3,
     ,					      cc, nn*9, wrkz, nz, ifail)
      if(ifail.ne.0) then
	write(*,*) 'splbcb: ifail ',ifail
	stop
      endif

      do ir=1,nr
	rgr(2*ir-1)=rgrs(ir)
      enddo
      nr2=2*nr-1
      do ir=2,nr2,2
	rgr(ir)=0.5*(rgr(ir-1)+rgr(ir+1))
      enddo
      do iz=1,nz
	zgr(2*iz-1)=zgrs(iz)
      enddo
      nz2=2*nz-1
      do iz=2,nz2,2
	zgr(iz)=0.5*(zgr(iz-1)+zgr(iz+1))
      enddo

      rgrss(1:nr2)=rgr(1:nr2)
      zgrss(1:nz2)=zgr(1:nz2)

      call evlbcs(rgrss, nr2, zgrss, nz2, pfmss, ngpr,
     ,            rgrs, nr, zgrs, nz, pfms, nr, cx, nn*3, cy, nn*3,
     ,                            cc, nn*9, wrkr, 4*nr, wrkz, nz, ifail)
      if(ifail.ne.0) then
	write(*,*) 'evlbcs: ifail ',ifail
	if(ifail.lt.0) stop
      end if

      nr=nr2
      nz=nz2
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
