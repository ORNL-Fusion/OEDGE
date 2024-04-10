      subroutine rdeqfu(lun,ngpr,ngpz,iret, nr,nz,
     ,           rmin,zmin,rmax,zmax,delr,delz,psilim,
     ,           btorc,rcntc,fg,pfm,rgr,zgr)
c=======================================================================
c*** where:
c***
c*** i)  nr, nz, rmin, rmax, delr, and delz define the rectangular
c***     mesh used to store the psi values via the functions rr & zz
c***     defined at the top of these code,
c***
c*** ii) psilim is the value at the separatrix,
c***
c*** iii) btorc is the toroidal magnetic field at a radius rcntc,
c***
c*** iv) the poloidal flux is pfm,
c***
c*** v) fg is the flux function R*Btor, rmax and zmax are the grid edges
c=======================================================================
cak-20051115: removed flux reduction to Wb/rad - it is already there
c
c  version : 15.11.2005 14:18
c
      implicit none
      integer lun,ngpr,ngpz,iret, nr,nz
      real*8 pfm(ngpr,ngpz),rgr(ngpr),zgr(ngpz)
      real*8 rmin,zmin,rmax,zmax,delr,delz,psilim,btorc,rcntc,fg
      integer i,j,k,l
c      real*8 u,r,z,p,rr,zz,pi
      real*8 u,r,z,p,rr,zz
      character*19 c,cc
c=======================================================================
      rr(i)=delr*(i-1)+rmin
      zz(i)=delz*(i-1)+zmin
c=======================================================================
c
      iret=0
      rewind lun
      read(lun,*)
      read(lun,*)
      read(lun,*) nr,nz,delr,delz,psilim
      read(lun,*)
      if(nr.gt.ngpr) then
          print *,'=== rdeqfu: nr > ngpr'
          iret=2
      end if
      if(nz.gt.ngpz) then
          print *,'=== rdeqfu: nz > ngpz'
          iret=2
      end if
      if(nr.le.0) then
          print *,'=== rdeqfu: nr < 1'
          iret=4
      end if
      if(nz.le.0) then
          print *,'=== rdeqfu: nz < 1'
          iret=4
      end if
      if(iret.ne.0) return

      read(lun,*)
      do i=1,nr !{
        rgr(i)=0
      end do !}
      do j=1,nz !{
        zgr(j)=0
      end do !}
c      pi=4.*atan(1.)
      do j=1,nz !{
        do i=1,nr !{
c          read(lun,700) c,r,z,p
          read(lun,*) c,r,z,p
          k=index(c,'(')
          l=index(c,')')
          if(k.eq.0 .or. l.eq.0 .or. l.le.k+2) go to 90
          cc=c(k+1:l-1)
          read(cc,*,err=90) k,l
          if(k.ne.i .or. l.ne.j) then !{
            write(*,*) '=== rdeqfu: wrong R,Z order : ',i,j
            stop
          end if !}
          if(j.eq.1) then !{
            rgr(i)=r
          else !}{
            if(rgr(i).ne.r) then !{
              write(*,*) '=== rdeqfu: different R : ',i,j
              stop
            end if !}
          end if !}
          if(i.eq.1) then !{
            zgr(j)=z
          else !}{
            if(zgr(j).ne.z) then !{
              write(*,*) '=== rdeqfu: different Z : ',i,j
              stop
            end if !}
          end if !}
c          pfm(i,j)=p/(2.*pi)
          pfm(i,j)=p
        end do !}
      end do !}
c      psilim=psilim/(2.*pi)
      rmin=rgr(1)
      zmin=zgr(1)
      rmax=rgr(nr)
      zmax=zgr(nz)
      delr=rgr(2)-rgr(1)
      delz=zgr(2)-zgr(1)
      return

 90   write(*,*) '=== rdeqfu: wrong line identification format',
     ,                                        ' - must be (i, j) : ',i,j
      stop
  700 format(a19,3e20.0)
c=======================================================================
      end
