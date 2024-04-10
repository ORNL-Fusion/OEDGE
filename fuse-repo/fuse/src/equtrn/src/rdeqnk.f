      subroutine rdeqnk(lun,ngpr,ngpz,iret, nr,nz,
     ,           rmin,zmin,rmax,zmax,delr,delz,psilim,
     ,           btorc,rcntc,fg,pfm,rgr,zgr)
c=====================================================
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
c=====================================================
c
c  version : 16.01.96 19:33
c
      real*8 pfm(ngpr,*),rgr(*),zgr(*)
      real*8 rmin,zmin,rmax,zmax,delr,delz,psilim,btorc,rcntc,fg,u
c=====================================================
      rr(i)=delr*(i-1)+rmin
      zz(i)=delz*(i-1)+zmin
c=====================================================
c
      iret=0
      rewind lun
      read(lun,*)
      read(lun,*)
      read(lun,*) u,u,u,fg
      read(lun,*)
      read(lun,*) nr,nz
      if(nr.gt.ngpr) then
          print *,'=== rdefit: nr > ngpr'
          iret=2
      end if
      if(nz.gt.ngpz) then
          print *,'=== rdefit: nz > ngpz'
          iret=2
      end if
      if(nr.le.0) then
          print *,'=== rdefit: nr < 1'
          iret=4
      end if
      if(nz.le.0) then
          print *,'=== rdefit: nz < 1'
          iret=4
      end if
      if(iret.ne.0) return
c
      read(lun,*)
      read(lun,*) rmin,rmax,zmin,zmax,delr,delz
      read(lun,*)
      read(lun,*)
      read(lun,*) rcntc,u,psilim
      read(lun,*)
      read(lun,*) ((pfm(i,j),i=1,nr),j=1,nz)
c
      btorc=fg/rcntc
      do i=1,nr
        rgr(i)=rr(i)
      enddo
      do i=1,nz
        zgr(i)=zz(i)
      enddo
c
      end
