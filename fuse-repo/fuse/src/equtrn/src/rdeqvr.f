      subroutine rdeqvr(lun,ngpr,ngpz,iret, nr,nz, btf,rtf,rgr,zgr,pfm)
c======================================================================
c*** Read the equilibrium data from TdeV
c***
c*** Input:
c***  lun     the logical unit number for the input
c***  ngpr    the maximum number of points in R direction
c***  ngpz    the maximum number of points in Z direction
c***
c*** Output:
c***  iret    return code (0 means OK)
c***  nr      the actual number of points in R direction
c***  nz      the actual number of points in Z direction
c***  btf     the toroidal magnetic field at the R=rtf
c***  rgr     the R values for the grid points
c***  zgr     the Z values for the grid points
c***  pfm     the values of the (poloidal flux - separatrix flux)
c***
c*** If the input file contains no data on the toroidal field (old
c*** version), then the field file must be pre-connected to LUN 3
c======================================================================
c
c  version : 06.02.99 20:52
c
      real*8 rgr(*), zgr(*), pfm(ngpr, *)
c... toroidal field in tesla, radius in m
      real*8 btf, rtf
      real ubtf, urtf
      integer i
      character ul*72,uch*4
c-----------------------------------------------------
c
c*** Read the plasma equilibrium ...
c
      rtf=0.
c
      uch='$r'
 10   read(lun,'(a72)',end=98) ul
      if(index(ul,'$r').le.0) go to 10
c
      uch='nr'
      read(lun,'(a72)',end=98) ul
      i=index(ul,'=')
      if(i.le.0) go to 99
      read(ul(i+1:72),*,err=99) nr
      if(nr.gt.ngpr) then
          print *,'==== rdeqvr: nr > ngpr ',nr
          iret=2
          return
      else if(nr.le.0) then
          print *,'=== rdeqvr: nr < 1'
          iret=4
          return
      end if
      uch='r'
      read(lun,*,err=99) (rgr(i),i=1,nr)
c
      uch='$z'
 15   read(lun,'(a72)',end=98) ul
      if(index(ul,'$z').le.0) go to 15
c
      uch='nz'
      read(lun,'(a72)',end=98) ul
      i=index(ul,'=')
      if(i.le.0) go to 99
      read(ul(i+1:72),*,err=99) nz
      if(nz.gt.ngpz) then
          print *,'=== rdeqvr: nz > ngpz'
          iret=2
          return
      else if(nz.le.0) then
          print *,'=== rdeqvr: nz < 1'
          iret=4
          return
      end if
      uch='z'
      read(lun,*,err=99) (zgr(i),i=1,nz)
c
      uch='$psi'
 30   read(lun,'(a72)',end=98) ul
      if(index(ul,'$psi').le.0) go to 30
c
      uch='psi'
      read(lun,*,err=99) ((pfm(i,j),i=1,nr),j=1,nz)
c
      if(rtf.le.0.) then
c
c*** Read the toroidal field and the corresponding radius
c*** from a separate file (for compatibility with old versions)
c
          uch='btf'
 40       continue
          read(3,*,err=40, end=99) btf
	  rtf=1.
          close(3)
      end if
      iret=0
      return
c-----------------------------------------------------
c
 98   print *,'==== rdeqvr: no label ',uch
      iret=8
      return
c
 99   print *,'==== rdeqvr: error by reading ',uch
      iret=8
c
      end
