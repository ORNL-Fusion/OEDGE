      subroutine rdeqdg(lun,ngpr,ngpz,iret, nr,nz,
     ,                                        psib,btf,rtf,rgr,zgr,pfm)
c=====================================================
c*** Read the equilibrium data written in the dg-compatible format.
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
c***  psib    the psi value at the separatrix
c***  btf     the toroidal magnetic field at the R=rtf
c***  rgr     the R values for the grid points
c***  zgr     the Z values for the grid points
c***  pfm     the values of the (poloidal flux - separatrix flux)
c***
c*** If the input file contains no data on the toroidal field (old
c*** version), then the field file must be pre-connected to LUN 3
c=====================================================
c
c  version : 05.08.95 16:50
c
      real*8 rgr(*), zgr(*), pfm(ngpr, *)
c... toroidal field in tesla, radius in m
      real*8 btf, rtf, psib
      real ubtf, urtf
c-----------------------------------------------------
c
c*** Read the plasma equilibrium ...
c
      iret=0
      call rdeqjh(lun,nr,nz,upsib,ubtf,urtf,*99)
      btf=ubtf
      rtf=urtf
      psib=upsib
      if(rtf.le.0.) then
c
c*** Read the toroidal field and the corresponding radius
c*** from a separate file (for compatibility with old versions)
c
          read(3,*,err=99)
          read(3,*,err=99) btf, rtf
          close(3)
      end if
      if(nr.gt.ngpr) then
          print *,'==== rdeqdg: nr > ngpr ',nr
          iret=2
      end if
      if(nz.gt.ngpz) then
          print *,'=== rdeqdg: nz > ngpz'
          iret=2
      end if
      if(nr.le.0) then
          print *,'=== rdeqdg: nr < 1'
          iret=4
      end if
      if(nz.le.0) then
          print *,'=== rdeqdg: nz < 1'
          iret=4
      end if
      if(iret.ne.0) return
c
      read(lun,8000) (rgr(i),i=1,nr)
      read(lun,*)
      read(lun,*)
      read(lun,8000) (zgr(i),i=1,nz)
      read(lun,*)
      read(lun,*)
      read(lun,8000) ((pfm(i,j),i=1,nr),j=1,nz)
 8000 format(5(3x,e14.8))
      iret=0
      return
c-----------------------------------------------------
c
 99   print *,'==== rdeqdg: error in the input files'
      iret=8
c
      end
