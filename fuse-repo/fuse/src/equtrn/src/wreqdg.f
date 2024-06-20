      subroutine wreqdg(lun,ngpr,iret, nr,nz, psib,
     ,                                             btf,rtf,rgr,zgr,pfm)
c=====================================================
c*** Write the equilibrium data in the dg-compatible format.
c***
c*** Input:
c***  lun     the logical unit number for the output
c***  ngpr    the maximum number of points in R direction
c***  nr      the actual number of points in R direction
c***  nz      the actual number of points in Z direction
c***  psib    the poloidal flux at the separatrix
c***  btf     the toroidal magnetic field at the R=rtf
c***  rgr     the R values for the grid points
c***  zgr     the Z values for the grid points
c***  pfm     the values of the poloidal flux
c***
c*** Output:
c***  iret    return code (0 means OK)
c=====================================================
c
c  version : 23.06.97 17:23
c
      real*8 rgr(*), zgr(*), pfm(ngpr, *)
c... toroidal field in tesla, radius in m
      real*8 btf, rtf, psib
c
c***  Write the toroidal field and the corresponding radius...
c*** -- now obsolette!
c
c      write(3,*,err=99) 'Toroidal field in Tesla, radius in m'
c      write(3,*,err=99) btf, rtf
c
c***  ... then the plasma equilibrium ...
c
      iret=0
      write(lun,*,err=99)
     /    '   jm   :=  no. of grid points in radial direction;'
      write(lun,*,err=99)
     /    '   km   :=  no. of grid points in vertical direction;'
      write(lun,*,err=99)
     /    '   r    :=  radial   co-ordinates of grid points  [m];'
      write(lun,*,err=99)
     /    '   z    :=  vertical co-ordinates of grid points  [m];'
      write(lun,*,err=99)
     /    '   psi  :=  flux per radiant at grid points     [Wb/rad];'
      write(lun,*,err=99)
     /    '   psib :=  psi at plasma boundary              [Wb/rad];'
      write(lun,*,err=99)
     /    '   btf  :=  toroidal magnetic field                  [t];'
      write(lun,*,err=99)
     /    '   rtf  :=  major radius at which btf is specified   [m];'
      write(lun,*,err=99)
      write(lun,*,err=99)
      write(lun,*,err=99) '   jm    = ', nr,';'
      write(lun,*,err=99) '   km    = ', nz,';'
      write(lun,*,err=99) '   psib  = ',psib,' Wb/rad;'
      write(lun,*,err=99) '   btf   = ',btf,' t;'
      write(lun,*,err=99) '   rtf   = ',rtf,' m;'
      write(lun,*,err=99)
      write(lun,*,err=99) '   r(1:jm);'
      write(lun,8000) (rgr(i),i=1,nr)
      write(lun,*)
      write(lun,*,err=99) '   z(1:km);'
      write(lun,8000) (zgr(i),i=1,nz)
      write(lun,*)
      write(lun,*) '     ((psi(j,k)-psib,j=1,jm),k=1,km)'
      write(lun,8000) ((pfm(i,j)-psib,i=1,nr),j=1,nz)
 8000 format(5(3x,e14.8))
      iret=0
      return
c-----------------------------------------------------
c
 99   print *,'==== wreqdg: error in the input files'
      iret=8
c
      end
