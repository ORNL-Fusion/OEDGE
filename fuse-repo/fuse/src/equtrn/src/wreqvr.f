      subroutine wreqvr(lun,ngpr,iret,nr,nz,rgr,zgr,pfm)
c
c  version : 18.01.96 15:02
c
c=====================================================
c*** Wrighting the CARRE equilibrium file
c***
c*** Input:
c***  lun     the logical unit number for the output
c***  ngpr    the maximum number of points in R direction
c***  nr      the actual number of points in R direction
c***  nz      the actual number of points in Z direction
c***  rgr     the R values for the grid points
c***  zgr     the Z values for the grid points
c***  pfm     the values of the poloidal flux
c***
c*** Output:
c***  iret    return code (0 means OK)
c=====================================================
      real*8 pfm(ngpr,*),rgr(ngpr),zgr(*)
c=====================================================
c
      write(lun,*,err=999)
      write(lun,'(i5,i6)',err=999) nr,nz
      write(lun,*,err=999) '$r'
      write(lun,710,err=999) 'nr',nr
      write(lun,720,err=999) (rgr(i),i=1,nr)
      write(lun,*,err=999) '$z'
      write(lun,710,err=999) 'nz',nz
      write(lun,720,err=999) (zgr(i),i=1,nz)
      write(lun,*,err=999)
      write(lun,*,err=999) '$psi'
      write(lun,720,err=999) ((pfm(i,j),i=1,nr),j=1,nz)
      iret=0
      return
c=====================================================
 999  iret=1
c=====================================================
  710 format(1x,a2,'=',i3)
  720 format(1p,1x,6e13.5)
c
      end
