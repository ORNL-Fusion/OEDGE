!------------------------------------------------------------
! File: save.f90
! Project: OSMtsolver
! Date: 02/04/2009
! Subroutines for saving data
!------------------------------------------------------------


!----------------------------------------
! Saves the data array to the wanted file
!----------------------------------------
subroutine saveDATA(Xfields,Nz)
   use prec_const
   implicit none

   real(float), dimension(Nz,4) :: Xfields
   integer :: iz, Nz

   open(101,file='Result/N.dat',status='unknown', access='append', form='unformatted')
   open(102,file='Result/V.dat',status='unknown', access='append', form='unformatted')
   open(103,file='Result/Te.dat',status='unknown', access='append', form='unformatted')
   open(104,file='Result/Ti.dat',status='unknown', access='append', form='unformatted')

   write(101) (Xfields(iz,1), iz=1,Nz)
   write(102) (Xfields(iz,2), iz=1,Nz)
   write(103) (Xfields(iz,3), iz=1,Nz)
   write(104) (Xfields(iz,4), iz=1,Nz)

   close(101)
   close(102)
   close(103)
   close(104)

end subroutine saveDATA



!----------------------------------------------------------
! Create (or clear) output files and saves initial
! conditions
!----------------------------------------------------------
subroutine saveDATA_first(Xfields,zmesh,Bamp,Nz)
   use prec_const
   implicit none

   real(float), dimension(Nz,4) :: Xfields
   real(float), dimension(Nz) :: zmesh, Bamp
   integer :: iz, Nz

   open(101,file='Result/N.dat',status='unknown', form='unformatted')
   open(102,file='Result/V.dat',status='unknown', form='unformatted')
   open(103,file='Result/Te.dat',status='unknown', form='unformatted')
   open(104,file='Result/Ti.dat',status='unknown', form='unformatted')
   open(105,file='Result/zmesh.dat',status='unknown', form='unformatted')
   open(106,file='Result/Bamp.dat',status='unknown', form='unformatted')

   write(101) (Xfields(iz,1), iz=1,Nz)
   write(102) (Xfields(iz,2), iz=1,Nz)
   write(103) (Xfields(iz,3), iz=1,Nz)
   write(104) (Xfields(iz,4), iz=1,Nz)
   write(105) (zmesh(iz), iz=1,Nz)
   write(106) (Bamp(iz), iz=1,Nz)

   close(101)
   close(102)
   close(103)
   close(104)
   close(105)
   close(106)

end subroutine saveDATA_first
