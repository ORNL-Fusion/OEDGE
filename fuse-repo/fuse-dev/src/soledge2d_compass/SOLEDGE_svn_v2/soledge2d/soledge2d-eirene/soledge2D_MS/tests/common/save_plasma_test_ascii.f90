subroutine save_plasma_test_ascii(level)
  use test_var
  use all_variables, only : zones, global_parameters, reference_parameters
  use hdf5
  implicit none
  integer*4,intent(in) :: level
  !other variable
  integer*4 :: Nx,Nz
  integer*4 :: k,n
  logical :: dir_e
#include "compile_opt.inc"
#if GFORTRAN==1
  inquire(File='Results',exist=dir_e)
#endif
#if GFORTRAN==0
  inquire(Directory='Results',exist=dir_e)
#endif
  if(.not.dir_e) then
     call system("mkdir Results")
  end if
  open(unit=80,file='Results/reference_parameters.txt',status='unknown')
  write(80,100) reference_parameters%fields%T0eV
  write(80,100) reference_parameters%fields%n0
  close(80)
  open(unit=80,file='Results/all_plasmas',status='unknown')
  do n=0,global_parameters%N_ions
     write(80,100) zones(1)%species(n)%var(1)%density(1,1)
100 format(512es15.7)
  end do
  close(80)
  open(unit=80,file='Results/rad_plasmas',status='unknown')
  do n=0,global_parameters%N_ions
     write(80,100) zones(1)%species(n)%sources%rad(1,1)
  end do
  close(80)
end subroutine save_plasma_test_ascii
