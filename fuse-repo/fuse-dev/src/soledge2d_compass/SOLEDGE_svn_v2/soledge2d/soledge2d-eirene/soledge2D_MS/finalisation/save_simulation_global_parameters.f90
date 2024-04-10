subroutine save_simulation_global_parameters()
  use all_variables, only : global_variables
  use hdf5
  implicit none
  !HDF5 variable
  integer(HID_T) :: file_id
  integer(HID_T) :: dataset_id
  integer(HID_T) :: dspace_id
  integer :: error
  character(50) :: filename
  integer(HSIZE_T) :: dim1d(1)
  !other variable
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
  call h5open_f(error)
  filename='Results/globals'
  call h5fcreate_f(trim(filename),H5F_ACC_TRUNC_F,file_id,error) 
  dim1d(1)=1
  call h5screate_simple_f(1,dim1d,dspace_id,error)
  call h5dcreate_f(file_id,"tempus",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
  call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,global_variables%tempus,dim1d,error)
  call h5dclose_f(dataset_id,error)
  call h5sclose_f(dspace_id,error)
  call h5fclose_f(file_id,error)
  call h5close_f(error)
end subroutine save_simulation_global_parameters
