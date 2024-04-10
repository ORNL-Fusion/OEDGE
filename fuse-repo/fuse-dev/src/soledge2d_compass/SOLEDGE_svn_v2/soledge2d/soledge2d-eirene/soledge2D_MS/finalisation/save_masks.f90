subroutine save_masks()
  use all_variables, only : zones, global_parameters
  use hdf5
  implicit none
  !HDF5 variable
  integer(HID_T) :: file_id
  integer(HID_T) :: group_id
  integer(HID_T) :: group2_id
  integer(HID_T) :: dataset_id
  integer(HID_T) :: dspace_id
  integer :: error
  character(50) :: filename
  character(50) :: groupname
  integer(HSIZE_T) :: dim1d(1)
  integer(HSIZE_T) :: dim2d(2)
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
  call h5open_f(error)
  filename="Results/masks"
  call h5fcreate_f(trim(filename),H5F_ACC_TRUNC_F,file_id,error) 
  do k=1,global_parameters%N_zones
     write(groupname,"(A4,I0)") "zone",k
     call h5gcreate_f(file_id,trim(groupname),group_id,error)
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     dim2d(1)=Nx+2
     dim2d(2)=Nz+2
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"chi1",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%masks%chi1,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"chi2",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%masks%chi2,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"chi3",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%masks%chi3,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"chi4",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%masks%chi4,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"chi5",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%masks%chi5,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"chi6",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%masks%chi6,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)

     call h5gclose_f(group_id,error)
  end do
  call h5fclose_f(file_id,error)
  call h5close_f(error)
end subroutine save_masks
