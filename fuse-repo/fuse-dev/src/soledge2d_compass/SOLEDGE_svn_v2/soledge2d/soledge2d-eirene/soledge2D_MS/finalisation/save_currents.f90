subroutine save_currents()
  use all_variables, only : zones, global_parameters, interp_data2, flags
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
  integer(HSIZE_T) :: dim3d(3)
  !other variable
  integer*4 :: Nx,Nz
  real*8,allocatable :: buffer(:,:)
  integer*4 :: k,n,ind_spec
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
  filename='Results/currents'
  call h5fcreate_f(trim(filename),H5F_ACC_TRUNC_F,file_id,error) 
   !Save zones
  do k=1,global_parameters%N_zones
     write(groupname,"(A4,I0)") "zone",k
     call h5gcreate_f(file_id,trim(groupname),group_id,error)
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     dim3d(1)=Nx
     dim3d(2)=Nz
     dim3d(3)=4
     call h5screate_simple_f(3,dim3d,dspace_id,error)
     call h5dcreate_f(group_id,"j_parallel",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%electric_fields(1)%j_parallel,dim3d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(3,dim3d,dspace_id,error)
     call h5dcreate_f(group_id,"j_adv_para_W",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%electric_fields(1)%j_para_adv_W,dim3d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(3,dim3d,dspace_id,error)
     call h5dcreate_f(group_id,"j_diff_W",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%electric_fields(1)%j_diff_W,dim3d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(3,dim3d,dspace_id,error)
     call h5dcreate_f(group_id,"j_pola",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%electric_fields(1)%j_perp,dim3d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(3,dim3d,dspace_id,error)
     call h5dcreate_f(group_id,"j_diam",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%species(1)%drifts%jdiam,dim3d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(3,dim3d,dspace_id,error)
     call h5dcreate_f(group_id,"j_adv_perp_W",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%electric_fields(1)%j_perp_adv_W,dim3d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5gclose_f(group_id,error)
  end do
  call h5fclose_f(file_id,error)

  call h5close_f(error)

end subroutine save_currents
