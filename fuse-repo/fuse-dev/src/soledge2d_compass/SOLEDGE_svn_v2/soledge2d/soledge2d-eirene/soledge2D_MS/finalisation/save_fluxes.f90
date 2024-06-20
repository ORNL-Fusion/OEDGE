subroutine save_fluxes()
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
  do n=0,global_parameters%N_ions
     write(filename,"(A15,I0)") "Results/fluxes_",n
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
        call h5dcreate_f(group_id,"fluxn",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%species(n)%fluxes%fluxn,dim3d,error)
        call h5dclose_f(dataset_id,error)
        call h5sclose_f(dspace_id,error)
        call h5screate_simple_f(3,dim3d,dspace_id,error)
        call h5dcreate_f(group_id,"fluxG",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%species(n)%fluxes%fluxG,dim3d,error)
        call h5dclose_f(dataset_id,error)
        call h5sclose_f(dspace_id,error)
        call h5screate_simple_f(3,dim3d,dspace_id,error)
        call h5dcreate_f(group_id,"fluxE",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%species(n)%fluxes%fluxE,dim3d,error)
        call h5dclose_f(dataset_id,error)
        call h5sclose_f(dspace_id,error)
        call h5gclose_f(group_id,error)
     end do
     call h5fclose_f(file_id,error)
  end do
  call h5close_f(error)
end subroutine save_fluxes
