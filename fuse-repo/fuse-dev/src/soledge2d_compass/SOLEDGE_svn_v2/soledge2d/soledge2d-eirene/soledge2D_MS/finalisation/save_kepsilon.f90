subroutine save_kepsilon()
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
  !other variable
  integer*4 :: Nx,Nz
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
  filename='Results/k_epsilon'
  call h5fcreate_f(trim(filename),H5F_ACC_TRUNC_F,file_id,error) 
  !Save zones
  do k=1,global_parameters%N_zones
     write(groupname,"(A4,I0)") "zone",k
     call h5gcreate_f(file_id,trim(groupname),group_id,error)
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     dim2d(1)=Nx+2
     dim2d(2)=Nz+2
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"k",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%kepsilon(1)%k,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"epsilon",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%kepsilon(1)%epsilon,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     dim2d(1)=Nx
     dim2d(2)=Nz
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"interchange",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%kepsilon(1)%interchange,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"kh",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%kepsilon(1)%kh,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"UEshear",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%kepsilon(1)%UE_shear,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)

     call h5gclose_f(group_id,error)
  end do

  if(flags%use_triangles) then
     groupname='triangles'
     call h5gcreate_f(file_id,trim(groupname),group_id,error)
     dim1d=interp_data2%N_knots
     call h5screate_simple_f(1,dim1d,dspace_id,error)
     call h5dcreate_f(group_id,"k",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,interp_data2%knots_k(:),dim1d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)        
     call h5screate_simple_f(1,dim1d,dspace_id,error)
     call h5dcreate_f(group_id,"epsilon",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,interp_data2%knots_epsilon(:),dim1d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)        
     call h5gclose_f(group_id,error)
  end if

  call h5fclose_f(file_id,error)
  call h5close_f(error)
end subroutine save_kepsilon
