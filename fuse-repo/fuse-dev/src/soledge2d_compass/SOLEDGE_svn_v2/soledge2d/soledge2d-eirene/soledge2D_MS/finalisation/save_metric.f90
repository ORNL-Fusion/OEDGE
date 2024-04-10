subroutine save_metric()
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
  filename="Results/metric"
  call h5fcreate_f(trim(filename),H5F_ACC_TRUNC_F,file_id,error) 
  do k=1,global_parameters%N_zones
     write(groupname,"(A4,I0)") "zone",k
     call h5gcreate_f(file_id,trim(groupname),group_id,error)
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     dim2d(1)=Nx+2
     dim2d(2)=Nz+2
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"cpp",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%metric_coefficients%cpp,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"cpt",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%metric_coefficients%cpt,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"ctt",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%metric_coefficients%ctt,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"c_tt",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%metric_coefficients%c_tt,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"c_pp",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%metric_coefficients%c_pp,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"c_pt",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%metric_coefficients%c_pt,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"Jac",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%metric_coefficients%jacobian,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"G",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%metric_coefficients%G,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     !jacobian
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"dPdR",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%jacobian%dPdR,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"dPdZ",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%jacobian%dPdZ,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"dTdR",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%jacobian%dTdR,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"dTdZ",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%jacobian%dTdZ,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"dRdP",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%jacobian%dRdP,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"dZdP",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%jacobian%dZdP,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"dRdT",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%jacobian%dRdT,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"dZdT",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%jacobian%dZdT,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     !surface and volume
     dim2d(1)=Nx
     dim2d(2)=Nz
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"dS_north",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%metric_coefficients%ds_north_dd,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"dS_south",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%metric_coefficients%ds_south_dd,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"dS_east",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%metric_coefficients%ds_east_dd,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"dS_west",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%metric_coefficients%ds_west_dd,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(2,dim2d,dspace_id,error)
     call h5dcreate_f(group_id,"dvol",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%metric_coefficients%dvol_dd,dim2d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5gclose_f(group_id,error)
  end do
  call h5fclose_f(file_id,error)
  call h5close_f(error)
end subroutine save_metric
