subroutine save_vorticity_test(level)
  use test_var
  use all_variables, only : zones, global_parameters
  use hdf5
  implicit none
  integer*4,intent(in) :: level
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
  do n=0,global_parameters%N_ions
     write(filename,"(I0,A17,I0)") level,"_electric_fields_",n
     call h5fcreate_f(trim("Results/"//filename),H5F_ACC_TRUNC_F,file_id,error) 
     do k=1,global_parameters%N_zones
        write(groupname,"(A4,I0)") "zone",k
        call h5gcreate_f(file_id,trim(groupname),group_id,error)
        Nx=zones(k)%mesh%Nx
        Nz=zones(k)%mesh%Nz
        dim2d(1)=Nx+2
        dim2d(2)=Nz+2
        call h5screate_simple_f(2,dim2d,dspace_id,error)
        call h5dcreate_f(group_id,"phi",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%electric_fields(1)%phi,dim2d,error)
        call h5dclose_f(dataset_id,error)
        call h5sclose_f(dspace_id,error)
        call h5screate_simple_f(2,dim2d,dspace_id,error)
        call h5dcreate_f(group_id,"vorticity",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%electric_fields(1)%vorticity,dim2d,error)
        call h5dclose_f(dataset_id,error)
        call h5sclose_f(dspace_id,error)
        call h5screate_simple_f(2,dim2d,dspace_id,error)
        call h5dcreate_f(group_id,"phi2",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%electric_fields(2)%phi,dim2d,error)
        call h5dclose_f(dataset_id,error)
        call h5sclose_f(dspace_id,error)
        call h5screate_simple_f(2,dim2d,dspace_id,error)
        call h5dcreate_f(group_id,"vorticity2",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%electric_fields(2)%vorticity,dim2d,error)
        call h5dclose_f(dataset_id,error)
        call h5sclose_f(dspace_id,error)
        call h5screate_simple_f(2,dim2d,dspace_id,error)
        call h5dcreate_f(group_id,"vorticity3",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,test_vort(k)%W,dim2d,error)
        call h5dclose_f(dataset_id,error)
        call h5sclose_f(dspace_id,error)
        call h5screate_simple_f(2,dim2d,dspace_id,error)
        call h5dcreate_f(group_id,"phi3",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,test_vort(k)%phi,dim2d,error)
        call h5dclose_f(dataset_id,error)
        call h5sclose_f(dspace_id,error)
        call h5screate_simple_f(2,dim2d,dspace_id,error)
        call h5dcreate_f(group_id,"r",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,test_vort(k)%r,dim2d,error)
        call h5dclose_f(dataset_id,error)
        call h5sclose_f(dspace_id,error)
        call h5screate_simple_f(2,dim2d,dspace_id,error)
        call h5dcreate_f(group_id,"theta",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,test_vort(k)%theta,dim2d,error)
        call h5dclose_f(dataset_id,error)
        call h5sclose_f(dspace_id,error)
        call h5gclose_f(group_id,error)
     end do
     call h5fclose_f(file_id,error)
  end do
  call h5close_f(error)
end subroutine save_vorticity_test
