subroutine read_chis()
  use all_variables, only : global_parameters,zones
  use Mzone
  use hdf5
  implicit none
    !HDF5 variable
  integer(HID_T) :: file_id
  integer(HID_T) :: group_id
  integer(HID_T) :: dataset_id
  integer :: error
  character(7),parameter :: filename="mesh.h5"
  character(50) :: groupname
  integer(HSIZE_T) :: dim2d(2)
  !other variable
  integer*4 :: k
  integer*4 :: Nx,Nz
  real*8,allocatable :: buffer(:,:)
  call h5open_f(error)
  call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,error)
  do k=1,global_parameters%N_zones
     write(groupname,"(A4,I0)") "zone",k
     call h5gopen_f(file_id,trim(groupname),group_id,error)
     dim2d(1)=zones(k)%mesh%Nx
     dim2d(2)=zones(k)%mesh%Nz
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     allocate(buffer(1:Nx,1:Nz))
     call h5dopen_f(group_id,"chi",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer,dim2d,error)
     zones(k)%masks%chi(1:Nx,1:Nz)=buffer
     call h5dclose_f(dataset_id,error)     
     call h5gclose_f(group_id,error) !close zone
     deallocate(buffer)
  end do
  call h5fclose_f(file_id,error) !close mesh.h5
  call h5close_f(error) ! close hdf5
end subroutine read_chis

