subroutine read_geometry()
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
     dim2d(1)=zones(k)%mesh%Nx+2
     dim2d(2)=zones(k)%mesh%Nz+2
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     allocate(buffer(0:Nx+1,0:Nz+1))
     call h5dopen_f(group_id,"Rgeom",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer,dim2d,error)
     zones(k)%mesh%Rgeom(0:Nx+1,0:Nz+1)=buffer
     call h5dclose_f(dataset_id,error)     
     call h5dopen_f(group_id,"Zgeom",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer,dim2d,error)
     zones(k)%mesh%Zgeom(0:Nx+1,0:Nz+1)=buffer
     call h5dclose_f(dataset_id,error)     
     call h5dopen_f(group_id,"Br",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer,dim2d,error)
     zones(k)%mesh%Br(0:Nx+1,0:Nz+1)=buffer
     call h5dclose_f(dataset_id,error)     
     call h5dopen_f(group_id,"Bz",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer,dim2d,error)
     zones(k)%mesh%Bz(0:Nx+1,0:Nz+1)=buffer
     call h5dclose_f(dataset_id,error)     
     call h5dopen_f(group_id,"Bphi",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer,dim2d,error)
     zones(k)%mesh%Bphi(0:Nx+1,0:Nz+1)=buffer
     call h5dclose_f(dataset_id,error)     
     deallocate(buffer)
     allocate(buffer(1:Nx+1,1:Nz+1))
     dim2d(1)=zones(k)%mesh%Nx+1
     dim2d(2)=zones(k)%mesh%Nz+1
     call h5dopen_f(group_id,"Rcorner",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer,dim2d,error)
     zones(k)%mesh%Rcorner(1:Nx+1,1:Nz+1)=buffer
     call h5dclose_f(dataset_id,error)     
     call h5dopen_f(group_id,"Zcorner",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer,dim2d,error)
     zones(k)%mesh%Zcorner(1:Nx+1,1:Nz+1)=buffer
     call h5dclose_f(dataset_id,error)     
     call h5gclose_f(group_id,error) !close zone
     deallocate(buffer)
     zones(k)%mesh%B=sqrt(Zones(k)%mesh%Br**2.D0+Zones(k)%mesh%Bz**2.D0+Zones(k)%mesh%Bphi**2.D0)

  end do
  call h5fclose_f(file_id,error) !close mesh.h5
  call h5close_f(error) ! close hdf5
end subroutine read_geometry

