subroutine load_ballooning_from_file()
  use all_variables, only : zones, global_parameters
  use hdf5
  implicit none
  !HDF5 variable
  integer(HID_T) :: file_id
  integer(HID_T) :: group_id
  integer(HID_T) :: group2_id
  integer(HID_T) :: dataset_id
  integer :: error
  character(7),parameter :: filename="mesh.h5"
  character(50) :: groupname
  integer(HSIZE_T) :: dim1d(1)
  integer(HSIZE_T) :: dim2d(2)
  !other variable
  integer*4 :: k,n
  integer*4 :: Nx,Nz
  real*8,allocatable :: buffer(:,:),ballD(:,:),ballNu(:,:),ballChi(:,:),ballChie(:,:)
  call h5open_f(error)
  call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,error)
  do k=1,global_parameters%N_zones
     write(groupname,"(A4,I0)") "zone",k
     call h5gopen_f(file_id,trim(groupname),group_id,error)
     dim2d(1)=zones(k)%mesh%Nx
     dim2d(2)=zones(k)%mesh%Nz
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     allocate(buffer(1:Nx,1:Nz),ballD(0:Nx+1,0:Nz+1),ballNu(0:Nx+1,0:Nz+1),ballChi(0:Nx+1,0:Nz+1),ballChie(0:Nx+1,0:Nz+1))
     call h5dopen_f(group_id,"ballooningD",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer,dim2d,error)
     ballD(1:Nx,1:Nz)=buffer
     call h5dopen_f(group_id,"ballooningNu",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer,dim2d,error)
     ballNu(1:Nx,1:Nz)=buffer
     call h5dopen_f(group_id,"ballooningChi",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer,dim2d,error)
     ballChi(1:Nx,1:Nz)=buffer
     call h5dopen_f(group_id,"ballooningChie",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer,dim2d,error)
     ballChie(1:Nx,1:Nz)=buffer
     call h5gclose_f(group_id,error) !close Zone
     zones(k)%species(0)%transport_perp%chi_p=zones(k)%species(0)%transport_perp%chi_p*ballChie
     zones(k)%species(0)%transport_perp%chi_t=zones(k)%species(0)%transport_perp%chi_t*ballChie
     do n=1,global_parameters%N_ions
        zones(k)%species(n)%transport_perp%D_p=zones(k)%species(n)%transport_perp%D_p*ballD
        zones(k)%species(n)%transport_perp%D_t=zones(k)%species(n)%transport_perp%D_t*ballD
        zones(k)%species(n)%transport_perp%nu_p=zones(k)%species(n)%transport_perp%nu_p*ballNu
        zones(k)%species(n)%transport_perp%nu_t=zones(k)%species(n)%transport_perp%nu_t*ballNu
        zones(k)%species(n)%transport_perp%chi_p=zones(k)%species(n)%transport_perp%chi_p*ballChi
        zones(k)%species(n)%transport_perp%chi_t=zones(k)%species(n)%transport_perp%chi_t*ballChi
     end do
     deallocate(buffer,ballD,ballNu,ballChi,ballChie)
  end do
  call h5fclose_f(file_id,error) !close file
  call h5close_f(error) !close hdf5
end subroutine load_ballooning_from_file
