subroutine load_plasma()
  use all_variables, only : zones, global_parameters
  use hdf5
  implicit none
  !HDF5 variable
  integer(HID_T) :: file_id
  integer(HID_T) :: group_id
  integer(HID_T) :: group2_id
  integer(HID_T) :: dataset_id
  integer :: error
  character(50) :: filename
  character(50) :: groupname
  integer(HSIZE_T) :: dim1d(1)
  integer(HSIZE_T) :: dim2d(2)
  !other variable
  integer*4 :: Nx,Nz
  real*8,allocatable :: buffer(:,:)
  integer*4 :: k,n
  call h5open_f(error)
  do n=0,global_parameters%N_ions
     write(filename,"(A15,I0)") "Results/plasma_",n
     call h5fopen_f(trim(filename),H5F_ACC_RDWR_F,file_id,error) 
     do k=1,global_parameters%N_zones
        write(groupname,"(A4,I0)") "zone",k
        call h5gopen_f(file_id,trim(groupname),group_id,error)
        Nx=zones(k)%mesh%Nx
        Nz=zones(k)%mesh%Nz
        allocate(buffer(0:Nx+1,0:Nz+1))
        dim2d(1)=Nx+2
        dim2d(2)=Nz+2
        call h5dopen_f(group_id,"density",dataset_id,error)
        call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer,dim2d,error)
        zones(k)%species(n)%var(1)%density=buffer
        call h5dclose_f(dataset_id,error)
        call h5dopen_f(group_id,"Gamma",dataset_id,error)
        call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer,dim2d,error)
        zones(k)%species(n)%var(1)%Gamma=buffer
        call h5dclose_f(dataset_id,error)
        call h5dopen_f(group_id,"temperature",dataset_id,error)
        call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer,dim2d,error)
        zones(k)%species(n)%var(1)%temperature=buffer
        call h5dclose_f(dataset_id,error)
        deallocate(buffer)
        allocate(buffer(1:Nx,1:Nz))
        dim2d(1)=Nx
        dim2d(2)=Nz
        call h5dopen_f(group_id,"alpham",dataset_id,error)
        call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer,dim2d,error)
        zones(k)%species(n)%penalisation_memories%alpham=buffer
        call h5dclose_f(dataset_id,error)
        call h5dopen_f(group_id,"alphap",dataset_id,error)
        call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer,dim2d,error)
        zones(k)%species(n)%penalisation_memories%alphap=buffer
        call h5dclose_f(dataset_id,error)
        deallocate(buffer)
        call h5gclose_f(group_id,error)
     end do
     call h5fclose_f(file_id,error)
  end do
  call h5close_f(error)
end subroutine load_plasma
