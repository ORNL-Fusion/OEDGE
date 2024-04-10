subroutine load_transport_coefficients()
  use all_variables, only : zones, global_parameters, flags, reference_parameters
  use hdf5
  use MradialFeedback
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
  integer*4 :: Nx,Nz,Nxshift
  real*8,allocatable :: buffer(:,:)
  integer*4 :: k,n,i
  real*8 :: D0 
  call h5open_f(error)
  call h5fopen_f("Results/transports_coefficients",H5F_ACC_RDWR_F,file_id,error) 
  do k=1,global_parameters%N_zones
     write(groupname,"(A4,I0)") "zone",k
     call h5gopen_f(file_id,trim(groupname),group_id,error)
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     allocate(buffer(0:Nx+1,0:Nz+1))
     dim2d(1)=Nx+2
     dim2d(2)=Nz+2
     call h5dopen_f(group_id,"D",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer,dim2d,error)
     do n=1,global_parameters%N_ions
        zones(k)%species(n)%transport_perp%D_p=buffer
        zones(k)%species(n)%transport_perp%nu_p=buffer
     end do
     call h5dclose_f(dataset_id,error)
     call h5dopen_f(group_id,"Chii",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer,dim2d,error)
     do n=1,global_parameters%N_ions
        zones(k)%species(n)%transport_perp%chi_p=buffer
     end do
     call h5dclose_f(dataset_id,error)
     call h5dopen_f(group_id,"Chie",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer,dim2d,error)
     zones(k)%species(0)%transport_perp%chi_p=buffer
     call h5dclose_f(dataset_id,error)
     deallocate(buffer)
     call h5gclose_f(group_id,error)
  end do
  call h5fclose_f(file_id,error)
  call h5close_f(error)
  if(flags%radialFeedback) then
     Nxshift=0
     do n=1,radialFeedbackData%Nzones
        k=radialFeedbackData%zone_profile(n)
        Nx=zones(k)%mesh%Nx
        !recompute new local diffusivity from flux and desired gradient
        Nz=radialFeedbackData%Nz
        D0=(reference_parameters%geometry%rs0**2.)/reference_parameters%fields%tau0
        do i=1,Nx
           radialFeedbackData%D(i+Nxshift)=zones(k)%species(1)%transport_perp%D_p(i,Nz)*D0
           radialFeedbackData%chie(i+Nxshift)=zones(k)%species(0)%transport_perp%chi_p(i,Nz)*D0
           radialFeedbackData%chii(i+Nxshift)=zones(k)%species(1)%transport_perp%chi_p(i,Nz)*D0
        end do
        Nxshift=Nxshift+Nx
     end do
  end if
end subroutine load_transport_coefficients
