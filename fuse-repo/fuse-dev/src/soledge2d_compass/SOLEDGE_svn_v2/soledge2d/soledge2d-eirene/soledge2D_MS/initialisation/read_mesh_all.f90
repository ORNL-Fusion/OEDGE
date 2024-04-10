subroutine read_mesh_all()
  use all_variables, only : global_parameters,zones,megazones
  use Mzone
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
  integer*4 :: k,i,j
  integer*4 :: Nx,Nz
  real*8,allocatable :: buffer(:,:)
  integer*4 :: isperiodic_int
  dim1d(1)=1
  call h5open_f(error)
  call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,error)
  !read N_zones
  call h5dopen_f(file_id,"NZones",dataset_id,error)
  call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,global_parameters%N_zones,dim1d,error)
  call h5dclose_f(dataset_id,error)
  allocate(zones(1:global_parameters%N_zones))

  !read N_megazones
  call h5dopen_f(file_id,"NMegazones",dataset_id,error)
  call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,global_parameters%N_megazones,dim1d,error)
  call h5dclose_f(dataset_id,error)
  allocate(megazones(1:global_parameters%N_megazones))

  do k=1,global_parameters%N_zones
     write(groupname,"(A4,I0)") "zone",k
     call h5gopen_f(file_id,trim(groupname),group_id,error)
     !read zone size
     call h5dopen_f(group_id,"Nx",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,zones(k)%mesh%Nx,dim1d,error)
     call h5dclose_f(dataset_id,error)
     call h5dopen_f(group_id,"Nz",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,zones(k)%mesh%Nz,dim1d,error)
     call h5dclose_f(dataset_id,error)
     call allocate_zone_mesh(zones(k)%mesh)
     !load mesh
     dim2d(1)=zones(k)%mesh%Nx
     dim2d(2)=zones(k)%mesh%Nz
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     allocate(buffer(1:Nx,1:Nz))
     call h5dopen_f(group_id,"x",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer,dim2d,error)
     zones(k)%mesh%x(1:Nx,1:Nz)=buffer
     call h5dclose_f(dataset_id,error)
     call h5dopen_f(group_id,"xm",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer,dim2d,error)
     zones(k)%mesh%x_minus_1half(1:Nx,1:Nz)=buffer
     call h5dclose_f(dataset_id,error)
     call h5dopen_f(group_id,"xp",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer,dim2d,error)
     zones(k)%mesh%x_plus_1half(1:Nx,1:Nz)=buffer
     call h5dclose_f(dataset_id,error)
     call h5dopen_f(group_id,"z",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer,dim2d,error)
     zones(k)%mesh%z(1:Nx,1:Nz)=buffer
     call h5dclose_f(dataset_id,error)
     call h5dopen_f(group_id,"zm",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer,dim2d,error)
     zones(k)%mesh%z_minus_1half(1:Nx,1:Nz)=buffer
     call h5dclose_f(dataset_id,error)
     call h5dopen_f(group_id,"zp",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer,dim2d,error)
     zones(k)%mesh%z_plus_1half(1:Nx,1:Nz)=buffer
     call h5dclose_f(dataset_id,error)
     !load mesh boundary
     call h5dopen_f(group_id,"xmin",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%mesh%xmin,dim1d,error)
     call h5dclose_f(dataset_id,error)
     call h5dopen_f(group_id,"xmax",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%mesh%xmax,dim1d,error)
     call h5dclose_f(dataset_id,error)
     call h5dopen_f(group_id,"zmin",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%mesh%zmin,dim1d,error)
     call h5dclose_f(dataset_id,error)
     call h5dopen_f(group_id,"zmax",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%mesh%zmax,dim1d,error)
     call h5dclose_f(dataset_id,error)

!!$     do i=1,Nx
!!$        do j=1,Nz
!!$           zones(k)%mesh%z(i,j)=j
!!$           zones(k)%mesh%x(i,j)=i
!!$           zones(k)%mesh%z_minus_1half(i,j)=j-0.5d0
!!$           zones(k)%mesh%z_plus_1half(i,j)=j+0.5d0
!!$           zones(k)%mesh%x_minus_1half(i,j)=i-0.5d0
!!$           zones(k)%mesh%x_plus_1half(i,j)=i+0.5d0
!!$        end do
!!$     end do
!!$     zones(k)%mesh%xmin=0.5d0
!!$     zones(k)%mesh%xmax=Nx+0.5
!!$     zones(k)%mesh%zmin=0.5d0
!!$     zones(k)%mesh%zmax=Nz+0.5

     !load neighbors info
     call h5gopen_f(group_id,'Neighbors',group2_id,error)
     call h5dopen_f(group2_id,"North",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,zones(k)%Neighbors(N_North),dim1d,error)
     call h5dclose_f(dataset_id,error)
     call h5dopen_f(group2_id,"South",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,zones(k)%Neighbors(N_South),dim1d,error)
     call h5dclose_f(dataset_id,error)
     call h5dopen_f(group2_id,"East",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,zones(k)%Neighbors(N_East),dim1d,error)
     call h5dclose_f(dataset_id,error)
     call h5dopen_f(group2_id,"West",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,zones(k)%Neighbors(N_West),dim1d,error)
     call h5dclose_f(dataset_id,error)
     call h5gclose_f(group2_id,error) !close Neighbors
     !load MagNeighbors info
     call h5gopen_f(group_id,'MagNeighbors',group2_id,error)
     call h5dopen_f(group2_id,"North",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,zones(k)%MagNeighbors(N_North),dim1d,error)
     call h5dclose_f(dataset_id,error)
     call h5dopen_f(group2_id,"South",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,zones(k)%MagNeighbors(N_South),dim1d,error)
     call h5dclose_f(dataset_id,error)
     call h5dopen_f(group2_id,"East",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,zones(k)%MagNeighbors(N_East),dim1d,error)
     call h5dclose_f(dataset_id,error)
     call h5dopen_f(group2_id,"West",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,zones(k)%MagNeighbors(N_West),dim1d,error)
     call h5dclose_f(dataset_id,error)
     call h5gclose_f(group2_id,error) !close MagNeighbors
     call h5gclose_f(group_id,error) !close Zone
     deallocate(buffer)
  end do

  do k=1,global_parameters%N_megazones
     write(groupname,"(A8,I0)") "megazone",k
     call h5gopen_f(file_id,trim(groupname),group_id,error)
     dim1d(1)=1
     call h5dopen_f(group_id,"size",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,megazones(k)%size,dim1d,error)
     call h5dclose_f(dataset_id,error)
     allocate(megazones(k)%zone_number(1:megazones(k)%size))
     dim1d(1)=megazones(k)%size
     call h5dopen_f(group_id,"configuration",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,megazones(k)%zone_number,dim1d,error)
     call h5dclose_f(dataset_id,error)
     dim1d(1)=1
     call h5dopen_f(group_id,"isperiodic",dataset_id,error)
     call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,isperiodic_int,dim1d,error)
     call h5dclose_f(dataset_id,error)
     if(isperiodic_int.eq.1) then
        megazones(k)%is_periodic=.true.
     else
        megazones(k)%is_periodic=.false.
     end if
     call h5gclose_f(group_id,error) !close Megazone
  end do

  call h5fclose_f(file_id,error) !close file
  call h5close_f(error) !close hdf5
end subroutine read_mesh_all


subroutine allocate_zone_mesh(mesh)
  use Mgeometry
  implicit none
  Type(Tmesh),intent(inout) :: mesh
  integer*4 ::  Nx,Nz
  Nx=mesh%Nx
  Nz=mesh%Nz
  allocate(mesh%x(0:Nx+1,0:Nz+1))
  allocate(mesh%x_minus_1half(0:Nx+1,0:Nz+1))
  allocate(mesh%x_plus_1half(0:Nx+1,0:Nz+1))
  allocate(mesh%z(0:Nx+1,0:Nz+1))
  allocate(mesh%z_minus_1half(0:Nx+1,0:Nz+1))
  allocate(mesh%z_plus_1half(0:Nx+1,0:Nz+1))
end subroutine allocate_zone_mesh
