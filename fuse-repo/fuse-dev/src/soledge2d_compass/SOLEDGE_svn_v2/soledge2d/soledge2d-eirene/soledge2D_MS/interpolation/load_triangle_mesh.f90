subroutine load_triangle_mesh()
  use all_variables, only : global_parameters, Interp_data2, zones
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
  integer(HSIZE_T) :: dim3d(3)
  !other variable
  integer*4 :: Nx,Nz
  real*8,allocatable :: buffer1(:)
  real*8,allocatable :: buffer2(:,:)
  integer*4,allocatable :: buffer1i(:)
  integer*4,allocatable :: buffer2i(:,:)
  integer*4,allocatable :: buffer3i(:,:,:)
  integer*4 :: k,nknot,n,ti,tj,tk,num_tri
  call h5open_f(error)
  filename='triangles.h5'
  dim1d=1
  call h5fopen_f(trim(filename),H5F_ACC_RDWR_F,file_id,error) 
  call h5dopen_f(file_id,"Ntriangles",dataset_id,error)
  call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,Interp_data2%N_triangles,dim1d,error)
  call h5dclose_f(dataset_id,error)
  call h5dopen_f(file_id,"Nknots",dataset_id,error)
  call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,Interp_data2%N_knots,dim1d,error)
  call h5dclose_f(dataset_id,error)
  call allocate_interpolated_data2()

  groupname='knots'
  allocate(buffer1(1:Interp_data2%N_knots))
  allocate(buffer1i(1:Interp_data2%N_knots))
  allocate(buffer2i(1:Interp_data2%N_knots,1:3))
  allocate(buffer3i(1:Interp_data2%N_knots,1:8,1:3))
  dim1d=Interp_data2%N_knots
  call h5gopen_f(file_id,trim(groupname),group_id,error)
  call h5dopen_f(group_id,"R",dataset_id,error)
  call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer1,dim1d,error)
  call h5dclose_f(dataset_id,error)
  do k=1,Interp_data2%N_knots
     Interp_data2%knots_R(k)=buffer1(k)
  end do
  call h5dopen_f(group_id,"Z",dataset_id,error)
  call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer1,dim1d,error)
  call h5dclose_f(dataset_id,error)
  do k=1,Interp_data2%N_knots
     Interp_data2%knots_Z(k)=buffer1(k)
  end do
  call h5dopen_f(group_id,"pass",dataset_id,error)
  call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,buffer1i,dim1d,error)
  call h5dclose_f(dataset_id,error)
  do k=1,Interp_data2%N_knots
     Interp_data2%knots_interp_points(k)%pass=buffer1i(k)
  end do
  call h5dopen_f(group_id,"nsol",dataset_id,error)
  call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,buffer1i,dim1d,error)
  call h5dclose_f(dataset_id,error)
  do k=1,Interp_data2%N_knots
     Interp_data2%knots_interp_points(k)%n_soledge=buffer1i(k)
  end do
  call h5dopen_f(group_id,"neir",dataset_id,error)
  call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,buffer1i,dim1d,error)
  call h5dclose_f(dataset_id,error)
  do k=1,Interp_data2%N_knots
     Interp_data2%knots_interp_points(k)%n_eirene=buffer1i(k)
  end do
  dim3d(1)=Interp_data2%N_knots
  dim3d(2)=8
  dim3d(3)=3
  call h5dopen_f(group_id,"sol",dataset_id,error)
  call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,buffer3i,dim1d,error)
  call h5dclose_f(dataset_id,error)
  do k=1,Interp_data2%N_knots
     Interp_data2%knots_interp_points(k)%sol=buffer3i(k,:,:)
  end do
  dim2d(1)=Interp_data2%N_knots
  dim2d(2)=3
  call h5dopen_f(group_id,"eir",dataset_id,error)
  call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,buffer2i,dim1d,error)
  call h5dclose_f(dataset_id,error)
  do k=1,Interp_data2%N_knots
     Interp_data2%knots_interp_points(k)%eir=buffer2i(k,:)
  end do
  call h5gclose_f(group_id,error)


  groupname='wall'
  call h5gopen_f(file_id,trim(groupname),group_id,error)
  dim1d=1
  call h5dopen_f(group_id,"Ntriangles",dataset_id,error)
  call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,Interp_data2%wall_data%N_triangles_on_wall,dim1d,error)
  call h5dclose_f(dataset_id,error)
  call allocate_wall_data()
  dim2d(1)=Interp_data2%wall_data%N_triangles_on_wall
  dim2d(2)=4
  call h5dopen_f(group_id,"back_interp",dataset_id,error)
  call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,Interp_data2%wall_data%back_interp_on_wall,dim2d,error)
  call h5dclose_f(dataset_id,error)
  dim1d=Interp_data2%wall_data%N_triangles_on_wall
  dim2d(1)=Interp_data2%wall_data%N_triangles_on_wall
  dim2d(2)=3
  deallocate(buffer1)
  deallocate(buffer2i)
  allocate(buffer1(1:Interp_data2%wall_data%N_triangles_on_wall))
  allocate(buffer2i(1:Interp_data2%wall_data%N_triangles_on_wall,1:3))
  call h5dopen_f(group_id,"north_s2d_to_use",dataset_id,error)
  call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,buffer2i,dim2d,error)
  call h5dclose_f(dataset_id,error)
  call h5dopen_f(group_id,"north_s2d_weight",dataset_id,error)
  call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer1,dim1d,error)
  call h5dclose_f(dataset_id,error)
  do k=1,Interp_data2%wall_data%N_triangles_on_wall
     do n=1,3
        Interp_data2%wall_data%s2d_to_use(k,1,n)=buffer2i(k,n)
     end do
     Interp_data2%wall_data%weight_s2d_to_use(k,1)=buffer1(k)
  end do
  call h5dopen_f(group_id,"south_s2d_to_use",dataset_id,error)
  call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,buffer2i,dim2d,error)
  call h5dclose_f(dataset_id,error)
  call h5dopen_f(group_id,"south_s2d_weight",dataset_id,error)
  call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer1,dim1d,error)
  call h5dclose_f(dataset_id,error)
  do k=1,Interp_data2%wall_data%N_triangles_on_wall
     do n=1,3
        Interp_data2%wall_data%s2d_to_use(k,2,n)=buffer2i(k,n)
     end do
     Interp_data2%wall_data%weight_s2d_to_use(k,2)=buffer1(k)
  end do
  call h5dopen_f(group_id,"east_s2d_to_use",dataset_id,error)
  call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,buffer2i,dim2d,error)
  call h5dclose_f(dataset_id,error)
  call h5dopen_f(group_id,"east_s2d_weight",dataset_id,error)
  call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer1,dim1d,error)
  call h5dclose_f(dataset_id,error)
  do k=1,Interp_data2%wall_data%N_triangles_on_wall
     do n=1,3
        Interp_data2%wall_data%s2d_to_use(k,3,n)=buffer2i(k,n)
     end do
     Interp_data2%wall_data%weight_s2d_to_use(k,3)=buffer1(k)
  end do
  call h5dopen_f(group_id,"west_s2d_to_use",dataset_id,error)
  call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,buffer2i,dim2d,error)
  call h5dclose_f(dataset_id,error)
  call h5dopen_f(group_id,"west_s2d_weight",dataset_id,error)
  call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,buffer1,dim1d,error)
  call h5dclose_f(dataset_id,error)
  do k=1,Interp_data2%wall_data%N_triangles_on_wall
     do n=1,3
        Interp_data2%wall_data%s2d_to_use(k,4,n)=buffer2i(k,n)
     end do
     Interp_data2%wall_data%weight_s2d_to_use(k,4)=buffer1(k)
  end do
  call h5gclose_f(group_id,error)

  groupname='triangles'
  call h5gopen_f(file_id,trim(groupname),group_id,error)
  dim2d(1)=Interp_data2%N_triangles
  dim2d(2)=3
  deallocate(buffer2i)
  allocate(buffer2i(1:Interp_data2%N_triangles,1:3))
  call h5dopen_f(group_id,"tri_knots",dataset_id,error)
  call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,Interp_data2%tri_knots,dim2d,error)
  call h5dclose_f(dataset_id,error)
  call h5dopen_f(group_id,"type_face",dataset_id,error)
  call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,Interp_data2%type_face,dim2d,error)
  call h5dclose_f(dataset_id,error)
  call h5dopen_f(group_id,"back_interp",dataset_id,error)
  call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,buffer2i,dim2d,error)
  call h5dclose_f(dataset_id,error)
  allocate(interp_data2%zones_data(1:global_parameters%N_zones))
  do k=1,global_parameters%N_zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     allocate(interp_data2%zones_data(k)%triangles(1:Nx,1:Nz,1:6))
     allocate(interp_data2%zones_data(k)%num_triangles(1:Nx,1:Nz))
     interp_data2%zones_data(k)%num_triangles=0
  end do
  do k=1,interp_data2%N_triangles
     ti=buffer2i(k,2)
     tj=buffer2i(k,3)
     tk=buffer2i(k,1)
     num_tri=interp_data2%zones_data(tk)%num_triangles(ti,tj)+1
     interp_data2%zones_data(tk)%triangles(ti,tj,num_tri)=k
     interp_data2%zones_data(tk)%num_triangles(ti,tj)=num_tri
  end do

  call h5gclose_f(group_id,error)

  call h5fclose_f(file_id,error)
  call H5close_f(error)
end subroutine load_triangle_mesh
