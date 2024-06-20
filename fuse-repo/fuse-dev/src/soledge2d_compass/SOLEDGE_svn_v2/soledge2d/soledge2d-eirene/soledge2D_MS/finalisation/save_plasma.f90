subroutine save_plasma()
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
  real*8,allocatable :: buffer(:,:)
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
  do n=0,global_parameters%N_ions
     write(filename,"(A15,I0)") "Results/plasma_",n
     call h5fcreate_f(trim(filename),H5F_ACC_TRUNC_F,file_id,error) 
     !save properties
     dim1d(1)=1
     call h5screate_simple_f(1,dim1d,dspace_id,error)
     call h5dcreate_f(file_id,"charge",H5T_NATIVE_INTEGER,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_INTEGER,zones(1)%species(n)%charge,dim1d,Error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     call h5screate_simple_f(1,dim1d,dspace_id,error)
     call h5dcreate_f(file_id,"mass",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
     call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(1)%species(n)%element%mass,dim1d,error)
     call h5dclose_f(dataset_id,error)
     call h5sclose_f(dspace_id,error)
     !Save zones
     do k=1,global_parameters%N_zones
        write(groupname,"(A4,I0)") "zone",k
        call h5gcreate_f(file_id,trim(groupname),group_id,error)
        Nx=zones(k)%mesh%Nx
        Nz=zones(k)%mesh%Nz
        allocate(buffer(0:Nx+1,0:Nz+1))
        dim2d(1)=Nx+2
        dim2d(2)=Nz+2
        call h5screate_simple_f(2,dim2d,dspace_id,error)
        call h5dcreate_f(group_id,"density",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%species(n)%var(1)%density,dim2d,error)
        call h5dclose_f(dataset_id,error)
        call h5sclose_f(dspace_id,error)
        call h5screate_simple_f(2,dim2d,dspace_id,error)
        call h5dcreate_f(group_id,"Gamma",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%species(n)%var(1)%Gamma,dim2d,error)
        call h5dclose_f(dataset_id,error)
        call h5sclose_f(dspace_id,error)
        call h5screate_simple_f(2,dim2d,dspace_id,error)
        call h5dcreate_f(group_id,"temperature",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%species(n)%var(1)%temperature,dim2d,error)
        call h5dclose_f(dataset_id,error)
        call h5sclose_f(dspace_id,error)
        call h5screate_simple_f(2,dim2d,dspace_id,error)
        call h5dcreate_f(group_id,"mach",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%species(n)%var(1)%Mach,dim2d,error)
        call h5dclose_f(dataset_id,error)
        call h5sclose_f(dspace_id,error)
        call h5screate_simple_f(2,dim2d,dspace_id,error)
        call h5dcreate_f(group_id,"vpinch",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%species(n)%transport_perp%v_pinch,dim2d,error)
        call h5dclose_f(dataset_id,error)
        call h5sclose_f(dspace_id,error)
        if(flags%neutral_model.eq.2) then ! fluid neutrals
           if(n.eq.1) then
              call h5screate_simple_f(2,dim2d,dspace_id,error)
              call h5dcreate_f(group_id,"Ndensity",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
              call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%neutrals%density,dim2d,error)
              call h5dclose_f(dataset_id,error)
              call h5sclose_f(dspace_id,error)
           end if
        end if
        dim2d(1)=Nx
        dim2d(2)=Nz
        call h5screate_simple_f(2,dim2d,dspace_id,error)
        call h5dcreate_f(group_id,"rad",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%species(n)%sources%rad,dim2d,error)
        call h5dclose_f(dataset_id,error)
        call h5sclose_f(dspace_id,error)
        call h5screate_simple_f(2,dim2d,dspace_id,error)
        call h5dcreate_f(group_id,"alpham",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%species(n)%penalisation_memories%alpham,dim2d,error)
        call h5dclose_f(dataset_id,error)
        call h5sclose_f(dspace_id,error)
        call h5screate_simple_f(2,dim2d,dspace_id,error)
        call h5dcreate_f(group_id,"alphap",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%species(n)%penalisation_memories%alphap,dim2d,error)
        call h5dclose_f(dataset_id,error)
        call h5sclose_f(dspace_id,error)
        deallocate(buffer)
        call h5gclose_f(group_id,error)
     end do
     if(flags%use_triangles) then
        groupname='triangles'
        call h5gcreate_f(file_id,trim(groupname),group_id,error)
        dim1d=interp_data2%N_knots
        call h5screate_simple_f(1,dim1d,dspace_id,error)
        call h5dcreate_f(group_id,"density",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,interp_data2%knots_density(:,n,1),dim1d,error)
        call h5dclose_f(dataset_id,error)
        call h5sclose_f(dspace_id,error)        
        call h5screate_simple_f(1,dim1d,dspace_id,error)
        call h5dcreate_f(group_id,"velocity",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,interp_data2%knots_velocity(:,n,1),dim1d,error)
        call h5dclose_f(dataset_id,error)
        call h5sclose_f(dspace_id,error)        
        call h5screate_simple_f(1,dim1d,dspace_id,error)
        call h5dcreate_f(group_id,"temperature",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,interp_data2%knots_temperature(:,n,1),dim1d,error)
        call h5dclose_f(dataset_id,error)
        call h5sclose_f(dspace_id,error)        
        call h5screate_simple_f(1,dim1d,dspace_id,error)
        call h5dcreate_f(group_id,"radiation",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,interp_data2%knots_radiation(:,n),dim1d,error)
        call h5dclose_f(dataset_id,error)
        call h5sclose_f(dspace_id,error)        
        if(n.eq.0) then
           !electrons
           dim1d=interp_data2%N_knots
           call h5screate_simple_f(1,dim1d,dspace_id,error)
           call h5dcreate_f(group_id,"Zeff",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
           call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,interp_data2%knots_Zeff(:),dim1d,error)
           call h5dclose_f(dataset_id,error)
           call h5sclose_f(dspace_id,error)
           dim1d=interp_data2%N_triangles
           call h5screate_simple_f(1,dim1d,dspace_id,error)
           call h5dcreate_f(group_id,"NRad",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
           call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,interp_data2%tri_Srad(:,0,1),dim1d,error)
           call h5dclose_f(dataset_id,error)
           call h5sclose_f(dspace_id,error)        
        else
           if(global_parameters%ions_list(n,2).eq.1) then
              ! charge = 1 ion
              ind_spec=global_parameters%ions_list(n,1)
              dim1d=interp_data2%N_triangles
              call h5screate_simple_f(1,dim1d,dspace_id,error)
              call h5dcreate_f(group_id,"Sn",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
              call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,interp_data2%tri_Sn(:,ind_spec,1),dim1d,error)
              call h5dclose_f(dataset_id,error)
              call h5sclose_f(dspace_id,error)        
              call h5screate_simple_f(1,dim1d,dspace_id,error)
              call h5dcreate_f(group_id,"Nn",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
              call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,interp_data2%tri_Nn(:,ind_spec,1),dim1d,error)
              call h5dclose_f(dataset_id,error)
              call h5sclose_f(dspace_id,error)        
              call h5screate_simple_f(1,dim1d,dspace_id,error)
              call h5dcreate_f(group_id,"Nm",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
              call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,interp_data2%tri_Nm(:,ind_spec,1),dim1d,error)
              call h5dclose_f(dataset_id,error)
              call h5sclose_f(dspace_id,error)        
              call h5screate_simple_f(1,dim1d,dspace_id,error)
              call h5dcreate_f(group_id,"Tn",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
              call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,interp_data2%tri_Tn(:,ind_spec,1),dim1d,error)
              call h5dclose_f(dataset_id,error)
              call h5sclose_f(dspace_id,error)        
              call h5screate_simple_f(1,dim1d,dspace_id,error)
              call h5dcreate_f(group_id,"Tm",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
              call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,interp_data2%tri_Tm(:,ind_spec,1),dim1d,error)
              call h5dclose_f(dataset_id,error)
              call h5sclose_f(dspace_id,error)        
              call h5screate_simple_f(1,dim1d,dspace_id,error)
              call h5dcreate_f(group_id,"Nrad",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
              call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,interp_data2%tri_Srad(:,ind_spec,1),dim1d,error)
              call h5dclose_f(dataset_id,error)
              call h5sclose_f(dspace_id,error)        
           end if
        end if
        call h5gclose_f(group_id,error)
     end if
     call h5fclose_f(file_id,error)
  end do
  call h5close_f(error)
end subroutine save_plasma
