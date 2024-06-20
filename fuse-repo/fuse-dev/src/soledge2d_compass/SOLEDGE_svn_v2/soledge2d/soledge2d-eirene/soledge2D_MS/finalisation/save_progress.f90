subroutine save_progress(n_ite)
  use all_variables, only : zones, global_parameters, interp_data2, flags, drift_flags
  use hdf5
  implicit none
  integer*4,intent(in) :: n_ite
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
  integer(HSIZE_T) :: dim3d(3)
  !other variable
  integer*4 :: Nx,Nz
  real*8,allocatable :: buffer(:,:)
  integer*4 :: k,n,ind_spec
  logical :: dir_e
  integer*4 :: nsave
#include "compile_opt.inc"
#if GFORTRAN==1
  inquire(File='Evolution',exist=dir_e)
#endif
#if GFORTRAN==0
  inquire(Directory='Evolution',exist=dir_e)
#endif
  if(.not.dir_e) then
     call system("mkdir Evolution")
  end if
  if(global_parameters%N_iterations.gt.10) then
     if(modulo(dble(n_ite),dble(global_parameters%N_iterations)/real(global_parameters%N_save)).eq.0) then
        nsave = n_ite/(floor(dble(global_parameters%N_iterations)/real(global_parameters%N_save)))
        call h5open_f(error)
        do n=0,global_parameters%N_ions
           write(filename,"(I0,A8,I0)") nsave,"_plasma_",n
           call h5fcreate_f(trim("Evolution/"//filename),H5F_ACC_TRUNC_F,file_id,error) 
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
              if(n.eq.0) then
                 !electrons
                 if (drift_flags%solve_phi) then
                 call h5screate_simple_f(2,dim2d,dspace_id,error)
                 call h5dcreate_f(group_id,"phi",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
                 call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%electric_fields(1)%phi,dim2d,error)
                 call h5dclose_f(dataset_id,error)
                 call h5sclose_f(dspace_id,error)
		end if
              end if
              if(drift_flags%solve_drift) then
                 call h5screate_simple_f(2,dim2d,dspace_id,error)
                 call h5dcreate_f(group_id,"uEp",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
                 call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%species(n)%drifts%uEp,dim2d,error)
                 call h5dclose_f(dataset_id,error)
                 call h5sclose_f(dspace_id,error)
                 call h5screate_simple_f(2,dim2d,dspace_id,error)
                 call h5dcreate_f(group_id,"uEt",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
                 call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%species(n)%drifts%uEt,dim2d,error)
                 call h5dclose_f(dataset_id,error)
                 call h5sclose_f(dspace_id,error)
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
              if(n.eq.0) then
                 !electrons
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

        !currents
        if (drift_flags%solve_phi) then
        write(filename,"(I0,A9)") nsave,"_currents"
        call h5fcreate_f(trim("Evolution/"//filename),H5F_ACC_TRUNC_F,file_id,error) 
        !Save zones
        do k=1,global_parameters%N_zones
           write(groupname,"(A4,I0)") "zone",k
           call h5gcreate_f(file_id,trim(groupname),group_id,error)
           Nx=zones(k)%mesh%Nx
           Nz=zones(k)%mesh%Nz
           dim3d(1)=Nx
           dim3d(2)=Nz
           dim3d(3)=4
           call h5screate_simple_f(3,dim3d,dspace_id,error)
           call h5dcreate_f(group_id,"j_parallel",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
           call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%electric_fields(1)%j_parallel,dim3d,error)
           call h5dclose_f(dataset_id,error)
           call h5sclose_f(dspace_id,error)
           call h5screate_simple_f(3,dim3d,dspace_id,error)
           call h5dcreate_f(group_id,"j_adv_W",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
           call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%electric_fields(1)%j_para_adv_W,dim3d,error)
           call h5dclose_f(dataset_id,error)
           call h5sclose_f(dspace_id,error)
           call h5screate_simple_f(3,dim3d,dspace_id,error)
           call h5dcreate_f(group_id,"j_diff_W",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
           call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%electric_fields(1)%j_diff_W,dim3d,error)
           call h5dclose_f(dataset_id,error)
           call h5sclose_f(dspace_id,error)
           call h5screate_simple_f(3,dim3d,dspace_id,error)
           call h5dcreate_f(group_id,"j_pola",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
           call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%electric_fields(1)%j_perp,dim3d,error)
           call h5dclose_f(dataset_id,error)
           call h5sclose_f(dspace_id,error)
           call h5screate_simple_f(3,dim3d,dspace_id,error)
           call h5dcreate_f(group_id,"j_diam",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
           call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%species(1)%drifts%jdiam,dim3d,error)
           call h5dclose_f(dataset_id,error)
           call h5sclose_f(dspace_id,error)
           call h5screate_simple_f(3,dim3d,dspace_id,error)
           call h5dcreate_f(group_id,"j_ExB",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
           call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%species(1)%drifts%jExB,dim3d,error)
           call h5dclose_f(dataset_id,error)
           call h5sclose_f(dspace_id,error)
           call h5screate_simple_f(3,dim3d,dspace_id,error)
           call h5dcreate_f(group_id,"j_BxdB",H5T_NATIVE_DOUBLE,dspace_id,dataset_id,error)
           call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,zones(k)%species(1)%drifts%jBxDB,dim3d,error)
           call h5dclose_f(dataset_id,error)
           call h5sclose_f(dspace_id,error)
           call h5gclose_f(group_id,error)
        end do
        call h5fclose_f(file_id,error)
        end if ! drift_flags%solve_phi

        call h5close_f(error)
     end if
  end if
end subroutine save_progress
