subroutine styx_init_sheath_model
  use eirmod_precision
  use eirmod_parmmod
  use eirmod_cpes
  use styx2eirene
  implicit none
  integer*4 :: ier,i

  include 'mpif.h'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! even if the PIC sheath model is not in use, incidence angles are calculated !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call mpi_barrier(mpi_comm_world,ier)
  call mpi_bcast(Nsou,1,mpi_integer,0,mpi_comm_world,ier)

  allocate(sheath1D(Nsou))
  do i=1,Nsou
    sheath1D(i)%alphaB=0._dp
    sheath1D(i)%tau=0._dp
    sheath1D(i)%ksi=0._dp
    sheath1D(i)%ksi0=0._dp
    sheath1D(i)%uparX=0._dp
    sheath1D(i)%uparY=0._dp
    sheath1D(i)%uparZ=0._dp
    sheath1D(i)%alpha_V=0._dp
    sheath1D(i)%beta_V=0._dp
    sheath1D(i)%hit_V=0
  enddo

  allocate(Bnu2xyz(3,3,Nsou),xyz2Bnu(3,3,Nsou))

  Bnu2xyz=0._dp
  xyz2Bnu=0._dp

! moved around -- check

  call mpi_barrier(mpi_comm_world,ier)
  call mpi_bcast(sheath_model,1,mpi_integer,0,mpi_comm_world,ier)

  if (sheath_model == 1) then
    call mpi_bcast(sheath1D_av%Ntau,1,mpi_integer,0,mpi_comm_world,ier)
    call mpi_bcast(sheath1D_av%Nksi,1,mpi_integer,0,mpi_comm_world,ier)
    call mpi_bcast(sheath1D_av%NalphaB,1,mpi_integer,0,mpi_comm_world,ier)

    if (my_pe /= 0) then
      do i=1,Nsou
        allocate(sheath1D(i)%E(sheath1D_av%Ntau,sheath1D_av%Nksi))
        allocate(sheath1D(i)%alpha(sheath1D_av%Ntau,sheath1D_av%Nksi))
        allocate(sheath1D(i)%beta(sheath1D_av%Ntau,sheath1D_av%Nksi))
      enddo
      allocate(sheath1D_av%tau(sheath1D_av%Ntau))
      allocate(sheath1D_av%ksi(sheath1D_av%Nksi))
    endif
  endif

  if (my_pe == 0) call styx_fill_sheath_data()

  call mpi_barrier(mpi_comm_world,ier)
       
  if (sheath_model == 1) then   
    do i=1,Nsou
      call mpi_bcast(sheath1D(i)%E,sheath1D_av%Ntau*sheath1D_av%Nksi,mpi_real8,0,mpi_comm_world,ier)
      call mpi_bcast(sheath1D(i)%alpha,sheath1D_av%Ntau*sheath1D_av%Nksi,mpi_real8,0,mpi_comm_world,ier)
      call mpi_bcast(sheath1D(i)%beta,sheath1D_av%Ntau*sheath1D_av%Nksi,mpi_real8,0,mpi_comm_world,ier)
    enddo
    call mpi_bcast(sheath1D_av%tau,sheath1D_av%Ntau,mpi_real8,0,mpi_comm_world,ier)
    call mpi_bcast(sheath1D_av%ksi,sheath1D_av%Nksi,mpi_real8,0,mpi_comm_world,ier)
  endif

  call mpi_bcast(Ntri_styx,1,mpi_integer,0,mpi_comm_world,ier)

  if (my_pe /= 0) then
    allocate(recsurfinv(3,Ntri_styx))
  endif

  call mpi_bcast(recsurfinv,Ntri_styx*3,mpi_integer,0,mpi_comm_world,ier)

  ! change of coordinates matrices (x,y,z <-> upar,n,uperp), for incidence angles
  call mpi_bcast(Bnu2xyz,9*nsou,mpi_real8,0,mpi_comm_world,ier)
  call mpi_bcast(xyz2Bnu,9*nsou,mpi_real8,0,mpi_comm_world,ier)
  

  do i=1,Nsou
    call mpi_bcast(sheath1D(i)%uparX,1,mpi_real8,0,mpi_comm_world,ier)
    call mpi_bcast(sheath1D(i)%uparY,1,mpi_real8,0,mpi_comm_world,ier)
    call mpi_bcast(sheath1D(i)%uparZ,1,mpi_real8,0,mpi_comm_world,ier)   
  enddo


end subroutine styx_init_sheath_model


