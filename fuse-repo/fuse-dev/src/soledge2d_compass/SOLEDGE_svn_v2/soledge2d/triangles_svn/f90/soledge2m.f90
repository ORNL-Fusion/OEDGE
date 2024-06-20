!##############################################################################!
!##    Program SOLEDGE2D                                                     ##!
!##    Version: 2.010                                                        ##!
!##    Date: 26, January 2012                                                ##!
!##    Author: Hugo Bufferand (M2P2 Laboratory, Marseille)                   ##!
!##    Contact: hugo.bufferand@gmail.com                                     ##!
!##############################################################################!

program soledge2m

  use domain_types
  use solver_types
  use N_solver
  use G_solver
  use Te_solver
  use Ti_solver
  use metrix
  use equimag
  use implicit
  use jacobs
  use initialisation
  use multidomain
  use saving
  use monitors
  use sources_perp
  use sources_para
  use derivatives
  use CFL_dt
  use center
  use neutrals
  use penalisation
  use RungeKutta
  use interpolation_soleir
  use hdf5 !DJM

  !############################################################################!
  !####             Variables declaration                        ##############!
  !############################################################################!

  implicit none
  ! Variables associated to each domain containing 
  ! mesh and plasma fields (n,G,Te,Ti)
  Type(Zone),allocatable:: Zones(:)   ! TYPE: see "domain_types.f90"

  ! Variables (solvers) containing perpendicular and parallel
  ! source terms for all equations and for each domain
  Type(dyn),allocatable:: dyns(:)   ! TYPE: see "solver_types.f90"

  ! Variables containing tridiagonal matrix and source
  ! vector for Te equation for each domain...
  Type(tridiagonal),allocatable:: dynTes(:) ! TYPE: see "solver_types.f90"
  ! ... and for each parallel topological structure
  Type(super_tridiagonal),allocatable :: dynsParae(:) ! TYPE: see "solver_types.f90"
  ! The same for Ti
  Type(tridiagonal),allocatable:: dynTis(:)
  Type(super_tridiagonal),allocatable :: dynsParai(:)
  ! The same for n
  Type(tridiagonal),allocatable:: dynNs(:)
  Type(super_tridiagonal),allocatable :: SdynNs(:)
  ! The same for Gamma
  Type(tridiagonal),allocatable:: dynGs(:)
  Type(super_tridiagonal),allocatable :: SdynGs(:)

  ! Variables containing Magnetic fields for each domain
  Type(MagField),allocatable :: MagFields(:)  !TYPE: see "equimag.f90"

  ! Variables containing Jacobians for each domain
  Type(Jacobian),allocatable :: Jacobians(:)  !TYPE: see "jacobs.f90"

  ! Variables containing Metric coefficients for each domain
  Type(metric),allocatable :: Metrics(:)     !TYPE: see "metrix.f90"

  ! Variables containing correction terms (mostly for Ti equation)
  Type(corrections),allocatable :: Cors(:)   !TYPE: see "solver_types.f90"

  ! Variables containing values for the mesh at core center
  ! (this is not fully implemented) ####CAREFUL####
  Type(Centre) :: centre_   !TYPE: see "domain_types.f90"

  ! Variable that contain information about multidomain
  ! topology in the parallel direction
  ! that is in fact the definition of parallel topological structures
  Type(config_parallel):: config_para     !TYPE: see "domain_types.f90"

  ! Sources used at the end of the 4 Runge-Kutta steps to compute fields at time n+1
  Type(Sources),allocatable :: RKsources(:)
  integer*4 RKstep  ! index of RK step (from 1 to 4)
  integer*4 isFinalRK ! boolean that tells whethet RK sources have been computed
  real*8 RKdt ! intermediate time step

  ! For neutrals field computation
  Type(Flux_n),allocatable :: Flux_ns(:)  !TYPE: see "neutrals.f90"

  integer*4 ntps,i,j,k,n !various index
  real*8 CFL,dt,tempus !CFL coefficient, time step, time
  integer*4 SolOK,normal_zone_number !various flags
  integer*4 limdir 
  real*8 resN,resG,resTe,resTi,resNn !residuals (absolute error)
  ! structure encompassing informations about the run such as
  ! number of iterations... contained in input.txt file
  Type(data) :: dat           !TYPE: see "domain_types.f90"
  integer*8 index1
  integer*4 Nx,Nz
  integer*4 flag_erre 
  real*8, dimension(:), ALLOCATABLE :: locdt !local CFL time step for each domain
  integer*4, dimension(:), ALLOCATABLE :: loclimdir !critical explicit term for time step
  real*8,allocatable :: Source(:)
  real*8 Stot

  ! DJM 27/1/12
  integer(HID_T) :: file_id            ! File identifier
  integer(HID_T) :: dset_id            ! Dataset identifier
  integer(HID_T) :: dspace_id          ! Dataspace identifier
  real*8, allocatable :: time_array(:) ! DJM 27/1/12, time array
  integer :: error ! Error flag
  integer(HSIZE_T), dimension(1) :: dim1d ! 1D dimensions
  character(len=9), parameter :: outfilename = "s2dout.h5" ! Output file name

  Type(InterpData1),allocatable::  InterpData1s(:)
  Type(InterpData2) InterpData2_

  !############################################################################!
  !####             Beginning of the main program                ##############!
  !############################################################################!


  ! Lecture du fichier d'input (input.txt)
  ! The information are stored in structure "dat" 
  ! that will be broadcast to each domain
  call load_param(dat)   !SR: see "initialisation.f90"

  ! DJM 27/1/12, iniitialise the time array
  allocate(time_array(1:dat%nite))

  ! allocation of the tables according to the number of domains
  allocate(Zones(dat%Nzones))  
  allocate(locdt(dat%Nzones))  
  allocate(loclimdir(dat%Nzones))  
  allocate(dyns(dat%Nzones))  
  allocate(dynTes(dat%Nzones))
  allocate(dynTis(dat%Nzones))
  allocate(dynNs(dat%Nzones))
  allocate(dynGs(dat%Nzones))
  allocate(MagFields(dat%Nzones))
  allocate(Jacobians(dat%Nzones))
  allocate(Metrics(dat%Nzones))
  allocate(Cors(dat%Nzones))
  allocate(RKsources(dat%Nzones))
  allocate(Flux_ns(dat%NZones))
  allocate(Source(dat%Nzones))
  allocate(InterpData1s(dat%Nzones))

  ! Initialisation of the parallel topology of the domain
  call init_config_para(config_para)  !SR: see "initialisation.f90"
  write(*,*) 'init config para OK'

  ! Initialisation of the zones including mesh
  call init_zones(Zones,dat%Nzones,centre_,config_para,dat)  !SR: see "initialisation.f90"
  write(*,*) 'init zone:   OK'

  call init_neutrals_flux(Zones,dat%NZones,Flux_ns)
  write(*,*) 'Init neutrals: OK'

  ! Initialisation of RK sources tables
  call init_RK(Zones,dat%NZones,RKsources)
  write(*,*) 'init RK sources OK'

  if(dat%isSLAB.eq.0) then 
     ! if the user does not choose Slab option
     ! Loading Magnetic field values  
     call init_MagFields(Zones,dat%Nzones,MagFields)   !SR: see "equimag.f90"
     call load_mag(Zones,dat%Nzones,MagFields,dat)     !SR: see "equimag.f90"
     write(*,*) 'init Mag:    OK'
     ! Calculation of Jacobians
     call init_jacobians(Jacobians,Zones,dat%Nzones)   !SR: see "jacobs.f90"
     call calcul_jacobians(Jacobians,Zones,dat%Nzones,MagFields)   !SR: see "jacobs.f90"
     !     call write_jac(JAcobians,Zones,dat%Nzones)
     write(*,*) 'init Jac:    OK'
     ! Computation of metric coefficients
     call calcul_metric(Metrics,Jacobians,MagFields,Zones,dat%Nzones,dat)   !SR: see "metrix.f90"
  else
     ! Metric coefficients are computed from 
     ! Slab data specified in input.txt file
     call calcul_metric_SLAB(Metrics,MagFields,Zones,dat%Nzones,dat%qref,dat)  !SR: see "metrix.f90"
  end if
  call write_met(Metrics,Zones,dat%Nzones)
  write(*,*) 'Metric      : OK'

  ! Computation of the plasma surface on the core side 
  ! (used to compute power flux from exhaust power)
  call compute_Score(Zones,Metrics,dat%NZones,dat,config_para)  !SR: see "multidomain.f90"

  ! Initialisation of the plasma / load of initial values
  call init_plasma(Zones,dat%Nzones,dat,centre_,Metrics,Flux_ns)   !SR: see "initialisation.f90"
  write(*,*) 'Init plasma : OK'

  ! Initialisation of the solvers, sources, matrix
  ! and corrections for each zone
  call init_solver(dyns,Zones,dat%Nzones,dat%eta)    !SR: see "initialisation.f90"
  write(*,*) 'Init solver : OK'
  call init_corrections(Zones,dat%Nzones,Cors)       !SR: see "initialisation.f90"
  !n
  call init_super_tridiag(config_para,SdynNs,Zones,dat%Nzones)    !SR: see "implicit.f90"
  call init_tridiag(dynNs,Zones,dat%Nzones)                       !SR: see "implicit.f90"
  !Gamma
  call init_super_tridiag(config_para,SdynGs,Zones,dat%Nzones)
  call init_tridiag(dynGs,Zones,dat%Nzones)
  !Te
  call init_super_tridiag(config_para,dynsParae,Zones,dat%Nzones)
  call init_tridiag(dynTes,Zones,dat%Nzones)  
  !Ti
  call init_super_tridiag(config_para,dynsParai,Zones,dat%Nzones)
  call init_tridiag(dynTis,Zones,dat%Nzones)  
  write(*,*) 'Init systems: OK'

  call init_interp(Zones,dat%NZones,InterpData1s,InterpData2_)

  call compute_Dmean(Zones,Metrics,dat%NZones,dat,config_para)

  ! broadcasting values on domain boundaries
  do k=1,dat%Nzones
     call load_neighboring_values(Zones(k),Zones,dat%Nzones,dat,centre_,Metrics)  !SR: see "multidomain.f90"
  end do

  ! Creation of the file where residues are stored
  open(unit=20,file='residus.txt',status='unknown')

  ! Creation of the file where errors are saved
  open(unit=40,file='soledge2m.log',status='unknown')

  write(*,212,advance='no')

  ! DJM 27/1/12, if starting from archive, set tempus to zero, otherwise set to
  ! the last time step of the last archive:
  if(dat%restart.eq.0) then
     tempus=0.
  else
     dim1d(1)=1
     call h5open_f(error)
     call h5fopen_f(outfilename,H5F_ACC_RDWR_F,file_id,error)
     call h5dopen_f(file_id, "endtime", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tempus, dim1d, error)
     call h5dclose_f(dset_id, error)
     call h5fclose_f(file_id, error)
     call h5close_f(error)
  end if
  CFL=dat%CFL0


  !$OMP PARALLEL private(k,ntps) default(SHARED)
  !############################################################################!
  !####             Beginning of the parallel loop               ##############!
  !############################################################################!
  do ntps=1,dat%nite
     ! Time step computation. For each zone, a time step is calculated considering
     ! CFL conditions for advection and a minimal time step for explicit perp
     ! diffusion. The minimum of these times is considered.
     !$OMP DO
     do k=1,dat%Nzones
        call calcul_dt(dyns(k),Zones,dat%Nzones,locdt(k),loclimdir(k),dat,Metrics)   !SR: see "CFL_dt.f90"
     end do
     !$OMP END DO
     !$OMP BARRIER
     !$OMP MASTER
     ! Looking for the minimum
     index1 = MINLOC(locdt,1)
     dt = locdt(index1)
     limdir = loclimdir(index1)
     ! Condition for no NaN result

     if(CFL.lt.dat%CFL0) then ! let's try with increased CFL
        CFL=CFL+(dat%CFL0-CFL)/20.d0
     end if
     !$OMP END MASTER
     !$OMP BARRIER
     ! One tries to solve the equation until one finds no NaN.
     ! One aborts if CFL is too small and reaches a critical value
     SolOK=0 !if one finds NaN in solution at the end, this flag remains to 0
     do while(SolOK.eq.0)
        !$OMP MASTER
        dt = dt*CFL     
        !$OMP END MASTER

        if(dat%isRK4.eq.1) then
           ! LOOP on Runge-Kutta steps (RK4)
           do RKstep=1,4
              !$OMP BARRIER
              !$OMP DO
              do k=1,dat%Nzones
                 ! Initializing fields and flux before intermediate time step iteration
                 ! At RKstep=4, isFinalRK flag is switched to 1
                 call preRK(Zones(k),RKsources(k),RKdt,dt,isFinalRK,RKstep)   !SR: see "RungeKutta.f90"
              end do
              !$OMP END DO

              !$OMP BARRIER
              !$OMP DO
              do k=1,dat%Nzones
                 ! Broadcast of values on domain boundary (for fluxes computation - plasma(2))
                 call load_neighboring_values(Zones(k),Zones,dat%Nzones,dat,centre_,Metrics)
                 ! Compute some derivatives
                 call calcul_deriv(dyns(k),Zones,dat%Nzones)  !(is it used ####CAREFUL####)   !SR: see "derivatives.f90"
              end do
              !$OMP END DO
              !$OMP BARRIER

              !$OMP MASTER
              ! Broadcast of the derivatives
              call load_NeighDerAll(dyns,Zones,dat%Nzones) !(is it used ####CAREFUL####)    !SR: see "derivatives.f90"

              !$OMP END MASTER
              !$OMP BARRIER

              ! ###RK4 Computing of the fluxes from intermediate time steps stored in plasma(2) ###
              ! The perpendicular explicit sources are calculated (for all equations)
              !$OMP DO 
              do k=1,dat%Nzones 
                 call calcul_SourcePerp(dyns(k),Zones,dat%Nzones,Metrics)  !SR: see "sources_perp.f90"
              end do
              !$OMP END DO
              !OMP BARRIER

              ! The parallel explicit sources are calculated (for all equations)
              !$OMP DO
              do k=1,dat%NZones
                 call calcul_Sourcepara(dyns(k),Zones,dat%NZones,Metrics,dat%isSlab,dat)  !SR: see "sources_para.f90"
              end do
              !$OMP END DO
              !$OMP BARRIER


              !$OMP DO
              do k=1,dat%NZones
                 ! Determining the value on the fluxes that must be taken into account for the RK4 step
                 call  computeFluxRK(Zones(k),dyns(k),RKsources(k),RKstep)  !SR: see "RungeKutta.f90"
              end do
              !$OMP END DO
              !$OMP BARRIER           


              ! Matrix and source vector computation for density and velocity equations
              !$OMP DO
              do k=1,dat%NZones
                 call calcul_matN(dyns(k),dynNs(k),Zones,dat%Nzones,RKdt,Metrics,dat%eta,dat)  !SR: see "N_solver.f90"
                 call calcul_sourceN(dyns(k),dynNs(k),Zones,dat%Nzones,RKdt,dat,dat%eta,RKsources(k),isFinalRK)           !SR: see "N_solver.f90"
                 call calcul_matG(dyns(k),dynGs(k),Zones,dat%Nzones,RKdt,Metrics,dat%eta)  !SR: see "G_solver.f90"
                 call calcul_sourceG(dyns(k),dynGs(k),Zones,dat%Nzones,RKdt,dat,RKsources(k),isFinalRK)           !SR: see "G_solver.f90"
              end do
              !$OMP END DO
              !$OMP BARRIER

              ! Solving density and velocity
              !$OMP DO
              do k=1,config_para%zone_para_number
                 ! --- Density ----
                 ! Concatenating zones matrix for solving parallel topological structures
                 call concat_N(dynNs,config_para,SdynNs(k),Zones,dat%Nzones,Metrics)                !SR: see "N_solver.f90"
                 ! Inverting the tridiagonal matrix fo Te
                 call invert_superT(SdynNs(k))                                                      !SR: see "implicit.f90"
                 ! Un concatenating to attribute Te to each Zone
                 call unconcat_N(SdynNs(k),config_para,dyns,Zones,dat%Nzones,flag_erre,Metrics)     !SR: see "N_solver.f90"
                 ! --- Velocity ---
                 ! Concatenating zones matrix for solving parallel topological structures
                 call concat_G(dynGs,config_para,SdynGs(k),Zones,dat%Nzones,Metrics)                !SR: see "G_solver.f90"
                 ! Inverting the tridiagonal matrix fo Ti
                 call invert_superT(SdynGs(k))                                                      !SR: see "implicit.f90"
                 ! Un concatenating to attribute Ti to each Zone
                 call unconcat_G(SdynGs(k),config_para,dyns,Zones,dat%Nzones,flag_erre,Metrics)     !SR: see "G_solver.f90"
              end do
              !$OMP END DO
              !$OMP BARRIER

              !$OMP DO 
              do k=1,dat%Nzones 
                 ! Computation of some convective term in the poloidal plane from implicitely solved density
                 ! This routine modifies the dyn%sources_perp tables 
                 call calc_source_perp_T_imp(dyns(k),Zones,dat%NZones,Metrics)                      !SR: see "sources_perp.f90"
              end do
              !$OMP END DO
              !OMP BARRIER

              ! Matrix and source vector computation for Te and Ti equations
              ! (this is optional: if flag isoT=1, the plasma is isothermal and Te=Te0, Ti=Ti0)
              if(dat%isoT.eq.0) then
                 !$OMP DO
                 do k=1,dat%Nzones
                    ! Computation of matrix for solving electronic temperature parallel diff
                    call calcul_matTe(dyns(k),dynTes(k),Zones,dat%Nzones,RKdt,Metrics,dat)                    !SR: see "Te_solver.f90"
                    ! Computation of source term for solving electronic temperature parallel diff
                    call calcul_sourceTe(dyns(k),dynTes(k),Zones,dat%Nzones,RKdt,dat,Cors(k),Metrics,dyns,RKsources(k),isFinalRK)    !SR: see "Te_solver.f90"
                    ! Computation of matrix for solving ionic temperature parallel diff
                    call calcul_matTi(dyns(k),dynTis(k),Zones,dat%Nzones,RKdt,Metrics,dat)                    !SR: see "Ti_solver.f90"
                    ! Computation of source term for solving ionic temperature parallel diff
                    ! call calul_source_perp_cor(Cors(k),dyns(k),Zones,dat%Nzones,Metrics)                  !SR: see "Ti_solver.f90"
                    call calcul_sourceTi(dyns(k),dynTis(k),Zones,dat%Nzones,RKdt,dat,Cors(k),RKsources(k),isFinalRK)                 !SR: see "Ti_solver.f90"
                 end do
                 !$OMP END DO
              end if
              !$OMP BARRIER

              ! Solving Te and Ti
              if(dat%isoT.eq.0) then
                 !$OMP DO
                 do k=1,config_para%zone_para_number
                    ! --- Te ----
                    ! Concatenating zones matrix for solving parallel megazones
                    call concat_Parae(dynTes,config_para,dynsParae(k),Zones,dat%Nzones,Metrics)            !SR: see "Te_solver.f90"
                    ! Inverting the tridiagonal matrix fo Te
                    call invert_superT(dynsParae(k))                                                       !SR: see "implicit.f90" 
                    ! Un concatenating to attribute Te to each Zone
                    call unconcat_Parae(dynsParae(k),config_para,dyns,Zones,dat%Nzones,flag_erre,Metrics)  !SR: see "Te_solver.f90"
                    ! --- Ti ----
                    ! Concatenating zones matrix for solving parallel megazones
                    call concat_Parai(dynTis,config_para,dynsParai(k),Zones,dat%Nzones)                    !SR: see "Ti_solver.f90"
                    ! Inverting the tridiagonal matrix fo Ti
                    call invert_superT(dynsParai(k))                                                       !SR: see "implicit.f90"
                    ! Un concatenating to attribute Ti to each Zone
                    call unconcat_Parai(dynsParai(k),config_para,dyns,Zones,dat%Nzones)                    !SR: see "Ti_solver.f90"
                 end do
                 !$OMP END DO
              end if
              !$OMP BARRIER

              ! Solving center mesh (averaging around)
              !$OMP MASTER
              !  call solve_center(centre_,Zones,dat%Nzones,Metrics,dt)
              !$OMP END MASTER



              !$OMP DO
              do k=1,dat%NZones
                 ! Updating the field that will be used at next RK step to compute fluxes
                 call postRK(Zones(k),dyns(k),RKstep)  !SR: see "RungeKutta.f90"
              end do
              !$OMP END DO
              !$OMP BARRIER
           end do ! end of Runge-Kutta loop for RK source computing
        else
           !$OMP BARRIER
           !$OMP DO
           do k=1,dat%Nzones
              ! Initializing fields and flux before intermediate time step iteration
              ! At RKstep=4, isFinalRK flag is switched to 1
              call preRK(Zones(k),RKsources(k),RKdt,dt,isFinalRK,1)   !SR: see "RungeKutta.f90"
           end do
           !$OMP END DO

           !$OMP BARRIER
           !$OMP DO
           do k=1,dat%Nzones
              ! Broadcast of values on domain boundary (for fluxes computation - plasma(2))
              call load_neighboring_values(Zones(k),Zones,dat%Nzones,dat,centre_,Metrics)
              ! Compute some derivatives
              call calcul_deriv(dyns(k),Zones,dat%Nzones)  !(is it used ####CAREFUL####)   !SR: see "derivatives.f90"
           end do
           !$OMP END DO
           !$OMP BARRIER

           !$OMP MASTER
           ! Broadcast of the derivatives
           call load_NeighDerAll(dyns,Zones,dat%Nzones) !(is it used ####CAREFUL####)    !SR: see "derivatives.f90"

           ! Computation of source due to neutral ####CAREFUL####
           if(dat%recyc.ne.0.d0) then
              ! calcul de la source de neutres
              call calcul_Sn(Zones,dat%Nzones,dat,Metrics)    !SR: see "neutrals.f90"
           end if
           !$OMP END MASTER
           !$OMP BARRIER

           ! ###RK4 Computing of the fluxes from intermediate time steps stored in plasma(2) ###
           ! The perpendicular explicit sources are calculated (for all equations)
           !$OMP DO 
           do k=1,dat%Nzones 
              call calcul_SourcePerp(dyns(k),Zones,dat%Nzones,Metrics)  !SR: see "sources_perp.f90"
           end do
           !$OMP END DO
           !OMP BARRIER

           ! The parallel explicit sources are calculated (for all equations)
           !$OMP DO
           do k=1,dat%NZones
              call calcul_Sourcepara(dyns(k),Zones,dat%NZones,Metrics,dat%isSlab,dat)  !SR: see "sources_para.f90"
           end do
           !$OMP END DO
           !$OMP BARRIER

           ! Matrix and source vector computation for density and velocity equations
           !$OMP DO
           do k=1,dat%NZones
              call calcul_matN(dyns(k),dynNs(k),Zones,dat%Nzones,RKdt,Metrics,dat%eta,dat)  !SR: see "N_solver.f90"
              call calcul_sourceN(dyns(k),dynNs(k),Zones,dat%Nzones,RKdt,dat,dat%eta,RKsources(k),0)           !SR: see "N_solver.f90"
              call calcul_matG(dyns(k),dynGs(k),Zones,dat%Nzones,RKdt,Metrics,dat%eta)  !SR: see "G_solver.f90"
              call calcul_sourceG(dyns(k),dynGs(k),Zones,dat%Nzones,RKdt,dat,RKsources(k),0)           !SR: see "G_solver.f90"
           end do
           !$OMP END DO
           !$OMP BARRIER

           ! Solving density and velocity
           !$OMP DO
           do k=1,config_para%zone_para_number
              ! --- Density ----
              ! Concatenating zones matrix for solving parallel topological structures
              call concat_N(dynNs,config_para,SdynNs(k),Zones,dat%Nzones,Metrics)                !SR: see "N_solver.f90"
              ! Inverting the tridiagonal matrix fo Te
              call invert_superT(SdynNs(k))                                                      !SR: see "implicit.f90"
              ! Un concatenating to attribute Te to each Zone
              call unconcat_N(SdynNs(k),config_para,dyns,Zones,dat%Nzones,flag_erre,Metrics)     !SR: see "N_solver.f90"
              ! --- Velocity ---
              ! Concatenating zones matrix for solving parallel topological structures
              call concat_G(dynGs,config_para,SdynGs(k),Zones,dat%Nzones,Metrics)                !SR: see "G_solver.f90"
              ! Inverting the tridiagonal matrix fo Ti
              call invert_superT(SdynGs(k))                                                      !SR: see "implicit.f90"
              ! Un concatenating to attribute Ti to each Zone
              call unconcat_G(SdynGs(k),config_para,dyns,Zones,dat%Nzones,flag_erre,Metrics)     !SR: see "G_solver.f90"
           end do
           !$OMP END DO
           !$OMP BARRIER

           !$OMP DO 
           do k=1,dat%Nzones 
              ! Computation of some convective term in the poloidal plane from implicitely solved density
              ! This routine modifies the dyn%sources_perp tables 
              call calc_source_perp_T_imp(dyns(k),Zones,dat%NZones,Metrics)                      !SR: see "sources_perp.f90"
           end do
           !$OMP END DO
           !OMP BARRIER

           ! Matrix and source vector computation for Te and Ti equations
           ! (this is optional: if flag isoT=1, the plasma is isothermal and Te=Te0, Ti=Ti0)
           if(dat%isoT.eq.0) then
              !$OMP DO
              do k=1,dat%Nzones
                 ! Computation of matrix for solving electronic temperature parallel diff
                 call calcul_matTe(dyns(k),dynTes(k),Zones,dat%Nzones,RKdt,Metrics,dat)                    !SR: see "Te_solver.f90"
                 ! Computation of source term for solving electronic temperature parallel diff
                 call calcul_sourceTe(dyns(k),dynTes(k),Zones,dat%Nzones,RKdt,dat,Cors(k),Metrics,dyns,RKsources(k),0)    !SR: see "Te_solver.f90"
                 ! Computation of matrix for solving ionic temperature parallel diff
                 call calcul_matTi(dyns(k),dynTis(k),Zones,dat%Nzones,RKdt,Metrics,dat)                    !SR: see "Ti_solver.f90"
                 ! Computation of source term for solving ionic temperature parallel diff
                 ! call calul_source_perp_cor(Cors(k),dyns(k),Zones,dat%Nzones,Metrics)                  !SR: see "Ti_solver.f90"
                 call calcul_sourceTi(dyns(k),dynTis(k),Zones,dat%Nzones,RKdt,dat,Cors(k),RKsources(k),0)                 !SR: see "Ti_solver.f90"
              end do
              !$OMP END DO
           end if
           !$OMP BARRIER

           ! Solving Te and Ti
           if(dat%isoT.eq.0) then
              !$OMP DO
              do k=1,config_para%zone_para_number
                 ! --- Te ----
                 ! Concatenating zones matrix for solving parallel megazones
                 call concat_Parae(dynTes,config_para,dynsParae(k),Zones,dat%Nzones,Metrics)            !SR: see "Te_solver.f90"
                 ! Inverting the tridiagonal matrix fo Te
                 call invert_superT(dynsParae(k))                                                       !SR: see "implicit.f90" 
                 ! Un concatenating to attribute Te to each Zone
                 call unconcat_Parae(dynsParae(k),config_para,dyns,Zones,dat%Nzones,flag_erre,Metrics)  !SR: see "Te_solver.f90"
                 if(flag_erre.eq.1) then
                    write(40,*) ntps,': Te boundary has been changed to avoid negative temperature in zone ', k
                 end if
                 ! --- Ti ----
                 ! Concatenating zones matrix for solving parallel megazones
                 call concat_Parai(dynTis,config_para,dynsParai(k),Zones,dat%Nzones)                    !SR: see "Ti_solver.f90"
                 ! Inverting the tridiagonal matrix fo Ti
                 call invert_superT(dynsParai(k))                                                       !SR: see "implicit.f90"
                 ! Un concatenating to attribute Ti to each Zone
                 call unconcat_Parai(dynsParai(k),config_para,dyns,Zones,dat%Nzones)                    !SR: see "Ti_solver.f90"
              end do
              !$OMP END DO
           end if
           !$OMP BARRIER

           ! Solving center mesh (averaging around)
           !$OMP MASTER
           !  call solve_center(centre_,Zones,dat%Nzones,Metrics,dt)
           !$OMP END MASTER

           !$OMP DO
           do k=1,dat%NZones
              ! Updating the field that will be used at next RK step to compute fluxes
              call postRK(Zones(k),dyns(k),1)  !SR: see "RungeKutta.f90"
           end do
           !$OMP END DO
           !$OMP BARRIER

        end if

        !$OMP DO
        do k=1,dat%Nzones
           ! Broadcast of values on domain boundary (for fluxes computation - plasma(2))
           call load_neighboring_values(Zones(k),Zones,dat%Nzones,dat,centre_,Metrics)
           ! Compute some derivatives
           call calcul_deriv(dyns(k),Zones,dat%Nzones)  !(is it used ####CAREFUL####)   !SR: see "derivatives.f90"
        end do
        !$OMP END DO
        !$OMP BARRIER

        !$OMP MASTER
        ! Broadcast of the derivatives
        call load_NeighDerAll(dyns,Zones,dat%Nzones) !(is it used ####CAREFUL####)    !SR: see "derivatives.f90"
        !$OMP END MASTER
        !$OMP BARRIER

        ! Computation of source due to neutral
        !$OMP DO
        do k=1,dat%Nzones
           if(dat%isNeutralON.eq.1) then
              ! calcul de la source de neutres
              call compute_outflux(Zones(k),Metrics(k),Flux_ns(k),dat,Source(k))
           end if
        end do
        !$OMP END DO
        !$OMP BARRIER

        !$OMP DO
        do k=1,dat%Nzones
           if(dat%isNeutralON.eq.1) then
              call evolN(Zones(k),dt,Metrics(k),Flux_ns(k),dyns(k),dat)
           end if
        end do
        !$OMP END DO
        !$OMP BARRIER
        !$OMP MASTER
        Stot=0.d0
        do k=1,dat%Nzones
           Stot=Stot+Source(k)
        end do
        !$OMP END MASTER
        !$OMP BARRIER

        ! cleaning and looking for errors (NaN, n<0, Te<0, Ti<0)
        !$OMP DO
        do k=1,dat%Nzones
           call clean_solution(dyns(k),Zones,dat%Nzones)            !SR: see "monitors.f90"
           call check_solution(dyns(k),Zones,dat%Nzones)            !SR: see "monitors.f90"
        end do
        !$OMP END DO
        !$OMP BARRIER

        ! deciding whether solution is satisfactory or not
        ! if not: SolOK=0 and back to beginning with smaller CFL
        ! if solution is OK, fields are updated and go to iteration n+1
        !$OMP MASTER
        normal_zone_number=0
        do k=1,dat%Nzones
           normal_zone_number=normal_zone_number+dyns(k)%isOK
        end do
        if(normal_zone_number.eq.dat%Nzones) then
           SolOK=1 ! Solution OK let us leave the loop

        else
           CFL=CFL*0.5d0 ! if Sol NOK, let us reduce CFL
           write(*,*) 'CFL reduced to ',CFL
           if(CFL.lt.0.01) then
              write(*,*) 'Soledge2D crashed at iteration ', ntps
              stop
           end if
        end if
        !$OMP END MASTER
        !$OMP BARRIER
     end do

     !$OMP DO
     do k=1,dat%Nzones
        ! The plasma parameters are updated
        call transfer_dyn(dyns(k),Zones,dat%Nzones,dat%isoT)      
     end do
     !$OMP END DO
     !$OMP BARRIER

     !$OMP MASTER
     tempus = tempus+dt

     ! The residues are calculated (based on max norm)
     call calcul_res(dyns,dat%Nzones,resN,resG,resTe,resTi,resNn,dt)    !SR: see "monitors.f90"
     if(dat%nite>=1000) then
        if(modulo(ntps,dat%nite/1000).eq.0) then    
           if(flag_erre.eq.1) then
              write(*,*) 'Negative Te on BC (mesh too scarce to handle Te gradients) - might cause solution to oscilate'
           end if
           write(*,211) real(ntps)/real(dat%nite)*100,dt ,limdir, index1
           call write_temporal_custom(Zones,dat%NZones,dat) ! writing Mach for monitoring with gnuplot
           write(20,600) resN, resG, resTe, resTi, resNn
           write(80,*) Stot
        endif
     else
        if(flag_erre.eq.1) then
           write(*,*) 'Negative Te on BC (mesh too scarce to handle Te gradients) - might cause solution to oscilate'
        end if
        write(*,211) real(ntps)/real(dat%nite)*100, dt, limdir, index1
        call write_temporal_custom(Zones,dat%NZones,dat) ! writing Mach for monitoring with gnuplot
        write(20,600) resN, resG, resTe, resTi, resNn
        write(80,*) Stot
     end if

     !Saving fields for video
     if(dat%fframe.gt.0) then
        if(dat%nite.gt.dat%fframe) then
           if(modulo(ntps,dat%fframe).eq.0) then
              !              call write_champs3(Centre_,Zones,dat%Nzones,ntps/dat%fframe)    !SR: see "saving.f90"
              write(168,102) (Zones(2)%plasma(1)%Mach(5,j),j=1,180)
              write(167,102) (Zones(2)%plasma(1)%density(5,j),j=1,180)
              write(169,102) (Zones(2)%plasma(1)%Te(5,j),j=1,180)
              write(170,102) (Zones(2)%plasma(1)%Ti(5,j),j=1,180)
102           format(512es15.7)
           end if
        end if
     end if

     ! DJM 27/1/12 save time in time array
     time_array(ntps)=tempus
     ! BEGIN DJM 19/1/12 save plasma parameters at required positions:
     do k=1,dat%Nzones
        do j=1,Zones(k)%ntracks
           Zones(k)%plasma(1)%neTime(ntps,j)=&
                Zones(k)%plasma(1)%density(Zones(k)%xtrackindex(j),Zones(k)%ztrackindex(j))
           Zones(k)%plasma(1)%MachTime(ntps,j)=&
                Zones(k)%plasma(1)%Mach(Zones(k)%xtrackindex(j),Zones(k)%ztrackindex(j))
           Zones(k)%plasma(1)%TeTime(ntps,j)=&
                Zones(k)%plasma(1)%Te(Zones(k)%xtrackindex(j),Zones(k)%ztrackindex(j))
           Zones(k)%plasma(1)%TiTime(ntps,j)=&
                Zones(k)%plasma(1)%Ti(Zones(k)%xtrackindex(j),Zones(k)%ztrackindex(j))
        end do
     end do

     !$OMP END MASTER
     !$OMP BARRIER
  end do
  !$OMP END PARALLEL
  close(20)
  close(40)
  close(70)

  ! Broadcast of boundary values before saving
  do k=1,dat%Nzones
     call load_neighboring_values(Zones(k),Zones,dat%Nzones,dat,centre_,Metrics)
  end do
  ! Saving of the results
  call write_champs4(Zones,dat%Nzones,dat,centre_,time_array)    !SR: see "saving.f90"
  call write_neutrals(Zones,Flux_ns,dat%NZones,centre_)   !SR: see "saving.f90"
  call save_neighbors(Zones,dat%Nzones)           !SR: see "saving.f90"

  call sol2eir(InterpData2_,dat%NZones,Zones,MagFields)
  call write_res_knots(InterpData2_)

  write(*,213) tempus

700 format( 512es15.7 )
211 format (F5.1,1X,'%',3X,es15.7,3X,I2,3X,I2)
212 format ('Avancement : 000.0 %')
213 format ('Temps final: ',F10.4)
600 format ( 5es15.7 )
605 format (es15.7,3x,es15.7,3x,es15.7,3x,I2,3X,I2)
end program soledge2m

