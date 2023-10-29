!
! jdemod - changed mod_divimp.f to mod_divimp.f90
!
! ======================================================================
! ======================================================================
!
MODULE mod_divimp_walldyn
  IMPLICIT none
  SAVE
  PUBLIC
  !
  !     ------------------------------------------------------------------
  ! 
  TYPE, PUBLIC :: type_walldyn

     REAL*4              :: n       ! neutral redistribution
     REAL*4              :: i       ! ion redistribution
     REAL*4, ALLOCATABLE :: iz (:)  ! ion redistribution by charge state
     REAL*4, ALLOCATABLE :: eiz(:)  ! energy of redistributed ions by charge state

  ENDTYPE type_walldyn

  INTEGER :: wdn_i
  TYPE(type_walldyn), ALLOCATABLE :: wdn(:,:)
  !
  !     ==================================================================
  !
CONTAINS
  !
  !     ------------------------------------------------------------------
  ! 
  SUBROUTINE divWdnAllocate(wallpts,nizs)
    IMPLICIT none

    INTEGER, INTENT(IN) :: wallpts,nizs

    INTEGER :: i,j

    IF (ALLOCATED(wdn)) THEN
       WRITE(0,*) 'ERROR divWdnAllocate: WDN array already allocated'
       STOP
    ENDIF

    wdn_i = 1

    ALLOCATE(wdn(wallpts+1,wallpts+1))

    DO i = 1, wallpts+1
       DO j = 1, wallpts+1
          ALLOCATE(wdn(i,j)%iz (nizs+1))
          ALLOCATE(wdn(i,j)%eiz(nizs+1))
          wdn(i,j)%n   = 0.0
          wdn(i,j)%i   = 0.0
          wdn(i,j)%iz  = 0.0
          wdn(i,j)%eiz = 0.0
       END DO
    END DO
    RETURN

  END SUBROUTINE divWdnAllocate
  !
  !     ==================================================================
  !
END MODULE mod_divimp_walldyn
!
! ======================================================================
! ======================================================================
!
MODULE mod_divimp
  IMPLICIT none
  SAVE
  PUBLIC
  !
  !     Parameters:
  !     ------------------------------------------------------------------
  INTEGER, PARAMETER :: MAXNLAUNCH = 1000,&  
                        MAXNBLK    =    1,&
                        MAXNATM    =    2,&
                        MAXNMOL    =    1,&
                        MAXNION    =    1,&
                        MAXNPHO    =    1,&
                        MAXNSRC    =    5  ! bulk, atoms, molecules, test ions, photons
  !
  !     Impurity-wall interaction structure:
  !     ------------------------------------------------------------------
  TYPE, PUBLIC :: type_wall_flux

     ! jdemod - added initialization of type elements so that this hopefully carries through when wall_flx is allocated(?)
     REAL*4 :: length
     REAL*4 :: area 

     REAL*4 :: in_par_blk(MAXNBLK,0:MAXNSRC)  ! *** NOTE *** if anything is changed here then the code in 
     REAL*4 :: in_ene_blk(MAXNBLK,0:MAXNSRC)  ! DIVSTORE and IOOUT needs to be update (slver=3.6)
     REAL*4 :: in_par_atm(MAXNATM,0:MAXNSRC)
     REAL*4 :: in_ene_atm(MAXNATM,0:MAXNSRC)
     REAL*4 :: in_par_mol(MAXNMOL,0:MAXNSRC)
     REAL*4 :: in_ene_mol(MAXNMOL,0:MAXNSRC)

     REAL*4 :: em_par_atm(MAXNATM,0:MAXNSRC) 
     REAL*4 :: em_ene_atm(MAXNATM,0:MAXNSRC)
     REAL*4 :: em_par_mol(MAXNMOL,0:MAXNSRC)

     REAL*4 :: launch(MAXNLAUNCH)
     REAL*4 :: prompt

  ENDTYPE type_wall_flux
  !
  !     -------------------------------------------------------------------     
  INTEGER :: wall_n , wall_nlaunch    

  TYPE(type_wall_flux), ALLOCATABLE :: wall_flx(:)

  !     DIVIMP input options:
  !     ------------------------------------------------------------------
  TYPE, PUBLIC :: type_options_divimp
     REAL*4    :: version = 1.0
     !...     Ribbon grid:
     INTEGER       :: rib_n
     INTEGER       :: rib_type
     INTEGER       :: rib_format
     CHARACTER*512 :: rib_file            ! file containing ribbon data
     INTEGER       :: rib_region          ! region of interest
     REAL          :: rib_r1              ! region of interest
     REAL          :: rib_z1            
     REAL          :: rib_r2  
     REAL          :: rib_z2  
     INTEGER       :: rib_rad_opt         ! radial refinement
     REAL          :: rib_rad_a      
     REAL          :: rib_rad_b
     REAL          :: rib_rad_c
     REAL          :: rib_rad_d
     INTEGER       :: rib_pol_opt         ! poloidal refinement
     INTEGER       :: rib_pol_n       
     REAL          :: rib_pol_a        
     REAL          :: rib_pol_b
     REAL          :: rib_pol_c
     REAL          :: rib_pol_d
     INTEGER       :: rib_pol_n_def       ! poloidal refinement
     REAL          :: rib_pol_a_def    
     REAL          :: rib_pol_b_def
     REAL          :: rib_pol_c_def
     REAL          :: rib_pol_d_def
     !...     Sputter data from other runs:
     !         INTEGER       :: nsputter
     !...     Store particle state data:
     INTEGER       :: pstate
     !...     Store particle state data:
     INTEGER       :: niter               ! number of DIVIMP iterations
  ENDTYPE type_options_divimp

  TYPE(type_options_divimp) :: opt_div
  INTEGER                   :: div_iter
  !
  !     ==================================================================
  !
  TYPE :: type_sputter_data
     REAL                :: version
     INTEGER             :: format
     INTEGER             :: data_type
     CHARACTER*256       :: case_name
     CHARACTER*256       :: extension
     REAL                :: fraction
     CHARACTER*256       :: tag
     REAL                :: absfac
     REAL                :: absfac_net
     INTEGER             :: atomic_number
     REAL                :: atomic_mass
     INTEGER             :: target_number
     INTEGER             :: charge
     INTEGER             :: charge_min
     INTEGER             :: charge_max
     INTEGER             :: nsegments
     INTEGER             :: type
     INTEGER,ALLOCATABLE :: ik(:)
     INTEGER,ALLOCATABLE :: ir(:)
     INTEGER,ALLOCATABLE :: id(:)
     REAL,ALLOCATABLE    :: r    (:)
     REAL,ALLOCATABLE    :: z    (:)
     REAL,ALLOCATABLE    :: dds  (:)
     REAL,ALLOCATABLE    :: te   (:)
     REAL,ALLOCATABLE    :: ti   (:)
     REAL,ALLOCATABLE    :: flux (:,:)
     REAL,ALLOCATABLE    :: e0   (:,:)
     REAL,ALLOCATABLE    :: yield(:,:)
     REAL,ALLOCATABLE    :: modifier(:)
     INTEGER             :: ncore
     REAL,ALLOCATABLE    :: core_rho (:)
     REAL,ALLOCATABLE    :: core_psin(:)
     REAL,ALLOCATABLE    :: core_ne  (:)
     REAL,ALLOCATABLE    :: core_te  (:)
     REAL,ALLOCATABLE    :: core_ti  (:)
     REAL,ALLOCATABLE    :: core_percent_nfrac(:)
     REAL,ALLOCATABLE    :: core_percent_efrac(:)
  ENDTYPE type_sputter_data

  TYPE(type_sputter_data), ALLOCATABLE :: sputter_data(:)
  INTEGER :: sputter_ndata,sputter_ilast
  REAL    :: sputter_nabsfac              ! temporary storage of nabsfac from divCompileSputteringYields
  !
  !     ------------------------------------------------------------------
  ! 
  INTEGER :: nymfs_global
  !
  !     ==================================================================
  !
CONTAINS
  !
  !     ------------------------------------------------------------------
  ! 
  SUBROUTINE divClean
    USE mod_divimp_walldyn

    IF (ALLOCATED(wall_flx)) DEALLOCATE(wall_flx)
    IF (ALLOCATED(wdn     )) DEALLOCATE(wdn)

  END SUBROUTINE divClean
  !
  !     ------------------------------------------------------------------
  ! 
  SUBROUTINE divInitializeOptions

    !      opt_div%nsputter = 0
    sputter_ndata = 0
    sputter_ilast = 0
    sputter_nabsfac = 0.0
    nymfs_global = 0

  END SUBROUTINE divInitializeOptions
  !
  !     ------------------------------------------------------------------
  ! 
  subroutine init_wall_flx(nwall)
    implicit none
    ! jdemod - this routine intiializes the wall_flx arrays to 0.0
    ! since this is not done by allocation itself
    integer :: in,nwall


    if (allocated(wall_flx)) then 
       do in = 1,nwall

          !     jdemod - added initialization of wall_flx
          wall_flx(in)%length  = 0.0
          wall_flx(in)%area    = 0.0

          wall_flx(in)%in_par_blk = 0.0 ! *** NOTE *** if anything is changed here then the code in 
          wall_flx(in)%in_ene_blk = 0.0 ! DIVSTORE and IOOUT needs to be update (slver=3.6)
          wall_flx(in)%in_par_atm = 0.0
          wall_flx(in)%in_ene_atm = 0.0
          wall_flx(in)%in_par_mol = 0.0
          wall_flx(in)%in_ene_mol = 0.0

          wall_flx(in)%em_par_atm = 0.0
          wall_flx(in)%em_ene_atm = 0.0
          wall_flx(in)%em_par_mol = 0.0

          wall_flx(in)%launch = 0.0
          wall_flx(in)%prompt = 0.0
       end do
    endif
  end subroutine init_wall_flx


END MODULE mod_divimp
!
! ======================================================================
! ======================================================================
!
MODULE mod_divimp_cneut
  IMPLICIT none
  SAVE
  PUBLIC
  !
  !     ------------------------------------------------------------------
  ! 
  INTEGER :: cneut_test
  !
  !     ==================================================================
  !
  !      CONTAINS

  !      SUBROUTINE     divClean
  !        IF (ALLOCATED(wall_flx)) DEALLOCATE(wall_flx)
  !      END SUBROUTINE divClean

END MODULE mod_divimp_cneut
!
! ======================================================================
! ======================================================================
!
MODULE mod_divimp_tdep
  IMPLICIT none
  SAVE
  PUBLIC
  !
  !     ------------------------------------------------------------------
  ! 
  TYPE, PUBLIC :: type_tdep_data
     REAL :: r
     REAL :: z
     REAL :: phi
     REAL :: s
     REAL :: cross
     REAL :: diag
     REAL :: temp
     REAL :: vel
     REAL :: charge
     REAL :: weight
  ENDTYPE type_tdep_data
  !
  !     ------------------------------------------------------------------
  ! 
  INTEGER :: dummy

  LOGICAL :: tdep_data_exists

  TYPE(type_tdep_data), ALLOCATABLE :: tdep_save(:)
  INTEGER :: tdep_save_n

  TYPE(type_tdep_data), ALLOCATABLE :: tdep_load(:)
  REAL    :: tdep_load_version
  INTEGER :: tdep_load_n
  REAL    :: tdep_load_absfac
  REAL    :: tdep_load_deltat
  REAL    :: tdep_load_qtim
  REAL    :: tdep_load_frac
  REAL    :: tdep_load_ions_injected
  REAL    :: tdep_load_ions_to_target

END MODULE mod_divimp_tdep
!
! ======================================================================
! ======================================================================
!
MODULE mod_divimp_3d
  IMPLICIT none
  SAVE
  PUBLIC
  !
  !     ------------------------------------------------------------------
  ! 
  INTEGER :: cneut_test

  !
  !     ==================================================================
  !
  !      CONTAINS
  !
  !      SUBROUTINE     divClean
  !        IF (ALLOCATED(wall_flx)) DEALLOCATE(wall_flx)
  !      END SUBROUTINE divClean
  !
END MODULE mod_divimp_3d
!
! ======================================================================
! ======================================================================
!
