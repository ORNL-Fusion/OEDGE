!     -*-Fortran-*-
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
      INTEGER, PARAMETER :: MAXNLAUNCH = 1000,  
     .                      MAXNBLK    =    1,
     .                      MAXNATM    =    2,
     .                      MAXNMOL    =    1,
     .                      MAXNION    =    1,
     .                      MAXNPHO    =    1,
     .                      MAXNSRC    =    5  ! bulk, atoms, molecules, test ions, photons
!
!     Impurity-wall interaction structure:
!     ------------------------------------------------------------------
      TYPE, PUBLIC :: type_wall_flux

         REAL*4 :: length
         REAL*4 :: area

         REAL*4 :: in_par_blk(MAXNBLK,0:MAXNSRC)
         REAL*4 :: in_ene_blk(MAXNBLK,0:MAXNSRC)
         REAL*4 :: in_par_atm(MAXNATM,0:MAXNSRC)
         REAL*4 :: in_ene_atm(MAXNATM,0:MAXNSRC)
         REAL*4 :: in_par_mol(MAXNMOL,0:MAXNSRC)
         REAL*4 :: in_ene_mol(MAXNMOL,0:MAXNSRC)

         REAL*4 :: em_par_atm(MAXNATM,0:MAXNSRC)
         REAL*4 :: em_ene_atm(MAXNATM,0:MAXNSRC)

         REAL*4 :: launch(MAXNLAUNCH)
         REAL*4 :: prompt

      ENDTYPE type_wall_flux
!
!     -------------------------------------------------------------------     
      INTEGER :: wall_n      ,  
     .           wall_nlaunch    

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
      ENDTYPE type_options_divimp
 
      TYPE(type_options_divimp) :: opt_div
!
!     ==================================================================
!
      CONTAINS

      SUBROUTINE     divClean
        IF (ALLOCATED(wall_flx)) DEALLOCATE(wall_flx)
      END SUBROUTINE divClean

      END MODULE mod_divimp
!
! ======================================================================
! ======================================================================
!
