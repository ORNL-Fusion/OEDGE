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
      INTEGER, PARAMETER :: MAXNLAUNCH = 10,  
     .                      MAXNBLK    =  1,
     .                      MAXNATM    =  2,
     .                      MAXNMOL    =  1,
     .                      MAXNION    =  1,
     .                      MAXNPHO    =  1,
     .                      MAXNSRC    =  5  ! bulk, atoms, molecules, test ions, photons
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
