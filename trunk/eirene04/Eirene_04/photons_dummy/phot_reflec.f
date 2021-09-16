      module phot_reflec
      use precision
      implicit none

      private

      public :: init_refl_hlm, reflect_hollmann

      contains
 
      subroutine init_refl_hlm()
      integer :: i
      
      write (6,*) ' init_refl_hlm: nothing done '
      write (6,*) ' photonic reflection not yet available '
      return

      end subroutine init_refl_hlm

  
! reflect
! Die Routine entscheidet, ob ein Photon mit Einfallwinkel theta_i und
! Wellenlaenge lambda_i reflektiert wird. Im Falle einer 
! Reflektion ist flag 1, ansonsten 0. 
      subroutine reflect_hollmann(theta_i, lambda_i, mat, flag, 
     .                   theta_out, alpha_out, rprob)
      real(dp), intent(in) :: theta_i, lambda_i
      integer, intent(in) :: mat
      real(dp), intent(out) :: theta_out, alpha_out, rprob
      integer, intent(out) :: flag

      theta_out = 0._dp
      alpha_out = 0._dp
      rprob = 0._dp
      flag = 0

      write (6,*) ' reflect_hollmann: nothing done '
      write (6,*) ' photonic reflection not yet available '
      return
     
      end subroutine reflect_hollmann

      end module phot_reflec




