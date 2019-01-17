      ! -*-Fortran-*-
      !
      ! The purpose of this module is to make the Fortran 77 fixed format
      ! parameter declarations available and usable in any future fortran 90
      ! source code without having to use fixed form source code. 
      !
      ! Depending on how well this works (or not) - common blocks may be
      ! converted to modules in the future. 
      !
      ! jde - June 20th, 2005
      !
      module mod_params

         implicit none
         include 'params'

      end module mod_params
