      ! -*- Fortran -*-
      !
      ! The purpose of this code is to make the far periphery common block 
      ! variables which are written in fixed format source code available in 
      ! free format source code fortran 90 modules that need them. 
      !

      module mod_fperiph
         use mod_params
         implicit none
 
         include 'fperiph_com'

      end module mod_fperiph
