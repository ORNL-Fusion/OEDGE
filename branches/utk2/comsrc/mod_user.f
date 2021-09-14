!     -*-Fortran-*-
      MODULE mod_user
      IMPLICIT none
      PUBLIC


      TYPE :: type_options_user
         INTEGER :: u_mom  ! User defined momentum volume source
      ENDTYPE type_options_user


      TYPE(type_options_user) :: opt_user


      REAL, ALLOCATABLE :: stuff(:)  ! <--- James' momentum source!



      END MODULE mod_user
