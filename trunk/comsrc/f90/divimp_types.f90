module divimp_types
  INTEGER, PARAMETER :: R8 = SELECTED_REAL_KIND (14)
  INTEGER, PARAMETER :: R4 = SELECTED_REAL_KIND (6)


!
! jdemod - type name changed to cell_type from type_cell due to potential for name
!          conflict of type type_cell with the same name in mod_eirene04.f90
!
      TYPE cell_type
        INTEGER :: index,ik,ir,nv,rzone,zzone,xpt,map
        integer :: mapik,mapir
        !integer :: targ,wall
        integer :: connect(4)
        REAL    :: rcen,zcen,bratio,rv(4),zv(4)
      ENDTYPE cell_type


end module divimp_types
