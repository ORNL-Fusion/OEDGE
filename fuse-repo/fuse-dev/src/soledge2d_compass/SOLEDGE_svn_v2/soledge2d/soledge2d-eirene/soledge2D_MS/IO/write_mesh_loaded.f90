  subroutine write_mesh_loaded()
    use all_variables, only : global_parameters,flags
    implicit none
    character(2) :: buffer
    write(buffer,'(I2)') global_parameters%N_zones
    write(*,*) 
    write(*,*) 'Mesh :'
    write(*,*) '~~~~~'
    write(*,*) 'Mesh loaded: '//trim(adjustl(buffer))//' zones'
    if(flags%is_Pen) then
       write(*,*) ' - Penalisation: ON'
    else
       write(*,*) ' - Penalisation: OFF'
    end if
    if(flags%is_SLAB) then
       write(*,*) ' - Simplified slab geometry'
    else
       write(*,*) ' - Full metric'
    end if
  end subroutine write_mesh_loaded
