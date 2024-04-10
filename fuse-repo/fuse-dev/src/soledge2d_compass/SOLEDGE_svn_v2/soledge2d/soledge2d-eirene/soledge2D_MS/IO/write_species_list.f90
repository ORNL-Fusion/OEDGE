  subroutine write_species_list()
    use all_variables, only : global_parameters
    implicit none
    character(2) :: buffer
    character(2) :: buffer2
    integer*4 :: n
    write(buffer,'(I2)') global_parameters%N_ions
    write(*,*) 'Total number of species to simulate: 1 + '//trim(adjustl(buffer))
    write(*,*) ' 0 - electrons (1-)'
    do n=1,global_parameters%N_ions
       write(buffer,'(I2)') n
       write(buffer2,'(I2)') global_parameters%ions_list(n,2)
       write(*,*) buffer//' - '//global_parameters%element_list(global_parameters%ions_list(n,1))%symbol&
            //' ('//trim(adjustl(buffer2))//'+)'
    end do
    open(unit=156,file='ions_list',status='unknown')
    write(156,*) 'e-'
    do n=1,global_parameters%N_ions
       write(buffer2,'(I2)') global_parameters%ions_list(n,2)
       write(156,*) trim(adjustl(global_parameters%element_list(global_parameters%ions_list(n,1))%symbol))&
            //trim(adjustl(buffer2))//'+'
    end do
    close(156)
  end subroutine write_species_list
