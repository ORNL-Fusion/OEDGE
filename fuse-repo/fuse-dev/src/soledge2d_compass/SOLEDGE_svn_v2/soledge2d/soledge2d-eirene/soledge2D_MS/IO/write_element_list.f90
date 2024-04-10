  subroutine write_element_list()
    use all_variables, only : global_parameters
    implicit none
    integer*4 i
    character(2) :: buffer
    write(*,*)
    write(*,*) 'Chemistry :'
    write(*,*) '~~~~~~~~~~'
    if(global_parameters%N_species.gt.1) then
       write(*,1) global_parameters%N_species
    else
       write(*,2) global_parameters%N_species
    end if
1 format('The plasma is made of ', I2 ,' elements:')
2 format('The plasma is made of ', I2 ,' element:')
    do i=1,global_parameters%N_species
       write(buffer,'(I2)') i 
       write(*,*) '- '//global_parameters%element_list(i)%name!,' (',trim(adjustl(buffer)),')'
    end do
  end subroutine write_element_list

