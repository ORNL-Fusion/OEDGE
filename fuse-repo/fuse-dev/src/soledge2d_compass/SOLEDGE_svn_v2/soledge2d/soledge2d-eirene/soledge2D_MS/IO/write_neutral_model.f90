  subroutine write_neutral_model()
    use all_variables, only : global_parameters,flags
    implicit none
    write(*,*) 
    write(*,*) 'Neutral model :'
    write(*,*) '~~~~~~~~~~~~~~'
    if(flags%neutral_model.eq.0) then
       write(*,*) ' - Type: None'
    else
       if(flags%neutral_model.eq.1) then
          write(*,*) ' - Type: EIRENE'
       else
          write(*,*) ' - Type: Fluid'
       end if
    end if
  end subroutine write_neutral_model
