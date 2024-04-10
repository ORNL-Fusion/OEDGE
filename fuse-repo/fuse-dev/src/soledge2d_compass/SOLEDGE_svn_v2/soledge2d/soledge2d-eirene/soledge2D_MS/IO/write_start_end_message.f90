  subroutine write_start_message()
    write(*,*) '   _____       ________    __          '
    write(*,*) '  / ___/____  / / ____/___/ /___ ____  '
    write(*,*) '  \__ \/ __ \/ / __/ / __  / __ `/ _ \ '
    write(*,*) ' ___/ / /_/ / / /___/ /_/ / /_/ /  __/ '
    write(*,*) '/____/\____/_/_____/\__,_/\__, /\___/  '
    write(*,*) '                         /____/        '
    write(*,*) 
  end subroutine write_start_message

  subroutine write_loop_begin_message()
    write(*,*)
    write(*,*) ' SolEdge2D simulation started:'
  end subroutine write_loop_begin_message

  subroutine write_end_message()
    write(*,*) 
    write(*,*) 'Soledge2D simulation finished'
  end subroutine write_end_message
