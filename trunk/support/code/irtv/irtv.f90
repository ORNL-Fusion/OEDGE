program irtv
  use process_irtv
  implicit none


  character*256 :: filename
  integer :: ierr,opt 
  logical :: finished


  finished = .false.


  do while (.not.finished)


     write (6,*) 'Options: 1-process  2-exit'
     read(5,*) opt

     if (opt.eq.1) then 

        write(6,*) 'Enter file name:'
        read(5,*) filename

        call open_irtv(filename,ierr)

        if (ierr.eq.0) then 

           call process_irtv_commands
        endif

     elseif (opt.eq.2) then 

        finished=.true.

     endif


  end do






end program irtv
