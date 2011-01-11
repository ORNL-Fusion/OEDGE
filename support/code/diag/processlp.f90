program processfs
  use analyse_lp

  implicit none

  character*100 :: filename,filenameout
  integer:: nlines,in
  real :: tstart, tend

  integer ::  shot 

  real :: minval,maxval

  integer :: nbins,printopt
  integer :: ierr
  
  logical :: first_run

  !
  ! Print opt = 0 ... gives baseline average evolution over time
  !
  ! Print opt = 1 ... gives all print data
  !
  !

  printopt = 0



  nbins  = 200



  do shot = 134582,134597

     write(filename,'(i6,a)') shot,'_tab.dat'
     write(0,'(a,a)') 'Filename:'//trim(filename)

     call read_lp(filename,ierr)

     if (ierr.ne.0) cycle

     call setup_lp(nbins)

     !write(0,'(a,a)') '3 Filename:'//trim(filename),' Channel: '//trim(fschan)

     call accumulate_data

     write(filenameout,'(i6,a)') shot,'_proc.dat'

     !write(0,'(a,a)') '4 Filename:'//trim(filename),' Channel: '//trim(fschan)

     call analyse_print_data(filenameout,printopt)

     write(0,'(a,a)') 'FilenameOut:'//trim(filenameout)

  end do



  first_run = .true.

  do shot = 134582,134597

     write(filename,'(i6,a)') shot,'_tab.dat'
     write(0,'(a,a)') 'Filename:'//trim(filename)

     call read_lp(filename,ierr)

     if (ierr.ne.0) cycle

     if (first_run) then 
        call setup_lp(nbins)
        first_run = .false.
     endif


     !write(0,'(a,a)') '3 Filename:'//trim(filename),' Channel: '//trim(fschan)

     call accumulate_data


  end do

     write(filenameout,'(a)') 'lp_proc_summary.dat'

     !write(0,'(a,a)') '4 Filename:'//trim(filename),' Channel: '//trim(fschan)

     call analyse_print_data(filenameout,printopt)

     write(0,'(a,a)') 'FilenameOut:'//trim(filenameout)





end program processfs
