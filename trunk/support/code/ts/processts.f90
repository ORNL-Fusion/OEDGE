program processts
  use analyse_ts 

  implicit none


  integer:: in


  character*512 :: infilename,ident,time,ofilename,elm_filename

  integer :: iunit, ounit,ierr

  !real,allocatable :: lp_data(:,:) ,lp_axis(:),lp_proc_data(:,:,:),lp_axis_psi(:,:)

  integer :: nlines,npts,ndata,ncols,nextra

  real :: tmin,tmax,deltar
  !character*5 :: tmins,tmaxs
  character*4 :: iform
  character*5 :: chisq
  integer :: exp_tmin,exp_tmax,nargs
  character*512 :: arg
  character*512 :: time_filename
  character*3 :: time_ext
  logical :: time_filt 

  real :: psi_min,psi_max
  integer :: npsi_bins,arg_cnt
  integer,external :: iargc


  logical :: elm_filt,remove_outlier

  !
  time_filt = .false.
  elm_filt = .false.
  remove_outlier = .true.

  !
  ! Read and assign command line arguments
  !
  ! Command line is:  'file name'   tmin    tmax   chi_lim -e 'elm file name'
  !
  nargs = iargc()

  if (nargs.lt.3) then 
     write(0,'(a)') 'procts command usage: procts <infilename> tmin  tmax  -e  <elmtimefilename> -t <time windows filename>'
     write(0,'(a)') 'The -e argument is optional'
     STOP 'Insufficient command line arguments'
  endif

  arg_cnt = 1

  do while (arg_cnt.le.nargs) 


    !call getarg(1,arg)

     ! Arguments are either in order or preceded by a flag
     ! First 4 arguments are in order : file tmin tmax chi_limit
     call getarg(arg_cnt,arg)
     arg_cnt = arg_cnt + 1

     write(0,'(a,i6,a,a,a)') 'ARG:',arg_cnt,':',trim(arg),':'

     ! arg index increased by 1 since arg_cnt is incremented after read

     if (arg_cnt.eq.2) then 
        infilename = trim(arg)

     elseif (arg_cnt.eq.3) then 
        read(arg,*) tmin
        write(0,*) 'TMIN=',tmin

     elseif (arg_cnt.eq.4) then 
        read(arg,*) tmax
        write(0,*) 'TMAX=',tmax


     elseif (arg.eq.'-e') then 



!  call getarg(4,arg)
!  if (arg.eq.'-e') then 
!     elm_filt = .true.
!     call getarg(5,arg)
!     elm_filename = trim(arg)
!     write(0,*) 'ELM=',trim(elm_filename)
!  endif

        elm_filt = .true.
        call getarg(arg_cnt,arg)
        arg_cnt = arg_cnt + 1

        elm_filename = trim(arg)
        write(0,*) 'ELM=',trim(elm_filename)

     elseif (arg.eq.'-t') then 
        time_ext='_tf'

        time_filt = .true.
        call getarg(arg_cnt,arg)
        arg_cnt = arg_cnt + 1

        time_filename = trim(arg)
        write(0,*) 'Time windows file=',trim(time_filename)

        call setup_time_filter(time_filt,time_filename)

    endif

  end do


  psi_min = 0.2
  psi_max = 1.3
  npsi_bins = int((psi_max-psi_min)/0.005)



  call setup_ts_bins(tmin,tmax,psi_min,psi_max,npsi_bins)

  write(0,*) 'Filename:',trim(infilename),':'

  call read_ts(infilename,nlines,elm_filt,elm_filename)

  call accumulate_ts_data(nlines,elm_filt)

  call analyse_print_ts_data(infilename,elm_filt)

  call deallocate_ts_data



end program processts
