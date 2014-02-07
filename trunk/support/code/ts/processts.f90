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

  real :: psi_min,psi_max
  integer :: npsi_bins
  integer,external :: iargc


  logical :: elm_filt,remove_outlier

  !
  elm_filt = .false.
  remove_outlier = .true.

  !
  ! Read and assign command line arguments
  !
  ! Command line is:  'file name'   tmin    tmax   chi_lim -e 'elm file name'
  !
  nargs = iargc()

  if (nargs.lt.3) then 
     write(0,'(a)') 'procts command usage: procts <infilename> tmin  tmax  -e  <elmtimefilename>'
     write(0,'(a)') 'The -e argument is optional'
     STOP 'Insufficient command line arguments'
  endif

  call getarg(1,arg)
  infilename = trim(arg)

  call getarg(2,arg)
  read(arg,*) tmin
  write(0,*) 'TMIN=',tmin

  call getarg(3,arg)
  read(arg,*) tmax
  write(0,*) 'TMAX=',tmax

  call getarg(4,arg)
  if (arg.eq.'-e') then 
     elm_filt = .true.
     call getarg(5,arg)
     elm_filename = trim(arg)
     write(0,*) 'ELM=',trim(elm_filename)
  endif

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
