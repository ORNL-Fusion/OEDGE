program proclp
  use utilities
  !use array_sort
  use langmuir_probes

  implicit none

  ! This routine reads in a typical Lanmgmuir probe data file from Jon and will 
  ! then bin average the data within a specified time window. 
  !
  ! Three distributions are possible 
  ! A) averaged by time window
  ! B) averaged by R-Rsep bins
  ! C) Binned by Jsat magnitude and R-Rsep
  !

  character*512 :: infilename,ident,time,ofilename

  integer :: iunit, ounit,ierr

  real,allocatable :: lp_data(:,:) ,lp_axis(:),lp_proc_data(:,:)

  integer :: nlines,npts,ndata,ncols

  real :: tmin,tmax,deltar,chisq_lim
  character*5 :: tmins,tmaxs
  character*4 :: iform
  character*5 :: chisq
  integer :: exp_tmin,exp_tmax
  character*256 :: arg



  ! Initialization for current LP data file format
  ncols =9

  !
  ! Read and assign command line arguments
  !
  ! Command line is:  'file name'   tmin    tmax
  !
  call getarg(1,arg)
  infilename = trim(arg)

  call getarg(2,arg)
  read(arg,*) tmin
  write(0,*) 'TMIN=',tmin

  call getarg(3,arg)
  read(arg,*) tmax
  write(0,*) 'TMAX=',tmax


  call getarg(4,arg)
  read(arg,*) chisq_lim
  write(chisq,'(f5.3)') chisq_lim
  write(0,*) 'Chisq:',trim(chisq),':'


  !infilename = '134083_tab.dat'
  !tmin = 2300.0
  !tmax = 4200.0
  !chisq_lim = 0.1

  deltar = 0.005
  


  call find_free_unit_number(iunit)

  open(iunit,file=trim(infilename),iostat=ierr)
  write(0,*) 'Opening file:',trim(infilename),' Status = ',ierr


  call read_lp_data_file(iunit,lp_data,nlines,ncols)

  ! analyse and bin the lp_data

  call bin_lp_data_r(lp_axis,lp_proc_data,npts,ndata,lp_data,nlines,ncols,deltar,tmin,tmax,chisq_lim)



  ! OUTPUT

  call find_free_unit_number(ounit)

  exp_tmin = int(log10(tmin))+1
  !write(iform,'(a2,i1,a1)') '(i',exp_tmin,')'
  !write(tmins(1:exp_tmin),form=iform) int(tmin)
  if (exp_tmin.eq.3) then 
     write(tmins(1:3),'(i3)') int(tmin)
  elseif (exp_tmin.eq.4) then
     write(tmins(1:4),'(i4)') int(tmin)
  endif


  write(0,*) 'FORM=',trim(iform),':TMINS=',trim(tmins),':'

  exp_tmax = int(log10(tmax))+1
  !write(iform,'(a2,i1,a1)') '(i',exp_tmax,')'
  !write(tmaxs(1:exp_tmax),form=iform) int(tmax)
  if (exp_tmax.eq.3) then 
     write(tmaxs(1:3),'(i3)') int(tmax)
  elseif (exp_tmax.eq.4) then
     write(tmaxs(1:4),'(i4)') int(tmax)
  endif

  !write(0,*) 'FORM=',trim(iform),':TMAXS=',trim(tmaxs),':'

  write(0,*) 'Times:',trim(tmins),':',trim(tmaxs),':'

  ofilename = ' '

  ofilename = 'bin_lp_'//trim(infilename)//'_'//tmins(1:exp_tmin)//'-'//tmaxs(1:exp_tmax)//'_'//trim(chisq)

  !ofilename = 'bin_lp_'//tmins(1:exp_tmin)
  !ofilename = trim(ofilename)//'-'
  !ofilename = trim(ofilename)//tmaxs(1:exp_tmax)
  !ofilename = trim(ofilename)//'_'//trim(infilename)

  write(0,*) 'Out file name:',trim(ofilename),':',ounit
  
  !ofilename = 'bin_lp_3000_3500_134083_tab.dat'
  
  open(ounit,file=ofilename,iostat=ierr)

  write(0,*) 'Out file name:',trim(ofilename),':',ierr

  write(time,'(1x,f10.2,a,f10.2)') tmin, ' TO ', tmax

  ident = trim(infilename)//' : '//trim(time)//' : '//'CHISQ < '//trim(chisq)

  call print_lp_bin_data(ounit,lp_axis,lp_proc_data,npts,ndata,ident)

  close(iunit)
  close(ounit)

end program proclp






