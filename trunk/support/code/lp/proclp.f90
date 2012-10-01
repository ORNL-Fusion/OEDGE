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

  character*512 :: infilename,ident,time,ofilename,elm_filename

  integer :: iunit, ounit,ierr

  real,allocatable :: lp_data(:,:) ,lp_axis(:),lp_proc_data(:,:,:),lp_axis_psi(:,:),lp_axis_r(:,:),elm_fractions(:,:)


  integer :: nlines,npts,ndata,ncols,nextra,n_elm_fractions,n_avs

  real :: tmin,tmax,deltar,chisq_lim
  character*5 :: tmins,tmaxs
  character*4 :: iform
  character*5 :: chisq
  integer :: exp_tmin,exp_tmax,nargs
  character*512 :: arg

  integer,external :: iargc

  integer :: arg_cnt
  logical :: elm_filt,remove_outlier
  real :: outlier_mult
  real :: elm_start_input,elm_end_input, elm_effect_start_input, elm_effect_end_input
  character*3 :: or_ext


  ! Set base number of average intervals to 1
  n_avs =1 
  n_elm_fractions = 0
  !
  elm_filt = .false.
  remove_outlier = .false.
  or_ext = ''

  ! default outlier detection condition > 3 x pre-average
  outlier_mult = 3.0


  ! Initialization for current LP data file format
  ! Add three columns for including additional ELM information
  ncols  = 9 
  nextra = 3
  !
  ! Read and assign command line arguments
  !
  ! Command line is:  'file name'   tmin    tmax   chi_lim -e 'elm file name'
  !
  nargs = iargc()

  if (nargs.lt.3) then 
     write(0,'(a)') 'proclp command usage: proclp <infilename> tmin  tmax  chi_lim -e  <elmtimefilename>'
     write(0,'(a)') 'The -e argument is optional'
     STOP 'Insufficient command line arguments'
  endif

  arg_cnt = 1

  do while (arg_cnt.le.nargs) 

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

     elseif (arg_cnt.eq.5) then 
        read(arg,*) chisq_lim
        write(chisq,'(f5.3)') chisq_lim
        write(0,*) 'Chisq:',trim(chisq),':'

     elseif (arg.eq.'-e') then 

        elm_filt = .true.
        call getarg(arg_cnt,arg)
        arg_cnt = arg_cnt + 1

        elm_filename = trim(arg)
        write(0,*) 'ELM=',trim(elm_filename)
        ! If ELM filtering is on - several more averaging intervals are calculated
        ! a) ELM vs. no-ELM  +2
        n_avs = n_avs + 2
        ! b) inter ELM fraction times  +n  (precoded as 5 for now)
        ! add bin averaging for fractions of ELM cycle
        n_elm_fractions = 5
        n_avs = n_avs+n_elm_fractions

        if (allocated(elm_fractions)) deallocate(elm_fractions)
        allocate(elm_fractions(n_elm_fractions,2),stat=ierr)
        if (ierr.ne.0) then 
           call errmsg('PROCLP','ERROR ALLOCATING ELM_FRACTIONS')
           stop
        endif

        ! Set array of ELM cycle averages
        ! These are manually set to allow for any fractions to be specified
        ! bin 1
        elm_fractions(1,1) = 0.1
        elm_fractions(1,2) = 0.3
        ! bin 2
        elm_fractions(2,1) = 0.2
        elm_fractions(2,2) = 0.4
        ! bin 3
        elm_fractions(3,1) = 0.4
        elm_fractions(3,2) = 0.6
        ! bin 4
        elm_fractions(4,1) = 0.6
        elm_fractions(4,2) = 0.8
        ! bin 
        elm_fractions(5,1) = 0.8
        elm_fractions(5,2) = 0.99

     elseif (arg.eq.'-or') then 

        remove_outlier = .true. 
        or_ext='_or'

     elseif (arg.eq.'-om') then 

        ! note -om implies -or so both are not needed
        remove_outlier = .true. 
        or_ext='_or'

        call getarg(arg_cnt,arg)
        arg_cnt = arg_cnt + 1
        read(arg,*) outlier_mult
        write(0,'(a,f6.2)') 'OULIER MULTIPLIER = ', outlier_mult

     endif


  end do


  ! Specify deltr for bin averaging 
  deltar = 0.005


  ! Set ELM time window characteristics
  ! ELM start and ELM end are used to try to pin down the ELM peak window
  ! ELM_effect start and end are used to remove data from the analysis that may be influenced by nearby ELMs
  ! try elm time window of -2 to +4
  elm_start_input =  0.5
  elm_end_input   =  1.5
  elm_effect_start_input = -1.0
  elm_effect_end_input   =  5.0



  call find_free_unit_number(iunit)

  open(iunit,file=trim(infilename),iostat=ierr)
  write(0,*) 'Opening file:',trim(infilename),' Status = ',ierr

  call read_lp_data_file(iunit,lp_data,nlines,ncols,nextra)

  ! Add ELM information to the tabulated LP data
  if (elm_filt) then 
     call filter_lp_data(elm_filename,lp_data,nlines,ncols,nextra,elm_start_input,elm_end_input,elm_effect_start_input,elm_effect_end_input)
  endif

  ! analyse and bin the lp_data

  call bin_lp_data_r(lp_axis,lp_axis_psi,lp_axis_r,lp_proc_data,npts,ndata,lp_data,nlines,ncols,nextra,deltar,tmin,tmax,chisq_lim,elm_filt,remove_outlier,outlier_mult,n_avs,n_elm_fractions,elm_fractions)


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

  if (remove_outlier) then 
     ofilename = 'bin_'//trim(infilename)//'_'//tmins(1:exp_tmin)//'-'//tmaxs(1:exp_tmax)//'_'//trim(chisq)//trim(or_ext)
  else
     ofilename = 'bin_'//trim(infilename)//'_'//tmins(1:exp_tmin)//'-'//tmaxs(1:exp_tmax)//'_'//trim(chisq)
  endif
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

  call print_lp_bin_data(ounit,lp_axis,lp_axis_r,lp_axis_psi,lp_proc_data,npts,ndata,ident,n_avs,elm_filt,n_elm_fractions,elm_fractions)

  close(iunit)
  close(ounit)


  if (remove_outlier) then 
     ofilename = 'lp_rev_'//tmins(1:exp_tmin)//'-'//tmaxs(1:exp_tmax)//'_'//trim(chisq)//trim(or_ext)//'_'//trim(infilename)
  else
     ofilename = 'lp_rev_'//tmins(1:exp_tmin)//'-'//tmaxs(1:exp_tmax)//'_'//trim(chisq)//'_'//trim(infilename)
  endif

  open(ounit,file=ofilename,iostat=ierr)

  call print_lp_data(ounit,lp_data,nlines,ncols,nextra,ident,tmin,tmax,chisq_lim)

  close(ounit)


end program proclp






