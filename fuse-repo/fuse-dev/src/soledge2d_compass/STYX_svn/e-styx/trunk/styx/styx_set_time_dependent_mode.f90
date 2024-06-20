subroutine styx_set_time_dependent_mode(timed)
  use eirmod_precision
  !use eirmod_parmmod
  use styx2eirene
  use eirmod_comprt, only : iunout 
  implicit none
  logical,intent(in) :: timed
  real*8 :: dum1(12,1),dum2(0:1)
  integer*4 :: dum3(9,1)
  real*8 :: fluxes
  integer*4 :: iprnl,npartt,mpartt,i,j

  if (timed) then
! time step (modified to the actual value later on)
    dt_eirene=1e-7_dp
    ntime_styx=1
! number of particles in time dependent array
!    nprnli_styx=Npart_Eirene(1)
! number of histories continued from the previous time cycle
!    nptst_styx=Npart_Eirene(1)

  if (nptst_styx>nprnli_styx) then
    write(iunout,*) 'Number of particles to ba launched > census array size, defaulted to the latter'
  endif


! total number of time steps for particle tracing 
    ntmstp_styx=1
! length of individual internal time steps
    dtimv_styx=dt_eirene
! initial time (NOT irrelevant, initialized to zero here updated in infcop with icalleir counter)
    time0_styx=0.d0
    icalleir=0
! number of snapshot tallies
    nsnvi_styx=0

! see parmod (variables in comprt)

    npartt=12
    mpartt=9
    iprnl=0
    fluxes=0.d0

! initialize census array = initial condition for neutral velocity distribution
! iprnl is the number of scores on the census array at the previous time step
! the flux is set to zero to ensure that no particles are launched

    open(unit=15,access='sequential',form='unformatted')
    rewind(15)
    write(15) iprnl,fluxes,dtimv_styx
    write(15) ((dum1(j,i),j=1,npartt),i=1,iprnl)
    write(15) (dum2(i)              ,i=0,iprnl)
    write(15) ((dum3(j,i),j=1,mpartt),i=1,iprnl)
    close(15)

    write(*,*) 'EIRENE run in time dependent mode, dt = ',dt_eirene
    write(*,*) 'initial condition = 0'
  else
    dt_eirene=-1.d0
    ntime_styx=0
    write(*,*) 'EIRENE run in stationnary mode'
  endif








end subroutine styx_set_time_dependent_mode
