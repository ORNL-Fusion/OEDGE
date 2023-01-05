module sic_2020_yield_data

use error_handling

implicit none

private


!
! Notes: 1) A value of -1 is inserted for data which is not present in the tabulated source data.
!        2) A value of 0 has been inserted in some cicumstances (where data is not available) in order
!           to facilitate the calculation of the angular weighted average yield. 
!

! Carbon (SiC) data
! D -> SiC

integer,parameter :: n_c_sic_d_e = 43

real :: c_sic_d_yields(n_c_sic_d_e)
real :: c_sic_d_e(n_c_sic_d_e)

data c_sic_d_e /12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 35, 40, 50, 55, 60, 65, 70, 75, 80, 90, 100, &
110, 120, 130, 140, 160, 180, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 1200, 1600, 2400, 3200, 6400, &
12000, 24000, 48000/

data c_sic_d_yields / &
       0,        0,        0,  1.00e-6, 8.00e-06, 5.80e-05, 1.52e-04, 3.44e-04, 7.17e-04, 1.33e-03, 2.64e-03, &
3.18e-03, 3.90e-03, 4.41e-03, 4.99e-03, 5.53e-03, 5.93e-03, 6.68e-03, 7.69e-03, 8.36e-03, 8.91e-03, 9.26e-03, &
9.67e-03, 1.06e-02, 1.12e-02, 1.18e-02, 1.25e-02, 1.30e-02, 1.33e-02, 1.36e-02, 1.34e-02, 1.33e-02, 1.32e-02, &
1.28e-02, 1.26e-02, 1.14e-02, 1.04e-02, 8.55e-03, 7.36e-03, 4.47e-03, 2.80e-03, 1.61e-03, 1.04e-03/  

! C -> SiC

integer,parameter :: n_c_sic_c_e = 36

real :: c_sic_c_yields(n_c_sic_c_e)
real :: c_sic_c_e(n_c_sic_c_e)

data c_sic_c_e /5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 110, 120, 130, 140, 150, 160, 180, 200, &
250, 300, 350, 400, 500, 600, 700, 800, 1200, 1600, 3200, 6400, 12000, 24000, 48000/

data c_sic_c_yields / &
       0, 8.00e-07, 7.56e-05, 5.58e-04, 1.70e-03, 3.68e-03, 6.34e-03, 9.64e-03, 1.72e-02, 2.54e-02, 3.36e-02, &
4.15e-02, 5.64e-02, 6.31e-02, 6.97e-02, 7.61e-02, 8.17e-02, 8.74e-02, 9.27e-02, 1.03e-01, 1.12e-01, 1.32e-01, &
1.49e-01, 1.62e-01, 1.74e-01, 1.93e-01, 2.09e-01, 2.21e-01, 2.30e-01, 2.55e-01, 2.68e-01, 2.79e-01, 2.63e-01, &
2.28e-01, 1.80e-01, 0/

! Si -> SiC

integer,parameter :: n_c_sic_si_e = 34

real :: c_sic_si_yields(n_c_sic_si_e)
real :: c_sic_si_e(n_c_sic_si_e)

data c_sic_si_e /5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 110, 120, 130, 140, 150, 160, 180, 200, &
250, 300, 350, 400, 500, 600, 700, 800, 1200, 1600, 3200, 6400, 12000/

data c_sic_si_yields / &
       0,        0, 1.00E-07, 4.00E-06, 4.33E-05, 1.77E-04, 4.98E-04, 1.03E-03, 2.97E-03, 6.14E-03, 1.03E-02, &
1.53E-02, 2.68E-02, 3.33E-02, 3.99E-02, 4.67E-02, 5.36E-02, 6.06E-02, 6.75E-02, 8.12E-02, 9.42E-02, 1.26E-01, &
1.53E-01, 1.79E-01, 2.01E-01, 2.40E-01, 2.73E-01, 3.02E-01, 3.26E-01, 3.99E-01, 4.49E-01, 5.51E-01, 6.19E-01, &
6.41E-01/


! Silicon (SiC) data
! D -> SiC

integer,parameter :: n_si_sic_d_e = 43

real :: si_sic_d_yields(n_si_sic_d_e)
real :: si_sic_d_e(n_si_sic_d_e)

data si_sic_d_e /12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 35, 40, 50, 55, 60, 65, 70, 75, 80, 90, 100, &
110, 120, 130, 140, 160, 180, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 1200, 1600, 2400, 3200, 6400, &
12000, 24000, 48000/

data si_sic_d_yields / &
       0,        0,        0,        0,        0, 2.00e-06, 8.00e-06, 1.10e-05, 5.50e-05, 1.59e-04, 5.83e-04, &
9.73e-04, 1.32e-03, 1.68e-03, 2.20e-03, 2.61e-03, 3.06e-03, 3.70e-03, 4.57e-03, 5.25e-03, 6.00e-03, 6.53e-03, &
7.24e-03, 8.00e-03, 8.58e-03, 9.34e-03, 1.04e-02, 1.10e-02, 1.16e-02, 1.17e-02, 1.19e-02, 1.18e-02, 1.18e-02, &
1.17e-02, 1.13e-02, 1.03e-02, 9.46e-03, 7.75e-03, 6.86e-03, 4.19e-03, 2.63e-03, 1.59e-03, 9.64e-04/

! Si -> SiC

integer,parameter :: n_si_sic_si_e = 34

real :: si_sic_si_yields(n_si_sic_si_e)
real :: si_sic_si_e(n_si_sic_si_e)

data si_sic_si_e /5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 110, 120, 130, 140, 150, 160, 180, 200, &
250, 300, 350, 400, 500, 600, 700, 800, 1200, 1600, 3200, 6400, 12000/

data si_sic_si_yields / &
       0,        0,        0, 6.00E-07, 6.00E-06, 2.89E-05, 1.20E-04, 3.02E-04, 1.12E-03, 2.70E-03, 5.19E-03, &
8.34E-03, 1.64E-02, 2.10E-02, 2.61E-02, 3.13E-02, 3.64E-02, 4.18E-02, 4.74E-02, 5.82E-02, 6.90E-02, 9.49E-02, &
1.18E-01, 1.40E-01, 1.59E-01, 1.94E-01, 2.22E-01, 2.47E-01, 2.68E-01, 3.32E-01, 3.74E-01, 4.63E-01, 5.22E-01, &
5.41E-01/

! C -> SiC

integer,parameter :: n_si_sic_c_e = 36

real :: si_sic_c_yields(n_si_sic_c_e)
real :: si_sic_c_e(n_si_sic_c_e)

data si_sic_c_e /5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 110, 120, 130, 140, 150, 160, 180, 200, &
250, 300, 350, 400, 500, 600, 700, 800, 1200, 1600, 3200, 6400, 12000, 24000, 48000/

data si_sic_c_yields / &
       0,        0, 4.80E-06, 9.90E-05, 4.97E-04, 1.42E-03, 2.96E-03, 5.06E-03, 1.04E-02, 1.68E-02, 2.32E-02, &
2.96E-02, 4.16E-02, 4.72E-02, 5.28E-02, 5.78E-02, 6.25E-02, 6.72E-02, 7.16E-02, 7.97E-02, 8.74E-02, 1.04E-01, &
1.18E-01, 1.29E-01, 1.39E-01, 1.56E-01, 1.69E-01, 1.80E-01, 1.88E-01, 2.10E-01, 2.21E-01, 2.32E-01, 2.20E-01, &
1.94E-01, 1.55E-01, 0/


public :: yld_sic, print_sic_yields

    integer,parameter ::  maxdata=100
    !real :: yield_self(maxdata),energy_self(maxdata)
    real :: yield_d(maxdata),energy_d(maxdata)

    real :: yield_si(maxdata),energy_si(maxdata)
    real :: yield_c(maxdata),energy_c(maxdata)

    integer :: loaded_mat = 0
    
    ! samod - Think loaded_opt was never actually defined, so doing that
    real :: loaded_opt = -1000.0

    !integer :: n_data_self = 0
    integer :: n_data_d = 0

    integer :: n_data_si = 0
    integer :: n_data_c = 0


    !save yield_self,energy_self,n_data_self
    save yield_d,energy_d,n_data_d
    save yield_c,energy_c,n_data_c
    save yield_si,energy_si,n_data_si

    save loaded_opt,loaded_mat
                                                                                                  
contains

  real function yld_sic(energy,matp,matt)
    implicit none
    integer,intent(in) :: matp,matt
    real energy

    ! MATT = 20 = C (SiC)
    ! MATT = 21 = Si (SiC)
    !
    ! MATP = 2    D
    ! MATP = 6    Self
    ! MATP = 5    C
    ! MATP = 8   Si
    !

    integer in,ipos
    external ipos
    !real lin_interp
    !external lin_interp

    ! Initialize
    !
    yld_sic = 0.0
    !
    ! Load data if not already loaded
    !
    if (loaded_mat.ne.matt) then 
       ! Load sputtering data
       !call get_yields(matt,maxdata,n_data_self,energy_self,yield_self,&
       !     n_data_d,energy_d,yield_d)

       call get_yields(matt,maxdata,n_data_c,energy_c,yield_c,n_data_si,energy_si,yield_si,n_data_d,energy_d,yield_d)

       ! Save loaded data option
       !write(0,*) 'Loading:',matt,n_data_self,n_data_d,loaded_opt,opt,loaded_mat
       loaded_mat = matt
    endif


    if (matp.eq.2) then 
       !
       ! Sputtering by D
       !
       yld_sic = lin_interp(energy,energy_d,yield_d,n_data_d,1)
    elseif (matp.eq.5) then 
!       write(0,*) 'Calculating sputtering yield due to C'
       !
       ! Sputtering by C
       !
       yld_sic = lin_interp(energy,energy_c,yield_c,n_data_c,1)
    elseif (matp.eq.8) then
!       write(0,*) 'Calculating sputtering yield due to Si'
       !
       ! Sputtering by Si
       !
       !write(0,*) 'Calculating Y using MATP=8...'
       yld_sic = lin_interp(energy,energy_si,yield_si,n_data_si,1)
    endif


    return
  end function yld_sic


  real function lin_interp(val,lookup,result_data,ndata,minopt)
    implicit none
    integer ndata
    integer, intent(in) ::  minopt
    real val,lookup(ndata),result_data(ndata)

    integer in,ipos
    external ipos

    if (val.lt.lookup(1)) then 
       if (minopt.eq.1) then 
          lin_interp = 0.0
       elseif (minopt.eq.2) then 
          lin_interp = result_data(1)
       endif
    elseif (val.eq.lookup(1)) then 
       lin_interp = result_data(1)
    elseif (val.ge.lookup(ndata)) then 
       lin_interp = result_data(ndata)
    else
       in = ipos(val,lookup,ndata)
       lin_interp = &
            (val-lookup(in-1))/(lookup(in)-lookup(in-1))* &
            (result_data(in)-result_data(in-1)) + &
             result_data(in-1)
    endif

    return
  end function lin_interp


  subroutine get_yields(matt,maxdata,n_data_c,energy_c,yield_c,n_data_si,energy_si,yield_si,n_data_d,energy_d,yield_d)

    implicit none
    integer :: matt
    integer :: n_data_d,maxdata,n_data_c,n_data_si
    !integer :: n_data_self

    real energy_d(maxdata),yield_d(maxdata)
    !real energy_self(maxdata),yield_self(maxdata)
    real energy_c(maxdata),yield_c(maxdata)
    real energy_si(maxdata),yield_si(maxdata)

    ! Yield data is taken from TRIM.SP runs (Bringuier, 2020).
    !
    ! MATT = 20    is C (SiC)
    ! MATT = 21    is Si (SiC)
    !
    !write(0,*) 'Get yield:',matt
    ! Load C (SiC) data 
    if (matt.eq.20) then 

       ! load self-sputtering data

       !write(0,*) 'SELF:',n_be_self_angle,n_be_self_e
       ! IPP/08 Krieger - SUN compiler hates more than 132 columns
       !call load_yield(n_data_self,maxdata,energy_self,yield_self,n_c_sic_self_e,c_sic_self_e,c_sic_self_yields)

       ! load C yields by D
       !write(0,*) 'D:',n_be_d_angle,n_be_d_e
       call load_yield(n_data_d,maxdata,energy_d,yield_d,n_c_sic_d_e,c_sic_d_e,c_sic_d_yields)

       ! load C yields by C
       call load_yield(n_data_c,maxdata,energy_c,yield_c,n_c_sic_c_e,c_sic_c_e,c_sic_c_yields)

       ! load C yields by Si
       call load_yield(n_data_si,maxdata,energy_si,yield_si,n_c_sic_si_e,c_sic_si_e,c_sic_si_yields)

    ! Load Si (SiC) data
    elseif (matt.eq.21) then 

       ! load self-sputtering data
       !call load_yield(n_data_self,maxdata,energy_self,yield_self,n_si_sic_self_e,si_sic_self_e,si_sic_self_yields)

       ! load Si yields by D
       call load_yield(n_data_d,maxdata,energy_d,yield_d,n_si_sic_d_e,si_sic_d_e,si_sic_d_yields)

       ! load Si yields by C
       call load_yield(n_data_c,maxdata,energy_c,yield_c,n_si_sic_c_e,si_sic_c_e,si_sic_c_yields)

       ! load Si yields by Si
       call load_yield(n_data_si,maxdata,energy_si,yield_si,n_si_sic_si_e,si_sic_si_e,si_sic_si_yields)

    endif

    return
  end subroutine get_yields


  subroutine load_yield(ndata,maxdata,energy,yield,n_src_e,src_e,src_yields)
    implicit none
    integer :: ndata,maxdata
    real    :: energy(maxdata), yield(maxdata)

    integer :: n_src_e
    real :: src_e(n_src_e)
    real :: src_yields(n_src_e)


    integer :: in,data_count

    ! load data for incident_angle from src_yields into yield array

    !write(0,'(a,3i5,3(1x,g12.5))') 'Loading yields:',ia,n_src_e,n_src_angle,incident_angle,src_angle(ia)

    ! Matching data not found - stop

    !
    ! Load requested data - leave out -1 entries - if average data is requested check for at least 5 good data values in average. 
    !

    data_count = 0

    do in = 1,n_src_e

          if (src_yields(in).gt.0.0) then 
             data_count = data_count + 1
             energy(data_count) = src_e(in)
             yield(data_count) = src_yields(in)
             !write(0,'(a,i6,1x,f10.2,1x,g12.5)') 'Loading:',data_count,energy(data_count),yield(data_count)
          endif
   
    end do

    ndata = data_count

    return
  end subroutine load_yield


  ! sazmod - Seems the _2020_ was not intended, so renaming
  !subroutine print_sic_2020_yields(ounit)
  subroutine print_sic_yields(ounit)
    implicit none
    !
    ! Prints the yield data currently in use to the specified unit
    !
    integer :: ounit

    integer :: in

    write(ounit,'(a)') '     2020 SiC SPUTTERING DATA:'
    write(ounit,'(a)') '     BASE SPUTTERING YIELDS BY D IN USE (MAY BE MULTIPLIED BY ANY SPECIFIED YIELD MULTIPLIERS)'
    write(ounit,'(25x,a,12x,a)')  'E (eV)','YIELD'
    do in = 1,n_data_d
       write(ounit,'(10x,i6,3x,f12.2,3x,g12.5)')  in,energy_d(in),yield_d(in)
    end do

    !write(ounit,'(a)')
    !write(ounit,'(a)') '     BASE SELF SPUTTERING YIELDS IN USE (MAY BE MULTIPLIED BY ANY SPECIFIED YIELD MULTIPLIERS)'
    !write(ounit,'(25x,a,12x,a)')  'E (eV)','YIELD'
    !do in = 1,n_data_self
    !   write(ounit,'(10x,i6,3x,f12.2,3x,g12.5)')  in,energy_self(in),yield_self(in)
    !end do

    write(ounit,'(a)')
    write(ounit,'(a)') '     BASE SPUTTERING YIELDS BY C IN USE (MAY BE MULTIPLIED BY ANY SPECIFIED YIELD MULTIPLIERS)'
    write(ounit,'(25x,a,12x,a)')  'E (eV)','YIELD'
    do in = 1,n_data_c
       write(ounit,'(10x,i6,3x,f12.2,3x,g12.5)')  in,energy_c(in),yield_c(in)
    end do

    write(ounit,'(a)')
    write(ounit,'(a)') '     BASE SPUTTERING YIELDS BY SI IN USE (MAY BE MULTIPLIED BY ANY SPECIFIED YIELD MULTIPLIERS)'
    write(ounit,'(25x,a,12x,a)')  'E (eV)','YIELD'
    do in = 1,n_data_si
       write(ounit,'(10x,i6,3x,f12.2,3x,g12.5)')  in,energy_si(in),yield_si(in)
    end do

  !end subroutine print_sic_2020_yields
  end subroutine print_sic_yields


end module sic_2020_yield_data
