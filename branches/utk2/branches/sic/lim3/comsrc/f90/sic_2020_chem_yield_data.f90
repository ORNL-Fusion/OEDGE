module sic_2020_chem_yield_data

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

integer,parameter :: n_c_sic_d_e = 6
integer,parameter :: n_c_sic_d_temp = 12

real :: c_sic_d_yields(n_c_sic_d_temp,n_c_sic_d_e)
real :: c_sic_d_temp(n_c_sic_d_temp)
real :: c_sic_d_e(n_c_sic_d_e)

data c_sic_d_e /15, 20, 30, 50, 100, 300/

data c_sic_d_temp /293, 323, 343, 373, 473, 573, 673, 873, 973, 1073, 1173, 1273/

data c_sic_d_yields / &
1.25e-03, 1.25e-03, 1.25e-03, 1.25e-03, 1.25e-03, 1.25e-03, 1.24e-03, 1.54e-03, 1.34e-03, 4.68e-03, 1.31e-04, 3.40e-05, &
1.49e-03, 1.49e-03, 1.49e-03, 1.49e-03, 1.49e-03, 1.49e-03, 1.49e-03, 1.76e-03, 1.46e-03, 5.08e-04, 1.42e-03, 3.68e-05, &
1.74e-03, 1.74e-03, 1.74e-03, 1.74e-03, 1.74e-03, 1.74e-03, 1.74e-03, 2.15e-03, 1.86e-03, 6.51e-04, 1.83e-04, 4.73e-05, &
1.77e-03, 1.77e-03, 1.77e-03, 1.77e-03, 1.77e-03, 1.77e-03, 1.77e-03, 2.64e-03, 2.63e-03, 9.39e-04, 2.64e-04, 6.85e-05, &
1.02e-03, 1.02e-03, 1.02e-03, 1.02e-03, 1.02e-03, 1.02e-03, 1.02e-03, 2.75e-03, 3.48e-03, 1.29e-03, 3.64e-04, 9.45e-05, &
9.77e-06, 9.77e-06, 9.77e-06, 9.77e-06, 9.77e-06, 9.80e-06, 1.46e-05, 2.52e-03, 4.04e-03, 1.53e-03, 4.35e-04, 1.13e-04/
  

! Silicon (SiC) data
! D -> SiC

integer,parameter :: n_si_sic_d_e = 6
integer,parameter :: n_si_sic_d_temp = 12

real :: si_sic_d_yields(n_si_sic_d_temp,n_si_sic_d_e)
real :: si_sic_d_temp(n_si_sic_d_temp)
real :: si_sic_d_e(n_si_sic_d_e)

data si_sic_d_e /15, 20, 30, 50, 100, 300/

data si_sic_d_temp /293, 323, 343, 373, 473, 573, 673, 873, 973, 1073, 1173, 1273/

data si_sic_d_yields / &
       0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0, &
       0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0, &
       0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0, &
       0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0, &
       0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0, &
       0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0/


public :: yld_sic_chem, print_sic_chem_yields

    integer,parameter ::  maxenergydata=100
    integer,parameter ::  maxtempdata=100
    real :: yield_d(maxtempdata,maxenergydata),energy_d(maxenergydata),temp_d(maxtempdata)

    integer :: loaded_mat = 0

    integer :: n_energy_data_d = 0
    integer :: n_temp_data_d = 0

    save yield_d,energy_d,temp_d,n_energy_data_d,n_temp_data_d
    save loaded_mat
                                                                                                  
contains


! Containing function for calculating the chemical sputtering yield from SiC
  real function yld_sic_chem(energy,matp,matt,temp)
    implicit none
    integer,intent(in) :: matp,matt
    real :: energy,temp

    ! MATT = 20 = C (SiC)
    ! MATT = 21 = Si (SiC)
    
    ! MATP = 2    D
    
    integer :: in

    ! Variables for linear interpolation
    real :: yield_temp_1,yield_temp_2
    real :: temp_real,energy_real
    logical :: flag_temp,flag_energy
    integer :: j,i_temp,i_energy

    ! Initialize
    !
    yld_sic_chem = 0.0
    !
    ! Load data if not already loaded
    !
    if (loaded_mat.ne.matt) then 
       ! Load sputtering data
       call get_yields(matt,maxtempdata,maxenergydata,n_temp_data_d,n_energy_data_d,temp_d,energy_d,yield_d)
       ! Save loaded data option
       !write(0,*) 'Loading:',matt,n_data_self,n_data_d,loaded_opt,opt,loaded_mat
       loaded_mat = matt
    endif

    ! Check value of yield
    !write(0,*) 'yield (yld_sic_chem) =', yield_d(4,4)
    !write(0,*) '# of energy data =', n_energy_data_d
    !write(0,*) '# of temp data =', n_temp_data_d

    if (matp.eq.2) then ! for chemical sputtering by D

       ! Linear interpolation function goes here

       ! Initialization
       flag_temp=.TRUE.
       flag_energy=.TRUE.
       temp_real=0.0
       energy_real=0.0
       j=1

       ! Interpolate over temperature
       if (temp.lt.temp_d(1)) then
          i_temp=1
          temp_real=temp          
          temp=temp_d(1)
       elseif (temp.gt.temp_d(n_temp_data_d)) then
          i_temp=n_temp_data_d-1
          temp_real=temp
          temp=temp_d(n_temp_data_d)
       else
          do while (flag_temp)
             if ((temp.ge.temp_d(j)).and.(temp.le.temp_d(j+1))) then
                i_temp=j
                flag_temp=.FALSE.
             else
                j=j+1 
             endif     
          enddo
       endif

       ! Interpolate over energy
       j=1
       if (energy.lt.energy_d(1)) then         
          i_energy=1
          energy_real=energy
          energy=energy_d(1)                    
       elseif (energy.gt.energy_d(n_energy_data_d)) then
          i_energy=n_energy_data_d-1
          energy_real=energy
          energy=energy_d(n_energy_data_d)
       else       
          do while (flag_energy)
             if ((energy.ge.energy_d(j)).and.(energy.le.energy_d(j+1))) then
                i_energy=j
                flag_energy=.FALSE.
             else
                j=j+1 
             endif     
          enddo
       endif

       ! Find interpolated yield
       yield_temp_1=(yield_d(i_temp+1,i_energy)-yield_d(i_temp,i_energy))/(temp_d(i_temp+1)- &
          temp_d(i_temp))*(temp-temp_d(i_temp))+yield_d(i_temp,i_energy)
       yield_temp_2=(yield_d(i_temp+1,i_energy+1)-yield_d(i_temp,i_energy+1))/(temp_d(i_temp+1)- &
          temp_d(i_temp))*(temp-temp_d(i_temp))+yield_d(i_temp,i_energy+1)
       yld_sic_chem=(yield_temp_2-yield_temp_1)/(energy_d(i_energy+1)-energy_d(i_energy))* &
          (energy-energy_d(i_energy))+yield_temp_1

       if (temp_real.ne.0.0) then
          temp=temp_real
       endif
       if (energy_real.ne.0.0) then
          energy=energy_real
       endif

       !write(0,*) 'Temperature =', temp
       !write(0,*) 'Energy =', energy
       !write(0,*) 'Chemical Yield =', yld_sic_chem

       !yld_sic_chem = lin_interp(energy,temp,energy_d,temp_d,yield_d,n_energy_data_d,n_temp_data_d,2)
    endif

    return
  end function yld_sic_chem


! Linear interpolation function - interpolates over both temperature and energy
  real function lin_interp(energy,temp,energy_table,temp_table,yield_table,n_energy_data,n_temp_data,minopt)
    implicit none
    integer n_energy_data,n_temp_data
    integer :: temp_index
    integer, intent(in) ::  minopt
    real :: energy,temp
    !real, dimension(n_energy_data) :: energy_table
    !real, dimension(n_temp_data) :: temp_table
    !real, dimension(n_temp_data,n_energy_data) :: yield_table
    real energy_table(n_energy_data)
    real temp_table(n_temp_data)
    real yield_table(n_temp_data,n_energy_data)
    real :: yield_temp_1,yield_temp_2
    real :: temp_real,energy_real
    logical :: flag_temp,flag_energy
    integer :: in,j,i_temp,i_energy

!   Initialization
    flag_temp=.TRUE.
    flag_energy=.TRUE.
    temp_real=0.0
    energy_real=0.0
    j=1

!    write(0,*) 'Temperature Table:', temp_table
!    write(0,*) 'Energy Table:', energy_table
!    write(0,*) 'Yield Table:', yield_table

       
!   interp(x_value,x_table,x_table_length,y_table,y_table_length)

!    if (temp.lt.temp_table(1)) then
!       if (energy.lt.energy_table(1)) then
!          lin_intep=yield_table(1,1)
!       elseif (energy.gt.energy_table(n_energy_data))
!          lin_interp=yield_table(1,n_energy_data)   
!
!       lin_interp=interp(energy,energy_table,n_energy_data,yield_table(1,1:n_energy_data),n_energy_data)

    ! Check value of yield
    !write(0,*) '# temp data', n_temp_data
    !write(0,*) '# energy data', n_energy_data
    !write(0,*) 'energy (lin_interp) =' , energy_table(4)
    !write(0,*) 'temp (lin_interp) =', temp_table(4)
    !write(0,*) 'yield (lin_interp) =', yield_table(4,4)    

!   Interpolate over temperature
    if (temp.lt.temp_table(1)) then
       !write(0,*) 'Temp is minimum'
       if (minopt.eq.2) then
          i_temp=1
          temp_real=temp          
          temp=temp_table(1)
       endif
    elseif (temp.gt.temp_table(n_temp_data)) then
       !write(0,*) 'Temp is maximum'
       i_temp=n_temp_data-1
       temp_real=temp
       temp=temp_table(n_temp_data)
    else
       !write(0,*) 'Temp is in the middle'
       do while (flag_temp)
          if ((temp.ge.temp_table(j)).and.(temp.le.temp_table(j+1))) then
             i_temp=j
             flag_temp=.FALSE.
          else
             j=j+1 
          endif     
       enddo
    endif

!   Interpolate over energy
    j=1
    if (energy.lt.energy_table(1)) then
       !write(0,*) 'Energy is minimum'
       if (minopt.eq.2) then
          i_energy=1
          energy_real=energy
          energy=energy_table(1)          
       endif
    elseif (energy.gt.energy_table(n_energy_data)) then
       !write(0,*) 'Energy is maximum'
       i_energy=n_energy_data-1
       energy_real=energy
       energy=energy_table(n_energy_data)
    else
       !write(0,*) 'Energy is in the middle'
       do while (flag_energy)
          if ((energy.ge.energy_table(j)).and.(energy.le.energy_table(j+1))) then
             i_energy=j
             flag_energy=.FALSE.
          else
             j=j+1 
          endif     
       enddo
    endif

!   Find interpolated yield
    yield_temp_1=(yield_table(i_temp+1,i_energy)-yield_table(i_temp,i_energy))/(temp_table(i_temp+1)- &
       temp_table(i_temp))*(temp-temp_table(i_temp))+yield_table(i_temp,i_energy)
    yield_temp_2=(yield_table(i_temp+1,i_energy+1)-yield_table(i_temp,i_energy+1))/(temp_table(i_temp+1)- &
       temp_table(i_temp))*(temp-temp_table(i_temp))+yield_table(i_temp,i_energy+1)
    lin_interp=(yield_temp_2-yield_temp_1)/(energy_table(i_energy+1)-energy_table(i_energy))* &
       (energy-energy_table(i_energy))+yield_temp_1

    if (temp_real.ne.0.0) then
       temp=temp_real
    endif
    if (energy_real.ne.0.0) then
       energy=energy_real
    endif

!   Write results to error file
    !write(0,*) 'Yield(373K,15eV) =', yield_table(4,1)
    !write(0,*) 'Yield(473K,15eV) =', yield_table(5,1)
    !write(0,*) 'Yield(373K,20eV) =', yield_table(4,2)
    !write(0,*) 'Yield(473K,20eV) =', yield_table(5,2)
    !write(0,*) 'Yield(373K,30eV) =', yield_table(4,3)
    !write(0,*) 'Yield(373K,50eV) =', yield_table(4,4)
    !write(0,*) 'Yield(373K,100eV) =', yield_table(4,5)
    !write(0,*) 'Temperature =', temp
    !write(0,*) 'Energy =', energy
    !write(0,*) 'Temperature Index =', i_temp
    !write(0,*) 'Energy Index =', i_energy
    !write(0,*) 'Chemical Yield =', lin_interp

    return
  end function lin_interp


! Subroutine used to pull appropriate tables for chemical sputtering yield calculation
  subroutine get_yields(matt,maxtempdata,maxenergydata,n_temp_data_d,n_energy_data_d,temp_d,energy_d,yield_d)
    implicit none
    integer :: matt
    integer :: n_temp_data_d,n_energy_data_d,maxtempdata,maxenergydata

    real :: energy_d(maxenergydata),temp_d(maxtempdata),yield_d(maxtempdata,maxenergydata)

    ! Yield data is taken from modified Roth model (Abrams, 2020).
    !
    ! MATT = 20    is C (SiC)
    ! MATT = 21    is Si (SiC)
    !
    !write(0,*) 'Get yield:',matt
    ! Load C (SiC) data 
    if (matt.eq.20) then 

       ! load deuterium yields
       !write(0,*) 'D:',n_be_d_angle,n_be_d_e
       call load_yield(n_temp_data_d,n_energy_data_d,maxtempdata,maxenergydata,temp_d,energy_d,yield_d, &
            n_c_sic_d_temp,n_c_sic_d_e,c_sic_d_temp,c_sic_d_e,c_sic_d_yields)

    ! Load Si (SiC) data
    elseif (matt.eq.21) then 

       ! load deuterium yields
       call load_yield(n_temp_data_d,n_energy_data_d,maxtempdata,maxenergydata,temp_d,energy_d,yield_d, &
            n_si_sic_d_temp,n_si_sic_d_e,si_sic_d_temp,si_sic_d_e,si_sic_d_yields)

    endif

    ! Check value of yield
    !write(0,*) 'yield (get_yields) =', yield_d(4,4)

    return
  end subroutine get_yields


  subroutine load_yield(ntempdata,nenergydata,maxtempdata,maxenergydata,temp_d,energy_d,yield_d,n_src_temp, &
       n_src_e,src_temp,src_e,src_yields)
    implicit none
    integer :: ntempdata,nenergydata,maxtempdata,maxenergydata
    real    :: temp_d(maxtempdata),energy_d(maxenergydata), yield_d(maxtempdata,maxenergydata)

    integer :: n_src_temp,n_src_e
    real :: src_temp(n_src_temp),src_e(n_src_e)
    real :: src_yields(n_src_temp,n_src_e)


    integer :: i,j,count_temp,count_energy

    ! load data for incident_angle from src_yields into yield array

    !write(0,'(a,3i5,3(1x,g12.5))') 'Loading yields:',ia,n_src_e,n_src_angle,incident_angle,src_angle(ia)

    ! Matching data not found - stop

    !
    ! Load requested data - leave out -1 entries - if average data is requested check for at least 5 good data values in average. 
    !

    count_temp = 0
    count_energy=0

    do i=1,n_src_e      
       count_energy=count_energy+1
       energy_d(i)=src_e(i)
       do j=1,n_src_temp
          temp_d(j)=src_temp(j)
          yield_d(j,i)=src_yields(j,i)

         ! Check output
          !write(0,*) 'i =', i
          !write(0,*) 'j =', j
          !write(0,*) 'yield =', yield_d(j,i)

          if (i.eq.n_src_e) then
             count_temp=count_temp+1
          endif   
          ! write(0,'(a,i6,1x,f10.2,1x,g12.5)') 'Loading:',data_count,energy(data_count),yield(data_count)         
       enddo

    enddo

    ntempdata=count_temp
    nenergydata=count_energy
    !ndata = count

    return
  end subroutine load_yield


  ! rea function interp(x_value,x_table,x_table_length,y_table,y_table_length)
  !  implicit none
  !  integer :: index,x_table_length,y_table_length
  !  real :: x_value,x_table(x_table_length),y_table(y_table_length)
  !  integer :: i

  !  do i=1,x_table_length
  !     if ((x_value.ge.x_table(i)).and.(x_value.le.x_table(i+1))) then
  !        index=i
  !     end if
  !  enddo

  !  interp=(y_table(index+1)-y_table(index))/(x_table(index+1)-x_table(index))*(x_value-x_table(index))+y_table(index)

  !end function interp


  ! sazm0d - Seems the _2020_ not what was meant to be written.
  !subroutine print_sic_2020_chem_yields(ounit)
  subroutine print_sic_chem_yields(ounit)
    implicit none
    !
    ! Prints the yield data currently in use to the specified unit
    !
    integer :: ounit

    integer :: in

    write(ounit,'(a)') '     2020 SiC CHEMICAL SPUTTERING DATA:'
    write(ounit,'(a)') '     BASE D SPUTTERING YIELDS IN USE (MAY BE MULTIPLIED BY ANY SPECIFIED YIELD MULTIPLIERS)'
    write(ounit,'(25x,a,12x,a)')  'E (eV)','YIELD'
    do in = 1,n_energy_data_d
       write(ounit,'(10x,i6,3x,f12.2,3x,g12.5)')  in,energy_d(in),yield_d(1,in)
    end do

  !end subroutine print_sic_2020_chem_yields
  end subroutine print_sic_chem_yields


end module sic_2020_chem_yield_data
