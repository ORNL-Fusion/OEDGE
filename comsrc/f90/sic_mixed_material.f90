module sic_mixed_material

implicit none
private

! This module implements the mixed material model for sputtering from
! SiC as described in Abrams Nuclear Fusion 2021. The incident angle is 
! assumed to be 45 degrees (versus normal as is often assumed in 3DLIM). 
!
! Y_Si,tot = sig_SiC * Y_SiC,Si + sig_Si * Y_Si
! Y_C,tot  = sig_SiC * Y_SiC,C  + sig_C  * Y_C
! 
! Where sig refers to the fracitonal abundances of Si, C and SiC within 
! some characteristic implantation depths of the surface.
! No other species, thus sig_Si + sig_C + sig_SiC = 1.
! Each yield is broken up into the constituient yields from each 
! impacting ion.
! For each yield calculation, see Eq's 1-6 in the paper. 


! D-->SiC, Si, physical. TrimSP calculated.
! This is the average of the yield for the crystollographic directions
! SiC(110), SiC(111) C-rich, and SiC(111) Si-rich.
integer,parameter :: n_d_sic_si = 23
real :: d_sic_si_yields(n_d_sic_si)
real :: d_sic_si_e(n_d_sic_si)
data d_sic_si_e / &
       25,    30,    35,    40,    50,    60,    70,    80,   100,   &
      120,   140,   160,   200,   300,   400,   500,   800,  1600,   &
     3200,  6400, 12000, 24000, 48000/
data d_sic_si_yields / &
  0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00,  &
  0.00E+00, 0.00E+00, 2.00E-07, 1.90E-06, 1.60E-05, 5.63E-05,  &
  2.50E-04, 1.10E-03, 1.92E-03, 2.57E-03, 3.81E-03, 4.73E-03,  &
  4.43E-03, 3.38E-03, 2.33E-03, 1.27E-03, 6.96E-04/  


! D-->SiC, Si, chemical
! zero

! C-->SiC, Si, physical
! This is the average of the yield for the crystollographic directions
! SiC(110), SiC(111) C-rich, and SiC(111) Si-rich.
integer,parameter :: n_c_sic_si = 21
real :: c_sic_si_yields(n_c_sic_si)
real :: c_sic_si_e(n_c_sic_si)
data c_sic_si_e / &
       25,    30,    35,    40,    50,    60,    70,    80,   100,   &
      120,   140,   160,   200,   300,   400,   500,   800,  1600,   &
     3200,  6400, 12000/
data c_sic_si_yields / &
  0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 3.33E-08, 2.97E-06,  &
  4.44E-05, 2.02E-04, 1.29E-03, 3.57E-03, 6.90E-03, 1.08E-02,  &
  1.97E-02, 4.15E-02, 6.00E-02, 7.52E-02, 1.08E-01, 1.52E-01,  &
  1.79E-01, 1.83E-01, 1.67E-01/

! D-->SiC, C, physical
! This is the average of the yield for the crystollographic directions
! SiC(110), SiC(111) C-rich, and SiC(111) Si-rich
integer,parameter :: n_d_sic_c = 23
real :: d_sic_c_yields(n_d_sic_c)
real :: d_sic_c_e(n_d_sic_c)
data d_sic_c_e /&
       25,    30,    35,    40,    50,    60,    70,    80,   100,   &
      120,   140,   160,   200,   300,   400,   500,   800,  1600,   &
     3200,  6400, 12000, 24000, 48000/
data d_sic_c_yields / &
  0.00E+00, 1.00E-07, 1.91E-05, 6.67E-05, 2.97E-04, 8.28E-04, &
  1.56E-03, 2.41E-03, 4.13E-03, 5.77E-03, 7.14E-03, 8.35E-03, &
  1.02E-02, 1.30E-02, 1.44E-02, 1.51E-02, 1.55E-02, 1.39E-02, &
  1.05E-02, 6.99E-03, 4.34E-03, 2.23E-03, 1.16E-03/

! D-->SiC, C, chemical
! 10% of the Roth yield.
integer,parameter :: n_d_sic_c_ch = 28
real :: d_sic_c_ch_yields(n_d_sic_c_ch)
real :: d_sic_c_ch_e(n_d_sic_c_ch)
data d_sic_c_ch_e / &
        5,    10,    15,    20,    25,    30,    35,    40,    50,  &
       60,    70,    80,   100,   120,   140,   160,   200,   300,  &
      400,   500,   600,   800,  1000,  1400,  1800,  2400,  3000,  &
     4000/
data d_sic_c_ch_yields /&
  2.74E-03, 8.64E-03, 1.25E-02, 1.49E-02, 1.64E-02, 1.74E-02,   &
  1.79E-02, 1.81E-02, 1.77E-02, 1.67E-02, 1.53E-02, 1.36E-02,   &
  1.02E-02, 7.13E-03, 4.75E-03, 3.06E-03, 1.19E-03, 9.77E-05,   &
  7.71E-06, 6.04E-07, 4.74E-08, 2.93E-10, 1.82E-12, 7.22E-17,   & 
  2.92E-21, 7.72E-28, 2.09E-34, 2.46E-45/

! C-->SiC, C, physical
integer,parameter :: n_c_sic_c = 21
real :: c_sic_c_yields(n_c_sic_c)
real :: c_sic_c_e(n_c_sic_c)
data c_sic_c_e / &
       25,    30,    35,    40,    50,    60,    70,    80,   100,   &
      120,   140,   160,   200,   300,   400,   500,   800,  1600,   &
     3200,  6400, 12000/
data c_sic_c_yields / &
  1.44E-04, 6.43E-04, 1.71E-03, 3.42E-03, 8.63E-03, 1.57E-02,   &
  2.39E-02, 3.24E-02, 4.97E-02, 6.57E-02, 8.08E-02, 9.43E-02,   &
  1.19E-01, 1.64E-01, 1.98E-01, 2.23E-01, 2.75E-01, 3.35E-01,   &
  3.64E-01, 3.53E-01, 3.09E-01/

! D-->Si, physical. TrimSP calculated.
integer,parameter :: n_d_si = 23
real :: d_si_yields(n_d_si)
real :: d_si_e(n_d_si)
data d_si_e / &
     25,    30,    35,    40,    50,    60,    70,    80,    100,  &
    120,   140,   160,   200,   300,   400,   500,   800,   1600,  &
    3200, 6400, 12000, 24000, 48000/
data d_si_yields / &
  7.40E-06, 2.11E-04, 8.11E-04, 1.78E-03, 4.52E-03, 7.76E-03, &
  1.12E-02, 1.45E-02, 2.04E-02, 2.57E-02, 3.04E-02, 3.42E-02, &
  4.08E-02, 5.14E-02, 5.70E-02, 6.00E-02, 6.25E-02, 5.70E-02, &
  4.50E-02, 3.09E-02, 1.85E-02, 9.25E-03, 5.18E-03/

! D-->Si, chemical (room temperature, 25C).
integer,parameter :: n_d_si_ch = 28
real :: d_si_ch_yields(n_d_si_ch)
real :: d_si_ch_e(n_d_si_ch)
data d_si_ch_e / &
      5,    10,    13,    16,    21,    28,    36,    48,    63,  &
     83,   110,   145,   191,   251,   331,   437,   575,   759,  &
   1000,  1318,  1738,  2291,  3020,  3981,  5248,  6918,  8318,  &
   9120/
data d_si_ch_yields / &
  7.34E-04, 2.21E-03, 2.85E-03, 3.27E-03, 3.62E-03, 2.44E-03,   &
  1.63E-03, 9.90E-04, 7.56E-04, 5.98E-04, 4.96E-04, 5.08E-04,   &
  5.20E-04, 5.32E-04, 5.44E-04, 5.55E-04, 5.67E-04, 5.79E-04,   &
  5.90E-04, 6.02E-04, 6.14E-04, 6.26E-04, 6.38E-04, 6.49E-04,   &
  6.61E-04, 6.73E-04, 6.80E-04, 6.84E-04/

! C-->Si, physical
integer,parameter :: n_c_si = 20
real :: c_si_yields(n_c_si)
real :: c_si_e(n_c_si)
data c_si_e / &
     25,    30,    35,    40,    50,    60,    70,    80,    100,  &
    120,   140,   160,   200,   300,   400,   500,   800,   1600,  &
    3200, 6400/
data c_si_yields / &
  1.14E-02, 2.50E-02, 4.27E-02, 6.25E-02, 1.06E-01, 1.51E-01,   &
  1.93E-01, 2.34E-01, 3.07E-01, 3.71E-01, 4.27E-01, 4.77E-01,   &
  5.63E-01, 7.22E-01, 8.33E-01, 9.16E-01, 1.07E+00, 1.23E+00,   &
  1.27E+00, 1.18E+00/

! D-->C, physical
integer,parameter :: n_d_c = 23
real :: d_c_yields(n_d_c)
real :: d_c_e(n_d_c)
data d_c_e / &
     25,    30,    35,    40,    50,    60,    70,    80,    100, &
    120,   140,   160,   200,   300,   400,   500,   800,   1600, &
   3200,  6400, 12000, 24000, 48000/
data d_c_yields / &
  9.70E-06, 2.05E-04, 7.92E-04, 1.79E-03, 4.65E-03, 7.97E-03, &
  1.15E-02, 1.47E-02, 2.07E-02, 2.58E-02, 3.01E-02, 3.36E-02, &
  3.91E-02, 4.64E-02, 4.92E-02, 5.00E-02, 4.82E-02, 3.91E-02, &
  2.73E-02, 1.67E-02, 9.85E-03, 5.25E-03, 2.84E-03/

! D-->C, chemical
! Roth equation
integer,parameter :: n_d_c_ch = 28
real :: d_c_ch_yields(n_d_c_ch)
real :: d_c_ch_e(n_d_c_ch)
data d_c_ch_e / &
        5,    10,    15,    20,    25,    30,    35,    40,    50,  &
       60,    70,    80,   100,   120,   140,   160,   200,   300,  &
      400,   500,   600,   800,  1000,  1400,  1800,  2400,  3000,  &
     4000/
data d_c_ch_yields / &
  2.74E-03, 8.64E-03, 1.25E-02, 1.49E-02, 1.64E-02, 1.74E-02,    &
  1.79E-02, 1.81E-02, 1.77E-02, 1.67E-02, 1.53E-02, 1.36E-02,    &
  1.02E-02, 7.13E-03, 4.75E-03, 3.06E-03, 1.19E-03, 9.77E-05,    &
  7.71E-06, 6.04E-07, 4.74E-08, 2.93E-10, 1.82E-12, 7.22E-17,    &
  2.92E-21, 7.72E-28, 2.09E-34, 2.46E-45/

! C-->C, physical
integer,parameter :: n_c_c = 23
real :: c_c_yields(n_c_c)
real :: c_c_e(n_c_c)
data c_c_e /& 
     25,    30,    35,    40,    50,    60,    70,    80,    100, &
    120,   140,   160,   200,   300,   400,   500,   800,   1600, &
   3200,  6400, 12000, 24000, 48000/
data c_c_yields / &
  7.38E-04, 3.07E-03, 7.74E-03, 1.48E-02, 3.54E-02, 6.09E-02,    &
  8.89E-02, 1.17E-01, 1.71E-01, 2.19E-01, 2.62E-01, 3.00E-01,    &
  3.64E-01, 4.81E-01, 5.62E-01, 6.20E-01, 7.28E-01, 8.35E-01,    &
  8.55E-01, 7.83E-01, 6.57E-01, 4.87E-01, 3.30E-01/
  
! Si --> C, physical
integer,parameter :: n_si_c = 20
real :: si_c_yields(n_si_c)
real :: si_c_e(n_si_c)
data si_c_e /& 
     25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 160, 200, 300, &
     400, 500, 800, 1600, 3200, 6400/
data si_c_yields / &
  8.30E-06, 1.46E-04, 8.10E-04, 2.55E-03, 1.06E-02, 2.53E-02, &
  4.56E-02, 7.03E-02, 1.27E-01, 1.85E-01, 2.44E-01, 2.99E-01, &
  3.98E-01, 5.93E-01, 7.38E-01, 8.52E-01, 1.09E+00, 1.42E+00, &
  1.68E+00, 1.85E+00/
  
! Si --> Si, physical
integer,parameter :: n_si_si = 20
real :: si_si_yields(n_si_si)
real :: si_si_e(n_si_si)
data si_si_e /& 
     25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 160, 200, 300, &
     400, 500, 800, 1600, 3200, 6400/
data si_si_yields / &
  1.23E-02, 2.77E-02, 4.84E-02, 7.29E-02, 1.29E-01, 1.87E-01, &
  2.45E-01, 3.01E-01, 4.03E-01, 4.96E-01, 5.79E-01, 6.54E-01, &
  7.88E-01, 1.05E+00, 1.24E+00, 1.40E+00, 1.74E+00, 2.23E+00, &
  2.64E+00, 2.88E+00/
  
! Si --> SiC, C, physical
integer,parameter :: n_si_sic_c = 31
real :: si_sic_c_yields(n_si_sic_c)
real :: si_sic_c_e(n_si_sic_c)
data si_sic_c_e /& 
     15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 110, 120, 130, 140, &
     150, 160, 180, 200, 250, 300, 350, 400, 500, 600, 700, 800, 1200, &
     1600, 3200, 6400/
data si_sic_c_yields / &
  1.06E-04, 1.27E-03, 4.69E-03, 1.07E-02, 1.87E-02, 2.84E-02, &
  5.15E-02, 7.77E-02, 1.05E-01, 1.34E-01, 1.88E-01, 2.14E-01, &
  2.39E-01, 2.62E-01, 2.85E-01, 3.07E-01, 3.28E-01, 3.67E-01, &
  4.03E-01, 4.83E-01, 5.51E-01, 6.10E-01, 6.62E-01, 7.49E-01, &
  8.20E-01, 8.80E-01, 9.32E-01, 1.08E+00, 1.19E+00, 1.39E+00, &
  1.52E+00/
  
! Si --> SiC, Si, physical
integer,parameter :: n_si_sic_si = 31
real :: si_sic_si_yields(n_si_sic_si)
real :: si_sic_si_e(n_si_sic_si)
data si_sic_si_e /& 
     15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 110, 120, 130, 140, &
     150, 160, 180, 200, 250, 300, 350, 400, 500, 600, 700, 800, 1200, &
     1600, 3200, 6400/
data si_sic_si_yields / &
  6.20E-06, 1.77E-04, 1.08E-03, 3.40E-03, 7.66E-03, 1.40E-02, &
  3.19E-02, 5.51E-02, 8.06E-02, 1.07E-01, 1.59E-01, 1.85E-01, &
  2.09E-01, 2.32E-01, 2.55E-01, 2.76E-01, 2.96E-01, 3.35E-01, &
  3.70E-01, 4.48E-01, 5.14E-01, 5.70E-01, 6.18E-01, 7.01E-01, &
  7.67E-01, 8.23E-01, 8.70E-01, 1.01E+00, 1.10E+00, 1.27E+00, &
  1.37E+00/
  
public :: yield_mm

!    integer,parameter ::  maxdata=100
!    real :: yield_d(maxdata), energy_d(maxdata)
!    real :: yield_si(maxdata), energy_si(maxdata)
!    real :: yield_c(maxdata), energy_c(maxdata)

!    integer :: loaded_mat = 0
!    real :: loaded_opt = -1000.0

!    integer :: n_data_d = 0
!    integer :: n_data_si = 0
!    integer :: n_data_c = 0

!    save yield_d, energy_d, n_data_d
!    save yield_c, energy_c, n_data_c
!    save yield_si, energy_si, n_data_si

!    save loaded_opt, loaded_mat
                                                                                                  
contains

!  real function yld_sic_mm(energy, matp, matt)
  
!    implicit none
!    integer, intent(in) :: matp, matt
!    real :: energy
    
!    ! MATT = 20 = C (SiC)
!    ! MATT = 21 = Si (SiC)
    
!    ! MATP = 2    D
!    ! MATP = 6    Self
!    ! MATP = 5    C
!    ! MATP = 8    Si
    
!    integer in,ipos
!    external ipos
!    !real lin_interp
!    !external lin_interp

!    ! Initialize
!    yld_sic = 0.0
   
!    ! Load data if not already loaded
!    if (loaded_mat.ne.matt) then 
!       ! Load sputtering data
!       !call get_yields(matt,maxdata,n_data_self,energy_self,yield_self,&
!       !     n_data_d,energy_d,yield_d)

!       call get_yields(matt, maxdata, n_data_c, energy_c, yield_c, &
!         n_data_si, energy_si, yield_si, n_data_d, energy_d, yield_d)

!       ! Save loaded data option
!       !write(0,*) 'Loading:',matt,n_data_self,n_data_d,loaded_opt,opt,loaded_mat
!       loaded_mat = matt
!    endif
  
!  end function yld_sic_mm
  
  real function lin_interp(val, lookup, result_data, ndata, minopt)
  
    ! Linear interpolation at val using the point-slope formula. If val
    ! is greater than the last lookup values, the last result_data is
    ! returned.
    ! 
    ! val (real) : The X value to get interpolated value at.
    ! lookup (real, ndata) : X coordinates of data to be interpolated.
    ! result_data (real, ndata) : Y coordinates
    ! ndata (int) : Number of coordinates
    ! minopt (int) : Option for how to treat values less than the 
    !   minimum in lookup. 0 = return zero. 1 = return the first 
    !   result_data value.
  
    implicit none
    integer ndata
    integer, intent(in) ::  minopt
    real val, lookup(ndata), result_data(ndata)

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
  
!  subroutine get_yields(matt, maxdata, n_data_c, energy_c, yield_c, &
!    n_data_si, energy_si, yield_si, n_data_d, energy_d, yield_d)

!    implicit none
!    integer :: matt
!    integer :: n_data_d, maxdata, n_data_c, n_data_si
    
!    ! For yields of X = Si or C, we need the following yields.
!    ! D->SiC,X   D->SiC,X,ch   C->SiC,X
!    ! D->X       D->X,ch       C->X
    

!   ! Arrays containing energies and yields for each yield needed in the
!   ! mixed material model. LEFT OFF HERE.
!    real energy_d(maxdata),    yield_d(maxdata)
!    real energy_d_ch(maxdata), yield_d_ch(maxdata)
!    real energy_c(maxdata),    yield_c(maxdata)
!    real energy_si(maxdata),   yield_si(maxdata)
    
!    ! MATT = 20 is C (SiC)
!    ! MATT = 21 is Si (SiC)
    
!    !write(0,*) 'Get yield:',matt
    
!    ! Load C related yields.
!    if (matt.eq.20) then 

!       ! Load SiC, C yields by D, physical
!       call load_yield(n_data_d, maxdata, energy_d, yield_d, &
!         n_d_sic_c_e, d_sic_c_e, d_sic_c_yields)
         
!       ! Load C yields by D, chemical

!       ! Load SiC, C yields by C
!       call load_yield(n_data_c, maxdata, energy_c, yield_c, &
!         n_c_sic_c_e, c_sic_c_e, c_sic_c_yields)

!       ! Load SiC, C yields by Si
!       !call load_yield(n_data_si, maxdata, energy_si, yield_si, &
!       !  n_c_sic_si_e, c_sic_si_e, c_sic_si_yields)
         
       

!    ! Load Si (SiC) data
!    elseif (matt.eq.21) then 

!       ! Load Si yields by D
!       call load_yield(n_data_d, maxdata, energy_d, yield_d,&
!         n_si_sic_d_e, si_sic_d_e, si_sic_d_yields)

!       ! Load Si yields by C
!       call load_yield(n_data_c, maxdata, energy_c, yield_c, &
!         n_si_sic_c_e, si_sic_c_e, si_sic_c_yields)

!       ! Load Si yields by Si
!       call load_yield(n_data_si, maxdata, energy_si, yield_si, &
!         n_si_sic_si_e, si_sic_si_e, si_sic_si_yields)

!    endif

!    return
!  end subroutine get_yields
  
!  subroutine load_yield(ndata, maxdata, energy, yield, n_src_e, src_e, &
!    src_yields)
    
!    ! Loads energies and yields into energy and yield.
    
!    implicit none
!    integer :: ndata, maxdata
!    real    :: energy(maxdata), yield(maxdata)

!    integer :: n_src_e
!    real :: src_e(n_src_e)
!    real :: src_yields(n_src_e)

!    integer :: in, data_count

!    data_count = 0

!    do in = 1, n_src_e

!          if (src_yields(in).gt.0.0) then 
!             data_count = data_count + 1
!             energy(data_count) = src_e(in)
!             yield(data_count) = src_yields(in)
!             !write(0,'(a,i6,1x,f10.2,1x,g12.5)') 'Loading:', &
!             !  data_count,energy(data_count),yield(data_count)
!          endif
   
!    end do

!    ndata = data_count

!    return
!  end subroutine load_yield
  

!  subroutine print_sic_yields_mm(ounit)
    
    
!    ! Prints the yield data currently in use to the specified unit.

!   implicit none
!    integer :: ounit, in

!    write(ounit,'(a)') '     2020 SiC SPUTTERING DATA:'
!    write(ounit,'(a)') '     BASE SPUTTERING YIELDS BY D IN USE (MAY BE MULTIPLIED BY ANY SPECIFIED YIELD MULTIPLIERS)'
!    write(ounit,'(25x,a,12x,a)')  'E (eV)','YIELD'
!    do in = 1,n_data_d
!       write(ounit,'(10x,i6,3x,f12.2,3x,g12.5)')  in,energy_d(in),yield_d(in)
!    end do

!    !write(ounit,'(a)')
!    !write(ounit,'(a)') '     BASE SELF SPUTTERING YIELDS IN USE (MAY BE MULTIPLIED BY ANY SPECIFIED YIELD MULTIPLIERS)'
!    !write(ounit,'(25x,a,12x,a)')  'E (eV)','YIELD'
!    !do in = 1,n_data_self
!    !   write(ounit,'(10x,i6,3x,f12.2,3x,g12.5)')  in,energy_self(in),yield_self(in)
!    !end do

!    write(ounit,'(a)')
!    write(ounit,'(a)') '     BASE SPUTTERING YIELDS BY C IN USE (MAY BE MULTIPLIED BY ANY SPECIFIED YIELD MULTIPLIERS)'
!    write(ounit,'(25x,a,12x,a)')  'E (eV)','YIELD'
!    do in = 1,n_data_c
!       write(ounit,'(10x,i6,3x,f12.2,3x,g12.5)')  in,energy_c(in),yield_c(in)
!    end do

!    write(ounit,'(a)')
!    write(ounit,'(a)') '     BASE SPUTTERING YIELDS BY SI IN USE (MAY BE MULTIPLIED BY ANY SPECIFIED YIELD MULTIPLIERS)'
!    write(ounit,'(25x,a,12x,a)')  'E (eV)','YIELD'
!    do in = 1,n_data_si
!       write(ounit,'(10x,i6,3x,f12.2,3x,g12.5)')  in,energy_si(in),yield_si(in)
!    end do

!  end subroutine print_sic_yields_mm
  
  real function yield_mm(energy, cion, frac_c, frac_si, mm_usage)
  
    ! We've already defined the energies and yields above, so why use
    ! get_yields like the other options if it just means duplicating
    ! arrays in memory?
    
    implicit none
    integer, intent(in) :: cion
    real :: energy, frac_c, frac_si, yield_sic_si, yield_sic_c, yield_si
    real :: yield_c, surf_si, surf_c, refl_coef
    integer :: mm_usage
    
    yield_mm = 0.0
    
    ! Hardcoding in value from Abrams paper.
    refl_coef = 0.1
    
    ! What is this, this isn't correct usage of this?
    ! Additional usage to just pull out the chemical sputtering yields 
    ! for D-->C.
    !if (mm_usage.eq.2) then
    !  yield_c = lin_interp(energy, d_c_ch_e, d_c_ch_yields, n_d_c_ch, 1)
    !  yield_mm = yield_c
    !  return
    !endif
    
    ! Total yield of C from C.
    yield_c = &
      lin_interp(energy, d_c_e, d_c_yields, n_d_c, 1) + &
      lin_interp(energy, d_c_ch_e, d_c_ch_yields, n_d_c_ch, 1) + &
      frac_c * lin_interp(energy, c_c_e, c_c_yields, n_c_c, 1) + &
      frac_si * lin_interp(energy, si_c_e, si_c_yields, n_si_c, 1) 
      
    ! Usage switch allows using just the graphite yields to facilitate 
    ! comparisons between SiC and C.
    if (mm_usage.eq.1) then
      yield_mm = yield_c
      return
    endif
      
    ! Total yield of Si from Si.
    yield_si = &
      lin_interp(energy, d_si_e, d_si_yields, n_d_si, 1) + &
      lin_interp(energy, d_si_ch_e, d_si_ch_yields, n_d_si_ch, 1) + &
      frac_c * lin_interp(energy, c_si_e, c_si_yields, n_c_si, 1) + &
      frac_si * lin_interp(energy, si_si_e, si_si_yields, n_si_si, 1)
        
    ! Usage switch allows using just the silicon yields to facilitate 
    ! comparisons between SiC and Si.
    if (mm_usage.eq.2) then
      yield_mm = yield_si
      return
    endif

    ! Total yield of C from SiC.
    yield_sic_c = &
    lin_interp(energy, d_sic_c_e, d_sic_c_yields, &
      n_d_sic_c, 1) + &
    lin_interp(energy, d_sic_c_ch_e, d_sic_c_ch_yields, &
      n_d_sic_c_ch, 1) + &
    frac_c * lin_interp(energy, c_sic_c_e, c_sic_c_yields, &
      n_c_sic_c, 1) + &
    frac_si * lin_interp(energy, si_sic_c_e, si_sic_c_yields, &
      n_si_sic_c, 1)

    ! Total yield of Si from SiC. Chemical sputtering of Si from SiC 
    ! is zero.
    yield_sic_si = &
      lin_interp(energy, d_sic_si_e, d_sic_si_yields, &
        n_d_sic_si, 1) + &
      frac_c * lin_interp(energy, c_sic_si_e, c_sic_si_yields, &
        n_c_sic_si, 1) + &
      frac_si * lin_interp(energy, si_sic_si_e, si_sic_si_yields, &
        n_si_sic_si, 1)
            
    ! Carbon and silicon surface abundances. SiC abundance is 
    ! just (1 - surf_c - surf_sic). If yields are zero (e.g. call
    ! function with energy = 0.0), then set C and Si surface abundances
    ! to zero.
    if (yield_c.eq.0.0) then
      surf_c = 0.0
    else
      surf_c = (1 - refl_coef) * frac_c / yield_c
    endif
    
    ! At low temperatures, surf_c can go above 1, but this is physically
    ! unrealistic. So we cap surf_c at 1.
    if (surf_c.gt.1) then
      surf_c = 1.0
    endif
    
    if ((yield_si + yield_sic_c - yield_sic_si).eq.0.0) then
      surf_si = 0.0
    else
      surf_si = (1 - surf_c) * (yield_sic_c - yield_sic_si) / &
        (yield_si + yield_sic_c - yield_sic_si)
    endif
          
    ! Finally, depending on the ion being simulated, add together.
    ! sazmod - Incredibly, I had these swapped in error... shit
    if (cion.eq.6) then
      yield_mm = (1 - surf_c - surf_si) * yield_sic_c + &
        surf_c * yield_c
    elseif (cion.eq.14) then
      yield_mm = (1 - surf_c - surf_si) * yield_sic_si + &
        surf_si * yield_si
    else
      write(0,*) 'yield_mm error! cion = ', cion
    endif
    
    ! Debugging print statement.
    !write(0,*) 'SiC MM yield: ',yield_mm,yield_sic_si,yield_sic_c, &
    !  surf_c,surf_si
    
    return
  
  end function yield_mm


end module sic_mixed_material
