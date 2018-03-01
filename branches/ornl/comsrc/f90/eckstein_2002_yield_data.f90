module eckstein_2002_yield_data

use error_handling

implicit none

private


!
! Notes: 1) A value of -1 is inserted for data which is not present in the tabulated source data.
!        2) A value of 0 has been inserted in some cicumstances (where data is not available) in order
!           to facilitate the calculation of the angular weighted average yield. 
!        3) Although included for completeness - the angle-averaged yields are not considered usable unless
!           a minimum of 5 data points contributed to the calculation of the average value. 
!



! Beryllium data from Eckstein 2002
! D -> Be ebd=3.38eV - Eckstein IPP report 9/132 June 2002 Page 16

integer,parameter :: n_be_d_e = 21
integer,parameter :: n_be_d_angle = 13

real :: be_d_yields(n_be_d_angle,n_be_d_e)
real :: be_d_angle(n_be_d_angle)
real :: be_d_e(n_be_d_e)


data be_d_e /10, 11, 12, 13, 14, 15, 17, 20, 25, 30, 40, 50, 70, 100, 140, 200, 300, 500, 1000, 2000, 3000/
data be_d_angle /0, 15 , 30 , 45 , 55 , 60 , 65 , 70 , 75 , 80 , 85 , 87, -1/

data be_d_yields / &
4.13e-6, 4.30e-6, 4.56e-6,   -1   ,   -1   ,   -1   ,   -1    ,  -1    ,   -1   ,   -1   ,  -1     , -1     , 2.14493e-06,&
2.61e-5, 3.05e-5, 3.15e-5, 2.69e-5, 1.86e-5,   -1   , 1.00e-5 ,  -1    , 3.67e-6, 8.00e-7, 5.00e-7 , -1     , 2.51589e-05,&
9.29e-5, 1.07e-4, 1.25e-4, 1.14e-4, 8.71e-5, 6.81e-5, 5.37e-5 ,  -1    , 2.31e-5, 1.39e-5, 9.62e-6 , -1     , 9.86572e-05,&
2.32e-4, 2.67e-4, 3.47e-4, 3.16e-4, 2.61e-4, 2.10e-4, 1.62e-4 ,  -1    , 7.67e-5, 4.78e-5, 3.18e-5 , -1     , 0.000266474,&
4.65e-4, 5.50e-4, 6.90e-4, 6.91e-4, 5.92e-4,   -1   , 3.92e-4 ,  -1    , 1.91e-4, 1.20e-4, 7.85e-5 , -1     , 0.000559756,&
8.10e-4, 9.42e-4, 1.29e-3, 1.29e-3, 1.12e-3, 9.82e-4, 7.97e-4 ,  -1    , 4.11e-4, 2.57e-4, 1.65e-4 , -1     , 0.0010252  ,&
1.73e-3, 2.02e-3, 2.60e-3, 3.03e-3, 2.92e-3,   -1   , 2.29e-3 ,  -1    , 1.25e-3, 7.67e-4, 4.62e-4 , -1     , 0.00231114 ,&
3.64e-3, 4.19e-3, 5.46e-3, 7.03e-3, 7.46e-3, 7.12e-3, 6.59e-3 ,  -1    , 3.76e-3, 2.14e-3, 1.16e-3 , -1     , 0.00525801 ,&
7.30e-3, 8.25e-3, 1.11e-2, 1.57e-2, 1.87e-2, 1.94e-2, 1.85e-2 ,  -1    , 1.08e-2, 5.51e-3, 2.51e-3 , -1     , 0.0117133  ,&
1.08e-2, 1.22e-2, 1.69e-2, 2.54e-2, 3.25e-2, 3.58e-2, 3.46e-2 ,  -1    , 2.07e-2, 1.02e-2, 3.82e-3 , -1     , 0.0190356  ,&
1.68e-2, 1.90e-2, 2.64e-2, 4.34e-2, 6.08e-2,   -1   , 7.10e-2 ,  -1    , 4.54e-2, 2.12e-2, 5.91e-3 , -1     , 0.03289    ,&
2.09e-2, 2.40e-2, 3.45e-2, 5.78e-2, 6.53e-2, 9.88e-2, 1.07e-1 ,  -1    , 7.40e-2,   -1   ,   -1    , -1     , 0.042348   ,&
  -1   ,   -1   ,   -1   ,   -1   ,   -1   ,   -1   , 1.60e-1 ,  -1    , 1.28e-1, 5.00e-2, 9.51e-3 , -1     , 0.0110313  ,&
3.10e-2, 3.59e-2, 5.42e-2, 9.78e-2, 1.49e-1, 1.80e-1, 2.13e-1 ,  -1    , 2.00e-1, 1.01e-1, 1.27e-2 , -1     , 0.079691   ,&
3.32e-2, 3.94e-2, 6.11e-2, 1.08e-1, 1.68e-1,   -1   , 2.48e-1 ,  -1    , 2.71e-1, 1.60e-1, 1.82e-2 , -1     , 0.0916089  ,&
3.51e-2, 4.03e-2, 6.32e-2, 1.14e-1, 1.77e-1, 2.20e-1, 2.68e-1 ,  -1    , 3.42e-1, 2.43e-1, 2.93e-2 , -1     , 0.0997961  ,&
3.44e-2, 4.20e-2, 6.34e-2, 1.12e-1, 1.75e-1, 2.12e-1, 2.72e-1 ,  -1    , 3.91e-1, 3.42e-1, 5.71e-2 , -1     , 0.102859   ,&
3.24e-2, 3.85e-2, 5.82e-2, 1.04e-1, 1.50e-1, 1.98e-1, 2.57e-1 ,  -1    , 4.08e-1, 4.39e-1, 1.38e-1 , -1     , 0.0983351  ,&
2.53e-2, 2.96e-2, 4.37e-2, 7.80e-2, 1.25e-1, 1.54e-1, 2.06e-1 ,  -1    , 3.65e-1, 4.59e-1, 3.34e-1 , -1     , 0.08083    ,&
1.76e-2,   -1   ,   -1   ,   -1   ,   -1   ,   -1   ,  -1     ,  -1    ,   -1   ,   -1   ,    -1   , -1     , 0.0176     ,&
1.25e-2,   -1   , 2.02e-2, 3.43e-2,   -1   , 7.29e-2,  -1     , 1.29e-1,   -1   , 2.89e-1, 4.59e-1 , 3.18e-1, 0.0414479   &
/  


! Be ->Be ebd=3.38eV - Eckstein IPP report 9/132 June 2002 Page 25


integer,parameter :: n_be_self_e = 26
integer,parameter :: n_be_self_angle = 10


real :: be_self_yields(n_be_self_angle,n_be_self_e)
real :: be_self_angle(n_be_self_angle)
real :: be_self_e(n_be_self_e)


data be_self_e /5, 6, 7, 8, 9, 10, 11, 12 , 13, 15, 17, 20, 25, 30, 40, 50, 70, 100, 200, 300, 500, 700, 1000, 2000, 3000, 5000/   
data be_self_angle /0 , 15, 30, 45, 60, 65, 75, 80, 85, -1/

data be_self_yields /& 
  0    ,    0   ,       0   ,   3.05e-6,   7.10e-6,   9.83e-6,   1.05e-5,   1.02e-5,   1.03e-5,  1.99517e-06,&     
  0    ,    0   ,    3.17e-6,   2.32e-5,   3.87e-5,   4.62e-5,   4.62e-5,   4.64e-5,   4.61e-5,  1.30308e-05,&     
  0    ,    0   ,    1.36e-5,   6.54e-5,   9.76e-5,   1.18e-4,   1.30e-4,   1.34e-4,   1.38e-4,  3.56533e-05,&     
  0    , 3.35e-6,    3.53e-5,   1.35e-4,   1.96e-4,   2.71e-4,   3.47e-4,   3.80e-4,   4.01e-4,  8.45538e-05,&     
  0    , 9.66e-6,    7.48e-5,   2.48e-4,   3.92e-4,   6.31e-4,   8.84e-4,   9.58e-4,   1.01e-3,  0.000181513,&      
  0    , 2.12e-5,    1.33e-4,   4.22e-4,   7.84e-4,   1.37e-3,   1.89e-3,   2.04e-3,   2.17e-3,  0.000358359,&      
6.78e-6, 4.03e-5,    2.21e-4,   7.39e-4,   1.51e-3,   2.62e-3,   3.47e-3,   3.77e-3,   3.74e-3,  0.000663614,&      
1.26e-5, 6.63e-5,    3.47e-4,   1.22e-3,   2.74e-3,   4.67e-3,   5.88e-3,   6.00e-3,   6.04e-3,  0.00113349 ,&      
2.14e-5, 1.01e-4,    5.39e-4,   2.07e-3,   4.48e-3,   7.29e-3,   8.69e-3,   8.73e-3,   8.84e-3,  0.00178831 ,&      
5.23e-5, 2.17e-4,    1.22e-3,   4.95e-3,   9.90e-3,   1.48e-2,   1.61e-2,   1.57e-2,   1.51e-2,  0.00377815 ,&      
1.15e-4, 4.40e-4,    2.39e-3,   9.79e-3,   1.79e-2,   2.48e-2,   2.50e-2,   2.35e-2,   2.22e-2,  0.00670206 ,&      
3.05e-4, 1.09e-3,    5.63e-3,   2.08e-2,   3.47e-2,   4.38e-2,   4.04e-2,   3.65e-2,   3.33e-2,  0.0129422  ,&      
1.09e-3, 3.32e-3,    1.54e-2,   4.77e-2,   7.14e-2,   8.21e-2,   6.87e-2,   5.85e-2,   5.09e-2,  0.0274842  ,&      
2.68e-3, 7.27e-3,    2.95e-2,   8.18e-2,   1.16e-1,   1.25e-1,   9.74e-2,   7.92e-2,   6.58e-2,  0.0458087  ,&      
8.41e-3, 1.93e-2,    6.54e-2,   1.58e-1,   2.12e-1,   2.15e-1,   1.50e-1,   1.12e-1,   8.65e-2,  0.0875607  ,&      
1.68e-2, 3.45e-2,    1.03e-1,   2.33e-1,   3.06e-1,   3.04e-1,   1.99e-1,   1.37e-1,   9.65e-2,  0.130115   ,&      
3.77e-2, 6.72e-2,    1.73e-1,   3.59e-1,   4.70e-1,   4.71e-1,   2.87e-1,   1.73e-1,   1.03e-1,  0.207418   ,&      
7.00e-2, 1.11e-1,    2.49e-1,   4.93e-1,   6.63e-1,   6.94e-1,   4.21e-1,   2.29e-1,   1.06e-1,  0.299338   ,&      
1.43e-1, 1.98e-1,    3.77e-1,   7.18e-1,   1.02e-0,   1.21e-0,   8.77e-1,   4.38e-1,   1.16e-1,  0.479921   ,&      
1.86e-1, 2.46e-1,    4.37e-1,   8.19e-1,   1.19e-0,   1.52e-0,   1.30e-0,   6.90e-1,   1.37e-1,  0.580202   ,&      
2.33e-1, 2.93e-1,    4.94e-1,   9.03e-1,   1.34e-0,   1.85e-0,   1.93e-0,   1.23e-0,   2.08e-1,  0.687856   ,&      
2.57e-1, 3.16e-1,    5.11e-1,   9.30e-1,   1.40e-0,   2.00e-0,   2.31e-0,   1.72e-0,   3.09e-1,  0.739884   ,&      
2.74e-1, 3.31e-1,    5.19e-1,   9.31e-1,   1.41e-0,   2.09e-0,   2.65e-0,   2.27e-0,   5.10e-1,  0.774371   ,&      
2.53e-1,   -1   ,      -1   ,     -1   ,     -1   ,     -1   ,     -1   ,     -1   ,     -1   ,  0.0253     ,&      
2.63e-1,   -1   ,      -1   ,     -1   ,     -1   ,     -1   ,   2.99e-0,   3.38e-0,   2.20e-0,  0.771981   ,&      
2.27e-1,   -1   ,      -1   ,     -1   ,     -1   ,     -1   ,     -1   ,     -1   ,     -1   ,  0.0227      &      
/

 

! Carbon data from Eckstein 2002 
! D -> C ebd=7.41eV - Eckstein IPP report 9/132 June 2002 Page 44

integer,parameter :: n_c_d_e = 12
integer,parameter :: n_c_d_angle = 10

real :: c_d_yields(n_c_d_angle,n_c_d_e)
real :: c_d_angle(n_c_d_angle)
real :: c_d_e(n_c_d_e)

data c_d_e /30, 33, 40, 50, 70, 100, 140, 200, 300, 500, 1000, 2000/ 

data c_d_angle/0, 15, 30, 45, 55, 65, 75, 80, 85, -1/

data c_d_yields/&
8.58e-5, 1.12e-4, 1.78e-4, 2.30e-4, 2.08e-4, 1.12e-4, 2.55e-5, 8.60e-6, 2.56e-6, 0.00014548 ,&  
2.16e-4,   -1   ,   -1   ,   -1   ,   -1   , 3.29e-4,   -1   ,   -1   ,    -1  , 0.000201391,&  
7.35e-4, 9.15e-4, 1.26e-3, 1.72e-3, 2.21e-3, 1.80e-3, 7.60e-4, 3.12e-4, 9.35e-5, 0.00126345 ,&  
1.96e-3, 2.33e-3, 3.09e-3, 4.58e-3, 5.87e-3, 5.95e-3, 4.15e-3, 1.76e-3, 4.45e-4, 0.00344837 ,&  
4.79e-3, 5.19e-3, 7.01e-3, 1.13e-2, 1.57e-2, 2.14e-2, 1.78e-2, 8.60e-3, 1.57e-3, 0.00916153 ,&  
8.18e-3, 8.77e-3, 1.21e-2, 1.96e-2, 3.08e-2, 4.59e-2, 4.78e-2, 2.72e-2, 4.04e-3, 0.0175714  ,&  
1.10e-2, 1.21e-2, 1.68e-2, 2.88e-2, 4.40e-2, 7.03e-2, 8.63e-2, 5.73e-2, 8.17e-3, 0.0261649  ,&  
1.32e-2, 1.48e-2, 2.41e-2, 3.56e-2, 5.93e-2, 9.10e-2, 1.21e-1, 9.76e-2, 1.41e-2, 0.0347071  ,&  
1.47e-2, 1.66e-2, 2.40e-2, 4.33e-2, 6.81e-2, 1.13e-1, 1.63e-1, 1.51e-1, 2.94e-2, 0.0413178  ,&  
1.44e-2, 1.72e-2, 2.72e-2, 4.58e-2, 7.39e-2, 1.20e-1, 1.89e-1, 2.12e-1, 7.29e-2, 0.0456091  ,&  
1.30e-2, 1.49e-2, 2.36e-2, 4.18e-2, 6.49e-2, 1.08e-1, 1.90e-1, 2.35e-1, 1.72e-1, 0.0423841  ,&  
1.02e-2,   -1   ,   -1   , 3.45e-2,  -1    ,    -1  , 1.88e-1, 2.22e-1, 2.45e-1, 0.0391297   &  
/                                                                                                  
                                                                                                     
! C -> C ebd=7.41eV - Eckstein IPP report 9/132 June 2002 Page 53                                    
                                                                                                     
integer,parameter :: n_c_self_e = 19                                                                 
integer,parameter :: n_c_self_angle = 10                                                              
                                                                                                     
real :: c_self_yields(n_c_self_angle,n_c_self_e)                                                           
real :: c_self_angle(n_c_self_angle)                                                                    
real :: c_self_e(n_c_self_e)



data c_self_e /8, 10, 12, 15, 17, 20, 25, 30, 40, 45, 50, 70, 100, 140, 200, 300, 500, 1000, 1200/ 
data c_self_angle /0, 15, 30, 45, 55, 65, 75, 80, 85, -1/

data c_self_yields /&
  0    ,    0   ,    0   ,    0   ,   0    , 5.00e-6, 7.69e-6, 8.99e-6, 9.65e-6, 6.17863e-07,&      
  0    ,    0   , 5.00e-7, 6.50e-6, 2.07e-5, 2.60e-5, 3.81e-5, 3.87e-5, 3.78e-5, 7.07043e-06,&      
  0    ,    0   , 3.21e-6, 2.36e-5, 5.03e-5, 7.39e-5, 8.42e-5, 8.64e-5, 8.58e-5, 1.89791e-05,&      
  0    ,    0   , 1.10e-5, 8.90e-5, 1.47e-4, 2.00e-4, 2.36e-4, 2.60e-4, 2.64e-4, 5.73842e-05,&      
 -1    ,   -1   ,   -1   ,   -1   ,   -1   , 3.95e-4,   -1   ,   -1   ,    -1  , 0.00016693 ,&       
  0    , 1.00e-5, 8.10e-5, 3.60e-4, 6.74e-4, 1.18e-3, 1.73e-3, 1.95e-3, 2.05e-3, 0.000325632,&       
5.00e-6, 3.10e-5, 2.76e-4, 1.28e-3, 2.90e-3, 5.02e-3, 6.61e-3, 7.13e-3, 7.12e-3, 0.00128819 ,&       
1.83e-5, 1.00e-4, 8.15e-4, 3.95e-3, 8.02e-3, 1.31e-2, 1.55e-2, 1.55e-2, 1.60e-2, 0.00341181 ,&       
1.35e-4, 5.96e-4, 3.73e-3, 1.56e-2, 2.85e-2, 4.06e-2, 4.16e-2, 3.40e-2, 3.85e-2, 0.0113394  ,&       
2.74e-4,    -1  ,   -1   ,   -1   ,   -1   ,   -1   ,   -1   ,   -1   ,   -1   , 0.000274   ,&       
5.21e-4, 1.89e-3, 1.04e-2, 3.58e-2, 5.81e-2, 7.50e-2, 7.57e-2, 7.03e-2, 6.37e-2, 0.0236474  ,&       
2.57e-3, 7.31e-3, 3.01e-2, 8.40e-2, 1.34e-1, 1.63e-1, 1.51e-1, 1.28e-1, 1.10e-1, 0.054879   ,&       
8.84e-3, 1.96e-2, 6.53e-2, 1.66e-1, 2.45e-1, 2.89e-1, 2.45e-1, 1.99e-1, 1.54e-1, 0.104272   ,&       
2.13e-2, 3.95e-2, 1.11e-1, 2.52e-1, 3.74e-1, 4.43e-1, 3.65e-1, 2.60e-1, 1.83e-1, 0.163937   ,&       
4.14e-2, 6.76e-2, 1.63e-1, 3.49e-1, 5.15e-1, 6.33e-1, 5.14e-1, 3.53e-1, 2.11e-1, 0.235121   ,&       
7.16e-2, 1.06e-1, 2.27e-1, 4.60e-1, 6.83e-1, 8.63e-1, 7.42e-1, 4.94e-1, 2.42e-1, 0.324181   ,&       
1.16e-1, 1.61e-1, 3.05e-1, 5.98e-1, 8.91e-1, 1.18e-0, 1.15e-0, 7.92e-1, 3.01e-1, 0.445642   ,&       
1.78e-1, 2.28e-1, 3.93e-1, 7.38e-1, 1.10e-0, 1.55e-0, 1.77e-0, 1.45e-0, 4.92e-1, 0.592167   ,&       
2.13e-1,   -1   ,   -1   ,   -1   ,   -1   ,   -1   ,   -1   ,   -1   ,   -1   , 0.213       &      
/                                                                                                    


public :: yld2002, print_eck2002_yields


    integer,parameter ::  maxdata=100
    real :: yield_self(maxdata),energy_self(maxdata)
    real :: yield_d(maxdata),energy_d(maxdata)

    real ::loaded_opt = -1000.0
    integer :: loaded_mat = 0

    integer :: n_data_self = 0
    integer :: n_data_d = 0


    save yield_self,energy_self,n_data_self
    save yield_d,energy_d,n_data_d
    save loaded_opt,loaded_mat

                                                                                                    

contains


  real function yld2002(energy,matp,matt,opt)
    implicit none
    integer,intent(in) :: matp,matt
    real energy
    real,intent(in) :: opt

    ! Opt specifies the requested incident angle : -1 for average.
    !
    ! MATT = 2 = Be
    ! MATT = 4 = C
    !
    ! MATP = 2    D
    ! MATP = 6    Self
    !

    integer in,ipos
    external ipos
    !real lin_interp
    !external lin_interp


    !
    ! Initialize
    !
    yld2002 = 0.0
    !
    ! Load data if not already loaded
    !
    if (loaded_opt.ne.opt.or.loaded_mat.ne.matt) then 
       ! Load sputtering data
       call get_yields(matt,maxdata,n_data_self,energy_self,yield_self,&
            n_data_d,energy_d,yield_d,opt)
       ! Save loaded data option
       !write(0,*) 'Loading:',matt,n_data_self,n_data_d,loaded_opt,opt,loaded_mat
       loaded_mat = matt
       loaded_opt = opt
    endif


    if (matp.eq.2) then 
       !
       ! D Sputtering
       !
       yld2002 = lin_interp(energy,energy_d,yield_d,n_data_d,1)
    elseif (matp.eq.6) then 
       !
       ! Self sputtering
       !
       yld2002 = lin_interp(energy,energy_self,yield_self,n_data_self,1)
    endif


    return
  end function yld2002




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



  subroutine get_yields(matt,maxdata,n_data_self,energy_self,yield_self,&
       n_data_d,energy_d,yield_d,opt)
    implicit none
    integer :: matt
    integer :: n_data_d,n_data_self,maxdata
    real opt
    real energy_d(maxdata),yield_d(maxdata)
    real energy_self(maxdata),yield_self(maxdata)

    !

    ! Opt specifies incident angle : use -1 for average : data must exist or code stops
    !
    ! Yield data is taken from Ref: W Eckstein, IPP 9/132 June 2002.
    !
    ! MATT = 2    is Be
    ! MATT = 4    is C
    !
    !write(0,*) 'Get yield:',matt
    ! Load Be data 
    if (matt.eq.2) then 

       ! load self-sputtering data

       !write(0,*) 'SELF:',n_be_self_angle,n_be_self_e
       ! IPP/08 Krieger - SUN compiler hates more than 132 columns
       call load_yield(n_data_self,maxdata,energy_self,yield_self,opt,n_be_self_e,be_self_e,n_be_self_angle,be_self_angle,&
                    &  be_self_yields)


       ! load deuterium yields
       !write(0,*) 'D:',n_be_d_angle,n_be_d_e
       call load_yield(n_data_d,maxdata,energy_d,yield_d,opt,n_be_d_e,be_d_e,n_be_d_angle,be_d_angle,be_d_yields)

       ! Load Carbon data
    elseif (matt.eq.4) then 


       ! load self-sputtering data

       call load_yield(n_data_self,maxdata,energy_self,yield_self,opt,n_c_self_e,c_self_e,n_c_self_angle,c_self_angle,c_self_yields)


       ! load deuterium yields
       call load_yield(n_data_d,maxdata,energy_d,yield_d,opt,n_c_d_e,c_d_e,n_c_d_angle,c_d_angle,c_d_yields)

    endif

    return
  end subroutine get_yields


  subroutine load_yield(ndata,maxdata,energy,yield,incident_angle,n_src_e,src_e,n_src_angle,src_angle,src_yields)
    implicit none
    integer :: ndata,maxdata
    real    :: energy(maxdata), yield(maxdata)

    real :: incident_angle

    integer :: n_src_e, n_src_angle
    real :: src_e(n_src_e)
    real :: src_angle(n_src_angle)
    real :: src_yields(n_src_angle,n_src_e)


    integer :: ia,in,data_count

    ! load data for incident_angle from src_yields into yield array

    ! First find if data for the specified incident_angle exists

    ia = -1

    do in = 1,n_src_angle

       if (src_angle(in).eq.incident_angle) then 
          ! found appropriate angle data
          ia = in
          exit
       endif

    end do

    !write(0,'(a,3i5,3(1x,g12.5))') 'Loading yields:',ia,n_src_e,n_src_angle,incident_angle,src_angle(ia)


    ! Matching data not found - stop
    if (ia.eq.-1) then
       call errmsg('ERROR IN LOAD_YIELD - SPUTTER DATA FOR SPECIFIED INCIDENT ANGLE NOT FOUND:',incident_angle)
       stop
    endif

    !
    ! Load requested data - leave out -1 entries - if average data is requested check for at least 5 good data values in average. 
    !

    data_count = 0

    do in = 1,n_src_e


       if (incident_angle.ge.0.0) then 
          ! specific angle

          if (src_yields(ia,in).gt.0.0) then 
             data_count = data_count + 1
             energy(data_count) = src_e(in)
             yield(data_count) = src_yields(ia,in)
             !write(0,'(a,i6,1x,f10.2,1x,g12.5)') 'Loading:',data_count,energy(data_count),yield(data_count)
          endif

       else
          ! average data
          if (src_yields(ia,in).gt.0.0.and.check_valid_average(src_yields,n_src_angle,n_src_e,in)) then 
             data_count = data_count + 1   
             energy(data_count) = src_e(in)
             yield(data_count) = src_yields(ia,in)
             !write(0,'(a,i6,1x,f10.2,1x,g12.5)') 'Loading AV:',data_count,energy(data_count),yield(data_count)
          endif

       endif
    end do

    ndata = data_count

    return
  end subroutine load_yield



  logical function check_valid_average(src_yields,n_src_a,n_src_e,ndata)
    implicit none
    integer :: n_src_e,n_src_a,ndata
    real :: src_yields(n_src_a,n_src_e)

    integer :: ia, no_data_count 

    no_data_count=0

    do ia = 1, n_src_a

       if (src_yields(ia,ndata).lt.0.0) then 
          no_data_count=no_data_count+1
       endif
    end do


    ! need at least 5 data points to consider it a good data point - note that .gt. is used since average is one of the
    ! values included in the n_src_a total. 

    if ((n_src_a - no_data_count).gt.5) then 
       check_valid_average = .true.
    else
       check_valid_average = .false.
    endif

    return
  end function check_valid_average


  subroutine print_eck2002_yields(ounit)
    implicit none
    !
    ! Prints the yield data currently in use to the specified unit
    !
    integer :: ounit

    integer :: in

    write(ounit,'(a)') '     2002 ECKSTEIN SPUTTERING DATA:'
    write(ounit,'(a)') '     BASE D SPUTTERING YIELDS IN USE (MAY BE MULTIPLIED BY ANY SPECIFIED YIELD MULTIPLIERS)'
    write(ounit,'(25x,a,12x,a)')  'E (eV)','YIELD'
    do in = 1,n_data_d
       write(ounit,'(10x,i6,3x,f12.2,3x,g12.5)')  in,energy_d(in),yield_d(in)
    end do

    write(ounit,'(a)')
    write(ounit,'(a)') '     BASE SELF SPUTTERING YIELDS IN USE (MAY BE MULTIPLIED BY ANY SPECIFIED YIELD MULTIPLIERS)'
    write(ounit,'(25x,a,12x,a)')  'E (eV)','YIELD'
    do in = 1,n_data_self
       write(ounit,'(10x,i6,3x,f12.2,3x,g12.5)')  in,energy_self(in),yield_self(in)
    end do




  end subroutine print_eck2002_yields



end module eckstein_2002_yield_data
