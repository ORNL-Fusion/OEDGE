module eckstein_2007_yield_data

  use error_handling

  implicit none

  private


  !

  ! Fitting Formula and Parameters used here are from:
  ! 
  ! Topics in Applied Physics 110
  !
  ! Behrisch and Eckstein (eds)
  !
  ! "Sputtering by Particle Bombardment, Experiments and Computer Calculations from Threshold to MeV Energies"
  ! 
  ! Chapter "Sputtering Yields" by W. Eckstein, Springer 2007, pgs 33 - 186
  !
  !
  ! Fitting Formula:
  !
  ! Y(E) = q S_n_KrC(El) [ E0/Eth -1 ]**mu / (  lamda/w(El) + (E0/Eth -1)**mu) 
  !
  ! S_n_KrC =  0.5 ln (1+1.2288 El) / w(El)
  !
  ! w(El) = El + 0.1728 * sqrt(El) + 0.008 * EL**0.1504
  !
  ! El = E0 (M2/(M1+M2)) * aL / (Z1*Z2*e**2) = E0 / Epsilon
  !
  ! aL = (9 * PI **2 / 128 )**(1/3) * aB (Z1**(2/3) + Z2**(2/3))**(-1/2) 
  !
  ! aB=0.0529177nm
  !  
  ! Z1, M1 = projectille atom
  ! Z2, M2 = target atom
  !
  ! q, lambda, mu = fitting paramters 
  !
  ! Epsilon is dependent only on the colliding materials and can be pre-calculated
  ! 
  ! 
  !
  ! The following is tabulated data from the report for the most common fusion materials. 
  ! Be, C, W and Mo 
  !
  ! Self-sputtering data has been moved to be the last entry in each table to facilitate 
  ! code use of the tables. 
  !
  !
  ! Table 1: W. Eckstein - pg 118
  !
  ! Be data
  !
  ! Ion   Target    Lambda       q       mu    Eth (eV)    Epsilon     Esb (eV)  Esb/Gamma
  !
  !  H     Be       0.8007    0.0564   1.5147   14.340    2.56510e2     3.38      9.32
  !  D     Be       1.7575    0.1044   1,9906   9.5059    2.82110e2     3.38      5.67
  !  T     Be       2.0794    0.1379   1.5660   9.4345    3.07966e2     3.38      4.49
  ! 3He    Be       0.7725    0.3310   1.6036  12.8963    6.65344e2     3.38      4.49
  ! 4He    Be       1.4745    0.3193   1.6989  12.3288    7.19545e2     3.38      3.97
  !  N     Be       5.2833    0.9334   2.5368  16.5425    5.46566e3     3.38      3.55 
  !  O     Be       1.2209    1.2024   1.6881  22.6648    6.97104e3     3.38      3.67
  ! Ne     Be       2.5474    1.8309   1.9400  22.7750    1.06588e4     3.38      3.96
  ! Ar     Be       0.8082    3.2032   1.5058  37.1816    3.68450e4     3.38      5.63 
  ! Kr     Be       0.3844    5.3588   1.9600  61.452     1.67028e5     3.38      9.64
  ! Xe     Be       0.4779    8.1740   1.8350  86.942     4.23834e5     3.38     14.06
  !  Be    Be       2.0334    0.8241   1.3437  16.9689    2.20796e3     3.38      3.38
  !
  ! C data
  !
  ! Ion   Target    Lambda       q       mu    Eth (eV)    Epsilon     Esb (eV)  Esb/Gamma
  !
  !  H     C        1.3533    0.0241   1.4103   38.630    4.14659e2     7.41     25.89
  !  D     C        1.2848    0.0539   1.1977   27.770    4.46507e2     7.41     15.08 
  !  T     C        1.9050    0.0718   1.1512   23.617    4.78673e2     7.41     11.54
  ! 3He    C        0.7341    0.2058   1.1956   29.883    1.02061e3     7.41     11.54
  ! 4He    C        4.5910    0.1951   1.7852   19.124    1.08716e3     7.41      9.88
  !  N     C        5.4288    0.7481   1.7701  34.9372    7.37899e3     7.41      7.45
  !  O     C        9.6110    1.0171   2.0102  34.1293    9.29758e3     7.41      7.56
  ! Ne     C        2.5015    1.1912   1.6551  46.6904    1.39308e4     7.41      7.92 
  ! Ar     C        1.2622    2.4576   1.3952  68.8460    4.57989e4     7.41     10.42 
  ! Kr     C        1.3628    3.4372   2.2366  88.2918    1.99609e5     7.41     16.90
  ! Xe     C        0.4408    4.3004   1.7734  145.4236   4.98349e5     7.41     24.13
  !  C     C       13.9666    0.7015   2.0947  21.4457    5.68684e3     7.41      7.41
  !
  !
  ! Table 5: W. Eckstein - pg 122
  !
  ! Mo Data
  !
  ! Ion   Target    Lambda       q       mu    Eth (eV)    Epsilon     Esb (eV)  Esb/Gamma
  !           
  !  H     Mo       0.5124    0.0114   1,1469  201.4886   4.71832e3     6.83      165.63 
  !  D     Mo       0.3241    0.0326   1.5410   97.7738   4.76698e3     6.83       84.95
  !  T     Mo       0.5078    0.0661   1.5955   67.1475   4.81614e3     6.83       57.71
  ! 3He    Mo       0.3541    0.1373   0.9926   75.3995   9.84614e3     6.83       57.71
  ! 4He    Mo       0.1537    0.1563   0.9989   59.3088   9.94365e3     6.83       44.44
  !  N     Mo       0.1157    1.7900   1.8032   28.2561   4.10879e4     6.83       15.36 
  !  O     Mo       0.1762    2.1069   2.4821   23.9169   4.83320e4     6.83       13.94
  ! Ne     Mo       0.2205    2.8995   2.6514   23.6170   6.38956e4     6.83       11.89
  ! Ar     Mo       0.1339    6.3606   1.9562   28.2149   1.43274e5     6.83        8.23
  ! Kr     Mo       0.1412   11.9419   1.6911   38.6337   4.17411e5     6.83        6.86
  ! Xe     Mo       0.2401   32.5719   1.6694   47.4030   8.47848e5     6.83        7.00
  ! Mo     Mo       0.3580   12.2715   2.1844   31.4737   5.33049e5     6.83        6.83
  !
  !
  ! Table 8: W. Eckstein - pg 125
  !
  ! W Data
  !
  ! Ion   Target    Lambda       q       mu    Eth (eV)    Epsilon     Esb (eV)  Esb/Gamma
  !           
  !  H     W        1.0087    0.0075   1.2046   457.42    9.86986e3     8.68      399.36  
  !  D     W        0.3583    0.0183   1.4410   228.84    9.92326e3     8.68      202.85 
  !  T     W        0.2870    0.0419   1.5802   153.8842  9.97718e3     8.68      136.48 
  ! 3He    W        0.2424    0.0884   1.2439   164.3474  2.02666e4     8.68      136.48
  ! 4He    W        0.1692    0.1151   1.7121   120.56    2.03728e4     8.68      104.13 
  !  N     W        0.0921    1.4389   2.0225   45.3362   7.90505e4     8.68       32.98
  !  O     W        0.0777    1.8824   1.7536   44.2135   9.19794e4     8.68       29.46
  ! Ne     W        0.0828    2.5520   1.9534   38.6389   1.19107e5     8.68       24.35
  ! Ar     W        0.2113    5.9479   2.3857   27.0503   2.46646e5     8.68       14.80
  ! Kr     W        0.1747   13.6917   2.5161   34.7592   6.36677e5     8.68       10.09
  ! Xe     W        0.1385   20.5321   2.0952   44.8701   1.18932e6     8.68        8.93
  !  W     W        2.2697   18.6006   3.1273   24.9885   1.99860e6     8.68        8.68
  !
  !


  !
  ! Set up named constants to reference the data and make the code easier to read
  !

  integer,parameter :: targ_be=1,targ_c=2,targ_mo=3,targ_w=4
  integer,parameter :: ion_h=1,ion_d=2,ion_t=3,ion_he3=4,ion_he4=5,ion_n=6,ion_o=7,ion_ne=8,&
       ion_ar=9,ion_kr=10,ion_xe=11,ion_self=12
  integer,parameter :: lambda_val=1,q_val=2,mu_val=3,eth_val=4,epsilon_val=5,esb_val=6

  !
  ! Set up the storage array for the parameter data 
  !

  integer,parameter :: n_target=4   ! 4 target materials - (Be, C, Mo, W)
  integer,parameter :: n_ion=12     ! 11 bombarding species + self-sputtering
  integer,parameter :: n_params=7   ! 7 tabulated parameters - 5 are used in the fitting formula

  real :: yield_parameters(n_params,n_ion,n_target)

  ! Target material 1: Be data
  ! Ion   Target    Lambda       q       mu    Eth (eV)    Epsilon     Esb (eV)  Esb/Gamma
  !
  !  H     Be       0.8007    0.0564   1.5147   14.340    2.56510e2     3.38      9.32
  !  D     Be       1.7575    0.1044   1,9906   9.5059    2.82110e2     3.38      5.67
  !  T     Be       2.0794    0.1379   1.5660   9.4345    3.07966e2     3.38      4.49
  ! 3He    Be       0.7725    0.3310   1.6036  12.8963    6.65344e2     3.38      4.49
  ! 4He    Be       1.4745    0.3193   1.6989  12.3288    7.19545e2     3.38      3.97
  !  N     Be       5.2833    0.9334   2.5368  16.5425    5.46566e3     3.38      3.55 
  !  O     Be       1.2209    1.2024   1.6881  22.6648    6.97104e3     3.38      3.67
  ! Ne     Be       2.5474    1.8309   1.9400  22.7750    1.06588e4     3.38      3.96
  ! Ar     Be       0.8082    3.2032   1.5058  37.1816    3.68450e4     3.38      5.63 
  ! Kr     Be       0.3844    5.3588   1.9600  61.452     1.67028e5     3.38      9.64
  ! Xe     Be       0.4779    8.1740   1.8350  86.942     4.23834e5     3.38     14.06
  !  Be    Be       2.0334    0.8241   1.3437  16.9689    2.20796e3     3.38      3.38

  data yield_parameters(:,:,1) /&
       0.8007,    0.0564,   1.5147,   14.340,    2.56510e2,     3.38,      9.32,  &
       1.7575,    0.1044,   1.9906,   9.5059,    2.82110e2,     3.38,      5.67,  &
       2.0794,    0.1379,   1.5660,   9.4345,    3.07966e2,     3.38,      4.49,  &
       0.7725,    0.3310,   1.6036,  12.8963,    6.65344e2,     3.38,      4.49,  &
       1.4745,    0.3193,   1.6989,  12.3288,    7.19545e2,     3.38,      3.97,  &
       5.2833,    0.9334,   2.5368,  16.5425,    5.46566e3,     3.38,      3.55,  &
       1.2209,    1.2024,   1.6881,  22.6648,    6.97104e3,     3.38,      3.67,  &
       2.5474,    1.8309,   1.9400,  22.7750,    1.06588e4,     3.38,      3.96,  &
       0.8082,    3.2032,   1.5058,  37.1816,    3.68450e4,     3.38,      5.63,  &
       0.3844,    5.3588,   1.9600,  61.452 ,    1.67028e5,     3.38,      9.64,  &
       0.4779,    8.1740,   1.8350,  86.942 ,    4.23834e5,     3.38,     14.06,  &
       2.0334,    0.8241,   1.3437,  16.9689,    2.20796e3,     3.38,      3.38   &
       /

  !
  ! Target Material 2: C Data
  !
  ! C data
  !
  ! Ion   Target    Lambda       q       mu    Eth (eV)    Epsilon     Esb (eV)  Esb/Gamma
  !
  !  H     C        1.3533    0.0241   1.4103   38.630    4.14659e2     7.41     25.89
  !  D     C        1.2848    0.0539   1.1977   27.770    4.46507e2     7.41     15.08 
  !  T     C        1.9050    0.0718   1.1512   23.617    4.78673e2     7.41     11.54
  ! 3He    C        0.7341    0.2058   1.1956   29.883    1.02061e3     7.41     11.54
  ! 4He    C        4.5910    0.1951   1.7852   19.124    1.08716e3     7.41      9.88
  !  N     C        5.4288    0.7481   1.7701  34.9372    7.37899e3     7.41      7.45
  !  O     C        9.6110    1.0171   2.0102  34.1293    9.29758e3     7.41      7.56
  ! Ne     C        2.5015    1.1912   1.6551  46.6904    1.39308e4     7.41      7.92 
  ! Ar     C        1.2622    2.4576   1.3952  68.8460    4.57989e4     7.41     10.42 
  ! Kr     C        1.3628    3.4372   2.2366  88.2918    1.99609e5     7.41     16.90
  ! Xe     C        0.4408    4.3004   1.7734  145.4236   4.98349e5     7.41     24.13
  !  C     C       13.9666    0.7015   2.0947  21.4457    5.68684e3     7.41      7.41
  !

  data yield_parameters(:,:,2) /&
       1.3533,    0.0241,   1.4103,   38.630 ,   4.14659e2,     7.41,     25.89,  &
       1.2848,    0.0539,   1.1977,   27.770 ,   4.46507e2,     7.41,     15.08,  & 
       1.9050,    0.0718,   1.1512,   23.617 ,   4.78673e2,     7.41,     11.54,  &
       0.7341,    0.2058,   1.1956,   29.883 ,   1.02061e3,     7.41,     11.54,  &
       4.5910,    0.1951,   1.7852,   19.124 ,   1.08716e3,     7.41,      9.88,  &
       5.4288,    0.7481,   1.7701,  34.9372 ,   7.37899e3,     7.41,      7.45,  &
       9.6110,    1.0171,   2.0102,  34.1293 ,   9.29758e3,     7.41,      7.56,  &
       2.5015,    1.1912,   1.6551,  46.6904 ,   1.39308e4,     7.41,      7.92,  & 
       1.2622,    2.4576,   1.3952,  68.8460 ,   4.57989e4,     7.41,     10.42,  & 
       1.3628,    3.4372,   2.2366,  88.2918 ,   1.99609e5,     7.41,     16.90,  &
       0.4408,    4.3004,   1.7734,  145.4236,   4.98349e5,     7.41,     24.13,  &
       13.9666,   0.7015,   2.0947,   21.4457,   5.68684e3,     7.41,      7.41   &
       /


  !
  ! Target Material 3: Mo Data
  !
  !
  ! Mo Data
  !
  ! Ion   Target    Lambda       q       mu    Eth (eV)    Epsilon     Esb (eV)  Esb/Gamma
  !           
  !  H     Mo       0.5124    0.0114   1,1469  201.4886   4.71832e3     6.83      165.63 
  !  D     Mo       0.3241    0.0326   1.5410   97.7738   4.76698e3     6.83       84.95
  !  T     Mo       0.5078    0.0661   1.5955   67.1475   4.81614e3     6.83       57.71
  ! 3He    Mo       0.3541    0.1373   0.9926   75.3995   9.84614e3     6.83       57.71
  ! 4He    Mo       0.1537    0.1563   0.9989   59.3088   9.94365e3     6.83       44.44
  !  N     Mo       0.1157    1.7900   1.8032   28.2561   4.10879e4     6.83       15.36 
  !  O     Mo       0.1762    2.1069   2.4821   23.9169   4.83320e4     6.83       13.94
  ! Ne     Mo       0.2205    2.8995   2.6514   23.6170   6.38956e4     6.83       11.89
  ! Ar     Mo       0.1339    6.3606   1.9562   28.2149   1.43274e5     6.83        8.23
  ! Kr     Mo       0.1412   11.9419   1.6911   38.6337   4.17411e5     6.83        6.86
  ! Xe     Mo       0.2401   32.5719   1.6694   47.4030   8.47848e5     6.83        7.00
  ! Mo     Mo       0.3580   12.2715   2.1844   31.4737   5.33049e5     6.83        6.83

  data yield_parameters(:,:,3) /&
       0.5124,    0.0114,   1.1469,  201.4886,   4.71832e3,     6.83,      165.63,  & 
       0.3241,    0.0326,   1.5410,   97.7738,   4.76698e3,     6.83,       84.95,  &
       0.5078,    0.0661,   1.5955,   67.1475,   4.81614e3,     6.83,       57.71,  &
       0.3541,    0.1373,   0.9926,   75.3995,   9.84614e3,     6.83,       57.71,  &
       0.1537,    0.1563,   0.9989,   59.3088,   9.94365e3,     6.83,       44.44,  &
       0.1157,    1.7900,   1.8032,   28.2561,   4.10879e4,     6.83,       15.36,  & 
       0.1762,    2.1069,   2.4821,   23.9169,   4.83320e4,     6.83,       13.94,  &
       0.2205,    2.8995,   2.6514,   23.6170,   6.38956e4,     6.83,       11.89,  &
       0.1339,    6.3606,   1.9562,   28.2149,   1.43274e5,     6.83,        8.23,  &
       0.1412,   11.9419,   1.6911,   38.6337,   4.17411e5,     6.83,        6.86,  &
       0.2401,   32.5719,   1.6694,   47.4030,   8.47848e5,     6.83,        7.00,  &
       0.3580,   12.2715,   2.1844,   31.4737,   5.33049e5,     6.83,        6.83   &
       /

  !
  ! Target Material 4: W Data
  !
  !
  ! W Data
  !
  ! Ion   Target    Lambda       q       mu    Eth (eV)    Epsilon     Esb (eV)  Esb/Gamma
  !           
  !  H     W        1.0087    0.0075   1.2046   457.42    9.86986e3     8.68      399.36  
  !  D     W        0.3583    0.0183   1.4410   228.84    9.92326e3     8.68      202.85 
  !  T     W        0.2870    0.0419   1.5802   153.8842  9.97718e3     8.68      136.48 
  ! 3He    W        0.2424    0.0884   1.2439   164.3474  2.02666e4     8.68      136.48
  ! 4He    W        0.1692    0.1151   1.7121   120.56    2.03728e4     8.68      104.13 
  !  N     W        0.0921    1.4389   2.0225   45.3362   7.90505e4     8.68       32.98
  !  O     W        0.0777    1.8824   1.7536   44.2135   9.19794e4     8.68       29.46
  ! Ne     W        0.0828    2.5520   1.9534   38.6389   1.19107e5     8.68       24.35
  ! Ar     W        0.2113    5.9479   2.3857   27.0503   2.46646e5     8.68       14.80
  ! Kr     W        0.1747   13.6917   2.5161   34.7592   6.36677e5     8.68       10.09
  ! Xe     W        0.1385   20.5321   2.0952   44.8701   1.18932e6     8.68        8.93
  !  W     W        2.2697   18.6006   3.1273   24.9885   1.99860e6     8.68        8.68

  data yield_parameters(:,:,4) /&
       1.0087,    0.0075,   1.2046,   457.42  ,  9.86986e3,     8.68,      399.36,  &  
       0.3583,    0.0183,   1.4410,   228.84  ,  9.92326e3,     8.68,      202.85,  & 
       0.2870,    0.0419,   1.5802,   153.8842,  9.97718e3,     8.68,      136.48,  & 
       0.2424,    0.0884,   1.2439,   164.3474,  2.02666e4,     8.68,      136.48,  &
       0.1692,    0.1151,   1.7121,   120.56  ,  2.03728e4,     8.68,      104.13,  & 
       0.0921,    1.4389,   2.0225,   45.3362 ,  7.90505e4,     8.68,       32.98,  &
       0.0777,    1.8824,   1.7536,   44.2135 ,  9.19794e4,     8.68,       29.46,  &
       0.0828,    2.5520,   1.9534,   38.6389 ,  1.19107e5,     8.68,       24.35,  &
       0.2113,    5.9479,   2.3857,   27.0503 ,  2.46646e5,     8.68,       14.80,  &
       0.1747,   13.6917,   2.5161,   34.7592 ,  6.36677e5,     8.68,       10.09,  &
       0.1385,   20.5321,   2.0952,   44.8701 ,  1.18932e6,     8.68,        8.93,  &
       2.2697,   18.6006,   3.1273,   24.9885 ,  1.99860e6,     8.68,        8.68   &
       /


  !
  ! The following arrays map the LIM/DIVIMP target and plasma materials to the indices used in this module - these lists
  ! must be updated by hand if the information changes in LIM/DIVIMP.
  !
  ! In addition, the setting of these is only performed on first execution 
  !
  !

  integer :: init_targ_mat, init_ion_mat,init_matt,init_matp
  logical :: eckstein2007_data_available


  integer,parameter :: targ_mats(19) = (/ &
       -1, &         !  MATT=1    Aluminum
       targ_be, &   !  MATT=2    Beryllium
       -1, &         !  MATT=3    Copper
       targ_c, &    !  MATT=4    Carbon/Graphite
       -1, &         !  MATT=5    Titanium
       -1, &         !  MATT=6    Iron
       -1, &         !  MATT=7    Nickel
       targ_mo, &   !  MATT=8    Molybdenum
       targ_w, &    !  MATT=9    Tungsten
       -1, &         !  MATT=10   Boron
       -1, &         !  MATT=11   Lithium
       -1, &         !  MATT=12   Chromium
       -1, &         !  MATT=13   "Deuterium"
       -1, &         !  MATT=14   "Helium"
       -1, &         !  MATT=15   "Neon"
       -1, &         !  MATT=16   "Argon"
       -1, &         !  MATT=17   "Oxygen"
       -1, &         !  MATT=18   "Chlorine"
       -1  &         !  MATT=19   "Nitrogen"
       /)

! slmod begin
  integer,parameter :: ion_mats(8) = (/ &
       ion_h, &    !  MATP=1    H
       ion_d, &    !  MATP=2    D
       ion_t, &    !  MATP=3    T
       ion_he4, &  !  MATP=4    HE4
       -1, &       !  MATP=5    C
       ion_self, & !  MATP=6    SELF
       ion_o,  &   !  MATP=7    O
       ion_ne &    !  MATP=8    Ne   ! Addition noted in SPUTTER.F as well. -SL, 10/05/2012 
       /) 
!
!  integer,parameter :: ion_mats(7) = (/ &
!       ion_h, &    !  MATP=1    H
!       ion_d, &    !  MATP=2    D
!       ion_t, &    !  MATP=3    T
!       ion_he4, &  !  MATP=4    HE4
!       -1, &       !  MATP=5    C
!       ion_self, & !  MATP=6    SELF
!       ion_o  &    !  MATP=7    O
!       /) 
! slmod end

  CHARACTER*18 :: TARMAT(19)= (/&
       ' ALUMINIUM       ',' BERYLLIUM       ',' COPPER          ',&
       ' GRAPHITE        ',' TITANIUM        ',' IRON            ',&
       ' NICKEL          ',' MOLYBDENUM      ',' TUNGSTEN        ',&
       ' BORON           ',' LITHIUM         ',' CHROMIUM        ',&
       ' "DEUTERIUM"     ',' "HELIUM"        ',' "NEON"          ',&
       ' "ARGON"         ',' "OXYGEN"        ',' "CHLORINE"      ',&
       ' "NITROGEN"      ' /)                                     

! slmod begin
  CHARACTER*6  :: PLAMAT(8) = (/ ' H    ',' D    ',' T    ',' HE4  ',' C    ',' SELF ',' O    ',' Ne   '  /)
!
!  CHARACTER*6  :: PLAMAT(7) = (/ ' H    ',' D    ',' T    ',' HE4  ',' C    ',' SELF ',' O    ',/)
! slmod end

! slmod begin
!  public :: yield_2007,init_eckstein_2007,print_eck2007_yields,eckstein2007_data_available,  &
!            get_target_index,get_plasma_index
!
  public :: yield_2007,init_eckstein_2007,print_eck2007_yields,eckstein2007_data_available
! slmod end

  save



contains
! slmod begin
!
!  integer function get_target_index(z)
!
!    integer, intent(in) :: z
!
!    get_target_index = -1
!
!    if (z.eq.4 ) get_target_index = 2  ! targ_be  MATT=2  Beryllium
!    if (z.eq.6 ) get_target_index = 4  ! targ_c   MATT=4  Carbon/Graphite
!    if (z.eq.42) get_target_index = 8  ! targ_mo  MATT=8  Molybdenum
!    if (z.eq.74) get_target_index = 9  ! targ_w   MATT=9  Tungsten
!
!    return 
!
!  end function get_target_index
!
!  integer function get_plasma_index(z,a)
!
!    integer, intent(in) :: z,a
!
!    get_plasma_index = -1
!
!    if (z.eq.1 .and.a.EQ.1) get_plasma_index =  1  !  ion_h     MATP=1    H
!    if (z.eq.1 .and.a.EQ.2) get_plasma_index =  2  !  ion_d     MATP=2    D
!    if (z.eq.1 .and.a.EQ.3) get_plasma_index =  3  !  ion_t     MATP=3    T
!    if (z.eq.2 .and.a.EQ.3) get_plasma_index = -1  !  ion_he3  
!    if (z.eq.2 .and.a.EQ.4) get_plasma_index =  4  !  ion_he4   MATP=4    HE4
!    if (z.eq.6            ) get_plasma_index =  5  !  ion_c     MATP=5
!    if (z.eq.7            ) get_plasma_index = -1  !  ion_n    
!    if (z.eq.8            ) get_plasma_index =  7  !  ion_o     MATP=7
!    if (z.eq.10           ) get_plasma_index = -1  !  ion_ne   
!    if (z.eq.18           ) get_plasma_index = -1  !  ion_ar   
!    if (z.eq.36           ) get_plasma_index = -1  !  ion_kr   
!    if (z.eq.54           ) get_plasma_index = -1  !  ion_xe   
!
!    if (z.EQ.-1           ) get_plasma_index =  6  !  ion_self  MATP=6    SELF
!
!    return 
!     
!  end function get_plasma_index
!
! slmod end

  subroutine init_eckstein_2007(matt,matp)

    integer :: matt,matp

    !
    ! This routine maps the LIM/DIVIMP indices to the internal routine indices
    !

    eckstein2007_data_available = .true.

    init_matt  = matt
    init_matp  = matp

    init_targ_mat = targ_mats(matt)
    init_ion_mat  = ion_mats(matp)

    if (init_targ_mat.eq.-1 .or.init_ion_mat.eq.-1) then 


       call errmsg('ECKSTEIN 2007 YIELD DATA','The combination of target ='//tarmat(matt)//&
            ' and bombarding ion = '//plamat(matp)//' is not supported in the 2007 data'//&
            ' currently available in the code. 1996 data will be used instead')

       eckstein2007_data_available = .false.

    endif

  end subroutine init_eckstein_2007





  real function yield_f(lambda,q,mu,eth,epsilon,e0)

    real,intent(in) :: lambda,q,mu,eth,epsilon,e0

    !
    ! Note: this function is small enough that it could be in-lined (in fact the compiler may do that
    !       anyway when optimization is on. 
    !

    !
    ! Calculate the yield using the Eckstein/Preuss formula and appropriate parameters
    !

    real :: el,w,s_n_krc

    !
    ! Check for impact energy below threshold
    !

    if (e0.le.eth) then 
       yield_f = 0.0
       return
    endif


    !
    ! Calculate the contributions to the equation
    !

    el = e0 / epsilon

    w = el + 0.1728 * sqrt(el) + 0.008 * el**0.1504

    s_n_krc = 0.5 * log(1.0 +1.2288 * el) / w


    !
    ! Calculate the yield
    !
    yield_f = q * s_n_krc * (e0/eth -1.0)**mu / (lambda/w + (e0/eth-1.0)**mu)



  end function yield_f


! slmod begin
  real function yield_2002_W_SS_60(e0) 
    real,intent(in) :: e0

    REAL :: result
    REAL :: yield(2,20) = (/  &
           15.00 ,  2.98E-6 ,  &
           20.00 ,  3.29E-5 ,  &
           25.00 ,  1.36E-4 ,  &
           30.00 ,  3.11E-4 ,  &
           40.00 ,  1.15E-3 ,  &
           50.00 ,  3.64E-3 ,  &
           60.00 ,  8.48E-3 ,  &
           70.00 ,  1.70E-2 ,  &
           80.00 ,  2.83E-2 ,  &
          100.00 ,  6.03E-2 ,  &
          120.00 ,  1.00E-1 ,  &
          140.00 ,  1.43E-1 ,  &
          200.00 ,  2.80E-1 ,  &
          300.00 ,  5.01E-1 ,  &
          350.00 ,  5.98E-1 ,  &
          500.00 ,  8.85E-1 ,  &
          800.00 ,  1.35E+0 ,  &
         1000.00 ,  1.62E+0 ,  &
         2000.00 ,  2.59E+0 ,  &
         2500.00 ,  2.98E+0  /)

      CALL Fitter(20,yield(1,1:20),yield(2,1:20),1,e0,result,'LINEAR')

      yield_2002_W_SS_60 = result

  end function yield_2002_W_SS_60
! slmod end

  real function yield_2007(matp,matt,e0) 
    integer, intent(in) :: matp,matt
    real,intent(in) :: e0

    !
    ! LIM and DIVIMP assign specific values to matt (target material) and matp (plasma or bombarding material) 
    ! These numbers have to be mapped to the indices used internally in this module. 
    !
    !

    real :: lambda,q,mu,eth,epsilon,yield2
    integer :: targ_mat, ion_mat

    targ_mat = targ_mats(matt)
    ion_mat  = ion_mats(matp)

    !

    lambda  = yield_parameters(lambda_val,ion_mat,targ_mat)
    q       = yield_parameters(q_val,ion_mat,targ_mat)
    mu      = yield_parameters(mu_val,ion_mat,targ_mat)
    eth     = yield_parameters(eth_val,ion_mat,targ_mat)
    epsilon = yield_parameters(epsilon_val,ion_mat,targ_mat)

    yield_2007 = yield_f(lambda,q,mu,eth,epsilon,e0)
! slmod begin

!    write(0,*) '  2007',matt,ion_mat,matp,targ_mat

    if (.false..and.ion_mat.eq.12.and.targ_mat.eq.4) then

      yield2 = yield_2002_W_SS_60(e0)
!      write(0,*) 'ss',e0,yield_2007,yield2
      yield_2007 = yield2
    endif

! slmod wend
  end function yield_2007


  subroutine print_eck2007_yields(ounit)
    implicit none
    !
    ! Prints the yield data currently in use to the specified unit
    !
    integer :: ounit
    integer :: in
    real :: e0,yieldp,yields

    if (.not.eckstein2007_data_available) return

    call prb
    write(ounit,'(a)')   '     2007 ECKSTEIN SPUTTERING DATA FIT PARAMETERS:'
    write(ounit,'(a,a)') '     TARGET MATERIAL IS:',tarmat(init_matt)
    write(ounit,'(a,a)') '     BOMBARDING ION IS :',plamat(init_matp)
    write(ounit,'(a)')   '     FIT PARAMETERS: '
    write(ounit,'(a)')   '            BOMBARDING ION      SELF-SPUTTERING '


    write(ounit,'(a,3x,2(g12.5,8x))') '  LAMBDA = ',yield_parameters(lambda_val,init_ion_mat,init_targ_mat) ,&
         yield_parameters(lambda_val,ion_self,init_targ_mat)    
    write(ounit,'(a,3x,2(g12.5,8x))') '  Q      = ',yield_parameters(q_val,init_ion_mat,init_targ_mat) ,&
         yield_parameters(q_val,ion_self,init_targ_mat)    
    write(ounit,'(a,3x,2(g12.5,8x))') '  MU     = ',yield_parameters(mu_val,init_ion_mat,init_targ_mat) ,&
         yield_parameters(mu_val,ion_self,init_targ_mat)    
    write(ounit,'(a,3x,2(g12.5,8x))') '  ETH    = ',yield_parameters(eth_val,init_ion_mat,init_targ_mat) ,&
         yield_parameters(eth_val,ion_self,init_targ_mat)    
    write(ounit,'(a,3x,2(g12.5,8x))') '  EPSILON= ',yield_parameters(epsilon_val,init_ion_mat,init_targ_mat) ,&
         yield_parameters(epsilon_val,ion_self,init_targ_mat)    

    call prb
    write(ounit,'(a)') '     YIELD FORMULA USED:'
    write(ounit,'(a)') '     Y(E) = q * S_n_KrC * [ E0/Eth -1 ]**mu / [  lambda/w + (E0/Eth -1)**mu] '
    write(ounit,'(a)') '     S_n_KrC =  0.5 ln (1+1.2288 EL) / w'
    write(ounit,'(a)') '     w = EL + 0.1728 * sqrt(EL) + 0.008 * EL**0.1504 '
    write(ounit,'(a)') '     EL = E0 (M2/(M1+M2)) * aL / (Z1*Z2*e**2) = E0 / Epsilon '
    write(ounit,'(a)') '     E0 is the impact energy of the bombarding species'
    call prb

    write(ounit,'(a)') '     SAMPLE YIELD DATA (NORMAL INCIDENCE): '
    write(ounit,'(a)') '     Energy (eV)       PRIMARY          SELF'

    do in = 0,14
       e0 = in * 10.0 + 10.0
       yieldp = yield_2007(init_matp,init_matt,e0)
       yields = yield_2007(6,init_matt,e0)
       write(ounit,'(5x,f10.2,2(5x,g12.5))') e0,yieldp,yields
    end do

    call prb

  end subroutine print_eck2007_yields




end module eckstein_2007_yield_data

