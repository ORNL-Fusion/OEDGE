module mod_soledge_input


  implicit none

  ! The purpose of this module is to split the input variables to the soledge
  ! module from the plasma calculation routines.
  ! This allows this module to be used in OUT without all the plasma code
  ! that is not required.

  ! jdemod - option to turn on/off use of SOL 12,13 etc
  integer,public :: soledge_opt
  
  integer,public :: cioptf_soledge


  integer,public :: csopt != 0
  integer,public :: cpopt != 0

  real,public :: csolls != 0.5
  real,public :: csollt != 0.5

  real,public :: csollr != 0.5  ! length of radiation source
  
  real,public :: csolfr != 0.0  ! strength of radiation as fraction of target power  (Popt 2 and 3)
  real,public :: csolpr != 0.0  ! strength of radiation source in absolute terms (popt 0 and 1) 
  
  real,public :: cfiz !=  0.0   ! fraction split between two ionization sources in sopt 4 and 5
  real,public :: sol13_pdist != 0.0
  real,public :: sol13_padd  != 0.0
 



  
contains

  



end module mod_soledge_input
