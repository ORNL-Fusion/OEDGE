! -*-Mode:f90-*-
! Storage_setup.f90
! Storage Setup File
!
! Adam McLean under Prof. Peter Stangeby
! Research Assistant: David Elder
! January, 2002
!
! This module sets up each of the common blocks
! and assigns data for cell-specific information.
!

Module HC_OUT_Storage_Setup

  Use ComHC ! Access to setup and input values.

  Implicit None

  ! Keep type declarations in memory.
  Save



  !Type HC_OUT_Data_Table_Type

  ! Contains all data read into OUT from the DIV binary output file.
  Character (Len=10), Dimension (Number_HC_Species) :: HC_State_List ! List of character names of each hydrocarbon species used in OUT processing.
  Integer, Dimension (Number_HC_Species+1) :: HC_Output_List ! Determines which hydrocarbons to include in stored binary file for OUT processing: 0=not included, 1=included.
  Real, Dimension (Max_Number_Walks,2) :: HC_Walks ! Follows position of particles as they travel for each timestep up to maxnws steps - typically 10,000.
  Real, Dimension (Number_HC_Species) :: HC_Factor_A ! Stores FACTA.
  Real, Dimension (Number_HC_Species) :: HC_Factor_B ! Stores FACTB.
  Real :: HC_FP ! Density normalization factor.
  Real :: HC_FT ! Density normalization factor.
  Real :: HC_ABSFAC ! Normalization for wall launches.
  !		Real :: FYTOT
  Integer :: HC_Particles_Launched ! Number of hydrocarbons launched.

  
  ! Allocatable arrays
  !Real, allocatable, Dimension (maxnks,maxnrs,Number_HC_Species) :: HC_Density ! Hydrocarbon density array. Records cell position and species of each fragment after each timestep.
  !Real, allocatable, Dimension (maxnks,maxnrs,Num_H_States) :: H_Density ! Hydrogen density array. Records cell position and species of each fragment at time of discard during breakup.
  !Real, allocatable, Dimension (maxnks,maxnrs,-1:0) :: HC_TIZS_CH ! CH ionization locations. (:,:,-1) for primary neutrals and (:,:,0) for those from primaries + reflections/self-sputters.
  Real, allocatable, Dimension (:,:,:) :: HC_Density ! Hydrocarbon density array. Records cell position and species of each fragment after each timestep.
  Real, allocatable, Dimension (:,:,:) :: H_Density ! Hydrogen density array. Records cell position and species of each fragment at time of discard during breakup.
  Real, allocatable, Dimension (:,:,:) :: HC_TIZS_CH ! CH ionization locations. (:,:,-1) for primary neutrals and (:,:,0) for those from primaries + reflections/self-sputters.
  !		Real, Dimension (maxnks,maxnrs,-1:0) :: HC_TIZS_C2 ! C2 ionization locations. (:,:,-1) for primary neutrals and (:,:,0) for those from primaries + reflections/self-sputters.

  ! Storage used solely in OUT.
  Real, allocatable, Dimension (:,:) :: HC_Trans_Prob ! Stores total transition probability to move out of one state given reactions possible and cell properties.
  Real, allocatable, Dimension (:,:) :: HC_KVALS
  Real, allocatable, Dimension (:,:) :: X_CH ! CH/CD excitation rate.
  Real, allocatable, Dimension (:,:) :: X_C2 ! C2 excitation rate.

  !End Type HC_Out_Data_Table_Type

  contains

    subroutine allocate_hc_out_storage
      use mod_params
      use allocate_arrays
      implicit none
      integer :: ierr
      
      call allocate_array(hc_density,maxnks,maxnrs,Number_HC_Species,'hc_density',ierr)
      call allocate_array(h_density,maxnks,maxnrs,Num_H_States,'h_density',ierr)
      call allocate_array(hc_tizs_ch,1,maxnks,1,maxnrs,-1,0,'hz_tizs_ch',ierr)

      ! Storage used solely in OUT.

      call allocate_array(hc_trans_prob,maxnks,maxnrs,'hc_trans_prob',ierr)
      call allocate_array(hc_kvals,maxnks,Number_HC_Species,'hc_kvals',ierr)
      call allocate_array(x_ch,maxnks,maxnrs,'x_ch',ierr)
      call allocate_array(x_c2,maxnks,maxnrs,'x_c2',ierr)

    end subroutine allocate_hc_out_storage

    subroutine deallocate_hc_out_storage
      implicit none

      if (allocated(hc_density)) deallocate(hc_density)
      if (allocated(h_density))  deallocate(h_density)
      if (allocated(hc_tizs_ch)) deallocate(hc_tizs_ch)
      if (allocated(hc_trans_prob)) deallocate(hc_trans_prob)
      if (allocated(hc_kvals)) deallocate(hc_kvals)
      if (allocated(x_ch)) deallocate(x_ch)
      if (allocated(x_c2)) deallocate(x_c2)
    
    end subroutine deallocate_hc_out_storage

  
End Module HC_OUT_Storage_Setup
