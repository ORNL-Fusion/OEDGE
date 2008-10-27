! -*-Mode:f90-*-
! HC_Release.f90
! Upon impact of a neutral with a first wall surface, this routine
! decides what hydrocarbon is released into the vessel by chemical
! sputtering.
!
! Adam McLean under Prof. Peter Stangeby
! Research Assistant: David Elder
! July, 2002
!
! This routine is called upon determination of species of HC
! initially launched.  Called in the HC_Wrapper program.	//depot/div6/src/f90/hc_inside_ion.d6b	# edit
 
 
Module HC_Release
 
  ! Required modules.
  Use HC_Init_Lib_Data ! Gain access to HC species designations.	
 
  ! Every good Fortran 90 program has...
  Implicit None	
 
Contains
 
  ! Routine to take care of HC released.
  Subroutine Release_HC (R_Location,Z_Location,Vessel_Segment,HC_Energy_In,HC_Angle_To_Normal_In,Sputtered_HC_Species)
 
    ! Required modules.
    Use ComHC ! Variables set by DIVIMP initialization.
    Use HC_Stack ! Gain access to carbon stack storage and stack operation subroutines.
    Use HC_Init_DIV_Data ! Global and cell data from DIVIMP common blocks.
    Use HC_Utilities ! Surface temp.
 
    ! Every good Fortran 90 program has...		
    Implicit None
 
    ! Define input/output variables
    Real, Intent (In) :: R_Location
    Real, Intent (In) :: Z_Location
    Integer, Intent (In) :: Vessel_Segment
    Real, Intent (In) :: HC_Energy_In
    Real, Intent (In) :: HC_Angle_To_Normal_In
    Integer, Intent (Out) :: Sputtered_HC_Species
 
    ! Define local variables that will be assigned values.
    Integer :: Local_HC_Sputtering_Model ! Sputtered release of preset HC or using Mech&Davis&Haasz.
    Integer :: Local_HC_Sputtered_HC_Species ! Preset HC to be released if HC_Sputter_Type set to default.
    Real :: Local_Back_Plasma_Mass ! Used in Mech/Haasz/Davis data.
    Real :: Surface_Temperature ! Also used in Mech/Haasz/Davis data.
 
    ! Set local variables from common blocks.
    Local_HC_Sputtering_Model = hc_sputtering_model ! Assign from ComHC.
    Local_HC_Sputtered_HC_Species = hc_sputtered_hc_species ! Assign from ComHC.
    Local_Back_Plasma_Mass =  Back_Plasma_Ion_Mass ! 2.0 for deuterium, 2.5 for 50/50 DT.
 
    ! Determine sputtered hydrocarbon molecule.
    If (Local_HC_Sputtering_Model .eq. 0) Then
       ! Preset HC used.
       Sputtered_HC_Species = Local_HC_Sputtered_HC_Species ! Typically Methane, CH4, species 10.
    ElseIf (Local_HC_Sputtering_Model .eq. 1) Then
       ! Use Mech, Davis, Haasz correlations for chemical sputtering release distribution.
       ! Find surface temperature at impinged wall/target segment.
       Surface_Temperature = Find_Wall_Temperature (Vessel_Segment)
 
       ! Return sputtered hydrocarbon species
       Call HC_Sputter_Mech_Haasz_Davis (HC_Energy_In, HC_Angle_To_Normal_In, Local_Back_Plasma_Mass, Surface_Temperature, &
&Sputtered_HC_Species)
    Else
       ! Other options not currently supported.
       Write(Output_Unit_HC_Alert,*)"Option",Local_HC_Sputtering_Model,"for Local_HC_Sputtering_Model not supported at this time."
       Write (Output_Unit_HC_Alert,*) "Program stop."
       Stop
    End If
 
  End Subroutine Release_HC
 
  Subroutine HC_Sputter_Mech_Haasz_Davis (HC_Energy_In, HC_Angle_To_Normal_In, Background_Plasma_Mass, Surface_Temperature, &
&Sputtered_HC_Species)
 
    ! Implements release model derived from data from Mech, Haasz, and Davis, 1998.
 
    ! Every good Fortran 90 program has...		
    Implicit None
 
    ! Define input/output variables
    Real, Intent (In) :: HC_Energy_In
    Real, Intent (In) :: HC_Angle_To_Normal_In
    Real, Intent (In) :: Background_Plasma_Mass
    Real, Intent (In) :: Surface_Temperature
    Integer, Intent (Out) :: Sputtered_HC_Species
 
    ! Declare local variables.
 
    ! Note, yhaasz97 info in neut.d6a is needed, along with parameterized hc species data.
 
  End Subroutine HC_Sputter_Mech_Haasz_Davis
 
End Module HC_Release
