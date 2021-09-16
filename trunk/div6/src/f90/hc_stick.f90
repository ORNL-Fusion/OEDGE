! -*-Mode:f90-*-
! HC_Stick.f90
! Find the reflection probability (1 - sticking probability)
! for a particular HC impinging on the wall/target at a
! specific energy, angle to the surface normal, and
! surface temperature.
! Based on sticking coefficients calculated by:
! -1) Abritrary default preset (typically 0.5).
! -2) Alman & Ruzic, 2002.
!
! Adam McLean under Prof. Peter Stangeby
! Research Assistant: David Elder
! July, 2002
!
! This routine is called upon determination of reflection
! for impinged HC.  Called in the HC_Wrapper program.
 
Module HC_Stick
 
  ! Required modules.
  Use HC_Init_Lib_Data ! Gain access to HC species designations.	
 
  ! Every good Fortran program has...
  Implicit None	
 
Contains
 
  ! Routine to take care of HC sticking coefficients.
  Subroutine HC_Sputtering_Coefficient (Vessel_Segment,HC_Species,HC_Energy,HC_Angle_To_Normal,Sputtering_Probability)
 
    ! Required modules.
    Use ComHC ! Access to global variables related to the hydrocarbon module.
    Use HC_Utilities ! Access hydrocarbon utilities.
    Use HC_Init_Lib_Data ! Get HC charge.
    Use HC_Get ! Wallpts.
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables
    Integer, Intent (In) :: Vessel_Segment
    Integer, Intent (In) :: HC_Species
    Real, Intent (In) :: HC_Energy
    Real, Intent (In) :: HC_Angle_To_Normal
    Real, Intent (Out) :: Sputtering_Probability
 
    ! Define local variables that will be assigned values.
    Integer Local_HC_Sticking_Coef_Model ! Pre-set sputtering probability, Alman and Ruzic, or Janev models.
    Real :: Surface_Temperature
 
    ! Set local variables from common blocks.
    Local_HC_Sticking_Coef_Model = hc_sticking_coef_model
    Surface_Temperature = Find_Wall_Temperature (Vessel_Segment)
 
    ! Find new hydrocarbon sticking probability.
    If (Local_HC_Sticking_Coef_Model .eq. 0) Then
       ! Use preset sputtering table to decide which single hydrocarbon species is emitted.
       If (Get_HC_Charge (HC_Species) .eq. 0) Then
          ! Assign neutral sputtering coefficient.
          Sputtering_Probability = 1.0 - hc_sticking_coef_preset
       ElseIf (Get_HC_Charge (HC_Species) .ge. 1) Then
          ! Assign ion sputtering coefficient.
          Sputtering_Probability = 1.0 - hc_sticking_coef_preset
       Else
          ! Error - give stern warning.
          Write (Output_Unit_HC_Alert,*) "Error: Hydocarbon charge cannot be less than zero."
          Write (Output_Unit_HC_Alert,*) "Program stopping."
          Stop
       End If
    ElseIf (Local_HC_Sticking_Coef_Model .eq. 1) Then
       ! Use Alman and Ruzic sticking data.
       Call Stick_Alman_Ruzic (Vessel_Segment,HC_Species,HC_Energy,HC_Angle_To_Normal,Surface_Temperature, Sputtering_Probability)
 
    ElseIf (Local_HC_Sticking_Coef_Model .eq. 2) Then
       ! Use Janev reflection data.
       Call Stick_Janev (Vessel_Segment,HC_Species,HC_Energy,HC_Angle_To_Normal,Surface_Temperature,Sputtering_Probability)
    Else
       ! Unsupported sputtering model selected.
       Write (Output_Unit_HC_Alert,*) "Unsupported sputtering model selected."
       Write (Output_Unit_HC_Alert,*) "Program Stopping."
       Stop
    End If
 
  End Subroutine HC_Sputtering_Coefficient
 
  Subroutine HC_Reflection_Coefficient (Vessel_Segment,HC_Species,HC_Energy,HC_Angle_To_Normal,Reflection_Probability)
 
    ! Required modules.
    Use ComHC ! Access to global variables related to the hydrocarbon module.
    Use HC_Utilities ! Access hydrocarbon utilities.
    Use HC_Init_Lib_Data ! Get HC charge.
    Use HC_Get ! Wallpts.
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables
    Integer, Intent (In) :: Vessel_Segment
    Integer, Intent (In) :: HC_Species
    Real, Intent (In) :: HC_Energy
    Real, Intent (In) :: HC_Angle_To_Normal
    Real, Intent (Out) :: Reflection_Probability
 
    ! Define local variables that will be assigned values.
    Integer Local_HC_Reflection_Coef_Model ! Pre-set reflection probability, Alman and Ruzic, or Janev models.
    Real :: Surface_Temperature
 
    ! Set local variables from common blocks.
    Local_HC_Reflection_Coef_Model = hc_reflection_coef_model
    Surface_Temperature = Find_Wall_Temperature (Vessel_Segment)
 
    ! Find new hydrocarbon reflection probability.
    If (Local_HC_Reflection_Coef_Model .eq. 0) Then
       ! Use preset reflection table to decide which single hydrocarbon species is emitted.
       If (Get_HC_Charge (HC_Species) .eq. 0) Then
          ! Assign neutral reflection coefficient.
          Reflection_Probability = hc_reflection_coef_preset
 
          !write (0,*) "Reflection coefficient:",Reflection_Probability
 
       ElseIf (Get_HC_Charge (HC_Species) .ge. 1) Then
          ! Assign ion reflection coefficient.
          Reflection_Probability = hc_reflection_coef_preset
       Else
          ! Error - give stern warning.
          Write (Output_Unit_HC_Alert,*) "Error: Hydocarbon charge cannot be less than zero."
          Write (Output_Unit_HC_Alert,*) "Program stopping."
          Stop
       End If
    ElseIf (Local_HC_Reflection_Coef_Model .eq. 1) Then
       ! Use Janev reflection data.
       !Call Reflect_Janev (Vessel_Segment, HC_Species, HC_Energy, HC_Angle_To_Normal, Surface_Temperature, Reflection_Probability)
 
    ElseIf (Local_HC_Reflection_Coef_Model .eq. 2) Then
       ! Use Alman and Ruzic reflection data.
       !Call Reflect_Alman_Ruzic (Vessel_Segment, HC_Species, HC_Energy, HC_Angle_To_Normal, Surface_Temperature, Reflection_Probability)
 
       ! Apply reflection coefficients from Alman and Ruzic, JNM 313-316, 2003
       ! for average energies of wall impact.
       If (HC_Species .eq. 10 .or. HC_Species .eq. 9) Then
          ! CH4 or CH4+, thermal launch assumed.
          Reflection_Probability = 1.0
       ElseIf (HC_Species .eq. 8 .or. HC_Species .eq. 7) Then
          ! CH3 or CH3+, average energy of ~0.5 eV.
          Reflection_Probability = 0.92
       ElseIf (HC_Species .eq. 6 .or. HC_Species .eq. 5) Then
          ! CH3 or CH3+, average energy of ~2.0 eV.
          Reflection_Probability = 0.50
       ElseIf (HC_Species .eq. 4 .or. HC_Species .eq. 3) Then
          ! CH3 or CH3+, average energy of ~3.5 eV.
          Reflection_Probability = 0.20
       ElseIf (HC_Species .eq. 2 .or. HC_Species .eq. 1) Then
          ! CH3 or CH3+, average energy of ~5.5 eV.
          Reflection_Probability = 0.25
       Else
          ! unsupported species.
          Write (0,*) "Unsupported species at reflection:",HC_Species
          Write (0,*) "Using reflection coefficient of 1.0."
          Reflection_Probability = 1.0
       End If
 
    ElseIf (Local_HC_Reflection_Coef_Model .eq. 3) Then
 
       ! Apply constant specified reflection coefficients
       ! for all HC fragments except CH4.
       If (HC_Species .eq. 10) Then
          ! Assign perfect reflection coefficient for neutral CH4.
          Reflection_Probability = 1.0
       ElseIf (HC_Species .eq. 2) Then
          ! Assign perfect absorption coefficient for neutral C.
          Reflection_Probability = 0.0
       Else
          ! Assign specified reflection coefficient.
          Reflection_Probability = hc_reflection_coef_preset
       End If
 
    ElseIf (Local_HC_Reflection_Coef_Model .eq. 4) Then
       ! Use data from Jacob - J.Nuc.Mat. 337-339 (2005) 839-846 : Table 1 pg 844
       ! The table gives sticking coefficients - converted here to reflection coefficients
       If (HC_Species .eq. 10 .or. HC_Species .eq. 9) Then
          ! CH4 or CH4+
          Reflection_Probability = 1.0
       ElseIf (HC_Species .eq. 8 .or. HC_Species .eq. 7) Then
          ! CH3 or CH3+
          Reflection_Probability = 1.0-1.0e-4
       ElseIf (HC_Species .eq. 6 .or. HC_Species .eq. 5) Then
          ! CH2 or CH2+
          Reflection_Probability = 1.0-0.025
       ElseIf (HC_Species .eq. 4 .or. HC_Species .eq. 3) Then
          ! CH or CH+
          Reflection_Probability = 0.0
       ElseIf (HC_Species .eq. 2 .or. HC_Species .eq. 1) Then
          ! C or C+
          Reflection_Probability = 0.0
       Else
          ! unsupported species.
          Write (0,*) "Unsupported species at reflection:",HC_Species
          Write (0,*) "Using reflection coefficient of 1.0."
          Reflection_Probability = 1.0
       End IF

    ElseIf (Local_HC_Reflection_Coef_Model .eq. 5) Then
       ! Reflection for all species except C and C+
       If (HC_Species .eq. 2 .or. HC_Species .eq. 1) Then
          ! C or C+
          Reflection_Probability = 0.0
       elseif (hc_species.ge.3.and.hc_species.le.10) then 
          ! All other HC species up to CH4
          Reflection_Probability = hc_reflection_coef_preset
       Else
          ! unsupported species.
          Write (0,*) "Unsupported species at reflection:",HC_Species
          Write (0,*) "Using reflection coefficient of 1.0."
          Reflection_Probability = 1.0
       End If

    Else
       ! Unsupported reflection model selected.
       Write (Output_Unit_HC_Alert,*) "Unsupported reflection model selected."
       Write (Output_Unit_HC_Alert,*) "Program Stopping."
       Write (0,*) "Unsupported reflection model selected."
       Write (0,*) "Program Stopping."
       Stop
    End If
 
  End Subroutine HC_Reflection_Coefficient
 
  Subroutine Stick_Alman_Ruzic(Vessel_Segment,HC_Species,HC_Energy,HC_Angle_To_Normal,Surface_Temperature,Reflection_Probability)
    ! Use Alman and Ruzic PSI2002 paper data for reflection species.
 
    Use HC_Utilities ! Access hydrocarbon utilities.
 
    ! Every good Fortran 90 program has...		
    Implicit None
 
    ! Define input/output variables
    Integer, Intent (In) :: Vessel_Segment
    Integer, Intent (In) :: HC_Species
    Real, Intent (In) :: HC_Energy
    Real, Intent (In) :: HC_Angle_To_Normal
    Real, Intent (In) :: Surface_Temperature
    Real, Intent (Out) :: Reflection_Probability
 
    ! Define local variables that will be assigned values.
 
 
  End Subroutine Stick_Alman_Ruzic
 
  Subroutine Stick_Janev (Vessel_Segment, HC_Species, HC_Energy, HC_Angle_To_Normal, Surface_Temperature, Reflection_Probability)
    ! Use Alman and Ruzic PSI2002 paper data for reflection species.
 
    Use HC_Utilities ! Access hydrocarbon utilities.
 
    ! Every good Fortran 90 program has...		
    Implicit None
 
    ! Define input/output variables
    Integer, Intent (In) :: Vessel_Segment
    Integer, Intent (In) :: HC_Species
    Real, Intent (In) :: HC_Energy
    Real, Intent (In) :: HC_Angle_To_Normal
    Real, Intent (In) :: Surface_Temperature
    Real, Intent (Out) :: Reflection_Probability
 
    ! Define local variables that will be assigned values.
 
 
  End Subroutine Stick_Janev
 
 
End Module HC_Stick
