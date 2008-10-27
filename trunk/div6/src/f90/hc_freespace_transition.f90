! -*-Mode:f90-*-
! hc_freespace_transition.f90
! Hydrocarbon Reaction Program
!
! Adam McLean under Prof. Peter Stangeby
! Research Assistant: David Elder
! October, 2002
!
! The following subroutines control the physics of transitioning from one
! hydrocarbon species to another.  They control direction, velocity and
! energy of the newly-formed hydrocarbon fragment in free space.  They
! are split into four main groups:
!	1) Neutral -> Neutral Transitions
!	2) Neutral -> Ion Transitions
!	3) Ion -> Ion Transitions
!	4) Ion -> Neutral Transitions
! Note that transitions occuring due to interaction with the reactor
! vessel (wall and targets) is handled by the hc_vessel_interact module.
 
Module HC_FreeSpace_Transition
 
  ! Gain access to external data and subroutines.
  Use ComHC
  Use HC_Init_Lib_Data
  Use HC_Init_DIV_Data
  Use HC_Get
  Use HC_Utilities
 
  ! Every good Fortran program has...
  Implicit None
 
Contains
 
  Subroutine Neut_Neut_FreeSpace_Angle (HC_Reaction,Followed_State,Current_Angle,Azimuthal_Angle,Seed,NRand)
 
    ! Every good Fortran program has...		
    Implicit None
 
    Integer, Intent (In) :: HC_Reaction
    Integer, Intent (In) :: Followed_State
    Real, Intent (Out) :: Current_Angle
    Real, Intent (Out) :: Azimuthal_Angle
    Double Precision, Intent (In) :: Seed
    Integer, Intent (InOut) :: NRand
 
    ! Declare local variables.
    Real :: Random_Value_1
    Real :: Random_Value_2
    Real :: Random_Value_3
    Real :: Beta, Psi, HC_Angle
 
    NRand = NRand + 1
    Call Surand2 (Seed, 1, Random_Value_1)
    NRand  = NRand + 1
    Call Surand2 (Seed, 1, Random_Value_2)
    NRand = NRand + 1
    Call Surand2 (Seed, 1, Random_Value_3)
 
    ! Calculate polar angle theta.
    !Current_Angle = 1.0 *  Pi_Value * Random_Value_1
    Current_Angle = ACOS (1.0-2.0*Random_Value_1)
 
    ! Project velocity into 3D space, phi.
    Azimuthal_Angle = 2.0 *  Pi_Value * Random_Value_2 ! 0 to 2PI.
 
    Current_Angle = ATAN (TAN (Current_Angle) * COS (Azimuthal_Angle))
 
    !Beta = 0.5 *  Pi_Value * Random_Value_1 ! 0 to PI/2.
    !Psi = 2.0 *  Pi_Value * Random_Value_2 ! 0 to 2PI.
    !HC_Angle = ATAN (TAN(Beta) * COS(Psi))
 
    !If (Random_Value_3 .le. 0.5) Then
    !	Current_Angle = HC_Angle
    !Else
    !	If (HC_Angle .ge. 0.0) Then
    !		Current_Angle = HC_Angle -  Pi_Value
    !	Else
    !		Current_Angle = HC_Angle +  Pi_Value
    !	End If
    !End If
 
    ! Line up with +Z axis in DIVIMP.
    !Current_Angle = Current_Angle +  Pi_Value / 2.0
    !If (Current_Angle .gt.  Pi_Value) Then
    !	Current_Angle = Current_Angle - 2.0 *  Pi_Value
    !End If
 
    Current_Angle = Current_Angle -  Pi_Value / 2.0
    If (Current_Angle .lt. - Pi_Value) Then
       Current_Angle = Current_Angle + 2.0 *  Pi_Value
    End If
 
    ! Angle is now between -2PI and 0.
    Current_Angle = SIGN(Current_Angle,Random_Value_3-0.5)
 
  End Subroutine Neut_Neut_FreeSpace_Angle
 
  Subroutine Neut_Neut_FreeSpace_Energy (HC_Reaction,Followed_State,H_Isotope_Composition,Current_Cell,Current_Ring,&
                                        &Current_Velocity,Current_Energy,Kin_Energy_Added)
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables.
    Integer, Intent (In) :: HC_Reaction
    Integer, Intent (In) :: Followed_State
    Integer, Dimension (Number_H_Species), Intent (In) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (In) :: H_Isotope_Composition
    Integer, Intent (In) :: Current_Cell
    Integer, Intent (In) :: Current_Ring
    Real, Intent (In) :: Current_Velocity
    Real, Intent (InOut) :: Current_Energy
    Real, Intent (Out) :: Kin_Energy_Added
 
    ! Find energy from current velocity.
    If (hc_energy_calc .eq. 0) Then
       Current_Energy = Find_HC_Mass (HC_State_Transform_Table (HC_Reaction) % End_C_States (Followed_State),H_Isotope_Composition)&
                            & * (Current_Velocity / 1.38E4) ** 2
    ElseIf (hc_energy_calc .eq. 1) Then
       Current_Energy = Current_Energy
    End If
 
    If (HC_State_Transform_Table (HC_Reaction) % HC_E_Type .eq. 0 .or. HC_Reaction_Kinetics .eq. 0) Then
       ! Do nothing to energy as found above.
       Kin_Energy_Added = 0.0
    ElseIf (HC_State_Transform_Table (HC_Reaction) % HC_E_Type .eq. 1) Then
       ! Add fixed amount of energy to hydrocarbon product.
       Current_Energy = Current_Energy + HC_State_Transform_Table (HC_Reaction) % HC_E
       Kin_Energy_Added = HC_State_Transform_Table (HC_Reaction) % HC_E
    ElseIf (HC_State_Transform_Table (HC_Reaction) % HC_E_Type .eq. 2) Then
       ! Add variable energy to hydrocarbon product corresponding to some fraction of the background energy.
       ! Note: Num_Electron_Reactions is found in HC_Init_Lib_Data.
       ! jdemod - removing reliance on some fixed number of electron reactions - code would never have worked for higher hydrocarbons
       if (hc_state_transform_table(hc_reaction)%reaction_type.eq.'p') then
       !If (HC_Reaction .gt. Num_Electron_Reactions) Then
          ! Use local electron plasma temperature.
          Current_Energy = Current_Energy + HC_State_Transform_Table (HC_Reaction) % HC_E * gktebs (Current_Cell,Current_Ring)
          Kin_Energy_Added = HC_State_Transform_Table (HC_Reaction) % HC_E * gktebs (Current_Cell,Current_Ring)
       Else
          ! Use local proton plasma temperature.
          Current_Energy = Current_Energy + HC_State_Transform_Table (HC_Reaction) % HC_E * gktibs (Current_Cell,Current_Ring)
          Kin_Energy_Added = HC_State_Transform_Table (HC_Reaction) % HC_E * gktibs (Current_Cell,Current_Ring)
       End If
    Else
       Write (Output_Unit_HC_Alert,*) "Error in Neut_Neut_Freespace_Energy: Unsupported energetic: ",HC_State_Transform_Table (&
                                          &HC_Reaction) % HC_E_Type
       Write (Output_Unit_HC_Alert,*) "Stopping program."
       Stop
    End If
 
  End Subroutine Neut_Neut_FreeSpace_Energy
 
  Subroutine Neut_Neut_FreeSpace_Velocity (HC_Reaction,Followed_State,H_Isotope_Composition,Current_Angle,Azimuthal_Angle,&
                                          &Current_Velocity,Current_Velocity_In_R,Current_Velocity_In_Z,Seed,NRand)
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables.
    Integer, Intent (In) :: HC_Reaction
    Integer, Intent (In) :: Followed_State
    Integer, Dimension (Number_H_Species), Intent (In) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (In) :: H_Isotope_Composition
    Real, Intent (In) :: Current_Angle
    Real, Intent (In) :: Azimuthal_Angle
    Real, Intent (InOut) :: Current_Velocity
    Real, Intent (Out) :: Current_Velocity_In_R
    Real, Intent (Out) :: Current_Velocity_In_Z
    Double Precision, Intent (In) :: Seed
    Integer, Intent (InOut) :: NRand
 
    ! Declare local variables.
    Real :: Random_Temp_1
    Real :: Energetic_Vel ! Velocity associated with hydrocarbon energetic.
 
    ! Find velocity associated with hydrocarbon energetic.
    Energetic_Vel = 1.38E4 * SQRT (HC_State_Transform_Table (HC_Reaction) % HC_E / Find_HC_Mass (HC_State_Transform_Table (&
                            &HC_Reaction) % End_C_States (Followed_State),H_Isotope_Composition)) ! Mass in amu, T in eV.
    Energetic_Vel = 0.0
    ! Find random values to determine which direction extra velocity acts.
    NRand  = NRand + 1
    Call Surand2 (Seed, 1, Random_Temp_1)
 
    ! Add extra velocity component.
    ! Note: SIGN (A,B) = Abs (A) * Sign (B)
    Current_Velocity = Current_Velocity + SIGN (0.5 * Energetic_Vel, Random_Temp_1 - 0.5)
 
    !write (0,*) "NEUTNEUTVEL",Current_Velocity,Energetic_Vel,Random_Temp_1,Azimuthal_Angle,Current_Velocity*COS(Azimuthal_Angle)
    ! Project velocity into 3D space.
    Current_Velocity = Current_Velocity * COS(Azimuthal_Angle)
    !Current_Velocity = 0.0 ! REMOVE AMMOD.
    ! Calculate R,Z components of velocity, probability of ionization, etc.
    Current_Velocity_In_R = Current_Velocity * COS (Current_Angle) *  Neutral_Time_Step
    Current_Velocity_In_Z = Current_Velocity * SIN (Current_Angle) *  Neutral_Time_Step
 
  End Subroutine Neut_Neut_FreeSpace_Velocity
 
  Subroutine Neut_Neut_FreeSpace_Location (HC_Reaction,Followed_State,Current_R,Current_Z,Current_Cell,Current_Ring)
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables.
    Integer, Intent (In) :: HC_Reaction
    Integer, Intent (In) :: Followed_State
    Real, Intent (InOut) :: Current_R
    Real, Intent (InOut) :: Current_Z
    Integer, Intent (InOut) :: Current_Cell
    Integer, Intent (InOut) :: Current_Ring
 
  End Subroutine Neut_Neut_FreeSpace_Location
 
  Subroutine Neut_Ion_FreeSpace_Angle (HC_Reaction,Followed_State,Current_Angle)
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables.
    Integer, Intent (In) :: HC_Reaction
    Integer, Intent (In) :: Followed_State
    Real, Intent (Out) :: Current_Angle
 
    ! Calculate angle from cell geometry.
    Current_Angle = 0.0				
 
  End Subroutine Neut_Ion_FreeSpace_Angle
 
  Subroutine Neut_Ion_FreeSpace_Energy (HC_Reaction,Followed_State,H_Isotope_Composition,Current_Cell,Current_Ring,&
                                       &Current_Velocity,Current_Energy,Kin_Energy_Added)
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables.
    Integer, Intent (In) :: HC_Reaction
    Integer, Intent (In) :: Followed_State
    Integer, Dimension (Number_H_Species), Intent (In) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (In) :: H_Isotope_Composition
    Integer, Intent (In) :: Current_Cell
    Integer, Intent (In) :: Current_Ring
    Real, Intent (In) :: Current_Velocity
    Real, Intent (InOut) :: Current_Energy
    Real, Intent (Out) :: Kin_Energy_Added
 
    If (hc_energy_calc .eq. 0) Then
       Current_Energy = Find_HC_Mass (HC_State_Transform_Table (HC_Reaction) % End_C_States (Followed_State),H_Isotope_Composition)&
                     & * (Current_Velocity / 1.38E4) ** 2
    ElseIf (hc_energy_calc .eq. 1) Then
       Current_Energy = Current_Energy
    End If
 
    If (HC_State_Transform_Table (HC_Reaction) % HC_E_Type .eq. 0 .or. HC_Reaction_Kinetics .eq. 0) Then
       ! Do nothing to energy as found above.
       Kin_Energy_Added = 0.0
    ElseIf (HC_State_Transform_Table (HC_Reaction) % HC_E_Type .eq. 1) Then
       ! Add fixed amount of energy to hydrocarbon product.
       Current_Energy = Current_Energy + HC_State_Transform_Table (HC_Reaction) % HC_E
       Kin_Energy_Added = HC_State_Transform_Table (HC_Reaction) % HC_E
    ElseIf (HC_State_Transform_Table (HC_Reaction) % HC_E_Type .eq. 2) Then
       ! Add variable energy to hydrocarbon product corresponding to some fraction of the background energy.
       ! Note: Num_Electron_Reactions is found in HC_Init_Lib_Data.
       ! jdemod - removing reliance on some fixed number of electron reactions - code would never have worked for higher hydrocarbons
       if (hc_state_transform_table(hc_reaction)%reaction_type.eq.'p') then
       ! If (HC_Reaction .gt. Num_Electron_Reactions) Then
          ! Use local electron plasma temperature.
          Current_Energy = Current_Energy + HC_State_Transform_Table (HC_Reaction) % HC_E * gktebs (Current_Cell,Current_Ring)
          Kin_Energy_Added = HC_State_Transform_Table (HC_Reaction) % HC_E * gktebs (Current_Cell,Current_Ring)
       Else
          Current_Energy = Current_Energy + HC_State_Transform_Table (HC_Reaction) % HC_E * gktibs (Current_Cell,Current_Ring)
          Kin_Energy_Added = HC_State_Transform_Table (HC_Reaction) % HC_E * gktibs (Current_Cell,Current_Ring)
       End If
    Else
       Write (Output_Unit_HC_Alert,*) "Error in Neut_Neut_Freespace_Energy: Unsupported energetic: ",HC_State_Transform_Table (&
                                           &HC_Reaction) % HC_E_Type
       Write (Output_Unit_HC_Alert,*) "Stopping program."
       Stop
    End If
 
    If (.not.  Grid_Error) Then
       ! IF SET TI=TB FOR STATE IZ APPLIES, BETTER DO IT
       If (Get_HC_Charge (HC_State_Transform_Table (HC_Reaction) % End_C_States (Followed_State)) .eq.  Equate_Ion_Temp_Charge) Then
          Current_Energy = MAX (Current_Energy, gktibs (Current_Cell,Current_Ring))
       End If
    Else
       Write (Output_Unit_HC_Alert,*) "Warning: Ion created off the grid:",HC_State_Transform_Table (HC_Reaction) % End_C_States (&
                                   &Followed_State),Current_Cell,Current_Ring
    End If
 
  End Subroutine Neut_Ion_FreeSpace_Energy
 
  Subroutine Neut_Ion_FreeSpace_Velocity (HC_Reaction,Followed_State,H_Isotope_Composition,Current_Ring,Current_S,Current_Velocity,&
                                         &Current_Velocity_In_S,Seed,NRand)
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables.
    Integer, Intent (In) :: HC_Reaction
    Integer, Intent (In) :: Followed_State
    Integer, Dimension (Number_H_Species), Intent (In) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (In) :: H_Isotope_Composition
    Integer, Intent (In) :: Current_Ring
    Real, Intent (In) :: Current_S
    Real, Intent (InOut) :: Current_Velocity ! Note that neutral Current_Velocity is set to 0.0.
    Real, Intent (Out) :: Current_Velocity_In_S
    Double Precision, Intent (In) :: Seed
    Integer, Intent (InOut) :: NRand
 
    ! Declare local variables.
    Real :: Random_Temp
    Real :: Energetic_Vel ! Velocity associated with hydrocarbon energetic.
 
    ! Add extra velocity associated with hydrocarbon energetic.
    !Energetic_Vel = 9.79E3 * SQRT (HC_State_Transform_Table (HC_Reaction) % HC_E / Find_HC_Mass 
    ! (HC_State_Transform_Table (HC_Reaction) % End_C_States (Followed_State),H_Isotope_Composition)) ! Mass in amu, T in eV.

    Energetic_Vel = 1.38E4 * SQRT (HC_State_Transform_Table (HC_Reaction) % HC_E / Find_HC_Mass (HC_State_Transform_Table (&
                   &HC_Reaction) % End_C_States (Followed_State),H_Isotope_Composition)) ! Mass in amu, T in eV.
 
    ! Initial ion velocity (H23).  Typically the same as CNEUTG DIVIMP input option.
    If (hc_neut_ion_velocity .eq. 0) Then
       ! Dead stop along S.
       Current_Velocity_In_S = 0.0
    ElseIf (hc_neut_ion_velocity .eq. 1 .or. hc_neut_ion_velocity .eq. 3) Then
       ! +/- 0.5 Vprevious along S.
       Random_Temp = geran (Seed)
       NRand  = NRand + 1
       Current_Velocity_In_S = SIGN (0.5 * Current_Velocity, Random_Temp - 0.5)
       If (hc_neut_ion_velocity .eq. 3) Then
          Random_Temp = geran (Seed)
          NRand = NRand + 1
          Current_Velocity_In_S = Current_Velocity_In_S * 2.0 * SQRT (Random_Temp)
       End If
 
    ElseIf (hc_neut_ion_velocity .eq. 2) Then
       Current_Velocity_In_S = SIGN (Current_Velocity, 0.5 * gksmaxs (Current_Ring) - Current_S)
    End If
 
    ! Find random values to determine which direction extra velocity acts.
    Random_Temp = geran (Seed)
    NRand  = NRand + 1
    ! Add extra velocity component.
    ! Note: SIGN (A,B) = Abs (A) * Sign (B)
    Current_Velocity_In_S = Current_Velocity_In_S + SIGN (0.5 * Energetic_Vel, Random_Temp - 0.5)
    !Current_Velocity_In_S = 0.0 ! REMOVE AMMOD, but not below.
    ! Reset neutral particle velocity.
    Current_Velocity = 0.0
 
  End Subroutine Neut_Ion_FreeSpace_Velocity
 
  Subroutine Neut_Ion_FreeSpace_Location (HC_Reaction,Followed_State,Current_R,Current_Z,Current_Cell,Current_Ring,STmp,Current_S,&
                                         &Current_Cross,Current_Theta,S_Start)
 
    use error_handling
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables.
    Integer, Intent (In) :: HC_Reaction
    Integer, Intent (In) :: Followed_State
    Real, Intent (In) :: Current_R
    Real, Intent (In) :: Current_Z
    Integer, Intent (In) :: Current_Cell
    Integer, Intent (In) :: Current_Ring
    Real, Intent (In) :: STmp
    Real, Intent (Out) :: Current_S
    Real, Intent (Out) :: Current_Cross
    Real, Intent (Out) :: Current_Theta
    Real, Intent (Out) :: S_Start
 
    ! SET Initial S and CROSS postion for particles.
    ! Note:  Check for control option is removed given that any HC can
    ! be launched with control option 0 so that should always be used.
    If ( Neutral_Init_Pos_Opt .eq. 0) Then
       Current_Cross = 0.0
       If (STmp .ne. 0.0 .and. &
            & ( Injection_Opt .eq.2 .or.  Injection_Opt .eq.3 .or.  Injection_Opt .eq. 5 .or.  Injection_Opt .eq. 6)) Then
          Current_S  = STmp
       Else
          Current_S  = gkss (Current_Cell,Current_Ring)
       End If
 
    ElseIf ( Neutral_Init_Pos_Opt .eq. 1) Then
 
       Current_Cross = 0.0
       If (STmp .ne. 0.0 .and. &
            & ( Injection_Opt .eq. 2 .or.  Injection_Opt .eq. 3 .or.  Injection_Opt .eq. 5 .or.  Injection_Opt .eq. 6)) Then
          Current_S  = STmp
       Else
          Call getscross_approx (Current_R, Current_Z, Current_S, Current_Cross, Current_Cell, Current_Ring)
       End If
    End If
 
    ! Record starting S-distance from nearest target.
    S_Start = MIN (Current_S, gksmaxs (Current_Ring) - Current_S)
 
    If (.not.  Grid_Error) Then
       ! Set Theta value for non-orthogonal transport
       If ( Non_Orthogonal_Grid_Opt .eq. 1 .or.  Non_Orthogonal_Grid_Opt .eq. 3) Then
          ! Non-orthogonal.
          If (Current_S .gt. gkss (Current_Cell, Current_Ring)) Then
             If (Current_Cell .lt. gnks (Current_Ring)) Then
                Current_Theta = gthetag (Current_Cell, Current_Ring) + (Current_S - gkss (Current_Cell, Current_Ring)) / gkfords (&
                     &Current_Cell, Current_Ring) * (gthetag (Current_Cell+1, Current_Ring) - gthetag (Current_Cell, Current_Ring))
             Else
                Current_Theta = gthetag (Current_Cell, Current_Ring) + (Current_S - gkss (Current_Cell, Current_Ring)) / gkfords (&
                          &Current_Cell, Current_Ring) * (gthetat (gidds (Current_Ring,1))- gthetag (Current_Cell, Current_Ring))
             End If
          Else
             If (Current_Cell .gt. 1) Then
                Current_Theta = gthetag (Current_Cell, Current_Ring) + (gkss (Current_Cell, Current_Ring) - Current_S) / gkbacds (&
                        &Current_Cell, Current_Ring)*(gthetag (Current_Cell, Current_Ring) - gthetag (Current_Cell-1,Current_Ring))
             Else
                Current_Theta = gthetag (Current_Cell, Current_Ring) + (gkss (Current_Cell, Current_Ring) - Current_S) / gkbacds (&
                         &Current_Cell, Current_Ring) * (gthetag (Current_Cell, Current_Ring) - gthetat (gidds (Current_Ring,2)))
             End If
          End If
       End If
    Else
       Write (Output_Unit_HC_Alert,*) "Warning: Ion created off the grid:",Current_Cell,Current_Ring
    End If
 
  End Subroutine Neut_Ion_FreeSpace_Location
 
  Subroutine Ion_Ion_FreeSpace_Angle (HC_Reaction,Followed_State,Current_Angle)
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables.
    Integer, Intent (In) :: HC_Reaction
    Integer, Intent (In) :: Followed_State
    Real, Intent (Out) :: Current_Angle
 
    ! Calculate angle from cell geometry.
    Current_Angle = 0.0				
 
  End Subroutine Ion_Ion_FreeSpace_Angle
 
  Subroutine Ion_Ion_FreeSpace_Energy (HC_Reaction,Followed_State,H_Isotope_Composition,Current_Cell,Current_Ring,&
                                     &Current_Velocity_In_S,Current_Energy,Kin_Energy_Added)
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables.
    Integer, Intent (In) :: HC_Reaction
    Integer, Intent (In) :: Followed_State
    Integer, Dimension (Number_H_Species), Intent (In) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (In) :: H_Isotope_Composition
    Integer, Intent (In) :: Current_Cell
    Integer, Intent (In) :: Current_Ring
    Real, Intent (In) :: Current_Velocity_In_S
    Real, Intent (InOut) :: Current_Energy
    Real, Intent (Out) :: Kin_Energy_Added
 
    ! Find temperature.
    If (hc_energy_calc .eq. 0) Then
       !Current_Energy = Find_HC_Mass (HC_State_Transform_Table (HC_Reaction) % End_C_States 
       ! (Followed_State),H_Isotope_Composition) * (Current_Velocity_In_S / 9.79E3) ** 2
       Current_Energy = Find_HC_Mass (HC_State_Transform_Table (HC_Reaction) % End_C_States (Followed_State),H_Isotope_Composition)&
                                      & * (Current_Velocity_In_S / 1.38E4) ** 2
    ElseIf (hc_energy_calc .eq. 1) Then
       Current_Energy = Current_Energy
    End If
 
    If (HC_State_Transform_Table (HC_Reaction) % HC_E_Type .eq. 0 .or. HC_Reaction_Kinetics .eq. 0) Then
       ! Do nothing to energy as found above.
       Kin_Energy_Added = 0.0
    ElseIf (HC_State_Transform_Table (HC_Reaction) % HC_E_Type .eq. 1) Then
       ! Add fixed amount of energy to hydrocarbon product.
       Current_Energy = Current_Energy + HC_State_Transform_Table (HC_Reaction) % HC_E
       Kin_Energy_Added = HC_State_Transform_Table (HC_Reaction) % HC_E
    ElseIf (HC_State_Transform_Table (HC_Reaction) % HC_E_Type .eq. 2) Then
       ! Add variable energy to hydrocarbon product corresponding to some fraction of the background energy.
       ! Note: Num_Electron_Reactions is found in HC_Init_Lib_Data.
       ! jdemod - removing reliance on some fixed number of electron reactions - code would never have worked for higher hydrocarbons
       if (hc_state_transform_table(hc_reaction)%reaction_type.eq.'p') then
       !If (HC_Reaction .gt. Num_Electron_Reactions) Then
          ! Use local electron plasma temperature.
          Current_Energy = Current_Energy + HC_State_Transform_Table (HC_Reaction) % HC_E * gktebs (Current_Cell,Current_Ring)
          Kin_Energy_Added = HC_State_Transform_Table (HC_Reaction) % HC_E * gktebs (Current_Cell,Current_Ring)
       Else
          Current_Energy = Current_Energy + HC_State_Transform_Table (HC_Reaction) % HC_E * gktibs (Current_Cell,Current_Ring)
          Kin_Energy_Added = HC_State_Transform_Table (HC_Reaction) % HC_E * gktibs (Current_Cell,Current_Ring)
       End If
    Else
       Write (Output_Unit_HC_Alert,*) "Error in Neut_Neut_Freespace_Energy: Unsupported energetic: ",HC_State_Transform_Table (&
                                    &HC_Reaction) % HC_E_Type
       Write (Output_Unit_HC_Alert,*) "Stopping program."
       Stop
    End If
 
    If (.not.  Grid_Error) Then
       ! IF SET TI=TB FOR STATE IZ APPLIES, BETTER DO IT
       If (Get_HC_Charge (HC_State_Transform_Table (HC_Reaction) % End_C_States (Followed_State)) .eq.  Equate_Ion_Temp_Charge) Then
          Current_Energy = MAX (Current_Energy, gktibs (Current_Cell,Current_Ring))
       End If
    Else
       Write (Output_Unit_HC_Alert,*) "Warning: Ion created off the grid:",HC_State_Transform_Table (HC_Reaction) % End_C_States (&
                              &Followed_State),Current_Cell,Current_Ring
    End If
 
  End Subroutine Ion_Ion_FreeSpace_Energy
 
  Subroutine Ion_Ion_FreeSpace_Velocity (HC_Reaction,Followed_State,H_Isotope_Composition,Current_Ring,Current_S,&
                                       &Current_Velocity_In_S,Seed,NRand)
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables.
    Integer, Intent (In) :: HC_Reaction
    Integer, Intent (In) :: Followed_State
    Integer, Dimension (Number_H_Species), Intent (In) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (In) :: H_Isotope_Composition
    Integer, Intent (In) :: Current_Ring
    Real, Intent (In) :: Current_S
    Real, Intent (InOut) :: Current_Velocity_In_S
    Double Precision, Intent (In) :: Seed
    Integer, Intent (InOut) :: NRand
 
    ! Declare local variables.
    Real :: Random_Temp
    Real :: Energetic_Vel ! Velocity associated with hydrocarbon energetic.
 
    ! Find velocity associated with hydrocarbon energetic.
    !Energetic_Vel = 9.79E3 * SQRT (HC_State_Transform_Table (HC_Reaction) % HC_E / Find_HC_Mass 
    ! (HC_State_Transform_Table (HC_Reaction) % End_C_States (Followed_State),H_Isotope_Composition)) ! Mass in amu, T in eV.
    Energetic_Vel = 1.38E4 * SQRT (HC_State_Transform_Table (HC_Reaction) % HC_E / Find_HC_Mass (HC_State_Transform_Table (&
                      &HC_Reaction) % End_C_States (Followed_State),H_Isotope_Composition)) ! Mass in amu, T in eV.
 
    ! Find random values to determine which direction extra velocity acts.
    Random_Temp = geran (Seed)
    NRand  = NRand + 1
 
    ! Add extra velocity component.
    ! Note: SIGN (A,B) = Abs (A) * Sign (B)
    Current_Velocity_In_S = Current_Velocity_In_S + SIGN (0.5 * Energetic_Vel, Random_Temp - 0.5)
    !Current_Velocity_In_S = 0.0 ! REMOVE AMMOD
 
  End Subroutine Ion_Ion_FreeSpace_Velocity
 
  Subroutine Ion_Ion_FreeSpace_Location (HC_Reaction,Followed_State,Current_Cell,Current_Ring,Current_S,Current_Cross,&
                                        &Current_Theta,S_Start)
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables.
    Integer, Intent (In) :: HC_Reaction
    Integer, Intent (In) :: Followed_State
    Integer, Intent (In) :: Current_Cell
    Integer, Intent (In) :: Current_Ring
    Real, Intent (InOut) :: Current_S
    Real, Intent (InOut) :: Current_Cross
    Real, Intent (Out) :: Current_Theta
    Real, Intent (Out) :: S_Start
 
    ! S, Cross, Cell, Ring are the same.
 
    ! Set Theta value for non-orthogonal transport
    If ( Non_Orthogonal_Grid_Opt .eq. 1 .or.  Non_Orthogonal_Grid_Opt .eq. 3) Then
       ! Non-orthogonal.
 
       If (Current_S .gt. gkss (Current_Cell, Current_Ring)) Then
          If (Current_Cell .lt. gnks (Current_Ring)) Then
             Current_Theta = gthetag (Current_Cell, Current_Ring) + (Current_S - gkss (Current_Cell, Current_Ring)) / gkfords (&
                  &Current_Cell, Current_Ring) * (gthetag (Current_Cell+1, Current_Ring) - gthetag (Current_Cell, Current_Ring))
          Else
             Current_Theta = gthetag (Current_Cell, Current_Ring) + (Current_S - gkss (Current_Cell, Current_Ring)) / gkfords (&
                         &Current_Cell, Current_Ring) * (gthetat (gidds (Current_Ring,1))- gthetag (Current_Cell, Current_Ring))
          End If
       Else
          If (Current_Cell .gt. 1) Then
             Current_Theta = gthetag (Current_Cell, Current_Ring) + (gkss (Current_Cell, Current_Ring)-Current_S) / gkbacds (&
                    &Current_Cell, Current_Ring)*(gthetag (Current_Cell, Current_Ring) - gthetag (Current_Cell-1,Current_Ring))
          Else
             Current_Theta = gthetag (Current_Cell, Current_Ring) + (gkss (Current_Cell, Current_Ring) - Current_S) / gkbacds (&
                        &Current_Cell, Current_Ring) * (gthetag (Current_Cell, Current_Ring) - gthetat (gidds (Current_Ring,2)))
          End If
       End If
    End If
 
    ! Record starting S-distance from nearest target.
    S_Start = MIN (Current_S, gksmaxs (Current_Ring) - Current_S)
 
  End Subroutine Ion_Ion_FreeSpace_Location
 
  Subroutine Ion_Neut_FreeSpace_Angle(HC_Reaction,Followed_State,Current_Cell,Current_Ring,Current_Angle,Azimuthal_Angle,Seed,NRand)
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables.
    Integer, Intent (In) :: HC_Reaction
    Integer, Intent (In) :: Followed_State
    Integer, Intent (In) :: Current_Cell
    Integer, Intent (In) :: Current_Ring
    Real, Intent (Out) :: Current_Angle
    Real, Intent (Out) :: Azimuthal_Angle
    Double Precision, Intent (In) :: Seed
    Integer, Intent (InOut) :: NRand
 
    ! Declare local variables.
    Real :: Random_Value_1
    Real :: Random_Value_2
    Real :: Beta, Psi, HC_Angle
 
    ! Find new angle for neutral.
    ! hc_ion_neut_angle ! ion->neutral transition angle distribution option: 0=isotropic, 1=sine biased forward, 2=S dir (H24).
 
    If (hc_ion_neut_angle .eq. 0) Then
 
       ! Isotropic distribution angle (DIVIMP default).
       ! Note: SIGN (A,B) returns the Abs(A) * sign(B)
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value_1)
       Current_Angle = 2.0 *  Pi_Value * Random_Value_1 -  Pi_Value ! -PI to PI.
 
       ! Project velocity into 3D space, phi.
       NRand  = NRand + 1
       Call Surand2 (Seed, 1, Random_Value_2)
       Azimuthal_Angle = 2.0 *  Pi_Value * Random_Value_2 -  Pi_Value ! -PI to PI.
 
       !NRand = NRand + 1
       !Call Surand2 (Seed, 1, Random_Value_2)
       !NRand = NRand + 1
       !Call Surand2 (Seed, 1, Random_Value_3)
 
       !Beta = 0.5 *  Pi_Value * Random_Value_1 ! 0 to PI/2.
       !Psi = 2.0 *  Pi_Value * Random_Value_2 ! 0 to 2PI.
       !HC_Angle = ATAN (TAN(Beta) * COS(Psi))
 
       !If (Random_Value_3 .le. 0.5) Then
       !	Current_Angle = HC_Angle
       !Else
       !	If (HC_Angle .ge. 0.0) Then
       !		Current_Angle = HC_Angle -  Pi_Value
       !	Else
       !		Current_Angle = HC_Angle +  Pi_Value
       !	End If
       !End If
 
       ! Line up with +Z axis in DIVIMP.
       !Current_Angle = Current_Angle +  Pi_Value / 2.0
       !If (Current_Angle .gt.  Pi_Value) Then
       !	Current_Angle = Current_Angle - 2.0 *  Pi_Value
       !End If
 
    ElseIf (hc_ion_neut_angle .eq. 1) Then
       ! Sine biased forward angle.
       ! Note: SIGN (A,B) returns the Abs(A) * sign(B)
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value_1)
       Current_Angle = SIGN (ACOS (Random_Value_1), Random_Value_2 - 0.5) + Find_Ion_Angle (Current_Cell, Current_Ring)
 
    ElseIf (hc_ion_neut_angle .eq. 2) Then
       ! Neutral gets last ion's magnetic field direction angle.
       Current_Angle = Find_Ion_Angle (Current_Cell, Current_Ring)
 
    End If
 
  End Subroutine Ion_Neut_FreeSpace_Angle
 
  Subroutine Ion_Neut_FreeSpace_Energy (HC_Reaction,Followed_State,H_Isotope_Composition,Current_Cell,Current_Ring,&
&Current_Velocity_In_S,Current_Energy,Kin_Energy_Added)
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables.
    Integer, Intent (In) :: HC_Reaction
    Integer, Intent (In) :: Followed_State
    Integer, Dimension (Number_H_Species), Intent (In) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (In) :: H_Isotope_Composition
    Integer, Intent (In) :: Current_Cell
    Integer, Intent (In) :: Current_Ring
    Real, Intent (In) :: Current_Velocity_In_S
    Real, Intent (InOut) :: Current_Energy
    Real, Intent (Out) :: Kin_Energy_Added
 
    ! Find temperature.
    If (hc_energy_calc .eq. 0) Then
       Current_Energy = Find_HC_Mass (HC_State_Transform_Table (HC_Reaction) % End_C_States (Followed_State),H_Isotope_Composition)&
                        & * (Current_Velocity_In_S / 1.38E4) ** 2
    ElseIf (hc_energy_calc .eq. 1) Then
       Current_Energy = Current_Energy
    End If
 
    If (HC_State_Transform_Table (HC_Reaction) % HC_E_Type .eq. 0 .or. HC_Reaction_Kinetics .eq. 0) Then
       ! Do nothing to energy as found above.
       Kin_Energy_Added = 0.0
    ElseIf (HC_State_Transform_Table (HC_Reaction) % HC_E_Type .eq. 1) Then
       ! Add fixed amount of energy to hydrocarbon product.
       Current_Energy = Current_Energy + HC_State_Transform_Table (HC_Reaction) % HC_E
       Kin_Energy_Added = HC_State_Transform_Table (HC_Reaction) % HC_E
    ElseIf (HC_State_Transform_Table (HC_Reaction) % HC_E_Type .eq. 2) Then
       ! Add variable energy to hydrocarbon product corresponding to some fraction of the background energy.
       ! Note: Num_Electron_Reactions is found in HC_Init_Lib_Data.
       ! jdemod - removing reliance on some fixed number of electron reactions - code would never have worked for higher hydrocarbons
       if (hc_state_transform_table(hc_reaction)%reaction_type.eq.'p') then
       !If (HC_Reaction .gt. Num_Electron_Reactions) Then
          ! Use local electron plasma temperature.
          Current_Energy = Current_Energy + HC_State_Transform_Table (HC_Reaction) % HC_E * gktebs (Current_Cell,Current_Ring)
          Kin_Energy_Added = HC_State_Transform_Table (HC_Reaction) % HC_E * gktebs (Current_Cell,Current_Ring)
       Else
          Current_Energy = Current_Energy + HC_State_Transform_Table (HC_Reaction) % HC_E * gktibs (Current_Cell,Current_Ring)
          Kin_Energy_Added = HC_State_Transform_Table (HC_Reaction) % HC_E * gktibs (Current_Cell,Current_Ring)
       End If
    Else
       Write (Output_Unit_HC_Alert,*) "Error in Neut_Neut_Freespace_Energy: Unsupported energetic: ",HC_State_Transform_Table (&
                  &HC_Reaction) % HC_E_Type
       Write (Output_Unit_HC_Alert,*) "Stopping program."
       Stop
    End If
 
  End Subroutine Ion_Neut_FreeSpace_Energy
 
  Subroutine Ion_Neut_FreeSpace_Velocity (HC_Reaction,Followed_State,H_Isotope_Composition,Current_Angle,Azimuthal_Angle,&
&Current_Cell,Current_Ring,HC_Temperature,Current_Velocity_In_S,Current_Velocity,Current_Velocity_In_R,Current_Velocity_In_Z,Seed,&
&NRand)
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables.
    Integer, Intent (In) :: HC_Reaction
    Integer, Intent (In) :: Followed_State
    Integer, Dimension (Number_H_Species), Intent (In) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (In) :: H_Isotope_Composition
    Real, Intent (In) :: Current_Angle
    Real, Intent (In) :: Azimuthal_Angle
    Integer, Intent (In) :: Current_Cell
    Integer, Intent (In) :: Current_Ring
    Real, Intent (In) :: HC_Temperature
    Real, Intent (InOut) :: Current_Velocity_In_S ! Note that neutral Current_Velocity_In_S is set to 0.0.
    Real, Intent (Out) :: Current_Velocity
    Real, Intent (Out) :: Current_Velocity_In_R
    Real, Intent (Out) :: Current_Velocity_In_Z
    Double Precision, Intent (In) :: Seed
    Integer, Intent (InOut) :: NRand
 
    ! Declare local variables.
    Real :: Random_Temp
    Real :: Energetic_Vel ! Velocity associated with hydrocarbon energetic.
 
    ! Initial neutral velocity.  Check for neutral heating option.
    If ( Impurity_Neutral_Vel_Opt .eq. 0) Then
       Current_Velocity = Current_Velocity_In_S *  Vel_Mult_Recomb_Neut
       !Current_Velocity = 1.38E4 * SQRT (HC_Temperature / Find_HC_Mass (HC_State_Transform_Table (HC_Reaction) % 
       ! End_C_States (Followed_State),H_Isotope_Composition)) *  Vel_Mult_Recomb_Neut
    Else If ( Impurity_Neutral_Vel_Opt .eq. 1 .or.  Impurity_Neutral_Vel_Opt .eq. 2) Then
       Current_Velocity = 1.38E4 * SQRT (gktibs (Current_Cell,Current_Ring) / Find_HC_Mass (HC_State_Transform_Table (HC_Reaction) &
                      &% End_C_States (Followed_State),H_Isotope_Composition)) *  Vel_Mult_Recomb_Neut
    End If
 
    ! Find velocity associated with hydrocarbon energetic.
    Energetic_Vel = 1.38E4 * SQRT (HC_State_Transform_Table (HC_Reaction) % HC_E / Find_HC_Mass (HC_State_Transform_Table (&
            &HC_Reaction) % End_C_States (Followed_State),H_Isotope_Composition)) ! Mass in amu, T in eV.
 
    ! Find random values to determine which direction extra velocity acts.
    NRand  = NRand + 1
    Call Surand2 (Seed, 1, Random_Temp)
 
    ! Add extra velocity component.
    ! Note: SIGN (A,B) = Abs (A) * Sign (B)
    Current_Velocity = Current_Velocity + SIGN (0.5 * Energetic_Vel, Random_Temp - 0.5)
 
    ! Project velocity into 3D space.
    Current_Velocity = Current_Velocity * COS(Azimuthal_Angle)
    !Current_Velocity = 0.0 ! REMOVE AMMOD
    ! Calculate R,Z components of velocity, probability of ionization, etc.
    Current_Velocity_In_R = Current_Velocity * COS (Current_Angle) *  Neutral_Time_Step
    Current_Velocity_In_Z = Current_Velocity * SIN (Current_Angle) *  Neutral_Time_Step
 
    ! Reset neutral particle velocity.
    Current_Velocity_In_S = 0.0
 
  End Subroutine Ion_Neut_FreeSpace_Velocity
 
  Subroutine Ion_Neut_FreeSpace_Location (HC_Reaction,Followed_State,Current_Cell,Current_Ring,Current_S,Current_Cross,Current_R,&
                                         &Current_Z,Current_Theta,S_Start)
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables
    Integer, Intent (In) :: HC_Reaction
    Integer, Intent (In) :: Followed_State
    Integer, Intent (In) :: Current_Cell
    Integer, Intent (In) :: Current_Ring
    Real, Intent (InOut) :: Current_S
    Real, Intent (InOut) :: Current_Cross
    Real, Intent (Out) :: Current_R
    Real, Intent (Out) :: Current_Z
    Real, Intent (Out) :: Current_Theta
    Real, Intent (Out) :: S_Start
 
    ! Update position of particle.  Note, all that is needed is R and Z.
    Call getrz (Current_Cell, Current_Ring, Current_S, Current_Cross, Current_R, Current_Z, RZ_Opt)
    ! Returns R, Z.
 
    ! Reset non-applicable variables.
    Current_S = 0.0
    Current_Cross = 0.0
    Current_Theta = 0.0
    S_Start = 0.0
 
  End Subroutine Ion_Neut_FreeSpace_Location
 
End Module HC_Freespace_Transition
