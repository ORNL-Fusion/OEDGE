! -*-Mode:f90-*-
! HC_Vessel_Interact.f90
! Upon HC impinging on the vessel wall and reflecting or sputtering,
! decides which HC species the reflecting HC becomes.
! Based on sticking and reflection coefficients calculated by:
! -1) Abritrary default preset (typically 0.5).
! -2) Alman & Ruzic, 2002.
!
! Adam McLean under Prof. Peter Stangeby
! Research Assistant: David Elder
! July, 2002
!
! This routine is called upon determination of reflection or sputtering
! for impinged HC.  Called in the HC_Wrapper program.
 
Module HC_Vessel_Interact
 
  ! Every good Fortran program has...
  Implicit None	
 
Contains
 
  Subroutine HC_Ion_Impact_ReLaunch (Type,Cur_HC_Spec,Last_HC_Species,H_Isotope_Composition,Last_H_Isotope_Composition,Launch_Reg, &
       & Sput_Weight,Neutral_Time_Step_Count,Ion_Time_Step_Count,Eq_Total_Ion_Time_Steps,Reflection_Counter,NeutType,Target_Index,&
       & Current_Theta, &
       & Current_R,Current_Z,HC_Temperature,Current_Velocity,Current_Velocity_In_S,Current_Angle,Seed,NRand,Current_Cell,&
       & Current_Ring, &
       & Current_Velocity_In_R,Current_Velocity_In_Z,Current_S,Current_Cross,hc_v)
    ! Particle is sputtered back into plasma.  Find new species and update all statistics.
 
    ! Required modules.
    Use ComHC ! Preset variables.
    Use HC_Init_DIV_Data ! Gain access to cell/global/geometry/current data structures.
    Use HC_Init_DIV_Diag ! Data diagnostics.
    Use HC_Init_Lib_Data ! Hydrocarbon data and utilities.
    Use HC_Utilities ! General purpose utilities.
    Use HC_Get ! gwallindex.
    use hc_velocity_type ! hc_v type
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables
    Character (len=*), Intent (In) :: Type
    Integer, Intent (InOut) :: Cur_HC_Spec
    Integer, Intent (Out) :: Last_HC_Species
    Integer, Dimension (Number_H_Species), Intent (InOut) :: H_Isotope_Composition
    Integer, Dimension (Number_H_Species), Intent (Out) :: Last_H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (InOut) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (Out) :: Last_H_Isotope_Composition
    Integer, Intent (InOut) :: Launch_Reg
    Real, Intent (InOut) :: Sput_Weight ! SPUTY
    Double Precision, Intent (InOut) :: Neutral_Time_Step_Count
    Double Precision, Intent (InOut) :: Ion_Time_Step_Count
    Double Precision, Intent (InOut) :: Eq_Total_Ion_Time_Steps
    Integer, Intent (InOut) :: Reflection_Counter
    Integer, Intent (In) :: NeutType
    ! jdemod - it appears that this should be intent(out)
    !Integer, Intent (In) :: Target_Index
    Integer, Intent (Out) :: Target_Index
    Real, Intent (InOut) :: Current_Theta
    Real, Intent (InOut) :: Current_R
    Real, Intent (InOut) :: Current_Z
    Real, Intent (InOut) :: HC_Temperature
    Real, Intent (InOut) :: Current_Velocity
    ! jdemod - it appears that current_velocity_is_s should be inout or out
    !Real, Intent (In) :: Current_Velocity_In_S
    Real, Intent (InOut) :: Current_Velocity_In_S
    Real, Intent (InOut) :: Current_Angle
    Double Precision, Intent (In) :: Seed
    Integer, Intent (InOut) :: NRand
    Integer, Intent (InOut) :: Current_Cell
    Integer, Intent (InOut) :: Current_Ring
    Real, Intent (InOut) :: Current_Velocity_In_R
    Real, Intent (InOut) :: Current_Velocity_In_Z
    Real, Intent (InOut) :: Current_S
    Real, Intent (InOut) :: Current_Cross
    
    type(hc_velocity_type1) :: hc_v
 
    ! Define local variables.
    Integer :: Orig_Launch_Reg
    !real :: temp_velocity_multiplier
 
    ! Record previous HC data.
    Last_HC_Species = Cur_HC_Spec
    Last_H_Isotope_Composition = H_Isotope_Composition
    Orig_Launch_Reg = Launch_Reg
 
    ! Record collision with target.
    HC_Num_Striking_Target (Cur_HC_Spec,Launch_Reg) =  HC_Num_Striking_Target (Cur_HC_Spec,Launch_Reg) + Sput_Weight
    Reflection_Counter = Reflection_Counter + 1
    Total_HC_Reflections (Cur_HC_Spec,Launch_Reg) =  Total_HC_Reflections (Cur_HC_Spec,Launch_Reg) + 1
    Max_HC_Reflections_Found ( Launch_HC_Species,Launch_Reg) = MAX0 (Reflection_Counter, Max_HC_Reflections_Found ( &
                                                               &Launch_HC_Species,Launch_Reg))
 
    ! Reflection-has-occurred indicator reset to false for each particle launched.
    HC_Hit_Vessel = .True.
 
    ! Reset time for new particle launch.
    Neutral_Time_Step_Count = 0.0
    Ion_Time_Step_Count = 0.0
    Eq_Total_Ion_Time_Steps = 0.0
 
    ! Newly launched particle properties.
    If (Type .eq. "Sputter") Then
       Call HC_Ion_Target_Location (Current_S,Current_Cross,Current_Cell,Current_Ring,Current_R,Current_Z,Target_Index)
       Call HC_Target_Sputter_Species (Cur_HC_Spec,Target_Index,Current_Theta,Current_R,Current_Z,HC_Temperature,Current_Velocity,&
                                      &Current_Velocity_In_S,Current_Angle)
       Call Get_HC_Isotope_Composition (Cur_HC_Spec,H_Isotope_Composition,Seed,NRand) ! Update new H content from target.
       Sput_Weight = Calc_Sputy (Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition))
       Call HC_Ion_Target_Sputter_Angle (Cur_HC_Spec,Current_Angle,HC_Temperature,Target_Index,Current_Angle,Seed,NRand)
       Call HC_Ion_Target_Sputter_Energy (Cur_HC_Spec,H_Isotope_Composition,Current_Cell,Current_Ring,Current_Angle,Target_Index,&
                                          &HC_Temperature)
       Call HC_Ion_Target_Sputter_Velocity (Cur_HC_Spec,H_Isotope_Composition,Current_Cell,Current_Ring,Current_Angle,Target_Index,&
                                   &HC_Temperature,Current_Velocity,Current_Velocity_In_R,Current_Velocity_In_Z,Seed,NRand,hc_v)

       ! Reset launch type for sputtering event.
       If (NeutType .eq. 2) Then
          ! Originally a target launch.
          If (Launch_Reg .eq. HC_Launch_Reg_Target_1_Dist .or. Launch_Reg .eq. HC_Launch_Reg_Target_1_Pin) Then
             Launch_Reg = HC_Launch_Reg_Sput_Target_1
          ElseIf (Launch_Reg .eq. HC_Launch_Reg_Target_2_Dist .or. Launch_Reg .eq. HC_Launch_Reg_Target_2_Pin) Then
             Launch_Reg = HC_Launch_Reg_Sput_Target_2
          End If
       ElseIf (NeutType .eq. 3) Then
          ! Originally a target self-sputter launch.
       ElseIf (NeutType .eq. 5) Then
          ! Originally a wall launch.
          Launch_Reg = HC_Launch_Reg_Sput_Wall
       ElseIf (NeutType .eq. 6) Then
          ! Originally a 2D neutral launch.
          ! No change.
       ElseIf (NeutType .eq. 7) Then
          ! Originally a reflected ion.
          ! No change.
       elseif (neuttype.eq.0) then 
          ! jdemod
          ! Freespace launch - for puffing ... - essentially equivalent to the 2D neutral launch except at one point
          ! No change - but need to check if neuttype is used in any arrays as an index
       Else
          ! Should not happen for HC launch.
          Write (Output_Unit_HC_Alert,*) "Error in HC_Ion_Impact_ReLaunch: Invalid original NeutType:",NeutType
          Write (Output_Unit_HC_Alert,*) "Program stopping."
          Stop
       End If


       !
       ! jdemod - update new particle launch for sputtered particle- NOTE: both sputtered and reflected particles are recorded below 
       !          since the launch region has changed. 
       !
       !   
       !HC_Sum_Fragments_Launched (Cur_HC_Spec,Launch_Reg) =  HC_Sum_Fragments_Launched (Cur_HC_Spec,Launch_Reg) + Sput_Weight

 
    ElseIf (Type .eq. "Reflect") Then
       Call HC_Ion_Target_Location (Current_S,Current_Cross,Current_Cell,Current_Ring,Current_R,Current_Z,Target_Index)
       Call HC_Target_Reflect_Species (Cur_HC_Spec,Target_Index,HC_Temperature,Current_Velocity_In_S,Current_Angle)
       Call Get_HC_Isotope_Composition (Cur_HC_Spec,H_Isotope_Composition,Seed,NRand) ! Update new H content from target.
       Sput_Weight = Calc_Sputy (Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition))
       Call HC_Ion_Target_Reflect_Angle (Cur_HC_Spec,Current_Angle,HC_Temperature,Target_Index,Current_Angle,Seed,NRand)
       Call HC_Ion_Target_Reflect_Energy (Cur_HC_Spec,H_Isotope_Composition,Current_Cell,Current_Ring,Current_Angle,Target_Index,&
                           &HC_Temperature)
       Call HC_Ion_Target_Reflect_Velocity (Cur_HC_Spec,H_Isotope_Composition,Current_Cell,Current_Ring,Current_Angle,Target_Index,&
                           &HC_Temperature,Current_Velocity,Current_Velocity_In_R,Current_Velocity_In_Z,Seed,NRand,hc_v)

       ! Reset launch type for reflection event.
       If (NeutType .eq. 2) Then
          ! Originally a target launch.
          If (Launch_Reg .eq. HC_Launch_Reg_Target_1_Dist .or. Launch_Reg .eq. HC_Launch_Reg_Target_1_Pin) Then
             Launch_Reg = HC_Launch_Reg_Refl_Target_1
          ElseIf (Launch_Reg .eq. HC_Launch_Reg_Target_2_Dist .or. Launch_Reg .eq. HC_Launch_Reg_Target_2_Pin) Then
             Launch_Reg = HC_Launch_Reg_Refl_Target_2
          End If
       ElseIf (NeutType .eq. 3) Then
          ! Originally a target self-sputter launch.
       ElseIf (NeutType .eq. 5) Then
          ! Originally a wall launch.
          Launch_Reg = HC_Launch_Reg_Refl_Wall
       ElseIf (NeutType .eq. 6) Then
          ! Originally a 2D neutral launch.
          ! No change.
       ElseIf (NeutType .eq. 7) Then
          ! Originally a reflected ion.
          ! No change.
       elseif (neuttype.eq.0) then 
          ! jdemod
          ! Freespace launch - for puffing ... 
          ! No change - but need to check if neuttype is used in any arrays as an index
       Else
          ! Should not happen for HC launch.
          Write (Output_Unit_HC_Alert,*) "Error in HC_Ion_Impact_ReLaunch: Invalid original NeutType:",NeutType
          Write (Output_Unit_HC_Alert,*) "Program stopping."
          Stop
       End If
 
    End If



    ! New position, energy/velocity, and direction are loaded.  Return to following code.
 
    ! Record sputtered particle relaunch position statistics.
    HC_Total_R_ReProd_Positions (Cur_HC_Spec,Orig_Launch_Reg) =  HC_Total_R_ReProd_Positions (Cur_HC_Spec,Orig_Launch_Reg) + &
                                                                &Current_R * Sput_Weight
    HC_Total_Z_ReProd_Positions (Cur_HC_Spec,Orig_Launch_Reg) =  HC_Total_Z_ReProd_Positions (Cur_HC_Spec,Orig_Launch_Reg) + &
                                                                &Current_Z * Sput_Weight
    HC_Total_ReProd_Angles (Cur_HC_Spec,Orig_Launch_Reg) =  HC_Total_ReProd_Angles (Cur_HC_Spec,Orig_Launch_Reg) + (Current_Angle) &
                                                          &* Sput_Weight
    HC_Sum_Fragments_ReLaunched (Cur_HC_Spec,gwallindex (Target_Index),Orig_Launch_Reg) =  HC_Sum_Fragments_ReLaunched (&
                                                        &Cur_HC_Spec,gwallindex (Target_Index),Orig_Launch_Reg) + Sput_Weight
 

    ! Sum velocity statistics for initial particle launch.
    ! jdemod - the individual particle velocity_multiplier is not visible in this scope - only the global multipler used 
    !          in certain options - for now - set temp_velocity_multipler to 1.0 
    ! Changed code to use Saved_velocity_mult which is recorded in hc_launch_velocity
    !temp_velocity_multiplier = 1.0
    HC_Total_ReProd_Vels_No_VMult (Cur_HC_Spec,Orig_Launch_Reg) =  HC_Total_ReProd_Vels_No_VMult (Cur_HC_Spec,Orig_Launch_Reg) + &
                                 & Current_Velocity /  Saved_Velocity_Mult * Sput_Weight
    HC_Max_ReProd_Vel_No_VMult (Cur_HC_Spec,Orig_Launch_Reg) = MAX ( HC_Max_ReProd_Vel_No_VMult (Cur_HC_Spec,Orig_Launch_Reg), &
                                                                    &Current_Velocity /  Saved_Velocity_Mult)
    !		 HC_Total_Vel_Ang_Mults (Cur_HC_Spec,Orig_Launch_Reg) =  HC_Total_Vel_Ang_Mults (Cur_HC_Spec,Orig_Launch_Reg) +  Velocity_Multiplier * Sput_Weight
    !		 HC_Max_Vel_Ang_Mults (Cur_HC_Spec,Orig_Launch_Reg) = MAX ( HC_Max_Vel_Ang_Mults (Cur_HC_Spec,Orig_Launch_Reg),  Velocity_Multiplier)

    HC_Total_ReProd_Velocities (Cur_HC_Spec,Orig_Launch_Reg) =  HC_Total_ReProd_Velocities (Cur_HC_Spec,Orig_Launch_Reg) + &
                                                               &Current_Velocity * Sput_Weight
    HC_Max_ReProd_Velocities (Cur_HC_Spec,Orig_Launch_Reg) = MAX ( HC_Max_ReProd_Velocities (Cur_HC_Spec,Orig_Launch_Reg), &
                                                                &Current_Velocity)
    HC_Tot_ReProd_Temperatures (Cur_HC_Spec,Orig_Launch_Reg)=  HC_Tot_ReProd_Temperatures (Cur_HC_Spec,Orig_Launch_Reg) + &
                                                              &HC_Temperature * Sput_Weight
 
    ! Record use of determined launch region for output.
    HC_Region_Used (Launch_Reg) =  HC_Region_Used (Launch_Reg) + 1
 
    ! Record sputtered particle launch statistics for new launch region.
    HC_Total_R_Prod_Positions(Cur_HC_Spec,Launch_Reg)=HC_Total_R_Prod_Positions(Cur_HC_Spec,Launch_Reg)+Current_R*Sput_Weight
    HC_Total_Z_Prod_Positions(Cur_HC_Spec,Launch_Reg)=HC_Total_Z_Prod_Positions(Cur_HC_Spec,Launch_Reg)+Current_Z*Sput_Weight
    HC_Total_S_Prod_Positions(Cur_HC_Spec,Launch_Reg)=HC_Total_S_Prod_Positions(Cur_HC_Spec,Launch_Reg)+Current_S*Sput_Weight
    HC_Total_Cross_Prod_Positions (Cur_HC_Spec,Launch_Reg) =  HC_Total_Cross_Prod_Positions (Cur_HC_Spec,Launch_Reg) + &
                                                            &Current_Cross * Sput_Weight
    HC_Total_Prod_Angles (Cur_HC_Spec,Launch_Reg) =  HC_Total_Prod_Angles (Cur_HC_Spec,Launch_Reg) + (Current_Angle) * Sput_Weight


    !
    ! jdemod - I think the following line was a bug - it is recording the launch fragment in the launch array - it has already been
    !          recorded in the relaunch array above - if it is recorded here then the particle accounting will not be correct since
    !          it is not a primary launched particle. Only SPUTTERED particles should increase the base particle count - not REFLECTED particles
    !          move this line to "SPUTTER" option.
    !
    !          Correction: the code above changes the launch_reg value to one indicating a particle reflection has occurred - it then increments
    !                      the particle count.It will duplicate data that is recorded in the hc_sum_fragments_relaunched array. I will have
    !                      to check and see whether the particle weight is used in the density normalization later in the code. It should not be
    !                      used since the counts resulting from this particle are due to the single primary particle being reflected. 
    !                      The code added to record sputtered particles only has been removed.   
    !
    HC_Sum_Fragments_Launched (Cur_HC_Spec,Launch_Reg) =  HC_Sum_Fragments_Launched (Cur_HC_Spec,Launch_Reg) + Sput_Weight
 
    ! Sum velocity statistics for particle launch in new launch region.
    ! Note:  This data will be added to by following DIVIMP self-sputtered launches.
    ! jdemod -
    ! see above temp_velocity_mult set to 1.0
    ! Changed code to use Saved_velocity_mult which is recorded in hc_launch_velocity
    !
    HC_Total_Prod_Vels_No_VMult (Cur_HC_Spec,Launch_Reg) =  HC_Total_Prod_Vels_No_VMult (Cur_HC_Spec,Launch_Reg) + &
                                                           &Current_Velocity /  Saved_Velocity_Mult * Sput_Weight
    HC_Max_Prod_Vel_No_VMult (Cur_HC_Spec,Launch_Reg) = MAX ( HC_Max_Prod_Vel_No_VMult (Cur_HC_Spec,Launch_Reg), Current_Velocity /&
                                                           &  Saved_Velocity_Mult)

    !		 HC_Total_Vel_Ang_Mults (Cur_HC_Spec,Launch_Reg) =  HC_Total_Vel_Ang_Mults (Cur_HC_Spec,Launch_Reg) +  Velocity_Multiplier * Sput_Weight
    !		 HC_Max_Vel_Ang_Mults (Cur_HC_Spec,Launch_Reg) = MAX ( HC_Max_Vel_Ang_Mults (Cur_HC_Spec,Launch_Reg),  Velocity_Multiplier)

    HC_Total_Prod_Velocities(Cur_HC_Spec,Launch_Reg)=HC_Total_Prod_Velocities(Cur_HC_Spec,Launch_Reg)+Current_Velocity*Sput_Weight
    HC_Max_Prod_Velocities (Cur_HC_Spec,Launch_Reg) = MAX ( HC_Max_Prod_Velocities (Cur_HC_Spec,Launch_Reg), Current_Velocity)
    HC_Tot_Prod_Temperatures(Cur_HC_Spec,Launch_Reg)=HC_Tot_Prod_Temperatures(Cur_HC_Spec,Launch_Reg)+HC_Temperature*Sput_Weight
 
    ! Remember original launch position.
    Relaunch_Cell = Current_Cell
    Relaunch_Ring = Current_Ring
    Relaunch_R = Current_R
    Relaunch_Z = Current_Z
    Relaunch_HC_Species = Cur_HC_Spec
 
    ! Particle initialization.
    HC_Has_Leaked = .False.
    HC_Has_Leaked_Core = .False.
    HC_Leak_Particles (Launch_Reg) = 1 ! CLEAKP
 
    ! Add ion sputter to total removal.
    ! Note:  An ion can only exit the grid by crossing a target segment, so the event should be added to NEROS (Target_Index,3) - total removal.
    HC_Erosion (Target_Index,3,Cur_HC_Spec) =  HC_Erosion (Target_Index,3,Cur_HC_Spec) + Sput_Weight ! NEROS.
 
  End Subroutine HC_Ion_Impact_ReLaunch
 
  Subroutine HC_Neutral_Impact_ReLaunch (Type,Cur_HC_Spec,Last_HC_Species,H_Isotope_Composition,Last_H_Isotope_Composition,&
       &Launch_Reg, &
       & Sput_Weight,Neutral_Time_Step_Count,Ion_Time_Step_Count,Eq_Total_Ion_Time_Steps,Reflection_Counter,NeutType,Wall_Index,&
       &Reflection_Angle,Segment_Normal_Angle,Current_Theta,S_Start, &
       & RNew,ZNew,Current_R,Current_Z,HC_Temperature,Current_Velocity,Current_Angle,Seed,NRand,Current_Cell,Current_Ring, &
       & Current_Velocity_In_R,Current_Velocity_In_Z,Current_S,Current_Cross,hc_v)
    ! Particle is reflected or sputtered back into plasma.  Find new species and update all statistics.
 
    ! Required modules.
    Use ComHC ! Preset variables.
    Use HC_Init_DIV_Data ! Gain access to cell/global/geometry/current data structures.
    Use HC_Init_DIV_Diag ! Data diagnostics.
    Use HC_Init_Lib_Data ! Hydrocarbon data and utilities.
    Use HC_Utilities ! General purpose utilities.
    Use HC_Get ! gwallindex.
    use hc_velocity_type ! HC velocity type declaration
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables
    Character (len=*), Intent (In) :: Type
    Integer, Intent (InOut) :: Cur_HC_Spec
    Integer, Intent (Out) :: Last_HC_Species
    Integer, Dimension (Number_H_Species), Intent (InOut) :: H_Isotope_Composition
    Integer, Dimension (Number_H_Species), Intent (Out) :: Last_H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (InOut) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (Out) :: Last_H_Isotope_Composition
    Integer, Intent (InOut) :: Launch_Reg
    Real, Intent (InOut) :: Sput_Weight ! SPUTY
    Double Precision, Intent (InOut) :: Neutral_Time_Step_Count
    Double Precision, Intent (InOut) :: Ion_Time_Step_Count
    Double Precision, Intent (InOut) :: Eq_Total_Ion_Time_Steps
    Integer, Intent (InOut) :: Reflection_Counter
    Integer, Intent (In) :: NeutType
    Integer, Intent (In) :: Wall_Index
    Real, Intent (Out) :: Reflection_Angle
    Real, Intent (IN) :: Segment_Normal_Angle
    Real, Intent (InOut) :: Current_Theta
    Real, Intent (InOut) :: S_Start
    Real, Intent (InOut) :: RNew
    Real, Intent (InOut) :: ZNew
    Real, Intent (InOut) :: Current_R
    Real, Intent (InOut) :: Current_Z
    Real, Intent (InOut) :: HC_Temperature
    Real, Intent (InOut) :: Current_Velocity
    Real, Intent (InOut) :: Current_Angle
    Double Precision, Intent (In) :: Seed
    Integer, Intent (InOut) :: NRand
    Integer, Intent (InOut) :: Current_Cell
    Integer, Intent (InOut) :: Current_Ring
    Real, Intent (InOut) :: Current_Velocity_In_R
    Real, Intent (InOut) :: Current_Velocity_In_Z
    Real, Intent (InOut) :: Current_S
    Real, Intent (InOut) :: Current_Cross

    type(hc_velocity_type1) :: hc_v ! HC velocity structure
 
    ! Define local variables.
    Integer :: Orig_Launch_Reg
    
    ! Record previous HC data.
    Last_HC_Species = Cur_HC_Spec
    Last_H_Isotope_Composition = H_Isotope_Composition
    Orig_Launch_Reg = Launch_Reg
 
    ! Check if particle struck target.
    If (INT (gwallpt (Wall_Index,18)) .gt. 0) Then
       ! Record collision with target.
       HC_Num_Striking_Target (Cur_HC_Spec,Launch_Reg) =  HC_Num_Striking_Target (Cur_HC_Spec,Launch_Reg) + Sput_Weight
    Else
       ! Record collision with wall.
       HC_Num_Reach_Wall (Cur_HC_Spec,Launch_Reg) =  HC_Num_Reach_Wall (Cur_HC_Spec,Launch_Reg) + Sput_Weight
    End If
    Reflection_Counter = Reflection_Counter + 1
    Total_HC_Reflections (Cur_HC_Spec,Launch_Reg) =  Total_HC_Reflections (Cur_HC_Spec,Launch_Reg) + 1
    Max_HC_Reflections_Found ( Launch_HC_Species,Launch_Reg) = MAX0 (Reflection_Counter, Max_HC_Reflections_Found ( &
                              &Launch_HC_Species,Launch_Reg))
 
    ! Rerelease-has-occurred indicator reset to false for each particle launched.
    HC_Hit_Vessel = .True.
 
    ! Reset time for new particle launch.
    Neutral_Time_Step_Count = 0.0
    Ion_Time_Step_Count = 0.0
    Eq_Total_Ion_Time_Steps = 0.0
 
    ! Newly launched particle properties.
    If (Type .eq. "Sputter") Then

       Call HC_Neut_Vessel_Sputter_Location (RNew,ZNew,Current_R,Current_Z,Current_Cell,Current_Ring,Current_S,Current_Cross,&
                                            &Current_Theta,S_Start)
       Call HC_Vessel_Sputter_Species (Cur_HC_Spec,H_Isotope_Composition,Wall_Index,Current_Theta,Current_R,Current_Z,&
                                      &HC_Temperature,Current_Velocity,0.0,Current_Angle)
       Call Get_HC_Isotope_Composition (Cur_HC_Spec,H_Isotope_Composition,Seed,NRand) ! Update new H content from wall.
       Sput_Weight = Calc_Sputy (Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition))
       Call HC_Neut_Vessel_Sputter_Angle (Cur_HC_Spec,Current_Angle,HC_Temperature,Wall_Index,Reflection_Angle,&
                                         &Segment_Normal_Angle,Current_Angle,Seed,NRand)
       
       Call HC_Neut_Vessel_Sputter_Energy (Cur_HC_Spec,H_Isotope_Composition,Wall_Index,Current_Velocity,HC_Temperature)
       Call HC_Neut_Vessel_Sputter_Velocity (Cur_HC_Spec,H_Isotope_Composition,Current_Cell,Current_Ring,Current_Angle,Wall_Index,&
                                            &HC_Temperature,Current_Velocity,Current_Velocity_In_R,Current_Velocity_In_Z,&
                                            &Seed,NRand,hc_v)
 
       ! Reset launch type for sputtering event.
       If (NeutType .eq. 2) Then
          ! Originally a target launch.
          If (Launch_Reg .eq. HC_Launch_Reg_Target_1_Dist .or. Launch_Reg .eq. HC_Launch_Reg_Target_1_Pin) Then
             Launch_Reg = HC_Launch_Reg_Sput_Target_1
          ElseIf (Launch_Reg .eq. HC_Launch_Reg_Target_2_Dist .or. Launch_Reg .eq. HC_Launch_Reg_Target_2_Pin) Then
             Launch_Reg = HC_Launch_Reg_Sput_Target_2
          End If
       ElseIf (NeutType .eq. 3) Then
          ! Originally a target self-sputter launch.
       ElseIf (NeutType .eq. 5) Then
          ! Originally a wall launch.
          Launch_Reg = HC_Launch_Reg_Sput_Wall
       ElseIf (NeutType .eq. 6) Then
          ! Originally a 2D neutral launch.
          ! No change.
       ElseIf (NeutType .eq. 7) Then
          ! Originally a reflected ion.
          ! No change.
       elseif (neuttype.eq.0) then 
          ! jdemod
          ! Freespace launch - for puffing ... 
          ! No change - but need to check if neuttype is used in any arrays as an index
       Else
          ! Should not happen for HC launch.
          Write (Output_Unit_HC_Alert,*) "Error in HC_Neutral_Impact_ReLaunch: Invalid original NeutType:",NeutType
          Write (Output_Unit_HC_Alert,*) "Program stopping."
          Stop
       End If


       !
       ! jdemod - update new particle launch for sputtered particle - see NOTES for ions - since launch region is changed this 
       !          is now done below for both sputtered and reflected particles. 
       !
       !
       !HC_Sum_Fragments_Launched (Cur_HC_Spec,Launch_Reg) =  HC_Sum_Fragments_Launched (Cur_HC_Spec,Launch_Reg) + Sput_Weight

 
    ElseIf (Type .eq. "Reflect") Then

       Call HC_Neut_Vessel_Reflect_Location (RNew,ZNew,Current_R,Current_Z,Current_Cell,Current_Ring,Current_S,Current_Cross,&
                                            &Current_Theta,S_Start)
       Call HC_Vessel_Reflect_Species (Cur_HC_Spec,Wall_Index,HC_Temperature,Current_Velocity,0.0,Current_Angle)
       Call Get_HC_Isotope_Composition (Cur_HC_Spec,H_Isotope_Composition,Seed,NRand) ! Update new H content from wall.
       Sput_Weight = Calc_Sputy (Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition))
       Call HC_Neut_Vessel_Reflect_Angle (Cur_HC_Spec,Current_Angle,HC_Temperature,Wall_Index,Reflection_Angle,&
                                         &Segment_Normal_Angle,Current_Angle,Seed,NRand)
       Call HC_Neut_Vessel_Reflect_Energy (Cur_HC_Spec,H_Isotope_Composition,Wall_Index,HC_Temperature)

       Call HC_Neut_Vessel_Reflect_Velocity (Cur_HC_Spec,H_Isotope_Composition,Current_Cell,Current_Ring,Current_Angle,Wall_Index,&
                                            &HC_Temperature,Current_Velocity,Current_Velocity_In_R,Current_Velocity_In_Z,&
                                            &Seed,NRand,hc_v)

       ! Reset launch type for reflection event.
       If (NeutType .eq. 2) Then
          ! Originally a target launch.
          If (Launch_Reg .eq. HC_Launch_Reg_Target_1_Dist .or. Launch_Reg .eq. HC_Launch_Reg_Target_1_Pin) Then
             Launch_Reg = HC_Launch_Reg_Refl_Target_1
          ElseIf (Launch_Reg .eq. HC_Launch_Reg_Target_2_Dist .or. Launch_Reg .eq. HC_Launch_Reg_Target_2_Pin) Then
             Launch_Reg = HC_Launch_Reg_Refl_Target_2
          End If
       ElseIf (NeutType .eq. 3) Then
          ! Originally a target self-sputter launch.
       ElseIf (NeutType .eq. 5) Then
          ! Originally a wall launch.
          Launch_Reg = HC_Launch_Reg_Refl_Wall
       ElseIf (NeutType .eq. 6) Then
          ! Originally a 2D neutral launch.
          ! No change.
       ElseIf (NeutType .eq. 7) Then
          ! Originally a reflected ion.
          ! No change.
       elseif (neuttype.eq.0) then 
          ! jdemod
          ! Freespace launch - for puffing ... 
          ! No change - but need to check if neuttype is used in any arrays as an index
       Else
          ! Should not happen for HC launch.
          Write (Output_Unit_HC_Alert,*) "Error in HC_Outside_Neutral: Invalid original NeutType:",NeutType
          Write (Output_Unit_HC_Alert,*) "Program stopping."
          Stop
       End If


 
    End If


    ! New position, energy/velocity, and direction are loaded.  Return to following code.
 
    ! Record sputtered particle relaunch position statistics.
    HC_Total_R_ReProd_Positions (Cur_HC_Spec,Orig_Launch_Reg) =  HC_Total_R_ReProd_Positions (Cur_HC_Spec,Orig_Launch_Reg) + &
                                                                &Current_R * Sput_Weight
    HC_Total_Z_ReProd_Positions (Cur_HC_Spec,Orig_Launch_Reg) =  HC_Total_Z_ReProd_Positions (Cur_HC_Spec,Orig_Launch_Reg) + &
                                                                &Current_Z * Sput_Weight
    HC_Total_ReProd_Angles (Cur_HC_Spec,Orig_Launch_Reg) =  HC_Total_ReProd_Angles (Cur_HC_Spec,Orig_Launch_Reg) + (Current_Angle) &
                                                        &* Sput_Weight
    HC_Sum_Fragments_ReLaunched (Cur_HC_Spec,Wall_Index,Orig_Launch_Reg) =  HC_Sum_Fragments_ReLaunched (Cur_HC_Spec,Wall_Index,&
                                                                           &Orig_Launch_Reg) + Sput_Weight
 
    ! Sum velocity statistics for initial particle launch.
    ! jdemod - change to Saved_velocity_mult from velocity_multiplier
    HC_Total_ReProd_Vels_No_VMult (Cur_HC_Spec,Orig_Launch_Reg) =  HC_Total_ReProd_Vels_No_VMult (Cur_HC_Spec,Orig_Launch_Reg) + &
                                                                 & Current_Velocity /  Saved_Velocity_Mult * Sput_Weight
    HC_Max_ReProd_Vel_No_VMult (Cur_HC_Spec,Orig_Launch_Reg) = MAX ( HC_Max_ReProd_Vel_No_VMult (Cur_HC_Spec,Orig_Launch_Reg), &
                                                             &Current_Velocity /  Saved_Velocity_Mult)
    !		 HC_Total_Vel_Ang_Mults (Cur_HC_Spec,Orig_Launch_Reg) =  HC_Total_Vel_Ang_Mults (Cur_HC_Spec,Orig_Launch_Reg) +  Velocity_Multiplier * Sput_Weight
    !		 HC_Max_Vel_Ang_Mults (Cur_HC_Spec,Orig_Launch_Reg) = MAX ( HC_Max_Vel_Ang_Mults (Cur_HC_Spec,Orig_Launch_Reg),  Velocity_Multiplier)
    HC_Total_ReProd_Velocities (Cur_HC_Spec,Orig_Launch_Reg) =  HC_Total_ReProd_Velocities (Cur_HC_Spec,Orig_Launch_Reg) + &
                                                              &Current_Velocity * Sput_Weight
    HC_Max_ReProd_Velocities (Cur_HC_Spec,Orig_Launch_Reg) = MAX ( HC_Max_ReProd_Velocities (Cur_HC_Spec,Orig_Launch_Reg), &
                                                                  &Current_Velocity)
    HC_Tot_ReProd_Temperatures (Cur_HC_Spec,Orig_Launch_Reg)=  HC_Tot_ReProd_Temperatures (Cur_HC_Spec,Orig_Launch_Reg) + &
                                                                    &HC_Temperature * Sput_Weight
 
    ! Record use of determined launch region for output.
    HC_Region_Used (Launch_Reg) =  HC_Region_Used (Launch_Reg) + 1
 
    ! Record sputtered particle launch statistics for new launch region.
    HC_Total_R_Prod_Positions(Cur_HC_Spec,Launch_Reg)=HC_Total_R_Prod_Positions(Cur_HC_Spec,Launch_Reg)+Current_R*Sput_Weight
    HC_Total_Z_Prod_Positions(Cur_HC_Spec,Launch_Reg)=HC_Total_Z_Prod_Positions(Cur_HC_Spec,Launch_Reg)+Current_Z*Sput_Weight
    HC_Total_S_Prod_Positions(Cur_HC_Spec,Launch_Reg)=HC_Total_S_Prod_Positions(Cur_HC_Spec,Launch_Reg)+Current_S*Sput_Weight
    HC_Total_Cross_Prod_Positions (Cur_HC_Spec,Launch_Reg) =  HC_Total_Cross_Prod_Positions (Cur_HC_Spec,Launch_Reg) + &
                                                                &Current_Cross * Sput_Weight
    HC_Total_Prod_Angles (Cur_HC_Spec,Launch_Reg) =  HC_Total_Prod_Angles (Cur_HC_Spec,Launch_Reg) + (Current_Angle) * Sput_Weight

    !
    ! jdemod - I think the following line was a bug - it is recording the launch fragment in the launch array - it has already been
    !          recorded in the relaunch array above - if it is recorded here then the particle accounting will not be correct since
    !          it is not a primary launched particle.
    !
    !          Correction: the code above changes the launch_reg value to one indicating a particle reflection has occurred - it then increments
    !                      the particle count.It will duplicate data that is recorded in the hc_sum_fragments_relaunched array. I will have
    !                      to check and see whether the particle weight is used in the density normalization later in the code. It should not be
    !                      used since the counts resulting from this particle are due to the single primary particle being reflected. 
    !                      The code added to record sputtered particles only has been removed.   
    !
    HC_Sum_Fragments_Launched (Cur_HC_Spec,Launch_Reg) =  HC_Sum_Fragments_Launched (Cur_HC_Spec,Launch_Reg) + Sput_Weight
 
    ! Sum velocity statistics for particle launch in new launch region.
    ! Note:  This data will be added to by following DIVIMP self-sputtered launches.
    ! jdemod - change to Saved_velocity_mult from velocity_multiplier
    HC_Total_Prod_Vels_No_VMult (Cur_HC_Spec,Launch_Reg) =  HC_Total_Prod_Vels_No_VMult (Cur_HC_Spec,Launch_Reg) + &
                                                           &Current_Velocity /  Saved_Velocity_Mult * Sput_Weight
    HC_Max_Prod_Vel_No_VMult (Cur_HC_Spec,Launch_Reg) = MAX ( HC_Max_Prod_Vel_No_VMult (Cur_HC_Spec,Launch_Reg), Current_Velocity /&
                                                          &  Saved_Velocity_Mult)
    !		 HC_Total_Vel_Ang_Mults (Cur_HC_Spec,Launch_Reg) =  HC_Total_Vel_Ang_Mults (Cur_HC_Spec,Launch_Reg) +  Velocity_Multiplier * Sput_Weight
    !		 HC_Max_Vel_Ang_Mults (Cur_HC_Spec,Launch_Reg) = MAX ( HC_Max_Vel_Ang_Mults (Cur_HC_Spec,Launch_Reg),  Velocity_Multiplier)
    HC_Total_Prod_Velocities(Cur_HC_Spec,Launch_Reg)=HC_Total_Prod_Velocities(Cur_HC_Spec,Launch_Reg)+Current_Velocity*Sput_Weight
    HC_Max_Prod_Velocities (Cur_HC_Spec,Launch_Reg) = MAX ( HC_Max_Prod_Velocities (Cur_HC_Spec,Launch_Reg), Current_Velocity)
    HC_Tot_Prod_Temperatures(Cur_HC_Spec,Launch_Reg)=HC_Tot_Prod_Temperatures(Cur_HC_Spec,Launch_Reg)+HC_Temperature*Sput_Weight
 
    ! Remember original launch position.
    Relaunch_Cell = Current_Cell
    Relaunch_Ring = Current_Ring
    Relaunch_R = Current_R
    Relaunch_Z = Current_Z
    Relaunch_HC_Species = Cur_HC_Spec
 
    ! Particle initialization.
    HC_Has_Leaked = .False.
    HC_Has_Leaked_Core = .False.
    HC_Leak_Particles (Launch_Reg) = 1 ! CLEAKP
 
    ! Modify commons for cell tau parallel, stopping, and heating for current cell and ring.
    Call Initialize_State_Prop_Data (Current_Cell,Current_Ring,Get_HC_Charge (Cur_HC_Spec))
    Call Modify_Taus (Cur_HC_Spec,H_Isotope_Composition,Sput_Weight,Current_Cell,Current_Ring)
    ! jdemod
    !Call Modify_Taus (Cur_HC_Spec,H_Isotope_Composition,Sput_Weight,HC_Temperature,Current_Cell,Current_Ring,Current_S,&
    !                             &Seed,Random_Numbers_Used,NRand)
 
    ! Add neutral sputter to total removal (NEROS(x,3)).
    If (INT (gwallpt (Wall_Index,18)) .gt. 0) Then
       HC_Erosion(INT(gwallpt(Wall_Index,18)),3,Cur_HC_Spec)=HC_Erosion(INT(gwallpt(Wall_Index,18)),3,Cur_HC_Spec)+Sput_Weight! NEROS.
    End If

 
  End Subroutine HC_Neutral_Impact_ReLaunch
 
  Subroutine HC_Vessel_Sputter_Species (HC_Species,H_Isotope_Composition,Vessel_Segment, &
       & Current_Theta,R_Position,Z_Position,HC_Temperature,HC_Velocity,HC_Velocity_In_S,HC_Angle)
 
    ! Required modules.
    Use HC_Init_Lib_Data ! Gain access to state and transition tables.
    Use HC_Stack ! Gain access to hydrocarbon storage stack and operations.
    Use ComHC ! Access to global variables related to the hydrocarbon module.
    Use HC_Utilities ! Access hydrocarbon utilities.
    Use HC_Get
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables
    Integer, Intent (InOut) :: HC_Species ! Note, only one species passed in and out, additional species loaded in the HC stack.
    Integer, Dimension (Number_H_Species), Intent (In) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (In) :: H_Isotope_Composition
    Integer, Intent (In) :: Vessel_Segment ! ID
    Real, Intent (In) :: Current_Theta ! THETA
    Real, Intent (In) :: R_Position ! R
    Real, Intent (In) :: Z_Position ! Z
    Real, Intent (InOut) :: HC_Temperature ! TEMI
    Real, Intent (InOut) :: HC_Velocity ! VEL
    Real, Intent (In) :: HC_Velocity_In_S ! VEL
    Real, Intent (InOut) :: HC_Angle
 
    ! Declare local variables.
    Integer :: i
 
    ! Proceed to sputter particle.
    ! First, find new hydrocarbon species (H28).
    If (hc_sputtering_species_model .eq. 0) Then
       ! Use preset sputtering table to decide which single hydrocarbon species is emitted.
       ! Note, this may change to multiple Carbon atoms if the primary channel of sputtering for C2 or higher is changed.
       HC_Species = HC_Sputtering_Table (HC_Species) % Sputtered_State_Number (1)
       Do i = 2,Highest_Carbon_Content,1 ! Up to 3 for C3Hx series of hydrocarbons.
          If (HC_Sputtering_Table (HC_Species) % Sputtered_State_Number (i) .ne. 0) Then
             ! Multiple sputtered hydocarbon species.  Store additional data on HC stack.
             ! Note that current S and Cross are given values of 0.0 because particles are always sputtered from the wall as neutrals.
             Call HC_Push (HC_Species,R_Position,Z_Position,0.0,0.0,HC_Angle,Current_Theta,HC_Velocity,HC_Velocity_In_S,&
                          &HC_Temperature)
          End If
       End Do
 
    ElseIf (hc_sputtering_species_model .eq. 1) Then
       ! Use Alman and Ruzic sputtering data.  Note the possibility of multiple hydrocarbons released and returned in HC_Species_Out.
       !Call Sputter_Species_Alman_Ruzic (HC_Species,Vessel_Segment,R_Position,Z_Position,HC_Temperature,HC_Velocity,HC_Angle_To_Normal,Surface_Temperature)
    ElseIf (hc_sputtering_species_model .eq. 2) Then
       ! Sputtered species equals incoming species. Do nothing.
 
    Else
       ! Unsupported sputtering model selected.
       Write (Output_Unit_HC_Alert,*) "Unsupported sputtering model selected:", hc_sputtering_species_model
       Write (Output_Unit_HC_Alert,*) "Program Stopping."
       Stop
    End If
 
  End Subroutine HC_Vessel_Sputter_Species
 
  Subroutine HC_Neut_Vessel_Sputter_Angle (HC_Species,Incoming_Angle, &
       & Incoming_Energy,Vessel_Segment,Sputtering_Angle,Segment_Normal_Angle,HC_Angle,Seed,NRand)
 
    ! Required modules.
    Use ComHC ! Access to global variables related to the hydrocarbon module.
    Use HC_Init_Div_Data ! Access to misc_data_table.
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables.
    Integer, Intent (In) :: HC_Species !
    Real, Intent (In) :: Incoming_Angle
    Real, Intent (In) :: Incoming_Energy ! TEMI
    Integer, Intent (In) :: Vessel_Segment ! ID
    Real, Intent (In) :: Sputtering_Angle ! TREF
    Real, Intent (In) :: Segment_Normal_Angle ! TNORM
    Real, Intent (Out) :: HC_Angle
    Double Precision, Intent (In) :: Seed
    Integer, Intent (InOut) :: NRand
 
    ! Declare local variables.
    Real :: HC_Angle_To_Normal
    Real :: Random_Value_1
    Real :: Random_Value_2
    Real :: Beta, Psi


    ! Find angle to normal given vessel segment and impact angle.
    HC_Angle_To_Normal = ABS (-Segment_Normal_Angle - Incoming_Angle)
    !write(0,*) "howyah",HC_Angle_To_Normal,hc_sputtering_angle_model
    ! Find direction of particle after reflection (H31).

    ! jdemod 
    ! This is being hard coded here to 5!?
    ! remove this since hc_sputtering_angle_model is input H57. It defaults to the value specified for DIVIMP particles
    ! but a different option can be selected for hydrocarbons
    ! hc_sputtering_angle_model = 5


    If (hc_sputtering_angle_model .eq. 0) Then
       ! This option should never happen.
       Write (Output_Unit_HC_Alert,*) "Error: Should not have entered function neut vessel hc_sputter for no sputter condition: &
&HC_Neut_Vessel_Sputter_Angle.  Note that the Neutral Wall Ref option should be set to 4."
       Write (Output_Unit_HC_Alert,*) "Program stopping."
       Stop
 
    ElseIf (hc_sputtering_angle_model .ge. 1 .and. hc_sputtering_angle_model .le. 4) Then
       !		ElseIf (hc_sputtering_angle_model .ge. 1 .or. hc_sputtering_angle_model .le. 4) Then
       ! Angle is supported by main REFANG DIVIMP code.
       ! =1, Specular sputtering from wall, bounds of +/-90 degrees from normal to vessel wall.
       ! =2, Isotropic sputtering with ACOS(x) distribution from vessel normal, bounds of +/-90 degrees from normal to vessel wall.
       ! =3, Isotropic sputtering with ASIN(x) distribution from vessel normal, bounds of +/-90 degrees from normal to vessel wall.
       ! =4, Isotropic sputtering with ASIN(SQRT(x)) distribution from vessel normal, bounds of +/-90 degrees from normal to vessel wall.
       HC_Angle = Sputtering_Angle
 
    ElseIf (hc_sputtering_angle_model .eq. 5) Then
 
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value_1)
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value_2)
 
       Beta = ASIN (SQRT(Random_Value_1))
       Psi = 2.0 *  Pi_Value * Random_Value_2

       ! jdemod - this formula is using "sputtering_angle" to calculate hc_angle. However, this value is found from the reflection_angle
       !          of the impacting species and gives a distribution around that which could potentially give an hc_angle that will point outside the
       !          vessel. I think this should be using segment_normal_angle instead
       !HC_Angle = Sputtering_Angle + ATAN (TAN (Beta) * COS (Psi))
       HC_Angle = Segment_normal_angle + ATAN (TAN (Beta) * COS (Psi))
       !write(0,*) "here",HC_Angle,Sputtering_Angle
    ElseIf (hc_sputtering_angle_model .eq. 10) Then
       ! Normal reflection from vessel wall
       HC_Angle = Segment_Normal_Angle
 
    ElseIf (hc_sputtering_angle_model .eq. 11) Then
       ! Alman and Ruzic angle data used.  Ignore Reflection_Angle (TREF) returned by INTCALC.
       !Call Sputter_Angle_Alman_Ruzic (HC_Species,Vessel_Segment,R_Position,Z_Position,HC_Temperature,HC_Velocity,HC_Angle_To_Normal,Surface_Temperature)
 
    Else		
       ! Unsupported reflection model selected.
       Write (Output_Unit_HC_Alert,*) "Unsupported reflection model selected: ", hc_sputtering_angle_model
       Write (Output_Unit_HC_Alert,*) "Program Stopping."
       Stop
    End If

  End Subroutine HC_Neut_Vessel_Sputter_Angle
 
  Subroutine HC_Ion_Target_Sputter_Angle (HC_Species,Incoming_Angle,Incoming_Energy,Target_Segment,Sputtering_Angle,Seed,NRand)
 
    ! Required modules.
    Use ComHC ! Access to global variables related to the hydrocarbon module.
    Use HC_Init_DIV_Data
    Use HC_Get
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables.
    Integer, Intent (In) :: HC_Species !
    Real, Intent (In) :: Incoming_Angle
    Real, Intent (In) :: Incoming_Energy ! TEMI
    Integer, Intent (In) :: Target_Segment ! ID
    Real, Intent (Out) :: Sputtering_Angle ! TREF
    Double Precision, Intent (In) :: Seed
    Integer, Intent (InOut) :: NRand
 
    ! Local variables.
    Real :: Random_Value_1
    Real :: Random_Value_2
 
    ! Declare local variables.
    Real :: HC_Angle_To_Normal
    Real :: Beta, Psi
 
    ! Find angle to normal given vessel segment and impact angle.
    HC_Angle_To_Normal = ABS (-gthetas (Target_Segment) - Incoming_Angle)
 
    ! Find direction of particle after reflection (H31).
    If (hc_sputtering_angle_model .eq. 0) Then
       ! This option should never happen.
       Write (Output_Unit_HC_Alert,*) "Error: Should not have entered function ion target hc_sputter for no sputter condition: &
&HC_Ion_Target_Sputter_Angle."
       Write (Output_Unit_HC_Alert,*) "Program Stopping."
       Stop
    ElseIf (hc_sputtering_angle_model .eq. 1) Then
       ! =1, Specular sputtering from wall, bounds of +/-90 degrees from normal to vessel wall.
       Write (Output_Unit_HC_Alert,*) "Specular ion sputtering not currently supported."
       Write (Output_Unit_HC_Alert,*) "Program Stopping."
       Stop
 
    ElseIf (hc_sputtering_angle_model .eq. 2) Then
       ! =2, Isotropic sputtering SINE distribution from vessel normal, bounds of +/-90 degrees from normal to vessel wall.
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value_1)
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value_2)
       Sputtering_Angle = gthetas (Target_Segment) + SIGN (ACOS (Random_Value_1), Random_Value_2 - 0.5)
 
    ElseIf (hc_sputtering_angle_model .eq. 3) Then
       ! =3, Isotropic sputtering COS distribution from vessel normal, bounds of +/-90 degrees from normal to vessel wall.
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value_1)
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value_2)
       Sputtering_Angle = gthetas (Target_Segment) + SIGN (ASIN (Random_Value_1), Random_Value_2 - 0.5)
 
    ElseIf (hc_sputtering_angle_model .eq. 4) Then
       ! =4, Isotropic sputtering COS distribution from vessel normal, bounds of +/-90 degrees from normal to vessel wall.
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value_1)
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value_2)
       Sputtering_Angle = gthetas (Target_Segment) + SIGN (ASIN (SQRT (Random_Value_1)), Random_Value_2 - 0.5)
 
    ElseIf (hc_sputtering_angle_model .eq. 5) Then
       ! =5, Isotropic sputtering COS distribution from vessel normal, bounds of +/-90 degrees from normal to vessel wall.
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value_1)
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value_2)
 
       Beta = ASIN (SQRT(Random_Value_1))
       Psi = 2.0 *  Pi_Value * Random_Value_2
       Sputtering_Angle = gthetas (Target_Segment) + ATAN (TAN (Beta) * COS (Psi))
 
    ElseIf (hc_sputtering_angle_model .eq. 10) Then
       ! Normal reflection from vessel wall
       Sputtering_Angle = gthetas (Target_Segment)
 
    ElseIf (hc_sputtering_angle_model .eq. 11) Then
       ! Alman and Ruzic angle data used.  Ignore Reflection_Angle (TREF) returned by INTCALC.
       !Call Sputter_Angle_Alman_Ruzic (HC_Species,Vessel_Segment,R_Position,Z_Position,HC_Temperature,HC_Velocity,HC_Angle_To_Normal,Surface_Temperature)
 
 
    Else		
       ! Unsupported reflection model selected.
       Write (Output_Unit_HC_Alert,*) "Unsupported reflection model selected: ", hc_sputtering_angle_model
       Write (Output_Unit_HC_Alert,*) "Program Stopping."	
       Stop
    End If
    !Write (Output_Unit_Scratch,*) "hc_ion_target_sputter_angle",Target_Segment,Incoming_Angle,Sputtering_Angle,gthetas(Target_Segment)
 
  End Subroutine HC_Ion_Target_Sputter_Angle
 
  Subroutine HC_Neut_Vessel_Sputter_Energy (HC_Species,H_Isotope_Composition,Vessel_Segment,HC_Velocity,Sputtered_Energy)
 
    ! Required modules.
    Use ComHC ! Access to global variables related to the hydrocarbon module.
    Use HC_Init_Lib_Data
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables.
    Integer, Intent (In) :: HC_Species !
    Integer, Dimension (Number_H_Species), Intent (In) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (In) :: H_Isotope_Composition
    Integer, Intent (In) :: Vessel_Segment ! IN
    Real, Intent (In) :: HC_Velocity ! VIN
    Real, Intent (Out) :: Sputtered_Energy ! TEMN
 
    !Sputtered_Energy = Find_HC_Mass (HC_Species,H_Isotope_Composition) * (HC_Velocity / 1.38E4) * (HC_Velocity / 1.38E4)
    Sputtered_Energy = hc_sput_energy_neutral_preset
 
  End Subroutine HC_Neut_Vessel_Sputter_Energy
 
  Subroutine HC_Ion_Target_Sputter_Energy (HC_Species,H_Isotope_Composition,Current_Cell,Current_Ring,Current_Angle,Target_Index,&
&Sputtered_Energy)
 
    ! Required modules.
    Use ComHC ! Access to global variables related to the hydrocarbon module.
    Use HC_Utilities
    Use HC_Init_Lib_Data
    Use HC_Init_DIV_Data
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables.
    Integer, Intent (In) :: HC_Species
    Integer, Dimension (Number_H_Species), Intent (In) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (In) :: H_Isotope_Composition
    Integer, Intent (In) :: Current_Cell
    Integer, Intent (In) :: Current_Ring
    Real, Intent (In) :: Current_Angle
    Integer, Intent (In) :: Target_Index ! IN
    Real, Intent (Out) :: Sputtered_Energy ! TEMN
 
    ! Declare local variables.
    Real :: HC_Velocity ! VIN, Used only for thermal temperature emission.
 
    ! Find energy and velocity of particle after reflection (H29, H30).
    If (hc_sputtering_energy_model .eq. 0) Then
       ! Use preset energy and calculate velocity.
       Sputtered_Energy = hc_sput_energy_ion_preset
 
    ElseIf (hc_sputtering_energy_model .eq. 1) Then
       ! Do not change energy at all - reflected particle energy is same as impact.
       ! Recalculate velocity for potential new mass.
 
    ElseIf (hc_sputtering_energy_model .eq. 2) Then
       ! Assign new temperature as well.
       HC_Velocity = SQRT (8.0 * 1.38066E-23 * Find_Target_Temperature (Target_Index) / ( Pi_Value * Find_HC_Mass (HC_Species,&
                          &H_Isotope_Composition)))
       Sputtered_Energy = Find_HC_Mass (HC_Species,H_Isotope_Composition) * (HC_Velocity / 1.38E4) * (HC_Velocity / 1.38E4)
 
    ElseIf (hc_sputtering_energy_model .eq. 3) Then
       ! Particle energy found with Alman & Ruzic data.
       !Call Sputter_EnergyVel_Alman_Ruzic (HC_Species,Target_Index,R_Position,Z_Position,HC_Temperature,HC_Velocity,HC_Angle_To_Normal,Surface_Temperature)
 
    Else		
       ! Unsupported reflection model selected.
       Write (Output_Unit_HC_Alert,*) "Unsupported sputtering energy/vel model selected: ", hc_sputtering_energy_model
       Write (Output_Unit_HC_Alert,*) "Program Stopping."		
       Stop
    End If
 
  End Subroutine HC_Ion_Target_Sputter_Energy
 
  Subroutine HC_Neut_Vessel_Sputter_Velocity (HC_Species,H_Isotope_Composition,Current_Cell,Current_Ring, &
       & Current_Angle,Vessel_Segment,HC_Temperature,HC_Velocity,Current_Velocity_In_R,Current_Velocity_In_Z,&
       & Seed,NRand,hc_v)
 
    ! Required modules.
    Use ComHC ! Access to global variables related to the hydrocarbon module.
    Use HC_Init_Lib_Data
    Use HC_Init_DIV_Data
    Use HC_Utilities ! MB velocity distribution.
    Use HC_Get
    use hc_kinetics
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables.
    Integer, Intent (In) :: HC_Species !
    Integer, Dimension (Number_H_Species), Intent (In) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (In) :: H_Isotope_Composition
    Integer, Intent (In) :: Current_Cell ! IK
    Integer, Intent (In) :: Current_Ring ! IR
    Real, Intent (In) :: Current_Angle
    Integer, Intent (In) :: Vessel_Segment ! IN
    Real, Intent (In) :: HC_Temperature
    Real, Intent (Out) :: HC_Velocity ! VIN
    Real, Intent (Out) :: Current_Velocity_In_R
    Real, Intent (Out) :: Current_Velocity_In_Z
    Double Precision, Intent (In) :: Seed
    Integer, Intent (InOut) :: NRand
 
    type(hc_velocity_type1) :: hc_v ! HC velocity structure



    ! Declare local variables.
    Real :: Random_Value
    Real :: Random_Angle
 
    ! jdemod - the code relating neutral temperature to the background plasma was taken from a DIVIMP recombined ion routine and not
    !          from a surface interaction routine - thus the option and the vel_mult_recomb_neut multiplier are not applicable here
    HC_Velocity = 1.38E4 * SQRT (HC_Temperature / Find_HC_Mass (HC_Species,H_Isotope_Composition))

    ! Check for neutral heating.
    !If ( Impurity_Neutral_Vel_Opt .eq. 0) Then
    !   HC_Velocity = 1.38E4 * SQRT (HC_Temperature / Find_HC_Mass (HC_Species,H_Isotope_Composition)) *  Vel_Mult_Recomb_Neut
    !Else If ( Impurity_Neutral_Vel_Opt .eq. 1 .or.  Impurity_Neutral_Vel_Opt .eq. 2) Then
    !   HC_Velocity = 1.38E4 * SQRT (gktibs (Current_Cell,Current_Ring) / Find_HC_Mass (HC_Species,H_Isotope_Composition)) *  &
    !                &Vel_Mult_Recomb_Neut
    !End If
 
    ! Adjust velocity for projection to 2-D poloidal plane.
    NRand = NRand + 1
    Call Surand2 (Seed, 1, Random_Value)
    Random_Angle =  Pi_Value * Random_Value -  Pi_Value / 2.0

    ! jdemod - limit random_angle to [PI/2-PI/180, -PI/2 + PI/180] to avoid V=0  
    if (abs(random_angle).gt.(pi_value/2.0 - Pi_value/180.0)) then 
       random_angle = sign(pi_value/2.0-pi_value/180.0,random_angle)
    endif

    ! jdemod
    !
    ! Remap the new velocity in 3D - include a toroidal component if one is required. 
    !
    if (hc_kinetics_opt.eq.1) then 

       if (debug_kinetics) &
          &  write(6,*) 'HC_Neut_Vessel_Sputter_Velocity'
       call reset_hc_velocity(hc_v,hc_velocity,current_angle,random_angle,hc_temperature)

    endif

    HC_Velocity = HC_Velocity * COS (Random_Angle)
 
    ! Calculate R,Z components of velocity, probability of ionization, etc.
    Current_Velocity_In_R = HC_Velocity * COS (Current_Angle) *  Neutral_Time_Step
    Current_Velocity_In_Z = HC_Velocity * SIN (Current_Angle) *  Neutral_Time_Step

  End Subroutine HC_Neut_Vessel_Sputter_Velocity
 
  Subroutine HC_Ion_Target_Sputter_Velocity (HC_Species,H_Isotope_Composition,Current_Cell,Current_Ring, &
       & Current_Angle,Target_Segment,HC_Temperature,HC_Velocity,Current_Velocity_In_R,Current_Velocity_In_Z,Seed,NRand,hc_v)
 
    ! Required modules.
    Use ComHC ! Access to global variables related to the hydrocarbon module.
    Use HC_Init_Lib_Data
    Use HC_Init_DIV_Data
    Use HC_Get
    use hc_kinetics

    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables.
    Integer, Intent (In) :: HC_Species !
    Integer, Dimension (Number_H_Species), Intent (In) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (In) :: H_Isotope_Composition
    Integer, Intent (In) :: Current_Cell ! IK
    Integer, Intent (In) :: Current_Ring ! IR
    Real, Intent (In) :: Current_Angle
    Integer, Intent (In) :: Target_Segment ! IN
    Real, Intent (In) :: HC_Temperature ! TEMI
    Real, Intent (Out) :: HC_Velocity ! VIN
    Real, Intent (Out) :: Current_Velocity_In_R
    Real, Intent (Out) :: Current_Velocity_In_Z
    Double Precision, Intent (In) :: Seed
    Integer, Intent (InOut) :: NRand

    type(hc_velocity_type1) :: hc_v ! HC velocity structure


 
    ! Declare local variables.
    Real :: Random_Value
    Real :: Random_Angle
 
    ! jdemod - the code relating neutral temperature to the background plasma was taken from a DIVIMP recombined ion routine and not
    !          from a surface interaction routine - thus the option and the vel_mult_recomb_neut multiplier are not applicable here
    HC_Velocity = 1.38E4 * SQRT (HC_Temperature / Find_HC_Mass (HC_Species,H_Isotope_Composition))

    ! Check for neutral heating.
    !If ( Impurity_Neutral_Vel_Opt .eq. 0) Then
    !   HC_Velocity = 1.38E4 * SQRT (HC_Temperature / Find_HC_Mass (HC_Species,H_Isotope_Composition)) *  Vel_Mult_Recomb_Neut
    !Else If ( Impurity_Neutral_Vel_Opt .eq. 1 .or.  Impurity_Neutral_Vel_Opt .eq. 2) Then
    !   HC_Velocity = 1.38E4 * SQRT (gktibs (Current_Cell,Current_Ring) / Find_HC_Mass (HC_Species,H_Isotope_Composition)) *  &
    !                &Vel_Mult_Recomb_Neut
    !End If
 
    ! Adjust velocity for projection to 2-D poloidal plane.
    NRand = NRand + 1
    Call Surand2 (Seed, 1, Random_Value)
    Random_Angle =  Pi_Value * Random_Value -  Pi_Value / 2.0

    ! jdemod - limit random_angle to [PI/2-PI/180, -PI/2 + PI/180] to avoid V=0  
    if (abs(random_angle).gt.(pi_value/2.0 - Pi_value/180.0)) then 
       random_angle = sign(pi_value/2.0-pi_value/180.0,random_angle)
    endif

    ! jdemod
    !
    ! Remap the new velocity in 3D - include a toroidal component if one is required. 
    !
    if (hc_kinetics_opt.eq.1) then 

       if (debug_kinetics) &
           & write(6,*) 'HC_Ion_Target_Sputter_Velocity'
       call reset_hc_velocity(hc_v,hc_velocity,current_angle,random_angle,hc_temperature)

    endif

 
    HC_Velocity = HC_Velocity * COS (Random_Angle)
 
    ! Calculate R,Z components of velocity, probability of ionization, etc.
    Current_Velocity_In_R = HC_Velocity * COS (Current_Angle) *  Neutral_Time_Step
    Current_Velocity_In_Z = HC_Velocity * SIN (Current_Angle) *  Neutral_Time_Step
 
    !Write (Output_Unit_Scratch,*) "HC_Ion_Target_Sputter_Velocity", Injection_Opt,Current_Cell,Current_Ring,Random_Angle,HC_Velocity

  End Subroutine HC_Ion_Target_Sputter_Velocity
 
  Subroutine HC_Neut_Vessel_Sputter_Location (NewR,NewZ,Current_R,Current_Z, &
       & Current_Cell,Current_Ring,Current_S,Current_Cross,Current_Theta,S_Start)
 
    Use ComHC
    Use HC_Init_DIV_Data
 
    ! Every good Fortran 90 program has...		
    Implicit None
 
    ! Define input/output variables
    Real, Intent (In) :: NewR
    Real, Intent (In) :: NewZ
    Real, Intent (Out) :: Current_R
    Real, Intent (Out) :: Current_Z
    Integer, Intent (Out) :: Current_Cell
    Integer, Intent (Out) :: Current_Ring
    Real, Intent (Out) :: Current_S
    Real, Intent (Out) :: Current_Cross
    Real, Intent (Out) :: Current_Theta
    Real, Intent (Out) :: S_Start
 
    ! Update position of particle.  Note, all that is needed is R and Z.
    Current_R = NewR
    Current_Z = NewZ
 
    Call gridpos (Current_Cell, Current_Ring, Current_R, Current_Z, .False.,  Grid_Error)
    ! Returns R, Z.
    !Write (Output_Unit_Scratch,*) "sputtering location:",Current_R,Current_Z, Grid_Error
    ! Reset non-applicable variables.
    Current_S = 0.0
    Current_Cross = 0.0
    Current_Theta = 0.0
    S_Start = 0.0
 
  End Subroutine HC_Neut_Vessel_Sputter_Location
 
  Subroutine Sputter_Species_Alman_Ruzic (HC_Species,Vessel_Segment,R_Position,Z_Position, &
       & HC_Temperature,HC_Velocity,HC_Velocity_In_S,HC_Angle_To_Normal,Surface_Temperature)
    ! Use Alman and Ruzic PSI2002 paper data for sputtering hydrocarbon species.
 
    Use HC_Utilities ! Access hydrocarbon utilities.
    Use HC_Stack ! Gain access to hydrocarbon storage stack and operations for hydrocarbon breakup into multiple species.
    Use ComHC
 
    ! Every good Fortran 90 program has...		
    Implicit None
 
    ! Define input/output variables
    Integer, Intent (InOut) :: HC_Species
    Integer, Intent (In) :: Vessel_Segment
    Real, Intent (In) :: R_Position
    Real, Intent (In) :: Z_Position
    Real, Intent (InOut) :: HC_Temperature
    Real, Intent (InOut) :: HC_Velocity
    Real, Intent (InOut) :: HC_Velocity_In_S
    Real, Intent (InOut) :: HC_Angle_To_Normal
    Real, Intent (In) :: Surface_Temperature
 
    ! Define local variables that will be assigned values.
    Integer, Dimension (Number_HC_Species,Highest_Carbon_Content) :: Alman_Ruzic_Sputtering_Table ! Currently 58 x 3 possibilities for reflect
 
 
    ! Push sputtered hydrocarbon(s) onto carbon storage stack.
    !Call HC_Push (HC_Species, R_Location, Z_Location, S_Location, Cross_Location, HC_Direction, Current_Theta, HC_Velocity, HC_Velocity_In_S, HC_Temperature)
 
  End Subroutine Sputter_Species_Alman_Ruzic
 
  Subroutine Sputter_Angle_Alman_Ruzic (HC_Species,Vessel_Segment,R_Position,Z_Position, &
       & HC_Temperature,HC_Velocity,HC_Velocity_In_S,HC_Angle_To_Normal,Surface_Temperature)
    ! Use Alman and Ruzic PSI2002 paper data for sputtering angle.
 
    Use HC_Utilities ! Access hydrocarbon utilities.
    Use HC_Stack ! Gain access to hydrocarbon storage stack and operations for hydrocarbon breakup into multiple species.
    Use ComHC
 
    ! Every good Fortran 90 program has...		
    Implicit None
 
    ! Define input/output variables
    Integer, Intent (InOut) :: HC_Species
    Integer, Intent (In) :: Vessel_Segment
    Real, Intent (In) :: R_Position
    Real, Intent (In) :: Z_Position
    Real, Intent (InOut) :: HC_Temperature
    Real, Intent (InOut) :: HC_Velocity
    Real, Intent (InOut) :: HC_Velocity_In_S
    Real, Intent (InOut) :: HC_Angle_To_Normal
    Real, Intent (In) :: Surface_Temperature
 
    ! Define local variables that will be assigned values.
    Integer, Dimension (Number_HC_Species,Highest_Carbon_Content) :: Alman_Ruzic_Reflection_Table ! Currently 58 x 3 possibilities for reflection table.
 
 
  End Subroutine Sputter_Angle_Alman_Ruzic
 
  Subroutine Sputter_EnergyVel_Alman_Ruzic (HC_Species,Vessel_Segment,R_Position, &
       & Z_Position,HC_Temperature,HC_Velocity,HC_Angle_To_Normal,Surface_Temperature)
    ! Use Alman and Ruzic PSI2002 paper data for sputtered energy and velocity.
 
    Use HC_Utilities ! Access hydrocarbon utilities.
    Use HC_Stack ! Gain access to hydrocarbon storage stack and operations for hydrocarbon breakup into multiple species.
    Use ComHC
 
    ! Every good Fortran 90 program has...		
    Implicit None
 
    ! Define input/output variables
    Integer, Intent (InOut) :: HC_Species
    Integer, Intent (In) :: Vessel_Segment
    Real, Intent (In) :: R_Position
    Real, Intent (In) :: Z_Position
    Real, Intent (InOut) :: HC_Temperature
    Real, Intent (InOut) :: HC_Velocity
    Real, Intent (InOut) :: HC_Angle_To_Normal
    Real, Intent (In) :: Surface_Temperature
 
    ! Define local variables that will be assigned values.
    Integer, Dimension (Number_HC_Species,Highest_Carbon_Content) :: Alman_Ruzic_Sputtering_Table ! Currently 58 x 3 possibilities for reflect
 
 
  End Subroutine Sputter_EnergyVel_Alman_Ruzic
 
  Subroutine HC_Target_Sputter_Species (HC_Species,Target_Segment,Current_Theta, &
       & R_Position,Z_Position,HC_Temperature,HC_Velocity,HC_Velocity_In_S,HC_Angle)
 
    ! Required modules.
    Use HC_Init_Lib_Data ! Gain access to state and transition tables.
    Use HC_Stack ! Gain access to hydrocarbon storage stack and operations.
    Use ComHC ! Access to global variables related to the hydrocarbon module.
    Use HC_Utilities ! Access hydrocarbon utilities.
    Use HC_Get
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables
    Integer, Intent (InOut) :: HC_Species ! Note, only one species passed in and out, additional species loaded in the HC stack.
    Integer, Intent (In) :: Target_Segment ! IN
    Real, Intent (In) :: Current_Theta ! THETA
    Real, Intent (In) :: R_Position ! R
    Real, Intent (In) :: Z_Position ! Z
    Real, Intent (InOut) :: HC_Temperature ! TEMI
    Real, Intent (InOut) :: HC_Velocity ! VEL
    Real, Intent (InOut) :: HC_Velocity_In_S ! VEL
    Real, Intent (InOut) :: HC_Angle
 
    ! Declare local variables.
    Integer :: i
 
    ! Proceed to sputter particle.
    ! First, find new hydrocarbon species (H28).
    If (hc_sputtering_species_model .eq. 0) Then
       ! Use preset sputtering table to decide which single hydrocarbon species is emitted.
       ! Note, this may change to multiple Carbon atoms if the primary channel of sputtering for C2 or higher is changed.
       HC_Species = HC_Sputtering_Table (HC_Species) % Sputtered_State_Number (1)
       Do i = 2, Highest_Carbon_Content, 1 ! Up to 3 for C3Hx series of hydrocarbons.
          If (HC_Sputtering_Table (HC_Species) % Sputtered_State_Number (i) .ne. 0) Then
             ! Multiple sputtered hydocarbon species.  Store additional data on HC stack.
             ! Note that current S and Cross are given values of 0.0 because particles are always sputtered from the wall as neutrals.
             Call HC_Push (HC_Species,R_Position,Z_Position,0.0,0.0,HC_Angle,Current_Theta,HC_Velocity,HC_Velocity_In_S,&
                          &HC_Temperature)
          End If
       End Do
 
    ElseIf (hc_sputtering_species_model .eq. 1) Then
       ! Use Alman and Ruzic sputtering data.  Note the possibility of multiple hydrocarbons released and returned in HC_Species_Out.
       !Call Sputter_Species_Alman_Ruzic (HC_Species,Target_Segment,R_Position,Z_Position,HC_Temperature,HC_Velocity,HC_Angle_To_Normal,Surface_Temperature)
    ElseIf (hc_sputtering_species_model .eq. 2) Then
       ! Sputtered species equals incoming species.  Do nothing.
    Else
       ! Unsupported sputtering model selected.
       Write (Output_Unit_HC_Alert,*) "Unsupported sputtering model selected:", hc_sputtering_species_model
       Write (Output_Unit_HC_Alert,*) "Program Stopping."
       Stop
    End If
 
  End Subroutine HC_Target_Sputter_Species
 
  Subroutine HC_Target_Reflect_Species (HC_Species,Target_Segment,HC_Temperature,HC_Velocity,HC_Angle)
 
    ! Required modules.
    Use HC_Init_Lib_Data ! Gain access to state and transition tables.
    Use ComHC ! Access to global variables related to the hydrocarbon module.
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables
    Integer, Intent (InOut) :: HC_Species ! Note, only one species passed in and out, additional species loaded in the HC stack.
    Integer, Intent (In) :: Target_Segment ! IN
    Real, Intent (In) :: HC_Temperature ! TEMI
    Real, Intent (In) :: HC_Velocity ! VEL
    Real, Intent (In) :: HC_Angle
 
    ! Find new hydrocarbon species (H28).
    If (hc_reflection_species_model .eq. 0) Then
       ! Use preset reflection table to decide which single hydrocarbon species is emitted.
       ! Note, this may change to multiple Carbon atoms if the primary channel of reflection for C2 or higher is changed.
       HC_Species = HC_Reflection_Table (HC_Species) % Reflected_State_Number
    ElseIf (hc_reflection_species_model .eq. 1) Then
       ! Use Alman and Ruzic reflection data.  Note the possibility of multiple hydrocarbons released and returned in HC_Species_Out.
       !Call Reflect_Species_Alman_Ruzic (HC_Species,Vessel_Segment,R_Position,Z_Position,HC_Temperature,HC_Velocity,HC_Angle_To_Normal,Surface_Temperature)
    Else
       ! Unsupported reflection model selected.
       Write (Output_Unit_HC_Alert,*) "Unsupported reflection model selected:", hc_reflection_species_model
       Write (Output_Unit_HC_Alert,*) "Program Stopping."
       Stop
    End If
 
  End Subroutine HC_Target_Reflect_Species
 
  Subroutine HC_Vessel_Reflect_Species (HC_Species,Vessel_Segment,HC_Temperature,HC_Velocity,HC_Velocity_In_S,HC_Angle)
 
    ! Required modules.
    Use HC_Init_Lib_Data ! Gain access to state and transition tables.
    Use ComHC ! Access to global variables related to the hydrocarbon module.
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables
    Integer, Intent (InOut) :: HC_Species ! Note, only one species passed in and out, additional species loaded in the HC stack.
    Integer, Intent (In) :: Vessel_Segment ! ID
    Real, Intent (In) :: HC_Temperature ! TEMI
    Real, Intent (In) :: HC_Velocity ! VEL
    Real, Intent (In) :: HC_Velocity_In_S ! VEL
    Real, Intent (In) :: HC_Angle
 
    ! Find new hydrocarbon species (H28).
    If (hc_reflection_species_model .eq. 0) Then
       ! Use preset reflection table to decide which single hydrocarbon species is emitted.
       ! Note, this may change to multiple Carbon atoms if the primary channel of reflection for C2 or higher is changed.
       HC_Species = HC_Reflection_Table (HC_Species) % Reflected_State_Number
    ElseIf (hc_reflection_species_model .eq. 1) Then
       ! Use Alman and Ruzic reflection data.  Note the possibility of multiple hydrocarbons released and returned in HC_Species_Out.
       !Call Reflect_Species_Alman_Ruzic (HC_Species,Vessel_Segment,R_Position,Z_Position,HC_Temperature,HC_Velocity,HC_Velocity_In_S,HC_Angle_To_Normal,Surface_Temperature)
    Else
       ! Unsupported reflection model selected.
       Write (Output_Unit_HC_Alert,*) "Unsupported reflection model selected:", hc_reflection_species_model
       Write (Output_Unit_HC_Alert,*) "Program Stopping."
       Stop
    End If
 
  End Subroutine HC_Vessel_Reflect_Species
 
  Subroutine HC_Neut_Vessel_Reflect_Angle (HC_Species,Incoming_Angle,Incoming_Energy, &
       & Vessel_Segment,Reflection_Angle,Segment_Normal_Angle,HC_Angle,Seed,NRand)
 
    ! Required modules.
    Use ComHC ! Access to global variables related to the hydrocarbon module.
    Use HC_Init_Div_Data ! Access to misc_data_table.
    Use HC_Get
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables.
    Integer, Intent (In) :: HC_Species !
    Real, Intent (In) :: Incoming_Angle
    Real, Intent (In) :: Incoming_Energy ! TEMI
    Integer, Intent (In) :: Vessel_Segment ! ID
    Real, Intent (In) :: Reflection_Angle ! TREF
    Real, Intent (In) :: Segment_Normal_Angle ! TNORM
    Real, Intent (Out) :: HC_Angle
    Double Precision, Intent (In) :: Seed
    Integer, Intent (InOut) :: NRand
 
    ! Declare local variables.
    Real :: HC_Angle_To_Normal
    Real :: Random_Value_1
    Real :: Random_Value_2
    Real :: Beta, Psi

    !write(output_unit_scratch,'(a,i8,10(1x,g20.12))') 'HC_NEUT_VESSEL_REFLECT_ANGLE1:', hc_reflection_angle_model,&
    !   &HC_Species,Incoming_Angle*raddeg,Incoming_Energy, &
    !   & Reflection_Angle*raddeg,Segment_Normal_Angle*raddeg,HC_Angle*raddeg

 
    ! Find angle to normal given vessel segment and impact angle.
    HC_Angle_To_Normal = ABS (-Segment_Normal_Angle - Incoming_Angle)
 
    ! Find direction of particle after reflection (H31).
    If (hc_reflection_angle_model .eq. 0) Then
       ! This option should never happen.
       Write (Output_Unit_HC_Alert,*) "Error: Should not have entered function neutral vessel hc_reflect for no reflect condition: &
&HC_Neut_Vessel_Reflect_Angle."
       Write (Output_Unit_HC_Alert,*) "Program stopping."
       Stop
    ElseIf (hc_reflection_angle_model .ge. 1 .and. hc_reflection_angle_model .le. 4) Then
       ! Angle is supported by main REFANG DIVIMP code.
       ! =1, Specular reflection from wall, bounds of +/-90 degrees from normal to vessel wall.
       ! =2, Isotropic reflection with ACOS(x) distribution from vessel normal, bounds of +/-90 degrees from normal to vessel wall.
       ! =3, Isotropic reflection with ASIN(x) distribution from vessel normal, bounds of +/-90 degrees from normal to vessel wall.
       ! =4, Isotropic reflection with ASIN(SQRT(x)) distribution from vessel normal, bounds of +/-90 degrees from normal to vessel wall.
       HC_Angle = Reflection_Angle
 
    ElseIf (hc_reflection_angle_model .eq. 5) Then
 
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value_1)
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value_2)
 
       Beta = ASIN (SQRT(Random_Value_1))
       Psi = 2.0 *  Pi_Value * Random_Value_2
       HC_Angle = Reflection_Angle + ATAN (TAN (Beta) * COS (Psi))
 
    ElseIf (hc_reflection_angle_model .eq. 10) Then
       ! Normal reflection from vessel wall
       HC_Angle = Segment_Normal_Angle
 
    ElseIf (hc_reflection_angle_model .eq. 11) Then
       ! Alman and Ruzic angle data used.  Ignore Reflection_Angle (TREF) returned by INTCALC.
       !Call Reflect_Angle_Alman_Ruzic (HC_Species,Vessel_Segment,R_Position,Z_Position,HC_Temperature,HC_Velocity,HC_Angle_To_Normal,Surface_Temperature)
 
    Else		
       ! Unsupported reflection model selected.
       Write (Output_Unit_HC_Alert,*) "Unsupported reflection model selected: ", hc_reflection_angle_model
       Write (Output_Unit_HC_Alert,*) "Program Stopping."	
    End If


    !write(output_unit_scratch,'(a,i8,10(1x,g20.12))') 'HC_NEUT_VESSEL_REFLECT_ANGLE1:', hc_reflection_angle_model,&
    !   &HC_Species,Incoming_Angle*raddeg,Incoming_Energy, &
    !   & Reflection_Angle*raddeg,Segment_Normal_Angle*raddeg,HC_Angle*raddeg


 
  End Subroutine HC_Neut_Vessel_Reflect_Angle
 
  Subroutine HC_Ion_Target_Reflect_Angle (HC_Species,Incoming_Angle,Incoming_Energy,Target_Segment,Reflection_Angle,Seed,NRand)
 
    ! Required modules.
    Use ComHC ! Access to global variables related to the hydrocarbon module.
    Use HC_Init_DIV_Data
    Use HC_Get
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables.
    Integer, Intent (In) :: HC_Species !
    Real, Intent (In) :: Incoming_Angle
    Real, Intent (In) :: Incoming_Energy ! TEMI
    Integer, Intent (In) :: Target_Segment ! ID
    Real, Intent (Out) :: Reflection_Angle ! TREF
    Double Precision, Intent (In) :: Seed
    Integer, Intent (InOut) :: NRand
 
    ! hc_reflection_angle_model ! Decide at what angle a reflected HC ejects the vessel wall at: -1=NRFOPT, 0=off, 1=specular, 2=isotropic, 3=normal, 4=Alman and Ruzic angle data (H33).	
 
    ! Local variables.
    Real :: Random_Value_1
    Real :: Random_Value_2
    Real :: Beta, Psi

    !write(output_unit_scratch,'(a,i8,10(1x,g20.12))') 'HC_ION_TARGET_REFLECT_ANGLE1:', hc_reflection_angle_model,target_segment,&
    !   &HC_Species,Incoming_Angle,Incoming_Energy, &
    !   & Reflection_Angle

 
    ! Find direction of particle after reflection (H31).
    If (hc_reflection_angle_model .eq. 0) Then
       ! This option should never happen.
       Write (Output_Unit_HC_Alert,*) "Error: Should not have entered function ion target hc_reflect for no reflect condition: &
&HC_Ion_Target_Reflect_Angle."
       Write (Output_Unit_HC_Alert,*) "Program stopping/"
       Stop
    ElseIf (hc_reflection_angle_model .eq. 1) Then
       ! =1, Specular reflection from wall, bounds of +/-90 degrees from normal to vessel wall.
       ! jdemod - Ion specular reflection doesn't make any sense - what is the incoming angle in this case?
       reflection_angle = 2.0*gthetas(target_segment) + sign ((PI-abs(incoming_angle)),-incoming_angle)
       Write (Output_Unit_HC_Alert,*) "Specular reflection not yet supported.  Program stopping."
       Stop
    ElseIf (hc_reflection_angle_model .eq. 2) Then
       ! =2, Isotropic reflection from vessel wall, bounds of +/-90 degrees from normal to vessel wall.
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value_1)
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value_2)
       Reflection_Angle = gthetas (Target_Segment) + SIGN (ACOS (Random_Value_1), Random_Value_2 - 0.5)
 
    ElseIf (hc_reflection_angle_model .eq. 3) Then
       ! =3, Isotropic sputtering COS distribution from vessel normal, bounds of +/-90 degrees from normal to vessel wall.
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value_1)
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value_2)
       Reflection_Angle = gthetas (Target_Segment) + SIGN (ASIN (Random_Value_1), Random_Value_2 - 0.5)
 
    ElseIf (hc_reflection_angle_model .eq. 4) Then
       ! =3, Isotropic sputtering SQRT SIN distribution from vessel normal, bounds of +/-90 degrees from normal to vessel wall.
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value_1)
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value_2)
       Reflection_Angle = gthetas (Target_Segment) + SIGN (ASIN (SQRT (Random_Value_1)), Random_Value_2 - 0.5)
 
    ElseIf (hc_reflection_angle_model .eq. 5) Then
       ! =5, Isotropic sputtering projected SQRT SIN distribution from vessel normal, bounds of +/-90 degrees from normal to vessel wall.
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value_1)
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value_2)
 
       Beta = ASIN (SQRT(Random_Value_1))
       Psi = 2.0 *  Pi_Value * Random_Value_2
       Reflection_Angle = gthetas (Target_Segment) + ATAN (TAN (Beta) * COS (Psi))
 
    ElseIf (hc_reflection_angle_model .eq. 10) Then
       ! Normal reflection from vessel wall
       Reflection_Angle = gthetas (Target_Segment)
 
    ElseIf (hc_reflection_angle_model .eq. 11) Then
       ! Alman and Ruzic angle data used.  Ignore Reflection_Angle (TREF) returned by INTCALC.
       !Call Reflect_Angle_Alman_Ruzic (HC_Species,Vessel_Segment,R_Position,Z_Position,HC_Temperature,HC_Velocity,HC_Angle_To_Normal,Surface_Temperature)
 
    Else		
       ! Unsupported reflection model selected.
       Write (Output_Unit_HC_Alert,*) "Unsupported reflection model selected: ", hc_reflection_angle_model
       Write (Output_Unit_HC_Alert,*) "Program Stopping."	
 
    End If

    

    !write(output_unit_scratch,'(a,i8,10(1x,g20.12))') 'HC_ION_TARGET_REFLECT_ANGLE2:', hc_reflection_angle_model,target_segment,&
    !   &HC_Species,Incoming_Angle,Incoming_Energy, &
    !   & Reflection_Angle


 
  End Subroutine HC_Ion_Target_Reflect_Angle
 
  Subroutine HC_Neut_Vessel_Reflect_Energy (HC_Species,H_Isotope_Composition,Vessel_Segment,Reflected_Energy)
 
    ! Required modules.
    Use ComHC ! Access to global variables related to the hydrocarbon module.
    Use HC_Init_Lib_Data
    Use HC_Utilities
    Use HC_Init_DIV_Data
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables.
    Integer, Intent (In) :: HC_Species !
    Integer, Dimension (Number_H_Species), Intent (In) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (In) :: H_Isotope_Composition
    Integer, Intent (In) :: Vessel_Segment ! IN
    Real, Intent (Out) :: Reflected_Energy ! TEMN
 
    ! Declare local variables.
    Real :: HC_Velocity ! VIN
 
    ! Find energy and velocity of particle after reflection (H29, H30).
    If (hc_reflection_energy_model .eq. 0) Then
       ! Use preset energy and calculate velocity.
       Reflected_Energy = hc_refl_energy_neutral_preset
 
    ElseIf (hc_reflection_energy_model .eq. 1) Then
       ! Do not change energy at all - reflected particle energy is same as impact.
       ! Recalculate velocity for potential new mass.
 
    ElseIf (hc_reflection_energy_model .eq. 2) Then
       ! Assign new temperature as well.			
       ! jdemod - I don't think the following expression can be correct since the 8kT/(PI m) term does not 
       !          convert the hc_mass to amu which is required given that k is defined as 1.38e-23)
       ! Assign new temperature as well
       ! jdemod - added AMU and removed 1.38e4 factor which is sqrt(2e/amu) - the expression seems almost correct - needed an AMU in the HC_velocity term
       !HC_Velocity = SQRT (8.0 * 1.38066E-23 * Find_Wall_Temperature (Vessel_Segment) / ( Pi_Value * Find_HC_Mass (HC_Species,&
       !             &H_Isotope_Composition)))
       !Reflected_Energy = Find_HC_Mass (HC_Species,H_Isotope_Composition) * (HC_Velocity / 1.38E4) * (HC_Velocity / 1.38E4)

       ! IPP/09 Krieger - shortened line (SUNWorkshop chokes over len>132)

       HC_Velocity = SQRT(8.0 * 1.38066E-23 *  Find_Wall_Temperature(Vessel_Segment) / (Pi_Value * AMU * Find_HC_Mass(HC_Species,&
                          &H_Isotope_Composition)))
       Reflected_Energy = 0.5* AMU / ech * Find_HC_Mass (HC_Species,H_Isotope_Composition) * HC_Velocity**2


 
    ElseIf (hc_reflection_energy_model .eq. 3) Then
       ! Particle energy found with Alman & Ruzic data.
       !Call Reflect_EnergyVel_Alman_Ruzic (HC_Species,Target_Segment,R_Position,Z_Position,HC_Temperature,HC_Velocity,HC_Angle_To_Normal,Surface_Temperature)
    Else
       ! Unsupported reflection model selected.
       Write (Output_Unit_HC_Alert,*) "Unsupported reflection energy/vel model selected: ", hc_reflection_energy_model
       Write (Output_Unit_HC_Alert,*) "Program Stopping."
       Stop
    End If
 
  End Subroutine HC_Neut_Vessel_Reflect_Energy
 
  Subroutine HC_Ion_Target_Reflect_Energy (HC_Species,H_Isotope_Composition, &
       & Current_Cell,Current_Ring,Current_Angle,Target_Index,Reflected_Energy)
 
    ! Required modules.
    Use ComHC ! Access to global variables related to the hydrocarbon module.
    Use HC_Init_Lib_Data
    Use HC_Utilities
    Use HC_Init_DIV_Data
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables.
    Integer, Intent (In) :: HC_Species
    Integer, Dimension (Number_H_Species), Intent (In) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (In) :: H_Isotope_Composition
    Integer, Intent (In) :: Current_Cell
    Integer, Intent (In) :: Current_Ring
    Real, Intent (In) :: Current_Angle
    Integer, Intent (In) :: Target_Index ! IN
    Real, Intent (Out) :: Reflected_Energy ! TEMN
 
    ! Declare local variables.
    Real :: HC_Velocity ! VIN
 
    ! Find energy and velocity of particle after reflection (H29, H30).
    If (hc_reflection_energy_model .eq. 0) Then
       ! Use preset energy and calculate velocity.
       Reflected_Energy = hc_refl_energy_ion_preset
 
    ElseIf (hc_reflection_energy_model .eq. 1) Then
       ! Do not change energy at all - reflected particle energy is same as impact.
       ! Recalculate velocity for potential new mass.
 
    ElseIf (hc_reflection_energy_model .eq. 2) Then
       ! Assign new temperature as well
       ! jdemod - added AMU and removed 1.38e4 factor which is sqrt(2e/amu) - the expression seems almost correct - needed an AMU in the HC_velocity term
       !HC_Velocity = SQRT (8.0 * 1.38066E-23 * Find_Target_Temperature (Target_Index) / ( Pi_Value  * Find_HC_Mass (HC_Species,&
       !                   &H_Isotope_Composition)))
       !Reflected_Energy = Find_HC_Mass (HC_Species,H_Isotope_Composition) * (HC_Velocity / 1.38E4) * (HC_Velocity / 1.38E4)

       ! IPP/09 Krieger - shortened line (SUNWorkshop chokes over len>132)

       HC_Velocity = SQRT(8.0 * 1.38066E-23 * Find_Target_Temperature(Target_Index) / (Pi_Value * AMU * Find_HC_Mass(HC_Species,&
                          &H_Isotope_Composition)))
       Reflected_Energy = 0.5* AMU / ech * Find_HC_Mass (HC_Species,H_Isotope_Composition) * HC_Velocity**2

 
    ElseIf (hc_reflection_energy_model .eq. 3) Then
       ! Particle energy found with Alman & Ruzic data.
       !Call Reflect_EnergyVel_Alman_Ruzic (HC_Species,Target_Index,R_Position,Z_Position,HC_Temperature,HC_Velocity,HC_Angle_To_Normal,Surface_Temperature)
    Else
       ! Unsupported reflection model selected.
       Write (Output_Unit_HC_Alert,*) "Unsupported reflection energy/vel model selected: ", hc_reflection_energy_model
       Write (Output_Unit_HC_Alert,*) "Program Stopping."
       Stop
    End If
 
  End Subroutine HC_Ion_Target_Reflect_Energy
 
  Subroutine HC_Neut_Vessel_Reflect_Velocity (HC_Species,H_Isotope_Composition,Current_Cell,Current_Ring, &
       & Current_Angle,Vessel_Segment,HC_Temperature,HC_Velocity,Current_Velocity_In_R,Current_Velocity_In_Z,Seed,NRand,hc_v)
 
    ! Required modules.
    Use ComHC ! Access to global variables related to the hydrocarbon module.
    Use HC_Init_Lib_Data
    Use HC_Init_DIV_Data
    Use HC_Get
    use hc_kinetics

    ! Every good Fortran program has...
    Implicit None
 
    ! Define input/output variables.
    Integer, Intent (In) :: HC_Species !
    Integer, Dimension (Number_H_Species), Intent (In) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (In) :: H_Isotope_Composition
    Integer, Intent (In) :: Current_Cell ! IK
    Integer, Intent (In) :: Current_Ring ! IR
    Real, Intent (In) :: Current_Angle
    Integer, Intent (In) :: Vessel_Segment ! IN
    Real, Intent (In) :: HC_Temperature ! TEMI
    Real, Intent (Out) :: HC_Velocity ! VIN
    Real, Intent (Out) :: Current_Velocity_In_R
    Real, Intent (Out) :: Current_Velocity_In_Z
    Double Precision, Intent (In) :: Seed
    Integer, Intent (InOut) :: NRand

    type(hc_velocity_type1) :: hc_v ! HC velocity structure

 
    ! Declare local variables.
    Real :: Random_Value
    Real :: Random_Angle
 
    ! jdemod - the code relating neutral temperature to the background plasma was taken from a DIVIMP recombined ion routine and not
    !          from a surface interaction routine - thus the option and the vel_mult_recomb_neut multiplier are not applicable here
    HC_Velocity = 1.38E4 * SQRT (HC_Temperature / Find_HC_Mass (HC_Species,H_Isotope_Composition)) 


    ! Check for neutral heating.
    !If ( Impurity_Neutral_Vel_Opt .eq. 0) Then
    !   HC_Velocity = 1.38E4 * SQRT (HC_Temperature / Find_HC_Mass (HC_Species,H_Isotope_Composition)) *  Vel_Mult_Recomb_Neut
    !Else If ( Impurity_Neutral_Vel_Opt .eq. 1 .or.  Impurity_Neutral_Vel_Opt .eq. 2) Then
    !   HC_Velocity = 1.38E4 * SQRT (gktibs (Current_Cell,Current_Ring) / Find_HC_Mass (HC_Species,H_Isotope_Composition)) *  &
    !                &Vel_Mult_Recomb_Neut
    !End If
 
    !write(0,'(a,3g18.6)') 'HC_REF_VEL:',hc_velocity, current_velocity_in_r,current_velocity_in_z

    if (hc_reflection_angle_model.eq.10) then 
       random_angle = 0.0
    else

       ! Adjust velocity for projection to 2-D poloidal plane.
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value)

       ! jdemod - random angle represents a projection angle into the toroidal direction
       !          it is assumed to be evenly distributed in the range -PI/2 to PI/2 
       !          the remaining projected HC_velocity in the R,Z plane is thus HC_Velocity * cos(random_angle)
       !          However, if random_angle is actually equal to +/- PI/2 the projected velocity -> 0 which
       !          results in a zero or near zero particle velocity at a material surface which is an issue. 
       !
       !          Fix this by:                  1) this maximum value for this is limited to 1 degree less than 
       !                                           PI/2 ... = [PI/2 - PI/180.0, -PI/2 + PI/180]
       !                                           
       Random_Angle =  Pi_Value * Random_Value -  Pi_Value / 2.0

       ! jdemod - limit random_angle to [PI/2-PI/180, -PI/2 + PI/180] to avoid V=0  
       if (abs(random_angle).gt.(pi_value/2.0 - Pi_value/180.0)) then 
          random_angle = sign(pi_value/2.0-pi_value/180.0,random_angle)
       endif

    endif

    !Write (Output_Unit_Scratch,'(a,2i8,4(1x,g20.12))') "neut_vessel_reflect_Velocity",&
    !           &Current_Cell,Current_Ring,Current_Angle*raddeg,HC_Velocity,Random_Angle,COS(Random_Angle)

    ! jdemod
    !
    ! Remap the new velocity in 3D - include a toroidal component if one is required. 
    !
    if (hc_kinetics_opt.eq.1) then 

       if (debug_kinetics) & 
            &  write(6,*) 'HC_Neut_Vessel_Reflect_Velocity'
       call reset_hc_velocity(hc_v,hc_velocity,current_angle,random_angle,hc_temperature)

    endif


    HC_Velocity = HC_Velocity * COS (Random_Angle)
 
    ! Calculate R,Z components of velocity, probability of ionization, etc.
    Current_Velocity_In_R = HC_Velocity * COS (Current_Angle) *  Neutral_Time_Step
    Current_Velocity_In_Z = HC_Velocity * SIN (Current_Angle) *  Neutral_Time_Step
 
    !if (current_velocity_in_z.eq.0.0) then 
    !   write(0,'(a,2i8,10(1x,g20.12))') 'ERROR in neutral vessel reflect velocity:',Current_Cell,Current_Ring,Current_Angle*raddeg,&
    ! & current_velocity_in_r,current_velocity_in_z,HC_Velocity,Random_Angle,COS(Random_Angle)
    !endif
    
    !write(0,'(a,10g18.6)') 'HC_REF_VEL:',hc_velocity, current_velocity_in_r,current_velocity_in_z,random_angle

    !if (hc_velocity.lt.100) stop 'hc_vel'

  End Subroutine HC_Neut_Vessel_Reflect_Velocity




  Subroutine HC_Surface_interaction_Velocity (HC_Species,H_Isotope_Composition,Current_Cell,Current_Ring, &
       & Current_Angle,Vessel_Segment,HC_Temperature,HC_Velocity,Current_Velocity_In_R,Current_Velocity_In_Z,Seed,NRand,hc_v)
 
    !
    ! jdemod - this code is the same for all of the various velocity routines - so I am looking at whether it can 
    !          be factored out into one routine. There may be some differences in the calling signatures and code for 
    !          ions vs. neutrals
    !

    ! Required modules.
    Use ComHC ! Access to global variables related to the hydrocarbon module.
    Use HC_Init_Lib_Data
    Use HC_Init_DIV_Data
    Use HC_Get
    use hc_kinetics
 
    ! Every good Fortran program has...
    Implicit None
 
    ! Define input/output variables.
    Integer, Intent (In) :: HC_Species !
    Integer, Dimension (Number_H_Species), Intent (In) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (In) :: H_Isotope_Composition
    Integer, Intent (In) :: Current_Cell ! IK
    Integer, Intent (In) :: Current_Ring ! IR
    Real, Intent (In) :: Current_Angle
    Integer, Intent (In) :: Vessel_Segment ! IN
    Real, Intent (In) :: HC_Temperature ! TEMI
    Real, Intent (Out) :: HC_Velocity ! VIN
    Real, Intent (Out) :: Current_Velocity_In_R
    Real, Intent (Out) :: Current_Velocity_In_Z
    Double Precision, Intent (In) :: Seed
    Integer, Intent (InOut) :: NRand

    type(hc_velocity_type1) :: hc_v ! HC velocity structure

 
    ! Declare local variables.
    Real :: Random_Value
    Real :: Random_Angle
 
    ! jdemod - the code relating neutral temperature to the background plasma was taken from a DIVIMP recombined ion routine and not
    !          from a surface interaction routine - thus the option and the vel_mult_recomb_neut multiplier are not applicable here
    HC_Velocity = 1.38E4 * SQRT (HC_Temperature / Find_HC_Mass (HC_Species,H_Isotope_Composition))

    ! Check for neutral heating.
    !If ( Impurity_Neutral_Vel_Opt .eq. 0) Then
    !   HC_Velocity = 1.38E4 * SQRT (HC_Temperature / Find_HC_Mass (HC_Species,H_Isotope_Composition)) *  Vel_Mult_Recomb_Neut
    !Else If ( Impurity_Neutral_Vel_Opt .eq. 1 .or.  Impurity_Neutral_Vel_Opt .eq. 2) Then
    !   HC_Velocity = 1.38E4 * SQRT (gktibs (Current_Cell,Current_Ring) / Find_HC_Mass (HC_Species,H_Isotope_Composition)) *  &
    !                &Vel_Mult_Recomb_Neut
    !End If
 
    ! Adjust velocity for projection to 2-D poloidal plane.
    NRand = NRand + 1
    Call Surand2 (Seed, 1, Random_Value)
    Random_Angle =  Pi_Value * Random_Value -  Pi_Value / 2.0

    ! jdemod - limit random_angle to [PI/2-PI/180, -PI/2 + PI/180] to avoid V=0  
    if (abs(random_angle).gt.(pi_value/2.0 - Pi_value/180.0)) then 
       random_angle = sign(pi_value/2.0-pi_value/180.0,random_angle)
    endif


    !Write (Output_Unit_Scratch,*) "neut_vessel_reflect_Velocity",Current_Cell,Current_Ring,Current_Angle,HC_Velocity,Random_Angle,COS(Random_Angle)

    ! jdemod
    !
    ! Remap the new velocity in 3D - include a toroidal component if one is required. 
    !
    if (hc_kinetics_opt.eq.1) then 

       if (debug_kinetics) & 
            &  write(6,*) 'HC_Surface_interaction_Velocity'
       call reset_hc_velocity(hc_v,hc_velocity,current_angle,random_angle,hc_temperature)

    endif


    HC_Velocity = HC_Velocity * COS (Random_Angle)
 
    ! Calculate R,Z components of velocity, probability of ionization, etc.
    Current_Velocity_In_R = HC_Velocity * COS (Current_Angle) *  Neutral_Time_Step
    Current_Velocity_In_Z = HC_Velocity * SIN (Current_Angle) *  Neutral_Time_Step

  End Subroutine HC_Surface_Interaction_Velocity


 
  Subroutine HC_Ion_Target_Reflect_Velocity (HC_Species,H_Isotope_Composition,Current_Cell,Current_Ring, &
       & Current_Angle,Target_Segment,HC_Temperature,HC_Velocity,Current_Velocity_In_R,Current_Velocity_In_Z,Seed,NRand,hc_v)
 
    ! Required modules.
    Use ComHC ! Access to global variables related to the hydrocarbon module.
    Use HC_Init_Lib_Data
    Use HC_Init_DIV_Data
    Use HC_Get
    use hc_kinetics
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables.
    Integer, Intent (In) :: HC_Species !
    Integer, Dimension (Number_H_Species), Intent (In) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (In) :: H_Isotope_Composition
    Integer, Intent (In) :: Current_Cell ! IK
    Integer, Intent (In) :: Current_Ring ! IR
    Real, Intent (In) :: Current_Angle
    Integer, Intent (In) :: Target_Segment ! IN
    Real, Intent (In) :: HC_Temperature ! TEMI
    Real, Intent (Out) :: HC_Velocity ! VIN
    Real, Intent (Out) :: Current_Velocity_In_R
    Real, Intent (Out) :: Current_Velocity_In_Z
    Double Precision, Intent (In) :: Seed
    Integer, Intent (InOut) :: NRand

    type(hc_velocity_type1) :: hc_v ! HC velocity structure
 
    ! Declare local variables.
    Real :: Random_Value
    Real :: Random_Angle
 




    ! jdemod - the code relating neutral temperature to the background plasma was taken from a DIVIMP recombined ion routine and not
    !          from a surface interaction routine - thus the option and the vel_mult_recomb_neut multiplier are not applicable here
    HC_Velocity = 1.38E4 * SQRT (HC_Temperature / Find_HC_Mass (HC_Species,H_Isotope_Composition))
    ! Check for neutral heating.
    !If ( Impurity_Neutral_Vel_Opt .eq. 0) Then
    !   HC_Velocity = 1.38E4 * SQRT (HC_Temperature / Find_HC_Mass (HC_Species,H_Isotope_Composition)) *  Vel_Mult_Recomb_Neut
    !Else If ( Impurity_Neutral_Vel_Opt .eq. 1 .or.  Impurity_Neutral_Vel_Opt .eq. 2) Then
    !   HC_Velocity = 1.38E4 * SQRT (gktibs (Current_Cell,Current_Ring) / Find_HC_Mass (HC_Species,H_Isotope_Composition)) *  &
    !                 &Vel_Mult_Recomb_Neut
    !End If
 

    ! normal launch does not need a 3D-> 2D correction
    if (hc_reflection_angle_model.eq.10) then 
       random_angle = 0.0
    else
       ! Adjust velocity for projection to 2-D poloidal plane.
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value)
       Random_Angle =  Pi_Value * Random_Value -  Pi_Value / 2.0

      ! jdemod - limit random_angle to [PI/2-PI/180, -PI/2 + PI/180] to avoid V=0  
      if (abs(random_angle).gt.(pi_value/2.0 - Pi_value/180.0)) then 
         random_angle = sign(pi_value/2.0-pi_value/180.0,random_angle)
      endif


    endif

    !Write (Output_Unit_Scratch,'(a,2i8,4(1x,g20.12))') "ion_target_reflect_Velocity",&
    !    &Current_Cell,Current_Ring,Current_Angle*raddeg,HC_Velocity,Random_Angle,COS(Random_Angle)
    ! jdemod
    !
    ! Remap the new velocity in 3D - include a toroidal component if one is required. 
    !
    if (hc_kinetics_opt.eq.1) then 

       if (debug_kinetics) & 
            &  write(6,*) 'HC_Ion_Target_Reflect_Velocity'
       call reset_hc_velocity(hc_v,hc_velocity,current_angle,random_angle,hc_temperature)

    endif

 
    HC_Velocity = HC_Velocity * COS (Random_Angle)
 
    ! Calculate R,Z components of velocity, probability of ionization, etc.
    Current_Velocity_In_R = HC_Velocity * COS (Current_Angle) *  Neutral_Time_Step
    Current_Velocity_In_Z = HC_Velocity * SIN (Current_Angle) *  Neutral_Time_Step

    !if (current_velocity_in_z.eq.0.0) then 
    !   write(0,'(a,2i8,10(1x,g20.12))') 'ERROR in Ion target reflect velocity:',Current_Cell,Current_Ring,Current_Angle*raddeg,&
    ! & current_velocity_in_r,current_velocity_in_z,HC_Velocity,Random_Angle,COS(Random_Angle)
    !endif

  End Subroutine HC_Ion_Target_Reflect_Velocity
 
  Subroutine HC_Neut_Vessel_Reflect_Location (NewR,NewZ,Current_R,Current_Z, &
       & Current_Cell,Current_Ring,Current_S,Current_Cross,Current_Theta,S_Start)
 
    Use ComHC
    Use HC_Init_DIV_Data
 
    ! Every good Fortran 90 program has...		
    Implicit None
 
    ! Define input/output variables
    Real, Intent (In) :: NewR
    Real, Intent (In) :: NewZ
    Real, Intent (Out) :: Current_R
    Real, Intent (Out) :: Current_Z
    Integer, Intent (Out) :: Current_Cell
    Integer, Intent (Out) :: Current_Ring
    Real, Intent (Out) :: Current_S
    Real, Intent (Out) :: Current_Cross
    Real, Intent (Out) :: Current_Theta
    Real, Intent (Out) :: S_Start
 
    ! Update position of particle.  Note, all that is needed is R and Z.
    Current_R = NewR
    Current_Z = NewZ
 
    Call gridpos (Current_Cell,Current_Ring,Current_R,Current_Z,.False., Grid_Error)
    ! Returns R, Z.
 
    ! Reset non-applicable variables.
    Current_S = 0.0
    Current_Cross = 0.0
    Current_Theta = 0.0
    S_Start = 0.0
 
  End Subroutine HC_Neut_Vessel_Reflect_Location
 
  Subroutine HC_Ion_Target_Location (Current_S,Current_Cross,Current_Cell,Current_Ring,Current_R,Current_Z,Target_Index)
 
    ! Required modules.
    Use ComHC ! Access to global variables related to the hydrocarbon module.
    Use HC_Init_Lib_Data
    Use HC_Init_DIV_Data
    Use HC_Get
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables.
    Real, Intent (InOut) :: Current_S
    Real, Intent (InOut) :: Current_Cross
    Integer, Intent (InOut) :: Current_Cell
    Integer, Intent (In) :: Current_Ring
    Real, Intent (Out) :: Current_R
    Real, Intent (Out) :: Current_Z
    Integer, Intent (Out) :: Target_Index
 
    ! Declare local variables.
    Integer, External :: Verify_ID
 
    If (Current_S .le. 0.0) Then
       Current_S = 0.0
       Current_Cell = 1
       ! Verify_ID takes IK,IR, and target to look at.
       Target_Index = Verify_ID (Current_Cell,Current_Ring,2)
    Else
       Current_S  = gksmaxs (Current_Ring)
       Current_Cell = gnks (Current_Ring)
       Target_Index = Verify_ID (Current_Cell,Current_Ring,1)
    End If
 
    ! Postion on target/initial position options.
    If ( Neutral_Init_Pos_Opt .eq. 0) Then
       Current_R = grp (Target_Index)
       Current_Z = gzp (Target_Index)
    ElseIf ( Neutral_Init_Pos_Opt .eq. 1) Then
       Call Position_On_Target (Current_R,Current_Z,Current_Cross,Target_Index)
    End If
    !Write (Output_Unit_Scratch,*) "hc_ion_target_location",Current_S,Current_Cross,Current_cell,Current_Ring,Target_Index,Current_R,Current_Z
    Current_Cross = 0.0
 
  End Subroutine HC_Ion_Target_Location
 
End Module HC_Vessel_Interact
