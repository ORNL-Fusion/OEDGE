! -*-Mode:f90-*-
! Check for PROMPT REDEPOSITION of the ORIGINAL ION
 
Module HC_Prompt
 
  Implicit None
 
Contains
 
  Subroutine HC_Prompt_Deposition (NProd,Current_Cell,Current_Ring,Current_R,Current_Z,Current_S,Current_Cross,Current_Theta, &
       & Cur_HC_Spec,H_Isotope_Composition,Current_Velocity_In_S,Current_Velocity,Current_Velocity_In_R,Current_Velocity_In_Z, &
       & HC_Temperature,Last_HC_Species,Last_H_Isotope_Composition,Launch_Reg,Reflection_Counter,NeutType,Max_Velocity_Randoms, &
       & Tries_Without_Change,Neutral_Time_Step_Count,Ion_Time_Step_Count,Eq_Total_Ion_Time_Steps,Seed,NRand,IFate,hc_v,Debug)
 
    Use HC_Init_DIV_Data ! Gain access to data structures.
    Use HC_Init_DIV_Diag
    Use HC_Get ! gknbs, glambda1 functions.
    Use HC_Put
    Use HC_Init_Lib_Data ! Get HC mass and charge.
    Use HC_Utilities ! Get SPUTY.
    Use HC_Vessel_Interact ! Reflection and sputtering routines.
    Use HC_Stick ! Refl. and Sput. coefficients.
 
    use hc_velocity_type

    ! Every good Fortran routine has...
    Implicit None
 
    ! Declare call-line variables.
    Integer, Intent (InOut) :: NProd ! NPROD
    Integer, Intent (InOut) :: Current_Cell ! IK
    Integer, Intent (InOut) :: Current_Ring ! IR
    Real, Intent (InOut) :: Current_R ! R
    Real, Intent (InOut) :: Current_Z ! Z
    Real, Intent (InOut) :: Current_S ! S
    Real, Intent (InOut) :: Current_Cross ! CROSS
    Real, Intent (InOut) :: Current_Theta
    Integer, Intent (InOut) :: Cur_HC_Spec
    Integer, Dimension (Number_H_Species), Intent (InOut) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (InOut) :: H_Isotope_Composition
    Real, Intent (InOut) :: Current_Velocity_In_S
    Real, Intent (Out) :: Current_Velocity ! VEL
    Real, Intent (Out) :: Current_Velocity_In_R ! XVELF
    Real, Intent (Out) :: Current_Velocity_In_Z ! YVELF
    Real, Intent (InOut) :: HC_Temperature ! TEMI
    Integer, Intent (Out) :: Last_HC_Species
! slmod begin
! This looks to have been changed in the calling routine in hc_inside_neutral.f90, so I am guessing
! that changing it here is the correct fix...
!
    Integer, Dimension (Number_H_Species), Intent (InOut) :: Last_H_Isotope_Composition
!
!    Integer, Dimension (Number_H_Products), Intent (InOut) :: Last_H_Isotope_Composition
! slmod end
    Integer, Intent (InOut) :: Launch_Reg ! M
    Integer, Intent (InOut) :: Reflection_Counter
    Integer, Intent (In) :: NeutType
    Real, Intent (InOut) :: Max_Velocity_Randoms ! RMAXS
    Integer, Intent (In) :: Tries_Without_Change
    Double Precision, Intent (InOut) :: Neutral_Time_Step_Count
    Double Precision, Intent (InOut) :: Ion_Time_Step_Count
    Double Precision, Intent (InOut) :: Eq_Total_Ion_Time_Steps ! CIST
    Double Precision, Intent (In) :: Seed ! SEED
    Integer, Intent (InOut) :: NRand ! NRAND
    Integer, Intent (Out) :: IFate ! IFATE
    Logical, Intent (In) :: Debug

    type(hc_velocity_type1) :: hc_v

 
    ! Declare local variables.
    Integer :: Target_Index ! ID
    Real :: Sput_Weight
    Real :: Dep_Energy
    Real :: Sheath_Fraction ! Output from call to promptdep.
    Integer :: Return_Code ! Output RC from call to promptdep.
    Real :: Reflection_Probability ! REFPROB
    Real :: Sputtering_Probability
    Real :: Random_Value ! RAN
    Real :: Random_Value_1 ! RAN1
    Real :: Random_Value_2 ! RAN2
    Real :: Best ! BEST
    Real :: DSQ ! DSQ
    Integer :: IM_Index ! IM
    Real :: Current_Angle
    Real :: Angle_Offset
    Real :: Reflection_Angle
    Real :: Sputtering_Angle
    Real :: HC_Angle_To_Normal
    Real :: Segment_Normal_Angle
    Logical :: Reflect_Ion
    Logical :: Sputter_Ion
 
    Integer :: Yield_Calc = 1
    Real :: SputNew
    Real :: RYield
    Real :: EMax
 
    Integer, External :: Verify_ID ! Finds closest target segment to R,Z position.
    Real, External :: YIELD ! Finds chemical sputtering yield.
 
    ! Assign local variables.
    Sput_Weight = Calc_Sputy (Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition))
    Current_Velocity = Current_Velocity_In_S		
    !write (0,*) "Trying to promptly redeposit:",Tries_Without_Change,Cur_HC_Spec
    ! Note, PROMPT is set to redepoosit ions only after their
    ! first timestep once the transition occurs.
    If (( Prompt_Deposition_Opt .eq. 1 .or. &
         &    Prompt_Deposition_Opt .eq. 2) .and. &
         &   Tries_Without_Change .eq. 0) Then
 
       ! Option activated.
       ! Check for prompt deposition
       Call promptdep (Current_Cell,Current_Ring,Target_Index,Current_R,Current_Z,REAL (Get_HC_Charge (Cur_HC_Spec)),Sput_Weight,&
                      &Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition),HC_Temperature,Sheath_Fraction,Return_Code)
 
       ! A return code of 1 indicates that prompt redeposition
       ! has occurred. The routine also returns the impact energy
       ! of the depositing ion.
 
       If (Return_Code .eq. 1) Then
 
          ! REFLECT ---- IONS may need to be reflected here
          ! Set default non reflection condition
 
          Reflect_Ion = .False.
          Sputter_Ion = .False.
          Reflection_Probability = 0.0
 
          ! Find closest vessel segment.
          If (Current_S .le. gksmaxs (Current_Ring) / 2.0) Then
             Current_S = 0.0
             Current_Cell = 1
             Target_Index = Verify_ID (Current_Cell, Current_Ring, 2)
          Else
             Current_S = gksmaxs (Current_Ring)
             Current_Cell = gnks (Current_Ring)
             Target_Index = Verify_ID (Current_Cell, Current_Ring, 1)
          End If
 
          ! Find angle of magnetic field.
          Current_Angle = Find_Ion_Angle (Current_Cell,Current_Ring)
 
          ! Check if reflection or sputtering is allowed.
          If (hc_ion_reflection_option .eq. 1 .or. hc_sputtering_option .eq. 1) Then
             ! Select random number.
             NRand = NRand + 1
             Call Surand2 (Seed, 1, Random_Value)
 
             ! Find angle to normal at impacting vessel segment.
             HC_Angle_To_Normal = Find_Angle_To_Normal (Current_Cell, Current_Ring, Current_Angle)
 
             ! Check to see how likely reflection is from HC_Stick module.
             Call HC_Reflection_Coefficient (Target_Index, Cur_HC_Spec, HC_Temperature, HC_Angle_To_Normal, Reflection_Probability)
             If (hc_ion_reflection_option .eq. 1 .and. Random_Value .le. Reflection_Probability) Then
                Reflect_Ion = .True.
             End If
 
             ! If reflection does not occur, check to see if HC sputters.
             Call HC_Sputtering_Coefficient (Target_Index, Cur_HC_Spec, HC_Temperature, HC_Angle_To_Normal, Sputtering_Probability)
 
             If (.not. Reflect_Ion) Then
                If (hc_sputtering_option .eq. 1) Then
                   If (Random_Value .le. (Reflection_Probability + Sputtering_Probability)) Then
                      Sputter_Ion = .True.
                   End If
                End If
             End If
          End If
 
          ! Test for reflection and sputtering.
          If (Reflect_Ion) Then
             ! Reflection from intersecting segment.
             ! Set various quantities related to reflection coefficients
             !write (0,*) "Reflecting prompt ion",NProd,Cur_HC_Spec,Current_R,Current_Z
 
             ! Particle is reflected back into plasma.  Find new species and update all statistics.
             Call HC_Ion_Impact_ReLaunch ("Reflect",Cur_HC_Spec,Last_HC_Species,H_Isotope_Composition,Last_H_Isotope_Composition,&
                  & Launch_Reg, &
                  & Sput_Weight,Neutral_Time_Step_Count,Ion_Time_Step_Count,Eq_Total_Ion_Time_Steps,Reflection_Counter,NeutType,&
                  & Target_Index,Current_Theta, &
                  & Current_R,Current_Z,HC_Temperature,Current_Velocity,Current_Velocity_In_S,Current_Angle,Seed,NRand,&
                  & Current_Cell,Current_Ring, &
                  & Current_Velocity_In_R,Current_Velocity_In_Z,Current_S,Current_Cross,hc_v)
 
             if (debug_hc) then
                Write (Output_Unit_Scratch,*) 'HC PROMPT ION REFLECTED:',Current_Cell,Current_Ring,Target_Index,Current_R,Current_Z,&
                  & Target_Index,ryield, &
                  & Target_Material,Dep_Energy,Sput_Weight,Segment_Normal_Angle,Current_Theta,HC_Temperature,Current_Velocity,&
                  & Current_Angle
             endif

             ! If output is requested, write to files.
             If (hc_evolve_print_option .eq. 1) Then
                Write (Output_Unit_Evolve,9500) "Prompt ion reflect. NeutType:",NeutType,"Launch:",Launch_Reg,"Vessel:",&
                     & Target_Index,"R:",Current_R,"Z:", &
                     & Current_Z,"Last HC:",Last_HC_Species,"New HC:",Cur_HC_Spec,"New Ang:",Current_Angle,"New Vel:",&
                     & Current_Velocity,"New T:",HC_Temperature
 9500           Format(3X,A,I2,1X,A,1X,I2,1X,A,1X,I3,1X,A,1X,F6.3,1X,A,1X,F6.3,&
                      &1X,A,1X,I2,1X,A,1X,I2,1X,A,1X,F6.3,1X,A,1X,F10.3,1X,A,1X,F6.3)
             End If
             If (hc_coord_print_option .eq. 1) Then
                Write (Output_Unit_Location,9501) "Prompt ion reflect. NeutType:",NeutType,"Launch:",Launch_Reg,"Vessel:",&
                     & Target_Index,"R:",Current_R,"Z:", &
                     & Current_Z,"Last HC:",Last_HC_Species,"New HC:",Cur_HC_Spec,"New Ang:",Current_Angle,"New Vel:",&
                     & Current_Velocity,"New T:",HC_Temperature
9501            Format(3X,A,I2,1X,A,1X,I2,1X,A,1X,I3,1X,A,1X,F6.3,1X,A,1X,F6.3,&
                      &1X,A,1X,I2,1X,A,1X,I2,1X,A,1X,F6.3,1X,A,1X,F10.3,1X,A,1X,F6.3)
             End If
 
             ! Continue following the newly reflected particle.
             Return
          Else
             ! ION is NOT reflected.  Record deposition and check for self-sputtering.
             ! Record all the regular loss statistics.
             HC_Num_Absorbed_TargWall (Cur_HC_Spec,Launch_Reg) =  HC_Num_Absorbed_TargWall (Cur_HC_Spec,Launch_Reg) + Sput_Weight
             HC_Min_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) = MIN ( HC_Min_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg), REAL (&
                                                               & Eq_Total_Ion_Time_Steps))
             HC_Max_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) = MAX ( HC_Max_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg), REAL (&
                                                               & Eq_Total_Ion_Time_Steps))
             HC_Tot_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) + &
                                                                & Eq_Total_Ion_Time_Steps * Sput_Weight
 
             ! Recording addition to total elapsed time in state Cur_HC_Spec.
             HC_Tot_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) + &
                                                                 & HC_Temperature * Sput_Weight
             HC_Tot_Velocity_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Velocity_At_Absorption (Cur_HC_Spec,Launch_Reg) + &
                                                                     & Current_Velocity * Sput_Weight
             HC_Tot_ABS_Vel_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_ABS_Vel_At_Absorption (Cur_HC_Spec,Launch_Reg) + ABS (&
                                                                    & Current_Velocity) * Sput_Weight
             HC_Tot_Elec_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Elec_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) + &
                                                                      & gktebs (Current_Cell, Current_Ring) * Sput_Weight
 
             IM_Index = MIN (INT (HC_Temperature / (0.2 *  Init_Electron_Temperature)) + 1, 10)
             ! jdemod - array use does not match declaration - add second index
             !HC_Ctexs (IM_index)  =  HC_Ctexs (IM_index) + HC_Temperature * Sput_Weight
             HC_Ctexs (IM_index,Launch_reg)  =  HC_Ctexs (IM_index,Launch_reg) + HC_Temperature * Sput_Weight
             HC_Num_Absorbed_Act_Target(Cur_HC_Spec,Launch_Reg)=HC_Num_Absorbed_Act_Target(Cur_HC_Spec,Launch_Reg)+Sput_Weight
             HC_RDep (Cur_HC_Spec,Launch_Reg) =  HC_RDep (Cur_HC_Spec,Launch_Reg) + Sput_Weight
 
             If (Target_Index .lt. 1 .or. Target_Index .gt.  Num_Target_Cells) Then
                Write (Output_Unit_HC_Alert,*) 'HC DEPS Error:',Target_Index,Cur_HC_Spec,Sput_Weight,Launch_Reg
                Write (Output_Unit_HC_Alert,*) "Program stopping."
                Stop
             Else
                HC_Deposit (Target_Index,Cur_HC_Spec) =  HC_Deposit (Target_Index,Cur_HC_Spec) + Sput_Weight
                HC_Erosion (Target_Index,1,Cur_HC_Spec) =  HC_Erosion (Target_Index,1,Cur_HC_Spec) + Sput_Weight
             End If
 
             ! Prompt deposition energy is returned by the promptdep
             ! subroutine - it is calculated as some fraction of the
             ! 3 kTe sheath/MPS drop - there are no contributions
             ! for ion velocity. This is because the ion should
             ! have had no time for energy transfer collisions and
             ! should simply have its' own initial energy.
             ! Note: V^2, so do not worry about -velocity.
             Dep_Energy = Sheath_Fraction * REAL (Get_HC_Charge (Cur_HC_Spec)) * gktebs(Current_Cell,Current_Ring) + 5.22E-9 * &
                        & Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) &
                        & * Current_Velocity * Current_Velocity + 2.0 * HC_Temperature
             !Dep_Energy = Sheath_Fraction * REAL (Get_HC_Charge (Cur_HC_Spec)) * gktebs(Current_Cell,Current_Ring) + 5.22E-9 * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) &
             !& * Current_Velocity /  Ion_Time_Step * Current_Velocity /  Ion_Time_Step + 2.0 * HC_Temperature
 
             !write (0,*) "PROMPT DATA:",Dep_Energy,Sheath_Fraction,Current_Velocity,HC_Temperature,Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition)
             ! Record average energy.
             HC_PromptDeps (Target_Index,5,Cur_HC_Spec) =  HC_PromptDeps (Target_Index,5,Cur_HC_Spec) + Sput_Weight * Dep_Energy
 
             If (Yield_Calc .eq. 1 .and. Sputter_Ion) Then
                ! Particle is sputtered back into plasma.  Find new species and update all statistics.
                Call HC_Ion_Impact_ReLaunch ("Sputter",Cur_HC_Spec,Last_HC_Species,H_Isotope_Composition,&
                     & Last_H_Isotope_Composition,Launch_Reg, &
                     & Sput_Weight,Neutral_Time_Step_Count,Ion_Time_Step_Count,Eq_Total_Ion_Time_Steps,Reflection_Counter,NeutType,&
                     & Target_Index,Current_Theta, &
                     & Current_R,Current_Z,HC_Temperature,Current_Velocity,Current_Velocity_In_S,Current_Angle,Seed,NRand,&
                     & Current_Cell,Current_Ring, &
                     & Current_Velocity_In_R,Current_Velocity_In_Z,Current_S,Current_Cross,hc_v)
 
                if (debug_hc) then 
                   Write (Output_Unit_Scratch,*) "HC PROMPT ION SPUTTERED:",Current_Cell,Current_Ring,Target_Index,Current_R,&
                     & Current_Z,gwallindex(Target_Index),ryield, &
                     & Target_Material,Dep_Energy,Sput_Weight,Segment_Normal_Angle,Current_Theta,HC_Temperature,Current_Velocity,&
                     & Current_Angle
                endif
 
                ! If output is requested, write to files.
                If (hc_evolve_print_option .eq. 1) Then
                   Write (Output_Unit_Evolve,9502) "Prompt ion sputter. NeutType:",NeutType,"Launch:",Launch_Reg,"Vessel:",&
                        & Target_Index,"R:",Current_R,"Z:", &
                        & Current_Z,"Last HC:",Last_HC_Species,"New HC:",Cur_HC_Spec,"New Ang:",Current_Angle,"New Vel:",&
                        & Current_Velocity,"New T:",HC_Temperature
9502               Format(3X,A,I2,1X,A,1X,I2,1X,A,1X,I3,1X,A,1X,F6.3,&
                         &1X,A,1X,F6.3,1X,A,1X,I2,1X,A,1X,I2,1X,A,1X,F6.3,1X,A,1X,F10.3,1X,A,1X,F6.3)
                End If
                If (hc_coord_print_option .eq. 1) Then
                   Write (Output_Unit_Location,9503) "Prompt Ion sputter. NeutType:",NeutType,"Launch:",Launch_Reg,"Vessel:",&
                        & Target_Index,"R:",Current_R,"Z:", &
                        & Current_Z,"Last HC:",Last_HC_Species,"New HC:",Cur_HC_Spec,"New Ang:",Current_Angle,"New Vel:",&
                        & Current_Velocity,"New T:",HC_Temperature
9503               Format(3X,A,I2,1X,A,1X,I2,1X,A,1X,I3,1X,A,1X,F6.3,&
                         &1X,A,1X,F6.3,1X,A,1X,I2,1X,A,1X,I2,1X,A,1X,F6.3,1X,A,1X,F10.3,1X,A,1X,F6.3)
                End If
 
                ! New position, energy/velocity, and direction are loaded.  Return to following code.
                Write (Output_Unit_Scratch,*) 'HC SPUTTERED from prompt deposit:',Current_Cell,Current_Ring,Target_Index,Current_R,&
                      & Current_Z,gikds(Target_Index),girds(Target_Index),ryield,gkmfss(Target_Index),&
                      & Target_Material,Dep_Energy,Sput_Weight
                Return
 
             ElseIf (Yield_Calc .eq. 2 .and. Sputter_Ion) Then
                ! DIVIMP yield determination.
                If (gkmfss(Target_Index) .ge. 0.0) Then
                   RYield = YIELD (6,  Target_Material,Dep_Energy,gktebs (Current_Cell,Current_Ring),gktibs (Current_Cell,&
                           &Current_Ring)) * gkmfss (Target_Index)
                ElseIf (gkmfss (Target_Index) .lt. 0.0 .and. gkmfss (Target_Index) .ge. -50.0) then
                   RYield = ABS (gkmfss (Target_Index))
                ElseIf (gkmfss(Target_Index) .le. -99.0) Then
                   RYield=YIELD(6,Target_Material,Dep_Energy,gktebs(Current_Cell,Current_Ring),gktibs(Current_Cell,Current_Ring))
                End If
 
                SputNew = Sput_Weight * RYield
                HC_YldTot (Cur_HC_Spec,Launch_Reg) =  HC_YldTot (Cur_HC_Spec,Launch_Reg) + SputNew
                HC_YldMax (Cur_HC_Spec,Launch_Reg) = MAX ( HC_YldMax (Cur_HC_Spec,Launch_Reg), SputNew)
 
                Write (Output_Unit_Scratch,*) 'HC SPUTTERED:',Current_Cell,Current_Ring,Target_Index,Current_R,Current_Z,gikds (&
                     & Target_Index),girds (Target_Index),ryield,gkmfss (Target_Index), &
                     & Target_Material,Dep_Energy,Sput_Weight,SputNew, HC_YldTot
 
                If (SputNew .gt.  Self_Sputter_Threshold) Then
                   HC_YThTot (Cur_HC_Spec,Launch_Reg) =  HC_YThTot (Cur_HC_Spec,Launch_Reg) + SputNew
                   NProd  = NProd + 1
                   If (hc_launch_angle_velocity .eq. 1 .or. hc_launch_angle_velocity .eq. 4 .or. hc_launch_angle_velocity .eq. 5 &
                       &.or.  Sputter_Opt .eq. 4) Then
                      EMax =  EMax_Factor * Dep_Energy
                      Max_Velocity_Randoms = 1.0 / (1.0 +  Target_Binding_Energy / EMax)**2.0
                   Else
                      Max_Velocity_Randoms = 1.0
                   End If
 
                   Call pxprods (NProd,Current_R)
                   Call pyprods (NProd,Current_Z)
                   Call psnews(NProd,SputNew)
 
                   ! For segments with a fixed sputtering yield - allow for
                   ! the energy of the sputtered particle to be set to a
                   ! specific value.
 
                   If ( Self_Sputter_Opt .eq. 2 .and. gkmfss (Target_Index) .lt. 0.0) Then
                      ! Re-emit with temperature on an ion impact re-release as only ions will be prompt-deposited.
                      Call peprods (NProd,hc_sput_energy_ion_preset)
                   Else
                      Call peprods (NProd,0.0)
                   End If
 
                   Call pidprods (NProd,Target_Index)
                   HC_Launchdat (NProd,2) = 1.0
 
                End If
 
             End If
 
             ! Exit due to prompt deposition.
             IFate = 26 ! HC Ion Prompt Deposition.
 
             Return
          End If
 
       End If
 
    End If
 
  End Subroutine HC_Prompt_Deposition
 
End Module HC_Prompt
