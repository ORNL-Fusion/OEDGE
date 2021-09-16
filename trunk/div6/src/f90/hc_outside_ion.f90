! -*-Mode:f90-*-
! Title: Go Molecular Transport Routine, Adam McLean, UTIAS, October 2001.
! Purpose: Transport a molecule within an arbitrary shape.
! Use: Used with the go executable shell script
! Notes:
 
Module HC_Outside_Ion
 
  ! Required modules.
 
  ! Every good Fortran program has...
  Implicit None	
 
Contains
 
  Subroutine HC_Trans_Outside_Ion (	&
       &	IProd,				&
       &	LPRod,				&
       &	NProd,				&
       &	Sput_Weight,			& ! SPUTY
       &	Seed,				& ! SEED
       &	NRand,  			& ! NRAND
       &	HC_Temperature,			& ! TEMN
       &	Launch_Reg,			& ! M
       &	NeutType,			&
       &	Max_Velocity_Randoms,		& ! RMAXS
       &	Tries_Without_Change,		&
       &	Eq_Total_Ion_Time_Steps,	& ! CIST
       &	Neutral_Time_Step_Count,	&
       &	Ion_Time_Step_Count,		&
       &	Reflection_Counter,		& ! NRFCNT
       &	Cur_HC_Spec,		&
       &	Last_HC_Species,		&
       &	H_Isotope_Composition,		&
       &	Last_H_Isotope_Composition,	&
       &	Current_R,			& ! R
       &	Current_Z,			& ! Z
       &	Current_Cell,			& ! IK
       &	Current_Ring,			& ! IR
       &	Last_R,				& ! ROLD
       &	Last_Z,				& ! ZOLD
       &	Current_S,			& ! S
       &	Current_Cross,			& ! CROSS
       &	Current_Theta,			& ! THETA
       &	Target_Index,			& ! ID
       &	Dir_Indicator,			& ! IS
       &	Current_Angle,			& ! TNEW
       &	Current_Velocity,		& ! VIN
       &	Current_Velocity_In_S,		& ! VINS
       &	Current_Velocity_In_R,		& ! XVELF
       &	Current_Velocity_In_Z,		& ! YVELF
       &	Random_Numbers_Used,		& ! KK
       &	IFate,				& ! IFATE
       &        hc_v,                           & ! HC Velocity structure
       &	Debug)
 
    ! External modules.
    Use ComHC ! Access user options, hydrocarbon common block.
    Use HC_Init_DIV_Data ! Gain access to cell/global/geometry/current data structures.
    Use HC_Init_DIV_Diag ! Data diagnostics.
    Use HC_Init_Lib_Data ! HC utilities.
    Use HC_Stack ! Save broken higher hydrocarbons upon reflection on stack.
    Use HC_Release ! Subroutine that finds the type of hydrocarbon to launch.
    Use HC_Put
    Use HC_Get
    Use HC_Vessel_Interact ! Sputtering and reflection coefficient calculation.
    Use HC_Utilities ! Find surface temperature.
    Use HC_Stick ! Refl. and Sput. coefficients.
    use hc_velocity_type ! jdemod - hc_v type 
 
    ! Every good Fortran routine should have...
    Implicit None
 
    ! Define input/output variables
    Integer, Intent (In) :: IProd
    Integer, Intent (In) :: LPRod
    Integer, Intent (InOut) :: NProd
    Real, Intent (InOut) :: Sput_Weight ! SPUTY
    Double Precision, Intent (In) :: Seed
    Integer, Intent (InOut) :: NRand
    Real, Intent (Out) :: HC_Temperature ! TEMN
    Integer, Intent (InOut) :: Launch_Reg ! M
    Integer, Intent (In) :: NeutType
    Real, Intent (InOut) :: Max_Velocity_Randoms ! RMAXS
    Integer, Intent (InOut) :: Tries_Without_Change
    Double Precision, Intent (InOut) :: Eq_Total_Ion_Time_Steps ! CIST
    Double Precision, Intent (InOut) :: Neutral_Time_Step_Count
    Double Precision, Intent (InOut) :: Ion_Time_Step_Count
    Integer, Intent (InOut) :: Reflection_Counter ! NRFCNT
    Integer, Intent (InOut) :: Cur_HC_Spec
    Integer, Intent (InOut) :: Last_HC_Species
    Integer, Dimension (Number_H_Species), Intent (InOut) :: H_Isotope_Composition
    Integer, Dimension (Number_H_Species), Intent (InOut) :: Last_H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (InOut) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (InOut) :: Last_H_Isotope_Composition
    Real, Intent (InOut) :: Current_R ! R
    Real, Intent (InOut) :: Current_Z ! Z
    Integer, Intent (InOut) :: Current_Cell ! IK
    Integer, Intent (InOut) :: Current_Ring ! IR
    Real, Intent (InOut) :: Last_R ! ROLD
    Real, Intent (InOut) :: Last_Z ! ZOLD
    Real, Intent (InOut) :: Current_S
    Real, Intent (InOut) :: Current_Cross
    Real, Intent (InOut) :: Current_Theta ! THETA
    Integer, Intent (Out) :: Target_Index ! ID
    Integer, Intent (InOut) :: Dir_Indicator ! IS
    Real, Intent (Out) :: Current_Angle ! TNEW
    Real, Intent (InOut) :: Current_Velocity ! VIN
    Real, Intent (InOut) :: Current_Velocity_In_S ! VINS
    Real, Intent (Out) :: Current_Velocity_In_R ! XVELF
    Real, Intent (Out) :: Current_Velocity_In_Z ! YVELF
    Integer, Intent (InOut) :: Random_Numbers_Used ! KK
    Integer, Intent (Out) :: IFate

    type(hc_velocity_type1) :: hc_v

    Logical, Intent (In) :: Debug
 
    ! Define local variables that will be assigned values.
    real :: smax_local

 
    Integer :: Ind
    Integer :: IndI ! Wall segment number crossed.
 
    Real :: Segment_Normal_Angle ! TNORM, Normal angle from the intersecting wall segment, returned by INTCALC.
    Real :: DSQ
    Real :: Best
 
    Integer :: Yield_Calc = 1
    Integer :: IM_Index ! IM
    Real :: RNew
    Real :: ZNew
    Real :: RTmp
    Real :: ZTmp
    Real :: Random_Value
    Real :: Random_Value_1
    Real :: Random_Value_2
    Real :: RefProb
    Real :: Move_Factor
    Real :: Dep_Energy
    Real :: SputNew
    Real :: RYield
    Real :: EMax
    Real :: Reflection_Probability ! Probability taken from wallpt(25,x).
    Real :: Sputtering_Probability
    Real :: Surface_Temperature ! Temperature taken from wallpt (19,x)
    Real :: HC_Angle_To_Normal ! Angle off the tangent of the surface element struck.
 
    Logical :: Reflect_Neut
    Logical :: Reflect_Ion
    Logical :: Sputter_Ion
 
    ! GA15 temporary variables.
    Integer :: GA15_IndWork (2, MaxPts) ! INDWORK
    Real :: GA15_Work (4 * MaxPts) ! WORK
    Real :: TDum (MaxPts)
    Real :: XDum (MaxPts)
    Real :: YDum (MaxPts)
 
    Logical :: Intersect ! SECT
    Logical :: Wall_Check_Result ! RESULT
 
    Real :: Adjust ! Not currently used but required to be declared.
    Real :: DCross(4)
 
    Integer, External :: Verify_ID
    Real, External :: YIELD
 
    ! Ion impact into walls or target.  Find the wall segment it left through and the approximate R,Z position where it crossed.
    ! Check if reached target.
    If (Current_S .le. 0.0 .or. Current_S .ge. gksmaxs (Current_Ring)) Then
       If ( Target_Mirror_Opt .eq. 0) Then
 
 
          !
          ! jdemod - What is this doing here? Why stop the code for S values > 100.0?????
          !
          !write (0,*) "HC_Ion_Target_Location",Current_S,Current_Cross,Current_Cell,Current_Ring,Current_R,Current_Z,Target_Index
          !If (Current_S .gt. 100.0 .or. Current_R .gt. 10.0) Then
          !	Write (Output_Unit_HC_Alert,*) "Serious positional error with S position in HC_Outside_Ion.  Stopping."
          !	Stop
          !End If
          Call HC_Ion_Target_Location (Current_S,Current_Cross,Current_Cell,Current_Ring,Current_R,Current_Z,Target_Index)
 
          ! Reflect first.
          Reflect_Ion = .False.
          Sputter_Ion = .False.
          Reflection_Probability = 0.0
 
          ! Find angle of magnetic field.
          Current_Angle = Find_Ion_Angle (Current_Cell,Current_Ring)
 
          ! Check if reflection is allowed.
          If (hc_ion_reflection_option .eq. 1 .or. hc_sputtering_option .eq. 1) Then
             ! Select random number.
             NRand = NRand + 1
             Call Surand2 (Seed, 1, Random_Value)
 
             ! Find angle to normal at impacting vessel segment.
             HC_Angle_To_Normal = Find_Angle_To_Normal (Current_Cell,Current_Ring,Current_Angle)
 
             ! Check to see how likely reflection is from HC_Stick module first.
             Call HC_Reflection_Coefficient (Target_Index, Cur_HC_Spec,HC_Temperature,HC_Angle_To_Normal,Reflection_Probability)
             If (hc_ion_reflection_option .eq. 1 .and. Random_Value .le. Reflection_Probability) Then
                Reflect_Ion = .True.
                Tries_Without_Change = 0
             Else
                ! If reflection does not occur, check to see if HC sputters.
                Call HC_Sputtering_Coefficient (Target_Index,Cur_HC_Spec,HC_Temperature,HC_Angle_To_Normal,Sputtering_Probability)
 
                If (.not. Reflect_Ion .and. hc_sputtering_option .eq. 1 .and. Random_Value .le. (Reflection_Probability + &
                    &Sputtering_Probability)) Then
                   Sputter_Ion = .True.
                   Tries_Without_Change = 0
                End If
             End If
          End If
          !Write (Output_Unit_Scratch,*) "VESSEL HC Ion:reflect",Reflect_Ion,"Sputter",Sputter_Ion,hc_neutral_reflection_option,hc_sputtering_option,&
          !& Random_Value,Reflection_Probability,Sputtering_Probability,Segment_Normal_Angle,HC_Angle_To_Normal
 
          ! Test for reflection and sputtering.
          If (Reflect_Ion) Then
             ! Particle is reflected.  Find new species (Note, reflected species is always neutral).
             !write (0,*) "DIVIMP-HC Reflecting ion",IProd,LPRod,NProd,Cur_HC_Spec,Current_R,Current_Z,Reflection_Probability,Random_Value
 
             ! Particle is reflected back into plasma.  Find new species and update all statistics.
             Call HC_Ion_Impact_ReLaunch ("Reflect",Cur_HC_Spec,Last_HC_Species,H_Isotope_Composition,Last_H_Isotope_Composition,&
                  & Launch_Reg, &
                  & Sput_Weight,Neutral_Time_Step_Count,Ion_Time_Step_Count,Eq_Total_Ion_Time_Steps,Reflection_Counter,NeutType,&
                  & Target_Index,Current_Theta, &
                  & Current_R,Current_Z,HC_Temperature,Current_Velocity,Current_Velocity_In_S,Current_Angle,Seed,NRand,&
                  & Current_Cell,Current_Ring, &
                  & Current_Velocity_In_R,Current_Velocity_In_Z,Current_S,Current_Cross,hc_v)
 
             !Write (Output_Unit_Scratch,*) 'HC ION REFLECTED AT TARGET:',Current_Cell,Current_Ring,Target_Index,Current_R,Current_Z,gwallindex(Target_Index),ryield, &
             !&  Target_Material,Dep_Energy,Sput_Weight,Segment_Normal_Angle,Current_Theta,HC_Temperature,Current_Velocity,Current_Angle
 
             ! If output is requested, open files.
             If (hc_evolve_print_option .eq. 1) Then
                Write (Output_Unit_Evolve,9505) "Ion reflection at target index:",Target_Index,"R:",Current_R,"Z:",Current_Z,"Last &
                                                &HC:",Last_HC_Species,"New HC:",Cur_HC_Spec
9505            Format (3X,A,1X,I3,1X,A,1X,F8.3,1X,A,1X,F8.3,1X,A,1X,I2,1X,A,1X,I2)
             End If
             If (hc_coord_print_option .eq. 1) Then
                Write (Output_Unit_Location,9506) "Ion reflection at target index:",Target_Index,"R:",Current_R,"Z:",Current_Z,&
                                                  &"Last HC:",Last_HC_Species,"New HC:",Cur_HC_Spec
9506            Format (3X,A,1X,I3,1X,A,1X,F8.3,1X,A,1X,F8.3,1X,A,1X,I2,1X,A,1X,I2)
             End If
 
             ! Continue following the newly reflected particle.

! ---- RETURN EXIT ----
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
                                                                     & Current_Velocity_In_S * Sput_Weight
             HC_Tot_ABS_Vel_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_ABS_Vel_At_Absorption (Cur_HC_Spec,Launch_Reg) + ABS (&
                                                                     & Current_Velocity_In_S) * Sput_Weight
             HC_Tot_Elec_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Elec_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) + &
                                                                      & gktebs (Current_Cell, Current_Ring) * Sput_Weight
 
             IM_Index = MIN (INT (HC_Temperature / (0.2 *  Init_Electron_Temperature)) + 1, 10)
             HC_Ctexs (IM_Index,launch_reg)  =  HC_Ctexs (IM_Index,launch_reg) + HC_Temperature * Sput_Weight
             HC_Num_Absorbed_Act_Target (Cur_HC_Spec,Launch_Reg) = HC_Num_Absorbed_Act_Target(Cur_HC_Spec,Launch_Reg) + Sput_Weight
             HC_RDep (Cur_HC_Spec,Launch_Reg) =  HC_RDep (Cur_HC_Spec,Launch_Reg) + Sput_Weight
 
             If (Target_Index .lt. 1 .or. Target_Index .gt.  Num_Target_Cells) Then
                Write (Output_Unit_HC_Alert,*) 'HC DEPS Error:',Target_Index,Cur_HC_Spec,Sput_Weight
                Write (Output_Unit_HC_Alert,*) "Program stopping."
                Stop
             Else
                HC_Deposit (Target_Index,Cur_HC_Spec) =  HC_Deposit (Target_Index,Cur_HC_Spec) + Sput_Weight ! DEPS.
                HC_Erosion (Target_Index,1,Cur_HC_Spec) =  HC_Erosion (Target_Index,1,Cur_HC_Spec) + Sput_Weight ! NEROS.
             End If
 
             Dep_Energy = 3.0 * REAL(Get_HC_Charge (Cur_HC_Spec)) * gktebs (Current_Cell,Current_Ring) + 5.22E-9 * Find_HC_Mass (&
                        & Cur_HC_Spec,H_Isotope_Composition) &
                        & * Current_Velocity /  Ion_Time_Step * Current_Velocity /  Ion_Time_Step + 2.0 * HC_Temperature
 
             If (Yield_Calc .eq. 1 .and. Sputter_Ion) Then
 
                ! Particle is sputtered back into plasma.  Find new species and update all statistics.
                Call HC_Ion_Impact_ReLaunch ("Sputter",Cur_HC_Spec,Last_HC_Species,H_Isotope_Composition,&
                     & Last_H_Isotope_Composition,Launch_Reg, &
                     & Sput_Weight,Neutral_Time_Step_Count,Ion_Time_Step_Count,Eq_Total_Ion_Time_Steps,Reflection_Counter,NeutType,&
                     & Target_Index,Current_Theta, &
                     & Current_R,Current_Z,HC_Temperature,Current_Velocity,Current_Velocity_In_S,Current_Angle,Seed,NRand,&
                     & Current_Cell,Current_Ring, &
                     & Current_Velocity_In_R,Current_Velocity_In_Z,Current_S,Current_Cross,hc_v)
 
                ! New position, energy/velocity, and direction are loaded.  Return to following code.
                !Write (Output_Unit_Scratch,*) 'HC ION SPUTTERED from target:',Current_Cell,Current_Ring,Current_Angle,Current_Velocity,Target_Index, &
                !& Current_R,Current_Z,gikds(Target_Index),girds(Target_Index),gkmfss(Target_Index), Target_Material,Dep_Energy,Sput_Weight
 
                ! If output is requested, open files.
                If (hc_evolve_print_option .eq. 1) Then
                   Write(Output_Unit_Evolve,9507)"Ion sputter.NeutType:",NeutType,"Launch:",Launch_Reg,"Target:",Target_Index,"R:",&
                        & Current_R,"Z:",Current_Z,"Last HC:",Last_HC_Species,"New HC:",Cur_HC_Spec,"New Ang:",Current_Angle,"New &
                        &Vel:",Current_Velocity,"New T:",HC_Temperature
9507               Format(3X,A,I2,1X,A,1X,I2,1X,A,1X,I3,1X,A,1X,F6.3,&
                         &1X,A,1X,F6.3,1X,A,1X,I2,1X,A,1X,I2,1X,A,1X,F6.3,1X,A,1X,F10.3,1X,A,1X,F6.3)
                End If
                If (hc_coord_print_option .eq. 1) Then
                   Write (Output_Unit_Location,9508) "    Ion sputter. NeutType:",NeutType,"Launch:",Launch_Reg,"Target:",&
                        & Target_Index,"R:", &
                        & Current_R,"Z:",Current_Z,"Last HC:",Last_HC_Species,"New HC:",Cur_HC_Spec,"New Ang:",Current_Angle,"New &
                        &Vel:",Current_Velocity,"New T:",HC_Temperature
9508               Format(3X,A,I2,1X,A,1X,I2,1X,A,1X,I3,1X,A,1X,F6.3,&
                         &1X,A,1X,F6.3,1X,A,1X,I2,1X,A,1X,I2,1X,A,1X,F6.3,1X,A,1X,F10.3,1X,A,1X,F6.3)
                End If
 
                ! Continue following the newly sputtered particle.

! ---- RETURN EXIT ----

                Return
 
             ElseIf (Yield_Calc .eq. 2 .and. Sputter_Ion) Then
                ! DIVIMP yield determination.
                If (gkmfss(Target_Index) .ge. 0.0) Then
                   RYield = YIELD (6,  Target_Material,Dep_Energy,gktebs(Current_Cell,Current_Ring),gktibs(Current_Cell,&
                                  &Current_Ring)) * gkmfss(Target_Index)
                ElseIf (gkmfss(Target_Index) .lt. 0.0 .and. gkmfss(Target_Index) .ge. -50.0) then
                   RYield = ABS(gkmfss(Target_Index))
                ElseIf (gkmfss(Target_Index) .le. -99.0) Then
                   RYield=YIELD(6,Target_Material,Dep_Energy,gktebs(Current_Cell,Current_Ring),gktibs(Current_Cell,Current_Ring))
                End If
 
                SputNew = Sput_Weight * RYield
                HC_YldTot =  HC_YldTot + SputNew
                HC_YldMax = MAX ( HC_YldMax, SputNew)
 
                !Write (Output_Unit_Scratch,*) 'HC SPUTTERED:',Current_Cell,Current_Ring,Target_Index,Current_R,Current_Z,gikds(Target_Index),girds(Target_Index), &
                !& RYield,gkmfss(Target_Index), Target_Material,Dep_Energy,Sput_Weight,SputNew, HC_YldTot
 
                If (SputNew .gt.  Self_Sputter_Threshold) Then
                   HC_YthTot =  HC_YthTot + SputNew
                   NProd = NProd + 1
                   Call psnews (NProd,SputNew)
                   If (hc_launch_angle_velocity .eq. 1 .or. hc_launch_angle_velocity .eq. 4 .or. hc_launch_angle_velocity .eq. 5 &
                       &.or.  Sputter_Opt .eq. 4) Then
                      EMax =  EMax_Factor * Dep_Energy
                      Max_Velocity_Randoms = 1.0 / (1.0 +  Target_Binding_Energy / EMax)**2.0
                   Else
                      Max_Velocity_Randoms = 1.0
                   End If
 
                   ! No adjustments for initial position or position on target options
                   ! since promptly redeposited particles are assumed to have not
                   ! travelled too far.
 
                   Current_R = grp(Target_Index)
                   Current_Z = gzp(Target_Index)
 
                   Call pxprods(NPROD,Current_R)
                   Call pyprods(NPROD,Current_Z)
 
                   ! For segments with a fixed sputtering yield - allow for
                   ! the energy of the sputtered particle to be set to a
                   ! specific value.
 
                   If ( Self_Sputter_Opt .eq. 2 .and. (gkmfss (Target_Index) .lt. 0.0 .and. gkmfss(Target_Index) .ge. -50.0)) Then
                      Call peprods (NProd, Init_Particle_Temperature)
                   Else
                      Call peprods (NProd,0.0)
                   End If
 
                   Call pidprods (NProd,Target_Index)
                   HC_Launchdat (NProd,2) = 0.0
 
                   ! Add to total removal.
                   HC_Erosion (Target_Index,3,Cur_HC_Spec) =  HC_Erosion (Target_Index,3,Cur_HC_Spec) + Sput_Weight ! NEROS.
                End If

! ---- RETURN EXIT ----
 
                Return
 
             End If
 
             ! Reflect ion and sputter ion must be false.  Ion collision with target.
             ! Add ion weight to wall element closest to grid departure.
             !
             ! jdemod - bug/error  - call uses the charge of the cur_hc_spec and not the actual species itself
             !                       which is what is used with all other references to hc_update_walldep
             !
             !Call HC_Update_WallDep (Current_Cell,Current_Ring,Get_HC_Charge (Cur_HC_Spec),Target_Index, 0,&
             !                      & Launch_Wall_Index,NeutType,Sput_Weight,Launch_Reg)
             !
             Call HC_Update_WallDep (Current_Cell,Current_Ring,Cur_HC_Spec,Target_Index, 0,&
                                   & Launch_Wall_Index,NeutType,Sput_Weight,Launch_Reg)
 

             ! Record collision with target.
             HC_Num_Striking_Target (Cur_HC_Spec,Launch_Reg) =  HC_Num_Striking_Target (Cur_HC_Spec,Launch_Reg) + Sput_Weight
             !write (0,*) "HERE: ION struck target",Current_R,Current_Z,Last_R,Last_Z
             IFate = 22 ! HC Ion Struck Target (DIV IFATE = 2).

! ---- RETURN EXIT ----

             Return
             ! Endif for test of target ion/neutral mirror.
 
          End If
 
          ! Mirror target specified
       ElseIf ( Target_Mirror_Opt .eq. 1) Then
          ! Reverse sign of the particle velocity
          ! to mimic target reflection/impact
          Current_Velocity_In_S = -Current_Velocity_In_S
          !write (0,*) "Mirroring1:",Current_S,Current_Cell,Current_Ring,-Current_Velocity_In_S,Current_Velocity_In_S,Reflection_Counter,Ion_Time_Step_Count
          If (Current_S .lt. 0) Then
             Current_S = -Current_S
          ElseIf (Current_S .gt. gksmaxs (Current_Ring)) Then
             Current_S = gksmaxs (Current_Ring) - (Current_S - gksmaxs (Current_Ring))
 
             ! For the mirror target - it is necessary to
             ! move the particle just slightly out from the
             ! target in order to avoid an infinite loop
             ! when S is exactly 0.0 or smax.
          ElseIf (Current_S .eq. 0.0) Then
             If (gkss (Current_Cell,Current_Ring) .eq. 0.0) Then
                Current_S = gkss (Current_Cell + 1,Current_Ring) / 100.0
             Else
                Current_S = gkss (Current_Cell,Current_Ring) / 100.0
             End If
          ElseIf (Current_S .eq. gksmaxs (Current_Ring)) Then
             If (gkss (Current_Cell,Current_Ring) .eq. gksmaxs (Current_Ring)) Then
                Current_S = gksmaxs (Current_Ring) - (gksmaxs (Current_Ring) - gkss (Current_Cell - 1, Current_Ring))/100.0
             Else
                Current_S = gksmaxs (Current_Ring) - (gksmaxs (Current_Ring) - gkss (Current_Cell,Current_Ring))/100.0
             End If
          End If
 
          ! Mirror target specified - confine particle - do not reflect velocity.
       ElseIf ( Target_Mirror_Opt .eq. 2) Then
          ! Do not reverse sign of the particle velocity to mimic target reflection/impact.
          !write (0,*) "Mirroring2:",Current_S,Current_Cell,Current_Ring,Current_Velocity_In_S,Current_Velocity_In_S,Reflection_Counter,Ion_Time_Step_Count
          If (Current_S .lt. 0) Then
             Current_S = -Current_S
          ElseIf (Current_S .gt. gksmaxs (Current_Ring)) Then
             Current_S = gksmaxs (Current_Ring) - (Current_S - gksmaxs (Current_Ring))
 
             ! For the mirror target - it is necessary to
             ! move the particle just slightly out from the
             ! target in order to avoid an infinite loop
             ! when S is exactly 0.0 or smax.
 
          ElseIf (Current_S .eq. 0.0) Then
             If (gkss (Current_Cell,Current_Ring) .eq. 0.0) Then
                Current_S = gkss(Current_Cell + 1, Current_Ring) / 100.0
             Else
                Current_S = gkss (Current_Cell,Current_Ring) / 100.0
             End If
          ElseIf (Current_S .eq. gksmaxs (Current_Ring)) Then
             If (gkss (Current_Cell,Current_Ring) .eq. gksmaxs (Current_Ring)) Then
                Current_S = gksmaxs (Current_Ring) - (gksmaxs (Current_Ring) - gkss (Current_Cell - 1,Current_Ring)) / 100.0
             Else
                Current_S = gksmaxs (Current_Ring) - (gksmaxs (Current_Ring) - gkss (Current_Cell,Current_Ring)) / 100.0
             End If
          End If
       End If
    Else
       ! This code should not be used if S=>0.0 or S<=SMAX.
       Write (Output_Unit_HC_Alert,*) "Error: HC_OUTSIDE_ION should not be used for particle within target boundaries:",Current_S,&
                                      &Current_Ring,gksmaxs (Current_Ring)
       Stop
    End If

    ! jdemod - NOTE: These only come into play for target mirror option cases! 
    !                It is questionable that update_cross should be run at this point. 
 
    ! Update the Cross-field transport term and make
    ! all other related position adjustments.
 
    If (Debug) Then
       Write(Output_Unit_Scratch,'(a,2i4,1p,4g12.5,l4)') 'UP_CROSS_HC:2A:',Current_Cell,Current_Ring,Current_S,Current_Theta,&
                                                         &Current_Cross,Ion_Time_Step_Count,Debug
    End If
 
    smax_local = gksmaxs (Current_Ring)
    Call Update_Cross (Current_Cell,Current_Ring, Launch_Cell, Launch_Ring,Random_Numbers_Used,Current_S,Current_Theta, &
         & Current_Cross,Adjust,Dcross, CKK_Minimum,smax_local, Local_K,nrand,iprod,eq_total_ion_time_steps,Debug)
 
    If (Debug) Then
       Write(Output_Unit_Scratch,'(a,2i4,1p,4g12.5,l4)') 'UP_CROSS:2B:',Current_Cell,Current_Ring,Current_S,Current_Theta,&
                                                         &Current_Cross,Ion_Time_Step_Count,Debug
    End If
 
    If (Current_S .lt. 0.0 .or. Current_S .gt. gksmaxs (Current_Ring)) Then
       ! This should be avoided by this code.
       Write (Output_Unit_HC_Alert,*) "Error: HC_OUTSIDE_ION should return the particle within the flux tube:",Current_S,&
                                     &Current_Ring,gksmaxs (Current_Ring)
    End If
 
    ! slmod begin.
    If ( RZ_Opt .gt. 0) Then
       Call GETRZ (Current_Cell,Current_Ring,Current_S,Current_Cross,Current_R,Current_Z, RZ_Opt)
    Else
       Current_R = grs (Current_Cell,Current_Ring)
       Current_Z = gzs (Current_Cell,Current_Ring)
    End If
    ! slmod end.
 
  End Subroutine HC_Trans_Outside_Ion
 
End Module HC_Outside_Ion
