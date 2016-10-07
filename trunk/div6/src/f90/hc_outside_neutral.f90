! -*-Mode:f90-*-
! Title: Go Molecular Transport Routine, Adam McLean, UTIAS, October 2001.
! Purpose: Transport a molecule within an arbitrary shape.
! Use: Used with the go executable shell script
! Notes:
 
Module HC_Outside_Neutral
 
  ! Required modules.
 
  ! Every good Fortran program has...
  Implicit None	
 
Contains
 
  Subroutine HC_Trans_Outside_Neutral (&
       &	IProd,				&
       &	LPRod,				&
       &	NProd,				&
       &	Seed,				&
       &	NRand,  			&
       &	Sput_Weight,			& ! SPUTY
       &	HC_Temperature,			& ! TEMN
       &	Launch_Reg,			& ! M
       &	NeutType,			&
       &	Tries_Without_Change,		&
       &	Eq_Total_Ion_Time_Steps,	& ! CIST
       &	Neutral_Time_Step_Count, 	&
       &	Ion_Time_Step_Count, 		&
       &	Reflection_Counter,		& ! NRFCNT
       &	MTC_Counter,			& ! MTCCNT
       &	All_Wall_R_Points,		& ! RW
       &	All_Wall_Z_Points,		& ! RZ
       &	Cur_HC_Spec,		&
       &	Last_HC_Species,		&
       &	H_Isotope_Composition,		&
       &	Last_H_Isotope_Composition,	&
       &	Current_R,			& ! R
       &	Current_Z,			& ! Z
       &	Current_S,			& ! S
       &	Current_Cross,			& ! CROSS
       &	Current_Cell,			& ! IK
       &	Current_Ring,			& ! IR
       &	Last_R,				& ! ROLD
       &	Last_Z,				& ! ZOLD
       &	S_Start,			& ! SSTART
       &	Max_Velocity_Randoms,		& ! RMAXS
       &	Current_Theta,			& ! THETA
       &	Current_Dir_Indicator,		& ! IS
       &	Current_Angle,			& ! TNEW
       &	Current_Velocity,		& ! VIN
       &	Current_Velocity_In_R,		& ! XVELF
       &	Current_Velocity_In_Z,		& ! YVELF
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
    Use HC_Get
    Use HC_Put
    Use HC_Vessel_Interact ! Sputtering and reflection coefficient calculation.
    Use HC_Utilities ! Find surface temperature.
    Use HC_Stick ! Refl. and Sput. coefficients.
    use hc_velocity_type ! HC velocity type
 
    ! Every good Fortran routine should have...
    Implicit None
 
    ! Define input/output variables
    Integer, Intent (In) :: IProd
    Integer, Intent (In) :: LPRod
    Integer, Intent (InOut) :: NProd
    Double Precision, Intent (In) :: Seed
    Integer, Intent (InOut) :: NRand
    Real, Intent (InOut) :: Sput_Weight
    Real, Intent (Out) :: HC_Temperature ! TEMN
    Integer, Intent (InOut) :: Launch_Reg ! M
    Integer, Intent (In) :: NeutType
    Integer, Intent (InOut) :: Tries_Without_Change
    Double Precision, Intent (InOut) :: Eq_Total_Ion_Time_Steps ! CIST
    Double Precision, Intent (InOut) :: Neutral_Time_Step_Count
    Double Precision, Intent (InOut) :: Ion_Time_Step_Count
    Integer, Intent (InOut) :: Reflection_Counter ! NRFCNT
    Integer, Intent (InOut) :: MTC_Counter
    Real, Intent (In) :: All_Wall_R_Points ( Max_Points)
    Real, Intent (In) :: All_Wall_Z_Points ( Max_Points)
    Integer, Intent (InOut) :: Cur_HC_Spec
    Integer, Intent (InOut) :: Last_HC_Species
    Integer, Dimension (Number_H_Species), Intent (InOut) :: H_Isotope_Composition
    Integer, Dimension (Number_H_Species), Intent (InOut) :: Last_H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (InOut) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (InOut) :: Last_H_Isotope_Composition
    Real, Intent (InOut) :: Current_R ! R
    Real, Intent (InOut) :: Current_Z ! Z
    Real, Intent (Out) :: Current_S
    Real, Intent (Out) :: Current_Cross
    Integer, Intent (InOut) :: Current_Cell ! IK
    Integer, Intent (InOut) :: Current_Ring ! IR
    Real, Intent (InOut) :: Last_R ! ROLD
    Real, Intent (InOut) :: Last_Z ! ZOLD
    Real, Intent (Out) :: S_Start ! SSTART
    Real, Intent (InOut) :: Max_Velocity_Randoms ! RMAXS
    Real, Intent (InOut) :: Current_Theta ! THETA
    Integer, Intent (InOut) :: Current_Dir_Indicator ! IS
    Real, Intent (Out) :: Current_Angle ! TNEW
    Real, Intent (InOut) :: Current_Velocity ! VIN
    Real, Intent (InOut) :: Current_Velocity_In_R ! XVELF
    Real, Intent (InOut) :: Current_Velocity_In_Z ! YVELF
    Integer, Intent (Out) :: IFate

    type (hc_velocity_type1) :: hc_v ! HC velocity structure

    ! jdemod - debug_hc_neut is a name conflict with the comhc module variable - changed name to just debug for now
    Logical, Intent (In) :: Debug
 
    ! Define local variables that will be assigned values.
 
    Integer :: ID
    Integer :: Ind
    Integer :: Wall_Index ! Wall segment number crossed.
    Integer :: Local_Neutral_Reflection_Opt ! Local NRFOPT.
    Integer :: Yield_Calc = 1
 
    Real :: Reflection_Angle ! TNEW, Reflection angle off normal returned from INTCALC.
    Real :: Segment_Normal_Angle ! TNORM, Normal angle from the intersecting wall segment, returned by INTCALC.
    Real :: DSQ
    Real :: Best
    Real :: RNew
    Real :: ZNew
    Real :: RTmp
    Real :: ZTmp
    Real :: Random_Value
    Real :: Random_Value_1
    Real :: Random_Value_2
    Real :: RefProb
    Real :: Move_Factor
    Real :: Reflection_Probability ! Probability taken from wallpt(25,x).
    Real :: Sputtering_Probability
    Real :: RYield
    Real :: Dep_Energy
    Real :: SputNew
    Real :: EMax
    Real :: Surface_Temperature ! Temperature taken from wallpt (19,x)
    Real :: HC_Angle_To_Normal ! Angle off the tangent of the surface element struck.
    Real :: Angle_Of_Impact ! Particle direction at moment of impact.
 
    Logical :: Reflect_Neut
    Logical :: Sputter_Neut
    Logical :: Intersect ! SECT
    integer :: intersect_result,intersect_index
 
    ! GA15 temporary variables.
    Integer :: GA15_IndWork (2,  Max_Points) ! INDWORK
    Real :: GA15_Work (4 *  Max_Points) ! WORK
    Real :: TDum ( Max_Points)
    Real :: XDum ( Max_Points)
    Real :: YDum ( Max_Points)
    Real :: Wall_Check_Result ! RESULT
    Real :: R_New
    Real :: Z_New
 
    Real, External :: ATAN2C ! Checks and calculates atan2 value.
    Real, External :: YIELD

    !
    ! step count - counts steps needed to see if we can get a neutral restart position in bounds
    !
    integer :: step_count
    logical :: debug_step_count = .true.

    ! iter_cnt is a debug quantity which can be uncommented and used to determine 
    ! on which iteration errors are being triggered for a given particle
    !integer,save :: iter_cnt = 0
    !iter_cnt = iter_cnt + 1

 
    If (Get_HC_Charge (Cur_HC_Spec) .eq. 0) Then
 
       ! Setup GA15 routines.
       Call GA15A ( Num_Boundary_Points,1,GA15_Work,4* Max_Points,GA15_IndWork, Max_Points,All_Wall_R_Points,All_Wall_Z_Points,&
&TDum,XDum,YDum,6)
 
       ! Neutral particle.
 
       If (Debug) Then
          Write (Output_Unit_HC_Alert,'(a,7(1x,g13.5))') 'DEBUG2:',Eq_Total_Ion_Time_Steps,Current_R,Current_Z,&
                 & Current_Velocity_In_R,Current_Velocity_In_Z,Current_Velocity,Wall_Check_Result
       End If
 
 
       ! Values of angle model of 0 - 4 coincide with NRFOPT.  If higher, need to alter for INTCALC but remembed later for reflection angle.
       If (hc_reflection_angle_model .le. 4) Then
          Local_Neutral_Reflection_Opt = hc_sputtering_angle_model
       Else
          ! Temporary assignment for use in INTCALC.
          Local_Neutral_Reflection_Opt = 1
       End If
 
 
       ! Particle is found to be outside the wall boundary.  Find the position and wall segment where it left.
       if (debug_hc) then 
          Write (Output_Unit_Scratch,'(a,2i4,2f12.7,1p,g12.4,l4)') 'WALL COLL:',Current_Cell,Current_Ring,Current_R,Current_Z,&
                     & Wall_Check_Result, Grid_Error
          Write (Output_Unit_Scratch,'(13x,i5,4f12.7,1p,7g12.4)') IProd,Last_R,Last_Z, Launch_R, Launch_Z,Current_Velocity_In_R,&
                     & Current_Velocity_In_Z,Eq_Total_Ion_Time_Steps
       endif

       ! write(6,'(a,i10,2i6,5g18.10)')  &
       !        'HC_OUTSIDE_NEUTRAL:FIND_WALL_INTERSESCTION-BEFORE:', &
       !        iter_cnt,current_cell,current_ring,current_r,current_z,last_r,last_z 

       ! jdemod - this routine handles the wall intersection calculations - since it was common to three different modules I pulled it out to a
       !          separate routine
       ! jdemod - add hc_cprint to the call to find_wall_intersection - this was modified in main code but not updated in the hc modules
       call find_wall_intersection(current_r,current_z,last_r,last_z,rnew,znew,reflection_angle,segment_normal_angle,&
                               &local_neutral_reflection_opt,wall_index,intersect_result,intersect,hc_cprint)


       !   write(6,'(a,i10,2i6,8g18.10,i10,l6)')  &
       !        'HC_OUTSIDE_NEUTRAL:FIND_WALL_INTERSESCTION-AFTER:', &
       !         iter_cnt,current_cell,current_ring,current_r,current_z,last_r,last_z,rnew,znew, &
       !         reflection_angle*raddeg,segment_normal_angle*raddeg,intersect_result,intersect
 
 
       ! Verify RNEW,ZNEW
       !			If ((ABS (Current_R - rnew) + ABS (Last_R - rnew)) .ne. ABS (Current_R - Last_R) .or. &
       !			& (ABS (Current_Z-znew) + ABS (Last_Z - znew)) .ne. ABS (Current_Z-Last_Z)) Then
       !
       !write (Output_Unit_Scratch,'(a,i5,l4,6(1x,g13.6))') 'Warning in HC_Outside_Neutral: POSSIBLE RNEW,ZNEW ERROR:',Wall_Index,Intersect,Current_R,Current_Z,rnew,znew,Last_R,Last_Z
       !write (0,'(a,i5,l4,6(1x,g13.6))') 'Warning in HC_Outside_Neutral: POSSIBLE RNEW,ZNEW ERROR:',Wall_Index,Intersect,Current_R,Current_Z,rnew,znew,Last_R,Last_Z
       !
       !				Call GA15B (Current_R,Current_Z,Wall_Check_Result, Num_Boundary_Points,1,GA15_Work,4* Max_Points, &
       !				& GA15_IndWork, Max_Points,All_Wall_R_Points,All_Wall_Z_Points,TDUM,XDUM,YDUM,6)
       ! Write (Output_Unit_Scratch,'(a,6(1x,g13.6))') 'Warning in HC_Outside_Neutral: Current_R,Current_Z,Wall_Check_Result:',Current_R,Current_Z,Wall_Check_Result
       !
       !				Call GA15B (Last_R,Last_Z,Wall_Check_Result, Num_Boundary_Points,1,GA15_Work,4* Max_Points, &
       !				& GA15_IndWork, Max_Points,All_Wall_R_Points,All_Wall_Z_Points,TDum,XDum,YDum,6)
       ! Write (Output_Unit_Scratch,'(a,6(1x,g13.6))') 'Warning in HC_Outside_Neutral: Last_R,Last_Z,Wall_Check_Result:',Current_R,Current_Z,Wall_Check_Result
       !
       !			End If
 
    Else ! Ion impact into target.  Find the wall segment it left through and the approximate R,Z position where it crossed.
 
       ! Dealt with by hc_outside_ion routine.
       Write (Output_Unit_HC_Alert,*) "Incorrect charge.  Program stopping."
       Stop
 
    End If
 
    ! Dealing with possible reflections
 
    ! Check for reflection - allowing for varying probabilities
    ! for each segment - if NRFOPT is on (ie. not zero).
 
    Reflect_Neut = .False.
    Sputter_Neut = .False.
 
    ! Check if reflection or sputtering is allowed.
    If (hc_neutral_reflection_option .eq. 1 .or. hc_sputtering_option .eq. 1) Then
       ! Select random number.
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value)
 
       ! Find angle to normal at impacting vessel segment.
       HC_Angle_To_Normal = Segment_Normal_Angle
 
       ! Check to see how likely reflection is from HC_Stick module first.
       Call HC_Reflection_Coefficient (Wall_Index,Cur_HC_Spec,HC_Temperature,HC_Angle_To_Normal,Reflection_Probability)
 
       !write (0,*) "Checking refl:",Cur_HC_Spec,hc_neutral_reflection_option,Random_Value,Reflection_Probability	
       !write (Output_Unit_Scratch,*) "Checking refl:",Cur_HC_Spec,hc_neutral_reflection_option,Random_Value,Reflection_Probability
 
       If (hc_neutral_reflection_option .eq. 1 .and. Random_Value .le. Reflection_Probability) Then
          Reflect_Neut = .True.
          Tries_Without_Change = 0
       Else
          ! If reflection does not occur, check to see if HC sputters.
          Call HC_Sputtering_Coefficient (Wall_Index,Cur_HC_Spec,HC_Temperature,HC_Angle_To_Normal,Sputtering_Probability)
 
          If (.not. Reflect_Neut .and. hc_sputtering_option .eq. 1 .and. Random_Value .le. (Reflection_Probability + &
              &Sputtering_Probability)) Then
             Sputter_Neut = .True.
             Tries_Without_Change = 0
          End If
       End If
    End If
    !write (Output_Unit_Scratch,*) "VESSEL HC neut: reflect",Reflect_Neut,"Sputter",Sputter_Neut,hc_neutral_reflection_option,hc_sputtering_option,&
    !& Random_Value,Reflection_Probability,Sputtering_Probability,Segment_Normal_Angle,HC_Angle_To_Normal
 
    ! Execute reflection or deposition.
    If (.not. Reflect_Neut) Then
 
       If (.not. Intersect) Then
          Write (Output_Unit_HC_Alert,*) 'Error in HC_Outside_Neutral: NO INTERSECTION POINT FOUND without HC neutral reflection - &
&ERROR',Current_R,Current_Z,reflect_neut
          Write (Output_Unit_Evolve,*) 'Error in HC_Outside_Neutral: NO INTERSECTION POINT FOUND without HC neutral reflection - &
&ERROR',Current_R,Current_Z,reflect_neut
          Write (Output_Unit_Location,*) 'Error in HC_Outside_Neutral: NO INTERSECTION POINT FOUND without HC neutral reflection - &
&ERROR',Current_R,Current_Z,reflect_neut
          !
          ! jdemod - remove stop condition - find wall intersection returns the nearest wall element in the case of an error
          !          so it should not be necessary to stop at this point - however, the wall intersection code needs to
          !          be revamped to reduce the chance of these errors occuring. 
          !
          ! Stop
       End If
 
       ! NEUTRAL is NOT reflected.  Record deposition and check for self-sputtering.
       ! Record all the regular loss statistics.
       HC_Num_Absorbed_TargWall (Cur_HC_Spec,Launch_Reg) =  HC_Num_Absorbed_TargWall (Cur_HC_Spec,Launch_Reg) + Sput_Weight
       HC_Min_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) = MIN ( HC_Min_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg), REAL (&
&Eq_Total_Ion_Time_Steps))
       HC_Max_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) = MAX ( HC_Max_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg), REAL (&
&Eq_Total_Ion_Time_Steps))
       HC_Tot_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) + &
&Eq_Total_Ion_Time_Steps * Sput_Weight
 
       ! Recording addition to total elapsed time in state Cur_HC_Spec.
       HC_Tot_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) + HC_Temperature * &
&Sput_Weight
       HC_Tot_Velocity_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Velocity_At_Absorption (Cur_HC_Spec,Launch_Reg) + &
&Current_Velocity * Sput_Weight
       HC_Tot_ABS_Vel_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_ABS_Vel_At_Absorption (Cur_HC_Spec,Launch_Reg) + ABS (&
&Current_Velocity) * Sput_Weight
       HC_Tot_Elec_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Elec_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) + gktebs &
&(Current_Cell, Current_Ring) * Sput_Weight
 
       !			IM_Index = MIN (INT (HC_Temperature / (0.2 *  Init_Electron_Temperature)) + 1, 10)
       !			 HC_Ctexs (IM_Index,Launch_Reg)  =  HC_Ctexs (IM_Index,Launch_Reg) + HC_Temperature * Sput_Weight
       HC_RDep (Cur_HC_Spec,Launch_Reg) =  HC_RDep (Cur_HC_Spec,Launch_Reg) + Sput_Weight
 
       ! Record neutral deposition only if particle strikes target region (NEROS (x,1)).
       If (INT (gwallpt (Wall_Index,18)) .gt. 0) Then
          HC_Erosion(INT(gwallpt(Wall_Index,18)),1,Cur_HC_Spec)=HC_Erosion(INT(gwallpt(Wall_Index,18)),1,Cur_HC_Spec)+Sput_Weight! NEROS.
          HC_Num_Absorbed_Act_Target (Cur_HC_Spec,Launch_Reg) =  HC_Num_Absorbed_Act_Target (Cur_HC_Spec,Launch_Reg) + Sput_Weight
       End If
 
       ! Check for sputtering.
       If (Yield_Calc .eq. 1 .and. Sputter_Neut) Then
 
          ! Particle is sputtered back into plasma.  Find new species and update all statistics.
          Call HC_Neutral_Impact_ReLaunch ("Sputter",Cur_HC_Spec,Last_HC_Species,H_Isotope_Composition,Last_H_Isotope_Composition,&
               &Launch_Reg, &
               & Sput_Weight,Neutral_Time_Step_Count,Ion_Time_Step_Count,Eq_Total_Ion_Time_Steps,Reflection_Counter,NeutType,&
               &Wall_Index,Reflection_Angle,Segment_Normal_Angle,Current_Theta,S_Start, &
               & RNew,ZNew,Current_R,Current_Z,HC_Temperature,Current_Velocity,Current_Angle,Seed,NRand,Current_Cell,Current_Ring, &
               & Current_Velocity_In_R,Current_Velocity_In_Z,Current_S,Current_Cross,hc_v)
 
          !Write (Output_Unit_Scratch,*) 'HC NEUTRAL SPUTTERED:',Current_Cell,Current_Ring,Wall_Index,Current_R,Current_Z,Wall_Index,ryield, &
          !&  Target_Material,Dep_Energy,Sput_Weight,Segment_Normal_Angle,Current_Theta,HC_Temperature,Current_Velocity,Current_Angle
 
          ! If output is requested, write to files.
          If (hc_evolve_print_option .eq. 1) Then
             Write (Output_Unit_Evolve,9500) "Neutral sputter. NeutType:",NeutType,"Launch:",Launch_Reg,"Vessel:",Wall_Index,"R:",&
                  &Current_R,"Z:", &
                  & Current_Z,"Last HC:",Last_HC_Species,"New HC:",Cur_HC_Spec,"New Ang:",Current_Angle,"New Vel:",&
                  &Current_Velocity,"New T:",HC_Temperature
9500 Format(3X,A,I2,1X,A,1X,I2,1X,A,1X,I3,1X,A,1X,F6.3,1X,A,1X,F6.3,1X,A,1X,I2,1X,A,1X,I2,1X,A,1X,F6.3,1X,A,1X,F10.3,1X,A,1X,F6.3)
          End If
          If (hc_coord_print_option .eq. 1) Then
             Write (Output_Unit_Location,9501) "Neutral sputter. NeutType:",NeutType,"Launch:",Launch_Reg,"Vessel:",Wall_Index,"R:&
&",Current_R,"Z:", &
                  & Current_Z,"Last HC:",Last_HC_Species,"New HC:",Cur_HC_Spec,"New Ang:",Current_Angle,"New Vel:",&
&Current_Velocity,"New T:",HC_Temperature
9501 Format(3X,A,I2,1X,A,1X,I2,1X,A,1X,I3,1X,A,1X,F6.3,1X,A,1X,F6.3,1X,A,1X,I2,1X,A,1X,I2,1X,A,1X,F6.3,1X,A,1X,F10.3,1X,A,1X,F6.3)
          End If
 
          ! Continue following the newly sputtered particle.
          Return
 
       ElseIf (Yield_Calc .eq. 2 .and. Sputter_Neut) Then
          ! DIVIMP yield determination.
          If (gkmfss(Wall_Index) .ge. 0.0) Then
             RYield = YIELD (6, Target_Material,Dep_Energy,gktebs (Current_Cell,Current_Ring),gktibs (Current_Cell,Current_Ring)) *&
& gkmfss (Wall_Index)
          ElseIf (gkmfss (Wall_Index) .lt. 0.0 .and. gkmfss (Wall_Index) .ge. -50.0) then
             RYield = ABS (gkmfss (Wall_Index))
          ElseIf (gkmfss(Wall_Index) .le. -99.0) Then
             RYield = YIELD (6, Target_Material,Dep_Energy,gktebs (Current_Cell,Current_Ring),gktibs (Current_Cell,Current_Ring))
          End If
 
          SputNew = Sput_Weight * RYield
          HC_YldTot =  HC_YldTot + SputNew
          HC_YldMax = MAX ( HC_YldMax, SputNew)
 
          !Write (Output_Unit_Scratch,*) 'HC SPUTTERED:',Current_Cell,Current_Ring,Wall_Index,Current_R,Current_Z,ryield,gkmfss (Wall_Index), Target_Material,Dep_Energy,Sput_Weight,SputNew, HC_YldTot
 
          If (SputNew .gt.  Self_Sputter_Threshold) Then
             ! Particle will be re-launched as a hydrocarbon in the
             ! self-sputtered group.  Do not include an addition to HC_Erosion
             ! since this will occur in the later call to Launch_HC.
 
             HC_YthTot (Cur_HC_Spec,Launch_Reg) =  HC_YthTot (Cur_HC_Spec,Launch_Reg) + SputNew
             NProd  = NProd + 1
             Call psnews (NPROD,SputNew)
             If (hc_launch_angle_velocity .eq. 1 .or. hc_launch_angle_velocity .eq. 4 .or. hc_launch_angle_velocity .eq. 5 .or.  &
&Sputter_Opt .eq. 4) Then
                EMax =  EMax_Factor * Dep_Energy
                Max_Velocity_Randoms = 1.0 / (1.0 +  Target_Binding_Energy / EMax)**2
             Else
                Max_Velocity_Randoms = 1.0
             End If
 
             Call pxprods (NProd,Current_R)
             Call pyprods (NProd,Current_Z)
 
             ! For segments with a fixed sputtering yield - allow for
             ! the energy of the sputtered particle to be set to a
             ! specific value.
 
             If ( Self_Sputter_Opt .eq. 2 .and. gkmfss (Wall_Index) .lt. 0.0) Then
                Call peprods (NProd,hc_sput_energy_neutral_preset)
             Else
                Call peprods (NProd,0.0)
             End If
 
             Call pidprods (NProd,Wall_Index)
             HC_Launchdat (NProd,2) = 1.0
 
          End If
 
          If (Debug) Then
             Write (Output_Unit_HC_Alert,*) 'FP HC RELAUNCH:',Current_Cell,Current_Ring,Wall_Index,Current_R,Current_Z,SputNew,&
&NProd,Dep_Energy
          End If
 
          Return
       End If
 
       ! No self-sputtering occurs.  Target collision
       If ((Wall_Index .gt.  Last_Wall_Index .and. Wall_Index .lt.  First_Trap_Index) .or. &
            &   (Wall_Index .gt.  Last_Trap_Index .and. Wall_Index .le.  Num_Wall_Points)) Then
 
          If (MTC_Counter .gt. 0) then
 
             HC_MTC_Striking_Target (Cur_HC_Spec,Launch_Reg) = &
                  &  HC_MTC_Striking_Target (Cur_HC_Spec,Launch_Reg) + Sput_Weight
          End If
 
          ! Record particles with invalid ID's in total
          If (Wall_Index .lt. 1 .or. Wall_Index .gt.  Num_Wall_Points) Then
             HC_WallsN ( Max_Points + 1,Cur_HC_Spec) = &
                  &  HC_WallsN ( Max_Points + 1,Cur_HC_Spec) + Sput_Weight
 
             If ( Launch_Wall_Index .ge. 1 .and.  Launch_Wall_Index .le.  Num_Wall_Points) Then
 
                HC_WTDep ( Launch_Wall_Index, Max_Points + 1,2) = &
                     &  HC_WTDep ( Launch_Wall_Index, Max_Points + 1,2) + Sput_Weight
 
             End If
          Else
             HC_WallsN (Wall_Index,Cur_HC_Spec) = &
                  &  HC_WallsN (Wall_Index,Cur_HC_Spec) + Sput_Weight
 
             If ( Launch_Wall_Index .ge. 1 .and.  Launch_Wall_Index .le.  Num_Wall_Points) Then
 
                HC_WTDep ( Launch_Wall_Index,Wall_Index, 2) = &
                     &  HC_WTDep ( Launch_Wall_Index,Wall_Index, 2) + Sput_Weight
             End If
          End If
 
          ! Record collision with target.
          HC_Num_Striking_Target (Cur_HC_Spec,Launch_Reg) =  HC_Num_Striking_Target (Cur_HC_Spec,Launch_Reg) + Sput_Weight
 
          IFate = 12
 
          !Write(Output_Unit_Scratch,'(a,5(1x,i6),4(1x,g12.5))') 'COLLISION WITH TARGET:',Wall_Index, &
          !&  Launch_Wall_Index,id,Current_Cell,Current_Ring,Current_R,Current_Z,gwallpt(Wall_Index,1),gwallpt(Wall_Index,2)
          Return
 
          ! Wall collision
       Else
 
          HC_Num_Reach_Wall (Cur_HC_Spec,Launch_Reg) =  HC_Num_Reach_Wall (Cur_HC_Spec,Launch_Reg) + Sput_Weight
          If (MTC_Counter.GT.0) Then
 
             HC_MTC_Reach_Wall (Cur_HC_Spec,Launch_Reg) =  HC_MTC_Reach_Wall (Cur_HC_Spec,Launch_Reg) + Sput_Weight
          End If
 
          ! Assign a value to ir that corresponds to the index of the wall segment crossed.
          If (Wall_Index .ge.  First_Wall_Index .and. Wall_Index .le.  Last_Wall_Index) Then
 
             Current_Ring =  Inner_Wall_Ring
 
          ElseIf (Wall_Index .ge.  First_Trap_Index .and. Wall_Index .le.  Last_Trap_Index) Then
 
             Current_Ring =  Upper_Trap_Ring
          Else
             Write (Output_Unit_HC_Alert,*) 'Neutral not on wall segment'//' at collision',Wall_Index,Current_Ring,Current_Cell,&
&Current_R,Current_Z
          End If
 
          HC_Walls (Current_Cell,Current_Ring,Cur_HC_Spec) = &
               &  HC_Walls (Current_Cell,Current_Ring,Cur_HC_Spec) + Sput_Weight
 
          ! Record particles with invalid ID's in total
          If (Wall_Index .lt. 1 .or. Wall_Index .gt.  Num_Wall_Points) Then
 
             HC_WallsN ( Max_Points + 1,Cur_HC_Spec) = &
                  &  HC_WallsN ( Max_Points + 1,Cur_HC_Spec) + Sput_Weight
 
             If ( Launch_Wall_Index .ge. 1 .and.  Launch_Wall_Index .le.  Num_Wall_Points) Then
 
                HC_WTDep ( Launch_Wall_Index, Max_Points + 1,2) = &
                     &  HC_WTDep ( Launch_Wall_Index, Max_Points + 1,2) + Sput_Weight
             End If
          Else
 
             HC_WallsN (Wall_Index,Cur_HC_Spec) = &
                  &  HC_WallsN (Wall_Index,Cur_HC_Spec) + Sput_Weight
 
             If ( Launch_Wall_Index .ge. 1 .and.  Launch_Wall_Index .le.  Num_Wall_Points) Then
 
                HC_WTDep ( Launch_Wall_Index,Wall_Index,2) = &
                     &  HC_WTDep ( Launch_Wall_Index,Wall_Index,2) + Sput_Weight
             End If
          End If
 
          IFate = 10
 
          !Write (Output_Unit_Scratch,*) "No-reflect Collision with wall"
          !Write (Output_Unit_Scratch,'(a,5(1x,i6),4(1x,g12.5))') ' COLLISION WITH WALL  :',Wall_Index, Launch_Wall_Index,id,Current_Cell,Current_Ring,Current_R,Current_Z
          Return
       End If
 
    ElseIf (Reflect_Neut) Then
 
       If (.not. Intersect) Then
 
          Write (Output_Unit_HC_Alert,'(a,2f14.8,l6)') 'NO INTERSECTION POINT FOUND with neutral HC reflection - ERROR',&
                       &Current_R,Current_Z,reflect_neut
          !Write (Output_Unit_Scratch,*) 'NO INTERSECTION POINT FOUND with neutral HC reflection - ERROR',Current_R,Current_Z,reflect_neut
 
          If ((Wall_Index .gt.  Last_Wall_Index .and. Wall_Index .lt.  First_Trap_Index) .or. &
               & (Wall_Index .gt.  Last_Trap_Index .and. Wall_Index .le.  Num_Wall_Points)) Then
             ! Target strike.
 
             If (MTC_Counter.gt.0) Then
 
                HC_MTC_Striking_Target (Cur_HC_Spec,Launch_Reg) =  HC_MTC_Striking_Target (Cur_HC_Spec,Launch_Reg) + Sput_Weight
             End If
 
             ! Record particles with invalid ID's in total
             If (Wall_Index .lt. 1 .or. Wall_Index .gt.  Num_Wall_Points) Then
 
                HC_WallsN ( Max_Points + 1,Cur_HC_Spec) = &
                     &  HC_WallsN ( Max_Points + 1,Cur_HC_Spec) + Sput_Weight
 
                If ( Launch_Wall_Index .ge. 1 .and.  Launch_Wall_Index .le.  Num_Wall_Points) Then
 
                   HC_WTDep ( Launch_Wall_Index, Max_Points + 1,2) = &
                        &  HC_WTDep ( Launch_Wall_Index, Max_Points + 1,2) + Sput_Weight
                End If
             Else
                HC_WallsN (Wall_Index,Cur_HC_Spec) = &
                     &  HC_WallsN (Wall_Index,Cur_HC_Spec) + Sput_Weight
 
                If ( Launch_Wall_Index .ge. 1 .and.  Launch_Wall_Index .le.  Num_Wall_Points) Then
 
                   HC_WTDep ( Launch_Wall_Index,Wall_Index,2) = &
                        &  HC_WTDep ( Launch_Wall_Index,Wall_Index,2) + Sput_Weight
                End If
             End If
             !Write (Output_Unit_Scratch,*) "Reflect Non-intersect Collision with target"
 
             ! Record collision with target.
             HC_Num_Striking_Target (Cur_HC_Spec,Launch_Reg) =  HC_Num_Striking_Target (Cur_HC_Spec,Launch_Reg) + Sput_Weight
 
             IFate = 12
             Return
          Else
             ! Wall strike.
 
             HC_Num_Reach_Wall (Cur_HC_Spec,Launch_Reg) =  HC_Num_Reach_Wall (Cur_HC_Spec,Launch_Reg) + Sput_Weight
             If (MTC_Counter .gt. 0) Then
 
                HC_MTC_Reach_Wall (Cur_HC_Spec,Launch_Reg) =  HC_MTC_Reach_Wall (Cur_HC_Spec,Launch_Reg) + Sput_Weight
             End If
 
             HC_Walls (Current_Cell,Current_Ring,Cur_HC_Spec) = &
                  &  HC_Walls (Current_Cell,Current_Ring,Cur_HC_Spec) + Sput_Weight
 
             ! Record particles with invalid ID's in total
             If (Wall_Index .lt. 1 .or. Wall_Index .gt.  Num_Wall_Points) Then
 
                HC_WallsN ( Max_Points + 1,Cur_HC_Spec) = &
                     &  HC_WallsN ( Max_Points + 1,Cur_HC_Spec) + Sput_Weight
 
                If ( Launch_Wall_Index .ge. 1 .and.  Launch_Wall_Index .le.  Num_Wall_Points) Then
 
                   HC_WTDep ( Launch_Wall_Index, Max_Points + 1,2) = &
                        &  HC_WTDep ( Launch_Wall_Index, Max_Points + 1,2) + Sput_Weight
                End If
             Else
                HC_WallsN (Wall_Index,Cur_HC_Spec) = &
                     &  HC_WallsN (Wall_Index,Cur_HC_Spec) + Sput_Weight
 
                If ( Launch_Wall_Index .ge. 1 .and.  Launch_Wall_Index .le.  Num_Wall_Points) Then
 
                   HC_WTDep ( Launch_Wall_Index,Wall_Index,2) = &
                        &  HC_WTDep ( Launch_Wall_Index,Wall_Index,2) + Sput_Weight
                End If
             End If
             !Write (Output_Unit_Scratch,*) "Reflect Non-intersect Collision with wall"
             IFate = 10
             Return
          End If
 
       Else ! Intersect.
          ! Reflection from intersecting segment.
          ! Set various quantities related to reflection coefficients
 
          !Write (0,*) "DIVIMP-HC Reflecting neutral",IProd,LPRod,NProd,Cur_HC_Spec,&
          !& Current_R,Current_Z,Reflection_Probability,Random_Value
 
          ! Particle is reflected back into plasma.  Find new species and update all statistics.
          Call HC_Neutral_Impact_ReLaunch ("Reflect",Cur_HC_Spec,Last_HC_Species,H_Isotope_Composition,Last_H_Isotope_Composition,&
               &Launch_Reg, &
               & Sput_Weight,Neutral_Time_Step_Count,Ion_Time_Step_Count,Eq_Total_Ion_Time_Steps,Reflection_Counter,NeutType,&
               &Wall_Index,Reflection_Angle,Segment_Normal_Angle,Current_Theta,S_Start, &
               & RNew,ZNew,Current_R,Current_Z,HC_Temperature,Current_Velocity,Current_Angle,Seed,NRand,Current_Cell,Current_Ring, &
               & Current_Velocity_In_R,Current_Velocity_In_Z,Current_S,Current_Cross,hc_v)
 
          !Write (Output_Unit_Scratch,*) 'HC NEUTRAL REFLECTED:',Current_Cell,Current_Ring,Wall_Index,Current_R,Current_Z,Wall_Index,ryield, &
          !&  Target_Material,Dep_Energy,Sput_Weight,Segment_Normal_Angle,Current_Theta,HC_Temperature,Current_Velocity,Current_Angle
 
          ! If output is requested, write to files.
          If (hc_evolve_print_option .eq. 1) Then
             Write (Output_Unit_Evolve,9502) "Neutral reflect. NeutType:",NeutType,"Launch:",Launch_Reg,"Vessel:",Wall_Index,"R:",&
                  &Current_R,"Z:", &
                  & Current_Z,"Last HC:",Last_HC_Species,"New HC:",Cur_HC_Spec,"New Ang:",Current_Angle,"New Vel:",&
                  &Current_Velocity,"New T:",HC_Temperature
9502 Format(3X,A,I2,1X,A,1X,I2,1X,A,1X,I3,1X,A,1X,F6.3,1X,A,1X,F6.3,1X,A,1X,I2,1X,A,1X,I2,1X,A,1X,F6.3,1X,A,1X,F10.3,1X,A,1X,F6.3)
          End If
          If (hc_coord_print_option .eq. 1) Then
             Write (Output_Unit_Location,9503) "Neutral reflect. NeutType:",NeutType,"Launch:",Launch_Reg,"Vessel:",&
                  & Wall_Index,"R:",Current_R,"Z:", &
                  & Current_Z,"Last HC:",Last_HC_Species,"New HC:",Cur_HC_Spec,"New Ang:",Current_Angle,"New Vel:",&
                  &Current_Velocity,"New T:",HC_Temperature
9503 Format(3X,A,I2,1X,A,1X,I2,1X,A,1X,I3,1X,A,1X,F6.3,1X,A,1X,F6.3,1X,A,1X,I2,1X,A,1X,I2,1X,A,1X,F6.3,1X,A,1X,F10.3,1X,A,1X,F6.3)
          End If
 
          ! Continue following the newly reflected particle.
          !Return
 

          !
          ! move_factor - controls size of corrective steps
          !
          Move_Factor = 1.0
          step_count = 0
          rtmp = rnew
          ztmp = znew

          !
          ! Check to see if the new particle position is inside 
          !

          Call GA15B (RTmp,ZTmp,Wall_Check_Result, Num_Boundary_Points,1,GA15_Work,4* Max_Points, &
                  & GA15_IndWork, Max_Points,All_Wall_R_Points,All_Wall_Z_Points,TDUM,XDUM,YDUM,6)

          !
          ! If it is not inside then iterate the particle position until it is - or we give up
          !

          if (wall_check_result.lt.0.0) then 
 
          Do
             step_count = step_count+1
             RTmp = Rtmp + Move_Factor * Current_Velocity_In_R
             ZTmp = Ztmp + Move_Factor * Current_Velocity_In_Z
 
             Call GA15B (RTmp,ZTmp,Wall_Check_Result, Num_Boundary_Points,1,GA15_Work,4* Max_Points, &
                  & GA15_IndWork, Max_Points,All_Wall_R_Points,All_Wall_Z_Points,TDUM,XDUM,YDUM,6)
             !write (0,*) "Wall_Check_Result",Wall_Check_Result
             !write (0,*) "RTMP",RNew,RTMP
             !write (0,*) "ZTMP",ZNew,ZTMP
 
             if (debug_step_count) then 
                write(6,'(a,i8,6g18.9)') 'STEP COUNT:',step_count,rtmp,ztmp,&
                                     &current_velocity_in_r,current_velocity_in_z,wall_check_result
             endif

             If (Wall_Check_Result .lt. 0) Then
 
                !Move_Factor = Move_Factor * 0.1
 
                If (step_count.gt.100) Then
 
                   Call GA15B(Last_R,Last_Z,Wall_Check_Result, Num_Boundary_Points,1,GA15_Work,4* Max_Points, &
                        & GA15_IndWork, Max_Points,All_Wall_R_Points,All_Wall_Z_Points,TDUM,XDUM,YDUM,6)
 
                   ! Write (Output_Unit_Scratch,'(a,8 (1x,g13.6))') 'WALL COLL: NEW Current_R,Current_Z OUTSIDE:', &
                   ! & RNew,ZNew,Current_R,Current_Z,Current_Velocity_In_R,Current_Velocity_In_Z,Current_Angle
                   ! Write (Output_Unit_Scratch,'(a,3i5,5 (1x,g13.6))') 'MORE1:',Current_Cell,Current_Ring, &
                   ! & Wall_Index,Last_R,Last_Z,Wall_Check_Result,gwallpt(Wall_Index,16),Current_Velocity
                   ! Write (Output_Unit_Scratch,'(a,8 (1x,g13.6))') 'MORE2:',gwallpt (Wall_Index,1),gwallpt (Wall_Index,2), &
                   ! & gwallpt (Wall_Index,20),gwallpt (Wall_Index,21),gwallpt (Wall_Index,22),gwallpt (Wall_Index,23), &
                   ! & gwallpt (Wall_Index,8),gwallpt (Wall_Index,9)
 
                   RTmp = RNew + Current_Velocity_In_R
                   ZTmp = ZNew + Current_Velocity_In_Z
 
                   Call GRIDPOS (Current_Cell,Current_Ring,Last_R,Last_Z,.false., Grid_Error)
 
                   CALL GA15B (RTmp,ZTmp,Wall_Check_Result, Num_Boundary_Points,1,GA15_Work,4* Max_Points, &
                        & GA15_IndWork, Max_Points,All_Wall_R_Points,All_Wall_Z_Points,TDUM,XDUM,YDUM,6)
 
                   Write (Output_Unit_Scratch,'(a,2i5,8 (1x,g13.6),l4)') 'MORE3:',&
                                        &Current_Cell,Current_Ring,rtmp,ztmp,Wall_Check_Result, Grid_Error
 
                   Call GRIDPOS (Current_Cell,Current_Ring,Current_R,Current_Z,.false., Grid_Error)
 
                   Call GA15B (Current_R,Current_Z,Wall_Check_Result, Num_Boundary_Points,1,GA15_Work,4* Max_Points, &
                        & GA15_IndWork, Max_Points,All_Wall_R_Points,All_Wall_Z_Points,TDUM,XDUM,YDUM,6)
 
                   Write (Output_Unit_Scratch,'(a,2i5,8 (1x,g13.6),l4)') 'MORE3:',&
                                          &Current_Cell,Current_Ring,Current_R,Current_Z,Wall_Check_Result, Grid_Error
 
                   HC_Num_Reach_Wall (Cur_HC_Spec,Launch_Reg) =  HC_Num_Reach_Wall (Cur_HC_Spec,Launch_Reg) + Sput_Weight
                   If (MTC_Counter .gt. 0) Then
                      HC_MTC_Reach_Wall (Cur_HC_Spec,Launch_Reg) =  HC_MTC_Reach_Wall (Cur_HC_Spec,Launch_Reg) + Sput_Weight
                   End If
 
                   If (Current_Ring.ne.Upper_Trap_Ring .and. Current_Ring.ne.Inner_Wall_Ring .and. (Current_Cell.ne.1.or.&
                                             &Current_Cell.ne.gnks(Current_Ring))) Then
                      Write (Output_Unit_Scratch,*) 'Neutral not in adjacent cell at collision',Current_Cell,Current_Ring
                   End If
 
                   HC_Walls (Current_Cell,Current_Ring,Cur_HC_Spec) = &
                        &  HC_Walls (Current_Cell,Current_Ring,Cur_HC_Spec) + Sput_Weight
                   ! Record particles with invalid ID's in total
 
                   If (Wall_Index .lt. 1 .or. Wall_Index .gt.  Num_Wall_Points) Then
 
                      HC_WallsN ( Max_Points + 1,Cur_HC_Spec) = &
                           &  HC_WallsN ( Max_Points + 1,Cur_HC_Spec) + Sput_Weight
 
                      If ( Launch_Wall_Index .ge. 1 .and.  Launch_Wall_Index .le.  Num_Wall_Points) Then
 
                         HC_WTDep ( Launch_Wall_Index, Max_Points + 1,2) = &
                              &  HC_WTDep ( Launch_Wall_Index, Max_Points + 1,2) + Sput_Weight
                      End If
 
                   Else
                      HC_WallsN (Wall_Index,Cur_HC_Spec) = &
                           &  HC_WallsN (Wall_Index,Cur_HC_Spec) + Sput_Weight
 
                      If ( Launch_Wall_Index .ge. 1 .and.  Launch_Wall_Index .le.  Num_Wall_Points) Then
 
                         HC_WTDep ( Launch_Wall_Index,Wall_Index,2) = &
                              &  HC_WTDep ( Launch_Wall_Index,Wall_Index,2) + Sput_Weight
                      End If
 
                   End If
                   Write (Output_Unit_HC_Alert,*) "Reflect intersect launch failure"
                   IFate = 10
                   Return
                End If
                Cycle
             Else
                Exit
             End If
          End Do

          endif
 
          !
          ! Update last particle position
          !

          last_r = current_r
          last_z = current_z

          !
          ! update to current particle position
          !
          Current_R = RTmp
          Current_Z = ZTmp
 
          ! Write (Output_Unit_Scratch,'(a,i5,1p,10g11.4)') 'NRF:',Wall_Index,RNEW,ZNEW,Last_R,Last_Z, &
          ! & Current_Angle,Current_R,Current_Z,Current_Velocity_In_R,Current_Velocity_In_Z,Eq_Total_Ion_Time_Steps
 
          If (Reflection_Counter .gt.  Max_Total_HC_Reflections) Then
 
             HC_Reflection_Loss (Cur_HC_Spec,Launch_Reg) =  HC_Reflection_Loss (Cur_HC_Spec,Launch_Reg) + Sput_Weight
 
             If ((Wall_Index .gt.  Last_Wall_Index .and. Wall_Index .lt.  First_Trap_Index) .or. &
                  & (Wall_Index .gt.  Last_Trap_Index .and. Wall_Index .le.  Num_Wall_Points)) Then
                ! Target strike.
 
                Write(Output_Unit_HC_Alert,*) 'Too many reflections:target:',Current_R,Current_Z,Launch_Reg, &
                     & Current_Cell,Current_Ring,Wall_Index,rnew,znew,Current_Velocity_In_R,Current_Velocity_In_Z
                Write(Output_Unit_HC_Alert,*)'wl:',First_Wall_Index,Last_Wall_Index,First_Trap_Index,Last_Trap_Index,Num_Wall_Points
                Write(Output_Unit_HC_Alert,*) 'wp:',(id,':',gwallpt (Wall_Index,id),',',id=1,13)
 
                Write(Output_Unit_Scratch,*) 'Too many reflections:target:',Current_R,Current_Z,Launch_Reg, &
                     & Current_Cell,Current_Ring,Wall_Index,rnew,znew,Current_Velocity_In_R,Current_Velocity_In_Z
                Write(Output_Unit_Scratch,*)'wl:',First_Wall_Index,Last_Wall_Index,First_Trap_Index,Last_Trap_Index,Num_Wall_Points
                Write(Output_Unit_Scratch,*) 'wp:',(id,':',gwallpt (Wall_Index,id),',',id=1,13)
 
                If (MTC_Counter.gt.0) Then
 
                   HC_MTC_Striking_Target (Cur_HC_Spec,Launch_Reg) =  HC_MTC_Striking_Target (Cur_HC_Spec,Launch_Reg) + Sput_Weight
                End If
 
                ! Record collision with target.
                HC_Num_Striking_Target (Cur_HC_Spec,Launch_Reg) =  HC_Num_Striking_Target (Cur_HC_Spec,Launch_Reg) + Sput_Weight
 
                IFate = 12
                Return
             Else
                ! Wall strike.
 
                Write (Output_Unit_HC_Alert,*) 'Too many reflections:wall:',Current_R,Current_Z,Launch_Reg,Current_Cell,&
&Current_Ring,Wall_Index,rnew,znew,Current_Velocity_In_R,Current_Velocity_In_Z
                !Write (Output_Unit_Scratch,*) 'Too many reflections:wall:',Current_R,Current_Z,Launch_Reg,Current_Cell,Current_Ring,Wall_Index,rnew,znew,Current_Velocity_In_R,Current_Velocity_In_Z
 
                HC_Num_Reach_Wall (Cur_HC_Spec,Launch_Reg) =  HC_Num_Reach_Wall (Cur_HC_Spec,Launch_Reg) + Sput_Weight
 
                IF (MTC_Counter.GT.0) Then
                   HC_MTC_Reach_Wall (Cur_HC_Spec,Launch_Reg) =  HC_MTC_Reach_Wall (Cur_HC_Spec,Launch_Reg) + Sput_Weight
                End If
 
                If (Current_Ring .ne.  Upper_Trap_Ring .and. Current_Ring .ne.  Inner_Wall_Ring) Then
 
                   Write (Output_Unit_HC_Alert,*) 'Neutral not in wall ring at collision',Current_Ring
 
                   If (Current_Ring .gt.  Upper_Trap_Ring) Then
                      Current_Ring =  Upper_Trap_Ring
                   Else
                      Current_Ring =  Inner_Wall_Ring
                   End If
                End If
 
                HC_Walls (Current_Cell,Current_Ring,Cur_HC_Spec) = &
                     &  HC_Walls (Current_Cell,Current_Ring,Cur_HC_Spec) + Sput_Weight
 
                ! Record particles with invalid ID's in total
                If (Wall_Index .lt. 1 .or. Wall_Index .gt.  Num_Wall_Points) Then
 
                   HC_WallsN ( Max_Points + 1,Cur_HC_Spec) = &
                        &  HC_WallsN ( Max_Points + 1,Cur_HC_Spec) + Sput_Weight
 
                   If ( Launch_Wall_Index.ge.1.and. Launch_Wall_Index.le. Num_Wall_Points) Then
 
                      HC_WTDep ( Launch_Wall_Index, Max_Points + 1,2) = &
                           &  HC_WTDep ( Launch_Wall_Index, Max_Points + 1,2) + Sput_Weight
 
                   End If
                Else
 
                   HC_WallsN (Wall_Index,Cur_HC_Spec) = &
                        &  HC_WallsN (Wall_Index,Cur_HC_Spec) + Sput_Weight
 
                   If ( Launch_Wall_Index .ge. 1 .and.  Launch_Wall_Index .le.  Num_Wall_Points) Then
 
                      HC_WTDep ( Launch_Wall_Index,Wall_Index,2) = &
                           &  HC_WTDep ( Launch_Wall_Index,Wall_Index,2) + Sput_Weight
                   End If
                End If
 
                IFate = 10
                Return
             End If
          End If
          !write (0,*) "HERE outside neutral"
       End If
    End If

 
  End Subroutine HC_Trans_Outside_Neutral
 
End Module HC_Outside_Neutral
