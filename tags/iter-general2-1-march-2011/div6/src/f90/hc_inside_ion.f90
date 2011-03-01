! -*-Mode:f90-*-
! Title: Go Molecular Transport Routine, Adam McLean, UTIAS, October 2001.
! Purpose: Transport a molecule within an arbitrary shape.
! Use: Used with the go executable shell script
! Notes:  Also takes care of far-periphery ion interaction.
!         For the far-periphery, no update of NEROS or DEPS is made.

Module HC_Inside_Ion

  ! Required modules.

  ! Every good Fortran 90 program has...
  Implicit None	

Contains

  Subroutine HC_Trans_Inside_Ion (	&
       &	IProd,				&
       &	LPRod,				&
       &	NProd,				&
       &	Sput_Weight,			&
       &	LatIZ,				&
       &	NatIZ,				&
       &	Seed,				&
       &	NRand,				&
       &	HC_Temperature,			& ! TEMN
       &	Kin_Energy_Added,		&
       &	Launch_Reg,			& ! M
       &	NeutType,			&
       &	Random_Numbers_Used,		& ! KK
       &	Max_Random_Values,		& ! KKLIM
       &	Tries_Without_Change,		&
       &	Eq_Total_Ion_Time_Steps,	& ! CIST
       &	Neutral_Time_Step_Count,	&
       &	Ion_Time_Step_Count,		&
       &	MTC_Time_Step_Count,		& ! MTCCIST
       &	MTC_Counter,			& ! MTCCNT
       &	S_Start,			& ! SSTART
       &	Max_Velocity_Randoms,		& ! RMAXS
       &	Reflection_Counter,		& ! NRFCNT
       &	Cur_HC_Spec,		&
       &	Last_HC_Species,		&
       &	Particle_Ionized,		&
       &	H_Isotope_Composition,		&
       &	Last_H_Isotope_Composition,	&
       &	Current_Velocity,		& ! VIN
       &	Current_Velocity_In_R,		& ! XVELF
       &	Current_Velocity_In_Z,  	& ! YVELF
       &	Current_Velocity_In_S,		& ! VINS
       &	Current_R,			& ! R
       &	Current_Z,			& ! Z
       &	Current_Cell,			& ! IK
       &	Current_Ring,			& ! IR
       &	Current_S,			& ! S
       &	Current_Cross,			& ! CROSS
       &	Current_Angle,  		& ! Polar ANGLE theta.
       &	Azimuthal_Angle,  		& ! Phi.
       &	Current_Theta,			& ! THETA
       &	Launch_Angle,			& ! ANGLAN
       &	Current_Tangent,		& ! TANGNT
       &	Tangent_Launch,			& ! TANLAN
       &	Last_R,				&
       &	Last_Z,				&
       &	Last_Cell,			& ! IKLAST, IKOLD
       &	Last_Ring,			& ! IRLAST, IROLD
       &	Min_Z,				& ! XXX
       &	Max_S,				& ! SSS
       &	Max_S_Ion_Removal,		& ! SSSS
       &	Time_Bin,			& ! IT
       &	TStepHC,			& ! TSTEPN
       &	IFate,				& ! IFATE
       &        hc_v,                           & ! HC velocity structure
                                ! jdemod - debug_hc_ion is a name collision with the global in hc_init_div_data with the same name - changed to debug here
                                !&	Debug_HC_Ion)
       &        Debug)

    ! External modules.
    Use HC_Init_DIV_Data ! Gain access to cell/global/geometry/current data structures.
    Use HC_Init_DIV_Diag ! Data diagnostics.
    Use ComHC ! Access user options, hydrocarbon common block.
    Use HC_Get ! External functions.
    Use HC_Put ! External functions.
    ! jdemod - no reason to need access to data load routines from here
    !Use HC_LdDta ! Access hydrocarbon storage data.
    Use HC_NewSt ! Subroutines for hydrocarbon transitions.
    Use HC_Utilities
    Use HC_Init_Lib_Data ! Hydrocarbon data and utilities.
    Use HC_Vessel_Interact ! Reflection and sputtering routines.
    Use HC_Stick ! Refl. and Sput. coefficients.
    Use HC_Prompt ! Prompt deposition.
    Use HC_Freespace_Transition ! Transition physics.
    Use HC_Stack ! Load stored hydrocarbons upon completion of primary molecule.

    use hc_kinetics      ! HC kinetics routines and options
    use bfield ! magnetic field information

    ! Every good Fortran routine should have...
    Implicit None

    ! Define input/output variables
    Integer, Intent (In) :: IProd
    Integer, Intent (In) :: LPRod
    Integer, Intent (InOut) :: NProd
    Real, Intent (InOut) :: Sput_Weight ! SPUTY
    Integer, Intent (In) :: LatIZ
    Integer, Intent (InOut) :: NatIZ
    Double Precision, Intent (In) :: Seed
    Integer, Intent (InOut) :: NRand
    Real, Intent (InOut) :: HC_Temperature ! TEMN
    Real, Intent (InOut) :: Kin_Energy_Added
    Integer, Intent (InOut) :: Launch_Reg ! M
    Integer, Intent (In) :: NeutType
    Integer, Intent (InOut) :: Random_Numbers_Used ! KK
    Integer, Intent (In) :: Max_Random_Values ! KKLIM
    Integer, Intent (InOut) :: Tries_Without_Change
    ! jdemod - eq_total_ion_time_steps - changed to InOut because something changes it
    Double Precision, Intent (InOut) :: Eq_Total_Ion_Time_Steps ! CIST
    Double Precision, Intent (InOut) :: Neutral_Time_Step_Count
    Double Precision, Intent (InOut) :: Ion_Time_Step_Count
    Double Precision, Intent (InOut) :: MTC_Time_Step_Count ! MTCCIST
    Integer, Intent (InOut) :: MTC_Counter
    Real, Intent (InOut) :: S_Start ! SSTART
    Real, Intent (InOut) :: Max_Velocity_Randoms ! RMAXS
    ! jdemod - reflection_counter - changed to InOut because something changes it
    Integer, Intent (InOut) :: Reflection_Counter ! NRFCNT
    Integer, Intent (InOut) :: Cur_HC_Spec
    Integer, Intent (InOut) :: Last_HC_Species
    Integer, Intent (InOut) :: Particle_Ionized
    Integer, Dimension (Number_H_Species), Intent (InOut) :: H_Isotope_Composition
    Integer, Dimension (Number_H_Species), Intent (InOut) :: Last_H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (InOut) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (InOut) :: Last_H_Isotope_Composition
    Real, Intent (InOut) :: Current_Velocity ! VIN
    Real, Intent (InOut) :: Current_Velocity_In_R ! XVELF
    Real, Intent (InOut) :: Current_Velocity_In_Z ! YVELF
    Real, Intent (InOut) :: Current_Velocity_In_S ! VINS
    Real, Intent (InOut) :: Current_R ! R
    Real, Intent (InOut) :: Current_Z ! Z
    Integer, Intent (InOut) :: Current_Cell ! IK
    Integer, Intent (InOut) :: Current_Ring ! IR
    Real, Intent (InOut) :: Current_S ! S
    Real, Intent (InOut) :: Current_Cross ! CROSS
    Real, Intent (InOut) :: Current_Angle ! Polar ANGLE theta.
    Real, Intent (InOut) :: Azimuthal_Angle ! Phi.
    Real, Intent (InOut) :: Current_Theta ! THETA
    Real, Intent (In) :: Launch_Angle ! ANGLAN
    Real, Intent (In) :: Current_Tangent ! TANGNT
    Real, Intent (In) :: Tangent_Launch ! TANLAN
    Real, Intent (InOut) :: Last_R !
    Real, Intent (InOut) :: Last_Z !
    Integer, Intent (InOut) :: Last_Cell ! IKLAST
    Integer, Intent (InOut) :: Last_Ring ! IRLAST
    Real, Dimension (Number_HC_Species,Number_Regions) , Intent (InOut) :: Min_Z ! XXX
    Real, Dimension (Number_HC_Species,Number_Regions) , Intent (InOut) :: Max_S ! SSS
    Real, Dimension (Number_HC_Species,Number_Regions) , Intent (InOut) :: Max_S_Ion_Removal ! SSSS
    Integer, Intent (InOut) :: Time_Bin ! IT
    Real, Intent (InOut) :: TStepHC ! TSTEPN
    Integer, Intent (Out) :: IFate ! IFATE

    type (hc_velocity_type1) :: hc_v   ! HC Velocity

    Logical, Intent (In) :: Debug


    ! Define local variables that will be assigned values.
    Integer :: ISOL ! Ionizdat first parameter ISOL (SOL,1 or Main,2).
    Integer :: IFP ! Ionizdat third parameter IFP (1 for C+ from regular launch, 2 for C+ from a far periphery relaunch.
    Integer :: IRFLCT ! Ionizdat fourth parameter IRFLCT (1 for C+ that has not been reflected, 2 for C+ that has been reflected).
    Integer :: INDI ! INDI, Temporary index value used to find wall segments.
    Integer :: IND ! IND, Temporary index value used to find wall segments.
    Integer :: Target_Index ! ID, wall/target index.
    Integer :: Temp_Cell ! JK
    Integer :: IM_Index ! IM
    Integer :: FP_Return ! RES, Return code from FPERIPH.
    Integer :: Followed_State ! Keep track of carbon-containing products in a breakup.
    Integer :: Add_Followed_State ! Keep track of additional carbon-containing products in a breakup.
    Integer  :: Old_Cell
    Integer :: Old_Ring
    Integer :: Yield_Calc = 1
    Integer :: Launch_Target

    Real :: Temp_Velocity_Factor ! TMPVELF
    Real :: Temp_Current_Velocity ! VINTMP
    Real :: Temp_Current_Velocity_In_R ! XVELFTMP
    Real :: Temp_Current_Velocity_In_Z ! YVELFTMP
    Real :: HC_Timestep ! Equals either DIVIMP ion or neutral timestep, depending on hydrocarbon.
    Real :: Random_Temp ! Temporary random number.
    Double Precision :: Total_Steps ! CIST
    Real :: Temp_Theta ! THETA1
    Real :: FP_Loss_Time ! fplosstim
    Real :: FP_Distance ! fpdist
    Double Precision :: FP_TimeSteps ! cistfp
    Real :: SputNew
    Real :: RYield
    Real :: EMax
    Real :: Reflection_Probability
    Real :: Sputtering_Probability
    Real :: Reflection_Angle
    Real :: Sputtering_Angle
    Real :: Angle_Offset
    Real :: Dep_Energy
    Real :: Segment_Normal_Angle
    Real :: HC_Angle_To_Normal
    Real :: Time_Tot
    Real :: Temp_Vel
    Real :: Old_S
    Real :: Old_Cross
    Real :: Old_R
    Real :: Old_Z
    Real :: Old_Theta
    Real :: Old_Angle
    Real :: Old_Velocity
    Real :: Old_Velocity_In_S
    Real :: Old_Temp
    Logical :: Reflect_Ion
    Logical :: Sputter_Ion
    Double Precision :: Double_Sput_Weight ! DSPUTY

    ! Variables for hydrocarbon change of state.
    Integer :: Counter
    Integer :: Max_Without_Change
    ! jdemod - only load states possible after the final state has been chosen
    !Integer, dimension (Max_Reactions_Per_State, Highest_Carbon_Content) :: States_Possible
    Integer, dimension (Highest_Carbon_Content) :: States_Possible
    Integer, dimension (Max_Reactions_Per_State) :: Reactions_Possible
    ! jdemod - this array total_probatilities is not used 
    !Real, dimension (Max_Reactions_Per_State) :: Total_Probabilities
    Real, dimension (Max_Reactions_Per_State) :: NS_Probabilities
    Real, dimension (Max_Reactions_Per_State) :: Reaction_Probabilities
    Real :: Random_Value
    Real :: Random_Value_1
    Real :: Random_Value_2
    Real :: Total_Probability

    ! jdemod - Variables supporting FP option 5
    real :: rsect,zsect,hc_mass,ctemav
    integer :: id_out,is_out,iz,istate

    ! jdemod - number of possible reactions for hc in the current state
    integer :: number_reactions
    integer :: in ! generic counting variable

    ! External functions.
    Integer, External :: FPERIPH
    Integer, External :: Verify_ID
    Integer, External :: IPOS
    Real, External :: YIELD


    Total_Steps = Neutral_Time_Step_Count + Ion_Time_Step_Count
    Double_Sput_Weight = DBLE (Sput_Weight)

    ! Take care of ion inside or moving inside main plasma.
    If (Current_Ring .lt.  Inner_SOL_Ring) Then

       ! Accumulate statistics if particle has not entered the core plasma previously.
       If (.not.  HC_InCore) Then
          HC_InCore = .True.
          HC_Num_Entered_Core (Cur_HC_Spec, Launch_Reg) =  HC_Num_Entered_Core (Cur_HC_Spec, Launch_Reg) + Sput_Weight
          !write (0,*) "LAUNCH IK,IR:", Launch_Cell,  Launch_Ring
          HC_Ion_Core_Density ( Launch_Cell,  Launch_Ring) = &
               &  HC_Ion_Core_Density ( Launch_Cell,  Launch_Ring) + Sput_Weight

          If (NeutType .gt. 0) Then
             ! Record occurance of leakage.
             HC_WTSource ( Launch_Wall_Index, Launch_Ring,3,NeutType) = &
                  &  HC_WTSource ( Launch_Wall_Index, Launch_Ring,3,NeutType) + Sput_Weight
             !write (Output_Unit_Scratch,*) "SOURCING4:", HC_WTSource ( Starting_Index, Launch_Ring,3,NeutType),NeutType,Cur_HC_Spec
             ! Record leakage location.
             HC_WTSource ( Launch_Wall_Index, Launch_Ring,4,NeutType) = &
                  &  HC_WTSource ( Launch_Wall_Index, Launch_Ring,4,NeutType) + Sput_Weight * S_Start
          End If
       End If

       ! Check if ion has moved into main plasma.
       If (Last_Ring .ge.  Inner_SOL_Ring) Then
          ! Record particle has entered core plasma - count it and save
          ! it's starting ionization position. (for possible scatter plot)
          If (.not.  HC_Has_Leaked_Core) Then
             HC_Has_Leaked_Core = .True.
             Number_Leaked_Core (Launch_Reg) =  Number_Leaked_Core (Launch_Reg) + 1
             Total_Leaked_Core (Launch_Reg) =  Total_Leaked_Core (Launch_Reg) + Sput_Weight
             If ( Number_Leaked_Core (Launch_Reg) .le. Max_Impurities) Then ! Note:  Max_Impurities is found in ComHC.
                HC_Leak_Position ( Number_Leaked_Core (Launch_Reg),Launch_Reg,1) =  Launch_R
                HC_Leak_Position ( Number_Leaked_Core (Launch_Reg),Launch_Reg,2) =  Launch_Z
             End If
          End If
       End If

       If ( HC_CflRin) Then
          ! Nonorth
          If (Debug) Then
             Write (Output_Unit_Scratch,9003) IProd,Total_Steps,Current_Cell,Current_Ring,Cur_HC_Spec,Current_R,Current_Z,&
                  &Current_S, Local_K,Current_Theta, &
                  & gksmaxs (Current_Ring),Current_Velocity,HC_Temperature,Current_Cross,Sput_Weight,Time_Bin,'HC entered main'
          End If
          ! Nonorth

          ! Current position statistics upon core entry.
          HC_Num_Orig_Enter_Core (Cur_HC_Spec,Launch_Reg)  =  HC_Num_Orig_Enter_Core (Cur_HC_Spec,Launch_Reg)  + Sput_Weight
          HC_Tot_Z_Orig_Enter_Core(Cur_HC_Spec,Launch_Reg)=HC_Tot_Z_Orig_Enter_Core(Cur_HC_Spec,Launch_Reg)+Current_Z*Sput_Weight
          HC_Tot_S_Orig_Enter_Core (Cur_HC_Spec,Launch_Reg) =  HC_Tot_S_Orig_Enter_Core (Cur_HC_Spec,Launch_Reg) + MIN (Current_S,&
               & gksmaxs (Current_Ring) - Current_S) * Sput_Weight

          ! Accumulate some statistics on different sources.
          ! Original Neutral from FP launch.

          If ( HC_LaunchDat (IProd, 2) .eq. 1.0) Then
             ! Original neutral refected.

             If ( HC_LaunchDat (IProd, 3) .ne. 0.0) Then
                HC_Tot_Core_From_FP_Ref (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Core_From_FP_Ref (Cur_HC_Spec,Launch_Reg) + Sput_Weight
             Else
                HC_Tot_Core_From_FP_NoRef(Cur_HC_Spec,Launch_Reg)=HC_Tot_Core_From_FP_NoRef(Cur_HC_Spec,Launch_Reg)+Sput_Weight
             End If
          Else
             ! Original neutral refected.
             If ( HC_LaunchDat (IProd, 3) .ne. 0) Then
                HC_Tot_Core_From_Reg_Ref (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Core_From_Reg_Ref (Cur_HC_Spec,Launch_Reg) + Sput_Weight
             Else
                HC_Tot_Core_From_Reg_NoRef(Cur_HC_Spec,Launch_Reg)=HC_Tot_Core_From_Reg_NoRef(Cur_HC_Spec,Launch_Reg)+Sput_Weight
             End If
          End If

          HC_Tot_Teq_Orig_Enter_Core (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Teq_Orig_Enter_Core (Cur_HC_Spec,Launch_Reg) +  &
               & Ion_Time_Step * Total_Steps * Sput_Weight
          HC_ELims (Current_Cell, 3, Cur_HC_Spec) =  HC_ELims (Current_Cell, 3, Cur_HC_Spec) + Sput_Weight

          HC_Tot_Teq_Orig_Enter_Main (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Teq_Orig_Enter_Main (Cur_HC_Spec,Launch_Reg) + Total_Steps &
               &  * Sput_Weight
          HC_Tot_Temp_Orig_Enter_Main (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Temp_Orig_Enter_Main (Cur_HC_Spec,Launch_Reg) + &
               & HC_Temperature * Sput_Weight
          HC_Min_Teq_Orig_Enter_Main (Cur_HC_Spec,Launch_Reg) = MIN ( HC_Min_Teq_Orig_Enter_Main (Cur_HC_Spec,Launch_Reg), Real(&
               &Total_Steps))
          HC_CflRin = .FALSE.

          ! Stop following particle in main - if option is set.
          If ( Stop_Ion_In_Core_Opt .eq. 1) Then
             HC_Stopped_Follow_Ion_In_Core(Cur_HC_Spec,Launch_Reg)=HC_Stopped_Follow_Ion_In_Core(Cur_HC_Spec,Launch_Reg)+Sput_Weight
             IFate = 21
             Return
          End If

          HC_Tot_Z_At_Main_Entry(Cur_HC_Spec,Launch_Reg)=HC_Tot_Z_At_Main_Entry(Cur_HC_Spec,Launch_Reg)+Current_Z*Sput_Weight
       End If

       ! Check for reflection off central mirror.
       If (Current_Ring .le.  Core_Ring) Then
          Current_Cross = -ABS (Current_Cross)
          !Write (Output_Unit_Scratch,*) 'WARNING: Adjusting HC particle found inside central mirror:',Current_Cell,Current_Ring,Current_Cross,Current_R,Current_Z,Current_S

          ! Set to one step outside core mirror ring.
          Current_Ring =  Core_Ring
          Current_Cross = -gdistout (Current_Cell,  Core_Ring) - gkperps (Current_Cell,  Core_Ring)
       End If
       If ( HC_CflRxa) Then
          HC_Num_Orig_Central_Reflect (Cur_HC_Spec,Launch_Reg) =  HC_Num_Orig_Central_Reflect (Cur_HC_Spec,Launch_Reg) + Sput_Weight
          HC_Tot_Teq_Orig_Central_Reflect (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Teq_Orig_Central_Reflect (Cur_HC_Spec,Launch_Reg) + &
               & Total_Steps * Sput_Weight
          HC_Tot_Temp_At_Central_Reflect (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Temp_At_Central_Reflect (Cur_HC_Spec,Launch_Reg) + &
               & HC_Temperature * Sput_Weight
          HC_Min_Teq_Orig_Central_Reflect (Cur_HC_Spec,Launch_Reg) = MIN ( HC_Min_Teq_Orig_Central_Reflect (Cur_HC_Spec,Launch_Reg)&
               &, Real(Total_Steps))
          HC_CflRxa = .False.
          If ( Debug) Then
             Write (Output_Unit_Scratch,9003) IProd,Total_Steps,Current_Cell,Current_Ring,Cur_HC_Spec,Current_R,Current_Z,&
                  & Current_S, Local_K, &
                  & Current_Theta,gksmaxs (Current_Ring),Current_Velocity,HC_Temperature,Current_Cross,Sput_Weight,Time_Bin,'HC &
                  &reflected'
          End If
       End If

       ! Take care of ion in SOL or trapped plasma.
    Else ! IR .ge. IRSEP.

       ! Accumulate some data if the particle enters various regions
       ! of the SOL or trapped plasma. (related to initial ionization
       ! positions of particles entering various regions).
       If (.not.  HC_InEdge) Then
          HC_InEdge = .True.
          HC_Ion_Edge_Density ( Launch_Cell, Launch_Ring) = &
               &  HC_Ion_Edge_Density ( Launch_Cell, Launch_Ring) + Sput_Weight
       End If

       If (.not.  HC_InTrap .and. Current_Ring .gt.  Inner_Wall_Ring .and. Current_Ring .le.  Num_Upper_Rings) Then
          ! TRAP region
          HC_InTrap = .True.
          HC_Ion_Trap_Density ( Launch_Cell, Launch_Ring) = &
               &  HC_Ion_Trap_Density ( Launch_Cell, Launch_Ring) + Sput_Weight
       ElseIf (.not.  HC_InDIV .and. Current_Ring .ge.  Inner_SOL_Ring .and. Current_Ring .le.  Inner_Wall_Ring .and. &
            & ((Current_Z .ge.  Z_X_Point .and.  ReFCT .eq. 1) .or. (Current_Z .le.  Z_X_Point .and.  ReFCT .eq. 0))) then
          ! Divertor region
          HC_InDIV = .True.
          HC_Ion_Divertor_Density ( Launch_Cell, Launch_Ring) = &
               &  HC_Ion_Divertor_Density ( Launch_Cell, Launch_Ring) + Sput_Weight
       ElseIf (.not.  HC_InMSOL .and. Current_Ring .ge.  Inner_SOL_Ring .and. Current_Ring .le.  Inner_Wall_Ring .and. &
            & ((Current_Z .lt.  Z_X_Point .and.  ReFCT .eq. 1) .or. (Current_Z .gt.  Z_X_Point .and.  ReFCT .eq. 0))) Then
          ! Main SOL Region.
          HC_InMSOL = .True.
          HC_Ion_MSOL_Density ( Launch_Cell, Launch_Ring) = &
               &  HC_Ion_MSOL_Density ( Launch_Cell, Launch_Ring) + Sput_Weight
       End If

       ! Check if HC reached walls.
       ! MODIFY TO SUPPORT ION REFLECTION INSTEAD OF ELIMINATION IF THE OPTION IS ACTIVE.
       ! OPTION 2 MOVES THE WALL TO THE OUTERMOST RING INSTEAD OF MIDWAY BETWEEN THE LAST AND NEXT TO LAST.

       If ( Non_Orthogonal_Grid_Opt .eq. 1 .or.  Non_Orthogonal_Grid_Opt .eq. 3) Then
          ! Non-orthogonal
          Last_Cell = Current_Cell
          Last_Ring = Current_Ring
       End If

       If (Current_Ring .eq.  Inner_Wall_Ring .or. Current_Ring .eq.  Upper_Trap_Ring .or. Current_Ring .eq.  Outer_Wall_Ring .or. &
            &Current_Ring .eq.  Lower_Trap_Ring) Then

          !write (0,*) "HC Ran into IRWALL", Inner_Wall_Ring, Outer_Wall_Ring, Far_Periphery_Opt,Current_Cell,Current_Ring,Current_s,Current_Cross

          If ((( Far_Periphery_Opt .eq. 0 .or. ( Far_Periphery_Opt .eq. 4 .and. Current_Ring .eq.  Inner_Wall_Ring)) .and. &
               & ( Ion_Wall_Opt .eq. 0 .or.  Ion_Wall_Opt .eq. 2)) .or. &
               & (( Far_Periphery_Opt .eq. 0 .or. ( Far_Periphery_Opt.eq.4.and.Current_Ring.eq. Inner_Wall_Ring)) .and.  &
               & Ion_Wall_Opt .eq. 1 .and. &
               & ((Current_Cross .le. 0.0 .and. (Current_Ring .eq.  Inner_Wall_Ring .or. Current_Ring .eq.  Outer_Wall_Ring) ).or. &
               &(Current_Cross.ge.0.0.and.(Current_Ring.eq.Upper_Trap_Ring.or.Current_Ring.eq.Lower_Trap_Ring)))))Then

             ! Record all the regular loss statistics.
             HC_Num_Absorbed_TargWall (Cur_HC_Spec,Launch_Reg) =  HC_Num_Absorbed_TargWall (Cur_HC_Spec,Launch_Reg) + Sput_Weight
             HC_Min_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) = MIN ( HC_Min_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg), Real(&
                  &Total_Steps))
             HC_Max_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) = MAX ( HC_Max_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg), Real(&
                  &Total_Steps))
             HC_Tot_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) + Total_Steps *&
                  & Sput_Weight

             ! Recording addition to total elapsed time in state Cur_HC_Spec.
             HC_Tot_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) + &
                  &HC_Temperature * Sput_Weight
             HC_Tot_Velocity_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Velocity_At_Absorption (Cur_HC_Spec,Launch_Reg) + &
                  &Current_Velocity_In_S * Sput_Weight
             HC_Tot_ABS_Vel_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_ABS_Vel_At_Absorption (Cur_HC_Spec,Launch_Reg) + ABS (&
                  &Current_Velocity_In_S) * Sput_Weight
             HC_Tot_Elec_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Elec_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) + &
                  &gktebs (Current_Cell, Current_Ring) * Sput_Weight

             IM_Index = MIN (INT (HC_Temperature / (0.2 *  Init_Electron_Temperature)) + 1, 10)
             HC_Ctexs (IM_Index,Launch_Reg)  =  HC_Ctexs (IM_Index,Launch_Reg) + HC_Temperature * Sput_Weight

             HC_RWall (Cur_HC_Spec,Launch_Reg) =  HC_RWall (Cur_HC_Spec,Launch_Reg) + Sput_Weight
             HC_Walls (Current_Cell,Current_Ring,Cur_HC_Spec) = &
                  &  HC_Walls (Current_Cell,Current_Ring,Cur_HC_Spec) + Sput_Weight

             ! Add ion weight to wall element closest to grid departure.
             Call HC_Update_WallDep (Current_Cell,Current_Ring,Cur_HC_Spec,0,0,Launch_Wall_Index,NeutType,Sput_Weight,Launch_Reg)

             IFate = 20 ! HC Ion Reached Vessel Wall (DIV IFATE = 1).
             Return

          ElseIf((Far_Periphery_Opt.eq.1.or.(Far_Periphery_Opt.eq.4.and.Current_Ring.eq.Upper_Trap_Ring)).and.Ion_Wall_Opt&
               &.eq.1.and.((Current_Cross.le.0.0.and.(Current_Ring.eq.Inner_Wall_Ring.or.Current_Ring.eq.Outer_Wall_Ring)).or.&
               & (Current_Cross .ge. 0.0 .and. (Current_Ring .eq.  Upper_Trap_Ring .or. Current_Ring .eq.  Lower_Trap_Ring)))) Then

             ! REFLECT ION AND CONTINUE.
             ! ALLOW ION TO CONTINUE IT WILL NEVER BE FARTHER THAN THE OUTERMOST RING.
             !write (0,*) "Reversing HC cross1 and continuing",-Current_Cross,Current_Cell,Current_Ring
             Current_Cross = -Current_Cross

          ElseIf (( Far_Periphery_Opt .eq. 1 .or. ( Far_Periphery_Opt .eq. 4 .and. Current_Ring .eq.  Upper_Trap_Ring)) &
               & .and. ( Ion_Wall_Opt .eq. 0 .or.  Ion_Wall_Opt .eq. 2)) Then

             ! REFLECT ION AND CONTINUE.
             ! ALLOW ION TO CONTINUE IT WILL NEVER BE FARTHER THAN
             ! THE INNER EDGE OF THE OUTERMOST RING REGION.
             ! FOR THETA ALLOW FOR (REMOTE) POSSIBILITY OF CROSSING
             ! SEPARATRIX IN THE INNER DIVERTOR REGION.
             !write (0,*) "Reversing HC cross2 and continuing",-Current_Cross,Current_Cell,Current_Ring
             Current_Cross = -Current_Cross
             If (Current_Ring .eq.  Inner_Wall_Ring) Then
                Temp_Cell = gikins (Current_Cell,Current_Ring)
                Current_Ring = girins (Current_Cell,Current_Ring)

                If ( Non_Orthogonal_Grid_Opt .eq. 1 .or.  Non_Orthogonal_Grid_Opt .eq. 3) Then
                   ! nonorth
                   If (Current_Ring .eq.  Num_Upper_Rings .and. Temp_Cell .ge.  Cell_TI) Then
                      Current_Theta = Current_Theta -  DThetaG
                   End If
                   ! nonorth
                End If

             ElseIf (Current_Ring .eq.  Upper_Trap_Ring) Then
                Temp_Cell = gikouts (Current_Cell,Current_Ring)
                Current_Ring = girouts (Current_Cell,Current_Ring)

                If ( Non_Orthogonal_Grid_Opt .eq. 1 .or.  Non_Orthogonal_Grid_Opt .eq. 3) Then
                   ! nonorth
                   If (Current_Ring .eq.  Inner_SOL_Ring .and. Temp_Cell .gt.  Cell_TO) Then
                      Current_Theta = Current_Theta +  DThetaG
                   End If
                   ! nonorth
                End If
             End If

             ! nonorth
             If ( Non_Orthogonal_Grid_Opt .eq. 0 .or.  Non_Orthogonal_Grid_Opt .eq. 2 .or. (( Non_Orthogonal_Grid_Opt .eq. 1 .or. &
                  &  Non_Orthogonal_Grid_Opt .eq. 3) .and. (gtagdv (Temp_Cell,Current_Ring) .eq. 0 .and. gtagdv (Last_Cell,&
                  &Last_Ring) .eq. 0))) Then
                ! BOTH POINTS UNSHIFTED FROM ORTHOGONAL GRID
                ! nonorth
                Current_Cell = Temp_Cell
                ! nonorth
             ElseIf ( Non_Orthogonal_Grid_Opt .eq. 1 .or.  Non_Orthogonal_Grid_Opt .eq. 3) Then
                ! EITHER POINT SHIFTED FROM ORTHOGONAL GRID
                If (Last_Ring .eq.  Inner_Wall_Ring) Then
                   Temp_Cell = giking (Current_Cell,Last_Ring)
                   Current_Ring = girins (Current_Cell,Last_Ring)
                ElseIf (Last_Ring .eq.  Upper_Trap_Ring) Then
                   Temp_Cell = gikoutg (Current_Cell,Last_Ring)
                   Current_Ring = girouts (Current_Cell,Last_Ring)
                End If

                ! Infinite loop bug - for temp_cell=1 and current_theta < thetag(temp_cell,current_ring)

                Do
                   If (Current_Theta .lt. gthetag (Temp_Cell,Current_Ring)) Then
                      If (Temp_Cell .gt. 1) Then
                         Temp_Theta = (gthetag (Temp_Cell,Current_Ring) - Current_Theta) / (gthetag (Temp_Cell,Current_Ring) - &
                              &gthetag (Temp_Cell - 1,Current_Ring))
                         If (Temp_Theta .gt. 0.5) Then
                            Temp_Cell = Temp_Cell - 1
                         Else
                            Exit
                         End If
                      else
                         exit
                      End If
                   Else
                      If (Temp_Cell .lt. gnks (Current_Ring)) Then
                         Temp_Theta = (Current_Theta - gthetag (Temp_Cell,Current_Ring)) / (gthetag (Temp_Cell + 1,Current_Ring) - &
                              &gthetag (Temp_Cell,Current_Ring))
                         If (Temp_Theta .ge. 0.5) Then
                            Temp_Cell = Temp_Cell + 1
                         Else
                            Exit
                         End If
                      else
                         exit
                      End If
                   End If
                End Do

                Current_Cell = Temp_Cell
             End If
             ! nonorth

             !
             ! jdemod - add support for far periphery option 5 - allowing for transport in FP
             !
          ElseIf((Far_Periphery_Opt.eq.3.or.far_periphery_opt.eq.5).and.&
               &((Ion_Wall_Opt.eq.0.or.Ion_Wall_Opt.eq.2).or.(Ion_Wall_Opt.eq.1.and.&
               & ((Current_Cross .le. 0.0 .and. (Current_Ring .eq.  Inner_Wall_Ring .or. Current_Ring .eq.  Outer_Wall_Ring)) .or. &
               & (Current_Cross .ge. 0.0 .and. (Current_Ring .eq.  Upper_Trap_Ring.or.Current_Ring.eq. Lower_Trap_Ring))) ) )) Then

             ! FPOPT 3 AND WALL OPT 1 OR 3:
             ! WALL AT LAST RING

             ! WALL IS REPLACED BY THE FAR PERIPHERY REGION.
             ! IONS HAVE 4 POSSIBLE FATES:
             ! 1) RE-ENTRY WHERE THEY EXITED - IF THEY DIFFUSE BACK.
             ! 2) DISCARDED BY REACHING TIME-LIMIT.
             ! 3) IMPACT WITH WALL AT DISTANCE FPXMAX.
             ! 4) FP LOSS TO TARGET PLATE IMPACT. USING
             ! CHARACTERISTIC TIME FPTIM.

             ! FPOPT 3 AND WALL OPT 0 OR 2:
             ! AS ABOVE BUT PLACING WALL/PERIPHERY RING BETWEEN
             ! THE LAST TWO RINGS. THIS IS NECESSARY BECAUSE THE
             ! JET SHOT DATA DOES NOT CONTAIN INFORMATION
             ! ON THE BACKGROUND CHARACTERISTICS FOR THE OUTERMOST
             ! RINGS - DESPITE CONTAINING ALL THE GEOMETRY DATA.
             !write (0,*) "Entering HC far periphery", Far_Periphery_Opt, Ion_Wall_Opt,Current_S,Current_CROSS,Total_Steps
             HC_FPEnt (Cur_HC_Spec,Launch_Reg) =  HC_FPEnt (Cur_HC_Spec,Launch_Reg) + Sput_Weight

             ! Outer or Inner FP
             If (Current_S .lt. (gksmaxs (Current_Ring) / 2.0)) Then
                ! OUTER
                FP_Distance =  hc_FPXMaxO
                FP_Loss_Time =  Ion_Time_Step /  FPTimeO
             Else
                ! INNER
                FP_Distance =  hc_FPXMaxI
                FP_Loss_Time =  Ion_Time_Step /  FPTimeI
             End If
             !write (0,*) "FP DATA HC:", FPXMaxI, FPXMaxO, Ion_Time_Step, FPTimeI, FPTimeO



             if (far_periphery_opt.eq.5) then 
                !
                !                For regular ions istate = iz - for hc fragments they differ
                !
                iz = Get_HC_Charge (Cur_HC_Spec)
                HC_mass = Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition)
                istate = Cur_HC_Spec
                ctemav = gctemav() 
                !
                !                call periphery transport routine

                ! Krieger IPP/07 - SUN compiler insists on 132 column limit
                if (debug_hc) & 
                     & write(6,'(a,i5,5g15.5)') 'HC ENTERING FP:',fp_return,hc_temperature,current_s,current_velocity_in_s, &
                     &                          total_steps

                !
                ! jdemode - note - current_velocity_in_s is not the same type of quantity as in 
                ! the main DIVIMP code - the timestep is NOT implicitly carried through and must
                ! be used explicitly here
                current_velocity_in_s = current_velocity_in_s * ion_time_step


                ! jdemod - at the present time the hc_v velocity structure is not passed into the far periphery - this only has
                !          an impact on the temperature of the particle and the evolution of the perpendicular velocity 
                !          component since the code does not assume that the parallel velocity remains constant - vpara is
                !          used to calculate the particle velocity when the hc kinetics are invoked. However, the question
                !          remains what to do about tperp, vperp and differing values of hc_temperature from time spent in the 
                !          fp - in particular the plasma conditions in the fp should probably not be represented by the 
                !          values from the outermost ring - otherwise the particles may heat too much.
                !          For the time being the best option is probably to stop tperp and vperp from evolving in the FP. 
                !          

                call fp_transport(iprod,current_cell,current_ring,iz,istate,current_s,current_theta,&
                     & current_cross,&
                     & current_velocity_in_s,HC_temperature,&
                     & HC_mass,nrand,total_steps,&
                     & fp_timesteps,Max_HC_Ion_Iter_To_TMax,ctemav,rsect,zsect,fp_return)


                ! jde - restore velocity to expected value in this routine
                current_velocity_in_s = current_velocity_in_s / ion_time_step

                ! Krieger IPP/07 - SUN compiler insists on 132 column limit
                if (debug_hc) & 
                     & write(6,'(a,i5,5g15.5)') 'HC EXITING FP :',fp_return,hc_temperature,current_s,current_velocity_in_s, &
                     &                                            total_steps



                !
                !                If a wall collision has occured then refine the rsect,zsect 
                !                values returned to actual wall locations and determine the ID
                !                and IS values related to that location. 
                !
                if (fp_return.eq.3) then 
                   call find_nearest_point_on_wall(rsect,zsect,id_out,is_out)
                endif

             else

                ! jdemod - add initial crossfield start position in FP to the FPERIPH call
                FP_Return = FPERIPH (Total_Steps, FP_TimeSteps, FP_Distance, FP_Loss_Time , &
                     &  Max_HC_Ion_Iter_To_TMax,NRand, FP_Diff_Rate,Seed,0.0)
             endif

             ! Update ion lifetimes for far periphery excursion
             Total_Steps  = Total_Steps + FP_TimeSteps
             Ion_Time_Step_Count = Ion_Time_Step_Count + FP_TimeSteps

             !write (0,*) "HC leaving FP:",Total_Steps,FP_TimeSteps,FP_Distance,FP_Loss_Time,FP_Return, Max_HC_Ion_Iter_To_TMax, FP_Diff_Rate

             If (Debug) Then
                Write (Output_Unit_Scratch,'(a,i5,1p,6g12.5)') 'FP HC:',FP_Return,Total_Steps, Max_HC_Ion_Iter_To_TMax, &
                     &  Ion_Time_Step,FP_Distance,FP_Loss_Time, FP_Diff_Rate
             End If

             If (FP_Return .eq. 1) Then
                ! SET CROSS TO ONE STEP INSIDE THE BOUNDARY.
                ! THIS IS THE POSITION IT WILL HAVE REACHED
                ! IN ORDER TO EXIT THE FPERIPH ROUTINE

                HC_FPExit (Cur_HC_Spec,Launch_Reg) =  HC_FPExit (Cur_HC_Spec,Launch_Reg) + Sput_Weight

                If ( Ion_Wall_Opt .eq. 1) Then
                   If (Current_Ring .eq.  Inner_Wall_Ring .or. Current_Ring .eq.  Outer_Wall_Ring) Then
                      Current_Cross = gkperps (Current_Cell,Current_Ring)
                   ElseIf (Current_Ring .eq.  Upper_Trap_Ring .or. Current_Ring .eq.  Lower_Trap_Ring) Then
                      Current_Cross = -gkperps (Current_Cell,Current_Ring)
                   End If
                ElseIf ( Ion_Wall_Opt .eq. 0 .or.  Ion_Wall_Opt .eq. 2) Then
                   If (Current_Ring .eq.  Inner_Wall_Ring .or. Current_Ring .eq.  Outer_Wall_Ring) Then
                      Temp_Cell = gikins (Current_Cell,Current_Ring)
                      Current_Ring = girins (Current_Cell,Current_Ring)

                      If ( Non_Orthogonal_Grid_Opt .eq. 1 .or.  Non_Orthogonal_Grid_Opt .eq. 3) Then
                         ! nonorth
                         If (Current_Ring .eq.  Num_Upper_Rings .and. Temp_Cell .ge.  Cell_TI) Then
                            Current_Theta = Current_Theta -  DThetaG
                         End If
                         ! nonorth
                      End If

                      Current_Cell = Temp_Cell
                      Current_Cross = -gdistout (Current_Cell,Current_Ring) + gkperps (Current_Cell,Current_Ring)

                   ElseIf (Current_Ring .eq.  Upper_Trap_Ring .or. Current_Ring .eq.  Lower_Trap_Ring) Then
                      Temp_Cell = gikouts (Current_Cell,Current_Ring)
                      Current_Ring = girouts (Current_Cell,Current_Ring)

                      If ( Non_Orthogonal_Grid_Opt .eq. 1 .or.  Non_Orthogonal_Grid_Opt .eq. 3) Then
                         ! nonorth
                         If (Current_Ring .eq.  Inner_SOL_Ring .and. Temp_Cell .gt.  Cell_TO) Then
                            Current_Theta = Current_Theta +  DThetaG
                         End If
                         ! nonorth
                      End If

                      Current_Cell = Temp_Cell
                      Current_Cross = gdistin (Current_Cell,Current_Ring) - gkperps (Current_Cell,Current_Ring)
                   End If

                   If ( Non_Orthogonal_Grid_Opt .eq. 1 .or.  Non_Orthogonal_Grid_Opt .eq. 3) Then
                      ! nonorth
                      If (gtagdv (Temp_Cell,Current_Ring) .eq. 0 .and. gtagdv (Last_Cell,Last_Ring) .eq. 0) Then
                         ! BOTH POINTS UNSHIFTED FROM ORTHOGONAL GRID
                         Current_Cell = Temp_Cell
                      Else
                         ! EITHER POINT SHIFTED FROM ORTHOGONAL GRID
                         If (Last_Ring .eq.  Inner_Wall_Ring) Then
                            Temp_Cell = giking (Current_Cell,Last_Ring)
                         ElseIf (Last_Ring .eq.  Upper_Trap_Ring) Then
                            Temp_Cell = gikoutg (Current_Cell,Last_Ring)
                         End If

                         ! Infinite loop bug - for temp_cell=1 and current_theta < thetag(temp_cell,current_ring)

                         Do
                            If (Current_Theta .lt. gthetag (Temp_Cell,Current_Ring)) Then
                               If (Temp_Cell .gt. 1) Then
                                  Temp_Theta = (gthetag (Temp_Cell,Current_Ring) - Current_Theta) / (gthetag (Temp_Cell,&
                                       &Current_Ring) - gthetag (Temp_Cell - 1,Current_Ring))
                                  If (Temp_Theta .gt. 0.5) Then
                                     Temp_Cell = Temp_Cell - 1
                                  Else
                                     Exit
                                  End If
                               else
                                  exit
                               End If
                            Else
                               If (Temp_Cell .lt. gnks (Current_Ring)) Then
                                  Temp_Theta = (Current_Theta - gthetag (Temp_Cell,Current_Ring)) / (gthetag (Temp_Cell+1,&
                                       &Current_Ring) - gthetag (Temp_Cell,Current_Ring))
                                  If (Temp_Theta .ge. 0.5) Then
                                     Temp_Cell = Temp_Cell + 1
                                  Else
                                     Exit
                                  End If
                               else
                                  exit
                               End If
                            End If
                         End Do
                         Current_Cell = Temp_Cell
                      End If

                      If (Last_Ring .eq.  Inner_Wall_Ring) Then
                         Current_Cross = -gdistout (Current_Cell,Current_Ring) + gkperps (Current_Cell,Current_Ring)
                      ElseIf (Last_Ring .eq.  Upper_Trap_Ring) Then
                         Current_Cross = gdistin (Current_Cell,Current_Ring) - gkperps (Current_Cell,Current_Ring)
                      End If
                      ! nonorth
                   End If
                End If
             ElseIf (FP_Return .eq. 2) Then
                ! CIST HAS EXCEEDED CSTMAX - GOTO 780 TO FINISH PROCESSING.
                HC_FPTTotal (Cur_HC_Spec,Launch_Reg) =  HC_FPTTotal (Cur_HC_Spec,Launch_Reg) + Sput_Weight

                ! jdemod - variable needs second index to be correct
                !HC_Num_Reach_Max_Iter (Cur_HC_Spec) =  HC_Num_Reach_Max_Iter (Cur_HC_Spec) + Sput_Weight
                HC_Num_Reach_Max_Iter (Cur_HC_Spec,launch_reg) =  HC_Num_Reach_Max_Iter (Cur_HC_Spec,launch_reg) + Sput_Weight

                ! jdemod - variable needs second index to be correct
                ! HC_Tot_Temp_Reach_Max_Iter (Cur_HC_Spec) = HC_Tot_Temp_Reach_Max_Iter(Cur_HC_Spec) + HC_Temperature * Sput_Weight
                HC_Tot_Temp_Reach_Max_Iter (Cur_HC_Spec,launch_reg) =  HC_Tot_Temp_Reach_Max_Iter (Cur_HC_Spec,launch_reg) &
                     & + HC_Temperature * Sput_Weight
                HC_Tot_Exist_Max_Iter (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Exist_Max_Iter (Cur_HC_Spec,Launch_Reg) + Sput_Weight
                IFate = 29 ! HC Ion Time = TMax.
                Return

             ElseIf (FP_Return .eq. 3) Then
                ! TREAT AS NORMAL WALL COLLISION
                If ( Far_Periphery_Recycle_Opt .eq. 0) Then

                   ! Normal collision - no recycling
                   HC_Num_Absorbed_TargWall(Cur_HC_Spec,Launch_Reg)=HC_Num_Absorbed_TargWall(Cur_HC_Spec,Launch_Reg)+Sput_Weight
                   HC_Min_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) = MIN ( HC_Min_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg), &
                        &Real(Total_Steps))
                   HC_Max_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) = MAX ( HC_Max_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg), &
                        &Real(Total_Steps))
                   HC_Tot_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) + &
                        &Total_Steps * Sput_Weight

                   ! Recording addition to total elapsed time in state Cur_HC_Spec.
                   HC_Tot_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) + &
                        & HC_Temperature * Sput_Weight
                   HC_Tot_Velocity_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Velocity_At_Absorption (Cur_HC_Spec,Launch_Reg)&
                        & + Current_Velocity_In_S * Sput_Weight
                   HC_Tot_ABS_Vel_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_ABS_Vel_At_Absorption (Cur_HC_Spec,Launch_Reg) +&
                        & ABS (Current_Velocity_In_S) * Sput_Weight
                   ! jdemod - this line was partially changed - variable needs second index to be correct
                   !HC_Tot_Elec_Temp_At_Absorption (Cur_HC_Spec) =  HC_Tot_Elec_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) + &
                   HC_Tot_Elec_Temp_At_Absorption (Cur_HC_Spec,launch_reg) = &
                        & HC_Tot_Elec_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) + &
                        & gktebs (Current_Cell, Current_Ring) * Sput_Weight

                   IM_Index = MIN (INT (HC_Temperature / (0.2 *  Init_Electron_Temperature)) + 1, 10)
                   ! jdemod - array use does not match declaration - add second index
                   !HC_Ctexs (IM_index)  =  HC_Ctexs (IM_index) + HC_Temperature * Sput_Weight
                   HC_Ctexs (IM_Index,Launch_reg)  =  HC_Ctexs (IM_Index,Launch_reg) + HC_Temperature * Sput_Weight

                   HC_RWall (Cur_HC_Spec,Launch_Reg) =  HC_RWall (Cur_HC_Spec,Launch_Reg)+ Sput_Weight

                   HC_Walls (Current_Cell,Current_Ring,Cur_HC_Spec) = &
                        &  HC_Walls (Current_Cell,Current_Ring,Cur_HC_Spec) + Sput_Weight

                   ! Add ion weight to wall element closest to grid departure.
                   Call HC_Update_WallDep(Current_Cell,Current_Ring,Cur_HC_Spec,0,id_out, &
                        & Launch_Wall_Index,NeutType,Sput_Weight,Launch_Reg)

                   IFate = 20 ! HC Ion Reached Vessel Wall (DIV IFATE = 1).
                   Return

                ElseIf ( Far_Periphery_Recycle_Opt .eq. 1) Then

                   ! Recycle particle from edge of nearest plate

                   ! For the outer wall the particle will
                   ! be launched from either ID=1 or ID=NDS
                   ! For the trap wall it will be launched
                   ! from ID = NDSIN or ID = NDSIN+1.

                   ! Unless DDS(ID) = 0 in which case it shifts
                   ! accross to the first target segment with non-zero size.

                   ! Note: To have reached this point IR must be equal to IRWALL or IRTRAP.

                   ! Record statistics as if for normal wall collision
                   !Write (Output_Unit_Scratch,*) 'FP HC RELAUNCH BW:',Current_Cell,Current_Ring,Current_R,Current_Z,Cur_HC_Spec

                   ! REFLECT ---- IONS may need to be reflected here
                   HC_Num_Absorbed_TargWall(Cur_HC_Spec,Launch_Reg)=HC_Num_Absorbed_TargWall(Cur_HC_Spec,Launch_Reg)+Sput_Weight
                   HC_Min_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) = MIN ( HC_Min_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg), &
                        &Real(Total_Steps))
                   HC_Max_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) = MAX ( HC_Max_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg), &
                        &Real(Total_Steps))
                   HC_Tot_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) + &
                        & Eq_Total_Ion_Time_Steps * Sput_Weight

                   ! Recording addition to total elapsed time in state Cur_HC_Spec.
                   HC_Tot_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) + &
                        &HC_Temperature * Sput_Weight
                   HC_Tot_Velocity_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Velocity_At_Absorption (Cur_HC_Spec,Launch_Reg)&
                        & + Current_Velocity_In_S * Sput_Weight
                   HC_Tot_ABS_Vel_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_ABS_Vel_At_Absorption (Cur_HC_Spec,Launch_Reg) +&
                        & ABS (Current_Velocity_In_S) * Sput_Weight
                   HC_Tot_Elec_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Elec_Temp_At_Absorption (Cur_HC_Spec,&
                        &Launch_Reg) + gktebs (Current_Cell, Current_Ring) * Sput_Weight

                   IM_Index = MIN (INT (HC_Temperature / (0.2 *  Init_Electron_Temperature)) + 1, 10)
                   ! jdemod - array use does not match declaration - add second index
                   !HC_Ctexs (IM_index)  =  HC_Ctexs (IM_index) + HC_Temperature * Sput_Weight
                   HC_Ctexs (IM_Index,Launch_reg)  =  HC_Ctexs (IM_Index,Launch_reg) + HC_Temperature * Sput_Weight

                   HC_RWall (Cur_HC_Spec,Launch_Reg) =  HC_RWall (Cur_HC_Spec,Launch_Reg) + Sput_Weight

                   HC_Walls (Current_Cell,Current_Ring,Cur_HC_Spec) = &
                        &  HC_Walls (Current_Cell,Current_Ring,Cur_HC_Spec) + Sput_Weight

                   ! Add ion weight to wall element closest to grid departure.
                   Call HC_Update_WallDep(Current_Cell,Current_Ring,Cur_HC_Spec,0,id_out,&
                        & Launch_Wall_Index,NeutType,Sput_Weight,Launch_Reg)

                   Reflect_Ion = .False.
                   Sputter_Ion = .False.
                   Reflection_Probability = 0.0

                   ! Find target segment for re-launch
                   If (Current_Cell .gt. gnks (Current_Ring) / 2) Then
                      Current_Cell = gnks (Current_Ring)
                      Target_Index = Verify_ID (Current_Cell,Current_Ring,1)
                   Else
                      Current_Cell = 1
                      Target_Index = Verify_ID (Current_Cell,Current_Ring,2)
                   End If

                   ! Postion on target/initial position options
                   If ( Neutral_Init_Pos_Opt .eq. 0) Then
                      Current_R = grp (Target_Index)
                      Current_Z = gzp (Target_Index)
                   ElseIf ( Neutral_Init_Pos_Opt .eq. 1) Then
                      Call Position_On_Target (Current_R,Current_Z,Current_Cross,Target_Index)
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
                      Call HC_Reflection_Coefficient (Target_Index, Cur_HC_Spec, HC_Temperature, HC_Angle_To_Normal, &
                           &Reflection_Probability)
                      If (hc_ion_reflection_option .eq. 1 .and. Random_Value .le. Reflection_Probability) Then
                         Reflect_Ion = .True.
                      End If

                      ! If reflection does not occur, check to see if HC sputters.
                      Call HC_Sputtering_Coefficient (Target_Index, Cur_HC_Spec, HC_Temperature, HC_Angle_To_Normal, &
                           &Sputtering_Probability)
                      If (.not. Reflect_Ion .and. hc_sputtering_option .eq. 1 .and. Random_Value .le. (Reflection_Probability + &
                           &Sputtering_Probability)) Then
                         Sputter_Ion = .True.
                      End If
                   End If

                   ! Assign sputtering and reflection angle from vessel segment, each +/-90 degrees isotropic from normal.
                   Segment_Normal_Angle = Find_Vessel_Segment_Normal (Current_Cell,Current_Ring)
                   NRand = NRand + 1
                   Call Surand2 (Seed, 1, Random_Value_1)
                   NRand = NRand + 1
                   Call Surand2 (Seed, 1, Random_Value_2)
                   Angle_Offset = SIGN (Random_Value_1 *  Pi_Value / 2.0, 0.5 - Random_Value_2)
                   Reflection_Angle = Segment_Normal_Angle + Angle_Offset
                   Sputtering_Angle = Segment_Normal_Angle + Angle_Offset

                   ! Test for reflection and sputtering.
                   If (Reflect_Ion) Then
                      ! Particle is reflected.  Find new species.
                      !write (0,*) "Reflecting FP31 ion",IProd,LPRod,NProd,Cur_HC_Spec,Current_R,Current_Z
                      ! Particle is reflected back into plasma.  Find new species and update all statistics.
                      Call HC_Ion_Impact_ReLaunch ("Reflect",Cur_HC_Spec,Last_HC_Species,H_Isotope_Composition,&
                           & Last_H_Isotope_Composition,Launch_Reg, &
                           & Sput_Weight,Neutral_Time_Step_Count,Ion_Time_Step_Count,Eq_Total_Ion_Time_Steps,Reflection_Counter,&
                           & NeutType,Target_Index,Current_Theta, &
                           & Current_R,Current_Z,HC_Temperature,Current_Velocity,Current_Velocity_In_S,Current_Angle,Seed,NRand,&
                           & Current_Cell,Current_Ring, &
                           & Current_Velocity_In_R,Current_Velocity_In_Z,Current_S,Current_Cross,hc_v)

                      !Write (Output_Unit_Scratch,*) 'Ion reflection in FP(3) at target:',Current_Cell,Current_Ring,Target_Index,Current_R,Current_Z,gwallindex(Target_Index),ryield, &
                      !&  Target_Material,Dep_Energy,Sput_Weight,Segment_Normal_Angle,Current_Theta,HC_Temperature,Current_Velocity,Current_Angle

                      ! If output is requested, open files.
                      If (hc_evolve_print_option .eq. 1) Then
                         Write (Output_Unit_Evolve,9510) "Ion reflection in FP(3) at target index:",Target_Index,"R:",Current_R,"Z:&
                              &",Current_R,"Last HC:",Last_HC_Species,"New HC:",Cur_HC_Spec
9510                     Format (3X,A,1X,I3,1X,A,1X,F8.3,1X,A,1X,F8.3,1X,A,1X,I2,1X,A,1X,I2)
                      End If
                      If (hc_coord_print_option .eq. 1) Then
                         Write (Output_Unit_Location,9511) "Ion reflection in FP(3) at target index:",Target_Index,"R:",Current_R,&
                              & "Z:",Current_R,"Last HC:",Last_HC_Species,"New HC:",Cur_HC_Spec
9511                     Format (3X,A,1X,I3,1X,A,1X,F8.3,1X,A,1X,F8.3,1X,A,1X,I2,1X,A,1X,I2)
                      End If

                      ! jdemod - this calculation of the R,Z velocity components is performed 
                      !          in hc_ion_target_reflect_velocity called from hc_ion_impact_relaunch 
                      !          for reflected particles and does not need to be repeated here. 
                      ! Calculate R,Z components of velocity, probability of ionization, etc.
                      !Current_Velocity_In_R = Current_Velocity * COS (Current_Angle) *  Neutral_Time_Step
                      !Current_Velocity_In_Z = Current_Velocity * SIN (Current_Angle) *  Neutral_Time_Step
                      !
                      ! New position, energy/velocity, and direction are loaded.  Return to following code.

                      Return
                   Else
                      ! ION is NOT reflected.  Record deposition and check for self-sputtering.
                      ! Record all the regular loss statistics.
                      HC_Num_Absorbed_TargWall(Cur_HC_Spec,Launch_Reg)=HC_Num_Absorbed_TargWall(Cur_HC_Spec,Launch_Reg)+Sput_Weight
                      HC_Min_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) = MIN ( HC_Min_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg), &
                           & REAL (Eq_Total_Ion_Time_Steps))
                      HC_Max_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) = MAX ( HC_Max_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg), &
                           & REAL (Eq_Total_Ion_Time_Steps))
                      HC_Tot_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) + &
                           & Eq_Total_Ion_Time_Steps * Sput_Weight

                      ! Recording addition to total elapsed time in state Cur_HC_Spec.
                      HC_Tot_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) + &
                           & HC_Temperature * Sput_Weight
                      HC_Tot_Velocity_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Velocity_At_Absorption (Cur_HC_Spec,&
                           & Launch_Reg) + Current_Velocity_In_S * Sput_Weight
                      HC_Tot_ABS_Vel_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_ABS_Vel_At_Absorption (Cur_HC_Spec,&
                           & Launch_Reg) + ABS (Current_Velocity_In_S) * Sput_Weight
                      HC_Tot_Elec_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Elec_Temp_At_Absorption (Cur_HC_Spec,&
                           & Launch_Reg) + gktebs (Current_Cell, Current_Ring) * Sput_Weight

                      IM_Index = MIN (INT (HC_Temperature / (0.2 *  Init_Electron_Temperature)) + 1, 10)
                      ! jdemod - array use does not match declaration - add second index
                      !HC_Ctexs (IM_index)  =  HC_Ctexs (IM_index) + HC_Temperature * Sput_Weight
                      HC_Ctexs (IM_Index,Launch_reg)  =  HC_Ctexs (IM_Index,Launch_reg) + HC_Temperature * Sput_Weight
                      HC_Num_Absorbed_Act_Target (Cur_HC_Spec,Launch_Reg) =  HC_Num_Absorbed_Act_Target (Cur_HC_Spec,Launch_Reg) + &
                           & Sput_Weight
                      HC_RDep (Cur_HC_Spec,Launch_Reg) =  HC_RDep (Cur_HC_Spec,Launch_Reg) + Sput_Weight

                      ! Do not record the deposit/erosion statistics of this as a standard relaunched particle.
                      Dep_Energy = 3.0 * Get_HC_Charge (Cur_HC_Spec) * gktebs (Current_Cell,Current_Ring) + 5.22E-9 * Find_HC_Mass &
                           & (Cur_HC_Spec,H_Isotope_Composition) * &
                           & Current_Velocity_In_S /  Ion_Time_Step * Current_Velocity_In_S /  Ion_Time_Step + 2.0 * HC_Temperature

                      If (Yield_Calc .eq. 1 .and. Sputter_Ion) Then

                         ! Particle is sputtered back into plasma.  Find new species and update all statistics.
                         Call HC_Ion_Impact_ReLaunch ("Sputter",Cur_HC_Spec,Last_HC_Species,H_Isotope_Composition,&
                              &Last_H_Isotope_Composition,Launch_Reg, &
                              & Sput_Weight,Neutral_Time_Step_Count,Ion_Time_Step_Count,Eq_Total_Ion_Time_Steps,Reflection_Counter,&
                              &NeutType,Target_Index,Current_Theta, &
                              & Current_R,Current_Z,HC_Temperature,Current_Velocity,Current_Velocity_In_S,Current_Angle,Seed,NRand,&
                              &Current_Cell,Current_Ring, &
                              & Current_Velocity_In_R,Current_Velocity_In_Z,Current_S,Current_Cross,hc_v)

                         ! New position, energy/velocity, and direction are loaded.  Return to following code.
                         !Write (Output_Unit_Scratch,*) 'HC SPUTTERED from ion impact FP=3:',Current_Cell,Current_Ring,Target_Index,Current_R,Current_Z,gikds (Target_Index),girds (Target_Index),ryield,gkmfss (Target_Index), Target_Material,Dep_Energy,Sput_Weight

                         ! If output is requested, open files.
                         If (hc_evolve_print_option .eq. 1) Then
                            Write (Output_Unit_Evolve,9520) "Ion sputter in FP(3) at target index:",Target_Index,"R:",Current_R,"Z:&
                                 &",Current_R,"Last HC:",Last_HC_Species,"New HC:",Cur_HC_Spec
9520                        Format (3X,A,1X,I3,1X,A,1X,F8.3,1X,A,1X,F8.3,1X,A,1X,I2,1X,A,1X,I2)
                         End If
                         If (hc_coord_print_option .eq. 1) Then
                            Write (Output_Unit_Location,9521) "Ion sputter in FP(3) at target index:",Target_Index,"R:",Current_R,&
                                 &"Z:",Current_R,"Last HC:",Last_HC_Species,"New HC:",Cur_HC_Spec
9521                        Format (3X,A,1X,I3,1X,A,1X,F8.3,1X,A,1X,F8.3,1X,A,1X,I2,1X,A,1X,I2)
                         End If

                         ! jdemod - this calculation of the R,Z velocity components is performed 
                         !          in hc_ion_target_sputter_velocity called from hc_ion_impact_relaunch 
                         !          for reflected particles and does not need to be repeated here. 
                         ! Calculate R,Z components of velocity, probability of ionization, etc.
                         !Current_Velocity_In_R = Current_Velocity * COS (Current_Angle) *  Neutral_Time_Step
                         !Current_Velocity_In_Z = Current_Velocity * SIN (Current_Angle) *  Neutral_Time_Step

                         Return

                      ElseIf (Yield_Calc .eq. 2 .and. Sputter_Ion) Then
                         ! DIVIMP yield determination.
                         If (gkmfss(Target_Index) .ge. 0.0) Then
                            RYield = YIELD (6,  Target_Material,Dep_Energy,gktebs (Current_Cell,Current_Ring),gktibs (Current_Cell,&
                                 &Current_Ring)) * gkmfss (Target_Index)
                         ElseIf (gkmfss (Target_Index) .lt. 0.0 .and. gkmfss (Target_Index) .ge. -50.0) then
                            RYield = ABS (gkmfss (Target_Index))
                         ElseIf (gkmfss(Target_Index) .le. -99.0) Then
                            RYield = YIELD (6, Target_Material,Dep_Energy,gktebs (Current_Cell,Current_Ring),gktibs (Current_Cell,&
                                 &Current_Ring))
                         End If

                         SputNew = Sput_Weight * RYield
                         HC_YldTot (Cur_HC_Spec,Launch_Reg) =  HC_YldTot (Cur_HC_Spec,Launch_Reg) + SputNew
                         HC_YldMax (Cur_HC_Spec,Launch_Reg) = MAX ( HC_YldMax (Cur_HC_Spec,Launch_Reg), SputNew)

                         !Write (Output_Unit_Scratch,*) 'HC SPUTTERED from ion impact:',Current_Cell,Current_Ring,Target_Index,Current_R,Current_Z,gikds (Target_Index),girds (Target_Index),ryield, &
                         !& gkmfss (Target_Index), Target_Material,Dep_Energy,Sput_Weight,SputNew, HC_YldTot

                         If (SputNew .gt.  Self_Sputter_Threshold) Then
                            HC_YthTot (Cur_HC_Spec,Launch_Reg) =  HC_YthTot (Cur_HC_Spec,Launch_Reg) + SputNew
                            NProd  = NProd + 1
                            Call psnews (NPROD,SputNew)
                            If (hc_launch_angle_velocity .eq. 1 .or. hc_launch_angle_velocity .eq. 4 .or. hc_launch_angle_velocity &
                                 &.eq. 5 .or.  Sputter_Opt .eq. 4) Then
                               EMax =  EMax_Factor * Dep_Energy
                               Max_Velocity_Randoms = 1.0 / (1.0 +  Target_Binding_Energy / EMax)**2.0
                            Else
                               Max_Velocity_Randoms = 1.0
                            End If

                            Call pxprods (NProd,Current_R)
                            Call pyprods (NProd,Current_Z)

                            ! For segments with a fixed sputtering yield - allow for
                            ! the energy of the sputtered particle to be set to a
                            ! specific value.

                            If ( Self_Sputter_Opt .eq. 2 .and. gkmfss (Target_Index) .lt. 0.0) Then
                               Call peprods (NProd,hc_sput_energy_ion_preset)
                            Else
                               Call peprods (NProd,0.0)
                            End If

                            Call pidprods (NProd,Target_Index)
                            HC_Launchdat (NProd,2) = 1.0

                         End If

                         If (Debug) Then
                            Write (Output_Unit_Scratch,*) 'FP HC RELAUNCH:',Current_Cell,Current_Ring,Target_Index,Current_R,&
                                 &Current_Z,gikds (Target_Index),girds (Target_Index),SputNew,NProd,Dep_Energy
                         End If

                         ! Record collision with target.
                         HC_Num_Striking_Target(Cur_HC_Spec,Launch_Reg)=HC_Num_Striking_Target(Cur_HC_Spec,Launch_Reg)+Sput_Weight

                         IFate = 22 ! HC Ion Struck Target.
                         Return
                      End If
                   End If
                End If

             ElseIf (FP_Return .eq. 4) Then
                ! FP TARGET IMPACT - RECORD AND GOTO NEXT ION
                ! TREAT AS NORMAL TARGET LOSS FOR STATISTICS BUT
                ! NOT SELF-SPUTTERING

                If ( Far_Periphery_Recycle_Opt .eq. 0) Then
                   HC_Num_Absorbed_TargWall(Cur_HC_Spec,Launch_Reg)=HC_Num_Absorbed_TargWall(Cur_HC_Spec,Launch_Reg)+Sput_Weight
                   HC_Min_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) = MIN ( HC_Min_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg), &
                        & Real(Total_Steps))
                   HC_Max_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) = MAX ( HC_Max_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg), &
                        & Real(Total_Steps))
                   HC_Tot_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) + &
                        & Eq_Total_Ion_Time_Steps * Sput_Weight

                   ! Recording addition to total elapsed time in state Cur_HC_Spec.
                   HC_Tot_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) + &
                        & HC_Temperature * Sput_Weight
                   HC_Tot_Velocity_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Velocity_At_Absorption (Cur_HC_Spec,Launch_Reg)&
                        & + Current_Velocity_In_S * Sput_Weight
                   HC_Tot_ABS_Vel_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_ABS_Vel_At_Absorption (Cur_HC_Spec,Launch_Reg) +&
                        & ABS (Current_Velocity_In_S) * Sput_Weight
                   HC_Tot_Elec_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Elec_Temp_At_Absorption (Cur_HC_Spec,&
                        & Launch_Reg) + gktebs (Current_Cell, Current_Ring) * Sput_Weight

                   IM_Index = MIN (INT (HC_Temperature / (0.2 *  Init_Electron_Temperature)) + 1, 10)
                   ! jdemod - array use does not match declaration - add second index
                   !HC_Ctexs (IM_index)  =  HC_Ctexs (IM_index) + HC_Temperature * Sput_Weight
                   HC_Ctexs (IM_Index,Launch_reg)  =  HC_Ctexs (IM_Index,Launch_reg) + HC_Temperature * Sput_Weight

                   HC_FPTarg (Cur_HC_Spec,Launch_Reg) =  HC_FPTarg (Cur_HC_Spec,Launch_Reg) + Sput_Weight
                   HC_FPTart (Cur_HC_Spec,Launch_Reg) =  HC_FPTart (Cur_HC_Spec,Launch_Reg) + Sput_Weight
                   HC_RFPTarg (Cur_HC_Spec,Launch_Reg) =  HC_RFPTarg (Cur_HC_Spec,Launch_Reg) + Sput_Weight

                   Reflect_Ion = .False.
                   Sputter_Ion = .False.
                   Reflection_Probability = 0.0

                   ! Find closest target segment for loss
                   If (Current_Cell .gt. gnks (Current_Ring) / 2) Then
                      Current_Cell = gnks (Current_Ring)
                      Target_Index = Verify_ID (Current_Cell,Current_Ring,1)
                   Else
                      Current_Cell = 1
                      Target_Index = Verify_ID (Current_Cell,Current_Ring,2)
                   End If

                   ! Add ion weight to wall element  closest to grid departure.
                   !Write (Output_Unit_Scratch,*) 'HC FPTARG:',Current_Cell,Current_Ring,Target_Index,gwallindex (Target_Index)

                   Call HC_Update_WallDep (Current_Cell,Current_Ring,Cur_HC_Spec,Target_Index,0, &
                        & Launch_Wall_Index,NeutType,Sput_Weight,Launch_Reg)

                   IFate = 25 ! HC Ion Hit Far-Periphery Target (DIV IFATE = 9).
                   Return

                ElseIf ( Far_Periphery_Recycle_Opt .eq. 1) Then
                   ! Recycle particle from edge of nearest plate.

                   ! Particles trigger a self-sputtering event upon recycling.

                   ! For the outer fp target the particle will
                   ! be launched from either ID=1 or ID=NDS
                   ! For the trap fp target it will be launched
                   ! from ID = NDSIN or ID = NDSIN+1.

                   ! Unless DDS(ID) = 0 in which case it shifts
                   ! accross to the first target segment with non-zero size.

                   ! Note: To have reached this point IR must be equal to IRWALL or IRTRAP.
                   ! Record statistics as if for normal FP-TARGET collision.

                   If (Debug) Then
                      Write (Output_Unit_Scratch,*) 'FP HC RELAUNCH BT:',Current_Cell,Current_Ring,Current_R,Current_Z,Cur_HC_Spec
                   End If

                   ! REFLECT ---- IONS may need to be reflected here
                   HC_Num_Absorbed_TargWall(Cur_HC_Spec,Launch_Reg)=HC_Num_Absorbed_TargWall(Cur_HC_Spec,Launch_Reg)+Sput_Weight
                   HC_Min_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) = MIN ( HC_Min_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg), &
                        & Real(Total_Steps))
                   HC_Max_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) = MAX ( HC_Max_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg), &
                        & Real(Total_Steps))
                   HC_Tot_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) + &
                        & Eq_Total_Ion_Time_Steps * Sput_Weight

                   ! Recording addition to total elapsed time in state Cur_HC_Spec.
                   HC_Tot_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) + &
                        & HC_Temperature * Sput_Weight
                   HC_Tot_Velocity_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Velocity_At_Absorption (Cur_HC_Spec,Launch_Reg)&
                        & + Current_Velocity_In_S * Sput_Weight
                   HC_Tot_ABS_Vel_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_ABS_Vel_At_Absorption (Cur_HC_Spec,Launch_Reg) +&
                        & ABS (Current_Velocity_In_S) * Sput_Weight
                   HC_Tot_Elec_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Elec_Temp_At_Absorption (Cur_HC_Spec,&
                        &Launch_Reg) + gktebs (Current_Cell, Current_Ring) * Sput_Weight

                   IM_Index = MIN (INT (HC_Temperature / (0.2 *  Init_Electron_Temperature)) + 1, 10)
                   ! jdemod - array use does not match declaration - add second index
                   !HC_Ctexs (IM_index)  =  HC_Ctexs (IM_index) + HC_Temperature * Sput_Weight
                   HC_Ctexs (IM_Index,launch_reg)  =  HC_Ctexs (IM_Index,launch_reg) + HC_Temperature * Sput_Weight

                   HC_FPTarg (Cur_HC_Spec,Launch_Reg) =  HC_FPTarg (Cur_HC_Spec,Launch_Reg) + Sput_Weight
                   HC_FPTart (Cur_HC_Spec,Launch_Reg) =  HC_FPTart (Cur_HC_Spec,Launch_Reg) + Sput_Weight
                   HC_RFPTarg (Cur_HC_Spec,Launch_Reg) =  HC_RFPTarg  (Cur_HC_Spec,Launch_Reg)+ Sput_Weight

                   Reflect_Ion = .False.
                   Sputter_Ion = .False.
                   Reflection_Probability = 0.0

                   ! Find target segment for re-launch
                   If (Current_Cell .gt. gnks (Current_Ring) / 2) Then
                      Current_Cell = gnks (Current_Ring)
                      Target_Index = Verify_ID (Current_Cell,Current_Ring,1)
                   Else
                      Current_Cell = 1
                      Target_Index = Verify_ID (Current_Cell,Current_Ring,2)
                   End If

                   ! Position on target/initial position options.
                   If ( Neutral_Init_Pos_Opt .eq. 0) Then
                      Current_R = grp (Target_Index)
                      Current_Z = gzp (Target_Index)
                   ElseIf ( Neutral_Init_Pos_Opt .eq. 1) Then
                      Call Position_On_Target (Current_R,Current_Z,Current_Cross,Target_Index)
                   End If

                   ! Find angle of magnetic field.
                   Current_Angle = Find_Ion_Angle (Current_Cell,Current_Ring)

                   ! Check if reflection is allowed.

                   ! jdemod - comment - the following code gets away with using one random number for
                   !                    either sputtering or reflection by summing the reflection and sputtering 
                   !                    probabilities in the second IF statement. 
                   If (hc_ion_reflection_option .eq. 1 .or. hc_sputtering_option .eq. 1) Then
                      ! Select random number.
                      NRand = NRand + 1
                      Call Surand2 (Seed, 1, Random_Value)

                      ! Find angle to normal at impacting vessel segment.
                      HC_Angle_To_Normal = Find_Angle_To_Normal (Current_Cell, Current_Ring, Current_Angle)

                      ! Check to see how likely reflection is from HC_Stick module.
                      Call HC_Reflection_Coefficient (Target_Index, Cur_HC_Spec, HC_Temperature, HC_Angle_To_Normal, &
                           & Reflection_Probability)
                      If (hc_ion_reflection_option .eq. 1 .and. Random_Value .le. Reflection_Probability) Then
                         Reflect_Ion = .True.
                      End If

                      ! If reflection does not occur, check to see if HC sputters.
                      Call HC_Sputtering_Coefficient (Target_Index, Cur_HC_Spec, HC_Temperature, HC_Angle_To_Normal, &
                           & Sputtering_Probability)
                      ! slmod begin
                      If ( .not.Reflect_Ion .and. hc_sputtering_option .eq. 1 .and. Random_Value .le. (&
                           & Reflection_Probability + Sputtering_Probability)) Then
                         !
                         !                      If (Reflect_Ion .eq. .False. .and. hc_sputtering_option .eq. 1 .and. Random_Value .le. (&
                         ! slmod end

                         Sputter_Ion = .True.
                      End If
                   End If

                   ! Test for reflection and sputtering.
                   If (Reflect_Ion) Then
                      ! Particle is reflected.  Find new species.
                      !write (0,*) "Reflecting FP41 ion",IProd,LPRod,NProd,Cur_HC_Spec,Current_R,Current_Z

                      ! Particle is reflected back into plasma.  Find new species and update all statistics.
                      Call HC_Ion_Impact_ReLaunch ("Reflect",Cur_HC_Spec,Last_HC_Species,H_Isotope_Composition,&
                           & Last_H_Isotope_Composition,Launch_Reg, &
                           & Sput_Weight,Neutral_Time_Step_Count,Ion_Time_Step_Count,Eq_Total_Ion_Time_Steps,Reflection_Counter,&
                           & NeutType,Target_Index,Current_Theta, &
                           & Current_R,Current_Z,HC_Temperature,Current_Velocity,Current_Velocity_In_S,Current_Angle,Seed,NRand,&
                           & Current_Cell,Current_Ring, &
                           & Current_Velocity_In_R,Current_Velocity_In_Z,Current_S,Current_Cross,hc_v)

                      !Write (Output_Unit_Scratch,*) 'Ion reflection in FP(4) at target:',Current_Cell,Current_Ring,Target_Index,Current_R,Current_Z,gwallindex(Target_Index),ryield, &
                      !&  Target_Material,Dep_Energy,Sput_Weight,Segment_Normal_Angle,Current_Theta,HC_Temperature,Current_Velocity,Current_Angle

                      ! If output is requested, open files.
                      If (hc_evolve_print_option .eq. 1) Then
                         Write (Output_Unit_Evolve,9512) "Ion reflection in FP(4) at target index:",Target_Index,"R:",Current_R,"Z:&
                              &",Current_R,"Last HC:",Last_HC_Species,"New HC:",Cur_HC_Spec
9512                     Format (3X,A,1X,I3,1X,A,1X,F8.3,1X,A,1X,F8.3,1X,A,1X,I2,1X,A,1X,I2)
                      End If
                      If (hc_coord_print_option .eq. 1) Then
                         Write (Output_Unit_Location,9513) "Ion reflection in FP(4) at target index:",Target_Index,"R:",Current_R,&
                              &"Z:",Current_R,"Last HC:",Last_HC_Species,"New HC:",Cur_HC_Spec
9513                     Format (3X,A,1X,I3,1X,A,1X,F8.3,1X,A,1X,F8.3,1X,A,1X,I2,1X,A,1X,I2)
                      End If

                      ! jdemod - this calculation of the R,Z velocity components is performed 
                      !          in hc_ion_target_reflect_velocity called from hc_ion_impact_relaunch 
                      !          for reflected particles and does not need to be repeated here. 
                      ! Calculate R,Z components of velocity, probability of ionization, etc.
                      !Current_Velocity_In_R = Current_Velocity * COS (Current_Angle) *  Neutral_Time_Step
                      !Current_Velocity_In_Z = Current_Velocity * SIN (Current_Angle) *  Neutral_Time_Step

                      ! Continue following the newly reflected particle.
                      Return
                   Else
                      ! ION is NOT reflected.  Record deposition and check for self-sputtering.
                      ! Record all the regular loss statistics.
                      HC_Num_Absorbed_TargWall(Cur_HC_Spec,Launch_Reg)=HC_Num_Absorbed_TargWall(Cur_HC_Spec,Launch_Reg)+Sput_Weight
                      HC_Min_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) = MIN ( HC_Min_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg), &
                           & REAL (Eq_Total_Ion_Time_Steps))
                      HC_Max_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) = MAX ( HC_Max_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg), &
                           & REAL (Eq_Total_Ion_Time_Steps))
                      HC_Tot_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Teq_To_Absorption (Cur_HC_Spec,Launch_Reg) + &
                           & Eq_Total_Ion_Time_Steps * Sput_Weight

                      ! Recording addition to total elapsed time in state Cur_HC_Spec.
                      HC_Tot_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) + &
                           & HC_Temperature * Sput_Weight
                      HC_Tot_Velocity_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Velocity_At_Absorption (Cur_HC_Spec,&
                           & Launch_Reg) + Current_Velocity_In_S * Sput_Weight
                      HC_Tot_ABS_Vel_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_ABS_Vel_At_Absorption (Cur_HC_Spec,&
                           & Launch_Reg) + ABS (Current_Velocity_In_S) * Sput_Weight
                      HC_Tot_Elec_Temp_At_Absorption (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Elec_Temp_At_Absorption (Cur_HC_Spec,&
                           & Launch_Reg) + gktebs (Current_Cell, Current_Ring) * Sput_Weight

                      IM_Index = MIN (INT (HC_Temperature / (0.2 *  Init_Electron_Temperature)) + 1, 10)
                      ! jdemod - array use does not match declaration - add second index
                      !HC_Ctexs (IM_index)  =  HC_Ctexs (IM_index) + HC_Temperature * Sput_Weight
                      HC_Ctexs(IM_Index,launch_reg)  =  HC_Ctexs(IM_Index,launch_reg) + HC_Temperature * Sput_Weight
                      HC_Num_Absorbed_Act_Target (Cur_HC_Spec,Launch_Reg) =  HC_Num_Absorbed_Act_Target (Cur_HC_Spec,Launch_Reg) + &
                           & Sput_Weight
                      HC_RDep (Cur_HC_Spec,Launch_Reg) =  HC_RDep (Cur_HC_Spec,Launch_Reg) + Sput_Weight

                      ! Add ion weight to wall element closest to grid departure.
                      Call HC_Update_WallDep (Current_Cell,Current_Ring,Cur_HC_Spec,Target_Index,0,&
                           & Launch_Wall_Index,NeutType,Sput_Weight,Launch_Reg)

                      Dep_Energy = 3.0 * Get_HC_Charge (Cur_HC_Spec) * gktebs (Current_Cell,Current_Ring) + 5.22E-9 * Find_HC_Mass &
                           &(Cur_HC_Spec,H_Isotope_Composition) * &
                           & Current_Velocity_In_S /  Ion_Time_Step * Current_Velocity_In_S /  Ion_Time_Step + 2.0 * HC_Temperature

                      If (Yield_Calc .eq. 1 .and. Sputter_Ion) Then

                         ! Particle is sputtered back into plasma.  Find new species and update all statistics.
                         Call HC_Ion_Impact_ReLaunch ("Sputter",Cur_HC_Spec,Last_HC_Species,H_Isotope_Composition,&
                              & Last_H_Isotope_Composition,Launch_Reg, &
                              & Sput_Weight,Neutral_Time_Step_Count,Ion_Time_Step_Count,Eq_Total_Ion_Time_Steps,Reflection_Counter,&
                              & NeutType,Target_Index,Current_Theta, &
                              & Current_R,Current_Z,HC_Temperature,Current_Velocity,Current_Velocity_In_S,Current_Angle,Seed,NRand,&
                              & Current_Cell,Current_Ring, &
                              & Current_Velocity_In_R,Current_Velocity_In_Z,Current_S,Current_Cross,hc_v)

                         ! New position, energy/velocity, and direction are loaded.  Return to following code.
                         !Write (Output_Unit_Scratch,*) 'HC SPUTTERED from ion impact FP=4:',Current_Cell,Current_Ring,Target_Index,Current_R,Current_Z,gikds(Target_Index),girds(Target_Index),ryield,gkmfss(Target_Index), Target_Material,Dep_Energy,Sput_Weight

                         ! If output is requested, open files.
                         If (hc_evolve_print_option .eq. 1) Then
                            Write (Output_Unit_Evolve,9514) "Ion sputter in FP(4) at target index:",Target_Index,"R:",Current_R,"Z:&
                                 &",Current_R,"Last HC:",Last_HC_Species,"New HC:",Cur_HC_Spec
9514                        Format (3X,A,1X,I3,1X,A,1X,F8.3,1X,A,1X,F8.3,1X,A,1X,I2,1X,A,1X,I2)
                         End If
                         If (hc_coord_print_option .eq. 1) Then
                            Write (Output_Unit_Location,9515) "Ion sputter in FP(4) at target index:",Target_Index,"R:",Current_R,&
                                 &"Z:",Current_R,"Last HC:",Last_HC_Species,"New HC:",Cur_HC_Spec
9515                        Format (3X,A,1X,I3,1X,A,1X,F8.3,1X,A,1X,F8.3,1X,A,1X,I2,1X,A,1X,I2)
                         End If

                         ! 2.
                         ! jdemod - this calculation of the R,Z velocity components is performed 
                         !          in hc_ion_target_sputter_velocity called from hc_ion_impact_relaunch 
                         !          for reflected particles and does not need to be repeated here. 
                         ! Calculate R,Z components of velocity, probability of ionization, etc.
                         !Current_Velocity_In_R = Current_Velocity * COS (Current_Angle) *  Neutral_Time_Step
                         !Current_Velocity_In_Z = Current_Velocity * SIN (Current_Angle) *  Neutral_Time_Step

                         ! Continue following the newly sputtered particle.
                         Return

                      ElseIf (Yield_Calc .eq. 2) Then
                         ! DIVIMP yield determination.
                         If (gkmfss(Target_Index) .ge. 0.0) Then
                            RYield = YIELD (6,  Target_Material,Dep_Energy,gktebs (Current_Cell,Current_Ring),gktibs (Current_Cell,&
                                 & Current_Ring)) * gkmfss (Target_Index)
                         ElseIf (gkmfss (Target_Index) .lt. 0.0 .and. gkmfss (Target_Index) .ge. -50.0) then
                            RYield = ABS (gkmfss (Target_Index))
                         ElseIf (gkmfss(Target_Index) .le. -99.0) Then
                            RYield = YIELD (6, Target_Material,Dep_Energy,gktebs (Current_Cell,Current_Ring),gktibs (Current_Cell,&
                                 & Current_Ring))
                         End If

                         SputNew = Sput_Weight * RYield
                         HC_YldTot (Cur_HC_Spec,Launch_Reg) =  HC_YldTot (Cur_HC_Spec,Launch_Reg) + SputNew
                         HC_YldMax (Cur_HC_Spec,Launch_Reg) = MAX ( HC_YldMax (Cur_HC_Spec,Launch_Reg), SputNew)

                         WRITE(Output_Unit_Scratch,*) 'HC SPUTTERED from ion impact:',Current_Cell,Current_Ring,Target_Index,&
                              & Current_R,Current_Z,gikds (Target_Index),girds (Target_Index),ryield,gkmfss (Target_Index), &
                              & Target_Material,Dep_Energy,Sput_Weight,SputNew, HC_YldTot

                         If (SputNew .gt.  Self_Sputter_Threshold) Then
                            HC_YthTot (Cur_HC_Spec,Launch_Reg) =  HC_YthTot (Cur_HC_Spec,Launch_Reg) + SputNew
                            NProd  = NProd + 1
                            Call psnews(NPROD,SputNew)
                            If (hc_launch_angle_velocity .eq. 1 .or. hc_launch_angle_velocity .eq. 4 .or. hc_launch_angle_velocity &
                                 &.eq. 5 .or.  Sputter_Opt .eq. 4) Then
                               EMax =  EMax_Factor * Dep_Energy
                               Max_Velocity_Randoms = 1.0 / (1.0 +  Target_Binding_Energy / EMax)**2.0
                            Else
                               Max_Velocity_Randoms = 1.0
                            End If

                            Call pxprods (NProd,Current_R)
                            Call pyprods (NProd,Current_Z)

                            ! For segments with a fixed sputtering yield - allow for
                            ! the energy of the sputtered particle to be set to a
                            ! specific value.

                            If ( Self_Sputter_Opt .eq. 2 .and. gkmfss (Target_Index) .lt. 0.0) Then
                               Call peprods (NProd,hc_sput_energy_ion_preset)
                            Else
                               Call peprods (NProd,0.0)
                            End If

                            Call pidprods (NProd,Target_Index)
                            HC_Launchdat (NProd,2) = 1.0

                         End If

                         If (Debug) Then
                            Write (Output_Unit_Scratch,*) 'FP HC RELAUNCH:',Current_Cell,Current_Ring,Target_Index,Current_R,&
                                 & Current_Z,gikds (Target_Index),girds (Target_Index),SputNew,NProd,Dep_Energy
                         End If

                         ! Record collision with target.
                         HC_Num_Striking_Target(Cur_HC_Spec,Launch_Reg)=HC_Num_Striking_Target(Cur_HC_Spec,Launch_Reg)+Sput_Weight

                         IFate = 22 ! HC Ion Struck Target.
                         Return
                      End If
                   End If
                End If
             End If

             ! ELSEIF ( Far_Periphery_Opt .eq. 2) Then

             ! THIS OPTION ALLOWS THE ION TO DRIFT OUT AS FAR AS
             ! IT CAN. IT WILL ALWAYS BE ASSOCIATED WITH THE
             ! CHARACTERISTICS OF THE OUTERMOST RING.

          End If
       End If

       If ( Non_Orthogonal_Grid_Opt .eq. 1 .or.  Non_Orthogonal_Grid_Opt .eq. 3) Then
          ! Nonorth
          If (Current_Ring .ne. Last_Ring) Then

             Local_K = gkks (Current_Ring)
             CKK_Maximum = MAX ( CKK_Maximum,  Local_K)

             !  Set range for poloidal drift velocity
             If ( Poloidal_Drift_Opt .eq. 1) Then
                SOL_Drift_Start =  SOL_Drift_Start * gksmaxs (Current_Ring)
                SOL_Drift_End =  SOL_Drift_End * gksmaxs (Current_Ring)
             End If

             If (gtagdv (Current_Cell, Current_Ring) .eq. 0 .and. gtagdv (Last_Cell,Last_Ring) .eq. 0) Then
                !  BOTH POINTS UNSHIFTED FROM ORTHOGONAL GRID.
                If (Current_S .gt. gkss (Last_Cell,Last_Ring)) Then
                   Current_S = gkss (Current_Cell, Current_Ring) + (Current_S - gkss (Last_Cell,Last_Ring)) * &
                        & (gkfords (Current_Cell, Current_Ring) / gkfords (Last_Cell,Last_Ring))
                Else
                   Current_S = gkss (Current_Cell, Current_Ring) - (gkss (Last_Cell,Last_Ring) - Current_S) * &
                        & (gkbacds (Current_Cell, Current_Ring) / gkbacds (Last_Cell,Last_Ring))
                End If
             Else
                ! EITHER POINT SHIFTED FROM ORTHOGONAL GRID
                If (Current_Theta .lt. gthetag (Current_Cell, Current_Ring)) Then
                   If (Current_Cell .eq. 1) Then
                      If (Current_Theta .le. gthetat (gidds (Current_Ring, 2))) Then
                         Current_S = 0.0
                      Else
                         Current_S = gkss (Current_Cell, Current_Ring) * (Current_Theta - gthetat (gidds (Current_Ring, 2))) &
                              & / (gthetag (Current_Cell, Current_Ring) - gthetat (gidds (Current_Ring, 2)))
                      End If
                   Else
                      Temp_Theta = (gthetag (Current_Cell,Current_Ring) - Current_Theta) &
                           & / (gthetag (Current_Cell, Current_Ring) - gthetag (Current_Cell - 1, Current_Ring))
                      Current_S = gkss (Current_Cell,Current_Ring) - gkbacds (Current_Cell, Current_Ring) * Temp_Theta
                   End If
                Else
                   If (Current_Cell .eq. gnks (Current_Ring)) Then
                      If (Current_Theta .ge. gthetat (gidds (Current_Ring, 1))) Then
                         Current_S = gksmaxs (Current_Ring)
                      Else
                         Current_S = gkss (Current_Cell, Current_Ring) + gkfords (Current_Cell, Current_Ring) &
                              & * (gthetat (gidds (Current_Ring, 1)) - Current_Theta) &
                              & / (gthetat (gidds (Current_Ring, 1)) - gthetag (Current_Cell, Current_Ring))
                      End If
                   Else
                      Temp_Theta = (Current_Theta - gthetag (Current_cell, Current_Ring)) &
                           & / (gthetag (Current_Cell + 1, Current_Ring) - gthetag (Current_Cell, Current_Ring))
                      Current_S = gkss (Current_Cell, Current_Ring) + gkfords (Current_Cell, Current_Ring) * Temp_Theta
                   End If
                End If
             End If
          End If
          ! Nonorth
       End If

       ! Check to see if particle has been moved outside edges.
       If (Current_S .le. 0.0 .or. Current_S .ge. gksmaxs (Current_Ring)) Then
          ! If ((Current_S .le. 0.0 .or. Current_S .ge. gksmaxs (Current_Ring)) .and. Current_Ring .ge.  Inner_SOL_Ring) Then
          Write (Output_Unit_HC_Alert,*) "Warning: HC ion has been artificially moved outside of S in HC_INSIDE_ION:",Current_R,&
               &Current_Z,Current_S,Current_Cell,Current_Ring,gksmaxs (Current_Ring)
       End If

    End If

    ! Exit from main plasma.
    If (Last_Ring .lt.  Inner_SOL_Ring .and. Current_Ring .ge.  Inner_SOL_Ring) Then
       If ( HC_CflRex) Then

          ! Record exit event.
          HC_Tot_Fragments_Exit_Main (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Fragments_Exit_Main (Cur_HC_Spec,Launch_Reg) + Sput_Weight

          If ( Debug) Then
             Write (Output_Unit_Scratch,9003) IProd,Total_Steps,Current_Cell,Current_Ring,Cur_HC_Spec,Current_R,Current_Z,&
                  & Current_S, Local_K,Current_Theta, &
                  & gksmaxs (Current_Ring),Current_Velocity,HC_Temperature,Current_Cross,Sput_Weight,Time_Bin,'HC exited main'
          End If
          HC_CflRex = .False.

          HC_Num_Orig_Enter_Main (Cur_HC_Spec,Launch_Reg)  =  HC_Num_Orig_Enter_Main (Cur_HC_Spec,Launch_Reg) + Sput_Weight
          HC_Tot_Z_Orig_Enter_Main(Cur_HC_Spec,Launch_Reg)=HC_Tot_Z_Orig_Enter_Main(Cur_HC_Spec,Launch_Reg)+Current_Z*Sput_Weight
          HC_Tot_S_Orig_Enter_Main (Cur_HC_Spec,Launch_Reg) =  HC_Tot_S_Orig_Enter_Main (Cur_HC_Spec,Launch_Reg) + MIN (Current_S, &
               & gksmaxs (Current_Ring) - Current_S) * Sput_Weight
          HC_ELims (Current_Cell,2,Cur_HC_Spec) =  HC_ELims (Current_Cell,2,Cur_HC_Spec) + Sput_Weight

          !				If ( HC_InMain) Then
          !					 HC_Cmmm (Cur_HC_Spec,Launch_Reg)  =  HC_Cmmm (Cur_HC_Spec,Launch_Reg) + Sput_Weight
          !					 HC_CmmmX (Cur_HC_Spec,Launch_Reg) =  HC_CmmmX (Cur_HC_Spec,Launch_Reg) + Current_Z * Sput_Weight
          !					 HC_CmmmS (Cur_HC_Spec,Launch_Reg) =  HC_CmmmS (Cur_HC_Spec,Launch_Reg) + MIN (Current_S, gksmaxs (Current_Ring) - Current_S) * Sput_Weight
          !					 HC_ELims (Current_Cell,2,Cur_HC_Spec) =  HC_ELims (Current_Cell,2,Cur_HC_Spec) + Sput_Weight
          !				Else
          !					 HC_Clll (Cur_HC_Spec,Launch_Reg)  =  HC_Clll (Cur_HC_Spec,Launch_Reg) + Sput_Weight
          !					 HC_ClllX (Cur_HC_Spec,Launch_Reg) =  HC_ClllX (Cur_HC_Spec,Launch_Reg) + Current_Z * Sput_Weight
          !					 HC_ClllS (Cur_HC_Spec,Launch_Reg) =  HC_ClllS (Cur_HC_Spec,Launch_Reg) + MIN (Current_S, gksmaxs (Current_Ring) - Current_S) * Sput_Weight
          !					 HC_ELims (Current_Cell,1,Cur_HC_Spec) =  HC_ELims (Current_Cell,1,Cur_HC_Spec) + Sput_Weight
          !				End If
       End If
    End If

    !  BOTH ROUTES CONTINUE HERE ...
    !write (Output_Unit_Scratch,*) "both routes continue here"
    ! Check for leakage - record particle once.
    If (gcheckleak() .and. Current_Ring .ge.  Inner_SOL_Ring &
         & .and. (.not.  HC_Has_Leaked .and. &
         & (Current_S .gt. gcleaks ( HC_Leak_Particles (Launch_Reg)) .and. Current_S .lt. (gksmaxs (Current_Ring) - &
         & gcleaks ( HC_Leak_Particles (Launch_Reg)))))) Then
       HC_Leak_Density ( HC_Leak_Particles (Launch_Reg), Cur_HC_Spec) = &
            &  HC_Leak_Density ( HC_Leak_Particles (Launch_Reg), Cur_HC_Spec) + Sput_Weight
       HC_Leak_Particles (Launch_Reg) =  HC_Leak_Particles (Launch_Reg) + 1
       If ( HC_Leak_Particles (Launch_Reg) .gt. gcleaksn ()) Then
          HC_Leak_Particles (Launch_Reg) = gcleaksn ()
          HC_Has_Leaked = .True.
          HC_Leak_Time (Launch_Reg) =  HC_Leak_Time (Launch_Reg) + Total_Steps *  Ion_Time_Step
       End If
    End If

    If (Random_Numbers_Used .gt. (Max_Random_Values + 10)) Then
       !Write (Output_Unit_Scratch,*) 'HC 3:',Random_Numbers_Used,Max_Random_Values
    End If

    ! CHECK FOR COLLISION.
    Random_Numbers_Used = Random_Numbers_Used + 1
    If ((HC_Temperature * granv (Random_Numbers_Used)) .le.  Local_HC_Change_State_Coll_Prob) Then
       HC_Tot_Ion_Collision (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Ion_Collision (Cur_HC_Spec,Launch_Reg) + Sput_Weight
    End If


    !
    ! Score particle in DDLIMS in Cur_HC_Spec position.
    ! If time point reached, score time position also, and increment.
    ! Note:  Do only for C neutral.
    !If (Cur_HC_Spec .eq. 2) Then
    !   HC_DDLims (Current_Cell,Current_Ring,Get_HC_Charge (Cur_HC_Spec)) =  HC_DDLims (Current_Cell,Current_Ring,Get_HC_Charge (&
    !                                                                      & Cur_HC_Spec)) + Double_Sput_Weight
    !   HC_DDts (Current_Cell,Current_Ring,Get_HC_Charge (Cur_HC_Spec)) =  HC_DDts (Current_Cell,Current_Ring,Get_HC_Charge (&
    !                                                                    & Cur_HC_Spec)) + Double_Sput_Weight * HC_Temperature
    !End If

    ! Do also for CH/CD neutral.
    !If (Cur_HC_Spec .eq. 4) Then
    !   HC_DDLims_CH (Current_Cell,Current_Ring,Get_HC_Charge (Cur_HC_Spec)) =  HC_DDLims_CH (Current_Cell,Current_Ring,&
    !                                                                         & Get_HC_Charge (Cur_HC_Spec)) + Double_Sput_Weight
    !   HC_DDts_CH (Current_Cell,Current_Ring,Get_HC_Charge (Cur_HC_Spec)) =  HC_DDts_CH (Current_Cell,Current_Ring,Get_HC_Charge (&
    !                                                                       & Cur_HC_Spec)) + Double_Sput_Weight * HC_Temperature
    !End If

    ! Determine the time bin for the particle
    If ( Number_Time_Steps .gt. 0) Then
       Time_Bin = IPOS (Total_Steps, gctimes (1,Cur_HC_Spec), Total_Steps + 1)
       ! Record time dependent contribution at each time step
       HC_Lims(Current_Cell,Current_Ring,Cur_HC_Spec,Time_Bin)=HC_Lims(Current_Cell,Current_Ring,Cur_HC_Spec,Time_Bin)+Sput_Weight
    End If

    ! Position extremes for current particle.
    Min_Z (Cur_HC_Spec, Launch_Reg) = MIN (Min_Z (Cur_HC_Spec, Launch_Reg), Current_Z)
    Max_S (Cur_HC_Spec, Launch_Reg) = MAX (Max_S (Cur_HC_Spec, Launch_Reg), MIN (Current_S,gksmaxs (Current_Ring) - Current_S))

    Last_Z = Current_Z

    If ( Ion_Diffuse) Then
       If ( Local_SPara .le. 0.0) Then
          HC_DParas (Cur_HC_Spec,1,Launch_Reg) =  HC_DParas (Cur_HC_Spec,1,Launch_Reg) + Double_Sput_Weight
          If (Current_Ring .ge.  Inner_SOL_Ring) Then
             HC_DParas (Cur_HC_Spec,5,Launch_Reg) =  HC_DParas (Cur_HC_Spec,5,Launch_Reg) + Double_Sput_Weight
          End If
       Else
          ! jdemod - need to add launch_reg to HC_Dparas
          !HC_DParas (Cur_HC_Spec, 2) =  HC_DParas (Cur_HC_Spec, 2) + Double_Sput_Weight
          HC_DParas (Cur_HC_Spec, 2,launch_reg) =  HC_DParas (Cur_HC_Spec, 2,launch_reg) + Double_Sput_Weight
          If (Current_Ring .ge.  Inner_SOL_Ring) Then
             HC_DParas (Cur_HC_Spec, 3,Launch_Reg) =  HC_DParas (Cur_HC_Spec, 3,Launch_Reg) + Double_Sput_Weight
             HC_DParas(Cur_HC_Spec,4,Launch_Reg)=HC_DParas(Cur_HC_Spec,4,Launch_Reg)+Double_Sput_Weight*DBLE(Local_SPara)
          End If
       End If
    End If

    ! IONISATION, RECOMBINATION and Hydrocarbon evolution.

    ! Check reaction possibilities for current state.
    If (Get_HC_Charge (Cur_HC_Spec) .eq. 0) Then
       HC_TimeStep =  Neutral_Time_Step
    Else
       HC_TimeStep =  Ion_Time_Step
    End If


    ! Zero out necessary arrays for all states posible.
    States_Possible = 0
    Reactions_Possible = 0
    !Total_Probabilities = 0.0
    NS_Probabilities = 0.0
    Reaction_Probabilities = 0.0
    Total_probability = 0.0
    number_reactions = 0

    ! Find total hydrocarbon state change probability for current species, timestep, and grid location.
    ! Note that E&L reaction rates are in cm^3/s.  Therefore, DIVIMP plasma
    ! densities must be reduced by 1E6 (from 1/m^3) for use.
    ! jdemod - removed states_possible from call
    Call HC_New_State (Current_Cell,Current_Ring,Cur_HC_Spec,HC_Data_Type,gknbs(Current_Cell,Current_Ring) / 1.0E6,gktebs(&
         &  Current_Cell,Current_Ring),gktibs(Current_Cell,Current_Ring), &
         &  Back_Plasma_Ion_Mass,HC_Timestep,Reactions_Possible,Reaction_Probabilities,NS_Probabilities,&
         &  Total_Probability,number_reactions)

    !Write (Output_Unit_Scratch,*) "NEWST ION :",Cur_HC_Spec,gknbs (Current_Cell,Current_Ring) / &
    !           & 1.0E6,gktebs (Current_Cell,Current_Ring),gktibs (Current_Cell,Current_Ring), &
    !           &  Back_Plasma_Ion_Mass,HC_Timestep,Reactions_Possible,States_Possible,&
    !           &  Reaction_Probabilities,NS_Probabilities,Total_Probability

    !Write (Output_Unit_Scratch,'(a,i6,2x,a,i6,g14.6)') "NEWST ION :",number_reactions,&
    !      & hc_state_table(cur_hc_spec)%state_name, Cur_HC_Spec,total_probability
    !do in = 1,number_reactions
    !   write(output_unit_scratch,'(a,i6,3(1x,g12.5))') &
    !      &hc_state_table(hc_state_transform_table(reactions_possible(in))%end_c_states(1))%state_name,in,&
    !      &Reactions_Possible(in),Reaction_Probabilities(in),NS_Probabilities(in)
    !end do

    ! jdemod - this should be an HC debugging print out - not a particle debug printout
    !If (Debug_hc) Then
    !   Write (Output_Unit_Scratch,'(a,i5,2x,a,g12.4,f10.3,f10.3,f8.5,f12.1)') "New State:",Cur_HC_Spec, HC_Data_Type, gknbs (Current_Cell, Current_Ring), &
    !        & gktebs(Current_Cell, Current_Ring), gktibs(Current_Cell, Current_Ring), total_probability, HC_Timestep
    !End If

    ! Check for allowance of hydrocarbon state transitions.
    If (hc_disable_transitions .eq. 1) Then
       ! Possibility of MTC collisions for neutrals results in differing treatment for ions vs. neutrals
       ! Transitions disabled.
       !Total_Probability = 0.0
       Return
    End If

    ! Multiply total probability of ion reactant transition to reduce time spent as ionized species (if multiplier >1.0).
    ! Note, HC_Ion_Reaction_Mult is set in subroutine HC_Begin, found in HC_Start.
    ! jdemod - NOTE: I'm not sure what this is used for - current value is set to 1.0 so it won't change anything
    !                It could be used to artificially speed the ion reaction processes ... ???
    If (Get_HC_Charge (Cur_HC_Spec) .ne. 0) Then
       Total_Probability = Total_Probability * HC_Ion_Reaction_Mult
       If (Total_Probability .gt. 1.0) Then
          Total_Probability = 1.0
       End If
    End If

    !write (0,*) "HERE",Ion_Time_Step_Count
    ! Force change of state at particular timestep.
    !If (Ion_Time_Step_Count .eq. 10000) then
    !Total_Probability = 1.0
    !Ion_Time_Step_COunt = Ion_Time_Step_Count + 1
    !else
    !Total_Probability = 0.0
    !endif
    !Total_Probability = Total_Probability / 5.0

    Random_Numbers_Used = Random_Numbers_Used + 1			
    If (granv (Random_Numbers_Used) .ge. Total_Probability) Then
       ! No reaction or MTC occurs.
       Tries_Without_Change = Tries_Without_Change + 1
       If (Tries_Without_Change .gt.  Max_HC_Ion_Iter_To_TMax) Then
          Write (Output_Unit_HC_Alert,*) "Error in HC_Inside_Ion:  Probability for reaction too low:", Total_Probability, &
               &Max_HC_Ion_Iter_To_TMax
          Write (Output_Unit_HC_Alert,*) "Program stopping."
          Stop
       End If
       ! Exit this routine if no change of state occurs ...
       Return
    End If

    ! Proceed to evolve hydrocarbon.
    ! Save previous state.
    Last_HC_Species = Cur_HC_Spec
    Last_H_Isotope_Composition = H_Isotope_Composition

    ! Get random number to decide which state is the destination.
    Random_Numbers_Used = Random_Numbers_Used + 1 ! Add to KK counter.
    Random_Temp = granv (Random_Numbers_Used)


    ! jdemod - why do this here? Why not perform this calculation in the probabilities routines?
    ! Convert normalized probability values to sums of state change given all previous reactions.
    ! jdemod - cumulative probabilities now returned by HC new state
    !Do Counter=2, size (NS_Probabilities), 1
    !   NS_Probabilities (Counter) = NS_Probabilities (Counter) + NS_Probabilities (Counter-1)
    !End Do

    ! Find which state we've moved into with another random value.
    ! jdemod - a linear search is not the most efficient
    Counter=1
    Do While (Random_Temp .gt. NS_Probabilities (Counter))
       ! Increment counter.
       Counter = Counter + 1
    End Do

    ! Reset tries without change counter.
    Tries_Without_Change = 0

    ! jdemod
    ! Load states possible from hc_state_transform_table - process selected is indicated by counter
    do in = 1, highest_carbon_content
       states_possible(in) = hc_state_transform_table(reactions_possible(counter))%end_c_states(in)
    end do

    ! Indicate changed state.  Note:  Always take the first carbon atomic/molecular product,
    ! then check for additional components and store in the stack.
    Followed_State = 1
    ! jdemod
    Cur_HC_Spec = States_Possible (Followed_State)
    !Cur_HC_Spec = States_Possible (Couter,Followed_State)

    !Cur_HC_Spec = 8
    !write (0,*) "Changing to 8"

    If (Cur_HC_Spec .eq. 0) Then
       Write (Output_Unit_HC_Alert,*) "Error in HC_Inside_Ion:",Reactions_Possible, States_Possible, Reaction_Probabilities, &
            & NS_Probabilities	
       Write (Output_Unit_HC_Alert,*) "Program stopping."
       Stop
    End If
    !write (Output_Unit_Scratch,*) "HC Change Of State from ion:",Last_HC_Species,Cur_HC_Spec,Current_R,Current_Z,Current_Cell,Current_Ring, HC_MTCProb (Current_Cell,Current_Ring,Cur_HC_Spec),Total_Probability

    ! Record left-over hydrogen based on background plasma/mass reduction and update H_Isotope_Composition.
    Call Record_Hydrogen_Release (Current_Cell,Current_Ring,Last_HC_Species,Reactions_Possible (Counter),Cur_HC_Spec,&
         &H_Isotope_Composition,Seed,NRand)

    ! Record reaction indicating a change of state.
    HC_Reaction_Count (Last_HC_Species,Reactions_Possible (Counter),Launch_Reg) =  HC_Reaction_Count (Last_HC_Species,&
         &Reactions_Possible (Counter),Launch_Reg) + 1.0

    Sput_Weight = Calc_Sputy (Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition))
    Double_Sput_Weight = DBLE (Sput_Weight)

    ! Modify commons for cell tau parallel, stopping, and heating for current cell.
    Call Initialize_State_Prop_Data (Current_Cell,Current_Ring,Get_HC_Charge (Cur_HC_Spec))
    Call Modify_Taus (Cur_HC_Spec,H_Isotope_Composition,Sput_Weight,Current_Cell,Current_Ring)
    ! jdemod
    !Call Modify_Taus (Cur_HC_Spec,H_Isotope_Composition,Sput_Weight,HC_Temperature,Current_Cell,Current_Ring,Current_S,Seed,&
    !                 &Random_Numbers_Used,NRand)

    ! Now, load additional particles from higher hydrocarbon ion breakup into stack storage.
    Add_Followed_State = Followed_State + 1

    Do While (Add_Followed_State .le. Highest_Carbon_Content)
       ! Check if molecule broke into more than one carbon-containing product.
       ! jdemod 
       !If (States_Possible (Counter, Add_Followed_State) .ne. 0) Then
       If (States_Possible (Add_Followed_State) .ne. 0) Then
          ! Found another carbon product.  Save required characteristics on the stack.
          Write (Output_Unit_HC_Alert,*) "SHOULD never be here for E&L database - humungous error"

          ! jdemod
          !Call HC_Push (States_Possible (Counter, Add_Followed_State), Current_R, Current_Z, Current_S, Current_Cross, &
          Call HC_Push (States_Possible (Add_Followed_State), Current_R, Current_Z, Current_S, Current_Cross, &
               &Current_Angle, Current_Theta, Current_Velocity, Current_Velocity_In_S, HC_Temperature)

          ! Increment counter as there may be another carbon-carrying component.
          Add_Followed_State = Add_Followed_State + 1
       Else
          ! No carbon-containing particles found.  Exit the loop.
          Exit
       End If

    End Do

    !---------------------------------------------------------------------------------------------
    !
    !------------------------ ION TO ION ---------------------------------------------------------
    !
    !---------------------------------------------------------------------------------------------

    ! ION TO ION.  Check for ionization and recombination and alter grid coordinates accordingly.
    If((Last_HC_Species.ne.Cur_HC_Spec).and.(Get_HC_Charge(Last_HC_Species).ne.0).and.(Get_HC_Charge(Cur_HC_Spec).ne.0))Then
       ! Particle began as ion, stayed an ion.
       Old_Cell = Current_Cell
       Old_Ring = Current_Ring
       Old_S = Current_S
       Old_Cross = Current_Cross
       Old_R = Current_R
       Old_Z = Current_Z
       Old_Theta = Current_Theta
       Old_Angle = Current_Angle
       Old_Velocity = Current_Velocity
       Old_Velocity_In_S = Current_Velocity_In_S
       Old_Temp = HC_Temperature

       ! jdemod
       ! Update the particle energy, velocity, temperature and trajectory as required - 
       !     option 0 is the original code
       !     Option 1 is the revised hc_kinetics code
       !

       if (hc_kinetics_opt.eq.0) then 

          Call Ion_Ion_FreeSpace_Location (Reactions_Possible (Counter),Followed_State,Current_Cell,Current_Ring,Current_S,&
               &Current_Cross,Current_Theta,S_Start)
          Call Ion_Ion_FreeSpace_Angle (Reactions_Possible (Counter),Followed_State,Current_Angle)
          ! Krieger IPP/07 - SUN compiler insists on 132 column limit
          Call Ion_Ion_FreeSpace_Energy (Reactions_Possible (Counter),Followed_State,H_Isotope_Composition,Current_Cell, &
               & Current_Ring,Current_Velocity_In_S,HC_Temperature,Kin_Energy_Added)
          Call Ion_Ion_FreeSpace_Velocity (Reactions_Possible (Counter),Followed_State,H_Isotope_Composition,Current_Ring, &
               & Current_S,Current_Velocity_In_S,Seed,NRand)

       elseif (hc_kinetics_opt.eq.1) then 

          if (debug_kinetics) &
               & write(6,'(a,10g12.5)') 'TP:INSIDE ION :',total_probability,hc_ion_reaction_mult

          call update_hc_kinetics(hc_v,ii_reaction,current_cell,current_ring,  & 
               & cur_hc_spec,reactions_possible(counter),followed_state,h_isotope_composition,  &
               & current_r,current_z,current_s, current_theta,current_cross,current_angle,  &
               & current_velocity_in_s,current_velocity_in_r,current_velocity_in_z, &
               & current_velocity,kin_energy_added,s_start, &
               & neutral_time_step,ion_time_step,sput_weight,iprod)

       endif

       ! Print hydrocarbon evolution data to file if option H31 selected in input file.
       If (HC_Evolve_Print_Option .eq. 1) Then
          ! Print option activated.
          Time_Tot = Neutral_Time_Step_Count *  Neutral_Time_Step + Ion_Time_Step_Count *  Ion_Time_Step
          Write (Output_Unit_Evolve,9300) Time_Tot,Neutral_Time_Step_Count,Ion_Time_Step_Count,Last_HC_Species,&
               & Last_H_Isotope_Composition(1),Last_H_Isotope_Composition(2), &
               & Last_H_Isotope_Composition(3),Cur_HC_Spec,H_Isotope_Composition(1),H_Isotope_Composition(2),H_Isotope_Composition(&
               & 3),Current_Cell,Current_Ring, &
               & Old_R,Current_R,Old_Z,Current_Z,Old_S,Current_S,Old_Cross,Current_Cross,Old_Angle,Current_Angle,Old_Velocity_In_S,&
               & Current_Velocity_In_S,Old_Temp,HC_Temperature
9300      Format (E10.4,F10.1,F10.1,I10,2X,I2,1X,I2,1X,I2,I10,2X,I2,1X,I2,1X,I2,I10,I10,F10.6,F10.6,F10.6,F10.6,F10.6,F10.6,F10.4,&
               &F10.4,F10.4,F10.4,F10.3,F10.3,F10.3,F10.3,F10.3,F10.3)
       End If





       ! Note, Launch_Reg is not updated here to keep it at the value where the original
       ! hydrocarbon was launched into the plasma.

       ! CALCULATE TIME AT WHICH DIFFUSION WILL BE FIRST APPLIED.	
       Local_RConst =  Calc_Hi

       If ( First_Diffuse_Opt .eq. 0) Then
          If ( Collision_Opt .ne. 1) Then
             Local_RConst = 0.0
          End If

       ElseIf ( First_Diffuse_Opt .eq. 1) Then
          If ( Local_HC_Change_State_Coll_Prob .gt. 0.0) Then
             NRand = NRand + 1
             Call Surand2 (Seed, 1, Random_Value)
             Local_RConst = -HC_Temperature * LOG (Random_Value) /  Local_HC_Change_State_Coll_Prob
          End If
       End If




       If (Debug) Then
          Write (Output_Unit_Scratch,9005)
          Write (Output_Unit_Scratch,9003) IProd,0.0,Current_Cell,Current_Ring,Cur_HC_Spec,Current_R,Current_Z,Current_S, Local_K,&
               & Current_Theta,gksmaxs (Current_Ring), &
               & Current_Velocity_In_S,HC_Temperature, Local_SPara,Current_Cross,Sput_Weight,Time_Bin,'HC Ion Appeared'
9003      FORMAT(1X,I5,F9.1,2I3,I2,2F9.5,F8.3,2F6.2,F8.3,1P,E15.8,0P,F7.1,1P,E8.1,0P,F8.5,F5.2,I2,:,1X,A,:,F8.5)
9005      FORMAT(1X,'--HC------TIME-IK-IR-HC','----R--------Z-------S------K---THETA--SMAX---','---DRIFTVEL------TEMI-PARADIFF-&
               &CROSS--FRAC-IT',14('-'))
       End If

       If (100*(IProd/100) .eq. IProd) Then
          Write (Output_Unit_Scratch,'('' DIV: HC Ion from molecule'',I6,'' STARTING'')') IProd
       End If

       HC_CflRxa = .True.
       HC_CflRex = .True.
       If (Current_Ring .ge.  Inner_SOL_Ring) Then
          HC_InMain = .False.
          HC_CflRin = .True.
       Else
          HC_InMain = .True.
          HC_CflRin = .False.
          HC_Tot_Start_Main =  HC_Tot_Start_Main + Sput_Weight
          If ( Stop_Ion_In_Core_Opt .eq. 1) Then
             HC_Stopped_Follow_Ion_In_Core(Cur_HC_Spec,Launch_Reg)=HC_Stopped_Follow_Ion_In_Core(Cur_HC_Spec,Launch_Reg)+Sput_Weight
             IFate = 21 ! HC Ion Reached Main Plasma.
             Return
          End If
       End If

       ! Calculate starting region and set appropriate logical
       ! variable and update the appropriate storage array. Then
       ! as the particle enters each of the other regions - set
       ! the logical variable to true and accumulate the count
       ! in the appropriate starting bin.
       Launch_Cell = Current_Cell
       Launch_Ring = Current_Ring

       !
       ! jdemod - the code appears to reinitilize the particle in grid region
       !          tracking toggles whenever a new state is entered. 
       !        - this is inconsistent with DIVIMP where these tracking values
       !          are set when the ion is created. 
       !        - the idea behind these variables was to track core leakage
       !          in terms of particles that crossed the separatrix at least 
       !          once and compare this to where these particles formed. 
       !        - For the CH module this should either be based on where the CH4 originates
       !          or solely on the C+ creation and not on the intermediate CH4 states. 
       !        - The simplest approach is probably to take the C+ creation point and use that
       !          as the starting point as is done with neut. The expectation would be that the 
       !          amount of CH4 fragments entering and leaving the confined plasma should be 
       !          relatively small. 
       !        - In any case, recording data for all of the different states without saving the 
       !          data by state is not useful since some particles may enter and leave the 
       !          core or other grid regions multiple times. 
       !        - Two fixes: 1) Do not copy HC_wtsource to wtsource
       !                     2) Track particle grid entry by CH4 injected and not by each breakup state. 
       !                        - comment out tracking toggle initialization here so that current
       !                        values carry over. 
       !
       !HC_InCore = .False.
       !HC_InEdge = .False.
       !HC_InMSOL = .False.
       !HC_InDIV = .False.
       !HC_InTrap = .False.

       If ( Launch_Ring .lt.  Inner_SOL_Ring) Then
          ! Particle leakage to core.
          If (.not.  HC_InCore) Then
             HC_InCore = .True.
             HC_Num_Entered_Core =  HC_Num_Entered_Core + Sput_Weight

             HC_Ion_Core_Density ( Launch_Cell, Launch_Ring) = &
                  &  HC_Ion_Core_Density ( Launch_Cell, Launch_Ring) + Sput_Weight

             If (NeutType .gt. 0) Then
                ! Record occurance of leakage.
                HC_WTSource ( Starting_Index, Launch_Ring,3,NeutType) = &
                     &  HC_WTSource ( Starting_Index, Launch_Ring,3,NeutType) + Sput_Weight

                ! Record leakage position.
                HC_WTSource ( Starting_Index, Launch_Ring,4,NeutType) = &
                     &  HC_WTSource ( Starting_Index, Launch_Ring,4,NeutType) + Sput_Weight * S_Start
             End If
          End If
       Else
          ! In the edge so set some sub-divisions
          HC_InEdge = .True.

          HC_Ion_Edge_Density ( Launch_Cell, Launch_Ring) = &
               &  HC_Ion_Edge_Density ( Launch_Cell, Launch_Ring) + Sput_Weight

          If (Current_Ring .gt.  Inner_Wall_Ring .and. Current_Ring .le.  Num_Upper_Rings) Then
             ! Particle starting in trap region
             HC_InTrap = .True.

             HC_Ion_Trap_Density ( Launch_Cell, Launch_Ring) = &
                  &  HC_Ion_Trap_Density ( Launch_Cell, Launch_Ring) + Sput_Weight

          ElseIf ((Current_Z .ge.  Z_X_Point .and.  ReFCT .eq. 1) .or. (Current_Z .le.  Z_X_Point .and.  ReFCT .eq. 0)) Then
             ! Divertor region
             HC_InDIV = .True.

             HC_Ion_Divertor_Density ( Launch_Cell, Launch_Ring) = &
                  &  HC_Ion_Divertor_Density ( Launch_Cell, Launch_Ring) + Sput_Weight

          Else
             ! Main SOL Region
             HC_InMSOL = .True.

             HC_Ion_MSOL_Density ( Launch_Cell, Launch_Ring) = &
                  &  HC_Ion_MSOL_Density ( Launch_Cell, Launch_Ring) + Sput_Weight

          End If
       End If

       ! Check for prompt deposition.
       ! jdemod - HC_prompt_deposition needs to reset last_hc_species as well as the current one or it can trigger
       !          change of state code found below - problem addressed by placing ion->neutral transition code in ELSEIF block opposite
       !          this code
       Call HC_Prompt_Deposition (		&
            &	NProd,				&
            &	Current_Cell,			&
            &	Current_Ring,			&
            &	Current_R,			&
            &	Current_Z,			&
            &	Current_S,			&
            &	Current_Cross,			&
            &	Current_Theta,			&
            &	Cur_HC_Spec,	           	&
            &	H_Isotope_Composition,		&
            &	Current_Velocity_In_S,		&
            &	Current_Velocity,		& ! VIN
            &	Current_Velocity_In_R,		& ! XVELF
            &	Current_Velocity_In_Z,  	& ! YVELF
            &	HC_Temperature,			&
            &	Last_HC_Species,		&
            &	Last_H_Isotope_Composition,	&
            &	Launch_Reg,			& ! M
            &	Reflection_Counter,		&
            &	NeutType,			&
            &	Max_Velocity_Randoms,		& ! RMAXS
            &	Tries_Without_Change,		& ! CISTIZ
            &	Neutral_Time_Step_Count,	&
            &	Ion_Time_Step_Count,		&
            &	Eq_Total_Ion_Time_Steps,	& ! CIST
            &	Seed,				&
            &	NRand,				&
            &	IFate,				&
            &   hc_v,                           & ! hc velocity data
            &	 Debug_HC_Prompt)				

       ! jdemod
       !End If

       !---------------------------------------------------------------------------------------------
       !
       !------------------------ ION TO NEUTRAL -----------------------------------------------------
       !
       !---------------------------------------------------------------------------------------------


       ! jdemod - make this part of an elseif in the above so that a promptly redeposited particle can't
       !          trigger this code
       ! ION TO NEUTRAL.  Check for recombination and alter grid coordinates accordingly.
    ElseIf((Last_HC_Species.ne.Cur_HC_Spec).and.(Get_HC_Charge(Last_HC_Species).ne.0).and.(Get_HC_Charge(Cur_HC_Spec).eq.0))Then
       ! Particle began as ion, neutralization occurred.
       Old_Cell = Current_Cell
       Old_Ring = Current_Ring
       Old_S = Current_S
       Old_Cross = Current_Cross
       Old_R = Current_R
       Old_Z = Current_Z
       Old_Theta = Current_Theta
       Old_Angle = Current_Angle
       Old_Velocity = Current_Velocity
       Old_Velocity_In_S = Current_Velocity_In_S
       Old_Temp = HC_Temperature

       ! jdemod - add new kinetics options
       if (hc_kinetics_opt.eq.0) then

          Call Ion_Neut_FreeSpace_Location (Reactions_Possible (Counter),Followed_State,Current_Cell,&
               &Current_Ring,Current_S,&
               &Current_Cross,Current_R,Current_Z,Current_Theta,S_Start)
          Call Ion_Neut_FreeSpace_Angle (Reactions_Possible (Counter),Followed_State,Current_Cell,&
               &Current_Ring,Current_Angle,&
               &Azimuthal_Angle,Seed,NRand)
          Call Ion_Neut_FreeSpace_Energy (Reactions_Possible (Counter),Followed_State,&
               &H_Isotope_Composition,Current_Cell,Current_Ring,&
               &Current_Velocity_In_S,HC_Temperature,Kin_Energy_Added)
          Call Ion_Neut_FreeSpace_Velocity (Reactions_Possible (Counter),Followed_State,&
               &H_Isotope_Composition,Current_Angle,&
               & Azimuthal_Angle, &
               & Current_Cell,Current_Ring,HC_Temperature,Current_Velocity_In_S,Current_Velocity,Current_Velocity_In_R,&
               & Current_Velocity_In_Z,Seed,NRand)

       elseif (hc_kinetics_opt.eq.1) then 

          call update_hc_kinetics(hc_v,in_reaction,current_cell,current_ring,  & 
               & cur_hc_spec,reactions_possible(counter),followed_state,h_isotope_composition,  &
               & current_r,current_z,current_s, current_theta,current_cross,current_angle,  &
               & current_velocity_in_s,current_velocity_in_r,current_velocity_in_z, &
               & current_velocity,kin_energy_added,s_start, &
               & neutral_time_step,ion_time_step,sput_weight,iprod)

       endif


       ! Print hydrocarbon evolution data to file if option H31 selected in input file.
       If (HC_Evolve_Print_Option .eq. 1) Then
          ! Print option activated.
          Time_Tot = Neutral_Time_Step_Count *  Neutral_Time_Step + Ion_Time_Step_Count *  Ion_Time_Step
          Write (Output_Unit_Evolve,9302) Time_Tot,Neutral_Time_Step_Count,Ion_Time_Step_Count,Last_HC_Species,&
               & Last_H_Isotope_Composition(1),Last_H_Isotope_Composition(2), &
               & Last_H_Isotope_Composition(3),Cur_HC_Spec,H_Isotope_Composition(1),H_Isotope_Composition(2),H_Isotope_Composition(&
               & 3),Current_Cell,Current_Ring, &
               & Old_R,Current_R,Old_Z,Current_Z,Old_S,Current_S,Old_Cross,Current_Cross,Old_Angle,Current_Angle,Old_Velocity_In_S,&
               & Current_Velocity,Old_Temp,HC_Temperature
9302      Format (E10.4,F10.1,F10.1,I10,2X,I2,1X,I2,1X,I2,I10,2X,I2,1X,I2,1X,I2,I10,I10,F10.6,F10.6,F10.6,F10.6,F10.6,F10.6,F10.4,&
               &F10.4,F10.4,F10.4,F10.3,F10.3,F10.3,F10.3,F10.3,F10.3)
       End If

       ! Update K, Cell S centre, and ring SMax.
       Call Initialize_Cell_Geom_Data (Current_Cell, Current_Ring)

       If (Debug) Then
          Write (Output_Unit_Scratch,9003) IProd,0.0,Current_Cell,Current_Ring,Cur_HC_Spec,Current_R,Current_Z,Current_S, Local_K,&
               & Current_Theta,gksmaxs (Current_Ring), &
               & Current_Velocity_In_S,HC_Temperature, Local_SPara,Current_Cross,Sput_Weight,Time_Bin,'HC Ion Appeared'
       End If

       If (100*(IProd/100) .eq. IProd) Then
          !Write (Output_Unit_Scratch,'('' DIV: HC Neutral from molecule'',I6,'' STARTING'')') IProd
       End If

       ! Finished transition, return to following code HC_Inside_Neutral.

    End If

    If (Cur_HC_Spec .eq. 1) Then ! Reached "C+".
       ! C+ production has occured.  Store particle details in arrays and totals.

       If (NeutType .eq. 2 .or. NeutType .eq. 5) Then
          HC_ChemIzs (Current_Cell,Current_Ring,Cur_HC_Spec) =  HC_ChemIzs (Current_Cell,Current_Ring,Cur_HC_Spec) + Sput_Weight
       End If

       NatIZ = NatIZ + 1
       HC_Num_Fragments_Reach_CIon (Launch_Reg) =  HC_Num_Fragments_Reach_CIon (Launch_Reg) + Sput_Weight ! HC routine variable in HC_Launch_Diag_Table, same as RatIZ / SatIZ.
       !write (0,*) "adding one to ion ionization:", HC_Num_Fragments_Reach_CIon (Launch_Reg),Launch_Reg
       Call pxatizs (LatIZ + NatIZ - 1,Current_R)
       Call pyatizs (LatIZ + NatIZ - 1,Current_Z)
       Call pkatizs (LatIZ + NatIZ - 1, Local_K)
       Call psatizs (LatIZ + NatIZ - 1,MIN (Current_S, gksmaxs (Current_Ring) - Current_S))
       Call ptemtizs (LatIZ + NatIZ - 1,HC_Temperature)
       Call psnews (LatIZ + NatIZ - 1,Sput_Weight)
       Call pcistizs (LatIZ + NatIZ - 1,REAL (Eq_Total_Ion_Time_Steps) *  Neutral_Time_Step /  Ion_Time_Step)

       ! Record starting wall or target element for neutral/ion.
       Call pidatizs (LatIZ + NatIZ - 1,1, Starting_Index)
       Call pidatizs (LatIZ + NatIZ - 1,2,NeutType)
       Call pidatizs (latIZ + NatIZ - 1,3, First_Ioniz_Ring) ! irstart.
       Call pidatizs (latIZ + NatIZ - 1,4, First_Ioniz_Cell) ! ikstart.
       Call ptravel_locations (latIZ + NatIZ - 1,1, HC_InCore) ! incore.
       Call ptravel_locations (latIZ + NatIZ - 1,2, HC_InEdge) ! inedge.
       Call ptravel_locations (latIZ + NatIZ - 1,3, HC_InMSOL) ! inmsol.
       Call ptravel_locations (latIZ + NatIZ - 1,4, HC_InDIV) ! indiv.
       Call ptravel_locations (latIZ + NatIZ - 1,5, HC_InTrap) ! intrap.

       HC_LaunchDat (LatIZ + NatIZ - 1,3) = REAL (Reflection_Counter)
       HC_LaunchDat (LatIZ + NatIZ - 1,2) =  HC_Launchdat (IProd,2)

       ! Note:  Already have changed velocity via the state change to C+.
       Call pvins (LatIZ + NatIZ - 1,Current_Velocity)

       ! NOTE:  David Elder,    1995  AUG 3
       !
       ! There were too many new variables being required
       ! to keep track of the various ionization data under
       ! all the varied circumstances - so the information
       ! was centralized into one array which contains
       ! all of the information - there are some other
       ! variables which accumulate the data for
       ! compatibility with the older code - but all of
       ! the new code uses the new setup.
       !
       ! Here is a summary of the array and its contents
       !
       ! IONIZDAT(2     ,2         ,2    ,2      ,5)
       !          isol  Launch_Reg   ifp   irflct  quant
       !
       ! isol = 1 for SOL information
       !      = 2 for MAIN plasma information
       !
       ! Launch_Target    = 1 for target 1	ik > nks(ir)/2
       !                  = 2 for target 2	ik =< nks(i2)/2
       !
       ! ifp  = 1 for a neutral resulting from a regular
       ! 	 launch
       !      = 2 for a neutral resulting from a Far Periphery
       ! 	 relaunch
       !
       ! irflct = 1 for a neutral that has NOT been reflected
       !       = 2 for a neutral that has been reflected
       !
       ! quant= 1 total weight of neutrals
       !      = 2 weighted R coordinate
       !      = 3 weighted Z coordinate
       !      = 4 weighted K value
       !      = 5 weighted S value
       Launch_Target = 0

       If (Launch_Reg .eq. HC_Launch_Reg_Target_1_Dist .or. &
            &   Launch_Reg .eq. HC_Launch_Reg_Target_1_Pin .or. &
            &   Launch_Reg .eq. HC_Launch_Reg_Sput_Target_1 .or. &
            &   Launch_Reg .eq. HC_Launch_Reg_Refl_Target_1) Then
          Launch_Target = 1
       ElseIf (Launch_Reg .eq. HC_Launch_Reg_Target_2_Dist .or. &
            &       Launch_Reg .eq. HC_Launch_Reg_Target_2_Pin .or. &
            &       Launch_Reg .eq. HC_Launch_Reg_Sput_Target_2 .or. &
            &       Launch_Reg .eq. HC_Launch_Reg_Refl_Target_2) Then
          Launch_Target = 2
       End If

       If (Launch_Target .ne. 0) Then
          If (Current_Ring .ge.  Inner_SOL_Ring) Then
             isol = 1
          Else
             isol = 2
          End If

          If ( Far_Periphery_Recycle_Opt .eq. 1 .and. glaunchdat(IProd,2) .eq. 1.0) Then
             ifp = 2
          Else
             ifp = 1
          End If

          If (( Neutral_Reflection_Opt .eq. 1 .or.  Neutral_Reflection_Opt .eq. 2) .and. Reflection_Counter .gt. 0) Then
             irflct = 2
          Else
             irflct = 1
          End If

          HC_LIonizDat (isol,Launch_Target,ifp,irflct,1) =  HC_LIonizDat (isol,Launch_Target,ifp,irflct,1) + Sput_Weight
          HC_LIonizDat (isol,Launch_Target,ifp,irflct,2) =  HC_LIonizDat (isol,Launch_Target,ifp,irflct,2) + Sput_Weight * Current_R
          HC_LIonizDat (isol,Launch_Target,ifp,irflct,3) =  HC_LIonizDat (isol,Launch_Target,ifp,irflct,3) + Sput_Weight * Current_Z
          HC_LIonizDat (isol,Launch_Target,ifp,irflct,4) =  HC_LIonizDat (isol,Launch_Target,ifp,irflct,4) + Sput_Weight *  Local_K
          HC_LIonizDat (isol,Launch_Target,ifp,irflct,5) =  HC_LIonizDat (isol,Launch_Target,ifp,irflct,5) + Sput_Weight * MIN (&
               & Current_S,gksmaxs (Current_Ring) - Current_S)
       End If

       IFate = 23 ! HC Ion Reduced to C+.
       Return

    End If

    ! ION REMOVAL
    Random_Numbers_Used = Random_Numbers_Used + 1
    If (granv (Random_Numbers_Used) .lt. gkplos (Current_Cell,Current_Ring,Cur_HC_Spec)) Then
       HC_Num_Ions_Lost (Cur_HC_Spec,Launch_Reg) =  HC_Num_Ions_Lost (Cur_HC_Spec,Launch_Reg) + Sput_Weight
       HC_Min_Teq_Ions_Lost (Cur_HC_Spec,Launch_Reg) = MIN ( HC_Min_Teq_Ions_Lost (Cur_HC_Spec,Launch_Reg), Real(Total_Steps))
       HC_Max_Teq_Ions_Lost (Cur_HC_Spec,Launch_Reg) = MAX ( HC_Max_Teq_Ions_Lost (Cur_HC_Spec,Launch_Reg), Real(Total_Steps))
       HC_Tot_Teq_Ions_Lost (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Teq_Ions_Lost (Cur_HC_Spec,Launch_Reg) + Total_Steps*Sput_Weight
       Max_S_Ion_Removal (Cur_HC_Spec, Launch_Reg) = MAX (Max_S_Ion_Removal (Cur_HC_Spec, Launch_Reg), MIN (Current_S, gksmaxs (&
            &Current_Ring) - Current_S))

       IFate = 24 ! HC Ion Removed.
       Return
    End If

  End Subroutine HC_Trans_Inside_Ion

End Module HC_Inside_Ion
