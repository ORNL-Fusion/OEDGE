! -*-Mode:f90-*-
! HC_Follow.f90
! Purpose:
! Wrapping routine called upon chemical sputtering event for:
! -selection of HC species to launch for individual hydrocarbon molecules.
! -transport, evolution, and recording of sputtered HC.
! -sticking and reflection of HC impact on vessel/divertor.
! All steps are accomplished in general for HC's up to multi-C
! in size with a stack for current C molecules as fragments of
! original sputtering event.
!
! Adam McLean under Prof. Peter Stangeby
! Research Assistant: David Elder
! July, 2002
 
Module HC_Follow
 
  ! Required modules.
 
  ! Every good Fortran 90 program has...
  Implicit None
 
Contains
 
  Subroutine HC_Transport (		&
       &	IProd,				& 	! Number of particle tracked (in).
       &	LProd,				&	! First number of HC fragments to start with (in).
       &	NProd_input,			&	! Continue to this point (in).
       &	LatIZ,				& 	! First position in ionization array to be filled (in).
       &	NatIZ,				& 	! Number of initial hydrocarbons that become C+ ions (out).
       &	Seed,				& 	! Random seed (in/out).
       &	NRand,				& 	! Count of total random numbers used (in/out).
       &	NeuTime,			& 	! Time spent tracking neutrals (in/out).
       &	IonTime,			&	! Time spent tracking ions for intermediate HC molecules (in/out).
       &	SFail,				& 	! Number of failed launches V>Vmax at least 1000 (out).
       &	Status,				& 	! Generation of launch =1, primary, =2, secondary, etc, 10=total (in).
       &	NeutType,			& 	! Target (1,2,3), wall (4,5), 2Dneutral (6), Refl Ion (7), or  physical (1,4) or chemical (2,5) sputter (in).
       &	HC_Porm,			&	! PORM
       &	Random_Numbers_Used,		& 	! KK
       &	Max_Random_Values, 		& 	! KKLIM
       &	Storage)				! Data storage for distribution analysis.
 
    ! Required modules.
    Use ComHC ! Hydrocarbon-related variables set in DIVIMP input file.
    Use HC_Utilities ! Gain access to hydrocarbon-specific utilities.
    Use HC_Get ! External functions.
    Use HC_Put ! Place data.
    Use HC_Init_Lib_Data ! HC mass and charge subroutines.
    Use HC_Stack ! Load stored hydrocarbons upon completion of primary molecule.
    Use HC_Inside_Ion ! Operate on particle inside wall boundary.
    Use HC_Inside_Neutral ! Operate on particle inside wall boundary.
    Use HC_Outside_Ion ! Operate on particle outside wall boundary.
    Use HC_Outside_Neutral ! Operate on particle outside wall boundary.
    Use HC_Release ! Find initially sputtered hydrocarbon.
    Use HC_Init_DIV_Data
    Use HC_Init_DIV_Diag
    Use HC_FreeSpace_Transition
    Use HC_Ion_Transport ! Move an ion.
    Use HC_Neutral_Transport ! Move a neutral.

    use hc_kinetics ! hc kinetics - routines, options and types
    ! use hc_velocity_type - included via the hc_kinetics routine
    use error_handling
    
    use subgrid_options ! Detailed subgrid for more accurate location reporting of density profiles ...
    use subgrid

    ! Every good Fortran 90 program has...		
    Implicit None
 
    ! Define input/output variables
    Integer, Intent (In) :: IProd
    Integer, Intent (In) :: LPRod
    Integer, Intent (In) :: NProd_input
    Integer, Intent (In) :: LatIZ
    Integer, Intent (Out) :: NatIZ
    Double Precision, Intent (In) :: Seed
    Integer, Intent (InOut) :: NRand
    Real, Intent (InOut) :: NeuTime
    Real, Intent (InOut) :: IonTime
    Real, Intent (Out) :: SFail
    Integer, Intent (In) :: Status
    ! jdemod - neuttype changed to INOUT since it is changed somewhere below
    Integer, Intent (InOut) :: NeutType
    Real, Intent (InOut) :: HC_Porm ! PORM
    Integer, Intent (InOut) :: Random_Numbers_Used ! KK
    Integer, Intent (In) :: Max_Random_Values ! KKLIM
    Real, Dimension (:,:), Intent (InOut) :: Storage
 
    ! Define local variables that will be assigned values.
    Integer :: MTC_Counter ! MTCCNT
    Integer :: IFate ! IFATE
    Integer :: Reflection_Counter ! NRFCNT
    ! jdemod - because of the way the HC code handles particle recycling - make 
    !          NPROD a local variable with its initial value set to that passed in 
    !          from neut - also, need to add array bounds checking to the pxprods and 
    !          pyprods routines 
    integer :: nprod
    integer :: nprod_orig

    !
    ! HC velocity type variable to hold various pieces of the particle velocity
    !
    type (hc_velocity_type1) :: hc_v
    
    !
    Real :: Current_R ! R
    Real :: Current_Z ! Z
    Real :: Current_S ! S
    Real :: Current_Cross ! CROSS
    Integer :: Current_Cell ! IK
    Integer :: Current_Ring ! IR
    Integer :: Segment_Index ! ID
    Real :: Current_Theta ! THETA
    Integer :: Current_Dir_Indicator ! IS
    Real :: Current_Energy !
    Real :: Current_Velocity ! VIN
    Real :: Current_Velocity_In_R ! XVELF
    Real :: Current_Velocity_In_Z ! YVELF
    Real :: Current_Velocity_In_S ! VINS
    Real :: Current_Angle ! ANGLE
    Real :: Azimuthal_Angle
    Real :: Current_Tangent ! TANGNT
    Integer :: Cur_HC_Spec
    Integer :: Last_HC_Species
    Integer :: Local_Launch_HC_Species
    Real :: Last_R ! ROLD
    Real :: Last_Z ! ZOLD
    Real :: Last_S ! SOLD
    Real :: Last_Cross ! CROSSOLD
    Integer :: Last_Cell ! IKLAST
    Integer :: Last_Ring ! IRLAST
    Real :: Last_Angle
    Real :: Last_Velocity_In_R ! XVELFTMP
    Real :: Last_Velocity_In_Z ! YVELFTMP
    Real :: Last_Velocity ! VINTMP
    Real :: Last_Velocity_In_S !
    Real :: Last_Temperature !
    Real :: TStepHC ! TSTEPN
    Real :: Sput_Weight ! SPUTY
    Integer :: Launch_Reg ! M
    ! Note:  Time counts are DP so they are good over 16,777,216.
    Double Precision :: Neutral_Time_Step_Count ! CIST during NEUT.
    Double Precision :: Ion_Time_Step_Count
    Double Precision :: MTC_Time_Step_Count ! MTCCIST
    Double Precision :: Time_Steps_In_State ! CISTIZ
    Double Precision :: Eq_Total_Ion_Time_Steps ! Same as CIST, after NEUT is done.
    Double Precision :: Total_Time_Step_Count
    Double Precision :: Total_Time_Step_Count_Limit
    Integer :: Tries_Without_Change
    Integer :: Time_Bin ! IT
    Real :: Launch_Angle ! ANGLAN
    Real :: Launch_Velocity ! VIN
    !jdemod - due to variables with the same name at the global level - velocity_multiplier is made local here
    !         the global saved_velocity_mult is used to keep the value of velocity_multiplier visible at a global
    !         scope
    real :: velocity_multiplier
    !
    Real :: Beta ! BETA
    Real :: Psi ! PSI
    Real :: True_Tangent ! TANTRU
    Real :: Tangent_Launch ! TANLAN
    Real :: Random_Temp
    Integer :: Temp_Counter
    Integer :: Size_Of_Stack ! Size of the current hydrocarbon stack.
    Real :: HC_Temperature ! TEMN
    Real :: Sheath_Fraction ! sheath_fraction
    Integer :: Return_Code ! rc
    Real :: Time_Value
    Real :: Ran1
    Real :: Ran2
    Integer :: Charge
    Real :: S_Start ! SSTART
    Real :: Max_Velocity_Randoms ! RMAXS
    Real :: STmp ! STMP
    Real, Dimension (Number_HC_Species,Number_Regions) :: Min_Z ! XXX
    Real, Dimension (Number_HC_Species,Number_Regions) :: Max_S ! SSS
    Real, Dimension (Number_HC_Species,Number_Regions) :: Max_S_Ion_Removal ! SSSS
    Real, Dimension (Number_HC_Species,Number_Regions) :: Particle_Reach_State ! Particle reached state.
    Integer :: Prev_Cell
    Integer :: Prev_Ring
    Real :: Particle_Distance ! Distance travelled as this particle.
    Real :: State_R ! Initial R position for current state.
    Real :: State_Z ! Initial Z position for current state.
    Real :: Tot_Velocity ! Velocity sum over many timesteps.
    Real :: Kin_Energy_Added
    ! jdemod - the h_isotope_composition records data by array index 1=H, 2=D, 3=T - so 2H's, 2D's and a T would be 2,2,1
    !        - h_isotope_composition was declared inconsistently - by number_h_species in some routines and number_h_products elsewhere 
    !        - there were more number_h_products declarations - however, looking at the usage is should have been declared number_h_species
    !        - this has been fixed throughout the code
    Integer, Dimension (Number_H_Species) :: H_Isotope_Composition ! Number of H/D/T isotopes in the HC.
    Integer, Dimension (Number_H_Species) :: Last_H_Isotope_Composition ! Previous number of H/D/T isotopes in the HC.
    !Integer, Dimension (Number_H_Products) :: H_Isotope_Composition ! Number of H/D/T isotopes in the HC.
    !Integer, Dimension (Number_H_Products) :: Last_H_Isotope_Composition ! Previous number of H/D/T isotopes in the HC.
    Integer :: Particle_Ionized
    Integer :: M
 
    ! GA15 temporary variables.
    Integer :: GA15_IndWork (2, MaxPts)
    Real :: GA15_Work (4 * MaxPts)
    Real :: TDum (MaxPts)
    Real :: XDum (MaxPts)
    Real :: YDum (MaxPts)
    Real :: All_Wall_R_Points (MaxPts)
    Real :: All_Wall_Z_Points (MaxPts)
    Real :: Wall_Check_Result ! RESULT
    Logical :: In_Main ! inmain
 
    Double Precision :: Double_Sput_Weight ! DSPUTY
 
    ! External function declarations.
    Integer, External :: IPOS ! IPOS(R,RS,NRS) Finds nearest higher value in RS array (length NRS) to given R.
    Real, External :: CELLWIDTH


    ! jdemod - initialization - set values of NPROD
    nprod = nprod_input
    nprod_orig = nprod_input
    
    ! jdemod - initialize kin_energy_added
    kin_energy_added = 0.0

    ! Call global data setup routines.
    ! jdemod - implim does not need to be reset for each particle - in addition - calling this 
    !          routine re-initializes all of the debug variables - turning them off - debugging 
    !          can not work with this call at this point
    !Call Initialize_Misc_Data_Table () ! Needed here to set IMPLIM at each new molecule follow.


    !if (debug_hc) &
    !     & write(0,'(a,i6)') 'Particle:',iprod


    ! jdemod - custom debugging
    !
    !if (iprod.eq.90.or.iprod.eq.34) then
    !   debug_HC_neutral=.true.
    !   debug_HC_ion=.true.
    !   hc_coord_print_option = 1
    !   hc_evolve_print_option = 1
    !   print_debug_x_cstep_neutral =1 
    !   print_debug_x_cstep_ion    =1 
    !else
    !   debug_HC_neutral=.false.
    !   debug_HC_ion=.false.
    !   hc_coord_print_option = 0
    !   hc_evolve_print_option = 0
    !   print_debug_x_cstep_neutral =0
    !   print_debug_x_cstep_ion    =0 
    !endif
 
    ! Set initial fate.
    IFate = 0
 
    ! Set initial time.
    Time_Value = 0.0
 
    ! Debug test.
    TStepHC =  Print_Debug_x_CStep_Ion ! TSTEPL=CSTEPL.  When to print debugging information.
 
    ! Zero momentum transfer collision (MTC) counter.
    MTC_Counter = 0
 
    ! Reflection counter.
    Reflection_Counter = 0
 
    ! Reflection-has-occurred indicator reset to false for each particle launched.
    HC_Hit_Vessel = .False.
 
    ! Initialize reached state variable.
    Particle_Reach_State = 0.0 ! Increment to 1.0 once particle reaches state.
 
    ! Set distance travelled by particle.
    Particle_Distance = 0.0
 
    ! jdemod - NOTE: Determines ONLY the HC species - many arguments NOT required
    ! Determine hydrocarbon species to be released.
    ! NOTE: Must find impinging energy and angle off normal of particle that began the chemical sputtering process here!!!
    Call Release_HC (Current_R,Current_Z, Launch_Wall_Index,0.0,0.0,Local_Launch_HC_Species) ! HC_Follow is iterated for neutral impacts.
 
    ! Save the originally launched hydrocarbon species.
    Cur_HC_Spec = Local_Launch_HC_Species
 
    ! Determine H isotope content given plasma background.  Note:  This routine should only be called once per HC launch.
    Call Get_HC_Isotope_Composition (Cur_HC_Spec,H_Isotope_Composition,Seed,NRand)
 
    ! Save composition to last value.
    Last_H_Isotope_Composition = H_Isotope_Composition
 
    ! Find launch location.
    If (Get_HC_Charge (Cur_HC_Spec) .eq. 0) Then
       ! Similar to CNEUTA = 0 in DIVIMP, Find current position given Launch_Option.  All 'current' and 'launch' variables are returned.
       Call HC_Starting_Position (IProd,LProd,NeutType,Current_R,Current_Z,Current_Cell,Current_Ring,Segment_Index,&
                                & Current_Dir_Indicator,Launch_Reg)
       ! Note, only required to be called once per particle as afterwards, impact position will always be known.
       ! Set Particle_Ionized to 0 to begin as neutral.
       Particle_Ionized = 0
    Else
       ! Similar to CNEUTA = 1 in DIVIMP, Launch position determined by INJECTION option.
       Call HC_Injection_Position (Seed,NRand,Current_R,Current_Z,Current_S,Current_Cross, &
            & NeutType,STmp,Current_Cell,Current_Ring,Segment_Index,Current_Dir_Indicator,Launch_Reg)
       ! Set Particle_Ionized to 1 to begin as neutral.
       Particle_Ionized = 1
    End If
    !Write (Output_Unit_Scratch,*) "R,Z,S,Cross",Current_R,Current_Z,Current_S,Current_Cross,Current_Cell,Current_Ring,Cur_HC_Spec,Get_HC_Charge (Cur_HC_Spec)
    !Write (0,*) "R,Z,S,Cross",Current_R,Current_Z,Current_S,Current_Cross,Current_Cell,Current_Ring,Cur_HC_Spec,Get_HC_Charge (Cur_HC_Spec)
 
    ! Remember original launch position.
    Launch_R = Current_R
    Launch_Z = Current_Z
    Launch_HC_Species = Cur_HC_Spec
 
    ! Initialize location data.- seems to save only KKS and KSMAXS values and some drift velocity limits?
    Call Initialize_Cell_Geom_Data (Current_Cell,Current_Ring)
 
    !
    ! jdemod - these should all be sput_weight - not 1.0 - a half-weight particle should only count as 1/2 a particle reaching a given state from the 
    !          perspective of generating averages and calculating code data - code moved to after sput_weight is calculated
    !
    !Particle_Reach_State (Cur_HC_Spec,Launch_Reg) = 1.0
    !HC_Num_Orig_Reach_State (Cur_HC_Spec,Launch_Reg) =  HC_Num_Orig_Reach_State (Cur_HC_Spec,Launch_Reg) + 1.0
    !HC_Tot_Reach_State (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Reach_State (Cur_HC_Spec,Launch_Reg) + 1.0
    !

    State_R = Current_R
    State_Z = Current_Z
 

    !
    ! jdemod - record start of particle walks
    !
    !write(6,'(a,3i5,5g16.8)') 'HC_WALKS:',current_cell,current_ring,Cur_HC_Spec, current_r,current_z
 
 
    If ( HC_Walk_Count .lt. Max_Number_Walks) Then
       ! write(0,*) 'WALKS:', HC_Walk_Count, Max_Number_Walks
       HC_Walks ( HC_Walk_Count, 1) = 100.0 *  Maximum_R
       HC_Walks ( HC_Walk_Count, 2) = 100.0 *  Maximum_Z
       HC_Walk_Count =  HC_Walk_Count + 1
       If ( HC_Walk_Count .lt. Max_Number_Walks) Then
          HC_Walks ( HC_Walk_Count, 1) = Current_R
          HC_Walks ( HC_Walk_Count, 2) = Current_Z
          HC_Walks ( HC_Walk_Count+1, 1) =  Calc_Hi
          HC_Walks ( HC_Walk_Count+1, 2) =  Calc_Hi
          HC_Walk_Count =  HC_Walk_Count + 1
       End If
    End If
 
 



    ! jdemod
    ! This call to TAUIN2 is a left over - since the charge is always 0 at this point - this code does nothing 
    ! In addition - Adam's routine called modify_taus simply recalculates the TAUS value from 
    ! scratch and doesn't use any code from TAU. 
    !Call TAUIN2 (Get_HC_Charge (Cur_HC_Spec))
 
    ! Assign last value to hydrocarbon species.
    Last_HC_Species = Cur_HC_Spec
    Last_H_Isotope_Composition = H_Isotope_Composition		
 
    Sput_Weight = Calc_Sputy (Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition))
    Double_Sput_Weight = Sput_Weight



    !
    ! jdemod - these should all be sput_weight - not 1.0 - a half-weight particle should only count as 1/2 a particle reaching a given state from the 
    !          perspective of generating averages and calculating code data.
    !
    !Particle_Reach_State (Cur_HC_Spec,Launch_Reg) = 1.0
    !HC_Num_Orig_Reach_State (Cur_HC_Spec,Launch_Reg) =  HC_Num_Orig_Reach_State (Cur_HC_Spec,Launch_Reg) + 1.0
    !HC_Tot_Reach_State (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Reach_State (Cur_HC_Spec,Launch_Reg) + 1.0
    !
    Particle_Reach_State (Cur_HC_Spec,Launch_Reg) = sput_weight
    HC_Num_Orig_Reach_State (Cur_HC_Spec,Launch_Reg) =  HC_Num_Orig_Reach_State (Cur_HC_Spec,Launch_Reg) + sput_weight
    HC_Tot_Reach_State (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Reach_State (Cur_HC_Spec,Launch_Reg) + sput_weight


 
    ! Set S, Cross, Theta, S_Start, value for non-orthogonal transport if beginning as an ion.
    If (Get_HC_Charge (Cur_HC_Spec) .ne. 0) Then
       Call Neut_Ion_FreeSpace_Location (0,1,Current_R,Current_Z,Current_Cell,Current_Ring,STmp,Current_S,Current_Cross,&
                                       & Current_Theta,S_Start)
    Else
       !
       ! jdemod - calling getscross_approx is not relevant for a neutral - in addition, the subsequent code here overwrites the 
       !          return values from the routine making it useless. Remove call.
       !          Also, since there has been no check to see if the particle is even on the grid or has a valid ik,ir index - 
       !          calling getscross_approx can only lead to problems here. 
       ! Note - all of the code here for launching ions might be deleted - however, if launching ionized HC fragments may be of 
       !        interest then the code should support ion launch - I think the code will need significant modifications to do this
       !        properly though.
       !
       !Call getscross_approx (Current_R,Current_Z,Current_S,Current_Cross,Current_Cell,Current_Ring)
       !
       ! jdemod - what is the following doing? I don't think these are correct for a neutral launch at all - S and CROSS only apply to ion transport.
       !        - consider it as initialization at this point for each new particle and remove later if needed. 
       !
 
       If (Launch_Reg .eq. 1) Then
          Current_S = gksmaxs (Current_Ring)
       ElseIf (Launch_Reg .eq. 2) Then
          Current_S = 0.0
       End If
       Current_Theta = 0.0
       Current_Cross = 0.0
       S_Start = Current_S ! SSTART
    End If

 
    ! Initialize timing counters.
    Tries_Without_Change = 0
    Neutral_Time_Step_Count = 0.0
    Ion_Time_Step_Count = 0.0
    MTC_Time_Step_Count = 0.0
    Time_Steps_In_State = 0.0
    Eq_Total_Ion_Time_Steps = 0.0
    Total_Time_Step_Count = 0.0
    Total_Time_Step_Count_Limit = 10000000
 
    ! Particle initialization.
    HC_Has_Leaked = .False.
    HC_Has_Leaked_Core = .False.
    HC_Leak_Particles (Launch_Reg) = 1 ! CLEAKP
    Min_Z = 0.0 ! XXX
    Max_S = 0.0 ! SSS
    Max_S_Ion_Removal = 0.0 ! SSSS
 
 
    !
    ! jdemod - time_bin can only be useful for time dependent analysis and density arrays must be available to support them - is that true?
    !
 
    ! Beginning time bin.	
    Time_Bin = IPOS (Neutral_Time_Step_Count, gctimes (1,0),  Number_Time_Steps + 1)
 
    ! Save erosion statistics
    ! Add to erosion.
    If (hc_launch_location .eq. 0 .or. hc_launch_location .eq. 3) Then
       ! Add to total removal from target.
       HC_Erosion (Segment_Index,3,Cur_HC_Spec) =  HC_Erosion (Segment_Index,3,Cur_HC_Spec) + Sput_Weight ! NEROS.
       ! Note:  NeutType = 3 indicates a target self-sputter launch which
       ! will only occur from the final call to LAUNCH from DIV.
       If (NeutType .ne. 3) Then
          ! Add to primary removal from target.
          HC_Erosion (Segment_Index,2,Cur_HC_Spec) =  HC_Erosion (Segment_Index,2,Cur_HC_Spec) + Sput_Weight ! NEROS.
       End If
       ! Add to removal from wall segment.
       If (gwallindex (Segment_Index) .ne. 0.0) Then
          HC_WallsE (gwallindex (Segment_Index),Cur_HC_Spec) =  HC_WallsE (gwallindex (Segment_Index),Cur_HC_Spec) + Sput_Weight
       Else
          Write (Output_Unit_HC_Alert,'(A,5I5,3(1X,G12.5))') "Error in HC_Follow: HC Wallse:Target?",Segment_Index,Current_Cell,&
                                                       & Current_Ring,IProd,LProd,gwallindex (Segment_Index)
          HC_WallsE ( Max_Points + 1,Cur_HC_Spec) =  HC_WallsE ( Max_Points + 1,Cur_HC_Spec) + Sput_Weight
       End If
    ElseIf (hc_launch_location .eq. 2 .or. hc_launch_location .eq. 4 .or. hc_launch_location .eq. 6) Then
       ! Add to primary removal from wall segment if it is also a target segment.
       ! Note: WALLPT is an array of type REAL.
       If (INT (gwallpt (Segment_Index,18)) .gt. 0) Then
          ! Note:  NeutType = 3 indicates a target self-sputter launch which
          ! will only occur from the final call to LAUNCH from DIV.
          If (NeutType .ne. 3) Then ! Note:  This should never occur.  Will need to change this value if we plan to support wall-launched self-sputtered particles.
             ! Add to primary removal from target.
             HC_Erosion (INT (gwallpt (Segment_Index,18)),2,Cur_HC_Spec) =  HC_Erosion (INT (gwallpt (Segment_Index,18)),2,&
                                                                          & Cur_HC_Spec) + Sput_Weight ! NEROS.
          End If
          ! Add to total removal from target.
          HC_Erosion (INT (gwallpt (Segment_Index,18)),3,Cur_HC_Spec) =  HC_Erosion (INT (gwallpt (Segment_Index,18)),3,&
                                                                      &  Cur_HC_Spec) + Sput_Weight ! NEROS.
       End If
       ! Add to removal from wall segment.
       If (Segment_Index .lt. 1 .or. Segment_Index .gt.  Num_Wall_Points) Then
          Write (Output_Unit_HC_Alert,'(A,6I5)') "Error in HC_Follow: HC Wallse:Target?:",Segment_Index,Current_Cell,Current_Ring,&
                                                & IProd,LProd,gwallindex (Segment_Index)
          HC_WallsE ( Max_Points + 1,Cur_HC_Spec) =  HC_WallsE ( Max_Points + 1,Cur_HC_Spec) + Sput_Weight
       Else
          HC_WallsE (Segment_Index,Cur_HC_Spec) =  HC_WallsE (Segment_Index,Cur_HC_Spec) + Sput_Weight
       End If
 
    End If


    ! Use the nrs+1 element to record the total source from
    ! all elements. The portion that is ionized is recorded
    ! in the rest of the array.   Do not do for ion injection
    ! Note: neuttype keeps the different indexing for target and wall sources separated
    If (Get_HC_Charge (Cur_HC_Spec) .eq. 0) Then
       If ( Starting_Index .ge. 1 .and.  Starting_Index .le.  Num_Wall_Points) Then
          If (NeutType .gt. 0) Then				
             HC_WTSource ( Starting_Index, Num_Upper_Rings + 1,1,NeutType) = &
                  &  HC_WTSource ( Starting_Index, Num_Upper_Rings + 1,1,NeutType) + Sput_Weight
          EndIf
          !Write (Output_Unit_Scratch,*) "SOURCING1:", HC_WTSource ( Starting_Index, Num_Upper_Rings + 1,1,NeutType),NeutType,Cur_HC_Spec
       Else
          Write (Output_Unit_HC_Alert,*) "Error in HC_Follow: Starting index is not within applicable range: ", Starting_Index, &
                                       & Num_Wall_Points
          Write (Output_Unit_HC_Alert,*) "Program stopping."
          Stop
       End If
    End If
 
    ! Load local copy of all wall points in RW and ZW for GA15.
    Do Temp_Counter = 1, MaxPts
       All_Wall_R_Points (Temp_Counter) = grw (Temp_Counter)
       All_Wall_Z_Points (Temp_Counter) = gzw (Temp_Counter)
    End Do
    !write (0,*) "ALLR",All_Wall_R_Points
    !write (0,*) "ALLZ",All_Wall_Z_Points
 
    HC_CflRin = .True.
    HC_CflRex = .True.
 


    ! Calculate angle of hydrcarbon release.
    Call HC_Launch_Angle (IProd,Beta,Psi,Launch_Angle,Seed,NRand) ! Called at each instance of reflection/sputtering as well.
    Azimuthal_Angle = Psi
 
    ! Calculate normal angle to target.
    Call Target_Normal (Segment_Index,Current_Dir_Indicator,Launch_Angle,True_Tangent,Tangent_Launch, &
                   &Current_Angle,Current_Tangent,&
                   & Status) ! Called at each instance of reflection as well.
    !write (Output_Unit_Scratch,*) "angle info:",Segment_Index,Current_Dir_Indicator,Launch_Angle,True_Tangent,Tangent_Launch,Current_Angle,Current_Tangent
    ! Find proper new angle.
    Current_Angle = Current_Angle + Current_Tangent ! Launch angle distributed around SIN (or other) added to tangent angle off starting wall/target segment.
    ! Artificial launch angle.
    !Current_Angle = 1.57




    ! jdemod - comment - the following code does nothing since max_velocity_randoms=1.0 is assigned below
    ! Find maximum random value for launch.
    If (HC_Launch_Location .eq. 0 .or. HC_Launch_Location .eq. 3 .or. HC_Launch_Location .eq. 6) Then
       ! Target launch.
       Max_Velocity_Randoms = gkrmax (Segment_Index) ! Here Segment_Index is the target index.
    ElseIf (HC_Launch_Location .eq. 2 .or. HC_Launch_Location .eq. 4) Then
       ! Wall launch.
       Max_Velocity_Randoms = gkrmaxw (Segment_Index) ! Here Segment_Index is the wall index.
    End If


    ! Ensure that there is no lower energy cuttoff for chemical sputtering.
    Max_Velocity_Randoms = 1.0




 
    ! Calculate velocity of hydrocarbon release.
    If (Get_HC_Charge (Cur_HC_Spec) .eq. 0) Then
       Call HC_Launch_Velocity	(IProd,LProd,Current_Cell,Current_Ring,Segment_Index,Cur_HC_Spec, &
            & H_Isotope_Composition,Sput_Weight, &
            & Max_Velocity_Randoms,Beta,Psi,Seed,NRand,Launch_Reg,HC_Temperature, Velocity_Multiplier,Launch_Velocity,IFate)
       ! Call at each instance of reflection as well.
       !HC_Temperature = 0.5
    Else
       ! Use ion injection parameters for launch.
       Call HC_Injection_Velocity (Current_Cell,Current_Ring,Cur_HC_Spec,H_Isotope_Composition,Sput_Weight,Seed,NRand,HC_Porm,&
                                   &HC_Temperature,Launch_Velocity)
       !HC_Temperature = 0.5
       !write (0,*) "IPROD:",IPROD
       !If (IPROD .gt. 500) Then
       !Launch_Velocity = 50000.0 ! REMOVE ME AMMOD
       !else
       !Launch_Velocity = -50000.0 ! REMOVE ME AMMOD
       !endif
       !Launch_Velocity = 0.0 ! REMOVE ME AMMOD
 
 
       ! Initialize state data and modify taus for starting state.
       Call Initialize_State_Prop_Data (Current_Cell,Current_Ring,Get_HC_Charge (Cur_HC_Spec)) ! Sets LFTS, LFPS, LFSS, LLLFPS, LTOLDS
       Call Modify_Taus (Cur_HC_Spec,H_Isotope_Composition,Sput_Weight,Current_Cell,Current_Ring)
       ! jdemod
       !Call Modify_Taus (Cur_HC_Spec,H_Isotope_Composition,Sput_Weight,HC_Temperature,Current_Cell,Current_Ring,Current_S,Seed,&
       !                & Random_Numbers_Used,NRand)
    End If


    ! Sum position statistics for initial particle launch.
    HC_Total_R_Prod_Positions(Cur_HC_Spec,Launch_Reg)=HC_Total_R_Prod_Positions(Cur_HC_Spec,Launch_Reg)+Current_R*Sput_Weight
    HC_Total_Z_Prod_Positions(Cur_HC_Spec,Launch_Reg)=HC_Total_Z_Prod_Positions(Cur_HC_Spec,Launch_Reg)+Current_Z*Sput_Weight
    HC_Total_S_Prod_Positions(Cur_HC_Spec,Launch_Reg)=HC_Total_S_Prod_Positions(Cur_HC_Spec,Launch_Reg)+Current_S*Sput_Weight
    HC_Total_Cross_Prod_Positions (Cur_HC_Spec,Launch_Reg) =  HC_Total_Cross_Prod_Positions (Cur_HC_Spec,Launch_Reg) + &
                                                            & Current_Cross * Sput_Weight
    HC_Total_Prod_Angles (Cur_HC_Spec,Launch_Reg) =  HC_Total_Prod_Angles (Cur_HC_Spec,Launch_Reg) + (Launch_Angle + &
                                                   & Tangent_Launch) * Sput_Weight
    HC_Sum_Fragments_Launched (Cur_HC_Spec,Launch_Reg) =  HC_Sum_Fragments_Launched (Cur_HC_Spec,Launch_Reg) + Sput_Weight
 

 
    ! Record creation of new species at the end of the previous timestep.
    If (Get_HC_Charge (Cur_HC_Spec) .eq. 0) Then
       Call Record_Evolve_Data (Cur_HC_Spec,Launch_Reg,Current_R,Current_Z,Current_Cell,Current_Ring,Current_S,Current_Cross,&
                              & Current_Angle,Launch_Velocity,HC_Temperature,Sput_Weight)
    Else
       Call Record_Evolve_Data (Cur_HC_Spec,Launch_Reg,Current_R,Current_Z,Current_Cell,Current_Ring,Current_S,Current_Cross,&
                              & Current_Angle,Launch_Velocity,HC_Temperature,Sput_Weight)
    End If

 
    ! Print status info for new particle.
    !Write (0,9292) "HC launched.  Species:",Cur_HC_Spec,"NeutType:",NeutType,"Particle:",IProd,"Launch:",hc_launch_location,"Vel/Ang:",hc_launch_angle_velocity,"R:",Current_R,"Z:",Current_Z,"Ang:",Current_Angle,"Vel:",Launch_Velocity,"Temp:",HC_Temperature
    !9292 Format (1X,A,1X,I2,1X,A,1X,I1,1X,A,1X,I7,1X,A,1X,I1,1X,A,1X,I1,1X,A,1X,F8.3,1X,A,1X,F8.3,1X,A,1X,F8.3,1X,A,1X,F10.3,1X,A,1X,F8.3)
 
    ! Write particle number to output files if selected.
    If (hc_coord_print_option .eq. 1) Then
       Write (Output_Unit_Location,*) ""
       Write(Output_Unit_Location,9300)"New particle launched.NeutType:",NeutType,"Particle:",IProd,"Launch:",hc_launch_location,&
                                & "Vel/Ang:",hc_launch_angle_velocity,"R:",Current_R,"Z:",Current_Z,"Ang:",Current_Angle,&
                                & "Vel:",Launch_Velocity,"Temp:",HC_Temperature
9300   Format (1X,A,1X,I1,1X,A,1X,I7,1X,A,1X,I1,1X,A,1X,I1,1X,A,1X,F8.3,1X,A,1X,F8.3,1X,A,1X,F8.3,1X,A,1X,F10.3,1X,A,1X,F8.3)
 
       write (Output_Unit_Location,9301) "Step Time","Neut","Ion","HC","Iso-H/D/T","Cell","Ring","R","Z","S","Cross","Angle","Vel",&
            & "Temp","Density", &
            & "E Temp","I Temp","E Field","LE Field","Alphs","Betas","Col Prob","Para In","Stop In","Heat In","SPara","VPara","R &
            &Const","R Gauss"
9301   Format (40A10)
    End If
    If (hc_evolve_print_option .eq. 1) Then
       Write (Output_Unit_Evolve,*) ""
       Write (Output_Unit_Evolve,9302) "New particle launched.  NeutType:",NeutType,"Particle:",IProd,"Launch:",hc_launch_location,&
           & "Vel/Ang:",hc_launch_angle_velocity,"R:",Current_R,"Z:",Current_Z,"Ang:",Current_Angle,"Vel:",Launch_Velocity,"Temp:",&
           & HC_Temperature
9302   Format (1X,A,1X,I1,1X,A,1X,I7,1X,A,1X,I1,1X,A,1X,I1,1X,A,1X,F8.3,1X,A,1X,F8.3,1X,A,1X,F8.3,1X,A,1X,F10.3,1X,A,1X,F8.3)
 
       Write (Output_Unit_Evolve,9303) "Event Time","Neut","Ion","Init HC","H/D/T B","Dest HC","H/D/T A","Cell","Ring","R Before",&
            & "R After","Z Before", &
            & "Z After","S Before","S After","Cross B","Cross A","Angle B","Angle A","Vel B","Vel A","Temp B","Temp A"
9303   Format (40A10)
    End If
 
    !storage (IProd,1) = Current_Angle
    !storage (IProd,2) = Launch_Velocity
    if (Launch_Velocity .ge. 5000.0 .or. Launch_Velocity .le. 0.0) Then
       write (Output_Unit_HC_Alert,*) "Check launch velocity:",Launch_Velocity
    end If
    !Write (Output_Unit_Scratch,*) "Returning"
    !Return	
 

    If (IFate .eq. 0) Then ! Proceed as normal

       !
       ! jdemod - Setup the structure that will contain the detail hc velocity information based on the launch data
       ! Assign the initial velocities for the HC fragment
       ! NOTE: UPDATE: The difference in treatment between ion and neutral launch velocity was a bug and is now fixed - launch
       !       velocity is now in m/s for both ions and neutrals. X launch velocity for neutrals is not scaled by the time step while launch_velocity for ions already has the 
       !       time step factored in. X
       ! NOTE: Current_angle lies in the R,Z plane while the azimuthal angle(psi) gives the out of plane 
       !       component. It is measured from psi = 0 lying in the R,Z plane anti-clockwise (looking down) from the current_angle vector. 
       call assign_hc_velocity(hc_v, launch_velocity, velocity_multiplier, current_angle, psi, hc_temperature, &
                           & current_velocity, current_velocity_in_s, current_velocity_in_r, current_velocity_in_z, &
                           & neutral_time_step, ion_time_step, get_hc_charge(cur_hc_spec), find_hc_mass(cur_hc_spec))


 
       ! Sum velocity statistics for initial particle launch.
       HC_Total_Prod_Vels_No_VMult (Cur_HC_Spec,Launch_Reg) =  HC_Total_Prod_Vels_No_VMult (Cur_HC_Spec,Launch_Reg) + &
                                                              &Launch_Velocity /  Velocity_Multiplier * Sput_Weight
       HC_Max_Prod_Vel_No_VMult (Cur_HC_Spec,Launch_Reg) = MAX ( HC_Max_Prod_Vel_No_VMult (Cur_HC_Spec,Launch_Reg), &
                                                         &Launch_Velocity /  Velocity_Multiplier)
       !			 HC_Total_Vel_Ang_Mults (Cur_HC_Spec,Launch_Reg) =  HC_Total_Vel_Ang_Mults (Cur_HC_Spec,Launch_Reg) +  Velocity_Multiplier * Sput_Weight
       !			 HC_Max_Vel_Ang_Mults (Cur_HC_Spec,Launch_Reg) = MAX ( HC_Max_Vel_Ang_Mults (Cur_HC_Spec,Launch_Reg),  Velocity_Multiplier)
       HC_Total_Prod_Velocities(Cur_HC_Spec,Launch_Reg)=HC_Total_Prod_Velocities(Cur_HC_Spec,Launch_Reg)+Launch_Velocity*Sput_Weight
       HC_Max_Prod_Velocities (Cur_HC_Spec,Launch_Reg) = MAX ( HC_Max_Prod_Velocities (Cur_HC_Spec,Launch_Reg), Launch_Velocity)



       HC_Tot_Prod_Temperatures(Cur_HC_Spec,Launch_Reg)=HC_Tot_Prod_Temperatures(Cur_HC_Spec,Launch_Reg)+HC_Temperature*Sput_Weight


       HC_Time_At_Production (IPROD,Cur_HC_Spec) = 0.0
       HC_Energy_At_Production (IPROD,Cur_HC_Spec) = HC_Temperature

       if (debug_hc) &
             &    write(6,'(a,2i5,1x,a,10(1x,g12.5))') 'START:',iprod,cur_hc_spec,hc_ident_species(cur_hc_spec), &
             &                                time_value,hc_temperature,kin_energy_added


 
       ! Record debugging information.			
       If ( Debug_HC_Neutral) Then
          ! jdemod - convert to common printout
          Write (Output_Unit_Scratch,9004) IProd, Neutral_Time_Step_Count , Ion_Time_Step_Count, Current_Cell,Current_Ring,&
                 & Cur_HC_Spec,Current_R,Current_Z,Current_S, Local_K, &
                 & Current_Theta,gksmaxs (Current_Ring),Current_Velocity,HC_Temperature, Local_SPara,Current_Cross,&
                 & Sput_Weight,Time_Bin,'Neutral HC Launch'

          !Write (Output_Unit_Scratch,9003) IProd, Neutral_Time_Step_Count, Current_Cell, Current_Ring,Current_R,Current_Z, Local_K,&
          !     & Current_Velocity, HC_Temperature, &
          !     & Sput_Weight,(Current_Angle) *  Rad_In_A_Deg,Time_Bin,'Neutral HC Launch'

       End If
 
       ! Check if hydrocarbon is going to strike target surface immediately.
       If ((hc_launch_location .ne. 1 .and. hc_launch_location .ne. 5) .and. &
            & ((Launch_Angle + Tangent_Launch .le. True_Tangent -  Pi_Value / 2.0) .or. &
            & (Launch_Angle + Tangent_Launch .ge. True_Tangent +  Pi_Value / 2.0))) Then
          Write (Output_Unit_HC_Alert,*) "Warning in HC_Follow: Out here:",Launch_Angle,Tangent_Launch,True_Tangent
          HC_Num_Striking_Target (Cur_HC_Spec,Launch_Reg) =  HC_Num_Striking_Target (Cur_HC_Spec,Launch_Reg) + Sput_Weight
 
          If (MTC_Counter .gt. 0) Then
             HC_MTC_Striking_Target (Cur_HC_Spec,Launch_Reg) =  HC_MTC_Striking_Target (Cur_HC_Spec,Launch_Reg) + Sput_Weight
          EndIf
          IFate = 12
 
       EndIf
 

       If (IFate .eq. 0) Then ! Continue as normal.
 
          If (IProd .eq. (IProd / 100) * 100) Then
             !Write (Output_Unit_Scratch,'(a,i6,3i4,5(1x,g13.5))') 'Launch:',IProd,Segment_Index,Current_Cell,Current_Ring,Current_R,Current_Z,Current_Velocity,Current_Angle
          EndIf
 
          ! Prepare GA15 boundary-checking routine.
          Call GA15A ( Num_Boundary_Points,1,GA15_Work,4* Max_Points,GA15_IndWork, Max_Points,All_Wall_R_Points,All_Wall_Z_Points,&
                     &TDum,XDum,YDum,6)
 
          ! Iterate around motion for each timestep until ionization occurs.
          ! Generation new set of randoms when old lot are used up.
 
          ! Set initial cell occupied.
          Last_Cell = Current_Cell
          Last_Ring = Current_Ring
          Prev_Cell = Current_Cell
          Prev_Ring = Current_Ring
          !write (0,*) "INIT TEMP:",HC_Temperature
          Last_Temperature = HC_Temperature
 
          ! Check if initial position is inside wall.
          Call GA15B (Current_R,Current_Z,Wall_Check_Result, Num_Boundary_Points,1,GA15_Work,4* Max_Points,GA15_IndWork, &
                     &Max_Points,All_Wall_R_Points,All_Wall_Z_Points,TDum,XDum,YDum,6)
 
          If (Wall_Check_Result .lt. 0.0) Then
             Write (Output_Unit_HC_Alert,'(a,3(1x,g16.8))') &
                  & 'ERROR: Particle initial position is'// &
                  & ' apparently outside wall:',Current_R,Current_Z,Wall_Check_Result
 
             Write (Output_Unit_HC_Alert,'(a,5(1x,i6),5(1x,g13.5))') &
                  & 'ERROR LAUNCH:',IProd, Launch_Wall_Index,Segment_Index,Current_Cell, &
                  & Current_Ring,Current_R,Current_Z,Current_Velocity,Current_Angle
 
             !  Move particle by one time step - it should only be outside for numerical reasons.
             Current_R = Current_R + Current_Velocity_In_R
             Current_Z = Current_Z + Current_Velocity_In_Z
             Neutral_Time_Step_Count = Neutral_Time_Step_Count + 1.0
             Eq_Total_Ion_Time_Steps = Eq_Total_Ion_Time_Steps +  Neutral_Time_Step /  Ion_Time_Step
 
             ! Check again if initial position is inside wall
             Call GA15B (Current_R,Current_Z,Wall_Check_Result, Num_Boundary_Points,1,GA15_Work,4* Max_Points,GA15_IndWork, &
                        &Max_Points,All_Wall_R_Points,All_Wall_Z_Points,TDum,XDum,YDum,6)
 
             If (Wall_Check_Result .ge. 0.0) Then
                Write (Output_Unit_HC_Alert,'(a,6(1x,g16.8))') 'ERROR FOLLOWUP: OK : Wall_Check_Result,Current_R,Current_Z',&
                                                             &Wall_Check_Result,Current_R,Current_Z
                HC_Launch_Grid_Error_Moved_Okay (Cur_HC_Spec,Launch_Reg) =  HC_Launch_Grid_Error_Moved_Okay (Cur_HC_Spec,&
                                                                            & Launch_Reg) + Sput_Weight
             Else
                Write (Output_Unit_HC_Alert,'(a,6(1x,g16.8))') 'ERROR FOLLOWUP: STILL OUTSIDE : Wall_Check_Result,Current_R,&
&Current_Z', Wall_Check_Result,Current_R,Current_Z
                HC_Launch_Grid_Error_Moved_Out (Cur_HC_Spec,Launch_Reg) =  HC_Launch_Grid_Error_Moved_Out (Cur_HC_Spec,Launch_Reg) &
                                                                       &  + Sput_Weight
             End If
 
          End If
 
          ! Write out diagnostic information
          If (Status .le. 10 .and. &
               & (((NProd - LProd + 1) .lt. 1000) .or. &
               & ((NProd - LProd + 1) .lt. 10000 .and. (IProd / 10) * 10.0 .eq. IProd) .or. &
               & ((IProd / 100) * 100.0 .eq. IProd))   ) Then
             ! jdemod - write to the scratch unit instead of the alert unit
             Write (Output_Unit_scratch,'(a,2i6,4i4,6(1x,g12.5))') 	'NEUT-A:',IProd,LProd+IProd-1,&
                  & Current_Cell,Current_Ring,Segment_Index,Current_Dir_Indicator, &
                  & Current_R,Current_Z,Current_Velocity,Current_Velocity_In_R,Current_Velocity_In_Z,Neutral_Time_Step_Count
          End If
 
 
          !--------------------------------------------------------- Loop Start ------------------------------------------------------------------------------------------
 
 
          ! Start loop for multiple C species from hydrocarbons larger than CH4.
          Do
 
             ! Reset velocity sum for particle thus far to zero for each C molecule.
             Tot_Velocity = 0.0
 
             ! Transport, evolve, and record each hydrocarbon.
             Do

 
                ! Check if we need more random numbers.
                If (Random_Numbers_Used .gt. Max_Random_Values) Then
                   Call peranv (Seed, Random_Numbers_Used) ! Fill RANV array externally.
                   NRand = NRand + Random_Numbers_Used
                   Random_Numbers_Used = 0
                End If
 
                ! Adjust impurity velocity for change of cell if option is set.
                If (Get_HC_Charge (Cur_HC_Spec) .eq. 0) Then
 
                   ! Neutral hydrocarbon species.
                   !Write (Output_Unit_Scratch,*) "Before neutral step:",Current_R,Current_Z,Current_Cell,Current_Ring,Current_Velocity,HC_Temperature,Current_Angle
                   Call Neutral_Move (		&
                        &	Current_R,		&
                        &	Current_Z,		&
                        &	Current_Cell,		&
                        &	Current_Ring,		&
                        &	Current_Velocity,	&
                        &	Current_Velocity_In_R,	&
                        &	Current_Velocity_In_Z,	&
                        &	Current_Angle,		&
                        &	Cur_HC_Spec,	&
                        &	H_Isotope_Composition,	&
                        &	Last_R,			&
                        &	Last_Z,			&
                        &	Last_Cell,		&
                        &	Last_Ring,		&
                        &	Last_Velocity,		&
                        &	Last_Velocity_In_R,	&
                        &	Last_Velocity_In_Z,	&
                        &	Last_Angle,		&
                        &	HC_Temperature,		&
                        &	 Debug_HC_Neutral)
                   !Write (Output_Unit_Scratch,*) "After neutral step:",Current_R,Current_Z,Current_Cell,Current_Ring,Current_Velocity,HC_Temperature,Current_Angle
 
                   ! Increase neutral time step count.
                   Neutral_Time_Step_Count = Neutral_Time_Step_Count + 1.0
                   Eq_Total_Ion_Time_Steps = Eq_Total_Ion_Time_Steps +  Neutral_Time_Step /  Ion_Time_Step
                   Total_Time_Step_Count = Total_Time_Step_Count + 1.0
 
                   ! Check if the neutral particle has gone outside the wall boundary (returns logical Wall_Check_Result).
                   Call GA15B(Current_R,Current_Z,Wall_Check_Result, Num_Boundary_Points,1,GA15_Work,4* Max_Points,GA15_IndWork, &
                             &Max_Points,All_Wall_R_Points,All_Wall_Z_Points,TDum,XDum,YDum,6)
 
                   ! Check if inside the WBC geometry boundary and count if not.
                   If (hc_wbc_comp_option .gt. 0) Then
                      !Write (Output_Unit_Scratch,'(a,2i,80e)') "HC neutboundcheck:",NPROD-LPROD+1,IPROD,HC_Temperature,Current_Velocity, &
                      !& Current_R,Current_R -  Launch_R, &
                      !& Current_Z,Current_Z -  Launch_Z
 
                      If (Current_R .gt. ( Launch_R + hc_WBC_hori_bound) .or. &
                           & Current_R .lt. ( Launch_R - hc_WBC_hori_bound) .or. &
                           & Current_Z .gt. ( Launch_Z + hc_WBC_vert_bound)) Then
 
                         ! Add to WBC counting and stop following.
                         IFate = 15 ! NEUT IFATE = 7, neutral hit WBC boundary.
 
                      End If
                   End If
 
                ElseIf (Get_HC_Charge (Cur_HC_Spec) .gt. 0) Then
                   ! Ionized hydrocarbon species.
 
                   ! Check to be sure particle is still within the grid before trying to move it.
                   ! Otherwise, go straight to HC_Outside_Ion routine.
                   If (Current_S .gt. 0.0 .or. Current_S .lt. gksmaxs (Current_Ring)) Then
 
                      !Write (Output_Unit_Scratch,*) "Before ion step",Current_R,Current_Z,Current_S,Current_Cross,Prev_Cell, &
                      !& Current_Cell,Last_Ring,Current_Ring,Current_Angle,Current_Velocity_In_S,HC_Temperature,Segment_Index
 
                      ! Save last position before updating.
                      Last_R = Current_R
                      Last_Z = Current_Z
                      Last_Temperature = HC_Temperature
 
                      ! Quickly adjust transport coefficients.
                      Call Adjust_Taus (Cur_HC_Spec,H_Isotope_Composition,Sput_Weight,HC_Temperature,Current_Cell,Current_Ring,&
                                       &Current_S,Launch_Reg,Eq_Total_Ion_Time_Steps,Seed,Random_Numbers_Used,NRand,hc_v)
 
                      Call Ion_Move (			&
                           &	Current_R,		&
                           &	Current_Z,		&
                           &	Current_S,		&
                           &	Current_Cross,		&
                           &	Current_Cell,		&
                           &	Current_Ring,		&
                           &	Cur_HC_Spec,	&
                           &	H_Isotope_Composition,	&
                           &	Sput_Weight,		&
                           &	Current_Angle,		&
                           &	Current_Velocity_In_S,	&
                           &	HC_Temperature,		&
                           &	Current_Theta,		&
                           &	Last_R,			&
                           &	Last_Z,			&
                           &	Last_S,			&
                           &	Last_Cross,		&
                           &	Last_Cell,		&
                           &	Last_Ring,		&
                           &	Last_Angle,		&
                           &	Last_Velocity_In_S,	&
                           &	Random_Numbers_Used,	&
                           &    nrand,                  &
                           &    iprod,                  &
                           &    eq_total_ion_time_steps,&
                           &    hc_v,                   & ! HC velocity and temperature structure
                           &    Debug_HC_Ion)
 
                      Call GETRZ (Current_Cell,Current_Ring,Current_S,Current_Cross,Current_R,Current_Z, RZ_Opt)
 
                      !Write (Output_Unit_Scratch,*) "after ion step",Current_R,Current_Z,Current_S,Current_Cross,Prev_Cell, &
                      !& Current_Cell,Last_Ring,Current_Ring,Current_Angle,Current_Velocity_In_S,HC_Temperature,Segment_Index,gksmaxs (Current_Ring)
 
                      !if (Current_Cell .ne. Prev_Cell) Write (Output_Unit_Scratch,'(A,2I3,2F10.1,2F10.2,4F10.6)') "HC MOVED CELL in transport!!!!",current_cell,prev_cell, &s
                      !& neutral_Time_Step_Count,Ion_Time_Step_Count, Current_Velocity,Current_Velocity_In_S,Current_R,Current_Z,Current_S,Current_Cross
                      !if (Current_Ring .ne. Prev_Ring) Write (Output_Unit_Scratch,'(A,2I3,2F10.1,2F10.2,4F10.6)') "HC MOVED RING in transport!!!!",current_ring,prev_ring, &
                      !& neutral_Time_Step_Count,Ion_Time_Step_Count,Current_Velocity,Current_Velocity_In_S,Current_R,Current_Z,Current_S,Current_Cross
 
                      ! Increase ion time step count.
                      Ion_Time_Step_Count = Ion_Time_Step_Count + 1.0
                      Eq_Total_Ion_Time_Steps = Eq_Total_Ion_Time_Steps + 1.0
                      Total_Time_Step_Count = Total_Time_Step_Count + 1.0
 
                      ! Check if ion is inside the WBC geometry boundary and count if not.
                      If (hc_wbc_comp_option .gt. 0) Then
                         ! Find approximate R,Z for current location.
                         Call GetRZ (Current_Cell,Current_Ring,Current_S,Current_Cross,Current_R,Current_Z, RZ_Opt)
 
                         !Write (Output_Unit_Scratch,'(a,2i,80e)') "HC ionboundcheck:",NPROD-LPROD+1,IPROD,HC_Temperature,Current_Velocity_In_S, &
                         !& Current_R,Current_R -  Launch_R, &
                         !& Current_Z,Current_Z -  Launch_Z, &
                         !& MIN(Current_S,gksmaxs (Current_Ring) - Current_S),Current_Cross, &
                         !& MIN(Current_S,gksmaxs (Current_Ring) - Current_S)*SIN (ASIN (1 / gkbfs(Current_Cell,Current_Ring)))
 
                         If (Current_R .gt. ( Launch_R + hc_WBC_hori_bound) .or. &
                              & Current_R .lt. ( Launch_R - hc_WBC_hori_bound) .or. &
                              & Current_Z .gt. ( Launch_Z + hc_WBC_vert_bound)) Then
                            ! Add to WBC counting and stop following.
                            IFate = 27 ! DIV IFATE = 11, Ion hit WBC boundary.
                         End If
                      End If
                   End If
                End If
                ! Ensure no motion takes place. AMMOD REMOVE
                !Current_Velocity = 0.0 ! AMMOD REMOVE
                !Current_Velocity_In_S = 0.0 ! AMMOD REMOVE
                !Current_Velocity_In_R = 0.0 ! AMMOD REMOVE
                !Current_Velocity_In_Z = 0.0 ! AMMOD REMOVE
 
                ! Record time in particular HC state.
                Time_Steps_In_State = Time_Steps_In_State + 1.0
 


                If (.not.  Grid_Error) Then
                   ! Assign new properties and geometry for new location if necessary.
                   ! If the hydrocarbon has moved into a new cell, update the cell properties table.
                   If (Current_Ring .ne. Prev_Ring .or. Current_Cell .ne. Prev_Cell) Then
                      Call Initialize_Cell_Geom_Data (Current_Cell,Current_Ring)
                      Call Initialize_State_Prop_Data (Current_Cell,Current_Ring,Get_HC_Charge (Cur_HC_Spec))
                      ! Modify taus for new cell.  Check if only required for ion transport???
                      Call Modify_Taus (Cur_HC_Spec,H_Isotope_Composition,Sput_Weight,Current_Cell,Current_Ring)
                      ! jdemod
                      !Call Modify_Taus (Cur_HC_Spec,H_Isotope_Composition,Sput_Weight,HC_Temperature,Current_Cell,Current_Ring,&
                      !                & Current_S,Seed,Random_Numbers_Used,NRand)
                      Prev_Ring = Current_Ring
                      Prev_Cell = Current_Cell
                   End If
                End If
 
                ! React and count hydrocarbon depending on whether it is a neutral or ion, inside or outside the calculation grid.
                If (Get_HC_Charge (Cur_HC_Spec) .eq. 0) Then
 
                   If (Wall_Check_Result .gt. 0.0) Then ! Neutral is inside boundary.

                      ! For every time step spent as a neutral C - check to see if line profile data is being collected and update that
                      ! information as necessary. This test is hard coded as state 2! If neutral C changes from state 2 then this would need
                      ! to be changed as well. 
                      
                      ! Krieger IPP/07 - SUN compiler insists on 132 column limit
                      if (cur_hc_spec.eq.2.and.(.not.grid_error)) then   ! neutral C on computational grid
                         call hc_update_line_profile(current_cell,current_ring,current_r,current_z,&
                                                   & current_velocity_in_r/neutral_time_step, &
                                                   & current_velocity_in_z/neutral_time_step,sput_weight)
                         ! Krieger IPP/07 - SUN compiler insists on 132 column limit
                         if (debug_hc) &
                              & write(6,'(a,2i6,10(1x,g12.5))') 'HC LP:',current_cell,current_ring,current_r,current_z,hc_v%v(1), &
                              &                          hc_v%v(2),current_velocity_in_r/neutral_time_step, &
                              &                          current_velocity_in_z/neutral_time_step

                         if ( (abs(current_velocity_in_r/neutral_time_step-hc_v%v(1)).gt.0.01).or.&
                             &(abs(current_velocity_in_z/neutral_time_step-hc_v%v(2)).gt.0.01)) then 

                               ! Krieger IPP/07 - SUN compiler insists on 132 column limit
                               write(6,'(a,2i6,10(1x,g12.5))') 'HC LP ERROR:',current_cell,current_ring,current_r,current_z, &
                              &                         hc_v%v(1),hc_v%v(2), current_velocity_in_r/neutral_time_step,  &
                              &                         current_velocity_in_z/neutral_time_step

                         endif
                     

                      endif
                         
 
                      ! Input the current hydrocarbon species, its characteristics, and output the final hydrocarbon and its characteristics.
                      Call HC_Trans_Inside_Neutral (		&
                           &	IProd,				&
                           &	LPRod,				&
                           &	NProd,				&
                           &	Sput_Weight,			&
                           &	Double_Sput_Weight,		&
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
                           &	Current_Velocity_In_S,		& ! VINS - only OUT
                           &	Current_R,			& ! R
                           &	Current_Z,			& ! Z
                           &	Current_Cell,			& ! IK
                           &	Current_Ring,			& ! IR
                           &	Current_S,			& ! S - only OUT
                           &	Current_Cross,			& ! CROSS - only OUT
                           &	Current_Theta,			& ! THETA
                           &	Current_Angle,  		& ! Polar ANGLE theta.
                           &	Azimuthal_Angle,  		& ! Phi.
                           &	Launch_Angle,			& ! ANGLAN
                           &	Current_Tangent,		& ! TANGNT
                           &	Tangent_Launch,			& ! TANLAN
                           &	Last_Cell,			& ! IKLAST
                           &	Last_Ring,			& ! IRLAST
                           &	Time_Bin,			& ! IT
                           &	TStepHC,			& ! TSTEPN
                           &	IFate,				& ! IFATE
                           &    hc_v,                           & ! HC Velocity 
                           &	 Debug_HC_Neutral)
 
                   ! jdemod - try strictly less than 0.0  (Leave as .le. ) 
                   !ElseIf (Wall_Check_Result .lt. 0.0) Then 

                   ElseIf (Wall_Check_Result .le. 0.0) Then 
                      !write(output_unit_scratch,'(a,5g20.10)') 'Wall_check_result:',wall_check_result,Wall_Check_Result .le. 0.0

                      ! Wall_Check_Result is less than or equal to 0.0, particle is outside wall boundary. 
                      ! See if it gets stuck or reflects, then find it's new position, direction, hydrocabron species, and energy.
                      ! Write (Output_Unit_Scratch,*) "Neutral gone outside vessel:",Current_R,Current_Z,Current_Cell,Current_Ring,Cur_HC_Spec
 
                      ! Set indicator for vessel interaction has occurred at least once.
                      HC_Hit_Vessel = .True.
                      
                      ! jdemod - if the particle state changes upon wall interactions then zero out the reaction transition data 
                      !          for the particle upon a wall interaction so that only the evolution since the 
                      !          wall interaction is recorded

                      HC_Time_At_Production (IPROD,1:number_hc_species-1) = 0.0
                      HC_Energy_At_Production (IPROD,1:number_hc_species-1) = 0.0
                      HC_Kin_E_Add_At_Production (IPROD,number_hc_species-1) = 0.0


                      Call HC_Trans_Outside_Neutral (		&
                           &	IProd,				&
                           &	LPRod,				&
                           &	NProd,				&
                           &	Seed,				&
                           &	NRand,				&
                           &	Sput_Weight,			&
                           &	HC_Temperature,			& ! TEMN
                           &	Launch_Reg,			& ! M
                           &	NeutType,			&
                           &	Tries_Without_Change,		&
                           &	Eq_Total_Ion_Time_Steps,	& ! CIST
                           &	Neutral_Time_Step_Count,	&
                           &	Ion_Time_Step_Count,		&
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
                           &	Current_Angle,  		& ! ANGLE
                           &	Current_Velocity,		& ! VIN
                           &	Current_Velocity_In_R,		& ! XVELF
                           &	Current_Velocity_In_Z,		& ! YVELF
                           &	IFate,				& ! IFATE
                           &    hc_v,                           & ! HC Velocity structure
                           &	Debug_HC_Neutral)
 
                   End If
                ElseIf (Get_HC_Charge (Cur_HC_Spec) .gt. 0) Then
                   ! Check for HC in core rings first, then in SOL/PP.
                   If(Current_Ring.lt.Inner_SOL_Ring.or.(Current_S.gt.0.0.and.Current_S.lt.gksmaxs(Current_Ring)))Then! Ion is not below the target boundaries.
 
                      ! Input the current hydrocarbon species, its characteristics, and output the final hydrocarbon and its characteristics.
                      Call HC_Trans_Inside_Ion (		&
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
                           &	Current_Velocity_In_S,		& ! VINS - only OUT
                           &	Current_R,			& ! R
                           &	Current_Z,			& ! Z
                           &	Current_Cell,			& ! IK
                           &	Current_Ring,			& ! IR
                           &	Current_S,			& ! S - only OUT
                           &	Current_Cross,			& ! CROSS - only OUT
                           &	Current_Angle,  		& ! Polar ANGLE theta.
                           &	Azimuthal_Angle,  		& ! Phi.
                           &	Current_Theta,			& ! THETA
                           &	Launch_Angle,			& ! ANGLAN
                           &	Current_Tangent,		& ! TANGNT
                           &	Tangent_Launch,			& ! TANLAN
                           &	Last_R,				&
                           &	Last_Z,				&
                           &	Last_Cell,			& ! IKLAST
                           &	Last_Ring,			& ! IRLAST
                           &	Min_Z,				& ! XXX
                           &	Max_S,				& ! SSS
                           &	Max_S_Ion_Removal,		& ! SSSS
                           &	Time_Bin,			& ! IT
                           &	TStepHC,			& ! TSTEPN
                           &	IFate,				& ! IFATE
                           &    hc_v,                           & ! HC Velocity structure
                           &	Debug_HC_Ion)
                   Else 
                      ! Particle is on grid boundary. See if it gets stuck or reflects, then find it's new position, direction, hydrocabron species, and energy.
                      !Write (Output_Unit_Scratch,*) "Ion gone outside grid at targets:",Current_R,Current_Z,Current_S,Current_Cell,Current_Ring,gksmaxs (Current_Ring),Cur_HC_Spec
 
                      ! Set indicator for vessel interaction has occurred at least once.
                      HC_Hit_Vessel = .True.
                      
                      ! jdemod - if the particle state changes upon wall interactions then zero out the reaction transition data 
                      !          for the particle upon a wall interaction so that only the evolution since the 
                      !          wall interation is recorded

                      HC_Time_At_Production (IPROD,1:number_hc_species-1) = 0.0
                      HC_Energy_At_Production (IPROD,1:number_hc_species-1) = 0.0
                      HC_Kin_E_Add_At_Production (IPROD,number_hc_species-1) = 0.0

 
                      Call HC_Trans_Outside_Ion (		&
                           &	IProd,				&
                           &	LPRod,				&
                           &	NProd,				&
                           &	Sput_Weight,			& ! SPUTY
                           &	Seed,				& ! SEED
                           &	NRand,  			& ! NRAND
                           &	HC_Temperature, 		& ! TEMN
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
                           &	Last_R, 			& ! ROLD
                           &	Last_Z, 			& ! ZOLD
                           &	Current_S,			& ! S
                           &	Current_Cross,			& ! CROSS
                           &	Current_Theta,			& ! THETA
                           &	Segment_Index,			& ! ID
                           &	Current_Dir_Indicator,  	& ! IS
                           &	Current_Angle,  		& ! ANGLE
                           &	Current_Velocity,		& ! VIN
                           &	Current_Velocity_In_S,  	& ! VINS
                           &	Current_Velocity_In_R,		& ! XVELF
                           &	Current_Velocity_In_Z,		& ! YVELF
                           &	Random_Numbers_Used,		& ! KK
                           &	IFate,  			& ! IFATE
                           &    hc_v,                           & ! HC Velocity 
                           &	Debug_HC_Ion)
                   End If
 
                End If

 
                ! Record end-of-timestep density in common block structure if particle is on the grid.
                !If (Current_Cell .ne. Prev_Cell) Write (Output_Unit_Scratch,'(A,2I3,2F10.1,2F10.2,4F10.6)') "HC MOVED CELL in process!!!!",current_cell,prev_cell, &
                !& neutral_Time_Step_Count,Ion_Time_Step_Count,Current_Velocity,Current_Velocity_In_S,Current_R,Current_Z,Current_S,Current_Cross
                !If (Current_Ring .ne. Prev_Ring) Write (Output_Unit_Scratch,'(A,2I3,2F10.1,2F10.2,4F10.6)') "HC MOVED RING in process!!!!",current_ring,prev_ring, &
                !& neutral_Time_Step_Count,Ion_Time_Step_Count,Current_Velocity,Current_Velocity_In_S,Current_R,Current_Z,Current_S,Current_Cross
! slmod begin
                If ( .not.Grid_Error ) Then
!
!               If ( Grid_Error .ne. .True.) Then
! slmod end
                   ! Assign new properties and geometry for new location if necessary.
                   ! If the hydrocarbon has moved into a new cell, update the cell properties table.
                   If (Current_Ring .ne. Prev_Ring .or. Current_Cell .ne. Prev_Cell) Then
                      Call Initialize_Cell_Geom_Data (Current_Cell,Current_Ring)
                      Call Initialize_State_Prop_Data (Current_Cell,Current_Ring,Get_HC_Charge (Cur_HC_Spec))
                      ! Modify taus for new cell.
                      Call Modify_Taus (Cur_HC_Spec,H_Isotope_Composition,Sput_Weight,Current_Cell,Current_Ring)
                      ! jdemod
                      !Call Modify_Taus (Cur_HC_Spec,H_Isotope_Composition,Sput_Weight,HC_Temperature,Current_Cell,Current_Ring,&
                      !                 &Current_S,Seed,Random_Numbers_Used,NRand)
                      Prev_Ring = Current_Ring
                      Prev_Cell = Current_Cell
                   End If
 
                   HC_Density (Current_Cell,Current_Ring,Cur_HC_Spec) = &
                        &  HC_Density (Current_Cell,Current_Ring,Cur_HC_Spec) + Sput_Weight

                   ! Record detailed data if the subgrid option is active
                   if (subgrid_opt.ne.0) call hc_update_subgrid(current_r,current_z,cur_hc_spec,sput_weight)

 
                   ! Record C+ and C data.
                   ! jdemod - this use of hc_ddts is incompatible with use elsewhere - it will also result in double accounting in some cases
                   !          I am removing the code elsewhere and making use of hc_ddts consistent with hc_density
                   !        - also need to add code to properly update hc_temperature
                   !If (Cur_HC_Spec .eq. 1 .or. Cur_HC_Spec .eq. 2) Then 
                      HC_DDTS (Current_Cell,Current_Ring,Cur_HC_Spec) = &
                           &  HC_DDTS (Current_Cell,Current_Ring,Cur_HC_Spec) + HC_Temperature * Sput_Weight
                   ! Record CH data.
                   !ElseIf (Cur_HC_Spec .eq. 4) Then
                   !   HC_DDTS (Current_Cell,Current_Ring,Cur_HC_Spec) = &
                   !        &  HC_DDTS (Current_Cell,Current_Ring,Cur_HC_Spec) + HC_Temperature * Sput_Weight
                   !End If
                End If
 
                ! Loop back if cutoff time not yet reached.
                If ( Debug_HC_Neutral .or.  Debug_HC_Ion) Then
                   ! jdemod - do not need this at each time step
                   !Write (Output_Unit_Scratch,*) "Local_debug is on"
                   Do
                      TStepHC = TStepHC +  Print_Debug_x_CStep_Ion
                      If (TStepHC .gt. (Neutral_Time_Step_Count + Ion_Time_Step_Count)) Then
                         Exit
                      End If
                   End Do
 
                   ! Non-orthogonal
                   ! Krieger IPP/07 - SUN compiler insists on 132 column limit
                   if (get_hc_charge(cur_hc_spec).eq.0) then 
                      Write (Output_Unit_Scratch,9004) IProd, Neutral_Time_Step_Count , Ion_Time_Step_Count, &
                        & Current_Cell,Current_Ring, Cur_HC_Spec,Current_R,Current_Z,Current_S, Local_K, &
                        & Current_Theta,gksmaxs (Current_Ring),Current_Velocity,HC_Temperature, Local_SPara,Current_Cross,&
                        & Sput_Weight,Time_Bin,'Update Neutral'
                   else
                      Write (Output_Unit_Scratch,9004) IProd, Neutral_Time_Step_Count , Ion_Time_Step_Count, &
                        & Current_Cell,Current_Ring,Cur_HC_Spec,Current_R,Current_Z,Current_S, Local_K, &
                        & Current_Theta,gksmaxs (Current_Ring),Current_Velocity,HC_Temperature, Local_SPara,Current_Cross,&
                        & Sput_Weight,Time_Bin,'Update Ion'
                   endif

9004   format(I5 ,2F9.1, 3I4, 2F10.6, F9.3, F8.1 ,F10.5 ,F8.2, g12.3, F8.3, f8.3,f8.3,f7.2,I4,:,1X,A)


                End If
 
                ! Find total time passed so far since the original particle launch.
                Time_Value = Neutral_Time_Step_Count *  Neutral_Time_Step + Ion_Time_Step_Count *  Ion_Time_Step
 
                ! Print hydrocarbon state, r, z data in output file if selected in DIVIMP input file, option H30.
                If (hc_coord_print_option .eq. 1) Then
                   ! Check for ion
                   If (Get_HC_Charge (Cur_HC_Spec) .gt. 0.0) Then
                      Call GETRZ (Current_Cell,Current_Ring,Current_S,Current_Cross,Current_R,Current_Z, RZ_Opt)
                   End If
 
                   ! Print to file.
                   If (Get_HC_Charge (Cur_HC_Spec) .eq. 0) Then ! Neutral HC.
                      Write (Output_Unit_Location,Fmt=9101) Time_Value,Neutral_Time_Step_Count,Ion_Time_Step_Count,Cur_HC_Spec, &
                           & H_Isotope_Composition(1),H_Isotope_Composition(2),H_Isotope_Composition(3),Current_Cell,Current_Ring,&
                           & Current_R, &
                           & Current_Z,Current_S,Current_Cross,Current_Angle,Current_Velocity,HC_Temperature,gknbs(Current_Cell,&
                           & Current_Ring), &
                           & gktebs (Current_Cell,Current_Ring),gktibs (Current_Cell,Current_Ring)
                   Else ! Charged HC.
                      Write (Output_Unit_Location,Fmt=9202) Time_Value,Neutral_Time_Step_Count,Ion_Time_Step_Count,Cur_HC_Spec, &
                           & H_Isotope_Composition(1),H_Isotope_Composition(2),H_Isotope_Composition(3),Current_Cell,Current_Ring,&
                           & Current_R, &
                           & Current_Z,Current_S,Current_Cross,Current_Angle,Current_Velocity_In_S,HC_Temperature,gknbs(&
                           & Current_Cell,Current_Ring), &
                           & gktebs (Current_Cell,Current_Ring),gktibs (Current_Cell,Current_Ring),gkes (Current_Cell,Current_Ring)&
                           &, Local_Electric_Field, &
                           &  Local_Alphs, Local_Betas,  Local_HC_Change_State_Coll_Prob, &
                           &  Local_HC_Tau_Parallel_Inv,  Local_HC_Tau_Stopping_Inv,  Local_HC_Tau_Heating_Inv,  Local_SPara, &
                           &  Local_VPara,  Local_RConst,  Local_RGauss
                   End If
9101               Format(F10.7,F10.1,F10.1,I10,2X,I2,1X,I2,1X,I2,I10,I10,F10.6,F10.6,&
                         &F10.3,F10.5,F10.3,F10.3,F10.4,E10.3,F10.3,F10.3)
9202               Format(F10.7,F10.1,F10.1,I10,2X,I2,1X,I2,1X,I2,I10,I10,F10.6,F10.6,&
                         &F10.3,F10.5,F10.3,F10.3,F10.4,E10.3,F10.3,F10.3,12E10.3)
                End If
 
                ! Record end of time step position.
 
                ! jdemod - debug output
                !write(6,'(a,3i5,5g16.8)') 'HC_WALKS:',current_cell,current_ring,Cur_HC_Spec, current_r,current_z
 
                If ( HC_Walk_Count .lt. Max_Number_Walks) Then
                   !write (0,*) "MODDING:",Total_Time_Step_Count,Walks_Record_Freq,MOD(INT(Total_Time_Step_Count),Walks_Record_Freq)
                   If (MOD(INT(Total_Time_Step_Count),Walks_Record_Freq).eq.0) Then
                      HC_Walks ( HC_Walk_Count, 1) = Current_R
                      HC_Walks ( HC_Walk_Count, 2) = Current_Z
                      HC_Walks ( HC_Walk_Count+1, 1) =  Calc_Hi
                      HC_Walks ( HC_Walk_Count+1, 2) =  Calc_Hi
                      HC_Walk_Count =  HC_Walk_Count + 1
                   End If
                End If
                !write (0,*) "MAXITER:", Max_HC_Ion_Iter_To_TMax
                If ((Neutral_Time_Step_Count + Ion_Time_Step_Count) .ge.  Max_HC_Ion_Iter_To_TMax) Then
                   ! Reached cutoff time.
                   HC_Num_Reach_Max_Iter (Cur_HC_Spec,Launch_Reg) =  HC_Num_Reach_Max_Iter (Cur_HC_Spec,Launch_Reg) + Sput_Weight
                   HC_Tot_Temp_Reach_Max_Iter (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Temp_Reach_Max_Iter (Cur_HC_Spec,Launch_Reg) + &
                                                             & HC_Temperature * Sput_Weight
                   HC_Tot_Exist_Max_Iter (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Exist_Max_Iter (Cur_HC_Spec,Launch_Reg) + Sput_Weight
                   IFate = 1
                End If
 
                If (Total_Time_Step_Count .ge. Total_Time_Step_Count_Limit) Then
                   ! Reached cutoff time.
                   !write (0,*) "Reached total time step count limit:",Total_Time_Step_Count_Limit
                   HC_Num_Reach_Max_Iter (Cur_HC_Spec,Launch_Reg) =  HC_Num_Reach_Max_Iter (Cur_HC_Spec,Launch_Reg) + Sput_Weight
                   HC_Tot_Temp_Reach_Max_Iter (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Temp_Reach_Max_Iter (Cur_HC_Spec,Launch_Reg) + &
                                      & HC_Temperature * Sput_Weight
                   HC_Tot_Exist_Max_Iter (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Exist_Max_Iter (Cur_HC_Spec,Launch_Reg) + Sput_Weight
                   IFate = 1
                End If
 
                ! Only add to distance for non-breakup timesteps or breakups
                ! between fragments of equal charges.  This is due to the
                ! approximation in location on the calculational grid made in the
                ! neutral - ion and vise-versa breakup steps.
                If (Cur_HC_Spec .eq. Last_HC_Species .or. Get_HC_Charge (Last_HC_Species) .eq. Get_HC_Charge (Cur_HC_Spec)) Then
                   !Write (Output_Unit_Scratch,*) "Distance:",Last_HC_Species,Launch_Reg, Launch_R, Launch_Z,Last_R,Last_Z, HC_Cumu_CroFly_Dist_Trav_By_Pos (Last_HC_Species,Launch_Reg)
                   ! Add to distance travelled based on the velocity over the timestep.
                   If (Get_HC_Charge (Last_HC_Species) .eq. 0) Then
                      ! Neutral.
                      HC_Cumu_Real_Dist_Trav_By_VT (Last_HC_Species,Launch_Reg) = &
                           &  HC_Cumu_Real_Dist_Trav_By_VT (Last_HC_Species,Launch_Reg) + &
                           & ABS(Last_Velocity) *  Neutral_Time_Step
                      Tot_Velocity = Tot_Velocity + ABS(Last_Velocity)
                   ElseIf (Get_HC_Charge (Last_HC_Species) .ne. 0) Then
                      ! Ion.
                      HC_Cumu_Real_Dist_Trav_By_VT (Last_HC_Species,Launch_Reg) = &
                           &  HC_Cumu_Real_Dist_Trav_By_VT (Last_HC_Species,Launch_Reg) + &
                           & ABS(Last_Velocity_In_S) *  Ion_Time_Step
                      Tot_Velocity = Tot_Velocity + ABS(Last_Velocity_In_S)
                   End If
                   !write (0,*) "Cumu Real By VT", HC_Cumu_Real_Dist_Trav_By_VT (Last_HC_Species,Launch_Reg),Get_HC_Charge (Last_HC_Species), &
                   !& Tot_Velocity,Last_Velocity,Last_Velocity_In_S
 
                   HC_Cumu_Real_Dist_Trav_By_Pos (Last_HC_Species,Launch_Reg) = &
                        &  HC_Cumu_Real_Dist_Trav_By_Pos (Last_HC_Species,Launch_Reg) + &
                        & SQRT ((Last_R - Current_R)**2 + (Last_Z - Current_Z)**2)
                   !write (0,*) "Cumu Real By POS",IPROD,LPROD,NPROD,SQRT((Last_R - Current_R)**2 + (Last_Z - Current_Z)**2),Last_R,Current_R,Last_Z, &
                   !& Current_Z,IFate,Current_Cell,Current_Ring,Current_Velocity,Current_Velocity_In_S,HC_Temperature,Current_S,Current_Cross,Last_HC_Species,Cur_HC_Spec
                End If


                ! Record average data for new state.
                If (Cur_HC_Spec .ne. Last_HC_Species) Then
                   storage (IProd,1) = Current_Angle
                   storage (IProd,2) = Launch_Velocity
 
                   ! Distance travelled in previous state by first and last position.
                   HC_State_Dist_Trav_By_Pos (Last_HC_Species,Launch_Reg) = &
                        &  HC_State_Dist_Trav_By_Pos (Last_HC_Species,Launch_Reg) + &
                        & SQRT ((Last_R - State_R)**2 + (Last_Z - State_Z)**2)
                   !write (0,*) "State Dist By Pos",SQRT ((Last_R - State_R)**2 + (Last_Z - State_Z)**2),Last_R,State_R,Last_Z,State_Z					
 
                   ! jdemod - possible bug - doesn't this calculation of the distance travelled depend on whether it is a neutral or an ion?
                   !                       - maybe it isn't used??
                   ! Record distance travelled in state.
                   Particle_Distance = Particle_Distance + ABS (Current_Velocity) *  Neutral_Time_Step
 
                   ! Distance travelled in previous state by average velocity * time in state.
                   ! Tot_Velocity will have one extra velocity added, so subtract last_velocity.
                   If (Get_HC_Charge (Last_HC_Species) .eq. 0) Then
                      ! Neutral.
                      HC_State_Dist_Trav_By_VT (Last_HC_Species,Launch_Reg) = &
                           &  HC_State_Dist_Trav_By_VT (Last_HC_Species,Launch_Reg) + &
                           & (Tot_Velocity-ABS(Last_Velocity)) *  Neutral_Time_Step
                   ElseIf (Get_HC_Charge (Last_HC_Species) .ne. 0) Then
                      ! Ion.
                      HC_State_Dist_Trav_By_VT (Last_HC_Species,Launch_Reg) = &
                           &  HC_State_Dist_Trav_By_VT (Last_HC_Species,Launch_Reg) + &
                           & (Tot_Velocity-ABS(Last_Velocity_In_S)) *  Ion_Time_Step
                   EndIf
                   !write (0,*) "State Dist By Vel",Tot_Velocity, HC_State_Dist_Trav_By_VT (Last_HC_Species,Launch_Reg),Get_HC_Charge (Last_HC_Species)
 
                   ! Record C2 band transition. C2->C2+.
                   !							If (Last_HC_Species .eq. 13 .and. Cur_HC_Spec .eq. 12) Then
                   ! Add to TIZS weighting for previous charge if ionization  by e or p impact of CH->CH+ occurred.
                   ! First, add for "all" neutrals in position -1 and 0.
                   !								 HC_TIZS_C2 (Current_Cell,Current_Ring,-1) =  HC_TIZS_C2 (Current_Cell,Current_Ring,-1) + Sput_Weight
                   !								 HC_TIZS_C2 (Current_Cell,Current_Ring,0) =  HC_TIZS_C2 (Current_Cell,Current_Ring,0) + Sput_Weight
                   ! Next, add for secondaries, etc. if already re-released.
                   !								If ( HC_Hit_Vessel) Then
                   !									 HC_TIZS_C2 (Current_Cell,Current_Ring,0) =  HC_TIZS_C2 (Current_Cell,Current_Ring,0) + Sput_Weight
                   !								End If
                   !							End If
 
                   

                   ! jdemod - the following code doesn't make any sense - it appears that hc_tizs and hc_tizs_ch are incremented by sput_weight twice in the 
                   !          current cell if the particle has had a vessel wall interaction but only once otherwise - that doesn't make any sense from what I 
                   !          can tell since hc_hit_vessel is set true for any number of vessel wall impacts. In addition, every time a particle reaches CH it would
                   !          be recorded here anyway - and would continue to be double recorded if it had struck the wall anytime in the past. 
                   !
                   ! Record CH/CD band transition.
                   If (Last_HC_Species .eq. 4 .and. Cur_HC_Spec .eq. 3) Then
                      ! Add to TIZS weighting for previous charge if ionization  by e or p impact of CH->CH+ occurred.
                      ! First, add for "all" neutrals in position -1 and 0.
                      HC_TIZS_CH (Current_Cell,Current_Ring,-1) =  HC_TIZS_CH (Current_Cell,Current_Ring,-1) + Sput_Weight
                      HC_TIZS_CH (Current_Cell,Current_Ring,0) =  HC_TIZS_CH (Current_Cell,Current_Ring,0) + Sput_Weight
                      ! jdemod - double recording doesn't seem reasonable to me 
                      ! Next, add for secondaries, etc. if already re-released.
                      !If ( HC_Hit_Vessel) Then
                      !   HC_TIZS_CH (Current_Cell,Current_Ring,0) =  HC_TIZS_CH (Current_Cell,Current_Ring,0) + Sput_Weight
                      !End If
                   End If
 
                   ! Record C band transition.
                   If (Last_HC_Species .eq. 2 .and. Cur_HC_Spec .eq. 1) Then
                      ! Add to TIZS weighting for previous charge if ionization  by e or p impact of C->C+ occurred.
                      HC_TIZS (Current_Cell,Current_Ring,-1) =  HC_TIZS (Current_Cell,Current_Ring,-1) + Sput_Weight
                      HC_TIZS (Current_Cell,Current_Ring,0) =  HC_TIZS (Current_Cell,Current_Ring,0) + Sput_Weight
                      ! jdemod - double recording doesn't seem reasonable to me 
                      ! Next, add for secondaries, etc. if already re-released.
                      !If ( HC_Hit_Vessel) Then
                      !   HC_TIZS (Current_Cell,Current_Ring,0) =  HC_TIZS_CH (Current_Cell,Current_Ring,0) + Sput_Weight
                      !End If
                   End If
                   ! New state position information.
                   State_R = Last_R
                   State_Z = Last_Z
 
                   !
                   ! jdemod - S and CROSS are not useful or significant for neutral particles - there is no need to execute this
                   !          code - in addition - since this code is after the transport step - there is a reasonable chance that
                   !          the current R,Z are not in the current_cell, current_ring that well calculated for the call to gridpos
                   !          prior to the transport step. 
                   !
                   ! If neutral->neutral reaction occurred, need to find S and K for the new location.
                   If ((Get_HC_Charge (Last_HC_Species) .eq. 0) .and. (Get_HC_Charge (Cur_HC_Spec) .eq. 0)) Then
                      ! jdemod
                      ! Calculate S, Cross for new location.
                      ! Only do this if grid_error is not set - if the particle is not on the grid then there is no point in
                      ! executing this code. 
                      ! I am still of the opinion that this is not useful for neutrals and if it was it should be calculated
                      ! at every time step - not just for neutral to neutral transitions - however, just in case the HC code 
                      ! collects some data on the S,CROSS position of neutrals at state change - might as well ensure these
                      ! have a reasonable value. 
                      if (.not.grid_error) then 
                         Last_S = Current_S
                         Last_Cross = Current_Cross

                         Call getscross_approx (Current_R, Current_Z, Current_S, Current_Cross, Current_Cell, Current_Ring)
                      else
                         current_s = 0.0
                         current_cross = 0.0
                      endif
                   End If
 
                   !Write (Output_Unit_Scratch,*) "Times:",Last_HC_Species,Time_Steps_In_State, Neutral_Time_Step, Ion_Time_Step
                   ! Note addition instead of assignment to account for breakaway HC's from larger hydrocarbon species accounted for together with their original molecule.
                   If (Get_HC_Charge (Last_HC_Species) .eq. 0) Then
                      ! Neutral.
                      HC_State_Time (Last_HC_Species,Launch_Reg) =  HC_State_Time (Last_HC_Species,Launch_Reg) +  &
                                                                  & Neutral_Time_Step * Time_Steps_In_State
                      HC_State_Timesteps (Last_HC_Species,Launch_Reg) =  HC_State_Timesteps (Last_HC_Species,Launch_Reg) + &
                                                                       & Time_Steps_In_State
                      HC_Time_To_Prod(Cur_HC_Spec,Launch_Reg)=HC_Time_To_Prod(Cur_HC_Spec,Launch_Reg)+(Time_Value-Neutral_Time_Step)
                      If ( Neutral_Time_Step * Time_Steps_In_State .gt.  HC_Max_Time_In_State (Last_HC_Species,Launch_Reg)) Then
                         HC_Max_TimeSteps (Last_HC_Species,Launch_Reg) = Time_Steps_In_State
                         HC_Max_Time_In_State (Last_HC_Species,Launch_Reg) = &
                              &  Neutral_Time_Step * Time_Steps_In_State
                      End If
                   Else
                      ! Ion.
                      HC_State_Time (Last_HC_Species,Launch_Reg) =  HC_State_Time (Last_HC_Species,Launch_Reg) +  Ion_Time_Step * &
                                                                  & Time_Steps_In_State
                      HC_State_Timesteps (Last_HC_Species,Launch_Reg) =  HC_State_Timesteps (Last_HC_Species,Launch_Reg) + &
                                                                       & Time_Steps_In_State
                      HC_Time_To_Prod(Cur_HC_Spec,Launch_Reg)=HC_Time_To_Prod(Cur_HC_Spec,Launch_Reg)+(Time_Value-Ion_Time_Step)
                      If ( Ion_Time_Step * Time_Steps_In_State .gt.  HC_Max_Time_In_State (Last_HC_Species,Launch_Reg)) Then
                         HC_Max_TimeSteps (Last_HC_Species,Launch_Reg) = Time_Steps_In_State
                         HC_Max_Time_In_State (Last_HC_Species,Launch_Reg) = &
                              &  Ion_Time_Step * Time_Steps_In_State
                      End If
                   End If
 
                   ! Record creation of new species at the end of the previous timestep.
                   If (Get_HC_Charge (Cur_HC_Spec) .eq. 0) Then
                      Call Record_Evolve_Data (Cur_HC_Spec,Launch_Reg,Current_R,Current_Z,Current_Cell,Current_Ring,Current_S,&
                                              &Current_Cross,Current_Angle,Current_Velocity,HC_Temperature,Sput_Weight)
                   Else
                      Call Record_Evolve_Data (Cur_HC_Spec,Launch_Reg,Current_R,Current_Z,Current_Cell,Current_Ring,Current_S,&
                                              &Current_Cross,Current_Angle,Current_Velocity_In_S,HC_Temperature,Sput_Weight)
                   End If
 
                   !Write (Output_Unit_Scratch,*) "Timesteps in state:",Time_Steps_In_State,Cur_HC_Spec,Last_HC_Species
 
                   ! Record time and energy at production of new species.
                   ! Note, do only if re-release has not occurred.
                   ! jdemod - kin_energy_added is never zeroed so it is a cumulative total of the energy added as we run 
                   !          through the particles. 
                   !        - added zeroing of kin_energy_added so it is the total gained by the particle

                   !
                   ! jdemod - evolution data is only recorded if the particle has not hit the vessel wall - this results in 
                   !          no data being stored for puffed particles which have several interactions before reaching the plasma.
                   !          This functionality is being changed to completely record the last evolution record - not the first before
                   !          a wall impact and thus - past evolution records for the current particle are zeroed out at wall interactions. 
                   !
                   !If (.not.  HC_Hit_Vessel) Then
                      HC_Time_At_Production (IPROD,Cur_HC_Spec) = Time_Value
                      HC_Energy_At_Production (IPROD,Cur_HC_Spec) = HC_Temperature
                      HC_Kin_E_Add_At_Production (IPROD,Cur_HC_Spec) = Kin_Energy_Added

                      if (debug_hc) &
                           &  write(6,'(a,2i5,1x,a,10(1x,g12.5))') 'PROD:',iprod,cur_hc_spec,hc_ident_species(cur_hc_spec), &
                           &                                time_value,hc_temperature,kin_energy_added
                      

                   !End If
 
                   ! Restart state counter if state change has occurred.
                   Time_Steps_In_State = 0.0
                   Tot_Velocity = 0.0
                   !Write (Output_Unit_Scratch,*) "Evolution details:",Cur_HC_Spec,Last_HC_Species,Particle_Reach_State (Last_HC_Species, Launch_Reg), &
                   !&  HC_Num_Orig_Reach_State (Cur_HC_Spec,Launch_Reg), HC_Tot_Reach_State (Cur_HC_Spec,Launch_Reg)
 
                   ! Increment state counter for number of times particle reaches state (adds at each change of state).
                   ! jdemod - this should be sput_weight not 1.0
                   HC_Tot_Reach_State (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Reach_State (Cur_HC_Spec,Launch_Reg) + sput_weight
 
                   ! Check to see if particle has made it to this species previously.
                   If (Particle_Reach_State (Cur_HC_Spec, Launch_Reg) .lt. 1.0) Then
                      ! Made it to particle for the first time.
                      !HC_Num_Orig_Reach_State (Cur_HC_Spec,Launch_Reg) =  HC_Num_Orig_Reach_State (Cur_HC_Spec,Launch_Reg) + 1.0
                      !Particle_Reach_State (Cur_HC_Spec, Launch_Reg) = Particle_Reach_State (Cur_HC_Spec, Launch_Reg) + 1.0
                      HC_Num_Orig_Reach_State (Cur_HC_Spec,Launch_Reg) =  HC_Num_Orig_Reach_State (Cur_HC_Spec,Launch_Reg) + sput_weight
                      Particle_Reach_State (Cur_HC_Spec, Launch_Reg) = Particle_Reach_State (Cur_HC_Spec, Launch_Reg) + sput_weight
                   End If
 
                   ! Reset Last_HC_Species to prevent from recording average data over again.
                   Last_HC_Species = Cur_HC_Spec
                   Last_H_Isotope_Composition = H_Isotope_Composition
 
                End If
 
                ! Check for particular number of timesteps to return at.
                !If (Neutral_Time_Step_Count .gt. 100000 .and.  HC_Hit_Vessel) Then
                ! Stop following reflected particle.
                !IFate = 1 ! Overall time limit.
                !Exit
                !End If
 
                !If (Neutral_Time_Step_Count .gt. 100000) Then
                ! Stop following particle.
                !IFate = 1 ! Overall time limit.
                !Exit
                !End If
 
                If (IFate .ne. 0) Then
                   !write (0,*) "FATE ne 0"
                   If (HC_WBC_Comp_Option .ne. 0) Then
                      !write (0,*) "WBC COMP ON",IFATE,Launch_Reg,Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition),Cur_HC_Spec,HC_Temperature
                      ! Record WBC boundary and target cases.
 
                      If (IFate .eq. 10) Then
                         ! Write WBC data for neutral redep.
                         ! jdemod - changed 2 to 2.0
                         HC_Tot_Temp_At_WBC_HC_Target(Cur_HC_Spec,Launch_Reg)=HC_Tot_Temp_At_WBC_HC_Target(Cur_HC_Spec,Launch_Reg)+&
                              & Sput_Weight * (0.5 * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) * 1.67E-27 * &
                              & Current_Velocity * Current_Velocity / 1.602E-19 + 2.0 * HC_Temperature)
                         ! jdemod - format specifier missing brackets
                         write (Output_Unit_HC_Alert,'(A,I6,4E12.4)') "WBC_HC_Event: ",Cur_HC_Spec,Current_Velocity,HC_temperature,&
                              & Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition),Sput_Weight * &
                              & (0.5 * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) * 1.67E-27 * Current_Velocity * &
                              & Current_Velocity / 1.602E-19 + 2.0 * HC_Temperature)
 
                      ElseIf (IFate .eq. 11) Then
                         ! Write WBC data for neutral no-redep.
                         HC_Tot_Temp_At_WBC_HC_Boundary (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Temp_At_WBC_HC_Boundary (Cur_HC_Spec,&
                              & Launch_Reg) + &
                              & Sput_Weight * (0.5 * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) * 1.67E-27 * &
                              & Current_Velocity * Current_Velocity / 1.602E-19 + 2 * HC_Temperature)
                         ! jdemod - format specifier missing brackets
                         write (Output_Unit_HC_Alert,'(A,I6,4E12.4)') "WBC_HC_Event: ",Cur_HC_Spec,Current_Velocity,HC_temperature,&
                              & Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition),Sput_Weight * &
                              & (0.5 * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) * 1.67E-27 * Current_Velocity * &
                              & Current_Velocity / 1.602E-19 + 2.0 * HC_Temperature)
 
                      ElseIf (IFate .eq. 12) Then
                         ! Write WBC data for neutral redep.
                         HC_Tot_Temp_At_WBC_HC_Target(Cur_HC_Spec,Launch_Reg)=HC_Tot_Temp_At_WBC_HC_Target(Cur_HC_Spec,Launch_Reg)+&
                              & Sput_Weight * (0.5 * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) * 1.67E-27 * &
                              & Current_Velocity * Current_Velocity / 1.602E-19 + 2 * HC_Temperature)
                         ! jdemod - format specifier missing brackets
                         write (Output_Unit_HC_Alert,'(A,I6,4E12.4)') "WBC_HC_Event: ",Cur_HC_Spec,Current_Velocity,HC_temperature,&
                              & Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition),Sput_Weight * &
                              & (0.5 * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) * 1.67E-27 * Current_Velocity * &
                              & Current_Velocity / 1.602E-19 + 2.0 * HC_Temperature)
 
                      ElseIf (IFate .eq. 14) Then
                         ! Write WBC data for neutral redep.
                         HC_Tot_Temp_At_WBC_HC_Target(Cur_HC_Spec,Launch_Reg)=HC_Tot_Temp_At_WBC_HC_Target(Cur_HC_Spec,Launch_Reg)+&
                              & Sput_Weight * (0.5 * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) * 1.67E-27 * &
                              & Current_Velocity * Current_Velocity / 1.602E-19 + 2 * HC_Temperature)
                         ! jdemod - format specifier missing brackets
                         write (Output_Unit_HC_Alert,'(A,I6,4E12.4)') "WBC_HC_Event: ",Cur_HC_Spec,Current_Velocity,HC_temperature,&
                              & Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition),Sput_Weight * &
                              & (0.5 * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) * 1.67E-27 * Current_Velocity * &
                              & Current_Velocity / 1.602E-19 + 2.0 * HC_Temperature)
 
                      ElseIf (IFate .eq. 15) Then
                         ! WBC comparison addition for neut out of bounds.
                         HC_Tot_Temp_At_WBC_HC_Boundary (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Temp_At_WBC_HC_Boundary (Cur_HC_Spec,&
                              & Launch_Reg) + &
                              & Sput_Weight * (0.5 * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) * 1.67E-27 * &
                              & Current_Velocity * Current_Velocity / 1.602E-19 + 2 * HC_Temperature)
                         ! jdemod - format specifier missing brackets
                         write (Output_Unit_HC_Alert,'(A,I6,4E12.4)') "WBC_HC_Bound: ",Cur_HC_Spec,Current_Velocity,HC_temperature,&
                              & Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition),Sput_Weight * &
                              & (0.5 * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) * 1.67E-27 * Current_Velocity * &
                              & Current_Velocity / 1.602E-19 + 2.0 * HC_Temperature)
 
 
                      ElseIf (IFate .eq. 20) Then
                         ! Write WBC data for ion redep.
                         HC_Tot_Temp_At_WBC_HC_Target(Cur_HC_Spec,Launch_Reg)=HC_Tot_Temp_At_WBC_HC_Target(Cur_HC_Spec,Launch_Reg)+&
                              & Sput_Weight * (0.5 * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) * 1.67E-27 * &
                              & Current_Velocity_In_S * Current_Velocity_In_S / 1.602E-19 + 2.0 * HC_Temperature)
                         ! jdemod - format specifier missing brackets
                         write (Output_Unit_HC_Alert,'(A,I6,4E12.4)') "WBC_HC_Event: ",Cur_HC_Spec,Current_Velocity_In_S,&
                              & HC_temperature,Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition),Sput_Weight * &
                              & (0.5 * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) * 1.67E-27 * Current_Velocity_In_S * &
                              & Current_Velocity_In_S / 1.602E-19 + 2.0 * HC_Temperature)
 
                      ElseIf (IFate .eq. 22) Then
                         ! Write WBC data for ion redep.
                         HC_Tot_Temp_At_WBC_HC_Target(Cur_HC_Spec,Launch_Reg)=HC_Tot_Temp_At_WBC_HC_Target(Cur_HC_Spec,Launch_Reg)+&
                              & Sput_Weight * (0.5 * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) * 1.67E-27 * &
                              & Current_Velocity_In_S * Current_Velocity_In_S / 1.602E-19 + 2 * HC_Temperature)
                         ! jdemod - format specifier missing brackets
                         write (Output_Unit_HC_Alert,'(A,I6,4E12.4)') "WBC_HC_Event: ",Cur_HC_Spec,Current_Velocity_In_S,&
                              & HC_temperature,Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition),Sput_Weight * &
                              & (0.5 * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) * 1.67E-27 * Current_Velocity_In_S * &
                              & Current_Velocity_In_S / 1.602E-19 + 2.0 * HC_Temperature)
                      ElseIf (IFate .eq. 25) Then
                         ! Write WBC data for ion redep.
                         HC_Tot_Temp_At_WBC_HC_Target(Cur_HC_Spec,Launch_Reg)=HC_Tot_Temp_At_WBC_HC_Target(Cur_HC_Spec,Launch_Reg)+&
                              & Sput_Weight * (0.5 * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) * 1.67E-27 * &
                              & Current_Velocity_In_S * Current_Velocity_In_S / 1.602E-19 + 2 * HC_Temperature)
                         ! jdemod - format specifier missing brackets
                         write (Output_Unit_HC_Alert,'(A,I6,4E12.4)') "WBC_HC_Event: ",Cur_HC_Spec,Current_Velocity_In_S,&
                              & HC_temperature,Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition),Sput_Weight * &
                              & (0.5 * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) * 1.67E-27 * Current_Velocity_In_S * &
                              & Current_Velocity_In_S / 1.602E-19 + 2.0 * HC_Temperature)
 
                      ElseIf (IFate .eq. 26) Then
                         ! Write WBC data for ion redep.
                         HC_Tot_Temp_At_WBC_HC_Target(Cur_HC_Spec,Launch_Reg)=HC_Tot_Temp_At_WBC_HC_Target(Cur_HC_Spec,Launch_Reg)+&
                              & Sput_Weight * (0.5 * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) * 1.67E-27 * &
                              & Current_Velocity_In_S * Current_Velocity_In_S / 1.602E-19 + 2.0 * HC_Temperature)
                         ! jdemod - format specifier missing brackets
                         write (Output_Unit_HC_Alert,'(A,I6,4E12.4)') "WBC_HC_Event: ",Cur_HC_Spec,Current_Velocity_In_S,&
                              & HC_temperature,Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition),Sput_Weight * &
                              & (0.5 * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) * 1.67E-27 * Current_Velocity_In_S * &
                              & Current_Velocity_In_S / 1.602E-19 + 2.0 * HC_Temperature)
 
                      ElseIf (IFate .eq. 21) Then
                         ! Write WBC data for ion no-redep.
                         HC_Tot_Temp_At_WBC_HC_Boundary (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Temp_At_WBC_HC_Boundary (Cur_HC_Spec,&
                              & Launch_Reg) + &
                              & Sput_Weight * (0.5 * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) * 1.67E-27 * &
                              & Current_Velocity_In_S * Current_Velocity_In_S / 1.602E-19 + 2.0 * HC_Temperature)
                         ! jdemod - format specifier missing brackets
                         write (Output_Unit_HC_Alert,'(A,I6,4E12.4)') "WBC_HC_Event: ",Cur_HC_Spec,Current_Velocity_In_S,&
                              & HC_temperature,Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition),Sput_Weight * &
                              & (0.5 * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) * 1.67E-27 * Current_Velocity_In_S * &
                              & Current_Velocity_In_S / 1.602E-19 + 2.0 * HC_Temperature)
 
                      ElseIf (IFate .eq. 27) Then
                         ! WBC comparison addition for ion out of bounds.
                         HC_Tot_Temp_At_WBC_HC_Boundary (Cur_HC_Spec,Launch_Reg) =  HC_Tot_Temp_At_WBC_HC_Boundary (Cur_HC_Spec,&
                              & Launch_Reg) + &
                              & Sput_Weight * (0.5 * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) * 1.67E-27 * &
                              & Current_Velocity_In_S * Current_Velocity_In_S / 1.602E-19 + 2.0 * HC_Temperature)
                         ! jdemod - format specifier missing brackets
                         write (Output_Unit_HC_Alert,'(A,I6,4E12.4)') "WBC_HC_Boundary: ",Cur_HC_Spec,Current_Velocity_In_S,&
                              & HC_temperature,Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition),Sput_Weight * &
                              & (0.5 * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) * 1.67E-27 * Current_Velocity * &
                              & Current_Velocity / 1.602E-19 + 2.0 * HC_Temperature)
 
                      End If
                   End If
 
                   ! Finished with current hydrocarbon.  Record timing data.
                   Time_Value = Neutral_Time_Step_Count *  Neutral_Time_Step + Ion_Time_Step_Count *  Ion_Time_Step
                   Time_All_HCs (Launch_Reg) =  Time_All_HCs (Launch_Reg) + Time_Value
                   If (Time_Value .gt.  Max_Time_Any_HC (Launch_Reg)) Then
                      Max_Time_Any_HC (Launch_Reg) = Time_Value
                   End If
 
                   ! Add remaining distance travelled in previous state by first and last position.
                   HC_State_Dist_Trav_By_Pos (Cur_HC_Spec,Launch_Reg) = &
                        &  HC_State_Dist_Trav_By_Pos (Cur_HC_Spec,Launch_Reg) + &
                        & SQRT ((Last_R - State_R)**2 + (Last_Z - State_Z)**2)
                   !write (0,*) "State Dist By Pos",SQRT ((Last_R - State_R)**2 + (Last_Z - State_Z)**2),Last_R,State_R,Last_Z,State_Z					
 
                   ! Add remaining distance travelled in previous state by average velocity * time in state.
                   ! Tot_Velocity will have one extra velocity added, so subtract last_velocity.
                   If (Get_HC_Charge (Cur_HC_Spec) .eq. 0) Then
                      ! Neutral.
                      HC_State_Dist_Trav_By_VT (Cur_HC_Spec,Launch_Reg) = &
                           &  HC_State_Dist_Trav_By_VT (Cur_HC_Spec,Launch_Reg) + &
                           & (Tot_Velocity) *  Neutral_Time_Step
                   ElseIf (Get_HC_Charge (Cur_HC_Spec) .ne. 0) Then
                      ! Ion.
 
                      HC_State_Dist_Trav_By_VT (Cur_HC_Spec,Launch_Reg) = &
                           &  HC_State_Dist_Trav_By_VT (Cur_HC_Spec,Launch_Reg) + &
                           & (Tot_Velocity) *  Ion_Time_Step
                   EndIf
                   !write (0,*) "State Dist By Vel",Tot_Velocity,  Ion_Time_Step, HC_State_Dist_Trav_By_VT (Cur_HC_Spec,Launch_Reg),Get_HC_Charge (Cur_HC_Spec)
 
                   ! Record distance travelled as the crow flies from the point of birth.
                   HC_Cumu_CroFly_Dist_Trav_By_Pos (Cur_HC_Spec,Launch_Reg) = &
                        &  HC_Cumu_CroFly_Dist_Trav_By_Pos (Cur_HC_Spec,Launch_Reg) + &
                        & SQRT ((Current_R -  Launch_R)**2 + (Current_Z -  Launch_Z)**2)
                   !write (0,*) "Crowfly from point of birth",SQRT ((Last_R -  Launch_R)**2 + (Last_Z -  Launch_Z)**2), &
                   !& Last_R, Launch_R,Last_Z, Launch_Z
 
                   ! Write final position data to C+ position stats file.
                   If (IFate .eq. 13 .or. IFate .eq. 23) Then
                      Write (Output_Unit_Cpos_Pos,9191) Time_Value,Neutral_Time_Step_Count,Ion_Time_Step_Count,Last_HC_Species,&
                           & Last_H_Isotope_Composition(1),Last_H_Isotope_Composition(2), &
                           & Last_H_Isotope_Composition(3),Cur_HC_Spec,H_Isotope_Composition(1),H_Isotope_Composition(2),&
                           & H_Isotope_Composition(3),Current_Cell,Current_Ring, &
                           & Last_R,Current_R,Last_Z,Current_Z,Last_S,Current_S,Last_Cross,Current_Cross,Last_Angle,Current_Angle,&
                           & Last_Velocity_In_S,Current_Velocity,Last_Temperature,HC_Temperature
9191                  Format (E10.4,F10.1,F10.1,I10,2X,I2,1X,I2,1X,I2,I10,2X,I2,1X,I2,1X,I2,I10,I10,F10.6,F10.6,F10.6,F10.6,F10.6,&
                             &F10.6,F10.4,F10.4,F10.4,F10.4,F10.3,F10.3,F10.3,F10.3,F10.3,F10.3)
                   EndIf
 
                   ! Exit current hydrocarbon loop.
                   Exit
                End If
 
             End Do
 
             ! Check for end of HC stack, exit if finished.
             Size_Of_Stack = HC_Stack_Size()
             If (Size_Of_Stack .eq. 0) Then
                ! Removed the last hydrocarbon for this neutral impact.
                Exit
             Else ! Original HC was larger than CH4 and has previously broken up into 2 or more C components.
                ! Check that proper data is available for this to occur.
                If (hc_evolution_model_primary .le. 1 .and. hc_evolution_model_secondary .le. 1) Then
                   Write (Output_Unit_HC_Alert,*) "WARNING: Data does not support multiple C molecules."
                   Write (Output_Unit_HC_Alert,*) "Program stopping."
                   Stop
                End If
                ! Pull a particle off the stack.
                Call HC_Pull (Cur_HC_Spec,Current_R,Current_Z,Current_S,Current_Cross,Current_Angle,Current_Theta,Current_Velocity,&
                             &Current_Velocity_In_S,HC_Temperature)
 
                ! Assign current cell and ring if particle is on grid, set grid_error if not.
                Call gridpos (Current_Cell, Current_Ring, Current_R, Current_Z,.False., Grid_Error)
 
                ! Initialize state data and modify taus for starting state.
                Call Initialize_State_Prop_Data (Current_Cell,Current_Ring,Get_HC_Charge (Cur_HC_Spec)) ! Sets LFTS, LFPS, LFSS, LLLFPS, LTOLDS
                Call Modify_Taus (Cur_HC_Spec,H_Isotope_Composition,Sput_Weight,Current_Cell,Current_Ring)
                ! jdemod
                !Call Modify_Taus (Cur_HC_Spec,H_Isotope_Composition,Sput_Weight,HC_Temperature,Current_Cell,Current_Ring,Current_S,&
                !                 &Seed,Random_Numbers_Used,NRand)
 
                ! If neutral
                If (Get_HC_Charge (Cur_HC_Spec) .eq. 0) Then
                   ! Calculate directional velocity.
                   Current_Velocity_In_R = Current_Velocity * COS (Current_Angle) *  Neutral_Time_Step
                   Current_Velocity_In_Z = Current_Velocity * SIN (Current_Angle) *  Neutral_Time_Step
 
                Else
                   ! Ionized remnant.
                   Last_S = Current_S
                   Last_Cross = Current_Cross

                   Call getscross_approx (Current_R, Current_Z, Current_S, Current_Cross, Current_Cell, Current_Ring)
 
                   ! Initial ion velocity (H23).  Typically the same as CNEUTG DIVIMP input option.
                   If (hc_neut_ion_velocity .eq. 0) Then
                      ! Dead stop along S.
                      Current_Velocity_In_S = 0.0
                   ElseIf (hc_neut_ion_velocity .eq. 1 .or. hc_neut_ion_velocity .eq. 3) Then
                      ! +/- 0.5 Vprevious along S.
                      Random_Temp = geran (Seed)
                      NRand  = NRand + 1
                      Current_Velocity_In_S = SIGN (Current_Velocity * 0.5, Random_Temp - 0.5)
                      If (hc_neut_ion_velocity .eq. 3) Then
                         Random_Temp = geran (Seed)
                         NRand = NRand + 1
                         Current_Velocity_In_S = Current_Velocity_In_S * 2.0 * SQRT (Random_Temp)
                      End If
                   ElseIf (hc_neut_ion_velocity .eq. 2) Then
                      Current_Velocity_In_S = SIGN (Current_Velocity, 0.5 * gksmaxs (Current_Ring) - Current_S)
                   End If
 
                End If
 
             End If
 
             ! Load next HC released in higher hydrocarbon breakup.
 
 
          End Do
 
       End If ! IFATE = 12, not equal to 0 from initial boundary strike check will end up here.
 
    End If ! IFATE = 14, not equal 0 from launch_velocity will end up here.
 
    ! Current ion or sub-ion finished with.  End of Do-loop.
 
    ! Record end of time step position.
    If ( HC_Walk_Count .lt. Max_Number_Walks) Then
       HC_Walks ( HC_Walk_Count, 1) = Current_R
       HC_Walks ( HC_Walk_Count, 2) = Current_Z
       HC_Walks ( HC_Walk_Count+1, 1) =  Calc_Hi
       HC_Walks ( HC_Walk_Count+1, 2) =  Calc_Hi
       HC_Walk_Count =  HC_Walk_Count + 1
    End If
 
    ! Record maximum timesteps in state.
    If (Time_Steps_In_State .gt.  HC_Max_TimeSteps (Last_HC_Species,Launch_Reg)) Then
       HC_Max_TimeSteps (Last_HC_Species,Launch_Reg) = Time_Steps_In_State
    End If
 
    ! Record total and maximum ion equivalent timesteps in state.
    HC_Total_Ion_Teq_Iter (Launch_Reg) =  HC_Total_Ion_Teq_Iter (Launch_Reg) + Eq_Total_Ion_Time_Steps
    !Write (Output_Unit_Scratch,*) "Total ion eq tss:", HC_Total_Ion_Teq_Iter,Eq_Total_Ion_Time_Steps
    If (Eq_Total_Ion_Time_Steps .gt.  HC_Max_Ion_Teq_Iter (Launch_Reg)) Then
       HC_Max_Ion_Teq_Iter (Launch_Reg) = Eq_Total_Ion_Time_Steps
    End If
 
    ! Note addition instead of assignment to account for breakaway HC's from larger hydrocarbon species accounted for together with their original molecule.
    If (Get_HC_Charge (Last_HC_Species) .eq. 0) Then
       HC_State_Time (Cur_HC_Spec,Launch_Reg) =  HC_State_Time (Cur_HC_Spec,Launch_Reg) +  Neutral_Time_Step * Time_Steps_In_State
       HC_State_Timesteps (Cur_HC_Spec,Launch_Reg) =  HC_State_Timesteps (Cur_HC_Spec,Launch_Reg) + Time_Steps_In_State
       If ( Neutral_Time_Step * Time_Steps_In_State .gt.  HC_Max_Time_In_State (Cur_HC_Spec,Launch_Reg)) Then
          HC_Max_Time_In_State (Cur_HC_Spec,Launch_Reg) =  HC_Max_Time_In_State (Cur_HC_Spec,Launch_Reg) +  Neutral_Time_Step * &
                                                         & Time_Steps_In_State
       End If
    Else
       HC_State_Time (Cur_HC_Spec,Launch_Reg) =  HC_State_Time (Cur_HC_Spec,Launch_Reg) +  Ion_Time_Step * Time_Steps_In_State
       HC_State_Timesteps (Cur_HC_Spec,Launch_Reg) =  HC_State_Timesteps (Cur_HC_Spec,Launch_Reg) + Time_Steps_In_State
       If ( Ion_Time_Step * Time_Steps_In_State .gt.  HC_Max_Time_In_State (Cur_HC_Spec,Launch_Reg)) Then
          HC_Max_Time_In_State (Cur_HC_Spec,Launch_Reg) =  HC_Max_Time_In_State (Cur_HC_Spec,Launch_Reg) +  Ion_Time_Step * &
                                                         & Time_Steps_In_State
       End If
    End If
 
    ! Record min/max position data.
    HC_Tot_Min_Z_Reach =  HC_Tot_Min_Z_Reach + Min_Z * Sput_Weight
    HC_Tot_Max_S_Reach =  HC_Tot_Max_S_Reach + Max_S * Sput_Weight
    HC_Tot_Max_S_Reach_Ion_Removal =  HC_Tot_Max_S_Reach_Ion_Removal + Max_S_Ion_Removal * Sput_Weight
 
    ! If output is requested, write to files.
    If (hc_evolve_print_option .eq. 1) Then
       Write (Output_Unit_Evolve,9400) "IFate:",IFate, Particle_Fate (IFate)
9400   Format (2X,A,1X,I2,1X,A)
    End If
 
    If (HC_Coord_Print_Option .eq. 1) Then
       Write (Output_Unit_Location,9401) "IFate:",IFate, Particle_Fate (IFate)
9401   Format (2X,A,1X,I2,1X,A)
    End If
 
    ! Non-orthogonal.
    If ( Debug_HC_Ion) Then

       Write (Output_Unit_Scratch,9004) IProd, Neutral_Time_Step_Count , Ion_Time_Step_Count, Current_Cell,Current_Ring,&
            & Cur_HC_Spec,Current_R,Current_Z,Current_S, Local_K, &
            & Current_Theta,gksmaxs (Current_Ring),Current_Velocity,HC_Temperature, Local_SPara,Current_Cross,&
            & Sput_Weight,Time_Bin,'FATE:'//Particle_Fate (IFate)

       ! jdemod - none of the debug output matches the formatting in 9003
       !Write (Output_Unit_Scratch,9003) IProd,Neutral_Time_Step_Count,Ion_Time_Step_Count,Current_Cell,Current_Ring,Cur_HC_Spec,&
       !     & Current_R,Current_Z,Current_S, Local_K, &
       !     & Current_Theta,gksmaxs (Current_Ring),Current_Velocity,HC_Temperature, Local_SPara,Current_Cross,Sput_Weight,Time_Bin,&
       !     & IFate

    End If
 
    ! Report and record IFate data.
    HC_IFate_Count (Cur_HC_Spec,Launch_Reg,Ifate) =  HC_IFate_Count (Cur_HC_Spec,Launch_Reg,IFate) + 1.0
    !Write (Output_Unit_Scratch,*) "IFATE info: IFATE, DIVIMPIFATE, SPECIES, LAUNCH REGION, SUM FOR FATE:", &
    !& IFate, Particle_Fate (IFate),Cur_HC_Spec,Launch_Reg, HC_IFate_Count (Cur_HC_Spec,Launch_Reg,Ifate)
 
    ! Record MTC event frequency.
    If (MTC_Counter .gt. 10) Then
       HC_MTCTotCnt (11,1) =  HC_MTCTotCnt (11,1) + Sput_Weight
    Else
       HC_MTCTotCnt (MTC_Counter,1) =  HC_MTCTotCnt (MTC_Counter,1) + Sput_Weight
    End If
 
    If ( Debug_HC_Neutral .or. (IFate .eq. 14)) Then

       Write (Output_Unit_Scratch,9004) IProd, Neutral_Time_Step_Count , Ion_Time_Step_Count, Current_Cell,Current_Ring,&
            & Cur_HC_Spec,Current_R,Current_Z,Current_S, Local_K, &
            & Current_Theta,gksmaxs (Current_Ring),Current_Velocity,HC_Temperature, Local_SPara,Current_Cross,&
            & Sput_Weight,Time_Bin,Particle_fate(ifate)

       !Write (Output_Unit_Scratch,9003) IProd,Neutral_Time_Step_Count,Current_Cell,Current_Ring,Current_R,Current_Z, Local_K,&
       !     & Current_Velocity, &
       !     & HC_Temperature,Sput_Weight,(Current_Angle) *  Rad_In_A_Deg,Time_Bin, Particle_Fate (IFate)


9003   FORMAT(1X,I5,F9.1,2I4,2F9.5,F9.4,F8.1,F7.2,F5.2,F8.3,I3,:,1X,A)
    End If
 
    ! Finished following original hydrocarbon and it's daughter particles.  Store counters in transport data block.
    HC_Ion_Timesteps_Count (Launch_Reg) =  HC_Ion_Timesteps_Count (Launch_Reg) + Ion_Time_Step_Count
    HC_Neutral_Timesteps_Count (Launch_Reg) =  HC_Neutral_Timesteps_Count (Launch_Reg) + Neutral_Time_Step_Count
 
    Neutime = Neutral_Time_Step_Count *  Neutral_Time_Step
    Iontime = Ion_Time_Step_Count *  Ion_Time_Step
 
    Open (unit=92,file="temphc.txt")
    ! Write rate HC_Density data for injir
    Write (92,*) "Neutral_Time_Step_Count",Neutral_Time_Step_Count
    Write (92,*) "Ion_Time_Step_Count",Ion_Time_Step_Count
    Write (92,*) "Time_Steps_In_State",Time_Steps_In_State
    Write (92,*) "Eq_Total_Ion_Time_Steps",Eq_Total_Ion_Time_Steps
    Write (92,*) "Diag_Table % HC_Neutral_Timesteps_Count", HC_Neutral_Timesteps_Count
    Write (92,*) "Diag_Table % HC_Ion_Timesteps_Count", HC_Ion_Timesteps_Count
    Write (92,*) "Diag_Table % Time_All_HCs", Time_All_HCs
    Write (92,*) "Diag_Table % Max_Time_Any_HC", Max_Time_Any_HC
    Write(92,*)"Diag_Table%HC_State_TimeSteps(Cur_HC_Spec,1:2)",HC_State_TimeSteps(Cur_HC_Spec,1),HC_State_TimeSteps(Cur_HC_Spec,2)
    Write(92,*)"Diag_Table%HC_Max_TimeSteps(Cur_HC_Spec,1:2)",HC_Max_TimeSteps(Cur_HC_Spec,1),HC_Max_TimeSteps(Cur_HC_Spec,2)
    write (92,*) "INJIR",  INJ_Ring_Number
    Do Temp_Counter = 1,gnks( INJ_Ring_Number)
       write (92,19) "cell",Temp_Counter,"counts", HC_Density (Temp_Counter, INJ_Ring_Number,Cur_HC_Spec),"length", &
            &     gksb(Temp_Counter, INJ_Ring_Number)-gksb(Temp_Counter-1, INJ_Ring_Number),"KBFS",gkbfs(Temp_Counter, &
            &     INJ_Ring_Number),"BRatio", &
            &     gbratio(Temp_Counter, INJ_Ring_Number),"Center S",gkss(Temp_Counter, INJ_Ring_Number), &
            &     "area",gkareas(Temp_Counter, INJ_Ring_Number), &
            &     "ratio L/A",(gksb(Temp_Counter, INJ_Ring_Number)-gksb(Temp_Counter-1, INJ_Ring_Number))&
            &                  /max(1e-8,gkareas(Temp_Counter,INJ_Ring_Number)), &
            &     "ratio C/A", HC_Density (Temp_Counter, INJ_Ring_Number,Cur_HC_Spec)&
            &                  /max(1e-8,gkareas(Temp_Counter, INJ_Ring_Number)), &
            &     "ratio C/L", HC_Density (Temp_Counter, INJ_Ring_Number,Cur_HC_Spec)/(gksb(Temp_Counter, INJ_Ring_Number)-gksb(&
            &     Temp_Counter-1, INJ_Ring_Number))
    End Do
 
    write (92,*) "By Ring:", Inner_Wall_Ring
    Do Temp_Counter=1, Inner_Wall_Ring-1
       ! Find cell at injir S.
 
       !       FOR NOW SIMPLY FIND GRID POINT CLOSEST TO DESIRED
       !       INJECTION POSITION. THIS IS ALL THAT IS CURRENTLY DONE.
       !       THE PARTIAL DISPLACEMENTS FROM THE GRID POINTS ARE IGNORED.
       !
       NRAND = NRAND + 1
       CALL SURAND2 (SEED,1,Random_Temp)
       STMP=Random_Temp*(INJ_Area_Upper_Bound-INJ_Area_Lower_Bound)*gKSMAXS(Temp_Counter)+INJ_Area_Lower_Bound*gKSMAXS(Temp_Counter)
       NRAND = NRAND + 1
       CALL SURAND2 (SEED, 1, Random_Temp)
       IF (Random_Temp.GT.0.5.and.( Injection_Opt.eq.2.or. Injection_Opt.eq.5)) THEN
          !
          !         OUTER PLATE - OTHERWISE INNER - for option 2
          !
          STMP = gKSMAXS(Temp_Counter) - STMP
       ENDIF
 
       ! FIND NEAREST IK CORRESPONDING TO DISTANCE S ALONG CONTOUR INJI
       M = 1
758    IF (M.LT.gNKS(Temp_Counter).AND.STMP.GT.gKSS(M,Temp_Counter)) THEN
          M = M + 1
          GOTO 758
       ENDIF
759    IF (M.GT.1.AND.STMP.LE.gKSS(M-1,Temp_Counter)) THEN
          M = M - 1
          GOTO 759
       ENDIF
       IF (M.GT.1.AND.(STMP-gKSS(M-1,Temp_Counter).LT.gKSS(M,Temp_Counter)-STMP)) M = M - 1
       If (Temp_Counter .lt.  Inner_SOL_Ring) Then
          M = 14
       ELSE
          M = 34
       EndIf
 
       !         write (92,18) "Ring",Temp_Counter,"cell",M,"counts", HC_Density(M,Temp_Counter,Cur_HC_Spec), &
       !	&	"length",gksb(M,Temp_Counter)-gksb(M-1,Temp_Counter),"Center S",gkss(M,Temp_Counter), &
       !	&	"r centre",grs(M,Temp_Counter),"z centre",gzs(M,Temp_Counter),"cell width",2.0*cellwidth(M,Temp_Counter,-1,2),"proj cell width",2.0*cellwidth(M,Temp_Counter,-1,2)*gkbfs(M,Temp_Counter), &
       !	&	"KBFS",gkbfs(M, INJ_Ring_Number),"BRatio",gbratio(M,Temp_Counter),"area",gkareas(M,Temp_Counter), &
       !	&	"ratio L/A",(gksb(M,Temp_Counter)-gksb(M-1, INJ_Ring_Number))/gkareas(M, INJ_Ring_Number), &
       !	&	"ratio C/A", HC_Density(M,Temp_Counter,Cur_HC_Spec)/gkareas(M,Temp_Counter), &
       !	&	"ratio C/L", HC_Density(M,Temp_Counter,Cur_HC_Spec)/(gksb(M,Temp_Counter)-gksb(M-1,Temp_Counter))
    End Do
17  Format (a,1x,e9.3,1x,a,1x,e9.3,1x,a,1x,e9.3,1x,a,1x,I3)
18  Format (2(A,1X,I3,1X),20(A,1X,E11.5,1X))
19  Format (1(A,1X,I3,1X),20(A,1X,E11.5,1X))
    Close (unit=92)	
 
  End Subroutine HC_Transport
 
End Module HC_Follow
