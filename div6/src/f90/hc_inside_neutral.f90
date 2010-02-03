! -*-Mode:f90-*-
! Title: HC Inside Grid Neutral Transport Routine, Adam McLean, UTIAS, October 2001.
! Purpose: Transport a neutral within the DIVIMP calculational grid.
! Use: Used with the HC Follow wrapping routine.
! Notes:

Module HC_Inside_Neutral

  ! Required modules.
  use mtc
  ! Every good Fortran program has...
  Implicit None	

Contains

  Subroutine HC_Trans_Inside_Neutral (	&
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
       &	Current_Velocity_In_S,		& ! VINS
       &	Current_R,			& ! R
       &	Current_Z,			& ! Z
       &	Current_Cell,			& ! IK
       &	Current_Ring,			& ! IR
       &	Current_S,			& ! S
       &	Current_Cross,			& ! CROSS
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
       &        hc_v,                           & ! HC Velocity 
       &	Debug)

    ! External modules.
    Use ComHC ! Access user options, hydrocarbon common block.
    Use HC_Init_DIV_Data ! Gain access to cell/global/geometry/current data structures.
    Use HC_Init_DIV_Diag ! Data diagnostics.
    Use HC_Put ! External data functions.
    Use HC_Get ! External data functions.
    Use HC_LdDta ! Access hydrocarbon storage data.
    Use HC_NewSt ! Subroutines for hydrocarbon transitions.
    Use HC_Utilities ! Counting arrays.
    Use HC_Init_Lib_Data ! Hydrocarbon transition data.
    Use HC_Prompt ! Contains prompt deposition routine.
    Use HC_Freespace_Transition ! Transition physics.
    Use HC_Stack ! Load stored hydrocarbons upon completion of primary molecule.

    use hc_kinetics ! HC kinetics routines, options and types
    ! use hc_velocity_type ! HC Velocity type definition - included via the hc_kinetics module
    use bfield ! Magnetic field information

    ! Every good Fortran routine should have...
    Implicit None

    ! Define input/output variables
    Integer, Intent (In) :: IProd
    Integer, Intent (In) :: LPRod
    Integer, Intent (InOut) :: NProd
    Real, Intent (InOut) :: Sput_Weight
    Double Precision, Intent (InOut) :: Double_Sput_Weight
    Integer, Intent (In) :: LatIZ
    Integer, Intent (InOut) :: NatIZ
    Double Precision, Intent (In) :: Seed
    Integer, Intent (InOut) :: NRand
    Real, Intent (InOut) :: HC_Temperature ! TEMN
    Real, Intent (Out) :: Kin_Energy_Added
    ! jdemod - made InOut since it may be changed in hc_prompt_deposition 
    Integer, Intent (InOut) :: Launch_Reg ! M - may be changed in hc_prompt_deposition 
    Integer, Intent (In) :: NeutType
    Integer, Intent (InOut) :: Random_Numbers_Used ! KK
    Integer, Intent (In) :: Max_Random_Values ! KKLIM
    Integer, Intent (InOut) :: Tries_Without_Change
    Double Precision, Intent (InOut) :: Eq_Total_Ion_Time_Steps ! CIST
    Double Precision, Intent (InOut) :: Neutral_Time_Step_Count
    Double Precision, Intent (InOut) :: Ion_Time_Step_Count
    Double Precision, Intent (InOut) :: MTC_Time_Step_Count ! MTCCIST
    Integer, Intent (InOut) :: MTC_Counter
    Real, Intent (InOut) :: S_Start ! SSTART
    Real, Intent (InOut) :: Max_Velocity_Randoms ! RMAXS
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
    Real, Intent (Out) :: Current_Velocity_In_S ! VINS
    Real, Intent (InOut) :: Current_R ! R
    Real, Intent (InOut) :: Current_Z ! Z
    Integer, Intent (InOut) :: Current_Cell ! IK - can be changed in hc_prompt_deposition
    Integer, Intent (InOut) :: Current_Ring ! IR 
    Real, Intent (Out) :: Current_S ! S
    Real, Intent (Out) :: Current_Cross ! CROSS
    Real, Intent (InOut) :: Current_Theta ! THETA
    Real, Intent (InOut) :: Current_Angle ! Polar ANGLE theta.
    Real, Intent (InOut) :: Azimuthal_Angle ! Phi.
    Real, Intent (In) :: Launch_Angle ! ANGLAN
    Real, Intent (In) :: Current_Tangent ! TANGNT
    Real, Intent (In) :: Tangent_Launch ! TANLAN
    Integer, Intent (In) :: Last_Cell ! IKLAST
    Integer, Intent (In) :: Last_Ring ! IRLAST
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
    Integer :: Followed_State ! Keep track of carbon-containing products in a breakup.
    Integer :: Add_Followed_State ! Keep track of additional carbon-containing products in a breakup.
    Integer :: Launch_Target

    Real :: Temp_Velocity_Factor ! TMPVELF
    Real :: Temp_Current_Velocity ! VINTMP
    Real :: Temp_Current_Velocity_In_R ! XVELFTMP
    Real :: Temp_Current_Velocity_In_Z ! YVELFTMP
    Real :: HC_Timestep ! Equals either DIVIMP ion or neutral timestep, depending on hydrocarbon.
    Real :: Random_Temp ! Temporary random number.
    Real :: Time_Tot
    Integer :: Old_Cell
    Integer :: Old_Ring
    Real :: Old_S
    Real :: Old_Cross
    Real :: Old_R
    Real :: Old_Z
    Real :: Old_Theta
    Real :: Old_Angle
    Real :: Old_Velocity
    Real :: Old_Velocity_In_S
    Real :: Old_Temp

    ! Variables for hydrocarbon change of state.
    Integer :: Max_Without_Change
    Integer :: Counter
    ! jdemod - determine states_possible AFTER reaction chosen
    Integer, dimension (Highest_Carbon_Content) :: States_Possible
    !Integer, dimension (Max_Reactions_Per_State, Highest_Carbon_Content) :: States_Possible
    Integer, dimension (Max_Reactions_Per_State) :: Reactions_Possible
    ! jdemod - removed unused array
    !Real, dimension (Max_Reactions_Per_State) :: Total_Probabilities
    Real, dimension (Max_Reactions_Per_State) :: NS_Probabilities
    Real, dimension (Max_Reactions_Per_State) :: Reaction_Probabilities
    Real :: Random_Value
    Real :: Random_Value_1
    Real :: Random_Value_2
    Real :: Total_Probability

    ! jdemod - added number_reactions returned from ns_new_state as number of processes possible
    integer :: number_reactions
    integer :: in ! generic counting variable

    real :: hc_mtc_prob  ! Probability of an HC MTC collision

    ! External functions.
    Integer, External :: IPOS


    !write(0,*) 'Launch reg:',launch_reg

    If ( Grid_Error) Then
       !write (Output_Unit_Location,*) "Off the grid" ! Causes large case file size as injected into SOL to blow up.
       ! Particle not in grid - record void region density data.
       ! Estimate which region the neutral lies in -

       If ( Z_X_Point .gt.  Z_Plasma_Centre) Then	
          ! X-point up configuration
          If (Current_Z .ge.  Z_X_Point) Then
             ! Divertor region.
             If (Current_R .gt. grp( Inner_Target_Points) .and. Current_R .lt. grp ( Inner_Target_Points + 1)) Then
                ! Private Plasma void
                HC_DDVoid (2,Cur_HC_Spec,Launch_Reg) =  HC_DDVoid (2,Cur_HC_Spec,Launch_Reg) + Double_Sput_Weight
             Else
                ! Other Divertor void regions.
                HC_DDVoid (3,Cur_HC_Spec,Launch_Reg) =  HC_DDVoid (3,Cur_HC_Spec,Launch_Reg) + Double_Sput_Weight
             End If

          Elseif (Current_Z .lt.  Z_X_Point) Then
             ! Main plasma void region
             HC_DDVoid (1,Cur_HC_Spec,Launch_Reg) =  HC_DDVoid (1,Cur_HC_Spec,Launch_Reg) + Double_Sput_Weight
          End If

       ElseIf ( Z_X_Point .le.  Z_Plasma_Centre) Then
          ! X-point down configuration.
          If (Current_Z .le.  Z_X_Point) Then
             ! Divertor region.
             If (Current_Z .lt. grp( Inner_Target_Points) .and. Current_R .gt. grp ( Inner_Target_Points + 1)) then
                ! Private Plasma void
                HC_DDVoid (2,Cur_HC_Spec,Launch_Reg) =  HC_DDVoid (2,Cur_HC_Spec,Launch_Reg) + Double_Sput_Weight
             Else
                ! Other Divertor void regions.
                HC_DDVoid (3,Cur_HC_Spec,Launch_Reg) =  HC_DDVoid (3,Cur_HC_Spec,Launch_Reg) + Double_Sput_Weight
             EndIf

          ElseIf (Current_Z .gt.  Z_X_Point) Then
             ! Main plasma void region
             HC_DDVoid (1,Cur_HC_Spec,Launch_Reg) =  HC_DDVoid (1,Cur_HC_Spec,Launch_Reg) + Double_Sput_Weight
          EndIf

       EndIf

    ElseIf (.not.  Grid_Error) Then

       ! Check if entered main plasma, reached centre, reached
       ! target plates, or survived to time cutoff point.

       If ( HC_CflRin .and. Current_Ring .lt.  Inner_SOL_Ring) Then

          HC_CflRin = .False.
          HC_Num_Enter_Main_Plasma (Cur_HC_Spec,Launch_Reg) =  HC_Num_Enter_Main_Plasma (Cur_HC_Spec,Launch_Reg) + Sput_Weight
          HC_Elims (Current_Cell,3,Cur_HC_Spec) =  HC_Elims (Current_Cell,3,Cur_HC_Spec) + Sput_Weight

       ElseIf ( HC_CflRex .and. Current_Ring .ge.  Inner_SOL_Ring) Then
          If (.not.  HC_CflRin) Then
             HC_CflRex = .False.
             HC_Num_Enter_Main_Plasma (Cur_HC_Spec,Launch_Reg) =  HC_Num_Enter_Main_Plasma (Cur_HC_Spec,Launch_Reg) + Sput_Weight
             HC_Elims (Current_Cell,1,Cur_HC_Spec) =  HC_Elims (Current_Cell,1,Cur_HC_Spec) + Sput_Weight
          End If
       End If

       If (Current_Ring .eq. 1) Then
          HC_Num_Reach_Centre (Cur_HC_Spec,Launch_Reg) =  HC_Num_Reach_Centre (Cur_HC_Spec,Launch_Reg) + Sput_Weight
          IFate = 11
          Return

       ElseIf (Eq_Total_Ion_Time_Steps .gt.  Max_HC_Neut_Iter_To_TMax) Then
          !
          ! jdemod - this array has the wrong number of dimensions - HC_Num_at_tmax(hc_species,launch_reg)
          !HC_Num_At_TMax (Launch_Reg) =  HC_Num_At_TMax (Launch_Reg) + Sput_Weight
          HC_Num_At_TMax (cur_hc_spec,Launch_Reg) =  HC_Num_At_TMax (cur_hc_spec,Launch_Reg) + Sput_Weight
          IFate = 19
          Return
       End If

       !
       ! jdemod - both hc_density and hc_ddts are scored in hc_follow - these variables seem to be trying to do the same 
       !          thing but use different and incorrect indexing at least in the case of hc_ddts - I am commenting out this code
       !
       ! Score particle in DDLIMS in Cur_HC_Spec position.
       ! If time point reached, score time position also, and increment.
       ! Note:  Do only for C neutral.
       !If (Cur_HC_Spec .eq. 2) Then
       !   HC_DDLims (Current_Cell,Current_Ring,Get_HC_Charge (Cur_HC_Spec)) =  HC_DDLims (Current_Cell,Current_Ring,Get_HC_Charge (&
       !                                                                      &Cur_HC_Spec)) + Double_Sput_Weight
       !   HC_DDts (Current_Cell,Current_Ring,Get_HC_Charge (Cur_HC_Spec)) =  HC_DDts (Current_Cell,Current_Ring,Get_HC_Charge (&
       !                                                                    &Cur_HC_Spec)) + Double_Sput_Weight * HC_Temperature
       !End If

       ! Do also for CH/CD neutral.
       !If (Cur_HC_Spec .eq. 4) Then
       !   HC_DDLims_CH (Current_Cell,Current_Ring,Get_HC_Charge (Cur_HC_Spec)) =  HC_DDLims_CH (Current_Cell,Current_Ring,&
       !                                                                  &Get_HC_Charge (Cur_HC_Spec)) + Double_Sput_Weight
       !   HC_DDts_CH (Current_Cell,Current_Ring,Get_HC_Charge (Cur_HC_Spec)) =  HC_DDts_CH (Current_Cell,Current_Ring,&
       !                                                     &Get_HC_Charge (Cur_HC_Spec)) + Double_Sput_Weight * HC_Temperature
       !End If

       ! Record density separately for any chemically sputtered component.
       If (NeutType .eq. 2 .or. NeutType .eq. 5) Then
          ! This should always be true for chemically launched HCs, but not for free space launches.
          HC_ChemDen (Current_Cell,Current_Ring,Cur_HC_Spec) =  HC_ChemDen (Current_Cell,Current_Ring,Cur_HC_Spec) + Sput_Weight
       End If

       ! Time lapse statistics.
       If ( Number_Time_Steps .gt. 0) Then
          ! Determine the time bin for the particle
          Time_Bin = IPOS (Eq_Total_Ion_Time_Steps, gctimes(1,0),  Number_Time_Steps + 1)

          If (Eq_Total_Ion_Time_Steps .ge. gctimes(Time_Bin,0)) Then

             HC_Lims (Current_Cell,Current_Ring,Cur_HC_Spec,Time_Bin) =  HC_Lims (Current_Cell,Current_Ring,Cur_HC_Spec,Time_Bin) +&
                  & Sput_Weight
             Time_Bin = Time_Bin + 1
          End If

       End If

       ! Set new ionization probability.
       ! Check for production of C+.  If it has not occured, jump back for another iteration.

       ! Print debugging information if applicable.
       If (Debug) Then
          If (Eq_Total_Ion_Time_Steps .ge. TStepHC) Then
             Do
                TStepHC = TStepHC +  Print_Debug_x_CStep_Neutral
                If (TStepHC .gt. Eq_Total_Ion_Time_Steps) Then
                   Exit
                End If
             End Do
             Write (Output_Unit_Scratch,9004) IProd, Time_Bin, Eq_Total_Ion_Time_Steps, Current_Cell,Current_Ring,&
                  &Current_R,Current_Z, Local_K,&
                  &Current_Velocity,HC_Temperature,Sput_Weight,(Current_Angle + Current_Tangent) *  Rad_In_A_Deg
9004         Format (2(1X,I5),1x,g18.8,2(1X,I5),10(1X,F12.5))
          End If
       End If

       ! Check reaction possibilities for current state.
       If (Get_HC_Charge (Cur_HC_Spec) .eq. 0) Then
          HC_TimeStep =  Neutral_Time_Step
       Else
          HC_TimeStep =  Ion_Time_Step
       End If


       !
       ! ******************************* NEW STATE TEST *******************************************
       ! 


       ! Zero out necessary arrays for all reactions possible per state.
       States_Possible = 0
       Reactions_Possible = 0
       ! jdemod - remove unused array
       !Total_Probabilities = 0.0
       NS_Probabilities = 0.0
       Reaction_Probabilities = 0.0
       Total_probability = 0.0

       ! Find total hydrocarbon state change probability for current species, timestep, and grid location.
       !Write (Output_Unit_Scratch,*) "in_neut",gknbs (Current_Cell,Current_Ring) / 1.0E6,gknbs (Current_Cell,Current_Ring)			

       ! Note that E&L reaction rates are in cm^3/s.  Therefore, DIVIMP plasma
       ! densities must be reduced by 1E6 (from 1/m^3) for use.
       ! jdemod - removed states_possible from call - added number_reactions
       Call HC_New_State (Current_Cell,Current_Ring,Cur_HC_Spec,HC_Data_Type,gknbs (Current_Cell,Current_Ring) / 1.0E6,gktebs (&
            & Current_Cell,Current_Ring),gktibs (Current_Cell,Current_Ring), &
            & Back_Plasma_Ion_Mass,HC_Timestep,Reactions_Possible,Reaction_Probabilities,NS_Probabilities,&
            & Total_Probability,number_reactions)

       !Write (Output_Unit_Scratch,*) "NEWST NEUT:",Cur_HC_Spec,gknbs (Current_Cell,Current_Ring) / &
       !        & 1.0E6,gktebs (Current_Cell,Current_Ring),gktibs (Current_Cell,Current_Ring), &
       !        &  Back_Plasma_Ion_Mass,HC_Timestep,Reactions_Possible,States_Possible,&
       !        &  Reaction_Probabilities,NS_Probabilities,Total_Probability

       !Write (Output_Unit_Scratch,'(a,i6,2x,a,i6,1x,g14.6)') "NEWST NEUT:",number_reactions,&
       !      & hc_state_table(cur_hc_spec)%state_name, Cur_HC_Spec,total_probability
       !do in = 1,number_reactions
       !   write(output_unit_scratch,'(a,i6,3(1x,g12.5))') &
       !      &hc_state_table(hc_state_transform_table(reactions_possible(in))%end_c_states(1))%state_name,in,&
       !      &Reactions_Possible(in),Reaction_Probabilities(in),NS_Probabilities(in)
       !end do

       ! jdemod - this should be an HC debugging print out - not a particle debug printout
       !If (Debug_hc) Then
       !   Write (Output_Unit_Scratch,'(a,i5,2x,a,g12.4,f10.3,f10.3,f8.5,f12.1)') "New State:",&
       !        & Cur_HC_Spec, HC_Data_Type, gknbs (Current_Cell, Current_Ring), &
       !        & gktebs (Current_Cell, Current_Ring), gktibs (Current_Cell, Current_Ring), total_probability,HC_Timestep
       !End If


       ! Check for allowance of hydrocarbon state transitions.
       If (hc_disable_transitions .eq. 1) Then
          ! Transitions disabled.
          Total_Probability = 0.0
       End If

       ! Check for particular number of timesteps to force change of state.
       !If (Neutral_Time_Step_Count .eq. 100) Then
       !Total_Probability = 1.0
       !Neutral_Time_Step_Count = Neutral_Time_Step_Count + 1
       !End If
       !Total_Probability = Total_Probability / 5.0

       hc_mtc_prob = mtc_prob(Back_Plasma_Ion_Mass,Find_HC_Mass(Cur_HC_Spec,H_Isotope_Composition),& 
            & current_cell,current_ring,neutral_time_step)

       Random_Numbers_Used = Random_Numbers_Used + 1	
       If (granv (Random_Numbers_Used) .ge. ( Hc_mtc_prob + Total_Probability)) Then
          ! No reaction or MTC occurs.
          Tries_Without_Change = Tries_Without_Change + 1
          If (Tries_Without_Change .gt.  Max_HC_Neut_Iter_To_TMax) Then
             Write (Output_Unit_HC_Alert,*) "Error:  Probability for reaction too low:", Total_Probability, Max_HC_Neut_Iter_To_TMax
          End If
          Return
       End If


       !
       ! jdemod - CODE RETURNS IF NO CHANGE OCCURS
       !
       !
       ! NOTE: at the present time MTC is not implemented for HC - it APPEARS as if it should be but the array HC_mtcprob is set to zero in the
       !       initialization stage and never set to a value later. So - MTC is not active for HC fragments at any time. 
       !       - need to implement a probability calculation
       !       - need to properly combine this with the change of state probability
       !       - need to reimplement the MTC collision code so that one block of code can be used anywhere in DIVIMP
       !       - need to implement an option to turn it on and off at a global level
       !       - HC_mtcprob is mass dependent and thus changes with BOTH cell and species - not efficient to pre-calculate it - I have created a mtc module to
       !         handle the mtc process in the hc code, neut and neutone. 

       ! State change event has occurred.
       !
       ! 1) hydrocarbon evolution
       ! 2) charge exchange collision (?) - not energetically possible
       ! 3) momentum transfer collision
       !
       ! Check for and deal with momentum transfer collision first -
       ! Since the particle will continue to be tracked as a neutral
       ! after the direction of it's velocity is adjusted.
       !
       ! Note: the same random number is used since it still decides
       ! randomly between these probabilities there is no
       ! need to draw another - simply check (for now) if
       ! it is greater than the ionization probability since
       ! the only additional state change currently implemented
       ! (11/14/97) is a momentum transfer collision.

       ! jdemod - if an event has occured and it is NOT a change of state then it has to be a momentum transfer collision. 
       !          This code was modified from the existing in order to follow the same rules as the DIVIMP code. 
       !        - Total_probability is only the total probability for an HC change of state returned by HC_newst - the mtc probability
       !          is saved in hc_mtc_prob

       If (granv(Random_Numbers_Used) .gt. total_probability) Then
          !If (granv(Random_Numbers_Used) .lt.  HC_MTCProb (Current_Cell,Current_Ring,Cur_HC_Spec)) Then

          !   call execute_mtc(1,MTC_Counter,MTC_Time_Step_Count,Eq_Total_Ion_Time_Steps,HC_MTCInf,& 
          !            & Current_Velocity_In_R,Current_Velocity_In_Z, &
          !            & Sput_Weight,Current_Velocity,HC_Temperature,Impurity_Neutral_Vel_Opt, &
          !            & neutral_time_step,Random_Numbers_Used,Find_HC_Mass(Cur_HC_Spec,H_Isotope_Composition),&
          !            & current_cell,current_ring,hc_v%v)
          call execute_mtc(1,MTC_Counter,MTC_Time_Step_Count,Eq_Total_Ion_Time_Steps,HC_MTCInf,& 
               & Current_Velocity_In_R,Current_Velocity_In_Z, &
               & Sput_Weight,Current_Velocity,HC_Temperature,Impurity_Neutral_Vel_Opt, &
               & neutral_time_step,Random_Numbers_Used,Find_HC_Mass(Cur_HC_Spec,H_Isotope_Composition),&
               & gktibs(current_cell,current_ring),hc_v%v)

          return

       endif


       ! Proceed to evolve hydrocarbon.
       ! Save previous state.
       Last_HC_Species = Cur_HC_Spec
       Last_H_Isotope_Composition = H_Isotope_Composition

       ! Get random number to decide which state is the destination.
       Random_Numbers_Used = Random_Numbers_Used + 1 ! Add to KK counter.
       Random_Temp = granv (Random_Numbers_Used)

       ! jdemod - cumulative probabilities now returned by hc_new_state
       ! Convert normalized probability values to sums of state change given all previous reactions.
       !Do Counter=2, size (NS_Probabilities), 1
       !   NS_Probabilities (Counter) = NS_Probabilities (Counter) + NS_Probabilities (Counter-1)
       !End Do

       ! Find which state we've moved into with another random value.
       Counter=1
       Do While (Random_Temp .gt. NS_Probabilities (Counter))
          ! Increment counter.
          Counter = Counter + 1
       End Do

       ! Reset tries without change counter.
       Tries_Without_Change = 0

       ! jdemod
       ! Copy states possible from hc_state_transform_table
       do in = 1, highest_carbon_content
          states_possible(in) = hc_state_transform_table(reactions_possible(counter))%end_c_states(in)
       end do

       ! Indicate changed state.  Note:  Always take the first carbon atomic/molecular product,
       ! then check for additional components and store in the stack.
       Followed_State = 1
       ! jdemod 
       Cur_HC_Spec = States_Possible (Followed_State)
       !Cur_HC_Spec = States_Possible (Counter, Followed_State)
       !Cur_HC_Spec = 8
       If (Cur_HC_Spec .eq. 0) Then
          Write (Output_Unit_HC_Alert,*) "Error in HC_Inside_Neutral: Evolution data:",Reactions_Possible, States_Possible, &
               &Reaction_Probabilities, NS_Probabilities, Total_Probability
          Write (Output_Unit_HC_Alert,*) "Program stopping."
          Stop
       End If
       !Write (Output_Unit_Scratch,*) "HC Change Of State from neutral:",Last_HC_Species,Cur_HC_Spec,Current_R,Current_Z,Current_Cell,Current_Ring, HC_MTCProb (Current_Cell,Current_Ring,Cur_HC_Spec),Total_Probability

       ! Record left-over hydrogen based on background plasma/mass reduction and update H_Isotope_Composition.
       Call Record_Hydrogen_Release (Current_Cell,Current_Ring,Last_HC_Species,Reactions_Possible (Counter),Cur_HC_Spec,&
            & H_Isotope_Composition,Seed,NRand)

       ! Record change of state.
       HC_Reaction_Count (Last_HC_Species,Reactions_Possible (Counter),Launch_Reg) =  HC_Reaction_Count (Last_HC_Species,&
            & Reactions_Possible (Counter),Launch_Reg) + 1.0

       Sput_Weight = Calc_Sputy (Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition))
       Double_Sput_Weight = DBLE (Sput_Weight)

       ! Modify commons for cell tau parallel, stopping, and heating for current cell.
       Call Initialize_State_Prop_Data (Current_Cell,Current_Ring,Get_HC_Charge (Cur_HC_Spec))
       Call Modify_Taus (Cur_HC_Spec,H_Isotope_Composition,Sput_Weight,Current_Cell,Current_Ring)
       ! jdemod
       !Call Modify_Taus (Cur_HC_Spec,H_Isotope_Composition,Sput_Weight,HC_Temperature,Current_Cell,Current_Ring,Current_S,Seed,&
       !                 & Random_Numbers_Used,NRand)

       ! Now, load additional particles from higher hydrocarbon neutral breakup into stack storage.
       Add_Followed_State = Followed_State + 1

       Do While (Add_Followed_State .le. Highest_Carbon_Content)
          ! Check if molecule broke into more than one carbon-containing product.
          !Write (Output_Unit_Scratch,*) "FOT:",Counter,Add_Followed_State,States_Possible (Counter, Add_Followed_State)
          ! jdemod
          If (States_Possible (Add_Followed_State) .ne. 0) Then
             !If (States_Possible (Counter, Add_Followed_State) .ne. 0) Then
             ! Found another carbon product.  Save required characteristics on the stack.
             Write (Output_Unit_Scratch,*) "SHOULD never be here for E&L database - humungous error"

             ! jdemod
             !Call HC_Push (States_Possible (Counter, Add_Followed_State), Current_R, Current_Z, Current_S, Current_Cross, &
             Call HC_Push (States_Possible (Add_Followed_State), Current_R, Current_Z, Current_S, Current_Cross, &
                  & Current_Angle, Current_Theta, Current_Velocity, Current_Velocity_In_S, HC_Temperature)

             ! Increment counter as there may be another carbon-carrying component.
             Add_Followed_State = Add_Followed_State + 1
          Else
             ! No carbon-containing particles found.  Exit the loop.
             Exit
          End If

       End Do

       !---------------------------------------------------------------------------------------------
       !
       !------------------------ NEUTRAL TO ION -----------------------------------------------------
       !
       !---------------------------------------------------------------------------------------------

       ! NEUTRAL TO ION. Check for ionization and recombination and alter grid coordinates accordingly.
       If((Last_HC_Species.ne.Cur_HC_Spec).and.(Get_HC_Charge(Last_HC_Species).eq.0).and.(Get_HC_Charge(Cur_HC_Spec).ne.0))Then
          ! Particle began as neutral, ionization occurred.
          ! Note: STmp = 0.0 here as during non-injection transitions, getsapprox should always be used.
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

          ! jdemod - add new hc kinetics option
          if (hc_kinetics_opt.eq.0) then 
             Call Neut_Ion_FreeSpace_Location (Reactions_Possible (Counter),Followed_State,&
                  &Current_R,Current_Z,Current_Cell,&
                  &Current_Ring,0.0,Current_S,Current_Cross,Current_Theta,S_Start)
             Call Neut_Ion_FreeSpace_Angle (Reactions_Possible (Counter),Followed_State,Current_Angle)

             ! Base new energy on neutral velocity, not new ion velocity.
             Call Neut_Ion_FreeSpace_Energy (Reactions_Possible (Counter),Followed_State,&
                  &H_Isotope_Composition,Current_Cell,&
                  &Current_Ring,Current_Velocity,HC_Temperature,Kin_Energy_Added)
             Call Neut_Ion_FreeSpace_Velocity (Reactions_Possible (Counter),Followed_State,&
                  &H_Isotope_Composition,Current_Ring,&
                  &Current_S,Current_Velocity,Current_Velocity_In_S, Seed, NRand)
          elseif (hc_kinetics_opt.eq.1) then 


             if (debug_kinetics) &
                  & write(6,'(a,10g12.5)') 'TP:INSIDE NEUT:',total_probability,hc_mtc_prob

             call update_hc_kinetics(hc_v,ni_reaction,current_cell,current_ring,  & 
                  & cur_hc_spec,reactions_possible(counter),followed_state,h_isotope_composition,  &
                  & current_r,current_z,current_s, current_theta,current_cross,current_angle,  &
                  & current_velocity_in_s,current_velocity_in_r,current_velocity_in_z, &
                  & current_velocity,kin_energy_added,s_start, &
                  & neutral_time_step,ion_time_step,sput_weight,iprod)



          endif


          ! Record ionization has occured at least once.
          Particle_Ionized = Particle_Ionized + 1

          ! jdemod - changed line order to group code dependent on neuttype.ne.0
          ! Record Ring and target index from which ion originated
          ! for target or wall launches.
          if (particle_ionized.eq.1) then 

             ! Record position of first ionization to pass back to DIVIMP.
             First_Ioniz_Ring = Old_Ring
             First_Ioniz_Cell = Old_Cell

             If (NeutType .ne. 0) Then

                HC_WTSource ( Starting_Index,Old_Ring,1,NeutType) = &
                     &  HC_WTSource ( Starting_Index,Old_Ring,1,NeutType) + Sput_Weight

                HC_WTSource ( Starting_Index,Old_Ring,2,NeutType) = &
                     &HC_WTSource(Starting_Index,Old_Ring,2,NeutType)+MIN(Current_S,gksmaxs(Current_Ring)-Current_S)*Sput_Weight
                !Write (Output_Unit_Scratch,*) "SOURCING2:", HC_WTSource ( Starting_Index,Old_Ring,1,NeutType),NeutType,Cur_HC_Spec

                ! Accumulate data on the ionization of eroded particles from each wall element,NeutType,Cur_HC_Spec
                HC_WallsE_I ( Launch_Wall_Index,Cur_HC_Spec) = &
                     &  HC_WallsE_I ( Launch_Wall_Index,Cur_HC_Spec) + Sput_Weight
             endif

          End If

          ! Print hydrocarbon evolution data to file if option H31 selected in input file.
          If (HC_Evolve_Print_Option .eq. 1) Then
             ! Print option activated.
             Time_Tot = Eq_Total_Ion_Time_Steps *  Ion_Time_Step
             Write (Output_Unit_Evolve,9300) Time_Tot,Neutral_Time_Step_Count,Ion_Time_Step_Count,Last_HC_Species,&
                  & Last_H_Isotope_Composition(1),Last_H_Isotope_Composition(2), &
                  & Last_H_Isotope_Composition(3),Cur_HC_Spec,H_Isotope_Composition(1),H_Isotope_Composition(2),&
                  & H_Isotope_Composition(3),Current_Cell,Current_Ring, &
                  & Old_R,Current_R,Old_Z,Current_Z,Old_S,Current_S,Old_Cross,Current_Cross,Old_Angle,Current_Angle,Old_Velocity,&
                  & Current_Velocity_In_S,Old_Temp,HC_Temperature
9300         Format (E10.4,F10.1,F10.1,I10,2X,I2,1X,I2,1X,I2,I10,2X,I2,1X,I2,1X,I2,I10,I10,F10.6,F10.6,F10.6,F10.6,F10.6,F10.6,&
                  &F10.4,F10.4,F10.4,F10.4,F10.3,F10.3,F10.3,F10.3,F10.3,F10.3)
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
             Write (Output_Unit_Scratch,9003) IProd,0.0,Current_Cell,Current_Ring,Cur_HC_Spec,Current_R,Current_Z,Current_S, &
                  & Local_K,Current_Theta,gksmaxs (Current_Ring), &
                  & Current_Velocity_In_S,HC_Temperature, Local_SPara,Current_Cross,Sput_Weight,Time_Bin,'HC Ion Appeared'
9003         FORMAT(1X,I5,F9.1,2I3,I2,2F9.5,F8.3,2F6.2,F8.3,1P,E15.8,0P,F7.1,1P,E8.1,0P,F8.5,F5.2,I2,:,1X,A,:,F8.5)
9005         FORMAT(1X,'--HC------TIME-IK-IR-HC','----R--------Z-------S------K---THETA--SMAX---','---DRIFTVEL------TEMI-PARADIFF-&
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
                HC_Stopped_Follow_Ion_In_Core =  HC_Stopped_Follow_Ion_In_Core + Sput_Weight
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

             If (.not.  HC_InCore) Then
                HC_InCore = .True.
                HC_Num_Entered_Core (Cur_HC_Spec,Launch_Reg) =  HC_Num_Entered_Core (Cur_HC_Spec,Launch_Reg) + Sput_Weight

                HC_Ion_Core_Density ( Launch_Cell, Launch_Ring) = &
                     &  HC_Ion_Core_Density ( Launch_Cell, Launch_Ring) + Sput_Weight

                If (NeutType .gt. 0) Then
                   HC_WTSource ( Starting_Index, Launch_Ring,3,NeutType) = &
                        &  HC_WTSource ( Starting_Index, Launch_Ring,3,NeutType) + Sput_Weight
                   !Write (Output_Unit_Scratch,*) "SOURCING3:", HC_WTSource ( Starting_Index, Launch_Ring,3,NeutType),NeutType,Cur_HC_Spec	

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
          Call HC_Prompt_Deposition (		&
               &	NProd,				&
               &	Current_Cell,			&
               &	Current_Ring,			&
               &	Current_R,			&
               &	Current_Z,			&
               &	Current_S,			&
               &	Current_Cross,			&
               &	Current_Theta,			&
               &	Cur_HC_Spec,		&
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
               &	Tries_Without_Change,		&
               &	Neutral_Time_Step_Count,	&
               &	Ion_Time_Step_Count,		&
               &	Eq_Total_Ion_Time_Steps,	& ! CIST
               &	Seed,				&
               &	NRand,				&
               &	IFate,				&
               &        hc_v,                           & ! hc_velocity data
               &	 Debug_HC_Prompt)				

          !End If

          !---------------------------------------------------------------------------------------------
          !
          !------------------------ NEUTRAL TO NEUTRAL -------------------------------------------------
          !
          !---------------------------------------------------------------------------------------------
          ! jdemod - placed this code in an ElseIf block since both this and the previous block should
          !          not both be executed. Change of states resulting from prompt deposition might possibly
          !          trigger this condition
          ! NEUTRAL TO NEUTRAL. Check for ionization and recombination and alter grid coordinates accordingly.
       ElseIf((Last_HC_Species.ne.Cur_HC_Spec).and.(Get_HC_Charge(Last_HC_Species).eq.0).and.(Get_HC_Charge(Cur_HC_Spec).eq.0))Then

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

             Call Neut_Neut_FreeSpace_Location (Reactions_Possible (Counter),Followed_State,&
                  & Current_R,Current_Z,Current_Cell,&
                  & Current_Ring)
             Call Neut_Neut_FreeSpace_Angle (Reactions_Possible (Counter),Followed_State,&
                  &Current_Angle,Azimuthal_Angle,Seed,NRand)
             !Write (Output_Unit_Scratch,*) "ANGLE",Old_Angle,Current_Angle,Reactions_Possible (Counter),Followed_State

             Call Neut_Neut_FreeSpace_Energy (Reactions_Possible (Counter),Followed_State,&
                  &H_Isotope_Composition,Current_Cell,&
                  &Current_Ring,Current_Velocity,HC_Temperature,Kin_Energy_Added)
             Call Neut_Neut_FreeSpace_Velocity (Reactions_Possible (Counter),Followed_State,&
                  &H_Isotope_Composition,Current_Angle,&
                  &Azimuthal_Angle,Current_Velocity,Current_Velocity_In_R,Current_Velocity_In_Z,Seed,NRand)

          elseif (hc_kinetics_opt.eq.1) then 

             call update_hc_kinetics(hc_v,nn_reaction,current_cell,current_ring,  & 
                  & cur_hc_spec,reactions_possible(counter),followed_state,h_isotope_composition,  &
                  & current_r,current_z,current_s, current_theta,current_cross,current_angle,  &
                  & current_velocity_in_s,current_velocity_in_r,current_velocity_in_z, &
                  & current_velocity,kin_energy_added,s_start, &
                  & neutral_time_step,ion_time_step,sput_weight,iprod)

          endif

          ! jdemod - not sure why this is here - on neutral to neutral transitions the value of S is not 
          !          significant - initialization to zero might only be useful to avoid use of 
          !          uninitialized data but that data should not be used anyway
          ! Note:  Must set current S to zero to account for an initial neutral to neutral step.
          Current_S = 0.0

          ! Print hydrocarbon evolution data to file if option H31 selected in input file.
          If (HC_Evolve_Print_Option .eq. 1) Then
             ! Print option activated.
             Time_Tot = Eq_Total_Ion_Time_Steps *  Ion_Time_Step
             Write (Output_Unit_Evolve,9301) Time_Tot,Neutral_Time_Step_Count,Ion_Time_Step_Count,Last_HC_Species,&
                  & Last_H_Isotope_Composition(1),Last_H_Isotope_Composition(2), &
                  & Last_H_Isotope_Composition(3),Cur_HC_Spec,H_Isotope_Composition(1),H_Isotope_Composition(2),&
                  & H_Isotope_Composition(3),Current_Cell,Current_Ring, &
                  & Old_R,Current_R,Old_Z,Current_Z,Old_S,Current_S,Old_Cross,Current_Cross,Old_Angle,Current_Angle,Old_Velocity,&
                  & Current_Velocity,Old_Temp,HC_Temperature
9301         Format (E10.4,F10.1,F10.1,I10,2X,I2,1X,I2,1X,I2,I10,2X,I2,1X,I2,1X,I2,I10,I10,F10.6,F10.6,F10.6,F10.6,F10.6,F10.6,&
                  &F10.4,F10.4,F10.4,F10.4,F10.3,F10.3,F10.3,F10.3,F10.3,F10.3)
          End If
          ! New neutral will not be C+ so return to hc_follow code.
          Return
       End If

       If (Cur_HC_Spec .eq. 1) Then ! Reached "C+".
          ! C+ production has occured.  Store particle details in arrays and totals.

          If (NeutType .eq. 2 .or. NeutType .eq. 5) Then
             HC_ChemIzs (Current_Cell,Current_Ring,Cur_HC_Spec) =  HC_ChemIzs (Current_Cell,Current_Ring,Cur_HC_Spec) + Sput_Weight
          End If

          NatIZ = NatIZ + 1
          HC_Num_Fragments_Reach_CIon (Launch_Reg) =  HC_Num_Fragments_Reach_CIon (Launch_Reg) + Sput_Weight ! HC routine variable in HC_Launch_Diag_Table, same as RatIZ / SatIZ.
          !write (0,*) "adding one to neutral ionization:", HC_Num_Fragments_Reach_CIon (Launch_Reg),Launch_Reg
          Call pxatizs (LatIZ + NatIZ - 1,Current_R)
          Call pyatizs (LatIZ + NatIZ - 1,Current_Z)
          Call pkatizs (LatIZ + NatIZ - 1, Local_K)
          Call psatizs (LatIZ + NatIZ - 1,MIN (Current_S, gksmaxs (Current_Ring) - Current_S))
          Call ptemtizs (LatIZ + NatIZ - 1,HC_Temperature)
          Call psnews (LatIZ + NatIZ - 1,Sput_Weight)
          Call pcistizs (LatIZ + NatIZ - 1,REAL (Eq_Total_Ion_Time_Steps))

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

          HC_LaunchDat (LatIZ + NatIZ - 1,3) =  HC_LaunchDat (LatIZ + NatIZ - 1,3) + REAL (Reflection_Counter)
          HC_LaunchDat (LatIZ + NatIZ - 1,2) =  HC_LaunchDat (LatIZ + NatIZ - 1,3) +  HC_Launchdat (IProd,2)

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
          ! Launch_Reg    = 1 for target 1	ik > nks(ir)/2
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

          ! jdemod - this needs to be set or the code can crash if launch_target is not initialized
          launch_target = 0

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

             If ( Far_Periphery_Recycle_Opt .eq. 1 .and.  HC_Launchdat (IProd,2) .eq. 1.0) Then
                ifp = 2
             Else
                ifp = 1
             End If

             If (( Neutral_Reflection_Opt .eq. 1 .or.  Neutral_Reflection_Opt .eq. 2) .and. Reflection_Counter .gt. 0) Then
                irflct = 2
             Else
                irflct = 1
             End If

             !write(0,*) 'bnds:',isol,launch_target,ifp,irflct,launch_reg

             HC_LIonizDat (isol,Launch_Target,ifp,irflct,1) =  HC_LIonizDat (isol,Launch_Target,ifp,irflct,1) + Sput_Weight
             HC_LIonizDat(isol,Launch_Target,ifp,irflct,2)=HC_LIonizDat(isol,Launch_Target,ifp,irflct,2)+Sput_Weight*Current_R
             HC_LIonizDat(isol,Launch_Target,ifp,irflct,3)=HC_LIonizDat(isol,Launch_Target,ifp,irflct,3)+Sput_Weight*Current_Z
             HC_LIonizDat(isol,Launch_Target,ifp,irflct,4)=HC_LIonizDat(isol,Launch_Target,ifp,irflct,4)+Sput_Weight*Local_K
             HC_LIonizDat (isol,Launch_Target,ifp,irflct,5) =  HC_LIonizDat (isol,Launch_Target,ifp,irflct,5) + Sput_Weight * MIN (&
                  &Current_S,gksmaxs (Current_Ring) - Current_S)
          End If

          IFate = 13 ! HC Neutral Reduced to C+.

          Return

       End If

    End If

  End Subroutine HC_Trans_Inside_Neutral



End Module HC_Inside_Neutral
