! -*-Mode:f90-*-
! HC_Start.f90
! Contains routines run before and after hydrocarbon code is
! run within DIVIMP.  Code here is only called and run once
! either at the beginning or ending of hydrocarbon code execution.
!
! Adam McLean under Prof. Peter Stangeby
! Research Assistant: David Elder
! March, 2003
 
Module HC_Start
 
  ! Every good Fortran 90 program has...
  Implicit None	
 
Contains
 
  ! Hydrocarbon setup.  Routine to perform all initialization that
  ! needs to be done only once throughout the entire case.  Note
  ! that this is done at the very beginning of hydrocarbon following.
  Subroutine HC_Begin ()
 
    ! Required modules.
    Use ComHC ! Contains DIVIMP input file hydrocarbon-related options.
    Use HC_Init_DIV_Data ! Contains global initialization routines.
    Use HC_Init_DIV_Diag ! HC diagnostic data.
    Use HC_Init_Lib_Data ! Hydrocarbon module initialization.
    Use HC_LdDta ! Access hydrocarbon storage data.
    use hc_put   ! Put data to DIVIMP arrays - included here to initialize error message channel 
    use hc_get   ! Get data from DIVIMP arrays - included here to initialize error message channel 
    use hc_diag_data ! Module to contain HC code diagnostic monitoring variables
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables.
    Integer :: IOStat
 
    ! Check if HCs are enabled.
    If (hc_follow_option .eq. 1) Then
 
       ! Perform initial file handling.
       ! Open hydrocarbon error file if HCs are enabled.
       Open (Unit = Output_Unit_HC_Alert, Action = "Write", File = "hc_alert.txt", IOSTat = IOStat)
       Write (Output_Unit_HC_Alert,'(A)') "Hydrocarbon Alert File"
       If (IOStat .gt. 0) Then
          Write (0,*) "Error in HC_Start:  Could not create HC Alert File:",IOStat
          Write (0,*) "Program stopping."
          Stop
       End If
       !
       ! jdemod 
       ! Initialize the error message channel in the put and get modules - this was the only value in comhc used
       ! in these modules - however since the modules used both comhc at a global level and params in each routine
       ! various error messages arose due to namespace collisions or in particlar multiple save statements for the same
       ! variable.
       !
       call hc_get_init_output(output_unit_hc_alert)
       call hc_put_init_output(output_unit_hc_alert)

 
       ! Open hydrocarbon data file if HCs are enabled.
       Open (Unit = Output_Unit_HC_Data, Action = "Write", File = "hc_output.txt", IOSTat = IOStat)
       Write (Output_Unit_HC_Data,'(A)') "Hydrocarbon Data File"
       If (IOStat .gt. 0) Then
          Write (Output_Unit_HC_Alert,*) "Error in HC_Start:  Could not create HC Data File:",IOStat
       End If
 
       ! Open hydrocarbon scratch file if HCs are enabled.
       Open (Unit = Output_Unit_Scratch, Action = "Write", File = "hc_scratch.txt", IOSTat = IOStat)
       Write (Output_Unit_Scratch,'(A)') "Hydrocarbon Scratch File"
       If (IOStat .gt. 0) Then
          Write (Output_Unit_HC_Alert,*) "Error in HC_Start:  Could not create HC Scratch File:",IOStat
       End If
 
       ! Open C+ initial position file if HCs are enabled.
       Open (Unit = Output_Unit_Cpos_Pos, Action = "Write", File = "hc_cpos_stat.txt", IOSTat = IOStat)
       Write (Output_Unit_Cpos_Pos,'(A)') "Initial C+ Position Data File"
       If (IOStat .gt. 0) Then
          Write (Output_Unit_HC_Alert,*) "Error in HC_Start:  Could not create C+ initial position file:",IOStat
       End If
       Write (Output_Unit_Cpos_Pos,*) " "
       Write (Output_Unit_Cpos_Pos,9190) "Event Time","Neut","Ion","Init HC","H/D/T B","Dest HC","H/D/T A","Cell","Ring","R &
&Before","R After","Z Before", &
            & "Z After","S Before","S After","Cross B","Cross A","Angle B","Angle A","Vel B","Vel A","Temp B","Temp A"
9190   Format (40A10)
 
       ! If additional output is requested, open files.
       If (hc_evolve_print_option .eq. 1) Then
          Open (Unit = Output_Unit_Evolve, Action = "Write", File ="hc_evolution.txt", IOSTat = IOStat)
          Write (Output_Unit_Evolve,'(A)') "Hydrocarbon Evolution Timestep Data"
9200      Format (40A10)
          If (IOStat .gt. 0) Then
             Write (Output_Unit_HC_Alert,*) "Error in HC_Start:  Could not create HC Evolution Timestep Data File:",IOStat
          End If
       End If


       If (hc_coord_print_option .eq. 1) Then
          Open (Unit = Output_Unit_Location, Action = "Write", File ="hc_location.txt", IOSTat = IOStat)
          Write (Output_Unit_Location,'(A)') "Hydrocarbon Location Timestep Data"
          If (IOStat .gt. 0) Then
             Write (Output_Unit_HC_Alert,*) "Error in HC_Start:  Could not create HC Location Timestep Data File:",IOStat
          End If
       End If
 
       ! Call global data setup routines.
       Call Initialize_Global_Geom_Data ()
       Call Initialize_Global_Prop_Data ()
       Call Initialize_DIVIMP_Options_Table ()
       Call Initialize_DIVIMP_HC_Data_Table ()
       Call Initialize_Misc_Data_Table ()
       ! jdemod - initialize the change of state routines - pass in the particle masses
       !Call Initialize_HC_Data (back_plasma_ion_mass,impurity_ion_mass)
       ! H isotope mass associated with hydrocarbons is now specified in input as input_HC_H_mass in the 
       ! comhc module - as a result it is assigned directly in the initialize_HC_data routine
       Call Initialize_HC_Data (impurity_ion_mass)
       Call Load_HC_Data ()
       Call Initialize_HC_Vessel_Interact () ! Setup reflection and sticking coefficient arrays.
 
       ! Set all launch diagnostics values for each target and sums to 0.0.
       ! All these are used in diagnostic printing after all HC particles have been followed.
       Call Initialize_HC_Launch_Diag_Table ()
       Call Initialize_Transport_Diag_Table ()
       Call Initialize_Evolve_Diag_Table ()
       Call Initialize_Rerelease_Diag_Table ()
       Call Initialize_VesselInt_Diag_Table ()
       Call Initialize_Death_Diag_Table ()
 
       ! Set debug status.
       Debug_HC_Ion = .false.
       Debug_HC_Neutral = .false.
 
       ! Set walk count in ComHC for HC particle tracks.
       HC_Walk_Count = 1 ! IW
 
       ! Set transport kinetics option found in ComHC.
       HC_Energy_Calc = 1
 
       ! Set self sputtering option found in ComHC.
       HC_Self_Sputter = 0
 
       ! Ion reactant total probability for reaction multiplier found in ComHC.
       ! 1.0 - no change in probabilities, otherwise a linear change in total reaction probability from an ionized reactant.
       HC_Ion_Reaction_Mult = 1.0
 
       ! Set flag to include reaction kinetics in post-transition reactant energy.
       ! 0 - no added energy to transition products, 1 - energy added as prescribed in available transition database(s).
       ! jdemod - this value has been moved to the input file as an optional input value - its default value is 1
       ! HC_Reaction_Kinetics = 1
 
       ! Set flag to indicate us of imrpvoed 3D reaction kinetics including
       ! influence on particle direction, velocity and energy.
       ! jdemod - if this operates as it sounds it should be implemented as a different hc_reaction_kinetics option value - not a 
       !          new option
       HC_Improved_Kinetics = 0
 


       ! Initialize the HC diagnostic data routines
       call init_hc_diag_data(number_hc_species)


    End If
 
  End Subroutine HC_Begin
 
  Subroutine HC_End (NProd,LProd,SExit,SMain,SHC,NeuTime,Charge_States)
 
    ! Required modules.
    Use HC_Output ! Contains results printing routine.
    Use HC_Put ! Access routines to DIVIMP data.
    Use HC_Utilities ! Contains Record_DIVIMP_Data.
    Use HC_WBC_Comp ! Contains Print_WBC_Table.
    Use ComHC ! Contains DIVIMP input file hydrocarbon-related options.
    use hc_diag_data ! Module to collect diagnostic data
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables
    Integer, Intent (In) :: NProd
    Integer, Intent (In) :: LProd
    Real, Intent (In) :: SExit
    Real, Intent (In) :: SMain
    Real, Intent (In) :: SHC
    Real, Intent (In) :: NeuTime
    Integer, Intent (In) :: Charge_States
 
    ! Complete HCs.  Write applicable impurity data back to DIVIMP common blocks.
    ! Note, do this before call to HC_Print_Output as some processing is done at printing
    ! that should not be done twice (also done in DIVIMP).
    Call Record_DIVIMP_Data ()
 
    ! Produce print-out of all recorded data relevant to hydrocarbons.
    Call HC_Print_Output (NProd,LProd,SExit,SMain,SHC,NeuTime)
 
    ! If output was requested, close files.
    If (hc_follow_option .eq. 1) Then
       Close (Unit = Output_Unit_HC_Data)
       Close (Unit = Output_Unit_Scratch)
       Close (Unit = Output_Unit_HC_Alert)
       Close (Unit = Output_Unit_Cpos_Pos)
       If (hc_evolve_print_option .eq. 1) Then
          Close (Unit = Output_Unit_Evolve)
       endif
       If (hc_coord_print_option .eq. 1) Then
          Close (Unit = Output_Unit_Location)
       End If
    End If
 
    ! Append remainder of HC_Case.out WBC table.
    If (HC_WBC_Comp_Option .gt. 0) Then
       Call Print_WBC_Comp (Charge_States)
    End If

    !
    ! Clean up HC diagnostic data
    !
    call end_hc_diag_data

  End Subroutine HC_End
 
End Module HC_Start
