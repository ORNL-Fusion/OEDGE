! -*-Mode:f90-*-
! Out_Utils.f90
! Contains generic free-format output utilities including:
! -HC_Trans_Prob: Calculates transition probability and saves to comhc.
!
! Adam McLean under Prof. Peter Stangeby
! Research Assistant: David Elder
! August, 2002
 
Module Out_Utils
 
  Use HC_Init_Lib_Data

  ! Every good Fortran program has...
  Implicit None	
 
Contains
 
  !
  ! jdemod - these functions are already available from HC_INIT_LIB_DATA
  ! 
!  Integer Function HC_Species_Ident (Species)
!    ! Return identification number for a hydrocarbon species.
! 
!    ! Required modules.
!    Use HC_Init_Lib_Data ! Gain access to carbon stack declaration.
! 
!    ! Every good Fortran program has...		
!    Implicit None
! 
!    ! Declare input/output variables.
!    Character (Len=*), Intent (In) :: Species
! 
!    ! Declare local variables.
!    Integer :: i,j
! 
!    ! Find species name in hydrocarbon state table.
!    Do j = 1, SIZE (HC_State_Table), 1
!       If (HC_State_Table (j) % State_Name .eq. Species) Then
!          ! Found the proper hydrocarbon name.  Assign its state number to the output array.
!          HC_Species_Ident = HC_State_Table (j) % State_Number
!          ! Finish looping around embedded Do.
!          Exit
!       Else
!          ! Check to see if we're at the end of the list of supported hydrocarbons.
!          If (j .eq. SIZE (HC_State_Table)) Then
!             Write (Output_Unit_HC_Alert,*) "Error in HC_Species_Ident: Unsupported hydrocarbon specified: ",Species
!             Write (Output_Unit_HC_Alert,*) "Program stopping."
!             Stop
!          Else
!             ! Try the next one.
!             Cycle
!          End If
!       End If
!    End Do
! 
!  End Function HC_Species_Ident
! 
!  Character (Len=10) Function HC_Ident_Species (Ident)
!    ! Return species names upon input of a hydrocarbon identification numbers.
! 
!    ! Required modules.
!    Use HC_Init_Lib_Data ! Gain access to carbon stack declaration.
! 
!    ! Every good Fortran 90 program has...		
!    Implicit None
! 
!    ! Declare input/output variables.
!    Integer, Intent (In) :: Ident
! 
!    ! Declare local variables.
!    Integer :: i
! 
!    ! Check to see if number is larger than supported hydrocarbons.
!    If (Ident .lt. 1 .or. ident .gt. Size (HC_State_Table)) Then
!       Write (Output_Unit_HC_Alert,*) "Error in HC_Ident_Species: Unsupported hydrocarbon specified. IDENT =",IDENT
!       Write (Output_Unit_HC_Alert,*) "Program stopping"
!       Stop
!    End If
!    HC_Ident_Species = HC_State_Table (Ident) % State_Name
! 
!  End Function HC_Ident_Species
! 
!  Integer Function Get_HC_Charge (State_Number)
! 
!    ! Required modules.
!    Use HC_Init_Lib_Data ! Gain access to carbon stack declaration.
! 
!    Implicit None
!    Integer, Intent (IN) :: State_Number
!    Integer :: Length
! 
!    ! Find de-spaced length of State_Name.
!    Length=0
!    Do
!       If (HC_State_Table (State_Number) % State_Name (Length+1:Length+1) .eq. "") Then
!          Exit
!       End If
!       Length = Length + 1
!    End Do
! 
!    ! Calculate molecular charge.
!    If (HC_State_Table (State_Number) % State_Name(Length:Length) .eq. "+") Then
!       Get_HC_Charge = 1
!       If (HC_State_Table (State_Number) % State_Name(Length - 1:Length - 1) .eq. "+") Then
!          Get_HC_Charge = 2
!       End If
!    Else
!       ! No positive charge.
!       Get_HC_Charge = 0
!    End If
! 
!  End Function Get_HC_Charge
! 

  Subroutine Fill_HC_Trans_Prob_Table (Start_HC_Species,End_HC_Species,HC_TimeStep)
    ! Fills the HC_Trans_Prob table
 
    Use ComHC ! Includes HC_Trans_Prob array.
    !		Use HC_Init_DIV_Data ! Load cell and global properties data structures.
    Use HC_Init_Out_Data ! HC diagnostics.
    Use HC_Get ! Get density, temperature functions.
    Use HC_NewSt ! Subroutines for hydrocarbon transitions.
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Declare input.output variables.
    Integer, Intent (In) :: Start_HC_Species
    Integer, Intent (In) :: End_HC_Species
    Real, Intent (Out) :: HC_TimeStep
 
    ! Local variables.
    Integer :: ik
    Integer :: ir
    Integer :: Reaction_Number
    Integer :: Multi_State
 
    ! Variables for hydrocarbon change of state.
    ! jdemod - only load states_possible after state selected 
    !          not needed at all for this application in OUT
    !Integer, Dimension (Highest_Carbon_Content) :: States_Possible
    Integer, Dimension (Max_Reactions_Per_State) :: Reactions_Possible
    Real, Dimension (Max_Reactions_Per_State) :: Total_Probabilities
    Real, Dimension (Max_Reactions_Per_State) :: NS_Probabilities
    Real, Dimension (Max_Reactions_Per_State) :: Reaction_Probabilities
    Real :: Total_Probability
    Character (len=6) :: Data_Type = "SigmaV"

    ! jdemod - added number_reactions to hc_new_state call
    integer :: number_reactions
    integer :: in

    ! jdemod - not needed - there were only 5 variables being accessed from this source - 2 of which
    !          were already available by another name in comhc - while the other 3 were more easily
    !          implemented as get functions
    ! Initialize global data.
    ! Call Initialize_Global_Prop_Data ()
    ! Call Initialize_Global_Geom_Data ()
 
    If (Get_HC_Charge (Start_HC_Species) .eq. 0) Then
       HC_TimeStep =  gfsrate()
    Else
       HC_TimeStep =  gqtim()
    End If
 
    ! Load up HC transition probability table.
    Do ir = 1,  maxnrs
       Do ik = 1,  maxnks
 
          ! Zero out comhc array from earlier possible run.
          HC_Trans_Prob (ik, ir) = 0.0
 
          ! Zero out necessary arrays for maximum reactions per state.
          !States_Possible = 0
          Reactions_Possible = 0
          Total_Probabilities = 0.0
          NS_Probabilities = 0.0
          Reaction_Probabilities = 0.0
 
          ! Calculate transition probabilities.
          ! Note that E&L reaction rates are in cm^3/s.  Therefore, DIVIMP plasma
          ! densities must be reduced by 1.0E6 (from 1/m^3) for use.
          ! jdemod - changed to remove states_possible and add number_reactions
          Call HC_New_State (ik,ir,Start_HC_Species,Data_Type,gknbs (ik, ir) / 1.0E6,gktebs(ik,ir),gktibs(ik,ir), &
               &  gplasma_mass(),HC_TimeStep,Reactions_Possible,Reaction_Probabilities,NS_Probabilities,&
               & Total_Probability,number_reactions)
 
          ! Record desired probability for final HC species.
          If (End_HC_Species .eq. -1) Then
             ! Record total probability.
             HC_Trans_Prob (ik, ir) = Total_Probability * 100
          Else
             ! Record individual destination probabilty.
             ! Note possibility of multiple reactions to same final state.
             ! jdemod - use number_reactions
             !Do Reaction_Number = 1, Max_Reactions_Per_State
             Do Reaction_Number = 1, number_reactions
                Do Multi_State = 1, Highest_Carbon_Content
                   ! jdemod - changed usage of states_possible - not needed now
                   if (hc_state_transform_table(reactions_possible(reaction_number))%end_c_states(multi_state)&
                       &.eq.end_hc_species) then 
                   !If (States_Possible (Reaction_Number,Multi_State) .eq. End_HC_Species) Then
                      ! Add this probability to the total.
                      HC_Trans_Prob (ik, ir) =  HC_Trans_Prob (ik, ir) + Reaction_Probabilities (Reaction_Number) * 100
                   End If
                End Do
             End Do
          End If
 
       End Do
    End Do
 
  End Subroutine Fill_HC_Trans_Prob_Table
 
  Subroutine Fill_HC_XCH (X_Method)
    ! Fills the X_CH table in HC_OUT_Data_Table.
 
    Use HC_Init_Out_Data ! HC diagnostics.
    !		Use HC_Init_DIV_Data ! Load cell and global properties data structures.
    Use HC_Get ! Get density, temperature functions.
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Declare input.output variables.
    Character (Len=*), Intent (In) :: X_Method
 
    ! Local variables.
    Integer :: ik
    Integer :: ir
    Real :: E_CD
    Real :: f_CD
    Real :: Local_E_Temp
    Real :: Local_E_Density
 
    E_CD = 12400.0/4308.0
    f_CD = 0.00516
 
    If (X_Method .eq. "NRL") Then
 
       Do ir = 1,  Maxnrs
          Do ik = 1,  Maxnks
 
             ! Zero out comhc array from earlier possible run.
             X_CH (ik,ir) = 0.0
 
             Local_E_Temp = gktebs (ik,ir)
             Local_E_Density = gknbs (ik,ir)
 
             ! NRL calculation.
             X_CH (ik,ir) = 1.0E-6 * 16.0E-6 * f_cd * EXP(-E_CD/Local_E_Temp) / (E_CD * SQRT(Local_E_Temp))
             ! Mewe correction.
             X_CH (ik,ir) =  X_CH (ik,ir) * (0.6+0.28*(LOG(E_CD/Local_E_Temp+1/(E_CD/Local_E_Temp))-0.4/(1+E_CD/Local_E_Temp)**2))
 
          End Do
       End Do
 
    ElseIf (X_Method .eq. "NAUJOKS") Then
 
       Do ir = 1,  Maxnrs
          Do ik = 1,  maxnks
 
             ! Zero out comhc array from earlier possible run.
             X_CH (ik,ir) = 0.0
 
             ! Naujoks calculation.
             X_CH (ik,ir) = 1.45E-14 * EXP(-3.19/gktebs(ik,ir))/(gktebs(ik,ir)**0.28)
 
          End Do
       End Do
 
    Else
       Write (0,*) "Option for Fill_HC_XCH not supported:",X_Method
       Stop
    End If
 
  End Subroutine Fill_HC_XCH
 
  Subroutine Fill_HC_XC2 (X_Method)
    ! Fills the X_C2 table in HC_OUT_Data_Table.
 
    Use HC_Init_Out_Data ! HC diagnostics.
    !		Use HC_Init_DIV_Data ! Load cell and global properties data structures.
    Use HC_Get ! Get density, temperature functions.
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Declare input.output variables.
    Character (Len=*), Intent (In) :: X_Method
 
    ! Local variables.
    Integer :: ik
    Integer :: ir
    Real :: E_CD
    Real :: f_CD
    Real :: Local_E_Temp
    Real :: Local_E_Density
 
    E_CD = 12400.0/5160.0
    f_CD = 0.03*0.73
 
    If (X_Method .eq. "NRL") Then
 
       Do ir = 1,  Maxnrs
          Do ik = 1,  Maxnks
 
             ! Zero out comhc array from earlier possible run.
             X_C2 (ik,ir) = 0.0
 
             Local_E_Temp = gktebs (ik,ir)
             Local_E_Density = gknbs (ik,ir)
 
             ! NRL calculation.
             X_C2 (ik,ir) = 1.0E-6 * 16.0E-6 * f_cd * EXP(-E_CD/Local_E_Temp) / (E_CD * SQRT(Local_E_Temp))
             ! Mewe correction.
             X_C2 (ik,ir) =  X_C2 (ik,ir) * (0.6+0.28*(LOG(E_CD/Local_E_Temp+1/(E_CD/Local_E_Temp))-0.4/(1+E_CD/Local_E_Temp)**2))
 
          End Do
       End Do
    Else
       Write (0,*) "Option for Fill_HC_XCH not supported:",X_Method
       Stop
    End If
 
  End Subroutine Fill_HC_XC2
 
End Module Out_Utils
