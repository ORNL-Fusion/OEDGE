! -*-Mode:f90-*-
! hc_new_state.f90 (current_state,density,temperature,target_energy,hydrogen_mass,timestep,output_flag,new_state)
! Hydrocarbon Reaction Program
!
! Adam McLean under Prof. Peter Stangeby
! Research Assistant: David Elder
! October, 1999
!
! Note: If this routine is called multiple times, be
!       sure to zero out input arrays between calls.
 
Module HC_NewSt
 
  ! Gain access to external data and subroutines.
  Use HC_Init_Lib_Data
  ! jdemod - put all the probability utility routines into one module - no need for 5 files
  use hc_probability
  !Use HC_RnLup
  !Use HC_RnPrb
  !Use HC_NSPrb
  !Use HC_StLup
  !Use HC_TlPrb
 
  ! Every good Fortran program should have...
  Implicit None
 
Contains
 
  Subroutine HC_New_State (Current_Cell,Current_Ring,Current_State,Data_Type,Plasma_Density,Plasma_Temperature,Target_Energy, &
       & Hydrogen_Mass,Timestep,Local_Reactions_Possible,Reaction_Probabilities,NS_Probabilities,Total_Probability,&
       & number_reactions)
    !Subroutine HC_New_State (Current_Cell,Current_Ring,Current_State,Data_Type,Plasma_Density,Plasma_Temperature,Target_Energy, &
    !     & Hydrogen_Mass,Timestep,Local_Reactions_Possible,States_Possible,Reaction_Probabilities,NS_Probabilities,Total_Probability,&
    !     & number_reactions)
    !
    ! jdemod - NOTE, WARNING!
    !          The plasma density currently passed into this routine is in units of cm-3 (not mks units as are used in the rest of the code)
    !          This is because both the EL coefficients and the JR tabulated data give sigma-V in terms of cm-2s-1
    !
    ! jdemod - interface changes 
    !          - remove states_possible - it is redundant once the new state is selected - only useful if hc_init_lib_data
    !          is not accessible from the transport routines - which is not the case - and even if it was - it would be more efficient to 
    !          add a routine here to return just those states leading from the selected process.
    !          - add number of reactions found - to make search code more efficient instead of just running through the entire array
    !          - change NS probability to be cumulative 
    !

    Implicit None

    ! Define externally important variables
    integer, intent(out) :: number_reactions
    Integer, Intent (IN) :: Current_Cell
    Integer, Intent (IN) :: Current_Ring
    Integer, Intent (IN) :: Current_State
    Character (len=*), Intent (IN) :: Data_Type
    Real, Intent (IN) :: Plasma_Density,Plasma_Temperature,Target_Energy,Hydrogen_Mass,Timestep
    Integer, Dimension (:), Intent (OUT) ::  Local_Reactions_Possible
    ! jdemod - remove states_possible
    !Integer, Dimension (:,:), Intent (Out) :: States_Possible
    Real, Dimension (:), Intent (OUT) :: Reaction_Probabilities,NS_Probabilities
    Real, Intent (OUT) :: Total_Probability

    integer,save :: probabilities_printed(number_hc_reactions) 
    real,save :: reaction_probability_data(number_hc_reactions)
    integer :: in,reaction_index
    integer,save :: init_newst_print = 1
    logical,save :: debug_newst = .false.

    if (debug_newst.and.(init_newst_print.eq.1)) then 
       probabilities_printed=0
       reaction_probability_data = 0.0
       init_newst_print = 2
       write(0,*) 'Initialize NEWST PRINT:',init_newst_print
    endif


    ! Fill array with possible transition states upon programmed reactions.
    Call HC_Reaction_Lookup (Current_State,Local_Reactions_Possible,number_reactions)

    ! Check plasma temperature is within the range of applicability for each reaction possible.
    if (hc_evolution_model_primary.eq.1) then 
       Call HC_Bounds_Check_Ehrhardt_Langer (Current_Cell,Current_Ring,Plasma_Temperature,Local_Reactions_Possible,Data_Type)
    endif

    ! jdemod - remove this - only get the states after reaction process has been selected - too much work otherwise
    ! Get list of final states possible 
    !Call HC_State_Table_Lookup (Local_Reactions_Possible,States_Possible)

    ! Calculate associated probability with each possible transition.
    Call HC_Reaction_Probability (Local_Reactions_Possible,Data_Type,Plasma_Density,Plasma_Temperature,Target_Energy,Hydrogen_Mass,&
         &Timestep,Reaction_Probabilities)

    ! jdemod - change NS_probability to cumulative instead of individual
    ! Find normalized state probability for individual transition given a transition.
    Call HC_NS_Probability (Reaction_Probabilities, NS_Probabilities,number_reactions)

    ! Find total probability for transition.
    Call HC_Total_Probability (Reaction_Probabilities,Total_Probability,number_reactions)

    !Total_Probability = Total_Probability / 3.0
    !If (Total_Probability .gt. 1.0) Then
    !   Total_Probability = 1.0
    !End If

    !if (Plasma_Temperature .ne. 20.0 .or. Plasma_Density .ne. 2.5E13) Then
    !	write (0,*) "HCing",Plasma_Temperature,Plasma_Density,Current_Cell,Current_Ring,Current_State
    !endif

    !write (0,*) "DIVIMP-HC:",Current_State,Local_Reactions_Possible,States_Possible,Plasma_Density,&
    !                     &Plasma_Temperature,Target_Energy,Hydrogen_Mass,Timestep,Reaction_Probabilities,NS_Probabilities
    !stop


    ! write out probabilities for the reaction

    if (debug_newst) then 

       do in = 1,number_reactions
          reaction_index = local_reactions_possible(in)
          if (probabilities_printed(reaction_index).eq.0) then 
             probabilities_printed(reaction_index) = 1
             reaction_probability_data(reaction_index) = reaction_probabilities(in)

             write(output_unit_hc_alert,'(a,2i8,a,a20,i6,10(1x,g12.5))') 'REACTION    :',current_state,reaction_index,':',&
                  &hc_state_transform_table(reaction_index)%reaction_desc(1:20),&
                  &in,plasma_density,plasma_temperature,reaction_probabilities(in),total_probability


          elseif (abs(reaction_probabilities(in)-reaction_probability_data(reaction_index)).gt.1e-8) then 

             write(output_unit_hc_alert,'(a,2i8,a,a20,i6,10(1x,g12.5))') 'REACTION ERR:',current_state,reaction_index,':',&
                  &hc_state_transform_table(reaction_index)%reaction_desc(1:20),&
                  &in,plasma_density,plasma_temperature,reaction_probabilities(in),total_probability


          endif


       end do

    endif

  End Subroutine HC_New_State
 
  Subroutine HC_Bounds_Check_Ehrhardt_Langer (Current_Cell,Current_Ring,Plasma_Temperature,Local_Reactions_Possible,Data_Type)
 
    ! Gain access to external data and subroutines.
    Use HC_Init_Lib_Data
 
    Implicit None
 
    ! Define externally important variables
    Integer, Intent (IN) :: Current_Cell
    Integer, Intent (IN) :: Current_Ring
    Real, Intent (IN) :: Plasma_Temperature
    Integer, Dimension (:), Intent (IN) ::  Local_Reactions_Possible
    Character (len=*), Intent (IN) :: Data_Type
 
    ! Define local variables.
    Integer :: i
    Logical :: Check_Temp
 
    ! Assign local variables.
    i = 1
    Check_Temp = .False.
    !Check_Temp = .True.
 
    ! Check to see if a temperature check is desired.
    If (.not. Check_Temp) Then
       Return
    End If
 
    ! Check that plasma temperature is higher than the minimum for each applicable reaction.
    If (Plasma_Temperature .ne. 0.0) Then
       Do While (Local_Reactions_Possible (i) .ne. 0)
 
          If (Data_Type .eq. 'Sigma') Then
             If (Plasma_Temperature .lt. HC_State_Transform_Table (i) % Sigma_Tmin_Limit) Then
                ! Print warning.
                Write (Output_Unit_HC_Alert,100) "Warning: Plasma temperature",Plasma_Temperature,"eV in Cell",Current_Cell,",&
                     &Ring",Current_Ring, &
                     & "is less than applicable for Sigma:",HC_State_Transform_Table (i) % Sigma_Tmin_Limit,"eV for reaction:",&
                     &Local_Reactions_Possible (i)
100             Format (A,1X,F6.3,1X,A,1X,I2,1X,A,1X,I2,1X,A,F6.3,1X,A,1X,I2)
             End If
          Else If (Data_Type .eq. 'SigmaV') Then
             If (Plasma_Temperature .lt. HC_State_Transform_Table (i) % SigmaV_Tmin_Limit) Then
                ! Print warning.
                Write (Output_Unit_HC_Alert,101) "Warning: Plasma temperature",Plasma_Temperature,"eV in Cell",Current_Cell,",&
                     &Ring",Current_Ring, &
                     & "is less than applicable for SigmaV:",HC_State_Transform_Table (i) % SigmaV_Tmin_Limit,"eV for reaction:",&
                     &Local_Reactions_Possible (i)
101             Format (A,1X,F6.3,1X,A,1X,I2,1X,A,1X,I2,1X,A,F6.3,1X,A,1X,I2)
             End If
          Else
             Write (Output_Unit_HC_Alert,*) 'Data type not defined in call to HC_Interp: ', Data_Type
          End If
 
          ! Also check upper limit.
          If (Plasma_Temperature .gt. 2000) Then
             ! Temperature too high for E&L data.
             Write (Output_Unit_HC_Alert,*) "Warning: Plasma temperature is higher than applicable using Ehrhardt and Langer data (&
                                             &2,000 eV):", Plasma_Temperature,"eV."
          End If
 
          ! Increase counter for number of reactions.
          i = i + 1
       End Do
    End If
 
  End Subroutine HC_Bounds_Check_Ehrhardt_Langer
 
End Module HC_NewSt
