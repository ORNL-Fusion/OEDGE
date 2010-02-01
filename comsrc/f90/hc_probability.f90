module hc_probability

!
! jdemod - this routine combines all of the separate small modules which 
!          implement the HC state transition probability code into one place
!
!          hc_rnprb
!          hc_tlprb
!          hc_nsprb
!          hc_rnlup
!
!

contains


  Subroutine HC_Reaction_Probability (Reactions_Possible, Data_Type, Plasma_Density, &
       & Plasma_Temperature, Target_Energy, Hydrogen_Mass, Timestep, Reaction_Probabilities)
    use hc_init_lib_data
    use hc_interp
    Implicit None

    ! Define externally read variables
    Integer, Dimension (:), Intent (IN) :: Reactions_Possible
    Character (len=*), Intent (IN) :: Data_Type
    Real, Intent (IN) :: Plasma_Density, Plasma_Temperature, Target_Energy, Hydrogen_Mass, Timestep
    Real, Dimension (:), Intent (OUT) :: Reaction_Probabilities

    ! Define internal variables
    Integer :: i
    Real :: Interpolant

    i = 1
    Do While (Reactions_Possible (i) .ne. 0)

       ! Find the cross section at specific temperature and target energy using specified data
       Call HC_Interpolate (Reactions_Possible (i),Data_Type,Plasma_Temperature,Target_Energy,Interpolant)
       !write (0,'(a,i5,a,5(1x,g12.5)') "REACTING",Reactions_Possible (i),Data_Type,Plasma_Temperature,Plasma_density,Target_Energy,Interpolant,Hydrogen_Mass

       ! jdemod - proton interaction rates may need to be adjusted for different masses of plasma particles. 
       ! Correct for hydrogen mass for proton/D/T reactions.
       ! Change this to check the reaction type of reactions_possible(i)
       ! This may not be valid - may need to change the mass dependence effects

       ! If (Reactions_Possible (i) .gt. Num_Electron_Reactions) Then
       if (hc_state_transform_table(Reactions_Possible (i))%reaction_type.eq.'p') then
          Interpolant = Interpolant / SQRT (Hydrogen_Mass)
       End If

       ! Calculate probability for individual reaction
       If (Data_Type .eq. 'Sigma') Then
          ! Data will give a cross section
          Reaction_Probabilities (i) = 1 - EXP ( - (Plasma_Density * Interpolant) * Timestep)
          If (Reaction_Probabilities (i) .gt. 0.99 .and. Reaction_Probabilities (i) .le. 1.0) Then
             Write (Output_Unit_HC_Alert,*) 'Warning:  A reaction probability greater than 99% was found.'// &
                                          & '  Try a shorter time step:',Reaction_Probabilities (i)
          End If
       Else If (Data_Type .eq. 'SigmaV') Then
          ! Data will give a reaction rate
          Reaction_Probabilities (i) = 1 - EXP ( - (Plasma_Density * Interpolant) * Timestep)
          If (Reaction_Probabilities (i) .gt. 0.99 .and. Reaction_Probabilities (i) .le. 1.0) Then
             Write (Output_Unit_HC_Alert,*) 'Warning:  A reaction probability greater than 99% was found.'// &
                                          & '  Try a shorter time step:',Reaction_Probabilities (i)
          End If
       Else
          Write (Output_Unit_HC_Alert,*) 'Unknown data type specified in HC_Reaction_Probability: ', Data_Type
          Stop
       End If
       !write (0,'(a,1x,i6,1x,a,10(1x,g12.5))') "REACTING",Reactions_Possible(i),Data_Type,Plasma_Density, &
       !& Plasma_Temperature,Target_Energy,Interpolant,Reaction_Probabilities(i),Hydrogen_Mass,Timestep
       ! Increase counter for number of reactions.
       i = i + 1
    End Do

  End Subroutine HC_Reaction_Probability




  Subroutine HC_Reaction_Lookup (Current_State, Reactions_Possible,number_reactions)
    use hc_init_lib_data
    Implicit None

    integer, intent(out) :: number_reactions
    Integer, Intent (In) :: Current_State
    Integer, Dimension (:), Intent (OUT) :: Reactions_Possible

    ! Begin reaction table position counter at 1.
    ! Integer :: j
    Integer :: k

    number_reactions = hc_reaction_table(current_state)%number_reactions

    ! jdemod - replaced with reaction table lookup
    do k = 1,hc_reaction_table(current_state)%number_reactions
       reactions_possible(k) = hc_reaction_table(current_state)%reaction(k)
    end do       

  End Subroutine HC_Reaction_Lookup


  Subroutine HC_State_Table_Lookup (Reactions_Possible, States_Possible)
    use hc_init_lib_data
    ! Purpose:  To return an array of states which the given reactions may possibly
    ! lead to.  This includes all states from breakup of higher hydrocarbons.

    ! All good Fortran programs have...
    Implicit None

    ! Declare input/output variables.
    Integer, Dimension (:), Intent (IN) :: Reactions_Possible
    Integer, Dimension (:,:), Intent (OUT) :: States_Possible
    Integer :: i,j
    i = 1

    Do While (Reactions_Possible (i) .ne. 0)
       ! Note:  For higher hydrocarbons, we ignore multiple breakup HC's for now, and load the
       ! additional ones into the stack structure within the HC_Inside routines.
       j = 1
       Do While (j .le. SIZE (States_Possible,2)) ! Currently, this will go to State_Possible (50,3).
          If (HC_State_Transform_Table (Reactions_Possible (i)) % End_C_States (j) .ne. 0) Then
             States_Possible (i,j) = HC_State_Transform_Table (Reactions_Possible (i)) % End_C_States (j)
             ! Do next possible end state if applicable.
          End If
          j = j + 1
       End Do

       ! Do next reaction.
       i = i + 1
    End Do

  End Subroutine HC_State_Table_Lookup



  Subroutine HC_NS_Probability (Reaction_Probabilities, NS_Probabilities,number_reactions)

    Implicit None

    integer, intent(IN) :: number_reactions
    Real, Dimension (:), Intent (IN) :: Reaction_Probabilities
    Real, Dimension (:), Intent (OUT) :: NS_Probabilities

    !
    ! jdemod - changed to return cumulative reaction probability
    !
    ! Take probabilities for individual reactions including lack of reaction
    ! and remove lack of reaction probability (probability given a reaction).

    !Real :: Sum_Probability
    Integer :: i 
    !Integer :: j

    !Sum_Probability = 0.0
    NS_Probabilities = 0.0

    ! Find sum of reaction probabilities with probability of lack of reaction.
    ! jdemod - pass in number of reactions to improve code efficiency
    Do i = 1, number_reactions
    !Do i = 1, size (Reaction_Probabilities), 1
       !Sum_Probability = Sum_Probability + Reaction_Probabilities (i)
       if (i.eq.1) then 
          NS_Probabilities(i) = Reaction_probabilities(i)
       else
          NS_Probabilities(i) = NS_Probabilities(i-1) + Reaction_probabilities(i)
       endif
    End Do


    !if (ns_probabilities(number_reactions).le.0.0) then 
    !   write(0,'(a,i7,1x,g12.5)') 'ERROR in HC_NS_Probability',number_reactions,ns_probabilities(number_reactions)
    !endif

    ! Find reaction probability given a reaction.
    ! Do j = 1, i - 1, 1

    ! jdemod
    !
    ! Only normalize if the total probability of a change of state is non-zero
    ! It is possible for a particle in a particular state to find itself in a part of the plasma with plasma conditions 
    ! such that the change of state probability in any given time step is numerically zero
    !
    if (ns_probabilities(number_reactions).gt.0) then 
       Do i = 1, number_reactions
          ! jdemod - renormalize the cumulative array
          NS_Probabilities (i) = NS_Probabilities(i) / NS_Probabilities(number_reactions)
         !NS_Probabilities (j) = Reaction_Probabilites (j) / Sum_Probability
       End Do
    endif

  End Subroutine HC_NS_Probability



  Subroutine HC_Total_Probability (Reaction_Probabilities, Total_Probability,number_reactions)

    Implicit None

    integer, intent(IN) :: number_reactions
    Real, Dimension (:), Intent (IN) :: Reaction_Probabilities
    Real, Intent (OUT) :: Total_Probability

    Integer :: i
    Real :: Total_Not_Probability

    i = 1
    Total_Not_Probability = 1.0

    ! jdemod - switch to number_reactions
    Do i = 1, number_reactions, 1
    !Do i = 1, size(Reaction_Probabilities), 1
       Total_Not_Probability = Total_Not_Probability * (1.0 - Reaction_Probabilities (i))
    End Do

    Total_Probability = 1.0 - Total_Not_Probability

  End Subroutine HC_Total_Probability



end module hc_probability
