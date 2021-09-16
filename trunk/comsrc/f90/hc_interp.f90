! -*-Mode:f90-*-
! HC Interpolation Routine.f90
! Hydrocarbon Reaction Program
!
! Adam McLean under Prof. Peter Stangeby
! Research Assistant: David Elder
! October, 1999
 
Module HC_Interp
 
Contains
 
  Subroutine HC_Interpolate (Reaction, Data_Type, Plasma_Temperature, Target_Energy, Interpolant)
 
    ! Gain access to data structures.
    Use HC_Init_Lib_Data
 
    ! Every good Fortran program has...
    Implicit None
 
    ! Define externally important variables.
    Integer, Intent (IN) :: Reaction
    Character (len=*), Intent (In) :: Data_Type
    Real, Intent (IN) :: Plasma_Temperature, Target_Energy
    Real, Intent (OUT) :: Interpolant
 
    ! Define internal veriables.
    Integer :: i, j, k
    Interpolant = 0
    If (Data_Type .eq. 'Sigma') Then
       ! Check transform table for target particle dependence and polynomial order.
       If (HC_State_Transform_Table (Reaction) % Sigma_TPD .eq. 1) Then
          ! Reaction is target particle independent.  Perform a single summation.
          Do j = 1, HC_State_Transform_Table (Reaction) % Sigma_Polynomial_Terms
             Interpolant = Interpolant + Sigma_Data (Reaction, j, 1) * (LOG (Plasma_Temperature))**(j-1)
          End Do
       Else If (HC_State_Transform_Table (Reaction) % Sigma_TPD .eq. 2) Then
          ! Reaction is target particle dependent.  Perform a double summation.
          Do j = 1, HC_State_Transform_Table (Reaction) % Sigma_Polynomial_Terms, 1
             Do k = 1, HC_State_Transform_Table (Reaction) % Sigma_Polynomial_Terms, 1
                Interpolant=Interpolant+Sigma_Data(Reaction,j,k)*((LOG(Target_Energy))**(j-1))*((LOG(Plasma_Temperature))**(k-1))
             End Do
          End Do
       elseif (HC_State_Transform_Table (Reaction) % Sigma_TPD .eq. 3) then 
          ! jdemod - added support for tabulated data and table interpolation
          !
          ! Data is tabulated sigma
          ! Perform table lookup - energy dependent reactions will have a p_energy array size greater than 1. 
          call table_interpolate(interpolant,Plasma_Temperature, Target_Energy,hc_state_transform_table(reaction)%t_index,&
                                &hc_state_transform_table(reaction)%e_index,hc_state_transform_table(reaction)%reaction_data)
       Else
          Write (Output_Unit_HC_Alert,*) 'Unknown target particle dependence value in HC_Interpolate: ', HC_State_Transform_Table (&
                                            &Reaction) % Sigma_TPD
          Write (Output_Unit_HC_Alert,*) "Program stopping."
          Stop
       End If
 
       ! jdemod - exponential is taken only for data types 1 and 2
       if (HC_State_Transform_Table (Reaction) % Sigma_TPD .eq. 1.or.HC_State_Transform_Table (Reaction) % Sigma_TPD .eq. 2) then 
          ! Take exponential of interpolant.
          Interpolant = EXP (Interpolant)
       endif


    Else If (Data_Type .eq. 'SigmaV') Then
       ! Interpolant will be a reaction rate.
       ! Check transform table for target particle dependence and polynomial order.
       If (HC_State_Transform_Table (Reaction) % SigmaV_TPD .eq. 1) Then
          ! Reaction is target particle independent.  Perform a single summation.
          Do j = 1, HC_State_Transform_Table (Reaction) % SigmaV_Polynomial_Terms, 1
             Interpolant = Interpolant + SigmaV_Data (Reaction, j, 1) * (LOG (Plasma_Temperature))**(j-1)
          End Do 
       Else If (HC_State_Transform_Table (Reaction) % SigmaV_TPD .eq. 2) Then
          ! Reaction is target particle dependent.  Perform a double summation.
          Do j = 1, HC_State_Transform_Table (Reaction) % SigmaV_Polynomial_Terms, 1
             Do k = 1, HC_State_Transform_Table (Reaction) % SigmaV_Polynomial_Terms, 1
                Interpolant=Interpolant+SigmaV_Data(Reaction,j,k)*((LOG(Target_Energy))**(j-1))*((LOG(Plasma_Temperature))**(k-1))
             End Do
          End Do
       elseif (HC_State_Transform_Table (Reaction) % SigmaV_TPD .eq. 3) then 
          ! jdemod - added support for tabulated data and table interpolation
          !
          ! Data is tabulated sigmav
          ! Perform table lookup - energy dependent reactions will have a p_energy array size greater than 1. 
          call table_interpolate(interpolant,Plasma_Temperature, Target_Energy,hc_state_transform_table(reaction)%t_index,&
                                &hc_state_transform_table(reaction)%e_index,hc_state_transform_table(reaction)%reaction_data)


       Else
          Write (Output_Unit_HC_Alert,*) 'Unknown target particle dependence value in HC_Interpolate: ', HC_State_Transform_Table (&
                                                                                  &Reaction) % SigmaV_TPD
          Write (Output_Unit_HC_Alert,*) "Program stopping."
          Stop
       End If

 
       ! jdemod - exponential is taken only for data types 1 and 2
       if (HC_State_Transform_Table (Reaction) % SigmaV_TPD .eq. 1.or.HC_State_Transform_Table (Reaction) % SigmaV_TPD .eq. 2) then
          ! Take exponential of interpolant.
          Interpolant = EXP (Interpolant)
       endif

    Else
       Write (Output_Unit_HC_Alert,*) 'Data type not defined in call to HC_Interp: ', Data_Type
       Write (Output_Unit_HC_Alert,*) "Program stopping."
       Stop
    End If
 

  End Subroutine HC_Interpolate



  subroutine table_interpolate(interpolant,plasma_temperature,target_energy,t_index,e_index,reaction_data)
    implicit none
    real :: interpolant,plasma_temperature,target_energy
    real,dimension(:) :: t_index,e_index
    real,dimension(:,:) :: reaction_data

    !
    ! Local variables
    !
    real :: interp_top, interp_bot
    integer :: t_num,e_num,t_top,e_top
    integer,external :: ipos
    logical :: t_low, t_high

    t_num = size(t_index)
    e_num = size(e_index)
    
    ! reaction data should have a size of (t_num,e_num)

    ! Need to find indices for interpolation

    !
    ! t_top and e_top are the table values needed for the interpolation
    !
    ! If e_top = 1 then no energy interpolation is required
    !

    t_top = ipos(plasma_temperature,t_index,t_num)

    ! Set logicals for temperature range - a good compiler should do something similar but this way it is explicit
    ! Not worth doing for energy since it is not used as frequently
    t_low  = plasma_temperature.lt.t_index(1)
    t_high = plasma_temperature.gt.t_index(t_num)


    ! Situation of only 1 energy bin - most common
    if (e_num.eq.1) then 

       if (t_low) then 
          interpolant = reaction_data(1,1)
       elseif (t_high) then
          interpolant = reaction_data(t_num,1)
       else
          interpolant = reaction_data(t_top-1,1) + &
                      & (plasma_temperature-t_index(t_top-1)) / (t_index(t_top) - t_index(t_top-1)) * &
                      & (reaction_data(t_top,1)-reaction_data(t_top-1,1))
       endif

    else
       ! Array has both temperature and energy dependence
       e_top = ipos(target_energy,e_index,e_num)
       !
       ! Need to double interpolate - first on energy then temperature
       !

       if (target_energy.lt.e_index(1)) then 
          if (t_low) then 
             interp_top = reaction_data(1,1)
             interp_bot = interp_top
          elseif (t_high) then 
             interp_top = reaction_data(t_num,1)
             interp_bot = interp_top
          else
             interp_top = reaction_data(t_top,1)
             interp_bot = reaction_data(t_top-1,1)
          endif
       elseif (target_energy.gt.e_index(e_num)) then 
          if (t_low) then 
             interp_top = reaction_data(1,e_num)
             interp_bot = interp_top
          elseif (t_high) then 
             interp_top = reaction_data(t_num,e_num)
             interp_bot = interp_top
          else
             interp_top = reaction_data(t_top,e_num)
             interp_bot = reaction_data(t_top-1,e_num)
          endif
       else
          if (t_low) then 
             interp_top = reaction_data(1,e_top-1) + &
                      & (target_energy-e_index(e_top-1)) / (e_index(e_top) - e_index(e_top-1)) * &
                      & (reaction_data(1,e_top)-reaction_data(1,e_top-1))
             interp_bot = interp_top
          elseif (t_high) then 
             interp_top = reaction_data(t_num,e_top-1) + &
                      & (target_energy-e_index(e_top-1)) / (e_index(e_top) - e_index(e_top-1)) * &
                      & (reaction_data(t_num,e_top)-reaction_data(t_num,e_top-1))
             interp_bot = interp_top
          else
             interp_top = reaction_data(t_top,e_top-1) + &
                      & (target_energy-e_index(e_top-1)) / (e_index(e_top) - e_index(e_top-1)) * &
                      & (reaction_data(t_top,e_top)-reaction_data(t_top,e_top-1))
             interp_bot = reaction_data(t_top-1,e_top-1) + &
                      & (target_energy-e_index(e_top-1)) / (e_index(e_top) - e_index(e_top-1)) * &
                      & (reaction_data(t_top-1,e_top)-reaction_data(t_top-1,e_top-1))
          endif
       endif

       ! Use temperature to interpolate the energy data
       

       if (t_low.or.t_high) then 

          interpolant = interp_top

       else

          interpolant = interp_bot + &
                      & (plasma_temperature-t_index(t_top-1)) / (t_index(t_top) - t_index(t_top-1)) * &
                      & (interp_top-interp_bot)

       endif


    endif

    !if (e_num.eq.1) then 
    !   write(6,'(a,4i6,2x,10g12.5)') 'Table interpolate e:', t_top,t_num,e_top,e_num,&
    !                                                   &t_index(t_top-1),plasma_temperature,t_index(t_top),target_energy,&
    !                                                   &reaction_data(t_top-1,1),&
    !                                                   &interpolant,reaction_data(t_top,1)
    !else
    !   write(6,'(a,4i6,2x,10g12.5)') 'Table interpolate p:', t_top,t_num,e_top,e_num,&
    !                                                   &t_index(t_top-1),plasma_temperature,t_index(t_top),&
    !                                                   &interp_bot,interpolant,interp_top,&
    !                                                   &e_index(e_top-1),target_energy,e_index(e_top)
    !endif
 
  end subroutine table_interpolate

End Module HC_Interp
