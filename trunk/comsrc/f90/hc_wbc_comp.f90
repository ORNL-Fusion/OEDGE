! -*-Mode:f90-*-
! HC_WBC_Comp.f90
! Provides all routines for DIVIMP-HC to record data and produce
! a table to compare HC and C transport and deposition between
! it and WBC (Jeff Brooks, ANL).
!
! Adam McLean under Prof. Peter Stangeby
! Research Assistant: David Elder
! January, 2003
 
Module HC_WBC_Comp
 
  ! Every good Fortran program has...
  Implicit None	
 
Contains
 
  Subroutine Record_WBC_Neut_Event (Mass,Velocity,Temperature,Mass_Fraction)
    ! Records particle mass and temperature data for purpose of comparison with WBC.
 
    ! Required modules.
    Use ComHC ! Input options.
    Use HC_Init_DIV_Data ! Time steps.
    Use HC_Init_DIV_Diag ! HC_Average data storage.
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Declare input/output variables.
    Real, Intent (In) :: Mass !
    Real, Intent (In) :: Velocity ! VEL
    Real, Intent (In) :: Temperature ! TEMI
    Real, Intent (In) :: Mass_Fraction ! SPUTY
 
    ! Check for applicability.
    If (hc_follow_option .ne. 0 .and. hc_wbc_comp_option .ne. 0) Then
       ! Add to statistics.
       ! Note:  Mass must be converted from AMU to kg, then E=1/2 mv^2 from J to eV.  T is in eV.
       HC_Tot_Temp_At_WBC_Target (0) =  HC_Tot_Temp_At_WBC_Target (0) + Mass_Fraction * (0.5 * Mass * 1.67E-27 * Velocity * &
                                      & Velocity / 1.602E-19 + 2.0 * Temperature)
       HC_Tot_At_WBC_Target (0) =  HC_Tot_At_WBC_Target (0) + Mass_Fraction
       ! write (0,*) "Record WBC neut:",velocity,temperature,mass_fraction	
       write (Output_Unit_HC_Alert,'(A20,I6,4E12.4)') "WBC_Neut_Event:","0",velocity /  Neutral_Time_Step,temperature,mass,&
                     & Mass_Fraction * (0.5 * Mass * 1.67E-27 * Velocity * Velocity / 1.602E-19 + 2.0 * Temperature)
    End If
  End Subroutine Record_WBC_Neut_Event
 
  Subroutine Record_WBC_Ion_Event (Charge_State,Mass,Velocity,Temperature,Mass_Fraction)
    ! Records particle mass and temperature data for purpose of comparison with WBC.
 
    ! Required modules.
    Use ComHC ! Input options.
    Use HC_Init_DIV_Data ! Time steps.
    Use HC_Init_DIV_Diag ! HC_Average data storage.		
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Declare input/output variables.
    Integer, Intent (In) :: Charge_State ! IZ
    Real, Intent (In) :: Mass ! VEL
    Real, Intent (In) :: Velocity ! VEL
    Real, Intent (In) :: Temperature ! TEMI
    Real, Intent (In) :: Mass_Fraction ! SPUTY
 
    ! Check for applicability.
    If (hc_follow_option .ne. 0 .and. hc_wbc_comp_option .ne. 0) Then
       ! Add to statistics.
       ! Note:  Mass must be converted from AMU to kg, then E=1/2 mv^2 from J to eV.  T is in eV.
       HC_Tot_Temp_At_WBC_Target (Charge_State) = 	HC_Tot_Temp_At_WBC_Target (Charge_State) + Mass_Fraction * (0.5 * Mass * 1.67E-&
&27 * Velocity * Velocity / 1.602E-19 + 2.0 * Temperature)
       HC_Tot_At_WBC_Target (Charge_State) =  HC_Tot_At_WBC_Target (Charge_State) + Mass_Fraction
       ! write (0,*) "Record WBC ion:",charge_state,velocity,temperature,mass_fraction
       write (Output_Unit_HC_Alert,'(A20,I6,4E12.4)') "WBC_ION_Event: ",charge_state,velocity /  Ion_Time_Step,temperature,mass,&
&Mass_Fraction * (0.5 * Mass * 1.67E-27 * Velocity * Velocity / 1.602E-19 + 2.0 * Temperature)
    End If
  End Subroutine Record_WBC_Ion_Event
 
  Integer Function Check_WBC_Neut_Position (Current_R,Current_Z,Mass,Velocity,Temperature,Mass_Fraction,Num_Pri_Launch,&
&Num_Sec_Launch,Number_Produced,Step_Count)
 
    ! Check if a neutral is inside the WBC geometry boundary and count if not.
    Use ComHC ! Input options.
    Use HC_Init_DIV_Diag ! HC_Average data storage.
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Declare input/output variables.
    Real, Intent (In) :: Current_R ! R
    Real, Intent (In) :: Current_Z ! Z
    Real, Intent (In) :: Mass ! CRMI
    Real, Intent (In) :: Velocity ! VEL
    Real, Intent (In) :: Temperature ! TEMI
    Real, Intent (In) :: Mass_Fraction ! SPUTY
    Integer, Intent (In) :: Num_Pri_Launch ! NPROD
    Integer, Intent (In) :: Num_Sec_Launch ! LPROD
    Integer, Intent (In) :: Number_Produced ! IPROD
    Real, Intent (In) :: Step_Count ! CIST
 
    Check_WBC_Neut_Position = 0
 
    If (hc_follow_option .ne. 0 .and. hc_wbc_comp_option .ne. 0) Then
 
       ! Write (0,'(a,2i,8e)') "Neut WBC Boundary check:",Num_Pri_Launch-Num_Sec_Launch+1,Number_Produced,Temperature,Velocity,SUM ( HC_Total_R_Prod_Positions) / &
       ! & (Num_Pri_Launch-Num_Sec_Launch+1),Current_R,Current_R - (Num_Pri_Launch-Num_Sec_Launch+1),SUM ( HC_Total_Z_Prod_Positions) / &
       ! & (Num_Pri_Launch-Num_Sec_Launch+1),Current_Z,Current_Z - SUM ( HC_Total_Z_Prod_Positions) / (Num_Pri_Launch-Num_Sec_Launch+1)
 
       If (Current_R .gt. (SUM ( HC_Total_R_Prod_Positions) / (Num_Pri_Launch-Num_Sec_Launch+1) + HC_WBC_Hori_Bound) .or. &
&Current_R .lt. (SUM ( HC_Total_R_Prod_Positions) / &
            & (Num_Pri_Launch-Num_Sec_Launch+1) - HC_WBC_Hori_Bound) .or. Current_Z .gt. (SUM ( HC_Total_Z_Prod_Positions) / (&
&Num_Pri_Launch-Num_Sec_Launch+1) + HC_WBC_Vert_Bound)) Then
 
          ! WBC comparison addition for neut out of bounds.
          HC_Tot_Temp_At_WBC_Boundary (0) =  HC_Tot_Temp_At_WBC_Boundary (0) + &
               & Mass_Fraction * (0.5 * Mass * 1.67E-27 * Velocity * Velocity / 1.602E-19 + 2.0 * Temperature)
          HC_Tot_At_WBC_Boundary (0) =  HC_Tot_At_WBC_Boundary (0) + Mass_Fraction
          Check_WBC_Neut_Position = 1
          write (Output_Unit_HC_Alert,'(A20,I6,4E12.4)') "WBC_Neut_Boundary:","0",velocity,temperature,mass,Mass_Fraction * (0.5 * &
&Mass * 1.67E-27 * Velocity * Velocity / 1.602E-19 + 2.0 * Temperature)
       End If
    End If
 
  End Function Check_WBC_Neut_Position
 
  Integer Function Check_WBC_Ion_Position (Current_Cell,Current_Ring,Current_S,Current_Cross,Charge_State,Mass,Velocity_In_S,&
&Temperature,Mass_Fraction,Num_Pri_Launch,Num_Sec_Launch,Number_Ionized,Impurity_Num)
 
    ! Check if an ion is inside the WBC geometry boundary and count if not.
    Use ComHC ! Input options.
    Use HC_Init_DIV_Diag ! HC_Average data storage.
    Use HC_Init_DIV_Data ! Gain access to cell/global/geometry/current data structures.
    Use HC_Get ! Contains gkbfs, gksmaxs.
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Declare input/output variables.
    Integer, Intent (In) :: Current_Cell ! IK
    Integer, Intent (In) :: Current_Ring ! IR
    Real, Intent (In) :: Current_S ! S
    Real, Intent (In) :: Current_Cross ! CROSS
    Integer, Intent (In) :: Charge_State ! IZ
    Real, Intent (In) :: Mass ! CRMI
    Real, Intent (In) :: Velocity_In_S ! VEL
    Real, Intent (In) :: Temperature ! TEMI
    Real, Intent (In) :: Mass_Fraction ! SPUTY
    Integer, Intent (In) :: Num_Pri_Launch ! NIMPS
    Integer, Intent (In) :: Num_Sec_Launch ! NIMPS2
    Integer, Intent (In) :: Number_Ionized ! NATIZ
    Integer, Intent (In) :: Impurity_Num ! IMP
 
    ! Declare local variables.
    Real :: R_Temp
    Real :: Z_Temp
 
    Check_WBC_Ion_Position = 0
 
    If (hc_follow_option .ne. 0 .and. hc_wbc_comp_option .ne. 0) Then
       ! Find approximate R,Z for current location.
       Call GETRZ (Current_Cell,Current_Ring,Current_S,Current_Cross,R_Temp,Z_Temp, RZ_Opt)
 
       ! Write (0,'(a,3i,11e)') "Ion WBC Boundary check:",Num_Pri_Launch+Num_Sec_Launch, Number_Ionized,Impurity_Num,Temperature,Velocity_In_S,R_Temp,R_Temp - SUM ( HC_Total_R_Prod_Positions) / (Num_Pri_Launch+Num_Sec_Launch),Z_Temp,Z_Temp - &
       ! & SUM ( HC_Total_Z_Prod_Positions)/ (Num_Pri_Launch+Num_Sec_Launch),gksmaxs(Current_Ring)-Current_S,Current_Cross,(gksmaxs(Current_Ring)-Current_S)*SIN(ASIN(1/gkbfs(Current_Cell,Current_Ring)))
 
       If (R_Temp .gt. (SUM ( HC_Total_R_Prod_Positions) / (Num_Pri_Launch+Num_Sec_Launch) + HC_WBC_Hori_Bound) &
            & .or. R_Temp .lt. (SUM ( HC_Total_R_Prod_Positions) / (Num_Pri_Launch+Num_Sec_Launch) - HC_WBC_Hori_Bound) &
            & .or. Z_Temp .gt. (SUM ( HC_Total_Z_Prod_Positions) / (Num_Pri_Launch+Num_Sec_Launch) + HC_WBC_Vert_Bound)) Then
 
          ! WBC comparison addition for neut out of bounds.
          HC_Tot_Temp_At_WBC_Boundary (Charge_State) =  HC_Tot_Temp_At_WBC_Boundary (Charge_State) + &
               & Mass_Fraction * (0.5 * Mass * 1.67E-27 * Velocity_In_S * Velocity_In_S / 1.602E-19 + 2.0 * Temperature)
          HC_Tot_At_WBC_Boundary (Charge_State) =  HC_Tot_At_WBC_Boundary (Charge_State) + Mass_Fraction
          Check_WBC_Ion_Position = 1
          Write (Output_Unit_HC_Alert,'(A20,I6,4E12.4)') "WBC_Ion_Boundary:",Charge_State,velocity_In_S,temperature,mass,&
&Mass_Fraction * (0.5 * Mass * 1.67E-27 * Velocity_In_S * Velocity_In_S / 1.602E-19 + 2.0 * Temperature)
       End If
    End If
 
  End Function Check_WBC_Ion_Position
 
  Subroutine Print_WBC_Comp (Charge_States)
    ! Prints the WBC comparison table as specified by Jeff Brooks.
 
    ! Required modules.
    Use ComHC ! Input options.
    Use HC_Init_Lib_Data ! Contains Get_HC_Charge function.
    Use HC_Init_DIV_Data ! Contains Get_HC_Charge function.
    Use HC_Init_DIV_Diag ! HC_Average data storage.
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Declare input/output variables.
    Integer, Intent (In) :: Charge_States
 
    ! Declare local variables.
    Integer :: Charge_State ! IZ
    Integer :: HC_Species
    Integer :: Num_HCs
    Integer :: Tot_Redep
    Integer :: Tot_Noredep
 
    ! Assign number of hydrocarbon species available.
    If (hc_higher_hcs_option .eq. 0) Then
       Num_HCs = 10
    Else
       Num_HCs = 58
    End If
 
    ! Re-open main HC output file to write WBC table at the end.
    Open (Unit = Output_Unit_HC_Data, File = "hc_output.txt", Status = 'Old', Position = 'Append')
 
    ! WBC comparison table title.
    Write (Output_Unit_HC_Data,*) "WBC Comparison Table"
    Write (Output_Unit_HC_Data,'(a7,4a10)') "Species","No. Redep","Avg E eV","# Not-Red","Avg E eV"
 
    ! Combine C neutral data from both HC following and regular DIVIMP following.
    HC_Species = 2
    Charge_State = Get_HC_Charge (HC_Species)
 
    Write (Output_Unit_HC_Data,'(a2,i1,4x,i10,f10.4,i10,f10.4)') "C+",Charge_State, &
         & INT ( HC_Tot_At_WBC_Target (Charge_State)) + &
         & SUM ( HC_IFate_Count (HC_Species,:,10)) + &
         & SUM ( HC_IFate_Count (HC_Species,:,12)) + &
         & SUM ( HC_IFate_Count (HC_Species,:,14)), &
 
         & ( HC_Tot_Temp_At_WBC_Target (Charge_State) + &
         & SUM ( HC_Tot_Temp_At_WBC_HC_Target (HC_Species,:))) / &
         & ( HC_Tot_At_WBC_Target (Charge_State) + &
         & SUM ( HC_IFate_Count (HC_Species,:,10)) + &
         & SUM ( HC_IFate_Count (HC_Species,:,12)) + &
         & SUM ( HC_IFate_Count (HC_Species,:,14)) +  Calc_Lo), &
 
         & INT ( HC_Tot_At_WBC_Boundary (Charge_State)) + &
         & SUM ( HC_IFate_Count (HC_Species,:,11)) + &
         & SUM ( HC_IFate_Count (HC_Species,:,15)), &
 
         & ( HC_Tot_Temp_At_WBC_Boundary (Charge_State) + &
         & SUM ( HC_Tot_Temp_At_WBC_HC_Boundary (HC_Species,:))) / &
         & (INT( HC_Tot_At_WBC_Boundary (Charge_State)) + &
         & SUM ( HC_IFate_Count (HC_Species,:,11)) + &
         & SUM ( HC_IFate_Count (HC_Species,:,15)) +  Calc_Lo)
 
    If (INT( HC_Tot_At_WBC_Boundary (Charge_State)) + SUM ( HC_IFate_Count (HC_Species,:,11)) + &
         & SUM ( HC_IFate_Count (HC_Species,:,15)) .gt. 0.0) Then
       write (Output_Unit_HC_Alert,*) "C0 nums", HC_Tot_At_WBC_Boundary (Charge_State),INT( HC_Tot_At_WBC_Boundary (Charge_State)),&
& HC_IFate_Count (HC_Species,1,11),&
            &  HC_IFate_Count (HC_Species,2,11), HC_IFate_Count (HC_Species,1,15), HC_IFate_Count (HC_Species,2,15)
       write (Output_Unit_HC_Alert,*) "C0 temps", HC_Tot_Temp_At_WBC_Boundary (Charge_State), HC_Tot_Temp_At_WBC_HC_Boundary (&
&HC_Species,1), HC_Tot_Temp_At_WBC_HC_Boundary (HC_Species,2)
    End If
 
    Do Charge_State = 1, Charge_States
       Write (Output_Unit_HC_Data,'(a2,i1,4x,i10,f10.4,i10,f10.4)') "C+",Charge_State, &
            & INT ( HC_Tot_At_WBC_Target (Charge_State)), HC_Tot_Temp_At_WBC_Target (Charge_State) / &
            & ( HC_Tot_At_WBC_Target (Charge_State) +  Calc_Lo), &
            & INT ( HC_Tot_At_WBC_Boundary (Charge_State)), HC_Tot_Temp_At_WBC_Boundary (Charge_State) / &
            & ( HC_Tot_At_WBC_Boundary (Charge_State) +  Calc_Lo)
    End Do
 
    !Write (0,*) "DATER:",SUM ( HC_IFate_Count (HC_Species,:,10)),SUM ( HC_Tot_Temp_At_WBC_HC_Target (HC_Species,:)),Num_HCs
 
    ! Note:  The C+ and C0 states (HC state 1 and 2) are reported above, and as soon as a HC reaches
    ! state 1 (C+), it is passed back to DIVIMP.  Therefore, reporting should start at 3.
    Do HC_Species = 3, Num_HCs
       If (Get_HC_Charge (HC_Species) .eq. 0) Then
          ! Redep: Sum of IFates 10,12,14.  No redep: Sum of IFates 11,15,19.  Except for C which is handled by DIVIMP.
          Write (Output_Unit_HC_Data,81) HC_State_Table (HC_Species) % State_Name, &
               & SUM ( HC_IFate_Count (HC_Species,:,10)) + &
               & SUM ( HC_IFate_Count (HC_Species,:,12)) + &
               & SUM ( HC_IFate_Count (HC_Species,:,14)), &
 
               & SUM ( HC_Tot_Temp_At_WBC_HC_Target (HC_Species,:)) / &
               & (SUM ( HC_IFate_Count (HC_Species,:,10)) + &
               & SUM ( HC_IFate_Count (HC_Species,:,12)) + &
               & SUM ( HC_IFate_Count (HC_Species,:,14)) +  Calc_Lo), &
 
               & SUM ( HC_IFate_Count (HC_Species,:,11)) + &
               & SUM ( HC_IFate_Count (HC_Species,:,15)), &
               & SUM ( HC_Tot_Temp_At_WBC_HC_Boundary (HC_Species,:)) / &
               & (SUM ( HC_IFate_Count (HC_Species,:,11)) + &
               & SUM ( HC_IFate_Count (HC_Species,:,15)) +  Calc_Lo)
 
          ! & (Tot_Reach_State_Per_Launch (HC_Species,1) + Tot_Reach_State_Per_Launch (HC_Species,2))
81        Format (a7,i10,f10.4,i10,f10.4)
       Else
          ! Redep: Sum of IFates 20,22,25,26.  No redep: Sum of IFates 21,27.  Except for C which is handled by DIVIMP.
          Write (Output_Unit_HC_Data,81) HC_State_Table (HC_Species) % State_Name, &
               & SUM ( HC_IFate_Count (HC_Species,:,20)) + &
               & SUM ( HC_IFate_Count (HC_Species,:,22)) + &
               & SUM ( HC_IFate_Count (HC_Species,:,25)) + &
               & SUM ( HC_IFate_Count (HC_Species,:,26)), &
 
               & SUM ( HC_Tot_Temp_At_WBC_HC_Target (HC_Species,:)) / &
               & (SUM ( HC_IFate_Count (HC_Species,:,20)) + &
               & SUM ( HC_IFate_Count (HC_Species,:,22)) + &
               & SUM ( HC_IFate_Count (HC_Species,:,25)) + &
               & SUM ( HC_IFate_Count (HC_Species,:,26)) +  Calc_Lo), &
 
               & SUM ( HC_IFate_Count (HC_Species,:,21)) + &
               & SUM ( HC_IFate_Count (HC_Species,:,27)), &
 
               & SUM ( HC_Tot_Temp_At_WBC_HC_Boundary (HC_Species,:)) / &
               & (SUM ( HC_IFate_Count (HC_Species,:,21)) + &
               & SUM ( HC_IFate_Count (HC_Species,:,27)) +  Calc_Lo)
       End If
    End Do
 
    Tot_Redep = SUM (INT ( HC_Tot_At_WBC_Target)) + &
         & SUM ( HC_IFate_Count (:,:,10)) + &
         & SUM ( HC_IFate_Count (:,:,12)) + &
         & SUM ( HC_IFate_Count (:,:,14)) + &
         & SUM ( HC_IFate_Count (:,:,20)) + &
         & SUM ( HC_IFate_Count (:,:,22)) + &
         & SUM ( HC_IFate_Count (:,:,25)) + &
         & SUM ( HC_IFate_Count (:,:,26))
 
    Tot_Noredep = SUM (INT ( HC_Tot_At_WBC_Boundary)) + &
         & SUM ( HC_IFate_Count (:,:,11)) + &
         & SUM ( HC_IFate_Count (:,:,15)) + &
         & SUM ( HC_IFate_Count (:,:,21)) + &
         & SUM ( HC_IFate_Count (:,:,27))
 
    Write (Output_Unit_HC_Data,'(a9,i8,2x,a10,i8,2x,a8,i8)') "Tot redep",Tot_Redep,"Tot no-red",Tot_Noredep,"Sum",Tot_Redep + &
&Tot_Noredep
 
    ! Close output file once again to finish.
    Close (Unit = Output_Unit_HC_Data)
 
  End Subroutine Print_WBC_Comp
 
End Module HC_WBC_Comp
