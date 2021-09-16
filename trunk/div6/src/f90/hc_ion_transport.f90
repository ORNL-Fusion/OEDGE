! -*-Mode:f90-*-
! HC_Ion_Transport.f90
! Recalculation of all Spitzer transport coefficients
! for particular hydrocarbon fragment mass.
!
! Adam McLean
 
Module Hc_Ion_Transport
 
  Implicit None
 
Contains
 
  Subroutine Ion_Move (Current_R,Current_Z,Current_S,Current_Cross,Current_Cell,Current_Ring, &
       & Cur_HC_Spec,H_Isotope_Composition,Sput_Weight,Current_Angle,Current_Velocity_In_S, &
       & Current_Temperature,Current_Theta,Last_R,Last_Z,Last_S,Last_Cross, &
       & Last_Cell,Last_Ring,Last_Angle,Last_Velocity_In_S,Random_Numbers_Used,nrand,iprod,&
       & eq_total_ion_time_steps,hc_v,Debug)
 
    Use ComHC ! Number_H_Species.
    Use HC_Init_DIV_Data ! Gain access to data structures.
    Use HC_Init_DIV_Diag ! HC diagnostics.
    Use HC_Get ! gknbs, glambda1 functions.
    Use HC_Init_Lib_Data ! Get HC mass and charge.
    Use HC_Utilities ! Get SPUTY.
    use hc_velocity_type ! HC velocity data definition
    use hc_kinetics_options ! kinetics options

    ! Every good Fortran routine has...
    Implicit None
 
    ! Declare call-line variables.
    Real, Intent (InOut) :: Current_R
    Real, Intent (InOut) :: Current_Z
    Real, Intent (InOut) :: Current_S
    Real, Intent (InOut) :: Current_Cross
    Integer, Intent (InOut) :: Current_Cell
    Integer, Intent (InOut) :: Current_Ring
    Integer, Intent (In) :: Cur_HC_Spec
    Integer, Dimension (Number_H_Species), Intent (In) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (In) :: H_Isotope_Composition
    Real, Intent (In) :: Sput_Weight
    Real, Intent (In) :: Current_Angle
    Real, Intent (InOut) :: Current_Velocity_In_S
    Real, Intent (InOut) :: Current_Temperature
    Real, Intent (InOut) :: Current_Theta
    Real, Intent (Out) :: Last_R
    Real, Intent (Out) :: Last_Z
    Real, Intent (Out) :: Last_S
    Real, Intent (Out) :: Last_Cross
    Integer, Intent (Out) :: Last_Cell
    Integer, Intent (Out) :: Last_Ring
    Real, Intent (Out) :: Last_Angle
    Real, Intent (Out) :: Last_Velocity_In_S
    Integer, Intent (InOut) :: Random_Numbers_Used ! KK
    Integer, Intent (InOut) :: Nrand ! Total count of random numbers - some functions do not take random numbers from ranv
    integer, Intent (In) :: iprod                  ! Index of current particle  
    real*8, intent(In) :: eq_total_ion_time_steps  ! Total time for the particle

    type(hc_velocity_type1) :: hc_v ! hc velocity data structure - including temperatures

    Logical, Intent (In) :: Debug

    ! jdemod - need a local smax variable since some of the code used here changes the value and passing in a function result just isn't a good idea
    real :: smax_local 
    
    Real :: Current_HC_Sputy ! SPUTY
    Double Precision :: DSputy ! DSPUTY
    Real :: FVEL ! Impurity ion fluid parallel velocity (m/s).
    Real :: FVH ! Background ion fluid parallel velocity (m/s).
    Real :: bg_drftvel ! Background drift velocity.
    Real :: imp_drftvel ! Impurity drift velocity.
    Real :: pol_drftv ! True poloidal drift velocity * ion timestep.
    Real :: FF ! Force of friction.
    Real :: FE ! Electric force.
    Real :: FIG ! Ion temperature gradient force.
    Real :: FEG ! Electron temperature gradient force.
    Real :: FVG ! Impurity pressure gradient force.
    Real :: QUANT ! Sum of all force velocities.
    Double Precision :: DVPara
    Double Precision :: DSPara
    Real :: CHIpara
    Real :: Lambda2
    Real :: XKpara
    Real :: XDparapara
    Real :: K11,K12,K13,D11,D12,D13
    Real :: Target
    Real :: Cells_To_Target
    Integer :: Velocity_Bin ! IN
    Integer :: Velocity_Cell ! IKV
 
    Real :: Adjust ! Not currently used but required to be declared.
    Real :: DCross(4)
    Real :: FTotal ! FTOTAL
    Real :: QTIM_Squared ! QTIM2
    Real :: CRMHC ! Local mass of hydrocarbon molecule.
    Real :: CRMI ! Mass of DIVIMP impurity.
    Real :: EFACTHC ! Electric field factor.
    Integer :: COption ! Returned from Decision
    !
    ! Local variables
    !
    real ds_dperpz,delta_s_dperpz,ds_pinch,ds_kpinchs
    external delta_s_dperpz,ds_kpinchs
 
 
    ! Initialize variables.
    Current_HC_Sputy = Calc_Sputy (Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition))
    DSputy = DBLE (Current_HC_Sputy)
    QUANT = 0.0
    QTIM_Squared =  Ion_Time_Step ** 2
    DCross = 0.0 ! Array length 4.
    CRMHC = Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition)
    CRMI =  Impurity_Ion_Mass
    EFACTHC =  Ion_Time_Step *  Ion_Time_Step *  ECH /  AMU / CRMHC
 
    ! Save previous values.
    Last_S = Current_S
    Last_Cross = Current_Cross
    Last_Cell = Current_Cell
    Last_Ring = Current_Ring
    Last_Angle = Current_Angle
    Last_Velocity_In_S = Current_Velocity_In_S
    Last_R = Current_R
    Last_Z = Current_Z
 
    ! Modify velocity to be distance per timestep.
    Current_Velocity_In_S = Current_Velocity_In_S *  Ion_Time_Step
 
    ! Assign quantities used in force calculations.
    ! Assume magnetic field is aligned parallel to z axis.
    FVEL = Current_Velocity_In_S
    !write (Output_Unit_Scratch,*) "ionhc",Current_S,Current_Cross

    !
    ! Allow for drifts to be applied inside the confined plasma
    ! (if an HC fragment actually manages to get there)
    !
    !If (Current_Ring .lt.  Inner_SOL_Ring) Then
    !   ! Ion in MAIN plasma.  No drift is overlaped on velocity.
    !   FVH = gkvhs (Current_Cell, Current_Ring)
    !Else ! Ion in SOL or TRAPPED plasma.	
    !
       ! Calculate poloidal drift V*T.
       !
       ! jdemod - call a subroutine to get the drift velocity for the ring
       !
       ! pol_drftv =  Background_Drift_Velocity *  Ion_Time_Step
       !
       ! SOL_Drift_Start and SOL_Drift_End are loaded from the data calculated
       ! in the main DIV module - maintaining local values is no longer needed
       !
       ! Note: Drifts may now be applied in the confined plasma (possibly to simulate
       !       toroidal rotation? Code modified here to support this)
       !
       ! pol_drftv is initialized with a zero velocity for all rings so that this 
       ! option will not affect things if the poloidal drift is not on. 
       !
       call get_drftv(pol_drftv,SOL_Drift_Start,SOL_Drift_End,current_ring)
       !
 
       ! Assign poloidal drift if applicable.
       If ( Poloidal_Drift_Opt .eq. 0) Then
          imp_drftvel = 0.0
          bg_drftvel = 0.0				
       ElseIf ( (Poloidal_Drift_Opt .eq. 1 .or.Poloidal_Drift_Opt.eq.3) .and. &
            & (Current_S .ge.  SOL_Drift_Start .and. Current_S .le.  SOL_Drift_End)) Then
          imp_drftvel = pol_drftv
          bg_drftvel = 0.0
       ElseIf ( Poloidal_Drift_Opt .eq. 2 .and. (Current_S .ge.  SOL_Drift_Start .and. Current_S .le.  SOL_Drift_End)) Then
          bg_drftvel = pol_drftv
          imp_drftvel = 0.0
       End If
 
       ! Find flow velocity with drift included.
       FVH = gkvhs (Current_Cell, Current_Ring) + bg_drftvel
    ! End If
    !write (Output_Unit_Scratch,*) "ionhc",Current_S,Current_Cross,FVEL, Inner_SOL_Ring,FVH,gkvhs (Current_Cell, Current_Ring),bg_drftvel
 
    ! Select Reiser or normal force calculations.
    If ( Reiser_Opt .eq. 1 .or.  Reiser_Opt .eq. 2) Then
       ! Call decision to examine gradients in the cell and determine
       ! if it is appropriate to use the Reiser formulation for the forces.
 
       Call Decision (COption,Current_Cell,Current_Ring)
 
       If (COption .ne. 0) Then
          ! If option 2 - update coeffiecients each time step
          If ( Reiser_Opt .eq. 2) then
             Call update_reiser_coeff (Current_Cell,Current_Ring,Get_HC_Charge (Cur_HC_Spec),Current_S,lambda2,FVH)
          Else
             Lambda2 = glambda1 (Get_HC_Charge (Cur_HC_Spec)) * gknbs (Current_Cell,Current_Ring)
          End If
 
          ! Calculate Reiser forces
 
          CHIpara = galphai (Current_Cell,Current_Ring) * (FVEL - FVH)
 
          Call COULOMB_COLL (XKpara,XDparapara,CHIpara,Lambda2,Current_Cell,Current_Ring,K11,K12,K13,D11,D12,D13)
 
          FF = K11 * gqtim2()
          FIG = K12 * gqtim2()
          FVG = K13 * gqtim2()
 
          Random_Numbers_Used = Random_Numbers_Used + 1
 
          ! Reset Vpara for Reiser fomulation.
          Local_VPara  = SQRT (XDparapara *  Ion_Time_Step) *  Local_RGauss *  Ion_Time_Step
          Local_SPara  = 0.0
 
       Else
          ! Default force behaviour - as standard.
          FF = gkfssmod (Current_Cell, Current_Ring)  *  Local_HC_Tau_Stopping_Inv * (FVH - FVEL)
          FIG =  Local_Betas * gkfigs (Current_Cell,Current_Ring) * CRMI / CRMHC
          FVG = 0.0
 
       End If
 
    Else
       ! Default force behaviour - as standard.
       FF = gkfssmod (Current_Cell,Current_Ring) *  Local_HC_Tau_Stopping_Inv * (FVH - FVEL)
       FIG =  Local_Betas * gkfigs (Current_Cell,Current_Ring) * CRMI / CRMHC
       FVG = 0.0
 
    End If
 
    FEG =  Local_Alphs * gkfegs (Current_Cell, Current_Ring) * CRMI / CRMHC
 
    ! Modify electric field strength if within applicable range of target.
    ! Check for applicability of sheath electric field modification.
    If (hc_presheath_efield .eq. 1) Then
 
       ! Do only for SOL rings (IR .ge. IRSEP).
       If (Current_Ring .ge.  Inner_SOL_Ring) Then
 
          ! Record cells to target.
          Cells_To_Target = MIN (Current_Cell, gnks (Current_Ring) - Current_Cell)
 
          ! Check if within magnetic pre-sheath.
          If (Cells_To_Target .le. hc_efield_cells) Then
             ! Add Brooks modification to local electric field strength.
             Local_Electric_Field = Sheath_E_Field (Current_Cell,Current_Ring,Current_S, &
                  &  Plate_Electron_Temperature,gktebs (Current_Cell,Current_Ring),gktibs (Current_Cell,Current_Ring), &
                  & gknbs (Current_Cell,Current_Ring),  Back_Plasma_Ion_Mass, Back_Plasma_Charge, &
                  & hc_efield_drop_fraction,gbts (Current_Cell,Current_Ring),gkes (current_Cell,Current_Ring))
             ! jdemod - *** BUG *** - since the sheath_e_field code can return gkes - it is incorrect to always scale the
             !          result by EFACTHC - check and scale appropriately
             if ( Local_Electric_Field.ne.gkes(current_Cell,Current_Ring)) then
                Local_Electric_Field =  Local_Electric_Field * EFACTHC
             else
                Local_Electric_Field =  Local_Electric_Field * CRMI/CRMHC
             endif
             !write (0,*) "BROOKS", Local_Electric_Field,EFACTHC
          Else
             Local_Electric_Field = gkes (Current_Cell,Current_Ring) * CRMI / CRMHC
             !write (0,*) "NOT BROOKS",gkes (Current_Cell,Current_Ring),CRMI,CRMHC
          End If
       End If
    Else
       ! Use standard DIVIMP electric field.
       Local_Electric_Field = gkes (Current_Cell,Current_Ring) * CRMI / CRMHC
       !write (0,*) "NOT BROOKS2",gkes (Current_Cell,Current_Ring),CRMI,CRMHC
    End If
 
    FE = Get_HC_Charge (Cur_HC_Spec) *  Local_Electric_Field
    !if (gkes(current_Cell,Current_Ring) .ne. 0.0) Then
    !write (0,*) "NOT ZZERO",gkes(current_Cell,Current_Ring),Current_Cell,Current_Ring
    !endif
 
    !write (0,*) "CHARGE,LEF",Get_HC_Charge (Cur_HC_Spec), Local_Electric_Field, &
    !& gkes (Current_Cell,Current_Ring),CRMI,CRMHC,hc_presheath_efield
 
    ! Calculate modifications to forces if any.
    If ( TEB_Coeff_Opt .eq. 3 .and. &
         & Current_S .gt.  TGrad_Zero_Dist * gksmaxs (Current_Ring) .and. &
         & Current_S .lt. gksmaxs (Current_Ring) * (1.0 -  TGrad_Zero_Dist)) Then
       FEG = 0.0
    End If
 
    If ( TIB_Coeff_Opt .eq. 3 .and. &
         & Current_S .gt.  TGrad_Zero_Dist * gksmaxs (Current_Ring) .and. &
         & Current_S .lt. gksmaxs (Current_Ring) * (1.0 -  TGrad_Zero_Dist)) Then
       FIG = 0.0
    End If
 
    ! Ion movement parallel to magnetic field.
    !FF = 0.0
    !FE = 0.0
    !FEG = 0.0
    !FIG = 0.0
    !FVG = 0.0
    QUANT = FF + FE + FEG + FIG + FVG
    Random_Numbers_Used = Random_Numbers_Used + 1
 
    DSPara = DBLE (SIGN ( Local_SPara, granv (Random_Numbers_Used) - 0.5))
    DVPara = DBLE (SIGN ( Local_VPara, granv (Random_Numbers_Used) - 0.5))
    !DSPara = 0.0
    !DVPara = 0.0
 
    !write (Output_Unit_Scratch,*) "QUANT HC",QUANT,FF,FE,FEG,FIG,FVG,FVH,FVEL, Local_VPara,dvpara,DSPara,REAL(Current_Ring), &
    !& REAL( Inner_SOL_Ring),Current_Velocity_In_S
 
    !write (0,'(A,20E12.4)') "QUANT HC",QUANT,FF,FE,FEG,FIG,FVG,FVH,FVEL, Local_VPara,dvpara, &
    !& DSPara,REAL(Current_Ring),REAL( Inner_SOL_Ring),Current_Velocity_In_S,Current_S
 
    !write (0,*) "More forces", Local_Electric_Field,gkes (Current_Cell,Current_Ring), &
    !& gbts (Current_Cell,Current_Ring),DSPara,DVPara, &
    !&  Local_SPara, Local_VPara,Current_S, &
    !& Current_Velocity_In_S,REAL(Current_Ring),REAL( Inner_SOL_Ring)
 
    !
    ! jdemod
    !
    ! Added a parallel transport effect resulting from cross-field diffusion along the axis
    ! perpendicular to both the poloidal and total field line directions.
    !
    ds_dperpz=delta_s_dperpz(current_cell,current_ring,Random_Numbers_Used)

    !
    ! jdemod - parallel motion due to radial drift - only non-zero when the option is active. 
    !
    ds_pinch = ds_kpinchs(current_cell,current_ring)
 
    ! Find new S position.
    !If (Current_Ring .lt.  Inner_SOL_Ring) Then
    ! Ion in main plasma.  No drift will be added.
    !   Current_S = Current_S + Current_Velocity_In_S + 0.5 * (QUANT + DVPara) + DSPara + ds_dperpz
    !Else ! Ion in SOL or TRAP.  Add drift velocity calculated above.

    !
    ! jdemod - the get_drftv rotuine above returns the appropriate value for the drift velocity for 
    !          any ring so splitting core/sol transport to avoid the drift is no longer necessary
    !
       Current_S = Current_S + Current_Velocity_In_S + 0.5 * (QUANT + DVPara) + DSPara + imp_drftvel + ds_dperpz + ds_pinch

    !End If
 
    !
    ! jdemod
    !
 
    !write (0,*) "forces",QUANT,DVPara,ff,fe,feg,fig,fvg,FVH,FVEL,gkfssmod (Current_Cell,Current_Ring), &
    !&  Local_HC_Tau_Stopping_Inv, Local_Betas, &
    !& gkfigs (Current_Cell,Current_Ring),CRMI,CRMHC
 
    Current_Velocity_In_S = Current_Velocity_In_S + QUANT + DVPara
    !write (Output_Unit_Scratch,*) "vels HC",Current_Velocity_In_S,QUANT,DVPara,Current_S,DSPara,imp_drftvel
 
    !
    ! Allow for diffusion of the perpendicular velocity if the corresponding vperp kinetics option is selected
    !
    if (hc_kinetics_opt.eq.1) then 
       if (hc_vperp_opt.eq.3) then 
          
          ! jdemod
          ! increment random numbers and change the perpendicular velocity by the same magnitude as the 
          ! parallel velocity. Note that the value of vperp is in actual m/s and not distance travelled in one
          ! time step - thus the change in vperp is modified by the ion_time_step which is included in local_vpara
          ! in earlier code - this option does not support spatial diffusion - i.e. changes in spara
          !
          random_numbers_used = random_numbers_used + 1
          hc_v%vperp = hc_v%vperp + sign(local_vpara/ion_time_step,granv (Random_Numbers_Used) - 0.5)

          if (debug_kinetics) &
               & write(6,'(a,10g12.5)') 'HC_VPERP_OPT_3:UPDATE VPERP:',hc_v%vperp,local_vpara/ion_time_step

       endif
    endif



    ! For collision option 7 - reverse sign of VEL under certain circumstances.
    If ( Collision_Opt .eq. 7 .and. DSPara .ne. 0.0) Then
       If (Current_Velocity_In_S .gt. 0.0 .and. DSPara .lt. 0.0) Then
          Current_Velocity_In_S = -Current_Velocity_In_S
       ElseIf (Current_Velocity_In_S .lt. 0.0 .and. DSPara .gt. 0.0) Then
          Current_Velocity_In_S = -Current_Velocity_In_S
       End If
    End If
 
    ! Advance ion specific counters.
    ! Note, Inner_SOL_Ring is equal to IRSEP in DIVIMP.
    If (Current_Ring .lt.  Inner_SOL_Ring) Then
       ! Ion in main plasma.
 
       ! Record average parallel steps for debugging purposes.
       If (DVPara .ne. 0.0) Then
          If (DVPara .lt. 0.0) Then
             HC_DVParaStep (5, Cur_HC_Spec) = &
                  &  HC_DVParaStep (5, Cur_HC_Spec) + DVPara
             HC_VParaStep (5, Cur_HC_Spec) = &
                  &  HC_VParaStep (5, Cur_HC_Spec) + ABS (Current_Velocity_In_S)
             HC_DVParaCnt (5, Cur_HC_Spec) = &
                  &  HC_DVParaCnt (5, Cur_HC_Spec) + 1.0
             HC_DVMaxV (5, Cur_HC_Spec) =  &
                  & MAX ( HC_DVMaxV (5, Cur_HC_Spec),DBLE (Current_Velocity_In_S))
             HC_DVMinV (5, Cur_HC_Spec) = &
                  & MIN ( HC_DVMinV (5, Cur_HC_Spec),DBLE (Current_Velocity_In_S))
          ElseIf (DVPara .gt. 0.0) Then
             HC_DVParaStep (6, Cur_HC_Spec) = &
                  &  HC_DVParaStep (6, Cur_HC_Spec) + DVPara
             HC_VParaStep (6, Cur_HC_Spec) = &
                  &  HC_VParaStep (6, Cur_HC_Spec) + ABS (Current_Velocity_In_S)
             HC_DVParaCnt (6, Cur_HC_Spec) = &
                  &  HC_DVParaCnt (6, Cur_HC_Spec) + 1.0
             HC_DVMaxV (6, Cur_HC_Spec) = &
                  & MAX ( HC_DVMaxV (6, Cur_HC_Spec),DBLE (Current_Velocity_In_S))
             HC_DVMinV (6, Cur_HC_Spec) = &
                  & MIN ( HC_DVMinV (6, Cur_HC_Spec),DBLE (Current_Velocity_In_S))
          End If
       End If
 
       ! Record some quantities in the core
       HC_CoreOuts (Cur_HC_Spec, 1) = &
            &  HC_CoreOuts (Cur_HC_Spec, 1) + DSputy
       HC_CoreOuts (Cur_HC_Spec, 2) = &
            &  HC_CoreOuts (Cur_HC_Spec, 2) + DSputy * Current_Temperature /  Local_HC_Change_State_Coll_Prob
       HC_CoreOuts (Cur_HC_Spec, 3) = &
            &  HC_CoreOuts (Cur_HC_Spec, 3) + DSputy /  Local_HC_Tau_Stopping_Inv
       If (Current_S .le. 0.5 * gksmaxs (Current_Ring)) Then
          HC_CoreOuts (Cur_HC_Spec, 4) = &
               &  HC_CoreOuts (Cur_HC_Spec, 4) + DSputy * FF
          HC_CoreOuts (Cur_HC_Spec, 5) = &
               &  HC_CoreOuts (Cur_HC_Spec, 5) + DSputy * FE
          HC_CoreOuts (Cur_HC_Spec, 6) = &
               &  HC_CoreOuts (Cur_HC_Spec, 6) + DSputy * FEG
          HC_CoreOuts (Cur_HC_Spec, 7) = &
               &  HC_CoreOuts (Cur_HC_Spec, 7) + DSputy * FIG
          HC_CoreOuts (Cur_HC_Spec, 8) = &
               &  HC_CoreOuts (Cur_HC_Spec, 8) + DSputy * FVEL
          HC_CoreOuts (Cur_HC_Spec, 9) = &
               &  HC_CoreOuts (Cur_HC_Spec, 9) + DSputy * FVH
       Else
          HC_CoreOuts (Cur_HC_Spec, 4) = &
               &  HC_CoreOuts (Cur_HC_Spec, 4) - DSputy * FF
          HC_CoreOuts (Cur_HC_Spec, 5) = &
               &  HC_CoreOuts (Cur_HC_Spec, 5) - DSputy * FE
          HC_CoreOuts (Cur_HC_Spec, 6) = &
               &  HC_CoreOuts (Cur_HC_Spec, 6) - DSputy * FEG
          HC_CoreOuts (Cur_HC_Spec, 7) = &
               &  HC_CoreOuts (Cur_HC_Spec, 7) - DSputy * FIG
          HC_CoreOuts (Cur_HC_Spec, 8) = &
               &  HC_CoreOuts (Cur_HC_Spec, 8) - DSputy * FVEL
          HC_CoreOuts (Cur_HC_Spec, 9) = &
               &  HC_CoreOuts (Cur_HC_Spec, 9) - DSputy * FVH
       EndIf
 
       ! Looping around main plasma contours.
       If (Current_S .lt. 0.0) Then
          Do
             Current_S = Current_S + gksmaxs (Current_Ring)
             If (Current_S .ge. 0.0) Then
                Exit
             End If
             Write (Output_Unit_Scratch,*) "Still here..adding.  Stopping"
             Stop
          End Do
       ElseIf (Current_S .gt. gksmaxs (Current_Ring)) Then
          Do
             Current_S = Current_S - gksmaxs (Current_Ring)
             If (Current_S .le. gksmaxs (Current_Ring)) Then
                Exit
             End If
             Write (Output_Unit_Scratch,*) "Still here...subtracting.  Stopping."
             Stop
          End Do
       End If
 
       ! Update the Cross-field transport term and make
       ! all other related position adjustments.
       If (debug) Then
          Write (Output_Unit_Scratch,'(a,2i4,1p,4g12.5,l4)') 'UP_CROSS:1A:',Current_Cell,Current_Ring,Current_S,Current_Theta,&
                                                             &Current_Cross,Debug
       End If
       !write (Output_Unit_Scratch,*) "hc updating cross1"		
 
       smax_local = gksmaxs(current_ring)

       Call update_cross (Current_Cell,Current_Ring,Last_Cell,Last_Ring,Random_Numbers_Used,Current_S,Current_Theta, &
            & Current_Cross,Adjust,DCross, CKK_Minimum,smax_local, Local_K,nrand,iprod,eq_total_ion_time_steps,Debug)
 
       If (Debug) Then
          Write (Output_Unit_Scratch,'(a,2i4,1p,4g12.5,l4)') 'UP_CROSS:1B:',Current_Cell,Current_Ring,Current_S,Current_Theta,&
                                                           &Current_Cross,Debug
       End If
 
    Else ! Ion in SOL or TRAP.
 
       ! Record average parallel steps for debugging purposes.
       If (DVPara .ne. 0.0) Then
          If (Current_S .lt. gksmaxs (Current_Ring) / 2.0) Then
             If (DVPara .lt. 0.0) Then
                HC_DVParaStep (1, Cur_HC_Spec) = &
                     &  HC_DVParaStep (1, Cur_HC_Spec) + DVPara
                HC_VParaStep (1, Cur_HC_Spec) = &
                     &  HC_VParaStep (1, Cur_HC_Spec) + ABS (Current_Velocity_In_S)
                HC_DVParaCnt (1, Cur_HC_Spec) = &
                     &  HC_DVParaCnt (1, Cur_HC_Spec) + 1.0
                HC_DVMaxV (1, Cur_HC_Spec) = &
                     & MAX ( HC_DVMaxV (1, Cur_HC_Spec), DBLE (Current_Velocity_In_S))
                HC_DVMinV (1, Cur_HC_Spec) = &
                     & MIN ( HC_DVMinV (1, Cur_HC_Spec), DBLE (Current_Velocity_In_S))
             ElseIf (DVPara .gt. 0.0) Then
                HC_DVParaStep (2, Cur_HC_Spec) = &
                     &  HC_DVParaStep (2, Cur_HC_Spec) + DVPara
                HC_VParaStep (2, Cur_HC_Spec) = &
                     &  HC_VParaStep (2, Cur_HC_Spec) + ABS (Current_Velocity_In_S)
                HC_DVParaCnt (2, Cur_HC_Spec) = &
                     &  HC_DVParaCnt (2, Cur_HC_Spec) + 1.0
                HC_DVMaxV (2, Cur_HC_Spec) = &
                     & MAX ( HC_DVMaxV (2, Cur_HC_Spec), DBLE (Current_Velocity_In_S))
                HC_DVMinV (2, Cur_HC_Spec) = &
                     & MIN ( HC_DVMinV (2, Cur_HC_Spec), DBLE (Current_Velocity_In_S))
             End If
          ElseIf (Current_S .gt. gksmaxs (Current_Ring) / 2.0) Then
             If (DVPara .lt. 0.0) Then
                HC_DVParaStep (3, Cur_HC_Spec) = &
                     &  HC_DVParaStep (3, Cur_HC_Spec) + DVPara
                HC_VParaStep (3, Cur_HC_Spec) = &
                     &  HC_VParaStep (3, Cur_HC_Spec) + ABS (Current_Velocity_In_S)
                HC_DVParaCnt (3, Cur_HC_Spec) = &
                     &  HC_DVParaCnt (3, Cur_HC_Spec) + 1.0
                HC_DVMaxV (3, Cur_HC_Spec) = &
                     & MAX ( HC_DVMaxV (3, Cur_HC_Spec), DBLE (Current_Velocity_In_S))
                HC_DVMinV (3, Cur_HC_Spec) = &
                     & MIN ( HC_DVMinV (3, Cur_HC_Spec), DBLE (Current_Velocity_In_S))
             ElseIf (DVPara .gt. 0.0) Then
                HC_DVParaStep (4, Cur_HC_Spec) = &
                     &  HC_DVParaStep (4, Cur_HC_Spec) + DVPara
                HC_VParaStep (4, Cur_HC_Spec) = &
                     &  HC_VParaStep (4, Cur_HC_Spec) + ABS (Current_Velocity_In_S)
                HC_DVParaCnt (4, Cur_HC_Spec) = &
                     &  HC_DVParaCnt (4, Cur_HC_Spec) + 1.0
                HC_DVMaxV (4, Cur_HC_Spec) = &
                     & MAX ( HC_DVMaxV (4, Cur_HC_Spec), DBLE (Current_Velocity_In_S))
                HC_DVMinV (4, Cur_HC_Spec) = &
                     & MIN ( HC_DVMinV (4, Cur_HC_Spec), DBLE (Current_Velocity_In_S))
             End If
          End If
       End If
 
       ! Spara - spatial diffusion.
 
       If (DSPara .ne. 0.0) Then
          If (Current_S .lt. gksmaxs (Current_Ring) / 2.0) Then
             If (DSPara.lt.0.0) Then
                HC_DSParaStep (1) =  HC_DSParaStep (1) + DSPara
                HC_DSParaCnt (1) =  HC_DSParaCnt (1) + 1.0
             ElseIf (DSPara.gt.0.0) Then
                HC_DSParaStep (2) =  HC_DSParaStep (2) + DSPara
                HC_DSParaCnt (2) =  HC_DSParaCnt (2) + 1.0
             End If
          ElseIf (Current_S .gt. gksmaxs (Current_Ring) / 2.0) Then
             If (DSPara .lt. 0.0) Then
                HC_DSParaStep (3) =  HC_DSParaStep (3) + DSPara
                HC_DSParaCnt (3) =  HC_DSParaCnt (3) + 1.0
             ElseIf (DSPara .gt. 0.0) Then
                HC_DSParaStep (4) =  HC_DSParaStep (4) + DSPara
                HC_DSParaCnt (4) =  HC_DSParaCnt (4) + 1.0
             End If
          End If
       End If
 
       ! Check if on grid.
       If (Current_S .le. 0.0 .or. Current_S .ge. gksmaxs (Current_Ring)) Then
          ! Off grid, do not update cross.
          ! This will lead to the hc_outside_ion routine.
       Else
          ! On grid, update cross.
 
          ! Print debug info if desired.
          If (Debug) Then
             Write (Output_Unit_Scratch,'(a,2i4,1p,4g12.5,L4)') 'Up_Cross:1A',Current_Cell,Current_Ring,Current_S,Current_Theta,&
                                                                &Current_Cross,Debug
          End If


 
          ! Update cross-field transport term.
          !write (0,*) "hc updating cross2",Current_Cell,Current_Ring,Current_S,Current_Theta,Current_Cross,Adjust,DCross
          !write (Output_Unit_Scratch,*) "hc updating cross2"

          smax_local = gksmaxs(current_ring)
          Call update_cross( 			& ! Call to update_cross found in cfield.d6a
               & Current_Cell,				& ! ik, input and output - ik index of particle location before and after update.
               & Current_Ring, 			& ! ir, input and output - ir index of particle location before and after update.
               & Last_Cell, 				& ! ikold, output - equals original particle ik location.
               & Last_Ring, 				& ! irold, output - equals original particle ir location.
               & Random_Numbers_Used, 			& ! kk, input - index into array of random numbers, ranv.
               & Current_S, 				& ! s, input - current s coordinate of the particle.
               & Current_Theta,		 	& ! theta, output - orthogonal coordinate along field line for given s value.
               & Current_Cross,		 	& ! cross, input - current position perpendicular to cell center line.
               & Adjust, 				& ! adjust, output calculated in adjust_cross - correction to cross-field coordinate when particle moves to another cell when expansion or compression may have occured.
               & DCross, 				& ! dcross, output - records statistics of number of cross-field steps in various grid regions.
               &  CKK_Minimum, 	& ! ckkmin, input and output - minimum k coordinate recorded by particle.
               & smax_local, 		& ! smax, input - s location of upper target (other target at s=0.0).
               &  Local_K, 		& ! k, input - k<1.0 in core, k=1.0 on separatrix, k>1.0 outside.
               & nrand,                  & ! total count of random numbers
               & iprod,                  & ! index of current particle
               & eq_total_ion_time_steps,& ! total ion equivalent time steps for particle
               & Debug) 				  ! debug, input - debugging information printed or not.
          
!
! jdemod - at some point in time need to update the transport code for HCs - problems exist at the moment in terms of code interdependencies
!
!      smax_local = gksmaxs(current_ring)
!      call do_crossfield_step(current_cell,current_ring,last_cell,last_ring,random_numbers_used,current_s,current_theta,current_cross,&
!                          &   last_theta,last_cross,&
!                          &   adjust,dcross,ckk_minimum,smax_local,local_k,debug,&
!                          &   seed,nrand,neutim,cist,imp,debug_all,&
!                          &   ifate)
!


          If (Debug) Then
             ! jdemod - format descriptor was 14 not l4 (one instead of "ELL")
             write (Output_Unit_Scratch,'(a,2i4,1p,4g12.5,l4)') 'Up_Cross:1B',Current_Cell,Current_Ring,Last_S,Current_Theta,&
                                                                &Current_Cross,Debug
          End If
       End If
 
    End If
 
    ! This set of operations done for both core and SOL ions.
    ! Check if on grid.
    If (Current_S .le. 0.0 .or. Current_S .ge. gksmaxs (Current_Ring)) Then
       ! Off grid, do not update cross.
       ! This will lead to the hc_outside_ion routine.
       ! jdemod - NOTE - this does not lead to the outside_ion routine - should it do so??
    Else
 
       ! jdemod GETRZ should be called for all values of rz_opt
       !If ( RZ_Opt .eq. 1) Then
          Call GETRZ (Current_Cell,Current_Ring,Current_S,Current_Cross,Current_R,Current_Z, RZ_Opt)
       !Else
       !   Current_R = grs (Current_Cell,Current_Ring)						
       !   Current_Z = gzs (Current_Cell,Current_Ring)						
       !End If
 
       ! psmod
       ! Record data on average forces acting on impurity particles.
       FTotal = ((FE + FEG + FIG + FF + FVG + DVPara) / QTIM_Squared) * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) *  AMU
 
       !  HC_FCell (Current_Cell,Current_Ring,Cur_HC_Spec) =  HC_FCell (Current_Cell,Current_Ring,Cur_HC_Spec) + FTotal * Sput_Weight
       ! jdemod - array was missing third dimension - also third dimension is scaled to maxizs - so is useless for 
       !          individual HC states - this code groups all HC fragments together by charge state - most of which are
       !          charge state 1. 
       HC_FCell (Current_Cell,Current_Ring,Get_HC_Charge (Cur_HC_Spec)) =  HC_FCell (Current_Cell,Current_Ring,& 
                                                                          & Get_HC_Charge (Cur_HC_Spec)) + FTotal * Sput_Weight
       ! Examine Frictional K11(FF), Thermal K12(FIG),
       ! and Background Velocity Gradient K13 (FVG) contributions seperately.
       HC_Ffi (Current_Cell,Current_Ring,Get_HC_Charge (Cur_HC_Spec)) = &
            &  HC_Ffi (Current_Cell,Current_Ring,Get_HC_Charge (Cur_HC_Spec)) + &
            & FF / QTIM_Squared * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) *  AMU * Sput_Weight
       HC_Fthi (Current_Cell,Current_Ring,Get_HC_Charge (Cur_HC_Spec)) = &
            &  HC_Fthi (Current_Cell,Current_Ring,Get_HC_Charge (Cur_HC_Spec)) + &
            & FIG / QTIM_Squared * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) *  AMU * Sput_Weight
       HC_Fvbg (Current_Cell,Current_Ring,Get_HC_Charge (Cur_HC_Spec)) = &
            &  HC_Fvbg (Current_Cell,Current_Ring,Get_HC_Charge (Cur_HC_Spec)) + &
            & FVG / QTIM_Squared * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) *  AMU * Sput_Weight
       HC_Force_Diff (Current_Cell,Current_Ring,Get_HC_Charge (Cur_HC_Spec)) = &
            &  HC_Force_Diff (Current_Cell, Current_Ring,Get_HC_Charge (Cur_HC_Spec)) + &
            & DVPara / QTIM_Squared * Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) *  AMU * Sput_Weight
       HC_IonVelAvg (Current_Cell,Current_Ring,Get_HC_Charge (Cur_HC_Spec)) = &
            &  HC_IonVelAvg (Current_Cell,Current_Ring,Get_HC_Charge (Cur_HC_Spec)) + &
            & Current_Velocity_In_S * Sput_Weight
       ! psmod
 
       ! (RIV)
       If ( Debug_HC_V) Then
          !  HC_SDVS (Current_Cell,Current_Ring,Cur_HC_Spec) =  HC_SDVS (Current_Cell,Current_Ring,Cur_HC_Spec)  + Sput_Weight * fvel
          !  HC_SDVS2 (Current_Cell,Current_Ring,Cur_HC_Spec) =  HC_SDVS2 (Current_Cell,Current_Ring,Cur_HC_Spec) + Sput_Weight * fvel**2.0
          HC_SDVS (Current_Cell,Current_Ring) =  HC_SDVS (Current_Cell,Current_Ring)  + Sput_Weight * fvel
          HC_SDVS2 (Current_Cell,Current_Ring) =  HC_SDVS2 (Current_Cell,Current_Ring) + Sput_Weight * fvel**2.0
 
          ! Note, SDVB needs to be recalculated here for CRMHC in the comparison.
          If(ABS(FVEL).gt.(9.79E3*SQRT(2.0*gktibs(Current_Cell,Current_Ring)/Back_Plasma_Ion_Mass)*Ion_Time_Step))Then
             HC_SDVS3 (Current_Cell,Current_Ring,1) =  HC_SDVS3 (Current_Cell,Current_Ring,1) + Sput_Weight * ABS (fvel)
             HC_SDVS3 (Current_Cell,Current_Ring,2) =  HC_SDVS3 (Current_Cell,Current_Ring,2) + Sput_Weight
             !write (0,*) "sdvb data:",ABS (FVEL),sput_Weight,9.79E3 * SQRT (2.0 * gktibs (Current_Cell,Current_Ring) / Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition)) *  Ion_Time_Step, &
             !&  HC_SDVS3 (Current_Cell,Current_Ring,1), HC_SDVS3 (Current_Cell,Current_Ring,2)
          End If
 
          ! Split up the velocity distribution and assign it to an appropriate bin.
          !Velocity_Bin = INT (FVEL / ( hc_VelSep *  hc_VelPlate))
          Velocity_Bin = INT (FVEL / ( hc_VelSep * 9.79E3 * SQRT (gktibs (1, Inner_SOL_Ring) / Find_HC_Mass (Cur_HC_Spec,&
&H_Isotope_Composition)) *  Ion_Time_Step))
          If (IABS (Velocity_Bin) .gt.  hc_NVel) Then
             Velocity_Bin = ISIGN ( hc_NVel,Velocity_Bin)
          End If
 
          If (Velocity_Bin .ge. 0 .and. FVEL .ge. 0.0) Then
             Velocity_Bin = Velocity_Bin + 1
          End If
 
          ! Set knot currently occupied.
          If ( Max_Velocity_Cells_Per_Ring .gt. gnks ( INJ_Ring_Number) .and. gnks ( INJ_Ring_Number) .gt. 0) Then
             Velocity_Cell = Current_Cell
          Else
             Velocity_Cell = 1
          End If
          HC_VelSpace (Velocity_Bin,Velocity_Cell) =  HC_VelSpace (Velocity_Bin,Velocity_Cell) + FVEL * Sput_Weight
          HC_VelWeight (Velocity_Bin,Velocity_Cell) =  HC_VelWeight (Velocity_Bin,Velocity_Cell) + Sput_Weight
       End If
 
       ! This will lead to the hc_inside_ion routine.
    End If
 
    ! Modify velocity to be true velocity, not distance per timestep.
    Current_Velocity_In_S = Current_Velocity_In_S /  Ion_Time_Step
 
  End Subroutine Ion_Move
 
End Module Hc_Ion_Transport
