! -*-Mode:f90-*-
! HC_Ouptut.f90
! Output routine specific to data from Hydrocarbon following.
!
! Adam McLean under Prof. Peter Stangeby
! Research Assistant: David Elder
! November, 2002
 
Module HC_Output
 
  ! Every good Fortran 90 program has...
  Implicit None	
 
Contains
 
  Subroutine HC_Print_Output (Pri_Neut_Launch,Sec_Neut_Launch,SExit,SMain,SHC,NeuTime)
 
    ! Required modules.
    Use ComHC ! Contains DIVIMP input file hydrocarbon-related options.
    Use HC_Init_DIV_Data ! Bloabl data.
    Use HC_Init_DIV_Diag ! Contains global initialization routines.
    Use HC_Init_Lib_Data ! Get HC mass and charge.
    Use HC_Utilities ! Contains data array initialization routine.
    Use HC_Get
    use hc_diag_data ! Diagnostic data from the hc routines - mostly reaction kinetics to start
 
    ! Every good Fortran program has...
    Implicit None
 
    ! Define input/output variables
    Integer, Intent (In) :: Pri_Neut_Launch
    Integer, Intent (In) :: Sec_Neut_Launch
    Real, Intent (In) :: SExit
    Real, Intent (In) :: SMain
    Real, Intent (In) :: SHC ! RNEUT, HCs launched.
    Real, Intent (In) :: NeuTime ! CPU time used to follow all HCs to completion.
 
    ! Define local variables that will be assigned values.
    Real, Dimension (Number_Regions) :: Particle_Launch_Sum
    Real, Dimension (Number_Regions) :: Particle_ReLaunch_Sum
    Real, Dimension (Number_Regions) :: Neutral_Launch_Sum
    Real, Dimension (Number_Regions) :: Ion_Launch_Sum
    Integer :: HC_Species
    Integer :: Current_HC_State
    Real :: Yield_Mult_Factor ! YMF, used to artificially decrease yield result.
    Real :: Temp_Sum ! TMPSUM
    Integer :: Launch_Reg ! M, Target (1-2), wall (3), FS PT (4), FS 2D (5), sputtered (6), reflected (7).
    Integer :: Current_Cell ! IK, Used to find average position of injected HCs.
    Integer :: Current_Ring ! IR, Used to find average position of injected HCs.
    Real :: EatIZTot ! EATIZTOT, Total temperature in both tagets.
    Real :: SatIZTot ! SATIZTOT, Total C+ particles produced.
    Real :: TemAV ! TEMAV, Average temperature of all C+ ions at production.
    Integer :: General_Counter
    Real :: HC_Factor ! FACT
    Double Precision :: Double_HC_Factor ! DACT
    Integer :: Temp_Counter
    Integer :: TimeStep ! IT
    Integer :: Target_Element ! ID
    ! jdemod - rename sdtimp to hc_sdtimp so the variable name doesn't overlap with ones in mod_diagvel
    Real, Dimension (maxnks,maxnrs) :: HC_SDTImp ! SDTIMP
    Integer :: Velocity_Cells ! IKV
    Integer :: Velocity_Bin ! IN
    Real :: Vel ! VEL
    Real :: Velocity_Total ! VELTOT
    Real, Dimension (- hc_NVel: hc_NVel+1) :: HC_VelCoord ! VELCOORD
    Integer :: Time_Bin ! IT
    Real :: DSum1
    Real :: DSum3
    Real :: DSum4
    Real, Dimension (-1:Number_HC_Species) :: SDtzs
    Real, Dimension (-1:Number_HC_Species) :: SDtzs2
    Integer :: Num_HCs ! Number of HCs being tracked.
    Integer, Dimension (20) :: Fates
    Integer :: Fate_Num
    Integer :: Reaction
    Real :: Reaction_Count_Sum
    Real, Dimension (Number_HC_Species,Number_HC_Species) :: HC_Transition_Table
    Integer :: Target_Species
    Integer :: Product_Count
    Integer :: All_Walls
    integer :: in ! generic counter for loops
    !
    Integer, Dimension (Number_H_Species) :: H_Isotope_Composition ! Note:  This is separate from the array throughout DIVIMP-HC in order to perform average background mass analysis in HC_Output.
    !Integer, Dimension (Number_H_Products) :: H_Isotope_Composition ! Note:  This is separate from the array throughout DIVIMP-HC in order to perform average background mass analysis in HC_Output.
    Real, Dimension (Number_HC_Species,Number_Regions) :: Tot_Reach_State_Per_Launch
    Character (Len = 6 * Highest_Carbon_Content) :: Full_Product_List ! Note: 5 characters is the current maximum for any single HC name.
    Integer :: HC_Charge
    Integer :: Particle
    Integer :: Launch_Pri_Sec
    Real :: Transitions
    Real :: Kin_Energy
    Real :: Area_Fact_CH
    !		Real :: Area_Fact_C2
    Real, Dimension (Number_HC_Species) :: Transition_Count_So_Far
    Real, Dimension (Number_HC_Species) :: Transition_Count_All_The_Way
    Real, Dimension (Number_HC_Species) :: Tot_Reach_State_All_The_Way
    Real, Dimension (Number_HC_Species) :: Time_All_The_Way
    Real, Dimension (Number_HC_Species) :: Energy_All_The_Way
    Real, Dimension (Number_HC_Species) :: Kin_Energy_All_The_Way
 
    ! Complete launch, follow and recording of all particles.  Print diagnostics.
 
    ! Find total of all HCs launched for all species and both targets.
    ! This figure is equivalent to TNEUT in DIVIMP.
    Particle_Launch_Sum = 0.0 ! Note: Array (Number_Regions).
    Neutral_Launch_Sum = 0.0
    Ion_Launch_Sum = 0.0
    H_Isotope_Composition = 0 ! Set to all 0's to perform plasma mass averaged analysis for mass.
 
    ! Assign number of hydrocarbon species available.
    num_hcs = number_hc_species

    write(0,*) 'H_Output: Assigning NUM_HCs =',number_hc_species

    ! jdemod - hard coding these is bad practice  
    !If (HC_Higher_HCs_Option .eq. 0) Then
    !   Num_HCs = 10
    !Else
    !   Num_HCs = 58
    !End If
 
    ! Perform all output for all regions which particles were launched from.
    Do Launch_Reg = 1,Number_Regions
       Do HC_Species = 1, Number_HC_Species
 
          If (Get_HC_Charge (HC_Species) .eq. 0) Then
             Neutral_Launch_Sum (Launch_Reg) = Neutral_Launch_Sum (Launch_Reg) +  HC_Sum_Fragments_Launched (HC_Species,Launch_Reg)
          Else
             Ion_Launch_Sum (Launch_Reg) = Ion_Launch_Sum (Launch_Reg) +  HC_Sum_Fragments_Launched (HC_Species,Launch_Reg)
          End If

          ! Sum up contribution from each target.
          Particle_Launch_Sum(Launch_Reg) = Particle_Launch_Sum(Launch_Reg) +  HC_Sum_Fragments_Launched(HC_Species,Launch_Reg)
          ! Krieger IPP/07 - SUN compiler insists on 132 column limit
          Particle_ReLaunch_Sum(Launch_Reg) = Particle_ReLaunch_Sum(Launch_Reg) + &
          &                                   SUM(HC_Sum_Fragments_ReLaunched(HC_Species,:,Launch_Reg))
 
          ! Check for zero.
          If ( HC_Sum_Fragments_Launched (HC_Species,Launch_Reg) .eq. 0.0) Then
             HC_Sum_Fragments_Launched (HC_Species,Launch_Reg) =  Calc_Lo
          End If
          Do All_Walls = 1, Max_Points
             If ( HC_Sum_Fragments_ReLaunched (HC_Species,All_Walls,Launch_Reg) .eq. 0.0) Then
                HC_Sum_Fragments_ReLaunched (HC_Species,All_Walls,Launch_Reg) =  Calc_Lo
             End If
          End Do
       End Do
 
       ! jdemod - why can't this be set to zero - if it causes a problem it might be hiding a bug
       ! Check for zero.
       If (Particle_Launch_Sum (Launch_Reg) .eq. 0.0) Then
          Particle_Launch_Sum (Launch_Reg) =  Calc_Lo * 1.0E6
       End If
 
    End Do
 
    ! Check that total particles launched equals sum of primary, secondary, and HC self-sputtered launched particles.
    ! jdemod - this usually does not seem to be consistent 
    If(INT(SUM(Particle_Launch_Sum)).ne.(Pri_Neut_Launch+Sec_Neut_Launch+INT(SUM(HC_Sum_Fragments_ReLaunched(:,:,:)))))Then
       ! Error - these should be equal.
       Write (Output_Unit_HC_Alert,*) "Warning in HC_Print_Output: Total particles launched is not consistent:", &
            & INT (SUM (Particle_Launch_Sum)),Pri_Neut_Launch,Sec_Neut_Launch,Pri_Neut_Launch+Sec_Neut_Launch,INT (SUM ( &
            & HC_Sum_Fragments_ReLaunched (:,:,:)))
    End If
 

    !
    ! jdemod - NOTE: since the HC state arrays are never indexed by either 0 or -1 these calculations are meainingless
    !
    ! Calculate FACTA and FACTB normalizations for bin sizes.
    !HC_Factor_A (0) = 0.0
    !If (SUM (Neutral_Launch_Sum) .gt. 0.0) Then
    !   HC_Factor_A (-1) = 1.0 / SUM (Neutral_Launch_Sum)
    !   HC_Factor_A (0) = 1.0 / SUM (Neutral_Launch_Sum)
    !ElseIf (SUM (Ion_Launch_Sum) .gt. 0.0) Then
    !   HC_Factor_A (-1) = 1.0 / SUM (Ion_Launch_Sum)
    !   HC_Factor_A(0) = 1.0 / SUM (Ion_Launch_Sum)
    !End If
    !HC_Factor_B (-1) =  HC_Factor_A (-1) *  Neutral_Time_Step
    !HC_Factor_B (0) =  HC_Factor_A (0) *  Neutral_Time_Step
    
    write(0,'(a,8(1x,g12.5))') 'HC particles launched:',sum(neutral_launch_sum), sum(ion_launch_sum),&
                                                      & sum(particle_launch_sum),sum(particle_relaunch_sum),&
                                                      & pri_neut_launch, sec_neut_launch
    !write(0,*) 'HC particle timesteps:',neutral_time_step,ion_time_step
    
 
    ! jdemod - neither of these values appears to be used or printed - in addition they 
    !          are incorrectly indexed using an index of 0 
    !        - the actual arrays appear to be indexed from -1 though I can find no use of this
    !          Also - what is the difference between entering "main" and entering "core" - often
    !          these mean the same thing and it isn't clear that both are needed here
    !
    !HC_Num_Orig_Enter_Main (0) = SExit
    !		 HC_Cmmm (0) = 0.0
    !HC_Num_Orig_Enter_Core (0) = SMain
 
    Do HC_Species = 1, Number_HC_Species
       ! Currently we do not discriminate between ions injected and neutrals launched.
       If (SUM (Ion_Launch_Sum) .gt. 0.0) Then
          HC_Factor_A (HC_Species) = 1.0 / SUM (Particle_Launch_Sum)
       Else
          HC_Factor_A (HC_Species) = 1.0 / SUM (Particle_Launch_Sum)
       End If
 
       ! Multiply by applicable timestep.
       If (Get_HC_Charge (HC_Species) .eq. 0.0) Then
          HC_Factor_B (HC_Species) =  HC_Factor_A (HC_Species) *  Neutral_Time_Step
       Else
          HC_Factor_B (HC_Species) =  HC_Factor_A (HC_Species) *  Ion_Time_Step
       End If
    End Do
 
    ! DEAL WITH ANOMALY FOR CONTINUOUS RINGS, WHERE TWO POINTS ON THE
    ! RING ARE COINCIDENT - COMBINE FIRST & LAST POINTS ON RING.
    Do HC_Species = 1, Number_HC_Species
       HC_ELims (1,3,HC_Species) =  HC_ELims (1,3,HC_Species) +  HC_ELims (gNKS ( Inner_SOL_Ring - 1), 3, HC_Species)
       HC_ELims (gNKS ( Inner_SOL_Ring - 1), 3, HC_Species) = 0.0
       Do Current_Ring = 1,  Inner_SOL_Ring - 1
 
          Do Launch_Pri_Sec = -1,0
             HC_TIZS (1,Current_Ring,Launch_Pri_Sec) =  HC_TIZS (1,Current_Ring,Launch_Pri_Sec) +  &
                                                      & HC_TIZS (gnks(Current_Ring),Current_Ring,Launch_Pri_Sec)
             HC_TIZS (1,Current_Ring,Launch_Pri_Sec) = 0.0
          End Do
          Do Launch_Pri_Sec = -1,0 ! -1=primaries only, 0=total.
             HC_TIZS_CH (1,Current_Ring,Launch_Pri_Sec) =  HC_TIZS_CH (1,Current_Ring,Launch_Pri_Sec) +  &
                                                         & HC_TIZS_CH (gnks(Current_Ring),Current_Ring,Launch_Pri_Sec)
             HC_TIZS_CH (1,Current_Ring,Launch_Pri_Sec) = 0.0
             !					 HC_TIZS_C2 (1,Current_Ring,Launch_Pri_Sec) =  HC_TIZS_C2 (1,Current_Ring,Launch_Pri_Sec) +  HC_TIZS_C2 (gnks(Current_Ring),Current_Ring,Launch_Pri_Sec)
             !					 HC_TIZS_C2 (1,Current_Ring,Launch_Pri_Sec) = 0.0
          End Do
 

          !
          ! jdemod - change the treatment of hc_ddts to match that of hc_density - this was inconsistent throughout the hc code 
          !          eliminate hc_ddlims
          !
          HC_DDts (1,Current_Ring,hc_species) =  HC_DDts (1,Current_Ring,hc_species) + &
                                              &  HC_DDts (gNKS (Current_Ring),Current_Ring,hc_species)
          !HC_DDLims (1,Current_Ring,0) =  HC_DDLims (1,Current_Ring,0) +  HC_DDLims (gNKS (Current_Ring),Current_Ring,0)
! slmod
! Apparently, HC_DDTs is declaired HC_DDts (gNKS (?,?,1) somewhere, and so can't take a
! zero in the 3rd index.  Sorry, can't figure this one out at the moment so commenting out...
!
!          HC_DDts (gNKS (Current_Ring),Current_Ring,0) = 0.0D0
! slmod end
          !HC_DDLims (gNKS (Current_Ring),Current_Ring,0) = 0.0D0
 
          HC_Density (1,Current_Ring,HC_Species) =  HC_Density (1,Current_Ring,HC_Species) + &
                                                  & HC_Density (gNKS (Current_Ring), Current_Ring,HC_Species)
          HC_Density (gNKS (Current_Ring),Current_Ring,HC_Species) = 0.0D0
 
          Do TimeStep = 1,  Number_Time_Steps
             HC_Lims (1, Current_Ring, HC_Species, TimeStep) =  HC_Lims (1, Current_Ring, HC_Species, TimeStep) + &
                                           & HC_Lims (gNKS (Current_Ring), Current_Ring, HC_Species, TimeStep)
             HC_Lims (gNKS (Current_Ring), Current_Ring, HC_Species, TimeStep) = 0.0
          End Do
 
          ! (RIV)
 
          HC_SDVS (1,Current_Ring)  =  HC_SDVS (1,Current_Ring) +  HC_SDVS (gnks (Current_Ring), Current_Ring)
          HC_SDVS (gnks (Current_Ring),Current_Ring) = 0.0
 
          If ( Debug_HC_V) Then
 
             HC_SDVS2 (1,Current_Ring)  =  HC_SDVS2 (1,Current_Ring) +  HC_SDVS2 (gnks (Current_Ring), Current_Ring)
             HC_SDVS2 (gnks (Current_Ring), Current_Ring) = 0.0
             HC_SDVS3 (1,Current_Ring,1) =  HC_SDVS3 (1,Current_Ring,1) +  HC_SDVS3 (gnks (Current_Ring), Current_Ring,1)
             HC_SDVS3 (gnks (Current_Ring), Current_Ring, 1) = 0.0
             HC_SDVS3 (1, Current_Ring, 2) =  HC_SDVS3 (1, Current_Ring, 2) +  HC_SDVS3 (gnks (Current_Ring), Current_Ring, 2)
             HC_SDVS3 (gnks (Current_Ring), Current_Ring, 2) = 0.0
 
          End If
       End Do
    End Do
 
    Do HC_Species = 1, Number_HC_Species
       Do Current_Ring = 1,  Num_Upper_Rings ! NRS.
          Do Current_Cell = 1, gNKS (Current_Ring)
 
             If (gkareas(Current_Cell, Current_Ring) .ne. 0.0) Then
                Double_HC_Factor = DBLE ( HC_Factor_B (HC_Species) / gkareas (Current_Cell, Current_Ring))
             Else
                Double_HC_Factor = 0.0
             End If
 
             !HC_DDLims (Current_Cell,Current_Ring,INT (Get_HC_Charge (HC_Species))) = Double_HC_Factor *  HC_DDLims (Current_Cell,&
             !                                                             &Current_Ring,INT(Get_HC_Charge (HC_Species)))
             HC_Density(Current_Cell,Current_Ring,HC_Species)=Double_HC_Factor*HC_Density(Current_Cell,Current_Ring,HC_Species)
 
             Do TimeStep = 1,  Number_Time_Steps
                HC_Lims (Current_Cell,Current_Ring,HC_Species,TimeStep) = Double_HC_Factor *  HC_Lims (Current_Cell,Current_Ring,&
&HC_Species,TimeStep) &
                     & / (gctimes (TimeStep,Get_HC_Charge (HC_Species))-gctimes (TimeStep-1,Get_HC_Charge (HC_Species)))
             End Do
 
             If (Get_HC_Charge (HC_Species) .eq. 0) Then
                HC_ChemDen(Current_Cell,Current_Ring,HC_Species)=Double_HC_Factor*HC_ChemDen(Current_Cell,Current_Ring,HC_Species)
             End If
          End Do
       End Do
    End Do
 
    ! Ionization of CH and C2 only.
    Do Launch_Pri_Sec = -1,0 ! Do for -1-> primaries only, and 0->including re-releases
       Do Current_Ring = 1,  Num_Upper_Rings
          Do Current_Cell = 1, gnks (Current_Ring)
             If (gkareas (Current_Cell,Current_Ring) .ne. 0.0) Then
                Area_Fact_CH =  HC_Factor_A (HC_Species_Ident("CH")) / gkareas (Current_Cell,Current_Ring)
                !						Area_Fact_C2 =  HC_Factor_A (HC_Species_Ident("C2")) / gkareas (Current_Cell,Current_Ring)
             Else
                Area_Fact_CH = 0.0
                !						Area_Fact_C2 = 0.0
             EndIf
 
             HC_TIZS_CH(Current_Cell,Current_Ring,Launch_Pri_Sec)=Area_Fact_CH*HC_TIZS_CH(Current_Cell,Current_Ring,Launch_Pri_Sec)
             !					 HC_TIZS_C2 (Current_Cell,Current_Ring,Launch_Pri_Sec) = Area_Fact_C2 *  HC_TIZS_C2 (Current_Cell,Current_Ring,Launch_Pri_Sec)
          End Do
       End Do
    End Do
 
    ! Modify leaked particle sum to account for start at value of 1.
    Do Launch_Reg = 1,Number_Regions
       If ( HC_Leak_Particles (Launch_Reg) .gt. 0) Then
          HC_Leak_Particles (Launch_Reg) =  HC_Leak_Particles (Launch_Reg) - 1
       End If
    End Do
 
    ! Modify states reached sum to prevent nans.
    Do Launch_Reg = 1,Number_Regions
       Do HC_Species = 1,Num_HCs
          Tot_Reach_State_Per_Launch (HC_Species,Launch_Reg) =  HC_Tot_Reach_State (HC_Species,Launch_Reg)
          If (Tot_Reach_State_Per_Launch (HC_Species,Launch_Reg) .eq. 0.0) Then
             Tot_Reach_State_Per_Launch (HC_Species,Launch_Reg) = 1E-7
          End If
       End Do
    End Do
 
    ! Print Table of Contents for HC data output.
    Write (Output_Unit_HC_Data,'(A)') ""
    Write (Output_Unit_HC_Data,'(A)') "Table of Contents"
    Write (Output_Unit_HC_Data,'(A)') "1) Particle LAUNCH Statistics"
    Write (Output_Unit_HC_Data,'(A)') "2) Particle TRANSPORT Statistics"
    Write (Output_Unit_HC_Data,'(A)') "3) Particle EVOLUTION Statistics"
    Write (Output_Unit_HC_Data,'(A)') "4) Particle RE-RELEASE Statistics"
    Write (Output_Unit_HC_Data,'(A)') "5) Particle VESSEL INTERACTION Statistics"
    Write (Output_Unit_HC_Data,'(A)') "6) Particle END-OF-LIFE Statistics"
 
    ! Print LAUNCH diagnostics.
    Write (Output_Unit_HC_Data,'(A)') ""
    Write (Output_Unit_HC_Data,'(A)') "1) Particle LAUNCH Statistics"
    Write (Output_Unit_HC_Data,'(A)') ""
    Write (Output_Unit_HC_Data,'(1X,A,19X,A,2X,200A12)') "Region specific launch statistics","Unit",( HC_Region_Names (Launch_Reg),&
& Launch_Reg=1,Number_Regions,1)
    Write(Output_Unit_HC_Data,27)"Number of particles launched","#",(Particle_Launch_Sum(Launch_Reg),Launch_Reg=1,Number_Regions,1)
    Write (Output_Unit_HC_Data,27) "Number of particles relaunched before reaching C+","#",(Particle_ReLaunch_Sum (Launch_Reg), &
&Launch_Reg=1,Number_Regions,1)
    Write (Output_Unit_HC_Data,28) "Fraction of particles launched","n/a",(Particle_Launch_Sum (Launch_Reg) / SUM (&
&Particle_Launch_Sum),Launch_Reg=1,Number_Regions,1)
    Write(Output_Unit_HC_Data,27)"Number of neutrals launched","#",(Neutral_Launch_Sum(Launch_Reg),Launch_Reg=1,Number_Regions,1)
    Write (Output_Unit_HC_Data,27) "Number of ions launched","#",(Ion_Launch_Sum (Launch_Reg), Launch_Reg=1,Number_Regions,1)
 
27  Format (1X,A50,2X,A4,2X,200F12.1) ! Good for values up to 999,999.9
28  Format (1X,A50,2X,A4,2X,200F12.6) ! Good for fractions up to 1.000000
 
    Do Launch_Reg = 1,Number_Regions
 
       ! Check that molecules were launched from this region.
       If ( HC_Region_Used (Launch_Reg) .gt. 0) Then
 
          ! Print headings.
          Write (Output_Unit_HC_Data,'(A)') ""
          Write (Output_Unit_HC_Data,'(1X,A,1X,A10,26X,A)') "Species and region specific launch statistics for launches:", &
&HC_Region_Names (Launch_Reg),"Species Launched"
          Write (Output_Unit_HC_Data,'(42X,A9,2X,A4,2X,100A10)') "Statistic","Unit",(HC_State_Table (HC_Species) % State_Name,&
&HC_Species = 1,Num_HCs,1)
 
          Write (Output_Unit_HC_Data,40) "Average R launch position","m",( HC_Total_R_Prod_Positions (HC_Species,Launch_Reg) / &
&Particle_Launch_Sum (Launch_Reg), HC_Species = 1,Num_HCs,1)
          Write (Output_Unit_HC_Data,40) "Average Z launch position","m",( HC_Total_Z_Prod_Positions (HC_Species,Launch_Reg) / &
&Particle_Launch_Sum (Launch_Reg), HC_Species = 1,Num_HCs,1)
          Write (Output_Unit_HC_Data,40) "Average S launch position","m",( HC_Total_S_Prod_Positions (HC_Species,Launch_Reg) / &
&Particle_Launch_Sum (Launch_Reg), HC_Species = 1,Num_HCs,1)
          Write (Output_Unit_HC_Data,40) "Average cross launch position","m",( HC_Total_Cross_Prod_Positions (HC_Species,&
&Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1,Num_HCs,1)
          Write (Output_Unit_HC_Data,40) "Average abs angle from normal at launch","deg",( HC_Total_Prod_Angles (HC_Species,&
&Launch_Reg) *raddeg / Particle_Launch_Sum (Launch_Reg), HC_Species = 1,Num_HCs,1)
          Write (Output_Unit_HC_Data,40) "Fraction of species launched","n/a",( HC_Sum_Fragments_Launched (HC_Species,Launch_Reg) /&
& Particle_Launch_Sum (Launch_Reg), HC_Species = 1,Num_HCs,1)
 
          Write (Output_Unit_HC_Data,42) "Average velocity at launch no vmult","m/s",( HC_Total_Prod_Vels_No_VMult (HC_Species,&
&Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1,Num_HCs,1)
          Write (Output_Unit_HC_Data,42) "Maximum velocity at launch no vmult","m/s",( HC_Max_Prod_Vel_No_VMult (HC_Species,&
&Launch_Reg), HC_Species = 1,Num_HCs,1)
          !
          !                               jdemod - the following line wasn't printing the right information in any case
          !
          !				Write (Output_Unit_HC_Data,42) "Average velocity/angle vmult","#",( HC_Total_Prod_Vels_No_VMult (HC_Species,Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1,Num_HCs,1)
          !				Write (Output_Unit_HC_Data,40) "Maximum velocity/angle vmult","#",( HC_Max_Vel_Ang_Mults (HC_Species,Launch_Reg), HC_Species = 1,Num_HCs,1)
          !
          Write (Output_Unit_HC_Data,42) "Average velocity at launch","m/s",( HC_Total_Prod_Velocities (HC_Species,Launch_Reg) / &
&Particle_Launch_Sum (Launch_Reg), HC_Species = 1,Num_HCs,1)
          Write (Output_Unit_HC_Data,42) "Maximum velocity at launch","m/s",( HC_Max_Prod_Velocities (HC_Species,Launch_Reg), &
&HC_Species = 1,Num_HCs,1)
 
          Write (Output_Unit_HC_Data,40) "Average temperature at launch","eV",( HC_Tot_Prod_Temperatures (HC_Species,Launch_Reg) / &
&Particle_Launch_Sum (Launch_Reg), HC_Species = 1,Num_HCs,1)
          Write (Output_Unit_HC_Data,42) "Number of failed launches","#",( HC_Num_Failed_Launches (HC_Species,Launch_Reg), &
&HC_Species = 1,Num_HCs,1)
          Write (Output_Unit_HC_Data,42) "Number launch velocities > maximum","#",( HC_Num_Vel_Greater_Max_Ran (HC_Species,&
&Launch_Reg), HC_Species = 1,Num_HCs,1)
          Write (Output_Unit_HC_Data,42) "Launch grid error moved okay","#",( HC_Launch_Grid_Error_Moved_Okay (HC_Species,&
&Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1,Num_HCs,1)
          Write (Output_Unit_HC_Data,42) "Launch grid error after 1 timestep","#",( HC_Launch_Grid_Error_Moved_Out (HC_Species,&
&Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1,Num_HCs,1)
          !				Write (Output_Unit_HC_Data,42) "Launchdat 1","#",( HC_LaunchDat (HC_Species,1) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1,Num_HCs,1)
          !				Write (Output_Unit_HC_Data,42) "Launchdat 2","#",( HC_LaunchDat (HC_Species,2) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1,Num_HCs,1)
          !				Write (Output_Unit_HC_Data,42) "Launchdat 3","#",( HC_LaunchDat (HC_Species,3) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1,Num_HCs,1)
          !				Write (Output_Unit_HC_Data,42) "Launchdat 4","#",( HC_LaunchDat (HC_Species,4) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1,Num_HCs,1)
       End If
 
40     Format (1X,A50,2X,A4,1X,100F10.4) ! Good for values up to 999.9999
41     Format (1X,A50,2X,A4,1X,100F10.2) ! Good for values up to 99,999.99
42     Format (1X,A50,2X,A4,1X,100F10.1) ! Good for values up to 999,999.9
43     Format (1X,A50,2X,A4,1X,100F10.6) ! Good for fractions up to 1.000000
 
    End Do
 
    ! Print TRANSPORT diagnostics.
    Write (Output_Unit_HC_Data,'(A)') ""
    Write (Output_Unit_HC_Data,'(A)') "2) Particle TRANSPORT Statistics"
    Write (Output_Unit_HC_Data,'(A)') ""
 
    ! Print single valued timing data.
    Write (Output_Unit_HC_Data,'(1X,A50,7X,F8.1,1X,A1)') "Total CPU time following all particles",NeuTime,"s"
    Write(Output_Unit_HC_Data,'(1X,A50,7X,F8.5,1X,A1)')"Average CPU time per particle launched",NeuTime/SUM(Particle_Launch_Sum),"s"
    Write (Output_Unit_HC_Data,'(A)') ""
 
    ! Print single-valued data for each target.
    ! Timing data.
    Write (Output_Unit_HC_Data,'(1X,A,19X,A,2X,200A12)') "Region specific launch statistics","Unit",( HC_Region_Names (Launch_Reg),&
& Launch_Reg=1,Number_Regions,1)
    Write (Output_Unit_HC_Data,156) "Average time of flight per particles","s",( Time_All_HCs (Launch_Reg) / Particle_Launch_Sum (&
&Launch_Reg), Launch_Reg=1,Number_Regions,1)
    Write (Output_Unit_HC_Data,156) "Max flight time for any single particle","s",( Max_Time_Any_HC (Launch_Reg), Launch_Reg=1,&
&Number_Regions,1)
    Write (Output_Unit_HC_Data,50) "Average Ion equivalent timesteps in lifetime","#",( HC_Total_Ion_Teq_Iter (Launch_Reg) / &
&Particle_Launch_Sum (Launch_Reg), Launch_Reg=1,Number_Regions,1)
    Write (Output_Unit_HC_Data,50) "Maximum Ion equivalent timesteps in lifetime","#",( HC_Max_Ion_Teq_Iter (Launch_Reg), &
&Launch_Reg=1,Number_Regions,1)
    Write (Output_Unit_HC_Data,154) "Average time to first diffusion","s",( HC_Tot_Time_First_Diff (Launch_Reg) / &
&Particle_Launch_Sum (Launch_Reg), Launch_Reg=1,Number_Regions,1)
    Write (Output_Unit_HC_Data,50) "Average neutral timesteps per particle","#",( HC_Neutral_Timesteps_Count (Launch_Reg) / &
&Particle_Launch_Sum (Launch_Reg), Launch_Reg=1,Number_Regions,1)
    Write (Output_Unit_HC_Data,50) "Average ion timesteps per particle","#",( HC_Ion_Timesteps_Count (Launch_Reg) / &
&Particle_Launch_Sum (Launch_Reg), Launch_Reg=1,Number_Regions,1)
    Write (Output_Unit_HC_Data,'(A)') ""
 
    ! Leakage data.
    Write (Output_Unit_HC_Data,'(1X,A,18X,A,2X,200A12)') "Region specific leakage statistics","Unit",( HC_Region_Names (Launch_Reg)&
&, Launch_Reg=1,Number_Regions,1)
    Write (Output_Unit_HC_Data,52) "Number of ions leaked into core","#",( Number_Leaked_Core (Launch_Reg), Launch_Reg=1,&
&Number_Regions,1)
    Write(Output_Unit_HC_Data,50)"Total ion leakage from core","#",(Total_Leaked_Core(Launch_Reg),Launch_Reg=1,Number_Regions,1)
    Write (Output_Unit_HC_Data,156) "Average time to ion leakage","s",( HC_Leak_Time (Launch_Reg) / Particle_Launch_Sum (&
&Launch_Reg),Launch_Reg=1,Number_Regions,1)
    Write (Output_Unit_HC_Data,51) "Fraction of leaked particles","n/a",( HC_Leak_Particles (Launch_Reg) / Particle_Launch_Sum (&
&Launch_Reg),Launch_Reg=1,Number_Regions,1)
    Write (Output_Unit_HC_Data,51) "Average launch R position of leaked particles","m",(SUM ( HC_Leak_Position (:,Launch_Reg,1)) / &
         & ( Number_Leaked_Core (Launch_Reg) +  Calc_Lo),Launch_Reg=1,Number_Regions,1)
    Write (Output_Unit_HC_Data,51) "Average launch Z position of leaked particles","m",(SUM ( HC_Leak_Position (:,Launch_Reg,2)) / &
         & ( Number_Leaked_Core (Launch_Reg) +  Calc_Lo),Launch_Reg=1,Number_Regions,1)
    Write (Output_Unit_HC_Data,'(A)') ""
 
! jdemod - added a few more decimal places for some and switched to exponential notation 
50  Format (1X,A50,2X,A4,2X,200E12.4) ! Good for values up to 999,999.9
51  Format (1X,A50,2X,A4,2X,200F12.6) ! Good for fractions/times up to 1.000000
52  Format (1X,A50,2X,A4,2X,200I12) ! Good for numbers up to 10,000,000
53  Format (1X,A50,2X,A4,2X,200F12.4) ! Good for reals up to 999.9999
154 Format (1X,A50,2X,A4,2X,200E12.4) ! Good for values up to 0.999E99
155 Format (1X,A50,2X,A4,2X,200F12.8) ! Good for fractions/times up to 1.00000000 
! jdemod - adding some formats that will do a better job on very small numbers
156 Format (1X,A50,2X,A4,2X,200(1x,E11.5)) ! Good for fractions/times up to 1.00000000 

    ! Print arrayed data.
    Do Launch_Reg = 1,Number_Regions
 
       ! Check that molecules were launched from this region.
       If ( HC_Region_Used (Launch_Reg) .gt. 0) Then
 
          ! Print headings.
          Write (Output_Unit_HC_Data,'(1X,A,1X,A10,25X,A)') "Species and region specific transport statistics for launches:", &
&HC_Region_Names (Launch_Reg),"Species Followed"
 
          ! Timing statistics.
          Write (Output_Unit_HC_Data,'(1X,A,35X,A,1X,100A11)') "Timing statistics","Unit",(HC_State_Table (HC_Species) % &
&State_Name,HC_Species = 1,Num_HCs,1)
          Write (Output_Unit_HC_Data,57) "Total particles that reach state","#",(Tot_Reach_State_Per_Launch (HC_Species,Launch_Reg)&
&, HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,541) "Average time spent in each state when reached","s",( HC_State_Time (HC_Species,&
&Launch_Reg) / (Tot_Reach_State_Per_Launch (HC_Species,Launch_Reg)), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,541) "Max time spent in each state when reached","s",( HC_Max_Time_In_State (HC_Species,&
&Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,57) "Average timesteps spent in each state when reached","s",( HC_State_TimeSteps (HC_Species,&
&Launch_Reg) / (Tot_Reach_State_Per_Launch (HC_Species,Launch_Reg)), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,57) "Max timesteps spent in each state","s",( HC_Max_TimeSteps (HC_Species,Launch_Reg), &
&HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,541) "Fraction of particles reaching time cutoff","n/a",( HC_Num_Reach_Max_Iter (HC_Species,&
&Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,57) "Particles existing at TMax","#",( HC_Tot_Exist_Max_Iter (HC_Species,Launch_Reg), &
&HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,541) "Total time spent in each state","s",( HC_State_Time (HC_Species,Launch_Reg), HC_Species &
&= 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,'(A)') ""
 
          ! State sums.
          Write (Output_Unit_HC_Data,'(1X,A,25X,A,X,100A11)') "State occurrence statistics","Unit",(HC_State_Table (HC_Species) % &
&State_Name,HC_Species = 1,Num_HCs,1)
          Write (Output_Unit_HC_Data,55) "Particle fraction reaching state at least once","n/a",( HC_Num_Orig_Reach_State (&
&HC_Species,Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,55) "Average times each particle reaches each state","#",( HC_Tot_Reach_State (HC_Species,&
&Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,'(A)') ""
 
          ! Position sums.
          Write (Output_Unit_HC_Data,'(1X,A,33X,A,X,100A11)') "Position statistics","Unit",(HC_State_Table (HC_Species) % &
&State_Name,HC_Species = 1,Num_HCs,1)
          ! jdemod - change to number from average
          !Write (Output_Unit_HC_Data,541) "Avg num of particles striking target per launch","n/a",( HC_Num_Striking_Target (&
          !       &HC_Species,Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,541) "Num of particles striking target per launch","n/a",( HC_Num_Striking_Target (&
                &HC_Species,Launch_Reg), HC_Species = 1, Num_HCs, 1)
          ! jdemod

          ! IPP/09 Krieger - shortened lines (SUNWorkshop chokes over len>132)

          Write(Output_Unit_HC_Data,541) "Avg num of particles enter main plasma per launch","n/a",( HC_Num_Enter_Main_Plasma (&
&HC_Species,Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write(Output_Unit_HC_Data,541) "Avg num of particles exit main plasma per launch","n/a",( HC_Tot_Fragments_Exit_Main (&
&HC_Species,Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write(Output_Unit_HC_Data,541) "Avg num of particles entering core per launch","n/a",( HC_Num_Entered_Core (HC_Species,&
&Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write(Output_Unit_HC_Data,541) "Avg num of particles reaching centre per launch","n/a",( HC_Num_Reach_Centre (HC_Species,&
&Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write(Output_Unit_HC_Data,541) "Avg num of ions stoppped follow in core per launch","n/a",( &
&HC_Stopped_Follow_Ion_In_Core (HC_Species,Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
 
          ! Void region statistics.
          Write (Output_Unit_HC_Data,57) "Avg TSs for neuts in main plasma void per launch","#",( HC_DDVoid (1,HC_Species,&
&Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,57) "Avg TSs for neuts in priv plasma void per launch","#",( HC_DDVoid (2,HC_Species,&
&Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,57) "Avg TSs for neuts in other div void per launch","#",( HC_DDVoid (3,HC_Species,Launch_Reg)&
& / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,'(A)') ""
 
          ! Distance sums.
          Write (Output_Unit_HC_Data,'(1X,A,14X,A,X,100A11)') "Distance travelled as state statistics","Unit",(HC_State_Table (&
&HC_Species) % State_Name,HC_Species = 1,Num_HCs,1)
          Write (Output_Unit_HC_Data,546) "Average cumulative real distance by position","m", &
               & ( HC_Cumu_Real_Dist_Trav_By_Pos (HC_Species,Launch_Reg) / ( HC_Tot_Reach_State (HC_Species,Launch_Reg) +  Calc_Lo)&
&, HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,546) "Average cumulative real distance by vel*time","m", &
               & ( HC_Cumu_Real_Dist_Trav_By_VT (HC_Species,Launch_Reg) / ( HC_Tot_Reach_State (HC_Species,Launch_Reg) +  Calc_Lo),&
& HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,546) "Average state distance by position","m", &
               & ( HC_State_Dist_Trav_By_Pos (HC_Species,Launch_Reg) / ( HC_Tot_Reach_State (HC_Species,Launch_Reg) +  Calc_Lo), &
&HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,546) "Average state distance by vel*time","m", &
               & ( HC_State_Dist_Trav_By_VT (HC_Species,Launch_Reg) / ( HC_Tot_Reach_State (HC_Species,Launch_Reg) +  Calc_Lo), &
&HC_Species = 1, Num_HCs, 1) 	
          Write (Output_Unit_HC_Data,57) "Number ending their life in this state","#",((SUM ( HC_IFate_Count (HC_Species,&
&Launch_Reg,:)) +  Calc_Lo), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,546) "Average as-the-crow-flies dist to death by pos","m", &
               & ( HC_Cumu_CroFly_Dist_Trav_By_Pos (HC_Species,Launch_Reg) / (SUM ( HC_IFate_Count (HC_Species,Launch_Reg,:)) +  &
&Calc_Lo), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,'(A)') ""
 
          ! Far-periphery statistics.
          Write (Output_Unit_HC_Data,'(1X,A,28X,A,X,100A11)') "Far-periphery statistics","Unit",(HC_State_Table (HC_Species) % &
&State_Name,HC_Species = 1,Num_HCs,1)
          Write (Output_Unit_HC_Data,57) "FPTarg","#",( HC_FPTarg (HC_Species,Launch_Reg) / Particle_Launch_Sum (Launch_Reg), &
&HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,54) "Fraction of ions lost to far periphery target","n/a",( HC_FPTart (HC_Species,Launch_Reg) &
&/ Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,54) "Fraction of ions plating on far periphery target","n/a",( HC_RFPTarg (HC_Species,&
&Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,54) "Fraction of ions entering far periphery","n/a",( HC_FPEnt (HC_Species,Launch_Reg) / &
&Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,54) "Fraction of ion exiting far periphery","n/a",( HC_FPExit (HC_Species,Launch_Reg) / &
&Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,54) "Fraction of ions existing to TMax in far periphery","n/a",( HC_FPTTotal (HC_Species,&
&Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,54) "Ion fraction entering main from FP w reflections","n/a",( HC_Tot_Core_From_FP_Ref (&
&HC_Species,Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,54) "Ion fraction entering main from FP w/o reflections","n/a",( HC_Tot_Core_From_FP_NoRef (&
&HC_Species,Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,'(A)') ""
 
          ! Movement indicators.
          Write (Output_Unit_HC_Data,'(1X,A,33X,A,X,100A11)') "Movement statistics","Unit",(HC_State_Table (HC_Species) % &
                & State_Name,HC_Species = 1,Num_HCs,1)
          Write (Output_Unit_HC_Data,54) "Fraction of orig particles entering main plasma","n/a",( HC_Num_Orig_Enter_Main (&
                & HC_Species,Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,55) "Average Z value of particles entering main plasma","m",( HC_Tot_Z_Orig_Enter_Main (&
                & HC_Species,Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,55) "Average S value of particles entering main plasma","m",( HC_Tot_S_Orig_Enter_Main (&
                & HC_Species,Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,54) "Fraction of particles entering core plasma","n/a",( HC_Num_Orig_Enter_Core (HC_Species,&
                & Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,55) "Average Z value of particles entering core plasma","m",( HC_Tot_Z_Orig_Enter_Core (&
                & HC_Species,Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,55) "Average S value of particles entering core plasma","m",( HC_Tot_S_Orig_Enter_Core (&
                & HC_Species,Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,57) "Average ion/neut eqiv ts's at orig core entry","#",( HC_Tot_Teq_Orig_Enter_Core (&
                & HC_Species,Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,55) "Average minimum Z value reached","m",( HC_Tot_Min_Z_Reach (HC_Species,Launch_Reg) / &
                & Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,55) "Average maximum S or SMAX-S value reached","m",( HC_Tot_Max_S_Reach (HC_Species,&
                & Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,55) "Average max S or SMAX-S reached for removed ions","m",( HC_Tot_Max_S_Reach_Ion_Removal (&
                & HC_Species,Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,54) "Regular launch fraction enter main w reflections","n/a",( HC_Tot_Core_From_Reg_Ref (&
                & HC_Species,Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,54) "Regular launch fraction enter main w/o reflections","n/a",( HC_Tot_Core_From_Reg_NoRef (&
                & HC_Species,Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,55) "Average R at species creation","m",( HC_Tot_R_Species_Creation (HC_Species,Launch_Reg) / &
                & Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,55) "Average Z at particle entry to main plasma","m",( HC_Tot_Z_At_Main_Entry (HC_Species,&
                & Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,54) "Fraction of launched particles reflecting off central mirror","n/a",( &
            & HC_Num_Orig_Central_Reflect (HC_Species,Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,57) "Average ion-eqivalent TSs to reflection off central mirror","#",( &
          & HC_Tot_Teq_Orig_Central_Reflect (HC_Species,Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,56) "Average temperature at central mirror reflection","eV",( HC_Tot_Temp_At_Central_Reflect (&
               & HC_Species,Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,57) "Minimum ion-eqivalent timesteps for reflection off central mirror","#",( &
           & HC_Min_Teq_Orig_Central_Reflect (HC_Species,Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,58) "Average number of ion collisions per launched particle","#",( HC_Tot_Ion_Collision (&
               & HC_Species,Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,54) "Fraction of particles starting in the main plasma","n/a",( HC_Tot_Start_Main (HC_Species,&
&Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,57) "Min ion-equivalent TSs to entry to main plasma","#",( HC_Min_Teq_Orig_Enter_Main (&
&HC_Species,Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,57) "Avg ion-equivalent TSs to entry to main plasma","#",( HC_Tot_Teq_Orig_Enter_Main (&
&HC_Species,Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,55) "Average temperature at entrance to main plasma","eV",( HC_Tot_Temp_Orig_Enter_Main (&
&HC_Species,Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,55) "Fraction of particles reaching state being lost","n/a",( HC_Num_Ions_Lost (HC_Species,&
&Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,57) "Minimum ion-equivalent timesteps to ion loss","#",( HC_Tot_Teq_Ions_Lost (HC_Species,&
&Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,57) "Maximum ion-equivalent timesteps to ion loss","#",( HC_Tot_Teq_Ions_Lost (HC_Species,&
&Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,57) "Average ion-equivalent timesteps to ion loss","#",( HC_Tot_Teq_Ions_Lost (HC_Species,&
&Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          ! Neutral MTC counting.
          Write (Output_Unit_HC_Data,57) "Momentum Transfer Collision struck","#",( HC_MTC_Striking_Target (HC_Species,Launch_Reg) &
&/ Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          ! Temperature statistics
          Write (Output_Unit_HC_Data,55) "Average temp of particles reaching cutoff time","eV",( HC_Tot_Temp_Reach_Max_Iter (&
&HC_Species,Launch_Reg) / Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,*) ""
 
       End If
 
54     Format (1X,A50,2X,A4,X,100F11.8) ! Good for fractions up to 1.00000
! jdemod - add exponential notation for some of these numbers
541    Format (1X,A50,2X,A4,X,100(1x,E10.5)) ! Good for fractions up to 1.00000
! jdemod - switched 545 to 546 using exponential output 
545    Format (1X,A50,2X,A4,X,100F11.7) ! Good for fractions up to 1.00000
!546    Format (1X,A50,2X,A4,X,100E11.6) ! Good for fractions up to 1.00000
546    Format (1X,A50,2X,A4,X,100(1x,E10.5)) ! Good for fractions up to 1.00000
55     Format (1X,A50,2X,A4,X,100F11.4) ! Good for values up to 999.9999
56     Format (1X,A50,2X,A4,X,100F11.2) ! Good for values up to 99,999.99
57     Format (1X,A50,2X,A4,X,100F11.1) ! Good for values up to 999,999.9
58     Format (1X,A50,2X,A4,X,100E11.2) ! Good for values up to 0.999E99
 
    End Do
 
    ! Print evolution statistics.
    Write (Output_Unit_HC_Data,*) "3) Hydrocarbon EVOLUTION Statistics"
    Write (Output_Unit_HC_Data,*) ""
 
    ! Print table of hydrocarbon reactions from each target.
    Do Launch_Reg = 1,Number_Regions
 
       ! Check that molecules were launched from this region.
       If ( HC_Region_Used (Launch_Reg) .gt. 0) Then
 
          If (hc_disable_transitions .eq. 1) Then
             ! Assign small value to reaction count to avoid 'nan's.
             Reaction_Count_Sum =  Calc_Lo
          Else
             Reaction_Count_Sum = SUM ( HC_Reaction_Count (:,:,Launch_Reg))
             If (Reaction_Count_Sum .eq. 0.0) Then
                Reaction_Count_Sum =  Calc_Lo
             End If
          End If
          Write (Output_Unit_HC_Data,'(1X,A,1X,A10)') "Reaction Table for launches:", HC_Region_Names (Launch_Reg)
          Write (Output_Unit_HC_Data,'(2X,A6,9X,A10,7X,A8,5X,A9,4X,A6,3X,A13)') "Source","Product(s)","Reaction","Raw Count","Tot &
&Fr","Fr per Launch"
 
          Do Reaction = 1,Number_HC_Reactions
             ! Note:  There is more space allocated than reactions loaded.
             ! Check to be sure the particular reaction has been assigned.
             If (HC_State_Transform_Table (Reaction) % Start_State .gt. 0) Then
 
                ! Construct full list of reaction products for benefit of each higher hydrocarbon.
                Product_Count = 1
 
                Do While (Product_Count .le. Highest_Carbon_Content)
                   If (HC_State_Transform_Table (Reaction) % End_C_States (Product_Count) .gt. 0) Then
                      ! Note:  CxHy+, 5 characters, is the longest currently supported name.
                      Full_Product_List (1 + 6 * (Product_Count - 1):6 * Product_Count) = HC_Ident_Species (&
&HC_State_Transform_Table (Reaction) % End_C_States (Product_Count))
                   Else
                      Full_Product_List (1 + 6 * (Product_Count - 1):6 * Product_Count) = "      "
                   End If
 
                   ! Move to the next product, if applicable.
                   Product_Count = Product_Count + 1
                End Do
 
                Write (Output_Unit_HC_Data,30) HC_Ident_Species (HC_State_Transform_Table (Reaction) % Start_State),"to", &
                     & Full_Product_List,"by", &
                     & HC_State_Transform_Table (Reaction) % Reaction_Type,"impact", &
                     & SUM ( HC_Reaction_Count (:,Reaction,Launch_Reg)), &
                     & SUM ( HC_Reaction_Count (:,Reaction,Launch_Reg)) / Reaction_Count_Sum, &
                     & SUM ( HC_Reaction_Count (:,Reaction,Launch_Reg)) / Particle_Launch_Sum (Launch_Reg)
30              Format (2X,A7,1X,A2,1X,A18,A2,1X,A1,1X,A6,7X,F7.1,2X,F8.4,4X,F8.4)
             End If
          End Do
          Write (Output_Unit_HC_Data,'(A)') ""
       End If
    End Do
 
    ! Print hydrocarbon reaction table for all HC's to all HC products.
    Do Launch_Reg = 1,Number_Regions
 
       ! Check that molecules were launched from this region.
       If ( HC_Region_Used (Launch_Reg) .gt. 0) Then
 
          ! Initialize transition table.
          HC_Transition_Table = 0.0
 
          ! Go through HC reaction table and load into local transition table.
          Do HC_Species = 1, Num_HCs, 1
             Do Reaction = 1,Number_HC_Reactions,1
                !Write (Output_Unit_Scratch,*) "Reaction",HC_Species,Reaction,HC_State_Transform_Table (Reaction) % End_C_States
                Product_Count = 1
                Do While (Product_Count .le. Highest_Carbon_Content)
                   If (HC_State_Transform_Table (Reaction) % End_C_States (Product_Count) .ne. 0) Then
                      HC_Transition_Table (HC_Species,HC_State_Transform_Table (Reaction) % End_C_States (Product_Count)) = &
                           & HC_Transition_Table (HC_Species,HC_State_Transform_Table (Reaction) % End_C_States (Product_Count)) + &
                           &  HC_Reaction_Count (HC_Species,Reaction,Launch_Reg)
                   End If
                   Product_Count = Product_Count + 1
                End Do
             End Do
          End Do
 
          Write (Output_Unit_HC_Data,'(1X,A,1X,A10)') "Transition table for launches:", HC_Region_Names (Launch_Reg)
          Write(Output_Unit_HC_Data,'(2X,A6,3X,100A8)')"Source",(HC_State_Table(HC_Species)%State_Name,HC_Species=1,Num_HCs,1)
          Do HC_Species = 1, Num_HCs
             Write (Output_Unit_HC_Data,'(2X,A6,3X,100F8.1)') HC_State_Table (HC_Species) % State_Name, (HC_Transition_Table (&
&HC_Species,Target_Species),Target_Species = 1, HC_Species - 1, 1)
          End Do
          Write (Output_Unit_HC_Data,'(A)') ""	
       End If
    End Do
 
    ! Print data for breakup chains that extend all the way to C+ only.
    ! Initialize sums.
    Transition_Count_So_Far = 0.0
    Time_All_The_Way = 0.0
    Energy_All_The_Way = 0.0
    Kin_Energy_All_The_Way = 0.0
    Transition_Count_All_The_Way = 0.0
    Tot_Reach_State_All_The_Way =  Calc_Lo
    

    write(6,'(a,i6)') 'Output:',max_impurities

    Do Particle = 1, Max_Impurities
 
       ! Particle was launched.
       If ( HC_Energy_At_Production (Particle,HC_Species_Ident ("CH4")) .ne. 0.0) Then
          
          ! Add up number of steps required to get to each state.
          Transitions = 0.0

          Do HC_Species = Num_HCs - 1, 1, -1
             If ( HC_Time_At_Production (Particle,HC_Species) .ne. 0.0) Then
                ! Add one to total transitions made to get to current state.
                !
                ! jdemod - transition_count_so_far is intended to determine the average number of
                !          transitions needed to reach this state from the original state 
                !          so the total count of transitions to reach the state is summed up as one
                !          loops through the array. 
                !

                Transitions = Transitions + 1.0
                Transition_Count_So_Far (HC_Species) = Transition_Count_So_Far (HC_Species) + Transitions
             End If
          End Do
 
          ! Check to see if particle evolved to C+ and record data for those particles only.
          If ( HC_Time_At_Production (Particle,HC_Species_Ident ("C+")) .ne. 0.0) Then
 
             ! Found a particle that evolved down to C+.
             Transitions = 0.0
             Kin_Energy = 0.0
             Do HC_Species = Num_HCs - 1, 1, -1
                If ( HC_Time_At_Production (Particle,HC_Species) .ne. 0.0) Then

                   !write (6,'(a,i6,a,a,a,f10.5,a,g12.5,a,g12.5)') "Particle:",particle," made it to ",HC_Ident_Species(hc_species),&
                   !     &" at time ", HC_Time_At_Production (Particle,HC_Species),&
                   !     &" with energy ", HC_Energy_At_Production (Particle,HC_Species),&
                   !     &" with kin added ",hc_kin_e_add_at_production(particle,hc_species)

                   Tot_Reach_State_All_The_Way (HC_Species) = Tot_Reach_State_All_The_Way (HC_Species) + 1.0
                   Transitions = Transitions + 1.0
                   Transition_Count_All_The_Way (HC_Species) = Transition_Count_All_The_Way (HC_Species) + Transitions
                End If
 
                Time_All_The_Way (HC_Species) = Time_All_The_Way (HC_Species) +  HC_Time_At_Production (Particle,HC_Species)
                Energy_All_The_Way (HC_Species) = Energy_All_The_Way (HC_Species) +  HC_Energy_At_Production (Particle,HC_Species)
                Kin_Energy = Kin_Energy +  HC_Kin_E_Add_At_Production (Particle,HC_Species)
                Kin_Energy_All_The_Way (HC_Species) = Kin_Energy_All_The_Way (HC_Species) + Kin_Energy
 
             End Do
 
          End If
 
       End If
 
    End Do
 
    !
    ! jdemod - some of these are wrong in that they divide data from the state evolution tables by totals for all particles
    !          reaching those states - which ignores the fact that some of these go through wall interactions and thus are 
    !          not properly recorded 
    !

    Write (Output_Unit_HC_Data,'(1X,A,47X,A)') "Species specific evolution statistics for all launches:","Species Produced"
    Write (Output_Unit_HC_Data,'(42X,A9,2X,A4,2X,100A11)') "Statistic","Unit",(HC_State_Table (HC_Species) % State_Name,HC_Species &
                        &= 1,Num_HCs,1)
    Write (Output_Unit_HC_Data,62) "Average number of transitions to get to species","#",(Transition_Count_So_Far (HC_Species) / &
                 &SUM (Tot_Reach_State_Per_Launch (HC_Species,:)), HC_Species = 1, Num_HCs, 1)
    Write (Output_Unit_HC_Data,62) "Average transitions only for particles reach C+","#",(Transition_Count_All_The_Way (HC_Species)&
                 & / Tot_Reach_State_All_The_Way (HC_Species), HC_Species = 1, Num_HCs, 1)
    Write (Output_Unit_HC_Data,60) "Average time at prod only for particles reach C+","s",(Time_All_The_Way (HC_Species) / &
                 &Tot_Reach_State_All_The_Way (HC_Species), HC_Species = 1, Num_HCs, 1)
    Write (Output_Unit_HC_Data,61) "Average energy at prod only for particles reach C+","eV",(Energy_All_The_Way (HC_Species)  / &
                 &Tot_Reach_State_All_The_Way (HC_Species), HC_Species = 1, Num_HCs, 1)
    Write (Output_Unit_HC_Data,62) "Average kin energy add at prod only for reach C+","eV",(Kin_Energy_All_The_Way (HC_Species)  / &
                 &Tot_Reach_State_All_The_Way (HC_Species), HC_Species = 1, Num_HCs, 1)
    Write (Output_Unit_HC_Data,*) ""
 
    ! Print arrayed evolution data.
    Do Launch_Reg = 1,Number_Regions
 
       ! Check that molecules were launched from this region.
       If ( HC_Region_Used (Launch_Reg) .gt. 0) Then		
 
          Write (Output_Unit_HC_Data,'(1X,A,1X,A10,29X,A)') "Species and region specific evolution statistics for launches:", &
&HC_Region_Names (Launch_Reg),"Species Produced"
          Write (Output_Unit_HC_Data,'(42X,A9,2X,A4,2X,100A11)') "Statistic","Unit",(HC_State_Table (HC_Species) % State_Name,&
&HC_Species = 1,Num_HCs,1)
 
          Write (Output_Unit_HC_Data,61) "Average R position at production of species","m",( HC_R_Pos_At_Prod (HC_Species,&
&Launch_Reg) / Tot_Reach_State_Per_Launch (HC_Species,Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,61) "Average Z position at production of species","m",( HC_Z_Pos_At_Prod (HC_Species,&
&Launch_Reg) / Tot_Reach_State_Per_Launch (HC_Species,Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,61) "Average S dist from start at production of species","m",( HC_S_Pos_At_Prod (HC_Species,&
&Launch_Reg) / Tot_Reach_State_Per_Launch (HC_Species,Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,61) "Average cross position at production of species","m",( HC_Cross_Pos_At_Prod (HC_Species,&
&Launch_Reg) / Tot_Reach_State_Per_Launch (HC_Species,Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,62) "Average angle at production of species (+x=0,+ccw)","deg",( HC_Angle_At_Prod (HC_Species,&
&Launch_Reg) *raddeg / Tot_Reach_State_Per_Launch (HC_Species,Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,63) "Average number of occurances of species per launch","#",( HC_Num_Fragments (HC_Species,&
&Launch_Reg) / Tot_Reach_State_Per_Launch (HC_Species,Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,63) "Average ABS(velocity) at production of species","m/s",( HC_Velocity_At_Prod (HC_Species,&
&Launch_Reg) / Tot_Reach_State_Per_Launch (HC_Species,Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,63) "Maximum ABS(velocity) at production of species","m/s",( HC_Max_Velocity_At_Prod (&
&HC_Species,Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,62) "Average temperature at production of species","eV",( HC_Temperature_At_Prod (HC_Species,&
&Launch_Reg) / Tot_Reach_State_Per_Launch (HC_Species,Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,60) "Average cumulative time at production of species","s",( HC_Time_To_Prod (HC_Species,&
&Launch_Reg) / Tot_Reach_State_Per_Launch (HC_Species,Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,*) ""
 
       End If
 
60     Format (1X,A50,2X,A4,2X,100F11.8) ! Good for fractions/times up to 1.00000000
605    Format (1X,A50,2X,A4,2X,100F11.7) ! Good for fractions/times up to 10.0000000
61     Format (1X,A50,2X,A4,2X,100F11.6) ! Good for locations up to 9999.999999
62     Format (1X,A50,2X,A4,2X,100F11.3) ! Good for counts up to 99,999.99
63     Format (1X,A50,2X,A4,2X,100F11.1) ! Good for counts up to 999,999.9
 
    End Do
 
    ! Print re-release statistics.
    Write (Output_Unit_HC_Data,*) "4) Hydrocarbon RE-RELEASE Statistics"
    Write (Output_Unit_HC_Data,*) ""
 
    Do Launch_Reg = 1,Number_Regions
 
       ! Check that molecules were launched from this region.
       If ( HC_Region_Used (Launch_Reg) .gt. 0) Then		
 
          Write (Output_Unit_HC_Data,'(1X,A,1X,A10,35X,A)') "Species and region specific re-release statistics for launches:", &
&HC_Region_Names (Launch_Reg),"Species"
          Write (Output_Unit_HC_Data,'(42X,A9,2X,A4,4X,100A11)') "Statistic","Unit",(HC_State_Table (HC_Species) % State_Name,&
&HC_Species = 1,Num_HCs,1)
 
          Write (Output_Unit_HC_Data,73) "Number of particles striking vessel","#",( HC_Num_Striking_Target (HC_Species,Launch_Reg)&
& +  HC_Num_Reach_Wall (HC_Species,Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,74) "Number of reflections per striking species","#",( Total_HC_Reflections (HC_Species,&
&Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,73) "Number of reflected species","#",(SUM( HC_Sum_Fragments_ReLaunched (HC_Species,:,&
&Launch_Reg)), HC_Species = 1, Num_HCs, 1)				
          Write (Output_Unit_HC_Data,71) "Average R position of re-released species","m",( HC_Total_R_ReProd_Positions (HC_Species,&
&Launch_Reg) / &
               & SUM ( HC_Sum_Fragments_ReLaunched (HC_Species,:,Launch_Reg)), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,71) "Average Z position of re-released species","m",( HC_Total_Z_ReProd_Positions (HC_Species,&
&Launch_Reg) / &
               & SUM ( HC_Sum_Fragments_ReLaunched (HC_Species,:,Launch_Reg)), HC_Species = 1, Num_HCs, 1)
          Write(Output_Unit_HC_Data,71)"Average angle of re-released species","deg",&
               &(HC_Total_ReProd_Angles(HC_Species,Launch_Reg)*raddeg /&
               & SUM ( HC_Sum_Fragments_ReLaunched (HC_Species,:,Launch_Reg)), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,72) "Average velocity of re-released species","m/s",( HC_Total_ReProd_Velocities (HC_Species,&
&Launch_Reg) / &
               & SUM ( HC_Sum_Fragments_ReLaunched (HC_Species,:,Launch_Reg)), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,72) "Maximum velocity of re-released species","m/s",( HC_Max_ReProd_Velocities (HC_Species,&
&Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,71) "Average temperature of re-released species","eV",( HC_Tot_ReProd_Temperatures (&
&HC_Species,Launch_Reg) / &
               & SUM ( HC_Sum_Fragments_ReLaunched (HC_Species,:,Launch_Reg)), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,73) "Number of relaunches outside wall but moved inside","#",( HC_ReLau_Grid_Error_Moved_Okay &
&(HC_Species,Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,73) "Number of relaunches outside wall not moved inside","#",( HC_ReLau_Grid_Error_Moved_Out (&
&HC_Species,Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,71) "Average number of reflections per region launches","#",( Total_HC_Reflections (&
&HC_Species,Launch_Reg) / &
               & Particle_Launch_Sum (Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,74) "Maximum re-release events for a single launch","#",( Max_HC_Reflections_Found (&
&HC_Species,Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,73) "Number of particles lost to reflection limit","#",( HC_Reflection_Loss (HC_Species,&
&Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,*) ""
 
       End If
 
70     Format (1X,A50,2X,A4,2X,100F11.5) ! Good for fractions/times up to 1.00000
71     Format (1X,A50,2X,A4,2X,100F11.4) ! Good for locations up to 999.9999
72     Format (1X,A50,2X,A4,2X,100F11.2) ! Good for counts up to 99,999.99
73     Format (1X,A50,2X,A4,2X,100F11.1) ! Good for counts up to 999,999.9
74     Format (1X,A50,2X,A4,2X,100I11) ! Good for counts up to 99,999,999,999
 
    End Do
 
    ! Print vessel interaction statistics.
    Write (Output_Unit_HC_Data,*) "5) Hydrocarbon VESSEL INTERACTION Statistics"
    Write (Output_Unit_HC_Data,*) ""
 
    Do Launch_Reg = 1,Number_Regions
 
       ! Check that molecules were launched from this region.
       If ( HC_Region_Used (Launch_Reg) .gt. 0) Then		
 
          Write (Output_Unit_HC_Data,'(1X,A,1X,A10,18X,A)') "Species and region specific vessel interaction statistics for &
&launches:", HC_Region_Names (Launch_Reg),"Species Interacting"
          Write (Output_Unit_HC_Data,'(42X,A9,2X,A4,2X,100A11)') "Statistic","Unit",(HC_State_Table (HC_Species) % State_Name,&
&HC_Species = 1,Num_HCs,1)
 
          Write (Output_Unit_HC_Data,73) "Number of ions interacting with far periphery","#",( HC_RWall (HC_Species,Launch_Reg), &
&HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,73) "Number of species interacting with walls+targets","#",( HC_RDep (HC_Species,Launch_Reg), &
&HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,73) "Number of MTC fragments interacting with walls","#",( HC_MTC_Reach_Wall (HC_Species,&
&Launch_Reg), HC_Species = 1, Num_HCs, 1)			
          Write (Output_Unit_HC_Data,*) ""
       End If
 
    End Do
 
    ! Print end-of-life statistics.
    Write (Output_Unit_HC_Data,*) "6) Hydrocarbon END-OF-LIFE Statistics"
    Write (Output_Unit_HC_Data,'(A)') ""
 
    ! Print table of particle fates for launches from each target.
    Do Launch_Reg = 1,Number_Regions
 
       ! Check that molecules were launched from this region.
       If ( HC_Region_Used (Launch_Reg) .gt. 0) Then		
 
          ! Print headings.
          Write (Output_Unit_HC_Data,'(1XA,1X,A10,37X,A)') "Particle Fate Table for launches:", HC_Region_Names (Launch_Reg),&
&"Fraction of Species Launched"
          Write (Output_Unit_HC_Data,'(1X,A4,32X,A9,1X,A7,4X,100A8)') "Fate","Raw Count","Tot Fr",(HC_State_Table (HC_Species) % &
&State_Name,HC_Species = 1,Num_HCs,1)
 
          ! Current IFATEs to print: 1,10,11,12,13,14,19,20,21,22,23,24,25,26,29
          If (hc_wbc_comp_option .gt. 0) Then
             Fates = (/ 1,10,11,12,13,14,15,19,20,21,22,23,24,25,26,27,29,0,0,0 /)
          Else
             Fates = (/ 1,10,11,12,13,14,19,20,21,22,23,24,25,26,29,0,0,0,0,0 /)
          End If
 
          Do Fate_Num = 1,Size(Fates)
             If (Fates (Fate_Num) .ne. 0) Then
                Write(Output_Unit_HC_Data,20)Particle_Fate(Fates(Fate_Num)),SUM(HC_IFate_Count(:,Launch_Reg,Fates(Fate_Num))),&
                     & SUM ( HC_IFate_Count (:,Launch_Reg,Fates (Fate_Num))) / Particle_Launch_Sum (Launch_Reg), &
                     & (REAL ( HC_IFate_Count (HC_Species,Launch_Reg,Fates (Fate_Num))) / &
                     & (Particle_Launch_Sum (Launch_Reg)),HC_Species = 1,Num_HCs,1)
20              Format (2X,A,T41,I6,2X,F6.4,2X,100F8.4)
             End If
          End Do
          Write (Output_Unit_HC_Data,'(A)') ""
       End If
    End Do
 
    Do Launch_Reg = 1,Number_Regions
 
       ! Check that molecules were launched from this region.
       If ( HC_Region_Used (Launch_Reg) .gt. 0) Then		
 
          Write (Output_Unit_HC_Data,'(1X,A,1X,A10,24X,A)') "Species and region specific end-of-life statistics for launches:", &
&HC_Region_Names (Launch_Reg),"Species At End-Of-Life"
          Write (Output_Unit_HC_Data,'(42X,A9,2X,A4,2X,100A11)') "Statistic","Unit",(HC_State_Table (HC_Species) % State_Name,&
&HC_Species = 1,Num_HCs,1)
 
          Write (Output_Unit_HC_Data,56) "Fraction of total launch absorbed on actual targs","#",( HC_Num_Absorbed_Act_Target (&
&HC_Species,Launch_Reg) / &
               & (Particle_Launch_Sum (Launch_Reg) +  Calc_Lo), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,56) "Fraction of total launch absorbed on walls+targs","#",( HC_Num_Absorbed_TargWall (&
&HC_Species,Launch_Reg) / &
               & (Particle_Launch_Sum (Launch_Reg) +  Calc_Lo), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,56) "Min ion-equivalent TSs to first fragment absorbed","#",( HC_Min_Teq_To_Absorption (&
&HC_Species,Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,56) "Max ion-equivalent TSs to last fragment absorbed","#",( HC_Max_Teq_To_Absorption (&
&HC_Species,Launch_Reg), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,56) "Average ion-equivalent TSs to fragment absorption","#",( HC_Tot_Teq_To_Absorption (&
&HC_Species,Launch_Reg) / &
               & ( HC_Num_Absorbed_TargWall (HC_Species,Launch_Reg) +  Calc_Lo), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,56) "Average fragment temp at fragment absorption","eV",( HC_Tot_Temp_At_Absorption (&
&HC_Species,Launch_Reg) / &
               & ( HC_Num_Absorbed_TargWall (HC_Species,Launch_Reg) +  Calc_Lo), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,56) "Average velocity at fragment absorption","m/s",( HC_Tot_Velocity_At_Absorption (&
&HC_Species,Launch_Reg) / &
               & ( HC_Num_Absorbed_TargWall (HC_Species,Launch_Reg) +  Calc_Lo), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,56) "Average ABS(velocity) at fragment absorption","m/s",( HC_Tot_ABS_Vel_At_Absorption (&
&HC_Species,Launch_Reg) / &
               & ( HC_Num_Absorbed_TargWall (HC_Species,Launch_Reg) +  Calc_Lo), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,56) "Average electron temp at location of absorption","eV",( HC_Tot_Elec_Temp_At_Absorption (&
&HC_Species,Launch_Reg) / &
               & ( HC_Num_Absorbed_TargWall (HC_Species,Launch_Reg) +  Calc_Lo), HC_Species = 1, Num_HCs, 1)
          !				Write (Output_Unit_HC_Data,55) "Average temperature per charge state at absorption","eV",( HC_Ctexs (HC_Species,Launch_Reg) / &
          !				& ( HC_Num_Absorbed_TargWall (HC_Species,Launch_Reg) +  Calc_Lo), HC_Species = 1, Num_HCs, 1)
          Write (Output_Unit_HC_Data,'(A)') ""
       End If
    End Do
 
    ! Recorded Ion Velocity (RIV).
    ! Note:  This is only carried out for the initial launch
    ! hydrocarbon species, hc_sputtered_hc_species.
    If ( Debug_HC_V .and. Get_HC_Charge (HC_Sputtered_HC_Species) .ne. 0) Then
 
       Do Current_Ring = 1,  Num_Upper_Rings
          Do Current_Cell = 1, gNKS (Current_Ring)
             If ( HC_Density (Current_Cell,Current_Ring,hc_sputtered_hc_species) .gt. 0.0D0) Then
                ! Also calculate the mean ion velocity.
                !write (0,*) "data:",Current_Cell,Current_Ring, HC_SDVS (Current_Cell,Current_Ring), HC_Density (Current_Cell,Current_Ring,hc_sputtered_hc_species), &
                !&  HC_SDVS (Current_Cell,Current_Ring) /  HC_Density (Current_Cell,Current_Ring,hc_sputtered_hc_species) /  Ion_Time_Step
                HC_SDVS (Current_Cell,Current_Ring) =  HC_SDVS (Current_Cell,Current_Ring) / &
                     &  HC_Density (Current_Cell,Current_Ring,hc_sputtered_hc_species) /  Ion_Time_Step
                HC_SDVS2 (Current_Cell,Current_Ring) =  HC_SDVS2 (Current_Cell,Current_Ring) / &
                     &  HC_Density (Current_Cell,Current_Ring,hc_sputtered_hc_species) /  Ion_Time_Step**2.0
 
                HC_SDTImp (Current_Cell,Current_Ring) = ( HC_SDVS2 (Current_Cell,Current_Ring) - &
                     &HC_SDVS(Current_Cell,Current_Ring)**2.0)/9.58084e7*Find_HC_Mass(hc_sputtered_hc_species,H_Isotope_Composition)
 
                HC_SDVS3 (Current_Cell,Current_Ring,1) =  HC_SDVS3 (Current_Cell,Current_Ring,1) / &
                     &  HC_SDVS3 (Current_Cell,Current_Ring,2) /  Ion_Time_Step
                HC_SDVS3 (Current_Cell,Current_Ring,2) =  HC_SDVS3 (Current_Cell,Current_Ring,2) / &
                     &  HC_Density (Current_Cell,Current_Ring,hc_sputtered_hc_species) *  Ion_Time_Step
 
             End If
          End Do
       End Do
 
       ! Print out debug information on cells and their contents
       ! where there are particles with Vz > Vb (local)
       Write (Output_Unit_HC_Data,*) 'SUMMARY OF HC ION VELOCITIES EXCEEDING LOCAL SOUND SPEED:'
 
       Write (Output_Unit_HC_Data,*) ' IK   IR   HC VB=(2kT/m)^0.5   AVERAGE VZ      FRACTION       NUMBER'
 
       Do Current_Ring = 1,  Num_Upper_Rings
          Do Current_Cell = 1, gNKS (Current_Ring)
             If ( HC_SDVS3 (Current_Cell,Current_Ring,2) .gt. 0.0) Then
                Write (Output_Unit_HC_Data,'(3i5,1x,g13.5,1x,g13.5,1x,2g16.8)') Current_Cell,Current_Ring,hc_sputtered_hc_species, &
                     & 9.79E3 * SQRT (2.0 * gktibs (Current_Cell,Current_Ring) / Find_HC_Mass (hc_sputtered_hc_species,&
&H_Isotope_Composition)), &
                     &  HC_SDVS3 (Current_Cell,Current_Ring,1), HC_SDVS3 (Current_Cell,Current_Ring,2), &
                     &  HC_SDVS3 (Current_Cell,Current_Ring,2) *  HC_Density (Current_Cell,Current_Ring,hc_sputtered_hc_species)
             End If
          End Do
       End Do
 
       ! Calculate the distribution of velocities
       If ( Max_Velocity_Cells_Per_Ring .gt. gnks ( INJ_Ring_Number) .and. gnks ( INJ_Ring_Number) .gt. 0) Then
          Velocity_Cells = gnks ( INJ_Ring_Number)
       Else
          Velocity_Cells = 1
       End If
 
       Do Current_Cell = 1,Velocity_Cells
          Velocity_Total = 0.0
 
          Do Velocity_Bin = - hc_NVel,  hc_NVel + 1
             If ( HC_VelSpace (Velocity_Bin,Current_Cell) .ne. 0.0) Then
                HC_VelSpace (Velocity_Bin,Current_Cell) =  HC_VelSpace (Velocity_Bin,Current_Cell) /  HC_VelWeight (Velocity_Bin,&
&Current_Cell) /  Ion_Time_Step
             End If
             Velocity_Total = Velocity_Total +  HC_VelWeight (Velocity_Bin,Current_Cell)
          End Do
 
          Do Velocity_Bin = - hc_NVel,  hc_NVel + 1
             If (Velocity_Total .ne. 0.0) Then
                HC_VelWeight (Velocity_Bin,Current_Cell) =  HC_VelWeight (Velocity_Bin,Current_Cell) / Velocity_Total
             End If
          End Do
 
       End Do
 
       Do Current_Cell = 1,Velocity_Cells
 
          If (Velocity_Cells .eq. 1) Then
             ! jdemod - usage of velplate inconsistent with code in hc_ion_transport
             Vel =  hc_VelPlate
             Write (Output_Unit_HC_Data,*) 'Velocity Weight Distribution: ',  hc_VelPlate
          Else
             Vel=(9.79E3*SQRT(gktibs(Current_Cell,INJ_Ring_Number)/Find_HC_Mass(hc_sputtered_hc_species,H_Isotope_Composition)))
             Write (Output_Unit_HC_Data,*) 'Velocity Weight Distribution: Knot=',Current_Cell,9.79E3 * SQRT(gktibs(Current_Cell, &
&INJ_Ring_Number)/Find_HC_Mass (hc_sputtered_hc_species,H_Isotope_Composition)), Vel
          End If
 
          Do Velocity_Bin = - hc_NVel,  hc_NVel + 1
             If (Velocity_Bin .le. 0) Then
                HC_VelCoord (Velocity_Bin) = (FLOAT (Velocity_Bin) - 0.5) *  hc_VelSep
             ElseIf (Velocity_Bin .gt. 0) Then
                HC_VelCoord (Velocity_Bin) = (FLOAT (Velocity_Bin - 1) + 0.5) *  hc_VelSep
             End If
             Write (Output_Unit_HC_Data,'(i4,1x,f8.3,1x,f14.6,12(1x,g16.6))') Velocity_Bin,HC_VelCoord (Velocity_Bin),HC_VelCoord (&
&Velocity_Bin) * Vel, &
                  &  HC_VelWeight (Velocity_Bin,Current_Cell), HC_VelSpace (Velocity_Bin,Current_Cell)					
          End Do
 
       End Do
 
       ! Print out velocity array.
       Write (Output_Unit_HC_Data,'(//1X,''TABLE OF VELOCITY VALUES:'')')
       Do Current_Ring = 1,  Num_Upper_Rings
          Write (Output_Unit_HC_Data,9031) 'PRIMARY','TOTNEUT','IONIZ 1','IONIZ 2','IONIZ 3','IONIZ 4','IONIZ 5','IONIZ 6'
          Write (Output_Unit_HC_Data,9032)
          Do Current_Cell = 1, gNKS (Current_Ring)
             Write (Output_Unit_HC_Data,9033) Current_Cell,Current_Ring,gRS (Current_Cell,Current_Ring),gZS (Current_Cell,&
&Current_Ring),( HC_SDVS (Current_Cell,Current_Ring))
          End Do
       End Do
 
       ! Print out impurity ion temperature array.
       Write (Output_Unit_HC_Data,'(//1X,''TABLE OF ION TEMPERATURE VALUES:'')')
       Do Current_Ring = 1,  Num_Upper_Rings
          Write (Output_Unit_HC_Data,9031) 'PRIMARY','TOTNEUT','IONIZ 1','IONIZ 2','IONIZ 3','IONIZ 4','IONIZ 5','IONIZ 6'
          Write (Output_Unit_HC_Data,9032)
          Do Current_Cell = 1, gNKS (Current_Ring)
             Write (Output_Unit_HC_Data,9033) Current_Cell,Current_Ring,gRS (Current_Cell,Current_Ring),gZS (Current_Cell,&
&Current_Ring),(HC_SDTImp (Current_Cell,Current_Ring))
          End Do
       End Do
 
       ! End of debugv.
    ElseIf ( Debug_HC_V .and. Get_HC_Charge (hc_sputtered_hc_species) .eq. 0) Then
       Write (Output_Unit_HC_Alert,*) "Warning in HC_Output:  Ion velocity information is not useful for a neutral launch particle."
       Write (Output_Unit_HC_Alert,*) "Continuing to next data."
    End If

    !
    ! Print out diagnostic data from hc_diag_data - pass unit number for output as an argument
    !
    call print_hc_diag_data(output_unit_hc_data)

    !


    !
    ! Print sheath efield diagnostics
    !
    write(output_unit_hc_data,*) "Summary of Sheath Efield statistics:"

    nsheath_total = sum(sheath_zbin_data(:,1))
    
    if (nsheath_total.gt.0.0) then 
       write(output_unit_hc_data,'(a,1x,g14.6)') "Total number of particles in sheath       = ",&
                                                     & nsheath_total
       write(output_unit_hc_data,'(a,1x,g14.6)') "Average particle Z position in Sheath (m) = ",&
                                                     & average_zsheath/nsheath_total
       write(output_unit_hc_data,'(a,1x,g14.6)') "Average particle S position in Sheath (m) = ",&
                                                    & average_ssheath/nsheath_total
       write(output_unit_hc_data,'(a,1x,g14.6)') "Average Sheath Efield                     = ",&
                                                    & average_sheath_efield/nsheath_total
       write(output_unit_hc_data,'(a,1x,g14.6)') "Average DIVIMP Efield                     = ",&
                                                    & average_div_efield/nsheath_total
       write(output_unit_hc_data,'(a,1x,g14.6)') "Average Sheath Bfield                     = ",&
                                                    & average_sheath_bfield/nsheath_total
       write(output_unit_hc_data,'(a,1x,g14.6)') "Average Sheath Bangle                     = ",&
                                                    & average_sheath_bangle/nsheath_total

       write(output_unit_hc_data,'(a)') "Sheath Efield profile:"

       do in = 1,nsheath_bins
          if (sheath_zbin_data(in,1).gt.0.0) then 
             write(output_unit_hc_data,'(i6,2(1x,g14.6))') in, (sheath_extent * (real(in)-0.5)),&
                  & sheath_zbin_data(in,2)/sheath_zbin_data(in,1)
          else
             write(output_unit_hc_data,'(i6,2(1x,g14.6))') in, (sheath_extent * (real(in)-0.5)),&
                  & 0.0
          endif
       end do

    endif


 
9031 FORMAT(/1X,' IK IR    R	  Z  ',12(2X,A7))
9032 FORMAT(1X,131('-'))
9033 FORMAT(1X,2I3,2F7.3,1P,12E9.2)
9034 FORMAT(39X , 1P , 12E9.2 )
 
    write (Output_Unit_Scratch,*) "HC output processing complete."
 
  End Subroutine HC_Print_Output
 
End Module HC_Output
