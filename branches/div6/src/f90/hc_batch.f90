! -*-Mode:f90-*-
! HC_Batch.f90
! Wrapping routine called upon chemical sputtering event for
! group of sputtered particles.  Loops to call hc_follow for
! each particle, keeping track of how many have been followed.
! This routine may be called multiple times, for target, wall,
! free-space and distrubted launch events.
!
! Adam McLean under Prof. Peter Stangeby
! Research Assistant: David Elder
! March, 2003

Module HC_Batch

  ! Every good Fortran 90 program has...
  Implicit None	

Contains

  ! Hydrocarbon Update.  Routine to take care of multiple sputtered particles.
  Subroutine HC_Launch (	&
       &	LProd,		&	! First number of HC fragments to start with (in).
       &	NProd,		&	! Continue to this point (in).
       &	LatIZ,		&	! First position in ionization array to be filled (in).
       &	NatIZ,		&	! Number of initial hydrocarbons that become C+ ions (out).
       &	SStruk,		&	! Total of HC fragments striking target (out).
       &	MTCStruk,	&	! Momentum transfer collision (out).
       &	SMain,		&	! Number of fragments entering main plasma, IR < IRSEP (out).
       &	SExit,		&	! HC Fragments exiting main plasma (out).
       &	SatIZ,		&	! Total of HC fragments ionized to C+ (out).
       &	SHC,		&	! Sum of HC fragments launched, SNEUT (out).
       &	SWallHC,	&	! Sum of HC fragments reaching the wall (out).
       &	MTCWallHC,	&	! Momentum transfer collision (out).
       &	SCent,		&	! Sum of HC fragments reaching IR=1 (out).
       &	StMax,		&	! Number of HC fragments existing at T=Tmax (out).
       &	Seed,		&	! Random seed (in/out).
       &	NRand,		&	! Count of total random numbers used (in/out).
       &	Neutime,	&	! Time spent tracking neutrals (in/out).
       &	SFail,		&	! Number of failed launches V>Vmax at least 1000 (out).
       &	Status,		&	! Generation of launch =1, primary, =2, secondary, etc, 10=total (in).
       &	matp,	&	! MATP, Background bombarding ion type, used for yields (but only printed currently) (in).
       &	matt,&	! MATT, Target type, used for yields (in).
       &	NeutType,	&	! Wall (4,5) or target (1,2,3), physical (1,4) or chemical (2,5) sputter (in).
       &	Local_CNeutB,	&	! Launch location option, overrides H21.
       &	Local_CNeutC	)	! Launch angle/velocity flag, overrides H22.

    ! Additional Input variables from common blocks:
    ! -XPRODS: Array of R coordinates of launch position for each particle set in div.d6a (in) /Commons/CNeut.
    ! -YPRODS: Array of Z coordinates of launch position for each particle set in div.d6a (in) /Commons/CNeut.
    ! -FSRATE: Neutral timestep (in) /Commons/Comtor, ComHC.
    ! -QTIM: Ion timestep (in) /Commons/Comtor, ComHC.
    ! -RMAXS: Maximum random numbers allowed in VIN calcs (in) /Commons/CNeut.
    !
    ! Additional Output variables from common blocks:
    ! -XATIZS: Array with R coordinates of injected particles at point of C+ production (out) /Commons/CNeut.
    ! -YATIZS: Array with Z coordinates of injected particles at point of C+ production (out) /Commons/CNeut.
    ! -KATIZS: Array with K values (found in KKS, set by tau.d6a) at point of C+ production (out) /Commons/CNeut.
    ! -SATIZS: Array with K values (found in KKS, set by tau.d6a) at point of C+ production (out) /Commons/CNeut.
    ! -VINS: Velocity of particle at point of C+ production (out) /Commons/CNeut.
    !
    ! -HC_RSTRUK: Total hydrocarbon fragments striking target (out).
    ! -HC_RATIZ: Total injected hydrocarbons that become C+ ions (out).
    ! -HC_RNEUT: Total hydrocarbons launched (out).
    ! -HC_RWALLN: Total hydrocarbons plating out on walls (out).
    ! -HC_RCENT: Total hydrocarbons reaching central mirror, IR=1 (out).
    ! -HC_RTMAX: Total hydrocarbons existing at T=Tmax (out).
    ! -HC_RFAIL: Total hydrocarbons with launch failures (out).
    ! -HC_RMAIN: Total hydrocarbons reaching main plasma, IR < IRSEP (out).
    !
    ! Additional Input/Output variables from common blocks:
    ! -SPUTYS: Array with fragment sizes of all neutrals (in/out) /Commons/CNeut.

    ! Required modules.
    Use ComHC ! Contains regional data.
    Use HC_Init_DIV_Data ! Contains global initialization routines.
    Use HC_Init_DIV_Diag ! HC diagnostic data.
    Use HC_Utilities ! Contains data array initialization routine.
    Use HC_Put ! Access routines to DIVIMP data.
    Use HC_Follow ! Contains HC_Follow subroutine.

    ! Every good Fortran program has...		
    Implicit None

    ! Define input/output variables
    Integer, Intent (In) :: LPRod
    Integer, Intent (In) :: NProd
    Integer, Intent (In) :: LatIZ
    Integer, Intent (Out) :: NatIZ ! Number of C+ atoms produced from bulk launch.
    Real, Intent (Out) :: SSTruk
    Real, Intent (Out) :: MTCStruk
    Real, Intent (Out) :: SMain
    Real, Intent (Out) :: SExit
    Real, Intent (Out) :: SatIZ
    Real, Intent (Out) :: SHC
    Real, Intent (Out) :: SWallHC
    Real, Intent (Out) :: MTCWallHC
    Real, Intent (Out) :: SCent
    Real, Intent (Out) :: STMax
    Double Precision, Intent (InOut) :: Seed
    Integer, Intent (InOut) :: NRand
    Real, Intent (InOut) :: Neutime
    Real, Intent (Out) :: SFail
    Integer, Intent (In) :: Status
    Integer, Intent (In) :: matp
    Integer, Intent (In) :: matt
    Integer, Intent (InOut) :: NeutType
    Integer, Intent (In) :: Local_CNeutB
    Integer, Intent (In) :: Local_CNeutC

    ! Define local variables that will be assigned values.
    Integer :: Random_Numbers_Used ! KK, Random numbers to be made available.
    Integer :: Max_Random_Values ! KKLIM, Limit to create new random number set and begin counting at 0.0 again.
    Integer :: IProd
    Integer :: Loop_Counter
    Integer :: HC_Species

    ! jdemod - integer for local percentage calculation
    integer :: perc

    Real :: IonTime
    Real :: HC_Porm ! PORM
    Real :: Starting_Time ! STATIM
    Real :: Time_Used ! TIMUSD
    Real :: Part_Time ! PARTIM
    Real :: Divide_Region
    Real :: Divide_Region_1
    Real :: Divide_Region_2
    Real :: Divide_Region_3
    Real :: Launch_Sum
    Real :: SFP ! Number of far-periphery launches.
    Real, Dimension (1000000,2) :: Storage

    ! External functions.
    Real, External :: ZA02AS ! Timing routine in SysPGI.u6a.

    real :: SHC_init, SStruk_init,MTCStruk_init,SMain_init, SExit_init,SatIZ_init
    real :: SWallHC_init,MTCWallHC_init,SCent_init,STMax_init, SFail_init, SFP_init

    ! jdemod - for debugging - these are the values originally calculated
    real :: SHC_org, SStruk_org,MTCStruk_org,SMain_org, SExit_org,SatIZ_org
    real :: SWallHC_org,MTCWallHC_org,SCent_org,STMax_org, SFail_org, SFP_org

    real :: SHC_alt, SStruk_alt,MTCStruk_alt,SMain_alt, SExit_alt,SatIZ_alt
    real :: SWallHC_alt,MTCWallHC_alt,SCent_alt,STMax_alt, SFail_alt, SFP_alt


    write(0,*) 'Running HC_Launch'
    ! Initialize storage.
    Storage = 0.0

    ! Initialize output variables.
    NatIZ = 0

    SHC = 0.0
    SSTruk = 0.0
    MTCStruk = 0.0
    SMain = 0.0
    SExit = 0.0
    SatIZ = 0.0
    SWallHC = 0.0
    MTCWallHC = 0.0
    SCent = 0.0
    STMax = 0.0
    SFail = 0.0
    SFP = 0.0

    ! jdemod - Incorrect statements below ... each different kind of launch has storage allocated for accumulating data
    !          The code then sums only the specific required pieces when scaling at the end. 
    !
    !        - the HC code uses the same arrays to accumulate all HC data regardless of the sources (the same way DIVIMP
    !          records impurity information.
    !        - However, the HC code uses the array data to determine the totals of particles lost to various
    !          mechanisms. Unfortunately, this makes the code non-reentrant unless the initial values of these are calculated
    !          and subtracted at the end since the calling code expects values ONLY for the current generation of particles. 

    SHC_init = SUM(HC_Sum_Fragments_Launched(:,1:4)) + SUM(HC_Sum_Fragments_Launched(:,5:8))
    SStruk_init = SUM ( HC_Num_Striking_Target (:,:))
    MTCStruk_init = SUM ( HC_MTC_Striking_Target (:,:))
    SMain_init = SUM ( HC_Num_Enter_Main_Plasma (:,:))
    SExit_init = SUM ( HC_Tot_Fragments_Exit_Main (:,:))
    SatIZ_init = SUM(HC_Num_Fragments_Reach_CIon (:))
    SWallHC_init = SUM ( HC_Num_Reach_Wall (:,:))
    MTCWallHC_init = SUM ( HC_MTC_Reach_Wall (:,:))
    SCent_init = SUM ( HC_Num_Reach_Centre (:,:))
    STMax_init = SUM ( HC_Num_At_TMax (:,:))
    SFail_init = SUM ( HC_Num_Failed_Launches (:,:))
    SFP_init = SUM ( HC_RFPTarg (:,:))



    ! Assign values to global variables.  Note these are set here
    ! because they are local variables in DIVIMP.
    Plasma_Type = matp ! MATP
    Target_Material = matt ! MATT

    ! jdemod - assigning a number prevents use as free space launch - though maybe the issue
    !          is the choice of value assigned - need to see if neuttype is used as an array index
    ! If NEUTTYPE equals 0, must assign a value.
    !If  (NeutType .eq. 0) Then
    !	NeutType = 1
    write (0,*) "Making NEUTTYPE:",NeutType		
    !EndIf

    ! Set overrides in ComHC.  This must be done for each
    ! group of hydrocarbons launched.
    If (hc_launch_location .ne. 6) Then
       hc_launch_location = Local_CNeutB
    End If
    hc_launch_angle_velocity = Local_CNeutC

    Write (Output_Unit_Scratch,9013) '*** Following Hydrocarbon Group ***'
9013 FORMAT(/1X,A)

    ! Check to see if printout of debug data was specified in DIVIMP input file.
    !
    ! NOTE: The HC code should really have its own debug flags - it isn't a good idea
    !       to piggyback this on the DIVIMP debug flags. Add a separate implementation to
    !       the TO-DO list. (jde Jan/2007)
    !
    If ( Print_Debug_x_CStep_Ion .gt. 0.0) Then
       ! jdemod - should print ion time step for ion debugging data
       !Write (Output_Unit_Scratch,9004) NINT ( Print_Debug_x_CStep_Ion),  Neutral_Time_Step
       Write (Output_Unit_Scratch,9004) NINT ( Print_Debug_x_CStep_Ion),  Ion_Time_Step
9004   FORMAT('DIV DEBUG: ION  DIAGNOSTICS TO BE PRINTED EVERY',I6,' TIMESTEPS  (DELTA T =',1P,G10.3,' SECONDS).',//)
       Debug_HC_Ion = .true. ! DEBUGL, Override input setting.
       !       Debug_HC_Ion = .false. ! DEBUGL, Override input setting.
    End If

    If ( Print_Debug_x_CStep_Neutral .gt. 0.0) Then
       Write (Output_Unit_Scratch,9005) NINT ( Print_Debug_x_CStep_Neutral),  Neutral_Time_Step
9005   FORMAT('DIV DEBUG: NEUT DIAGNOSTICS TO BE PRINTED EVERY',I6,' TIMESTEPS  (DELTA T =',1P,G10.3,' SECONDS).',//)
       Debug_HC_Neutral = .true. ! DEBUGN, Override input setting.
       !       Debug_HC_Neutral = .false. ! DEBUGN, Override input setting.
    End If

    write (0,*) 'Debugging:',debug_hc_ion,debug_hc_neutral


    ! Begin timing for this group of particles.
    Starting_Time = ZA02AS (1)

    ! Create random values that will be used later.
    ! Note:  these are not loaded before HC is called for an ion injection case
    ! and will offset the number of calls to SURAND by 3 in comparison to a non-HC case.
    !Call peranva (Seed,NProd - LProd + 1)
    !Call peranvb (Seed,NProd - LProd + 1)
    !Call peranvc (Seed,NProd - LProd + 1)
    !NRand = NRand + 3 * (NProd - LProd + 1)
    Random_Numbers_Used = 1000 *  Num_Sectors ! KK = 1000 * ISECT
    Max_Random_Values = Random_Numbers_Used - 100 ! KKLIM = KK - 10
    HC_Porm = -1.0 ! PORM, used for direction of ion injection (INJECTION options only).

    ! Loop for each hydrocarbon fragment to be launched in this group.
    Do IProd = 1, NProd - LProd + 1, 1

       ! jdemod - changed to print progress only every 10%
       !if (iprod/100.0.eq.real(int(iprod/100))) write(0,*) 'Following Particle:',iprod,nprod,lprod
       if (mod(iprod,(nprod-lprod+1)/10).eq.0) then 
            perc = int((iprod*10)/((nprod-lprod+1)/10))
            write(0,'(a,i3,a,i8)')   'HC Following: ',perc,' % complete. Particle # =',iprod
       endif


       !Do IProd = 1, 1000000, 1
       if (debug_hc) then 
          Write (Output_Unit_Scratch,'(a,1x,i6,1x,a,1x,i6,1x,a,1x,i9)') "Particle #",IProd,"Total particles",NProd-LProd+1,"Starting &
            &step", HC_Walk_Count
       endif

       ! Particle transport (?)


       Call HC_Transport (	&
            &	IProd,		&
            &	LPRod,		&
            &	NProd,		&
            &	LatIZ,		&
            &	NatIZ,		&
            &	Seed,		&
            &	NRand,		&
            &	Neutime,	&
            &	IonTime,	&
            &	SFail,		&
            &	Status,		&
            &	NeutType,	&
            &	HC_Porm,        & ! PORM
            &	Random_Numbers_Used,	& ! KK
            &	Max_Random_Values, 	& ! KKLIM
            &	Storage)


       ! Note, at this point, IonTime is not returned to the call to HC_Launch in Neut.d6a.

       ! Use timing implementation from DIVIMP, originally used in LIM.
       ! Check for CPU time limit before moving to next particle.
       ! THE FOLLOWING TIME LIMIT CODE IS TAKEN DIRECTLY FROM LIM WHERE
       ! IT WORKS WITH REASONABLE SUCCESS. THE SAME IS EXPECTED HERE. IT
       ! IS NECESSARY TO PREVENT RUN-AWAY CONDITIONS FROM DEVELOPING. ON
       ! THE WORKSTATION LONG RUNS ARE ACCEPTABLE, SO A TIME LIMIT OF
       ! 36000.0 SECONDS MAY BE USED.
       ! SEE IF TEST OF CPU TIME USED IS DUE TRAP CASE WHERE TIMUSD=0 OCCURS.

       If (IProd .ge.  Impurity_Limit) Then ! Note, in DIV, the comparison is to IMP, the number of ionized components, not IPROD, the number of original neutrals being followed.
          Time_Used = ZA02AS(1) - Starting_Time
          If (Time_Used .gt. 0.0) Then
             Part_Time = ( CPU_Time_Limit - Neutime - Iontime) / Time_Used
          Else
             Part_Time = 10.0
          End If
          !Write (Output_Unit_Scratch,'a,1x,i6,1x,a,1x,f10.3,1x,a,1x,f10.5,1x,a,1x,f12.9,1x,a,1x,f12.9') "IProd",IProd,"Time Used",Time_Used,"Part Time",Part_Time,"Neut Time",Neutime,"Ion Time",Iontime
          !If (Part_Time .ge. 1.05) Then
          ! jdemod - what is the point of changing the value from 1.05 to -100000? 
          If (Part_Time .ge. -100000) Then
             Impurity_Limit = INT (REAL (IProd) * MIN (4.0, 0.25 + 0.75 * Part_Time))
          Else

             ! HAVE RUN OUT OF CPU TIME, STOP ITERATION
             ! THERE ARE SO MANY COUNTERS IT IS VIRTUALLY IMPOSSIBLE TO
             ! WIND UP THE ROUTINE CLEANLY.  JUST WORK OUT HOW MANY IONS
             ! HAVE BEEN FOLLOWED IN THE TIME ALLOTTED AND CALL IT A DAY.

             Impurity_Limit = IProd
             !					Do Loop_Counter = IProd+1,  HC_Num_Fragments_Reach_CIon
             !						 HC_Num_Fragments_Reach_CIon =  HC_Num_Fragments_Reach_CIon - gsputys (Loop_Counter)
             !					End Do
             ! NATIZ = IMP ! Not required for multi-generation following model use for hydrocarbon transport.
             Write (Output_Unit_HC_Alert,'('' Error:  CPU time limit reached during hydrocarbon following.'')')
             Write (Output_Unit_HC_Alert,'('' Number of hydrocarbons reduced to'',I5)') NINT ( HC_Num_Fragments_Reach_CIon)
             Call PRB
             Call PRC ('Error:  CPU time limit reached during hydrocarbon following.')
             Call PRI ('Number of hydrocarbons reduced to ', NINT ( HC_Num_Fragments_Reach_CIon))
             Call PRC ('Increase HC CPU time or increase quantum timestep.')
             Call PRC ('Results for hydrocarbons will be printed but they should be treated with caution...')
             Exit ! Exit from particle following DO loop.
          End If
       End If
    End Do

    ! Finished following particles
    ! Print diagnostic and summary information

    ! --------------------------------------------------------------------------------------------------------------------------------------------------------

    write(0,*) 'HC_LAUNCH Printing:',neuttype

    ! Print out statistical data recorded for this particle group.

    Call Bin_Test (storage (:,1),-2* Pi_Value,2* Pi_Value,500)
    Call Bin_Test (storage (:,2),0.0,5000.0,100)


    !
    !
    !     jdemod
    !
    !     The following are the definitions of neuttype from neut.f
    !
    !     The HC code is usually only called for chemical sputtering options though it could 
    !     potentially be called for a freespace HC launch as well. 
    !
    !     The HC code is invoked if it's option is turned on and either chemical sputtering is called for
    !     or a freespace launch (cneutb=1) is specified.
    !
    !     In addition - all sputtering and reflection by hydrocarbons is handled internally - this code
    !     does not create additional batchs or change the value of neuttype based on the specific neutral - thus
    !     the summarized numbers here MUST include the reflected or sputtered contributions as well. These summations are
    !     independent of neuttype. In fact neuttype = 3 is misinterpreted in the existing code - it applies to an entire batch
    !     of self-sputtered particles since this is how DIVIMP was implemented originally. 
    !
    !
    !     Free Space= 0 - most data not recorded - cneutb=1
    !     Targ+phys = 1
    !     Targ+chem = 2
    !     Targ-self = 3
    !     Wall+phys = 4
    !     Wall+chem = 5
    !     2Dneutral = 6
    !     Refl Ion  = 7
    !
    !     A problem may arise since some of the launch regions share the same index - on the other hand this 
    !     may be due to a misunderstanding of the DIVIMP launch options. For example,DIVIMP launch options 0 and 3 are
    !     identical - both target physically sputtered. Historically, option 0 did not have self-sputtering while option 3
    !     did - however as more launch options were added the self-sputtering was split to a separate input cselfs and 
    !     launch options 0 and 3 became the same. The situation is similar for launch options 2 and 4. 
    !
    !


    ! Sum values for both targets, all hydrocarbon species for this group launch.
    ! jdemod - these number leave out a lot of the sources from different regions
    !          
    ! Notes: 
    ! 1) All particle sources need to be added together including sputtered and free space launches  
    !    - reflected particles are not a particle source
    ! 2) All particle fates need to be summed including reflected particles and all particle sources
    ! 3) These summaries should be the same no matter which neut type is launched. 
    !
    ! HC_Launch_Reg_Target_1_Dist = 1	! DIVIMP launch option 0. NeutType 2.
    ! HC_Launch_Reg_Target_1_Pin = 1		! DIVIMP launch option 3. NeutType 2.
    ! HC_Launch_Reg_Target_2_Dist = 2	! DIVIMP launch option 0. NeutType 2.
    ! HC_Launch_Reg_Target_2_Pin = 2		! DIVIMP launch option 3. NeutType 2.
    ! HC_Launch_Reg_Wall_Homo = 3		! DIVIMP launch option 2. NeutType 5.
    ! HC_Launch_Reg_Wall_Dist = 3		! DIVIMP launch option 4. NeutType 5.
    !
    ! HC_Launch_Reg_Free_Space_PT = 4	! DIVIMP launch option 1. NeutType 6.
    ! HC_Launch_Reg_Free_Space_2D = 4	! DIVIMP launch option 5. NeutType 6.
    !
    ! HC_Launch_Reg_Sput_Target_1 = 5	! Data specific for sputtered particles. NeutType 3.
    ! HC_Launch_Reg_Sput_Target_2 = 6	! Data specific for sputtered particles. NeutType 3.
    ! HC_Launch_Reg_Sput_Wall = 7		! Data specific for sputtered particles. NeutType 3.
    ! HC_Launch_Reg_Sput_FS_PT = 8		! Data specific for sputtered particles. NeutType 3.
    ! HC_Launch_Reg_Sput_FS_2D = 8		! Data specific for sputtered particles. NeutType 3.
    ! 
    ! HC_Launch_Reg_Refl_Target_1 = 9	! Data specific for reflected particles. NeutType 7.
    ! HC_Launch_Reg_Refl_Target_2 = 10	! Data specific for reflected particles. NeutType 7.
    ! HC_Launch_Reg_Refl_Wall = 11		! Data specific for reflected particles. NeutType 7.
    !
    ! HC_Launch_Reg_Refl_FS_PT = 12		! Data specific for reflected particles. NeutType 7.
    ! HC_Launch_Reg_Refl_FS_2D = 12		! Data specific for reflected particles. NeutType 7.


    !
    ! All of these arrays contain data for each HC state - 10 in all from CH4 to C+
    ! Thus - summming over all states may include some particle more than once in the totals 


    ! Total number of fragments launched is the total from launch_regions 1 to 4
    ! plus sputtered particles 
    SHC_alt = SUM(HC_Sum_Fragments_Launched(:,1:4)) + SUM(HC_Sum_Fragments_Launched(:,5:8))

    ! These variables are supposed to record final fate statistics - meaning what happens to the
    ! particle at the end of its lifetime. This appears to be the case even though the launch_region
    ! is changed while following the particle. 

    ! Number reaching the wall
    SWallHC_alt = SUM ( HC_Num_Reach_Wall (:,:))

    ! Number reaching the wall after a Momentum Transfer Collision
    MTCWallHC_alt = SUM ( HC_MTC_Reach_Wall (:,:))

    ! Number reaching the centre of the plasma
    SCent_alt = SUM ( HC_Num_Reach_Centre (:,:))

    ! Number at Tmax
    STMax_alt = SUM ( HC_Num_At_TMax (:,:))

    ! Number striking target 
    SStruk_alt = SUM ( HC_Num_Striking_Target (:,:))

    ! Number striking target after an MTC 
    MTCStruk_alt = SUM ( HC_MTC_Striking_Target (:,:))

    ! Number of failed launches
    SFail_alt = SUM ( HC_Num_Failed_Launches (:,:))

    ! Number of fragments reaching C+
    SatIZ_alt = SUM(HC_Num_Fragments_Reach_CIon (:))


    !write(0,'(a,10(1x,g12.5))') 'SATIZ:',HC_Num_Fragments_Reach_CIon (HC_Launch_Reg_Target_1_Dist),&
    !     &	 HC_Num_Fragments_Reach_CIon (HC_Launch_Reg_Target_1_Pin), &
    !     &	 HC_Num_Fragments_Reach_CIon (HC_Launch_Reg_Target_2_Dist), &
    !     &	 HC_Num_Fragments_Reach_CIon (HC_Launch_Reg_Target_2_Pin),Divide_Region,satiz

    ! Number entered main plasma
    SMain_alt = SUM ( HC_Num_Enter_Main_Plasma (:,:))

    ! Number exiting main plasma
    SExit_alt = SUM ( HC_Tot_Fragments_Exit_Main (:,:))

    ! Number reaching Far Periphery target
    SFP_alt = SUM ( HC_RFPTarg (:,:))




    ! jdemod - replacing the code below with one block caused errors. Some launch_reg 
    !          are used for recording reflected and other data that should not be included
    !          with original sources or losses from original sources. 
    !
    ! Note: NeutType 1 and 4 are physically sputtered and should never call HC_Launch.
    If (NeutType .eq. 2) Then
       ! Chemically sputtered from target.
       If (HC_Launch_Reg_Target_1_Dist .eq. HC_Launch_Reg_Target_1_Pin .and. &
            &   HC_Launch_Reg_Target_2_Dist .eq. HC_Launch_Reg_Target_2_Pin .and. &
            &   HC_Launch_Reg_Target_1_Dist .eq. HC_Launch_Reg_Target_2_Dist) Then
          Divide_Region = 4.0
       ElseIf (HC_Launch_Reg_Target_1_Dist .eq. HC_Launch_Reg_Target_1_Pin .or. &
            &   HC_Launch_Reg_Target_2_Dist .eq. HC_Launch_Reg_Target_2_Pin) Then
          Divide_Region = 2.0
       Else
          Divide_Region = 1.0
       End If
    
       !write(0,*) 'HC_LAUNCH Printing:divide_region:',divide_region,HC_Launch_Reg_Target_1_Dist,&
       !     &HC_Launch_Reg_Target_1_Pin,HC_Launch_Reg_Target_2_Dist,HC_Launch_Reg_Target_2_Pin
        
          !jdemod - original code
    
           SHC_org = (SUM ( HC_Sum_Fragments_Launched (:,HC_Launch_Reg_Target_1_Dist))+ &
                &     SUM ( HC_Sum_Fragments_Launched (:,HC_Launch_Reg_Target_1_Pin))+ &
                &     SUM ( HC_Sum_Fragments_Launched (:,HC_Launch_Reg_Target_2_Dist))+ &
                &     SUM ( HC_Sum_Fragments_Launched (:,HC_Launch_Reg_Target_2_Pin)))/Divide_Region
    
           SWallHC_org = (SUM ( HC_Num_Reach_Wall (:,HC_Launch_Reg_Target_1_Dist))+ &
                &	  SUM ( HC_Num_Reach_Wall (:,HC_Launch_Reg_Target_1_Pin))+ &
                &	  SUM ( HC_Num_Reach_Wall (:,HC_Launch_Reg_Target_2_Dist))+ &
                &	  SUM ( HC_Num_Reach_Wall (:,HC_Launch_Reg_Target_2_Pin)))/Divide_Region
    
           MTCWallHC_org = (SUM ( HC_MTC_Reach_Wall (:,HC_Launch_Reg_Target_1_Dist))+ &
                &	    SUM ( HC_MTC_Reach_Wall (:,HC_Launch_Reg_Target_1_Pin))+ &
                &	    SUM ( HC_MTC_Reach_Wall (:,HC_Launch_Reg_Target_2_Dist))+ &
                &	    SUM ( HC_MTC_Reach_Wall (:,HC_Launch_Reg_Target_2_Pin)))/Divide_Region
    
           SCent_org = (SUM ( HC_Num_Reach_Centre (:,HC_Launch_Reg_Target_1_Dist))+ &
                &	SUM ( HC_Num_Reach_Centre (:,HC_Launch_Reg_Target_1_Pin))+ &
                &	SUM ( HC_Num_Reach_Centre (:,HC_Launch_Reg_Target_2_Dist))+ &
                &	SUM ( HC_Num_Reach_Centre (:,HC_Launch_Reg_Target_2_Pin)))/Divide_Region
    
           STMax_org = (SUM ( HC_Num_At_TMax (:,HC_Launch_Reg_Target_1_Dist))+ &
                &	SUM ( HC_Num_At_TMax (:,HC_Launch_Reg_Target_1_Pin))+ &
                &	SUM ( HC_Num_At_TMax (:,HC_Launch_Reg_Target_1_Dist))+ &
                &	SUM ( HC_Num_At_TMax (:,HC_Launch_Reg_Target_1_Pin)))/Divide_Region
    
           SStruk_org = (SUM ( HC_Num_Striking_Target (:,HC_Launch_Reg_Target_1_Dist))+ &           
                &	 SUM ( HC_Num_Striking_Target (:,HC_Launch_Reg_Target_1_Pin))+ &            
                &	 SUM ( HC_Num_Striking_Target (:,HC_Launch_Reg_Target_2_Dist))+ &           
                &	 SUM ( HC_Num_Striking_Target (:,HC_Launch_Reg_Target_2_Pin)))/Divide_Region

           write(output_unit_scratch,'(a,4i8,20(1x,g18.8))') 'SSTRUK_org:',&
                & HC_Launch_Reg_Target_1_Dist,HC_Launch_Reg_Target_1_Pin,HC_Launch_Reg_Target_2_Dist,HC_Launch_Reg_Target_2_Pin, divide_region, &
                &  SUM ( HC_Num_Striking_Target (:,HC_Launch_Reg_Target_1_Dist)),&
                &  SUM ( HC_Num_Striking_Target (:,HC_Launch_Reg_Target_1_Pin )),&
                &  SUM ( HC_Num_Striking_Target (:,HC_Launch_Reg_Target_2_Dist)),&
                &  SUM ( HC_Num_Striking_Target (:,HC_Launch_Reg_Target_2_Pin)),&
                &  (SUM ( HC_Num_Striking_Target (:,HC_Launch_Reg_Target_1_Dist))+ &           
                &    SUM ( HC_Num_Striking_Target (:,HC_Launch_Reg_Target_1_Pin))+ &            
                &    SUM ( HC_Num_Striking_Target (:,HC_Launch_Reg_Target_2_Dist))+ &           
                &    SUM ( HC_Num_Striking_Target (:,HC_Launch_Reg_Target_2_Pin)))/Divide_Region,&
                &  SUM ( HC_Num_Striking_Target (:,:))


    
           MTCStruk_org = (SUM ( HC_MTC_Striking_Target (:,HC_Launch_Reg_Target_1_Dist))+ &
                &	   SUM ( HC_MTC_Striking_Target (:,HC_Launch_Reg_Target_1_Pin))+ &
                &	   SUM ( HC_MTC_Striking_Target (:,HC_Launch_Reg_Target_2_Dist))+ &
                &	   SUM ( HC_MTC_Striking_Target (:,HC_Launch_Reg_Target_2_Pin)))/Divide_Region
    
           SFail_org = (SUM ( HC_Num_Failed_Launches (:,HC_Launch_Reg_Target_1_Dist))+ &
                &	SUM ( HC_Num_Failed_Launches (:,HC_Launch_Reg_Target_1_Pin))+ &
                &	SUM ( HC_Num_Failed_Launches (:,HC_Launch_Reg_Target_2_Dist))+ &
                &	SUM ( HC_Num_Failed_Launches (:,HC_Launch_Reg_Target_2_Pin)))/Divide_Region
    
           SatIZ_org = ( HC_Num_Fragments_Reach_CIon (HC_Launch_Reg_Target_1_Dist)+ &
                &	 HC_Num_Fragments_Reach_CIon (HC_Launch_Reg_Target_1_Pin)+ &
                &	 HC_Num_Fragments_Reach_CIon (HC_Launch_Reg_Target_2_Dist)+ &
                &	 HC_Num_Fragments_Reach_CIon (HC_Launch_Reg_Target_2_Pin))/Divide_Region
    
           !write(0,'(a,10(1x,g12.5))') 'SATIZ:',HC_Num_Fragments_Reach_CIon (HC_Launch_Reg_Target_1_Dist),&
           !     &	 HC_Num_Fragments_Reach_CIon (HC_Launch_Reg_Target_1_Pin), &
           !     &	 HC_Num_Fragments_Reach_CIon (HC_Launch_Reg_Target_2_Dist), &
           !     &	 HC_Num_Fragments_Reach_CIon (HC_Launch_Reg_Target_2_Pin),Divide_Region,satiz
                
           SMain_org = (SUM ( HC_Num_Enter_Main_Plasma (:,HC_Launch_Reg_Target_1_Dist))+ &
                &	SUM ( HC_Num_Enter_Main_Plasma (:,HC_Launch_Reg_Target_1_Pin))+ &
                &	SUM ( HC_Num_Enter_Main_Plasma (:,HC_Launch_Reg_Target_2_Dist))+ &
                &	SUM ( HC_Num_Enter_Main_Plasma (:,HC_Launch_Reg_Target_2_Pin)))/Divide_Region
           SExit_org = (SUM ( HC_Tot_Fragments_Exit_Main (:,HC_Launch_Reg_Target_1_Dist))+ &
                &	SUM ( HC_Tot_Fragments_Exit_Main (:,HC_Launch_Reg_Target_1_Pin))+ &
                &	SUM ( HC_Tot_Fragments_Exit_Main (:,HC_Launch_Reg_Target_2_Dist))+ &
                &	SUM ( HC_Tot_Fragments_Exit_Main (:,HC_Launch_Reg_Target_2_Pin)))/Divide_Region
           SFP_org = (SUM ( HC_RFPTarg (:,HC_Launch_Reg_Target_1_Dist))+ &
                &      SUM ( HC_RFPTarg (:,HC_Launch_Reg_Target_1_Pin))+ &
                &      SUM ( HC_RFPTarg (:,HC_Launch_Reg_Target_2_Dist))+ &
                &      SUM ( HC_RFPTarg (:,HC_Launch_Reg_Target_2_Pin)))/Divide_Region
        ElseIf (NeutType .eq. 3) Then
           ! NeutType indicates target launched and later self-sputtered.  For DIVIMP-HC,
           ! this would be all sputtered particles including those originally from the wall or in FS.
     
           ! Divide_Region_1 = Target chemical launches
           ! Divide_Region_2 = Wall chemical launches
           ! Divide_Region_3 = FS chemical launches
           Divide_Region = 1.0
           Divide_Region_1 = 1.0
           Divide_Region_2 = 1.0
           Divide_Region_3 = 1.0
     
           If (HC_Launch_Reg_Sput_Target_1 .eq. HC_Launch_Reg_Sput_Target_2 .and. &
                &   HC_Launch_Reg_Sput_Target_1 .eq. HC_Launch_Reg_Sput_Wall .and. &
                &   HC_Launch_Reg_Sput_Target_1 .eq. HC_Launch_Reg_Sput_FS_Pt .and. &
                &   HC_Launch_Reg_Sput_Target_1 .eq. HC_Launch_Reg_Sput_FS_2D) Then
              Divide_Region = 5.0
           ElseIf (HC_Launch_Reg_Sput_Target_1 .eq. HC_Launch_Reg_Sput_Target_2) Then
              Divide_Region_1 = 2.0
           ElseIf (HC_Launch_Reg_Sput_FS_PT .eq. HC_Launch_Reg_Sput_FS_2D) Then
              Divide_Region_3 = 2.0
           End If
     
           SHC_org = ((SUM ( HC_Sum_Fragments_Launched (:,HC_Launch_Reg_Sput_Target_1))+ &
                &      SUM ( HC_Sum_Fragments_Launched (:,HC_Launch_Reg_Sput_Target_2)))/Divide_Region_1+ &
                &      SUM ( HC_Sum_Fragments_Launched (:,HC_Launch_Reg_Sput_Wall))/Divide_Region_2+ &
                &      (SUM ( HC_Sum_Fragments_Launched (:,HC_Launch_Reg_Sput_FS_Pt))+ &
                &      SUM ( HC_Sum_Fragments_Launched (:,HC_Launch_Reg_Sput_FS_2D)))/Divide_Region_3)/Divide_Region
           SWallHC_org = ((SUM ( HC_Num_Reach_Wall (:,HC_Launch_Reg_Sput_Target_1))+ &
                &	   SUM ( HC_Num_Reach_Wall (:,HC_Launch_Reg_Sput_Target_2)))/Divide_Region_1+ &
                &	   SUM ( HC_Num_Reach_Wall (:,HC_Launch_Reg_Sput_Wall))/Divide_Region_2+ &
                &	   (SUM ( HC_Num_Reach_Wall (:,HC_Launch_Reg_Sput_FS_Pt))+ &
                &	   SUM ( HC_Num_Reach_Wall (:,HC_Launch_Reg_Sput_FS_2D)))/Divide_Region_3)/Divide_Region
           MTCWallHC_org = ((SUM ( HC_MTC_Reach_Wall (:,HC_Launch_Reg_Sput_Target_1))+ &
                &	     SUM ( HC_MTC_Reach_Wall (:,HC_Launch_Reg_Sput_Target_2)))/Divide_Region_1+ &
                &	     SUM ( HC_MTC_Reach_Wall (:,HC_Launch_Reg_Sput_Wall))/Divide_Region_2+ &
                &	     (SUM ( HC_MTC_Reach_Wall (:,HC_Launch_Reg_Sput_FS_Pt))+ &
                &	     SUM ( HC_MTC_Reach_Wall (:,HC_Launch_Reg_Sput_FS_2D)))/Divide_Region_3)/Divide_Region
           SCent_org = ((SUM ( HC_Num_Reach_Centre (:,HC_Launch_Reg_Sput_Target_1))+ &
                &	 SUM ( HC_Num_Reach_Centre (:,HC_Launch_Reg_Sput_Target_2)))/Divide_Region_1+ &
                &	 SUM ( HC_Num_Reach_Centre (:,HC_Launch_Reg_Sput_Wall))/Divide_Region_2+ &
                &	 (SUM ( HC_Num_Reach_Centre (:,HC_Launch_Reg_Sput_FS_Pt))+ &
                &	 SUM ( HC_Num_Reach_Centre (:,HC_Launch_Reg_Sput_FS_2D)))/Divide_Region_3)/Divide_Region
           STMax_org = ((SUM ( HC_Num_At_TMax (:,HC_Launch_Reg_Sput_Target_1))+ &
                &	 SUM ( HC_Num_At_TMax (:,HC_Launch_Reg_Sput_Target_2)))/Divide_Region_1+ &
                &	 SUM ( HC_Num_At_TMax (:,HC_Launch_Reg_Sput_Wall))/Divide_Region_2+ &
                &	 (SUM ( HC_Num_At_TMax (:,HC_Launch_Reg_Sput_FS_Pt))+ &
                &	 SUM ( HC_Num_At_TMax (:,HC_Launch_Reg_Sput_FS_2D)))/Divide_Region_3)/Divide_Region
           SStruk_org = ((SUM ( HC_Num_Striking_Target (:,HC_Launch_Reg_Sput_Target_1))+ &
                &	  SUM ( HC_Num_Striking_Target (:,HC_Launch_Reg_Sput_Target_2)))/Divide_Region_1+ &
                &	  SUM ( HC_Num_Striking_Target (:,HC_Launch_Reg_Sput_Wall))/Divide_Region_2+ &
                &	  (SUM ( HC_Num_Striking_Target (:,HC_Launch_Reg_Sput_FS_Pt))+ &
                &	  SUM ( HC_Num_Striking_Target (:,HC_Launch_Reg_Sput_FS_2D)))/Divide_Region_3)/Divide_Region
           MTCStruk_org = ((SUM ( HC_MTC_Striking_Target (:,HC_Launch_Reg_Sput_Target_1))+ &
                &	    SUM ( HC_MTC_Striking_Target (:,HC_Launch_Reg_Sput_Target_2)))/Divide_Region_1+ &
                &	    SUM ( HC_MTC_Striking_Target (:,HC_Launch_Reg_Sput_Wall))/Divide_Region_2+ &
                &	    (SUM ( HC_MTC_Striking_Target (:,HC_Launch_Reg_Sput_FS_Pt))+ &
                &	    SUM ( HC_MTC_Striking_Target (:,HC_Launch_Reg_Sput_FS_2D)))/Divide_Region_3)/Divide_Region
           SFail_org = ((SUM ( HC_Num_Failed_Launches (:,HC_Launch_Reg_Sput_Target_1))+ &
                &	 SUM ( HC_Num_Failed_Launches (:,HC_Launch_Reg_Sput_Target_2)))/Divide_Region_1+ &
                &	 SUM ( HC_Num_Failed_Launches (:,HC_Launch_Reg_Sput_Wall))/Divide_Region_2+ &
                &	 (SUM ( HC_Num_Failed_Launches (:,HC_Launch_Reg_Sput_FS_Pt))+ &
                &	 SUM ( HC_Num_Failed_Launches (:,HC_Launch_Reg_Sput_FS_2D)))/Divide_Region_3)/Divide_Region
           SatIZ_org = (( HC_Num_Fragments_Reach_CIon (HC_Launch_Reg_Sput_Target_1)+ &
                &	  HC_Num_Fragments_Reach_CIon (HC_Launch_Reg_Sput_Target_2))/Divide_Region_1+ &
                &	  HC_Num_Fragments_Reach_CIon (HC_Launch_Reg_Sput_Wall)/Divide_Region_2+ &
                &	 ( HC_Num_Fragments_Reach_CIon (HC_Launch_Reg_Sput_FS_Pt)+ &
                &	  HC_Num_Fragments_Reach_CIon (HC_Launch_Reg_Sput_FS_2D))/Divide_Region_3)/Divide_Region
           SMain_org = ((SUM ( HC_Num_Enter_Main_Plasma (:,HC_Launch_Reg_Sput_Target_1))+ &
                &	 SUM ( HC_Num_Enter_Main_Plasma (:,HC_Launch_Reg_Sput_Target_2)))/Divide_Region_1+ &
                &	 SUM ( HC_Num_Enter_Main_Plasma (:,HC_Launch_Reg_Sput_Wall))/Divide_Region_2+ &
                &	 (SUM ( HC_Num_Enter_Main_Plasma (:,HC_Launch_Reg_Sput_FS_Pt))+ &
                &	 SUM ( HC_Num_Enter_Main_Plasma (:,HC_Launch_Reg_Sput_FS_2D)))/Divide_Region_3)/Divide_Region
           SExit_org = ((SUM ( HC_Tot_Fragments_Exit_Main (:,HC_Launch_Reg_Sput_Target_1))+ &
                &	 SUM ( HC_Tot_Fragments_Exit_Main (:,HC_Launch_Reg_Sput_Target_2)))/Divide_Region_1+ &
                &	 SUM ( HC_Tot_Fragments_Exit_Main (:,HC_Launch_Reg_Sput_Wall))/Divide_Region_2+ &
                &	 (SUM ( HC_Tot_Fragments_Exit_Main (:,HC_Launch_Reg_Sput_FS_Pt))+ &
                &	 SUM ( HC_Tot_Fragments_Exit_Main (:,HC_Launch_Reg_Sput_FS_2D)))/Divide_Region_3)/Divide_Region
           SFP_org = ((SUM ( HC_RFPTarg (:,HC_Launch_Reg_Sput_Target_1))+ &
                &      SUM ( HC_RFPTarg (:,HC_Launch_Reg_Sput_Target_2)))/Divide_Region_1+ &
                &      SUM ( HC_RFPTarg (:,HC_Launch_Reg_Sput_Wall))/Divide_Region_2+ &
                &      (SUM ( HC_RFPTarg (:,HC_Launch_Reg_Sput_FS_Pt))+ &
                &      SUM ( HC_RFPTarg (:,HC_Launch_Reg_Sput_FS_2D)))/Divide_Region_3)/Divide_Region
        ElseIf (NeutType .eq. 5) Then
           ! Chemically sputtered from wall.
           If (HC_Launch_Reg_Wall_Homo .eq. HC_Launch_Reg_Wall_Dist) Then
              Divide_Region = 2.0
           Else
              Divide_Region = 1.0
           End If
     
           SHC_org = (SUM ( HC_Sum_Fragments_Launched (:,HC_Launch_Reg_Wall_Homo))+ &
                &     SUM ( HC_Sum_Fragments_Launched (:,HC_Launch_Reg_Wall_Dist)))/Divide_Region
           SWallHC_org = (SUM ( HC_Num_Reach_Wall (:,HC_Launch_Reg_Wall_Homo))+ &
                &	  SUM ( HC_Num_Reach_Wall (:,HC_Launch_Reg_Wall_Dist)))/Divide_Region
           MTCWallHC_org = (SUM ( HC_MTC_Reach_Wall (:,HC_Launch_Reg_Wall_Homo))+ &
                &	    SUM ( HC_MTC_Reach_Wall (:,HC_Launch_Reg_Wall_Dist)))/Divide_Region
           SCent_org = (SUM ( HC_Num_Reach_Centre (:,HC_Launch_Reg_Wall_Homo))+ &
                &	SUM ( HC_Num_Reach_Centre (:,HC_Launch_Reg_Wall_Dist)))/Divide_Region
           STMax_org = (SUM ( HC_Num_At_TMax (:,HC_Launch_Reg_Wall_Homo))+ &
                &	SUM ( HC_Num_At_TMax (:,HC_Launch_Reg_Wall_Homo)))/Divide_Region
           SStruk_org = (SUM ( HC_Num_Striking_Target (:,HC_Launch_Reg_Wall_Homo))+ &
                &	 SUM ( HC_Num_Striking_Target (:,HC_Launch_Reg_Wall_Dist)))/Divide_Region
           MTCStruk_org = (SUM ( HC_MTC_Striking_Target (:,HC_Launch_Reg_Wall_Homo))+ &
                &	   SUM ( HC_MTC_Striking_Target (:,HC_Launch_Reg_Wall_Dist)))/Divide_Region
           SFail_org = (SUM ( HC_Num_Failed_Launches (:,HC_Launch_Reg_Wall_Homo))+ &
                &	SUM ( HC_Num_Failed_Launches (:,HC_Launch_Reg_Wall_Dist)))/Divide_Region

           SatIZ_org = (SUM ( HC_Num_Fragments_Reach_CIon (:HC_Launch_Reg_Wall_Homo))+ &
                &	SUM ( HC_Num_Fragments_Reach_CIon (:HC_Launch_Reg_Wall_Dist)))/Divide_Region

           !SatIZ_org = (SUM ( HC_Num_Fragments_Reach_CIon (:))+ &
           !     &	SUM ( HC_Num_Fragments_Reach_CIon (:)))/Divide_Region

           SMain_org = (SUM ( HC_Num_Enter_Main_Plasma (:,HC_Launch_Reg_Wall_Homo))+ &
                &	SUM ( HC_Num_Enter_Main_Plasma (:,HC_Launch_Reg_Wall_Dist)))/Divide_Region
           SExit_org = (SUM ( HC_Tot_Fragments_Exit_Main (:,HC_Launch_Reg_Wall_Homo))+ &
                &	SUM ( HC_Tot_Fragments_Exit_Main (:,HC_Launch_Reg_Wall_Dist)))/Divide_Region
           SFP_org = (SUM ( HC_RFPTarg (:,HC_Launch_Reg_Wall_Homo))+ &
                &      SUM ( HC_RFPTarg (:,HC_Launch_Reg_Wall_Dist)))/Divide_Region
     
           !write (0,*) "NEUTTYPE HC DATA",HC_Launch_Reg_Wall_Homo,HC_Launch_Reg_Wall_Dist,SatIZ, HC_Num_Fragments_Reach_CIon,Divide_Region
     
        !ElseIf (NeutType .eq. 0) Then
           ! jdemod - neuttype=0 represents freespace launches - try treating as neuttype=6 for code purposes as seems to be done elsewhere. 
           ! Ion injection case.
           !STMax = NProd - LProd + 1
           !Write (0,*) "Ion Injection:",STMax
        ! jdemod
        ElseIf (NeutType .eq. 6.or.neuttype.eq.0) Then
           ! 2D neutral launched particles in DIVIMP.  For DIVIMP-HC, this includes both
           ! 2D and point launched particles since it is unlikely that both will be used
           ! simultaneously.
           If (HC_Launch_Reg_Free_Space_PT .eq. HC_Launch_Reg_Free_Space_2D) Then
              Divide_Region = 2.0
           Else
              Divide_Region = 1.0
           End If
     
           SHC_org = (SUM ( HC_Sum_Fragments_Launched (:,HC_Launch_Reg_Free_Space_PT))+ &
                &     SUM ( HC_Sum_Fragments_Launched (:,HC_Launch_Reg_Free_Space_2D)))/Divide_Region
           SWallHC_org = (SUM ( HC_Num_Reach_Wall (:,HC_Launch_Reg_Free_Space_PT))+ &
                &	  SUM ( HC_Num_Reach_Wall (:,HC_Launch_Reg_Free_Space_2D)))/Divide_Region
           MTCWallHC_org = (SUM ( HC_MTC_Reach_Wall (:,HC_Launch_Reg_Free_Space_PT))+ &
                &	    SUM ( HC_MTC_Reach_Wall (:,HC_Launch_Reg_Free_Space_2D)))/Divide_Region
           SCent_org = (SUM ( HC_Num_Reach_Centre (:,HC_Launch_Reg_Free_Space_PT))+ &
                &	SUM ( HC_Num_Reach_Centre (:,HC_Launch_Reg_Free_Space_2D)))/Divide_Region
           STMax_org = (SUM ( HC_Num_At_TMax (:,HC_Launch_Reg_Free_Space_PT))+ &
                &	SUM ( HC_Num_At_TMax (:,HC_Launch_Reg_Free_Space_2D)))/Divide_Region
           SStruk_org = (SUM ( HC_Num_Striking_Target (:,HC_Launch_Reg_Free_Space_PT))+ &
                &	 SUM ( HC_Num_Striking_Target (:,HC_Launch_Reg_Free_Space_2D)))/Divide_Region
           MTCStruk_org = (SUM ( HC_MTC_Striking_Target (:,HC_Launch_Reg_Free_Space_PT))+ &
                &	   SUM ( HC_MTC_Striking_Target (:,HC_Launch_Reg_Free_Space_2D)))/Divide_Region
           SFail_org = (SUM ( HC_Num_Failed_Launches (:,HC_Launch_Reg_Free_Space_PT))+ &
                &	SUM ( HC_Num_Failed_Launches (:,HC_Launch_Reg_Free_Space_2D)))/Divide_Region
           SatIZ_org = ( HC_Num_Fragments_Reach_CIon (HC_Launch_Reg_Free_Space_PT)+ &
                &	 HC_Num_Fragments_Reach_CIon (HC_Launch_Reg_Free_Space_2D))/Divide_Region
           SMain_org = (SUM ( HC_Num_Enter_Main_Plasma (:,HC_Launch_Reg_Free_Space_PT))+ &
                &	SUM ( HC_Num_Enter_Main_Plasma (:,HC_Launch_Reg_Free_Space_2D)))/Divide_Region
           SExit_org = (SUM ( HC_Tot_Fragments_Exit_Main (:,HC_Launch_Reg_Free_Space_PT))+ &
                &	SUM ( HC_Tot_Fragments_Exit_Main (:,HC_Launch_Reg_Free_Space_2D)))/Divide_Region
           SFP_org = (SUM ( HC_RFPTarg (:,HC_Launch_Reg_Free_Space_PT))+ &
                &      SUM ( HC_RFPTarg (:,HC_Launch_Reg_Free_Space_2D)))/Divide_Region
        ElseIf (NeutType .eq. 7) Then
           ! NeutType indicates target launched and later reflected.  For DIVIMP-HC,
           ! this would be all reflected particles including those from the wall.
     
           ! Divide_Region_1 = Target chemical launches
           ! Divide_Region_2 = Wall chemical launches
           ! Divide_Region_3 = FS chemical launches
           Divide_Region = 1.0
           Divide_Region_1 = 1.0
           Divide_Region_2 = 1.0
           Divide_Region_3 = 1.0
     
           If (HC_Launch_Reg_Refl_Target_1 .eq. HC_Launch_Reg_Refl_Target_2 .and. &
                &   HC_Launch_Reg_Refl_Target_1 .eq. HC_Launch_Reg_Refl_Wall .and. &
                &   HC_Launch_Reg_Refl_Target_1 .eq. HC_Launch_Reg_Refl_FS_Pt .and. &
                &   HC_Launch_Reg_Refl_Target_1 .eq. HC_Launch_Reg_Refl_FS_2D) Then
              Divide_Region = 5.0
           ElseIf (HC_Launch_Reg_Refl_Target_1 .eq. HC_Launch_Reg_Refl_Target_2) Then
              Divide_Region_1 = 2.0
           ElseIf (HC_Launch_Reg_Refl_FS_PT .eq. HC_Launch_Reg_Refl_FS_2D) Then
              Divide_Region_3 = 2.0
           End If
     
           SHC_org = ((SUM ( HC_Sum_Fragments_Launched (:,HC_Launch_Reg_Refl_Target_1))+ &
                &      SUM ( HC_Sum_Fragments_Launched (:,HC_Launch_Reg_Refl_Target_2)))/Divide_Region_1+ &
                &      SUM ( HC_Sum_Fragments_Launched (:,HC_Launch_Reg_Refl_Wall))/Divide_Region_2+ &
                &      (SUM ( HC_Sum_Fragments_Launched (:,HC_Launch_Reg_Refl_FS_Pt))+ &
                &      SUM ( HC_Sum_Fragments_Launched (:,HC_Launch_Reg_Refl_FS_2D)))/Divide_Region_3)/Divide_Region
           SWallHC_org = ((SUM ( HC_Num_Reach_Wall (:,HC_Launch_Reg_Refl_Target_1))+ &
                &	   SUM ( HC_Num_Reach_Wall (:,HC_Launch_Reg_Refl_Target_2)))/Divide_Region_1+ &
                &	   SUM ( HC_Num_Reach_Wall (:,HC_Launch_Reg_Refl_Wall))/Divide_Region_2+ &
                &	   (SUM ( HC_Num_Reach_Wall (:,HC_Launch_Reg_Refl_FS_Pt))+ &
                &	   SUM ( HC_Num_Reach_Wall (:,HC_Launch_Reg_Refl_FS_2D)))/Divide_Region_3)/Divide_Region
           MTCWallHC_org = ((SUM ( HC_MTC_Reach_Wall (:,HC_Launch_Reg_Refl_Target_1))+ &
                &	     SUM ( HC_MTC_Reach_Wall (:,HC_Launch_Reg_Refl_Target_2)))/Divide_Region_1+ &
                &	     SUM ( HC_MTC_Reach_Wall (:,HC_Launch_Reg_Refl_Wall))/Divide_Region_2+ &
                &	     (SUM ( HC_MTC_Reach_Wall (:,HC_Launch_Reg_Refl_FS_Pt))+ &
                &	     SUM ( HC_MTC_Reach_Wall (:,HC_Launch_Reg_Refl_FS_2D)))/Divide_Region_3)/Divide_Region
           SCent_org = ((SUM ( HC_Num_Reach_Centre (:,HC_Launch_Reg_Refl_Target_1))+ &
                &	 SUM ( HC_Num_Reach_Centre (:,HC_Launch_Reg_Refl_Target_2)))/Divide_Region_1+ &
                &	 SUM ( HC_Num_Reach_Centre (:,HC_Launch_Reg_Refl_Wall))/Divide_Region_2+ &
                &	 (SUM ( HC_Num_Reach_Centre (:,HC_Launch_Reg_Refl_FS_Pt))+ &
                &	 SUM ( HC_Num_Reach_Centre (:,HC_Launch_Reg_Refl_FS_2D)))/Divide_Region_3)/Divide_Region
           STMax_org = ((SUM ( HC_Num_At_TMax (:,HC_Launch_Reg_Refl_Target_1))+ &
                &	 SUM ( HC_Num_At_TMax (:,HC_Launch_Reg_Refl_Target_2)))/Divide_Region_1+ &
                &	 SUM ( HC_Num_At_TMax (:,HC_Launch_Reg_Refl_Wall))/Divide_Region_2+ &
                &	 (SUM ( HC_Num_At_TMax (:,HC_Launch_Reg_Refl_FS_Pt))+ &
                &	 SUM ( HC_Num_At_TMax (:,HC_Launch_Reg_Refl_FS_2D)))/Divide_Region_3)/Divide_Region
           SStruk_org = ((SUM ( HC_Num_Striking_Target (:,HC_Launch_Reg_Refl_Target_1))+ &
                &	  SUM ( HC_Num_Striking_Target (:,HC_Launch_Reg_Refl_Target_2)))/Divide_Region_1+ &
                &	  SUM ( HC_Num_Striking_Target (:,HC_Launch_Reg_Refl_Wall))/Divide_Region_2+ &
                &	  (SUM ( HC_Num_Striking_Target (:,HC_Launch_Reg_Refl_FS_Pt))+ &
                &	  SUM ( HC_Num_Striking_Target (:,HC_Launch_Reg_Refl_FS_2D)))/Divide_Region_3)/Divide_Region
           MTCStruk_org = ((SUM ( HC_MTC_Striking_Target (:,HC_Launch_Reg_Refl_Target_1))+ &
                &	    SUM ( HC_MTC_Striking_Target (:,HC_Launch_Reg_Refl_Target_2)))/Divide_Region_1+ &
                &	    SUM ( HC_MTC_Striking_Target (:,HC_Launch_Reg_Refl_Wall))/Divide_Region_2+ &
                &	    (SUM ( HC_MTC_Striking_Target (:,HC_Launch_Reg_Refl_FS_Pt))+ &
                &	    SUM ( HC_MTC_Striking_Target (:,HC_Launch_Reg_Refl_FS_2D)))/Divide_Region_3)/Divide_Region
           SFail_org = ((SUM ( HC_Num_Failed_Launches (:,HC_Launch_Reg_Refl_Target_1))+ &
                &	 SUM ( HC_Num_Failed_Launches (:,HC_Launch_Reg_Refl_Target_2)))/Divide_Region_1+ &
                &	 SUM ( HC_Num_Failed_Launches (:,HC_Launch_Reg_Refl_Wall))/Divide_Region_2+ &
                &	 (SUM ( HC_Num_Failed_Launches (:,HC_Launch_Reg_Refl_FS_Pt))+ &
                &	 SUM ( HC_Num_Failed_Launches (:,HC_Launch_Reg_Refl_FS_2D)))/Divide_Region_3)/Divide_Region
           SatIZ_org = (( HC_Num_Fragments_Reach_CIon (HC_Launch_Reg_Refl_Target_1)+ &
                &	  HC_Num_Fragments_Reach_CIon (HC_Launch_Reg_Refl_Target_2))/Divide_Region_1+ &
                &	  HC_Num_Fragments_Reach_CIon (HC_Launch_Reg_Refl_Wall)/Divide_Region_2+ &
                &	 ( HC_Num_Fragments_Reach_CIon (HC_Launch_Reg_Refl_FS_Pt)+ &
                &	  HC_Num_Fragments_Reach_CIon (HC_Launch_Reg_Refl_FS_2D))/Divide_Region_3)/Divide_Region
           SMain_org = ((SUM ( HC_Num_Enter_Main_Plasma (:,HC_Launch_Reg_Refl_Target_1))+ &
                &	 SUM ( HC_Num_Enter_Main_Plasma (:,HC_Launch_Reg_Refl_Target_2)))/Divide_Region_1+ &
                &	 SUM ( HC_Num_Enter_Main_Plasma (:,HC_Launch_Reg_Refl_Wall))/Divide_Region_2+ &
                &	 (SUM ( HC_Num_Enter_Main_Plasma (:,HC_Launch_Reg_Refl_FS_Pt))+ &
                &	 SUM ( HC_Num_Enter_Main_Plasma (:,HC_Launch_Reg_Refl_FS_2D)))/Divide_Region_3)/Divide_Region
           SExit_org = ((SUM ( HC_Tot_Fragments_Exit_Main (:,HC_Launch_Reg_Refl_Target_1))+ &
                &	 SUM ( HC_Tot_Fragments_Exit_Main (:,HC_Launch_Reg_Refl_Target_2)))/Divide_Region_1+ &
                &	 SUM ( HC_Tot_Fragments_Exit_Main (:,HC_Launch_Reg_Refl_Wall))/Divide_Region_2+ &
                &	 (SUM ( HC_Tot_Fragments_Exit_Main (:,HC_Launch_Reg_Refl_FS_Pt))+ &
                &	 SUM ( HC_Tot_Fragments_Exit_Main (:,HC_Launch_Reg_Refl_FS_2D)))/Divide_Region_3)/Divide_Region
           SFP_org = ((SUM ( HC_RFPTarg (:,HC_Launch_Reg_Refl_Target_1))+ &
                &      SUM ( HC_RFPTarg (:,HC_Launch_Reg_Refl_Target_2)))/Divide_Region_1+ &
                &      SUM ( HC_RFPTarg (:,HC_Launch_Reg_Refl_Wall))/Divide_Region_2+ &
                &      (SUM ( HC_RFPTarg (:,HC_Launch_Reg_Refl_FS_Pt))+ &
                &      SUM ( HC_RFPTarg (:,HC_Launch_Reg_Refl_FS_2D)))/Divide_Region_3)/Divide_Region
        Else			
           ! NeutType not allowed.
           Write (Output_Unit_HC_Alert,*) "Error in HC_Batch: NeutType value not allowed in DIVIMP-HC:",NeutType
           Write (Output_Unit_HC_Alert,*) "Program stopping."
           !Stop
        End If
    ! 
    ! jdemod - end 
    !



    ! This will give the increment from the current generation. 
    ! NOTE: Need to verify that the code doesn't double count in some circumstances since there is the 
    !       definition of Divide_region in the original code. 

    ! jdemod - assign appropriate values to these before paassing back to calling routine - use the _org values
    !          since they appear to be correctly calculated ... the other methods include reflected and other particles
    !          that appear to be treated as separate launch regions. 

    !SHC      = SHC      - SHC_init
    !SSTruk   = SSTruk   - SSTruk_init   
    !MTCStruk = MTCStruk - MTCStruk_init 
    !SMain    = SMain    - SMain_init    
    !SExit    = SExit    - SExit_init    
    !SatIZ    = SatIZ    - SatIZ_init    
    !SWallHC  = SWallHC  - SWallHC_init  
    !MTCWallHC= MTCWallHC- MTCWallHC_init
    !SCent    = SCent    - SCent_init    
    !STMax    = STMax    - STMax_init    
    !SFail    = SFail    - SFail_init    
    !SFP      = SFP      - SFP_init      

    SHC      =  SHC_org
    SSTruk   =  SSTruk_org   
    MTCStruk =  MTCStruk_org 
    SMain    =  SMain_org    
    SExit    =  SExit_org    
    SatIZ    =  SatIZ_org    
    SWallHC  =  SWallHC_org  
    MTCWallHC=  MTCWallHC_org
    SCent    =  SCent_org    
    STMax    =  STMax_org    
    SFail    =  SFail_org    
    SFP      =  SFP_org      


    ! Output sums as calculated.
    If (NeutType .ne. 3) Then
       Write (Output_Unit_Scratch,'(A,2I3)') "Chemically sputtered HC launch complete. NEUTTYPE, STATUS:",NeutType,Status
    Else
       Write (Output_Unit_Scratch,'(A,2I3)') "Self-sputtered HC launch complete. NEUTTYPE, STATUS:",NeutType,Status
    End If
    Write (Output_Unit_Scratch,*) "Prod:",nprod-lprod+1
    Write (Output_Unit_Scratch,'(a,5(1x,g18.8),l6)') "SHC:",SHC, SUM(HC_Sum_Fragments_Launched(:,1:4)), SUM(HC_Sum_Fragments_Launched(:,5:8)),SHC_init, SHC_org, SHC.eq.SHC_org
    Write (Output_Unit_Scratch,'(a,3(1x,g18.8),l6)') "SStruk:",SStruk,        SSTruk_init   ,  SSTruk_org   , SSTruk.eq.SSTruk_org  
    Write (Output_Unit_Scratch,'(a,3(1x,g18.8),l6)') "MTCStruk:",MTCStruk,    MTCStruk_init ,  MTCStruk_org , MTCStruk .eq. MTCStruk_org 
    Write (Output_Unit_Scratch,'(a,3(1x,g18.8),l6)') "SMain:",SMain,          SMain_init    ,  SMain_org    , SMain    .eq. SMain_org    
    Write (Output_Unit_Scratch,'(a,3(1x,g18.8),l6)') "SExit:",SExit,          SExit_init    ,  SExit_org    , SExit    .eq. SExit_org    
    Write (Output_Unit_Scratch,'(a,3(1x,g18.8),l6)') "SatIZ:",SatIZ,          SatIZ_init    ,  SatIZ_org    , SatIZ    .eq. SatIZ_org    
    Write (Output_Unit_Scratch,'(a,3(1x,g18.8),l6)') "SWallHC:",SWallHC,      SWallHC_init  ,  SWallHC_org  , SWallHC  .eq. SWallHC_org  
    Write (Output_Unit_Scratch,'(a,3(1x,g18.8),l6)') "MTCWallHC:",MTCWallHC,  MTCWallHC_init,  MTCWallHC_org, MTCWallHC.eq. MTCWallHC_org
    Write (Output_Unit_Scratch,'(a,3(1x,g18.8),l6)') "SCent:",SCent,          SCent_init    ,  SCent_org    , SCent    .eq. SCent_org    
    Write (Output_Unit_Scratch,'(a,3(1x,g18.8),l6)') "STMax:",STMax,          STMax_init    ,  STMax_org    , STMax    .eq. STMax_org    
    Write (Output_Unit_Scratch,'(a,3(1x,g18.8),l6)') "SFail:",SFail,          SFail_init    ,  SFail_org    , SFail    .eq. SFail_org    
    Write (Output_Unit_Scratch,'(a,3(1x,g18.8),l6)') "SFP:",SFP,              SFP_init      ,  SFP_org      , SFP      .eq. SFP_org      


    ! Check that all sums are still whole numbers.  If not, there was likely a divide error somewhere.
    ! jdemod - if partial weight particles are in use from sputtering then these may not be whole numbers
    !
    If (MOD (SHC,1.0) .ne. 0.0) Then
       Write (Output_Unit_HC_Alert,*) "Warning in HC_Batch:  Sum of SHC not a whole number:",SHC
    ElseIf (MOD (SWallHC,1.0) .ne. 0.0) Then
       Write (Output_Unit_HC_Alert,*) "Warning in HC_Batch:  Sum of SWallHC not a whole number:",SWallHC	
    ElseIf (MOD (MTCWallHC,1.0) .ne. 0.0) Then
       Write (Output_Unit_HC_Alert,*) "Warning in HC_Batch:  Sum of MTCWallHC not a whole number:",MTCWallHC
    ElseIf (MOD (SCent,1.0) .ne. 0.0) Then
       Write (Output_Unit_HC_Alert,*) "Warning in HC_Batch:  Sum of SCent not a whole number:",SCent
    ElseIf (MOD (STMax,1.0) .ne. 0.0) Then
       Write (Output_Unit_HC_Alert,*) "Warning in HC_Batch:  Sum of STMax not a whole number:",STMax
    ElseIf (MOD (SStruk,1.0) .ne. 0.0) Then
       Write (Output_Unit_HC_Alert,*) "Warning in HC_Batch:  Sum of SStruk not a whole number:",SStruk
    ElseIf (MOD (MTCStruk,1.0) .ne. 0.0) Then
       Write (Output_Unit_HC_Alert,*) "Warning in HC_Batch:  Sum of MTCStruk not a whole number:",MTCStruk
    ElseIf (MOD (SFail,1.0) .ne. 0.0) Then
       Write (Output_Unit_HC_Alert,*) "Warning in HC_Batch:  Sum of SFail not a whole number:",SFail
    ElseIf (MOD (SatIZ,1.0) .ne. 0.0) Then
       Write (Output_Unit_HC_Alert,*) "Warning in HC_Batch:  Sum of SatIZ not a whole number:",SatIZ
    ElseIf (MOD (SMain,1.0) .ne. 0.0) Then
       Write (Output_Unit_HC_Alert,*) "Warning in HC_Batch:  Sum of SMain not a whole number:",SMain
    ElseIf (MOD (SExit,1.0) .ne. 0.0) Then
       Write (Output_Unit_HC_Alert,*) "Warning in HC_Batch:  Sum of SExit not a whole number:",SExit
    ElseIf (MOD (SFP,1.0) .ne. 0.0) Then
       Write (Output_Unit_HC_Alert,*) "Warning in HC_Batch:  Sum of SFP not a whole number:",SFP
    End If

    ! Note:  This sum will not add up for sputtered and reflected sums as the statistics may
    ! be added to by the DIVIMP-HC code during other types of launches.  Also, more than one
    ! group of self-sputtered particles may be launches (2nd + 3rd +... generations) so
    ! sputtered statistics will sum from generation to generation.
    If (NeutType .ne. 3 .and. NeutType .ne. 7) Then

       ! Check that the sum of all sums is the same as the total for particles launched.
       If (SHC .ne. REAL (NProd - LProd + 1)) Then
          Write (Output_Unit_HC_Alert,*) "Warning 1 in HC_Batch: Total particle sums launched from HC_Transport not consistent: ",&
               &REAL (NProd - LProd + 1),SHC
       End If

       ! Check that the sum of particles returned is the same as the total for particles lauched
       Launch_Sum = SWallHC + MTCWallHC + SCent + STMax + SStruk + MTCStruk + SFail + SatIZ + SMain + SExit + SFP
       If (Launch_Sum .ne. REAL (NProd - LProd + 1)) Then
          Write (Output_Unit_HC_Alert,*) "Warning 2 in HC_Batch: Total particle sums returned from HC_Transport not consistent: ",&
               &REAL (NProd - LProd + 1),Launch_Sum
          Write (Output_Unit_HC_Alert,'(a,20(1x,g15.5))') "Warning 2 Details:",real(NProd-LProd+1),SWallHC, MTCWallHC, SCent, STMax, SStruk,&
                                                                             & MTCStruk, SFail, SatIZ, SMain, SExit, SFP, SHC
       End If
    End If

    ! Finish timing for this group of hydrocarbons launched.
    NeuTime = (NeuTime + IonTime) + ZA02AS(1) - Starting_Time ! Hydrocarbon 'neutim' is time spent following all HC charge states until C+ production.

    !NATIZ = 0 ! AMMOD REMOVE



    write(0,*) 'END Of HC_LAUNCH:',nprod,natiz


  End Subroutine HC_Launch

End Module HC_Batch
