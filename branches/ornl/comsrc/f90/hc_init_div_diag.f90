! -*-Mode:f90-*-
! Initialize_Diagnostics.f90
! Diagnostic Data Storage Initialization File
!
! Adam McLean under Prof. Peter Stangeby
! Research Assistant: David Elder
! January, 2002
!
! This module initializes each of the required data structures for
! all cell properties.  It acts as the data-passing front-end between
! original DIVIMP code and all the added transport and hydrocarbon
! routines added.

Module HC_Init_DIV_Diag

  ! Use statements.
  Use HC_Storage_Setup ! Access to derived data types.	

  ! Every good Fortran 90 program has...
  Implicit None

  ! Define data tables.
  !Type (HC_Launch_Diag_Table_Type) :: HC_Launch_Diag_Table
  !Type (HC_Transport_Diag_Table_Type) :: HC_Transport_Diag_Table
  !Type (HC_Evolve_Diag_Table_Type) :: HC_Evolve_Diag_Table
  !Type (HC_Rerelease_Diag_Table_Type) :: HC_Rerelease_Diag_Table
  !Type (HC_Vessel_Int_Diag_Table_Type) :: HC_Vessel_Int_Diag_Table
  !Type (HC_Death_Diag_Table_Type) :: HC_Death_Diag_Table

Contains

  Subroutine Initialize_HC_Launch_Diag_Table

    ! Required modules.
    Use ComHC ! Contains Number_HC_Species.
    Use HC_Init_DIV_Data ! Contains Global_Prop_Table.
    Use HC_Init_Lib_Data ! Contains HC_State_Table.

    ! Every good Fortran 90 program has...
    Implicit None

    ! Declare local variables.
    Integer :: i

    ! Note: Arrays are defined in a data-parallel way (i.e. the unary
    ! operator is assigned to all array positions at once).

    HC_Total_R_Prod_Positions = 0.0 ! XTOT 
    HC_Total_Z_Prod_Positions = 0.0 ! YTOT
    HC_Total_S_Prod_Positions = 0.0 ! 
    HC_Total_Cross_Prod_Positions = 0.0 ! 
    HC_Total_Prod_Angles = 0.0 ! ATOT
    HC_Sum_Fragments_Launched = 0.0 ! RNEUT
    HC_Total_Prod_Vels_No_VMult = 0.0 ! VTOTA
    HC_Max_Prod_Vel_No_VMult = 0.0 ! VTOTAM
    !		 HC_Total_Vel_Ang_Mults = 0.0 ! VMULTT
    !		 HC_Max_Vel_Ang_Mults = 0.0 ! VMULTM
    HC_Total_Prod_Velocities = 0.0 ! VTOT
    HC_Max_Prod_Velocities = 0.0 ! VTOTM
    HC_Tot_Prod_Temperatures = 0.0 ! ETOT
    HC_Num_Failed_Launches = 0.0 ! RFAIL 
    HC_Num_Vel_Greater_Max_Ran = 0.0 ! 
    HC_Launch_Grid_Error_Moved_Okay = 0.0 ! ERR_OK
    HC_Launch_Grid_Error_Moved_Out = 0.0 ! ERR_OUT
    HC_LaunchDat = 0.0 ! LAUNCHDAT 

    ! Fill the state table in the storage area to pass to OUT.
    Do i = 1,Number_HC_Species,1
       HC_State_List (i) = HC_State_Table (i) % State_Name
    End Do

    ! Fill in HC_Output_List.


  End Subroutine Initialize_HC_Launch_Diag_Table

  Subroutine Initialize_Transport_Diag_Table

    ! Set RSTMax: Max number of iterations up to time = TMax (0.1S)
    ! (differs from CSTMax which applies to ion,  by FSRate/QTim and
    ! because TMax is 0.1 here but 10 for ions).
    ! Set NatIZ : Number of hydrocarbon "fragments" that evolve to C+.
    ! Set RProd : Sum of fragments produced.

    ! Required modules.
    Use HC_Init_DIV_Data ! Contains Global_Prop_Table.
    Use HC_Get

    ! Every good Fortran 90 program has...
    Implicit None

    ! Declare local variables.
    Integer :: Current_Ring
    Integer :: Current_Cell

    ! Timing data.
    HC_Walk_Count = 1
    HC_Neutral_Timesteps_Count = 0
    HC_Ion_Timesteps_Count = 0
    Time_All_HCs = 0.0 ! CISTOT
    Max_Time_Any_HC = 0.0 ! CISMAX
    HC_State_Time = 0.0
    HC_Max_Ion_Teq_Iter = 0.0
    HC_Total_Ion_Teq_Iter = 0.0
    HC_Max_Time_In_State = 0.0
    HC_State_TimeSteps = 0 ! CIEIZS
    HC_Max_TimeSteps = 0
    HC_Tot_Time_First_Diff = 0.0 ! RDIFFT
    HC_Num_Reach_Max_Iter = 0.0 ! CICUTS
    HC_Tot_Exist_Max_Iter = 0.0 ! TCUT

    ! State sums.
    HC_Num_Orig_Reach_State = 0.0 ! CITIZS
    HC_Tot_Reach_State = 0.0 ! 

    ! Movement sums.
    HC_Num_Striking_Target = 0.0 ! RSTRUK
    HC_Num_Enter_Main_Plasma = 0.0 ! RMAIN 
    HC_Tot_Fragments_Exit_Main = 0.0 ! REXIT
    HC_Num_Entered_Core = 0.0 ! NUM_ENTERED_CORE
    HC_Num_Reach_Centre = 0.0 ! RCENT
    HC_Stopped_Follow_Ion_In_Core = 0.0 ! STOPPED_FOLLOW

    ! Position indicators.
    HC_Walks = 0.0 ! WALKS
    HC_InMain = .False.
    HC_InCore = .False.
    HC_InEdge = .False.
    HC_InMSol = .False.
    HC_InDiv = .False.
    HC_InTrap = .False.
    HC_CflRxa = .False.
    HC_CflRex = .False.
    HC_CflRin = .False.

    ! Distance totals.
    HC_Cumu_CroFly_Dist_Trav_By_Pos = 0.0
    HC_Cumu_Real_Dist_Trav_By_Pos = 0.0
    HC_Cumu_Real_Dist_Trav_By_VT = 0.0		   
    HC_State_Dist_Trav_By_Pos = 0.0
    HC_State_Dist_Trav_By_VT = 0.0

    ! Density recording by current position.
    HC_Density = 0.0
    H_Density = 0.0
    HC_ChemDen = 0.0 ! CHEMDEN
    HC_ChemIzs = 0.0 ! CHEMIZS
    HC_DDVoid = 0.0 ! DDVOID
    HC_ELims = 0.0 ! ELIMS
    HC_Lims = 0.0 ! LIMS
    HC_DParas = 0 ! DPARAS

    ! C and CH/CD neutral data for spectroscopic analysis.
    HC_TIZS = 0.0 ! TIZS
    HC_TIZS_CH = 0.0 ! TIZS
    !                HC_TIZS_C2 = 0.0 ! C2 TIZS
    ! jdemod - removed all except hc_ddts 
    !HC_DDLims = 0.0 ! DDLIMS for C band emission.
    !HC_DDLims_CH = 0.0 ! DDLIMS for CH/CD band emission.
    HC_DDts = 0.0 ! DDTS
    !HC_DDts_CH = 0.0 ! DDTS

    ! Density recording by starting position.
    HC_Ion_Core_Density = 0.0 ! NCORE
    HC_Ion_Edge_Density = 0.0 ! NEDGE
    HC_Ion_Trap_Density = 0.0 ! NTRAP
    HC_Ion_Divertor_Density = 0.0 ! NDIVERT
    HC_Ion_MSOL_Density = 0.0 ! NMSOL

    ! MTC, Momentum Transfer Collisions.
    HC_MTC_Striking_Target = 0.0 ! MTC_RSTRUK
    HC_MTCInf = 0.0 ! MTCINF
    HC_MTCTotCnt = 0.0 ! MTCTOTCNT

    ! Leakage statistics.
    Number_Leaked_Core = 0 ! NLEAKCORE
    Total_Leaked_Core = 0.0 ! TOTLEAKCORE
    HC_Has_Leaked = .False. ! HASLEAKED
    HC_Has_Leaked_Core = .False. ! HASLEAKEDCORE
    HC_Leak_Position = 0.0 ! CLEAKPOS
    HC_Leak_Time = 0.0 ! CLEAKT
    HC_Leak_Density = 0.0 ! CLEAKN
    HC_Leak_Particles = 0 ! CLEAKP

    ! Far-periphery counting.
    HC_FPTarg = 0.0 ! FPTARG
    HC_FPTart = 0.0 ! FPTART
    HC_RFPTarg = 0.0 ! RFPTARG
    HC_FPEnt = 0.0 ! FPEnt
    HC_FPExit = 0.0 ! FPExit
    HC_FPTTotal = 0.0 ! FPTTOT
    HC_Tot_Core_From_FP_Ref = 0.0 ! CVVFPREF
    HC_Tot_Core_From_FP_NoRef = 0.0 ! CVVFPNRF

    ! Counters used during ion transport.
    HC_CoreOuts = 0 ! COREOUTS
    HC_DSParaNorm = 0 ! DSPARANORM
    HC_DVParaNorm = 0 ! DVPARANORM
    HC_DSParaStep = 0 ! DSPARASTEP
    HC_DVParaStep = 0 ! DVPARASTEP
    HC_VParaStep = 0 ! VPARASTEP
    HC_VParaNorm = 0 ! VPARASTEP
    HC_DSParaCnt = 0 ! DSPARACNT
    HC_DVParaCnt = 0 ! DVPARACNT
    HC_DVMinV = 0 ! DVMINV
    HC_DVMaxV = 0 ! DVMAXV

    ! Force counters for ion transport.
    HC_FCell = 0.0 ! FCELL
    HC_Ffi = 0.0 ! FFI
    HC_Fthi = 0.0 ! FTHI
    HC_Fvbg = 0.0 ! FVBG
    HC_Force_Diff = 0.0 ! FORCE_DIFF
    HC_IonVelAvg = 0.0 ! IONVELAVG

    ! Velocity statistics.
    HC_SDVS = 0.0
    HC_SDVS2 = 0.0
    HC_SDVS3 = 0.0
    HC_VelSpace = 0.0
    HC_VelWeight = 0.0

    ! Temperature statistics.
    HC_Tot_Temp_Reach_Max_Iter = 0.0 ! CRTRCS

    ! Movement indicators.
    HC_Num_Orig_Enter_Main = 0.0 ! CLLL
    HC_Tot_Z_Orig_Enter_Main = 0.0 ! CLLLX
    HC_Tot_S_Orig_Enter_Main = 0.0 ! CLLLS
    HC_Num_Orig_Enter_Core = 0.0 ! CNNN
    HC_Tot_Z_Orig_Enter_Core = 0.0 ! CNNNX
    HC_Tot_S_Orig_Enter_Core = 0.0 ! CNNNS
    HC_Tot_Teq_Orig_Enter_Core = 0.0 ! CNNNT
    HC_Tot_Min_Z_Reach = 0.0 ! CXXX 
    HC_Tot_Max_S_Reach = 0.0 ! CSSS
    HC_Tot_Max_S_Reach_Ion_Removal = 0.0 ! CSSSS
    HC_Tot_Core_From_Reg_Ref = 0.0 ! CVVFPREF
    HC_Tot_Core_From_Reg_NoRef = 0.0 ! CVVFPNRF
    HC_Tot_R_Species_Creation = 0.0 ! CVVXC
    HC_Tot_Z_At_Main_Entry = 0.0 ! CVVXE
    HC_Num_Orig_Central_Reflect = 0.0 ! CICRXA
    HC_Tot_Teq_Orig_Central_Reflect = 0.0 ! CISRXA
    HC_Tot_Temp_At_Central_Reflect = 0.0 ! CITRXA
    HC_Min_Teq_Orig_Central_Reflect = 1.0 ! CIFRXA
    HC_Tot_Ion_Collision = 0.0 ! CICCOL
    HC_Tot_Start_Main = 0.0 ! CICRNJ
    HC_Min_Teq_Orig_Enter_Main = 1.0 ! CIFRIN
    HC_Tot_Teq_Orig_Enter_Main = 0.0 ! CISRIN
    HC_Tot_Temp_Orig_Enter_Main = 0.0 ! CITRIN
    HC_Num_Ions_Lost = 0.0 ! CICLOS
    HC_Min_Teq_Ions_Lost = 1.0 ! CIFLOS
    HC_Max_Teq_Ions_Lost = 0.0 ! CILLOS
    HC_Tot_Teq_Ions_Lost = 0.0 ! CISLOS

  End Subroutine Initialize_Transport_Diag_Table

  Subroutine Initialize_Evolve_Diag_Table

    ! Every good Fortran 90 program has...
    Implicit None

    HC_Reaction_Count = 0.0 ! Reaction statistics counter.
    HC_LIonizDat = 0 ! LIONIZDAT
    HC_Num_Fragments_Reach_CIon = 0.0 ! NATIZ

    HC_Factor_A = 0.0 ! FACTA
    HC_Factor_B = 0.0 ! FACTB

    ! Assign arrays for all species and target of origin.
    HC_R_Pos_At_Prod = 0.0 ! XATIZ
    HC_Z_Pos_At_Prod = 0.0 ! YATIZ
    HC_S_Pos_At_Prod = 0.0 ! 
    HC_Cross_Pos_At_Prod = 0.0 ! 
    HC_Angle_At_Prod = 0.0 ! ATOT
    HC_Num_Fragments = 0.0 ! RNEUT
    HC_Velocity_At_Prod = 0.0 ! VATIZ
    HC_Max_Velocity_At_Prod = 0.0 ! VATIZM
    HC_Temperature_At_Prod = 0.0 ! EATIZ
    HC_Time_To_Prod = 0.0 ! TATIZ	
    HC_Time_At_Production = 0.0
    HC_Energy_At_Production = 0.0

  End Subroutine Initialize_Evolve_Diag_Table

  Subroutine Initialize_Rerelease_Diag_Table

    ! Every good Fortran 90 program has...
    Implicit None

    HC_Sum_Fragments_ReLaunched = 0.0 ! RNEUT

    ! Assign arrays for all species and target of origin.
    HC_Total_R_ReProd_Positions = 0.0 ! 
    HC_Total_Z_ReProd_Positions = 0.0 ! 
    HC_Total_ReProd_Angles = 0.0 ! 
    HC_Total_ReProd_Vels_No_VMult = 0.0 ! 
    HC_Max_ReProd_Vel_No_VMult = 0.0 ! 
    !		 HC_Total_Vel_Ang_Mults = 0.0 ! 
    !		 HC_Max_Vel_Ang_Mults = 0.0 ! 
    HC_Total_ReProd_Velocities = 0.0 ! 
    HC_Max_ReProd_Velocities = 0.0 ! 
    HC_Tot_ReProd_Temperatures = 0.0 ! 
    HC_ReLau_Grid_Error_Moved_Okay = 0.0 ! 
    HC_ReLau_Grid_Error_Moved_Out = 0.0 ! 

    Total_HC_Reflections = 0 ! TOTRF
    Max_HC_Reflections_Found = 0 ! MAXRF
    HC_Reflection_Loss = 0.0 ! NRFLOSS

    HC_YldTot = 0.0 ! YLDTOT
    HC_YThTot = 0.0 ! YTHTOT
    HC_YldMax = 0.0 ! YLDMAX

  End Subroutine Initialize_Rerelease_Diag_Table

  Subroutine Initialize_VesselInt_Diag_Table

    ! Every good Fortran 90 program has...
    Implicit None

    HC_Deposit = 0.0
    HC_Erosion = 0.0
    HC_PromptDeps = 0.0
    HC_Walls = 0.0
    HC_WallsE = 0.0
    HC_WallsE_I = 0.0
    HC_WallsI = 0.0
    HC_WallsN = 0.0
    HC_WTDep = 0.0
    HC_WTSource = 0.0

    ! Assign arrays for all species and target of origin.
    HC_RWall = 0.0 ! RWALL
    HC_RDep = 0.0 ! RDEP
    HC_Num_Reach_Wall = 0.0 ! RWALLN
    HC_MTC_Reach_Wall = 0.0 ! MTC_RWALLN

  End Subroutine Initialize_VesselInt_Diag_Table

  Subroutine Initialize_Death_Diag_Table

    ! Every good Fortran 90 program has...
    Implicit None

    HC_IFate_Count = 0 ! Record particle IFate data.
    HC_Wall_Deposit_Count = 0 ! Records deposit events for each wall, each HC species.
    HC_Carbon_Wall_Deposit_Count = 0 ! Total number of carbon atoms depositing on each wall segment.

    ! jdemod - HC_NUM_at_Tmax initialized twice - second just below - restore initialization of HC_Ctexs
    ! HC_Num_At_TMax = 0.0 ! RTMAX
    HC_Ctexs = 0.0 ! CTEXS

    ! Assign arrays for all species and target of origin.
    HC_Num_At_TMax = 0.0 ! RTMAX
    HC_Num_Absorbed_TargWall = 0.0 ! CICABS
    HC_Min_Teq_To_Absorption = 1.0 ! CIFABS
    HC_Max_Teq_To_Absorption = 0.0 ! CILABS
    HC_Tot_Teq_To_Absorption = 0.0 ! CISABS
    HC_Tot_Temp_At_Absorption = 0.0 ! CRTABS
    HC_Tot_Velocity_At_Absorption = 0.0 ! CRVABS
    HC_Tot_ABS_Vel_At_Absorption = 0.0 ! CRAVAV
    HC_Tot_Elec_Temp_At_Absorption = 0.0 ! CTBS
    HC_Num_Absorbed_Act_Target = 0.0 ! ACTTARG

    HC_Tot_At_WBC_Boundary = 0.0
    HC_Tot_Temp_At_WBC_Boundary = 0.0
    HC_Tot_At_WBC_Target = 0.0
    HC_Tot_Temp_At_WBC_Target = 0.0

  End Subroutine Initialize_Death_Diag_Table

End Module HC_Init_DIV_Diag
