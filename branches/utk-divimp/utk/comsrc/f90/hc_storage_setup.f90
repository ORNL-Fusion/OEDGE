! -*-Mode:f90-*-
! Storage_setup.f90
! Storage Setup File
!
! Adam McLean under Prof. Peter Stangeby
! Research Assistant: David Elder
! January, 2002
!
! This module sets up each of the common blocks
! and assigns data for cell-specific information.
!

Module HC_Storage_Setup

  Use ComHC ! Access to setup and input values.

  Implicit None

  ! Keep type declarations in memory.
  Save

  !Type Cell_Geom_Table_Type

  ! Cell geometry dependent only.  Values that are the same throughout the cell.
  Real :: Local_K ! K
  Real :: Local_Ring_SMax ! SMAX

  !End Type Cell_Geom_Table_Type

  !Type State_Prop_Table_Type

  ! Contains all data that is hydrocarbon-state
  ! and/or position dependent.

  ! Cell properties dependent only.
  Real :: Local_Alphs ! ALPHS
  Real :: Local_Betas ! BETAS
  ! Cell position and HC Mass and/or charge dependent.		
  Real :: Local_HC_Change_State_Coll_Prob ! KFPS
  Real :: Local_HC_Tau_Parallel_Inv ! KKKFPS
  Real :: Local_HC_Tau_Stopping_Inv ! KFSS
  Real :: Local_HC_Tau_Heating_Inv ! KFTS
  Real :: Local_HC_Temp ! KTOLDS
  Real :: Local_Electric_Field ! KES
  Real :: Local_VPara ! VPARA
  Real :: Local_SPara ! SPARA
  Real :: Local_RConst ! RCONST
  Real :: Local_RGauss ! RGAUSS

  !End Type State_Prop_Table_Type

  !Type Global_Prop_Table_Type

  ! Contains all global properties which do not
  ! change in the course of hydrocarbon following.

  Integer :: Equate_Ion_Temp_Charge ! CIZDSRT
  Integer :: Impurity_Atomic_Num ! CION
  Integer :: INJ_Ring_Number ! INJIR
  Integer :: Max_Total_HC_Reflections ! MAXNRFCNT
  Integer :: Number_Injected_Particles ! INJNUM
  Integer :: Number_Time_Steps ! NTS
  Integer,parameter :: hc_NVel=nvel ! NVEL
  Integer :: Max_Ionization_States ! MAXIZS
  Integer :: Target_Material ! MATT
  Integer :: Plasma_Type ! MATP
  Integer :: Max_Time_Slices ! MAXNTS
  Real :: Ion_Time_Step ! QTIM
  Real :: Neutral_Time_Step ! FSRATE
  Real :: Z_Enh_Factor ! CZENH	
  Real :: Back_Plasma_Ion_Mass ! CRMB
  Real :: Back_Plasma_Charge ! RIZB
  Real :: Impurity_Ion_Mass ! CRMI
  Real :: Vel_Mult_Recomb_Neut ! CVRMULT
  Real :: Target_Binding_Energy ! CEBD
  Real :: Init_Particle_Temperature ! CTEM1
  Real :: Alt_Init_Particle_Temperature ! CTEM2
  Real :: Init_Electron_Temperature ! CTEB0
  Real :: Plate_Electron_Temperature ! CTEBP
  Real :: Max_HC_Ion_Iter_To_TMax ! CSTMAX
  Real :: Max_HC_Neut_Iter_To_TMax ! RSTMAX.
  Real :: CKK_Minimum ! CKKMIN
  Real :: CKK_Maximum ! CKKMAX
  Real :: Background_Drift_Velocity ! CDRFTV
  Real :: SOL_Drift_Start ! SDRFT_START
  Real :: SOL_Drift_End ! SDRFT_END
  Real :: SOL_Drift_Vel_Start ! CDRFTV_START
  Real :: SOL_Drift_Vel_End ! CDRFTV_END
  Real :: TGrad_Zero_Dist ! CSTGRAD
  Real :: EMax_Factor ! CEMAXF
  Real :: Self_ZEff ! CIZEFF
  Real :: Self_Sputter_Threshold ! CTRESH
  Real :: DThetaG ! DTHETG
  Real :: FPTimeO ! FPTIM0
  Real :: FPTimeI ! FPTIMI
  Real :: hc_FPXMaxO ! FPXMAX0
  Real :: hc_FPXMaxI ! FPXMAXI
  Real :: FP_Diff_Rate ! DIFFR
  Real,parameter :: hc_VelSep=velsep ! VELSEP
  Real :: hc_VelPlate ! VELPLATE
  Real :: INJ_Area_Lower_Bound ! INJF1
  Real :: INJ_Area_Upper_Bound ! INJF2
  Real :: Z_Charge ! CHZO
  Real :: ZO_Temp_Grad_Parameter ! CZO
  Real :: Cosine_Dist_Power !  CNIN
  Real :: Const_Velocity_Mult ! CVAMULT
  Real :: Saved_Velocity_Mult ! Velocity multiplier/particle used in HC_Launch_velocity
  Logical :: Ion_Diffuse ! DIFFUS

  ! Array data definitions.
  Real, Dimension (maxnks,maxnrs,Number_HC_Species) :: HC_MTCProb ! 		

  !End Type Global_Prop_Table_Type

  !Type Global_Geom_Table_Type

  ! Contains all geometry data relevent to hydrocarbon
  ! following that does not change in the course of
  ! the case.

  Integer :: Cell_TI ! IKTI
  Integer :: Cell_TO ! IKTO
  Integer :: Core_Ring ! IRCORE
  Integer :: Inner_SOL_Ring ! IRSEP
  Integer :: Inner_Target_Points ! NDSIN
  Integer :: Inner_Wall_Ring ! IRWALL
  Integer :: Max_Rings ! MAXNRS
  Integer :: Max_Cells_Per_Ring ! MAXNKS
  Integer :: Max_Points ! MAXPTS
  Integer :: Max_Target_Cells ! MAXNDS
  Integer :: Max_Velocity_Cells_Per_Ring ! MAXVNKS
  Integer :: Num_Boundary_Points ! PCNT
  Integer :: Num_Sectors ! ISECT
  Integer :: Num_Upper_Rings ! NRS
  Integer :: Num_Lower_Rings ! NRS2
  Integer :: Num_Wall_Points !  WALLPTS
  ! jdemod - needed to make the name different from the variable in comtor
  Integer :: hc_ReFCT ! REFCT
  ! jdemod - switch to using xpoint_up to define grid geometry from refct
  Logical :: hc_xpoint_up ! xpoint_up
  Integer :: Special_Ring ! IRSPEC
  Integer :: Outer_SOL_Ring ! IRSEP2
  Integer :: Outer_Wall_Ring ! IRWALL2
  Integer :: Upper_Trap_Ring ! IRTRAP
  Integer :: Lower_Trap_Ring ! IRTRAP2
  Integer :: First_Wall_Index ! WLWALL1
  Integer :: Last_Wall_Index ! WLWALL2
  Integer :: First_Trap_Index ! WLTRAP1
  Integer :: Last_Trap_Index ! WLTRAP2
  Integer :: Num_Target_Cells ! NDS
  Integer :: End_Target_Point ! NDSIN3
  Real :: R_Plasma_Centre ! R0
  Real :: Z_Plasma_Centre ! Z0
  Real :: R_X_Point ! RZP
  Real :: Z_X_Point ! ZXP
  Real :: Normal_Measure_To_X_Axis ! CSNORM
  Real :: Maximum_R ! RMAX
  Real :: Maximum_Z ! ZMAX
  Logical :: Grid_Error ! GRIDERR

  !End Type Global_Geom_Table_Type

  !Type DIVIMP_Options_Table_Type

  ! Contains all DIVIMP input options relevent to
  ! hydrocarbon transport.

  Integer :: Control_Switch ! CNEUTA
  Integer :: Impurity_Neutral_Vel_Opt ! CNEUTVEL
  Integer :: Collision_Opt ! CIOPTB
  Integer :: Friction_Opt ! CIOPTC
  Integer :: Heating_Opt ! CIOPTD
  Integer :: Injection_Opt ! CIOPTE
  Integer :: SOL_Opt ! CIOPTF
  Integer :: TEB_Coeff_Opt ! CIOPTM
  Integer :: TIB_Coeff_Opt ! CIOPTN
  Integer :: Reiser_Opt ! CIOPTR
  Integer :: Neutral_Reflection_Opt ! NRFOPT
  Integer :: Neutral_Mom_Coll_Opt ! MTCOPT
  Integer :: Far_Periphery_Opt ! FPOPT
  Integer :: Far_Periphery_Recycle_Opt ! FPROPT
  Integer :: Normal_Measure_Opt ! CNEUTE
  Integer :: Target_Position_Opt ! CTARGOPT
  Integer :: Non_Orthogonal_Grid_Opt ! NORTHOPT
  Integer :: Prompt_Deposition_Opt ! PROMPT_DEPOPT
  Integer :: Stop_Ion_In_Core_Opt ! CSTOP
  Integer :: First_Diffuse_Opt ! CDIFOP
  Integer :: Poloidal_Drift_Opt ! CPDRFT
  Integer :: RZ_Opt ! RZOPT
  Integer :: Neutral_Init_Pos_Opt ! INIT_POS_OPT
  Integer :: Self_Sputter_Opt ! CSELFS
  Integer :: Sputter_Opt ! CNEUTD
  Integer :: Target_Mirror_Opt ! CMIROPT
  Integer :: Ion_Wall_Opt ! CIONR
  Integer :: Grid_Opt ! CGRIDOPT

  !End Type DIVIMP_Options_Table_Type

  !Type DIVIMP_HC_Data_Table_Type

  ! Contains all DIVIMP input options relevent to
  ! hydrocarbon transport.

  Integer, Dimension (Number_Regions) :: HC_Region_Used ! Set if particle from a particular launch location was used (0-not used, 1-used).
  Character (Len=50), Dimension (Number_HC_IFates) :: Particle_Fate ! FATE, Character descriptions of each possible fate.
  Character (Len=10), Dimension (Number_Regions) :: HC_Region_Names ! Character description of launch regions.

  !End Type DIVIMP_HC_Data_Table_Type

  !Type Particle_Location_Table_Type

  Logical :: HC_Hit_Vessel ! Record the first time the particle hits the wall or target.

  ! Starting position indicators for initial launch.
  ! These do not change throughout each particle
  ! following.
  Real :: Launch_R ! ORGR
  Real :: Launch_Z ! ORGZ
  Real :: Launch_S ! 
  Real :: Launch_Cross ! 
  Integer :: Launch_HC_Species
  Integer :: Launch_Cell ! IKORG
  Integer :: Launch_Ring ! IRORG
  Integer :: First_Ioniz_Cell ! IKSTART in DIV.
  Integer :: First_Ioniz_Ring ! IRSTART in DIV.
  Integer :: Starting_Index ! IDSTART
  Integer :: Launch_Wall_Index ! IWSTART	
  Integer :: Launch_Type ! IDTYPE in DIV, NEUTTYPE in NEUT
  ! IDTYPE, defining the type of the launch - Wall or Target - Physical or Chemically sputtered.
  ! Targ+phys = 1, Targ+chem = 2, Targ-self = 3, Wall+phys = 4, Wall+chem = 5, 2Dneutral = 6, Refl Ion  = 7

  ! Position indicators specific to reflection/sputtering.
  ! These are initialized at 0 and are updated at each
  ! successful reflection/sputtering event.
  Real :: Relaunch_R ! Re-ORGR
  Real :: Relaunch_Z ! Re-ORGZ
  Real :: Relaunch_S ! 
  Real :: Relaunch_Cross ! 
  Integer :: Relaunch_HC_Species
  Integer :: Relaunch_Cell ! Re-IKORG
  Integer :: Relaunch_Ring ! Re-IRORG
  Integer :: Relaunch_Target_Index ! Re-IDSTART
  Integer :: Relaunch_Wall_Index ! Re-IWSTART
  Integer :: Relaunch_Type ! Re-IDTYPE

  !End Type Particle_Location_Table_Type

  ! Begin diagnostic grouping definitions applicable to all launches.
  ! Contains all internal HC specific counters with data that can be
  ! sent to OUT file, and/or printed with case summary.  Input
  ! from the DIVIMP case file for HC options is stored in the ComHC 
  ! common block.

  !Type HC_Launch_Diag_Table_Type

  ! Default output recorded for OUT processing.
  Character(Len=10), Dimension (Number_HC_Species) :: HC_State_List ! List of character names of each hydrocarbon species used in OUT processing.
  Integer, Dimension (Number_HC_Species + 1) :: HC_Output_List ! Determines which hydrocarbons to include in stored binary file for OUT processing: 0=not included, 1=included.

  ! Contains all data stored at particle launch.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Total_R_Prod_Positions ! XTOT
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Total_Z_Prod_Positions ! YTOT
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Total_S_Prod_Positions ! 
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Total_Cross_Prod_Positions ! 
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Total_Prod_Angles ! ATOT
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Sum_Fragments_Launched ! RNEUT
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Total_Prod_Vels_No_VMult
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Max_Prod_Vel_No_VMult
  !		Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Total_Vel_Ang_Mults
  !		Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Max_Vel_Ang_Mults
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Total_Prod_Velocities ! VTOT
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Max_Prod_Velocities ! VTOTM
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Prod_Temperatures ! ETOT
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Num_Failed_Launches ! RFAIL
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Num_Vel_Greater_Max_Ran ! REJECT
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Launch_Grid_Error_Moved_Okay ! ERR_OK, Begins just outside the wall boundary, but one step moved it in okay.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Launch_Grid_Error_Moved_Out ! ERR_OUT, Begins just outside the wall boundary and moving it does not help.	

  Real, Dimension (Max_Impurities,4) :: HC_LaunchDat ! maximp found in comhc.

  !End Type HC_Launch_Diag_Table_Type

  !Type HC_Transport_Diag_Table_Type

  ! Contains all data stored at each timestep during transport.

  ! Timing data.
  Integer :: HC_Walk_Count ! IW, Used for storing walk position data.
  Double Precision, Dimension (Number_Regions) :: HC_Neutral_Timesteps_Count ! Number of timesteps spent as a neutral.
  Double Precision, Dimension (Number_Regions) :: HC_Ion_Timesteps_Count ! Number of timesteps spent as an ion.
  Double Precision :: Time_All_HCs (Number_Regions) ! CISTOT, Total time steps spent following all HCs launched.
  Double Precision :: Max_Time_Any_HC (Number_Regions) ! CISMAX, Maximum time steps for any single HC following.
  Double Precision, Dimension (Number_HC_Species,Number_Regions) :: HC_State_Time ! CIEIZS, Time spent in each state from each target.
  Double Precision, Dimension (Number_Regions) :: HC_Max_Ion_Teq_Iter ! Maximum ion eqivalent timesteps through lifetime.
  Double Precision, Dimension (Number_Regions) :: HC_Total_Ion_Teq_Iter ! Total ion eqivalent timesteps through lifetime.
  Double Precision, Dimension (Number_HC_Species,Number_Regions) :: HC_Max_Time_In_State ! Maximum time spent in each state from each target.
  Double Precision, Dimension (Number_HC_Species,Number_Regions) :: HC_State_TimeSteps ! Timesteps spent in each state from each target.
  Double Precision, Dimension (Number_HC_Species,Number_Regions) :: HC_Max_TimeSteps ! Maximum number of timesteps spent in each state from each target.
  Double Precision :: HC_Tot_Time_First_Diff (Number_Regions) ! RDIFFT, Total time to first diffusion.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Num_Reach_Max_Iter ! CICUTS, Particles reaching time cutoff.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Exist_Max_Iter ! TCUT, HCs existing at TMax (CSTMAX).

  ! State sums.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Num_Orig_Reach_State ! CITIZS, Number of original particles that reach each state.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Reach_State ! Total number of times the particle reaches each state (fraction of particles may be >1.0 due to reflections at higher HC state).

  ! Movement sums.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Num_Striking_Target ! RSTRUK
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Num_Enter_Main_Plasma ! RMAIN
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Fragments_Exit_Main ! REXIT
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Num_Entered_Core ! NUM_ENTERED_CORE, Number of ions entering core plasma.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Num_Reach_Centre ! RCENT
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Stopped_Follow_Ion_In_Core ! STOPPED_FOLLOW, when ion within inner SOL ring and option CSTOP is set.

  ! Position indicators.
  Real, Dimension (Max_Number_Walks,2) :: HC_Walks ! Follows position of particles as they travel for each timestep up to maxnws steps - typically 10,000.
  Logical :: HC_InMain ! INMAIN, Indicates if the particle begin in the main plasma.
  Logical :: HC_InCore ! INCORE, Indicates if the particle has entered the core plasma in its lifetime.
  Logical :: HC_InEdge ! INEDGE, Indicates if the particle has entered the edge region in its lifetime.
  Logical :: HC_InMSOL ! INMSOL, Indicates if the particle has entered the main SOL region in its lifetime.
  Logical :: HC_InDIV ! INDIV, Indicates if the particle has entered the divertor region in its lifetime.
  Logical :: HC_InTrap ! INTRAP, Indicates if the particle has entered the trap region in its lifetime.
  Logical :: HC_CflRxa ! CFLRXA, Indicates particle has not reflected off central mirror at least once (particle is not counted as re-reflecting if it does so more than once).
  Logical :: HC_CflRex ! CFLREX, Indicates particle has not exited main plasma at least once (particle is not counted as re-exiting if it crosses more than once).
  Logical :: HC_CflRin ! CFLRIN, Indicates particle has not entered main plasma at least once (particle is not counted as re-entering if it crosses more than once).
  ! Note:  Entering statistics are not recorded if particle begins within the main plasma.

  ! Distance sums.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Cumu_CroFly_Dist_Trav_By_Pos ! Cumulative distance travelled from birth as the crow flies calculated by r,z position.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Cumu_Real_Dist_Trav_By_Pos ! Cumulative total distance travelled as calculated by r,z position.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Cumu_Real_Dist_Trav_By_VT ! Cumulative total distance travelled as calculated by velocity * timestep.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_State_Dist_Trav_By_Pos ! State distance travelled as calculated by r,z position.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_State_Dist_Trav_By_VT ! State distance travelled as calculated by velocity * timestep.

  ! Density recording by current position.
  Double Precision, Dimension (maxnks,maxnrs,Number_HC_Species) :: HC_Density ! Hydrocarbon density array. Records cell position and species of each fragment after each timestep.
  Double Precision, Dimension (maxnks,maxnrs,Num_H_States) :: H_Density ! Hydrogen isotope density array. Records cell position and species of each H released from a breakup at the point it is created.
  Real, Dimension (maxnks,maxnrs,Number_HC_Species) :: HC_ChemDen ! CHEMDEN
  Real, Dimension (maxnks,maxnrs,Number_HC_Species) :: HC_ChemIZS ! CHEMIZS
  Real, Dimension (3,Number_HC_Species,Number_Regions) :: HC_DDVoid ! DDVOID
  Real, Dimension (maxnks,3,-1:Number_HC_Species) :: HC_ELims ! ELIMS
  Real, Dimension (maxnks,maxnrs,-1:Number_HC_Species,maxnts) :: HC_Lims ! LIMS
  Double Precision, Dimension (Number_HC_Species,6,Number_Regions) :: HC_DParas ! DPARAS

  ! C and CH/CD neutral data for spectroscopic analysis.
  Real, Dimension (maxnks,maxnrs,-1:maxizs) :: HC_TIZS ! TIZS, Should contain only C (0) and C+ (1) data.
  Real, Dimension (maxnks,maxnrs,-1:0) :: HC_TIZS_CH ! TIZS, For CH/CD band line emission.
  !		Real, Dimension (maxnks,maxnrs,-1:0) :: HC_TIZS_C2 ! TIZS, For C2 band line emission, useful for higher hydrocarbons.
  ! jdemod - corrected use of hc_ddts and eliminated unused/misused variables
  !Double Precision, Dimension (maxnks,maxnrs,0:maxizs) :: HC_DDLims ! DDLIMS, Should contain only C (0) and C+ (1) data.
  !Double Precision, Dimension (maxnks,maxnrs,0:maxizs) :: HC_DDLims_CH ! DDLIMS, For CH/CD band line emission.
  !Double Precision, Dimension (maxnks,maxnrs,0:maxizs) :: HC_DDts ! DDTS, Should contain only C (0) and C+ (1) data.
  Double Precision, Dimension (maxnks,maxnrs,number_hc_species) :: HC_DDts ! DDTS - should contain temperatures for all HC species
  !Double Precision, Dimension (maxnks,maxnrs,0:maxizs) :: HC_DDts_CH ! DDTS, For CH/CD band line emission.

  ! Density recording, recorded by starting ik,ir.
  Real, Dimension (maxnks,maxnrs) :: HC_Ion_Core_Density ! NCORE
  Real, Dimension (maxnks,maxnrs) :: HC_Ion_Edge_Density ! NEDGE
  Real, Dimension (maxnks,maxnrs) :: HC_Ion_Trap_Density ! NTRAP
  Real, Dimension (maxnks,maxnrs) :: HC_Ion_Divertor_Density ! NDIVERT
  Real, Dimension (maxnks,maxnrs) :: HC_Ion_MSOL_Density ! NMSOL

  ! MTC, Momentum Transfer Collisions.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_MTC_Striking_Target ! MTC_RSTRUK
  Real, Dimension (7,3) :: HC_MTCInf ! MTCInf, Records momentum transfer collision information.
  Real, Dimension (0:11,3) :: HC_MTCTotCnt ! MTCTOTCNT

  ! Leakage statistics.
  Integer :: Number_Leaked_Core (Number_Regions) ! NLEAKCORE
  Real :: Total_Leaked_Core (Number_Regions) ! TOTLEAKCORE
  Logical :: HC_Has_Leaked ! HASLEAKED
  Logical :: HC_Has_Leaked_Core ! HASLEAKEDCORE
  Real, Dimension (Max_Impurities,Number_Regions,2) :: HC_Leak_Position ! CLEAKPOS
  Real :: HC_Leak_Time (Number_Regions) ! CLEAKT
  Real, Dimension (maxpts,Number_HC_Species+1) :: HC_Leak_Density ! CLEAKN
  Integer :: HC_Leak_Particles (Number_Regions) ! CLEAKP

  ! Far-periphery counting.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_FPTarg ! FPTARG
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_FPTart ! FPTART, Number lost to far periphery target.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_RFPTarg ! RFPTARG, Number of HC ions plating on FP target.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_FPEnt ! FPENT, Number entering far periphery.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_FPExit ! FPEXIT, Number exiting far periphery.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_FPTTotal ! FPTTOT, Number existing to TMax in far periphery.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Core_From_FP_Ref ! CVVFPREF, Number of particles entering main from far-periphery with reflections.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Core_From_FP_NoRef ! CVVFPNRF, Number of particles entering main from far-periphery without reflections.

  ! Counters used during ion transport.
  Double Precision, Dimension (Number_HC_Species, 9) :: HC_CoreOuts ! Force and time data for ion transport.
  Double Precision, Dimension (6) :: HC_DSParaNorm ! DSPARANORM
  Double Precision, Dimension (6, Number_HC_Species + 1) :: HC_DVParaNorm ! DVPARANORM
  Double Precision, Dimension (6) :: HC_DSParaStep ! DSPARASTEP
  Double Precision, Dimension (6, Number_HC_Species + 1) :: HC_DVParaStep ! DVPARASTEP
  Double Precision, Dimension (6, Number_HC_Species + 1) :: HC_VParaStep ! VPARASTEP
  Double Precision, Dimension (6, Number_HC_Species + 1) :: HC_VParaNorm ! VPARANORM
  Double Precision, Dimension (6) :: HC_DSParaCnt ! DSPARACNT
  Double Precision, Dimension (6, Number_HC_Species + 1) :: HC_DVParaCnt ! DVPARACNT
  Double Precision, Dimension (6, Number_HC_Species + 1) :: HC_DVMinV ! DVMINV
  Double Precision, Dimension (6, Number_HC_Species + 1) :: HC_DVMaxV ! DVMAXV

  ! Force counters for ion transport.
  Real, Dimension (maxnks,maxnrs,maxizs) :: HC_FCell ! FCELL
  Real, Dimension (maxnks,maxnrs,maxizs) :: HC_Ffi ! FFI
  Real, Dimension (maxnks,maxnrs,maxizs) :: HC_Fthi ! FTHI
  Real, Dimension (maxnks,maxnrs,maxizs) :: HC_Fvbg ! FVBG
  Real, Dimension (maxnks,maxnrs,maxizs) :: HC_Force_Diff ! DIFF
  Real, Dimension (maxnks,maxnrs,maxizs) :: HC_IonVelAvg ! VELAVG

  ! Velocity statistics.
  ! Note, we expect to trave only one charge state.
  Real, Dimension (maxnks,maxnrs) :: HC_SDVS ! SVDS
  Real, Dimension (maxnks,maxnrs) :: HC_SDVS2 ! SVDS2
  Real, Dimension (maxnks,maxnrs,2) :: HC_SDVS3 ! SVDS3
  Real, Dimension (-nvel:nvel+1,maxvnks) :: HC_VelSpace ! VELSPACE
  Real, Dimension (-nvel:nvel+1,maxvnks) :: HC_VelWeight ! VELWEIGHT

  ! Temperature statistics.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Temp_Reach_Max_Iter ! CRTRCS, Temperature of particles reaching cutoff time.

  ! Movement indicators.
  Real, Dimension (-1:Number_HC_Species,Number_Regions) :: HC_Num_Orig_Enter_Main ! CLLL, Number of original particles entering main plasma (particles not counted >1 if enter main plasma >1).
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Z_Orig_Enter_Main ! CLLLX, Total Z value of all original main plasma entering events.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_S_Orig_Enter_Main ! CLLLS, Total MIN(S,S-SMAX) value of all original main plasma entering events.
  Real, Dimension (-1:Number_HC_Species,Number_Regions) :: HC_Num_Orig_Enter_Core ! CNNN, Number of original particles entering core plasma.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Z_Orig_Enter_Core ! CNNNX, Total Z value of all original core entering events.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_S_Orig_Enter_Core ! CNNNS, Total MIN(S,S-SMAX) value of all original core entering events.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Teq_Orig_Enter_Core ! CNNNT, Total ion/neutral equivalent timesteps at time of original core entry event.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Min_Z_Reach ! CXXX, Total minimum Z value reached.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Max_S_Reach ! CSSS, Total maximum S or SMAX-S value reached.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Max_S_Reach_Ion_Removal ! CSSSS, Total maximum S or SMAX-S value reached for removed ions.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Core_From_Reg_Ref ! CVVREFM, Number of particles entering main from regular launch with reflections.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Core_From_Reg_NoRef ! CVVNRFM, Number of particles entering main from regular launch without reflections.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_R_Species_Creation ! CVVXC, Total R position at species creation.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Z_At_Main_Entry ! CVVXE, Total Z position at particle entry into main plasma.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Num_Orig_Central_Reflect ! CICRXA, Total original particles reflecting off central mirror.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Teq_Orig_Central_Reflect ! CISRXA, Total ion-equivalent timesteps for original particles to reflect off central mirror.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Temp_At_Central_Reflect ! CITRXA, Total temperatures when original particles reflect off central mirror.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Min_Teq_Orig_Central_Reflect ! CIFRXA, Minimum ion-equivalent timesteps for reflection off central mirror.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Ion_Collision ! CICCOL, Total number of ion collisions.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Start_Main ! CICRNJ, Number of particles started in main plasma.
  ! Note: CIFRIN and CISRIN are recorded in REAL, not double presision - may be a problem if timesteps go over 16 million.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Min_Teq_Orig_Enter_Main ! CIFRIN, Minimum number of ion-equivalent timesteps to original particles entering main plasma.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Teq_Orig_Enter_Main ! CISRIN, Total ion-equivalent timesteps to original particles entering main plasma.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Temp_Orig_Enter_Main ! CITRIN, Total temperature of original particles entering main plasma.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Num_Ions_Lost ! CICLOS, Number of ions lost.
  ! Note: CIFLOS, CILLOS, CISLOS are recorded in REAL, not double presision - may be a problem if timesteps go over 16 million.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Min_Teq_Ions_Lost ! CIFLOS, Minimum time to ion loss.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Max_Teq_Ions_Lost ! CILLOS, Maximum time to ion loss.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Teq_Ions_Lost ! CISLOS, Total time to ion loss.

  !		Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Orig_Enter_Main ! CICRIN, Total original particles entering main plasma.
  !		Real, Dimension (-1:Number_HC_Species,Number_Regions) :: HC_Cmmm ! CMMM, Note: No distinction between originally ionized in main or SOL is required for HCs.
  !		Real, Dimension (Number_HC_Species,Number_Regions) :: HC_CmmmX ! CMMMX, 
  !		Real, Dimension (Number_HC_Species,Number_Regions) :: HC_CmmmS ! CMMMS, 
  !		Real, Dimension (Number_HC_Species,Number_Regions) :: HC_CNorgS ! CNORGS, 
  !		Real, Dimension (Number_HC_Species,Number_Regions) :: HC_CNorgR ! CNORGR, 
  !		Real, Dimension (Number_HC_Species,Number_Regions) :: HC_CNorgZ ! CNORGZ, 
  !		Real, Dimension (Number_HC_Species,Number_Regions) :: HC_CvvRM ! CVVRM, 
  !		Real, Dimension (Number_HC_Species,Number_Regions) :: HC_CvvZM ! CVVZM, 
  !		Real, Dimension (Number_HC_Species,Number_Regions) :: HC_CvvSM ! CVVSM, 

  !End Type HC_Transport_Diag_Table_Type

  !Type HC_Evolve_Diag_Table_Type

  ! Contains all data stored during state transitions.
  Real, Dimension (Number_HC_Species,Number_HC_Reactions,Number_Regions) :: HC_Reaction_Count ! Count for each available transition for each species.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_R_Pos_At_Prod ! XATIZ
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Z_Pos_At_Prod ! YATIZ
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_S_Pos_At_Prod ! 
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Cross_Pos_At_Prod ! 
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Angle_At_Prod ! ATOT
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Num_Fragments ! RNEUT
  Real :: HC_Num_Fragments_Reach_CIon (Number_Regions) ! Number of particles that reach C+.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Velocity_At_Prod ! VATIZ
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Max_Velocity_At_Prod ! VATIZM
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Temperature_At_Prod ! EATIZ
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Time_To_Prod ! TATIZ

  Real, Dimension (Max_Impurities,Number_HC_Species) :: HC_Time_At_Production
  Real, Dimension (Max_Impurities,Number_HC_Species) :: HC_Energy_At_Production
  Real, Dimension (Max_Impurities,Number_HC_Species) :: HC_Kin_E_Add_At_Production

  Real, Dimension (2,2,2,2,5) :: HC_LIonizDat ! lionizdat

  ! jdemod - removed 0 and -1 indices since they don't apply to the HC code
  ! Calculated for out.
  Real, Dimension (Number_HC_Species) :: HC_Factor_A ! Species dependent normalization per unit particle density.
  Real, Dimension (Number_HC_Species) :: HC_Factor_B ! Species dependent normalization per unit particle density per unit time (ion or neutral timestep).
  !Real, Dimension (-1:Number_HC_Species) :: HC_Factor_A ! Species dependent normalization per unit particle density.
  !Real, Dimension (-1:Number_HC_Species) :: HC_Factor_B ! Species dependent normalization per unit particle density per unit time (ion or neutral timestep).

  !End Type HC_Evolve_Diag_Table_Type

  !Type HC_Rerelease_Diag_Table_Type

  ! Contains all data stored during reflection/sputtering.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Total_R_ReProd_Positions ! XTOT
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Total_Z_ReProd_Positions ! YTOT
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Total_ReProd_Angles ! ATOT
  Real, Dimension (Number_HC_Species+1,maxpts+1,Number_Regions) :: HC_Sum_Fragments_ReLaunched ! RNEUT
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Total_ReProd_Vels_No_VMult
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Max_ReProd_Vel_No_VMult
  !		Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Total_Vel_Ang_Mults
  !		Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Max_Vel_Ang_Mults
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Total_ReProd_Velocities ! VTOT
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Max_ReProd_Velocities ! VTOTM
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_ReProd_Temperatures ! ETOT
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_ReLau_Grid_Error_Moved_Okay ! ERR_OK, Begins just outside the wall boundary, but one step moved it in okay.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_ReLau_Grid_Error_Moved_Out ! ERR_OUT, Begins just outside the wall boundary and moving it does not help.	

  Integer, Dimension (Number_HC_Species,Number_Regions) :: Total_HC_Reflections ! TOTRF, Total reflections following all HC derivatives.
  Integer, Dimension (Number_HC_Species,Number_Regions) :: Max_HC_Reflections_Found ! MAXRF, Maximum reflections per particle as calculated while following particles.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Reflection_Loss ! NRFLOSS, Number of HC's not counted because we went over the maximum allowable reflection limit.

  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_YldTot ! YLDTOT, Total self-sputtering yield.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_YThTot ! YTHTOT, Total number of yields above threshold yield.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_YldMax ! YLDMAX, Maximum yield for any neutral HC fragment.

  !End Type HC_Rerelease_Diag_Table_Type

  !Type HC_Vessel_Int_Diag_Table_Type

  Real, Dimension (maxnds,Number_HC_Species) :: HC_Deposit ! DEPS, Deposited particles. ! Note maxnds is use associated in ComHC by params.
  Real, Dimension (maxnds,5,Number_HC_Species) :: HC_Erosion ! NEROS, Number of particles eroded from targets.
  Real, Dimension (maxnds,6,Number_HC_Species) :: HC_PromptDeps ! PROMPTDEPS, Prompt deposited particles.
  Real, Dimension (maxnks,maxnrs,0:Number_HC_Species+1) :: HC_Walls ! WALLS
  Real, Dimension (maxpts+1,Number_HC_Species) :: HC_WallsE ! WALLSE, Wall erosion.
  Real, Dimension (maxpts+1,Number_HC_Species) :: HC_WallsE_I ! WALLSE_I, Wall erosion, ionized particles.
  Real, Dimension (maxpts+1,Number_HC_Species) :: HC_WallsI ! WALLSI, Wall erosion, ionized particles.
  Real, Dimension (maxpts+1,Number_HC_Species) :: HC_WallsN ! WALLSN, Wall erosion, neutral particles.
  Real, Dimension (maxpts,maxpts+1,3) :: HC_WTDep ! WTDEP
  Real, Dimension (maxpts,maxnrs,4,6) :: HC_WTSource ! WTSOURCE, Source origins.

  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_RWall ! RWALL, Number of HCs plating out on walls (compare to RWALLN which includes neutrals).
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_RDep ! RDEP, Number of HCs plating out on targets.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Num_Reach_Wall ! RWALLN, Total HC fragments plating out on walls.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_MTC_Reach_Wall ! MTC_RWALLN			

  !End Type HC_Vessel_Int_Diag_Table_Type

  !Type HC_Death_Diag_Table_Type

  ! Contains all data stored at completion of particle following.
  Integer, Dimension (Number_HC_Species,Number_Regions,Number_HC_IFates) :: HC_IFate_Count ! Record particle IFate data.
  Integer, Dimension (maxpts+1,Number_HC_Species+1) :: HC_Wall_Deposit_Count ! Records deposit events for each wall, each HC species.
  Integer, Dimension (maxpts+1) :: HC_Carbon_Wall_Deposit_Count ! Total number of carbon atoms depositing on each wall segment.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Num_At_TMax ! RTMAX

  ! Counters for end of lifetime.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Num_Absorbed_TargWall ! CICABS, Number of fragments absorbed on targets and walls.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Min_Teq_To_Absorption ! CIFABS, Ion-eqivalent timesteps to first particle absorbed.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Max_Teq_To_Absorption ! CILABS, Ion-eqivalent timesteps to last fragment absorbed.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Teq_To_Absorption ! CISABS, Total time to fragment absorption.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Temp_At_Absorption ! CRTABS, Total temperature at fragment absorption.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Velocity_At_Absorption ! CRVABS, Total velocity at fragment absorption.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_ABS_Vel_At_Absorption ! CRAVAV, Total ABS(Velocity at fragment absorption.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Elec_Temp_At_Absorption ! CTBS, Total electron temperature at fragment absorption.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Num_Absorbed_Act_Target ! ACTTARG, Number of fragments absorbed on actual target.
  Real, Dimension (10,Number_Regions) :: HC_Ctexs ! CTEXS, Total temperature for each charge state at fragment absorption.

  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Temp_At_WBC_HC_Target ! Storage for target/wall collision temperature in HC code.
  Real, Dimension (Number_HC_Species,Number_Regions) :: HC_Tot_Temp_At_WBC_HC_Boundary ! Storage for freespace boundary collision temperature in HC code.
  Real, Dimension (0:MAXIZS) :: HC_Tot_At_WBC_Boundary ! Storage for freespace boundary collision count for comparison with WBC.
  Real, Dimension (0:MAXIZS) :: HC_Tot_Temp_At_WBC_Boundary ! Storage for freespace boundary collision for comparison with WBC.
  Real, Dimension (0:MAXIZS) :: HC_Tot_At_WBC_Target ! Storage for target collision count for comparison with WBC.
  Real, Dimension (0:MAXIZS) :: HC_Tot_Temp_At_WBC_Target ! Storage for freespace boundary collision for comparison with WBC.

  !End Type HC_Death_Diag_Table_Type

  !Type Misc_Data_Table_Type

  ! Contains all global miscellaneous data.
  Real :: Calc_Hi ! 
  Real :: Calc_Lo ! 
  Real :: Pi_Value !
  Real :: Root_2 ! 
  !		Real :: AMU ! 
  !		Real :: ECH ! 
  !		Real :: KBoltz ! 
  Real :: Rad_In_A_Deg ! 
  Logical :: Debug_HC_Ion ! DEBUGL
  Logical :: Debug_HC_Neutral ! DEBUGN
  Logical :: Debug_HC_V ! DEBUGV
  Logical :: Debug_HC_Prompt ! DEBUG_PROMPT
  Real :: Print_Debug_x_CStep_Ion ! CSTEPL
  Real :: Print_Debug_x_CStep_Neutral ! CSTEPN
  Integer :: Impurity_Limit ! IMPLIM
  Real :: CPU_Time_Limit ! CPULIM
  !End Type Misc_Data_Table_Type


End Module HC_Storage_Setup
