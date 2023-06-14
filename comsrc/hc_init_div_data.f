! -*-Fortran-*-
! Initialize_Cell_Properties.f90
! Cell Properties Initialization File
!
! Adam McLean under Prof. Peter Stangeby
! Research Assistant: David Elder
! January, 2002
!
! This module initializes each of the required data structures for
! all cell properties.  It acts as the data-passing front-end between
! original DIVIMP code and all the added transport and hydrocarbon
! following routines added.

      Module HC_Init_DIV_Data

	! Use statements.
	Use HC_Storage_Setup ! Access to derived data types.	

	! Every good Fortran 90 program has...
	Implicit None

	! Define data tables.
	!Type (Cell_Geom_Table_Type) :: Cell_Geom_Table
	!Type (State_Prop_Table_Type) :: State_Prop_Table
	!Type (Global_Prop_Table_Type) :: Global_Prop_Table
	!Type (Global_Geom_Table_Type) :: Global_Geom_Table
	!Type (DIVIMP_Options_Table_Type) :: DIVIMP_Options_Table
	!Type (DIVIMP_HC_Data_Table_Type) :: DIVIMP_HC_Data_Table
	!Type (Particle_Location_Table_Type) :: Particle_Location_Table
	!Type (Misc_Data_Table_Type) :: Misc_Data_Table

      Contains

	Subroutine Initialize_Cell_Geom_Data (ik,ir)
	
		! Required modules.
		Use HC_Get ! gkks, gksmaxs.
		
		! Every good Fortran 90 program has...
		Implicit None
		
		! Declare subroutine variables.
		Integer, Intent (In) :: ik ! Cell.
		Integer, Intent (In) :: ir ! Ring.

		! Geometric properties.
		! Cell position dependent only.
		 Local_K = gkks (ir)
		 Local_Ring_SMax = gksmaxs (ir)

		! Update drfit values with new SMAX.
		 SOL_Drift_Start = 
     >             SOL_Drift_Vel_Start *  
     >            Local_Ring_SMax
		 SOL_Drift_End = 
     >             SOL_Drift_Vel_End *  
     >            Local_Ring_SMax
		 CKK_Maximum = MAX (
     >             CKK_Maximum,  Local_K)

	End Subroutine Initialize_Cell_Geom_Data

	Subroutine Initialize_State_Prop_Data (ik,ir,iz)
			
		! Required modules.
		Use HC_Get ! Contains external functions to interact with DIVIMP.

		! Every good Fortran 90 program has...
		Implicit None
		
		! Declare subroutine variables.
		Integer, Intent(In) :: ik ! Cell.
		Integer, Intent(In) :: ir ! Ring.
		Integer, Intent(In) :: iz ! Charge.
		
		If (iz .gt. 0) Then
		   Local_Alphs = gkalphs (iz) ! Used to find FEG (electron temperature gradient force).
		   Local_Betas = gkbetas (iz) ! Used to find FIG (ion temperature gradient force).
!write (0,*) "Updating state prop ",ik,ir,iz
		
                  ! Cell position and HC Mass and/or charge dependent.  Note that these are only correct after modify_taus.f90 has been run.
		   Local_HC_Change_State_Coll_Prob = 
     >              MAX ( Calc_Lo, gkfps (ik,ir,iz)) ! LFPS, Change of state collision prob 0<= kfps <= 1
		   Local_HC_Tau_Parallel_Inv = 
     >              gkkkfps (ik,ir,iz) ! 
		   Local_HC_Tau_Stopping_Inv = 
     >              MAX ( Calc_Lo, gkfss (ik,ir,iz)) ! LFSS, Local friction force susceptablility.
		   Local_HC_Tau_Heating_Inv = 
     >              MAX ( Calc_Lo, gkfts (ik,ir,iz)) ! LFTS
		   Local_HC_Temp = gctemav() ! LTOLDS
		
		  ! Note that the e-field may be updated at each timestep with the enhanced sheath model of Brooks.
		   Local_Electric_Field = gkes (ik,ir) ! KES
                End If

		! Setting of spara and vpara currently correspond to collision option 13.
		 Local_SPara = 0.0 ! 
		 Local_VPara = 0.0 ! 
		 Local_RConst = 0.0 ! 
		 Local_RGauss = 0.0 ! 

	End Subroutine Initialize_State_Prop_Data

	Subroutine Initialize_Global_Prop_Data ()

		! Every good Fortran 90 program has...
      use mod_params
      use mod_comtor
      use mod_dynam4
      use mod_commv
      use mod_promptdep
      use mod_reiser_com
      use mod_cgeom
      use mod_fperiph_com
      use mod_diagvel
      use mod_driftvel
		Implicit None

		! Included common blocks.
c	Include 'params' ! Required by comtor
c	Include 'comtor' ! Contains czenh,crmb,rizb,qtim,fsrate,crmi,cion,cneutvel,mtcopt,fpropt(COMTOR2),ctargopt,northopt,cdrftv,cdrftv_start,cdrftv_end,cteb0,ctebp,czo,chzo
c	Include 'dynam4' ! Contains nts,cstmax
c	Include 'commv' ! Contains ckkmin,ckkmax
c	Include 'promptdep' ! Contains prompt_depopt
c	Include 'reiser_com' ! Contains coptr
c	Include 'cgeom' ! Contains ikti,ikto,dthetg
c	Include 'fperiph_com' ! Contains fptimi,fpxmaxi
		!Include 'slcom' ! Contains fptimi,fpxmaxi
c	Include 'diagvel' ! Contains nvel,velplate,velsep
                !
                ! jdemod           
                !
c               include 'driftvel'
                
		! Global properties.
		 Ion_Time_Step = qtim ! Ion time step (s).
		 Neutral_Time_Step = fsrate ! Neutral time step (s).
		 Number_Time_Steps = nts ! Number of time steps (typically 1).
		 Z_Enh_Factor = czenh ! Charge enhancement factor.
		 Back_Plasma_Ion_Mass = crmb ! Background plasma ion mass (2=D, 2.5=DT).
		 Back_Plasma_Charge = rizb ! Background plasma ion charge (normally 1).
		 Impurity_Ion_Mass = crmi ! Impurity ion mass (12.0=C).
		 Impurity_Atomic_Num = cion ! Impurity atomic number (6=C).
		 Vel_Mult_Recomb_Neut = cvrmult ! Velocity multiplier for recombined neutrals (normally 1.0).
		 Max_Ionization_States = maxizs ! Maximum number of impurity charge states handled by DIVIMP (normally 6).
		 Max_Time_Slices = maxnts ! Number of time frames to model in time (normally 1).
		 Target_Binding_Energy = cebd ! Target material surface binding energy (eV).
		 Init_Particle_Temperature = ctem1 ! Initial ion, neutral, gas, etc. temperature (replaces Tin, Ein, Tg and converted to velocity for case of ion injection)
		 Alt_Init_Particle_Temperature = 
     >            ctem2 ! Alternate inital particle temperature used for Vel/Ang flag 9.
		 Init_Electron_Temperature = cteb0 ! Temperature of electrons at t=0.
		 Plate_Electron_Temperature = ctebp ! Temperature of electrons at plates.
		 Max_Total_HC_Reflections = 500 ! MAXNRFCNT, Maximum number of HC reflections off the wall or targets.
		 Max_HC_Ion_Iter_To_TMax = 
     >            cstmax ! CSTMAX, Maximum timesteps before we conclude no reaction will occur (typically 1E7).
!     >            10.0 /  Ion_Time_Step ! CSTMAX, Maximum timesteps before we conclude no reaction will occur (typically 1E7).
!     >          10000 ! CSTMAX, Maximum timesteps before we conclude no reaction will occur.
		 Max_HC_Neut_Iter_To_TMax = 
     >            0.1 /  Neutral_Time_Step ! RSTMAX, Local to NEUT/NEUTONE.
		 CKK_Minimum = ckkmin ! Used for update_cross in ion transport.
		 CKK_Maximum = ckkmax ! Used for update_cross in ion transport.
		 Background_Drift_Velocity = cdrftv ! Poloidal drift velocity in the SOL.
		 SOL_Drift_Vel_Start = cdrftv_start ! Starting velocity for SOL poloidal drift.
		 SOL_Drift_Vel_End = cdrftv_end ! Ending velocity for SOL poloidal drift.
		 SOL_Drift_Start = 0.0 ! sdrft_start, Starting position for SOL poloidal drift.
		 SOL_Drift_End = 0.0 ! sdrft_end, Ending position for SOL poloidal drift.
		 TGrad_Zero_Dist = cstgrad ! Distance where Tgrad forces become 0.
		 EMax_Factor = cemaxf ! EMax-factor.
		 Ion_Diffuse = .False. ! DIFFUS, Ion diffusion selector.
		 Self_ZEff = cizeff ! Self Z effective.
		 Self_Sputter_Threshold = ctresh ! Self-sputtering threshold (eV).
		 Equate_Ion_Temp_Charge = cizset ! Charge at which Ti is set to Tb.
		 Cosine_Dist_Power = cnin ! Cosine distribution power from DIVIMP input file.
		 Const_Velocity_Mult = cvamult ! Velocity multiplier used for Vel/Ang flags 14 and 15.
		
		 DThetaG = dthetg ! DTHETG
		 FPTimeO = fptimO ! FPTIMO
		 FPTimeI = fptimI ! FPTIMI
		 hc_FPXMaxO = fpxmaxO ! FPXMAXO
		 hc_FPXMaxI = fpxmaxI ! FPXMAXI
		 FP_Diff_Rate = SQRT (2.0 * cdperpfp * qtim) ! DIFFR
!		 NVel = nvel ! NVEL
!		 VelSep = velsep ! VELSEP 
		 hc_VelPlate = velplate ! VELPLATE
		 INJ_Ring_Number = injir ! Ring number for particle injection.
		 INJ_Area_Lower_Bound = injf1 ! Injection area lower bound.
		 INJ_Area_Upper_Bound = injf2 ! Injection area upper bound.
		 Number_Injected_Particles = injnum ! Number of injected particles.
		 Z_Charge = chzo ! CHZO
		 ZO_Temp_Grad_Parameter = czo ! CZO

		! Array data.
		 HC_MTCProb = 0 ! Note: Size maxnks*maxnrs*Number_HC_Species.

		! Note:  MATT and MATP are set in HC_Follow because they
		!        are local variables in DIVIMP.

                 ! jdemod - unstructured input quantity which adjusts the collisionality - it was passed as an argument
                 ! rather than used as a global value so HC code needs adjustment
                 hc_sf_tau = sf_tau
                 

                 
	End Subroutine Initialize_Global_Prop_Data

	Subroutine Initialize_Global_Geom_Data ()
	
		! Every good Fortran 90 program has...
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_diagvel
		Implicit None
		
		! Include common blocks.
c	Include 'params' ! Contains maxnrs,maxnks,maxpts,maxnds,isect
c	Include 'cgeom' ! Contains nds,nrs,nrs2,rmax,zmax,refct
c	Include 'comtor' ! Contains irspec,pcnt(COMTOR2),cgridopt(COMTOR2)
c	Include 'diagvel' ! Contains maxvnks

		! Global geometry properties.
		 Cell_TI = ikti ! IKTI, Last cell counting from the inner target.
		 Cell_TO = ikto ! IKTO, First cell counting from the inner target.
		 Max_Rings = maxnrs ! Maximum rings in the grid.
		 Max_Cells_Per_Ring = maxnks ! Maximum cells per ring in the grid.
		 Max_Points = maxpts ! Maximum points in the wall definition.
		 Max_Target_Cells = maxnds ! Maximum cells along the target.
		 Num_Sectors = isect ! Number of sectors.
		 Num_Upper_Rings = nrs ! Number of rings in grid.
		 Num_Lower_Rings = nrs2 ! Number of rings in grid.
		 Special_Ring = irspec ! Ring index to begin hi-res study.
		 Core_Ring = ircore ! Core ring.
		 Inner_SOL_Ring = irsep ! Main plasma ring.
		 Outer_SOL_Ring = irsep2 ! Main plasma ring.
		 Inner_Wall_Ring = irwall ! Wall ring.
		 Outer_Wall_Ring = irwall2 ! Wall ring.
		 Upper_Trap_Ring = irtrap ! Upper trap rings.
		 Lower_Trap_Ring = irtrap2 ! Lower trap rings.
		 First_Wall_Index = wlwall1 ! First outer wall index.
		 Last_Wall_Index = wlwall2 ! Last outer wall index.
		 Num_Wall_Points = wallpts ! Number of wall index points.
		 First_Trap_Index = wltrap1 ! First trap index.
		 Last_Trap_Index = wltrap2 ! Last trap index.
		 Num_Target_Cells = nds ! Number of target cells.
		 Inner_Target_Points = ndsin ! Number of inner target points.
		 End_Target_Point = ndsin3 ! 
		 R_Plasma_Centre = R0 ! X coordinate of plasma centre.
		 Z_Plasma_Centre = Z0 ! Y coordinate of plasma centre.
		 R_X_Point = rxp ! X coordinate of X point.
		 Z_X_Point = zxp ! Y coordinate of X point.
		 Normal_Measure_To_X_Axis = csnorm ! Measure theta from this many degrees from x=0.
		 Num_Boundary_Points = pcnt ! Number of total boundary points.
		 Maximum_R = rmax ! Maximum R on grid.
		 Maximum_Z = zmax ! Maximum Z on grid.
		 Grid_Error = .false. ! GRIDERR, False means no grid error (written by gridpos).

                 !
                 ! jdemod - refct isn't reliable to indicate grid geometry and is being replaced by the xpoint_up
                 !          value calculated based on grid geometry
                 !
                 hc_xpoint_up = xpoint_up ! set the overall grid geometry
                 ! jdemod
                 ! note: this doesn't work since the variable names are the same ... the one in scope is 
                 !       the one found in comtor
		 hc_ReFCT = refct ! Reflection, set in TAU from GRID2D geometry.
		 Max_Velocity_Cells_Per_Ring = maxnks ! switch to maxnks - MAXVNKS, Used for ion velocity counting.
		
	End Subroutine Initialize_Global_Geom_Data

	Subroutine Initialize_DIVIMP_Options_Table

		! Every good Fortran 90 program has...
      use mod_cgeom
      use mod_comtor
      use mod_diagvel
      use mod_reiser_com
      use mod_promptdep
      use mod_fperiph_com
      use mod_driftvel
		Implicit None
		

		! DIVIMP input options related to global properties.
		 Control_Switch = cneuta ! Begin particles as neutrals or ions.
		 Impurity_Neutral_Vel_Opt = cneutvel ! Impurity neutral velocity type option.
		 Collision_Opt = cioptb ! Collision option.
		 Friction_Opt = cioptc ! Friction option.
		 Heating_Opt = cioptd ! Heating option.
		 Injection_Opt = ciopte ! Injection option.
		 SOL_Opt = cioptf ! SOL option.
		 TEB_Coeff_Opt = cioptm ! TEB gradient coefficient option.
		 TIB_Coeff_Opt = cioptn ! TIB gradient coefficient option.
		 Reiser_Opt = cioptr ! Reiser collision option.
		 Neutral_Reflection_Opt = nrfopt ! Neutral reflection option.
		 Neutral_Mom_Coll_Opt = mtcopt ! Neutral momentum collision option from DIVIMP input file.
		 Far_Periphery_Opt = fpopt ! Far periphery option.
		 Far_Periphery_Recycle_Opt = fpropt ! Far periphery recycle as neutrals relaunched.
		 Normal_Measure_Opt = cneute ! Normal option from DIVIMP (0-2).
		 Target_Position_Opt = ctargopt ! Position of target (0-6).
		 Non_Orthogonal_Grid_Opt = northopt ! NONORTH, Non-orthogonal treatment of grid option (0-3).
		 Prompt_Deposition_Opt = 
     >            prompt_depopt ! Prompt deposition option (0-1).
		 Stop_Ion_In_Core_Opt = cstop ! Stop following ions when they reach the core (0-1).
		 First_Diffuse_Opt = cdifop ! Change calculation of SPARA depending on time (0-3).
		 Poloidal_Drift_Opt = cpdrft ! SOL poloidal drift option (0-1).
		 RZ_Opt = rzopt ! Calculate actual R,Z coordinates (0-1).
		 Neutral_Init_Pos_Opt = init_pos_opt ! Determines neutral position after neutralization (0-1).
		 Self_Sputter_Opt = cselfs ! Self-sputter option.
		 Sputter_Opt = cneutd ! Chemical sputter option.
		 Target_Mirror_Opt = cmiropt ! Target mirror option.
		 Ion_Wall_Opt =  cionr ! Ion wall option.
		! DIVIMP input options related to global geometry.
		 Grid_Opt = cgridopt ! Type of standard grid to use from the DIVIMP input file.
	End Subroutine Initialize_DIVIMP_Options_Table
		
	Subroutine Initialize_DIVIMP_HC_Data_Table ()

		! Required modules.
		!Use ComHC ! Includes all launch region numbers.

		! Every good Fortran 90 program has...
		Implicit None

		! Included common blocks.

		! Set all regions to not-used.  Note:  These are updated
		! in HC_Starting_Position as new molecules are launched.
		 HC_Region_Used = 0 ! Set all regions to non-contributing (0) state.

		! Common hydrocarbon particle FATE definition.
		 Particle_Fate (1) = 
     >            'HC Time = TMax'				! IFATE = 1
		! Neutral particle FATE definition.
		 Particle_Fate (10) = 
     >            'HC Neutral Reached Vessel Wall'		! IFATE = 10, NEUT IFATE=1
		 Particle_Fate (11) = 
     >            'HC Neutral Reached Centre Plasma' 		! IFATE = 11, NEUT IFATE=2
		 Particle_Fate (12) = 
     >            'HC Neutral Struck Target' 			! IFATE = 12, NEUT IFATE=4
		 Particle_Fate (13) = 
     >            'HC Neutral Reduced to C+' 			! IFATE = 13, NEUT IFATE=5
		 Particle_Fate (14) = 
     >            'HC Neutral Failed Launch' 			! IFATE = 14, NEUT IFATE=6
		 Particle_Fate (15) = 
     >            'HC Neutral Hit WBC Boundary'			! IFATE = 15, NEUT IFATE=7 (added by AM)
		 Particle_Fate (19) = 
     >            'HC Neutral Time = TMax'			! IFATE = 19, NEUT IFATE=3
		! Ion particle FATE definition.
		 Particle_Fate (20) = 
     >            'HC Ion Reached Vessel Wall'			! IFATE = 20, DIV IFATE=1
		 Particle_Fate (21) = 
     >            'HC Ion Reached Main Plasma'			! IFATE = 21, DIV IFATE=7
		 Particle_Fate (22) = 
     >            'HC Ion Struck Target'			! IFATE = 22, DIV IFATE=2
		 Particle_Fate (23) = 
     >            'HC Ion Reduced to C+'			! IFATE = 23, DIV IFATE=5
		 Particle_Fate (24) = 
     >            'HC Ion Removed'				! IFATE = 24, DIV IFATE=8
		 Particle_Fate (25) = 
     >            'HC Ion Hit Far-Periphery Target'		! IFATE = 25, DIV IFATE=9
		 Particle_Fate (26) = 
     >            'HC Ion Prompt Deposition' 			! IFATE = 26, DIV IFATE=10
		 Particle_Fate (27) = 
     >            'HC Ion Hit WBC Boundary'  			! IFATE = 27, DIV IFATE=11 (added by AM)
		 Particle_Fate (29) = 
     >            'HC Ion Time = TMax'				! IFATE = 29, DIV IFATE=3

		! Name all launch regions for output purposes.
                If (HC_Launch_Reg_Target_1_Dist .eq.
     >            HC_Launch_Reg_Target_1_Pin) Then
		    HC_Region_Names
     >               (HC_Launch_Reg_Target_1_Dist) = "Target 1"
		    HC_Region_Names
     >               (HC_Launch_Reg_Target_1_Pin) = "Target 1"
                Else
		    HC_Region_Names
     >               (HC_Launch_Reg_Target_1_Dist) = "Tar1 Dist"
		    HC_Region_Names
     >               (HC_Launch_Reg_Target_1_Pin) = "Tar1 Pin"
                End If
		   
                If (HC_Launch_Reg_Target_2_Dist .eq.
     >            HC_Launch_Reg_Target_2_Pin) Then
		    HC_Region_Names
     >               (HC_Launch_Reg_Target_2_Dist) = "Target 2"
		    HC_Region_Names
     >               (HC_Launch_Reg_Target_2_Pin) = "Target 2"
                Else
		    HC_Region_Names
     >               (HC_Launch_Reg_Target_2_Dist) = "Tar2 Dist"
		    HC_Region_Names
     >               (HC_Launch_Reg_Target_2_Pin) = "Tar2 Pin"
                End If

                If (HC_Launch_Reg_Wall_Homo .eq.
     >            HC_Launch_Reg_Wall_Dist) Then
		    HC_Region_Names
     >               (HC_Launch_Reg_Wall_Homo) = "Wall"
		    HC_Region_Names
     >               (HC_Launch_Reg_Wall_Dist) = "Wall"
                Else
		    HC_Region_Names
     >               (HC_Launch_Reg_Wall_Homo) = "Wall Homo"
		    HC_Region_Names
     >               (HC_Launch_Reg_Wall_Dist) = "Wall Dist"
                End If

                If (HC_Launch_Reg_Free_Space_PT .eq.
     >            HC_Launch_Reg_Free_Space_2D) Then
                    HC_Region_Names
     >               (HC_Launch_Reg_Free_Space_PT) = "Free Sp"
		    HC_Region_Names
     >               (HC_Launch_Reg_Free_Space_2D) = "Free Sp"
		Else
                    HC_Region_Names
     >               (HC_Launch_Reg_Free_Space_PT) = "FS Point"
		    HC_Region_Names
     >               (HC_Launch_Reg_Free_Space_2D) = "FS 2D"
		End If
		
                If (HC_Launch_Reg_Sput_Target_1 .eq.
     >              HC_Launch_Reg_Sput_Target_2 .and.
     >		    HC_Launch_Reg_Sput_Target_2 .eq.
     >		    HC_Launch_Reg_Sput_Wall .and.
     >              HC_Launch_Reg_Sput_Wall .eq.
     >              HC_Launch_Reg_Sput_FS_PT .and.
     >              HC_Launch_Reg_Sput_FS_PT .eq.
     >              HC_Launch_Reg_Sput_FS_2D) Then
		    HC_Region_Names
     >               (HC_Launch_Reg_Sput_Target_1) = "Sputtered"
		    HC_Region_Names
     >               (HC_Launch_Reg_Sput_Target_2) = "Sputtered"
		    HC_Region_Names
     >               (HC_Launch_Reg_Sput_Wall) = "Sputtered"
		    HC_Region_Names
     >               (HC_Launch_Reg_Sput_FS_PT) = "Sputtered"
		    HC_Region_Names
     >               (HC_Launch_Reg_Sput_FS_2D) = "Sputtered"
                ElseIf (HC_Launch_Reg_Sput_Target_1 .eq.
     >              HC_Launch_Reg_Sput_Target_2 .and.
     >		    HC_Launch_Reg_Sput_Target_2 .eq.
     >		    HC_Launch_Reg_Sput_Wall .and.
     >              HC_Launch_Reg_Sput_FS_PT .eq.
     >              HC_Launch_Reg_Sput_FS_2D) Then
		    HC_Region_Names
     >               (HC_Launch_Reg_Sput_Target_1) = "Sput W+T"
		    HC_Region_Names
     >               (HC_Launch_Reg_Sput_Target_2) = "Sput W+T"
		    HC_Region_Names
     >               (HC_Launch_Reg_Sput_Wall) = "Sput W+T"
		    HC_Region_Names
     >               (HC_Launch_Reg_Sput_FS_PT) = "Sput FS"
		    HC_Region_Names
     >               (HC_Launch_Reg_Sput_FS_2D) = "Sput FS"
		ElseIf (HC_Launch_Reg_Sput_Target_1 .eq.
     >              HC_Launch_Reg_Sput_Target_2 .and.
     >		    HC_Launch_Reg_Sput_Target_2 .ne.
     >		    HC_Launch_Reg_Sput_Wall .and.
     >              HC_Launch_Reg_Sput_FS_PT .eq.
     >              HC_Launch_Reg_Sput_FS_2D) Then
		    HC_Region_Names
     >               (HC_Launch_Reg_Sput_Target_1) = "Sput Targ"
		    HC_Region_Names
     >               (HC_Launch_Reg_Sput_Target_2) = "Sput Targ"
		    HC_Region_Names
     >               (HC_Launch_Reg_Sput_Wall) = "Sput Wall"
		    HC_Region_Names
     >               (HC_Launch_Reg_Sput_FS_PT) = "Sput FS"
		    HC_Region_Names
     >               (HC_Launch_Reg_Sput_FS_2D) = "Sput FS"
                ElseIf (HC_Launch_Reg_Sput_Target_1 .ne.
     >              HC_Launch_Reg_Sput_Target_2 .and.
     >		    HC_Launch_Reg_Sput_Target_2 .ne.
     >		    HC_Launch_Reg_Sput_Wall .and.
     >              HC_Launch_Reg_Sput_FS_PT .eq.
     >              HC_Launch_Reg_Sput_FS_2D) Then
		    HC_Region_Names
     >               (HC_Launch_Reg_Sput_Target_1) = "Sput Tar 1"
		    HC_Region_Names
     >               (HC_Launch_Reg_Sput_Target_2) = "Sput Tar 2"
		    HC_Region_Names
     >               (HC_Launch_Reg_Sput_Wall) = "Sput Wall"
		    HC_Region_Names
     >               (HC_Launch_Reg_Sput_FS_PT) = "Sput FS"
		    HC_Region_Names
     >               (HC_Launch_Reg_Sput_FS_2D) = "Sput FS"
     		Else
		    HC_Region_Names
     >               (HC_Launch_Reg_Sput_Target_1) = "Sput Tar 1"
		    HC_Region_Names
     >               (HC_Launch_Reg_Sput_Target_2) = "Sput Tar 2"
		    HC_Region_Names
     >               (HC_Launch_Reg_Sput_Wall) = "Sput Wall"
		    HC_Region_Names
     >               (HC_Launch_Reg_Sput_FS_PT) = "Sput FSPT"
		    HC_Region_Names
     >               (HC_Launch_Reg_Sput_FS_2D) = "Sput FS2D"
		End If		

                If (HC_Launch_Reg_Refl_Target_1 .eq.
     >              HC_Launch_Reg_Refl_Target_2 .and.
     >		    HC_Launch_Reg_Refl_Target_2 .eq.
     >		    HC_Launch_Reg_Refl_Wall .and.
     >              HC_Launch_Reg_Refl_Wall .eq.
     >              HC_Launch_Reg_Refl_FS_PT .and.
     >              HC_Launch_Reg_Refl_FS_PT .eq.
     >              HC_Launch_Reg_Refl_FS_2D) Then
		    HC_Region_Names
     >               (HC_Launch_Reg_Refl_Target_1) = "Reflected"
		    HC_Region_Names
     >               (HC_Launch_Reg_Refl_Target_2) = "Reflected"
		    HC_Region_Names
     >               (HC_Launch_Reg_Refl_Wall) = "Reflected"
		    HC_Region_Names
     >               (HC_Launch_Reg_Refl_FS_PT) = "Reflected"
		    HC_Region_Names
     >               (HC_Launch_Reg_Refl_FS_2D) = "Reflected"
                ElseIf (HC_Launch_Reg_Refl_Target_1 .eq.
     >              HC_Launch_Reg_Refl_Target_2 .and.
     >		    HC_Launch_Reg_Refl_Target_2 .eq.
     >		    HC_Launch_Reg_Refl_Wall .and.
     >              HC_Launch_Reg_Refl_FS_PT .eq.
     >              HC_Launch_Reg_Refl_FS_2D) Then
		    HC_Region_Names
     >               (HC_Launch_Reg_Refl_Target_1) = "Refl W+T"
		    HC_Region_Names
     >               (HC_Launch_Reg_Refl_Target_2) = "Refl W+T"
		    HC_Region_Names
     >               (HC_Launch_Reg_Refl_Wall) = "Refl W+T"
		    HC_Region_Names
     >               (HC_Launch_Reg_Refl_FS_PT) = "Refl FS"
		    HC_Region_Names
     >               (HC_Launch_Reg_Refl_FS_2D) = "Refl FS"
		ElseIf (HC_Launch_Reg_Refl_Target_1 .eq.
     >              HC_Launch_Reg_Refl_Target_2 .and.
     >		    HC_Launch_Reg_Refl_Target_2 .ne.
     >		    HC_Launch_Reg_Refl_Wall .and.
     >              HC_Launch_Reg_Refl_FS_PT .eq.
     >              HC_Launch_Reg_Refl_FS_2D) Then
		    HC_Region_Names
     >               (HC_Launch_Reg_Refl_Target_1) = "Refl Targ"
		    HC_Region_Names
     >               (HC_Launch_Reg_Refl_Target_2) = "Refl Targ"
		    HC_Region_Names
     >               (HC_Launch_Reg_Refl_Wall) = "Refl Wall"
		    HC_Region_Names
     >               (HC_Launch_Reg_Refl_FS_PT) = "Refl FS"
		    HC_Region_Names
     >               (HC_Launch_Reg_Refl_FS_2D) = "Refl FS"
                ElseIf (HC_Launch_Reg_Refl_Target_1 .ne.
     >              HC_Launch_Reg_Refl_Target_2 .and.
     >		    HC_Launch_Reg_Refl_Target_2 .ne.
     >		    HC_Launch_Reg_Refl_Wall .and.
     >              HC_Launch_Reg_Refl_FS_PT .eq.
     >              HC_Launch_Reg_Refl_FS_2D) Then
		    HC_Region_Names
     >               (HC_Launch_Reg_Refl_Target_1) = "Refl Tar 1"
		    HC_Region_Names
     >               (HC_Launch_Reg_Refl_Target_2) = "Refl Tar 2"
		    HC_Region_Names
     >               (HC_Launch_Reg_Refl_Wall) = "Refl Wall"
		    HC_Region_Names
     >               (HC_Launch_Reg_Refl_FS_PT) = "Refl FS"
		    HC_Region_Names
     >               (HC_Launch_Reg_Refl_FS_2D) = "Refl FS"
     		Else
		    HC_Region_Names
     >               (HC_Launch_Reg_Refl_Target_1) = "Refl Tar 1"
		    HC_Region_Names
     >               (HC_Launch_Reg_Refl_Target_2) = "Refl Tar 2"
		    HC_Region_Names
     >               (HC_Launch_Reg_Refl_Wall) = "Refl Wall"
		    HC_Region_Names
     >               (HC_Launch_Reg_Refl_FS_PT) = "Refl FSPT"
		    HC_Region_Names
     >               (HC_Launch_Reg_Refl_FS_2D) = "Refl FS2D"
		End If		

	End Subroutine Initialize_DIVIMP_HC_Data_Table

	Subroutine Initialize_Misc_Data_Table ()
	
		! Every good Fortran program has...
      use mod_params
      use mod_cneut
      use mod_cgeom
      use mod_comtor
      use mod_diagvel
		Implicit None
		
		! Include common blocks.
c	Include 'params' ! Contains hi, lo, pi, raddeg
c	Include 'cneut' ! Contains xprods, yprods
c	Include 'cgeom' ! Contains rmax
c	Include 'comtor' ! Contains debugn, cstepn
c	Include 'diagvel' ! Contains debugv
		
		 Calc_Hi = HI ! Machine hi value (very big number).
		 Calc_Lo = LO ! Machine lo value (very small number).
		 Pi_Value = PI ! Preset value of PI.
		 Root_2 = root2 ! SQRT(2).
!		 AMU = amu ! Atomic mass unit.
!		 ECH = ech ! Electric constant.
!		 KBoltz = kboltz ! Boltzmann's constant.
		 Rad_In_A_Deg = RADDEG ! Radians in a degree.		
		 Debug_HC_Neutral = debugn ! DEBUGN
		 Debug_HC_Ion = debugl ! DEBUGL
		 Debug_HC_V = debugv ! DEBUGV
		 Debug_HC_Prompt = debugl ! 
		 Print_Debug_x_CStep_Ion = cstepl ! CSTEPL, Print diagnostics every x timesteps.
		 Print_Debug_x_CStep_Neutral = cstepn ! CSTEPN, Print diagnostics every x timesteps.
		 Impurity_Limit = 4 ! IMPLIM, Consider moving IMPLIM setting in DIV.d6a to a common block.
		 CPU_Time_Limit = 100.0 ! CPULIM, Consider moving CPULIM to a common block.

	End Subroutine Initialize_Misc_Data_Table
	
	Subroutine Initialize_HC_Vessel_Interact ()
	
		! Use blocks.
		Use ComHC ! Includes maxnds by association with 'params'.
	
		! Every good Fortran rogram has...
		Implicit None
		
		! Declare local variables.
		Integer :: Vessel_Segment
		
		! First, fill reflection table.
		Do Vessel_Segment = 1, Max_Target_Cells
			hc_reflection_coefs (Vessel_Segment) = 
     >                    hc_reflection_coef_preset
		End Do
		
		! Note that if the sticking coef. model is not preset, the coefficient is
		! calculated based on local conditions and hydrocarbon properties.
		If (hc_sticking_coef_model .eq. 0) Then
			Do Vessel_Segment = 1, 
     >                    Max_Target_Cells
				hc_sticking_coefs (Vessel_Segment) = 
     >                            hc_sticking_coef_preset
			End Do
		End If
	
	End Subroutine Initialize_HC_Vessel_Interact
	
      End Module HC_Init_DIV_Data
