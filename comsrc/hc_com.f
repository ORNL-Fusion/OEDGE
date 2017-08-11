! -*-Fortran-*-
! Hydrocarbon Transport and Evolution Common Block
! By Adam McLean
! Under Dr. Peter Stangeby, Research Assistant David Elder
! April, 2002
!
! Common block declaration for hydrocarbon following.
! Note: Written in fixed format for compatibility with main DIVIMP common block 'params'.
! Note: Initialization is done in setup.d6a for all variables and arrays declared here, lines 2000-2200.

      Module ComHC

        ! Every good Fortran program has...
        Implicit none        

        !
        ! Make sure that the values assigned to these variables are retained across uses
        !
        logical, parameter :: debug_hc = .false.


	! Required include files.
        Include 'params' ! Contains maxnds, maxnrs, maxnks, maxpts, maxnws.
	Include 'diagvel' ! Contains nvel required by init_diag.
	
	! Define local HC constants.
	Integer, Parameter :: Number_HC_Species = 10		! Handles all HCs from C3H8 to C+. ! jdemod - set to 10 for now to conserve storage (CH4)
	Integer, Parameter :: Number_H_Species = 3		! H, H+, H2. ! jdemod - I think this comment should be H,D,T to avoid confusion
	Integer, Parameter :: Num_H_States = 4*Number_H_Species ! Note: 12 - H+, H, H2+, H2, D+, D, D2+, D2, T+, T, T2+, T2.
        !
        ! jdemod - increase number_h_products to 4 since CH4 could release 4 separate H's 
	Integer, Parameter :: Number_H_Products = 4		! Handles up to 3 H products per reaction.
        ! jdemod - changing highest_carbon_content to 1 does change some of the output formatting
        Integer, Parameter :: Highest_Carbon_Content = 3        ! Greatest number of carbon atoms supported per HC molecule.
	Integer, Parameter :: Number_HC_Reactions = 150		! Number of reactions available from database.
	Integer, Parameter :: Number_HC_IFates = 29		! Number of particle fates.
	Integer, Parameter :: Output_Unit_HC_Data = 63		! Output unit for main data file.
	Integer, Parameter :: Output_Unit_Evolve = 64		! Output unit for detailed HC evolution file.
	Integer, Parameter :: Output_Unit_Location = 65		! Output unit for detailed timestep file.
	Integer, Parameter :: Output_Unit_Scratch = 66		! Output unit for detailed timestep file.

        ! jdemod - alerts open a file with a specific name so unit number is not that important - however, some code was outputting to 
        !          hc_alert before the file was open (code in pr_hydrocarbon_options) which required moving the call to hc_begin
	Integer, Parameter :: Output_Unit_HC_Alert = 67	! Output unit for warning/error file.

	Integer, Parameter :: Output_Unit_Cpos_Pos = 68		! Output unit for initial C+ positions.
	Integer, Parameter :: Max_Number_Walks = 10000          ! MAXNWS, Number of r,z pairs to store for plotting.  jdemod - reduce from 1000000 in case of memory limits(?)
!	Integer, Parameter :: Walks_Record_Freq = 100		! Frequency of r,z pairs to store for plotting. 1=record all points.  100=record every 100th point.
!       jdemod - walks are useless for analysis unless recorded at every timestep
	Integer, Parameter :: Walks_Record_Freq = 1		! Frequency of r,z pairs to store for plotting. 1=record all points.  100=record every 100th point.
!       jdemod
	Integer, Parameter :: Max_Impurities  = maximp		! MAXIMP, Maximum number of impurity particles allowed.
	Real, Parameter :: HC_WBC_Hori_Bound = 0.50		! 50 cm in either lateral direction from launch position.
	Real, Parameter :: HC_WBC_Vert_Bound = 0.05		! 5.0 cm up from launch position.
	Character (Len=6), Parameter :: HC_Data_Type = "SigmaV"	! Decide to use reaction rates or cross sections.
	
	! HC launch region beginning points for data storage and output.  Note:  This allows
	! independance in data reporting for particles launched from different areas in different
	! ways, and alternatively combining data collection for specific areas and types of launch.
	! These are used in HC_Starting_Position to define the region at the point of launch,
	! then may be updated for sputtered or reflected hydrocarbons.
	Integer, Parameter :: Number_Regions = 12		! 2 targets, plus the wall as one region, plus free-space point and 2D launches, sputtered, reflected particles.
	Integer, Parameter :: HC_Launch_Reg_Target_1_Dist = 1	! DIVIMP launch option 0. NeutType 2.
	Integer, Parameter :: HC_Launch_Reg_Target_1_Pin = 1		! DIVIMP launch option 3. NeutType 2.
	Integer, Parameter :: HC_Launch_Reg_Target_2_Dist = 2	! DIVIMP launch option 0. NeutType 2.
	Integer, Parameter :: HC_Launch_Reg_Target_2_Pin = 2		! DIVIMP launch option 3. NeutType 2.
	Integer, Parameter :: HC_Launch_Reg_Wall_Homo = 3		! DIVIMP launch option 2. NeutType 5.
	Integer, Parameter :: HC_Launch_Reg_Wall_Dist = 3		! DIVIMP launch option 4. NeutType 5.
	Integer, Parameter :: HC_Launch_Reg_Free_Space_PT = 4	! DIVIMP launch option 1. NeutType 6.
	Integer, Parameter :: HC_Launch_Reg_Free_Space_2D = 4	! DIVIMP launch option 5. NeutType 6.
	Integer, Parameter :: HC_Launch_Reg_Sput_Target_1 = 5	! Data specific for sputtered particles. NeutType 3.
	Integer, Parameter :: HC_Launch_Reg_Sput_Target_2 = 6	! Data specific for sputtered particles. NeutType 3.
	Integer, Parameter :: HC_Launch_Reg_Sput_Wall = 7		! Data specific for sputtered particles. NeutType 3.
	Integer, Parameter :: HC_Launch_Reg_Sput_FS_PT = 8		! Data specific for sputtered particles. NeutType 3.
	Integer, Parameter :: HC_Launch_Reg_Sput_FS_2D = 8		! Data specific for sputtered particles. NeutType 3.
	Integer, Parameter :: HC_Launch_Reg_Refl_Target_1 = 9	! Data specific for reflected particles. NeutType 7.
	Integer, Parameter :: HC_Launch_Reg_Refl_Target_2 = 10	! Data specific for reflected particles. NeutType 7.
	Integer, Parameter :: HC_Launch_Reg_Refl_Wall = 11		! Data specific for reflected particles. NeutType 7.
	Integer, Parameter :: HC_Launch_Reg_Refl_FS_PT = 12		! Data specific for reflected particles. NeutType 7.
	Integer, Parameter :: HC_Launch_Reg_Refl_FS_2D = 12		! Data specific for reflected particles. NeutType 7.

        ! Hydrocarbon characteristics - default mass of hydrogen isotope
        real ::        input_HC_H_mass          ! default mass of hydrogen isotopes associated with the hydrocarbons in the simulation - can be different from the background plasma         
	
	! Hydrocarbon specific runtime options.
	Integer        HC_Follow_Option		! Determines if Hydrocarbon following should take place: 0=not activated (default), 1=activated (H15).
	Integer        HC_Higher_HCs_Option	! Determines if Hydrocarbons past CH4 are used in the code: 0=not activated (default), 1=activated (H16).
	Integer        HC_WBC_Comp_Option	! Determines if the WBC geometry model should be used for particle following and counting: 0=not activated (default), 1=activated (H17).
	
	! Hydrocarbon launch options.
	Integer        HC_Sputtering_Model		! Type of sputter release: 0-preset (default), 1-Mech&Davis&Haasz (H20).
	Integer        HC_Sputtered_HC_Species		! Preset HC species (from species table) to be released from wall/divertor; default - Methane (H21).
	Integer        HC_Evolution_Model_Primary	! Primary reaction rate data source; available: 1=E&L (default), 2=Alman/Ruzic/Brooks, 3=Janev (H22).
	Integer        HC_Evolution_Model_Secondary	! Secondary reaction rate data source; available: 0-none, 1=E&L, 2=Alman/Ruzic/Brooks, 3=Janev (default - 0, none) (H23).
	Integer        HC_Launch_Location		! Same as CNEUTB.  Available: -1 to 6, default -1=CNEUTB (H24).

        Integer        HC_Launch6_wall_index            ! jdemod - wall index for HC launch option 6 - this needs to be coordinated with DIVIMP launch options
        Logical        HC_Launch6_wall_index_set        ! jdemod - wall index must be set if option 6 is selected - warning issued if not

	Integer        HC_Launch_Angle_Velocity		! Launch initial velocity/angle options, same as CNEUTC; default -1=CNEUTC (H25).
	Integer        HC_Launch_Velocity_Model		! No Maxwellian, single, or double Maxwellian (H26).
	Real           HC_Dual_MB_Pri_Vel_Flux		! Primary contribution to dual MB velocity flux (H27).
	Real           HC_Dual_MB_Sec_Mean_Temp		! Mean temperature of hotter secondary Maxwellian (H28).
	Integer        HC_Energy_Calc			! Transition kinetics energy determination.
	Integer        HC_Self_Sputter			! Self-sputtering of HCs from C impact.

	! Hydrocarbon flight options.
	Integer        HC_Neut_Ion_Velocity	! neutral->ion initial velocity option, same as CNEUTG; default -1=CNEUTG (H30).
	Integer        HC_Ion_Neut_Angle	! ion->neutral transition angle distribution option: 0=isotropic, 1=sine biased forward, 2=S dir (H31).
	Integer        HC_Ion_Neut_Velocity	! ion->neutral transition velocity distribution option: 0=equal to ion, 1=? (H32).
	Integer	       HC_Lambda_Calc		! Use lambda=15 everywhere, or improved calculation of Sivukhin: 0=off, 1=on (H33).
	Integer        HC_Disable_Transitions	! Turn off ability for initial hydrocarbon to evolve: 0=off, 1=on (H34).
	Integer        HC_Presheath_Efield	! Use improved model of electric field of Brooks: 0=off, 1=on (H35).
	Real           HC_Efield_Drop_Fraction	! Fraction fD of potential drop in Debye region: 0.0-1.0, 0.25 typical (H36).
	Integer        HC_Efield_Cells		! Number of cells from target to apply improved sheath model (H37).
	Integer        HC_Ion_Reaction_Mult	! Multiplier for ion reactant transition probabilities.
	Integer        HC_Reaction_Kinetics     ! Option to include transition energy contribution.
	Integer        HC_Improved_Kinetics     ! Option to include 3D effects in transition direction, velocity and energy.

	! Hydrocarbon reflection options.
	Integer        HC_Neutral_Reflection_Option	! HC neutral reflection: 0=off, 1=on (H40).
	Integer        HC_Ion_Reflection_Option		! HC ion reflection: 0=off, 1=on (H41)
	Integer        HC_Reflection_Coef_Model		! Method to calculate reflection coefficient: 0=preset, 1=Alman&Ruzic, 2=Janev (H42).
	Real           HC_Reflection_Coef_Preset	! Preset default reflection coefficient: 0.0-1.0 (H43)
	Integer        HC_Reflection_species_Model	! Decide what HC species a reflected HC becomes: 0=preset reflection table, 1=Alman and Ruzic reflection data  (H44).
	Integer        HC_Reflection_Energy_Model	! Decide what HC energy a reflected HC gains: 0=preset energy (H46/H47), 1=impact energy, 2=thermal, 3=Alman and Ruzic reflection data (H45).
	Real	       HC_Refl_Energy_Neutral_Preset	! HC energy a reflected HC leaves with after a neutral particle impact (eV) (H46).
	Real	       HC_Refl_Energy_Ion_Preset	! HC energy a reflected HC leaves with after an ionized particle impact (eV) (H47).
	Integer        HC_Reflection_Angle_Model	! Decide at what angle a reflected HC ejects the vessel wall at: -1=NRFOPT, 0=off, 1=specular, 2=isotropic, 3=normal, 4=Alman and Ruzic angle data (H48).

	! Hydrocarbon sputtering options.
	Integer        HC_Sputtering_Option		! HC sputtering option: 0=off, 1=on (H50).
	Integer        HC_Sticking_Coef_Model		! Method to calculate sticking coefficient: 0=preset, 1=Alman&Ruzic, 2=Janev (H51).
	Real	       HC_Sticking_Coef_Preset		! Preset default sticking coefficient: -1.0=CTRESH, 0.0-1.0 (H52).
	Integer        HC_Sputtering_Species_Model	! Decide what HC species a sputtered HC becomes: 0=preset reflection table, 1=Alman and Ruzic reflection data (H53).
	Integer        HC_Sputtering_Energy_Model	! Decide what HC energy a sputtered HC gains: 0=preset energy (H55/H56), 1=impact energy, 2=thermal, 3=Alman and Ruzic reflection data (H54).
	Real	       HC_Sput_energy_Neutral_Preset	! HC energy a sputtered HC leaves with after a neutral particle impact (eV) (H55).
	Real	       HC_Sput_energy_Ion_Preset	! HC energy a sputtered HC leaves with after an ionized particle impact (eV) (H56).
	Integer        HC_Sputtering_Angle_Model	! Decide at what angle a sputtered HC ejects the vessel wall at: -1=NRFOPT, 0=off, 1=specular, 2=isotropic, 3=normal, 4=Alman and Ruzic angle data (H57).
	
	! Input arrays.
	Real, Dimension (maxnds) :: HC_Reflection_Coefs	! Probability of reflection. Can be individually set or set all to same preset value with H28.
	Real, Dimension (maxnds) :: HC_Sticking_Coefs	! Probability of sticking (equals 1.0 - Prob. of sputtering). Can be individually set or set all to same preset value with H36.
	
	! Additional optional output recorded or saved to file.
	Integer        HC_Coord_Print_Option	! Prints r,z position data and hydrocarbon species at each timestep: 0=off, 1=on (H50).
	Integer        HC_Evolve_Print_Option	! Prints r,z position data and hydrocarbon transition when they occur: 0=off, 1=on (H51).

        ! jdemod - add equivalent to cprint option in DIVIMP - should have the same value loaded in global_hc_assign_inputs
        integer :: hc_cprint ! General level of diagnostic printing matching the cprint value in the main code


        ! Temporary variables for debugging the sheath E-field
        integer, parameter :: nsheath_bins = 100
        real, parameter :: sheath_extent = 0.00005
        integer :: sheath_bin
        real :: sheath_zbin_data(nsheath_bins,2)
        real :: average_sheath_efield = 0.0
        real :: average_sheath_bfield = 0.0
        real :: average_sheath_bangle = 0.0
        real :: average_zsheath = 0.0
        real :: average_ssheath = 0.0
        real :: average_div_efield = 0.0
        real :: ntimes_div, ntimes_sheath
        real :: nsheath_total = 0.0

        save


      End Module ComHC
