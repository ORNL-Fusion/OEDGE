! -*-Mode:f90-*-
! Hydrocarbon_initialize.f90
! Hydrocarbon Calculation Initialization File
!
! Contains generic hydrocarbon functions including:
! -HC_Species_Ident: Returns hydrocarbon numeric identifier with its species name.
! -HC_Ident_Species: Returns hydrocarbon species name with its numeric identifier.
 
! Adam McLean under Prof. Peter Stangeby
! Research Assistant: David Elder
! October, 1999
!
! This module initializes each of the required data structures.
!
 
Module HC_Init_Lib_Data
 
  use error_handling
  Use HC_Lib_Setup
  Use ComHC ! Access to global variables related to the hydrocarbon module.
 
  ! Every good Fortran 90 program has...
  Implicit None
 
  ! Set all integer parameter variables.
  Integer, Parameter :: Num_States = Number_HC_Species ! Current 58 including doubly charged HCs over CH4 (eg. C2H++).
  ! jdemod - generalizing counts and storage - no longer separates e and p reactions
  !Integer, Parameter :: Num_Electron_Reactions = 35
  !Integer, Parameter :: Num_Proton_Reactions = 10
  integer, parameter :: max_EL_reactions = 45
  
  Integer, Parameter :: Max_Coeff_Temp_Depd = 9
  Integer, Parameter :: Max_Coeff_Energy_Depd = 9
 
  ! Define State Table
  Type (State_Entry) :: HC_State_Table (Num_States)
  Type (State_Entry) :: H_State_Table (Num_H_States)
  Type (State_Transform_Table_Type) :: HC_State_Transform_Table (Number_HC_Reactions) ! Number of reactions is from ComHC.
  Type (Reaction_Table_Type) :: HC_Reaction_Table (Num_States)
  Type (Reflection_Table_Type) :: HC_Reflection_Table (Num_States)
  Type (Sputtering_Table_Type) :: HC_Sputtering_Table (Num_States)
 
  ! Define storage for polynomial coefficient values
  ! Sigma data is as (Num_Electron_Reactions + Num_Proton_Reactions) X Max_Coeff_Temp_Depd X Max_Coeff_Energy_Depd.
  ! jdemod
  !Real, Dimension (Num_Electron_Reactions + Num_Proton_Reactions,Max_Coeff_Temp_Depd,Max_Coeff_Energy_Depd) :: Sigma_Data
  Real, Dimension (max_EL_reactions,Max_Coeff_Temp_Depd,Max_Coeff_Energy_Depd) :: Sigma_Data
 
  ! SigmaV data is as Proton Reactions X Number of polynomials X Number of cofficients per polynomial.
  ! jdemod
  !Real, Dimension (Num_Electron_Reactions + Num_Proton_Reactions,Max_Coeff_Temp_Depd,Max_Coeff_Energy_Depd) :: SigmaV_Data
  Real, Dimension (max_EL_reactions,Max_Coeff_Temp_Depd,Max_Coeff_Energy_Depd) :: SigmaV_Data

  ! jdemod - Store the masses in local variables - then this module becomes independent of hc_init_div_data and all the other 
  ! baggage that goes with it - these can be set in the initialization routine. 
  ! jdemod - back_plasma_ion_mass can be locally set to the hydrogen isotope associated with the methane. Which can then be
  !          different from the dominant plasma isotope. The mass of the dominant plasma isotope is used in the "p" interaction
  !          HC interactions with a first order correction changing the effective interaction velocity and thus sigmaV by sqrt(H_mass)
  ! Mote: changed back_plasma_ion_mass to HC_H_mass to make the distinction clear
  real,private :: HC_H_mass,impurity_ion_mass

  integer :: count_e_reactions  ! Total number of e-reactions 
  integer :: count_p_reactions  ! Total number of p-reactions 
  integer :: count_hc_reactions ! Total number of HC state transitions

 
Contains
 
  Subroutine Initialize_HC_Data (crmi)
 
    Use ComHC ! Access to global variables related to the hydrocarbon module.
    !Use HC_Init_DIV_Data ! Data arrays for output processing.
 
    ! Every good Fortran 90 program has...
    Implicit None
 
    ! jdemod - pass in the masses for the simulation
    real :: crmi

    ! Declare local variables.
    Integer, Dimension (Number_H_Species) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products) :: H_Isotope_Composition
 
    ! Some counters to get around the matrixes.
    Integer :: i, j, k
 

    ! set module variables
    ! jdemod - pass in the masses for the simulation
    ! NOTE: input_HC_H_mass is now in comhc and is read in using the unstructured input routines
    HC_H_mass = input_HC_H_mass
    impurity_ion_mass = crmi


    ! Initialize arrays.
    ! Completely fill state table.
    HC_State_Table = (/		State_Entry (1,'C+',0,0,-1,0.0),&
         State_Entry (2,'C',0,0,-1,0.0),		&
         State_Entry (3,'CH+',0,0,-1,0.0),	&
         State_Entry (4,'CH',0,0,-1,0.0),		&
         State_Entry (5,'CH2+',0,0,-1,0.0),	&
         State_Entry (6,'CH2',0,0,-1,0.0),	&
         State_Entry (7,'CH3+',0,0,-1,0.0),	&
         State_Entry (8,'CH3',0,0,-1,0.0),	&
         State_Entry (9,'CH4+',0,0,-1,0.0),	&
         State_Entry (10,'CH4',0,0,-1,0.0) /)	
 
    !						State_Entry (11,'C2++',0,0,0,0.0),	&
    !						State_Entry (12,'C2+',0,0,0,0.0),	&
    !						State_Entry (13,'C2',0,0,0,0.0),	&
    !						State_Entry (14,'C2H++',0,0,0,0.0),	&
    !						State_Entry (15,'C2H+',0,0,0,0.0),	&
    !						State_Entry (16,'C2H',0,0,0,0.0),	&
    !						State_Entry (17,'C2H2++',0,0,0,0.0),	&
    !						State_Entry (18,'C2H2+',0,0,0,0.0),	&
    !						State_Entry (19,'C2H2',0,0,0,0.0),	&
    !						State_Entry (20,'C2H3++',0,0,0,0.0),	&
    !						State_Entry (21,'C2H3+',0,0,0,0.0),	&
    !						State_Entry (22,'C2H3',0,0,0,0.0),	&
    !						State_Entry (23,'C2H4++',0,0,0,0.0),	&
    !						State_Entry (24,'C2H4+',0,0,0,0.0),	&
    !						State_Entry (25,'C2H4',0,0,0,0.0),	&
    !						State_Entry (26,'C2H5++',0,0,0,0.0),	&
    !						State_Entry (27,'C2H5+',0,0,0,0.0),	&
    !						State_Entry (28,'C2H5',0,0,0,0.0),	&
    !						State_Entry (29,'C2H6++',0,0,0,0.0),	&
    !						State_Entry (30,'C2H6+',0,0,0,0.0),	&
    !						State_Entry (31,'C2H6',0,0,0,0.0),	&
    !						State_Entry (32,'C3++',0,0,0,0.0),	&sigmav_jr_CHy_DE.dat
    !						State_Entry (33,'C3+',0,0,0,0.0),	&
    !						State_Entry (34,'C3',0,0,0,0.0),	&
    !						State_Entry (35,'C3H++',0,0,0,0.0),	&
    !						State_Entry (36,'C3H+',0,0,0,0.0),	&
    !						State_Entry (37,'C3H',0,0,0,0.0),	&
    !						State_Entry (38,'C3H2++',0,0,0,0.0),	&
    !						State_Entry (39,'C3H2+',0,0,0,0.0),	&
    !						State_Entry (40,'C3H2',0,0,0,0.0),	&
    !						State_Entry (41,'C3H3++',0,0,0,0.0),	&
    !						State_Entry (42,'C3H3+',0,0,0,0.0),	&
    !						State_Entry (43,'C3H3',0,0,0,0.0),	&
    !						State_Entry (44,'C3H4++',0,0,0,0.0),	&
    !						State_Entry (45,'C3H4+',0,0,0,0.0),	&
    !						State_Entry (46,'C3H4',0,0,0,0.0),	&
    !						State_Entry (47,'C3H5++',0,0,0,0.0),	&
    !						State_Entry (48,'C3H5+',0,0,0,0.0),	&
    !						State_Entry (49,'C3H5',0,0,0,0.0),	&
    !						State_Entry (50,'C3H6++',0,0,0,0.0),	&
    !						State_Entry (51,'C3H6+',0,0,0,0.0),	&
    !						State_Entry (52,'C3H6',0,0,0,0.0),	&
    !						State_Entry (53,'C3H7++',0,0,0,0.0),	&
    !						State_Entry (54,'C3H7+',0,0,0,0.0),	&
    !						State_Entry (55,'C3H7',0,0,0,0.0),	&
    !						State_Entry (56,'C3H8++',0,0,0,0.0),	&
    !						State_Entry (57,'C3H8+',0,0,0,0.0),	&
    !						State_Entry (58,'C3H8',0,0,0,0.0)   /)	&
 
 
    H_State_Table = (/		State_Entry (1,'H+',0,1,1,1.0),		      &
         State_Entry (2,'H',0,1,0,1.0),  	      &
         State_Entry (3,'H2+',0,2,1,2.0 * 1.0),        &
         State_Entry (4,'H2',0,2,0,2.0 * 1.0),	      &
         State_Entry (5,'D+',0,1,1,2.0), 	      &
         State_Entry (6,'D',0,1,0,2.0),  	      &
         State_Entry (7,'D2+',0,2,1,2.0 * 2.0),        &
         State_Entry (8,'D2',0,2,0,2.0 * 2.0),	      &
         State_Entry (9,'T+',0,1,1,3.0), 	      &
         State_Entry (10,'T',0,1,0,3.0), 	      &
         State_Entry (11,'T2+',0,2,1,2.0 * 3.0),       & 
         State_Entry (12,'T2',0,2,0,2.0 * 3.0) /)
 
    ! Initialize H_Isotope for background-averaged plasma analysis.
    H_Isotope_Composition = 0
 
    ! Fill the state table in with each hydrocarbon mass and charge.
    Do i = 1,Num_States,1
       Call Get_HC_Atomic_Concentration (i,HC_State_Table (i) % Carbon_Content,HC_State_Table (i) % Hydrogen_Content)
       HC_State_Table (i) % State_Charge = Get_HC_Charge (i)
       ! jdemod - NOTE: this call MUST have H_isotope_composition present and set to zero - this is used to initialize the value that
       !          find_hc_mass will return later if h_isotope_composition is not present
       HC_State_Table (i) % State_Mass = Find_HC_Mass (i,H_Isotope_Composition) ! Note:  This HC mass is the plasma background-averaged mass.
       !write (6,*) "MASSING:",i,HC_State_Table (i) % State_Name,HC_State_Table (i) % State_Mass,HC_State_Table (i) % State_Charge
 
    End Do
 
    ! jdemod - do NOT initialize the reaction type - this should be assigned when the reaction is read in and then 
    !          later used to differentiate between the treatment of different reactions
    ! jdemod - also need to initialize the entire state transform table - not just specific pieces 
    ! Fill state transform table with reaction number, type (e or p with proton
    ! reactions at the end of the array and zeros in start and end state.

    ! Initialize the state transform table

    Do i = 1,size(hc_state_transform_table)
       HC_State_Transform_Table (i) % Reaction_Number = i
       HC_State_Transform_Table (i) % Start_State = 0
       HC_State_Transform_Table (i) % End_C_States = 0
       HC_State_Transform_Table (i) % End_H_States = 0
       HC_State_Transform_Table (i) % Reaction_Identifier = ''
       HC_State_Transform_Table (i) % Reaction_Type = ''
       HC_State_Transform_Table (i) % Reaction_desc = ''
       HC_State_Transform_Table (i) % Reaction_Product_Mass = 0
       HC_State_Transform_Table (i) % Charged_Reaction_Product_Mass = 0
       HC_State_Transform_Table (i) % HC_E_Type = 0
       HC_State_Transform_Table (i) % HC_E = 0.0
       HC_State_Transform_Table (i) % H_E_Type = 0
       HC_State_Transform_Table (i) % H_E = 0.0
       HC_State_Transform_Table (i) % Sigma_TPD = 0
       HC_State_Transform_Table (i) % Sigma_Polynomial_Terms = 0
       HC_State_Transform_Table (i) % SigmaV_TPD = 0
       HC_State_Transform_Table (i) % SigmaV_Polynomial_Terms = 0
    End Do
    

    !Do i = 1,Num_Electron_Reactions,1
    !   HC_State_Transform_Table (i) % Reaction_Number = i
    !   HC_State_Transform_Table (i) % Start_State = 0
    !   HC_State_Transform_Table (i) % End_C_States = 0
    !   HC_State_Transform_Table (i) % End_H_States = 0
    !   HC_State_Transform_Table (i) % Reaction_Identifier = ''
    !   HC_State_Transform_Table (i) % Reaction_Type = 'e'
    !   HC_State_Transform_Table (i) % HC_E_Type = 0
    !   HC_State_Transform_Table (i) % HC_E = 0.0
    !   HC_State_Transform_Table (i) % H_E_Type = 0
    !   HC_State_Transform_Table (i) % H_E = 0.0
    !   HC_State_Transform_Table (i) % Sigma_TPD = 0
    !   HC_State_Transform_Table (i) % Sigma_Polynomial_Terms = 0
    !   HC_State_Transform_Table (i) % SigmaV_TPD = 0
    !   HC_State_Transform_Table (i) % SigmaV_Polynomial_Terms = 0
    !End Do
    !Do i = Num_Electron_Reactions + 1,Num_Electron_Reactions + Num_Proton_Reactions,1
    !   HC_State_Transform_Table (i) % Reaction_Number = i
    !   HC_State_Transform_Table (i) % Start_State = 0
    !   HC_State_Transform_Table (i) % End_C_States = 0
    !   HC_State_Transform_Table (i) % End_H_States = 0
    !   HC_State_Transform_Table (i) % Reaction_Identifier = ''
    !   HC_State_Transform_Table (i) % Reaction_Type = 'p'
    !   HC_State_Transform_Table (i) % HC_E_Type = 0
    !   HC_State_Transform_Table (i) % HC_E = 0.0
    !   HC_State_Transform_Table (i) % H_E_Type = 0
    !   HC_State_Transform_Table (i) % H_E = 0.0			
    !   HC_State_Transform_Table (i) % Sigma_TPD = 0
    !   HC_State_Transform_Table (i) % Sigma_Polynomial_Terms = 0
    !   HC_State_Transform_Table (i) % SigmaV_TPD = 0
    !   HC_State_Transform_Table (i) % SigmaV_Polynomial_Terms = 0
    !End Do
 
    ! Fill reaction table with reaction number and zeros.
    Do i = 1,Num_States,1
       HC_Reaction_Table (i) % State_Number = i
       HC_Reaction_Table (i) % Number_Reactions = 0
       Do j=1,size(hc_reaction_table(i)%reaction)
          HC_Reaction_Table (i) % Reaction (j) = 0
       End Do
    End Do
 
    ! Fill Sigma/SigmaV data tables with zeros.
    Sigma_Data = 0.0
    SigmaV_Data = 0.0
 
    ! Fill the reflection table with reflected state information.
    ! Note, this table as used assumes that reflection has occurred (the HC is not stuck)
    ! and hence, their energy is relatively low (high energy -> embedding of HC in vessel
    ! or divertor), but high enough to initiate chemical sputtering (~1-5eV typically).
    ! Note, currently the table is set to reflect all HCs as neutrals, fully stripped
    ! particles as themselves, and all others with equal carbon atoms, in their volatile
    ! form (volatile HC's: CH4, C2H2, C2H4, C2H6, C3H6, C3H8).
 
    Do i = 1, Num_States,1 ! 1 to 58.
       HC_Sputtering_Table (i) % State_Number = i
       Do j = 1,Highest_Carbon_Content,1 ! Should go from 1 to 3.
          HC_Sputtering_Table (i) % Sputtered_State_Number (j) = 0
       End Do
    End Do
 
    Do i = 1, Num_States,1 ! 1 to 58.
       HC_Reflection_Table (i) % State_Number = i
       HC_Reflection_Table (i) % Reflected_State_Number = 0
    End Do
 
    ! Hydrocarbon species filled to next highest volatile molecule.
    HC_Sputtering_Table (1) % Sputtered_State_Number (1) = 2
    HC_Sputtering_Table (2) % Sputtered_State_Number (1) = 2
    HC_Sputtering_Table (3) % Sputtered_State_Number (1) = 10
    HC_Sputtering_Table (4) % Sputtered_State_Number (1) = 10
    HC_Sputtering_Table (5) % Sputtered_State_Number (1) = 10
    HC_Sputtering_Table (6) % Sputtered_State_Number (1) = 10
    HC_Sputtering_Table (7) % Sputtered_State_Number (1) = 10
    HC_Sputtering_Table (8) % Sputtered_State_Number (1) = 10
    HC_Sputtering_Table (9) % Sputtered_State_Number (1) = 10
    HC_Sputtering_Table (10) % Sputtered_State_Number (1) = 10
    !		HC_Sputtering_Table (11) % Sputtered_State_Number (1) = 13
    !		HC_Sputtering_Table (12) % Sputtered_State_Number (1) = 13
    !		HC_Sputtering_Table (13) % Sputtered_State_Number (1) = 13
    !		HC_Sputtering_Table (14) % Sputtered_State_Number (1) = 19
    !		HC_Sputtering_Table (15) % Sputtered_State_Number (1) = 19
    !		HC_Sputtering_Table (16) % Sputtered_State_Number (1) = 19
    !		HC_Sputtering_Table (17) % Sputtered_State_Number (1) = 19
    !		HC_Sputtering_Table (18) % Sputtered_State_Number (1) = 19
    !		HC_Sputtering_Table (19) % Sputtered_State_Number (1) = 19
    !		HC_Sputtering_Table (20) % Sputtered_State_Number (1) = 25
    !		HC_Sputtering_Table (21) % Sputtered_State_Number (1) = 25
    !		HC_Sputtering_Table (22) % Sputtered_State_Number (1) = 25
    !		HC_Sputtering_Table (23) % Sputtered_State_Number (1) = 25
    !		HC_Sputtering_Table (24) % Sputtered_State_Number (1) = 25
    !		HC_Sputtering_Table (25) % Sputtered_State_Number (1) = 25
    !		HC_Sputtering_Table (26) % Sputtered_State_Number (1) = 31
    !		HC_Sputtering_Table (27) % Sputtered_State_Number (1) = 31
    !		HC_Sputtering_Table (28) % Sputtered_State_Number (1) = 31
    !		HC_Sputtering_Table (29) % Sputtered_State_Number (1) = 31
    !		HC_Sputtering_Table (30) % Sputtered_State_Number (1) = 31
    !		HC_Sputtering_Table (31) % Sputtered_State_Number (1) = 31
    !		HC_Sputtering_Table (32) % Sputtered_State_Number (1) = 34
    !		HC_Sputtering_Table (33) % Sputtered_State_Number (1) = 34
    !		HC_Sputtering_Table (34) % Sputtered_State_Number (1) = 34
    !		HC_Sputtering_Table (35) % Sputtered_State_Number (1) = 52
    !		HC_Sputtering_Table (36) % Sputtered_State_Number (1) = 52
    !		HC_Sputtering_Table (37) % Sputtered_State_Number (1) = 52
    !		HC_Sputtering_Table (38) % Sputtered_State_Number (1) = 52
    !		HC_Sputtering_Table (39) % Sputtered_State_Number (1) = 52
    !		HC_Sputtering_Table (40) % Sputtered_State_Number (1) = 52
    !		HC_Sputtering_Table (41) % Sputtered_State_Number (1) = 52
    !		HC_Sputtering_Table (42) % Sputtered_State_Number (1) = 52
    !		HC_Sputtering_Table (43) % Sputtered_State_Number (1) = 52
    !		HC_Sputtering_Table (44) % Sputtered_State_Number (1) = 52
    !		HC_Sputtering_Table (45) % Sputtered_State_Number (1) = 52
    !		HC_Sputtering_Table (46) % Sputtered_State_Number (1) = 52
    !		HC_Sputtering_Table (47) % Sputtered_State_Number (1) = 52
    !		HC_Sputtering_Table (48) % Sputtered_State_Number (1) = 52
    !		HC_Sputtering_Table (49) % Sputtered_State_Number (1) = 52
    !		HC_Sputtering_Table (50) % Sputtered_State_Number (1) = 52
    !		HC_Sputtering_Table (51) % Sputtered_State_Number (1) = 52
    !		HC_Sputtering_Table (52) % Sputtered_State_Number (1) = 52
    !		HC_Sputtering_Table (53) % Sputtered_State_Number (1) = 58
    !		HC_Sputtering_Table (54) % Sputtered_State_Number (1) = 58
    !		HC_Sputtering_Table (55) % Sputtered_State_Number (1) = 58
    !		HC_Sputtering_Table (56) % Sputtered_State_Number (1) = 58
    !		HC_Sputtering_Table (57) % Sputtered_State_Number (1) = 58
    !		HC_Sputtering_Table (58) % Sputtered_State_Number (1) = 58
 
    ! Simple neutralization in most cases.
    HC_Reflection_Table (1) % Reflected_State_Number = 2 ! 2
    HC_Reflection_Table (2) % Reflected_State_Number = 2 ! 2
    HC_Reflection_Table (3) % Reflected_State_Number = 10 ! 4
    HC_Reflection_Table (4) % Reflected_State_Number = 10 ! 4
    HC_Reflection_Table (5) % Reflected_State_Number = 10 ! 6
    HC_Reflection_Table (6) % Reflected_State_Number = 10 ! 6
    HC_Reflection_Table (7) % Reflected_State_Number = 10 ! 8
    HC_Reflection_Table (8) % Reflected_State_Number = 10 ! 8
    HC_Reflection_Table (9) % Reflected_State_Number = 10 ! 10
    HC_Reflection_Table (10) % Reflected_State_Number = 10 ! 10
    !		HC_Reflection_Table (11) % Reflected_State_Number = 13 ! 13
    !		HC_Reflection_Table (12) % Reflected_State_Number = 13 ! 13
    !		HC_Reflection_Table (13) % Reflected_State_Number = 13 ! 13
    !		HC_Reflection_Table (14) % Reflected_State_Number = 19 ! 16
    !		HC_Reflection_Table (15) % Reflected_State_Number = 19 ! 16
    !		HC_Reflection_Table (16) % Reflected_State_Number = 19 ! 16
    !		HC_Reflection_Table (17) % Reflected_State_Number = 19 ! 19
    !		HC_Reflection_Table (18) % Reflected_State_Number = 19 ! 19
    !		HC_Reflection_Table (19) % Reflected_State_Number = 19 ! 19
    !		HC_Reflection_Table (20) % Reflected_State_Number = 25 ! 22
    !		HC_Reflection_Table (21) % Reflected_State_Number = 25 ! 22
    !		HC_Reflection_Table (22) % Reflected_State_Number = 25 ! 22
    !		HC_Reflection_Table (23) % Reflected_State_Number = 25 ! 25
    !		HC_Reflection_Table (24) % Reflected_State_Number = 25 ! 25
    !		HC_Reflection_Table (25) % Reflected_State_Number = 25 ! 25
    !		HC_Reflection_Table (26) % Reflected_State_Number = 31 ! 28
    !		HC_Reflection_Table (27) % Reflected_State_Number = 31 ! 28
    !		HC_Reflection_Table (28) % Reflected_State_Number = 31 ! 28
    !		HC_Reflection_Table (29) % Reflected_State_Number = 31 ! 31
    !		HC_Reflection_Table (30) % Reflected_State_Number = 31 ! 31
    !		HC_Reflection_Table (31) % Reflected_State_Number = 31 ! 31
    !		HC_Reflection_Table (32) % Reflected_State_Number = 34 ! 34
    !		HC_Reflection_Table (33) % Reflected_State_Number = 34 ! 34
    !		HC_Reflection_Table (34) % Reflected_State_Number = 34 ! 34
    !		HC_Reflection_Table (35) % Reflected_State_Number = 52 ! 37
    !		HC_Reflection_Table (36) % Reflected_State_Number = 52 ! 37
    !		HC_Reflection_Table (37) % Reflected_State_Number = 52 ! 37
    !		HC_Reflection_Table (38) % Reflected_State_Number = 52 ! 40
    !		HC_Reflection_Table (39) % Reflected_State_Number = 52 ! 40
    !		HC_Reflection_Table (40) % Reflected_State_Number = 52 ! 40
    !		HC_Reflection_Table (41) % Reflected_State_Number = 52 ! 43
    !		HC_Reflection_Table (42) % Reflected_State_Number = 52 ! 43
    !		HC_Reflection_Table (43) % Reflected_State_Number = 52 ! 43
    !		HC_Reflection_Table (44) % Reflected_State_Number = 52 ! 46
    !		HC_Reflection_Table (45) % Reflected_State_Number = 52 ! 46
    !		HC_Reflection_Table (46) % Reflected_State_Number = 52 ! 46
    !		HC_Reflection_Table (47) % Reflected_State_Number = 52 ! 49
    !		HC_Reflection_Table (48) % Reflected_State_Number = 52 ! 49
    !		HC_Reflection_Table (49) % Reflected_State_Number = 52 ! 49
    !		HC_Reflection_Table (50) % Reflected_State_Number = 52 ! 52
    !		HC_Reflection_Table (51) % Reflected_State_Number = 52 ! 52
    !		HC_Reflection_Table (52) % Reflected_State_Number = 52 ! 52
    !		HC_Reflection_Table (53) % Reflected_State_Number = 58 ! 55
    !		HC_Reflection_Table (54) % Reflected_State_Number = 58 ! 55
    !		HC_Reflection_Table (55) % Reflected_State_Number = 58 ! 55
    !		HC_Reflection_Table (56) % Reflected_State_Number = 58 ! 58
    !		HC_Reflection_Table (57) % Reflected_State_Number = 58 ! 58
    !		HC_Reflection_Table (58) % Reflected_State_Number = 58 ! 58
 
  End Subroutine Initialize_HC_Data
 
  Subroutine Get_HC_Isotope_Composition (State,H_Isotope_Composition,Seed,NRand)
 
    ! Purpose:  To return the initial launch mass of a hydrocarbon molecule
    !           given a particular distribution of plasma background (i.e.
    !           hydrogen isotope concentration).
 
    ! Required external modules.
    !Use HC_Init_DIV_Data ! Contains background plasma mass.
 
    Implicit None
    Integer, Intent (In) :: State
    Integer, Dimension (Number_H_Species), Intent (Out) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (Out) :: H_Isotope_Composition
    Double Precision, Intent (In) :: Seed
    Integer, Intent (InOut) :: NRand
 
    ! Declare local variables.
    Real :: Random_Value
    Real :: D_Fraction
    Real :: D_Mobility_Enh
    Integer :: C_Num
    Integer :: H_Num
    Integer :: Number
    Real :: Avg_H_Mass ! Background average plasma mass.		
 
    ! Reset mass to begin summing components.
 
    H_Isotope_Composition = 0
    H_Num = 0
    Number = 0
 
    ! jdemod - need to initialize all variables that are used
    !        - I think this value needs to be set to 1 - a value of 0 falls through to all T in the test below
    !D_Mobility_Enh = 0.0
    D_Mobility_Enh = 1.0
 
    ! Set masses of molecular components.
    Avg_H_Mass =  HC_H_mass ! 1 for H, 2 for D, 2.5 for 50:50 DT.
 
    ! Determine fraction of Deuterium in a D/T plasma.
    D_Fraction = -1.0
    If (Avg_H_Mass .ge. 2.0 .and. Avg_H_Mass .le. 3.0) Then
       ! Plasma contains some tritium.
       D_Fraction = 1.0 - MOD (Avg_H_Mass, 2.0)
    End If
 
    ! Determine number of carbon and hydrogen atoms in the given hydrocarbon.
    H_Num = HC_State_Table (State) % Hydrogen_Content
 
    ! Determine mass of each hydrogen atom in current hydrocarbon.
    Do
       ! Exit immediately for C with no H.
       If (Number .eq. H_Num) Then
          Exit
       End If
 
       Number = Number + 1
 
       ! Call new random number and add to total.
       CALL SURAND2 (SEED, 1, Random_Value)
       NRand = NRand + 1
       If (D_Fraction .eq. -1.0) Then
          ! H
          H_Isotope_Composition (1) = H_Isotope_Composition (1) + 1
       ElseIf (Random_Value .le. (D_Fraction+D_Mobility_Enh)) Then
          ! D
          H_Isotope_Composition (2) = H_Isotope_Composition (2) + 1
       Else
          ! T
          H_Isotope_Composition (3) = H_Isotope_Composition (3) + 1
       End If
    End Do
 
  End Subroutine Get_HC_Isotope_Composition
 
  Subroutine Get_HC_Atomic_Concentration (State,Carbon_Content,Hydrogen_Content)
 
    ! Purpose:  To return the number of carbon and hydrogen atoms in a
    !           given hydrocarbon molecule given it's assigned state
    !           number via the HC_State_Table.
 
    ! Required external modules.
    Use ComHC ! Highest_Carbon_Content.
 
    Implicit None
    Integer, Intent (In) :: State
    Integer, Intent (Out) :: Carbon_Content
    Integer, Intent (Out) :: Hydrogen_Content
 
    ! Reset summing components.
    Carbon_Content = 0
    Hydrogen_Content = 0
 
    ! First, add up C mass.
    If (HC_State_Table (State) % State_Name (2:2) .ne. "H" .and. HC_State_Table (State) % State_Name (2:2) .ne. "h" .and. &
         &   HC_State_Table (State) % State_Name (2:2) .ne. "D" .and. HC_State_Table (State) % State_Name (2:2) .ne. "d" .and. &
         &   HC_State_Table (State) % State_Name (2:2) .ne. "T" .and. HC_State_Table (State) % State_Name (2:2) .ne. "t" .and. &
         &   HC_State_Table (State) % State_Name (2:2) .ne. "+" .and. HC_State_Table (State) % State_Name (2:2) .ne. " ") Then
       ! Hydrocarbon has more than one C atom (but must be less than 3).
       Read (HC_State_Table (State) % State_Name (2:2), FMT='(I1)') Carbon_Content
       If (Carbon_Content .gt. Highest_Carbon_Content) Then
          Write (Output_Unit_HC_Alert,*) "Error in Get_HC_Atomic_Concentration: HC contains more C than supported:",&
                                         &Highest_Carbon_Content
          Write (Output_Unit_HC_Alert,*) "Program stopping."
          Stop
       End If
 
       ! Next, add up H mass for hydrocarbon with more than one C.
       If (HC_State_Table (State) % State_Name (3:3) .ne. " " .and. HC_State_Table (State) % State_Name (3:3) .ne. "+") Then
          If (HC_State_Table (State) % State_Name (4:4) .ne. " " .and. HC_State_Table (State) % State_Name (4:4) .ne. "+") Then
             ! Hydrocarbon has more than one H
             Read (HC_State_Table (State) % State_Name (4:4), FMT='(I1)') Hydrogen_Content
          Else
             ! Hydrocarbon has one H.
             Hydrogen_Content = Hydrogen_Content + 1
          End If
       Else
          ! Hydrocarbon has no H, only possibly charge states.
       End If
 
    Else
       ! One C atom in HC.
       Carbon_Content = 1
 
       ! Next, add up H mass
       If (HC_State_Table (State) % State_Name (3:3) .ne. "") Then
          If (HC_State_Table (State) % State_Name (3:3) .ne. "+") Then
             ! Hydrocarbon has more than one H
             Read (HC_State_Table (State) % State_Name (3:3), FMT='(I1)') Hydrogen_Content
          Else
             ! Hydrocarbon has one H.
             Hydrogen_Content = Hydrogen_Content + 1
          End If
       Else
          ! Hydrocarbon has zero or one H.
          If (HC_State_Table (State) % State_Name (2:2) .ne. "H") Then
             ! Hydrocarbon has no H.
          Else
             ! Hydrocarbon has one H.
             Hydrogen_Content = Hydrogen_Content + 1
          End If
       End If
    End If
 
  End Subroutine Get_HC_Atomic_Concentration
 
  Real Function Find_HC_Mass (State_Number,H_Isotope_Composition)

    ! Purpose:  To return the mass of a particle given a particular state
    !           and number of each hydrogen isotope in the molecule.

    ! Required external modules.
    !Use HC_Init_DIV_Data ! Contains background plasma mass.

    Implicit None
    Integer, Intent (In) :: State_Number ! State number from 1 to 58.
    Integer, Dimension (Number_H_Species), Intent (In),optional :: H_Isotope_Composition ! Number of H/D/T's in the HC.
    !Integer, Dimension (Number_H_Products), Intent (In) :: H_Isotope_Composition ! Number of H/D/T's in the HC.
    Real :: C_Mass ! Impurity mass.
    Real :: Avg_H_Mass ! Background average plasma mass.
    Real :: H1_Mass ! H
    Real :: H2_Mass ! D
    Real :: H3_Mass ! T

    if (present(h_isotope_composition)) then 
       !
       ! If an h_isotope_composition has been provided
       !


       ! Reset mass to begin summing components.
       Find_HC_Mass = 0.0

       ! Set masses of molecular components.
       C_Mass =  Impurity_Ion_Mass ! Typically 12 for C.
       Avg_H_Mass =  HC_H_mass ! 1 for H, 2 for D, 2.5 for 50:50 DT.

       H1_Mass = 1.0
       !
       ! jdemod - these calculations of H2 and H3 mass appear incorrect 
       !        - if the isotope composition is specified then these should be 1,2,3
       !        - if one is using the background average then either isotope composition
       !          is not specified or it should be zero - in either case these should still be 
       !          assigned 1,2,3 values. 
       !
       !H2_Mass = Avg_H_Mass - MOD (Avg_H_Mass, 2.0) ! 2.0
       !H3_Mass = Avg_H_Mass + (1.0 - MOD (Avg_H_Mass, 2.0)) ! 3.0

       H2_Mass = 2.0  ! 2.0
       H3_Mass = 3.0  ! 3.0

       ! Determine type of mass to calculate.
       If (SUM (H_Isotope_Composition) .eq. 0) Then
          ! Base mass on average background.
          Find_HC_Mass = C_Mass * REAL (HC_State_Table (State_Number) % Carbon_Content) + Avg_H_Mass * REAL (HC_State_Table (&
               &State_Number) % Hydrogen_Content)
       Else
          ! Perform exact mass analysis.
          Find_HC_Mass = C_Mass * REAL (HC_State_Table (State_Number) % Carbon_Content) + &
               & H1_Mass * REAL (H_Isotope_Composition (1)) +&
               & H2_Mass * REAL (H_Isotope_Composition (2)) +&
               & H3_Mass * REAL (H_Isotope_Composition (3))
       End If

       !write(0,'(a,5i5,10f8.2)') 'H_mass:',state_number,sum(h_isotope_composition),&
       !   &    h_isotope_composition(1),h_isotope_composition(2),h_isotope_composition(3),&
       !   &    avg_h_mass,HC_H_mass,h1_mass,h2_mass,h3_mass,find_hc_mass,c_mass,REAL (HC_State_Table (State_Number) % Carbon_Content),&
       !   &    mod(avg_h_mass,2.0),avg_h_mass-mod(avg_h_mass,2.0)

    else
       ! 
       ! Note: the form without the h_isotope_composition should only be used AFTER initialization is completed
       !
       ! If a specific composition is not provided then return the mass based on the background plasma ion mass
       !
       find_hc_mass = hc_state_table(state_number)%state_mass

    endif

    !if (present(h_isotope_composition)) then 
    !   write(0,'(a,2i5,6f8.2)') 'H_mass:',state_number,sum(h_isotope_composition),&
    !      &    avg_h_mass,HC_H_mass,find_hc_mass
    !else
    !   write(0,'(a,i5,5x,6f8.2)') 'H_mass:',state_number,&
    !      &    avg_h_mass,HC_H_mass,find_hc_mass
    !endif



  End Function Find_HC_Mass
 
  Integer Function Get_HC_Charge (State_Number)

    Implicit None
    Integer, Intent (IN) :: State_Number
    Integer :: Length

    !
    ! jdemod - improve code efficiency - use precalculated values for state charge unless
    !          the values have not been initialized - in which case calculate the charge
    !          based on the state name. 
    !        - Initialization of the state_charge field occurs in the initialize_hc_data routine
    !
    if (state_number.lt.1) then 
       call errmsg('ERROR: Get_hc_charge - invalid state requested',state_number)
       get_hc_charge = 0
    else
       get_hc_charge = hc_state_table(state_number)%state_charge
    endif

    if (get_hc_charge.eq.-1) then 

       ! Find de-spaced length of State_Name.
       ! jdemod - change to using len_trim 
       length = len_trim(HC_State_Table (State_Number) % State_Name )


       !write(0,'(a,i5,6a)') 'CHARGE:',length,':',hc_state_table(state_number)%state_name(1:length),&
       !              &':',HC_State_Table (State_Number) % State_Name(Length:Length),':'
       !Length=0
       !Do
       !   If (HC_State_Table (State_Number) % State_Name (Length+1:Length+1) .eq. "") Then
       !      Exit
       !   End If
       !   Length = Length + 1
       !End Do

       ! Calculate molecular charge.
       If (HC_State_Table (State_Number) % State_Name(Length:Length) .eq. "+") Then
          Get_HC_Charge = 1
          If (HC_State_Table (State_Number) % State_Name(Length - 1:Length - 1) .eq. "+") Then
             Get_HC_Charge = 2
          End If
       Else
          ! No positive charge.
          Get_HC_Charge = 0
       End If

    endif

  End Function Get_HC_Charge
 
  Integer Function HC_Species_Ident (Species)
    ! Return identification number for a hydrocarbon species.
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Declare input/output variables.
    Character (Len=*), Intent (In) :: Species
 
    ! Declare local variables.
    Integer :: i,j

    ! jdemod - corrected possible problems for strings of unequal lengths
    integer :: targlen, speclen

    speclen = len_trim(species)
 
    ! Find species name in hydrocarbon state table.
    Do j = 1, SIZE (HC_State_Table), 1
       ! jdemod - find length of string
       targlen = len_trim(HC_State_Table (j) % State_Name)
       If (HC_State_Table (j) % State_Name(1:targlen) .eq. Species(1:speclen)) Then
          ! Found the proper hydrocarbon name.  Assign its state number to the output array.
          HC_Species_Ident = HC_State_Table (j) % State_Number
          ! Finish looping around embedded Do.
          Exit
       Else
          ! Check to see if we're at the end of the list of supported hydrocarbons.
          If (j .eq. SIZE (HC_State_Table)) Then
             Write (Output_Unit_HC_Alert,*) "Error in HC_Species_Ident: Unsupported hydrocarbon specified: ",Species(1:speclen),':'
             Write (Output_Unit_HC_Alert,*) "Program stopping."
             Write (0,*) "Error in HC_Species_Ident: Unsupported hydrocarbon specified: ",Species(1:speclen),':'
             Write (0,*) "Program stopping."
             Stop
          Else
             ! Try the next one.
             Cycle
          End If
       End If
    End Do
 
  End Function HC_Species_Ident
 
  Integer Function H_Species_Ident (Species)
    ! Return identification number for a hydrogen species.
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Declare input/output variables.
    Character (Len=*), Intent (In) :: Species
 
    ! Declare local variables.
    Integer :: i,j

    ! jdemod - corrected possible problems for strings of unequal lengths
    integer :: targlen, speclen

    speclen = len_trim(species)
 
    ! Find species name in hydrogen state table.
    Do j = 1, SIZE (H_State_Table), 1
       targlen = len_trim(H_State_Table (j) % State_Name)
       If (H_State_Table (j) % State_Name(1:targlen) .eq. Species(1:speclen)) Then
          ! Found the proper hydrogen name.  Assign its state number to the output array.
          H_Species_Ident = H_State_Table (j) % State_Number
          ! Finish looping around embedded Do.
          Exit
       Else
          ! Check to see if we're at the end of the list of supported hydrocarbons.
          If (j .eq. SIZE (H_State_Table)) Then
             Write (Output_Unit_HC_Alert,*) "Error in H_Species_Ident: Unsupported hydrogen specified: ",Species(1:speclen),':'
             Write (Output_Unit_HC_Alert,*) "Program stopping."
             Write (0,*) "Error in H_Species_Ident: Unsupported hydrogen specified: ",Species(1:speclen),':'
             Write (0,*) "Program stopping."
             Stop
          Else
             ! Try the next one.
             Cycle
          End If
       End If
    End Do
 
  End Function H_Species_Ident
 
  Character (Len=10) Function HC_Ident_Species (Ident)
    ! Return species names upon input of a hydrocarbon identification numbers.
 
    ! Every good Fortran 90 program has...		
    Implicit None
 
    ! Declare input/output variables.
    Integer, Intent (In) :: Ident
 
    ! Declare local variables.
    Integer :: i
 
    ! Check to see if number is larger than supported hydrocarbons.
    If (Ident .lt. 1 .or. Ident .gt. SIZE (HC_State_Table)) Then
       Write (Output_Unit_HC_Alert,*) "Error in HC_Ident_Species: Unsupported hydrocarbon specified. IDENT =",IDENT
       Write (Output_Unit_HC_Alert,*) "Program stopping"
       Stop
    End If
    HC_Ident_Species = HC_State_Table (Ident) % State_Name
 
  End Function HC_Ident_Species
 
  Character (Len=3) Function H_Ident_Species (Ident)
    ! Return species names upon input of a hydrogen identification numbers.
 
    ! Every good Fortran 90 program has...		
    Implicit None
 
    ! Declare input/output variables.
    Integer, Intent (In) :: Ident
 
    ! Declare local variables.
    Integer :: i
 
    ! Check to see if number is larger than supported hydrocarbons.
    If (Ident .lt. 1 .or. Ident .gt. SIZE (H_State_Table)) Then
       Write (Output_Unit_HC_Alert,*) "Error in H_Ident_Species: Unsupported hydrogen specified. IDENT =",IDENT
       Write (Output_Unit_HC_Alert,*) "Program stopping"
       Stop
    End If
    H_Ident_Species = H_State_Table (Ident) % State_Name
 
  End Function H_Ident_Species
 
  Subroutine Parse_Reaction (Reaction_Number, Eqn_Inputs, Eqn_Outputs)
 
    ! Required modules.
    Use ComHC ! Gain access to HC common block.
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Declare input/output variables.
    Integer, Intent (In) :: Reaction_Number
    Character (Len=*), Intent (In) :: Eqn_Inputs
    Character (Len=*), Intent (In) :: Eqn_Outputs
 
    ! Declare local variables.
    !Integer :: Starting_Position
    !Integer :: Position_Counter
    Integer :: Number_Of_Carbons
    Integer :: Number_Of_Hydrogens

    integer :: start_state,end_state
    integer :: start_pos,end_pos,plus_pos
    integer :: char_cnt,end_len
    integer :: product_count,in

    !write(0,'(a,a,a,i6)') 'EQN I:',eqn_inputs(1:len_trim(eqn_inputs)),':',len_trim(eqn_inputs)
    !write(0,'(a,a,a,i6)') 'EQN O:',eqn_outputs(1:len_trim(eqn_outputs)),':',len_trim(eqn_outputs)

    ! jdemod - this routine needs to be re-written to handle more complex reaction listings 
    !          1) A number specifier indicates a number of identical product species
    !          2) The order of C and H products is not enforced - H may precede C
    !          3) In addition make the identification of the source state more robust
    ! This routine expects the state formats to be listed without super or subscript characters ( _ or ^ ) 
    ! As well - a + sign will indicate either a charged state or the next product 
    !         - a single + indicates the next reaction product
    !         - a double + indicates the last state is charged and the next reaction product

    ! Source state - assume source state runs to end of first string

    start_pos = scan(eqn_inputs,'Cc')
    end_pos=scan(eqn_inputs,'CcHcDdTtpe+012345678',.true.)
 
    ! Assign starting state, which is unique in all cases.
    !write (0,*) "PARSING",Eqn_Inputs,"A",Reaction_Number,Eqn_Outputs
    ! Note that Eqn_Inputs (1-2) are "e+" or "p+".

    start_state = HC_Species_Ident (Eqn_Inputs (start_pos:end_pos)) 
    HC_State_Transform_Table (Reaction_Number) % Start_State = start_state

    if (scan(eqn_inputs(1:1),'eE').ne.0) then
       HC_State_Transform_Table(reaction_number) % Reaction_Type = 'e'
    elseif (scan(eqn_inputs(1:1),'pP').ne.0) then 
       HC_State_Transform_Table(reaction_number) % Reaction_Type = 'p'
    else
       write(0,*) 'ERROR: HC Parse_reaction : Unknown HC reaction type - not e or p: type=',eqn_inputs(1:1)
       write(Output_Unit_HC_Alert,*) 'ERROR: HC Parse_reaction : Unknown HC reaction type - not e or p: type=',eqn_inputs(1:1)
       stop
    endif


    !write(0,'(a,2i4,1x,a)') 'State::',start_pos,start_state,Eqn_Inputs (start_pos:end_pos)
 

    ! Now - scan the eqn_outpus to determine the H and C states present

    ! loop through the string eating up the characters
    char_cnt = 1 

    ! due to the possibility of non-printing characters in the reaction product string - we use a different method to 
    ! find the last character than len_trim
    !end_len = len_trim(eqn_outputs)
    end_len=scan(eqn_outputs,'CcHcDdTtpe+012345678',.true.)

    ! Zero counters
    Number_Of_Carbons = 0
    Number_Of_Hydrogens = 0

    product_count = 1

    ! Loop through the eqn_outputs string
    ! Sum up the mass of the reaction products for use later
    do 

       ! first character should be either a C for a hydrocarbon or an H for a hydrogenic species - this should always be the case
       ! the code to process each species will eat up all the other characters

       if (scan(eqn_outputs(char_cnt:char_cnt),'Cc').ne.0) then 
          ! found a C product

          start_pos = char_cnt

          ! Find the end of the product name - This is either the next "+" or the end of the string

          if (char_cnt.eq.end_len) then 

             end_pos = end_len

          else

             plus_pos = char_cnt + scan(eqn_outputs(char_cnt+1:),'+') 


             ! single "+" or no "+" is at end of string
             if (plus_pos.eq.end_len.or.plus_pos.eq.char_cnt) then 

                end_pos = end_len
             
             else

                ! check to see if this is a charged product
                if (eqn_outputs(plus_pos+1:plus_pos+1).eq.'+') then 
                   ! first plus sign found is charge indicator
                   end_pos = plus_pos

                else
                   ! plus sign found is product connector
                   end_pos = plus_pos -1

                endif
             endif

          endif

          ! need to add 2 to skip the "+" sign
          char_cnt = end_pos + 2

          end_state = HC_Species_Ident (Eqn_Outputs(start_pos:end_pos))

          !write(0,'(a,5i5,1x,a)') 'HC state:',product_count,end_state,start_pos,end_pos,plus_pos,Eqn_Outputs(start_pos:end_pos)

          do in = 1,product_count
             Number_Of_Carbons = Number_Of_Carbons + 1
             HC_State_Transform_Table (Reaction_Number) % End_C_States (Number_Of_Carbons) = end_state
             !
             ! jdemod - Note - these reaction_product_mass calculations assume the default H isotope composition
             !
             hc_state_transform_table(reaction_number)%reaction_product_mass = & 
                   &  hc_state_transform_table(reaction_number)%reaction_product_mass + &
                   &  hc_state_table(end_state)%state_mass
             !
             ! Add to mass of charged reaction products
             !
             if (hc_state_table(end_state)%state_charge .gt.0) then 
                hc_state_transform_table(reaction_number)%charged_reaction_product_mass = & 
                   &  hc_state_transform_table(reaction_number)%charged_reaction_product_mass + &
                   &  hc_state_table(end_state)%state_mass
             endif


          end do

          ! reset product count
          product_count = 1

       elseif (scan(eqn_outputs(char_cnt:char_cnt),'HhDdTt').ne.0) then 
          ! Found a Hydrogen product


          start_pos = char_cnt

          ! Find the end of the product name - This is either the next "+" or the end of the string

          if (char_cnt.eq.end_len) then 

             end_pos = end_len

          else

             plus_pos = char_cnt + scan(eqn_outputs(char_cnt+1:),'+') 


             ! single "+" or no "+" is at end of string
             if (plus_pos.eq.end_len.or.plus_pos.eq.char_cnt) then 

                end_pos = end_len
             
             else

                ! check to see if this is a charged product
                if (eqn_outputs(plus_pos+1:plus_pos+1).eq.'+') then 
                   ! first plus sign found is charge indicator
                   end_pos = plus_pos

                else
                   ! plus sign found is product connector
                   end_pos = plus_pos -1

                endif
             endif

          endif

          ! need to add 2 to skip the "+" sign
          char_cnt = end_pos + 2

          end_state = H_Species_Ident (Eqn_Outputs(start_pos:end_pos))

          !write(0,'(a,5i5,1x,a,a)') 'H  state:',product_count,end_state,start_pos,end_pos,plus_pos,&
          !                                     &Eqn_Outputs(start_pos:end_pos),':'

          do in = 1,product_count
             Number_Of_Hydrogens = Number_Of_Hydrogens + 1
             HC_State_Transform_Table (Reaction_Number) % End_H_States (Number_Of_Hydrogens) = end_state
             !
             ! Note - all Hydrogenic components are designated with an H in the reaction listings - this will 
             !        give these a mass of 1.0 - even if it happens to be D or T injected. In order to compensate
             !        for this - the state_mass is multiplied by the actual dominant hydrogenic isotope mass.
             !
             hc_state_transform_table(reaction_number)%reaction_product_mass = & 
                   &  hc_state_transform_table(reaction_number)%reaction_product_mass + &
                   &  h_state_table(end_state)%state_mass * HC_H_mass
             if (h_state_table(end_state)%state_charge .gt.0) then 
                 hc_state_transform_table(reaction_number)%charged_reaction_product_mass = & 
                     &  hc_state_transform_table(reaction_number)%charged_reaction_product_mass + &
                     &  h_state_table(end_state)%state_mass * HC_H_mass
             endif
          end do

          ! reset product count
          product_count = 1

       elseif (scan(eqn_outputs(char_cnt:char_cnt),'1234').ne.0) then 
          ! Found a product count
          ! The next product starts with a product count - read the count, save it and apply to the next state read
          read(eqn_outputs(char_cnt:char_cnt),'(i1)') product_count

          char_cnt = char_cnt + 1

       elseif (eqn_outputs(char_cnt:char_cnt).eq.'e') then 
          ! discard single electrons
          ! check to see if this is the end of the string or not - if it is the end - do nothing

          ! Find the end of the product name - This is either the next "+" or the end of the string

          if (char_cnt.eq.end_len) then 

             end_pos = end_len

          else

             plus_pos = char_cnt + scan(eqn_outputs(char_cnt+1:),'+') 


             ! single "+" or no "+" is at end of string
             if (plus_pos.eq.end_len.or.plus_pos.eq.char_cnt) then 

                end_pos = end_len
             
             else

                ! check to see if this is a charged product
                if (eqn_outputs(plus_pos+1:plus_pos+1).eq.'+') then 
                   ! first plus sign found is charge indicator
                   end_pos = plus_pos

                else
                   ! plus sign found is product connector
                   end_pos = plus_pos -1

                endif
             endif

          endif

          ! need to add 2 to skip the "+" sign

          char_cnt = end_pos + 2

       elseif (eqn_outputs(char_cnt:char_cnt).eq.'+') then 
          ! This is an error condition - something is wrong parsing the string - "+" should get used up
          Write(0,'(a)') 'ERROR: Problem in HC routine PARSE REACTION:'
          write(0,'(a)') '       Unexpected "+" in reaction product listing'
          write(0,'(a,a,a)') 'EQN OUTPUTS:',eqn_outputs(1:len_trim(eqn_outputs)),':'
          char_cnt=char_cnt + 1
       elseif (eqn_outputs(char_cnt:char_cnt).eq.' ') then 
          ! This is an error condition - something is wrong parsing the string - shouldn't be spaces
          Write(0,'(a)') 'ERROR: Problem in HC routine PARSE REACTION:'
          write(0,'(a)') '       Unexpected space in reaction product listing'
          write(0,'(a,a,a)') 'EQN OUTPUTS:',eqn_outputs(1:len_trim(eqn_outputs)),':'
          char_cnt=char_cnt + 1
       endif

       ! loop exit condition
       if (char_cnt.gt.end_len) exit

    end do

!    ! Assign target states, which there may be up to 3 carbon and 3 hydrogen atoms or molecules.
!    Starting_Position = 1
!    Position_Counter = 1
! 
!    ! Check that first character is 'C' or 'c'.
!    If(Eqn_Outputs(Position_Counter:Position_Counter).ne."C".and.Eqn_Outputs(Position_Counter:Position_Counter).ne."c")Then
!       Write (Output_Unit_HC_Alert,*) "Error in HC_LdDta: First character of output states should be C or c: ",Eqn_Outputs
!       Write (Output_Unit_HC_Alert,*) "Please check hc_xsec_ehrhardt_langer.dat."
!       Stop
!    End If
!     
!    Number_Of_Carbons = 0
!    Number_Of_Hydrogens = 0
!    Do
!       Position_Counter = Position_Counter + 1 ! First iteration will make this 2.
!       If (Eqn_Outputs (Position_Counter:Position_Counter)   .eq. " "  .or. &
!            &   Eqn_Outputs (Position_Counter:Position_Counter + 1) .eq. "+C" .or. Eqn_Outputs (Position_Counter:Position_Counter +&
!                                                               & 1) .eq. "+c" .or. &
!            &   Eqn_Outputs (Position_Counter:Position_Counter + 1) .eq. "+H" .or. Eqn_Outputs (Position_Counter:Position_Counter +&
!                                                               & 1) .eq. "+h" .or. &
!            &   Eqn_Outputs (Position_Counter:Position_Counter + 1) .eq. "+D" .or. Eqn_Outputs (Position_Counter:Position_Counter +&
!                                                               & 1) .eq. "+d" .or. &
!            &   Eqn_Outputs (Position_Counter:Position_Counter + 1) .eq. "+T" .or. Eqn_Outputs (Position_Counter:Position_Counter +&
!                                                               & 1) .eq. "+t") Then
!     
!          ! Come up to either a new carbon atom/molecule, or a new hydrogen atom/molecule product.  Record the previous spaces as a carbon atom/molecule.
!          Number_Of_Carbons = Number_Of_Carbons + 1
!          HC_State_Transform_Table (Reaction_Number) % End_C_States (Number_Of_Carbons) = HC_Species_Ident (Eqn_Outputs (&
!                                                  &Starting_Position:Position_Counter - 1))
!          Starting_Position = Position_Counter + 1
! 
!          ! Check for " ", "+H" or "+h" to go to hydrogen product loading code.
!          If (Eqn_Outputs (Position_Counter:Position_Counter + 1) .ne. "+C" .and. Eqn_Outputs (Position_Counter:Position_Counter + &
!                                                              &1) .ne. "+c") Then
!             Exit
!          End If
!          ! There must be another C.  Continue around...						
!       End If
!    End Do
! 
!    If (Eqn_Outputs (Position_Counter:Position_Counter) .ne. " ") Then
!       ! Move in front of the next "+".
!       Position_Counter = Position_Counter + 1
! 
!       Do
!          ! Move on top of the next H.
!          Position_Counter = Position_Counter + 1
! 
!          If (Eqn_Outputs (Position_Counter:Position_Counter)     .eq. " "  .or. &
!               &   Eqn_Outputs (Position_Counter:Position_Counter + 1) .eq. "+H" .or. Eqn_Outputs (Position_Counter:&
!                                                &Position_Counter + 1) .eq. "+h" .or. &
!               &   Eqn_Outputs (Position_Counter:Position_Counter + 1) .eq. "+D" .or. Eqn_Outputs (Position_Counter:&
!                                                &Position_Counter + 1) .eq. "+d" .or. &
!               &   Eqn_Outputs (Position_Counter:Position_Counter + 1) .eq. "+T" .or. Eqn_Outputs (Position_Counter:&
!                                                &Position_Counter + 1) .eq. "+t") Then
! 
!             ! Come up to either a new carbon atom/molecule, or a new hydrogen atom/molecule product.  Record the previous spaces as a carbon atom/molecule.
!             HC_State_Transform_Table (Reaction_Number) % End_H_States (Number_Of_Hydrogens + 1) = H_Species_Ident (Eqn_Outputs (&
!                                                                      &Starting_Position:Position_Counter - 1))
!             Number_Of_Hydrogens = Number_Of_Hydrogens + 1
!             Starting_Position = Position_Counter + 1
! 
!             ! Check that next position is not another hydrogen to finish.
!             If (Eqn_Outputs (Position_Counter:Position_Counter) .eq. " ") Then
!                Exit
!             End If
!             ! There must be another H.  Continue around...						
!          End If
!       End Do
!    End If
! 
  End Subroutine Parse_Reaction
 
End Module HC_Init_Lib_Data
