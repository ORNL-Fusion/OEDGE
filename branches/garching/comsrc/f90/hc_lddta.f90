! -*-Mode:f90-*-
! Hydrocarbon_data_read.f90
! Hydrocarbon Data File Reading Module
!
! Adam McLean under Prof. Peter Stangeby
! Research Assistant: David Elder
! October, 1999
!
! Module to define and read sigma and sigma-v  from DIVIMP
! hydrocarbon data file.
!
 
! Define record structure to be read into from data file.
Module HC_LdDta
 
  ! Main routine to call proper data file depending on input
  ! selection of evolution model and dataset.
  Use ComHC ! Gain access to HC common block.
  Use HC_Init_Lib_Data ! Call routine to setup data structures.
 
  ! Every good Fortran program has...
  Implicit None
 
  Integer, Parameter :: Max_Lines_Poly_Definition = 10
 
  ! Define all labels that will be encountered in data file.
  Character *12 :: Reaction_label = 'Reaction:   '
  Character *12 :: Input_label = 'Inputs:     '
  Character *12 :: Output_label = 'Outputs:    '
  Character *12 :: SigmaFT_label = 'SigmaFT:    '
  Character *12 :: SigmaVFT_label = 'SigmaVFT:   '
  Character *12 :: HC_E_Type_label = 'HC_E_Type:  '
  Character *12 :: HC_E_label = 'HC_E:       '
  Character *12 :: H_E_Type_label = 'H_E_Type:   '
  Character *12 :: H_E_label = 'H_E:        '
  Character *12 :: Datatype_label = 'Datatype:   '
  Character *12 :: Data_label = 'Data:       '
  Character *12 :: Limits_label = 'Limits:     '
  Character *12 :: Error_label = 'Error:      '
 
  ! Expected data format structure for up to both background plasma temperature
  ! dependent and target particle energy dependance.
  Structure / Data_Set /
     Character *20		Reaction
     Character *10		Eqn_Inputs ! Note:  Longest HC reaction reactants will be 'e+CxHy+', or 7 characters.
     Character *40		Eqn_Outputs ! Note:  Longest HC reaction products will be 'CHy++CHy++CHy++Hz++Hz++Hz+', or 26 characters.
     Integer *1		SigmaFT
     Integer *1		SigmaVFT
     Integer *1		Data_Type_HC_Release_Energetic
     Real			HC_Energy_Energetic
     Integer *1		Data_Type_H_Release_Energetic
     Real			H_Energy_Energetic
     Character *10		Data_Type_SigmaFT
     Integer *2		Data_Poly_SigmaFT
     Integer *2		Data_Lines_SigmaFT
     Integer *2		Data_Values_Per_Line_SigmaFT
     Integer *2		Data_Poly_Temp_SigmaFT
     Character *200		Poly_Coeff_Sigma (Max_Lines_Poly_Definition)
     Real			Limit_Sigma_Emin
     Real			Limit_Sigma_ValueEMin
     Real			Limit_Sigma_ValueEMax
     Real			Limit_Sigma_Energy_Error
     Real			Error_Energy
     Character *10		Data_Type_SigmaVFT
     Integer *2		Data_Poly_SigmaVFT
     Integer *2		Data_Lines_SigmaVFT
     Integer *2		Data_Values_Per_Line_SigmaVFT
     Integer *2		Data_Poly_Temp_SigmaVFT
     Character *200		Poly_Coeff_SigmaV (Max_Lines_Poly_Definition)
     Real			Limit_SigmaV_Tmin
     Real			Limit_SigmaV_ValueTMin
     Real			Limit_SigmaV_ValueTMax
     Real			Limit_SigmaV_Temp_Error
     Real			Error_Temperature
  End Structure
 
! Define variables used in module.
 
Integer :: Data_file_unit = 10, mult
Integer :: Temp_count = 1, i = 0, y = 0, num_elec_so_far = 0, num_prot_so_far = 0
Logical :: Unit_open
 
Logical :: End_of_file
Character *200 :: Buffer
 
! Initialize data structure.
Record /Data_Set/ Data_Structure
 
Contains
 
Subroutine Load_HC_Data ()
 
  ! Every good Fortran program has...
  Implicit None
 
  ! Local variables.
  Integer :: i

  ! initialize counters - removing unusual e and p reaction count dependence - these counters are in hc_init_lib_data
  count_e_reactions = 0
  count_p_reactions = 0
  count_hc_reactions = 0
  !

 
  If (hc_evolution_model_primary .eq. 1) Then
     Call Load_HC_Data_Ehrhardt_Langer ()
     Write (Output_Unit_HC_Data,*) "Finished loading hydrocarbon transition database: Ehrhardt and Langer"
  Else If (hc_evolution_model_primary .eq. 2) Then
     Call Load_HC_Data_Alman_Ruzic_Brooks ()
     Write (Output_Unit_HC_Data,*) "Finished loading hydrocarbon transition database: Alman and Ruzic"
  Else If (hc_evolution_model_primary .eq. 3) Then
     Call Load_HC_Data_Janev_reiter ()
     Write (Output_Unit_HC_Data,*) "Finished loading hydrocarbon transition database: Janev and Reiter"
  End If
 
!
! jdemod - the hc_evolution_model_secondary code does not work and has never been tested so it is 
!          being commented out for now. 
!
!  If (hc_evolution_model_secondary .eq. 1) Then
!     ! Only load E&L data if not already loaded.
!     If (hc_evolution_model_primary .ne. 1) Then
!        Call Load_HC_Data_Ehrhardt_Langer ()
!     End If
!  Else If (hc_evolution_model_secondary .eq. 2) Then
!     ! Only load A&R&B data if not already loaded.
!     If (hc_evolution_model_primary .ne. 2) Then
!        Call Load_HC_Data_Alman_Ruzic_Brooks ()
!     End If
!  Else If (hc_evolution_model_secondary .eq. 3) Then
!     ! Only load Janev data if not already loaded.
!     If (hc_evolution_model_primary .ne. 3) Then
!        Call Load_HC_Data_Janev ()
!     End If
!  End If


!  stop "Debugging: End of HC data load"

 
End Subroutine Load_HC_Data
 
Subroutine Load_HC_Data_Ehrhardt_Langer ()
 
  ! Every good Fortran 90 program has...
  Implicit None
 
  ! Some counters to get around the matrixes.
  Integer :: i, j, k
  Integer :: Io_Result = 0
  Integer :: C_Product_Counter
  Integer :: Current_Reaction

  character*256 :: divhome
  character*256 :: el_data_filename


  ! get divimp directory so we can look up data files

  call GetEnv('DIVHOME',divhome)

  write(0,'(a,a,a)') 'Load_HC_Data_Erhardt_Langer:DIVIMP HOME:',trim(divhome),':'

  ! set datafile name

  el_data_filename=trim(divhome)//'/hc_data/hc_xsec_ehrhardt_langer.dat'

  ! initialize
  current_reaction = 0

  call find_free_unit_number(data_file_unit)

  ! Open data file for reading.
  Open (Unit = Data_file_unit, File = el_data_filename, Status = 'Old',IOSTAT=Io_Result)
  If (Io_Result .ne. 0) then
     Write (Output_Unit_HC_Alert,*) 'ERROR in HC_LdDta: File not found ... Name = ',el_data_filename
     Write (0,*) 'ERROR in HC_LdDta: File not found ... Name = ',el_data_filename
     Write (6,*) 'ERROR in HC_LdDta: File not found ... Name = ',el_data_filename
     Stop
  End If
 
  ! Read through data file to the end.
  Do While (Io_Result .eq. 0)
 
     ! Read a line from the data file.
     Read (Unit = Data_file_unit, FMT = '(A200)',IOSTAT=Io_Result) Buffer
 
     ! Read until first character is not '!' or '#' or ' '.
     Do While ((Buffer(1:1) .eq. '!') .or. (Buffer(1:1) .eq. ' ') .or. (Buffer(1:1) .eq. '#') .or. len_trim(Buffer).le.1)
        Read (Unit = Data_file_unit, FMT = '(A200)',IOSTAT=Io_Result) Buffer
     End Do
 
     ! Interpret the line as a component of defined data structure.
     ! Check for Reaction label.
 
     If (Buffer (1:12) .eq. Reaction_label) Then
        Read (Buffer (13:), FMT = '(A10)') Data_Structure.Reaction
        Cycle
        ! Check for Inputs.
     Else If (Buffer (1:12) .eq. Input_label) Then
        Read (Buffer (13:), FMT = '(A10)') Data_Structure.Eqn_Inputs
        Cycle
        ! Check for Outputs.
     Else If (Buffer (1:12) .eq. Output_label) Then
        Read (Buffer (13:), FMT = '(A40)') Data_Structure.Eqn_Outputs
        Cycle
        ! Check for SigmaFT.
     Else If (Buffer (1:12) .eq. SigmaFT_label) Then
        Read (Buffer (13:), FMT = '(I1)') Data_Structure.SigmaFT
        Cycle
        ! Check for SigmaVFT.
     Else If (Buffer (1:12) .eq. SigmaVFT_label) Then
        Read (Buffer (13:), FMT = '(I1)')Data_Structure.SigmaVFT
        Cycle
        ! Check for hydrocarbon energetic type.
     Else If (Buffer (1:12) .eq. HC_E_Type_label) Then
        Read (Buffer (13:), FMT = '(I1)') Data_Structure.Data_Type_HC_Release_Energetic
        If (Data_Structure.Data_Type_HC_Release_Energetic .gt. 2) Then
           ! Error:  Data listed was not what was intended.
           Write (Output_Unit_HC_Alert,*) 'Error HC_Lddta:  HC energy release type not supported: ',&
                 &Data_Structure.Data_Type_HC_Release_Energetic
           Write (Output_Unit_HC_Alert,*) 'Program stopping.'
           Stop
        End If
        If (Data_Structure.Data_Type_HC_Release_Energetic .gt. 0) Then
           ! Check for hydrocarbon energy added to product if applicable.
           ! Read next line...
           Read (Unit = Data_file_unit, FMT = '(A200)',IOSTAT=io_result) Buffer	
           If (Buffer (1:12) .eq. HC_E_label) Then
              Read (Buffer (13:), FMT = *) Data_Structure.HC_Energy_Energetic
           End If
        End If
        Cycle
        ! Check for hydrogen energetic type.
     Else If (Buffer (1:12) .eq. H_E_Type_label) Then
        Read (Buffer (13:), FMT = '(I1)') Data_Structure.Data_Type_H_Release_Energetic
        If (Data_Structure.Data_Type_HC_Release_Energetic .gt. 2) Then
           ! Error:  Data listed was not what was intended.
           Write (Output_Unit_HC_Alert,*) 'Error HC_Lddta:  H energy release type not supported: ',&
                                            &Data_Structure.Data_Type_H_Release_Energetic
           Write (Output_Unit_HC_Alert,*) 'Program stopping.'
           Stop
        End If
        If (Data_Structure.Data_Type_H_Release_Energetic .gt. 0) Then
           ! Check for hydrogen energy added to product if applicable.
           ! Read next line...
           Read (Unit = Data_file_unit, FMT = '(A200)',IOSTAT=io_result) Buffer	
           If (Buffer (1:12) .eq. H_E_label) Then
              Read (Buffer (13:), FMT = *) Data_Structure.H_Energy_Energetic
           End If
        End If
        ! Check for Datatype.  Do not cycle before this check.
     End If
 
     Do While ((len_trim(Buffer) .gt. 1) .and. (io_result .eq. 0))
        ! Read next line...
        Read (Unit = Data_file_unit, FMT = '(A200)',IOSTAT=io_result) Buffer	
 
        If ((Len_Trim(Buffer) .le. 1) .or. (io_result .ne. 0)) Then
           Cycle
        End If
 
        ! Check for data type
        If (Buffer (1:12) .eq. Datatype_label) Then
           ! Assume one information line (found) followed by data lines and one limit/error line.
           ! jdemod - SCAN is not correct in this context - needs to use INDEX - scan reports any character in the set not the string
           !If ((Scan (Buffer, 'Sigma') .ne. 0) .and. (Scan (Buffer, '-V') .eq. 0)) Then
           If ((Index (Buffer, 'Sigma') .ne. 0) .and. (Index (Buffer, '-V') .eq. 0)) Then
 
              ! Found a Sigma data set.  Check to see if it's already been defined.
              If (Data_Structure.Data_Type_SigmaFT .ne. '') Then
                 ! Data type is not defined yet.  Read the next line...
                 Read (Buffer (13:), FMT = '(A10)') Data_Structure.Data_Type_SigmaFT
                 If (Data_Structure.SigmaFT .eq. 1) Then
                    ! Found a target particle independent sigma data set.  Read in Data line.
                    Read (Unit = Data_file_unit, FMT = '(A12,I2,I2,I2)') Data_label,Data_Structure.Data_Poly_SigmaFT,&
                          &Data_Structure.Data_Lines_SigmaFT,Data_Structure.Data_Values_Per_Line_SigmaFT
                 Else If (Data_Structure.SigmaFT .eq. 2) Then
                    ! Found a target particle dependent sigma data set.  Read in Data lines.
                    Read (Unit = Data_file_unit, FMT = '(A12,I2,I2,I2,I2)') Data_label, Data_Structure.Data_Poly_SigmaFT,&
            &Data_Structure.Data_Lines_SigmaFT,Data_Structure.Data_Values_Per_Line_SigmaFT,Data_Structure.Data_Poly_Temp_SigmaFT
                 Else
                    ! Error:  Data listed was not what was intended.
                    Write (Output_Unit_HC_Alert,*) 'Error HC_Lddta:  This data has not been specified: ',Buffer
                    Stop
                 End If
                 ! Read data lines
                 Do i=1,Data_Structure.Data_Lines_SigmaFT,1
                    Read (Unit = Data_file_unit, FMT = '(A200)') Data_Structure.Poly_Coeff_Sigma (i)
                 End Do
              Else
                 Write (Output_Unit_HC_Alert,*) 'Error HC_Lddta:  This value has already been defined: ',Buffer
                 Exit
              End If
              ! Read next line...
 
              Read (Unit = Data_file_unit, FMT = '(A200)',IOSTAT=io_result) Buffer
              ! Check for Limits or Errors
              If (Buffer (1:12) .eq. Limits_label) Then
                 ! Check to make sure we've saved a single-dependence series.
                 If (Data_Structure.SigmaFT .eq. 1) Then
                    Read (Buffer (13:), FMT = *) Data_Structure.Limit_Sigma_Emin,Data_Structure.Limit_Sigma_ValueEMin,&
                        &Data_Structure.Limit_Sigma_ValueEMax,Data_Structure.Limit_Sigma_Energy_Error
                 Else
                    Write (Output_Unit_HC_Alert,*) 'Error HC_Lddta:  Limit information found though no single dependence data has &
                                                   &been written: ',Buffer
                    Exit
                 End If
                 ! Check for Errors
              Else If (Buffer (1:12) .eq. Error_label) Then
                 ! Check to make sure we've saved a double-dependence series.
                 If (Data_Structure.SigmaFT .eq. 2) Then
                    Read (Buffer (13:), FMT = *) Data_Structure.Error_Energy
                 Else
                    Write (Output_Unit_HC_Alert,*) 'Error HC_Lddta:  Sigma information should be in Limit form: ', Buffer
                    Exit
                 End If
              Else
                 Write (Output_Unit_HC_Alert,*) 'Data Sigma must be followed by limit information ', Buffer
                 Exit
              End If
              Cycle
           ! jdemod - SCAN is not correct in this context - needs to use INDEX - scan reports any character in the set not the string
           !Else If (Scan (Buffer, 'Sigma-V') .ne. 0) Then
           Else If (Index (Buffer, 'Sigma-V') .ne. 0) Then
              ! Found a SigmaV data set.   Check to see if it's already been defined.
              If (Data_Structure.Data_Type_SigmaVFT .ne. '') Then
                 ! Data type is not defined yet.  Read the next line...
                 Read (Buffer (13:), FMT = '(A10)') Data_Structure.Data_Type_SigmaVFT
                 If (Data_Structure.SigmaVFT .eq. 1) Then
                    ! Found a target particle independent sigmaV data set.  Read in Data line.
                    Read (Unit = Data_file_unit, FMT = '(A12,I2,I2,I2)') Data_label,Data_Structure.Data_Poly_SigmaVFT,&
                         &Data_Structure.Data_Lines_SigmaVFT,Data_Structure.Data_Values_Per_Line_SigmaVFT
                 Else If (Data_Structure.SigmaVFT .eq. 2) Then
                    ! Found a target particle dependent sigmaV data set.  Read in Data lines.
                    Read (Unit = Data_file_unit, FMT = '(A12,I2,I2,I2,I2)') Data_label,Data_Structure.Data_Poly_SigmaVFT,&
                          &Data_Structure.Data_Lines_SigmaVFT,Data_Structure.Data_Values_Per_Line_SigmaVFT,&
                          &Data_Structure.Data_Poly_Temp_SigmaVFT
                 Else
                    ! Error:  Data listed was not what was intended.
                    Write (Output_Unit_HC_Alert,*) 'Error HC_Lddta:  This data has not been specified: ',Buffer
                    Exit
                 End If
                 ! Read data lines
                 Do i=1,Data_Structure.Data_Lines_SigmaVFT,1
                    Read (Unit = Data_file_unit, FMT = '(A200)') Data_Structure.Poly_Coeff_SigmaV (i)
                 End Do
              Else
                 Write (Output_Unit_HC_Alert,*) 'SigmaV has already been defined: ',Buffer
                 Exit
              End If
              Read (Unit = Data_file_unit, FMT = '(A200)',IOSTAT=io_result) Buffer
              ! Check for Limits or Errors
              If (Buffer (1:12) .eq. Limits_label) Then
                 ! Check to make sure we've saved a single-dependence series.
                 If (Data_Structure.SigmaVFT .eq. 1) Then
                    Read (Buffer (13:), FMT = *) Data_Structure.Limit_SigmaV_Tmin,Data_Structure.Limit_SigmaV_ValueTMin,&
                                           &Data_Structure.Limit_SigmaV_ValueTMax,Data_Structure.Limit_SigmaV_Temp_Error
                 Else
                    Write (Output_Unit_HC_Alert,*) 'Error HC_Lddta:  Limit information found though no single dependence data has &
                                                    &been written: ',Buffer
                    Exit
                 End If
                 ! Check for Errors
              Else If (Buffer (1:12) .eq. Error_label) Then
                 ! Check to make sure we've saved a double-dependence series.
                 If (Data_Structure.SigmaVFT .eq. 2) Then
                    Read (Buffer (13:), FMT = *) Data_Structure.Error_Temperature
                 Else
                    Write (Output_Unit_HC_Alert,*) 'Error HC_Lddta:  Sigma-V information should be in Error form: ', Buffer
                    Exit
                 End If
              Else
                 Write (Output_Unit_HC_Alert,*) 'Error HC_Lddta:  Data SigmaV must be followed by limit information ', Buffer
                 Exit
              End If
              Cycle
           Else
              Write (Output_Unit_HC_Alert,*) 'Error in HC_LdDta: Unknown data type found: ', Buffer
              Stop
           End If
        Else
           Write (Output_Unit_HC_Alert,*) 'Error in HC_Lddta:  Data type definition should follow: ', Buffer
           Stop
        End If
        Read (Unit = Data_file_unit, FMT = '(A200)',IOSTAT=io_result) Buffer
     End Do
 
     ! Process data in Data_Structure.
 
     ! Read transition states (initial and final) into state table.
     If (Data_Structure.Eqn_Inputs (1:2) .eq. 'e+' .or. Data_Structure.Eqn_Inputs (1:2) .eq. 'E+') Then
        ! Increase electron reaction number.
        ! jdemod - store all reactions consecutively 
        num_elec_so_far = num_elec_so_far + 1
        !Current_Reaction = num_elec_so_far
        !current_reaction = current_reaction + 1
 
        ! Load in HC_State_Transform_Table % Start_State and End_States.
        !Call Parse_Reaction (current_reaction, Data_Structure.Eqn_Inputs, Data_Structure.Eqn_Outputs)
 
        !HC_State_Transform_Table (current_reaction) % Reaction_Identifier = Data_Structure.Reaction
        !HC_State_Transform_Table (current_reaction) % HC_E_Type = Data_Structure.Data_Type_HC_Release_Energetic
        !HC_State_Transform_Table (current_reaction) % HC_E = Data_Structure.HC_Energy_Energetic
        !HC_State_Transform_Table (current_reaction) % H_E_Type = Data_Structure.Data_Type_H_Release_Energetic
        !HC_State_Transform_Table (current_reaction) % H_E = Data_Structure.H_Energy_Energetic
 
     Else If (Data_Structure.Eqn_Inputs (1:2) .eq. 'p+' .or. Data_Structure.Eqn_Inputs (1:2) .eq. 'P+' .or. &
          &        Data_Structure.Eqn_Inputs (1:2) .eq. 'h+' .or. Data_Structure.Eqn_Inputs (1:2) .eq. 'H+' .or. &
          &        Data_Structure.Eqn_Inputs (1:2) .eq. 'h+' .or. Data_Structure.Eqn_Inputs (1:2) .eq. 'H+' .or. &
          &        Data_Structure.Eqn_Inputs (1:2) .eq. 't+' .or. Data_Structure.Eqn_Inputs (1:2) .eq. 'T+') Then
 
        ! Increase proton reaction number.
        num_prot_so_far = num_prot_so_far + 1
        !Current_Reaction = Num_Electron_Reactions + num_prot_so_far
        !current_reaction = current_reaction + 1
 
        ! Load in HC_State_Transform_Table % Start_State and End_States.
        !Call Parse_Reaction (current_reaction, Data_Structure.Eqn_Inputs, Data_Structure.Eqn_Outputs)
 
        !HC_State_Transform_Table (current_reaction) % Reaction_Identifier = Data_Structure.Reaction
        !HC_State_Transform_Table (current_reaction) % HC_E_Type=Data_Structure.Data_Type_HC_Release_Energetic
        !HC_State_Transform_Table (current_reaction) % HC_E = Data_Structure.HC_Energy_Energetic
        !HC_State_Transform_Table (current_reaction) % H_E_Type=Data_Structure.Data_Type_H_Release_Energetic
        !HC_State_Transform_Table (current_reaction) % H_E = Data_Structure.H_Energy_Energetic

        !HC_State_Transform_Table (Num_Electron_Reactions + num_prot_so_far) % Reaction_Identifier = Data_Structure.Reaction
        !HC_State_Transform_Table (Num_Electron_Reactions + num_prot_so_far)%HC_E_Type=Data_Structure.Data_Type_HC_Release_Energetic
        !HC_State_Transform_Table (Num_Electron_Reactions + num_prot_so_far) % HC_E = Data_Structure.HC_Energy_Energetic
        !HC_State_Transform_Table (Num_Electron_Reactions + num_prot_so_far)%H_E_Type=Data_Structure.Data_Type_H_Release_Energetic
        !HC_State_Transform_Table (Num_Electron_Reactions + num_prot_so_far) % H_E = Data_Structure.H_Energy_Energetic
 
     Else
        Write (Output_Unit_HC_Alert,*) 'Unknown reaction!: ',Buffer
        Exit
     End If

     ! increment reaction count 
     current_reaction = current_reaction + 1
 
     ! Load in HC_State_Transform_Table % Start_State and End_States.
     ! jdemod - parse reaction now also determines the reaction_type
     Call Parse_Reaction (current_reaction, Data_Structure.Eqn_Inputs, Data_Structure.Eqn_Outputs)
 
     HC_State_Transform_Table (current_reaction) % Reaction_Identifier = Data_Structure.Reaction
     HC_State_Transform_Table (current_reaction) % HC_E_Type=Data_Structure.Data_Type_HC_Release_Energetic
     HC_State_Transform_Table (current_reaction) % HC_E = Data_Structure.HC_Energy_Energetic
     HC_State_Transform_Table (current_reaction) % H_E_Type=Data_Structure.Data_Type_H_Release_Energetic
     HC_State_Transform_Table (current_reaction) % H_E = Data_Structure.H_Energy_Energetic
 
     ! Read reaction number into reaction table and increment number of reactions possible for that initial state by one.
     HC_Reaction_Table (HC_Species_Ident (Data_Structure.Eqn_Inputs (3:))) % Number_Reactions = HC_Reaction_Table (&
                       &HC_Species_Ident (Data_Structure.Eqn_Inputs (3:))) % Number_Reactions + 1
     If ((HC_Reaction_Table (HC_Species_Ident (Data_Structure.Eqn_Inputs (3:))) % Number_Reactions) > 10) Then
        Print *,'Error:  Increase number of possible reactions from ',HC_Reaction_Table (HC_Species_Ident (&
                                                            &Data_Structure.Eqn_Inputs (3:))) % Number_Reactions
        Exit
     End If
     ! Add current reaction number to reaction table.
     HC_Reaction_Table (HC_Species_Ident (Data_Structure.Eqn_Inputs (3:))) % Reaction (HC_Reaction_Table (HC_Species_Ident (&
                                         &Data_Structure.Eqn_Inputs (3:))) % Number_Reactions) = Current_Reaction

 
     ! Read in polynomial coefficients .
     ! jdemod - SCAN is not correct in this context - needs to use INDEX - scan reports any character in the set not the string
     ! If((Scan(Data_Structure.Data_Type_SigmaFT,'Sigma').ne.0).and.(Scan(Data_Structure.Data_Type_SigmaFT,'-V').eq.0))Then
     If((index(Data_Structure.Data_Type_SigmaFT,'Sigma').ne.0).and.(index(Data_Structure.Data_Type_SigmaFT,'-V').eq.0))Then
        ! Load data in transform table for use in interpolation routine.
 
        ! jdemod - these blocks of code are the same - only difference is the reaction index
        !If (Scan (Data_Structure.Eqn_Inputs, 'e') .ne. 0) Then

        HC_State_Transform_Table (current_reaction) % Sigma_TPD = Data_Structure.SigmaFT
        HC_State_Transform_Table (current_reaction) % Sigma_Polynomial_Terms = Data_Structure.Data_Poly_SigmaFT
        HC_State_Transform_Table (current_reaction) % Sigma_Tmin_Limit = Data_Structure.Limit_Sigma_Emin
        HC_State_Transform_Table (current_reaction) % Sigma_ValueEMin_Limit = Data_Structure.Limit_Sigma_ValueEMin
        HC_State_Transform_Table (current_reaction) % Sigma_ValueEMax_Limit = Data_Structure.Limit_Sigma_ValueEMax
        HC_State_Transform_Table (current_reaction) % Sigma_Energy_Error_Limit = Data_Structure.Limit_Sigma_Energy_Error

        !Else If (Scan (Data_Structure.Eqn_Inputs, 'p') .ne. 0) Then
        !   HC_State_Transform_Table (Num_Electron_Reactions + num_prot_so_far) % Sigma_TPD = Data_Structure.SigmaFT
        !   HC_State_Transform_Table(Num_Electron_Reactions+num_prot_so_far)%Sigma_Polynomial_Terms=Data_Structure.Data_Poly_SigmaFT
        !   HC_State_Transform_Table (Num_Electron_Reactions + num_prot_so_far) % Sigma_Tmin_Limit = Data_Structure.Limit_Sigma_Emin
        !   HC_State_Transform_Table (Num_Electron_Reactions + num_prot_so_far) % Sigma_ValueEMin_Limit = &
        !                                   &Data_Structure.Limit_Sigma_ValueEMin
        !   HC_State_Transform_Table (Num_Electron_Reactions + num_prot_so_far) % Sigma_ValueEMax_Limit = &
        !                                   &Data_Structure.Limit_Sigma_ValueEMax
        !   HC_State_Transform_Table (Num_Electron_Reactions + num_prot_so_far) % Sigma_Energy_Error_Limit = &
        !                                   &Data_Structure.Limit_Sigma_Energy_Error
        !End If


        ! Read in Sigma data.  Check if it's target particle dependent.
        If (Data_Structure.SigmaFT .eq. 1) Then
           ! Sigma data is target particle indepedent (little array).
           Do i = 1,Data_Structure.Data_Lines_SigmaFT,1
              ! Store line data in temporary array of length Max number of values per line.
              Read (Data_Structure.Poly_Coeff_Sigma (i), FMT = *) &
                      & Sigma_Data (current_reaction, 1 + (i-1) * Data_Structure.Data_Values_Per_Line_SigmaFT : (i-1) * &
                             &Data_Structure.Data_Values_Per_Line_SigmaFT + Data_Structure.Data_Values_Per_Line_SigmaFT, 1)


              !If (Scan (Data_Structure.Eqn_Inputs, 'e') .ne. 0) Then
              !   ! jdemod - changed to current_reaction
              !   Read (Data_Structure.Poly_Coeff_Sigma (i), FMT = *) &
              !        & Sigma_Data (current_reaction, 1 + (i-1) * Data_Structure.Data_Values_Per_Line_SigmaFT : (i-1) * &
              !               &Data_Structure.Data_Values_Per_Line_SigmaFT + Data_Structure.Data_Values_Per_Line_SigmaFT, 1)
              !   !Read (Data_Structure.Poly_Coeff_Sigma (i), FMT = *) &
              !   !     & Sigma_Data (num_elec_so_far, 1 + (i-1) * Data_Structure.Data_Values_Per_Line_SigmaFT : (i-1) * &
              !   !            &Data_Structure.Data_Values_Per_Line_SigmaFT + Data_Structure.Data_Values_Per_Line_SigmaFT, 1)
              !Else If (Scan (Data_Structure.Eqn_Inputs, 'p') .ne. 0) Then
              !   ! jdemod - changed to current_reaction
              !   Read (Data_Structure.Poly_Coeff_Sigma (i), FMT = *) &
              !        & Sigma_Data (current_reaction, 1 + (i-1) * &
              !        &Data_Structure.Data_Values_Per_Line_SigmaFT : (i-1) * Data_Structure.Data_Values_Per_Line_SigmaFT + &
              !        &Data_Structure.Data_Values_Per_Line_SigmaFT, 1)
              !   !Read (Data_Structure.Poly_Coeff_Sigma (i), FMT = *) &
              !   !     & Sigma_Data (Num_Electron_Reactions + num_prot_so_far, 1 + (i-1) * &
              !   !     &Data_Structure.Data_Values_Per_Line_SigmaFT : (i-1) * Data_Structure.Data_Values_Per_Line_SigmaFT + &
              !   !     &Data_Structure.Data_Values_Per_Line_SigmaFT, 1)
              !End If

           End Do
        Else If (Data_Structure.SigmaFT .eq. 2) Then
           ! Sigma data is target particle dependent (big array!).
           Do i = 1,Data_Structure.Data_Lines_SigmaFT,1
              Do y = 1, Data_Structure.Data_Poly_Temp_SigmaFT,1

                 ! Store line data in temporary array of length Max number of values per line.
                 ! jdemod - changed to current_reaction
                 Read (Data_Structure.Poly_Coeff_Sigma (i), FMT = *) &
                         & Sigma_Data (current_reaction, y : y + Data_structure.Data_Values_Per_Line_SigmaFT - 1, i)


                 !If (Scan (Data_Structure.Eqn_Inputs, 'e') .ne. 0) Then
                 !   ! jdemod - changed to current_reaction
                 !   Read (Data_Structure.Poly_Coeff_Sigma (i), FMT = *) &
                 !        & Sigma_Data (current_reaction, y : y + Data_structure.Data_Values_Per_Line_SigmaFT - 1, i)
                 !   !Read (Data_Structure.Poly_Coeff_Sigma (i), FMT = *) &
                 !   !     & Sigma_Data (num_elec_so_far, y : y + Data_structure.Data_Values_Per_Line_SigmaFT - 1, i)
                 !Else If (Scan (Data_Structure.Eqn_Inputs, 'p') .ne. 0) Then
                 !   ! jdemod - changed to current_reaction
                 !   Read (Data_Structure.Poly_Coeff_Sigma (i), FMT = *) &
                 !        &Sigma_Data(current_reaction,y:y+Data_Structure.Data_Values_Per_Line_SigmaFT-1,i)
                 !   !Read (Data_Structure.Poly_Coeff_Sigma (i), FMT = *) &
                 !   !     &Sigma_Data(Num_Electron_Reactions+num_prot_so_far,y:y+Data_Structure.Data_Values_Per_Line_SigmaFT-1,i)
                 !End If

              End Do
           End Do
        Else
           Write (Output_Unit_HC_Alert,*) 'Error:  Sigma data format not defined'
           Exit
        End If
     End If

     ! jdemod - SCAN is not correct in this context - needs to use INDEX - scan reports any character in the set not the string
     !If (Scan (Data_Structure.Data_Type_SigmaVFT, 'Sigma-V') .ne. 0) Then
     If (index (Data_Structure.Data_Type_SigmaVFT, 'Sigma-V') .ne. 0) Then
        ! Load data in transform table for use in interpolation routine.

        ! jdemod - use current_reaction - note the "p" code left out a number of values
        HC_State_Transform_Table (current_reaction) % SigmaV_TPD = Data_Structure.SigmaVFT
        HC_State_Transform_Table (current_reaction) % SigmaV_Polynomial_Terms = Data_Structure.Data_Poly_SigmaVFT
        HC_State_Transform_Table (current_reaction) % SigmaV_Tmin_Limit = Data_Structure.Limit_SigmaV_Tmin
        HC_State_Transform_Table (current_reaction) % SigmaV_ValueEMin_Limit = Data_Structure.Limit_SigmaV_ValueTMin
        HC_State_Transform_Table (current_reaction) % SigmaV_ValueEMax_Limit = Data_Structure.Limit_SigmaV_ValueTMax
        HC_State_Transform_Table (current_reaction) % SigmaV_Energy_Error_Limit = Data_Structure.Limit_SigmaV_Temp_Error


        !If (Scan (Data_Structure.Eqn_Inputs, 'e') .ne. 0) Then
        !   HC_State_Transform_Table (num_elec_so_far) % SigmaV_TPD = Data_Structure.SigmaVFT
        !   HC_State_Transform_Table (num_elec_so_far) % SigmaV_Polynomial_Terms = Data_Structure.Data_Poly_SigmaVFT
        !   HC_State_Transform_Table (num_elec_so_far) % SigmaV_Tmin_Limit = Data_Structure.Limit_SigmaV_Tmin
        !   HC_State_Transform_Table (num_elec_so_far) % SigmaV_ValueEMin_Limit = Data_Structure.Limit_SigmaV_ValueTMin
        !   HC_State_Transform_Table (num_elec_so_far) % SigmaV_ValueEMax_Limit = Data_Structure.Limit_SigmaV_ValueTMax
        !   HC_State_Transform_Table (num_elec_so_far) % SigmaV_Energy_Error_Limit = Data_Structure.Limit_SigmaV_Temp_Error
        !Else If (Scan (Data_Structure.Eqn_Inputs, 'p') .ne. 0) Then
        !   HC_State_Transform_Table (Num_Electron_Reactions + num_prot_so_far) % SigmaV_TPD = Data_Structure.SigmaVFT
        !   HC_State_Transform_Table (Num_Electron_Reactions + num_prot_so_far) % SigmaV_Polynomial_Terms = &
        !                                                                 &Data_Structure.Data_Poly_SigmaVFT
        !   HC_State_Transform_Table (num_elec_so_far) % SigmaV_Energy_Error_Limit = Data_Structure.Error_Temperature
        !End If


        ! Read in SigmaV data.  Check if it's target particle dependent.
        If (Data_Structure.SigmaVFT .eq. 1) Then
           ! SigmaV data is target particle indepedent (little array).
           Do i = 1,Data_Structure.Data_Lines_SigmaVFT,1
              ! Store line data in temporary array of length Max number of values per line.
              ! jdemod - use current_reaction
              Read (Data_Structure.Poly_Coeff_SigmaV (i), FMT = *) &
                   &SigmaV_Data (current_reaction, 1 + (i-1) * Data_Structure.Data_Values_Per_Line_SigmaVFT : (i-1) * &
                   &Data_Structure.Data_Values_Per_Line_SigmaVFT + Data_Structure.Data_Values_Per_Line_SigmaVFT, 1)

              !If (Scan (Data_Structure.Eqn_Inputs, 'e') .ne. 0) Then
              !   Read (Data_Structure.Poly_Coeff_SigmaV (i), FMT = *) &
              !        &SigmaV_Data (num_elec_so_far, 1 + (i-1) * Data_Structure.Data_Values_Per_Line_SigmaVFT : (i-1) * &
              !        &Data_Structure.Data_Values_Per_Line_SigmaVFT + Data_Structure.Data_Values_Per_Line_SigmaVFT, 1)
              !Else If (Scan (Data_Structure.Eqn_Inputs, 'p') .ne. 0) Then
              !   Read (Data_Structure.Poly_Coeff_SigmaV (i), FMT = *) &
              !        &SigmaV_Data (Num_Electron_Reactions + num_prot_so_far, 1 + (i-1) * &
              !        &Data_Structure.Data_Values_Per_Line_SigmaVFT : (i-1) * Data_Structure.Data_Values_Per_Line_SigmaVFT + &
              !        &Data_Structure.Data_Values_Per_Line_SigmaVFT, 1)
              !End If

           End Do
        Else If (Data_Structure.SigmaVFT .eq. 2) Then
           ! SigmaV data is target particle dependent (big array!).
           Do i = 1,Data_Structure.Data_Lines_SigmaVFT,1
              Do y = 1, Data_Structure.Data_Poly_Temp_SigmaVFT,Data_structure.Data_Values_Per_Line_SigmaVFT
                 ! jdemod - use current_reaction

                 Read (Data_Structure.Poly_Coeff_SigmaV (i), FMT = *) &
                         &SigmaV_Data (current_reaction, y : y + Data_structure.Data_Values_Per_Line_SigmaVFT-1,i)

                 ! Store line data in temporary array of length Max number of values per line.
                 !If (Scan (Data_Structure.Eqn_Inputs, 'e') .ne. 0) Then
                 !   Read (Data_Structure.Poly_Coeff_SigmaV (i), FMT = *) &
                 !        &SigmaV_Data (num_elec_so_far, y : y + Data_structure.Data_Values_Per_Line_SigmaVFT-1,i)
                 !Else If (Scan (Data_Structure.Eqn_Inputs, 'p') .ne. 0) Then
                 !   Read (Data_Structure.Poly_Coeff_SigmaV (i), FMT = *) &
                 !        &SigmaV_Data(Num_Electron_Reactions+num_prot_so_far,y:y+Data_Structure.Data_Values_Per_Line_SigmaVFT-1,i)
                 !End If

              End Do
           End Do
        Else
           Write (Output_Unit_HC_Alert,*) 'Error:  SigmaV data format not defined'
           Exit
        End If
     End If
  End Do
  Close (Unit = Data_file_unit)


 
  count_e_reactions = num_elec_so_far
  count_p_reactions = num_prot_so_far
  count_hc_reactions = current_reaction

  write(0,*) 'EL Reaction counts:',num_elec_so_far,num_prot_so_far,count_hc_reactions

 
End Subroutine Load_HC_Data_Ehrhardt_Langer
 
Subroutine Load_HC_Data_Janev_reiter ()

  use hc_read_janev_reiter
  ! Every good Fortran 90 program has...
  Implicit None

  !
  ! Filenames to be read are hard coded at the moment
  !
  !
  integer,parameter :: max_files = 10
  character*256 :: divhome

  character*256 :: files(max_files)

  integer :: filecount
  integer :: strlen,in,i,rc
  integer :: nt,ne,it

  !
  ! Total number of files to read
  !
  filecount = 10

  call GetEnv('DIVHOME',divhome)

  write(0,'(a,a,a)') 'Load_HC_Data_Janev_Reiter:DIVIMP HOME:',trim(divhome),':'

  !
  ! e-reactions
  !

  files(1) = trim(divhome)//'/hc_data/hc_sigmav_jr_CHy_DE.dat'
  files(2) = trim(divhome)//'/hc_data/hc_sigmav_jr_CHy+_DE.dat'
  files(3) = trim(divhome)//'/hc_data/hc_sigmav_jr_CHy+_DI.dat'
  files(4) = trim(divhome)//'/hc_data/hc_sigmav_jr_CHy+_DR.dat'
  files(5) = trim(divhome)//'/hc_data/hc_sigmav_jr_CI_DI.dat'
  files(6) = trim(divhome)//'/hc_data/hc_sigmav_jr_CH_I_DI.dat'
  files(7) = trim(divhome)//'/hc_data/hc_sigmav_jr_CH2_I_DI.dat'
  files(8) = trim(divhome)//'/hc_data/hc_sigmav_jr_CH3_I_DI.dat'
  files(9) = trim(divhome)//'/hc_data/hc_sigmav_jr_CH4_I_DI.dat'

  !
  ! p-reactions
  !
  files(10) = trim(divhome)//'/hc_data/hc_sigmav_jr_CHy_CX.dat'

  call init_load_jr

  do in = 1,filecount

     strlen = len_trim(files(in))
     call load_jr(files(in)(1:strlen),rc)

  end do

 
  call set_jr_reaction_counts(count_hc_reactions,count_e_reactions,count_p_reactions)


  do in = 1,count_hc_reactions

     write(6,'(a,i4,1x,a,a4,3i6,2(1x,g12.5),a)') 'Processes:',in,hc_state_transform_table(in)%reaction_identifier,&
                       & hc_state_transform_table(in)%reaction_type,&
                       & hc_state_transform_table(in)%reaction_number,&
                       & hc_state_transform_table(in)%start_state,&
                       & hc_state_transform_table(in)%end_c_states(1),&
                       & hc_state_transform_table(in)%reaction_product_mass,&
                       & hc_state_transform_table(in)%charged_reaction_product_mass,&
                       & hc_state_transform_table(in)%reaction_desc(1:len_trim(hc_state_transform_table(in)%reaction_desc))
                                      
  end do


  write(6,*) 'Reaction table listing:'
  do in = 1,number_hc_species
     write(6,*) 'start state:',in,':',hc_state_table(in)%state_name
     write(6,*) 'reactions from state:',(hc_reaction_table(in)%reaction(i),i=1,hc_reaction_table(in)%number_reactions)
  end do   

  do in = 1,count_hc_reactions

     nt = size(hc_state_transform_table(in)%t_index)
     ne = size(hc_state_transform_table(in)%e_index)
     
     write(6,'(a,3i6,1x,a)') 'Process:',in,nt,ne,hc_state_transform_table(in)%reaction_desc

     do it = 1,nt
        write(6,'(2(1x,g18.8))') hc_state_transform_table(in)%t_index(it),&
               hc_state_transform_table(in)%reaction_data(it,ne)
     end do 
     
  end do
  



  write(6,*) 'Total reactions:', count_hc_reactions,count_e_reactions,count_p_reactions

  !stop "Debugging: End of jr data load"

 
 
  ! Some counters to get around the matrixes.
  !Integer :: i, j, k
  !Integer :: Io_Result = 0
 
  ! Check for unit number assignment.
  !Do While (Unit_open)
  !   Data_file_unit = Data_file_unit + 1
  !   Inquire (Unit = Data_file_unit, Opened = Unit_open)
  !End Do
 
  ! Open data file for reading.
  !Open (Unit = Data_file_unit, File = "hc_xsec_janev.dat", Status = 'Old')
  !If (Io_Result .ne. 0) then
  !   Write (Output_Unit_HC_Alert,*) 'ERROR: File not found ... Name = ','hc_xsec_janev.dat'
  !   Stop
  !End If
 
  ! Read through data file to the end.
  !Do While (Io_Result .eq. 0)
 
     ! Read a line from the data file.
     !Read (Unit = Data_file_unit, FMT = '(A200)',IOSTAT=Io_Result) Buffer
 
     ! Read until first character is not '!' or '#' or ' '.
     !Do While ((Buffer(1:1) .eq. '!') .or. (Buffer(1:1) .eq. ' ') .or. (Buffer(1:1) .eq. '#') .or. len_trim(Buffer).le.1)
     !   Read (Unit = Data_file_unit, FMT = '(A200)',IOSTAT=Io_Result) Buffer
     !End Do
 
     ! Interpret the line as a component of defined data structure.
     ! Check for Reaction label.
 
  !End Do
  !Close (Unit = Data_file_unit)
 
End Subroutine Load_HC_Data_Janev_reiter
 
Subroutine Load_HC_Data_Alman_Ruzic_Brooks ()
 
  ! Every good Fortran 90 program has...
  Implicit None
 
  ! Some counters to get around the matrixes.
  Integer :: i, j, k
  Integer :: Io_Result = 0
 
  ! Check for unit number assignment.
  Do While (Unit_open)
     Data_file_unit = Data_file_unit + 1
     Inquire (Unit = Data_file_unit, Opened = Unit_open)
  End Do
 
  ! Open data file for reading.
  Open (Unit = Data_file_unit, File = "hc_xsec_alman_ruzic_brooks.dat", Status = 'Old')
  If (Io_Result .ne. 0) then
     Write (Output_Unit_HC_Alert,*) 'ERROR: File not found ... Name = ','hc_xsec_alman_ruzic_brooks.dat'
     Stop
  End If
 
  ! Read through data file to the end.
  Do While (Io_Result .eq. 0)
 
     ! Read a line from the data file.
     Read (Unit = Data_file_unit, FMT = '(A200)',IOSTAT=Io_Result) Buffer
 
     ! Read until first character is not '!' or '#' or ' '.
     Do While ((Buffer(1:1) .eq. '!') .or. (Buffer(1:1) .eq. ' ') .or. (Buffer(1:1) .eq. '#') .or. len_trim(Buffer).le.1)
        Read (Unit = Data_file_unit, FMT = '(A200)',IOSTAT=Io_Result) Buffer
     End Do
 
     ! Interpret the line as a component of defined data structure.
     ! Check for Reaction label.
 
  End Do
  Close (Unit = Data_file_unit)
 
End Subroutine Load_HC_Data_Alman_Ruzic_Brooks
 
End Module HC_LdDta
