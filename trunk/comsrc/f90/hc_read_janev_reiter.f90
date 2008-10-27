module hc_read_janev_reiter

  Type process_Table_Type
     Integer :: Reaction_Number
     Character (len = 10) :: Start_state
     ! Note: when dealing with hydrocarbons higher than methane - the end_state quantity will need to be an array and
     !       the read code will need
     !       to be modified accordingly by scanning for all product states and recording them all. 
     character (len = 10) :: Eqn_inputs
     character (len = 40) :: Eqn_outputs
     character (len = 100) :: reaction_desc
     integer :: Reaction_Type
  End Type process_Table_Type


  integer,private,parameter :: max_energies    = 20
  integer,private,parameter :: max_temperature = 100
  integer,private,parameter :: max_reaction_string_len = 100
  character*7,private,parameter :: ProcessID = 'Process' 
  character*7,private,parameter :: ProzessID = 'Prozess' 
  character*7,private,parameter :: E_ID = '      E' 


  integer, private, save :: total_reactions_read 
  integer, private, save :: num_e_reactions,num_p_reactions

  integer,private :: jr_unit_num
  integer,private :: process_count,temperature_count,energy_count

  integer,private :: reaction_type
  integer,private,parameter :: e_reaction =1 
  integer,private,parameter :: p_reaction =2


  type(process_table_type),private, allocatable, save, dimension(:) :: process_table
  real*8,private,allocatable, save, dimension(:,:,:) :: reaction_data
  real*8,private,allocatable, save, dimension(:) :: t_index, e_index
  

  integer,private,parameter :: max_line_length = 500
  character*(max_line_length),private :: buffer
! slmod begin
!...bug
  character*6,parameter :: buff_fmt = '(a500)'
!
!  character*5,parameter :: buff_fmt = '(a500)'
! slmod end


contains


  subroutine init_load_jr
    implicit none
    !
    ! Initialization code for the Janev-Reiter database reading routines
    !
    total_reactions_read = 0
    num_e_reactions = 0
    num_p_reactions = 0

  end subroutine init_load_jr


  subroutine load_jr(filename,rc)
    implicit none
    character*(*) :: filename
    integer rc


    !
    ! Local variables
    !

    integer ios


    ! open file

    call find_free_unit_number(jr_unit_num)

    open(unit=jr_unit_num,file=filename,status='old',iostat=ios)

    if (ios.ne.0) then 
       write(0,*) 'ERROR Opening Janev-Reiter rate file:',filename
       rc = ios
       return
    endif

    ! analyse the input data file to determine data type and size
    call analyse_file


    ! Allocate temporary storage for reactions and reaction data
    call allocate_storage


    ! - add the reaction to the hc_reaction_table for the appropriate species
    ! - allocate storage and pointers for the indexing arrays and the data value 
    !   array
    ! - Second dimension is 1 for e reactions
    !
    ! - Note - the number of energy and temperature elements to the data set is not
    !          known a priori though this data could be added to the output file 
    !        - as a result, this code assumes a maximum of 200 temperature elements
    !          and 50 energy elements to the data set. Note that e- reactions only 
    !          have an electron temperature dependence
    ! 


    !
    ! read in and parse the reaction list  - storing each reaction to a 
    ! separate entry in the hc_state_transform_table
    !
    call read_process_data_block


    ! loop over processes - assign process number - start state and end state numbers
    ! assign whether it is an electron or proton reaction. 
    call read_sigmav_data


    !
    ! Clean up after reading in the data
    !
    call finished_reading_file


    ! 
    ! Assign data to the data structures used in the calling code
    !
    call assign_jr_data


    !
    ! Record totals for all processes read so far by this module
    !

    total_reactions_read = total_reactions_read + process_count

    if (reaction_type.eq.e_reaction) then 
       num_e_reactions = num_e_reactions + process_count
    elseif (reaction_type.eq.p_reaction) then 
       num_p_reactions = num_p_reactions + process_count
    endif

    !
    ! Deallocate local storage
    !
    call deallocate_storage


  end subroutine load_jr


  subroutine read_process_data_block
    implicit none

    logical :: reading_processes

    character*100 :: reaction_desc
    character*100 :: eqn_inputs
    character*100 :: eqn_outputs
    character*10 :: start_state_string
    integer :: proc_len
    integer :: process_num
    integer :: pos1,pos2,pos3,pos4
    integer :: ios

    integer current_process_count

    reading_processes = .false.
    proc_len = len_trim(ProcessID)

    do
       read(jr_unit_num,buff_fmt,iostat=ios) buffer

       if (ios.ne.0) then 
          exit
 
       elseif (buffer(1:proc_len).eq.ProcessID) then 

          if (.not.reading_processes) then 
             reading_processes = .true.
          endif

          !
          ! Interpret the process line
          !
          
          !
          ! Find the ":" so as to read the reaction description as 
          ! well as the process number
          !
          pos1 = proc_len
          pos2 = scan(buffer,':')

          read(buffer(pos1+1:pos2-1),*) process_num

          read(buffer(pos2+1:),'(a100)') reaction_desc

          ! Process reaction_desc into Eqn_inputs and Eqn_outputs

          pos3 = index(reaction_desc,'->') 
          eqn_inputs = reaction_desc(1:pos3-1)
          pos4 = index(reaction_desc,'Eopt=') 
          ! if reaction kinetic data is available then truncate the outputs where the Energy option 
          ! starts
          !write(0,'(2a)') reaction_desc(1:len_trim(reaction_desc)),':'
          if (pos4.eq.0) then 
             eqn_outputs = reaction_desc(pos3+2:)
             !write(0,'(a,2i6,a)') 'eqn_outputs:',pos3+2,pos4-1,eqn_outputs(1:len_trim(eqn_outputs))
          else
             eqn_outputs = reaction_desc(pos3+2:pos4-1)
             !write(0,'(a,2i6,a)') 'eqn_outputs:',pos3+2,pos4-1,eqn_outputs(1:len_trim(eqn_outputs))
          endif
          call mod_state_desc(eqn_inputs)
          call mod_state_desc(eqn_outputs)

          !
          ! Find start state
          !
          pos1 = scan(reaction_desc,'C') 
          pos2 = scan(reaction_desc(pos1+1:),' ')
          start_state_string = reaction_desc(pos1:pos1+pos2-1)

          !
          ! Substitute symbols and remove "_" from state description
          !
          call mod_state_desc(start_state_string)

          !
          ! Find end state - scan from end
          !
          !pos1 = scan(reaction_desc,'C',.true.) 
          !pos2 = scan(reaction_desc(pos1+1:),' ')
          !end_state_string = reaction_desc(pos1:pos1+pos2-1)
          !write(0,*) 'End  0:',end_state_string
          !
          ! Substitute symbols and remove "_" from state description
          !
          !call mod_desc_name(end_state_string)
          !write(0,*) 'End  1:',end_state_string


          ! Edit the inputs and outputs strings to match the conventions that Adam has
          ! previously implemented,
          ! In particular - strip +e from the end of outputs 
          !                 if outputs starts with H+ - change this to end with +H instead
          
          !call mod_outputs(eqn_outputs)


          process_table(process_num)%reaction_number = process_num
          process_table(process_num)%start_state = start_state_string
          !process_table(process_num)%end_state = end_state_string
          process_table(process_num)%eqn_inputs = eqn_inputs(1:len_trim(eqn_inputs))
          process_table(process_num)%eqn_outputs = eqn_outputs(1:len_trim(eqn_outputs))
          process_table(process_num)%reaction_type = reaction_type
          process_table(process_num)%reaction_desc = reaction_desc

       elseif (reading_processes) then 
          ! If the code is reading processes and finds a line which 
          ! does not contain the process keyword then the process block
          ! is finished
          exit
       endif

    end do 

  end subroutine read_process_data_block


  subroutine read_sigmav_data
    implicit none
    ! This routine reads in the sigmav data from the file into the data array

    logical :: read_finished
    integer :: in,cnt,ios
    integer :: processes_read
    integer :: proc_len,e_len

    ! Useful string lengths
    proc_len = len_trim(ProcessID)
    e_len = len_trim(E_ID)

    !
    ! Continue through the file to the reaction data 
    !
    ! 
    read_finished = .false.
    ios = 0

    if (reaction_type.eq.e_reaction) then 


       do

          read(jr_unit_num,buff_fmt,iostat=ios) buffer
          
          if(buffer(1:e_len).eq.E_ID) then 

             do cnt = 1,temperature_count

                read(jr_unit_num,buff_fmt,iostat=ios) buffer
                read(buffer,*) t_index(cnt),(reaction_data(in,cnt,1),in=1,process_count)

                !write(0,'(a,a)') 'Buff:',buffer(1:len_trim(buffer))
                !write(0,'(a,f10.5,2x,20g12.5)') 'Data:',t_index(cnt),(reaction_data(in,cnt,1),in=1,process_count)

             end do

             read_finished = .true.

          endif

          if (ios.ne.0.or.read_finished) exit

       end do


    elseif (reaction_type.eq.p_reaction) then 

       processes_read = 0

       do 

          read(jr_unit_num,buff_fmt,iostat=ios) buffer
          !write(0,'(a,a,a)') 'BUFF:',buffer(1:len_trim(buffer)),':'
          !
          ! Need to read one block for each process - look for ProcessID
          !
          if (buffer(1:proc_len).eq.ProcessID.or.buffer(1:proc_len).eq.ProzessID) then

             ! Discard 2 lines between process ID and reaction data 
             read(jr_unit_num,buff_fmt,iostat=ios) buffer
             read(jr_unit_num,buff_fmt,iostat=ios) buffer

             ! Read in energy indices - should be the same for each data set
             ! write(0,'(a,a,a)') 'EBUFF:',buffer(1:len_trim(buffer)),':'
             read(buffer,*) (e_index(in),in=1,energy_count)
             
             processes_read = processes_read+1
             
             ! Read in reaction data
             do cnt = 1,temperature_count
                read(jr_unit_num,buff_fmt,iostat=ios) buffer
                read(buffer,*) t_index(cnt),(reaction_data(processes_read,cnt,in),in = 1,energy_count) 
             end do

             if (processes_read.eq.process_count) then 
                read_finished = .true.
             endif

          end if

          if (ios.ne.0.or.read_finished) exit


       end do


    endif
       

  end subroutine read_sigmav_data


  subroutine finished_reading_file
    implicit none


    !
    ! Clean up after reading in the data - at the moment - just close the file
    !

    close(jr_unit_num)

  end subroutine finished_reading_file


  subroutine assign_jr_data

    use hc_init_lib_data
    implicit none
    
    !
    ! Local variables
    !
    integer :: in,it,ie
    integer :: process_number
    integer :: state,ierr

    !
    ! Reaction energy options
    !
    character*(5),parameter :: eopt_string = 'Eopt='
    character*(3),parameter :: ek_string   = 'Ek='
    integer :: eopt_loc,ek_loc
    integer :: eopt
    real    :: ek

    !
    ! This routine assigns the jr data read in from the file to the reaction table used in the rest of
    ! the HC module. 
    !
    ierr = 0
    !
    ! Loop through reactions from total_reactions_read+1 to total_reactions_read + process_count
    !
    ! Load data into the hc_state_transform_table and set values appropriately - allocate pointers and load data
    !

    do in = 1,process_count

       process_number = total_reactions_read + in

       if (in.eq.1) then 
          ! the temperature and energy index arrays need only be allocated once since they are repeated for
          ! all the data in the file - since this is the case - the pointers are set to point to the first allocated. 
          ! Note: these data are not deallocated at any time before the end of execution
          allocate(hc_state_transform_table(process_number)%t_index(temperature_count),stat=ierr)
          if (ierr.ne.0) write(0,*) 'Error allocating t_index:',process_number,temperature_count

          allocate(hc_state_transform_table(process_number)%e_index(energy_count),stat=ierr)
          if (ierr.ne.0) write(0,*) 'Error allocating e_index:',process_number,energy_count
          
          !
          ! Assign the electron temperature and proton energy dependencies for the cross sections to the newly allocated arrays
          !
          !do in = 1,temperature_count
          !   hc_state_transform_table(process_number)%t_index(in) = t_index(in)
          !enddo

          hc_state_transform_table(process_number)%t_index = real(t_index(:))
          hc_state_transform_table(process_number)%e_index = real(e_index(:))
             
       else
          !
          ! Assign pointers to the allocated index arrays
          !
          hc_state_transform_table(process_number)%t_index => hc_state_transform_table(total_reactions_read+1)%t_index
          hc_state_transform_table(process_number)%e_index => hc_state_transform_table(total_reactions_read+1)%e_index
       endif

       ! 
       ! Allocate the sigmav data storage and copy over the sigmav data from the local module storage
       !
       
       allocate(hc_state_transform_table(process_number)%reaction_data(temperature_count,energy_count),stat=ierr)
       if (ierr.ne.0) write(0,*) 'Error allocating reaction_data:',temperature_count,energy_count

       !
       ! Copy the cross-section of the array required into the allocated target
       !
       hc_state_transform_table(process_number)%reaction_data = reaction_data(in,:,:)

       !write(6,*) 'HC Reaction Data:',process_number,in
       !do it = 1,temperature_count
       !   write(6,'(f12.5,4(1x,g18.8)') hc_state_transform_table(process_number)%t_index(it),&
       !         hc_state_transform_table(process_number)%reaction_data(it,1),reaction_data(in,it,1)
       !end do 

       !
       ! Set data block flags
       !

       hc_state_transform_table(process_number)%reaction_number = process_number
        
       !
       ! Use existing code to assign start and end states
       !
       
       call parse_reaction(process_number,process_table(in)%eqn_inputs,process_table(in)%eqn_outputs)



       !HC_State_Transform_Table(process_number)%Reaction_Identifier = Data_Structure.Reaction
       !
       ! The following options are all related to the energy balance of particles in the reaction in
       ! terms of whether the hydrocarbon fragments or the hydrogen fragments pick up energy when the 
       ! process occurs and if so - how much. Set to zeroes as default values for now. 
       !
       !
       ! The reaction kinetic energies and options for their application have been added to the process
       ! lines of the JR data files. The options for application is preceded by the string Eopt= and the 
       ! energy value is preceded by the indicator Ek= - if these are not present values of zero are assigned.
       ! These data are extracted from the reaction_desc string of the local process table. 
       !

       if (index(process_table(in)%reaction_desc,eopt_string).ne.0) then 
          !
          ! reaction description includes energy and option data 
          ! 1) extract energy and option value
          !
          eopt_loc = index(process_table(in)%reaction_desc,eopt_string)
          ek_loc   = index(process_table(in)%reaction_desc,ek_string)
          read(process_table(in)%reaction_desc(eopt_loc+len_trim(eopt_string):),*) eopt
          read(process_table(in)%reaction_desc(ek_loc+len_trim(ek_string):),*) ek

          !write(0,'(a,i6,g12.5)') 'Ek:',eopt,ek

          !
          ! Note - set these values to 0.0 for now until the hc_freespace_transition code is properly implemented
          !
          ! eopt = 3 - energy divided among all products
          ! eopt = 4 - energy divided among charged products Ek_neut = 0.0
          !
          !if (eopt.eq.3.or.eopt.eq.4) then 
          !   HC_State_Transform_Table(process_number)%HC_E_Type=  0.0
          !   HC_State_Transform_Table(process_number)%HC_E = 0.0
          !else
          HC_State_Transform_Table(process_number)%HC_E_Type=  eopt
          HC_State_Transform_Table(process_number)%HC_E = ek
          !endif

       else
          HC_State_Transform_Table(process_number)%HC_E_Type=  0
          HC_State_Transform_Table(process_number)%HC_E = 0.0
          !
          HC_State_Transform_Table(process_number)%H_E_Type= 0
          HC_State_Transform_Table(process_number)%H_E = 0.0
       endif

       !
       ! Add a routine to load specific tabulated energy data for specific processes
       !

       !
       ! Moved the reaction_type determination to the parse reaction routine - based on leading e or p of the 
       ! start state string. 
       !
       !if (reaction_type.eq.e_reaction) then
       !   HC_State_Transform_Table(process_number) % Reaction_Type = 'e'
       !elseif (reaction_type.eq.p_reaction) then 
       !   HC_State_Transform_Table(process_number) % Reaction_Type = 'p'
       !endif

       HC_State_Transform_Table(process_number) % reaction_desc = process_table(in)%reaction_desc
       
       
       !
       ! Set quantities related to sigmav tables
       !
       ! Set sigmaV type to 3 - meaning table interpolation
       !
       HC_State_Transform_Table (process_number) % SigmaV_TPD = 3
       !
       ! Use t_index and e_index max and min to specify max and min temperatures and energies
       !
       
       HC_State_Transform_Table (process_number) % SigmaV_Tmin_Limit = t_index(1)
       HC_State_Transform_Table (process_number) % SigmaV_ValueEMin_Limit = reaction_data(in,1,1)
       HC_State_Transform_Table (process_number) % SigmaV_ValueEMax_Limit = reaction_data(in,temperature_count,energy_count)
       HC_State_Transform_Table (process_number) % SigmaV_Energy_Error_Limit = 0.0


       !
       ! Zero out sigma related entries in the transform table
       !
       HC_State_Transform_Table (process_number) % Sigma_TPD = 0
       HC_State_Transform_Table (process_number) % Sigma_Polynomial_Terms = 0
       HC_State_Transform_Table (process_number) % Sigma_Tmin_Limit = 0
       HC_State_Transform_Table (process_number) % Sigma_ValueEMin_Limit = 0
       HC_State_Transform_Table (process_number) % Sigma_ValueEMax_Limit = 0 
       HC_State_Transform_Table (process_number) % Sigma_Energy_Error_Limit = 0


       !
       ! Update reaction table entries
       !

       state = HC_species_Ident(process_table(in)%Eqn_Inputs (3:))

       !write(0,'(a,a,i6)') 'IDENT:',process_table(in)%Eqn_Inputs (3:),state

       ! Read reaction number into reaction table and increment number of reactions possible for that initial state by one.
       HC_Reaction_Table (state) % Number_Reactions = HC_Reaction_Table (state) % Number_Reactions + 1

       If ((HC_Reaction_Table (state) % Number_Reactions) > max_reactions_per_state) Then
          write(0,'(a,2i6)') 'Error:  Increase number of possible reactions for state (max_reactions_per_state: ',&
                                                &state,max_reactions_per_state
          stop
       End If
       ! Add current reaction number to reaction table.
       HC_Reaction_Table (state) % Reaction (HC_Reaction_Table (state) % Number_Reactions) = process_number
 

    end do


  end subroutine assign_jr_data




  subroutine mod_state_desc(state_string)
    implicit none
    character*(*) :: state_string

    integer :: strlen,in,in2,in_next

    !
    ! Remove _ ^ and space from state names - replace ^ with + if required
    !

    !
    ! May need to rewrite p reactions to place the H at the end ... i.e. start all OUTPUT reaction sections with C
    !

    strlen = len_trim(state_string)

    do in = 1,strlen

       if (state_string(in:in).eq.'^') then 
          ! Check to see if next character is a "+" 
          if (in.ne.strlen.and.state_string(in+1:in+1).eq."+") then 
             ! Delete the ^
             do in2 = in+1, strlen
                state_string(in2-1:in2-1) = state_string(in2:in2)
             end do
             ! Replace last character with a space
             state_string(strlen:strlen) = ' '
          else
             ! replace ^ with + for use in DIVIMP and compatibility with internal state names
             state_string(in:in) = '+'
          endif

       elseif (state_string(in:in).eq.'_') then 
          ! Delete underscores from the state name
          do in2 = in+1, strlen
             state_string(in2-1:in2-1) = state_string(in2:in2)
          end do
          ! Replace last character with a space
          state_string(strlen:strlen) = ' '
       elseif (state_string(in:in).eq.' '.and.in.ne.strlen) then 
           ! Delete the space by moving text up unless we are at the last character of the string
           ! find the next non-space character
           in_next = 0
           do in2 = in+1,strlen
              if (state_string(in2:in2).ne.' ') then
                 in_next = in2
                 exit
              endif
           end do 
           
           if (in_next.gt.0) then
              do in2 = in_next, strlen
                 state_string(in2-(in_next-in):in2-(in_next-in)) = state_string(in2:in2)
              end do
              ! Replace last character with a space
              do in2 = strlen-(in_next-in)+1,strlen
                 state_string(in2:in2) = ' '
              end do
           else
              ! If there are only spaces left in the string then the processing is completed
              exit

           endif
       endif

    end do

  end subroutine mod_state_desc

  subroutine mod_outputs(outputs)
    implicit none
    character*(*) :: outputs

    ! Edit the outputs string so it does not end with a +e (lists only H and C products)
    ! Edit the outputs string so it does not begin with H+ - change to +H at end 
    
    integer :: l_out,c_start,copy_len,in
    character*(100) :: out_copy
    
    l_out = len_trim(outputs)

    ! remove trailing +e

    if (outputs(l_out-1:l_out).eq.'+e') then 
       outputs(l_out-1:l_out) = '  '
    endif

    ! Move leading terms if any to end

    c_start = scan(outputs,'Cc')

    if (c_start.ne.1) then 
       out_copy = outputs(c_start:)
       copy_len = len_trim(out_copy)
       ! Add connection character
       out_copy(copy_len+1:copy_len+1) = '+'
       ! copy leading section up to just before the + before the C
       do in = 1,c_start-2
          out_copy(copy_len+1+in:copy_len+1+in) = outputs(in:in)
       end do 
       ! assign copy back 
       copy_len = len_trim(out_copy)
       outputs = out_copy(1:copy_len)
    end if

  end subroutine mod_outputs


  subroutine analyse_file
    implicit none
    !
    ! This routine scans through the input file to determine certain properties needed
    ! for storage allocation.
    ! 1) Number of processes in file
    ! 2) "e" or "p" processes
    ! 3) Number of temperature values in the table
    ! 4) Number of energy values in the cross-section tables if it is a "p" reaction
    !
    ! At the end it rewinds the file so it is ready for further processing
    ! It assumes that all "p" tables in the file have the same format as the first


    logical :: reading_processes,process_block_finished,analysis_finished
    character*100 :: reaction_desc
    integer :: proc_len,e_len
    integer :: pos1,pos2,pos3
    integer :: char_cnt
    integer :: in
    integer :: extstr,start_pos,strlen
    external extstr
    
    integer ios

    ! Initialization
    process_count = 0
    temperature_count = 0
    energy_count = 0

    ! Useful string lengths
    proc_len = len_trim(ProcessID)
    e_len = len_trim(E_ID)

    ! flags for read process
    reading_processes=.false.
    process_block_finished = .false.
    analysis_finished = .false.

    ! I/O error flag
    ios = 0


    do 

       ! Read in the file
       read(jr_unit_num,buff_fmt,iostat=ios) buffer

       !
       ! check for process block and read it in - determine if the file is e or p processes - assume
       ! only one type of process in each file. 
       !

       if ((.not.process_block_finished).and.buffer(1:proc_len).eq.ProcessID) then 

          process_count = process_count + 1

          if (process_count.eq.1) then 

             reading_processes = .true.

             pos1 = scan(buffer,':')
             read(buffer(pos1+1:),'(a100)') reaction_desc

             strlen = extstr(reaction_desc,start_pos)

             ! Determine type of processes in file
             
             if (reaction_desc(start_pos:start_pos).eq.'e') then
                reaction_type = e_reaction
             elseif (reaction_desc(start_pos:start_pos).eq.'p') then
                reaction_type = p_reaction
             endif

          endif
       elseif (reading_processes.and.(.not.process_block_finished)) then 
          ! Note finished reading the process block - this should be on the blank line after the process block
          process_block_finished = .true.

       elseif (process_block_finished.and.reaction_type.eq.e_reaction.and.buffer(1:e_len).eq.E_ID) then 
          ! found reaction table for e-reactions - should be process_count + 1 columns of data
          ! read data block and count lines until blank line or eof is reached
          temperature_count = 0

          do 
             read(jr_unit_num,buff_fmt,iostat=ios) buffer
             if (len_trim(buffer).eq.0.or.ios.ne.0) exit
             temperature_count = temperature_count +1

          end do

          energy_count = 1

          analysis_finished = .true.

       elseif (process_block_finished.and.reaction_type.eq.p_reaction.and.(buffer(1:proc_len).eq.ProcessID.or.&
              &buffer(1:proc_len).eq.ProzessID)) then 
          ! need to determine BOTH how many temperatures as well as how many energies
          
          ! read in and discard blank line 
          read(jr_unit_num,buff_fmt,iostat=ios) buffer

          ! Read in line with energies listed
          read(jr_unit_num,buff_fmt,iostat=ios) buffer

          ! read in energies into real array

          ! Need to determine how many energy values are listed - loop through counting blocks of data on line
          ! Assume that the line will contain only numbers, +, -, ., and E/e or D/d and spaces
          ! Data is separated by one or more white space
          
          char_cnt = 0
          energy_count = 0

          do
             pos1 = char_cnt + scan(buffer(char_cnt+1:),'0123456789+-.EeDd')

             if (pos1.eq.char_cnt) then 
                exit
             else
                energy_count = energy_count + 1
                ! find end of item
                char_cnt = pos1 + scan(buffer(pos1+1:),' ')
                if (char_cnt.eq.pos1) then
                   ! if there are no more separators then there are no more data values
                   exit
                endif
             endif

          end do 

          ! Now read in the rest of the table to determine number of temperatures - assume that all the p-reactions in the 
          ! file will have the same number of lines of temperature data

          temperature_count = 0

          do 
             read(jr_unit_num,buff_fmt,iostat=ios) buffer
             if (len_trim(buffer).eq.0.or.ios.ne.0) exit
             temperature_count = temperature_count +1

          end do

          analysis_finished = .true.

       endif
       

       if (ios.ne.0.or.analysis_finished) exit


     end do

     ! restore input file

     rewind(jr_unit_num)

  end subroutine analyse_file



  subroutine allocate_storage
    implicit none
    ! This routine allocates storage for the reaction processes that are read in as well as the reaction data 
    integer ierr

    ierr = 0 

    allocate(reaction_data(process_count,temperature_count,energy_count),stat=ierr)
    if (ierr.ne.0) write(0,*) 'Error allocating reaction data:'

    allocate(t_index(temperature_count),stat=ierr)
    if (ierr.ne.0) write(0,*) 'Error allocating t-index:'

    allocate(e_index(energy_count),stat=ierr)
    if (ierr.ne.0) write(0,*) 'Error allocating e-index:'

    ! Allocate storage to store process information

    allocate(process_table(process_count),stat=ierr)
    if (ierr.ne.0) write(0,*) 'Error allocating process table:'

    ! Initialize allocated numerical storage
    reaction_data = 0.0
    t_index = 0.0
    e_index = 0.0

  end subroutine allocate_storage


  subroutine deallocate_storage
    implicit none
    ! This routine deallocates storage for the reaction processes that are read in as well as the reaction data 

    deallocate(reaction_data)
    deallocate(t_index)
    deallocate(e_index)

    ! DeAllocate storage to store process information

    deallocate(process_table)

  end subroutine deallocate_storage

  subroutine set_jr_reaction_counts(count_hc_reactions,count_e_reactions,count_p_reactions)
    implicit none
    integer :: count_hc_reactions,count_e_reactions,count_p_reactions

    count_hc_reactions=total_reactions_read
    count_e_reactions=num_e_reactions
    count_p_reactions=num_p_reactions

  end subroutine set_jr_reaction_counts



end module hc_read_janev_reiter
