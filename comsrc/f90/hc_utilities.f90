! -*-Mode:f90-*-
! HC_Utilities.f90
! Contains generic hydrocarbon utilities including:
! -HC_Species_Ident: Returns hydrocarbon numeric identifier with its species name.
! -HC_Ident_Species: Returns hydrocarbon species name with its numeric identifier.
!
! Adam McLean under Prof. Peter Stangeby
! Research Assistant: David Elder
! July, 2002
 
Module HC_Utilities
 
  use error_handling
  ! Every good Fortran program has...
  Implicit None	
 
Contains
 
  Real Function Calc_Sputy (HC_Mass)
    ! Returns the number of carbons in the hydrocarbon molecule given its mass.
 
    ! Required modules.
    Use ComHC ! Contains Output_Unit_HC_Alert.
 
    ! Every good Fortran program has...
    Implicit None
 
    ! Declare input/output variables.
    Real, Intent (In) :: HC_Mass
 
    ! jdemod
    !
    ! This code will not work correctly for CT4 or any heavier hydrocarbons with a mix of species. This routine should be
    ! modified to return the number of carbons associated with the current hc species. 
    ! 
    ! i.e. sputy = hc_state_table(cur_hc_spec)%carbon_content * 1.0
    !
    ! However, even this is a bit of an open question since I think a C3H8 should still only have a weight of one for the molecule - even
    ! though it will split into 3 equally weighted C's at some point.
    !
    !

    If (HC_Mass .lt. 12) Then

       ! Error, this should never happen.

       call errmsg("Error in Calc_Sputy: Hydrocarbon with mass less than 12 AMU should never occur:",HC_Mass)

       Write (Output_Unit_HC_Alert,*) "Error in Calc_Sputy: Hydrocarbon with mass less than 12 AMU should never occur:",HC_Mass
       Write (Output_Unit_HC_Alert,*) "Program stopping."
       Stop
    ElseIf (HC_Mass .ge. 12 .and. HC_Mass .lt. 24) Then
       ! Single carbon atom in the hydrocarbon.
       Calc_Sputy = 1
    ElseIf (HC_Mass .ge. 24 .and. HC_Mass .lt. 36) Then
       ! C2 molecule.
       Calc_Sputy = 2
    ElseIf (HC_Mass .ge. 36 .and. HC_Mass .lt. 48) Then
       ! C3 molecule.
       Calc_Sputy = 3
    Else
       ! Error, this should never happen.
       Write (Output_Unit_HC_Alert,*) "Error in Calc_Sputy: Hydrocarbon with mass greater than 48 AMU is not supported at this &
&time:",HC_Mass
       Write (Output_Unit_HC_Alert,*) "Program stopping."
       Stop
    End If
 
  End Function Calc_Sputy

 Real Function Calc_Sputy2 (cur_hc_spec)
    ! Returns the number of carbons in the hydrocarbon molecule given its mass.
 
    ! Required modules.
    Use ComHC ! Contains Output_Unit_HC_Alert.
    use hc_init_lib_data
 
    ! Every good Fortran program has...
    Implicit None
 
    ! Declare input/output variables.
    integer, Intent (In) :: cur_hc_spec
 
    !
    ! Revised calc_sputy which is based on the number of C in the species
    !

    calc_sputy2 = hc_state_table(cur_hc_spec)%carbon_content
 
  End Function Calc_Sputy2
 

 
  Real Function Find_Wall_Temperature (Wall_Segment)
 
    ! Find surface temperature (in K) for use with models for hydrocarbon energy/species release.
 
    Use HC_Get ! Access vessel geometry data.
 
    ! Every good Fortran 90 program has...		
    Implicit None
 
    ! Declare input/output variables.
    Integer, Intent (In) :: Wall_Segment
 
    ! Grab surface temperature from DIVIMP WALLPT array (ind,19)
    Find_Wall_Temperature = gwallpt (Wall_Segment,19)
 
  End Function Find_Wall_Temperature
 
  Real Function Find_Target_Temperature (Target_Segment)
 
    ! Find surface temperature (in K) for use with models for hydrocarbon energy/species release.
 
    Use HC_Get ! Access vessel geometry data.
 
    ! Every good Fortran 90 program has...		
    Implicit None
 
    ! Declare input/output variables.
    Integer, Intent (In) :: Target_Segment
 
    ! Grab surface temperature from DIVIMP TEMDS array (maxnds).
    Find_Target_Temperature = gtempds (Target_Segment)
 
  End Function Find_Target_Temperature
 
  Subroutine HC_Starting_Position (IProd,LProd,NeutType,R_Launch,Z_Launch,HC_Cell,HC_Ring,Segment_Index,HC_Dir_Indicator,Launch_Reg)
 
    ! Determine starting position based on DIVIMP input options.
 
    Use ComHC
    Use HC_Get
    Use HC_Init_DIV_Data
 
    ! Every good Fortran 90 program has...
    Implicit None
 
    ! Declare input/output variables.
    Integer, Intent (In) :: IProd
    Integer, Intent (In) :: LProd
    Integer, Intent (InOut) :: NeutType
    Real, Intent (Out) :: R_Launch
    Real, Intent (Out) :: Z_Launch
    Integer, Intent (Out) :: HC_Cell ! IK
    Integer, Intent (Out) :: HC_Ring ! IR
    Integer, Intent (Out) :: Segment_Index ! ID
    Integer, Intent (Out) :: HC_Dir_Indicator ! IS
    Integer, Intent (Out) :: Launch_Reg ! M
 
    ! Declare local variables.
    Integer :: Launch_Option
    Integer :: Grid_Option
 
    If (IProd + LProd -1 .le. Max_Impurities) Then
       R_Launch = gxprods (IProd + LProd - 1)
       Z_Launch = gyprods (IProd + LProd - 1)
    End If

    !write(0,*) 'INIT:',r_launch,z_launch
    !write(output_unit_scratch,*) 'INIT:',r_launch,z_launch
    
 
    ! Set IK, IR to 0,0.
    HC_Cell = 0
    HC_Ring = 0
 
    ! Assign local variables.
    Launch_Option = hc_launch_location ! CNEUTB.
    Grid_Option =  Grid_Opt
 
    If (Launch_Option .eq. 0 .or. Launch_Option .eq. 3) Then
       ! TARGET launch.
       ! ik, ir - indices of nearest bin centre.
       ! id - target index from which the launch is occuring.
       HC_Dir_Indicator = 0
       Segment_Index = gidprods (IProd + LProd - 1) ! ID
       HC_Cell = gikds (Segment_Index) ! IK
       HC_Ring = girds (Segment_Index) ! IR
       Grid_Error = .false.
 
       ! Find which target M to source from.
       ! Check for
       If (Grid_Option .eq. 0 .or. Grid_Option .eq. 1 .or. Grid_Option .eq. 3) Then
 
          ! Assume M = 2.
          If (Launch_Option .eq. 0) Then
             If (NeutType .ne. 3) Then
                ! Primary launch event.
                Launch_Reg = HC_Launch_Reg_Target_2_Dist
             Else
                ! Target self-sputtering event.
                Launch_Reg = HC_Launch_Reg_Sput_Target_2
             End If
          ElseIf (Launch_Option .eq. 3) Then
             If (NeutType .ne. 3) Then
                ! Primary launch event.
                Launch_Reg = HC_Launch_Reg_Target_2_Pin
             Else
                ! Target self-sputtering event.
                Launch_Reg = HC_Launch_Reg_Sput_Target_2
             End If
          End If
 
          ! Check to see if M not 2.
          If ((Launch_Option .eq. 0) .and. Segment_Index .le.  Inner_Target_Points) Then
             ! Distributed target launch.
             If (NeutType .ne. 3) Then
                ! Primary launch event.
                Launch_Reg = HC_Launch_Reg_Target_1_Dist
             Else
                ! Target self-sputtering event.
                Launch_Reg = HC_Launch_Reg_Sput_Target_1
             End If
          ElseIf ((Launch_Option .eq. 3) .and. Segment_Index .le.  Inner_Target_Points) Then
             ! EIRENE-distributed target launch.
             If (NeutType .ne. 3) Then
                ! Primary launch event.
                Launch_Reg = HC_Launch_Reg_Target_1_Pin
             Else
                ! Target self-sputtering event.
                Launch_Reg = HC_Launch_Reg_Sput_Target_1
             End If
          End If
 
       ElseIf (Grid_Option .eq. 2) Then
 
          ! Assume M = 2.
          If (Launch_Option .eq. 0) Then
             If (NeutType .ne. 3) Then
                ! Primary launch event.
                Launch_Reg = HC_Launch_Reg_Target_2_Dist
             Else
                ! Target self-sputtering event.
                Launch_Reg = HC_Launch_Reg_Sput_Target_2
             End If
          ElseIf (Launch_Option .eq. 3) Then
             If (NeutType .ne. 3) Then
                ! Primary launch event.
                Launch_Reg = HC_Launch_Reg_Target_2_Pin
             Else
                ! Target self-sputtering event.
                Launch_Reg = HC_Launch_Reg_Sput_Target_2
             End If
          End If
 
          ! Check to see if M not 2.
          If((Launch_Option.eq.0).and.(Segment_Index.gt.Inner_Target_Points.and.Segment_Index.le.End_Target_Point))Then
             ! Distributed target launch.
             If (NeutType .ne. 3) Then
                ! Primary launch event.
                Launch_Reg = HC_Launch_Reg_Target_1_Dist
             Else
                ! Target self-sputtering event.
                Launch_Reg = HC_Launch_Reg_Sput_Target_1
             End If
          ElseIf((Launch_Option.eq.3).and.(Segment_Index.gt.Inner_Target_Points.and.Segment_Index.le.End_Target_Point))Then
             ! EIRENE-distributed target launch.
             If (NeutType .ne. 3) Then
                ! Primary launch event.
                Launch_Reg = HC_Launch_Reg_Target_1_Pin
             Else
                ! Target self-sputtering event.
                Launch_Reg = HC_Launch_Reg_Sput_Target_1
             End If
          End If
       End If
 
    ElseIf (Launch_Option .eq. 1 .or. Launch_Option .eq. 5) Then
       ! Free space launch.
       ! IK, IR are the indices of the closest bin centre.
       ! ID, IS contain the value 0 and are not used.
       Segment_Index = gidprods (IProd + LProd - 1)
       HC_Dir_Indicator = 0
       Call gridpos (HC_Cell, HC_Ring, R_Launch, Z_Launch, .true.,  Grid_Error)
 
       ! Set launch region for statistics collection.
       If (Launch_Option .eq. 1) Then
          If (NeutType .ne. 3) Then
             Launch_Reg = HC_Launch_Reg_Free_Space_PT
          Else
             Launch_Reg = HC_Launch_Reg_Sput_FS_PT
          End If
       ElseIf (Launch_Option .eq. 5) Then
          If (NeutType .ne. 3) Then
             ! jdemod - I think this should be launch_reg
             !Launch_Option = HC_Launch_Reg_Free_Space_2D
             Launch_Reg = HC_Launch_Reg_Free_Space_2D
          Else
             ! jdemod - I think this should be launch_reg
             Launch_Reg = HC_Launch_Reg_Sput_FS_2D
          End If
       End If

    ElseIf (Launch_Option .eq. 2 .or. Launch_Option .eq. 4) Then
       ! WALL launch.
       Segment_Index = gidprods (IProd + LProd - 1)
       HC_Dir_Indicator = gisprods (IProd + LProd - 1)
       Call gridpos (HC_Cell, HC_Ring, R_Launch, Z_Launch, .true.,  Grid_Error)
 
       ! Set launch region for statistics collection.
       If (Launch_Option .eq. 2) Then
          If (NeutType .ne. 3) Then
             Launch_Reg = HC_Launch_Reg_Wall_Homo
          Else
             Launch_Reg = HC_Launch_Reg_Sput_Wall
          End If
       ElseIf (Launch_Option .eq. 4) Then
          If (NeutType .ne. 3) Then
             ! jdemod - I think this should be launch_reg
             !Launch_Option = HC_Launch_Reg_Wall_Dist
             Launch_Reg = HC_Launch_Reg_Wall_Dist
          Else
             ! jdemod - I think this should be launch_reg
             Launch_Reg = HC_Launch_Reg_Sput_Wall
          End If
       End If
 
    ElseIf (Launch_Option .eq. 6) Then
       ! Start from a single WALL or TARGET element - index needs to be appropriately assigned 
       ! Note that the issue with detection of particles on cell boundaries needs to be resolved. 
       Segment_Index = HC_launch6_wall_index
       HC_Dir_Indicator = 0 ! Should be 0/1 equally in random distribution.
       !R_Launch = gwallpt (gwallindex (Segment_Index), 1)
       !Z_Launch = gwallpt (gwallindex (Segment_Index), 2)
       R_Launch = gwallpt (Segment_Index, 1)
       Z_Launch = gwallpt (Segment_Index, 2)
       Call gridpos (HC_Cell,HC_Ring,R_Launch,Z_Launch,.true., Grid_Error)
 
       ! Adjust NeutType for target/wall location.
       ! Assume target 1 for now.
       !write (0,*) "aherenow:",NeutType,HC_Ring,HC_Cell,R_Launch,Z_Launch,gidds (HC_Ring,1), Grid_Error
! slmod begin
       If ( .not.Grid_Error ) Then
!
!       If ( Grid_Error .ne. .true.) Then
! slmod end
          ! Starting on a target.
          ! Change neuttype to indicate a target launch.
          If (NeutType .eq. 4) Then
             NeutType = 1
          ElseIf (NeutType .eq. 5) Then
             NeutType = 2
          End If
 
          ! Set launch region for statistics collection.
          If (NeutType .ne. 3) Then
             If (gwallpt (Segment_Index,16) .eq. 4) Then
                ! Outer target
                Launch_Reg = HC_Launch_Reg_Target_1_Dist
             ElseIf (gwallpt (Segment_Index,16) .eq. 1) Then
                ! Index target
                Launch_Reg = HC_Launch_Reg_Target_2_Dist
             End If
          Else
             If (gwallpt (Segment_Index,16) .eq. 4) Then
                ! Outer target
                Launch_Reg = HC_Launch_Reg_Sput_Target_1
             ElseIf (gwallpt (Segment_Index,16) .eq. 1) Then
                ! Index target
                Launch_Reg = HC_Launch_Reg_Sput_Target_2
             End If
          End If
       Else
          ! Starting on a wall.
          ! Change neuttype to indicate a wall launch.
          If (NeutType .eq. 1) Then
             NeutType = 4
          ElseIf (NeutType .eq. 2) Then
             NeutType = 5
          End If
 
          ! Set launch region for statistics collection.
          If (NeutType .ne. 3) Then
             Launch_Reg = HC_Launch_Reg_Wall_Dist
          Else
             Launch_Reg = HC_Launch_Reg_Sput_Wall
          End If
       End If
 
    End If
 
    ! Assign starting index for either target or wall coordinate of launch.
    Starting_Index = Segment_Index ! Saved as either target or wall index.

    
    !
    ! This should not be here - starting index should be set correctly above for wall/target launches
    !
    !If (Launch_Option .eq. 6) Then
    !   ! Assume target 1 for now.
    !   Starting_Index = gidds (HC_Ring,1)
    !End If
 
    ! Assign starting index for wall.
 
    ! Krieger IPP/07 - SUN compiler insists on 132 column limit
    if (debug_hc) & 
         &  write (6,'(a,2i8,2g12.5,3i8)') &
           & "HC_DEBUG: NEUTRAL INITIAL POSITION",NeutType,Segment_Index,R_Launch,Z_Launch,HC_Cell,HC_Ring, &
           & Inner_Target_Points
 
    If (NeutType .eq. 1 .or. NeutType .eq. 2 .or. NeutType .eq. 3) Then
       Launch_Wall_Index = gwallindex ( Starting_Index)
    ElseIf (NeutType .eq. 4 .or. NeutType .eq. 5) Then
       Launch_Wall_Index =  Starting_Index
    Else
       Launch_Wall_Index = -1
    End If
 
    ! Record use of determined launch region for output.
    HC_Region_Used (Launch_Reg) =  HC_Region_Used (Launch_Reg) + 1
 
    !write (0,*) "endstuff",Segment_Index, Starting_Index,&
    !&  Launch_Wall_Index,Launch_Reg,Grid_Option,Launch_Option,R_Launch, &
    !&  R_Plasma_Centre, Lower_Trap_Ring, &
    !&  Inner_SOL_Ring,gkss (HC_Cell, HC_Ring),gksmaxs (HC_Ring) / 2.0, &
    !&  Num_Upper_Rings, End_Target_Point
    !write (0,*) "POS Data",Launch_option, Grid_option, Neuttype, Launch_Reg
 
  End Subroutine HC_Starting_Position
 
  Subroutine HC_Injection_Position (Seed,NRand,R_Launch,Z_Launch,Current_S,Current_Cross, &
       & NeutType,STmp,HC_Cell,HC_Ring,Segment_Index,HC_Dir_Indicator,Launch_Reg)
 
    ! Determine starting position based on DIVIMP input options.
 
    Use ComHC
    Use HC_Init_DIV_Data
    Use HC_Get ! Contains gcxsc,gcysc
 
    ! Every good Fortran 90 program has...
    Implicit None
 
    ! Declare input/output variables.
    Double Precision, Intent (In) :: Seed ! Use for new random numbers.
    Integer, Intent (InOut) :: NRand ! NRAND		
    Real, Intent (Out) :: R_Launch ! R
    Real, Intent (Out) :: Z_Launch ! Z
    Real, Intent (Out) :: Current_S ! S
    Real, Intent (Out) :: Current_Cross ! Cross
    Integer, Intent (InOut) :: NeutType		
    Real, Intent (Out) :: STmp ! Calculated S.
    Integer, Intent (Out) :: HC_Cell ! IK
    Integer, Intent (Out) :: HC_Ring ! IR
    Integer, Intent (Out) :: Segment_Index ! ID
    Integer, Intent (Out) :: HC_Dir_Indicator ! IS
    Integer, Intent (Out) :: Launch_Reg ! M
 
    ! Declare local variables.
    Real, Dimension ( Max_Cells_Per_Ring *  Max_Rings) :: Local_InjProb
    Real :: Random_Value
    Integer :: Injection_Bin
    Integer :: Index_Value
 
    Real, External :: IPOS
    STmp = 0.0
 
    ! Set launch region for statistics collection.
    Launch_Reg = HC_Launch_Reg_Free_Space_2D
 
    If ( Injection_Opt .eq. 1) Then ! CIOPTE
       R_Launch = gcxsc ()
       Z_Launch = gcysc ()
 
       ! Assign starting cell and ring.
       Call gridpos (HC_Cell,HC_Ring,R_Launch,Z_Launch,.true., Grid_Error)
 
       ! Set launch region for statistics collection.
       Launch_Reg = HC_Launch_Reg_Free_Space_PT
 
       Launch_Cell = HC_Cell
       Launch_Ring = HC_Ring
 
    ElseIf ( Injection_Opt .eq. 2 .or.  Injection_Opt .eq. 3 .or.  Injection_Opt .eq. 5 .or.  Injection_Opt .eq. 6) Then
       ! FOR NOW SIMPLY FIND GRID POINT CLOSEST TO DESIRED
       ! INJECTION POSITION. THIS IS ALL THAT IS CURRENTLY DONE.
       ! THE PARTIAL DISPLACEMENTS FROM THE GRID POINTS ARE IGNORED.
 
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value)
       STmp = Random_Value * ( INJ_Area_Upper_Bound -  INJ_Area_Lower_Bound) * gksmaxs ( INJ_Ring_Number) +  INJ_Area_Lower_Bound *&
& gksmaxs ( INJ_Ring_Number)
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value)
       If (Random_Value .gt. 0.5 .and. ( Injection_Opt .eq. 2 .or.  Injection_Opt .eq. 5)) Then
 
          ! OUTER PLATE - OTHERWISE INNER - for option 2
          STmp = gksmaxs ( INJ_Ring_Number) - STmp
       End If
 
       ! FIND R,Z COORDINATES OF NEAREST GRID POINT
       ! FIND NEAREST IK CORRESPONDING TO DISTANCE S ALONG CONTOUR INJI
       HC_Cell = 1
       Do
          If (HC_Cell .lt. gnks ( INJ_Ring_Number) .and. STmp .gt. gkss (HC_Cell, INJ_Ring_Number)) Then
             HC_Cell = HC_Cell + 1
          Else
             Exit
          End If
       End Do
 
       Do
          If (HC_Cell .gt. 1 .and. STmp .le. gkss (HC_Cell - 1,  INJ_Ring_Number)) Then
             HC_Cell = HC_Cell - 1
          Else
             Exit
          End If
       End Do
 
       If (HC_Cell .gt. 1 .and. (STmp - gkss (HC_Cell - 1,  INJ_Ring_Number) .lt. gkss (HC_Cell,  INJ_Ring_Number) - STmp)) Then
          HC_Cell = HC_Cell - 1
       End If
       !HC_Cell = 34
       HC_Ring =  INJ_Ring_Number
 
       Launch_Cell = HC_Cell
       Launch_Ring = HC_Ring
 
       If ( RZ_Opt .eq. 1) Then
          Call GETRZ (HC_Cell,HC_Ring,STmp,Current_Cross,R_Launch,Z_Launch, RZ_Opt)
       Else
          R_Launch = grs (HC_Cell,HC_Ring)					
          Z_Launch = gzs (HC_Cell,HC_Ring)					
       End If
 
       !R_Launch = grs (HC_Cell,HC_Ring)
       !Z_Launch = gzs (HC_Cell,HC_Ring)
 
       Write (6,*) 'HC stmp:',STmp,HC_Cell,HC_Ring, INJ_Ring_Number,gkss (HC_Cell, INJ_Ring_Number),R_Launch,Z_Launch,gksmaxs ( &
&INJ_Ring_Number)
       !Write (0,*) 'HC stmp:',STmp,HC_Cell,HC_Ring, INJ_Ring_Number,gkss (HC_Cell, INJ_Ring_Number),R_Launch,Z_Launch,gksmaxs ( INJ_Ring_Number)
 
    ElseIf ( Injection_Opt .eq. 4 .or.  Injection_Opt .eq. 7) Then
       ! Find the initial injection bin
       ! Note: this option works only if DIVIMP Control Option is set to 1 to launch ions.
 
       ! Read all external injprob values into local injprob
       Index_Value = 1
       Do
          Local_InjProb (Index_Value) = ginjprob (Index_Value)
          Index_Value = Index_Value + 1
          If (Index_Value .gt.  Max_Cells_Per_Ring *  Max_Rings) Then
             Exit
          End If
       End Do
 
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Random_Value)
       !write (0,*) "Data:",Random_Value, Number_Injected_Particles,SUM(Local_InjProb)
 
       Injection_Bin = IPOS (Random_Value,Local_InjProb, Number_Injected_Particles)
       HC_Cell = ginjkind (Injection_Bin)
       HC_Ring = ginjrind (Injection_Bin)
 
       Launch_Cell = HC_Cell
       Launch_Ring = HC_Ring
 
       R_Launch = grs (HC_Cell, HC_Ring)
       Z_Launch = gzs (HC_Cell, HC_Ring)
 
    End If
 
    HC_Dir_Indicator = 0.0
    Segment_Index = 1
 
    ! Assign starting index for either target or wall coordinate of launch.
    Starting_Index = Segment_Index ! Saved as either target or wall index.
    If (hc_launch_location .eq. 6) Then
       ! Assume target 1 for now.
       Starting_Index = gidds (HC_Ring,1)
    End If
 
    ! Assign starting index for wall.
    !write (0,*) "stuff",NeutType,Segment_Index,R_Launch,Z_Launch,HC_Cell,HC_Ring, Inner_Target_Points, Starting_Index,gwallindex ( Starting_Index)
    If (NeutType .eq. 1 .or. NeutType .eq. 2 .or. NeutType .eq. 3) Then
       Launch_Wall_Index = gwallindex ( Starting_Index)
       Launch_Wall_Index = 1
    ElseIf (NeutType .eq. 4 .or. NeutType .eq. 5) Then
       Launch_Wall_Index =  Starting_Index
    Else
       Launch_Wall_Index = -1
    End If
 
    ! Record use of determined launch region for output.
    HC_Region_Used (Launch_Reg) =  HC_Region_Used (Launch_Reg) + 1
 
    !write (0,*) "endstuff",Segment_Index, Starting_Index,&
    !&  Launch_Wall_Index,Launch_Reg,Grid_Option,Launch_Option,R_Launch, &
    !&  R_Plasma_Centre, Lower_Trap_Ring, &
    !&  Inner_SOL_Ring,gkss (HC_Cell, HC_Ring),gksmaxs (HC_Ring) / 2.0, &
    !&  Num_Upper_Rings, End_Target_Point
    !write (0,*) "POS Data",Launch_option, Grid_option, Neuttype, Launch_Reg
 
  End Subroutine HC_Injection_Position
 
  Subroutine HC_Launch_Angle (IProd,Beta,Psi,Launch_Angle,Seed,NRand)
 
    ! Hydrocarbon launch angle with all available options from the launch subroutine in neut.d6a.
    ! Finds a random value from -90 to + 90 degrees (-pi/2 to pi/2).
 
    Use ComHC
    Use HC_Init_DIV_Data
    Use HC_Get
 
    ! Every good Fortran 90 program has...
    Implicit None
 
    ! Declare input/output variables.
    Integer, Intent (In) :: IProd
    Real, Intent (Out) :: Beta
    Real, Intent (Out) :: Psi
    Real, Intent (Out) :: Launch_Angle
    Double Precision, Intent (In) :: Seed ! Use for new random numbers.
    Integer, Intent (InOut) :: NRand ! NRAND		
 
    ! Declare local variables.
    Integer :: Velocity_Angle_Option
    Real :: Random_Value_1
    Real :: Random_Value_2
    real :: temp
    save random_value_1,random_value_2,temp
    
    ! Initialize output variables.
    Launch_Angle = 0.0
    Beta = 0.0
    Psi = 0.0
 
    ! Assign local variables from modules.
    Velocity_Angle_Option = hc_launch_angle_velocity ! Set to CNEUTC or user specified for independent HC use.
 
    ! Changed all access to RANV arrays to local random call.
    NRand = NRand + 1
    Call Surand2 (Seed, 1, Random_Value_1)
    NRand = NRand + 1
    Call Surand2 (Seed, 1, Random_Value_2)
 
    If (Velocity_Angle_Option .eq. 0) Then
       Launch_Angle = SIGN (ASIN (Random_Value_1), Random_Value_2-0.5)
    ElseIf (Velocity_Angle_Option .eq. 3 .or. Velocity_Angle_Option .eq. 5 .or. Velocity_Angle_Option .eq. 9) Then
       ! Note, this is recommended for most HC launches.
       !Launch_Angle = SIGN (ASIN (SQRT(Random_Value_1)), Random_Value_2-0.5)
       !write (0,'(a,i6,2g12.5,i6,5g12.5)') "Launch angle is:",Velocity_Angle_Option,Beta,Psi,IProd,Launch_Angle,&
       !     & random_value_1,random_value_2,sqrt(random_value_1),asin(sqrt(random_value_1))
       !
       ! NOTE: The following repeat of the same two identical lines is a workaround for a PGI compiler bug of some
       !       sort - still in the process of tracking it down - jde - May 15 2007
       !
       !Beta = ASIN (sqrt(random_value_1))
       !write(0,*) 'beta1:',beta,random_value_1

       Beta = ASIN (sqrt(random_value_1))
       !write(0,*) 'beta2:',beta,random_value_1
       Psi = 2.0 *  Pi_Value * Random_Value_2

       Launch_Angle = ATAN (TAN (Beta) * COS (Psi))

       !write (6,'(a,i6,2g12.5,i6,5g12.5)') "Launch angle is:",Velocity_Angle_Option,Beta,Psi,IProd,Launch_Angle,&
       !  & random_value_1,random_value_2,sqrt(random_value_1),asin(sqrt(random_value_1))

    ElseIf (Velocity_Angle_Option .eq. 1 .or. Velocity_Angle_Option .eq. 2 .or. Velocity_Angle_Option .eq.4) Then
       Beta = ASIN (SQRT(Random_Value_1))
       Psi = 2.0 *  Pi_Value * Random_Value_2
       Launch_Angle = ATAN (TAN (Beta) * COS (Psi))
    ElseIf (Velocity_Angle_Option .eq. 6) Then
       Launch_Angle = 0.0
    ElseIf (Velocity_Angle_Option .eq. 7 .or. Velocity_Angle_Option .eq. 11) Then
       Launch_Angle = SIGN (ACOS ((1.0-Random_Value_1)**(1.0/3.0)),Random_Value_2-0.5)
    ElseIf (Velocity_Angle_Option .eq. 8 .or. Velocity_Angle_Option .eq. 15) Then
       Launch_Angle = 2.0 *  Pi_Value * Random_Value_1 -  Pi_Value
    ElseIf (Velocity_Angle_Option .eq. 10) Then
       Beta = ACOS ((1.0 - Random_Value_1) ** (1.0/3.0))
       PSI = 2.0 *  Pi_Value * Random_Value_2
       Launch_Angle = Beta
       If (Psi .lt.  Pi_Value / 2.0 .or. Psi .gt. 3.0 *  Pi_Value / 2.0) Then
          Launch_Angle = -Beta
       End If
    ElseIf (Velocity_Angle_Option .eq. 12 .or. Velocity_Angle_Option .eq. 13) Then
       Beta = ACOS ((1.0 - Random_Value_1)**(1.0 / ( Cosine_Dist_Power + 1.0)))
       PSI = 2.0 *  Pi_Value * Random_Value_2
       Launch_Angle = ATAN (TAN (Beta) * COS (Psi))
    End If

    !write (6,'(a,i6,2g12.5,i6,5g12.5)') "Launch angle is:",Velocity_Angle_Option,Beta,Psi,IProd,Launch_Angle,&
    !     & random_value_1,random_value_2,sqrt(random_value_1),asin(sqrt(random_value_1))
    !write (0,'(a,i6,2g12.5,i6,5g12.5)') "Launch angle is:",Velocity_Angle_Option,Beta,Psi,IProd,Launch_Angle,&
    !     & random_value_1,random_value_2,sqrt(random_value_1),asin(sqrt(random_value_1))

  End Subroutine HC_Launch_Angle
 
  Subroutine Target_Normal (Launch_Segment,Dir_Indicator,Launch_Angle,True_Tangent,Tangent_Launch,Current_Angle,Current_Tangent,&
&Status)
 
    ! Find normal angle to specified wall segment and add to input angle.
    ! For side puff gas injection, puff along Y=0, ie. reset tanlan to 0.
    ! This may result in some neutrals being launched straight at the
    ! target surface - count these but do not follow (reflect) them.
    ! Variables ANGLAN and TANLAN refer to +Y region only.  True values
    ! used are ANGLE and TANGNT which correct for -Y or +Y regions.
 
    Use ComHC
    Use HC_Init_DIV_Data
    Use HC_Get
 
    ! Every good Fortran 90 program has...
    Implicit None
 
    ! Declare input/output variables.
    Integer, Intent (In) :: Launch_Segment ! ID
    Integer, Intent (In) :: Dir_Indicator ! IS
    Real, Intent (In) :: Launch_Angle ! ANGLAN
    Real, Intent (Out) :: True_Tangent ! TANTRU
    Real, Intent (Out) :: Tangent_Launch ! TANLAN
    Real, Intent (Out) :: Current_Angle ! ANGLE
    Real, Intent (Out) :: Current_Tangent ! TANGNT
    Integer, Intent (In) :: Status
 
    ! Declare local variables.
    Integer :: Launch_Option
 
    ! Assign local variables from modules.
    Launch_Option = hc_launch_location ! Set to CNEUTC or user specified for independent HC use.
 
    If (Launch_Option .eq. 1 .or. Launch_Option .eq. 5) Then
       ! Launch from given R,Z or from cell centres.  No need to adjust launch angle.
       True_Tangent = 0.0
    ElseIf (Launch_Option .eq. 2 .or. Launch_Option .eq. 4 .or. Launch_Option .eq. 6) Then
       ! Launch from a wall segment.  "Launch_Segment" is a WALL.
       ! The IS variable contains the information about which
       ! side of a wall launch point the particle is launched from.
       If (Dir_Indicator .eq. 0) Then
          ! Anti-Clockwise
          True_Tangent = gwallpt (Launch_Segment,8) + 0.5 *  Pi_Value
       ElseIf (Dir_Indicator .eq. 1) Then
          ! Clockwise
          True_Tangent = gwallpt (Launch_Segment,9) - 0.5 *  Pi_Value
       End If
 
    ElseIf (Launch_Option .eq. 0 .or. Launch_Option .eq. 3) Then
       ! Launch from target.  "Launch_Segment" is a TARGET.
       ! Adjust launch angle by the target's segment angle.
       True_Tangent = gthetas (Launch_Segment)
    End If
 
    !  Normal_Measure_Opt = CNEUTE in init_data.
    If ( Normal_Measure_Opt .eq. 0 .or. ( Normal_Measure_Opt .eq. 1 .and. Status .gt. 1)) Then
       Tangent_Launch = True_Tangent
    Else
       Tangent_Launch =  Normal_Measure_To_X_Axis
    End If
 
    Current_Angle = Launch_Angle
    Current_Tangent = Tangent_Launch	
    !write (0,*) "TARGET_NORMAL",Launch_Segment,Dir_Indicator,Launch_Angle,True_Tangent,Tangent_Launch,Current_Angle,Current_Tangent,Status
  End Subroutine Target_Normal
 
  Subroutine HC_Launch_Velocity (IProd,LProd,Current_Cell,Current_Ring,Segment_Index,HC_Species,H_Isotope_Composition, &
       & Sput_Weight,Max_Velocity_Randoms,Beta,PSI,Seed,NRand,Launch_Reg,HC_Temperature,Velocity_Multiplier,Launch_Velocity,IFate)
 
    ! Calculate launch velocity of hydrocarbon species.  First select random number in the
    ! range 0 < RAN < RANMAX.  If not found, reject number and find another one.
    ! To prevent a potentially infinite loop, count number of rejected velocities in
    ! NREJEC and limit to (say) 1000.  If this limit is exceeded, then the launch is
    ! called a "failed launch".
 
    ! Use required modules.
    Use HC_Init_DIV_Data
    Use HC_Init_DIV_Diag
    Use HC_Get
    Use ComHC ! Number_H_Species.
    Use HC_Init_Lib_Data
 
    ! Every good Fortran program has...
    Implicit None
 
    ! Declare input/output variables.
    Integer, Intent (In) :: IProd
    Integer, Intent (In) :: LProd
    Integer, Intent (In) :: Current_Cell ! IK
    Integer, Intent (In) :: Current_Ring ! IR
    Integer, Intent (In) :: Segment_Index ! ID (wall or target).
    Integer, Intent(In) :: HC_Species
    Integer, Dimension (Number_H_Species), Intent (In) :: H_Isotope_Composition ! Number of H/D/T's in the HC.
    !Integer, Dimension (Number_H_Products), Intent (In) :: H_Isotope_Composition ! Number of H/D/T's in the HC.
    Real, Intent (In) :: Sput_Weight
    Real, Intent (In) :: Max_Velocity_Randoms ! RMAXS
    Real, Intent (In) :: Beta ! BETA, required to find VMULT.
    Real, Intent (In) :: Psi ! PSI, required to find VMULT.
    Double Precision, Intent (In) :: Seed ! Use for new random numbers.
    Integer, Intent (InOut) :: NRand ! NRAND
    Integer, Intent (In) :: Launch_Reg ! Used for statistics.
    Real, Intent (Out) :: HC_Temperature ! TEMN
    Real, Intent (Out) :: Velocity_Multiplier ! VMULT
    Real, Intent (Out) :: Launch_Velocity ! VIN
    Integer, Intent (InOut) :: IFate !
 
    ! Declare local variables.
    Integer :: NRejec
    Integer :: Velocity_Angle_Option
    Integer :: Maxwellian_Option
    Character (len=12) :: Input_Type_Dist_1
    Real :: Energy_Temp_Dist_1
    Character (len=12) :: Input_Type_Dist_2
    Real :: Energy_Temp_Dist_2
    Real :: Primary_Contribution
    Real :: HC_Mass
    Real :: Random_Value
 
    ! Assign local variables.
    Velocity_Angle_Option = HC_Launch_Angle_Velocity ! Set to CNEUTC or user specified for independent HC use.
    Maxwellian_Option = HC_Launch_Velocity_Model ! Launch velocity using a single/dual Maxwellian distribution.
 
    ! Find mass of hydrocarbon.
    HC_Mass = Find_HC_Mass (HC_Species,H_Isotope_Composition)

    !
    ! jdemod - initialize the velocity multiplier so that is it set even for a Maxwellian distribution - otherwise
    !          the code elsewhere won't work for that option - the velocity multiplier used is saved at the end of this
    !          routine in saved_velocity_mult
    !
    velocity_multiplier=1.0

 
    ! Check for Maxwellian velocity distribution.
    If (Maxwellian_Option .eq. 0) Then
       ! Velocity calculation done as a constant.
       ! Calculate velocity multiplier VMULT.
       If (Velocity_Angle_Option .eq. 1 .or. Velocity_Angle_Option .eq. 4 .or. Velocity_Angle_Option .eq. 12 .or. &
          &Velocity_Angle_Option .eq. 13 .or. Velocity_Angle_Option .eq. 14) Then
          Velocity_Multiplier = SQRT (ABS(COS(Beta)**2+(SIN(Beta)**2*COS(Psi)**2)))
       Else
          Velocity_Multiplier = 1.0
       End If
 
       If (IProd .le. Max_Impurities) Then
 
          If (geprods (IProd) .gt. 0.0) Then
             !HC_Temperature  = geprods(IProd)
             If (Get_HC_Charge (HC_Species) .eq. 0) Then
                HC_Temperature = HC_Sput_Energy_Neutral_Preset
             Else
                HC_Temperature = HC_Sput_Energy_Ion_Preset
             End If
             Launch_Velocity = 1.38E4 * SQRT (HC_Temperature / HC_Mass)
          Else
             !Random_Value = granvc(IProd)
             NRand = NRand + 1
             Call Surand2 (Seed, 1, Random_Value)
 
             NRejec = 0
             !write (0,*) "velocity const",Max_Velocity_Randoms
             Do
                ! Check for an applicable lower energy cutoff.  None should currently apply for chemical sputtering.
                If (Random_Value .gt. Max_Velocity_Randoms) Then
                   HC_Num_Vel_Greater_Max_Ran(HC_Species,Launch_Reg)=HC_Num_Vel_Greater_Max_Ran(HC_Species,Launch_Reg)+Sput_Weight
                   NRejec = NRejec + 1
                   Call Surand2 (Seed, 1, Random_Value)
                   NRand = NRand + 1
                   If (NRejec .lt. 1000) Then
                      Cycle
                   End If
 
                   ! Write WBC data for neutral redep.
                   If (hc_wbc_comp_option .ne. 0) Then
                      HC_Tot_Temp_At_WBC_HC_Target(HC_Species,Launch_Reg)=HC_Tot_Temp_At_WBC_HC_Target(HC_Species,Launch_Reg)+&
                           &Sput_Weight*(0.5*HC_Mass*1.67E-27*Launch_Velocity*Launch_Velocity/1.602E-19+2*HC_Temperature)
                   End If
 
                   Launch_Velocity = 0.0
                   HC_Temperature  = 0.0
                   HC_Num_Failed_Launches (HC_Species,Launch_Reg) =  HC_Num_Failed_Launches (HC_Species,Launch_Reg) + Sput_Weight
                   IFate = 14 ! Hydrocarbon failed launch.  Record as an MTC event.
                   Return
                End If
                ! Exit do loop.
                Exit
             End Do
 
             ! HC_Temperature =  Init_Particle_Temperature
             If (Get_HC_Charge (HC_Species) .eq. 0) Then
                HC_Temperature = HC_Sput_Energy_Neutral_Preset
             Else
                HC_Temperature = HC_Sput_Energy_Ion_Preset
             End If
 
             If (Velocity_Angle_Option .eq. 0 .or. Velocity_Angle_Option .eq. 1 .or. Velocity_Angle_Option .eq. 4 .or. &
                &Velocity_Angle_Option .eq. 5) Then
                Launch_Velocity=1.38E4*SQRT(Target_Binding_Energy/(1.0/SQRT(Random_Value)-1.0)/HC_Mass)*Velocity_Multiplier
             ElseIf (Velocity_Angle_Option .eq. 2) Then
                If (Random_Value .ge. 1.0) Then
                   Random_Value = 0.999999
                End If
                Launch_Velocity = 1.38E4 * SQRT (HC_Temperature * ABS (LOG (1.0-Random_Value)) / HC_Mass)
             ElseIf (Velocity_Angle_Option .eq. 3 .or. Velocity_Angle_Option .eq. 6 .or. Velocity_Angle_Option .eq. 7 .or. &
                    &Velocity_Angle_Option .eq. 8 .or. Velocity_Angle_Option .eq. 10 .or. Velocity_Angle_Option .eq. 11) Then
                Launch_Velocity = 1.38E4 * SQRT (HC_Temperature / HC_Mass)
 
                ! Note 156 VEL/ANG 9 Flag. Notice that launches alternatively at
                ! EIN1, EIN2,  but we also have the +/-Y thing.  Hence here
                ! we need to launch 2 particles at EIN1, followed by 2 at EIN2,
                ! etc to ensure that some of each are launched on each side of Y=0.
             ElseIf (Velocity_Angle_Option .eq. 9) Then
                If (2 * (IPRod / 4) .eq. IPRod / 2) Then
                   Launch_Velocity = 1.38E4 * SQRT (HC_Temperature / HC_Mass)
                Else
                   Launch_Velocity = 1.38E4 * SQRT ( Alt_Init_Particle_Temperature / HC_Mass)
                End If
             ElseIf (Velocity_Angle_Option .eq. 12) Then
                Launch_Velocity = 1.38E4 * SQRT (HC_Temperature / HC_Mass) * Velocity_Multiplier
             ElseIf (Velocity_Angle_Option .eq. 13) Then
                If (Random_Value .ge. 1.0) Then
                   Random_Value = 0.999999
                End If
                Launch_Velocity = 1.38E4 * SQRT (HC_Temperature * ABS (LOG (1.0 - Random_Value)) / HC_Mass) * Velocity_Multiplier
             ElseIf (Velocity_Angle_Option .eq. 14) then
                Launch_Velocity = 1.38E4 * SQRT (gktids (Segment_Index) / HC_Mass) * Velocity_Multiplier * Const_velocity_mult
             ElseIf (Velocity_Angle_Option .eq. 15) then
                Launch_Velocity=1.38E4*SQRT(gktibs(Current_Cell,Current_Ring)/HC_Mass)*Velocity_Multiplier*Const_velocity_mult
             End If
          End If
       Else
          ! Default to constant velocity if IProd greater than allowed impurities.
          !HC_Temperature  = geprods(IProd)
          If (Get_HC_Charge (HC_Species) .eq. 0) Then
             HC_Temperature = HC_Sput_Energy_Neutral_Preset
          Else
             HC_Temperature = HC_Sput_Energy_Ion_Preset
          End If
          Launch_Velocity = 1.38E4 * SQRT (HC_Temperature / HC_Mass)
 
       End If
    Else
       ! Maxwellian velocity option is selected.
       ! Call one random value to determine velocity of particle.
       ! Note: SIGN (A,B) returns the Abs(A) * sign(B).
       Call SURAND2 (Seed, 1, Random_Value)
 
       ! Note, this is no longer used.  Particle temperature is determined by the surface temperature.
       If (Get_HC_Charge (HC_Species) .eq. 0) Then
          HC_Temperature = HC_Sput_Energy_Neutral_Preset
       Else
          HC_Temperature = HC_Sput_Energy_Ion_Preset
       End If
       ! Convert to temperature.
       Energy_Temp_Dist_1 = HC_Temperature *  ECH /  KBoltz
       Input_Type_Dist_1 = "Temperature"
 
       ! Acquire segment temperature to decide thermal launch distribution.
       If (gwallpt (Segment_Index,16) .eq. 1 .or. gwallpt (Segment_Index,16) .eq. 4) Then
          ! Particle is launched from an inner or outer target segment (1 and 4 respectively).
          ! Note: Find_Target_Temperature assumes input of a target segment value, not wall segment.
          Energy_Temp_Dist_1 = Find_Target_Temperature (INT(gwallpt(Segment_Index,18)))
          Input_Type_Dist_1 = "Temperature"
       Else
          ! If particle is launched from a wall (7), baffle (9) or private plasma zone segment (8).
          Energy_Temp_Dist_1 = Find_Wall_Temperature (Segment_Index)
          Input_Type_Dist_1 = "Temperature"
       End If
 
       If (hc_dual_mb_sec_mean_temp .gt. 100) Then
          ! Input is most likely in degrees K.
          Input_Type_Dist_2 = "Temperature"
       Else
          ! Input is most likely in eV.
          Input_Type_Dist_2 = "Energy"
       End If
       Energy_Temp_Dist_2 = hc_dual_mb_sec_mean_temp
       Primary_Contribution = hc_dual_mb_pri_vel_flux
 
       If (Maxwellian_Option .eq. 1) Then
          ! Override secondary distribution with 0.0 contribution.
          Primary_Contribution = 1.0			
       End If
 
       Call Maxwellian_Launch_Velocity_v3 (Input_Type_Dist_1,Energy_Temp_Dist_1,Input_Type_Dist_2,Energy_Temp_Dist_2,&
&Primary_Contribution,HC_Mass,Random_Value,Launch_Velocity)
       !write (0,*) "surtemp:",Energy_Temp_Dist_1,Energy_Temp_Dist_2,Segment_Index,Find_Target_Temperature (Segment_Index),Find_Wall_Temperature (Segment_Index),gwallpt (Segment_Index,16),INT(gwallpt(Segment_Index,18)),Launch_Velocity
 
    End If
    !Write (0,'(a,5(1x,g12.5))') "Launch Velocity is:",Launch_Velocity,HC_mass,HC_temperature,Random_Value



    ! jdemod - record the velocity_multiplier used for this particle launch
    saved_velocity_mult = velocity_multiplier

 
  End Subroutine HC_Launch_Velocity
 
  Subroutine HC_Injection_Velocity (Current_Cell,Current_Ring,HC_Species, &
       & H_Isotope_Composition,Sput_Weight,Seed,NRand,HC_Porm,HC_Temperature,Launch_Velocity)
 
    ! Calculate launch velocity of hydrocarbon species for an ion injection case.
    ! jdemod - NOTE!!: The HC code expects the ion velocity in terms of m/s - it is only converted to
    !                  dist/timestep in the ion_move routine in hc_ion_transport. As a result, it is bug
    !                  here to scale the initial velocity by the ion_time_step. 
 
    ! Use required modules.
    Use HC_Init_DIV_Data
    Use HC_Get
    Use ComHC
    Use HC_Init_Lib_Data
 
    ! Every good Fortran program has...
    Implicit None
 
    ! Declare input/output variables.
    Integer, Intent (In) :: Current_Cell ! IK
    Integer, Intent (In) :: Current_Ring ! IR
    Integer, Intent(In) :: HC_Species
    Integer, Dimension (Number_H_Species), Intent (In) :: H_Isotope_Composition ! Number of H/D/T's in the HC.
    !Integer, Dimension (Number_H_Products), Intent (In) :: H_Isotope_Composition ! Number of H/D/T's in the HC.
    Real, Intent (In) :: Sput_Weight
    Double Precision, Intent (In) :: Seed ! Use for new random numbers.
    Integer, Intent (InOut) :: NRand ! NRAND		
    Real, Intent (InOut) :: HC_Porm ! PORM
    Real, Intent (Out) :: HC_Temperature ! TEMN
    Real, Intent (Out) :: Launch_Velocity ! VIN
 
    ! Declare local variables.
    Real :: Ran1
    Real :: Ran2
 
    HC_Porm = -1.0 * HC_Porm ! Direction modifier.
 
    If ( Injection_Opt .eq. 1) Then

       !jdemod - remove time step
       !Launch_Velocity=9.79e3*SQRT(Init_Particle_Temperature/Find_HC_Mass(HC_Species,H_Isotope_Composition))*HC_Porm*Ion_Time_Step

       Launch_Velocity=9.79e3*SQRT(Init_Particle_Temperature/Find_HC_Mass(HC_Species,H_Isotope_Composition))*HC_Porm
 
    ElseIf ( Injection_Opt .eq. 2 .or.  Injection_Opt .eq. 3) Then
       !jdemod - remove time step
       !Launch_Velocity =9.79e3*SQRT(Init_Particle_Temperature/Find_HC_Mass(HC_Species,H_Isotope_Composition))*HC_Porm*Ion_Time_Step
       Launch_Velocity =9.79e3*SQRT(Init_Particle_Temperature/Find_HC_Mass(HC_Species,H_Isotope_Composition))*HC_Porm
 
    ElseIf ( Injection_Opt .eq. 5 .or.  Injection_Opt .eq. 6) Then
       Do
          NRand = NRand + 1
          Call Surand2 (Seed, 1, Ran1)
          If (Ran1 .ne. 0.0) Then
             Exit
          End If
       End Do
       Do
          NRand = NRand + 1
          Call Surand2 (Seed, 1, Ran2)
          If (Ran2 .ne. 0.0) Then
             Exit
          End If
       End Do
       Local_RGauss = SQRT (-2.0 * LOG (Ran1)) * COS (2.0 *  Pi_Value * Ran2)
       !jdemod - remove time step
       !Launch_Velocity = 9.79e3 * SQRT ( Init_Particle_Temperature / Find_HC_Mass (HC_Species,H_Isotope_Composition)) *  &
       !                & Ion_Time_Step *  Local_RGauss
       Launch_Velocity = 9.79e3 * SQRT ( Init_Particle_Temperature / Find_HC_Mass (HC_Species,H_Isotope_Composition)) *  &
                       & Local_RGauss
 
    ElseIf ( Injection_Opt .eq. 4 .or.  Injection_Opt .eq. 7) Then
       !jdemod - remove time step
       !Launch_Velocity = 9.79E3 * SQRT (2.0 * gpinenz (Current_Cell,Current_Ring) / Find_HC_Mass (HC_Species,H_Isotope_Composition)&
       !                & ) * 0.5 * HC_Porm *  Ion_Time_Step
       Launch_Velocity = 9.79E3 * SQRT (2.0 * gpinenz (Current_Cell,Current_Ring) / Find_HC_Mass (HC_Species,H_Isotope_Composition)&
                       & ) * 0.5 * HC_Porm
 
    End If
 
    ! Set initial ion temperature.
    If ( Injection_Opt .eq. 7 .or.  Injection_Opt .eq. 4) Then
       HC_Temperature = gpinenz (Current_Cell,Current_Ring)
    Else
       HC_Temperature =  Init_Particle_Temperature
    End If
 
  End Subroutine HC_Injection_Velocity
 
  Subroutine Maxwellian_Launch_Velocity (Input_Type_Dist_1,Energy_Temp_Dist_1,Input_Type_Dist_2, &
       & Energy_Temp_Dist_2,Primary_Contribution,Mass,Random_Value,Launch_Velocity)
 
    ! Used to determine launch velocity of released atoms and molecules
    ! under single or double Maxwellian distribution functions.  Implements
    ! the dual Maxwell-Boltzmann distribution model for hydrocarbon release
    ! of E. Vietzke described in "Energy distributions of CD4 and CD3
    ! chemically released from graphite by D+ and D0/Ne+ impact", Journal
    ! of Nuclear Materials 290-293 (2001) 158-161.  MB distributions are
    ! integrated with a variable-subinterval Simpson's solver over the
    ! velocity space from 0 m/s to z * max (avg. velocity associated with
    ! the two input MB distributions) where z is typically 5.  If a single
    ! MB distribution is desired, a primary contribution value of 1.0 can
    ! be used to negate the secondary distribution.
 
    ! Use required modules.
    Use ComHC
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Declare input/output variables.
    Character (len=*), Intent (In) :: Input_Type_Dist_1
    Real, Intent (In) :: Energy_Temp_Dist_1
    Character (len=*), Intent (In) :: Input_Type_Dist_2
    Real, Intent (In) :: Energy_Temp_Dist_2
    Real, Intent (In) :: Primary_Contribution
    Real, Intent (In) :: Mass
    Real, Intent (In) :: Random_Value
    Real, Intent (Out) :: Launch_Velocity
 
    ! Declare local constants.
    Real :: Mass_1_AMU = 1.67E-27
    Real :: Electronic_Charge = 1.60219E-19
    Real :: Boltzmann_Constant = 1.38066E-23
    Integer :: Subintervals = 10000
    Integer :: Upper_Velocity_Multiplier = 6
    Real :: Lower_Velocity_Limit = 0.0
 
    ! Declare local variables.
    Real :: Beta_Dist_1 ! Beta = m/2kT1.
    Real :: Beta_Dist_2 ! Beta = m/2kT2.
    Real :: Avg_Velocity_Dist_1 ! Average velocity asociated with energy of MB distribution 1.
    Real :: Avg_Velocity_Dist_2 ! Average velocity asociated with energy of MB distribution 2.
    Real :: Upper_Velocity_Limit ! Upper limit to check against.
    Real :: Delta_V ! Velocity step for each interval.
    Real :: Interval_Sum ! Ongoing sum of integrated area.
    Real :: Current_Velocity ! Velocity for integration routine.
    Real :: Current_Velocity_Squared ! V^2 calculated once to save computation.
    Real :: Simpson_Weight ! Weighting factor.
    Real :: Dist_1 ! Evaluation of primary MB function.
    Real :: Dist_2 ! Evaluation of secondary MB function.
    Integer :: Integration_Steps ! Tracks number of steps taken in Simpson's integration.
 
    ! Read in input energy (typically used for energetic launch) or temperature
    ! (typically used for thermal launch) and calculate Beta constant.
    If (Input_Type_Dist_1 .eq. "Temperature") Then
       Beta_Dist_1 = (Mass * Mass_1_AMU) / (2.0 * Energy_Temp_Dist_1 * Boltzmann_Constant)
       Avg_Velocity_Dist_1 = SQRT (8.0 * Boltzmann_Constant * Energy_Temp_Dist_1 / (Pi * Mass * Mass_1_AMU))
    ElseIf (Input_Type_Dist_1 .eq. "Energy") Then
       Beta_Dist_1 = (Mass * Mass_1_AMU) / (2.0 * Energy_Temp_Dist_1 * Electronic_Charge)
       Avg_Velocity_Dist_1 = SQRT (8.0 * Electronic_Charge * Energy_Temp_Dist_1 / (Pi * Mass * Mass_1_AMU))
    Else
       Write (Output_Unit_HC_Alert,*) "Error in Maxwellian_Launch_Velocity: Incorrect input energy or temperature to &
&Maxwelliam_Launch_Velocity primary distribution."
       Write (Output_Unit_HC_Alert,*) "Please use 'energy' or 'temperature' as input parameters.  Program stopping."
       Stop
    End If
 
    If (Input_Type_Dist_2 .eq. "Temperature") Then
       Beta_Dist_2 = (Mass * Mass_1_AMU) / (2.0 * Energy_Temp_Dist_2 * Boltzmann_Constant)
       Avg_Velocity_Dist_2 = SQRT (8.0 * Boltzmann_Constant * Energy_Temp_Dist_2 / (Pi * Mass * Mass_1_AMU))
    ElseIf (Input_Type_Dist_2 .eq. "Energy") Then
       Beta_Dist_2 = (Mass * Mass_1_AMU) / (2.0 * Energy_Temp_Dist_2 * Electronic_Charge)
       Avg_Velocity_Dist_2 = SQRT (8.0 * Electronic_Charge * Energy_Temp_Dist_2 / (Pi * Mass * Mass_1_AMU))
    Else
       Write (Output_Unit_HC_Alert,*) "Error in Maxwellian_Launch_Velocity: Incorrect input energy or temperature to &
&Maxwelliam_Launch_Velocity secondary distribution."
       Write (Output_Unit_HC_Alert,*) "Please use 'energy' or 'temperature' as input parameters.  Program stopping."
       Stop
    End If
 
    ! Check that contribution is between 0.0 and 1.0.
    If (Primary_Contribution .lt. 0.0 .or. Primary_Contribution .gt. 1.0) Then
       Write (Output_Unit_HC_Alert,*) "Error in Maxwellian_Launch_Velocity:  Contribution to primary distribution must lay between &
&0.0 and 1.0."
       Write (Output_Unit_HC_Alert,*) "Program stopping."
       Stop
    End If
 
    ! Check that mass is greater than 0.0.
    If (Mass .lt. 0.0) Then
       Write (Output_Unit_HC_Alert,*) "Error in Maxwellian_Launch_Velocity:  Particle mass must be greater than 0.0."
       Write (Output_Unit_HC_Alert,*) "Program stopping."
       Stop
    End If
 
    ! Check that number of subintervals is reasonable.
    If (Subintervals .lt. 10) Then
       Write (Output_Unit_HC_Alert,*) "Warning in Maxwellian_Launch_Velocity: Number of subintervals is <10.  A poor approximation &
&may result."
    ElseIf (Subintervals .gt. 1000000) Then
       Write (Output_Unit_HC_Alert,*) "Warning in Maxwellian_Launch_Velocity: Number of subintervals is >10,000.  A high &
&calculation time may result."
    End If
 
    ! Check that number of subintervals is an even value.
    If (MOD (Subintervals, 2) .gt. 0.0) Then
       Write (Output_Unit_HC_Alert,*) "Error in Maxwellian_Launch_Velocity: Number of subintervals must be an even number."
       Write (Output_Unit_HC_Alert,*) "Program stopping."
       Stop
    End If
 
    ! Check divergence of two distributions.
    If (Primary_Contribution .ne. 0.0 .and. Primary_Contribution .ne. 1.0) Then
       If (Energy_Temp_Dist_1 / Energy_Temp_Dist_2 .gt. (Subintervals / 2.0) .or. Energy_Temp_Dist_1 / Energy_Temp_Dist_2 .lt. (&
&2.0 / Subintervals)) Then
          ! Two distributions will lead to unbalanced result.  Abort.
          Write (Output_Unit_HC_Alert,*) "Error in Maxwellian_Launch_Velocity:  Two distributions are too disjointed.  Resulting &
&integration"
          Write (Output_Unit_HC_Alert,*) "will have too few points in prime MB distribution space."
          Write (Output_Unit_HC_Alert,*) "Program stopping."
          Stop
       End If
    End If
 
    ! Set upper integration limits.
    Upper_Velocity_Limit = Upper_Velocity_Multiplier * MAX (Avg_Velocity_Dist_1, Avg_Velocity_Dist_2)
 
    ! Find velocity intervals.
    Delta_V = (Upper_Velocity_Limit - Lower_Velocity_Limit) / Subintervals
 
    ! Reset sum.
    Integration_Steps = 0
    Launch_Velocity = 0.0
 
    ! Start integration.  Sum up intervals until total reaches the input Random_Value.
    Do Current_Velocity = Lower_Velocity_Limit, Upper_Velocity_Limit, Delta_V
 
       Integration_Steps = Integration_Steps + 1
       Current_Velocity_Squared = Current_Velocity**2
 
       Interval_Sum=Primary_Contribution*(1-(1+Beta_Dist_1*Current_Velocity_Squared)*EXP(-Beta_Dist_1*Current_Velocity_Squared))+&
            &(1-Primary_Contribution)*(1-(1+Beta_Dist_2*Current_Velocity_Squared)*EXP(-Beta_Dist_2*Current_Velocity_Squared))
 
       !write (0,*) "int",Random_Value,Primary_Contribution,Current_Velocity,Interval_Sum
 
       ! Check against input value.
       If (Interval_Sum .gt. Random_Value) Then
          ! Passed the velocity interval.  Finish up.
          Launch_Velocity = (Current_Velocity + (Current_Velocity - Delta_V)) / 2.0
          Exit
       End If
    End Do
 
    ! Check that a value was assigned.  If not, Random_Value is higher than the sum of all intervals.
    ! This implies that the velocity to return must be higher than the maximum anticipated by Upper_Velocity_Multiplier.
    If (Launch_Velocity .eq. 0.0) Then
       ! Check if Interval_Sum is within 1% of Random_Value.
       If (ABS (Interval_Sum - Random_Value) .lt. 0.01) Then
          Write(Output_Unit_HC_Alert,*)"Warning in Maxwellian_Launch_Velocity:Launch velocity is higher than:",Upper_Velocity_Limit
          Launch_Velocity = (Upper_Velocity_Multiplier + 1) * MAX (Avg_Velocity_Dist_1, Avg_Velocity_Dist_2)
          Write (Output_Unit_HC_Alert,*) "Using:", Launch_Velocity
       Else
          ! Must raise maximum velocity for integration - not enough velocity space is being considered.
          Write (Output_Unit_HC_Alert,*) "Error in Maxwellian_Launch_Velocity: Upper velocity multiplier is too low:", &
&Upper_Velocity_Multiplier, Upper_Velocity_Limit
          Write (Output_Unit_HC_Alert,*) "or number of intervals is too low:", Subintervals, Interval_Sum
          Write (Output_Unit_HC_Alert,*) "Program stopping."
          Stop
       End If
    End If
 
  End Subroutine Maxwellian_Launch_Velocity
 
  Subroutine Maxwellian_Launch_Velocity_v3 (Input_Type_Dist_1,Energy_Temp_Dist_1,Input_Type_Dist_2, &
       & Energy_Temp_Dist_2,Primary_Contribution,Mass,Random_Value,Launch_Velocity)
 
    ! Used to determine launch velocity of released atoms and molecules
    ! under single or double Maxwellian distribution functions.  Impliments
    ! the dual Maxwell-Boltzman distribution model for hydorcarbon relase
    ! of E. Vietzke described in "Energy distributions of CD4 and CD3
    ! chemically released from graphite by D+ and D0/Ne+ impact", Journal
    ! of Nuclear Materials 290-293 (2001) 158-161.  MB distributions are
    ! integrated with a variable-subinterval Simpson's solver over the
    ! velocity space from 0 m/s to z * max (avg. velocity associated with
    ! the two input MB distributions) where z is typically 5.  If a single
    ! MB distribution is desired, a primary contribution value of 1.0 can
    ! be used to negate the secondary distribution.
 
    ! External modules
    Use ComHC ! Output_Unit_HC_Alert.
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Declare input/output variables.
    Character (len=*), Intent (In) :: Input_Type_Dist_1
    Real, Intent (In) :: Energy_Temp_Dist_1
    Character (len=*), Intent (In) :: Input_Type_Dist_2
    Real, Intent (In) :: Energy_Temp_Dist_2
    Real, Intent (In) :: Primary_Contribution
    Real, Intent (In) :: Mass
    Real, Intent (In) :: Random_Value
    Real, Intent (Out) :: Launch_Velocity
 
    ! Declare local constants.
    Real :: Mass_1_AMU = 1.67E-27
    Real :: Electronic_Charge = 1.60219E-19
    Real :: Boltzmann_Constant = 1.38066E-23
    Integer :: Subintervals = 10000
    Integer :: Upper_Velocity_Multiplier = 6
    Real :: Lower_Velocity_Limit = 0.0
 
    ! Declare local variables.
    Real :: Beta_Dist_1 ! Beta = m/2kT1.
    Real :: Beta_Dist_1_Three_Halves ! Beta1^1.5.
    Real :: Beta_Dist_2 ! Beta = m/2kT2.
    Real :: Beta_Dist_2_Three_Halves ! Beta2^1.5.
    Real :: Avg_Velocity_Dist_1 ! Average velocity associated with energy of MB distribution 1.
    Real :: Avg_Velocity_Dist_2 ! Average velocity associated with energy of MB distribution 2.
    Real :: Upper_Velocity_Limit ! Upper limit to check against.
    Real :: Delta_V ! Velocity step for each interval.
    Real :: Interval_Sum ! Ongoing sum of integrated area.
    Real :: Current_Velocity ! Velocity for integration routine.
    Real :: Current_Velocity_Squared ! V^2 calculated once to save computation.
    Real :: Simpson_Weight ! Weighting factor.
    Real :: Dist_1 ! Evaluation of primary MB function.
    Real :: Dist_2 ! Evaluation of secondary MB function.
    Real :: SQRT_Pi ! Square root of PI.
    Integer :: Integration_Steps ! Tracks number of steps taken in Simpson's integration.
 
    SQRT_Pi = SQRT (Pi)
 
    ! Check that input energy/temperatures are >0.0.
    If (Energy_Temp_Dist_1 .le. 0.0) Then
       Write (Output_Unit_HC_Alert,*) "Error in Maxwellian_Launch_Velocity_v3:  Input energy/temperature 1 must be >0.0."
       Write (Output_Unit_HC_Alert,*) "Program stopping."
       Stop
    ElseIf (Energy_Temp_Dist_2 .le. 0.0) Then
       Write (Output_Unit_HC_Alert,*) "Error in Maxwellian_Launch_Velocity_v3:  Input energy/temperature 2 must be >0.0."
       Write (Output_Unit_HC_Alert,*) "Program stopping."
       Stop
    End If
 
    ! Read in input energy (typically used for energetic launch) or temperature
    ! (typically used for thermal launch) and calculate Beta constant.
    If (Input_Type_Dist_1 .eq. "Temperature") Then
       Beta_Dist_1 = (Mass * Mass_1_AMU) / (2.0 * Energy_Temp_Dist_1 * Boltzmann_Constant)
       Avg_Velocity_Dist_1 = SQRT (8.0 * Boltzmann_Constant * Energy_Temp_Dist_1 / (Pi * Mass * Mass_1_AMU))
    ElseIf (Input_Type_Dist_1 .eq. "Energy") Then
       Beta_Dist_1 = (Mass * Mass_1_AMU) / (2.0 * Energy_Temp_Dist_1 * Electronic_Charge)
       Avg_Velocity_Dist_1 = SQRT (8.0 * Electronic_Charge * Energy_Temp_Dist_1 / (Pi * Mass * Mass_1_AMU))
    Else
       Write (Output_Unit_HC_Alert,*) "Error in Maxwellian_Launch_Velocity_v3: Incorrect input energy or temperature to &
&Maxwelliam_Launch_Velocity primary distribution."
       Write (Output_Unit_HC_Alert,*) "Please use 'energy' or 'temperature' as input parameters.  Program stopping."
       Stop
    End If
    Beta_Dist_1_Three_Halves = Beta_Dist_1**(1.5)
 
    If (Input_Type_Dist_2 .eq. "Temperature") Then
       Beta_Dist_2 = (Mass * Mass_1_AMU) / (2.0 * Energy_Temp_Dist_2 * Boltzmann_Constant)
       Avg_Velocity_Dist_2 = SQRT (8.0 * Boltzmann_Constant * Energy_Temp_Dist_2 / (Pi * Mass * Mass_1_AMU))
    ElseIf (Input_Type_Dist_2 .eq. "Energy") Then
       Beta_Dist_2 = (Mass * Mass_1_AMU) / (2.0 * Energy_Temp_Dist_2 * Electronic_Charge)
       Avg_Velocity_Dist_2 = SQRT (8.0 * Electronic_Charge * Energy_Temp_Dist_2 / (Pi * Mass * Mass_1_AMU))
    Else
       Write (Output_Unit_HC_Alert,*) "Error in Maxwellian_Launch_Velocity_v3: Incorrect input energy or temperature to &
&Maxwelliam_Launch_Velocity secondary distribution."
       Write (Output_Unit_HC_Alert,*) "Please use 'energy' or 'temperature' as input parameters.  Program stopping."
       Stop
    End If
    Beta_Dist_2_Three_Halves = Beta_Dist_2**(1.5)
    !write (0,*) "beta1",Beta_Dist_1,"avg1",Avg_Velocity_Dist_1
    !write (0,*) "beta2",Beta_Dist_2,"avg2",Avg_Velocity_Dist_2
 
    ! Check that contribution is between 0.0 and 1.0.
    If (Primary_Contribution .lt. 0.0 .or. Primary_Contribution .gt. 1.0) Then
       Write (Output_Unit_HC_Alert,*) "Error in Maxwellian_Launch_Velocity_v3:  Contribution to primary distribution must lay &
&between 0.0 and 1.0."
       Write (Output_Unit_HC_Alert,*) "Program stopping."
       Stop
    End If
 
    ! Check that mass is greater than 0.0.
    If (Mass .lt. 0.0) Then
       Write (Output_Unit_HC_Alert,*) "Error in Maxwellian_Launch_Velocity_v3:  Particle mass must be greater than 0.0."
       Write (Output_Unit_HC_Alert,*) "Program stopping."
       Stop
    End If
 
    ! Check that number of subintervals is reasonable.
    If (Subintervals .lt. 10) Then
       Write (Output_Unit_HC_Alert,*) "Warning in Maxwellian_Launch_Velocity_v3: Number of subintervals is <10.  A poor &
&approximation may result."
    ElseIf (Subintervals .gt. 100000) Then
       Write (Output_Unit_HC_Alert,*) "Warning in Maxwellian_Launch_Velocity_v3: Number of subintervals is >10,000.  A high &
&calculation time may result."
    End If
 
    ! Check that number of subintervals is an even value.
    If (MOD (Subintervals, 2) .gt. 0.0) Then
       Write (Output_Unit_HC_Alert,*) "Error in Maxwellian_Launch_Velocity_v3: Number of subintervals must be an even number."
       Write (Output_Unit_HC_Alert,*) "Program stopping."
       Stop
    End If
 
    ! Check divergence of two distributions.
    If (Primary_Contribution .ne. 0.0 .and. Primary_Contribution .ne. 1.0) Then
       If (Energy_Temp_Dist_1 / Energy_Temp_Dist_2 .gt. (Subintervals / 2.0) .or. Energy_Temp_Dist_1 / Energy_Temp_Dist_2 .lt. (&
&2.0 / Subintervals)) Then
          ! Two distributions will lead to unbalanced result.  Abort.
          Write (Output_Unit_HC_Alert,*) "Error in Maxwellian_Launch_Velocity_v3:  Two distributions are too disjointed.  &
&Resulting integration"
          Write (Output_Unit_HC_Alert,*) "will have too few points in prime MB distribution space."
          Write (Output_Unit_HC_Alert,*) "Program stopping."
          Stop
       End If
    End If
 
    ! Set upper integration limits.
    Upper_Velocity_Limit = Upper_Velocity_Multiplier * MAX (Avg_Velocity_Dist_1, Avg_Velocity_Dist_2)
    !Upper_Velocity_Limit = 5000.0
    ! Find velocity intervals.
    Delta_V = (Upper_Velocity_Limit - Lower_Velocity_Limit) / Subintervals
 
    ! Reset sum.
    Integration_Steps = 0
    Launch_Velocity = 0.0
    Interval_Sum = 0.0
 
    ! Start integration.  Sum up intervals until total reaches the input Random_Value.
    Do Current_Velocity = Lower_Velocity_Limit, Upper_Velocity_Limit, Delta_V
 
       Integration_Steps = Integration_Steps + 1
       Current_Velocity_Squared = Current_Velocity ** 2.0
 
       Interval_Sum = Primary_Contribution * (1.0 - (1.0 + Beta_Dist_1 * Current_Velocity_Squared) * EXP (-Beta_Dist_1 * &
            &Current_Velocity_Squared)) + &
            &(1.0-Primary_Contribution)*(1.0-(1.0+Beta_Dist_2*Current_Velocity_Squared)*EXP(-Beta_Dist_2*Current_Velocity_Squared))
 
       !write (0,*) "int",Random_Value,Primary_Contribution,Beta_Dist_1,Current_Velocity,Current_Velocity_Squared,Interval_Sum
 
       ! Check against input value.
       If (Interval_Sum .gt. Random_Value) Then
          ! Passed the velocity interval.  Finish up.
          Launch_Velocity = (Current_Velocity + (Current_Velocity - Delta_V)) / 2.0
          Exit
       End If
    End Do
 
    ! Check that a value was assigned.  If not, Random_Value is higher than the sum of all intervals.
    ! This implies that the velocity to return must be higher than the maximum anticipated by Upper_Velocity_Multiplier.
    If (Launch_Velocity .eq. 0.0) Then
       ! Check if Interval_Sum is within 1% of Random_Value.
       If (ABS (Interval_Sum - Random_Value) .lt. 0.01) Then
          Write (Output_Unit_HC_Alert,*) "Warning in Maxwellian_Launch_Velocity_v3:  Launch velocity is higher than:",&
                             &Upper_Velocity_Limit
          Launch_Velocity = (Upper_Velocity_Multiplier + 1) * MAX (Avg_Velocity_Dist_1, Avg_Velocity_Dist_2)
          Write (Output_Unit_HC_Alert,*) "Using:", Launch_Velocity
       Else
          ! Must raise maximum velocity for integration - not enough velocity space is being considered.
          Write (Output_Unit_HC_Alert,*) "Error in Maxwellian_Launch_Velocity_v3: Upper velocity multiplier is too low:", &
                                   &Upper_Velocity_Multiplier, Upper_Velocity_Limit
          Write (Output_Unit_HC_Alert,*) "or number of intervals is too low:", Subintervals, Interval_Sum
          Write (Output_Unit_HC_Alert,*) "Program stopping."
          Stop
       End If
    End If
    !write (0,*) "MAXS ",Input_Type_Dist_1, Energy_Temp_Dist_1, Input_Type_Dist_2, Energy_Temp_Dist_2, &
    !& Primary_Contribution, Mass, Random_Value,Launch_Velocity,Upper_Velocity_Limit, &
    !& Lower_Velocity_Limit,Delta_V,Upper_Velocity_Multiplier,Integration_Steps,Interval_Sum,Random_Value
 
  End Subroutine Maxwellian_Launch_Velocity_v3
 
  Subroutine Modify_Taus (Cur_HC_Spec,H_Isotope_Composition,Sput_Weight,Current_Cell,Current_Ring)
!      jdemod - no random numbers used in this routine - removed other unused arguments
!       & HC_Temperature,Current_Cell,Current_Ring,Current_S,Seed,Random_Numbers_Used,NRand)
 
    ! Required modules.
    Use ComHC ! Number_H_Species.
    Use HC_Init_DIV_Data ! Load cell and global properties data structures.
    Use HC_Init_Lib_Data ! Get HC mass and charge functions.
    Use HC_Get
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Declare input/output variables.
    Integer, Intent (In) :: Cur_HC_Spec
    Integer, Dimension (Number_H_Species), Intent (In) :: H_Isotope_Composition ! Number of H/D/T's in the HC.
    !Integer, Dimension (Number_H_Products), Intent (In) :: H_Isotope_Composition ! Number of H/D/T's in the HC.
    Real, Intent (In) :: Sput_Weight
    !Real, Intent (InOut) :: HC_Temperature
    Integer, Intent (In) :: Current_Cell
    Integer, Intent (In) :: Current_Ring
    !Real, Intent (In) :: Current_S
    !Double Precision, Intent (In) :: Seed
    !Integer, Intent (InOut) :: Random_Numbers_Used ! KK
    !Integer, Intent (InOut) :: NRand
 
    ! Define local variables that will be assigned values.
    Real :: Current_HC_Mass
    Integer :: Current_HC_Charge
    Real :: RIZB ! Background charge.
    Real :: CRMHC ! Local mass of hydrocarbon molecule.
    Real :: CRMI ! Mass of DIVIMP impurity.
    Real :: CRMB ! Background plasma ion mass.
    Integer :: IrSpec ! Special ring where high-resolution measurement begins.
    Real :: CTEMAV ! Cell average temperature.
    Real :: QTIM ! Ion time step.
    Real :: CIOPTB ! Collision option.
    Real :: CIOPTC ! Friction option.
    Real :: CIOPTD ! Heating option.
    Real :: CZENH ! Z enhancement factor.
    Integer :: CDIFOP ! First diffuse option.
    ! Define local variables that will contain calculated values.		
    Real :: Tau_Mass_Factor ! Sqrt (CRMI/CRMHC) multiplication factor.
    Real :: C215A ! Intermediate step for TAU Heating calc.
    Real :: C350A ! Intermediate step for TAU Stopping calc.
    Real :: C350B ! Intermediate step for TAU Parallel calc.
    Real :: RIZSQR ! Current_HC_Charge^2.
    Real :: ROOTTT ! Square root of CTEMAV.
    Real :: ROOTMI ! Square root of CRMHC.
    Real :: TAU ! Intermediate step for TAU Heating calc.
    Real :: STAU ! Intermediate step for TAU Heating calc.
    Real :: FTAU ! Intermediate step for TAU Heating calc.
    Real :: FTAUP ! Intermediate step for TAU Heating calc.
    Real :: FTAUS ! Intermediate step for TAU Heating calc.
    Real :: FTAUT ! Intermediate step for TAU Heating calc.
    Real :: Ran1 ! Temporary random number.
    Real :: Ran2 ! Temporary random number.
    Real :: Ratio1
    Real :: Ratio2
    Real :: Lambda
    Real :: Reduced_Mass ! MU
    Logical :: Debug_Modify_Taus
 
    ! Assign local variables.
    Current_HC_Charge = Get_HC_Charge (Cur_HC_Spec)
    Current_HC_Mass = Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition)		
    Debug_Modify_Taus = .False.
 
    ! Assign local variables from global data tables.
    RIZB =  Back_Plasma_Charge
    CRMHC = Current_HC_Mass
    CRMI =  Impurity_Ion_Mass
    CRMB =  Back_Plasma_Ion_Mass
    IRSpec =  Special_Ring
    CTEMAV = gctemav ()
    QTIM =  Ion_Time_Step
    CIOPTB =  Collision_Opt
    CIOPTC =  Friction_Opt
    CIOPTD =  Heating_Opt
    CZENH =  Z_Enh_Factor
    CDIFOP =  First_Diffuse_Opt
 
    ! Calculation for Lambda - see Dolan, page 40.
    If (hc_lambda_calc .eq. 0) Then
       ! Use constant for lambda.
       Lambda = 15.0
    Else If (hc_lambda_calc .eq. 1) Then
       ! Use calculation for Lambda found in Dolan for ion-collision dominated plasma.
       ! Originally in Sivukhin, D.V., Coulomb collisions in a fully ionized plasma in
       ! Review of Plasma Physics (Consultation Bureau, New York, 1966) Vol. 4, p.88.
       Lambda = 30.0 - 0.5 * LOG (gknbs (Current_Cell, Current_Ring)) + 1.5 * LOG (gktibs (Current_Cell, Current_Ring))
    Else
       ! Unsupported hc_lambda_calc option.
       Write (Output_Unit_HC_Alert,*) "Error in Modify_Taus: Unsupported option for hc_lambda_calc:",hc_lambda_calc
       Write (Output_Unit_HC_Alert,*) "Program stopping."
       Stop
    End If
 
    FTAU = CZENH * SQRT (CRMB) * RIZB * RIZB * Lambda * QTIM
    FTAUP = FTAU * 6.8E-14
    FTAUS = FTAU * 6.8E-14 * (1.0 + CRMB / CRMHC)
    FTAUT = FTAU * 1.4E-13
    C215A = FTAU * 1.4E-13 * SQRT (CRMHC)
    C350A = 9.0E-14 * (1.0 + CRMHC / CRMB) * RIZB * RIZB * CZENH * Lambda * QTIM / SQRT (CRMHC)
    C350B = 9.0E-14 * SQRT (CRMHC) * RIZB * RIZB * CZENH * Lambda * QTIM / CRMB
    RIZSQR = REAL (Current_HC_Charge) * REAL (Current_HC_Charge)
    ROOTTT = SQRT (CTEMAV)
    ROOTMI = SQRT (CRMHC)
    STAU = gknbs (Current_Cell,Current_Ring) * RIZSQR / (CRMHC * gktibs (Current_Cell,Current_Ring)**1.5)
    ! First, modify tau collision probability and tau parallel.
    ! KFPS(IK,IR,IZ) = STAU * KTIBS(IK,IR) * FTAUP * 2.0
    Local_HC_Change_State_Coll_Prob = STAU * gktibs (Current_Cell,Current_Ring) * FTAUP * 2.0
    !write (0,*) "Data modify taus",RIZB,CRMHC,CRMI,CRMB,IRSPEC,QTIM,CZENH,LAMBDA, Local_HC_Change_State_Coll_Prob
    !write (0,*) "More data",FTAU,FTAUP,FTAUS,FTAUT,C215A,C350A,C350B,RIZSQR,ROOTTT,ROOTMI,STAU,gknbs (Current_Cell,Current_Ring), &
    !& gktibs (Current_Cell,Current_Ring),Current_HC_Charge,Cur_HC_Spec
 
    ! Check for options effecting tau coll prob that depend on mass used in DIVIMP input.
    If (CIOPTB .eq. 1 .and. Current_Ring .ge. IRSpec) Then
       ! KFPS(IK,IR,IZ) = 0.0
       Local_HC_Change_State_Coll_Prob = 0.0
    ElseIf (CIOPTB .eq. 2 .and. Current_Ring .ge. IRSpec) Then
       ! KFPS(IK,IR,IZ) = 2.0 * KNBS(IK,IR) * 6.8E-14 * RIZB * REAL (CIZEFF) * RIZSQR * LAMBDA / (ROOTMI * ROOTTT) * QTIM
       Local_HC_Change_State_Coll_Prob = 2.0 * gknbs (Current_Cell,Current_Ring) * 6.8E-14 * RIZB * REAL ( Self_ZEff) * RIZSQR * &
&Lambda / (ROOTMI * ROOTTT) *  Ion_Time_Step
    ElseIf (CIOPTB .eq. 4 .and. Current_Ring .ge. IRSpec .and. CTEMAV .gt. gktibs (Current_Cell,Current_Ring) * CRMHC / CRMB) Then
       ! KFPS(IK,IR,IZ) = C350B * RIZSQR * KTIBS(IK,IR) * KNBS (IK,IR) / (CTEMAV * ROOTTT)
       Local_HC_Change_State_Coll_Prob = C350B * RIZSQR * gktibs (Current_Cell,Current_Ring) * gknbs (Current_Cell,Current_Ring) / &
&(CTEMAV * ROOTTT)
    End If
 
    ! From TAUIN2, line 90.
    If ( Local_HC_Change_State_Coll_Prob .eq. 0.0) Then
       ! KKKFPS(IK,IR,IZ) = 0.0
       Local_HC_Tau_Parallel_Inv = 0.0
    ElseIf ((CIOPTB .eq. 3 .or. CIOPTB .eq. 7) .and. Current_Ring .ge. IRSpec) Then
       ! KKKFPS(IK,IR,IZ) = SQRT (9.76E8 * CTEMAV / CRMI) * QTIM / KFPS(IK,IR,IZ)
       Local_HC_Tau_Parallel_Inv = SQRT (9.76E8 * CTEMAV / CRMHC) * QTIM /  Local_HC_Change_State_Coll_Prob
    ElseIf (CIOPTB .eq. 5) Then
       ! KKKFPS(IK,IR,IZ) = 1.56e4 * SQRT(CTEMAV/CRMI) * QTIM
       Local_HC_Tau_Parallel_Inv = 1.56e4 * SQRT (CTEMAV / CRMHC) * QTIM
    ElseIf (CIOPTB .eq. 6) Then
       ! KKKFPS(IK,IR,IZ) = 1.56e4 * SQRT(1.0/CRMI * KFPS(ik,ir,iz)/2.0) * QTIM
       Local_HC_Tau_Parallel_Inv = 1.56e4 * SQRT (1.0 / CRMHC *  Local_HC_Change_State_Coll_Prob / 2.0) * QTIM
    ElseIf (CIOPTB .eq. 8) Then
       ! KKKFPS(IK,IR,IZ) = QTIM * SQRT ((PI/4.0 * 4.88E8) /(KFPS(IK,IR,IZ)*CRMI))
       Local_HC_Tau_Parallel_Inv = QTIM * SQRT (( Pi_Value / 4.0 * 4.88E8) / ( Local_HC_Change_State_Coll_Prob * CRMHC))
    ElseIf (CIOPTB .eq. 9.and. Current_Ring .ge. IRSpec) Then
       ! KKKFPS(IK,IR,IZ)= SQRT(((PI/4.0)* 9.76E8)*CTEMAV/CRMI)* QTIM /  Local_HC_Change_State_Coll_Prob
       Local_HC_Tau_Parallel_Inv = SQRT ((( Pi_Value / 4.0) * 9.76E8) * CTEMAV / CRMHC)* QTIM /  Local_HC_Change_State_Coll_Prob
    ElseIf (CIOPTB .eq. 10) Then
       ! KKKFPS(IK,IR,IZ) = 1.56e4 * SQRT(PI/4.0 * CTEMAV/CRMI) * QTIM
       Local_HC_Tau_Parallel_Inv = 1.56e4 * SQRT ( Pi_Value / 4.0 * CTEMAV / CRMHC) * QTIM
    ElseIf (CIOPTB .eq. 11 .or. CIOPTB .eq. 12) Then
       ! KKKFPS(IK,IR,IZ) = 1.56e4 * SQRT(PI/4.0 * 1.0/CRMI *kfps(ik,ir,iz)/2.0) * QTIM
       Local_HC_Tau_Parallel_Inv = 1.56e4 * SQRT ( Pi_Value / 4.0 * 1.0 / CRMHC *  Local_HC_Change_State_Coll_Prob / 2.0) * QTIM
    ElseIf (CIOPTB .eq. 13) Then
       ! KKKFPS(IK,IR,IZ) = 1.56e4 * SQRT(PI/4.0 * 1.0/CRMI * (kfps(ik,ir,iz)*(1.0+CRMB/CRMI)) /2.0) * QTIM
       Local_HC_Tau_Parallel_Inv = 1.56e4 * SQRT &
            & ( Pi_Value / 4.0 * 1.0 / CRMHC * &
            & ( Local_HC_Change_State_Coll_Prob * &
            & (1.0 + CRMB / CRMHC)) / 2.0) * QTIM
 
       !write (0,*) "COLL opt 13", Local_HC_Tau_Parallel_Inv, Local_HC_Change_State_Coll_Prob
    Else
       ! KKKFPS(IK,IR,IZ) = QTIM * SQRT (4.88E8 /(KFPS(IK,IR,IZ)*CRMI))
       Local_HC_Tau_Parallel_Inv = QTIM * SQRT (4.88E8 /( Local_HC_Change_State_Coll_Prob * CRMHC))
    End If
 
    ! Next, modify tau stopping.
    TAU = STAU * FTAUS
    !write (0,*) "tau stopping",TAU,STAU,FTAUS
    If (TAU .gt.1.0E-3) Then
       ! KFSS(IK,IR,IZ) = 1.0-EXP(-TAU)
       Local_HC_Tau_Stopping_Inv = 1.0 - EXP (-TAU)
    Else
       ! KFSS(IK,IR,IZ) = TAU
       Local_HC_Tau_Stopping_Inv = TAU
    End If
    ! Check for options effecting tau stopping used in DIVIMP input.
    If (CIOPTC .eq. 1 .and. Current_Ring .ge. IRSpec) Then
       ! KFSS(IK,IR,IZ) = 0.0
       Local_HC_Tau_Stopping_Inv = 0.0
    ElseIf (CIOPTC .eq. 2 .and. Current_Ring .ge. IRSpec) Then
       TAU =  Local_HC_Change_State_Coll_Prob / (2.0 * CTEMAV)
       If (TAU .gt. 1.0E-3) Then
          ! KFSS(IK,IR,IZ) = 1.0-EXP(-TAU)
          Local_HC_Tau_Stopping_Inv = 1.0 - EXP (-TAU)
       Else
          ! KFSS(IK,IR,IZ) = TAU
          Local_HC_Tau_Stopping_Inv = TAU
       End If
    ElseIf (CIOPTC .eq. 3 .and. Current_Ring .ge. IRSpec .and. CTEMAV .gt. gktibs (Current_Cell,Current_Ring)*CRMHC / CRMB) Then
       TAU = C350A * RIZSQR * gknbs (Current_Cell,Current_Ring) / (CTEMAV * ROOTTT)
       If (TAU.Gt.1.0E-3) Then
          ! KFSS(IK,IR,IZ) = 1.0-EXP(-TAU)
          Local_HC_Tau_Stopping_Inv = 1.0 - EXP (-TAU)
       Else
          ! KFSS(IK,IR,IZ) = TAU
          Local_HC_Tau_Stopping_Inv = TAU
       End If
    End If
 
    ! Finally, modify tau heating.
    TAU = STAU * FTAUT
    !write (0,*) "tau h",TAU,STAU,FTAUT
    If (TAU .gt. 1.0E-3) Then
       ! KFTS(IK,IR,IZ) = 1.0-EXP(-TAU)
       Local_HC_Tau_Heating_Inv = 1.0-EXP(-TAU)
    Else
       ! KFTS(IK,IR,IZ) = TAU
       Local_HC_Tau_Heating_Inv = TAU
    End If
 
    ! Check for options effecting tau heating used in DIVIMP input.
    If (CIOPTD .eq. 1 .and. Current_Ring .ge. IRSpec) Then
       ! KFTS(IK,IR,IZ) = 0.0
       Local_HC_Tau_Heating_Inv = 0.0
    ElseIf (CIOPTD .eq. 2 .and. Current_Ring .ge. IRSpec) Then
       ! KFTS(IK,IR,IZ) = 1.0
       Local_HC_Tau_Heating_Inv = 1.0
    ElseIf (CIOPTD .eq. 3) Then
       TAU=C215A*RIZSQR*gknbs(Current_Cell,Current_Ring)/((CRMHC*gktibs(Current_Cell,Current_Ring)+CRMB*CTEMAV)**1.5)
       If (TAU.Gt.1.0E-3) Then
          ! KFTS(IK,IR,IZ) = 1.0-EXP(-TAU)
          !write (0,*) "-1ep(-Tauing):",1.0-EXP(-TAU)
          Local_HC_Tau_Heating_Inv = 1.0-EXP(-TAU)
       Else
          ! KFTS(IK,IR,IZ) = TAU
          !write (0,*) "Tauing:",tau
          Local_HC_Tau_Heating_Inv = TAU
       End If
    End If
 
    ! Note, ALPHAS requires no modifications (only charge dependent).
    ! Make modifications to BETAS for specific cell and impurity.
    If ( TIB_Coeff_Opt .eq. 0) Then
       Local_Betas = 0.0
    ElseIf ( TIB_Coeff_Opt .eq. 1 .or.  TIB_Coeff_Opt .eq. 3) Then
       Reduced_Mass = CRMHC / (CRMHC+CRMB)
       Local_Betas = -3.0*(1.0-Reduced_Mass-5.0* Root_2*(1.1*Reduced_Mass**2.5-0.35*Reduced_Mass**1.5)* RIZSQR) / (2.6 - 2.0*&
&Reduced_Mass + 5.4*Reduced_Mass*Reduced_Mass)
    ElseIf ( TIB_Coeff_Opt .eq. 2) Then
       Local_Betas =  Z_Charge * RIZSQR / (REAL ( ZO_Temp_Grad_Parameter) + SQRT (0.5*(1.0+CRMB/CRMHC)))
    EndIf
    !write (0,*) "ALPHS", TIB_Coeff_Opt, Local_Betas
  End Subroutine Modify_Taus
 
  Subroutine Adjust_Taus (Cur_HC_Spec,H_Isotope_Composition,Sput_Weight,HC_Temperature, &
       & Current_Cell,Current_Ring,Current_S,Launch_Reg,Eq_Total_Ion_Time_Steps,Seed,Random_Numbers_Used,NRand,hc_v)
 
    ! Required modules.
    Use ComHC ! Number_H_Species.
    Use HC_Init_DIV_Data ! Load cell and global properties data structures.
    Use HC_Init_DIV_Diag
    Use HC_Init_Lib_Data ! Get HC mass and charge functions.
    Use HC_Get
    use hc_kinetics
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Declare input.output variables.
    Integer, Intent (In) :: Cur_HC_Spec
    Integer, Dimension (Number_H_Species), Intent (In) :: H_Isotope_Composition ! Number of H/D/T's in the HC.
    !Integer, Dimension (Number_H_Products), Intent (In) :: H_Isotope_Composition ! Number of H/D/T's in the HC.
    Real, Intent (In) :: Sput_Weight
    Real, Intent (InOut) :: HC_Temperature
    Integer, Intent (In) :: Current_Cell
    Integer, Intent (In) :: Current_Ring
    Real, Intent (In) :: Current_S
    Integer, Intent (In) :: Launch_Reg ! M
    Double Precision, Intent (In) :: Eq_Total_Ion_Time_Steps
    Double Precision, Intent (In) :: Seed
    Integer, Intent (InOut) :: Random_Numbers_Used ! KK
    Integer, Intent (InOut) :: NRand

    type(hc_velocity_type1) :: hc_v
 
    ! Define local variables that will be assigned values.
    Real :: Current_HC_Mass
    Integer :: Current_HC_Charge
    Real :: RIZB ! Background charge.
    Real :: CRMHC ! Local mass of hydrocarbon molecule.
    Real :: CRMI ! Mass of DIVIMP impurity.
    Real :: CRMB ! Background plasma ion mass.
    Integer :: IrSpec ! Special ring where high-resolution measurement begins.
    Real :: CTEMAV ! Cell average temperature.
    Real :: QTIM ! Ion time step.
    Real :: CIOPTB ! Collision option.
    Real :: CIOPTC ! Friction option.
    Real :: CIOPTD ! Heating option.
    Real :: CZENH ! Z enhancement factor.
    Integer :: CDIFOP ! First diffuse option.
    ! Define local variables that will contain calculated values.		
    Real :: Tau_Mass_Factor ! Sqrt (CRMI/CRMHC) multiplication factor.
    Real :: C215A ! Intermediate step for TAU Heating calc.
    Real :: C350A ! Intermediate step for TAU Stopping calc.
    Real :: C350B ! Intermediate step for TAU Parallel calc.
    Real :: RIZSQR ! Current_HC_Charge^2.
    Real :: TAU ! Intermediate step for TAU Heating calc.
    Real :: Ran1 ! Temporary random number.
    Real :: Ran2 ! Temporary random number.
    Real :: Ratio1
    Real :: Ratio2
    Logical :: Debug_Adjust_Taus
 
    ! Assign local variables.
    Current_HC_Charge = Get_HC_Charge (Cur_HC_Spec)
    Current_HC_Mass = Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition)		
    Debug_Adjust_Taus = .False.
 
    ! Assign local variables from global data tables.
    RIZB =  Back_Plasma_Charge
    CRMHC = Current_HC_Mass
    CRMI =  Impurity_Ion_Mass
    CRMB =  Back_Plasma_Ion_Mass
    IRSpec =  Special_Ring
    CTEMAV = gctemav ()
    QTIM =  Ion_Time_Step
    CIOPTB =  Collision_Opt
    CIOPTC =  Friction_Opt
    CIOPTD =  Heating_Opt
    CZENH =  Z_Enh_Factor
    CDIFOP =  First_Diffuse_Opt
 
    ! PARALLEL DIFFUSION AND ION TEMPERATURE
 
    ! QUICKLY ADJUST RELEVANT CHARACTERISTIC TIMES
    ! FOR NON-STANDARD PLASMA OPTIONS.  BECAUSE OF THE SQUARE ROOTS
    ! THESE CALCULATIONS ARE DONE SPARINGLY TO PREVENT DRAMATIC
    ! SLOWING OF PROGRAM EXECUTION - FOR EXAMPLE, WHENEVER THE
    ! TEMPERATURE CHANGES BY 20% OR MORE.
 
    If (((CIOPTB .eq. 2 .or. CIOPTB .eq. 3 .or. CIOPTB .eq. 4 .or. CIOPTB .eq. 9 .or. CIOPTC .eq. 2 .or. CIOPTC .eq. 3) .and. &
&Current_Ring .ge. IRSpec) .or. CIOPTD .ge. 3 .or. CIOPTB .eq. 5 .or. CIOPTB .eq. 10) Then
       ! Warning:  These calculations depend on CTEMAV which is not dependable.
       Write (Output_Unit_HC_Alert,*) "Warning in Adjust_Taus:  These parallel transport calculations depend on CTEMAV which is &
&not dependable."
       If (HC_Temperature .gt. 1.2 * CTEMAV .or. HC_Temperature .lt. 0.8 * CTEMAV) Then
          If (CIOPTB .eq. 2) Then
             RATIO1 = SQRT (CTEMAV / HC_Temperature)
             RATIO2 = 1.0 / SQRT (RATIO1)
             Local_HC_Change_State_Coll_Prob = RATIO1 *  Local_HC_Change_State_Coll_Prob
             Local_HC_Tau_Parallel_Inv = RATIO2 *  Local_HC_Tau_Parallel_Inv
          ElseIf (CIOPTB .eq. 3 .or. CIOPTB .eq. 7 .or. CIOPTB .eq. 9 .or. CIOPTB .eq. 5 .or. CIOPTB .eq. 10) Then
             RATIO2 = SQRT (HC_Temperature / CTEMAV)
             Local_HC_Tau_Parallel_Inv = RATIO2 *  Local_HC_Tau_Parallel_Inv
          ElseIf (CIOPTB .eq. 4) Then
             If (HC_Temperature .gt. (gktibs (Current_Cell,Current_Ring) * CRMHC / CRMB)) Then
                Local_HC_Change_State_Coll_Prob = C350B * REAL ((Get_HC_Charge (Cur_HC_Spec))**2) * gktibs (Current_Cell,&
&Current_Ring) * gknbs (Current_Cell,Current_Ring) / ((HC_Temperature) ** 1.5)
                Local_HC_Tau_Parallel_Inv =  Ion_Time_Step * SQRT (4.88E8 / ( Local_HC_Change_State_Coll_Prob * CRMHC))
             End If
          End If
 
          If (CIOPTC .eq. 2) Then
             TAU =  Local_HC_Change_State_Coll_Prob / (2.0 * HC_Temperature)
             If (TAU .gt. 1.0E-3) Then
                Local_HC_Tau_Stopping_Inv = 1.0 - EXP (-TAU)
             Else
                Local_HC_Tau_Stopping_Inv = TAU
             End If
          ElseIf (CIOPTC .eq. 3) Then
             If (HC_Temperature .gt. gktibs (Current_Cell,Current_Ring) * CRMHC / CRMB) Then
                TAU = C350A * RIZSQR * gknbs (Current_Cell,Current_Ring) / (HC_Temperature ** 1.5)
                If (TAU .gt. 1.0E-3) Then
                   Local_HC_Tau_Stopping_Inv = 1.0 - EXP (-TAU)
                Else
                   Local_HC_Tau_Stopping_Inv = TAU
                End If
             End If
          End If
 
          If (CIOPTD .eq. 3) Then
             TAU=C215A*RIZSQR*gknbs(Current_Cell,Current_Ring)/((CRMHC*gktibs(Current_Cell,Current_Ring)+CRMB*HC_Temperature)**1.5)
             If (TAU .gt. 1.0E-3) Then
                Local_HC_Tau_Heating_Inv = 1.0 - EXP (-TAU)
             Else
                Local_HC_Tau_Heating_Inv = TAU
             End If
          End If
       End If
    End If
 
    ! Calculate parallel diffusion factor SPARA.
    If ( Ion_Diffuse) Then
       Local_SPara = HC_Temperature *  Local_HC_Tau_Parallel_Inv
    Else
       If (CDIFOP .eq. 2) Then
          If ( Local_HC_Change_State_Coll_Prob .gt. 0.0) Then
             Local_RConst = 2.0 * HC_Temperature /  Local_HC_Change_State_Coll_Prob
          Else
             Local_RConst =  Calc_Hi
          End If
       End If
       If (Eq_Total_Ion_Time_Steps .ge.  Local_RConst .or.  Local_RConst .lt. 1.0) Then
          HC_Tot_Time_First_Diff(Launch_Reg)=HC_Tot_Time_First_Diff(Launch_Reg)+Eq_Total_Ion_Time_Steps*Ion_Time_Step*Sput_Weight
          Ion_Diffuse = .True.
          Local_SPara  = HC_Temperature *  Local_HC_Tau_Parallel_Inv
       Else
          Local_SPara  = 0.0
       End If
    End If
    !write (0,*) "ion diffuse", Ion_Diffuse, Local_SPara,HC_Temperature, Local_HC_Tau_Parallel_Inv, &
    !&  Local_HC_Change_State_Coll_Prob, HC_Tot_Time_First_Diff (Launch_Reg),Launch_Reg
    ! Collision options 3,4,7,9: reset SPARA if necessary.
    If (Current_Ring .ge. IRSpec) Then
       If (CIOPTB .eq. 3 .or. CIOPTB .eq. 4 .or. CIOPTB .eq. 7 .or. CIOPTB .eq. 9) Then
          Random_Numbers_Used = Random_Numbers_Used + 1
          If (granv (Random_Numbers_Used) .gt.  Local_HC_Change_State_Coll_Prob / (2.0 * HC_Temperature)) Then
             Local_SPara = 0.0
          End If
       End If
    End If
 
    ! Collision option 5 - set parallel diffusive velocity step
    If (CIOPTB .eq. 5 .or. CIOPTB .eq. 10) Then
       Local_VPara = 0.0
       Local_SPara  = 0.0
       Random_Numbers_Used = Random_Numbers_Used + 1
       If (granv (Random_Numbers_Used) .le.  Local_HC_Change_State_Coll_Prob / (2.0*HC_Temperature)) Then
          Local_VPara =  Local_HC_Tau_Parallel_Inv
       End If
 
    ElseIf (CIOPTB .eq. 6 .or. CIOPTB .eq. 11) Then
       ! Collision option 6 ... changes on each time-step
       Local_SPara  = 0.0
       Local_VPara =  Local_HC_Tau_Parallel_Inv * SQRT ( Local_HC_Change_State_Coll_Prob / (2.0 * HC_Temperature))
 
    ElseIf (cioptb .eq. 12 .or. cioptb .eq. 13) Then
       ! Collision option 12 ... changes on each time-step
       ! Collision option 13 ... tau para / (1+Mb/Mi) ! ???
       Local_SPara  = 0.0
       Local_VPara =  Local_HC_Tau_Parallel_Inv ! ???
       Do
          NRand = NRand + 1
          Call Surand2 (Seed, 1, Ran1)
          If (Ran1 .ne. 0.0) Then
             Exit
          End If
       End Do
       NRand = NRand + 1
       Call Surand2 (Seed, 1, Ran2)
       Local_RGauss = SQRT (-2.0 * LOG (Ran1)) * COS (2.0 *  Pi_Value * Ran2)
       Local_VPara =  Local_VPara *  Local_RGauss
       !write (Output_Unit_Scratch,*) "paras hc", Local_HC_Tau_Parallel_Inv, Local_RGauss, Local_VPara
 
       !write (0,*) "paras hc", Local_HC_Tau_Parallel_Inv, Local_RGauss, Local_VPara,Ran1,Ran2
 
    ElseIf(CIOPTB.eq.14.and.(Current_S.gt.TGrad_Zero_Dist.and.Current_S.lt.(gksmaxs(Current_Ring)-TGrad_Zero_Dist)))Then
       Local_SPara  = 0.0
       Local_VPara  = 0.0
    End If
 
    If (Debug_Adjust_Taus) Then
       ! Call TAUCHK (Current_Cell,Current_Ring,Get_HC_Charge(Cur_HC_Spec), Local_SPara,HC_Temperature)
    End If


    ! Krieger IPP/07 - SUN compiler insists on 132 column limit
    if (debug_kinetics) &
      & write(6,'(a,i5,10g12.5)') &
	  &   'UPDATE HC TPERP1:',Get_HC_Charge (Cur_HC_Spec),hc_temperature,hc_v%t,hc_v%tperp,gktibs(Current_Cell,Current_Ring),&
      &   gktibs(Current_Cell,Current_Ring)-HC_v%tperp,Local_HC_Tau_Heating_Inv

 
    ! Increment temperature based on energy differential.
    !write (0,*) "Before:", Local_HC_Tau_Heating_Inv,HC_Temperature,gktibs (Current_Cell,Current_Ring)
    HC_Temperature=MAX(Calc_Lo,HC_Temperature+(gktibs(Current_Cell,Current_Ring)-HC_Temperature)*Local_HC_Tau_Heating_Inv)
    !write (0,*) "After:",HC_Temperature

    !
    ! Update particle temperatures at this point if required
    !
    if (hc_kinetics_opt.eq.1) then

       ! overall particle temperature - should be the same as hc_temperature

       hc_v%t = MAX(Calc_Lo,hc_v%t+(gktibs(Current_Cell,Current_Ring)-HC_v%t)*Local_HC_Tau_Heating_Inv)

       ! perpendicular particle temperature - evolves the same as the overall temperature

       hc_v%tperp = MAX(Calc_Lo,hc_v%tperp+(gktibs(Current_Cell,Current_Ring)-HC_v%tperp)*Local_HC_Tau_Heating_Inv)

       ! Krieger IPP/07 - SUN compiler insists on 132 column limit
       if (debug_kinetics) &
            & write(6,'(a,i5,10g12.5)')  & 
              & 'UPDATE HC TPERP2:',Get_HC_Charge (Cur_HC_Spec),hc_temperature,hc_v%t,hc_v%tperp,gktibs(Current_Cell,Current_Ring),&
              & gktibs(Current_Cell,Current_Ring)-HC_v%tperp,Local_HC_Tau_Heating_Inv

    endif


    ! Reiser code ion transport - gaussian value.
    If ( Reiser_Opt .eq. 1 .or.  Reiser_Opt .eq. 2) Then
       Do
          NRand = NRand + 1
          Call Surand2 (Seed,1,Ran1)
          If (Ran1 .ne. 0.0) Then
             Exit
          End If
       End Do
       NRand = NRand + 1
       Call Surand2 (Seed,1,Ran2)
       Local_RGauss = SQRT (-2.0 * LOG (Ran1)) * COS (2.0 *  Pi_Value * Ran2)
    End If
 
  End Subroutine Adjust_Taus
 
  Real Function Sheath_E_Field (Current_Cell,Current_Ring,S_Distance,TeBP,Te,Ti,n,miB,ZiB,fD,B_field,EDiv)
 
    ! Use required modules.
    Use ComHC ! Contains Output_Unit_HC_Alert.
    Use HC_Get ! Contains gkbfs, gksmaxs.
 
    ! Purpose: To calculate the near-target electric field (V) and
    ! overwrite existing DIVIMP calculated field if returned value
    ! is greater than the DIVIMP value.  This E-field magnitude is
    ! then used in the calculation for electric force acting on the
    ! impurity ion directed toward the target.  Model based on
    ! Brooks
 
    ! Declare input/output variables.
    Integer, Intent (In) :: Current_Cell ! Current location in ring, ik.
    Integer, Intent (In) :: Current_Ring ! Current location, ir.
    Real, Intent (In) :: S_Distance ! Distance to position of particle along S (m).
    Real, Intent (In) :: TeBP ! Background plasma electron temperature at plates (eV).
    Real, Intent (In) :: Te ! Background plasma electron temperature (eV).
    Real, Intent (In) :: Ti ! Background plasma ion temperature (eV).
    Real, Intent (In) :: n ! Background plasma density (1/m^3).
    Real, Intent (In) :: miB ! Background plasma ion mass (amu).
    Real, Intent (In) :: ZiB ! Background plasma ion charge (e).
    Real, Intent (In) :: fD ! Fraction of voltage drop in Debye sheath.
    Real, Intent (In) :: B_field ! Magnetic field (T)
    Real, Intent (In) :: Ediv ! DIVIMP calculated E (V/m).
 
 
    ! jdemod - NEVER use single character variable names except in trivially small subroutines - they are very difficult to find
    ! ALSO   - comhc includes the params common block which includes most of these constants
    !          ech, amu, eps0
    !
    ! Declare local variables.
    !Real :: e0 = 8.85E-12
    !Real :: e = 1.602E-19
    !Real :: mAMU = 1.67E-27
 
    Real :: Sheath_Voltage_Multiplier = -3.0 ! Multiple for target electron temperature.
    Real :: S_Distance_To_Target ! Distance from target along S factoring in which target it came from (m).
    Real :: Z_Distance ! Perpendicular distance to target (m).
    Real :: B_Angle_To_Target ! Angle between B and the target (degrees).
    Real :: Debye_Dist ! Background plasma Debye length (m).
    Real :: Larmor_Radius ! Background plasma Larmor radius (m).
    Real :: phi0 ! Sheath voltage drop (V).
    Real :: phi1 ! Phi coefficient.
    Real :: phi2 ! Phi coefficient.
    Real :: EOut ! E field strength calculated by Brooks model and compared to EDiv.

    ! Calculate S distance to closest target.
    S_Distance_To_Target = MIN (S_Distance, gksmaxs (Current_Ring) - S_Distance)
 
    ! Check that S is not too large for applicability.
    If (S_Distance_To_Target .gt. 1.0) Then
       Write (Output_Unit_HC_Alert, *) "Warning in Sheath_E_Field: S is too large (>1.0 m) for small angle assumption:", &
            & S_Distance_To_Target,Current_Cell,Current_Ring,S_Distance
    End If
 
    ! Check plasma background density.
    If (n .lt. 1.0E17 .or. n .gt. 1.0E25) Then
       Write (Output_Unit_HC_Alert, *) "Warning in Sheath_E_Field: Plasma density is out of range (<1E17 or >1E25):", n
    End If
 
    ! Check fraction of voltage drop in Debye sheath.
    If (fD .lt. 0.0 .or. fD .gt. 1.0) Then
       Write (Output_Unit_HC_Alert, *) "Error in Sheath_E_Field: fD is out of range (<0.0 or >1.0):", fD
       Write (Output_Unit_HC_Alert, *) "Program stopping."
       Stop
    End If
 
    ! Find angle between magnetic field B and target surface.  Note: ASIN returns a value in radians.
    B_Angle_To_Target = ASIN (1 / gkbfs (Current_Cell, Current_Ring))
 
    ! Check for a realistic angle.  Note: 0.5236 radians = 30 degrees.
    If (B_Angle_To_Target .lt. 0.0 .or. B_Angle_To_Target .gt. 0.5236) Then
       Write (Output_Unit_HC_Alert, *) "Warning in Sheath_E_Field: B_Angle_To_Target is out of range (<0.0deg or >30deg):", &
&B_Angle_To_Target
    End If
 
    ! Find perpendicular distance to target (in Z).
    Z_Distance = S_Distance_To_Target * SIN (B_Angle_To_Target)
 
    ! Check that Z is not too large.
    If (Z_Distance .gt. 0.2) Then
       Write (Output_Unit_HC_Alert, *) "Warning in Sheath_E_Field: Z is too large (>0.2 m) for small angle assumption:", Z_Distance
    End If
 
    ! Calculate Debye length.
    Debye_Dist = SQRT (2.0*eps0*ech*Te*Ti / (n*(Te+Ti)*ech**2.0))
 
    ! Calculate Larmor radius.
    Larmor_Radius = SQRT (ech*Ti / (miB*AMU)) / (ZiB*ech*B_field / (miB*AMU))
 
    ! Calculate sheath voltage drop.
    phi0 = Sheath_Voltage_Multiplier * TeBP
 
    ! Calculate phi1.
    phi1 = phi0*(EXP (-6.0*Debye_Dist/Larmor_Radius) - (1.0-fD)) / (EXP (-6.0*Debye_Dist/Larmor_Radius) - EXP (-3.0))
 
    ! Calculate phi2.
    phi2 = 	phi0 * (1.0-fD-EXP (-3.0))/ (EXP (-6.0*Debye_Dist/Larmor_Radius)-EXP (-3.0))
 
    ! Find E.
    EOut  = -phi1 / (2.0*Debye_Dist)*EXP (-Z_Distance / (2.0*Debye_Dist)) - phi2/Larmor_Radius * EXP (-Z_Distance/Larmor_Radius)
 
    ! Project E to B-field.
    EOut = EOut * SIN (B_Angle_To_Target)
    !write (0,*) "Efield:",Z_Distance,S_Distance_To_Target,gksmaxs(Current_Ring),B_Angle_To_Target,EOut
 
    !
    ! jdemod - there is a bug here since the SIGN code won't work when the DIVIMP E-field is zero. Need to find a better
    !          way to properly assign the sign for the electric field force. 
    !        - the sign convention in DIVIMP is the following - 
    !          forces towards the IK=1 target are NEGATIVE 
    !          forces towards the IK=NKS(IR) target are POSITIVE
    !
    !

    Eout = sign(Eout,real(current_cell-gnks(current_ring)/2))


    !
    ! Compare to E-field calculated by DIVIMP.  Only assign Brooks calc if value is > ABS (DIVIMP E-field).
    !
    If (abs(EOut) > ABS (EDiv)) Then
       Sheath_E_Field = Eout
       ntimes_sheath = ntimes_sheath + 1.0
    Else
       ! Use DIVIMP value for e-field instead of sheath value.
       Sheath_E_Field = EDiv
       ntimes_div = ntimes_div + 1.0
    End If


    !
    ! jdemod - temporary for debugging the sheath E-field
    !
 
    sheath_bin = int(Z_Distance/sheath_extent)+1
    sheath_bin = min(nsheath_bins,max(1,sheath_bin))

    sheath_zbin_data(sheath_bin,1) = sheath_zbin_data(sheath_bin,1) + 1.0
    sheath_zbin_data(sheath_bin,2) = sheath_zbin_data(sheath_bin,2) + abs(Eout)

    average_zsheath = average_zsheath + z_distance
    average_ssheath = average_ssheath + s_distance_to_target
    average_sheath_efield = average_sheath_efield + sheath_e_field
    average_sheath_bfield = average_sheath_bfield + b_field
    average_sheath_bangle = average_sheath_bangle + b_angle_to_target*raddeg
    average_div_efield = average_div_efield + ediv

  End Function Sheath_E_Field
 
  Real Function Find_Ion_Angle (Current_Cell,Current_Ring)
    ! Return angle from positive X axis (clockwise-positive) for the
    ! magnetic field in the particular cell.
 
    ! Required modules.
    Use HC_Get ! gkorpg.
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Declare input/output variables.
    Integer, Intent (In) :: Current_Cell
    Integer, Intent (In) :: Current_Ring
 
    ! Declare local variables.
    Real :: R_Start
    Real :: R_End
    Real :: Z_Start
    Real :: Z_End
    Integer :: Cell_IN
 
    Real, External :: ATAN2C
 
    ! Assign beginning and end points, depending on location.
    ! Starting location is always furthest from the nearest target.
    ! End location is always closest to the nearest target.
 
    Cell_IN = gkorpg (Current_Cell,Current_Ring)
 
    If (Current_Cell .gt. gnks (Current_Ring) / 2.0) Then
       ! Outer leg of DIIID grid.
       ! Starting point is between vertices 1 and 2.
       R_Start = 0.5 * (grvertp (1,Cell_IN) + grvertp (2,Cell_IN))
       Z_Start = 0.5 * (gzvertp (1,Cell_IN) + gzvertp (2,Cell_IN))
 
       ! End point is between vertices 3 and 4.
       R_End = 0.5 * (grvertp (3,Cell_IN) + grvertp (4,Cell_IN))
       Z_End = 0.5 * (gzvertp (3,Cell_IN) + gzvertp (4,Cell_IN))
 
    Else
       ! Inner leg sign of DIIID grid.
       ! Starting point is between vertices 3 and 4.
       R_Start = 0.5 * (grvertp (3,Cell_IN) + grvertp (4,Cell_IN))
       Z_Start = 0.5 * (gzvertp (3,Cell_IN) + gzvertp (4,Cell_IN))
 
       ! End point is between vertices 1 and 2.
       R_End = 0.5 * (grvertp (1,Cell_IN) + grvertp (2,Cell_IN))
       Z_End = 0.5 * (gzvertp (1,Cell_IN) + gzvertp (2,Cell_IN))
    End If
 
    ! Find angle.
    Find_Ion_Angle = atan2c (Z_End - Z_Start, R_End - R_Start)
 
  End Function Find_Ion_Angle
 
  Real Function Find_Vessel_Segment_Normal (Current_Cell,Current_Ring)
    ! Return angle of the normal off a target-facing vessel segment.
 
    ! Required modules.
    Use HC_Get ! gkorpg.
    Use HC_Init_DIV_Data ! PI.
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Declare input/output variables.
    Integer, Intent (In) :: Current_Cell
    Integer, Intent (In) :: Current_Ring
 
    ! Declare local variables.
    Real :: R_Start
    Real :: R_End
    Real :: Z_Start
    Real :: Z_End
    Integer :: Cell_IN
    Integer :: Vessel_Segment_Index
 
    Real, External :: ATAN2C
 
    ! Find index of cell.
    Vessel_Segment_Index = gkorpg (Current_Cell, Current_Ring)
 
    ! Check cell.
    If (Current_Cell .ne. 1 .and. Current_Cell .ne. gnks (Current_Ring)) Then
       Write (Output_Unit_HC_Alert,*) "Warning in Find_Vessel_Segment_Normal: Cell is not against a target location:",Current_Cell,&
&Current_Ring,gnks (Current_Ring)
    End If
 
    ! Find normal angle to vessel segment.
    If (Current_Cell .gt. gnks (Current_Ring) / 2.0) Then
       ! Outer leg of DIIID grid.
       ! Starting point is vertices 4.
       R_Start = grvertp (4,Vessel_Segment_Index)
       Z_Start = gzvertp (4,Vessel_Segment_Index)
 
       ! End point is vertices 3.
       R_End = grvertp (3,Vessel_Segment_Index)
       Z_End = gzvertp (3,Vessel_Segment_Index)	
 
       Find_Vessel_Segment_Normal = atan2c (Z_End - Z_Start, R_End - R_Start) +  Pi_Value / 2.0
    Else
       ! Inner leg sign of DIIID grid.
       ! Starting point is between vertices 3 and 4.
       R_Start = grvertp (2,Vessel_Segment_Index)
       Z_Start = gzvertp (2,Vessel_Segment_Index)
 
       ! End point is between vertices 1 and 2.
       R_End = grvertp (1,Vessel_Segment_Index)
       Z_End = gzvertp (1,Vessel_Segment_Index)
 
       Find_Vessel_Segment_Normal = atan2c (Z_End - Z_Start, R_End - R_Start) -  Pi_Value / 2.0
    End If
 
    !write (0,*) "VESSEL",Current_Cell,Current_Ring,Vessel_Segment_Index,R_Start,Z_Start,R_End,Z_End,Find_Vessel_Segment_Normal,gwallpt(Vessel_Segment_Index,8),gwallpt(Vessel_Segment_Index,9)
 
  End Function Find_Vessel_Segment_Normal
 
  Real Function Find_Angle_To_Normal (Current_Cell,Current_Ring,Current_Angle)
    ! Return angle of incidence off the normal of an
    ! incoming particle to a vessel segment.
 
    ! Required modules.
    Use HC_Get ! gnks.
    Use HC_Init_DIV_Data ! PI.
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Declare input/output variables.
    Integer, Intent (In) :: Current_Cell
    Integer, Intent (In) :: Current_Ring
    Real, Intent (In) :: Current_Angle
 
    ! Declare local variables.
    Real :: Wall_Segment_Normal_Angle ! TNORM
 
    ! Check cell.
    If (Current_Cell .ne. 1 .and. Current_Cell .ne. gnks (Current_Ring)) Then
       Write (Output_Unit_HC_Alert,*) "Warning in Find_Angle_To_Normal: Cell is not against a target location:",Current_Cell,&
&Current_Ring,gnks (Current_Ring)
    End If
 
    Wall_Segment_Normal_Angle = Find_Vessel_Segment_Normal (Current_Cell,Current_Ring)
 
    ! Calculate angle difference.
    Find_Angle_To_Normal = ABS (-Wall_Segment_Normal_Angle - Current_Angle)
 
  End Function Find_Angle_To_Normal
 
  Subroutine HC_Update_Walldep (IK,IR,Cur_HC_Spec,IDT,IDW,IWStart,IDType,Sput_Weight,Launch_Reg)
    ! jdemod - allow for the specification of pre-calculated wall or target indexes as appropriate (instead of just target idn->idt)
    ! Required modules.
    Use ComHC ! Contains Output_Unit_HC_Alert for printed warnings.
    Use HC_Get
    Use HC_Init_DIV_Data ! Global DIVIMP data.
    Use HC_Init_DIV_Diag ! HC diagnostic data.
 
    Implicit None
 
    ! UPDATE_WALLDEP:
 
    ! This routine records the ion impact with the wall segment
    ! closest to where the particle left the grid. If ID is not
    ! equal to zero this is treated as a specified target index
    ! and is mapped to the equivalent wall index through the
    ! wallindex array. If it is zero - the R,Z of the cell centre
    ! specified by the ik,ir indices is used to find the shortest
    ! distance to the nearest wall segment centre - this is then
    ! used for recording the exit wall segment for the ion.
 
    ! David Elder	Nov 5, 1998
    ! Modified by Adam McLean, November, 2002 for use with HC following.
 
    ! Declare input variables.
    Integer, Intent (In) :: IK
    Integer, Intent (In) :: IR
    Integer, Intent (In) :: Cur_HC_Spec
    Integer, Intent (In) :: IDT,idw ! jdemod - add idw
    Integer, Intent (In) :: IWStart
    Integer, Intent (In) :: IDType ! Same as NEUTTYPE in NEUT.
    Real, Intent (In) :: Sput_Weight
    Integer, Intent (In) :: Launch_Reg
 
    ! Declare local variables.
    Real Best,DSQ,R,Z
    Integer Ind,ID


      write(6,'(a,8i6,4g12.5)') 'HC_Update_walldep:',ik,ir,cur_hc_spec,idt,idw,&
     &         iwstart,idtype,launch_reg,sput_weight


 
    ! If a target segment is not specified find the wall segment centre
    ! closest to the cell centre of the particle exit.
    ! jdemod - added idw

    If (IDT .eq. 0.and.idw.eq.0) Then  
 
       R = gRS (IK,IR)
       Z = gZS (IK,IR)
 
       BEST = HI
       DSQ  = HI
       IND = 1
       Do ID = 1, Num_Wall_Points
          DSQ = (gWALLPT (ID,1) - R) ** 2 + (gWALLPT (ID,2) - Z) ** 2
          If (DSQ .lt. BEST) Then
             BEST = DSQ
             IND = ID
          End If
       End Do
 
       If (IND .lt. 1 .or. IND .gt.  Num_Wall_Points) Then
 
          Write (Output_Unit_HC_Alert,*) "Warning in HC_Update_WallDep: No wall found:",idt,ind,ik,ir
          HC_WallsI ( Num_Wall_Points + 1,Cur_HC_Spec) = &
               &  HC_WallsI ( Num_Wall_Points + 1,Cur_HC_Spec) + Sput_Weight
 
          If (IWStart .ge. 1 .and. IWStart .le.  Num_Wall_Points) Then
 
             HC_WTDep (IWStart,  Num_Wall_Points + 1, 1) = &
                  &  HC_WTDep (IWStart,  Num_Wall_Points + 1, 1) + Sput_Weight
          End If
       Else
          HC_WallsI (ind,Cur_HC_Spec) = &
               &  HC_WallsI (ind,Cur_HC_Spec) + Sput_Weight
 
          If (IWStart .ge. 1 .and. IWStart .le.  Num_Wall_Points) Then
 
             HC_WTDep (iwstart, ind, 1) = &
                  &  HC_WTDep (iwstart, ind, 1) + Sput_Weight
          End If
 
       End If
 
    Elseif (idt.ne.0) then 
      ! Target segment specified
       If (gwallindex (idt) .ne. 0) Then
 
          HC_WallsI (gwallindex (idt),Cur_HC_Spec) =  HC_WallsI (gwallindex (idt),Cur_HC_Spec) + Sput_Weight
 
          If (IWStart .ge. 1 .and. IWStart .le.  Num_Wall_Points) Then
 
             HC_WTDep (IWStart, gwallindex (idt),1) =   HC_WTDep (IWStart, gwallindex (idt),1) + Sput_Weight
          End If
       Else
          Write (Output_Unit_HC_Alert,'(a,3i5)') "Warning in HC_Update_WallDep: Wallsi: target?:",idt, gwallindex (idt)
          HC_WallsI ( Num_Wall_Points + 1,Cur_HC_Spec) = &
               &  HC_WallsI ( Num_Wall_Points + 1,Cur_HC_Spec) + Sput_Weight
 
          If (IWStart .ge. 1 .and. IWStart .le.  Num_Wall_Points) Then
 
             HC_WTDep (iwstart,  Num_Wall_Points + 1,1) = &
                  &  HC_WTDep (iwstart,  Num_Wall_Points + 1,1) + Sput_Weight
          End If
 
       End If
 

    elseif (idw.ge.1.and.idw.le.num_wall_points) then 
!
!      Wall segment specified 
!

             hc_wallsi(idw,cur_hc_spec) = hc_wallsi(idw,cur_hc_spec)+ sput_weight

             if (iwstart.ge.1.and.iwstart.le.num_wall_points) then 
                hc_wtdep(iwstart,idw,1) =  hc_wtdep(iwstart,idw,1) + sput_weight
             endif

    else 

          Write (Output_Unit_HC_Alert,'(a,3i5)') "Warning in HC_Update_WallDep: Wallsi: idw > num_wall_points?:",idw

          hc_wallsi(num_wall_points+1,cur_hc_spec) = hc_wallsi(num_wall_points+1,cur_hc_spec) + sput_weight

          if (iwstart.ge.1.and.iwstart.le.num_wall_points) then 
             hc_wtdep(iwstart,num_wall_points+1,1) =  hc_wtdep(iwstart,num_wall_points+1,1) + sput_weight
          endif

    endif 



    Return
 
  End Subroutine HC_Update_Walldep
 
  Subroutine Record_Evolve_Data (Cur_HC_Spec,Launch_Reg,Current_R,Current_Z,Current_Cell, &
       & Current_Ring,Current_S,Current_Cross,Launch_Angle,Current_Velocity,HC_Temperature,Sput_Weight)
    ! Record average data for an evolution event.
 
    ! Required modules.
    Use ComHC ! Input options.
    Use HC_Init_DIV_Diag ! HC_Average data storage.
    Use HC_Init_DIV_Data ! PI.
    Use HC_Init_Lib_Data ! Contains Get_HC_Charge function.
    Use HC_Get ! Contains gksmaxs.
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Declare input/output variables.
    Integer, Intent (In) :: Cur_HC_Spec ! This is the new species, after transition.
    Integer, Intent (In) :: Launch_Reg
    Real, Intent (In) :: Current_R
    Real, Intent (In) :: Current_Z
    Integer, Intent (In) :: Current_Cell
    Integer, Intent (In) :: Current_Ring
    Real, Intent (In) :: Current_S
    Real, Intent (In) :: Current_Cross
    Real, Intent (In) :: Launch_Angle
    Real, Intent (In) :: Current_Velocity
    Real, Intent (In) :: HC_Temperature
    Real, Intent (In) :: Sput_Weight
 
    ! Define local variables.
    Real :: Local_Angle ! Store altered angle.
 
    ! Note, HC_AvTPos and HC_AvNPos are determined with the average position
    ! found after following all molecules and reported in HC_Start.
    HC_R_Pos_At_Prod (Cur_HC_Spec,Launch_Reg) =  HC_R_Pos_At_Prod (Cur_HC_Spec,Launch_Reg) + Current_R
    HC_Z_Pos_At_Prod (Cur_HC_Spec,Launch_Reg) =  HC_Z_Pos_At_Prod (Cur_HC_Spec,Launch_Reg) + Current_Z
    !write (0,*) "Launch_Reg",Launch_Reg,Cur_HC_Spec,gksmaxs (Current_Ring),gksmaxs (Current_Ring) - Current_S,Current_S
    If (Get_HC_Charge (Cur_HC_Spec) .ne. 0) Then
       If (Launch_Reg .eq. HC_Launch_Reg_Target_1_Dist .or. Launch_Reg .eq. HC_Launch_Reg_Target_1_Pin .or. &
            Launch_Reg .eq. HC_Launch_Reg_Sput_Target_1 .or. Launch_Reg .eq. HC_Launch_Reg_Refl_Target_1) Then
          ! HC produced at outboard target.
          HC_S_Pos_At_Prod(Cur_HC_Spec,Launch_Reg)=HC_S_Pos_At_Prod(Cur_HC_Spec,Launch_Reg)+(gksmaxs(Current_Ring)-Current_S)
       ElseIf (Launch_Reg .eq. HC_Launch_Reg_Target_2_Dist .or. Launch_Reg .eq. HC_Launch_Reg_Target_2_Pin .or. &
            Launch_Reg .eq. HC_Launch_Reg_Sput_Target_2 .or. Launch_Reg .eq. HC_Launch_Reg_Refl_Target_2) Then
          ! HC produced at inboard target.
          HC_S_Pos_At_Prod (Cur_HC_Spec,Launch_Reg) =  HC_S_Pos_At_Prod (Cur_HC_Spec,Launch_Reg) + Current_S
       End If
    End If
    HC_Cross_Pos_At_Prod (Cur_HC_Spec,Launch_Reg) =  HC_Cross_Pos_At_Prod (Cur_HC_Spec,Launch_Reg) + Current_Cross
    If (Launch_Angle .lt. 0.0) Then
       ! Keep all angles positive, +x axis = 0.0, increasing ccw.
       Local_Angle = Launch_Angle + 2.0 *  Pi_Value
    Else
       Local_Angle = Launch_Angle
    End If
    HC_Angle_At_Prod (Cur_HC_Spec,Launch_Reg) =  HC_Angle_At_Prod (Cur_HC_Spec,Launch_Reg) + Local_Angle
    HC_Num_Fragments (Cur_HC_Spec,Launch_Reg) =  HC_Num_Fragments (Cur_HC_Spec,Launch_Reg) + 1.0
    HC_Velocity_At_Prod (Cur_HC_Spec,Launch_Reg) =  HC_Velocity_At_Prod (Cur_HC_Spec,Launch_Reg) + ABS (Current_Velocity)
    HC_Max_Velocity_At_Prod(Cur_HC_Spec,Launch_Reg)=MAX(HC_Max_Velocity_At_Prod(Cur_HC_Spec,Launch_Reg),ABS(Current_Velocity))
    HC_Temperature_At_Prod (Cur_HC_Spec,Launch_Reg) =  HC_Temperature_At_Prod (Cur_HC_Spec,Launch_Reg) + HC_Temperature
 
  End Subroutine Record_Evolve_Data
 
  Subroutine Record_Hydrogen_Release (H_Cell,H_Ring,Last_HC_Species,Reaction_Number,Cur_HC_Spec,H_Isotope_Composition,Seed,NRand)
 
    ! Purpose:  To record the position of any lost hydrogen during a hydrocarbon transition.
    !           This is based on the hydrogen atoms/molecules released during the reaction
    !           that takes place as loaded.
    !           Note:  The probability of D vs. T loss is currently independent of background
    !                  plasma, only the ratio of D to T remaining in the HC molecule and
    !                  an additional enhancement to encourage tritium release over D due to
    !                  its additional mass.
 
    ! Required modules.
    Use ComHC ! Input options.
    Use HC_Init_DIV_Diag ! H_Density array.
    Use HC_Init_DIV_Data ! HC_Average data storage.
    Use HC_Init_Lib_Data ! Contains Get_HC_Charge function.
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Declare input/output variables.
    Integer, Intent (In) :: H_Cell
    Integer, Intent (In) :: H_Ring
    Integer, Intent (In) :: Last_HC_Species
    Integer, Intent (In) :: Reaction_Number
    Integer, Intent (In) :: Cur_HC_Spec
    Integer, Dimension (Number_H_Species), Intent (InOut) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (InOut) :: H_Isotope_Composition
    Double Precision, Intent (In) :: Seed
    Integer, Intent (InOut) :: NRand
 
    ! Declare local variables.
    Integer :: Product_Number
    Integer :: H_Sum
    Integer :: H_State_Num
    Integer :: H_Content
    Integer :: H_Counter
    Integer :: H_Multiplier
    Integer :: H_Size
    Integer :: H_State
    Real :: Avg_H_Mass
    Real :: Mass_Of_H ! Must be 1.0, 2.0 or 3.0 (H/D/T).
    Real :: Random_Value
    Real :: Prob_For_Remove
    Real :: D_Fraction
 
    Real :: T_Release_Enh_Factor

    integer :: in
 
    ! Find total H's in the HC molecule.
    H_Sum = SUM (H_Isotope_Composition)
 
    ! Set number of each isotope states loaded in HC_Init_Lib_Data.
    H_Size =  INT (Num_H_States / Number_H_Species) ! Typically 4.
 
    ! Set Tritium release enhancement factor.
    T_Release_Enh_Factor = 2.0 ! Should possibly be set in ComHC or in the input file.
 
    ! Determine type of H that was lost given known
    ! background and current molecular composition.
    ! jdemod - the mass of the HC H isotopes is now input and can differ from the 
    !          background plasma - so use the input value here instead
    !Avg_H_Mass =  Back_Plasma_Ion_Mass
    Avg_H_Mass =  input_HC_H_mass

    ! First, check for added mass via H impact reactions.
    ! jdemod - check to see if reaction_type is p 
    !If (Reaction_Number .gt. Num_Electron_Reactions) Then
    if (hc_state_transform_table(reaction_number)%reaction_type.eq.'p') then 
       ! H reaction.  Add background plasma H isotope mass to hydrocarbon.
       ! First, check for mixture of H isotopes.
       If (MOD (Avg_H_Mass, 1.0) .ne. 0.0) Then
          ! Record H isotope production from the given reaction.
          ! Note:  Only dealing with D/T now.
          If (Avg_H_Mass .ge. 2.0 .and. Avg_H_Mass .le. 3.0) Then
             ! Plasma contains some tritium.
             ! Note:  This value is currently not used.
             D_Fraction = 1.0 - MOD (Avg_H_Mass, 2.0)
          Else
             ! Fractional H not between 2 and 3 are not currently supported.
             Write (Output_Unit_HC_Alert,*) "Error in Record_Hydrogen_Release: Fractional H mass not between mass 2 and 3 are not &
                                            &currently supported.",Avg_H_Mass,MOD (Avg_H_Mass, 1.0)
             Write (Output_Unit_HC_Alert,*) "Program stopping."
             Stop
          End If
 
          ! Call new random number and add to total.
          Call SURAND2 (Seed, 1, Random_Value)
          NRand = NRand + 1
 
          ! Deal with combination of both D and T.
          If (Random_Value .lt. D_Fraction) Then
             ! Remove a D.
             Mass_Of_H = 2.0
          Else
             ! Remove a T.
             Mass_Of_H = 3.0
          End If
       Else
          ! 100% background species.  No isotope determination required.
          Mass_Of_H = Avg_H_Mass
       End If
 
       ! Update H_Isotope_Composition.			
       H_Isotope_Composition (INT (Mass_Of_H)) = H_Isotope_Composition (INT (Mass_Of_H)) + 1.0
 
    End If
 
    ! Determine what H isotopes were lost given plasma background.
    Product_Number = 0
    Do
       Product_Number = Product_Number + 1
       ! Check number of hydrogen products (atoms or molecules).
       ! jdemod - bug - if there are 3 hydrogen products from the reaction then this code checks for
       !          a product number beyond the end of the array when it should exit instead. 
       ! jdemod - do not want to execute the test with a possible out of bounds array subscript
       if (product_number.gt.number_h_products) exit  
       If (HC_State_Transform_Table (Reaction_Number) % End_H_States (Product_Number) .ne. 0) Then

          !write(6,'(a,4i6,1x,3a)') 'Reaction:',reaction_number,product_number,&
          !             &HC_State_Transform_Table (Reaction_Number) % End_H_States (Product_Number),&
          !       &H_State_Table (HC_State_Transform_Table (Reaction_Number) % End_H_States (Product_Number)) % Hydrogen_Content,&
          !       &':',hc_state_transform_table(reaction_number)%reaction_desc,':'
 
          ! Determine hydrogen atom content in released H.
          H_Content = H_State_Table (HC_State_Transform_Table (Reaction_Number) % End_H_States (Product_Number)) % Hydrogen_Content
 
          ! Do loop for all H's in H atom/molecule released (typically up to a mass of 2 - H2/D2/DT/T2 are all possible).
          H_Counter = 0
          Do
             H_Counter = H_Counter + 1
 
             ! Check for mixture of H isotopes.
             If (MOD (Avg_H_Mass, 1.0) .ne. 0.0) Then
 
                ! Record H isotope production from the given reaction.
                ! Note:  Only dealing with D/T now.
                If (Avg_H_Mass .ge. 2.0 .and. Avg_H_Mass .le. 3.0) Then
                   ! Plasma contains some tritium.
                   ! Note:  This value is currently not used.
                   D_Fraction = 1.0 - MOD (Avg_H_Mass, 2.0)
                Else
                   ! Fractional H not between 2 and 3 are not currently supported.
                   Write (Output_Unit_HC_Alert,*) "Error in Record_Hydrogen_Release: Fractional H mass not between mass 2 and 3 &
                                               &are not currently supported.",Avg_H_Mass,MOD (Avg_H_Mass, 1.0)
                   Write (Output_Unit_HC_Alert,*) "Program stopping."
                   Stop
                End If
 
                ! Call new random number and add to total.
                Call SURAND2 (Seed, 1, Random_Value)
                NRand = NRand + 1
 
                ! Check that both isotopes exist in current molecule.
                If (H_Isotope_Composition (2) .eq. 0) Then
                   ! No D left, must be T.
                   Mass_Of_H = 3.0
                ElseIf (H_Isotope_Composition (3) .eq. 0) Then
                   ! No D left, must be D.
                   Mass_Of_H = 2.0
                ElseIf (H_Isotope_Composition (2) .eq. 0 .and. H_Isotope_Composition (3) .eq. 0) Then
                   ! Both D and T are zero.  Should definitely not happen.
                   Write (Output_Unit_HC_Alert,*) "Error in Record_Hydrogen_Release: Total of both D and T are zero when a &
                                      &reduction should take place: ",H_Isotope_Composition (2),H_Isotope_Composition (3)
                   Write (Output_Unit_HC_Alert,*) "Program stopping."
                   Stop
                Else			
                   ! Deal with combination of both D and T.
                   ! Note:  May want to use D_Fraction here in the future.
                   Prob_For_Remove = 1 / (H_Isotope_Composition (2) + T_Release_Enh_Factor * H_Isotope_Composition (3))
                   If (Random_Value .lt. H_Isotope_Composition (2) * Prob_For_Remove) Then
                      ! Remove a D.
                      Mass_Of_H = 2.0
                   Else
                      ! Remove a T.
                      Mass_Of_H = 3.0
                   End If
                End If
             Else
                ! 100% background species.  No isotope determination required.
                Mass_Of_H = Avg_H_Mass
             End If
 
             ! Update H_Isotope_Composition.			
             H_Isotope_Composition (INT (Mass_Of_H)) = H_Isotope_Composition (INT (Mass_Of_H)) - 1.0
 
             !write(6,'(a,g12.5,1x,3i6)') 'H-isotope:',mass_of_h,(h_isotope_composition(in),in=1,3)
             

             ! Check for negative total for isotope.
             If (H_Isotope_Composition (INT (Mass_Of_H)) .lt. 0) Then
                ! Serious error in H isotope determination.
                Write (Output_Unit_HC_Alert,*) "Error in Record_Hydrogen_Release: H Isotope total reduced below 0: ",INT (&
                                 &Mass_Of_H),H_Isotope_Composition (INT (Mass_Of_H))
                Write (Output_Unit_HC_Alert,*) "Program stopping."
                Stop
             End If
 
             ! Check for end of current atom/molecule released.
             If (H_Counter .eq. H_Content) Then
                Exit
             End If
          End Do
 
          ! Record hydrogen loss in H_Density array.
          If (Mass_Of_H .eq. 1.0) Then
             H_Multiplier = 1
          ElseIf (Mass_Of_H .eq. 2.0) Then
             H_Multiplier = 2
          ElseIf (Mass_Of_H .eq. 3.0) Then
             H_Multiplier = 3
          Else
             Write (Output_Unit_HC_Alert,*) "Error in Record_Hydrogen_Release: H isotope mass above 3 are not supported: ",Mass_Of_H
             Write (Output_Unit_HC_Alert,*) "Program stopping."
             Stop
          End If
 
          ! H states from 1 to 4.
          If (H_Counter .eq. 1) Then
             ! H or H+.
             If(H_State_Table(HC_State_Transform_Table(Reaction_Number)%End_H_States(Product_Number))%State_Charge.eq.0)Then
                ! H.
                H_State = 2 + H_Size * (H_Multiplier - 1)
             ElseIf(H_State_Table(HC_State_Transform_Table(Reaction_Number)%End_H_States(Product_Number))%State_Charge.eq.1)Then
                ! H+.
                H_State = 1 + H_Size * (H_Multiplier - 1)
             Else
                ! Higher H charges not possible.
                Write (Output_Unit_HC_Alert,*) "Error in Record_Hydrogen_Release: H Isotope charge cannot be above 1: ", &
                     & INT (Mass_Of_H),H_State_Table (HC_State_Transform_Table (Reaction_Number) % End_H_States (Product_Number)) &
                                           &% State_Charge
                Write (Output_Unit_HC_Alert,*) "Program stopping."
                Stop 				
             End If
          ElseIf (H_Counter .eq. 2) Then
             ! H2 or H2+.
             If(H_State_Table(HC_State_Transform_Table(Reaction_Number)%End_H_States(Product_Number))%State_Charge.eq.0)Then
                ! H2.
                H_State = 4 + H_Size * (H_Multiplier - 1)
             ElseIf(H_State_Table(HC_State_Transform_Table(Reaction_Number)%End_H_States(Product_Number))%State_Charge.eq.1)Then
                ! H2+.
                H_State = 3 + H_Size * (H_Multiplier - 1)
             Else
                ! Higher H2 charges not possible.
                Write (Output_Unit_HC_Alert,*) "Error in Record_Hydrogen_Release: H Isotope charge cannot be above 1: ", &
                     & INT (Mass_Of_H),H_State_Table (HC_State_Transform_Table (Reaction_Number) % End_H_States (Product_Number)) &
                                                         &% State_Charge
                Write (Output_Unit_HC_Alert,*) "Program stopping."
                Stop 				
             End If
          Else
             ! Higher H3 or higher not possible.
             Write(Output_Unit_HC_Alert,*)"Error in Record_Hydrogen_Release:H molecular mass above 2 are not possible:",H_Counter
             Write (Output_Unit_HC_Alert,*) "Program stopping."
             Stop 				
          End If
 
          ! Check that H_State is not above that allowed.
          If (H_State .gt. Num_H_States) Then
             Write (Output_Unit_HC_Alert,*) "Error in Record_Hydrogen_Release: H state is too high: ",H_State,Num_H_States
             Write (Output_Unit_HC_Alert,*) "Program stopping."
             Stop 				
          End If
 
          ! Finally, record H loss to plasma.
          !Write (Output_Unit_Scratch,'(A10,10I4)') "H_DENS:",H_Cell,H_Ring,H_State,Last_HC_Species,Reaction_Number,Cur_HC_Spec
          H_Density (H_Cell,H_Ring,H_State) =  H_Density (H_Cell,H_Ring,H_State) + 1.0
 
       Else
          ! No more H products for this reaction.
          Exit
       End If
    End Do
 
  End Subroutine Record_Hydrogen_Release
 
  Subroutine Record_DIVIMP_Data ()
 
    ! Purpose:  To write all information from HC following routine
    ! relevant to DIVIMP impurity particle following.
 
    ! Required modules.
    Use HC_Put ! Put/Add subroutines.
    Use HC_Init_DIV_Data ! DIVIMP specific data.
    Use HC_Init_DIV_Diag ! Diagnostic data.
    Use HC_Init_Lib_Data ! Diagnostic data.
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Declare local variables.
    ! jdemod 
    integer :: index_cnt
    Integer :: i
    Integer :: j
    Integer :: k
    Integer :: l
    Integer :: m
 
    ! DIVIMP variables written which are not IZ dependent in DIVIMP.
    ! Variables written to OUT in divstore:
    ! Vessel_Int - HC_Erosion -> NEROS.  Required for plots 83,89.
    ! Vessel_Int - HC_WallsE -> WALLSE.  Required for plots 815.
    ! Vessel_Int - HC_WallsE_I -> WALLSE_I.  Required for plots 815.
    ! Vessel_Int - HC_WallsI -> WALLSI.  Required for plots 815.
    ! Vessel_Int - HC_WallsN -> WALLSN.  Required for plots 815.
    ! Vessel_Int - HC_WtDep - > WTDEP.  Required for plots 817.
    ! Vessel_Int - HC_WtSource -> WTSOURCE.  Required for plots 98,817.
    ! Transport - HC_Leak_Density -> CLEAKN.  Required for plots 96.
    ! Transport - HC_Leak_Position -> CLEAKPOS.  Required for plots 96.
    ! Transport - HC_Ion_Core_Density -> NCORE.  Required for plots 161-170.
    ! Transport - HC_Ion_Trap_Density -> NTRAP.  Required for plots 161-170.
    ! Transport - HC_Ion_Edge_Denisty -> NEDGE.  Required for plots 161-170.
    ! Transport - HC_Ion_Divertor_Density -> NDIVERT.  Required for plots 161-170.
    ! Transport - HC_Ion_MSOL_Density -> NMSOL.  Required for plots 161-170.
 
    ! Transport - HC_ChemDen -> CHEMDEN.  Not currently used in OUT.
    ! Transport - HC_ChemIZS -> CHEMIZS.  Not currently used in OUT.
    ! Vessel_Int - HC_PromptDeps -> PROMPTDEPS.  Not currently used in OUT.
    ! DWELFS.  Required for plots 101-104.
 
    ! HC_Erosion,HC_PromptDeps.
    Do i = 1, Max_Target_Cells
       ! Do for each target element.
       ! Copy all of NEROS(i,3) to NEROS(i,2)
       Do j = 1,3,1
          ! Do for ions at NEROS(:,1) and neutrals at NEROS(:,3).
          ! Note: Results for all hydrocarbon species (mass SPUTY) are added at once.
          Call paneros (i,j,SUM ( HC_Erosion (i,j,:)))
       End Do
       !Write (0,*) "SAVING HC_NEROS:",i,SUM( HC_Erosion (i,1,:)),SUM( HC_Erosion (i,2,:)),SUM( HC_Erosion (i,3,:))
       Do j = 1,6,1
          Call papromptdeps (i,j,SUM ( HC_PromptDeps (i,j,:)))
       End Do
    End Do
    !Write (0,*) "TOTAL NEROS1:",SUM ( HC_Erosion (:,1,:))
    !Write (0,*) "TOTAL NEROS3:",SUM ( HC_Erosion (:,3,:))
 
    ! HC_WallsE, HC_WallsE_I, HC_WtDep.
    Do i = 1, Max_Points + 1
       ! Do for each wall element.
       Call pawallse (i,SUM ( HC_WallsE (i,:)))
       Call pawallse_i (i,SUM ( HC_WallsE_I (i,:)))
       Call pawallsi (i,SUM ( HC_WallsI (i,:)))
       Call pawallsn (i,SUM ( HC_WallsN (i,:)))
 
       Do j = 1, Max_Points
          Do k = 1,3,1
             Call pawtdep (j,i,k, HC_WtDep (j,i,k))
          End Do
       End Do
    End Do
 
    ! jdemod - base wtsource data in DIVIMP on only the C+ ion transport as is done with neut data
    !        - the HC code was collecting wtsource data for each state instead of particle so there
    !          were issues with data adding up - changed to once/particle
    !
    ! HC_WtSource.
    Do i = 1, Max_Points
       Do j = 1, Max_Rings
          Do k = 1,4,1
             Do l = 1,6,1
                Call pawtsource (i,j,k,l, HC_WtSource (i,j,k,l))
             End Do
          End Do
       End Do
    End Do
 
    ! HC_Ion_Core_Density,HC_Ion_Trap_Density,HC_Ion_Edge_Density,HC_Ion_Divert_Density,HC_Ion_MSOL_Density,HC_ChemDen,HC_ChemIzs.
    Do i = 1, Max_Cells_Per_Ring
       Do j = 1, Max_Rings
          Call pancore (i,j, HC_Ion_Core_Density (i,j))
          Call pantrap (i,j, HC_Ion_Trap_Density (i,j))
          Call panedge (i,j, HC_Ion_Edge_Density (i,j))
          Call pandivert (i,j, HC_Ion_Divertor_Density (i,j))
          Call panmsol (i,j, HC_Ion_MSOL_Density (i,j))
          Call pachemden (i,j,SUM ( HC_ChemDen (i,j,:)))
          Call pachemizs (i,j,SUM ( HC_ChemIzs (i,j,:)))	
       End Do
    End Do
 
    ! DIVIMP variables written which are IZ dependent in DIVIMP.
    ! Variables written to OUT in divstore:
    ! Transport - HC_DDLims -> DDLIMS.  Required for plots 89,91-95,101-104,105-106.
    ! Transport - HC_DDts -> DDts.  Required for plots 109-110.
    ! Vessel_Int - HC_Deposit -> DEPS.  Required for plots 82.
    ! Transport - HC_ELims -> ELIMS.  Required for plots 85.
    ! Transport - HC_Lims -> LIMS.  Required for plots 101-104.
    ! Vessel_Int - HC_Walls -> WALLS.  Required for plots 81,90.
    ! Transport - HC_Leak_Density -> CLEAKN.  Required for plots 96.
    ! Transport - HC_Ffi -> FFI.  Required for plots 575-580.
    ! Transport - HC_FCell -> FCELL.  Required for plots 575-580.
    ! Transport - HC_Fthi -> FTHI.  Required for plots 575-580.
    ! Transport - HC_Fvbg -> FVBG.  Required for plots 575-580.
    ! Transport - HC_Force_Diff -> DIFF.  Required for plots 575-580.
    ! Transport - HC_IonVelAvg -> VELAVG.  Required for plots 575-580.
    ! Transport - HC_TIZS.  Required for plots 131-132.
    ! Evolve - HC_IonizDat.  Required for DAT output file.
 
    ! ZEFFS.  Required for plots 91-95. Depends only on DDLIMS - post-processed by DIVIMP only.
    ! DWELTS.  Required for plots 101-104. Not used - okay.
    ! HPOWLS.  Required for plots 141-150.
    ! POWLS.  Required for plots 141-150.
    ! HLINES.  Required for plots 151-156.
    ! LINES.  Required for plots 151-156.
 
    ! HC_Walls,HC_DDLims,HC_ELims,HC_Lims,HC_TIZS,HC_DDts.
    Do i = 1, Max_Cells_Per_Ring
       Do j = 1, Max_Rings
          Do k = -1,0 ! Do only for primaries (-1) and primaries+secondaries+etc. (0).
             Call patizs (i,j,k, HC_TIZS (i,j,k))
          End Do
          !
          ! jdemod - I don't undestand why hc_ddts is being written to all of the charge states unless maximum_ionization_states 
          !          only applies to the hydrocarbons ... even then - hc_ddts must be organized by charge state
          !
          ! hc_ddts seems to be used one way in hc_follow and with different indexing in hc_inside_ion and neutral
          ! treat hc_ddts the same as hc_density

          !Do k = 0, Max_Ionization_States,1
             !Call paddts (i,j,k, HC_DDts (i,j,k)) ! Required for calculation of ZEFFS in DIVIMP.
             !Call paddlims (i,j,k, HC_DDLims (i,j,k))
          !End Do
          Do k = 1,Number_HC_Species
             l = Get_HC_Charge (k)
             Call pawalls (i,j,l, HC_Walls (i,j,k))
             Do m = 1, Max_Time_Slices
                ! Do for all time slices.
                Call palims (i,j,l,m, HC_Lims (i,j,k,m))
             End Do
          End Do
 
          ! Write charge state 0 and 1 data for carbon impurity to DDLIMS. Note: C+ = state 1, C0 = state 2.
          Call paddlims (i,j,0, HC_Density (i,j,2))
          Call paddlims (i,j,1, HC_Density (i,j,1))
          ! Write temperature data for state 0 and 1 for carbon impurity to DDLIMS. Note: C+ = state 1, C0 = state 2.
          Call paddts (i,j,0, HC_Density (i,j,2))
          Call paddts (i,j,1, HC_Density (i,j,1))
       End Do
       Do j = 1,3,1
          Do k = 1,Number_HC_Species
             l = Get_HC_Charge (k)
             Call paelims (i,j,l, HC_ELims (i,j,k))
          End Do
       End Do
    End Do
 
    ! Deposition: HC_Deposit.
    Do i = 1, Max_Target_Cells
       ! Do for each target element.
       Do k = 1,Number_HC_Species
          l = Get_HC_Charge (k)
          If (l .gt. 0) Then
             Call padeps (i,l, HC_Deposit (i,k))
          End If
       End Do
    End Do
 
    ! jdemod - this is wrong since the hc_leak_position is indexed by launch_region as well. 
    !
    ! Leakage: HC_Leak_Position.
    ! Note, HASLEAKED must also be saved and returned to DIVIMP for each particle launched.
    !Do i = 1,Max_Impurities ! Note:  Max_Impurities is found in ComHC.
    !   Do j = 1,2 ! Done for particles from both targets.
    !      Call pacleakpos (i,j, HC_Leak_Position (i,j))
    !   End Do
    !End Do
    index_cnt = 0  
    Do i = 1,number_regions ! Note:  number_regions is found in ComHC.
       Do j = 1,number_leaked_core(i) ! Done for particles from both targets.
          index_cnt= index_cnt+1
          do k = 1,2
             ! jdemod - i,j index was reversed - corrected to j,i
             Call pacleakpos (index_cnt,k, HC_Leak_Position (j,i,k))
          end do  
       End Do
    End Do
   
 
    ! Leakage: HC_Leak_Density.
    Do i = 1, Max_Points
       ! Do for each wall element.
       Do k = 1,Number_HC_Species
          l = Get_HC_Charge (k)
          If (l .gt. 0) Then
             Call pacleakn (i,l, HC_Leak_Density (i,k))
          End If
       End Do
    End Do
 
    ! jdemod
    ! These should NOT be loaded into the DIVIMP arrays since the forces
    ! will be overwritten or added to by the ion force calculations in 
    ! DIVIMP - this is possibly useful information but only in the context
    ! of the HC module itself - also - the code storing these quantities 
    ! is based on the charge state of the hydrocarbon fragment - thus is almost
    ! always 1 - and groups the forces for all ionized fragments together no 
    ! matter what their mass
    !
    ! Ion forces: HC_Ffi,HC_FCell,HC_Fthi,HC_Fvbg,HC_Force_Diff,HC_IonVelAvg.
    !Do i = 1, Max_Cells_Per_Ring
    !   Do j = 1, Max_Rings
    !      Do k = 1, Max_Ionization_States
    !         Call paffi (i,j,k, HC_Ffi (i,j,k))
    !         Call pafcell (i,j,k, HC_FCell (i,j,k))
    !         Call pafthi (i,j,k, HC_Fthi (i,j,k))
    !         Call pafvbg (i,j,k, HC_Fvbg (i,j,k))
    !         Call padiff (i,j,k, HC_Force_Diff (i,j,k))
    !         Call paionvelavg (i,j,k, HC_IonVelAvg (i,j,k))
    !      End Do
    !   End Do
    !End Do
 
    Do i = 1,2,1
       Do j = 1,2,1
          Do k = 1,2,1
             Do l = 1,2,1
                Do m = 1,5,1
                   Call paionizdat (i,j,k,l,m, HC_LIonizDat (i,j,k,l,m))
                End Do
             End Do
          End Do
       End Do
    End Do
 
    ! Variables not written to OUT in divstore, not IZ dependent in DIVIMP:
    ! Transport - HC_DDVoid -> DDVOID
    ! Transport - HC_MTCInf -> MTCINF
    ! Transport - HC_MTCTotCnt -> MTCTOTCNT
    ! Global_Prop_Table - Ion_Diffuse -> DIFFUS
 
    ! Variables not written to OUT in divstore, IZ dependent in DIVIMP:
    ! Transport - HC_CicUts -> CICUTS ! Number of each HC reaching cutoff time
    ! Transport - HC_CieIzs -> CIEIZS ! Time spent in each state.
    ! Transport - HC_CitIzs -> CITIZS ! Timesteps spent in each state.
    ! Transport - HC_Clll
    ! Transport - HC_ClllS
    ! Transport - HC_ClllX
    ! Transport - HC_Cmmm
    ! Transport - HC_CmmmS
    ! Transport - HC_CmmmX
    ! Transport - HC_Cnnn
    ! Transport - HC_CnnnK
    ! Transport - HC_CnnnKT
    ! Transport - HC_CnnnS
    ! Transport - HC_CnnnT
    ! Transport - HC_CnnnX
    ! Transport - HC_CnorgR
    ! Transport - HC_CnorgS
    ! Transport - HC_CnorgZ
    ! Transport - HC_CrtrcS
    ! Transport - HC_Csss
    ! Transport - HC_csssS
    ! Transport - HC_cxxx
    ! Evolve - HC_IonizDat -> IONIZDAT
 
  End Subroutine Record_DIVIMP_Data
 
  Subroutine Bin_Test (Array,Lower,Upper,Bins)
 
    ! Required modules.
    Use ComHC ! Contains Output_Unit_HC_Alert.
 
    ! Every good Fortran program has...		
    Implicit None
 
    ! Define input/output variables
    Real, Dimension (:), Intent (In) :: Array
    Real, Intent (In) :: Lower
    Real, Intent (In) :: Upper
    Integer, Intent (In) :: Bins
 
    ! Declare local variables.
    Integer :: Number ! Number of values in the array.
    Integer :: i,j ! Loop counters.
    Integer, Dimension (Bins) :: Distrib ! Contains summed occurances of each range.
    Real :: Width ! Numeric range for each bin.
    Real :: Sum_Up ! Total.
 
    ! Format statement for statistics print out.
10  Format (a,1x,f8.3,1x,a,1x,f8.3,1x,a,1x,i6,1x,a,1x,f8.1)
 
    Number = Size (Array)
    Width = (Upper - Lower) / Bins
 
    ! Find first non-zero element in input array.
    Do i = Number,1,-1
       If (Array (i) .ne. 0.0) Then
          Number = i
          Exit
       End If
    End Do
 
    Write (Output_Unit_HC_Alert,*) "Array Distribution"
    Write (Output_Unit_HC_Alert,*) "Number:",Number,"Bins:",Bins
 
    ! Begin loop to differentiate each array value in location 'bins' for successive divisions.
    i = 0
    Do
       i = i + 1
       Distrib (i) = 0
       ! Begin loop to search through every range to see if it is within the division 'bin'.
       j = 0
       Do
          j = j + 1
          If (Array(j) .le. (Lower+i*Width) .and. Array(j) .gt. (Lower+(i-1)*Width)) Then
             Distrib(i) = Distrib(i) + 1
          End If
          If (j .ge. Number) Then
             Exit
          End If
       End Do
       Write (Output_Unit_HC_Alert,10) "Lower:",Lower+(i-1)*Width,"Upper:",Lower+i*Width,"Number:",Distrib(i),"Percent:",REAL(&
&Distrib(i))/REAL(Number)*100.0
       If (i .ge. Bins) Then
          Exit
       End If
    End Do
    ! Sum up to be sure total is the number of runs.
    i = 0
    Sum_Up = 0
    Do
       i = i + 1
       Sum_Up = Sum_Up + Distrib(i)
       If (i .ge. Bins) Then
          Exit
       End If
    End Do
    Write (Output_Unit_HC_Alert,*) "Sum:",Sum_Up
 
  End Subroutine Bin_Test
 
End Module HC_Utilities
