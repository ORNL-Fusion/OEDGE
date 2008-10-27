! -*-Mode:f90-*-
! HC_HC_Neutral_Transport.f90
! Go neutral tansport routine
! Move a particle around a box, bouncing off the walls and continuing with
! original velocity.
!
! Adam McLean
 
Module HC_Neutral_Transport
 
  ! Every good Fortran program has...
  Implicit None
 
Contains
 
  Subroutine Neutral_Move(Current_R,Current_Z,Current_Cell,Current_Ring,Current_Velocity, &
       & Current_Velocity_In_R,Current_Velocity_In_Z,Current_Angle,Cur_HC_Spec, &
       & H_Isotope_Composition,Last_R,Last_Z,Last_Cell,Last_Ring,Last_Velocity, &
       & Last_Velocity_In_R,Last_Velocity_In_Z,Last_Angle,HC_Temperature,Debug)
 
    Use ComHC ! Number_H_Species.
    Use HC_Init_DIV_Data ! Contains global geometry data.
    Use HC_Get ! Access external DIVIMP data.
    Use HC_Init_Lib_Data ! Get HC mass and charge functions.
 
    ! Every good Fortran program has...
    Implicit None
 
    ! Declare variables in the subroutine call and their intent.
    Real, Intent (InOut) :: Current_R
    Real, Intent (InOut) :: Current_Z
    Integer, Intent (InOut) :: Current_Cell
    Integer, Intent (InOut) :: Current_Ring
    Real, Intent (InOut) :: Current_Velocity
    Real, Intent (InOut) :: Current_Velocity_In_R
    Real, Intent (InOut) :: Current_Velocity_In_Z
    Real, Intent (InOut) :: Current_Angle
    Integer, Intent (In) :: Cur_HC_Spec
    Integer, Dimension (Number_H_Species), Intent (In) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (In) :: H_Isotope_Composition
    Real, Intent (Out) :: Last_R
    Real, Intent (Out) :: Last_Z
    Integer, Intent (Out) :: Last_Cell
    Integer, Intent (Out) :: Last_Ring
    Real, Intent (Out) :: Last_Velocity
    Real, Intent (Out) :: Last_Velocity_In_R
    Real, Intent (Out) :: Last_Velocity_In_Z	
    Real, Intent (Out) :: Last_Angle
    Real, Intent (InOut) :: HC_Temperature
    Logical, Intent (In) :: Debug
 
    Last_Angle = Current_Angle ! Note not change in angle is made while moving a neutral.
 
    ! Assign last values.
    Last_R = Current_R
    Last_Z = Current_Z
    Last_Cell = Current_Cell
    Last_Ring = Current_Ring
    Last_Velocity = Current_Velocity
 
    ! Update (R,Z) coordinates and array indicies.
    Current_R = Current_R + Current_Velocity_In_R ! Note velocities already multiplied by FSRATE upon launch and calculation of R/Z components.
    Current_Z = Current_Z + Current_Velocity_In_Z ! Note velocities already multiplied by FSRATE upon launch and calculation of R/Z components.
 
    ! Check to see if grid position has changed (returns new Cell and Ring for given R,Z).
    Call gridpos (Current_Cell,Current_Ring,Current_R,Current_Z,.false., Grid_Error)
 
    ! Check for neutral heating.
    If ( Impurity_Neutral_Vel_Opt .eq. 1 .and. (Current_Cell .ne. Last_Cell .or. Current_Ring .ne. Last_Ring)) Then
 
       ! Assign last value.
       Last_Velocity = Current_Velocity
 
       ! Assign new velocity
       Current_Velocity = 1.38e4 * SQRT (gktibs (Current_Cell,Current_Ring) / Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition))
 
       ! Assign new temperature
       HC_Temperature = Find_HC_Mass (Cur_HC_Spec,H_Isotope_Composition) * (Current_Velocity / 1.38E4) * (Current_Velocity / 1.38E4)
 
       ! Assign last values.
       Last_Velocity_In_R = Current_Velocity_In_R ! Note no change if the particle has not moved into another cell.
       Last_Velocity_In_Z = Current_Velocity_In_Z ! Note no change if the particle has not moved into another cell.
 
       ! Reset velocity components
       Current_Velocity_In_R = Last_Velocity_In_R * Current_Velocity / Last_Velocity
       Current_Velocity_In_Z = Last_Velocity_In_Z * Current_Velocity / Last_Velocity
 
       If (Debug) Then
          Write (Output_Unit_Scratch,'(a,2(1x,g13.5),2i5,2(1x,g13.5))') 'Debug neutral transport:',Current_Velocity_In_R,&
&Current_Velocity_In_Z,Current_Cell,Current_Ring,gktibs(Current_Cell,Current_Ring),Current_Velocity
       End If
 
    End If
 
  End Subroutine Neutral_Move
 
End Module HC_Neutral_Transport
