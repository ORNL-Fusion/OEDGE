! -*-Mode:f90-*-
! HC_Stack.f90
! Contains description, declaration, and operations for
! the stack implementation of carbon following during
! higher hydrocarbon breakup.
! 
! Adam McLean under Prof. Peter Stangeby
! Research Assistant: David Elder
! July, 2002

Module HC_Stack

  Use ComHC ! Find out how big of molecules are supported.

  ! Every good Fortran 90 program has...
  Implicit None	

  ! Keep type declarations in memory.
  Save

  ! Set integer parameter variables,
  Integer, Parameter :: Stack_Length = Highest_Carbon_Content * 10 ! Currently equals 30.
  ! One location for each carbon in the highest possible hydrocarbon 
  ! supported by data in the code, minus one which is handled by default.

  ! Declare stack entry table format.
  Type Carbon_Stack
     Integer :: HC_Species
     Real :: HC_R_Position
     Real :: HC_Z_Position
     Real :: HC_S_Position
     Real :: HC_Cross_Position
     Real :: HC_Energy
     Real :: HC_Velocity
     Real :: HC_Velocity_In_S
     Real :: HC_Theta
     Real :: HC_Direction
  End Type Carbon_Stack

  ! Initialize carbon stack.
  Type (Carbon_Stack) :: Stack (Stack_Length)

Contains

  Integer Function HC_Stack_Size ()
    ! Find and return current length of hydrocarbon stack.

    ! Every good Fortran 90 program has...		
    Implicit None

    ! Declare local variables.
    Integer i

    ! Initialize variables.
    HC_Stack_Size = 0

    Do i=1, Size(Stack), 1
       If (Stack (i) % HC_Species .ne. 0) Then
          HC_Stack_Size = HC_Stack_Size + 1
       Else
          ! Once found the first 0 entry, we're done.
          Exit
       End If
    End Do

  End Function HC_Stack_Size


  Subroutine HC_Push (Species,R_Position,Z_Position,S_Position,Cross_Position, &
       & Direction,Theta,Velocity,Velocity_In_S,Energy)
    ! Save single hydrocarbon species on the stack.

    ! Every good Fortran 90 program has...		
    Implicit None

    ! Declare input/output variables.
    Integer, Intent (In) :: Species
    Real, Intent (In) :: R_Position ! X position.
    Real, Intent (In) :: Z_Position ! Z position.
    Real, Intent (In) :: S_Position ! S position.
    Real, Intent (In) :: Cross_Position ! Cross position.
    Real, Intent (In) :: Direction ! Angle direction (degrees), 0 oriented on the +ve X axis.
    Real, Intent (In) :: Theta ! Theta value for non orthogonal transport.
    Real, Intent (In) :: Velocity ! Hydrocarbon velocity (m/s).
    Real, Intent (In) :: Velocity_In_S ! Hydrocarbon ion velocity (m/s).
    Real, Intent (In) :: Energy ! Total hydrcarbon energy (eV).

    ! Declare local variables.
    Integer :: i, j

    ! Check if stack has room for the number of species being added.
    If (Stack (Size(Stack)) % HC_Species .ne. 0) Then
       ! Stack is full.  Report an error and stop the program.
       Write (Output_Unit_HC_Alert,*) "Error: Carbon stack does not have &
            & sufficient space for species being added."
       Write (Output_Unit_HC_Alert,*) "Program stopping."
    End If

    ! Position HC to be pushed in order of size.

    Do i = 1, Size(Stack), 1
       If ((Stack (i) % HC_Species .ne. 0) .And. (Stack (i) % HC_Species .lt. Species)) Then
          ! Space is full and current HC in the space is smaller 
          ! than the one being pushed. Move to next position.
          Cycle
       Else ! Save the HC in the current spot
          ! Move existing HC to next available spot.
          Do j = Size(Stack), i+1, -1
             ! Work backwards from the end to the current spot.
             Stack (j) % HC_Species = Stack (j-1) % HC_Species
             Stack (j) % HC_R_Position = Stack (j-1) % HC_R_Position
             Stack (j) % HC_Z_Position = Stack (j-1) % HC_Z_Position
             Stack (j) % HC_S_Position = Stack (j-1) % HC_S_Position
             Stack (j) % HC_Cross_Position = Stack (j-1) % HC_Cross_Position
             Stack (j) % HC_Direction = Stack (j-1) % HC_Direction
             Stack (j) % HC_Theta = Stack (j-1) % HC_Theta
             Stack (j) % HC_Velocity = Stack (j-1) % HC_Velocity
             Stack (j) % HC_Velocity_In_S = Stack (j-1) % HC_Velocity_In_S
             Stack (j) % HC_Energy = Stack (j-1) % HC_Energy
          End Do
          ! Save HC.
          Stack (i) % HC_Species = Species
          Stack (i) % HC_R_Position = R_Position
          Stack (i) % HC_Z_Position = Z_Position
          Stack (i) % HC_S_Position = S_Position
          Stack (i) % HC_Cross_Position = Cross_Position
          Stack (i) % HC_Direction = Direction
          Stack (i) % HC_Theta = Theta
          Stack (i) % HC_Velocity = Velocity
          Stack (i) % HC_Energy = Energy
          Exit
       End If
    End Do

  End Subroutine HC_Push

  Subroutine HC_Pull (Species, R_Position, Z_Position, S_Position, Cross_Position, &
       & Direction, Theta, Velocity, Velocity_In_S, Energy)
    ! Load hydrocarbon species from the stack.  
    ! Note, only one hydrocarbon will be retrieved for any call.

    ! Every good Fortran 90 program has...		
    Implicit None

    ! Declare input/output variables.
    Integer, Intent (Out) :: Species
    Real, Intent (Out) :: R_Position ! X position.
    Real, Intent (Out) :: Z_Position ! Z position.
    Real, Intent (Out) :: S_Position ! S position.
    Real, Intent (Out) :: Cross_Position ! Cross position.
    Real, Intent (Out) :: Direction ! Angle direction (degrees), 0 oriented on the +ve X axis.
    Real, Intent (Out) :: Theta ! Theta value for non-orthogonal transport.
    Real, Intent (Out) :: Velocity ! Hydrocarbon velocity (m/s).
    Real, Intent (Out) :: Velocity_In_S ! Hydrocarbon ion velocity (m/s).
    Real, Intent (Out) :: Energy ! Total hydrcarbon energy or temperature (eV).

    ! Declare local variables.
    Integer :: i

    ! Check if stack is empty.
    If (Stack (1) % HC_Species .eq. 0) Then
       ! Stack is empty.  Report an error and stop the program.
       Write (Output_Unit_HC_Alert,*) "Error: Carbon stack is already empty."
       Write (Output_Unit_HC_Alert,*) "Program stopping."
    End If

    ! Pull First storage location into output variables.
    Species = Stack (1) % HC_Species
    R_Position = Stack (1) % HC_R_Position
    Z_Position = Stack (1) % HC_Z_Position
    S_Position = Stack (1) % HC_S_Position
    Cross_Position = Stack (1) % HC_Cross_Position
    Direction = Stack (1) % HC_Direction
    Theta = Stack (1) % HC_Theta
    Velocity = Stack (1) % HC_Velocity
    Velocity_In_S = Stack (1) % HC_Velocity_In_S
    Energy = Stack (1) % HC_Energy

    ! Move remaining hydrocarbons down to fill the space.
    Do i = 1, Size(Stack) - 1, 1
       If (Stack (i+1) % HC_Species .ne. 0) Then

          ! Found a filled position
          Stack (i) % HC_Species = Stack (i+1) % HC_Species
          Stack (i) % HC_R_Position = Stack (i+1) % HC_R_Position
          Stack (i) % HC_Z_Position = Stack (i+1) % HC_Z_Position
          Stack (i) % HC_S_Position = Stack (i+1) % HC_S_Position
          Stack (i) % HC_Cross_Position = Stack (i+1) % HC_Cross_Position
          Stack (i) % HC_Direction = Stack (i+1) % HC_Direction
          Stack (i) % HC_Theta = Stack (i+1) % HC_Theta
          Stack (i) % HC_Velocity = Stack (i+1) % HC_Velocity
          Stack (i) % HC_Velocity_In_S = Stack (i+1) % HC_Velocity_In_S
          Stack (i) % HC_Energy = Stack (i+1) % HC_Energy

          ! Zero out copied location.
          Stack (i+1) % HC_Species = 0
          Stack (i+1) % HC_R_Position = 0.0
          Stack (i+1) % HC_Z_Position = 0.0
          Stack (i+1) % HC_S_Position = 0.0
          Stack (i+1) % HC_Cross_Position = 0.0
          Stack (i+1) % HC_Direction = 0.0
          Stack (I+1) % HC_Theta = 0.0
          Stack (i+1) % HC_Velocity = 0.0
          Stack (i+1) % HC_Velocity_In_S = 0.0
          Stack (i+1) % HC_Energy = 0.0	

       Else
          ! Now zero out previous location.
          Stack (i) % HC_Species = 0
          Stack (i) % HC_R_Position = 0.0
          Stack (i) % HC_Z_Position = 0.0
          Stack (i) % HC_S_Position = 0.0
          Stack (i) % HC_Cross_Position = 0.0
          Stack (i) % HC_Direction = 0.0
          Stack (i) % HC_Theta = 0.0
          Stack (i) % HC_Velocity = 0.0
          Stack (i) % HC_Velocity_In_S = 0.0
          Stack (i) % HC_Energy = 0.0	

       End If

    End Do

  End Subroutine HC_Pull

End Module HC_Stack
