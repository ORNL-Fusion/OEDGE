!------------------------------------------------------------
! File: OSMtsolver.f90
! Project: OSMtsolver
! Date: 31/03/2009
! Main file of OSM time solver for density, parallel momentum,
!	ion and electron energies parallel balance equation
!------------------------------------------------------------


!-----------------------------------
! Main time solver subroutine that is
! called by OSM.
! Inputs:
!	- zmesh(Nz): parallel mesh array
!	- Bamp(Nz): amplitude of the magnetic field
!	- Nz: number of points in zmesh
!	- dt: time step
!	- mi: ion mass
!	- Zi: ion charge
!	- Nfields: number of fields (and equations)
!		we want to evolve. If Nfields=2, a
!		fixed T version is used (equations
!		on N and V only, Te and Ti profiles
!		not evolving); if Nfields = 4, all
!		the fields are evolved
!	- Dcoef(Nz,Nfields):
!		diffusion coefficients for each of
!		the unknowns
!	- Sources(Nz,Nfields):
!		source terms for each of
!		the unknowns
!	- BCs(2,Nfields):
!		boundary conditions on each side
!		of the field line for each of the
!		unknowns 
!	- period: the period of the mesh if it is
!		periodical, 0 otherwise
!	- Xfield_prev(Nz,4): the four unknown fields
!		(N, V, Te, Ti) at the previous time step.
!		 Te and Ti are just a parameter if Nfields=2
!	- Xfield_next(Nz,Nfields): the unknown fields
!		(N, V, Te, Ti) at the next time step.
!		 Te and Ti are not returned if Nfields=2
!-----------------------------------
subroutine OSMtsolver(zmesh, Bamp, Nz, dt, mi, Zi, Nfields, &
		Dcoef, Sources, BCs, &
		period, Xfield_prev, Xfield_next)
   use prec_const
   implicit none

   integer :: Nz, Zi, Nfields
   real(float) :: dt, mi
   real(float), dimension(Nz) :: zmesh, Bamp
   real(float), dimension(Nz,Nfields) :: Dcoef, Sources
   real(float), dimension(2,Nfields) :: BCs
   real(float), dimension(Nz,4) :: Xfield_prev
   real(float), dimension(Nz,Nfields) :: Xfield_next
   real(float) :: period

   real(float), parameter :: e = 1.602E-19
   real(float), dimension(:,:), allocatable :: Linop
   real(float), dimension(Nfields*Nz) :: RHS
   real(float), dimension(Nz) :: dzBamp, dzVoverB
   integer, dimension(Nfields*Nz) :: ipiv
   integer :: ifield, info, iz, KL, KU
   real(float), dimension(Nfields*Nz,2*Nfields**2) :: Uwoodbury
   real(float), dimension(2*Nfields**2,Nfields*Nz) :: Vwoodbury
   real(float), dimension(2*Nfields**2,2*Nfields**2) :: Hwoodbury
   real(float), dimension(2*Nfields**2) :: tmpwoodbury
   integer, dimension(2*Nfields**2) :: ipiv2

! slmod begin
   real(float), dimension(Nz) :: dummy
! slmod end

   ! Size of the band matrix
   KL = 2*Nfields-1
   KU = 2*Nfields-1
   allocate(Linop(2*KL+KU+1,Nfields*Nz))
   Uwoodbury = ZE
   Vwoodbury = ZE
   Hwoodbury = ZE
   tmpwoodbury = ZE

   ! A quick check to see if we have Nfields = 2 or 4
   if (.NOT.((Nfields.EQ.2).OR.(Nfields.EQ.4))) then
      print*, 'OSMtsolver received illegal Nfields argument &
		different from 2 or 4. Exiting...'
      stop
   endif

   ! A quick check to see if the user wants a diffusion term on density
   if (.NOT.(sum(abs(Dcoef(1:Nz,1))).EQ.0)) then
      print*, 'Warning: diffusion coefficient  for density set to &
		non zero value. Inconsistency may exist in the model.'
   endif

   ! Zeroing operators and right-hand side
   Linop = ZE
   RHS = ZE

   ! We change for the variable set that is actually used in the solver
   !	N=>N, V=>Gamma, Te=>Pe, Ti=>Pi
   Xfield_prev(1:Nz,2) = Xfield_prev(1:Nz,2)*Xfield_prev(1:Nz,1)
   Xfield_prev(1:Nz,3) = Zi*Xfield_prev(1:Nz,3)*Xfield_prev(1:Nz,1)
   Xfield_prev(1:Nz,4) = Xfield_prev(1:Nz,4)*Xfield_prev(1:Nz,1)

   ! Derivative of the magnetic field (we'll need it several times)
   call calc_deriv(Bamp, Nz, zmesh, dzBamp, period)

   ! Let's start with the density conservation equation
   !***************************************************
   call add_dtterm(Xfield_prev(1:Nz,1),Nz,Nfields,1,Linop,KL,KU,RHS,dt)
   call add_diffterm(Xfield_prev(1:Nz,1),Dcoef(1:Nz,1),zmesh(1:Nz),Bamp(1:Nz), &
			Nz,Nfields,1,1,Linop,KL,KU,period,Uwoodbury,Vwoodbury)
   do iz = 1, Nz
      RHS(Nfields*(iz-1)+1) = RHS(Nfields*(iz-1)+1) &
		+Sources(iz,1)
   end do

   ! Momentum Balance equation
   !**************************
   call add_dtterm(Xfield_prev(1:Nz,2),Nz,Nfields,2,Linop,KL,KU,RHS,dt)

   call add_diffterm(Xfield_prev(1:Nz,2),Dcoef(1:Nz,2),zmesh(1:Nz),Bamp(1:Nz), &
			Nz,Nfields,2,2,Linop,KL,KU,period,Uwoodbury,Vwoodbury)

! slmod begin
   dummy = -Dcoef(1:Nz,2)*Xfield_prev(1:Nz,2)/Xfield_prev(1:Nz,1)
   call add_diffterm(Xfield_prev(1:Nz,2),dummy, &
			zmesh(1:Nz),Bamp(1:Nz),Nz,Nfields,2,1,Linop,KL,KU,period,Uwoodbury,Vwoodbury)
!
!   call add_diffterm(Xfield_prev(1:Nz,2),-Dcoef(1:Nz,2) &
!						*Xfield_prev(1:Nz,2)/Xfield_prev(1:Nz,1), &
!			zmesh(1:Nz),Bamp(1:Nz),Nz,Nfields,2,1,Linop,KL,KU,period,Uwoodbury,Vwoodbury)
! slmod end
   do iz = 1, Nz
      RHS(Nfields*(iz-1)+2) = RHS(Nfields*(iz-1)+2) &
		!-e/mi*(Xfield_prev(iz,3)+Xfield_prev(iz,4))*dzBamp(iz)/Bamp(iz) &
		+Sources(iz,2)/mi
   end do

   ! Pressure Balance equations
   !***************************
   if (Nfields.EQ.4) then

      ! A derivative we gonna need
! slmod begin
      dummy = Xfield_prev(1:Nz,2)/(Bamp(1:Nz)*Xfield_prev(1:Nz,1))
      call calc_deriv(dummy, Nz, zmesh, dzVoverB, period)
!
!      call calc_deriv(Xfield_prev(1:Nz,2)/(Bamp(1:Nz)*Xfield_prev(1:Nz,1)), Nz, zmesh, dzVoverB, period)
! slmod end

      ! Equation on Te (here on pe)
      call add_dtterm(Xfield_prev(1:Nz,3),Nz,Nfields,3,Linop,KL,KU,RHS,dt)

! slmod begin
      dummy = TW/TH*Dcoef(1:Nz,3)
      call add_diffterm(Xfield_prev(1:Nz,3),dummy,zmesh(1:Nz),Bamp(1:Nz), &
			Nz,Nfields,3,3,Linop,KL,KU,period,Uwoodbury,Vwoodbury)
!
!      call add_diffterm(Xfield_prev(1:Nz,3),TW/TH*Dcoef(1:Nz,3),zmesh(1:Nz),Bamp(1:Nz), &
!			Nz,Nfields,3,3,Linop,KL,KU,period,Uwoodbury,Vwoodbury)
! slmod end

! slmod begin
      dummy = -TW/TH*Dcoef(1:Nz,3)*Xfield_prev(1:Nz,3)/(Zi*Xfield_prev(1:Nz,1))
      call add_diffterm(Xfield_prev(1:Nz,3),dummy, &
 			zmesh(1:Nz),Bamp(1:Nz),Nz,Nfields,3,1,Linop,KL,KU,period,Uwoodbury,Vwoodbury)
!
!      call add_diffterm(Xfield_prev(1:Nz,3),-TW/TH*Dcoef(1:Nz,3) &
! 						*Xfield_prev(1:Nz,3)/(Zi*Xfield_prev(1:Nz,1)), &
! 			zmesh(1:Nz),Bamp(1:Nz),Nz,Nfields,3,1,Linop,KL,KU,period,Uwoodbury,Vwoodbury)
! slmod end
      do iz = 1, Nz
         RHS(Nfields*(iz-1)+3) = RHS(Nfields*(iz-1)+3) &
		-TW/TH*Xfield_prev(iz,3)*Bamp(iz)*dzVoverB(iz) &
		+TW/TH*Sources(iz,3)
      end do

      ! Equation on Ti (here on pi)
      call add_dtterm(Xfield_prev(1:Nz,4),Nz,Nfields,4,Linop,KL,KU,RHS,dt)
! slmod begin
      dummy = TW/TH*Dcoef(1:Nz,4)
      call add_diffterm(Xfield_prev(1:Nz,4),dummy,zmesh(1:Nz),Bamp(1:Nz), &
			Nz,Nfields,4,4,Linop,KL,KU,period,Uwoodbury,Vwoodbury)
!
!      call add_diffterm(Xfield_prev(1:Nz,4),TW/TH*Dcoef(1:Nz,4),zmesh(1:Nz),Bamp(1:Nz), &
!			Nz,Nfields,4,4,Linop,KL,KU,period,Uwoodbury,Vwoodbury)
! slmod end

! slmod begin
      dummy = -TW/TH*Dcoef(1:Nz,4)*Xfield_prev(1:Nz,4)/Xfield_prev(1:Nz,1)
      call add_diffterm(Xfield_prev(1:Nz,4),dummy, &
			zmesh(1:Nz),Bamp(1:Nz),Nz,Nfields,4,1,Linop,KL,KU,period,Uwoodbury,Vwoodbury)
!
!      call add_diffterm(Xfield_prev(1:Nz,4),-TW/TH*Dcoef(1:Nz,4) &
!						*Xfield_prev(1:Nz,4)/Xfield_prev(1:Nz,1), &
!			zmesh(1:Nz),Bamp(1:Nz),Nz,Nfields,4,1,Linop,KL,KU,period,Uwoodbury,Vwoodbury)
! slmod end
      do iz = 1, Nz
         RHS(Nfields*(iz-1)+4) = RHS(Nfields*(iz-1)+4) &
		-TW/TH*Xfield_prev(iz,4)*Bamp(iz)*dzVoverB(iz) &
		+TW/TH*Sources(iz,4)
      end do

   endif

   ! MUSCL scheme for advection terms
   !*********************************
    call add_MUSCL_advect(Xfield_prev(1:Nz,1:4),RHS, &
 			zmesh,Bamp,Nz,Nfields,mi,BCs,period)

   ! Adding boundary conditions (if non periodic, ie in the SOL)
   !***************************
   if (period.EQ.ZE) then
      ! Boundary condition on the flux
         ! iz = 1
         Linop(:,2) = ZE
         Linop(KL+KU+1,2) = ON
         Linop(KL+KU+1-1,2) = -BCs(1,2)
         RHS(2) = ZE
         ! iz = Nz
         Linop(:,Nfields*(Nz-1)+2) = ZE
         Linop(KL+KU+1,Nfields*(Nz-1)+2) = ON
         Linop(KL+KU+1-1,Nfields*(Nz-1)+2) = -BCs(2,2)
         RHS(Nfields*(Nz-1)+2) = ZE
      ! Boundary conditions on Te and Ti (here on p_e and p_i)
      if (Nfields.EQ.4) then
         ! iz = 1
         Linop(:,3) = ZE
         Linop(KL+KU+1,3) = ON
         Linop(KL+KU+1-2,3) = -BCs(1,3)
         RHS(3) = ZE
         Linop(:,4) = ZE
         Linop(KL+KU+1,4) = ON
         Linop(KL+KU+1-3,4) = -BCs(1,4)
         RHS(4) = ZE
         ! iz = Nz
         Linop(:,Nfields*(Nz-1)+3) = ZE
         Linop(KL+KU+1,Nfields*(Nz-1)+3) = ON
         Linop(KL+KU+1-2,Nfields*(Nz-1)+3) = -BCs(2,3)
         RHS(Nfields*(Nz-1)+3) = ZE
         Linop(:,Nfields*(Nz-1)+4) = ZE
         Linop(KL+KU+1,Nfields*(Nz-1)+4) = ON
         Linop(KL+KU+1-3,Nfields*(Nz-1)+4) = -BCs(2,4)
         RHS(Nfields*(Nz-1)+4) = ZE
      endif
   endif

   ! Solving all the system
   !***********************
   call DGBTRF(Nz*Nfields,Nz*Nfields,KL,KU,Linop,2*KL+KU+1,ipiv,info)
   if (period.EQ.0) then ! Straight resolution
      call DGBTRS('T',Nz*Nfields,KL,KU,1,Linop,2*KL+KU+1,ipiv,RHS,Nfields*Nz,info)
   else ! Resolution using Woodbury formula
      call DGBTRS('T',Nz*Nfields,KL,KU,2*Nfields**2,Linop,2*KL+KU+1,ipiv,UWoodbury,Nfields*Nz,info)
      call DGBTRS('T',Nz*Nfields,KL,KU,1,Linop,2*KL+KU+1,ipiv,RHS,Nfields*Nz,info)
      do ifield = 1, 2*Nfields**2
         Hwoodbury(ifield,ifield) = ON
      end do
      Hwoodbury = Hwoodbury + matmul(Vwoodbury,Uwoodbury)
      call DGETRF(2*Nfields**2,2*Nfields**2,Hwoodbury,2*Nfields**2,ipiv2,info)
      tmpwoodbury = matmul(Vwoodbury,RHS)
      call DGETRS('N',2*Nfields**2,1,Hwoodbury,2*Nfields**2, &
			ipiv2,tmpwoodbury,2*Nfields**2,info)
      RHS = RHS - matmul(Uwoodbury,tmpwoodbury)
   endif

   ! Switching back to the N, V, Te, Ti unknown system and storing result
   do iz = 1, Nz
      Xfield_next(iz,1) = RHS(Nfields*(iz-1)+1)
      Xfield_next(iz,2) = RHS(Nfields*(iz-1)+2)/RHS(Nfields*(iz-1)+1)
   end do
   if (Nfields.EQ.4) then
      do iz = 1, Nz
         Xfield_next(iz,3) = RHS(Nfields*(iz-1)+3)/(Zi*RHS(Nfields*(iz-1)+1))
         Xfield_next(iz,4) = RHS(Nfields*(iz-1)+4)/RHS(Nfields*(iz-1)+1)
      end do
   else
      Xfield_prev(1:Nz,3) = Xfield_prev(1:Nz,3)/(Zi*Xfield_prev(1:Nz,1))
      Xfield_prev(1:Nz,4) = Xfield_prev(1:Nz,4)/Xfield_prev(1:Nz,1)
   endif

end subroutine







!--------------------------------------------------
! Subroutine that adds diffusion terms to the
! linear system of the model, using a fully
! implicit scheme
!	- Xfield(Nz): the current value of the
!		field to be advanced
!	- D_X(Nz): the diffusion coefficient
!	- zmesh(Nz): array of cells positions
!	- Nz: number of cells
!	- Nfields: number of fields in the system
!		(important to fill the operator correctly)
!	- ieq: equation number to which the diffusion
!		term is added
!	- ifield: field number to which the term applies
!	- Linop(2*KL+KU+1,Nfields*Nz): operator in
!		band storage
!	- KL: number of lower sub-diagonals
!	- KU: number of upper sub-diagonals
!	- period: the period of the mesh if it is
!		periodical, 0 otherwise
!	- Uwoodbury(Nfields*Nz,2*Nfields**2): first
!		correction matrix for Woodbury
!		formula for non-band elements in case
!		of a periodical mesh (see Numerical Recipes)
!	- Vwoodbury(2*Nfields**2,Nfields*Nz): second
!		correction matrix for Woodbury
!		formula for non-band elements in case
!		of a periodical mesh (see Numerical Recipes)
!--------------------------------------------------
subroutine add_dtterm(Xfield, Nz, Nfields, ifield, Linop, KL, KU, RHS, dt)
   use prec_const
   implicit none

   integer :: Nz, Nfields, ifield, KL, KU
   real(float), dimension(Nz) :: Xfield
   real(float), dimension(2*KL+KU+1,Nfields*Nz) :: Linop
   real(float), dimension(Nfields*Nz) :: RHS
   real(float) :: dt

   integer :: iz, diagind,lineind

   diagind = KL+KU+1

   ! Loop on inner points
   do iz = 1, Nz
      lineind = Nfields*(iz-1)+ifield
      Linop(diagind,lineind) = Linop(diagind,lineind)+ON/dt
      RHS(lineind) = RHS(lineind)+Xfield(iz)/dt
   end do

end subroutine add_dtterm





!-----------------------------------
! Subroutine that adds the boundary
! conditions
!************************
! NOT USED CURRENTLY
!************************
!-----------------------------------
subroutine add_BCs(Linop,KL,KU,RHS,BCs,Nz,Nfields,ifield)
   use prec_const
   implicit none

   integer :: Nz, Nfields, ifield, KL, KU
   real(float), dimension(2,Nfields) :: BCs
   real(float), dimension(2*KL+KU+1,Nfields*Nz) :: Linop
   real(float), dimension(Nfields*Nz) :: RHS

   integer :: iz, diagind, lineind

   diagind = KL+KU+1

   ! iz = 1
   lineind = ifield
   Linop(:,lineind) = ZE
   Linop(diagind,lineind) = ON
   RHS(lineind) = BCs(1,ifield)

   ! iz = Nz
   lineind = Nfields*(Nz-1)+ifield
   Linop(:,lineind) = ZE
   Linop(diagind,lineind) = ON
   RHS(lineind) = BCs(2,ifield)

end subroutine add_BCs










!--------------------------------------------------
! Subroutine that adds diffusion terms to the
! linear system of the model, using a fully
! implicit scheme
!	- Xfield(Nz): the current value of the
!		field to be advanced
!	- D_X(Nz): the diffusion coefficient
!	- zmesh(Nz): array of cells positions
!	- Bamp(Nz): amplitude of the magnetic field
!	- Nz: number of cells
!	- Nfields: number of fields in the system
!		(important to fill the operator correctly)
!	- ieq: equation number to which the diffusion
!		term is added
!	- ifield: field number to which the term applies
!	- Linop(2*KL+KU+1,Nfields*Nz): operator in
!		band storage
!	- KL: number of lower sub-diagonals
!	- KU: number of upper sub-diagonals
!	- period: the period of the mesh if it is
!		periodical, 0 otherwise
!	- Uwoodbury(Nfields*Nz,2*Nfields**2): first
!		correction matrix for Woodbury
!		formula for non-band elements in case
!		of a periodical mesh (see Numerical Recipes)
!	- Vwoodbury(2*Nfields**2,Nfields*Nz): second
!		correction matrix for Woodbury
!		formula for non-band elements in case
!		of a periodical mesh (see Numerical Recipes)
!--------------------------------------------------
subroutine add_diffterm(Xfield, D_X, zmesh, Bamp, Nz, Nfields, ieq, ifield, Linop, KL, KU, &
		period, Uwoodbury, Vwoodbury)

   use prec_const
   implicit none

   integer :: Nz, Nfields, ieq, ifield, KL, KU
   real(float), dimension(Nz) :: Xfield, D_X, zmesh, Bamp
   real(float), dimension(2*KL+KU+1,Nfields*Nz) :: Linop
   real(float), dimension(Nfields*Nz,2*Nfields**2) :: Uwoodbury
   real(float), dimension(2*Nfields**2,Nfields*Nz) :: Vwoodbury
   real(float) :: period

   integer :: iz, lineind, diagind
   real(float), dimension(1:Nz+1) :: Bamp_edge
   real(float) :: tempvol

   diagind = KL+KU+1

   ! Value of the magnetic field at cell edges
   do iz = 2, Nz
      Bamp_edge(iz) = HF*(Bamp(iz-1)+Bamp(iz))
   end do
   if (period.EQ.ZE) then
      Bamp_edge(1) = HF*(TH*Bamp(1)-Bamp(2))
      Bamp_edge(Nz+1) = HF*(TH*Bamp(Nz)-Bamp(Nz-1))
   else
      Bamp_edge(1) = HF*(Bamp(1)+Bamp(Nz))
      Bamp_edge(Nz+1) = Bamp_edge(1)
   endif

   ! Loop on inner points
   do iz = 2, Nz-1
      lineind = Nfields*(iz-1)+ieq
      tempvol = Bamp(iz)*(Bamp_edge(iz+1)/Bamp_edge(iz)-ON)/log(Bamp_edge(iz+1)/Bamp_edge(iz))
      Linop(diagind+ifield-ieq-Nfields,lineind) = &
		Linop(diagind+ifield-ieq-Nfields,lineind) &
		-HF*tempvol/Bamp_edge(iz)*(D_X(iz)+D_X(iz-1))/(zmesh(iz)-zmesh(iz-1))
      Linop(diagind+ifield-ieq,lineind) = Linop(diagind+ifield-ieq,lineind) &
		+HF*tempvol/Bamp_edge(iz+1)*(D_X(iz+1)+D_X(iz))/(zmesh(iz+1)-zmesh(iz)) &
		+HF*tempvol/Bamp_edge(iz)*(D_X(iz)+D_X(iz-1))/(zmesh(iz)-zmesh(iz-1))
      Linop(diagind+ifield-ieq+Nfields,lineind) = &
		Linop(diagind+ifield-ieq+Nfields,lineind) &
		-HF*tempvol/Bamp_edge(iz+1)*(D_X(iz+1)+D_X(iz))/(zmesh(iz+1)-zmesh(iz))
   end do

   ! Preparing treatment of periodic boundary conditions with Woodbury formula
   if (period.NE.ZE) then
      ! iz = 1 side
      iz = 1
      lineind = ieq
      tempvol = Bamp(iz)*(Bamp_edge(iz+1)/Bamp_edge(iz)-ON)/log(Bamp_edge(iz+1)/Bamp_edge(iz))
      Linop(diagind+ifield-ieq,lineind) = Linop(diagind+ifield-ieq,lineind) &
		+HF*tempvol/Bamp_edge(iz+1)*(D_X(iz+1)+D_X(iz))/(zmesh(iz+1)-zmesh(iz)) &
		+HF*tempvol/Bamp_edge(iz)*(D_X(Nz)+D_X(1))/modulo(zmesh(1)-zmesh(Nz),period)
      Linop(diagind+ifield-ieq+Nfields,lineind) = &
		Linop(diagind+ifield-ieq+Nfields,lineind) &
		-HF*tempvol/Bamp_edge(iz+1)*(D_X(iz+1)+D_X(iz))/(zmesh(iz+1)-zmesh(iz))
      Uwoodbury(lineind,(ieq-1)*Nfields+ifield) = ON
      Vwoodbury((ieq-1)*Nfields+ifield,(Nz-1)*Nfields+ifield) = &
		-HF*tempvol/Bamp_edge(iz)*(D_X(Nz)+D_X(1))/modulo(zmesh(1)-zmesh(Nz),period)
      ! iz = Nz side
      iz = Nz
      lineind = Nfields*(Nz-1)+ieq
      tempvol = Bamp(iz)*(Bamp_edge(iz+1)/Bamp_edge(iz)-ON)/log(Bamp_edge(iz+1)/Bamp_edge(iz))
      Linop(diagind+ifield-ieq-Nfields,lineind) = &
		Linop(diagind+ifield-ieq-Nfields,lineind) &
		-HF*tempvol/Bamp_edge(iz)*(D_X(iz)+D_X(iz-1))/(zmesh(iz)-zmesh(iz-1))
      Linop(diagind+ifield-ieq,lineind) = Linop(diagind+ifield-ieq,lineind) &
		+HF*tempvol/Bamp_edge(iz+1)*(D_X(Nz)+D_X(1))/modulo(zmesh(1)-zmesh(Nz),period) &
		+HF*tempvol/Bamp_edge(iz)*(D_X(iz)+D_X(iz-1))/(zmesh(iz)-zmesh(iz-1))
      Uwoodbury(lineind,Nfields**2+(ieq-1)*Nfields+ifield) = ON
      Vwoodbury(Nfields**2+(ieq-1)*Nfields+ifield,ifield) =  &
		-HF*tempvol/Bamp_edge(iz+1)*(D_X(Nz)+D_X(1))/modulo(zmesh(1)-zmesh(Nz),period)
   endif

end subroutine add_diffterm






!-----------------------------------------------------------------
! This subroutine is specific to the OSMtsolver model
! It computes the contribution to the right-hand side of the
! linear system of the advection terms using an explicit
! TVD MUSCL scheme (flux-limiter method)
!	- Xfield(1:Nz,1:4): the four fields of the model
!		(density, flux, e pressure, i pressure)
!	- RHS(Nfields*Nz): the right-hand side of the linear
!		system to be completed with advection terms
!	- zmesh(Nz): the cell points in the parallel direction
!	- Bamp(Nz): the amplitude of the Magnetic field
!	- Nz: number of cell points
!	- Nfields: number of fields to be advanced (either 2=
!		isothermal or 4=non-isothermal)
!	- mi: mass of the considered ion specie
!	- BCs: the boundary conditions array
!	- period: period of the mesh if closed field line, 0 if open
!-----------------------------------------------------------------
subroutine add_MUSCL_advect(Xfield,RHS,zmesh,Bamp,Nz,Nfields,mi,BCs,period)
   use prec_const
   implicit none

   integer :: Nz, Nfields
   real(float), dimension(1:Nz) :: zmesh, Bamp
   real(float), dimension(1:Nz,1:4) :: Xfield
   real(float), dimension(Nfields*Nz) :: RHS
   real(float) :: mi
   real(float), parameter :: e = 1.602E-19
   real(float), dimension(2,Nfields) :: BCs
   real(float) :: period

   integer :: ifield, iz
   real(float), dimension(0:Nz+1) :: zmesh_extr
   real(float), dimension(1:Nz+1) :: Bamp_edge
   real(float), dimension(0:Nz+1,4) :: Xfield_extr
   real(float), dimension(1:Nz+1,1:4) :: Xfield_left,Xfield_right
   real(float), dimension(1:Nz+1,1:Nfields) :: flux_left,flux_right
   real(float), dimension(1:Nz+1,1:Nfields) :: flux_KT
   real(float), dimension(1:Nz,1:4) :: flux_lim
   real(float), dimension(1:Nz+1) :: spectral_rad

   ! Extrapolated fields (depends on whether the BCs are periodic or not)
   zmesh_extr(1:Nz) = zmesh(1:Nz)
   Xfield_extr(1:Nz,1:4) = Xfield(1:Nz,1:4)
   if (period.EQ.ZE) then
      zmesh_extr(0) = TW*zmesh(1)-zmesh(2)
      zmesh_extr(Nz+1) = TW*zmesh(Nz)-zmesh(Nz-1)
      Xfield_extr(0,1:4) = TW*Xfield(1,1:4)-Xfield(2,1:4)
      Xfield_extr(Nz+1,1:4) = TW*Xfield(Nz,1:4)-Xfield(Nz-1,1:4)
   else
      zmesh_extr(0) = zmesh(Nz)-period
      zmesh_extr(Nz+1) = zmesh(1)+period
      Xfield_extr(0,1:4) = Xfield(Nz,1:4)
      Xfield_extr(Nz+1,1:4) = Xfield(1,1:4)
   endif

   ! Value of the magnetic field at cell edges
   do iz = 2, Nz
      Bamp_edge(iz) = HF*(Bamp(iz-1)+Bamp(iz))
   end do
   if (period.EQ.ZE) then
      Bamp_edge(1) = HF*(TH*Bamp(1)-Bamp(2))
      Bamp_edge(Nz+1) = HF*(TH*Bamp(Nz)-Bamp(Nz-1))
   else
      Bamp_edge(1) = HF*(Bamp(1)+Bamp(Nz))
      Bamp_edge(Nz+1) = Bamp_edge(1)
   endif

   ! Computation of flux limiter condition
   do ifield = 1, 4
      do iz = 1, Nz
         if (Xfield_extr(iz+1,ifield).EQ.Xfield_extr(iz,ifield)) then
            if (Xfield_extr(iz,ifield).GT.Xfield_extr(iz-1,ifield)) then
               flux_lim(iz,ifield) = 1.E20
            elseif (Xfield_extr(iz,ifield).LT.Xfield_extr(iz-1,ifield)) then
               flux_lim(iz,ifield) = -1.E20
            else
               flux_lim(iz,ifield) = ON
            endif
         else
            flux_lim(iz,ifield) = max((Xfield_extr(iz,ifield)-Xfield_extr(iz-1,ifield)) &
			/(Xfield_extr(iz+1,ifield)-Xfield_extr(iz,ifield)),ZE)
         endif
      end do
   end do

   ! Computation of the flux limiter: Ospre limiter function
   flux_lim = (ON+HF)*(flux_lim**2+flux_lim)/(flux_lim**2+flux_lim+ON)

   ! What are the left and right extrapolated values in the MUSCL scheme?
   do ifield = 1, 4
      do iz = 2, Nz+1
         Xfield_left(iz,ifield) = Xfield_extr(iz-1,ifield) &
		+HF*flux_lim(iz-1,ifield)*(Xfield_extr(iz,ifield)-Xfield_extr(iz-1,ifield))
      end do
      do iz = 1, Nz
         Xfield_right(iz,ifield) = Xfield_extr(iz,ifield) &
		-HF*flux_lim(iz,ifield)*(Xfield_extr(iz+1,ifield)-Xfield_extr(iz,ifield))
      end do
   end do
   if (period.EQ.ZE) then
      do ifield = 1, 4
         Xfield_left(1,ifield) = Xfield_right(1,ifield)
         Xfield_right(Nz+1,ifield) = Xfield_left(Nz+1,ifield)
      end do
   else
      do ifield = 1, 4
         Xfield_left(1,ifield) = Xfield_extr(Nz,ifield) &
		+HF*flux_lim(Nz,ifield)*(Xfield_extr(1,ifield)-Xfield_extr(Nz,ifield))
         Xfield_right(Nz+1,ifield) = Xfield_extr(1,ifield) &
		-HF*flux_lim(1,ifield)*(Xfield_extr(2,ifield)-Xfield_extr(1,ifield))
      end do
   endif

   ! Fluxes at cell edges
   flux_left(1:Nz+1,1) = Xfield_left(1:Nz+1,2)
   flux_right(1:Nz+1,1) = Xfield_right(1:Nz+1,2)
   flux_left(1:Nz+1,2) = Xfield_left(1:Nz+1,2)**2/Xfield_left(1:Nz+1,1) &
		+(Xfield_left(1:Nz+1,3)+Xfield_left(1:Nz+1,4))*e/mi
   flux_right(1:Nz+1,2) = Xfield_right(1:Nz+1,2)**2/Xfield_right(1:Nz+1,1) &
		+(Xfield_right(1:Nz+1,3)+Xfield_right(1:Nz+1,4))*e/mi
   if (Nfields.EQ.4) then
      flux_left(1:Nz+1,3) = Xfield_left(1:Nz+1,3)*Xfield_left(1:Nz+1,2)/Xfield_left(1:Nz+1,1)
      flux_right(1:Nz+1,3) = Xfield_right(1:Nz+1,3)*Xfield_right(1:Nz+1,2)/Xfield_right(1:Nz+1,1)
      flux_left(1:Nz+1,4) = Xfield_left(1:Nz+1,4)*Xfield_left(1:Nz+1,2)/Xfield_left(1:Nz+1,1)
      flux_right(1:Nz+1,4) = Xfield_right(1:Nz+1,4)*Xfield_right(1:Nz+1,2)/Xfield_right(1:Nz+1,1)
   endif

   ! Computation of spectral radii
   do iz = 1, Nz+1
      spectral_rad(iz) = max(abs(Xfield_extr(iz,2)/Xfield_extr(iz,1) &
		-sqrt((Xfield_extr(iz,3)+Xfield_extr(iz,4))/(mi/e*Xfield_extr(iz,1)))), &
	abs(Xfield_extr(iz,2)/Xfield_extr(iz,1) &
		+sqrt((Xfield_extr(iz,3)+Xfield_extr(iz,4))/(mi/e*Xfield_extr(iz,1)))), &
	abs(Xfield_extr(iz-1,2)/Xfield_extr(iz-1,1) &
		-sqrt((Xfield_extr(iz-1,3)+Xfield_extr(iz-1,4))/(mi/e*Xfield_extr(iz-1,1)))), &
	abs(Xfield_extr(iz-1,2)/Xfield_extr(iz-1,1) &
		+sqrt((Xfield_extr(iz-1,3)+Xfield_extr(iz-1,4))/(mi/e*Xfield_extr(iz-1,1)))))
   end do

   ! Computation of fluxes in Kurganov and Tadmor scheme
   do ifield = 1, Nfields
      flux_KT(1:Nz+1,ifield) = HF*((flux_left(1:Nz+1,ifield)+flux_right(1:Nz+1,ifield)) &
		-spectral_rad(1:Nz+1)*(Xfield_right(1:Nz+1,ifield)-Xfield_left(1:Nz+1,ifield)))
   end do
   if (period.EQ.ZE) then
      flux_KT(1,1) = BCs(1,2)*Xfield_right(1,1)
      flux_KT(Nz+1,1) = BCs(2,2)*Xfield_left(Nz+1,1)
      flux_KT(1,2) = BCs(1,2)*Xfield_right(1,2)+(Xfield_right(1,3)+Xfield_right(1,4))*e/mi
      flux_KT(Nz+1,2) = BCs(2,2)*Xfield_left(Nz+1,2)+(Xfield_left(Nz+1,3)+Xfield_left(Nz+1,4))*e/mi
   endif

   ! We can now define the right-hand side of the system... pfuiii!
   do ifield = 1, Nfields
      do iz = 1, Nz
         if (Bamp_edge(iz+1).NE.Bamp_edge(iz)) then
            RHS(Nfields*(iz-1)+ifield) = RHS(Nfields*(iz-1)+ifield)+ &
! slmod begin
! gfortran not happy about the -Bamp -- is it legit?
		(-Bamp(iz))*(Bamp_edge(iz+1)/Bamp_edge(iz)-ON)/log(Bamp_edge(iz+1)/Bamp_edge(iz))* &
!
!		-Bamp(iz)*(Bamp_edge(iz+1)/Bamp_edge(iz)-ON)/log(Bamp_edge(iz+1)/Bamp_edge(iz))* &
! slmod end
		(flux_KT(iz+1,ifield)/Bamp_edge(iz+1)-flux_KT(iz,ifield)/Bamp_edge(iz)) &
		/(HF*(zmesh_extr(iz+1)-zmesh_extr(iz-1)))
         else
            RHS(Nfields*(iz-1)+ifield) = RHS(Nfields*(iz-1)+ifield)+ &
! slmod begin
		(-Bamp(iz))* &
!
!  		-Bamp(iz)* &
! slmod end
		(flux_KT(iz+1,ifield)/Bamp_edge(iz+1)-flux_KT(iz,ifield)/Bamp_edge(iz)) &
		/(HF*(zmesh_extr(iz+1)-zmesh_extr(iz-1)))
         endif
      end do
   end do

end subroutine add_MUSCL_advect






!-------------------------------------------------------------------------
! A very simple subroutine to compute the 2nd order first derivative of
! a 1D array with respect to another 1D array
!	- field(Nz): the field we want do derivate
!	- Nz: size of the array
!	- zmesh(Nz): the array with respect to which the field is
!		derivated
!	- dzfield(Nz): the derivative on output
!	- period: real variable characterizing the periodicity
!		if period=0: the direction is non periodic and uncentered
!			1st order discretization is used at the boundaries
!		if period#0: the direction is periodic with period being
!			the period
!-------------------------------------------------------------------------
subroutine calc_deriv(field, Nz, zmesh, dzfield, period)
   use prec_const
   implicit none

   integer :: Nz
   real(float), dimension(Nz) :: field, dzfield, zmesh
   real(float) :: period

   integer :: iz

   ! Out of boundaries
   do iz = 2, Nz-1
      dzfield(iz) = (field(iz+1)-field(iz-1))/(zmesh(iz+1)-zmesh(iz-1))
   end do

   ! At the boundaries
   if (period.EQ.ZE) then ! non periodic
      dzfield(1) = (field(2)-field(1))/(zmesh(2)-zmesh(1))
      dzfield(Nz) = (field(Nz)-field(Nz-1))/(zmesh(Nz)-zmesh(Nz-1))
   else ! periodic
      dzfield(1) = (field(2)-field(Nz))/MODULO(zmesh(2)-zmesh(Nz),period)
      dzfield(Nz) = (field(1)-field(Nz-1))/MODULO(zmesh(1)-zmesh(Nz-1),period)
   endif

end subroutine calc_deriv
