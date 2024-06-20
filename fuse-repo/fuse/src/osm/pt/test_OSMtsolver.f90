!-----------------------------------------------------------------------
! Test code for OSM time solver
!-----------------------------------------------------------------------

PROGRAM FRONTS

   use prec_const
   implicit none

   integer :: it, itDIAG, iz
   integer, parameter :: Ntdiag = 80
   integer, parameter :: dtdiag = 100
   integer, parameter :: Nz = 200
   integer, parameter :: Zi = 1
   integer, parameter :: Nfields = 4
   real(float), parameter :: dt = 1.E-6
   real(float), parameter :: e = 1.602E-19
   real(float), parameter :: mi = 3.345E-27
   real(float), dimension(Nz) :: zmesh, Bamp
   real(float), dimension(Nz,Nfields) :: Dcoef, Sources
   real(float), dimension(Nz,4) :: Xfield_prev
   real(float), dimension(Nz,Nfields) :: Xfield_next
   real(float) :: period
   real(float), dimension(2,Nfields) :: BCs

   !----------------------------------------
   ! Initializes arrays
   !----------------------------------------
   zmesh = (/((iz-1)*100._float/(Nz-1),iz=1,Nz)/)
   period = TW*zmesh(Nz)-zmesh(Nz-1)
   Bamp = 0.585+0.3*cos(FO*PI*zmesh/period+ON)
   Xfield_prev(1:Nz,1) = 1.E18!*(ON+HF*cos(FO*PI*zmesh/max(zmesh(Nz),period)))
   !Xfield_prev(95:106,1) = 2.E18_float
   Xfield_prev(1:Nz,3) = 10._float
   !Xfield_prev(114:120,3) = 20._float
   Xfield_prev(1:Nz,4) = 10._float
   !Xfield_prev(95:106,4) = 5._float
   Xfield_prev(1:Nz,2) = ZE
   !Xfield_prev(95:106,2) = 1.E4_float

   period = 0.

   call saveDATA_first(Xfield_prev,zmesh,Bamp,Nz)
! slmod begin
   call saveDATA_ascii(0,Xfield_prev,zmesh,Bamp,Nz)
! slmod end

   !----------------------------------------
   ! Solver loop
   !----------------------------------------
   itDIAG = 0
   print*,'Starting simulation'
   DO it = 1, Ntdiag*dtdiag
      Dcoef(1:Nz,1) = ZE
      Dcoef(1:Nz,2) = ZE
      Sources(1:Nz,1) = 1.E21_float
      Sources(1:Nz,2) = ZE
      BCs(1,1) = 1.E18_float
      BCs(2,1) = 1.E18_float
      BCs(1,2) = -sqrt((Xfield_prev(1,3)+Xfield_prev(1,4))/(mi/e))
      BCs(2,2) = sqrt((Xfield_prev(Nz,3)+Xfield_prev(Nz,4))/(mi/e))
      if (Nfields.EQ.4) then
         Dcoef(1:Nz,3) = 4.E6_float
         Dcoef(1:Nz,4) = 4.E5_float
         Sources(1:Nz,3) = 4.E22_float
         Sources(1:Nz,4) = 4.E22_float
         BCs(1,3) = 10._float
         BCs(2,3) = 10._float
         BCs(1,4) = 10._float
         BCs(2,4) = 10._float
      endif

      call OSMtsolver(zmesh, Bamp, Nz, dt, mi, Zi, Nfields, &
		Dcoef, Sources, BCs, &
		period, Xfield_prev, Xfield_next)

      Xfield_prev(1:Nz,1:Nfields) = Xfield_next(1:Nz,1:Nfields)

!   if (itDIAG.eq.799) then
!     write(0,*) 'nz',Nz
!     write(0,*) Xfield_prev(1:10,1)
!     stop
!   endif

      IF (mod(it,dtdiag) .EQ. 0) THEN
	 itDIAG = itDIAG+1
         print*,' itDIAG = ',itDIAG,' out of ',Ntdiag,' (',Ntdiag*dtdiag,')'
	 call saveDATA(Xfield_prev,Nz)
	 call saveDATA_ascii(it,Xfield_prev,zmesh,Bamp,Nz)
      ENDIF

   END DO

   write(0,*) 'nz',Nz
   write(0,*) Xfield_prev(1:10,1)

   print*,'Simulation successfull... use Matlab to visualize'

   stop
end

