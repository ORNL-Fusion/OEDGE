c     -*-Fortran-*-
c
c ======================================================================
c
c input.f
c
c ======================================================================
c
      SUBROUTINE User_InitializeOptions
      USE mod_sol28_global
      USE mod_user
      IMPLICIT none

      opt_user%u_mom = 0

      


      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE User_AllocateArrays
      USE mod_sol28_global
      USE mod_user
      IMPLICIT none

c      WRITE(0,*) 'NCELL IN SETUPARRAYS:',ncell
      ALLOCATE(stuff(ncell))      


      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE User_DeallocateArrays
      USE mod_sol28_global
      USE mod_user
      IMPLICIT none


      DEALLOCATE(stuff)      


      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE User_LoadOptions(fp,itag,buffer)
      USE mod_sol28_global
      USE mod_user
      IMPLICIT none

      INTEGER   fp,itag
      CHARACTER buffer*(*)

      SELECTCASE (buffer(3:itag-1))

        CASE('S U_MOM')
          CALL ReadOptionI(buffer,2,opt_user%u_mom)


        CASE DEFAULT
          CALL ER('User_LoadOptions','Unrecognized input tag',*99)
      ENDSELECT 

      RETURN
 99   WRITE(0,*) 'TAG: >'//buffer(3:itag-1)//'<'
      STOP
      END
c
c ======================================================================
c
c output.f
c
c ======================================================================
c
      SUBROUTINE User_GenerateOutputFiles
      USE mod_interface
      USE mod_sol28_params
      USE mod_sol28_global
      USE mod_user
      IMPLICIT none

      INTEGER ion,itube,icell,ic1,ic2

      INTEGER i1,ival(10)
      REAL val(10)

      INTEGER, ALLOCATABLE :: it(:),ic(:)

      ion = 1

      ALLOCATE(it(ncell))
      ALLOCATE(ic(ncell))

      CALL inOpenInterface('osm.test.idl')

      DO itube = 1, ntube
        ic1 = tube(itube)%cell_index(LO)
        ic2 = tube(itube)%cell_index(HI)
        DO icell = ic1, ic2
          it(icell) = itube
          ic(icell) = icell
        ENDDO

c....   Low index target:
        CALL inPutData(itube,'TUBE','none')
        CALL inPutData(0    ,'CELL','none')
        CALL inPutData(0.0,'s','m')
        CALL inPutData(tube(itube)%te(LO),'Te','eV')
        CALL inPutData(tube(itube)%ne(LO),'ne','m-3')
c....   Cells centres:
        CALL inPutData(it(ic1:ic2),'TUBE','none')
        CALL inPutData(ic(ic1:ic2),'CELL','none')
        CALL inPutData(cell(ic1:ic2)%s,'s','m')
        CALL inPutData(fluid(ic1:ic2,ion)%te,'Te','eV')
        CALL inPutData(fluid(ic1:ic2,ion)%ne,'ne','m-3')        
c....   High index target:
        CALL inPutData(itube,'TUBE','none')
        CALL inPutData(ic2+1,'CELL','none')
        CALL inPutData(tube(itube)%smax,'s','m')
        CALL inPutData(tube(itube)%te(HI),'Te','eV')
        CALL inPutData(tube(itube)%ne(HI),'ne','m-3')
      ENDDO

      CALL inCloseInterface 

      DEALLOCATE(it)
      DEALLOCATE(ic)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c source.f
c
c ======================================================================
c
      SUBROUTINE User_VolumeParIonSource(target,source)
      USE mod_sol28_solver
      USE mod_user
      IMPLICIT none
 
      INTEGER,INTENT(IN)  :: target
      REAL*8 ,INTENT(OUT) :: source(*)

      INTEGER ic1,ic2

      ic1 = icbnd1(target)
      ic2 = icbnd2(target)

      source(ic1:ic2) = 0.0D0

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE User_VolumeParRecSource(target,source)
      USE mod_sol28_solver
      USE mod_user
      IMPLICIT none
 
      INTEGER,INTENT(IN)  :: target
      REAL*8 ,INTENT(OUT) :: source(*)

      INTEGER ic1,ic2

      ic1 = icbnd1(target)
      ic2 = icbnd2(target)

      source(ic1:ic2) = 0.0D0

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE User_VolumeParAnoSource(target,source)
      USE mod_sol28_solver
      USE mod_user
      IMPLICIT none
 
      INTEGER,INTENT(IN)  :: target
      REAL*8 ,INTENT(OUT) :: source(*)

      INTEGER ic1,ic2

      ic1 = icbnd1(target)
      ic2 = icbnd2(target)

      source(ic1:ic2) = 0.0D0

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE User_ParticleSource(ion,target,source)
      USE mod_sol28_solver
      USE mod_user
      IMPLICIT none
 
      INTEGER,INTENT(IN)  :: ion,target
      REAL*8 ,INTENT(OUT) :: source(*)

      INTEGER ic1,ic2

      ic1 = icbnd1(target)
      ic2 = icbnd2(target)

      source(ic1:ic2) = 0.0D0

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE User_MomentumSource(ion,target,source)
      USE mod_sol28_solver
      USE mod_user
      IMPLICIT none
 
      INTEGER,INTENT(IN)  :: ion,target
      REAL*8 ,INTENT(OUT) :: source(*)

      INTEGER ic1,ic2

      ic1 = icbnd1(target)  ! These set the cell index bounds for half flux-tube...
      ic2 = icbnd2(target)

c      STOP 'HOLY CRAP'

      ! DO ic = ic1, ic2
      !   R = cell(ic)%cencar(1)
      !   Z = cell(ic)%cencar(2)
      !   S = cell(ic)%s
      !   V = cell(ic)%vol
      ! ENDDO  

      source(ic1:ic2) = 0.0D0

      !ico = tube%cell_index(1) - 1
      ! stuff(ic1+ico:ic2+ico) = source(ic1:ic2)

      RETURN
 99   STOP
      END
c ======================================================================
c
      SUBROUTINE User_VolumeEneRecSource(target,source)
      USE mod_sol28_solver
      USE mod_user
      IMPLICIT none
 
      INTEGER,INTENT(IN)  :: target
      REAL*8 ,INTENT(OUT) :: source(*)

      INTEGER ic1,ic2

      ic1 = icbnd1(target)
      ic2 = icbnd2(target)

      source(ic1:ic2) = 0.0D0

      RETURN
 99   STOP
      END
c ======================================================================
c
      SUBROUTINE User_VolumeEneIonSource(target,source)
      USE mod_sol28_solver
      USE mod_user
      IMPLICIT none
 
      INTEGER,INTENT(IN)  :: target
      REAL*8 ,INTENT(OUT) :: source(*)

      INTEGER ic1,ic2

      ic1 = icbnd1(target)
      ic2 = icbnd2(target)

      source(ic1:ic2) = 0.0D0

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE User_DumpData(fp,fname,title)
      USE mod_sol28_global
      USE mod_user
      IMPLICIT none

      INTEGER       fp
      CHARACTER*(*) fname,title 

      RETURN
 99   STOP
      END
c ======================================================================
c
c main.f
c
c ======================================================================
c
      SUBROUTINE User_SetupSolverOptions(count)
      USE mod_sol28_solver

      INTEGER, INTENT(IN):: count

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE User_MainLoop(condition,count)
      USE mod_sol28_global
      USE mod_user
      IMPLICIT none

      LOGICAL             :: condition
      INTEGER, INTENT(IN) :: count
      
      RETURN
 99   STOP
      END
