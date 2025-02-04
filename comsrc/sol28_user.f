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
      IMPLICIT none


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
      ! jdemod - ion argument not currently used
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
c
c ======================================================================
c
c solver_main.f
c
c ======================================================================
c
      SUBROUTINE User_SetupSolverOptions(count)
      USE mod_sol28_solver
      implicit none
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
      
c      CALL CalcCFGradients
c      CALL CalcPolGradients
c      CALL CalcDriftSources

      RETURN
 99   STOP
      END
