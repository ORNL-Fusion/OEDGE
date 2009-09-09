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

      REAL GetCs2

      INTEGER ion,itube,icell,ic1,ic2,itube1,itube2

      INTEGER i1,ival(10),target
      REAL val(10),rdum

      INTEGER, ALLOCATABLE :: it(:),ic(:)
      REAL   , ALLOCATABLE :: cs(:)

      ion = 1

      itube1 = 1
      itube2 = ntube
c      itube1 = 4
c      itube2 = 4  ! ntube

      ALLOCATE(it(ncell))
      ALLOCATE(ic(ncell))
      ALLOCATE(cs(ncell))

c...  ------------------------------------------------------------------
      CALL inOpenInterface('osm.idl.plasma')

      DO itube = itube1, itube2
        ic1 = tube(itube)%cell_index(LO)
        ic2 = tube(itube)%cell_index(HI)
        DO icell = ic1, ic2
          it(icell) = itube
          ic(icell) = icell
          cs(icell) = GetCs2(fluid(icell,ion)%te,fluid(icell,ion)%ti)
        ENDDO

c....   Low index target:
        target = LO
        CALL inPutData(itube,'TUBE','none')
        CALL inPutData(-1   ,'CELL','none')
        CALL inPutData(0.0,'s','m')
        CALL inPutData(-1.0,'sbnd1','m')
        CALL inPutData(-1.0,'sbnd2','m')
        CALL inPutData(tube(itube)%ne(target)    ,'ne','m-3')
        CALL inPutData(tube(itube)%vi(target,ion),'vi','m s-1')
        rdum = GetCs2(tube(itube)%te(target),tube(itube)%ti(target,ion))
        CALL inPutData(rdum                      ,'cs','m s-1')
        CALL inPutData(tube(itube)%te(target)    ,'te','eV')
        CALL inPutData(tube(itube)%ti(target,ion),'ti','eV')
c....   Cells centres:
        CALL inPutData(it(ic1:ic2),'TUBE','none')
        CALL inPutData(ic(ic1:ic2),'CELL','none')
        CALL inPutData(cell(ic1:ic2)%s,'s','m')
        CALL inPutData(cell(ic1:ic2)%sbnd(1),'sbnd1','m')
        CALL inPutData(cell(ic1:ic2)%sbnd(2),'sbnd2','m')
        CALL inPutData(fluid(ic1:ic2,ion)%ne,'ne','m-3')        
        CALL inPutData(fluid(ic1:ic2,ion)%vi,'vi','m s-1')        
        CALL inPutData(cs   (ic1:ic2)       ,'cs','m s-1')
        CALL inPutData(fluid(ic1:ic2,ion)%te,'te','eV')
        CALL inPutData(fluid(ic1:ic2,ion)%ti,'ti','eV')
c....   High index target:
        target = HI
        CALL inPutData(itube,'TUBE','none')
        CALL inPutData(-1   ,'CELL','none')
        CALL inPutData(tube(itube)%smax          ,'s' ,'m')
        CALL inPutData(-1.0,'sbnd1','m')
        CALL inPutData(-1.0,'sbnd2','m')
        CALL inPutData(tube(itube)%ne(target)    ,'ne','m-3')
        CALL inPutData(tube(itube)%vi(target,ion),'vi','m s-1')
        rdum = GetCs2(tube(itube)%te(target),tube(itube)%ti(target,ion))
        CALL inPutData(rdum                      ,'cs','m s-1')
        CALL inPutData(tube(itube)%te(target)    ,'te','eV')
        CALL inPutData(tube(itube)%ti(target,ion),'ti','eV')
      ENDDO
      CALL inCloseInterface 

c...  ------------------------------------------------------------------
      CALL inOpenInterface('osm.idl.sources')
      DO itube = itube1, itube2
        ic1 = tube(itube)%cell_index(LO)
        ic2 = tube(itube)%cell_index(HI)
c....   Cells centres:
        CALL inPutData(it(ic1:ic2),'TUBE','none')
        CALL inPutData(ic(ic1:ic2),'CELL','none')
        CALL inPutData(cell(ic1:ic2)%s,'s','m')
        CALL inPutData(fluid(ic1:ic2,ion)%parsrc,'parsrc','???')        
        CALL inPutData(fluid(ic1:ic2,ion)%momsrc,'momsrc','???')        
        CALL inPutData(fluid(ic1:ic2,ion)%enesrc,'enesrc','???')
        CALL inPutData(fluid(ic1:ic2,ion)%eneion,'eneion','???')
        CALL inPutData(fluid(ic1:ic2,ion)%eneano,'eneano','???')
        CALL inPutData(fluid(ic1:ic2,ion)%enisrc,'enisrc','???')
      ENDDO

      CALL inCloseInterface 

c...  ------------------------------------------------------------------
      CALL inOpenInterface('osm.idl.params')
      CALL inPutData(2.0,'flupar mass','amu')
      CALL inCloseInterface 


      DEALLOCATE(it)
      DEALLOCATE(ic)
      DEALLOCATE(cs)

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
c
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
