c     -*-Fortran-*-
c
c ======================================================================
c
c
      SUBROUTINE ExportTetrahedrons(fname)
      USE mod_interface
      USE mod_geometry
      USE mod_eirene06_locals
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'
      INCLUDE 'slout'

      CHARACTER, INTENT(IN) :: fname*(*)

      INTEGER GetNumberOfObjects

      INTEGER status,ndat,ik,ir,iobj,ion
      REAL*8  p(3)
      REAL, ALLOCATABLE :: tdata(:,:)

      ion = 1

      CALL LoadObjects(fname(1:LEN_TRIM(fname)),status)
      IF (status.EQ.-1) THEN
        WRITE(0,*) 'MESSAGE ExportTetrahedrons: File not found'
        RETURN
      ENDIF

c      CALL LoadGrid('osm.raw')

      ndat = GetNumberOfObjects('default')
      ALLOCATE(tdata(ndat,2))
      CALL LoadTriangleData(2,1,1,1,tdata(1,1),'default')  ! Atom density (from plot 987)
      CALL LoadTriangleData(3,1,1,1,tdata(1,2),'default')  ! Mol. density (from plot 987)

      CALL inOpenInterface('idl.tet_centroid',ITF_WRITE)
      DO iobj = 1, nobj
        CALL CalcCentroid(iobj,2,p)
        CALL inPutData(SNGL(p(1)),'X','m')                     
        CALL inPutData(SNGL(p(2)),'Y','m')                     
        CALL inPutData(SNGL(p(3)),'Z','m')                     
        IF (grp(obj(iobj)%group)%origin.EQ.GRP_MAGNETIC_GRID) THEN
          ik = obj(iobj)%index(IND_IK)
          ir = obj(iobj)%index(IND_IR)
c          WRITE(0,*) 'ind:',iobj,i
          CALL inPutData(knbs (ik,ir),'NE','m-3')        
          CALL inPutData(ktebs(ik,ir),'TE','eV')
        ELSE
          CALL inPutData(0.0         ,'NE','m-3')        
          CALL inPutData(0.0         ,'TE','eV')
        ENDIF
        CALL inPutData(tdata(iobj,1),'N_D' ,'m-3')
        CALL inPutData(tdata(iobj,2),'N_D2','m-3')
      ENDDO
      CALL inCloseInterface

      DEALLOCATE(tdata)
      CALL geoClean
c      CALL osmClean

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c
      SUBROUTINE GenerateEIRENEDataFiles
      USE mod_interface
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'
      INCLUDE 'slout'

      INTEGER   ik,ir,ike,index(MAXNKS),pos(MAXNKS),tube(MAXNKS)
      CHARACTER unit*10

      WRITE(0,*) 'IDL EIRENE DATA FILES'

      index = -1
      DO ik = 1, MAXNKS
        pos(ik) = ik
      ENDDO

      unit = 'ph m-3 s-1'

      CALL inOpenInterface('idl.fluid_eirene',ITF_WRITE)
      DO ir = 2, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        ike = nks(ir)
        IF (ir.LT.irsep) ike = ike - 1
        tube = ir - 1                      ! TUBE is set to the OSM fluid grid system, where
        IF (ir.GT.irwall) tube = tube - 2  ! the boundary rings are not present
        CALL inPutData(index(1:ike),'INDEX','N/A')                     
        CALL inPutData(pos  (1:ike),'POS'  ,'N/A')                     
        CALL inPutData(tube (1:ike),'TUBE' ,'N/A')                     
        CALL inPutData(kss(1:ike,ir),'S','m')
        CALL inPutData(pinion (1:ike,ir),'ION_NET','m-3 s-1')        
        CALL inPutData(pinrec (1:ike,ir),'REC_NET','m-3 s-1')        
        CALL inPutData(pinmp  (1:ike,ir),'MOM_NET','?')
        CALL inPutData(pinqe  (1:ike,ir),'QE_NET' ,'?')
        CALL inPutData(pinqi  (1:ike,ir),'QI_NET' ,'?')
        CALL inPutData(pinatom(1:ike,ir),'ATM_DENS','m-3')
        CALL inPutData(pinmol (1:ike,ir),'MOL_DENS','m-3')
        CALL inPutData(pinline(1:ike,ir,6,H_BALPHA),'BALMER_ALPHA',unit)
        CALL inPutData(pinline(1:ike,ir,6,H_BGAMMA),'BALMER_GAMMA',unit)
      ENDDO
      CALL inCloseInterface
      RETURN
 99   STOP
      END
c
c ======================================================================
c
c
      SUBROUTINE GenerateDIVIMPDataFiles(nizs,cizsc,crmi,cion,absfac,
     .                                   title)
      USE mod_interface
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slout'
      INCLUDE 'comgra'
      INCLUDE 'cgeom'
      INCLUDE 'walls_com'
      INCLUDE 'dynam2'
      INCLUDE 'dynam3'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER  , INTENT(IN) :: nizs,cizsc,cion
      REAL     , INTENT(IN) :: crmi,absfac
      CHARACTER, INTENT(IN) :: title*(*)

      REAL GetFlux

      INTEGER   ik,ir,iz,id,in,status,ike,target,fp,
     .          index(MAXNKS),pos(MAXNKS),tube(MAXNKS)
      REAL      totfypin,impact_energy,pos1,pos2
      CHARACTER tag*64

      WRITE(0,*) 'IDL DIVIMP DATA FILES'

      index = -1
      DO ik = 1, MAXNKS
        pos(ik) = ik
      ENDDO

c...  Calculate the impurity source that came back from EIRENE:
c       (from code in divoutput.f)
      totfypin= 0.0
      DO id = 1, wallpts
         in = wallpt(id,17)
         totfypin = totfypin + flxhw3(in) * wallpt(id,7)  ! per meter toroidally s-1
      ENDDO 
      IF (totfypin.EQ.0.0) totfypin = 1.0

c...  Dump impurity data:
      CALL inOpenInterface('idl.divimp_imp_density',ITF_WRITE)
      CALL inPutData(absfac  ,'DIV_IMPURITY_INFLUX','m-1 s-1')
      CALL inPutData(totfypin,'EIR_IMPURITY_INFLUX','m-1 s-1')
      CALL inPutData(cizsc,'IMP_INITIAL_IZ','N/A')
      CALL inPutData(nizs ,'IMP_MAX_IZ'    ,'N/A')
      CALL inPutData(SNGL(cion),'IMP_Z'         ,'N/A')
      CALL inPutData(crmi ,'IMP_A'         ,'N/A')
      CALL inPutData(irsep-1 ,'GRID_ISEP' ,'N/A')  ! Just passing these as a check when
      CALL inPutData(irtrap-2,'GRID_IPFZ' ,'N/A')  ! plotting with the grid geometry 
      CALL inPutData(eirtorfrac,'TOROIDAL_FRACTION' ,'N/A')  
      DO ir = 2, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        ike = nks(ir)
        IF (ir.LT.irsep) ike = ike - 1
        tube = ir - 1                      ! TUBE is set to the OSM fluid grid system, where
        IF (ir.GT.irwall) tube = tube - 2  ! the boundary rings are not present
        CALL inPutData(index(1:ike)   ,'INDEX' ,'N/A')                     
        CALL inPutData(pos  (1:ike)   ,'POS'   ,'N/A')                     
        CALL inPutData(tube (1:ike)   ,'TUBE'  ,'N/A')                     
        CALL inPutData(kss  (1:ike,ir),'S'     ,'m')                     
        CALL inPutData(kps  (1:ike,ir),'P'     ,'m')
        CALL inPutData(kvols(1:ike,ir),'VOLUME','m-3')
      ENDDO
      DO iz = 0, MIN(nizs,cion)
        WRITE(tag,'(A,I0.2)') 'IMP_DENS_',iz
        DO ir = 2, nrs
          IF (idring(ir).EQ.BOUNDARY) CYCLE
          ike = nks(ir)
          IF (ir.LT.irsep) ike = ike - 1
          CALL inPutData(sdlims(1:ike,ir,iz),TRIM(tag),'m-3 s-1')                     
        ENDDO
      ENDDO
      CALL inCloseInterface 

      WRITE(0,*) 'IDL DIVIMP DATA FILES 2'

      CALL inOpenInterface('idl.divimp_imp_ionisation',ITF_WRITE)
      CALL inPutData(absfac  ,'DIV_IMPURITY_INFLUX','m-1 s-1')
      CALL inPutData(totfypin,'EIR_IMPURITY_INFLUX','m-1 s-1')
      CALL inPutData(cizsc   ,'IMP_INITIAL_IZ'     ,'N/A')
      CALL inPutData(nizs    ,'IMP_MAX_IZ'         ,'N/A')
      CALL inPutData(cion    ,'IMP_Z'              ,'N/A')
      CALL inPutData(crmi    ,'IMP_A'              ,'N/A')
      CALL inPutData(irsep-1 ,'GRID_ISEP'          ,'N/A')  ! Just passing these as a check when
      CALL inPutData(irtrap-2,'GRID_IPFZ'          ,'N/A')  ! plotting with the grid geometry 
      DO ir = 2, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        ike = nks(ir)
        IF (ir.LT.irsep) ike = ike - 1
        tube = ir - 1                      ! TUBE is set to the OSM fluid grid system, where
        IF (ir.GT.irwall) tube = tube - 2  ! the boundary rings are not present
        CALL inPutData(index(1:ike)   ,'INDEX','N/A')                     
        CALL inPutData(pos  (1:ike)   ,'POS'  ,'N/A')                     
        CALL inPutData(tube (1:ike)   ,'TUBE' ,'N/A')                     
        CALL inPutData(kss  (1:ike,ir),'S'    ,'N/A')                     
      ENDDO
      DO iz = 0, MIN(nizs,cion)
        WRITE(tag,'(A,I0.2)') 'IMP_IONIZ_',iz
        DO ir = 2, nrs
          IF (idring(ir).EQ.BOUNDARY) CYCLE
          ike = nks(ir)
          IF (ir.LT.irsep) ike = ike - 1
          CALL inPutData(tizs(1:ike,ir,iz),TRIM(tag),'m-3 s-1')                     
        ENDDO
      ENDDO
      CALL inCloseInterface 
c
c     Target fluxes
c
c From mom.f in div6:
c CICABS - total ion flux to target
c CRAVAV - absolute velocity as the ion is lost
c CRTBS  - average background ele temperature at loss
c CRTABS - average background ion temperature at loss
c     RENEGY = 3.0 * REAL(IZ) * CTBS(IZ) / CICABS(IZ) +
c              5.22E-9 * CRMI * VEXIT * VEXIT +
c              2.0 * CRTABS(IZ) / CICABS(IZ)

c...  Just missing at the moment: velocity of the ion as it enters the sheath, 
c     which I'm leaving off for now...

      WRITE(0,*) 'IDL DIVIMP DATA FILES 3',nds,MIN(nizs,cion)


      CALL inOpenInterface('idl.divimp_flux_target',ITF_WRITE)
      CALL inPutData(absfac  ,'DIV_IMPURITY_INFLUX','m-1 s-1')
      CALL inPutData(totfypin,'EIR_IMPURITY_INFLUX','m-1 s-1')
      CALL inPutData(cizsc   ,'IMP_INITIAL_IZ'     ,'N/A')
      CALL inPutData(nizs    ,'IMP_MAX_IZ'         ,'N/A')
      CALL inPutData(SNGL(cion),'IMP_Z'              ,'N/A')
      CALL inPutData(crmi    ,'IMP_A'              ,'N/A')
      CALL inPutData(irsep-1 ,'GRID_ISEP'          ,'N/A')  ! Just passing these as a check when
      CALL inPutData(irtrap-2,'GRID_IPFZ'          ,'N/A')  ! plotting with the grid geometry 
      DO id = 1, nds
        ir = irds(id)
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        tube(1) = ir - 1                         ! TUBE is set to the OSM fluid grid system, where
        IF (ir.GT.irwall) tube(1) = tube(1) - 2  ! the boundary rings are not present
        target = IKHI
        IF (ikds(id).EQ.1) target = IKLO
        CALL inPutData(id             ,'INDEX'     ,'N/A')                     
        CALL inPutData(target         ,'TARGET'    ,'N/A')                     
        CALL inPutData(tube(1)        ,'TUBE'      ,'N/A')                     
        CALL inPutData(wallindex(id)  ,'INDEX_WALL','N/A')                     
      ENDDO
      DO id = 1, nds
        ir = irds(id)
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        CALL inPutData(rp(id)   ,'R_CEN' ,'m')
        CALL inPutData(zp(id)   ,'Z_CEN' ,'m') 
        CALL inPutData(dds(id)  ,'LENGTH','m') 
        CALL inPutData(kteds(id),'TE'    ,'eV')                     
        CALL inPutData(ktids(id),'TI'    ,'eV')                     
        DO iz = 1, MIN(nizs,cion)
          WRITE(tag,'(A,I0.2,A)') 'IMP_',iz,'_'
          impact_energy = 3.0 * kteds(id) * REAL(iz) +   ! Missing contribution from ion velocity at sheath entrance...
     .                    2.0 * ktids(id) 
c
c  NEROS(,1) - deposition, see NEUT.F, ION_PARALLEL_TRANSPORT.F, ION_TRANSPORT.F
c  NEROS(,2) - same as (,3), see NEUT.F
c  NEROS(,3) - erosion, see DIV.F
c  NEROS(,4) - net (,1) + (,3)
c  NEROS(,5) - 
c          FACT = 0.0
c          IF (TDEP.GT.0.0) FACT = TNEUT / TDEP
c          NEROS(ID,5) = FACT * NEROS(ID,1) + NEROS(ID,3)
c
c  From OUT000.F:
c        WRITE(ELABS(1),'(A,F8.4)')'    TOTAL DEPOSITION =',SUM(1)+SUM(6)
c        WRITE(ELABS(2),'(A,F8.4)')'    PRIMARY REMOVAL  =',SUM(2)+SUM(7)
c        WRITE(ELABS(3),'(A,F8.4)')'    TOTAL REMOVAL    =',SUM(3)+SUM(8)
c        WRITE(ELABS(4),'(A,2F7.4)')'    NET EROSION=',     SUM(4),SUM(9)
c        WRITE(ELABS(5),'(A,2F7.4)')'    NENNL      =',    SUM(5),SUM(10)
c
          CALL inPutData(deps(id,iz)  ,TRIM(tag)//'FLUX','m-2 s-1')                     
          CALL inPutData(impact_energy,TRIM(tag)//'E0'  ,'eV')                     
          CALL inPutData(0.0          ,TRIM(tag)//'VI'  ,'m s-1')                     
        ENDDO      
      ENDDO  
      CALL inCloseInterface 

c...  ASCII data file for Sophie and MatLab:
      fp = 99
      OPEN (UNIT=fp,FILE='mlb.erosion',ACCESS='SEQUENTIAL',
     .      STATUS='REPLACE')
      WRITE(fp,'(A)') '* DIVIMP data for MatLab - target erosion'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A)') '* Title: '//TRIM(title)
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A)') '* r      = rho in ribbon grid land'
      WRITE(fp,'(A)') '* z      = distance along the field line, s'
      WRITE(fp,'(A)') '* dist1  = position along the wall of the '//
     .                'start of the target segment'
      WRITE(fp,'(A)') '* dist2  = end of the target segment'
      WRITE(fp,'(A)') '*   So, centre of segment is (dis1+dist2)/2. '//
     .                ' The origin of this distance-along-the-wall'
      WRITE(fp,'(A)') '*   is at the upper left corner of the grid '//
     .                'where rho=0.0 and s=max(s), with the '
      WRITE(fp,'(A)') '*   distance then proceeding clockwise.'
      WRITE(fp,'(A)') '* theta   = angle between field line and '//
     .                'the target in degrees'
      WRITE(fp,'(A)') '* D+ flux = flux densiy on target relative '//
     .                'to the surface normal (not the field line)'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A)') '* The erosion/deposition units are just '//
     .                'funny DIVIMP units at the moment.  The'
      WRITE(fp,'(A)') '* ratio of erosion to D+ flux should give '//
     .                'a relative yield for each segment.'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A1,A5,3A6,4A10,3A12,A8,A10)')
     .  '*','wall','targ','cell','ring','r (m)','z (m)',
     .  'dist1 (m)','dist2 (m)','erosion','deposition','net',
     .  'theta','D+ flux'
      pos1 = 0.0
      pos2 = 0.0
      DO id = 1, wallpts
        in = NINT(wallpt(id,18))
        ik = ikds(MAX(1,in))
        ir = irds(MAX(1,in))
        pos2 = pos1 + wallpt(id,7)
        IF (in.NE.0.AND.ik.NE.0.AND.ir.NE.0) THEN
          IF (ik.EQ.0.OR.ir.EQ.0) CYCLE
          WRITE(fp,'(4I6,4F10.5,1P,3E12.2,0P,F8.2,1P,E10.2,0P)')
     .      id,in,
     .      ikds(in),irds(in),
     .      rp(in),zp(in),pos1,pos2,
     .      -neros(in,3),-neros(in,1),-neros(in,4),
     .      90.0-ACOS(costet(in))*180.0/PI,
     .      knds(in) * ABS(kvds(in)) * costet(in) * bratio(ik,ir) 
        ENDIF
        pos1 = pos2
      ENDDO
c
c
c
      WRITE(0,*) 'IDL DIVIMP DATA FILES 4'

      CALL inOpenInterface('idl.divimp_flux_wall',ITF_WRITE)
      CALL inPutData(absfac  ,'DIV_IMPURITY_INFLUX','m-1 s-1')
      CALL inPutData(totfypin,'EIR_IMPURITY_INFLUX','m-1 s-1')
      CALL inPutData(cizsc   ,'IMP_INITIAL_IZ'     ,'N/A')
      CALL inPutData(nizs    ,'IMP_MAX_IZ'         ,'N/A')
      CALL inPutData(SNGL(cion),'IMP_Z'              ,'N/A')
      CALL inPutData(crmi    ,'IMP_A'              ,'N/A')
      CALL inPutData(irsep-1 ,'GRID_ISEP'          ,'N/A')  ! Just passing these as a check when
      CALL inPutData(irtrap-2,'GRID_IPFZ'          ,'N/A')  ! plotting with the grid geometry 

c     FLUXHW - FLUX OF HYDROGEN (ATOMS AND MOLECULES) TO THE WALL
c     FLXHW2 - FLUX OF HYDROGEN (ATOMS AND IONS) TO THE WALL
c     FLXHW3 - FLUX OF IMPURITIES SPUTTERED FROM THE WALL (N/A)
c     FLXHW4 - FLUX OF IMPURITIES REDEPOSITED ONTO THE WALL (N/A)  --- *HACK* AVERAGE IMPURITY LAUNCH ENERGY
c     FLXHW5 - AVERAGE ENERGY OF ATOMS HITTING THE WALL (EV)
c     FLXHW6 - FLUX OF HYDROGEN ATOMS TO THE WALL
c     FLXHW7 - AVERAGE ENERGY OF MOLECULES HITTING THE WALL (eV)
c     FLXHW8 - EIRENE REPORTED HYDROGEN ION FLUXES TO THE WALL 

C     WALLPT (IND,1) = R
C     WALLPT (IND,2) = Z
C     WALLPT (IND,3) = WEIGHT FACTOR FOR ANTI-CLOCKWISE
C     WALLPT (IND,4) = WEIGHT FACTOR FOR CLOCKWISE
C     WALLPT (IND,5) = LENGTH OF 1/2 SEGMENT ANTI-CLOCKWISE
C     WALLPT (IND,6) = LENGTH OF 1/2 SEGMENT CLOCKWISE
C     WALLPT (IND,7) = TOTAL LENGTH OF LAUNCH SEGMENT
C     WALLPT (IND,8) = ANGLE FOR ANTI-CLOCKWISE LAUNCH
C     WALLPT (IND,9) = ANGLE FOR CLOCKWISE LAUNCH
C     WALLPT (IND,10) = NET PROBABILITY ANTI-CLOCKWISE
C     WALLPT (IND,11) = NET PROBABILITY CLOCKWISE
C     WALLPT (IND,12) = NET PROBABILITY FOR ENTIRE SEGMENT
C     WALLPT (IND,13) = FINAL PROBABILITY FOR SEGMENT
c
c     wallpt (ind,16) = TYPE OF WALL SEGMENT
c                       1 = Outer Target (JET) - inner for Xpt down
c                       4 = Inner Target (JET) - outer      "
c                       7 = Main Wall
c                       8 = Private Plasma Wall
c
c                       9 = Baffle Segment
c
c                       These are similar to the quantity in the JVESM
c                       array associated with the NIMBUS wall
c                       specification. The difference is that the
c                       Main Wall is split into Inner and Outer Divertor
c                       Wall as well as the Main (SOL) Wall - this
c                       is not done here.
c
c     WALLPT (ind,17) = INDEX into the NIMBUS flux data returned
c                       for each wall segment - ONLY if the NIMBUS
c                       wall option has been specified. NOTE: if
c                       the NIMBUS wall has been specified - it is
c                       still combined with the DIVIMP target polygon
c                       corners because rounding errors may result in
c                       small discrepancies between the coordinates.
c
c     WALLPT (IND,18) = Index of corresponding target segment if the wall
c                       segment is also a target segment.
c
c     WALLPT (IND,19) = Temperature of wall segment in Kelvin (K)
c
c     WALLPT (IND,20) = RSTART
c     WALLPT (IND,21) = ZSTART
c     WALLPT (IND,22) = REND
c     WALLPT (IND,23) = ZEND
c
c     wallpt (ind,24) = Used for additional indexing information - used
c                       as IK knot number for wall and trap wall option 7
c
c     wallpt (ind,25) = Value of reflection coefficient - if reflection
c                       for this segment is turned off the value here
c                       will be zero. If a positive value is specified
c                       then regular reflection occurs. If it is negative
c                       then a PTR (prompt thermal re-emission) type
c                       reflection is used. The value for this is
c                       set with the individual YMF's and is read from
c                       the CYMFS array.
c
c     wallpt (ind,26) = IK value of nearest plasma cell to wall segment
c     wallpt (ind,27) = IR value of nearest plasma cell to wall segment
c     wallpt (ind,28) = Minimum distance to outermost ring
c     wallpt (ind,29) = Plasma Te at wall segment - Temporary storage for RI
c     wallpt (ind,30) = Plasma Ti at wall segment - Temporary storage for ZI
c     wallpt (ind,31) = Plasma density at wall segment

      DO id = 1, wallpts
        CALL inPutData(id             ,'INDEX'      ,'N/A')                     
        CALL inPutData(wallpt(id,1)   ,'R_CEN'      ,'m')  
        CALL inPutData(wallpt(id,2)   ,'Z_CEN'      ,'m')                     
        CALL inPutData(wallpt(id,20)  ,'R_VERTEX1'  ,'m')                     
        CALL inPutData(wallpt(id,21)  ,'Z_VERTEX1'  ,'m')                     
        CALL inPutData(wallpt(id,22)  ,'R_VERTEX2'  ,'m')                     
        CALL inPutData(wallpt(id,23)  ,'Z_VERTEX2'  ,'m')                     
        CALL inPutData(wallpt(id,7)   ,'LENGTH'     ,'m')                     
        CALL inPutData(wallpt(id,19)  ,'TEMPERATURE','K')                     
        in = wallpt(id,17)
        CALL inPutData(in             ,'INDEX_PIN'    ,'N/A')                     
        CALL inPutData(flxhw6(in)     ,'ATOM_PAR_FLUX'  ,'m-2 s-1')                     
        CALL inPutData(flxhw5(in)     ,'ATOM_AVG_ENERGY','eV')                     
        CALL inPutData(fluxhw(in)-flxhw6(in),'MOL_PAR_FLUX'  ,'m-2 s-1')                     
        CALL inPutData(flxhw7(in)           ,'MOL_AVG_ENERGY','eV')                     
      ENDDO
      CALL inCloseInterface 

      WRITE(0,*) 'DONE'

c     FLUXHW - FLUX OF HYDROGEN (ATOMS AND MOLECULES) TO THE WALL
c     FLXHW2 - FLUX OF HYDROGEN (ATOMS AND IONS) TO THE WALL
c     FLXHW3 - FLUX OF IMPURITIES SPUTTERED FROM THE WALL (N/A)
c     FLXHW4 - FLUX OF IMPURITIES REDEPOSITED ONTO THE WALL (N/A)  --- *HACK* AVERAGE IMPURITY LAUNCH ENERGY
c     FLXHW5 - AVERAGE ENERGY OF ATOMS HITTING THE WALL (EV)
c     FLXHW6 - FLUX OF HYDROGEN ATOMS TO THE WALL
c     FLXHW7 - AVERAGE ENERGY OF MOLECULES HITTING THE WALL (eV)
c     FLXHW8 - EIRENE REPORTED HYDROGEN ION FLUXES TO THE WALL 


c...  ASCII data file for Sophie and MatLab:
      fp = 99
      OPEN (UNIT=fp,FILE='mlb.plasma',ACCESS='SEQUENTIAL',
     .      STATUS='REPLACE')
      WRITE(fp,'(A)') '* DIVIMP data for MatLab - background plasma'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A)') '* Title: '//TRIM(title)
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A)') '* r = rho in ribbon grid land'
      WRITE(fp,'(A)') '* z = distance along the field line, s'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A1,A5,A6,7A10)') 
     .  '*','cell','ring','r (m)','z (m)','n (m-3)','v (m s-1)',
     .  'Te (eV)','Ti (eV)'
      DO ir = 2, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        ike = nks(ir) 
        IF (ir.LT.irsep) ike = ike - 1
        DO ik = 1, ike
          WRITE(fp,'(2I6,2F10.5,1P,2E10.2,0P,2F10.2)')
     .      ik,ir,
     .      rs(ik,ir),zs(ik,ir),
     .      knbs(ik,ir),kvhs(ik,ir),ktebs(ik,ir),ktibs(ik,ir)
        ENDDO
      ENDDO

c...  ASCII data file for Sophie and MatLab:
      fp = 99
      OPEN (UNIT=fp,FILE='mlb.impurities',ACCESS='SEQUENTIAL',
     .      STATUS='REPLACE')
      WRITE(fp,'(A)') '* DIVIMP data for MatLab - impurity densities'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A)') '* Title: '//TRIM(title)
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A)') '* r = rho in ribbon grid land'
      WRITE(fp,'(A)') '* z = distance along the field line, s'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A)') '* The 0,1,2,... are the charge states of the '//
     .                'impurity ions.  The units are absolute densities'
      WRITE(fp,'(A)') '* in particle/m^3, I think... actually, not '//
     .                'completely sure for ribbon grids...'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A1,A5,A6,100A10)') 
     .  '*','cell','ring','r (m)','z (m)','0','1','2','3','4','5','6'
      DO ir = 2, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        ike = nks(ir) 
        IF (ir.LT.irsep) ike = ike - 1
        DO ik = 1, ike
          WRITE(fp,'(2I6,2F10.5,1P,100E10.2,0P)')
     .      ik,ir,
     .      rs(ik,ir),zs(ik,ir),
     .      (sdlims(ik,ir,iz)*absfac,iz=0,MAXIZS)
        ENDDO
      ENDDO


      RETURN
 99   STOP
      END
c
c ======================================================================
c
c
      SUBROUTINE GenerateOSMDataFiles
      USE mod_geometry
      IMPLICIT none

      INTEGER status

      CALL LoadGrid('osm.raw')
      CALL LoadObjects('osm_geometry.raw',status)

      CALL GenerateOutputFiles



      CALL SaveFluidGridGeometry       

      CALL osmClean
      CALL geoClean
 
      RETURN
 99   STOP
      END
c
c ======================================================================
c
c
      SUBROUTINE FluxSurfacesInTheSOL
      USE mod_interface
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slout'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      INTEGER ik,ir,id

      WRITE(6,'(A)') '* FLUX SURFACES FOR SPENCER'
      WRITE(6,'(A)') '* <rho (m)>'
      WRITE(6,'(A)') '* <number of points>'
      WRITE(6,'(A)') '* <index> <R (m)> <Z (m)>'
      WRITE(6,'(A)') '*'

      ir = irsep
      WRITE(6,*) rho(ir,IN14)
      WRITE(6,*) nks(ir)+1
      DO ik = 1, nks(ir)
        id = korpg(ik,ir)
        WRITE(6,'(I6,2F12.6)') ik,rvertp(1,id),zvertp(1,id)
      ENDDO
      WRITE(6,'(I6,2F12.6)') ik,rvertp(4,id),zvertp(4,id)

      DO ir = irsep+10, 82, 10
        IF (rho(ir,IN14).EQ.0.0) CYCLE
        WRITE(6,*) rho(ir,IN14)
        WRITE(6,*) nks(ir)+1
        DO ik = 1, nks(ir)
          id = korpg(ik,ir)
          WRITE(6,'(I6,2F12.6)') ik,rvertp(1,id),zvertp(1,id)
        ENDDO
        WRITE(6,'(I6,2F12.6)') ik,rvertp(4,id),zvertp(4,id)
      ENDDO

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c
      SUBROUTINE DumpDataToIDL
      USE mod_interface
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slout'
      INCLUDE 'comgra'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      REAL GetJsat,GetCs,CalcPressure

      INTEGER n,ik,ike,ir,i1,i2,id,itube,ipos
      REAL x1,x2,y1,y2,t,mach,pe,p,jsat
      REAL   , ALLOCATABLE :: x(:),y(:),v(:),s(:)
      CHARACTER   tag_x*11,tag_y*11,file*512
      CHARACTER*7 tag
      CHARACTER*2 target_tag(2)

      WRITE(0,*) 'IDL DUMP DATA FILES'

c...  Dump data for processing in IDL:
      file = 'osm.idl'
      WRITE(6,*) '999: Dumping OSM data to interface file'
      WRITE(6,*) '     FILE = >',TRIM(file),'<'
      CALL inOpenInterface(file,ITF_WRITE)
      n = 10  ! Resolution parameter for all cells... need something more refined...
      CALL inPutData(irsep,'grid_irsep','none')            
      CALL inPutData(nrs  ,'grid_nrs  ','none')            
      DO ir = 2, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        ike = nks(ir)
        IF (ir.LT.irsep) ike = nks(ir) - 1
        DO ik = 1, ike
          id = korpg(ik,ir)
c         Plasma:
          CALL inPutData(knbs    (ik,ir),'cell_ne'    ,'m-3')
          CALL inPutData(ktebs   (ik,ir),'cell_te'    ,'eV')
          CALL inPutData(pinalpha(ik,ir),'cell_dalpha','ph m-3 s-1')            
c         Geometry:
          CALL inPutData(kss(ik,ir)     ,'cell_s'     ,'m')
          CALL inPutData(rs(ik,ir)      ,'cell_x'     ,'m')            
          CALL inPutData(zs(ik,ir)      ,'cell_y'     ,'m')            
          CALL inPutData(ik             ,'cell_ik'    ,'none')            
          CALL inPutData(ir             ,'cell_ir'    ,'none')            
          CALL inPutData(nvertp(id)     ,'cell_npts'  ,'none')            
          DO i2 = 1, nvertp(id)
            WRITE(tag_x,'(A,I1)') 'cell_vtx_x',i2
            WRITE(tag_y,'(A,I1)') 'cell_vtx_y',i2
            CALL inPutData(rvertp(i2,id),tag_x,'m')            
            CALL inPutData(zvertp(i2,id),tag_y,'m')            
          ENDDO
        ENDDO
      ENDDO
c...  Setup interpolation data for 2D interpolation in IDL:
      n = 10  ! Resolution parameter for all cells... need something more refined...
      DO ir = 1, nrs
        ike = nks(ir)
        IF (ir.LT.irsep) ike = nks(ir) - 1
        ALLOCATE(x(ike*n))
        ALLOCATE(y(ike*n))
        ALLOCATE(s(ike*n))
        ALLOCATE(v(ike*n))
        i1 = 0
        DO ik = 1, ike
          IF (idring(ir).EQ.BOUNDARY) THEN
            IF     (ir.EQ.1.OR.ir.EQ.irtrap) THEN
              id = korpg(ikouts(ik,ir),irouts(ik,ir))
              x1 = rvertp(1,id)
              y1 = zvertp(1,id)
              x2 = rvertp(4,id)
              y2 = zvertp(4,id)
            ELSEIF (ir.EQ.irwall) THEN
              id = korpg(ikins(ik,ir),irins(ik,ir))
              x1 = rvertp(2,id)
              y1 = zvertp(2,id)
              x2 = rvertp(3,id)
              y2 = zvertp(3,id)
            ELSE
              CALL ER('DumpDataToIDL','Unrecognized boundary ring',*99)
            ENDIF
          ELSE
            id = korpg(ik,ir)
            x1 = 0.5 * (rvertp(1,id) + rvertp(2,id))
            y1 = 0.5 * (zvertp(1,id) + zvertp(2,id))
            x2 = 0.5 * (rvertp(3,id) + rvertp(4,id))
            y2 = 0.5 * (zvertp(3,id) + zvertp(4,id))
          ENDIF
          DO t = 0.0, 0.99, 1.0/REAL(n)
            i1 = i1 + 1
            s(i1) = (1.0 - t) * ksb(ik-1,ir) + t * ksb(ik,ir)
            x(i1) = (1.0 - t) * x1           + t * x2
            y(i1) = (1.0 - t) * y1           + t * y2
          ENDDO
        ENDDO
        CALL inPutData(x,'x','m')            
        CALL inPutData(y,'y','m')            
c...    Interpolate along the field line to improve spatial resolution:
        CALL Fitter(ike,kss(1,ir),pinalpha(1,ir),ike*n,s,v,'LINEAR')
        IF (ir.GE.irsep) THEN
          v(1    ) = 0.0
          v(ike*n) = 0.0
        ENDIF
        CALL inPutData(v,'inter_dalpha','ph m-3 s-1')            
        CALL Fitter(ike,kss(1,ir),knbs    (1,ir),ike*n,s,v,'LINEAR')
        IF (idring(ir).EQ.BOUNDARY.OR.ir.LT.irsep) v = 0.0
        IF (ir.GE.irsep) THEN
          v(1    ) = 0.0
          v(ike*n) = 0.0
        ENDIF
        CALL inPutData(v,'inter_ne','m-3')            
        CALL Fitter(ike,kss(1,ir),ktebs   (1,ir),ike*n,s,v,'LINEAR')
        IF (idring(ir).EQ.BOUNDARY.OR.ir.LT.irsep) v = 0.0
        IF (ir.GE.irsep) THEN
          v(1    ) = 0.0
          v(ike*n) = 0.0
        ENDIF
        CALL inPutData(v,'inter_te','eV')            
        DEALLOCATE(x)
        DEALLOCATE(y)
        DEALLOCATE(s)
        DEALLOCATE(v)
      ENDDO
      CALL inCloseInterface
c
c     ----------------------------------------------------------------------
c     Write out target data:
c
c     If changing anything here, need to change it in GenerateOutputFiles in
c     sol28_output.f as well, so that the OSM and OUT generated 
c     idl.fluid_targets files remain in sync.
c
      file = 'osm.idl.fluid_targets'
      target_tag(IKLO) = 'LO'
      target_tag(IKHI) = 'HI'
      WRITE(6,*) '999: Dumping OSM data to interface file - targets'
      WRITE(6,*) '     FILE = >',TRIM(file),'<'
      CALL inOpenInterface(file,ITF_WRITE)
      ir = irtrap
      IF (nopriv) ir = irsep
      DO WHILE(ir.NE.irwall-1)
        ir = ir + 1
        IF (ir.EQ.nrs+1) ir = irsep
        itube = ir
        IF (ir.LT.irwall) itube = itube - 1
        IF (ir.GT.irwall) itube = itube - 3
        CALL inPutData(itube        ,'TAR_TUBE','none')                    
        CALL inPutData(ir           ,'TAR_RING','none')                    
        CALL inPutData(psitarg(ir,2),'TAR_PSIN','none')                    
        CALL inPutData(rho(ir,CELL1),'TAR_RHO' ,'m') 
        DO ipos = IKLO, IKHI  
          IF (ipos.EQ.IKLO) THEN
            id = idds(ir,2)
          ELSE
            id = idds(ir,1)
          ENDIF
          tag = 'TAR_'//target_tag(ipos)//'_'
          IF (ipos.EQ.IKLO) THEN
            CALL inPutData(0.0       ,tag//'S','m')                    
            CALL inPutData(0.0       ,tag//'P','m') 
          ELSE
            CALL inPutData(ksmaxs(ir),tag//'S','m')                    
            CALL inPutData(-1.0      ,tag//'P','m')                    
          ENDIF
          pe = knds(id) * kteds(id)
          p  = CalcPressure(knds(id),kteds(id),ktids(id),kvds(id))
          jsat = GetJsat(kteds(id),ktids(id),knds(id),kvds(id))
          mach = kvds(id) / GetCs(kteds(id),ktids(id))
         CALL inPutData(-1         ,tag//'TARGET_INDEX','none')     
         CALL inPutData(-1         ,tag//'LOCATION'    ,'none')                    
          CALL inPutData(jsat      ,tag//'JSAT'  ,'Amps')                    
          CALL inPutData(knds(id)  ,tag//'NE'    ,'m-3')                    
          CALL inPutData(kvds(id)  ,tag//'V'     ,'m s-1')                    
          CALL inPutData(mach      ,tag//'MACHNO','none')                    
          CALL inPutData(   pe *ECH,tag//'PE'    ,'Pa')
          CALL inPutData((p-pe)*ECH,tag//'PI'    ,'Pa')                    
          CALL inPutData(kteds(id) ,tag//'TE'    ,'eV')                    
          CALL inPutData(ktids(id) ,tag//'TI'    ,'eV')                    
        ENDDO                      
      ENDDO
      CALL inCloseInterface

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c
      SUBROUTINE CoreRadialParticleTransport
      USE mod_interface
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slout'
      INCLUDE 'comgra'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER ik,ir1,ir2,ikmid1,ikmid2
      REAL    flux,deltax,deltan,gradn,sion,dperp,svol

      flux = 0.0

      DO ir1 = 1, irsep-1
        IF (idring(ir1).EQ.BOUNDARY) CYCLE

c...    Find outer midplane cells:
        ikmid1 = 0
        DO ik = 1, nks(ir1)-2
          IF (rs(ik,ir1).GT.r0.AND.
     .        ((zs(ik,ir1).GE.0.0.AND.zs(ik+1,ir1).LT.0.0).OR.
     .         (zs(ik,ir1).LT.0.0.AND.zs(ik+1,ir1).GE.0.0))) 
     .      ikmid1 = ik
        ENDDO
        IF (ikmid1.EQ.0)
     .    CALL ER('CoreRadialParticleTransport','Outer midplane '//
     .            'cell IKMID1 not identified',*99)

        ikmid2 = 0
        ir2 = ir1 + 1
        DO ik = 1, nks(ir2)-2
          IF (rs(ik,ir2).GT.r0.AND.
     .        ((zs(ik,ir2).GE.0.0.AND.zs(ik+1,ir2).LT.0.0).OR.
     .         (zs(ik,ir2).LT.0.0.AND.zs(ik+1,ir2).GE.0.0))) 
     .      ikmid2 = ik
        ENDDO
        IF (ikmid2.EQ.0)
     .    CALL ER('CoreRadialParticleTransport','Outer midplane '//
     .            'cell IKMID2 not identified',*99)

c...    Density gradient at outer boundary:
        gradn = 0.0
        SELECTCASE (0)
          CASE (0)  ! Transport at outer midplane
            deltax = rho(ir2,CELL1) - rho(ir1,CELL1)
            deltan = knbs(ikmid2,ir2) - knbs(ikmid1,ir1)
            gradn = deltan / deltax
c          CASE (1)  ! Volume averaged
          CASE DEFAULT
            CALL ER('CoreRadialParticleTransport','Unknown density '//
     .              'gradient option',*99)
        ENDSELECT
c        WRITE(0,*) '     :',ikmid1,ir1,ikmid2,ir2
c        WRITE(0,*) '     :',deltax,deltan
c        WRITE(0,*) 'GRADN:',gradn

c...    Area of outer boundary:

c...    Ionisation source on the ring:
        SELECTCASE (1)
          CASE (0)  ! Outer midplane
            deltax = rho(ir1,OUT23) - rho(ir1,IN14)
            sion = pinion(ikmid1,ir1) * deltax
          CASE (1)  ! Volume averaged source, transport at outer midplane
            deltax = rho(ir1,OUT23) - rho(ir1,IN14)
            sion = 0.0
            svol = 0.0
            DO ik = 1, nks(ir1)-1
              sion = sion + pinion(ik,ir1) * kvols(ik,ir1)
              svol = svol + kvols (ik,ir1)
            ENDDO
            sion = sion / svol * deltax
          CASE DEFAULT
            CALL ER('CoreRadialParticleTransport','Unknown density '//
     .              'gradient option',*99)
        ENDSELECT

c        WRITE(0,*) 'SION :',deltax,pinion(ikmid1,ir1),sion


        flux = flux + sion

        dperp = -1.0 * flux / gradn

c        WRITE(0,*) 'FLUX : ',flux
c        WRITE(0,*) 'DPERP: ',dperp

        WRITE(0,'(A,3F10.6,1P,5E10.2,0P,F10.3) ')
     .    'DPERP:',psitarg(ir1,2),rho(ir1,CELL1),rs(ikmid1,ir1),
     .     knbs(ikmid1,ir1),sion/deltax,sion,flux,gradn,dperp

      ENDDO


      RETURN
 99   STOP
      END
c
c ======================================================================
c taken from tau.d6a
c
       subroutine calc_wallprad(nizs)
       implicit none
c
       integer nizs 
c
       include 'params'
       include 'cgeom'
       include 'comtor'
       include 'dynam3'
       include 'printopt'
c
c      CALC_WALLPRAD:
c
c      This routine calculates the total radiated power onto each wall segment.
c      It also calculates the total radiated H and IMP power. 
c      In order to calculate the amount on each wall segment it has to 
c      allocate the amount from each cell that would strike each wall elememt.
c      This routine uses the data in HPOWLS and POWLS to make these calculations.
C      The POWLS data is multiplied by the ABSFAC quantity to get absolute scalings
c      unless ABSFAC is zero. 
c
       integer in,ir,ik,iz,it
       real*8 imp_prad,h_prad,tot_prad
       real*8 rcent,zcent
c
       integer n_int,max_int,backward,forward,ivert
       parameter(max_int=10,backward=-1,forward=1)
       real*8 r_int(max_int),z_int(max_int)
       integer n_intersections(maxpts+1),next_vert
       external next_vert 
c
       real*8 datan2c
       real*8 wall_angle(maxpts+1),net_angle(maxpts),result_angle
       external datan2c
c
       integer nstart,nend,ntest,segtype
       real*8 rtest,ztest
c
c      Initialization 
c
       REAL wallprad2(maxpts+7,3)

       wallprad2 = 0.0

c      call rzero(wallprad,(maxpts+7)*3)
c
c      For every cell on the grid - this code must loop through every element of 
c      the wall.  
c       
       do ir = 1,nrs
         do ik = 1,nks(ir)
c
c           Sum up impurity and hydrogen radiation contributions
c             
c           Impurity
c
            imp_prad = 0.0
c 
            do iz = 0,nizs
c
               imp_prad = imp_prad + powls(ik,ir,iz)
c
            end do
c
            if (absfac.gt.0.0) imp_prad = imp_prad * absfac
c
c           Hydrogen
c
            h_prad = 0.0
c 
            do iz = 0,1
c
               h_prad = h_prad + hpowls(ik,ir,iz)
c
            end do
c
c           Scale for cell area and major radius - converting to watts
c
            imp_prad = imp_prad * kareas(ik,ir) * 2.0 * PI * rs(ik,ir) 
            h_prad   = h_prad   * kareas(ik,ir) * 2.0 * PI * rs(ik,ir) 
            tot_prad = h_prad + imp_prad
c
c           Don't perform geometry calculations unless there is a radiation contribution 
c
            if (tot_prad.eq.0.0) cycle
c
c           Record totals of radiated source in highest array element - these
c           will be compared to totals integrated over all wall elements
c           and should be the same.
c
            wallprad2(maxpts+7,1) = wallprad2(maxpts+7,1) + imp_prad
            wallprad2(maxpts+7,2) = wallprad2(maxpts+7,2) + h_prad
            wallprad2(maxpts+7,3) = wallprad2(maxpts+7,3) + tot_prad
c
c           Calculate the geometry setup information - can each wall vertex be 
c           seen from the current cell and what are the angles to all of the 
c           vertices.
c 
            rcent = rs(ik,ir)
            zcent = zs(ik,ir)  
c
            do in = 1,wallpts
c
c              Calculate initial intersections and whether LOS to each
c              vertex is clear. Last point and first point of the wall
c              should be identical
c
               n_intersections(in) = 0
c
               rtest = wallpt(in,20)
               ztest = wallpt(in,21)
c
               call calc_wall_intersections(ntest,max_int,r_int,z_int,
     >                      rcent,zcent,rtest,ztest,.true.)
c            
               n_intersections(in) = ntest 
c
c              Calculate wall_angle array
c
c              Get angles to the start and end of the wall segment
c 
               wall_angle(in) = datan2c(ztest-zcent,rtest-rcent)
c
c              Convert angles to 0 to 2 PI range  from -PI to PI
c
               if (wall_angle(in).lt.0.0)
     >             wall_angle(in)=wall_angle(in)+2.0*PI
c
            end do
c
c           Assign values to complete the wall
c
            n_intersections(wallpts+1) = n_intersections(1)
            wall_angle(wallpts+1) = wall_angle(1)
c
c           Calculate the net angle for each segment of the wall
c
            do in = 1,wallpts
c
               net_angle(in) = wall_angle(in+1)- wall_angle(in)
c
               if (net_angle(in).gt.PI) 
     >             net_angle(in) = net_angle(in) - 2.0 * PI
               if (net_angle(in).lt.-PI) 
     >             net_angle(in) = net_angle(in) + 2.0 * PI 
c
            end do 
c
c  
c           Basic LOS rules:
c           
c           1) Wall is listed in a basically clockwise order so assuming that 
c              all angles are positive in the range 0.0 to 2 PI then the angle
c              to the leading point of a wall segment will always be less 
c              than the angle to the second point for a forward oriented wall 
c              segment. 
c           2) All wall elements that have some exposure to radiation from the
c              cell being examined will be forward going assuming all radiation 
c              is from within the closed vessel. 
c           3) If the LOS from the cell to the ends of the wall segment cross the
c              the wall at any point - then the view to that wall segment is 
c              assumed to be blocked.
c           4) If the LOS from the cell to the ends of the wall segment does
c              not cross the wall then the LOS is considered unblocked and the
c              entire segment is visible to the source. 
c           5) If either end of a backward going wall segment is blocked then 
c              none of that line segment can be seen from the cell. 
c           6) If the wall segment is forward going but one vertex is blocked then
c              the wall element is partially obscured and the proportion of the
c              segment exposed to the source cell is approximately calculated.  
c           7) Only clockwise oriented segments can receive radiation (net_angle<0)
c           8) The total exposure of a partially blocked segment is limited by
c              the visible vertex and either the next or last visible vertex along 
c              the wall. 
c 
c

            do in = 1,wallpts
c
               if (net_angle(in).lt.0.0) then 
c
c                 If net_angle is less than zero then the wall element is 
c                 forward going or clockwise relative to the observation point.
c             
c                 Check number of intersections for the wall element vertices  
c
c                 Totals only need updating for unobstructed and partially 
c                 obstructed views.
c
c             
c                 Wall segment is unobstructed 
c
                  if (n_intersections(in).eq.0.and.
     >                n_intersections(in+1).eq.0) then 
c
                     result_angle = -net_angle(in) 
c
c
c                 Clockwise wall segment with one blocked vertex.
c                 Need to calculate the required angle
c                 Need angle of last unobstructed vertex
c
                  elseif (n_intersections(in).gt.0.and.
     >                   n_intersections(in+1).eq.0) then 
c      	    
                     ivert = next_vert(n_intersections,
     >                                 wallpts+1,in,backward)
c      	    
                     result_angle = wall_angle(in+1)- wall_angle(ivert)
c      	    
                     if (result_angle.gt.PI) 
     >                   result_angle = result_angle - 2.0 * PI
                     if (result_angle.lt.-PI) 
     >                   result_angle = result_angle + 2.0 * PI 

                     result_angle = -result_angle
c
c
c                 Need angle of next unobstructed vertex
c
                  elseif (n_intersections(in).eq.0.and.
     >                   n_intersections(in+1).gt.0) then
       	  
                     ivert = next_vert(n_intersections,
     >                                 wallpts+1,in+1,forward)
c      	  
                     result_angle = wall_angle(ivert)- wall_angle(in)
c      	  
                     if (result_angle.gt.PI) 
     >                   result_angle = result_angle - 2.0 * PI
                     if (result_angle.lt.-PI) 
     >                   result_angle = result_angle + 2.0 * PI 

                     result_angle = -result_angle
c
c
c                 LOS obstructed
c
                  elseif (n_intersections(in).gt.0.and.
     >                    n_intersections(in+1).gt.0) then
c
                      result_angle = 0.0 
c
                  endif
c
c                 Calculate contributions to this wall segment
c
                  wallprad2(in,1) = wallprad2(in,1) 
     >                        + imp_prad * result_angle/(2.0*PI)
                  wallprad2(in,2) = wallprad2(in,2) 
     >                        + h_prad * result_angle/(2.0*PI)
                  wallprad2(in,3) = wallprad2(in,3) 
     >                        + tot_prad * result_angle/(2.0*PI)
               end if

            end do
c
        end do
c
      end do  
c
c     Print outs of results
c
      do in = 1,wallpts
c
c       The following grand totals are calculated on a complete torus
c       basis factoring in the major radius. 
c
c
c       Check the tag for the wall element to see which type of 
c       segment it is:
c
c       1 = Target 1 (Outer target for X-point up grids)
c       4 = Target 2 (Inner target for X-point up grids)
c       7 = (and others 2,3) Main Vessel wall
c       8 = PFZ wall
c       9,10 = baffle segments - added to PFZ only applies 
c                                for some JET grids 
c
        segtype = wallpt(in,16)
c
c       First Target (Outer for X-point up - inner for down)
c  
        if (segtype.eq.1) then          
         wallprad2(maxpts+1,1) = wallprad2(maxpts+1,1) + wallprad2(in,1)
         wallprad2(maxpts+1,2) = wallprad2(maxpts+1,2) + wallprad2(in,2)
         wallprad2(maxpts+1,3) = wallprad2(maxpts+1,3) + wallprad2(in,3)
c
c       Second Target (Inner for X-point up - Outer for down)
c  
        elseif (segtype.eq.4) then          
         wallprad2(maxpts+2,1) = wallprad2(maxpts+2,1) + wallprad2(in,1)
         wallprad2(maxpts+2,2) = wallprad2(maxpts+2,2) + wallprad2(in,2)
         wallprad2(maxpts+2,3) = wallprad2(maxpts+2,3) + wallprad2(in,3)
c
c       Main Vessel wall elements
c
        elseif (segtype.eq.7.or.segtype.eq.2.or.segtype.eq.3) then
         wallprad2(maxpts+3,1) = wallprad2(maxpts+3,1) + wallprad2(in,1)
         wallprad2(maxpts+3,2) = wallprad2(maxpts+3,2) + wallprad2(in,2)
         wallprad2(maxpts+3,3) = wallprad2(maxpts+3,3) + wallprad2(in,3)
c
c       Private Flux Zone wall elements   
c
        elseif (segtype.eq.8.or.segtype.eq.9.or.segtype.eq.10) then
         wallprad2(maxpts+4,1) = wallprad2(maxpts+4,1) + wallprad2(in,1)
         wallprad2(maxpts+4,2) = wallprad2(maxpts+4,2) + wallprad2(in,2)
         wallprad2(maxpts+4,3) = wallprad2(maxpts+4,3) + wallprad2(in,3)
        else
         wallprad2(maxpts+5,1) = wallprad2(maxpts+5,1) + wallprad2(in,1)
         wallprad2(maxpts+5,2) = wallprad2(maxpts+5,2) + wallprad2(in,2)
         wallprad2(maxpts+5,3) = wallprad2(maxpts+5,3) + wallprad2(in,3)
        endif
c
      end do 
c
c     Sum up grand totals
c
      do in = 1,5
         do it = 1,3
            wallprad2(maxpts+6,it) = wallprad2(maxpts+6,it) 
     >                            + wallprad2(maxpts+in,it)
         end do
      end do
c
c     Convert all of the wall element data to W/m2
      do in = 1,wallpts
         do it = 1,3
            if (wallpt(in,7).gt.0.0.and.wallpt(in,1).gt.0.0) then
               wallprad2(in,it) = wallprad2(in,it)
     >                           /(wallpt(in,7)*2.0*PI*wallpt(in,1))
            endif
         end do
      end do
c
c
c     Print summary to unit 6 
      write(6,*) 
      write(6,'(a)') 'WALL Radiation Flux Summary [W/m2]:'
      write(6,'(2(a4,2X),4A12)') 'ID','TYPE','IMPURITY','HYDROGEN',
     .                           'TOTAL'
c     Wall elements
      do in = 1,wallpts
         write(6,'(2(i4,2x),4(1x,g12.5))') in,int(wallpt(in,16)),
     >             (wallprad2(in,it),it=1,3)
      end do
c
c
      write(6,*) 
      write(6,*) 'Regional Radiation Totals in [W]'
      write(6,*) 
      write(6,80) 'TOT  '//inner,(wallprad2(maxpts+1,it),it=1,3)
      write(6,80) 'TOT  '//outer,(wallprad2(maxpts+2,it),it=1,3)
      write(6,80) 'TOT   MAIN'  ,(wallprad2(maxpts+3,it),it=1,3)
      write(6,80) 'TOT    PFZ'  ,(wallprad2(maxpts+4,it),it=1,3)
      write(6,80) 'TOT   MISC'  ,(wallprad2(maxpts+5,it),it=1,3)
      write(6,*) 
      write(6,80) 'TOTAL  SEG'  ,(wallprad2(maxpts+6,it),it=1,3)
      write(6,80) 'TOTAL  SRC'  ,(wallprad2(maxpts+7,it),it=1,3)
      write(6,*) 

 80   FORMAT(a10,1P,4(1x,e12.5),0P)

      return
      end
c
c
c
      integer function next_vert(n_intersections,maxp,startin,step)
      implicit none
      integer startin,maxp,step
      integer n_intersections(maxp)
c
c     NEXT_VERT: This routine finds the next entry in n_intersections
c                which is equal to zero - moving in the direction
c                defined by step and starting at startin. 
c
c
      integer current
c
      current = startin
c  
      do while (n_intersections(current).ne.0)
c
         current = current + step
c          
         if (current.gt.maxp) current = 1
         if (current.lt.1) current = maxp
c
      end do
c 
      next_vert = current
c
      return
      end
c
c 
c
      subroutine calc_wall_intersections(n_int,max_int,r_int,z_int,
     >                       rstart,zstart,rend,zend,ignore_end)
      implicit none
      integer n_int,max_int
      real*8 r_int(max_int),z_int(max_int)
      real*8 rstart,zstart,rend,zend
      logical ignore_end
c
      include 'params'
      include 'comtor'
c
c     CALC_WALL_INTERSECTIONS:
c
c     This routine loops through the entire wall an reports all
c     of the intersections of the given line (rstart,zstart) (rend,zend) 
c     with the wall. The ignore_end flag instructs the code to not
c     calculate any intersections for the end-point since it is known
c     to be on the wall. 
c
      integer in
      real*8 rnew,znew,tnew,tnorm
      logical sect
c
      real eps
      parameter (eps=1.0e-6)
c
      n_int = 0
      call rzero(r_int,max_int)
      call rzero(z_int,max_int)
c
      do in = 1,wallpts
c
c         Initialize for call - most not needed
c
          rnew = 0.0
          znew = 0.0         
          tnew = 0.0
          tnorm= 0.0
          sect = .false.
c
          CALL INTCALCDP(Rend,Zend,Rstart,Zstart,
     >                 dble(WALLPT(IN,1)),dble(WALLPT(IN,2)),
     >                 dble(WALLPT(IN,8)),dble(WALLPT(IN,9)),
     >                 dble(WALLPT(IN,5)),dble(WALLPT(IN,6)),
     >                 RNEW,ZNEW,TNEW,tnorm,
     >                 SECT,nrfopt)
c
c         Record any intersections found
c       
          if (sect) then 
c
c            If not (endpoint ignored and this is the endpoint)
c
             if (.not.(ignore_end.and.
     >             (abs(rnew-rend).lt.eps.and.
     >              abs(znew-zend).lt.eps))) then 
c
c               Record point
c
                if (n_int.lt.max_int) then 
                   n_int = n_int + 1
                   r_int(n_int) = rnew
                   z_int(n_int) = znew   
                endif 
c
             endif
c
          endif  
c
      end do
c  
      return
      end
c
c ======================================================================
c
c
      SUBROUTINE OutputDivertorProfiles
      USE mod_eirene04
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slout'
      INCLUDE 'comgra'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'


      INTEGER ik,ir,fp,ix,iy,nx,ny,ngrid,i1,n,i2
      REAL*8  dx,dy,x1,x2,y1,y2,x(100000),y(100000),v(100000),
     .        xgrid(10000),ygrid(10000),vgrid(10000,5),
     .        x9(10),y9(10),v9(10)


      WRITE(0,*) '8:OUTPUT DIVERTOR PROFILES'

      ngrid = 0

c...  
      fp = 99
      OPEN (UNIT=fp,FILE='div.random',ACCESS='SEQUENTIAL',
     .      STATUS='REPLACE')
      WRITE(fp,'(A)') '* OSM-EIRENE profiles for 990429019 at 950 ms '//
     .                'on irregular X,Y mesh'
      WRITE(fp,'(A,5(A15))') '*','x','y','n_e','T_e','Ly_alpha'
      WRITE(fp,'(A,5(A15))') '*','(m)','(m)','(m-3)','(eV)',
     .                       '(ph m-3 s-1)'
      DO ir = 1, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        DO ik = 1, nks(ir)      
          IF (zs(ik,ir).GT.-0.2) CYCLE

          ngrid = ngrid + 1
          xgrid(ngrid) = rs(ik,ir)
          ygrid(ngrid) = zs(ik,ir)
          vgrid(ngrid,1) = knbs    (ik,ir)
          vgrid(ngrid,2) = ktebs   (ik,ir)
          vgrid(ngrid,3) = pinalpha(ik,ir)

          WRITE(fp,'(1X,2F15.6,1P,E15.4,0P,F15.2,1P,E15.4,0P)')
     .      xgrid(ngrid),ygrid(ngrid),vgrid(ngrid,1:3)
c     .      rs(ik,ir),zs(ik,ir),knbs(ik,ir),ktebs(ik,ir),pinalpha(ik,ir)

        ENDDO
      ENDDO
      CLOSE(fp)



      IF (.FALSE.) THEN
        DO i1 = 1, 10
          DO i2 = 1, 10
            x(i2+(i1-1)*10) = REAL(i2)
            y(i2+(i1-1)*10) = SQRT(REAL(i1))
            v(i2+(i1-1)*10) = REAL(i2)
          ENDDO
          x9(i1) = REAL(i1)+0.5
          y9(i1) = 5.0
        ENDDO
        WRITE(0,*) 'x:',x(1:10)
        WRITE(0,*) 'y:',y(1:10)
        WRITE(0,*) 'v:',v(1:10)
        CALL Interpolate(100,x,y,v,9,x9,y9,v9,0)
        WRITE(0,*) 'x1:',x9
        WRITE(0,*) 'y1:',y9
        WRITE(0,*) 'v1:',v9

        RETURN
      ENDIF
   




c...  
      x1 =  0.42
      x2 =  0.72
      y1 = -0.63
      y2 = -0.23
      nx = 30
      ny = 40
      n = nx*ny
      dx = (x2 - x1) / REAL(nx)
      dy = (y2 - y1) / REAL(ny)
      DO iy = 1, ny
        DO ix = 1, nx
          x(ix+(iy-1)*nx) = x1 + (0.5+REAL(ix-1))*dx
          y(ix+(iy-1)*nx) = y1 + (0.5+REAL(iy-1))*dy
        ENDDO
      ENDDO

c      DO iy = 1, ny
c        DO ix = 1, nx
c          WRITE(0,*) ix,x(ix+(iy-1)*ny),y(ix+(iy-1)*ny)
c        ENDDO
c      ENDDO

      fp = 99
      OPEN (UNIT=fp,FILE='div.regular',ACCESS='SEQUENTIAL',
     .      STATUS='REPLACE')
      WRITE(fp,'(A)') '* OSM-EIRENE profiles for 990429019 at 950 ms '//
     .                'on regular X,Y mesh'
      WRITE(fp,'(A,5(A15))') '*','x','y','n_e','T_e','Ly_alpha'
      WRITE(fp,'(A,5(A15))') '*','(m)','(m)','(m-3)','(eV)',
     .                       '(ph m-3 s-1)'

c     CALL Interpolate(n,x,y,v,1,x(1),y(1),v(1),0)
c      CALL Interpolate(10,xgrid(1:10),ygrid(1:10),vgrid(1:10,1),
c     .                 10,x    (1:10),y    (1:10),v    (1:10),0)
c      CALL Interpolate(ngrid,xgrid,ygrid,vgrid(1,1),n,x,y,v,0)

      DO i1 = 1, n
        WRITE(fp,'(1X,2F15.6,1P,E15.4,0P,F15.2,1P,E15.4,0P)')
     .    x(i1),y(i1),0.0,0.0,0.0
      ENDDO
 
      CLOSE(fp)


      RETURN
99    STOP
      END
c
c ======================================================================
c
c
      SUBROUTINE CoreProfileAnalysis(nizs,cizsc,crmi,cion,absfac)
      USE mod_eirene04
      USE mod_interface
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slout'
      INCLUDE 'comgra'
      INCLUDE 'cgeom'
      INCLUDE 'walls_com'
      INCLUDE 'dynam2'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER, INTENT(IN) :: nizs,cizsc,cion
      REAL   , INTENT(IN) :: crmi,absfac

      REAL CalcPressure, GetCs

      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost

      INTEGER ik,ir,npro,ikmid,i1,i2,fp,nre,itri,id,ncol,iz,iz1(200),
     .        ikmid1(200),ring(200),in
      CHARACTER dummy*256
      REAL    midpro(200,10),avolpro(200,0:100),psin(200),r(200),pmax,
     .        tvolp(200,0:100),deltar(200),dpsin(200),frac,p(200),
     .        dp(200),t_D(MAXNKS,MAXNRS),t_D2(MAXNKS,MAXNRS),
     .        r1,r2,z1,z2,rho1(200),scale,impurity_influx,rsmax,rsmin,
     .        pressure,machno,eirene_influx,
     .        totfpin,totfypin,tmpy,totzfpin,tothpin,totfapin,totmhpin

      REAL, ALLOCATABLE :: tdata(:)

c...  Needed:
c     -volume averaged + outer midplane profiles + poloidal distributions v theta v normalized (origin at midplane), radial averages
c     -ionisation,D density/temperature,D2 density/temperature,Dalpha
c     -also ne,Te,Ti
c     -eventually: vb, toroidal torque

c      CALL LoadTriangles
c      ALLOCATE(tdata (ntri))
c      ALLOCATE(tdata1(ntri))


      fp = 98
      OPEN (UNIT=fp,FILE='core_analysis.dat',ACCESS='SEQUENTIAL',
     .      STATUS='REPLACE')


c     ------------------------------------------------------------------
c     IMPURITIES
c     ------------------------------------------------------------------

c...  Calculate the impurity source that came back from EIRENE:
c       (from code in divoutput.f)
      totfpin = 0.0
      totfypin= 0.0
      totzfpin= 0.0
      tothpin = 0.0
      totmhpin= 0.0
      totfapin= 0.0
      DO id = 1, wallpts
         in = wallpt(id,17)
         totfpin  = totfpin  + flxhw2(in) * wallpt(id,7)
         totfypin = totfypin + flxhw3(in) * wallpt(id,7)
         totzfpin = totzfpin + flxhw4(in) * wallpt(id,7)
         totfapin = totfapin + flxhw6(in) * wallpt(id,7)
         tothpin  = tothpin  + flxhw6(in) * flxhw5(in) * wallpt(id,7)
         totmhpin = totmhpin + (fluxhw(in)-flxhw6(in))* kboltz
     .                        * wallpt(id,19) * wallpt(id,7)
      ENDDO 
      impurity_influx = absfac  ! per meter toroidally s-1  changed on 02/04/2010 -SL
      eirene_influx = totfypin  ! per meter toroidally s-1
      IF (totfypin.EQ.0.0) totfypin = 1.0

c...  Outer midplane profiles:
      npro = 0
      midpro = 0.0  
      scale = 1.0
      DO ir = 2, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        IF (ir.EQ.irwall.OR.ir.EQ.irtrap) CYCLE   ! Not sure why, but seems necessary...
        ikmid = 0
        DO ik = 1, nks(ir)-2
          IF (rs(ik,ir).GT.r0.AND.
c     .        ((zs(ik,ir).GE.0.64.AND.zs(ik+1,ir).LT.0.64).OR.
c     .         (zs(ik,ir).LT.0.64.AND.zs(ik+1,ir).GE.0.64))) THEN
     .        ((zs(ik,ir).GE.z0.AND.zs(ik+1,ir).LT.z0).OR.
     .         (zs(ik,ir).LT.z0.AND.zs(ik+1,ir).GE.z0))) THEN
c     .        ((zs(ik,ir).GE.0.0.AND.zs(ik+1,ir).LT.0.0).OR.     ! changed 05/01/2010
c     .         (zs(ik,ir).LT.0.0.AND.zs(ik+1,ir).GE.0.0))) THEN
            ikmid = ik
          ENDIF
        ENDDO
        IF (ikmid.EQ.0) CYCLE  ! Flux tube does not cross the y-axis on the LFS
        rsmax = 0.0
        DO ik = 1, nks(ir)-2
          IF (rs(ik,ir).GT.rsmax.AND.ABS(zs(ik,ir)).LT.0.3*r0) THEN
            ikmid = ik
            rsmax = rs(ik,ir)
          ENDIF
        ENDDO
        npro = npro + 1
        ring(npro) = ir
        rho1(npro) = rho(ir,CELL1) 
        r(npro) = rs(ikmid,ir)
        psin(npro) = psitarg(ir,2)
        ncol = 0
        DO iz = 5, 25, 5
          ncol = ncol + 1
          IF (ir.EQ.2) iz1(ncol) = iz
          IF (iz.EQ.5) ikmid1(npro) = ikmid
          midpro(npro,ncol) = sdlims(ikmid,ir,iz) * scale
c          WRITE(0,*) ir,iz,sdlims(:,ir,iz)
        ENDDO
      ENDDO
      WRITE(fp,*) '* Outer midplane profiles (max of bulge):'
      WRITE(fp,'(A4,8X,4A9,2X,A9,A10,A6,A9,2A9,10I10)') 
     .  '*','r','rho','psin','L','ne','vb','M','P',
     .  'Te','Ti',(iz1(i1),i1=1,ncol)
      WRITE(fp,*) npro
      DO i1 = 1, npro
        machno = kvhs(ikmid1(i1),ring(i1)) / qt /
     .           GetCs(ktebs(ikmid1(i1),ring(i1)),
     .                 ktibs(ikmid1(i1),ring(i1)))
        pressure = CalcPressure(knbs (ikmid1(i1),ring(i1)),
     .                          ktebs(ikmid1(i1),ring(i1)),
     .                          ktibs(ikmid1(i1),ring(i1)),
     .                          kvhs (ikmid1(i1),ring(i1))/qt)* ECH
        WRITE(fp,'(3I4,3F9.5,F9.2,2X,'//
     .           '1P,E9.2,E10.2,0P,F6.2,1P,E9.2,0P'//
     .           '2F9.2,1P,9E10.2,2X,I4)') 
     .    i1,ikmid1(i1),ring(i1),
     .    r(i1),rho1(i1),psin(i1),ksmaxs(ring(i1)),
     .    knbs(ikmid1(i1),ring(i1)),
     .    kvhs(ikmid1(i1),ring(i1)) / qt,
     .    machno,pressure,
     .    ktebs(ikmid1(i1),ring(i1)),ktibs(ikmid1(i1),ring(i1)),
     .    (midpro(i1,i2),i2=1,ncol)
      ENDDO

      WRITE(0,*) 'IDL CORE DATA FILES'

      CALL inOpenInterface('osm.idl.midplane',ITF_WRITE)
c      CALL inPutData(ring  (     1:npro ),'MID_IMPURITY_SOURCE','m-2 s-1')
      CALL inPutData(ring  (     1:npro ),'MID_RING','none')
      CALL inPutData(r     (     1:npro ),'MID_R'   ,'m')
      CALL inPutData(rho1  (     1:npro ),'MID_RHO' ,'m')
      CALL inPutData(psin  (     1:npro ),'MID_PSIN','none')
      CALL inPutData(ksmaxs(ring(1:npro)),'MID_L'   ,'m')
      DO i1 = 1, npro
        machno = kvhs(ikmid1(i1),ring(i1)) / qt / 
     .           GetCs(ktebs(ikmid1(i1),ring(i1)),
     .                 ktibs(ikmid1(i1),ring(i1)))
        pressure = CalcPressure(knbs (ikmid1(i1),ring(i1)),
     .                          ktebs(ikmid1(i1),ring(i1)),
     .                          ktibs(ikmid1(i1),ring(i1)),
     .                          kvhs (ikmid1(i1),ring(i1))/qt)* ECH
        CALL inPutData(knbs (ikmid1(i1),ring(i1))   ,'MID_NE','m-3')
        CALL inPutData(kvhs (ikmid1(i1),ring(i1))/qt,'MID_VB','m s-1')
        CALL inPutData(machno                       ,'MID_M' ,'none')
        CALL inPutData(pressure                     ,'MID_P' ,'Pa')
        CALL inPutData(ktebs(ikmid1(i1),ring(i1))   ,'MID_TE','eV')
        CALL inPutData(ktibs(ikmid1(i1),ring(i1))   ,'MID_TI','eV')
      ENDDO
      CALL inCloseInterface

c...  Volume averaged radial profiles in the core:
      WRITE(0 ,*) 'NIZS             =',nizs
      WRITE(0 ,*) 'IMPURITY_INFLUX  =',impurity_influx
      WRITE(0 ,*) 'TOROIDAL_FRACTION=',1.0
      WRITE(fp,*) 'NIZS=             ',nizs
      WRITE(fp,*) 'IMPURITY_INFLUX=  ',impurity_influx
      WRITE(fp,*) 'TOROIDAL_FRACTION=',1.0
      impurity_influx = impurity_influx * 1.0
      avolpro = 0.0  ! Full volume averaged radial profiles
      tvolp = 0.0  ! Totals
      npro = 0
      DO ir = 2, irsep-1
        npro = npro + 1
        ikmid = -1
        DO ik = 1, nks(ir)-2
          IF (rs(ik,ir).GT.r0.AND.
c     .        ((zs(ik,ir).GE.0.64.AND.zs(ik+1,ir).LT.0.64).OR.
c     .         (zs(ik,ir).LT.0.64.AND.zs(ik+1,ir).GE.0.64))) ikmid = ik
     .        ((zs(ik,ir).GE.z0.AND.zs(ik+1,ir).LT.z0).OR.
     .         (zs(ik,ir).LT.z0.AND.zs(ik+1,ir).GE.z0))) ikmid = ik
        ENDDO
        ikmid1(npro) = ikmid
        ring  (npro) = ir
        r     (npro) = rs     (ikmid,ir)
        rho1  (npro) = rho    (ir,CELL1) 
        psin  (npro) = psitarg(ir,2)
        DO ik = 1, nks(ir)-1
         tvolp(npro,0) = tvolp(npro,0)+kvols(ik,ir)                  ! Volume
         ncol = 0
         DO iz = 1, nizs
           tvolp(npro,iz)=tvolp(npro,iz)+kvols(ik,ir)*sdlims(ik,ir,iz)
         ENDDO
         tvolp(npro,nizs+1)=tvolp(npro,nizs+1)+kvols(ik,ir)*knbs (ik,ir)  ! ne
         tvolp(npro,nizs+2)=tvolp(npro,nizs+2)+kvols(ik,ir)*ktebs(ik,ir)  ! Te
         tvolp(npro,nizs+3)=tvolp(npro,nizs+3)+kvols(ik,ir)*ktibs(ik,ir)  ! D+ temperature
        ENDDO
        DO iz = 1, nizs
          tvolp(npro,nizs+4)=tvolp(npro,nizs+4)+tvolp(npro,iz)            ! Total amount of impurity on ring
          tvolp(npro,nizs+5)=tvolp(npro,nizs+5)+tvolp(npro,iz)*REAL(iz)   ! Total amount * charge on ring
        ENDDO
        avolpro(npro,1:nizs+5) = tvolp(npro,1:nizs+5) / tvolp(npro,0)     ! Divide by the total volume on the ring to get average values
      ENDDO
      avolpro(:,1     :nizs  )=avolpro(:,1     :nizs  )*impurity_influx
      avolpro(:,nizs+4:nizs+5)=avolpro(:,nizs+4:nizs+5)*impurity_influx
      avolpro(1:npro,nizs+6) = avolpro(1:npro,nizs+4) /  ! Fraction of number
     .                         avolpro(1:npro,nizs+1)  
      avolpro(1:npro,nizs+7) = avolpro(1:npro,nizs+5) /  ! Fraction of postive charge 
     .                         avolpro(1:npro,nizs+1)  
c...  Zeff:
      npro = 0
      DO ir = 2, irsep-1
        npro = npro + 1
        avolpro(npro,nizs+8) = avolpro(npro,nizs+1)  ! Values for Z=1
        avolpro(npro,nizs+9) = avolpro(npro,nizs+1)
        DO iz = 1, nizs                                              ! This isn quite right since Z=1 for impurities should be folded
          avolpro(npro,nizs+8) = avolpro(npro,nizs+8) +              ! into the above initial value, but since Z=1 for impurities
     .                           avolpro(npro,iz    ) * REAL(iz)**2  ! is zero in the core (!) it doesn't matter...
          avolpro(npro,nizs+9) = avolpro(npro,nizs+9) +
     .                           avolpro(npro,iz    ) * REAL(iz)  
        ENDDO
        avolpro(npro,nizs+8) = avolpro(npro,nizs+8)/avolpro(npro,nizs+9) 
      ENDDO

      WRITE(fp,*) '* Volume averaged core impurity radial profiles:'
      WRITE(fp,'(A4,3A10,2X,5A10,2X,80I10)') 
     .  '*','rho','psin','vol','frac_%','frac_e_%','Zeff','ne^20','Te',
     .  (i1,i1=1,nizs)
      WRITE(fp,*) npro
      DO i1 = 1, npro
        WRITE(fp,'(I4,3F10.5,2X,3F10.6,2F10.2,1P,2X,80E10.2)') 
     .    i1,rho1(i1),psin(i1),tvolp(i1,0),
     .    avolpro(i1,nizs+6)*100.0,avolpro(i1,nizs+7)*100.0,
     .    avolpro(i1,nizs+8),
     .    avolpro(i1,nizs+1)*1.0E-20,avolpro(i1,nizs+2),
     .    (avolpro(i1,iz),iz=1,nizs),
     .    avolpro(i1,nizs+4)
      ENDDO

      CALL inOpenInterface('osm.idl.core_impurities',ITF_WRITE)
 
      CALL inPutData(absfac       ,'DIV_IMPURITY_INFLUX','m-1 s-1')
      CALL inPutData(eirene_influx,'EIR_IMPURITY_INFLUX','m-1 s-1')
      CALL inPutData(cizsc,'IMP_INITIAL_IZ','NA')
      CALL inPutData(nizs ,'IMP_MAX_IZ'    ,'NA')
      CALL inPutData(cion ,'IMP_Z'         ,'NA')
      CALL inPutData(crmi ,'IMP_A'         ,'NA')
      CALL inPutData(ring  (1:npro)                 ,'RING'   ,'NA')
      CALL inPutData(r     (1:npro)                 ,'MID_R'  ,'m' )
      CALL inPutData((r    (1:npro)-r0)/(r(npro)-r0),'MID_R/A','m' )  ! Not quite right but almost...
      CALL inPutData(rho1  (1:npro)                 ,'RHO'    ,'m' )
      CALL inPutData(psin  (1:npro)                 ,'PSIN'   ,'NA')
      CALL inPutData(ksmaxs(ring(1:npro))           ,'L'      ,'m' )

      CALL inPutData(avolpro(1:npro,nizs+1),'AVG_NE','m-3')
      CALL inPutData(avolpro(1:npro,nizs+2),'AVG_TE','eV')
      CALL inPutData(avolpro(1:npro,nizs+3),'AVG_TI','eV')

      CALL inPutData(tvolp  (1:npro,0     )      ,'IMP_VOL'     ,'m-3')
      CALL inPutData(avolpro(1:npro,nizs+6)*100.0,'IMP_FRAC_%'  ,'NA')
      CALL inPutData(avolpro(1:npro,nizs+7)*100.0,'IMP_FRAC_E_%','NA')
      CALL inPutData(avolpro(1:npro,nizs+8)      ,'IMP_ZEFF'    ,'NA')
      DO iz = 1, nizs
       WRITE(dummy,'(A,I0.2)') 'IMP_AVG_DENS_',iz
       CALL inPutData(avolpro(1:npro,iz),dummy(1:LEN_TRIM(dummy)),'m-3')
      ENDDO

      DO i1 = 1, npro
        machno = kvhs(ikmid1(i1),ring(i1)) / qt / 
     .           GetCs(ktebs(ikmid1(i1),ring(i1)),
     .                 ktibs(ikmid1(i1),ring(i1)))
        pressure = CalcPressure(knbs (ikmid1(i1),ring(i1)),
     .                          ktebs(ikmid1(i1),ring(i1)),
     .                          ktibs(ikmid1(i1),ring(i1)),
     .                          kvhs (ikmid1(i1),ring(i1))/qt)* ECH
        CALL inPutData(knbs (ikmid1(i1),ring(i1))   ,'MID_NE','m-3')
        CALL inPutData(kvhs (ikmid1(i1),ring(i1))/qt,'MID_VB','m s-1')
        CALL inPutData(machno                       ,'MID_M' ,'none')
        CALL inPutData(pressure                     ,'MID_P' ,'Pa')
        CALL inPutData(ktebs(ikmid1(i1),ring(i1))   ,'MID_TE','eV')
        CALL inPutData(ktibs(ikmid1(i1),ring(i1))   ,'MID_TI','eV')
      ENDDO

      CALL inCloseInterface
c
c     ------------------------------------------------------------------
c...  Inner midplane profiles:
      npro = 0
      midpro = 0.0  
      scale = 1.0
      DO ir = 2, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        IF (ir.EQ.irwall.OR.ir.EQ.irtrap) CYCLE   ! Not sure why, but seems necessary...
        ikmid = 0
        DO ik = 1, nks(ir)-2
          IF (rs(ik,ir).LT.r0.AND.
     .        ((zs(ik,ir).GE.0.0.AND.zs(ik+1,ir).LT.0.0).OR.
     .         (zs(ik,ir).LT.0.0.AND.zs(ik+1,ir).GE.0.0))) THEN
            ikmid = ik
          ENDIF
        ENDDO
        IF (ikmid.EQ.0) CYCLE  ! Flux tube does not cross the y-axis on the HFS
c        rsmin = 0.0
c        DO ik = 1, nks(ir)-2
c          IF (rs(ik,ir).GT.rsmin.AND.ABS(zs(ik,ir)).LT.0.3*r0) THEN  ! Poorly defined
c            ikmid = ik
c            rsmin = rs(ik,ir)
c          ENDIF
c        ENDDO
        npro = npro + 1
        ring(npro) = ir
        rho1(npro) = rho(ir,CELL1) 
        r(npro) = rs(ikmid,ir)
        psin(npro) = psitarg(ir,2)
        ncol = 0
        DO iz = 5, 25, 5
          ncol = ncol + 1
          IF (ir.EQ.2) iz1(ncol) = iz
          IF (iz.EQ.5) ikmid1(npro) = ikmid
          midpro(npro,ncol) = sdlims(ikmid,ir,iz) * scale
c          WRITE(0,*) ir,iz,sdlims(:,ir,iz)
        ENDDO
      ENDDO
      WRITE(fp,*) '* Inner midplane profiles (Z approx. 0.0):'
      WRITE(fp,'(A4,8X,4A9,2X,A9,A10,A6,A9,2A9,10I10)') 
     .  '*','r','rho','psin','L','ne','vb','M','P',
     .  'Te','Ti',(iz1(i1),i1=1,ncol)
      WRITE(fp,*) npro
      DO i1 = 1, npro
        machno = kvhs(ikmid1(i1),ring(i1)) / qt /
     .           GetCs(ktebs(ikmid1(i1),ring(i1)),
     .                 ktibs(ikmid1(i1),ring(i1)))
        pressure = CalcPressure(knbs (ikmid1(i1),ring(i1)),
     .                          ktebs(ikmid1(i1),ring(i1)),
     .                          ktibs(ikmid1(i1),ring(i1)),
     .                          kvhs (ikmid1(i1),ring(i1))/qt)* ECH
        WRITE(fp,'(3I4,3F9.5,F9.2,2X,'//
     .           '1P,E9.2,E10.2,0P,F6.2,1P,E9.2,0P'//
     .           '2F9.2,1P,9E10.2,2X,I4)') 
     .    i1,ikmid1(i1),ring(i1),
     .    r(i1),rho1(i1),psin(i1),ksmaxs(ring(i1)),
     .    knbs(ikmid1(i1),ring(i1)),
     .    kvhs(ikmid1(i1),ring(i1)) / qt,
     .    machno,pressure,
     .    ktebs(ikmid1(i1),ring(i1)),ktibs(ikmid1(i1),ring(i1)),
     .    (midpro(i1,i2),i2=1,ncol)
      ENDDO



      WRITE(fp,*)
      IF (fp.NE.0) CLOSE(fp)


   
c
c     *** STOPPING HERE FOR NOW! ***
c
      RETURN



c...  Get T_D and T_D2 (crap):
      CALL LoadTriangles
      ALLOCATE(tdata(ntri))

      t_d = 0.0
      t_d2 = 0.0

      GOTO 10  ! *** Having to skip this block for now, crashing...

      tdata = 0.0
      CALL LoadTriangleData(2,1,7,0,tdata,'default')  
      DO itri = 1, ntri
        IF (tri(itri)%type.NE.MAGNETIC_GRID) CYCLE        
        ik = tri(itri)%index(1)
        ir = tri(itri)%index(2)
        WRITE(0,*) 'TRI:',itri,ik,ir
        t_d(ik,ir) = t_d(ik,ir) + tdata(itri)
      ENDDO
      t_d = t_d / 2.0

      tdata = 0.0
      CALL LoadTriangleData(3,1,7,0,tdata,'default')  
      DO itri = 1, ntri
        IF (tri(itri)%type.NE.MAGNETIC_GRID) CYCLE        
        ik = tri(itri)%index(1)
        ir = tri(itri)%index(2)
        t_d2(ik,ir) = t_d2(ik,ir) + tdata(itri)
      ENDDO
      t_d2 = t_d2 / 2.0

 10   CONTINUE

c...  Inner midplane profiles:
      npro = 0
      midpro = 0.0  
      DO ir = 2, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        ikmid = 0
        DO ik = 1, nks(ir)-2
          IF (rs(ik,ir).LT.r0.AND.
     .        ((zs(ik,ir).GE.0.0.AND.zs(ik+1,ir).LT.0.0).OR.
     .         (zs(ik,ir).LT.0.0.AND.zs(ik+1,ir).GE.0.0))) THEN
            ikmid = ik
            IF (ABS(zs(ik,ir)).GT.ABS(zs(ik+1,ir))) ikmid = ik + 1
          ENDIF
        ENDDO
        IF (ikmid.NE.0) THEN
          npro = npro + 1
          psin(npro) = psitarg(ir,2)
          r(npro) = rs(ikmid,ir)
          midpro(npro,1) = pinion  (ikmid,ir)  ! Ionisation
          midpro(npro,2) = pinalpha(ikmid,ir)  ! Dalpha
          midpro(npro,3) = pinatom (ikmid,ir)  ! D density
          midpro(npro,4) = t_D     (ikmid,ir)  ! D temperature
          midpro(npro,5) = pinmol  (ikmid,ir)  ! D2 density
          midpro(npro,6) = t_D2    (ikmid,ir)  ! D2 temperature
          midpro(npro,7) = knbs    (ikmid,ir)  ! ne
          midpro(npro,8) = ktebs   (ikmid,ir)  ! Te
          midpro(npro,9) = ktibs   (ikmid,ir)  ! D+ temperature
        ENDIF
      ENDDO

      WRITE(fp,*) '* Inner midplane profiles:'
      WRITE(fp,'(A4,2A8,9A10)') 
     .  '*','r','psin','Ionis.','Dalpha','n_D','T_D','n_D2','T_D2',
     .  'n_e','T_e','T_D+'
      WRITE(fp,*) npro
      DO i1 = 1, npro
        WRITE(fp,'(I4,2F8.5,1P,9E10.2)') 
     .    i1,r(i1),psin(i1),midpro(i1,1:9)
      ENDDO

c...  Outer midplane profiles:
      npro = 0
      midpro = 0.0  
      DO ir = 2, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        IF (ir.EQ.irwall.OR.ir.EQ.irtrap) CYCLE   ! Not sure why, but seems necessary...
        ikmid = 0
        DO ik = 1, nks(ir)-2
          IF (rs(ik,ir).GT.r0.AND.
     .        ((zs(ik,ir).GE.0.0.AND.zs(ik+1,ir).LT.0.0).OR.
     .         (zs(ik,ir).LT.0.0.AND.zs(ik+1,ir).GE.0.0))) THEN
            ikmid = ik
            IF (ABS(zs(ik,ir)).GT.ABS(zs(ik+1,ir))) ikmid = ik + 1
          ENDIF
        ENDDO
        IF (ikmid.NE.0) THEN
          npro = npro + 1
          psin(npro) = psitarg(ir,2)
          r(npro) = rs(ikmid,ir)
          midpro(npro,1) = pinion  (ikmid,ir)  ! Ionisation
          midpro(npro,2) = pinalpha(ikmid,ir)  ! Dalpha
          midpro(npro,3) = pinatom (ikmid,ir)  ! D density
          midpro(npro,4) = t_D     (ikmid,ir)  ! D temperature
          midpro(npro,5) = pinmol  (ikmid,ir)  ! D2 density
          midpro(npro,6) = t_D2    (ikmid,ir)  ! D2 temperature
          midpro(npro,7) = knbs    (ikmid,ir)  ! ne
          midpro(npro,8) = ktebs   (ikmid,ir)  ! Te
          midpro(npro,9) = ktibs   (ikmid,ir)  ! D+ temperature
        ENDIF
      ENDDO
      WRITE(fp,*) '* Outer midplane profiles:'
      WRITE(fp,'(A4,2A8,9A10)') 
     .  '*','r','psin','Ionis.','Dalpha','n_D','T_D','n_D2','T_D2',
     .  'n_e','T_e','T_D+'
      WRITE(fp,*) npro
      DO i1 = 1, npro
        WRITE(fp,'(I4,2F8.5,1P,9E10.2)') 
     .    i1,r(i1),psin(i1),midpro(i1,1:9)
      ENDDO

c...  Volume averaged radial profiles:
c      IF (nrs.EQ.65) nre = 27  ! sonnet_15169_259
c      IF (nrs.EQ.72) nre = 36  ! 
      nre = irsep - 1
      npro = nre - 1 
      avolpro = 0.0  ! Full volume averaged radial profiles
      tvolp = 0.0  ! Totals
      DO ir = 2, nre
        psin(ir-1) = psitarg(ir,2)
        deltar(ir-1) = rho(ir,OUT23) - rho(ir,IN14)
c        deltar(ir-1) = (r  (ir        ) - r  (ir-1    )) /   ! Check that rho is okay...
c     .                 (rho(ir+1,CELL1) - rho(ir,CELL1))
        frac = (rho(ir  ,OUT23) - rho(ir,CELL1)) /
     .         (rho(ir+1,CELL1) - rho(ir,CELL1))
        dpsin(ir-1) = (psitarg(ir+1,2) - psitarg(ir,2)) * frac * 2.0
        DO ik = 1, nks(ir)-1
          IF (rs(ik,ir).GT.rxp) CYCLE
          tvolp(ir-1,0) = tvolp(ir-1,0)+kvols(ik,ir)                  ! Volume
          tvolp(ir-1,1) = tvolp(ir-1,1)+kvols(ik,ir)*pinion  (ik,ir)  ! Ionisation
          tvolp(ir-1,2) = tvolp(ir-1,2)+kvols(ik,ir)*pinalpha(ik,ir)  ! Dalpha
          tvolp(ir-1,3) = tvolp(ir-1,3)+kvols(ik,ir)*pinatom (ik,ir)  ! D density
          tvolp(ir-1,4) = tvolp(ir-1,4)+kvols(ik,ir)*t_D     (ik,ir)  ! D temperature
          tvolp(ir-1,5) = tvolp(ir-1,5)+kvols(ik,ir)*pinmol  (ik,ir)  ! D2 density
          tvolp(ir-1,6) = tvolp(ir-1,6)+kvols(ik,ir)*t_D2    (ik,ir)  ! D2 temperature
          tvolp(ir-1,7) = tvolp(ir-1,7)+kvols(ik,ir)*knbs    (ik,ir)  ! ne
          tvolp(ir-1,8) = tvolp(ir-1,8)+kvols(ik,ir)*ktebs   (ik,ir)  ! Te
          tvolp(ir-1,9) = tvolp(ir-1,9)+kvols(ik,ir)*ktibs   (ik,ir)  ! D+ temperature
        ENDDO
        avolpro(ir-1,1:9) = tvolp(ir-1,1:9) / tvolp(ir-1,0)
      ENDDO
      WRITE(fp,*) '* Total inner radial core profiles:'
      WRITE(fp,'(A4,4A8,10A10)') 
     .  '*','r','psin','dr','dpsin','Vol','Ionis.','Dalpha','n_D','T_D',
     .  'n_D2','T_D2','n_e','T_e','T_D+'
      WRITE(fp,*) npro
      DO i1 = 1, npro
        WRITE(fp,'(I4,4F8.5,1P,10E10.2)') 
     .    i1,r(i1),psin(i1),deltar(i1),dpsin(i1),tvolp(i1,0:9)
      ENDDO
      WRITE(fp,'(4X,32X,1P,10E10.2)') 
     .  (SUM(tvolp(1:npro,i1)),i1=0,9)

      avolpro = 0.0  ! Full volume averaged radial profiles
      tvolp = 0.0  ! Totals
      DO ir = 2, nre
        DO ik = 1, nks(ir)-1
          IF (rs(ik,ir).LT.rxp) CYCLE
          tvolp(ir-1,0) = tvolp(ir-1,0)+kvols(ik,ir)                  ! Volume
          tvolp(ir-1,1) = tvolp(ir-1,1)+kvols(ik,ir)*pinion  (ik,ir)  ! Ionisation
          tvolp(ir-1,2) = tvolp(ir-1,2)+kvols(ik,ir)*pinalpha(ik,ir)  ! Dalpha
          tvolp(ir-1,3) = tvolp(ir-1,3)+kvols(ik,ir)*pinatom (ik,ir)  ! D density
          tvolp(ir-1,4) = tvolp(ir-1,4)+kvols(ik,ir)*t_D     (ik,ir)  ! D temperature
          tvolp(ir-1,5) = tvolp(ir-1,5)+kvols(ik,ir)*pinmol  (ik,ir)  ! D2 density
          tvolp(ir-1,6) = tvolp(ir-1,6)+kvols(ik,ir)*t_D2    (ik,ir)  ! D2 temperature
          tvolp(ir-1,7) = tvolp(ir-1,7)+kvols(ik,ir)*knbs    (ik,ir)  ! ne
          tvolp(ir-1,8) = tvolp(ir-1,8)+kvols(ik,ir)*ktebs   (ik,ir)  ! Te
          tvolp(ir-1,9) = tvolp(ir-1,9)+kvols(ik,ir)*ktibs   (ik,ir)  ! D+ temperature
        ENDDO
        avolpro(ir-1,1:9) = tvolp(ir-1,1:9) / tvolp(ir-1,0)
      ENDDO
      WRITE(fp,*) '* Total outer radial core profiles:'
      WRITE(fp,'(A4,4A8,10A10)') 
     .  '*','r','psin','dr','psin','Vol','Ionis.','Dalpha','n_D','T_D',
     .  'n_D2','T_D2','n_e','T_e','T_D+'
      WRITE(fp,*) npro
      DO i1 = 1, npro
        WRITE(fp,'(I4,4F8.5,1P,10E10.2)') 
     .    i1,r(i1),psin(i1),deltar(i1),dpsin(i1),tvolp(i1,0:9)
      ENDDO
      WRITE(fp,'(4X,32X,1P,10E10.2)') 
     .  (SUM(tvolp(1:npro,i1)),i1=0,9)

c...  Volume averaged poloidal profiles:
      nre = irsep - 1
      npro = nks(nre) - 1
      avolpro = 0.0  ! Full volume averaged radial profiles
      tvolp = 0.0  ! Totals
      DO ik = 1, nks(nre) - 1
        id = korpg(ik,nre)
        r1 = 0.5 * (rvertp(1,id) + rvertp(2,id))
        r2 = 0.5 * (rvertp(3,id) + rvertp(4,id))
        z1 = 0.5 * (zvertp(1,id) + zvertp(2,id))
        z2 = 0.5 * (zvertp(3,id) + zvertp(4,id))
        dp(ik) = SQRT((r2-r1)**2 + (z2-z1)**2)
        IF (ik.EQ.1) THEN
          p(ik) = 0.5 * dp(ik)
        ELSE
          p(ik) = p(ik-1) + 0.5 * dp(ik)
        ENDIF
        DO ir = 2, irsep - 1
          tvolp(ik,0) = tvolp(ik,0)+kvols(ik,ir)                  ! Volume
          tvolp(ik,1) = tvolp(ik,1)+kvols(ik,ir)*pinion  (ik,ir)  ! Ionisation
          tvolp(ik,2) = tvolp(ik,2)+kvols(ik,ir)*pinalpha(ik,ir)  ! Dalpha
          tvolp(ik,3) = tvolp(ik,3)+kvols(ik,ir)*pinatom (ik,ir)  ! D density
          tvolp(ik,4) = tvolp(ik,4)+kvols(ik,ir)*t_D     (ik,ir)  ! D temperature
          tvolp(ik,5) = tvolp(ik,5)+kvols(ik,ir)*pinmol  (ik,ir)  ! D2 density
          tvolp(ik,6) = tvolp(ik,6)+kvols(ik,ir)*t_D2    (ik,ir)  ! D2 temperature
          tvolp(ik,7) = tvolp(ik,7)+kvols(ik,ir)*knbs    (ik,ir)  ! ne
          tvolp(ik,8) = tvolp(ik,8)+kvols(ik,ir)*ktebs   (ik,ir)  ! Te
          tvolp(ik,9) = tvolp(ik,9)+kvols(ik,ir)*ktibs   (ik,ir)  ! D+ temperature
        ENDDO
        avolpro(ir-1,1:9) = tvolp(ir-1,1:9) / tvolp(ir-1,0)
      ENDDO
      pmax = p(ik-1) + 0.5 * dp(ik-1)
      WRITE(fp,*) '* Total poloidal core profiles:'
      WRITE(fp,'(A4,3A8,10A10)') 
     .  '*','p','dp','fp','Vol','Ionis.','Dalpha','n_D','T_D',
     .  'n_D2','T_D2','n_e','T_e','T_D+'
      WRITE(fp,*) npro
      DO i1 = 1, npro
        WRITE(fp,'(I4,3F8.5,1P,10E10.2)') 
     .    i1,p(i1),dp(i1),dp(i1)/pmax,tvolp(i1,0:9)
      ENDDO
      WRITE(fp,'(4X,24X,1P,10E10.2)') 
     .  (SUM(tvolp(1:npro,i1)),i1=0,9)

c...  Volume averaged radial profiles:  *** THIS APPEARS VERY SIMILAR TO THE SIMILARLY LABELLED CODE THAT APPEARS ABOVE ***
      nre = irsep - 1
      npro = nre - 1 
      avolpro = 0.0  ! Full volume averaged radial profiles
      tvolp = 0.0  ! Totals
      DO ir = 2, nre
        rho1(ir-1) = rho(ir,CELL1) 
        psin(ir-1) = psitarg(ir,2)
        deltar(ir-1) = rho(ir,OUT23) - rho(ir,IN14)
        frac = (rho(ir  ,OUT23) - rho(ir,CELL1)) /
     .         (rho(ir+1,CELL1) - rho(ir,CELL1))
        dpsin(ir-1) = (psitarg(ir+1,2) - psitarg(ir,2)) * frac * 2.0
        DO ik = 1, nks(ir)-1
          tvolp(ir-1,0) = tvolp(ir-1,0)+kvols(ik,ir)                  ! Volume
          tvolp(ir-1,1) = tvolp(ir-1,1)+kvols(ik,ir)*pinion  (ik,ir)  ! Ionisation
          tvolp(ir-1,2) = tvolp(ir-1,2)+kvols(ik,ir)*pinalpha(ik,ir)  ! Dalpha
          tvolp(ir-1,3) = tvolp(ir-1,3)+kvols(ik,ir)*pinatom (ik,ir)  ! D density
          tvolp(ir-1,4) = tvolp(ir-1,4)+kvols(ik,ir)*t_D     (ik,ir)  ! D temperature
          tvolp(ir-1,5) = tvolp(ir-1,5)+kvols(ik,ir)*pinmol  (ik,ir)  ! D2 density
          tvolp(ir-1,6) = tvolp(ir-1,6)+kvols(ik,ir)*t_D2    (ik,ir)  ! D2 temperature
          tvolp(ir-1,7) = tvolp(ir-1,7)+kvols(ik,ir)*knbs    (ik,ir)  ! ne
          tvolp(ir-1,8) = tvolp(ir-1,8)+kvols(ik,ir)*ktebs   (ik,ir)  ! Te
          tvolp(ir-1,9) = tvolp(ir-1,9)+kvols(ik,ir)*ktibs   (ik,ir)  ! D+ temperature
        ENDDO
        avolpro(ir-1,1:9) = tvolp(ir-1,1:9) / tvolp(ir-1,0)
      ENDDO
      WRITE(fp,*) '* Total inner radial core profiles:'
      WRITE(fp,'(A4,5A8,10A10)') 
     .  '*','r','dr','psin','dpsin','rho','Vol','Ionis.','Dalpha',
     .  'n_D','T_D','n_D2','T_D2','n_e','T_e','T_D+'
      WRITE(fp,*) npro
      DO i1 = 1, npro
        WRITE(fp,'(I4,5F8.5,1P,10E10.2)') 
     .    i1,r(i1),deltar(i1),psin(i1),dpsin(i1),rho1(i1),
     .    tvolp(i1,0:9)
      ENDDO
      WRITE(fp,'(4X,40X,1P,10E10.2)') 
     .  (SUM(tvolp(1:npro,i1)),i1=0,9)

c...  Clear arrays:
      CALL DEALLOC_VERTEX  
      CALL DEALLOC_SURFACE   
      CALL DEALLOC_TRIANGLE
      IF (ALLOCATED(tdata)) DEALLOCATE(tdata)


      RETURN
 99   STOP
      END

c
c ======================================================================
c
c
      SUBROUTINE DTSAnalysis
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slout'
      INCLUDE 'comgra'

      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost

c      INTEGER MAXGXS,MAXNGS

      INTEGER, PARAMETER :: MAXPTS_T = 50, MAXPTS_X = 25, MAXCHANNEL = 3

      INTEGER ic,ix,it,nt,nx,nc,nc1,nc2,count,interpolation
      REAL    a_n(0:MAXCHANNEL,6),a_t(0:MAXCHANNEL,6),
     .        t(MAXPTS_T),deltat,t1,t2,
     .        x(MAXPTS_X),deltax,x1,x2,
     .        n_e(MAXPTS_X,MAXPTS_T),t_e(MAXPTS_X,MAXPTS_T),
     .        sig(MAXPTS_X,MAXPTS_T),sig_ex(MAXPTS_X,MAXPTS_T),
     .        c_n(0:MAXCHANNEL,MAXPTS_T),c_t(0:MAXCHANNEL,MAXPTS_T),
     .        c_x(0:MAXCHANNEL),qmin,qmax,
     .        w(MAXPTS_X,MAXPTS_T),wsum,avg_n_e,avg_t_e,rv(4),zv(4),
     .        ms_avg_sig(MAXPTS_X),navg,tavg,ms_n_e(MAXPTS_X),
     .        ms_t_e(MAXPTS_X),ms_tmin,ms_tval,ms_sig(MAXPTS_X),
     .        ms_sig_ex(MAXPTS_X)

      CHARACTER xlabel*36,ylabel*64,label*5

c...  ADAS:
      CHARACTER dummy*80,adasid*80,adasex*3
      INTEGER   adasyr,isele,iselr,iselx,iseld,ierr,ircode
      REAL      pecvals(MAXGXS,MAXNGS),wlngth,
     .          tevals(MAXGXS,MAXNGS),nevals(MAXGXS,MAXNGS)



      DATA a_n /  0.0        ,  0.0        , -1.356E+18,  0.0,  
     .           -2.96446E+18, -5.15030E+18,  1.756E+19,  0.0,  
     .            3.47025E+19,  7.13752E+19, -6.441E+19, -1.6953E+18,
     .           -1.20367E+20, -2.94906E+20,  5.406E+19,  3.5908E+19,
     .            9.98718E+19,  3.01845E+20, -7.449E+19, -1.8084E+20,
     .            1.00793E+20,  3.97685E+20,  9.150E+20,  6.4927E+20 /

      DATA a_t / 0.0,  0.01518,  0.00335,  0.0,
     .           0.0, -0.24115, -0.05540,  0.0,
     .           0.0,  1.38015,  0.33977,  0.04335,
     .           0.0, -3.29071, -0.86869, -0.27682, 
     .           0.0,  2.56044,  0.58925,  0.07477,
     .           0.6,  1.86015,  1.49758,  2.34373 /

c      DATA a_n / 1.0E+20, 1.1E+20, 1.2E+20, 1.3E+20, 1.4E+20,
c     .           2.0E+20, 1.0E+20, 2.2E+20, 2.3E+20, 2.4E+20,
c     .           3.0E+20, 3.1E+20, 1.0E+20, 3.3E+20, 3.4E+20,
c     .           2.0E+20, 2.1E+20, 2.2E+20, 1.0E+20, 2.4E+20,
c     .           1.0E+20, 2.0E+20, 3.0E+20, 4.0E+20, 5.0E+20,
c     .           5.0E+20, 5.0E+20, 5.0E+20, 5.0E+20, 6.0E+20 /

c      DATA a_t / 1.0, 1.1, 1.2, 1.3, 1.4,
c     .           2.0, 2.0, 2.2, 2.3, 2.4,
c     .           3.0, 3.1, 3.0, 3.3, 3.4,
c     .           2.0, 2.1, 2.2, 4.0, 2.4,
c     .           1.0, 1.0, 1.0, 1.0, 1.0,
c     .           1.0, 1.0, 1.0, 1.0, 1.0 /

      DATA c_x / -1.346, -1.306, -1.266, -1.210 /


c...  Time series:
      nx = MAXPTS_X ! 40
      x1 = -1.346
      x2 = -1.210
      deltax = ABS(x2 - x1) / REAL(nx - 1)
      x(1) = x1
      DO ix = 2, nx
        x(ix) = x(ix-1) + deltax          
      ENDDO

c...  Space series:
      nt = MAXPTS_T ! 46
      t1 = 0.0
      t2 = 4.6
      deltat = ABS(t2 - t1) / REAL(nt - 1)
      t(1) = t1
      DO it = 2, nt
        t(it) = t(it-1) + deltat          
      ENDDO

c...  n(t) and T(t) at DTS points:
      nc1 = 0
      nc2 = 3
      DO ic = nc1, nc2
        DO it = 1, nt        
c          c_n(ic,it) = (t2 - t(it)) / (t2 - t1) * a_n(ic,5) +
c     .                 (t(it) - t1) / (t2 - t1) * a_n(ic,6)
c          c_t(ic,it) = (t2 - t(it)) / (t2 - t1) * a_t(ic,5) +
c     .                 (t(it) - t1) / (t2 - t1) * a_t(ic,6)

          c_n(ic,it) = a_n(ic,1) * t(it)**5 + a_n(ic,2) * t(it)**4 + 
     .                 a_n(ic,3) * t(it)**3 + a_n(ic,4) * t(it)**2 + 
     .                 a_n(ic,5) * t(it)    + a_n(ic,6)
          c_t(ic,it) = a_t(ic,1) * t(it)**5 + a_t(ic,2) * t(it)**4 + 
     .                 a_t(ic,3) * t(it)**3 + a_t(ic,4) * t(it)**2 + 
     .                 a_t(ic,5) * t(it)    + a_t(ic,6)
        ENDDO
      ENDDO

c      c_n = c_n * 0.5


c      DO it = 1, nt
c        WRITE(0,*) 'TIME:',t(it)
c        DO ic = nc1, nc2
c          WRITE(0,*) c_x(ic),c_n(ic,it),c_t(ic,it)
c        ENDDO
c      ENDDO
c      STOP 'HAlting..'

c      WRITE(0,*) 'T:'
c      WRITE(0,*) t

c      WRITE(0,*) 'A_N:'
c      WRITE(0,*) a_n(nc1:nc2,6)

c      WRITE(0,*) 'CN0:'
c      WRITE(0,*) c_n(0,1:nt)
c      WRITE(0,*) c_n(1,1:nt)
c      WRITE(0,*) c_n(2,1:nt)
c      WRITE(0,*) c_n(3,1:nt)
c      WRITE(0,*) c_n(4,1:nt)

c...  Interpolate to get n(x,t) and T(x,t):
      interpolation = 3
      SELECTCASE (interpolation)
        CASE (1)      ! Linear average
          DO it = 1, nt
c            WRITE(0,*) 'FIT:'
c            WRITE(0,*) c_x(nc1:nc2)
c            WRITE(0,*) c_n(nc1:nc2,it)
            CALL Fitter(nc2-nc1+1,c_x(nc1),c_n(nc1,it),
     .                  nx,x,n_e(1,it),'LINEAR')
            CALL Fitter(nc2-nc1+1,c_x(nc1),c_t(nc1,it),
     .                  nx,x,t_e(1,it),'LINEAR')
          ENDDO
        CASE (2)      ! Spline average
          DO it = 1, nt
            CALL Fitter(nc2-nc1+1,c_x(nc1),c_n(nc1,it),
     .                  nx,x,n_e(1,it),'SPLINE')
            CALL Fitter(nc2-nc1+1,c_x(nc1),c_t(nc1,it),
     .                  nx,x,t_e(1,it),'SPLINE')
          ENDDO
        CASE (3)      ! Linear, linear oscillation
          DO it = 1, nt
            CALL Fitter(nc2-nc1+1,c_x(nc1),c_n(nc1,it),
     .                  nx,x,n_e(1,it),'SPLINE')
            CALL Fitter(nc2-nc1+1,c_x(nc1),c_t(nc1,it),
     .                  nx,x,t_e(1,it),'LINEAR')
          ENDDO
        CASE (4)      ! Spline, linear oscillation
          DO it = 1, nt
            CALL Fitter(nc2-nc1+1,c_x(nc1),c_n(nc1,it),
     .                  nx,x,n_e(1,it),'LINEAR')
            CALL Fitter(nc2-nc1+1,c_x(nc1),c_t(nc1,it),
     .                  nx,x,t_e(1,it),'SPLINE')
          ENDDO
        CASE (5)      ! Linear, sinusoidal oscillation
        CASE (6)      ! Spline, sinusoidal oscillation
        CASEDEFAULT
      ENDSELECT

c      WRITE(0,*) 'X:'
c      WRITE(0,*) x

c      WRITE(0,*) 'C_X:'
c      WRITE(0,*) c_x

c      WRITE(0,*) 'NE0:'
c      DO ix = 1, nx
c        WRITE(0,*) (n_e(ix,it),it=1,nt)
c      ENDDO
c      WRITE(0,*) 'TE0:'
c      DO ix = 1, nx
c        WRITE(0,*) (t_e(ix,it),it=1,nt)
c      ENDDO


c...  Get ADAS data:
c     '000   ADAS' '*' 96  'pju'  12   78   0  0
c     '000   Refl'  4  0.4   -1.0 -1.0
      IF (.TRUE.) THEN

        CALL RDG1(dummy,adasid,adasyr,adasex,
     .            isele,iselr,iselx,iseld,ierr)
        IF (ierr.NE.0) 
     .    CALL ER('DTSAnalysis','ADAS probelm',*99)

        DO it = 1, nt        
        
          nevals(1:nx,1) = n_e(1:nx,it)
          tevals(1:nx,1) = t_e(1:nx,it)
c          nevals(1:nx,1) = n_e((it-1)*nx+1:it*nx)
c          tevals(1:nx,1) = t_e((it-1)*nx+1:it*nx)
          pecvals = 0.0

          adasid(1:1) = '*'

c          WRITE(0,*) 'ADAS:',adasyr,adasex,isele,iselr,iselx
c          WRITE(0,*) 'ADAS:',nx,tevals(1:nx,1),nevals(1:nx,1)

c...      Recombination:
          pecvals = 0.0
          CALL LdADAS2(1,0,adasid,adasyr,adasex,isele,iselr,iselx,  ! PEC data returned    1,1 = cz,iz (hydrogen)
     .                 nx,tevals,1,nevals,2,pecvals,wlngth,ircode)  ! 2 = recombination

          sig   (1:nx,it) = (pecvals(1:nx,1) * nevals(1:nx,1)) * 
     .                                         nevals(1:nx,1)

c...      Excitation:
          pecvals = 0.0
          CALL LdADAS2(1,0,adasid,adasyr,adasex,isele,iselr,iselx,  ! PEC data returned    1,1 = cz,iz (hydrogen)
     .                 nx,tevals,1,nevals,1,pecvals,wlngth,ircode)  ! 2 = recombination

          sig_ex(1:nx,it) = (pecvals(1:nx,1) * nevals(1:nx,1)) * 
     .                                         nevals(1:nx,1)


c          WRITE(0,*) 'ADAS:',adasid
c          WRITE(0,*) 'PEC:'
c          WRITE(0,*) pecvals(1:nx,1)
          IF (it.EQ.1.OR.ircode.NE.0) WRITE(0,*) 'WL:',wlngth,ircode
        ENDDO

      ENDIF

c      sig = 1.0

c      WRITE(0,*) 'SIG:'
c      DO ix = 1, nx
c        WRITE(0,*) (sig(ix,it),it=1,nt)
c      ENDDO


c...  Weight function:
      wsum = 0.0
      DO it = 1, nt
        w(1:nx,it) = sig(1:nx,it) * deltax * deltat
        DO ix = 1, nx
          wsum = wsum + w(ix,it)
        ENDDO
      ENDDO
      w(1:nx,1:nt) = w(1:nx,1:nt) / wsum


      wsum = 0.0
      DO it = 1, nt
        DO ix = 1, nx
          wsum = wsum + w(ix,it)
        ENDDO
      ENDDO
c      WRITE(0,*) 'W:',wsum
c      DO ix = 1, nx
c        WRITE(0,*) (w(ix,it),it=1,nt)
c      ENDDO

c...  Average density and temperature:
      avg_n_e = 0.0
      avg_t_e = 0.0
      DO it = 1, nt
        DO ix = 1, nx
          avg_n_e = avg_n_e + n_e(ix,it) * w(ix,it)
          avg_t_e = avg_t_e + t_e(ix,it) * w(ix,it)
        ENDDO
      ENDDO

c...  Average emission signal vs x:
      ms_tmin = 100.0
      ms_tval = -1.0
      tavg = 0.8 ! 4.25 ! 0.8 ! 3.0
      navg = 0.0
      ms_avg_sig = 0.0
      DO it = 1, nt
        IF (t(it).GE.tavg-0.5.AND.t(it).LE.tavg+0.5) THEN
          WRITE(0,*) 'AVG TIME:',t(it)
          navg = navg + 1.0
          ms_avg_sig(1:nx) = ms_avg_sig(1:nx) + sig(1:nx,it)
        ENDIF
        IF (ABS(t(it)-tavg).LT.ms_tmin) THEN
          ms_n_e   (1:nx) = n_e   (1:nx,it)
          ms_t_e   (1:nx) = t_e   (1:nx,it)
          ms_sig   (1:nx) = sig   (1:nx,it)
          ms_sig_ex(1:nx) = sig_ex(1:nx,it)
          ms_tval = t(it)
          ms_tmin = ABS(t(it)-tavg)
        ENDIF
      ENDDO
      ms_avg_sig(1:nx) = ms_avg_sig(1:nx) / navg

      WRITE(0,*) 'WAVELENGTH, INTER.:',wlngth,interpolation
      WRITE(0,*) 'MS AVG:',tavg,navg,ms_tval
      DO ix = 1, nx
        WRITE(0,'(F8.3,1P,2E10.2,0P,F7.2,1P,2E10.2,0P,F7.2)')
     .    x(ix),ms_avg_sig(ix),ms_n_e(ix),ms_t_e(ix),
     .    ms_sig(ix),ms_sig_ex(ix),
     .    MIN(100.0,ms_sig(ix)/ms_sig_ex(ix))
      ENDDO

      WRITE(0,*) 'AVG_N, AVG_T:',avg_n_e,avg_t_e




c...  Plots:
      slopt = 1
      plottype(1) = 2

      IF (.TRUE.) THEN
        slopt2 = 1
        iopt_ghost = 1 ! 2


c...    ???
        cxmin = t1
        cxmax = t2

        count = 0

        DO ic = nc1, nc2
c        DO ix = 1, nx

c          IF (MOD(REAL(x(ix)),1.0).NE.0.0) CYCLE

          count = count + 1

          map1x = 0.10             ! Add sensitivity to nxbin, nybin...
          map2x = map1x + 0.50 
          map1y = REAL(count-1) * 0.10 + 0.05
          map2y = map1y + 0.10

          cymin =  HI
          cymax = -HI
          DO it = 1, nt
            cymin = MIN(cymin,c_n(ic,it))
            cymax = MAX(cymax,c_n(ic,it))
c            cymin = MIN(cymin,n_e(ix,it))
c            cymax = MAX(cymax,n_e(ix,it))
          ENDDO
          xlabel = 't (s)     '
          IF (count.NE.1) xlabel = 'none      ' 
          ylabel = 'n_e (m-3) '
          CALL GRTSET_TRIM(' ',' ',' ',' ',' ',      ! TITLE,REF,nVIEW,PLANE,glabel,
     .                     cxmin,cxmax,cymin,cymax,
     .                     ' ',xlabel,ylabel,        ! TABLE,XLAB,YLAB,
     .                     0,' ',0,' ',1)            ! 0,smooth,0,ANLY,1)
          WRITE(label,'(I4,1X)') ic
          CALL GRTRAC(t(1:nt),c_n(ic,1:nt),nt,label,'LINE',1)
c          CALL GRTRAC(t(1:nt),n_e(ix,1:nt),nt,label,'LINE',1)
          CALL DrawFrame

          map1x = 0.80             ! Add sensitivity to nxbin, nybin...
          map2x = map1x + 0.50 
          cymin =  HI
          cymax = -HI
          DO it = 1, nt
            cymin = MIN(cymin,c_t(ic,it))
            cymax = MAX(cymax,c_t(ic,it))
          ENDDO
          ylabel = 'T_e (eV) '
          CALL GRTSET_TRIM(' ',' ',' ',' ',' ',      ! TITLE,REF,nVIEW,PLANE,glabel,
     .                     cxmin,cxmax,cymin,cymax,
     .                     ' ',xlabel,ylabel,        ! TABLE,XLAB,YLAB,
     .                     0,' ',0,' ',1)            ! 0,smooth,0,ANLY,1)
          WRITE(label,'(I4,1X)') ic
          CALL GRTRAC(t(1:nt),c_t(ic,1:nt),nt,label,'LINE',1)
          CALL DrawFrame

c          WRITE(0,*) 'DATA:',c_n(ic,1:nt)
c          WRITE(0,*) 'DATA:',n_e(ix,1:nt)

        ENDDO



        CALL FRAME

      ENDIF

      IF (.TRUE.) THEN
c...    ???
        cxmin = t1
        cxmax = t2

        count = 0

        DO ix = 1, nx, nx / 8

          count = count + 1

          map1x = 0.10             ! Add sensitivity to nxbin, nybin...
          map2x = map1x + 0.35 
          map1y = REAL(count-1) * 0.10 + 0.05
          map2y = map1y + 0.10

          cymin =  HI
          cymax = -HI
          DO it = 1, nt
            cymin = MIN(cymin,n_e(ix,it))
            cymax = MAX(cymax,n_e(ix,it))
          ENDDO
          xlabel = 't (s)     '
          IF (count.NE.1) xlabel = 'none      ' 
          ylabel = 'n_e (m-3) '
          CALL GRTSET_TRIM(' ',' ',' ',' ',' ',      ! TITLE,REF,nVIEW,PLANE,glabel,
     .                     cxmin,cxmax,cymin,cymax,
     .                     ' ',xlabel,ylabel,        ! TABLE,XLAB,YLAB,
     .                     0,' ',0,' ',1)            ! 0,smooth,0,ANLY,1)
          WRITE(label,'(F5.2,1X)') x(ix)
          CALL GRTRAC(t(1:nt),n_e(ix,1:nt),nt,label,'LINE',1)
          CALL DrawFrame

          map1x = map2x + 0.10     ! Add sensitivity to nxbin, nybin...
          map2x = map1x + 0.35 
          cymin =  HI
          cymax = -HI
          DO it = 1, nt
            cymin = MIN(cymin,t_e(ix,it))
            cymax = MAX(cymax,t_e(ix,it))
          ENDDO
          ylabel = 'T_e (eV) '
          CALL GRTSET_TRIM(' ',' ',' ',' ',' ',      ! TITLE,REF,nVIEW,PLANE,glabel,
     .                     cxmin,cxmax,cymin,cymax,
     .                     ' ',xlabel,ylabel,        ! TABLE,XLAB,YLAB,
     .                     0,' ',0,' ',1)            ! 0,smooth,0,ANLY,1)
          CALL GRTRAC(t(1:nt),t_e(ix,1:nt),nt,label,'LINE',1)
          CALL DrawFrame

          map1x = map2x + 0.10     ! Add sensitivity to nxbin, nybin...
          map2x = map1x + 0.35 
          cymin =  HI
          cymax = -HI
          DO it = 1, nt
            cymin = MIN(cymin,sig(ix,it))
            cymax = MAX(cymax,sig(ix,it))
          ENDDO
          ylabel = 'sig '
          CALL GRTSET_TRIM(' ',' ',' ',' ',' ',      ! TITLE,REF,nVIEW,PLANE,glabel,
     .                     cxmin,cxmax,cymin,cymax,
     .                     ' ',xlabel,ylabel,        ! TABLE,XLAB,YLAB,
     .                     0,' ',0,' ',1)            ! 0,smooth,0,ANLY,1)
          CALL GRTRAC(t(1:nt),sig(ix,1:nt),nt,label,'LINE',1)
          CALL DrawFrame

        ENDDO

        CALL FRAME

      ENDIF

      IF (.TRUE.) THEN
c...    ???
        cxmin = x1
        cxmax = x2

        count = 0

        DO it = 1, nt, nt / 8

          count = count + 1

          map1x = 0.10             ! Add sensitivity to nxbin, nybin...
          map2x = map1x + 0.35
          map1y = REAL(count-1) * 0.10 + 0.05
          map2y = map1y + 0.10

          cymin =  HI
          cymax = -HI
          DO ix = 1, nx
            cymin = MIN(cymin,n_e(ix,it))
            cymax = MAX(cymax,n_e(ix,it))
          ENDDO
          xlabel = 'x (m)     '
          IF (count.NE.1) xlabel = 'none      ' 
          ylabel = 'n_e (m-3) '
          CALL GRTSET_TRIM(' ',' ',' ',' ',' ',      ! TITLE,REF,nVIEW,PLANE,glabel,
     .                     cxmin,cxmax,cymin,cymax,
     .                     ' ',xlabel,ylabel,        ! TABLE,XLAB,YLAB,
     .                     0,' ',0,' ',1)            ! 0,smooth,0,ANLY,1)
          WRITE(label,'(F5.2,1X)') t(it)
          CALL GRTRAC(x(1:nx),n_e(1:nx,it),nx,label,'LINE',1)
          CALL DrawFrame

          map1x = map2x + 0.10
          map2x = map1x + 0.35 
          cymin =  HI
          cymax = -HI
          DO ix = 1, nx
            cymin = MIN(cymin,t_e(ix,it))
            cymax = MAX(cymax,t_e(ix,it))
          ENDDO
          ylabel = 'T_e (eV) '
          CALL GRTSET_TRIM(' ',' ',' ',' ',' ',      ! TITLE,REF,nVIEW,PLANE,glabel,
     .                     cxmin,cxmax,cymin,cymax,
     .                     ' ',xlabel,ylabel,        ! TABLE,XLAB,YLAB,
     .                     0,' ',0,' ',1)            ! 0,smooth,0,ANLY,1)
          CALL GRTRAC(x(1:nx),t_e(1:nx,it),nx,label,'LINE',1)
          CALL DrawFrame

          map1x = map2x + 0.10
          map2x = map1x + 0.35 
          cymin =  HI
          cymax = -HI
          DO ix = 1, nx
            cymin = MIN(cymin,sig(ix,it))
            cymax = MAX(cymax,sig(ix,it))
          ENDDO
          ylabel = 'sig '
          CALL GRTSET_TRIM(' ',' ',' ',' ',' ',      ! TITLE,REF,nVIEW,PLANE,glabel,
     .                     cxmin,cxmax,cymin,cymax,
     .                     ' ',xlabel,ylabel,        ! TABLE,XLAB,YLAB,
     .                     0,' ',0,' ',1)            ! 0,smooth,0,ANLY,1)
          CALL GRTRAC(x(1:nx),sig(1:nx,it),nx,label,'LINE',1)
          CALL DrawFrame

        ENDDO

        CALL FRAME

      ENDIF

c...  2D plots:
      IF (.TRUE.) THEN

        cxmin = t1 - 0.5 * deltat
        cxmax = t2 + 0.5 * deltat
        cymin = x1 - 0.5 * deltax
        cymax = x2 + 0.5 * deltax

        map1x = 0.10          
        map2x = map1x + 0.55
        map2y = 0.95
        map1y = map2y - 0.20
        xlabel = 'none  '
        ylabel = 'Z (m) '
        CALL GRTSET_TRIM(' ',' ',' ',' ',' ',      ! TITLE,REF,nVIEW,PLANE,glabel,
     .                   cxmin,cxmax,cymin,cymax,
     .                   ' ',xlabel,ylabel,        ! TABLE,XLAB,YLAB,
     .                   0,' ',0,' ',1)            ! 0,smooth,0,ANLY,1)
        CALL PSPACE(map1x,map2x,map1y,map2y)       ! Don't know why this is necessary
        CALL MAP   (cxmin,cxmax,cymin,cymax)
        qmin = 0.0
        qmax = 0.0
        DO it = 1, nt
          DO ix = 1, nx
            qmax = MAX(qmax,n_e(ix,it))
          ENDDO
        ENDDO             
        DO it = 1, nt
          DO ix = 1, nx
            rv(1) = t(it) - 0.5 * deltat
            zv(1) = x(ix) + 0.5 * deltax
            rv(2) = t(it) - 0.5 * deltat
            zv(2) = x(ix) - 0.5 * deltax
            rv(3) = t(it) + 0.5 * deltat
            zv(3) = x(ix) - 0.5 * deltax
            rv(4) = t(it) + 0.5 * deltat
            zv(4) = x(ix) + 0.5 * deltax
            CALL SetCol255_04(2,n_e(ix,it),qmin,qmax) 
            CALL FILCOL(255)                      
            CALL LINCOL(255)                      
            CALL PTPLOT(rv(1),zv(1),1,4,1)
          ENDDO
        ENDDO
        CALL DrawColourScale(1,2,qmin,qmax,'none')
        CALL DrawFrame

        map1x = 0.10          
        map2x = map1x + 0.55
        map2y = 0.75
        map1y = map2y - 0.20
        xlabel = 'none  '
        ylabel = 'Z (m) '
        CALL GRTSET_TRIM(' ',' ',' ',' ',' ',      ! TITLE,REF,nVIEW,PLANE,glabel,
     .                   cxmin,cxmax,cymin,cymax,
     .                   ' ',xlabel,ylabel,        ! TABLE,XLAB,YLAB,
     .                   0,' ',0,' ',1)            ! 0,smooth,0,ANLY,1)
        CALL PSPACE(map1x,map2x,map1y,map2y)
        CALL MAP   (cxmin,cxmax,cymin,cymax)
        qmin = 0.0
        qmax = 0.0
        DO it = 1, nt
          DO ix = 1, nx
            qmax = MAX(qmax,t_e(ix,it))
          ENDDO
        ENDDO             
        DO it = 1, nt
          DO ix = 1, nx
            rv(1) = t(it) - 0.5 * deltat
            zv(1) = x(ix) + 0.5 * deltax
            rv(2) = t(it) - 0.5 * deltat
            zv(2) = x(ix) - 0.5 * deltax
            rv(3) = t(it) + 0.5 * deltat
            zv(3) = x(ix) - 0.5 * deltax
            rv(4) = t(it) + 0.5 * deltat
            zv(4) = x(ix) + 0.5 * deltax
            CALL SetCol255_04(2,t_e(ix,it),qmin,qmax) 
            CALL FILCOL(255)                      
            CALL LINCOL(255)                      
            CALL PTPLOT(rv(1),zv(1),1,4,1)
          ENDDO
        ENDDO
        CALL DrawFrame

        map1x = 0.10          
        map2x = map1x + 0.55
        map2y = 0.55
        map1y = map2y - 0.20
        xlabel = 't (s) '
        ylabel = 'Z (m) '
        CALL GRTSET_TRIM(' ',' ',' ',' ',' ',      ! TITLE,REF,nVIEW,PLANE,glabel,
     .                   cxmin,cxmax,cymin,cymax,
     .                   ' ',xlabel,ylabel,        ! TABLE,XLAB,YLAB,
     .                   0,' ',0,' ',1)            ! 0,smooth,0,ANLY,1)
        CALL PSPACE(map1x,map2x,map1y,map2y)
        CALL MAP   (cxmin,cxmax,cymin,cymax)
        qmin = 0.0
        qmax = 0.0
        DO it = 1, nt
          DO ix = 1, nx
            qmax = MAX(qmax,sig(ix,it))
          ENDDO
        ENDDO             
        DO it = 1, nt
          DO ix = 1, nx
            rv(1) = t(it) - 0.5 * deltat
            zv(1) = x(ix) + 0.5 * deltax
            rv(2) = t(it) - 0.5 * deltat
            zv(2) = x(ix) - 0.5 * deltax
            rv(3) = t(it) + 0.5 * deltat
            zv(3) = x(ix) - 0.5 * deltax
            rv(4) = t(it) + 0.5 * deltat
            zv(4) = x(ix) + 0.5 * deltax
            CALL SetCol255_04(2,sig(ix,it),qmin,qmax) 
            CALL FILCOL(255)                      
            CALL LINCOL(255)                      
            CALL PTPLOT(rv,zv,1,4,1)
          ENDDO
        ENDDO
        CALL DrawFrame




        CALL FRAME
 
      ENDIF



      RETURN
 99   STOP
      END
c
c ======================================================================
c
c sburoutine: CalcCrossFieldFlux
c
c 
      SUBROUTINE CalcRadialFlux(osmtmp)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'
      INCLUDE 'slout'

      INTEGER ik,ir
      REAL    t1,t2,frac,tgrad,area,vperp,osmtmp(MAXNKS,MAXNRS)


      osmtmp = 0.0

      IF (nopriv) RETURN

      DO ir = irsep, irsep
c      DO ir = irsep, 19
        DO ik = nks(ir), ikto, -1
c        DO ik = nks(ir), 2, -1
c        DO ik = nks(ir), nks(ir)/2, -1

C          IF (zs(ik,ir).GT.zxp) EXIT

          IF (.TRUE.) THEN
c...        Temperature gradient calculated along center of cell:          

            IF     (ik.EQ.nks(ir)) THEN
              t1 = kteds(idds(ir,1))
            ELSEIF (ik.EQ.1) THEN
              t1 = kteds(idds(ir,2))
            ELSE
              frac = (ksb(ik  ,ir) - kss(ik,ir)) / 
     .               (kss(ik+1,ir) - kss(ik,ir))
              t1 = ktebs(ik,ir) * (1.0 - frac) + ktebs(ik+1,ir) * frac
            ENDIF

            frac = (ksb(ik-1,ir) - kss(ik-1,ir)) / 
     .             (kss(ik  ,ir) - kss(ik-1,ir))
            t2 = ktebs(ik-1,ir) * (1.0 - frac) + ktebs(ik,ir) * frac

            tgrad = (t1 - t2) / (kpb(ik,ir) - kpb(ik-1,ir))

            vperp = 0.71 * tgrad / 5.0 

c...        Approximation:
            area = 2.0 * PI * rs(ik,ir) * (kpb(ik,ir) - kpb(ik-1,ir)) * 
     .             eirtorfrac

            osmtmp(ik,ir) = vperp * knbs(ik,ir) * area / kvols(ik,ir)

c          WRITE(0,'(A,2I6,1P,4E10.2,0P,4X,3F10.2)') 
c     .      'FLUX:',ik,ir,vperp,tgrad,area,osmtmp(ik,ir),t1,t2,
c     .      ktebs(ik,ir)
          ELSE
            CALL ER('CalcRadialFlux','Invalid option',*99)
          ENDIF

        ENDDO
      ENDDO

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: CheckImaginaryDensity
c
c 
      SUBROUTINE CheckImaginaryDensity
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'
      INCLUDE 'slout'

      REAL    GetCs,CalcPressure

      INTEGER i1,ik,ir,id0,id6
      REAL    isat0,isat6,ion1x,rec1x,cfp1x,flx0,flx6,flx1x,
     .        p0,p6,mom1x,p1x,fact1,fact2

      DO i1 = 1, 4

        IF (i1.EQ.1) ir = 71
        IF (i1.EQ.2) ir = 18
        IF (i1.EQ.3) ir = 19
        IF (i1.EQ.4) ir = 22


        id0 = idds(ir,2)
        id6 = idds(ir,1)

        isat0 = knds(id0) * ECH * GetCs(kteds(id0),ktids(id0))
        isat6 = knds(id6) * ECH * GetCs(kteds(id6),ktids(id6))

        flx0 = -isat0 / ECH
        flx6 =  isat6 / ECH

        p0 = CalcPressure(knds(id0),kteds(id0),ktids(id0),kvds(id0))*ECH
        p6 = CalcPressure(knds(id6),kteds(id6),ktids(id6),kvds(id6))*ECH

        DO ik = nks(ir)/2, nks(ir)

          CALL CalcIntegral3(osmion,1,ik,ir,ion1x,8)
          CALL CalcIntegral3(osmrec,1,ik,ir,rec1x,8)
          CALL CalcIntegral3(osmcfp,1,ik,ir,cfp1x,8)
          flx1x = ion1x + flx0 - rec1x + cfp1x
          CALL CalcIntegral3(osmmp,1,ik,ir,mom1x,8)
          p1x = p0 + mom1x
          fact1 = (1.0*p1x)**2 / ((4.0*2.0*1.67E-27 * flx1x) * flx1x) /
     .            ((ktebs(ik,ir) + ktibs(ik,ir)) * ECH)

          CALL CalcIntegral3(osmion,ik,nks(ir),ir,ion1x,10)
          CALL CalcIntegral3(osmrec,ik,nks(ir),ir,rec1x,10)
          CALL CalcIntegral3(osmcfp,ik,nks(ir),ir,cfp1x,10)
          flx1x = -ion1x + flx6 + rec1x - cfp1x
          flx1x = -ion1x + flx6 + rec1x
          flx1x = -ion1x + flx6 
          CALL CalcIntegral3(osmmp,ik,nks(ir),ir,mom1x,10)
          p1x = p6 - mom1x
          fact2 = (1.5*p1x)**2 / ((4.0*2.0*1.67E-27 * flx1x) * flx1x) /
     .            ((ktebs(ik,ir) + ktibs(ik,ir)) * ECH)


          WRITE(6,'(A,2I6,2F10.2,1P,2X,2E10.2,0P,2F10.2,1P,E10.2,0P)') 
     .      'img:',ik,ir,fact1,fact2,flx1x,p1x,
     .      ktebs(ik,ir),ktibs(ik,ir),knbs(ik,ir)

        ENDDO
      ENDDO




      RETURN
 99   STOP
      END
c
c ======================================================================
c
c sburoutine: CalcPerpVel
c
c 
      SUBROUTINE CalcPerpVel
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'
      INCLUDE 'slout'

      REAL GetFlux

      INTEGER ik,ir
      REAL solsrc,solsnk,solden,solflx,
     .     pfzsrc,pfzsnk,pfzden,pfzflx,
     .     cfdvel(MAXNRS),cfdflx,area


c...  Find flux into PFZ from outer SOL:
      solsrc = 0.0
      solsnk = 0.0
      DO ir = irsep, irsep+4
c...    Ion source:
        DO ik = nks(ir)/2+1, nks(ir)
          solsrc = solsrc + kvols(ik,ir) * pinion(ik,ir) / eirtorfrac
        ENDDO
c...    Ion sink:
        solsnk = solsnk + ABS(GetFlux(IKHI,ir)) / eirtorfrac
      ENDDO

      solflx = 0.9 * (solsrc - solsnk) 

c...  Find average density for outer sol:
      solden = 0.0
      area = 0.0
      ir = irsep
      DO ik = ikti, nks(ir)
c...    Area:
        area = 2.0 * PI * rs(ik,ir) * (kpb(ik,ir) - kpb(ik-1,ir))
        solden = solden + (kpb(ik,ir) - kpb(ik-1,ir)) * knbs(ik,ir)
      ENDDO
      solden = solden / (kpb(nks(ir),ir) - kpb(ikti-1,ir))

      cfdvel(ir) = solflx / (solden * area)

      WRITE(0,*) 'OUT:',ir,solflx,solden,cfdvel(ir)

c... Find PFZ sink:
      pfzsrc = 0.0
      pfzsnk = 0.0
      DO ir = irtrap+1, nrs
        DO ik = 1, nks(ir)
          pfzsrc = pfzsrc + kvols(ik,ir) * pinion(ik,ir) / eirtorfrac
          pfzsnk = pfzsnk + kvols(ik,ir) * pinrec(ik,ir) / eirtorfrac
        ENDDO
        pfzsnk = pfzsnk + (ABS(GetFlux(IKLO,ir)) + 
     .                     ABS(GetFlux(IKHI,ir))) / eirtorfrac
      ENDDO

      WRITE(0,*) 'PFZ:',pfzsrc-pfzsnk,pfzsrc/pfzsnk


c...  Calculate source and sink on ring and excess flux to 
c     be transported cross-field:
      DO ir = nrs, irtrap+1, -1

        IF (ir.EQ.nrs) THEN
          cfdflx = solflx
        ELSE
          cfdflx = pfzflx
        ENDIF

c...    Find source and sink on ring:
        pfzsrc = 0.0
        pfzsnk = 0.0
        DO ik = 1, nks(ir)
          pfzsrc = pfzsrc + kvols(ik,ir) * pinion(ik,ir) / eirtorfrac
          pfzsnk = pfzsnk + kvols(ik,ir) * pinrec(ik,ir) / eirtorfrac
        ENDDO
        pfzsnk = pfzsnk + (ABS(GetFlux(IKLO,ir)) + 
     .                     ABS(GetFlux(IKHI,ir))) / eirtorfrac

        pfzflx = cfdflx + pfzsrc - pfzsnk

c...    Find representative density:
        pfzden = 0.0
        DO ik = 1, nks(ir)
          pfzden = pfzden + (kpb(ik,ir) - kpb(ik-1,ir)) * knbs(ik,ir)
        ENDDO
        pfzden = pfzden / (kpb(nks(ir),ir) - kpb(0,ir))
        pfzden = pfzden * 0.1

        cfdvel(ir) = pfzflx / pfzden        

        WRITE(0,*) 'OUT:',ir,pfzflx,pfzden,cfdvel(ir)

      ENDDO

      RETURN
99    STOP
      END
c
c ======================================================================
c
c sburoutine: CalcCrossFieldFlux
c
c 
      SUBROUTINE CalcCrossFieldFlux
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'
      INCLUDE 'slout'


      INTEGER ik,ir
      REAL    t1,t2,tgrad,area,flux,frac,sum,vperp,vperpavg,totdist,
     .        densityavg,totarea


      DO ir = irsep, irsep
c      DO ir = irsep, 19

        sum = 0.0
        totdist = 0.0
        vperpavg = 0.0
        densityavg = 0.0
        totarea = 0.0


c        DO ik = nks(ir), ikti, -1
        DO ik = nks(ir), nks(ir)/2, -1
          IF (zs(ik,ir).GT.zxp) EXIT

c...      Temperature gradient:          
          IF (ik.EQ.nks(ir)) THEN
            t1 = kteds(idds(ir,1))
          ELSE
            frac = (ksb(ik  ,ir) - kss(ik,ir)) / 
     .             (kss(ik+1,ir) - kss(ik,ir))
            t1 = ktebs(ik,ir) * (1.0 - frac) + ktebs(ik+1,ir) * frac
          ENDIF

          frac = (ksb(ik-1,ir) - kss(ik-1,ir)) / 
     .           (kss(ik  ,ir) - kss(ik-1,ir))
          t2 = ktebs(ik-1,ir) * (1.0 - frac) + ktebs(ik,ir) * frac

          tgrad = (t1 - t2) / (kpb(ik,ir) - kpb(ik-1,ir))

          vperp = 0.71 * tgrad / 5.0 

          vperpavg = vperpavg + vperp * (kpb(ik,ir) - kpb(ik-1,ir))
          densityavg = densityavg + knbs(ik,ir) * 
     .                              (kpb(ik,ir) - kpb(ik-1,ir))
          totdist = totdist + (kpb(ik,ir) - kpb(ik-1,ir))

          area = 2.0 * PI * rs(ik,ir) * (kpb(ik,ir) - kpb(ik-1,ir))

          totarea = totarea + area

          flux = vperp * area * knbs(ik,ir)

          sum = sum + flux

          WRITE(0,'(A,2I6,1P,4E10.2,0P,4X,3F10.2)') 
     .      'FLUX:',ik,ir,vperp,tgrad,area,flux,t1,t2,ktebs(ik,ir)


        ENDDO

        WRITE(0,*) 'sum:',sum,vperpavg/totdist,densityavg/totdist,
     .             totarea    ,ikti

      ENDDO

      RETURN
99    STOP
      END
c
c ======================================================================
c
c sburoutine: AnalyseStrata
c
c Calculate neutral production and destruction based on strata:
c
      SUBROUTINE AnalyseStrata
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'
      INCLUDE 'slout'

      REAL GetFlux,ATAN2C,CalcWidth,GetCs

      INTEGER ik,ir,i1,fp,numstrata,ir1,id,ik1,optflow,ikm,irdiv1,irdiv2
      REAL    ion(MAXSTRATA,MAXNRS,4),rec(MAXNRS),srcstr,tarflx,
     .        flux(2,MAXNRS),sumion(0:MAXSTRATA,4),sumstr(0:7),intion,
     .        intrec,deltar,deltaz,alpha,beta,cost,width,ionflx,flxvel,
     .        sumtarflx,sumintrec,sumintion,rdum1,rdum2,rdum3,rdum4,
     .        temin,sum1,sum2
     
      REAL*8 a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd
      CHARACTER*256 cdum1,dummy

      numstrata = eirnstrai

c      irdiv1 = 20  ! hds
c      irdiv2 = 22
c      irdiv1 = 21  ! mds
c      irdiv2 = 21
      irdiv1 = irwall - 1
      irdiv2 = irwall - 1

      fp = 99
      OPEN (UNIT=fp,FILE='outdata.dat',ACCESS='SEQUENTIAL',
     .      STATUS='REPLACE')


      sum1 = 0.0
      ir = irsep
      DO ik = 55, nks(ir)
        sum1 = sum1 + pinion(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)) 
      ENDDO
      sum2 = 0.0
      DO ik = nks(ir)/2, nks(ir)
        sum2 = sum2 + pinion(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)) 
      ENDDO
      WRITE(fp,*) 'SUM1/SUM2=',sum1,sum2,sum1/sum2
      WRITE(fp,*)



c...  Magnitude of each stratum:
      WRITE(fp,*) 
      WRITE(fp,'(A)') 'Stratum source strength:'

      WRITE(fp,*) 
      sumstr(1) = 0.0
      DO ir = irsep, irdiv1
        flux(IKLO,ir) = ABS(GetFlux(IKLO,ir)) / eirtorfrac
        sumstr(1) = sumstr(1) + flux(IKLO,ir)
        WRITE(fp,11) 'IS ring',ir,flux(IKLO,ir),sumstr(1)
      ENDDO      
      WRITE(fp,10) 'Inner SOL ',sumstr(1)

      WRITE(fp,*) 
      sumstr(2) = 0.0
      DO ir = irtrap+1, nrs
        flux(IKLO,ir) = ABS(GetFlux(IKLO,ir)) / eirtorfrac
        sumstr(2) = sumstr(2) + flux(IKLO,ir)
        WRITE(fp,11) 'IP ring',ir,flux(IKLO,ir),sumstr(2) 
      ENDDO      
      WRITE(fp,10) 'Inner PFZ ',sumstr(2)

      WRITE(fp,*) 
      sumstr(3) = 0.0
      DO ir = irsep, irdiv2
        flux(IKHI,ir) = ABS(GetFlux(IKHI,ir)) / eirtorfrac
        sumstr(3) = sumstr(3) + flux(IKHI,ir)
        WRITE(fp,11) 'OS ring',ir,flux(IKHI,ir),sumstr(3)
      ENDDO      
      WRITE(fp,10) 'Outer SOL ',sumstr(3)

      WRITE(fp,*) 
      sumstr(4) = 0.0
      DO ir = irtrap+1, nrs
        flux(IKHI,ir) = ABS(GetFlux(IKHI,ir)) / eirtorfrac
        sumstr(4) = sumstr(4) + flux(IKHI,ir)
        WRITE(fp,11) 'OP ring',ir,flux(IKHI,ir),sumstr(4) 
      ENDDO      
      WRITE(fp,10) 'Outer PFZ ',sumstr(4)

c...  Inner SOL recombination:
      WRITE(fp,*) 
      sumstr(5) = 0.0
      DO ir = irsep, irdiv1
c      DO ir = irsep, irwall-1
        rec(ir) = 0.0
        DO ik = 1, nks(ir)/2
          rec(ir) = rec(ir) + kvols(ik,ir) * pinrec(ik,ir) / eirtorfrac
        ENDDO
        sumstr(5) = sumstr(5) + rec(ir)
        WRITE(fp,11) 'SO ring',ir,rec(ir),sumstr(5)
      ENDDO
      WRITE(fp,10) 'Inner SOL rec   ',sumstr(5)

c...  Outer SOL recombination:
      WRITE(fp,*) 
      sumstr(6) = 0.0
      DO ir = irsep, irdiv2
c      DO ir = irsep, irwall-1
        rec(ir) = 0.0
        DO ik = nks(ir)/2 + 1, nks(ir)
          rec(ir) = rec(ir) + kvols(ik,ir) * pinrec(ik,ir) / eirtorfrac
        ENDDO
        sumstr(6) = sumstr(6) + rec(ir)
        WRITE(fp,11) 'SO ring',ir,rec(ir),sumstr(6)
      ENDDO
      WRITE(fp,10) 'Outer SOL rec   ',sumstr(6)

c...  PFZ recombination:
      WRITE(fp,*) 
      sumstr(7) = 0.0
      DO ir = irtrap+1, nrs
        rec(ir) = 0.0
        DO ik = 1, nks(ir)
          rec(ir) = rec(ir) + kvols(ik,ir) * pinrec(ik,ir) / eirtorfrac
        ENDDO
        sumstr(7) = sumstr(7) + rec(ir)
        WRITE(fp,11) 'PF ring',ir,rec(ir),sumstr(7)
      ENDDO
      WRITE(fp,10) 'PFZ rec   ',sumstr(7)

c...  Total neutral source:
      WRITE(fp,*) 
      sumstr(0) = 0.0
      DO i1 = 1, 7
        sumstr(0) = sumstr(0) + sumstr(i1)
      ENDDO


c...  Add puffing sources:
      WRITE(0,*)
      DO i1 = 1, eirnpuff
        WRITE(0 ,80) 'PUFF',i1,eirpuff(i1,MAXPUFF)/ECH/eirtorfrac,
     .                     NINT(eirpuff(i1,MAXPUFF-1))
        WRITE(fp,80) 'PUFF',i1,eirpuff(i1,MAXPUFF)/ECH/eirtorfrac,
     .                     NINT(eirpuff(i1,MAXPUFF-1))
        sumstr(0) = sumstr(0) + eirpuff(i1,MAXPUFF) / ECH / eirtorfrac
      ENDDO
80    FORMAT(A4,I3,1P,3X,E10.2,0P,I6)

      WRITE(fp,*)
      WRITE(fp,10) 'Tot source',sumstr(0)
      WRITE(0 ,10) 'Tot source',sumstr(0)




c...  PFZ recombination:
      WRITE(fp,*) 
      rdum3 = 0.0
      rdum4 = 0.0
      DO ir = irtrap+1, nrs
        ikm = -1
        temin = -HI
        DO ik = 1, nks(ir)
          IF (ktebs(ik,ir).GT.temin) THEN
            temin = ktebs(ik,ir)
            ikm = ik
          ENDIF
        ENDDO
c        WRITE(0,*) 'IKMID:',ikm,ir

        rdum1= 0.0
        DO ik = 1, ikm
          rdum1 = rdum1 + kvols(ik,ir) * pinrec(ik,ir) / eirtorfrac
        ENDDO
        rdum2= 0.0
        DO ik = ikm+1, nks(ir)
          rdum2 = rdum2 + kvols(ik,ir) * pinrec(ik,ir) / eirtorfrac
        ENDDO
        WRITE(fp,'(A,I6,1P,2E10.2,0P)') 'PF temp',ir,rdum1,rdum2
        rdum3 = rdum3 + rdum1
        rdum4 = rdum4 + rdum2
      ENDDO
      WRITE(fp,'(A,1P,2E10.2,0P)') 'PFZ i/o   ',rdum3,rdum4


c...  Calculate the contribution to ionisation in each region per stratum:
      WRITE(fp,*) 
      WRITE(fp,'(A)') 'Ionisation per stratum for each region:'

      CALL RZero(ion,MAXSTRATA*MAXNRS*3)

c...  Inner SOL:
      WRITE(fp,*) 
      CALL RZero(sumion(0,1),MAXSTRATA+1)
      DO ir = irsep, irdiv1
        DO ik = 1, nks(ir)/2
          sumion(0,1) = sumion(0,1) + kvols(ik,ir) * pinion(ik,ir) / 
     .                                eirtorfrac
          DO i1 = 1, numstrata
            intion = kvols(ik,ir) * pindata(ik,ir,H_ION1+i1-1) / 
     .               eirtorfrac
            ion(i1,ir,1) = ion(i1,ir,1) + intion
            sumion(i1,1) = sumion(i1,1) + intion
          ENDDO
        ENDDO
        WRITE(fp,11) 'I ring',ir,(ion(i1,ir,1),i1=1,numstrata)
      ENDDO
      WRITE(fp,10) 'Inner SOL',(sumion(i1,1),i1=1,numstrata),sumion(0,1)
10    FORMAT(A,1P,12(E10.2))
11    FORMAT(A,I3,1P,12(E10.2))

c...  Outer SOL:
      WRITE(fp,*) 
      CALL RZero(sumion(0,2),MAXSTRATA+1)
      DO ir = irsep, irdiv2
        DO ik = nks(ir)/2+1, nks(ir)
          sumion(0,2) = sumion(0,2) + kvols(ik,ir) * pinion(ik,ir) / 
     .                                eirtorfrac
          DO i1 = 1, numstrata
            intion = kvols(ik,ir) * pindata(ik,ir,H_ION1+i1-1) / 
     .               eirtorfrac
            ion(i1,ir,2) = ion(i1,ir,2) + intion
            sumion(i1,2) = sumion(i1,2) + intion
          ENDDO
        ENDDO
        WRITE(fp,11) 'O ring',ir,(ion(i1,ir,2),i1=1,numstrata)
      ENDDO
      WRITE(fp,10) 'Outer SOL',(sumion(i1,2),i1=1,numstrata),sumion(0,2)

c...  PFZ:
      WRITE(fp,*) 
      CALL RZero(sumion(0,3),MAXSTRATA+1)
      DO ir = irtrap+1, nrs
c        CALL RZero(ion(1,1,3),MAXSTRATA)
        DO ik = 1, nks(ir)
          sumion(0,3) = sumion(0,3) + kvols(ik,ir) * pinion(ik,ir) / 
     .                                eirtorfrac
          DO i1 = 1, numstrata
            intion = kvols(ik,ir) * pindata(ik,ir,H_ION1+i1-1) / 
     .               eirtorfrac
            ion(i1,ir,3) = ion(i1,ir,3) + intion
            sumion(i1,3) = sumion(i1,3) + intion
          ENDDO
        ENDDO
        WRITE(fp,11) 'P ring',ir,(ion(i1,ir,3),i1=1,numstrata)
      ENDDO
      WRITE(fp,10) 'Total PFZ',(sumion(i1,3),i1=1,numstrata),sumion(0,3)

c...  Core:
      WRITE(fp,*) 
      CALL RZero(sumion(0,4),MAXSTRATA+1)
      DO ir = 2, irsep-1
        DO ik = 1, nks(ir)
          sumion(0,4) = sumion(0,4) + kvols(ik,ir) * pinion(ik,ir) / 
     .                                eirtorfrac
          DO i1 = 1, numstrata
            intion = kvols(ik,ir) * pindata(ik,ir,H_ION1+i1-1) / 
     .               eirtorfrac
            ion(i1,ir,4) = ion(i1,ir,4) + intion
            sumion(i1,4) = sumion(i1,4) + intion
          ENDDO
        ENDDO
        WRITE(fp,11) 'C ring',ir,(ion(i1,ir,4),i1=1,numstrata)
      ENDDO
      WRITE(fp,10) 'Total COR',(sumion(i1,4),i1=1,numstrata),sumion(0,4)

      WRITE(fp,*)
      WRITE(fp,10) 'TOTAL IONISATION',sumion(0,1)+sumion(0,2)+
     .                                sumion(0,3)+sumion(0,4)



c...  Some fun sums for window-frame:
      rdum1 = 0.0
      DO ir = 2, nrs
        IF (idring(ir).EQ.-1) CYCLE
        DO ik = 1, nks(ir)
          rdum1 = rdum1 + kvols(ik,ir) * pinion(ik,ir)
        ENDDO
      ENDDO
      WRITE(0,*) 'TOTAL IONISATION:',rdum1

      rdum1 = 0.0
      DO ir = 35, 39
        DO ik = 1, nks(ir)
          rdum1 = rdum1 + kvols(ik,ir) * pinion(ik,ir)
        ENDDO
      ENDDO
      WRITE(0,*) 'WINDOW IONISATION:',rdum1

      rdum1 = 0.0
      DO ir = 2,irsep-1
        DO ik = 1, 11
          rdum1 = rdum1 + kvols(ik,ir) * pinion(ik,ir)
        ENDDO
      ENDDO
      WRITE(0,*) 'CORE IONISATION INSIDE SOURCE:',rdum1

      rdum1 = 0.0
      DO ir = 2,irsep-1
        DO ik = 12, nks(ir)-1
          rdum1 = rdum1 + kvols(ik,ir) * pinion(ik,ir)
        ENDDO
      ENDDO
      WRITE(0,*) 'CORE IONISATION OUTSIDE SOURCE:',rdum1

      rdum1 = 0.0
      DO ir = 2,irsep-1
        DO ik = 1, nks(ir)-1
          IF (zs(ik,ir).GT.0.0) CYCLE
          rdum1 = rdum1 + kvols(ik,ir) * pinion(ik,ir)
        ENDDO
      ENDDO
      WRITE(0,*) 'CORE IONISATION BELOW MIDPLANE:',rdum1

      rdum1 = 0.0
      DO ir = 2,irsep-1
        DO ik = 1, nks(ir)-1
          IF (rs(ik,ir).GT.rxp) CYCLE
          rdum1 = rdum1 + kvols(ik,ir) * pinion(ik,ir)
        ENDDO
      ENDDO
      WRITE(0,*) 'CORE IONISATION INSIDE X-POINT:',rdum1


      rdum1 = 0.0
      DO ir = irsep, 34
        rdum1 = rdum1 + flux(IKHI,ir)
      ENDDO
      WRITE(0,*) 'FLUX OUTER DIVERTOR:',rdum1

      rdum1 = 0.0
      DO ir = 35, 39
        rdum1 = rdum1 + flux(IKLO,ir) + flux(IKHI,ir)
      ENDDO
      WRITE(0,*) 'FLUX WINDOW:',rdum1








c...  Normalized ionisation as a function of ring/half-ring and region:
      WRITE(fp,*) 
      WRITE(fp,'(A)') 'Normalized ionisation:'

c...  Inner SOL:
      WRITE(fp,*) 
      DO ir = irsep, irwall-1
        srcstr = flux(IKLO,ir) + rec(ir)
        WRITE(fp,12) 'I ring',ir,(ion(i1,ir,1)/srcstr,i1=1,numstrata)
      ENDDO
      srcstr = sumstr(1) + sumstr(5)
      WRITE(fp,13) 'Inner SOL',(sumion(i1,1)/srcstr,i1=1,numstrata),
     .                         sumion(0,1)/srcstr
12    FORMAT(A,I3,12(F10.3:))
13    FORMAT(A,3X,12(F10.3:))

c...  Outer SOL:
      WRITE(fp,*) 
      DO ir = irsep, irwall-1
        srcstr = flux(IKHI,ir)
        WRITE(fp,12) 'O ring',ir,(ion(i1,ir,2)/srcstr,i1=1,numstrata)
      ENDDO
      srcstr = sumstr(3) + sumstr(6)
      WRITE(fp,13) 'Outer SOL',(sumion(i1,2)/srcstr,i1=1,numstrata),
     .                         sumion(0,2)/srcstr

c...  PFZ:
      WRITE(fp,*) 
      DO ir = irtrap+1, nrs
        srcstr = flux(IKLO,ir) + flux(IKHI,ir) + rec(ir)
        WRITE(fp,12) 'P ring',ir,(ion(i1,ir,3)/srcstr,i1=1,numstrata)
      ENDDO
      srcstr = sumstr(2) + sumstr(4) + sumstr(7)
      WRITE(fp,13) 'Total PFZ',(sumion(i1,3)/srcstr,i1=1,numstrata),
     .                         sumion(0,3)/srcstr
      




c...  Particle balance in the divertor:


c...  Check if parallel plasma flow should be over-written:
      READ(5,'(A80)') dummy
      IF (dummy(8:11).EQ.'Flow'.OR.dummy(8:11).EQ.'FLOW'.OR.
     .    dummy(8:11).EQ.'flow') THEN
        READ(dummy,*) cdum1,optflow
        WRITE(0,*) 'CALCULATING FLOW!'
        CALL CalcFlow(optflow)
      ELSE
      
        WRITE(0,*) 'SCALING BACK SOURCES IN ANALYSESTRATA'
      
        DO ir = 1, nrs
          DO ik = 1, nks(ir)       
            pinion(ik,ir) = pinion(ik,ir) / eirsrcmul
            pinrec(ik,ir) = pinrec(ik,ir) / eirsrcmul
          ENDDO
        ENDDO
      
        BACKSPACE 5
      ENDIF


c...  Inside IR < 20:
c      ir1 = 24
      ir1 = 19
 
      a1 =  0.50811D0
      a2 = -0.39848D0
      b1 =  DBLE(rxp)
      b2 =  DBLE(zxp)

      sumintion = 0.0
      sumintrec = 0.0
      sumtarflx = 0.0

c      DO ir = ir1, ir1

      DO ir = irsep, ir1

        DO ik = 1, nks(ir)/2
          id = korpg(ik,ir)
          c1 = DBLE(0.5 * (rvertp(1,id) + rvertp(2,id)))
          c2 = DBLE(0.5 * (zvertp(1,id) + zvertp(2,id)))
          d1 = DBLE(0.5 * (rvertp(3,id) + rvertp(4,id)))
          d2 = DBLE(0.5 * (zvertp(3,id) + zvertp(4,id)))
          CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
c          WRITE(0,'(2I6,8F12.4,2F14.3)') 
c     .      ir,ik,a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd
          IF (tab.GE.0.0.AND.tab.LE.1.0.AND.tcd.GE.0.0.AND.tcd.LE.1.0) 
     .      EXIT
        ENDDO
        IF (ik.EQ.nks(ir)/2+1) THEN
          CALL ER('AnalyseStrata','Inner divertor boundary not found',
     .            *99)
        ENDIF
c        ik1 = ik + 1
        ik1 = ik

c...    Neutral soure and sink in inner SOL:
        intrec = 0.0
        intion = 0.0
        DO ik = 1, ik1
          intrec = intrec + kvols(ik,ir) * pinrec(ik,ir) / eirtorfrac
          intion = intion + kvols(ik,ir) * pinion(ik,ir) / eirtorfrac
        ENDDO
c...    Neutral production at target:
        tarflx = flux(IKLO,ir)

        sumintrec = sumintrec + intrec
        sumintion = sumintion + intion
        sumtarflx = sumtarflx + tarflx

c...    Ion flux across boundary:
        deltar = SNGL(c1 - d1)
        deltaz = SNGL(c2 - d2)

          alpha  = ATAN2C(deltar,deltaz)
          deltar = SNGL(a1 - b1)
          deltaz = SNGL(a2 - b2)
          beta   = ATAN2C(deltar,deltaz) - alpha
          cost   = COS(PI / 2.0 - beta)
          IF (cost.LT.0.0) THEN
            STOP 'ERROR'
          ENDIF

c          width = CalcWidth(1,ir,CENTER,TOTAL)
          width = CalcWidth(ik,ir,CENTER,TOTAL)

        ionflx = knbs(ik1,ir) * kvhs(ik1,ir) / qt *
     .           width * 2.0 * PI * rs(ik1,ir) *
     .           eirtorfrac / kbfs(ik,ir)

          flxvel = (tarflx + intrec) - intion

          flxvel = flxvel / 
     .             (knbs(ik1,ir) * width * 2.0 * PI * rs(ik1,ir) *
     .              eirtorfrac / kbfs(ik,ir))
     
          WRITE(0,'(A,I6,3F10.4)') 'COST:',ir,-flxvel,
     .         flxvel/GetCs(ktebs(ik1,ir),ktibs(ik1,ir)),
     .         kvhs(ik1,ir)/qt

      ENDDO

      WRITE(0,*) 'INNER SOL REC:',sumintrec
      WRITE(0,*) 'INNER SOL ION:',-sumintion
      WRITE(0,*) 'INNER SOL FLX:',sumtarflx


c...  Outside IR < 21:
c      ir1 = 24
      ir1 = 20

      a1 =  0.61703
      a2 = -0.43445
      b1 =  rxp
      b2 =  zxp


      intion = 0.0
      intrec = 0.0
      tarflx = 0.0
      DO ir = irsep, ir1
        DO ik = nks(ir), nks(ir)/2, -1
          id = korpg(ik,ir)
          c1 = DBLE(0.5 * (rvertp(1,id) + rvertp(2,id)))
          c2 = DBLE(0.5 * (zvertp(1,id) + zvertp(2,id)))
          d1 = DBLE(0.5 * (rvertp(3,id) + rvertp(4,id)))
          d2 = DBLE(0.5 * (zvertp(3,id) + zvertp(4,id)))
          CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
c          WRITE(0,'(2I6,8F12.4,2F14.3)') 
c     .      ir,ik,a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd
          IF (tab.GE.0.0.AND.tab.LE.1.0.AND.tcd.GE.0.0.AND.tcd.LE.1.0) 
     .      EXIT
        ENDDO
        IF (ik.EQ.nks(ir)/2-1) THEN
          CALL ER('AnalyseStrata','Outer divertor boundary not found',
     .            *99)
        ENDIF
        ik1 = ik - 1

c...    Ion flux across boundary:




c...    Neutral soure and sink in inner SOL:
        DO ik = nks(ir), ik1, -1
          intrec = intrec + kvols(ik,ir) * pinrec(ik,ir) / eirtorfrac
          intion = intion + kvols(ik,ir) * pinion(ik,ir) / eirtorfrac
        ENDDO
c...    Neutral production at target:
        tarflx = tarflx + flux(IKHI,ir)
      ENDDO

      WRITE(0,*) 'OUTER SOL REC:',intrec
      WRITE(0,*) 'OUTER SOL ION:',-intion
      WRITE(0,*) 'OUTER SOL FLX:',tarflx





c...  PFZ (calculated above):


      WRITE(0,*) '      PFZ REC:',sumstr(7)
      WRITE(0,*) '      PFZ ION:',-sumion(0,3)
      WRITE(0,*) '      PFZ FLX:',(sumstr(2) + sumstr(4))

c... Puff:
c...  Add puffing sources:
      DO i1 = 1, eirnpuff
        IF (NINT(eirpuff(i1,MAXPUFF-1)).EQ.1) THEN
          WRITE(0,*) '         PUFF:',eirpuff(i1,MAXPUFF)/ECH,i1
        ENDIF
      ENDDO


      CLOSE(fp)



      RETURN
99    STOP
      END
c
c ======================================================================
c
c sburoutine: CalcFlow
c
c Calculate the flow pattern from unscaled sources and sinks using the 
c calculated background plasma.
c
c
c
      SUBROUTINE CalcFlow(optflow)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'
      INCLUDE 'slout'

      INTEGER optflow,ik,ik1,ir
      REAL    ionsrc,tarflx,tarflx1,tarflx2,recsrc,cfpsrc,deltas,parsrc,
     .        pfzdef,totvol,dummy(MAXNKS),tdummy(MAXNDS),pfzsrc,tmpsum,
     .        scale,discfp(MAXNKS),
     .        dftsrc,dist(MAXNKS,MAXNRS),sum,rdum1


      REAL GetFlux


      WRITE(char(26),'(A,I2)') 'Flow calculation option =',optflow

      IF (.FALSE..AND.optflow.EQ.3) THEN
c...    Limit the amount of over/underionisation at the low index target:

        DO ir = irsep, irwall-1

c          tarflx = ABS(knds(idds(ir,1)) * kvds(idds(ir,1)))
c          CALL CalcIntegral4(pinrec,nks(ir)/2+1,nks(ir),ir,recsrc,2)

          tarflx = ABS(knds(idds(ir,2)) * kvds(idds(ir,2)))
          CALL CalcIntegral4(pinrec,nks(ir)/2+1,nks(ir),ir,recsrc,2)


          CALL CalcIntegral3(pinion,1,nks(ir)/2,ir,ionsrc,2)

          scale = ionsrc / (recsrc + tarflx)

          WRITE(0,*) 'scale:',ir,scale
      
          IF     (scale.GT.1.20) THEN
            scale = scale / 1.20      
          ELSEIF (scale.LT.0.80) THEN
            scale = scale / 0.80      
          ELSE  
            scale = 1.0
          ENDIF
       
          DO ik = 1, nks(ir)/2
            pinion(ik,ir) = pinion(ik,ir) / scale
          ENDDO         

        ENDDO

      ENDIF

      IF (optflow.EQ.2.OR.optflow.EQ.3) THEN      
c...    Calculate ion deficit for the PFZ:
        tarflx = 0.0
        recsrc = 0.0
        ionsrc = 0.0
        DO ir = 2, irsep-1
c...      Calculate total volume recombination and ionisation 
c         sources (on the standard grid):
          DO ik = 1, nks(ir)-1
            recsrc = recsrc + pinrec(ik,ir) * kvols(ik,ir)
            ionsrc = ionsrc + pinion(ik,ir) * kvols(ik,ir)
          ENDDO
        ENDDO

        pfzdef = ionsrc - recsrc - tarflx

        tmpsum = pfzdef

c        WRITE(0,*) 'TARFLX=',tarflx
c        WRITE(0,*) 'RECSRC=',recsrc
c        WRITE(0,*) 'IONSRC=',ionsrc
c        WRITE(0,*) 'CORE EXCESS:',pfzdef,tmpsum

        tarflx = 0.0
        recsrc = 0.0
        ionsrc = 0.0
        DO ir = irsep, irwall-1
c...      Calculate the total flux to the targets:     
          tarflx = tarflx + ABS(GetFlux(IKLO,ir))
c...      Calculate total volume recombination and ionisation 
c         sources (on the standard grid):
          DO ik = 1, nks(ir)/2
            recsrc = recsrc + pinrec(ik,ir) * kvols(ik,ir)
            ionsrc = ionsrc + pinion(ik,ir) * kvols(ik,ir)
          ENDDO
        ENDDO

        pfzdef = ionsrc - recsrc - tarflx

        tmpsum = tmpsum + pfzdef

c        WRITE(0,*) 'TARFLX DIVERTOR ONLY!=',tarflx
c        WRITE(0,*) 'RECSRC DIVERTOR ONLY!=',recsrc
c        WRITE(0,*) 'IONSRC DIVERTOR ONLY!=',ionsrc
c        WRITE(0,*) 'INNER EXCESS:',pfzdef,tmpsum

        tarflx = 0.0
        recsrc = 0.0
        ionsrc = 0.0
        DO ir = irsep, irwall-1
c...      Calculate the total flux to the targets:     
          tarflx = tarflx + ABS(GetFlux(IKHI,ir))
c...      Calculate total volume recombination and ionisation 
c         sources (on the standard grid):
          DO ik = nks(ir)/2+1, nks(ir)
            recsrc = recsrc + pinrec(ik,ir) * kvols(ik,ir)
            ionsrc = ionsrc + pinion(ik,ir) * kvols(ik,ir)
          ENDDO
        ENDDO

        pfzdef = ionsrc - recsrc - tarflx

        tmpsum = tmpsum + pfzdef

c        WRITE(0,*) 'TARFLX=',tarflx
c        WRITE(0,*) 'RECSRC=',recsrc
c        WRITE(0,*) 'IONSRC=',ionsrc
c        WRITE(0,*) 'OUTER EXCESS:',pfzdef,tmpsum

        tarflx = 0.0
        recsrc = 0.0
        ionsrc = 0.0
        DO ir = irsep, irsep
c...      Calculate the total flux to the targets:     
c          tarflx = tarflx + ABS(GetFlux(IKLO,ir))
          tarflx = tarflx + ABS(GetFlux(IKHI,ir))
c...      Calculate total volume recombination and ionisation 
c         sources (on the standard grid):
          DO ik = nks(ir)/2+1, nks(ir)
            recsrc = recsrc + pinrec(ik,ir) * kvols(ik,ir)
            ionsrc = ionsrc + pinion(ik,ir) * kvols(ik,ir)
          ENDDO
        ENDDO

        pfzdef = ionsrc - recsrc - tarflx

c        WRITE(0,*) 'TARFLX=',tarflx
c        WRITE(0,*) 'RECSRC=',recsrc
c        WRITE(0,*) 'IONSRC=',ionsrc
c        WRITE(0,*) 'SEPARATRIX EXCESS:',pfzdef,tmpsum

        tarflx = 0.0
        recsrc = 0.0
        ionsrc = 0.0
        DO ir = irtrap+1, nrs
c...      Calculate the total flux to the targets:     
          tarflx = tarflx + ABS(GetFlux(IKLO,ir))
          tarflx = tarflx + ABS(GetFlux(IKHI,ir))
c...      Calculate total volume recombination and ionisation 
c         sources (on the standard grid):
          DO ik = 1, nks(ir)
            recsrc = recsrc + pinrec(ik,ir) * kvols(ik,ir)
            ionsrc = ionsrc + pinion(ik,ir) * kvols(ik,ir)
          ENDDO
        ENDDO

        pfzdef = ionsrc - recsrc - tarflx - 3.0E+20

        WRITE(0,*) 'TARFLX=',tarflx
        WRITE(0,*) 'RECSRC=',recsrc
        WRITE(0,*) 'IONSRC=',ionsrc
        WRITE(0,*) 'PFZ EXCESS:',pfzdef,tmpsum

      ELSEIF (optflow.EQ.1.OR.optflow.EQ.4.OR.optflow.EQ.5.OR.
     .        optflow.EQ.6) THEN

c...   Scale back everything:

       WRITE(0,*) 'SCALING BACK SOURCES IN CALCFLOW'

       DO ir = 1, nrs
         DO ik = 1, nks(ir)       
           pinion(ik,ir) = pinion(ik,ir) / eirsrcmul
           pinrec(ik,ir) = pinrec(ik,ir) / eirsrcmul
         ENDDO
       ENDDO
      ELSE
        CALL ER('CalcFlow','Invalid OPTFLOW value',*99)
      ENDIF

      DO ir = irsep, irwall-1
c...  

c...    Magnitude of net cross-field flux (uniform distribution
c       over the entire ring):

        tarflx1 = ABS(knds(idds(ir,2)) * kvds(idds(ir,2)))
        tarflx2 = ABS(knds(idds(ir,1)) * kvds(idds(ir,1)))
        CALL CalcIntegral4(pinrec,1,nks(ir),ir,recsrc,2)
        CALL CalcIntegral3(pinion,1,nks(ir),ir,ionsrc,2)
        CALL CalcIntegral3(osmcfpflx(1,1,1),1,nks(ir),ir,dftsrc,2)

        CALL RZero(osmcfp(1,ir),MAXNKS)

        IF (optflow.EQ.3) THEN
c...      Assign cross-field flux between the target and the x-point that prevents
c         flow reversal:

          ik1 = 0
          DO ik = nks(ir)/2, nks(ir)
            IF (zs(ik,ir).GT.zxp) ik1 = ik + 1
          ENDDO
          IF (ik1.EQ.0) STOP 'IK TROUBLE'

          totvol = 0.0
          DO ik = ik1, nks(ir)
            totvol = totvol + kvols(ik,ir)
          ENDDO

          pfzsrc = pfzdef / eirtorfrac / totvol 

          tarflx2 = ABS(knds(idds(ir,1)) * kvds(idds(ir,1)))
          CALL CalcIntegral4(pinrec,nks(ir)/2,nks(ir),ir,recsrc,4)
          CALL CalcIntegral3(pinion,nks(ir)/2,nks(ir),ir,ionsrc,4)
          cfpsrc = -(ionsrc - recsrc - tarflx2)
          IF (cfpsrc.LT.0.0) THEN
            WRITE(0,*) 'FLOW REVERSAL COMPENSATION:',ir,cfpsrc/pfzsrc
            WRITE(6,*) 'FLOW REVERSAL COMPENSATION:',ir,cfpsrc/pfzsrc
            WRITE(char(ir-irsep+11),'(I2,F10.4)') ir,cfpsrc/pfzsrc

            pfzsrc = cfpsrc
          ELSE
            pfzsrc = 0.0
          ENDIF

          tarflx1 = ABS(knds(idds(ir,2)) * kvds(idds(ir,2)))
          CALL CalcIntegral4(pinrec,1,nks(ir)/2,ir,recsrc,4)
          CALL CalcIntegral3(pinion,1,nks(ir)/2,ir,ionsrc,4)

c...      Calculate the distribution of the cross-field particle source:
          cfpsrc = -(ionsrc - recsrc - tarflx1)
          CALL GetDist(ikbound(ir,IKLO),nks(ir)/2,ir,dummy,tdummy,
     .                 discfp,1,.TRUE.,'CFP')
          DO ik = 1, nks(ir)/2
            osmcfp(ik,ir) = 0.0
          ENDDO
          DO ik = ikbound(ir,IKLO), nks(ir)/2
            osmcfp(ik,ir) = discfp(ik) * cfpsrc 
          ENDDO

c...      Check:
          tarflx1 = ABS(knds(idds(ir,2)) * kvds(idds(ir,2)))
          tarflx2 = ABS(knds(idds(ir,1)) * kvds(idds(ir,1)))
          CALL CalcIntegral4(pinrec,1,nks(ir)/2,ir,recsrc,2)
          CALL CalcIntegral3(pinion,1,nks(ir)/2,ir,ionsrc,2)
          CALL CalcIntegral3(osmcfp,1,nks(ir)/2,ir,cfpsrc,4)          
c          WRITE(0,*) 'CRAP 1:',ir, ionsrc - recsrc - tarflx1,
c     .                cfpsrc

          tarflx2 = ABS(knds(idds(ir,1)) * kvds(idds(ir,1)))
          CALL CalcIntegral4(pinrec,nks(ir)/2,nks(ir),ir,recsrc,4)
          CALL CalcIntegral3(pinion,nks(ir)/2,nks(ir),ir,ionsrc,4)

c...      Calculate the distribution of the cross-field particle source:
          cfpsrc = -(ionsrc - recsrc - tarflx2 + pfzsrc)
          CALL GetDist(nks(ir)/2+1,nks(ir),ir,dummy,tdummy,
     .                 discfp,1,.TRUE.,'CFP')
          DO ik = nks(ir)/2+1, nks(ir)
            osmcfp(ik,ir) = discfp(ik) * cfpsrc 
          ENDDO

c...      Add on cf transport of particles to PFZ:
          cfpsrc = pfzsrc / (ksmaxs(ir) - ksb(ik1-1,ir))
          DO ik = ik1, nks(ir)
            osmcfp(ik,ir) = osmcfp(ik,ir) + cfpsrc
          ENDDO

c...      Check:
          tarflx1 = ABS(knds(idds(ir,2)) * kvds(idds(ir,2)))
          tarflx2 = ABS(knds(idds(ir,1)) * kvds(idds(ir,1)))
          CALL CalcIntegral4(pinrec,1,nks(ir),ir,recsrc,2)
          CALL CalcIntegral3(pinion,1,nks(ir),ir,ionsrc,2)
          CALL CalcIntegral3(osmcfp,1,nks(ir),ir,cfpsrc,2)          
c          WRITE(0,*) 'CRAP A:',ir,ionsrc - recsrc - tarflx1 - tarflx2,
c     .                cfpsrc

        ELSEIF (optflow.EQ.2) THEN
c...      Assign cross-field flux into PFZ:

          ik1 = 0
          DO ik = nks(ir)/2, nks(ir)
            IF (zs(ik,ir).GT.zxp) ik1 = ik + 1
          ENDDO
          IF (ik1.EQ.0) STOP 'IK TROUBLE'

          totvol = 0.0
          DO ik = ik1, nks(ir)
            totvol = totvol + kvols(ik,ir)
          ENDDO

          pfzsrc = pfzdef / eirtorfrac / totvol 

          IF     (ir.EQ.14) THEN
            pfzsrc = pfzsrc * 0.35
          ELSEIF (ir.EQ.15) THEN
            pfzsrc = pfzsrc * 0.10
          ELSEIF (ir.EQ.16) THEN
            pfzsrc = pfzsrc * 0.02
          ELSE
            pfzsrc = 0.0
          ENDIF

c...RUN A RANGE OF DISTRIBUTIONS TO SEE WHICH GIVES A GOOD MATCH TO THE VSP DATA?
c...BASE THE DISTRIBUTION OF THE CF FLUX ON T GRADIENT?  OTHER DEPENDENCIES?  HOW SENSTIVIE IS THIS?  SCALE INSIDE TO LIMIT FLOW REVERSAL THERE SO THAT FLOW AT THE VSP LOCATION IS NOT DOMINATED BY THE INNER OVERIONISTAION?
c...THIS MASSIVE FLOW REVERSAL FROM THE INSIDE...? FIX?  (CF DISTRIBUTIONS? MOVE PEAK NEXT TO CORE?  REDUCE PEAK HEIGHT?)
c...DO NEED TO ITERATE LONGER?  I SAW THAT THE "REDUCTION" BIT WAS NOT STABLE

c...      Calculate the distribution of the cross-field particle source:
          cfpsrc = -(ionsrc - recsrc - tarflx1 - tarflx2 + pfzsrc)
          CALL GetDist(ikbound(ir,IKLO),nks(ir),ir,dummy,tdummy,
c     .                 osmcfp(1,ir),17,.TRUE.,'CFP')
     .                 osmcfp(1,ir),1,.TRUE.,'CFP')
          DO ik = ikbound(ir,IKLO), nks(ir)
            osmcfp(ik,ir) = osmcfp(ik,ir) * cfpsrc 
          ENDDO

c...      Add on cf transport of particles to PFZ:
          cfpsrc = pfzsrc / (ksmaxs(ir) - ksb(ik1-1,ir))
          DO ik = ik1, nks(ir)
            osmcfp(ik,ir) = osmcfp(ik,ir) + cfpsrc
          ENDDO

c...      Check:
          CALL CalcIntegral3(osmcfp,1,nks(ir),ir,cfpsrc,2)          
c          WRITE(0,*) 'CRAP B:',ionsrc - recsrc - tarflx1 - tarflx2,
c     .                cfpsrc

        ELSEIF (optflow.EQ.4.OR.optflow.EQ.5) THEN

          CALL CalcIntegral4(pinrec,nks(ir)/2,nks(ir),ir,recsrc,2)
          CALL CalcIntegral3(pinion,nks(ir)/2,nks(ir),ir,ionsrc,2)

          cfpsrc = -(ionsrc - recsrc - tarflx2) 

c...      Calculate the distribution of the cross-field particle source:
c          CALL GetDist(nks(ir)/2,nks(ir),ir,dummy,tdummy,
c     .                 osmcfp(1,ir),9,.TRUE.,'CFP')
          CALL GetDist(nks(ir)/2,nks(ir),ir,dummy,tdummy,
     .                 osmcfp(1,ir),1,.TRUE.,'CFP')

c          CALL GetDist(ikbound(ir,IKLO),nks(ir),ir,dummy,tdummy,
c     .                 osmcfp(1,ir),1,.TRUE.,'CFP')

          DO ik = nks(ir)/2, nks(ir)
            osmcfp(ik,ir) = osmcfp(ik,ir) * cfpsrc 
          ENDDO


        ELSEIF (optflow.EQ.6) THEN
c...
          DO ik = 1, nks(ir)
            dist(ik,ir) = 0.0
            osmcfp(ik,ir) = 0.0
          ENDDO

          sum = 0.0
          DO ik = 1, nks(ir)
            dist(ik,ir) = 1.0 
c            sum = sum + dist(ik,ir)
          ENDDO
  
          DO ik = 1, nks(ir)
            deltas = ksb(ik,ir) - ksb(ik-1,ir)
c            dist(ik,ir) = dist(ik,ir) / sum / deltas
            sum = sum + deltas          
          ENDDO

          DO ik = 1, nks(ir)
            dist(ik,ir) = dist(ik,ir) / sum
          ENDDO

          CALL CalcIntegral3(dist,1,nks(ir),ir,rdum1,8)
c          WRITE(0,*) 'RDUM!:',rdum1

          dftsrc = 0.0
          cfpsrc = -(ionsrc - recsrc - dftsrc + tarflx1 - tarflx2) 

          DO ik = 1, nks(ir)
            osmcfp(ik,ir) = dist(ik,ir) * cfpsrc
          ENDDO

          CALL CalcIntegral3(osmcfp,1,nks(ir),ir,rdum1,8)
          WRITE(0,*) 'RDUM!:',rdum1,cfpsrc

c          STOP 'sdfsdf'


        ELSE
          cfpsrc = -(ionsrc - recsrc - tarflx1 - tarflx2) 

c...      Calculate the distribution of the cross-field particle source:
c          CALL GetDist(ikbound(ir,IKLO),nks(ir),ir,dummy,tdummy,
cc     .                 osmcfp(1,ir),17,.TRUE.,'CFP')
c     .                 osmcfp(1,ir),1,.TRUE.,'CFP')
     .                 
c          WRITE(0,'(A,I6,1P,4E10.2,0P)') 
c     .      'FLOW:',ir,ionsrc,tarflx1,tarflx2,cfpsrc
c          WRITE(0,*) 'IKBOUND:',ikbound(ir,IKLO)
          DO ik = 1, nks(ir)
            osmcfp(ik,ir) = 1.0
          ENDDO
       
          DO ik = ikbound(ir,IKLO), nks(ir)
            osmcfp(ik,ir) = osmcfp(ik,ir) * cfpsrc 
          ENDDO
        ENDIF      

c...    Inner SOL:
        DO ik = 1, nks(ir)/2
          CALL CalcIntegral4(pinrec,1,ik,ir,recsrc,2)
          CALL CalcIntegral3(pinion,1,ik,ir,ionsrc,2)
          CALL CalcIntegral3(osmcfpflx(1,1,1),1,ik,ir,dftsrc,2)
          dftsrc = 0.0
          CALL CalcIntegral3(osmcfp,1,ik,ir,cfpsrc,2)
          parsrc = -(tarflx1 + recsrc - ionsrc - dftsrc - cfpsrc)
          kvhs(ik,ir) = parsrc / knbs(ik,ir) * qt
        ENDDO

c...    Outer SOL:
        DO ik = nks(ir)/2, nks(ir)
          CALL CalcIntegral4(pinrec,ik,nks(ir),ir,recsrc,2)
          CALL CalcIntegral3(pinion,ik,nks(ir),ir,ionsrc,2)
          CALL CalcIntegral3(osmcfpflx(1,1,1),ik,nks(ir),ir,dftsrc,2)
          dftsrc = 0.0
          CALL CalcIntegral3(osmcfp,ik,nks(ir),ir,cfpsrc,2)
          parsrc = tarflx2 + recsrc - ionsrc - dftsrc - cfpsrc
          kvhs(ik,ir) = parsrc / knbs(ik,ir) * qt
        ENDDO
 
      ENDDO

      RETURN
99    STOP
      END
c
c ======================================================================
c
c
      SUBROUTINE Development(iopt,nizs2,cizsc2,crmi2,cion2,absfac2,
     .                       title)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'ppplas'
      INCLUDE 'slout'


      INTEGER   iopt,nizs2,ik,ir,i1,cizsc2,cion2
      REAL      array(MAXNKS,MAXNRS),crmi2,absfac2
      CHARACTER title*(*)

      IF     (iopt.EQ.2) THEN
        RETURN
      ELSEIF (iopt.EQ.3) THEN
        WRITE(0,*) 'OUT TERMINATED BY PLOT 999, OPTION 3'
        WRITE(6,*) 'OUT TERMINATED BY PLOT 999, OPTION 3'
        STOP
      ELSEIF (iopt.EQ.4) THEN
        CALL CheckImaginaryDensity
        RETURN
      ELSEIF (iopt.EQ.5) THEN
        CALL CalcPotential
        RETURN
      ELSEIF (iopt.EQ.6) THEN
        CALL DTSanalysis
c        CALL DTSanalysis(MAXGXS,MAXNGS)
        RETURN
      ELSEIF (iopt.EQ.7) THEN
        CALL CoreProfileAnalysis(nizs2,cizsc2,crmi2,cion2,absfac2)
        RETURN
      ELSEIF (iopt.EQ.8) THEN
        CALL OutputDivertorProfiles
        RETURN
      ELSEIF (iopt.EQ.9) THEN
        CALL CoreRadialParticleTransport
        RETURN
      ELSEIF (iopt.EQ.10) THEN
        WRITE(0,*) 'DUMING DATA TO IDL'
        CALL DumpDataToIDL
        RETURN
      ELSEIF (iopt.EQ.11) THEN
        WRITE(0,*) 'FLUX SURFACES FOR SPENCER'
        CALL FluxSurfacesInTheSOL
        RETURN
      ELSEIF (iopt.EQ.12) THEN
        CALL GenerateOSMDataFiles
        RETURN
      ELSEIF (iopt.EQ.13) THEN
        CALL calc_wallprad(nizs2)
        RETURN
      ELSEIF (iopt.EQ.14) THEN
        CALL GenerateDIVIMPDataFiles(nizs2,cizsc2,crmi2,cion2,
     .                               absfac2,title)
        RETURN
      ELSEIF (iopt.EQ.15) THEN
        CALL GenerateEIRENEDataFiles
        RETURN
      ELSEIF (iopt.EQ.16) THEN
        CALL ExportTetrahedrons('tetrahedrons.raw')
        RETURN
      ENDIF



      RETURN
     


      IF (nrmindex.GT.0) THEN

        WRITE(6,*) 
        WRITE(6,*) 'Calculated norms:'
        WRITE(6,*) 

        DO i1 = 1, nrmindex
        WRITE(6,'(I6,F12.6,I6,5X,A)') 
     .    nrmnum(i1),nrmvalue(i1),nrmstep(i1),
     .    nrmcomment(i1)(1:LEN_TRIM(nrmcomment(i1)))
        ENDDO

      ENDIF

      CALL AnalyseStrata

      CALL CalcCrossFieldFlux

      CALL CalcPerpVel

      CALL GREND
      STOP 'HALTING OUT FROM PLOT 999'

      CALL ProbePath2

      CALL AnalyseSolution(6)

c      CALL LoadEIRENEAtomicData

c      CALL CalcRadiatedPower(array,1)

      RETURN
99    STOP
      END
c
c ======================================================================
c
c
c
c
c
c
c
c
c
c
c
c
c =========
c
c subroutine: CalculateRho
c
c
      SUBROUTINE CalculateRho

      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      INTEGER GetModel
      INTEGER ik,ir,ir1,iki,iko,id1,id2,id3,ii,id,in,midnks,i1,
     .        ikto3,ikti3

      REAL rhozero

      DOUBLE PRECISION a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd


      WRITE(0,*) 'CALCULATING RHO IN OUT'

      a1 = r0
      a2 = 0.0
      b1 = r0 + 100.0
      b2 = 0.0

      DO ir = 2, irwall-1
        DO ik = 1, nks(ir)
          id = korpg(ik,ir)

          c1 = rvertp(1,id)
          c2 = zvertp(1,id)
          d1 = rvertp(4,id)
          d2 = zvertp(4,id)

          CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)

          IF (tab.GE.0.0.AND.tab.LE.1.0.AND.tcd.GE.0.0.AND.tcd.LE.1.0)
     .      rho(ir,IN14) = r0 + tab * 100.0

          c1 = 0.5 * (rvertp(1,id) + rvertp(2,id))
          c2 = 0.5 * (zvertp(1,id) + zvertp(2,id))
          d1 = 0.5 * (rvertp(3,id) + rvertp(4,id))
          d2 = 0.5 * (zvertp(3,id) + zvertp(4,id))

          CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)

          IF (tab.GE.0.0.AND.tab.LE.1.0.AND.tcd.GE.0.0.AND.tcd.LE.1.0)
     .      rho(ir,CELL1) = r0 + tab * 100.0

          c1 = rvertp(2,id)
          c2 = zvertp(2,id)
          d1 = rvertp(3,id)
          d2 = zvertp(3,id)

          CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)

          IF (tab.GE.0.0.AND.tab.LE.1.0.AND.tcd.GE.0.0.AND.tcd.LE.1.0)
     .      rho(ir,OUT23) = r0 + tab * 100.0
        ENDDO
      ENDDO

      rhozero = rho(irsep,IN14)

      DO ir = 2, irwall-1
        rho(ir,IN14)  = rho(ir,IN14)  - rhozero
        rho(ir,CELL1)  = rho(ir,CELL1)  - rhozero
        rho(ir,OUT23) = rho(ir,OUT23) - rhozero
      ENDDO

c This sucks... it is good to get rho from grid.. what if not a CMOD
c grid being used...?
      rho(nrs,IN14)  = 2 * rho(irsep,IN14) - rho(irsep,OUT23)
      rho(nrs,OUT23) = rho(irsep,IN14)
      rho(nrs,CELL1)  = 0.5 * (rho(nrs,IN14) + rho(nrs,OUT23))

      DO ir = nrs-1, irtrap+1, -1
        rho(ir,OUT23) = rho(ir+1,IN14)
        rho(ir,IN14)  = rho(ir+1,IN14) - 
     .                  (rho(ir+1,OUT23) - rho(ir+1,IN14))
c        rho(ir,IN14)  = 2 * rho(ir+1,IN14) - rho(ir,OUT23)
        rho(ir,CELL1)  = 0.5 * (rho(ir,IN14) + rho(ir,OUT23))
      ENDDO


      DO ir = 1, nrs
        WRITE(6,'(A,I4,1P,3E15.7)')
     .    'IR RHO = ',ir,(rho(ir,in),in=1,3)
      ENDDO

      RETURN
99    STOP 
      END






c
c ======================================================================
c
c subroutine: ProbePath2
c
      SUBROUTINE ProbePath2
      IMPLICIT   none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      INTEGER mark,i1,ik,ir,ikprb,id
      REAl    s,z
      REAL*8  a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd,tmin,tprb
c
c     Assign starting point point and ending point for probe trajectory
c     and assume that the probe moves in a straight line:
c
      a1 = 1.9400
      a2 = 0.5831 + 1.60
      
      b1 = 1.9400
      b2 = 0.9408 + 1.60

      DO ir = 2, irsep-1
        tmin  = HI
        ikprb = NULL
        DO ik = 1, nks(ir)-1
          id = korpg(ik,ir)
          DO i1 = 1, 2
            IF (i1.EQ.1) THEN
              c1 = 0.5 * (rvertp(1,id) + rvertp(2,id))
              c2 = 0.5 * (zvertp(1,id) + zvertp(2,id))
              d1 = rs(ik,ir)
              d2 = zs(ik,ir)
            ELSE
              c1 = rs(ik,ir)
              c2 = zs(ik,ir)
              d1 = 0.5 * (rvertp(3,id) + rvertp(4,id))
              d2 = 0.5 * (zvertp(3,id) + zvertp(4,id))
            ENDIF

            CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)

c            WRITE(6,'(A,3I4,1P,2E15.7)') 'path = ',ik,i1,ir,tab,tcd

c...        Assumes probe trajectory only cuts the grid once:
            IF (tab.GT.0.0.AND.tab.LT.1.0.AND.
     .          tcd.GT.0.0.AND.tcd.LT.1.0) THEN
              ikprb = ik
              tprb  = tcd
              tmin  = ABS(tab)
              IF (i1.EQ.1) THEN
                s = ksb(ik-1,ir) + tcd * (kss(ik,ir) - ksb(ik-1,ir))
                z = c2 + tcd * (d2 - c2)
              ELSE
                s = kss(ik  ,ir) + tcd * (ksb(ik,ir) - kss(ik  ,ir))
                z = c2 + tcd * (d2 - c2)
              ENDIF
              WRITE(6,'(A,3I4,1P,2E15.7,0P,2F10.4)') 'path = ',
     .          ik,i1,ir,tab,tcd,s,z
            ENDIF

          ENDDO
        ENDDO
      ENDDO

      RETURN
99    STOP
      END
c
c
c

      SUBROUTINE slCheckLOS(xb,wb,theta,dist,ind,x3,y3)

      IMPLICIT none

      REAL    xb(2),wb(2),dist,theta

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'
      INCLUDE 'slout'

      REAL    LONGWAY
      INTEGER NUMWPT
      PARAMETER (LONGWAY=1.0E+3, NUMWPT=20)

      INTEGER in1,in2,ncross,idum1,icross(NUMWPT),ind
      REAL*8  ar,az,br,bz,cr,cz,er,ez,tab,tcd
      REAL    tcross(NUMWPT),rdum1,x3,y3
c
c     Find points where LOS intersects the wall:
c
      dist   = HI
      ncross = 0
      ind    = 0
      x3     = 0.0
      y3     = 0.0

      ar = xb(1)
      az = xb(2)
      br = wb(1) * LONGWAY
      bz = wb(2) * LONGWAY

      DO in1 = 1, wallpts

        cr = wallpt(in1,20)
        cz = wallpt(in1,21)
        er = wallpt(in1,22)
        ez = wallpt(in1,23)

        CALL CalcInter(ar,az,br,bz,cr,cz,er,ez,tab,tcd)

        IF (tab.GT.0.0.AND.tab.LT.1.0.AND.
     .      tcd.GT.0.0.AND.tcd.LT.1.0) THEN
          ncross         = ncross + 1
          tcross(ncross) = tab
          icross(ncross) = in1
        ENDIF
      ENDDO

      IF (ncross.GT.NULL) THEN
        IF (ncross.GT.NUMWPT) STOP 'ERROR CHECKLOS: Array out of bounds'

10      DO in1 = 1, ncross-1
          IF (tcross(in1).GT.tcross(in1+1)) THEN
            rdum1          = tcross (in1)
            idum1          = icross (in1)
            tcross (in1)   = tcross (in1+1)
            icross (in1)   = icross (in1+1)
            tcross (in1+1) = rdum1
            icross (in1+1) = idum1
            GOTO 10
          ENDIF
        ENDDO


        IF     (slopt5.EQ.2.AND.ncross.GE.4) THEN
c...      See notes in PLOT980 routine where
c         slopt5 set equal to 2:
c          IF (ncross.LT.3) STOP 'NO WORK'
          in1 = 4
        ELSEIF (MOD(ncross,2).EQ.NULL) THEN
          in1 = 2
        ELSE
          in1 = 1
        ENDIF

        br = ar + tcross(in1) * (br - ar)
        bz = az + tcross(in1) * (bz - az)

        dist = SQRT((ar - br)**2.0 + (az - bz)**2.0)
        ind  = icross(in1)
        x3   = SNGL(br)
        y3   = SNGL(bz)
      ENDIF

      dist = dist * 1.001

c     WRITE(6,'(A,4F8.3,2F8.4,F6.1,2X,I4,F12.4,0P,I4,1P,10(E11.3:))')
c    . 'CHECKLOS: ',ar,az,br,bz,wb(1),wb(2),theta*180.0/PI,
c    . in1,dist,
c    . ncross,
c    . (tcross(in1),in1=1,ncross)

      RETURN
      END

      SUBROUTINE slLOSINT (TVALS,TOUTS,TWIDS,NUMTHE,ROBS,ZOBS,AVPTS,VS,
     .                     DUM1)
      IMPLICIT NONE
C
C  *********************************************************************
C  *                                         WRITE(0,*) tvals2(1:numth2,1)                                    *
C  *  LOSINT:  INTEGRATE THE VARIABLE VS ALONG A FAN OF SIGHT LINES    *
C  *           FROM A COMMON OBSERVATION POINT (ROBS,ZOBS).  SPATIAL   *
C  *           RESOLUTION IS SIMULATED BY AN AVPTS-POINT AVERAGE OF    *
C  *           A SET OF CHORDS SPANNING THE INTERVAL TWIDS.  AT THE    *
C  *           MOMENT THERE ARE NO WEIGHTS ON THIS AVERAGING WHICH     *
C  *           IMPLIES THAT WE ARE ASSUMING A RECTANGULAR, RATHER      *
C  *           THAN A CIRCULAR, VIEWING CONE.  THE INTEGRAL IS         *
C  *           PERFORMED BY FINDING THE PATH LENGTHS OF THE LOS IN     *
C  *           EACH PLASMA CELL.                                       *
C  *                                                                   *
C  *            CHRIS FARRELL  (HUNTERSKIL)  MARCH 1989                *
C  *            LORNE HORTON   (JET)         JULY  1993                *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE "PARAMS"
      include 'params'
C     INCLUDE "CGEOM"
      include 'cgeom'
      INTEGER NUMTHE,AVPTS
      REAL    TVALS(MAXTHE),TOUTS(MAXTHE),TWIDS(MAXTHE),
     >        ROBS,ZOBS,VS(MAXNKS,MAXNRS)
C
      INTEGER I,J,K,IK,IR,NINT,SIDE(2)
      REAL    THETA,XB(2),WB(2),DIST(2)
c slmod begin
      INCLUDE 'slout'

      INTEGER    MAXIND     ,MAXVIEW
      PARAMETER (MAXIND=1000,MAXVIEW=1000)

      INTEGER IN,IN1,IN2,nview,idum1,refopt
      REAL    WDIST,DUM1,SNGL1(MAXIND),WGHT1(0:MAXIND),
     .                   SNGL2(MAXIND),WGHT2(0:MAXIND),
     .        xb2(2*MAXVIEW),wb2(2*MAXVIEW),theta2(MAXVIEW),
     .        weight2(MAXVIEW),rdum1,rdum2,refweight,wght3,delta,frac,
     .        lastval,lastval2,dummy1
      CHARACTER*128 graph6,cdum1


      nshow = 0

      refopt = 0
      READ(5,'(A128)',END=10) graph6
      IF (graph6(8:11).EQ.'Refl'.OR.graph6(8:11).EQ.'REFL'.OR.
     .    graph6(8:11).EQ.'refl') THEN
        READ(graph6,*) cdum1,refopt,refweight,rdum1,rdum2
        WRITE(6,*) 'REFOPT,REFWEIGHT=',refopt,refweight
        WRITE(char(29),'(A,I6  )') 'REFLECTION OPTION=',refopt
        WRITE(char(30),'(A,F6.2)') 'REFLECTION WEIGHT=',refweight
      ELSE
        BACKSPACE 5
      ENDIF
10    CONTINUE
c slmod end
C
      XB(1) = ROBS
      XB(2) = ZOBS
C
C  LOOP OVER SIGHT LINES
C

c...WHAT AN IDIOT -- HOW COULD I POSSIBLY LEAVE THIS IN HERE? -SL, FEB 21, 2004
c      DO ir = irtrap+1, nrs
c        DO ik = 1, nks(ir)
c          IF (rs(ik,ir).GT.0.58) CYCLE
c          wght0(ik,ir) = 0.0
c        ENDDO          
c      ENDDO

      DO 200 I = 1, NUMTHE
        TVALS(I) = 0.0
c slmod begin
        IN2      = 0
        WGHT2(0) = 0.0
        WGHT1(I) = 0.0
c slmod end
C
C  LOOP OVER CHORDS FOR AVERAGING
C
        DO 100 J = 1, AVPTS

          IF (AVPTS.EQ.1) THEN
            THETA = TOUTS(I)
          ELSE
            THETA = TOUTS(I) - 0.5*TWIDS(I) + (J-1)*TWIDS(I)/(AVPTS-1)
          ENDIF
          THETA = THETA*DEGRAD
          WB(1) = COS(THETA)
          WB(2) = SIN(THETA)

c slmod begin
c...      For each view, specify an array of reflected views according
c         to some selected distribution.  Load these reflected views
c         into an array along with assigned weights.  The parent view
c         is the first element in the array and has a weight of 1:

c...      Assign element from parent view:

          nview = 1

          xb2(1)     = xb(1)
          xb2(2)     = xb(2)
          wb2(1)     = wb(1)
          wb2(2)     = wb(2)
          theta2(1)  = theta
          weight2(1) = 1.0

c slmod end
C
C  LOOP OVER PLASMA CELLS
C
c slmod begin
          IF (LOSOPT.GT.0) THEN

            IF (refopt.NE.0) THEN

c...          Generate child views:
              CALL GetReflections(nview,xb2,wb2,theta2,weight2,MAXVIEW,
     .                            refopt,refweight,rdum1,rdum2)   



              DO in1 = 1, nview


                IF (in1.GT.1) THEN
                  frac = 1.0 / REAL(nview - 1)
                ELSE
                  frac = 1.0            
                ENDIF


                dummy1 = 0.0

                in2 = 2 * (in1 - 1) + 1
                CALL slCheckLOS(xb2(in2),wb2(in2),theta2(in1),wdist,
     .                          idum1,rdum1,rdum2)

                DO ir = 1, nrs
                  

                  IF (ir.EQ.1.OR.ir.EQ.irwall.OR.ir.EQ.irtrap) CYCLE
                  DO ik = 1, nks(ir)

                    k = korpg(ik,ir)

                    CALL Inters(nvertp(k),rvertp(1,k),zvertp(1,k),
     .                          xb2(in2),wb2(in2),nint,dist,side)

                    IF (nint.EQ.2.AND.dist(1).LT.wdist.AND.
     .                                dist(2).GT.0.0  ) THEN
                      IF (DIST(1).LT.0.0) THEN
                        DELTA = DIST(2)
                      ELSE
                        DELTA = DIST(2) - DIST(1)
                      ENDIF

                      WGHT3= WGHT0(IK,IR) * DELTA * weight2(in1) * frac

                      dummy1 = dummy1 + vs(ik,ir) * delta * weight2(in1)

                      TVALS(I) = TVALS(I) + WGHT3 * VS(IK,IR)
                      WGHT1(I) = WGHT1(I) + WGHT3

                    ENDIF
                  ENDDO
                ENDDO
              
c             IF (in1.EQ.1) WRITE(0,*) 'WGHT 1:',I,THETA/DEGRAD,WGHT1(I)
c             IF (in1.EQ.2) WRITE(0,*) 'WGHT 2:',I,THETA/DEGRAD,WGHT1(I)


                IF (i.EQ.5) THEN
  
                  IF (in1.EQ.1) THEN
                    lastval = 0.0
                    lastval2 = 0.0
                  ELSE
                  ENDIF

c                  WRITE(0,*)  '  - ',in1,wght1(i) - lastval,
c     .                   tvals(i)-lastval2,dummy1

                    lastval = wght1(i)
                    lastval2 = tvals(i)
                ENDIF

              ENDDO

            ELSE

c              WRITE(0,*) 'VIEW:',i,touts(i)
c              WRITE(6,*) 'VIEW:',i,touts(i)

              CALL slCheckLOS(xb,wb,theta,wdist,idum1,rdum1,rdum2)

              IN1      = 0
              DO IR = 1, NRS
                IF (ir.EQ.1.OR.ir.EQ.irwall.OR.ir.EQ.irtrap) CYCLE
                DO IK = 1, NKS(IR)
                  K = KORPG(IK,IR)
                  CALL INTERS(NVERTP(K),RVERTP(1,K),ZVERTP(1,K),
     .                        XB,WB,NINT,DIST,SIDE)
                  IF (NINT.EQ.2.AND.DIST(1).LT.WDIST.AND.
     .                              DIST(2).GT.0.0  ) THEN
                    IF (DIST(1).LT.0.0) THEN
                      DELTA = DIST(2)
                    ELSE
                      DELTA = DIST(2) - DIST(1)
                    ENDIF
 
                    WGHT3= WGHT0(IK,IR) * DELTA

                    TVALS(I) = TVALS(I) + WGHT3 * VS(IK,IR)
                    WGHT1(I) = WGHT1(I) + WGHT3

                    IF (.FALSE..AND.I.Eq.20) THEN
                      WRITE(6,'(A,3I6,F10.3,1P,5E12.4,0P)') 
     .                  '-',i,ik,ir,vs(ik,ir),wght0(ik,ir),
     .                  (DIST(2) - MAX(0.0,DIST(1))),tvals(i),wght3,
     .                  wght1(i)
c     .                  (DIST(2) - MAX(0.0,DIST(1))),wght3,wght1(i)
                    ENDIF

                  ENDIF
                ENDDO
              ENDDO

            ENDIF

          ELSE

            IF (refopt.NE.0) THEN 
c...          Generate child views:
              CALL GetReflections(nview,xb2,wb2,theta2,weight2,MAXVIEW,
     .                            refopt,refweight,rdum1,rdum2)   
              DO in1 = 1, nview

                IF (in1.GT.1) THEN
                  frac = 1.0 / REAL(nview - 1)
                ELSE
                  frac = 1.0            
                ENDIF

                in2 = 2 * (in1 - 1) + 1
                CALL slCheckLOS(xb2(in2),wb2(in2),theta2(in1),wdist,
     .                          idum1,rdum1,rdum2)
                DO ir = 1, nrs
                  IF (ir.EQ.1.OR.ir.EQ.irwall.OR.ir.EQ.irtrap) CYCLE
                  DO ik = 1, nks(ir)
                    k = korpg(ik,ir)
                    CALL Inters(nvertp(k),rvertp(1,k),zvertp(1,k),
     .                          xb2(in2),wb2(in2),nint,dist,side)
                    IF (nint.EQ.2.AND.dist(1).LT.wdist.AND.
     .                                dist(2).GT.0.0  ) THEN
                      IF (dist(1).LT.0.0) THEN
                        tvals(i) = tvals(i) + vs(ik,ir) * weight2(in1) *
     .                                        dist(2) * frac
                      ELSE
                        tvals(i) = tvals(i) + vs(ik,ir) * weight2(in1) * 
     .                                       (dist(2) - dist(1)) * frac
                      ENDIF
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO


            ELSE


              CALL slCheckLOS(xb,wb,theta,wdist,idum1,rdum1,rdum2)
              DO IR = 1, NRS
c                IF (ir.GE.14.AND.ir.LE.18) CYCLE
                DO IK = 1, NKS(IR)
                  K = KORPG(IK,IR)
                  IF (K.EQ.0) CYCLE
                  CALL INTERS(NVERTP(K),RVERTP(1,K),ZVERTP(1,K),
     .                        XB,WB,NINT,DIST,SIDE)
                   IF (NINT.EQ.2.AND.DIST(1).LT.WDIST.AND.
     .                               DIST(2).GT.0.0  ) THEN
                     IF (DIST(1).LT.0.0) THEN
                       TVALS(I) = TVALS(I) + VS(IK,IR)*DIST(2)
                     ELSE
                       TVALS(I) = TVALS(I) + VS(IK,IR)*(DIST(2)-DIST(1))
                    ENDIF
                  ENDIF
                ENDDO
              ENDDO


c              WRITE(0,*) 'VIEW:',i,touts(i),tvals(i),tvals(MAX(1,i-1))


            ENDIF
c
c            DO IR = 1, NRS
c              DO IK = 1, NKS(IR)
c                K = KORPG(IK,IR)
c                CALL INTERS(NVERTP(K),RVERTP(1,K),ZVERTP(1,K),
c     .                      XB,WB,NINT,DIST,SIDE)
c                IF (NINT.EQ.2.AND.DIST(1).LT.WDIST.AND.
c     .                            DIST(2).GT.0.0  ) THEN
c                  IF (DIST(1).LT.0.0) THEN
c                    TVALS(I) = TVALS(I) + VS(IK,IR)*DIST(2)
c                  ELSE
c                    TVALS(I) = TVALS(I) + VS(IK,IR)*(DIST(2)-DIST(1))
c                  ENDIF
c                ENDIF
c              ENDDO
c            ENDDO
c slmod end
          ENDIF

  100   CONTINUE

        IF (LOSOPT.GT.0) THEN
          TVALS(I) = TVALS(I) / WGHT1(I)
        ELSE
          TVALS(I) = TVALS(I) / AVPTS
        ENDIF
c
c          DO 20 IR = 1, NRS
c            DO 10 IK = 1, NKS(IR)
c              K = KORPG(IK,IR)
c              CALL INTERS(NVERTP(K),RVERTP(1,K),ZVERTP(1,K),
c     >                    XB,WB,NINT,DIST,SIDE)
c              IF (NINT.EQ.2 .AND. DIST(2).GT.0.0) THEN
c                IF (DIST(1).LT.0.0) THEN
c                  TVALS(I) = TVALS(I) + VS(IK,IR)*DIST(2)
c                ELSE
c                  TVALS(I) = TVALS(I) + VS(IK,IR)*(DIST(2)-DIST(1))
c                ENDIF
c
c                write(6,*) 'i:',i,ik,ir,tvals(i),vs(ik,ir),
c     >                     dist(2),dist(1)
c              ENDIF
c   10       CONTINUE
c   20     CONTINUE
c  100   CONTINUE
c        TVALS(I) = TVALS(I) / AVPTS
c slmod end



  200 CONTINUE
C
      RETURN
      END


c
c EMPTY
c
      SUBROUTINE DUM
      END





      SUBROUTINE GetReflections(nview,xb2,wb2,theta2,weight2,MAXVIEW,
     .                          refopt,refweight,rdum1,rdum2)   
      IMPLICIT none


      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slout'

      INTEGER nview,in1,MAXVIEW,i1,refopt,iview
      REAL    WDIST,x0,x1,x2,x3,x4,y0,y1,y2,y3,y4,m,
     .        xb2(2*MAXVIEW),wb2(2*MAXVIEW),theta2(MAXVIEW),
     .        weight2(MAXVIEW),refweight,dangle4,wgtsum,
     .        a,b,c,angle1,angle2,angle3,angle4,rdum1,rdum2,dangle

      LOGICAL message
      DATA message /.TRUE./
      SAVE


c      x0 = xb2(1)
c      y0 = xb2(2)

      angle1 = theta2(1)


c...  Find surface normal:

      CALL slCheckLOS(xb2(1),wb2(1),theta2(1),wdist,in1,x3,y3)

      x1 = wallpt(in1,20)
      y1 = wallpt(in1,21)
      x2 = wallpt(in1,22)
      y2 = wallpt(in1,23)

c      WRITE(0,*) '0:',theta2(1)/degrad
c      WRITE(0,*) '1:',x1,y1,x2,y2
c      WRITE(0,*) '2:',x3,y3,in1

c...  Find angle of surface segment:
      a = SQRT((x2 -         x1)**2.0 + (y2 - y1)**2.0)      
      b = SQRT((x2 - 100.0 * x1)**2.0 + (y2 - y1)**2.0)     
      c = SQRT((x1 - 100.0 * x1)**2.0 + (y1 - y1)**2.0)      

      angle2 = ACOS((b**2.0 - a**2.0 - c**2.0) / (-2.0 * a * c))

      IF (y2.LT.y1) angle2 = 2.0 * PI - angle2

c...  Angle of normal is angle of surface segment less 90 degrees:
      angle2 = angle2 - (PI / 2.0)

c.... Rotate parent angle as necessary:
      IF (ABS(angle1-angle2).GT.(PI/2.0)) angle1 = angle1 - PI       
      IF (angle1.LT.0.0   ) angle1 = angle1 + 2.0 * PI
      IF (angle1.GE.2.0*PI) angle1 = angle1 - 2.0 * PI

c      WRITE(0,*) 'B:',angle1/degrad,angle2/degrad,angle3/degrad

      IF     (ABS(refopt).EQ.1) THEN
c....   Specular reflection:

c....   Calcualte outgoing angle:
        angle3 = (angle2 - angle1) + angle2
        IF (angle3.LT.0.0   ) angle3 = angle3 + 2.0 * PI
        IF (angle3.GE.2.0*PI) angle3 = angle3 - 2.0 * PI

        nview = nview + 1
        xb2(2*nview-1) = x3
        xb2(2*nview  ) = y3
        wb2(2*nview-1) = COS(angle3)
        wb2(2*nview  ) = SIN(angle3)
        theta2(nview)  = angle3 
        weight2(nview) = refweight


         

      ELSEIF (ABS(refopt).EQ.2) THEN
c....   Diffuse reflection:

        DO dangle = -89.0, 89.0, 1.0
c        DO dangle = 1.0, 89.0, 20.0

c....     Calcualte outgoing angle:
          angle3 = angle2 + dangle * PI / 180.0
c          angle3 = (angle2 - (angle2 + dangle)) + angle2
          IF (angle3.LT.0.0   ) angle3 = angle3 + 2.0 * PI
          IF (angle3.GE.2.0*PI) angle3 = angle3 - 2.0 * PI

          nview = nview + 1
          xb2(2*nview-1) = x3
          xb2(2*nview  ) = y3
          wb2(2*nview-1) = COS(angle3)
          wb2(2*nview  ) = SIN(angle3)
          theta2(nview)  = angle3 
          weight2(nview) = refweight

        ENDDO

      ELSEIF (ABS(refopt).EQ.3) THEN
c....   Specular/diffuse reflection:

c....   Calcualte outgoing angle:
        angle3 = (angle2 - angle1) + angle2
        IF (angle3.LT.0.0   ) angle3 = angle3 + 2.0 * PI
        IF (angle3.GE.2.0*PI) angle3 = angle3 - 2.0 * PI

        IF (message) THEN
          WRITE(0,*) 
          WRITE(0,*) '********************************'
          WRITE(0,*) ' HIGH-RES DIFFUSE REFLECTION ON'
          WRITE(0,*) '********************************'
          WRITE(0,*) 
          message = .FALSE.
        ENDIF

        DO dangle = -85.0, 85.0, 1.0
c        DO dangle = -85.0, 85.0, 1.0

c....     Calcualte outgoing angle:
          angle4 = angle2 + dangle * PI / 180.0
          IF (angle4.LT.0.0   ) angle4 = angle4 + 2.0 * PI
          IF (angle4.GE.2.0*PI) angle4 = angle4 - 2.0 * PI

          nview = nview + 1
          xb2(2*nview-1) = x3
          xb2(2*nview  ) = y3
          wb2(2*nview-1) = COS(angle4)
          wb2(2*nview  ) = SIN(angle4)
          theta2(nview)  = angle4 
          IF (ABS(angle3-angle4).LT.(PI/180.0)*(5.0/2.0)) THEN
            weight2(nview) = 0.10
            WRITE(0,'(A,3F10.5,I6)') 'SPEC:',angle2*180/PI,
     .                 angle3*180/PI,angle4*180/PI,nview
          ELSE
            weight2(nview) = refweight
          ENDIF

        ENDDO

      ELSEIF (ABS(refopt).EQ.4) THEN
c....   Proper diffuse reflection:

       DO dangle = -89.0, 89.0, 5.0
c        DO dangle = 1.0, 89.0, 5.0

c....     Calcualte outgoing angle:
          angle3 = angle2 + dangle * PI / 180.0
c          angle3 = (angle2 - (angle2 + dangle)) + angle2
          IF (angle3.LT.0.0   ) angle3 = angle3 + 2.0 * PI
          IF (angle3.GE.2.0*PI) angle3 = angle3 - 2.0 * PI

          nview = nview + 1
          xb2(2*nview-1) = x3
          xb2(2*nview  ) = y3
          wb2(2*nview-1) = COS(angle3)
          wb2(2*nview  ) = SIN(angle3)
          theta2(nview)  = angle3 
          weight2(nview) = refweight * ABS(SIN(dangle*PI/180.0))
        ENDDO


      ELSEIF (ABS(refopt).EQ.5) THEN
c....   50% specular, 50% diffuse:

c...    Specular reflection angle:
        angle4 = (angle2 - angle1) + angle2
        IF (angle4.LT.0.0   ) angle4 = angle4 + 2.0 * PI
        IF (angle4.GE.2.0*PI) angle4 = angle4 - 2.0 * PI
        dangle4 = HI

        DO dangle = -89.0, 89.0, 1.0

c....     Calcualte outgoing angle:
          angle3 = angle2 + dangle * PI / 180.0
          IF (angle3.LT.0.0   ) angle3 = angle3 + 2.0 * PI
          IF (angle3.GE.2.0*PI) angle3 = angle3 - 2.0 * PI

          nview = nview + 1
          xb2(2*nview-1) = x3
          xb2(2*nview  ) = y3
          wb2(2*nview-1) = COS(angle3)
          wb2(2*nview  ) = SIN(angle3)
          theta2(nview)  = angle3 
          weight2(nview) = refweight * ABS(SIN(dangle*PI/180.0))

          IF (ABS(angle3-angle4).LT.dangle4) THEN
            dangle4 = ABS(angle3-angle4)
            iview = nview
          ENDIF

        ENDDO

c...    Sum the weights of the diffuse reflections and assign it to the
c       specular reflection chord:        
        wgtsum = 0.0
        DO i1 = 1, nview
          IF (i1.EQ.iview) CYCLE
          wgtsum = wgtsum + weight2(i1)
        ENDDO
        weight2(iview) = 0.3 * wgtsum


c        WRITE(0,*) 'ANDLGE:',theta2(iview),angle4,angle2,angle1
c        STOP 'dsfsd'

      ELSE
        CALL ER('GetReflections','Sorry, option not supported',*99)
      ENDIF


c...  Output:
      IF (refopt.LT.0.AND.nshow+nview.LT.MAXSHOW) THEN

        DO i1 = 1, nview

          IF (refopt.EQ.-5.AND.i1.GT.1.AND.i1.NE.iview) CYCLE

          nshow = nshow + 1
          rshow(nshow) = xb2(2*i1-1)
          zshow(nshow) = xb2(2*i1  )
          ashow(nshow) = theta2(i1)
        ENDDO
      ENDIF



      RETURN
99    STOP
      END


c          nview = nview + 1
c          xb2(3)     = xb2(1)     
c          xb2(4)     = xb2(2)    
c          wb2(3)     = wb2(1)    
c          wb2(4)     = wb2(2)    
c          theta2(2)  = theta2(1) 
c          weight2(2) = 1.0








c
c
c ROUTINES FROM SLOUT.O5A
c
c




c
c ======================================================================
c
c subroutine: SetupSourcePlot
c
      SUBROUTINE SetupSourcePlot(iref,graph,mode)
      IMPLICIT   none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'
      INCLUDE 'slout'

      COMMON /OLDPLASMA/ oldknbs ,oldktebs ,oldktibs ,oldkvhs ,
     .                   oldknbs2,oldktebs2,oldktibs2,oldkvhs2
      REAL
     .     oldknbs  (MAXNKS,MAXNRS),oldktebs (MAXNKS,MAXNRS),
     .     oldktibs (MAXNKS,MAXNRS),oldkvhs  (MAXNKS,MAXNRS),
     .     oldktebs2(MAXNKS,MAXNRS),oldktibs2(MAXNKS,MAXNRS),
     .     oldknbs2 (MAXNKS,MAXNRS),oldkvhs2 (MAXNKS,MAXNRS)


      INTEGER   iref,mode
      CHARACTER graph*(*)

      INTEGER i1,iter,error1,error2,error3,status,load,ik,ir,
     .        u1,u2,u3

      INTEGER   NUMPLTS
      PARAMETER (NUMPLTS=30)
      INTEGER   plotflag(NUMPLTS),plotload(NUMPLTS)
      DATA      (plotflag(i1),i1=1,NUMPLTS)
     .  /  711, 701, 703, 723, 705, 725, 726, 000, 000, 000,
     .     219, 012, 606, 605, 618, 000, 966, 960, 959, 958,
     .     970, 974, 978, 984, 980, 353, 000, 999, 011, 012/
      DATA      (plotload(i1),i1=1,NUMPLTS)
     .  /    3,   3,   3,   3,   3,   3,   3,   0,   0,   0,
     .       3,   4,   3,   3,   0,   0,   4,   3,   3,   3,
     .       3,   3,   4,   3,   3,   4,   0,   3,   4,   4/

      error1 = 0
      error2 = 0
      error3 = 0
      status = 0

c      WRITE(0,*) '    ',graph(1:15)

      IF     (stepopt.EQ.0) THEN
        READ(graph( 7: 9),'(I3)') iref
        READ(graph(11:12),'(I2)') iter
      ELSEIF (nsteplist.EQ.99) THEN
        READ(graph( 7: 9),'(I3)') iref
        iter = stepopt
      ELSE
        READ(graph( 7: 9),'(I3)') iref
        iter = steplist(stepopt)
      ENDIF

      loadstep = iter

c
c
c
      DO i1 = 1, NUMPLTS
        IF (plotflag(i1).EQ.iref) status = i1
      ENDDO

      IF (status.EQ.0) THEN
        WRITE(0,'(A,I3,A)')
     .    'WARNING SetupSourcePlot: Plot ',iref,' untested'
        load = 3
      ELSE
        load = plotload(status)
      ENDIF

      u1 = 40
      u2 = 41
      u3 = 42

      OPEN(u1,FILE='source.dat',FORM='UNFORMATTED',STATUS='OLD',ERR=99)
      OPEN(u2,FILE='plasma.dat',FORM='UNFORMATTED',STATUS='OLD',ERR=99)
      OPEN(u3,FILE='geomty.dat',FORM='UNFORMATTED',STATUS='OLD',ERR=99)

c      OPEN(PINOUT2,STATUS='UNKNOWN',FORM='UNFORMATTED')
c      OPEN(PINOUT3,STATUS='UNKNOWN',FORM='UNFORMATTED')
c      OPEN(PINOUT4,STATUS='UNKNOWN',FORM='UNFORMATTED')
c
c      WRITE(0,*) 'IREF,ITER: ',iref,iter

      DO i1 = 0, iter
c        WRITE(0,'(A,I3,A)') 'Reading iteration ',i1,'...'


        IF (load.EQ.1.OR.load.EQ.3.OR.load.EQ.4)
     .    CALL ReadSources (u1,error1)
c     .    CALL ReadSources (PINOUT2,error1)
        IF (load.EQ.2.OR.load.EQ.3.OR.load.EQ.4)
     .    CALL ReadPlasma  (u2,error2)
c     .    CALL ReadPlasma  (PINOUT3,error2)
        IF (rel_step.EQ.0.OR.adp_opt.GT.0.OR.load.EQ.4) THEN
c          WRITE(0,*) 'MAK: REL_STEP= ',rel_step,adp_opt,load
c          WRITE(0,*) 'LOADING GEOMETRY DATA  STEP= ',rel_step
          CALL ReadGeometry(u3,error3)
c          CALL ReadGeometry(PINOUT4,error3)
          IF (error3.NE.0) THEN
            CALL WN('SetupSourcePlot','Geometry data not found')
            error3 = 0
          ENDIF
        ENDIF

c        WRITE(0,*) 'MARK: IKBOUNDS= ',ikbound(3,IKLO),ikbound(3,IKHI)

c          WRITE(0,*) 'NBR -:',nbr

        IF (error1.NE.0.OR.error2.NE.0.OR.error3.NE.0) THEN
          IF (nsteplist.NE.99) 
     .      WRITE(0,'(A)') 'ERROR SetupSourcePlot: Source data not '//
     .                     'found'
     .
 

         IF (error3.NE.0)
     .      CALL WN('SSP','Problem with geometry data')
          iref =  0
          mode = -1


          CLOSE(u1)
          CLOSE(u2)
          CLOSE(u3)
c          CLOSE(PINOUT2)
c          CLOSE(PINOUT3)
c          CLOSE(PINOUT4)

          IF (nsteplist.EQ.99) nsteplist = 0



          RETURN

        ENDIF
      ENDDO

      mode = 1

c     WRITE(0,'(A,I3,A,3I3)')
c    .  'Substituting 954 for ',iref,'  STEP plot,STEP,ITER =',
c    .  iter,rel_step,rel_iter
      WRITE(6,'(A,I3,A,I3)')
     .  'Substituting 954 for ',iref,', iteration ',iter

c      CLOSE(PINOUT2)
c      CLOSE(PINOUT3)
c      CLOSE(PINOUT4)
          CLOSE(u1)
          CLOSE(u2)
          CLOSE(u3)

c...  Scaling the flow velocity:
      DO ir = 1, MAXNRS
        DO ik = 1, MAXNKS
          kvhs(ik,ir) = kvhs(ik,ir) * qt
        ENDDO
      ENDDO

      RETURN
99    WRITE(0,'(5X,A)') graph
      WRITE(0,'(5X,A,2I4)') 'IREF ITER     = ',iref,iter
      WRITE(0,'(5X,A,4I4)') 'LOAD ERROR1-3 = ',
     .  load,error1,error2,error3
      iref =  0
      mode = -1

c      CLOSE(PINOUT2)
c      CLOSE(PINOUT3)
c      CLOSE(PINOUT4)
          CLOSE(u1)
          CLOSE(u2)
          CLOSE(u3)

      DO ir = 1, MAXNKS
        DO ik = 1, MAXNKS
          oldktebs(ik,ir) = ktebs(ik,ir)
          oldktibs(ik,ir) = ktibs(ik,ir)
          oldknbs (ik,ir) = knbs (ik,ir)
          oldkvhs (ik,ir) = kvhs (ik,ir)
        ENDDO
      ENDDO

      RETURN
      END
c
c ======================================================================
c
c function: RingNo
c
c This was copied from the DIVIMP source module plasma.d6a so that 
c the ShiftTargetData subroutine would compile with OUT.  It would be 
c better if RingNo were moved to utility.u6a, but I don't want to 
c make such an extensive change (well, its not that big) to the main 
c DIVIMP modules at this time (SL - Mar 29, 2000):
c
      INTEGER FUNCTION RINGNO(IR,LPDAT,NLPDAT,MAXN1,MAXN2,ierr)
      implicit none
      INTEGER IR,NLPDAT
      INTEGER MAXN1,MAXN2,ierr
      REAL LPDAT(MAXN1,MAXN2)
C
      INTEGER I
C
C     THIS ROUTINE CHECKS THROUGH THE DATA ARRAY TO SEE IF THE
C     FIRST VALUE MATCHES THE RING NUMBER BEING SOUGHT. IF IT MATCHES
C     IT RETURNS THE INDEX INTO THE ARRAY, IF IT DOESN'T MATCH IT RETURN
C     THE INDEX OF THE LAST SET OF DATA.
C
      ierr = 0
c
      DO 100 I = 1,NLPDAT
         IF (INT(LPDAT(I,1)).EQ.IR) THEN
            GOTO 200
         ENDIF
100   CONTINUE
      I=NLPDAT
      WRITE(6,*) 'DATA NOT FOUND FOR RING - ',IR,' USED ',
     >           INT(LPDAT(I,1)),' INSTEAD'
      ierr = 1
200   RINGNO = I
      RETURN
      END
c
c
c


c
c
c dummy 
c
c

      SUBROUTINE OutputData(fp,comment)
      IMPLICIT none

      INTEGER   fp
      CHARACTER comment*(*)

      RETURN
      END
