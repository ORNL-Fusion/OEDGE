c     -*-Fortran-*-
c ======================================================================
c ======================================================================
c
c
c ======================================================================
c ======================================================================
c
c block: 
c
c PROCESSTRIANGLES
c READTRIANGLES
c DUMPTRIANGLES
c
c ======================================================================
c
c
c
c ======================================================================
c
c subroutine: ProcessTriangles
c
      SUBROUTINE ProcessTriangles
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'


      REAL GetEAD

      INTEGER    TRI_BX  ,TRI_BY  ,TRI_BZ  
      PARAMETER (TRI_BX=1,TRI_BY=2,TRI_BZ=3)
      
c...  LOAD DENSITIES AND PLOT THE DIFFERENCE, PERHAPS ALSO THE INCREASE IN THE
c     IONISATION RATE...


      INTEGER vern,trin,coln

      INTEGER, ALLOCATABLE :: trivi(:,:),colt(:)
      REAL   , ALLOCATABLE :: verx(:),very(:),verz(:),tridata(:,:),
     .                        tricenx(:),triceny(:),tricenz(:)
      

      INTEGER fp,i1,i2,idum1
      REAL*8  a1,a2,b1,b2,c1,c2,d1,d2,t1,t2






      WRITE(0,*) 'HERE IN PROCESS TRIANGLES'


c...  Load in vertex data:
      fp = 99
      OPEN(UNIT=fp,FILE='objects.npco_char',ACCESS='SEQUENTIAL',
c      OPEN(UNIT=fp,FILE='triageom.npco_char',ACCESS='SEQUENTIAL',
     .     STATUS='OLD',ERR=96)      
      READ(fp,*) vern
      ALLOCATE(verx(vern))
      ALLOCATE(very(vern))
      ALLOCATE(verz(vern))
      DO i1 = 1, vern
        READ(fp,*,ERR=95) idum1,verx(i1),very(i1)
      ENDDO
      CLOSE (fp)

c...  Load in triangle vertex indices:
      fp = 99
      OPEN(UNIT=fp,FILE='objects.elemente',ACCESS='SEQUENTIAL',
c      OPEN(UNIT=fp,FILE='triageom.elemente',ACCESS='SEQUENTIAL',
     .     STATUS='OLD',ERR=96)      
      READ(fp,*) trin
      ALLOCATE(trivi(trin,3))
      DO i1 = 1, trin
        READ(fp,*,ERR=95) idum1,(trivi(i1,i2),i2=1,3)
      ENDDO
      CLOSE (fp)

c...  Load in triangle vertex indices:
      fp = 99
      OPEN(UNIT=fp,FILE='objects.plasma',ACCESS='SEQUENTIAL',
c      OPEN(UNIT=fp,FILE='triageom.plasma',ACCESS='SEQUENTIAL',
     .     STATUS='OLD',ERR=96)      
      READ(fp,*) ! Has to match the number of header lines in eirene06...
      READ(fp,*)  
      READ(fp,*)
      READ(fp,*)
      READ(fp,*)
      READ(fp,*)
      READ(fp,*) trin

c...THIS NEEDS TO BE READ FROM DATA FILE!
      coln = 11
      ALLOCATE(colt(coln))
      colt(TRI_BX) = 7
      colt(TRI_BY) = 8
      colt(TRI_BZ) = 9

      ALLOCATE(tridata(trin,coln))
      DO i1 = 1, trin
        READ(fp,*,ERR=95) idum1,(tridata(i1,i2),i2=1,coln)
      ENDDO
      CLOSE (fp)



c...  Calculate center of each triangle:
      ALLOCATE(tricenx(trin))
      ALLOCATE(triceny(trin))
      ALLOCATE(tricenz(trin))
      DO i1 = 1, trin
        a1 = DBLE(verx(trivi(i1,1)))
        a2 = DBLE(very(trivi(i1,1)))
        b1 = DBLE(0.5*(verx(trivi(i1,2))+verx(trivi(i1,3))))
        b2 = DBLE(0.5*(very(trivi(i1,2))+very(trivi(i1,3))))
        c1 = DBLE(verx(trivi(i1,2)))
        c2 = DBLE(very(trivi(i1,2)))
        d1 = DBLE(0.5*(verx(trivi(i1,1))+verx(trivi(i1,3))))
        d2 = DBLE(0.5*(very(trivi(i1,1))+very(trivi(i1,3))))

        CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,t1,t2)
     
        tricenx(i1) = a1 + t1 * (b1 - a1)
        triceny(i1) = a2 + t1 * (b2 - a2)
      ENDDO


c...  Output an R,Z list:
      fp = 99
      OPEN(UNIT=fp,FILE='celldata.dat',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
      WRITE(fp,'(A)') ' SHOT: 990429019  TIME: 950'
      WRITE(fp,'(A)') '    IK    IR     R (m)     Z (m)'
      DO i1 = 1, trin
        WRITE(fp,'(2I6,2F10.6)') i1,0,tricenx(i1)/100.0,
     .                                triceny(i1)/100.0
      ENDDO
      CLOSE (fp)        


c...  Output triangles.dat for plotting:
      fp = 99
      OPEN(UNIT=fp,FILE='triangles.dat',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
      WRITE(fp,*) trin
      DO i1 = 1, trin
        WRITE(fp,'(I6,3(2F10.6,2X),2X,2F10.6)') 
     .    i1,
     .    (verx(trivi(i1,i2))/100.0,very(trivi(i1,i2))/100.0,i2=1,3),
     .    tricenx(i1)/100.0,triceny(i1)/100.0  
      ENDDO
      CLOSE (fp)        

      CALL DumpGrid('PROCESSING TRAINGLES')


c SHOT: 116253   TIME: 3000
c    IK    IR     R (m)     Z (m)
c     1     2  1.570436 -0.787203
c     2     2  1.562577 -0.785282
c     3     2  1.544741 -0.779894
c     4     2  1.517825 -0.768564
c     5     2  1.482221 -0.747810





      RETURN
 95   STOP 'FILE ACCESS ERROR'
 96   STOP 'PROBLEM OPENING FILE'
 99   STOP
      END

c ======================================================================
c
c subroutine: ReadTriangles
c
      SUBROUTINE ReadTriangles
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      REAL GetEAD

      INTEGER    NCOLS  ,COLVOL  ,COLDEN  ,COLIOT  ,COLIOO  ,COLREC
     .                  ,COLATM  ,COLPHO  ,COLNET  ,COLEMI  ,COLABS
      PARAMETER (NCOLS=9,COLVOL=3,COLDEN=2,COLIOT=4,COLIOO=6,COLREC=5,
     .                   COLATM=7,COLPHO=8,COLNET=9,COLEMI=10,COLABS=11)
      
c...  LOAD DENSITIES AND PLOT THE DIFFERENCE, PERHAPS ALSO THE INCREASE IN THE
c     IONISATION RATE...

      REAL, ALLOCATABLE :: tridata(:,:)
      

      INTEGER fp,i1,i2,ik,ike,ir,idum1,trin,count
      REAL    sumvol,sumden,sumrec,sumiont,sumiono,sumpho,
     .        rdum1,rdum2,rdum3



      ALLOCATE(tridata(MAXNKS*MAXNRS,NCOLS))

      WRITE(0,*) 'HERE IN READTRIANGLES'

      fp = 99
      OPEN(UNIT=fp,FILE='facelift.dat',ACCESS='SEQUENTIAL',
     .     STATUS='OLD',ERR=96)      

      READ(fp,*) trin
      IF (.FALSE..AND.trin.NE.eirtrin) 
     .  CALL ER('ReadTriangles','Grid incompatability',*97)

      READ(fp,*) 
      DO i1 = 1, trin
        READ(fp,*,ERR=95) idum1,(tridata(i1,i2),i2=1,6)
c        WRITE(0,'(I6,2X,1P,6E12.4)') idum1,(tridata(i1,i2),i2=1,6)
      ENDDO

c...  Process triangle data for each cell:
c      DO ir = irtrap+1, nrs
      DO ir = 2, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        ike = nks(ir)
        IF (ir.LT.irsep) ike = ike - 1
        DO ik = 1, ike
          count = 0
          sumvol = 0.0
          sumden = 0.0
          sumrec = 0.0
          sumiont = 0.0
          sumiono = 0.0
          DO i1 = 1, trin
            IF (ik.EQ.eirtriik(i1).AND.ir.EQ.eirtriir(i1)) THEN
              count = count + 1
              sumvol = sumvol + tridata(i1,COLVOL)  
              sumden = sumden + tridata(i1,COLDEN)*tridata(i1,COLVOL)  
              sumrec = sumrec + tridata(i1,COLREC)*tridata(i1,COLVOL)
              sumiont = sumiont + tridata(i1,COLIOT)*tridata(i1,COLVOL)
              sumiono = sumiono + tridata(i1,COLIOO)*tridata(i1,COLVOL)
            ENDIF
          ENDDO

          sumden = sumden / sumvol
          sumrec = sumrec / sumvol
          sumiont = sumiont / sumvol
          sumiono = sumiono / sumvol

c...      Change from EIRENE to DIVIMP units:
          sumvol = sumvol * 1.0E-06 * eirtorfrac
          sumden = sumden * 1.0E+06
          sumrec = sumrec 
          sumiont = sumiont
          sumiono = sumiono

          rdum1 = GetEAD(ktebs(ik,ir),knbs(ik,ir),1,'H.4 ')
          rdum2 = GetEAD(ktebs(ik,ir),knbs(ik,ir),3,'H.4 ')

c          WRITE(0,'(A,2I4,I5,F10.2,1P,4(2E10.2,1X),2E14.4,0P,2F10.4)') 
c     .      'TRI:',ik,ir,count,ktebs(ik,ir),
c     .      sumvol,kvols(ik,ir),
c     .      sumden,knbs(ik,ir),
c     .      sumiont,rdum1,
c     .      sumrec,rdum2,
c     .      sumiono,sumiono/rdum1
          

          mulion(ik,ir) = sumiono / rdum1



          IF (ABS(sumvol-kvols(ik,ir)).GT.0.01*sumvol) 
     .      STOP 'VOLUME ERROR'
c          IF (ABS(sumden-knbs(ik,ir) ).GT.0.01*sumden)
c     .      STOP 'DENSITY ERROR'
          
        ENDDO
      ENDDO

      CLOSE (fp)


      fp = 99
      OPEN(UNIT=fp,FILE='facelift.dat.2',ACCESS='SEQUENTIAL',
     .     STATUS='OLD',ERR=96)      
      READ(fp,*) trin
      IF (.FALSE..AND.trin.NE.eirtrin) 
     .  CALL ER('ReadTriangles','Grid incompatability',*97)
      READ(fp,*) 
      tridata = 0.0
      DO i1 = 1, trin
        READ(fp,*,ERR=95) idum1,(tridata(i1,i2),i2=1,NCOLS)
      ENDDO
c...  Process triangle data for each cell:
      DO ir = 2, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        ike = nks(ir)
        IF (ir.LT.irsep) ike = ike - 1
        DO ik = 1, ike
          count = 0
          sumvol = 0.0
          sumpho = 0.0
          DO i1 = 1, trin
            IF (ik.EQ.eirtriik(i1).AND.ir.EQ.eirtriir(i1)) THEN
              count = count + 1
              sumvol = sumvol + tridata(i1,COLVOL)  
              sumpho  = sumpho + tridata(i1,COLDEN)*tridata(i1,COLVOL)
c              WRITE(0,*) 'SUMPHO:',sumpho,sumvol
            ENDIF
          ENDDO

          sumpho = sumpho / sumvol

c...      Change from EIRENE to DIVIMP units:
          sumvol = sumvol * 1.0E-06 * eirtorfrac
          sumpho = sumpho * 1.0E+06

          eirpho1(ik,ir) = sumpho 

c          WRITE(0,'(A,2I4,I5,F10.2,1P,2(2E10.2,1X)')
c     .      'TRI:',ik,ir,count,ktebs(ik,ir),
c     .      sumvol,kvols(ik,ir),
c     .      eirpho1(ik,ir),0.0

          IF (ABS(sumvol-kvols(ik,ir)).GT.0.01*sumvol) 
     .      STOP 'VOLUME ERROR'
          
        ENDDO
      ENDDO
      CLOSE (fp)


      fp = 99
      OPEN(UNIT=fp,FILE='facelift.dat.3',ACCESS='SEQUENTIAL',
     .     STATUS='OLD',ERR=96)      
      READ(fp,*) trin
      IF (.FALSE..AND.trin.NE.eirtrin) 
     .  CALL ER('ReadTriangles','Grid incompatability',*97)
      READ(fp,*) 
      tridata = 0.0
      DO i1 = 1, trin
        READ(fp,*,ERR=95) idum1,(tridata(i1,i2),i2=1,NCOLS)
      ENDDO
c...  Process triangle data for each cell:
      DO ir = 2, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        ike = nks(ir)
        IF (ir.LT.irsep) ike = ike - 1
        DO ik = 1, ike
          count = 0
          sumvol = 0.0
          sumpho = 0.0
          DO i1 = 1, trin
            IF (ik.EQ.eirtriik(i1).AND.ir.EQ.eirtriir(i1)) THEN
              count = count + 1
              sumvol = sumvol + tridata(i1,COLVOL)  
              sumpho  = sumpho + tridata(i1,COLDEN)*tridata(i1,COLVOL)
            ENDIF
          ENDDO

          sumpho = sumpho / sumvol

c...      Change from EIRENE to DIVIMP units:
          sumvol = sumvol * 1.0E-06 * eirtorfrac
          sumpho = sumpho * 1.0E+06

          eirpho2(ik,ir) = sumpho 

          IF (ABS(sumvol-kvols(ik,ir)).GT.0.01*sumvol) 
     .      STOP 'VOLUME ERROR'
          
        ENDDO
      ENDDO
      CLOSE (fp)


      DEALLOCATE(tridata)

      RETURN
 95   WRITE(0     ,*) 'ERROR READING TRIANGLE DATA FILE'
      STOP
 96   WRITE(0     ,*) 'TRIANGLE INPUT FILE NOT FOUND'
      WRITE(PINOUT,*) 'TRIANGLE INPUT FILE NOT FOUND'
      RETURN
 97   WRITE(0,*) 'TRIN=',trin,eirtrin
 99   STOP
      END
c
c ======================================================================
c ======================================================================
c
c block: Structuring grid before passing to EIRENE
c
c STRUCTUREGRID
c
c
c ======================================================================
c
c subroutine: StructureGrid
c
      SUBROUTINE StructureGrid(mode)
      IMPLICIT none

c     Input:
      INTEGER mode

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      REAL       TOL
      PARAMETER (TOL=1.0E-07)

      INTEGER i1,ik1,ik2,ir,ir1,maxcells,numadd,id1,id2,ik

      WRITE(0,*) 'STRUCTURING GRID - NEW ROUTINE'
      write(50,*)
      write(50,'(A)') 'Structuring the grid'
      write(50,*)

      CALL OutputData(87,'Before balancing grid')

c...  Need to make sure that all core cells have the same number of
c     cells:
      DO ir = 3, irsep-1
        IF (nks(ir).NE.nks(2))
     .    CALL ER('StructureGrid','The core rings do not all have '//
     .            'the same number of cells',*99)
      ENDDO


c...  Find the longest ring:
      maxcells = 0
      DO ir = 2, irwall-1
        IF (ir.LT.irsep) THEN
          maxcells = MAX(maxcells,nks(ir)-1+nks(ir+irtrap-1))
        ELSE
          maxcells = MAX(maxcells,nks(ir))
        ENDIF
      ENDDO
      WRITE(SLOUT,*) 'MAXCELLS=',maxcells



      IF (.TRUE.) THEN

c...  Add virtual cells to "pad" the grid, so that each ring has the
c     same number of cells:
      DO ir = 2, irwall-1
        ir1 = 0
        IF (ir.LT.irsep.AND.nks(ir)-1+nks(ir+irtrap-1).LT.maxcells) THEN
          ir1 = ir + irtrap - 1
          numadd = maxcells - (nks(ir) - 1 + nks(ir+irtrap-1))
        ELSEIF (ir.GE.irsep.AND.nks(ir).LT.maxcells) THEN
          ir1 = ir 
          numadd = maxcells - nks(ir)
        ENDIF
        WRITE(SLOUT,*) 'ADDING CELLS:',ir,ir1,numadd

        IF (ir1.GT.0) THEN
c...      Alternate adding cells to the beginning and ends of the ring:
          DO i1 = 1, numadd, 2
            IF (i1  .LE.numadd) THEN
              CALL InsertCell(1       ,ir1,BEFORE,VIRTUAL)
            ENDIF
            IF (i1+1.LE.numadd) THEN
              CALL InsertCell(nks(ir1),ir1,AFTER ,VIRTUAL)
            ENDIF
          ENDDO
        ENDIF

      ENDDO

      ENDIF




c...  Are the virtual cells fixed up?
      CALL OutputData(87,'Grid balanced')

      CALL BuildMap
      CALL SetupGrid



c      STOP 'STRUCTURING GRID'

c...  Check the integrity of the grid:
      DO ir = 2, nrs
        IF (idring(ir).EQ.-1) CYCLE
        DO ik = 1, nks(ir)-1
          id1 = korpg(ik,ir)
          id2 = korpg(ik+1,ir)
          IF (.NOT.(ABS(rvertp(3,id1)-rvertp(2,id2)).LT.TOL.AND.
     .              ABS(zvertp(3,id1)-zvertp(2,id2)).LT.TOL.AND.
     .              ABS(rvertp(4,id1)-rvertp(1,id2)).LT.TOL.AND.
     .              ABS(zvertp(4,id1)-zvertp(1,id2)).LT.TOL)) THEN     
            WRITE(0,*) 'PROBLEM WITH GRID:',ik,ir,'. HALTING PROGRAM.'
  
            CALL OutputData(85,'Error detected when structuring grid')
            STOP
          ENDIF
        ENDDO
      ENDDO



      RETURN
99    WRITE(EROUT,'(5X,A,I4)') ' IR = ',ir
      CALL SetupGrid
      CALL OutputGrid(87,'Severe problems when structuring grid')
      STOP ' '
      END
c
c ========================================================================
c
c
c subroutine: UnstructureGrid
c
c Makes sure that the grid appears to be structured.
c
c Assume the core is fine... add later...
c
      SUBROUTINE UnstructureGrid
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      INTEGER ik,ir
c
      WRITE(0,*) 'UNSTRUCTURING GRID'

      DO ir = irsep, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE

        DO ik = 1, nks(ir)
          IF (virtag(1,ir).EQ.1) THEN
            CALL DeleteCell(1,ir)
          ELSE
            GOTO 10
          ENDIF
        ENDDO
10      CONTINUE

        DO ik = nks(ir), 1, -1
          IF (virtag(nks(ir),ir).EQ.1) THEN
            CALL DeleteCell(nks(ir),ir)
          ELSE
            GOTO 20
          ENDIF
        ENDDO
20      CONTINUE

      ENDDO

      vpolyp = vpolmin

      CALL SetupGrid
      CALL BuildMap


      CALL OutputData(85,"Done unstructuring grid")
      WRITE(0,*) 'DONE'

      RETURN
      END
c
c ======================================================================
c ======================================================================
c
c block: Neutral wall construction
c
c BUILDNEUTRALWALL
c
c
c ======================================================================
c
c subroutine: AssignNIMBUSWall
c
c
c
c
      SUBROUTINE AssignNIMBUSWall
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER ir,iv,i1,id,i2,ndivadsur,in

      REAL       TOL
      PARAMETER (TOL=1.0E-6)

      ndivadsur = 0

      CALL OutputGrid(87,'Before assigning NIMBUS wall')

      IF (stopopt.EQ.121) WRITE(0,*) 'BUILDING NIMBUS WALL'


c...  STILL NEED TO ASSIGN JVESM!

      IF     ((grdnmod.NE.0.AND.grdmod(1,1).NE.887.0).OR.
     .        iflexopt(8).EQ.11) THEN
c...                          980116023:
 

        IF (.TRUE..OR.eirspdatmode.EQ.3) THEN
c...      Match up NIMBUS and DIVIMP wall indicies (xVESM and WALLPT arrays):
          i1 = 1
c...      Both the surface properities and puff specification options 
c         need to be up to date:
          IF (eirpmode.GT.0.AND.eirpmode.LT.4) 
     .      CALL ER('AssignNIMBUSWall','*E16 and *070 input revisions'//
     .              ' not compatible',*99)
        ELSE
c...      Start the NIMBUS wall at IRTRAP on the "outer" target -- this
c         is how the wall is organized when calling NIMBUS:
          i1 = 0
          DO i2 = wallpts, 1, -1
            IF (wallpt(i2,16).EQ.1.0) i1 = i2
          ENDDO
          IF (i1.EQ.0) CALL ER('AssignNIMBUSWall','Bad wall data',*99)
        ENDIF

        nvesm = 0
        DO WHILE (nvesm.LT.wallpts)
          nvesm = nvesm + 1
          rvesm(nvesm,1) = wallpt(i1,20)
          zvesm(nvesm,1) = wallpt(i1,21)
          rvesm(nvesm,2) = wallpt(i1,22)
          zvesm(nvesm,2) = wallpt(i1,23)

          jvesm(nvesm) = NINT(wallpt(i1,16))

          wallpt(i1,17) = REAL(nvesm)

          IF (i1.EQ.wallpts) THEN
            i1 = 1
          ELSE
            i1 = i1 + 1
          ENDIF
        ENDDO    

c...    Assign NIMINDEX:
        DO in = 1, nds
          id = korpg(ikds(in),irds(in))
          nimindex(in) = 0
          DO i1 = 1, nvesm
            IF (ikds(in).EQ.1) THEN
c...          Low IK index target:
              IF (ABS(rvertp(1,id)-rvesm(i1,1)).LT.TOL.AND.
     .            ABS(zvertp(1,id)-zvesm(i1,1)).LT.TOL.AND.
     .            ABS(rvertp(2,id)-rvesm(i1,2)).LT.TOL.AND.
     .            ABS(zvertp(2,id)-zvesm(i1,2)).LT.TOL) THEN
                IF (nimindex(in).EQ.0) THEN
                  nimindex(in) = i1
                ELSE
                  CALL ER('AssignNIMBUSWall','Multiple wall indecies '//
     .                    'identified',*99)
                ENDIF 
              ENDIF
            ELSE
c...          High IK index target:
              IF (ABS(rvertp(3,id)-rvesm(i1,1)).LT.TOL.AND.
     .            ABS(zvertp(3,id)-zvesm(i1,1)).LT.TOL.AND.
     .            ABS(rvertp(4,id)-rvesm(i1,2)).LT.TOL.AND.
     .            ABS(zvertp(4,id)-zvesm(i1,2)).LT.TOL) THEN
                IF (nimindex(in).EQ.0) THEN
                  nimindex(in) = i1
                ELSE
                  CALL ER('AssignNIMBUSWall','Multiple wall indecies '//
     .                    'identified',*99)
                ENDIF 
              ENDIF
            ENDIF
          ENDDO
        ENDDO

c...    Avoid the rest of this routine:
        RETURN

      ELSEIF (thesis) THEN
c...    C-Mod modelling:

        STOP 'CODE TAGGED FOR DELETION IN ASSIGNNIMBUS WALL'
        WRITE(0,*) 'JVESM NOT ASSIGNED!'

        nvesm = wallpts
        DO i1 = wallpts, 1, -1
          rvesm(i1,1) = wallpt2(i1,1)
          zvesm(i1,1) = wallpt2(i1,2)
          IF (i1.EQ.wallpts) THEN
            i2 = 1
          ELSE
            i2 = i1 + 1
          ENDIF
          rvesm(i1,2) = wallpt2(i2,1)
          zvesm(i1,2) = wallpt2(i2,2)
        ENDDO

      ELSEIF (stopopt.EQ.121.AND.wallpts.EQ.0) THEN

        WRITE(0,*) 'JVESM NOT ASSIGNED!'

c...    Broken grid, non-C-Mod modelling:
c
c       Assign neutral wall from the WALLCO array.  This is necessary
c       if the WALLPT array has not been assigned yet -- this method will
c       fail (as is) when assigning WALLPT(x,17) below:
        nvesm = 0
        i1    = 0
        DO WHILE (nvesm.LT.nwall-1)
          nvesm = nvesm + 1
          i1    = i1    + 1
          rvesm(nvesm,1) = wallco(i1  ,1)
          zvesm(nvesm,1) = wallco(i1  ,2)
          rvesm(nvesm,2) = wallco(i1+1,1)
          zvesm(nvesm,2) = wallco(i1+1,2)
        ENDDO
        i1 = 0
        DO WHILE (i1.LT.nwall2-1)
          nvesm = nvesm + 1
          i1    = i1    + 1
          rvesm(nvesm,1) = wallco2(i1  ,1)
          zvesm(nvesm,1) = wallco2(i1  ,2)
          rvesm(nvesm,2) = wallco2(i1+1,1)
          zvesm(nvesm,2) = wallco2(i1+1,2)
        ENDDO

      ELSE

c...    Non-broken grid, non-thesis modelling:

        IF (eirspdatmode.EQ.3) THEN
c...      Match up NIMBUS and DIVIMP wall indicies (xVESM and WALLPT arrays):
          i1 = 1
          
c...      Both the surface properities and puff specification options 
c         need to be up to date:
          IF (eirpmode.GT.0.AND.eirpmode.LT.4) 
     .      CALL ER('AssignNIMBUSWall','*E16 and *070 input revisions'//
     .              ' not compatible',*99)
        ELSE
c...      Start the NIMBUS wall at IRTRAP on the "outer" target -- this
c         is how the wall is organized when calling NIMBUS:
          i1 = 0
          DO i2 = wallpts, 1, -1
            IF (wallpt(i2,16).EQ.1.0) i1 = i2
          ENDDO
          IF (i1.EQ.0) CALL ER('AssignNIMBUSWall','Bad wall data',*99)
        ENDIF

        nvesm = 0
        DO WHILE (nvesm.LT.wallpts)
          nvesm = nvesm + 1
          rvesm(nvesm,1) = wallpt(i1,20)
          zvesm(nvesm,1) = wallpt(i1,21)
          rvesm(nvesm,2) = wallpt(i1,22)
          zvesm(nvesm,2) = wallpt(i1,23)
c
c jdemod - added calculation of JVESM for this grid option - needs to be added
c          for other cases.  
c
c 'OT'  Oueter Target       1
c 'OC'  Outer Corner        2
c 'OD'  Outer Divertor Wall 3   (Vessel wall roughly below X-point I think)
c 'IT'  Inner Target        4
c 'IC'  Inner Corner        5
c 'ID'  Inner Divertor Wall 6
c 'MS'  Main Vessel Wall    7
c 'PV'  Private Plasma Wall 8
c 'CB'  Compound Baffle     9
c 'BA'  Baffle Segment      10
c       Other               0
c
c         Use the simple solution for now - extra information useful but not required - 
c         the critical items are OUTER target, INNER target, MAIN wall and PFZ wall - which
c         are all setup in the walls routine - copy this to jvesm. 
c
          jvesm(nvesm) = wallpt(i1,16)
c
c          IF (wallpt(i1,16).EQ.7) THEN
c            IF ((zvesm(nvesm,2).GE.0.5*(zxp+z0).AND.
c     .           zvesm(nvesm,1).LT.0.5*(zxp+z0)).OR.
c     .          (zvesm(nvesm,2).LE.0.5*(zxp+z0).AND.
c     .           zvesm(nvesm,1).GT.0.5*(zxp+z0))) THEN
c              IF     (cvesm.EQ.3) THEN
c                jvesm(nvesm) = cvesm
c                cvesm = 7
c              ELSEIF (cvesm.EQ.7) THEN
c                cvesm = 6
c                jvesm(nvesm) = cvesm
c              ELSE
c                IF (stopopt.EQ.121) THEN
c
c                ELSE
c                  CALL ER('WriteBlock03b','Unable to assign JVESM',*99)
c                ENDIF
c              ENDIF
c            ELSE
c             jvesm(nvesm) = cvesm
c            ENDIF
c          ELSE
c            jvesm(nvesm) = wallpt(i1,16)
c          ENDIF
c
c
c jdemod         
c
          IF (i1.EQ.wallpts) THEN
            i1 = 1
          ELSE
            i1 = i1 + 1
          ENDIF
c
        ENDDO

c        nvesm = wallpts
c        DO i1 = wallpts, 1, -1
c          rvesm(i1,1) = wallpt(i1,20)
c          zvesm(i1,1) = wallpt(i1,21)
c          rvesm(i1,2) = wallpt(i1,22)
c          zvesm(i1,2) = wallpt(i1,23)
c        ENDDO

      ENDIF

c...  Additional surfaces specified in DIVIMP:
      DO i1 = 1, eirnasdat
c...    Only want to include 2D additional surfaces:
        IF (eirasdat(i1,1).NE.1.0) CYCLE

        ndivadsur = ndivadsur + 1

        IF     (eirasdat(i1,2).EQ.1.0) THEN
c...      Data point specified in DIVIMP input file:
          rvesm(i1+nvesm,1) = eirasdat(i1,3)
          zvesm(i1+nvesm,1) = eirasdat(i1,4)
        ELSEIF (eirasdat(i1,2).EQ.2.0) THEN
c...      Data point references grid:
          ir = INT(ABS(eirasdat(i1,3)))
          iv = INT(    eirasdat(i1,4) )
          IF   (eirasdat(i1,3).LT.0.0) THEN
            id = korpg(1      ,ir)
          ELSE
            id = korpg(nks(ir),ir)
          ENDIF
          rvesm(i1+nvesm,1) = rvertp(iv,id)
          zvesm(i1+nvesm,1) = zvertp(iv,id)
        ELSE
          CALL ER('AssignNIMBUSWall','Unsupported vertex code',*99)
        ENDIF

        IF     (eirasdat(i1,5).EQ.1.0) THEN
          rvesm(i1+nvesm,2) = eirasdat(i1,6)
          zvesm(i1+nvesm,2) = eirasdat(i1,7)
        ELSEIF (eirasdat(i1,5).EQ.2.0) THEN
          ir = INT(ABS(eirasdat(i1,6)))
          iv = INT(    eirasdat(i1,7) )
          IF   (eirasdat(i1,6).LT.0.0) THEN
            id = korpg(1      ,ir)
          ELSE
            id = korpg(nks(ir),ir)
          ENDIF
          rvesm(i1+nvesm,2) = rvertp(iv,id)
          zvesm(i1+nvesm,2) = zvertp(iv,id)
        ENDIF
      ENDDO
      nvesp = ndivadsur
c      nvesm = nvesm + eirnasdat

      IF (stopopt.EQ.121) THEN

c...    Processing of the DIVIMP specified additional surfaces
c       is not required here (at the moment) because the additional
c       surfaces are included in the WALLPT array when 
c       BuildNeutralWall is called.

c...    Argh! This a real mess now, since Dave wants these
c       additional surfaces to be accounted for using NVESP, so
c       that all the neutral wall options -- which only makes sense 
c       of course, just more work for the thesis code...

        WRITE(0,*)
        WRITE(0,*) ' *** SETTING NVESP TO ZERO *** '
        WRITE(0,*)

        nvesp = 0
      ENDIF

      CALL OutputData(85,'Mapping to PIN wall')


      IF (grdnmod.NE.0.AND.grdmod(1,1).NE.887.0) THEN



      ELSEIF (thesis) THEN
c      IF (thesis.OR.stopopt.EQ.121) THEN
c...    Map WALLPT(x,17) to xVESM arrays:
        DO i1 = 1, wallpts
          wallpt(i1,17) = -1.0
          DO i2 = 1, nvesm + nvesp
            IF (ABS(wallpt(i1,20)-rvesm(i2,1)).LT.TOL.AND.
     .          ABS(wallpt(i1,21)-zvesm(i2,1)).LT.TOL.AND.
     .          ABS(wallpt(i1,22)-rvesm(i2,2)).LT.TOL.AND.       
     .          ABS(wallpt(i1,23)-zvesm(i2,2)).LT.TOL)
     .        wallpt(i1,17) = REAL(i2)
          ENDDO
          IF (wallpt(i1,17).EQ.-1.0)
     .      CALL ER('AssignNIMBUSWall','Cannot map DIVIMP wall to '//
     .                                 'NIMBUS wall - A',*99)
        ENDDO
      ELSEIF (stopopt.EQ.121.AND.wallpts.EQ.0) THEN
        WRITE(0,*) 'FIX WALLPT(x,17) ASSIGNMENT'
      ELSE
        DO i1 = 1, wallpts
          wallpt(i1,17) = -1.0
          DO i2 = 1, nvesm + nvesp
            IF (ABS(wallpt(i1,20)-rvesm(i2,1)).LT.TOL.AND.
     .          ABS(wallpt(i1,21)-zvesm(i2,1)).LT.TOL.AND.
     .          ABS(wallpt(i1,22)-rvesm(i2,2)).LT.TOL.AND.       
     .          ABS(wallpt(i1,23)-zvesm(i2,2)).LT.TOL) 
     .         wallpt(i1,17) = REAL(i2)
          ENDDO
          WRITE(SLOUT,'(A,4I6)') 'CRAP: ',i1,i2,wallpts,nvesm
          IF (wallpt(i1,17).EQ.-1.0)
     .      CALL ER('AssignNIMBUSWall','Cannot map DIVIMP wall to '//
     .                                 'PIN wall - B',*99)
        ENDDO
      ENDIF    

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: BuildTargets
c
c Generalized grid replacement for DOTARG.  Target indexing proceeds from
c the high IK index target at IRWALL and proceeds clockwise.
c
c ctargoptg = 6 required.
c
c Quantities assigned:
c
c     IDDS     -
c     IKDS     -
c     IRDS     -
c     RP       - from PLATCO?
c     ZP       - from PLATCO?
c     NDS      -
c     NDSIN    -
c     NDSIN2   -
c     NDSIN3   -
c     THETAS   -
c     DDS      -
c     THETAS2  -
c     DDS2     -
c     COSTET   -
c     SEPDIST  - distance from separatrix (very crude at present)
c     SEPDIST2 - distance from separatrix (not so bad)
c     RSPDIST  - good?
c     ZSPDIST  -
c  
c Need special consideration for FRC grid (sonnet_grid_sub_type=1), and 
c ITER grid (proper double-null).
c
c
      SUBROUTINE BuildTargets
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER FindNum
      REAL    ATAN3C

      INTEGER treg(MAXNRS),ntreg,i1,i2,i3,in,id,ir,
     .        index(0:MAXNRS)
      LOGICAL cont
      REAL    deltar,deltaz,angle(0:MAXNRS)

c      WRITE(0,*) 'HERE IN BUILDTARGETS'

c..ADD CHECK TO MAKE SURE THAT NO RINGS ARE MOVED HERE!
      CALL SequenceGrid

c...  Check sonnet_grid_sub_type.eq.1 - modifications need to be made in SequenceGrid, or 
c     can it just be done here...?

c...  Check for ITER grid - where should that be accommodated?

c...  Check the CTARGIOT=6:

c...  Establish sequence of target regions to be processed:


      DO i1 = 1, grdntreg(IKLO)
        ir = grdtseg(MAX(grdntseg(i1,IKLO)/2,1),i1,IKLO)
        id = korpg(1,ir)
        deltar = 0.5 * (rvertp(1,id) + rvertp(2,id)) - r0
        deltaz = 0.5 * (zvertp(1,id) + zvertp(2,id)) - z0
        angle(i1) = ATAN3C(deltaz,deltar)
        index(i1) = -i1
      ENDDO
      DO i1 = 1, grdntreg(IKHI)
        ir = grdtseg(MAX(grdntseg(i1,IKHI)/2,1),i1,IKHI)
        id = korpg(nks(ir),ir)
        deltar = 0.5 * (rvertp(3,id) + rvertp(4,id)) - r0
        deltaz = 0.5 * (zvertp(3,id) + zvertp(4,id)) - z0
        angle(i1+grdntreg(IKLO)) = ATAN3C(deltaz,deltar)
c        WRITE(0,*) 'HI IR=',i1,ir,ATAN3C(deltaz,deltar)
c        WRITE(0,*) 'HI IR=',deltaz,deltar
        index(i1+grdntreg(IKLO)) = i1
      ENDDO

      DO i1 = 1, grdntreg(IKLO)+grdntreg(IKHI)
c        WRITE(0,*) 'ANGLES:',angle(i1),angle(grdntreg(IKLO)+1)
        IF (angle(i1).GT.angle(grdntreg(IKLO)+1)) 
     .    angle(i1) = angle(i1) - 360.0
c        WRITE(0,*) 'ANGLES:',angle(i1),index(i1)
      ENDDO
   
c...  Sort indecies:
      cont = .TRUE.
      DO WHILE (cont)
        cont = .FALSE.
        DO i1 = 1, grdntreg(IKLO)+grdntreg(IKHI)-1
          IF (angle(i1).LT.angle(i1+1)) THEN
            angle(0) = angle(i1)
            index(0) = index(i1)
            angle(i1) = angle(i1+1)
            index(i1) = index(i1+1)
            angle(i1+1) = angle(0)
            index(i1+1) = index(0)
            cont = .TRUE.
          ENDIF
        ENDDO
      ENDDO

      ntreg = grdntreg(IKLO) + grdntreg(IKHI)
c      WRITE(0,*)
      DO i1 = 1, ntreg
        treg(i1) = index(i1)
c        WRITE(0,*) 'ANGLES:',angle(i1),index(i1)
      ENDDO


      in = 0
      idds = 0
      DO i1 = 1, ntreg
        IF     (treg(i1).LT.0) THEN
c...      Process low IK index target:
          i2 = -treg(i1)
          DO i3 = 1, grdntseg(i2,IKLO)
            in = in + 1
            ir = grdtseg(i3,i2,IKLO)
            idds(ir,2) = in 
            ikds(in) = 1
            irds(in) = ir
            id = FindNum(platco,nplat,ir)
            IF (id.EQ.0) 
     .        CALL ER('BuildTargets','Unable to find low IK '//
     .                'index target coordinates',*99)
            rp(in) = platco(id,4)
            zp(in) = platco(id,5)

c            WRITE(0,*) 'IN  -->',ir,idds(ir,2)
          ENDDO

        ELSEIF (treg(i1).GT.0) THEN
c...      Process high IK index target:
          i2 = treg(i1)
          DO i3 = grdntseg(i2,IKHI), 1, -1
            in = in + 1
            ir = grdtseg(i3,i2,IKHI)
            idds(ir,1) = in 
            ikds(in) = nks(ir)
            irds(in) = ir
            id = FindNum(platco,nplat,ir)
            IF (id.EQ.0) 
     .        CALL ER('BuildTargets','Unable to find high IK '//
     .                'index target coordinates',*99)
            rp(in) = platco(id,2)
            zp(in) = platco(id,3)

c            WRITE(0,*) 'OUT -->',ir,idds(ir,1)
          ENDDO
c...      Set NDSIN:
          IF (i2.EQ.grdntreg(IKHI)) ndsin = in
        ENDIF
      ENDDO
c...  Assign total number of targets:
      nds = in
      ndsin2 = nds
      ndsin3 = ndsin


c      WRITE(0,*) 'DONE IN BUILDTARGETS'

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: BuildNeutralWall
c
c The assumption is made that the neutral wall is well defined, so that no
c processing is required in order for wall segments to be cut to fit the 
c target segments.  These kinds of modifications are made via the ...
c array:
c
      SUBROUTINE BuildNeutralWall
      USE mod_geometry
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      REAL       TOL
      PARAMETER (TOL=5.0E-07)

      REAL   ATAN2C 
      REAL*8 SideLength

      INTEGER in,ik,ir,i1,i2,i3,id,id1,id2,wdat(0:8,6),ndat,
     .        imin,walln,ii,iv,istart,iend,irstart,ir1,ir2,
     .        wallc(2*MAXPTS+1),wallt(MAXPTS),holdc,holdt,ikn,irn,
     .        ik1,ik2,cnt,walln_2
      LOGICAL terminate
      REAL    holdr1,holdz1,holdr2,holdz2,rstart,zstart,
     .        wallr1(MAXPTS+1,2),wallz1(MAXPTS+1,2),
     .        r1,z1,r2,z2,r3,z3,dist,ri,zi,t,rvector,zvector
      REAL*8  d_wallr1(2*MAXPTS+1,2),d_wallz1(2*MAXPTS+1,2)
c
      terminate = .FALSE.

      IF (stopopt.EQ.14) RETURN

      CALL DB('Entering BUILDNEUTRALWALL')
      WRITE(88,*) 'ENTERING BUILDNEUTRALWALL'

      CALL OutputGrid(85,'Before building neutral wall')

      in    = 0
      walln = 0
      wallc = 0
     
      nvesm = 0
      nvesp = 0

c...  Orientation of neutral wall is reveresed, since the DIVIMP
c     specification (CLOCKWISE) is opposite EIRENE's.
c
c     FIX! A check should be made for the neutral wall option, not for 
c     whether or not there is neutral wall data in the input file:
c      IF (sloutput) WRITE(0,*) 'ISSUE WARNING ON WALL OPTION'
c      CALL WN('BuildNeutralWall','Not checking wall option')

c...  Load main wall data:
      IF (iflexopt(8).EQ.11.OR.iflexopt(8).EQ.12.OR.
     .    (cneur.EQ.2.AND.ctrap.EQ.3)) THEN
c...                 980116023:           990429019 new:

        DO i1 = 1, nwall-1
          IF (wallco(i1,1).NE.wallco(i1+1,1).OR.
     .        wallco(i1,2).NE.wallco(i1+1,2)) THEN
      
            walln = walln + 1
            wallr1(walln,1) = wallco(i1,1)
            wallz1(walln,1) = wallco(i1,2)
            wallr1(walln,2) = wallco(i1+1,1)
            wallz1(walln,2) = wallco(i1+1,2)
          ENDIF
        ENDDO

c...    Load private plasma wall:
        DO i1 = 1, nwall2-1
          IF (wallco2(i1,1).NE.wallco2(i1+1,1).OR.
     .        wallco2(i1,2).NE.wallco2(i1+1,2)) THEN
      
            walln = walln + 1
            wallr1(walln,1) = wallco2(i1  ,1)
            wallz1(walln,1) = wallco2(i1  ,2)
            wallr1(walln,2) = wallco2(i1+1,1)
            wallz1(walln,2) = wallco2(i1+1,2)
      
          ENDIF
        ENDDO

      ELSEIF (cneur.EQ.4.AND.ctrap.EQ.4) THEN
c...    The neutral wall is contained in the xVES arrays, which are 
c       typically assigned from the grid file.  It is assumed that the
c       neutral wall specification proceeds clockwise around the vessel:

        DO i1 = 1, nves-1
          IF (rves(i1).EQ.rves(i1+1).AND.
     .        zves(i1).EQ.zves(i1+1)) CYCLE

          walln = walln + 1
          wallr1(walln,1) = rves(i1)
          wallz1(walln,1) = zves(i1)
          wallr1(walln,2) = rves(i1+1)
          wallz1(walln,2) = zves(i1+1)
c
c         jdemod - turn on debugging
c
        ENDDO

      ELSE
        CALL ER('BuildNeutralWall','Unrecognized wall options',*99)
      ENDIF


c
c         jdemod - turn on debugging
c
      write(6,*) 'Debug information is in <casename>.src - pinout'
      do i1 = 1,walln
          WRITE(pinout,'(A,2I8,10(1X,G18.8))') 'WALL_Start:',i1,walln,
     >     wallr1(i1,1),wallz1(i1,1),wallr1(i1,2),wallz1(i1,2)
      end do

          
c      WRITE(0,*) 'FAILURE FOR FLOATING WALL SURFACES'

      write(pinout,*) 'Eirnasdat:',eirnasdat,
     >                 (eirasdat(i1,1),i1=1,eirnasdat)
c
     
      DO i1 = 1, eirnasdat
c       ----------------------------------------------------------------
        IF (eirasdat(i1,1).EQ.1.0.OR.eirasdat(i1,1).EQ.98.0) THEN
          IF (eirasdat(i1,1).EQ.98.0) THEN
            istart = 1
            iend = 1
c...        Count number of data lines in input file:
            DO WHILE (i1+iend+1.LE.eirnasdat) 
              IF (eirasdat(i1+iend+1,1).GE.0.0) EXIT
              iend = iend + 1
            ENDDO
          ELSE
c...        Insert a 2D additional surface (toroidally infinite):
            istart = 0
            iend   = 0            
          ENDIF
          DO i3 = istart, iend
c...        This somewhat awkward loop is due to a (despirate it seems) attempt
c           to reuse some code.  I3=0 for the case where a segment is being inserted,
c           so that the vertex data is just taken from the data line referred
c           to by the I1 loop.  For modifications to the vertex data of existing segments,
c           that were from standard neutral wall input, the I1 data line
c           just tells how many subsequent datalines are present for modifications:
            IF (eirasdat(i1,1).EQ.98.0) THEN
              i2 = NINT(ABS(eirasdat(i1+i3,1)))
            ELSE
              walln = walln + 1
              i2    = walln
            ENDIF
c...        First vertex:
            IF     (eirasdat(i1+i3,2).EQ.1.0) THEN
c...          Data point specified in DIVIMP input file:
              wallr1(i2,1) = eirasdat(i1+i3,3)
              wallz1(i2,1) = eirasdat(i1+i3,4)
              wallc (i2)   = -999
            ELSEIF (eirasdat(i1+i3,2).EQ.2.0) THEN
c...          Data point references grid:
              ir = INT(ABS(eirasdat(i1+i3,3)))
              iv = INT(    eirasdat(i1+i3,4) )
              IF   (eirasdat(i1+i3,3).LT.0.0) THEN
                id = korpg(1      ,ir)
              ELSE
                id = korpg(nks(ir),ir)
              ENDIF
              wallr1(i2,1) = rvertp(iv,id)
              wallz1(i2,1) = zvertp(iv,id)
              wallc (i2)   = -999
            ENDIF
c...        Second vertex:
            IF     (eirasdat(i1+i3,5).EQ.1.0) THEN
              wallr1(i2,2) = eirasdat(i1+i3,6)
              wallz1(i2,2) = eirasdat(i1+i3,7)
              wallc (i2)   = -999
            ELSEIF (eirasdat(i1+i3,5).EQ.2.0) THEN
              ir = INT(ABS(eirasdat(i1+i3,6)))
              iv = INT(    eirasdat(i1+i3,7) )
              IF   (eirasdat(i1+i3,6).LT.0.0) THEN
                id = korpg(1      ,ir)
              ELSE
                id = korpg(nks(ir),ir)
              ENDIF
              wallr1(i2,2) = rvertp(iv,id)
              wallz1(i2,2) = zvertp(iv,id)
              wallc (i2)   = -999
            ENDIF
          ENDDO
c       ----------------------------------------------------------------
        ELSEIF (eirasdat(i1,1).EQ.96.0.OR.
     .          eirasdat(i1,1).EQ.97.0) THEN
c...      Add wall segments that are associated with a specified region of
c         the main SOL virtual ring (IRWALL):
          cnt = 0
          ik1 = 1
          ir1 = irins(1,irwall)

          DO ik = 2, nks(irwall)+1
            IF (ir1.NE.irins(ik,irwall).OR.ik.EQ.nks(irwall)+1) THEN
              cnt = cnt + 1
              ik2 = ik - 1
c
              IF (cnt.GE.NINT(eirasdat(i1,2)).AND.
     .            cnt.LE.NINT(eirasdat(i1,3))) THEN
                id = korpg(ikins(ik1,irwall),ir1)
                walln = walln + 1
                wallr1(walln,1) = rvertp(2,id)
                wallz1(walln,1) = zvertp(2,id)
                IF (eirasdat(i1,1).EQ.96.0) THEN
c...              Construct wall segments by projecting the fluid grid
c                 cells radially outward at the ends of the IRWALL group:
                  t = 0.001
                  wallr1(walln,2) = rvertp(2,id) + t
                  wallz1(walln,2) = zvertp(2,id)
                  id = korpg(ikins(ik2,irwall),ir1)
                  walln = walln + 1
                  wallr1(walln,1) = wallr1(walln-1,2)
                  wallz1(walln,1) = wallz1(walln-1,2)
                  wallr1(walln,2) = rvertp(3,id) + t
                  wallz1(walln,2) = zvertp(3,id)
                  walln = walln + 1
                  wallr1(walln,1) = wallr1(walln-1,2)
                  wallz1(walln,1) = wallz1(walln-1,2)
                  wallr1(walln,2) = rvertp(3,id)
                  wallz1(walln,2) = zvertp(3,id)

c                  DO i2 = walln-2,walln
c                    WRITE(0,*) '  wall:',wallr1(i2,1),
c     .                                   wallz1(i2,1)
c                    WRITE(0,*) '      :',wallr1(i2,2),
c     .                                   wallz1(i2,2)
c                  ENDDO

                ELSE
c...              Build wall segments by projecting radially outward along
c                 the fluid cell poloidal boundaries:
                  t = 2.0
                  DO i2 = ik1, ik2
                    id = korpg(ikins(i2,irwall),irins(i2,irwall))               
                    wallr1(walln,2) = (1.-t)*rvertp(1,id)+t*rvertp(2,id)
                    wallz1(walln,2) = (1.-t)*zvertp(1,id)+t*zvertp(2,id)
                    walln = walln + 1
                    wallr1(walln,1) = wallr1(walln-1,2)
                    wallz1(walln,1) = wallz1(walln-1,2)
                  ENDDO              
                  wallr1(walln,2) = (1.-t)*rvertp(4,id) + t*rvertp(3,id)
                  wallz1(walln,2) = (1.-t)*zvertp(4,id) + t*zvertp(3,id)
                ENDIF
                walln = walln + 1
                wallr1(walln,1) = wallr1(walln-1,2)
                wallz1(walln,1) = wallz1(walln-1,2)
                wallr1(walln,2) = rvertp(3,id)
                wallz1(walln,2) = zvertp(3,id)
              ENDIF
              ik1 = ik
              ir1 = irins(ik1,irwall)
            ENDIF
          ENDDO
c       ----------------------------------------------------------------
        ELSEIF (eirasdat(i1,1).EQ.99.0) THEN
c...      Delete specified segment index from neutral wall specification:
c          WRITE(0,*) 'TRYING TO DELETE:',NINT(eirasdat(i1,2)),
c     .                                   NINT(eirasdat(i1,5))
          DO i2 = NINT(eirasdat(i1,2)), MIN(walln,NINT(eirasdat(i1,5)))
c
c=======
c          DO i2 = NINT(eirasdat(i1,2)), 
c     .            NINT(MAX(eirasdat(i1,3),eirasdat(i1,5)))  ! Added option to use 3rd entry -SL, 30/09/2010
c>>>>>>> .merge-right.r477
c
            wallr1(i2,1) = r0
            wallz1(i2,1) = z0
            wallr1(i2,2) = r0
            wallz1(i2,2) = z0
          ENDDO
c       ----------------------------------------------------------------
        ELSEIF (eirasdat(i1,1).EQ.999.0) THEN
c...      Debug: Save data and avoid attempted wall sequencing:
          nvesm = walln
          DO i2 = 1, walln
            jvesm(i2) = 8
            rvesm(i2,1) = wallr1(i2,1)
            zvesm(i2,1) = wallz1(i2,1)
            rvesm(i2,2) = wallr1(i2,2)
            zvesm(i2,2) = wallz1(i2,2)
          ENDDO 
          terminate = .TRUE.
c...      Leave the loop:
          EXIT
c       ----------------------------------------------------------------
        ELSEIF (eirasdat(i1,1).EQ.998.0) THEN
c...      Setup OSM geometry:
          CALL MapRingstoTubes
          CALL DumpData_OSM('output.trouble1','trouble1')
c...      Automated clipping:

          if (cprint.eq.3.or.cprint.eq.9) then 
            DO i2 = 1, walln
               WRITE(pinout,'(A,I6,2(2F14.7,2X))') 'WALLN, SENT    : ',
     .        i2,wallr1(i2,1),wallz1(i2,1),wallr1(i2,2),wallz1(i2,2)
               WRITE(6,'(A,I6,2(2F14.7,2X))') 'WALLN, SENT    : ',
     .        i2,wallr1(i2,1),wallz1(i2,1),wallr1(i2,2),wallz1(i2,2)
            ENDDO
          endif 
c
c          d_wallr1 = DBLE(wallr1)
c          d_wallz1 = DBLE(wallz1)
c          CALL ClipWallToGrid(walln,d_wallr1,d_wallz1,MAXPTS+1,.TRUE.)
c          wallr1 = SNGL(d_wallr1)
c          wallz1 = SNGL(d_wallz1)           
          d_wallr1(1:walln,:) = DBLE(wallr1(1:walln,:))
          d_wallz1(1:walln,:) = DBLE(wallz1(1:walln,:))
          walln_2 = walln
          CALL osmClipWallToGrid
     .           (walln_2,wallc,d_wallr1,d_wallz1,2*MAXPTS+1)
          walln = walln_2
          wallr1(1:walln,:) = SNGL(d_wallr1(1:walln,:))
          wallz1(1:walln,:) = SNGL(d_wallz1(1:walln,:))           
c
          if (cprint.eq.3.or.cprint.eq.9) then 
             DO i2 = 1, walln
               WRITE(pinout,'(A,I6,2(2F14.7,2X))') 'WALLN, RETURNED: ',
     .           i2,wallr1(i2,1),wallz1(i2,1),wallr1(i2,2),wallz1(i2,2)
               WRITE(6,'(A,I6,2(2F14.7,2X))') 'WALLN, RETURNED: ',
     .           i2,wallr1(i2,1),wallz1(i2,1),wallr1(i2,2),wallz1(i2,2)
             ENDDO
          endif
c...      Wipe the geometry arrays:
          CALL geoClean
          CALL osmClean
        ENDIF
c       ----------------------------------------------------------------
      ENDDO

c...  Process neutral wall data and remove deleted segments:

      write(pinout,'(a,i8)') 'Delete zero length elements:',walln

      i3 = walln
      DO i1 = i3, 1, -1
        IF (wallr1(i1,1).EQ.wallr1(i1,2).AND.
     .      wallz1(i1,1).EQ.wallz1(i1,2)) THEN
c
c     jdemod
c
          WRITE(pinout,'(a,i8,10(1x,g18.8))') 'DELETING:',i1,
     >         wallr1(i1,1),wallr1(i1,2),wallz1(i1,1),wallz1(i1,2)
c
          DO i2 = i1, walln-1
            wallr1(i2,1) = wallr1(i2+1,1)
            wallz1(i2,1) = wallz1(i2+1,1)
            wallr1(i2,2) = wallr1(i2+1,2)
            wallz1(i2,2) = wallz1(i2+1,2)             
          ENDDO
          walln = walln - 1
        ENDIF
      ENDDO

      IF (terminate) THEN
        CALL OutputData(85,'Checking neutral wall')
        CALL DumpGrid('BUILDING WALL SEGMENT LIST')
      ENDIF


c     Add target segments and tag them by assigning WALLC
c     a positive integer value corresponding to the ring number:
      DO ir = irsep, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE

        walln = walln + 1
        wallr1(walln,1) = rvertp(1,korpg(1,ir))
        wallz1(walln,1) = zvertp(1,korpg(1,ir))
        wallr1(walln,2) = rvertp(2,korpg(1,ir))
        wallz1(walln,2) = zvertp(2,korpg(1,ir))
        wallc (walln  ) = ir
        wallt (walln  ) = 1

        walln = walln + 1
        wallr1(walln,1) = rvertp(3,korpg(nks(ir),ir))
        wallz1(walln,1) = zvertp(3,korpg(nks(ir),ir))
        wallr1(walln,2) = rvertp(4,korpg(nks(ir),ir))
        wallz1(walln,2) = zvertp(4,korpg(nks(ir),ir))
        wallc (walln  ) = ir
        wallt (walln  ) = 4
      ENDDO

      do i1 = 1,walln
          WRITE(pinout,'(A,4I8,10(1X,G18.8))') 'WALL_Mid  :',i1,walln,
     >     wallc(i1),wallt(i1),
     >     wallr1(i1,1),wallz1(i1,1),wallr1(i1,2),wallz1(i1,2)
      end do


c...  Sequence neutral wall segments (clockwise), starting
c     at the highest IR index of the low IK index targets that
c     are still in the divertor proper:

      write(pinout,*) 'GRDNMOD:',grdnmod

      IF (grdnmod.NE.0) THEN
        DO ir1 = irsep, irwall-1
          id1 = korpg(1,ir1)
          irstart = 0
          DO ir2 = irsep, irwall-1
            id2 = korpg(1,ir2)
            IF (ABS(rvertp(2,id1)-rvertp(1,id2)).LT.TOL.AND.
     .          ABS(zvertp(2,id1)-zvertp(1,id2)).LT.TOL) THEN
c...          Target segment found to be continuous with a neighbouring
c             segment, so flag ISTART so that we don't start here:
              irstart = ir 
              EXIT
            ENDIF
          ENDDO
          IF (irstart.EQ.0) THEN
c...        No continuous neighbouring segment found, start here:
            irstart = ir1
            EXIT
          ENDIF
        ENDDO
        IF (irstart.EQ.0) 
     .    CALL ER('BuildNeutralWall','Unable to find starting '//
     .            'point',*99)
      ELSE
        irstart = 0
        DO ir = irsep, irwall-1
          IF (idring(ir).EQ.TARTOWAL.OR.idring(ir).EQ.TARTOTAR)
     .      irstart = ir
        ENDDO
      ENDIF

      write(pinout,*) 'IRSTART:',irstart

      DO i1 = 1, walln-1
        IF (SQRT((wallr1(i1,1)-wallr1(i1,2))**2 + 
     .           (wallz1(i1,1)-wallz1(i1,2))**2).LT.1.0E-6) 
     .   CALL WN('BuildNeutralWall','Very short segment detected')
      ENDDO

      IF (cgridopt.EQ.LINEAR_GRID.OR.cgridopt.EQ.RIBBON_GRID) THEN
c...    For linear grids the 'wall' array does not go all the way around
c       the simulation domain/volume since the center of the plasma is at
c       R=0.0:
        rstart = rvertp(1,korpg(1,1))
        zstart = zvertp(1,korpg(1,1))
      ELSE
        rstart = rvertp(2,korpg(1,irstart))
        zstart = zvertp(2,korpg(1,irstart))
      ENDIF

      ! jdemod 
      write(pinout,'(a)') 'WALL BEFORE SEQUENCING:'
      do i1 = 1,walln
         write(pinout,'(i6,4(1x,g12.5))') i1,wallr1(i1,1),wallz1(i1,1),
     >                                 wallr1(i1,2),wallz1(i1,2)
      end do


      write(pinout,*) 'Rstart,zstart:',rstart,zstart

      i2 = 0
      DO i1 = 1, walln-1

        IF (i1.EQ.1) THEN
          WRITE(pinout,*) '   SEARCH:',i1,walln,rstart,zstart,irstart
        ELSE
          WRITE(pinout,*) '   SEARCH:',i1-1,walln,wallr1(MAX(1,i1-1),2),
     .                                 wallz1(MAX(1,i1-1),2)
        ENDIF

        DO i2 = 1, 2
          DO i3 = i1, walln
            IF ((i2.EQ.1.AND.wallc(i3).LE.0).OR.       ! check target segments on first pass
     .          (i2.EQ.2.AND.wallc(i3).GT.0)) CYCLE


        IF (i1.EQ.1) THEN
          WRITE(pinout,*) '  COMPARE:',i1,i3,
     >            ABS(wallr1(i3,1)-rstart),ABS(wallz1(i3,1)-zstart),
     >           ABS(wallr1(i3,1)-rstart               ).LT.TOL.AND. 
     .           ABS(wallz1(i3,1)-zstart               ).LT.TOL
        ELSE
          WRITE(pinout,*) '  COMPARE:',i1-1,i3,
     .           ABS(wallr1(i3,1)-wallr1(MAX(1,i1-1),2)),
     .           ABS(wallz1(i3,1)-wallz1(MAX(1,i1-1),2)),
     .           ABS(wallr1(i3,1)-wallr1(MAX(1,i1-1),2)).LT.TOL.AND.
     .           ABS(wallz1(i3,1)-wallz1(MAX(1,i1-1),2)).LT.TOL

        ENDIF


            IF ((i1.EQ.1.AND.
     .           ABS(wallr1(i3,1)-rstart               ).LT.TOL.AND. 
     .           ABS(wallz1(i3,1)-zstart               ).LT.TOL).OR.
     .          (i1.GT.1.AND.
     .           ABS(wallr1(i3,1)-wallr1(MAX(1,i1-1),2)).LT.TOL.AND.
     .           ABS(wallz1(i3,1)-wallz1(MAX(1,i1-1),2)).LT.TOL)) THEN
	  
              WRITE(pinout,*) '    MATCH:',i3,wallr1(i3,1),wallz1(i3,1)
              WRITE(pinout,*) '         :',i3,wallr1(i3,2),wallz1(i3,2)
	  
              holdr1 = wallr1(i1,1) 
              holdz1 = wallz1(i1,1) 
              holdr2 = wallr1(i1,2) 
              holdz2 = wallz1(i1,2) 
              holdc  = wallc (i1)
              holdt  = wallt (i1)
	  
              wallr1(i1,1) = wallr1(i3,1)
              wallz1(i1,1) = wallz1(i3,1)
              wallr1(i1,2) = wallr1(i3,2)
              wallz1(i1,2) = wallz1(i3,2)
              wallc (i1)   = wallc (i3)
              wallt (i1)   = wallt (i3)
	  
              wallr1(i3,1) =  holdr1
	      wallz1(i3,1) =  holdz1
	      wallr1(i3,2) =  holdr2
	      wallz1(i3,2) =  holdz2
              wallc (i3)   =  holdc
	      wallt (i3)   =  holdt
 	  
              EXIT
            ENDIF
	  
          ENDDO
          IF (i3.LT.walln+1) EXIT

        ENDDO

        IF     (i2.EQ.3) THEN
c...      The next wall segment could not be found:
          CALL ER('BuildNeutralWall','Unable to sequence wall',*99)
        ELSEIF (wallr1(1,1).EQ.wallr1(i1,2).AND.
     .          wallz1(1,1).EQ.wallz1(i1,2)) THEN
c...      Wall is closed:
          walln = i1
          EXIT
        ENDIF

      ENDDO

c
c     jdemod
c     - close the wall specification by connecting the first to the 
c       last wall points - the specification given here leaves out the
c       core boundary on ribbon grids. 
c     - I'm not sure if linear grids need the same fix
c
      if (cgridopt.EQ.LINEAR_GRID.OR.cgridopt.eq.RIBBON_GRID) then 
         walln = walln+1
         wallr1(walln,1) = wallr1(walln-1,2)
         wallz1(walln,1) = wallz1(walln-1,2)
         wallr1(walln,2) = wallr1(1,1)
         wallz1(walln,2) = wallz1(1,1)
         wallc(walln) = 0
         wallt(walln) = 0
      endif

c
c
c
      pcnt    = 0
      wallpts = 0

      ! jdemod - sizes wrong in calls to rzero and with new fortran syntax it is
      !          safer and easier to properly initialize the arrays with an assignment
      wallpt = 0.0
      wallpt2 = 0.0
      rw = 0.0
      zw = 0.0
      !CALL RZero(wallpt ,MAXPTS*19)
      !CALL RZero(wallpt2,MAXPTS*2)
      !CALL RZero(rw     ,MAXPTS)
      !CALL RZero(zw     ,MAXPTS)

c...  Assign quantities in the WALLPT array (the principal wall data 
c     array).  This is usually done in the DOWALL routine:
c
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

      wallpts = walln

c      write(0,*) 'BUILDNEUTRALWALL:WALLN:',walln

      if (cprint.eq.3.or.cprint.eq.9) then 
        write(6,*) 'BUILDNEUTRALWALL:WALLN:',walln

         do in = 1,walln

c         write(0,'(a,3i8,10(1x,g18.8))') 'BNW:',in,
c     >             wallt(in),wallc(in),
c     >             wallr1(in,1),wallz1(in,1),
c     >             wallr1(in,2),wallz1(in,2)
            write(6,'(a,3i8,10(1x,g18.8))') 'BNW:',in,
     >             wallt(in),wallc(in),
     >             wallr1(in,1),wallz1(in,1),
     >             wallr1(in,2),wallz1(in,2)
         end do
      ENDIF
c
      DO in = 1, walln
        r1 = wallr1(in,1)
        z1 = wallz1(in,1)
        r3 = wallr1(in,2)
        z3 = wallz1(in,2)
        r2 = 0.5 * (r1 + r3)
        z2 = 0.5 * (z1 + z3)

        wallpt(in,1) = r2
        wallpt(in,2) = z2

        wallpt2(in,1) = r1
        wallpt2(in,2) = z1

        wallpt(in,3) = 1.0
        wallpt(in,4) = 1.0
        wallpt(in,5) = SQRT((r1 - r2)**2 + (z1 - z2)**2)
        wallpt(in,6) = SQRT((r2 - r3)**2 + (z2 - z3)**2)
        wallpt(in,7) = wallpt(in,5) + wallpt(in,6)
        wallpt(in,8) = ATAN2C(REAL(z1-z2),REAL(r1-r2))
        wallpt(in,9) = ATAN2C(REAL(z3-z2),REAL(r3-r2))

c...  WALLPT(x,16) (limited assignment):
c
c   * 'OT'  Outer Target        1
c     'OC'  Outer Corner        2
c     'OD'  Outer Divertor Wall 3   (Vessel wall roughly below X-point I think)
c   * 'IT'  Inner Target        4
c     'IC'  Inner Corner        5
c     'ID'  Inner Divertor Wall 6
c   * 'MS'  Main Vessel Wall    7
c   * 'PV'  Private Plasma Wall 8
c     'CB'  Compound Baffle     9
c     'BA'  Baffle Segment      10
c           Other               0
c
        wallpt(in,16)= wallt(in)

        IF (wallc(in).GT.0) THEN
          IF (wallt(in).EQ.1) THEN
            wallpt(in,18) = idds(wallc(in),2)
          ELSE
            wallpt(in,18) = idds(wallc(in),1)
          ENDIF
        ELSE
          wallpt(in,18)= 0.0
        ENDIF

        wallpt(in,20) = r1
        wallpt(in,21) = z1
        wallpt(in,22) = r3
        wallpt(in,23) = z3

        IF (wallt(in).EQ.0) THEN
c...      This block of code from DOWALL. Apply to all non-target wall 
c         segments:
          CALL Find_Nearest_Boundary_Cell(ikn,irn,r2,z2,dist,ri,zi)  
          IF     (ctrap.EQ.7) THEN
            WRITE(0,*) 'PROBLEM: CODE DEVELOPMENT REQUIRED IN '//
     .                 'BUILDNEUTRALWALL IF TRAP SEGMENTS ARE '//
     .                 'TO BE IDENTIFIED. HALTING DIVIMP.'
            STOP
          ELSEIF (cneur.NE.7) then 
c...        Record nearest plasma cell:
            wallpt(in,26) = REAL(ikn)
            wallpt(in,27) = REAL(irn)
          ENDIF
c...      Save distance and intersection point - intersection is temporary(eh?):
          wallpt(in,28) = dist
          wallpt(in,29) = ri
          wallpt(in,30) = zi
        ENDIF

        IF (.TRUE.) THEN  
c...      Assign RW and ZW:
          IF (pcnt.GT.MAXPTS-3) 
     .      CALL ER('BuildNeutralWall','Too many wall points, '//
     .              'increase MAXNKS',*98)
          pcnt = pcnt + 1
          rw(pcnt) = r1
          zw(pcnt) = z1

          if (in.eq.walln) then 
             ! jdemod - close the figure describing the neutral wall
             pcnt = pcnt + 1
             rw(pcnt) = r3
             zw(pcnt) = z3
             pcnt = pcnt + 1
             rw(pcnt) = wallpt(1,20)
             zw(pcnt) = wallpt(1,21)
          endif
c
c         This is a little different than in DOWALL... make sure it won't
c         cause a problem...  
c           The above is an old comment, and I don't recall what the difference is,
c           but not too worried since RW and ZW don't appear to be used in anger 
c           anywhere. -SL, 30/09/2010
c
c         jdemod
c           -rw and zw are the points that define the boundary for neutral transport
c           -these are used extensively in neut and neutone in calls to GA15A
c           -However, all that is required is that the points included should
c            properly specify the neutral bounding surface. If wallpts has been set
c            up correctly then assigning just r1,z1 and closing the loop should be sufficient
c
c
c         jdemod - remove additional duplicate vertices
c          IF (wallt(in).EQ.1.OR.wallt(in).EQ.4) THEN
c            pcnt = pcnt + 1
c            rw(pcnt) = r2
c            zw(pcnt) = z2
c          ENDIF
c          IF (wallt(in).NE.wallt(in+1)) THEN
c            pcnt = pcnt + 1
c            rw(pcnt) = r3
c            zw(pcnt) = z3
c          ENDIF
c
        ENDIF

      ENDDO


      IF (pcnt.GT.MAXPTS) 
     .  CALL ER('BuildNeutralWall','RZ,ZW out of bounds. Increase '//
     .          'MAXPTS.',*99)

c...  Assign WLWALL1,2 and WLTRAP1,2:
      wlwall1 = 0
      wlwall2 = 0
      wltrap1 = 0
      wltrap2 = 0
      DO i1 = 1, wallpts
        in = MAX(1,NINT(wallpt(i1,18)))
        ik = ikds(in)
        ir = irds(in)
        IF (iflexopt(8).EQ.11) THEN
c...      980116023 -- assumption here is that WALLPT starts with the 
c         neutral wall near IRWALL on the high field side:
          IF (wlwall1.EQ.0.AND.wallpt(i1,18).EQ.0.0) wlwall1 = i1
          IF (wlwall2.EQ.0.AND.wallpt(i1,18).NE.0.0) wlwall2 = i1 - 1
          IF (wltrap1.EQ.0.AND.wlwall2.NE.0.AND.
     .                         wallpt(i1,18).EQ.0.0) wltrap1 = i1
          IF (wltrap2.EQ.0.AND.wltrap1.NE.0.AND.
     .                         wallpt(i1,18).NE.0.0) wltrap2 = i1 - 1
        ELSE
          IF (wlwall1.EQ.0.AND.ik.EQ.1      .AND.
     .                         ir.EQ.irwall-1) wlwall1 = i1 + 1
          IF (wlwall1.NE.0.AND.ik.EQ.nks(ir).AND.
     .                         ir.EQ.irwall+1) wlwall2 = i1 - 1
          IF (wltrap1.EQ.0.AND.ik.EQ.nks(ir).AND.
     .                         ir.EQ.irtrap+1) wltrap1 = i1 + 1
          IF (wltrap1.NE.0.AND.ik.EQ.1      .AND.
     .                         ir.EQ.irtrap+1) wltrap2 = i1 - 1
        ENDIF
      ENDDO
c...  Blank WLWALL1,2 for a generalized geometry, since the main chamber
c     wall is no longer well defined:
      IF (grdnmod.NE.0.AND.iflexopt(8).NE.11) THEN
c...                       980116023:
        wlwall1 = 0
        wlwall2 = 0
      ENDIF
      IF (cgridopt.EQ.LINEAR_GRID.OR.irtrap.GT.nrs.OR.
     .    cgridopt.EQ.RIBBON_GRID) THEN      
c...    No PFZ here:
        IF (wallpts.LT.MAXPTS) THEN
          wltrap1 = wallpts + 1
          wltrap2 = wallpts + 1
        ELSE
          CALL ER('BuildNeutralWall','Sorry, need bigger MAXPTS',*99)
        ENDIF
      ENDIF
 
c      IF (sloutput) 
c     .  WRITE(0,*) 'DONE BUILDING NEUTRAL WALL -- NEED TO ASSIGN '//
c     .             'MORE WALLPT ENTRIES'

c...  Assign WALLINDEX:
      DO i1 = 1, wallpts
        in = NINT(wallpt(i1,18))
        IF (in.NE.0) wallindex(in) = i1
      ENDDO


c
c
c      write(6,'(a,2i8)') 'Wallpts:',wallpts,pcnt
c      do i1 = 1,wallpts
c         write(6,'(10(1x,g18.8))') wallpt(i1,20),wallpt(i1,21),
c     >              wallpt(i1,22),wallpt(i1,23),rw(i1),zw(i1)
c      end do




c      CALL OutputData(85,'DONE BUILDING NEUTRAL WALL')
c      STOP 'sdfsddsf'

c      IF (terminate) THEN
c        CALL OutputData(85,'Checking neutral wall')
c        CALL DumpGrid('BUILDING WALL SEGMENT LIST')
c      ENDIF


      RETURN
98    WRITE(0,*) '  MAXPTS = ',MAXPTS
      WRITE(0,*) '  NWALL  = ',walln
      WRITE(0,*) '  GUESS  = ',walln-COUNT(wallt.NE.0)+
     .                         2*COUNT(wallt.NE.0)
      CALL DumpGrid('TOO MANY NEUTRAL WALL POINTS')
99    CONTINUE
      WRITE(0,*) '   I1   :',i1
      WRITE(0,*) '   R1,Z1:',wallr1(i1,2),wallz1(i1,2)
      WRITE(pinout,*)
      WRITE(pinout,*) '   R1,Z1:',wallr1(i1,2),wallz1(i1,2)
      DO i1 = 1, walln
        WRITE(pinout,*)  i1,(wallr1(i1,i2),wallz1(i1,i2),i2=1,2)
      ENDDO
      nvesm = i1-1
      DO i2 = 1, i1-1
        jvesm(i2) = 8
        rvesm(i2,1) = wallr1(i2,1)
        zvesm(i2,1) = wallz1(i2,1)
        rvesm(i2,2) = wallr1(i2,2)
        zvesm(i2,2) = wallz1(i2,2)
      ENDDO 
      CALL DumpGrid('UNABLE TO SEQUENCE NEUTRAL WALL')
      STOP
      END
c
c ======================================================================
c
c subroutine: BuildGridPolygons
c
c
c
      SUBROUTINE BuildGridPolygons
      use error_handling
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'grbound'
      INCLUDE 'slcom'

      REAL*8, PARAMETER :: TOL = 3.0E-7

      INTEGER id,walln,i1,i2,ik,ir,kind,fp
      REAL    holdr(2),holdz(2),wallr1(MAXPTS+1,2),wallz1(MAXPTS+1,2)


      IF (cionr.NE.2) 
     .  CALL ER('BuildGridPolygons','CIONR=2 is this only support '//
     .          'ion wall option for generalized grids',*99)

      fp = PINOUT

      walln = 0


c...  Add target segments:
      DO ir = irsep, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        walln = walln + 1
        if (walln.gt.maxpts) then 
           ! jdemod - add some bounds checking
           call errmsg('BuildGridPolygons','NWALL > MAXPTS')
           write(0,*) 'NWALL = ',nwall,' MAXPTS = ',maxpts
           stop 'NWALL > MAXPTS'
        endif
        id = korpg(1,ir)
        wallr1(walln,1) = rvertp(1,id)
        wallz1(walln,1) = zvertp(1,id)
        wallr1(walln,2) = rvertp(2,id)
        wallz1(walln,2) = zvertp(2,id)
        walln = walln + 1
        if (walln.gt.maxpts) then 
           ! jdemod - add some bounds checking
           call errmsg('BuildGridPolygons','NWALL > MAXPTS')
           stop 'NWALL > MAXPTS'
           write(0,*) 'NWALL = ',nwall,' MAXPTS = ',maxpts
        endif
        id = korpg(nks(ir),ir)
        wallr1(walln,1) = rvertp(3,id)
        wallz1(walln,1) = zvertp(3,id)
        wallr1(walln,2) = rvertp(4,id)
        wallz1(walln,2) = zvertp(4,id)
      ENDDO
c...  Add IRWALL boundary ring segments:
      ir = irwall
      DO ik = 1, nks(ir)
        id = korpg(ikins(ik,ir),irins(ik,ir))
        walln = walln + 1
        if (walln.gt.maxpts) then 
           ! jdemod - add some bounds checking
           call errmsg('BuildGridPolygons','NWALL > MAXPTS')
           stop 'NWALL > MAXPTS'
        endif
        IF (irouts(ikins(ik,ir),irins(ik,ir)).EQ.irwall) THEN
          wallr1(walln,1) = rvertp(2,id)
          wallz1(walln,1) = zvertp(2,id)
          wallr1(walln,2) = rvertp(3,id)
          wallz1(walln,2) = zvertp(3,id)
        ELSE
c          WRITE(0,*) '----> ',ikins(ik,ir),irins(ik,ir)
          wallr1(walln,1) = rvertp(4,id)
          wallz1(walln,1) = zvertp(4,id)
          wallr1(walln,2) = rvertp(1,id)
          wallz1(walln,2) = zvertp(1,id)
        ENDIF
      ENDDO
c...  Add IRTRAP boundary ring segments:
      IF (cgridopt.NE.LINEAR_GRID.AND.cgridopt.NE.RIBBON_GRID) THEN
        ir = irtrap
        DO ik = 1, nks(ir)
          id = korpg(ikouts(ik,ir),irouts(ik,ir))
          walln = walln + 1
        if (walln.gt.maxpts) then 
           ! jdemod - add some bounds checking
           call errmsg('BuildGridPolygons','NWALL > MAXPTS')
           stop 'NWALL > MAXPTS'
        endif
          wallr1(walln,1) = rvertp(4,id)
          wallz1(walln,1) = zvertp(4,id)
          wallr1(walln,2) = rvertp(1,id)
          wallz1(walln,2) = zvertp(1,id)
        ENDDO
      ENDIF

c      WRITE(0,*) 'WALLN=',walln
      WRITE(6,*) 'WALLN=',walln

      DO i1 = 2, walln-1
        DO i2 = i1, walln
          IF (ABS(wallr1(i2,1)-wallr1(i1-1,2)).LT.TOL.AND.
     .        ABS(wallz1(i2,1)-wallz1(i1-1,2)).LT.TOL) THEN
c...        Need to check that they are exactly the same, and if not
c           then make them exactly the same:  
            holdr(    1:2) = wallr1(i1,1:2) 
            holdz(    1:2) = wallz1(i1,1:2) 
            wallr1(i1,1:2) = wallr1(i2,1:2)
            wallz1(i1,1:2) = wallz1(i2,1:2)
            wallr1(i2,1:2) = holdr(    1:2)
            wallz1(i2,1:2) = holdz(    1:2)
            EXIT
          ENDIF
        ENDDO
c...    The next wall segment could not be found:
        IF (i2.EQ.walln+1) 
     .    CALL ER('BuildGridPolygons','Unable to polygon grid',*99)
      ENDDO

c...  Assign:
      ionwpts = walln + 1
      riw(1:ionwpts-1) = wallr1(1:walln,1)
      ziw(1:ionwpts-1) = wallz1(1:walln,1)
      riw(ionwpts) = riw(1)
      ziw(ionwpts) = ziw(1)


c      write(0,*) 'BGP:IONWPTS:',ionwpts
      write(6,*) 'BGP:IONWPTS:',ionwpts

      if (cprint.eq.3.or.cprint.eq.9) then 
         do i1 = 1,ionwpts
            !write(0,'(a,i8,10(1x,g18.8))') 'IONW:',i1,riw(i1),ziw(i1)
            write(6,'(a,i8,10(1x,g18.8))') 'IONW:',i1,riw(i1),ziw(i1)
         end do
      endif

c...  Not really sure what this does, but found in IONWALL in WALLS.F, 
c     seems to setup some work arrays:
      kind = 1
      CALL GA15A(IONWPTS,KIND,iwWORK,4*MAXPTS,iwINDW,MAXPTS,
     >           RIW,ZIW,iwTDUM,iwXDUM,iwYDUM,6)
      WRITE(SLOUT,'(A,I10)') 'RETURN FROM GA15A, TAU =',iwINDW(2,1)

c...  Core boundary polygon:
      !
      ! jdemod
      ! the description used for the core boundary would seem to leave out 
      ! the edge of the last polygon - this implicitly assumes that the 
      ! first and last cells of a core ring are the same - which is not the 
      ! case for a ribbon grid. 
      ! ALSO - use of ring #2 here is also based on standard grid assumptions
      !        though this works for a ribbon grid. 
      !
      ioncpts = nks(2) 
      DO ik = 1, nks(2) - 1
        id = korpg(ik,2)
        rcw(ik) = rvertp(1,id)
        zcw(ik) = zvertp(1,id)
      enddo
      if (cgridopt.eq.RIBBON_GRID) then 
         ! add the information for the last polygon - both vertices 
         ! for a ribbon grid
         ioncpts = ioncpts + 1
        id = korpg(nks(ir),2)
        rcw(ioncpts) = rvertp(1,id)
        zcw(ioncpts) = zvertp(1,id)
         ioncpts = ioncpts + 1
        rcw(ioncpts) = rvertp(4,id)
        zcw(ioncpts) = zvertp(4,id)
      endif

      ! Close the figure
      rcw(ioncpts) = rcw(1)
      zcw(ioncpts) = zcw(1)


c      write(0,*) 'BGP:IONCPTS:',ioncpts
      write(6,*) 'BGP:IONCPTS:',ioncpts

      if (cprint.eq.3.or.cprint.eq.9) then 
         do i1 = 1,ioncpts
            !write(0,'(a,i8,10(1x,g18.8))') 'CORW:',i1,rcw(i1),zcw(i1)
            write(6,'(a,i8,10(1x,g18.8))') 'CORW:',i1,rcw(i1),zcw(i1)
         end do
      endif

      CALL GA15A(IONCPTS,KIND,icWORK,4*MAXPTS,icINDW,MAXPTS,
     >             RCW,ZCW,icTDUM,icXDUM,icYDUM,6)

c...  Some stuff from the bottom of the IONWALL routine in WALLS.F:
      if     (xygrid.eq.1) then
        CALL ER('BuildGridPolygons','Sorry, XYGRID=1 not supported '//
     .          'at the moment',*99)
      elseif (xygrid.eq.0) then
         DR = (RMAX-RMIN) / REAL(maxgxs-4)
         DZ = (ZMAX-ZMIN) / REAL(maxgys-4)
         RMIN = RMIN - 2.0 * DR
         RMAX = RMAX + 2.0 * DR
         ZMIN = ZMIN - 2.0 * DZ
         ZMAX = ZMAX + 2.0 * DZ
      endif

      RETURN
99    WRITE(fp,*) '  I1,2=',i1,i2
      WRITE(fp,*) '  WALL1:',wallr1(i1-1,2),wallz1(i1-1,2)
      WRITE(fp,*) '  WALL2:'
      DO i2 = 2, walln
        WRITE(fp,'(1X,A,I4,4F16.8,4X,2F16.8)')
     .     '       :',i2,wallr1(i2,1),wallz1(i2,1),
     .                   wallr1(i2,2),wallz1(i2,2),
     .     wallr1(i2,1)-wallr1(i2-1,2),
     .     wallz1(i2,1)-wallz1(i2-1,2)
      ENDDO
      CALL OutputData(85,'TRYING TO POLYGON THE GRID')
      CALL DumpGrid('TRYING TO POLYGON THE GRID')
      STOP
      END
c
c ======================================================================
c ======================================================================
c
c block: Low level grid manipulation
c
c COPYRING
c INSERTRING
c DELETERING
c MOVECELL
c
c ======================================================================
c
c subroutine: MergeRings
c
c 
c
c
      SUBROUTINE MergeRings(ir)
      USE mod_grid_divimp
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER ir

      INTEGER ik,id1,id2,i1

c...  Check:
      IF (nks(ir).NE.nks(ir+1)) 
     .  CALL ER('MergeRings','Rings to be merged do not have the '//
     .          'same numbers of cells',*99)

      DO ik = 1, nks(ir)
c...    Only a small number of cell quantities need to be specified, 
c       in the ends of the rings are to be cut shortly:

        id1 = korpg(ik,ir)
        id2 = korpg(ik,ir+1)

        rvertp(2,id1)=rvertp(2,id2)
        zvertp(2,id1)=zvertp(2,id2)
        rvertp(3,id1)=rvertp(3,id2)
        zvertp(3,id1)=zvertp(3,id2)
        IF (ALLOCATED(d_rvertp)) THEN
          d_rvertp(2,id1) = d_rvertp(2,id2)
          d_zvertp(2,id1) = d_zvertp(2,id2)
          d_rvertp(3,id1) = d_rvertp(3,id2)
          d_zvertp(3,id1) = d_zvertp(3,id2)
        ENDIF 
        rs(ik,ir) = 0.0
        zs(ik,ir) = 0.0
        DO i1 = 1, 4
          rs(ik,ir+1) = rs(ik,ir+1) + 0.25 * rvertp(i1,id2)
          zs(ik,ir+1) = zs(ik,ir+1) + 0.25 * zvertp(i1,id2)
        ENDDO
      ENDDO

      CALL DeleteRing(ir+1)

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: AddOuterRing
c
c 
c
c
      SUBROUTINE AddOuterRing(ir,frac1)
      USE mod_grid_divimp
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER ir
      REAL    frac1

      INTEGER ik,id1,id2,i1
      REAL    frac

      WRITE(0,*) 'DATA:',ir,irsep,irwall,frac1

      IF (ir.GE.irsep.AND.ir.LT.irwall) THEN

c...    Adding a ring outside the ..
c       
        CALL InsertRing(ir,AFTER,PERMANENT)

        WRITE(0,*) 'DATA:',ir,irsep,irwall,frac1

c...    Setup geometric quantities for the new ring:
        nks   (ir+1) = nks(ir)
        irorg2(ir+1) = ir
        idring(ir+1) = TARTOTAR 

        frac = frac1 + 1.0

        DO ik = 1, nks(ir+1)

c...      Copy cell geometry data for the ring.  Only a limited
c         set of cell/ring data needs to be specified because
c         the ring is being added shortly after the grid is read in,
c         and before it is processed by OEDGE:

          bratio(ik,ir+1) = bratio(ik,ir)
          kbfs  (ik,ir+1) = kbfs  (ik,ir)                
          rs    (ik,ir+1) = rs    (ik,ir)
          zs    (ik,ir+1) = zs    (ik,ir)

          id1 = korpg(ik,ir)
          id2 = korpg(ik,ir+1)
          nvertp(id2) = 4           ! Who are we kidding, it will always be 4...

          rvertp(1,id2)=rvertp(2,id1)
          zvertp(1,id2)=zvertp(2,id1)
          rvertp(2,id2)=rvertp(1,id1)+frac*(rvertp(2,id1)-rvertp(1,id1))
          zvertp(2,id2)=zvertp(1,id1)+frac*(zvertp(2,id1)-zvertp(1,id1))
          rvertp(3,id2)=rvertp(4,id1)+frac*(rvertp(3,id1)-rvertp(4,id1))
          zvertp(3,id2)=zvertp(4,id1)+frac*(zvertp(3,id1)-zvertp(4,id1))
          rvertp(4,id2)=rvertp(3,id1)
          zvertp(4,id2)=zvertp(3,id1)
          IF (ALLOCATED(d_rvertp)) THEN
            d_rvertp(1,id2)=d_rvertp(2,id1)
            d_zvertp(1,id2)=d_zvertp(2,id1)
            d_rvertp(2,id2)=d_rvertp(1,id1)+DBLE(frac)*(d_rvertp(2,id1)-
     .                                                  d_rvertp(1,id1))
            d_zvertp(2,id2)=d_zvertp(1,id1)+DBLE(frac)*(d_zvertp(2,id1)-
     .                                                  d_zvertp(1,id1))
            d_rvertp(3,id2)=d_rvertp(4,id1)+DBLE(frac)*(d_rvertp(3,id1)-
     .                                                  d_rvertp(4,id1))
            d_zvertp(3,id2)=d_zvertp(4,id1)+DBLE(frac)*(d_zvertp(3,id1)-
     .                                                  d_zvertp(4,id1))
            d_rvertp(4,id2)=d_rvertp(3,id1)
            d_zvertp(4,id2)=d_zvertp(3,id1)
          ENDIF

          rs(ik,ir+1) = 0.0
          zs(ik,ir+1) = 0.0
          DO i1 = 1, 4
            rs(ik,ir+1) = rs(ik,ir+1) + 0.25 * rvertp(i1,id2)
            zs(ik,ir+1) = zs(ik,ir+1) + 0.25 * zvertp(i1,id2)
          ENDDO
        ENDDO
      ELSE
        CALL ER('AddOuterRing','Core and PFZ cut require work',*99)
      ENDIF

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: DupeRing
c
c 
c
c
      SUBROUTINE DupeRing(ir)
      USE mod_grid_divimp
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER ir

      INTEGER ik,id1,id2,i1,irset

      IF (idring(ir).EQ.BOUNDARY) 
     .  CALL ER('TailorGrid','Trying to duplicate a boundary ring',*99)

      IF (ir.GE.irsep.AND.ir.LE.nrs) THEN
c      IF (ir.GE.irsep.AND.ir.LT.irwall) THEN
c...    This ring is to be cut in the middle, so add another ring to the grid
c       at that will serve as the outer half of the ring being cut:
        IF (ir.LT.irwall) THEN
          CALL InsertRing(irwall,BEFORE,PERMANENT)
          irset = irwall - 1
        ELSE
          CALL InsertRing(nrs,AFTER,PERMANENT)
          irset = nrs
        ENDIF

c...    Setup geometric quantities for the new ring:
        nks   (irset) = nks(ir)
        irorg2(irset) = ir
        idring(irset) = TARTOTAR - 100
        psitarg(irset,1:2) = psitarg(ir,1:2)

        DO ik = 1, nks(ir)  
c...      Copy cell geometry data for the ring.  Only a limited
c         set of cell/ring data needs to be specified because
c         the ring is being added shortly after the grid is read in,
c         and before it is processed by OEDGE:
          bratio(ik,irset) = bratio(ik,ir)
          kbfs  (ik,irset) = kbfs  (ik,ir)                
          rs    (ik,irset) = rs    (ik,ir)
          zs    (ik,irset) = zs    (ik,ir)

          id1 = korpg(ik,ir)
          id2 = korpg(ik,irset)
          nvertp(id2) = nvertp(id1)
          DO i1 = 1, nvertp(id2)
            rvertp(i1,id2) = rvertp(i1,id1)
            zvertp(i1,id2) = zvertp(i1,id1)
            IF (ALLOCATED(d_rvertp)) THEN
              d_rvertp(i1,id2) = d_rvertp(i1,id1)
              d_zvertp(i1,id2) = d_zvertp(i1,id1)
            ENDIF
          ENDDO
        ENDDO
      ELSE 
        CALL ER('DupeRing','Core ring duplication requires work',*99)
      ENDIF

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: ExpandGrid
c
      SUBROUTINE ExpandGrid(ndupe,size_frac,ir_reference) 
      USE mod_grid_divimp
      IMPLICIT none
  
      INTEGER, INTENT(IN) :: ndupe,ir_reference
      REAL   , INTENT(IN) :: size_frac

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER ik,id,idupe,irset,irref
      REAL*8  frac

      IF (idring(ir_reference).EQ.BOUNDARY) 
     .  CALL ER('ExpandGrid','Trying to expand grid with a boundary '//
     .          'ring',*99)

      IF (.NOT.ALLOCATED(d_rvertp)) 
     .  CALL ER('ExpandGrid','Expecting double precision vertex '//
     .          'arrays to be allocated',*99)

      IF (size_frac.LE.0.0) 
     .  CALL ER('ExpandGrid','Scaing fraction .LE. 0.0',*99)


      irref = ir_reference
      frac  = DBLE(size_frac + 1.0)


      WRITE(0,*) 'FRAC:',frac

      DO idupe = 1, ndupe

        CALL DupeRing(irref)

        IF (irref.LT.irwall) THEN 
          irset = irwall - 1  ! The convention in DupeRing
          DO ik = 1, nks(irset) 
            id = korpg(ik,irset)
          
            d_rvertp(2,id) =         d_rvertp(1,id) + 
     .                       frac * (d_rvertp(2,id) - d_rvertp(1,id))
            d_zvertp(2,id) =         d_zvertp(1,id) + 
     .                       frac * (d_zvertp(2,id) - d_zvertp(1,id))
          
            d_rvertp(3,id) =         d_rvertp(4,id) + 
     .                       frac * (d_rvertp(3,id) - d_rvertp(4,id))
            d_zvertp(3,id) =         d_zvertp(4,id) + 
     .                       frac * (d_zvertp(3,id) - d_zvertp(4,id))
          
            d_rvertp(1,id) = d_rvertp(2,korpg(ik,irref))
            d_zvertp(1,id) = d_zvertp(2,korpg(ik,irref))
            d_rvertp(4,id) = d_rvertp(3,korpg(ik,irref))
            d_zvertp(4,id) = d_zvertp(3,korpg(ik,irref))

            rvertp(1:4,id) = SNGL(d_rvertp(1:4,id))
            zvertp(1:4,id) = SNGL(d_zvertp(1:4,id))
          ENDDO
        ELSE
          STOP 'PFZ GRID EXTENSION NEEDS WORK'
        ENDIF

        idring(irset) = idring(irref)

        irref = irset

      ENDDO 

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: CopyRing
c
c Copy ring and cell quantities from one ring to another (existing) ring.
c
c
      SUBROUTINE CopyRing(ir1,ir2)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER ir1,ir2
      INTEGER ik

      IF (ir1.GT.MAXNRS.OR.ir2.GT.MAXNRS)
     .  CALL ER('InsertRing','Ring index is out of bounds',*99)

c...  Copy ring:
      nks   (ir2)   = nks   (ir1)
      irorg2(ir2)   = irorg2(ir1)
      idring(ir2)   = idring(ir1)
      ikto2 (ir2)   = ikto2 (ir1)
      ikti2 (ir2)   = ikti2 (ir1)
      idds  (ir2,1) = idds  (ir1,1)
      idds  (ir2,2) = idds  (ir1,2)
      ksmaxs(ir2)   = ksmaxs(ir1)
      psitarg(ir2,1) = psitarg(ir1,1)
      psitarg(ir2,2) = psitarg(ir1,2)
      DO ik = 1, nks(ir2)
        CALL MoveCell(ik,ir2,ik,ir1)
      ENDDO

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: InsertRing
c
c Target quantities are not adjusted, as it is assumed that this
c grid manipluation is done either at the beginning of the DIVIMP run,
c immediately after the geometry data has been read in, or in
c conjunction with ring deletion, so that the changes are removed before
c calling other DIVIMP routines (as would be the case when writing the
c EIRENE plasma file).
c
c
c
      SUBROUTINE InsertRing(irref,mode,type)
      USE mod_grid_divimp
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER FetchKORPG

      INTEGER irref,mode,type
      INTEGER ik,ir,nshift,irend,ir1,ik1

      nshift =  1
      ir     = -1
c
c     Checks:
c
      IF (irref.LT.1.OR.irref.GT.nrs)
     .  CALL ER('InsertRing','Reference ring is out of bounds',*99)

      IF (nrs+nshift.GT.MAXNRS)
     .  CALL ER('InsertRing','MAXNRS exceeded',*99)
c
c
c
      IF     (mode.EQ.BEFORE) THEN
        irend = irref + 1
      ELSEIF (mode.EQ.AFTER ) THEN
        irend = irref + 2
      ELSE
        CALL ER('InsertRing','Invalid mode',*99)
      ENDIF
c
c     Make room for new ring:
c
      DO ir = nrs+1, irend, -1
c Account for new virtual cell related quantities...
        nks   (ir)   = nks   (ir-nshift)
        irorg2(ir)   = irorg2(ir-nshift)
        idring(ir)   = idring(ir-nshift)
        ikto2 (ir)   = ikto2 (ir-nshift)
        ikti2 (ir)   = ikti2 (ir-nshift)
        idds  (ir,1) = idds  (ir-nshift,1)
        idds  (ir,2) = idds  (ir-nshift,2)
        ksmaxs(ir)   = ksmaxs(ir-nshift)
        psitarg(ir,1) = psitarg(ir-nshift,1)
        psitarg(ir,2) = psitarg(ir-nshift,2)
        DO ik = 1, nks(ir)
          CALL MoveCell(ik,ir,ik,ir-nshift)
        ENDDO
      ENDDO
c
c     Setup ring and cell quantities on new ring:
c
      ir = irend - 1

      nks   (ir)   = nks (irref)
c...bug : Allowed non-zero target fluxes for virtual rings. - Nov 30, 1999
      idds  (ir,1) = MAXNDS
      idds  (ir,2) = MAXNDS
c      idds  (ir,1) = idds(irref,1)
c      idds  (ir,2) = idds(irref,2)

      irorg2(ir)   = -1
c      idring(ir)   = -1
      idring(ir) = idring(irref)
      ksmaxs(ir)   = REAL(nks(ir))

      ksb(0,ir) = 0.0

      knds(idds(ir,1)) = 1.0E+12
      knds(idds(ir,2)) = 1.0E+12
      kteds(idds(ir,1)) = 5.0
      kteds(idds(ir,2)) = 5.0
      ktids(idds(ir,1)) = 5.0
      ktids(idds(ir,2)) = 5.0

c...BUG! Of sorts.  This makes sure that enough pointers are setup for each ring:
      DO ik = 1, MAXNKS
c      DO ik = 1, nks(ir)

c...    Assign a KORPG index:
        IF (type.EQ.PERMANENT) THEN
          korpg(ik,ir) = FetchKORPG(ik,ir)
        ELSE
          korpg(ik,ir) = MAXNRS*MAXNKS
        ENDIF

        nvertp(korpg(ik,ir)) = 0
        rvertp(:,korpg(ik,ir)) = 0.0  ! BUG Need the assignment of getting random values, which caused problem
        zvertp(:,korpg(ik,ir)) = 0.0  ! for grid grid_iter_10d -SL, 22/02/2012
c
c       jdemod - assigning arbitrary values of -1.0 to the cell centers 
c                causes problems with various pieces of code elsewhere
c                e.g. rmin,rmax are calculated as min(rs(ik,ir)) and max(rs(ik,ir)) which
c                     fails if one of these has been set to -1.0
c                Alternate idea - set rs(ik,ir) and zs(ik,ir) to the values in the 
c                reference ring. 
c
c        rs    (ik,ir) = -1.0
c        zs    (ik,ir) = -1.0
c
        rs    (ik,ir) = rs (ik,irref)
        zs    (ik,ir) = zs (ik,irref) 
c
        kbfs  (ik,ir) =  1.0
        bratio(ik,ir) =  1.0
        knbs  (ik,ir) =  1.0E+12
        knes  (ik,ir) =  1.0E+12
        kvhs  (ik,ir) =  0.0
        ktibs (ik,ir) =  5.0
        ktebs (ik,ir) =  5.0
        kss   (ik,ir) =  REAL(ik) - 0.5
        ksb   (ik,ir) =  REAL(ik)

        IF (ALLOCATED(divimp_ik)) THEN
          divimp_ik(ik,ir) = divimp_ik(ik,irref)
          divimp_ir(ik,ir) = divimp_ir(ik,irref)
        ENDIF
      ENDDO
c
c     Update global grid parameters:
c
      IF (mode.EQ.BEFORE) THEN
        IF (irref.LE.irsep  ) irsep   = irsep   + nshift
        IF (irref.LE.irsep2 ) irsep2  = irsep2  + nshift
        IF (stopopt.EQ.121) STOP 'PFZ BREAK B'
        IF (irref.LE.irbreak) irbreak = irbreak + nshift
        IF (irref.LE.irwall ) irwall  = irwall  + nshift
        DO ir = irref, nrs + nshift
          IF (irorg2(ir).GT.0) irorg2(ir) = irorg2(ir) + nshift
        ENDDO
      ELSE
        IF (irref.LT.irsep  ) irsep   = irsep   + nshift
        IF (irref.LT.irsep2 ) irsep2  = irsep2  + nshift
c        IF (stopopt.EQ.121) STOP 'PFZ BREAK C'
        IF (irref.LT.irbreak) irbreak = irbreak + nshift
        IF (irref.LE.irwall ) irwall  = irwall  + nshift
        DO ir = irref+2, nrs + nshift
          IF (irorg2(ir).GT.0) irorg2(ir) = irorg2(ir) + nshift
        ENDDO
      ENDIF

      irtrap = irwall + nshift
      nrs = nrs + nshift

      irtrap2 = irtrap
      irwall2 = irwall

c      WRITE(0,*) 'INSERTRING:',irref,nrs

      RETURN
99    WRITE(0,*) 'IR   =',ir
      WRITE(0,*) 'IRREF=',irref,nrs
      STOP 'InsertRing'
      END
c
c ======================================================================
c
c subroutine: SplitRing
c 
c
c
      SUBROUTINE SplitRing(ir,sposition)
      USE mod_grid_divimp
      IMPLICIT none

      INTEGER ir
      REAL    sposition

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER ik,ii,id,in1,in2,ir1,ir2,ikmp1,ikmp2
      LOGICAL status
      REAL*8  r(6),z(6),spos,xval1(MAXNRS),xval2(MAXNRS),
     .                       yval1(MAXNRS),yval2(MAXNRS)


      IF (ir.EQ.1.OR.ir.EQ.irwall.OR.ir.EQ.irtrap)
     .  CALL ER('SplitRing','Trying to split a boundary ring',*99)

      CALL InsertRing(ir,AFTER,PERMANENT)

      spos = DBLE(sposition)

c...  Initialize the cells on the new ring:
      DO ik = 1, nks(ir+1)
         
c        IF (nvertp(korpg(ik,ir)).NE.4) 
c     .    CALL ER('SplitRing','Cell does not have 4 verticies',*99)

c...    Assume cells have 4 verticies:
        id = korpg(ik,ir)
        IF (ALLOCATED(d_rvertp)) THEN
          r(1:4) = d_rvertp(1:4,id)
          z(1:4) = d_zvertp(1:4,id)
        ELSE
          r(1:4) = rvertp(1:4,id)
          z(1:4) = zvertp(1:4,id)
        ENDIF
c        DO ii = 1, 4
c          r(ii) = rvertp(ii,id)
c          z(ii) = zvertp(ii,id)
c        ENDDO

        r(5) = r(1) + spos * (r(2) - r(1))
        z(5) = z(1) + spos * (z(2) - z(1))
        r(6) = r(4) + spos * (r(3) - r(4))
        z(6) = z(4) + spos * (z(3) - z(4))

c        WRITE(0,'(A,2I6,12F6.2)') 'lfsgk',ik,ir,(r(ii),z(ii),ii=1,6)

c...    Assign a KORPG index:
c        IF (npolyp+1.LT.vpolmin) THEN
c          npolyp = npolyp + 1
c          id     = npolyp
c        ELSE
c          CALL ER('SplitRing','Cell index out of bounds',*99)
c        ENDIF
c
c        korpg(ik,ir+1) = id

c...    Build the cell geometry:
        id = korpg(ik,ir+1)
        nvertp(id)   = 4
        IF (ALLOCATED(d_rvertp)) THEN
          d_rvertp(1,id) = r(5)
          d_zvertp(1,id) = z(5)
          d_rvertp(2,id) = r(2)
          d_zvertp(2,id) = z(2)
          d_rvertp(3,id) = r(3)
          d_zvertp(3,id) = z(3)
          d_rvertp(4,id) = r(6)
          d_zvertp(4,id) = z(6)      
        ENDIF
        rvertp(1,id) = SNGL(r(5))
        zvertp(1,id) = SNGL(z(5))
        rvertp(2,id) = SNGL(r(2))
        zvertp(2,id) = SNGL(z(2))
        rvertp(3,id) = SNGL(r(3))
        zvertp(3,id) = SNGL(z(3))
        rvertp(4,id) = SNGL(r(6))
        zvertp(4,id) = SNGL(z(6))      
        rs(ik,ir+1) = 0.25 * (rvertp(1,id) + rvertp(2,id) +
     .                        rvertp(3,id) + rvertp(4,id))
        zs(ik,ir+1) = 0.25 * (zvertp(1,id) + zvertp(2,id) +
     .                        zvertp(3,id) + zvertp(4,id))

c...    Should perhaps do something more complicated here, 
c       but this will do to start:
        bratio(ik,ir+1) = bratio(ik,ir)
        kbfs  (ik,ir+1) = kbfs  (ik,ir) 
      ENDDO

c...  Resize the cells on the pre-existing ring:
      DO ik = 1, nks(ir)
c...    Build the cell geometry:
        id = korpg(ik,ir)
        IF (ALLOCATED(d_rvertp)) THEN
          d_rvertp(2,id) = d_rvertp(1,korpg(ik,ir+1))
          d_zvertp(2,id) = d_zvertp(1,korpg(ik,ir+1))
          d_rvertp(3,id) = d_rvertp(4,korpg(ik,ir+1))
          d_zvertp(3,id) = d_zvertp(4,korpg(ik,ir+1))
        ENDIF
        rvertp(2,id) = rvertp(1,korpg(ik,ir+1))
        zvertp(2,id) = zvertp(1,korpg(ik,ir+1))
        rvertp(3,id) = rvertp(4,korpg(ik,ir+1))
        zvertp(3,id) = zvertp(4,korpg(ik,ir+1))
        rs(ik,ir) = 0.25 * (rvertp(1,id) + rvertp(2,id) +
     .                      rvertp(3,id) + rvertp(4,id))
        zs(ik,ir) = 0.25 * (zvertp(1,id) + zvertp(2,id) +
     .                      zvertp(3,id) + zvertp(4,id))
      ENDDO

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: DeleteRing
c
c
      SUBROUTINE DeleteRing(irref)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER irref, nshift
      INTEGER ik,ir

      nshift = 1
c
c     Check:
c
      IF (irref.LT.1.OR.irref.GT.nrs)
     .  CALL ER('DeleteRing','Ring to be deleted is outside grid',*99)
c
c     Move rings IRREF through NRS:
c
      DO ir = irref, nrs
        nks   (ir)   = nks   (ir+nshift)
        irorg2(ir)   = irorg2(ir+nshift)
        idring(ir)   = idring(ir+nshift)
        ikto2 (ir)   = ikto2 (ir+nshift)
        ikti2 (ir)   = ikti2 (ir+nshift)
        idds  (ir,1) = idds  (ir+nshift,1)
        idds  (ir,2) = idds  (ir+nshift,2)
        ksmaxs(ir)   = ksmaxs(ir+nshift)
        psitarg(ir,1) = psitarg(ir+nshift,1)
        psitarg(ir,2) = psitarg(ir+nshift,2)
        DO ik = 1, nks(ir)
          CALL MoveCell(ik,ir,ik,ir+nshift)
        ENDDO
      ENDDO
c
c...  Update target quantities:
c      DO in = 1, 



c
c     Update global grid quantities:
c
      nrs  = nrs  - nshift
      nrs2 = nrs2 - nshift

      IF (irref.LT.irsep)   irsep   = irsep   - nshift
      IF (irref.LE.irsep2)  irsep2  = irsep2  - nshift
      IF (irref.LT.irbreak) irbreak = irbreak - nshift
      IF (irref.LE.irtrap ) irtrap  = irtrap  - nshift
      IF (irref.LE.irtrap2) irtrap2 = irtrap2 - nshift
      IF (irref.LE.irwall)  irwall  = irwall  - nshift
      IF (irref.LE.irwall2) irwall2 = irwall2 - nshift
      DO ir = irref, nrs 
        IF (irorg2(ir).GT.0) irorg2(ir) = irorg2(ir) - nshift
      ENDDO

      irtrap2 = irtrap
      irwall2 = irwall

      RETURN
99    STOP
      END
c
c ======================================================================
c
c
      SUBROUTINE ResetRing(ir,irref)
      IMPLICIT none

      INTEGER ir,irref

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER FetchKORPG

      INTEGER ik,ii,id,ik1,ir1,id1
      LOGICAL valid

      IF (irref.NE.-1) nks(ir) = nks(irref)

      ikmids(ir) = nks(ir) / 2

c Need to recalculate quantities... i.e. kss ...
c Can't reset the boundary ring cells until cell center values can
c be recalculated... at the moment this seems to be being done in FindLink...
c move it here...?
      DO ik = 1, nks(ir)
c        CALL MoveCell(ik,ir,MAXNKS,MAXNRS)
        kbfs  (ik,ir) = 1.0
        bratio(ik,ir) = 0.0
        ktebs (ik,ir) = 1.0
        ktibs (ik,ir) = 1.0
        knbs  (ik,ir) = 1.0
        knes  (ik,ir) = 1.0
        kvhs  (ik,ir) = 0.0
        korpg (ik,ir) = FetchKORPG(ik,ir)
        id = korpg(ik,ir)
        DO ii = 1, nvertp(id)
          rvertp(ii,id) = 0.0
          zvertp(ii,id) = 0.0
        ENDDO
        nvertp(id) = 0
      ENDDO

      IF (irref.EQ.-1) THEN
c...    Assumes that connection map is built:
        DO ik = 1, nks(ir)
          thetag(ik,ir) = thetag(ikins(ik,ir),irins(ik,ir))
          virtag(ik,ir) = virtag(ikins(ik,ir),irins(ik,ir))
        ENDDO
        IF (ir.GT.irsep) THEN
          thetat(idds(ir,1)) = thetat(idds(irins(nks(ir),ir),1))
          thetat(idds(ir,2)) = thetat(idds(irins(1      ,ir),2))
        ENDIF
      ELSE
        DO ik = 1, nks(ir)
          thetag(ik,ir) = thetag(ik,irref)
          virtag(ik,ir) = virtag(ik,irref)
        ENDDO
        IF (ir.GT.irsep) THEN
          thetat(idds(ir,1)) = thetat(idds(irref,1))
          thetat(idds(ir,2)) = thetat(idds(irref,2))
        ENDIF
      ENDIF

      RETURN
 99   WRITE(0,*) 'IK,IR,NPOLYP:',ik,ir,npolyp,vpolmin
      STOP
      END
c
c ======================================================================
c
c subroutine: InsertCell
c
c At the moment, this routine does not adjust CUTPT1 or CUTPT2.
c
c The adjusting the target variables is not a concern since the cells
c added to the ends of a ring have zero volume.  The IDDS reference to
c calculated target quantities is a function of IR, not IK or NKS(IR)
c (which is modified here).
c
c have to get rid of the ikto/ikti restrictions - in DeleteCell also
c
c fix npolyp problem - adding real cells after adding virtual cells
c will screw things up if npolyp is not incremented with each
c virtual cell - be that is shitty for multiple calls to EIRENE
c where the korpg index could overflow - need to set aside a bank for
c for virtual cells that can be reset, and allow enough room for all
c the real cells that are going to be added
c
c
c
c
c The following cell quantities are accounted for:
c
c bratio,kfbs,kvhs,knbs,ktibs,ktebs,kss,ksb,rs,zs,korpg,rvertp,zvertp
c
c
      SUBROUTINE InsertCell(ikcell,ir,mode,type)
      USE mod_grid_divimp
      IMPLICIT none
c
c Input:
c
      INTEGER ikcell,ir,mode,type

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      REAL GetArea

      INTEGER ii,i1,iri,ik,ik1,in,id,ikstart,ikend,ikinf,midnks
c
c Cells cannot be deleted from core rings:
c
      IF (ir.LT.irsep.AND.ikcell.GE.nks(ir).AND.mode.EQ.AFTER)
     .  CALL ER('InsertCell','Cannot add end cell to core ring',*99)

c
c Check to see if there is room on the ring for a new cell:
c
      IF ((nks(ir)+1).GT.MAXNKS-5)
     .  CALL ER('InsertCell','Ring array space full',*99)

c
c Find index:
c
      IF (type.EQ.VIRTUAL) THEN
        IF (vpolyp+1.LE.MAXNKS*MAXNRS) THEN
          in = vpolyp + 1
        ELSE
          WRITE(0,*) 'ERROR (InsertCell): Virtual cell index is ',
     .               'out of bounds (VPOLYP = ',vpolyp,')'
          STOP
        ENDIF
      ELSEIF (type.EQ.NONVIRTUAL) THEN
        IF (grdnmod.EQ.0) THEN
          IF (npolyp+1.LT.vpolmin) THEN
            in = npolyp + 1
          ELSE
            WRITE(0,*) 'ERROR (InsertCell): Nonvirtual cell index is ',
     .                 'out of bounds (NPOLYP = ',npolyp,'  VPOLMIN',
     .                 ' = ',vpolmin,'  IKCELL = ',ikcell,'  IR = ',
     .                 ir,')'
            STOP
          ENDIF        
        ELSE
          IF (npolyp+1.GE.MAXNKS*MAXNRS) THEN
            CALL ER('InsertCell','NONVIRTUAL cell index out of bounds',
     .              *99)
          ELSE
            in = npolyp + 1
            IF (in.GE.vpolmin) THEN
              vpolmin = vpolmin + 1
              vpolyp = vpolmin
            ENDIF
          ENDIF        
        ENDIF
      ELSE
        WRITE(0,*) 'ERROR (InsertCell): Unsupported option (TYPE = ',
     .             type,')'
        STOP
      ENDIF
c
c Add a cell before cell IKCELL,IR:
c
      IF (mode.EQ.BEFORE) THEN
        ikstart = ikcell + 1
        ik      = ikcell
        ikinf   = ikcell + 1
      ELSEIF (mode.EQ.AFTER) THEN
        ikstart = ikcell + 2
        ik      = ikcell + 1
        ikinf   = ikcell
      ELSE
        WRITE(0,*) 'ERROR (InsertCell): Unsupported option (MODE = ',
     .             mode,')'
        STOP
      ENDIF

      IF (ikcell.EQ.nks(ir).AND.mode.EQ.AFTER) THEN
        ikend = 0
      ELSE
        ikend = nks(ir) + 1
      ENDIF
c
c Shift cell quantities along ring:
c
c      WRITE(0,*) 'STATUS: IK1 IR IKEND IKSTART = ',ik1,ir,ikend,ikstart

      DO ik1 = ikend, ikstart, -1
        CALL MoveCell(ik1,ir,ik1-1,ir)
      ENDDO
c
c Assign cell quantities:
c
      id = korpg(ikinf,ir)

      IF (type.EQ.VIRTUAL) THEN
        IF (mode.EQ.BEFORE) THEN
          rvertp(3,in) = rvertp(2,id)
          zvertp(3,in) = zvertp(2,id)
          rvertp(4,in) = rvertp(1,id)
          zvertp(4,in) = zvertp(1,id)
          rvertp(1,in) = rvertp(4,in)
          zvertp(1,in) = zvertp(4,in)
          rvertp(2,in) = rvertp(3,in)
          zvertp(2,in) = zvertp(3,in)
        ELSE
          rvertp(1,in) = rvertp(4,id)
          zvertp(1,in) = zvertp(4,id)
          rvertp(2,in) = rvertp(3,id)
          zvertp(2,in) = zvertp(3,id)
          rvertp(3,in) = rvertp(2,in)
          zvertp(3,in) = zvertp(2,in)
          rvertp(4,in) = rvertp(1,in)
          zvertp(4,in) = zvertp(1,in)
        ENDIF
      ELSE
        rvertp(:,in) = 0.0
        zvertp(:,in) = 0.0
        IF (ALLOCATED(d_rvertp)) THEN
          d_rvertp(:,in) = 0.0D0
          d_zvertp(:,in) = 0.0D0
        ENDIF
      ENDIF

      nvertp(in) = 4

      rs(ik,ir) = 0.25 * (rvertp(1,in) + rvertp(2,in) +
     +                    rvertp(3,in) + rvertp(4,in))
      zs(ik,ir) = 0.25 * (zvertp(1,in) + zvertp(2,in) +
     +                    zvertp(3,in) + zvertp(4,in))

c Should all of these quantities be assigned in this way...?
      IF (stopopt.EQ.121) THEN
        bratio(ik,ir) = -1.0
      ELSE
        bratio(ik,ir) = bratio(ikinf,ir)        
      ENDIF
      kbfs  (ik,ir) = kbfs  (ikinf,ir)
      kvhs  (ik,ir) = kvhs  (ikinf,ir)
      knbs  (ik,ir) = knbs  (ikinf,ir)
      knes  (ik,ir) = knes  (ikinf,ir)
      ktibs (ik,ir) = ktibs (ikinf,ir)
      ktebs (ik,ir) = ktebs (ikinf,ir)
      kes   (ik,ir) = kes   (ikinf,ir)
      korpg (ik,ir) = in
      kvols (ik,ir) = GetArea(ik,ir)

      pinion  (ik,ir) = pinion  (ikinf,ir)
      pinatom (ik,ir) = pinatom (ikinf,ir)
      pinmol  (ik,ir) = pinmol  (ikinf,ir)
      pinena  (ik,ir) = pinena  (ikinf,ir)
      pinenm  (ik,ir) = pinenm  (ikinf,ir)
      pinalpha(ik,ir) = pinalpha(ikinf,ir)
      DO i1 = 1, 6
        pinline (ik,ir,i1,H_BALPHA) = pinline (ikinf,ir,i1,H_BALPHA)
        pinline (ik,ir,i1,H_BGAMMA) = pinline (ikinf,ir,i1,H_BGAMMA)
      ENDDO
      pinz0   (ik,ir) = pinz0   (ikinf,ir)
      pinionz (ik,ir) = pinionz (ikinf,ir)
      pinenz  (ik,ir) = pinenz  (ikinf,ir)
      pinrec  (ik,ir) = pinrec  (ikinf,ir)
      pinqi   (ik,ir) = pinqi   (ikinf,ir)
      pinqe   (ik,ir) = pinqe   (ikinf,ir)
      pinmp   (ik,ir) = pinmp   (ikinf,ir)
c      pinvdist(ik,ir) = pinvdist(ikinf,ir)
      osmpei  (ik,ir) = osmpei  (ikinf,ir)
      osmpmk  (ik,ir) = osmpmk  (ikinf,ir)
      osmmp   (ik,ir) = osmmp   (ikinf,ir)
      osmqe   (ik,ir) = osmqe   (ikinf,ir)

      pinrec  (ik,ir) = pinrec  (ikinf,ir)

      DO in = 1, eirnsdtor
        DO i1 = 1, MAXBGK
          pinbgk(ik,ir,i1+(in-1)*MAXBGK) = 
     .      pinbgk(ikinf,ir,i1+(in-1)*MAXBGK)
        ENDDO
      ENDDO

      IF (ALLOCATED(divimp_ik)) THEN
        divimp_ik(ik,ir) = divimp_ik(ikinf,ir)
        divimp_ir(ik,ir) = divimp_ir(ikinf,ir)
      ENDIF

c
c *** LOOK HERE *** This may be trouble... but should work for the moment.
c
      kss(ik,ir) = 0.0
c      kss(ik,ir) = kss(ikinf,ir)

c Not sure if all this will still work...?

c *** FIX *** - have to catch ksb(0,ir) when setting the first cell -- should
c recalculate ksb and kps after adjusting the grid... not such a big deal
c if the grid is being modified just after begin read in... bigger deal later
c when the grid will be used and ksb has already been calculated...

      IF (ik.EQ.1.AND.ksb(0,ir).NE.0.0) THEN
        ksb(1,ir) = ksb(0,ir)
        ksb(0,ir) = 0.0
      ELSE
        IF (type.EQ.VIRTUAL) ksb(ik,ir) = 0.0
      ENDIF
c
c
c

      nks(ir) = nks(ir) + 1

c
c Check to see if IKTO and IKTI should be modified:
c
      IF (nopriv) THEN
        ikto2 = 0
        ikti2 = nks(irsep) + 1
      ELSE
        IF (ir.EQ.irsep.OR.ir.GT.irtrap.OR.
     .      (ir.GT.irsep.AND.ir.LT.irwall.AND.type.EQ.VIRTUAL)) THEN

          IF ((ikcell.LE.ikto2(ir).AND.mode.EQ.BEFORE).OR.
     .        (ikcell.LT.ikto2(ir).AND.mode.EQ.AFTER))
     .      ikto2(ir) = ikto2(ir) + 1

          IF ((ikcell.LE.ikti2(ir).AND.mode.EQ.BEFORE).OR.
     .        (ikcell.LT.ikti2(ir).AND.mode.EQ.AFTER))
     .      ikti2(ir) = ikti2(ir) + 1

        ELSEIF (ir.LT.irwall) THEN
          midnks = ikti - ikto - 1

          ikto2(ir) = MAX(1      ,(nks(ir) - midnks) / 2)
          ikti2(ir) = MIN(nks(ir),ikto2(ir) + midnks + 1)
        ENDIF
      ENDIF

      ikto = ikto2(irsep)
      ikti = ikti2(irsep)
c
c Distinquish between a real and virtual cell:
c
      IF (type.EQ.NONVIRTUAL) THEN
        virtag(ik,ir) = 0
        npolyp        = npolyp + 1
      ELSE
        virtag(ik,ir) = 1
        vpolyp        = vpolyp + 1
      ENDIF

      IF ((ikcell.LE.ikmids(ir).AND.mode.EQ.BEFORE).OR.
     .    (ikcell.LT.ikmids(ir).AND.mode.EQ.AFTER ))
     .  ikmids(ir) = ikmids(ir) + 1


      RETURN
99    CALL OutputData(85,'Error in InsertCell')
      WRITE(EROUT,'(5X,A,4I4)')
     .  ' IKCELL,IR MODE TYPE = ',ikcell,ir,mode,type
      STOP
      END
c
c ======================================================================
c
c FUNCION: GETAREA
c
      REAL FUNCTION GETAREA(IK,IR)

      INCLUDE 'params'
      INCLUDE 'cgeom'

      DOUBLE PRECISION a,b,base,c1,c2,area,r(4),z(4),rc,zc,height

      integer l,lp1,kp

      INTEGER IK,IR,I,IK1,id
c      REAL    A,B,BASE,C1,C2

      IF (NVERTP(KORPG(IK,IR)).NE.4) THEN
         WRITE(0,*) 'Warning (GETAREA): Cell does not'
     +              ,' have 4 verticies:',IK,IR
         GETAREA=0
         RETURN
      ENDIF

      area = 0.0

      kp = korpg(ik,ir)

      DO L = 1, NVERTP(KP)
        LP1 = L + 1
        IF (L.EQ.NVERTP(KP)) LP1 = 1
        area = area + (RVERTP(LP1,KP)*ZVERTP(L  ,KP) -
     +                 RVERTP(L  ,KP)*ZVERTP(LP1,KP))
      ENDDO

      GetArea = 0.5 * abs(area)

      RETURN
      END
c
c ======================================================================
c
c subroutine: DeleteCell
c
c KSS and KSB do not have to be recalculated if DeleteCell is being used
c to remove cells that were added using InsertCell, because all such cells
c have ...
c
c At the moment, this routine does not adjust CUTPT1 or CUTPT2.
c
c The following cell quantities are accounted for:
c
c bratio,kfbs,kvhs,knbs,ktibs,ktebs,kss,ksb,rs,zs,korpg,rvertp,zvertp
c
c ======================================================================

      SUBROUTINE DeleteCell(ikcell,ir)

      INTEGER ikcell,ir

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER ik,id1,id2
c
c Cells cannot be deleted from core rings:
c
      IF (ir.LT.irsep) THEN
        WRITE(0,*) 'ERROR (DeleteCell): Cells cannot be ',
     .             'deleted from core rings (IR = ',ir,')'
        STOP
      ENDIF


c      IF (ir.EQ.irsep.AND.
c     .    (ikcell.EQ.ikto2(ir).OR.ikcell.EQ.ikti2(ir)-1))
c     .  CALL ER('DeleteCell','Cannot delete cells near hard'//
c     .                       ' cut points',*99)
c      IF (ir.GT.irtrap.AND.
c     .    (ikcell.EQ.ikto2(ir).OR.ikcell.EQ.ikti2(ir)-1)) THEN
c        CALL WN('DeleteCell','Trying to delete cell near hard'//
c     .          'cut point in PFZ')
c        WRITE(0,*) '    IK,IR=',ikcell,ir
c      ENDIF

       IF ((ir.EQ.irsep.OR.ir.EQ.nrs).AND.
c       IF ((ir.EQ.irsep.OR.ir.GT.irtrap).AND.
     .    (ikcell.EQ.ikto2(ir).OR.ikcell.EQ.ikti2(ir)-1))
     .  CALL ER('DeleteCell','Cannot delete cells near hard'//
     .                       ' cut points',*99)


      korpg (ikcell,ir)        = MAXNKS*MAXNRS
      nvertp(korpg(ikcell,ir)) = 0


c
c     Delete cell IKCELL,IR:
c
      DO ik = ikcell, nks(ir)-1
        CALL MoveCell(ik,ir,ik+1,ir)
      ENDDO
c
c     Check to see if IKTO and IKTI should be modified:
c
      IF (ikcell.LT.ikto2(ir)) ikto2(ir) = ikto2(ir) - 1
      IF (ikcell.LT.ikti2(ir)) ikti2(ir) = ikti2(ir) - 1

      IF (ikcell.LE.ikmids(ir)) ikmids(ir) = ikmids(ir) - 1

      IF (ir.EQ.irsep.AND.ikcell.LT.ikto) ikto = ikto - 1
      IF (ir.EQ.irsep.AND.ikcell.LT.ikti) ikti = ikti - 1

      nks(ir) = nks(ir) - 1


      RETURN
99    WRITE(EROUT,*) 'IK IR = ',ik,ir
      STOP
      END
c
c ======================================================================

c Make sure that a virtual cell isn't being split...
c
      SUBROUTINE SplitCell(ikcell,ircell,splitpos,code)
      USE mod_grid_divimp
      IMPLICIT none
c
c Input:
c
      INTEGER ikcell,ircell
      REAL*8  splitpos

c     Output:
      INTEGER code

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER ik,ir,ir1,ii,id,ik1,ik2,id1,id2,nadd,hadd,ladd
      REAL*8  r(6),z(6),spos,len,deltar(4),deltaz(4)
      REAL    brat,newksb,ksb2,newkpb,kpb2,r1,z1

      code = 1

      ik = ikcell
      ir = ircell
      spos = splitpos
c      spos = DBLE(splitpos)
      id = korpg(ik,ir)
c
c
c
c      IF (nks(ir)+1.EQ.MAXNKS) THEN
      IF (nks(ir)+1.EQ.MAXNKS-5) THEN
        CALL MS('SplitCell','Ring array bound reached')
        RETURN
      ENDIF

      IF (ir.EQ.irsep.OR.(ir.GE.irtrap.AND.ir.LE.nrs)) THEN
        ladd = ikto
        hadd = nks(irsep) - ikti + 1

        DO ir1 = irtrap + 1, nrs
          ladd = MAX(ladd,ikto2(ir1))
          hadd = MAX(hadd,nks(ir1)-ikti2(ir1)+1)
        ENDDO

        WRITE(SLOUT,'(A,8I6)') 'SPLITCELL: STRUCT',ik,ir,ladd,hadd,
     .                ikto2(ir),nks(ir)-ikti2(ir)+1,
     .                ladd+hadd+(ikti-ikto-1)+1,nks(ir)

        IF (ik.LE.ikto2(ir).AND.ikto2(ir).GE.ladd.AND.
c     .      ladd+hadd+(ikti-ikto-1)+1.GT.MAXNKS) THEN
     .      ladd+hadd+(ikti-ikto-1)+1.GT.MAXNKS-10) THEN
          CALL WN('SplitCell','Structurability violation A')
          WRITE(EROUT,'(2I4)') ik,ir
          WRITE(0    ,'(2I4)') ik,ir
          RETURN
        ENDIF

        IF (ik.GE.ikti2(ir).AND.nks(ir)-ikti2(ir)+1.GE.hadd.AND.
c     .      ladd+hadd+(ikti-ikto-1)+1.GT.MAXNKS) THEN
     .      ladd+hadd+(ikti-ikto-1)+1.GT.MAXNKS-10) THEN
          CALL WN('SplitCell','Structurability violation B')
          WRITE(EROUT,'(3I4)') ik,ir,hadd
          WRITE(0    ,'(3I4)') ik,ir,hadd
          RETURN
        ENDIF



c        IF (ladd+hadd+(ikti-ikto-1)+1.GT.MAXNKS-12) THEN
c          CALL MS('SplitCell','Structurability violation C')
c          WRITE(EROUT,'(2I4)') ik,ir
c          RETURN
c        ENDIF

      ENDIF
c
c     Make sure there is room to construct outer wall ring:
c
c      IF (nbr.GT.0.AND.ir.EQ.nrs+1) THEN

      IF (grdnmod.EQ.0.AND.nbr.GT.0.AND.ir.EQ.irwall-1) THEN
        IF (stopopt.EQ.121) THEN

        ELSE
          r1 = rvertp(1,korpg(1,irbreak))
          z1 = zvertp(1,korpg(1,irbreak))

          ik2 = 0
          DO ik1 = 1, nks(irbreak-1)
            id1 = korpg(ik1,irbreak-1)

            IF (rvertp(3,id1).EQ.r1.AND.
     .          zvertp(3,id1).EQ.z1) ik2 = ik1
          ENDDO

c          WRITE(0,*) '   IK2 CELLS = ',ik2,ik2+nks(ir)+1

          IF (ik2.EQ.0) CALL ER('SplitCell','Cannot locate inner '//
     .                                      'wall contact A',*99)

          IF (ik2+nks(ir)+1.GE.MAXNKS-5) THEN
            CALL WN('SplitCell','Outer ring full')
            RETURN
          ENDIF
        ENDIF
      ENDIF


      IF (nbr.GT.0.AND.ir.EQ.irbreak-1) THEN
        IF (stopopt.EQ.121) THEN

        ELSEIF (grdnmod.NE.0) THEN
c...      Not sure why this was necessary, but I am turning it 
c         off for now.   

        ELSE
          r1 = rvertp(1,korpg(1,irbreak))
          z1 = zvertp(1,korpg(1,irbreak))
          ik2 = 0
          DO ik1 = 1, nks(irbreak-1)
            id1 = korpg(ik1,irbreak-1)
            IF (rvertp(3,id1).EQ.r1.AND.
     .          zvertp(3,id1).EQ.z1) ik2 = ik1
          ENDDO

          IF (ik2.EQ.0) CALL ER('SplitCell','Cannot locate inner '//
     .                                      'wall contact B',*99)

          IF (ik.LE.ik2.AND.ik2+nks(ir)+1.GE.MAXNKS-5) THEN
            CALL WN('SplitCell','Leaves no space for outer ring')
            RETURN
          ENDIF
        ENDIF
      ENDIF





      code = -1

      IF (splitpos.LE.0.0.OR.splitpos.GE.1.0)
     .  CALL ER('SplitCell','Invalid split position',*99)

      IF (ir.LT.irsep.AND.ik.GE.nks(ir))
     .  CALL ER('SplitCell','Cannot split last core cell',*99)

      IF (nvertp(id).NE.4)
     .  CALL ER('SplitCell','Cell must have 4 verticies',*99)

c
c     Determine vertices for new cells:
c
      IF (ALLOCATED(d_rvertp)) THEN
        r(1:4) = d_rvertp(1:4,id)
        z(1:4) = d_zvertp(1:4,id)
      ELSE
        r(1:4) = DBLE(rvertp(1:4,id))
        z(1:4) = DBLE(zvertp(1:4,id))
      ENDIF

      r(5) = r(1) + spos * (r(4) - r(1))
      z(5) = z(1) + spos * (z(4) - z(1))
      r(6) = r(2) + spos * (r(3) - r(2))
      z(6) = z(2) + spos * (z(3) - z(2))

c     Check that the length of the poloidal sides of the new cells
c     is sufficient:
c
      deltar(1) = r(1) - r(5)
      deltar(2) = r(2) - r(6)
      deltar(3) = r(5) - r(4)
      deltar(4) = r(6) - r(3)
      deltaz(1) = z(1) - z(5)
      deltaz(2) = z(2) - z(6)
      deltaz(3) = z(5) - z(4)
      deltaz(4) = z(6) - z(3)

      DO ii = 1, 4
c...bug
        len = DSQRT(deltar(ii)**2.0D0 + deltaz(ii)**2)
c        len = SQRT(deltar(ii)**2.0 + deltaz(ii)**2.0)

c        WRITE(SLOUT,*) 'MARK: SPLIT ',ikcell,ircell,len,deltar(ii),
c     .                 deltaz(ii),grd_minpl
c        WRITE(SLOUT,*) spos,r(1),r(4),r(5),z(1),z(4),z(5)

        IF (grdnmod.EQ.0.AND.len.LT.grd_minpl) THEN
          CALL MS('SplitCell','Side too short')
          RETURN
        ENDIF
      ENDDO

c Gross assumption...?
      brat   = bratio(ik,ir)

      newksb = ksb(ik-1,ir) + spos * (ksb(ik,ir) - ksb(ik-1,ir))
      ksb2   = ksb(ik,ir)

      newkpb = kpb(ik-1,ir) + spos * (kpb(ik,ir) - kpb(ik-1,ir))
      kpb2   = kpb(ik,ir)
c
c     Add a cell to the ring:
c
c       WRITE(0,*) 'NKS=',nks(ir),ir
      CALL InsertCell(ik,ir,AFTER,NONVIRTUAL)
c       WRITE(0,*) 'NKS=',nks(ir),ir

c This has to be added because of the way that InsertCell adjusts ikto...
      IF (ir.EQ.irsep.AND.ik.EQ.ikto.OR.ik.EQ.ikto2(ir)) THEN
        ikto2(ir) = ikto2(ir) + 1
        ikto      = ikto2(ir)
      ENDIF

      IF (ik.EQ.ikmids(ir)) ikmids(ir) = ikmids(ir) + 1

      ik1 = ik
      id1 = korpg(ik1,ir)

      ik2 = ik + 1
      id2 = korpg(ik2,ir)
c
c Assign cell quantities:
c
      rvertp(3,id2) = rvertp(3,id1)
      zvertp(3,id2) = zvertp(3,id1)
      rvertp(4,id2) = rvertp(4,id1)
      zvertp(4,id2) = zvertp(4,id1)
      IF (ALLOCATED(d_rvertp)) THEN
        d_rvertp(3,id2) = d_rvertp(3,id1)
        d_zvertp(3,id2) = d_zvertp(3,id1)
        d_rvertp(4,id2) = d_rvertp(4,id1)
        d_zvertp(4,id2) = d_zvertp(4,id1)
      ENDIF

c      rvertp(1,id1) = SNGL(r(1))
c      zvertp(1,id1) = SNGL(z(1))
c      rvertp(2,id1) = SNGL(r(2))
c      zvertp(2,id1) = SNGL(z(2))
      rvertp(3,id1) = SNGL(r(6))
      zvertp(3,id1) = SNGL(z(6))
      rvertp(4,id1) = SNGL(r(5))
      zvertp(4,id1) = SNGL(z(5))
      IF (ALLOCATED(d_rvertp)) THEN
        d_rvertp(3,id1) = r(6)
        d_zvertp(3,id1) = z(6)
        d_rvertp(4,id1) = r(5)
        d_zvertp(4,id1) = z(5)
      ENDIF

      rvertp(1,id2) = SNGL(r(5))
      zvertp(1,id2) = SNGL(z(5))
      rvertp(2,id2) = SNGL(r(6))
      zvertp(2,id2) = SNGL(z(6))
      IF (ALLOCATED(d_rvertp)) THEN
        d_rvertp(1,id2) = r(5)
        d_zvertp(1,id2) = z(5)
        d_rvertp(2,id2) = r(6)
        d_zvertp(2,id2) = z(6)
      ENDIF
c      rvertp(3,id2) = SNGL(r(3))
c      zvertp(3,id2) = SNGL(z(3))
c      rvertp(4,id2) = SNGL(r(4))
c      zvertp(4,id2) = SNGL(z(4))

      rs(ik1,ir) = 0.25 * (rvertp(1,id1) + rvertp(2,id1) +
     .                     rvertp(3,id1) + rvertp(4,id1))
      zs(ik1,ir) = 0.25 * (zvertp(1,id1) + zvertp(2,id1) +
     .                     zvertp(3,id1) + zvertp(4,id1))

      rs(ik2,ir) = 0.25 * (rvertp(1,id2) + rvertp(2,id2) +
     .                     rvertp(3,id2) + rvertp(4,id2))
      zs(ik2,ir) = 0.25 * (zvertp(1,id2) + zvertp(2,id2) +
     .                     zvertp(3,id2) + zvertp(4,id2))

c The magnetic field should be lineraly interpolated...
      bratio(ik1,ir) = brat
      bratio(ik2,ir) = brat

      IF (brat.NE.0.0) THEN
        kbfs(ik1,ir) = 1.0 / brat
        kbfs(ik2,ir) = 1.0 / brat
      ELSE
        kbfs(ik1,ir) = 1.0
        kbfs(ik2,ir) = 1.0
      ENDIF

      ksb(ik1,ir) = newksb
      ksb(ik2,ir) = ksb2
      kss(ik1,ir) = 0.5 * (ksb(ik1,ir) + ksb(ik1-1,ir))
      kss(ik2,ir) = 0.5 * (ksb(ik2,ir) + ksb(ik1  ,ir))

      kpb(ik1,ir) = newkpb
      kpb(ik2,ir) = kpb2
      kps(ik1,ir) = 0.5 * (kpb(ik1,ir) + kpb(ik1-1,ir))
      kps(ik2,ir) = 0.5 * (kpb(ik2,ir) + kpb(ik1  ,ir))
c
c     Adjust closing cell on core rings if necessary:
c
      IF (ir.LT.irsep.AND.ik.EQ.1) THEN
        rs    (nks(ir),ir) = rs    (1,ir)
        zs    (nks(ir),ir) = zs    (1,ir)
        korpg (nks(ir),ir) = korpg (1,ir)
        bratio(nks(ir),ir) = bratio(1,ir)
        kbfs  (nks(ir),ir) = kbfs  (1,ir)
      ENDIF

c      WRITE(SLOUT,'(A,3I4,2X,2I5,4F10.4)')
c     .  'SPLITCELL IK1,IK2 NEWKSB KSS2 = ',
c     .  ik1,ik2,nks(ir),korpg(ik1,ir),korpg(ik2,ir),
c     .  ksb(ik1,ir),ksb(ik2,ir),kss(ik1,ir),kss(ik2,ir)
c      WRITE(SLOUT,'(5X,6F9.6)') (r(ii),ii=1,6)
c      WRITE(SLOUT,'(5X,6F9.6)') (z(ii),ii=1,6)

      code = 0

      RETURN
99    WRITE(EROUT,'(5X,A,2I4)')
     .  'IK,IR = ',ik,ir
      END
c
c
c ======================================================================
c
c subroutine: MoveCell
c
c
      SUBROUTINE MoveCell(ik1,ir1,ik2,ir2)
      USE mod_grid_divimp
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER ik1,ik2,ir1,ir2,i1,i2


      IF (ik1.LT.1.OR.ik1.GT.nks(ir1)+1.OR.
     .    ik2.LT.1.OR.ik2.GT.nks(ir2)+1)
     .  CALL ER('MoveCell','Index out of bounds',*99)

      korpg (ik1,ir1) = korpg (ik2,ir2)
      rs    (ik1,ir1) = rs    (ik2,ir2)
      zs    (ik1,ir1) = zs    (ik2,ir2)
      bratio(ik1,ir1) = bratio(ik2,ir2)
      kbfs  (ik1,ir1) = kbfs  (ik2,ir2)
      kvhs  (ik1,ir1) = kvhs  (ik2,ir2)
      knbs  (ik1,ir1) = knbs  (ik2,ir2)
      knes  (ik1,ir1) = knes  (ik2,ir2)
      ktibs (ik1,ir1) = ktibs (ik2,ir2)
      ktebs (ik1,ir1) = ktebs (ik2,ir2)
      kes   (ik1,ir1) = kes   (ik2,ir2)
      kss   (ik1,ir1) = kss   (ik2,ir2)
      ksb   (ik1,ir1) = ksb   (ik2,ir2)
      kps   (ik1,ir1) = kps   (ik2,ir2)
      kpb   (ik1,ir1) = kpb   (ik2,ir2)
      kfords(ik1,ir1) = kfords(ik2,ir2)
      kbacds(ik1,ir1) = kbacds(ik2,ir2)
      kvols (ik1,ir1) = kvols (ik2,ir2)
      virtag(ik1,ir1) = virtag(ik2,ir2)
      thetag(ik1,ir1) = thetag(ik2,ir2)

      pinion  (ik1,ir1) = pinion  (ik2,ir2)
      pinatom (ik1,ir1) = pinatom (ik2,ir2)
      pinmol  (ik1,ir1) = pinmol  (ik2,ir2)
      pinena  (ik1,ir1) = pinena  (ik2,ir2)
      pinenm  (ik1,ir1) = pinenm  (ik2,ir2)
      pinalpha(ik1,ir1) = pinalpha(ik2,ir2)
      DO i1 = 1, 6
        pinline (ik1,ir1,i1,H_BALPHA) = pinline(ik2,ir2,i1,H_BALPHA)
        pinline (ik1,ir1,i1,H_BGAMMA) = pinline(ik2,ir2,i1,H_BGAMMA)
      ENDDO
      pinz0   (ik1,ir1) = pinz0   (ik2,ir2)
      pinionz (ik1,ir1) = pinionz (ik2,ir2)
      pinenz  (ik1,ir1) = pinenz  (ik2,ir2)
      pinrec  (ik1,ir1) = pinrec  (ik2,ir2)
      pinqi   (ik1,ir1) = pinqi   (ik2,ir2)
      pinqe   (ik1,ir1) = pinqe   (ik2,ir2)
      pinmp   (ik1,ir1) = pinmp   (ik2,ir2)
c...  TO BE REMOVED:
      DO i1 = 1, MAXDATA
        pindata(ik1,ir1,i1) = pindata(ik2,ir2,i1)
      ENDDO
      DO i1 = 1, MAXSTRATA
        DO i2 = 1, 3
          pinstrata(ik1,ir1,i2,i1) = pinstrata(ik2,ir2,i2,i1)
        ENDDO
      ENDDO
      DO i1 = 1, NMOMCHA
        pinploss(ik1,ir2,i1) = pinploss(ik2,ir2,i1)
      ENDDO
      DO i1 = 1, 3
        pinioncomp(ik1,ir2,i1) = pinioncomp(ik2,ir2,i1)
      ENDDO

c      pinvdist(ik1,ir) =  pinvdist(ik2,irref)
      osmpei  (ik1,ir1) = osmpei  (ik2,ir2)
      osmcfp  (ik1,ir1) = osmcfp  (ik2,ir2)
      osmcfe  (ik1,ir1) = osmcfe  (ik2,ir2)
      osmcfi  (ik1,ir1) = osmcfi  (ik2,ir2)
      osmmp   (ik1,ir1) = osmmp   (ik2,ir2)
      osmqe   (ik1,ir1) = osmqe   (ik2,ir2)

      mulrec  (ik1,ir1) = mulrec  (ik2,ir2)
      mulion  (ik1,ir1) = mulion  (ik2,ir2)
      mulqer  (ik1,ir1) = mulqer  (ik2,ir2)
      mulqei  (ik1,ir1) = mulqei  (ik2,ir2)
      DO i2 = 1, eirnsdtor
        DO i1 = 1, MAXBGK
          pinbgk(ik1,ir1,i1+(i2-1)*MAXBGK) =
     .      pinbgk(ik2,ir2,i1+(i2-1)*MAXBGK)
        ENDDO
      ENDDO

      pinior  (ik1,ir1) = pinior  (ik2,ir2)
      pinmpr  (ik1,ir1) = pinmpr  (ik2,ir2)
      pinqir  (ik1,ir1) = pinqir  (ik2,ir2)
      pinqer  (ik1,ir1) = pinqer  (ik2,ir2)

      IF (ALLOCATED(divimp_ik)) THEN
        divimp_ik(ik1,ir1) = divimp_ik(ik2,ir2)
        divimp_ir(ik1,ir1) = divimp_ir(ik2,ir2)
      ENDIF

      RETURN
99    WRITE(EROUT,'(5X,A,4I5)') 'IK1 IR1 IK2 IR2 = ',ik1,ir1,ik2,ir2
      WRITE(0    ,'(5X,A,4I5)') 'IK1 IR1 IK2 IR2 = ',ik1,ir1,ik2,ir2
      STOP
      END
c
c ======================================================================
c


c
c
c
c
c
c ======================================================================
c ======================================================================
c
c block: Building magnetic grid connection map
c
c
c ======================================================================
c
      INTEGER FUNCTION FetchKORPG(ik1,ir1)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      INTEGER ik,ir,ik1,ir1,id,lastid
      LOGICAL valid

      DATA lastid /1/
      SAVE

c...  Assign KORPG in sequence until the end of the allocated index space:
      IF (npolyp.LT.vpolmin-1) THEN
        FetchKORPG = npolyp + 1
        npolyp = npolyp + 1
c        WRITE(0,*) 'KORPG FOUND:',ik1,ir1,'FAST'
        RETURN
      ENDIF

c...  Search for a free KORPG index:
 10   DO id = lastid, vpolmin-1
        valid = .TRUE.
        ir = 0
        DO WHILE (valid.AND.ir.LT.MAXNRS)
          ir = ir + 1
          DO ik = 1, MAXNKS
            IF (ik.EQ.ik1.AND.ir.EQ.ir1) CYCLE
            IF (korpg(ik,ir).EQ.id) THEN
              valid = .FALSE.
              EXIT
            ENDIF
          ENDDO
        ENDDO
        IF (valid) THEN
          FetchKORPG = id
          lastid = id          
c          WRITE(0,*) 'KORPG FOUND:',ik1,ir1,korpg(ik1,ir1).EQ.id
          RETURN
        ENDIF
      ENDDO
      IF (lastid.NE.1) THEN
c      IF (ir.EQ.nrs.AND.lastid.NE.1) THEN
c...    One final search all the way through:
        lastid = 1
        GOTO 10
c        ir = 0
      ENDIF

      CALL ER('FetchKORPG','Unable to find index',*99)
      RETURN
 99   CALL OutputData(86,'Unable to find KORPG index')
      CALL DumpGrid('Unable to find KORPG index')
      STOP
      END
c
c
c ======================================================================
c
c subroutine: GenWallRing
c
      SUBROUTINE GenWallRing
      IMPLICIT   none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

c     Functions:
      INTEGER FetchKORPG
      REAL atan3c


c     Local variables:
      INTEGER IK,IR,STATUS

      INTEGER ii,ncell,wallik(2*MAXNKS),wallir(2*MAXNKS),tempik,tempir,
     .        i,in,ik1,ir1,id1,iavg,istart,i1,i2
      LOGICAL valid,debug
      REAL    rcen,zcen,cellr(2*MAXNKS),cellz(2*MAXNKS),ravg,zavg,
     .        wallth(2*MAXNKS),xpth,deltaz,deltar,tempr,tempz,tempth
c
c Initialize variables:
c
      debug = .TRUE.
      STATUS = 0
      i = 0
c
c Find center of main plasma... check if this is done elsewhere...
c make sure it is done for all grids...
c

c...C-MOD
c      rcen = rxp
c      zcen = zxp 
c...DIII-D
c      rcen = r0
c      zcen = z0

      rcen = 0.5 * (r0 + rxp)
      zcen = 0.5 * (z0 + zxp)
     
c    
c     Find cells without a neighbour in the 'outward' (increasing
c     IR) direction:
c    
      CALL OutputGrid(87,'Before generating wall ring')
     
      ncell = 0

      DO ir = irsep, nrs
        IF (idring(ir).EQ.-1) CYCLE
        DO ik = 1, nks(ir)
          IF (virtag(ik,ir).EQ.0.AND.
     .        (irouts(ik,ir).EQ.ir.OR.
     .         irins (ik,ir).EQ.ir)) THEN
c     .         irins (ik,ir).EQ.ir.AND.ir.LT.irwall)) THEN  ! Why the ir.LT.irwall..?
     .        
c          IF (irouts(ik,ir).EQ.ir.AND.virtag(ik,ir).EQ.0) THEN
     
            in = korpg(ik,ir)
     
c            write(0,*) '      ncell',ncell,2*MAXNKS

            IF (ncell+1.EQ.2*MAXNKS)
     .        CALL ER('GenWallRing','Array bound violation',*99)
     
            ncell = ncell + 1
c            write(0,*) 'ik,ir,ncell',ik,ir,ncell
             
            wallik(ncell) = ik
            wallir(ncell) = ir

            cellr(ncell) = 0.5 * (rvertp(2,in) + rvertp(3,in))
            cellz(ncell) = 0.5 * (zvertp(2,in) + zvertp(3,in))
c    
c           Calculate poloidal angle...
c    
            deltar = cellr(ncell) - rcen
            deltaz = cellz(ncell) - zcen

c Check this...
            wallth(ncell) = atan3c(deltaz,deltar)

            IF (grdnmod.NE.0) THEN
              wallth(ncell) = wallth(ncell) + 90.0
              IF (wallth(ncell).GT.0.0) 
     .          wallth(ncell) = wallth(ncell) - 360.0
            ENDIF

          ENDIF
        ENDDO
      ENDDO
      
      WRITE(0,*) 'GenWallRing: NCELL=',ncell
      IF (ncell+1.GT.2*MAXNKS)
     .  CALL ER('GenWallRing','Array bound violation, increase '//
     .          'MAXNKS',*99)

c...  Adjust WALLTH based on x-point location:
      IF (cgridopt.EQ.LINEAR_GRID.OR.cgridopt.EQ.RIBBON_GRID) THEN
        wallth(1:ncell) = cellz(1:ncell)
      ELSE
        deltar = rxp - rcen
        deltaz = zxp - zcen
        xpth = atan3c(deltaz,deltar)

        IF (debug) 
     .    WRITE(0,*) 'BROKEN MAP (NCELL,RCEN,ZCEN,XPTH): ',
     .      ncell,rcen,zcen,xpth
        
        write(88,*) 'what the fuck',rxp,rcen
        write(88,*) '             ',zxp,zcen
        write(88,*) '             ',xpth
        DO ii = 1, ncell
          write(88,*) 'wallth:',ii,wallth(ii),wallik(ii),wallir(ii)
        ENDDO
c
        DO in = 1, ncell
          IF (wallth(in).GT.xpth) wallth(in) = wallth(in) - 360.0
        ENDDO
        IF (grdnmod.NE.0.AND.irbreak.GT.irsep) THEN
          DO in = 1, ncell
            IF (wallik(in).EQ.1.AND.
     .          (wallir(in).EQ.grdtseg(grdntseg(1,IKLO),1,IKLO).OR.
     .           wallir(in).EQ.irwall-1)) THEN
              istart = in
              EXIT
            ENDIF
          ENDDO
          DO in = 1, ncell
            IF (wallth(in).GT.wallth(istart)) wallth(in)=wallth(in)-360.0
          ENDDO
        ENDIF
      ENDIF

      write(88,*) '             '
      DO ii = 1, ncell
        write(88,*) 'wallth:',ii,wallth(ii),wallik(ii),wallir(ii)
      ENDDO


c...  Sort wall segments:

10    status = 0
      DO ii = 1, ncell-1
        IF (wallth(ii).LT.wallth(ii+1)) THEN
          status = 1
          tempik = wallik(ii)
          tempir = wallir(ii)
          tempr  = cellr (ii)
          tempz  = cellz (ii)
          tempth = wallth(ii)
          wallik(ii) = wallik(ii+1)
          wallir(ii) = wallir(ii+1)
          cellr (ii) = cellr (ii+1)
          cellz (ii) = cellz (ii+1)
          wallth(ii) = wallth(ii+1)
          wallik(ii+1) = tempik
          wallir(ii+1) = tempir
          cellr (ii+1) = tempr
          cellz (ii+1) = tempz
          wallth(ii+1) = tempth
        ENDIF
      ENDDO
      IF (status.EQ.1) GOTO 10

c...  Make sure segments are sorted properly for each
c     sub-section of the wall ring that are associated
c     with a particular non-boundary ring:

      IF (.TRUE.) THEN
c...    Make sure that cells for particular IK regions are 
c       clustered together (the WALLTH ordering gets things
c       close, but not always quite right):
        DO i1 = 1, ncell-1          
          DO i2 = i1+1, ncell
            IF ((wallir(i2)  .EQ.wallir(i1)  .AND.
     .           wallik(i2)  .EQ.wallik(i1)+1).OR.
     .          (wallir(i2)  .EQ.wallir(i1)  .AND.
     .           wallir(i1+1).NE.wallir(i1)  .AND.
     .           .TRUE.)) THEN
              tempik = wallik(i1+1)
              tempir = wallir(i1+1)
              tempr  = cellr (i1+1)
              tempz  = cellz (i1+1)
              tempth = wallth(i1+1)
              wallik(i1+1) = wallik(i2)
              wallir(i1+1) = wallir(i2)
              cellr (i1+1) = cellr (i2)
              cellz (i1+1) = cellz (i2)
              wallth(i1+1) = wallth(i2)
              wallik(i2) = tempik
              wallir(i2) = tempir
              cellr (i2) = tempr
              cellz (i2) = tempz
              wallth(i2) = tempth
              EXIT
            ENDIF
          ENDDO
        ENDDO
c...    Sort each region in ascending IK order:
        DO i1 = 1, ncell-1          
          DO i2 = i1+1, ncell
            IF (wallir(i2).EQ.wallir(i1)  .AND.
     .          wallik(i2).LT.wallik(i1)) THEN
     
              tempik = wallik(i1)
              tempir = wallir(i1)
              tempr  = cellr (i1)
              tempz  = cellz (i1)
              tempth = wallth(i1)
     
              wallik(i1) = wallik(i2)
              wallir(i1) = wallir(i2)
              cellr (i1) = cellr (i2)
              cellz (i1) = cellz (i2)
              wallth(i1) = wallth(i2)
     
              wallik(i2) = tempik
              wallir(i2) = tempir
              cellr (i2) = tempr
              cellz (i2) = tempz
              wallth(i2) = tempth
            ENDIF
          ENDDO
        ENDDO
      ELSE
c        status = 0
c        DO WHILE (status.EQ.0) 
c          status = 1
c          DO ii = 1, ncell-1
c            IF (wallir(ii).EQ.wallir(ii+1).AND.
c     .          wallik(ii).GT.wallik(ii+1)) THEN
c              status = 0
c     
c              tempik = wallik(ii)
c              tempir = wallir(ii)
c              tempr  = cellr (ii)
c              tempz  = cellz (ii)
c              tempth = wallth(ii)
c     
c              wallik(ii) = wallik(ii+1)
c              wallir(ii) = wallir(ii+1)
c              cellr (ii) = cellr (ii+1)
c              cellz (ii) = cellz (ii+1)
c              wallth(ii) = wallth(ii+1)
c     
c              wallik(ii+1) = tempik
c              wallir(ii+1) = tempir
c              cellr (ii+1) = tempr
c              cellz (ii+1) = tempz
c              wallth(ii+1) = tempth
c            ENDIF
c          ENDDO
c        ENDDO
      ENDIF     

      IF (debug) THEN
        DO ii = 1, ncell
          WRITE(50,*) 'IRCELLS:',wallik(ii),wallir(ii),wallth(ii)
        ENDDO
      ENDIF


c      CALL DumpGrid('BUMMER MAN')


c...  Assign cell quantities for IR=IRWALL:
      ir = irwall
      nks(ir) = ncell
      DO ik = 1, ncell
        rs    (ik,ir) = cellr (ik)
        zs    (ik,ir) = cellz (ik)
        ikins (ik,ir) = wallik(ik)
        irins (ik,ir) = wallir(ik)
        ikouts(ik,ir) = ik
        irouts(ik,ir) = irwall
        IF (irins(wallik(ik),wallir(ik)).EQ.wallir(ik)) THEN
c...      For the "secondary PFZ" of double-null extended grids:
          ikins(wallik(ik),wallir(ik)) = ik
          irins(wallik(ik),wallir(ik)) = irwall
        ELSE     
          ikouts(wallik(ik),wallir(ik)) = ik
          irouts(wallik(ik),wallir(ik)) = irwall
        ENDIF    
        thetag(ik,ir) = thetag(wallik(ik),wallir(ik))
        virtag(ik,ir) = virtag(wallik(ik),wallir(ik))
      ENDDO
      thetat(idds(ir,2)) = thetag(nks(ir),ir)
      thetat(idds(ir,1)) = thetag(nks(ir),ir)
      DO ik = 1, nks(ir)
        kbfs  (ik,ir) = 1.0
        bratio(ik,ir) = 0.0
        korpg (ik,ir) = FetchKORPG(ik,ir)
        ktebs (ik,ir) = 1.0
        ktibs (ik,ir) = 1.0
        knbs  (ik,ir) = 1.0
        knes  (ik,ir) = 1.0
        kvhs  (ik,ir) = 0.0
        thetag(ik,ir) = thetag(ikins(ik,ir),irins(ik,ir))
        virtag(ik,ir) = virtag(ikins(ik,ir),irins(ik,ir))
      ENDDO
      thetat(idds(ir,1)) = thetat(idds(irins(nks(ir),ir),1))
      thetat(idds(ir,2)) = thetat(idds(irins(1      ,ir),2))
c...  Scan wall ring and make sure that the cells are ordered properly.  This 
c     should not be a problem, but you never know:
      IF (grdnmod.NE.0) THEN
        ir = irwall
        DO i1 = 1, nks(ir)-1
          DO i2 = i1+1, nks(ir)
            IF (irins(i1,ir).EQ.irins(i2,ir).AND.
     .          ikins(i1,ir).GE.ikins(i2,ir)) 
     .        CALL ER('GenWallRing','IRWALL cells out of order',*99)
          ENDDO
        ENDDO
      ENDIF

      RETURN
99    CALL OutputData(85,'Error in GenWallRing')
      WRITE(0,*) 'IK,NKS=',i1,i2,nks(ir)
      WRITE(0,*) 'IKINS,IRINS 1:',ikins(i1,ir),irins(i1,ir)
      WRITE(0,*) 'IKINS,IRINS 2:',ikins(i2,ir),irins(i2,ir)
      STOP
      END
c
c ======================================================================
c
c subroutine: FindLink
c
c Check if sending a reference to a virtual cell...?
c
c Currently assumes that core rings are structured...
c
      SUBROUTINE FindLink(ikcell,ircell,side,iklink,irlink)
      USE mod_grid_divimp
      IMPLICIT none

c     Input:
      INTEGER ikcell,ircell,side

c     Output:
      INTEGER iklink,irlink

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      REAL*8     TOL
      PARAMETER (TOL = 1.0E-07)

      REAL     CenLen
      INTEGER  CalcPoint

      INTEGER ik,ir,i1,in,incell,ikmin,irmin,ikend,cpc,cpd,id,
     .        nmatch,ikmatch(MAXNKS),irmatch(MAXNRS)
      LOGICAL clean,side_reset,cont
      REAL    dist,distmin
      REAL*8  a1,a2,b1,b2,c1,c2,d1,d2,t1,t2

      side_reset = .FALSE.

c
c     Initialization:
c
      CALL IZero(ikmatch,MAXNKS)

      incell = korpg(ikcell,ircell)

c...  Handle the special case of boundary rings, and set up vertex
c     information for standard cases:
      IF (virtag(ikcell,ircell).GT.0) THEN
        irlink = ircell
        iklink = ikcell
        RETURN
      ENDIF

      IF (side.EQ.SIDE14) THEN
        IF (ircell.EQ.1.OR.ircell.EQ.irtrap) THEN
          irlink = ircell
          iklink = ikcell
          in = korpg(ikcell,ircell+1)
          rs(ikcell,ircell) = 0.5 * (rvertp(2,in) + rvertp(3,in))
          zs(ikcell,ircell) = 0.5 * (zvertp(2,in) + zvertp(3,in))
          RETURN
        ELSEIF (ircell.EQ.irtrap+1.OR.
     .          ircell.EQ.2.OR.
     .          ircell.EQ.irwall) THEN
          irlink = ircell - 1
          iklink = ikcell
          RETURN
        ELSE
          IF (ALLOCATED(d_rvertp)) THEN
            c1 = d_rvertp(1,incell)
            c2 = d_zvertp(1,incell)
            d1 = d_rvertp(4,incell)
            d2 = d_zvertp(4,incell)
          ELSE 
            c1 = DBLE(rvertp(1,incell))
            c2 = DBLE(zvertp(1,incell))
            d1 = DBLE(rvertp(4,incell))
            d2 = DBLE(zvertp(4,incell))
          ENDIF
          IF (ircell.EQ.irsep.AND.
     .        (ikcell.LE.ikto.OR.ikcell.GE.ikti)) THEN
            ir = nrs
          ELSE
            ir = ircell - 1
          ENDIF
        ENDIF
      ELSEIF (side.EQ.SIDE23) THEN
        IF (ircell.EQ.irwall) THEN
          irlink = ircell
          iklink = ikcell
          in = korpg(ikcell,ircell-1)
          rs(ikcell,ircell) = 0.5 * (rvertp(1,in) + rvertp(4,in))
          zs(ikcell,ircell) = 0.5 * (zvertp(1,in) + zvertp(4,in))
          RETURN
c Assume that the core rings are sturctured...
        ELSEIF (ircell.EQ.1.OR.ircell.EQ.irtrap.OR.
     .          (ircell.EQ.irwall-1.AND.irbreak.EQ.0)) THEN
          irlink = ircell + 1
          iklink = ikcell
          RETURN
        ELSE
          IF (ALLOCATED(d_rvertp)) THEN
            c1 = d_rvertp(2,incell)
            c2 = d_zvertp(2,incell)
            d1 = d_rvertp(3,incell)
            d2 = d_zvertp(3,incell)
          ELSE 
            c1 = DBLE(rvertp(2,incell))
            c2 = DBLE(zvertp(2,incell))
            d1 = DBLE(rvertp(3,incell))
            d2 = DBLE(zvertp(3,incell))
          ENDIF
          IF (ircell.EQ.nrs) THEN
            ir = irsep
          ELSE
            ir = ircell + 1
          ENDIF
        ENDIF
      ELSE
        CALL ER('FindLink','Illegal side option',*99)
      ENDIF

c      IF (ircell.EQ.388.AND.side.EQ.SIDE14) THEN
c        WRITE(0,*) 'BUILDMAP 388:',ir
c      ENDIF

c...  Loop over rings to find neighbouring cells:
      nmatch = 0
      cont = .TRUE.
      DO WHILE (cont)
        cont  = .FALSE.
        clean = .FALSE.
        IF     (ir.EQ.ircell.OR.idring(ir).EQ.BOUNDARY) THEN
          ikend = 0
        ELSEIF (ir.LT.irsep) THEN
          ikend = nks(ir) - 1
        ELSE
          ikend = nks(ir)
        ENDIF
        DO ik = 1, ikend
c         Exclude virtual cells, unless specified cell is also virtual:
          IF (virtag(ik,ir).NE.0) CYCLE
          in = korpg(ik,ir)
          IF (side.EQ.SIDE23) THEN
            IF (ALLOCATED(d_rvertp)) THEN
              a1 = d_rvertp(1,in)
              a2 = d_zvertp(1,in)
              b1 = d_rvertp(4,in)
              b2 = d_zvertp(4,in)
            ELSE 
              a1 = DBLE(rvertp(1,in))
              a2 = DBLE(zvertp(1,in))
              b1 = DBLE(rvertp(4,in))
              b2 = DBLE(zvertp(4,in))
            ENDIF
          ELSE
            IF (ALLOCATED(d_rvertp)) THEN
              a1 = d_rvertp(2,in)
              a2 = d_zvertp(2,in)
              b1 = d_rvertp(3,in)
              b2 = d_zvertp(3,in)
            ELSE 
              a1 = DBLE(rvertp(2,in))
              a2 = DBLE(zvertp(2,in))
              b1 = DBLE(rvertp(3,in))
              b2 = DBLE(zvertp(3,in))
            ENDIF
          ENDIF
          cpc = CalcPoint(a1,a2,b1,b2,c1,c2,t1)
          cpd = CalcPoint(a1,a2,b1,b2,d1,d2,t2)
c          IF (side.EQ.SIDE14.AND.
c     .        ircell.EQ.388.AND.ir.EQ.2) THEN
c            WRITE(88,'(A,4I6,2X,2I6,2F15.8)') 
c     .        'LINK 388:',ikcell,ircell,ik,ir,cpc,cpd,t1,t2
c            WRITE(88,'(A,4F15.8)')
c     .        '        :',a1,a2,b1,b2
c            WRITE(88,'(A,4F15.8)')
c     .        '        :',c1,c2,d1,d2
c          ENDIF
          IF ((cpc.EQ.1.AND.t1.LT.1.0D0).AND.
     .        (cpd.EQ.1.AND.t2.GT.0.0D0)) THEN
c...        Both end points are connected to a neighbouring cell.  This is 
c           the "classic" assumption when building the connection map, so 
c           exit the search (this could be improved in the future):
            nmatch = 1
            iklink = ik
            irlink = ir
            clean = .TRUE.
            cont  = .FALSE.
            EXIT  ! Leave IK loop
          ELSEIF ((cpc.EQ.1.AND.t1.LT.1.0D0).OR.
     .            (cpd.EQ.1.AND.t2.GT.0.0D0)) THEN
c...        Only one end point is connected with a neighbouring cell, so
c           store the cell index in a list of candidates:
            nmatch          = nmatch + 1
            ikmatch(nmatch) = ik
            irmatch(nmatch) = ir
            cont = .TRUE.
          ENDIF
        ENDDO

        IF (nmatch.EQ.0.OR..NOT.clean) THEN
c...      No connection was identified, so check the entire grid:
          cont = .TRUE.
          IF (side.EQ.SIDE14) THEN
            IF (side_reset) THEN
              ir = ir - 1
            ELSE
              side_reset = .TRUE.
              ir = nrs
            ENDIF
          ELSE
            IF (side_reset) THEN
              ir = ir + 1
           ELSE
              side_reset = .TRUE.
              ir = 3
            ENDIF
          ENDIF

          IF (side.EQ.SIDE23.AND.
     .        ikcell.EQ.26.AND.ircell.EQ.120) THEN
             IF (sloutput) THEN
               WRITE(88,*) 'LINK PROBLEM: SEARCHING',ir
             ENDIF 
          ENDIF

c...      Exit when a virtual ring is encountered:
          IF ((side.EQ.SIDE14.AND.ir.EQ.1    ).OR.
     .        (side.EQ.SIDE23.AND.ir.EQ.nrs+1)) cont = .FALSE.


c          IF (ircell.EQ.388.AND.side.EQ.SIDE14) THEN
c            WRITE(0,*) 'BUILDMAP 388 searching:',ir
c          ENDIF


        ENDIF
      ENDDO

      IF (.NOT.clean.AND.nmatch.LE.1) THEN
c...    If no match was found, so point the cell to itself (for
c       later processing):
        iklink = ikcell
        irlink = ircell
      ELSEIF (clean) THEN
c...    Neighbouring cell was found:
      ELSE
c...    More than one match was found, so select the one that gives the shortest
c       distance between cell centres:
        distmin = HI
        DO i1 = 1, nmatch
          ik = ikmatch(i1)
          ir = irmatch(i1)
          dist = CenLen(ikcell,ircell,ik,ir)
          IF (dist.LT.distmin) THEN
            distmin = dist
            ikmin   = ik
            irmin   = ir
          ENDIF
        ENDDO
        iklink = ikmin
        irlink = irmin
      ENDIF

      RETURN
99    STOP
      END
c
c ======================================================================
c
c SUBROUTINE: BuildMap
c
c There is a problem here due to the assignment of ikouts to the wall
c ring before the check is made for a triangular cell.  As a result,
c a cell is added to the constructed wall ring that shouldn't be there,
c because the parent cell is triangular.
c
c jul 9, 97 - need to account for virtual cells with NO POLYGONS
c
      SUBROUTINE BuildMap
      USE mod_grid_divimp
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'
c
c Declare functions:
c
      REAL    CalcWidth,atan3c

c
c Declare local variables:
c
      INTEGER IK,IR,iki,iko,iri,iro,irs,ire,idb,id,id1,id2,ik1,ir1
c
c     Update boundary ring data in case cells have been
c     inserted/deleted from neighbouring rings:
c
c Make sure the cell centers for bounary rings are recalculated at the
c correct time...
c      WRITE(0,*) 'Resetting boundary rings'

      IF (ALLOCATED(d_rvertp)) THEN
        rvertp = SNGL(d_rvertp)
        zvertp = SNGL(d_zvertp)
      ENDIF


c      ! jdemod
c      if (cgridopt.ne.RIBBON_GRID) then 
       CALL ResetRing(1     ,2)
       CALL ResetRing(irwall,irwall-1)
       CALL ResetRing(irtrap,irtrap+1)
c      endif
c
c
c      WRITE(0,*) 'MESSAGE (BuildMap): Generating connection map'

      WRITE(50,*)
      WRITE(50,*) 'Building connection map:'
      WRITE(50,*)

c      WRITE(0,*) 'BUILDMAP: FINDLINK ASSUMES CORE IS OKAY - FIX'

      CALL DB('Calling FindLink')


      CALL OutputData(85,'Before calling FindLink')
c
c Build the connection map:
c
      DO ir = 1, nrs
        DO ik = 1, nks(ir)
          CALL FindLink(ik,ir,SIDE14,ikins (ik,ir),irins (ik,ir))
          CALL FindLink(ik,ir,SIDE23,ikouts(ik,ir),irouts(ik,ir))
        ENDDO
      ENDDO
c
c Handle wall ring for broken grids:
c
c      WRITE(0,*) 'IRBREAK:',irbreak
c      CALL DumpGrid('GEN WALL RING')


      CALL DB('Calling GenWallRing')
      IF (nbr.GT.0) CALL GenWallRing



c
c Calculate KINDS and KOUTDS, the distances between cell centers:
c
      CALL MS('BuildMap','*** IMPROVE KINDS/FINDS ***')

      DO ir = 1, nrs
        DO ik = 1, nks(ir)
          IF (ir.EQ.1.OR.ir.EQ.irtrap) THEN
            kinds(ik,ir) = 0.0
          ELSE
            iki = ikins(ik,ir)
            iri = irins(ik,ir)
            kinds(ik,ir) = SQRT((rs(ik,ir) - rs(iki,iri))**2.0 +
     .                          (zs(ik,ir) - zs(iki,iri))**2.0)
          ENDIF
          IF (ir.EQ.irwall) THEN
            koutds(ik,ir) = 0.0
          ELSE
            iko = ikouts(ik,ir)
            iro = irouts(ik,ir)
            koutds(ik,ir) = SQRT((rs(ik,ir) - rs(iko,iro))**2.0 +
     .                           (zs(ik,ir) - zs(iko,iro))**2.0)
          ENDIF
        ENDDO
      ENDDO
c
c     IKBREAK:
c

      CALL IZero(ikbreak,MAXNRS)

      IF (nbr.GT.0) THEN

        IF (irbreak.LE.irwall) THEN
          irs = irbreak - 1
          ire = irsep
        ELSE
          irs = irbreak - 1
          ire = irtrap  + 1
        ENDIF

        IF (irbreak.EQ.irsep) STOP 'FOLLOWING CODE WILL FAIL'

        DO ir = irs, ire, -1
          IF (stopopt.EQ.121.OR.stopopt.EQ.123) THEN   ! 123 for last gasp of thesis code...
c...        DIII-D plenum field line extensions:
            DO ik = nks(ir), 1, -1
              IF (ir.EQ.irbreak-1) THEN
                idb = korpg(nks(irbreak),irbreak)
                id  = korpg(ik          ,ir     )
                IF (zvertp(3,id).EQ.zvertp(4,idb).AND.
     .              rvertp(3,id).EQ.rvertp(4,idb)) ikbreak(ir) = ik
              ELSE
                IF (ikouts(ik,ir).GE.ikbreak(ir+1)) ikbreak(ir) = ik
              ENDIF
            ENDDO
          ELSE
            DO ik = 1, nks(ir)
              IF (ir.EQ.irbreak-1) THEN
                IF (irouts(ik,ir).EQ.irwall       ) ikbreak(ir) = ik + 1
              ELSE
                IF (ikouts(ik,ir).LE.ikbreak(ir+1)) ikbreak(ir) = ik
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF

c...  Clear connection map values for virtual cells:
      DO ir = 1, nrs
        DO ik = 1, nks(ir)
          IF (virtag(ik,ir).EQ.1) THEN
            irins (ik,ir) = -1
            irouts(ik,ir) = -1
          ENDIF
        ENDDO
      ENDDO

      CALL DB('Done in BuildMap - short return')

c.... Assign polygons for virtual rings:

c...  Do me a little scan to be sure that there is no KORPG overlap:
      DO ir = 1, nrs
        DO ik = 1, nks(ir)
          DO ir1 = 1, nrs
            DO ik1 = 1, nks(ir1)
              IF (korpg(ik,ir).EQ.korpg(ik1,ir1).AND.
     .            .NOT.(ir.EQ.ir1.AND.ik.EQ.ik1).AND.
     .            .NOT.(ir.LT.irsep.AND.(ik.EQ.1.OR.
     .                                   ik.EQ.nks(ir)))) THEN
                CALL ER('BuildMap','KORPG overlap detected',*99)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      
c      ! jdemod - do not need to modify polygons on the ribbon grid
c      !        - however, it may be a good idea to insert one virtual ring
c      !          for ir=1 in case too much code assumes IR=1 is virtual in core
c      !
c      if (cgridopt.ne.RIBBON_GRID) then 

      ir = 1
      DO ik = 1, nks(ir)
        IF (virtag(ik,ir).EQ.1) CYCLE
        id1 = korpg(ik,ir)
        id2 = korpg(ikouts(ik,ir),irouts(ik,ir))
        nvertp(id1) = 4
        rvertp(1,id1) = rvertp(1,id2)
        zvertp(1,id1) = zvertp(1,id2)
        rvertp(2,id1) = rvertp(1,id2)
        zvertp(2,id1) = zvertp(1,id2)
        rvertp(3,id1) = rvertp(4,id2)
        zvertp(3,id1) = zvertp(4,id2)
        rvertp(4,id1) = rvertp(4,id2)
        zvertp(4,id1) = zvertp(4,id2)
        rs(ik,ir) = 0.5 * (rvertp(1,id1) + rvertp(3,id1))
        zs(ik,ir) = 0.5 * (zvertp(1,id1) + zvertp(3,id1))
c        WRITE(0,*) 'NVERTP',ik,ir,id1,nvertp(id1)
      ENDDO

      ir = irwall
      DO ik = 1, nks(ir)
        IF (virtag(ik,ir).EQ.1) CYCLE
        id1 = korpg(ik,ir)
        id2 = korpg(ikins(ik,ir),irins(ik,ir))
        nvertp(id1) = 4
        IF (irins(ikins(ik,ir),irins(ik,ir)).EQ.irwall) THEN
c...      For the "secondary PFZ" of double-null extended grids:
          rvertp(1,id1) = rvertp(1,id2)
          zvertp(1,id1) = zvertp(1,id2)
          rvertp(2,id1) = rvertp(1,id2)
          zvertp(2,id1) = zvertp(1,id2)
          rvertp(3,id1) = rvertp(4,id2)
          zvertp(3,id1) = zvertp(4,id2)
          rvertp(4,id1) = rvertp(4,id2)
          zvertp(4,id1) = zvertp(4,id2)
        ELSE
          rvertp(1,id1) = rvertp(2,id2)
          zvertp(1,id1) = zvertp(2,id2)
          rvertp(2,id1) = rvertp(2,id2)
          zvertp(2,id1) = zvertp(2,id2)
          rvertp(3,id1) = rvertp(3,id2)
          zvertp(3,id1) = zvertp(3,id2)
          rvertp(4,id1) = rvertp(3,id2)
          zvertp(4,id1) = zvertp(3,id2)
        ENDIF
        rs(ik,ir) = 0.5 * (rvertp(1,id1) + rvertp(3,id1))
        zs(ik,ir) = 0.5 * (zvertp(1,id1) + zvertp(3,id1))
      ENDDO

      ir = irtrap
      DO ik = 1, nks(ir)
        IF (virtag(ik,ir).EQ.1) CYCLE
        id1 = korpg(ik,ir)
        id2 = korpg(ikouts(ik,ir),irouts(ik,ir))
        nvertp(id1) = 4
        rvertp(1,id1) = rvertp(1,id2)
        zvertp(1,id1) = zvertp(1,id2)
        rvertp(2,id1) = rvertp(1,id2)
        zvertp(2,id1) = zvertp(1,id2)
        rvertp(3,id1) = rvertp(4,id2)
        zvertp(3,id1) = zvertp(4,id2)
        rvertp(4,id1) = rvertp(4,id2)
        zvertp(4,id1) = zvertp(4,id2)
        rs(ik,ir) = 0.5 * (rvertp(1,id1) + rvertp(3,id1))
        zs(ik,ir) = 0.5 * (zvertp(1,id1) + zvertp(3,id1))
      ENDDO

c      ! jdemod - end of ribbon grid IF
c      endif

      CALL DB('Done in BuildMap')
c      CALL OutputData(87,'END OF BUILDMAP')
c      STOP 'sdfsd'


      RETURN
 99   WRITE(0,*) 'DATA:',ik,ir,ik1,ir1
      CALL DumpGrid('CHECKING KORPG OVERLAP IN BUILDMAP')
      STOP
      END
c
c
c
c
c
c
c
c ========================================================================
c
c ========================================================================
c
c block: Grid Shaping
c
c SEGCHK
c SHAPETARGET
c TAILORGRID
c
c ======================================================================
c
c function: SegChk
c
c 
      LOGICAL FUNCTION SegChk(ik,ir,rvp,zvp,side,nseg,rseg,zseg,
     .                        r,z)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

c...  Input:
      INTEGER ik,ir,side,nseg
      REAL*8  rvp(5,MAXNKS*MAXNRS),zvp(5,MAXNKS*MAXNRS),
     .        rseg(2*nseg),zseg(2*nseg)

c...  Output:
      REAL*8  r,z

      REAL*8     DTOL
      PARAMETER (DTOL = 1.0D-07)

      INTEGER in,id
      REAL*8  a1,a2,b1,b2,c1,c2,d1,d2,f,t1,t2,mint1

      SegChk = .FALSE.

c...  Validate input:
      id = korpg(ik,ir)
      IF     (side.EQ.SIDE14) THEN
        a1 = rvp(1,id)
        a2 = zvp(1,id)
        b1 = rvp(4,id)
        b2 = zvp(4,id)
      ELSEIF (side.EQ.SIDE23) THEN
        a1 = rvp(2,id)
        a2 = zvp(2,id)
        b1 = rvp(3,id)
        b2 = zvp(3,id)
      ELSEIF (side.EQ.SIDE34) THEN
        a1 = rvp(3,id)
        a2 = zvp(3,id)
        b1 = rvp(4,id)
        b2 = zvp(4,id)
      ELSE
        WRITE(0,*) 'ERROR (WallChk): Illegal SIDE option ',
     .             'specified (IK = ',ik,'  IR = ',ir,
     .             '  SIDE = ',side,')'
        STOP
      ENDIF

c...  Search for intesections between the grid cells and the line segment
c     specified in the OEDGE input file, returning the first intersection:
      mint1 = 1.0D+30
      DO in = 1, 2*nseg, 2
        c1 = rseg(in)
        c2 = zseg(in)
        d1 = rseg(in+1)
        d2 = zseg(in+1)
c...    Find the intersection points:
        CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,t1,t2)
c...    This is necessary because the DIII-D SONNET grids are already refined
c       at the target, while the C-Mod grids are not:
c        f = 10000.0D0
c...    C-Mod:
        f = 1.0D0
        IF ((t1+DTOL*f.GE.0.0D0.AND.t1-DTOL*f.LE.1.0D0).AND.
     .      (t2+DTOL*f.GE.0.0D0.AND.t2-DTOL*f.LE.1.0D0).AND.
     .      t1.LT.mint1) THEN
          SegChk = .TRUE.
          mint1 = t1
          r = a1 + t1 * (b1 - a1)
          z = a2 + t1 * (b2 - a2)
          WRITE(50,'(3X,A,2I4,1X,I3,2X,I3,2X,2F11.7,L4)')
     .      'SegChk IK,IR,T,S WS T1,T2 WC = ',
     .      ik,ir,side,in,t1,t2,SegChk
          RETURN
        ENDIF
      ENDDO

      RETURN
      END
c
c ======================================================================
c
      SUBROUTINE MorphGrid(mode,index)
      USE mod_grid_divimp
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom' 

      INTEGER mode,index,i1,ik,ir,id,ike,ind1
      REAL*8  deltax,deltay,rv,zv,fx,fy,x1,x2,y1,y2


      IF     (mode.GE.2.AND.mode.LE.7) THEN

c...    Lower divertor stretch:

        ind1 = index + 1

        DO WHILE (grdmod(ind1,1).EQ.0.0.AND.ind1.LT.grdnmod) 

          IF (grdmod(ind1+1,1).NE.0.0) 
     .      CALL ER('MorphGrid','Invalid GRDMOD data',*99)
      
          x1 = DBLE(grdmod(ind1,2))   
          x2 = DBLE(grdmod(ind1,3))   
          y1 = DBLE(grdmod(ind1,4))   
          y2 = DBLE(grdmod(ind1,5))
          deltax = DBLE(grdmod(ind1+1,2))   
          deltay = DBLE(grdmod(ind1+1,3))

          WRITE(0,*) 'MODE:',ind1,mode
          WRITE(0,*) 'MODE:',x1,x2
          WRITE(0,*) 'MODE:',y1,y2
          WRITE(0,*) 'MODE:',deltax,deltay

          DO ir = 1, nrs
            ike = nks(ir)
            IF (ir.LT.irsep) ike = ike - 1
            DO ik = 1, ike
              id = korpg(ik,ir)
              DO i1 = 1, nvertp(id)
                IF (ALLOCATED(d_rvertp)) THEN
                  rv = d_rvertp(i1,id)
                  zv = d_zvertp(i1,id)
                ELSE
                  rv = rvertp(i1,id)
                  zv = zvertp(i1,id)
                ENDIF
                IF (rv.LT.MIN(x1,x2).OR.rv.GT.MAX(x1,x2).OR.
     .              zv.LT.MIN(y1,y2).OR.zv.GT.MAX(y1,y2)) CYCLE

                SELECTCASE (mode)
                  CASE (2)
                    fx = 1.0
                    fy = 1.0 - DABS((zv - 0.5*(y1+y2)) / (0.5*(y2-y1)))
                    rvertp(i1,id) = rv 
                    zvertp(i1,id) = zv + deltay * fx * fy
                    IF (ALLOCATED(d_rvertp)) THEN
                      d_rvertp(i1,id) = rv 
                      d_zvertp(i1,id) = zv + deltay * fx * fy
                    ELSE
                      rvertp(i1,id) = rv 
                      zvertp(i1,id) = zv + deltay * fx * fy
                    ENDIF
                  CASE (3)
                    fx = DABS((zv - y1) / (y2 - y1))**1.0
                    fy = DABS((rv - x1) / (x2 - x1))**0.5
                    IF (ALLOCATED(d_rvertp)) THEN
                      d_rvertp(i1,id) = rv + deltax * fx * fy
                      d_zvertp(i1,id) = zv
                    ELSE
                      rvertp(i1,id) = rv + deltax * fx * fy
                      zvertp(i1,id) = zv
                    ENDIF
                  CASE (4)
c...                General adjustment which goes to zero at the domain boundary:
                    fx = 1.0 - DABS((rv - 0.5*(x1+x2)) / (0.5*(x2-x1)))
                    fy = 1.0 - DABS((zv - 0.5*(y1+y2)) / (0.5*(y2-y1)))
                    IF (ALLOCATED(d_rvertp)) THEN
                      d_rvertp(i1,id) = rv + deltax * fx * fy
                      d_zvertp(i1,id) = zv + deltay * fx * fy
                    ELSE
                      rvertp(i1,id) = rv + deltax * fx * fy
                      zvertp(i1,id) = zv + deltay * fx * fy
                    ENDIF
                  CASE (5)
c...                Horizontal/vertical adjustment that goes to zero at Y1,Y2/X1,X2:
                    IF (deltax.NE.0.0) THEN
                      fx = 1.0 - DABS((rv-0.5*(x1+x2)) / (0.5*(x2-x1)))
                      IF (ALLOCATED(d_rvertp)) THEN
                        d_rvertp(i1,id) = rv + deltax * fx
                      ELSE
                        rvertp(i1,id) = rv + deltax * fx
                      ENDIF
                    ENDIF
                    IF (deltay.NE.0.0) THEN
                      fy = 1.0 - DABS((zv-0.5*(y1+y2)) / (0.5*(y2-y1)))
                      IF (ALLOCATED(d_rvertp)) THEN
                        d_zvertp(i1,id) = zv + deltay * fy
                      ELSE
                        zvertp(i1,id) = zv + deltay * fy
                      ENDIF
                    ENDIF
                  CASE (6)
c...                Horizonta/vertical adjustment that goes to zero at Y1,1:
                    IF (deltax.NE.0.0) THEN
                      fx = DABS((zv - y1) / (y2 - y1))**0.5
                      fy = DABS((rv - x1) / (x2 - x1))**1.0
                      IF (ALLOCATED(d_rvertp)) THEN
                        d_rvertp(i1,id) = rv + deltax * fx * fy
                      ELSE
                        rvertp(i1,id) = rv + deltax * fx * fy
                      ENDIF
                    ENDIF
                    IF (deltay.NE.0.0) THEN
                      fx = DABS((zv - y1) / (y2 - y1))**1.0
                      fy = DABS((rv - x1) / (x2 - x1))**0.5
                      IF (ALLOCATED(d_rvertp)) THEN
                        d_zvertp(i1,id) = zv + deltay * fx * fy
                      ELSE
                        zvertp(i1,id) = zv + deltay * fx * fy
                      ENDIF  
                    ENDIF
                  CASE (7)
c...                Yank:
                    IF (deltax.NE.0.0) THEN
                      fx = 1.0 - DABS((zv-0.5*(y1+y2)) / (0.5*(y2-y1)))
                      IF (ALLOCATED(d_rvertp)) THEN
                        d_rvertp(i1,id) = rv + deltax * fx 
                      ELSE
                        rvertp(i1,id) = rv + deltax * fx 
                      ENDIF
                    ENDIF
                    IF (deltay.NE.0.0) THEN
                      fy = 1.0 - DABS((rv-0.5*(x1+x2)) / (0.5*(x2-x1)))
                      IF (ALLOCATED(d_rvertp)) THEN
                        d_zvertp(i1,id) = zv + deltay * fy
                      ELSE
                        zvertp(i1,id) = zv + deltay * fy
                      ENDIF
                    ENDIF

                  CASE DEFAULT
                    CALL ER('MorphGrid','Unknown MODE',*99)
                ENDSELECT

              ENDDO
            ENDDO
          ENDDO      

          ind1 = ind1 + 2
    
        ENDDO

      ELSEIF (mode.EQ.1) THEN
c...    hidelta10:

c...    Lower core/SOL stretch:
        x1 =  0.1
        x2 =  1.5
        y1 = -0.1
        y2 = -1.9
c        y1 = -1.9
c        y2 = -0.1

        deltax =  0.0
        deltay = -0.2

        DO ir = 1, nrs
          ike = nks(ir)
          IF (ir.LT.irsep) ike = ike - 1
          DO ik = 1, ike
            id = korpg(ik,ir)
            DO i1 = 1, nvertp(id)
              rv = rvertp(i1,id)
              zv = zvertp(i1,id)
              IF (rv.GE.x1.AND.rv.LE.x2.AND.
     .            zv.GE.MIN(y1,y2).AND.zv.LE.MAX(y1,y2)) THEN
                fx = 1.0
                fy = 1.0 - DABS((zv-0.5*(y1+y2)) / (0.5*(y2-y1)))
                IF (ALLOCATED(d_rvertp)) THEN
                  d_rvertp(i1,id) = rv 
                  d_zvertp(i1,id) = zv + deltay * fx * fy
                ELSE
                  rvertp(i1,id) = rv 
                  zvertp(i1,id) = zv + deltay * fx * fy
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO      

c...    Upper core/SOL stretch:
        x1 =  0.1
        x2 =  1.5
        y1 =  0.1
        y2 =  1.9
        deltax =  0.0
        deltay = +0.2
        DO ir = 1, nrs
          ike = nks(ir)
          IF (ir.LT.irsep) ike = ike - 1
          DO ik = 1, ike
            id = korpg(ik,ir)
            DO i1 = 1, nvertp(id)
              rv = rvertp(i1,id)
              zv = zvertp(i1,id)
              zv = zvertp(i1,id)
              IF (rv.GE.x1.AND.rv.LE.x2.AND.
     .            zv.GE.MIN(y1,y2).AND.zv.LE.MAX(y1,y2)) THEN
c              IF (rv.GE.x1.AND.rv.LE.x2.AND.
c     .            zv.GE.y1.AND.zv.LE.y2) THEN
                fx = 1.0
                fy = 1.0 - DABS((zv - 0.5*(y1+y2)) / (0.5*(y2-y1)))
                IF (ALLOCATED(d_rvertp)) THEN
                  d_rvertp(i1,id) = rv 
                  d_zvertp(i1,id) = zv + deltay * fx * fy
                ELSE
                  rvertp(i1,id) = rv 
                  zvertp(i1,id) = zv + deltay * fx * fy
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO      

c...    Lower divertor stretch:
        x1 =  0.5
        x2 =  1.7
        y1 = -1.1
        y2 = -2.0
c        y1 = -2.0
c        y2 = -1.1
        deltax = -0.4
        deltay =  0.0
        DO ir = 1, nrs
          ike = nks(ir)
          IF (ir.LT.irsep) ike = ike - 1
          DO ik = 1, ike
            id = korpg(ik,ir)
            DO i1 = 1, nvertp(id)
              rv = rvertp(i1,id)
              zv = zvertp(i1,id)
              IF (rv.GE.x1.AND.rv.LE.x2.AND.
     .            zv.GE.MIN(y1,y2).AND.zv.LE.MAX(y1,y2)) THEN
                fx = DABS((zv - y1) / (y2 - y1))**1.0
c                fx = ABS((zv - y2) / (y2 - y1))**1.0
                fy = DABS((rv - x1) / (x2 - x1))**0.5
                IF (ALLOCATED(d_rvertp)) THEN
                  d_rvertp(i1,id) = rv + deltax * fx * fy
                  d_zvertp(i1,id) = zv
                ELSE
                  rvertp(i1,id) = rv + deltax * fx * fy
                  zvertp(i1,id) = zv
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO      

c...    Upper divertor stretch:
        x1 =  0.5
        x2 =  1.7
        y1 =  1.1
        y2 =  2.0
        deltax = -0.4
        deltay =  0.0
        DO ir = 1, nrs
          ike = nks(ir)
          IF (ir.LT.irsep) ike = ike - 1
          DO ik = 1, ike
            id = korpg(ik,ir)
            DO i1 = 1, nvertp(id)
              rv = rvertp(i1,id)
              zv = zvertp(i1,id)
              IF (rv.GE.x1.AND.rv.LE.x2.AND.
     .            zv.GE.MIN(y1,y2).AND.zv.LE.MAX(y1,y2)) THEN
c     .            zv.GE.y1.AND.zv.LE.y2) THEN
                fx = DABS((zv - y1) / (y2 - y1))**1.0
                fy = DABS((rv - x1) / (x2 - x1))**0.5
                IF (ALLOCATED(d_rvertp)) THEN
                  d_rvertp(i1,id) = rv + deltax * fx * fy
                  d_zvertp(i1,id) = zv
                ELSE
                  rvertp(i1,id) = rv + deltax * fx * fy
                  zvertp(i1,id) = zv
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO      


      ELSE
        CALL ER('MorphGrid','Unrecognized mode',*99)
      ENDIF

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: ShapeTarget
c
c
c
      SUBROUTINE ShapeTarget(dataindex,rvp,zvp)
      USE mod_grid_divimp
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

c     Input:
      INTEGER dataindex
      REAL*8  rvp(5,MAXNKS*MAXNRS),zvp(5,MAXNKS*MAXNRS)


      COMMON /GRID/ iktop,irout,irin
      INTEGER       iktop(MAXNRS),irout(MAXNRS),irin(MAXNRS)

      INTEGER WallChk
      LOGICAL SegChk

      INTEGER ik,ik1,ik2,ir,ii,id,id1,id2,status1,status2,irs,ire,i1,i2,
     .        i3,nseg,type,mode
      REAL*8  len(MAXNKS),totlen,t,rseg(2*MAXNRS),zseg(2*MAXNRS)
      REAL*8  r1,z1,r2,z2,dist1,dist2
      REAL*8  minlength

c...  Maintain compatability with thesis version:
      IF (grdnmod.EQ.0) THEN
        STOP 'ShapeTarget_Old NO LONGER SUPPORTED'
        CALL ShapeTarget_Old(rvp,zvp)
        RETURN
      ENDIF

c...    This is necessary because the DIII-D grids are already refined near the
c       targets.  However, this scale limit is still arbitrary and may not
c       jive with future DIII-D grids:
c      minlength = 1.0D-04
c...    C-Mod:
      minlength = 2.0D-03

c...  Check GRDMOD array, which is loaded in the OEDGE input file and contains
c     the parameters for modifying the grid, for rings that need to be 
c     modified:
      i1 = dataindex

c...  Assing cut parameters:
      type = NINT(grdmod(i1,1))
      mode = NINT(grdmod(i1,2))
      irs  = NINT(grdmod(i1,4))
      ire  = NINT(grdmod(i1,5))

      IF (irs.EQ.0.AND.ire.EQ.0) RETURN

      WRITE(SLOUT,'(A,4I6)') 
     .  'SHAPETARGET - MODE,IRS,IRE:',mode,irs,ire,i1
      IF (irs.LT.irsep.OR.ire.GT.nrs.OR.irs.GT.ire)
     .  CALL ER('ShapeTarget','Invalid ring indexing',*99)
c...  Load line segment data:
      nseg = -1
      rseg = 0.0D0
      zseg = 0.0D0
      i2 = 1
c...  Advance the index pointer (I1+I2) to the start of the 
c     line segment data:
      DO WHILE (grdmod(i1+i2,1).NE.0.0.AND.i1+i2.LE.grdnmod) 
        i2 = i2 + 1
      ENDDO
      IF (i1+i2.GT.grdnmod) 
     .  CALL ER('ShapeTargets','Line segment data not found',*99)
c...  Load line segment data:
      DO WHILE (grdmod(i1+i2,1).EQ.0.0.AND.i1+i2.LE.grdnmod) 
        nseg = nseg + 2
        rseg(nseg  ) = DBLE(grdmod(i1+i2,2))
        zseg(nseg  ) = DBLE(grdmod(i1+i2,3))
        rseg(nseg+1) = DBLE(grdmod(i1+i2,4))
        zseg(nseg+1) = DBLE(grdmod(i1+i2,5))
        i2 = i2 + 1
      ENDDO
      IF (nseg.EQ.-1) 
     .  CALL ER('ShapeTarget','Line segment data not found',*99)

c...  Need to make room for a new ring if one of the original
c     targets is not discarded:
      IF (type.EQ.1.AND.mode.EQ.3) THEN
c...    For at least one case (likely when keeping the outer target), a 
c       ring is going to be added:
c
c...    Cases: 1 - Cut inner target 
c              2 - Cut outer target
c              3 - Cut mid-ring, keeping original inner and outer targets
c                  (which requires inserting a new ring for the outer half)
c
c...    Duplicate the ring and insert before IRWALL:
c
c...    Move grid vertex data into RVERTP and ZVERTP:
c        DO i2 = 1, MAXNKS*MAXNRS
c          DO i3 = 1, 4
c            rvertp(i3,i2) = rvp(i3,i2)
c            zvertp(i3,i2) = zvp(i3,i2)
c          ENDDO
c        ENDDO
        d_rvertp = rvp
        d_zvertp = zvp
        DO ir = irs, ire
          CALL DupeRing(ir)
        ENDDO
        rvp = d_rvertp
        zvp = d_zvertp
c        DO i2 = 1, MAXNKS*MAXNRS
c          DO i3 = 1, 4
c            rvp(i3,i2) = rvertp(i3,i2)
c            zvp(i3,i2) = zvertp(i3,i2)
c          ENDDO
c        ENDDO
      ENDIF

      irbreak = MIN(irbreak,irs)
     
      IF (type.EQ.1) THEN
        DO ir = irs, ire
          IF (ir.EQ.irwall.OR.ir.EQ.irtrap) CYCLE

c...  Need to set origin ring as well!
          id = idring(ir)
c...      Correct the IDRING specification for rings with a type 3 cut:
          IF (id.LT.-50) id = id + 100

          idring(ir) = -1
c
c TARTOTAR  - 1
c TARTOWAL  - 2
c WALTOTAR  - 3
c WALTOWAL  - 4
c
          IF (mode.EQ.1.AND.id.EQ.TARTOTAR) idring(ir) = WALTOTAR
          IF (mode.EQ.1.AND.id.EQ.TARTOWAL) idring(ir) = WALTOWAL
          IF (mode.EQ.1.AND.id.EQ.WALTOTAR) idring(ir) = WALTOTAR
          IF (mode.EQ.1.AND.id.EQ.WALTOWAL) idring(ir) = WALTOWAL

          IF (mode.EQ.2.AND.id.EQ.TARTOTAR) idring(ir) = TARTOWAL
          IF (mode.EQ.2.AND.id.EQ.WALTOTAR) idring(ir) = WALTOWAL
          IF (mode.EQ.2.AND.id.EQ.WALTOWAL) idring(ir) = WALTOWAL

          IF (mode.EQ.3.AND.id.EQ.TARTOTAR) idring(ir) = TARTOWAL
          IF (mode.EQ.3.AND.id.EQ.TARTOWAL) idring(ir) = TARTOWAL
          IF (mode.EQ.3.AND.id.EQ.WALTOTAR) idring(ir) = WALTOWAL
c...      IDRING unchanged:
          IF (mode.EQ.3.AND.id.EQ.WALTOWAL) idring(ir) = WALTOWAL

          IF (idring(ir).EQ.-1) 
     .      CALL ER('ShapeTarget','Unable to set IDRING',*99)                
        ENDDO
      ENDIF  
c
c...  Loop over rings:
c
      DO ir = irs, ire
        IF (ir.EQ.irwall.OR.ir.EQ.irtrap) CYCLE

        IF (mode.EQ.1) THEN
c
c...      Low index target:
c
          DO ik1 = nks(ir), 1, -1
            IF (SegChk(ik1,ir,rvp,zvp,SIDE14,nseg,rseg,zseg,r1,z1))
     .        GOTO 120
          ENDDO
          CALL ER('ShapeTarget','No intersection-low side 1-4',*99)
120       CONTINUE
          DO ik2 = nks(ir), 1, -1
            IF (SegChk(ik2,ir,rvp,zvp,SIDE23,nseg,rseg,zseg,r2,z2))
     .        GOTO 130
          ENDDO
          CALL ER('ShapeTarget','No intersection-low side 2-3',*99)
130       CONTINUE

          totlen = 0.0D0
          t      = 0.0D0

          id1 = korpg(ik1,ir)
          id2 = korpg(ik2,ir)
c
c...      Check to see if the length of the cropped cell side will be too
c         short.  If so, then move the vertex of the cell side being
c         cut over to the intersection (between the cell side and the
c         wall) point:
          dist1 = DSQRT((rvp(4,id1) - r1)**2 + (zvp(4,id1) - z1)**2)
          dist2 = DSQRT((rvp(3,id2) - r2)**2 + (zvp(3,id2) - z2)**2)

          WRITE(SLOUT,'(A,3I4,2F10.6)')
     .      'Checking side lengths LOW INDEX:',ir,ik1,ik2,dist1,dist2

          IF (dist1.LT.minlength) THEN
c
c Do the global search and update for all verticies so that
c MAX(1.0,REAL(ik2-ik1+1)) can be used to scale min side length...
c
            rvp(4,id1) = r1
            zvp(4,id1) = z1

            ik1 = ik1 + 1
            id1 = korpg(ik1,ir)
            rvp(1,id1) = r1
            zvp(1,id1) = z1
 
            WRITE(SLOUT,'(A,2I4,F10.6)') 'Move vertex 4 :',ik1-1,ir,
     .                                   dist1
          ENDIF

          IF (dist2.LT.minlength) THEN
            rvp(3,id2) = r2
            zvp(3,id2) = z2
            ik2 = ik2 + 1
            id2 = korpg(ik2,ir)
            rvp(2,id2) = r2
            zvp(2,id2) = z2
            WRITE(SLOUT,'(A,2I4,F10.6)') 'Move vertex 3 :',ik2-1,ir,
     .                                   dist2
          ENDIF
c
          IF (ik1.EQ.ik2) THEN
            rvp(1,id1) = r1
            zvp(1,id1) = z1
            rvp(2,id2) = r2
            zvp(2,id2) = z2
          ELSEIF (ik1.GT.ik2) THEN
            rvp(1,id2) = r1
            zvp(1,id2) = z1
            rvp(2,id2) = r2
            zvp(2,id2) = z2
            DO ik = ik2, ik1
              id      = korpg(ik,ir)
              len(ik) = DSQRT((rvp(2,id) - rvp(3,id))**2 +
     .                        (zvp(2,id) - zvp(3,id))**2)
              totlen  = totlen + len(ik)
            ENDDO
            DO ik = ik1, ik2+1, -1
              id = korpg(ik,ir)
              t  = t + len(ik) / totlen
              rvp(1,id) = rvp(4,id1) + t * (r1 - rvp(4,id1))
              zvp(1,id) = zvp(4,id1) + t * (z1 - zvp(4,id1))
              rvp(4,korpg(ik-1,ir)) = rvp(1,id)
              zvp(4,korpg(ik-1,ir)) = zvp(1,id)
              WRITE(50,*) 'splitting side: ',ik,ir,'low ik target'
            ENDDO
          ELSEIF (ik2.GT.ik1) THEN
            rvp(1,id1) = r1
            zvp(1,id1) = z1
            rvp(2,id1) = r2
            zvp(2,id1) = z2
            DO ik = ik1, ik2
              id      = korpg(ik,ir)
              len(ik) = DSQRT((rvp(1,id) - rvp(4,id))**2 +
     .                        (zvp(1,id) - zvp(4,id))**2)
              totlen  = totlen + len(ik)
            ENDDO
            DO ik = ik2, ik1+1, -1
              id = korpg(ik,ir)
              t  = t + len(ik) / totlen
              rvp(2,id) = rvp(3,id2) + t * (r2 - rvp(3,id2))
              zvp(2,id) = zvp(3,id2) + t * (z2 - zvp(3,id2))
              rvp(3,korpg(ik-1,ir)) = rvp(2,id)
              zvp(3,korpg(ik-1,ir)) = zvp(2,id)
              WRITE(50,*) 'splitting side: ',ik,ir,t,
     .                    'low index target'
            ENDDO
          ENDIF
c         Delete excess cells:
          IF (ik1.GT.1.AND.ik2.GT.1) THEN
            DO ii = 1, MIN(ik1,ik2)-1
              iktop(ir) = iktop(ir) - 1
              CALL DeleteCell(1,ir)
              WRITE(50,'(3X,2A,I3,A,I3,A,I4,A,I4)') 'ShapeTarget:',
     .        ' DELETING CELL  ring ',ir,' :  target ',1,
     .        ' : korpg ',id,' :  nks ',nks(ir)
            ENDDO
          ENDIF
        ELSEIF (mode.EQ.2.OR.mode.EQ.3) THEN
c
c...      High index target:
c
          WRITE(SLOUT,*) 'SHAPE TARGET: IR= ',ir
          WRITE(SLOUT,*) '1-4'
          DO ik1 = 1, nks(ir)
            IF (SegChk(ik1,ir,rvp,zvp,SIDE14,nseg,rseg,zseg,r1,z1))
     .        GOTO 140
          ENDDO
          CALL ER('ShapeTarget','No intersection- high side 1-4',*99)
140       CONTINUE
          WRITE(SLOUT,*) '2-3'
          DO ik2 = 1, nks(ir)
            IF (SegChk(ik2,ir,rvp,zvp,SIDE23,nseg,rseg,zseg,r2,z2))
     .        GOTO 150
          ENDDO
          WRITE(SLOUT,*) 'NKS =', nks(ir)
          CALL ER('ShapeTarget','No intersection- high side 2-3',*99)
150       CONTINUE
          totlen = 0.0D0
          t      = 0.0D0
          id1 = korpg(ik1,ir)
          id2 = korpg(ik2,ir)
c
          dist1 = DSQRT((rvp(1,id1)-r1)**2 + (zvp(1,id1)-z1)**2)
          dist2 = DSQRT((rvp(2,id2)-r2)**2 + (zvp(2,id2)-z2)**2)
            WRITE(SLOUT,'(A,3I4,2F10.6)')
     .        'Checking side lengths HIGH INDEX:',ir,ik1,ik2,dist1,
     .        dist2
          IF (dist1.LT.minlength) THEN
            rvp(1,id1) = r1
            zvp(1,id1) = z1
            ik1 = ik1 - 1
            id1 = korpg(ik1,ir)
            rvp(4,id1) = r1
            zvp(4,id1) = z1
            WRITE(SLOUT,'(A,2I4,F10.6)') 'Move vertex 1 :',ik1+1,ir,
     .                                   dist1
          ENDIF
 
          IF (dist2.LT.minlength) THEN
            rvp(2,id2) = r2
            zvp(2,id2) = z2
            ik2 = ik2 - 1
            id2 = korpg(ik2,ir)
            rvp(3,id2) = r2
            zvp(3,id2) = z2
            WRITE(SLOUT,'(A,2I4,F10.6)') 'Move vertex 2 :',ik2+1,ir,
     .                                   dist2
          ENDIF
c
          IF (ik1.EQ.ik2) THEN
            rvp(4,id1) = r1
            zvp(4,id1) = z1
            rvp(3,id2) = r2
            zvp(3,id2) = z2
            WRITE(SLOUT,*) 'IK1=IK2:',ik1,ik2
          ELSEIF (ik1.GT.ik2) THEN
            WRITE(SLOUT,*) 'IK1>IK2:',ik1,ik2
            rvp(4,id1) = r1
            zvp(4,id1) = z1
            rvp(3,id1) = r2
            zvp(3,id1) = z2
            DO ik = ik2, ik1
              id      = korpg(ik,ir)
              len(ik) = DSQRT((rvp(1,id) - rvp(4,id))**2 +
     .                        (zvp(1,id) - zvp(4,id))**2)
              totlen  = totlen + len(ik)
            ENDDO
            DO ik = ik2, ik1-1
              id = korpg(ik,ir)
              t  = t + len(ik) / totlen
              rvp(3,id) = rvp(2,id2) + t * (r2 - rvp(2,id2))
              zvp(3,id) = zvp(2,id2) + t * (z2 - zvp(2,id2))
              rvp(2,korpg(ik+1,ir)) = rvp(3,id)
              zvp(2,korpg(ik+1,ir)) = zvp(3,id)
              WRITE(50,'(5X,2A,2I4,F10.4,A)') '   :',
     .          ' Splitting side: ',ik,ir,t,' outer target'
            ENDDO
          ELSEIF (ik2.GT.ik1) THEN
            WRITE(SLOUT,*) 'IK1<IK2:',ik1,ik2
            rvp(4,id2) = r1
            zvp(4,id2) = z1
            rvp(3,id2) = r2
            zvp(3,id2) = z2
            DO ik = ik1, ik2
              id      = korpg(ik,ir)
              len(ik) = DSQRT((rvp(2,id) - rvp(3,id))**2 +
     .                        (zvp(2,id) - zvp(3,id))**2)
              totlen  = totlen + len(ik)
            ENDDO
            DO ik = ik1, ik2-1
              id = korpg(ik,ir)
              t  = t + len(ik) / totlen
              rvp(4,id) = rvp(1,id1) + t * (r1 - rvp(1,id1))
              zvp(4,id) = zvp(1,id1) + t * (z1 - zvp(1,id1))
              rvp(1,korpg(ik+1,ir)) = rvp(4,id)
              zvp(1,korpg(ik+1,ir)) = zvp(4,id)
              WRITE(50,'(5X,2A,2I4,F10.4,A)') '   :',
     .          ' Splitting side: ',ik,ir,t,' outer target'
            ENDDO
          ENDIF
c         Delete excess cells:
          IF (ik1.LT.nks(ir).AND.ik2.LT.nks(ir)) THEN
            DO ii = nks(ir), MAX(ik1,ik2)+1, -1
              CALL DeleteCell(nks(ir),ir)
              WRITE(50,'(3X,2A,I3,A,I3,A,I4,A,I4)') 'ShapeTarget:',
     .        ' DELETING CELL  ',ir,' :  target ',2,
     .        ' :  korpg ',id,' :  nks ',nks(ir)
            ENDDO
          ENDIF
        ELSE
          CALL ER('ShapeTarget','Unrecognized mode',*99)
        ENDIF
c Sep 22, 97 - Should recaluate BRATIO here as well, based
c on linear extrapolation...
      ENDDO


      RETURN
99    CONTINUE
      WRITE(0,*) 'IR,IDRING=',ir,id,mode
      WRITE(0,*) 'IRS,IRE,IRSEP,NRS=',irs,ire,irsep,nrs
      DO i1 = 1, 4
        rvertp(i1,:) = SNGL(rvp(i1,:))
        zvertp(i1,:) = SNGL(zvp(i1,:))
      ENDDO
      CALL DumpGrid('ERROR IN ShapeTarget')
      STOP
      END
c
c ========================================================================
c
c subroutine: SequenceTargets
c
      SUBROUTINE SequenceTargets(nregion1,nlist1,ilist1,
     .                           nregion2,nlist2,ilist2)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER ir,i1,i2,ir1,ir2,id1,id2,id,
     .        nregion1,nregion2,
     .        nlist1(0:MAXNRS),nlist2(0:MAXNRS),
     .        ilist1(MAXNRS,0:MAXNRS),ilist2(MAXNRS,0:MAXNRS)
      LOGICAL cont,found,next

      REAL       TOL
c
c     jdemod - change tol to 1.0e-5 for some high resolution grids I was using (sub-mm at midplane)
c
      PARAMETER (TOL=1.0E-05)
c      PARAMETER (TOL=1.0E-04)


c      CALL DUMPGRID('BUUMMMER')

      nregion1    = 1
      nlist1(1)   = 1 
      ilist1(1,1) = irsep
      ir1 = irsep
      ir2 = irsep
      cont = .TRUE.
      DO WHILE(cont)
        cont  = .FALSE.
        found = .FALSE.

c...    Low index targets:

c...    Look for a neighbouring target segment on vertex 2 side
c       of the low IK target segment:
        DO ir = irsep, nrs
          IF (ir1.EQ.ir.OR.idring(ir).EQ.BOUNDARY.OR.
     .        ir1.LT.irwall.AND.ir.GT.irtrap.OR.
     .        ir1.GT.irtrap.AND.ir.LT.irwall) CYCLE
          id1 = korpg(1,ir1)
          id2 = korpg(1,ir )
c          WRITE(0,*) 'VERT:',ir1,rvertp(2,id1),ir,rvertp(1,id2)
c          WRITE(0,*) '    :',ir1,zvertp(2,id1),ir,zvertp(1,id2)
          IF (ABS(rvertp(2,id1)-rvertp(1,id2)).LT.TOL.AND.
     .        ABS(zvertp(2,id1)-zvertp(1,id2)).LT.TOL) THEN
            nlist1(nregion1) = nlist1(nregion1) + 1
            ilist1(nlist1(nregion1),nregion1) = ir
            ir1 = ir
            found = .TRUE.
c            WRITE(0,*) 'IR1:',ir1,nregion1
          ENDIF
        ENDDO
c...    Look for a neighbouring target segment on the vertex 1 side 
c       of the low IK target segment:
        IF (cgridopt.NE.LINEAR_GRID.AND.cgridopt.NE.RIBBON_GRID) THEN
          DO ir = nrs, irsep, -1
            IF (ir2.EQ.ir.OR.idring(ir).EQ.BOUNDARY.OR.
     .          (ir2.LT.irwall.AND.ir.GT.irtrap.OR.
     .           ir2.GT.irtrap.AND.ir.LT.irwall)) CYCLE
            id1 = korpg(1,ir2)
            id2 = korpg(1,ir)
            IF (ir2.EQ.irsep.and.(ir.EQ.nrs.or.ir.eq.irsep+1)) THEN
c              WRITE(0,'(A,2F10.6,3I6)') 
c     .   'X-->',ABS(rvertp(1,id1)-rvertp(2,id2)),
c     .         ABS(zvertp(1,id1)-zvertp(2,id2)), 
c     .         nlist1(nregion1),ir2,ir
            ENDIF
            IF (ir2.EQ.irsep+1.and.(ir.EQ.irsep+2.or.ir.eq.irsep)) THEN
c              WRITE(0,'(A,2F10.6,3I6)') 
c     .   'Y-->',ABS(rvertp(1,id1)-rvertp(2,id2)),
c     .         ABS(zvertp(1,id1)-zvertp(2,id2)), 
c     .         nlist1(nregion1),ir2,ir
            ENDIF
            IF (ABS(rvertp(1,id1)-rvertp(2,id2)).LT.TOL.AND.
     .          ABS(zvertp(1,id1)-zvertp(2,id2)).LT.TOL) THEN
c...          Make room:
              DO i1 = nlist1(nregion1), 1, -1
                ilist1(i1+1,nregion1) = ilist1(i1,nregion1)
              ENDDO
              nlist1(nregion1) = nlist1(nregion1) + 1
              ilist1(1,nregion1) = ir
c              WRITE(0,*) 'FOUND',ir,ir2
              ir2 = ir
              found = .TRUE.
c              WRITE(0,*) 'FOUND',ir,ir2
            ENDIF
          ENDDO
        ENDIF
c            WRITE(0,*) 'DONE LOOP'
        IF (found) THEN
c...      Check if a connecting target segment was found (but was
c         not sequential in terms of the target index):
          cont = .TRUE.          
        ELSE
c...      Find the innermost target segment that is not yet included
c         in the identified target regions:
          DO ir = irsep, nrs          
            IF (idring(ir).EQ.BOUNDARY) CYCLE
            next = .TRUE.
            DO i1 = 1, nregion1
              DO i2 = 1, nlist1(i1)
                IF (ilist1(i2,i1).EQ.ir) next = .FALSE.
              ENDDO
            ENDDO
            IF (next) THEN
              ir1 = ir
              ir2 = ir
              nregion1 = nregion1 + 1
              nlist1(nregion1) = 1
              ilist1(1,nregion1) = ir1
              cont = .TRUE.
c              WRITE(0,*) 'NEXT',ir1,nregion1
              EXIT
            ENDIF 
          ENDDO
        ENDIF

      ENDDO

c      WRITE(0,*) 'SEQUENCE LOW :'
c      DO i1 = 1, nregion1
c          WRITE(0,'(I6,I6,50(I6:))') 
c     .  i1,nlist1(i1),(ilist1(i2,i1),i2=1,nlist1(i1))
c      ENDDO


c...  HIGH INDEX TARGETS:

c          WRITE(0,*) 'ID1:'
c          WRITE(0,*) 'ID1:',nrs,nks(irtrap+1),irtrap+1

      nregion2    = 1
      nlist2(1)   = 1 
      ilist2(1,1) = irsep
      ir1 = irsep
      ir2 = irsep
      cont = .TRUE.
      DO WHILE(cont)
        cont  = .FALSE.
        found = .FALSE.

        DO ir = irsep, nrs
c        DO ir = irsep, irwall-1 !     Not including PFZ for now:
          IF (ir1.EQ.ir.OR.idring(ir).EQ.BOUNDARY.OR.
     .        ir1.LT.irwall.AND.ir.GT.irtrap.OR.
     .        ir1.GT.irtrap.AND.ir.LT.irwall) CYCLE

c          WRITE(0,*) 'ID1:',nks(ir1),ir1,irwall

          id1 = korpg(nks(ir1),ir1)
          id2 = korpg(nks(ir),ir)
          IF (ABS(rvertp(3,id1)-rvertp(4,id2)).LT.TOL.AND.
     .        ABS(zvertp(3,id1)-zvertp(4,id2)).LT.TOL) THEN
            nlist2(nregion2) = nlist2(nregion2) + 1
            ilist2(nlist2(nregion2),nregion2) = ir
            ir1 = ir
            found = .TRUE.
          ENDIF
        ENDDO
        IF (cgridopt.NE.LINEAR_GRID.AND.cgridopt.NE.RIBBON_GRID) THEN
          DO ir = nrs, irsep, -1
c          DO ir = irwall-1, irsep, -1
            IF (ir2.EQ.ir.OR.idring(ir).EQ.BOUNDARY.OR.
     .          ir2.LT.irwall.AND.ir.GT.irtrap.OR.
     .          ir2.GT.irtrap.AND.ir.LT.irwall) CYCLE
            id1 = korpg(nks(ir2),ir2)
            id2 = korpg(nks(ir),ir)
            IF (ABS(rvertp(4,id1)-rvertp(3,id2)).LT.TOL.AND.
     .          ABS(zvertp(4,id1)-zvertp(3,id2)).LT.TOL) THEN
c...          Make room in the list:
              DO i1 = nlist2(nregion2), 1, -1
                ilist2(i1+1,nregion2) = ilist2(i1,nregion2)
              ENDDO
              nlist2(nregion2) = nlist2(nregion2) + 1
              ilist2(1,nregion2) = ir
              ir2 = ir
              found = .TRUE.
            ENDIF
          ENDDO
        ENDIF
        IF (found) THEN
c...      Check if a connecting target segment was found (but was
c         not sequential in terms of the target index):
          cont = .TRUE.          
        ELSE
c...      Find the innermost target segment that is not yet included
c         in the list of target regions:
          DO ir = irsep, nrs
            IF (idring(ir).EQ.BOUNDARY) CYCLE
c          DO ir = irsep, irwall-1          
            next = .TRUE.
            DO i1 = 1, nregion2
              DO i2 = 1, nlist2(i1)
                IF (ilist2(i2,i1).EQ.ir) next = .FALSE.
              ENDDO
            ENDDO
            IF (next) THEN
              ir1 = ir
              ir2 = ir
              nregion2 = nregion2 + 1
              nlist2(nregion2) = 1
              ilist2(1,nregion2) = ir1
              cont = .TRUE.
              EXIT
            ENDIF 
          ENDDO
        ENDIF

      ENDDO

c      WRITE(0,*) 'SEQUENCE HIGH:'
c      DO i1 = 1, nregion2
c          WRITE(0,'(I6,I6,50(I6:))') 
c     .  i1,nlist2(i1),(ilist2(i2,i1),i2=1,nlist2(i1))
c      ENDDO

      RETURN
99    STOP
      END
c
c ========================================================================
c
c subroutine: SequenceGrid
c
      SUBROUTINE SequenceGrid
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      REAL ATAN3C

      INTEGER ir,i1,i2,i3,i4,ir1,ir2,nlist,ilist(0:MAXNRS),id1,id2,id,
     .        nregion1,nregion2,nlist1(0:MAXNRS),nlist2(0:MAXNRS),
     .        ilist1(MAXNRS,0:MAXNRS),ilist2(MAXNRS,0:MAXNRS),fp
      LOGICAL cont,status,found,debug
      REAL    deltar,deltaz,angle1(0:MAXNRS),angle2(0:MAXNRS)

      REAL       TOL
c
c     jdemod - changed tolerance to 1.0e-5 since I was using some grids with sub-mm resolution at mid-plane
c
      PARAMETER (TOL=1.0E-05)
C      PARAMETER (TOL=1.0E-04)
c      PARAMETER (TOL=1.0E-03)

      debug = .TRUE.
      fp    = PINOUT
 
c...  Remove small gaps from between adjacent targets, assuming that 
c     such small gaps are the result of grid irregularities, rather than 
c     representing discontinuities in the targets:

      DO ir1 = irsep, nrs
        IF (idring(ir1).EQ.BOUNDARY) CYCLE
        DO ir2 = irsep, nrs
          IF (idring(ir2).EQ.BOUNDARY.OR.ir1.EQ.ir2) CYCLE
c...      Check low IK index target:
          id1 = korpg(1,ir1)
          id2 = korpg(1,ir2)          
          IF (rvertp(1,id1).NE.rvertp(2,id2).AND.
     .        zvertp(1,id1).NE.zvertp(2,id2).AND.
     .        ABS(rvertp(1,id1)-rvertp(2,id2)).LT.TOL.AND.
     .        ABS(zvertp(1,id1)-zvertp(2,id2)).LT.TOL) THEN 
            rvertp(1,id1) = rvertp(2,id2)
            zvertp(1,id1) = zvertp(2,id2)
            WRITE(0,*) 'DEBUG: PROBLEM LOW  ',ir1,ir2,id1,id2
          ENDIF
c...      Check high IK index target:
          id1 = korpg(nks(ir1),ir1)
          id2 = korpg(nks(ir2),ir2)          
          IF (rvertp(4,id1).NE.rvertp(3,id2).AND.
     .        zvertp(4,id1).NE.zvertp(3,id2).AND.
     .        ABS(rvertp(4,id1)-rvertp(3,id2)).LT.TOL.AND.
     .        ABS(zvertp(4,id1)-zvertp(3,id2)).LT.TOL) THEN
            rvertp(4,id1) = rvertp(3,id2)
            zvertp(4,id1) = zvertp(3,id2)
            WRITE(0,*) 'DEBUG: PROBLEM HIGH ',ir1,ir2
          ENDIF
        ENDDO
      ENDDO

c...  Find contiguous targets (this may not be the most straightfoward method of 
c     doing this, but I am looking for generality at present):
      CALL SequenceTargets(nregion1,nlist1,ilist1,
     .                     nregion2,nlist2,ilist2)

c      WRITE(0,*) 'PST:',nrs,irtrap+1,nks(irtrap+1)

      IF (debug) THEN
        WRITE(fp,*)
        WRITE(fp,*) 'Target grouping A:'
        DO i1 = 1, nregion1
          WRITE(fp,'(I6,I6,50(I6:))') 
     .      i1,nlist1(i1),(ilist1(i2,i1),i2=1,nlist1(i1))
        ENDDO
        WRITE(fp,*)
        DO i1 = 1, nregion2
          WRITE(fp,'(I6,I6,50(I6:))') 
     .      i1,nlist2(i1),(ilist2(i2,i1),i2=1,nlist2(i1))
        ENDDO
      ENDIF

c...  Now, determine the preferred order of the target regions, so that
c     the rings can be sorted properly:

c...  Find a poloidal angular coordinate for each low index target region, so that
c     they can be properly ordered:
      DO i1 = 1, nregion1
        ir = ilist1(MAX(nlist1(i1)/2,1),i1)
        id = korpg(1,ir)
        IF (cgridopt.EQ.LINEAR_GRID.OR.cgridopt.EQ.RIBBON_GRID) THEN 
          angle1(i1) = zvertp(1,id)
        ELSE
          deltar = 0.5 * (rvertp(1,id) + rvertp(2,id)) - rxp
          deltaz = 0.5 * (zvertp(1,id) + zvertp(2,id)) - zxp
c          deltar = 0.5 * (rvertp(1,id) + rvertp(2,id)) - r0
c          deltaz = 0.5 * (zvertp(1,id) + zvertp(2,id)) - z0
          angle1(i1) = ATAN3C(deltaz,deltar)
          IF (i1.GT.1.AND.angle1(i1).GT.angle1(1)) 
     .      angle1(i1) = angle1(i1) - 360.0
        ENDIF
      ENDDO
c...  Sort target regions clockwise, with the first region containing
c     the separatrix:
      cont = .TRUE.
      DO WHILE(cont)
        cont = .FALSE.
        DO i1 = 2, nregion1-1
          IF (angle1(i1).LT.angle1(i1+1)) THEN
            nlist1(0) = nlist1(i1)
            angle1(0) = angle1(i1) 
            DO i2 = 1, nlist1(i1)
              ilist1(i2,0) = ilist1(i2,i1)
            ENDDO
            nlist1(i1) = nlist1(i1+1)
            angle1(i1) = angle1(i1+1) 
            DO i2 = 1, nlist1(i1+1)
              ilist1(i2,i1) = ilist1(i2,i1+1)
            ENDDO
            nlist1(i1+1) = nlist1(0)
            angle1(i1+1) = angle1(0) 
            DO i2 = 1, nlist1(0)
              ilist1(i2,i1+1) = ilist1(i2,0)
            ENDDO
            cont = .TRUE.
          ENDIF
        ENDDO
      ENDDO

      IF (debug) THEN
        WRITE(fp,*)
        WRITE(fp,*) 'Target grouping B:'
        DO i1 = 1, nregion1
          WRITE(fp,'(I6,I6,50(I6:))') 
     .      i1,nlist1(i1),(ilist1(i2,i1),i2=1,nlist1(i1))
        ENDDO
        WRITE(fp,*)
        DO i1 = 1, nregion2
          WRITE(fp,'(I6,I6,50(I6:))') 
     .      i1,nlist2(i1),(ilist2(i2,i1),i2=1,nlist2(i1))
        ENDDO
      ENDIF

c...  Move rings around:
      IF (.TRUE.) THEN
c...    Sort rings based on the low IK index regions:
        ir1 = irsep
        DO i1 = 1, nregion1
          DO i2 = 1, nlist1(i1)
            ir2 = ilist1(i2,i1) 
            IF     (ir1.LT.irwall.AND.ir2.GT.irtrap.OR.
     .              ir1.GT.irtrap.AND.ir2.LT.irwall) THEN
c...          No SOL/PFZ mixing:
            ELSEIF (ir1.NE.ir2) THEN
c              IF (sloutput) THEN
c                WRITE(SLOUT,*) 'SWAPPING RINGS IN SEQUENCEGRID',ir1,ir2
c                WRITE(0    ,*) 'SWAPPING RINGS IN SEQUENCEGRID',ir1,ir2
c              ENDIF
              CALL CopyRing(ir1,nrs+1)
              CALL CopyRing(ir2,ir1)            
              CALL CopyRing(nrs+1,ir2)
c...          Replace IR1 in the list with IR2:
              DO i3 = 1, nregion1
                DO i4 = 1, nlist1(i3)
                  IF     (ilist1(i4,i3).EQ.ir1) THEN
                    ilist1(i4,i3) = ir2
                  ELSEIF (ilist1(i4,i3).EQ.ir2) THEN
                    ilist1(i4,i3) = ir1
                  ENDIF
                ENDDO
              ENDDO
            ENDIF

            ir1 = ir1 + 1 ! DEV
            DO WHILE(idring(ir1).EQ.BOUNDARY)
              ir1 = ir1 + 1
            ENDDO
          ENDDO
        ENDDO
      ENDIF



c      WRITE(0,*) 'PST:',nrs,irtrap+1,nks(irtrap+1)

c...  Identify target regions again, with the reorganized grid:
      CALL SequenceTargets(nregion1,nlist1,ilist1,
     .                     nregion2,nlist2,ilist2)

c...  Sort high IK (outer) index target regions counter-clockwise.  The low IK
c     index (inner) target regions should be fine, since the rings were sorted
c     based on the inner target regions:

c...  Find a poloidal angular coordinate for each high index target region:
      DO i1 = 1, nregion2
        ir = ilist2(MAX(nlist2(i1)/2,1),i1)
        id = korpg(nks(ir),ir)
        IF (cgridopt.EQ.LINEAR_GRID.OR.cgridopt.EQ.RIBBON_GRID) THEN 
          angle2(i1) = zvertp(1,id)
        ELSE
          deltar = 0.5 * (rvertp(3,id) + rvertp(4,id)) - rxp
          deltaz = 0.5 * (zvertp(3,id) + zvertp(4,id)) - zxp
c          deltar = 0.5 * (rvertp(3,id) + rvertp(4,id)) - r0
c          deltaz = 0.5 * (zvertp(3,id) + zvertp(4,id)) - z0
          angle2(i1) = ATAN3C(deltaz,deltar)
          DO WHILE (i1.NE.1.AND.angle2(i1).LT.angle2(1))
            angle2(i1) = angle2(i1) + 360.0
          ENDDO
c          IF (angle2(i1).LT.angle2(1)) angle2(i1) = angle2(i1) + 360.0
c          WRITE(0,*) 'ANGLE:',i1,angle2(i1)
        ENDIF
      ENDDO
c...  Sort the region indeces:
      cont = .TRUE.
      DO WHILE(cont)
        cont = .FALSE.
        DO i1 = 2, nregion2-1
          IF (angle2(i1).GT.angle2(i1+1)) THEN
            nlist2(0) = nlist2(i1)
            angle2(0) = angle2(i1) 
            DO i2 = 1, nlist2(i1)
              ilist2(i2,0) = ilist2(i2,i1)
            ENDDO
            nlist2(i1) = nlist2(i1+1)
            angle2(i1) = angle2(i1+1) 
            DO i2 = 1, nlist2(i1+1)
              ilist2(i2,i1) = ilist2(i2,i1+1)
            ENDDO
            nlist2(i1+1) = nlist2(0)
            angle2(i1+1) = angle2(0) 
            DO i2 = 1, nlist2(0)
              ilist2(i2,i1+1) = ilist2(i2,0)
            ENDDO
            cont = .TRUE.
          ENDIF
        ENDDO
      ENDDO

      IF (debug) THEN
        WRITE(fp,*)
        WRITE(fp,*) 'Target grouping C:'
        DO i1 = 1, nregion1
          WRITE(fp,'(I6,I6,50(I6:))') 
     .      i1,nlist1(i1),(ilist1(i2,i1),i2=1,nlist1(i1))
        ENDDO
        WRITE(fp,*)
        DO i1 = 1, nregion2
          WRITE(fp,'(I6,I6,50(I6:))') 
     .      i1,nlist2(i1),(ilist2(i2,i1),i2=1,nlist2(i1))
        ENDDO
      ENDIF

c...  Add IRWALL:
c      WRITE(0,*) 'IRWALL:',irwall,irtrap

      found = .FALSE.
      DO i1 = 1, nregion1
        i2 = nlist1(i1)
        IF (ilist1(i2,i1).EQ.irwall-1) THEN
          nlist1(i1) = nlist1(i1) + 1
          ilist1(i2+1,i1) = irwall
          found = .TRUE.
c          WRITE(0,*) 'FOUND1:',i1
        ENDIF
      ENDDO
      IF (.NOT.found) 
     .  CALL ER('SequenceGrid','Unable to add IRWALL 1',*98)

      found = .FALSE.
      DO i1 = 1, nregion2
        i2 = nlist2(i1)
        IF (ilist2(i2,i1).EQ.irwall-1) THEN
          nlist2(i1) = nlist2(i1) + 1
          ilist2(i2+1,i1) = irwall
          found = .TRUE.
c          WRITE(0,*) 'FOUND2:',i1
        ENDIF
      ENDDO
      IF (.NOT.found) 
     .  CALL ER('SequenceGrid','Unable to add IRWALL 2',*98)

c...  Add TRAP:
      found = .FALSE.
      DO i1 = 1, nregion1
        IF (ilist1(1,i1).EQ.irtrap+1) THEN
          DO i2 = nlist1(i1), 1, -1 
            ilist1(i2+1,i1) = ilist1(i2,i1)
          ENDDO
          nlist1(i1) = nlist1(i1) + 1
          ilist1(1,i1) = irtrap
          found = .TRUE.
        ENDIF
      ENDDO
      IF (.NOT.found.AND.irtrap.LT.nrs)
     .  CALL ER('SequenceGrid','Unable to add IRTRAP 1',*99)

      found = .FALSE.
      DO i1 = 1, nregion2
        IF (ilist2(1,i1).EQ.irtrap+1) THEN
          DO i2 = nlist2(i1), 1, -1 
            ilist2(i2+1,i1) = ilist2(i2,i1)
          ENDDO
          nlist2(i1) = nlist2(i1) + 1
          ilist2(1,i1) = irtrap
          found = .TRUE.
        ENDIF
      ENDDO
      IF (.NOT.found.AND.irtrap.LT.nrs) 
     .  CALL ER('SequenceGrid','Unable to add IRTRAP 1',*99)


      IF (sloutput) THEN
        WRITE(pinout,*)
        WRITE(pinout,*) 'Final target grouping:'
        DO i1 = 1, nregion1
          WRITE(pinout,'(I6,I6,50(I6:))') 
     .      i1,nlist1(i1),(ilist1(i2,i1),i2=1,nlist1(i1))
        ENDDO
        WRITE(pinout,*)
        DO i1 = 1, nregion2
          WRITE(pinout,'(I6,I6,50(I6:))') 
     .      i1,nlist2(i1),(ilist2(i2,i1),i2=1,nlist2(i1))
        ENDDO
        WRITE(pinout,*) 'SHOULD TARGET-TO-TARGET GRID REGIONS BE '//
     .                  'ASSIGNED?'
      ENDIF


c...  *** Can't do this because it doesn't work for a break in the PFZ itself, a la C-Mod... ***
c...  Quick check to make sure that the PFZ (primary) rings in REGION 2 (high
c     IK targets) are actually in the last group:
c      DO i1 = 1, nregion2 
c        DO i2 = 1, nlist2(i1)
c          IF ((i1.LT.nregion2.AND.ilist2(i2,i1).GT.irtrap).OR.
c     .        (i1.EQ.nregion2.AND.ilist2(i2,i1).LT.irtrap)) 
c     .      CALL ER('SequenceGrid','High IK target PFZ rings not in '//
c     .              'last group',*98)
c        ENDDO
c      ENDDO


      IF (.TRUE.) THEN
        WRITE(PINOUT,*)
        WRITE(PINOUT,*) 'TARGET GROUPS:'
        WRITE(PINOUT,*)
        DO i1 = 1, nregion1
            WRITE(PINOUT,'(I6,I6,50(I6:))') 
     .    i1,nlist1(i1),(ilist1(i2,i1),i2=1,nlist1(i1))
        ENDDO
        WRITE(PINOUT,*)
        DO i1 = 1, nregion2
            WRITE(PINOUT,'(I6,I6,50(I6:))') 
     .    i1,nlist2(i1),(ilist2(i2,i1),i2=1,nlist2(i1))
        ENDDO
      ENDIF

c...  Copy results into global arrays:
      DO i1 = 1, nregion1
        grdntreg(IKLO) = nregion1      
        grdntseg(i1,IKLO) = nlist1(i1)      
        DO i2 = 1, nlist1(i1)
          grdtseg(i2,i1,IKLO) = ilist1(i2,i1)
        ENDDO
      ENDDO
      DO i1 = 1, nregion2
        grdntreg(IKHI) = nregion2
        grdntseg(i1,IKHI) = nlist2(i1)      
        DO i2 = 1, nlist2(i1)
          grdtseg(i2,i1,IKHI) = ilist2(i2,i1)
        ENDDO
      ENDDO


c      STOP 'YEAH!'



      RETURN
 98   WRITE(0,*)
      WRITE(0,*) '  GRID SEQUENCING:'
      DO i1 = 1, nregion1
          WRITE(0,'(I6,I6,50(I6:))') 
     .  i1,nlist1(i1),(ilist1(i2,i1),i2=1,nlist1(i1))
      ENDDO
      WRITE(0,*)
      DO i1 = 1, nregion2
          WRITE(0,'(I6,I6,50(I6:))') 
     .  i1,nlist2(i1),(ilist2(i2,i1),i2=1,nlist2(i1))
      ENDDO
99    CALL DumpGrid('PROBLEMS SEQUENCING GRID')
      STOP
      END
c
c ======================================================================
c
c
c
c
c
      SUBROUTINE LoadGridData(mode)
      IMPLICIT none

      INTEGER mode

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER   fp,ik,ir,ik1,ir1
      CHARACTER buffer*1024
      REAL      rdum1,rdum2,rdum3,rdum4,diff1,diff2,diff3

      REAL, ALLOCATABLE :: bratio2(:,:), psin2(:,:)

      ALLOCATE(bratio2(MAXNKS,MAXNRS))
      ALLOCATE(psin2  (MAXNKS,MAXNRS))

      fp = 98
      SELECTCASE(mode)
        CASE (0,2,3,10)
          OPEN(UNIT=fp,FILE='grid.dat'  ,ACCESS='SEQUENTIAL',
     .         STATUS='OLD',ERR=96)      
        CASE (1)
          OPEN(UNIT=fp,FILE='grid.dat.1',ACCESS='SEQUENTIAL',
     .         STATUS='OLD',ERR=96)      
        CASE DEFAULT
          CALL ER('LoadGridData','Unknown MODE',*99)
      ENDSELECT

c...  Clear header:
      SELECTCASE (mode)
        CASE (10)
          READ(fp,*) 
          READ(fp,*) 
          DO ir = 1, nrs
            IF (idring(ir).EQ.BOUNDARY) CYCLE
            DO ik = 1, nks(ir)
              READ(fp,*,END=98) rdum1,rdum2,rdum3,rdum4,ik1,ir1
              IF (ik1.EQ.ik.AND.ir1.EQ.ir) THEN
                bratio2(ik,ir) = rdum3        
              ELSE
                CALL ER('LoadGridData','IK/IR index does not match',*99)
              ENDIF
            ENDDO
          ENDDO
        CASE (0:1)
          READ(fp,*) 
          READ(fp,*) 
          DO WHILE(.TRUE.)
            READ(fp,*,END=10) rdum1,rdum2,rdum3,rdum4,ik,ir
            bratio2(ik,ir) = rdum3        
            psin2(ik,ir) = rdum4
          ENDDO
 10       CONTINUE
        CASE (2)
          WRITE(0,*) 'MAKE SURE THE .SUP FILE HAS A VERSION NUMBER'
          WRITE(buffer,'(1024X)')
          READ(fp,*) 
          buffer(1:1) = '*'
          DO WHILE(buffer(1:1).EQ.'*')
            READ(fp,'(A1024)') buffer
c            WRITE(0,*) 'WHAT'//buffer(1:1)//'<'
          ENDDO
          BACKSPACE(fp)
        CASE (3)
          READ(fp,*) 
          WRITE(buffer,'(1024X)')
          buffer(1:1) = '*'
          DO WHILE(buffer(1:1).EQ.'*')
            READ(fp,'(A1024)') buffer
          ENDDO
          BACKSPACE(fp)
      ENDSELECT

c...  Compare:
c      WRITE(SLOUT,*)
c      WRITE(SLOUT,*) 'BRATIO RECACULATION:'
c      DO ir = 1, nrs
c        IF (idring(ir).EQ.BOUNDARY) CYCLE
c        DO ik = 1, nks(ir)
c          diff1 = (bratio2(ik,ir) - bratio(ik,ir)) / bratio(ik,ir)*100.0
c          diff2 = ((psin2(ik,ir)-1.0) - (psitarg(ir,2)-1.0)) / 
c     .            (psitarg(ir,2)-1.0)*100.0
c          diff3 = ((psin2(ik,ir)-1.0) - (psitarg(ir,1)-1.0)) / 
c     .            (psitarg(ir,1)-1.0)*100.0
c          WRITE(SLOUT,'(2I6,2F14.6,F8.3,A,2X,3F10.4,2(F8.3,A))') 
c     .      ik,ir,bratio(ik,ir),bratio2(ik,ir),diff1,'%',
c     .      psin2(ik,ir),psitarg(ir,2),psitarg(ir,1),
c     .      diff2,'%  ',diff3,'%'
c        ENDDO
c      ENDDO

c...  Assign grid values:
      psitarg = 0.0
c      bratio = 0.0

      SELECTCASE (mode)
        CASE (10) 
          DO WHILE (.TRUE.)
            READ(fp,*,END=20) rdum1,rdum2,rdum3,rdum4,ik,ir
            IF (ir.LT.irsep) THEN
              psitarg(ir,1) = rdum4
              psitarg(ir,2) = psitarg(ir,1)
            ELSE
              IF (ik.EQ.1      ) psitarg(ir,2) = rdum4
              IF (ik.EQ.nks(ir)) psitarg(ir,1) = rdum4
            ENDIF
          ENDDO
 20       CONTINUE
          DO ir = 1, nrs
            IF (idring(ir).EQ.BOUNDARY) CYCLE
            DO ik = 1, nks(ir)
              IF (bratio2(ik,ir).LE.0.0) 
     .          CALL ER('LoadGridData','BRATIO.LE.0.0',*99)
            ENDDO
            bratio(1:nks(ir),ir) = bratio2(1:nks(ir),ir)
          ENDDO
        CASE (0:1)
          DO ir = 1, nrs
            IF (idring(ir).EQ.BOUNDARY) CYCLE
            psitarg(ir,2) = psin2(1,ir)
            psitarg(ir,1) = psin2(nks(ir),ir)
            DO ik = 1, nks(ir)
              IF (bratio2(ik,ir).NE.0.0) bratio(ik,ir) = bratio2(ik,ir)
            ENDDO
          ENDDO
        CASE (2) 
          DO ir = 1, nrs
            IF (idring(ir).EQ.BOUNDARY) CYCLE
            DO ik = 1, nks(ir)
              READ(fp,'(A1024)',END=25) buffer
c              WRITE(0,*) 'BUFFER',buffer(1:100)
              READ(buffer,*) ik1,ir1,rdum1
              IF (ik.NE.ik1.OR.ir.NE.ir1) 
     .          CALL ER('LoadGridData','Cell index mismatch',*99)
              IF (ik.EQ.1      ) psitarg(ir,2) = rdum1
              IF (ik.EQ.nks(ir)) psitarg(ir,1) = rdum1
            ENDDO
          ENDDO
 25       CONTINUE
        CASE (3) 
          DO ir = 1, nrs
            IF (idring(ir).EQ.BOUNDARY) CYCLE
            DO ik = 1, nks(ir)
              READ(fp,'(A1024)',END=30) buffer
              READ(buffer,*) ik1,ir1,rdum1,rdum2
              IF (ik.NE.ik1.OR.ir.NE.ir1) 
     .          CALL ER('LoadGridData','Cell index mismatch',*99)
              IF (ik.EQ.1      ) psitarg(ir,2) = rdum1
              IF (ik.EQ.nks(ir)) psitarg(ir,1) = rdum1
              bratio(ik,ir) = rdum2
              kbfs  (ik,ir) = 1.0 / (rdum2 + 1.0E-10)
            ENDDO
          ENDDO
 30       CONTINUE
      ENDSELECT

c...  Close stream:
      CLOSE (fp)
c...  Clear arrays:
      DEALLOCATE(bratio2)
      DEALLOCATE(psin2)

      RETURN
 96   WRITE(0,*) 'LoadGridData: Cannot find data file'
      GOTO 99
 98   CALL ER('LoadGridData','Error reading supplimental grid data '//
     .        'file',*99)
 99   WRITE(0,*) 'MODE   =',mode
      WRITE(0,*) 'IK,IR  =',ik ,ir
      WRITE(0,*) 'IK1,IR1=',ik1,ir1
      STOP
      END
c
c ======================================================================
c
      SUBROUTINE DumpGrid(note)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
 
      CHARACTER*(*) note

c...  For call to STORE:
      CHARACTER title*174,desc*1024,job*72,equil*60
      REAL      facta(-1:MAXIZS),factb(-1:MAXIZS)

      INTEGER i1,i2,ik,ir,id,in,iwall

c...  Recalculate cell centers:
c      DO ir = 1, nrs
c        DO ik = 1, nks(ir)
c          id = korpg(ik,ir)
c          rs(ik,ir) = 0.0
c          zs(ik,ir) = 0.0
c          DO in = 1, nvertp(id)
c            rs(ik,ir) = rs(ik,ir) + 1.0/REAL(nvertp(id))*rvertp(in,id)
c            zs(ik,ir) = zs(ik,ir) + 1.0/REAL(nvertp(id))*zvertp(in,id)
c          ENDDO
c        ENDDO
c      ENDDO

c      CALL SetupGrid

      IF (nvesm.EQ.0) THEN
        nvesm = nwall
        DO iwall = 1, nwall
          jvesm(iwall) = 8
          rvesm(iwall,1) = wallco(iwall ,1)
          zvesm(iwall,1) = wallco(iwall ,2)
          rvesm(iwall,2) = wallco(iwall+1,1)
          zvesm(iwall,2) = wallco(iwall+1,2)
        ENDDO 
      ENDIF

      CALL SaveSolution
      CALL OutputData(86,note(1:LEN_TRIM(note)))

      title = note(1:LEN_TRIM(note))
      desc  = 'Call to STORE from DumpGrid'
      job   = 'Call to STORE from DumpGrid'
      equil = 'Call to STORE from DumpGrid'

      WRITE(0,*) 'CALLING STORE'
      CALL Store(title,desc,1,job,equil,facta,factb,1,1)

      WRITE(0,*) 'HALTING CODE FROM DUMPGRID WHILE '//
     .           note(1:LEN_TRIM(note))
      STOP

      RETURN
99    STOP
      END
c     
c ========================================================================
c
c subroutine: CalcPolDist
c
c 
c
c
c
      SUBROUTINE CalcPolDist(ir,pdist,ploc)
      USE mod_grid_divimp
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER ik,ir,id,ilocmax,ikend
      REAL    pdist(MAXNKS),ploc(0:MAXNKS),r1,r2,z1,z2,plocshift

      ikend = nks(ir)
      IF (ir.LT.irsep) ikend = ikend - 1

      DO ik = 1, ikend
        id = korpg(ik,ir)
        IF (ALLOCATED(d_rvertp)) THEN
          r1 = 0.5 * SNGL(d_rvertp(1,id) + d_rvertp(2,id))
          z1 = 0.5 * SNGL(d_zvertp(1,id) + d_zvertp(2,id))
          r2 = 0.5 * SNGL(d_rvertp(3,id) + d_rvertp(4,id))
          z2 = 0.5 * SNGL(d_zvertp(3,id) + d_zvertp(4,id))
        ELSE
          r1 = 0.5 * (rvertp(1,id) + rvertp(2,id))
          z1 = 0.5 * (zvertp(1,id) + zvertp(2,id))
          r2 = 0.5 * (rvertp(3,id) + rvertp(4,id))
          z2 = 0.5 * (zvertp(3,id) + zvertp(4,id))
        ENDIF
        pdist(ik) = SQRT((r1 - r2)**2 + (z1 - z2)**2)
        IF (ik.EQ.1) THEN
          ploc(ik) = 0.5 * pdist(ik)
        ELSE
          ploc(ik) = ploc(ik-1) + 0.5 * (pdist(ik-1) + pdist(ik))
        ENDIF
      ENDDO
      ploc(0) = ploc(ikend) + 0.5 * pdist(nks(ir))
c...  Set coordinate so that it is increasing towards the center from
c     both targets:
      DO ik = 1, ikend
        ploc(ik) = ABS(ABS(ploc(ik) - 0.5 * ploc(0)) - 0.5 * ploc(0))
      ENDDO
c...  Adjust poloidal coordinate so that the 1st and last cells are at 
c     position 0.0:
      plocshift = ploc(1)
      DO ik = 1, ikend-1
        IF (ploc(ik).LT.ploc(ik+1)) ploc(ik) = ploc(ik) - plocshift
      ENDDO
      plocshift = ploc(ikend)
      DO ik = ikend, 2, -1
        IF (ploc(ik).LT.ploc(ik-1)) ploc(ik) = ploc(ik) - plocshift
      ENDDO
      ploc(0      ) = 0.0
      ploc(ikend+1) = 0.0

      RETURN
 99   STOP
      END
c     
c ========================================================================
c ========================================================================
c
c subroutine: PoloidalRefinement
c
c 
c
c
c
      SUBROUTINE PoloidalRefinement(ir1,mode,param)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER, INTENT(IN) :: ir1,mode
      REAL   , INTENT(IN) ::  param

      INTEGER i1,ik,ir,status,nks2,iks,ike,ikto3,ikti3,lastmode,lastir,
     .        ikmark1(MAXNRS),ikmark2(MAXNRS),tmpnmod
      LOGICAL cont
      REAL    pdist(MAXNKS),ploc(0:MAXNKS),tloc,tdist,lastparam

      DATA lastmode,lastir,lastparam /0,0,0.0/

      SAVE

      ir = ir1

      IF (mode.GE.1.AND.mode.LE.3) THEN
c...    Simple refinement of low IK index near target region:
    
        cont = .TRUE.
        DO WHILE (cont)
          cont = .FALSE.
          nks2 = nks(ir)
          IF (ir.LT.irsep) nks2 = nks2 - 1
          DO ik = nks2, 1, -1
c...        Calculate poloidal distribution of cells:
            CALL CalcPolDist(ir,pdist,ploc)

            IF (ploc(ik)/param.GT.5.0) THEN
              tdist = 200.0
            ELSE
              tdist = 1.5 * grd_minpl * EXP(ploc(ik)/param)
            ENDIF

            WRITE(88,'(A,2I6,3F12.5,L2)') 
     .        'TDIST:',ik,ir1,tdist,pdist(ik),grd_minpl,
     .         (tdist.LT.pdist(ik).AND.
     .          (mode.EQ.3.OR.
     .           mode.EQ.1.AND.ploc(ik).LT.ploc(ik+1).OR.
     .           mode.EQ.2.AND.ploc(ik).LT.ploc(ik-1)))
                     
            IF (tdist.LT.pdist(ik).AND.
     .          (mode.EQ.3.OR.
     .           mode.EQ.1.AND.ploc(ik).LT.ploc(ik+1).OR.
     .           mode.EQ.2.AND.ploc(ik).LT.ploc(ik-1))) THEN

              CALL SplitCell(ik,ir,0.5D0,status)
c              WRITE(0,*) '  POLOIDAL:',ik,ir,status
              IF (status.EQ.-1) THEN
                CALL ER('PoloidalRefinement','Unable to refine grid',
     .                  *99)
              ELSEIF (status.EQ.1) THEN
c...            Soft fail:
                WRITE(0,*) 'POLOIDAL REFINEMENT: SOFT FAIL'                
              ELSE
                cont = .TRUE.
              ENDIF

            ENDIF
c            STOP 'sdfsd'

          ENDDO    
 
        ENDDO
c     -------------------------------------------------------------------
      ELSEIF ((mode.GE.4.AND.mode.LE.9).OR.mode.EQ.11) THEN
c...    Refine the divertor region:

c...    Triggers for resetting grid parameters:
        IF ((ir.NE.lastir+1.AND.ir.NE.irtrap+1).OR.
     .      mode.NE.lastmode.OR.lastparam.NE.param) THEN
          IF (mode.EQ.4.OR.mode.EQ.5) THEN
            ikto3 = NINT(REAL(ikto           )*param)
            ikti3 = NINT(REAL(nks(irsep)-ikti)*param)
          ELSE
            ikto3 = ikto
            ikti3 = nks(irsep) - ikti
          ENDIF
        ENDIF

        lastir = ir
        lastmode = mode
        lastparam = param

        IF     (mode.EQ.4) THEN
          iks = 1
          ike = ikto3
        ELSEIF (mode.EQ.5) THEN
          IF (ir.LT.irwall) THEN
            iks = nks(ir) - ikti3  
            ike = nks(ir)
          ELSE
            iks = nks(ir) - ikti3
            ike = nks(ir)
          ENDIF
        ELSEIF (mode.EQ.6) THEN
          iks = NINT(REAL(ikto3-1)*(1.0-param)) + 1
          ike = ikto3
        ELSEIF (mode.EQ.7) THEN
          iks = nks(ir) - ikti3
          ike = nks(ir) - NINT(REAL(ikti3)*(1.0-param))
        ELSEIF (mode.EQ.8) THEN
c...      Near target refinement in the inner divertor (outer on JET):
          iks = 1
          ike = MIN(nks(ir)/2,NINT(param))
        ELSEIF (mode.EQ.9) THEN
c...      Near target refinement in the outer divertor (inner on JET):
          iks = MAX(nks(ir)/2,nks(ir) - NINT(param) + 1)
          ike = nks(ir)
c        ELSEIF (mode.EQ.11) THEN
cc...      The whole ring:
c          iks = 1
c          ike = nks(ir)
        ENDIF

c        WRITE(0,*) 'IKS,IKE=',ir,mode,iks,ike

        CALL CalcPolDist(ir,pdist,ploc)
        DO ik = ike, iks, -1
c          IF (ir.EQ.irsep.AND.ik.LE.ikto3) ikto3 = ikto3 + 1
c          IF (ir.EQ.irsep.AND.ik.LT.ikti3) ikti3 = ikti3 + 1
          IF (grd_minpl.LT.pdist(ik)) THEN
            CALL SplitCell(ik,ir,0.5D0,status)
c            WRITE(0,*) '  POLOIDAL:',ik,ir,status
            IF (status.NE.0) 
     .        CALL ER('PoloidalRefinement','Unable to refine grid',
     .                *99)
          ENDIF
        ENDDO
c     -------------------------------------------------------------------
      ELSEIF (mode.EQ.10) THEN
c...    Refine the x-point:

c...    Hack:
        IF (grdmod(1,1).EQ.887.0) THEN
          tmpnmod = grdnmod
          grdnmod = 0
        ENDIF

        CALL BuildMap

        IF (grdmod(1,1).EQ.887.0) grdnmod = tmpnmod

        ikmark1(irsep) = ikti2(irsep)
        ikmark2(irsep) = ikto2(irsep)
        ikmark1(irsep+1) = ikouts(ikmark1(irsep),irsep)
        ikmark2(irsep+1) = ikouts(ikmark2(irsep),irsep)
        ikmark1(irsep+2) = ikouts(ikmark1(irsep+1),irsep+1)
        ikmark2(irsep+2) = ikouts(ikmark2(irsep+1),irsep+1)

        DO ir = irsep, irsep + 1
c        DO ir = irsep, irsep + 2
c...      Outer SOL (inside on JET):
          iks = ikmark1(ir) - NINT(param)
          ike = ikmark1(ir) + NINT(param) - 1
          DO ik = ike, iks, -1
            CALL SplitCell(ik,ir,0.5D0,status)
            IF (status.NE.0) 
     .        CALL ER('PoloidalRefinement','Unable to refine gridA',*99)
          ENDDO
c...      Inner SOL (outside on JET):
          iks = ikmark2(ir) - NINT(param) + 1
          ike = ikmark2(ir) + NINT(param)
          DO ik = ike, iks, -1
            CALL SplitCell(ik,ir,0.5D0,status)
            IF (status.NE.0) 
     .        CALL ER('PoloidalRefinement','Unable to refine gridB',*99)
          ENDDO
        ENDDO

        DO ir = nrs-2, nrs
c...      PFZ:
          iks = ikto2(ir) - NINT(param) + 1
          ike = ikti2(ir) + NINT(param) - 1
          DO ik = ike, iks, -1
            CALL SplitCell(ik,ir,0.5D0,status)
            IF (status.NE.0) 
     .        CALL ER('PoloidalRefinement','Unable to refine grid',*99)
          ENDDO
        ENDDO
c     -------------------------------------------------------------------
      ELSEIF (mode.EQ.12) THEN
c...    Double the number of cells on the ring (split each cell in half):
        DO i1 = 1, NINT(param)
          nks2 = nks(ir)
          DO ik = nks2, 1, -1
            CALL SplitCell(ik,ir,0.5D0,status)
            IF (status.EQ.-1)
     .        CALL ER('PoloidalRefinement','Unable to refine grid',*99)
          ENDDO    
        ENDDO
c     -------------------------------------------------------------------
      ELSEIF (mode.EQ.14) THEN
c...    Inner midplane:
        DO ik = nks(ir), 1, -1
          IF (rs(ik,ir).LT.r0.AND.ABS(zs(ik,ir)-z0).LT.param) THEN
            CALL SplitCell(ik,ir,0.5D0,status)
            IF (status.EQ.-1)
     .        CALL ER('PoloidalRefinement','Unable to refine grid',*99)
          ENDIF
        ENDDO
c     -------------------------------------------------------------------
      ELSEIF (mode.EQ.15) THEN
c...    Outer midplane:
        DO ik = nks(ir), 1, -1
          IF (rs(ik,ir).GT.r0.AND.ABS(zs(ik,ir)-z0).LT.param) THEN
            CALL SplitCell(ik,ir,0.5D0,status)
            IF (status.EQ.-1)
     .        CALL ER('PoloidalRefinement','Unable to refine grid',*99)
          ENDIF
        ENDDO
c     -------------------------------------------------------------------
      ELSE
        CALL ER('PoloidalRefinement','Invalid MODE',*99)
      ENDIF


      RETURN
 99   WRITE(0,'(A,2I6,F10.2)') 'DATA=',ir,mode,param
      STOP
      END
c
c ======================================================================
c
      REAL*8 FUNCTION VertexDisplacement(id1,iv1,id2,iv2)
      USE mod_grid_divimp
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'

      INTEGER, INTENT(IN) :: id1,id2,iv1,iv2

      REAL*8  a1,a2,b1,b2

      IF (ALLOCATED(d_rvertp)) THEN
        a1 = d_rvertp(iv1,id1)
        a2 = d_zvertp(iv1,id1)
        b1 = d_rvertp(iv2,id2)
        b2 = d_zvertp(iv2,id2)
      ELSE
        a1 = DBLE(rvertp(iv1,id1))
        a2 = DBLE(zvertp(iv1,id1))
        b1 = DBLE(rvertp(iv2,id2))
        b2 = DBLE(zvertp(iv2,id2))
      ENDIF

      VertexDisplacement = DSQRT((a1 - b1)**2 + (a2 - b2)**2)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      REAL*8 FUNCTION SideLength(id,iside)
      USE mod_grid_divimp
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'

      INTEGER, INTENT(IN) :: id,iside 

      INTEGER i1,i2
      REAL*8  a1,a2,b1,b2

      i1 = iside
      i2 = iside + 1
      IF (i2.EQ.nvertp(id)+1) i2 = 1 

      IF (ALLOCATED(d_rvertp)) THEN
        a1 = d_rvertp(i1,id)
        a2 = d_zvertp(i1,id)
        b1 = d_rvertp(i2,id)
        b2 = d_zvertp(i2,id)
      ELSE
        a1 = DBLE(rvertp(i1,id))
        a2 = DBLE(zvertp(i1,id))
        b1 = DBLE(rvertp(i2,id))
        b2 = DBLE(zvertp(i2,id))
      ENDIF

      SideLength = DSQRT((a1 - b1)**2 + (a2 - b2)**2)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE TightenGrid
      USE mod_grid_divimp
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER CalcPoint
      LOGICAL PointOnLine
      REAL*8  SideLength,VertexDisplacement

      INTEGER ik,ir,id,ik1,ir1,id1,ike,ike1,iside,v1,v2,v3,v4,fp,
     .        ik2,ir2,ik3,status,id2,id3,id4,i1,iv

      REAL*8  dist,tol,length1,length2,x(3),y(3),s,t,
     .        a1,a2,b1,b2,c1,c2,d1,d2,t1

      REAL*8, PARAMETER :: DTOL = 1.0D-06


      fp = 88

      CALL BuildMap
      CALL SetupGrid

      CALL OutputData(85,'Tightening the grid')
              
      IF (cgridopt.NE.LINEAR_GRID.AND.irsep2.NE.irsep.AND.
     .    cgridopt.NE.RIBBON_GRID) 
     .  CALL WN('TightenGrid','Not ready for double null grids '//
     .          'due to the secondary PFR')  !,*99)

      id = -1 ! TEMP

      DO ir = 2, irwall-1
        IF (ringtype(ir).EQ.PFZ) CYCLE  ! Avoiding the SOL PFZ problem for double-null grids, rather than fixing it...
        DO ik = 1, nks(ir)+1
          id1 = korpg(MIN(ik,nks(ir)),ir)
          IF (ik.LT.nks(ir)) THEN
            ik1 = ik + 1
            ir1 = ir
            id2 = korpg(ik1,ir1)
            d_rvertp(1,id2) = d_rvertp(4,id1)
            d_zvertp(1,id2) = d_zvertp(4,id1)
            d_rvertp(2,id2) = d_rvertp(3,id1)
            d_zvertp(2,id2) = d_zvertp(3,id1)
          ENDIF
          ir1 = irouts(MIN(ik,nks(ir)),ir)
          IF (idring(ir1).EQ.BOUNDARY) CYCLE
          IF (ik.EQ.nks(ir)+1) THEN
            id2 = id1
            a1 = d_rvertp(3,id1)
            a2 = d_zvertp(3,id1)
            b1 = d_rvertp(2,id1)
            b2 = d_zvertp(2,id1)
          ELSE
            id2 = korpg(MAX(1,ik-1),ir)
            a1 = d_rvertp(2,id1)
            a2 = d_zvertp(2,id1)
            b1 = d_rvertp(3,id1)
            b2 = d_zvertp(3,id1)
          ENDIF
          DO ik1 = 1, nks(ir1)+1
            id3 = korpg(MIN(ik1,nks(ir1)),ir1)
            IF (ik1.EQ.nks(ir1)+1) THEN
              id4 = id3
              c1 = d_rvertp(4,id3)
              c2 = d_zvertp(4,id3)
              d1 = d_rvertp(1,id3)
              d2 = d_zvertp(1,id3)
            ELSE
              id4 = korpg(MAX(1,ik1-1),ir1)
              c1 = d_rvertp(1,id3)
              c2 = d_zvertp(1,id3)
              d1 = d_rvertp(4,id3)
              d2 = d_zvertp(4,id3)
            ENDIF
            dist    = DSQRT((a1 - c1)**2 + (a2 - c2)**2)
            length1 = MIN(SideLength(id1,2),SideLength(id2,2)) 
            length2 = MIN(SideLength(id3,4),SideLength(id4,4)) 
            tol = MIN(1.0D-3, 0.05D0*MIN(length1,length2))
c            IF (.FALSE..AND.ik.EQ.2.and.ir.eq.24) THEN
            IF (.FALSE..AND.ik.EQ.30.and.ir.eq.119) THEN
              WRITE(0,*) '---',ik,ir,ik1,ir1
              WRITE(0,*) '---',dist,tol
              WRITE(0,*) '---',d_rvertp(2,id1),d_rvertp(1,id3)
              WRITE(0,*) '---',d_rvertp(3,id1),d_rvertp(4,id3)
              WRITE(0,*) '---',d_zvertp(2,id1),d_zvertp(1,id3)
              WRITE(0,*) '---',d_zvertp(3,id1),d_zvertp(4,id3)
            ENDIF
            IF (dist.GT.0.0D0.AND.dist.LT.tol.AND.
     .          (a1.NE.c1.OR.a2.NE.c2)) THEN
              WRITE(fp,*) 'TIGHTENING SOL:',ik,ir,ik1,ir1
              WRITE(fp,*) '              :',dist,tol      
              IF (ik1.EQ.nks(ir1)+1) THEN
                d_rvertp(4,id3) = a1
                d_zvertp(4,id3) = a2
              ELSE
                d_rvertp(1,id3) = a1
                d_zvertp(1,id3) = a2
                IF (ik1.GT.1) THEN
                  id3 = korpg(ik1-1,ir1)
                  d_rvertp(4,id3) = a1
                  d_zvertp(4,id3) = a2
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      DO ir = irsep, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        IF (ir.LT.irwall.AND.ringtype(ir).EQ.PFZ) CYCLE
c...    Low index target:
        ik1 = 1
        ir1 = ir
        ik2 = ikins(ik1,ir1)
        ir2 = irins(ik1,ir1)
        IF (.FALSE..AND.idring(ir2).EQ.BOUNDARY.AND.ir1.NE.ir2) THEN
          id1 = korpg(ik1,ir1)
          id2 = korpg(ik2,ir2)
          id3 = korpg(1  ,ir2)
          dist    = VertexDisplacement(id1,1,id3,2)
          length1 = MIN(SideLength(id1,1),SideLength(id1,4))
          length2 = MIN(SideLength(id3,1),SideLength(id3,2))
          tol = MIN(1.0D-4, 0.05D0*MIN(length1,length2))
          IF (dist.LT.tol) THEN
c          IF     (DABS(d_rvertp(1,id1)-d_rvertp(2,id3)).LT.DTOL.AND.   
c                  DABS(d_zvertp(1,id1)-d_zvertp(2,id3)).LT.DTOL) THEN  
c...        Make sure the points are exactly the same, since small errors
c           can creep in when cutting the grid:
            d_rvertp(1,id1) = d_rvertp(2,id3)
            d_zvertp(1,id1) = d_zvertp(2,id3)
          ELSEIF (irbreak.GT.0.AND.ir.GE.irbreak-1) THEN
            dist    = VertexDisplacement(id1,1,id2,2)
            length1 = MIN(SideLength(id1,1),SideLength(id1,4))
            length2 = MIN(SideLength(id2,1),SideLength(id2,2))
            tol = MIN(1.0D-4, 0.05D0*MIN(length1,length2))
            IF (dist.LT.tol) THEN
              d_rvertp(1,id1) = d_rvertp(2,id2)
              d_zvertp(1,id1) = d_zvertp(2,id2)
            ELSE
c            IF (DABS(d_rvertp(1,id1)-d_rvertp(2,id2)).GT.DTOL.OR.
c                DABS(d_zvertp(1,id1)-d_zvertp(2,id2)).GT.DTOL) THEN
c...          A step is evident between neighbouring target segments,
c             which indicates a discontinuity.  To make everything tidy,
c             split the cell so the end of the target segment that is
c             part way up the neighbouring ring sits on top of a cell 
c             vertex:
              a1 = d_rvertp(2,id2)
              a2 = d_zvertp(2,id2)
              b1 = d_rvertp(3,id2)
              b2 = d_zvertp(3,id2)
              c1 = d_rvertp(1,id1)
              c2 = d_zvertp(1,id1)
              IF (CalcPoint(a1,a2,b1,b2,c1,c2,t1).EQ.1) THEN
      WRITE(0,*) 'A1,2=',a1,a2
      WRITE(0,*) 'B1,2=',b1,b2
      WRITE(0,*) 'C1,2=',c1,c2
      WRITE(0,*) '     ',dist,tol
                WRITE(0,*) 'SPLITTING TARGET NEIGHBOUR A',t1
                CALL SplitCell(ik2,ir2,t1,status)
                CALL BuildMap
              ELSE
c               *** NOTE *** Some code added to the - B occurence below (high index target), 
c               so look there if a problem shows up here... -SL, 29/09/2010
                CALL ER('TightenGrid','Target vertex does not match '//
     .                  'boundary - A',*99)
              ENDIF
            ENDIF
          ELSE
            CALL ER('TightenGrid','Target vertex mismatch on the '//
     .              'standard grid region - A',*99)
          ENDIF
        ENDIF  
    
        ik2 = ikouts(ik1,ir1)
        ir2 = irouts(ik1,ir1)
        id2 = korpg (ik2,ir2)
        IF (.FALSE..AND.idring(ir2).NE.BOUNDARY.AND.ir1.NE.ir2) THEN
          IF (DABS(d_rvertp(2,id1)-d_rvertp(1,id2)).GT.DTOL.OR.
     .        DABS(d_zvertp(2,id1)-d_zvertp(1,id2)).GT.DTOL) THEN
            STOP 'CODE NOT READY A'
          ENDIF
        ENDIF
      ENDDO

c...  Check that target segments are properly connected to vertices on 
c     neighbouring rings, which can be a problem because ... I can't honestly
c     remember, but there must have been a problem at some point in the past
c     for me to have written this code: -SL, 29/09/2010
      DO ir = irsep, nrs
        IF (ir.LT.irwall.AND.ringtype(ir).EQ.PFZ) CYCLE
c        DO ir = irbreak-1, nrs
c        DO ir = irbreak, nrs
        IF (idring(ir).EQ.BOUNDARY.OR.(ir.EQ.irsep.AND.nopriv)) CYCLE  ! NOPRIV check necessary because of the 
c...    High index target:                                             ! extra cell added to the core rings
        ik1 = nks(ir)
        ir1 = ir
        ik2 = ikins(ik1,ir1)
        ir2 = irins(ik1,ir1)
        IF (idring(ir2).NE.BOUNDARY.AND.ir1.NE.ir2) THEN
          id1 = korpg(ik1     ,ir1)
          id2 = korpg(ik2     ,ir2)
          id3 = korpg(nks(ir2),ir2)
          dist    = VertexDisplacement(id1,4,id3,3)
          length1 = MIN(SideLength(id1,3),SideLength(id1,4))
          length2 = MIN(SideLength(id3,3),SideLength(id3,2))
          tol = MIN(1.0D-4, 0.05D0*MIN(length1,length2))
          IF (dist.LT.tol) THEN
c          IF     (DABS(d_rvertp(4,id1)-d_rvertp(3,id3)).LT.DTOL.AND.     ! NEED TO FIX THIS PROPERLY-MAR 15
c     .            DABS(d_zvertp(4,id1)-d_zvertp(3,id3)).LT.DTOL) THEN
c...        Make sure the points are exactly the same, since small errors
c           can creep in when cutting the grid: 
            d_rvertp(4,id1) = d_rvertp(3,id3)
            d_zvertp(4,id1) = d_zvertp(3,id3)
          ELSEIF (irbreak.GT.0.AND.ir.GE.irbreak-1) THEN
            dist    = VertexDisplacement(id1,4,id2,3)
            length1 = MIN(SideLength(id1,3),SideLength(id1,4))
            length2 = MIN(SideLength(id2,3),SideLength(id2,2))
            tol = MIN(1.0D-4, 0.05D0*MIN(length1,length2))
c...        Check for matching vertex on neighbour:
            IF (dist.LT.tol) THEN
              d_rvertp(4,id1) = d_rvertp(3,id2)
              d_zvertp(4,id1) = d_zvertp(3,id2)
            ELSE
c...          Brute force check to make sure there's no vertex on the
c             neighbouring ring that matches the target vertex:
              DO ik3 = 1, nks(ir2)
                id3 = korpg(ik3,ir2)
                length1 = MIN(SideLength(id1,3),SideLength(id1,4))
                length2 = MIN(SideLength(id3,3),SideLength(id3,2))
                tol = MIN(1.0D-4, 0.05D0*MIN(length1,length2))
                dist = VertexDisplacement(id1,4,id3,3)                  
                IF (dist.LT.tol) THEN
                  d_rvertp(4,id1) = d_rvertp(3,id3)
                  d_zvertp(4,id1) = d_zvertp(3,id3)
                  EXIT
                ENDIF
              ENDDO
c...          No vertex found, so try splitting the neighbouring cell
c             to create a matching vertex:
              IF (ik3.EQ.nks(ir2)+1) THEN
                a1 = d_rvertp(2,id2)
                a2 = d_zvertp(2,id2)
                b1 = d_rvertp(3,id2)
                b2 = d_zvertp(3,id2)
                c1 = d_rvertp(4,id1)
                c2 = d_zvertp(4,id1)
                IF (CalcPoint(a1,a2,b1,b2,c1,c2,t1).EQ.1) THEN
                  WRITE(0,*) 'SPLITTING TARGET NEIGHBOUR B',t1
                  CALL SplitCell(ik2,ir2,t1,status)
                  CALL BuildMap
                ELSE
                  CALL ER('TightenGrid','Target vertex does not '//
     .                    'match boundary - B',*99)
                ENDIF
              ENDIF
            ENDIF
          ELSE
            CALL ER('TightenGrid','Target vertex mismatch on the '//
     .              'standard grid region - B',*99)
          ENDIF
        ENDIF

        ir2 = ikouts(ik1,ir1)
        ir2 = irouts(ik1,ir1)
        id2 = korpg (ik2,ir2)
        IF (idring(ir2).NE.BOUNDARY.AND.ir1.NE.ir2) THEN
          IF (DABS(d_rvertp(3,id1)-d_rvertp(4,id2)).GT.DTOL.OR.
     .        DABS(d_zvertp(3,id1)-d_zvertp(4,id2)).GT.DTOL) THEN
c...        Find and split the appropriate cell:
            DO ik2 = nks(ir2)/2, nks(ir2)
              id2 = korpg(ik2,ir2)
              a1 = d_rvertp(1,id2)
              a2 = d_zvertp(1,id2)
              b1 = d_rvertp(4,id2)
              b2 = d_zvertp(4,id2)
              c1 = d_rvertp(3,id1)
              c2 = d_zvertp(3,id1)
              IF (CalcPoint(a1,a2,b1,b2,c1,c2,t1).EQ.1) THEN
                IF (t1.GT.0.0D0+DTOL.AND.t1.LT.1.0D0-DTOL) THEN
                  WRITE(0,*) 'SPLITTING TARGET NEIGHBOUR C',t1
                  CALL SplitCell(ik2,ir2,t1,status)
                  CALL BuildMap
                  EXIT
                ENDIF
              ELSE
              ENDIF
            ENDDO
          ENDIF
        ENDIF
      ENDDO

c...  This is for really messed up grids, where targets from opposite ends
c     of different rings may be touching.  This happens with grids that go
c     all the way out to the wall, creating many private flux regions:
c      DO ir1 = irsep, nrs
c        IF (idring(ir1).EQ.BOUNDARY) CYCLE
c        DO ir2 = irsep, nrs
c          IF (idring(ir2).EQ.BOUNDARY.OR.ir1.EQ.ir2) CYCLE
c          id1 = korpg(nks(ir1),ir1)
c          id2 = korpg(1       ,ir2)
c          dist    = VertexDisplacement(id1,4,id2,1)
c          length1 = MIN(SideLength(id1,3),SideLength(id1,4))
c          length2 = MIN(SideLength(id2,1),SideLength(id2,2))
c          tol = MIN(1.0D-4, 0.05D0*MIN(length1,length2))
c          IF (dist.GT.0.0D0.AND.dist.LT.tol) THEN
cc...        Make sure the points are exactly the same, since small errors
cc           can creep in when cutting the grid: 
c            WRITE(0,*) 'ZIPPING GRID',ir1,ir2
c            d_rvertp(1,id2) = d_rvertp(4,id1)
c            d_zvertp(1,id2) = d_zvertp(4,id1)
c          ENDIF
c        ENDDO
c      ENDDO

      DO ir1 = irsep, nrs
        IF (idring(ir1).EQ.BOUNDARY) CYCLE
        IF (ir1.LT.irwall.AND.ringtype(ir1).EQ.PFZ) CYCLE
        id1 = korpg(1       ,ir1)
        id2 = korpg(nks(ir1),ir1)
        DO ir2 = irsep, nrs
          IF (idring(ir2).EQ.BOUNDARY.OR.ir1.EQ.ir2) CYCLE
          DO ik = 1, nks(ir2)+1
            id = korpg(ik,ir2)
            IF (ik.EQ.nks(ir2)+1) THEN
              iv = 3
            ELSE
              iv = 2
            ENDIF
            dist    = VertexDisplacement(id1,1,id,iv)
            length1 = MIN(SideLength(id1,1),SideLength(id1,4))
            length2 = MIN(SideLength(id ,1),SideLength(id ,2))
            tol = MIN(1.0D-4, 0.05D0*MIN(length1,length2))

c            IF (ir1.EQ.128.AND.ir2.EQ.120) THEN
c              WRITE(0,*) '>>>',ir1,ir2
c              WRITE(0,*) '---',dist,tol
c              WRITE(0,*) '---',d_rvertp(1 ,id1),d_rvertp(1 ,id1)
c              WRITE(0,*) '---',d_rvertp(iv,id ),d_rvertp(iv,id )
c            ENDIF        

            IF (dist.GT.0.0D0.AND.dist.LT.tol) THEN
c...          Make sure the points are exactly the same, since small errors
c             can creep in when cutting the grid: 
              WRITE(0,*) 'MASSIVE GRID ZIPPING',ir1,ir2
              d_rvertp(1,id1) = d_rvertp(iv,id)
              d_zvertp(1,id1) = d_zvertp(iv,id)
            ENDIF
            dist    = VertexDisplacement(id2,4,id,iv)
            length1 = MIN(SideLength(id2,3),SideLength(id2,4))
            length2 = MIN(SideLength(id ,1),SideLength(id ,2))
            tol = MIN(1.0D-4, 0.05D0*MIN(length1,length2))
            IF (dist.GT.0.0D0.AND.dist.LT.tol) THEN
              WRITE(0,*) 'MASSIVE GRID ZIPPING',ir1,ir2
              d_rvertp(4,id2) = d_rvertp(iv,id)
              d_zvertp(4,id2) = d_zvertp(iv,id)
            ENDIF
            IF (ik.EQ.nks(ir)+1) THEN
              iv = 4
            ELSE
              iv = 1
            ENDIF
            dist    = VertexDisplacement(id1,2,id,iv)
            length1 = MIN(SideLength(id1,1),SideLength(id1,2))
            length2 = MIN(SideLength(id ,3),SideLength(id ,4))
            tol = MIN(1.0D-4, 0.05D0*MIN(length1,length2))
            IF (dist.GT.0.0D0.AND.dist.LT.tol) THEN
              WRITE(0,*) 'MASSIVE GRID ZIPPING',ir1,ir2
              d_rvertp(2,id1) = d_rvertp(iv,id)
              d_zvertp(2,id1) = d_zvertp(iv,id)
            ENDIF
            dist    = VertexDisplacement(id2,3,id,iv)
            length1 = MIN(SideLength(id2,3),SideLength(id2,2))
            length2 = MIN(SideLength(id ,3),SideLength(id ,4))
            tol = MIN(1.0D-4, 0.05D0*MIN(length1,length2))
            IF (dist.GT.0.0D0.AND.dist.LT.tol) THEN
              WRITE(0,*) 'MASSIVE GRID ZIPPING',ir1,ir2
              d_rvertp(3,id2) = d_rvertp(iv,id)
              d_zvertp(3,id2) = d_zvertp(iv,id)
            ENDIF


          ENDDO
        ENDDO
      ENDDO

c...  Clean up grid as necessary:
      CALL BuildMap
      DO ir = 2, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        IF (ir.LT.irwall.AND.ringtype(ir).EQ.PFZ) CYCLE
        ike = nks(ir)
        IF (ir.LT.irsep) ike = nks(ir) - 1
        DO ik = 1, ike
          id = korpg(ik,ir)
          DO iside = 1, 2
            IF (iside.EQ.1) THEN
c...          Side 14:
              v1 = 1
              v2 = 4
              v3 = 2
              v4 = 3
              ir1 = irins(ik,ir)  
            ELSE
c...          Side 23:
              v1 = 2
              v2 = 3
              v3 = 1
              v4 = 4
              ir1 = irouts(ik,ir)  
            ENDIF
            x(1) = d_rvertp(v1,id)
            y(1) = d_zvertp(v1,id)
            x(2) = d_rvertp(v2,id)
            y(2) = d_zvertp(v2,id)
            ike1 = nks(ir1) + 1
            IF (ir1.LT.irsep) ike1 = nks(ir1)
            DO ik1 = 1, ike1
              id1 = korpg(ik1,ir1)
              IF (ik1.EQ.nks(ir1)+1) THEN
                x(3) = d_rvertp(v4,id1)
                y(3) = d_zvertp(v4,id1)
              ELSE
                x(3) = d_rvertp(v3,id1)
                y(3) = d_zvertp(v3,id1)
              ENDIF

c              IF (ik.EQ.30.AND.ir.EQ.119.AND.iside.EQ.1) THEN
c                STOP 'REMOVE THIS EXCEPTION HANDLE'
c                IF (.NOT.PointOnLine(x,y,s,t,3,.TRUE.)) CYCLE
c              ELSE
                IF (.NOT.PointOnLine(x,y,s,t,3,.FALSE.)) CYCLE
c              ENDIF

              IF (DABS(s-t).GT.1.0D-8.AND.DABS(s-t).LT.1.0D0) THEN
c                IF (sloutput) THEN
c                  WRITE(0,*) '************************************'
c                  WRITE(0,*) '  PROBLEM: IK,IR,ISIDE=',ik,ir,iside
c                  WRITE(0,*) '  PROBLEM: IK,IR,ISIDE=',s,t
c                  WRITE(0,*) '  DOING NOTHING...'
c                  WRITE(0,*) '************************************'
c                  CYCLE
c                ENDIF

                IF (ik1.EQ.nks(ir1)+1) THEN
                  d_rvertp(v4,id1) = x(1) + s * (x(2) - x(1))
                  d_zvertp(v4,id1) = y(1) + s * (y(2) - y(1))
                ELSE
                  d_rvertp(v3,id1) = x(1) + s * (x(2) - x(1))
                  d_zvertp(v3,id1) = y(1) + s * (y(2) - y(1))
                ENDIF
                IF (ik1.GT.1.AND.ik1.LT.nks(ir1)+1) THEN
                  id1 = korpg(ik1-1,ir1)
                  d_rvertp(v4,id1) = x(1) + s * (x(2) - x(1))
                  d_zvertp(v4,id1) = y(1) + s * (y(2) - y(1))
                ENDIF
              ELSE
                WRITE(fp,'(A,2I6,2F15.8)') 'PASSED THE TEST',ik,ir,s,t
              ENDIF

            ENDDO
          ENDDO

        ENDDO
      ENDDO

      rvertp = SNGL(d_rvertp)
      zvertp = SNGL(d_zvertp)

      RETURN
 99   WRITE(0,*) 'IR1,2=',ir1,ir2
      WRITE(0,*) 'A1,2 =',a1,a2
      WRITE(0,*) 'B1,2 =',b1,b2
      WRITE(0,*) 'C1,2 =',c1,c2
      WRITE(0,*) 'D1,2 =',d1,d2
      rvertp = SNGL(d_rvertp)
      zvertp = SNGL(d_zvertp)
      CALL OutputGrid(85,'PROBLEM IN TIGHTENGRID')
      CALL DumpGrid('PROBLEM IN TIGHTENGRID')
      STOP
      END
c
c ========================================================================
c
c subroutine: TailorGrid
c
c Note - there are no virtual rings at this point -- should add them first?
c
cc
      SUBROUTINE TailorGrid
      USE mod_grid_divimp
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      COMMON /GRID/ iktop,irout,irin


      INTEGER iktop(MAXNRS),irout(MAXNRS),irin(MAXNRS)

      INTEGER id1,id2,id3,ir1,ir2,irs,ire,nlist,ilist(0:MAXNRS),ike,ike1


      INTEGER ii,i1,i2,i3,ir,ik,id,ik1,ik2,status,in,in1,tmpnmod,fp
      REAL*8  r1,z1,r2,z2,rvp(5,MAXNKS*MAXNRS),zvp(5,MAXNKS*MAXNRS)
c     .        a1,a2,b1,b2,c1,c2,t1
      REAL    maxz,spos

      INTEGER count
      CHARACTER buffer*512


      fp = pinout

      IF (sloutput) WRITE(fp,*) 'HERE IN TAILORGRID'


      IF (ALLOCATED(d_rvertp)) THEN
        DEALLOCATE(d_rvertp)
        DEALLOCATE(d_zvertp)
      ENDIF
      ALLOCATE(d_rvertp(5,MAXNKS*MAXNRS))
      ALLOCATE(d_zvertp(5,MAXNKS*MAXNRS))
      d_rvertp = DBLE(rvertp)
      d_zvertp = DBLE(zvertp)

c      WRITE(0,*) 'here we go...'
c      CALL SetupGrid
c      CALL BuildMap   
c      CALL TightenGrid 

c      WRITE(fp,*) 'IRBREAK:',irbreak

      irbreak = 0

c.... Delete and split rings:
      DO i1 = 1, grdnmod
 
c        WRITE(fp,*) 'GRDMOD:',i1,grdmod(i1,1)

        IF     (grdmod(i1,1).EQ.3.0) THEN
c...      Type = 1 - cut a ring
c                2 - extend the end of a ring
c                3 - delete a ring
c                4 - split a ring
 
c...      Assing cut parameters:
          irs  = NINT(grdmod(i1,4))
          ire  = NINT(grdmod(i1,5))
          DO ir = MAX(irs,ire), MIN(irs,ire), -1
            CALL DeleteRing(ir)
          ENDDO

        ELSEIF (grdmod(i1,1).EQ.4.0) THEN
c...      Assigning cut parameters:
          irs  = NINT(grdmod(i1,4))
          ire  = NINT(grdmod(i1,5))
          spos = grdmod(i1,3)
          DO ir = MAX(irs,ire), MIN(irs,ire), -1
            CALL SplitRing(ir,spos)
          ENDDO

        ELSEIF (grdmod(i1,1).EQ.5.0) THEN
c...      Duplicate a ring:
          irs  = NINT(grdmod(i1,4))
          ire  = NINT(grdmod(i1,5))
          DO ir = irs, ire
            CALL DupeRing(ir)
          ENDDO

        ELSEIF (grdmod(i1,1).EQ.9.0) THEN
c...      Extend grid:
          CALL AddOuterRing(NINT(grdmod(i1,4)),grdmod(i1,3)) 

        ELSEIF (grdmod(i1,1).EQ.10.0) THEN
c...      Squish 2 rings together:
          irs  = NINT(grdmod(i1,4))
          ire  = NINT(grdmod(i1,5))
          DO ir = irs, ire
            CALL MergeRings(ir) 
          ENDDO
c          CALL MergeRings(NINT(grdmod(i1,4))) 

        ELSEIF (grdmod(i1,1).EQ.11.0) THEN
c...      Create a new ring by expanding a ring out radially:
          CALL ExpandGrid(NINT(grdmod(i1,2)),grdmod(i1,3),
     .                    NINT(grdmod(i1,4)))

        ELSEIF (grdmod(i1,1).EQ.700.0) THEN
c...      Morph grid:
          CALL MorphGrid(NINT(grdmod(i1,2)),i1)

        ELSEIF (grdmod(i1,1).EQ.1.0.OR.grdmod(i1,1).EQ.2.0) THEN
c...      Shape the grid targets:
          rvp = d_rvertp
          zvp = d_zvertp
          IF (irbreak.EQ.0) irbreak = MAXNRS
          idring(1)      = BOUNDARY
          idring(irwall) = BOUNDARY
          idring(irtrap) = BOUNDARY
          CALL ShapeTarget(i1,rvp,zvp)
c...      Move grid vertex data into RVERTP and ZVERTP:
          d_rvertp = rvp
          d_zvertp = zvp
          rvertp = SNGL(d_rvertp)
          zvertp = SNGL(d_zvertp)

        ELSEIF (grdmod(i1,1).EQ.999.0) THEN
c...      Stop grid modifications and dump current grid state
c         for plotting in OUT:
          IF (irsep-2.NE.nrs-irtrap) THEN
            WRITE(fp,*)
            WRITE(fp,*) 'WARNING: GRID UNBALANCED',irsep-2,nrs-irtrap
            WRITE(fp,*) 
          ENDIF
          rvertp = SNGL(d_rvertp)
          zvertp = SNGL(d_zvertp)
          CALL SetupGrid
          CALL BuildMap
          CALL DumpGrid('MODIFYING RINGS AND TARGETS')

        ELSEIF (grdmod(2,1).NE.886.0.AND.
     .          (grdmod(i1,1).EQ.888.0.OR.grdmod(i1,1).GE.6.0)) THEN
c...      Moving on to polodial refinement:
          EXIT

        ELSEIF (grdmod(i1,1).EQ.885.0) THEN
c...      Hack for thesis version:
          WRITE(fp,*) 'MOVING TO POLOIDAL REFINEMENT CODE',irbreak
          tmpnmod = grdnmod
          grdnmod = 0
          EXIT

        ENDIF

      ENDDO

      CALL SetupGrid

c...  The grid is most likely a mess, fix it:
      CALL SequenceGrid

c...  Check if the core and PFZ are balanced (NOTE: this will not be an issue when
c     the EIRENE interface is converted to triangles):
      IF (sloutput.AND.pincode.LE.3.AND.irsep-2.NE.nrs-irtrap) THEN
        WRITE(fp,*)
        WRITE(fp,*) 'WARNING: GRID UNBALANCED',irsep-2,nrs-irtrap
        WRITE(fp,*) 
      ENDIF

c...  Assign IRBREAK:
      i1 = 1


      CALL FindGridBreak


      IF (sloutput) WRITE(fp,*) 'IRBREAK,NBR=',irbreak,nbr

c      CALL BuildMap

c...  Split cells along broken target so that the nearest neighbour of the
c     target cell is well defined:


      CALL TightenGrid

c      WRITE(fp,*) 'IRBREAK:',irbreak
c      CALL DumpGrid('BUILDING WALL RING')
c      STOP


c...  Reset IRSEP2 for connected double-null grids, which may have happened
c     in the call to SequenceGrid:
      IF (connected) irsep2 = irouts(ikins(ikti-1,irsep)+1,irsep-1)


      IF (sloutput) WRITE(fp,*) 'DONE IN BUILDMAP'


      IF (grdmod(1,1).EQ.887.0) grdnmod = tmpnmod

      rvertp = SNGL(d_rvertp)
      zvertp = SNGL(d_zvertp)

      DO i1 = 1, grdnmod

        IF     (grdmod(i1,1).EQ.6.0) THEN
          IF (grdmod(i1,2).EQ.-1.0) THEN
c...        Change minium poloidal cell size:
            grd_minpl = grdmod(i1,3)
          ELSE
c...        Poloidal refinement:
            IF (.FALSE..AND.grdmod(i1,2).EQ.5.0) THEN
              CALL OUTPUTDATA(85,'sdfsd')
              STOP 'sdfsd'
            ENDIF

            CALL SetupGrid

            irs = NINT(grdmod(i1,4))
            ire = NINT(grdmod(i1,5))
            IF (irs.EQ.-99) irs = irsep
            IF (ire.EQ.-99) ire = nrs
            DO ir = irs, ire
             IF (ir.NE.0.AND.idring(MAX(1,ir)).EQ.BOUNDARY) CYCLE
             CALL PoloidalRefinement(ir,NINT(grdmod(i1,2)),grdmod(i1,3))
c             WRITE(0,*) 'POLOIDAL REF:',ir,nks(ir)
            ENDDO

c            WRITE(0,*) 'here we go...'
c            CALL BuildMap
c            CALL TightenGrid
          ENDIF

c          rvertp = SNGL(d_rvertp)
c          zvertp = SNGL(d_zvertp)

        ELSEIF (grdmod(i1,1).EQ.7.0) THEN
c...      Dump grid data to an external file:
          OPEN (UNIT=98,FILE='cell.dat',ACCESS='SEQUENTIAL',
     .          STATUS='NEW')      
          IF (.FALSE.) THEN
            WRITE(98,*) 'SHOT: 119919   TIME: 3000'
            WRITE(98,'(2A6,2A10)') 'IK','IR','R (m)','Z (m)'
            DO ir = 1, nrs
              IF (idring(ir).EQ.BOUNDARY) CYCLE
              DO ik = 1, nks(ir)
                WRITE(98,'(2I6,2F10.6)') ik,ir,rs(ik,ir),zs(ik,ir)
              ENDDO
            ENDDO
          ELSEIF (.TRUE.) THEN
            WRITE(98,'(A)') '* Grid geometry file from DIVIMP'
            WRITE(98,'(A,1X,2A6,6A12)') 
     .        '*','IK','IR','Rcen (m)','Zcen (m)','Rmid1 (m)',
     .        'Zmid1 (m)','R1 (m)','Z1 (m)'
            WRITE(98,*) 3  ! Format code
            DO ir = 2, nrs
              IF (idring(ir).EQ.BOUNDARY) CYCLE
c...          Points in core near the x-point (more accurate due to flux expansion):
              DO ik = 1, nks(ir)
                id = korpg(ik,ir)
                WRITE(98,'(2X,2I6,6F12.7)') ik,ir,
     .            rs(ik,ir),zs(ik,ir),
     .            0.5*(rvertp(1,id)+rvertp(2,id)),
     .            0.5*(zvertp(1,id)+zvertp(2,id)),
     .            rvertp(1,id),
     .            zvertp(1,id)
              ENDDO
            ENDDO
          ELSEIF (.FALSE.) THEN
            WRITE(98,*) 'SHOT: 119919   TIME: 3000'
            WRITE(98,'(2A6,2A10)') 'IK','IR','R (m)','Z (m)'
            DO ir = 1, irsep-1
              IF (idring(ir).EQ.BOUNDARY) CYCLE
c...          Points in core near the x-point (more accurate due to flux expansion):
              ik = 1
              id = korpg(ik,ir)
              WRITE(98,'(2I6,2F10.6)') ik,ir,
     .          0.5*(rvertp(1,id)+rvertp(2,id)),
     .          0.5*(zvertp(1,id)+zvertp(2,id))
            ENDDO
            DO ir = 1, nrs
              IF (idring(ir).EQ.BOUNDARY) CYCLE
c...          Points at the centers of inner targets:
              IF (ir.GE.irsep) THEN
                ik = 1
                id = korpg(ik,ir)
                WRITE(98,'(2I6,2F10.6)') ik,ir,
     .            0.5*(rvertp(1,id)+rvertp(2,id)),
     .            0.5*(zvertp(1,id)+zvertp(2,id))
              ENDIF
            ENDDO
            DO ir = 1, nrs
              IF (idring(ir).EQ.BOUNDARY) CYCLE
c...          Points near the outer midplane:
              IF (ir.LE.irwall) THEN
                DO ik = 1, nks(ir)
                  IF (rs(ik,ir).GT.r0.AND.zs(ik  ,ir).GT.z0.AND.      ! I think RS/ZS are properly defined here...
     .                                    zs(ik+1,ir).LT.z0) THEN
                    id = korpg(ik,ir)
                    WRITE(98,'(2I6,2F10.6)') ik,ir,
     .                0.5*(rvertp(3,id)+rvertp(4,id)),
     .                0.5*(zvertp(3,id)+zvertp(4,id))
                  ENDIF
                ENDDO
              ENDIF
            ENDDO
            DO ir = 1, nrs
              IF (idring(ir).EQ.BOUNDARY) CYCLE
c...          Points at the centers of outer targets:
              IF (ir.GE.irsep) THEN
                ik = nks(ir)
                id = korpg(ik,ir)
                WRITE(98,'(2I6,2F10.6)') ik,ir,
     .            0.5*(rvertp(3,id)+rvertp(4,id)),
     .            0.5*(zvertp(3,id)+zvertp(4,id))
              ENDIF
            ENDDO
          ENDIF
          CLOSE (98)

        ELSEIF (grdmod(i1,1).EQ.8.0) THEN

          IF     (grdmod(i1,2).EQ. 0.0.OR.
     .            grdmod(i1,2).EQ. 1.0.OR.
     .            grdmod(i1,2).EQ. 2.0.OR.
     .            grdmod(i1,2).EQ. 3.0.OR.
     .            grdmod(i1,2).EQ.10.0) THEN
c...        Load Bratio and PSIn data from external file:
            CALL LoadGridData(NINT(grdmod(i1,2)))

          ELSEIF (grdmod(i1,2).EQ.99.0) THEN
c *HACK*
            REWIND(4)
            DO WHILE (buffer(1:16).NE.'PSI-DOUBLE-NULLd') 
              READ(4,'(A512)') buffer
            ENDDO
            READ(buffer(17:),*) i2
            psitarg = 0.0
            DO i3 = 1, i2
              READ(4,*) ir,psitarg(ir,2),psitarg(ir,1)
            ENDDO
c...        Quick check:
            DO ir = 2, nrs
              IF (idring(ir).EQ.BOUNDARY) CYCLE
              IF (psitarg(ir,2).EQ.0.0.OR.psitarg(ir,1).EQ.0.0)
     .          CALL ER('TailorGrid','Bad PSITARG hack',*99)
            ENDDO
            IF (sloutput) THEN
              WRITE(fp,*)
              WRITE(fp,*) '----------------------'
              WRITE(fp,*) 'PSITARG -- SUPER-HACK!'
              WRITE(fp,*) '----------------------'
              WRITE(fp,*)
            ENDIF 
          ELSE
            CALL ER('TailorGrid','Invalid .dat grid option',*99)
          ENDIF

        ELSEIF (grdmod(i1,1).EQ.11.0) THEN
c...      Another chance to delete rings, which are ordered now:
          irs  = NINT(grdmod(i1,4))
          ire  = NINT(grdmod(i1,5))
          DO ir = MAX(irs,ire), MIN(irs,ire), -1
            CALL DeleteRing(ir)
          ENDDO

        ELSEIF (grdmod(i1,1).EQ.998.0) THEN
c...      Pretend the grid wasn't modified:
          irbreak = 0
          nbr = 0

          CALL BuildMap
          CALL SetupGrid

          grdnmod = 0

          cutpt1   = ikto + 1
          cutpt2   = ikti + 1
          maxkpts  = nks(irsep) + 2  

c          CALL DumpGrid('TRICY POLOIDAL DISTRIBTION')
c          RETURN
          EXIT

        ELSEIF (grdmod(i1,1).EQ.999.0) THEN
c...      Stop grid modifications and dump current grid state
c         for plotting in OUT:
          CALL BuildMap
          CALL SetupGrid
          CALL DumpGrid('MODIFYING POLOIDAL DISTRIBUTION')
        ENDIF

      ENDDO

c...  Hack: Turn off new generalized grid code:
      IF (grdmod(1,1).EQ.887.0) grdnmod = 0

c...  Calculate some derived grid quantities:
      CALL SetupGrid
      CALL TightenGrid

      IF (nopriv) THEN
        rxp = rvertp(1,korpg(1,irsep))
        zxp = zvertp(1,korpg(1,irsep))
      ELSE
        IF (sloutput) WRITE(fp,*) 'ADJUSTING X-POINT'
        rxp = rvertp(4,korpg(ikto,irsep))
        zxp = zvertp(4,korpg(ikto,irsep))
      ENDIF

c      DO ik = 1, 10
c        id = korpg(ik,irsep)
c        WRITE(0,*) 'CHECK: ',rvertp(1,id),d_rvertp(1,id)
c        WRITE(0,*) '     : ',zvertp(1,id),d_zvertp(1,id)
c        WRITE(0,*) '     : ',rvertp(4,id),d_rvertp(4,id)
c        WRITE(0,*) '     : ',zvertp(4,id),d_zvertp(4,id)
c      ENDDO

c      DEALLOCATE(d_rvertp)
c      DEALLOCATE(d_zvertp)

      RETURN
99    CONTINUE
c      WRITE(0,*) 'IR=',ir,t1,id1
c      WRITE(0,*) 'IK1,IR1,D1=',ik1,ir1,id1
c      WRITE(0,*) 'IR=',c1,c2
c      WRITE(0,*) 'IK2,IR1,D2=',ik2,ir2,id2
c      WRITE(0,*) 'A1,2=',a1,a2
c      WRITE(0,*) 'B1,2=',b1,b2
c      WRITE(0,*) 'C1,2=',c1,c2
      CALL OutputGrid(85,'PROBLEM IN TAILORGRID')
      CALL DumpGrid('PROBLEM IN TAILORGRID')
      STOP
      END
c
c ========================================================================
c ========================================================================
c
c 
c
c ========================================================================
c
      SUBROUTINE FindKnot(nknot,knot,NUMZONE,izone,condition,
     .                    index1,index2)
      IMPLICIT none
      INCLUDE 'params'  ! for SLOUTPUT
c
c     jdemod
c
c     Condition=1 finds two different knots which have the same 
c                 R,Z values for knot number 1. This condition
c                 should only occur at an Xpoint for two cells in the 
c                 inner/outer SOL
c
c     NOTE: The cells sharing vertex 1 or 4 will be in the inner and outer SOL
c           The cells sharing vertex 2 or 3 will be in the core and PFZ
c
c           There seems to be a bug in the code when finding the innermost
c           core ring - using vertex 1 finds a cell in the main SOL adjacent
c           to the core. Stepping inward will usually work if the cell 
c           next to the core is chosen - however, if the separatrix is unusual
c           it could be that the Z coordinate of the cell adjacent to the PFZ
c           will be closer to the center of the plasma and thus the wrong cell
c           will be chosen to find the center of the grid. In order to fix this, 
c           the test vertex for the XPoint should be set to 2 or 3 so that 
c           cells inside and outside the core are chosen - more clearly above 
c           and below the Xpoints. I'll have to check and see if this has an
c           impact on the calculations for double null grids. 
c
c
c     Condition=2 finds knot 1=2 and knot 4=3 for test cell vs. other cells
c
c     Condition=3 finds knot 3=2 and knot 4=1 for test cell vs. other cells
c
c     Condition=4 finds knot 2=1 and knot 3=4 for test cell vs. other cells
c
c     Condition=5 finds knot 1=4 and knot 2=3 for test cell vs. other cells
c
c     Condition=6 finds two different knots which have the same 
c                 R,Z values for knot number 2. This condition
c                 should only occur at an Xpoint for two cells in the core/pfz
c
c

      TYPE type_cell
        INTEGER :: index,ik,ir,nv,rzone,zzone,xpt,map
        REAL    :: rcen,zcen,bratio,rv(4),zv(4)
      ENDTYPE type_cell

      INTEGER nknot,index1,index2,NUMZONE,izone(NUMZONE+1,NUMZONE),
     .        condition
      TYPE(type_cell) :: knot(0:nknot)      

      REAL*8, PARAMETER :: DTOL=1.0D-06

      INTEGER i1,i2,r1,z1

c
c     jdemod - Logical for debugging
c
      logical output
c
      output = .false.

      i1 = index1

      index2 = -1

      if (condition.ne.1.and.output) then 

         write(6,'(a,100i8)') 'FINDKNOT1:',
     >      nknot,NUMZONE,condition,index1,index2
         write(6,'(a,100i8)') 'FINDKNOT:',
     >      izone
c
         write(0,'(a,100i8)') 'FINDKNOT1:',
     >      nknot,NUMZONE,condition,index1,index2
c      write(0,'(a,100i8)') 'FINDKNOT:',
c     >      izone

      endif

      DO z1 = knot(i1)%zzone-1, knot(i1)%zzone+1
        IF (z1.LT.1.OR.z1.GT.NUMZONE) CYCLE
        DO r1 = knot(i1)%rzone-1, knot(i1)%rzone+1
          IF (r1.LT.1.OR.r1.GT.NUMZONE) CYCLE  

c          DO i2 = 1, nknot
          DO i2 = izone(r1,z1), izone(r1+1,z1)-1
             
            if (output) then 
               write(6,'(a,15i8)') 'IZONE:',i1,z1,r1,i2,izone(r1,z1),
     >                izone(r1+1,z1)-1
            endif

            IF (i1.EQ.i2) CYCLE
c...
            IF     (condition.EQ.1.AND.
     .              ABS(knot(i1)%rv(1)-knot(i2)%rv(1)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(1)-knot(i2)%zv(1)).LT.DTOL) THEN

              IF (sloutput) WRITE(0,*) 'XPOINT SOL:',i1,i2
              if (output) then
                 WRITE(0,*) 'XPOINT SOL:',i1,i2
                 WRITE(6,*) 'XPOINT SOL:',i1,i2
              endif
                 
              index2 = i2
              RETURN
c
c     jdemod - Add vertex 2 for Xpoint detection for core rings so that the cells found
c              will be clearly above and below the Xpoint
c

            ELSEIF     (condition.EQ.6.AND.
     .              ABS(knot(i1)%rv(2)-knot(i2)%rv(2)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(2)-knot(i2)%zv(2)).LT.DTOL) THEN

              IF (sloutput) WRITE(0,*) 'XPOINT CORE/PFZ:',i1,i2
              if (output) then
                 WRITE(0,*) 'XPOINT CORE/PFZ:',i1,i2
                 WRITE(6,*) 'XPOINT CORE/PFZ:',i1,i2
              endif
                 
              index2 = i2
              RETURN

            ELSEIF (condition.EQ.2.AND.
     .              ABS(knot(i1)%rv(1)-knot(i2)%rv(2)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(1)-knot(i2)%zv(2)).LT.DTOL.AND.
     .              ABS(knot(i1)%rv(4)-knot(i2)%rv(3)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(4)-knot(i2)%zv(3)).LT.DTOL) THEN

              if (output) then 
                WRITE(0,*) 'SIDE INWARD 41:',i1,i2
                WRITE(6,*) 'SIDE INWARD 41:',i1,i2
              endif

              index2 = i2
              RETURN

            ELSEIF (condition.EQ.3.AND.
     .              ABS(knot(i1)%rv(3)-knot(i2)%rv(2)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(3)-knot(i2)%zv(2)).LT.DTOL.AND.
     .              ABS(knot(i1)%rv(4)-knot(i2)%rv(1)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(4)-knot(i2)%zv(1)).LT.DTOL) THEN

              if (output) then 
                WRITE(0,*) 'SIDE UP 34:',i1,i2
                WRITE(6,*) 'SIDE UP 34:',i1,i2
              endif

              index2 = i2
              RETURN

            ELSEIF (condition.EQ.4.AND.
     .              ABS(knot(i1)%rv(2)-knot(i2)%rv(1)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(2)-knot(i2)%zv(1)).LT.DTOL.AND.
     .              ABS(knot(i1)%rv(3)-knot(i2)%rv(4)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(3)-knot(i2)%zv(4)).LT.DTOL) THEN


               if (output) then 
                 WRITE(0,*) 'SIDE OUTWARD 23:',i1,i2
                 WRITE(6,*) 'SIDE OUTWARD 23:',i1,i2
               endif

              index2 = i2
              RETURN

            ELSEIF (condition.EQ.5.AND.
     .              ABS(knot(i1)%rv(1)-knot(i2)%rv(4)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(1)-knot(i2)%zv(4)).LT.DTOL.AND.
     .              ABS(knot(i1)%rv(2)-knot(i2)%rv(3)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(2)-knot(i2)%zv(3)).LT.DTOL) THEN
c...          Matching sides 12 and 34:

               if (output) then
                 WRITE(0,*) 'SIDE DOWN 12:',i1,i2
                 WRITE(6,*) 'SIDE DOWN 12:',i1,i2
               endif

              index2 = i2
              RETURN



            ENDIF

          ENDDO

        ENDDO
      ENDDO

c
c     If the code reaches here it has failed to find a match in the zoning system - try scanning all cells
c
      if (condition.ne.1) then 

      if (output) then 
         write(0,*) 
     >     'FAILED TO FIND KNOT IN ZONED SYSTEM: CHECKING ALL CELLS'
         write(6,*) 
     >     'FAILED TO FIND KNOT IN ZONED SYSTEM: CHECKING ALL CELLS'
      endif
c
c      DO z1 = knot(i1)%zzone-1, knot(i1)%zzone+1
c        IF (z1.LT.1.OR.z1.GT.NUMZONE) CYCLE
c        DO r1 = knot(i1)%rzone-1, knot(i1)%rzone+1
c          IF (r1.LT.1.OR.r1.GT.NUMZONE) CYCLE  

        DO i2 = 1, nknot
c          DO i2 = izone(r1,z1), izone(r1+1,z1)-1
             
            write(6,'(a,15i8)') 'IZONE:',i1,z1,r1,i2,izone(r1,z1),
     >                izone(r1+1,z1)-1

            IF (i1.EQ.i2) CYCLE
c...
            IF     (condition.EQ.1.AND.
     .              ABS(knot(i1)%rv(1)-knot(i2)%rv(1)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(1)-knot(i2)%zv(1)).LT.DTOL) THEN

              IF (sloutput) WRITE(0,*) 'XPOINT:',i1,i2
              if (output) then
                 write(0,*) 'ZONING ERROR: CELL FOUND:',i2,
     >                        knot(i2)%rzone,knot(i2)%zzone
                 write(6,*) 'ZONING ERROR: CELL FOUND:',i2,
     >                        knot(i2)%rzone,knot(i2)%zzone
                 WRITE(0,*) 'XPOINT:',i1,i2
                 WRITE(6,*) 'XPOINT:',i1,i2
              endif
                 
              index2 = i2
              RETURN

            ELSEIF (condition.EQ.2.AND.
     .              ABS(knot(i1)%rv(1)-knot(i2)%rv(2)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(1)-knot(i2)%zv(2)).LT.DTOL.AND.
     .              ABS(knot(i1)%rv(4)-knot(i2)%rv(3)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(4)-knot(i2)%zv(3)).LT.DTOL) THEN

              if (output) then 
                write(0,*) 'ZONING ERROR: CELL FOUND:',i2,
     >                        knot(i2)%rzone,knot(i2)%zzone
                write(6,*) 'ZONING ERROR: CELL FOUND:',i2,
     >                        knot(i2)%rzone,knot(i2)%zzone
                WRITE(0,*) 'SIDE INWARD 41:',i1,i2
                WRITE(6,*) 'SIDE INWARD 41:',i1,i2
              endif

              index2 = i2
              RETURN

            ELSEIF (condition.EQ.3.AND.
     .              ABS(knot(i1)%rv(3)-knot(i2)%rv(2)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(3)-knot(i2)%zv(2)).LT.DTOL.AND.
     .              ABS(knot(i1)%rv(4)-knot(i2)%rv(1)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(4)-knot(i2)%zv(1)).LT.DTOL) THEN

              if (output) then 
                write(0,*) 'ZONING ERROR: CELL FOUND:',i2,
     >                        knot(i2)%rzone,knot(i2)%zzone
                write(6,*) 'ZONING ERROR: CELL FOUND:',i2,
     >                        knot(i2)%rzone,knot(i2)%zzone
                WRITE(0,*) 'SIDE UP 34:',i1,i2
                WRITE(6,*) 'SIDE UP 34:',i1,i2
              endif

              index2 = i2
              RETURN

            ELSEIF (condition.EQ.4.AND.
     .              ABS(knot(i1)%rv(2)-knot(i2)%rv(1)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(2)-knot(i2)%zv(1)).LT.DTOL.AND.
     .              ABS(knot(i1)%rv(3)-knot(i2)%rv(4)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(3)-knot(i2)%zv(4)).LT.DTOL) THEN


               if (output) then 
                 write(0,*) 'ZONING ERROR: CELL FOUND:',i2,
     >                        knot(i2)%rzone,knot(i2)%zzone
                 write(6,*) 'ZONING ERROR: CELL FOUND:',2,
     >                        knot(i2)%rzone,knot(i2)%zzone
                 WRITE(0,*) 'SIDE OUTWARD 23:',i1,i2
                 WRITE(6,*) 'SIDE OUTWARD 23:',i1,i2
               endif

              index2 = i2
              RETURN

            ELSEIF (condition.EQ.5.AND.
     .              ABS(knot(i1)%rv(1)-knot(i2)%rv(4)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(1)-knot(i2)%zv(4)).LT.DTOL.AND.
     .              ABS(knot(i1)%rv(2)-knot(i2)%rv(3)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(2)-knot(i2)%zv(3)).LT.DTOL) THEN
c...          Matching sides 12 and 34:

               if (output) then
                 write(0,*) 'ZONING ERROR: CELL FOUND:',i2,
     >                        knot(i2)%rzone,knot(i2)%zzone
                 write(6,*) 'ZONING ERROR: CELL FOUND:',i2,
     >                        knot(i2)%rzone,knot(i2)%zzone
                 WRITE(0,*) 'SIDE DOWN 12:',i1,i2
                 WRITE(6,*) 'SIDE DOWN 12:',i1,i2
               endif

              index2 = i2
              RETURN



            ENDIF

          ENDDO

c        ENDDO
c      ENDDO
      endif



      RETURN
 99   STOP
      END

c
c ========================================================================
c
      SUBROUTINE MoveKnot(knot1,knot)
      IMPLICIT none

      TYPE type_cell
        INTEGER :: index,ik,ir,nv,rzone,zzone,xpt,map
        REAL    :: rcen,zcen,bratio,rv(4),zv(4)
      ENDTYPE type_cell

      TYPE(type_cell) :: knot1,knot      

      INTEGER i1

      knot%index  = knot1%index
      knot%ik     = knot1%ik
      knot%ir     = knot1%ir
      knot%rzone  = knot1%rzone
      knot%zzone  = knot1%zzone
      knot%xpt    = knot1%xpt
      knot%rcen   = knot1%rcen
      knot%zcen   = knot1%zcen
      knot%bratio = knot1%bratio
      DO i1 = 1, 4
        knot%rv(i1) = knot1%rv(i1)
        knot%zv(i1) = knot1%zv(i1)
      ENDDO

      RETURN
 99   STOP
      END
c
c ========================================================================
c
      SUBROUTINE ReadGeneralisedGrid(gridunit,ik,ir,
     .                               rshift,zshift,indexiradj)
      USE mod_grid_divimp
      IMPLICIT none

      INTEGER gridunit,ik,ir,indexiradj
      REAL    rshift,zshift

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'
      INCLUDE 'pindata'

c..TMP
      CHARACTER title*174,desc*1024,job*72,equil*60
      REAL      facta(-1:MAXIZS),factb(-1:MAXIZS)

      INTEGER, PARAMETER :: NUMZONE = 5
      REAL*8,  PARAMETER :: TOL = 1.0D-06

      INTEGER   nknot,i1,i2,z1,r1,kind,nxpt,ixpt(0:2),cxpt(0:2),i3,
     .          izone(NUMZONE+1,NUMZONE),newi1,icore(0:2),id,tmpnks,
     .          ikmax,irmax,ir1,istart,
     .          numpsi,ikpsi(MAXNRS),irpsi(MAXNRS)
      LOGICAL   cont,deleteknot,output,swap
      REAL      vrmin,vzmin,vrmax,vzmax,rspan,zspan,area,valpsi(MAXNRS)
      REAL*8    rvdp(4),zvdp(4),areadp
      CHARACTER buffer*1000

      INTEGER, ALLOCATABLE :: imap(:,:)

      TYPE type_cell
        INTEGER :: index,ik,ir,nv,rzone,zzone,xpt,map
        REAL    :: rcen,zcen,bratio,rv(4),zv(4)
      ENDTYPE type_cell
 
      TYPE(type_cell),ALLOCATABLE :: knot(:)

c
c     jdemod - flag for turning on and off debugging output
c
      output = .TRUE.

      ALLOCATE(knot(0:MAXNKS*MAXNRS))
      ALLOCATE(imap(MAXNKS,0:MAXNRS))

c...  Find the start of the cell/knot information in the grid file:
      WRITE(buffer,'(1000X)')
      DO WHILE (buffer(4:8).NE.'=====')
        READ(gridunit,'(A10)',END=98) buffer
      ENDDO

c...  Read the knot data:
      nknot = 0
      DO WHILE(nknot.EQ.0.OR.buffer(4:10).EQ.'Element')
        nknot = nknot + 1
        READ(gridunit,80,END=97) knot(nknot)%index,
     .                           knot(nknot)%ik   ,knot(nknot)%ir,
     .                           knot(nknot)%rv(2),knot(nknot)%zv(2),
     .                           knot(nknot)%rv(3),knot(nknot)%zv(3)
        READ(gridunit,81,END=97) knot(nknot)%bratio,
     .                           knot(nknot)%rcen ,knot(nknot)%zcen
        READ(gridunit,82,END=97) knot(nknot)%rv(1),knot(nknot)%zv(1),
     .                           knot(nknot)%rv(4),knot(nknot)%zv(4)
        knot(nknot)%nv = 4
c...    Dividing line:       
        READ(gridunit,*)
        READ(gridunit,'(A10)',END=20) buffer
        BACKSPACE(gridunit)
      ENDDO
 80   FORMAT(10X,I5,4X,I3,1x,i3,4x,e17.10e2,1x,e17.10e2,8x,e17.10e2,1x,
     .       E17.10E2)
 81   FORMAT(18x,e17.10e2,14x,e17.10e2,1x,e17.10e2)
 82   FORMAT(30x,e17.10e2,1x,e17.10e2,8x,e17.10e2,1x,e17.10e2)
c...  End of file continuation:
 20   CONTINUE



c...  Delete zero volume cells:

c...  Strip those boundary cells:
c     jdemod
      IF (output) then
         write(0,*) 'CELLS READ IN:', nknot,
     >            knot(nknot)%ir,knot(nknot)%ik
         write(6,*) 'CELLS READ IN:', nknot,
     >            knot(nknot)%ir,knot(nknot)%ik
         WRITE(0,*) 'STRIPPING...'
         WRITE(6,*) 'STRIPPING...'
      endif
c
      ikmax = 0
      irmax = 0
      DO i1 = 1, nknot
        IF (knot(i1)%ik.GT.ikmax) ikmax = knot(i1)%ik
        IF (knot(i1)%ir.GT.irmax) irmax = knot(i1)%ir
      ENDDO

c     jdemod
      IF (output) then
         write(0,*) 'IKMAX,IRMAX:', ikmax,irmax
         write(6,*) 'IKMAX,IRMAX:', ikmax,irmax
      endif


c...  (This is not the most efficient way of doing things -- certaintly it would
c      be better just to avoid storing the cells as the grid file is read in --    <---WHAT?
c      but it is general, which is the goal here, and makes minimal assumptions 
c      about the structure of the grid file):
c
c     jdemod - the following code loops through the grid and removes all cells for which the 
c              coordinates of vertex 3=2 and 4=1 - these are the parallel to the field line 
c              vertices - this will remove cells at the targets but not any boundary rings
c
      i1 = 1
      DO WHILE(i1.LE.nknot)
        IF     (.FALSE..AND.
     .          (knot(i1)%ik.EQ.0.OR.knot(i1)%ik.EQ.ikmax.OR.     ! virtual cells on end of rings
     .           knot(i1)%ir.EQ.0.OR.knot(i1)%ir.EQ.irmax)) THEN  ! virtual rings
c...      Remove boundary knots:
          deleteknot = .TRUE.
        ELSEIF (ABS(knot(i1)%rv(3)-knot(i1)%rv(2)).LT.TOL.AND.
     .          ABS(knot(i1)%zv(3)-knot(i1)%zv(2)).LT.TOL.AND.
     .          ABS(knot(i1)%rv(4)-knot(i1)%rv(1)).LT.TOL.AND.
     .          ABS(knot(i1)%zv(4)-knot(i1)%zv(1)).LT.TOL) THEN
c...      Also get rid of zero volume cells, which can be present in UEDGE
c         double null grids.  The above condition is the best identifier
c         for these (for grids generated with UEDGE anyway):
          deleteknot = .TRUE.
        ELSE
c...      Cell to be kept, advance index:
          deleteknot = .FALSE.
          i1 = i1 + 1
        ENDIF
c...    Delete the knot:
        IF (deleteknot) THEN
          IF (i1.LT.nknot) THEN
            DO i2 = i1, nknot
              CALL MoveKnot(knot(i2+1),knot(i2))
            ENDDO
          ENDIF
          nknot = nknot - 1
c          jdemod
          IF (output) then
             WRITE(0,*) 'GONE',i1
             WRITE(6,*) 'GONE',i1
          endif

        ENDIF
      ENDDO
c
c     jdemod
      IF (output) then 
         WRITE(0,*) 'DONE',nknot
         WRITE(6,*) 'DONE',nknot
      endif


c...  Search the grid for remaining virtual cells (typically zero-volume):
c      DO i1 = 1, nknot
c        DO i2 = 1, knot(i1)%nv
c          rvdp(i2) = DBLE(knot(i1)%rv(i2))
c          zvdp(i2) = DBLE(knot(i1)%zv(i2))
c        ENDDO
c        areadp = 0.0
c        DO i2 = 1, knot(i1)%nv
c          i3 = i2 + 1
c          IF (i2.EQ.knot(i1)%nv) i3 = 1
c          areadp = areadp + (rvdp(i3) * zvdp(i2) -
c     .                       rvdp(i2) * zvdp(i3))
c        ENDDO
c        area = 0.5 * SNGL(DABS(areadp))
c
c         cont = .FALSE.
c          IF (DABS(rvdp(3)-rvdp(2)).LT.TOL.AND.
c     .        DABS(zvdp(3)-zvdp(2)).LT.TOL.AND.
c     .        DABS(rvdp(4)-rvdp(1)).LT.TOL.AND.
c     .        DABS(zvdp(4)-zvdp(1)).LT.TOL) cont = .TRUE.
c
c
c        IF (area.LT.1.0E-08) 
c     .    WRITE(0,*) 'I!,AREA:',knot(i1)%index,area,cont
c      ENDDO
c
c      STOP 'sdfsd'

c...  R,Z shifts:

c...  Assign knot sector (for improved efficiency in the search routines):
      vrmin =  HI
      vrmax = -HI
      vzmin =  HI
      vzmax = -HI
      DO i1 = 1, nknot      
        DO i2 = 1, knot(i1)%nv
          IF (knot(i1)%rv(i2).LT.vrmin) vrmin = knot(i1)%rv(i2)
          IF (knot(i1)%rv(i2).GT.vrmax) vrmax = knot(i1)%rv(i2)
          IF (knot(i1)%zv(i2).LT.vzmin) vzmin = knot(i1)%zv(i2)
          IF (knot(i1)%zv(i2).GT.vzmax) vzmax = knot(i1)%zv(i2)
        ENDDO
      ENDDO
      vrmin = vrmin - 0.001
      vrmax = vrmax + 0.001
      vzmin = vzmin - 0.001
      vzmax = vzmax + 0.001
c
c     jdemod
      IF (output) then
         WRITE(0,*) 'MIN,MAX:',vrmin,vrmax,vzmin,vzmax
         WRITE(6,*) 'MIN,MAX:',vrmin,vrmax,vzmin,vzmax
      endif
c      IF (output) WRITE(0,*) 'MIN,MAX:',vrmin,vrmax,vzmin,vzmax
c...  Assign knot to search zone:
      rspan = (vrmax - vrmin) / REAL(NUMZONE)
      zspan = (vzmax - vzmin) / REAL(NUMZONE)
      DO i1 = 1, nknot
        knot(i1)%rzone = INT( (knot(i1)%rcen - vrmin) / rspan ) + 1
        knot(i1)%zzone = INT( (knot(i1)%zcen - vzmin) / zspan ) + 1
c        WRITE(0,*) 'SPAN:',i1,knot(i1)%rzone,knot(i1)%zzone
      ENDDO

c...  Index cells by zone, work now for saved time later:
      kind = 1
      DO z1 = 1, NUMZONE
        DO r1 = 1, NUMZONE
          izone(r1,z1) = kind
          DO i1 = kind, nknot
            IF     (knot(i1)%rzone.EQ.r1.AND.knot(i1)%zzone.EQ.z1) THEN
              IF (i1.EQ.kind) THEN
c...            Do nothing:
              ELSE
c...            Swap knots:
                CALL MoveKnot(knot(kind),knot(0)   )
                CALL MoveKnot(knot(i1)  ,knot(kind))
                CALL MoveKnot(knot(0)   ,knot(i1)  )
c                WRITE(0,*) 'ZONING SWAP:',kind,i1
              ENDIF
              kind = kind + 1
            ENDIF
          ENDDO
        ENDDO
      ENDDO
c...  
      DO i1 = 1, NUMZONE-1
        izone(NUMZONE+1,i1) = izone(1,i1+1)
      ENDDO
      izone(NUMZONE+1,NUMZONE) = nknot + 1

c      DO i1 = 1, nknot
c        WRITE(0,*) 'ZONED:',i1,knot(i1)%rzone,knot(i1)%zzone
c      ENDDO
c      DO z1 = 1, NUMZONE
c        DO r1 = 1, NUMZONE
c          WRITE(0,*) 'IZONE:',izone(r1,z1) 
c        ENDDO
c      ENDDO

c...  Search for x-point(s):
c
c     jdemod - it appears that the X-point finding algorithm used is the following:
c            - the only cells on the grid which can have an identical vertex - both index and value 
c              and not be the same cell must occur at the Xpoint - the shared vertex in this 
c              case IS the Xpoint. This is the condition tested for when FindKnot is called with 
c              a 1. Search efficiency has been enhanced by using the zones set up above. 
c            - Zero volume cells would be an issue with this algorithm - code above this has 
c              apparently removed zero volume cells where vertices 3=2 and 1=4 - however, it would 
c              appear that boundary rings around the plasma and in the core/PFZ have not been 
c              removed.
c
      nxpt = 0
      DO i1 = 1, nknot
        IF (knot(i1)%xpt.NE.0) CYCLE

        IF (nxpt.EQ.2) THEN
          WRITE(0,*)
          WRITE(0,*) '--------------------------------------------'
          WRITE(0,*) '-   MORE THAN 2 XPTS FOUND, IGNORING...    -'
          WRITE(0,*) '--------------------------------------------'
          WRITE(0,*)
          EXIT
        ENDIF

        CALL FindKnot(nknot,knot,NUMZONE,izone,1,i1,i2)

        IF (i2.NE.-1) THEN
c
c     jdemod - the code appears to assume that the midplane is at 0.0 - this 
c              should probably be replaced with the zc value defining the 
c              center of the confined plasma. 
c
c...      Select the appropriate cell, whichever is closest to the midplane will
c         be the cell associated with the core (which we want to build first):
          IF    (knot(i1)%zcen.LT.0.0.AND.knot(i2)%zcen.LT.0.0) THEN
c
c           jdemod
            IF (output) then 
               WRITE(0,'(a,2i6,4(1x,g12.5))') 'XPT DN:',i1,i2,
     >            knot(i1)%rv(1),knot(i1)%zv(1),
     >            knot(i1)%zcen,knot(i2)%zcen
               WRITE(6,'(a,2i6,4(1x,g12.5))') 'XPT DN:',i1,i2,
     >            knot(i1)%rv(1),knot(i1)%zv(1),
     >            knot(i1)%zcen,knot(i2)%zcen
            endif
c            IF (output) WRITE(0,*) '  >> ',i1,i2
c            IF (output) WRITE(0,*) '     ',knot(i1)%zcen,knot(i2)%zcen
            nxpt = nxpt + 1
            IF (knot(i1)%zcen.GT.knot(i2)%zcen) THEN
              ixpt(nxpt) = i1
              knot(i1)%xpt = i2
              knot(i2)%xpt = i1
            ELSE
              ixpt(nxpt) = i2
              knot(i1)%xpt = i2
              knot(i2)%xpt = i1
            ENDIF
          ELSEIF(knot(i1)%zcen.GT.0.0.AND.knot(i2)%zcen.GT.0.0) THEN
            nxpt = nxpt + 1
c
c           jdemod
            IF (output) then 
               WRITE(0,'(a,2i6,4(1x,g12.5))') 'XPT UP:',i1,i2,
     >            knot(i1)%rv(1),knot(i1)%zv(1),
     >            knot(i1)%zcen,knot(i2)%zcen
               WRITE(6,'(a,2i6,4(1x,g12.5))') 'XPT UP:',i1,i2,
     >            knot(i1)%rv(1),knot(i1)%zv(1),
     >            knot(i1)%zcen,knot(i2)%zcen
            endif
c            IF (output) WRITE(0,*) '  >> ',i1,i2
c            IF (output) WRITE(0,*) '     ',knot(i1)%zcen,knot(i2)%zcen
            IF (knot(i1)%zcen.LT.knot(i2)%zcen) THEN
              ixpt(nxpt) = i1
              knot(i1)%xpt = i2
              knot(i2)%xpt = i1
            ELSE
              ixpt(nxpt) = i2
              knot(i1)%xpt = i2
              knot(i2)%xpt = i1
            ENDIF
          ELSE
            CALL ER('Readgeneralisedgrid','Unrecognized '//
     .              'x-point configuration',*99)
          ENDIF
        ENDIF
      ENDDO

c...  No x-points found, for some reason:
      IF (nxpt.EQ.0) 
     .  CALL ER('Readgeneralisedgrid','No x-points found',*99)

      IF (output) THEN
        DO i1 = 1, nxpt
          WRITE(0,'(a,5i6)') 'XPTS:',i1,
     .               ixpt(i1),knot(ixpt(i1))%index,
     .               knot(ixpt(i1))%xpt,
     .               knot(knot(ixpt(i1))%xpt)%index
          WRITE(6,'(a,5i6)') 'XPTS:',i1,
     .               ixpt(i1),knot(ixpt(i1))%index,
     .               knot(ixpt(i1))%xpt,
     .               knot(knot(ixpt(i1))%xpt)%index
        ENDDO
      ENDIF
c
c     jdemod - the Xpoint finding algorithm returns the two cells that share the 
c              Xpoint vertex as index 1. One of these cells should be the cell
c              in the SOL adjacent to the first cell on the core ring at the
c              Xpoint. The second of these cells is below the Xpoint adjacent to 
c              PFZ. By using the cell "closer to the midplane" it should choose
c              the cell adjacent to the confined plasma - the cell sharing the side
c              with vertices where 1 = 2 and 4 = 3 should be the first cell on the 
c              core ring. 
c
c              The code then walks inward from the Xpoint finding the cell on the
c              innermost core ring corresponding to the first cell on the ring at the 
c              Xpoint. The variable cxpt records the number of rings from the Xpoint
c              to the innermost ring. If this value is the same for two different
c              Xpoints then the grid is a connected double null. If not - the difference
c              in the two values should define the number of rings in the secondary 
c              plasma between the two Xpoints for the DDN plasma configuration. 
c
c...  Searching for the start of the core center ring:
c
c     jdemod - when using cells in the main SOL for finding the center of the core there is
c              a problem when the center point of the cell adjacent to the PFZ is inverted 
c              relative to the Xpoint location. This can happen for grids with a strikepoint that
c              is closer to Z=0.0 than the Xpoint. However, using cells that are actually in 
c              the core and pfz by using a different vertex for Xpoint detection allows this 
c              problem to be avoided. The definition of cxpt needs to be revised to accomodate this.
c


c
      cxpt = 0
      icore = 0
      DO i3 = 1, nxpt
        newi1 = ixpt(i3)
        cont = .TRUE.
        DO WHILE (cont)
          i1 = newi1
          cxpt(i3) = cxpt(i3) + 1
          cont = .FALSE.
c
c     jdemod - when the search has moved all the way inward and 
c              can no longer find a cell with a matching side then
c              i2=-1 is returned and the code moves onto any other
c              Xpoints. 
c
          CALL FindKnot(nknot,knot,NUMZONE,izone,2,i1,i2)

          IF (i2.NE.-1) THEN
            cont = .TRUE.
            newi1 = i2 
            icore(i3) = i2
          ENDIF
        ENDDO
c
c       jdemod
        IF (output) then 
           WRITE(0,*) 'CORE xpt:',i3,ixpt(i3),cxpt(i3),icore(i3)
           WRITE(0,'(a,i6,10(1x,g12.5))') 'CELL1:',ixpt(i3),
     >           knot(ixpt(i3))%rcen,  knot(ixpt(i3))%zcen          
           WRITE(0,'(a,i6,10(1x,g12.5))') 'CELL2:',icore(i3),
     >           knot(icore(i3))%rcen,  knot(icore(i3))%zcen          
           WRITE(6,*) 'CORE xpt:',i3,ixpt(i3),cxpt(i3),icore(i3)
           WRITE(6,'(a,i6,10(1x,g12.5))') 'CELL1:',ixpt(i3),
     >           knot(ixpt(i3))%rcen,  knot(ixpt(i3))%zcen          
           WRITE(6,'(a,i6,10(1x,g12.5))') 'CELL2:',icore(i3),
     >           knot(icore(i3))%rcen,  knot(icore(i3))%zcen          
        endif
c        IF (output) WRITE(0,*) 'Cxpt:',ixpt(i3),cxpt(i3),icore(i3)
      ENDDO

      IF (nxpt.GT.1) THEN
c...    Check that the x-points are ordered properly, with the primary x-point
c       at index 1, and whether or not the double-null grid is connected:
        swap = .FALSE.
        IF     (nxpt.GT.1.AND.cxpt(1).EQ.cxpt(2)) THEN
c...      Connected:
          IF (knot(ixpt(1))%zcen.GT.0.0) swap = .TRUE.
          connected = .TRUE.
c
c         jdemod
          IF (output) then 
             WRITE(0,*) 'CONNECTED DN DETECTED'
             WRITE(6,*) 'CONNECTED DN DETECTED'
          endif
c          IF (output) WRITE(0,*) 'CONNECTED DN DETECTED'
        ELSEIF (nxpt.GT.1.AND.cxpt(1).GT.cxpt(2)) THEN
          swap = .TRUE.
          connected = .FALSE.
        ENDIF
        IF (swap) THEN
c
c         jdemod
          IF (output) then 
             WRITE(0,*) 'SWAPPING X-POINTS'
             WRITE(6,*) 'SWAPPING X-POINTS'
          endif
c          IF (output) WRITE(0,*) 'SWAPPING X-POINTS'
          ixpt (0) = ixpt (1)
          ixpt (1) = ixpt (2)
          ixpt (2) = ixpt (0)
          cxpt (0) = cxpt (1)
          cxpt (1) = cxpt (2)
          cxpt (2) = cxpt (0)
          icore(0) = icore(1)
          icore(1) = icore(2)
          icore(2) = icore(0)
        ENDIF
      ENDIF

c...  Make sure that x-point knot indeces are ordered properly, with the inner (lower radius)
c     of each pair listed in IXPT:
c      DO i1 = 1, nxpt
c        IF (knot(ixpt(i1))%rcen.GT.knot(knot(ixpt(i1))%xpt)%rcen) THEN
c          jdemod
c          IF (output) then 
c              WRITE(0,*) 'SWAPPING X-POINT PAIR',i1
c              WRITE(6,*) 'SWAPPING X-POINT PAIR',i1
c          endif       
c          IF (output) WRITE(0,*) 'SWAPPING X-POINT PAIR',i1
c          ixpt(i1) = knot(ixpt(i1))%xpt
c        ENDIF
c      ENDDO
 
c...  Location of the primary separatrix is known:
      irsep  = cxpt(1)
      irsep2 = irsep

c...  Build the grid:
c
c     jdemod
      IF (output) then 
         WRITE(0,*) 'PROCESSING CORE RINGS'
         WRITE(6,*) 'PROCESSING CORE RINGS'
      endif
c      IF (output) WRITE(0,*) 'PROCESSING CORE RINGS'

c...  Start with the core rings:
      ik = 1
      ir = 1
      i1 = icore(1) 
      DO WHILE(ir.LT.irsep)
        imap(1,ir) = i1
        cont = .TRUE.
        DO WHILE(cont)
          cont = .FALSE.
c...      Move along the ring:
          CALL FindKnot(nknot,knot,NUMZONE,izone,3,i1,i2)
          IF (i2.NE.-1) THEN
            IF (i2.NE.imap(1,ir)) THEN
              i1 = i2 
              ik = ik + 1
              imap(ik,ir) = i1
              cont = .TRUE.
c
c             jdemod
              IF (output) then 
                 WRITE(0,'(a,4i6,5(1x,g12.5))') 'CORE MAP:',
     >                ik,ir,i1,imap(1,ir),
     >                knot(i1)%rcen,knot(i1)%zcen
                 WRITE(6,'(a,4i6,5(1x,g12.5))') 'CORE MAP:',
     >                ik,ir,i1,imap(1,ir),
     >                knot(i1)%rcen,knot(i1)%zcen
              endif
c              IF (output) WRITE(0,*) 'CORE MAP:',ik,ir,i1
            ENDIF
          ELSE
            CALL ER('Readgeneralisedgrid','Bad IK step',*99)
          ENDIF

        ENDDO
        nks(ir) = ik
c...    Step outward, still in the core:        
        CALL FindKnot(nknot,knot,NUMZONE,izone,4,imap(1,ir),i2)
        IF (i2.NE.-1) THEN        
          i1 = i2
          ik = 1
          ir = ir + 1
        ELSE
          CALL ER('Readgeneralisedgrid','Bad IR step',*99)
        ENDIF
      ENDDO

c...  SOL rings: 
c
c     jdemod
      IF (output) then 
         WRITE(0,*) 'PROCESSING SOL RINGS'
         WRITE(6,*) 'PROCESSING SOL RINGS'
      endif
c      IF (output) WRITE(0,*) 'PROCESSING SOL RINGS'
c
c     jdemod - Doesn't this assume an Xpoint down configuration? At least as far as the 
c              "high field side" reference goes? I think the code itself still works. 
c

c...  Step out of the core on the high field side:
      i1 = imap(1,irsep-1)
      CALL FindKnot(nknot,knot,NUMZONE,izone,4,i1,i2)
c
      write(0,*) 'SOL RINGS:',irsep,i1,i2,knot(i2)%rcen,knot(i2)%zcen
c
      IF (i2.NE.-1) THEN  
c...    Move down to the target:
        i1 = i2
        cont = .TRUE.
        DO WHILE(cont)
          cont = .FALSE.
          CALL FindKnot(nknot,knot,NUMZONE,izone,5,i1,i2)
          IF (i2.NE.-1) THEN 
            i1 = i2
            cont = .TRUE.
          ENDIF
        ENDDO
      ELSE
        CALL ER('Readgeneralisedgrid','Bad IR step to SOL',*99)
      ENDIF
c...  Target located, start mapping the SOL:
      ik = 1
      ir = irsep
      imap(ik,ir) = i1
      cont = .TRUE.
      DO WHILE(cont)
        cont = .FALSE.
c...    Move along the ring:
        CALL FindKnot(nknot,knot,NUMZONE,izone,3,i1,i2)
        IF (i2.NE.-1) THEN
          i1 = i2 
          ik = ik + 1
          imap(ik,ir) = i1
          cont = .TRUE.
c
c         jdemod
          IF (output) then 
             WRITE(0,'(a,4i6,5(1x,g12.5))') 'INNER SOL MAP:',
     >            ik,ir,i1,imap(1,ir),
     >            knot(i1)%rcen,knot(i1)%zcen
             WRITE(6,'(a,4i6,5(1x,g12.5))') 'INNER SOL MAP:',
     >            ik,ir,i1,imap(1,ir),
     >            knot(i1)%rcen,knot(i1)%zcen
          endif
c          IF (output) WRITE(0,*) 'INNER SOL MAP:',ik,ir,i1
        ENDIF
c...    Step radially outward if ring is finished:
        IF (.NOT.cont) THEN
          nks(ir) = ik
          i1 = imap(1,ir)
          CALL FindKnot(nknot,knot,NUMZONE,izone,4,i1,i2)          
          IF (i2.NE.-1) THEN
            i1 = i2
            ik = 1
            ir = ir + 1
            imap(ik,ir) = i1
            cont = .TRUE.
c
c           jdemod
            IF (output) then 
             WRITE(0,'(a,4i6,5(1x,g12.5))') 'INNER SOL MAP:NEW RING:',
     >            ik,ir,i1,imap(1,ir),
     >            knot(i1)%rcen,knot(i1)%zcen
             WRITE(6,'(a,4i6,5(1x,g12.5))') 'INNER SOL MAP:NEW RING:',
     >            ik,ir,i1,imap(1,ir),
     >            knot(i1)%rcen,knot(i1)%zcen
            endif
c            IF (output) WRITE(0,*) 'INNER SOL MAP NEW RING:',ik,ir,i1
          ENDIF
        ENDIF
      ENDDO
      irwall = ir
      irtrap = ir + 1
      nrs = ir

      IF (nxpt.GT.1) THEN
c
c       jdemod
        IF (output) then
           WRITE(0,*) 'PROCESSING LOW FIELD SOL'
           WRITE(6,*) 'PROCESSING LOW FIELD SOL'
        endif
c        IF (output) WRITE(0,*) 'PROCESSING LOW FIELD SOL'

c...    Process the low field side looking for any rings that
c       were not processed when looking around the high field side (which
c       is usually the case for double-nulls):

c...    Register all knots that have been mapped to the grid:
        DO ir = 1, nrs
          DO ik = 1, nks(ir)
            knot(imap(ik,ir))%map = 1
          ENDDO
        ENDDO
c
c     jdemod - why not use - i1=imap(nks(irsep-1),irsep-1) ?
c
c...    Start with the first cell on the outer-most core ring and move to 
c       the last cell on the same ring:
        i1 = imap(1,irsep-1)
        CALL FindKnot(nknot,knot,NUMZONE,izone,5,i1,i2)
        IF (i2.EQ.-1) 
     .    CALL ER('Readgeneralisedgrid','Core map problems',*99)
        i1 = i2
c...    Move into the SOL:
        CALL FindKnot(nknot,knot,NUMZONE,izone,4,i1,i2)
        IF (i2.NE.-1) THEN
          i1 = i2
          IF (knot(i2)%map.EQ.1) THEN  
c...        Keep moving outward until a cell with no assigned mapping is
c           found:
            cont = .TRUE.
            DO WHILE(cont)
              cont = .FALSE.
              CALL FindKnot(nknot,knot,NUMZONE,izone,4,i1,i2)
              IF (i2.NE.-1) THEN 
                i1 = i2
                IF (knot(i1)%map.NE.0) cont = .TRUE.
              ELSE
c...            Either an error or a single-null grid is being tested:
                STOP 'SINGLE NULL GRID BEING TESTED OR ERROR?'
              ENDIF
            ENDDO
          ENDIF
        ELSE
          CALL ER('Readgeneralisedgrid','Bad IR step to SOL',*99)
        ENDIF

c...    An unmapped cell has been found, proceed to target:
        cont = .TRUE.
        DO WHILE(cont)
          cont = .FALSE.
          CALL FindKnot(nknot,knot,NUMZONE,izone,5,i1,i2)
          IF (i2.NE.-1) THEN 
            i1 = i2
            cont = .TRUE.
          ENDIF
        ENDDO

c...    Target located, start mapping the low field SOL:
        ik = 1
        ir = nrs + 1
        IF (connected) irsep2 = ir
        imap(ik,ir) = i1
        cont = .TRUE.
        DO WHILE(cont)
          cont = .FALSE.
c...      Move along the ring:
          CALL FindKnot(nknot,knot,NUMZONE,izone,3,i1,i2)
          IF (i2.NE.-1) THEN
            i1 = i2 
            ik = ik + 1
            imap(ik,ir) = i1
            cont = .TRUE.
c
c           jdemod
            IF (output) then 
               WRITE(0,'(a,4i6,5(1x,g12.5))') 'OUTER SOL MAP:',
     >            ik,ir,i1,imap(1,ir),
     >            knot(i1)%rcen,knot(i1)%zcen
               WRITE(6,'(a,4i6,5(1x,g12.5))') 'OUTER SOL MAP:',
     >            ik,ir,i1,imap(1,ir),
     >            knot(i1)%rcen,knot(i1)%zcen
            endif
c            IF (output) WRITE(0,*) 'OUTER SOL MAP:',ik,ir,i1
          ENDIF
c...      Step radially outward if ring is finished:
          IF (.NOT.cont) THEN
c
c           jdemod
            IF (output) then 
               WRITE(0,*) 'STEPPING OUT:',ik,ir,i1
               WRITE(6,*) 'STEPPING OUT:',ik,ir,i1
            endif
c            IF (output) WRITE(0,*) 'STEPPING OUT:',ik,ir,i1
            nks(ir) = ik
            i1 = imap(1,ir)
            CALL FindKnot(nknot,knot,NUMZONE,izone,4,i1,i2)          
            IF (i2.NE.-1) THEN
              i1 = i2
              ik = 1
              ir = ir + 1
              imap(ik,ir) = i1
              cont = .TRUE.
c
c             jdemod
              IF (output) then 
                 WRITE(0,*) 'OUTER SOL MAP NEW RING:',ik,ir,i1
                 WRITE(6,*) 'OUTER SOL MAP NEW RING:',ik,ir,i1
              endif
c              IF (output) WRITE(0,*) 'OUTER SOL MAP NEW RING:',ik,ir,i1
            ELSE
c...          Assume the outer boundary of the grid:
c
c             jdemod
              IF (output) then 
                 WRITE(0,*) 'ASSUMING OUTER GRID BOUNDARY'
                 WRITE(6,*) 'ASSUMING OUTER GRID BOUNDARY'
              endif
c              IF (output) WRITE(0,*) 'ASSUMING OUTER GRID BOUNDARY'
            ENDIF
          ENDIF
        ENDDO
        irwall = ir
        irtrap = ir + 1
        nrs = ir

c...    Register all knots that have been mapped to the grid:
        DO ir = 1, nrs
          DO ik = 1, nks(ir)
            knot(imap(ik,ir))%map = 1
          ENDDO
        ENDDO

c
c       jdemod
        IF (output) then 
           WRITE(0,*) 'PROCESSING SECONDARY PFZ'
           WRITE(6,*) 'PROCESSING SECONDARY PFZ'
        endif
c        IF (output) WRITE(0,*) 'PROCESSING SECONDARY PFZ'
c...    Process the secondary x-point PFR, which is just considered part of the
c       SOL for generalized grids:
        i1 = knot(ixpt(2))%xpt
c        i1 = ixpt(2)
c...    Move into the PFR:
        CALL FindKnot(nknot,knot,NUMZONE,izone,2,i1,i2)
        IF (i2.EQ.-1) 
     .    CALL ER('Readgeneralisedgrid','PFR2 problems',*99)
        i1 = i2
        ik = 1
        ir = nrs + 1
        imap(ik,ir) = i1

c...    Proceed to target:
c
c       jdemod
        IF (output) then 
           WRITE(0,*) 'LOOKING FOR TARGET'
           WRITE(6,*) 'LOOKING FOR TARGET'
        endif
c        IF (output) WRITE(0,*) 'LOOKING FOR TARGET'
        cont = .TRUE.
        DO WHILE(cont)
          cont = .FALSE.
          CALL FindKnot(nknot,knot,NUMZONE,izone,5,i1,i2)
c
c         jdemod
          IF (output) then 
             WRITE(0,*) 'MOVING',i1,i2,istart
             WRITE(6,*) 'MOVING',i1,i2,istart
          endif
c          IF (output) WRITE(0,*) 'MOVING',i1,i2,istart
          IF (i2.NE.-1) THEN 
            i1 = i2
            ik = ik + 1
            imap(ik,ir) = i1
            cont = .TRUE.
          ENDIF
        ENDDO
 
c...    Target located, start mapping the secondary PFR:
c
c       jdemod
        IF (output) then 
           WRITE(0,*) 'TARGET LOCATED'
           WRITE(6,*) 'TARGET LOCATED'
        endif
c        IF (output) WRITE(0,*) 'TARGET LOCATED'
        ik = 1
        ir = nrs + 1
        imap(ik,ir) = i1
        cont = .TRUE.
        DO WHILE(cont)
          cont = .FALSE.
c...      Move along the ring:
          CALL FindKnot(nknot,knot,NUMZONE,izone,3,i1,i2)
c
c         jdemod
          IF (output) then 
             WRITE(0,*) 'MOVING',i1,i2
             WRITE(6,*) 'MOVING',i1,i2
          endif
c          IF (output) WRITE(0,*) 'MOVING',i1,i2
          IF (i2.NE.-1) THEN
            i1 = i2 
            ik = ik + 1
            imap(ik,ir) = i1
            cont = .TRUE.
c
c           jdemod
            IF (output) then 
               WRITE(0,*) '2ND PFR MAP:',ik,ir,i1,knot(i1)%zcen
               WRITE(6,*) '2ND PFR MAP:',ik,ir,i1,knot(i1)%zcen
            endif
c            IF (output) WRITE(0,*) '2ND PFR MAP:',ik,ir,i1,knot(i1)%zcen
          ENDIF
c...      Step radially outward if ring is finished:
          IF (.NOT.cont) THEN
            nks(ir) = ik
            i1 = imap(1,ir)
            CALL FindKnot(nknot,knot,NUMZONE,izone,2,i1,i2)          
            IF (i2.NE.-1) THEN
              i1 = i2
              ik = 1
              ir = ir + 1
              imap(ik,ir) = i1
              cont = .TRUE.
c
c             jdemod
              IF (output) then  
                WRITE(0,*) '2ND PFR MAP NEW RING:',ik,ir,i1,
     >             knot(i1)%zcen
                WRITE(6,*) '2ND PFR MAP NEW RING:',ik,ir,i1,
     >             knot(i1)%zcen
              endif
c              IF (output) 
c     .          WRITE(0,*) '2ND PFR MAP NEW RING:',ik,ir,i1,
c     .          knot(i1)%zcen
            ELSE
c...          Assume the outer boundary of the grid:
            ENDIF
          ENDIF
        ENDDO
        irwall = ir
        irtrap = ir + 1
        nrs = ir

      ENDIF ! Done processing double-null rings





      IF (.FALSE.) THEN   ! DEBUG
        nks(ir) = ik
        irwall = ir
        irtrap = ir + 1
        nrs = ir
        id = 0
        DO ir = 1, nrs
          DO ik = 1, nks(ir)        
            i1 = imap(ik,ir)
            rs(ik,ir) = knot(i1)%rcen
            zs(ik,ir) = knot(i1)%zcen
            bratio(ik,ir) = knot(i1)%bratio
            id = id + 1
            korpg(ik,ir) = id
            nvertp(id) = knot(i1)%nv
            DO i2 = 1, nvertp(id)
              rvertp(i2,id) = knot(i1)%rv(i2)
              zvertp(i2,id) = knot(i1)%zv(i2)
            ENDDO
          ENDDO
        ENDDO
        ikto = -1
        ikti = -1
        DO ik = 1, nks(irsep)
          IF (connected) THEN
            IF (imap(ik,irsep).EQ.ixpt(1)          ) ikto = ik    ! Not sure this will always work... 
            IF (imap(ik,irsep).EQ.knot(ixpt(2))%xpt) ikti = ik -1
c            IF (imap(ik,irsep).EQ.ixpt(1)) ikto = ik - 1  ! Not sure this will always work...
c            IF (imap(ik,irsep).EQ.ixpt(2)) ikti = ik     
          ENDIF
        ENDDO
        CALL SaveSolution
        CALL OutputData(86,'MAST!')
        title = '...'
        desc  = 'Call to STORE from DumpGrid'
        job   = 'Call to STORE from DumpGrid'
        equil = 'Call to STORE from DumpGrid'
        WRITE(0,*) 'CALLING STORE'
        CALL Store(title,desc,1,job,equil,facta,factb,1,1)
        WRITE(0,*) 'FUN WITH MAST GRIDS!'
        STOP
      ENDIF


c
c     jdemod
      IF (output) then 
         WRITE(0,*) 'PROCESSING PRIMARY PFZ'
         WRITE(6,*) 'PROCESSING PRIMARY PFZ'
      endif   
c      IF (output) WRITE(0,*) 'PROCESSING PRIMARY PFZ'

c...  Process the primary x-point PFR:
      i1 = knot(ixpt(1))%xpt
c...  Move into the PFR:
      CALL FindKnot(nknot,knot,NUMZONE,izone,2,i1,i2)
      IF (i2.EQ.-1) 
     .  CALL ER('Readgeneralisedgrid','PFR1 problems',*99)
c
c     jdemod
      IF (output) then
         WRITE(0,*) '  INTO PFZ',knot(i1)%index,knot(i2)%index
         WRITE(6,*) '  INTO PFZ',knot(i1)%index,knot(i2)%index
      endif
c      IF (output) 
c     .  WRITE(0,*) '  INTO PFZ',knot(i1)%index,knot(i2)%index
      i1 = i2
c...  Proceed to target:
      cont = .TRUE.
      DO WHILE(cont)
        cont = .FALSE.
        CALL FindKnot(nknot,knot,NUMZONE,izone,5,i1,i2)
c        jdemod
c        IF (output) then 
c           WRITE(0,*) '  TO TARGET',i1,i2
c           WRITE(0,*) '  TO TARGET',knot(i1)%index
c           WRITE(0,*) '  TO TARGET',knot(i2)%index
c           WRITE(6,*) '  TO TARGET',i1,i2
c           WRITE(6,*) '  TO TARGET',knot(i1)%index
c           WRITE(6,*) '  TO TARGET',knot(i2)%index
c        endif
c        IF (output) WRITE(0,*) '  TO TARGET',i1,i2
c        IF (output) WRITE(0,*) '  TO TARGET',knot(i1)%index
c        IF (output) WRITE(0,*) '  TO TARGET',knot(i2)%index
c        STOP 'tet'
        IF (i2.NE.-1) THEN 
          i1 = i2
          cont = .TRUE.
        ENDIF
      ENDDO
c...  Target located, start mapping the primary PFR:
c     
c     jdemod
      IF (output) then 
         WRITE(0,*) 'PROCESSING PRIMARY PFZ, TARGET LOCATED'
         WRITE(6,*) 'PROCESSING PRIMARY PFZ, TARGET LOCATED'
      endif
c      IF (output) WRITE(0,*) 'PROCESSING PRIMARY PFZ, TARGET LOCATED'
      ik = 1
      ir = nrs + 1
      imap(ik,ir) = i1
      cont = .TRUE.
      DO WHILE(cont)
        cont = .FALSE.
c...    Move along the ring:
        CALL FindKnot(nknot,knot,NUMZONE,izone,3,i1,i2)
        IF (i2.NE.-1) THEN
          i1 = i2 
          ik = ik + 1
          imap(ik,ir) = i1
          cont = .TRUE.
c          jdemod
c          IF (output) then 
c             WRITE(0,*) '1ST PFR MAP:',ik,ir,i1,knot(i1)%zcen
c             WRITE(6,*) '1ST PFR MAP:',ik,ir,i1,knot(i1)%zcen
c          endif
c          IF (output) WRITE(0,*) '1ST PFR MAP:',ik,ir,i1,knot(i1)%zcen
        ENDIF
c...    Step radially outward if ring is finished:
        IF (.NOT.cont) THEN
          nks(ir) = ik
          i1 = imap(1,ir)
          CALL FindKnot(nknot,knot,NUMZONE,izone,2,i1,i2)          
          IF (i2.NE.-1) THEN
            i1 = i2
            ik = 1
            ir = ir + 1
            imap(ik,ir) = i1
            cont = .TRUE.
c
c           jdemod
            IF (output) then 
             WRITE(0,'(a,4i6,5(1x,g12.5))') 'PRIM PFZ MAP:',
     >            ik,ir,i1,imap(1,ir),
     >            knot(i1)%rcen,knot(i1)%zcen
             WRITE(6,'(a,4i6,5(1x,g12.5))') 'PRIM PFZ MAP:',
     >            ik,ir,i1,imap(1,ir),
     >            knot(i1)%rcen,knot(i1)%zcen
            endif
c            IF (output) 
c     .        WRITE(0,*) '1ST PFR MAP NEW RING:',ik,ir,i1,knot(i1)%zcen
          ELSE
c...        Assume the outer boundary of the grid:
          ENDIF
c         jdemod
          IF (output) then 
             WRITE(0,*) 'PROCESSING PRIMARY PFZ, BUZZING...'
             WRITE(6,*) 'PROCESSING PRIMARY PFZ, BUZZING...'
          endif
c          IF (output) WRITE(0,*) 'PROCESSING PRIMARY PFZ, BUZZING...'
        ENDIF
      ENDDO
      nrs = ir

c...  Need to reorder the rings in the primary PFZ:
      DO i1 = 0, (nrs-irtrap+1)/2-1
c
c       jdemod
        IF (output) then 
           WRITE(0,*) 'I1???=',i1
           WRITE(6,*) 'I1???=',i1
        endif
c        IF (output) WRITE(0,*) 'I1???=',i1
        tmpnks = nks(irtrap+i1)
        DO ik = 1, nks(irtrap+i1)
          imap(ik,0) = imap(ik,irtrap+i1)
        ENDDO
        nks(irtrap+i1) = nks(nrs-i1)
        DO ik = 1, nks(nrs-i1)
          imap(ik,irtrap+i1) = imap(ik,nrs-i1)
        ENDDO
        nks(nrs-i1) = tmpnks
        DO ik = 1, tmpnks
          imap(ik,nrs-i1) = imap(ik,0)
        ENDDO
      ENDDO

c...  Put grid together:
c     jdemod
      IF (output) then 
         WRITE(0,*) 'PUTTING GRID TOGETHER'
         WRITE(6,*) 'PUTTING GRID TOGETHER'
      endif
c      IF (output) WRITE(0,*) 'PUTTING GRID TOGETHER'
      id = 0
      DO ir = 1, nrs
        DO ik = 1, nks(ir)        
          i1 = imap(ik,ir)
          rs(ik,ir) = knot(i1)%rcen
          zs(ik,ir) = knot(i1)%zcen
          bratio(ik,ir) = knot(i1)%bratio
          id = id + 1
          korpg(ik,ir) = id
          nvertp(id) = knot(i1)%nv
          DO i2 = 1, nvertp(id)
            rvertp(i2,id) = knot(i1)%rv(i2)
            zvertp(i2,id) = knot(i1)%zv(i2)
          ENDDO
c...      Store these in case B2 data from Rhozansky is being loaded:
C         IPP/11 - Karl: put in check for divimp_ik allocated
          IF (ALLOCATED(divimp_ik)) THEN
            divimp_ik(ik,ir) = knot(i1)%ik 
            divimp_ir(ik,ir) = knot(i1)%ir
          ENDIF
        ENDDO
      ENDDO

c...  Set NPOLYP:
      npolyp  = id
      vpolmin = (MAXNKS*MAXNRS - npolyp) / 2 + npolyp
      vpolyp  = vpolmin

c...  Find IKTO,IKTI:
      ikto = -1
      ikti = -1
      DO ik = 1, nks(irsep)
        IF (connected) THEN
          IF (imap(ik,irsep).EQ.ixpt(1)          ) ikto = ik    ! Not sure this will always work...
          IF (imap(ik,irsep).EQ.knot(ixpt(2))%xpt) ikti = ik - 1
c          IF (imap(ik,irsep).EQ.ixpt(1)) ikto = ik      ! Not sure this will always work...
c          IF (imap(ik,irsep).EQ.ixpt(2)) ikti = ik - 1
c          IF (imap(ik,irsep).EQ.knot(ixpt(2))%xpt) ikti = ik - 1  ! Not sure this will always work...
        ELSE
          IF (imap(ik,irsep).EQ.ixpt(1)          ) ikto = ik - 1
          IF (imap(ik,irsep).EQ.knot(ixpt(1))%xpt) ikti = ik 
        ENDIF
      ENDDO


      IF (ikto.EQ.-1.OR.ikti.EQ.-1)
     .  CALL ER('Readgeneralisedgrid','IKTI or IKTO not found',*99)

c...  Add virtual rings 1 (core boundary), IRWALL (SOL) and IRTRAP (PFZ):
      CALL InsertRing(1,BEFORE,PERMANENT)
      CALL InsertRing(nrs-irsep+2,AFTER,PERMANENT)
      CALL InsertRing(nrs-irsep+3,BEFORE,PERMANENT)

c     CUTPT1
c     CUTPT2
c     MAXKPTS
c     MAXRINGS
c     CUTRING
c
      nopriv = .FALSE.

      cutpt1 = ikto
      cutpt2 = ikti             ! These are semi-bogus for a connected double-null...?
      cutring = irsep - 1
      maxkpts = nks(irsep)
      maxrings = irwall
      indexiradj = 1

c...TMP
c      ik = 1
c      ir = nrs + 1
c      nks(ir) = 1 
c      imap(ik,ir) = i2
c      irwall = ir
c      irtrap = ir + 1
c      nrs = ir 
c      WRITE(0,*) 'i2:',i2,knot(i2)%map

      IF (.NOT..TRUE.) THEN
c        id = 0
c        DO ir = 1, nrs
c          DO ik = 1, nks(ir)        
c            i1 = imap(ik,ir)
c            rs(ik,ir) = knot(i1)%rcen
c            zs(ik,ir) = knot(i1)%zcen
c            bratio(ik,ir) = knot(i1)%bratio
c            id = id + 1
c            korpg(ik,ir) = id
c            nvertp(id) = knot(i1)%nv
c            DO i2 = 1, nvertp(id)
c              rvertp(i2,id) = knot(i1)%rv(i2)
c              zvertp(i2,id) = knot(i1)%zv(i2)
c            ENDDO
c          ENDDO
c        ENDDO
        CALL SaveSolution
        CALL OutputData(86,'MAST!')
        title = '...'
        desc  = 'Call to STORE from DumpGrid'
        job   = 'Call to STORE from DumpGrid'
        equil = 'Call to STORE from DumpGrid'
        WRITE(0,*) 'CALLING STORE'
        CALL Store(title,desc,1,job,equil,facta,factb,1,1)
        WRITE(0,*) 'FUN WITH MAST GRIDS!'
        STOP
      ENDIF

c...  Add virtual boundary cells, which will be stripped off later:
      IF (ctargopt.EQ.0.OR.ctargopt.EQ.1.OR.ctargopt.EQ.2.OR.
     .    ctargopt.EQ.3.OR.ctargopt.EQ.6) 
     .  CALL AddPoloidalBoundaryCells      

c...  Look for PSIn data for full double null grids (code mostly 
c     from tau.d6a):
      READ(gridunit,'(A)',END=25) buffer
      IF (buffer(1:16).EQ.'PSI-DOUBLE-NULLd') THEN   ! direct assignment to each ring, post mortem...
        READ(buffer(17:),*) numpsi
c...    The PSI values are to be loaded in TailorGrid...glorious hack!
        DO i1 = 1, numpsi
          READ(gridunit,*,END=97)
        ENDDO
      ELSEIF (buffer(1:15).EQ.'PSI-DOUBLE-NULL') THEN
        READ(buffer(16:),*) numpsi
c...    The PSI values are listed with one on each line
c       indexed by knot and ring index based on the SONNET 
c       grid coordinates:
        DO i1 = 1, numpsi
          READ(gridunit,*,END=97) ikpsi(i1),irpsi(i1),valpsi(i1)
        ENDDO
c...    Assign to grid rings:
        DO ir = 2, irwall-1
c...      (Need the "-1" because a virtual core ring has been added to the grid)
          ir1 = knot(imap(1,ir-1))%ir 
          DO i1 = 1, numpsi
            IF (irpsi(i1).EQ.ir1) EXIT
          ENDDO
          IF (i1.EQ.numpsi+1) 
     .      CALL ER('Readgeneralisedgrid','Problem with PSIn',*99)
c          WRITE(0,*) '--',ir,valpsi(i1)
          psitarg(ir,1) = valpsi(i1)
          psitarg(ir,2) = valpsi(i1)
          IF (ir.LT.irsep) THEN
            psitarg(ir-1+irtrap,1) = valpsi(i1)
            psitarg(ir-1+irtrap,2) = valpsi(i1)
          ENDIF
        ENDDO
        WRITE(0,*)
        WRITE(0,*) '--------------------------------------------------'
        WRITE(0,*) 'BOGUS PSITARG -- ALSO USING INNER TARGET DATA ONLY'
        WRITE(0,*) '--------------------------------------------------'
        WRITE(0,*)
      ELSE
        BACKSPACE gridunit
      ENDIF
 25   CONTINUE

      DEALLOCATE(knot)
      DEALLOCATE(imap)

c      IF (nrs.EQ.60) THEN
c        WRITE(0,*)
c        WRITE(0,*) '--------------------------------------------------'
c        WRITE(0,*) 'HARDCODING IRSEP2 = 30 (FOR IR = 60) '
c        WRITE(0,*) '--------------------------------------------------'
c        WRITE(0,*)
c        irsep2 = 35
c      ENDIF

c...  For consistency with original SONNET code in tau.d6a:
      ir = maxrings
      ik = maxkpts

      RETURN
 97   CALL ER('Readgeneralisedgrid','Unexpected end-of-file',*99)
 98   CALL ER('Readgeneralisedgrid','Problem accessing grid file',*99)
 99   WRITE(0,*) 'IXPT:',ixpt(1),ixpt(2)
      STOP  
      END
c
c
c ========================================================================
c
c jdemod - nopriv moved to common cgeom
c
c      SUBROUTINE ReadQuasiDoubleNull(gridunit,nopriv,ik,ir,
c     .                               rshift,zshift,indexiradj)
      SUBROUTINE ReadQuasiDoubleNull(gridunit,ik,ir,
     .                               rshift,zshift,indexiradj)
      IMPLICIT none

      INTEGER gridunit,ik,ir,indexiradj
c      LOGICAL nopriv

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'
      INCLUDE 'pindata'

      REAL*8 TOL
      PARAMETER (TOL=1.0D-06)

      INTEGER i1,i2,ik1,ik2,ik3,ir2,ikmax,irmax,id1,id2,in,id,
     .        tmpnks,iks,ike,count
      REAL    rshift,zshift
      CHARACTER*100 buffer

      REAL*8, ALLOCATABLE :: rcen(:,:),zcen(:,:),bval(:,:),
     .                       rvp(:,:,:),zvp(:,:,:)

 
      quasidn = .TRUE.

c...  Allocate the temporary storage arrays for the grid:
      ALLOCATE(rvp(0:MAXNKS,0:MAXNRS,5))
      ALLOCATE(zvp(0:MAXNKS,0:MAXNRS,5))
      ALLOCATE(rcen(0:MAXNKS,0:MAXNRS))
      ALLOCATE(zcen(0:MAXNKS,0:MAXNRS))
      ALLOCATE(bval(0:MAXNKS,0:MAXNRS))

c...  Process the grid file:
      ikmax = 0
      irmax = 0

      DO WHILE (.TRUE.)
        READ(gridunit,'(A100)',END=290) buffer

        IF     (buffer(1:6).EQ.'SHIFT:'.OR.
     .          buffer(1:6).EQ.'Shift:'.OR.
     .          buffer(1:6).EQ.'shift:') THEN

c...      Get data for shifting the grid: 
          READ (buffer(7:),*) rshift,zshift

        ELSEIF (buffer(1:5).EQ.'GEOM:'.OR.
     .          buffer(1:5).EQ.'Geom:'.OR.
     .          buffer(1:5).EQ.'geom:') THEN

c...      Get grid index data: 
          READ (buffer(6:),*) irmax,irsep,irsep2,ikmax,cutpt1,cutpt2

        ELSEIF (buffer(4:8).EQ.'====='.OR.buffer(4:8).EQ.'-----') THEN

c...      Check that cell data exits on the next line:
          READ(gridunit,'(A100)',END=290) buffer
          BACKSPACE(gridunit)
          IF (buffer(1:10).NE.'   Element') EXIT

c...      Read in cell data:
          READ(gridunit,9000) in,ik,ir,
     .                        rvp(ik,ir,2),zvp(ik,ir,2),
     .                        rvp(ik,ir,3),zvp(ik,ir,3)
          READ(gridunit,9001) bval(ik,ir),rcen(ik,ir),zcen(ik,ir)
          READ(gridunit,9002) rvp(ik,ir,1),zvp(ik,ir,1),
     .                        rvp(ik,ir,4),zvp(ik,ir,4)
c          ikmax = MAX(ik,ikmax)
c          irmax = MAX(ir,irmax)
        ENDIF
      ENDDO
290   CONTINUE      

      nrs = irmax
      nks = ikmax

c...  Shift everything by RSHIFT and ZSHIFT:
      DO ir = 1, irmax
        DO ik = 1, ikmax       
          rcen(ik,ir) = rcen(ik,ir) + DBLE(rshift)
          zcen(ik,ir) = zcen(ik,ir) + DBLE(zshift)
          DO i1 = 1, 4
            rvp(ik,ir,i1) = rvp(ik,ir,i1) + DBLE(rshift)
            zvp(ik,ir,i1) = zvp(ik,ir,i1) + DBLE(zshift)
          ENDDO 
        ENDDO
      ENDDO

c...  Need to clean up the cells along the midplane, since the "top" and "bottom"
c     of the UEDGE grid do not quite match up (a very small shift is required).  Note
c     that the ZSHIFT applied to the grid has to be reasonably accurate in order for
c     the cell boundaries to actually lie along the Z=0.0 line, since the grid is
c     not made with Z=0.0 at the midplane.  The 2 cm bound used below is arbitrary, and 
c     is representative of the accuracy with which ZSHIFT moves the appropriate cell
c     boundaries to the midplane:
      DO ir = 1, irmax
        count = 0
        DO ik = 1, ikmax-1       
          IF (DABS(zvp(ik,ir,3)).LT.5.0D-02.AND.
     .        DABS(zvp(ik,ir,4)).LT.5.0D-02) THEN
            rvp(ik+1,ir,2) = rvp(ik,ir,3)
            zvp(ik+1,ir,2) = zvp(ik,ir,3)
            rvp(ik+1,ir,1) = rvp(ik,ir,4)
            zvp(ik+1,ir,1) = zvp(ik,ir,4)
            count = count + 1
          ENDIF
        ENDDO
        IF (count.NE.2) THEN
          WRITE(0,*) 'DEBUG:',ir,irsep
          CALL ER('ReadQuasiDoubleNull','Unable to correct midplane '//
     .            'cells',*99)
        ENDIF
      ENDDO

c...  Extract cells in the divertor and adjacent to the core in the SOL 
c     classic (between IRSEP and IRSEP2).  Reorganize rings so that 
c     neighbouring poloidal cells are also next to each other in index space:
      DO ir = irsep, nrs
        DO ik1 = 1, ikmax-1
          DO ik2 = ik1+1, ikmax
            IF (DABS(rvp(ik1,ir,3)-rvp(ik2,ir,2)).LT.TOL.AND.
     .          DABS(zvp(ik1,ir,3)-zvp(ik2,ir,2)).LT.TOL.AND.
     .          DABS(rvp(ik1,ir,4)-rvp(ik2,ir,1)).LT.TOL.AND.
     .          DABS(zvp(ik1,ir,4)-zvp(ik2,ir,1)).LT.TOL) THEN
c...          Found adjacent cell along the ring.  Move it next to the 
c             current focus cell (IK1,IR), if necessary:
              IF (ik1+1.NE.ik2) THEN
                rcen(0,0) = rcen(ik1+1,ir)
                zcen(0,0) = zcen(ik1+1,ir)
                bval(0,0) = bval(ik1+1,ir)
                rcen(ik1+1,ir) = rcen(ik2,ir)
                zcen(ik1+1,ir) = zcen(ik2,ir)
                bval(ik1+1,ir) = bval(ik2,ir)
                rcen(ik2,ir) = rcen(0,0)
                zcen(ik2,ir) = zcen(0,0)
                bval(ik2,ir) = bval(0,0)
                DO i1 = 1, 4
                  rvp(0,0,i1) = rvp(ik1+1,ir,i1)
                  zvp(0,0,i1) = zvp(ik1+1,ir,i1)
                  rvp(ik1+1,ir,i1) = rvp(ik2,ir,i1)
                  zvp(ik1+1,ir,i1) = zvp(ik2,ir,i1)
                  rvp(ik2,ir,i1) = rvp(0,0,i1)
                  zvp(ik2,ir,i1) = zvp(0,0,i1)
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDDO
c...    Set NKS for the ring by cutting the ring at the first virtual cell:
        DO ik = 1, ikmax
          IF (DABS(rvp(ik,ir,3)-rvp(ik,ir,2)).LT.TOL.AND.
     .        DABS(zvp(ik,ir,3)-zvp(ik,ir,2)).LT.TOL.AND.
     .        DABS(rvp(ik,ir,4)-rvp(ik,ir,1)).LT.TOL.AND.
     .        DABS(zvp(ik,ir,4)-zvp(ik,ir,1)).LT.TOL) THEN
            nks(ir) = ik - 1
            EXIT
          ENDIF
        ENDDO
      ENDDO

c...  Assign to OEDGE grid geometry arrays for "standard" rings:
      id = 0
      DO ir = 1, nrs 
        DO ik = 1, nks(ir)
          rs(ik,ir) = SNGL(rcen(ik,ir))
          zs(ik,ir) = SNGL(zcen(ik,ir))
          bratio(ik,ir) = SNGL(bval(ik,ir))
          id = id + 1
          korpg(ik,ir) = id
          nvertp(id) = 4
          DO i1 = 1, 4
            rvertp(i1,id) = SNGL(rvp(ik,ir,i1))
            zvertp(i1,id) = SNGL(zvp(ik,ir,i1))
          ENDDO
        ENDDO
      ENDDO

      DO ir = 1, nrs
        irorg2(ir) = ir 
      ENDDO

c...  Add outer half of rings with IR greater than IRSEP2+1 to the SOL.  Identify 
c     the virtual cells bounding the poloidal region of interest:
      DO ir = irsep2+1, nrs       
        iks = 0
        ike = 0
        DO ik = nks(ir)+2, ikmax
          IF (DABS(rvp(ik,ir,3)-rvp(ik,ir,2)).LT.TOL.AND.
     .        DABS(zvp(ik,ir,3)-zvp(ik,ir,2)).LT.TOL.AND.
     .        DABS(rvp(ik,ir,4)-rvp(ik,ir,1)).LT.TOL.AND.
     .        DABS(zvp(ik,ir,4)-zvp(ik,ir,1)).LT.TOL) THEN
            IF (iks.EQ.0) iks = ik + 1     
            IF (iks.NE.0) ike = ik - 1       
          ENDIF
        ENDDO

        ir2 = nrs + ir - irsep2
        irorg2(ir2) = ir
        nks(ir2) = ike - iks + 1
        DO ik = iks, ike
          ik2 = ik - iks + 1
          rs(ik2,ir2) = SNGL(rcen(ik,ir))
          zs(ik2,ir2) = SNGL(zcen(ik,ir))
          bratio(ik2,ir2) = SNGL(bval(ik,ir))
          id = id + 1
          korpg(ik2,ir2) = id
          nvertp(id) = 4
          DO i1 = 1, 4
            rvertp(i1,id) = SNGL(rvp(ik,ir,i1))
            zvertp(i1,id) = SNGL(zvp(ik,ir,i1))
          ENDDO
        ENDDO
      ENDDO
      nrs = 2 * nrs - irsep2

c...  Add (classic) PFZ rings:
      DO ir = 1, irsep-1
        ir2 = nrs + ir
        nks(ir2) = cutpt1 + nks(ir) - cutpt2
        irorg2(ir2) = ir2
        ik2 = 0
        DO ik = 1, nks(ir)
          IF (ik.GT.cutpt1.AND.ik.LT.cutpt2) CYCLE
          ik2 = ik2 + 1
          rs(ik2,ir2) = SNGL(rcen(ik,ir))
          zs(ik2,ir2) = SNGL(zcen(ik,ir))
          bratio(ik2,ir2) = SNGL(bval(ik,ir))
          id = id + 1
          korpg(ik2,ir2) = id
          nvertp(id) = 4
          DO i1 = 1, 4
            rvertp(i1,id) = SNGL(rvp(ik,ir,i1))
            zvertp(i1,id) = SNGL(zvp(ik,ir,i1))
          ENDDO
c          WRITE(0,*) 'ADDING:',ik2,ir2
        ENDDO
      ENDDO
      nrs = nrs + irsep - 1

c...  Adjust core rings:
      DO ir = 1, irsep-1
        DO ik = 1, nks(ir) - cutpt1
          rs    (ik,ir) = rs    (ik+cutpt1,ir)
          zs    (ik,ir) = zs    (ik+cutpt1,ir)
          bratio(ik,ir) = bratio(ik+cutpt1,ir)
          korpg (ik,ir) = korpg (ik+cutpt1,ir)
        ENDDO
        nks(ir) = cutpt2 - cutpt1 - 1

c...    Make core ring continuous:
        DO ik1 = 1, nks(ir)-1
          id1 = korpg(ik1,ir)
          DO ik2 = ik1+1, nks(ir)
            id2 = korpg(ik2,ir)
            IF (ABS(rvertp(3,id1)-rvertp(2,id2)).LT.SNGL(TOL).AND.
     .          ABS(zvertp(3,id1)-zvertp(2,id2)).LT.SNGL(TOL).AND.
     .          ABS(rvertp(4,id1)-rvertp(1,id2)).LT.SNGL(TOL).AND.
     .          ABS(zvertp(4,id1)-zvertp(1,id2)).LT.SNGL(TOL)) THEN
c...          Found adjacent cell along the ring.  Move it next to the 
c             current focus cell (IK1,IR), if necessary:
              rs(MAXNKS,MAXNRS) = rs(ik1+1,ir)
              zs(MAXNKS,MAXNRS) = zs(ik1+1,ir)
              rs(ik1+1 ,ir)     = rs(ik2,ir)
              zs(ik1+1 ,ir)     = zs(ik2,ir)
              rs(ik2   ,ir)     = rs(MAXNKS,MAXNRS)
              zs(ik2   ,ir)     = zs(MAXNKS,MAXNRS)
              bratio(MAXNKS,MAXNRS) = bratio(ik1+1,ir)
              bratio(ik1+1 ,ir)     = bratio(ik2,ir)
              bratio(ik2   ,ir)     = bratio(MAXNKS,MAXNRS)
              korpg(MAXNKS,MAXNRS) = korpg(ik1+1,ir)
              korpg(ik1+1 ,ir)     = korpg(ik2,ir)
              korpg(ik2   ,ir)     = korpg(MAXNKS,MAXNRS)
              EXIT
            ENDIF
          ENDDO
c          WRITE(0,*) 'SORTING:',ik1,ik2,nks(ir),tmpnks
        ENDDO

c...    Set NKS by stopping ring at the first virtual cell, or at the
c       first discontinuity:
        tmpnks = nks(ir)
        nks(ir) = 0
        DO ik = 1, tmpnks
          id1 = korpg(ik,ir)
          id2 = korpg(MIN(ik+1,tmpnks),ir)  
          IF (ABS(rvertp(3,id1)-rvertp(2,id1)).LT.SNGL(TOL).AND.
     .        ABS(zvertp(3,id1)-zvertp(2,id1)).LT.SNGL(TOL).AND.
     .        ABS(rvertp(4,id1)-rvertp(1,id1)).LT.SNGL(TOL).AND.
     .        ABS(zvertp(4,id1)-zvertp(1,id1)).LT.SNGL(TOL)) THEN
            nks(ir) = ik - 1
            EXIT
          ENDIF
          IF (.NOT.(ik.LT.tmpnks.AND.
     .        ABS(rvertp(3,id1)-rvertp(2,id2)).LT.SNGL(TOL).AND.
     .        ABS(zvertp(3,id1)-zvertp(2,id2)).LT.SNGL(TOL).AND.
     .        ABS(rvertp(4,id1)-rvertp(1,id2)).LT.SNGL(TOL).AND.
     .        ABS(zvertp(4,id1)-zvertp(1,id2)).LT.SNGL(TOL))) THEN     
            nks(ir) = ik 
            EXIT
          ENDIF
        ENDDO

c...    Close core ring by adding cell to end of ring:
c        nks(ir) = nks(ir) + 1
c        rs(nks(ir),ir) = rs(1,ir)
c        zs(nks(ir),ir) = zs(1,ir)
c        bratio(nks(ir),ir) = bratio(1,ir)
c        korpg(nks(ir),ir) = korpg(1,ir)
      ENDDO

c...  Set NPOLYP:
      npolyp = 0
      DO ir = 1, nrs
        DO ik = 1, nks(ir)
          IF (idring(ir).EQ.-1) CYCLE
          npolyp = MAX(npolyp,korpg(ik,ir))
        ENDDO
      ENDDO
      vpolmin = (MAXNKS*MAXNRS - npolyp) / 2 + npolyp
      vpolyp  = vpolmin

      CALL OutputData(86,'Before adding virtual rings')

c...  Add virtual rings 1 (core boundary), IRWALL (SOL) and IRTRAP (PFZ):
      CALL InsertRing(1,BEFORE,PERMANENT)
      CALL InsertRing(nrs-irsep+2,AFTER,PERMANENT)
      CALL InsertRing(nrs-irsep+3,BEFORE,PERMANENT)

c      CALL ResetRing(1,2)
c      CALL ResetRing(irwall,irwall-1)
c      CALL ResetRing(irtrap,irtrap+1)


c...  Set grid quantities:
c
c     IDRING
c     IRWALL
c     IRTRAP
c     NOPRIV
c     IKTO
c     IKTI
c     NPOLYP
c   

      irwall = nrs - irsep + 1 
      irwall2 = irwall
      irtrap = irwall + 1
      irtrap2 = irtrap

      IF (grdnmod.EQ.0) grdnmod = -1
      irbreak =  irsep2 + 1
      nbr = irwall - irbreak 

      idring = TARTOTAR
      idring(1) = -1
      idring(irwall) = -1
      idring(irtrap) = -1
      i1 = (irwall - 1 - irsep2) / 2
      DO ir = irsep2+1, irwall-1
        IF (ir.LE.irsep2+i1) THEN
          idring(ir) = TARTOWAL
        ELSE
          idring(ir) = WALTOTAR
        ENDIF
      ENDDO

c...  Find IKTO and IKTI from the separatrix ring:
      ir = irsep
      DO ik1 = 1, nks(ir)-2
        id1 = korpg(ik1,ir)
        DO ik2 = ik1+2, nks(ir)
          id2 = korpg(ik2,ir)
          IF (ABS(rvertp(4,id1)-rvertp(1,id2)).LT.TOL.AND.
     .        ABS(zvertp(4,id1)-zvertp(1,id2)).LT.TOL) THEN
            ikto = ik1
            ikti = ik2
          ENDIF
        ENDDO
      ENDDO
 
c...  Build connection map so that the boundary rings
c     are assigned (zero volume) polygons:
c      CALL BuildMap 

c     CUTPT1
c     CUTPT2
c     MAXKPTS
c     MAXRINGS
c     CUTRING
c
      nopriv = .FALSE.

      cutpt1 = ikto
      cutpt2 = ikti
      cutring = irsep - 1
      maxkpts = nks(irsep)
      maxrings = irwall
      indexiradj = 1

c...  For consistency with original code:
      ir = maxrings
      ik = maxkpts

      DEALLOCATE(rvp)
      DEALLOCATE(zvp)
      DEALLOCATE(rcen)
      DEALLOCATE(zcen)
      DEALLOCATE(bval)

      CALL OutputData(85,'Working on it')
c      CALL SaveSolution
c      STOP 'WORKING ON IT'

      WRITE(0,*) '********************'
      WRITE(0,*) 'SETTING IRSEP2=IRSEP'
      WRITE(0,*) '********************'
      irsep2 = irsep

      RETURN
 9000 format(3x,'Element',i5,' = (',i3,',',i3,'): (',
     >       e17.10e2,',',e17.10e2,')',6x,'(',e17.10e2,',',e17.10e2,')')
 9001 format(3x,'Field ratio  = ',e17.10e2,13x,
     >       '(',e17.10e2,',',e17.10e2,')')
 9002 format(3x,26x,'(',e17.10e2,',',e17.10e2,')',6x,
     >       '(',e17.10,',',e17.10,')')
99    STOP
      END
c
c
c
      SUBROUTINE PrepQuasiDoubleNull
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'
      INCLUDE 'pindata'

      INTEGER ik,ir,id
c      INTEGER ik,ir,id,ikmp2(MAXNRS),ikmp1(MAXNRS)
      REAL    dist


c KILL    
c...  Map PSIN data to correct rings:
      
c...  Find inner midplane and outer midplane cell indecies.  If inner
c     midplane index does not exist (since the rings is only on the outside) 
c     then use the low field side data only.  Same for the inside.  It is
c     safe to assume that the grid is structured at this point, since no
c     grid modifications have been done yet (this could change in the future
c     if the grid generators become more flexible, but that likely won't be
c     for a while):
      
c...  1 cm accuracy arbitrary here - fix:

c      CALL outputdata(85,'1')

c      IF (cgridopt.EQ.0) 
c     .  CALL ER('PrepQuasiDoubleNull',
c     .          'Check that the following logic works for JET grids',
c     .          *99)


c      DO ir = irsep, irwall-1
c        ikmp2(ir) = 0
c        ikmp1(ir) = 0
c        DO ik = 1, nks(ir)       
c          id = korpg(ik,ir)
c          IF ((zvertp(1,id).LE.0.0.AND.zvertp(4,id).GT.0.0).OR.
c     .        (zvertp(1,id).GT.0.0.AND.zvertp(4,id).LE.0.0)) THEN
c            IF (rvertp(1,id).LT.r0) THEN
c              IF (ikmp2(ir).EQ.0) THEN
c                ikmp2(ir) = ik
c              ELSE
c                CALL ER('PrepQuasiDoubleNull',
c     .                  'Multiple low IK midplane cells found',*98)
c              ENDIF
c            ENDIF
c            IF (rvertp(1,id).GT.r0) THEN
c              IF (ikmp1(ir).EQ.0) THEN
c                ikmp1(ir) = ik
c              ELSE
c                CALL ER('PrepQuasiDoubleNull',
c     .                  'Multiple high IK midplane cells found',*98)
c              ENDIF
c            ENDIF
c          ENDIF
c        ENDDO
c...    This will trigger when rings in the upper part of the MC 
c       are included -- development required:
c        IF (ikmp2(ir).EQ.0.AND.ikmp1(ir).EQ.0) 
c     .    CALL ER('PrepQuasiDoubleNull','Unable to identify midplane '//
c     .            'cells',*99)
cc...    This is all sort of a strange business because the PSIN values
cc       should be the same for both the inner and outer targets, and 
cc       yet 2 values are stored... not sure why:
c        IF (ikmp2(ir).GT.0) psitarg(ir,2) = psitarg(irorg2(ir),2)
c        IF (ikmp2(ir).EQ.0) psitarg(ir,2) = psitarg(irorg2(ir),1)
c        IF (ikmp1(ir).GT.0) psitarg(ir,1) = psitarg(irorg2(ir),1)
c        IF (ikmp1(ir).EQ.0) psitarg(ir,1) = psitarg(irorg2(ir),2)
c      
cc...    Load PSIN data arrays, which are used later to load PSIN data on the 
cc       modified grid.  Note that the fact that the UEDGE double null grids have 
cc       cell boundaries along the midplane, which makes getting the geometry
cc       data (RHO on the outside) particularly easy, and may not be the case for 
cc       other grids.  This also breaks down for rings that do not pass the inner
cc       or outer midplanes:
c      
c        IF (ikmp2(ir).GT.0) THEN
c          id = korpg(ikmp2(ir),ir)
c          dist = 0.5 * SQRT((rvertp(3,id)-rvertp(4,id))**2+
c     .                      (zvertp(3,id)-zvertp(4,id))**2)
c          psindat(2) = psindat(2) + 1
c          psidat(psindat(2)  ,3) = psidat(psindat(2),3) + dist
c          psidat(psindat(2)+1,3) = psidat(psindat(2),3) + dist
c          psidat(psindat(2)  ,4) = psitarg(ir,2)
c        ENDIF
c        IF (ikmp1(ir).GT.0) THEN
c          id = korpg(ikmp1(ir),ir)
c          dist = 0.5 * SQRT((rvertp(3,id)-rvertp(4,id))**2+
c     .                      (zvertp(3,id)-zvertp(4,id))**2)
c          psindat(1) = psindat(1) + 1
c          psidat(psindat(1)  ,1) = psidat(psindat(1),1) + dist
c          psidat(psindat(1)+1,1) = psidat(psindat(1),1) + dist
c          psidat(psindat(1)  ,2) = psitarg(ir,1)
c        ENDIF
c      ENDDO

c      CALL outputdata(86,'2')
c      STOP 'sdfsdf'
      
c...  Add virtual boundary cells, which will be stripped off later:
      IF (ctargopt.EQ.0.OR.ctargopt.EQ.1.OR.ctargopt.EQ.2.OR.
     .    ctargopt.EQ.3.OR.ctargopt.EQ.6) 
     .  CALL AddPoloidalBoundaryCells      

      RETURN
c98    WRITE(0,*) 'IR=',ik,ir,ikmp2(ir),ikmp1(ir)
      CALL OutputData(85,'Problems in PrepQuasiDoubleNull')
99    STOP
      END
c ======================================================================
c
c    
c
      SUBROUTINE AddPoloidalBoundaryCells
      USE mod_grid_divimp
      IMPLICIT none

      include 'params'
      include 'cgeom'

      INTEGER ik,ir

      DO ir = irsep, nrs

        IF (nks(ir)+2.GT.MAXNKS) 
     .    CALL ER('AddPoloidalBoundaryCells','MAXNKS exceeded',*99)

        DO ik=nks(ir)+1,2,-1
          rs    (ik,ir) = rs    (ik-1,ir)
          zs    (ik,ir) = zs    (ik-1,ir)
          bratio(ik,ir) = bratio(ik-1,ir)
          kbfs  (ik,ir) = kbfs  (ik-1,ir)
          korpg (ik,ir) = korpg (ik-1,ir)
          IF (ALLOCATED(divimp_ik)) THEN
            divimp_ik(ik,ir) = divimp_ik(ik-1,ir)
            divimp_ir(ik,ir) = divimp_ir(ik-1,ir)
          ENDIF
        ENDDO
        bratio(1,ir) = bratio(2,IR)
        kbfs  (1,ir) = kbfs  (2,IR)
        rs    (1,ir) = 0.5 * (RVERTP(1,KORPG(2,IR)) *
     .                        RVERTP(2,KORPG(2,IR)))
        zs    (1,ir) = 0.5 * (ZVERTP(1,KORPG(2,IR)) *
     .                        ZVERTP(2,KORPG(2,IR)))
      
        BRATIO(NKS(IR)+2,IR) =BRATIO(NKS(IR)-1,IR)
        KBFS  (NKS(IR)+2,IR) =KBFS  (NKS(IR)-1,IR)
        RS    (NKS(IR)+2,IR) =0.5 * (RVERTP(3,KORPG(NKS(IR)+1,IR))*
     +                               RVERTP(4,KORPG(NKS(IR)+1,IR)))
        ZS    (NKS(IR)+2,IR) =0.5 * (ZVERTP(3,KORPG(NKS(IR)+1,IR))*
     +                               ZVERTP(4,KORPG(NKS(IR)+1,IR)))
        NKS(IR) = NKS(IR) + 2
      ENDDO
      ikti = ikti + 1
      ikto = ikto + 1

      RETURN
99    STOP
      END



c
c ======================================================================
c
c subroutine: DumpQuadrangles
c
      SUBROUTINE DumpQuadrangles
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'


      REAL       TOL
      PARAMETER (TOL=1.0E-07)

      REAL GetMach,GetJsat,GetFlux 

      TYPE type_quad
        INTEGER :: ik(0:4),ir(0:4),map(4)
      ENDTYPE type_quad

      INTEGER fp,fp2,ik1,ik2,ik3,ik4,ir1,ir2,ir3,ir4,i1,i2,ik,ir,in,
     .        nquad
      LOGICAL status
      REAL    x1,x2,y1,y2,deltax,deltay,Bfrac,fact,
     .        tarte,tarti,tarne,tarv,tarflux,tarisat,tarM
      REAL*8  Bx,By,Bz,beta,brat

      INTEGER, ALLOCATABLE :: quadmap(:,:)

      TYPE(type_quad),ALLOCATABLE :: quad(:)


      fp  = 99
      fp2 = 98

      ALLOCATE(quadmap(MAXNKS,MAXNRS))
      ALLOCATE(quad(MAXNKS*MAXNRS))

      nquad = 0      
      DO ir = 1, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE

        DO ik = 1, nks(ir)
          IF (ir.LT.irsep.AND.ik.EQ.nks(ir)) CYCLE

          IF (ir.LT.irsep) THEN
            ik1 = ik - 1
            ir1 = ir 
            IF (ik.EQ.1        ) ik1 = nks(ir) - 1
            ik2 = ikouts(ik,ir)
            ir2 = irouts(ik,ir) 
            ik3 = ik + 1
            ir3 = ir 
            IF (ik.EQ.nks(ir)-1) ik3 = 1
            ik4 = ikins(ik,ir)
            ir4 = irins(ik,ir)
            IF (ir4.LE.1) THEN
              ik4 = -1
              ir4 = -1
            ENDIF
          ELSE
            ik1 = ik - 1
            ir1 = ir
            IF (ik.EQ.1      ) THEN
              ik1 = -4
              ir1 = -4
            ENDIF
            ik2 = ikouts(ik,ir)
            ir2 = irouts(ik,ir)
            IF (irouts(ik,ir).EQ.irwall) THEN
              ik2 = -3
              ir2 = -3
            ENDIF
            ik3 = ik + 1
            ir3 = ir
            IF (ik.EQ.nks(ir)) THEN
              ik3 = -5
              ir3 = -5 
            ENDIF
            ik4 = ikins(ik,ir)
            ir4 = irins(ik,ir)
            IF (irins(ik,ir).EQ.irtrap) THEN
              ik4 = -2
              ir4 = -2
            ENDIF
          ENDIF

          nquad = nquad + 1
          quad(nquad)%ik(0) = ik
          quad(nquad)%ik(1) = ik1          
          quad(nquad)%ik(2) = ik2
          quad(nquad)%ik(3) = ik3
          quad(nquad)%ik(4) = ik4
          quad(nquad)%ir(0) = ir
          quad(nquad)%ir(1) = ir1
          quad(nquad)%ir(2) = ir2
          quad(nquad)%ir(3) = ir3
          quad(nquad)%ir(4) = ir4

          quadmap(ik,ir) = nquad
        ENDDO
      ENDDO

c...  Convert connection map to 1D indexing from IK,IR indexing:
      DO i1 = 1, nquad
        DO i2 = 1, 4
          IF (quad(i1)%ik(i2).LT.0) THEN
            quad(i1)%map(i2) = quad(i1)%ik(i2)
          ELSE
            quad(i1)%map(i2) = quadmap(quad(i1)%ik(i2),quad(i1)%ir(i2))
          ENDIF
        ENDDO
      ENDDO


c...  Dump triangles:
      OPEN(UNIT=fp,FILE='quadrangles.geometry',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      

c...  Header:
      WRITE(fp,'(A)') '* VERSION 1.0 OSM GEOMETRY FILE FOR QUADRANGLE'//
     .                ' GRID'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A)') '* CELL DATA'
      WRITE(fp,'(A)') '*'

c          WRITE(fp,'(2I4,2X,2F6.2,2X,I3,1X,4(3I4,2X),2X,4(2F10.6,2X:))')

      WRITE(fp,'(A7,2X,2A6,2X,A3,2X,4(2A4,2X),2X,4(2A10,2X:))')
     .  '* INDEX','RS','ZS','NP',
     .  'IN1','IS1','IN2','IS2','IN3','IS3',
     .  'IN4','IS4','RV1','ZV1','RV2','ZV2','RV3','ZV3',
     .  'RV4','ZV4'
     .  
      WRITE(fp,'(A7,2X,2A6,2X,3X,2X,4(10X),2X,4(2A10,2X:))')
     .  '*      ','(m)','(m)','(m)','(m)','(m)','(m)','(m)','(m)',
     .            '(m)','(m)'

      WRITE(fp,*) nquad
      
      DO i1 = 1, nquad
        ik = quad(i1)%ik(0)
        ir = quad(i1)%ir(0)
        in = korpg(ik,ir)

        WRITE(fp,'(I7,2X,2F6.2,2X,I3,2X,4(2I4,2X),2X,4(2F10.6,2X),
     .             2X,2I4)')
     .    i1,
     .    rs(ik,ir),zs(ik,ir),nvertp(in),
     .    quad(i1)%map(1),3,quad(i1)%map(2),4,
     .    quad(i1)%map(3),1,quad(i1)%map(4),2,
     .    (rvertp(i2,in),zvertp(i2,in),i2=1,nvertp(in)),
     .    ik,ir
      ENDDO
      CLOSE(fp)      


      fact = qtim * qtim * emi / crmi

c...  Dump plasma data:
      OPEN(UNIT=fp ,FILE='quadrangles.plasma',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
      OPEN(UNIT=fp2,FILE='quadrangles.efield',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      

c...  Header:
      WRITE(fp,'(A)') '* VERSION 1.0 OSM PLASMA FILE FOR QUADRANGLE '//
     .                'GRID'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A)') '* BULK PLASMA DATA'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A7,2A8,2X,A10,2X,3A10,2X,3A12)') 
     .  '* Index','Te','Ti','ne','vx',
     .  'vy','vz','Bx','By','Bz'
      WRITE(fp,'(A7,2A8,2X,A10,2X,3A10,2X,3A12)') 
     .  '*      ','(eV)','(eV)','(m-3)','(m s-1)',
     .  '(m s-1)','(m s-1)','(Tesla)','(Tesla)','(Tesla)'


c...  Header:
      WRITE(fp2,'(A)') '* VERSION 1.0 OSM PLASMA FILE FOR QUADRANGLE '//
     .                'GRID'
      WRITE(fp2,'(A)') '*'
      WRITE(fp2,'(A)') '* BULK PLASMA DATA (PART DEUX)'
      WRITE(fp2,'(A)') '*'
      WRITE(fp2,'(A7,3A12)') 
     .  '* Index','Ex','Ey','Ez'
      WRITE(fp2,'(A7,3A12)')
     .  '*      ','V m-1','V m-1','V m-1'

      WRITE(fp ,*) nquad
      WRITE(fp2,*) nquad
      
      DO i1 = 1, nquad
        ik = quad(i1)%ik(0)
        ir = quad(i1)%ir(0)
        in = korpg(ik,ir)

        x1 = 0.5 * (rvertp(1,in) + rvertp(2,in))
        y1 = 0.5 * (zvertp(1,in) + zvertp(2,in))
        x2 = 0.5 * (rvertp(3,in) + rvertp(4,in))
        y2 = 0.5 * (zvertp(3,in) + zvertp(4,in))
        deltax = (x2 - x1)
        deltay = (y2 - y1)
        brat = DBLE(bratio(ik,ir))
        beta = DBLE(deltax / deltay)
        Bz = DSQRT(1.0 - brat**2)
        By = brat * DSQRT(1.0D0 / (1.0D0 + beta**2)) * 
     .       DBLE(SIGN(1.0,deltay))
        Bx = beta * By
c...    CBPHI is the on-axis B-field value specified in the OSM input file:
        Bfrac = cbphi * r0 / rs(ik,ir) 

        WRITE(fp,'(I7,2F8.2,2X,1P,E10.2,2X,3E10.2,2X,3E12.4,0P,2X,2I4)')
     .    i1,
     .    ktebs(ik,ir),ktibs(ik,ir),knbs(ik,ir),
     .    SNGL(Bx)*kvhs(ik,ir)/qtim,
     .    SNGL(By)*kvhs(ik,ir)/qtim,
     .    SNGL(Bz)*kvhs(ik,ir)/qtim,
     .    SNGL(Bx)*Bfrac,SNGL(By)*Bfrac,SNGL(Bz)*Bfrac,
     .    ik,ir

c...    Dump efield data:
        WRITE(fp2,'(I7,1P,3E12.4,2X,2I4)') 
     .    i1,
     .    SNGL(Bx)*kes(ik,ir)/fact,
     .    SNGL(By)*kes(ik,ir)/fact,
     .    SNGL(Bz)*kes(ik,ir)/fact,
     .    ik,ir
      ENDDO

c...  All done:
      CLOSE(fp)      
      CLOSE(fp2)      

      DEALLOCATE(quadmap)
      DEALLOCATE(quad)

      RETURN
96    WRITE(0,*) 'DUMPQUADRANGLES: PROBLEMS WITH FILE ACCESS'
      STOP
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
      SUBROUTINE BuildLinearGrid
      IMPLICIT none
      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'


      INTEGER id,ik,ir,i1,grid_option,nrings_inner,nrings_outer
      REAL*8  r,delr,L,r1,r2,z1,z2,frac1,frac2,
     .        vessel_radius,brat,frac,r_inner,r_outer,delta
     .        , delz ! jdemod

      grid_option = 8 ! 3  ! 7  ! 6

      brat = 1.0

      r0 = 0.0001D0  ! Need this tiny displacement to keep EIRENE04 from falling over 
c      r0 = 0.0000001D0  ! Need this tiny displacement to keep EIRENE04 from falling over 

      SELECTCASE (grid_option)
        CASE (1)  ! Full vessel, mirrored
          vessel_radius = 0.02D0
          L = 3.6D0                   ! Total length of mirrored plasma column (m)
          r = 0.015D0                 ! Plasma radius (m)
          z0 = L / 2.0D0              ! Height of the centre of the plasma (m)
          delr = (vessel_radius - r)  ! Distance from plasma to outer wall (m)
          maxrings = 10               ! Number of flux tubes (if changed, also need to change triangle grid in input file)
          nks(1:maxrings) = 100       ! Number of cells on each tube
        CASE (2)  ! Full vessel
          vessel_radius = 0.02D0
          L = 1.8D0
          r = 0.015D0          
          z0 = L / 2.0D0              
          delr = (vessel_radius - r)  
          maxrings = 10      
          nks(1:maxrings) = 100         
        CASE (3) ! Target chamber
          brat = 0.05 ! 0.985 ! 0.5
  
          vessel_radius = 0.02D0 ! 0.05D0 ! 0.02D0
          L = 0.55D0  ! 0.56D0
          r = 0.015D0 ! 0.03D0  ! 0.015D0
          z0 = L / 2.0D0 ! 0.0 ! L / 2.0D0              
          delr = (vessel_radius - r)  
          maxrings = 10      
          nks(1:maxrings) = 20  ! 50  ! 175
        CASE (4) ! Target chamber: fancy
          vessel_radius = 0.16D0 
          L = 0.55D0
          r = 0.08D0          
          z0 = L / 2.0D0              
          delr = (vessel_radius - r)  
          maxrings = 50      
          nks(1:maxrings) = 150
        CASE (5) ! Target chamber: fancy #2, full vessel and small volume 
          vessel_radius = 0.16D0 
          L = 0.55D0
          r = 0.03D0          
          z0 = L / 2.0D0              
          delr = (vessel_radius - r)  
          maxrings = 50      
          nks(1:maxrings) = 150
        CASE (6) ! Target chamber: fancy #3, small volume 
          vessel_radius = 0.05D0 
          L = 0.55D0
          r = 0.03D0          
          z0 = L / 2.0D0              
          delr = (vessel_radius - r)  
          maxrings = 50      
          nks(1:maxrings) = 150
        CASE (7) ! Grid with 2 radial regions
          vessel_radius = 0.05D0 
          L = 0.55D0
          z0 = L / 2.0D0              
          nrings_inner = 30
          nrings_outer = 20
          maxrings = nrings_inner + nrings_outer
          r_inner = 0.01D0
          r_outer = 0.03D0     
          delr = (vessel_radius - r_outer)  
          nks(1:maxrings) = 150
          ! jdemod
        CASE (8)  ! Rectangular grid
          vessel_radius = 0.02D0
          L = 50.0D0
          r = 0.20D0          
          z0 = L / 2.0D0              
          delr = (vessel_radius - r)  
          maxrings = 50      
          nks(1:maxrings) = 100         
      ENDSELECT

      id = 0

      DO ir = 1, maxrings
        IF (.TRUE.) THEN
          IF     (grid_option.EQ.4) THEN
            frac1 = (DBLE(ir-1) / DBLE(maxrings))**0.7
            frac2 = (DBLE(ir  ) / DBLE(maxrings))**0.7
            r1 = frac1 * r
            r2 = frac2 * r
          ELSEIF (grid_option.EQ.7) THEN
            IF (ir.LE.nrings_inner) THEN
              frac  = DBLE(ir-1) / DBLE(nrings_inner)
              delta = r_inner / DBLE(nrings_inner)
              r1 = frac * r_inner
              r2 = r1 + delta
            ELSE
              frac  = DBLE(ir-nrings_inner-1) / DBLE(nrings_outer)
              delta = (r_outer - r_inner) / DBLE(nrings_outer)
              r1 = r_inner + frac * (r_outer - r_inner)
              r2 = r1 + delta
            ENDIF
          ELSEif (grid_option.eq.8) then 
            ! jdemod
            frac  = DBLE(ir-1) / DBLE(maxrings)
            delta = r / DBLE(maxrings)
            r1 = frac * r 
            r2 = r1 + delta
          ELSE
            frac  = DBLE(ir-1) / DBLE(maxrings)
            delta = r / DBLE(maxrings)
            r1 = frac * r 
            r2 = r1 + delta
          ENDIF
          IF (ir.EQ.1) r1 = r1 + r0
        ENDIF
c
c        delz = L/nks(ir)
c
C       Krieger IPP/07 - SUN compiler does not know SNGL, replaced by REAL  -strange since SNGL is used elsewhere... -SL
C       psitarg(ir,1) = ABS(0.5*(SNGL(r1+r2)))
C       psitarg(ir,2) = ABS(0.5*(SNGL(r1+r2)))
        psitarg(ir,1) = ABS(0.5*(REAL(r1+r2)))
        psitarg(ir,2) = ABS(0.5*(REAL(r1+r2)))
        idring(ir) = TARTOTAR

c        WRITE(0,*) 'IR:',ir,psitarg(ir,1)

c       nks(ir) = 100
c
        DO ik = 1, nks(ir)

          SELECTCASE (grid_option)
            CASE (1)  ! Full vessel, mirrored
              IF (.TRUE.) THEN
                frac = DBLE(ik-1) / DBLE(nks(ir)) 
                delta = L / DBLE(nks(ir)) 
                z1 = (0.5 - frac) * L
                z2 = z1 - delta
              ENDIF
            CASE (2)  ! Full vessel
              IF (.TRUE.) THEN
                frac = DBLE(ik-1) / DBLE(nks(ir)) 
                delta = L / DBLE(nks(ir)) 
                z1 = (1.0 - frac) * L
                z2 = z1 - delta       
              ENDIF
            CASE (3)  ! Target chamber
              IF (.TRUE.) THEN
                frac = DBLE(ik-1) / DBLE(nks(ir)) 
                delta = L / DBLE(nks(ir)) 
                z1 = (0.5 - frac) * L + z0
c                z1 = (1.0 - frac) * L 
c
                z2 = z1 - delta       
c
              ENDIF
            CASE (4:7) ! Target chamber: fancy
              frac1 = DBLE(ik-1) / DBLE(nks(ir)) 
              frac2 = DBLE(ik  ) / DBLE(nks(ir)) 
              frac1 = SIGN(.5D0,frac1-.5D0)*(ABS(frac1-.5)/0.5)**1.0+0.5
              frac2 = SIGN(.5D0,frac2-.5D0)*(ABS(frac2-.5)/0.5)**1.0+0.5
              z1 = (1.0 - frac1) * L
              z2 = (1.0 - frac2) * L     
            CASE (8)  
               ! jdemod
              IF (.TRUE.) THEN
                frac = DBLE(ik-1) / DBLE(nks(ir)) 
                delta = L / DBLE(nks(ir)) 
                z1 = (1.0 - frac) * L
                z2 = z1 - delta       
              ENDIF
          ENDSELECT

c          frac = ((ABS(0.5 * (z1 + z2) - z0) + 0.001) / L * 2.0)**0.05
c          IF (ir.EQ.2) WRITE(0,*) frac
c          bratio(ik,ir) = SNGL(brat * frac)
          bratio(ik,ir) = SNGL(brat)
          kbfs  (ik,ir) = 1.0 / brat
          bts   (ik,ir) = cbphi 

          id = id + 1

          korpg(ik,ir) = id

          nvertp(id) = 4

!          frac = 1.0D0 + 1.0D0 * DBLE(ik-1) / DBLE(nks(ir) - 1)
          frac = 1.0D0

          IF (ik.EQ.1) THEN
            rvertp(1,id) = SNGL(r1)
            rvertp(2,id) = SNGL(r2)
            zvertp(1,id) = SNGL(z1)
            zvertp(2,id) = SNGL(z1)
          ELSE
            rvertp(1,id) = rvertp(4,id-1)
            rvertp(2,id) = rvertp(3,id-1)
            zvertp(1,id) = zvertp(4,id-1)
            zvertp(2,id) = zvertp(3,id-1)
          ENDIF

          rvertp(3,id) = SNGL(r2 * frac)
          rvertp(4,id) = SNGL(r1 * frac)
          zvertp(3,id) = SNGL(z2)
          zvertp(4,id) = SNGL(z2)

          rs(ik,ir) = 0.0
          zs(ik,ir) = 0.0
          DO i1 = 1, nvertp(id)
            rs(ik,ir) = rs(ik,ir) + rvertp(1,id)
            zs(ik,ir) = zs(ik,ir) + zvertp(1,id)
          ENDDO
          rs(ik,ir) = rs(ik,ir) / REAL(nvertp(id))
          zs(ik,ir) = zs(ik,ir) / REAL(nvertp(id))

        ENDDO

      ENDDO

      npolyp  = id
      vpolmin = (MAXNKS*MAXNRS - npolyp) / 2 + npolyp
      vpolyp  = vpolmin

      ikto = 2
      ikti = 3

      rves = 0
      rvesm = 0

      irsep  = 1
      irwall = maxrings 
      irtrap = irwall
      nrs    = irwall
      nbr    = 0

c      WRITE(0,*) 'NVERT:',nvertp(5)

      CALL InsertRing(1         ,BEFORE,PERMANENT)
      CALL InsertRing(maxrings+1,AFTER ,PERMANENT)

c      WRITE(0,*) 'NVERT:',nvertp(5)

c...  Necessary..? 
      cutring = 1
      cutpt1 = ikto
      cutpt2 = ikti

      idring(1) = -1
      idring(nrs) = -1

c...  Modify the grid based on entries in the GRDMOD array assigned 
c     from the input file:
c      IF (grdnmod.NE.0) CALL TailorGrid

      rmin = HI
      rmax = LO
      zmin = HI
      zmax = LO
      DO ir = 1, nrs
        DO ik = 1, nks(ir)
          rmin = MIN(rmin,rs(ik,ir))
          rmax = MAX(rmax,rs(ik,ir))
          zmin = MIN(zmin,zs(ik,ir))
          zmax = MAX(zmax,zs(ik,ir))
        ENDDO
      ENDDO

      rxp =  0.0 
      zxp =  0.25 * zmin + 0.75 * zmax
c      zxp = -0.25 * L

c...  Neutral wall
      IF (cneur.EQ.4) THEN

         SELECTCASE (grid_option)
           CASE (1)  ! Full vessel, mirrored
             nves = 20
             ir = irwall-1
             r1 = DBLE(rvertp(2,korpg(1      ,ir)))
             r2 = r1 + delr
             z1 = DBLE(zvertp(2,korpg(1      ,ir)))
             z2 = DBLE(zvertp(3,korpg(nks(ir),ir)))
  
             rves(1)  =  SNGL(r1)
             zves(1)  =  SNGL(z1)
             rves(2)  =  SNGL(r2)
             zves(2)  =  SNGL(z1)
  
             rves(3)  =  SNGL(r2)
             zves(3)  =  SNGL(L) / 2.0 - 0.56
             rves(4)  =  SNGL(r1) + 0.0101 
             zves(4)  =  SNGL(L) / 2.0 - 0.56
             rves(5)  =  SNGL(r1) + 0.0001
             zves(5)  =  SNGL(L) / 2.0 - 0.57
             rves(6)  =  SNGL(r2)
             zves(6)  =  SNGL(L) / 2.0 - 0.57
  
             rves(7)  =  SNGL(r2)
             zves(7)  =  0.03
             rves(8)  =  SNGL(r1) + 0.0001
             zves(8)  =  0.03
             rves(9)  =  SNGL(r1) + 0.0001
             zves(9)  =  0.02
             rves(10) =  SNGL(r2)
             zves(10) =  0.02
  
             rves(11) =  r2
             zves(11) = -0.02
             rves(12) =  r1 + 0.0001
             zves(12) = -0.02
             rves(13) =  r1 + 0.0001
             zves(13) = -0.03
             rves(14) =  r2
             zves(14) = -0.03
  
             rves(15) =  r2
             zves(15) = -L / 2.0 + 0.56
             rves(16) =  r1 + 0.0001
             zves(16) = -L / 2.0 + 0.56
             rves(17) =  r1 + 0.0101
             zves(17) = -L / 2.0 + 0.57
             rves(18) =  r2
             zves(18) = -L / 2.0 + 0.57
  
             rves(19) =  r2
             zves(19) =  z2
             rves(20) =  r1
             zves(20) =  z2
           CASE (2)  ! Full vessel
             nves = 12
             ir = irwall-1
             r1 = rvertp(2,korpg(1      ,ir))
             r2 = r1 + delr
             z1 = zvertp(2,korpg(1      ,ir))
             z2 = zvertp(3,korpg(nks(ir),ir))
  
             rves(1)  =  r1
             zves(1)  =  z1
             rves(2)  =  r2
             zves(2)  =  z1
  
             rves(3)  =  r2
             zves(3)  =  L - 0.56
             rves(4)  =  r1 + 0.0101 
             zves(4)  =  L - 0.56
             rves(5)  =  r1 + 0.0001 
             zves(5)  =  L - 0.57
             rves(6)  =  r2
             zves(6)  =  L - 0.57
  
             rves(7)  =  r2
             zves(7)  =  0.03
             rves(8)  =  r1 + 0.0001
             zves(8)  =  0.03
             rves(9)  =  r1 + 0.0001
             zves(9)  =  0.02
             rves(10) =  r2
             zves(10) =  0.02
  
             rves(11) =  r2
             zves(11) =  z2
             rves(12) =  r1
             zves(12) =  z2
           CASE (3)  ! Target chamber
             nves = 7
             ir = irwall-1
             r1 = rvertp(2,korpg(1      ,ir)) - 0.0001 ! So that the clipping code is required / activated
             r2 = r1 + delr
             z1 = zvertp(2,korpg(1      ,ir))
             z2 = zvertp(3,korpg(nks(ir),ir)) 
  
             rves(1) =  r1
             zves(1) =  z1
             rves(2) =  r2
             zves(2) =  z1

             rves(3) =  r2
             zves(3) =  0.55 * z1 + 0.45 * z2
             rves(4) =  r2
             zves(4) =  0.50 * z1 + 0.50 * z2
             rves(5) =  r2
             zves(5) =  0.45 * z1 + 0.55 * z2
  
             rves(6) =  r2
             zves(6) =  z2 
             rves(7) =  rvertp(3,korpg(nks(ir),ir)) - 0.0001 ! r1
             zves(7) =  z2
           CASE (4:5)  ! Target chamber: fancy
             nves = 10
             ir = irwall-1
             r1 = rvertp(2,korpg(1      ,ir))
             r2 = r1 + delr
             z1 = zvertp(2,korpg(1      ,ir))
             z2 = zvertp(3,korpg(nks(ir),ir))
  
             rves(1) =  r1
             zves(1) =  z1
             rves(2) =  r2
             zves(2) =  z1

             rves(3) =  r2
             zves(3) =  0.11

             rves(4) =  0.21
             zves(4) =  0.11

             rves(5) =  0.21
             zves(5) =  0.06

             rves(6) =  0.212
             zves(6) =  0.06

             rves(7) =  0.212
             zves(7) =  0.05

             rves(8) =  0.21
             zves(8) =  0.05

             rves(9) =  0.21
             zves(9) =  z2
  
             rves(10) =  r1
             zves(10) =  z2
           CASE (6:7)  ! Target chamber: fancy #3, small volume
             nves = 10
             ir = irwall-1
             r1 = rvertp(2,korpg(1      ,ir))
             r2 = r1 + delr
             z1 = zvertp(2,korpg(1      ,ir))
             z2 = zvertp(3,korpg(nks(ir),ir))
  
             rves(1) =  r1
             zves(1) =  z1
             rves(2) =  r2
             zves(2) =  z1

             rves(3) =  r2
             zves(3) =  0.11

             rves(4) =  r2 + 0.01
             zves(4) =  0.11

             rves(5) =  r2 + 0.01
             zves(5) =  0.06

             rves(6) =  r2 + 0.015
             zves(6) =  0.06

             rves(7) =  r2 + 0.015
             zves(7) =  0.05

             rves(8) =  r2 + 0.01
             zves(8) =  0.05

             rves(9) =  r2 + 0.01
             zves(9) =  z2
  
             rves(10) =  r1
             zves(10) =  z2
        ENDSELECT
      ENDIF


      IF (.TRUE.) THEN
        nvesm = nves - 1
        DO i1 = 1, nves-1
          rvesm(i1,1) = rves(i1)
          zvesm(i1,1) = zves(i1)
          rvesm(i1,2) = rves(i1+1)
          zvesm(i1,2) = zves(i1+1)
        ENDDO
      ENDIF
 
c      CALL DumpGrid('BUILDING LINEAR GRID')

      IF (grdnmod.GT.0) CALL TailorGrid


      CALL OutputData(85,'Linear')


c...  Add virtual boundary cells, which will be stripped off later:
      IF (CTARGOPT.EQ.0.OR.CTARGOPT.EQ.1.OR.CTARGOPT.EQ.2.OR.
     .    CTARGOPT.EQ.3.OR.CTARGOPT.EQ.6) 
     .   CALL AddPoloidalBoundaryCells

c      STOP 'WHA-WHO!'

      RETURN
 99   STOP
      END


c
c ======================================================================
c
c     
c     
c     -------------------------------------------------------------------
c     
c     The following code puts together the ITER ribbon grid. It 
c     uses the stand alone grid generator which uses CASTEM data
c     and calls it in-line.
c     
c     
      SUBROUTINE BuildRibbonGrid
      use ribbon_grid_options
      use error_handling
      use castem_field_line_data
      IMPLICIT none
      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'


      character*512 :: ident_file,intersection_file
      integer :: ierr
      integer :: in,ik,ir,it,is
      character*512 :: cmd,source_dir

      ! set grid run descriptor
      crun = 'ITER FIRST WALL RIBBON GRID'

      write(0,*) 'Building ribbon grid:'
      write(0,*) 'RG_CASTEM_DATA:',trim(rg_castem_data),':'

      ierr = 0

      ! get path to the data directory
      
      call get_div_data_dir(source_dir,ierr)
      
      if (ierr.ne.0) then 
         ! error getting data directory
         call errmsg('BuildRibbonGrid',
     >            'Error obtaining data directory from environment')
         stop 'Build Ribbon Grid 1'
      endif

      write(0,*) 

      if (ribbon_input_format_opt.eq.0) then 
         ! CASTEM formatted input file

         !ident_file = 'DATA_IDENTIFIER_260410.txt'
         ident_file = 'DATA_IDENTIFIER_'//trim(rg_castem_data)//'.txt'

         !intersection_file = 'DATA_RHO_S_260410.txt'
         intersection_file = 'DATA_RHO_S_'//trim(rg_castem_data)//'.txt'

      
         cmd = 'cp '//trim(source_dir)//'/'//trim(ident_file)//' .'
         call run_system_command(cmd,ierr)
 
         if (ierr.ne.0) then 
            ! error copying ident file
            call errmsg('BuildRibbonGrid',
     >               'Error copying ident file cmd='//trim(cmd))
            stop 'Build Ribbon Grid 2a'
         endif

         cmd = 'cp '//trim(source_dir)//'/'//
     >                trim(intersection_file)//' .'
         call run_system_command(cmd,ierr)

         if (ierr.ne.0) then 
            ! error copying ident file
            call errmsg('BuildRibbonGrid',
     >               'Error copying intersection file cmd='//trim(cmd))
            stop 'Build Ribbon Grid 3'
         endif


         call read_identifier_data(ident_file,ierr)

         if (ierr.ne.0) then 
            call errmsg('Error reading IDENTIFIER data',ierr)
            return
         endif


         call read_castem_intersection_data(intersection_file,ierr)

         if (ierr.ne.0) then 
            call errmsg('Error reading INTERSECTION data',ierr)
            return
         endif

         call calculate_castem_limiter_surface


      elseif (ribbon_input_format_opt.eq.1) then
         ! RAY formatted input file

         !intersection_file = 'DATA_RHO_S_260410.txt'
         intersection_file = trim(rg_castem_data)

      
         cmd = 'cp '//trim(source_dir)//'/'//
     >                trim(intersection_file)//' .'

c         write(0,*) 'cmd:',trim(cmd)

         call run_system_command(cmd,ierr)
 
         if (ierr.ne.0) then 
            ! error copying ident file
            call errmsg('BuildRibbonGrid',
     >               'Error copying ident file cmd='//trim(cmd))
            stop 'Build Ribbon Grid 2b'
         endif

         call read_ray_intersection_data(intersection_file,ierr)

         if (ierr.ne.0) then 
            call errmsg('Error reading INTERSECTION data',ierr)
            return
         endif

         call calculate_ray_limiter_surface

      endif



      call print_field_line_summary

      call generate_grid

      call write_grid

      call assign_grid_to_divimp(maxnrs,maxnks,mves,nrs,nks,
     >     nves,rves,zves,
     >     npolyp,korpg,
     >     nvertp,rvertp,zvertp,
     >     rs,zs)

      call deallocate_castem_storage

      write(0,*) 'Completed ribbon grid generation:'

      ! Assign Bratio = 1.0 to start - may need a more appropriate value. 
      ! Assign these arrays constant values for now
      bratio = 1.0
      kbfs = 1.0 / bratio
      bts  = cbphi
 
      ! Assign PSITARG
      do ir = 1,nrs
                                ! assign R value to PSITARG for these cases
                                ! target 1 is at the nks(ir) end of the ring
                                ! target 2 is at the ik=1 end of the ring
         ik = nks(ir)
         in = korpg(ik,ir)
                                ! target is between vertices 3,4 at the UP end of the ring
         psitarg(ir,1) = (rvertp(3,in) + rvertp(4,in))/2.0

         ik = 1
                                ! target is between vertices 1,2 at the DOWN end of the ring
         psitarg(ir,1) = (rvertp(1,in) + rvertp(2,in))/2.0


         ! assign IDRING as TARTOTAR for all rings to start since this is true for a ribbon grid
         idring(ir) = TARTOTAR

      end do


      vpolmin = (MAXNKS*MAXNRS - npolyp) / 2 + npolyp
      vpolyp  = vpolmin

      ikto = 0
      ikti = maxnks

      irsep = 1
      irwall = nrs
      irtrap = nrs
      nbr = 0

c
c     Flag all cells are non-orthogonal
c
      tagdv = 1.0
c

      ! insert zero volume boundary rings? - boundary rings would need be added to every PFZ section - not really feasible?
      

c     WRITE(0,*) 'NVERT:',nvertp(5)
c      CALL InsertRing(1         ,BEFORE,PERMANENT)
c      CALL InsertRing(maxrings+1,AFTER ,PERMANENT)
c     WRITE(0,*) 'NVERT:',nvertp(5)

                                ! these should not be needed
      cutring = 1
      cutpt1 = 0
      cutpt2 = 0

c...  Modify the grid based on entries in the GRDMOD array assigned 
c     from the input file:
c     IF (grdnmod.NE.0) CALL TailorGrid

      ! these should be based on polygons not cell centers
      rmin = HI
      rmax = LO
      zmin = HI
      zmax = LO
      DO ir = 1, nrs
         DO ik = 1, nks(ir)
            rmin = MIN(rmin,rs(ik,ir))
            rmax = MAX(rmax,rs(ik,ir))
            zmin = MIN(zmin,zs(ik,ir))
            zmax = MAX(zmax,zs(ik,ir))
         ENDDO
      ENDDO

      ! no Xpoint or limiter tip
      rxp = 0.0
      zxp = 0.0
      
      r0 = 0.0
      z0 = 0.0

      
         nvesm = nves - 1
         DO in = 1, nves-1
            rvesm(in,1) = rves(in)
            zvesm(in,1) = zves(in)
            rvesm(in,2) = rves(in+1)
            zvesm(in,2) = zves(in+1)
         ENDDO

! possibly need to add poloidal boundary cells at end of rings - perhaps avoid for now by choosing appropriate target option


c...  Add virtual boundary cells, which will be stripped off later:
c      IF (CTARGOPT.EQ.0.OR.CTARGOPT.EQ.1.OR.CTARGOPT.EQ.2.OR.
c     .     CTARGOPT.EQ.3.OR.CTARGOPT.EQ.6) 
c     .     CALL AddPoloidalBoundaryCells



      RETURN
      END
c
c ======================================================================
c
