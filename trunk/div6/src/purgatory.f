c
c ======================================================================
c
c subroutine: WriteInputFile
c
c

      SUBROUTINE WriteInputFile_97

      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER          ik,ik1,ik2,ir,i1,i2,fp1,fp2,nstsi,walln
      INTEGER          mode
      REAL             material(4)
      DOUBLE PRECISION wallr(MAXPTS,2),wallz(MAXPTS,2)
      CHARACTER        buffer*200

      DATA material / 9642., 1206., 18474., 904./
c
c     Check whether DIVIMP input option requests EIRENE data file:
c
      IF (eirdata.NE.1) RETURN
c
c     Initialization:
c

c
c     Set surface temperatures (in eV):
c
      eirtemp1 = ctargt * 1.38E-23 / ECH
      eirtemp2 = cwallt * 1.38E-23 / ECH


      walln = 0

      DO ir = irsep, nrs
        virloc(ir,IKLO) = 1
        virloc(ir,IKHI) = nks(ir)
      ENDDO


      WRITE(SLOUT,*)

      DO i1 = 1, wallpts

        WRITE(SLOUT,'(2F15.7,F10.3)')
     .    wallpt2(i1,1),wallpt2(i1,2),wallpt(i1,18)

      ENDDO

      WRITE(SLOUT,*)

      DO i1 = 1, wallpts

        WRITE(SLOUT,'(2F15.7,F10.3)')
     .    wallpt(i1,20),wallpt(i1,21),wallpt(i1,18)

      ENDDO


      WRITE(SLOUT,*)

      DO i1 = 1, pcnt

        WRITE(SLOUT,'(2F15.7)')
     .    rw(i1),zw(i1)

        IF (rw(i1).EQ.rw(i1+1).AND.zw(i1).EQ.zw(i1+1))
     .    WRITE(SLOUT,*) 'DUP'

      ENDDO

      fp1 = EIRIN
      fp2 = EIROUT

      OPEN(UNIT=fp1,FORM='FORMATTED',ERR=95,STATUS='OLD')
      OPEN(UNIT=fp2,FORM='FORMATTED',ERR=95,STATUS='REPLACE')

      nvesm = wallpts

      DO i1 = wallpts, 1, -1
        rvesm(i1,1) = wallpt(i1,20)
        zvesm(i1,1) = wallpt(i1,21)
c        rvesm(i1,1) = wallpt2(i1,1)
c        zvesm(i1,1) = wallpt2(i1,2)

        IF (i1.EQ.wallpts) THEN
          i2 = 1
        ELSE
          i2 = i1 + 1
        ENDIF

        rvesm(i1,2) = wallpt(i2,20)
        zvesm(i1,2) = wallpt(i2,21)
c        rvesm(i1,2) = wallpt2(i2,1)
c        zvesm(i1,2) = wallpt2(i2,2)

        IF (i1.EQ.wallpts) THEN
          i2 = 1
        ELSE
          i2 = i1 + 1
        ENDIF

        IF (wallpt(i1,18).EQ.0) THEN
          walln = walln + 1

          wallr(walln,1) = wallpt(i2,20)
          wallz(walln,1) = wallpt(i2,21)
          wallr(walln,2) = wallpt(i1,20)
          wallz(walln,2) = wallpt(i1,21)
c          wallr(walln,1) = wallpt2(i2,1)
c          wallz(walln,1) = wallpt2(i2,2)
c          wallr(walln,2) = wallpt2(i1,1)
c          wallz(walln,2) = wallpt2(i1,2)
        ENDIF
      ENDDO
c
c     This is a somewhat convoluted loop at the moment, which reads
c     through an existing EIRENE input data file that serves as a
c     template for the new data file being written.  Grid and neutral
c     wall specific data are substituted into the template, as well
c     as any EIRENE options/settings that are specifiable
c     from DIVIMP (such as EIRENE execution time):
c
10    CONTINUE

      CALL ReadLine(fp1,buffer,1,*50,*98)

20    CONTINUE

      IF (buffer(1:6).EQ.'*** 0.') THEN


c Need to remove the requirement that the template file have an intitial
c seciton labelled *** 0...

c
c       This section has been added to EIRENE and contains options
c       for the new EIRENE code that are related to the
c       generalization of the grid:
c
        WRITE(fp2,'(A)') '*** 0. DIVIMP RELATED SETUP DATA (DIVIMP)'
        WRITE(fp2,'(2A,I6)') '''Geometry option  (GEOMOPT)  ',
     .    '0-standard   1-from DIVIMP    ''',eirgeom
        WRITE(fp2,'(2A,I6)') '''Grid option      (GRIDOPT)  ',
     .    '0-structured 1-generalized    ''',eirgrid
c        WRITE(fp2,'(2A,I6)') '''AddUsr option    (ADDOPT)   ',
c     .    '0-execute    1-do not execute ''',eiradd
        WRITE(fp2,'(2A,I6)') '''Wall data option (NEUTOPT)  ',
     .    '0-standard   1-accurate       ''',eirneut
        WRITE(fp2,'(2A,I6)') '''Debug option     (DEBUGOPT) ',
     .    '0-off                         ''',eirdebug
c
c       Advance the template file until the next section tag
c       is found:
c
21      CALL ReadLine(fp1,buffer,1,*97,*98)
        IF (buffer(1:3).NE.'***') GOTO 21
        BACKSPACE fp1

      ELSEIF (buffer(1:6).EQ.'*** 1.') THEN
c
c       Set EIRENE execution time:
c
        WRITE(fp2,'(A)') '*** 1. DATA FOR OPERATING MODE (DIVIMP)'
        CALL UpdateLine1I(fp1,fp2,buffer,3,eirtime)
        CALL TransferLine(fp1,fp2,buffer,1)

      ELSEIF (buffer(1:6).EQ.'*** 2.') THEN
c
c       Specify grid size (number of rings and cells per ring):
c
        WRITE(fp2,'(A)') '*** 2. DATA FOR STANDARD MESH (DIVIMP)'

23      CALL TransferLine(fp1,fp2,buffer,1)
        IF (buffer(1:2).EQ.'* ') GOTO 23

        CALL TransferLine(fp1,fp2,buffer,2)
        CALL UpdateLine2I(fp1,fp2,buffer,1,3,irwall-1,nks(irsep)+3)
        CALL TransferLine(fp1,fp2,buffer,3)
        CALL UpdateLine1I(fp1,fp2,buffer,1,nks(irsep)+3)
        CALL TransferLine(fp1,fp2,buffer,10)

      ELSEIF (buffer(1:6).EQ.'*** 3a') THEN
c
c       Configure standard surfaces in the grid, such as targets,
c       cut points and core/wall/pfz boundaries.  The settings
c       for neutral production at the targets are currently
c       hard coded and so are not copied from the template file:
c
        WRITE(fp2,'(A)')
     .    '*** 3a. DATA FOR NON DEFAULT STANDARD SURFACES (DIVIMP)'

        IF (nbr.GT.0) THEN
c         This will have to be improved for more general wall
c         contact which has more than one broken region...
          WRITE(fp2,'(I6)') 15
        ELSE
          WRITE(fp2,'(I6)') 14
        ENDIF

        nstsi = 1
        WRITE(fp2,'(A,I2,A)')
     .    '* RADIAL FLUXSURFACE ',1,' ABSORBING (PLASMA CORE)'
        WRITE(fp2,'(5I6)')
     .    nstsi,1,1,ikto+1,ikti+2
        WRITE(fp2,'(A)')
     .    '     2    -2     0     0     0     6     0     0'

        nstsi = nstsi + 1
        WRITE(fp2,'(A,I2,A)')
     .    '* RADIAL FLUXSURFACE ',1,' TRANSPARENT (PFZ, LEFT TARGET)'
        WRITE(fp2,'(5I6)')
     .    nstsi,1,1,virloc(irtrap+1,IKLO),ikto+1
        WRITE(fp2,'(A)')
     .    '    -1     0010220     0     0     6     0001001'

        nstsi = nstsi + 1
        WRITE(fp2,'(A,I2,A)')
     .    '* RADIAL FLUXSURFACE ',1,' TRANSPARENT (PFZ, RIGHT TARGET)'
        WRITE(fp2,'(5I6)')
     .    nstsi,1,1,ikti+2,virloc(irtrap+1,IKHI)+(ikti-ikto)+2
        WRITE(fp2,'(A)')
     .    '    -1     0010220     0     0     6     0001001'

        IF (nbr.GT.0) THEN
c
c         Partial radial surfaces due to broken rings:
c
          DO ir = irbreak-1, irwall-2
            ik1 = 0
            ik2 = 0
            DO ik = virloc(ir,IKLO), virloc(ir,IKHI)
              IF (irouts(ik,ir).EQ.irwall.AND.ik1.EQ.0) ik1 = ik

              IF (irouts(ik,ir).NE.irwall.AND.ik1.NE.0) THEN
                ik2 = ik - 1

                IF (ik1.GT.ikto) ik1 = ik1 + 1
                IF (ik1.GT.ikti) ik1 = ik1 + 1

                IF (ik2.GT.ikto) ik2 = ik2 + 1
                IF (ik2.GT.ikti) ik2 = ik2 + 1

                nstsi = nstsi + 1
                WRITE(fp2,'(A,I2,A)')
     .            '* RADIAL FLUXSURFACE ',ir,' TRANSPARENT (WALL)'
                WRITE(fp2,'(5I6)')
     .            nstsi,1,ir,ik1,ik2+1
                WRITE(fp2,'(A)')
     .            '    -1     0020110     0     0     6     0001001'

                ik1 = 0
                ik2 = 0
              ENDIF
            ENDDO
          ENDDO
        ENDIF

        nstsi = nstsi + 1
        WRITE(fp2,'(A,I2,A)')
     .    '* RADIAL FLUXSURFACE ',irwall-1,', TRANSPARENT (OUTER WALL)'
        WRITE(fp2,'(5I6)')
     .    nstsi,1,irwall-1,virloc(irwall-1,IKLO),virloc(irwall-1,IKHI)+3
        WRITE(fp2,'(A)')
     .    '    -1     0020110     0     0     6     0001001'

        nstsi = nstsi + 1
        WRITE(fp2,'(A,I2,2A)')
     .    '* POLOIDAL FLUXSURFACE ',1,', SOURCE/REFLECTING ',
     .    '(LEFT TARGET)'
        WRITE(fp2,'(5I6)')
     .    nstsi,2,1,1,irwall-1
        WRITE(fp2,'(6I6)') 1,0,0,0,0,6
        WRITE(fp2,'(2I6)') 1,2

c Replace with an array...
        WRITE(fp2,'(1P,2E12.4)') material(eirmat1),eirtemp1

        WRITE(fp2,'(2A)')
     .    '  1.0000E 00  1.0000E 00  1.0000E 00  1.0000E 00',
     .    '  0.5000E 00  1.0000E 00'

        nstsi = nstsi + 1
        WRITE(fp2,'(A,I2,A)')
     .    '* POLOIDAL FLUXSURFACE ',nks(irsep)+3,
     .    ', SOURCE/REFLECTING (RIGHT TARGET)'
        WRITE(fp2,'(5I6)')
     .    nstsi,2,nks(irsep)+3,1,irwall-1
        WRITE(fp2,'(6I6)') 1,0,0,0,0,6
        WRITE(fp2,'(2I6)') 1,2

c Replace with an array...
        WRITE(fp2,'(1P,2E12.4)') material(eirmat1),eirtemp1

        WRITE(fp2,'(2A)')
     .    '  1.0000E 00  1.0000E 00  1.0000E 00  1.0000E 00',
     .    '  0.5000E 00  1.0000E 00'

        nstsi = nstsi + 1
        WRITE(fp2,'(A,I2,A)')
     .    '* POLOIDAL FLUXSURFACE ',ikto+1,', TRANSPARENT LEFT PFZ/CORE'
        WRITE(fp2,'(5I6)') nstsi,2,ikto+1,1,irsep-1
        WRITE(fp2,'(6I6)') -3,0,0,ikto+2,0,3

        nstsi = nstsi + 1
        WRITE(fp2,'(A,I2,A)')
     .    '* POLOIDAL FLUXSURFACE ',ikto+2,', TRANSPARENT LEFT PFZ/CORE'
        WRITE(fp2,'(5I6)') nstsi,2,ikto+2,1,irsep-1
        WRITE(fp2,'(6I6)') -3,0,0,ikti+1,0,3

        nstsi = nstsi + 1
        WRITE(fp2,'(A,I2,A)')
     .  '* POLOIDAL FLUXSURFACE ',ikti+1,', TRANSPARENT RIGHT PFZ/CORE'
        WRITE(fp2,'(5I6)') nstsi,2,ikti+1,1,irsep-1
        WRITE(fp2,'(6I6)') -3,0,0,ikti+2,0,3

        nstsi = nstsi + 1
        WRITE(fp2,'(A,I2,2A)')
     .    '* POLOIDAL FLUXSURFACE ',ikti+2,', TRANSPARENT RIGHT',
     .    ' PFZ/CORE'
        WRITE(fp2,'(5I6)') nstsi,2,ikti+2,1,irsep-1
        WRITE(fp2,'(6I6)') -3,0,0,ikto+1,0,3

        nstsi = nstsi + 1
        WRITE(fp2,'(A,I2,A)')
     .    '* POLOIDAL FLUXSURFACE ',ikto+1,', TRANSPARENT LEFT SOL'
        WRITE(fp2,'(5I6)') nstsi,2,ikto+1,irsep-1,irwall-1
        WRITE(fp2,'(6I6)') -3,0,0,ikto+2,0,4

        nstsi = nstsi + 1
        WRITE(fp2,'(A,I2,A)')
     .    '* POLOIDAL FLUXSURFACE ',ikto+2,', TRANSPARENT LEFT SOL'
        WRITE(fp2,'(5I6)') nstsi,2,ikto+2,irsep-1,irwall-1
        WRITE(fp2,'(6I6)') -3,0,0,ikti+1,0,4

        nstsi = nstsi + 1
        WRITE(fp2,'(A,I2,A)')
     .    '* POLOIDAL FLUXSURFACE ',ikti+1,', TRANSPARENT RIGHT SOL'
        WRITE(fp2,'(5I6)') nstsi,2,ikti+1,irsep-1,irwall-1
        WRITE(fp2,'(6I6)') -3,0,0,ikti+2,0,4

        nstsi = nstsi + 1
        WRITE(fp2,'(A,I2,A)')
     .    '* POLOIDAL FLUXSURFACE ',ikti+1,', TRANSPARENT RIGHT SOL'
        WRITE(fp2,'(5I6)') nstsi,2,ikti+2,irsep-1,irwall-1
        WRITE(fp2,'(6I6)') -3,0,0,ikto+1,0,4

30      CALL ReadLine(fp1,buffer,1,*97,*98)
        IF (buffer(1:3).NE.'***') GOTO 30
        BACKSPACE fp1

      ELSEIF ((buffer(1:6).EQ.'*** 3b'.OR.
     .         buffer(1:6).EQ.'*** 3B').AND.eirneut.EQ.1) THEN
c
c       This section details the neutral wall segments (referred to
c       as 'additional surfaces' in EIRENE).  The geometry
c       data for these surfaces comes from either the neutral wall
c       data in the DIVIMP input file or the grid file (JET grids).
c
c       The ADDUSR routine in EIRENE clips surfaces to fit the targets,
c       but this has already been done in DIVIMP so the ADDUSR
c       routine should not be executed (see EIREADD DIVIMP option):
c
        WRITE(fp2,'(A)')
     .    '*** 3b. DATA FOR ADDITIONAL SURFACES (DIVIMP)'
        WRITE(fp2,'(I6)') walln

        DO i1 = 1, walln
          WRITE(fp2,'(A,I4,A)')
     .      '*',i1,' :'
          WRITE(fp2,'(A)')
     .      ' 2.00000E+00 1.00000E+00 1.00000E-05'
          WRITE(fp2,'(A)')
     .      '     1     2     0     0     0     1     0     0     0'
          WRITE(fp2,'(6E14.7)')
     .      wallr(i1,1) * 100.0, wallz(i1,1) * 100.0,-1.0E+20,
     .      wallr(i1,2) * 100.0, wallz(i1,2) * 100.0, 1.0E+20
          WRITE(fp2,'(A)')
     .      '     1     0     0     0'

c Replace with an array...
          WRITE(fp2,'(1P,2E12.4)') material(eirmat2),eirtemp2

          WRITE(fp2,'(2A)')
     .      ' 1.00000E+00 1.00000E+00 0.00000E+00 1.00000E+00',
     .      ' 5.00000E-01 1.00000E+00'
        ENDDO
c
40      CALL ReadLine(fp1,buffer,1,*97,*98)
        IF (buffer(1:3).NE.'***') GOTO 40
        BACKSPACE fp1

      ELSEIF (buffer(1:6).EQ.'*** 7.') THEN
        WRITE(fp2,'(A)')
     .    '*** 7. DATA FOR PRIMARY SOURCES OF NEUTRALS (DIVIMP)'

        CALL TransferLine(fp1,fp2,buffer,3)

        CALL ReadLine(fp1,buffer,1,*97,*98)
        WRITE(fp2,'(A,I3)') '* POLOIDAL POLYGON ',1

        CALL TransferLine(fp1,fp2,buffer,7)
        CALL UpdateLine1I(fp1,fp2,buffer,5,irwall-1)
        CALL TransferLine(fp1,fp2,buffer,5)

        CALL ReadLine(fp1,buffer,1,*97,*98)
        WRITE(fp2,'(A,I3)') '* POLOIDAL POLYGON ',nks(irsep)+3

        CALL TransferLine(fp1,fp2,buffer,7)
        CALL UpdateLine2I(fp1,fp2,buffer,3,5,nks(irsep)+3,irwall-1)
        CALL TransferLine(fp1,fp2,buffer,5)

        WRITE(fp2,'(A)') '* VOLUME SOURCE '
        CALL TransferLine(fp1,fp2,buffer,13)


      ELSEIF (buffer(1:6).EQ.'*** 11') THEN
        WRITE(fp2,'(A)')
     .    '*** 11. DATA FOR NUMERICAL AND GRAPHICAL OUTPUT (DIVIMP)'

        CALL TransferLine(fp1,fp2,buffer,14)
        CALL UpdateLine1I(fp1,fp2,buffer,2,nks(irsep)+2)
        CALL TransferLine(fp1,fp2,buffer,13)

      ELSEIF (buffer(1:6).EQ.'*** 14') THEN
        WRITE(fp2,'(A)')
     .    '*** 14. DATA FOR INTERFACING (DIVIMP)'

        CALL TransferLine(fp1,fp2,buffer,4)
        CALL UpdateLine2I(fp1,fp2,buffer,1,2,nks(irsep)+2,irwall-2)
        CALL TransferLine(fp1,fp2,buffer,2)
        CALL UpdateLine1I(fp1,fp2,buffer,6,irwall-1)
        CALL UpdateLine2I(fp1,fp2,buffer,2,6,nks(irsep)+2,irwall-1)
        CALL TransferLine(fp1,fp2,buffer,8)

        CALL ReadLine(fp1,buffer,1,*50,*98)
        CALL ER("WriteInputFile","Expected end of file",*99)

      ELSE
c
c     Section identifier is not recognized, so
c     copy the entire section as is:
c
        CALL WriteLine(fp2,buffer)
        GOTO 10
      ENDIF

      CALL ReadLine(fp1,buffer,1,*97,*98)
      IF (buffer(1:3).NE.'***') THEN
        CALL ER('WriteInputFile','Invalid template format',*99)
      ENDIF

      GOTO 20

50    CONTINUE

      CLOSE (fp1)
      CLOSE (fp2)

      IF (eirneut.EQ.0) THEN
        nvesm = 0
        write(0,*)
        write(0,*) ' TEMPORARY BLANKING OF NEUTRAL WALL DATA! '
        write(0,*)
      ENDIF

      RETURN
c
c     Error code:
c
95    CALL ER('WriteInputFile','File error',*99)
97    CALL ER('WriteInputFile','Unexpected end of file',*99)
98    CALL ER('WriteInputFile','Problems reading template file',*99)
99    WRITE(EROUT,*) '  Last line read: '
      WRITE(EROUT,*) '  "',buffer,'"'
      STOP
      END
c
c ======================================================================
c
c subroutine: WriteGeometryFile
c
c
c
c
c
      SUBROUTINE WriteGeometryFile_97
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER nnks,dimxh,dimyh,nncut,nxcut1,nxcut2

      INTEGER           ik,ik1,ir1,ir,id
      DOUBLE PRECISION  x1,x2,x3,x4,y1,y2,y3,y4
c
c     Check to see if DIVIMP should write the geometry
c     file for EIRENE:
c
      IF (eirgeom.EQ.0) RETURN
c
c     Initialisation:
c
      nnks   = nks(irsep)
      dimxh  = nnks
      dimyh  = irwall-2
      nncut  = 2
      nxcut1 = ikto
      nxcut2 = ikti

      OPEN(UNIT=52,ACCESS='SEQUENTIAL',STATUS='REPLACE')
c
c     Write the header to the geometry file:
c
      WRITE(52,'(A)') 'DIVIMP/EIRENE Geometry Data'
      WRITE(52,*)
      WRITE(52,'(A30,I6)') 'dimhx                         ',dimxh
      WRITE(52,'(A30,I6)') 'dimyh                         ',dimyh
      WRITE(52,'(A30,I6)') 'nncut                         ',nncut
      WRITE(52,'(A30,I6)') 'nxcut1+1                      ',nxcut1+1
      WRITE(52,'(A30,I6)') 'nxcut2                        ',nxcut2
      WRITE(52,'(A30,I6)') '[dummy]                       ',0
      WRITE(52,'(A30,I6)') '[dummy]                       ',0
      WRITE(52,'(A30,I6)') '[dummy]                       ',0
      WRITE(52,*)

      DO ik = 1, nnks
        DO ir = 2, irwall
c
c         Determine the cell index:
c
          IF (ir.LT.irsep) THEN
            IF (ik.LE.nxcut1) THEN
              ik1 = ik
              ir1 = irtrap+ir-1
            ELSEIF (ik.GE.nxcut2) THEN
              ik1 = ik-(nxcut2-nxcut1-1)
              ir1 = irtrap+ir-1
            ELSE
              ik1 = ik-nxcut1
              ir1 = ir
            ENDIF
          ELSEIF (ir.EQ.irwall) THEN
            ik1 = ik
            ir1 = ir-1
          ELSE
            ik1 = ik
            ir1 = ir
          ENDIF

          id = korpg(ik1,ir1)
c
c         Write geometry data:
c
          x1 = rvertp(1,id)
          y1 = zvertp(1,id)
          x2 = rvertp(2,id)
          y2 = zvertp(2,id)

          IF (ir.EQ.irwall) THEN
            WRITE(52,3334)
     .        x2,y2,0.0,0.0,ir1,ik1,' n'
          ELSE
            WRITE(52,3334)
     +        x1,y1,x2,y2,ir1,ik1,' n'
          ENDIF
        ENDDO
c
c       Add in the extra cell information around cut points:
c
        IF (ik.EQ.nxcut1.OR.ik.EQ.nxcut2-1.OR.ik.EQ.nnks) THEN
          DO ir = 2, irwall
c
c           Determine the cell index:
c
            IF (ir.LT.irsep) THEN
              IF (ik.EQ.nxcut1) THEN
                ik1 = ik
                ir1 = irtrap+ir-1
              ELSEIF (ik.EQ.nxcut2-1) THEN
                ik1 = nks(ir)-1
                ir1 = ir
              ELSEIF (ik.EQ.nnks) THEN
                ik1 = nks(irtrap+ir-1)
                ir1 = irtrap+ir-1
              ENDIF
            ELSEIF (ir.EQ.irwall) THEN
              ik1 = ik
              ir1 = ir - 1
            ELSE
              ik1 = ik
              ir1 = ir
            ENDIF

            id = korpg(ik1,ir1)

            x3 = rvertp(3,id)
            y3 = zvertp(3,id)
            x4 = rvertp(4,id)
            y4 = zvertp(4,id)

            IF ((ir.EQ.irwall.AND.ik.EQ.nxcut1  ).OR.
     .          (ir.EQ.irwall.AND.ik.EQ.nxcut2-1).OR.
     .          (ir.EQ.irwall.AND.ik.EQ.nnks    )) THEN
              WRITE(52,3334)
     .          x3,y3,0.0,0.0,ir1,ik1,' p'
            ELSE
              WRITE(52,3334)
     .          x4,y4,x3,y3,ir1,ik1,' p'
            ENDIF
          ENDDO
        ENDIF

      ENDDO
      CLOSE(52)

      RETURN
3334  FORMAT(4E15.7,5X,2I6,2A)
      END


c
c ======================================================================
c
c subroutine: BuildTriangles
c
      SUBROUTINE BuildTriangles
      USE mod_eirene04
      IMPLICIT none

c DEFUNCT

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'


      REAL       TOL
      PARAMETER (TOL=1.0E-07)

      REAL GetMach,GetJsat,GetFlux 

      INTEGER fp,fp2,ik1,ik2,ir1,ir2,i1,i2,trin,id,nl,v1,v1n,v2,v2n,
     .        vern,in,tarside,region
      LOGICAL status
      REAL    x1,x2,y1,y2,deltax,deltay,Bfrac,fact,
     .        tarte,tarti,tarne,tarv,tarflux,tarisat,tarM
      REAL*8  Bx,By,Bz,beta,brat

      INTEGER, ALLOCATABLE :: nlink1(:,:)    ,nlink2(:,:),
     .                        llink1(:,:,:,:),llink2(:,:,:,:),
     .                        triik(:),triir(:),trisf(:),
     .                        trimap(:,:),triside(:,:),trivert(:,:),
     .                        trisurface(:,:)
      REAL, ALLOCATABLE :: trix(:,:),verx(:),verz(:),
     .                     triy(:,:),very(:)


c      RETURN

      WRITE(0,*) 'BUILDING TRIANGLES'

      CALL BuildMap

      CALL OutputData(85,'Building triangles')


      ALLOCATE(nlink1(MAXNKS,MAXNRS))
      ALLOCATE(nlink2(MAXNKS,MAXNRS))
      ALLOCATE(llink1(MAXNKS,MAXNRS,100,2))
      ALLOCATE(llink2(MAXNKS,MAXNRS,100,2))

      ALLOCATE(triik(MAXNKS*MAXNRS))
      ALLOCATE(triir(MAXNKS*MAXNRS))
      ALLOCATE(trisf(MAXNKS*MAXNRS))
      ALLOCATE(trix(MAXNKS*MAXNRS,3))
      ALLOCATE(triy(MAXNKS*MAXNRS,3))
      ALLOCATE(trimap(MAXNKS*MAXNRS,3))
      ALLOCATE(triside(MAXNKS*MAXNRS,3))
      ALLOCATE(trisurface(MAXNKS*MAXNRS,3))
      ALLOCATE(trivert(MAXNKS*MAXNRS,3))

      ALLOCATE(verx(MAXNKS*MAXNRS))
      ALLOCATE(very(MAXNKS*MAXNRS))
      ALLOCATE(verz(MAXNKS*MAXNRS))

      vern = 0
      verx = 0.0
      very = 0.0
      verz = 0.0


      fp = 99
      fp2 = 98


c...  Build a list of reference quadrilaterals, and a list of connection
c     map references to each cell:    

      nlink1 = 0
      nlink2 = 0
      DO ir1 = 2, nrs
        IF (idring(ir1).EQ.-1) CYCLE
        DO ik1 = 1, nks(ir1)
          DO ir2 = 1, nrs
            IF (idring(ir2).EQ.BOUNDARY) CYCLE
            DO ik2 = 1, nks(ir2)
              IF (ir2.LT.irsep.AND.ik2.EQ.nks(ir2)) CYCLE
              IF     (irouts(ik2,ir2).EQ.ir1.AND.
     .                ikouts(ik2,ir2).EQ.ik1) THEN

                nlink1(ik1,ir1) = nlink1(ik1,ir1) + 1
                llink1(ik1,ir1,nlink1(ik1,ir1),1) = ik2
                llink1(ik1,ir1,nlink1(ik1,ir1),2) = ir2
              ELSEIF (irins (ik2,ir2).EQ.ir1.AND.
     .                ikins (ik2,ir2).EQ.ik1) THEN
                nlink2(ik1,ir1) = nlink2(ik1,ir1) + 1
                llink2(ik1,ir1,nlink2(ik1,ir1),1) = ik2
                llink2(ik1,ir1,nlink2(ik1,ir1),2) = ir2
              ENDIF
            ENDDO
          ENDDO
c...      
          IF (nlink1(ik1,ir1).LE.1) nlink1(ik1,ir1) = 0
          IF (nlink2(ik1,ir1).LE.1) nlink2(ik1,ir1) = 0
c...      Identify cells that adjacent to boundary rings:
          IF (idring(irins (ik1,ir1)).EQ.-1) nlink1(ik1,ir1) = -1
          IF (idring(irouts(ik1,ir1)).EQ.-1) nlink2(ik1,ir1) = -1
        ENDDO
c        DO ik1 = 1, nks(ir1)     
c         IF (ir1.EQ.irsep.AND.ik1.EQ.25) THEN
c           WRITE(0,*)  'DUMP:',ik1,nlink1(ik1,ir1),nlink2(ik1,ir1)
c         ENDIF
c        ENDDO
      ENDDO 

c...  Construct 2 triangles from each regular OEDGE cell:
c
c     counter-clockwise
c     Triangle 1 contains OEDGE cell verticies 1 and 2, and 
c     triangle 2 contains verticies 3 and 4.
c
c
      DO ir1 = 2, nrs
        IF (idring(ir1).EQ.-1) CYCLE

c        IF (ir1.NE.irtrap+1.AND.ir1.NE.irtrap+2) CYCLE
c        IF (ir1.NE.irsep.AND.ir1.NE.nrs.AND.ir1.NE.irsep+1.AND.
c     .      ir1.NE.irsep-1) CYCLE
c        IF (ir1.LT.irtrap) CYCLE

        DO ik1 = 1, nks(ir1)
          IF (ir1.LT.irsep.AND.ik1.EQ.nks(ir1)) CYCLE
 
          id = korpg(ik1,ir1)

          trin = trin + 1
          triik(trin) = ik1
          triir(trin) = ir1
          trisf(trin) = 23
          trix(trin,1) = rvertp(1,id)         
          triy(trin,1) = zvertp(1,id) 
          trix(trin,2) = rvertp(3,id) 
          triy(trin,2) = zvertp(3,id) 
          trix(trin,3) = rvertp(2,id) 
          triy(trin,3) = zvertp(2,id) 

          trin = trin + 1
          triik(trin) = ik1
          triir(trin) = ir1
          trisf(trin) = 14
          trix(trin,1) = rvertp(3,id) 
          triy(trin,1) = zvertp(3,id) 
          trix(trin,2) = rvertp(1,id) 
          triy(trin,2) = zvertp(1,id) 
          trix(trin,3) = rvertp(4,id) 
          triy(trin,3) = zvertp(4,id) 
        ENDDO
      ENDDO

c.... Slice triangles with more than one neighbour:
      i1 = 0
      status = .TRUE.
      DO WHILE (status)
        i1 = i1 + 1

        ik1 = triik(i1)
        ir1 = triir(i1)

        IF     (trisf(i1).EQ.23) THEN
          IF (nlink2(ik1,ir1).LE.0) CYCLE

          nl = nlink2(ik1,ir1)

c...      Make room for new triangles:
          DO i2 = trin, i1+1, -1
            triik(i2+nl-1)   = triik(i2)
            triir(i2+nl-1)   = triir(i2)
            trisf(i2+nl-1)   = trisf(i2)
            trix (i2+nl-1,1) = trix (i2,1)
            triy (i2+nl-1,1) = triy (i2,1)
            trix (i2+nl-1,2) = trix (i2,2)
            triy (i2+nl-1,2) = triy (i2,2)
            trix (i2+nl-1,3) = trix (i2,3)
            triy (i2+nl-1,3) = triy (i2,3)
          ENDDO
          trin = trin + nl - 1
  
c...      Duplicate triangle:
          DO i2 = i1+1, i1+nl-1
            triik(i2) =  triik(i1)
            triir(i2) =  triir(i1)
            trisf(i2) = -trisf(i1)
            trix(i2,1) = trix(i1,1)
            triy(i2,1) = triy(i1,1)
            trix(i2,2) = trix(i1,2)
            triy(i2,2) = triy(i1,2)
            trix(i2,3) = trix(i1,3)
            triy(i2,3) = triy(i1,3)
          ENDDO

c...      Update trailing vertex from connection map data:
          DO i2 = 2, nl
            ik2 = llink2(ik1,ir1,i2,1)
            ir2 = llink2(ik1,ir1,i2,2)
            id = korpg(ik2,ir2)
            trix(i1+i2-1,3) = rvertp(1,id) 
            triy(i1+i2-1,3) = zvertp(1,id) 
          ENDDO

c...      Make triangles in this series consistent with each subsequent
c         triangle in the series:
          DO i2 = 1, nl-1
            trix(i1+i2-1,2) = trix(i1+i2,3)  
            triy(i1+i2-1,2) = triy(i1+i2,3) 
          ENDDO

        ELSEIF (trisf(i1).EQ.14) THEN
          IF (nlink1(ik1,ir1).LE.0) CYCLE

          nl = nlink1(ik1,ir1)

c...      Make room for new triangles:
          DO i2 = trin, i1+1, -1
            triik(i2+nl-1)   = triik(i2)
            triir(i2+nl-1)   = triir(i2)
            trisf(i2+nl-1)   = trisf(i2)
            trix (i2+nl-1,1) = trix (i2,1)
            triy (i2+nl-1,1) = triy (i2,1)
            trix (i2+nl-1,2) = trix (i2,2)
            triy (i2+nl-1,2) = triy (i2,2)
            trix (i2+nl-1,3) = trix (i2,3)
            triy (i2+nl-1,3) = triy (i2,3)
          ENDDO
          trin = trin + nl - 1
  
c...      Duplicate triangle:
          DO i2 = i1+1, i1+nl-1
            triik(i2) =  triik(i1)
            triir(i2) =  triir(i1)
            trisf(i2) = -trisf(i1)
            trix(i2,1) = trix(i1,1)
            triy(i2,1) = triy(i1,1)
            trix(i2,2) = trix(i1,2)
            triy(i2,2) = triy(i1,2)
            trix(i2,3) = trix(i1,3)
            triy(i2,3) = triy(i1,3)
          ENDDO

c...      Update trailing vertex from connection map data:
          DO i2 = 2, nl
            ik2 = llink1(ik1,ir1,i2,1)
            ir2 = llink1(ik1,ir1,i2,2)
            id = korpg(ik2,ir2)
            trix(i1+i2-1,2) = rvertp(2,id) 
            triy(i1+i2-1,2) = zvertp(2,id) 
            IF (ik1.EQ.25.AND.
     .        ir1.EQ.irsep) 
     .        WRITE(0,*) 'ik,ir2=',ik2,ir2,rvertp(2,id),zvertp(2,id),nl
          ENDDO

c...      Make triangles in this series consistent with each subsequent
c         triangle in the series:
          DO i2 = 1, nl-1
            trix(i1+i2-1,3) = trix(i1+i2,2)  
            triy(i1+i2-1,3) = triy(i1+i2,2) 
          ENDDO

        ENDIF

        status = .FALSE.
        IF (i1.LT.trin) status = .TRUE.
      ENDDO


c         DO i2 = 1, trin
c           WRITE(0,'(A,I4,3F10.5)') '->',i2,trix(i2,1),triy(i2,1)
c           WRITE(0,'(A,I4,3F10.5)') '  ',i2,trix(i2,2),triy(i2,2)
c           WRITE(0,'(A,I4,3F10.5)') '  ',i2,trix(i2,3),triy(i2,3)
c         ENDDO
 


c...  Build connection map:

c...  THIS CAN BE MORE EFFICIENT:
      DO i1 = 1, trin
        DO v1 = 1, 3
          v1n = v1 + 1
          IF (v1n.EQ.4) v1n = 1
          DO i2 = 1, trin
            IF (i1.EQ.i2) CYCLE
            DO v2 = 1, 3    
              v2n = v2 - 1
              IF (v2n.EQ.0) v2n = 3
             
              IF (ABS(trix(i1,v1 )-trix(i2,v2 )).LT.TOL.AND.
     .            ABS(triy(i1,v1 )-triy(i2,v2 )).LT.TOL.AND.
     .            ABS(trix(i1,v1n)-trix(i2,v2n)).LT.TOL.AND.
     .            ABS(triy(i1,v1n)-triy(i2,v2n)).LT.TOL) THEN
                trimap (i1,v1) = i2
                triside(i1,v1) = v2n
              ENDIF

            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO i2 = 1, trin
        DO v1 = 1, 3
          IF (trimap(i2,v1).EQ.0) THEN
c...        Identify neighbourless surfaces and assign the appropriate mapping
c           to the appropriate non-default standard surface:
            ik1 = triik(i2)
            ir1 = triir(i2)
            IF     (v1.EQ.3.AND.ik1.EQ.1       .AND.ir1.GE.irsep) THEN
c...          Low IK index target segment, EIRENE poloidal suface 4:
              trisurface(i2,v1) = 2                            
c              trisurface(i2,v1) = 4                            
            ELSEIF (v1.EQ.3.AND.ik1.EQ.nks(ir1).AND.ir1.GE.irsep) THEN
c...          High IK index target segment, EIRENE poloidal surface 5:
              trisurface(i2,v1) = 3
c              trisurface(i2,v1) = 5                            
            ELSEIF (v1.EQ.2.AND.irins(ik1,ir1) .EQ.1) THEN
c...          Core boundary, EIRENE radial surface 1:
              trisurface(i2,v1) = 1                            
            ELSEIF ((v1.EQ.2.AND.irouts(ik1,ir1).EQ.irwall).OR.
c...          Added for MAST-DN secondary PFZ:
     .              (v1.EQ.2.AND.irins (ik1,ir1).EQ.irwall)) THEN
c...          SOL boundary, EIRENE radial surface 4:
              trisurface(i2,v1) = 4                            
c              trisurface(i2,v1) = 2                            
            ELSEIF (v1.EQ.2.AND.irins(ik1,ir1) .EQ.irtrap) THEN
c...          PFZ boundary, EIRENE radial surface 5:
              trisurface(i2,v1) = 5
c              trisurface(i2,v1) = 3                            
            ELSEIF (v1.EQ.1.OR.v1.EQ.3) THEN
              WRITE(0,*) 'PROBLEM 3:',i2,v1,triik(i2),triir(i2)
            ENDIF
             
          ELSEIF (.FALSE..AND.v1.EQ.2.AND.ABS(trisf(i2)).EQ.23.AND.  
     .            (triir(i2).EQ.irsep-1.OR.triir(i2).EQ.nrs)) THEN
c...        Label "right" side of separatrix as EIRENE non-standard surface 6:
            trisurface(i2,v1) = 6
          ELSEIF (.FALSE..AND.v1.EQ.2.AND.ABS(trisf(i2)).EQ.14.AND.
     .            triir(i2).EQ.irsep) THEN
c...        Label "left" side of separatrix as EIRENE non-standard surface 7:
            trisurface(i2,v1) = 7

          ELSEIF (trimap(trimap(i2,v1),triside(i2,v1)).NE.i2) THEN
            WRITE(0,*) 'PROBLEM 4:',i2,v1,triik(i2),triir(i2)

          ELSEIF (triside(trimap(i2,v1),triside(i2,v1)).NE.v1) THEN
            WRITE(0,*) 'PROBLEM 5:',i2,v1,triik(i2),triir(i2)

          ENDIF
        ENDDO
      ENDDO

c...  All is well, so generate a master list of points with pointers to the respective
c     triangle verticies:
      DO i1 = 1, trin
        DO v1 = 1, 3

c...      Scan list of verticies and assign the relevant vertex
c         number to the triangle:
          status = .FALSE.
          DO i2 = 1, vern
            IF (ABS(verx(i2)-trix(i1,v1)).LT.TOL.AND.
     .          ABS(very(i2)-triy(i1,v1)).LT.TOL) THEN
              IF (trivert(i1,v1).EQ.0) THEN
                trivert(i1,v1) = i2
                status = .TRUE.
              ELSE
                CALL ER('BuildTriangles','Redundant vertex',*99)
              ENDIF
            ENDIF
          ENDDO
c...      Vertex was not found, so add it to the list:
          IF (.NOT.status) THEN
            vern = vern + 1
            verx(vern) = trix(i1,v1)
            very(vern) = triy(i1,v1)
            trivert(i1,v1) = vern
          ENDIF

        ENDDO
      ENDDO

c      WRITE(0,*) 'vern:',vern,trin






c...  Dump triangles:
c      OPEN(UNIT=fp,FILE='triangles.points',ACCESS='SEQUENTIAL',
c     .     STATUS='REPLACE',ERR=96)      
      nver = vern
c      WRITE(fp,*) vern
      DO i1 = 1, vern
        ver(i1,1) = verx(i1)
        ver(i1,2) = very(i1)
c        WRITE(fp,'(I6,3F12.6))') i1,verx(i1)*100.0,very(i1)*100.0,0.0
      ENDDO
c      CLOSE(fp)      

c...  Dump triangles:
c      OPEN(UNIT=fp,FILE='triangles.sides',ACCESS='SEQUENTIAL',
c     .     STATUS='REPLACE',ERR=96)      
      ntri = trin
c      WRITE(fp,*) trin
      DO i1 = 1, trin
        DO v1 = 1, 3
          tri(i1)%ver(v1) = trivert(i1,v1)
        ENDDO  
c        WRITE(fp,'(I6,4X,3I6))') i1,(trivert(i1,v1),v1=1,3)
      ENDDO
c      CLOSE(fp)      

c...  Dump triangles:
c      OPEN(UNIT=fp,FILE='triangles.map',ACCESS='SEQUENTIAL',
c     .     STATUS='REPLACE',ERR=96)      
c      WRITE(fp,*) trin
      DO i1 = 1, trin
        tri(i1)%sideindex(1,2) = ABS(trisf(i1))
        DO v1 = 1, 3
          tri(i1)%map(v1) = 0 ! trimap(i1,v1)
          tri(i1)%sid(v1) = triside(i1,v1)
          tri(i1)%sur(v1) = 0 !trisurface(i1,v1)
          IF (trisurface(i1,v1).EQ.2) tri(i1)%sideindex(2,v1) = IKLO
          IF (trisurface(i1,v1).EQ.3) tri(i1)%sideindex(2,v1) = IKHI
        ENDDO  
c        WRITE(fp,'(I6,4X,3(3I6,4X),2I6)') i1,
c     .    (trimap(i1,v1),triside(i1,v1),trisurface(i1,v1),v1=1,3),
c     .    triik(i1),triir(i1)
      ENDDO
c      CLOSE(fp)      


      fact = qtim * qtim * emi / crmi

c...  Dump plasma data:
c      OPEN(UNIT=fp ,FILE='triangles.plasma',ACCESS='SEQUENTIAL',
c     .     STATUS='REPLACE',ERR=96)      
c      OPEN(UNIT=fp2,FILE='triangles.efield',ACCESS='SEQUENTIAL',
c     .     STATUS='REPLACE',ERR=96)      

c...  Header:
c      WRITE(fp,'(A)') '* VERSION 1.0 OSM PLASMA FILE FOR TRIANGULAR '//
c     .                'GRID'
c      WRITE(fp,'(A)') '*'
c      WRITE(fp,'(A)') '* BULK PLASMA DATA'
c      WRITE(fp,'(A)') '*'
c      WRITE(fp,'(A7,2A8,2X,A10,2X,3A10,2X,3A12)') 
c     .  '* Index','Te','Ti','ne','vx',
c     .  'vy','vz','Bx','By','Bz'
c      WRITE(fp,'(A7,2A8,2X,A10,2X,3A10,2X,3A12)') 
c     .  '*      ','(eV)','(eV)','(cm-3)','(cm s-1)',
c     .  '(cm s-1)','(cm s-1)','(Tesla)','(Tesla)','(Tesla)'


c...  Header:
c      WRITE(fp2,'(A)') '* VERSION 1.0 OSM PLASMA FILE FOR TRIANGULAR '//
c     .                'GRID'
c      WRITE(fp2,'(A)') '*'
c      WRITE(fp2,'(A)') '* BULK PLASMA DATA (PART DEUX)'
c      WRITE(fp2,'(A)') '*'
c      WRITE(fp2,'(A7,3A12)') 
c     .  '* Index','Ex','Ey','Ez'
c      WRITE(fp2,'(A7,3A12)')
c     .  '*      ','V m-1','V m-1','V m-1'

      
c      WRITE(fp,*) trin
      DO i1 = 1, trin
        ik1 = triik(i1)
        ir1 = triir(i1)
        id = korpg(ik1,ir1)
        x1 = 0.5 * (rvertp(1,id) + rvertp(2,id))
        y1 = 0.5 * (zvertp(1,id) + zvertp(2,id))
        x2 = 0.5 * (rvertp(3,id) + rvertp(4,id))
        y2 = 0.5 * (zvertp(3,id) + zvertp(4,id))
        deltax = (x2 - x1)
        deltay = (y2 - y1)
        brat = DBLE(bratio(ik1,ir1))
        beta = DBLE(deltax / deltay)
        Bz = DSQRT(1.0 - brat**2)
        By = brat * DSQRT(1.0D0 / (1.0D0 + beta**2)) * 
     .       DBLE(SIGN(1.0,deltay))
        Bx = beta * By
c...    CBPHI is the on-axis B-field value specified in the OSM input file:
        Bfrac = cbphi * r0 / rs(ik1,ir1) 

        tri(i1)%type = MAGNETIC_GRID

        tri(i1)%index(1) = ik1
        tri(i1)%index(2) = ir1

        tri(i1)%plasma(1) = ktebs(ik1,ir1)  ! 5.0
        tri(i1)%plasma(2) = ktibs(ik1,ir1)  ! 5.0
        tri(i1)%plasma(3) = knbs (ik1,ir1)  ! 1.0E+12
        tri(i1)%plasma(4) = SNGL(Bx)*kvhs(ik1,ir1)
        tri(i1)%plasma(5) = SNGL(By)*kvhs(ik1,ir1)
        tri(i1)%plasma(6) = SNGL(Bz)*kvhs(ik1,ir1)

        tri(i1)%bfield(1) = SNGL(Bx)*Bfrac
        tri(i1)%bfield(2) = SNGL(By)*Bfrac
        tri(i1)%bfield(3) = SNGL(Bz)*Bfrac

        tri(i1)%efield(1) = SNGL(Bx)*kes(ik1,ir1)/fact
        tri(i1)%efield(2) = SNGL(By)*kes(ik1,ir1)/fact
        tri(i1)%efield(3) = SNGL(Bz)*kes(ik1,ir1)/fact

c        WRITE(fp,'(I7,2F8.2,2X,1P,E10.2,2X,3E10.2,2X,
c     .             3E12.4,0P,6X,2I4)') i1,
c     .         ktebs(ik1,ir1),ktibs(ik1,ir1),knbs(ik1,ir1)*1.0E-06,
c     .         SNGL(Bx)*kvhs(ik1,ir1)*100.0,  !/qtim,
c     .         SNGL(By)*kvhs(ik1,ir1)*100.0,  !/qtim,
c     .         SNGL(Bz)*kvhs(ik1,ir1)*100.0,  !/qtim,
c     .         SNGL(Bx)*Bfrac,SNGL(By)*Bfrac,SNGL(Bz)*Bfrac,
c     .         ik1,ir1

c...    Dump efield data:
c        WRITE(fp2,'(I7,1P,3E12.4,10X,0P,2I4,1P,E12.4)') i1,
c     .         SNGL(Bx)*kes(ik1,ir1)/fact,
c     .         SNGL(By)*kes(ik1,ir1)/fact,
c     .         SNGL(Bz)*kes(ik1,ir1)/fact,
c     .         ik1,ir1,kes(ik1,ir1)/fact

      ENDDO



c...  All done:
c      CLOSE(fp)      
c      CLOSE(fp2)      

c...  Store triangle list in case opacity data is read from a separate
c     EIRENE run:
c      IF (trin.GT.MAXNKS*MAXNRS) 
c     .  CALL ER('DumpTrinagles','Arrays bounds error',*99)
c      eirtrin = trin
c      DO i1 = 1, trin
c        eirtriik(i1) = triik(i1)
c        eirtriir(i1) = triir(i1)
c      ENDDO

c...  mod_eirene04:
c      ntri = trin
c      nver = vern


      DEALLOCATE(nlink1)
      DEALLOCATE(nlink2)
      DEALLOCATE(llink1)
      DEALLOCATE(llink2)

      DEALLOCATE(triik)
      DEALLOCATE(triir)
      DEALLOCATE(trix)
      DEALLOCATE(triy)
      DEALLOCATE(trimap)
      DEALLOCATE(triside)
      DEALLOCATE(trisurface)
      DEALLOCATE(trivert)

      DEALLOCATE(verx)
      DEALLOCATE(very)
      DEALLOCATE(verz)

c      STOP 'DUMPING'
      WRITE(0,*) 'DONE'

      RETURN
96    WRITE(0,*) 'BUILDTRIANGES: PROBLEMS WITH FILE ACCESS'
      STOP
99    STOP
      END
c
