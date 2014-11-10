c slmod begin
c
c ======================================================================
c
      SUBROUTINE USRTIM

      USE PARMMOD
      USE COMPRT
      USE CCONA

      IMPLICIT NONE

      INTEGER i1,i2,fp,idum1
      REAL*8 totbinc
      CHARACTER*1024 cdum1

c...  Set up parameters:
      ENTRY USRTI0

      fp = 31

      READ(fp,'(A)',ERR=97,END=98) cdum1

      IF (cdum1(1:10).NE.'[TIME-TO-I') THEN
        WRITE(0,*) 'ERROR: TIME-TO-IONISATION PARAMETERS EXPECTED '//
     .             'BUT NOT FOUND'
        WRITE(0,*) 'CDUM1:',cdum1
        STOP
      ENDIF

      READ(fp,'(A)',ERR=97,END=98) cdum1
      DO WHILE(cdum1(1:1).EQ.'*')
        READ(fp,'(A)',ERR=97,END=98) cdum1     
      ENDDO  
      BACKSPACE fp

c...  Load the number of regions for recording time-to-ionisation
c     information:
      READ(fp,*) timnum

      DO i1 = 1, timnum

c...    Read bounding rectangle(s):
        READ(fp,*) idum1,timx1(i1),timx2(i1),timy1(i1),timy2(i1),
     .             timnbin(i1)
        IF (timnbin(i1).NE.0) THEN
c...      Load bin information:
          timnbin(i1) = timnbin(i1) + 2
          DO i2 = 1, timnbin(i1)+1
            READ(fp,*) timbin(i1,i2)             
          ENDDO      
        ENDIF
c...    Scale units from [m] to [cm]:
        timx1(i1) = timx1(i1) * 100.0D0
        timx2(i1) = timx2(i1) * 100.0D0
        timy1(i1) = timy1(i1) * 100.0D0
        timy2(i1) = timy2(i1) * 100.0D0
      ENDDO

      RETURN

c...  Process atom ionisation events: 
      ENTRY USRTI1

c...  Check time recording regions to see if particle
c     qualifies:
      DO i1 = 1, timnum
        IF (x0.GT.timx1(i1).AND.x0.LT.timx2(i1).AND.
     .      y0.GT.timy1(i1).AND.y0.LT.timy2(i1)) THEN

          WRITE(6,*) '--TIME-->',npanu,i1


c...      Update the particle count and log time:
          timcnt(i1) = timcnt(i1) + weight
          timavg(i1) = timavg(i1) + weight * time

c...      Bin the particle time:
          DO i2 = 1, timnbin(i1)
            IF (time.GE.timbin(i1,i2).AND.time.LT.timbin(i1,i2+1)) 
     .        timbinc(i1,i2) = timbinc(i1,i2) + weight               
          ENDDO

        ENDIF

      ENDDO

      RETURN


c...  Dump results to the transfer file:
      ENTRY USRTI2

      fp = 32

      WRITE(fp,'(A)') '[TIME-TO-IONISATION STATISTICS]'

      WRITE(fp,'(A8,4A8,2X,2A14)') 'REGION','X1','X2','Y1','Y2','COUNT',
     .                            'AVERAGE'
      DO i1 = 1, timnum
        timavg(i1) = timavg(i1) / (timcnt(i1) + EPS10)

        WRITE(fp,'(I8,4F8.3,2X,F14.4,1P,E14.4,0P)') 
     .    i1,timx1(i1)*0.01D0,timx2(i1)*0.01D0,
     .       timy1(i1)*0.01D0,timy2(i1)*0.01D0,
     .       timcnt(i1),timavg(i1)

      ENDDO

      DO i1 = 1, timnum
        IF (timnbin(i1).EQ.0) CYCLE

        totbinc = 0.0D0
        DO i2 = 1, timnbin(i1)
          totbinc = totbinc + timbinc(i1,i2)
        ENDDO

        WRITE(fp,'(2A8,4A14)') 'REGION','BIN','T1','T2',
     .                         'FRACTION','COUNT'
        DO i2 = 1, timnbin(i1)
          WRITE(fp,'(2I8,1P,2E12.4,0P,2F14.4)')
     .      i1,i2,timbin(i1,i2),timbin(i1,i2+1),
     .      timbinc(i1,i2)/(totbinc+EPS10)*100.0D0,timbinc(i1,i2)
        ENDDO

      ENDDO

      RETURN
97    WRITE(0,*) 'ERROR: ERROR READING TRANSFER FILE'
      GOTO 99
98    WRITE(0,*) 'ERROR: UNEXPECTED EOF ON TRANSFER FILE'
99    STOP
      END
c
c ======================================================================
c
      SUBROUTINE BGKUSR

c NEED TO DUPE A NEW ONE -- ACTUALLY, IT MAY BE BETTER TO PROCESS THE OLD
c ONE BEFORE IT IS OUTPUT TO DIVIMP, SO THAT IT IS READY TO GO ON THE
C NEXT EIRENE CALL

       RETURN
       END
c
c ======================================================================
c
      SUBROUTINE JUMUSR

      USE PARMMOD
      USE CADGEO
      USE CLGIN
      USE CGEOM
      USE COMUSR
      USE CGRID

      IMPLICIT NONE

      REAL*8 V1(3),V2(3),V3(3),V4(3),XMIN,XMAX,YMIN,YMAX

      INTEGER ICOUNT,NMATCH,ASCCODE1,NLIMI1,NSBOX1,NSBOX2,FP,
     .        JSCRTCH(NLIM),I,NC,J,I1,I2,I3,I4,I5,NCELL

c     ASSIGN IGJUM3

      REAL*8 TOL
      PARAMETER (TOL=1.0D-03)

c...  F90 short-hand for this?
      WRITE(0,*) 'ASSIGNING IGJUM4 DEFAULT',nrad
c      STOP 'CHECKING NRAD'
      DO NCELL=1,NRAD
        IGJUM4(NCELL,0)=1
        DO J=1,8
          IGJUM4(NCELL,J)=0
        ENDDO
      ENDDO

c...  Don't proceed if there isn't a vacuum grid:
      IF (SBGKI.GT.EBGKI) THEN
        WRITE(0,*) 'NOT SPECIFYING IGJUM4'
        RETURN
      ENDIF

c...  Attempt to read stored IGJUM3:
c      FP=98
c      OPEN(UNIT=FP,FILE='jummap.dat',FORM='UNFORMATTED',
c     .     STATUS='OLD',ERR=10)
c      READ(FP) ASCCODE1,NLIMI1,NSBOX1,NSBOX2
c      IF (ASCCODE1.EQ.ASCCODE.AND.NLIMI1.EQ.(NLIMI/NOPTM1+1).AND.
c     .    NSBOX1  .EQ.MIN(NOPTIM,NSBOX).AND.NSBOX2.EQ.NSBOX) THEN
c        WRITE(0,*) 'LOADING IGJUM3 FROM DATA FILE'
c        READ(FP) ((IGJUM3(J,I),I=1,NLIMI/NOPTM1+1),
c     .                         J=1,MIN(NOPTIM,NSBOX))
c        READ(FP) ((IGJUM4(J,I),I=0,8),J=1,NSBOX)
c        CLOSE(FP)
c        WRITE(0,*) 'DONE'
c        RETURN
c      ENDIF
c10    CONTINUE
c      CLOSE(FP)
c      WRITE(0,*) 'GENERATING IGJUM3'


      DO J=1,NLIM
        JSCRTCH(J) = 0
      ENDDO
 

c ASSUME ADDITIONAL SURFACES DO NOT OVERLAP THE STANDARD GRID
c SO DO NOT CHECK ADDITIONAL SURFACES WHEN A PARTICLE IS ON THE
c STANDARD GRID.  ONLY TURN OFF POLYGON SURFACES THAT ARE 
c INFINITE IN THE Z DIRECTION, OR THAT ARE PART OF THE VACUUM GRID.

c      WRITE(0,*) 'MARK: JUM PARAMS= ',NSURF,NSURF+ASCNCELL

c      WRITE(0,*) 'WARNING: IGJUM4 SET TO CHECK ALL ADDITIONAL '//
c     .           'SURFACES FROM STANDARD GRID'

      DO NCELL=1,NSURF
        IGJUM4(NCELL,0) =  0
        IF (HSTDI.NE.HADDI) THEN
          IGJUM4(NCELL,0)               =  IGJUM4(NCELL,0) + 1
          IGJUM4(NCELL,IGJUM4(NCELL,0)) = -3
        ENDIF
        IGJUM4(NCELL,0)               =  IGJUM4(NCELL,0) + 1
        IGJUM4(NCELL,IGJUM4(NCELL,0)) = -4

c...    Check the standard additional surfaces as well:
c        IGJUM4(NCELL,0)               =  IGJUM4(NCELL,0) + 1
c        IGJUM4(NCELL,IGJUM4(NCELL,0)) = -1

        IF (NCELL.LE.NOPTIM) THEN
          DO J=1,NLIMI
            IF (P1(3,J).EQ.-1.0D+20.AND.P2(3,J).EQ.1.0D+20) THEN 
              IF (NLIMPB >= NLIMPS) THEN
                IGJUM3(NCELL,J)=1
              ELSE
                CALL BITSET (IGJUM3,0,NOPTIM,NCELL,J,1,NBITS)
              ENDIF
            ENDIF
          ENDDO
c VACUUM GRID CELLS
          DO J=SBGKI,EBGKI
            IF (NLIMPB >= NLIMPS) THEN
              IGJUM3(NCELL,J)=1
            ELSE
              CALL BITSET (IGJUM3,0,NOPTIM,NCELL,J,1,NBITS)
            ENDIF
          ENDDO
        ENDIF

      ENDDO

c...  LGJUM1:
      IF (NLIMPB >= NLIMPS) THEN

      ELSE
        WRITE(0,*) 'IGJUM1 MEMORY SAVING OPTION NOT COMPLETE'
      ENDIF







c LOOP OVER ADDITIONAL SURFACES, AND DECIDE IF THEY ARE 
c COINCIDENT WITH THE POLYGON BOUNDARY OF AN ADDITIONAL CELL

c      DO I1=1,1
      I1 = 0
      NC = 1
      DO I=1,ASCNCELL*ASCNCUT

        I1 = I1 + 1
        IF (I1.EQ.ASCNCELL+1) THEN
          I1 = 1
          NC = NC + 1
        ENDIF

        NCELL=ASCCELL(I1)+(NC-1)*ASCNCELL+NSURF

        IGJUM4(NCELL,0) = 0

        ICOUNT=0

        IF (output)  WRITE(0,*) 'MARK PROGRESS... ',I1,ASCNCELL

c LOOP OVER ALL ADDITIONAL CELLS (VACUUM SURFACES)

        DO J=SBGKI,CBGKI

C
          IF (JSCRTCH(J).EQ.2) CYCLE

c IF AN ADDITIONAL SURFACE IS NOT A SIMPLE LINE SEGMENT WITH
c INFINITE EXTENT, THEN ALWAYS CHECK IT
c          IF (P1(3,J).NE.-1.0D+20.AND.P2(3,J).NE.1.0D+20) CYCLE

c TURN OFF ADDITIONAL SURFACE BY DEFAULT

          IF (NCELL.LE.NOPTIM) THEN
            IF (NLIMPB >= NLIMPS) THEN
              IGJUM3(NCELL,J)=1
            ELSE
              CALL BITSET (IGJUM3,0,NOPTIM,NCELL,J,1,NBITS)
            ENDIF
          ENDIF
c LOOP OVER ALL ADDITIONAL SURFACES AND COMPARE WITH
c ADDITIONAL CELL POLYGON SIDES

          IF (ASC3DMODE.EQ.0.OR.ASC3DMODE.EQ.2) THEN

            DO I3=1,ASCNVERTEX(I1)

              IF (I3.EQ.ASCNVERTEX(I1)) THEN
                I4=1
              ELSE
                I4=I3+1
              ENDIF
              IF ((DABS(P1(1,J)-ASCVERTEX(I3*2-1,I1)).LT.TOL.AND.
     .             DABS(P1(2,J)-ASCVERTEX(I3*2  ,I1)).LT.TOL.AND.
     .             DABS(P2(1,J)-ASCVERTEX(I4*2-1,I1)).LT.TOL.AND.
     .             DABS(P2(2,J)-ASCVERTEX(I4*2  ,I1)).LT.TOL).OR.
     .            (DABS(P2(1,J)-ASCVERTEX(I3*2-1,I1)).LT.TOL.AND.
     .             DABS(P2(2,J)-ASCVERTEX(I3*2  ,I1)).LT.TOL.AND.
     .             DABS(P1(1,J)-ASCVERTEX(I4*2-1,I1)).LT.TOL.AND.
     .             DABS(P1(2,J)-ASCVERTEX(I4*2  ,I1)).LT.TOL)) THEN

                IF (NCELL.LE.NOPTIM) THEN           
                  IF (NLIMPB >= NLIMPS) THEN
                    IGJUM3(NCELL,J)=0
                  ELSE
                    CALL BITSET (IGJUM3,0,NOPTIM,NCELL,J,0,NBITS)
                  ENDIF
                ENDIF
                ICOUNT=ICOUNT+1

c...            Add surface to short list of surfaces to check when in 
c               this cell:
                IGJUM4(NCELL,0)               =  IGJUM4(NCELL,0) + 1
                IGJUM4(NCELL,IGJUM4(NCELL,0)) =  J

              ENDIF

c            WRITE(6,'(A,3I6,2F14.8,2X,2F14.8,4I6)') '---->',
c     .        I1,j,I3,
c     .        P1(1,J),P1(2,J),
c     .        ASCVERTEX(I3*2-1,I1),ASCVERTEX(I3*2,I1),
c     .        ncell,j,nsurf,igjum3(ncell,j)
c            WRITE(6,'(A,18X,2F14.8,2X,2F14.8)') '---->',
c     .        P2(1,J),P2(2,J),
c     .        ASCVERTEX(I4*2-1,I1),ASCVERTEX(I4*2,I1)
         
            ENDDO

          ELSEIF (ASC3DMODE.EQ.1) THEN

            STOP 'NO LONGER IN USE'

c SEARCH EACH SIDE OF THE CUBE.  ASSUMPTIONS ARE MADE ABOUT THE ORDER
c OF THE VERTICIES BASED ON HOW THE 3D VACUUM MESH IS WRITTEN
c BY DIVIMP.
            DO I3=1,ASCNVERTEX(I1)
              IF (I3.EQ.ASCNVERTEX(I1)) THEN
                I4=1
              ELSE
                I4=I3+1
              ENDIF

              V1(1) = ASCVERTEX(2*I4-1,I1)
              V1(2) = ASCVERTEX(2*I4  ,I1)
              V1(3) = ASCZMIN3D(       I1)
              V2(1) = ASCVERTEX(2*I4-1,I1)
              V2(2) = ASCVERTEX(2*I4  ,I1)
              V2(3) = ASCZMAX3D(       I1)
              V4(1) = ASCVERTEX(2*I3-1,I1)
              V4(2) = ASCVERTEX(2*I3  ,I1)
              V4(3) = ASCZMAX3D(       I1)
              V3(1) = ASCVERTEX(2*I3-1,I1)
              V3(2) = ASCVERTEX(2*I3  ,I1)
              V3(3) = ASCZMIN3D(       I1)

              NMATCH = 0
              DO I5 = 1, 3
c                WRITE(0,'(A,1P,4E12.5)') 
c     .            'V:',V1(I5),V2(I5),V3(I5),V4(I5)
                IF ((DABS(P1(I5,J)-V1(I5)).LT.TOL.AND.
     .               DABS(P2(I5,J)-V2(I5)).LT.TOL.AND.
     .               DABS(P3(I5,J)-V3(I5)).LT.TOL.AND.
     .               DABS(P4(I5,J)-V4(I5)).LT.TOL).OR.
     .              (DABS(P1(I5,J)-V4(I5)).LT.TOL.AND.
     .               DABS(P2(I5,J)-V3(I5)).LT.TOL.AND.
     .               DABS(P3(I5,J)-V2(I5)).LT.TOL.AND.
     .               DABS(P4(I5,J)-V1(I5)).LT.TOL))  NMATCH = NMATCH + 1
              ENDDO

c MATCH FOUND FOR CELL SIDE
              IF (NMATCH.EQ.3) THEN
c                WRITE(0,*) 'MATCH 0:',i1,j

c...            Add surface to short list of surfaces to check when in 
c               this cell:
                IGJUM4(NCELL,0)               =  IGJUM4(NCELL,0) + 1
                IGJUM4(NCELL,IGJUM4(NCELL,0)) =  J

                JSCRTCH(J) = JSCRTCH(J) + 1

                IF (NCELL.LE.NOPTIM) THEN
                  IF (NLIMPB >= NLIMPS) THEN
                    IGJUM3(NCELL,J)=0
                  ELSE
                    CALL BITSET (IGJUM3,0,NOPTIM,NCELL,J,0,NBITS)
                  ENDIF
                ENDIF
c                IGJUM3(NCELL,J)=0

                ICOUNT=ICOUNT+1
              ENDIF

            ENDDO

c CELL END SURFACE 5:
            XMIN =  1.0D+20
            XMAX = -1.0D+20
            YMIN =  1.0D+20
            YMAX = -1.0D+20
            DO I3=1,ASCNVERTEX(I1)
              XMIN=MIN(XMIN,ASCVERTEX(2*I3-1,I1))
              XMAX=MAX(XMAX,ASCVERTEX(2*I3-1,I1))
              YMIN=MIN(YMIN,ASCVERTEX(2*I3  ,I1))
              YMAX=MAX(YMAX,ASCVERTEX(2*I3  ,I1))
            ENDDO

            V1(1) = XMAX
            V1(2) = YMIN
            V1(3) = ASCZMIN3D(I1)
            V2(1) = XMIN      
            V2(2) = YMIN      
            V2(3) = ASCZMIN3D(I1)
            V3(1) = XMAX      
            V3(2) = YMAX      
            V3(3) = ASCZMIN3D(I1)
            V4(1) = XMIN      
            V4(2) = YMAX      
            V4(3) = ASCZMIN3D(I1)
	    		      
            NMATCH = 0
            DO I5 = 1, 3
c              WRITE(6,'(A,1P,4E12.5,I6)') 
c     .          'P:',V1(I5,J),V2(I5,J),V3(I5,J),V4(I5,J),J
c              WRITE(6,'(A,1P,4E12.5)') 
c     .          'V:',V1(I5),V2(I5),V3(I5),V4(I5)
              IF (DABS(P1(I5,J)-V1(I5)).LT.TOL.AND.
     .            DABS(P2(I5,J)-V2(I5)).LT.TOL.AND.
     .            DABS(P3(I5,J)-V3(I5)).LT.TOL.AND.
     .            DABS(P4(I5,J)-V4(I5)).LT.TOL) NMATCH = NMATCH + 1
            ENDDO
            IF (NMATCH.EQ.3) THEN

c...          Add surface to short list of surfaces to check when in 
c             this cell:
              IGJUM4(NCELL,0)               =  IGJUM4(NCELL,0) + 1
              IGJUM4(NCELL,IGJUM4(NCELL,0)) =  J

              JSCRTCH(J) = JSCRTCH(J) + 1

c               WRITE(0,*) 'MATCH 5:',i1,j
              IF (NCELL.LE.NOPTIM) THEN
                IF (NLIMPB >= NLIMPS) THEN
                  IGJUM3(NCELL,J)=0
                ELSE
                  CALL BITSET (IGJUM3,0,NOPTIM,NCELL,J,0,NBITS)
                ENDIF
              ENDIF
c              IGJUM3(NCELL,J)=0
c              ICOUNT=ICOUNT+1
            ENDIF

            V1(3) = ASCZMAX3D(I1)
            V2(3) = ASCZMAX3D(I1)
            V4(3) = ASCZMAX3D(I1)
            V3(3) = ASCZMAX3D(I1)
	    		      
            NMATCH = 0
            DO I5 = 1, 3
c              WRITE(0,'(A,1P,4E12.5)') 
c     .          'V:',V1(I5),V2(I5),V3(I5),V4(I5)
              IF (DABS(P1(I5,J)-V1(I5)).LT.TOL.AND.
     .            DABS(P2(I5,J)-V2(I5)).LT.TOL.AND.
     .            DABS(P3(I5,J)-V3(I5)).LT.TOL.AND.
     .            DABS(P4(I5,J)-V4(I5)).LT.TOL) NMATCH = NMATCH + 1
            ENDDO
            IF (NMATCH.EQ.3) THEN

c...          Add surface to short list of surfaces to check when in 
c             this cell:
              IGJUM4(NCELL,0)               =  IGJUM4(NCELL,0) + 1
              IGJUM4(NCELL,IGJUM4(NCELL,0)) =  J

              JSCRTCH(J) = JSCRTCH(J) + 1

c               WRITE(0,*) 'MATCH 6:',i1,j
              IF (NCELL.LE.NOPTIM) THEN
                IF (NLIMPB >= NLIMPS) THEN
                  IGJUM3(NCELL,J)=0
                ELSE
                  CALL BITSET (IGJUM3,0,NOPTIM,NCELL,J,0,NBITS)
                ENDIF
              ENDIF
c              IGJUM3(NCELL,J)=0
c              ICOUNT=ICOUNT+1
            ENDIF

          ENDIF

        ENDDO

c        WRITE(6,*) 'ICOUNT=',ncell,icount,ascnvertex(i1)
      
        IF (ICOUNT.EQ.ASCNVERTEX(I1)) THEN
c        IF (ASC3DMODE.EQ.0.AND.ICOUNT.EQ.ASCNVERTEX(I1).OR.
c     .      ASC3DMODE.EQ.1.AND.ICOUNT.EQ.ASCNVERTEX(I1)) THEN
c     .      ASC3DMODE.EQ.1.AND.ICOUNT.EQ.ASCNVERTEX(I1)+2)) THEN
          DO J=1,SBGKI-1
            IF (P1(3,J).EQ.-1.0D+20.AND.P2(3,J).EQ.1.0D+20) THEN
              IF (NCELL.LE.NOPTIM) THEN
                IF (NLIMPB >= NLIMPS) THEN
                  IGJUM3(NCELL,J)=1
                ELSE
                  CALL BITSET (IGJUM3,0,NOPTIM,NCELL,J,1,NBITS)
                ENDIF
              ENDIF
c             IGJUM3(NCELL,J)=1
            ENDIF
          ENDDO
        ELSE
c        ELSEIF (ASC3DMODE.EQ.1.OR.ASC3DMODE.EQ.2) THEN
c...      Make sure that all the wall surfaces are checked:
          IGJUM4(NCELL,0)               =  IGJUM4(NCELL,0) + 1
          IGJUM4(NCELL,IGJUM4(NCELL,0)) =  -1
        ENDIF

c...    Add appropriate toroidal region boundary surfaces:
        IF (ASC3DMODE.EQ.2) THEN
          IF (NC.GT.1      ) THEN
            IGJUM4(NCELL,0)               =  IGJUM4(NCELL,0) + 1
            IGJUM4(NCELL,IGJUM4(NCELL,0)) =  CBGKI+NC-1
          ENDIF
          IF (NC.LT.ASCNCUT) THEN
            IGJUM4(NCELL,0)               =  IGJUM4(NCELL,0) + 1
            IGJUM4(NCELL,IGJUM4(NCELL,0)) =  CBGKI+NC
          ENDIF
        ENDIF                         

c...    Always check the manditory surfaces:      
c        IF (ASC3DMODE.EQ.1.OR.ASC3DMODE.EQ.2) THEN
          IF (HADDI.NE.EBGKI) THEN
            IGJUM4(NCELL,0)               =  IGJUM4(NCELL,0) + 1
            IGJUM4(NCELL,IGJUM4(NCELL,0)) =  -2
          ENDIF
          IGJUM4(NCELL,0)               =  IGJUM4(NCELL,0) + 1
          IGJUM4(NCELL,IGJUM4(NCELL,0)) =  -4
c        ENDIF
      ENDDO

c...  For additional cell 1:
      IGJUM4(NSURF+1,0) =  3
      IGJUM4(NSURF+1,1) = -1
      IGJUM4(NSURF+1,2) = -2
      IGJUM4(NSURF+1,3) = -4


c      DO NCELL=1,NSBOX
c        WRITE(6,'(A,3I6,3X,10I6)') 'MARK: IGJUM4=',
c     .    NCELL,NSURF,MAX(0,NCELL-NSURF),
c     .    (IGJUM4(NCELL,I1),I1=0,9)
c      ENDDO


c...  Dump IJGUM3 to a file:
c      WRITE(0,*) 'WRITING IGJUM3 TO FILE',ASCCODE
c      FP=98
c      OPEN(UNIT=FP,FILE='jummap.dat',FORM='UNFORMATTED',
c     .     STATUS='REPLACE',ERR=98)
c      WRITE(FP) ASCCODE,NLIMI/NOPTM1+1,MIN(NOPTIM,NSBOX),NSBOX
c      WRITE(FP) ((IGJUM3(J,I),I=1,NLIMI/NOPTM1+1),
c     .                        J=1,MIN(NOPTIM,NSBOX))
c      WRITE(FP) ((IGJUM4(J,I),I=0,8),J=1,NSBOX)
c      CLOSE(FP)
c      WRITE(0,*) 'DONE'

      RETURN
 98   WRITE(0,*) 'UNABLE TO OPEN IGJUM3 STORAGE FILE'
 99   STOP
      END
c
c ======================================================================
c
      INTEGER FUNCTION CheckCell(nr,np,xval,yval,cp)

      USE PARMMOD
      USE COMPRT
      USE CGRID
      USE CGEOM
      USE CLOGAU
      USE CPOLYG
      USE CCONA
      USE CLGIN
      USE CUPD
      USE CTRIG

      IMPLICIT NONE

c
c Input:
c
      INTEGER          nr,np
      REAL*8 xval,yval,cp(4)

      REAL*8 TOL
      PARAMETER (TOL = -1.0D-03)

c      REAL*8 vx(4),vy(4),ax,ay,bx,by,x,y
c sl*16
      REAL*8 vx(4),vy(4),ax,ay,bx,by,x,y,tx(4),ty(4)
      INTEGER          ii,i1,i2,icnt

      CheckCell = 0

      DO ii = 1,4
        cp(ii) = -2.0
      ENDDO
c
c
c
      IF (nr.LT.1.OR.nr.GE.nr1st) THEN
        CheckCell = -1
        RETURN
      ENDIF
c
c
c
      IF (np.LT.NPOINT(1,1).OR.np.EQ.NPOINT(2,1).OR.
     .    np.EQ.NPOINT(2,2).OR.np.GT.NPOINT(2,3)-1) THEN
        CheckCell = -1
        RETURN
      ENDIF
c
c
c
      IF (gridopt.EQ.1.AND.rvrtag(nr,np).EQ.1) RETURN

      IF (printopt.EQ.99)
     .   WRITE(6,'(A,2I4)') 'CHECKCELL: (NR,NP)',
     .    nr,np

      x = xval
      y = yval

      IF (gridopt.EQ.1) THEN
        vx(1) = xvert(nr,np  ,1)
        vx(2) = xvert(nr,np  ,2)
        vx(3) = xvert(nr,np+1,2)
        vx(4) = xvert(nr,np+1,1)

        vy(1) = yvert(nr,np  ,1)
        vy(2) = yvert(nr,np  ,2)
        vy(3) = yvert(nr,np+1,2)
        vy(4) = yvert(nr,np+1,1)
      ELSE
        STOP 'ERROR (CheckCell): Unsupported geometry option'
      ENDIF

      DO ii = 1,4
        cp(ii) = -1.0
      ENDDO

      DO ii = 1, 4
        IF (ii.EQ.4) THEN
          i1 = 4
          i2 = 1
        ELSE
          i1 = ii
          i2 = ii + 1
        ENDIF

        ax = x - vx(i1)
        ay = y - vy(i1)

        bx = vx(i2) - vx(i1)
        by = vy(i2) - vy(i1)

        cp(ii) = ax * by - ay * bx

        IF (bx.NE.0.0) THEN
          tx(ii) = (x - vx(i1)) / bx
        ELSE
          tx(ii) = 0.0
        ENDIF

        IF (by.NE.0.0) THEN
          ty(ii) = (y - vy(i1)) / by
        ELSE
          ty(ii) = 0.0
        ENDIF

        IF (printopt.EQ.99)
     .    WRITE(6,'(3X,3I2,2E10.3,6E10.3,2X,E11.4)')
     .      ii,i1,i2,vx(ii),vy(ii),x,y,ax,ay,bx,by,cp(ii)
      ENDDO

      CheckCell = 1

      icnt = 0

      DO ii = 1,4
        IF (cp(ii).LT.TOL) CheckCell = 0
        IF (cp(ii).LT.0.0) icnt = icnt + 1
      ENDDO

      IF (icnt.GT.1) CheckCell = 0

      IF (CheckCell.EQ.0.AND.icnt.EQ.1) THEN
        icnt = 0
        DO ii = 1, 4
          IF (tx(ii).GE.0.0.AND.tx(ii).LE.1.0.AND.
     .        abs(tx(ii)-ty(ii)).LT.TOL)
     .      icnt = icnt + 1

          IF (printopt.EQ.99)
     .      WRITE(6,'(A,1P,2F15.7)')
     .        'CHECKCELL: ',tx(ii),ty(ii)
        ENDDO

c Not sure if this is a good plan or not... this should not be necessary...
c I think the problem is in DIVIMP and in the grid refinement code...
        IF (icnt.EQ.1) THEN
          CheckCell = 1

          IF (printopt.EQ.99)
     .      WRITE(6,'(A,2I5)')
     .        'CHECKCELL: Found particle on surface ',np,nr
        ENDIF
      ENDIF

      RETURN
      END
c
c ======================================================================
c
      SUBROUTINE FindPoloidalCell(nr,x,y)

      USE PARMMOD
      USE COMPRT
      USE CUPD
      USE CGEOM
      USE CPOLYG

      IMPLICIT NONE

      INTEGER NR
      REAL*8 X, Y

      INTEGER PolPos,CheckCell,ip1,ip2,ii, ITEMP1, II1, II2
      REAL*8 cp(4)

        ip1 = PolPos(NR,X,Y,cp,'FOLNEUT 2')

        IF (printopt.GE.1.AND.printopt.LE.10) THEN
          WRITE(6,'(3X,A,I4,3X,A,I4,2F10.3)')
     .      'POLCELL: PolPos = ',ip1,
     .      ' (NR X,Y) ',nr,x,y
          WRITE(6,'(10X,4F14.8)') cp(1),cp(2),cp(3),cp(4)
        ENDIF



c        ip2 = PolPos(nr-nincx,x,Y,cp,'FOLNEUT 3')
c
c        IF (printopt.GE.1.AND.printopt.LE.10) THEN
c          WRITE(6,'(3X,A,I4,3X,A,I4,2F10.3)')
c     .      'FOLNEUT: PolPos = ',ip2,
c     .      ' (NR-NINCX x,Y) ',nr-nincx,x,y
c            WRITE(6,'(10X,4F14.8)') cp(1),cp(2),cp(3),cp(4)
c        ENDIF

        IF (ip1.NE.npcell) THEN
          IF (ip1.GT.0) THEN
            npcell = ip1
            ipolgn = ip1
          ELSEIF (ip1.LT.0) THEN
c            WRITE(0,*) 'ERROR (FindPolCell): Unable to find cell',
c     .        ' ( NR = ',nr,')'
            WRITE(6,*) 'ERROR (FindPolCell): Unable to find cell',
     .        ' ( NR = ',nr,')'

            ip1 = PolPos(nr,x,Y,cp,'FOLNEUT 4')

            itemp1 = printopt
            printopt = 99

c            DO ii1 = 1, npplg
c              DO ii2 = npoint(1,ii1),npoint(2,ii1)-1
c                ip1 = CheckCell(nr-nincx,ii2,x,y,cp)
c              ENDDO
c            ENDDO

            DO ii1 = 1, npplg
              DO ii2 = npoint(1,ii1),npoint(2,ii1)-1
                ip1 = CheckCell(nr,ii2,x,y,cp)
              ENDDO
            ENDDO

            printopt = itemp1

c            WRITE(0,*) 'NPANU = ',npanu
            WRITE(6,*) 'NPANU = ',npanu

            WRITE(6,*) 'PolPos routine - error'
            npcell = -1

c
c
c Can't use LEARC1 at the moment, because there is no differention
c of surface geometry based on which side the particle is on (or coming
c from or going to).  It is okay for particles outside the grid,
c becasue there is crude decision making in LEARC1 based on the ring index.
c
c
c            isl1 = -1
c            isl2 = -1
c            isl1 = LEARC1(x,y,0.0D0,isl2,nr,nr,
c     .                    .TRUE.,.FALSE.,
c     .                    npanu,'TEST 02     ')
c            WRITE(6,'(3X,A,I5)') 'FINDPOLOIDALCELL: LEARC1 = ',isl2
c            npcell = isl2

            RETURN
          ELSE
            WRITE(6,'(3X,A)') 'POLCELL:'
            WRITE(6,'(3X,A)') 'POLCELL: ERROR'
            WRITE(6,'(3X,A)') 'POLCELL:'
          ENDIF
        ENDIF

        IF (printopt.GE.1.AND.printopt.LE.10) THEN
          WRITE(6,'(3X,A)')
     +      'POLCELL: Output (NR,NP IR,IP MRSURF X0,Y0,Z0 X,Y,Z00)'
          WRITE(6,'(5X,A,2I4,A,2I4,A,I4,6F10.4)')
     +      '(',NR,NPCELL,') (',IRCELL,IPCELL,')',MRSURF,
     +      X0,Y0,Z0,x,y,Z00
        ENDIF

c        WRITE(6,'(10X,A,I4,2X,2I4,2X,2I4,2X,2F12.6)')
c     .    'NR NPCELL,IPCELL IP1,IP2 x,Y',
c     .    nr,npcell,ipcell,ip1,ip2,x,y

c      DO ii = 1, npplg
c        DO ip = npoint(1,ii), npoint(2,ii)
c          WRITE(6,'(3X,A,I4,3X,A,2I4,2F10.3)')
c     .      'FOLNEUT: CheckCell = ',CheckCell(nr,ipx00,y,cp),
c     .      ' (NR,IP X,Y) ',nr,ip,x,y
c          WRITE(6,'(10X,4F12.6)') cp(1),cp(2),cp(3),cp(4)
c         ENDDO
c      ENDDO


c        ip1=LEARC2(X00,Y00,NR,NPANU,'Check   1    ')
c        Ip2=LEARC2(X00,Y00,NR-nincx,NPANU,'Check   2    ')

c        WRITE(6,'(10X,A,I4,2X,2I4,2X,2I4,2X,2F12.6)')
c     .    'NRCELL NPCELL,IPCELL IP1,IP2 X00,Y00',
c     .    nrcell,npcell,ipcell,ip1,ip2,x00,y00




      RETURN
      END
c
c ======================================================================
c
      INTEGER FUNCTION PolPos(nr,x,y,cp,callid)

      USE PARMMOD
      USE COMPRT
      USE CGRID
      USE CGEOM
      USE CLOGAU
      USE CPOLYG
      USE CCONA
      USE CLGIN
      USE CUPD
      USE CTRIG

      IMPLICIT NONE

      INTEGER          nr
      REAL*8 x,y,cp(4)
      CHARACTER callid*(*)

      REAL*8 TOL
      PARAMETER (TOL = -1.0E-05)


      INTEGER CheckCell

      INTEGER          ii,i2,ip,status,cpminip(200),chkcell
      REAL*8 tx,ty,x1,y1,cp1(4),cpmin(200),cpmax
c
c Initialization:
c
      status = 0
      PolPos = -1
c
c Check if NR is outside standard grid:
c
      IF (nr.LT.1.OR.nr.GE.nr1st) THEN
        PolPos = 0
        RETURN
      ENDIF

      cp(1) = -1.0
      cp(2) = -1.0
      cp(3) = -1.0
      cp(4) = -1.0

      DO ii = 1, npplg
        DO ip = npoint(1,ii), npoint(2,ii)-1

           chkcell = CheckCell(nr,ip,x,y,cp1)

           IF (chkcell.EQ.1) THEN
             status = status + 1

             cp(1)  = cp1(1)
             cp(2)  = cp1(2)
             cp(3)  = cp1(3)
             cp(4)  = cp1(4)

             cpmin(status) = 1.0
             DO i2 = 1, 4
               cpmin(status) = MIN(cpmin(status),cp(i2))
             ENDDO
             cpminip(status) = ip

             PolPos = ip
           ELSEIF (chkcell.EQ.-1) THEN
             WRITE(0,*) 'ERROR (CheckCell): Violation ( NRCELL = ',
     .         nr,'  IP = ',ip,'  id = ',callid,')'
             STOP
           ENDIF
         ENDDO
      ENDDO

      IF (status.EQ.1) THEN
        IF (cpmin(1).LT.-1.0D-8) THEN
          IF (printopt.GE.1.AND.printopt.LE.10)
     .      WRITE(6,'(3X,A,1P,E14.7)')
     .        ' POLPOS: *** Questionable *** (CPMAX) ',cpmin(1)
        ENDIF
      ELSEIF (status.LT.10) THEN
        cpmax = -1.0
        DO ii = 1, status
          IF (cpmin(ii).GT.cpmax) THEN
            cpmax  = cpmin  (ii)
            PolPos = cpminip(ii)
          ENDIF




          IF (printopt.GE.1.AND.printopt.LE.10)
     .      WRITE(6,'(3X,A,I2,A,I3,A,E14.7,3A)')
     .        ' POLPOS: Ambiguity (NRCELL = ',
     .        nr,'  NPCELL = ',cpminip(ii),'  CPMIN = ',cpmin(ii),
     .        '  id = ',callid,')'
        ENDDO

        IF (cpmax.LT.TOL)
c        IF (printopt.GE.1.AND.printopt.LE.10.AND.cpmax.LT.-1.0D-8)
     .    WRITE(6,'(3X,A,2I5,1P,E15.7)')
     .      ' POLPOS: Ambiquity POLPOS,NR CPMAX =  ',PolPos,nr,cpmax


      ELSEIF (status.GE.10) THEN

         WRITE(6,'(3X,A)')
     +    ' POLPOS: Output (NR,NP IR,IP MRSURF X0,Y0,Z0 X00,Y00,Z00)'
         WRITE(6,'(5X,A,2I4,A,2I4,A,I4,6F10.4,1X,A)')
     +    '(',NRCELL,NPCELL,') (',IRCELL,IPCELL,')',MRSURF,
     +    X0,Y0,Z0,X00,Y00,Z00,callid

c
c
c
        status = 0
        printopt = 99

        DO ii = 1, npplg
          DO ip = npoint(1,ii), npoint(2,ii)-1
             chkcell = CheckCell(nr,ip,x,y,cp1)
             IF (chkcell.EQ.1) WRITE(6,*) '   MATCH'
           ENDDO
        ENDDO

c
c
c



        WRITE(0,*) 'NPANU = ',npanu
        STOP 'ERROR (PolPos): More than 3 cells identified'
      ENDIF

      RETURN
      END
c
c ======================================================================
c
      SUBROUTINE CHKSTD(Z02,CELL)

      USE PARMMOD
      USE CLOGAU
      USE CGRID
      USE CGEOM
      USE CCONA
      USE COMPRT
 
      IMPLICIT NONE

      INTEGER CELL
      REAL*8 Z02, Z03

      IF (NLTRA) THEN
        IF (NLTOR) THEN
          STOP 'NLTOR=.TRUE. NOT SUPPORTED IN CHKVAC'
        ELSE
c...      Set the toroidal direction comparison variable to PHI instead
c         of Z0:
          PHISEG=2.0D0*PIA/DBLE(NTTRA-1)
          Z03=(PHI+DBLE(NTRSEG-1)*PHISEG)/DEGRAD
c          WRITE(0,*) 'SUPERPHI=',Z03,ntrseg
c          WRITE(0,*) 'SUPERPHI=',phi,degrad,phiseg

        ENDIF
      ELSE
        Z03=Z02
      ENDIF

      DO CELL=1,EIRNSDTOR-1

        IF (Z03.GE.EIRSDTOR(CELL).AND.Z03.LT.EIRSDTOR(CELL+1)) EXIT

      ENDDO

      IF (CELL.EQ.EIRNSDTOR) THEN
        WRITE(6,*) 'CANNOT FIND TOROIDAL REGION ON STANDARD GRID',NPANU
        WRITE(0,*) 'CANNOT FIND TOROIDAL REGION ON STANDARD GRID',NPANU
        WRITE(0,*) 'Z0=',Z03

        DO CELL=1,EIRNSDTOR-1

          WRITE(0,*) EIRSDTOR(CELL),EIRSDTOR(CELL+1)

        ENDDO

        CELL=0
      ENDIF

     
      RETURN
99    STOP
      END
c
c ======================================================================
c
      LOGICAL FUNCTION CHKTRA(Z02,MSURF2)

      USE PARMMOD
      USE CLOGAU
      USE CGRID
      USE CGEOM
      USE CCONA
      USE COMPRT

      IMPLICIT NONE

      INTEGER MSURF2, I1
      REAL*8 Z02, Z03
     

      IF (NLTRA) THEN
        IF (NLTOR) THEN
          STOP 'NLTOR=.TRUE. NOT SUPPORTED IN CHKVAC'
        ELSE
c...      Set the toroidal direction comparison variable to PHI instead
c         of Z0:
          PHISEG=2.0D0*PIA/DBLE(NTTRA-1)
          Z03=(PHI+DBLE(NTRSEG-1)*PHISEG)/DEGRAD
c          PHISEG=2.0D0*PIA/DBLE(NTTRA-1)
c          Z03=(PHI+DBLE(NTRSEG-1)*PHISEG)/DEGRAD
          
        ENDIF
      ELSE
        Z03 = Z02
      ENDIF

      CHKTRA=.FALSE.

      I1=0
      DO WHILE(.NOT.CHKTRA.AND.I1.LT.EIRNTRANS)
        I1=I1+1
        IF (MSURF2.EQ.NLIM+IDNINT(EIRTRANS(I1,1)).OR.
     .      MSURF2.EQ.-1) THEN
          IF (Z03.GT.EIRTRANS(I1,2).AND.
     .        Z03.LT.EIRTRANS(I1,3)) THEN
            CHKTRA=.TRUE.
          ENDIF
        ENDIF

c          WRITE(0,'(A,3F12.4,I6,L6,2I6)') 
c     .      '-->',z03,EIRTRANS(I1,2),EIRTRANS(I1,3),
c     .      MSURF2,CHKTRA,I1,EIRNTRANS

      ENDDO
      RETURN
      END
c
c ======================================================================
c
      INTEGER FUNCTION CHKVAC(NEWCELL,OLDCELL,X0,Y0,Z0,ILEARN,NPANU)

      USE PARMMOD
      USE CGEOM

      IMPLICIT NONE

      INTEGER NEWCELL,OLDCELL,NPANU
      REAL*8  X0,Y0,Z0

      INTEGER I1,I2,I3,FP,INIT,INDEX,ILEARN,
     .        IMATCH,CUT1,OLDCELL2,ICOUNT,LASTCL,CELLREGION
      LOGICAL NEW,OLD,OUTPUT1,SWITCH,STATUS
      REAL*8  X1,X2,TX,Y1,Y2,TY,XMIN,XMAX,YMIN,YMAX, TOR

      REAL*8     TOL1,TOL2

      INTEGER :: ISMART(0:1000,5)

      DATA INIT,ICOUNT /0,0/
      SAVE

      OLDCELL2 = OLDCELL

      ICOUNT = ICOUNT + 1

      LASTCL = -1

      IF (INIT.EQ.0) THEN
        INIT=1
        ISMART=0
        WRITE(0,*) 'ZEROING ISMART:',nlim
      ENDIF

      IF (ILEARN.LT.1.OR.ILEARN.GT.4) THEN      
        WRITE(0,*) 'LEARNING INDEX OUT OF RANGE',ILEARN
        STOP 'STOP IN CHKVAC'
      ENDIF

      IF (ISMART(0,ILEARN  ).EQ.NLIM.OR.
     .    ISMART(0,ILEARN+1).EQ.NLIM) THEN
        WRITE(0,*) 'ISMART INDEX OUT OF BOUNDS'
        STOP 'STOP IN CHKVAC'
      ENDIF

c      TOL1 = 2.0D-06
      TOL1 = 2.0D-05
c      TOL2 = 1.0D-05
      TOL2 = 1.0D-04

      NEW = .FALSE.
      OLD = .FALSE.

c  SOME OPTIMIZATION REQUIRED!

      output1 = .FALSE.
c...OUTPUT ON ERRORS IS TURNED OFF!
c      IF (optvac.EQ.2) output1 = .TRUE.

      CELLREGION=-1

      CHKVAC = -1
      IMATCH = -1

      STATUS = .TRUE.

      fp = 6

      TOR=Z0

c      IF (NLTRA) THEN
c        IF (NLTOR) THEN
c          STOP 'NLTOR=.TRUE. NOT SUPPORTED IN CHKVAC'
c        ELSE
c...      Set the toroidal direction comparison variable to PHI instead
c         of Z0:
c          PHISEG=2.0D0*PIA/DBLE(NTTRA-1)
c          TOR=(PHI+DBLE(NTRSEG-1)*PHISEG)*RADDEG
c          WRITE(0,*) 'TOR:',torc
c        ENDIF
c      ENDIF

c...  Find which toroidal block the particle is in, and adjust OLDCELL2
c     accordingly:
       IF (ASC3DMODE.EQ.2.AND.OLDCELL.GT.-1) THEN
        DO CUT1=1,ASCNCUT
          IF (TOR.GE.ASCZMIN3D(CUT1).AND.
     .        TOR.LT.ASCZMAX3D(CUT1)) THEN
            OLDCELL2 = OLDCELL2 - (CUT1 - 1) * ASCNCELL
            EXIT
          ENDIF
        ENDDO
      ELSEIF (OLDCELL.LT.-1) THEN
        OLDCELL2=-1
C RESTRICT ADDITIONAL CELL SEARCH TO
C VACUUM GRID REGION -OLDCELL+1
        CELLREGION=-OLDCELL-1
      ENDIF

c      IF (npanu.EQ.821) THEN
c        DO i1 = 1, ISMART(0,ILEARN)
c          WRITE(0,*) 'SMART:',ilearn,i1,ISMART(i1,ILEARN)
c        ENDDO
c      ENDIF


      DO WHILE(STATUS) 

        IF (output1) WRITE(fp,*) 'CHKVAC: ',tol2

        DO INDEX = -ISMART(0,ILEARN)+1, ASCNCELL

C...      Learning opportunity:
          IF (INDEX.LE.0) THEN
c            WRITE(0,*) 'ISMART   : ',ilearn,-INDEX+1,
c     .        ISMART(-INDEX+1,ILEARN),ismart(0,ilearn)
            I1 = ISMART(-INDEX+1,ILEARN)
          ELSE
            I1 = INDEX
          ENDIF

          IF (output1) WRITE(fp,'(A,4I6,2F10.5)') 
     .      'CHKVAC: ',I1,newcell,oldcell,chkvac,asczmin3d(i1),
     .      asczmax3d(i1)

          IF (CELLREGION.NE.-1.AND.ASCREGION(I1).NE.CELLREGION) CYCLE

          DO I2 = 1, ASCNVERTEX(I1)
            IF (I2.EQ.ASCNVERTEX(I1)) THEN
              I3 = 1
            ELSE
              I3 = I2 + 1
            ENDIF
            X1 = ASCVERTEX(2*I2-1,I1)
            Y1 = ASCVERTEX(2*I2  ,I1)
            X2 = ASCVERTEX(2*I3-1,I1)
            Y2 = ASCVERTEX(2*I3  ,I1)

c            IF (i1.EQ.ASCNCELL.and.ilearn.eq.4) THEN
c              WRITE(0,'(A,4F12.4)') '-SEARCH->',x1,y1,x2,y2
c            ENDIF

            SWITCH = .FALSE.
            IF (DABS(X1-X2).LT.1.0D-08) THEN
              IF (DABS(X0-X1).LT.TOL1) THEN
                TX = +99.0D0
              ELSE
                TX = -99.0D0
              ENDIF
            ELSE
              TX = (X0 - X1) / (X2 - X1)
              IF (TX.LT.0.25) THEN
                TX = (X0 - X2) / (X1 - X2)
                SWITCH = .TRUE.
              ENDIF
            ENDIF
            IF (DABS(Y1-Y2).LT.1.0D-08) THEN
c            IF (Y1.EQ.Y2) THEN
              IF (DABS(Y0-Y1).LT.TOL1) THEN
                TY = +99.0D0
              ELSE
                TY = -99.0D0
              ENDIF
            ELSE
              IF (SWITCH) THEN
                TY = (Y0 - Y2) / (Y1 - Y2)
              ELSE
                TY = (Y0 - Y1) / (Y2 - Y1)
              ENDIF
            ENDIF

            IF (((DABS(TX-TY).LT.TOL2.AND.
     .           TX.GE.0.0D0.AND.TX.LE.1.0D0.AND.
     .           TY.GE.0.0D0.AND.TY.LE.1.0D0).OR.
     .          (TX.EQ.+99.0D0.AND.TY.GE.0.0D0.AND.TY.LE.1.0D0).OR.
     .          (TY.EQ.+99.0D0.AND.TX.GE.0.0D0.AND.TX.LE.1.0D0)).AND.
     .          (ASC3DMODE.EQ.0.OR.ASC3DMODE.EQ.2.OR.
     .           (TOR.GT.ASCZMIN3D(I1).AND.TOR.LT.ASCZMAX3D(I1)))) THEN



              IF     (NEWCELL.EQ. 0) THEN
c...            Find additional cell that particle is in:
                CHKVAC = ASCCELL(I1)   
                IMATCH = I1
c...bug?
              ELSEIF (NEWCELL.EQ.-1.AND.OLDCELL2.NE.ASCCELL(I1)) THEN
c              ELSEIF (NEWCELL.EQ.-1.AND.OLDCELL.NE.ASCCELL(I1)) THEN
c...            Find additional cell that particle is in:
                IF     (CHKVAC.EQ.-1) THEN
                  CHKVAC = ASCCELL(I1)   
                  IMATCH = I1
                  LASTCL = CHKVAC
                ELSEIF (ASCCELL(I1).NE.LASTCL) THEN
c                ELSEIF (CHKVAC.NE.ASCCELL(I1)) THEN
                  CHKVAC = -2
                  IMATCH = -2
                ENDIF
              ENDIF


c              IF (npanu.EQ.821) THEN
c                WRITE(0,'(A,5I6,4F10.6,5I6,F10.7)') 
c     .            '-->',npanu,oldcell-1,oldcell2-1,asccell(i1)-1,
c     .                  i2,x0,x1,x2,TOR,icount,chkvac,lastcl,index,i1,
c     .                  tol2
c              ENDIF

c...          Find which toroidal block the particle is in, and increment
c             cell index accordingly:
c...bug?
              IF (ASC3DMODE.EQ.2.AND.CHKVAC.GE.1.AND.
     .            OLDCELL2.NE.ASCCELL(I1)) THEN
c              IF (ASC3DMODE.EQ.2.AND.CHKVAC.GE.1) THEN
                DO CUT1=1,ASCNCUT
                  IF (TOR.GE.ASCZMIN3D(CUT1).AND.
     .                TOR.LT.ASCZMAX3D(CUT1)) THEN
                    CHKVAC = CHKVAC + ASCNCELL * (CUT1 - 1)
                    EXIT
                  ENDIF
                ENDDO
              ENDIF

            ENDIF

            IF (output1) 
     .        WRITE(fp,'(A,2I6,2F12.6,2X,3F10.5,I6,4F12.6)') 
     .          '   ->',I1,I2,TX,TY,X0,Y0,TOR,CHKVAC,X1,X2,Y1,Y2

          ENDDO

c...      Forced exit if using learning curve:
c          IF (INDEX.LE.0) WRITE(0,*) 'TRYING:',i1,index,imatch

          IF (INDEX.LE.0.AND.CHKVAC.GT.0) THEN
c            IF (npanu.EQ.821) WRITE(0,*) 'TRYING TO BOOT',INDEX
            EXIT
          ENDIF

        ENDDO

c       IF (npanu.EQ.821) WRITE(0,*) 'BOOTED',index,i1

        IF (CHKVAC.LT.0.AND.TOL2.LE.1.0D-02) THEN
          TOL2 = TOL2 * 10.0D0
        ELSE

c...      Add cell to learning curve:
          IF (INDEX.GT.0.AND.imatch.GT.0) THEN
c...        About 10% performance improvement...or so...
            ISMART(0               ,ILEARN) = ISMART(0,ILEARN) + 1
            ISMART(ISMART(0,ILEARN),ILEARN) = imatch
c            WRITE(0,*) 'ISMART NEW:',ilearn,INDEX,imatch,
c     .                 ISMART(0,ILEARN)
          ENDIF

          STATUS = .FALSE.
        ENDIF

      ENDDO



c...  Trying to save the day by looking to see if the particle was launched
c     inside an additional cell, rather than on its boundary.  This can 
c     happen if a polygon surface is not aligned perfectly with the
c     the vacuum grid:


      IF (ILEARN.EQ.4.AND.CHKVAC.LT.0) THEN      



        CHKVAC = -1
        DO INDEX = -ISMART(0,ILEARN+1)+1, ASCNCELL
          IF (INDEX.LE.0) THEN
            I1 = ISMART(-INDEX+1,ILEARN+1)
          ELSE
            I1 = INDEX
          ENDIF
          IF (CELLREGION.NE.-1.AND.ASCREGION(I1).NE.CELLREGION) CYCLE
          XMIN =  1.0D+10
          XMAX = -1.0D+10
          YMIN =  1.0D+10
          YMAX = -1.0D+10
          DO I2 = 1, ASCNVERTEX(I1)
            XMIN = MIN(XMIN,ASCVERTEX(2*I2-1,I1))
            XMAX = MAX(XMAX,ASCVERTEX(2*I2-1,I1))
            YMIN = MIN(YMIN,ASCVERTEX(2*I2  ,I1))
            YMAX = MAX(YMAX,ASCVERTEX(2*I2  ,I1))
          ENDDO

          IF (X0.GT.XMIN.AND.X0.LT.XMAX.AND.
     .        Y0.GT.YMIN.AND.Y0.LT.YMAX) THEN
            IF (CHKVAC.EQ.-1) THEN
              CHKVAC = ASCCELL(I1)

              IF (INDEX.GT.0) THEN
                ISMART(0                 ,ILEARN+1)=ISMART(0,ILEARN+1)+1
                ISMART(ISMART(0,ILEARN+1),ILEARN+1)=I1
                WRITE(0,*) 'NEW DESPIRATE ISMART:',
     .            ILEARN+1,ISMART(0,ILEARN+1)
              ENDIF
            ELSE
              CHKVAC = -2
            ENDIF
          ENDIF
c...logic is flawed here:
          IF (INDEX.LE.0.AND.CHKVAC.GE.1) EXIT
        ENDDO
        IF (ASC3DMODE.EQ.2.AND.CHKVAC.GE.1) THEN
          DO CUT1=1,ASCNCUT
            IF (TOR.GE.ASCZMIN3D(CUT1).AND.
     .          TOR.LT.ASCZMAX3D(CUT1)) THEN
              CHKVAC = CHKVAC + ASCNCELL * (CUT1 - 1)
              EXIT
            ENDIF
          ENDDO
        ENDIF

      ENDIF



      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: CHECKNUM
c
      SUBROUTINE CHECKNUM(TAG,I1,I2,X1,X2)

      USE PARMMOD

      IMPLICIT NONE

      INTEGER :: I1, I2, CHKCNT, CHKERR
      REAL*8 X1, X2

      CHARACTER TAG*(*)



      RETURN


c      IF (printopt.GE.1.AND.printopt.LE.10) THEN
c        WRITE(6,'(10X,2A,2I4,2G15.7,$)')
c     +    'CHECKNUM: ',TAG,I1,I2,X1,X2
c      ENDIF

      CHKCNT = CHKCNT + 1

      IF (X1.EQ.X2) THEN

c        IF (printopt.GE.1.AND.printopt.LE.10) WRITE(6,*) ' '
c        IF (i1.EQ.25.AND.printopt.GE.1.AND.printopt.LE.10) THEN
c          WRITE(6,'(10X,A,A6,2I4,2G15.7)')
c     +      'CHECKNUM: ',TAG,I1,I2,X1,X2
c        ENDIF



      ELSE
c        IF (printopt.EQ.2.OR.debugopt.EQ.-2) THEN
c          WRITE(6,'(10X,2A,2I4,2G15.7,A)')
c     +      'CHECKNUM: ',TAG,I1,I2,X1,X2,' ERROR'
c        ENDIF

          WRITE(6,'(10X,2A,2I4,2G15.7,A)')
     +      'CHECKNUM: ',TAG,I1,I2,X1,X2,' ERROR'

        CHKERR = CHKERR + 1
      ENDIF

      RETURN
      END
c
c ======================================================================
c
c subroutine: MODUSR
c
      SUBROUTINE MODUSR

      USE PARMMOD
      USE COMUSR
      USE CPOLYG
      USE CCOUPL
      USE CGEOM
      USE CGRID

      IMPLICIT NONE

      INTEGER IT,IX,IY,IN,IFL,IPLS,NRED,NDX2

      REAL*8 DUMMY(0:NDXP,0:NDYP),SNID(0:NDXP,0:NDYP,NFL)

      IF (NITER.GT.0) CALL MODBGK

c      CALL ASDUSR

c      WRITE(0,*) 'NOTHING IN MODUSR'

c       ===============================================================
c
c       Output data for BGK iterations.  Some assumptions about the
c       number of ion species is being made, namely that there are 6:
c       D+, C+, D, D2, DD2, D2D.  Only the last 4 are to be saved here
c       and subsequently read by DIVIMP (they will not be processed to
c       the degree that the other PIN quantities are, in order to keep
c       things simple).
c
c
c...use this as the flag to identify when BGK is required -- good?
        IF (NITER.GT.0.AND.IITER.EQ.NITER) THEN

          IF (output) WRITE(0,*) 'WRITING BGK BACKGROUND PLASMA DATA'

          DO IT=1,NBMLT

c
c           ION DENSITY
c
            IF (output) WRITE(0 ,*) '[BGK: ION DENSITY]'

            SNID=0.0D0
            DO IFL=1,NFLA
              DO IPLS=1,NPLSI
                IF (IFLB(IPLS).NE.IFL.AND.IFLB(IPLS)+30.NE.IFL) CYCLE
                IF (output) WRITE(0,*) 'IFL,IPLS= ',IFL,IPLS
                DO IX=1,NDXA
                  DO IY=1,NDYA
                    IN=IY+(IX-1)*NR1ST+(IT-1)*NSTRD
                    SNID(IX,IY,IFL)=DIIN(IPLS,IN)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
            IF (NCUTL.NE.NCUTB)
     .        CALL INDMPI (SNID,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .                     NCUTB,NCUTL,NPOINT,NPPLG,1,1)

            NRED=(NPPLG-1)*(NCUTL-NCUTB)
            NDX2=NDXA-NRED
            WRITE(32,'(A,I3,A1)') '[BGK: ION DENSITY, REGION ',IT,']'
            CALL NEUTR(32,NDX2,NDYA,NFLA,SNID,NDX,NDY,NFL,1,1)
c
c           ION TEMPERATURE
c
            IF (output) WRITE(0 ,*) '[BGK: ION TEMPERATURE]'

            SNID=0.0D0
            DO IFL=1,NFLA
              DO IPLS=1,NPLSI
                IF (IFLB(IPLS).NE.IFL.AND.IFLB(IPLS)+30.NE.IFL) CYCLE
                DO IX=1,NDXA
                  DO IY=1,NDYA
                    IN=IY+(IX-1)*NR1ST+(IT-1)*NSTRD
                    SNID(IX,IY,IFL)=TIIN(IPLS,IN)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
            IF (NCUTL.NE.NCUTB)
     .        CALL INDMPI (SNID,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .                     NCUTB,NCUTL,NPOINT,NPPLG,1,1)
            NRED=(NPPLG-1)*(NCUTL-NCUTB)
            NDX2=NDXA-NRED
            WRITE(32,'(A,I3,A)') '[BGK: ION TEMPERATURE, REGION ',IT,']'
            CALL NEUTR(32,NDX2,NDYA,NFLA,SNID,NDX,NDY,NFL,1,1)
c
c           X-VELOCITY
c
            IF (output) WRITE(0 ,*) '[BGK: X-VELOCITY]'

            SNID=0.0D0
            DO IFL=1,NFLA
              DO IPLS=1,NPLSI
                IF (IFLB(IPLS).NE.IFL.AND.IFLB(IPLS)+30.NE.IFL) CYCLE
                DO IX=1,NDXA
                  DO IY=1,NDYA
                    IN=IY+(IX-1)*NR1ST+(IT-1)*NSTRD
                    SNID(IX,IY,IFL)=VXIN(IPLS,IN)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
            IF (NCUTL.NE.NCUTB)
     .        CALL INDMPI (SNID,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .                     NCUTB,NCUTL,NPOINT,NPPLG,1,1)
            NRED=(NPPLG-1)*(NCUTL-NCUTB)
            NDX2=NDXA-NRED
            WRITE(32,'(A,I3,A)') '[BGK: X-VELOCITY, REGION ',IT,']'
            CALL NEUTR(32,NDX2,NDYA,NFLA,SNID,NDX,NDY,NFL,1,1)
c
c           Y-VELOCITY
c
            IF (output) WRITE(0 ,*) '[BGK: Y-VELOCITY]'

            SNID=0.0D0
            DO IFL=1,NFLA
              DO IPLS=1,NPLSI
                IF (IFLB(IPLS).NE.IFL.AND.IFLB(IPLS)+30.NE.IFL) CYCLE
                DO IX=1,NDXA
                  DO IY=1,NDYA
                    IN=IY+(IX-1)*NR1ST+(IT-1)*NSTRD
                    SNID(IX,IY,IFL)=VYIN(IPLS,IN)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
           IF (NCUTL.NE.NCUTB)
     .        CALL INDMPI (SNID,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .                     NCUTB,NCUTL,NPOINT,NPPLG,1,1)
            NRED=(NPPLG-1)*(NCUTL-NCUTB)
            NDX2=NDXA-NRED
            WRITE(32,'(A,I3,A)') '[BGK: Y-VELOCITY, REGION ',IT,']'
            CALL NEUTR(32,NDX2,NDYA,NFLA,SNID,NDX,NDY,NFL,1,1)
c
c           Z-VELOCITY
c
            IF (output) WRITE(0 ,*) '[BGK: Z-VELOCITY]'

            SNID=0.0D0
            DO IFL=1,NFLA
              DO IPLS=1,NPLSI
                IF (IFLB(IPLS).NE.IFL.AND.IFLB(IPLS)+30.NE.IFL) CYCLE
                DO IX=1,NDXA
                  DO IY=1,NDYA
                    IN=IY+(IX-1)*NR1ST+(IT-1)*NSTRD
                    SNID(IX,IY,IFL)=VZIN(IPLS,IN)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
            IF (NCUTL.NE.NCUTB)
     .        CALL INDMPI (SNID,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .                     NCUTB,NCUTL,NPOINT,NPPLG,1,1)
            NRED=(NPPLG-1)*(NCUTL-NCUTB)
            NDX2=NDXA-NRED
            WRITE(32,'(A,I3,A)') '[BGK: Z-VELOCITY, REGION ',IT,']'
            CALL NEUTR(32,NDX2,NDYA,NFLA,SNID,NDX,NDY,NFL,1,1)

          ENDDO

          IF (output) WRITE(0,*) 'DONE'

        ELSE

          IF (output) WRITE(0,*) 'NOT WRITING BGK DATA TO DIVIMP'

        ENDIF
c
c      IF (NFILEL.EQ.1) THEN
c        WRITE(0,*) 'MODIFYING NFILEL'
c        WRITE(6,*) 'MODIFYING NFILEL'
c
c        NFILEL = 2
c      ENDIF
c slmod end
      RETURN
      END

c ======================================================================
c
c subroutine: ASDUSR
c
c Outputs additional surface data.
c
c
c
c     DIIN
c     TIIN
c     VXIN
c     VYIN
c     VZIN
c
c     PDENA Neutral atom density
c     EDENA    "     "   energy
c     PDENM Neutral molecule density
c     EDENM    "       "     energy
c     PDENI Ionised molecule density
c     EDENI    "       "     energy
c
c
c
      SUBROUTINE ASDUSR

      USE PRECISION
      USE PARMMOD
      USE COMSOU
      USE COMUSR
      USE COUTAU
      USE CESTIM
      USE CGRID
 
      IMPLICIT NONE

      INTEGER FP,IN,COUNT,I1,I2,ISTR,ISTRA,IATM,IMOL

      DATA COUNT /0/
      SAVE

      COUNT = COUNT + 1

      IF (output)
     .WRITE(0,*) 'WRITING ASD'

      FP=98
      OPEN(UNIT=FP,FILE='addsur.dat',ACCESS='SEQUENTIAL',
     .     STATUS='UNKNOWN',POSITION='APPEND')

      WRITE(FP,90) '*'
      WRITE(FP,94) '* DATA FOR ADDITIONAL EIRENE SURFACES ',COUNT
      WRITE(FP,90) '*                                     '
      WRITE(FP,90) '* NRADD (NUMBER OF SURFACES)=         '
      WRITE(FP,92) NRADD

C SURFACE DATA
      WRITE(FP,90) '* SURFACE DATA '
      WRITE(FP,91) 0,'IN','VOLUME (CM-3)'
      DO IN=NSURF+1,NSURF+NRADD
        WRITE(FP,92) IN-NSURF,IN,VOL(IN)
      ENDDO
C NEUTRAL ATOM DATA

      IF (NATMI.GT.2) THEN
        WRITE(0,*) 'ONLY PASSING 1ST AND 2ND ATOMIC INDEX TO DIVIMP'
        WRITE(6,*) 'ONLY PASSING 1ST AND 2ND ATOMIC INDEX TO DIVIMP'
      ENDIF

      WRITE(FP,90) '* NEUTRAL ATOM SPECIES       '
      WRITE(FP,90) '* NATMI (NUMBER OF SPECIES)= '
      WRITE(FP,92) MIN(2,NATMI)
      DO I1=1,MIN(2,NATMI)
        WRITE(FP,91) I1,'IN','PDENA (    )','EDENA (   )'
        DO IN=NSURF+1,NSURF+NRADD
          WRITE(FP,92) IN-NSURF,IN,PDENA(I1,IN),EDENA(I1,IN)
        ENDDO
      ENDDO
C NEUTRAL MOLECULE DATA
      WRITE(FP,90) '* NEUTRAL MOLECULE SPECIES   '
      WRITE(FP,90) '* NMOLI (NUMBER OF SPECIES)= '
      WRITE(FP,92) NMOLI
      DO I1=1,NMOLI
        WRITE(FP,91) I1,'IN','PDENM (    )','EDENM (   )'
        DO IN=NSURF+1,NSURF+NRADD
c...bug - Nov 22, 1999
          WRITE(FP,92) IN-NSURF,IN,PDENM(I1,IN),EDENM(I1,IN)
c          WRITE(FP,92) IN-NSURF,IN,PDENA(I1,IN),EDENA(I1,IN)
        ENDDO
      ENDDO
C TEST ION DATA
      WRITE(FP,90) '* TEST ION SPECIES           '
      WRITE(FP,90) '* NIONI (NUMBER OF SPECIES)= '
      WRITE(FP,92) NIONI
      DO I1=1,NIONI
        WRITE(FP,91) I1,'IN','PDENI (    )','EDENI (   )'
        DO IN=NSURF+1,NSURF+NRADD
          WRITE(FP,92) IN-NSURF,IN,PDENI(I1,IN),EDENI(I1,IN)
        ENDDO
      ENDDO
C BULK ION DATA (FOR BGK ROUTINES)
      WRITE(FP,90) '* BULK ION SPECIES           '
      WRITE(FP,90) '* NPLSI (NUMBER OF SPECIES)= '
      WRITE(FP,92) NPLSI
      DO I1=1,NPLSI
        WRITE(FP,91) I1,'IN','DIIN','TIIN','VXIN','VYIN','VZIN','EDRIFT'
        DO IN=NSURF+1,NSURF+NRADD
          WRITE(FP,92) IN-NSURF,IN,DIIN(I1,IN),TIIN(I1,IN),
     .                 VXIN(I1,IN),VYIN(I1,IN),VZIN(I1,IN),EDRIFT(I1,IN)
        ENDDO
      ENDDO

c RESIDUALS CALCULATED IN MODBGK ROUTINE
      IF (NITER.GE.1) THEN
        WRITE(FP,90) 'BGK RESIDUALS'
        WRITE(FP,*) NPLSI
        DO I1 = 1, NPLSI
          WRITE(FP,93) I1,(EIRRES(I2,I1),I2=1,7)
        ENDDO
      ENDIF

C STRATA DATA
      WRITE(FP,90) 'STRATUM DATA FOR SELECTED ADDITIONAL CELLS'
      WRITE(FP,*) NSTRDAT,NSTRAI,NATMI,NMOLI
      DO ISTR=1,NSTRDAT
        I1=STRDAT(ISTR)
        WRITE(FP,*) STRDAT(ISTR)
        DO ISTRA=1,NSTRAI
          DO IATM=1,NATMI
            WRITE(FP,'(2I6,1P,2E12.4)') ISTRA,IATM, 
     .        PSTRDATA(IATM,ISTR,ISTRA),ESTRDATA(IATM,ISTR,ISTRA)
          ENDDO
          DO IMOL=1,NMOLI
            WRITE(FP,'(2I6,1P,2E12.4)') ISTRA,IMOL, 
     .        PSTRDATM(IMOL,ISTR,ISTRA),ESTRDATM(IMOL,ISTR,ISTRA)
          ENDDO
        ENDDO
c...    Totals:
        DO IATM=1,NATMI
          WRITE(FP,'(2I6,1P,2E12.4)') 0,IATM, 
     .      PDENA(IATM,NSURF+I1),EDENA(IATM,NSURF+I1)
        ENDDO
        DO IMOL=1,NMOLI
          WRITE(FP,'(2I6,1P,2E12.4)') 0,IMOL, 
     .      PDENM(IMOL,NSURF+I1),EDENM(IMOL,NSURF+I1)
        ENDDO
      ENDDO

c...  SNEAK IN THE FLUX FOR EACH STRATUM AS WELL
      WRITE(FP,'(1P,100(E12.4:)') (FLUXT(ISTRA),ISTRA=1,NSTRAI)



c...TEMP
      WRITE(6,*) 'TESTING STRATA DATA'
      DO ISTR=1,NSTRDAT
        I1=STRDAT(ISTR)
        WRITE(6,*)
        WRITE(6,*)  ISTR,I1,':'

        DO ISTRA=1,NSTRAI
          DO IATM=1,NATMI
            WRITE(6,'(2I6,1P,2E12.4)') ISTRA,IATM, 
     .        PSTRDATA(IATM,ISTR,ISTRA),
     .        ESTRDATA(IATM,ISTR,ISTRA)
          ENDDO
          DO IMOL=1,NMOLI
            WRITE(6,'(2I6,1P,2E12.4)') ISTRA,IMOL, 
     .        PSTRDATM(IMOL,ISTR,ISTRA),
     .        ESTRDATM(IMOL,ISTR,ISTRA)
          ENDDO
        ENDDO
 
        WRITE(6,*)
        DO IATM=1,NATMI
          WRITE(6,'(I6,1P,2E12.4)') IATM, 
     .      PDENA(IATM,NSURF+I1),
     .      EDENA(IATM,NSURF+I1)
        ENDDO
        DO IMOL=1,NMOLI
          WRITE(6,'(I6,1P,2E12.4)') IMOL, 
     .      PDENM(IMOL,NSURF+I1),
     .      EDENM(IMOL,NSURF+I1)
        ENDDO
      ENDDO

      CLOSE (FP)

      RETURN
90    FORMAT(A)
91    FORMAT(I3,A5,7X,20(A15:))
92    FORMAT(3X,I5,I7,1P,20(E15.7:),0P,I6)
93    FORMAT(1P,I4,7E12.4,0P)
94    FORMAT(A,I6)
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
      SUBROUTINE ReadData(zeile)

c put in subroutine... add a check to see if *** 0 exists...
c
      USE PRECISION
      USE PARMMOD
      USE CGEOM
      USE COMPRT

      IMPLICIT NONE

      INTEGER i1,i2
      CHARACTER zeile*(*)


      geomopt  = 0
      gridopt  = 0
      addopt   = 0
      neutopt  = 0
      debugopt = 0
      gchkopt  = 0
      printopt = 0
      optuser  = 0
      optzmotion = 0

      eirntrans = 0


       WRITE(6,'(1X,A72)') ZEILE

      WRITE(6,*) 'READING OPTIONS FOR EIRENE02'

        READ (IUNIN,*) comment,geomopt
        READ (IUNIN,*) comment,gridopt
c        READ (IUNIN,*) comment,addopt
        READ (IUNIN,*) comment,neutopt
        READ (IUNIN,*) comment,debugopt
        READ (IUNIN,*) comment,cxd2opt      
        READ (IUNIN,*) comment,optuser      
        READ (IUNIN,*) comment,eirntrans
        IF (eirntrans.GT.0) THEN
          DO i1 = 1, eirntrans
            READ(IUNIN,*) (eirtrans(i1,i2),i2=1,3)
          ENDDO
        ENDIF
        READ (IUNIN,*) comment,nstrdat
        READ (IUNIN,'(6I6)') (strdat(i1),i1=1,nstrdat)        
        READ (IUNIN,*) comment,volcor2

       IF (debugopt.GE.1.AND.debugopt.LE.20) THEN
         printopt = debugopt
       ELSEIF (debugopt.EQ.-41) THEN
         optzmotion = 1
         debugopt = 0
       ELSEIF (debugopt.EQ.-42) THEN
         optzmotion = 2
         debugopt = 0
       ELSEIF (debugopt.EQ.-43) THEN
         optzmotion = 3
         debugopt = 0
       ELSEIF (debugopt.EQ.-44) THEN
         optzmotion = 4
         debugopt = 0
       ELSEIF (debugopt.EQ.-20) THEN
         opttest=-(debugopt+19)
         debugopt=0
         WRITE(0,*) 'OPTTEST= ',opttest
       ELSE
         printopt = 0
       ENDIF

      IF (debugopt.NE.0.OR.optzmotion.NE.0) THEN
        WRITE(0,*)
        WRITE(0,*) '=================================='
        WRITE(0,*)
        WRITE(0,'(1X,A,I6)') 'grid option            ',gridopt
        WRITE(0,'(1X,A,I6)') 'geometry source option ',geomopt
        WRITE(0,'(1X,A,I6)') 'geometry check option  ',gchkopt
c        WRITE(0,'(1X,A,I6)') 'addusr option          ',addopt
        WRITE(0,'(1X,A,I6)') 'neutral wall option    ',neutopt
        WRITE(0,'(1X,A,I6)') 'debug option           ',debugopt
        WRITE(0,'(1X,A,I6)') 'print option           ',printopt
        WRITE(0,'(1X,A,I6)') 'CX D2+ production      ',cxd2opt      
        WRITE(0,'(1X,A,I6)') 'Z motion               ',optzmotion
        WRITE(0,'(1X,A,I6)')
        WRITE(0,*) '=================================='
        WRITE(0,*)
      ENDIF

c Read next 2 lines:
c       READ (IUNIN,'(A72)') ZEILE
c       WRITE(6,'(1X,A72)') ZEILE
c       READ (IUNIN,'(A72)') ZEILE


       RETURN
       END



c
c ======================================================================
c
c
c
c
c
      SUBROUTINE GEODIV(NDXA,NDYA,NPLP,NR1ST,
     .                   PUX,PUY,PVX,PVY,NBMLT)


      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CGEOM
      USE CLOGAU
      IMPLICIT NONE

       

      INTEGER, INTENT(INOUT) :: NDXA,NDYA,NPLP,NR1ST,NBMLT
      REAL(DP), INTENT(OUT) :: PUX(*),PUY(*),PVX(*),PVY(*)


      CHARACTER*80 LINE,LINE2
      INTEGER DIMXH,DIMYH,NNCUT,NXCUT1,NXCUT2,d1,d2,d3,dummy
      CHARACTER c1*12
      REAL*8 MERK(NDY)
      INTEGER k


      INTEGER fp,i1,i2,idum1,np

C  GEOMETRY DATA: CELL VERTICES (LINDA ---> EIRENE)
      REAL(DP) ::
     R  X1(NDX),Y1(NDX),X2(NDX),Y2(NDX),X3(NDX),Y3(NDX),
     R  X4(NDX),Y4(NDX)
      INTEGER :: IX, IY, J, I



C  ACTUAL MESH USED IN THIS RUN
C
       REWIND 30

3366  FORMAT(/)



c Temporary...
c       TRCHST = .FALSE.


c
c Read file header:
      IF (debugopt.NE.0)
     .  WRITE(0,*) 'GEODIV: Reading geometry file'

      READ(30,'(A60)'   ) comment
      READ(30,*)
      READ(30,'(A30,I6)') comment,dimxh
      READ(30,'(A30,I6)') comment,dimyh
      READ(30,'(A30,I6)') comment,nncut
      READ(30,'(A30,I6)') comment,nxcut1
      READ(30,'(A30,I6)') comment,nxcut2
      READ(30,'(A30,I6)') comment,dummy
      READ(30,'(A30,I6)') comment,dummy
      READ(30,'(A30,I6)') comment,dummy
      READ(30,*)


      NPLP = 3


      NPOINT(1,1) = 1
      NPOINT(2,1) = nxcut1
      NPOINT(1,2) = nxcut1 + 1
      NPOINT(2,2) = nxcut2 + 1
      NPOINT(1,3) = nxcut2 + 2
      NPOINT(2,3) = dimxh + 3
c
c     Read in geometry data:
c
      DO ix= 1,dimxh+3
        DO iy = 1, dimyh+1

          READ(30,*,END=25,ERR=20)
     +      xvert(iy,ix,1),yvert(iy,ix,1),
     +      xvert(iy,ix,2),yvert(iy,ix,2),
     +      d1,d2

c          READ(30,3334,END=25,ERR=20)
c     +      xvert(iy,ix,1),yvert(iy,ix,1),
c     +      xvert(iy,ix,2),yvert(iy,ix,2),
c     +      d1,d2,c1

            xpol(iy,ix) = xvert(iy,ix,1)
            ypol(iy,ix) = yvert(iy,ix,1)
        ENDDO
      ENDDO

3334  FORMAT(4D15.7,5X,2I6,A2)
      goto 30

20    CONTINUE
      write(0,*) 'GEODIV: Error reading geometry file',
     .  ' (IX,IY) ',ix,iy
      GOTO 30

25    CONTINUE
      write(0,*) 'GEODIV: Error EOF while reading geometry file',
     .  ' (IX,IY) ',ix,iy
30    CONTINUE


      IF (GRIDOPT.EQ.1) THEN
c
c Calculate PVRTAG:
c
        DO iy = 1, dimyh+1
          DO ix = 1, dimxh+3 - 1
              pvrtag(iy,ix) = 0
              rvrtag(iy,ix) = 0
          ENDDO

          DO ix = NPOINT(1,1),NPOINT(2,1)-1
            IF (xvert(iy,ix,1).EQ.xvert(iy,ix+1,1).AND.
     .          yvert(iy,ix,1).EQ.yvert(iy,ix+1,1).AND.
     .          xvert(iy,ix,2).EQ.xvert(iy,ix+1,2).AND.
     .          yvert(iy,ix,2).EQ.yvert(iy,ix+1,2)) THEN
              pvrtag(iy,ix+1) = 1
            ENDIF
          ENDDO

          DO ix = NPOINT(1,3),NPOINT(2,3)-1
            IF (xvert(iy,ix,1).EQ.xvert(iy,ix+1,1).AND.
     .          yvert(iy,ix,1).EQ.yvert(iy,ix+1,1).AND.
     .          xvert(iy,ix,2).EQ.xvert(iy,ix+1,2).AND.
     .          yvert(iy,ix,2).EQ.yvert(iy,ix+1,2)) THEN
              pvrtag(iy,ix) = 1
            ENDIF
          ENDDO

        ENDDO
c
c Calcule RVRTAG:
c
        DO iy = 1, dimyh+1
          DO ix = NPOINT(1,1),NPOINT(2,1)-1
            IF (xvert(iy,ix,1).EQ.xvert(iy,ix+1,1).AND.
     .          yvert(iy,ix,1).EQ.yvert(iy,ix+1,1).AND.
     .          xvert(iy,ix,2).EQ.xvert(iy,ix+1,2).AND.
     .          yvert(iy,ix,2).EQ.yvert(iy,ix+1,2)) THEN
              rvrtag(iy,ix) = 1
            ENDIF
          ENDDO

          DO ix = NPOINT(1,3),NPOINT(2,3)-1
            IF (xvert(iy,ix,1).EQ.xvert(iy,ix+1,1).AND.
     .          yvert(iy,ix,1).EQ.yvert(iy,ix+1,1).AND.
     .          xvert(iy,ix,2).EQ.xvert(iy,ix+1,2).AND.
     .          yvert(iy,ix,2).EQ.yvert(iy,ix+1,2)) THEN
              rvrtag(iy,ix) = 1
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      WRITE(60,*)
      WRITE(60,*) 'Mesh projection:'
      WRITE(60,*)

      DO 1015 IY=1,NDYA
        WRITE(60,*)
        DO 1014 IX=1,NDXA
          IF (gridopt.EQ.1) THEN
            X1(IX)=XVERT(IY,IX  ,1)
            Y1(IX)=YVERT(IY,IX  ,1)
            X2(IX)=XVERT(IY,IX+1,1)
            Y2(IX)=YVERT(IY,IX+1,1)
            X3(IX)=XVERT(IY,IX  ,2)
            Y3(IX)=YVERT(IY,IX  ,2)
            X4(IX)=XVERT(IY,IX+1,2)
            Y4(IX)=YVERT(IY,IX+1,2)
          ELSE

          ENDIF

          WRITE(60,'(2I4,8E11.4)')
     +     iy,ix,
     +     x1(ix),y1(ix),x2(ix),y2(ix),x3(ix),y3(ix),x4(ix),y4(ix)

1014    CONTINUE
        CALL MSHPROJ (X1,Y1,X2,Y2,X3,Y3,X4,Y4,PUX,PUY,PVX,PVY,NDXA,
     .                NR1ST,IY)
1015  CONTINUE
C
c
c       Output geometry data for debugging:
c
        WRITE(60,*)
        WRITE(60,*) 'Geometry data:'

        DO ix = 1, dimyh+1
          WRITE(60,*)
          DO iy = 1, ndya+1
            WRITE(60,'(4F13.8,1X,2I3,5X,2I4,2X,2I4)')
     +      xvert(iy,ix,1),yvert(iy,ix,1),
     +      xvert(iy,ix,2),yvert(iy,ix,2),
     +      pvrtag(iy,ix),rvrtag(iy,ix),iy,ix,iy+1,ix
          ENDDO
        ENDDO


        WRITE(60,*)
        WRITE(60,*) 'Geometry data (DIVIMP style):'

        DO iy = 1, ndya+1
         WRITE(60,*)
         DO ix = 1, dimxh+3
            WRITE(60,'(4F13.8,1X,2I3,5X,2I4,2X,2I4)')
     +      xvert(iy,ix,1),yvert(iy,ix,1),
     +      xvert(iy,ix,2),yvert(iy,ix,2),
     +      pvrtag(iy,ix),rvrtag(iy,ix),ix,iy,ix,iy+1
          ENDDO
        ENDDO




      NP=NPOINT(2,NPLP)
      DO 1020 J=1,NDYA+1
        DO 1020 I=1,NP
          XPOL(J,I)=XPOL(J,I)*100.
          YPOL(J,I)=YPOL(J,I)*100.

          IF (gridopt.EQ.1) THEN
            DO k = 1, 2
              xvert(j,i,k) = xvert(j,i,k) * 100.0
              yvert(j,i,k) = yvert(j,i,k) * 100.0

c              IF (ABS(xvert(j,i,k)).LT.5.D-5) xvert(j,i,k) = 0.0
c              IF (ABS(yvert(j,i,k)).LT.5.D-5) yvert(j,i,k) = 0.0
            ENDDO
          ENDIF

1020  CONTINUE


      close(60)

      WRITE(61,*) 'GEOMETRY DATA - GEODIV'
      WRITE(61,*)
      WRITE(61,*) 'ndxa   ',ndxa
      WRITE(61,*) 'ndya   ',ndya
      WRITE(61,*) 'nplp   ',nplp
      WRITE(61,*)
      WRITE(61,*) '1: ',NPOINT(1,1),NPOINT(2,1)
      WRITE(61,*) '2: ',NPOINT(1,2),NPOINT(2,2)
      WRITE(61,*) '3: ',NPOINT(1,3),NPOINT(2,3)

      DO ix = 1, NPOINT(2,nplp)
         WRITE(61,*)
         DO iy = 1, ndya+1
          WRITE(61,'(2F10.4,2I6)') xpol(iy,ix),ypol(iy,ix),iy,ix
        ENDDO
      ENDDO
      CLOSE(61)

c...  READ VACUUM GRID POLYGONS:

      optvac = 0

      fp = 98
      OPEN(UNIT=fp,FILE='vacuum_grid.dat',ACCESS='SEQUENTIAL',
     .     STATUS='OLD',ERR=99)

      optvac = 1

      READ(fp,*) 
      READ(fp,*) ascncell,idum1,asc3dmode,asccode,ascncut
c      WRITE(0,*) 'MARK: ASC3DMODE=',asc3dmode,ascncell,idum1
      DO i1 = 1, ascncell       
        READ(fp,*) idum1,asccell(i1),ascnvertex(i1),ascregion(i1)
        READ(fp,*) asczmin3D(i1),asczmax3D(i1)
        DO i2 = 1, ascnvertex(i1)
          READ(fp,*) ascvertex(2*i2-1,i1),ascvertex(2*i2,i1)
        ENDDO
c        IF (i1.EQ.ascncell) THEN
c          DO i2 = 1, ascnvertex(i1)
c            WRITE(6,*) '-->',ascvertex(2*i2-1,i1),ascvertex(2*i2,i1)
c          ENDDO
c        ENDIF
      ENDDO

c...  Attempt to read in the location of the toroidal surfaces
c     that do NBLOCK switching:

      IF (NLMLT) THEN
        READ(fp,*) 
        READ(fp,*) EIRNSDTOR
        DO I1=1,EIRNSDTOR
          READ(fp,*) EIRSDTOR(I1)
        ENDDO

        WRITE(6,*) 'TOROIDAL NBLOCK SWITHCING SURFACES'
c        WRITE(0,*) 'TOROIDAL NBLOCK SWITHCING SURFACES'
        DO I1=1,EIRNSDTOR
          WRITE(6,*) I1,EIRSDTOR(I1)
c          WRITE(0,*) I1,EIRSDTOR(I1)
        ENDDO

        IF (EIRNSDTOR.NE.NBMLT+1) THEN
          WRITE(0,*) 'WARNING: EIRNSDTOR.NE.NBMLT+1'
          WRITE(6,*) 'WARNING: EIRNSDTOR.NE.NBMLT+1'
        ENDIF
      ENDIF

      CLOSE(FP)

      IF (EIRNSDTOR.GT.0.AND.EIRNTRANS.EQ.0) THEN
        WRITE(0,*) 'ERROR: OLD TRANSPAREND ND-SURFACE SYSTEM'
        WRITE(6,*) 'ERROR: OLD TRANSPAREND ND-SURFACE SYSTEM'
        CALL EXIT
      ENDIF


c      WRITE(6,*)
c      WRITE(6,*) 'VACUUM GRID POLYGON DATA'
c      WRITE(6,*) '  NUMBER OF CELLS= ',ascncell
c      DO i1 = 1, ascncell       
c        WRITE(6,*) i1,asccell(i1),ascnvertex(i1)
c        DO i2 = 1, ascnvertex(i1)
c          WRITE(6,'(2F14.7)') ascvertex(2*i2-1,i1),ascvertex(2*i2,i1)
c        ENDDO
c      ENDDO



99    CONTINUE

      IF (optvac.GT.0) THEN
        WRITE(6,*) 
        WRITE(6,*) '********************'
        WRITE(6,*) 
        WRITE(6,*) 'OPTVAC= ',optvac
        WRITE(6,*) 
        WRITE(6,*) '********************'
        WRITE(6,*) 
      ENDIF


      IF (debugopt.NE.0)
     .  WRITE(0,*) 'GEODIV: Done'

      RETURN



*//END GEOMD//
      END







c slmod end

      SUBROUTINE BROAD_USR
      IMPLICIT NONE
      RETURN
      END
C
*//GEOMD//
C=======================================================================
C          S U B R O U T I N E   G E O M D
C=======================================================================
      SUBROUTINE GEOMD(NDXA,NDYA,NPLP,NR1ST,
     .                 PUX,PUY,PVX,PVY)
C
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CGEOM
      IMPLICIT NONE
C
      REAL(DP), INTENT(OUT) :: PUX(*),PUY(*),PVX(*),PVY(*)
      INTEGER, INTENT(INOUT) :: NDXA,NDYA,NPLP,NR1ST

      character(110) :: zeile
      REAL(DP) :: br(0:ndxp,0:ndyp,4),bz(0:ndxp,0:ndyp,4)
C
C  GEOMETRY DATA: CELL VERTICES (LINDA ---> EIRENE)
      REAL(DP) ::
     R  X1(NDX),Y1(NDX),X2(NDX),Y2(NDX),X3(NDX),Y3(NDX),
     R  X4(NDX),Y4(NDX)
      REAL(DP) :: DX, DY
      INTEGER :: I0, I0E, IX, IY, I1, I2, I3, I4, IPART, I, J

      ndxa=0
      ndya=0

      read (30,*)
      read (30,*)
      read (30,*)
      read (30,*)

1     continue
      read (30,'(a110)',end=99) zeile
      i0=index(zeile,'(')
      i0e=index(zeile,')')
      read (zeile(i0+1:i0e-1),*) ix,iy
      ndxa=max(ndxa,ix)
      ndya=max(ndya,iy)
      i1=index(zeile,': (')
      i2=index(zeile(i1+3:),')')+i1+2
      read (zeile(i1+3:i2-1),*) br(ix,iy,4),bz(ix,iy,4)
      i3=index(zeile(i2+1:),'(')+i2
      i4=i3+index(zeile(i3+1:),')')
      read (zeile(i3+1:i4-1),*) br(ix,iy,3),bz(ix,iy,3)

      read (30,'(a110)') zeile

      read (30,'(a110)') zeile
      i1=index(zeile,'(')
      i2=index(zeile,')')
      read (zeile(i1+1:i2-1),*) br(ix,iy,1),bz(ix,iy,1)
      i3=i2+index(zeile(i2+1:),'(')
      i4=i2+index(zeile(i2+1:),')')
      read (zeile(i3+1:i4-1),*) br(ix,iy,2),bz(ix,iy,2)

      read (30,*)
      goto 1


99    continue
      ndxa=ndxa-1
      ndya=ndya-1
C
      DO 1015 IY=1,NDYA
        DO 1014 IX=1,NDXA
          X1(IX)=br(ix,iy,1)
          Y1(IX)=bz(ix,iy,1)
          X2(IX)=br(ix,iy,2)
          Y2(IX)=bz(ix,iy,2)
          X3(IX)=br(ix,iy,4)
          Y3(IX)=bz(ix,iy,4)
          X4(IX)=br(ix,iy,3)
          Y4(IX)=bz(ix,iy,3)
1014    CONTINUE
        CALL MSHPROJ (X1,Y1,X2,Y2,X3,Y3,X4,Y4,PUX,PUY,PVX,PVY,NDXA,
     .                NR1ST,IY)
1015  CONTINUE
C
C SEARCH FOR THE CUTS
C
      IPART=1
      NPOINT(1,IPART)=1
      IY=1
      DO IX=1,NDXA
        DX=BR(IX+1,IY,1)-BR(IX,IY,2)
        DY=BZ(IX+1,IY,1)-BZ(IX,IY,2)
        IF (DX*DX+DY*DY.GT.EPS10) THEN
C CUT GEFUNDEN
          NPOINT(2,IPART)=IX+1+IPART-1
          IPART=IPART+1
          NPOINT(1,IPART)=IX+1+IPART-1
        ENDIF
      ENDDO
      NPOINT(2,IPART)=NDXA+IPART
C
      NPLP=IPART
C
      DO IY=1,NDYA
        DO IPART=1,NPLP
          DO IX=NPOINT(1,IPART),NPOINT(2,IPART)-1
            XPOL(IY,IX)=BR(IX-(IPART-1),IY,1)
            YPOL(IY,IX)=BZ(IX-(IPART-1),IY,1)
          ENDDO
          XPOL(IY,NPOINT(2,IPART))=BR(NPOINT(2,IPART)-IPART,IY,2)
          YPOL(IY,NPOINT(2,IPART))=BZ(NPOINT(2,IPART)-IPART,IY,2)
        ENDDO
      ENDDO
C INTRODUCE OUTERMOST RADIAL POLYGON
      DO IPART=1,NPLP
        DO IX=NPOINT(1,IPART),NPOINT(2,IPART)-1
          XPOL(NDYA+1,IX)=BR(IX-(IPART-1),NDYA,4)
          YPOL(NDYA+1,IX)=BZ(IX-(IPART-1),NDYA,4)
        ENDDO
        XPOL(NDYA+1,NPOINT(2,IPART))=BR(NPOINT(2,IPART)-IPART,NDYA,3)
        YPOL(NDYA+1,NPOINT(2,IPART))=BZ(NPOINT(2,IPART)-IPART,NDYA,3)
      ENDDO
C
      DO J=1,NDYA+1
        DO I=1,NPOINT(2,NPLP)
          XPOL(J,I)=XPOL(J,I)*100.
          YPOL(J,I)=YPOL(J,I)*100.
        END DO
      END DO
C
      ndxa=npoint(2,nplp)-1
c
C     do j=1,ndya+1
C       write (6,*)
C       write (6,*) 'in geomd polygon ',j
C       write (6,'(1p,6e12.4)') (xpol(j,i),ypol(j,i),i=1,npoint(2,nplp))
C     enddo
C
      RETURN
*//END GEOMD//
      END
C
C
C
      SUBROUTINE GEOUSR
C
C   PREPARE DATA FOR LIMITER-SURFACES
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CADGEO
      USE COMPRT
      USE CTRCEI
      USE CCONA
      USE CGEOM
      USE CGRID
      USE CLGIN
      USE CINIT
      USE CPOLYG
      IMPLICIT NONE
      CHARACTER(80) :: ZEILE
      REAL(DP) :: XCOOR, YCOOR, ZCOOR, xan, xen, yan, yen, zan, zen 
      INTEGER :: NADMOD, NASMOD, NORMOD, NRS, IPUNKT, I,NSSIR, NSSIP,
     .           IDIR, IR, IP, IT, IC, IN, NAS
      logical :: lx1,lx2,lx3,lx4,ly1,ly2,ly3,ly4
C
C MODIFY GEOMETRY
C
C
c slmod begin 
      IF (neutopt.EQ.1) RETURN
c
c It is unnecessary to alter geometry data in this
c routine because...
c
c
c slmod end
      READ (IUNIN,'(A80)') ZEILE
      READ (IUNIN,'(3I6)') NADMOD,NASMOD,NORMOD

      DO I=1,NADMOD
        READ (IUNIN,'(2I6,3E12.4)') NRS,IPUNKT,XCOOR,YCOOR,ZCOOR

        SELECT CASE(IPUNKT)

        CASE DEFAULT
           WRITE (6,*) 'WRONG POINTNUMBER IN ADDUSR '
           WRITE (6,*) 'INPUT LINE READING'
           WRITE (6,'(2I6,1P,3E12.4)') NRS,IPUNKT,XCOOR,YCOOR,ZCOOR
           WRITE (6,*) ' IS IGNORED '

        CASE (1)
           P1(1,NRS)=XCOOR
           P1(2,NRS)=YCOOR
           P1(3,NRS)=ZCOOR

        CASE (2)
           P2(1,NRS)=XCOOR
           P2(2,NRS)=YCOOR
           P2(3,NRS)=ZCOOR

        CASE (3)
           P3(1,NRS)=XCOOR
           P3(2,NRS)=YCOOR
           P3(3,NRS)=ZCOOR

        CASE (4)
           P4(1,NRS)=XCOOR
           P4(2,NRS)=YCOOR
           P4(3,NRS)=ZCOOR

        CASE (5)
           P5(1,NRS)=XCOOR
           P5(2,NRS)=YCOOR
           P5(3,NRS)=ZCOOR

        CASE (6)
           P6(1,NRS)=XCOOR
           P6(2,NRS)=YCOOR
           P6(3,NRS)=ZCOOR

        END SELECT
      ENDDO

      DO I=1,NASMOD
        READ (IUNIN,'(5I6)') NAS,IPUNKT,NSSIR,NSSIP
        IF (IPUNKT.EQ.1) THEN
          P1(1,NAS)=XPOL(NSSIR,NSSIP)
          P1(2,NAS)=YPOL(NSSIR,NSSIP)
        ELSEIF (IPUNKT.EQ.2) THEN
          P2(1,NAS)=XPOL(NSSIR,NSSIP)
          P2(2,NAS)=YPOL(NSSIR,NSSIP)
        ELSE
          WRITE (6,*) 'WRONG POINTNUMBER IN ADDUSR '
          WRITE (6,*) 'INPUT LINE READING'
          WRITE (6,'(5I6)') NAS,IPUNKT,NSSIR,NSSIP
          WRITE (6,*) ' IS IGNORED '
        ENDIF
      ENDDO

      DO I = 1, NORMOD
        READ (IUNIN,'(5I6)') IDIR,IR,IP
        IF (IDIR == 1) THEN
          PLNX(IR,IP) = -PLNX(IR,IP)
          PLNY(IR,IP) = -PLNY(IR,IP)
        ELSE IF (IDIR == 2) THEN
          PPLNX(IR,IP) = -PPLNX(IR,IP)
          PPLNY(IR,IP) = -PPLNY(IR,IP)
        ELSE
          WRITE (6,*) ' IDIR =',IDIR,' NOT FORESEEN IN GEOUSR '
          WRITE (6,*) IDIR,IR,IP
          WRITE (6,*) ' IS IGNORED '
        END IF
      END DO
C
C
C  ABSCHALTEN NICHT ERREICHBARER ODER DOPPELT VORHANDENER FLAECHEN
C
C   HIERHER: LGJUM1, LGJUM2 SETZEN ZUR BESCHLEUNIGUNG (NICHT UNBEDINGT
C   NOETIG)
C
C   LGJUM1(J,I)=.TRUE. :
C   ABSCHALTEN DER FLAECHE I, FALLS TEILCHEN AUF J SITZT
C
C   LGJUM2(J,I)=.TRUE. :
C   ABSCHALTEN DES ERSTEN SCHNITTPUNKTES MIT FLAECHE I, FALLS
C   TEILCHEN AUF J SITZT (FALLS I EINE FLAECHE ZWEITER ORDNUNG IST)
C
C   DEFAULTS: LGJUM1(J,J)=.TRUE. FUER EBENE FLAECHEN,
C             LGJUM2(J,J)=.TRUE. FUER FLAECHEN ZWEITER ORDNUNG
C
C
C
C  SET SOME VOLUMES EXPLIZIT
C
C
C  MODIFY REFLECTION MODEL AT TARGET PLATES
C
C     do i=1,nlimps
C       do isp=1,natmi+nmoli+nioni
C         recyct(isp,i)=1.
C       enddo
C     enddo

      RETURN
      END


      SUBROUTINE iniUSR
      IMPLICIT NONE
      RETURN
      END


      FUNCTION LEAUSR(A,B,C)
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: A, B, C
      INTEGER :: LEAUSR
      LEAUSR=1
      RETURN
      END
c
c
c slmod begin
c      subroutine modusr
c      return
c      end
c slmod end
C
C
      SUBROUTINE MSHADJ (X1,Y1,X2,Y2,XPLG,YPLG,NPLP,NPOINT,M1,NDX,IR)
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: X1(*),Y1(*),X2(*),Y2(*)
      INTEGER, INTENT(IN) ::M1, NDX, IR
      REAL(DP), INTENT(INOUT) :: XPLG(M1,*),YPLG(M1,*)
      INTEGER, INTENT(INOUT) :: NPOINT(2,*), NPLP
      REAL(DP) :: EPS, D, D1
      INTEGER :: NPRT, NP, I
      LOGICAL :: LDUMCL,LPRT
C  LPRT=.TRUE.  : VALID PART
C  LDUMCL=.TRUE.: DUMMY CELL
C
      EPS=1.E-20
      NPRT=0
      NP=0
      LDUMCL=.FALSE.
      LPRT=.FALSE.
C
      DO 1 I=1,NDX
        D=ABS((X1(I)-X2(I))**2+(Y1(I)-Y2(I))**2)
        IF (D.LE.EPS) THEN
          IF (LPRT) THEN
C  ENDE DER VALID ZELLEN
            NPOINT(2,NPRT)=NP
            LPRT=.FALSE.
            LDUMCL=.TRUE.
          ENDIF
        ELSEIF (D.GT.EPS) THEN
C  STARTPUNKT DIESES POLYGONS = ERSTER VALID PUNKT?
          IF (.NOT.LPRT.AND..NOT.LDUMCL) THEN
            LPRT=.TRUE.
            NPRT=NPRT+1
            NP=NP+1
            XPLG(IR,NP)=X1(I)
            YPLG(IR,NP)=Y1(I)
            NPOINT(1,NPRT)=NP
          ELSEIF (.NOT.LPRT.AND.LDUMCL) THEN
            D1=ABS((X1(I+1)-X2(I+1))**2+(Y1(I+1)-Y2(I+1))**2)
            IF (D1.GT.EPS) THEN
C  ENDE DER DUMMYZELLEN: SUCHE NAECHSTE ZELLE MIT D1.GT.0.
              LDUMCL=.FALSE.
              LPRT=.TRUE.
              NPRT=NPRT+1
              NPOINT(1,NPRT)=NP
            ENDIF
          ENDIF
        ENDIF
        NP=NP+1
        XPLG(IR,NP)=X2(I)
        YPLG(IR,NP)=Y2(I)
1     CONTINUE
      NPOINT(2,NPRT)=NP
      NPLP=NPRT
      RETURN
      END


      SUBROUTINE PLAUSR
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CSTEP
      USE CGRID
      USE CINIT
      USE COMSOU
      IMPLICIT NONE
      REAL(DP) :: FACTOR, STEP
      INTEGER :: IG, IR, IT, IN, IP, IPLS, NRWL
      REAL(DP) :: XMACH


      RETURN
      END
C
C
      SUBROUTINE PLTUSR(PLABLE,J)
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CADGEO
      USE CCONA
      USE CLGIN
      IMPLICIT NONE
      LOGICAL, INTENT(INOUT) :: PLABLE
      INTEGER, INTENT(IN) :: J
      INTEGER :: IB, MERK, IER
      REAL(DP) :: XS, YS, ZX0, ZY0, ZZ0, CX, CY, CZ, RZYLB, B0B, B1B,
     .          B2B, B3B, F0B, F1B, F2B, F3B, X0, Y0, B0, B1, B2, Z1, Z2
      INTEGER :: IN
      RETURN
      END
CDK USER
C
C   USER SUPPLIED SUBROUTINES
C
C           ************
C           *  TEXTOR  *  (HELIUM INCLUDED)
C           ************
C
      SUBROUTINE PROUSR (PRO,INDX,P0,P1,P2,P3,P4,P5,PROVAC,N)
C
C     P0 : CENTRAL VALUE
C     P1 : STARTING RADIUS FOR POLYNOMIAL
C     P2 : SWITCH FOR PHASE 1: OH-PHASE
C                           2: NI-PHASE
C     P3 : FACTOR FOR TI: TI=K*TE
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CGRID
      USE CGEOM
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: P0, P1, P2, P3, P4, P5, PROVAC
      REAL(DP), INTENT(OUT) :: PRO(*)
      INTEGER, INTENT(IN) :: INDX, N

      RETURN
      END


      SUBROUTINE REFUSR
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: XMW,XCW,XMP,XCP,ZCOS,ZSIN,EXPI,RPROB,
     .                        E0TERM
      INTEGER, INTENT(IN) :: IGASF,IGAST
      ENTRY RF0USR
      ENTRY SPTUSR
      ENTRY SP0USR
      ENTRY SP1USR
      ENTRY RF1USR (XMW,XCW,XMP,XCP,IGASF,IGAST,ZCOS,ZSIN,EXPI,
     .              RPROB,E0TERM,*,*,*,*)
      RETURN
      END
c
c
      subroutine retusr(sig)
      USE PRECISION
      implicit none
      real(dp), intent(in) :: sig
      return
      end
C
C
      SUBROUTINE SAMUSR (NLSF,X0,Y0,Z0,
     .              SORAD1,SORAD2,SORAD3,SORAD4,SORAD5,SORAD6,
     .              IRUSR,IPUSR,ITUSR,IAUSR,IBUSR,
     .              TEWL,TIWL,DIWL,VXWL,VYWL,VZWL,WEISPZ)
C
C  SAMPLE INITAL COORDIANTES X,Y,Z ON ADDITIONAL SURFACE NLLI
C
      USE PRECISION
      USE PARMMOD
      USE CADGEO
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: SORAD1,SORAD2,SORAD3,SORAD4,SORAD5,SORAD6
      REAL(DP), INTENT(OUT) :: X0,Y0,Z0,TEWL,TIWL,DIWL,VXWL,VYWL,VZWL,
     .                       WEISPZ
      INTEGER, INTENT(IN) :: NLSF,is1, is2
      INTEGER, INTENT(OUT) :: IRUSR, IPUSR, ITUSR, IAUSR, IBUSR
      REAL(DP) :: X, Y, T, B0, B1, B2, Z1, Z2, RANF
      INTEGER :: IER

      entry sm0usr (is1,is2,sorad1,sorad2,sorad3,sorad4,sorad5,sorad6)
      return

      entry SM1USR (NLSF,X0,Y0,Z0,
     .              SORAD1,SORAD2,SORAD3,SORAD4,SORAD5,SORAD6,
     .              IRUSR,IPUSR,ITUSR,IAUSR,IBUSR,
     .              TEWL,TIWL,DIWL,VXWL,VYWL,VZWL,WEISPZ)

      RETURN
      END


      SUBROUTINE SIGUSR(IFIRST,JJJ,ZDS,DUMMY1,PSIG,DUMMY2,ARGST,
     .                  XD0,YD0,ZD0,XD1,YD1,ZD1)
C
C  INPUT:
C          IFIRST: FLAG FOR INITIALISATION
C          NCELL:   INDEX IN TALLY ARRAYS FOR CURRENT ZONE
C          JJJ:    INDEX OF SEGMENT ALONG CHORD
C          ZDS:    LENGTH OF SEGMENT NO. JJJ
C  OUTPUT: CONTRIB. FROM CELL NCELL AND CHORD SEGMENT JJJ TO:
C          THE H ALPHA FLUX PSIG(I),I=0,4 CONTRIBUTIONS
C          FROM ATOMS, MOLECULES, TEST IONS AND BULK IONS
C          THE INTEGRANT ARGST SUCH THAT INTEGR.(ARGST*DL) = PSIG
C
      USE PRECISION
      USE PARMMOD
      USE CESTIM
      USE COMUSR
      USE COMPRT
      USE CSPEI
      USE COMSOU
      USE CLOGAU
      USE COMXS
      USE COMSIG
      USE CUPD
      USE CZT1
      USE CTRCEI
      USE CGRID
      USE CCONA
      USE CADGEO
      USE CLGIN
      USE CSDVI
      USE COUTAU
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: JJJ
      INTEGER, INTENT(INOUT) :: IFIRST
      REAL(DP), INTENT(INOUT) :: PSIG(0:NSPZ),ARGST(0:NSPZ,NRAD)
      REAL(DP), INTENT(IN) :: ZDS,DUMMY1,DUMMY2,XD0,YD0,ZD0,XD1,YD1,ZD1
      RETURN
      END
c
c
      subroutine talusr (ICOUNT,VECTOR,TALTOT,TALAV,
     .              TXTTL,TXTSP,TXTUN,ILAST,*)
      USE PRECISION
      USE PARMMOD
      USE CGRID
      USE CGEOM
      implicit NONE
      integer, intent(in) :: icount
      integer, intent(out) :: ilast
      real(dp), intent(in) :: vector(*), TALTOT, TALAV
      character(len=*) :: txttl,txtsp,txtun
      integer :: i

      return 1
      end


      SUBROUTINE TIMUSR(N,X,Y,Z,VX,VY,VZ,N1,N2,T,IC,IE,NP,NL)
      USE PRECISION
      USE PARMMOD
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: X,Y,Z,VX,VY,VZ,T,cx,cy,cz,sc
      INTEGER, INTENT(IN) :: N, N1, N2, IC, IE, NP, IS, NRCELL
      LOGICAL :: NL

      ENTRY NORUSR(is,x,y,z,cx,cy,cz,sc,VX,VY,VZ,NRCELL)

      RETURN
      END


      SUBROUTINE TMSUSR (T0)
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: T0
      RETURN
      END
C
C
      SUBROUTINE UPCUSR(WS,IND)
C
C  USER SUPPLIED COLLISION ESTIMATOR, VOLUME AVERAGED
C
      USE PRECISION
      USE PARMMOD
      USE CESTIM
      USE COMUSR
      USE COMPRT
      USE COMXS
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: WS
      INTEGER, INTENT(IN) :: IND
C
C     WS=WEIGHT/SIGTOT=WEIGHT/(VEL*ZMFPI)=WEIGHT/(VEL*SIGMA,MACR.)
C
C  FOR PARTICLE DENSITY IN CELL NO. NCELL
C     COLV(1,NCELL)=COLV(1,NCELL)+WS
C
      RETURN
      END
c
c
      subroutine upnusr
      return
      end
C
C
      SUBROUTINE UPSUSR(WT,IND)
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: WT
      INTEGER, INTENT(IN) :: IND
      RETURN
      END
C
C
      SUBROUTINE UPTUSR(XSTOR2,XSTORV2,WV)
C
C  USER SUPPLIED TRACKLENGTH ESTIMATOR, VOLUME AVERAGED
C
      USE PRECISION
      USE PARMMOD
      USE CESTIM
      USE COMUSR
      USE COMPRT
      USE CUPD
      USE COMXS
      USE CSPEZ
      USE CGRID
      USE CLOGAU
      USE CCONA
      USE CPOLYG
      USE CZT1
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: XSTOR2(MSTOR1,MSTOR2,N2ND+N3RD),
     .                         XSTORV2(NSTORV,N2ND+N3RD), WV
      REAL(DP) :: CNDYNA(NATM),CNDYNP(NPLS)
CDR
      REAL(DP) :: VPX(NRAD),VPY(NRAD),VRX(NRAD),VRY(NRAD)
CDR
      INTEGER :: IFIRST, IAT, IPL, I, IR, IP, IRD, IA1, IA2, IA3, NA4,
     .           INDEXM, INDEXF
      DATA IFIRST/0/
      SAVE

      IF (IFIRST.EQ.0) THEN
        IFIRST=1
        DO IAT=1,NATMI
          CNDYNA(IAT)=1.D3*AMUA*RMASSA(IAT)
        END DO
        DO IPL=1,NPLSI
          CNDYNP(IPL)=1.D3*AMUA*RMASSP(IPL)
        END DO
C
CDR
CDR  PROVIDE A RADIAL UNIT VECTOR PER CELL
CDR  VPX,VPY,  NEEDED FOR PROJECTING PARTICLE VELOCITIES
CDR  SAME FOR POLOIDAL UNIT VECTOR VRX,VRY
C
        DO I=1,NRAD
          VPX(I)=0.
          VPY(I)=0.
          VRX(I)=0.
          VRY(I)=0.
        END DO
        DO IR=1,NR1STM
          DO IP=1,NP2NDM
            IRD=IR+(IP-1)*NR1P2
            VPX(IRD)=PLNX(IR,IP)
            VPY(IRD)=PLNY(IR,IP)
            VRX(IRD)=PPLNX(IR,IP)
            VRY(IRD)=PPLNY(IR,IP)
          END DO
        END DO
        IA1=NATMI+NMOLI
        IA2=2*IA1
        IA3=3*IA1
        NA4=4*IA1
        INDEXM=NPLSI
        INDEXF=2*NPLSI
      ENDIF
      RETURN
      END





      FUNCTION VDION (I)
      USE PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I
      REAL(DP) :: VDION
      VDION=0.
      RETURN
      END


      SUBROUTINE VECUSR (I,VX,VY,VZ,IPLS)
      USE PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I, IPLS
      REAL(DP), INTENT(IN) :: VX,VY,VZ
      RETURN
      END


      SUBROUTINE VOLUSR(N,A)
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: A(*)
      INTEGER, INTENT(IN) :: N
      RETURN
      END
c
c
      subroutine w7xint
      return
      end
