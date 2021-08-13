      SUBROUTINE MKCENS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'COMNNL'
      INCLUDE 'COMSOU'
      INCLUDE 'COMUSR'
      INCLUDE 'CTEXT'
      INCLUDE 'CINIT'
      SAVE
      DATA IFIRST/0/
C
C  SET DEFAULTS FOR SOURCE DUE TO INITIAL CONDITION, VALID ONLY FOR
C  FIRST TIMESTEP. MODIFIED FOR LATER TIMESTEPS IN SUBR. TMSTEP
C
C  DEFINE ONE MORE STRATUM
      IF (IFIRST.EQ.0) THEN
        IFIRST=1
        NSTRAI=NSTRAI+1
C  CHECK STORAGE
        IF (NSTRAI.GT.NSTRA) THEN
          CALL MASPRM('NSTRA',5,NSTRA,'NSTRAI',6,NSTRAI,IERROR)
          CALL EXIT
        ENDIF
C
        TXTSOU(NSTRAI)='SOURCE DUE TO INITIAL CONDITION          '
C  SOURCE DISTRIBUTION SAMPLED FROM CENSUS ARRAYS RPARTC,IPARTC
        NLCNS(NSTRAI)=.TRUE.
        NLPNT(NSTRAI)=.FALSE.
        NLLNE(NSTRAI)=.FALSE.
        NLSRF(NSTRAI)=.FALSE.
        NLVOL(NSTRAI)=.FALSE.
C  DO NOT CALL IF2COP(NSTRAI)
        INDSRC(NSTRAI)=-1
C
        NLAVRP(NSTRAI)=.FALSE.
        NLAVRT(NSTRAI)=.FALSE.
        NLSYMP(NSTRAI)=.FALSE.
        NLSYMT(NSTRAI)=.FALSE.
        NPTS(NSTRAI)=0
        NINITL(NSTRAI)=2000*NINITL(NSTRAI-1)+1
        NEMODS(NSTRAI)=1
        NAMODS(NSTRAI)=1
        FLUX(NSTRAI)=0.
        NLATM(NSTRAI)=.FALSE.
        NLMOL(NSTRAI)=.FALSE.
        NLION(NSTRAI)=.FALSE.
        NLPLS(NSTRAI)=.FALSE.
        NSPEZ(NSTRAI)=0
        NSRFSI(NSTRAI)=0
C
        SORENI(NSTRAI)=0.
        SORENE(NSTRAI)=0.
        SORVDX(NSTRAI)=0.
        SORVDY(NSTRAI)=0.
        SORVDZ(NSTRAI)=0.
        SORCOS(NSTRAI)=0.
        SORMAX(NSTRAI)=0.
        SORCTX(NSTRAI)=0.
        SORCTY(NSTRAI)=0.
        SORCTZ(NSTRAI)=0.
C
      ENDIF
C
C  READ INITIAL POPULATION FROM PREVIOUS RUN, OVERWRITE DEFAULTS
C
      IF (NFILEJ.EQ.2.OR.NFILEJ.EQ.3) THEN

C  NEW TIMESTEP
        IF (DTIMVN.LE.0.D0) THEN
          DTIMVN=DTIMV
C       ELSE
C         DTIMVN=DTIMVN
        ENDIF
C
C  READ CENSUS ARRAY FROM OLD TIMESTEP
        CALL RSNAP
        DTIMVO=DTIMV
C
        WRITE (6,*) 'INITIAL POPULATION FOR FIRST TIMESTEP'
        WRITE (6,*) 'READ FROM FILE FT 15 '
        WRITE (6,*) 'PARTICLES AND FLUX STORED FOR INITIAL '
        WRITE (6,*) 'DISTRIBUTION IN PREVIOUS RUN '
        CALL MASJ1('IPRNL   ',IPRNL)
        CALL MASR1('FLUX    ',FLUX(NSTRAI))
C
        IF (DTIMVN.NE.DTIMVO) THEN
          FLUX(NSTRAI)=FLUX(NSTRAI)*DTIMVO/DTIMVN
C
          WRITE (6,*) 'FLUX IS RESCALED BY DTIMV_OLD/DTIMV_NEW '
          CALL MASR1('FLUX    ',FLUX(NSTRAI))
          CALL LEER(1)
        ENDIF
C
        DTIMV=DTIMVN
C
        IF (NPTST.LE.0) THEN
          NPTS(NSTRAI)=IPRNL
        ELSEIF (NPTST.GT.0) THEN
          NPTS(NSTRAI)=NPTST
        ENDIF
C
        IF (NPTS(NSTRAI).GT.0.AND.FLUX(NSTRAI).GT.0) THEN
          NSRFSI(NSTRAI)=1
          SORWGT(1,NSTRAI)=1.D0
        ENDIF
C
        CALL LEER(2)
        DO I=1,IPRNL
          RPARTC(I,10)=TIME0
        ENDDO
        WRITE (6,*) 'PARTICLE CLOCK RESET TO TIME0'
        WRITE (6,*) 'FIRST TIMESTEP RUNS FROM TIM1 TO TIM2:  '
          CALL MASR2('TIM1, TIM2      ',TIME0,TIME0+DTIMV)
          CALL LEER(2)
C
      ENDIF
      RETURN
      END
C
      SUBROUTINE MULTI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'CGRID'
      INCLUDE 'COMUSR'
      INCLUDE 'CGEOM'
      INCLUDE 'CINIT'
C
C  GEOMETRY DATA
C
      ENTRY MULTIG
C
C  ZONE VOLUMES, KNOWN IN ZONE 1 TO NSTRD
      DO 130 J=2,NBMLT
        DO 120 I=1,NSTRD
          VOL(I+(J-1)*NSTRD)=VOL(I)*VOLCOR(J)
c slmod begin - tr
c...      VOL2 is the toroidal volume and is used for the
c         calculation of the volume recombination source:
          VOL2(I+(J-1)*NSTRD)=VOL2(I)*VOLCOR(J)
c slmod end
120     CONTINUE
130   CONTINUE
      DO 140 I=1,NSTRD
        VOL(I)=VOL(I)*VOLCOR(1)
c slmod begin - tr
        VOL2(I)=VOL2(I)*VOLCOR(1)
c slmod end
140   CONTINUE
C
      RETURN
C
C  PLASMA DATA
C
      ENTRY MULTIP
C
C  INDPRO.LE.4: ONLY RADIAL PLASMA PROFILES GIVEN
C  RADIAL PLASMA PROFILES, KNOWN IN ZONES 1 TO NR1ST
      DO 210 J=2,NBLCKS
        IF (INDPRO(1).LE.4) THEN
          DO 201 I=1,NR1ST
            TEIN(I+(J-1)*NR1ST)=TEIN(I)
201       CONTINUE
        ENDIF
        IF (INDPRO(2).LE.4) THEN
          DO 202 K=1,NPLSI
            DO 202 I=1,NR1ST
              TIIN(K,I+(J-1)*NR1ST)=TIIN(K,I)
202       CONTINUE
        ENDIF
        IF (INDPRO(3).LE.4) THEN
          DO 204 K=1,NPLSI
            DO 204 I=1,NR1ST
              DIIN(K,I+(J-1)*NR1ST)=DIIN(K,I)
204       CONTINUE
        ENDIF
        IF (INDPRO(4).LE.4) THEN
          DO 205 K=1,NPLSI
            DO 205 I=1,NR1ST
              VXIN(K,I+(J-1)*NR1ST)=VXIN(K,I)
              VYIN(K,I+(J-1)*NR1ST)=VYIN(K,I)
              VZIN(K,I+(J-1)*NR1ST)=VZIN(K,I)
205       CONTINUE
        ENDIF
        IF (INDPRO(5).LE.4) THEN
          DO 206 I=1,NR1ST
            BXIN(I+(J-1)*NR1ST)=BXIN(I)
            BYIN(I+(J-1)*NR1ST)=BYIN(I)
            BZIN(I+(J-1)*NR1ST)=BZIN(I)
206       CONTINUE
        ENDIF
        IF (INDPRO(6).LE.4) THEN
          DO 207 K=1,NAINI
          DO 207 I=1,NR1ST
            ADIN(K,I+(J-1)*NR1ST)=ADIN(K,I)
207       CONTINUE
        ENDIF
210   CONTINUE
C
C  INDPRO.GT.4: ONLY NSTRD=NR1ST*NP2ND*NT3RD PLASMA DATA GIVEN
      DO 310 J=2,NBMLT
        IF (INDPRO(1).GT.4) THEN
          DO 301 I=1,NSTRD
            TEIN(I+(J-1)*NSTRD)=TEIN(I)
301       CONTINUE
        ENDIF
        IF (INDPRO(2).GT.4) THEN
          DO 302 K=1,NPLSI
            DO 302 I=1,NSTRD
              TIIN(K,I+(J-1)*NSTRD)=TIIN(K,I)
302       CONTINUE
        ENDIF
        IF (INDPRO(3).GT.4) THEN
          DO 304 K=1,NPLSI
            DO 304 I=1,NSTRD
              DIIN(K,I+(J-1)*NSTRD)=DIIN(K,I)
304       CONTINUE
        ENDIF
        IF (INDPRO(4).GT.4) THEN
          DO 305 K=1,NPLSI
            DO 305 I=1,NSTRD
              VXIN(K,I+(J-1)*NSTRD)=VXIN(K,I)
              VYIN(K,I+(J-1)*NSTRD)=VYIN(K,I)
              VZIN(K,I+(J-1)*NSTRD)=VZIN(K,I)
305       CONTINUE
        ENDIF
        IF (INDPRO(5).GT.4) THEN
          DO 306 I=1,NSTRD
            BXIN(I+(J-1)*NSTRD)=BXIN(I)
            BYIN(I+(J-1)*NSTRD)=BYIN(I)
            BZIN(I+(J-1)*NSTRD)=BZIN(I)
306       CONTINUE
        ENDIF
        IF (INDPRO(6).LE.4) THEN
          DO 307 K=1,NAINI
          DO 307 I=1,NSTRD
            ADIN(K,I+(J-1)*NSTRD)=ADIN(K,I)
307       CONTINUE
        ENDIF
310   CONTINUE
C
c slmod begin - tr
c...  Duplicate local multipliers passed from DIVIMP:
      DO IRRC=1,NREC
        DO J=2,NBMLT
          DO I=1,NSTRD
            TABRCM(IRRC,I+(J-1)*NSTRD)=TABRCM(IRRC,I)
            TABDSM(IRRC,I+(J-1)*NSTRD)=TABDSM(IRRC,I)
            TABREM(IRRC,I+(J-1)*NSTRD)=TABREM(IRRC,I)
            TABDEM(IRRC,I+(J-1)*NSTRD)=TABDEM(IRRC,I)
          ENDDO
        ENDDO
      ENDDO
c slmod end
      RETURN
      END
C
C
C FROM DETLEV, FEB 23, 2000
      SUBROUTINE SETEQ_NEW
c slmod begin - not tr
c      USE COMUSR
c slmod end
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
c slmod begin - not tr
      INCLUDE 'COMUSR'
c slmod end
      INCLUDE 'CADGEO'
      INCLUDE 'CTRCEI'
      INCLUDE 'CLGIN'
      LOGICAL BITGET
C
C  FIRST: ADDITIONAL SURFACES
C
      DO 1 J=1,NLIMI
        IEQ=IABS(ILEQUI(J))
        IF (IEQ.EQ.0) GOTO 1
!PB     IF (LGJUM0(J).OR.LGJUM0(IEQ)) GOTO 1
        IF (IGJUM0(J)+IGJUM0(IEQ) .NE. 0) GOTO 1
        SG=ISIGN(1,ILEQUI(J))
C   COEFFICIENTS OF SURFACE J ARE SET EQUAL TO THOSE OF SURFACE IEQ
        A0LM(J)=SG*A0LM(IEQ)
        A1LM(J)=SG*A1LM(IEQ)
        A2LM(J)=SG*A2LM(IEQ)
        A3LM(J)=SG*A3LM(IEQ)
        A4LM(J)=SG*A4LM(IEQ)
        A5LM(J)=SG*A5LM(IEQ)
        A6LM(J)=SG*A6LM(IEQ)
        A7LM(J)=SG*A7LM(IEQ)
        A8LM(J)=SG*A8LM(IEQ)
        A9LM(J)=SG*A9LM(IEQ)
        IF (JUMLIM(J).EQ.0) THEN
          ALM(J)=SG*ALM(IEQ)
          BLM(J)=SG*BLM(IEQ)
          CLM(J)=SG*CLM(IEQ)
        ELSE
          ALM(J)=ALM(IEQ)
          BLM(J)=BLM(IEQ)
          CLM(J)=CLM(IEQ)
        ENDIF
        IF (JUMLIM(J).NE.0) THEN
          IF (NLIMPB >= NLIMPS) THEN
            IGJUM1(J,IEQ)=1
            IGJUM1(IEQ,J)=1
          ELSE
            CALL BITSET (IGJUM1,0,NLIMPS,J,IEQ,1,NBITS)
            CALL BITSET (IGJUM1,0,NLIMPS,IEQ,J,1,NBITS)
          END IF
        ELSE
          IF (NLIMPB >= NLIMPS) THEN
            IGJUM2(J,IEQ)=1
            IGJUM2(IEQ,J)=1
          ELSE
            CALL BITSET (IGJUM2,0,NLIMPS,J,IEQ,1,NBITS)
            CALL BITSET (IGJUM2,0,NLIMPS,IEQ,J,1,NBITS)
          END IF
        ENDIF
C
        IF (TRCSUR) THEN
          WRITE (6,*) 'EQUATION OF SURFACE NO. ',J,' IS SET EQUAL TO '
          WRITE (6,*) 'EQUATION OF SURFACE NO. ',IEQ
          CALL LEER(1)
        ENDIF
1     CONTINUE
C
C  NEXT: NON DEFAULT STANDARD SURFACES
C  JUMLIM=0 FOR STANDARD SURFACES, ONLY LGJUM2 IS USED IN TIME-ROUTINES
C
      DO 10 J=NLIM+1,NLIM+NSTSI
        IF (ILEQUI(J).EQ.0.OR.IGJUM0(J).NE.0) GOTO 10
        IF (INUMP(J-NLIM,1).NE.0) IDIMP=1
        IF (INUMP(J-NLIM,2).NE.0) IDIMP=2
        IF (INUMP(J-NLIM,3).NE.0) IDIMP=3
C  NEGATIVE SIGN OF ILEQUI IS MEANINGLESS FOR STANDARD SURFACES
C  BECAUSE THEIR ORIENTATION CAN NOT BE CHANGED
        IEQ1=IABS(ILEQUI(J))
C  FIND INDEX IEQ OF SURFACE IEQ1, TO BE EQUALLED WITH SURFACE J
C  SEARCH ONLY IN THE SAME (RADIAL, POLOIDAL OR TOROIDAL) GRID
        IEQ=0
        DO 11 I=1,NSTSI
          IF (INUMP(I,IDIMP).EQ.IEQ1) IEQ=NLIM+I
11      CONTINUE
        IF (NLIMPB >= NLIMPS) THEN
          IGJUM2(J,IEQ)=1
          IGJUM2(IEQ,J)=1
        ELSE
          CALL BITSET (IGJUM2,0,NLIMPS,J,IEQ,1,NBITS)
          CALL BITSET (IGJUM2,0,NLIMPS,IEQ,J,1,NBITS)
        END IF
        IF (TRCSUR) THEN
          WRITE (6,*) 'EQUATION OF SURFACE NO. ',J,' IS SET EQUAL TO '
          WRITE (6,*) 'EQUATION OF SURFACE NO. ',IEQ
          CALL LEER(1)
        ENDIF
10    CONTINUE
C
C  LGJUM1: COMMUTATIVE
C  LGJUM2: ASSOCIATIVE
      IF (NLIMPB >= NLIMPS) THEN
        DO J=1,NLIMI
          IF (SUM(IGJUM1(J,1:NLIMI)) > 1 ) THEN
            DO I=1,NLIMI
              IF (IGJUM1(J,I) .NE. 0) IGJUM1(I,J)=1
            ENDDO
          ENDIF
          IF (JUMLIM(J).EQ.0) THEN
            IF (SUM(IGJUM2(J,1:NLIMI)) > 1 ) THEN
              DO I=1,NLIMI
                DO K=1,NLIMI
                  IF (IGJUM2(I,J)*IGJUM2(J,K) .NE. 0) THEN
                    IGJUM2(I,K)=1
                    IGJUM2(K,I)=1
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ELSE
        NLBT = NLIMI/NBITS+1
        DO J=1,NLIMI
          IF (SUM(IGJUM1(J,1:NLBT)) > 1 ) THEN
            DO I=1,NLIMI
              IF (BITGET(IGJUM1,0,NLIMPS,J,I,NBITS))
     .          CALL BITSET (IGJUM1,0,NLIMPS,I,J,1,NBITS)
            END DO
          END IF
          IF (JUMLIM(J).EQ.0) THEN
            IF (SUM(IGJUM2(I,1:NLBT)) > 1 ) THEN
              DO I=1,NLIMI
                IF (BITGET(IGJUM2,0,NLIMPS,I,J,NBITS)) THEN
                  DO K=1,NLIMI
                    IF (BITGET(IGJUM1,0,NLIMPS,J,K,NBITS)) THEN
                      CALL BITSET (IGJUM2,0,NLIMPS,I,K,1,NBITS)
                      CALL BITSET (IGJUM2,0,NLIMPS,K,I,1,NBITS)
                    END IF
                  END DO
                END IF
              END DO
            END IF
          END IF
        END DO
      END IF
      DO 200 J=NLIM+1,NLIM+NSTSI
        DO 200 I=NLIM+1,NLIM+NSTSI
        IF (NLIMPB >= NLIMPS) THEN
          IF (IGJUM1(J,I).NE.0) IGJUM1(I,J)=1
          DO K=NLIM+1,NLIM+NSTSI
            IF (IGJUM2(I,J)*IGJUM2(J,K) .NE. 0) THEN
              IGJUM2(I,K)=1
              IGJUM2(K,I)=1
            ENDIF
          END DO
        ELSE
          IF (BITGET(IGJUM1,0,NLIMPS,J,I,NBITS))
     .      CALL BITSET (IGJUM1,0,NLIMPS,I,J,1,NBITS)
          DO K=NLIM+1,NLIM+NSTSI
            IF (BITGET(IGJUM2,0,NLIMPS,I,J,NBITS) .AND.
     .          BITGET(IGJUM1,0,NLIMPS,J,K,NBITS)) THEN
              CALL BITSET (IGJUM2,0,NLIMPS,I,K,1,NBITS)
              CALL BITSET (IGJUM2,0,NLIMPS,K,I,1,NBITS)
            END IF
          END DO
        END IF
200   CONTINUE
      RETURN
      END
C
      SUBROUTINE SETEQ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'CADGEO'
      INCLUDE 'COMUSR'
      INCLUDE 'CTRCEI'
      INCLUDE 'CLGIN'
      LOGICAL BITGET
C
C  FIRST: ADDITIONAL SURFACES
C
      IF (output) WRITE(0,*) 'MARK: 1ST LOOP'
      DO 1 J=1,NLIMI
        IEQ=IABS(ILEQUI(J))
        IF (IEQ.EQ.0) GOTO 1
!PB     IF (LGJUM0(J).OR.LGJUM0(IEQ)) GOTO 1
        IF (IGJUM0(J)+IGJUM0(IEQ) .NE. 0) GOTO 1
        SG=ISIGN(1,ILEQUI(J))
C   COEFFICIENTS OF SURFACE J ARE SET EQUAL TO THOSE OF SURFACE IEQ
        A0LM(J)=SG*A0LM(IEQ)
        A1LM(J)=SG*A1LM(IEQ)
        A2LM(J)=SG*A2LM(IEQ)
        A3LM(J)=SG*A3LM(IEQ)
        A4LM(J)=SG*A4LM(IEQ)
        A5LM(J)=SG*A5LM(IEQ)
        A6LM(J)=SG*A6LM(IEQ)
        A7LM(J)=SG*A7LM(IEQ)
        A8LM(J)=SG*A8LM(IEQ)
        A9LM(J)=SG*A9LM(IEQ)
        IF (JUMLIM(J).EQ.0) THEN
          ALM(J)=SG*ALM(IEQ)
          BLM(J)=SG*BLM(IEQ)
          CLM(J)=SG*CLM(IEQ)
        ELSE
          ALM(J)=ALM(IEQ)
          BLM(J)=BLM(IEQ)
          CLM(J)=CLM(IEQ)
        ENDIF
        IF (JUMLIM(J).NE.0) THEN
          IF (NLIMPB >= NLIMPS) THEN
            IGJUM1(J,IEQ)=1
            IGJUM1(IEQ,J)=1
          ELSE
            CALL BITSET (IGJUM1,0,NLIMPS,J,IEQ,1,NBITS)
            CALL BITSET (IGJUM1,0,NLIMPS,IEQ,J,1,NBITS)
          END IF
        ELSE
          IF (NLIMPB >= NLIMPS) THEN
            IGJUM2(J,IEQ)=1
            IGJUM2(IEQ,J)=1
          ELSE
            CALL BITSET (IGJUM2,0,NLIMPS,J,IEQ,1,NBITS)
            CALL BITSET (IGJUM2,0,NLIMPS,IEQ,J,1,NBITS)
          END IF
        ENDIF
C
        IF (TRCSUR) THEN
          WRITE (6,*) 'EQUATION OF SURFACE NO. ',J,' IS SET EQUAL TO '
          WRITE (6,*) 'EQUATION OF SURFACE NO. ',IEQ
          CALL LEER(1)
        ENDIF
1     CONTINUE
C
C  NEXT: NON DEFAULT STANDARD SURFACES
C  JUMLIM=0 FOR STANDARD SURFACES, ONLY LGJUM2 IS USED IN TIME-ROUTINES
C
      IF (output) WRITE(0,*) 'MARK: SETEQ 10 LOOP'
      DO 10 J=NLIM+1,NLIM+NSTSI
        IF (ILEQUI(J).EQ.0.OR.IGJUM0(J).NE.0) GOTO 10
        IF (INUMP(J-NLIM,1).NE.0) IDIMP=1
        IF (INUMP(J-NLIM,2).NE.0) IDIMP=2
        IF (INUMP(J-NLIM,3).NE.0) IDIMP=3
C  NEGATIVE SIGN OF ILEQUI IS MEANINGLESS FOR STANDARD SURFACES
C  BECAUSE THEIR ORIENTATION CAN NOT BE CHANGED
        IEQ1=IABS(ILEQUI(J))
C  FIND INDEX IEQ OF SURFACE IEQ1, TO BE EQUALLED WITH SURFACE J
C  SEARCH ONLY IN THE SAME (RADIAL, POLOIDAL OR TOROIDAL) GRID
        IEQ=0
        DO 11 I=1,NSTSI
          IF (INUMP(I,IDIMP).EQ.IEQ1) IEQ=NLIM+I
11      CONTINUE
        IF (NLIMPB >= NLIMPS) THEN
          IGJUM2(J,IEQ)=1
          IGJUM2(IEQ,J)=1
        ELSE
          CALL BITSET (IGJUM2,0,NLIMPS,J,IEQ,1,NBITS)
          CALL BITSET (IGJUM2,0,NLIMPS,IEQ,J,1,NBITS)
        END IF
        IF (TRCSUR) THEN
          WRITE (6,*) 'EQUATION OF SURFACE NO. ',J,' IS SET EQUAL TO '
          WRITE (6,*) 'EQUATION OF SURFACE NO. ',IEQ
          CALL LEER(1)
        ENDIF
10    CONTINUE
C
C  LGJUM1: COMMUTATIVE
C  LGJUM2: ASSOCIATIVE
      IF (output) WRITE(0,*) 'MARK: SETEQ 100 LOOP'
      DO 100 J=1,NLIMI
        DO 100 I=1,NLIMI
        IF (NLIMPB >= NLIMPS) THEN
          IF (IGJUM1(J,I).NE.0) IGJUM1(I,J)=1
          IF (JUMLIM(J).EQ.0) THEN
            DO K=1,NLIMI
              IF (IGJUM2(I,J)*IGJUM2(J,K) .NE. 0) THEN
                IGJUM2(I,K)=1
                IGJUM2(K,I)=1
              ENDIF
            END DO
          ENDIF
        ELSE
          IF (BITGET(IGJUM1,0,NLIMPS,J,I,NBITS))
     .      CALL BITSET (IGJUM1,0,NLIMPS,I,J,1,NBITS)
c slmod begin - not tr (substantial changes to EIRENE02 code here)
          IF (JUMLIM(J).EQ.0) THEN
            DO K=1,NLIMI
              IF (BITGET(IGJUM2,0,NLIMPS,I,J,NBITS) .AND.
     .            BITGET(IGJUM1,0,NLIMPS,J,K,NBITS)) THEN
                CALL BITSET (IGJUM2,0,NLIMPS,I,K,1,NBITS)
                CALL BITSET (IGJUM2,0,NLIMPS,K,I,1,NBITS)
              END IF
            END DO
          ENDIF
c
c          DO K=1,NLIMI
c            IF (BITGET(IGJUM2,0,NLIMPS,I,J,NBITS) .AND.
c     .          BITGET(IGJUM1,0,NLIMPS,J,K,NBITS)) THEN
c              CALL BITSET (IGJUM2,0,NLIMPS,I,K,1,NBITS)
c              CALL BITSET (IGJUM2,0,NLIMPS,K,I,1,NBITS)
c            END IF
c          END DO
c slmod end
        END IF
100   CONTINUE
      IF (output) WRITE(0,*) 'MARK: SETEQ 200 LOOP'
      DO 200 J=NLIM+1,NLIM+NSTSI
        DO 200 I=NLIM+1,NLIM+NSTSI
        IF (NLIMPB >= NLIMPS) THEN
          IF (IGJUM1(J,I).NE.0) IGJUM1(I,J)=1
          DO K=NLIM+1,NLIM+NSTSI
            IF (IGJUM2(I,J)*IGJUM2(J,K) .NE. 0) THEN
              IGJUM2(I,K)=1
              IGJUM2(K,I)=1
            ENDIF
          END DO
        ELSE
          IF (BITGET(IGJUM1,0,NLIMPS,J,I,NBITS))
     .      CALL BITSET (IGJUM1,0,NLIMPS,I,J,1,NBITS)
          DO K=NLIM+1,NLIM+NSTSI
            IF (BITGET(IGJUM2,0,NLIMPS,I,J,NBITS) .AND.
     .          BITGET(IGJUM1,0,NLIMPS,J,K,NBITS)) THEN
              CALL BITSET (IGJUM2,0,NLIMPS,I,K,1,NBITS)
              CALL BITSET (IGJUM2,0,NLIMPS,K,I,1,NBITS)
            END IF
          END DO
        END IF
200   CONTINUE
      RETURN
      END
C
      SUBROUTINE SETFIT(TRCSUR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'CADGEO'
      INCLUDE 'CLGIN'
      INCLUDE 'CCONA'
C
      DIMENSION IEQ(2),XS(8),YS(8),DST(8,2)
      DIMENSION XA(14),XB(14)
      LOGICAL TRCSUR
      LOGICAL LINFX,LINFY,LINFZ
C
      DO 1 I=1,NLIMI
C  IS ILFIT OPTION IN USE?
        IF (ILFIT(I).EQ.0.OR.IGJUM0(I).NE.0) GOTO 1
C  SELECT THE SURFACE NUMBERS OF THE NEIGHBORING SURFACES
        IEQ(1)=ILFIT(I)/1000
        IEQ(2)=ILFIT(I)-IEQ(1)*1000
C  SURFACE I MUST BE GIVEN BY TWO POINT OPTION WITH ONE
C  IGNORABLE CO-ORDINATE
C  USE THIRD POINT FOR IDENTIFICATION OF 2-POINT INPUT OPTION
C  BECAUSE RLB HAS ALREADY BEEN OVERWRITTEN
        IF (P3(1,I).LT.1.D50.AND.P3(2,I).LT.1.D50.AND.
     .      P3(3,I).LT.1.D50) GOTO 991
C  WHICH CO-ORDINATE OF SURFACE NO. I IS IGNORABLE?
        LINFX=P3(1,I).GT.1.D50
        LINFY=P3(2,I).GT.1.D50
        LINFZ=P3(3,I).GT.1.D50
        IPNT1=0
        IPNT2=0
        DO 2 J=1,2
          IF (IEQ(J).EQ.0) GOTO 2
C  CONNECT SURFACE I WITH SURFACE NO. IE
          IE=IEQ(J)
C  SURFACE IE MUST BE GIVEN BY RLB=1 OR RLB=1.5 OPTION, WITH SAME
C  IGNORABLE CO-ORDINATE AS SURFACE I
          IF ((RLB(IE).NE.1..AND.RLB(IE).NE.1.5).OR.
     .        (IGJUM0(IE).NE.0)) GOTO 991
C
C  IDENTIFY IGNORABLE CO-ORDINATE OF SURFACE IE
          IF (LINFX) THEN
            IF (A1LM(IE).NE.0..OR.A4LM(IE).NE.0..OR.
     .          A7LM(IE).NE.0..OR.A8LM(IE).NE.0.D0) GOTO 5
C  X IS IGNORABLE, IN BOTH SURFACES: I AND IE
            A0=A0LM(IE)
            A1=A2LM(IE)
            A2=A3LM(IE)
            A3=A5LM(IE)
            A4=A6LM(IE)
            A5=A9LM(IE)
            XL=YLIMS1(IE,1)
            XR=YLIMS2(IE,1)
            YL=ZLIMS1(IE,1)
            YR=ZLIMS2(IE,1)
            IF (JUMLIM(I).NE.0) THEN
              XP1=P1(2,I)
              XP2=P2(2,I)
              YP1=P1(3,I)
              YP2=P2(3,I)
            ELSE
              B0=A0LM(I)
              B1=A2LM(I)
              B2=A3LM(I)
              B3=A5LM(I)
              B4=A6LM(I)
              B5=A9LM(I)
              XMIN=YLIMS1(I,1)
              XMAX=YLIMS2(I,1)
              YMIN=ZLIMS1(I,1)
              YMAX=ZLIMS2(I,1)
            ENDIF
            GOTO 100
          ENDIF
5         CONTINUE
          LINFX=.FALSE.
          IF (LINFY) THEN
            IF (A2LM(IE).NE.0..OR.A5LM(IE).NE.0..OR.
     .          A7LM(IE).NE.0..OR.A9LM(IE).NE.0.D0) GOTO 6
C  Y IS IGNORABLE, IN BOTH SURFACES: I AND IE
            A0=A0LM(IE)
            A1=A1LM(IE)
            A2=A3LM(IE)
            A3=A4LM(IE)
            A4=A6LM(IE)
            A5=A8LM(IE)
            XL=XLIMS1(IE,1)
            XR=XLIMS2(IE,1)
            YL=ZLIMS1(IE,1)
            YR=ZLIMS2(IE,1)
            IF (JUMLIM(I).NE.0) THEN
              XP1=P1(1,I)
              XP2=P2(1,I)
              YP1=P1(3,I)
              YP2=P2(3,I)
            ELSE
              B0=A0LM(I)
              B1=A1LM(I)
              B2=A3LM(I)
              B3=A4LM(I)
              B4=A6LM(I)
              B5=A8LM(I)
              XMIN=XLIMS1(I,1)
              XMAX=XLIMS2(I,1)
              YMIN=ZLIMS1(I,1)
              YMAX=ZLIMS2(I,1)
            ENDIF
            GOTO 100
          ENDIF
6         CONTINUE
          LINFY=.FALSE.
          IF (LINFZ) THEN
            IF (A3LM(IE).NE.0..OR.A6LM(IE).NE.0..OR.
     .          A8LM(IE).NE.0..OR.A9LM(IE).NE.0.D0) GOTO 7
C  Z IS IGNORABLE, IN BOTH SURFACES: I AND IE
            A0=A0LM(IE)
            A1=A1LM(IE)
            A2=A2LM(IE)
            A3=A4LM(IE)
            A4=A5LM(IE)
            A5=A7LM(IE)
            XL=XLIMS1(IE,1)
            XR=XLIMS2(IE,1)
            YL=YLIMS1(IE,1)
            YR=YLIMS2(IE,1)
            IF (JUMLIM(I).NE.0) THEN
              XP1=P1(1,I)
              XP2=P2(1,I)
              YP1=P1(2,I)
              YP2=P2(2,I)
            ELSE
              B0=A0LM(I)
              B1=A1LM(I)
              B2=A2LM(I)
              B3=A4LM(I)
              B4=A5LM(I)
              B5=A7LM(I)
              XMIN=XLIMS1(I,1)
              XMAX=XLIMS2(I,1)
              YMIN=YLIMS1(I,1)
              YMAX=YLIMS2(I,1)
            ENDIF
            GOTO 100
          ENDIF
7         CONTINUE
          LINFZ=.FALSE.
          GOTO 990
100       CONTINUE
C
C
          IS=0
C  SCHNITTPUNKTE MIT X=XL AND X=XR
          IF (ABS(A4).LE.EPS12) THEN
            YT=-(A0+A1*XL+A3*XL*XL)/(A2+A5*XL)
            IF (YT.GE.YL.AND.YT.LE.YR) THEN
              IS=IS+1
              XS(IS)=XL
              YS(IS)=YT
            ENDIF
            YT=-(A0+A1*XR+A3*XR*XR)/(A2+A5*XR)
            IF (YT.GE.YL.AND.YT.LE.YR) THEN
              IS=IS+1
              XS(IS)=XR
              YS(IS)=YT
            ENDIF
          ELSE
            AH=0.5*(A2+A5*XL)/A4
            B=(A0+A1*XL+A3*XL*XL)/A4
            RAD=AH*AH-B
            IF (RAD.GE.0.D0) THEN
              SRAD=SQRT(RAD)
              YT=-AH+SRAD
              IF (YT.GE.YL.AND.YT.LE.YR) THEN
                IS=IS+1
                XS(IS)=XL
                YS(IS)=YT
              ENDIF
              YT=-AH-SRAD
              IF (YT.GE.YL.AND.YT.LE.YR) THEN
                IS=IS+1
                XS(IS)=XL
                YS(IS)=YT
              ENDIF
            ENDIF
            AH=0.5*(A2+A5*XR)/A4
            B=(A0+A1*XR+A3*XR*XR)/A4
            RAD=AH*AH-B
            IF (RAD.GE.0.D0) THEN
              SRAD=SQRT(RAD)
              YT=-AH+SRAD
              IF (YT.GE.YL.AND.YT.LE.YR) THEN
                IS=IS+1
                XS(IS)=XR
                YS(IS)=YT
              ENDIF
              YT=-AH-SRAD
              IF (YT.GE.YL.AND.YT.LE.YR) THEN
                IS=IS+1
                XS(IS)=XR
                YS(IS)=YT
              ENDIF
            ENDIF
          ENDIF
C  SCHNITTPUNKTE MIT Y=YL AND Y=YR
          IF (ABS(A3).LE.EPS12) THEN
            XT=-(A0+A2*YL+A4*YL*YL)/(A1+A5*YL)
            IF (XT.GE.XL.AND.XT.LE.XR) THEN
              IS=IS+1
              YS(IS)=YL
              XS(IS)=XT
            ENDIF
            XT=-(A0+A2*YR+A4*YR*YR)/(A1+A5*YR)
            IF (XT.GE.XL.AND.XT.LE.XR) THEN
              IS=IS+1
              YS(IS)=YR
              XS(IS)=XT
            ENDIF
          ELSE
            AH=0.5*(A1+A5*YL)/A3
            B=(A0+A2*YL+A4*YL*YL)/A3
            RAD=AH*AH-B
            IF (RAD.GE.0.D0) THEN
              SRAD=SQRT(RAD)
              XT=-AH+SRAD
              IF (XT.GE.XL.AND.XT.LE.XR) THEN
                IS=IS+1
                YS(IS)=YL
                XS(IS)=XT
              ENDIF
              XT=-AH-SRAD
              IF (XT.GE.XL.AND.XT.LE.XR) THEN
                IS=IS+1
                YS(IS)=YL
                XS(IS)=XT
              ENDIF
            ENDIF
            AH=0.5*(A1+A5*YR)/A3
            B=(A0+A2*YR+A4*YR*YR)/A3
            RAD=AH*AH-B
            IF (RAD.GE.0.D0) THEN
              SRAD=SQRT(RAD)
              XT=-AH+SRAD
              IF (XT.GE.XL.AND.XT.LE.XR) THEN
                IS=IS+1
                YS(IS)=YR
                XS(IS)=XT
              ENDIF
              XT=-AH-SRAD
              IF (XT.GE.XL.AND.XT.LE.XR) THEN
                IS=IS+1
                YS(IS)=YR
                XS(IS)=XT
              ENDIF
            ENDIF
          ENDIF
C
C  SELECT THE DISTANCES TO THE INTERSECTION POINTS
          DO 3 K=1,IS
            DST(K,1)=1.D60
            DST(K,2)=1.D60
            IF (IPNT1.EQ.0)
     .      DST(K,1)=SQRT((XS(K)-XP1)**2+(YS(K)-YP1)**2)
            IF (IPNT2.EQ.0)
     .      DST(K,2)=SQRT((XS(K)-XP2)**2+(YS(K)-YP2)**2)
3         CONTINUE
C
          IMIN1=1
          IMIN2=1
          DO 4 K=2,IS
            IF (DST(K,1).LT.DST(IMIN1,1)) IMIN1=K
            IF (DST(K,2).LT.DST(IMIN2,2)) IMIN2=K
4         CONTINUE
C
          IF (TRCSUR) THEN
            WRITE (6,*) 'SETFIT, I,IE ',I,IE
            DO 4711 K=1,IS
              WRITE (6,*) 'K,XS,YS,DIST1,DIST2 ',
     .                     K,XS(K),YS(K),DST(K,1),DST(K,2)
4711        CONTINUE
          ENDIF
C
C  SET THE SELECTED POINT
          IF (DST(IMIN1,1).LT.DST(IMIN2,2)) THEN
            XP1=XS(IMIN1)
            YP1=YS(IMIN1)
            IPNT1=1
            IF (TRCSUR)
     .      WRITE (6,*) 'PNT.1 REPLACED BY ',XP1,YP1,' FOR SURF. ',I
          ELSE
            XP2=XS(IMIN2)
            YP2=YS(IMIN2)
            IPNT2=1
            IF (TRCSUR)
     .      WRITE (6,*) 'PNT.2 REPLACED BY ',XP2,YP2,' FOR SURF. ',I
          ENDIF
C
C  RESET THE POINTS ON THE ARRAYS
          IF (TRCSUR) WRITE (6,*) ' LINFXYZ ',LINFX,LINFY,LINFZ
          IF (LINFX) THEN
            P1(2,I)=XP1
            P1(3,I)=YP1
            P2(2,I)=XP2
            P2(3,I)=YP2
            DZ3=(P1(3,I)-P2(3,I))
            DY2=(P1(2,I)-P2(2,I))
            A0LM(I)=DZ3*P1(2,I)-DY2*P1(3,I)
            A1LM(I)=0.
            A2LM(I)=-DZ3
            A3LM(I)=DY2
            YLIMS1(I,1)=MIN(P1(2,I),P2(2,I))
            YLIMS2(I,1)=MAX(P1(2,I),P2(2,I))
            ZLIMS1(I,1)=MIN(P1(3,I),P2(3,I))
            ZLIMS2(I,1)=MAX(P1(3,I),P2(3,I))
            IF (P1(2,I).EQ.P2(2,I)) THEN
              YLIMS1(I,1)=YLIMS1(I,1)-0.1
              YLIMS2(I,1)=YLIMS2(I,1)+0.1
            ENDIF
            IF (P1(3,I).EQ.P2(3,I)) THEN
              ZLIMS1(I,1)=ZLIMS1(I,1)-0.1
              ZLIMS2(I,1)=ZLIMS2(I,1)+0.1
            ENDIF
          ELSEIF (LINFY) THEN
            P1(1,I)=XP1
            P1(3,I)=YP1
            P2(1,I)=XP2
            P2(3,I)=YP2
            DZ3=(P1(3,I)-P2(3,I))
            DX1=(P1(1,I)-P2(1,I))
            A0LM(I)=DZ3*P1(1,I)-DX1*P1(3,I)
            A1LM(I)=-DZ3
            A2LM(I)=0.
            A3LM(I)=DX1
            XLIMS1(I,1)=MIN(P1(1,I),P2(1,I))
            XLIMS2(I,1)=MAX(P1(1,I),P2(1,I))
            ZLIMS1(I,1)=MIN(P1(3,I),P2(3,I))
            ZLIMS2(I,1)=MAX(P1(3,I),P2(3,I))
            IF (P1(1,I).EQ.P2(1,I)) THEN
              XLIMS1(I,1)=XLIMS1(I,1)-0.1
              XLIMS2(I,1)=XLIMS2(I,1)+0.1
            ENDIF
            IF (P1(3,I).EQ.P2(3,I)) THEN
              ZLIMS1(I,1)=ZLIMS1(I,1)-0.1
              ZLIMS2(I,1)=ZLIMS2(I,1)+0.1
            ENDIF
          ELSEIF (LINFZ) THEN
            P1(1,I)=XP1
            P1(2,I)=YP1
            P2(1,I)=XP2
            P2(2,I)=YP2
            DY2=(P1(2,I)-P2(2,I))
            DX1=(P1(1,I)-P2(1,I))
            A0LM(I)=DY2*P1(1,I)-DX1*P1(2,I)
            A1LM(I)=-DY2
            A2LM(I)=DX1
            A3LM(I)=0.
            XLIMS1(I,1)=MIN(P1(1,I),P2(1,I))
            XLIMS2(I,1)=MAX(P1(1,I),P2(1,I))
            YLIMS1(I,1)=MIN(P1(2,I),P2(2,I))
            YLIMS2(I,1)=MAX(P1(2,I),P2(2,I))
            IF (P1(1,I).EQ.P2(1,I)) THEN
              XLIMS1(I,1)=XLIMS1(I,1)-0.1
              XLIMS2(I,1)=XLIMS2(I,1)+0.1
            ENDIF
            IF (P1(2,I).EQ.P2(2,I)) THEN
              YLIMS1(I,1)=YLIMS1(I,1)-0.1
              YLIMS2(I,1)=YLIMS2(I,1)+0.1
            ENDIF
          ENDIF
C
          AT=MAX(ABS(A1LM(I)),ABS(A2LM(I)),ABS(A3LM(I)))
          IF (ABS(A1LM(I)).EQ.AT) JUMLIM(I)=1
          IF (ABS(A2LM(I)).EQ.AT) JUMLIM(I)=2
          IF (ABS(A3LM(I)).EQ.AT) JUMLIM(I)=3
          XNORM=SQRT(A1LM(I)*A1LM(I)+A2LM(I)*A2LM(I)+A3LM(I)*A3LM(I))
          IF (XNORM.LE.EPS60) GOTO 993
          A0LM(I)=A0LM(I)/XNORM
          A1LM(I)=A1LM(I)/XNORM
          A2LM(I)=A2LM(I)/XNORM
          A3LM(I)=A3LM(I)/XNORM
          JUM=JUMLIM(I)
          GOTO (91,92,93),JUM
91          ALM(I)=-A0LM(I)/A1LM(I)
            BLM(I)=-A2LM(I)/A1LM(I)
            CLM(I)=-A3LM(I)/A1LM(I)
          GOTO 97
92          ALM(I)=-A0LM(I)/A2LM(I)
            BLM(I)=-A1LM(I)/A2LM(I)
            CLM(I)=-A3LM(I)/A2LM(I)
          GOTO 97
93          ALM(I)=-A0LM(I)/A3LM(I)
            BLM(I)=-A1LM(I)/A3LM(I)
            CLM(I)=-A2LM(I)/A3LM(I)
97        CONTINUE
C
          IF (TRCSUR) THEN
            WRITE (6,*) ' A0-A3 ',A0LM(I),A1LM(I),A2LM(I),A3LM(I)
            WRITE (6,*) ' XLIMS ',XLIMS1(I,1),XLIMS2(I,1)
            WRITE (6,*) ' YLIMS ',YLIMS1(I,1),YLIMS2(I,1)
            WRITE (6,*) ' ZLIMS ',ZLIMS1(I,1),ZLIMS2(I,1)
          ENDIF
C
C  TWO POINT OPTION FINISHED. P1,P2 REDEFINED
C  ALL OTHER SURFACE COEFFICIENTS ALSO REDEFINED

C  NOW: GENERAL SECOND ORDER EQUATION, RLB=1., FOR SURFACE I
C
2       CONTINUE
1     CONTINUE
C
      RETURN
990   CONTINUE
      WRITE (6,*) 'ERROR IN SUBR. SETFIT '
      WRITE (6,*) 'INCONSISTENCY IN IGNORABLE CO-ORDINATES DETECTED'
      WRITE (6,*) 'BETWEEN REQUESTING SURFACE I= ',I,' AND IE= ',IE
      WRITE (6,*) 'JUMLIM(I),LINFX,LINFY,LINFZ ',
     .             JUMLIM(I),LINFX,LINFY,LINFZ
      WRITE (6,*) 'EXIT CALLED '
      CALL EXIT
991   CONTINUE
      WRITE (6,*) 'ERROR IN SUBR. SETFIT '
      WRITE (6,*) 'FIT OPTION FOR RLB(IE) = ',RLB(IE),' NOT FORESEEN'
      WRITE (6,*) 'REQUEST FROM SURFACE NO. ',I
      WRITE (6,*) 'IE = ',IE,' EXIT CALLED '
      CALL EXIT
992   CONTINUE
      WRITE (6,*) 'ERROR IN SUBR. SETFIT '
      WRITE (6,*) 'THE VALID AREAS DO NOT INTERSECT'
      WRITE (6,*) 'I,IE ',I,IE
      CALL EXIT
993   CONTINUE
      WRITE (6,*) 'ERROR IN SUBR. SETFIT '
      WRITE (6,*) 'STRAIGHT LINE NO I= ',I,' COLLAPSED TO A POINT'
      WRITE (6,*) 'SURFACE NO. I IS REDUNDANT. USE CH0 I/I OPTION '
      CALL EXIT
      END
C
      SUBROUTINE GRID (IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'CGEOM'
      INCLUDE 'CGRID'
      INCLUDE 'CPOLYG'
      INCLUDE 'CTRIG'
      INCLUDE 'COMUSR'
      INCLUDE 'CLOGAU'
      INCLUDE 'CADGEO'
      INCLUDE 'CTRCEI'
      INCLUDE 'CINIT'
      INCLUDE 'CCONA'
      INCLUDE 'CUPD'
      INCLUDE 'CLGIN'
      DIMENSION POLY(NPMAX),TRIAN(NTMAX)
      EQUIVALENCE (POLY(1),XPOL(1,1))
      EQUIVALENCE (TRIAN(1),XTRIAN(1))
C STATEMENT FUNCTION FOR GRID PARAMETERS FOR LEVGEO=2 OPTION
      ELPARM(R,PIN,POUT,EX1)=(PIN-POUT)*(1.-R**EX1)**1.+POUT
C
      GOTO(100,200,300),IND
C
C   RADIAL GRID
C
100   CONTINUE
C
      IF (NR1ST.LT.2) RETURN
C
      IF (LEVGEO.EQ.1) THEN
C
C  GRID DATA GENERATION FOR LEVGEO.EQ.1
C
        IF (INDGRD(IND).LE.4) THEN
C** USE ONE OF THE EIRENE DEFAULT GRID OPTIONS
          IF (RRA.GT.RAA) THEN
            NLOCAL=NR1STM
            RSURF(NR1ST)=RRA
          ELSE
            NLOCAL=NR1ST
          ENDIF
C
          ND=NLOCAL-NRSEP+1
          DO 101 J=1,ND
            RSURF(J)=RIA+DBLE((J-1))/DBLE(ND-1)*(RGA-RIA)
101       CONTINUE
          DO 102 J=ND+1,NLOCAL
            RSURF(J)=RGA+DBLE(J-ND)/DBLE(NLOCAL-ND)*(RAA-RGA)
102       CONTINUE
        ELSEIF (INDGRD(IND).EQ.5) THEN
C** TAKE RADIAL GRID DATA FROM USER SUPPLIED SUBROUTINE
C         CALL PROUSR (RSURF,2+4*NPLS+3,0.D0,0.D0,0.D0,0.D0,
C    .                 0.D0,0.D0,0.D0,NR1ST)
        ELSEIF (INDGRD(IND).EQ.6) THEN
C** TAKE RADIAL GRID DATA FROM INTERFACE
C         CALL PROFR (RSURF,6+5*NPLS+NAIN,1,1,NR1ST)
        ENDIF
C  SET DERIVED GRID DATA FOR LEVGEO = 1 OPTION
C
C  NOTHING TO BE DONE HERE
C
        IF (TRCGRD) THEN
          CALL LEER(1)
          WRITE (6,*) 'GRIDPOINTS IN X DIRECTION '
          CALL LEER(1)
          CALL MASRR1('  N, RSURF ',RSURF,NR1ST,3)
          CALL LEER(2)
        ENDIF
C
      ELSEIF (LEVGEO.EQ.2) THEN
C
C  GRID DATA GENERATION FOR LEVGEO.EQ.2
C
        IF (INDGRD(IND).LE.4) THEN
C** USE ONE OF THE EIRENE DEFAULT GRID OPTIONS
          IF (RRA.GT.RAA) THEN
            NLOCAL=NR1STM
            RSURF(NR1ST)=RRA
            ELL(NR1ST)=ELLCH
            EP1(NR1ST)=EP1CH
            TRI(NR1ST)=TRICH
          ELSE
            NLOCAL=NR1ST
          ENDIF
C
          IF (INDGRD(IND).EQ.1) THEN
            ND=NLOCAL-NRSEP+1
            DO 105 J=1,ND
              RSURF(J)=RIA+DBLE(J-1)/DBLE(ND-1)*(RGA-RIA)
105         CONTINUE
            DO 106 J=ND+1,NLOCAL
              RSURF(J)=RGA+DBLE(J-ND)/DBLE(NLOCAL-ND)*(RAA-RGA)
106         CONTINUE
          ELSEIF (INDGRD(IND).EQ.2) THEN
C** RADIAL GRID WITH CONSTANT AREA
            IF (NLCRC) THEN
              GESFL=(RAA*RAA-RIA*RIA)*PIA
              FRING=GESFL/(NLOCAL-1)
              RSURF(1)=RIA
              RSURF(NLOCAL)=RAA
              DO 103 J=2,NLOCAL-1
                RSURF(J)=RIA+SQRT((J-1)*FRING*PIAI)
103           CONTINUE
            ELSEIF (NLELL) THEN
              GESFL=(RAA*RAA*ELLOT-RIA*RIA*ELLIN)*PIA
              FRING=GESFL/(NLOCAL-1)
              RSURF(1)=RIA
              RSURF(NLOCAL)=RAA
              DO 104 J=2,NLOCAL-1
C  SOLVE RSURF**2*ELL(RSURF)-(J-1)*FRING/PIA=0., RSURF(J-1)<RSURF<RAA
                CONST=(J-1)*FRING*PIAI
                RL=RSURF(J-1)+EPS30
                RR=RAA
108             RRL=(RL-RSURF(1))/(RSURF(NLOCAL)-RSURF(1))
109             RRR=(RR-RSURF(1))/(RSURF(NLOCAL)-RSURF(1))
                FL=RL**2*ELPARM(RRL,ELLIN,ELLOT,EXELL)-CONST
                FR=RR**2*ELPARM(RRR,ELLIN,ELLOT,EXELL)-CONST
                QUOTI=(RR-RL)/(FR-FL)
                RN=-FL*QUOTI+RL
                RRN=(RN-RSURF(1))/(RSURF(NLOCAL)-RSURF(1))
                FN=RN**2*ELPARM(RRN,ELLIN,ELLOT,EXELL)-CONST
                IF (ABS(FN/CONST).LT.EPS6) THEN
                  RSURF(J)=RN
                  GOTO 104
                ELSEIF (FN.LT.0.D0) THEN
                  RL=RN
                  GOTO 108
                ELSEIF (FN.GT.0.D0) THEN
                  RR=RN
                  GOTO 109
                ENDIF
104           CONTINUE
            ELSEIF (NLTRI) THEN
              GESFL=((RAA*RAA-2.*TRIOT*TRIOT)*ELLOT-
     .               (RIA*RIA-2.*TRIIN*TRIIN)*ELLIN)*PIA
              FRING=GESFL/(NLOCAL-1)
              RSURF(1)=RIA
              RSURF(NLOCAL)=RAA
              WRITE (6,*) 'NLTRI NOT READY IN GRID'
              CALL EXIT
            ENDIF
          ELSE
            WRITE (6,*) 'INVALID OPTION ENCOUNTERED IN SUBR. GRID'
            WRITE (6,*) 'INDGRD(IND),NLCRC ',INDGRD(IND),NLCRC
            WRITE (6,*) 'EXIT CALLED FROM SUBR. GRID'
            CALL EXIT
          ENDIF
        ELSEIF (INDGRD(IND).EQ.5) THEN
C** TAKE RADIAL GRID DATA FROM USER SUPPLIED SUBROUTINE
          CALL PROUSR (RSURF,2+4*NPLS+3,0.D0,0.D0,0.D0,0.D0,
     .                 0.D0,0.D0,0.D0,NR1ST)
        ELSEIF (INDGRD(IND).EQ.6) THEN
C** TAKE RADIAL GRID DATA FROM INTERFACE
C         CALL PROFR (RSURF,2+4*NPLS+3,1,1,NR1ST)
        ENDIF
C
C  SET DATA FOR GRID OF ELLIPSES
C
        IF (INDGRD(IND).LE.2) THEN
          IF (NLCRC) THEN
            DO 111 J=1,NR1ST
              EP1(J)=0.
              ELL(J)=1.
              TRI(J)=0.
111         CONTINUE
          ELSEIF (NLELL) THEN
            EP1(1)=EP1IN
            ELL(1)=ELLIN
            DO 112 J=2,NLOCAL
              RR=(RSURF(J)-RSURF(1))/(RSURF(NLOCAL)-RSURF(1))
              EP1(J)=ELPARM(RR,EP1IN,EP1OT,EXEP1)
              ELL(J)=ELPARM(RR,ELLIN,ELLOT,EXELL)
112         CONTINUE
            DO 113 J=1,NR1ST
              TRI(J)=0.
113         CONTINUE
          ELSEIF (NLTRI) THEN
            EP1(1)=EP1IN
            ELL(1)=ELLIN
            TRI(1)=TRIIN
            IF (ABS(TRI(1)/(RSURF(1)+EPS30)).GT.EPS30) THEN
              WRITE (6,*) 'FROM SUBR. GRID: '
              WRITE (6,*) 'ERROR IN TRIANGULARITY PARAMETERS '
              CALL EXIT
            ENDIF
            DO 114 J=2,NLOCAL
              RR=(RSURF(J)-RSURF(1))/(RSURF(NLOCAL)-RSURF(1))
              EP1(J)=ELPARM(RR,EP1IN,EP1OT,EXEP1)
              ELL(J)=ELPARM(RR,ELLIN,ELLOT,EXELL)
              TRI(J)=ELPARM(RR,TRIIN,TRIOT,EXTRI)
114         CONTINUE
          ENDIF
        ELSEIF (INDGRD(IND).EQ.5) THEN
C** TAKE ELLIP. GRID DATA FROM USER SUPPLIED SUBROUTINE
C         CALL PROUSR (EP1,2+4*NPLS+?,0.,0.,0.,0.,0.,0.,0.,NR1ST)
C         CALL PROUSR (ELL,2+4*NPLS+?,0.,0.,0.,0.,0.,0.,0.,NR1ST)
C         CALL PROUSR (TRI,2+4*NPLS+?,0.,0.,0.,0.,0.,0.,0.,NR1ST)
        ELSEIF (INDGRD(IND).EQ.6) THEN
C** TAKE ELLIP. GRID DATA FROM INTERFACE
C         CALL PROFR (EP1,2+4*NPLS+?,1,1,NR1ST)
C         CALL PROFR (ELL,2+4*NPLS+?,1,1,NR1ST)
C         CALL PROFR (TRI,2+4*NPLS+?,1,1,NR1ST)
        ENDIF
C
C  SET DERIVED GRID DATA FOR LEVGEO = 2 OPTION
C  (SAME FOR ALL INDGRD OPTIONS)
C
        DO 115 J=1,NR1ST
          RQ(J)=RSURF(J)*RSURF(J)
          ELLQ(J)=ELL(J)*ELL(J)
115     CONTINUE
C
        IF (TRCGRD) THEN
          CALL MASRR4('  N, RSURF,EP1,ELL,TRI',
     .                      RSURF,EP1,ELL,TRI,NR1ST)
          CALL LEER(2)
        ENDIF
C
      ELSEIF (LEVGEO.EQ.3) THEN
C
C  GRID DATA GENERATION FOR LEVGEO.EQ.3
C
        IF (INDGRD(IND).LE.4) THEN
C  ALL POLYGON DATA HAVE BEEN READ FROM INPUT FILE, NOTHING ELSE
C  TO BE DONE HERE
        ELSEIF (INDGRD(IND).EQ.5) THEN
C*** POLYGON DATA NOT JET AVAILABLE FROM PROUSR (INDGRD.EQ.5 OPTION)
        ELSEIF (INDGRD(IND).EQ.6) THEN
c slmod begin - not tr
c...note: PLOY(1) is equivalenced to XPOL(1,1), and PROFR assigns
c         PLOY from the RWK array.
c slmod end
          CALL PROFR(POLY,6+5*NPLS+NAIN,1,1,NPMAX)
          DO 131 I=1,2
          DO 131 J=1,NPPART
            NPOINT(I,J)=XPOINT(I,J)
131       CONTINUE
        ENDIF
C
        IF (PLREFL.GT.0.D0) THEN
C  SET ONE ADDITIONAL OUTERMOST POLYGON, DISTANCE PLREFL (CM) FROM THE
C  OUTERMOST POLYGON SPECIFIED SO FAR
          I=NR1STM
          DO 133 J=1,NPPLG
            DO 132 K=NPOINT(1,J),NPOINT(2,J)
c slmod begin - grid - tr
              IF (GRIDOPT.EQ.1) THEN
                VPXX=XVERT(I-1,K,2)-XVERT(I-1,K,1)
                VPYY=YVERT(I-1,K,2)-YVERT(I-1,K,1)
              ELSE
                VPXX=XPOL(I,K)-XPOL(I-1,K)
                VPYY=YPOL(I,K)-YPOL(I-1,K)
              ENDIF
c
c              VPXX=XPOL(I,K)-XPOL(I-1,K)
c              VPYY=YPOL(I,K)-YPOL(I-1,K)
c slmod end
              XNORM=SQRT(VPXX**2+VPYY**2)
              VPX=VPXX/XNORM*PLREFL
              VPY=VPYY/XNORM*PLREFL
              XPOL(I+1,K)=XPOL(I,K)+VPX
              YPOL(I+1,K)=YPOL(I,K)+VPY
c slmod begin - grid - tr
              XVERT(I+1,K,1)=XVERT(I,K,1)+VPX
              XVERT(I+1,K,2)=XVERT(I,K,2)+VPX
              YVERT(I+1,K,1)=YVERT(I,K,1)+VPY
              YVERT(I+1,K,2)=YVERT(I,K,2)+VPY
c slmod end
132         CONTINUE
133       CONTINUE
        ENDIF
C
        IF (XPCOR.NE.0.D0) THEN
C  SHIFT WHOLE POLYGON MESH IN X DIRECTION BY XPCOR (CM)
C
          DO 135 I=1,NR1ST
            DO 135 J=1,NPPLG
              DO 135 K=NPOINT(1,J),NPOINT(2,J)
                XPOL(I,K)=XPOL(I,K)+XPCOR
c slmod begin - grid - tr
                XVERT(1,K,1)=XVERT(I,K,1)+XPCOR
                XVERT(1,K,2)=XVERT(I,K,2)+XPCOR
c slmod end
135       CONTINUE
        ENDIF
C
        IF (YPCOR.NE.0.D0) THEN
C  SHIFT WHOLE POLYGON MESH IN Y DIRECTION BY YPCOR (CM)
C
          DO 136 I=1,NR1ST
            DO 136 J=1,NPPLG
              DO 136 K=NPOINT(1,J),NPOINT(2,J)
                YPOL(I,K)=YPOL(I,K)+YPCOR
c slmod begin - grid - tr
                YVERT(I,K,1)=YVERT(I,K,1)+YPCOR
                YVERT(I,K,2)=YVERT(I,K,2)+YPCOR
c slmod end
136       CONTINUE
        ENDIF
C
C  SET DERIVED GRID DATA FOR LEVGEO = 3 OPTION
C  (SAME FOR ALL INDGRD OPTIONS)
C
        DO 140 I=1,NR1ST
          DO 140 J=1,NPPLG
            DO 140 K=NPOINT(1,J),NPOINT(2,J)-1
c slmod begin - grid - tr
              IF (GRIDOPT.EQ.1) THEN
                VPLX(I,K)=XVERT(I,K+1,1)-XVERT(I,K,1)
                VPLY(I,K)=YVERT(I,K+1,1)-YVERT(I,K,1)
              ELSE
                VPLX(I,K)=XPOL(I,K+1)-XPOL(I,K)
                VPLY(I,K)=YPOL(I,K+1)-YPOL(I,K)
              ENDIF
c
c              VPLX(I,K)=XPOL(I,K+1)-XPOL(I,K)
c              VPLY(I,K)=YPOL(I,K+1)-YPOL(I,K)
c slmod end
140     CONTINUE
        DO 141 I=1,NR1ST
          DO 141 J=1,NPPLG
            IF (J.EQ.1) THEN
              DO 142 K=1,NPOINT(1,1)
                BGL(I,K)=0.
142           CONTINUE
            ELSE
              DO 143 K=NPOINT(2,J-1),NPOINT(1,J)
                BGL(I,K)=BGL(I,NPOINT(2,J-1))
143           CONTINUE
            ENDIF
            DO 141 K=NPOINT(1,J)+1,NPOINT(2,J)
              BGL(I,K)=BGL(I,K-1)+SQRT(VPLX(I,K-1)**2+VPLY(I,K-1)**2)
141     CONTINUE
C
C   CALCULATE THE OUTER NORMALS OF POLYGONS
C
        DO 144 I=1,NR1ST
          DO 144 J=1,NPPLG
            DO 144 K=NPOINT(1,J),NPOINT(2,J)-1
              PLABS=SQRT(VPLX(I,K)**2+VPLY(I,K)**2)
              PLNX(I,K)=VPLY(I,K)/(PLABS+EPS60)
              PLNY(I,K)=-VPLX(I,K)/(PLABS+EPS60)
144     CONTINUE
C
        DO 147 I=1,NR1ST
          IUP=I+1
          IDN=I
          IF (IUP.GT.NR1ST) THEN
            IUP=I
            IDN=I-1
          ENDIF
          DO 147 J=1,NPPLG
          DO 147 K=NPOINT(1,J),NPOINT(2,J)-1
c slmod begin - grid - tr
 146        IF (GRIDOPT.EQ.1) THEN
              XD=XVERT(IDN,K+1,2)-XVERT(IDN,K+1,1)
              YD=YVERT(IDN,K+1,2)-YVERT(IDN,K+1,1)
            ELSE
              XD=XPOL(IUP,K+1)-XPOL(IDN,K+1)
              YD=YPOL(IUP,K+1)-YPOL(IDN,K+1)
            ENDIF
c
c146         XD=XPOL(IUP,K+1)-XPOL(IDN,K+1)
c            YD=YPOL(IUP,K+1)-YPOL(IDN,K+1)
c slmod end
            IF (XD*XD+YD*YD.LT.EPS30) THEN
              IF (IUP.LT.NR1ST) THEN
                IUP=IUP+1
              ELSE
                IDN=IDN-1
              ENDIF
              GOTO 146
            ENDIF
            XS=SIGN(1.D0,XD*PLNX(I,K)+YD*PLNY(I,K))
            PLNX(I,K)=PLNX(I,K)*XS
            PLNY(I,K)=PLNY(I,K)*XS
147     CONTINUE
C
        IF (TRCGRD) THEN
          WRITE (6,*) ' NO. OF VALID PARTS = ',NPPLG
          DO 155 J=1,NR1ST
            WRITE (6,*) ' POLYGON NO. J = ',J
            DO 156 K=1,NPPLG
              WRITE (6,*) 'IA = ',NPOINT(1,K),' IE = ',NPOINT(2,K)
              WRITE (6,'(/1X,1P,6E12.4)') (XPOL(J,I),YPOL(J,I),
     .                                   I=NPOINT(1,K),NPOINT(2,K))
156         CONTINUE
155       CONTINUE
          CALL LEER(2)
          WRITE (6,*) 'ARCLENGTH OF RADIAL SURFACES AT Z=0.'
          DO 153 I=1,NR1ST
            WRITE (6,*) 'I = ',I
            WRITE (6,*) ' BGL(I,K),K=1,NP2ND ',(BGL(I,K),K=1,NP2ND)
            CALL LEER(1)
153       CONTINUE
        ENDIF
C
        CALL SNEIGH
C
      ELSEIF (LEVGEO.EQ.4) THEN
C
C  GRID DATA GENERATION FOR LEVGEO.EQ.4
C
        IF (INDGRD(1).NE.6) THEN
          WRITE (6,*) ' WRONG GRID OPTION SPECIFIED FOR FEM GRID '
          CALL EXIT
        ENDIF
        CALL PROFR(TRIAN,6+5*NPLS+NAIN,1,1,NTMAX)
        DO 164 I=1,3
        DO 164 J=1,NTRII
          NECKE(I,J)=XECKE(I,J)
          NCHBAR(I,J)=XNCHBR(I,J)
          NSEITE(I,J)=SEITE(I,J)
164     CONTINUE
C
C  SET DERIVED GRID DATA FOR LEVGEO = 4 OPTION
C  (SAME FOR ALL INDGRD OPTIONS)
C
        DO 165 I=1,NTRII
          VTRIX(1,I)=XTRIAN(NECKE(2,I))-XTRIAN(NECKE(1,I))
          VTRIY(1,I)=YTRIAN(NECKE(2,I))-YTRIAN(NECKE(1,I))
          VTRIX(2,I)=XTRIAN(NECKE(3,I))-XTRIAN(NECKE(2,I))
          VTRIY(2,I)=YTRIAN(NECKE(3,I))-YTRIAN(NECKE(2,I))
          VTRIX(3,I)=XTRIAN(NECKE(1,I))-XTRIAN(NECKE(3,I))
          VTRIY(3,I)=YTRIAN(NECKE(1,I))-YTRIAN(NECKE(3,I))
165     CONTINUE
C
C
C   CALCULATE THE OUTER NORMALS OF TRIANGLES
C
        DO 161 I=1,NTRII
          PLABS1=SQRT(VTRIX(1,I)**2+VTRIY(1,I)**2)
          PLABS2=SQRT(VTRIX(2,I)**2+VTRIY(2,I)**2)
          PLABS3=SQRT(VTRIX(3,I)**2+VTRIY(3,I)**2)
          PTRIX(1,I)=VTRIY(1,I)/(PLABS1+EPS60)
          PTRIX(2,I)=VTRIY(2,I)/(PLABS2+EPS60)
          PTRIX(3,I)=VTRIY(3,I)/(PLABS3+EPS60)
          PTRIY(1,I)=-VTRIX(1,I)/(PLABS1+EPS60)
          PTRIY(2,I)=-VTRIX(2,I)/(PLABS2+EPS60)
          PTRIY(3,I)=-VTRIX(3,I)/(PLABS3+EPS60)
161     CONTINUE
C
        DO 162 I=1,NTRII
          XD1=VTRIX(1,I)+PTRIX(1,I)
          YD1=VTRIY(1,I)+PTRIY(1,I)
          XS1=SIGN(1.D0,XD1*VTRIX(1,I)+YD1*VTRIY(1,I))
          PTRIX(1,I)=PTRIX(1,I)*XS1
          PTRIY(1,I)=PTRIY(1,I)*XS1
          XD2=VTRIX(2,I)+PTRIX(2,I)
          YD2=VTRIY(2,I)+PTRIY(2,I)
          XS2=SIGN(1.D0,XD2*VTRIX(2,I)+YD2*VTRIY(2,I))
          PTRIX(2,I)=PTRIX(2,I)*XS2
          PTRIY(2,I)=PTRIY(2,I)*XS2
          XD3=VTRIX(3,I)+PTRIX(3,I)
          YD3=VTRIY(3,I)+PTRIY(3,I)
          XS3=SIGN(1.D0,XD3*VTRIX(3,I)+YD3*VTRIY(3,I))
          PTRIX(3,I)=PTRIX(3,I)*XS3
          PTRIY(3,I)=PTRIY(3,I)*XS3
162     CONTINUE
C
        IF (TRCGRD) THEN
          WRITE (6,*) ' NUMBER OF TRIANGLES = ',NTRII
          WRITE (6,*) ' I,(XTRIAN(J),YTRIAN(J),J=1,3) '
          DO 163 I=1,NTRII
            WRITE (6,'(/1X,I4,1X,1P,6E12.4)')
     .                               I,(XTRIAN(NECKE(J,I)),
     .                                  YTRIAN(NECKE(J,I)),J=1,3)
163       CONTINUE
          CALL LEER(2)
        ENDIF
C
      ELSEIF (LEVGEO.EQ.5) THEN
C
C  GENERAL GEOMETRY OPTION: NOTHING TO DONE HERE
C
      ENDIF
C
C  SET GEOMETRICAL CONTANTS FOR Y OR POLOIDAL CO-ORDINATE
C
      IF (LEVGEO.EQ.1) THEN
        YDF=YAA-YIA
        PSURF(1)=YIA
      ELSEIF (LEVGEO.EQ.2) THEN
        YDF=(YAA-YIA)*DEGRAD
        PSURF(1)=YIA
      ELSEIF (LEVGEO.EQ.3) THEN
        YDF=1.
        PSURF(1)=YIA
      ELSEIF (LEVGEO.EQ.4) THEN
        YDF=1.
        PSURF(1)=YIA
      ELSEIF (LEVGEO.EQ.5) THEN
        YDF=1.
C
C  GENERAL GEOMETRY OPTION: NOTHING TO BE DONE HERE
C
      ENDIF
      IF (YDF.LE.0.D0) GOTO 992
C
      IF (TRCGRD.AND..NOT.NLPOL) THEN
        CALL LEER(2)
        WRITE (6,*) 'CONSTANTS FOR POLOIDAL OR Y DIRECTION'
        CALL MASR3('YDF,YIA,YAA=            ',YDF,YIA,YAA)
        CALL LEER(1)
      ENDIF
C
C  SET GEOMETRICAL CONTANTS FOR Z CO-ORDINATE
C
C  A) IN TOROIDAL APPROXIMATION:
C
      IF (NLTRA) THEN
        IF (NTTRA.LE.3.OR.ROA.LE.0.D0) GOTO 991
        IF (LEVGEO.EQ.1) THEN
          XDIFF=0.
        ELSEIF (LEVGEO.EQ.2) THEN
          XDIFF=EP1OT
        ELSEIF (LEVGEO.EQ.3) THEN
          XDIFF=0.
        ELSEIF (LEVGEO.EQ.4) THEN
          XDIFF=0.
        ELSEIF (LEVGEO.EQ.5) THEN
C
C  GENERAL GEOMETRY OPTION: NOTHING TO BE DONE HERE
C
        ENDIF
C
        ZDF=(ZAA-ZIA)*DEGRAD
C
C  ALPHA: HALF OF THE ANGLE INCREMENT IN EQUIDISTANT TOROIDAL ANGLE GRID
        ALPHA=0.5*(ZDF/DBLE(NTTRAM))
        TANAL=TAN(ALPHA)
        SINAL=SIN(2.*ALPHA)
        COSAL=COS(2.*ALPHA)
C
        DPHI=1./(2.*ALPHA)
C
C  ROA IS THE LARGE RADIUS OF THE TORUS
C  (RMTOR,0,0) IS THE ORIGIN OF LOCAL CO-ORDINATE SYSTEM IN
C              EACH TOROIDAL CELL
C  RMTOR SUCH THAT VOLUME OF TORUS = VOLUME OF THE NTTRAM SEGMENTS
C  AT PRESENT: FULLFILLED AT SURFACE DEFINED BY (RAA,EP1OT,ELLOT)
C              OR AT A POLYGON WITH XDIFF=0. (IF THERE IS ONE)
        RMTOR=(ROA+XDIFF)*ALPHA/TANAL-XDIFF
        ZHALF=ALPHA
        ZFULL=ZHALF*2.
C  SET ZSURF EVEN IF NLTOR=FALSE, FOR 3D GEOMETRY PLOTS
        DO 170 J=1,NTTRA
          ZSURF(J)=ZIA*DEGRAD+(J-1)/DPHI
170     CONTINUE
        DO 172 J=1,NTTRAM
          ZZONE(J)=0.5*(ZSURF(J)+ZSURF(J+1))
172     CONTINUE
        RORIG=RMTOR
C
C  B) IN CYLIND. APPROXIMATION:
C     ROA AND RMTOR ARE IRRELEVANT IN THIS CASE, AND ARE NOT DEFINED
C
      ELSEIF (NLTRZ) THEN
        ZDF=ZAA-ZIA
        IF (ZDF.LE.0.D0) GOTO 991
        ZSURF(1)=ZIA
        ZZONE(1)=(ZAA+ZIA)*0.5
        DPHI=1./ZDF
C
        RORIG=0.
C
C  C) IN TORUS CO-ORDINATES
C
      ELSEIF (NLTRT) THEN
        ZDF=(ZAA-ZIA)*DEGRAD
        IF (ZDF.LE.0.D0) GOTO 991
        ZSURF(1)=ZIA
        ZZONE(1)=(ZAA+ZIA)*0.5
        DPHI=1./ZDF
C
        RORIG=0.
C
      ELSE
        WRITE (6,*) ' ERROR IN INPUT DATA! '
        WRITE (6,*) ' NLTRA OR NLTRZ OR NLTRT MUST BE .TRUE. '
        CALL EXIT
      ENDIF
C
      IF (TRCGRD) THEN
        CALL LEER(2)
        IF (.NOT.NLTOR) THEN
          WRITE (6,*) 'CONSTANTS FOR TOROIDAL OR Z DIRECTION'
          CALL MASR3('ZDF,ZIA,ZAA=            ',ZDF,ZIA,ZAA)
        ENDIF
        IF (NLTRA) THEN
          WRITE (6,*) 'ROA,RMTOR= ',ROA,RMTOR
C         IF (.NOT.NLTOR) THEN
C           CALL MASRR1 (' N,  ZSURF ',ZSURF,NTTRA,3)
C           CALL MASRR1 (' N,  ZZONE ',ZZONE,NTTRAM,3)
C         ENDIF
          CALL LEER(2)
        ENDIF
        CALL LEER(1)
      ENDIF
C
C  SET SURFACE AREA OF NON DEFAULT STANDARD SURFACES
C
      DO 180 ISTS=1,NSTSI
        IF (INUMP(ISTS,1).NE.0) THEN
          NLJ=NLIM+ISTS
          IF (NLTRZ) THEN
            SAREA(NLJ)=YDF*ZDF
          ELSEIF (NLTRA) THEN
          ELSEIF (NLTRT) THEN
          ENDIF
        ENDIF
180   CONTINUE
C
      RETURN
C
C   POLOIDAL OR Y-GRID
C
200   CONTINUE
C
C  IF NLSYMP, Y-GRID MUST BE SYMMETRIC
C
      IF (LEVGEO.EQ.1) THEN
C   Y-GRID
        IF (INDGRD(IND).EQ.1) THEN
          ND=NP2ND-NPSEP+1
          DO 205 J=1,ND
            PSURF(J)=YIA+DBLE((J-1))/DBLE(ND-1)*(YGA-YIA)
205       CONTINUE
          DO 206 J=ND+1,NP2ND
            PSURF(J)=YGA+DBLE(J-ND)/DBLE(NP2ND-ND)*(YAA-YGA)
206       CONTINUE
        ELSEIF (INDGRD(IND).EQ.2) THEN
          DO 207 J=1,NP2ND
            PSURF(J)=YIA+(J-1)/DBLE(NP2NDM)*(YAA-YIA)
207       CONTINUE
C       ELSEIF (INDGRD(IND).EQ.5) THEN
C** TAKE Y GRID DATA FROM USER SUPPLIED SUBROUTINE
C TO BE WRITTEN
C         CALL PROUSR (PSURF,2+4*NPLS+3,0.,0.,0.,0.,0.,0.,0.,NP2ND)
C       ELSEIF (INDGRD(IND).EQ.6) THEN
C** TAKE Y GRID DATA FROM INTERFACE
C TO BE WRITTEN
C         CALL PROFR (PSURF,2+4*NPLS+5,1,1,NP2ND)
        ENDIF
C
        DO 210 J=1,NP2NDM
210       PHZONE(J)=(PSURF(J)+PSURF(J+1))/2.
C
        IF (TRCGRD) THEN
          CALL LEER(1)
          WRITE (6,*) 'GRIDPOINTS IN Y DIRECTION '
          CALL LEER(1)
          CALL MASRR1('  N, PSURF ',PSURF,NP2ND,3)
          CALL LEER(2)
        ENDIF
C
      ELSEIF (LEVGEO.EQ.2) THEN
C
        IF (INDGRD(IND).EQ.1) THEN
          ND=NP2ND-NPSEP+1
          DO 211 J=1,ND
            PSURF(J)=(YIA+DBLE((J-1))/DBLE(ND-1)*(YGA-YIA))*DEGRAD
211       CONTINUE
          DO 212 J=ND+1,NP2ND
            PSURF(J)=(YGA+DBLE(J-ND)/DBLE(NP2ND-ND)*(YAA-YGA))*DEGRAD
212       CONTINUE
        ELSEIF (INDGRD(IND).EQ.2) THEN
          DO 213 J=1,NP2ND
            PSURF(J)=(YIA+(J-1)/DBLE(NP2NDM)*(YAA-YIA))*DEGRAD
213       CONTINUE
        ENDIF
        DO 215 J=1,NP2ND
          COSPH(J)=COS(PSURF(J))
          SINPH(J)=SIN(PSURF(J))
215     CONTINUE
C
        IF (TRCGRD) THEN
          CALL LEER(1)
          WRITE (6,*) 'GRIDPOINTS IN POLOIDAL DIRECTION '
          CALL LEER(1)
          CALL MASRR1('  N, PSURF ',PSURF,NP2ND,3)
          CALL LEER(2)
        ENDIF
C
C
        NPPLG=1
        NPOINT(1,1)=1
        NPOINT(2,1)=NP2ND
        IFLAG=2
c slmod begin - grid - not tr (XPOL and YPOL references are absent in EIRENE02)
        STOP 'GRID: XPOL and YPOL are assigned - needs development'
c slmod end
        DO 1240 IR=1,NR1STM
          IRP=IR+1
          DO 1250 IP=1,NP2NDM
            IPP=IP+1
            CALL ARELLP(EP1(IRP),EP1(IR),ELL(IRP),ELL(IR),
     .                  TRI(IRP),TRI(IR),
     .                  RSURF(IRP),RSURF(IR),PSURF(IPP),PSURF(IP),IFLAG,
     .                  AELL,SX,SY,X1,Y1,X2,Y2,X3,Y3,X4,Y4)
C
            XPOL(IRP,IPP)=X1
            YPOL(IRP,IPP)=Y1
            XPOL(IR,IPP)=X2
            YPOL(IR,IPP)=Y2
            XPOL(IRP,IP)=X3
            YPOL(IRP,IP)=Y3
            XPOL(IR,IP)=X4
            YPOL(IR,IP)=Y4
1250      CONTINUE
1240    CONTINUE
C
      ENDIF
C
      IF (LEVGEO.EQ.2.OR.LEVGEO.EQ.3) THEN
C
        IF (TRCGRD) THEN
          DO 219 I=1,NP2ND
            WRITE (6,*) ' PERP. POLYGON NO. I = ',I
            WRITE (6,*) ' JA = ',1,' JE = ',NR1ST
            WRITE (6,'(/1X,1P,6E12.4)') (XPOL(K,I),YPOL(K,I),
     .             K=1,NR1ST)
219       CONTINUE
        ENDIF
C
c slmod begin - tr
        WRITE(6,*) 'NOTE: EXTENDED GEOMETRY INFORMATION NOT '//
     .             'LISTED'
c slmod end
        DO 220 K=1,NP2ND
          DO 220 I=1,NR1STM
c slmod begin - grid - tr
            IF (GRIDOPT.EQ.1) THEN
              VVTX(I,K)=XVERT(I,K,2)-XVERT(I,K,1)
              VVTY(I,K)=YVERT(I,K,2)-YVERT(I,K,1)
            ELSE
              VVTX(I,K)=XPOL(I+1,K)-XPOL(I,K)
              VVTY(I,K)=YPOL(I+1,K)-YPOL(I,K)
            ENDIF
c
c            VVTX(I,K)=XPOL(I+1,K)-XPOL(I,K)
c            VVTY(I,K)=YPOL(I+1,K)-YPOL(I,K)
c slmod end
220     CONTINUE
C
        DO 221 K=1,NP2ND
          BGLP(1,K)=0.
          DO 222 I=1,NR1STM
            BGLP(I+1,K)=BGLP(I,K)+SQRT(VVTX(I,K)**2+VVTY(I,K)**2)
222       CONTINUE
221     CONTINUE
C
        IF (TRCGRD) THEN
          CALL LEER(2)
          WRITE (6,*) 'ARCLENGTH OF POLOIDAL SURFACES AT Z=0.'
          DO 223 K=1,NP2ND
            WRITE (6,*) 'K = ',K
            WRITE (6,*) ' BGLP(I,K),I=1,NR1ST ',(BGLP(I,K),I=1,NR1ST)
            CALL LEER(1)
223       CONTINUE
        ENDIF
C
C   CALCULATE THE OUTER NORMALS OF POLYGONS
C
        DO 224 K=1,NP2ND
          DO 225 I=1,NR1STM
            IF (ABS(VVTY(I,K)).LT.EPS12) THEN
              PPLNX(I,K)=0.
              PPLNY(I,K)=1.
            ELSE
              PPLNX(I,K)=1.
              PPLNY(I,K)=-VVTX(I,K)/VVTY(I,K)
            ENDIF
225       CONTINUE
          DO 224 I=1,NR1STM
            PLABS=SQRT(PPLNX(I,K)**2+PPLNY(I,K)**2)
            PPLNX(I,K)=PPLNX(I,K)/PLABS
            PPLNY(I,K)=PPLNY(I,K)/PLABS
224     CONTINUE
C
        DO 227 I=1,NR1STM
        DO 227 J=1,NPPLG
          DO 227 K=NPOINT(1,J),NPOINT(2,J)
            KUP=K+1
            KDN=K
            IF (KUP.GT.NPOINT(2,J)) THEN
              KUP=K
              KDN=K-1
            ENDIF
c slmod begin - grid - tr
 226        IF (GRIDOPT.EQ.1) THEN
c             IF (i.EQ.23.AND.K.EQ.NPOINT(2,J)) THEN
c               WRITE(6,*) 'SEARCH:',k,kup,kdn
c             ENDIF
              XD=XVERT(I,KUP,2)-XVERT(I,KDN,2)
              YD=YVERT(I,KUP,2)-YVERT(I,KDN,2)
            ELSE
              XD=XPOL(I+1,KUP)-XPOL(I+1,KDN)
              YD=YPOL(I+1,KUP)-YPOL(I+1,KDN)
            ENDIF
c
c226         XD=XPOL(I+1,KUP)-XPOL(I+1,KDN)
c            YD=YPOL(I+1,KUP)-YPOL(I+1,KDN)
c slmod end
            IF (XD*XD+YD*YD.LT.EPS30) THEN
c slmod begin - new
c SPECIAL CASE - CRAP, I NEED TO GET NPPLG=1 AND AVOID THE NEED
c FOR THIS PATCH.
              IF     (J.EQ.1    .AND.KUP.GE.NPOINT(2,J)-1) THEN
                IF (KUP.EQ.NPOINT(2,J)-1) THEN
                  KUP=NPOINT(1,J+1)
                ELSE
                  KUP=KUP+1
                ENDIF
              ELSEIF (J.EQ.NPPLG.AND.KUP.EQ.NPOINT(2,J)) THEN   
                IF (KDN.EQ.NPOINT(2,J-1)) THEN
                  KDN=NPOINT(2,J-1)-1
                ELSE
                  KDN=KDN-1
                ENDIF
              ELSEIF (KUP.LT.NPOINT(2,J)) THEN
c
c              IF (KUP.LT.NPOINT(2,J)) THEN
c slmod end
                KUP=KUP+1
              ELSE
                KDN=KDN-1
              ENDIF
              GOTO 226
            ENDIF
            XS=SIGN(1.D0,XD*PPLNX(I,K)+YD*PPLNY(I,K))
            PPLNX(I,K)=PPLNX(I,K)*XS
            PPLNY(I,K)=PPLNY(I,K)*XS
227     CONTINUE
C
C  IDENTIFY DEAD CELLS IN GRID CUTS
C
        IT=1
        DO 229 IR=1,NR1ST-1
        DO 229 J=1,NPPLG-1
          DO 229 IP=NPOINT(2,J),NPOINT(1,J+1)-1
            NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
            NSTGRD(NCELL)=2
229     CONTINUE
C
      ELSEIF (LEVGEO.GT.3) THEN
C
        WRITE (6,*) 'ERROR EXIT FROM GRID. NLPOL ',LEVGEO
      ENDIF
C
C  1ST AND 2ND GRID DEFINED
C  SET SURFACE AREA OF NON DEFAULT STANDARD SURFACES
C
C  RADIAL (1ST GRID) SURFACES. OVERWRITE EARLIER VALUES FROM
C  CALL GRID(1)
C
      IF (LEVGEO.EQ.1) THEN
        DO 280 ISTS=1,NSTSI
          IF (INUMP(ISTS,1).NE.0) THEN
            IR=INUMP(ISTS,1)
            NLJ=NLIM+ISTS
            IF (NLTRZ) THEN
              SAREA(NLJ)=PSURF(IRPTE(ISTS,2))-PSURF(IRPTA(ISTS,2))
              SAREA(NLJ)=SAREA(NLJ)*ZDF
            ELSEIF (NLTRA) THEN
            ELSEIF (NLTRT) THEN
            ENDIF
          ENDIF
280     CONTINUE
      ELSEIF (LEVGEO.EQ.2.OR.LEVGEO.EQ.3) THEN
        DO 290 ISTS=1,NSTSI
          IF (INUMP(ISTS,1).NE.0) THEN
            IR=INUMP(ISTS,1)
            NLJ=NLIM+ISTS
            IF (NLTRZ) THEN
              SAREA(NLJ)=BGL(IR,IRPTE(ISTS,2))-BGL(IR,IRPTA(ISTS,2))
              SAREA(NLJ)=SAREA(NLJ)*ZDF
c slmod begin - debug - not tr


c slmod end
            ELSEIF (NLTRA) THEN
              SAREA(NLJ)=0.
              DO 291 IP=IRPTA(ISTS,2),IRPTE(ISTS,2)-1
c slmod begin - grid - tr
                IF (GRIDOPT.EQ.1) THEN
                  XS=(XVERT(IR,IP+1,1)+XVERT(IR,IP,1))*0.5+RMTOR
                ELSE
                  XS=((XPOL(IR,IP+1)+XPOL(IR,IP))*0.5)+RMTOR
                ENDIF
c
c                XS=((XPOL(IR,IP+1)+XPOL(IR,IP))*0.5)+RMTOR
c slmod end
c slmod begin - not tr
c BUG! THIS LINE WAS LEFT IN BY ACCIDENT - FEB 12, 2003
c                XS=((XPOL(IR,IP+1)+XPOL(IR,IP))*0.5)+RMTOR
c slmod end
                SAREA(NLJ)=SAREA(NLJ)+(BGL(IR,IP+1)-BGL(IR,IP))*XS
291           CONTINUE
              SAREA(NLJ)=SAREA(NLJ)*TANAL/ALPHA*PI2A
            ENDIF
          ENDIF
290     CONTINUE
      ELSE
C TO BE WRITTEN
      ENDIF
C
C  POLOIDAL (2ND GRID) SURFACES.
C
      IF (LEVGEO.EQ.1) THEN
        DO 285 ISTS=1,NSTSI
          IF (INUMP(ISTS,2).NE.0) THEN
            IP=INUMP(ISTS,2)
            NLJ=NLIM+ISTS
            IF (NLTRZ) THEN
              SAREA(NLJ)=RSURF(IRPTE(ISTS,1))-RSURF(IRPTA(ISTS,1))
              SAREA(NLJ)=SAREA(NLJ)*ZDF
            ELSEIF (NLTRA) THEN
              SAREA(NLJ)=0.
              DO 286 IR=IRPTA(ISTS,1),IRPTE(ISTS,1)-1
                XS=(RSURF(IR+1)+RSURF(IR))*0.5+RMTOR
                SAREA(NLJ)=SAREA(NLJ)+(RSURF(IR+1)-RSURF(IR))*XS
286           CONTINUE
              SAREA(NLJ)=SAREA(NLJ)*TANAL/ALPHA*PI2A
            ELSEIF (NLTRT) THEN
              SAREA(NLJ)=0.
              DO 287 IR=IRPTA(ISTS,1),IRPTE(ISTS,1)-1
                XS=(RSURF(IR+1)+RSURF(IR))*0.5*2.*PIA
                SAREA(NLJ)=SAREA(NLJ)+(RSURF(IR+1)-RSURF(IR))*XS
287           CONTINUE
            ENDIF
          ENDIF
285     CONTINUE
      ELSEIF (LEVGEO.EQ.2.OR.LEVGEO.EQ.3) THEN
        DO 295 ISTS=1,NSTSI
          IF (INUMP(ISTS,2).NE.0) THEN
            IP=INUMP(ISTS,2)
            NLJ=NLIM+ISTS
            IF (NLTRZ) THEN
              SAREA(NLJ)=BGLP(IRPTE(ISTS,1),IP)-BGLP(IRPTA(ISTS,1),IP)
              SAREA(NLJ)=SAREA(NLJ)*ZDF
            ELSEIF (NLTRA) THEN
              SAREA(NLJ)=0.
              DO 296 IR=IRPTA(ISTS,1),IRPTE(ISTS,1)-1
c slmod begin - grid - tr
                IF (GRIDOPT.EQ.1) THEN
                  XS=(XVERT(IR,IP,2)+XVERT(IR,IP,1))*0.5+RMTOR
                ELSE
                  XS=(XPOL(IR+1,IP)+XPOL(IR,IP))*0.5+RMTOR
                ENDIF
c
c                XS=(XPOL(IR+1,IP)+XPOL(IR,IP))*0.5+RMTOR
c slmod end
                SAREA(NLJ)=SAREA(NLJ)+(BGLP(IR+1,IP)-BGLP(IR,IP))*XS
296           CONTINUE
              SAREA(NLJ)=SAREA(NLJ)*TANAL/ALPHA*PI2A
            ELSEIF (NLTRT) THEN
              SAREA(NLJ)=0.
              DO 297 IR=IRPTA(ISTS,1),IRPTE(ISTS,1)-1
c slmod begin - grid - tr
                IF (GRIDOPT.EQ.1) THEN
                  XS=(XVERT(IR,IP,2)+XVERT(IR,IP,1))*PIA
                ELSE
                  XS=(XPOL(IR+1,IP)+XPOL(IR,IP))*0.5*2.*PIA
                ENDIF
c
c                XS=(XPOL(IR+1,IP)+XPOL(IR,IP))*0.5*2.*PIA
c slmod end
                SAREA(NLJ)=SAREA(NLJ)+(BGLP(IR+1,IP)-BGLP(IR,IP))*XS
297           CONTINUE
            ENDIF
          ENDIF
295     CONTINUE
      ELSE
C TO BE WRITTEN
      ENDIF
C
      RETURN
C
C   TOROIDAL OR Z-GRID
C
300   CONTINUE
C
C  IF NLSYMT, Z-GRID MUST BE SYMMETRIC
C
      IF (NLTRZ) THEN
C   Z-GRID
        IF (INDGRD(IND).EQ.1) THEN
          ND=NT3RD-NTSEP+1
          DO 305 J=1,ND
            ZSURF(J)=ZIA+DBLE((J-1))/DBLE(ND-1)*(ZGA-ZIA)
305       CONTINUE
          DO 306 J=ND+1,NT3RD
            ZSURF(J)=ZGA+DBLE(J-ND)/DBLE(NT3RD-ND)*(ZAA-ZGA)
306       CONTINUE
        ELSEIF (INDGRD(IND).EQ.2) THEN
          DO 307 J=1,NT3RD
            ZSURF(J)=ZIA+(J-1)/DBLE(NT3RDM)*(ZAA-ZIA)
307       CONTINUE
C       ELSEIF (INDGRD(IND).EQ.5) THEN
C** TAKE Z GRID DATA FROM USER SUPPLIED SUBROUTINE
C TO BE WRITTEN
C         CALL PROUSR (ZSURF,2+4*NPLS+3,0.,0.,0.,0.,0.,0.,0.,NT3RD)
C       ELSEIF (INDGRD(IND).EQ.6) THEN
C** TAKE Z GRID DATA FROM INTERFACE
C TO BE WRITTEN
C         CALL PROFR (ZSURF,2+4*NPLS+5,1,1,NT3RD)
        ENDIF
C
        DO 310 J=1,NT3RDM
310       ZZONE(J)=(ZSURF(J)+ZSURF(J+1))/2.
C
C     ELSEIF (NLTRA) THEN
C   GRID FOR TOROIDAL APPROXIMATION OF CYLINDER: ALREADY DONE
C
      ENDIF
C
      IF (TRCGRD) THEN
        CALL LEER(1)
        WRITE (6,*) 'GRIDPOINTS IN Z DIRECTION'
        CALL LEER(1)
        CALL MASRR1 (' N,  ZSURF ',ZSURF,NT3RD,3)
        CALL LEER(2)
      ENDIF
C
      RETURN
C
991   CONTINUE
      WRITE (6,*) 'GRID DATA INCONSISTENCY: 3RD GRID.  ZAA > ZIA ?'
      WRITE (6,*) 'ZIA,ZAA,NTTRA,ROA= ',ZIA,ZAA,NTTRA,ROA
      CALL EXIT
992   CONTINUE
      WRITE (6,*) 'GRID DATA INCONSISTENCY: 2ND GRID.  YAA > YIA ?'
      WRITE (6,*) 'YIA,YAA = ',YIA,YAA
      CALL EXIT
993   CONTINUE
      WRITE (6,*) 'GRID DATA INCONSISTENCY: 1ST GRID.  RAA > RIA ?'
      WRITE (6,*) 'RIA,RAA = ',RIA,RAA
      CALL EXIT
      END
C
      SUBROUTINE PLASMA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'CGRID'
      INCLUDE 'COMUSR'
      INCLUDE 'CINIT'
      INCLUDE 'CCONA'
      DIMENSION PLPRM(NPLPRM)
      EQUIVALENCE (PLPRM(1),TEIN(1))
      DIMENSION HELP(NRAD)
C
      DO 5 J=1,NPLPR1
        PLPRM(J)=0.
5     CONTINUE
C
C  SET EIRENE VACUUM MODEL DATA. I.E. IF ALL TEMPERATURES ARE
C  LESS THAN TVAC OR THE BACKGROUND DENSITY IS LESS THAN DVAC,
C  THEN THIS ZONE IS CONSIDERED TO BE AN "EIRENE VACUUM ZONE",
C  PARTICLE MEAN FREE PATHES IN SUCH ZONES ARE
C  EQUAL TO 1.D10 (CM) AND REACTION RATES ARE ZERO
      TVAC=0.01
      DVAC=1.D2
      VVAC=0.
C
C     SET DENSITY, TEMPERATURE AND MACH NUMBER PROFILES
C     ON MESH "RHOZNE(J)"
C
C  ELECTRON TEMPERATURE
      IND=INDPRO(1)
      GOTO (101,102,103,104,105,106,107),IND
101     CALL PROFN (TEIN,TE0,TE1,TE2,TE3,TE4,TE5,TVAC)
        GOTO 110
102     CALL PROFE (TEIN,TE0,TE1,TE2,TE4,TE5,TVAC)
        GOTO 110
103     CALL PROFS (TEIN,TE0,TE1,TE5,TVAC)
        GOTO 110
C  INDPRO=4 NOT IN USE
104     CONTINUE
        GOTO 110
105     CALL PROUSR (TEIN,0,TE0,TE1,TE2,TE3,TE4,TE5,TVAC,NSURF)
        GOTO 110
106     CALL PROFR (TEIN,0,1,1,NSURF)
        GOTO 110
107     CALL PROFR (TEIN,0,1,1,NSBOX)
        GOTO 110
110   CONTINUE
C  ION TEMPERATURE
      IND=INDPRO(2)
      DO 120 IPLS=1,NPLSI
        GOTO (111,112,113,114,115,116,117),IND
111     CALL PROFN (HELP,TI0(IPLS),TI1(IPLS),TI2(IPLS),
     .                   TI3(IPLS),TI4(IPLS),TI5(IPLS),TVAC)
        CALL RESETP (TIIN,HELP,IPLS,1,NPLS,NR1ST)
        GOTO 120
112     CALL PROFE (HELP,TI0(IPLS),TI1(IPLS),TI2(IPLS),
     .                   TI4(IPLS),TI5(IPLS),TVAC)
        CALL RESETP (TIIN,HELP,IPLS,1,NPLS,NR1ST)
        GOTO 120
113     CALL PROFS (HELP,TI0(IPLS),TI1(IPLS),TI5(IPLS),TVAC)
        CALL RESETP (TIIN,HELP,IPLS,1,NPLS,NR1ST)
        GOTO 120
C  INDPRO=4 NOT IN USE
114     CONTINUE
        GOTO 120
115     CALL PROUSR (HELP,1+0*NPLS,TI0(IPLS),TI1(IPLS),TI2(IPLS),
     .                      TI3(IPLS),TI4(IPLS),TI5(IPLS),TVAC,NSURF)
        CALL RESETP (TIIN,HELP,IPLS,1,NPLS,NSURF)
        GOTO 120
120     CONTINUE
        GOTO 1120
116     CALL PROFR (TIIN,1+0*NPLS,NPLSI,NPLS,NSURF)
        GOTO 1120
117     CALL PROFR (TIIN,1+0*NPLS,NPLSI,NPLS,NSBOX)
        GOTO 1120
1120  CONTINUE
C  ION DENSITY
      IND=INDPRO(3)
      DO 130 K=1,NPLSI
        GOTO (121,122,123,124,125,126,127),IND
121       CALL PROFN (HELP,DI0(K),DI1(K),DI2(K),DI3(K),DI4(K),DI5(K),
     .                DVAC)
          CALL RESETP (DIIN,HELP,K,1,NPLS,NR1ST)
          GOTO 130
122       CALL PROFE (HELP,DI0(K),DI1(K),DI2(K),DI4(K),DI5(K),
     .                DVAC)
          CALL RESETP (DIIN,HELP,K,1,NPLS,NR1ST)
          GOTO 130
123       CALL PROFS (HELP,DI0(K),DI1(K),DI5(K),DVAC)
          CALL RESETP (DIIN,HELP,K,1,NPLS,NR1ST)
          GOTO 130
C  INDPRO=4 NOT IN USE
124     CONTINUE
          GOTO 130
125       CALL PROUSR (HELP,1+1*NPLS,DI0(K),DI1(K),DI2(K),DI3(K),DI4(K),
     .                 DI5(K),DVAC,NSURF)
          CALL RESETP (DIIN,HELP,K,1,NPLS,NSURF)
          GOTO 130
130     CONTINUE
        GOTO 1130
126     CALL PROFR (DIIN,1+1*NPLS,NPLSI,NPLS,NSURF)
        GOTO 1130
127     CALL PROFR (DIIN,1+1*NPLS,NPLSI,NPLS,NSBOX)
        GOTO 1130
1130  CONTINUE
C  DRIFT VELOCITY
      IND=INDPRO(4)
      DO 140 IPLS=1,NPLSI
        GOTO (131,132,133,134,135,136,137),IND
131     CALL PROFN (HELP,VX0(IPLS),VX1(IPLS),VX2(IPLS),
     .                   VX3(IPLS),VX4(IPLS),VX5(IPLS),VVAC)
        CALL RESETP (VXIN,HELP,IPLS,1,NPLS,NR1ST)
        CALL PROFN (HELP,VY0(IPLS),VY1(IPLS),VY2(IPLS),
     .                   VY3(IPLS),VY4(IPLS),VY5(IPLS),VVAC)
        CALL RESETP (VYIN,HELP,IPLS,1,NPLS,NR1ST)
        CALL PROFN (HELP,VZ0(IPLS),VZ1(IPLS),VZ2(IPLS),
     .                   VZ3(IPLS),VZ4(IPLS),VZ5(IPLS),VVAC)
        CALL RESETP (VZIN,HELP,IPLS,1,NPLS,NR1ST)
        GOTO 140
132     CALL PROFE (HELP,VX0(IPLS),VX1(IPLS),VX2(IPLS),
     .                   VX4(IPLS),VX5(IPLS),VVAC)
        CALL RESETP (VXIN,HELP,IPLS,1,NPLS,NR1ST)
        CALL PROFE (HELP,VY0(IPLS),VY1(IPLS),VY2(IPLS),
     .                   VY4(IPLS),VY5(IPLS),VVAC)
        CALL RESETP (VYIN,HELP,IPLS,1,NPLS,NR1ST)
        CALL PROFE (HELP,VZ0(IPLS),VZ1(IPLS),VZ2(IPLS),
     .                   VZ4(IPLS),VZ5(IPLS),VVAC)
        CALL RESETP (VZIN,HELP,IPLS,1,NPLS,NR1ST)
        GOTO 140
133     CALL PROFS (HELP,VX0(IPLS),VX1(IPLS),VX5(IPLS),VVAC)
        CALL RESETP (VXIN,HELP,IPLS,1,NPLS,NR1ST)
        CALL PROFS (HELP,VY0(IPLS),VY1(IPLS),VY5(IPLS),VVAC)
        CALL RESETP (VYIN,HELP,IPLS,1,NPLS,NR1ST)
        CALL PROFS (HELP,VZ0(IPLS),VZ1(IPLS),VZ5(IPLS),VVAC)
        CALL RESETP (VZIN,HELP,IPLS,1,NPLS,NR1ST)
        GOTO 140
C  INDPRO=4 NOT IN USE
134     CONTINUE
        GOTO 140
135     CALL PROUSR (HELP,1+2*NPLS,VX0(IPLS),VX1(IPLS),VX2(IPLS),
     .                        VX3(IPLS),VX4(IPLS),VX5(IPLS),VVAC,NSURF)
        CALL RESETP (VXIN,HELP,IPLS,1,NPLS,NSURF)
        CALL PROUSR (HELP,1+3*NPLS,VY0(IPLS),VY1(IPLS),VY2(IPLS),
     .                        VY3(IPLS),VY4(IPLS),VY5(IPLS),VVAC,NSURF)
        CALL RESETP (VYIN,HELP,IPLS,1,NPLS,NSURF)
        CALL PROUSR (HELP,1+4*NPLS,VZ0(IPLS),VZ1(IPLS),VZ2(IPLS),
     .                        VZ3(IPLS),VZ4(IPLS),VZ5(IPLS),VVAC,NSURF)
        CALL RESETP (VZIN,HELP,IPLS,1,NPLS,NSURF)
        GOTO 140
140   CONTINUE
C  SCALE FROM MACH NUMBER PROFILE TO CM/SEC PROFILE?
      IF (NLMACH) THEN
        DO 1141 IPLS=1,NPLSI
          DO 1142 ICELL=1,NSURF
            FACT=CVEL2A*SQRT((TIIN(IPLS,ICELL)+
     .                        TEIN(ICELL))/RMASSP(IPLS))
            VXIN(IPLS,ICELL)=VXIN(IPLS,ICELL)*FACT
            VYIN(IPLS,ICELL)=VYIN(IPLS,ICELL)*FACT
            VZIN(IPLS,ICELL)=VZIN(IPLS,ICELL)*FACT
1142      CONTINUE
1141    CONTINUE
      ENDIF
      GOTO 1140
136   CALL PROFR (VXIN,1+2*NPLS,NPLSI,NPLS,NSURF)
      CALL PROFR (VYIN,1+3*NPLS,NPLSI,NPLS,NSURF)
      CALL PROFR (VZIN,1+4*NPLS,NPLSI,NPLS,NSURF)
      GOTO 1140
137   CALL PROFR (VXIN,1+2*NPLS,NPLSI,NPLS,NSBOX)
      CALL PROFR (VYIN,1+3*NPLS,NPLSI,NPLS,NSBOX)
      CALL PROFR (VZIN,1+4*NPLS,NPLSI,NPLS,NSBOX)
      GOTO 1140
1140  CONTINUE
C  MAGNETIC FIELD UNIT VECTOR
      IND=INDPRO(5)
      GOTO (141,141,141,141,145,146,147),IND
C
C  DEFAULT: POSITIVE Z DIRECTION
141   CONTINUE
      DO 142 J=1,NSURF
        BXIN(J)=0.
        BYIN(J)=0.
        BZIN(J)=1.
142   CONTINUE
      GOTO 150
145   CONTINUE
      CALL PROUSR (BXIN,1+5*NPLS,BX0,BX1,BX2,BX3,BX4,BX5,0.D0,NSURF)
      CALL PROUSR (BYIN,2+5*NPLS,BY0,BY1,BY2,BY3,BY4,BY5,0.D0,NSURF)
      CALL PROUSR (BZIN,3+5*NPLS,BZ0,BZ1,BZ2,BZ3,BZ4,BZ5,1.D0,NSURF)
      GOTO 150
146   CALL PROFR (BXIN,1+5*NPLS,1,1,NSURF)
      CALL PROFR (BYIN,2+5*NPLS,1,1,NSURF)
      CALL PROFR (BZIN,3+5*NPLS,1,1,NSURF)
      GOTO 150
147   CALL PROFR (BXIN,1+5*NPLS,1,1,NSBOX)
      CALL PROFR (BYIN,2+5*NPLS,1,1,NSBOX)
      CALL PROFR (BZIN,3+5*NPLS,1,1,NSBOX)
      GOTO 150
150   CONTINUE
C
C  CHECK FOR ZERO MAGNETIC FIELD IN STANDARD MESH
      DO 153 IR=1,NR1STM
      DO 153 IP=1,NP2NDM
      DO 153 IT=1,NT3RDM
        JJ=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
        IF (BXIN(JJ)**2+BYIN(JJ)**2+BZIN(JJ)**2.LE.EPS30) THEN
          WRITE (6,*) 'ZERO B-FIELD IN STANDARD CELL JJ= ',JJ
          CALL EXIT
        ENDIF
153   CONTINUE
C
C  ADDITIONAL INPUT TALLIES
      IND=INDPRO(6)
      DO 160 K=1,NAINI
        GOTO (151,151,151,151,155,156,157),IND
C  DEFAULT: ZERO
151     CONTINUE
        DO 1151 J=1,NR1ST
          ADIN(K,J)=0.
1151    CONTINUE
        GOTO 160
155     CALL PROUSR (HELP,6+5*NPLS,BD,BD,BD,BD,BD,BD,0.D0,NSURF)
        CALL RESETP (ADIN,HELP,K,1,NAIN,NSURF)
        GOTO 160
160   CONTINUE
      GOTO 1160
156   CALL PROFR (ADIN,6+5*NPLS,NAINI,NAIN,NSURF)
      GOTO 1160
157   CALL PROFR (ADIN,6+5*NPLS,NAINI,NAIN,NSBOX)
      GOTO 1160
1160  CONTINUE
C
C
C   SET VACUUM DATA IN ADDITIONAL REGIONS OUTSIDE THE
C   THE STANDARD MESH
C
      IF (INDPRO(1).LE.6) THEN
        DO J=NSURF+1,NSURF+NRADD
          TEIN(J)=TVAC
        ENDDO
      ENDIF
      IF (INDPRO(2).LE.6) THEN
        DO J=NSURF+1,NSURF+NRADD
          DO 17 IPLS=1,NPLSI
c slmod begin - not tr
c            IF (IPLS.GT.2.AND.J.GT.NSURF+1) THEN
c              TIIN(IPLS,J)=2.0
c            ELSE
c              TIIN(IPLS,J)=TVAC
c            ENDIF
c
            TIIN(IPLS,J)=TVAC
c slmod end
17        CONTINUE
        ENDDO
      ENDIF
      IF (INDPRO(3).LE.6) THEN
        DO J=NSURF+1,NSURF+NRADD
          DO 18 IPLS=1,NPLSI
c slmod begin - not tr
c            IF (IPLS.GT.2.AND.J.GT.NSURF+1) THEN
c              DIIN(IPLS,J)=1.0E+12
c            ELSE
c              DIIN(IPLS,J)=DVAC
c            ENDIF
c
            DIIN(IPLS,J)=DVAC
c slmod end
18        CONTINUE
        ENDDO
      ENDIF
      IF (INDPRO(4).LE.6) THEN
        DO J=NSURF+1,NSURF+NRADD
          DO 19 IPLS=1,NPLSI
            VXIN(IPLS,J)=VVAC
            VYIN(IPLS,J)=VVAC
            VZIN(IPLS,J)=VVAC
19        CONTINUE
        ENDDO
      ENDIF
      IF (INDPRO(5).LE.6) THEN
        DO J=NSURF+1,NSURF+NRADD
          BXIN(J)=0.
          BYIN(J)=0.
          BZIN(J)=1.
        ENDDO
      ENDIF
      IF (INDPRO(6).LE.6) THEN
        DO J=NSURF+1,NSURF+NRADD
          DO 20 IAIN=1,NAINI
            ADIN(IAIN,J)=0.
20        CONTINUE
        ENDDO
      ENDIF
C
      RETURN
      END
