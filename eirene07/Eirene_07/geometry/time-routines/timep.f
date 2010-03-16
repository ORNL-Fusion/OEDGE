C
C
      SUBROUTINE TIMEP (ZRAD)
C
C   CALCULATE TIME SEGMENTS IN Y- OR POLOIDAL MESH
C
C   INPUT:
C
C       X00,Y00,Z00: STARTING POINT FOR THIS TRACK
C       ZRAD   =
C       NRCELL = RADIAL CELL NUMBER, IN WHICH THIS TRACK OF LENGTH ZRAD
C                                    IS PERFORMED
C              = 0, IF TRACK OUTSIDE RADIAL MESH
C                   IF MRSURF .NE. 0 THEN
C                   REENTRY FOUND AT RADIAL SURFACE MRSURF.
C                   OTHERWISE: REENTRY AT NON DEFAULT
C                   POLOIDAL SURFACES IS SEARCHED FOR IN THIS CALL
C       NTCELL =
C       ITCELL =
C       NPCELL = POLOIDAL CELL INDEX OF POINT X00,Y00,Z00
C                AT WHICH THIS TRACK OF LENGTH ZRAD STARTS
C  WARNING: IF NLSRFY, NPCELL MAY BE WRONG
C       IPCELL =
C
C       NCOUT,BLPD,NCOUNT
C
C   OUTPUT:
C
C       X00,Y00,Z00: END POINT FOR THIS TRACK
C       NCOU   = TOTAL NUMBER OF CELLS, WHICH THE CURRENT TRACK OF
C                LENGTH ZRAD HAS CROSSED
C       CLPD(I)= LENGTH OF THE PART NO. I (I=1,NCOU)
C       JUPC(I)= CELL NUMBER IN P-GRID FOR TRACK I
C       LUPC(I)= SURFACE NUMBER IN P-GRID AN THE END OF TRACK I
C       MUPC(I)= ORIENTATION OF TRACK I
C       NPCELL = LAST POLOIDAL CELL INDEX OF X00,Y00,Z00 POINT,
C                ON WHICH THIS TRACK OF LENGTH ZRAD ENDS
C       IPCELL = NUMBER OF FINAL POLOIDAL CELL NPCELL
C                IF TRACK ORIGINATED INSIDE STANDARD MESH
C              = NUMBER OF POLOIDAL CELL, AT WHICH REENTRY WAS FOUND
C                ON RADIAL OR TOROIDAL SURFACE
C                IF TRACK ORIGINATED OUTSIDE  STANDARD MESH
C       IRCELL = NUMBER OF RADIAL CELL, AT WHICH REENTRY AT POLOIDAL
C                SURFACE MPSURF WAS FOUND. OTHERWISE: UNMODIFIED.
C
C              = 0  IF TRACK COMPLETELY OUTSIDE STANDARD MESH
C
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CLOGAU
      USE CUPD
      USE CPOLYG
      USE CGRID
      USE CGEOM
      USE COMPRT
      USE COMSOU
      USE CLGIN

      IMPLICIT NONE

      REAL(DP), INTENT(INOUT) :: ZRAD
      REAL(DP) :: HELP, GS, GC, F, XNEN, DUM, SUM, BB, DY, X0TEST,
     .          Y0TEST, V1, V2, T1, ZRADS, ZRD, WIN1, X000, Z0T,
     .          TWIN1, XT, Y000, TSS, X0T, Y0T
      INTEGER :: ILLZ, IANP, I, IENP, JHELP, ISW, JJC, J1, NRCLLP,
     .           LHELP, J2, INCY, JN, ICOU, NYSAVE, NHELP,
     .           LEARCA, I1, MPTEST, IR, IRSAVE, NCOUPE, ISTS, IADD,
     .           ICOUT, NCPEN, IPOLGS, NCPAN, J, LEARC2, NJC, LEARC1,
     .           IST, ITEST, JSH, NN ,IN
      INTEGER :: NCOUNS(N2ND+N3RD)
      LOGICAL :: LCUTY(N2NDPLG), LCUTX(N1ST), lnincz
      SAVE
C
C     IF (NLTRC) THEN
C       CALL LEER(1)
C       IF (NRCELL.GT.0) THEN
C         WRITE (iunout,*) 'TIMEP FROM INSIDE, NPANU ', NPANU
C         WRITE (iunout,*) 'ZRAD,NRCELL,NPCELL '
C         WRITE (iunout,*)  ZRAD,NRCELL,NPCELL
C       ELSE
C         WRITE (iunout,*) 'TIMEP FROM OUTSIDE, NPANU ', NPANU
C         WRITE (iunout,*) 'MRSURF,MTSURF,ZRAD '
C         WRITE (iunout,*)  MRSURF,MTSURF,ZRAD
C       ENDIF
C     ENDIF
C
      ICOUT=1
      IADD=0
      ZRADS=ZRAD
      IPOLGS=IPOLGN
      ncpan=0
      ncpen=0
C
      IF (NLTOR) THEN
C       NCOUT=NCOUT
        ZRAD=BLPD(1)
        NTCELL=NCOUNT(1)
        IF (NCOUT.GT.1) IPOLGN=0
C       IF (NLTRC) 
C     .   WRITE (iunout,*) 'WG. NLTOR: ZRAD,NTCELL ',ZRAD,NTCELL
      ELSE
        NCOUT=1
C       ZRAD=ZRAD
C       NTCELL=1
      ENDIF
C
10000 CONTINUE
C
CDR NOV.99: ADDED BECAUSE OF FOLION OPTION, WITH DEFAULT B-FIELD (IN Z DIRECTION
C
      IF (ABS(VELZ).EQ.1.D0) THEN
        NCOUP=1
        ALPD(NCOUP)=ZRAD
        JUPC(NCOUP)=NPCELL
C       X00=X00+ZRAD*VELX
C       Y00=Y00+ZRAD*VELY
        IPCELL=NPCELL
        MPSURF=0
        GOTO 5000
      ENDIF
C
      IF (LEVGEO.EQ.3) THEN
C
        IF (NRCELL.GT.0) GOTO 10
C
C  PARTICLE OUTSIDE STANDARD MESH
C  CHECK AT NON DEFAULT POLOIDAL SURFACES
C
        NCOUP=0
        ZRD=ZRAD
        DO 3 ISTS=1,NSTSI
          MPTEST=INUMP(ISTS,2)
          IF (MPTEST.NE.0) THEN
C  TEST POLOIDAL SURFACE NO. MPTEST FOR REENTRY
            DO 2 IR=1,NR1STM
              I1=IR+1
              V1=(YPOL(IR,MPTEST)-Y00)*VELX-(XPOL(IR,MPTEST)-X00)*VELY
              V2=(YPOL(I1,MPTEST)-Y00)*VELX-(XPOL(I1,MPTEST)-X00)*VELY
              LCUTX(IR)=V1*V2.LE.0.
2           CONTINUE
            DO 4 IR=1,NR1STM
              IF (LCUTX(IR)) THEN
                T1=((XPOL(IR,MPTEST)-X00)*VVTY(IR,MPTEST)-
     .              (YPOL(IR,MPTEST)-Y00)*VVTX(IR,MPTEST))
     .             /(VELX*VVTY(IR,MPTEST)-VELY*VVTX(IR,MPTEST)+EPS60)
C               IF (NLTRC) WRITE (iunout,*) 'IR,MPTEST,T1 ',IR,MPTEST,T1
                IF (T1.LT.0..OR.T1.GE.ZRD) GOTO 4
C               IF (NLTRC) 
C     .           WRITE (iunout,*) 'VALID INTERSECTION AT T1= ',T1
                NCOUP=1
                LUPC(NCOUP)=MPTEST
                IRSAVE=IR
                MUPC(NCOUP)=SIGN(1._DP,VELX*PPLNX(IR,MPTEST)+
     .                                 VELY*PPLNY(IR,MPTEST))
                JUPC(NCOUP)=1
                ZRD=T1
              ENDIF
4           CONTINUE
          ENDIF
3       CONTINUE
        IF (NCOUP.GT.0.AND.LUPC(MAX(1,NCOUP)).NE.MPSURF) THEN
C  REENTRY FOUND, REDUCE ZRAD TO T1
C  NCOUP=1 AT THIS POINT
          NCOUPE=1
          MPSURF=LUPC(NCOUPE)
          IPOLGN=LUPC(NCOUPE)
          NINCY=MUPC(NCOUPE)
          IRCELL=IRSAVE
          ZRAD=ZRD
          ISRFCL=0
          ALPD(NCOUPE)=ZRAD
          MRSURF=0
          MTSURF=0
          MASURF=0
          NINCX=0
          NINCZ=0
C         IF (NLTRC) THEN
C           WRITE (iunout,*) 'REENTRY FOUND, MPSURF,ZRAD = ',MPSURF,ZRAD
C           WRITE (iunout,*) 'IRCELL ',IRCELL
C         ENDIF
          GOTO 31
        ELSE
C  NO REENTRY FOUND
          NCOUP=1
          JUPC(1)=1
          ALPD(1)=ZRAD
          IPCELL=IPOLGN
          MPSURF=0
C         IF (NLTRC) THEN
C           WRITE (iunout,*) 'NO REENTRY INTO POLOIDAL GRID FOUND '
C         ENDIF
          X00=X00+ZRAD*VELX
          Y00=Y00+ZRAD*VELY
          GOTO 5000
        ENDIF
C
10      CONTINUE
C
C  PARTICLE INSIDE STANDARD MESH, RADIAL CELL NO. NRCELL
C
        IF (NCOUP.EQ.0) THEN
          WRITE (iunout,*) 'ERROR IN TIMEP: NCOUP=0'
          RETURN
        ENDIF
C       IF (NLTRC) THEN
C         WRITE (iunout,*) ' TIMEP IN NEIGHBOR PART '
C         WRITE (iunout,*) ' ALPD ',(ALPD(J),J=1,NCOUP)
C         WRITE (iunout,*) ' JUPC ',(JUPC(J),J=1,NCOUP)
C         WRITE (iunout,*) ' LUPC ',(LUPC(J),J=1,NCOUP)
C         WRITE (iunout,*) ' MUPC ',(MUPC(J),J=1,NCOUP)
C       ENDIF
        NLSRFY=.FALSE.
        IF ((icout == 1) .and.
     .    (SQRT((X0-X00)**2+(Y0-Y00)**2).GT.EPS10)) THEN
          DO J=1,NCOUP
            ALPD(J)=ALPD(J)-ZT
          ENDDO
        ENDIF
C
C  ACCOUNT FOR OTHER SURFACES INSIDE POLOIDAL MESH
C
        lnincz=.false.
        if (ityp==3) lnincz=(zrad < alpd(1))
        IF (ALPD(NCOUP).GT.ZRAD) THEN
          DO J=1,NCOUP
            IF (ALPD(J).GT.ZRAD) THEN
              ncpan=j+1
              ncpen=ncoup+1
              do jsh=ncoup,j+1,-1
                alpd(jsh+1)=alpd(jsh)-zrad
                jupc(jsh+1)=jupc(jsh)
                lupc(jsh+1)=lupc(jsh)
                mupc(jsh+1)=mupc(jsh)
              end do
              alpd(ncpan)=alpd(j)-zrad
              jupc(ncpan)=jupc(j)
              lupc(ncpan)=lupc(j)
              mupc(ncpan)=mupc(j)
              ALPD(J)=ZRAD
              IPOLGN=JUPC(J)
              NCOUP=J
              LUPC(NCOUP)=0
              MUPC(NCOUP)=0
              GOTO 4711
            ENDIF
          ENDDO
4711      CONTINUE
        ENDIF
C
C   ADJUST ALPD AND ACCOUNT FOR "NONDEFAULT" POLOIDAL SURFACES
C
        TSS=0.
        IST=0
        DO J=1,NCOUP
          ALPD(J)=ALPD(J)-TSS
          TSS=TSS+ALPD(J)
C
          ITEST=INMP2I(NRCELL,LUPC(J),0)
          IN=ITEST+NLIM
!pb          IF (ITEST.NE.0.AND.ILIIN(IN).NE.0) THEN
          IF ((ityp==3).or.(ITEST.NE.0.AND.ILIIN(IN).NE.0)) THEN
C
C  TRACK ENDS ON ONE OF THE NON DEFAULT POLOIDAL SURFACES
C
C           IF (NLTRC) THEN
C             WRITE (iunout,*) ' TRACK TERMINATED'
C             WRITE (iunout,*) ' ITEST,ILIIN ',ITEST,ILIIN(IN)
C           ENDIF
            NCOUPE=J
            MPSURF=LUPC(NCOUPE)
            IPOLGN=LUPC(NCOUPE)
            NINCY=MUPC(NCOUPE)
            ZRAD=TSS
            if ((ITEST.NE.0.AND.ILIIN(IN).NE.0).or.
     .          (1._dp-zrad/(zrads+eps60) > eps10)) then
              ISRFCL=0
              MASURF=0
            end if
!pb
            if (lnincz) then         ! toroidal surface
              nincx=0
              nincy=0
              mrsurf=0
              mpsurf=0
            elseif ((ityp==3).and.(nincx.ne.0)) then   ! radial surface
              nincy=0
              nincz=0
              mpsurf=0
              mtsurf=0
            else                     ! poloidal surface
              nincx=0
              nincz=0
              mrsurf=0
              mtsurf=0
            end if
!pb
            NCOUP=NCOUPE
            GOTO 311
          ENDIF
C
        ENDDO
C
C  LAST CELL, TRACK DOES NOT END ON A POLOIDAL SURFACE
C
C  INDEX IPOLGN OF LAST CELL NOT KNOWN ?
        IF (IPOLGN.EQ.0) THEN
          X0T=X00+VELX*ZRAD
          Y0T=Y00+VELY*ZRAD
          NN=LEARC1(X0T,Y0T,Z0T,IPOLGN,
     .              NRCELL,NRCELL,.FALSE.,.FALSE.,
     .              NPANU,'TIMEP       ')
        ENDIF
C
        MPSURF=0
        NINCY=0
C
311     CONTINUE
C       IF (NLTRC) THEN
C         WRITE (iunout,*) ' ALPD ',(ALPD(J),J=1,NCOUP)
C         WRITE (iunout,*) ' JUPC ',(JUPC(J),J=1,NCOUP)
C         WRITE (iunout,*) ' LUPC ',(LUPC(J),J=1,NCOUP)
C         WRITE (iunout,*) ' MUPC ',(MUPC(J),J=1,NCOUP)
C       ENDIF
C
        X00=X00+ZRAD*VELX
        Y00=Y00+ZRAD*VELY
        NPCELL=JUPC(NCOUP)
        IPCELL=NPCELL
        GOTO 5000
C
C
30      CONTINUE
C
C  LAST CELL, TRACK DOES NOT END ON A POLOIDAL SURFACE
C
C  INDEX OF LAST CELL NOT KNOWN (E.G. DUE TO ADD. SURFACE) ?
CDR  ERROR: FALLS NLSRFY, GGFLS NPCELL FALSCH
        IF (IPOLGN.EQ.0) THEN
CDR       IF (NCOUP.EQ.0) THEN
CDR         IPOLGN=NPCELL
CDR       ELSE
            X0T=X00+VELX*ZRAD
            Y0T=Y00+VELY*ZRAD
            IPOLGN=LEARC2(X0T,Y0T,NRCELL,NPANU,'TIMEP       ')
CDR       ENDIF
        ENDIF
C
        NCOUPE=NCOUP+1
        JUPC(NCOUPE)=IPOLGN
        LUPC(NCOUPE)=0
        MUPC(NCOUPE)=0
        ALPD(NCOUPE)=ZRAD-TSS
        MPSURF=0
        NINCY=0
C
31      CONTINUE
        NCOUP=NCOUPE
C       IF (NLTRC) THEN
C         WRITE (iunout,*) ' ALPD ',(ALPD(J),J=1,NCOUP)
C         WRITE (iunout,*) ' JUPC ',(JUPC(J),J=1,NCOUP)
C         WRITE (iunout,*) ' LUPC ',(LUPC(J),J=1,NCOUP)
C         WRITE (iunout,*) ' MUPC ',(MUPC(J),J=1,NCOUP)
C       ENDIF
C
        X00=X00+ZRAD*VELX
        Y00=Y00+ZRAD*VELY
        NPCELL=JUPC(NCOUP)
        IPCELL=NPCELL
        GOTO 5000
C
      ELSEIF (LEVGEO.EQ.2) THEN
        IF (NLCRC) THEN
C
C
C  DISTANCE TO NEXT X- OR RADIAL SURFACE KNOWN?
C
C       IF (TS.GE.1.D30) GOTO 1000
C
C  YES! ZRAD IS THE DISTANCE TRAVELLED IN X- OR RADIAL CELL NO. NRCELL
C
          NCOUP=1
          IF (NRCELL.LT.1.OR.NRCELL.GT.NR1STM) THEN
            X00=X00+VELX*ZRAD
            Y00=Y00+VELY*ZRAD
            WIN1=MOD(ATAN2(Y00,X00)+PI2A-PSURF(1),PI2A)+PSURF(1)
            NPCELL=WIN1/YDF*DBLE(NP2NDM)+1.
            GOTO 5000
          ENDIF
C
C  THE OLD POLOIDAL CELL INDEX IS: NPCELL
C  FIND THE NEW CELL INDEX : NJC
C
          X000=X00+VELX*ZRAD
          Y000=Y00+VELY*ZRAD
          WIN1=MOD(ATAN2(Y000,X000)+PI2A-PSURF(1),PI2A)+PSURF(1)
          NJC=WIN1/YDF*DBLE(NP2NDM)+1.
C
          TWIN1=0.
          NCOUP=0
          IF (NJC.EQ.NPCELL) GOTO 150
C
C   FIND ORIENTATION IN THETA-GRID
C
          XT=-Y00*VELX+X00*VELY
          NINCY=1
          IF (XT.LT.0.) NINCY=-1
C
C   CONTRIBUTION TO EACH THETA-CELL
C   NPCELL : STARTINDEX
C   JJC    : SURFACEINDEX
C   J1     : CELLINDEX
C   NJC    : ENDINDEX
          J1=NPCELL
100       JJC=J1
          IF (NINCY.EQ.1) JJC=JJC+1
C   TIMESTEP FROM X00,Y00 TO THETA-SURFACE, THETA=WIN
          GS=SINPH(JJC)
          GC=COSPH(JJC)
          XNEN=VELX*GS-VELY*GC
          F=(Y00*GC-X00*GS)/(XNEN+EPS60)
          NCOUP=NCOUP+1
          JUPC(NCOUP)=J1
          LUPC(NCOUP)=JJC
          ALPD(NCOUP)=F-TWIN1
          TWIN1=F
C
          J1=J1+NINCY
          IF (J1.EQ.0) J1=NP2NDM
          IF (J1.EQ.NP2ND) J1=1
!pb          IF (J1.NE.NJC) GOTO 100
          IF ((ityp.ne.3).and.(J1.NE.NJC)) GOTO 100
C
C   LAST THETA-CELL
C
150       NCOUP=NCOUP+1
          JUPC(NCOUP)=NJC
          ALPD(NCOUP)=ZRAD-TWIN1
          X00=X000
          Y00=Y000
          NPCELL=NJC
          MPSURF=0
          NINCY=0
          GOTO 5000
C
        ELSE
C
C  NEW PART: NOT NLCRC
C  PARTICLE INSIDE STANDARD MESH
C
          NCOUP=0
          NRCLLP=NRCELL+1
C
C   SEARCH FOR ALL POSSIBLE INTERSECTIONS WITHIN THE RADIAL CELL NRCELL
C
          DO 111 I=1,NP2ND
111         LCUTY(I)=.FALSE.
C
          DO 112 J=1,NP2ND
            V1=(YPOL(NRCELL,J)-Y00)*VELX-(XPOL(NRCELL,J)-X00)*VELY
            V2=(YPOL(NRCLLP,J)-Y00)*VELX-(XPOL(NRCLLP,J)-X00)*VELY
            LCUTY(J)=V1*V2.LE.0.
C           IF (NLTRC) THEN
C             IF (LCUTY(J)) WRITE (iunout,*) 'LCUTY=TRUE FOR ',J
C           ENDIF
112       CONTINUE
          IF (NLSRFY) THEN
            LCUTY(MPSURF)=.FALSE.
            NLSRFY=.FALSE.
C  PSURF(1)=PSURF(NP2ND)
            IF (MPSURF.EQ.1) LCUTY(NP2ND)=.FALSE.
            IF (MPSURF.EQ.NP2ND) LCUTY(1)=.FALSE.
          ENDIF
C
          IANP=ILLZ(NP2ND,LCUTY,1)+1
          IENP=NP2ND-ILLZ(NP2ND,LCUTY,-1)
C
C   COMPUTE THE FLIGHT TIMES TO THE INTERSECTION POINTS
C
        DO 114 I=IANP,IENP
          IF (LCUTY(I)) THEN
            T1=((XPOL(NRCELL,I)-X00)*VVTY(NRCELL,I)-
     .          (YPOL(NRCELL,I)-Y00)*VVTX(NRCELL,I))
     .         /(VELX*VVTY(NRCELL,I)-VELY*VVTX(NRCELL,I)+EPS60)
C           IF (NLTRC) WRITE (iunout,*) 'I,T1 ',I,T1
            IF (T1.LT.0..OR.T1.GE.ZRAD) GOTO 114
C           IF (NLTRC) WRITE (iunout,*) 'VALID INTERSECTION AT T1= ',T1
            NCOUP=NCOUP+1
            LUPC(NCOUP)=I
            MUPC(NCOUP)=SIGN(1._DP,VELX*PPLNX(NRCELL,I)+
     .                             VELY*PPLNY(NRCELL,I))
            IF (MUPC(NCOUP).EQ.-1) THEN
              JUPC(NCOUP)=I
              ALPD(NCOUP)=T1
              IF (I.EQ.NP2ND) NCOUP=NCOUP-1
            ELSEIF (MUPC(NCOUP).EQ.1) THEN
              JUPC(NCOUP)=I-1
              ALPD(NCOUP)=T1
              IF (I.EQ.1) NCOUP=NCOUP-1
            ENDIF
          ENDIF
114     CONTINUE
C
C   REARRANGE THE FLIGHT TIMES IN ASCENDING ORDER
C
115     ISW=0
        DO 120 J=1,NCOUP-1
          IF (ALPD(J).GT.ALPD(J+1)) THEN
            ISW=ISW+1
            HELP=ALPD(J)
            ALPD(J)=ALPD(J+1)
            ALPD(J+1)=HELP
            JHELP=JUPC(J)
            JUPC(J)=JUPC(J+1)
            JUPC(J+1)=JHELP
            LHELP=LUPC(J)
            LUPC(J)=LUPC(J+1)
            LUPC(J+1)=LHELP
            NHELP=MUPC(J)
            MUPC(J)=MUPC(J+1)
            MUPC(J+1)=NHELP
          ENDIF
120     CONTINUE
        IF (ISW.GT.0.AND.NCOUP.GT.2) GOTO 115
C       IF (NLTRC.AND.NCOUP.GT.0) THEN
C         WRITE (iunout,*) ' NACH SORTIEREN '
C         WRITE (iunout,*) ' ALPD ',(ALPD(J),J=1,NCOUP)
C         WRITE (iunout,*) ' JUPC ',(JUPC(J),J=1,NCOUP)
C         WRITE (iunout,*) ' LUPC ',(LUPC(J),J=1,NCOUP)
C         WRITE (iunout,*) ' MUPC ',(MUPC(J),J=1,NCOUP)
C       ENDIF
C
        DO 125 J=1,NCOUP-1
          IF (ABS(ALPD(J+1)-ALPD(J)).LE.EPS30) THEN
            IF (JUPC(J).LE.0.OR.JUPC(J).GE.NP2ND) THEN
C             IF (NLTRC) THEN
C               WRITE (iunout,*) ' VERTAUSCHE ALPD(',J,') UND (',J+1,')'
C               WRITE (iunout,*) ' ALPD = ',ALPD(J),ALPD(J+1)
C               WRITE (iunout,*) ' JUPC = ',JUPC(J),JUPC(J+1)
C             ENDIF
              HELP=ALPD(J)
              ALPD(J)=ALPD(J+1)
              ALPD(J+1)=HELP
              JHELP=JUPC(J)
              JUPC(J)=JUPC(J+1)
              JUPC(J+1)=JHELP
              LHELP=LUPC(J)
              LUPC(J)=LUPC(J+1)
              LUPC(J+1)=LHELP
              NHELP=MUPC(J)
              MUPC(J)=MUPC(J+1)
              MUPC(J+1)=NHELP
            ENDIF
          ENDIF
125     CONTINUE
C
C   ADJUST ALPD AND ACCOUNT FOR "NONDEFAULT" POLOIDAL SURFACES
C
        TSS=0.
        IST=0
        DO 130 J=1,NCOUP
          ALPD(J)=ALPD(J)-TSS
          TSS=TSS+ALPD(J)
C
          ITEST=INMP2I(NRCELL,LUPC(J),0)
          IN=ITEST+NLIM
!pb          IF (ITEST.NE.0.AND.ILIIN(IN).NE.0) THEN
          IF ((ityp==3).or.(ITEST.NE.0.AND.ILIIN(IN).NE.0)) THEN
C
C  TRACK ENDS ON ONE OF THE NON DEFAULT POLOIDAL SURFACES
C
C           IF (NLTRC) THEN
C             WRITE (iunout,*) ' TRACK TERMINATED'
C             WRITE (iunout,*) ' ITEST,ILIIN ',ITEST,ILIIN(IN)
C           ENDIF
            NCOUPE=J
            MPSURF=LUPC(NCOUPE)
            IPOLGN=LUPC(NCOUPE)
            NINCY=MUPC(NCOUPE)
            ZRAD=TSS
            ISRFCL=0
            NINCX=0
            NINCZ=0
            MRSURF=0
            MTSURF=0
            MASURF=0
            GOTO 131
          ENDIF
C
130     CONTINUE
C
C  LAST CELL, TRACK DOES NOT END ON A POLOIDAL SURFACE
C
C  INDEX OF LAST CELL NOT KNOWN (E.G. DUE TO ADD. SURFACE) ?
CDR  ERROR: FALLS NLSRFY, GGFLS NPCELL FALSCH
        IF (IPOLGN.EQ.0) THEN
CDR       IF (NCOUP.EQ.0) THEN
CDR         IPOLGN=NPCELL
CDR       ELSE
            X0T=X00+VELX*ZRAD
            Y0T=Y00+VELY*ZRAD
            IPOLGN=LEARC2(X0T,Y0T,NRCELL,NPANU,'TIMEP       ')
CDR       ENDIF
        ENDIF
C
        NCOUPE=NCOUP+1
        JUPC(NCOUPE)=IPOLGN
        LUPC(NCOUPE)=0
        MUPC(NCOUPE)=0
        ALPD(NCOUPE)=ZRAD-TSS
        MPSURF=0
        NINCY=0
C
131       CONTINUE
          NCOUP=NCOUPE
C         IF (NLTRC) THEN
C           WRITE (iunout,*) ' ALPD ',(ALPD(J),J=1,NCOUP)
C           WRITE (iunout,*) ' JUPC ',(JUPC(J),J=1,NCOUP)
C           WRITE (iunout,*) ' LUPC ',(LUPC(J),J=1,NCOUP)
C           WRITE (iunout,*) ' MUPC ',(MUPC(J),J=1,NCOUP)
C         ENDIF
C
          X00=X00+ZRAD*VELX
          Y00=Y00+ZRAD*VELY
          NPCELL=JUPC(NCOUP)
          IPCELL=NPCELL
          GOTO 5000
C
        ENDIF
C
      ELSEIF (LEVGEO.EQ.1) THEN
C
C  IDENTICAL, UP TO NAMES, TO TOROIDAL GRID PART, NLTRZ BLOCK
C
        IF (NRCELL.GT.0) GOTO 2900
C
C  PARTICLE OUTSIDE STANDARD MESH
C  CHECK AT NON DEFAULT POLOIDAL SURFACES
C
        NCOUP=0
        ZRD=ZRAD
        BB=VELY+EPS60
        NYSAVE=1
        IF (VELY.LT.0.D0) NYSAVE=-1
        DO 2903 ISTS=1,NSTSI
          MPTEST=INUMP(ISTS,2)
          IF (MPTEST.NE.0) THEN
C  TEST POLOIDAL SURFACE NO. MPTEST FOR REENTRY
C  TIME FROM Y00 TO PSURF
            DY=PSURF(MPTEST)-Y00
            F=DY/BB
C           IF (NLTRC) WRITE (iunout,*) 'MPTEST,F,DY ',MPTEST,F,DY
            IF (F.LE.ZRD.AND.F.GT.0.D0) THEN
              X0TEST=X00+VELX*F
              IF (X0TEST.GE.RSURF(1).AND.X0TEST.LE.RSURF(NR1ST)) THEN
                IRSAVE=LEARCA(X0TEST,RSURF,1,NR1ST,1,'TIMEP 1    ')
                NCOUP=1
                JUPC(NCOUP)=1
                LUPC(NCOUP)=MPTEST
                MUPC(NCOUP)=NYSAVE
                ZRD=F
              ENDIF
            ENDIF
          ENDIF
2903    CONTINUE
        IF (NCOUP.GT.0.AND.LUPC(MAX(1,NCOUP)).NE.MPSURF) THEN
C  REENTRY FOUND, REDUCE ZRAD TO F
C  NCOUP=1 AT THIS POINT
          NCOUP=1
          MPSURF=LUPC(NCOUP)
          NINCY=MUPC(NCOUP)
          IRCELL=IRSAVE
          ZRAD=ZRD
          ISRFCL=0
          ALPD(NCOUP)=ZRAD
          MRSURF=0
          MTSURF=0
          MASURF=0
          NINCX=0
          NINCZ=0
C         IF (NLTRC) THEN
C           WRITE (iunout,*) 'REENTRY FOUND, MPSURF,ZRAD = ',MPSURF,ZRAD
C           WRITE (iunout,*) 'IRCELL ',IRCELL
C         ENDIF
          Y00=Y00+VELY*ZRD
          NPCELL=JUPC(NCOUP)
          IPCELL=NPCELL
          GOTO 5000
        ELSE
C  NO REENTRY FOUND
          NCOUP=1
          JUPC(1)=1
          ALPD(1)=ZRAD
          IF (MRSURF.GT.0) THEN
C  CHECK VALID RANGE ON MRSURF
            Y0TEST=Y00+VELY*ZRAD
C           IF (NLTRC) 
C     .       WRITE (iunout,*) 'CHECK VALID RANGE: Y0TEST ',Y0TEST
            IF (Y0TEST.GE.PSURF(1).AND.Y0TEST.LE.PSURF(NP2ND)) THEN
              IPCELL=LEARCA(Y0TEST,PSURF,1,NP2ND,1,'TIMEP 2     ')
            ELSE
              MRSURF=0
              MTSURF=0
              NINCX=0
              NINCZ=0
            ENDIF
          ENDIF
          MPSURF=0
          NINCY=0
C         IF (NLTRC) THEN
C           WRITE (iunout,*) 'NO REENTRY FOUND '
C         ENDIF
          Y00=Y00+ZRD*VELY
          GOTO 5000
        ENDIF
C
C  PARTICLE IN STANDARD MESH, RADIAL CELL NRCELL
C
2900    CONTINUE
        Y000=Y00+VELY*ZRAD
C       IF (NLTRC) WRITE (iunout,*) 'Y000,Y00,ZRAD ',Y000,Y00,ZRAD
C
        DUM=0.
        NCOUP=1
C
C  J1: CELL INDEX
C  J2: SURFACE INDEX
        J1=NPCELL
        IF (VELY.LT.0.) THEN
          INCY=0
          NINCY=-1
          IF (NLSRFY) J1=MPSURF-1
        ELSE
          INCY=1
          NINCY=1
          IF (NLSRFY) J1=MPSURF
        ENDIF
        J2=J1+INCY
C
        NLSRFY=.FALSE.
C
3000    CONTINUE
        IF (J2.LE.0.OR.J2.GT.NP2ND) THEN
          WRITE (iunout,*) 'ERROR IN TIMEP ',J2,J1,VELY
          CALL EXIT_OWN(1)
        ENDIF
C  TIME FROM Y00 TO PSURF
        IF (MPSURF.EQ.J2) THEN
          J1=J1+NINCY
          J2=J1+INCY
          GOTO 3000
        ENDIF
        DY=(PSURF(J2)-Y00)
!pb reduce zrad if the trajectory comes too close to the next standard surface
!pb without hitting it
        IF (ABS(ABS(DY/ZRAD)-1._DP) <= EPS10) THEN
          ZRAD=(1._DP-EPS6)*ZRAD
          ZRADS=ZRAD
        END IF
        BB=VELY+EPS60
        F=DY/BB
C       IF (NLTRC) WRITE (iunout,*) 'J2,F,DY ',J2,F,DY
        IF (F.LE.ZRAD) THEN
          JUPC(NCOUP)=J1
          LUPC(NCOUP)=J2
          ALPD(NCOUP)=F-DUM
          DUM=F
C  STOP HISTORY AT NON DEFAULT STANDARD SURFACE J2
          ITEST=INMP2I(0,J2,0)
          IN=ITEST+NLIM
!pb          IF (ITEST.NE.0.AND.ILIIN(IN).NE.0) THEN
          IF ((ityp==3).or.(ITEST.NE.0.AND.ILIIN(IN).NE.0)) THEN
            ZRAD=F
            ISRFCL=0
            NPCELL=J1
            MPSURF=J2
            MRSURF=0
            NINCX=0
            NINCZ=0
            MTSURF=0
            MASURF=0
            IPCELL=NPCELL
            Y00=PSURF(J2)
            GOTO 5000
          ENDIF
C  NEXT CELL
          J1=J1+NINCY
          J2=J1+INCY
          NCOUP=NCOUP+1
          GOTO 3000
        ENDIF
C
C  LAST CELL
C
3100    CONTINUE
C       IF (NLTRC) WRITE (iunout,*) 'LAST CELL ',ZRAD,DUM
        IF (MPSURF.EQ.0) THEN
          NPCELL=J1
        ELSE
          Y0T=Y00+VELY*ZRAD
          NPCELL=LEARCA(Y0T,PSURF,1,NP2ND,1,'TIMEP 3     ')
        ENDIF
        MPSURF=0
        JUPC(NCOUP)=J1
        ALPD(NCOUP)=ZRAD-DUM
        IPCELL=NPCELL
        Y00=Y000
        GOTO 5000
C
      ENDIF
C
5000  CONTINUE
      DO 5100 J=1,NCOUP
        CLPD(IADD+J)=ALPD(J)
        NUPC(IADD+J)=(JUPC(J)-1)+(NTCELL-1)*NP2T3
        NCOUNP(IADD+J)=JUPC(J)
        NCOUNS(IADD+J)=NTCELL
        IF (ALPD(J).LE.0..OR.JUPC(J).LE.0.OR.JUPC(J).GE.NP2ND) THEN
          WRITE (iunout,*) 'ERROR DETECTED IN TIMEP '
          WRITE (iunout,*) 'J,IADD+J,ALPD,JUPC ',
     .                      J,IADD+J,ALPD(J),JUPC(J)
        ENDIF
C       IF (NLTRC) WRITE (iunout,*) 'TIMEP ',
C     .     J+IADD,CLPD(J+IADD),NUPC(J+IADD),NCOUNP(J+IADD)
5100  CONTINUE
C
      IF (ICOUT.LT.NCOUT.AND.MPSURF.EQ.0) THEN
        ICOUT=ICOUT+1
        IPOLGN=0
        IF (ICOUT.EQ.NCOUT) IPOLGN=IPOLGS
        ZRAD=BLPD(ICOUT)
        if (ncpan.gt.0) then
          jn=0
          do j=ncpan,ncpen
            jn=jn+1
            alpd(jn)=alpd(j)
            jupc(jn)=jupc(j)
            lupc(jn)=lupc(j)
            mupc(jn)=mupc(j)
          end do
        endif
        NTCELL=NCOUNT(ICOUT)
        IADD=IADD+NCOUP
        ncoup=jn
C       IF (NLTRC)
C    .  WRITE (iunout,*) 'NEXT TOR. CELL: ZRAD,NTCELL,IADD ',
C    .                                    ZRAD,NTCELL,IADD
        GOTO 10000
      ENDIF
C
      NCOU=IADD+NCOUP
C
      SUM=0.
      DO 5110 ICOU=1,NCOU
        SUM=SUM+CLPD(ICOU)
        NCOUNT(ICOU)=NCOUNS(ICOU)
C       WRITE (iunout,*) 'ICOU,NCOUNT,CLPD ',
C    .                    ICOU,NCOUNT(ICOU),CLPD(ICOU)
5110  CONTINUE
      IF (MPSURF.EQ.0.AND.ABS(SUM-ZRADS).GT.EPS10) THEN
        WRITE (iunout,*) 'ERROR IN TIMEP: NPANU,SUM,ZRADS ',
     .                               NPANU,SUM,ZRADS
        RETURN ! changed to avoid crash of long runs
!        CALL EXIT_OWN(1)
      ENDIF
C
      ZRAD=SUM
      RETURN
C
991   CONTINUE
      WRITE (iunout,*) 
     .  'REENTRANCE FROM VACUUM REGION AND LEVGEO=1 IN TIMEP'
      CALL EXIT_OWN(1)
      END
