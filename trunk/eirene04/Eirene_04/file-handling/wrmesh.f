c  nov 03:  use relative distances to find neighbor segment,
c           otherwise sometimes problems with non-closing polygons encountered.
      SUBROUTINE WRMESH
      USE PRECISION
      USE PARMMOD
      USE CADGEO
      USE CPLOT
      USE CPOLYG
      USE CGEOM
      USE CLGIN
      USE CGRPTL
      USE CTRIG
      USE CGRID
      IMPLICIT NONE


      INTEGER, PARAMETER :: MAXPOIN=500
      REAL(DP) :: partcont(maxpoin,2,2), maxlen
      REAL(DP) :: XPE, YPE, HELP, XT, YT, PHI1, X1, X2, Y1, Y2, PHI2
      REAL(DP) :: DISTQI,DISTQJ1,DISTQJ2
      INTEGER  :: MAXCONT, ICONT, IPOIN, IWST, IWEN, IWL, IWP,
     .            IWAN, IMN, I, NCONT, J, IUHR, ISTORE, IP, IH, IFOUND,
     .            ICO, IPO, IN, IS, IS1, ITRI
      INTEGER  :: IDIAG(MAXPOIN),irip(maxpoin,2)
      REAL(SP) :: xmin,xmax,ymin,ymax,deltax,deltay,delta,xcm,ycm
      REAL(SP) :: XP,YP

C INITIALISIERUNG DER PLOTDATEN
      xmin = CH2X0-CH2MX
      ymin = CH2Y0-CH2MY
      xmax = CH2X0+CH2MX
      ymax = CH2Y0+CH2MY
      deltax = abs(xmax-xmin)
      deltay = abs(ymax-ymin)
      delta = max(deltax,deltay)
      xcm = 24. * deltax/delta
      ycm = 24. * deltay/delta

C ANZAHL DER KONTOUREN BESTIMMEN
C ILPLG WIRD IM INPUT-BLOCK 3 EINGELESEN
      CALL LEER(2)
      WRITE (6,*) 'SUBROUTINE WRMESH CALLED '
      CALL LEER(1)
      NCONT = 0
      DO I=1,NLIMI
        NCONT = MAX(NCONT,ABS(ILPLG(I)))
      ENDDO
      DO I=NLIM+1,NLIM+NSTSI
        NCONT = MAX(NCONT,ABS(ILPLG(I)))
      ENDDO

      WRITE (6,*) 'NUMBER OF CONTOURS FOR FEM MESH: ',NCONT
      CALL LEER(1)
      if (ncont == 0) return

      ALLOCATE (NCONPOINT(NCONT))
      ALLOCATE (XCONTOUR(MAXPOIN,NCONT)) 
      ALLOCATE (YCONTOUR(MAXPOIN,NCONT)) 
      NCONPOINT = 0
      ICO = 0

      call grnxtf
      call grsclc(3.,3.,3.+real(xcm,kind(1.e0)),3.+real(ycm,kind(1.e0)))
      call grsclv(real(xmin,kind(1.e0)),real(ymin,kind(1.e0)),
     .            real(xmax,kind(1.e0)),real(ymax,kind(1.e0)))

      DO ICONT = 1,NCONT
        IPOIN = 0
        MAXLEN = 0.
        irip=0
C AKTUELLE KONTOUR BESTIMMEN, STUECKE MIT ILPLG=ICONT GEHOEREN ZUR
C AKTUELLEN CONTOUR, ANFANGS UND ENDPUNKT DIESES STUECKES WERDEN AUF
C PARTCONT GESPEICHERT
        DO I=1,NLIMI
          IF (ABS(ILPLG(I)) .EQ. ICONT) THEN
            IUHR=ILPLG(I)
C 0 < RLB(I) < 2
C 2-PUNKT OPTION WIRD IM TIMEA0 AUF RLB=1 ZURUECKGEFUEHRT
            IF ((RLB(I) .GT. 0.) .AND. (RLB(I) .LT. 2.) .AND.
     >          (P3(1,I) .EQ. 1.D55 .OR. P3(2,I) .EQ. 1.D55
     >          .OR. P3(3,I) .EQ. 1.D55)) THEN
              IPOIN = IPOIN + 1
              IF (A3LM(I) .EQ. 0.) THEN
C               X,Y-KOORDINATEN
                PARTCONT(IPOIN,1,1) = P1(1,I)
                PARTCONT(IPOIN,1,2) = P1(2,I)
                PARTCONT(IPOIN,2,1) = P2(1,I)
                PARTCONT(IPOIN,2,2) = P2(2,I)
                idiag(ipoin)=i
              ELSEIF (A2LM(I) .EQ. 0.) THEN
C               X,Z-KOORDINATEN
                PARTCONT(IPOIN,1,1) = P1(1,I)
                PARTCONT(IPOIN,1,2) = P1(3,I)
                PARTCONT(IPOIN,2,1) = P2(1,I)
                PARTCONT(IPOIN,2,2) = P2(3,I)
                idiag(ipoin)=i
              ELSEIF (A1LM(I) .EQ. 0.) THEN
C               Y,Z-KOORDINATEN
                PARTCONT(IPOIN,1,1) = P1(2,I)
                PARTCONT(IPOIN,1,2) = P1(3,I)
                PARTCONT(IPOIN,2,1) = P2(2,I)
                PARTCONT(IPOIN,2,2) = P2(3,I)
                idiag(ipoin)=i
              ENDIF
              maxlen = maxlen +
     >               sqrt((partcont(ipoin,1,1)-partcont(ipoin,2,1))**2
     >                   +(partcont(ipoin,1,2)-partcont(ipoin,2,2))**2)
            ELSE
C  ERROR
              WRITE(6,*) 'FALSCHE ANGABE FUER RLB, RLB = ',RLB(I),
     >                    ILPLG(I),I
            ENDIF
          ENDIF
        ENDDO

        IF (LEVGEO == 3) THEN
        DO I=1,NSTSI
          IF (ABS(ILPLG(NLIM+I)) .EQ. ICONT) THEN
            IUHR=ILPLG(NLIM+I)
            IF (INUMP(I,2) .NE. 0) THEN
C             POLOIDAL
              DO J=IRPTA(I,1),IRPTE(I,1)-1
                IF ((XPOL(J,INUMP(I,2)) .NE. XPOL(J+1,INUMP(I,2))) .OR.
     >              (YPOL(J,INUMP(I,2)) .NE. YPOL(J+1,INUMP(I,2)))) THEN
                  IPOIN = IPOIN + 1
                  PARTCONT(IPOIN,1,1) = XPOL(J,INUMP(I,2))
                  PARTCONT(IPOIN,1,2) = YPOL(J,INUMP(I,2))
                  PARTCONT(IPOIN,2,1) = XPOL(J+1,INUMP(I,2))
                  PARTCONT(IPOIN,2,2) = YPOL(J+1,INUMP(I,2))
                idiag(ipoin)=-i
                irip(ipoin,1)=j
                irip(ipoin,2)=INUMP(I,2)
              maxlen = maxlen +
     >               sqrt((partcont(ipoin,1,1)-partcont(ipoin,2,1))**2
     >                   +(partcont(ipoin,1,2)-partcont(ipoin,2,2))**2)
                ENDIF
              ENDDO
            ELSEIF (INUMP(I,1) .NE. 0) THEN
C             RADIAL
              DO J=IRPTA(I,2),IRPTE(I,2)-1
                IF ((XPOL(INUMP(I,1),J) .NE. XPOL(INUMP(I,1),J+1)) .OR.
     >              (YPOL(INUMP(I,1),J) .NE. YPOL(INUMP(I,1),J+1))) THEN
                  IPOIN = IPOIN + 1
                  PARTCONT(IPOIN,1,1) = XPOL(INUMP(I,1),J)
                  PARTCONT(IPOIN,1,2) = YPOL(INUMP(I,1),J)
                  PARTCONT(IPOIN,2,1) = XPOL(INUMP(I,1),J+1)
                  PARTCONT(IPOIN,2,2) = YPOL(INUMP(I,1),J+1)
                  idiag(ipoin)=-i
                  irip(ipoin,1)=INUMP(I,1)
                  irip(ipoin,2)=j
                  maxlen = maxlen +
     >               sqrt((partcont(ipoin,1,1)-partcont(ipoin,2,1))**2
     >                   +(partcont(ipoin,1,2)-partcont(ipoin,2,2))**2)
                ENDIF
              ENDDO
            ELSE
C  ERROR
              WRITE(6,*) 'CASE NOT FORESEEN: INUMP: ',
     >                     (INUMP(I,J),J=1,3)
            ENDIF
          ENDIF
        ENDDO
        ELSEIF (LEVGEO == 4) THEN
          DO ITRI = 1, NTRII
            DO IS = 1, 3
              IN=INMTI(IS,ITRI)
              IF (IN /= 0) THEN
                IF (ABS(ILPLG(IN)) == ICONT) THEN
                  IUHR=ILPLG(IN)
                  IS1 = IS+1
                  IF (IS1 > 3) IS1=1
                  IPOIN = IPOIN + 1
                  PARTCONT(IPOIN,1,1) = XTRIAN(NECKE(IS,ITRI))
                  PARTCONT(IPOIN,1,2) = YTRIAN(NECKE(IS,ITRI))
                  PARTCONT(IPOIN,2,1) = XTRIAN(NECKE(IS1,ITRI))
                  PARTCONT(IPOIN,2,2) = YTRIAN(NECKE(IS1,ITRI))
                  idiag(ipoin)=IN
                  irip(ipoin,1)=itri
                  irip(ipoin,2)=is
                  maxlen = maxlen +
     >               sqrt((partcont(ipoin,1,1)-partcont(ipoin,2,1))**2
     >                   +(partcont(ipoin,1,2)-partcont(ipoin,2,2))**2)
                  
                END IF
              END IF
            END DO
          END DO
        END IF

        IF (IPOIN.LE.0) THEN
          WRITE(6,*) 'CONTOUR ',ICONT,' NOT FOUND'
          GOTO 1000
        ENDIF

        call grnwpn(icont)
C STUECKE DER AKTUELLEN KONTOUR WERDEN SORTIERT
        DO I=1,IPOIN-1
          XPE = PARTCONT(I,2,1)
          YPE = PARTCONT(I,2,2)
          DISTQI=(PARTCONT(I,2,1)-PARTCONT(I,1,1))**2+
     .           (PARTCONT(I,2,2)-PARTCONT(I,1,2))**2
          IFOUND=0
          DO J=I+1,IPOIN
            DISTQJ1=(XPE-PARTCONT(J,1,1))**2+
     .              (YPE-PARTCONT(J,1,2))**2
            DISTQJ2=(XPE-PARTCONT(J,2,1))**2+
     .              (YPE-PARTCONT(J,2,2))**2
CDR         IF ((XPE .EQ. PARTCONT(J,1,1)) .AND.
CDR  >          (YPE .EQ. PARTCONT(J,1,2))) THEN
            IF (DISTQJ1/DISTQI.LE.1.D-10) THEN
              IFOUND=1
              HELP = PARTCONT(I+1,1,1)
              PARTCONT(I+1,1,1) = PARTCONT(J,1,1)
              PARTCONT(J,1,1) = HELP
              HELP = PARTCONT(I+1,1,2)
              PARTCONT(I+1,1,2) = PARTCONT(J,1,2)
              PARTCONT(J,1,2) = HELP

              HELP = PARTCONT(I+1,2,1)
              PARTCONT(I+1,2,1) = PARTCONT(J,2,1)
              PARTCONT(J,2,1) = HELP
              HELP = PARTCONT(I+1,2,2)
              PARTCONT(I+1,2,2) = PARTCONT(J,2,2)
              PARTCONT(J,2,2) = HELP

              ih=idiag(i+1)
              idiag(i+1)=idiag(j)
              idiag(j)=ih

              ih=irip(i+1,1)
              irip(i+1,1)=irip(j,1)
              irip(j,1)=ih
              ih=irip(i+1,2)
              irip(i+1,2)=irip(j,2)
              irip(j,2)=ih
            ELSEIF (DISTQJ2/DISTQI.LE.1.D-10) THEN
CDR         ELSEIF ((XPE .EQ. PARTCONT(J,2,1)) .AND.
CDR  >              (YPE .EQ. PARTCONT(J,2,2))) THEN
              IFOUND=1
              HELP = PARTCONT(J,1,1)
              PARTCONT(J,1,1) = PARTCONT(J,2,1)
              PARTCONT(J,2,1) = HELP
              HELP = PARTCONT(J,1,2)
              PARTCONT(J,1,2) = PARTCONT(J,2,2)
              PARTCONT(J,2,2) = HELP

              HELP = PARTCONT(I+1,1,1)
              PARTCONT(I+1,1,1) = PARTCONT(J,1,1)
              PARTCONT(J,1,1) = HELP
              HELP = PARTCONT(I+1,1,2)
              PARTCONT(I+1,1,2) = PARTCONT(J,1,2)
              PARTCONT(J,1,2) = HELP

              HELP = PARTCONT(I+1,2,1)
              PARTCONT(I+1,2,1) = PARTCONT(J,2,1)
              PARTCONT(J,2,1) = HELP
              HELP = PARTCONT(I+1,2,2)
              PARTCONT(I+1,2,2) = PARTCONT(J,2,2)
              PARTCONT(J,2,2) = HELP

              ih=idiag(i+1)
              idiag(i+1)=idiag(j)
              idiag(j)=ih

              ih=irip(i+1,1)
              irip(i+1,1)=irip(j,1)
              irip(j,1)=ih
              ih=irip(i+1,2)
              irip(i+1,2)=irip(j,2)
              irip(j,2)=ih
            ENDIF
          ENDDO
          IF (IFOUND.EQ.0) THEN
            WRITE (6,*) 'NO MATCHING POINT FOUND FOR CONTOUR ',ICONT
            write(6,*) i,idiag(i),irip(i,1),irip(i,2),
     >                   partcont(i,1,1),partcont(i,1,2),
     >                   partcont(i,2,1),partcont(i,2,2)
            WRITE (6,*) 'USE NEXT POINT '
            IP=I+1
            write(6,*) iP,idiag(iP),irip(ip,1),irip(ip,2),
     >                    partcont(iP,1,1),partcont(iP,1,2),
     >                    partcont(iP,2,1),partcont(iP,2,2)
          ENDIF
        ENDDO
        IF ((PARTCONT(1,1,1) .NE. PARTCONT(IPOIN,2,1)) .OR.
     >      (PARTCONT(1,1,2) .NE. PARTCONT(IPOIN,2,2))) THEN
          WRITE(6,*) 'CONTOUR ',ICONT,' IS NOT CLOSED'
        ELSE
          WRITE(6,*) 'CLOSED CONTOUR ',ICONT
        ENDIF
        do i=1,ipoin
          write(6,*) i,idiag(i),irip(i,1),irip(i,2),
     >                 partcont(i,1,1),partcont(i,1,2),
     >                 partcont(i,2,1),partcont(i,2,2)
        enddo
        XP = PARTCONT(1,1,1)
        YP = PARTCONT(1,1,2)
        call grjmp(REAL(XP,KIND(1.E0)),REAL(YP,KIND(1.E0)))
        DO I=2,IPOIN
          XP = PARTCONT(I,1,1)
          YP = PARTCONT(I,1,2)
          call grdrw(REAL(XP,KIND(1.E0)),REAL(YP,KIND(1.E0)))
        ENDDO
        XP = PARTCONT(1,1,1)
        YP = PARTCONT(1,1,2)
        call grDRW(REAL(XP,KIND(1.E0)),REAL(YP,KIND(1.E0)))

C  BERECHNUNG VON DELTA ALS MITTLERE LAENGE DER TEILSTUECKE
C  DELTA IST MASS FUER DIE GROESSE DER DREIECKE
        IMN=0
        IF (ICONT .EQ. 1) THEN
          maxlen = maxlen / REAL(IPOIN,KIND(1.D0))
          WRITE(78,*) maxlen
          WRITE(78,*)
        endif

C BESTIMMUNG DES UHRZEIGERSINNS DER KONTOUR
C KONTOUR MUSS FUER DIE TRIANGULIERUNG FOLGENDERMASSEN AUSGEGEBEN
C WERDEN:
C  - IM UHRZEIGERSINN FUER INNERE BEGRENZUNGEN DES GEBIETES (POSITIV)
C  - GEGEN UHRZEIGERSINN FUER AEUSSERE BEGRENZUNGEN DES GEBIETES (NEGATIV)
        YMIN=PARTCONT(1,1,2)
        DO I=1,IPOIN
          IF (PARTCONT(I,2,2) .LT. YMIN) THEN
            YMIN = PARTCONT(I,2,2)
            IMN=I
          ENDIF
        ENDDO
 
        IF (IMN .EQ. 0) THEN
          XT = PARTCONT(1,1,1)
          YT = PARTCONT(1,1,2)
C  PUNKT, DER IM UMLAUF DER VORHERGEHENDE IST
          X1 = PARTCONT(IPOIN,1,1)
          Y1 = PARTCONT(IPOIN,1,2)
C  PUNKT, DER IM UMLAUF DER NAECHSTE IST
          X2 = PARTCONT(1,2,1)
          Y2 = PARTCONT(1,2,2)
        ELSE
C  SONDERFALL IMN=IPOIN ENTFAELLT, DA ERSTER PUNKT GLEICH LETZTER
C  PUNKT GILT
          XT = PARTCONT(IMN,2,1)
          YT = PARTCONT(IMN,2,2)
C  PUNKT, DER IM UMLAUF DER VORHERGEHENDE IST
          X1 = PARTCONT(IMN,1,1)
          Y1 = PARTCONT(IMN,1,2)
C  PUNKT, DER IM UMLAUF DER NAECHSTE IST
          X2 = PARTCONT(IMN+1,2,1)
          Y2 = PARTCONT(IMN+1,2,2)

        ENDIF
 
C  BESTIMME POLARWINKEL VON (X1,Y1) UND (X2,Y2) MIT (XT,YT) ALS URSPRUNG
        PHI1 = ATAN2 (Y1-YT,X1-XT)
        PHI2 = ATAN2 (Y2-YT,X2-XT)

        IF (PHI2 .GT. PHI1) THEN
C  ABSPEICHERUNG ERFOLGTE IM UHRZEIGERSINN
          ISTORE = 1
        ELSE
          ISTORE = -1
        ENDIF
C  IUHR=ILPLG > 0 ==> IM UHRZEIGERSINN AUSGEBEN
C  IUHR=ILPLG < 0 ==> ENTGEGEN DEM UHRZEIGERSINN AUSGEBEN
        IWAN=1
        IWEN=IPOIN
        IWST=1
        IWP=1
        IWL=2
        IF (ISTORE*IUHR .LT. 0.) THEN
          IWAN=IPOIN
          IWEN=1
          IWST=-1
          IWP=2
          IWL=1
        ENDIF

        WRITE(78,*) IPOIN+1
        IF (IUHR > 0) THEN
          ICO=ICO+1
          NCONPOINT(ICO)=IPOIN+1
          IPO=0
        END IF
        DO I=iwan,iwen,iwst
          WRITE(78,'(1P,2(2X,E21.14))')
     >          PARTCONT(I,IWP,1),PARTCONT(I,IWP,2)
          IF (IUHR > 0) THEN
            IPO=IPO+1
            XCONTOUR(IPO,ICO) = PARTCONT(I,IWP,1)
            YCONTOUR(IPO,ICO) = PARTCONT(I,IWP,2)
          END IF
        ENDDO
        WRITE(78,'(1P,2(2X,E21.14))') PARTCONT(IWEN,IWL,1),
     >                               PARTCONT(IWEN,IWL,2)
        IF (IUHR > 0) THEN
          IPO=IPO+1
          XCONTOUR(IPO,ICO) = PARTCONT(IWEN,IWL,1)
          YCONTOUR(IPO,ICO) = PARTCONT(IWEN,IWL,2)
        END IF
1000    CONTINUE
      ENDDO
      NCONTOUR=ICO
C  END OF NCONT LOOP
      call leer(1)
      write (6,*) 'input file fort.78 for FEM mesh generator written '
      call leer(2)
      call grnwpn(1)
      call grnxtf
      END
