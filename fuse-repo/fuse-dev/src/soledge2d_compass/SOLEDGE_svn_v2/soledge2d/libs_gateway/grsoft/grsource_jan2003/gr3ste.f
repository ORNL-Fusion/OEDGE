C@PROCESS OPT(3) NOSDUMP NOGOSTMT FIPS(F)
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C UPDATE:  19. 7.91  GROTEN
C Busch 9. 6.93 OPEN,CLOSE macht Probleme auf DEC Workstation
C               fuer GR3STE war der interaktive IO ausserdem falsch,
C               da dieses Program auch im Batch auf der CRAY laufen
C               koennte
c
c     debug subchk
c     end debug
      SUBROUTINE GR3STE(AR,IER,CENTR,ART,IOCULI)
c     .. scalar arguments ..
      REAL              CENTR
      INTEGER           IER,IOCULI
      CHARACTER(len=*)  ART
c     ..
c     .. array arguments ..
      REAL              AR(*)
c     ..
c     .. local scalars ..
      REAL              BILD11,BILD12,BILD13,BILD14,BILD21,BILD22,
     $                  BILD23,BILD24,PP1,PP14,PP2,PP3,PP4,RANDX,RANDY,
     $                  VERS,VERSFA,XMAX,XMIN,YMAX,YMIN
c     ..
c     .. local arrays ..
      REAL              AXT(3,3),EXT(3,3),FEST(2,2),SCHUB(3)
c     ..
c     .. external subroutines ..
      EXTERNAL          GR3CEN,GR3EXT,GR3PLO,GR3TRL,GR3WIN,GRCHRC,
     $                  GRSCLC,GRSCLV
c     ..
c     .. intrinsic functions ..
      INTRINSIC         MAX,MIN
c     ..
c     .. common blocks ..
      COMMON            /GRPIC/FLPIC,NSCLC,NSCLV,NSCLP,RAHMEN,XMAXCM,
     $                  YMAXCM,XDCPIC,YDCPIC
CDEC$ PSECT /GRPIC/ NOSHR
      COMMON            /GRPP/PP
CDEC$ PSECT /GRPP/ NOSHR
      REAL              XDCPIC,XMAXCM,YDCPIC,YMAXCM
      INTEGER           FLPIC,NSCLC,NSCLP,NSCLV,RAHMEN
      REAL              PP(18)
      SAVE /GRPIC/,/GRPP/
c     ..
c     .. data statements ..
      DATA              SCHUB/0.,0.,0./
c     ..
c     .. executable statements ..
c     write(93,*) 'gr3ste- pp'
c     write(93,*) pp(1),pp(2),pp(3),pp(4)
c     write(93,*) pp(5),pp(6),pp(7),pp(8)
*
      IF (IOCULI.EQ.1) THEN
         VERSFA = 0.05
      ELSEIF (IOCULI.EQ.2) THEN
         VERSFA = 0.1
      ELSEIF (IOCULI.EQ.799) THEN
         VERSFA = 0.1
      ELSE
         VERSFA = 0.05
      END IF
*
      PP1    = PP(1)
      PP2    = PP(2)
      PP3    = PP(3)
      PP4    = PP(4)
      PP14   = PP(14)
c---- groesse fuer die 2 bilder
      RANDX  = XMAXCM/300.
      RANDY  = YMAXCM/300.
      BILD11 = RANDX
      BILD12 = RANDY
      BILD13 = XMAXCM/2. - RANDX
      BILD14 = YMAXCM - RANDY

      BILD21 = XMAXCM/2. + RANDX
      BILD22 = RANDY
      BILD23 = XMAXCM - RANDX
      BILD24 = YMAXCM - RANDY

c---- schriftgroesse kleiner machen
      CALL GRCHRC(PP14*0.6,0.,0)
c---- versetzung berechnen, verschieben, zentralprojektion, extrema
      CALL GR3EXT(AR,IER,EXT)
      VERS   = (EXT(1,3)-EXT(1,1))*VERSFA
      SCHUB(1) = VERS
      CALL GR3TRL(AR,IER,SCHUB)
      CALL GR3CEN(AR,IER,CENTR)
      CALL GR3EXT(AR,IER,EXT)
c---- aufheben der zentralprojektion, zur anderen seite schieben,
c---- zentalprojektion, extrema
      CALL GR3CEN(AR,IER,0.)
      SCHUB(1) = -2.*VERS
      CALL GR3TRL(AR,IER,SCHUB)
      CALL GR3CEN(AR,IER,CENTR)
      CALL GR3EXT(AR,IER,AXT)
c---- festes projektionsfeld
      XMIN   = MIN(EXT(1,1),AXT(1,1))
      XMAX   = MAX(EXT(1,3),AXT(1,3))
      YMIN   = MIN(EXT(2,1),AXT(2,1))
      YMAX   = MAX(EXT(2,3),AXT(2,3))
      FEST(1,1) = XMIN - (XMAX-XMIN)*0.01
      FEST(1,2) = XMAX + (XMAX-XMIN)*0.01
      FEST(2,1) = YMIN - (YMAX-YMIN)*0.01
      FEST(2,2) = YMAX + (YMAX-YMIN)*0.01
      CALL GR3WIN(FEST,IER)
c---- ausgabe des zweiten, rechten bildes
c     write(93,*) 'gr3ste - 1'
c     write(93,*)  bild21,bild22,bild23,bild24
c     write(93,*)  bild21,bild22,bild23,bild24
      CALL GRSCLC(BILD21,BILD22,BILD23,BILD24)
      CALL GRSCLV(BILD21,BILD22,BILD23,BILD24)
c---- orientierungshilfe zur justierung der stereobrille
c     call grjmps(20.,20.,4)
C Busch 9. 6.93 OPEN,CLOSE macht Probleme auf DEC Workstation
C Busch 9. 6.93 I/O darf nicht sein, da GR3STE im Batch auf CRAY erlaubt
      IER    = 1
      CALL GR3PLO(AR,IER,ART)
c--------nicht nur test: nicht loeschen !    bremse
cori  IF (IER.GE.10) THEN
cori     READ (5,FMT='(A1)',END=10) STR
c  10    CONTINUE
CORI     CLOSE (5)
CORI     OPEN (5)
cori  END IF
c-----end nicht nur test
c---- aufheben der zentralprojektion, verschieben nach rechts
c---- erneut zentralprojektion, ausgabe linkes bild
      CALL GR3CEN(AR,IER,0.)
      SCHUB(1) = -SCHUB(1)
      CALL GR3TRL(AR,IER,SCHUB)
      CALL GR3CEN(AR,IER,CENTR)
c     write(93,*) 'gr3ste - 2'
c     write(93,*)  bild11,bild12,bild13,bild14
c     write(93,*)  bild11,bild12,bild13,bild14
      CALL GRSCLC(BILD11,BILD12,BILD13,BILD14)
      CALL GRSCLV(BILD11,BILD12,BILD13,BILD14)
c---- orientierungshilfe zur justierung der stereobrille
c     call grjmps(0.25,20.,5)
      IER    = 1
      CALL GR3PLO(AR,IER,ART)
c--------nicht nur test: nicht loeschen !    bremse
C Busch 9. 6.93 OPEN,CLOSE macht Probleme auf DEC Workstation
C Busch 9. 6.93 I/O darf nicht sein, da GR3STE im Batch auf CRAY erlaubt
c     IF (IER.GE.10) THEN
c        READ (5,FMT='(A1)',END=20) STR
c  20    CONTINUE
CORI     CLOSE (5)
CORI     OPEN (5)
c     END IF
c---- aufheben der zentralprojektion, schieben in ursprungslage
c---- normale zentrierung und zeichengroesse
      CALL GR3CEN(AR,IER,0.)
      SCHUB(1) = -VERS
      CALL GR3TRL(AR,IER,SCHUB)
      FEST(1,1) = 0.
      FEST(1,2) = 0.
      FEST(2,1) = 0.
      FEST(2,2) = 0.
      CALL GR3WIN(FEST,IER)
      CALL GRCHRC(PP14,0.,0)

      END
