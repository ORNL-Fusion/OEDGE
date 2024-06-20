C@PROCESS OPT(3) NOSDUMP NOGOSTMT
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C-----------------------------------------------------------------------
C     GRDN
C     GR-SOFTWARE : DIN WERTE WERDEN ZUGEWIESEN
C     AUTHOR: Marlene Busch
C     Date:   14. 12. 1990
C     UPDATE:    8.  6. 1991
c
C     INPUT:
C     INTEGER IDIN    0,1,2,3,4
C     Output:
C     REAL X,Y
C     Default: Rueckgabe A0 Format, falls IDIN nicht gueltigen
C     Wert hatte
C     Es wird entgegen der Norm immer ein liegendes Rechteck in cm
C     zurueckgegeben. (entsprechend dem Laenge/Breite Verhaeltnis
C     von Bildschirm bzw. A0 Format des Plotters)
C     Die Werte koennen dann zum Aufruf von GRSCLP benutzt werden.
C
C-----------------------------------------------------------------------
      SUBROUTINE GRDN  (IDIN, X, Y)

      INTEGER IDIN
      REAL    X,Y

      IF ( IDIN.LT.0. OR.IDIN.GT.4) IDIN=0
      IF (IDIN.EQ.0) THEN
c        Y=83.1
         Y=83.095
c        X=117.9
         X=117.895
         ELSEIF ( IDIN .EQ. 1) THEN
            Y=58.4
            X=83.1
            ELSEIF ( IDIN .EQ. 2) THEN
               X=58.4
               Y=41.0
               ELSEIF ( IDIN .EQ. 3) THEN
                  X=41.0
                  Y=28.7
                  ELSEIF ( IDIN .EQ. 4) THEN
                     X=28.7
                     Y=20.0
      ENDIF
      END
