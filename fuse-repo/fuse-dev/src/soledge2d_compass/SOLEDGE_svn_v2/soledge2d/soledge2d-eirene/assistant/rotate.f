 
C
C*DK ROTATE
      SUBROUTINE EIRENE_ROTATE(VLABX,VLABY,VLABZ,VLOCX,VLOCY,VLOCZ,
     .                  CX,CY,CZ,CS)
C
C   VLAB IST DIE RICHTUNG DES EINFLIEGENDEN TEILCHENS IM LABORSYSTEM
C   VLOC IST DIE RICHTUNG DES REFLEKTIERTEN TEILCHENS IM LOKALEN SYSTEM
C   C IST DIE POSITIVE X RICHTUNG IM LOKALEN SYSTEM
C   ES WIRD VLOC INS LABORSYSTEM ZURUECKTRANSFORMIERT UND ALS VLAB
C   ZURUECKGEGEBEN
C   DAS LOKALE SYSTEM WIRD SO BESTIMMT, DASS DAS TEILCHEN IN SEINER
C   X-Z-EBENE EINFLIEGT, MIT POSITIVER Z-GESCHWINDIGKEIT.
C   ES WIRD VORAUSGESETZT, DASS CS = COS(C,VLAB) POSITIV
C   UND DASS VLOCX NEGATIV (D.H. BEIM OUTPUT COS(C,VLAB).LE.0.)
C
      USE EIRMOD_PRECISION
      IMPLICIT NONE
 
      REAL(DP), INTENT(INOUT) :: VLABX, VLABY, VLABZ, VLOCX, VLOCY,
     .                           VLOCZ
      REAL(DP), INTENT(IN) :: CX, CY, CZ, CS
      REAL(DP) :: SS, SSI, A2, A3, B2, B3, C2, C3
 
      SS=SQRT(1.-CS*CS)
      SSI=1./SS
C
C     A1=CX
C     B1=CY
C     C1=CZ
C
      A3=(VLABX-CS*CX)*SSI
      B3=(VLABY-CS*CY)*SSI
      C3=(VLABZ-CS*CZ)*SSI
C   (A2,B2,C2)=(A3,B3,C3) KREUZ (A1,B1,C1)
C
      A2=B3*CZ-C3*CY
      B2=C3*CX-A3*CZ
      C2=A3*CY-B3*CX
C
C     ROTATE WITH MATRIX/A1 A2 A3/
C                       /B1 B2 B3/
C                       /C1 C2 C3/
      VLABX=CX*VLOCX+A2*VLOCY+A3*VLOCZ
      VLABY=CY*VLOCX+B2*VLOCY+B3*VLOCZ
      VLABZ=CZ*VLOCX+C2*VLOCY+C3*VLOCZ
      RETURN
C
      ENTRY EIRENE_ROTATF(VLABX,VLABY,VLABZ,VLOCX,VLOCY,VLOCZ,
     .             CX,CY,CZ)
C  HIER IST ENTWEDER CS=1., D.H. SENKRECHTER EINFLUG, ODER
C  DIE ORIENTIERUNG DER Y-Z-ACHSEN IM LOCALEN SYSTEM SPIELT
C  WEGEN DER SYMMETRIE DER VERTEILUNG DES REFLEXIONSWINKELS
C  KEINE ROLLE, IST INSB. UNABHAENGIG VON VLAB WAEHLBAR
C
C  1. FALL:  ABS(CZ).NE.1.
C
      IF (ABS(CZ).GE.0.99999) GOTO 1
C
      SS=SQRT(CY*CY+CX*CX)
      SSI=1./SS
C
      A2=-CY*SSI
      B2=CX*SSI
C     C2=0.
C
      A3=-CZ*B2
      B3=CZ*A2
C     C3=SS
C     ROTATE WITH MATRIX/A1 A2 A3/
C                       /B1 B2 B3/
C                       /C1 C2 C3/
      VLABX=CX*VLOCX+A2*VLOCY+A3*VLOCZ
      VLABY=CY*VLOCX+B2*VLOCY+B3*VLOCZ
      VLABZ=CZ*VLOCX+         SS*VLOCZ
      RETURN
C
C  2. FALL: CZ=1. ODER CZ=-1., D.H. CX=CY=0.
C
 1    CONTINUE
C     A2=0.
C     B2=-CZ
C     C2=0.
C
C     A3=1.=CZ*CZ
C     B3=0.
C     C3=0.
C     ROTATE WITH MATRIX/A1 A2 A3/
C                       /B1 B2 B3/
C                       /C1 C2 C3/
      VLABX=                     VLOCZ
      VLABY=        -CZ*VLOCY
      VLABZ=CZ*VLOCX
      RETURN
C
      ENTRY EIRENE_ROTATI(VLABX,VLABY,VLABZ,VLOCX,VLOCY,VLOCZ,
     .             CX,CY,CZ)
C  WIE BEI ENTRY ROTATF, ABER ES WIRD MIT INVERSER
C  (=TRANSPONIERTER) MATRIX GEDREHT.
C
C  1. FALL:  ABS(CZ).NE.1.
C
      IF (ABS(CZ).GE.0.99999) GOTO 2
C
      SS=SQRT(CY*CY+CX*CX)
      SSI=1./SS
C
      A2=-CY*SSI
      B2=CX*SSI
C     C2=0.
C
      A3=-CZ*B2
      B3=CZ*A2
C     C3=SS
C     ROTATE WITH MATRIX/A1 A2 A3/
C                       /B1 B2 B3/
C                       /C1 C2 C3/
      VLOCX=CX*VLABX+CY*VLABY+CZ*VLABZ
      VLOCY=A2*VLABX+B2*VLABY
      VLOCZ=A3*VLABX+B3*VLABY+SS*VLABZ
      RETURN
C
C  2. FALL: CZ=1. ODER CZ=-1., D.H. CX=CY=0.
C
 2    CONTINUE
C     A2=0.
C     B2=-CZ
C     C2=0.
C
C     A3=1.=CZ*CZ
C     B3=0.
C     C3=0.
C     ROTATE WITH MATRIX/A1 A2 A3/
C                       /B1 B2 B3/
C                       /C1 C2 C3/
      VLOCX=                  CZ*VLABZ
      VLOCY=        -CZ*VLABY
      VLOCZ=   VLABX
      RETURN
      END
