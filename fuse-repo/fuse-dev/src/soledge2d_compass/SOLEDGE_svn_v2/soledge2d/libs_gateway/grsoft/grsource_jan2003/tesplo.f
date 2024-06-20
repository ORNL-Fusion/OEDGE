      PROGRAM TESPLO

      INTEGER, PARAMETER :: NXD=30,NYD=30
      REAL*8 :: TAB(0:NXD+1,0:NYD+1), TMIN, TMAX
      CHARACTER(40) :: NAM

      REAL ROTX,ROTZ,ROTY,angx,angz,plxa,plxe,plya,plye
      INTEGER IXA,IXE,IYA,IYE,NCONT 
      LOGICAL LCONT,lshad,L3D,LARBI3D,LARBISO,LCUT,LRAPS,lgrid
      COMMON /COMMPL/ ROTX,ROTZ,ROTY,angx,angz,plxa,plxe,plya,plye,
     .  IXA,IXE,IYA,IYE,NCONT,lshad,L3D,LCONT,LARBI3D,LARBISO,LCUT,
     .  LRAPS,lgrid

      TMIN=666.D0
      TMAX=666.D0

      TAB=0.D0
      DO IX=1,20
        DO IY=1,20
          TAB(IX,IY)=IX*IY
        END DO
      END DO

      NAM = 'TAB                          '
      IXA=1
      IXE=20
      IYA=1
      IYE=20
      L3D=.TRUE.
      LCONT=.TRUE.
      NCONT=18
      ANGX=25.
      ANGZ=210.
      ROTX = -60.
      ROTY = 0.
      ROTZ = 30.

      CALL GRSTRT (35,8)
      CALL GRSCLP (39.5,39.5,0)
      CALL PLRECT (TAB,NAM,666.D0,666.D0,.FALSE.,0)
      call grnxtf
      CALL GRSCLP (39.5,39.5,0)
      CALL PLRECT (TAB,NAM,666.D0,666.D0,.FALSE.,0)
      CALL GRENDE

      STOP

      END PROGRAM TESPLO
*//PLRECT//
*=======================================================================
*                 S U B R O U T I N E  P L R E C T
*=======================================================================
*
*
*  1. PURPOSE
*
*     PLRECT PROVIDES GRAPHICAL OUTPUT ON A RECTANGULAR GRID
*
*  2. SPECIFICATIONS
*      TAB: A REAL (NDIMX X NDIMY)-MATRIX, WHICH CONTAINS THE VALUES
*           OF THE FUNCTION THAT SHOULD BE PLOTTED
*      IXA, IXE : INTEGERS, BEGIN AND END VALUES ON THE X-AXIS
*      IYA, IYE : INTEGERS, BEGIN AND END VALUES ON THE Y-AXIS
*      X : REAL ARRAY WITH DIMENSION NDIMX, X-VALUES
*      Y : REAL ARRAY WITH DIMENSION NDIMY, Y-VALUES
*      IOPT : INTEGER TO SPECIFY PLOTTYPE
*         IOPT=1 3-DIMENSIONAL PLOT
*         IOPT=2 CONTOUR PLOT
*         IOPT=3 BOTH PLOTS ON ONE PAGE
*      NCONT : INTEGER, NUMBER OF CONTOUR LINES
*      NAME : CHARACTER*40, HEADER OF THE PLOT
*
*=======================================================================
*//DECLARATIONS//

      SUBROUTINE PLRECT(TABIN,NAME,MINI,MAXI,LLOG,IS)
      IMPLICIT NONE

      INTEGER, PARAMETER :: NXD=30,NYD=30

      REAL ROTX,ROTZ,ROTY,angx,angz,plxa,plxe,plya,plye
      INTEGER IXA,IXE,IYA,IYE,NCONT 
      LOGICAL LCONT,lshad,L3D,LARBI3D,LARBISO,LCUT,LRAPS,lgrid
      COMMON /COMMPL/ ROTX,ROTZ,ROTY,angx,angz,plxa,plxe,plya,plye,
     .  IXA,IXE,IYA,IYE,NCONT,lshad,L3D,LCONT,LARBI3D,LARBISO,LCUT,
     .  LRAPS,lgrid

* -- INPUT AND LOCAL VARIABLES

      INTEGER MAXCONT
      PARAMETER (MAXCONT=50)

      INTEGER I,J,IS
      INTEGER ICOL(MAXCONT)
      LOGICAL LLOG

      REAL DRDMPA(15),XMIN,XMAX,YMIN,YMAX,TABMIN,TABMAX, DUMMY
      REAL X(NXD),Y(NYD),TAB(NXD,NYD),WERT(MAXCONT)
      INTEGER ILEN,IDUMMY,ILENGTH
      REAL*8 TABIN(0:NXD+1,0:NYD+1),MAXI,MINI
      REAL*8 HUGEREAL

      CHARACTER*8  OPT
      CHARACTER*12 DIST
      CHARACTER*40 NAME
      CHARACTER*4 CCONT
      CHARACTER*53 AUSNAME
      CHARACTER*10 FORM
      DATA OPT/'HREIK'/

      HUGEREAL = 10.**RANGE(1.E0)
      WRITE(6,*) HUGEREAL,HUGE(1.)
      DO I=1,NXD
        X(I) = REAL(I)
        DO J=1,NYD
          TAB(I,J)=MAX(MIN(TABIN(I,J),HUGEREAL),-HUGEREAL)
        ENDDO
      ENDDO
      DO J=1,NYD
        Y(J) = REAL(J)
      ENDDO
      IF (IS .GT. 0) THEN
        ILEN = ILENGTH(NAME)
        FORM = '(A  ,A,I2)'
        WRITE(FORM(3:4),'(I2)') ILEN
        WRITE(6,*) 'FORM:',FORM
        WRITE(AUSNAME,FORM) NAME,' for fluid ',IS
      ELSE
        AUSNAME = NAME
      ENDIF
C
C     SEARCHING FOR THE MAX AND THE MIN OF THE MATRIX TAB
C
      IF (ABS(MINI-666.) .GT. 1.E-6) THEN
        DO I=IXA,IXE
          DO J=IYA,IYE
            IF (TAB(I,J) .LT. REAL(MINI)) TAB(I,J)=REAL(MINI)
          ENDDO
        ENDDO
      ENDIF
      IF (ABS(MAXI-666.) .GT. 1.E-6) THEN
        DO I=IXA,IXE
          DO J=IYA,IYE
            IF (TAB(I,J) .GT. REAL(MAXI)) TAB(I,J)=REAL(MAXI)
          ENDDO
        ENDDO
      ENDIF

      TABMIN=MINVAL(TAB(IXA:IXE,IYA:IYE))
      TABMAX=MAXVAL(TAB(IXA:IXE,IYA:IYE))

      IF (LLOG) THEN
        TABMIN=LOG10(MAX(TABMIN,1.E-30))
        TABMAX=LOG10(MAX(TABMAX,1.E-30))
        DO I=IXA,IXE
          DO J=IYA,IYE
            TAB(I,J)=LOG10(MAX(TAB(I,J),1.E-30))
          ENDDO
        ENDDO
      ENDIF

      IF (L3D) THEN
        XMIN=MINVAL(X(IXA:IXE))
        XMAX=MAXVAL(X(IXA:IXE))
        YMIN=MINVAL(Y(IYA:IYE))
        YMAX=MAXVAL(Y(IYA:IYE))
        DRDMPA(1)=0.0546875
        DRDMPA(2)=ANGZ
        DRDMPA(3)=ANGX
        DRDMPA(4)=27.
        DRDMPA(5)=0.
        DRDMPA(6)=XMIN
        IF (MOD(XMAX-XMIN,10.) .EQ. 0.) THEN
          DRDMPA(7)=XMAX
        ELSE
          DRDMPA(7)=XMIN+(INT((XMAX-XMIN)/10.)+1)*10.
        ENDIF
        DRDMPA(8)=YMIN
        IF (MOD(YMAX-YMIN,5.) .EQ. 0.) THEN
          DRDMPA(9)=YMAX
        ELSE
          DRDMPA(9)=YMIN+(INT((YMAX-YMIN)/5.)+1)*5.
        ENDIF
        DRDMPA(10)=TABMIN
        DRDMPA(11)=TABMAX
        DRDMPA(12)=IXA
        DRDMPA(13)=IXE
        DRDMPA(14)=IYA
        DRDMPA(15)=IYE
        CALL GRFRBN(2,1,1,2,1)
      ENDIF

      IF (LCONT) THEN
        IF (NCONT.EQ.0) NCONT=18
        WRITE (CCONT,'(I4)') NCONT
C       CALCULATE THE VALUES TO PLOT THE CONTOUR LINES
        DO I=1,NCONT
          ICOL(I)=1+(I-1)/3
          WERT(I)=TABMIN+(I-1)*(TABMAX-TABMIN)/(NCONT-1)
        ENDDO
        WRITE (DIST,'(1P,E12.4)') ABS(WERT(2)-WERT(1))
      ENDIF

      IF (L3D .AND. LCONT) THEN
        CALL GR90DG
        CALL GRCHRC(.4,0.,16)
C
C       PLOT THE FUNCTION AND THE CONTOUR PLOTS
C
        DRDMPA(4)=18
        CALL GRSCLC(5.,22.,16.,38.)
        CALL GRDRDM(DRDMPA,NXD,TAB,X,Y)
        DRDMPA(1)=0.
        CALL GRDSH(0.1,0.2,0.1)
        CALL GRDRDU(DRDMPA,NXD,TAB,X,Y)
        CALL GRDSH(1.,0.,1.)
C
C       PLOT THE CONTOUR PLOTS
C
        CALL GRSCLC(3.,5.,18.,20.)
        CALL GRHHNL(NXD,TAB,IXE-IXA+1,X,IYE-IYA+1,Y,NCONT,WERT,ICOL,
     .              OPT,IXA,1,IYA,1)
C
C       PLOT THE DISCRIPTION OF THE PICTURE
C
        CALL GRCHRC(.4,0.,16)
        CALL GRTXT(2.,-1.,-1,AUSNAME)
        CALL GRTXT(2.,-2.,-1,'Number of contour lines:')
        CALL GRTXTC(-1,CCONT)
        CALL GRTXT(2.,-3.,-1,'Difference between the contour lines:')
        CALL GRTXTC(-1,DIST)
        CALL GRSCLV(REAL(IXA),REAL(IYA),REAL(IXE),REAL(IYE))
        CALL GRAXS(-1,'X=2,Y=3',-1,' ',-1,' ')
        CALL GRNXTF
      ELSEIF (L3D) THEN
C
C       PLOT THE FUNCTION
C
        CALL GRSCLC(6.1,3.,22.1,28.7)
        CALL GRCHRC(.4,0.,16)
        CALL GRDRDM(DRDMPA,NXD,TAB,X,Y)
        DRDMPA(1)=0.
        CALL GRDSH(0.1,0.2,0.1)
        CALL GRDRDU(DRDMPA,NXD,TAB,X,Y)
        CALL GRDSH(1.,0.,1.)
        CALL GRTXT(9.,-1.,-1,AUSNAME)
        CALL GRSCLC(0.,0.,28.7,28.7)
        CALL GRSCLV(0.,0.,28.7,28.7)
        CALL GRCHRC(0.3,35.,IDUMMY)
C       CALL GRTXT(11.,24.5,-1,'INNER TARGET')
C       CALL GRTXT(25.7,20.5,-1,'OUTER TARGET')
        CALL GRCHRC(0.3,-15.,IDUMMY)
C       CALL GRTXT(23.,27.,-1,'WALL')
C       CALL GRTXT(8.,20.,-1,'PF')
C       CALL GRTXT(15.,20.2,-1,'CORE')
C       CALL GRTXT(13.,18.,-1,'PF')
        CALL GRCHRC(0.3,0.,IDUMMY)
        CALL GRNXTF
      ELSEIF (LCONT) THEN
        CALL GRSCLC(3.,5.,24.,26.)
        CALL GRCHRC(.4,0.,16)
        CALL GRHHNL(NXD,TAB,IXE,X,IYE,Y,NCONT,WERT,ICOL,
     .              OPT,IXA,1,IYA,1)
C
C       PLOT THE DESCRIPTION OF THE GRAPHIC
C
        CALL GRSCLV(3.,5.,24.,26.)
        CALL GRTXT(4.,3.,-1,AUSNAME)
        CALL GRTXT(4.,2.,-1,'Number of contour lines:')
        CALL GRTXTC(-1,CCONT)
        CALL GRTXT(4.,1.,-1,'Difference between the contour lines:')
        CALL GRTXTC(-1,DIST)
        CALL GRSCLV(REAL(IXA),REAL(IYA),REAL(IXE),REAL(IYE))
        CALL GRAXS(-1,'X=2,Y=3',-1,' ',-1,' ')
        CALL GRNXTF
      ENDIF
      END

      INTEGER FUNCTION ILENGTH(TXT)
      IMPLICIT NONE
      CHARACTER*(*) TXT
      INTEGER I

      I = LEN(TXT)
10    IF (TXT(I:I) .EQ. ' ') THEN
        I=I-1
        IF (I .GT. 0) GOTO 10
      ENDIF
      ILENGTH=I
      END
