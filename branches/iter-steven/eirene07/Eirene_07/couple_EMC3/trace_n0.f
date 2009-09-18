CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     NAME: TIMUSR                                                C
C FUNCTION: follow a particle untill it leaves the cell IRGEN     C
C     Y. FENG                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TIMUSR(IRGEN,X0,Y0,Z0,VELX,VELY,VELZ,NEW,
     .                  ICNXT,TIMET,ICOS,IER,NPANU,SURF_P)
      USE GEOMETRY_PL 
      USE PHYSICAL_CELL
      USE NEUTRAL_TRANSPORT

      IMPLICIT NONE 

      INTEGER ::  IRGEN,NEW,ICNXT,ICOS,IER,NPANU
      REAL*8  ::  X0,Y0,Z0,VELX,VELY,VELZ,TIMET
      LOGICAL ::  SURF_P
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 1.   description of the variables                                  C
C Input:                                                             C
C    IRGEN: initial cell number                                      C
C X0,Y0,Z0: initial position                                         C
C VELX,Y,Z: velocity vector                                          C
C      NEW: = 0, new particle                                        C
C          <=>0, following particle, initially on a surface          C
C   SURF_P: .TRUE. a surface particle                                C
C    NPANU: number of the MC particle                                C
C Output:                                                            C
C    TIMET: time for the particle until it leaves this cell          C
C           (also the travelling distance upto now because VL=1)     C
C   ICNXT : next cell number                                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER,SAVE :: IERR,IRNXT,IPNXT,ITNXT
      INTEGER      :: ISR1,ISTS,I,K
      REAL*8       :: XNORM,YNORM,ZNORM
      REAL*8       :: SCOS,VL
      LOGICAL,SAVE :: DIAGNO=.FALSE.
c     LOGICAL,SAVE :: DIAGNO=.TRUE. 
      LOGICAL      :: OUT_OF_RANGE
C--------------------------------------------------------------------
C begin
      IER = 0
      IERR = 0
c     DIAGNO = NPANU== 139
c     if(NPANU==2) stop

      IF( OUT_OF_RANGE(IC0_N0,1,NCELL_N2) ) THEN
        IF(DIAGNO) CALL WRMESS('IC0.LE.0 .OR. IC0.GT.NCELL_N2')
        IERR = 1
        TIMET=-1.
        WRITE(6,*)'*** ERROR IC0_N0 OUT OF RANGE'             
        WRITE(6,'(4I5)')NZ0_N0,JR0_N0,JP0_N0,JT0_N0
        RETURN
      ENDIF
C New particle
      IF(NEW.EQ.0) THEN
        IRGEN = IC0_N0
        TIMET = 0.

C replace the initial coordinates and velocity
        XP0 = X0
        YP0 = Y0
        ZP0 = Z0

        VX = VELX
        VY = VELY
        VZ = VELZ

C the new particle located just inside a cell
        IF(.NOT.SURF_P) THEN
          ISRF = 0
          ITRIA=0
        ENDIF
C maximal travellin distance to a plate
        TLMAX = 1.D30
        IF(NTRIANG.GT.0) CALL TLMAX_TRAVEL

        XFIN = X0
        YFIN = Y0
        ZFIN = Z0

        IF(DIAGNO) THEN
         WRITE(6,*)'==>TIMUSR FOLLOWS A NEW PARTICLE'
         WRITE(6,*)'==>NPANU:',NPANU                   
         WRITE(6,*)'ISRF  NZ0   JR0  JP0  JT0   IRGEN'
         WRITE(6,'(5I5,I8)')ISRF,NZ0_N0,JR0_N0,JP0_N0,JT0_N0,IRGEN
         WRITE(6,*)'INITIAL COORDINATS FROM EIRENE P0 AND V0:'
         WRITE(6,'(1P,6E12.4)')X0,Y0,Z0,VELX,VELY,VELZ
         WRITE(6,*)'Maximal distance:',TLMAX
         IF(ITRIA.NE.0)  THEN 
           WRITE(6,*)'To the Plate:',IPLATE,' Element:',ITRIA
         ENDIF  
         WRITE(6,'(A53)')
     &   '  IZ0  JR0  JP0  JT0  M_SF IS0-->IS1    TIMET     IC0'
        ENDIF
        NEW = 1
      ELSE
C Old particle:
C the particle is still located in the last cell, on the surface ISRF.
C the particle entries the cell IRGEN
        CALL INTO_CELL()
        IF(IERR.NE.0) THEN
         IF(DIAGNO) CALL WRMESS('ERROR IN INTO_CELL')
         TIMET=-1.
         WRITE(6,*)'*** ERROR FROM INTO_CELL WITH IERR=',IERR
         WRITE(6,'(4I5)')NZ0_N0,JR0_N0,JP0_N0,JT0_N0
         RETURN
        ENDIF
        IRGEN = IC0_N0
      ENDIF
c
c 4. follow particle
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c 5. determine the intersection point
C INPUT: JR0_N0,JP0_N0,JT0_N0,NZ0_N0: CELL ADDRESS
C        XFIN,YFIN,ZFIN   : INITIAL POSITION
C        VX,VY,VZ         : VELOCITY
C        ISRF             : INITILALLY LOCATED SURFACE NUMBER
C        JDRJ,JDPJ,JDTJ   : BOUNDARY SURFACE SEQUENCE
c  definition of the points in the corners of the cell:
       CALL CORNER_POINTS()

c initial surface number
       ISR1 = ISRF
       CALL INTERSECT()
COUTPUT: XFIN,YFIN,ZFIN   : NEW POSITION
C        ISRF             : SURFACE NUMBER OF INTERSECTION
C        TFIN             : TRAVELED DISTANCE
C        IERR >0          : ERROR
       IF(IERR.NE.0) THEN
        IF(DIAGNO) THEN 
          CALL WRMESS('ERROR IN INTERSECT')
          WRITE(6,*)'IERR=',IERR
          WRITE(6,'(7I5,E12.4,I6)')
     +    NZ0_N0,JR0_N0,JP0_N0,JT0_N0,MSURF,ISR1,ISRF,TIMET,IRGEN
        ENDIF
        TIMET=-1.
        WRITE(6,*)'*** ERROR FROM INTERSECT WITH IERR=',IERR
        WRITE(6,'(4I5)')NZ0_N0,JR0_N0,JP0_N0,JT0_N0
        RETURN
       ENDIF

      TIMET  = TFIN + TIMET
C Intersect an additional surface
      IF(TIMET.GT.TLMAX) THEN
         XFIN = XP0 + VX*TLMAX
         YFIN = YP0 + VY*TLMAX
         ZFIN = ZP0 + VZ*TLMAX

         TIMET = TLMAX
         IF(DIAGNO) THEN
          WRITE(6,*)'The Particle reaches a limiter'
          WRITE(6,*)'Distance traveled till now:',TIMET
          WRITE(6,'(A15,3f10.4)')'Final position:',XFIN,YFIN,ZFIN
         ENDIF
         ICNXT=ISR_TYPE(IPLATE)
         ICOS=1
         ISRF = 0
         RETURN
      ENDIF

      IRNXT =-I_JUMP(ISRF)*IDRJ
      IPNXT =-J_JUMP(ISRF)*IDPJ
      ITNXT =-K_JUMP(ISRF)*IDTJ

      CALL CHECK_SF_N0()
      IF(DIAGNO) WRITE(6,'(7I5,E12.4,I6)')
     +NZ0_N0,JR0_N0,JP0_N0,JT0_N0,MSURF,ISR1,ISRF,TIMET,IRGEN

      IF    (MSURF == 0) THEN
C NORMAOL
        ICNXT = IRGEN + 0
      ELSEIF(MSURF <  0) THEN
C Nontransparent surface
        ICNXT = MSURF
        ICOS=1
        ITRIA = 0
        
        IF(DIAGNO) THEN
          IF    (IRNXT /= 0) THEN 
            WRITE(6,'(A31,I4)')
     &      'A NON-TRANSP. R-SURF. WITH IND=',MSURF
            WRITE(6,'(A11,I4)')'SURFACE NR:',JR0_N0+INSRF(IRNXT)
          ELSEIF(IPNXT /= 0) THEN
            WRITE(6,'(A31,I4)')
     &      'A NON-TRANSP. P-SURF. WITH IND=',MSURF
            WRITE(6,'(A11,I4)')'SURFACE NR:',JP0_N0+INSRF(IPNXT)
          ELSE
            WRITE(6,'(A31,I4)')
     &      'A NON-TRANSP. T-SURF. WITH IND=',MSURF
            WRITE(6,'(A11,I4)')'SURFACE NR:',JT0_N0+INSRF(ITNXT)
          ENDIF
        ENDIF

      ELSE
C Non-default surfaces 
        ICNXT = -1
        ICOS=1
        ITRIA = 0

        IF(DIAGNO) THEN
          IF    (IRNXT /= 0) THEN 
            WRITE(6,'(A31,I4)')
     &      'A NON-DEFAULT R-SURF. WITH IND=',MSURF
            WRITE(6,'(A11,I4)')'SURFACE NR:',JR0_N0+INSRF(IRNXT)
          ELSEIF(IPNXT /= 0) THEN
            WRITE(6,'(A31,I4)')
     &      'A NON-DEFAULT P-SURF. WITH IND=',MSURF
            WRITE(6,'(A11,I4)')'SURFACE NR:',JP0_N0+INSRF(IPNXT)
          ELSE
            WRITE(6,'(A31,I4)')
     &      'A NON-DEFAULT T-SURF. WITH IND=',MSURF
            WRITE(6,'(A11,I4)')'SURFACE NR:',JT0_N0+INSRF(ITNXT)
          ENDIF
        ENDIF
      ENDIF

      RETURN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      NAME: NORUSR                                                    C
C  FUNCTION: DETERMINE THE NORMAL VECTOR OF THE REFLECTING SURFACE     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ENTRY NORUSR(ISTS,X0,Y0,Z0,XNORM,YNORM,ZNORM,SCOS,
     .             VELX,VELY,VELZ,IRGEN)
cbp: here: ists, scos not used
C Additional surface ?
      IF(ITRIA /= 0) THEN
        XNORM = V_TRIA_X(ITRIA)
        YNORM = V_TRIA_Y(ITRIA)
        ZNORM = V_TRIA_Z(ITRIA)
C CALLED BY SAMUSR?
        IF(NSR_N0.NE.0) THEN
         XNORM= XNORM*NSR_N0
         YNORM= YNORM*NSR_N0
         ZNORM= ZNORM*NSR_N0

         NSR_N0 = 0
        ELSEIF(XNORM*VX+YNORM*VY+ZNORM*VZ.LT.0.) THEN
         XNORM=-XNORM
         YNORM=-YNORM
         ZNORM=-ZNORM
        ENDIF

        IF(DIAGNO) write(6,*)VELX,VELY,VELZ,XNORM,YNORM,ZNORM
      ELSEIF(MSURF < 0) THEN
C REFLECTING SURFACE ---------------------------------------------
c surface number:
         i=isrf
         do 2 k=1,3
            x3eck(k) = xmy(ipunkt(i,k))
            y3eck(k) = ymy(ipunkt(i,k))
            z3eck(k) = zmy(ipunkt(i,k))
 2       continue
c
C NORMAL VECTOR
         XNORM =(Y3ECK(3)-Y3ECK(1))*(Z3ECK(2)-Z3ECK(1)) -
     +          (Z3ECK(3)-Z3ECK(1))*(Y3ECK(2)-Y3ECK(1))
         YNORM =(Z3ECK(3)-Z3ECK(1))*(X3ECK(2)-X3ECK(1)) -
     +          (X3ECK(3)-X3ECK(1))*(Z3ECK(2)-Z3ECK(1))
         ZNORM =(X3ECK(3)-X3ECK(1))*(Y3ECK(2)-Y3ECK(1)) -
     +          (Y3ECK(3)-Y3ECK(1))*(X3ECK(2)-X3ECK(1))

         VL = DSQRT(XNORM**2+YNORM**2+ZNORM**2)
         XNORM=XNORM/VL
         YNORM=YNORM/VL
         ZNORM=ZNORM/VL

         IF(XNORM*VX+YNORM*VY+ZNORM*VZ.LT.0.) THEN
         XNORM=-XNORM
         YNORM=-YNORM
         ZNORM=-ZNORM
         ENDIF

         IF(DIAGNO) THEN
          WRITE(6,*)'===> NORUSR CALCULATES THE NORMAOL VECTOR'
          WRITE(6,*)'COS=',XNORM*VELX+YNORM*VELY+ZNORM*VELZ
         ENDIF

      ELSEIF( MSURF >  0) THEN
C Periodic surface, leading to a new particle
C Depending on MSURF:

         CALL INTO_CELL()
         IRGEN = IC0_N0
         IF(IERR.NE.0) THEN
            IF(DIAGNO) CALL WRMESS('ERROR IN INTO_CELL')
            TIMET=-1.
            RETURN
         ENDIF

         VELX = VX
         VELY = VY
         VELZ = VZ

         X0   = XFIN
         Y0   = YFIN
         Z0   = ZFIN
      ENDIF

      RETURN

      CONTAINS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE CORNER_POINTS()

       IMPLICIT NONE

C INPUT: JR0_N0,JP0_N0,JT0_N0,NZ0_N0: FINE CELL ADDRESS
C        JDRJ,JDPJ,JDTJ : CORNER POINTS SELECTION SEQUENCE
C
COUTPUT: XMY,YMY,ZMY(8): X- Y- Z-COORDINATES OF CORNER POINTS
C-----------------------------------------------------------------
C      1     I+N_JUMP(IDRJ)    J+N_JUMP(IDPJ)   K+N_JUMP(IDTJ)
C      2     I+N_JUMP(IDRJ)    J+M_JUMP(IDPJ)   K+N_JUMP(IDTJ)
C      3     I+M_JUMP(IDRJ)    J+M_JUMP(IDPJ)   K+N_JUMP(IDTJ)
C      4     I+M_JUMP(IDRJ)    J+N_JUMP(IDPJ)   K+N_JUMP(IDTJ)
C      5     I+N_JUMP(IDRJ)    J+N_JUMP(IDPJ)   K+M_JUMP(IDTJ)
C      6     I+N_JUMP(IDRJ)    J+M_JUMP(IDPJ)   K+M_JUMP(IDTJ)
C      7     I+M_JUMP(IDRJ)    J+M_JUMP(IDPJ)   K+M_JUMP(IDTJ)
C      8     I+M_JUMP(IDRJ)    J+N_JUMP(IDPJ)   K+M_JUMP(IDTJ)
C------------------------------------------------------------------
C  where ID*J=-1,1, indicating the points sequence.
c  N_JUMP(L) = 0,1 and M_JUMP=1,0 when L=-1,1.
c  Standard case: IDRJ=IDPJ=IDTJ= -1.
c
c  The cell is surrounded by 12 triangles defined IPUNKT(I,J)
c  I: surface number (i=1,12)
c  J: point number   (J=1,3)
C
      INTEGER,DIMENSION(8) :: IG_P
      INTEGER :: I,I1,I2,J1,J2,K1,K2,ISTEPJ,ISTEPK,IG_BKGRD



      ISTEPJ=                  SRF_RADI(NZ0_N0)
      ISTEPK= SRF_POLO(NZ0_N0)*SRF_RADI(NZ0_N0)  
c two toroidal planes with M1 and M2
      I1 = JR0_N0 + N_JUMP(IDRJ)
      I2 = JR0_N0 + M_JUMP(IDRJ)

      J1 = JP0_N0 + N_JUMP(IDPJ)
      J2 = JP0_N0 + M_JUMP(IDPJ)

      K1 = JT0_N0 + N_JUMP(IDTJ)
      K2 = JT0_N0 + M_JUMP(IDTJ)

      IG_BKGRD= K1*ISTEPK + GRID_P_OS(NZ0_N0)
      IG_P(1) = I1 + J1*ISTEPJ + IG_BKGRD
      IG_P(2) = I1 + J2*ISTEPJ + IG_BKGRD
      IG_P(3) = I2 + J2*ISTEPJ + IG_BKGRD
      IG_P(4) = I2 + J1*ISTEPJ + IG_BKGRD

      IG_BKGRD= K2*ISTEPK + GRID_P_OS(NZ0_N0)
      IG_P(5) = I1 + J1*ISTEPJ + IG_BKGRD
      IG_P(6) = I1 + J2*ISTEPJ + IG_BKGRD
      IG_P(7) = I2 + J2*ISTEPJ + IG_BKGRD
      IG_P(8) = I2 + J1*ISTEPJ + IG_BKGRD

      DO I=1,4
         XMY(I) = RG(IG_P(I))*COSPHI(K1+PHI_PL_OS(NZ0_N0))
         YMY(I) = RG(IG_P(I))*SINPHI(K1+PHI_PL_OS(NZ0_N0))
         ZMY(I) = ZG(IG_P(I))
      ENDDO
      DO I=5,8
         XMY(I) = RG(IG_P(I))*COSPHI(K2+PHI_PL_OS(NZ0_N0))
         YMY(I) = RG(IG_P(I))*SINPHI(K2+PHI_PL_OS(NZ0_N0))
         ZMY(I) = ZG(IG_P(I))
      ENDDO
      RETURN
      END SUBROUTINE CORNER_POINTS

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     NAME: INTERSECT                                            C
C FUNCTION: DETERMINE THE INTERSECTION POINT OF FLIGHT OF PARTICEC
C           WITH A CELL BOUNDARY                                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INTERSECT()

      IMPLICIT NONE              
C INPUT: JR0_N0,JP0_N0,JT0_N0,NZ0_N0: CELL ADDRESS
C INPUT: XFIN,YFIN,ZFIN : INITIAL POSITION
C INPUT: VX,VY,VZ       : VELOCITY
C INPUT: ISRF           : INITILALLY LOCATED SURFACE NUMBER
C INPUT: JDRJ,JDPJ,JDTJ : BOUNDARY SURFACE SEQUENCE

COUTPUT: XFIN,YFIN,ZFIN : NEW POSITION
C        ISRF           : SURFACE NUMBER OF INTERSECTION
C        TFIN           : TRAVELED DISTANCE
C        IERR >0        : ERROR

C  parameters used in this program only
      REAL*8, PARAMETER :: TOLER=1.E-12
      REAL*8, DIMENSION(12) :: 
     R        xy,       xz,       yz 
      INTEGER,DIMENSION(12) ::
     I        J1,       J2,       J3
      INTEGER :: NMIN,N_CROSS,ICHECK,indmin,MINI,IN
      REAL*8  :: zeitmn,XCR,YCR,ZCR,CHECK_1,CHECK_2
C--------------------------------------------------------------------
      ierr = 0

c side surfaces:
      DO I=1,12
        J1(I)=ipunkt(i,1)
        J2(I)=ipunkt(i,2)
        J3(I)=ipunkt(i,3)
      ENDDO

      DO i=1,12

        xy(i) = (xmy(j2(i))-xmy(j1(i)))*(ymy(j3(i))-ymy(j1(i)))
     +-         (ymy(j2(i))-ymy(j1(i)))*(xmy(j3(i))-xmy(j1(i)))

        xz(i) = (zmy(j2(i))-zmy(j1(i)))*(xmy(j3(i))-xmy(j1(i)))
     +-         (xmy(j2(i))-xmy(j1(i)))*(zmy(j3(i))-zmy(j1(i)))

        yz(i) = (ymy(j2(i))-ymy(j1(i)))*(zmy(j3(i))-zmy(j1(i)))
     +-         (zmy(j2(i))-zmy(j1(i)))*(ymy(j3(i))-ymy(j1(i)))
      ENDDO
C
c DISTAN: Distance of the point (xfin,yfin,zfin) to the plane
c SINFI = SIN(FI): FI= the incidence angle
      DO i=1,12
       DISTAN(I)=- (xfin-xmy(j1(I)))*yz(I)
     &           - (yfin-ymy(j1(I)))*xz(I)
     &           - (zfin-zmy(j1(I)))*xy(I)
       SINFI(I) = vx*yz(I) + vy*xz(I) + vz*xy(I)
      ENDDO

c IF DABS(SINFI).LT.TOLER: either the plane does not exist or
c                    v parallel to the plane
      nmin=0

C     IF(DIAGNO)write(6,'(6E12.4)')SINFI(1:12),DISTAN(1:12)
C     IF(DIAGNO)write(6,'(3E12.4)')VX,VY,VZ 
      IF(ISRF/=0) DISTAN(ISRF) = 0.
      IF( ISRF >= 9) THEN
        IF    (ISRF == 9) THEN
          DISTAN(10) = 0.
        ELSEIF(ISRF ==10) THEN
          DISTAN( 9) = 0.
        ELSEIF(ISRF ==11) THEN
          DISTAN(12) = 0.
        ELSE
          DISTAN(11) = 0.
        ENDIF   
      ENDIF

      DO I=1,12
       IF(SINFI(I)*DISTAN(I) > 0.) THEN
        TCROSS(I) = DISTAN(I)/SINFI(I)
        nmin = nmin + 1
        iebene(nmin) = i
       ENDIF
      ENDDO
c=================================================================
c define the smallest time:
c-----------------------------------------------------------------
c check whether the point is "real" or "virtual":
C TOTAL NUMBER OF THE INTERSECTION POINTS: N_CROSS = nmin
      N_CROSS = nmin
      IF(N_CROSS.eq. 0) THEN
        IERR = 1
        RETURN
      ENDIF

      DO 100 ICHECK=1,N_CROSS

      indmin = iebene(1)
      zeitmn = tcross(indmin)
      mini   = 1

      do i=2,nmin
         in = iebene(i)
         if(tcross(in).lt.zeitmn) then
           zeitmn = tcross(in)
           indmin = in
           mini   = i
         endif
      enddo

      xcr = xfin + vx*zeitmn
      ycr = yfin + vy*zeitmn
      zcr = zfin + vz*zeitmn
c define a triangle:
      do i=1,3
      x3eck(i) = xmy(ipunkt(indmin,i)) - xcr
      y3eck(i) = ymy(ipunkt(indmin,i)) - ycr
      z3eck(i) = zmy(ipunkt(indmin,i)) - zcr
      enddo
c
      vpr1(1) = y3eck(1)*z3eck(2) - z3eck(1)*y3eck(2)
      vpr1(2) = z3eck(1)*x3eck(2) - x3eck(1)*z3eck(2)
      vpr1(3) = x3eck(1)*y3eck(2) - y3eck(1)*x3eck(2)
c
      vpr2(1) = y3eck(2)*z3eck(3) - y3eck(3)*z3eck(2)
      vpr2(2) = z3eck(2)*x3eck(3) - z3eck(3)*x3eck(2)
      vpr2(3) = x3eck(2)*y3eck(3) - x3eck(3)*y3eck(2)
c
c first check:
      CHECK_1=vpr1(1)*vpr2(1)+vpr1(2)*vpr2(2)+vpr1(3)*vpr2(3)
c     IF(DIAGNO) WRITE(6,'(A8,E12.4)')'CHECK_1=',CHECK_1
      IF(CHECK_1.GT.0.) THEN
c
        vpr3(1) = y3eck(3)*z3eck(1) - z3eck(3)*y3eck(1)
        vpr3(2) = z3eck(3)*x3eck(1) - x3eck(3)*z3eck(1)
        vpr3(3) = x3eck(3)*y3eck(1) - y3eck(3)*x3eck(1)
c
c last check:
        CHECK_2=vpr3(1)*vpr2(1)+vpr3(2)*vpr2(2)+vpr3(3)*vpr2(3)
c       IF(DIAGNO )WRITE(6,'(A8,E12.4)')'CHECK_2=',CHECK_2
        IF(CHECK_2.GT.0.) THEN

         XFIN = XCR
         YFIN = YCR
         ZFIN = ZCR 
         TFIN = ZEITMN
         ISRF = INDMIN
         RETURN
        ENDIF
      ENDIF

      nmin=nmin-1
      do i=mini,nmin
        iebene(I) = iebene(I+1)
      enddo

100   CONTINUE
      IERR = 2
      RETURN
      END SUBROUTINE INTERSECT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     NAME: CHECK_SF_N0                                          C
C FUNCTION: CHECK THE BOUNDARY SURFACE FOR NEUTRAL               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE CHECK_SF_N0()
       USE SURFACE_PL

       IMPLICIT NONE             

       INTEGER :: ISURF
C INPUT: JR0_N0,JP0_N0,JT0_N0,NZ0_N0 : CELL ADDRESS
C        IRNXT,IPNXT,ITNXT: R,P,T JUMP INDEX

COUTPUT: MSURF          =0: NORMAL SURAFCE, PARTICLE CONTINUES
C                       <0: NON-TRANSPARENT SURFACE, |MSURF| SURFACE NUMBER
C-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C R1: Radial surface
      IF(IRNXT.NE.0) THEN
        ISURF= JR0_N0 + INSRF(IRNXT) 
     .+       (JP0_N0+JT0_N0*ZON_POLO(NZ0_N0))*SRF_RADI(NZ0_N0) 
     .+        NRS_OFF(NZ0_N0)
        MSURF = IDSURR(ISURF)
C P1: Poloidal surface
      ELSEIF(IPNXT.NE.0) THEN
        ISURF= JR0_N0 + (JP0_N0+INSRF(IPNXT)
     .+        JT0_N0*SRF_POLO(NZ0_N0))*ZON_RADI(NZ0_N0) 
     .+        NPS_OFF(NZ0_N0)

        MSURF = IDSURP(ISURF)
C T1: Toroidal surface
      ELSEIF(ITNXT.NE.0) THEN
        ISURF= JR0_N0 + (JP0_N0+(JT0_N0+INSRF(ITNXT))
     .*        ZON_POLO(NZ0_N0))*ZON_RADI(NZ0_N0) 
     .+        NTS_OFF(NZ0_N0)

        MSURF = IDSURT(ISURF)
      ENDIF
      RETURN
      END SUBROUTINE CHECK_SF_N0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     NAME: INTO_CELL                                            C
C FUNCTION: DETERMINE THE CELL NUMBER INTO WHICH A PARTICLE      C
C           TRAVELS. IF NECESSARY, THE FLIGHT CAN BE CHANGED,    C
C           POSITION CAN BE SHIFTED ACCORDING TO THE SYMMETY OF  C
C           THE GIVEN GEOMETRY                                   C
C           THIS SUBROUTINE IS STANDARD FOR W7-AS AND -X CONF.   C
C           IN OTHER CASES, IT MUST BY REPLACED BY THE USER.     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE INTO_CELL()
       USE MAPPING_TOROIDAL

       IMPLICIT NONE               
C INPUT: JR0_N0,JP0_N0,JT0_N0,NZ0_N0  : CELL ADDRESS
C        XFIN,YFIN,ZFIN   : POSITION
C        VX,VY,VZ         : VELOCITY
C        MSURF            : Surface Type (see CHECK_SF_N0)
C        ISRF             : LOCATED SURFACE NUMBER
C        JDRJ,JDPJ,JDTJ   : BOUNDARY SURFACE SEQUENCE
C        IRNXT,IPNXT,ITNXT: R,P,T JUMP INDEX

COUTPUT: XFIN,YFIN,ZFIN   : NEW POSITION
C        ISRF             : NEW SURFACE NUMBER
C        JR0,JP0,JT0,NZ0  : NEW CELL ADDRESS
C        VX,VY,VZ         : NEW VELOCITY
C        JDRJ,JDPJ,JDTJ   : NEW BOUNDARY SURFACE SEQUENCE
C        IERR >0          : ERROR
C                           THE SYMMETRY OF THE GIVEN GEOMETRY
       INTEGER :: IG,IS,MAP_SF,PAI_SF,TSURF,IPHI
       REAL*8  :: R0,Z0,V_R,V_PHI,RJ,PJ,TJ
C-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
      IERR = 0
C Case0: Normal surface
      IF    (MSURF.EQ.0) THEN
        JR0_N0 = JR0_N0 + IRNXT
        JP0_N0 = JP0_N0 + IPNXT
        JT0_N0 = JT0_N0 + ITNXT
        IG     = JR0_N0 + (JP0_N0+JT0_N0*ZON_POLO(NZ0_N0))
     .*          ZON_RADI(NZ0_N0) + MESH_P_OS(NZ0_N0)
        IC0_N0 = IDCELL(IG)
        IF( OUT_OF_RANGE(IC0_N0,1,NCELL_N2) ) IERR = 1
        
        ISRF  = NEWSRF(ISRF)
        RETURN
      ENDIF 
         
      IF    ( IRNXT /=0 ) THEN
        CALL R_SF_JUMP(MSURF,IRNXT
     .,      NZ0_N0,JR0_N0,JP0_N0,JT0_N0,0.,0.,0. )
        ISRF  = NEWSRF(ISRF)
      ELSEIF( IPNXT /=0 ) THEN
        CALL P_SF_JUMP(MSURF,IPNXT
     .,      NZ0_N0,JR0_N0,JP0_N0,JT0_N0,0.,0.,0. )
        ISRF  = NEWSRF(ISRF)
      ELSEIF( ITNXT /=0 ) THEN
        TSURF = JT0_N0 + INSRF(ITNXT)

        IF    ( MSURF == 1 ) THEN 
C Periodic
          R0    = SQRT(XFIN**2+YFIN**2)
          IPHI  = TSURF  + PHI_PL_OS(NZ0_N0)
          V_R   = VX*COSPHI(IPHI) + VY*SINPHI(IPHI)
          V_PHI =-VX*SINPHI(IPHI) + VY*COSPHI(IPHI)

          IF(DIAGNO) THEN
           WRITE(6,*)'**>Periodic toroidal surface:'
           WRITE(6,*)'   IPHI  XFIN  YFIN  ZFIN  VX  VY  VZ'
           WRITE(6,'(A4,I4,6F9.3)')'old:'
     .,                  IPHI,XFIN,YFIN,ZFIN,VX,VY,VZ
          ENDIF


          CALL T_SF_JUMP(MSURF,ITNXT
     .,        NZ0_N0,JR0_N0,JP0_N0,JT0_N0,0.,0.,0. )

          TSURF = JT0_N0 + INSRF(-ITNXT)
          IPHI  = TSURF  + PHI_PL_OS(NZ0_N0)

          XFIN  = R0*COSPHI(IPHI)
          YFIN  = R0*SINPHI(IPHI)

          VX   = V_R*COSPHI(IPHI) - V_PHI*SINPHI(IPHI)
          VY   = V_R*SINPHI(IPHI) + V_PHI*COSPHI(IPHI)

          IF(DIAGNO) THEN
           WRITE(6,'(A4,I4,6F9.3)')'new:'
     .,                  IPHI,XFIN,YFIN,ZFIN,VX,VY,VZ
          ENDIF
          ISRF  = NEWSRF(ISRF)
        ELSEIF( MSURF == 2 ) THEN
C Up/down symmetric surface
          IPHI  = TSURF  + PHI_PL_OS(NZ0_N0)
          V_R   = VX*COSPHI(IPHI) + VY*SINPHI(IPHI)
          V_PHI =-VX*SINPHI(IPHI) + VY*COSPHI(IPHI)

          IF(DIAGNO) THEN
           WRITE(6,*)'**>Asymmetric toroidal surface:'
           WRITE(6,*)'   IPHI  XFIN  YFIN  ZFIN  VX  VY  VZ'
           WRITE(6,'(A4,I4,6F9.3)')'old:'
     .,                  IPHI,XFIN,YFIN,ZFIN,VX,VY,VZ
          ENDIF
          ZFIN  =-ZFIN
          V_PHI =-V_PHI

          VX    = V_R*COSPHI(IPHI) - V_PHI*SINPHI(IPHI)
          VY    = V_R*SINPHI(IPHI) + V_PHI*COSPHI(IPHI)
          VZ    =-VZ

          IF(DIAGNO) THEN
           WRITE(6,'(A4,I4,6F9.3)')'new:'
     .,                  IPHI,XFIN,YFIN,ZFIN,VX,VY,VZ
          ENDIF

          CALL T_SF_JUMP(MSURF,ITNXT
     .,        NZ0_N0,JR0_N0,JP0_N0,JT0_N0,0.,0.,0. )

        ELSEIF( MSURF == 3 ) THEN 
C Mapping 
          MAP_SF = 0
          DO IS=1,TOTAL_MAP_SF_T
            IF(  NZ0_N0 == ZONE_NR_MAP_T(IS) .AND. 
     .           TSURF  == T_SF_NR_MAP_T(IS) ) THEN 
              MAP_SF = IS
              EXIT
            ENDIF
          ENDDO
          IF( MAP_SF == 0 ) THEN
            WRITE(6,*)'Mapping surface not found'
            CALL STOP_ALL('Stopped in INTO_CELL')
          ENDIF 

          IF(DIAGNO) THEN
           WRITE(6,*)'**>Mapping toroidal surface:'
           WRITE(6,*)'MAP NR.  XFIN  YFIN  ZFIN  VX  VY  VZ'
           WRITE(6,'(A4,I4,6F9.3)')'old:'
     .,                MAP_SF,XFIN,YFIN,ZFIN,VX,VY,VZ
          ENDIF

          PAI_SF = MAP_PAIR_SF_T(MAP_SF)
          JT0_N0 = T_SF_NR_MAP_T(PAI_SF)

          IPHI  = TSURF + PHI_PL_OS(NZ0_N0)
          V_R   = VX*COSPHI(IPHI) + VY*SINPHI(IPHI)
          V_PHI =-VX*SINPHI(IPHI) + VY*COSPHI(IPHI)

          R0 = SQRT(XFIN**2+YFIN**2)
          Z0 = ZFIN
          CALL PERFORM_MAPPING_TOROIDAL_N0
     .    (MAP_SF,NZ0_N0,JR0_N0,JP0_N0,RJ,PJ,R0,Z0,IERR)

          IF(IERR > 0 .OR. IERR==-3) THEN  
            IERR = 2
            RETURN
          ELSEIF(IERR < 0) THEN
            RJ = -0.99999
            IF(IERR == -2) RJ=-RJ
            TJ = JT0_N0
            IPHI  = JT0_N0 + PHI_PL_OS(NZ0_N0)

            CALL RZ_REAL_COORDINATES
     .          (NZ0_N0,JR0_N0,JP0_N0,RJ,PJ,TJ,R0,Z0)

            XFIN  = R0*COSPHI(IPHI)
            YFIN  = R0*SINPHI(IPHI)
            ZFIN  = Z0
            IERR  = 0
          ENDIF

C Mapping on itself
          IF( PAI_SF == MAP_SF ) THEN
            ZFIN  =-ZFIN
            V_PHI =-V_PHI

            VX    = V_R*COSPHI(IPHI) - V_PHI*SINPHI(IPHI)
            VY    = V_R*SINPHI(IPHI) + V_PHI*COSPHI(IPHI)
            VZ    =-VZ
C Mapping on the last toroidal surfaces
          ELSEIF((JT0_N0 == 0                .AND. NZ0_N0 ==  0 ) .OR.
     .           (JT0_N0 == ZON_TORO(NZ0_N0) .AND. NZ0_N0 == NZONET-1)
     .          )THEN
            IPHI  = JT0_N0 + PHI_PL_OS(NZ0_N0)

            XFIN  = R0*COSPHI(IPHI)
            YFIN  = R0*SINPHI(IPHI)

            VX   = V_R*COSPHI(IPHI) - V_PHI*SINPHI(IPHI)
            VY   = V_R*SINPHI(IPHI) + V_PHI*COSPHI(IPHI)
            ISRF  = NEWSRF(ISRF)
          ELSE  
            ISRF  = NEWSRF(ISRF)
          ENDIF
          IF(JT0_N0 == ZON_TORO(NZ0_N0)) JT0_N0 = JT0_N0 - 1

          IF(DIAGNO) THEN
           WRITE(6,'(A4,I4,6F9.3)')'new:'
     .,                PAI_SF,XFIN,YFIN,ZFIN,VX,VY,VZ
          ENDIF
        ENDIF 

      ENDIF

      IG     = JR0_N0 + (JP0_N0+JT0_N0*ZON_POLO(NZ0_N0))
     .*        ZON_RADI(NZ0_N0) + MESH_P_OS(NZ0_N0)
      IC0_N0 = IDCELL(IG)
      IF( OUT_OF_RANGE(IC0_N0,1,NCELL_N2) ) IERR = 1

      RETURN
      END SUBROUTINE INTO_CELL

      END SUBROUTINE TIMUSR
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     NAME: TLMAX_TRAVEL                                         C
C FUNCTION: DETERMINE THE the travel distance of a particle to   C
C           a plate                                              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TLMAX_TRAVEL
      USE NEUTRAL_TRANSPORT 
      
      IMPLICIT NONE            
C INPUT: XP0,YP0,ZP0 : Particle location
C        VX,VY,VZ    : Velocity
      INTEGER :: IP,N_PL,IPM,MI,NMIN,N_CROSS,INDMIN,MINI,IN,ICHECK
     I,          I,J,L,LL,N,NY,N_YA,NY_B,L1,L2,IPOSI,ITR
      REAL*8  :: T_MIN,T_MI,T_MAX,ZEITMN,CHECK_1,CHECK_2,XCR,YCR,ZCR
     R,          T1,T2,DX1,DX2,DY1,DY2,DZ1,DZ2,XR1,XR2,YR1,YR2,ZR1,ZR2 

C Defaul value
      TLMAX = 1.D30

      IF(NTRIANG.LT.10) THEN
        CALL TLMAX_TRAVEL1
      ELSE
        CALL TLMAX_TRAVEL2
      ENDIF
      RETURN

      CONTAINS 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TLMAX_TRAVEL1
      IMPLICIT NONE            
C INPUT: XP0,YP0,ZP0 : Particle location
C        VX,VY,VZ    : Velocity
C--------------------------------------------------------------------
C parameters used in this program only

C case1: if NTRIANG is mall, scan all the triangles
c DISTAN: Distance of the point (xfin,yfin,zfin) to the plane
c SINFI = SIN(FI): FI= the incidence angle
      DO I=1,NTRIANG
       DISTAN(I)=- (XP0-X_TRIA(1,I))*V_TRIA_X(I)
     &           - (YP0-Y_TRIA(1,I))*V_TRIA_Y(I)
     &           - (ZP0-Z_TRIA(1,I))*V_TRIA_Z(I)

       SINFI(I) =   VX*V_TRIA_X(I)
     &           +  VY*V_TRIA_Y(I)
     &           +  VZ*V_TRIA_Z(I)
      ENDDO

c IF SINFI=0.: either the plane does not exist or
c              v parallel to the plane
      nmin=0
      DO I=1,NTRIANG
       IF(SINFI(I).NE.0. .AND. I.NE.ITRIA) THEN
        TCROSS(I) = DISTAN(I)/SINFI(I)
        IF(TCROSS(I).GT.0.) THEN
         nmin = nmin + 1
         iebene(nmin) = i
        ENDIF
       ENDIF
      ENDDO
c=================================================================
c define the smallest time:
c-----------------------------------------------------------------
c check whether the point is "real" or "virtual":
C TOTAL NUMBER OF THE INTERSECTION POINTS: N_CROSS = nmin
      N_CROSS = nmin
      IPLATE= 0
      ITRIA = 0
      IF(N_CROSS == 0) RETURN

      CHECK_LOOP : DO ICHECK=1,N_CROSS

      indmin = iebene(1)
      zeitmn = tcross(indmin)
      mini   = 1

      do i=2,nmin
         in = iebene(i)
         if(tcross(in).lt.zeitmn) then
           zeitmn = tcross(in)
           indmin = in
           mini   = i
         endif
      enddo

      xcr = XP0 + vx*zeitmn
      ycr = YP0 + vy*zeitmn
      zcr = ZP0 + vz*zeitmn
c define a triangle:
      do i=1,3
      x3eck(i) = X_TRIA(I,indmin) - xcr
      y3eck(i) = Y_TRIA(I,indmin) - ycr
      z3eck(i) = Z_TRIA(I,indmin) - zcr
      enddo
c
      vpr1(1) = y3eck(1)*z3eck(2) - z3eck(1)*y3eck(2)
      vpr1(2) = z3eck(1)*x3eck(2) - x3eck(1)*z3eck(2)
      vpr1(3) = x3eck(1)*y3eck(2) - y3eck(1)*x3eck(2)
c
      vpr2(1) = y3eck(2)*z3eck(3) - y3eck(3)*z3eck(2)
      vpr2(2) = z3eck(2)*x3eck(3) - z3eck(3)*x3eck(2)
      vpr2(3) = x3eck(2)*y3eck(3) - x3eck(3)*y3eck(2)
c
c first check:
      CHECK_1=vpr1(1)*vpr2(1)+vpr1(2)*vpr2(2)+vpr1(3)*vpr2(3)
      IF(CHECK_1.GE.0.) THEN
c
        vpr3(1) = y3eck(3)*z3eck(1) - z3eck(3)*y3eck(1)
        vpr3(2) = z3eck(3)*x3eck(1) - x3eck(3)*z3eck(1)
        vpr3(3) = x3eck(3)*y3eck(1) - y3eck(3)*x3eck(1)
c
c last check:
        CHECK_2=vpr3(1)*vpr2(1)+vpr3(2)*vpr2(2)+vpr3(3)*vpr2(3)
        IF(CHECK_2.GE.0.) THEN

         TLMAX  = zeitmn
         ITRIA  = indmin
         IPLATE = NP_NUM(indmin)
         EXIT CHECK_LOOP
        ENDIF
      ENDIF

      nmin=nmin-1
      do i=mini,nmin
        iebene(I) = iebene(I+1)
      enddo

      ENDDO CHECK_LOOP 

      RETURN
      END SUBROUTINE TLMAX_TRAVEL1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TLMAX_TRAVEL2
      
      IMPLICIT NONE            
C INPUT: XP0,YP0,ZP0 : Particle location
C        VX,VY,VZ    : Velocity
C--------------------------------------------------------------------
C parameters used in this program only

      DO IP=1,NPSORT
       DISTAN(IP)=-(XP0-X_SORT(IP))*V_SORT_X(IP)
     &            -(YP0-Y_SORT(IP))*V_SORT_Y(IP)
     &            -(ZP0-Z_SORT(IP))*V_SORT_Z(IP)

       SINFI(IP) =  VX*V_SORT_X(IP)
     &           +  VY*V_SORT_Y(IP)
     &           +  VZ*V_SORT_Z(IP)
      ENDDO

c IF SINFI=0.: either the plane does not exist or
c              v parallel to the plane
      N_PL=0
      DO IP=1,NPSORT
       IF(SINFI(IP).NE.0.) THEN
        TCR_PL(1,IP)  = (DISTAN(IP)+D_MINI(IP))/SINFI(IP)
        TCR_PL(2,IP)  = (DISTAN(IP)+D_MAXI(IP))/SINFI(IP)
        IF(TCR_PL(1,IP) > 0. .OR. TCR_PL(2,IP) > 0.) THEN
         N_PL         = N_PL + 1
         IEB_PL(N_PL) = IP
        ENDIF
       ENDIF
      ENDDO
c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*
c     write(6,*)'N_PL=',N_PL
      IF(N_PL.EQ.0) THEN 
        IPLATE= 0
        ITRIA = 0
        RETURN 
      ENDIF
C smallst distance
250   IPM   = IEB_PL(1)
      T_MIN = MIN(TCR_PL(1,IPM),TCR_PL(2,IPM))
      MI    = 1

      DO I=2,N_PL
         IP = IEB_PL(I)
         T_MI = MIN(TCR_PL(1,IP),TCR_PL(2,IP))
         IF(T_MI.LT.T_MIN) THEN
           T_MIN = T_MI
           IPM   = IP
           MI    = I
         ENDIF
      ENDDO

c     write(6,*) MI,IPM,T_MIN

      N_YA  = 1
      IEB_YA(1) = MI

      T_MAX     = MAX(TCR_PL(1,IPM),TCR_PL(2,IPM))
      DO I=1,N_PL
        IF(I.NE.MI) THEN
         IP = IEB_PL(I)
         T_MI = MIN(TCR_PL(1,IP),TCR_PL(2,IP))

         IF(T_MI.LE.T_MAX) THEN
           N_YA = N_YA + 1
           IEB_YA(N_YA) = I
         ENDIF
        ENDIF
      ENDDO

      LL = 0
      DO 300 N=1,N_YA
       IPOSI = IEB_PL(IEB_YA(N))
       T1 = TCR_PL(1,IPOSI)
       T2 = TCR_PL(2,IPOSI)

       DX1 =  VX*T1
       DY1 =  VY*T1
       DZ1 =  VZ*T1
       DX2 =  VX*T2
       DY2 =  VY*T2
       DZ2 =  VZ*T2

       XR1 = XP0 + MAX(DX1,DX2)
       YR1 = YP0 + MAX(DY1,DY2)
       ZR1 = ZP0 + MAX(DZ1,DZ2)
       XR2 = XP0 + MIN(DX1,DX2)
       YR2 = YP0 + MIN(DY1,DY2)
       ZR2 = ZP0 + MIN(DZ1,DZ2)

       IF( (XR1.LT.X_SMIN(IPOSI) .OR. XR2.GT.X_SMAX(IPOSI))
     & .OR.(YR1.LT.Y_SMIN(IPOSI) .OR. YR2.GT.Y_SMAX(IPOSI))
     & .OR.(ZR1.LT.Z_SMIN(IPOSI) .OR. ZR2.GT.Z_SMAX(IPOSI))
     &   ) GOTO 300

       L1 = LL + 1
       DO I=NPS_GRP(IPOSI-1)+1,NPS_GRP(IPOSI)
         LL        = LL + 1
         MTRIA(LL) = NPN_GRP(I)
       ENDDO
       L2 = LL
       NYBACK(L1:L2) = N

       DO I=L1,L2
         DISTAN(I)=- (XP0-X_TRIA(1,MTRIA(I)))*V_TRIA_X(MTRIA(I))
     &             - (YP0-Y_TRIA(1,MTRIA(I)))*V_TRIA_Y(MTRIA(I))
     &             - (ZP0-Z_TRIA(1,MTRIA(I)))*V_TRIA_Z(MTRIA(I))

          SINFI(I)=  VX*V_TRIA_X(MTRIA(I))
     &            +  VY*V_TRIA_Y(MTRIA(I))
     &            +  VZ*V_TRIA_Z(MTRIA(I))
       ENDDO
300   CONTINUE

      NY_B=0

      IF(LL.EQ.0) GOTO 900

c IF SINFI=0.: either the plane does not exist or
c              v parallel to the plane

       NMIN = 0
       DO I=1,LL
        IF(SINFI(I).NE.0. .AND. MTRIA(I).NE.ITRIA) THEN
        TCROSS(I) = DISTAN(I)/SINFI(I)
        IF(TCROSS(I).GT.0.) THEN
         NMIN         = NMIN + 1
         IEBENE(NMIN) = I
        ENDIF
       ENDIF
      ENDDO
      IF(NMIN.EQ.0) GOTO 900

260   indmin = iebene(1)
      zeitmn = tcross(indmin)
      mini   = 1

      do i=2,nmin
         in = iebene(i)
         if(tcross(in).lt.zeitmn) then
           zeitmn = tcross(in)
           indmin = in
           mini   = i
         endif
      enddo
      ITR = MTRIA(indmin)

      xcr = XP0 + vx*zeitmn
      ycr = YP0 + vy*zeitmn
      zcr = ZP0 + vz*zeitmn
c define a triangle:
      do i=1,3
      x3eck(i) = X_TRIA(I,ITR) - xcr
      y3eck(i) = Y_TRIA(I,ITR) - ycr
      z3eck(i) = Z_TRIA(I,ITR) - zcr
      enddo
c
      vpr1(1) = y3eck(1)*z3eck(2) - z3eck(1)*y3eck(2)
      vpr1(2) = z3eck(1)*x3eck(2) - x3eck(1)*z3eck(2)
      vpr1(3) = x3eck(1)*y3eck(2) - y3eck(1)*x3eck(2)
c
      vpr2(1) = y3eck(2)*z3eck(3) - y3eck(3)*z3eck(2)
      vpr2(2) = z3eck(2)*x3eck(3) - z3eck(3)*x3eck(2)
      vpr2(3) = x3eck(2)*y3eck(3) - x3eck(3)*y3eck(2)
c
c first check:
      CHECK_1=vpr1(1)*vpr2(1)+vpr1(2)*vpr2(2)+vpr1(3)*vpr2(3)
      IF(CHECK_1.GE.0.) THEN
c
        vpr3(1) = y3eck(3)*z3eck(1) - z3eck(3)*y3eck(1)
        vpr3(2) = z3eck(3)*x3eck(1) - x3eck(3)*z3eck(1)
        vpr3(3) = x3eck(3)*y3eck(1) - y3eck(3)*x3eck(1)
c
c last check:
        CHECK_2=vpr3(1)*vpr2(1)+vpr3(2)*vpr2(2)+vpr3(3)*vpr2(3)
        IF(CHECK_2.GE.0.) THEN

         IF(zeitmn.GT.T_MAX) THEN
           NY_B = NY_B + 1
           NYSTOR(NY_B) = NYBACK(indmin)
         ELSE
           TLMAX  = zeitmn
           ITRIA  = ITR
           IPLATE = NP_NUM(ITRIA)
           return
         ENDIF
        ENDIF
      ENDIF
      nmin=nmin-1
      do i=mini,nmin
        iebene(I) = iebene(I+1)
      enddo
      IF(nmin.GT.0) GOTO 260

c900   IF(DIAGNO) THEN
c      WRITE(6,'(A12,I4,4X,16I3)')'Initial: NU=',N_PL,IEB_PL(1:N_PL)
c      WRITE(6,'(A12,I4,4X,16I3)')' Checks: NU=',N_YA,IEB_YA(1:N_YA)
c      WRITE(6,'(A12,I4,4X,16I3)')' virtau: NU=',NY_B,NYSTOR(1:NY_B)
c      ENDIF
900   DO 901 NY=1,N_YA
       DO L=1,NY_B
       IF(NY.EQ.NYSTOR(L)) GOTO 901
       ENDDO
       I=IEB_YA(NY)
       IEB_PL(I) = 0
901   CONTINUE

c      IF(DIAGNO)
c    &WRITE(6,'(A12,I4,4X,16I3)')'Correct: NU=',N_PL,IEB_PL(1:N_PL)

      I=1
902   IF(IEB_PL(I).EQ.0) THEN
         N_PL=N_PL-1
         IF(N_PL.LT.I) GOTO 903
         DO J=I,N_PL
          IEB_PL(J) = IEB_PL(J+1)
         ENDDO
      ELSE
         I = I+1
         IF(I.GT.N_PL) GOTO 903
      ENDIF
      GOTO 902

c903   IF(DIAGNO) THEN
c      WRITE(6,'(A12,I4,4X,16I3)')'  Final: NU=',N_PL,IEB_PL(1:N_PL)
c      WRITE(6,*)'zeitmn,T_MAX=',zeitmn,T_MAX
c      stop
c      ENDIF
903   IF(N_PL.EQ.0) THEN
        IPLATE= 0
        ITRIA = 0
        RETURN
      ELSE
        GOTO 250
      ENDIF
      END SUBROUTINE TLMAX_TRAVEL2

      END SUBROUTINE TLMAX_TRAVEL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      NAME: LEAUSR                                                    C
C  FUNCTION: DETERMINE THE CELL NUMBERS OF A GIVEN POINT(X0,Y0,Z0)     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION LEAUSR(X0,Y0,Z0)

      USE GEOMETRY_PL
      USE PHYSICAL_CELL
      USE NEUTRAL_TRANSPORT

      IMPLICIT NONE            
      INTEGER ::  LEAUSR
     I,           J,IRGEN,NEW,ICNXT,NPANU,ICOS,IER,ISTS
     I,           ITR_STO,NTR_STO
     I,           IG,IG1,IG2,JP,JP1,JP2,JPD,L1,L2,ID,JD
      REAL*8  ::  X0,Y0,Z0,VELX,VELY,VELZ,TIMET
     R,           R_0,Z_0,PHI,PHI1,PHI2,RJ,PJ,TJ,DIST,DISTA
     R,           R,Z,XNORM,YNORM,ZNORM,SCOS
     R,           R_CENTER,Z_CENTER,D1,D2,D3,A,B,C1,C2
      LOGICAL ::  SURF_P

      LOGICAL:: DIAGNO=.FALSE.
      SAVE
C-----------------------------------------------------------------------
C  INPUT:
C     (X0,Y0,Z0 ):  A GIVEN POINT P0
C OUTPUT:
C          LEAUSR:  CELL NUMBER
C-----------------------------------------------------------------------

      LEAUSR = 0

      IF(DIAGNO) THEN
        WRITE(6,*)'====> LEAUSR DETERMINES CELL NUMBER'
        WRITE(6,'(A7,1p,3E12.4)')' INPUT:',X0,Y0,Z0
      ENDIF

      IF( N0S_DEF == 0) GOTO 200

      R_0  = SQRT(X0**2 + Y0**2)
      Z_0  = Z0

C  Phi-plane on which the particle is located
      IER = 1
      PHI = ATAN2(Y0,X0)
      LOOP_NZ0 : DO NZ0_N0=0,NZONET-1
      LOOP_JT0 : DO JT0_N0=0,ZON_TORO(NZ0_N0)-1
        PHI1 = PHI_PLANE(PHI_PL_OS(NZ0_N0)+JT0_N0  )
        PHI2 = PHI_PLANE(PHI_PL_OS(NZ0_N0)+JT0_N0+1)
        IF( PHI>= PHI1 .AND. PHI<= PHI2) THEN 
          IER = 0
          EXIT LOOP_NZ0
        ENDIF 
      ENDDO LOOP_JT0
      ENDDO LOOP_NZ0

      IF( IER/=0) THEN
        IF(DIAGNO)
     .  WRITE(6,*)'PHI,PHI_0,PHI_1=',PHI,PHI_PLANE(0)
     .,            PHI_PLANE(PHI_PL_OS(NZONET))
        RETURN
      ENDIF 

      JR0_N0 = (ZON_RADI(NZ0_N0)-1)/2
      
      IG1=  JT0_N0*SRF_POLO(NZ0_N0)*SRF_RADI(NZ0_N0)
     .+     GRID_P_OS(NZ0_N0)
      IG2= IG1 + SRF_RADI(NZ0_N0)*ZON_POLO(NZ0_N0)/2

      R_CENTER = 0.5*(RG(IG1)+RG(IG2))
      Z_CENTER = 0.

      A =  Z_0 - Z_CENTER
      B =-(R_0 - R_CENTER)

      JPD = (ZON_POLO(NZ0_N0)-1)/6
      JP1 = 0

      IG1 = JR0_N0 + (JP1+JT0_N0*SRF_POLO(NZ0_N0))*SRF_RADI(NZ0_N0)
     .+               GRID_P_OS(NZ0_N0)
      D1 = (RG(IG1)-R_CENTER)*A + (ZG(IG1)-Z_CENTER)*B
      C1 =-(RG(IG1)-R_CENTER)*B + (ZG(IG1)-Z_CENTER)*A

      LOOP_SREACH_JP1_2: DO L1=0,ZON_POLO(NZ0_N0)
        JP2 = JP1 + JPD
        IF(JP2 > ZON_POLO(NZ0_N0)) JP2 = ZON_POLO(NZ0_N0)
        
        IG2 = JR0_N0+(JP2+JT0_N0*SRF_POLO(NZ0_N0))*SRF_RADI(NZ0_N0)
     .+               GRID_P_OS(NZ0_N0)
        D2  = (RG(IG2)-R_CENTER)*A + (ZG(IG2)-Z_CENTER)*B
        C2  =-(RG(IG2)-R_CENTER)*B + (ZG(IG2)-Z_CENTER)*A

      
        IF( D1*D2.LE.0. .AND. (C1>0. .OR. C2>0.) ) THEN

         IG = IG1
         LOOP_SREACH_JP0_N0 : DO 
          IF(JP2-JP1<=1) THEN
           IF(-(RG(IG)-R_CENTER)*B+(ZG(IG)-Z_CENTER)*A > 0.) THEN
             JP0_N0 = JP1
             EXIT LOOP_SREACH_JP1_2 
           ELSE
             EXIT LOOP_SREACH_JP0_N0
           ENDIF
          ENDIF 

          JP= (JP1+JP2)/2
          IG= JR0_N0+(JP+JT0_N0*SRF_POLO(NZ0_N0))*SRF_RADI(NZ0_N0)
     .+                GRID_P_OS(NZ0_N0)
          D3= (RG(IG)-R_CENTER)*A + (ZG(IG)-Z_CENTER)*B
          IF    (D1*D3 <= 0.) THEN 
             JP2 = JP
             D2  = D3
          ELSEIF(D2*D3 <= 0.) THEN
             JP1 = JP
             D1  = D3
          ELSE  
             IER = 1
             EXIT LOOP_SREACH_JP1_2
          ENDIF
         ENDDO LOOP_SREACH_JP0_N0

        ELSEIF(  JP2==ZON_POLO(NZ0_N0) ) THEN  
             IER = 2
             EXIT LOOP_SREACH_JP1_2
        ENDIF 
        D1 = D2
        C1 = C2
        IG1=IG2
        JP1=JP2
      ENDDO LOOP_SREACH_JP1_2

      IF( IER/=0) RETURN 

C Now, we have JR0,JP0,JT0
C Find the grid point close to the (X0,Y0,Z0)
      IG= JR0_N0+(JP0_N0+JT0_N0*SRF_POLO(NZ0_N0))*SRF_RADI(NZ0_N0)
     .+           GRID_P_OS(NZ0_N0)

      DISTA = (RG(IG)-R_0)**2 + (ZG(IG)-Z_0)**2

      ID = 1
      JD = 1
      LOOP_SMALLEST_D: DO

       IF(JR0_N0 == 0 .OR. JR0_N0 == ZON_RADI(NZ0_N0)-1 
     ..OR.JP0_N0 == 0 .OR. JP0_N0 == ZON_POLO(NZ0_N0)-1) 
     . EXIT LOOP_SMALLEST_D

       IG1 = IG

       IG2 = IG + ID
       DIST= (RG(IG2)-R_0)**2 + (ZG(IG2)-Z_0)**2
       IF(DIST < DISTA) THEN
         DISTA = DIST
         JR0_N0= JR0_N0 + ID
         IG    = IG2
       ELSE
         IG2 = IG - ID
         DIST= (RG(IG2)-R_0)**2 + (ZG(IG2)-Z_0)**2
         IF(DIST < DISTA) THEN
           DISTA = DIST
           JR0_N0= JR0_N0 - ID
           IG    = IG2
           ID    =-ID
         ENDIF
       ENDIF

       IG2 = IG + SRF_RADI(NZ0_N0)*JD
       DIST= (RG(IG2)-R_0)**2 + (ZG(IG2)-Z_0)**2
       IF(DIST < DISTA) THEN
         DISTA = DIST
         JP0_N0= JP0_N0 + JD
         IG    = IG2
       ELSE
         IG2 = IG - SRF_RADI(NZ0_N0)*JD
         DIST= (RG(IG2)-R_0)**2 + (ZG(IG2)-Z_0)**2
         IF(DIST < DISTA) THEN
           DISTA = DIST
           JP0_N0= JP0_N0 - JD
           IG    = IG2
           JD    =-JD
         ENDIF
       ENDIF
       
       IF(IG == IG1) EXIT LOOP_SMALLEST_D
      ENDDO LOOP_SMALLEST_D
C*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
200   NTR_STO = NTRIANG
      ITR_STO = ITRIA
      NTRIANG = 0

      RJ = 0.
      PJ = 0.
      TJ = JT0_N0 + 0.00001               

      CALL RZ_REAL_COORDINATES
     .     (NZ0_N0,JR0_N0,JP0_N0,RJ,PJ,TJ,R,Z)
      IG1 = PHI_PL_OS(NZ0_N0) + JT0_N0
      IG2 = IG1 + 1
      PHI = PHI_PLANE(IG1)+(TJ-JT0_N0)*(PHI_PLANE(IG2)-PHI_PLANE(IG1))

      XP0 = R*COS(PHI)
      YP0 = R*SIN(PHI)
      ZP0 = Z
  
      DISTA= SQRT((X0-XP0)**2+(Y0-YP0)**2+(Z0-ZP0)**2)
      VELX = (X0-XP0)/DISTA
      VELY = (Y0-YP0)/DISTA
      VELZ = (Z0-ZP0)/DISTA

      SURF_P = .FALSE.
      IC0_N0 = IDCELL(JR0_N0+(JP0_N0+JT0_N0*ZON_POLO(NZ0_N0))
     .*               ZON_RADI(NZ0_N0) + MESH_P_OS(NZ0_N0)   )
    
      NEW = 0

      IF(DIAGNO) THEN
        WRITE(6,*)'DISTANCE TO (X0,Y0,Z0):',DISTA,' CM'
        WRITE(6,*)'NZ0   JR0  JP0  JT0     DIST(CM)'
      ENDIF

      ICOS  = 1
      
      IF(DIAGNO) ICOS = 2
      TRACE_LOOP : DO 
            CALL TIMUSR(IRGEN,XP0,YP0,ZP0,VELX,VELY,VELZ,NEW,
     .                  ICNXT,TIMET,ICOS,IER,IRGEN,SURF_P)

                 IF(DIAGNO) WRITE(6,'(4I5,F12.6)')
     .           NZ0_N0,JR0_N0,JP0_N0,JT0_N0,TIMET

            IF(TIMET > DISTA) THEN
              LEAUSR = IC0_N0
              IF(DIAGNO) WRITE(6,*) '****LEAUSR=',IC0_N0
              ISRF  = 0
              EXIT TRACE_LOOP
            ENDIF

            IF(TIMET <= 0.) THEN
              WRITE(6,*)'IERR In CALLING TIMUSR'
              
              EXIT TRACE_LOOP
            ELSEIF(ICNXT == -1) THEN
              CALL NORUSR(ISTS,XP0,YP0,ZP0,XNORM,YNORM,ZNORM
     .,                   SCOS,VELX,VELY,VELZ,IRGEN)
              DISTA = DISTA - TIMET
              SURF_P=.TRUE.
              NEW   = 0
            ELSEIF(ICNXT <0 ) THEN
              WRITE(6,*)'A non-transparent surface'
              WRITE(6,*)'MSURF=',MSURF
              EXIT TRACE_LOOP
            ENDIF 
      END DO TRACE_LOOP
      NTRIANG = NTR_STO 
      ITRIA   = ITR_STO
      RETURN
      END FUNCTION LEAUSR
