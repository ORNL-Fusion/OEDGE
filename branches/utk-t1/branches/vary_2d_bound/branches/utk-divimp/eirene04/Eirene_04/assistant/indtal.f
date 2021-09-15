C
C
C*DK INDTAL
      SUBROUTINE INDTAL (IND,M,NX,NY,NZ,NB)
C
C   called from subr. STATIS
C               subr. STATIS_BGK
C               subr. STATIS_COP
C   SIMILAR TO SUBR. INTTAL:
C   PROVIDE ARRAY IND(J,K), J=1,NRAD, K=1,8
C   THE VALUES IN = IND(J,K) , FOR EACH CELL J,
C   ARE THE INDICES OF THOSE BINS TO WHICH CELL J CONTRIBUTES
C   FOR EACH EIRENE CELL JE: IND(JE,1)=JE,
C
C   2ND INDEX .LE. 8, BECAUSE J MAY CONTRIBUTE TO
C
C    1     J              (FUNCTION OF X,Y,Z)
C    2     X      AVERAGE (FUNCTION OF Y,Z)
C    3     Y      AVERAGE (FUNCTION OF X,Z)
C    4     Z      AVERAGE (FUNCTION OF X,Y)
C    5     X-Y    AVERAGE (FUNCTION OF Z)
C    6     X-Z    AVERAGE (FUNCTION OF Y)
C    7     Y-Z    AVERAGE (FUNCTION OF X)
C    8     X-Y-Z  AVERAGE   --
C
C   M: LEADING DIMENSION OF IND IN CALLING PROGRAM
C  FURTHER (TO SPEED UP SUBR.STATIS)
C
C
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: M, NX, NY, NZ, NB
      INTEGER, INTENT(OUT) :: IND(M,8)
      INTEGER :: IX, IY, IZ, NIR, NIRT, NIRTP, IADD, K, IR, IJ, IK,
     .           NIT, NIPR, NIP, NITP, NS, I, NXM, NYM, NZM, IB, J,
     .           N1DEL, N2DEL
C
C NS: NUMBER OF STANDARD MESH CELLS
      NS=NX*NY*NZ*NB
C
      DO 10 I=1,M
        DO 10 J=1,8
          IND(I,J)=0
10    CONTINUE
C
C  additional cells contribute only to themselves
C
      DO 11 I=NS+1,M
        IND(I,1)=I
11    CONTINUE
C
      IF (NS.EQ.0) RETURN
C
      N1DEL=0
      IF (NY.GT.1.OR.NZ.GT.1) N1DEL=NX
      N2DEL=0
      IF (NZ.GT.1) N2DEL=NY
C
      NXM=MAX(1,NX-1)
      NYM=MAX(1,NY-1)
      NZM=MAX(1,NZ-1)

C  LOOP OVER STANDARD MESH BLOCKS
      DO 1000 IB=1,NB
      IADD=(IB-1)*NX*NY*NZ
C
      NIRTP=NX+((NY-1)+(NZ-1)*N2DEL)*N1DEL+IADD
      DO 100 IY=1,NYM
        NIRT=NX+((IY-1)+(NZ-1)*N2DEL)*N1DEL+IADD
        DO 101 IZ=1,NZM
          NIR=NX+((IY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
          DO 102 IX=1,NXM
            K=IX+((IY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
            IND(K,1)=K
            IND(K,2)=NIR
            IND(K,6)=NIRT
            IND(K,8)=NIRTP
102       CONTINUE
101     CONTINUE
100   CONTINUE
C
      DO 200 IZ=1,NZM
        NIPR=NX+((NY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
        DO 201 IX=1,NXM
          NIP=IX+((NY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
          DO 202 IY=1,NYM
            K=IX+((IY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
            IND(K,1)=K
            IND(K,3)=NIP
            IND(K,5)=NIPR
202       CONTINUE
201     CONTINUE
200   CONTINUE
C
      DO 300 IX=1,NXM
        NITP=IX+((NY-1)+(NZ-1)*N2DEL)*N1DEL+IADD
        DO 301 IY=1,NYM
          NIT=IX+((IY-1)+(NZ-1)*N2DEL)*N1DEL+IADD
          DO 302 IZ=1,NZM
            K=IX+((IY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
            IND(K,1)=K
            IND(K,4)=NIT
            IND(K,7)=NITP
302       CONTINUE
301     CONTINUE
300   CONTINUE
C
1000  CONTINUE
C  teste: gibt es ein ir, fuer verschiedene j aber gleiche ind(ir,j)
C         dann wuerde naemlich in statis doppelt gezaehlt
C         falls weniger als 3 dimensionen, kann sowas vorkommen.
      DO 2000 J=1,8
2000  continue
      do 2001 ir=1,ns
        do j=1,8
          ij=ind(ir,j)
          do k=1,8
            ik=ind(ir,k)
            if (ik.ne.0.and.ik.eq.ij.and.k.ne.j) then
              if (k.gt.j) ind(ir,k)=0
              if (j.gt.k) ind(ir,j)=0
            endif
          enddo
        enddo
2001  continue
C
      do 2002 ir=1,ns
        do j=1,8
          ij=ind(ir,j)
        enddo
2002  continue

C
      RETURN
      END
