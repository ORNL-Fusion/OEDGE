C
C
C
      SUBROUTINE
     .  EIRENE_PRTTAL(T1,T2,T3,PROF,NR,NP,NT,NB,NTT,IFLAG,IFILE)
C
C  PRINT VOLUME AVERAGED TALLIES
C  IFLAG=-1:  ONLY HEADER IS PRINTED
C  IFLAG= 0:  ONLY MEAN VALUES IN EACH BLOCK
C  IFLAG= 1:  ADDITIONALLY: 1D PROFILES (THESE MAY BE AVERAGES)
C  IFLAG= 2:  ADDITIONALLY: 2D PROFILES (THESE MAY BE AVERAGES)
C  IFLAG= 3:  ADDITIONALLY: 3D PROFILES
C  IFLAG> 3:  ONLY FULL PROFILES, NO AVERAGES
C
C  IFILE> 0:  WRITE FULL TALLY ONTO STREAM FORT.IFILE
C
      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
 
      CHARACTER(*), INTENT(IN) :: T1, T2, T3
      REAL(DP), INTENT(IN) :: PROF(*)
      INTEGER, INTENT(IN) :: NR, NP, NT, NB, NTT, IFLAG, IFILE
      REAL(DP) :: H(6)
      INTEGER :: K(6)
      INTEGER :: JR, JP, JT, IJ, N1DEL, N2DEL, IADD, JA, IA, IB, I,
     .           IC, IT, IP, NRM, NS, NTM, NPM, IRAD, IST, NCOL, IR
      CHARACTER(1) :: TL(72)
 
      DATA TL/72*'='/
      SAVE
 
      CALL EIRENE_LEER(3)
      WRITE (iunout,'(72A1)') TL
      WRITE (iunout,'(72A1)') TL
      WRITE (iunout,'(A9,A72)') 'TALLY:   ',T1
      WRITE (iunout,'(A9,A24)') 'SPECIES: ',T2
      WRITE (iunout,'(A9,A24)') 'UNITS:   ',T3
      WRITE (iunout,'(72A1)') TL
      WRITE (iunout,'(72A1)') TL
      CALL EIRENE_LEER(1)
      IF (IFILE.GT.0) THEN
        OPEN (UNIT=IFILE,position='APPEND')
        WRITE (IFILE,*) TL
        WRITE (IFILE,*) TL
        WRITE (IFILE,*) T1
        WRITE (IFILE,*) T2
        WRITE (IFILE,*) T3
        WRITE (IFILE,*) TL
        WRITE (IFILE,*) TL
        WRITE (IFILE,*) NR,NP,NT,NB,NTT
        DO IRAD=1,NTT,5
          WRITE (IFILE,*) (PROF(IR),IR=IRAD,MIN(IRAD+4,NTT))
        ENDDO
        CLOSE (UNIT=IFILE)
      ENDIF
 
11111 IF (IFLAG.LT.0) RETURN
C  NCOL: NUMBER OF PRINTED DATA PER LINE, .LE.6
      NCOL=5
C
      NS=NR*NP*NT*NB
      NRM=MAX(1,NR-1)
      NPM=MAX(1,NP-1)
      NTM=MAX(1,NT-1)
      N1DEL=0
      IF (NP.GT.1.OR.NT.GT.1) N1DEL=NR
      N2DEL=0
      IF (NT.GT.1) N2DEL=NP
C
C LOOP OVER STANDARD MESH BLOCKS
C
      IF (NS.EQ.0) GOTO 50000
      DO 10000 IB=1,NB
      IF (NB.GT.1) THEN
        WRITE (iunout,*) TL
        WRITE (iunout,777) IB
        WRITE (iunout,*) TL
      ENDIF
      IADD=(IB-1)*NR*NP*NT
C
      IF (IFLAG.LE.2) GOTO 1000
      IF (NR.LE.1.OR.NP.LE.1.OR.NT.LE.1) GOTO 1000
C
C  3 D PROFILES
C
      WRITE (iunout,*) TL
      WRITE (iunout,81)
      WRITE (iunout,*) TL
      CALL EIRENE_LEER(1)
      DO 1 JT=1,NTM
        WRITE (iunout,77) JT
        CALL EIRENE_LEER(1)
        DO 11 JP=1,NPM
          IJ=1
          IR=0
          WRITE (iunout,7) JP
110       DO 111 JR=IJ,NRM
            IC=JR+((JP-1)+(JT-1)*N2DEL)*N1DEL+IADD
            IR=IR+1
            IJ=IJ+1
            K(IR)=JR
            H(IR)=PROF(IC)
            IF (IR.GE.NCOL) GOTO 112
111       CONTINUE
112       CONTINUE
          WRITE (iunout,6) (K(I),H(I),I=1,IR)
          IR=0
          IF (IJ.LE.NRM) GOTO 110
C  NEXT SEGMENT
          CALL EIRENE_LEER(2)
11      CONTINUE
        WRITE (iunout,*) TL
1     CONTINUE
      WRITE (iunout,*) TL
      IF (IFLAG.GT.3) GOTO 10000
C
C  2 D PROFILES
C
1000  CONTINUE
C
      IF (IFLAG.LE.1) GOTO 2000
      IF (NR.LE.1.AND.NP.LE.1) GOTO 2000
      IF (NP.LE.1.AND.NT.LE.1) GOTO 2000
      IF (NR.LE.1.AND.NT.LE.1) GOTO 2000
C
C  POLOIDAL AND TOROIDAL PROFILE, RADIALLY AVERAGED
C
      IF (NTM.GT.1.AND.NPM.GT.1) THEN
        WRITE (iunout,82)
        IF (NR.GT.1) WRITE (iunout,881)
        IF (NR.GT.1) CALL EIRENE_LEER(1)
        DO 2 JT=1,NTM
          WRITE (iunout,77) JT
          IJ=1
          IP=0
220       DO 222 JP=IJ,NPM
            IC=NR+((JP-1)+(JT-1)*N2DEL)*N1DEL+IADD
            IP=IP+1
            IJ=IJ+1
            K(IP)=JP
            H(IP)=PROF(IC)
            IF (IP.GE.NCOL) GOTO 223
222       CONTINUE
223       CONTINUE
          WRITE (iunout,6) (K(I),H(I),I=1,IP)
          IP=0
          IF (IJ.LE.NPM) GOTO 220
C  NEXT SEGMENT
          CALL EIRENE_LEER(2)
2       CONTINUE
        WRITE (iunout,*) TL
      ENDIF
C
C  RADIAL AND POLOIDAL PROFILE, TOROIDALLY AVERAGED
C
      IF (NRM.GT.1.AND.NPM.GT.1) THEN
        WRITE (iunout,81)
        IF (NT.GT.1) WRITE (iunout,883)
        IF (NT.GT.1) CALL EIRENE_LEER(1)
        DO 3 JP=1,NPM
          WRITE (iunout,7) JP
          IJ=1
          IR=0
330       DO 333 JR=IJ,NRM
            IC=JR+((JP-1)+(NT-1)*N2DEL)*N1DEL+IADD
            IR=IR+1
            IJ=IJ+1
            K(IR)=JR
            H(IR)=PROF(IC)
            IF (IR.GE.NCOL) GOTO 334
333       CONTINUE
334       CONTINUE
          WRITE (iunout,6) (K(I),H(I),I=1,IR)
          IR=0
          IF (IJ.LE.NRM) GOTO 330
C  NEXT SEGMENT
          CALL EIRENE_LEER(2)
3       CONTINUE
        WRITE (iunout,*) TL
      ENDIF
C
C  RADIAL AND TOROIDAL PROFILE, POLOIDALLY AVERAGED
C
      IF (NRM.GT.1.AND.NTM.GT.1) THEN
        WRITE (iunout,81)
        IF (NP.GT.1) WRITE (iunout,882)
        IF (NP.GT.1) CALL EIRENE_LEER(1)
        DO 4 JT=1,NTM
          WRITE (iunout,77) JT
          IJ=1
          IR=0
440       DO 444 JR=IJ,NRM
            IC=JR+((NP-1)+(JT-1)*N2DEL)*N1DEL+IADD
            IR=IR+1
            IJ=IJ+1
            K(IR)=JR
            H(IR)=PROF(IC)
            IF (IR.GE.NCOL) GOTO 445
444       CONTINUE
445       CONTINUE
          WRITE (iunout,6) (K(I),H(I),I=1,IR)
          IR=0
          IF (IJ.LE.NRM) GOTO 440
C  NEXT SEGMENT
          CALL EIRENE_LEER(2)
4       CONTINUE
        WRITE (iunout,*) TL
      ENDIF
      IF (IFLAG.GT.3) GOTO 10000
C
C  1 D PROFILES
C
2000  CONTINUE
C
      IF (IFLAG.LE.0) GOTO 3000
C  RADIAL PROFILE, POLOIDALLY AND TOROIDALLY AVERAGED
C
      IF (NRM.GT.1) THEN
        WRITE (iunout,81)
        IF (NP.GT.1.AND.NT.EQ.1) WRITE (iunout,882)
        IF (NP.EQ.1.AND.NT.GT.1) WRITE (iunout,883)
        IF (NP.GT.1.AND.NT.GT.1) WRITE (iunout,8883)
        IJ=1
        IR=0
1110    DO 1111 JR=IJ,NRM
          IC=JR+((NP-1)+(NT-1)*N2DEL)*N1DEL+IADD
          IR=IR+1
          IJ=IJ+1
          K(IR)=JR
          H(IR)=PROF(IC)
          IF (IR.GE.NCOL) GOTO 1112
1111    CONTINUE
1112    CONTINUE
        WRITE (iunout,6) (K(I),H(I),I=1,IR)
        IR=0
        IF (IJ.LE.NRM) GOTO 1110
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) TL
      ENDIF
C
C  POLOIDAL PROFILE, RADIALLY AND TOROIDALLY AVERAGED
C
      IF (NPM.GT.1) THEN
        WRITE (iunout,82)
        IF (NR.GT.1.AND.NT.EQ.1) WRITE (iunout,881)
        IF (NR.EQ.1.AND.NT.GT.1) WRITE (iunout,883)
        IF (NR.GT.1.AND.NT.GT.1) WRITE (iunout,8882)
        IJ=1
        IP=0
1220    DO 1222 JP=IJ,NPM
          IC=NR+((JP-1)+(NT-1)*N2DEL)*N1DEL+IADD
          IP=IP+1
          IJ=IJ+1
          K(IP)=JP
          H(IP)=PROF(IC)
          IF (IP.GE.NCOL) GOTO 1223
1222    CONTINUE
1223    CONTINUE
        WRITE (iunout,6) (K(I),H(I),I=1,IP)
        IP=0
        IF (IJ.LE.NPM) GOTO 1220
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) TL
      ENDIF
C
C  TOROIDAL PROFILE, RADIALLY AND POLOIDALLY AVERAGED
C
      IF (NTM.GT.1) THEN
        IJ=1
        IT=0
        WRITE (iunout,83)
        IF (NR.GT.1.AND.NP.EQ.1) WRITE (iunout,881)
        IF (NR.EQ.1.AND.NP.GT.1) WRITE (iunout,882)
        IF (NR.GT.1.AND.NP.GT.1) WRITE (iunout,8881)
1330    DO 1333 JT=IJ,NTM
          IC=NR+((NP-1)+(JT-1)*N2DEL)*N1DEL+IADD
          IT=IT+1
          IJ=IJ+1
          K(IT)=JT
          H(IT)=PROF(IC)
          IF (IT.GE.NCOL) GOTO 1334
1333    CONTINUE
1334    CONTINUE
        WRITE (iunout,6) (K(I),H(I),I=1,IT)
        IT=0
        IF (IJ.LE.NTM) GOTO 1330
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) TL
      ENDIF
      IF (IFLAG.GT.3) GOTO 10000
C
3000  CONTINUE
      IC=NR+((NP-1)+(NT-1)*N2DEL)*N1DEL+IADD
      WRITE (iunout,8888) PROF(IC)
      WRITE (iunout,*) TL
      CALL EIRENE_LEER(4)
C
10000 CONTINUE
C
50000 CONTINUE
C  ADDITIONAL CELLS
      IF (NTT.GT.NS) THEN
        WRITE (iunout,*) TL
        WRITE (iunout,7777)
        WRITE (iunout,*) TL
      ENDIF
      IJ=NS+1
      IA=0
550   DO 555 JA=IJ,NTT
        IC=JA
        IA=IA+1
        IJ=IJ+1
        K(IA)=JA-NS
        H(IA)=PROF(IC)
        IF (IA.GE.NCOL) GOTO 556
555   CONTINUE
556   CONTINUE
      WRITE (iunout,6) (K(IC),H(IC),IC=1,IA)
      IA=0
      IF (IJ.LE.NTT) GOTO 550
      CALL EIRENE_LEER(2)
C
6     FORMAT (1X,6(I4,2X,1PE12.4,2X))
7     FORMAT (1X,'Y- OR POLOIDAL SEGMENT NUMBER ',I4)
77    FORMAT (1X,'Z- OR TOROIDAL SEGMENT NUMBER ',I4)
777   FORMAT (1X,'STANDARD MESH BLOCK NUMBER ',I4)
7777  FORMAT (1X,'ADDITIONAL CELLS ')
81    FORMAT (1X,'X- OR RADIAL PROFILE ')
82    FORMAT (1X,'Y- OR POLOIDAL PROFILE ')
83    FORMAT (1X,'Z- OR TOROIDAL PROFILE ')
881   FORMAT (1X,'X- OR RADIAL AVERAGE ')
882   FORMAT (1X,'Y- OR POLOIDAL AVERAGE ')
883   FORMAT (1X,'Z- OR TOROIDAL AVERAGE ')
8881  FORMAT (1X,'X- OR RAD. AND Y- OR POL. AVERAGE ',1PE12.4)
8882  FORMAT (1X,'X- OR RAD. AND Z- OR TOR. AVERAGE ',1PE12.4)
8883  FORMAT (1X,'Y- OR POL. AND Z- OR TOR. AVERAGE ',1PE12.4)
8888  FORMAT (1X,'BLOCK AVERAGE ',1PE12.4)
      RETURN
      END
