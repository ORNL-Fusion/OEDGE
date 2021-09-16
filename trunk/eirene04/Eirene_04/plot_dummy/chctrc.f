c------------------------------------------------------------------------
      SUBROUTINE CHCTRC(XPLO,YPLO,ZPLO,IFLAG,ISYM)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CLOGAU
      USE CTRCEI
      USE COMPRT
      USE COMSOU
      USE CCONA
      USE CUPD

      IMPLICIT NONE
C
      INTEGER,PARAMETER :: NTXHST=18

      REAL(DP), INTENT(IN) :: XPLO, YPLO, ZPLO 
      INTEGER, INTENT(IN) :: IFLAG, ISYM
      INTEGER :: IGENCS, IGENCV, NLLI
      CHARACTER(20) :: TXTHST(NTXHST)

      DATA TXTHST
     .           /'LOCATE(1)           ',
     .            'ELECTR. IMPACT(2)   ',
     .            'HEAVY PAR. IMPACT(3)',
     .            'PHOTON IMPACT(4)    ',
     .            'ELASTIC COLL.(5)    ',
     .            'CHARGE EXCHANGE(6)  ',
     .            'FOKKER PLANCK(7)    ',
     .            'SURFACE(8)          ',
     .            'SPLITTING(9)        ',
     .            'RUSSIAN ROULETTE(10)',
     .            'PERIODICITY(11)     ',
     .            'RESTART:A. SPLT.(12)',
     .            'SAVE:COND. EXP.(13) ',
     .            'RESTART:COND EXP(14)',
     .            'TIME LIMIT(15)      ',
     .            'GENERATION LIMIT(16)',
     .            'FLUID LIMIT(17)     ',
     .            'ERROR DETECTED      '/
C
C  WRITE TRACK DATA
C
      IF (TRCHST) THEN
        CALL LEER(1)
        WRITE (6,*) TXTHST(ISYM)
        IF (ISYM.EQ.1)  CALL MASJ1('NPANU   ',NPANU)
        WRITE (6,'(1X,A8)') TEXTS(ISPZ)
        CALL MASJ4 ('ITIME,IFPATH,IUPDTE,ICOL        ',
     .               ITIME,IFPATH,IUPDTE,ICOL)
        CALL MASR3 ('X0,Y0,Z0                ',XPLO,YPLO,ZPLO)
        CALL MASR6 ('VELX,VELY,VELZ,VEL,E0,WEIGHT                    ',
     .               VELX,VELY,VELZ,VEL,E0,WEIGHT)
        CALL MASR1 ('TIME    ',TIME)
        CALL MASJ2 ('IGENCV,IGENCS   ',IGENCV,IGENCS)
        CALL MASJ4 ('NRCELL,IPOLG,NACELL,NBLOCK      ',
     .               NRCELL,IPOLG,NACELL,NBLOCK)
        IF (NLTOR) THEN
          CALL MASJ1 ('NTCELL  ',NTCELL)
        ENDIF
        IF (NLTRA) THEN
          CALL MASJ1R ('NNTCLL,PHI      ',NNTCLL,PHI/DEGRAD)
        ENDIF
        IF (NLPOL) THEN
          CALL MASJ1 ('NPCELL  ',NPCELL)
        ENDIF
        IF ((ISYM.GE.6.AND.ISYM.LE.10).OR.
     .      (ISYM.EQ.1.AND.NLSRF(ISTRA))) THEN
          CALL MASJ5 ('MRSURF,MPSURF,MTSURF,MASURF,NLLI        ',
     .                 MRSURF,MPSURF,MTSURF,MASURF,NLLI)
          CALL MASR1 ('SCOS    ',SCOS)
        ENDIF
      ENDIF

      RETURN
      END
