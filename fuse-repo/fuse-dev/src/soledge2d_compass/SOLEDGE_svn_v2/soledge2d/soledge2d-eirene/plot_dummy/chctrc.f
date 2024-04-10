c------------------------------------------------------------------------
      SUBROUTINE EIRENE_CHCTRC(XPLO,YPLO,ZPLO,IFLAG,ISYM)
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CLOGAU
      USE EIRMOD_CTRCEI
      USE EIRMOD_COMPRT
      USE EIRMOD_COMSOU
      USE EIRMOD_CCONA
      USE EIRMOD_CUPD
 
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
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) TXTHST(ISYM)
        IF (ISYM.EQ.1)  CALL EIRENE_MASJ1('NPANU   ',NPANU)
        WRITE (iunout,'(1X,A8)') TEXTS(ISPZ)
        CALL EIRENE_MASJ4 ('ITIME,IFPATH,IUPDTE,ICOL        ',
     .               ITIME,IFPATH,IUPDTE,ICOL)
        CALL EIRENE_MASR3 ('X0,Y0,Z0                ',XPLO,YPLO,ZPLO)
        CALL EIRENE_MASR6
     .  ('VELX,VELY,VELZ,VEL,E0,WEIGHT                    ',
     .               VELX,VELY,VELZ,VEL,E0,WEIGHT)
        CALL EIRENE_MASR1 ('TIME    ',TIME)
        CALL EIRENE_MASJ2 ('IGENCV,IGENCS   ',IGENCV,IGENCS)
        CALL EIRENE_MASJ4 ('NRCELL,IPOLG,NACELL,NBLOCK      ',
     .               NRCELL,IPOLG,NACELL,NBLOCK)
        IF (NLTOR) THEN
          CALL EIRENE_MASJ1 ('NTCELL  ',NTCELL)
        ENDIF
        IF (NLTRA) THEN
          CALL EIRENE_MASJ1R ('NNTCLL,PHI      ',NNTCLL,PHI/DEGRAD)
        ENDIF
        IF (NLPOL) THEN
          CALL EIRENE_MASJ1 ('NPCELL  ',NPCELL)
        ENDIF
        IF ((ISYM.GE.6.AND.ISYM.LE.10).OR.
     .      (ISYM.EQ.1.AND.NLSRF(ISTRA))) THEN
          CALL EIRENE_MASJ5 ('MRSURF,MPSURF,MTSURF,MASURF,NLLI        ',
     .                 MRSURF,MPSURF,MTSURF,MASURF,NLLI)
          CALL EIRENE_MASR1 ('SCOS    ',SCOS)
        ENDIF
      ENDIF
 
      RETURN
      END
