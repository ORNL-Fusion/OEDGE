!pb  25.10.06:  format specifications corrected
 
      SUBROUTINE EIRENE_OUTSPEC
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CCONA
      USE EIRMOD_CTRCEI
      USE EIRMOD_CTEXT
      USE EIRMOD_CSDVI
 
      IMPLICIT NONE
      INTEGER :: IADTYP(0:4)
      INTEGER :: IOUT, ISPC, I, IT, IE
      REAL(DP) :: EN
      CHARACTER(10) :: TEXTYP(0:4)
 
C  SPECTRA
 
      IOUT = 20
      OPEN (UNIT=IOUT,FILE='spectra.out')
 
      TEXTYP(0) = 'PHOTONS   '
      TEXTYP(1) = 'ATOMS     '
      TEXTYP(2) = 'MOLECULES '
      TEXTYP(3) = 'TEST IONS '
      TEXTYP(4) = 'BULK IONS '
      IADTYP(0:4) = (/ 0, NSPH, NSPA, NSPAM, NSPAMI /)
 
      DO ISPC=1,NADSPC
        I = ESTIML(ISPC)%PSPC%ISPCSRF
        IT = ESTIML(ISPC)%PSPC%ISPCTYP
 
        WRITE (IOUT,*)
        WRITE (IOUT,*)
        WRITE (IOUT,*)
 
        IF (ESTIML(ISPC)%PSPC%ISRFCLL == 0)  THEN
          IF (I > NLIM) THEN
            WRITE (IOUT,'(A,A,I6)') ' SPECTRUM CALCULATED FOR',
     .                     ' NONDEFAULT STANDARD SURFACE ',I-NLIM
          ELSE
            WRITE (IOUT,'(A,A,I6)') ' SPECTRUM CALCULATED FOR',
     .                     ' ADDITIONAL SURFACE ',I
          END IF
          IF (ESTIML(ISPC)%PSPC%IDIREC > 0) THEN
            WRITE (iunout,'(A,3(ES12.4,A1))')
     .      ' IN DIRECTION (',ESTIML(ISPC)%PSPC%SPCVX,',',
     .      ESTIML(ISPC)%PSPC%SPCVY,',',ESTIML(ISPC)%PSPC%SPCVZ,')'
          END IF
          IF (IT == 1) THEN
            WRITE (IOUT,'(A,A)') ' TYPE OF SPECTRUM : ',
     .                'INCIDENT PARTICLE FLUX IN AMP/BIN(EV)   '
          ELSE IF (IT == 2) THEN
            WRITE (IOUT,'(A,A)') ' TYPE OF SPECTRUM : ',
     .                'INCIDENT ENERGY FLUX IN WATT/BIN(EV)    '
          END IF
 
        ELSE IF (ESTIML(ISPC)%PSPC%ISRFCLL == 1)  THEN
          WRITE (IOUT,'(A,A,I6)') ' SPECTRUM CALCULATED FOR',
     .                   ' SCORING CELL ',I
          IF (ESTIML(ISPC)%PSPC%IDIREC > 0) THEN
            WRITE (iunout,'(A,3(ES12.4,A1))')
     .      ' IN DIRECTION (',ESTIML(ISPC)%PSPC%SPCVX,',',
     .      ESTIML(ISPC)%PSPC%SPCVY,',',ESTIML(ISPC)%PSPC%SPCVZ,')'
          END IF
          IF (IT == 1) THEN
            WRITE (iunout,'(A20,A)') ' TYPE OF SPECTRUM : ',
     .        'SPECTRAL PARTICLE DENSITY IN #/CM**3/BIN(EV)   '
          ELSEIF (IT == 2) THEN
            WRITE (iunout,'(A20,A)') ' TYPE OF SPECTRUM : ',
     .        'SPECTRAL ENERGY DENSITY IN EV/CM**3/BIN(EV)    '
          ELSEIF (IT == 3) THEN
            WRITE (iunout,'(A20,A)') ' TYPE OF SPECTRUM : ',
     .        'SPECTRAL MOMENTUM DENSITY IN (G*CM/S)/CM**3/BIN(EV)    '
          END IF
 
        ELSE IF (ESTIML(ISPC)%PSPC%ISRFCLL == 2)  THEN
          WRITE (IOUT,'(A,A)') ' SPECTRUM CALCULATED FOR',
     .                   ' GEOMETRICAL CELL ',I
          IF (ESTIML(ISPC)%PSPC%IDIREC > 0) THEN
            WRITE (iunout,'(A,3(ES12.4,A1))')
     .      ' IN DIRECTION (',ESTIML(ISPC)%PSPC%SPCVX,',',
     .      ESTIML(ISPC)%PSPC%SPCVY,',',ESTIML(ISPC)%PSPC%SPCVZ,')'
          END IF
          IF (IT == 1) THEN
            WRITE (iunout,'(A20,A)') ' TYPE OF SPECTRUM : ',
     .        'SPECTRAL PARTICLE DENSITY IN #/CM**3/BIN(EV)   '
          ELSEIF (IT == 2) THEN
            WRITE (iunout,'(A20,A)') ' TYPE OF SPECTRUM : ',
     .        'SPECTRAL ENERGY DENSITY IN EV/CM**3/BIN(EV)    '
          ELSEIF (IT == 3) THEN
            WRITE (iunout,'(A20,A)') ' TYPE OF SPECTRUM : ',
     .        'SPECTRAL MOMENTUM DENSITY IN (G*CM/S)/CM**3/BIN(EV)    '
          END IF
        END IF
 
        WRITE (IOUT,'(A20,A9)') ' TYPE OF PARTICLE : ',
     .         TEXTYP(ESTIML(ISPC)%PSPC%IPRTYP)
        IF (ESTIML(ISPC)%PSPC%IPRSP == 0) THEN
          WRITE (IOUT,'(A10,10X,A16)') ' SPECIES :',
     .                'SUM OVER SPECIES'
        ELSE
          WRITE (IOUT,'(A10,10X,A8)') ' SPECIES :',
     .          TEXTS(IADTYP(ESTIML(ISPC)%PSPC%IPRTYP)+
     .          ESTIML(ISPC)%PSPC%IPRSP)
        END IF
 
        WRITE (IOUT,'(A15,5X,ES12.4)') ' MINIMAL ENERGY ',
     .         ESTIML(ISPC)%PSPC%SPCMIN
        WRITE (IOUT,'(A15,5X,ES12.4)') ' MAXIMAL ENERGY ',
     .         ESTIML(ISPC)%PSPC%SPCMAX
        WRITE (IOUT,'(A16,4x,I6)') ' NUMBER OF BINS ',
     .         ESTIML(ISPC)%PSPC%NSPC
        WRITE (IOUT,*)
        IF (ESTIML(ISPC)%PSPC%SPCINT > EPS60) THEN
          IF (NSIGI_SPC == 0) THEN
            DO IE=1, ESTIML(ISPC)%PSPC%NSPC
              EN = ESTIML(ISPC)%PSPC%SPCMIN +
     .             (IE-0.5)*ESTIML(ISPC)%PSPC%SPCDEL
              WRITE (IOUT,'(I6,2ES12.4)') IE,EN,
     .               ESTIML(ISPC)%PSPC%SPC(IE)
            END DO
          ELSE
            DO IE=1, ESTIML(ISPC)%PSPC%NSPC
              EN = ESTIML(ISPC)%PSPC%SPCMIN +
     .             (IE-0.5)*ESTIML(ISPC)%PSPC%SPCDEL
              WRITE (IOUT,'(I6,3ES12.4)') IE,EN,
     .               ESTIML(ISPC)%PSPC%SPC(IE),
     .               ESTIML(ISPC)%PSPC%SDV(IE)
            END DO
          END IF
        ELSE
          WRITE (IOUT,'(A)') ' SPECTRUM IDENTICAL 0 '
        END IF
        WRITE (IOUT,*)
        WRITE (IOUT,'(A,ES12.4)') ' INTEGRAL OF SPECTRUM ',
     .         ESTIML(ISPC)%PSPC%SPCINT
        IF (NSIGI_SPC > 0)
     .    WRITE (IOUT,'(A,ES12.4)') ' STANDARD DEVIATION  ',
     .                   ESTIML(ISPC)%PSPC%SGMS
      END DO
 
      RETURN
      END SUBROUTINE EIRENE_OUTSPEC
