c  sept. 05:  five more tallies added to step function, see also CSTEP.f
c  nov.  05:  add eltot and ve to step function data
 
C  write plasma (background) data, source distribution and atomic data
C  on unit 13.
C
c  at entry RPLAM:
C  read plasma (background) data, source distribution and atomic data
C  from unit 13.
C
C  trcfle:  confirm writing on printout on unit IUNOUT
C  IFLG  :
 
      SUBROUTINE EIRENE_WRPLAM(TRCFLE,IFLG)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CZT1
      USE EIRMOD_COMSOU
      USE EIRMOD_CSTEP
      USE EIRMOD_COMXS
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IFLG
      LOGICAL TRCFLE
      REAL(DP), ALLOCATABLE :: RDUM(:)
      INTEGER, ALLOCATABLE :: IDUM(:)
      LOGICAL, ALLOCATABLE :: LDUM(:)
      INTEGER :: NRDUM, NIDUM, NLDUM
C
      OPEN (UNIT=13,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 13
      WRITE (13)
C  REAL
     R           TEIN,TIIN,DEIN,DIIN,VXIN,VYIN,VZIN,
     R           BXIN,BYIN,BZIN,BFIN,ADIN,EDRIFT,
     R           VOL,WGHT,BXPERP,BYPERP,
     R           FLXOUT,SAREA,
     R           TEINL,TIINL,DEINL,DIINL,BVIN,PARMOM,
     R           RMASSI,RMASSA,RMASSM,RMASSP,
     R           DIOD,DATD,DMLD,DPLD,DPHD,
     R           DION,DATM,DMOL,DPLS,DPHOT,
     R           TVAC,DVAC,VVAC,ALLOC,
     T           TEXTS,
C  MUSR, INTEGER
     I           NSPH  ,NPHOTI,NPHOTIM,NFOLPH,NGENPH,
     I           NSPA  ,NATMI,NATMIM,NMASSA,NCHARA,NFOLA,NGENA,
     I           NSPAM ,NMOLI,NMOLIM,NMASSM,NCHARM,NFOLM,NGENM,
     I           NSPAMI,NIONI,NIONIM,NMASSI,NCHARI,NCHRGI,NFOLI,NGENI,
     I           NSPTOT,NPLSI,NPLSIM,NMASSP,NCHARP,NCHRGP,NBITS,
     I           NSNVI,NCPVI,NADVI,NBGVI,NALVI,NCLVI,NADSI,NALSI,NAINI,
     I           NPRT,ISPEZ,ISPEZI,
C  LUSR, LOGICAL
     L           LGVAC,LGDFT
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: module EIRMOD_COMUSR.f '
      CALL EIRENE_WRITE_CMDTA
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: RCMDTA,ICMDTA'
      CALL EIRENE_WRITE_CMAMF
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: RCMAMF,ICMAMF'
      WRITE (13) RCZT1,RCZT2,ZT1,ZRG
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: RCZT1,RCZT2,ZT1,ZRG'
      WRITE (13) RCMSOU,SREC,EIO,EEL,
     .           ICMSOU,INGRDA,INGRDE,NSTRAI,
     .           LCMSOU,NLSYMP,NLSYMT
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: RCMSOU,ICMSOU,LCMSOU'
      IF (ALLOCATED(FLSTEP))
     .  WRITE (13) FLSTEP,ELSTEP,FLTOT,ELTOT,VF,VE,
     .             QUOT,ADD,QUOTI,ADDIV,
     .             TESTEP,TISTEP,RRSTEP,VXSTEP,VYSTEP,VZSTEP,DISTEP,
     .             FESTEP,FISTEP,SHSTEP,VPSTEP,MCSTEP,
     .             IRSTEP,IPSTEP,ITSTEP,IASTEP,IBSTEP,IGSTEP,
     .             ISTUF,NSMAX,NSPSTI,NSPSTE
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: module EIRMOD_CSTEP.f'
      CLOSE (UNIT=13)
      RETURN
C
      ENTRY EIRENE_RPLAM(TRCFLE,IFLG)
      OPEN (UNIT=13,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 13
      READ (13)
C  REAL
     R           TEIN,TIIN,DEIN,DIIN,VXIN,VYIN,VZIN,
     R           BXIN,BYIN,BZIN,BFIN,ADIN,EDRIFT,
     R           VOL,WGHT,BXPERP,BYPERP,
     R           FLXOUT,SAREA,
     R           TEINL,TIINL,DEINL,DIINL,BVIN,PARMOM,
     R           RMASSI,RMASSA,RMASSM,RMASSP,
     R           DIOD,DATD,DMLD,DPLD,DPHD,
     R           DION,DATM,DMOL,DPLS,DPHOT,
     R           TVAC,DVAC,VVAC,ALLOC,
     T           TEXTS,
C  MUSR, INTEGER
     I           NSPH  ,NPHOTI,NPHOTIM,NFOLPH,NGENPH,
     I           NSPA  ,NATMI,NATMIM,NMASSA,NCHARA,NFOLA,NGENA,
     I           NSPAM ,NMOLI,NMOLIM,NMASSM,NCHARM,NFOLM,NGENM,
     I           NSPAMI,NIONI,NIONIM,NMASSI,NCHARI,NCHRGI,NFOLI,NGENI,
     I           NSPTOT,NPLSI,NPLSIM,NMASSP,NCHARP,NCHRGP,NBITS,
     I           NSNVI,NCPVI,NADVI,NBGVI,NALVI,NCLVI,NADSI,NALSI,NAINI,
     I           NPRT,ISPEZ,ISPEZI,
C  LUSR, LOGICAL
     L           LGVAC,LGDFT
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: module EIRMOD_COMUSR.f '
      CALL EIRENE_READ_CMDTA
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: RCMDTA,ICMDTA'
      CALL EIRENE_READ_CMAMF
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: RCMAMF,ICMAMF'
      READ (13) RCZT1,RCZT2,ZT1,ZRG
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: RCZT1,RCZT2,ZT1,ZRG'
      IF (IFLG == 0) THEN
        READ (13) RCMSOU,SREC,EIO,EEL,
     .            ICMSOU,INGRDA,INGRDE,NSTRAI,
     .            LCMSOU,NLSYMP,NLSYMT
      ELSE
        NRDUM = SIZE(RCMSOU) + SIZE(SREC) + SIZE(EIO) + SIZE(EEL)
        NIDUM = SIZE(ICMSOU) + SIZE(INGRDA) + SIZE(INGRDE) + 1
        NLDUM = SIZE(LCMSOU) + SIZE(NLSYMP) + SIZE(NLSYMT)
        ALLOCATE (RDUM(NRDUM))
        ALLOCATE (IDUM(NIDUM))
        ALLOCATE (lDUM(NLDUM))
        READ (13) RDUM, IDUM, LDUM
        DEALLOCATE (RDUM)
        DEALLOCATE (IDUM)
        DEALLOCATE (lDUM)
      END IF
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: RCMSOU,ICMSOU,LCMSOU'
      IF (ALLOCATED(FLSTEP))
     .   READ (13) FLSTEP,ELSTEP,FLTOT,ELTOT,VF,VE,
     .             QUOT,ADD,QUOTI,ADDIV,
     .             TESTEP,TISTEP,RRSTEP,VXSTEP,VYSTEP,VZSTEP,DISTEP,
     .             FESTEP,FISTEP,SHSTEP,VPSTEP,MCSTEP,
     .             IRSTEP,IPSTEP,ITSTEP,IASTEP,IBSTEP,IGSTEP,
     .             ISTUF,NSMAX,NSPSTI,NSPSTE
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: module EIRMOD_CSTEP.f '
      CLOSE (UNIT=13)
      RETURN
      END
