C
C
      SUBROUTINE SETAMD(ICAL)
C
C  SET ATOMIC AND MOLECULAR DATA: DRIVER
C
      USE PRECISION
      USE PARMMOD
      USE COMXS
      USE COMSOU
      USE CZT1
      USE PHOTON

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ICAL
      INTEGER :: I, IRPI, IRDS

      IF (ICAL == 0) THEN
        NRCX=0
        NREL=0
        NRPI=0
        NRDS=0
        NREC=0
        NBGV=0
        NROT=0
        CALL XSECTA_PARAM
        CALL XSECTM_PARAM
        CALL XSECTI_PARAM
        CALL XSECTP_PARAM
        CALL XSECTPH_PARAM
        
        NRCX=MAX(1,NRCX)
        NREL=MAX(1,NREL)
        NRPI=MAX(1,NRPI)
        NRDS=MAX(1,NRDS)
        NREC=MAX(1,NREC)
        NBGV=MAX(1,NBGV)
        NROT=MAX(1,PHV_NROTA+PHV_NROTPH)
        
        CALL SET_PARMMOD(2)
        CALL ALLOC_COMXS(2)
        CALL ALLOC_COMSOU(2)
        CALL ALLOC_CZT1(2)
        RETURN
      ELSE
        CALL INIT_CMDTA(2)
      END IF

      NRCXI=0
      NRELI=0
      NRPII=0
      NREII=0
      NRRCI=0
      NRBGI=0
      CALL XSECTA
      CALL XSECTM
      CALL XSECTI
      CALL XSECTP
      CALL XSECTPH
      CALL CONDENSE

      IPATDS = 0
      IPMLDS = 0
      IPIODS = 0 
      IPPLDS = 0
      DO IRDS=1,NRDS
        ipatds(IRDS,0)=COUNT(PATDS(IRDS,1:) > 0)
        IF (ipatds(IRDS,0).GT.0) THEN          ! inserted by Derek Harting 26.03.2004
             IPATDS(IRDS,1:ipatds(IRDS,0))=PACK( (/ (i,i=1,natm) /),
     .                                     PATDS(IRDS,1:) > 0)
        ENDIF
        ipmlds(IRDS,0)=COUNT(PMLDS(IRDS,1:) > 0)
        IF (ipmlds(IRDS,0).GT.0) THEN          ! inserted by Derek Harting 26.03.2004
             IPMLDS(IRDS,1:ipmlds(IRDS,0))=PACK( (/ (i,i=1,nmol) /),
     .                                     PMLDS(IRDS,1:) > 0)
        ENDIF
        ipiods(IRDS,0)=COUNT(PIODS(IRDS,1:) > 0)

        IF (ipiods(IRDS,0).GT.0) THEN          ! inserted by Derek Harting 26.03.2004
             IPIODS(IRDS,1:ipiods(IRDS,0))=PACK( (/ (i,i=1,nion) /),
     .                                     PIODS(IRDS,1:) > 0)
        ENDIF
        ipplds(IRDS,0)=COUNT(PPLDS(IRDS,1:) > 0)
        IF (ipplds(IRDS,0).GT.0) THEN          ! inserted by Derek Harting 26.03.2004
             IPPLDS(IRDS,1:ipplds(IRDS,0))=PACK( (/ (i,i=1,npls) /),
     .                                     PPLDS(IRDS,1:) > 0)
        ENDIF
      END DO

      IPATPI = 0
      IPMLPI = 0
      IPIOPI = 0
      IPPLPI = 0
      DO IRPI=1,NRPI
        ipatpi(IRPI,0)=COUNT(PATPI(IRPI,1:) > 0)
        IF (ipatpi(IRPI,0).GT.0) THEN         ! inserted by Derek Harting 26.03.2004
             IPATPI(IRPI,1:ipatpi(IRPI,0))=PACK( (/ (i,i=1,natm) /),
     .                                     PATPI(IRPI,1:) > 0)
        ENDIF
        ipmlpi(IRPI,0)=COUNT(PMLPI(IRPI,1:) > 0)
        IF (ipmlpi(IRPI,0).GT.0) THEN         ! inserted by Derek Harting 26.03.2004
             IPMLPI(IRPI,1:ipmlpi(IRPI,0))=PACK( (/ (i,i=1,nmol) /),
     .                                     PMLPI(IRPI,1:) > 0)
        ENDIF
        ipiopi(IRPI,0)=COUNT(PIOPI(IRPI,1:) > 0)
        IF (ipiopi(IRPI,0).GT.0) THEN         ! inserted by Derek Harting 26.03.2004
             IPIOPI(IRPI,1:ipiopi(IRPI,0))=PACK( (/ (i,i=1,nion) /),
     .                                     PIOPI(IRPI,1:) > 0)
        ENDIF
        ipplpi(IRPI,0)=COUNT(PPLPI(IRPI,1:) > 0)
        IF (ipplpi(IRPI,0).GT.0) THEN         ! inserted by Derek Harting 26.03.2004
             IPPLPI(IRPI,1:ipplpi(IRPI,0))=PACK( (/ (i,i=1,npls) /),
     .                                     PPLPI(IRPI,1:) > 0)
        ENDIF
      END DO

      RETURN
      END
