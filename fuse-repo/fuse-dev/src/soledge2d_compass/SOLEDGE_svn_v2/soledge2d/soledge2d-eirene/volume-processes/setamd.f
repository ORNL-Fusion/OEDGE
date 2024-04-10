C 27.6.05:  PHV_NROTA, PHV_NROTPH REMOVED
C
      SUBROUTINE EIRENE_SETAMD(ICAL)
C
C  SET ATOMIC AND MOLECULAR DATA: DRIVER
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMXS
      USE EIRMOD_COMSOU
      USE EIRMOD_COMUSR
      USE EIRMOD_CZT1
      USE EIRMOD_PHOTON
 
      IMPLICIT NONE
 
      real(dp) :: tpb1, tpb2, second_own
      INTEGER, INTENT(IN) :: ICAL
      INTEGER :: I, IRPI, IRDS
 
!pb      tpb1 = second_own()
 
      IF (ICAL == 0) THEN
        NRCX=0
        NREL=0
        NRPI=0
        NRDS=0
        NREC=0
        NBGV=0
        NROT=0
        CALL EIRENE_XSECTA_PARAM
        CALL EIRENE_XSECTM_PARAM
        CALL EIRENE_XSECTI_PARAM
        CALL EIRENE_XSECTP_PARAM
        CALL EIRENE_XSECTPH_PARAM
 
        NRCX=MAX(1,NRCX)
        NREL=MAX(1,NREL)
        NRPI=MAX(1,NRPI)
        NRDS=MAX(1,NRDS)
        NREC=MAX(1,NREC)
        NBGV=MAX(1,NBGV)
        NROT=MAX(1,NROT)
 
        CALL EIRENE_SET_PARMMOD(2)
        CALL EIRENE_ALLOC_COMXS(2)
        CALL EIRENE_ALLOC_COMSOU(2)
        CALL EIRENE_ALLOC_CZT1(2)

        MAXSPC(0:4) = (/ NPHOTI,NATMI,NMOLI,NIONI,NPLSI /)
 
!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for setamd(0) ',tpb2-tpb1
!pb        tpb1 = tpb2
 
        RETURN
      ELSE
        CALL EIRENE_INIT_CMDTA(2)
 
!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for init_cmdta ',tpb2-tpb1
!pb        tpb1 = tpb2
 
      END IF
 
      NRCXI=0
      NRELI=0
      NRPII=0
      NREII=0
      NRRCI=0
      NRBGI=0
      CALL EIRENE_XSECTA
 
!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for xsecta ',tpb2-tpb1
!pb        tpb1 = tpb2
 
      CALL EIRENE_XSECTM
 
!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for xsectm ',tpb2-tpb1
!pb        tpb1 = tpb2
 
      CALL EIRENE_XSECTI
 
!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for xsecti ',tpb2-tpb1
!pb        tpb1 = tpb2
 
      CALL EIRENE_XSECTP
 
!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for xsectp ',tpb2-tpb1
!pb        tpb1 = tpb2
 
      CALL EIRENE_XSECTPH
 
!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for xsectph ',tpb2-tpb1
!pb        tpb1 = tpb2
 
      CALL EIRENE_CONDENSE
 
!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for condense ',tpb2-tpb1
!pb        tpb1 = tpb2
 
 
      IPATDS = 0
      IPMLDS = 0
      IPIODS = 0
      IPPLDS = 0
      DO IRDS=1,NRDS
        ipatds(IRDS,0)=COUNT(PATDS(IRDS,1:) > 0)
        IF (ipatds(IRDS,0).GT.0) THEN          ! inserted by Derek Harting 26.03
             IPATDS(IRDS,1:ipatds(IRDS,0))=PACK( (/ (i,i=1,natm) /),
     .                                     PATDS(IRDS,1:) > 0)
        END IF
        ipmlds(IRDS,0)=COUNT(PMLDS(IRDS,1:) > 0)
        IF (ipmlds(IRDS,0).GT.0) THEN          ! inserted by Derek Harting 26.03
             IPMLDS(IRDS,1:ipmlds(IRDS,0))=PACK( (/ (i,i=1,nmol) /),
     .                                     PMLDS(IRDS,1:) > 0)
        END IF
        ipiods(IRDS,0)=COUNT(PIODS(IRDS,1:) > 0)
        IF (ipiods(IRDS,0).GT.0) THEN         ! inserted by Derek Harting 26.03.
             IPIODS(IRDS,1:ipiods(IRDS,0))=PACK( (/ (i,i=1,nion) /),
     .                                     PIODS(IRDS,1:) > 0)
        END IF
        ipplds(IRDS,0)=COUNT(PPLDS(IRDS,1:) > 0)
        IF (ipplds(IRDS,0).GT.0) THEN         ! inserted by Derek Harting 26.03.
             IPPLDS(IRDS,1:ipplds(IRDS,0))=PACK( (/ (i,i=1,npls) /),
     .                                     PPLDS(IRDS,1:) > 0)
        END IF
      END DO
 
      IPATPI = 0
      IPMLPI = 0
      IPIOPI = 0
      IPPLPI = 0
      DO IRPI=1,NRPI
        ipatpi(IRPI,0)=COUNT(PATPI(IRPI,1:) > 0)
        IF (ipatpi(IRPI,0).GT.0)              ! inserted by Derek Harting 26.03.
     &       IPATPI(IRPI,1:ipatpi(IRPI,0))=PACK( (/ (i,i=1,natm) /),
     .                                     PATPI(IRPI,1:) > 0)
        ipmlpi(IRPI,0)=COUNT(PMLPI(IRPI,1:) > 0)
        IF (ipmlpi(IRPI,0).GT.0)              ! inserted by Derek Harting 26.03.
     &       IPMLPI(IRPI,1:ipmlpi(IRPI,0))=PACK( (/ (i,i=1,nmol) /),
     .                                     PMLPI(IRPI,1:) > 0)
        ipiopi(IRPI,0)=COUNT(PIOPI(IRPI,1:) > 0)
        IF (ipiopi(IRPI,0).GT.0)              ! inserted by Derek Harting 26.03.
     &       IPIOPI(IRPI,1:ipiopi(IRPI,0))=PACK( (/ (i,i=1,nion) /),
     .                                     PIOPI(IRPI,1:) > 0)
        ipplpi(IRPI,0)=COUNT(PPLPI(IRPI,1:) > 0)
        IF (ipplpi(IRPI,0).GT.0)              ! inserted by Derek Harting 26.03.
     &       IPPLPI(IRPI,1:ipplpi(IRPI,0))=PACK( (/ (i,i=1,npls) /),
     .                                     PPLPI(IRPI,1:) > 0)
      END DO
 
!pb        tpb2 = second_own()
!pb        write (6,*) ' cpu time for packs und counts ',tpb2-tpb1
!pb        tpb1 = tpb2
 
 
      RETURN
      END
