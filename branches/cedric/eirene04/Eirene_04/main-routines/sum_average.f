C
      SUBROUTINE SUM_AVERAGE (ISTRA)
C
C  MONTE CARLO CALCULATION
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COUTAU
      USE CSPEZ
      USE CESTIM
      USE CLGIN

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ISTRA
      INTEGER :: J, IATM, IMOL, IION, IPLS, IPHOT, ISPZ, IADV, ICLV,
     .           ISNV, ICPV, IBGV
C
C   SUM OVER SURFACE INDEX
C   IN THE SURFACE AVERAGED ESTIMATORS
C
      DO 632 IMOL=1,NMOLI
        IF (.NOT.LOGMOL(IMOL,ISTRA)) CYCLE
        DO 633 J=1,NLIMPS
          IF (ILIIN(J).LE.0) CYCLE
          PRFAMI(IMOL,ISTRA)=PRFAMI(IMOL,ISTRA)+PRFAML(IMOL,J)
          PRFMMI(IMOL,ISTRA)=PRFMMI(IMOL,ISTRA)+PRFMML(IMOL,J)
          PRFIMI(IMOL,ISTRA)=PRFIMI(IMOL,ISTRA)+PRFIML(IMOL,J)
          PRFPHMI(IMOL,ISTRA)=PRFPHMI(IMOL,ISTRA)+PRFPHML(IMOL,J)
          PRFPMI(IMOL,ISTRA)=PRFPMI(IMOL,ISTRA)+PRFPML(IMOL,J)
          POTMLI(IMOL,ISTRA)=POTMLI(IMOL,ISTRA)-POTML(IMOL,J)
          ERFAMI(IMOL,ISTRA)=ERFAMI(IMOL,ISTRA)+ERFAML(IMOL,J)
          ERFMMI(IMOL,ISTRA)=ERFMMI(IMOL,ISTRA)+ERFMML(IMOL,J)
          ERFIMI(IMOL,ISTRA)=ERFIMI(IMOL,ISTRA)+ERFIML(IMOL,J)
          ERFPHMI(IMOL,ISTRA)=ERFPHMI(IMOL,ISTRA)+ERFPHML(IMOL,J)
          ERFPMI(IMOL,ISTRA)=ERFPMI(IMOL,ISTRA)+ERFPML(IMOL,J)
          EOTMLI(IMOL,ISTRA)=EOTMLI(IMOL,ISTRA)-EOTML(IMOL,J)
          SPTMLI(IMOL,ISTRA)=SPTMLI(IMOL,ISTRA)+SPTML(IMOL,J)
633     CONTINUE
632   CONTINUE
C
      DO 635 IATM=1,NATMI
        IF (.NOT.LOGATM(IATM,ISTRA)) CYCLE
        DO 636 J=1,NLIMPS
          IF (ILIIN(J).LE.0) CYCLE
          PRFAAI(IATM,ISTRA)=PRFAAI(IATM,ISTRA)+PRFAAT(IATM,J)
          PRFMAI(IATM,ISTRA)=PRFMAI(IATM,ISTRA)+PRFMAT(IATM,J)
          PRFIAI(IATM,ISTRA)=PRFIAI(IATM,ISTRA)+PRFIAT(IATM,J)
          PRFPHAI(IATM,ISTRA)=PRFPHAI(IATM,ISTRA)+PRFPHAT(IATM,J)
          PRFPAI(IATM,ISTRA)=PRFPAI(IATM,ISTRA)+PRFPAT(IATM,J)
          POTATI(IATM,ISTRA)=POTATI(IATM,ISTRA)-POTAT(IATM,J)
          ERFAAI(IATM,ISTRA)=ERFAAI(IATM,ISTRA)+ERFAAT(IATM,J)
          ERFMAI(IATM,ISTRA)=ERFMAI(IATM,ISTRA)+ERFMAT(IATM,J)
          ERFIAI(IATM,ISTRA)=ERFIAI(IATM,ISTRA)+ERFIAT(IATM,J)
          ERFPHAI(IATM,ISTRA)=ERFPHAI(IATM,ISTRA)+ERFPHAT(IATM,J)
          ERFPAI(IATM,ISTRA)=ERFPAI(IATM,ISTRA)+ERFPAT(IATM,J)
          EOTATI(IATM,ISTRA)=EOTATI(IATM,ISTRA)-EOTAT(IATM,J)
          SPTATI(IATM,ISTRA)=SPTATI(IATM,ISTRA)+SPTAT(IATM,J)
636     CONTINUE
635   CONTINUE
C
      DO 638 IION=1,NIONI
        IF (.NOT.LOGION(IION,ISTRA)) CYCLE
        DO 639 J=1,NLIMPS
          IF (ILIIN(J).LE.0) CYCLE
          PRFAII(IION,ISTRA)=PRFAII(IION,ISTRA)+PRFAIO(IION,J)
          PRFMII(IION,ISTRA)=PRFMII(IION,ISTRA)+PRFMIO(IION,J)
          PRFIII(IION,ISTRA)=PRFIII(IION,ISTRA)+PRFIIO(IION,J)
          PRFPHII(IION,ISTRA)=PRFPHII(IION,ISTRA)+PRFPHIO(IION,J)
          PRFPII(IION,ISTRA)=PRFPII(IION,ISTRA)+PRFPIO(IION,J)
          POTIOI(IION,ISTRA)=POTIOI(IION,ISTRA)-POTIO(IION,J)
          ERFAII(IION,ISTRA)=ERFAII(IION,ISTRA)+ERFAIO(IION,J)
          ERFMII(IION,ISTRA)=ERFMII(IION,ISTRA)+ERFMIO(IION,J)
          ERFIII(IION,ISTRA)=ERFIII(IION,ISTRA)+ERFIIO(IION,J)
          ERFPHII(IION,ISTRA)=ERFPHII(IION,ISTRA)+ERFPHIO(IION,J)
          ERFPII(IION,ISTRA)=ERFPII(IION,ISTRA)+ERFPIO(IION,J)
          EOTIOI(IION,ISTRA)=EOTIOI(IION,ISTRA)-EOTIO(IION,J)
          SPTIOI(IION,ISTRA)=SPTIOI(IION,ISTRA)+SPTIO(IION,J)
639     CONTINUE
638   CONTINUE
C
      DO IPHOT=1,NPHOTI
        IF (.NOT.LOGPHOT(IPHOT,ISTRA)) CYCLE
        DO J=1,NLIMPS
          IF (ILIIN(J).LE.0) CYCLE
          PRFAPHTI(IPHOT,ISTRA)=PRFAPHTI(IPHOT,ISTRA)+PRFAPHT(IPHOT,J)
          PRFMPHTI(IPHOT,ISTRA)=PRFMPHTI(IPHOT,ISTRA)+PRFMPHT(IPHOT,J)
          PRFIPHTI(IPHOT,ISTRA)=PRFIPHTI(IPHOT,ISTRA)+PRFIPHT(IPHOT,J)
          PRFPHPHTI(IPHOT,ISTRA)=PRFPHPHTI(IPHOT,ISTRA)+
     .                           PRFPHPHT(IPHOT,J)
          PRFPPHTI(IPHOT,ISTRA)=PRFPPHTI(IPHOT,ISTRA)+PRFPPHT(IPHOT,J)
          POTPHTI(IPHOT,ISTRA)=POTPHTI(IPHOT,ISTRA)-POTPHT(IPHOT,J)
          ERFAPHTI(IPHOT,ISTRA)=ERFAPHTI(IPHOT,ISTRA)+ERFAPHT(IPHOT,J)
          ERFMPHTI(IPHOT,ISTRA)=ERFMPHTI(IPHOT,ISTRA)+ERFMPHT(IPHOT,J)
          ERFIPHTI(IPHOT,ISTRA)=ERFIPHTI(IPHOT,ISTRA)+ERFIPHT(IPHOT,J)
          ERFPHPHTI(IPHOT,ISTRA)=ERFPHPHTI(IPHOT,ISTRA)+
     .                           ERFPHPHT(IPHOT,J)
          ERFPPHTI(IPHOT,ISTRA)=ERFPPHTI(IPHOT,ISTRA)+ERFPPHT(IPHOT,J)
          EOTPHTI(IPHOT,ISTRA)=EOTPHTI(IPHOT,ISTRA)-EOTPHT(IPHOT,J)
          SPTPHTI(IPHOT,ISTRA)=SPTPHTI(IPHOT,ISTRA)+SPTPHT(IPHOT,J)
        END DO
      END DO
C
      DO 641 IPLS=1,NPLSI
        IF (.NOT.LOGPLS(IPLS,ISTRA)) CYCLE
        DO 642 J=1,NLIMPS
          IF (ILIIN(J).LE.0) CYCLE
          POTPLI(IPLS,ISTRA)=POTPLI(IPLS,ISTRA)-POTPL(IPLS,J)
          EOTPLI(IPLS,ISTRA)=EOTPLI(IPLS,ISTRA)-EOTPL(IPLS,J)
          SPTPLI(IPLS,ISTRA)=SPTPLI(IPLS,ISTRA)+SPTPL(IPLS,J)
642     CONTINUE
641   CONTINUE
C
      DO 651 ISPZ=1,NSPTOT
        DO 652 J=1,NLIMPS
          IF (ILIIN(J).LE.0) CYCLE
            SPUMPI(ISPZ,ISTRA)=SPUMPI(ISPZ,ISTRA)+SPUMP(ISPZ,J)
652     CONTINUE
651   CONTINUE
C
C  SUM OVER SPECIES INDEX FOR INTEGRATED VOLUME AVERAGED TALLIES
C                         AND INTEGRATED SURFACE AVERAGED TALLIES
C
      DO 661 IPLS=1,NPLSI
        PAPLI(0,ISTRA)=PAPLI(0,ISTRA)+PAPLI(IPLS,ISTRA)
        PMPLI(0,ISTRA)=PMPLI(0,ISTRA)+PMPLI(IPLS,ISTRA)
        PIPLI(0,ISTRA)=PIPLI(0,ISTRA)+PIPLI(IPLS,ISTRA)
        PPHPLI(0,ISTRA)=PPHPLI(0,ISTRA)+PPHPLI(IPLS,ISTRA)
        PPPLI(0,ISTRA)=PPPLI(0,ISTRA)+PPPLI(IPLS,ISTRA)
        POTPLI(0,ISTRA)=POTPLI(0,ISTRA)+POTPLI(IPLS,ISTRA)
        EOTPLI(0,ISTRA)=EOTPLI(0,ISTRA)+EOTPLI(IPLS,ISTRA)
        SPTPLI(0,ISTRA)=SPTPLI(0,ISTRA)+SPTPLI(IPLS,ISTRA)
661   CONTINUE
      DO 662 IION=1,NIONI
        PDENII(0,ISTRA)=PDENII(0,ISTRA)+PDENII(IION,ISTRA)
        EDENII(0,ISTRA)=EDENII(0,ISTRA)+EDENII(IION,ISTRA)
        PPIOI (0,ISTRA)=PPIOI (0,ISTRA)+PPIOI (IION,ISTRA)
        PAIOI (0,ISTRA)=PAIOI (0,ISTRA)+PAIOI (IION,ISTRA)
        PMIOI (0,ISTRA)=PMIOI (0,ISTRA)+PMIOI (IION,ISTRA)
        PIIOI (0,ISTRA)=PIIOI (0,ISTRA)+PIIOI (IION,ISTRA)
        PPHIOI(0,ISTRA)=PPHIOI(0,ISTRA)+PPHIOI(IION,ISTRA)
        POTIOI(0,ISTRA)=POTIOI(0,ISTRA)+POTIOI(IION,ISTRA)
        PRFAII(0,ISTRA)=PRFAII(0,ISTRA)+PRFAII(IION,ISTRA)
        PRFMII(0,ISTRA)=PRFMII(0,ISTRA)+PRFMII(IION,ISTRA)
        PRFIII(0,ISTRA)=PRFIII(0,ISTRA)+PRFIII(IION,ISTRA)
        PRFPHII(0,ISTRA)=PRFPHII(0,ISTRA)+PRFPHII(IION,ISTRA)
        PRFPII(0,ISTRA)=PRFPII(0,ISTRA)+PRFPII(IION,ISTRA)
        EOTIOI(0,ISTRA)=EOTIOI(0,ISTRA)+EOTIOI(IION,ISTRA)
        ERFAII(0,ISTRA)=ERFAII(0,ISTRA)+ERFAII(IION,ISTRA)
        ERFMII(0,ISTRA)=ERFMII(0,ISTRA)+ERFMII(IION,ISTRA)
        ERFIII(0,ISTRA)=ERFIII(0,ISTRA)+ERFIII(IION,ISTRA)
        ERFPHII(0,ISTRA)=ERFPHII(0,ISTRA)+ERFPHII(IION,ISTRA)
        ERFPII(0,ISTRA)=ERFPII(0,ISTRA)+ERFPII(IION,ISTRA)
        SPTIOI(0,ISTRA)=SPTIOI(0,ISTRA)+SPTIOI(IION,ISTRA)
        PGENII(0,ISTRA)=PGENII(0,ISTRA)+PGENII(IION,ISTRA)
        EGENII(0,ISTRA)=EGENII(0,ISTRA)+EGENII(IION,ISTRA)
        VGENII(0,ISTRA)=VGENII(0,ISTRA)+VGENII(IION,ISTRA)
        EELFI (0,ISTRA)=EELFI (0,ISTRA)+EELFI (IION,ISTRA)
662   CONTINUE
      DO 663 IMOL=1,NMOLI
        PDENMI(0,ISTRA)=PDENMI(0,ISTRA)+PDENMI(IMOL,ISTRA)
        EDENMI(0,ISTRA)=EDENMI(0,ISTRA)+EDENMI(IMOL,ISTRA)
        PPMLI (0,ISTRA)=PPMLI (0,ISTRA)+PPMLI (IMOL,ISTRA)
        PAMLI (0,ISTRA)=PAMLI (0,ISTRA)+PAMLI (IMOL,ISTRA)
        PMMLI (0,ISTRA)=PMMLI (0,ISTRA)+PMMLI (IMOL,ISTRA)
        PIMLI (0,ISTRA)=PIMLI (0,ISTRA)+PIMLI (IMOL,ISTRA)
        PPHMLI(0,ISTRA)=PPHMLI(0,ISTRA)+PPHMLI(IMOL,ISTRA)
        POTMLI(0,ISTRA)=POTMLI(0,ISTRA)+POTMLI(IMOL,ISTRA)
        PRFAMI(0,ISTRA)=PRFAMI(0,ISTRA)+PRFAMI(IMOL,ISTRA)
        PRFMMI(0,ISTRA)=PRFMMI(0,ISTRA)+PRFMMI(IMOL,ISTRA)
        PRFIMI(0,ISTRA)=PRFIMI(0,ISTRA)+PRFIMI(IMOL,ISTRA)
        PRFPHMI(0,ISTRA)=PRFPHMI(0,ISTRA)+PRFPHMI(IMOL,ISTRA)
        PRFPMI(0,ISTRA)=PRFPMI(0,ISTRA)+PRFPMI(IMOL,ISTRA)
        EOTMLI(0,ISTRA)=EOTMLI(0,ISTRA)+EOTMLI(IMOL,ISTRA)
        ERFAMI(0,ISTRA)=ERFAMI(0,ISTRA)+ERFAMI(IMOL,ISTRA)
        ERFMMI(0,ISTRA)=ERFMMI(0,ISTRA)+ERFMMI(IMOL,ISTRA)
        ERFIMI(0,ISTRA)=ERFIMI(0,ISTRA)+ERFIMI(IMOL,ISTRA)
        ERFPHMI(0,ISTRA)=ERFPHMI(0,ISTRA)+ERFPHMI(IMOL,ISTRA)
        ERFPMI(0,ISTRA)=ERFPMI(0,ISTRA)+ERFPMI(IMOL,ISTRA)
        SPTMLI(0,ISTRA)=SPTMLI(0,ISTRA)+SPTMLI(IMOL,ISTRA)
        PGENMI(0,ISTRA)=PGENMI(0,ISTRA)+PGENMI(IMOL,ISTRA)
        EGENMI(0,ISTRA)=EGENMI(0,ISTRA)+EGENMI(IMOL,ISTRA)
        VGENMI(0,ISTRA)=VGENMI(0,ISTRA)+VGENMI(IMOL,ISTRA)
663   CONTINUE
      DO 664 IATM=1,NATMI
        PDENAI(0,ISTRA)=PDENAI(0,ISTRA)+PDENAI(IATM,ISTRA)
        EDENAI(0,ISTRA)=EDENAI(0,ISTRA)+EDENAI(IATM,ISTRA)
        PPATI (0,ISTRA)=PPATI (0,ISTRA)+PPATI (IATM,ISTRA)
        PAATI (0,ISTRA)=PAATI (0,ISTRA)+PAATI (IATM,ISTRA)
        PMATI (0,ISTRA)=PMATI (0,ISTRA)+PMATI (IATM,ISTRA)
        PIATI (0,ISTRA)=PIATI (0,ISTRA)+PIATI (IATM,ISTRA)
        PPHATI(0,ISTRA)=PPHATI(0,ISTRA)+PPHATI(IATM,ISTRA)
        POTATI(0,ISTRA)=POTATI(0,ISTRA)+POTATI(IATM,ISTRA)
        PRFAAI(0,ISTRA)=PRFAAI(0,ISTRA)+PRFAAI(IATM,ISTRA)
        PRFMAI(0,ISTRA)=PRFMAI(0,ISTRA)+PRFMAI(IATM,ISTRA)
        PRFIAI(0,ISTRA)=PRFIAI(0,ISTRA)+PRFIAI(IATM,ISTRA)
        PRFPHAI(0,ISTRA)=PRFPHAI(0,ISTRA)+PRFPHAI(IATM,ISTRA)
        PRFPAI(0,ISTRA)=PRFPAI(0,ISTRA)+PRFPAI(IATM,ISTRA)
        EOTATI(0,ISTRA)=EOTATI(0,ISTRA)+EOTATI(IATM,ISTRA)
        ERFAAI(0,ISTRA)=ERFAAI(0,ISTRA)+ERFAAI(IATM,ISTRA)
        ERFMAI(0,ISTRA)=ERFMAI(0,ISTRA)+ERFMAI(IATM,ISTRA)
        ERFIAI(0,ISTRA)=ERFIAI(0,ISTRA)+ERFIAI(IATM,ISTRA)
        ERFPHAI(0,ISTRA)=ERFPHAI(0,ISTRA)+ERFPHAI(IATM,ISTRA)
        ERFPAI(0,ISTRA)=ERFPAI(0,ISTRA)+ERFPAI(IATM,ISTRA)
        SPTATI(0,ISTRA)=SPTATI(0,ISTRA)+SPTATI(IATM,ISTRA)
        PGENAI(0,ISTRA)=PGENAI(0,ISTRA)+PGENAI(IATM,ISTRA)
        EGENAI(0,ISTRA)=EGENAI(0,ISTRA)+EGENAI(IATM,ISTRA)
        VGENAI(0,ISTRA)=VGENAI(0,ISTRA)+VGENAI(IATM,ISTRA)
664   CONTINUE
      DO IPHOT=1,NPHOTI
        PDENPHI (0,ISTRA)=PDENPHI (0,ISTRA)+PDENPHI (IPHOT,ISTRA)
        EDENPHI (0,ISTRA)=EDENPHI (0,ISTRA)+EDENPHI (IPHOT,ISTRA)
        PPPHTI  (0,ISTRA)=PPPHTI  (0,ISTRA)+PPPHTI  (IPHOT,ISTRA)
        PAPHTI  (0,ISTRA)=PAPHTI  (0,ISTRA)+PAPHTI  (IPHOT,ISTRA)
        PMPHTI  (0,ISTRA)=PMPHTI  (0,ISTRA)+PMPHTI  (IPHOT,ISTRA)
        PIPHTI  (0,ISTRA)=PIPHTI  (0,ISTRA)+PIPHTI  (IPHOT,ISTRA)
        PPHPHTI (0,ISTRA)=PPHPHTI (0,ISTRA)+PPHPHTI (IPHOT,ISTRA)
        POTPHTI (0,ISTRA)=POTPHTI (0,ISTRA)+POTPHTI (IPHOT,ISTRA)
        PRFAPHTI(0,ISTRA)=PRFAPHTI(0,ISTRA)+PRFAPHTI(IPHOT,ISTRA)
        PRFMPHTI(0,ISTRA)=PRFMPHTI(0,ISTRA)+PRFMPHTI(IPHOT,ISTRA)
        PRFIPHTI(0,ISTRA)=PRFIPHTI(0,ISTRA)+PRFIPHTI(IPHOT,ISTRA)
        PRFPHPHTI(0,ISTRA)=PRFPHPHTI(0,ISTRA)+PRFPHPHTI(IPHOT,ISTRA)
        PRFPPHTI(0,ISTRA)=PRFPPHTI(0,ISTRA)+PRFPPHTI(IPHOT,ISTRA)
        EOTPHTI (0,ISTRA)=EOTPHTI (0,ISTRA)+EOTPHTI (IPHOT,ISTRA)
        ERFAPHTI(0,ISTRA)=ERFAPHTI(0,ISTRA)+ERFAPHTI(IPHOT,ISTRA)
        ERFMPHTI(0,ISTRA)=ERFMPHTI(0,ISTRA)+ERFMPHTI(IPHOT,ISTRA)
        ERFIPHTI(0,ISTRA)=ERFIPHTI(0,ISTRA)+ERFIPHTI(IPHOT,ISTRA)
        ERFPHPHTI(0,ISTRA)=ERFPHPHTI(0,ISTRA)+ERFPHPHTI(IPHOT,ISTRA)
        ERFPPHTI(0,ISTRA)=ERFPPHTI(0,ISTRA)+ERFPPHTI(IPHOT,ISTRA)
        SPTPHTI (0,ISTRA)=SPTPHTI (0,ISTRA)+SPTPHTI (IPHOT,ISTRA)
        PGENPHI (0,ISTRA)=PGENPHI (0,ISTRA)+PGENPHI (IPHOT,ISTRA)
        EGENPHI (0,ISTRA)=EGENPHI (0,ISTRA)+EGENPHI (IPHOT,ISTRA)
        VGENPHI (0,ISTRA)=VGENPHI (0,ISTRA)+VGENPHI (IPHOT,ISTRA)
      END DO
      DO 665 IADV=1,NADVI
        ADDVI(0,ISTRA)=ADDVI(0,ISTRA)+ADDVI(IADV,ISTRA)
665   CONTINUE
      DO 666 ICLV=1,NCLVI
        COLVI(0,ISTRA)=COLVI(0,ISTRA)+COLVI(ICLV,ISTRA)
666   CONTINUE
      DO 667 ISNV=1,NSNVI
        SNAPVI(0,ISTRA)=SNAPVI(0,ISTRA)+SNAPVI(ISNV,ISTRA)
667   CONTINUE
      DO 668 ICPV=1,NCPVI
        COPVI(0,ISTRA)=COPVI(0,ISTRA)+COPVI(ICPV,ISTRA)
668   CONTINUE
      DO 669 IBGV=1,NBGVI
        BGKVI(0,ISTRA)=BGKVI(0,ISTRA)+BGKVI(IBGV,ISTRA)
669   CONTINUE
C
      RETURN
      END SUBROUTINE SUM_AVERAGE