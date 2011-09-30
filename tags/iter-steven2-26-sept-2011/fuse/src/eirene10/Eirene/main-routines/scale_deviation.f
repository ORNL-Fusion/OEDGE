      SUBROUTINE SCALE_DEVIATION (ISTR, ZWW, ZW, ZVOLNT, ZVOLWT,
     .                     ZVOLIN, ZVOLIW, SCLTAL, N1DIM)

      USE PRECISION
      USE PARMMOD
      USE CSDVI
      USE CGRID
      USE COUTAU

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ISTR, N1DIM
      REAL(DP), INTENT(IN) :: ZWW, ZW, ZVOLNT, ZVOLWT
      REAL(DP), INTENT(IN) :: SCLTAL(N1DIM,*)
      REAL(DP), INTENT(IN) :: ZVOLIN(*), ZVOLIW(*)

      INTEGER :: ISDV, ITAL, ICELL, ISPZ

      DO 900 ISDV=1,NSIGCI
        DO 901 ITAL=1,NTALR
          IF (IIHC(1,ISDV).NE.ITAL) GOTO 901
          ISPZ=MAX(1,IGHC(1,ISDV))
          IF (SCLTAL(ISPZ,ITAL).EQ.1) THEN
            DO 902 ICELL=1,NSBOX_TAL
              SIGMAC(1,ISDV,ICELL)=SIGMAC(1,ISDV,ICELL)*ZVOLIN(ICELL)*
     .                             XMCP(ISTR)
              SIGMAC(0,ISDV,ICELL)=SIGMAC(0,ISDV,ICELL)*ZVOLIN(ICELL)*
     .                             XMCP(ISTR)
902         CONTINUE
            SGMCS(1,ISDV)=SGMCS(1,ISDV)*ZVOLNT*XMCP(ISTR)
            SGMCS(0,ISDV)=SGMCS(0,ISDV)*ZVOLNT*XMCP(ISTR)
          ELSEIF (SCLTAL(ISPZ,ITAL).EQ.2) THEN
            DO 903 ICELL=1,NSBOX_TAL
              SIGMAC(1,ISDV,ICELL)=SIGMAC(1,ISDV,ICELL)*ZW*
     .                             XMCP(ISTR)
              SIGMAC(0,ISDV,ICELL)=SIGMAC(0,ISDV,ICELL)*ZW*
     .                             XMCP(ISTR)
903         CONTINUE
            SGMCS(1,ISDV)=SGMCS(1,ISDV)*ZW*XMCP(ISTR)
            SGMCS(0,ISDV)=SGMCS(0,ISDV)*ZW*XMCP(ISTR)
          ELSEIF (SCLTAL(ISPZ,ITAL).EQ.3) THEN
            DO 904 ICELL=1,NSBOX_TAL
              SIGMAC(1,ISDV,ICELL)=SIGMAC(1,ISDV,ICELL)*ZVOLIW(ICELL)*
     .                             XMCP(ISTR)
              SIGMAC(0,ISDV,ICELL)=SIGMAC(0,ISDV,ICELL)*ZVOLIW(ICELL)*
     .                             XMCP(ISTR)
904         CONTINUE
            SGMCS(1,ISDV)=SGMCS(1,ISDV)*ZVOLWT*XMCP(ISTR)
            SGMCS(0,ISDV)=SGMCS(0,ISDV)*ZVOLWT*XMCP(ISTR)
          ELSEIF (SCLTAL(ISPZ,ITAL).EQ.4) THEN
            DO 905 ICELL=1,NSBOX_TAL
              SIGMAC(1,ISDV,ICELL)=SIGMAC(1,ISDV,ICELL)*ZWW*
     .                             XMCP(ISTR)
              SIGMAC(0,ISDV,ICELL)=SIGMAC(0,ISDV,ICELL)*ZWW*
     .                             XMCP(ISTR)
905         CONTINUE
            SGMCS(1,ISDV)=SGMCS(1,ISDV)*ZWW*XMCP(ISTR)
            SGMCS(0,ISDV)=SGMCS(0,ISDV)*ZWW*XMCP(ISTR)
          ENDIF
901     CONTINUE
C
        DO 911 ITAL=1,NTALR
          IF (IIHC(2,ISDV).NE.ITAL) GOTO 911
          ISPZ=MAX(1,IGHC(2,ISDV))
          IF (SCLTAL(ISPZ,ITAL).EQ.1) THEN
            DO 912 ICELL=1,NSBOX_TAL
              SIGMAC(2,ISDV,ICELL)=SIGMAC(2,ISDV,ICELL)*ZVOLIN(ICELL)*
     .                             XMCP(ISTR)
              SIGMAC(0,ISDV,ICELL)=SIGMAC(0,ISDV,ICELL)*ZVOLIN(ICELL)*
     .                             XMCP(ISTR)
912         CONTINUE
            SGMCS(2,ISDV)=SGMCS(2,ISDV)*ZVOLNT*XMCP(ISTR)
            SGMCS(0,ISDV)=SGMCS(0,ISDV)*ZVOLNT*XMCP(ISTR)
          ELSEIF (SCLTAL(ISPZ,ITAL).EQ.2) THEN
            DO 913 ICELL=1,NSBOX_TAL
              SIGMAC(2,ISDV,ICELL)=SIGMAC(2,ISDV,ICELL)*ZW*
     .                             XMCP(ISTR)
              SIGMAC(0,ISDV,ICELL)=SIGMAC(0,ISDV,ICELL)*ZW*
     .                             XMCP(ISTR)
913         CONTINUE
            SGMCS(2,ISDV)=SGMCS(2,ISDV)*ZW*XMCP(ISTR)
            SGMCS(0,ISDV)=SGMCS(0,ISDV)*ZW*XMCP(ISTR)
          ELSEIF (SCLTAL(ISPZ,ITAL).EQ.3) THEN
            DO 914 ICELL=1,NSBOX_TAL
              SIGMAC(2,ISDV,ICELL)=SIGMAC(2,ISDV,ICELL)*ZVOLIW(ICELL)*
     .                             XMCP(ISTR)
              SIGMAC(0,ISDV,ICELL)=SIGMAC(0,ISDV,ICELL)*ZVOLIW(ICELL)*
     .                             XMCP(ISTR)
914         CONTINUE
            SGMCS(2,ISDV)=SGMCS(2,ISDV)*ZVOLWT*XMCP(ISTR)
            SGMCS(0,ISDV)=SGMCS(0,ISDV)*ZVOLWT*XMCP(ISTR)
          ELSEIF (SCLTAL(ISPZ,ITAL).EQ.4) THEN
            DO 915 ICELL=1,NSBOX_TAL
              SIGMAC(2,ISDV,ICELL)=SIGMAC(2,ISDV,ICELL)*ZWW*
     .                             XMCP(ISTR)
              SIGMAC(0,ISDV,ICELL)=SIGMAC(0,ISDV,ICELL)*ZWW*
     .                             XMCP(ISTR)
915         CONTINUE
            SGMCS(2,ISDV)=SGMCS(2,ISDV)*ZWW*XMCP(ISTR)
            SGMCS(0,ISDV)=SGMCS(0,ISDV)*ZWW*XMCP(ISTR)
          ENDIF
911     CONTINUE
900   CONTINUE

      RETURN

      END SUBROUTINE SCALE_DEVIATION
