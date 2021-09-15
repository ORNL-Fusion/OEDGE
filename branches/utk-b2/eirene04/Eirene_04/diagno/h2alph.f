C
C
      FUNCTION H2ALPH(TE,LREACT)
C
C  DEGAS MODEL FOR H ALPHA PRODUCTION VIA H2 DISSOCIATION
C  INPUT: TE IN EV
C
      USE PRECISION

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: TE
      INTEGER, INTENT(IN) :: LREACT

      INTEGER :: ITEMPS(2)
      REAL(DP) :: ZA(6,2), ZB(6,2), ZTEMPS(6,2)
      REAL(DP) :: H2ALPH, ZTE
      INTEGER :: IT, IINDEX, I
      SAVE
      DATA ITEMPS/6,4/
C                          LREACT=1: E + H20 -> 2H0
      DATA (ZTEMPS(I,1),I=1,6)/
     .      2.7,4.8,10.0,25.0,100.0,1.0E10/
      DATA (ZA(I,1),I=1,6)/
     .      5.1E-5,2.57E-4,5.52E-4,1.7E-3,3.57E-3,1.06E-2/
      DATA (ZB(I,1),I=1,6)/
     .      2.994,1.368,0.881,0.389,0.161,-0.076/
C                          LREACT=2: E + H2+ -> H0 + H+ +E
      DATA (ZTEMPS(I,2),I=1,4)/
     .      3.5,6.3,30.0,1.0E10/
      DATA (ZA(I,2),I=1,4)/
     .      2.1E-2,1.07E-2,3.2E-3,1.2E-2/
      DATA (ZB(I,2),I=1,4)/
     .     -0.822,-0.243,0.392,1.0E-10/
C
      ZTE=MAX(1._DP,TE)
      ZTE=MIN(ZTE,1.0E10_DP)
      IINDEX=0
      DO 15 IT=1,ITEMPS(LREACT)
        IINDEX=IINDEX+1
        IF (ZTE.LE.ZTEMPS(IT,LREACT)) GOTO 20
15    CONTINUE
C
20    H2ALPH=ZA(IINDEX,LREACT)*ZTE**ZB(IINDEX,LREACT)
      RETURN
      END
