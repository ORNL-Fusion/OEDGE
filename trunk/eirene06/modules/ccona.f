C   6.12.05   AU_TO_CM2 added here (and removed from fpatha, veloel)
C   1.01.06   hplnk_bar = hplnk/2Pi added here (and set in setcon.f)
      MODULE CCONA

      USE PRECISION

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: ALLOC_CCONA

      REAL(DP), PUBLIC, TARGET, SAVE :: RCONA(43)

      REAL(DP), PUBLIC, POINTER, SAVE ::
     R EPS60,   EPS30,   EPS12,    EPS10,    EPS6,      EPS5,
     R PMASSA,  PMASSE,  AMUA,
     R PIA,     PI2A,    PIHA,     PIQU,     PISQ,     PIAI,      SQ2,
     R SQ2I,    DEGRAD,  RADDEG,
     R CVELI2,  CVELAA,  CVEL2A,   ELCHA,    EFACT,     EFCT23,   EVKEL,
     R EIONH,   EIONH2,  EIONHE,
     R EV_TO_J, J_TO_EV, J_TO_ERG, ERG_TO_J, EV_TO_ERG, ERG_TO_EV,
     R HPLANCK, CLIGHT, MUB, HPLNK, HPCL, EV2HZ,
     R AU_TO_CM2, HPLNK_BAR

      CONTAINS

      SUBROUTINE ALLOC_CCONA

      EPS60     => RCONA(1)
      EPS30     => RCONA(2)
      EPS12     => RCONA(3)
      EPS10     => RCONA(4)
      EPS6      => RCONA(5)
      EPS5      => RCONA(6)
      PMASSA    => RCONA(7)
      PMASSE    => RCONA(8)
      AMUA      => RCONA(9)
      PIA       => RCONA(10)
      PI2A      => RCONA(11)
      PIHA      => RCONA(12)
      PISQ      => RCONA(13)
      PIAI      => RCONA(14)
      SQ2       => RCONA(15)
      SQ2I      => RCONA(16)
      DEGRAD    => RCONA(17)
      RADDEG    => RCONA(18)
      CVELI2    => RCONA(19)
      CVELAA    => RCONA(20)
      CVEL2A    => RCONA(21)
      ELCHA     => RCONA(22)
      EFACT     => RCONA(23)
      EFCT23    => RCONA(24)
      EVKEL     => RCONA(25)
      EIONH     => RCONA(26)
      EIONH2    => RCONA(27)
      EIONHE    => RCONA(28)
      EV_TO_J   => RCONA(29)
      J_TO_EV   => RCONA(30)
      J_TO_ERG  => RCONA(31)
      ERG_TO_J  => RCONA(32)
      EV_TO_ERG => RCONA(33)
      ERG_TO_EV => RCONA(34)
      HPLANCK   => RCONA(35)
      CLIGHT    => RCONA(36)
      MUB       => RCONA(37)
      PIQU      => RCONA(38)
      HPLNK     => RCONA(39)
      HPCL      => RCONA(40)
      EV2HZ     => RCONA(41)
      AU_TO_CM2 => RCONA(42)
      HPLNK_BAR => RCONA(43)

      RETURN
      END SUBROUTINE ALLOC_CCONA

      END MODULE CCONA
