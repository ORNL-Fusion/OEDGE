!pb  30.11.06: ltstv added for testing purposes
!    20.06.07: constant NLOGAU = number of logicals introduced

      MODULE CLOGAU

      USE PRECISION

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: ALLOC_CLOGAU

      INTEGER, PUBLIC, PARAMETER :: NLOGAU=31

      LOGICAL, PUBLIC, TARGET, SAVE :: LLOGAU(NLOGAU)

      LOGICAL, PUBLIC, POINTER, SAVE ::
     L NLSCL,  NLDRFT, NLCRR,  NLTEST, NLANA,  NLERG,  NLMOVIE, NLPLAS,
     L NLRAD,  NLSLB,  NLCRC,  NLELL,  NLTRI,  NLPLG,  NLGEN,
     L NLFEM,  NLTET,
     L NLPOL,  NLPLY,  NLPLA,  NLPLP,
     L NLTOR,  NLTRZ,  NLTRA,  NLTRT,
     L NLMLT,  NLADD,
     L NLTRIM, LHABER, NLONE,  LTSTV

      CONTAINS


      SUBROUTINE ALLOC_CLOGAU

      NLSCL   => LLOGAU(1)
      NLDRFT  => LLOGAU(2)
      NLCRR   => LLOGAU(3)
      NLTEST  => LLOGAU(4)
      NLANA   => LLOGAU(5)
      NLERG   => LLOGAU(6)
      NLMOVIE => LLOGAU(7)
      NLPLAS  => LLOGAU(8)
      NLRAD   => LLOGAU(9)
      NLSLB   => LLOGAU(10)
      NLCRC   => LLOGAU(11)
      NLELL   => LLOGAU(12)
      NLTRI   => LLOGAU(13)
      NLPLG   => LLOGAU(14)
      NLGEN   => LLOGAU(15)
      NLFEM   => LLOGAU(16)
      NLTET   => LLOGAU(17)
      NLPOL   => LLOGAU(18)
      NLPLY   => LLOGAU(19)
      NLPLA   => LLOGAU(20)
      NLPLP   => LLOGAU(21)
      NLTOR   => LLOGAU(22)
      NLTRZ   => LLOGAU(23)
      NLTRA   => LLOGAU(24)
      NLTRT   => LLOGAU(25)
      NLMLT   => LLOGAU(26)
      NLADD   => LLOGAU(27)
      NLTRIM  => LLOGAU(28)
      LHABER  => LLOGAU(29)
      NLONE   => LLOGAU(30)
      LTSTV   => LLOGAU(31)

      LLOGAU = .FALSE.

      END SUBROUTINE ALLOC_CLOGAU

      END MODULE CLOGAU
