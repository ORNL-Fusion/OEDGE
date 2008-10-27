      MODULE CADGEO

      USE PRECISION
      USE PARMMOD

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: ALLOC_CADGEO, DEALLOC_CADGEO, INIT_CADGEO

      INTEGER, PUBLIC, SAVE ::
     I NPLIM,  NADGEO, MADGEO

      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     R        RADGEO(:,:)

      INTEGER, PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     I         IADGEO(:,:)

C NADGEO, REAL
      REAL(DP), PUBLIC, POINTER, SAVE ::
     R A0LM(:),     A1LM(:),     A2LM(:),     A3LM(:),     A4LM(:),
     R A5LM(:),     A6LM(:),     A7LM(:),     A8LM(:),     A9LM(:),
     R ALIMS(:,:),  XLIMS(:,:),  YLIMS(:,:),  ZLIMS(:,:),
     R ALIMS0(:,:), XLIMS1(:,:), YLIMS1(:,:), ZLIMS1(:,:),
     R XLIMS2(:,:), YLIMS2(:,:), ZLIMS2(:,:),
     R XLIMS3(:,:), YLIMS3(:,:), ZLIMS3(:,:),
     R RLB(:),      
     R P1(:,:),     P2(:,:),     P3(:,:),     P4(:,:),    P5(:,:),
     R P6(:,:),
     R ALM(:),      BLM(:),      CLM(:),
     R PS13(:,:),   PS23(:,:),
     R PS24(:,:),   PS34(:,:),   PS35(:,:),   PS45(:,:),
     R P1A(:),      P2A(:),      P1B(:),      P2B(:),
     R P1C(:),      P2C(:)

C MADGEO, INTEGER
      INTEGER, PUBLIC, POINTER, SAVE ::
     I ILIN(:),     ISCN(:)

      INTEGER, PUBLIC, SAVE :: NLIMI

      LOGICAL, PUBLIC, ALLOCATABLE, SAVE ::
     L RLBNOT(:)

      CONTAINS


      SUBROUTINE ALLOC_CADGEO

      IF (ALLOCATED(RADGEO)) RETURN

      NPLIM=10+14*9+10+12*3
      NADGEO=NLIM*NPLIM
      MADGEO=NLIM*2

      ALLOCATE (RADGEO(NPLIM,NLIM))
      ALLOCATE (IADGEO(2,NLIM))
      ALLOCATE (RLBNOT(NLIM))

      WRITE (55,*) ' CADGEO ',NLIM*(NPLIM+1)*8 + 2*NLIM*4

      A0LM => RADGEO(1,:)
      A1LM => RADGEO(2,:)
      A2LM => RADGEO(3,:)
      A3LM => RADGEO(4,:)
      A4LM => RADGEO(5,:)
      A5LM => RADGEO(6,:)
      A6LM => RADGEO(7,:)
      A7LM => RADGEO(8,:)
      A8LM => RADGEO(9,:)
      A9LM => RADGEO(10,:)
      ALIMS => RADGEO(11:19,:)
      XLIMS => RADGEO(20:28,:)
      YLIMS => RADGEO(29:37,:)
      ZLIMS => RADGEO(38:46,:)
      ALIMS0 => RADGEO(47:55,:)
      XLIMS1 => RADGEO(56:64,:)
      YLIMS1 => RADGEO(65:73,:)
      ZLIMS1 => RADGEO(74:82,:)
      XLIMS2 => RADGEO(83:91,:)
      YLIMS2 => RADGEO(92:100,:)
      ZLIMS2 => RADGEO(101:109,:)
      XLIMS3 => RADGEO(110:118,:)
      YLIMS3 => RADGEO(119:127,:)
      ZLIMS3 => RADGEO(128:136,:)
      RLB => RADGEO(137,:)
      P1 => RADGEO(138:140,:)
      P2 => RADGEO(141:143,:)
      P3 => RADGEO(144:146,:)
      P4 => RADGEO(147:149,:)
      P5 => RADGEO(150:152,:)
      P6 => RADGEO(153:155,:)
      ALM => RADGEO(156,:)
      BLM => RADGEO(157,:)
      CLM => RADGEO(158,:)
      PS13 => RADGEO(159:161,:)
      PS23 => RADGEO(162:164,:)
      PS24 => RADGEO(165:167,:)
      PS34 => RADGEO(168:170,:)
      PS35 => RADGEO(171:173,:)
      PS45 => RADGEO(174:176,:)
      P1A => RADGEO(177,:)
      P2A => RADGEO(178,:)
      P1B => RADGEO(179,:)
      P2B => RADGEO(180,:)
      P1C => RADGEO(181,:)
      P2C => RADGEO(182,:)

      ILIN => IADGEO(1,:)
      ISCN => IADGEO(2,:)

      CALL INIT_CADGEO

      RETURN
      END SUBROUTINE ALLOC_CADGEO

      SUBROUTINE DEALLOC_CADGEO

      IF (.NOT.ALLOCATED(RADGEO)) RETURN

      DEALLOCATE (RADGEO)
      DEALLOCATE (IADGEO)
      DEALLOCATE (RLBNOT)

      RETURN
      END SUBROUTINE DEALLOC_CADGEO
      

      SUBROUTINE INIT_CADGEO

      RADGEO = 0.D0
      IADGEO = 0
      RLBNOT = .FALSE.

      RETURN
      END SUBROUTINE INIT_CADGEO

      END MODULE CADGEO


