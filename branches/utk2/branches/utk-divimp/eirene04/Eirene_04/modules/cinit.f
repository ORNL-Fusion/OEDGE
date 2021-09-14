      MODULE CINIT

      USE PRECISION
      USE PARMMOD

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: ALLOC_CINIT, DEALLOC_CINIT, INIT_CINIT,
     .          TDENMODEL, TDENMODAR

      INTEGER, PUBLIC, SAVE ::
     I         NCINIT, MCINIT, LCINIT

      REAL(DP),  PUBLIC, TARGET, ALLOCATABLE, SAVE :: RCINIT(:)
      INTEGER, PUBLIC, TARGET, ALLOCATABLE, SAVE :: ICINIT(:)
      LOGICAL, PUBLIC, TARGET, SAVE :: LCNIT(2)

      REAL(DP), PUBLIC, POINTER, SAVE ::
     R EP1IN,  EP1CH,  EP1OT,  EXEP1,
     R ELLIN,  ELLOT,  ELLCH,  EXELL,
     R TRIIN,  TRIOT,  TRICH,  EXTRI,
     R TE0,    TE1,    TE2,    TE3,    TE4,   TE5,
     R BP0,    BP1,    BP2,    BP3,    BP4,   BP5,
     R B0,     B1,     B2,     B3,     B4,    B5

      REAL(DP), PUBLIC, POINTER, SAVE ::
     R TI0(:),TI1(:),TI2(:),TI3(:),TI4(:),TI5(:),
     R VX0(:),VX1(:),VX2(:),VX3(:),VX4(:),VX5(:),
     R VY0(:),VY1(:),VY2(:),VY3(:),VY4(:),VY5(:),
     R VZ0(:),VZ1(:),VZ2(:),VZ3(:),VZ4(:),VZ5(:),
     R DI0(:),DI1(:),DI2(:),DI3(:),DI4(:),DI5(:)

      INTEGER, PUBLIC, POINTER, SAVE ::
     I INDPRO(:),INDGRD(:),INDSRC(:)

      LOGICAL, PUBLIC, POINTER, SAVE ::
     L NLMACH,NLMLTI

      CHARACTER(66), PUBLIC, SAVE :: CASENAME
      CHARACTER(10), PUBLIC, ALLOCATABLE, SAVE :: CDENMODEL(:) 

      TYPE TDENMODEL
        REAL(DP) :: G_BOLTZ, DELTAE, A_CORONA
        INTEGER :: NRE
        INTEGER, POINTER :: ISP(:)
        CHARACTER(8), POINTER :: FNAME(:)
        CHARACTER(4), POINTER :: H2(:)
        CHARACTER(9), POINTER :: REACTION(:)
        CHARACTER(3), POINTER :: CR(:)
      END TYPE

      TYPE TDENMODAR
        TYPE(TDENMODEL), POINTER :: TDM
      END TYPE
      
      TYPE(TDENMODAR), PUBLIC, ALLOCATABLE, SAVE :: TDMPAR(:) 

      CONTAINS

      SUBROUTINE ALLOC_CINIT

      IF (ALLOCATED(RCINIT)) RETURN

      NCINIT=30+30*NPLS
      MCINIT=12+3+NSTRA
      LCINIT=2

      ALLOCATE (RCINIT(NCINIT))
      ALLOCATE (ICINIT(MCINIT))
      ALLOCATE (CDENMODEL(NPLS))
      ALLOCATE (TDMPAR(NPLS))

      WRITE (55,'(A,T25,I15)')
     .      ' CINIT ',NCINIT*8 + MCINIT*4 +
     .                NPLS*LEN(CDENMODEL(1))

      EP1IN  => RCINIT(1)
      EP1CH  => RCINIT(2)
      EP1OT  => RCINIT(3)
      EXEP1  => RCINIT(4)
      ELLIN  => RCINIT(5)
      ELLOT  => RCINIT(6)
      ELLCH  => RCINIT(7)
      EXELL  => RCINIT(8)
      TRIIN  => RCINIT(9)
      TRIOT  => RCINIT(10)
      TRICH  => RCINIT(11)
      EXTRI  => RCINIT(12)
      TE0    => RCINIT(13)
      TE1    => RCINIT(14)
      TE2    => RCINIT(15)
      TE3    => RCINIT(16)
      TE4    => RCINIT(17)
      TE5    => RCINIT(18)
      BP0    => RCINIT(19)
      BP1    => RCINIT(20)
      BP2    => RCINIT(21)
      BP3    => RCINIT(22)
      BP4    => RCINIT(23)
      BP5    => RCINIT(24)
      B0     => RCINIT(25)
      B1     => RCINIT(26)
      B2     => RCINIT(27)
      B3     => RCINIT(28)
      B4     => RCINIT(29)
      B5     => RCINIT(30)
      TI0    => RCINIT(31         : 30+ 1*NPLS)
      TI1    => RCINIT(31+ 1*npls : 30+ 2*NPLS)
      TI2    => RCINIT(31+ 2*npls : 30+ 3*NPLS)
      TI3    => RCINIT(31+ 3*npls : 30+ 4*NPLS)
      TI4    => RCINIT(31+ 4*npls : 30+ 5*NPLS)
      TI5    => RCINIT(31+ 5*npls : 30+ 6*NPLS)
      VX0    => RCINIT(31+ 6*npls : 30+ 7*NPLS)
      VX1    => RCINIT(31+ 7*npls : 30+ 8*NPLS)
      VX2    => RCINIT(31+ 8*npls : 30+ 9*NPLS)
      VX3    => RCINIT(31+ 9*npls : 30+10*NPLS)
      VX4    => RCINIT(31+10*npls : 30+11*NPLS)
      VX5    => RCINIT(31+11*npls : 30+12*NPLS)
      VY0    => RCINIT(31+12*npls : 30+13*NPLS)
      VY1    => RCINIT(31+13*npls : 30+14*NPLS)
      VY2    => RCINIT(31+14*npls : 30+15*NPLS)
      VY3    => RCINIT(31+15*npls : 30+16*NPLS)
      VY4    => RCINIT(31+16*npls : 30+17*NPLS)
      VY5    => RCINIT(31+17*npls : 30+18*NPLS)
      VZ0    => RCINIT(31+18*npls : 30+19*NPLS)
      VZ1    => RCINIT(31+19*npls : 30+20*NPLS)
      VZ2    => RCINIT(31+20*npls : 30+21*NPLS)
      VZ3    => RCINIT(31+21*npls : 30+22*NPLS)
      VZ4    => RCINIT(31+22*npls : 30+23*NPLS)
      VZ5    => RCINIT(31+23*npls : 30+24*NPLS)
      DI0    => RCINIT(31+24*npls : 30+25*NPLS)
      DI1    => RCINIT(31+25*npls : 30+26*NPLS)
      DI2    => RCINIT(31+26*npls : 30+27*NPLS)
      DI3    => RCINIT(31+27*npls : 30+28*NPLS)
      DI4    => RCINIT(31+28*npls : 30+29*NPLS)
      DI5    => RCINIT(31+29*npls : 30+30*NPLS)

      INDPRO => ICINIT( 1 : 12)
      INDGRD => ICINIT(13 : 15)
      INDSRC => ICINIT(16 : 15+NSTRA)

      NLMACH => LCNIT(1)
      NLMLTI => LCNIT(2)

      CALL INIT_CINIT

      RETURN

      END SUBROUTINE ALLOC_CINIT


      SUBROUTINE DEALLOC_CINIT

      IF (.NOT.ALLOCATED(RCINIT)) RETURN

      DEALLOCATE (RCINIT)
      DEALLOCATE (ICINIT)
      DEALLOCATE (CDENMODEL)

      RETURN
      END SUBROUTINE DEALLOC_CINIT


      SUBROUTINE INIT_CINIT

      RCINIT = 0._DP
      ICINIT = 0
      LCNIT  = .FALSE.
!pb      CDENMODEL = REPEAT(' ',LEN(CDENMODEL(1)))
      CDENMODEL = REPEAT(' ',LEN(CDENMODEL))
C
C  INITIALIZE GEOMETRY DATA
C
      ELLIN=1._DP
      ELLOT=1._DP
      ELLCH=1._DP
      TRIIN=1._DP
      TRIOT=1._DP
      TRICH=1._DP

      RETURN
      END SUBROUTINE INIT_CINIT

      END MODULE CINIT
