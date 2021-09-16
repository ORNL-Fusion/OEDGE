C -*-Fortran-*-
C/ MODULE A1
C
         SUBROUTINE AATOM(PMASS,KIMP,KOUT,KFAIL)
         implicit none
         REAL PMASS(1)
         integer kimp,kout,kfail(1)
C
C A.1 INITIALISE ARRAYS FOR NON-CORONAL IMPURITY RADIATION
C
C VERSION 05/12/85 AVMA JET/OXFORD
C
C-----------------------------------------------------------------------
CL                  C1.1     IMPURITY IDENTIFICATION
C VERSION 07/09/84 AVMA JET/OXFORD
       COMMON/COMAT1/
     I   NORDER,   NIMP  ,   NSTATE
       INTEGER
     I   NORDER(3)       ,   NSTATE(3), NIMP
C-----------------------------------------------------------------------
CL                  C1.2     INTERPOLATION VARIABLES
C VERSION 07/09/84 AVMA JET/OXFORD
       COMMON/COMAT2/
     R   TE    ,   TELOG ,   DTLOG ,   TMINLG,   TMAXLG,   DENE  ,
     R   DENLOG,   DNLOG ,   DMINLG,   DMAXLG,
     I   NTE   ,   NNE
       integer nte,nne
       real dtlog,tminlg,tmaxlg,dnlog,dminlg,dmaxlg
       REAL
     R   TE(34)          ,   TELOG(34)       ,   DENE(18)        ,
     R   DENLOG(18)
C-----------------------------------------------------------------------
CL                  C1.3     INTERPOLATED RATE COEFFICIENTS
C VERSION 07/09/84 AVMA JET/OXFORD
       COMMON/COMAT3/
     R   SAL   ,   RAL   ,   DAL
       REAL
     R   SAL(34,28,3)    ,   RAL(18,34,28,3) ,   DAL(18,34,28,3)
C-----------------------------------------------------------------------
CL                  C1.4     INTERPOLATED RADIATION COEFFICIENTS
C VERSION 20/11/85 AVMA JET/OXFORD
       COMMON/COMAT4/
     R   RLINE ,   BREMS ,   RAD   ,   DIE
       REAL
     R   RLINE(34,29,3)  ,   BREMS(34,28,3)  ,   RAD(18,34,28,3),
     R   DIE(18,34,28,3)
C-----------------------------------------------------------------------
CL                  C1.5     INTERPOLATED IONISATION LOSS COEFFICIENTS
C VERSION 20/11/85 AVMA JET/OXFORD
       COMMON/COMAT5/
     R   RION  ,   P
c
c     jdemod - add save statement to help ensure data is retained
c
       save

       REAL
     R   RION(34,28,3)   ,   P(28,3)
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
CL                  C1.5     LOTZ IONISATION
C
       REAL
     .   ZA1(3,1)  ,             ZA3(3,2)  , ZA4(3,4)  , ZA5(3,6)  ,
     .   ZA6(3,7)  , ZA7(3,8)  , ZA8(3,10) , ZA9(3,13) , ZA10(3,14),
     .   ZA11(3,17), ZA12(3,18), ZA13(3,22), ZA14(3,24), ZA15(3,26)
C
       REAL
     .   ZB1(3,1)  ,             ZB3(3,2)  , ZB4(3,4)  , ZB5(3,6)  ,
     .   ZB6(3,7)  , ZB7(3,8)  , ZB8(3,10) , ZB9(3,13) , ZB10(3,14),
     .   ZB11(3,17), ZB12(3,18), ZB13(3,22), ZB14(3,24), ZB15(3,26)
C
       REAL
     .   ZC1(3,1)  ,             ZC3(3,2)  , ZC4(3,4)  , ZC5(3,6)  ,
     .   ZC6(3,7)  , ZC7(3,8)  , ZC8(3,10) , ZC9(3,13) , ZC10(3,14),
     .   ZC11(3,17), ZC12(3,18), ZC13(3,22), ZC14(3,24), ZC15(3,26)
C
       REAL
     .   ZP1(3,1)  ,             ZP3(3,2)  , ZP4(3,4)  , ZP5(3,6)  ,
     .   ZP6(3,7)  , ZP7(3,8)  , ZP8(3,10) , ZP9(3,13) , ZP10(3,14),
     .   ZP11(3,17), ZP12(3,18), ZP13(3,22), ZP14(3,24), ZP15(3,26),
     .   ZP16(3,28)
C
       INTEGER
     .   IQ1(3,1)  ,             IQ3(3,2)  , IQ4(3,4)  , IQ5(3,6)  ,
     .   IQ6(3,7)  , IQ7(3,8)  , IQ8(3,10) , IQ9(3,13) , IQ10(3,14),
     .   IQ11(3,17), IQ12(3,18), IQ13(3,22), IQ14(3,24), IQ15(3,26)
C
       REAL
     .   ZA(3,28,3),ZB(3,28,3),ZC(3,28,3),ZP(3,28,3)
       INTEGER IQ(3,28,3)
c
c       DIMENSION
c     .   ZA1(3,1)  ,             ZA3(3,2)  , ZA4(3,4)  , ZA5(3,6)  ,
c     .   ZA6(3,7)  , ZA7(3,8)  , ZA8(3,10) , ZA9(3,13) , ZA10(3,14),
c     .   ZA11(3,17), ZA12(3,18), ZA13(3,22), ZA14(3,24), ZA15(3,26)
C
c       DIMENSION
c     .   ZB1(3,1)  ,             ZB3(3,2)  , ZB4(3,4)  , ZB5(3,6)  ,
c     .   ZB6(3,7)  , ZB7(3,8)  , ZB8(3,10) , ZB9(3,13) , ZB10(3,14),
c     .   ZB11(3,17), ZB12(3,18), ZB13(3,22), ZB14(3,24), ZB15(3,26)
C
c       DIMENSION
c     .   ZC1(3,1)  ,             ZC3(3,2)  , ZC4(3,4)  , ZC5(3,6)  ,
c     .   ZC6(3,7)  , ZC7(3,8)  , ZC8(3,10) , ZC9(3,13) , ZC10(3,14),
c     .   ZC11(3,17), ZC12(3,18), ZC13(3,22), ZC14(3,24), ZC15(3,26)
C
c       DIMENSION
c     .   ZP1(3,1)  ,             ZP3(3,2)  , ZP4(3,4)  , ZP5(3,6)  ,
c     .   ZP6(3,7)  , ZP7(3,8)  , ZP8(3,10) , ZP9(3,13) , ZP10(3,14),
c     .   ZP11(3,17), ZP12(3,18), ZP13(3,22), ZP14(3,24), ZP15(3,26),
c     .   ZP16(3,28)
C
c       DIMENSION
c     .   IQ1(3,1)  ,             IQ3(3,2)  , IQ4(3,4)  , IQ5(3,6)  ,
c     .   IQ6(3,7)  , IQ7(3,8)  , IQ8(3,10) , IQ9(3,13) , IQ10(3,14),
c     .   IQ11(3,17), IQ12(3,18), IQ13(3,22), IQ14(3,24), IQ15(3,26)
C
c       DIMENSION
c     .   ZA(3,28,3),ZB(3,28,3),ZC(3,28,3),ZP(3,28,3),IQ(3,28,3)
C
C-----------------------------------------------------------------------
CL                  C1.6     ATOMIC STRUCTURE
C
       INTEGER
     .   IZ1(1)  ,           IZ3(2)  , IZ4(4)  , IZ5(6)  ,
     .   IZ6(7)  , IZ7(8)  , IZ8(10) , IZ9(13) , IZ10(14),
     .   IZ11(17), IZ12(18), IZ13(22), IZ14(24), IZ15(26),
     .   IZ16(28)
C
       INTEGER
     .   IK1(1)  ,           IK3(2)  , IK4(4)  , IK5(6)  ,
     .   IK6(7)  , IK7(8)  , IK8(10) , IK9(13) , IK10(14),
     .   IK11(17), IK12(18), IK13(22), IK14(24), IK15(26),
     .   IK16(28)
C
       INTEGER
     .   IZ(28,3),IK(28,3)
c
c       DIMENSION
c     .   IZ1(1)  ,           IZ3(2)  , IZ4(4)  , IZ5(6)  ,
c     .   IZ6(7)  , IZ7(8)  , IZ8(10) , IZ9(13) , IZ10(14),
c     .   IZ11(17), IZ12(18), IZ13(22), IZ14(24), IZ15(26),
c     .   IZ16(28)
C
c       DIMENSION
c     .   IK1(1)  ,           IK3(2)  , IK4(4)  , IK5(6)  ,
c     .   IK6(7)  , IK7(8)  , IK8(10) , IK9(13) , IK10(14),
c     .   IK11(17), IK12(18), IK13(22), IK14(24), IK15(26),
c     .   IK16(28)
C
c       DIMENSION
c     .   IZ(28,3),IK(28,3)
C
C-----------------------------------------------------------------------
CL                  C1.7     BURGESS-CHIDICHIMO IONISATION
C
       REAL
     .   ZCON5(6)     ,  ZCON6(7)     ,  ZCON7(8)    ,  ZCON12(18)  ,
     .   ZCON15(26)   ,  ZCON16(28)  
       INTEGER
     .   KSI5(6,6)    ,  KSI6(6,7)    ,  KSI7(6,8)   ,  KSI12(6,18) ,
     .   KSI15(6,26)  ,  KSI16(6,28) 
       REAL 
     .   POT5(6,6)    ,  POT6(6,7)    ,  POT7(6,8)   ,  POT12(6,18) ,
     .   POT15A(6,13) ,  POT15B(6,13) ,  POT16A(6,14),  POT16B(6,14)
C
       REAL     ZCON(28,3) 
       INTEGER  KSI(6,28,3)
       REAL     POT(6,28,3)
c
c       DIMENSION
c     .   ZCON5(6)     ,  ZCON6(7)     ,  ZCON7(8)    ,  ZCON12(18)  ,
c     .   ZCON15(26)   ,  ZCON16(28)   ,
c     .   KSI5(6,6)    ,  KSI6(6,7)    ,  KSI7(6,8)   ,  KSI12(6,18) ,
c     .   KSI15(6,26)  ,  KSI16(6,28)  ,
c     .   POT5(6,6)    ,  POT6(6,7)    ,  POT7(6,8)   ,  POT12(6,18) ,
c     .   POT15A(6,13) ,  POT15B(6,13) ,  POT16A(6,14),  POT16B(6,14)
C
c       DIMENSION
c     .   ZCON(28,3),  KSI(6,28,3),  POT(6,28,3)
C
C-----------------------------------------------------------------------
CL                  C1.8     TRANSITION DATA
C
       REAL
     .   ZE1(3,1)  ,             ZE3(3,2)  , ZE4(3,4)  , ZE5(3,6)  ,
     .   ZE6(3,7)  , ZE7(3,8)  , ZE8(3,10) , ZE9(3,13) , ZE10(3,14),
     .   ZE11(3,17), ZE12(3,18), ZE13(3,22), ZE14(3,24), ZE15(3,26),
     .   ZE16(3,28)
C
       REAL
     .   ZF1(3,1)  ,             ZF3(3,2)  , ZF4(3,4)  , ZF5(3,6)  ,
     .   ZF6(3,7)  , ZF7(3,8)  , ZF8(3,10) , ZF9(3,13) , ZF10(3,14),
     .   ZF11(3,17), ZF12(3,18), ZF13(3,22), ZF14(3,24), ZF15(3,26),
     .   ZF16(3,28)
C
       REAL
     .   ZG1(3,1)  ,             ZG3(3,2)  , ZG4(3,4)  , ZG5(3,6)  ,
     .   ZG6(3,7)  , ZG7(3,8)  , ZG8(3,10) , ZG9(3,13) , ZG10(3,14),
     .   ZG11(3,17), ZG12(3,18), ZG13(3,22), ZG14(3,24), ZG15(3,26),
     .   ZG16(3,28)
C
       REAL
     .   DE(3,28,3)   ,  F(3,28,3),   G(3,28,3)
c
c       DIMENSION
c     .   ZE1(3,1)  ,             ZE3(3,2)  , ZE4(3,4)  , ZE5(3,6)  ,
c     .   ZE6(3,7)  , ZE7(3,8)  , ZE8(3,10) , ZE9(3,13) , ZE10(3,14),
c     .   ZE11(3,17), ZE12(3,18), ZE13(3,22), ZE14(3,24), ZE15(3,26),
c     .   ZE16(3,28)
C
c       DIMENSION
c     .   ZF1(3,1)  ,             ZF3(3,2)  , ZF4(3,4)  , ZF5(3,6)  ,
c     .   ZF6(3,7)  , ZF7(3,8)  , ZF8(3,10) , ZF9(3,13) , ZF10(3,14),
c     .   ZF11(3,17), ZF12(3,18), ZF13(3,22), ZF14(3,24), ZF15(3,26),
c     .   ZF16(3,28)
C
c       DIMENSION
c     .   ZG1(3,1)  ,             ZG3(3,2)  , ZG4(3,4)  , ZG5(3,6)  ,
c     .   ZG6(3,7)  , ZG7(3,8)  , ZG8(3,10) , ZG9(3,13) , ZG10(3,14),
c     .   ZG11(3,17), ZG12(3,18), ZG13(3,22), ZG14(3,24), ZG15(3,26),
c     .   ZG16(3,28)
C
c       DIMENSION
c     .   DE(3,28,3)   ,  F(3,28,3),   G(3,28,3)
C
C-----------------------------------------------------------------------
CL                  C1.9     IMPURITY IDENTIFICATION
C
       INTEGER
     . INAM1(3,1)  ,INAM2(3,1)  , INAM3(3,2) , INAM4(3,4) ,INAM5(3,6)  ,
     . INAM6(3,7)  ,INAM7(3,8)  , INAM8(3,10), INAM9(3,13),INAM10(3,14),
     . INAM11(3,17),INAM12(3,18),INAM13(3,22),INAM14(3,24),INAM15(3,26),
     . INAM16(3,28)
C
       INTEGER
     .   INAME(3,28,3),  ZMASS(16),   ISTATE(16)
c
c       DIMENSION
c     . INAM1(3,1)  ,INAM2(3,1)  , INAM3(3,2) , INAM4(3,4) ,INAM5(3,6)  ,
c     . INAM6(3,7)  ,INAM7(3,8)  , INAM8(3,10), INAM9(3,13),INAM10(3,14),
c     . INAM11(3,17),INAM12(3,18),INAM13(3,22),INAM14(3,24),INAM15(3,26),
c     . INAM16(3,28)
C
c       DIMENSION
c     .   INAME(3,28,3),  ZMASS(16),   ISTATE(16)
C
C-----------------------------------------------------------------------
CL                  C1.10    VAINSHTEIN LINE RADIATION
C
       REAL ZDE6(10,7),ZBE6(10,7),ZCHI6(10,7)
C
       REAL ZDE(10,10,3),ZCHI(10,10,3),ZBE(10,10,3)
C
C-----------------------------------------------------------------------
C
       REAL ZDIEL(3,26,51,28,3)
c
c       DIMENSION ZDE6(10,7),ZBE6(10,7),ZCHI6(10,7)
C
c       DIMENSION ZDE(10,10,3),ZCHI(10,10,3),ZBE(10,10,3)
C
C-----------------------------------------------------------------------
C
c       DIMENSION ZDIEL(3,26,51,28,3)
C
C ----------------------------------------------------------------------
CL              1.     INITIALISE ARRAYS
C
CL     1.1 HYDROGEN AND 1.2 DEUTERIUM
C     DESCRIPTION
       DATA ZMASS(1),ISTATE(1) / 1.0 ,  1/
       DATA ZMASS(2),ISTATE(2) / 2.0 ,  1/
       DATA INAM1 / 4HH  I, 4H    , 4H    /
       DATA INAM2 / 4HD  I, 4H    , 4H    /
C
C     COEFFICIENT A FROM BEHRINGER 14/11/85
       DATA ZA1/
     . 4.0 , 0.0 , 0.0 /
C
C     COEFFICIENT B FROM BEHRINGER 14/11/85
       DATA ZB1/
     . 0.6  , 0.0  , 0.0  /
C
C     COEFFICIENT C FROM BEHRINGER 14/11/85
       DATA ZC1/
     . 0.56 , 0.0  , 0.0  /
C
C     COEFFICIENT P FROM BEHRINGER 14/11/85
       DATA ZP1/
     .    13.6  ,     1.0 ,     1.0 /
C
C     NUMBER OF EQUIVALENT ELECTRONS Q
       DATA IQ1/
     . 1 , 0 , 0 /
C
C     PRINCIPAL QUANTUM NUMBER NZ OF VALENCE SHELL
       DATA IZ1/
     . 1 /
C
C     NUMBER OF ELECTRONS IN VALENCE SHELL
       DATA IK1/
     .  0 /
C
C     EXCITATION ENERGY DE FROM BEHRINGER 14/11/85
       DATA ZE1/
     .      0.0 ,      0.0 ,      0.0 /
C
C     OSCILLATOR STRENGTH F FROM BEHRINGER 14/11/85
       DATA ZF1/
     . 0.0 , 0.0 , 0.0 /
C
C     GAUNT FACTOR G FROM BEHRINGER 14/11/85
       DATA ZG1/
     . 0.0 , 0.0 , 0.0 /
C
C-----------------------------------------------------------------------
C
CL     1.3 HELIUM
C     IMPURITY DESCRIPTION
       DATA ZMASS(3),ISTATE(3) / 4.0 ,  2/
       DATA INAM3 / 4HHE I, 4H    , 4H    , 4HHE I, 4HI   , 4H    /
C
C     COEFFICIENT A FROM BEHRINGER 14/11/85
       DATA ZA3/
     . 4.0 , 0.0 , 0.0 , 4.4 , 0.0 , 0.0 /
C
C     COEFFICIENT B FROM BEHRINGER 14/11/85
       DATA ZB3/
     . 0.75 , 0.0  , 0.0  , 0.38 , 0.0  , 0.0  /
C
C     COEFFICIENT C FROM BEHRINGER 14/11/85
       DATA ZC3/
     . 0.46 , 0.0  , 0.0  , 0.6  , 0.0  , 0.0  /
C
C     COEFFICIENT P FROM BEHRINGER 14/11/85
       DATA ZP3/
     .    24.6  ,     1.0 ,     1.0 ,    54.4   ,     1.0 ,     1.0 /
C
C     NUMBER OF EQUIVALENT ELECTRONS Q
       DATA IQ3/
     . 2 , 0 , 0 , 1 , 0 , 0 /
C
C     PRINCIPAL QUANTUM NUMBER NZ OF VALENCE SHELL
       DATA IZ3/
     . 1 , 1 /
C
C     NUMBER OF ELECTRONS IN VALENCE SHELL
       DATA IK3/
     .  1 ,  0 /
C
C     EXCITATION ENERGY DE FROM BEHRINGER 14/11/85
       DATA ZE3/
     .    21.2  ,   23.1  ,    0.0  ,    40.8  ,   48.4  ,    0.0   /
C
C     OSCILLATOR STRENGTH F FROM BEHRINGER 14/11/85
       DATA ZF3/
     . 0.6  , 0.15 , 0.0  , 0.42 , 0.15 , 0.0  /
C
C     GAUNT FACTOR G FROM BEHRINGER 14/11/85
       DATA ZG3/
     .  -0.8  , -0.25 ,  0.0  , -0.42 , -0.35 ,  0.0  /
C
C-----------------------------------------------------------------------
C
CL     1.4 BERYLLIUM
C     IMPURITY DESCRIPTION
       DATA ZMASS(4),ISTATE(4) / 9.0 ,  4/
       DATA INAM4 / 4HBE I, 4H    , 4H    , 4HBE I, 4HI   , 4H    ,
     .              4HBE I, 4HII  , 4H    , 4HBE I, 4HV   , 4H    /
C
C     COEFFICIENT A FROM BEHRINGER 01/08/85
       DATA ZA4/
     . 4.0 , 4.4 , 0.0 , 4.4 , 4.5 , 0.0 , 4.5 , 0.0 , 0.0 ,
     . 4.5 , 0.0 , 0.0 /
C
C     COEFFICIENT B FROM BEHRINGER 01/08/85
       DATA ZB4/
     . 0.7 , 0.6 , 0.0 , 0.0 , 0.4 , 0.0 , 0.3 , 0.0 , 0.0 ,
     . 0.0 , 0.0 , 0.0 /
C
C     COEFFICIENT C FROM BEHRINGER 01/08/85
       DATA ZC4/
     . 0.5 , 0.6 , 0.0 , 0.0 , 0.6 , 0.0 , 0.6 , 0.0 , 0.0 ,
     . 0.0 , 0.0 , 0.0 /
C
C     COEFFICIENT P FROM BEHRINGER 01/08/85
       DATA ZP4/
     .     9.32 ,   115.0 ,     1.0 ,    18.2   ,   125.0 ,     1.0 ,
     .   154.0  ,     1.0 ,     1.0 ,   218.0   ,     1.0 ,     1.0 /
C
C     NUMBER OF EQUIVALENT ELECTRONS Q
       DATA IQ4/
     . 2 , 2 , 0 , 1 , 2 , 0 , 2 , 0 , 0 , 1 , 0 , 0 /
C
C     PRINCIPAL QUANTUM NUMBER NZ OF VALENCE SHELL
       DATA IZ4/
     .  2 , 2 , 1 , 1 /
C
C     NUMBER OF ELECTRONS IN VALENCE SHELL
       DATA IK4/
     .  1 ,  0 ,  1 ,  0 /
C
C     EXCITATION ENERGY DE FROM BEHRINGER 01/08/85
       DATA ZE4/
     .   5.3 ,   7.5 , 0.0 ,   4.0 ,  12.0 , 0.0 ,
     . 124.0 , 140.0 , 0.0 , 163.0 , 193.0 , 0.0 /
C
C     OSCILLATOR STRENGTH F FROM BEHRINGER 01/08/85
       DATA ZF4/
     .  1.36 , 0.6  , 0.0 , 0.51 , 0.83 , 0.0 ,
     .  0.6  , 0.15 , 0.0 , 0.42 , 0.15 , 0.0 /
C
C     GAUNT FACTOR G FROM BEHRINGER 01/08/85
       DATA ZG4/
     .  1.0  ,-0.25 , 0.0 , 1.0  ,-0.2  , 0.0 ,
     . -0.4  ,-0.25 , 0.0 ,-0.42 ,-0.35 , 0.0 /
C
C     COEFFICIENT A (LOTZ) BITC FEBRUARY 1984
C      DATA ZA4/
C    . 4.0 , 4.2 , 0.0 , 4.4 , 4.4 , 0.0 , 4.5 , 0.0 , 0.0 ,
C    . 4.5 , 0.0 , 0.0 /
C
C     COEFFICIENT B (LOTZ) BITC FEBRUARY 1984
C      DATA ZB4/
C    . 0.7  , 0.6  , 0.0  , 0.0  , 0.4  , 0.0  , 0.3  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT C (LOTZ) BITC FEBRUARY 1984
C      DATA ZC4/
C    . 0.5  , 0.6  , 0.0  , 0.0  , 0.6  , 0.0  , 0.6  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT P (LOTZ) BITC FEBRUARY 1984
C      DATA ZP4/
C    .     9.3  ,   115.0 ,     1.0 ,    18.2   ,   126.0 ,     1.0 ,
C    .   153.9  ,     1.0 ,     1.0 ,   217.7   ,     1.0 ,     1.0 /
C
C     CHI (SOBELMAN) BITC FEBRUARY 1984
C      DATA ZE4/
C    . 0.06 , 0.17 , 0.0  , 0.95 , 0.0  , 0.0  , 0.72 , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  /
C
C     Q (SOBELMAN) BITC FEBRUARY 1984
C      DATA ZF4/
C    . 1.0 , 1.0 , 0.0 , 2.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 ,
C    . 0.0 , 0.0 , 0.0 /
C
C     A (SOBELMAN) BITC FEBRUARY 1984
C      DATA ZG4/
C    . 16.6 ,  5.4 ,  0.0 , 25.7 ,  0.0 ,  0.0 , 36.7 ,  0.0 ,  0.0 ,
C    .  0.0 ,  0.0 ,  0.0 /
C
C     TRANSITION ENERGY (VAINSHTEIN) BITC FEBRUARY 1984
C      DATA ZDE4/
C    .   5.28 ,   7.46 ,   7.46 ,   2.73 ,   7.46 ,
C    .   7.46 ,   7.46 ,   0.0  ,   0.0  ,   0.0  ,
C    .   3.96 ,  11.96 ,  10.94 ,  12.16 ,   0.0  ,
C    .   0.0  ,   0.0  ,   0.0  ,   0.0  ,   0.0  ,
C    . 121.65 , 123.67 , 139.82 , 140.4  , 140.28 ,
C    . 118.59 , 121.92 , 139.01 , 139.89 , 141.14 ,
C    . 163.29 , 163.29 , 193.53 , 193.53 , 193.53 ,
C    . 204.11 , 204.11 , 204.11 ,   0.0  ,   0.0  /
C
C     B (VAINSHTEIN) BITC FEBRUARY 1984
C      DATA ZBE4/
C    . 42.8   ,  6.94  , 15.56  ,  3.72  ,  1.22  ,
C    .  4.59  , 26.7   ,  0.0   ,  0.0   ,  0.0   ,
C    .  8.29  ,  0.72  ,  3.81  ,  7.52  ,  0.0   ,
C    .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
C    .  5.98  , 21.4   ,  4.52  , 15.98  ,  0.88  ,
C    .  2.655 , 16.65  ,  2.805 , 17.4   ,  4.23  ,
C    .  5.23  , 20.4   ,  4.29  , 16.2   ,  1.23  ,
C    .  4.03  , 15.2   ,  1.51  ,  0.0   ,  0.0   /
C
C     CHI (VAINSHTEIN) BITC FEBRUARY 1984
C      DATA ZCHI4/
C    .  0.24  ,  0.77  ,  1.14  , -0.22  , -0.27  ,
C    . -0.35  , -0.8   ,  0.0   ,  0.0   ,  0.0   ,
C    .  0.389 ,  0.084 ,  0.722 ,  0.484 ,  0.0   ,
C    .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
C    .  0.55  ,  0.048 ,  0.51  ,  0.039 ,  0.025 ,
C    . -0.522 , -0.877 , -0.592 , -0.97  , -1.45  ,
C    .  0.89  ,  0.288 ,  0.911 ,  0.323 ,  0.388 ,
C    .  0.917 ,  0.335 ,  0.397 ,  0.0   ,  0.0   /
C
C-----------------------------------------------------------------------
C
CL     1.5 CARBON
C     IMPURITY DESCRIPTION
       DATA ZMASS(5),ISTATE(5) /12.0 ,  6/
       DATA INAM5 / 4HC  I, 4H    , 4H    , 4HC  I, 4HI   , 4H    ,
     .              4HC  I, 4HII  , 4H    , 4HC  I, 4HV   , 4H    ,
     .              4HC  V, 4H    , 4H    , 4HC  V, 4HI   , 4H    /
C
C     COEFFICIENT A FROM BITC+IMPUDI FOR 0+, 1+, 5+ ONLY 02/09/85
       DATA ZA5/
     . 3.5 , 4.0 , 0.0 , 4.2 , 4.4 , 0.0 , 0.0 , 0.0 , 0.0 ,
     . 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 4.5 , 0.0 , 0.0 /
C
C     COEFFICIENT B FROM BITC+IMPUDI FOR 0+, 1+, 5+ ONLY 02/09/85
       DATA ZB5/
     . 0.7  , 0.7  , 0.0  , 0.4  , 0.4  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT C FROM BITC+IMPUDI FOR 0+, 1+, 5+ ONLY 02/09/85
       DATA ZC5/
     . 0.4  , 0.5  , 0.0  , 0.6  , 0.6  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT P FROM BITC+IMPUDI 14/11/84
       DATA ZP5/
     .    11.3  ,    16.6 ,     1.0 ,    24.4   ,    30.9 ,     1.0 ,
     .    47.9  ,   325.0 ,     1.0 ,    64.5   ,   343.0 ,     1.0 ,
     .   392.0  ,     1.0 ,     1.0 ,   490.0   ,     1.0 ,     1.0 /
C
C     EQUIVALENT ELECTRONS FROM BITC+IMPUDI FOR 0+, 1+, 5+ ONLY 02/09/85
       DATA IQ5/
     . 2 , 2 , 0 , 1 , 2 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
     . 1 , 0 , 0 /
C
C     KSI (BURGESS-CHIDICHIMO) FOR 2+, 3+, 4+ ONLY 02/09/85
       DATA KSI5/
     .  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     .  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,  2 ,  1 ,  0 ,  0 ,  0 ,  0 ,
     .  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /
C
C     POTENTIAL (BURGESS-CHIDICHIMO) FOR 2+, 3+, 4+ ONLY 02/09/85
       DATA POT5/
     .     1.0 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .     1.0 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .    47.9 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .   300.0 ,   64.5 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .   392.1 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .     1.0 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   /
C
C     CONSTANT C (BURGESS-CHIDICHIMO) FOR 2+, 3+, 4+ ONLY 02/09/85
       DATA ZCON5/
     .  0.0  , 0.0  , 2.56 , 1.82 , 2.28 , 0.0  /
C
C     PRINCIPAL QUANTUM NUMBER NZ OF VALENCE SHELL BEHRINGER 14/11/84
       DATA IZ5/
     .  2 , 2 , 2 , 2 , 1 , 1 /
C
C     NUMBER OF ELECTRONS IN VALENCE SHELL BEHRINGER 14/11/84
       DATA IK5/
     .  3 ,  2 ,  1 ,  0 ,  1 ,  0 /
C
C     EXCITATION ENERGY DE FROM BEHRINGER+IMPUDI 14/11/84
       DATA ZE5/
     .     0.0  ,     0.0  ,     0.0  ,     9.3  ,    13.7  ,     0.0  ,
     .    13.0  ,    32.0  ,     0.0  ,     8.0  ,    40.0  ,     0.0  ,
     .   308.0  ,   354.0  ,     0.0  ,   367.0  ,   435.0  ,     0.0  /
C
C     OSCILLATOR STRENGTH F FROM BEHRINGER+IMPUDI 14/11/84
       DATA ZF5/
     .  0.0   ,  0.0   ,  0.0   ,  0.27  ,  0.58  ,  0.0   ,
     .  0.81  ,  0.27  ,  0.0   ,  0.29  ,  0.24  ,  0.0   ,
     .  0.65  ,  0.14  ,  0.0   ,  0.42  ,  0.15  ,  0.0   /
C
C     GAUNT FACTOR G FROM BEHRINGER (TYPING ERROR CORRECTED 14/11/84)
       DATA ZG5/
     . 0.0  , 0.0  , 0.0  , 0.8  , 0.75 , 0.0  ,
     . 0.8  ,-0.25 , 0.0  , 1.0  ,-0.1  , 0.0  ,
     .-0.8  ,-0.25 , 0.0  ,-0.42 ,-0.35 , 0.0  /
C
C     COEFFICIENT A (LOTZ) AS USED IN BITC AND IMPUDI 14/11/84
C      DATA ZA5/
C    . 3.5 , 4.0 , 0.0 , 4.2 , 4.4 , 0.0 , 4.5 , 4.5 , 0.0 ,
C    . 4.5 , 4.5 , 0.0 , 4.5 , 0.0 , 0.0 , 4.5 , 0.0 , 0.0 /
C
C     COEFFICIENT B (LOTZ) AS USED IN BITC AND IMPUDI 14/11/84
C      DATA ZB5/
C    . 0.7  , 0.7  , 0.0  , 0.4  , 0.4  , 0.0  , 0.2  , 0.2  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT C (LOTZ) AS USED IN BITC AND IMPUDI 14/11/84
C      DATA ZC5/
C    . 0.4  , 0.5  , 0.0  , 0.6  , 0.6  , 0.0  , 0.6  , 0.6  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT Q (LOTZ) AS USED IN BITC AND IMPUDI 14/11/84
C      DATA IQ5/
C    . 2 , 2 , 0 , 1 , 2 , 0 , 2 , 2 , 0 , 1 , 2 , 0 , 2 , 0 , 0 ,
C    . 1 , 0 , 0 /
C
C     EXCITATION ENERGY (BEIGMAN) (NOT TESTED)
C      DATA ZE5/
C    .     0.0  ,     0.0  ,     0.0  ,    11.0  ,     0.0  ,     0.0  ,
C    .    12.7  ,     0.0  ,     0.0  ,     8.0  ,     0.0  ,     0.0  ,
C    .   307.0  ,     0.0  ,     0.0  ,   366.0  ,     0.0  ,     0.0  ,
C
C     M (BEIGMAN) (NOT TESTED)
C      DATA ZF5/
C    . 0.0   , 0.0   , 0.0   , 1.666 , 0.0   , 0.0   ,
C    . 2.0   , 0.0   , 0.0   , 1.0   , 0.0   , 0.0   ,
C    . 2.0   , 0.0   , 0.0   , 1.0   , 0.0   , 0.0   /
C
C     C (BEIGMAN) (NOT TESTED)
C      DATA ZG5/
C    . 0.0   , 0.0   , 0.0   , 0.82  , 0.0   , 0.0   ,
C    . 1.025 , 0.0   , 0.0   , 1.06  , 0.0   , 0.0   ,
C    . 0.052 , 0.0   , 0.0   , 0.072 , 0.0   , 0.0   /
C
C     TRANSITION ENERGY (VAINSHTEIN) (NOT TESTED) (CHANGE ARRAY)
C      DATA ZDE5/
C    .   0.0  ,   0.0  ,   0.0  ,   0.0  ,   0.0  ,
C    .   0.0  ,   0.0  ,   0.0  ,   0.0  ,   0.0  ,
C    .  11.0  ,   0.0  ,   0.0  ,   0.0  ,   0.0  ,
C    .   0.0  ,   0.0  ,   0.0  ,   0.0  ,   0.0  ,
C    .  12.7  ,  32.0  ,  30.5  ,  34.3  ,   6.5  ,
C    .  32.3  ,  29.5  ,  33.5  ,   0.0  ,   0.0  ,
C    .   8.0  ,  39.5  ,  37.5  ,  40.0  ,   0.0  ,
C    .   0.0  ,   0.0  ,   0.0  ,   0.0  ,   0.0  ,
C    . 307.0  , 307.0  , 355.0  , 355.0  , 355.0  ,
C    . 304.0  , 304.0  , 355.0  , 355.0  , 355.0  ,
C    . 367.2  , 367.2  , 435.6  , 435.6  , 435.6  ,
C    . 459.0  , 459.0  , 459.0  ,   0.0  ,   0.0  /
C
C     B (VAINSHTEIN) (NOT TESTED) (CHANGE ARRAY)
C      DATA ZBE5/
C    .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
C    .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
C    . 16.6   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
C    .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
C    . 16.0   ,  1.8   ,  8.7   , 16.0   ,  1.8   ,
C    .  2.1   ,  1.57  , 11.25  ,  0.0   ,  0.0   ,
C    .  6.2   ,  1.5   ,  4.5   ,  9.2   ,  0.0   ,
C    .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
C    .  5.98  , 21.4   ,  4.52  , 15.98  ,  0.88  ,
C    .  2.655 , 16.65  ,  2.805 , 17.4   ,  4.23  ,
C    .  5.23  , 20.4   ,  4.29  , 16.2   ,  1.23  ,
C    .  4.03  , 15.2   ,  1.51  ,  0.0   ,  0.0   /
C
C     CHI (VAINSHTEIN) (NOT TESTED) (CHANGE ARRAY)
C      DATA ZCHI5/
C    .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
C    .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
C    .  0.37  ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
C    .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
C    .  0.55  ,  0.005 ,  0.78  ,  0.525 , -0.2   ,
C    . -0.45  , -0.67  , -0.78  ,  0.0   ,  0.0   ,
C    .  0.73  ,  0.01  ,  0.82  ,  0.59  ,  0.0   ,
C    .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
C    .  0.55  ,  0.048 ,  0.51  ,  0.039 ,  0.025 ,
C    . -0.522 , -0.877 , -0.592 , -0.97  , -1.45  ,
C    .  0.89  ,  0.288 ,  0.911 ,  0.323 ,  0.388 ,
C    .  0.917 ,  0.335 ,  0.397 ,  0.0   ,  0.0   /
C
C-----------------------------------------------------------------------
C
CL     1.6 NITROGEN
C     IMPURITY DESCRIPTION
       DATA ZMASS(6),ISTATE(6) /14.0 ,  7/
       DATA INAM6 / 4HN  I, 4H    , 4H    , 4HN  I, 4HI   , 4H    ,
     .              4HN  I, 4HII  , 4H    , 4HN  I, 4HV   , 4H    ,
     .              4HN  V, 4H    , 4H    , 4HN  V, 4HI   , 4H    ,
     .              4HN  V, 4HII  , 4H    /
C
C     COEFFICIENT A FROM BITC FOR 0+, 1+, 6+ ONLY 20/11/84
       DATA ZA6/
     . 3.2 , 4.0 , 0.0 , 3.9 , 4.4 , 0.0 , 0.0 , 0.0 , 0.0 ,
     . 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 ,
     . 4.5 , 0.0 , 0.0 /
C
C     COEFFICIENT B FROM BITC FOR 0+, 1+, 6+ ONLY 20/11/84
       DATA ZB6/
     . 0.83 , 0.7  , 0.0  , 0.46 , 0.4  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT C FROM BITC FOR 0+, 1+, 6+ ONLY 20/11/84
       DATA ZC6/
     . 0.22 , 0.5  , 0.0  , 0.62 , 0.6  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT P FROM BITC 16/11/84
       DATA ZP6/
     .    14.5  ,    20.3 ,     1.0 ,    29.6   ,    36.7 ,     1.0 ,
     .    47.4  ,    55.8 ,     1.0 ,    77.5   ,   471.0 ,     1.0 ,
     .    97.9  ,   494.0 ,     1.0 ,   552.1   ,     1.0 ,     1.0 ,
     .   667.0  ,     1.0 ,     1.0 /
C
C     EQUIVALENT ELECTRONS Q FROM BITC FOR 0+, 1+, 6+ ONLY 20/11/84
       DATA IQ6/
     . 3 , 2 , 0 , 2 , 2 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
     . 0 , 0 , 0 , 1 , 0 , 0 /
C
C     KSI (BURGESS-CHIDICHIMO) FOR 2+, 3+, 4+, 5+ ONLY 19/11/84
       DATA KSI6/
     .  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     .  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     .  2 ,  1 ,  0 ,  0 ,  0 ,  0 ,  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     .  0 ,  0 ,  0 ,  0 ,  0 ,  0 /
C
C     POTENTIAL (BURGESS-CHIDICHIMO) FOR 2+, 3+, 4+, 5+ ONLY 19/11/84
       DATA POT6/
     .     1.0 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .     1.0 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .    47.4 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .    77.5 ,   64.5 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .   420.0 ,   97.9 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .   552.1 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .     1.0 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   /
C
C     CONSTANT C (BURGESS-CHIDICHIMO) FOR 2+, 3+, 4+, 5+ ONLY 19/11/84
       DATA ZCON6/
     .  0.0  , 0.0  , 2.18 , 2.44 , 2.38 , 3.28 , 0.0  /
C
C     PRINCIPAL QUANTUM NUMBER NZ OF VALENCE SHELL FROM BITC 16/11/84
       DATA IZ6/
     . 2 , 2 , 2 , 2 , 2 , 1 , 1 /
C
C     NUMBER OF ELECTRONS IN VALENCE SHELL FROM BITC, 20/11/84
       DATA IK6/
     .  4 ,  3 ,  2 ,  1 ,  0 ,  1 ,  0 /
C
C     CHI (SOBELMAN) BITC 20/11/84
       DATA ZE6/
     .  0.19  , 0.44  , 0.29  , 0.091 , 0.29  , 0.21  ,
     .  0.08  , 0.19  , 0.0   , 0.032 , 0.15  , 0.0   ,
     .  0.79  , 0.0   , 0.0   , 0.68  , 0.0   , 0.0   ,
     .  0.0   , 0.0   , 0.0   /
C
C     Q (SOBELMAN) BITC 20/11/84
       DATA ZF6/
     .  1.333 , 2.0   , 0.667 , 1.667 , 2.0   , 0.333 ,
     .  2.0   , 2.0   , 0.0   , 1.0   , 1.0   , 0.0   ,
     .  2.0   , 0.0   , 0.0   , 1.0   , 0.0   , 0.0   ,
     .  0.0   , 0.0   , 0.0   /
C
C     A (SOBELMAN) BITC 20/11/84
       DATA ZG6/
     . 21.9  ,  1.42 , 41.1  , 12.88 ,  3.16 , 70.31 ,
     . 16.8  ,  7.0  ,  0.0  ,  4.81 ,  6.76 ,  0.0  ,
     . 19.5  ,  0.0  ,  0.0  , 28.3  ,  0.0  ,  0.0  ,
     .  0.0  ,  0.0  ,  0.0  /
C
C     TRANSITION ENERGY (VAINSHTEIN) BITC 21/11/84
       DATA ZDE6/
     .  10.9  ,  10.3  ,  12.0  ,  13.0  ,   0.0  ,
     .   0.0  ,   0.0  ,   0.0  ,   0.0  ,   0.0  ,
     .  14.5  ,  18.5  ,  18.5  ,  21.2  ,  21.3  ,
     .  23.3  ,   0.0  ,   0.0  ,   0.0  ,   0.0  ,
     .  16.1  ,  27.4  ,  30.5  ,  33.1  ,   0.0  ,
     .   0.0  ,   0.0  ,   0.0  ,   0.0  ,   0.0  ,
     .  16.2  ,  48.2  ,  46.8  ,  50.2  ,  50.3  ,
     .  53.2  ,  52.1  ,   0.0  ,   0.0  ,   0.0  ,
     .  10.0  ,  56.6  ,  59.2  ,  60.1  ,   0.0  ,
     .   0.0  ,   0.0  ,   0.0  ,   0.0  ,   0.0  ,
     . 426.4  , 419.8  , 430.7  , 426.3  , 496.7  ,
     . 494.9  , 498.0  , 496.7  , 497.6  , 497.6  ,
     . 500.3  , 500.3  , 592.9  , 593.0  , 593.0  ,
     . 625.4  , 625.4  , 625.4  ,   0.0  ,   0.0  /
C
C     B (VAINSHTEIN) BITC 21/11/84
       DATA ZBE6/
     .  6.3   ,  2.62  , 11.7   ,  6.09  ,  0.0   ,
     .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
     .  9.07  ,  1.2   ,  0.12  ,  7.13  ,  0.13  ,
     . 13.33  ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
     .  8.83  ,  0.4   ,  4.23  , 10.9   ,  0.0   ,
     .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
     .  6.6   ,  9.08  ,  0.98  ,  2.32  ,  1.61  ,
     . 17.04  ,  8.01  ,  0.0   ,  0.0   ,  0.0   ,
     .  1.8   ,  4.68  ,  2.17  ,  9.28  ,  0.0   ,
     .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
     .  5.98  ,  2.65  , 21.4   , 16.6   ,  4.52  ,
     .  2.81  , 15.98  , 17.4   ,  0.88  ,  4.23  ,
     .  5.23  , 20.4   ,  4.29  , 16.2   ,  1.23  ,
     .  4.03  , 15.2   ,  1.51  ,  0.0   ,  0.0   /
C
C     CHI (VAINSHTEIN) BITC 21/11/84
       DATA ZCHI6/
     .100.43  ,  0.12  ,  0.882 ,  0.131 ,  0.0   ,
     .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
     .100.64  ,  0.14  , -0.17  ,  0.66  , -0.04  ,
     .  0.21  ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
     .100.76  ,  0.13  ,  0.73  ,  0.31  ,  0.0   ,
     .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
     .100.96  ,  0.81  , -0.47  ,  0.005 , -0.39  ,
     .  0.56  , -0.63  ,  0.0   ,  0.0   ,  0.0   ,
     .101.0   ,  0.83  ,  0.013 ,  0.61  ,  0.0   ,
     .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
     .  0.55  , -0.522 ,  0.048 , -0.877 ,  0.51  ,
     . -0.592 ,  0.039 , -0.97  ,  0.025 , -1.45  ,
     .  0.89  ,  0.288 ,  0.911 ,  0.323 ,  0.388 ,
     .  0.917 ,  0.335 ,  0.397 ,  0.0   ,  0.0   /
C
C     COEFFICIENT A (LOTZ) (BITC, 16/11/84)
C      DATA ZA6/
C    . 3.2 , 4.0 , 0.0 , 3.9 , 4.4 , 0.0 , 4.5 , 4.5 , 0.0 ,
C    . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 0.0 , 0.0 ,
C    . 4.5 , 0.0 , 0.0 /
C
C     COEFFICIENT B (LOTZ) (BITC, 16/11/84)
C      DATA ZB6/
C    . 0.83 , 0.7  , 0.0  , 0.46 , 0.4  , 0.0  , 0.2  , 0.2  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT C (LOTZ) (BITC, 16/11/84)
C      DATA ZC6/
C    . 0.22 , 0.5  , 0.0  , 0.62 , 0.6  , 0.0  , 0.6  , 0.6  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT Q (LOTZ) (BITC, 16/11/84)
C      DATA IQ6/
C    . 3 , 2 , 0 , 2 , 2 , 0 , 1 , 2 , 0 , 2 , 2 , 0 , 1 , 2 , 0 ,
C    . 2 , 0 , 0 , 1 , 0 , 0 /
C
C-----------------------------------------------------------------------
C
CL     1.7 OXYGEN
C     IMPURITY DESCRIPTION
       DATA ZMASS(7),ISTATE(7) /16.0 ,  8/
       DATA INAM7 / 4HO  I, 4H    , 4H    , 4HO  I, 4HI   , 4H    ,
     .              4HO  I, 4HII  , 4H    , 4HO  I, 4HV   , 4H    ,
     .              4HO  V, 4H    , 4H    , 4HO  V, 4HI   , 4H    ,
     .              4HO  V, 4HII  , 4H    , 4HO  V, 4HIII , 4H    /
C
C     COEFFICIENT A FOR 0+ ,1+ ,6+ ,7+ ONLY 19/11/84
       DATA ZA7/
     . 2.8 , 4.0 , 0.0 , 3.7 , 4.4 , 0.0 , 0.0 , 0.0 , 0.0 ,
     . 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 ,
     . 4.5 , 0.0 , 0.0 , 4.5 , 0.0 , 0.0 /
C
C     COEFFICIENT B FOR 0+, 1+, 6+, 7+ ONLY 19/11/84
       DATA ZB7/
     . 0.74 , 0.7  , 0.0  , 0.6  , 0.4  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT C FOR 0+, 1+, 6+, 7+ ONLY 19/11/84
       DATA ZC7/
     . 0.24 , 0.5  , 0.0  , 0.5  , 0.6  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT P
       DATA ZP7/
     .    13.6  ,    28.5 ,     1.0 ,    35.1   ,    42.6 ,     1.0 ,
     .    54.9  ,    63.8 ,     1.0 ,    77.4   ,    87.6 ,     1.0 ,
     .   114.0  ,   644.0 ,     1.0 ,   138.0   ,   670.0 ,     1.0 ,
     .   739.0  ,     1.0 ,     1.0 ,   871.0   ,     1.0 ,     1.0 /
C
C     EQUIVALENT ELECTRONS Q FOR 0+, 1+, 6+, 7+ ONLY 19/11/84
       DATA IQ7/
     . 4 , 2 , 0 , 3 , 2 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
     . 0 , 0 , 0 , 2 , 0 , 0 , 1 , 0 , 0 /
C
C     KSI (BURGESS-CHIDICHIMO) FOR 2+, 3+, 4+, 5+ ONLY 19/11/84
       DATA KSI7/
     .  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     .  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     .  2 ,  2 ,  0 ,  0 ,  0 ,  0 ,  2 ,  1 ,  0 ,  0 ,  0 ,  0 ,
     .  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /
C
C     POTENTIAL (BURGESS-CHIDICHIMO) FOR 2+, 3+, 4+, 5+ ONLY 19/11/84
       DATA POT7/
     .     1.0 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .     1.0 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .    54.9 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .    77.4 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .   550.0 ,  113.9 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .   530.0 ,  138.1 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .     1.0 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .     1.0 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   /
C
C     CONSTANT C (BURGESS-CHIDICHIMO) FOR 2+, 3+, 4+, 5+ ONLY 19/11/84
       DATA ZCON7/
     .  0.0  , 0.0  , 2.36 , 2.25 , 2.87 , 2.61 , 0.0  , 0.0  /
C
C     PRINCIPAL QUANTUM NUMBER NZ OF VALENCE SHELL
       DATA IZ7/
     .  2 , 2 , 2 , 2 , 2 , 2 , 1 , 1 /
C
C     NUMBER OF ELECTRONS IN VALENCE SHELL
       DATA IK7/
     .  5 ,  4 ,  3 ,  2 ,  1 ,  0 ,  1 ,  0 /
C
C     EXCITATION ENERGY DE FROM IMPUDI+OI FROM OTTAVIANI 14/11/85
       DATA ZE7/
     .    12.0  ,     0.0  ,     0.0  ,    15.0  ,    26.0  ,     0.0  ,
     .    16.0  ,    40.0  ,     0.0  ,    21.0  ,    51.0  ,     0.0  ,
     .    20.0  ,    72.0  ,     0.0  ,    12.0  ,    83.0  ,     0.0  ,
     .   574.0  ,   666.0  ,     0.0  ,   653.0  ,   774.0  ,     0.0  /
C
C     OSCILLATOR STRENGTH F FROM IMPUDI+OI FROM OTTAVIANI 14/11/85
       DATA ZF7/
     .  0.14  ,  0.0   ,  0.0   ,  0.43  ,  0.43  ,  0.0   ,
     .  0.4   ,  0.65  ,  0.0   ,  0.63  ,  0.55  ,  0.0   ,
     .  0.53  ,  0.59  ,  0.0   ,  0.2   ,  0.26  ,  0.0   ,
     .  0.7   ,  0.15  ,  0.0   ,  0.42  ,  0.15  ,  0.0   /
C
C     GAUNT FACTOR G FROM IMPUDI+OI FROM OTTAVIANI 14/11/85
       DATA ZG7/
     .-0.3  , 0.0  , 0.0  , 0.75 ,-0.25 , 0.0  ,
     . 0.75 ,-0.25 , 0.0  , 0.8  ,-0.25 , 0.0  ,
     . 0.8  ,-0.25 , 0.0  , 1.1  ,-0.1  , 0.0  ,
     .-0.8  ,-0.25 , 0.0  ,-0.42 ,-0.35 , 0.0  /
C
C     COEFFICIENT A (LOTZ)
C      DATA ZA7/
C    . 2.8 , 4.0 , 0.0 , 3.7 , 4.4 , 0.0 , 4.5 , 4.5 , 0.0 ,
C    . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
C    . 4.5 , 0.0 , 0.0 , 4.5 , 0.0 , 0.0 /
C
C     COEFFICIENT B (LOTZ)
C      DATA ZB7/
C    . 0.74 , 0.7  , 0.0  , 0.6  , 0.4  , 0.0  , 0.3  , 0.2  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT C (LOTZ)
C      DATA ZC7/
C    . 0.24 , 0.5  , 0.0  , 0.5  , 0.6  , 0.0  , 0.6  , 0.6  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT Q (LOTZ)
C      DATA IQ7/
C    . 4 , 2 , 0 , 3 , 2 , 0 , 2 , 2 , 0 , 1 , 2 , 0 , 2 , 2 , 0 ,
C    . 1 , 2 , 0 , 2 , 0 , 0 , 1 , 0 , 0 /
C
C     EXCITATION ENERGY (OTTAVIANI)
C      DATA ZE7/
C    .    12.0  ,     0.0  ,     0.0  ,    14.9  ,    27.3  ,     0.0  ,
C    .    19.3  ,    40.7  ,     0.0  ,    20.5  ,    52.0  ,     0.0  ,
C    .    19.7  ,    72.0  ,     0.0  ,    12.0  ,    82.6  ,     0.0  ,
C    .   573.9  ,     0.0  ,     0.0  ,   652.8  ,     0.0  ,     0.0  /
C
C     OSCILLATOR STRENGTH (OTTAVIANI)
C      DATA ZF7/
C    .  0.14  ,  0.0   ,  0.0   ,  0.43  ,  0.43  ,  0.0   ,
C    .  0.52  ,  0.57  ,  0.0   ,  0.63  ,  0.5   ,  0.0   ,
C    .  0.53  ,  0.59  ,  0.0   ,  0.196 ,  0.262 ,  0.0   ,
C    .  0.694 ,  0.0   ,  0.0   ,  0.416 ,  0.0   ,  0.0   /
C
C     GAUNT FACTOR (OTTAVIANI)
C      DATA ZG7/
C    .-0.3  , 0.0  , 0.0  , 1.0  ,-0.3  , 0.0  ,
C    . 1.0  ,-0.3  , 0.0  , 1.0  ,-0.3  , 0.0  ,
C    . 1.0  ,-0.3  , 0.0  , 1.0  ,-0.3  , 0.0  ,
C    .-0.3  , 0.0  , 0.0  ,-0.3  , 0.0  , 0.0  /
C
C-----------------------------------------------------------------------
C
CL     1.8 NEON
C     IMPURITY DESCRIPTION
       DATA ZMASS(8),ISTATE(8) /20.0 , 10/
       DATA INAM8 / 4HNE I, 4H    , 4H    , 4HNE I, 4HI   , 4H    ,
     .              4HNE I, 4HII  , 4H    , 4HNE I, 4HV   , 4H    ,
     .              4HNE V, 4H    , 4H    , 4HNE V, 4HI   , 4H    ,
     .              4HNE V, 4HII  , 4H    , 4HNE V, 4HIII , 4H    ,
     .              4HNE I, 4HX   , 4H    , 4HNE X, 4H    , 4H    /
C
C     COEFFICIENT A FROM BITC, 16/11/84
       DATA ZA8/
     . 2.6 , 4.0 , 0.0 , 3.2 , 4.4 , 0.0 , 4.2 , 4.5 , 0.0 ,
     . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
     . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 0.0 , 0.0 ,
     . 4.5 , 0.0 , 0.0 /
C
C     COEFFICIENT B FROM BITC, 16/11/84
       DATA ZB8/
     . 0.92 , 0.7  , 0.0  , 0.83 , 0.4  , 0.0  , 0.5  , 0.2  , 0.0  ,
     . 0.2  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT C FROM BITC, 16/11/84
       DATA ZC8/
     . 0.19 , 0.5  , 0.0  , 0.48 , 0.6  , 0.0  , 0.6  , 0.6  , 0.0  ,
     . 0.6  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT P FROM BEHRINGER REINSTATED 12/08/85
       DATA ZP8/
     .    21.6  ,    48.5 ,     1.0 ,    41.1   ,    66.4 ,     1.0 ,
     .    63.5  ,    86.2 ,     1.0 ,    97.1   ,   108.0 ,     1.0 ,
     .   126.0  ,   139.0 ,     1.0 ,   158.0   ,   172.0 ,     1.0 ,
     .   207.0  ,  1073.0 ,     1.0 ,   239.0   ,  1107.0 ,     1.0 ,
     .  1196.0  ,     1.0 ,     1.0 ,  1362.0   ,     1.0 ,     1.0 /
C
C     NUMBER OF EQUIVALENT ELECTRONS Q
       DATA IQ8/
     . 6 , 2 , 0 , 5 , 2 , 0 , 4 , 2 , 0 , 3 , 2 , 0 , 2 , 2 , 0 ,
     . 1 , 2 , 0 , 2 , 2 , 0 , 1 , 2 , 0 , 2 , 0 , 0 , 1 , 0 , 0 /
C
C     PRINCIPAL QUANTUM NUMBER NZ OF VALENCE SHELL
       DATA IZ8/
     .  2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 1 , 1 /
C
C     NUMBER OF ELECTRONS IN VALENCE SHELL
       DATA IK8/
     .  7 ,  6 ,  5 ,  4 ,  3 ,  2 ,  1 ,  0 ,  1 ,  0 /
C
C     EXCITATION ENERGY DE FROM BEHRINGER REINSTATED 12/08/85
       DATA ZE8/
     .    16.84 ,     0.0 ,     0.0 ,    26.9   ,     0.0 ,     0.0 ,
     .    25.4  ,    50.0 ,     0.0 ,    22.8   ,    62.0 ,     0.0 ,
     .    26.0  ,    87.0 ,     0.0 ,    27.5   ,   101.0 ,     0.0 ,
     .    26.7  ,   127.0 ,     0.0 ,    16.0   ,   141.0 ,     0.0 ,
     .   921.0  ,  1072.0 ,     0.0 ,  1022.0   ,  1211.0 ,     0.0 /
C
C     OSCILLATOR STRENGTH F FROM BEHRINGER REINSTATED 12/08/85
       DATA ZF8/
     .  0.174 ,  0.0   ,  0.0   ,  0.33  ,  0.0   ,  0.0   ,
     .  0.26  ,  0.71  ,  0.0   ,  0.34  ,  0.6   ,  0.0   ,
     .  0.4   ,  0.82  ,  0.0   ,  0.43  ,  0.56  ,  0.0   ,
     .  0.57  ,  0.5   ,  0.0   ,  0.152 ,  1.0   ,  0.0   ,
     .  0.723 ,  0.149 ,  0.0   ,  0.42  ,  0.08  ,  0.0   /
C
C     GAUNT FACTOR G FROM BEHRINGER REINSTATED 12/08/85
       DATA ZG8/
     . -0.25 ,  0.0  ,  0.0  ,  0.8  ,  0.0  ,  0.0  ,
     .  0.8  , -0.25 ,  0.0  ,  0.8  , -0.25 ,  0.0  ,
     .  0.8  , -0.25 ,  0.0  ,  0.8  , -0.25 ,  0.0  ,
     .  0.8  , -0.25 ,  0.0  ,  0.8  , -0.25 ,  0.0  ,
     . -1.0  , -0.25 ,  0.0  , -0.5  , -0.25 ,  0.0  /
C
C     COEFFICIENT P (LOTZ) (BITC, 16/11/84)
C      DATA ZP8/
C    .    21.6  ,    48.5 ,     1.0 ,    41.0   ,    66.4 ,     1.0 ,
C    .    63.4  ,    86.2 ,     1.0 ,    97.1   ,   108.1 ,     1.0 ,
C    .   126.2  ,   138.6 ,     1.0 ,   157.9   ,   171.7 ,     1.0 ,
C    .   207.3  ,  1074.0 ,     1.0 ,   239.1   ,  1108.0 ,     1.0 ,
C    .  1195.8  ,     1.0 ,     1.0 ,  1362.2   ,     1.0 ,     1.0 /
C
C     CHI (SOBELMAN) BITC 22/11/84
C      DATA ZE8/
C    .  0.39  , 0.73  , 0.39  , 0.18  , 0.46  , 0.30  ,
C    .  0.1   , 0.32  , 0.24  , 0.06  , 0.24  , 0.2   ,
C    .  0.04  , 0.19  , 0.17  , 0.04  , 0.16  , 0.0   ,
C    .  0.02  , 0.14  , 0.0   , 0.72  , 0.0   , 0.0   ,
C    .  0.65  , 0.0   , 0.0   , 0.0   , 0.0   , 0.0   /
C
C     Q (SOBELMAN) BITC 22/11/84
C      DATA ZF8/
C    .  0.333 , 2.0   , 1.667 , 0.667 , 2.0   , 1.333 ,
C    .  1.0   , 2.0   , 1.0   , 1.333 , 2.0   , 0.667 ,
C    .  1.667 , 2.0   , 0.333 , 2.0   , 2.0   , 0.0   ,
C    .  1.0   , 1.0   , 0.0   , 2.0   , 0.0   , 0.0   ,
C    .  1.0   , 0.0   , 0.0   , 0.0   , 0.0   , 0.0   /
C
C     A (SOBELMAN) BITC 22/11/84
C      DATA ZG8/
C    . 35.0  ,  1.0  , 20.5  , 21.4  ,  2.19 , 40.7  ,
C    . 13.5  ,  3.31 , 69.2  ,  8.91 ,  4.71 , 83.3  ,
C    .  5.75 ,  5.8  , 96.6  ,  6.7  ,  8.4  ,  0.0  ,
C    .  2.19 ,  6.87 ,  0.0  , 15.9  ,  0.0  ,  0.0  ,
C    . 21.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  /
C
C     TRANSITION ENERGY (VAINSHTEIN) BITC 22/11/84
C      DATA ZDE8A/
C    .    0.0  ,    0.0  ,    0.0  ,    0.0  ,    0.0  ,
C    .    0.0  ,    0.0  ,    0.0  ,    0.0  ,    0.0  ,
C    .   26.9  ,   29.1  ,   27.2  ,   32.9  ,   30.8  ,
C    .   36.0  ,    0.0  ,    0.0  ,    0.0  ,    0.0  ,
C    .   25.4  ,   43.0  ,   41.5  ,   47.0  ,   46.0  ,
C    .   52.0  ,    0.0  ,    0.0  ,    0.0  ,    0.0  ,
C    .   22.9  ,   59.5  ,   60.6  ,   66.1  ,   66.0  ,
C    .   71.8  ,    0.0  ,    0.0  ,    0.0  ,    0.0  ,
C    .   27.1  ,   74.0  ,   75.0  ,   80.8  ,   75.0  ,
C    .   86.7  ,    0.0  ,    0.0  ,    0.0  ,    0.0  /
C
C      DATA ZDE8B/
C    .   27.8  ,   89.6  ,   95.7  ,  101.2  ,    0.0  ,
C    .    0.0  ,    0.0  ,    0.0  ,    0.0  ,    0.0  ,
C    .   26.7  ,  123.8  ,  121.3  ,  127.2  ,  127.5  ,
C    .  132.9  ,  130.8  ,    0.0  ,    0.0  ,    0.0  ,
C    .   16.0  ,  136.4  ,  140.7  ,  142.3  ,    0.0  ,
C    .    0.0  ,    0.0  ,    0.0  ,    0.0  ,    0.0  ,
C    .  916.3  ,  905.0  ,  922.0  ,  914.8  , 1071.8  ,
C    . 1069.1  , 1073.7  , 1071.8  , 1074.0  , 1073.3  ,
C    . 1021.5  , 1021.7  , 1210.8  , 1210.9  , 1211.0  ,
C    . 1277.0  , 1277.1  , 1277.1  ,    0.0  ,    0.0  /
C
C     B (VAINSHTEIN) BITC 22/11/84
C      DATA ZBE8A/
C    .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
C    .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
C    .  2.43  ,  4.02  ,  0.49  , 17.33  ,  0.81  ,
C    . 18.1   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
C    .  3.95  ,  2.16  ,  0.48  , 16.41  ,  0.22  ,
C    . 27.47  ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
C    .  3.16  ,  1.2   ,  0.53  , 13.23  ,  0.41  ,
C    . 30.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
C    .  4.45  ,  0.66  ,  0.25  ,  9.47  ,  0.44  ,
C    . 25.13  ,  0.0   ,  0.0   ,  0.0   ,  0.0   /
C
C      DATA ZBE8B/
C    .  5.05  ,  0.28  ,  4.79  , 14.37  ,  0.0   ,
C    .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
C    .  4.48  ,  9.54  ,  0.78  ,  5.46  ,  1.57  ,
C    . 18.4   ,  6.05  ,  0.0   ,  0.0   ,  0.0   ,
C    .  1.89  ,  4.83  ,  3.17  ,  9.55  ,  0.0   ,
C    .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
C    .  5.98  ,  2.65  , 21.4   , 16.6   ,  4.52  ,
C    .  2.81  , 15.98  , 17.4   ,  0.88  ,  4.23  ,
C    .  5.23  , 20.4   ,  4.29  , 16.2   ,  1.23  ,
C    .  4.03  , 15.2   ,  1.51  ,  0.0   ,  0.0   /
C
C     CHI (VAINSHTEIN) BITC 22/11/84
C      DATA ZCHI8A/
C    .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
C    .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
C    .100.48  ,  0.18  , -0.30  ,  0.68  , -0.26  ,
C    .  0.17  ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
C    .100.78  ,  0.15  , -0.24  ,  0.74  , -0.05  ,
C    .  0.26  ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
C    .100.98  ,  0.14  , -0.43  ,  0.76  , -0.04  ,
C    .  0.32  ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
C    .100.96  ,  0.14  , -0.46  ,  0.77  , -0.17  ,
C    .  0.37  ,  0.0   ,  0.0   ,  0.0   ,  0.0   /
C
C      DATA ZCHI8B/
C    .100.96  ,  0.14  ,  0.78  ,  0.4   ,  0.0   ,
C    .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
C    .100.99  ,  0.8   , -0.41  ,  0.03  , -0.41  ,
C    .  0.62  , -0.55  ,  0.0   ,  0.0   ,  0.0   ,
C    .101.0   ,  0.85  ,  0.03  ,  0.64  ,  0.0   ,
C    .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
C    .  0.55  , -0.522 ,  0.048 , -0.877 ,  0.51  ,
C    . -0.592 ,  0.039 , -0.97  ,  0.025 , -1.45  ,
C    .  0.89  ,  0.288 ,  0.911 ,  0.323 ,  0.388 ,
C    .  0.917 ,  0.335 ,  0.397 ,  0.0   ,  0.0   /
C
C-----------------------------------------------------------------------
C
CL     1.9 ALUMINIUM
C     IMPURITY DESCRIPTION
       DATA ZMASS(9),ISTATE(9) /27.0 , 13/
       DATA INAM9 / 4HAL I, 4H    , 4H    , 4HAL I, 4HI   , 4H    ,
     .              4HAL I, 4HII  , 4H    , 4HAL I, 4HV   , 4H    ,
     .              4HAL V, 4H    , 4H    , 4HAL V, 4HI   , 4H    ,
     .              4HAL V, 4HII  , 4H    , 4HAL V, 4HIII , 4H    ,
     .              4HAL I, 4HX   , 4H    , 4HAL X, 4H    , 4H    ,
     .              4HAL X, 4HI   , 4H    , 4HAL X, 4HII  , 4H    ,
     .              4HAL X, 4HIII , 4H    /
C
C     COEFFICIENT A (LOTZ)
       DATA ZA9/
     . 4.0 , 4.0 , 3.0 , 4.4 , 3.7 , 4.4 , 4.5 , 4.2 , 4.5 ,
     . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
     . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
     . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 0.0 , 0.0 ,
     . 4.5 , 0.0 , 0.0 /
C
C     COEFFICIENT B (LOTZ)
       DATA ZB9/
     . 0.3  , 0.4  , 0.9  , 0.2  , 0.8  , 0.4  , 0.0  , 0.6  , 0.2  ,
     . 0.3  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT C (LOTZ)
       DATA ZC9/
     . 0.6  , 0.6  , 0.2  , 0.6  , 0.4  , 0.6  , 0.0  , 0.5  , 0.6  ,
     . 0.6  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT P (LOTZ)
       DATA ZP9/
     .     5.99 ,    10.6 ,    77.0 ,    18.8   ,    90.0 ,   134.0 ,
     .    28.4  ,   103.0 ,   148.0 ,   120.0   ,   164.5 ,     1.0 ,
     .   153.8  ,   193.9 ,     1.0 ,   190.5   ,   225.2 ,     1.0 ,
     .   241.4  ,   257.8 ,     1.0 ,   284.6   ,   302.5 ,     1.0 ,
     .   330.2  ,   349.5 ,     1.0 ,   398.6   ,  1922.0 ,     1.0 ,
     .   441.9  ,  1968.0 ,     1.0 ,  2086.0   ,     1.0 ,     1.0 ,
     .  2304.1  ,     1.0 ,     1.0 /
C
C     COEFFICIENT Q (LOTZ)
       DATA IQ9/
     . 1 , 2 , 6 , 2 , 6 , 2 , 1 , 6 , 2 , 6 , 2 , 0 , 5 , 2 , 0 ,
     . 4 , 2 , 0 , 3 , 2 , 0 , 2 , 2 , 0 , 1 , 2 , 0 , 2 , 2 , 0 ,
     . 1 , 2 , 0 , 2 , 0 , 0 , 1 , 0 , 0 /
C
C     PRINCIPAL QUANTUM NUMBER NZ OF VALENCE SHELL
       DATA IZ9/
     .  3 , 3 , 3 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 1 , 1 /
C
C     NUMBER OF ELECTRONS IN VALENCE SHELL
       DATA IK9/
     .  2 ,  1 ,  0 ,  7 ,  6 ,  5 ,  4 ,  3 ,  2 ,  1 ,
     .  0 ,  1 ,  0 /
C
C     EXCITATION ENERGY
       DATA ZE9/
     .     0.0  ,     0.0 ,     0.0 ,     0.0   ,     0.0 ,     0.0 ,
     .     0.0  ,     0.0 ,     0.0 ,     0.0   ,     0.0 ,     0.0 ,
     .     0.0  ,     0.0 ,     0.0 ,     0.0   ,     0.0 ,     0.0 ,
     .     0.0  ,     0.0 ,     0.0 ,     0.0   ,     0.0 ,     0.0 ,
     .     0.0  ,     0.0 ,     0.0 ,     0.0   ,     0.0 ,     0.0 ,
     .     0.0  ,     0.0 ,     0.0 ,     0.0   ,     0.0 ,     0.0 ,
     .     0.0  ,     0.0 ,     0.0 /
C
C     OSCILLATOR STRENGTH
       DATA ZF9/
     .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
     .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
     .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
     .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
     .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
     .  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,
     .  0.0   ,  0.0   ,  0.0   /
C
C     GAUNT FACTOR
       DATA ZG9/
     .  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,
     .  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,
     .  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,
     .  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,
     .  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,
     .  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,
     .  0.0  ,  0.0  ,  0.0  /
C
C-----------------------------------------------------------------------
C
CL     1.10 SILICON
C     IMPURITY DESCRIPTION
       DATA ZMASS(10),ISTATE(10) /28.0 , 14/
       DATA INAM10/ 4HSI I, 4H    , 4H    , 4HSI I, 4HI   , 4H    ,
     .              4HSI I, 4HII  , 4H    , 4HSI I, 4HV   , 4H    ,
     .              4HSI V, 4H    , 4H    , 4HSI V, 4HI   , 4H    ,
     .              4HSI V, 4HII  , 4H    , 4HSI V, 4HIII , 4H    ,
     .              4HSI I, 4HX   , 4H    , 4HSI X, 4H    , 4H    ,
     .              4HSI X, 4HI   , 4H    , 4HSI X, 4HII  , 4H    ,
     .              4HSI X, 4HIII , 4H    , 4HSI X, 4HIV  , 4H    /
C
C     COEFFICIENT A FROM BEHRINGER 29/11/85
       DATA ZA10/
     . 4.0 , 4.0 , 0.0 , 4.4 , 4.4 , 0.0 , 4.5 , 4.5 , 0.0 ,
     . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
     . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
     . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
     . 4.5 , 0.0 , 0.0 , 4.5 , 0.0 , 0.0 /
C
C     COEFFICIENT B FROM BEHRINGER 29/11/85
       DATA ZB10/
     . 0.3  , 0.4  , 0.0  , 0.2  , 0.2  , 0.0  , 0.0  , 0.4  , 0.0  ,
     . 0.0  , 0.2  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT C FROM BEHRINGER 29/11/85
       DATA ZC10/
     . 0.6  , 0.6  , 0.0  , 0.6  , 0.6  , 0.0  , 0.0  , 0.6  , 0.0  ,
     . 0.0  , 0.6  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
C
C-----------------------------------------------------------------------
C
C     IONISATION POTENTIALS FOR FIRST SHELL ARE FROM BURKE ET AL.
C
C-----------------------------------------------------------------------
C
C     COEFFICIENT P FROM BEHRINGER 14/11/85
       DATA ZP10/
     .     8.151,    13.5 ,     1.0 ,    16.35  ,    22.9 ,     1.0 ,
     .    33.49 ,   133.0 ,     1.0 ,    45.14  ,   148.0 ,     1.0 ,
     .   166.77 ,   217.0 ,     1.0 ,   205.27  ,   250.0 ,     1.0 ,
     .   246.5  ,   285.0 ,     1.0 ,   303.54  ,   321.0 ,     1.0 ,
     .   351.13 ,   371.0 ,     1.0 ,   401.38  ,   423.0 ,     1.0 ,
     .   476.35 ,  2259.0 ,     1.0 ,   523.43  ,  2309.0 ,     1.0 ,
     .  2437.65 ,     1.0 ,     1.0 ,  2673.2   ,     1.0 ,     1.0 /
C
C
C     COEFFICIENT Q FROM BEHRINGER 29/11/85
       DATA IQ10/
     . 2 , 2 , 0 , 1 , 2 , 0 , 2 , 6 , 0 , 1 , 6 , 0 , 6 , 2 , 0 ,
     . 5 , 2 , 0 , 4 , 2 , 0 , 3 , 2 , 0 , 2 , 2 , 0 , 1 , 2 , 0 ,
     . 2 , 2 , 0 , 1 , 2 , 0 , 2 , 0 , 0 , 1 , 0 , 0 /
C
C     PRINCIPAL QUANTUM NUMBER NZ OF VALENCE SHELL BEHRINGER 29/11/85
       DATA IZ10/
     .  3 , 3 , 3 , 3 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 1 , 1 /
C
C     NUMBER OF EQUIVALENT ELECTRONS FROM BEHRINGER 29/11/85
       DATA IK10/
     .  3 ,  2 ,  1 ,  0 ,  7 ,  6 ,  5 ,  4 ,  3 ,  2 ,
     .  1 ,  0 ,  1 ,  0 /
C
C     EXCITATION ENERGY DE FROM BEHRINGER 14/11/85
       DATA ZE10/
     .     5.3  ,     5.1 ,     0.0 ,     9.9   ,    15.0 ,     0.0 ,
     .    10.3  ,    21.9 ,     0.0 ,     8.9   ,    27.0 ,     0.0 ,
     .   105.0  ,   128.0 ,     0.0 ,    50.0   ,   150.0 ,     0.0 ,
     .    45.0  ,   170.0 ,     0.0 ,    39.0   ,   191.0 ,     0.0 ,
     .    44.0  ,   200.0 ,     0.0 ,    48.0   ,   245.0 ,     0.0 ,
     .    41.0  ,   284.0 ,     0.0 ,    24.5   ,   303.0 ,     0.0 ,
     .  1863.0  ,  1991.0 ,     0.0 ,  2004.0   ,  2375.0 ,     0.0 /
C
C     OSCILLATOR STRENGTH F FROM BEHRINGER 14/11/85
       DATA ZF10/
     .  0.23  ,  0.002 ,  0.0   ,  2.4   ,  0.22  ,  0.0   ,
     .  1.7   ,  0.05  ,  0.0   ,  0.81  ,  0.04  ,  0.0   ,
     .  0.21  ,  1.04  ,  0.0   ,  0.11  ,  1.0   ,  0.0   ,
     .  0.19  ,  1.0   ,  0.0   ,  0.22  ,  0.4   ,  0.0   ,
     .  0.28  ,  0.65  ,  0.0   ,  0.29  ,  0.63  ,  0.0   ,
     .  0.27  ,  0.65  ,  0.0   ,  0.11  ,  0.26  ,  0.0   ,
     .  0.7   ,  0.15  ,  0.0   ,  0.42  ,  0.15  ,  0.0   /
C
C     GAUNT FACTOR G FROM BEHRINGER 14/11/85
       DATA ZG10/
     .  0.75 , -0.25 ,  0.0  ,  0.75 , -0.25 ,  0.0  ,
     .  0.75 , -0.25 ,  0.0  ,  0.75 , -0.25 ,  0.0  ,
     . -0.25 , -0.25 ,  0.0  ,  0.75 , -0.25 ,  0.0  ,
     .  0.75 , -0.25 ,  0.0  ,  0.75 , -0.25 ,  0.0  ,
     .  0.75 , -0.25 ,  0.0  ,  0.75 , -0.25 ,  0.0  ,
     .  0.75 , -0.25 ,  0.0  ,  1.1  , -0.1  ,  0.0  ,
     . -0.8  , -0.25 ,  0.0  , -0.42 , -0.35 ,  0.0  /
C
C-----------------------------------------------------------------------
C
CL     1.11 CHLORINE
C     IMPURITY DESCRIPTION
       DATA ZMASS(11),ISTATE(11) /35.5 , 17/
       DATA INAM11/ 4HCL I, 4H    , 4H    , 4HCL I, 4HI   , 4H    ,
     .              4HCL I, 4HII  , 4H    , 4HCL I, 4HV   , 4H    ,
     .              4HCL V, 4H    , 4H    , 4HCL V, 4HI   , 4H    ,
     .              4HCL V, 4HII  , 4H    , 4HCL V, 4HIII , 4H    ,
     .              4HCL I, 4HX   , 4H    , 4HCL X, 4H    , 4H    ,
     .              4HCL X, 4HI   , 4H    , 4HCL X, 4HII  , 4H    ,
     .              4HCL X, 4HIII , 4H    , 4HCL X, 4HIV  , 4H    ,
     .              4HCL X, 4HV   , 4H    , 4HCL X, 4HVI  , 4H    ,
     .              4HCL X, 4HVII , 4H    /
C
C     COEFFICIENT A FROM THE TFR GROUP, JULY 1984
       DATA ZA11/
     . 4.0   , 4.0   , 4.5   , 4.4   , 4.4   , 4.5   ,
     . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
     . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
     . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
     . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
     . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
     . 3.375 , 3.375 , 3.375 , 3.375 , 3.375 , 0.0   ,
     . 3.375 , 3.375 , 0.0   , 4.5   , 0.0   , 0.0   ,
     . 4.5   , 0.0   , 0.0   /
C
C     COEFFICIENT B FROM THE TFR GROUP, JULY 1984
       DATA ZB11/
     . 0.6  , 0.4  , 0.0  , 0.3  , 0.2  , 0.0  , 0.2  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT C FROM BEHRINGER AND TFR GROUP, JULY 1984
       DATA ZC11/
     . 0.5  , 0.6  , 0.0  , 0.6  , 0.6  , 0.0  , 0.6  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT P FROM THE TFR GROUP, JULY 1984
       DATA ZP11/
     .    13.0  ,    24.5 ,   206.0 ,    23.8   ,    36.0 ,   224.0 ,
     .    39.9  ,    49.0 ,   243.0 ,    53.5   ,    64.0 ,   263.0 ,
     .    67.6  ,    80.0 ,   283.0 ,    97.0   ,   303.0 ,   371.0 ,
     .   114.0  ,   322.0 ,   390.0 ,   348.0   ,   417.0 ,  3064.0 ,
     .   400.0  ,   461.0 ,  3126.0 ,   456.0   ,   507.0 ,  3190.0 ,
     .   529.0  ,   553.0 ,  3252.0 ,   592.0   ,   617.0 ,  3315.0 ,
     .   656.0  ,   683.0 ,  3377.0 ,   750.0   ,  3437.0 ,     1.0 ,
     .   809.0  ,  3499.0 ,     1.0 ,  3658.0   ,     1.0 ,     1.0 ,
     .  3946.0  ,     1.0 ,     1.0 /
C
C     NUMBER OF EQUIVALENT ELECTRONS FROM THE TFR GROUP, JULY 1984
       DATA IQ11/
     . 5 , 2 , 6 , 4 , 2 , 6 , 3 , 2 , 6 , 2 , 2 , 6 , 1 , 2 , 6 ,
     . 2 , 6 , 2 , 1 , 6 , 2 , 6 , 2 , 2 , 5 , 2 , 2 , 4 , 2 , 2 ,
     . 3 , 2 , 2 , 2 , 2 , 2 , 1 , 2 , 2 , 2 , 2 , 0 , 1 , 2 , 0 ,
     . 2 , 0 , 0 , 1 , 0 , 0 /
C
C     PRINCIPAL QUANTUM NUMBER NZ OF VALENCE SHELL
       DATA IZ11/
     .  3 , 3 , 3 , 3 , 3 , 3 , 3 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 ,
     .  1 , 1 /
C
C     NUMBER OF ELECTRONS IN VALENCE SHELL
       DATA IK11/
     .  6 ,  5 ,  4 ,  3 ,  2 ,  1 ,  0 ,  7 ,  6 ,  5 ,
     .  4 ,  3 ,  2 ,  1 ,  0 ,  1 ,  0 /
C
C     EXCITATION ENERGY DE FROM BEHRINGER 04/12/85
       DATA ZE11/
     .     0.0  ,     0.0 ,     0.0 ,    15.0   ,    20.0 ,     0.0 ,
     .    12.0  ,    22.0 ,     0.0 ,    19.0   ,    27.0 ,     0.0 ,
     .    20.0  ,    38.0 ,     0.0 ,    18.5   ,     0.0 ,     0.0 ,
     .    15.4  ,    63.0 ,     0.0 ,   210.0   ,   248.0 ,     0.0 ,
     .    68.0  ,   260.0 ,     0.0 ,    60.0   ,   280.0 ,     0.0 ,
     .    52.0  ,   325.0 ,     0.0 ,    50.0   ,   350.0 ,     0.0 ,
     .    46.0  ,   393.0 ,     0.0 ,    52.0   ,   441.0 ,     0.0 ,
     .    32.3  ,   465.0 ,     0.0 ,  2789.0   ,  3271.0 ,     0.0 ,
     .  2960.0  ,  3508.0 ,     0.0 /
C
C     OSCILLATOR STRENGTH F FROM BEHRINGER 04/12/85
       DATA ZF11/
     .  0.0   ,  0.0   ,  0.0   ,  0.3   ,  1.0   ,  0.0   ,
     .  1.2   ,  1.0   ,  0.0   ,  1.5   ,  0.4   ,  0.0   ,
     .  1.5   ,  0.4   ,  0.0   ,  1.3   ,  0.0   ,  0.0   ,
     .  0.63  ,  0.12  ,  0.0   ,  0.15  ,  2.2   ,  0.0   ,
     .  0.09  ,  1.0   ,  0.0   ,  0.11  ,  1.2   ,  0.0   ,
     .  0.2   ,  1.2   ,  0.0   ,  0.25  ,  1.0   ,  0.0   ,
     .  0.26  ,  0.7   ,  0.0   ,  0.22  ,  0.7   ,  0.0   ,
     .  0.09  ,  0.5   ,  0.0   ,  0.75  ,  0.23  ,  0.0   ,
     .  0.42  ,  0.15  ,  0.0   /
C
C     GAUNT FACTOR G FROM BEHRINGER 04/12/85
       DATA ZG11/
     .  0.0  ,  0.0  ,  0.0  ,  0.75 , -0.2  ,  0.0  ,
     .  0.75 ,  0.75 ,  0.0  ,  0.75 , -0.2  ,  0.0  ,
     .  0.75 , -0.2  ,  0.0  ,  0.75 ,  0.0  ,  0.0  ,
     .  0.75 , -0.2  ,  0.0  , -0.25 , -0.2  ,  0.0  ,
     .  0.75 , -0.2  ,  0.0  ,  0.75 , -0.2  ,  0.0  ,
     .  0.75 , -0.2  ,  0.0  ,  0.75 , -0.2  ,  0.0  ,
     .  0.75 , -0.2  ,  0.0  ,  0.75 , -0.2  ,  0.0  ,
     .  1.1  , -0.1  ,  0.0  , -0.2  , -0.2  ,  0.0  ,
     . -0.2  , -0.35 ,  0.0  /
C
C     COEFFICIENT A (LOTZ) (BEHRINGER, JULY 1984)
C      DATA ZA11/
C    . 4.0 , 4.0 , 0.0 , 4.4 , 4.4 , 0.0 , 4.5 , 4.5 , 0.0 ,
C    . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
C    . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
C    . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
C    . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
C    . 4.5 , 0.0 , 0.0 , 4.5 , 0.0 , 0.0 /
C
C     COEFFICIENT B (LOTZ) (BEHRINGER, JULY 1984)
C      DATA ZB11/
C    . 0.6  , 0.4  , 0.0  , 0.3  , 0.2  , 0.0  , 0.2  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT P (LOTZ) (BEHRINGER, JULY 1984)
C      DATA ZP11/
C    .    13.0  ,    24.5 ,     1.0 ,    23.8   ,    36.0 ,     1.0 ,
C    .    39.9  ,    48.9 ,     1.0 ,    53.5   ,    64.1 ,     1.0 ,
C    .    67.6  ,    79.8 ,     1.0 ,    97.0   ,   303.0 ,     1.0 ,
C    .   114.0  ,   322.0 ,     1.0 ,   348.0   ,   417.0 ,     1.0 ,
C    .   400.0  ,   461.0 ,     1.0 ,   456.0   ,   507.0 ,     1.0 ,
C    .   529.0  ,   553.0 ,     1.0 ,   592.0   ,   617.0 ,     1.0 ,
C    .   656.0  ,   683.0 ,     1.0 ,   750.0   ,  3437.0 ,     1.0 ,
C    .   809.0  ,  3499.0 ,     1.0 ,  3658.0   ,     1.0 ,     1.0 ,
C    .  3946.0  ,     1.0 ,     1.0 /
C
C     COEFFICIENT Q (LOTZ) (BEHRINGER, JULY 1984)
C      DATA IQ11/
C    . 5 , 2 , 0 , 4 , 2 , 0 , 3 , 2 , 0 , 2 , 2 , 0 , 1 , 2 , 0 ,
C    . 2 , 6 , 0 , 1 , 6 , 0 , 6 , 2 , 0 , 5 , 2 , 0 , 4 , 2 , 0 ,
C    . 3 , 2 , 0 , 2 , 2 , 0 , 1 , 2 , 0 , 2 , 2 , 0 , 1 , 2 , 0 ,
C    . 2 , 0 , 0 , 1 , 0 , 0 /
C
C-----------------------------------------------------------------------
C
CL     1.12 ARGON
C     IMPURITY DESCRIPTION
       DATA ZMASS(12),ISTATE(12) /39.9 , 18/
       DATA INAM12/ 4HAR I, 4H    , 4H    , 4HAR I, 4HI   , 4H    ,
     .              4HAR I, 4HII  , 4H    , 4HAR I, 4HV   , 4H    ,
     .              4HAR V, 4H    , 4H    , 4HAR V, 4HI   , 4H    ,
     .              4HAR V, 4HII  , 4H    , 4HAR V, 4HIII , 4H    ,
     .              4HAR I, 4HX   , 4H    , 4HAR X, 4H    , 4H    ,
     .              4HAR X, 4HI   , 4H    , 4HAR X, 4HII  , 4H    ,
     .              4HAR X, 4HIII , 4H    , 4HAR X, 4HIV  , 4H    ,
     .              4HAR X, 4HV   , 4H    , 4HAR X, 4HVI  , 4H    ,
     .              4HAR X, 4HVII , 4H    , 4HAR X, 4HVIII, 4H    /
C
C     COEFFICIENT A FROM TFR 07/84 NOT FOR 1+,2+,3+,4+,5+ 12/08/85
       DATA ZA12/
     . 4.0   , 4.0   , 0.0   , 0.0   , 0.0   , 0.0   ,
     . 0.0   , 0.0   , 0.0   , 0.0   , 0.0   , 0.0   ,
     . 0.0   , 0.0   , 0.0   , 0.0   , 0.0   , 0.0   ,
     . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
     . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
     . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
     . 4.5   , 4.5   , 4.5   , 3.375 , 3.375 , 3.375 ,
     . 3.375 , 3.375 , 0.0   , 3.375 , 3.375 , 0.0   ,
     . 4.5   , 0.0   , 0.0   , 4.5   , 0.0   , 0.0   /
C
C     COEFFICIENT B FROM TFR 07/84 NOT FOR 1+,2+,3+,4+,5+  12/08/85
       DATA ZB12/
     . 0.62 , 0.4  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT C FROM TFR 07/84 NOT FOR 1+,2+,3+,4+,5+  12/08/85
       DATA ZC12/
     . 0.4  , 0.6  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT Q FROM TFR 07/84 NOT FOR 1+,2+,3+,4+,5+  12/08/85
       DATA IQ12/
     . 6 , 2 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
     . 0 , 0 , 0 , 2 , 6 , 2 , 1 , 6 , 2 , 6 , 2 , 2 , 5 , 2 , 2 ,
     . 4 , 2 , 2 , 3 , 2 , 2 , 2 , 2 , 2 , 1 , 2 , 2 , 2 , 2 , 0 ,
     . 1 , 2 , 0 , 2 , 0 , 0 , 1 , 0 , 0 /
C
C     COEFFICIENT P (LOTZ)
       DATA ZP12/
     .    15.8  ,    29.2 ,     1.0 ,    27.6   ,    30.6 ,     1.0 ,
     .    40.6  ,    44.5 ,     1.0 ,    59.7   ,    59.8 ,   250.0 ,
     .    75.2  ,    75.0 ,   250.0 ,    91.2   ,    91.0 ,   250.0 ,
     .   125.0  ,   275.0 ,   348.0 ,   143.0   ,   253.0 ,   326.0 ,
     .   423.0  ,   498.0 ,  3488.0 ,   479.0   ,   545.0 ,  3554.0 ,
     .   539.0  ,   594.0 ,  3621.0 ,   618.0   ,   644.0 ,  3688.0 ,
     .   686.0  ,   713.0 ,  3755.0 ,   755.0   ,   784.0 ,  3821.0 ,
     .   855.0  ,  3885.0 ,     1.0 ,   918.0   ,  3951.0 ,     1.0 ,
     .  4121.0  ,     1.0 ,     1.0 ,  4426.0   ,     1.0 ,     1.0 /
C
C     KSI (BURGESS-CHIDICHIMO) FOR 1+, 2+, 3+, 4+, 5+ 20/11/84
       DATA KSI12/
     .  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  2 ,  5 ,  0 ,  0 ,  0 ,  0 ,
     .  2 ,  4 ,  0 ,  0 ,  0 ,  0 ,  8 ,  5 ,  0 ,  0 ,  0 ,  0 ,
     .  8 ,  4 ,  0 ,  0 ,  0 ,  0 ,  8 ,  3 ,  0 ,  0 ,  0 ,  0 ,
     .  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     .  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     .  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     .  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     .  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     .  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /
C
C     POTENTIAL (BURGESS-CHIDICHIMO) FOR 1+, 2+, 3+, 4+, 5+ 20/11/84
       DATA POT12/
     .     1.0 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .    30.6 ,   27.6 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .    44.5 ,   40.7 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .   250.0 ,   59.8 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .   250.0 ,   75.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .   250.0 ,   91.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .     1.0 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .     1.0 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .     1.0 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .     1.0 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .     1.0 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .     1.0 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .     1.0 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .     1.0 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .     1.0 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .     1.0 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .     1.0 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .     1.0 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   /
C
C     CONSTANT C (BURGESS-CHIDICHIMO) FOR 1+, 2+, 3+, 4+, 5+ 20/11/84
       DATA ZCON12/
     .  0.0  , 1.86 , 2.40 , 2.11 , 2.40 ,
     .  2.72 , 0.0  , 0.0  , 0.0  , 0.0  ,
     .  0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     .  0.0  , 0.0  , 0.0  /
C
C     PRINCIPAL QUANTUM NUMBER NZ OF VALENCE SHELL
       DATA IZ12/
     .  3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 2 , 2 , 2 , 2 , 2 , 2 , 2 ,
     .  2 , 1 , 1 /
C
C     NUMBER OF ELECTRONS IN VALENCE SHELL
       DATA IK12/
     .  7 ,  6 ,  5 ,  4 ,  3 ,  2 ,  1 ,  0 ,  7 ,  6 ,
     .  5 ,  4 ,  3 ,  2 ,  1 ,  0 ,  1 ,  0 /
C
C     EXCITATION ENERGY
       DATA ZE12/
     .     0.0  ,     0.0 ,     0.0 ,    23.0   ,    13.5 ,     0.0 ,
     .    27.0  ,    14.0 ,     0.0 ,    27.5   ,    15.0 ,     0.0 ,
     .    27.5  ,    20.0 ,     0.0 ,    28.0   ,    18.0 ,     0.0 ,
     .    21.0  ,    70.0 ,     0.0 ,    17.5   ,    77.0 ,     0.0 ,
     .     0.0  ,   260.0 ,     0.0 ,    75.0   ,   330.0 ,     0.0 ,
     .    65.0  ,   360.0 ,     0.0 ,    55.0   ,   400.0 ,     0.0 ,
     .    60.0  ,   430.0 ,     0.0 ,    60.0   ,   460.0 ,     0.0 ,
     .    56.0  ,   500.0 ,     0.0 ,    35.5   ,   520.0 ,     0.0 ,
     .     0.0  ,  3140.0 ,     0.0 ,     0.0   ,  3300.0 ,     0.0 /
C
C     OSCILLATOR STRENGTH
       DATA ZF12/
     .  0.0   ,  0.0   ,  0.0   ,  4.0   ,  0.25  ,  0.0   ,
     .  3.0   ,  0.44  ,  0.0   ,  2.8   ,  0.7   ,  0.0   ,
     .  1.2   ,  0.85  ,  0.0   ,  0.5   ,  1.0   ,  0.0   ,
     .  1.2   ,  0.2   ,  0.0   ,  0.6   ,  0.12  ,  0.0   ,
     .  0.0   ,  2.0   ,  0.0   ,  0.085 ,  3.4   ,  0.0   ,
     .  0.15  ,  1.2   ,  0.0   ,  0.175 ,  1.6   ,  0.0   ,
     .  0.22  ,  1.1   ,  0.0   ,  0.22  ,  0.63  ,  0.0   ,
     .  0.2   ,  0.65  ,  0.0   ,  0.085 ,  0.35  ,  0.0   ,
     .  0.0   ,  0.77  ,  0.0   ,  0.0   ,  0.42  ,  0.0   /
C
C     GAUNT FACTOR, NO LINE RADIATION 12/08/85
       DATA ZG12/
     .  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,
     .  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,
     .  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,
     .  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,
     .  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,
     .  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,
     .  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,
     .  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,
     .  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  /
C
C     COEFFICIENT A (LOTZ) WITH CORRECTIONS FROM TFR GROUP (JULY 1984)
C      DATA ZA12/
C    . 4.0   , 4.0   , 0.0   , 2.814 , 2.948 , 0.0   ,
C    . 3.915 , 3.915 , 0.0   , 3.42  , 3.42  , 3.42  ,
C    . 3.915 , 3.915 , 3.915 , 4.41  , 4.41  , 4.41  ,
C    . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
C    . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
C    . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
C    . 4.5   , 4.5   , 4.5   , 3.375 , 3.375 , 3.375 ,
C    . 3.375 , 3.375 , 0.0   , 3.375 , 3.375 , 0.0   ,
C    . 4.5   , 0.0   , 0.0   , 4.5   , 0.0   , 0.0   /
C
C     COEFFICIENT A (LOTZ)
C      DATA ZA12/
C    . 4.0 , 4.0 , 0.0 , 4.2 , 4.4 , 0.0 , 4.5 , 4.5 , 0.0 ,
C    . 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 ,
C    . 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 ,
C    . 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 ,
C    . 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 0.0 ,
C    . 4.5 , 4.5 , 0.0 , 4.5 , 0.0 , 0.0 , 4.5 , 0.0 , 0.0 /
C
C     COEFFICIENT B (LOTZ)
C      DATA ZB12/
C    . 0.62 , 0.4  , 0.0  , 0.3  , 0.2  , 0.0  , 0.2  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT C (LOTZ)
C      DATA ZC12/
C    . 0.4  , 0.6  , 0.0  , 0.6  , 0.6  , 0.0  , 0.6  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT Q (LOTZ)
C      DATA IQ12/
C    . 6 , 2 , 0 , 5 , 2 , 0 , 4 , 2 , 0 , 3 , 2 , 8 , 2 , 2 , 8 ,
C    . 1 , 2 , 8 , 2 , 6 , 2 , 1 , 6 , 2 , 6 , 2 , 2 , 5 , 2 , 2 ,
C    . 4 , 2 , 2 , 3 , 2 , 2 , 2 , 2 , 2 , 1 , 2 , 2 , 2 , 2 , 0 ,
C    . 1 , 2 , 0 , 2 , 0 , 0 , 1 , 0 , 0 /
C
C     GAUNT FACTOR
C      DATA ZG12/
C    .  0.0  ,  0.0  ,  0.0  ,  1.0  , -1.0  ,  0.0  ,
C    .  1.0  , -1.0  ,  0.0  ,  1.0  , -1.0  ,  0.0  ,
C    .  1.0  , -1.0  ,  0.0  ,  1.0  , -1.0  ,  0.0  ,
C    .  1.0  , -1.0  ,  0.0  ,  1.0  , -1.0  ,  0.0  ,
C    .  0.0  , -1.0  ,  0.0  ,  1.0  , -1.0  ,  0.0  ,
C    .  1.0  , -1.0  ,  0.0  ,  1.0  , -1.0  ,  0.0  ,
C    .  1.0  , -1.0  ,  0.0  ,  1.0  , -1.0  ,  0.0  ,
C    .  1.0  , -1.0  ,  0.0  ,  1.0  , -1.0  ,  0.0  ,
C    .  0.0  , -1.0  ,  0.0  ,  0.0  , -1.0  ,  0.0  /
C
C-----------------------------------------------------------------------
C
CL     1.13 TITANIUM
C     IMPURITY DESCRIPTION
       DATA ZMASS(13),ISTATE(13) /48.0 , 22/
       DATA INAM13/ 4HTI I, 4H    , 4H    , 4HTI I, 4HI   , 4H    ,
     .              4HTI I, 4HII  , 4H    , 4HTI I, 4HV   , 4H    ,
     .              4HTI V, 4H    , 4H    , 4HTI V, 4HI   , 4H    ,
     .              4HTI V, 4HII  , 4H    , 4HTI V, 4HIII , 4H    ,
     .              4HTI I, 4HX   , 4H    , 4HTI X, 4H    , 4H    ,
     .              4HTI X, 4HI   , 4H    , 4HTI X, 4HII  , 4H    ,
     .              4HTI X, 4HIII , 4H    , 4HTI X, 4HIV  , 4H    ,
     .              4HTI X, 4HV   , 4H    , 4HTI X, 4HVI  , 4H    ,
     .              4HTI X, 4HVII , 4H    , 4HTI X, 4HVIII, 4H    ,
     .              4HTI X, 4HIX  , 4H    , 4HTI X, 4HX   , 4H    ,
     .              4HTI X, 4HXI  , 4H    , 4HTI X, 4HXII , 4H    /
C
C     COEFFICIENT A FROM BEHRINGER 14/11/85
       DATA ZA13/
     . 4.0 , 3.5 , 0.0 , 4.4 , 3.9 , 0.0 , 4.5 , 4.5 , 0.0 ,
     . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
     . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
     . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
     . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
     . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
     . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 0.0 , 0.0 ,
     . 4.5 , 0.0 , 0.0 /
C
C     COEFFICIENT B FROM BEHRINGER 14/11/85
       DATA ZB13/
     . 0.4  , 0.7  , 0.0  , 0.0  , 0.5  , 0.0  , 0.3  , 0.2  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT C FROM BEHRINGER 14/11/85
       DATA ZC13/
     . 0.6  , 0.4  , 0.0  , 0.0  , 0.6  , 0.0  , 0.6  , 0.6  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT P FROM BEHRINGER 14/11/85
       DATA ZP13/
     .     6.82 ,     8.0 ,     1.0 ,    13.6   ,    15.0 ,     1.0 ,
     .    27.7  ,    68.0 ,     1.0 ,    43.2   ,    83.0 ,     1.0 ,
     .    99.9  ,   124.0 ,     1.0 ,   120.0   ,   144.0 ,     1.0 ,
     .   141.0  ,   165.0 ,     1.0 ,   170.0   ,   188.0 ,     1.0 ,
     .   193.0  ,   213.0 ,     1.0 ,   217.0   ,   238.0 ,     1.0 ,
     .   266.0  ,   723.0 ,     1.0 ,   292.0   ,   753.0 ,     1.0 ,
     .   789.0  ,   890.8 ,     1.0 ,   863.0   ,   951.0 ,     1.0 ,
     .   942.0  ,  1015.0 ,     1.0 ,  1044.0   ,  1077.0 ,     1.0 ,
     .  1131.2  ,  1165.0 ,     1.0 ,  1221.0   ,  1256.0 ,     1.0 ,
     .  1346.0  ,  5959.0 ,     1.0 ,  1425.0   ,  6041.0 ,     1.0 ,
     .  6249.0  ,     1.0 ,     1.0 ,  6626.0   ,     1.0 ,     1.0 /
C
C     NUMBER OF EQUIVALENT ELECTRONS FROM BEHRINGER 14/11/85
       DATA IQ13/
     . 2 , 2 , 0 , 1 , 2 , 0 , 2 , 6 , 0 , 1 , 6 , 0 , 6 , 2 , 0 ,
     . 5 , 2 , 0 , 4 , 2 , 0 , 3 , 2 , 0 , 2 , 2 , 0 , 1 , 2 , 0 ,
     . 2 , 6 , 0 , 1 , 6 , 0 , 6 , 2 , 0 , 5 , 2 , 0 , 4 , 2 , 0 ,
     . 3 , 2 , 0 , 2 , 2 , 0 , 1 , 2 , 0 , 2 , 2 , 0 , 1 , 2 , 0 ,
     . 2 , 0 , 0 , 1 , 0 , 0 /
C
C     PRINCIPAL QUANTUM NUMBER NZ OF VALENCE SHELL BEHRINGER 14/11/85
       DATA IZ13/
     .  4 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 2 , 2 , 2 ,
     .  2 , 2 , 2 , 2 , 2 , 1 , 1 /
C
C     NUMBER OF ELECTRONS IN VALENCE SHELL 14/11/85
       DATA IK13/
     .  1 , 10 ,  9 ,  8 ,  7 ,  6 ,  5 ,  4 ,  3 ,  2 ,
     .  1 ,  0 ,  7 ,  6 ,  5 ,  4 ,  3 ,  2 ,  1 ,  0 ,
     .  1 ,  0 /
C
C     EXCITATION ENERGY DE FROM BEHRINGER 14/11/85
       DATA ZE13/
     .     3.5  ,     0.0 ,     0.0 ,     3.5   ,     0.0 ,     0.0 ,
     .     9.6  ,     0.0 ,     0.0 ,    40.0   ,    16.0 ,     0.0 ,
     .    49.0  ,     0.0 ,     0.0 ,    13.0   ,    50.0 ,     0.0 ,
     .    25.0  ,    48.0 ,     0.0 ,    25.0   ,    46.1 ,     0.0 ,
     .    25.0  ,    45.0 ,     0.0 ,    26.3   ,    42.8 ,     0.0 ,
     .    32.1  ,   140.0 ,     0.0 ,    26.5   ,   151.0 ,     0.0 ,
     .   525.0  ,   642.0 ,     0.0 ,    88.0   ,   520.0 ,     0.0 ,
     .    88.0  ,   610.0 ,     0.0 ,    75.0   ,   645.0 ,     0.0 ,
     .    70.0  ,   680.0 ,     0.0 ,    81.0   ,   800.0 ,     0.0 ,
     .    73.2  ,   816.0 ,     0.0 ,    46.1   ,   814.0 ,     0.0 ,
     .  4730.0  ,  5560.0 ,     0.0 ,  4970.0   ,  5890.0 ,     0.0 /
C
C     OSCILLATOR STRENGTH F FROM BEHRINGER 14/11/85
       DATA ZF13/
     .  1.2   ,  0.0   ,  0.0   ,  1.0   ,  0.0   ,  0.0   ,
     .  0.24  ,  0.0   ,  0.0   ,  3.4   ,  0.05  ,  0.0   ,
     .  5.3   ,  0.0   ,  0.0   ,  0.2   ,  3.0   ,  0.0   ,
     .  0.3   ,  2.2   ,  0.0   ,  0.5   ,  1.9   ,  0.0   ,
     .  0.6   ,  1.0   ,  0.0   ,  0.1   ,  0.79  ,  0.0   ,
     .  0.96  ,  1.2   ,  0.0   ,  0.48  ,  0.9   ,  0.0   ,
     .  2.9   ,  1.0   ,  0.0   ,  0.05  ,  8.0   ,  0.0   ,
     .  0.1   ,  1.6   ,  0.0   ,  0.11  ,  1.6   ,  0.0   ,
     .  0.16  ,  1.2   ,  0.0   ,  0.16  ,  0.75  ,  0.0   ,
     .  0.18  ,  2.25  ,  0.0   ,  0.07  ,  1.1   ,  0.0   ,
     .  0.8   ,  0.16  ,  0.0   ,  0.42  ,  0.08  ,  0.0   /
C
C     GAUNT FACTOR G FROM BEHRINGER 14/11/85
       DATA ZG13/
     .  0.8  ,  0.0  ,  0.0  ,  0.8  ,  0.0  ,  0.0  ,
     .  0.8  ,  0.0  ,  0.0  ,  0.65 , -0.25 ,  0.0  ,
     .  0.65 ,  0.0  ,  0.0  ,  0.8  ,  0.65 ,  0.0  ,
     .  0.77 ,  0.8  ,  0.0  ,  0.8  ,  0.8  ,  0.0  ,
     .  0.8  ,  0.8  ,  0.0  ,  0.8  ,  0.8  ,  0.0  ,
     .  0.8  , -0.25 ,  0.0  ,  0.8  , -0.25 ,  0.0  ,
     . -0.25 , -0.25 ,  0.0  ,  0.8  , -0.25 ,  0.0  ,
     .  0.8  , -0.25 ,  0.0  ,  0.8  , -0.25 ,  0.0  ,
     .  0.8  , -0.25 ,  0.0  ,  0.8  , -0.25 ,  0.0  ,
     .  0.8  , -0.25 ,  0.0  ,  0.8  , -0.25 ,  0.0  ,
     . -1.0  , -0.25 ,  0.0  , -0.5  , -0.25 ,  0.0  /
C
C-----------------------------------------------------------------------
C
CL     1.14 CHROMIUM
C     IMPURITY DESCRIPTION
       DATA ZMASS(14),ISTATE(14) /52.0 , 24/
       DATA INAM14/ 4HCR I, 4H    , 4H    , 4HCR I, 4HI   , 4H    ,
     .              4HCR I, 4HII  , 4H    , 4HCR I, 4HV   , 4H    ,
     .              4HCR V, 4H    , 4H    , 4HCR V, 4HI   , 4H    ,
     .              4HCR V, 4HII  , 4H    , 4HCR V, 4HIII , 4H    ,
     .              4HCR I, 4HX   , 4H    , 4HCR X, 4H    , 4H    ,
     .              4HCR X, 4HI   , 4H    , 4HCR X, 4HII  , 4H    ,
     .              4HCR X, 4HIII , 4H    , 4HCR X, 4HIV  , 4H    ,
     .              4HCR X, 4HV   , 4H    , 4HCR X, 4HVI  , 4H    ,
     .              4HCR X, 4HVII , 4H    , 4HCR X, 4HVIII, 4H    ,
     .              4HCR X, 4HIX  , 4H    , 4HCR X, 4HX   , 4H    ,
     .              4HCR X, 4HXI  , 4H    , 4HCR X, 4HXII , 4H    ,
     .              4HCR X, 4HXIII, 4H    , 4HCR X, 4HXIV , 4H    /
C
C     COEFFICIENT A FROM TFR 05/12/85
       DATA ZA14/
     . 4.0 , 2.7 , 4.0 , 3.4 , 4.4 , 0.0 , 4.2 , 4.5 , 4.5 ,
     . 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 ,
     . 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 ,
     . 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 ,
     . 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 ,
     . 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 ,
     . 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 0.0 ,
     . 4.5 , 4.5 , 0.0 , 4.5 , 0.0 , 0.0 , 4.5 , 0.0 , 0.0 /
C
C     COEFFICIENT B FROM TFR 05/12/85
       DATA ZB14/
     . 0.0  , 0.9  , 0.6  , 0.8  , 0.3  , 0.0  , 0.5  , 0.2  , 0.0  ,
     . 0.2  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT C FROM TFR 05/12/85
       DATA ZC14/
     . 0.0  , 0.2  , 0.4  , 0.5  , 0.6  , 0.0  , 0.6  , 0.6  , 0.0  ,
     . 0.6  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
C
C-----------------------------------------------------------------------
C
C     IONISATION POTENTIALS FOR FIRST SHELL ARE FROM BURKE ET AL.
C
C-----------------------------------------------------------------------
C
C     COEFFICIENT P FROM TFR 05/12/85
       DATA ZP14/
     .     6.766,     8.25,    48.0 ,    16.498 ,    67.0 ,    16.498,
     .    30.96 ,    86.0 ,   117.0 ,    49.1   ,   104.0 ,   136.0  ,
     .    69.46 ,   123.0 ,   154.0 ,    90.636 ,   142.0 ,   172.0  ,
     .   160.18 ,   191.0 ,   744.0 ,   184.7   ,   214.0 ,   775.0  ,
     .   209.3  ,   238.0 ,   807.0 ,   244.4   ,   266.0 ,   840.0  ,
     .   270.8  ,   294.0 ,   873.0 ,   298.0   ,   323.0 ,   906.0 ,
     .   354.8  ,   940.0 ,  1050.0 ,   384.2   ,   974.0 ,  1083.0 ,
     .  1010.6  ,  1129.0 ,  6621.0 ,  1097.0   ,  1197.0 ,  6711.0 ,
     .  1185.0  ,  1267.0 ,  6802.0 ,  1299.0   ,  1336.0 ,  6893.0 ,
     .  1396.0  ,  1434.0 ,  6984.0 ,  1496.0   ,  1535.0 ,  7074.0 ,
     .  1634.0  ,  7164.0 ,     1.0 ,  1721.4   ,  7254.0 ,     1.0 ,
     .  7482.0  ,     1.0 ,     1.0 ,  7894.8   ,     1.0 ,     1.0 /
C
C     NUMBER OF EQUIVALENT ELECTRONS FROM TFR 05/12/85
       DATA IQ14/
     . 1 , 5 , 6 , 5 , 6 , 0 , 4 , 6 , 2 , 3 , 6 , 2 , 2 , 6 , 2 ,
     . 1 , 6 , 2 , 6 , 2 , 6 , 5 , 2 , 6 , 4 , 2 , 6 , 3 , 2 , 8 ,
     . 2 , 2 , 8 , 1 , 2 , 8 , 2 , 6 , 2 , 1 , 6 , 2 , 6 , 2 , 2 ,
     . 5 , 2 , 2 , 4 , 2 , 2 , 3 , 2 , 2 , 2 , 2 , 2 , 1 , 2 , 2 ,
     . 2 , 2 , 0 , 1 , 2 , 0 , 2 , 0 , 0 , 1 , 0 , 0 /
C
C     PRINCIPAL QUANTUM NUMBER NZ OF VALENCE SHELL
       DATA IZ14/
     .  4 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 2 ,
     .  2 , 2 , 2 , 2 , 2 , 2 , 2 , 1 , 1 /
C
C     NUMBER OF ELECTRONS IN VALENCE SHELL
       DATA IK14/
     . 13 , 12 , 11 , 10 ,  9 ,  8 ,  7 ,  6 ,  5 ,  4 ,
     .  3 ,  2 ,  1 ,  0 ,  7 ,  6 ,  5 ,  4 ,  3 ,  2 ,
     .  1 ,  0 ,  1 ,  0 /
C
C     EXCITATION ENERGY (TFR, 13 SEPTEMBER 1984)
       DATA ZE14/
     .     0.0  ,     0.0 ,     0.0 ,    10.0   ,     7.0 ,     0.0 ,
     .    15.0  ,    15.0 ,     0.0 ,    20.0   ,    20.0 ,     0.0 ,
     .    20.0  ,    30.0 ,     0.0 ,    40.0   ,    50.0 ,     0.0 ,
     .    62.0  ,   110.0 ,     0.0 ,    60.0   ,    30.0 ,     0.0 ,
     .    58.0  ,    30.0 ,     0.0 ,    56.0   ,    30.0 ,     0.0 ,
     .    52.0  ,    34.0 ,     0.0 ,    50.0   ,    32.0 ,     0.0 ,
     .    37.0  ,   185.0 ,     0.0 ,    31.0   ,   200.0 ,     0.0 ,
     .   650.0  ,     0.0 ,     0.0 ,   110.0   ,   700.0 ,     0.0 ,
     .   100.0  ,   775.0 ,     0.0 ,    85.0   ,   800.0 ,     0.0 ,
     .    90.0  ,   850.0 ,     0.0 ,    85.0   ,   875.0 ,     0.0 ,
     .    83.0  ,   950.0 ,     0.0 ,    50.0   ,   990.0 ,     0.0 ,
     .  5700.0  ,     0.0 ,     0.0 ,  5900.0   ,     0.0 ,     0.0 /
C
C     OSCILLATOR STRENGTH (TFR, 13 SEPTEMBER 1984)
       DATA ZF14/
     .  0.0   ,  0.0   ,  0.0   ,  1.8   ,  0.25  ,  0.0   ,
     .  2.0   ,  0.25  ,  0.0   ,  2.3   ,  0.31  ,  0.0   ,
     .  2.45  ,  0.31  ,  0.0   ,  2.6   ,  0.26  ,  0.0   ,
     .  4.0   ,  0.15  ,  0.0   ,  3.0   ,  0.185 ,  0.0   ,
     .  2.3   ,  0.33  ,  0.0   ,  2.1   ,  0.5   ,  0.0   ,
     .  0.9   ,  0.65  ,  0.0   ,  0.38  ,  0.77  ,  0.0   ,
     .  0.9   ,  0.36  ,  0.0   ,  0.45  ,  0.21  ,  0.0   ,
     .  2.2   ,  0.0   ,  0.0   ,  0.065 ,  4.6   ,  0.0   ,
     .  0.11  ,  1.6   ,  0.0   ,  0.13  ,  1.6   ,  0.0   ,
     .  0.165 ,  1.2   ,  0.0   ,  0.16  ,  0.65  ,  0.0   ,
     .  0.17  ,  0.63  ,  0.0   ,  0.07  ,  0.375 ,  0.0   ,
     .  0.78  ,  0.0   ,  0.0   ,  0.42  ,  0.0   ,  0.0   /
C
C     GAUNT FACTOR G NOT AVAILABLE 14/11/85
       DATA ZG14/
     .  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,
     .  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,
     .  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,
     .  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,
     .  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,
     .  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,
     .  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,
     .  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,
     .  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,
     .  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,
     .  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,
     .  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  /
C
C-----------------------------------------------------------------------
C
CL     1.15 IRON
C     IMPURITY DESCRIPTION
       DATA ZMASS(15),ISTATE(15) /56.0 , 26/
       DATA INAM15/ 4HFE I, 4H    , 4H    , 4HFE I, 4HI   , 4H    ,
     .              4HFE I, 4HII  , 4H    , 4HFE I, 4HV   , 4H    ,
     .              4HFE V, 4H    , 4H    , 4HFE V, 4HI   , 4H    ,
     .              4HFE V, 4HII  , 4H    , 4HFE V, 4HIII , 4H    ,
     .              4HFE I, 4HX   , 4H    , 4HFE X, 4H    , 4H    ,
     .              4HFE X, 4HI   , 4H    , 4HFE X, 4HII  , 4H    ,
     .              4HFE X, 4HIII , 4H    , 4HFE X, 4HIV  , 4H    ,
     .              4HFE X, 4HV   , 4H    , 4HFE X, 4HVI  , 4H    ,
     .              4HFE X, 4HVII , 4H    , 4HFE X, 4HVIII, 4H    ,
     .              4HFE X, 4HIX  , 4H    , 4HFE X, 4HX   , 4H    ,
     .              4HFE X, 4HXI  , 4H    , 4HFE X, 4HXII , 4H    ,
     .              4HFE X, 4HXIII, 4H    , 4HFE X, 4HXIV , 4H    ,
     .              4HFE X, 4HXV  , 4H    , 4HFE X, 4HXVI , 4H    /
C
C     COEFFICIENT A FROM TFR FOR 0+,1+,2+ ONLY 11/11/85
       DATA ZA15/
     . 4.0   , 2.6   , 4.0   , 4.4   , 3.3   , 4.4   ,
     . 4.0   , 4.5   , 4.5   , 0.0   , 0.0   , 0.0   ,
     . 0.0   , 0.0   , 0.0   , 0.0   , 0.0   , 0.0   ,
     . 0.0   , 0.0   , 0.0   , 0.0   , 0.0   , 0.0   ,
     . 0.0   , 0.0   , 0.0   , 0.0   , 0.0   , 0.0   ,
     . 0.0   , 0.0   , 0.0   , 0.0   , 0.0   , 0.0   ,
     . 0.0   , 0.0   , 0.0   , 0.0   , 0.0   , 0.0   ,
     . 0.0   , 0.0   , 0.0   , 0.0   , 0.0   , 0.0   ,
     . 0.0   , 0.0   , 0.0   , 0.0   , 0.0   , 0.0   ,
     . 0.0   , 0.0   , 0.0   , 0.0   , 0.0   , 0.0   ,
     . 0.0   , 0.0   , 0.0   , 0.0   , 0.0   , 0.0   ,
     . 0.0   , 0.0   , 0.0   , 0.0   , 0.0   , 0.0   ,
     . 0.0   , 0.0   , 0.0   , 0.0   , 0.0   , 0.0   /
C
C     COEFFICIENT B FROM TFR FOR 0+,1+,2+ ONLY 12/11/85
       DATA ZB15/
     . 0.4  , 0.92 , 0.6  , 0.0  , 0.85 , 0.3  , 0.6  , 0.2  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT C FROM TFR FOR 0+,1+,2+ ONLY 12/11/85
       DATA ZC15/
     . 0.6  , 0.19 , 0.4  , 0.0  , 0.4  , 0.6  , 0.5  , 0.6  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT Q FROM TFR FOR 0+,1+,2+ ONLY 11/11/85
       DATA IQ15/
     . 2 , 6 , 6 , 1 , 6 , 6 , 6 , 6 , 2 , 0 , 0 , 0 , 0 , 0 , 0 ,
     . 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
     . 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
     . 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
     . 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
     . 0 , 0 , 0 /
C
C-----------------------------------------------------------------------
C
C    IONISATION POTENTIALS FOR FIRST SHELL ARE FROM BURKE ET AL.
C
C-----------------------------------------------------------------------
C
C     COEFFICIENT P (FOR 2ND AND 3RD SHELL TFR DATA, JULY 1984)
       DATA ZP15/
     .     7.87 ,     9.0 ,    59.0 ,    16.188 ,    17.5 ,    81.0 ,
     .    30.652,   103.0 ,   141.0 ,    54.8   ,   125.0 ,   162.0 ,
     .    75.0  ,   147.0 ,   184.0 ,    99.1   ,   169.0 ,   205.0 ,
     .   124.98 ,   190.0 ,   227.0 ,   151.06  ,   213.0 ,   249.0 ,
     .   233.6  ,   271.0 ,   965.0 ,   262.1   ,   297.0 ,  1000.0 ,
     .   290.3  ,   324.0 ,  1036.0 ,   330.8   ,   356.0 ,   756.0 ,
     .   361.0  ,   388.0 ,   756.0 ,   392.2   ,   421.0 ,   756.0 ,
     .   457.0  ,   756.0 ,   888.0 ,   489.26  ,   721.0 ,   854.0 ,
     .  1262.2  ,  1397.0 ,  7891.0 ,  1362.0   ,  1471.0 ,  7989.0 ,
     .  1469.0  ,  1548.0 ,  8088.0 ,  1582.0   ,  1622.0 ,  8187.0 ,
     .  1689.0  ,  1731.0 ,  8286.0 ,  1799.0   ,  1842.0 ,  8384.0 ,
     .  1958.6  ,  8482.0 ,     1.0 ,  2045.8   ,  8580.0 ,     1.0 ,
     .  8828.14 ,     1.0 ,     1.0 ,  9227.65  ,     1.0 ,     1.0 /
C
C     KSI (BURGESS-CHIDICHIMO) FOR STATES HIGHER THAN 2+ 12/11/85
       DATA KSI15/
     . 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
     . 0 , 0 , 0 , 0 , 0 , 0 , 5 , 6 , 0 , 0 , 0 , 0 ,
     . 4 , 6 , 0 , 0 , 0 , 0 , 3 , 6 , 2 , 0 , 0 , 0 ,
     . 2 , 6 , 2 , 0 , 0 , 0 , 1 , 6 , 2 , 0 , 0 , 0 ,
     . 6 , 2 , 6 , 0 , 0 , 0 , 5 , 2 , 6 , 0 , 0 , 0 ,
     . 4 , 2 , 6 , 0 , 0 , 0 , 3 , 2 , 6 , 0 , 0 , 0 ,
     . 2 , 2 , 6 , 0 , 0 , 0 , 1 , 2 , 6 , 0 , 0 , 0 ,
     . 2 , 6 , 2 , 0 , 0 , 0 , 1 , 6 , 2 , 0 , 0 , 0 ,
     . 6 , 2 , 2 , 0 , 0 , 0 , 5 , 2 , 2 , 0 , 0 , 0 ,
     . 4 , 2 , 2 , 0 , 0 , 0 , 3 , 2 , 2 , 0 , 0 , 0 ,
     . 2 , 2 , 2 , 0 , 0 , 0 , 1 , 2 , 2 , 0 , 0 , 0 ,
     . 2 , 2 , 0 , 0 , 0 , 0 , 1 , 2 , 0 , 0 , 0 , 0 ,
     . 2 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 /
C
C     POTENTIAL (BURGESS-CHIDICHIMO) FOR STATES HIGHER THAN 2+ 11/11/85
       DATA POT15A/
     .     1.0 ,    1.0 ,    1.0 ,    1.0 ,    1.0 ,    1.0 ,
     .     1.0 ,    1.0 ,    1.0 ,    1.0 ,    1.0 ,    1.0 ,
     .     1.0 ,    1.0 ,    1.0 ,    1.0 ,    1.0 ,    1.0 ,
     .    54.8 ,   54.8 ,    1.0 ,    1.0 ,    1.0 ,    1.0 ,
     .    75.0 ,   75.0 ,    1.0 ,    1.0 ,    1.0 ,    1.0 ,
     .    99.0 ,   99.0 ,   99.0 ,    1.0 ,    1.0 ,    1.0 ,
     .   125.0 ,  125.0 ,  125.0 ,    1.0 ,    1.0 ,    1.0 ,
     .   151.0 ,  151.0 ,  151.0 ,    1.0 ,    1.0 ,    1.0 ,
     .   235.0 ,  235.0 ,  782.0 ,    1.0 ,    1.0 ,    1.0 ,
     .   262.0 ,  262.0 ,  730.0 ,    1.0 ,    1.0 ,    1.0 ,
     .   290.0 ,  290.0 ,  738.0 ,    1.0 ,    1.0 ,    1.0 ,
     .   331.0 ,  331.0 ,  746.0 ,    1.0 ,    1.0 ,    1.0 ,
     .   361.0 ,  361.0 ,  742.0 ,    1.0 ,    1.0 ,    1.0 /
C
       DATA POT15B/
     .   392.0 ,  392.0 ,  743.0 ,    1.0 ,    1.0 ,    1.0 ,
     .   457.0 ,  745.0 ,  847.0 ,    1.0 ,    1.0 ,    1.0 ,
     .   489.3 ,  709.0 ,  820.0 ,    1.0 ,    1.0 ,    1.0 ,
     .  1265.0 , 1265.0 , 7348.0 ,    1.0 ,    1.0 ,    1.0 ,
     .  1358.0 , 1358.0 , 6626.0 ,    1.0 ,    1.0 ,    1.0 ,
     .  1456.0 , 1456.0 , 6631.0 ,    1.0 ,    1.0 ,    1.0 ,
     .  1582.0 , 1582.0 , 6632.0 ,    1.0 ,    1.0 ,    1.0 ,
     .  1689.0 , 1689.0 , 6605.0 ,    1.0 ,    1.0 ,    1.0 ,
     .  1799.0 , 1799.0 , 6586.0 ,    1.0 ,    1.0 ,    1.0 ,
     .  1950.0 , 6585.0 ,    1.0 ,    1.0 ,    1.0 ,    1.0 ,
     .  2045.0 , 6532.0 ,    1.0 ,    1.0 ,    1.0 ,    1.0 ,
     .  8828.0 ,    1.0 ,    1.0 ,    1.0 ,    1.0 ,    1.0 ,
     .  9278.0 ,    1.0 ,    1.0 ,    1.0 ,    1.0 ,    1.0 /
C
C     PARAMETER C BURGESS-CHIDICHIMO FOR STATES HIGHER THAN 2+ 12/11/85
       DATA ZCON15/
     .  0.0 , 0.0 , 0.0 , 2.3 , 2.3 , 2.3 , 2.3 , 2.3 , 2.3 , 2.3 ,
     .  2.3 , 2.3 , 2.3 , 2.3 , 2.3 , 2.3 , 2.3 , 2.3 , 2.3 , 2.3 ,
     .  2.3 , 2.3 , 2.3 , 2.3 , 2.3 , 2.3 /
C
C     PRINCIPAL QUANTUM NUMBER NZ OF VALENCE SHELL
       DATA IZ15/
     .  4 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 ,
     .  3 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 1 , 1 /
C
C     NUMBER OF ELECTRONS IN VALENCE SHELL
       DATA IK15/
     .  1 , 14 , 13 , 12 , 11 , 10 ,  9 ,  8 ,  7 ,  6 ,
     .  5 ,  4 ,  3 ,  2 ,  1 ,  0 ,  7 ,  6 ,  5 ,  4 ,
     .  3 ,  2 ,  1 ,  0 ,  1 ,  0 /
C
C     EXCITATION ENERGY DE FROM BEHRINGER + FEI FROM OTTAVIANI 15/08/85
       DATA ZE15/
     .     2.5  ,     5.0 ,     0.0 ,     5.0   ,     0.0 ,     0.0 ,
     .    12.0  ,     0.0 ,     0.0 ,    24.0   ,     0.0 ,     0.0 ,
     .    32.0  ,     0.0 ,     0.0 ,    42.0   ,     0.0 ,     0.0 ,
     .    50.0  ,     0.0 ,     0.0 ,    72.5   ,     0.0 ,     0.0 ,
     .    72.5  ,     0.0 ,     0.0 ,    35.0   ,    71.0 ,     0.0 ,
     .    35.0  ,    69.0 ,     0.0 ,    34.5   ,    64.5 ,     0.0 ,
     .    34.5  ,    62.0 ,     0.0 ,    35.0   ,    57.0 ,     0.0 ,
     .    43.5  ,   245.0 ,     0.0 ,    36.5   ,   260.0 ,     0.0 ,
     .   820.0  ,   900.0 ,     0.0 ,   130.0   ,   820.0 ,     0.0 ,
     .   110.0  ,   915.0 ,     0.0 ,   100.0   ,   970.0 ,     0.0 ,
     .   105.0  ,  1000.0 ,     0.0 ,    95.0   ,  1050.0 ,     0.0 ,
     .    97.0  ,  1075.0 ,     0.0 ,    57.5   ,  1150.0 ,     0.0 ,
     .  6700.0  ,  7900.0 ,     0.0 ,     0.0   ,     0.0 ,     0.0 /
C
C     OSCILLATOR STRENGTH F FROM BEHRINGER + FEI FROM OTTAVIANI 15/08/85
       DATA ZF15/
     .  2.0   ,  3.0   ,  0.0   ,  1.0   ,  0.0   ,  0.0   ,
     .  1.0   ,  0.0   ,  0.0   ,  1.0   ,  0.0   ,  0.0   ,
     .  1.0   ,  0.0   ,  0.0   ,  1.0   ,  0.0   ,  0.0   ,
     .  1.0   ,  0.0   ,  0.0   ,  3.4   ,  0.0   ,  0.0   ,
     .  3.73  ,  0.0   ,  0.0   ,  0.17  ,  2.8   ,  0.0   ,
     .  0.3   ,  2.15  ,  0.0   ,  0.47  ,  2.1   ,  0.0   ,
     .  0.6   ,  1.0   ,  0.0   ,  0.71  ,  0.35  ,  0.0   ,
     .  0.8   ,  1.23  ,  0.0   ,  0.4   ,  0.92  ,  0.0   ,
     .  2.84  ,  2.12  ,  0.0   ,  0.06  ,  8.2   ,  0.0   ,
     .  0.1   ,  1.6   ,  0.0   ,  0.12  ,  1.6   ,  0.0   ,
     .  0.15  ,  1.2   ,  0.0   ,  0.15  ,  0.7   ,  0.0   ,
     .  0.16  ,  2.3   ,  0.0   ,  0.07  ,  1.2   ,  0.0   ,
     .  0.8   ,  0.16  ,  0.0   ,  0.0   ,  0.0   ,  0.0   /
C
C     GAUNT FACTOR G FROM BEHRINGER + FEI FROM OTTAVIANI  15/08/85
       DATA ZG15/
     .  1.0  , -0.3  ,  0.0  ,  0.8  ,  0.0  ,  0.0  ,
     .  0.7  ,  0.0  ,  0.0  , -0.25 ,  0.0  ,  0.0  ,
     .  0.7  ,  0.0  ,  0.0  ,  0.7  ,  0.0  ,  0.0  ,
     . -0.25 ,  0.0  ,  0.0  ,  0.65 ,  0.0  ,  0.0  ,
     .  0.65 ,  0.0  ,  0.0  ,  0.77 ,  0.65 ,  0.0  ,
     .  0.77 ,  0.77 ,  0.0  ,  0.77 ,  0.77 ,  0.0  ,
     .  0.77 ,  0.77 ,  0.0  ,  0.77 ,  0.7  ,  0.0  ,
     .  0.77 , -0.25 ,  0.0  ,  0.77 , -0.25 ,  0.0  ,
     . -0.25 , -0.25 ,  0.0  ,  0.77 , -0.25 ,  0.0  ,
     .  0.77 , -0.25 ,  0.0  ,  0.77 , -0.25 ,  0.0  ,
     .  0.77 , -0.25 ,  0.0  ,  0.77 , -0.25 ,  0.0  ,
     .  0.77 , -0.25 ,  0.0  ,  0.77 , -0.25 ,  0.0  ,
     . -1.0  , -0.25 ,  0.0  ,  0.0  ,  0.0  ,  0.0  /
C
C     COEFFICIENT A FROM TFR JULY 1984
C      DATA ZA15/
C    . 4.0   , 2.6   , 4.0   , 4.4   , 3.3   , 4.4   ,
C    . 4.0   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
C    . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
C    . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
C    . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
C    . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
C    . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
C    . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
C    . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
C    . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
C    . 4.5   , 4.5   , 4.5   , 3.375 , 3.375 , 3.375 ,
C    . 3.375 , 3.375 , 0.0   , 3.375 , 3.375 , 0.0   ,
C    . 4.5   , 0.0   , 0.0   , 4.5   , 0.0   , 0.0   /
C
C     COEFFICIENT A (LOTZ) (BEHRINGER, JULY 1984)
C      DATA ZA15/
C    . 4.0 , 2.6 , 0.0 , 4.4 , 3.3 , 0.0 , 4.0 , 4.5 , 0.0 ,
C    . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
C    . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
C    . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
C    . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
C    . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
C    . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
C    . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
C    . 4.5 , 0.0 , 0.0 , 4.5 , 0.0 , 0.0 /
C
C     COEFFICIENT A (LOTZ) (OTTAVIANI, OCT 1983)
C      DATA ZA15/
C    . 4.0   , 2.6   , 4.0   , 4.5   , 4.5   , 4.5   ,
C    . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
C    . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
C    . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
C    . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
C    . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
C    . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
C    . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
C    . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
C    . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
C    . 4.5   , 4.5   , 4.5   , 4.5   , 4.5   , 4.5   ,
C    . 4.5   , 4.5   , 0.0   , 4.5   , 4.5   , 0.0   ,
C    . 4.5   , 0.0   , 0.0   , 4.5   , 0.0   , 0.0   /
C
C     COEFFICIENT B FROM TFR JULY 1984
C      DATA ZB15/
C    . 0.4  , 0.92 , 0.6  , 0.0  , 0.85 , 0.3  , 0.6  , 0.2  , 0.0  ,
C    . 0.3  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT B (LOTZ) (BEHRINGER, JULY 1984)
C      DATA ZB15/
C    . 0.4  , 0.92 , 0.0  , 0.0  , 0.85 , 0.0  , 0.6  , 0.2  , 0.0  ,
C    . 0.3  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT B (LOTZ) (OTTAVIANI, OCT 1983)
C      DATA ZB15/
C    . 0.4  , 0.92 , 0.6  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT C FROM TFR JULY 1984
C      DATA ZC15/
C    . 0.6  , 0.19 , 0.4  , 0.0  , 0.4  , 0.6  , 0.5  , 0.6  , 0.0  ,
C    . 0.6  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT C (LOTZ) (BEHRINGER, JULY 1984)
C      DATA ZC15/
C    . 0.6  , 0.19 , 0.0  , 0.0  , 0.4  , 0.0  , 0.5  , 0.6  , 0.0  ,
C    . 0.6  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT C (LOTZ) (OTTAVIANI, OCT 1983)
C      DATA ZC15/
C    . 0.6  , 0.19 , 0.4  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT P (LOTZ) (BEHRINGER, JULY 1984)
C      DATA ZP15/
C    .     7.87 ,     9.0 ,     1.0 ,    16.2   ,    17.5 ,     1.0 ,
C    .    30.7  ,   103.0 ,     1.0 ,    54.8   ,   125.0 ,     1.0 ,
C    .    75.0  ,   147.0 ,     1.0 ,    99.0   ,   169.0 ,     1.0 ,
C    .   125.0  ,   190.0 ,     1.0 ,   151.0   ,   213.0 ,     1.0 ,
C    .   235.0  ,   271.0 ,     1.0 ,   262.0   ,   297.0 ,     1.0 ,
C    .   290.0  ,   324.0 ,     1.0 ,   331.0   ,   356.0 ,     1.0 ,
C    .   361.0  ,   388.0 ,     1.0 ,   392.0   ,   421.0 ,     1.0 ,
C    .   457.0  ,  1185.0 ,     1.0 ,   490.0   ,  1223.0 ,     1.0 ,
C    .  1265.0  ,  1397.0 ,     1.0 ,  1358.0   ,  1471.0 ,     1.0 ,
C    .  1456.0  ,  1548.0 ,     1.0 ,  1582.0   ,  1622.0 ,     1.0 ,
C    .  1689.0  ,  1731.0 ,     1.0 ,  1799.0   ,  1842.0 ,     1.0 ,
C    .  1950.0  ,  8482.0 ,     1.0 ,  2045.0   ,  8580.0 ,     1.0 ,
C    .  8828.0  ,     1.0 ,     1.0 ,  9278.0   ,     1.0 ,     1.0 /
C
C     COEFFICIENT P (LOTZ) (FOR 2ND AND 3RD SHELL OTTAVIANI, OCT 1983)
C      DATA ZP15/
C    .     7.87 ,     9.0 ,    59.0 ,    16.188 ,    17.5 ,    81.0 ,
C    .    30.652,   103.0 ,   141.0 ,    54.8   ,   125.0 ,   162.0 ,
C    .    75.0  ,   147.0 ,   184.0 ,    99.1   ,   169.0 ,   205.0 ,
C    .   124.98 ,   190.0 ,   227.0 ,   151.06  ,   213.0 ,   249.0 ,
C    .   233.6  ,   270.8 ,   965.0 ,   262.1   ,   297.2 ,  1000.0 ,
C    .   290.3  ,   324.4 ,  1036.0 ,   330.8   ,   355.8 ,  1073.0 ,
C    .   361.0  ,   387.8 ,  1110.0 ,   392.2   ,   421.0 ,  1147.0 ,
C    .   457.0  ,  1185.0 ,  1307.0 ,   489.26  ,  1223.0 ,  1344.0 ,
C    .  1262.2  ,  1397.0 ,     1.0 ,  1362.0   ,  1471.0 ,     1.0 ,
C    .  1469.0  ,  1548.0 ,     1.0 ,  1582.0   ,  1622.0 ,     1.0 ,
C    .  1689.0  ,  1783.0 ,     1.0 ,  1799.0   ,  1842.0 ,     1.0 ,
C    .  1958.6  ,  8482.0 ,     1.0 ,  2045.8   ,  8580.0 ,     1.0 ,
C    .  8828.14 ,     1.0 ,     1.0 ,  9227.65  ,     1.0 ,     1.0 /
C
C     COEFFICIENT Q FROM TFR JULY 1984
C      DATA IQ15/
C    . 2 , 6 , 6 , 1 , 6 , 6 , 6 , 6 , 2 , 5 , 6 , 2 , 4 , 6 , 2 ,
C    . 3 , 6 , 2 , 2 , 6 , 2 , 1 , 6 , 2 , 6 , 2 , 6 , 5 , 2 , 6 ,
C    . 4 , 2 , 6 , 3 , 2 , 8 , 2 , 2 , 8 , 1 , 2 , 8 , 2 , 6 , 2 ,
C    . 1 , 6 , 2 , 6 , 2 , 2 , 5 , 2 , 2 , 4 , 2 , 2 , 3 , 2 , 2 ,
C    . 2 , 2 , 2 , 1 , 2 , 2 , 2 , 2 , 0 , 1 , 2 , 0 , 2 , 0 , 0 ,
C    . 1 , 0 , 0 /
C
C     COEFFICIENT Q (LOTZ) (BEHRINGER, JULY 1984)
C      DATA IQ15/
C    . 2 , 6 , 0 , 1 , 6 , 0 , 6 , 6 , 0 , 5 , 6 , 0 , 4 , 6 , 0 ,
C    . 3 , 6 , 0 , 2 , 6 , 0 , 1 , 6 , 0 , 6 , 2 , 0 , 5 , 2 , 0 ,
C    . 4 , 2 , 0 , 3 , 2 , 0 , 2 , 2 , 0 , 1 , 2 , 0 , 2 , 6 , 0 ,
C    . 1 , 6 , 0 , 6 , 2 , 0 , 5 , 2 , 0 , 4 , 2 , 0 , 3 , 2 , 0 ,
C    . 2 , 2 , 0 , 1 , 2 , 0 , 2 , 2 , 0 , 1 , 2 , 0 , 2 , 0 , 0 ,
C    . 1 , 0 , 0 /
C
C     COEFFICIENT Q (LOTZ) (OTTAVIANI, OCT 1983)
C      DATA IQ15/
C    . 2 , 6 , 6 , 1 , 6 , 6 , 6 , 6 , 2 , 5 , 6 , 2 , 4 , 6 , 2 ,
C    . 3 , 6 , 2 , 2 , 6 , 2 , 1 , 6 , 2 , 6 , 2 , 6 , 5 , 2 , 6 ,
C    . 4 , 2 , 6 , 3 , 2 , 6 , 2 , 2 , 6 , 1 , 2 , 6 , 2 , 6 , 2 ,
C    . 1 , 6 , 2 , 6 , 2 , 0 , 5 , 2 , 0 , 4 , 2 , 0 , 3 , 2 , 0 ,
C    . 2 , 2 , 0 , 1 , 2 , 0 , 2 , 2 , 0 , 1 , 2 , 0 , 2 , 0 , 0 ,
C    . 1 , 0 , 0 /
C
C     GAUNT FACTOR (BEHRINGER)
C      DATA ZG15/
C    .  1.0  , -0.3  ,  0.0  ,  0.8  ,  0.0  ,  0.0  ,
C    .  0.7  ,  0.0  ,  0.0  ,  0.25 ,  0.0  ,  0.0  ,
C    .  0.7  ,  0.0  ,  0.0  ,  0.7  ,  0.0  ,  0.0  ,
C    .  0.25 ,  0.0  ,  0.0  ,  0.65 ,  0.0  ,  0.0  ,
C    .  0.65 ,  0.0  ,  0.0  ,  0.77 ,  0.65 ,  0.0  ,
C    .  0.77 ,  0.77 ,  0.0  ,  0.77 ,  0.77 ,  0.0  ,
C    .  0.77 ,  0.77 ,  0.0  ,  0.77 ,  0.7  ,  0.0  ,
C    .  0.77 ,  0.25 ,  0.0  ,  0.77 ,  0.25 ,  0.0  ,
C    .  0.25 ,  0.25 ,  0.0  ,  0.77 ,  0.25 ,  0.0  ,
C    .  0.77 ,  0.25 ,  0.0  ,  0.77 ,  0.25 ,  0.0  ,
C    .  0.77 ,  0.25 ,  0.0  ,  0.77 ,  0.25 ,  0.0  ,
C    .  0.77 ,  0.25 ,  0.0  ,  0.77 ,  0.25 ,  0.0  ,
C    .  1.0  ,  0.25 ,  0.0  ,  0.0  ,  0.0  ,  0.0  /
C
C     EXCITATION ENERGY (OTTAVIANI) (ONLY FIRST 8 STATES CORRECT)
C      DATA ZE15/
C    .     2.5  ,     5.0 ,     0.0 ,     5.2   ,    10.0 ,     0.0 ,
C    .    10.0  ,    15.0 ,     0.0 ,    18.0   ,    23.5 ,     0.0 ,
C    .    32.0  ,    40.0 ,     0.0 ,    42.0   ,    60.0 ,     0.0 ,
C    .    56.4  ,    82.2 ,     0.0 ,    64.0   ,   108.0 ,     0.0 ,
C    .    72.5  ,   149.5 ,     0.0 ,    35.0   ,    71.0 ,   161.0 ,
C    .    35.0  ,    68.0 ,   171.0 ,    35.0   ,    64.0 ,   188.0 ,
C    .    39.0  ,    62.0 ,   199.0 ,    42.0   ,    59.0 ,   210.0 ,
C    .    43.6  ,   234.0 ,     0.0 ,    36.0   ,   246.0 ,     0.0 ,
C    .   814.0  ,     0.0 ,     0.0 ,   132.0   ,   861.0 ,     0.0 ,
C    .   118.0  ,   919.0 ,     0.0 ,    98.0   ,   966.0 ,     0.0 ,
C    .   105.0  ,  1000.0 ,     0.0 ,    97.5   ,  1054.0 ,     0.0 ,
C    .    93.0  ,  1125.0 ,     0.0 ,    57.5   ,  1167.0 ,     0.0 ,
C    .  6687.0  ,  8000.0 ,     0.0 ,  6892.0   ,  8400.0 ,     0.0 /
C
C     OSCILLATOR STRENGTH (OTTAVIANI) (ONLY FIRST 8 STATES CORRECT)
C      DATA ZF15/
C    .  2.0   ,  3.0   ,  0.0   ,  1.0   ,  3.0   ,  0.0   ,
C    .  1.0   ,  3.0   ,  0.0   ,  1.5   ,  2.5   ,  0.0   ,
C    .  2.0   ,  2.0   ,  0.0   ,  2.0   ,  1.5   ,  0.0   ,
C    .  3.0   ,  2.5   ,  0.0   ,  3.0   ,  2.0   ,  0.0   ,
C    .  2.65  ,  2.0   ,  0.0   ,  0.17  ,  2.5   ,  2.0   ,
C    .  0.3   ,  2.0   ,  2.0   ,  0.5   ,  2.0   ,  2.0   ,
C    .  0.6   ,  1.0   ,  1.5   ,  0.7   ,  0.5   ,  1.0   ,
C    .  0.8   ,  1.5   ,  0.0   ,  0.4   ,  1.0   ,  0.0   ,
C    .  3.7   ,  0 0   ,  0.0   ,  0.05  ,  2.5   ,  0.0   ,
C    .  0.1   ,  2.5   ,  0.0   ,  0.12  ,  2.0   ,  0.0   ,
C    .  0.15  ,  1.5   ,  0.0   ,  0.15  ,  1.0   ,  0.0   ,
C    .  0.15  ,  1.0   ,  0.0   ,  0.07  ,  0.8   ,  0.0   ,
C    .  0.77  ,  0.3   ,  0.0   ,  0.42  ,  0.14  ,  0.0   /
C
C     GAUNT FACTOR (OTTAVIANI) (ONLY FIRST 8 STATES CORRECT)
C      DATA ZG15/
C    .  1.0  , -0.3  ,  0.0  ,  1.0  , -0.3  ,  0.0  ,
C    .  1.0  , -0.3  ,  0.0  ,  1.0  , -0.3  ,  0.0  ,
C    .  1.0  , -0.3  ,  0.0  ,  1.0  , -0.3  ,  0.0  ,
C    .  1.0  , -0.3  ,  0.0  ,  1.0  , -0.3  ,  0.0  ,
C    .  1.0  ,  0.3  ,  0.0  ,  1.0  ,  1.0  ,  0.3  ,
C    .  1.0  ,  1.0  ,  0.3  ,  1.0  ,  1.0  ,  0.3  ,
C    .  1.0  ,  1.0  ,  0.3  ,  1.0  ,  1.0  ,  0.3  ,
C    .  1.0  ,  0.3  ,  0.0  ,  1.0  ,  0.3  ,  0.0  ,
C    .  0.3  ,  0.0  ,  0.0  ,  1.0  ,  0.3  ,  0.0  ,
C    .  1.0  ,  0.3  ,  0.0  ,  1.0  ,  0.3  ,  0.0  ,
C    .  1.0  ,  0.3  ,  0.0  ,  1.0  ,  0.3  ,  0.0  ,
C    .  1.0  ,  0.3  ,  0.0  ,  1.0  ,  0.3  ,  0.0  ,
C    .  0.3  ,  0.3  ,  0.0  ,  0.3  ,  0.3  ,  0.0  /
C
C-----------------------------------------------------------------------
C
CL     1.16 NICKEL
C     IMPURITY DESCRIPTION
       DATA ZMASS(16),ISTATE(16) /59.0 , 28/
       DATA INAM16/ 4HNI I, 4H    , 4H    , 4HNI I, 4HI   , 4H    ,
     .              4HNI I, 4HII  , 4H    , 4HNI I, 4HV   , 4H    ,
     .              4HNI V, 4H    , 4H    , 4HNI V, 4HI   , 4H    ,
     .              4HNI V, 4HII  , 4H    , 4HNI V, 4HIII , 4H    ,
     .              4HNI I, 4HX   , 4H    , 4HNI X, 4H    , 4H    ,
     .              4HNI X, 4HI   , 4H    , 4HNI X, 4HII  , 4H    ,
     .              4HNI X, 4HIII , 4H    , 4HNI X, 4HIV  , 4H    ,
     .              4HNI X, 4HV   , 4H    , 4HNI X, 4HVI  , 4H    ,
     .              4HNI X, 4HVII , 4H    , 4HNI X, 4HVIII, 4H    ,
     .              4HNI X, 4HIX  , 4H    , 4HNI X, 4HX   , 4H    ,
     .              4HNI X, 4HXI  , 4H    , 4HNI X, 4HXII , 4H    ,
     .              4HNI X, 4HXIII, 4H    , 4HNI X, 4HXIV , 4H    ,
     .              4HNI X, 4HXV  , 4H    , 4HNI X, 4HXVI , 4H    ,
     .              4HNI X, 4HXVII, 4H    , 4HNI X, 4HXVII, 4HI   /
C
C-----------------------------------------------------------------------
C
C    IONISATION POTENTIALS FOR FIRST SHELL ARE FROM BURKE ET AL.
C    IONISATION POTENTIALS FOR SECOND AND THIRD SHELL ARE NOT USED
C
C-----------------------------------------------------------------------
C
C     COEFFICIENT P (LOTZ) (FOR 2ND AND 3RD SHELL TFR GROUP, JULY 1984)
       DATA ZP16/
     .     7.638,    10.0 ,    73.0 ,    18.169 ,    97.0 ,    18.169 ,
     .    35.32 ,   122.0 ,   166.0 ,    54.9   ,   146.0 ,   190.0   ,
     .    76.1  ,   171.0 ,   215.0 ,   108.0   ,   196.0 ,   239.0   ,
     .   133.0  ,   221.0 ,   264.0 ,   162.0   ,   246.0 ,   288.0 ,
     .   193.0  ,   271.0 ,   313.0 ,   224.6   ,   296.0 ,   338.0 ,
     .   321.0  ,   363.0 ,  1214.0 ,   352.0   ,   393.0 ,  1253.0 ,
     .   384.0  ,   423.0 ,  1293.0 ,   430.0   ,   458.0 ,   917.0 ,
     .   464.0  ,   494.0 ,   917.0 ,   499.0   ,   531.0 ,   917.0 ,
     .   571.3  ,   917.0 ,  1066.0 ,   607.02  ,   868.0 ,  1017.0 ,
     .  1541.0  ,  1694.0 ,  9275.0 ,  1648.0   ,  1775.0 ,  9381.0 ,
     .  1756.0  ,  1858.0 ,  9488.0 ,  1894.0   ,  1938.0 ,  9595.0 ,
     .  2011.0  ,  2056.0 ,  9702.0 ,  2131.0   ,  2178.0 ,  9808.0 ,
     .  2295.0  ,  9914.0 ,     1.0 ,  2399.0   , 10020.0 ,     1.0 ,
     . 10289.5  ,     1.0 ,     1.0 , 10775.33  ,     1.0 ,     1.0 /
C
C     EFFECTIVE NUMBER OF ELECTRONS IN A SHELL (BURGESS-CHIDICHIMO)
       DATA KSI16/
     .  2 ,  2 ,  6 ,  2 ,  6 , 10 ,  2 ,  2 ,  6 ,  2 ,  6 ,  9 ,
     .  2 ,  2 ,  6 ,  2 ,  6 ,  8 ,  2 ,  2 ,  6 ,  2 ,  6 ,  7 ,
     .  2 ,  2 ,  6 ,  2 , 12 ,  0 ,  2 ,  2 ,  6 , 13 ,  0 ,  0 ,
     .  2 ,  2 ,  6 , 12 ,  0 ,  0 ,  2 ,  2 ,  6 , 11 ,  0 ,  0 ,
     .  2 ,  2 ,  6 , 10 ,  0 ,  0 ,  2 ,  2 ,  6 ,  9 ,  0 ,  0 ,
     .  2 ,  2 ,  6 ,  8 ,  0 ,  0 ,  2 ,  2 ,  6 ,  7 ,  0 ,  0 ,
     .  2 ,  2 ,  6 ,  6 ,  0 ,  0 ,  2 ,  2 ,  6 ,  5 ,  0 ,  0 ,
     .  2 ,  2 ,  6 ,  4 ,  0 ,  0 ,  2 ,  2 ,  6 ,  3 ,  0 ,  0 ,
     .  2 ,  2 ,  6 ,  2 ,  0 ,  0 ,  2 ,  2 ,  6 ,  1 ,  0 ,  0 ,
     .  2 ,  8 ,  0 ,  0 ,  0 ,  0 ,  2 ,  7 ,  0 ,  0 ,  0 ,  0 ,
     .  2 ,  6 ,  0 ,  0 ,  0 ,  0 ,  2 ,  5 ,  0 ,  0 ,  0 ,  0 ,
     .  2 ,  4 ,  0 ,  0 ,  0 ,  0 ,  2 ,  3 ,  0 ,  0 ,  0 ,  0 ,
     .  2 ,  2 ,  0 ,  0 ,  0 ,  0 ,  2 ,  1 ,  0 ,  0 ,  0 ,  0 ,
     .  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,  1 ,  0 ,  0 ,  0 ,  0 ,  0 /
C
C     EFFECTIVE POTENTIAL IN A SHELL (BURGESS-CHIDICHIMO)
       DATA POT16A/
     .  8350.0 , 1020.0 ,  860.0 , 111.0 , 66.0 ,  7.637 ,
     .  8300.0 , 1020.0 ,  860.0 , 110.0 , 66.0 , 18.17  ,
     .  8300.0 , 1030.0 ,  860.0 , 109.0 , 66.0 , 35.32  ,
     .  8300.0 , 1030.0 ,  860.0 , 109.0 , 66.0 , 54.9   ,
     .  8300.0 , 1030.0 ,  860.0 , 108.0 , 76.1 ,  1.0   ,
     .  8300.0 , 1040.0 ,  860.0 , 108.0 ,  1.0 ,  1.0   ,
     .  8300.0 , 1040.0 ,  860.0 , 133.0 ,  1.0 ,  1.0   ,
     .  8300.0 , 1040.0 ,  860.0 , 162.0 ,  1.0 ,  1.0   ,
     .  8300.0 , 1050.0 ,  860.0 , 193.0 ,  1.0 ,  1.0   ,
     .  8300.0 , 1050.0 ,  860.0 , 224.5 ,  1.0 ,  1.0   ,
     .  8300.0 , 1050.0 ,  860.0 , 321.0 ,  1.0 ,  1.0   ,
     .  8300.0 , 1060.0 ,  860.0 , 352.0 ,  1.0 ,  1.0   ,
     .  8300.0 , 1060.0 ,  860.0 , 384.0 ,  1.0 ,  1.0   ,
     .  8300.0 , 1060.0 ,  860.0 , 430.0 ,  1.0 ,  1.0   /
C
       DATA POT16B/
     .  8300.0 , 1070.0 ,  860.0 , 464.0 ,  1.0 ,  1.0   ,
     .  8300.0 , 1070.0 ,  860.0 , 499.0 ,  1.0 ,  1.0   ,
     .  8300.0 , 1070.0 ,  860.0 , 571.0 ,  1.0 ,  1.0   ,
     .  8300.0 , 1029.0 ,  863.0 , 607.0 ,  1.0 ,  1.0   ,
     .  8300.0 , 1541.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .  7800.0 , 1648.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .  7800.0 , 1756.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .  7800.0 , 1894.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .  7800.0 , 2011.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .  7800.0 , 2131.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .  7800.0 , 2295.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     .  7649.0 , 2399.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     . 10290.0 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   ,
     . 10775.0 ,    1.0 ,    1.0 ,   1.0 ,  1.0 ,  1.0   /
C
C     FITTING CONSTANT C (BURGESS-CHIDICHIMO)
       DATA ZCON16/
     .  1.2 , 1.7 , 2.3 , 2.3 , 2.3 , 2.3 , 2.3 , 2.3 , 2.3 , 2.3 ,
     .  2.3 , 2.3 , 2.3 , 2.3 , 2.3 , 2.3 , 2.3 , 2.3 , 2.3 , 2.3 ,
     .  2.3 , 2.3 , 2.3 , 2.3 , 2.3 , 2.3 , 2.3 , 2.3 /
C
C     PRINCIPAL QUANTUM NUMBER NZ OF VALENCE SHELL
       DATA IZ16/
     .  3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 ,
     .  3 , 3 , 3 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 1 , 1 /
C
C     NUMBER OF ELECTRONS IN VALENCE SHELL
       DATA IK16/
     . 17 , 16 , 15 , 14 , 13 , 12 , 11 , 10 ,  9 ,  8 ,
     .  7 ,  6 ,  5 ,  4 ,  3 ,  2 ,  1 ,  0 ,  7 ,  6 ,
     .  5 ,  4 ,  3 ,  2 ,  1 ,  0 ,  1 ,  0 /
C
C     EXCITATION ENERGY (BEHRINGER, 6 SEPTEMBER 1984)
       DATA ZE16/
     .    0.0 ,    0.0 ,    0.0 ,   10.0 ,   15.0 ,    0.0 ,
     .   15.0 ,   25.0 ,    0.0 ,   20.0 ,   30.0 ,    0.0 ,
     .   30.0 ,   45.0 ,    0.0 ,   40.0 ,   60.0 ,    0.0 ,
     .   50.0 ,   70.0 ,    0.0 ,   60.0 ,  100.0 ,    0.0 ,
     .   70.0 ,  121.0 ,    0.0 ,   85.0 ,  137.0 ,    0.0 ,
     .   84.0 ,  177.0 ,    0.0 ,   78.0 ,  165.0 ,    0.0 ,
     .   56.0 ,  172.0 ,    0.0 ,   69.0 ,  188.0 ,    0.0 ,
     .   62.0 ,  200.0 ,    0.0 ,   58.0 ,  250.0 ,    0.0 ,
     .   50.0 ,  288.0 ,    0.0 ,   41.0 ,  302.0 ,    0.0 ,
     .  997.0 , 1233.0 ,    0.0 ,  149.0 , 1042.0 ,    0.0 ,
     .  129.0 , 1033.0 ,    0.0 ,  112.0 , 1100.0 ,    0.0 ,
     .  119.0 , 1127.0 ,    0.0 ,  108.0 , 1200.0 ,    0.0 ,
     .  105.0 , 1327.0 ,    0.0 ,   67.7 , 1367.0 ,    0.0 ,
     . 7804.0 , 9170.0 ,    0.0 , 8049.0 , 9578.0 ,    0.0 /
C
C     OSCILLATOR STRENGTH (BEHRINGER, 6 SEPTEMBER 1984)
       DATA ZF16/
     .  0.0  ,  0.0  ,  0.0  ,  2.0  ,  1.0  ,  0.0  ,
     .  2.0  ,  1.0  ,  0.0  ,  2.0  ,  1.0  ,  0.0  ,
     .  2.0  ,  1.0  ,  0.0  ,  2.0  ,  1.0  ,  0.0  ,
     .  2.0  ,  1.0  ,  0.0  ,  2.0  ,  1.0  ,  0.0  ,
     .  2.0  ,  1.1  ,  0.0  ,  2.0  ,  1.3  ,  0.0  ,
     .  2.3  ,  1.1  ,  0.0  ,  2.3  ,  0.8  ,  0.0  ,
     .  1.1  ,  0.4  ,  0.0  ,  1.7  ,  0.3  ,  0.0  ,
     .  2.0  ,  0.4  ,  0.0  ,  1.1  ,  0.4  ,  0.0  ,
     .  0.77 ,  0.45 ,  0.0  ,  0.37 ,  0.3  ,  0.0  ,
     .  3.5  ,  1.1  ,  0.0  ,  0.06 ,  1.3  ,  0.0  ,
     .  0.09 ,  1.6  ,  0.0  ,  0.11 ,  1.6  ,  0.0  ,
     .  0.16 ,  1.4  ,  0.0  ,  0.21 ,  1.0  ,  0.0  ,
     .  0.15 ,  0.87 ,  0.0  ,  0.07 ,  0.56 ,  0.0  ,
     .  0.77 ,  0.24 ,  0.0  ,  0.42 ,  0.15 ,  0.0  /
C
C     GAUNT FACTOR (BEHRINGER, 6 SEPTEMBER 1984)
       DATA ZG16/
     . 0.0  , 0.0  , 0.0  , 0.8  ,-0.25 , 0.0  ,
     . 0.8  ,-0.25 , 0.0  , 0.8  ,-0.25 , 0.0  ,
     . 0.8  ,-0.25 , 0.0  , 0.8  ,-0.25 , 0.0  ,
     . 0.8  ,-0.25 , 0.0  , 0.8  ,-0.25 , 0.0  ,
     . 0.8  ,-0.25 , 0.0  , 0.8  ,-0.25 , 0.0  ,
     . 0.8  ,-0.25 , 0.0  , 0.8  ,-0.25 , 0.0  ,
     . 0.8  ,-0.25 , 0.0  , 0.8  ,-0.25 , 0.0  ,
     . 0.8  ,-0.25 , 0.0  , 0.8  ,-0.25 , 0.0  ,
     . 0.8  ,-0.25 , 0.0  , 0.8  ,-0.25 , 0.0  ,
     .-0.25 ,-0.25 , 0.0  , 0.8  ,-0.25 , 0.0  ,
     . 0.8  ,-0.25 , 0.0  , 0.8  ,-0.25 , 0.0  ,
     . 0.8  ,-0.25 , 0.0  , 0.8  ,-0.25 , 0.0  ,
     . 0.8  ,-0.25 , 0.0  , 0.8  ,-0.25 , 0.0  ,
     .-0.8  ,-0.25 , 0.0  ,-0.42 ,-0.35 , 0.0  /
C
C     COEFFICIENT A (LOTZ) (BEHRINGER, JULY 1984)
C      DATA ZA16/
C    . 4.0 , 2.4 , 0.0 , 3.0 , 0.0 , 0.0 , 3.8 , 4.5 , 0.0 ,
C    . 4.4 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
C    . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
C    . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
C    . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
C    . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
C    . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
C    . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 ,
C    . 4.5 , 4.5 , 0.0 , 4.5 , 4.5 , 0.0 , 4.5 , 0.0 , 0.0 ,
C    . 4.5 , 0.0 , 0.0 /
C
C     COEFFICIENT A (LOTZ) (TFR GROUP, JULY 1984)
C      DATA ZA16/
C    . 4.0 , 2.4 , 4.0 , 3.0 , 4.4 , 0.0 , 3.8 , 4.5 , 4.5 ,
C    . 4.4 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 ,
C    . 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 ,
C    . 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 ,
C    . 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 ,
C    . 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 ,
C    . 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 ,
C    . 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 4.5 , 3.375,3.375,3.375,
C    . 3.375,3.375,0.0 , 3.375,3.375,0.0 , 4.5 , 0.0 , 0.0 ,
C    . 4.5 , 0.0 , 0.0 /
C
C     COEFFICIENT B (LOTZ) (BEHRINGER, JULY 1984)
C      DATA ZB16/
C    . 0.4  , 0.94 , 0.0  , 0.9  , 0.0  , 0.0  , 0.7  , 0.2  , 0.0  ,
C    . 0.4  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT B (LOTZ) (TFR GROUP, JULY 1984)
C      DATA ZB16/
C    . 0.4  , 0.94 , 0.6  , 0.9  , 0.3  , 0.0  , 0.7  , 0.2  , 0.0  ,
C    . 0.4  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT C (LOTZ) (BEHRINGER, JULY 1984)
C      DATA ZC16/
C    . 0.6  , 0.17 , 0.0  , 0.25 , 0.0  , 0.0  , 0.5  , 0.6  , 0.0  ,
C    . 0.6  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT C (LOTZ) (TFR GROUP, JULY 1984)
C      DATA ZC16/
C    . 0.6  , 0.17 , 0.4  , 0.25 , 0.6  , 0.0  , 0.5  , 0.6  , 0.0  ,
C    . 0.6  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
C    . 0.0  , 0.0  , 0.0  /
C
C     COEFFICIENT P (LOTZ) (BEHRINGER, JULY 1984)
C      DATA ZP16/
C    .     8.68 ,    10.0 ,     1.0 ,    18.2   ,    18.2 ,     1.0 ,
C    .    35.2  ,   122.0 ,     1.0 ,    54.9   ,   146.0 ,     1.0 ,
C    .    75.5  ,   171.0 ,     1.0 ,   108.0   ,   196.0 ,     1.0 ,
C    .   133.0  ,   221.0 ,     1.0 ,   162.0   ,   246.0 ,     1.0 ,
C    .   193.0  ,   271.0 ,     1.0 ,   225.0   ,   296.0 ,     1.0 ,
C    .   321.0  ,   363.0 ,     1.0 ,   352.0   ,   393.0 ,     1.0 ,
C    .   384.0  ,   423.0 ,     1.0 ,   430.0   ,   458.0 ,     1.0 ,
C    .   464.0  ,   494.0 ,     1.0 ,   499.0   ,   531.0 ,     1.0 ,
C    .   571.0  ,  1458.0 ,     1.0 ,   608.0   ,  1500.0 ,     1.0 ,
C    .  1546.0  ,  1694.0 ,     1.0 ,  1648.0   ,  1775.0 ,     1.0 ,
C    .  1756.0  ,  1858.0 ,     1.0 ,  1894.0   ,  1938.0 ,     1.0 ,
C    .  2011.0  ,  2056.0 ,     1.0 ,  2131.0   ,  2178.0 ,     1.0 ,
C    .  2295.0  ,  9914.0 ,     1.0 ,  2399.0   , 10020.0 ,     1.0 ,
C    . 10290.0  ,     1.0 ,     1.0 , 10775.0   ,     1.0 ,     1.0 /
C
C     COEFFICIENT Q (LOTZ) (BEHRINGER, JULY 1984)
C      DATA IQ16/
C    . 2 , 8 , 0 , 9 , 0 , 0 , 8 , 6 , 0 , 7 , 6 , 0 , 6 , 6 , 0 ,
C    . 5 , 6 , 0 , 4 , 6 , 0 , 3 , 6 , 0 , 2 , 6 , 0 , 1 , 6 , 0 ,
C    . 6 , 2 , 0 , 5 , 2 , 0 , 4 , 2 , 0 , 3 , 2 , 0 , 2 , 2 , 0 ,
C    . 1 , 2 , 0 , 2 , 6 , 0 , 1 , 6 , 0 , 6 , 2 , 0 , 5 , 2 , 0 ,
C    . 4 , 2 , 0 , 3 , 2 , 0 , 2 , 2 , 0 , 1 , 2 , 0 , 2 , 2 , 0 ,
C    . 1 , 2 , 0 , 2 , 0 , 0 , 1 , 0 , 0 /
C
C     COEFFICIENT Q (LOTZ) (TFR GROUP, JULY 1984)
C      DATA IQ16/
C    . 2 , 8 , 6 , 9 , 6 , 0 , 8 , 6 , 2 , 7 , 6 , 2 , 6 , 6 , 2 ,
C    . 5 , 6 , 2 , 4 , 6 , 2 , 3 , 6 , 2 , 2 , 6 , 2 , 1 , 6 , 2 ,
C    . 6 , 2 , 6 , 5 , 2 , 6 , 4 , 2 , 6 , 3 , 2 , 8 , 2 , 2 , 8 ,
C    . 1 , 2 , 8 , 2 , 6 , 2 , 1 , 6 , 2 , 6 , 2 , 2 , 5 , 2 , 2 ,
C    . 4 , 2 , 2 , 3 , 2 , 2 , 2 , 2 , 2 , 1 , 2 , 2 , 2 , 2 , 0 ,
C    . 1 , 2 , 0 , 2 , 0 , 0 , 1 , 0 , 0 /
C
C     EXCITATION ENERGY (TFR, 6 SEPTEMBER 1984) (TESTED)
C      DATA ZE16/
C    .    0.0 ,    0.0 ,    0.0 ,   17.0 ,    9.0 ,    0.0 ,
C    .   12.0 ,   16.0 ,    0.0 ,   35.0 ,   22.0 ,    0.0 ,
C    .   40.0 ,   35.0 ,    0.0 ,   40.0 ,   50.0 ,    0.0 ,
C    .   40.0 ,   60.0 ,    0.0 ,   45.0 ,   70.0 ,    0.0 ,
C    .   50.0 ,  100.0 ,    0.0 ,   85.0 ,  130.0 ,    0.0 ,
C    .   84.0 ,  200.0 ,    0.0 ,   80.0 ,   41.0 ,    0.0 ,
C    .   80.0 ,   41.0 ,    0.0 ,   73.0 ,   40.0 ,    0.0 ,
C    .   70.0 ,   40.0 ,    0.0 ,   65.0 ,   40.0 ,    0.0 ,
C    .   50.0 ,  290.0 ,    0.0 ,   40.0 ,  300.0 ,    0.0 ,
C    .    0.0 , 1000.0 ,    0.0 ,  140.0 , 1050.0 ,    0.0 ,
C    .  125.0 , 1100.0 ,    0.0 ,  115.0 , 1150.0 ,    0.0 ,
C    .  120.0 , 1250.0 ,    0.0 ,  110.0 , 1200.0 ,    0.0 ,
C    .  102.0 , 1300.0 ,    0.0 ,   65.0 , 1350.0 ,    0.0 ,
C    .    0.0 , 7750.0 ,    0.0 ,    0.0 , 8000.0 ,    0.0 /
C
C     OSCILLATOR STRENGTH (TFR, 6 SEPTEMBER 1984) (TESTED)
C      DATA ZF16/
C    .  0.0  ,  0.0  ,  0.0  ,  0.4  ,  0.23 ,  0.0  ,
C    .  0.7  ,  0.29 ,  0.0  ,  1.0  ,  0.35 ,  0.0  ,
C    .  1.4  ,  0.33 ,  0.0  ,  1.5  ,  0.32 ,  0.0  ,
C    .  1.7  ,  0.34 ,  0.0  ,  1.9  ,  0.36 ,  0.0  ,
C    .  2.0  ,  0.3  ,  0.0  ,  2.15 ,  0.275,  0.0  ,
C    .  3.4  ,  0.2  ,  0.0  ,  2.6  ,  0.155,  0.0  ,
C    .  2.0  ,  0.28 ,  0.0  ,  1.8  ,  0.43 ,  0.0  ,
C    .  0.78 ,  0.55 ,  0.0  ,  3.2  ,  0.66 ,  0.0  ,
C    .  0.75 ,  0.45 ,  0.0  ,  0.35 ,  0.25 ,  0.0  ,
C    .  0.0  ,  2.4  ,  0.0  ,  0.055,  5.0  ,  0.0  ,
C    .  0.095,  1.6  ,  0.0  ,  0.11 ,  1.6  ,  0.0  ,
C    .  0.14 ,  1.2  ,  0.0  ,  0.14 ,  0.675,  0.0  ,
C    .  0.15 ,  0.6  ,  0.0  ,  0.065,  0.38 ,  0.0  ,
C    .  0.0  ,  0.79 ,  0.0  ,  0.0  ,  0.42 ,  0.0  /
C
C     GAUNT FACTOR (TFR, 6 SEPTEMBER 1984) (ASK ELISABETH KALLNE)
C      DATA ZG16/
C    . 0.0  , 0.0  , 0.0  , 0.8  ,-0.25 , 0.0  ,
C    . 0.8  ,-0.25 , 0.0  , 0.8  ,-0.25 , 0.0  ,
C    . 0.8  ,-0.25 , 0.0  , 0.8  ,-0.25 , 0.0  ,
C    . 0.8  ,-0.25 , 0.0  , 0.8  ,-0.25 , 0.0  ,
C    . 0.8  ,-0.25 , 0.0  , 0.8  ,-0.25 , 0.0  ,
C    . 0.8  ,-0.25 , 0.0  , 0.8  ,-0.25 , 0.0  ,
C    . 0.8  ,-0.25 , 0.0  , 0.8  ,-0.25 , 0.0  ,
C    . 0.8  ,-0.25 , 0.0  , 0.8  ,-0.25 , 0.0  ,
C    . 0.8  ,-0.25 , 0.0  , 0.8  ,-0.25 , 0.0  ,
C    . 0.0  ,-0.25 , 0.0  , 0.8  ,-0.25 , 0.0  ,
C    . 0.8  ,-0.25 , 0.0  , 0.8  ,-0.25 , 0.0  ,
C    . 0.8  ,-0.25 , 0.0  , 0.8  ,-0.25 , 0.0  ,
C    . 0.8  ,-0.25 , 0.0  , 0.8  ,-0.25 , 0.0  ,
C    . 0.0  ,-0.25 , 0.0  , 0.0  ,-0.35 , 0.0  /
C
C-----------------------------------------------------------------------
c
c      Local variables:
c
       integer j3,idial,j2,j1,j213,j214,j,istat,istatm,istatp
       integer j2sq,inz,inz1,j0,int,jj,j2p1,j2p12,j2p13,j2p1sq
       integer imid,imid2,idn,istar,ids
c
       real zexp1,zexp2,zexp3,zq,zbeta,zte,zsa,zpot,z,zinv
       real ze,zw,zpp,zcc,zpc,zs,zt,zterm,zmu,zcns,zeht,zppl
       real zs1,zehj2,zcons,zne,zsum,zpp2,zdc,zeht15,zdb
       real zpp1,zehj2p,zde1,zde2,zde3,zy1,zy2,zy3,zgg1,zgg2,zgg3
       real zza11,zza12,zza13,zza21,zza22,zza23,zaa1,zaa2,zaa3
       real zd1,zd2,zd3,zcor,zzn,zn1,zn2,zn3,zero,zline,zpe,zpeh
       real zdel,zd,zfac,zsqrte,zpelog,zdie,z2,zgaunt
c
C-----------------------------------------------------------------------
CL              2.     DETERMINE DIALLING ORDER
C
C   SCAN OVER MAXIMUM OF THREE IMPURITIES
         NIMP=MIN0(3,KIMP)
         DO 499 J3=1,NIMP
         IDIAL=0
         KFAIL(J3)=1
C
C   FIND MATCH WITH INTERNAL LIST OF MASSES
         DO 210 J2=1,16
         IF ( ABS(PMASS(J3)-ZMASS(J2)) .GT. 1.0E-05 ) GO TO 210
         IDIAL=J2
         KFAIL(J3)=0
 210     CONTINUE
C
C-----------------------------------------------------------------------
CL              3.     FILL COMMON BLOCKS IN RIGHT ORDER
C
C   FILL ARRAYS FOR MATCHING IMPURITY ONLY
         NORDER(J3)=IDIAL
         IF ( IDIAL .EQ. 0 ) GO TO 499
         NSTATE(J3)=ISTATE(IDIAL)
C
C     GET VALUES FOR THE DIALLED IMPURITY
         GO TO  (310,315,320,325,330,335,340,345,350,355,
     .           360,365,375,380,385,395), IDIAL
C
C     ELEMENT J3 IS HYDROGEN
 310     CONTINUE
         DO 312 J2=1,NSTATE(J3)
         DO 311 J1=1,3
         ZA(J1,J2,J3)   =ZA1(J1,J2)
         ZB(J1,J2,J3)   =ZB1(J1,J2)
         ZC(J1,J2,J3)   =ZC1(J1,J2)
         ZP(J1,J2,J3)   =ZP1(J1,J2)
         IQ(J1,J2,J3)   =IQ1(J1,J2)
         DE(J1,J2,J3)   =ZE1(J1,J2)
         F(J1,J2,J3)    =ZF1(J1,J2)
         G(J1,J2,J3)    =ZG1(J1,J2)
         INAME(J1,J2,J3)=INAM1(J1,J2)
 311     CONTINUE
         IZ(J2,J3)=IZ1(J2)
         IK(J2,J3)=IK1(J2)
         P(J2,J3) =ZP1(1,J2)
 312     CONTINUE
         GO TO 499
C
C     ELEMENT J3 IS DEUTERIUM
 315     CONTINUE
         DO 317 J2=1,NSTATE(J3)
         DO 316 J1=1,3
         ZA(J1,J2,J3)   =ZA1(J1,J2)
         ZB(J1,J2,J3)   =ZB1(J1,J2)
         ZC(J1,J2,J3)   =ZC1(J1,J2)
         ZP(J1,J2,J3)   =ZP1(J1,J2)
         IQ(J1,J2,J3)   =IQ1(J1,J2)
         DE(J1,J2,J3)   =ZE1(J1,J2)
         F(J1,J2,J3)    =ZF1(J1,J2)
         G(J1,J2,J3)    =ZG1(J1,J2)
         INAME(J1,J2,J3)=INAM2(J1,J2)
 316     CONTINUE
         IZ(J2,J3)=IZ1(J2)
         IK(J2,J3)=IK1(J2)
         P(J2,J3) =ZP1(1,J2)
 317     CONTINUE
         GO TO 499
C
C     ELEMENT J3 IS HELIUM
 320     CONTINUE
         DO 322 J2=1,NSTATE(J3)
         DO 321 J1=1,3
         ZA(J1,J2,J3)   =ZA3(J1,J2)
         ZB(J1,J2,J3)   =ZB3(J1,J2)
         ZC(J1,J2,J3)   =ZC3(J1,J2)
         ZP(J1,J2,J3)   =ZP3(J1,J2)
         IQ(J1,J2,J3)   =IQ3(J1,J2)
         DE(J1,J2,J3)   =ZE3(J1,J2)
         F(J1,J2,J3)    =ZF3(J1,J2)
         G(J1,J2,J3)    =ZG3(J1,J2)
         INAME(J1,J2,J3)=INAM3(J1,J2)
 321     CONTINUE
         IZ(J2,J3)=IZ3(J2)
         IK(J2,J3)=IK3(J2)
         P(J2,J3) =ZP3(1,J2)
 322     CONTINUE
         GO TO 499
C
C     ELEMENT J3 IS BERYLLIUM
 325     CONTINUE
         DO 328 J2=1,NSTATE(J3)
         DO 326 J1=1,3
         ZA(J1,J2,J3)   =ZA4(J1,J2)
         ZB(J1,J2,J3)   =ZB4(J1,J2)
         ZC(J1,J2,J3)   =ZC4(J1,J2)
         ZP(J1,J2,J3)   =ZP4(J1,J2)
         IQ(J1,J2,J3)   =IQ4(J1,J2)
         DE(J1,J2,J3)   =ZE4(J1,J2)
         F(J1,J2,J3)    =ZF4(J1,J2)
         G(J1,J2,J3)    =ZG4(J1,J2)
         INAME(J1,J2,J3)=INAM4(J1,J2)
 326     CONTINUE
C        DO 327 J1=1,10
C        ZDE(J1,J2,J3) =ZDE4(J1,J2)
C        ZBE(J1,J2,J3) =ZBE4(J1,J2)
C        ZCHI(J1,J2,J3)=ZCHI4(J1,J2)
C327     CONTINUE
         IZ(J2,J3)=IZ4(J2)
         IK(J2,J3)=IK4(J2)
         P(J2,J3) =ZP4(1,J2)
 328     CONTINUE
         GO TO 499
C
C     ELEMENT J3 IS CARBON
 330     CONTINUE
         DO 333 J2=1,NSTATE(J3)
         DO 331 J1=1,3
         ZA(J1,J2,J3)   =ZA5(J1,J2)
         ZB(J1,J2,J3)   =ZB5(J1,J2)
         ZC(J1,J2,J3)   =ZC5(J1,J2)
         ZP(J1,J2,J3)   =ZP5(J1,J2)
         IQ(J1,J2,J3)   =IQ5(J1,J2)
         DE(J1,J2,J3)   =ZE5(J1,J2)
         F(J1,J2,J3)    =ZF5(J1,J2)
         G(J1,J2,J3)    =ZG5(J1,J2)
         INAME(J1,J2,J3)=INAM5(J1,J2)
 331     CONTINUE
         DO 332 J1=1,6
         KSI(J1,J2,J3)=KSI5(J1,J2)
         POT(J1,J2,J3)=POT5(J1,J2)
 332     CONTINUE
         IZ(J2,J3)  =IZ5(J2)
         IK(J2,J3)  =IK5(J2)
         P(J2,J3)   =ZP5(1,J2)
         ZCON(J2,J3)=ZCON5(J2)
 333     CONTINUE
         GO TO 499
C
C     ELEMENT J3 IS NITROGEN
 335     CONTINUE
         DO 339 J2=1,NSTATE(J3)
         DO 336 J1=1,3
         ZA(J1,J2,J3)   =ZA6(J1,J2)
         ZB(J1,J2,J3)   =ZB6(J1,J2)
         ZC(J1,J2,J3)   =ZC6(J1,J2)
         ZP(J1,J2,J3)   =ZP6(J1,J2)
         IQ(J1,J2,J3)   =IQ6(J1,J2)
         DE(J1,J2,J3)   =ZE6(J1,J2)
         F(J1,J2,J3)    =ZF6(J1,J2)
         G(J1,J2,J3)    =ZG6(J1,J2)
         INAME(J1,J2,J3)=INAM6(J1,J2)
 336     CONTINUE
         DO 337 J1=1,6
         KSI(J1,J2,J3)=KSI6(J1,J2)
         POT(J1,J2,J3)=POT6(J1,J2)
         ZCON(J2,J3)=ZCON6(J2)
 337     CONTINUE
         DO 338 J1=1,10
         ZDE(J1,J2,J3) =ZDE6(J1,J2)
         ZBE(J1,J2,J3) =ZBE6(J1,J2)
         ZCHI(J1,J2,J3)=ZCHI6(J1,J2)
 338     CONTINUE
         IZ(J2,J3)=IZ6(J2)
         IK(J2,J3)=IK6(J2)
         P(J2,J3) =ZP6(1,J2)
 339     CONTINUE
         GO TO 499
C
C     ELEMENT J3 IS OXYGEN
 340     CONTINUE
         DO 343 J2=1,NSTATE(J3)
         DO 341 J1=1,3
         ZA(J1,J2,J3)   =ZA7(J1,J2)
         ZB(J1,J2,J3)   =ZB7(J1,J2)
         ZC(J1,J2,J3)   =ZC7(J1,J2)
         ZP(J1,J2,J3)   =ZP7(J1,J2)
         IQ(J1,J2,J3)   =IQ7(J1,J2)
         DE(J1,J2,J3)   =ZE7(J1,J2)
         F(J1,J2,J3)    =ZF7(J1,J2)
         G(J1,J2,J3)    =ZG7(J1,J2)
         INAME(J1,J2,J3)=INAM7(J1,J2)
 341     CONTINUE
         DO 342 J1=1,6
         KSI(J1,J2,J3)=KSI7(J1,J2)
         POT(J1,J2,J3)=POT7(J1,J2)
 342     CONTINUE
         IZ(J2,J3)=IZ7(J2)
         IK(J2,J3)=IK7(J2)
         P(J2,J3) =ZP7(1,J2)
         ZCON(J2,J3)=ZCON7(J2)
 343     CONTINUE
         GO TO 499
C
C     ELEMENT J3 IS NEON
 345     CONTINUE
         DO 347 J2=1,NSTATE(J3)
         DO 346 J1=1,3
         ZA(J1,J2,J3)   =ZA8(J1,J2)
         ZB(J1,J2,J3)   =ZB8(J1,J2)
         ZC(J1,J2,J3)   =ZC8(J1,J2)
         ZP(J1,J2,J3)   =ZP8(J1,J2)
         IQ(J1,J2,J3)   =IQ8(J1,J2)
         DE(J1,J2,J3)   =ZE8(J1,J2)
         F(J1,J2,J3)    =ZF8(J1,J2)
         G(J1,J2,J3)    =ZG8(J1,J2)
         INAME(J1,J2,J3)=INAM8(J1,J2)
 346     CONTINUE
         IZ(J2,J3)=IZ8(J2)
         IK(J2,J3)=IK8(J2)
         P(J2,J3) =ZP8(1,J2)
 347     CONTINUE
C        DO 349 J2=1,5
C        DO 348 J1=1,10
C        ZDE(J1,J2,J3) =ZDE8A(J1,J2)
C        ZBE(J1,J2,J3) =ZBE8A(J1,J2)
C        ZCHI(J1,J2,J3)=ZCHI8A(J1,J2)
C348     CONTINUE
C349     CONTINUE
C        DO 3492 J2=6,10
C        DO 3491 J1=1,10
C        J25=J2-5
C        ZDE(J1,J2,J3) =ZDE8B(J1,J25)
C        ZBE(J1,J2,J3) =ZBE8B(J1,J25)
C        ZCHI(J1,J2,J3)=ZCHI8B(J1,J25)
C3491    CONTINUE
C3492    CONTINUE
         GO TO 499
C
C     ELEMENT J3 IS ALUMINIUM
 350     CONTINUE
         DO 352 J2=1,NSTATE(J3)
         DO 351 J1=1,3
         ZA(J1,J2,J3)   =ZA9(J1,J2)
         ZB(J1,J2,J3)   =ZB9(J1,J2)
         ZC(J1,J2,J3)   =ZC9(J1,J2)
         ZP(J1,J2,J3)   =ZP9(J1,J2)
         IQ(J1,J2,J3)   =IQ9(J1,J2)
         DE(J1,J2,J3)   =ZE9(J1,J2)
         F(J1,J2,J3)    =ZF9(J1,J2)
         G(J1,J2,J3)    =ZG9(J1,J2)
         INAME(J1,J2,J3)=INAM9(J1,J2)
 351     CONTINUE
         IZ(J2,J3)=IZ9(J2)
         IK(J2,J3)=IK9(J2)
         P(J2,J3) =ZP9(1,J2)
 352     CONTINUE
         GO TO 499
C
C     ELEMENT J3 IS SILICON
 355     CONTINUE
         DO 357 J2=1,NSTATE(J3)
         DO 356 J1=1,3
         ZA(J1,J2,J3)   =ZA10(J1,J2)
         ZB(J1,J2,J3)   =ZB10(J1,J2)
         ZC(J1,J2,J3)   =ZC10(J1,J2)
         ZP(J1,J2,J3)   =ZP10(J1,J2)
         IQ(J1,J2,J3)   =IQ10(J1,J2)
         DE(J1,J2,J3)   =ZE10(J1,J2)
         F(J1,J2,J3)    =ZF10(J1,J2)
         G(J1,J2,J3)    =ZG10(J1,J2)
         INAME(J1,J2,J3)=INAM10(J1,J2)
 356     CONTINUE
         IZ(J2,J3)=IZ10(J2)
         IK(J2,J3)=IK10(J2)
         P(J2,J3) =ZP10(1,J2)
 357     CONTINUE
         GO TO 499
C
C     ELEMENT J3 IS CHLORINE
 360     CONTINUE
         DO 362 J2=1,NSTATE(J3)
         DO 361 J1=1,3
         ZA(J1,J2,J3)   =ZA11(J1,J2)
         ZB(J1,J2,J3)   =ZB11(J1,J2)
         ZC(J1,J2,J3)   =ZC11(J1,J2)
         ZP(J1,J2,J3)   =ZP11(J1,J2)
         IQ(J1,J2,J3)   =IQ11(J1,J2)
         DE(J1,J2,J3)   =ZE11(J1,J2)
         F(J1,J2,J3)    =ZF11(J1,J2)
         G(J1,J2,J3)    =ZG11(J1,J2)
         INAME(J1,J2,J3)=INAM11(J1,J2)
 361     CONTINUE
         IZ(J2,J3)=IZ11(J2)
         IK(J2,J3)=IK11(J2)
         P(J2,J3) =ZP11(1,J2)
 362     CONTINUE
         GO TO 499
C
C     ELEMENT J3 IS ARGON
 365     CONTINUE
         DO 368 J2=1,NSTATE(J3)
         DO 366 J1=1,3
         ZA(J1,J2,J3)   =ZA12(J1,J2)
         ZB(J1,J2,J3)   =ZB12(J1,J2)
         ZC(J1,J2,J3)   =ZC12(J1,J2)
         ZP(J1,J2,J3)   =ZP12(J1,J2)
         IQ(J1,J2,J3)   =IQ12(J1,J2)
         DE(J1,J2,J3)   =ZE12(J1,J2)
         F(J1,J2,J3)    =ZF12(J1,J2)
         G(J1,J2,J3)    =ZG12(J1,J2)
         INAME(J1,J2,J3)=INAM12(J1,J2)
 366     CONTINUE
         DO 367 J1=1,6
         KSI(J1,J2,J3)=KSI12(J1,J2)
         POT(J1,J2,J3)=POT12(J1,J2)
 367     CONTINUE
         IZ(J2,J3)=IZ12(J2)
         IK(J2,J3)=IK12(J2)
         P(J2,J3) =ZP12(1,J2)
         ZCON(J2,J3)=ZCON12(J2)
 368     CONTINUE
         GO TO 499
C
C     ELEMENT J3 IS TITANIUM
 375     CONTINUE
         DO 377 J2=1,NSTATE(J3)
         DO 376 J1=1,3
         ZA(J1,J2,J3)   =ZA13(J1,J2)
         ZB(J1,J2,J3)   =ZB13(J1,J2)
         ZC(J1,J2,J3)   =ZC13(J1,J2)
         ZP(J1,J2,J3)   =ZP13(J1,J2)
         IQ(J1,J2,J3)   =IQ13(J1,J2)
         DE(J1,J2,J3)   =ZE13(J1,J2)
         F(J1,J2,J3)    =ZF13(J1,J2)
         G(J1,J2,J3)    =ZG13(J1,J2)
         INAME(J1,J2,J3)=INAM13(J1,J2)
 376     CONTINUE
         IZ(J2,J3)=IZ13(J2)
         IK(J2,J3)=IK13(J2)
         P(J2,J3) =ZP13(1,J2)
 377     CONTINUE
         GO TO 499
C
C     ELEMENT J3 IS CHROMIUM
 380     CONTINUE
         DO 382 J2=1,NSTATE(J3)
         DO 381 J1=1,3
         ZA(J1,J2,J3)   =ZA14(J1,J2)
         ZB(J1,J2,J3)   =ZB14(J1,J2)
         ZC(J1,J2,J3)   =ZC14(J1,J2)
         ZP(J1,J2,J3)   =ZP14(J1,J2)
         IQ(J1,J2,J3)   =IQ14(J1,J2)
         DE(J1,J2,J3)   =ZE14(J1,J2)
         F(J1,J2,J3)    =ZF14(J1,J2)
         G(J1,J2,J3)    =ZG14(J1,J2)
         INAME(J1,J2,J3)=INAM14(J1,J2)
 381     CONTINUE
         IZ(J2,J3)=IZ14(J2)
         IK(J2,J3)=IK14(J2)
         P(J2,J3) =ZP14(1,J2)
 382     CONTINUE
         GO TO 499
C
C     ELEMENT J3 IS IRON
 385     CONTINUE
         DO 388 J2=1,NSTATE(J3)
         DO 386 J1=1,3
         ZA(J1,J2,J3)   =ZA15(J1,J2)
         ZB(J1,J2,J3)   =ZB15(J1,J2)
         ZC(J1,J2,J3)   =ZC15(J1,J2)
         ZP(J1,J2,J3)   =ZP15(J1,J2)
         IQ(J1,J2,J3)   =IQ15(J1,J2)
         DE(J1,J2,J3)   =ZE15(J1,J2)
         F(J1,J2,J3)    =ZF15(J1,J2)
         G(J1,J2,J3)    =ZG15(J1,J2)
         INAME(J1,J2,J3)=INAM15(J1,J2)
 386     CONTINUE
         DO 387 J1=1,6
         KSI(J1,J2,J3)  =KSI15(J1,J2)
 387     CONTINUE
         IZ(J2,J3)=IZ15(J2)
         IK(J2,J3)=IK15(J2)
         P(J2,J3) =ZP15(1,J2)
         ZCON(J2,J3)=ZCON15(J2)
 388     CONTINUE
         DO 390 J2=1,13
         DO 389 J1=1,6
         POT(J1,J2,J3)=POT15A(J1,J2)
 389     CONTINUE
 390     CONTINUE
         DO 392 J2=14,26
         J213=J2-13
         DO 391 J1=1,6
         POT(J1,J2,J3)=POT15B(J1,J213)
 391     CONTINUE
 392     CONTINUE
         GO TO 499
C
C     ELEMENT J3 IS NICKEL
 395     CONTINUE
         DO 398 J2=1,NSTATE(J3)
         DO 396 J1=1,3
         ZP(J1,J2,J3)   =ZP16(J1,J2)
         DE(J1,J2,J3)   =ZE16(J1,J2)
         F(J1,J2,J3)    =ZF16(J1,J2)
         G(J1,J2,J3)    =ZG16(J1,J2)
         INAME(J1,J2,J3)=INAM16(J1,J2)
 396     CONTINUE
         DO 397 J1=1,6
         KSI(J1,J2,J3)  =KSI16(J1,J2)
 397     CONTINUE
         IZ(J2,J3)  =IZ16(J2)
         IK(J2,J3)  =IK16(J2)
         P(J2,J3)   =ZP16(1,J2)
         ZCON(J2,J3)=ZCON16(J2)
 398     CONTINUE
         DO 400 J2=1,14
         DO 399 J1=1,6
         POT(J1,J2,J3)=POT16A(J1,J2)
 399     CONTINUE
 400     CONTINUE
         DO 402 J2=15,28
         J214=J2-14
         DO 401 J1=1,6
         POT(J1,J2,J3)=POT16B(J1,J214)
 401     CONTINUE
 402     CONTINUE
 499     CONTINUE
C
C-----------------------------------------------------------------------
CL              4.     INITIALISE TEMPERATURE AND DENSITY ARRAYS
C
C     CONSTANTS
         ZEXP1= 1.0/17.0
         ZEXP2= 2.0/17.0
         ZEXP3=12.0/17.0
C
C     DATA FOR LOGARITHMIC INTERPOLATION OF TEMPERATURE
         DTLOG=0.15
         TMINLG=0.0
         TMAXLG=5.0
         NTE=(TMAXLG-TMINLG)/DTLOG+1.5
C
C     DATA FOR LOGARITHMIC INTERPOLATION OF DENSITY
         DNLOG= 0.3
         DMINLG=10.0
         DMAXLG=15.0
         NNE=(DMAXLG-DMINLG)/DNLOG+1.5
C
C     LOGARITHMIC  AND REAL TEMPERATURE ARRAY
         DO 410 J=1,NTE
         TELOG(J)=TMINLG+(J-1)*DTLOG
         TE(J)=10.0**TELOG(J)
 410     CONTINUE
C
C     LOGARITHMIC AND REAL DENSITY ARRAY
         DO 420 J=1,NNE
         DENLOG(J)=DMINLG+(J-1)*DNLOG
         DENE(J)=10.0**DENLOG(J)
 420     CONTINUE
C
C-----------------------------------------------------------------------
CL              5.     RATE COEFFICIENTS FOR IONISATION
C
C     LOOP OVER IMPURITY SPECIES
         DO 1230 J3=1,NIMP
         IF ( NORDER(J3) .EQ. 0 ) GO TO 1230
         ISTAT =NSTATE(J3)
         ISTATM=ISTAT-1
         ISTATP=ISTAT+1
C
CL              5.1    IONISATION FOLLOWING BURGESS-CHIDICHIMO
C
C     USED FOR NICKEL
         IF ( NORDER(J3) .EQ. 16 ) THEN
C
C     LOOP OVER STATES OF IONISATION
           DO 515 J2=1,ISTAT
           ZQ=J2-1
           ZBETA=(SQRT((100.0*ZQ+91.0)/(4.0*ZQ+3.0))-5.0)/4.0
C
C     LOOP OVER TEMPERATURE
           DO 510 J1=1,NTE
           ZTE=TE(J1)
           ZSA=0.0
C
C     LOOP OVER SIX SHELLS
           DO 505 J=1,6
           ZPOT=(13.6/POT(J,J2,J3))**1.5
           Z=POT(J,J2,J3)/ZTE
           ZINV=1.0+1.0/Z
           ZE=EXP(-Z)*(ALOG(ZINV)-0.4/(1.0+Z)**2)
           ZW=ALOG(ZINV)**(ZBETA/ZINV)
           ZSA=ZSA+KSI(J,J2,J3)*ZPOT*SQRT(Z)*ZE*ZW
 505       CONTINUE
C
C     FILL ARRAY SAL WITH LOGARITHM OF IONISATION RATE
           ZSA=AMAX1(2.1715E-08*ZCON(J2,J3)*ZSA,1.0E-35)
           SAL(J1,J2,J3)=ALOG10(ZSA)
 510       CONTINUE
 515       CONTINUE
C
CL              5.2    IONISATION FOLLOWING BURGESS-CHIDICHIMO AND LOTZ
C
C     USED FOR CARBON, NITROGEN, OXYGEN, ARGON AND IRON
         ELSE IF ( ( NORDER(J3) .EQ. 5  )
     .       .OR.  ( NORDER(J3) .EQ. 6  )
     .       .OR.  ( NORDER(J3) .EQ. 7  )
     .       .OR.  ( NORDER(J3) .EQ. 12 )
     .       .OR.  ( NORDER(J3) .EQ. 15 ) ) THEN
 
C     LOOP OVER STATES OF IONISATION
           DO 530 J2=1,ISTAT
           ZQ=J2-1
           ZBETA=(SQRT((100.0*ZQ+91.0)/(4.0*ZQ+3.0))-5.0)/4.0
C
C     LOOP OVER TEMPERATURE
           DO 525 J1=1,NTE
           ZTE=TE(J1)
           ZSA=0.0
C
C     LOOP OVER SIX SHELLS
           DO 520 J=1,6
           ZPOT=(13.6/POT(J,J2,J3))**1.5
           Z=POT(J,J2,J3)/ZTE
           ZINV=1.0+1.0/Z
           ZE=EXP(-Z)*(ALOG(ZINV)-0.4/(1.0+Z)**2)
           ZW=ALOG(ZINV)**(ZBETA/ZINV)
           ZSA=ZSA+KSI(J,J2,J3)*ZPOT*SQRT(Z)*ZE*ZW
 520     CONTINUE
C
C     FILL ARRAY SAL WITH LOGARITHM OF IONISATION RATE
           SAL(J1,J2,J3)=2.1715E-08*ZCON(J2,J3)*ZSA
 525       CONTINUE
 530       CONTINUE
C
C     LOOP OVER STATES OF IONISATION
           DO 555 J2=1,ISTAT
C
C     LOOP OVER TEMPERATURE
           DO 550 J1=1,NTE
           ZTE=TE(J1)
           ZSA=0.0
C
C     LOOP OVER THREE SHELLS
           DO 545 J=1,3
           ZPP=ZP(J,J2,J3)/ZTE
           ZCC=ZC(J,J2,J3)
           ZPC=ZPP+ZCC
C
C     IONISATION GIVEN BY LOTZ, Z.PHYS. 216(1968) 241
           ZS=EXP(-ZPP)*(ALOG(1.0+1.0/ZPP)-0.4/(1.0+ZPP)**2)
           ZT=EXP(-ZPC)*(ALOG(1.0+1.0/ZPC)-0.4/(1.0+ZPC)**2)
           ZTERM=ZA(J,J2,J3)*IQ(J,J2,J3)*
     .            (ZS/ZPP-ZB(J,J2,J3)*EXP(ZCC)*ZT/ZPC)
           ZSA=ZSA+ZTERM
 545       CONTINUE
C
C     FILL ARRAY SAL WITH LOGARITHM OF IONISATION RATE
           ZSA=6.7E-07*ZSA/(ZTE*SQRT(ZTE))
           SAL(J1,J2,J3)=ALOG10(AMAX1(ZSA+SAL(J1,J2,J3),1.0E-35))
 550       CONTINUE
 555       CONTINUE
C
CL              5.3    IONISATION FOLLOWING LOTZ
C
C     USED FOR ALL OTHER ELEMENTS
         ELSE
C
C     LOOP OVER STATES OF IONISATION
           DO 580 J2=1,ISTAT
C
C     LOOP OVER TEMPERATURE
           DO 575 J1=1,NTE
           ZTE=TE(J1)
           ZSA=0.0
C
C     LOOP OVER THREE SHELLS
           DO 570 J=1,3
           ZPP=ZP(J,J2,J3)/ZTE
           ZCC=ZC(J,J2,J3)
           ZPC=ZPP+ZCC
C
C     IONISATION GIVEN BY LOTZ, Z.PHYS. 216(1968) 241
           ZS=EXP(-ZPP)*(ALOG(1.0+1.0/ZPP)-0.4/(1.0+ZPP)**2)
           ZT=EXP(-ZPC)*(ALOG(1.0+1.0/ZPC)-0.4/(1.0+ZPC)**2)
           ZTERM=ZA(J,J2,J3)*IQ(J,J2,J3)*
     .            (ZS/ZPP-ZB(J,J2,J3)*EXP(ZCC)*ZT/ZPC)
           ZSA=ZSA+ZTERM
 570       CONTINUE
C
C     FILL ARRAY SAL WITH LOGARITHM OF IONISATION RATE
           ZSA=AMAX1(6.7E-07*ZSA/(ZTE*SQRT(ZTE)),1.0E-35)
           SAL(J1,J2,J3)=ALOG10(ZSA)
 575       CONTINUE
 580       CONTINUE
         ENDIF
C
C-----------------------------------------------------------------------
CL              6.     RATE COEFFICIENTS FOR RADIATIVE RECOMBINATION
C
C     LOOP OVER STATES OF IONISATION
         DO 640 J2=1,ISTAT
         J2SQ=J2*J2
C
C     OPEN PLACES IN THE VALENCE SHELL - MU/NZ**3
         INZ=IZ(J2,J3)
         ZMU=(2.0*INZ**2-IK(J2,J3))/INZ**3
C
C     CONSTANT IN MAXIMUM SHELL CONTRIBUTING TO RADIATIVE RECOMBINATION
         ZCNS = 126.0*FLOAT(J2)**ZEXP3/13.6**ZEXP1
C
C     LOWER INDEX OF SUMMATION
         INZ1 = INZ+1.01
C
C     LOOP OVER TEMPERATURE
         DO 630 J1=1,NTE
         ZTE  = TE(J1)
         ZEHT = 13.6/ZTE
         ZPP1 = ZP(1,J2,J3)/ZTE
C
C     RADIATIVE CONTRIBUTION FROM THE VALENCE SHELL - VON GOELER
         ZS1=SQRT(ZEHT)*ZMU*J2SQ*ZPP1
     .      *(ALOG(1.0+1.0/ZPP1)-0.4/(1.0+ZPP1)**2)
C
C     HYDROGEN LIKE ENERGY LEVELS - Z**2*EH/TE
         ZEHJ2 = ZEHT*J2SQ
C
C     DENSITY INDEPENDENT PART OF CONTRIBUTIONS FROM HIGHER SHELLS
         ZCONS= 2.0*J2SQ**2*ZEHT**1.5
C
C     LOOP OVER DENSITIES
         DO 620 J0=1,NNE
         ZNE = DENE(J0)
C
C     HIGHEST CONTRIBUTING SHELL NEGLECTING EXPONENTIAL TERM IN NT
         INT = ZCNS*ZTE**ZEXP1/ZNE**ZEXP2+1.000001
         INT = MIN0(INT,1000)
C
C     RESET SUM TO FIRST CONTRIBUTION
         ZSUM = ZS1
C
C     CONTRIBUTION FROM HIGHER SHELLS
         DO 610 JJ=INZ1,INT
         ZPP2 = ZEHJ2/JJ**2
         ZSUM=ZSUM+ZCONS/JJ**3*(ALOG(1.0+1.0/ZPP2)-0.4/(1.0+ZPP2)**2)
 610     CONTINUE
C
C     FILL ARRAY RAL WITH LOGARITHM OF RADIATIVE RECOMBINATION RATE
         ZSUM=2.6E-14*ZSUM
         RAL(J0,J1,J2,J3) = ALOG10(ZSUM)
 620     CONTINUE
 630     CONTINUE
 640     CONTINUE
C
C-----------------------------------------------------------------------
CL              7.     DIELECTRONIC RECOMBINATION
C
CL              7.1        DIELECTRONIC RECOMBINATION FOLLOWING SOBELMAN
C
C     USED FOR NITROGEN
         IF ( NORDER(J3) .EQ. 6 ) THEN
C
C     LOOP OVER STATES OF IONISATION
         DO 715 J2=1,ISTATM
         J2P1  =J2+1
         J2P12 =J2P1*J2P1
         J2P13 =J2P12*J2P1
C
         ZDC = 1.0E-13*J2P13
C
C     LOOP OVER TEMPERATURE
         DO 710 J1=1,NTE
         ZEHT = 13.6/TE(J1)
         ZEHT15=ZEHT**1.5
C
C     LOOP OVER DENSITIES, BUT NO REAL DENSITY DEPENDENCE IN THIS CASE
         DO 705 J0=1,NNE
         ZDB=0.0
C
C     LOOP OVER TRANSITIONS
         DO 704 J=1,3
         ZS=EXP(-ZEHT*J2P12*DE(J,J2,J3))
         ZDIEL(J,J0,J1,J2,J3)=F(J,J2,J3)*G(J,J2,J3)*ZS*ZDC*ZEHT15
         ZDB=ZDB+ZDIEL(J,J0,J1,J2,J3)
 704     CONTINUE
C
C     FILL ARRAY RAL WITH LOGARITHM OF DIELECTRONIC RECOMBINATION RATE
         DAL(J0,J1,J2,J3)=ALOG10(AMAX1(ZDB,1.0E-35))
 705     CONTINUE
 710     CONTINUE
 715     CONTINUE
C
C     FULLY STRIPPED ION
         J2=ISTAT
C
C     LOOP OVER TEMPERATURE
         DO 725 J1=1,NTE
C
C     LOOP OVER DENSITIES
         DO 720 J0=1,NNE
C
C     FILL ARRAY DAL WITH "ZERO" DIELECTRONIC RECOMBINATION RATE
         ZDIEL(1,J0,J1,J2,J3)=1.0E-35
         ZDIEL(2,J0,J1,J2,J3)=1.0E-35
         ZDIEL(3,J0,J1,J2,J3)=1.0E-35
         DAL(J0,J1,J2,J3)=-75.0
 720     CONTINUE
 725     CONTINUE
C
CL              7.2        DIELECTRONIC RECOMBINATION FOLLOWING BURGESS
C
         ELSE
C
C     LOOP OVER STATES OF IONISATION
         DO 740 J2=1,ISTATM
         J2P1  =J2+1
         J2SQ  =J2*J2
         J2P1SQ=J2P1*J2P1
         ZEHJ2P=J2P1*13.6
C
C     TRANSITION ENERGY
         ZDE1= DE(1,J2P1,J3)
         ZDE2= DE(2,J2P1,J3)
         ZDE3= DE(3,J2P1,J3)
C
C     FACTOR Y
         ZY1 = ZDE1/ZEHJ2P
         ZY2 = ZDE2/ZEHJ2P
         ZY3 = ZDE3/ZEHJ2P
C
C     MULTIPLICATION FACTOR TO GET ORIGINAL BURGESS FORMULA
         ZGG1= 0.0
         ZGG2= 0.0
         ZGG3= 0.0
C
C     MULTIPLICATION FACTOR TO GET MODIFICATION BY MERTS
         IF ( G(1,J2P1,J3) .LT. 0.0 ) ZGG1= 1.0
         IF ( G(2,J2P1,J3) .LT. 0.0 ) ZGG2= 1.0
         IF ( G(3,J2P1,J3) .LT. 0.0 ) ZGG3= 1.0
C
C     FACTOR A IN ORIGINAL BURGESS
         ZZA11= 1.0+0.105*ZY1+0.015*ZY1**2
         ZZA12= 1.0+0.105*ZY2+0.015*ZY2**2
         ZZA13= 1.0+0.105*ZY3+0.015*ZY3**2
C
C     FACTOR A IN MERTS MODIFICATION
         ZZA21= 2.0+0.420*ZY1+0.060*ZY1**2
         ZZA22= 2.0+0.420*ZY2+0.060*ZY2**2
         ZZA23= 2.0+0.420*ZY3+0.060*ZY3**2
C
C     TOTAL FACTOR
         ZAA1= (1.0-ZGG1)*ZZA11+ZGG1*ZZA21
         ZAA2= (1.0-ZGG2)*ZZA12+ZGG2*ZZA22
         ZAA3= (1.0-ZGG3)*ZZA13+ZGG3*ZZA23
C
         ZD1 = F(1,J2P1,J3)*SQRT(ZY1)/ZAA1
         ZD2 = F(2,J2P1,J3)*SQRT(ZY2)/ZAA2
         ZD3 = F(3,J2P1,J3)*SQRT(ZY3)/ZAA3
C
         ZDC = 2.4E-09*J2P1SQ*SQRT(J2*J2P1/(J2**2+13.4))
C
C     AVERAGE TRANSITION ENERGY
         ZCOR= 1.0+0.015*J2*J2SQ/J2P1SQ
C
C     MODIFICATION IN EXPONENTIAL
         ZDE1= ZDE1/ZCOR
         ZDE2= ZDE2/ZCOR
         ZDE3= ZDE3/ZCOR
C
C     LOOP OVER TEMPERATURE
         DO 735 J1=1,NTE
         ZTE  = TE(J1)
         ZEHT = 13.6/ZTE
C
C     LOOP OVER DENSITIES
         DO 730 J0=1,NNE
         ZNE = DENE(J0)
         ZZN = 1.51E+17/ZNE * J2**6 * SQRT(ZTE)
         ZZN = ZZN**(1.0/7.0)
C
C     DENSITY DEPENDENCE
         ZN1 = (1.0-ZGG1)/(1.0+200.0/ZZN)+ZGG1/(1.0+666.7/(J2P1*ZZN)**2)
         ZN2 = (1.0-ZGG2)/(1.0+200.0/ZZN)+ZGG2/(1.0+666.7/(J2P1*ZZN)**2)
         ZN3 = (1.0-ZGG3)/(1.0+200.0/ZZN)+ZGG3/(1.0+666.7/(J2P1*ZZN)**2)
C
         ZDIEL(1,J0,J1,J2,J3)=ZN1/SQRT(ZTE**3)*ZD1*EXP(-ZDE1/ZTE)*ZDC
         ZDIEL(2,J0,J1,J2,J3)=ZN2/SQRT(ZTE**3)*ZD2*EXP(-ZDE2/ZTE)*ZDC
         ZDIEL(3,J0,J1,J2,J3)=ZN3/SQRT(ZTE**3)*ZD3*EXP(-ZDE3/ZTE)*ZDC
C
C     FILL ARRAY RAL WITH LOGARITHM OF DIELECTRONIC RECOMBINATION RATE
      ZDB=ZDIEL(1,J0,J1,J2,J3)+ZDIEL(2,J0,J1,J2,J3)+ZDIEL(3,J0,J1,J2,J3)
         DAL(J0,J1,J2,J3)=ALOG10(AMAX1(ZDB,1.0E-35))
 730     CONTINUE
 735     CONTINUE
 740     CONTINUE
C
C     FULLY STRIPPED ION
         J2=ISTAT
C
C     LOOP OVER TEMPERATURE
         DO 750 J1=1,NTE
C
C     LOOP OVER DENSITIES
         DO 745 J0=1,NNE
C
C     FILL ARRAY DAL WITH "ZERO" DIELECTRONIC RECOMBINATION RATE
         ZERO=1.0E-35/3.0
         ZDIEL(1,J0,J1,J2,J3)=ZERO
         ZDIEL(2,J0,J1,J2,J3)=ZERO
         ZDIEL(3,J0,J1,J2,J3)=ZERO
         DAL(J0,J1,J2,J3)=-75.0
 745     CONTINUE
 750     CONTINUE
         ENDIF
C
C-----------------------------------------------------------------------
CL              8.     LINE RADIATION LOSSES
C
C     USED FOR NITROGEN
         IF ( NORDER(J3) .EQ. 6 ) THEN
C
CL              8.1        RADIATION FOLLOWING VAINSHTEIN
C
C     NEUTRAL STATE
         J2=1
C
C     LOOP OVER TEMPERATURE
         DO 820 J1=1,NTE
         ZTE=TE(J1)
         ZLINE=0.0
         ZPE =ZP(1,J2,J3)
         ZPEH=(13.6/ZPE)**1.5
C
C     LOOP OVER TRANSITIONS
         DO 810 J=1,10
         IF ( ZCHI(J,J2,J3) .GT. 100.0 ) THEN
           ZDEL=ZDE(J,J2,J3)
           ZD=((ZPE-ZDEL)/ZDEL)**1.5
           ZBETA=ZDEL/ZTE
           ZFAC=SQRT(ZBETA+3.0)*ALOG(16.0+1.0/ZBETA)
           ZLINE=ZLINE+ZD*ZDEL*ZBE(J,J2,J3)*ZFAC*EXP(-ZBETA)
     .                *SQRT(ZBETA)/(ZBETA+ZCHI(J,J2,J3)-100.0)
         ELSE IF ( ZCHI(J,J2,J3) .GT. 0.0 ) THEN
           ZDEL=ZDE(J,J2,J3)
           ZD=((ZPE-ZDEL)/ZDEL)**1.5
           ZBETA=ZDEL/ZTE
           ZFAC=SQRT(ZBETA+1.0)
           ZLINE=ZLINE+ZD*ZDEL*ZBE(J,J2,J3)*ZFAC*EXP(-ZBETA)
     .                *SQRT(ZBETA)/(ZBETA+ZCHI(J,J2,J3))
         ELSE IF ( ZCHI(J,J2,J3) .LT. 0.0 ) THEN
           ZDEL=ZDE(J,J2,J3)
           ZD=((ZPE-ZDEL)/ZDEL)**1.5
           ZBETA=ZDEL/ZTE
           ZFAC=ZBETA/SQRT(ZBETA+1.0)
           ZLINE=ZLINE+ZD*ZDEL*ZBE(J,J2,J3)*ZFAC*EXP(-ZBETA)
     .                *SQRT(ZBETA)/(ZBETA-ZCHI(J,J2,J3))
         ENDIF
 810     CONTINUE
C
C     FILL ARRAY RLINE WITH LOGARITHM LINE RADIATION
         RLINE(J1,J2,J3)=ALOG10(AMAX1(ZLINE*ZPEH*1.6E-27,1.0E-35))
 820     CONTINUE
C
C     LOOP OVER STATES OF IONISATION
         DO 850 J2=2,ISTAT
C
C     LOOP OVER TEMPERATURE
         DO 840 J1=1,NTE
         ZTE=TE(J1)
         ZLINE=0.0
         ZPE =ZP(1,J2,J3)
         ZPEH=(13.6/ZPE)**1.5
C
C     LOOP OVER TRANSITIONS
         DO 830 J=1,10
         IF ( ZCHI(J,J2,J3) .GT. 100.0 ) THEN
           ZDEL=ZDE(J,J2,J3)
           ZD=((ZPE-ZDEL)/ZDEL)**1.5
           ZBETA=ZDEL/ZTE
           ZFAC=SQRT(ZBETA+3.0)*ALOG(16.0+1.0/ZBETA)
           ZLINE=ZLINE+ZD*ZDEL*ZBE(J,J2,J3)*ZFAC*EXP(-ZBETA)
     .                *SQRT(ZBETA)/(ZBETA+ZCHI(J,J2,J3)-100.0)
         ELSE IF ( ZCHI(J,J2,J3) .GT. 0.0 ) THEN
           ZDEL=ZDE(J,J2,J3)
           ZD=((ZPE-ZDEL)/ZDEL)**1.5
           ZBETA=ZDEL/ZTE
           ZFAC=ZBETA+1.0
           ZLINE=ZLINE+ZD*ZDEL*ZBE(J,J2,J3)*ZFAC*EXP(-ZBETA)
     .                *SQRT(ZBETA)/(ZBETA+ZCHI(J,J2,J3))
         ELSE IF ( ZCHI(J,J2,J3) .LT. 0.0 ) THEN
           ZDEL=ZDE(J,J2,J3)
           ZD=((ZPE-ZDEL)/ZDEL)**1.5
           ZBETA=ZDEL/ZTE
           ZFAC=ZBETA
           ZLINE=ZLINE+ZD*ZDEL*ZBE(J,J2,J3)*ZFAC*EXP(-ZBETA)
     .                *SQRT(ZBETA)/(ZBETA-ZCHI(J,J2,J3))
         ENDIF
 830     CONTINUE
C
C     FILL ARRAY RLINE WITH LOGARITHM LINE RADIATION
         RLINE(J1,J2,J3)=ALOG10(AMAX1(ZLINE*ZPEH*1.6E-27,1.0E-35))
 840     CONTINUE
 850     CONTINUE
C
C     FULLY STRIPPED ION
         J2=ISTATP
C
C     LOOP OVER TEMPERATURE
         DO 860 J1=1,NTE
C
C     FILL ARRAY RLINE WITH "ZERO" LINE RADIATION
         RLINE(J1,J2,J3)=-75.0
 860     CONTINUE
C
         ELSE
CL              8.2        LINE RADIATION FOLLOWING VAN REGEMORTER
C
C     LOOP OVER STATES OF IONISATION
         DO 880 J2=1,ISTAT
C
C     LOOP OVER TEMPERATURE
         DO 870 J1=1,NTE
         ZTE=TE(J1)
         ZSQRTE=SQRT(ZTE)
         ZLINE=0.0
C
C     LOOP OVER TRANSITIONS
         DO 865 J=1,3
         ZT=F(J,J2,J3)*ABS(G(J,J2,J3))*EXP(-DE(J,J2,J3)/ZTE)
         ZLINE=ZLINE+ZT
 865     CONTINUE
C
C     FILL ARRAY RLINE WITH LOGARITHM OF LINE RADIATION
         RLINE(J1,J2,J3)=ALOG10(AMAX1(ZLINE*2.53E-24/ZSQRTE,1.0E-35))
 870     CONTINUE
 880     CONTINUE
C
C     FULLY STRIPPED ION
         J2=ISTATP
C
C     LOOP OVER TEMPERATURE
         DO 890 J1=1,NTE
C
C     FILL ARRAY RLINE WITH "ZERO" LINE RADIATION
         RLINE(J1,J2,J3)=-75.0
 890     CONTINUE
         ENDIF
C
C-----------------------------------------------------------------------
CL              9.         RADIATIVE RECOMBINATION LOSSES
C
C     LOOP OVER STATES OF IONISATION
         DO 930 J2=1,ISTAT
C
C     LOOP OVER TEMPERATURE
         DO 920 J1=1,NTE
         ZPELOG=ALOG10(ZP(1,J2,J3))+ALOG10(1.6)-19.0
C
C     LOOP OVER DENSITIES
         DO 910 J0=1,NNE
         RAD(J0,J1,J2,J3)=RAL(J0,J1,J2,J3)+ZPELOG
 910     CONTINUE
 920     CONTINUE
 930     CONTINUE
C
C-----------------------------------------------------------------------
CL             10.         DIELECTRONIC RECOMBINATION LOSSES
C
C     LOOP OVER STATES OF IONISATION
         DO 1040 J2=1,ISTATM
         J2P1=J2+1
         ZPE=ZP(1,J2,J3)
C
C     LOOP OVER TEMPERATURE
         DO 1030 J1=1,NTE
C
C     LOOP OVER DENSITIES
         DO 1020 J0=1,NNE
         ZDIE=0.0
C
C     LOOP OVER TRANSITIONS
         DO 1010 J=1,3
         IF ( NORDER(J3) .EQ. 6 ) THEN
         ZDIE=ZDIE+ZDIEL(J,J0,J1,J2,J3)*(ZPE+DE(J,J2,J3)*13.6*J2P1*J2P1)
         ELSE
         ZDIE=ZDIE+ZDIEL(J,J0,J1,J2,J3)*(ZPE+DE(J,J2P1,J3))
         ENDIF
 1010    CONTINUE
         DIE(J0,J1,J2,J3)=ALOG10(AMAX1(ZDIE*1.6E-19,1.0E-35))
 1020    CONTINUE
 1030    CONTINUE
 1040    CONTINUE
C
C     FULLY STRIPPED ION
         J2=ISTAT
C
C     LOOP OVER TEMPERATURE
         DO 1060 J1=1,NTE
C
C     LOOP OVER DENSITIES
         DO 1050 J0=1,NNE
         DIE(J0,J1,J2,J3)=-75.0
 1050    CONTINUE
 1060    CONTINUE
C
C-----------------------------------------------------------------------
CL             11.         BREMSSTRAHLUNG LOSSES
C
C     LOOP OVER STATES OF IONISATION
         DO 1120 J2=1,ISTAT
         Z2=1.542E-32*J2**2
C
C     LOOP OVER TEMPERATURE
         DO 1110 J1=1,NTE
         ZSQRTE=SQRT(TE(J1))
         BREMS(J1,J2,J3)=ALOG10(AMAX1(Z2*ZSQRTE,1.0E-35))
 1110    CONTINUE
 1120    CONTINUE
C
C-----------------------------------------------------------------------
CL             12.         IONISATION ENERGY LOSSES
C
C
C     LOOP OVER STATES OF IONISATION
         DO 1220 J2=1,ISTAT
C
C     LOOP OVER TEMPERATURE
         DO 1210 J1=1,NTE
         ZPELOG=ALOG10(ZP(1,J2,J3))+ALOG10(1.6)-19.0
         RION(J1,J2,J3)=SAL(J1,J2,J3)+ZPELOG
 1210    CONTINUE
 1220    CONTINUE
 1230    CONTINUE
C
C-----------------------------------------------------------------------
CL             13.         WRITE ATOMIC DATA TO FILE KOUT
C
C     CHECK IF OUTPUT IS WANTED
         IF ( KOUT .GT. 0 ) THEN
C
C     SCAN OVER SPECIES
         DO 1390 J3=1,NIMP
         IF ( NORDER(J3) .EQ. 0 ) GO TO 1390
         IDIAL=NSTATE(J3)
C
C     USED FOR NICKEL
         IF ( IDIAL .EQ. 28 ) THEN
C
C     WRITE TITLE FOR BURGESS-CHIDICHIMO COEFFICIENTS
           WRITE (KOUT,1)
C
C     SCAN OVER STATES OF IONISATION
           DO 1315 J2=1,IDIAL
C
C     DETERMINE ON WHICH LINE ION IDENTIFICATION WILL BE WRITTEN
           IMID=0
           DO 1305 J=1,6
           IF ( KSI(J,J2,J3) .GT. 0 ) THEN
             IMID=IMID+1
           ENDIF
 1305      CONTINUE
           IMID2=IMID/2
           IF ( IMID2*2 .EQ. IMID ) THEN
             IMID=IMID2
           ELSE
             IMID=IMID2+1
           ENDIF
C
C     WRITE BURGESS-CHIDICHIMO COEFFICIENTS
           IF ( IMID .EQ. 0 ) GO TO 1315
           WRITE (KOUT,20)
           DO 1310 J=1,6
           IF ( KSI(J,J2,J3) .GT. 0 ) THEN
             IF ( J .NE. IMID ) THEN
               WRITE (KOUT,2) POT(J,J2,J3),KSI(J,J2,J3)
             ELSE
               WRITE (KOUT,3) ( INAME(JJ,J2,J3), JJ=1,3 ),ZP(1,J2,J3),
     .                        POT(J,J2,J3),KSI(J,J2,J3),ZCON(J2,J3)
             ENDIF
           ENDIF
 1310      CONTINUE
           WRITE (KOUT,4)
 1315      CONTINUE
C
C     USED FOR CARBON, NITROGEN, OXYGEN, ARGON, IRON
         ELSE IF ( ( IDIAL .EQ. 6 )
     .      .OR.   ( IDIAL .EQ. 7 )
     .      .OR.   ( IDIAL .EQ. 8 )
     .      .OR.   ( IDIAL .EQ. 18)
     .      .OR.   ( IDIAL .EQ. 26) ) THEN
C
C     WRITE TITLE FOR BURGESS-CHIDICHIMO COEFFICIENTS
           WRITE (KOUT,1)
C
C     SCAN OVER STATES OF IONISATION
           DO 1330 J2=1,IDIAL
C
C     DETERMINE ON WHICH LINE ION IDENTIFICATION WILL BE WRITTEN
           IMID=0
           DO 1320 J=1,6
           IF ( KSI(J,J2,J3) .GT. 0 ) THEN
             IMID=IMID+1
           ENDIF
 1320      CONTINUE
           IMID2=IMID/2
           IF ( IMID2*2 .EQ. IMID ) THEN
             IMID=IMID2
           ELSE
             IMID=IMID2+1
           ENDIF
C
C     WRITE BURGESS-CHIDICHIMO COEFFICIENTS
           IF ( IMID .EQ. 0 ) GO TO 1330
           WRITE (KOUT,20)
           DO 1325 J=1,6
           IF ( KSI(J,J2,J3) .GT. 0 ) THEN
             IF ( J .NE. IMID ) THEN
               WRITE (KOUT,2) POT(J,J2,J3),KSI(J,J2,J3)
             ELSE
               WRITE (KOUT,3) ( INAME(JJ,J2,J3), JJ=1,3 ),ZP(1,J2,J3),
     .                        POT(J,J2,J3),KSI(J,J2,J3),ZCON(J2,J3)
             ENDIF
           ENDIF
 1325      CONTINUE
           WRITE (KOUT,4)
 1330      CONTINUE
C
C     WRITE TITLE FOR LOTZ COEFFICIENTS
           WRITE (KOUT,5)
C
C     SCAN OVER STATES OF IONISATION
           DO 1345 J2=1,IDIAL
C
C     DETERMINE ON WHICH LINE ION IDENTIFICATION WILL BE WRITTEN
           IMID=0
           DO 1335 J=1,3
           IF ( IQ(J,J2,J3) .GT. 0 ) THEN
             IMID=IMID+1
           ENDIF
 1335      CONTINUE
           IMID2=IMID/2
           IF ( IMID2*2 .EQ. IMID ) THEN
             IMID=IMID2
           ELSE
             IMID=IMID2+1
           ENDIF
C
C     WRITE LOTZ COEFFICIENTS
           IF ( IMID .EQ. 0 ) GO TO 1345
           WRITE (KOUT,18)
           DO 1340 J=1,3
           IF ( IQ(J,J2,J3) .GT. 0 ) THEN
             IF ( J .NE. IMID ) THEN
             WRITE (KOUT,6) ZP(J,J2,J3),ZA(J,J2,J3),ZB(J,J2,J3),
     .                      ZC(J,J2,J3),IQ(J,J2,J3)
             ELSE
             WRITE (KOUT,7) ( INAME(JJ,J2,J3), JJ=1,3 ),
     .                      ZP(J,J2,J3),ZA(J,J2,J3),ZB(J,J2,J3),
     .                      ZC(J,J2,J3),IQ(J,J2,J3)
             ENDIF
           ENDIF
 1340      CONTINUE
           WRITE (KOUT,8)
 1345      CONTINUE
C
C     USED FOR ALL OTHER ELEMENTS
         ELSE
C
C     WRITE TITLE FOR LOTZ COEFFICIENTS
           WRITE (KOUT,5)
C
C     SCAN OVER STATES OF IONISATION
           DO 1360 J2=1,IDIAL
C
C     DETERMINE ON WHICH LINE ION IDENTIFICATION WILL BE WRITTEN
           IMID=0
           DO 1350 J=1,3
           IF ( IQ(J,J2,J3) .GT. 0 ) THEN
             IMID=IMID+1
           ENDIF
 1350      CONTINUE
           IMID2=IMID/2
           IF ( IMID2*2 .EQ. IMID ) THEN
             IMID=IMID2
           ELSE
             IMID=IMID2+1
           ENDIF
C
C     WRITE LOTZ COEFFICIENTS
           IF ( IMID .EQ. 0 ) GO TO 1360
           WRITE (KOUT,18)
           DO 1355 J=1,3
           IF ( IQ(J,J2,J3) .GT. 0 ) THEN
             IF ( J .NE. IMID ) THEN
               WRITE (KOUT,6) ZP(J,J2,J3),ZA(J,J2,J3),ZB(J,J2,J3),
     .                      ZC(J,J2,J3),IQ(J,J2,J3)
             ELSE
               WRITE (KOUT,7) ( INAME(JJ,J2,J3), JJ=1,3 ),
     .                        ZP(J,J2,J3),ZA(J,J2,J3),ZB(J,J2,J3),
     .                        ZC(J,J2,J3),IQ(J,J2,J3)
             ENDIF
           ENDIF
 1355      CONTINUE
           WRITE (KOUT,8)
 1360      CONTINUE
         ENDIF
C
C     WRITE TITLE FOR SOBELMAN TRANSITION DATA
         IF ( IDIAL .EQ.  7 ) THEN
           WRITE (KOUT,9)
C
C     SCAN OVER STATES OF IONISATION
           DO 1375 J2=1,IDIAL
C
C     DETERMINE ON WHICH LINE ION IDENTIFICATION WILL BE WRITTEN
           IMID=0
           DO 1365 J=1,3
           IF ( F(J,J2,J3) .GT. 0.0 ) THEN
             IMID=IMID+1
           ENDIF
 1365      CONTINUE
           IMID2=IMID/2
           IF ( IMID2*2 .EQ. IMID ) THEN
             IMID=IMID2
           ELSE
             IMID=IMID2+1
           ENDIF
C
C     WRITE COEFFICIENTS
           IF ( IMID .EQ. 0 ) GO TO 1375
           WRITE (KOUT,19)
           DO 1370 J=1,3
           IF ( F(J,J2,J3) .GT. 0.0 ) THEN
             IF ( J .NE. IMID ) THEN
               WRITE (KOUT,11) DE(J,J2,J3),F(J,J2,J3),G(J,J2,J3)
             ELSE
               WRITE (KOUT,12) ( INAME(JJ,J2,J3), JJ=1,3 ),
     .                         IZ(J2,J3),IK(J2,J3),DE(J,J2,J3),
     .                         F(J,J2,J3),G(J,J2,J3)
             ENDIF
           ENDIF
 1370      CONTINUE
           WRITE (KOUT,13)
 1375      CONTINUE
C
         ELSE
C
           WRITE (KOUT,10)
C
C     SCAN OVER STATES OF IONISATION
           DO 1354 J2=1,IDIAL
           WRITE (KOUT,22)
C
C     DETERMINE ON WHICH LINE ION IDENTIFICATION WILL BE WRITTEN
           IMID=0
           DO 1351 J=1,3
           IF ( F(J,J2,J3) .GT. 0.0 ) THEN
             IMID=IMID+1
           ENDIF
 1351      CONTINUE
           IMID2=IMID/2
           IF ( IMID2*2 .EQ. IMID ) THEN
             IMID=IMID2
           ELSE
             IMID=IMID2+1
           ENDIF
C
C     WRITE TRANSITION DATA
           DO 1352 J=1,3
           IF ( F(J,J2,J3) .GT. 0.0 ) THEN
             IDN=1
             ZGAUNT=ABS(G(J,J2,J3))
             IF ( ZGAUNT .EQ. G(J,J2,J3) ) IDN=0
             IF ( J .NE. IMID ) THEN
               IF ( IDN .EQ. 0 ) THEN
                 WRITE (KOUT,23) DE(J,J2,J3),F(J,J2,J3),ZGAUNT,IDN
               ELSE
                 IDN=0
                 WRITE (KOUT,26) DE(J,J2,J3),F(J,J2,J3),ZGAUNT
               ENDIF
             ELSE
               IF ( IDN .EQ. 0 ) THEN
                 WRITE (KOUT,24) ( INAME(JJ,J2,J3), JJ=1,3 ),
     .                           IZ(J2,J3),IK(J2,J3),DE(J,J2,J3),
     .                           F(J,J2,J3),ZGAUNT,IDN
               ELSE
                 IDN=0
                 WRITE (KOUT,27) ( INAME(JJ,J2,J3), JJ=1,3 ),
     .                           IZ(J2,J3),IK(J2,J3),DE(J,J2,J3),
     .                           F(J,J2,J3),ZGAUNT
               ENDIF
             ENDIF
           ELSE IF ( J .EQ. 1 ) THEN
             WRITE (KOUT,28) ( INAME(JJ,J2,J3), JJ=1,3 ),
     .                       IZ(J2,J3),IK(J2,J3)
           ENDIF
 1352      CONTINUE
           WRITE (KOUT,25)
 1354      CONTINUE
         ENDIF
C
C     WRITE TITLE FOR VAINSHTEIN LINE RADIATION DATA
         IF ( ( IDIAL .EQ.  7 ) ) THEN
           WRITE (KOUT,14)
C
C     SCAN OVER STATES OF IONISATION
           ISTAR=0
           DO 1380 J2=1,IDIAL
           WRITE (KOUT,21)
C
C     DETERMINE ON WHICH LINE ION IDENTIFICATION WILL BE WRITTEN
           IMID=0
           DO 1361 J=1,10
           IF ( ZCHI(J,J2,J3) .GT. 0.0 ) THEN
             IMID=IMID+1
           ELSE IF ( ZCHI(J,J2,J3) .LT. 0.0 ) THEN
             IMID=IMID+1
           ENDIF
 1361      CONTINUE
           IMID2=IMID/2
           IF ( IMID2*2 .EQ. IMID ) THEN
             IMID=IMID2
           ELSE
             IMID=IMID2+1
           ENDIF
C
C     WRITE COEFFICIENTS
           DO 1371 J=1,10
           IF ( ZCHI(J,J2,J3) .GT. 100.0 ) THEN
             ZCHI(J,J2,J3)=ZCHI(J,J2,J3)-100.0
             ISTAR=1
             IDS=0
             IF ( J .NE. IMID ) THEN
               WRITE (KOUT,29) ZDE(J,J2,J3),ZBE(J,J2,J3),
     .                         ZCHI(J,J2,J3),IDS
             ELSE
               WRITE (KOUT,30) ( INAME(JJ,J2,J3), JJ=1,3 ),
     .                         ZDE(J,J2,J3),ZBE(J,J2,J3),
     .                         ZCHI(J,J2,J3),IDS
             ENDIF
           ELSE IF ( ZCHI(J,J2,J3) .GT. 0.0 ) THEN
             IDS=0
             IF ( J .NE. IMID ) THEN
               WRITE (KOUT,15) ZDE(J,J2,J3),ZBE(J,J2,J3),
     .                         ZCHI(J,J2,J3),IDS
             ELSE
               WRITE (KOUT,16) ( INAME(JJ,J2,J3), JJ=1,3 ),
     .                         ZDE(J,J2,J3),ZBE(J,J2,J3),
     .                         ZCHI(J,J2,J3),IDS
             ENDIF
           ELSE IF ( ZCHI(J,J2,J3) .LT. 0.0 ) THEN
             IDS=1
             IF ( J .NE. IMID ) THEN
               WRITE (KOUT,15) ZDE(J,J2,J3),ZBE(J,J2,J3),
     .                         ABS(ZCHI(J,J2,J3)),IDS
             ELSE
               WRITE (KOUT,16) ( INAME(JJ,J2,J3), JJ=1,3 ),
     .                         ZDE(J,J2,J3),ZBE(J,J2,J3),
     .                         ABS(ZCHI(J,J2,J3)),IDS
             ENDIF
           ENDIF
 1371      CONTINUE
           WRITE (KOUT,17)
 1380      CONTINUE
         ENDIF
         IF ( ISTAR .EQ. 1 ) WRITE(KOUT,31)
 1390    CONTINUE
         ENDIF
C
C     LIST OF FORMATS
 1       FORMAT(//9X,1HI,13(1H-),1HI,39(1H-),1HI/9X,1HI,13X,1HI,39X,1HI/
     .    9X,1HI,13X,1HI,6X,27HCOEFFICIENTS FOR IONISATION,6X,1HI/
     .    9X,1HI,13X,1HI,39X,1HI/
     .    9X,1HI,13(1H-),2(1HI,12(1H-)),2(1HI,6(1H-)),1HI/
     .    9X,1HI,13X,1HI,2(12X,1HI),2(6X,1HI)/
     .    9X,1HI,5X,3HION,5X,1HI,3X,6HP (EV),3X,1HI,3X,6HE (EV),3X,1HI,
     .    1X,4HZETA,1X,1HI,2X,1HC,3X,1HI/
     .    9X,1HI,13X,1HI,2(12X,1HI),2(6X,1HI)/
     .    9X,1HI,13(1H-),2(1HI,12(1H-)),2(1HI,6(1H-)),1HI)
C
 20      FORMAT(9X,1HI,13X,1HI,2(12X,1HI),2(6X,1HI))
C
 2       FORMAT(9X,1HI,13X,1HI,12X,1HI,F10.3,2X,1HI,I4,2X,1HI,6X,1HI)
C
 3       FORMAT(9X,1HI,1X,3A4,1HI,F10.3,2X,1HI,F10.3,2X,1HI,I4,2X,1HI,
     .    F5.2,1X,1HI)
C
 4       FORMAT(9X,1HI,13X,1HI,2(12X,1HI),2(6X,1HI)/
     .    9X,1HI,13(1H-),2(1HI,12(1H-)),2(1HI,6(1H-)),1HI)
C
 5       FORMAT(//4X,1HI,13(1H-),1HI,49(1H-),1HI/4X,1HI,13X,1HI,49X,1HI/
     .    4X,1HI,13X,1HI,9X,32HLOTZ COEFFICIENTS FOR IONISATION,8X,1HI/
     .    4X,1HI,13X,1HI,49X,1HI/
     .    4X,1HI,13(1H-),2(1HI,12(1H-)),2(1HI,8(1H-)),1HI,5(1H-),1HI/
     .    4X,1HI,13X,2(1HI,12X),2(1HI,8X),1HI,5X,1HI/
     .    4X,1HI,5X,3HION,5X,1HI,3X,6HP (EV),3X,1HI,1X,10HA (CM2EV2),
     .    1X,1HI,4X,1HB,3X,1HI,4X,1HC,3X,1HI,2X,1HQ,2X,1HI/
     .    4X,1HI,13X,2(1HI,12X),2(1HI,8X),1HI,5X,1HI/
     .    4X,1HI,13(1H-),2(1HI,12(1H-)),2(1HI,8(1H-)),1HI,5(1H-),1HI)
C
 18      FORMAT(4X,1HI,13X,2(1HI,12X),2(1HI,8X),1HI,5X,1HI)
C
 6       FORMAT(4X,1HI,13X,1HI,F10.3,2X,1HI,F8.3,4X,1HI,F6.2,2X,1HI,
     .    F6.2,2X,1HI,I3,2X,1HI)
C
 7       FORMAT(4X,1HI,1X,3A4,1HI,F10.3,2X,1HI,F8.3,4X,1HI,F6.2,2X,1HI,
     .    F6.2,2X,1HI,I3,2X,1HI)
C
 8       FORMAT(4X,1HI,13X,2(1HI,12X),2(1HI,8X),1HI,5X,1HI/
     .    4X,1HI,13(1H-),2(1HI,12(1H-)),2(1HI,8(1H-)),1HI,5(1H-),1HI)
C
 9       FORMAT(//6X,1HI,13(1H-),1HI,12(1H-),1HI,32(1H-),1HI/
     .    6X,1HI,13X,1HI,5X,3HION,4X,1HI,32X,1HI/6X,1HI,13X,1HI,
     .    2X,9HSTRUCTURE,1X,1HI,9X,15HTRANSITION DATA,8X,1HI/
     .    6X,1HI,13X,1HI,12X,1HI,32X,1HI/
     .    6X,1HI,13(1H-),1HI,5(1H-),1HI,6(1H-),1HI,3(10(1H-),1HI)/
     .    6X,1HI,13X,1HI,5X,1HI,6X,1HI,3(10X,1HI)/
     .    6X,1HI,5X,3HION,5X,1HI,2X,2HNZ,1X,1HI,2X,3HKSI,1X,1HI,
     .    4X,3HCHI,3X,1HI,5X,1HQ,4X,1HI,5X,1HA,4X,1HI/
     .    6X,1HI,13X,1HI,5X,1HI,6X,1HI,3(10X,1HI)/
     .    6X,1HI,13(1H-),1HI,5(1H-),1HI,6(1H-),1HI,3(10(1H-),1HI))
C
 19      FORMAT(6X,1HI,13X,1HI,5X,1HI,6X,1HI,3(10X,1HI))
C
 10      FORMAT(//1X,1HI,13(1H-),1HI,12(1H-),1HI,41(1H-),1HI/
     .    1X,1HI,13X,1HI,5X,3HION,4X,1HI,41X,1HI/1X,1HI,13X,1HI,
     .    2X,9HSTRUCTURE,1X,1HI,13X,15HTRANSITION DATA,13X,1HI/
     .    1X,1HI,13X,1HI,12X,1HI,41X,1HI/
     .    1X,1HI,13(1H-),1HI,5(1H-),1HI,6(1H-),1HI,3(10(1H-),1HI),
     .    8(1H-),1HI/
     .    1X,1HI,13X,1HI,5X,1HI,6X,1HI,3(10X,1HI),8X,1HI/
     .    1X,1HI,5X,3HION,5X,1HI,2X,2HNZ,1X,1HI,2X,3HKSI,1X,1HI,
     .    4X,2HDE,4X,1HI,5X,1HF,4X,1HI,5X,1HG,4X,1HI,3X,2HDN,3X,1HI/
     .    1X,1HI,13X,1HI,5X,1HI,6X,1HI,3(10X,1HI),8X,1HI/
     .    1X,1HI,13(1H-),1HI,5(1H-),1HI,6(1H-),1HI,3(10(1H-),1HI),
     .    8(1H-),1HI)
C
 11      FORMAT(6X,1HI,13X,1HI,5X,1HI,6X,1HI,F8.3,2X,1HI,2(F8.2,2X,1HI))
C
 12      FORMAT(6X,1HI,1X,3A4,1HI,I3,2X,1HI,I4,2X,1HI,F8.3,2X,1HI,
     .    2(F8.2,2X,1HI))
C
 13      FORMAT(6X,1HI,13X,1HI,5X,1HI,6X,1HI,3(10X,1HI)/
     .    6X,1HI,13(1H-),1HI,5(1H-),1HI,6(1H-),1HI,3(10(1H-),1HI))
C
 14      FORMAT(//4X,1HI,13(1H-),1HI,49(1H-),1HI/4X,1HI,13X,1HI,49X,1HI/
     .    4X,1HI,13X,1HI,4X,42HVAINSHTEIN COEFFICIENTS FOR LINE RADIATIO
     .N,  3X,1HI/4X,1HI,13X,1HI,49X,1HI/
     .    4X,1HI,13(1H-),2(1HI,12(1H-)),2(1HI,11(1H-)),1HI/
     .    4X,1HI,13X,2(1HI,12X),2(1HI,11X),1HI/
     .    4X,1HI,5X,3HION,5X,1HI,3X,7HDE (EV),2X,1HI,6X,1HB,
     .    5X,1HI,4X,3HCHI,4X,1HI,5X,2HDS,4X,1HI/
     .    4X,1HI,13X,2(1HI,12X),2(1HI,11X),1HI/
     .    4X,1HI,13(1H-),2(1HI,12(1H-)),2(1HI,11(1H-)),1HI)
C
 15      FORMAT(4X,1HI,13X,1HI,1X,F8.2,3X,1HI,F9.3,3X,1HI,F9.3,2X,
     .    1HI,4X,I2,5X,1HI)
C
 16      FORMAT(4X,1HI,1X,3A4,1HI,1X,F8.2,3X,1HI,F9.3,3X,1HI,F9.3,2X,
     .    1HI,4X,I2,5X,1HI)
C
 29      FORMAT(4X,1HI,13X,1HI,1X,F8.2,3X,1HI,F9.3,3X,1HI,1X,1H*,1X,
     .    F6.3,2X,1HI,4X,I2,5X,1HI)
C
 30      FORMAT(4X,1HI,1X,3A4,1HI,1X,F8.2,3X,1HI,F9.3,3X,1HI,1X,1H*,1X,
     .    F6.3,2X,1HI,4X,I2,5X,1HI)
C
 31      FORMAT(//9X,48H*: SEE 'A PACKAGE FOR NON-CORONAL IMPURITY DATA'
     .    /12X,49HPART I, SECTION 2.4, BY A.E.P.M. ABELS-VAN MAANEN//)
C
 17      FORMAT(4X,1HI,13X,2(1HI,12X),2(1HI,11X),1HI/
     .    4X,1HI,13(1H-),2(1HI,12(1H-)),2(1HI,11(1H-)),1HI)
C
 21      FORMAT(4X,1HI,13X,2(1HI,12X),2(1HI,11X),1HI)
C
 22      FORMAT(1X,1HI,13X,1HI,5X,1HI,6X,1HI,3(10X,1HI),8X,1HI)
C
 23      FORMAT(1X,1HI,13X,1HI,5X,1HI,6X,1HI,1X,F8.3,1X,1HI,F8.3,2X,1HI,
     .    F8.2,2X,1HI,4X,I1,3X,1HI)
C
 24      FORMAT(1X,1HI,1X,3A4,1HI,I3,2X,1HI,I4,2X,1HI,1X,F8.3,1X,1HI,
     .    F8.3,2X,1HI,F8.2,2X,1HI,4X,I1,3X,1HI)
C
 25      FORMAT(1X,1HI,13X,1HI,5X,1HI,6X,1HI,3(10X,1HI),8X,1HI/
     .    1X,1HI,13(1H-),1HI,5(1H-),1HI,6(1H-),1HI,3(10(1H-),1HI),
     .    8(1H-),1HI)
C
 26      FORMAT(1X,1HI,13X,1HI,5X,1HI,6X,1HI,1X,F8.3,1X,1HI,F8.3,2X,1HI,
     .    F8.2,2X,1HI,4X,1H1,3X,1HI)
C
 27      FORMAT(1X,1HI,1X,3A4,1HI,I3,2X,1HI,I4,2X,1HI,1X,F8.3,1X,1HI,
     .    F8.3,2X,1HI,F8.2,2X,1HI,4X,1H1,3X,1HI)
C
 28      FORMAT(1X,1HI,1X,3A4,1HI,I3,2X,1HI,I4,2X,1HI,3(10X,1HI),8X,1HI)
C
         RETURN
         END
C/ MODULE A2
C
         SUBROUTINE RRATES(PTE,PNE,POUT,KMESH)
         implicit none
C
      integer kmesh
      REAL  PTE(KMESH), PNE(KMESH), POUT(4,28,3,KMESH)
C
C A.2  CALCULATE RATE COEFFICIENTS
C
C VERSION 02/12/85 AVMA JET/OXFORD
C
C-----------------------------------------------------------------------
CL                  C1.1     IMPURITY IDENTIFICATION
C VERSION 07/09/84 AVMA JET/OXFORD
       COMMON/COMAT1/
     I   NORDER,   NIMP  ,   NSTATE
       integer
     I   NORDER(3)       ,   NSTATE(3), nimp
C-----------------------------------------------------------------------
CL                  C1.2     INTERPOLATION VARIABLES
C VERSION 07/09/84 AVMA JET/OXFORD
       COMMON/COMAT2/
     R   TE    ,   TELOG ,   DTLOG ,   TMINLG,   TMAXLG,   DENE  ,
     R   DENLOG,   DNLOG ,   DMINLG,   DMAXLG,
     I   NTE   ,   NNE
       integer nte,nne
       real dtlog,tminlg,tmaxlg,dnlog,dminlg,dmaxlg
       REAL
     R   TE(34)          ,   TELOG(34)       ,   DENE(18)        ,
     R   DENLOG(18)
C-----------------------------------------------------------------------
CL                  C1.3     INTERPOLATED RATE COEFFICIENTS
C VERSION 07/09/84 AVMA JET/OXFORD
       COMMON/COMAT3/
     R   SAL   ,   RAL   ,   DAL
       REAL
     R   SAL(34,28,3)    ,   RAL(18,34,28,3) ,   DAL(18,34,28,3)
c
c      Local Variables:
c
       integer j,j3,j2,jj,ite,ite1,ine,ine1,istate

       real zte,ztelog,zne,znelog,zint,zind,zslg,zra1,zra2,zrlg
       real z1,zda1,zda2,zdlg,z2
c
C-----------------------------------------------------------------------
c
c      DIMENSION  PTE(KMESH), PNE(KMESH), POUT(4,28,3,KMESH)
C
C-----------------------------------------------------------------------
CL              1.         CALCULATE RATE COEFFICIENTS
C
C     LOOP OVER MESH POINTS
         DO 150 J=1,KMESH
C
C     ELECTRON TEMPERATURE IN EV
         ZTE=PTE(J)
C
C     ELECTRON TEMPERATURE OUTSIDE INTERVAL
         IF ( ( ZTE .LT. 1.0 ) .OR. ( ZTE .GT. 1.0E+05 ) ) THEN
           DO 30 J3=1,NIMP
           DO 20 J2=1,NSTATE(J3)
           DO 10 JJ=1,4
           POUT(JJ,J2,J3,J)=-1.0
 10        CONTINUE
 20        CONTINUE
 30        CONTINUE
         ELSE
C
C     ELECTRON TEMPERATURE INSIDE INTERVAL
           ZTELOG=ALOG10(ZTE)
C
C     ELECTRON DENSITY IN CM**-3
           ZNE=PNE(J)
C
C     ELECTRON DENSITY OUTSIDE INTERVAL
           IF ( ( ZNE .LT. 1.0E+10 ) .OR. ( ZNE .GT. 1.0E+15 ) ) THEN
             DO 60 J3=1,NIMP
             DO 50 J2=1,NSTATE(J3)
             DO 40 JJ=1,4
             POUT(JJ,J2,J3,J)=-1.0
 40          CONTINUE
 50          CONTINUE
 60          CONTINUE
           ELSE
C
C     ELECTRON DENSITY INSIDE INTERVAL
             ZNELOG=ALOG10(ZNE)
C
C     INTERVAL FOR TEMPERATURE INTERPOLATION
             ITE = (ZTELOG-TMINLG)/DTLOG + 1.99999
             ! jdemod - need to change from max0 - it implies a certain integer kind these days
             ITE = MAX (ITE,2)
             ITE = MIN (ITE,NTE)
             ITE1= ITE-1
c             write(0,*) 'RRATES:TELOG:',ite,ite1
             ZINT= (TELOG(ITE)-ZTELOG)/(TELOG(ITE)-TELOG(ITE1))
C
C     INTERVAL FOR DENSITY INTERPOLATION
             INE = (ZNELOG-DMINLG)/DNLOG + 1.99999
             ! jdemod - need to change from max0 - it implies a certain integer kind these days
             INE = MAX (INE,2)
             INE = MIN (INE,NNE)
             INE1= INE-1
             ZIND= (DENLOG(INE)-ZNELOG)/(DENLOG(INE)-DENLOG(INE1))
C
C     SCAN OVER IMPURITY SPECIES
             DO 140 J3=1,NIMP
             ISTATE=NSTATE(J3)
C
C     ZERO RATES FOR UNDEFINED IMPURITY MASS
             IF ( NORDER(J3) .EQ. 0 ) THEN
               DO 120 J2=1,ISTATE
               DO 110 JJ=1,4
               POUT(JJ,J2,J3,J)=0.0
 110           CONTINUE
 120           CONTINUE
             ELSE
C
C     RATE COEFFICIENT FOR IONISATION
               DO 130 J2=1,ISTATE
               ZSLG=SAL(ITE1,J2,J3)*ZINT+SAL(ITE,J2,J3)*(1.0-ZINT)
               POUT(1,J2,J3,J)=10.0**ZSLG
C
C     RATE COEFFICIENT FOR RADIATIVE RECOMBINATION
               ZRA1=RAL(INE1,ITE1,J2,J3)*ZINT
     .             +RAL(INE1,ITE,J2,J3)*(1.0-ZINT)
               ZRA2=RAL(INE ,ITE1,J2,J3)*ZINT
     .             +RAL(INE ,ITE,J2,J3)*(1.0-ZINT)
               ZRLG=ZRA1*ZIND+ZRA2*(1.0-ZIND)
               Z1=10.0**ZRLG
               POUT(3,J2,J3,J)=Z1
C
C     RATE COEFFICIENT FOR DIELECTRONIC RECOMBINATION
               ZDA1=DAL(INE1,ITE1,J2,J3)*ZINT
     .             +DAL(INE1,ITE,J2,J3)*(1.0-ZINT)
               ZDA2=DAL(INE ,ITE1,J2,J3)*ZINT
     .             +DAL(INE ,ITE,J2,J3)*(1.0-ZINT)
               ZDLG=ZDA1*ZIND+ZDA2*(1.0-ZIND)
               Z2=10.0**ZDLG
               POUT(4,J2,J3,J)=Z2
C
C     TOTAL RATE COEFFICIENT FOR RECOMBINATION
               POUT(2,J2,J3,J)=Z1+Z2
 130           CONTINUE
             ENDIF
 140         CONTINUE
           ENDIF
         ENDIF
 150     CONTINUE
C
         RETURN
         END
C/ MODULE A3
C
         SUBROUTINE ABUND(PIN,POUT,KMESH,NLNEUT)
         implicit none
c
         integer kmesh
C
C A.3  CALCULATE ABUNDANCE IN CORONAL EQUILIBRIUM
C
C VERSION 02/12/85 AVMA JET/OXFORD
C
C-----------------------------------------------------------------------
CL                  C1.1     IMPURITY IDENTIFICATION
C VERSION 07/09/84 AVMA JET/OXFORD
       COMMON/COMAT1/
     I   NORDER,   NIMP  ,   NSTATE
       integer
     I   NORDER(3)       ,   NSTATE(3),   nimp
C-----------------------------------------------------------------------
C
       real PIN(4,28,3,KMESH), POUT(29,3,KMESH), ZLOG(29)
c       DIMENSION PIN(4,28,3,KMESH), POUT(29,3,KMESH), ZLOG(29)
C
       LOGICAL NLNEUT
c
c      Local variables:
c
       integer j,j3,istate,istat1,j2,jp1
       real    zfact,zdenom
C
C-----------------------------------------------------------------------
CL              1.         CALCULATE FRACTION OF EACH STATE
C
C     LOOP OVER MESH POINTS
         DO 190 J=1,KMESH
C
C     SCAN OVER IMPURITY SPECIES
         DO 180 J3=1,NIMP
         ISTATE=NSTATE(J3)
         ISTAT1=ISTATE+1
C
C     ZERO FRACTIONS FOR UNDEFINED IMPURITY MASS
         IF ( NORDER(J3) .EQ. 0 ) THEN
           DO 110 J2=1,ISTAT1
           POUT(J2,J3,J)=0.0
 110       CONTINUE
         ELSE
C
C     INITIALISE NUMBERS
           ZFACT  =0.0
           ZDENOM =0.0
C
C     FRACTIONS OF TOTAL IMPURITY DENSITY
           IF ( NLNEUT ) THEN
             ZLOG(1)=0.0
C
C     FIRST SCAN OVER IMPURITY STATES
             DO 120 J2=1,ISTATE
             JP1=J2+1
             ZLOG(JP1)=ZLOG(J2)+ALOG10(AMAX1(PIN(1,J2,J3,J),1.0E-35))
     .                         -ALOG10(AMAX1(PIN(2,J2,J3,J),1.0E-35))
             ZFACT=AMAX1(ZFACT,ZLOG(JP1))
 120         CONTINUE
C
C     SECOND SCAN OVER IMPURITY STATES
             DO 130 J2=1,ISTAT1
             ZLOG(J2)=AMAX1(ZLOG(J2)-ZFACT,-10.0)
             ZLOG(J2)=10.0**(ZLOG(J2))
             ZDENOM=ZDENOM+ZLOG(J2)
 130         CONTINUE
C
C     FRACTIONAL ABUNDANCE FOR EACH STATE
             DO 140 J2=1,ISTAT1
             POUT(J2,J3,J)=ZLOG(J2)/ZDENOM
 140         CONTINUE
C
C     FRACTIONS OF TOTAL IMPURITY ION DENSITY
           ELSE
             ZLOG(2)=0.0
C
C     FIRST SCAN OVER IMPURITY STATES
             DO 150 J2=2,ISTATE
             JP1=J2+1
             ZLOG(JP1)=ZLOG(J2)+ALOG10(AMAX1(PIN(1,J2,J3,J),1.0E-35))
     .                         -ALOG10(AMAX1(PIN(2,J2,J3,J),1.0E-35))
             ZFACT=AMAX1(ZFACT,ZLOG(JP1))
 150         CONTINUE
C
C     SECOND SCAN OVER IMPURITY STATES
             DO 160 J2=2,ISTAT1
             ZLOG(J2)=AMAX1(ZLOG(J2)-ZFACT,-10.0)
             ZLOG(J2)=10.0**(ZLOG(J2))
             ZDENOM=ZDENOM+ZLOG(J2)
 160         CONTINUE
C
C     FRACTIONAL ABUNDANCE FOR EACH STATE
             POUT(1,J3,J)=0.0
             DO 170 J2=2,ISTAT1
             POUT(J2,J3,J)=ZLOG(J2)/ZDENOM
 170         CONTINUE
           ENDIF
         ENDIF
 180     CONTINUE
 190     CONTINUE
C
         RETURN
         END
C/ MODULE A4
C
       SUBROUTINE RDSHRT(PTE,PNE,PNZ,PDVOL,POUT,KMESH)
       implicit none
c
       integer kmesh
C
C A.4  CALCULATE IMPURITY RADIATION (SHORT VERSION)
C
C VERSION 02/12/85 AVMA JET/OXFORD
C
C-----------------------------------------------------------------------
CL                  C1.1     IMPURITY IDENTIFICATION
C VERSION 07/09/84 AVMA JET/OXFORD
       COMMON/COMAT1/
     I   NORDER,   NIMP  ,   NSTATE
       integer
     I   NORDER(3)       ,   NSTATE(3),  nimp
C-----------------------------------------------------------------------
CL                  C1.2     INTERPOLATION VARIABLES
C VERSION 07/09/84 AVMA JET/OXFORD
       COMMON/COMAT2/
     R   TE    ,   TELOG ,   DTLOG ,   TMINLG,   TMAXLG,   DENE  ,
     R   DENLOG,   DNLOG ,   DMINLG,   DMAXLG,
     I   NTE   ,   NNE
       integer nte,nne
       real dtlog,tminlg,tmaxlg,dnlog,dminlg,dmaxlg
       real
     R   TE(34)          ,   TELOG(34)       ,   DENE(18)        ,
     R   DENLOG(18)
C-----------------------------------------------------------------------
CL                  C1.3     INTERPOLATED RATE COEFFICIENTS
C VERSION 07/09/84 AVMA JET/OXFORD
       COMMON/COMAT3/
     R   SAL   ,   RAL   ,   DAL
       real
     R   SAL(34,28,3)    ,   RAL(18,34,28,3) ,   DAL(18,34,28,3)
C-----------------------------------------------------------------------
CL                  C1.4     INTERPOLATED RADIATION COEFFICIENTS
C VERSION 20/11/85 AVMA JET/OXFORD
       COMMON/COMAT4/
     R   RLINE ,   BREMS ,   RAD   ,   DIE
       real
     R   RLINE(34,29,3)  ,   BREMS(34,28,3)  ,   RAD(18,34,28,3),
     R   DIE(18,34,28,3)
C-----------------------------------------------------------------------
CL                  C1.5     INTERPOLATED IONISATION LOSS COEFFICIENTS
C VERSION 20/11/85 AVMA JET/OXFORD
       COMMON/COMAT5/
     R   RION  ,   P
       real
     R   RION(34,28,3)   ,   P(28,3)
C-----------------------------------------------------------------------
C
      real PTE(1), PNE(1), PNZ(29,3,KMESH), PDVOL(1)
C
      real POUT(4,3,KMESH)
C
c
c     Local Variables:
c
      integer j,j3,jj,ite,ite1,ine,ine1,istate,j2,j2m1

      real zte,ztelog,zne,znelog,zdvol,zint,zind,znz,zllg,zl,z
      real zrad,zion,znz1,zra1,zra2,zrlg,zr,zda1,zda2,zdlg,zd,zblg,zb
      real zilg,zi,z1,z2,zz1,zz2
c
C-----------------------------------------------------------------------
CL              1.         CALCULATE RADIATION
C
C     LOOP OVER MESH POINTS
         DO 140 J=1,KMESH
C
C     ELECTRON TEMPERATURE IN EV
         ZTE=PTE(J)
C
C     ELECTRON TEMPERATURE OUTSIDE INTERVAL
         IF ( ( ZTE .LT. 1.0 ) .OR. ( ZTE .GT. 1.0E+05 ) ) THEN
           DO 20 J3=1,NIMP
           DO 10 JJ=1,4
           POUT(JJ,J3,J)=-1.0
 10        CONTINUE
 20        CONTINUE
         ELSE
C
C     ELECTRON TEMPERATURE INSIDE INTERVAL
           ZTELOG=ALOG10(ZTE)
C
C     ELECTRON DENSITY IN CM**-3
           ZNE=PNE(J)
C
C     ELECTRON DENSITY OUTSIDE INTERVAL
           IF ( ( ZNE .LT. 1.0E+10 ) .OR. ( ZNE .GT. 1.0E+15 ) ) THEN
             DO 40 J3=1,NIMP
             DO 30 JJ=1,4
             POUT(JJ,J3,J)=-1.0
 30          CONTINUE
 40          CONTINUE
           ELSE
C
C     ELECTRON DENSITY INSIDE INTERVAL
             ZNELOG=ALOG10(ZNE)
C
C     VOLUME FACTOR DVOL IN CM**3
             ZDVOL=PDVOL(J)
C
C     INTERVAL FOR TEMPERATURE INTERPOLATION
             ITE = (ZTELOG-TMINLG)/DTLOG + 1.99999
             ! jdemod - need to change from max0 - it implies a certain integer kind these days
             ITE = MAX (ITE,2)
             ITE = MIN (ITE,NTE)
             ITE1= ITE-1

             ZINT= (TELOG(ITE)-ZTELOG)/(TELOG(ITE)-TELOG(ITE1))
C
C     INTERVAL FOR DENSITY INTERPOLATION
             INE = (ZNELOG-DMINLG)/DNLOG + 1.99999
             ! jdemod - need to change from max0 - it implies a certain integer kind these days
             INE = MAX (INE,2)
             INE = MIN (INE,NNE)
             INE1= INE-1
             ZIND= (DENLOG(INE)-ZNELOG)/(DENLOG(INE)-DENLOG(INE1))
C
C     SCAN OVER IMPURITY SPECIES
             DO 130 J3=1,NIMP
C
C     ZERO ENERGY TERMS FOR UNDEFINED IMPURITY MASS
             IF ( NORDER(J3) .EQ. 0 ) THEN
               DO 110 JJ=1,4
               POUT(JJ,J3,J)=0.0
 110           CONTINUE
             ELSE
C
C     NEUTRAL IMPURITY
               ISTATE=NSTATE(J3)+1
               ZNZ=PNZ(1,J3,J)
C
C     TOTAL RADIATION AND ENERGY LOSS
               ZLLG=RLINE(ITE1,1,J3)*ZINT+RLINE(ITE,1,J3)*(1.0-ZINT)
               ZL=10.0**ZLLG*ZNZ*ZNE
               POUT(1,J3,J)=ZL
               POUT(2,J3,J)=ZL
C
C     VOLUME INTEGRATED TOTAL RADIATION AND ENERGY LOSS
               Z=ZDVOL*ZL
               ZRAD=Z
               ZION=Z
C
C     SCAN OVER IONISED IMPURITY STATES
               DO 120 J2=2,ISTATE
               J2M1=J2-1
               ZNZ =PNZ(J2  ,J3,J)
               ZNZ1=PNZ(J2M1,J3,J)
C
C     LINE RADIATION
               ZLLG=RLINE(ITE1,J2,J3)*ZINT+RLINE(ITE,J2,J3)*(1.0-ZINT)
               ZL=10.0**ZLLG*ZNZ*ZNE
C
C     RADIATIVE RECOMBINATION RADIATION
               ZRA1=RAD(INE1,ITE1,J2M1,J3)*ZINT
     .             +RAD(INE1,ITE,J2M1,J3)*(1.0-ZINT)
               ZRA2=RAD(INE ,ITE1,J2M1,J3)*ZINT
     .             +RAD(INE ,ITE,J2M1,J3)*(1.0-ZINT)
               ZRLG=ZRA1*ZIND+ZRA2*(1.0-ZIND)
               ZR=10.0**ZRLG*ZNZ*ZNE
C
C     DIELECTRONIC RECOMBINATION RADIATION
               ZDA1=DIE(INE1,ITE1,J2M1,J3)*ZINT
     .             +DIE(INE1,ITE,J2M1,J3)*(1.0-ZINT)
               ZDA2=DIE(INE ,ITE1,J2M1,J3)*ZINT
     .             +DIE(INE ,ITE,J2M1,J3)*(1.0-ZINT)
               ZDLG=ZDA1*ZIND+ZDA2*(1.0-ZIND)
               ZD=10.0**ZDLG*ZNZ*ZNE
C
C     BREMSSTRAHLUNG
             ZBLG=BREMS(ITE1,J2M1,J3)*ZINT+BREMS(ITE,J2M1,J3)*(1.0-ZINT)
               ZB=10.0**ZBLG*ZNZ*ZNE
C
C     IONISATION ENERGY LOSS
               ZILG=RION(ITE1,J2M1,J3)*ZINT+RION(ITE,J2M1,J3)*(1.0-ZINT)
               ZI=10.0**ZILG*ZNZ1*ZNE
C
C     RATE COEFFICIENT FOR RADIATIVE RECOMBINATION
               ZRA1=RAL(INE1,ITE1,J2M1,J3)*ZINT
     .             +RAL(INE1,ITE,J2M1,J3)*(1.0-ZINT)
               ZRA2=RAL(INE ,ITE1,J2M1,J3)*ZINT
     .             +RAL(INE ,ITE,J2M1,J3)*(1.0-ZINT)
               ZRLG=ZRA1*ZIND+ZRA2*(1.0-ZIND)
               Z1=10.0**ZRLG
C
C     RATE COEFFICIENT FOR DIELECTRONIC RECOMBINATION
               ZDA1=DAL(INE1,ITE1,J2M1,J3)*ZINT
     .             +DAL(INE1,ITE,J2M1,J3)*(1.0-ZINT)
               ZDA2=DAL(INE ,ITE1,J2M1,J3)*ZINT
     .             +DAL(INE ,ITE,J2M1,J3)*(1.0-ZINT)
               ZDLG=ZDA1*ZIND+ZDA2*(1.0-ZIND)
               Z2=10.0**ZDLG
C
C     TOTAL RADIATION AND ENERGY LOSS
               ZZ1=ZL+ZR+ZD+ZB
               ZZ2=ZZ1+ZI-1.6E-19*P(J2M1,J3)*ZNZ*ZNE*(Z1+Z2)
               POUT(1,J3,J)=POUT(1,J3,J)+ZZ1
               POUT(2,J3,J)=POUT(2,J3,J)+ZZ2
 120           CONTINUE
C
C     VOLUME INTEGRATED TOTAL RADIATION AND ENERGY LOSS
               ZRAD=ZRAD+ZDVOL*POUT(1,J3,J)
               ZION=ZION+ZDVOL*POUT(2,J3,J)
               POUT(3,J3,J)=ZRAD
               POUT(4,J3,J)=ZION
             ENDIF
 130         CONTINUE
           ENDIF
         ENDIF
 140     CONTINUE
C
         RETURN
         END
C/ MODULE A5
C
       SUBROUTINE RDLONG(PTE,PNE,PNZ,PDVOL,POUT,KMESH)
       implicit none
c
       integer kmesh
C
C A.5  CALCULATE IMPURITY RADIATION
C
C VERSION 02/12/85 AVMA JET/OXFORD
C
C-----------------------------------------------------------------------
CL                  C1.1     IMPURITY IDENTIFICATION
C VERSION 07/09/84 AVMA JET/OXFORD
       COMMON/COMAT1/
     I   NORDER,   NIMP  ,   NSTATE
       integer
     I   NORDER(3)       ,   NSTATE(3),   nimp
C-----------------------------------------------------------------------
CL                  C1.2     INTERPOLATION VARIABLES
C VERSION 07/09/84 AVMA JET/OXFORD
       COMMON/COMAT2/
     R   TE    ,   TELOG ,   DTLOG ,   TMINLG,   TMAXLG,   DENE  ,
     R   DENLOG,   DNLOG ,   DMINLG,   DMAXLG,
     I   NTE   ,   NNE
       integer nte,nne
       real dtlog,tminlg,tmaxlg,dnlog,dminlg,dmaxlg
       real
     R   TE(34)          ,   TELOG(34)       ,   DENE(18)        ,
     R   DENLOG(18)
C-----------------------------------------------------------------------
CL                  C1.3     INTERPOLATED RATE COEFFICIENTS
C VERSION 07/09/84 AVMA JET/OXFORD
       COMMON/COMAT3/
     R   SAL   ,   RAL   ,   DAL
       real
     R   SAL(34,28,3)    ,   RAL(18,34,28,3) ,   DAL(18,34,28,3)
C-----------------------------------------------------------------------
CL                  C1.4     INTERPOLATED RADIATION COEFFICIENTS
C VERSION 20/11/85 AVMA JET/OXFORD
       COMMON/COMAT4/
     R   RLINE ,   BREMS ,   RAD   ,   DIE
       real
     R   RLINE(34,29,3)  ,   BREMS(34,28,3)  ,   RAD(18,34,28,3),
     R   DIE(18,34,28,3)
C-----------------------------------------------------------------------
CL                  C1.5     INTERPOLATED IONISATION LOSS COEFFICIENTS
C VERSION 20/11/85 AVMA JET/OXFORD
       COMMON/COMAT5/
     R   RION  ,   P
       real
     R   RION(34,28,3)   ,   P(28,3)
C-----------------------------------------------------------------------
C
      real PTE(1), PNE(1), PNZ(29,3,KMESH), PDVOL(1)
C
      real POUT(14,29,3,KMESH)
C
      real ZRAD(7)
c
c     Local Variables:
c
       integer j,j3,j2,jj,ite,ite1,ine,ine1,j2m1,istate,jj7

       real zte,ztelog,zne,znelog,zdvol,zint,zind,znz,zllg,zl,z
       real znz1,zra1,zra2,zrlg
       real zr,zda1,zda2,zdlg,zd,zblg,zb,zilg,zi,z1,z2,zz1,zz2
C
C-----------------------------------------------------------------------
CL              1.         CALCULATE RADIATION
C
C     LOOP OVER MESH POINTS
         DO 160 J=1,KMESH
C
C     ELECTRON TEMPERATURE IN EV
         ZTE=PTE(J)
C
C     ELECTRON TEMPERATURE OUTSIDE INTERVAL
         IF ( ( ZTE .LT. 1.0 ) .OR. ( ZTE .GT. 1.0E+05 ) ) THEN
           DO 30 J3=1,NIMP
           ISTATE=NSTATE(J3)+1
           DO 20 J2=1,ISTATE
           DO 10 JJ=1,14
           POUT(JJ,J2,J3,J)=-1.0
 10        CONTINUE
 20        CONTINUE
 30        CONTINUE
         ELSE
C
C     ELECTRON TEMPERATURE INSIDE INTERVAL
           ZTELOG=ALOG10(ZTE)
C
C     ELECTRON DENSITY IN CM**-3
           ZNE=PNE(J)
C
C     ELECTRON DENSITY OUTSIDE INTERVAL
           IF ( ( ZNE .LT. 1.0E+10 ) .OR. ( ZNE .GT. 1.0E+15 ) ) THEN
             DO 60 J3=1,NIMP
             ISTATE=NSTATE(J3)+1
             DO 50 J2=1,ISTATE
             DO 40 JJ=1,14
             POUT(JJ,J2,J3,J)=-1.0
 40          CONTINUE
 50          CONTINUE
 60          CONTINUE
           ELSE
C
C     ELECTRON DENSITY INSIDE INTERVAL
             ZNELOG=ALOG10(ZNE)
C
C     VOLUME FACTOR DVOL IN CM**3
             ZDVOL=PDVOL(J)
C
C     INTERVAL FOR TEMPERATURE INTERPOLATION
             ITE = (ZTELOG-TMINLG)/DTLOG + 1.99999
             ! jdemod - need to change from max0 - it implies a certain integer kind these days
             !ITE = MAX0 (ITE,2)
             !ITE = MIN0 (ITE,NTE)
             ITE = MAX (ITE,2)
             ITE = MIN (ITE,NTE)
             ITE1= ITE-1

             if (ite.lt.2) then 
                write (0,*) 'RDLONG:TELOG:ITE,ITE1:',ite,ite1,nte,
     >                     max(ite,2),min(ite,nte)
             endif

             ZINT= (TELOG(ITE)-ZTELOG)/(TELOG(ITE)-TELOG(ITE1))
C
C     INTERVAL FOR DENSITY INTERPOLATION
             INE = (ZNELOG-DMINLG)/DNLOG + 1.99999
             ! jdemod - need to change from max0 - it implies a certain integer kind these days
             INE = MAX (INE,2)
             INE = MIN (INE,NNE)
             INE1= INE-1
             ZIND= (DENLOG(INE)-ZNELOG)/(DENLOG(INE)-DENLOG(INE1))
C
C     SCAN OVER IMPURITY SPECIES
             DO 150 J3=1,NIMP
             ISTATE=NSTATE(J3)+1
C
C     ZERO ENERGY TERMS FOR UNDEFINED IMPURITY MASS
             IF ( NORDER(J3) .EQ. 0 ) THEN
               DO 120 J2=1,ISTATE
               DO 110 JJ=1,14
               POUT(JJ,J2,J3,J)=0.0
 110           CONTINUE
 120           CONTINUE
             ELSE
C
C     NEUTRAL IMPURITY
               ZNZ=PNZ(1,J3,J)
C
C     INDIVIDUAL ENERGY TERMS
               ZLLG=RLINE(ITE1,1,J3)*ZINT+RLINE(ITE,1,J3)*(1.0-ZINT)
               ZL=10.0**ZLLG*ZNZ*ZNE
               POUT(1,1,J3,J)=ZL
               POUT(2,1,J3,J)=ZL
               POUT(3,1,J3,J)=ZL
               POUT(4,1,J3,J)=0.0
               POUT(5,1,J3,J)=0.0
               POUT(6,1,J3,J)=0.0
               POUT(7,1,J3,J)=0.0
C
C     VOLUME INTEGRATED TERMS INITIALISED
               Z=ZDVOL*ZL
               POUT(8,1,J3,J) =Z
               POUT(9,1,J3,J) =Z
               POUT(10,1,J3,J)=Z
               POUT(11,1,J3,J)=0.0
               POUT(12,1,J3,J)=0.0
               POUT(13,1,J3,J)=0.0
               POUT(14,1,J3,J)=0.0
C
               ZRAD(1)=Z
               ZRAD(2)=Z
               ZRAD(3)=Z
               ZRAD(4)=0.0
               ZRAD(5)=0.0
               ZRAD(6)=0.0
               ZRAD(7)=0.0
C
C     SCAN OVER IONISED IMPURITY STATES
               DO 140 J2=2,ISTATE
               J2M1=J2-1
               ZNZ =PNZ(J2  ,J3,J)
               ZNZ1=PNZ(J2M1,J3,J)
C
C     LINE RADIATION
               ZLLG=RLINE(ITE1,J2,J3)*ZINT+RLINE(ITE,J2,J3)*(1.0-ZINT)
               ZL=10.0**ZLLG*ZNZ*ZNE
               POUT(3,J2,J3,J)=ZL
C
C     RADIATIVE RECOMBINATION RADIATION
               ZRA1=RAD(INE1,ITE1,J2M1,J3)*ZINT
     .             +RAD(INE1,ITE,J2M1,J3)*(1.0-ZINT)
               ZRA2=RAD(INE ,ITE1,J2M1,J3)*ZINT
     .             +RAD(INE ,ITE,J2M1,J3)*(1.0-ZINT)
               ZRLG=ZRA1*ZIND+ZRA2*(1.0-ZIND)
               ZR=10.0**ZRLG*ZNZ*ZNE
               POUT(4,J2,J3,J)=ZR
C
C     DIELECTRONIC RECOMBINATION RADIATION
               ZDA1=DIE(INE1,ITE1,J2M1,J3)*ZINT
     .             +DIE(INE1,ITE,J2M1,J3)*(1.0-ZINT)
               ZDA2=DIE(INE ,ITE1,J2M1,J3)*ZINT
     .             +DIE(INE ,ITE,J2M1,J3)*(1.0-ZINT)
               ZDLG=ZDA1*ZIND+ZDA2*(1.0-ZIND)
               ZD=10.0**ZDLG*ZNZ*ZNE
               POUT(5,J2,J3,J)=ZD
C
C     BREMSSTRAHLUNG
               ZBLG=BREMS(ITE1,J2M1,J3)*ZINT
     .             +BREMS(ITE,J2M1,J3)*(1.0-ZINT)
               ZB=10.0**ZBLG*ZNZ*ZNE
               POUT(6,J2,J3,J)=ZB
C
C     IONISATION ENERGY LOSS
               ZILG=RION(ITE1,J2M1,J3)*ZINT+RION(ITE,J2M1,J3)*(1.0-ZINT)
               ZI=10.0**ZILG*ZNZ1*ZNE
               POUT(7,J2,J3,J)=ZI
C
C     RATE COEFFICIENT FOR RADIATIVE RECOMBINATION
               ZRA1=RAL(INE1,ITE1,J2M1,J3)*ZINT
     .             +RAL(INE1,ITE,J2M1,J3)*(1.0-ZINT)
               ZRA2=RAL(INE ,ITE1,J2M1,J3)*ZINT
     .             +RAL(INE ,ITE,J2M1,J3)*(1.0-ZINT)
               ZRLG=ZRA1*ZIND+ZRA2*(1.0-ZIND)
               Z1=10.0**ZRLG
C
C     RATE COEFFICIENT FOR DIELECTRONIC RECOMBINATION
               ZDA1=DAL(INE1,ITE1,J2M1,J3)*ZINT
     .             +DAL(INE1,ITE,J2M1,J3)*(1.0-ZINT)
               ZDA2=DAL(INE ,ITE1,J2M1,J3)*ZINT
     .             +DAL(INE ,ITE,J2M1,J3)*(1.0-ZINT)
               ZDLG=ZDA1*ZIND+ZDA2*(1.0-ZIND)
               Z2=10.0**ZDLG
C
C     TOTAL RADIATION AND ENERGY LOSS
               ZZ1=ZL+ZR+ZD+ZB
               ZZ2=ZZ1+ZI-1.6E-19*P(J2M1,J3)*ZNZ*ZNE*(Z1+Z2)
               POUT(1,J2,J3,J)=ZZ1
               POUT(2,J2,J3,J)=ZZ2
C
C     VOLUME INTEGRATED TERMS
               DO 130 JJ=1,7
               JJ7=JJ+7
               ZRAD(JJ)=ZRAD(JJ)+ZDVOL*POUT(JJ,J2,J3,J)
               POUT(JJ7,J2,J3,J)=ZRAD(JJ)
 130           CONTINUE
 140           CONTINUE
             ENDIF
 150         CONTINUE
           ENDIF
         ENDIF
 160     CONTINUE
C
         RETURN
         END
