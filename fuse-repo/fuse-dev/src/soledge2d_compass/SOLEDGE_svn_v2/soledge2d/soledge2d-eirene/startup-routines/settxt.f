c  bug fix:  text(71-13) --> text(71)
c  17.3.06: txttal and txttlw added for additional tallies
C
      SUBROUTINE EIRENE_SETTXT
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CTEXT
      USE EIRMOD_COUTAU
 
      IMPLICIT NONE
 
      INTEGER :: IATM, IION, IPLS, IMOL, ISPZ, IPHOT, I, J, N1,
     .           N2, N3, N4, N5, N6, N7, N8, N9, N10, N11
      CHARACTER(24) :: TEXT24
      CHARACTER(72) :: TEXT72
 
      TXTTAL(1,1)='PARTICLE DENSITY (ATOMS)                         '
      TXTTAL(1,2)='PARTICLE DENSITY (MOLECULES)                     '
      TXTTAL(1,3)='PARTICLE DENSITY (TEST IONS)                     '
      TXTTAL(1,4)='PARTICLE DENSITY (PHOTONS)                       '
      TXTTAL(1,5)='ENERGY DENSITY (ATOMS)                           '
      TXTTAL(1,6)='ENERGY DENSITY (MOLECULES)                       '
      TXTTAL(1,7)='ENERGY DENSITY (TEST IONS)                       '
      TXTTAL(1,8)='ENERGY DENSITY (PHOTONS)                         '
      TXTTAL(1,9)=
     . 'PARTICLE SOURCE (ELECTRONS) FROM ATOM-PLASMA INTERACTION    '
      TXTTAL(1,10)=
     . 'PARTICLE SOURCE (ATOMS) FROM ATOM-PLASMA INTERACTION        '
      TXTTAL(1,11)=
     . 'PARTICLE SOURCE (MOLECULES) FROM ATOM-PLASMA INTERACTION    '
      TXTTAL(1,12)=
     . 'PARTICLE SOURCE (TEST IONS) FROM ATOM-PLASMA INTERACTION    '
      TXTTAL(1,13)=
     . 'PARTICLE SOURCE (PHOTONS) FROM ATOM-PLASMA INTERACTION    '
      TXTTAL(1,14)=
     . 'PARTICLE SOURCE (BULK IONS) FROM ATOM-PLASMA INTERACTION    '
      TXTTAL(1,15)=
     . 'PARTICLE SOURCE (ELECTRONS) FROM MOLECULE-PLASMA INTERACTION'
      TXTTAL(1,16)=
     . 'PARTICLE SOURCE (ATOMS) FROM MOLECULE-PLASMA INTERACTION    '
      TXTTAL(1,17)=
     . 'PARTICLE SOURCE (MOLECULES) FROM MOLECULE-PLASMA INTERACTION'
      TXTTAL(1,18)=
     . 'PARTICLE SOURCE (TEST IONS) FROM MOLECULE-PLASMA INTERACTION'
      TXTTAL(1,19)=
     . 'PARTICLE SOURCE (PHOTONS) FROM MOLECULE-PLASMA INTERACTION'
      TXTTAL(1,20)=
     . 'PARTICLE SOURCE (BULK IONS) FROM MOLECULE-PLASMA INTERACTION'
      TXTTAL(1,21)=
     . 'PARTICLE SOURCE (ELECTRONS) FROM TEST ION-PLASMA INTERACTION'
      TXTTAL(1,22)=
     . 'PARTICLE SOURCE (ATOMS) FROM TEST ION-PLASMA INTERACTION    '
      TXTTAL(1,23)=
     . 'PARTICLE SOURCE (MOLECULES) FROM TEST ION-PLASMA INTERACTION'
      TXTTAL(1,24)=
     . 'PARTICLE SOURCE (TEST IONS) FROM TEST ION-PLASMA INTERACTION'
      TXTTAL(1,25)=
     . 'PARTICLE SOURCE (PHOTONS) FROM TEST ION-PLASMA INTERACTION  '
      TXTTAL(1,26)=
     . 'PARTICLE SOURCE (BULK IONS) FROM TEST ION-PLASMA INTERACTION'
      TXTTAL(1,27)=
     . 'PARTICLE SOURCE (ELECTRONS) FROM PHOTON-PLASMA INTERACTION  '
      TXTTAL(1,28)=
     . 'PARTICLE SOURCE (ATOMS) FROM PHOTON-PLASMA INTERACTION      '
      TXTTAL(1,29)=
     . 'PARTICLE SOURCE (MOLECULES) FROM PHOTON-PLASMA INTERACTION  '
      TXTTAL(1,30)=
     . 'PARTICLE SOURCE (TEST IONS) FROM PHOTON-PLASMA INTERACTION  '
      TXTTAL(1,31)=
     . 'PARTICLE SOURCE (PHOTONS) FROM PHOTON-PLASMA INTERACTION    '
      TXTTAL(1,32)=
     . 'PARTICLE SOURCE (BULK IONS) FROM PHOTON-PLASMA INTERACTION  '
      TXTTAL(1,33)=
     . 'ENERGY SOURCE (ELECTRONS) FROM ATOM-PLASMA INTERACTION      '
      TXTTAL(1,34)=
     . 'ENERGY SOURCE (ATOMS) FROM ATOM-PLASMA INTERACTION          '
      TXTTAL(1,35)=
     . 'ENERGY SOURCE (MOLECULES) FROM ATOM-PLASMA INTERACTION      '
      TXTTAL(1,36)=
     . 'ENERGY SOURCE (TEST IONS) FROM ATOM-PLASMA INTERACTION      '
      TXTTAL(1,37)=
     . 'ENERGY SOURCE (PHOTONS) FROM ATOM-PLASMA INTERACTION        '
      TXTTAL(1,38)=
     . 'ENERGY SOURCE (BULK IONS) FROM ATOM-PLASMA INTERACTION      '
      TXTTAL(1,39)=
     . 'ENERGY SOURCE (ELECTRONS) FROM MOLECULE-PLASMA INTERACTION  '
      TXTTAL(1,40)=
     . 'ENERGY SOURCE (ATOMS) FROM MOLECULE-PLASMA INTERACTION      '
      TXTTAL(1,41)=
     . 'ENERGY SOURCE (MOLECULES) FROM MOLECULE-PLASMA INTERACTION  '
      TXTTAL(1,42)=
     . 'ENERGY SOURCE (TEST IONS) FROM MOLECULE-PLASMA INTERACTION  '
      TXTTAL(1,43)=
     . 'ENERGY SOURCE (PHOTONS) FROM MOLECULE-PLASMA INTERACTION    '
      TXTTAL(1,44)=
     . 'ENERGY SOURCE (BULK IONS) FROM MOLECULE-PLASMA INTERACTION  '
      TXTTAL(1,45)=
     . 'ENERGY SOURCE (ELECTRONS) FROM TEST ION-PLASMA INTERACTION  '
      TXTTAL(1,46)=
     . 'ENERGY SOURCE (ATOMS) FROM TEST ION-PLASMA INTERACTION      '
      TXTTAL(1,47)=
     . 'ENERGY SOURCE (MOLECULES) FROM TEST ION-PLASMA INTERACTION  '
      TXTTAL(1,48)=
     . 'ENERGY SOURCE (TEST IONS) FROM TEST ION-PLASMA INTERACTION  '
      TXTTAL(1,49)=
     . 'ENERGY SOURCE (PHOTONS) FROM TEST ION-PLASMA INTERACTION    '
      TXTTAL(1,50)=
     . 'ENERGY SOURCE (BULK IONS) FROM TEST ION-PLASMA INTERACTION  '
      TXTTAL(1,51)=
     . 'ENERGY SOURCE (ELECTRONS) FROM PHOTON-PLASMA INTERACTION    '
      TXTTAL(1,52)=
     . 'ENERGY SOURCE (ATOMS) FROM PHOTON-PLASMA INTERACTION        '
      TXTTAL(1,53)=
     . 'ENERGY SOURCE (MOLECULES) FROM PHOTON-PLASMA INTERACTION    '
      TXTTAL(1,54)=
     . 'ENERGY SOURCE (TEST IONS) FROM PHOTON-PLASMA INTERACTION    '
      TXTTAL(1,55)=
     . 'ENERGY SOURCE (PHOTONS) FROM PHOTON-PLASMA INTERACTION      '
      TXTTAL(1,56)=
     . 'ENERGY SOURCE (BULK IONS) FROM PHOTON-PLASMA INTERACTION    '
C  TALLY NTALA=57 (SEE PARMMOD.F)
C        ADDIT. TRACKLENGTH ESTIMATED TALLIES
C        TXTTAL IS OVERWRITTEN BY INPUT BLOCK 10A
      TXTTAL(1,NTALA)=
     . 'ADDITIONAL TALLIES, TRACKLENGTH ESTIMATOR, SUBR. UPTUSR.F   '
C  TALLY NTALC=58 (SEE PARMMOD.F)
C        ADDIT. COLLISION ESTIMATED TALLIES
C        TXTTAL IS OVERWRITTEN BY INPUT BLOCK 10B
      TXTTAL(1,NTALC)=
     . 'ADDITIONAL TALLIES, COLLISION ESTIMATOR, SUBR. UPCUSR.F     '
C  TALLY NTALT=59 (SEE PARMMOD.F)
C        ADDIT. SNAPSHOT ESTIMATED TALLIES
C        TXTTAL IS OVERWRITTEN BY INPUT BLOCK 13B
      TXTTAL(1,NTALT)=
     . 'ADDITIONAL TALLIES, SNAPSHOT ESTIMATOR, SUBR. UPNUSR.F      '
C  TALLY NTALM=60 (SEE PARMMOD.F)
C        ADDIT. TALLIES FOR INTERFACING TO OTHER CODES
C        TXTTAL MAY BE OVERWRITTEN IN SUBR. INFCOP
      TXTTAL(1,NTALM)=
     . 'ADDITIONAL TALLIES FOR INTERFACING, SUBR. INFCOP.F          '
C  TALLY NTALB=61 (SEE PARMMOD.F)
C        ADDIT. TALLIES FOR ITERATIVE MODE (BGK-ITERATION)
      TXTTAL(1,NTALB)=
     . 'ADDITIONAL TALLIES FOR ITERATIVE MODE, SUBR. UPTBGK.F       '
C  TALLY NTALB=62 (SEE PARMMOD.F)
C        ADDIT. TALLIES, ALGEBRAIC EXPRESSION IN EXISTING TALLIES
C        TXTTAL IS OVERWRITTEN BY INPUT BLOCK 10C
      TXTTAL(1,NTALR)=
     . 'ADDITIONAL TALLIES, ALGEBRAIC EXPRESSIONS, INPUT BLOCK 10C  '
      TXTTAL(1,63)=
     . 'PARTICLE SINK (ATOMS) DUE TO GENERATION LIMIT               '
      TXTTAL(1,64)=
     . 'PARTICLE SINK (MOLECULES) DUE TO GENERATION LIMIT           '
      TXTTAL(1,65)=
     . 'PARTICLE SINK (TEST IONS) DUE TO GENERATION LIMIT           '
      TXTTAL(1,66)=
     . 'PARTICLE SINK (PHOTONS) DUE TO GENERATION LIMIT             '
      TXTTAL(1,67)=
     . 'ENERGY SINK (ATOMS) DUE TO GENERATION LIMIT                 '
      TXTTAL(1,68)=
     . 'ENERGY SINK (MOLECULES) DUE TO GENERATION LIMIT             '
      TXTTAL(1,69)=
     . 'ENERGY SINK (TEST IONS) DUE TO GENERATION LIMIT             '
      TXTTAL(1,70)=
     . 'ENERGY SINK (PHOTONS) DUE TO GENERATION LIMIT               '
      TXTTAL(1,71)=
     . 'MOMENTUM SINK (ATOMS) DUE TO GENERATION LIMIT               '
      TXTTAL(1,72)=
     . 'MOMENTUM SINK (MOLECULES) DUE TO GENERATION LIMIT           '
      TXTTAL(1,73)=
     . 'MOMENTUM SINK (TEST IONS) DUE TO GENERATION LIMIT           '
      TXTTAL(1,74)=
     . 'MOMENTUM SINK (PHOTONS) DUE TO GENERATION LIMIT             '
      TXTTAL(1,75)=
     . 'PRIMARY PARTICLE SOURCE (ATOMS) FROM PLASMA INTERACTIONS    '
      TXTTAL(1,76)=
     . 'PRIMARY PARTICLE SOURCE (MOLECULES) FROM PLASMA INTERACTIONS'
      TXTTAL(1,77)=
     . 'PRIMARY PARTICLE SOURCE (TEST IONS) FROM PLASMA INTERACTIONS'
      TXTTAL(1,78)=
     . 'PRIMARY PARTICLE SOURCE (PHOTONS) FROM PLASMA INTERACTIONS  '
      TXTTAL(1,79)=
     . 'PRIMARY PARTICLE SOURCE (BULK IONS) FROM PLASMA INTERACTIONS'
      TXTTAL(1,80)=
     . 'PRIMARY ENERGY SOURCE (ATOMS) FROM PLASMA INTERACTIONS      '
      TXTTAL(1,81)=
     . 'PRIMARY ENERGY SOURCE (MOLECULES) FROM PLASMA INTERACTIONS  '
      TXTTAL(1,82)=
     . 'PRIMARY ENERGY SOURCE (TEST IONS) FROM PLASMA INTERACTIONS  '
      TXTTAL(1,83)=
     . 'PRIMARY ENERGY SOURCE (PHOTONS) FROM PLASMA INTERACTIONS    '
      TXTTAL(1,84)=
     . 'PRIMARY ENERGY SOURCE (BULK IONS) FROM PLASMA INTERACTIONS  '
      TXTTAL(1,85)=
     . 'MOMENTUM DENSITY, X-DIRECTION (ATOMS)                       '
      TXTTAL(1,86)=
     . 'MOMENTUM DENSITY, X-DIRECTION (MOLECULES)                   '
      TXTTAL(1,87)=
     . 'MOMENTUM DENSITY, X-DIRECTION (TEST IONS)                   '
      TXTTAL(1,88)=
     . 'MOMENTUM DENSITY, X-DIRECTION (PHOTONS)                     '
      TXTTAL(1,89)=
     . 'MOMENTUM DENSITY, Y-DIRECTION (ATOMS)                       '
      TXTTAL(1,90)=
     . 'MOMENTUM DENSITY, Y-DIRECTION (MOLECULES)                   '
      TXTTAL(1,91)=
     . 'MOMENTUM DENSITY, Y-DIRECTION (TEST IONS)                   '
      TXTTAL(1,92)=
     . 'MOMENTUM DENSITY, Y-DIRECTION (PHOTONS)                     '
      TXTTAL(1,93)=
     . 'MOMENTUM DENSITY, Z-DIRECTION (ATOMS)                       '
      TXTTAL(1,94)=
     . 'MOMENTUM DENSITY, Z-DIRECTION (MOLECULES)                   '
      TXTTAL(1,95)=
     . 'MOMENTUM DENSITY, Z-DIRECTION (TEST IONS)                   '
      TXTTAL(1,96)=
     . 'MOMENTUM DENSITY, Z-DIRECTION (PHOTONS)                     '
      TXTTAL(1,97)=
     . 'MOMENTUM SOURCE (BULK IONS) FROM ATOM-PLASMA INTERACTION    '
      TXTTAL(1,98)=
     . 'MOMENTUM SOURCE (BULK IONS) FROM MOLECULE-PLASMA INTERACTION'
      TXTTAL(1,99)=
     . 'MOMENTUM SOURCE (BULK IONS) FROM TEST ION-PLASMA INTERACTION'
      TXTTAL(1,100)=
     . 'MOMENTUM SOURCE (BULK IONS) FROM PHOTON-PLASMA INTERACTION  '
C
      DO 1 J=1,NTALV
        DO 1 I=2,N1MX
          TEXT72=TXTTAL(1,J)
          TXTTAL(I,J)=TEXT72
1     CONTINUE
C
      TXTUNT(1,1)='CM**-3                  '
      TXTUNT(1,2)='CM**-3                  '
      TXTUNT(1,3)='CM**-3                  '
      TXTUNT(1,4)='CM**-3                  '
      TXTUNT(1,5)='EV*CM**-3               '
      TXTUNT(1,6)='EV*CM**-3               '
      TXTUNT(1,7)='EV*CM**-3               '
      TXTUNT(1,8)='EV*CM**-3               '
      TXTUNT(1,9)='AMP*CM**-3              '
      TXTUNT(1,10)='AMP*CM**-3              '
      TXTUNT(1,11)='AMP*CM**-3              '
      TXTUNT(1,12)='AMP*CM**-3              '
      TXTUNT(1,13)='AMP*CM**-3              '
      TXTUNT(1,14)='AMP*CM**-3              '
      TXTUNT(1,15)='AMP*CM**-3              '
      TXTUNT(1,16)='AMP*CM**-3              '
      TXTUNT(1,17)='AMP*CM**-3              '
      TXTUNT(1,18)='AMP*CM**-3              '
      TXTUNT(1,19)='AMP*CM**-3              '
      TXTUNT(1,20)='AMP*CM**-3              '
      TXTUNT(1,21)='AMP*CM**-3              '
      TXTUNT(1,22)='AMP*CM**-3              '
      TXTUNT(1,23)='AMP*CM**-3              '
      TXTUNT(1,24)='AMP*CM**-3              '
      TXTUNT(1,25)='AMP*CM**-3              '
      TXTUNT(1,26)='AMP*CM**-3              '
      TXTUNT(1,27)='AMP*CM**-3              '
      TXTUNT(1,28)='AMP*CM**-3              '
      TXTUNT(1,29)='AMP*CM**-3              '
      TXTUNT(1,30)='AMP*CM**-3              '
      TXTUNT(1,31)='AMP*CM**-3              '
      TXTUNT(1,32)='AMP*CM**-3              '
      TXTUNT(1,33)='WATT*CM**-3             '
      TXTUNT(1,34)='WATT*CM**-3             '
      TXTUNT(1,35)='WATT*CM**-3             '
      TXTUNT(1,36)='WATT*CM**-3             '
      TXTUNT(1,37)='WATT*CM**-3             '
      TXTUNT(1,38)='WATT*CM**-3             '
      TXTUNT(1,39)='WATT*CM**-3             '
      TXTUNT(1,40)='WATT*CM**-3             '
      TXTUNT(1,41)='WATT*CM**-3             '
      TXTUNT(1,42)='WATT*CM**-3             '
      TXTUNT(1,43)='WATT*CM**-3             '
      TXTUNT(1,44)='WATT*CM**-3             '
      TXTUNT(1,45)='WATT*CM**-3             '
      TXTUNT(1,46)='WATT*CM**-3             '
      TXTUNT(1,47)='WATT*CM**-3             '
      TXTUNT(1,48)='WATT*CM**-3             '
      TXTUNT(1,49)='WATT*CM**-3             '
      TXTUNT(1,50)='WATT*CM**-3             '
      TXTUNT(1,51)='WATT*CM**-3             '
      TXTUNT(1,52)='WATT*CM**-3             '
      TXTUNT(1,53)='WATT*CM**-3             '
      TXTUNT(1,54)='WATT*CM**-3             '
      TXTUNT(1,55)='WATT*CM**-3             '
      TXTUNT(1,56)='WATT*CM**-3             '
      TXTUNT(1,NTALA)='TO BE READ              '
      TXTUNT(1,NTALC)='TO BE READ              '
      TXTUNT(1,NTALM)='TO BE DEFINED IN INFCOP '
      TXTUNT(1,NTALR)='TO BE READ              '
      TXTUNT(1,NTALB)='TO BE DEFINED IN BGK    '
C  GENERATION LIMIT TALLIES
      TXTUNT(1,63)='AMP*CM**-3              '
      TXTUNT(1,64)='AMP*CM**-3              '
      TXTUNT(1,65)='AMP*CM**-3              '
      TXTUNT(1,66)='AMP*CM**-3              '
      TXTUNT(1,67)='EV*CM**-3               '
      TXTUNT(1,68)='EV*CM**-3               '
      TXTUNT(1,69)='EV*CM**-3               '
      TXTUNT(1,70)='EV*CM**-3               '
      TXTUNT(1,71)='CM/S*CM**-3             '
      TXTUNT(1,72)='CM/S*CM**-3             '
      TXTUNT(1,73)='CM/S*CM**-3             '
      TXTUNT(1,74)='CM/S*CM**-3             '
      TXTUNT(1,75)='AMP*CM**-3              '
      TXTUNT(1,76)='AMP*CM**-3              '
      TXTUNT(1,77)='AMP*CM**-3              '
      TXTUNT(1,78)='AMP*CM**-3              '
      TXTUNT(1,79)='AMP*CM**-3              '
      TXTUNT(1,80)='WATT*CM**-3             '
      TXTUNT(1,81)='WATT*CM**-3             '
      TXTUNT(1,82)='WATT*CM**-3             '
      TXTUNT(1,83)='WATT*CM**-3             '
      TXTUNT(1,84)='WATT*CM**-3             '
      TXTUNT(1,85)='G*CM/SEC*CM**-3         '
      TXTUNT(1,86)='G*CM/SEC*CM**-3         '
      TXTUNT(1,87)='G*CM/SEC*CM**-3         '
      TXTUNT(1,88)='G*CM/SEC*CM**-3         '
      TXTUNT(1,89)='G*CM/SEC*CM**-3         '
      TXTUNT(1,90)='G*CM/SEC*CM**-3         '
      TXTUNT(1,91)='G*CM/SEC*CM**-3         '
      TXTUNT(1,92)='G*CM/SEC*CM**-3         '
      TXTUNT(1,93)='G*CM/SEC*CM**-3         '
      TXTUNT(1,94)='G*CM/SEC*CM**-3         '
      TXTUNT(1,95)='G*CM/SEC*CM**-3         '
      TXTUNT(1,96)='G*CM/SEC*CM**-3         '
      TXTUNT(1,97)='AMP*CM**-3              '
      TXTUNT(1,98)='AMP*CM**-3              '
      TXTUNT(1,99)='AMP*CM**-3              '
      TXTUNT(1,100)='AMP*CM**-3             '
      DO 2 J=1,NTALV
        DO 2 I=2,N1MX
          TEXT24=TXTUNT(1,J)
          TXTUNT(I,J)=TEXT24
2     CONTINUE
C
      TXTTLW(1,1)='PARTICLE FLUX, INCIDENT, ATOMS                   '
      TXTTLW(1,2)='PARTICLE FLUX, EMITTED, ATS. => ATOMS            '
      TXTTLW(1,3)='PARTICLE FLUX, EMITTED, MLS. => ATOMS            '
      TXTTLW(1,4)='PARTICLE FLUX, EMITTED, T.I. => ATOMS            '
      TXTTLW(1,5)='PARTICLE FLUX, EMITTED, PHS. => ATOMS            '
      TXTTLW(1,6)='PARTICLE FLUX, EMITTED, B.I. => ATOMS            '
      TXTTLW(1,7)='PARTICLE FLUX, INCIDENT, MOLECULES               '
      TXTTLW(1,8)='PARTICLE FLUX, EMITTED, ATS. => MOLECULES        '
      TXTTLW(1,9)='PARTICLE FLUX, EMITTED, MLS. => MOLECULES        '
      TXTTLW(1,10)='PARTICLE FLUX, EMITTED, T.I. => MOLECULES        '
      TXTTLW(1,11)='PARTICLE FLUX, EMITTED, PHS. => MOLECULES        '
      TXTTLW(1,12)='PARTICLE FLUX, EMITTED, B.I. => MOLECULES        '
      TXTTLW(1,13)='PARTICLE FLUX, INCIDENT, TEST IONS               '
      TXTTLW(1,14)='PARTICLE FLUX, EMITTED, ATS. => TEST IONS        '
      TXTTLW(1,15)='PARTICLE FLUX, EMITTED, MLS. => TEST IONS        '
      TXTTLW(1,16)='PARTICLE FLUX, EMITTED, T.I. => TEST IONS        '
      TXTTLW(1,17)='PARTICLE FLUX, EMITTED, PHS. => TEST IONS        '
      TXTTLW(1,18)='PARTICLE FLUX, EMITTED, B.I. => TEST IONS        '
      TXTTLW(1,19)='PARTICLE FLUX, INCIDENT, PHOTONS                 '
      TXTTLW(1,20)='PARTICLE FLUX, EMITTED, ATS. => PHOTONS          '
      TXTTLW(1,21)='PARTICLE FLUX, EMITTED, MLS. => PHOTONS          '
      TXTTLW(1,22)='PARTICLE FLUX, EMITTED, T.I. => PHOTONS          '
      TXTTLW(1,23)='PARTICLE FLUX, EMITTED, PHS. => PHOTONS          '
      TXTTLW(1,24)='PARTICLE FLUX, EMITTED, B.I. => PHOTONS          '
      TXTTLW(1,25)='PARTICLE FLUX, INCIDENT, BULK IONS               '
 
      TXTTLW(1,26)='ENERGY FLUX, INCIDENT, ATOMS                     '
      TXTTLW(1,27)='ENERGY FLUX, EMITTED, ATS. => ATOMS              '
      TXTTLW(1,28)='ENERGY FLUX, EMITTED, MLS. => ATOMS              '
      TXTTLW(1,29)='ENERGY FLUX, EMITTED, T.I. => ATOMS              '
      TXTTLW(1,30)='ENERGY FLUX, EMITTED, PHS. => ATOMS              '
      TXTTLW(1,31)='ENERGY FLUX, EMITTED, B.I. => ATOMS              '
      TXTTLW(1,32)='ENERGY FLUX, INCIDENT, MOLECULES                 '
      TXTTLW(1,33)='ENERGY FLUX, EMITTED, ATS. => MOLECULES          '
      TXTTLW(1,34)='ENERGY FLUX, EMITTED, MLS. => MOLECULES          '
      TXTTLW(1,35)='ENERGY FLUX, EMITTED, T.I. => MOLECULES          '
      TXTTLW(1,36)='ENERGY FLUX, EMITTED, PHS. => MOLECULES          '
      TXTTLW(1,37)='ENERGY FLUX, EMITTED, B.I. => MOLECULES          '
      TXTTLW(1,38)='ENERGY FLUX, INCIDENT, TEST IONS                 '
      TXTTLW(1,39)='ENERGY FLUX, EMITTED, ATS. => TEST IONS          '
      TXTTLW(1,40)='ENERGY FLUX, EMITTED, MLS. => TEST IONS          '
      TXTTLW(1,41)='ENERGY FLUX, EMITTED, T.I. => TEST IONS          '
      TXTTLW(1,42)='ENERGY FLUX, EMITTED, PHS. => TEST IONS          '
      TXTTLW(1,43)='ENERGY FLUX, EMITTED, B.I. => TEST IONS          '
      TXTTLW(1,44)='ENERGY FLUX, INCIDENT, PHOTONS                   '
      TXTTLW(1,45)='ENERGY FLUX, EMITTED, ATS. => PHOTONS            '
      TXTTLW(1,46)='ENERGY FLUX, EMITTED, MLS. => PHOTONS            '
      TXTTLW(1,47)='ENERGY FLUX, EMITTED, T.I. => PHOTONS            '
      TXTTLW(1,48)='ENERGY FLUX, EMITTED, PHS. => PHOTONS            '
      TXTTLW(1,49)='ENERGY FLUX, EMITTED, B.I. => PHOTONS            '
      TXTTLW(1,50)='ENERGY FLUX, INCIDENT, BULK IONS                 '
 
      TXTTLW(1,51)='SPUTTERED FLUX BY INCIDENT ATOMS                 '
      TXTTLW(1,52)='SPUTTERED FLUX BY INCIDENT MOLECULES             '
      TXTTLW(1,53)='SPUTTERED FLUX BY INCIDENT TEST IONS             '
      TXTTLW(1,54)='SPUTTERED FLUX BY INCIDENT PHOTONS               '
      TXTTLW(1,55)='SPUTTERED FLUX BY INCIDENT BULK IONS             '
      TXTTLW(1,56)='SPUTTERED FLUX, TOTAL                            '
C  TALLY NTLSA=57 (SEE PARMMOD.F)
C   ADDIT. TALLIES, SUBR. UPSUSR.F
C   TXTTLW IS OVERWRITTEN BY INPUT BLOCK 10D
      TXTTLW(1,57)='ADDITIONAL SURFACE TALLY, SUBR. UPSUSR.F         '
C  TALLY NTLSR=58 (SEE PARMMOD.F)
C   ADDIT. TALLIES, ALGEBRAIC EXPRESSION IN EXISTING TALLIES
C   TXTTLW IS OVERWRITTEN BY INPUT BLOCK 10E
      TXTTLW(1,58)='ALGEBRAIC EXPRESSION IN SURFACE AVERAGED TALLIES '
      TXTTLW(1,59)='PUMPED FLUX BY SPECIES                           '
C
      DO J=1,NTALS
        TXTTLW(2:N2MX,J)=TXTTLW(1,J)
      END DO
C
      TXTUNW(1,1)='AMP                     '
      TXTUNW(1,2)='AMP                     '
      TXTUNW(1,3)='AMP                     '
      TXTUNW(1,4)='AMP                     '
      TXTUNW(1,5)='AMP                     '
      TXTUNW(1,6)='AMP                     '
      TXTUNW(1,7)='AMP                     '
      TXTUNW(1,8)='AMP                     '
      TXTUNW(1,9)='AMP                     '
      TXTUNW(1,10)='AMP                     '
      TXTUNW(1,11)='AMP                     '
      TXTUNW(1,12)='AMP                     '
      TXTUNW(1,13)='AMP                     '
      TXTUNW(1,14)='AMP                     '
      TXTUNW(1,15)='AMP                     '
      TXTUNW(1,16)='AMP                     '
      TXTUNW(1,17)='AMP                     '
      TXTUNW(1,18)='AMP                     '
      TXTUNW(1,19)='AMP                     '
      TXTUNW(1,20)='AMP                     '
      TXTUNW(1,21)='AMP                     '
      TXTUNW(1,22)='AMP                     '
      TXTUNW(1,23)='AMP                     '
      TXTUNW(1,24)='AMP                     '
      TXTUNW(1,25)='AMP                     '
 
      TXTUNW(1,26)='WATT                    '
      TXTUNW(1,27)='WATT                    '
      TXTUNW(1,28)='WATT                    '
      TXTUNW(1,29)='WATT                    '
      TXTUNW(1,30)='WATT                    '
      TXTUNW(1,31)='WATT                    '
      TXTUNW(1,32)='WATT                    '
      TXTUNW(1,33)='WATT                    '
      TXTUNW(1,34)='WATT                    '
      TXTUNW(1,35)='WATT                    '
      TXTUNW(1,36)='WATT                    '
      TXTUNW(1,37)='WATT                    '
      TXTUNW(1,38)='WATT                    '
      TXTUNW(1,39)='WATT                    '
      TXTUNW(1,40)='WATT                    '
      TXTUNW(1,41)='WATT                    '
      TXTUNW(1,42)='WATT                    '
      TXTUNW(1,43)='WATT                    '
      TXTUNW(1,44)='WATT                    '
      TXTUNW(1,45)='WATT                    '
      TXTUNW(1,46)='WATT                    '
      TXTUNW(1,47)='WATT                    '
      TXTUNW(1,48)='WATT                    '
      TXTUNW(1,49)='WATT                    '
      TXTUNW(1,50)='WATT                    '
 
      TXTUNW(1,51)='AMP                     '
      TXTUNW(1,52)='AMP                     '
      TXTUNW(1,53)='AMP                     '
      TXTUNW(1,54)='AMP                     '
      TXTUNW(1,55)='AMP                     '
      TXTUNW(1,56)='AMP                     '
      TXTUNW(1,57)='TO BE READ              '
      TXTUNW(1,58)='TO BE READ              '
      TXTUNW(1,59)='AMP                     '
      DO J=1,NTALS
        TXTUNW(2:N2MX,J)=TXTUNW(1,J)
      END DO
C
      TXTPLS(1,1)='PLASMA TEMPERATURE                               '
      TXTPLS(1,2)='PLASMA TEMPERATURE                               '
      TXTPLS(1,3)='PLASMA DENSITY (BULK PARTICLES)                  '
      TXTPLS(1,4)='PLASMA DENSITY (BULK PARTICLES)                  '
      TXTPLS(1,5)='DRIFT VELOCITY IN X-DIRECTION (BULK IONS)        '
      TXTPLS(1,6)='DRIFT VELOCITY IN Y-DIRECTION (BULK IONS)        '
      TXTPLS(1,7)='DRIFT VELOCITY IN Z-DIRECTION (BULK IONS)        '
      TXTPLS(1,8)='MAGN. FIELD UNIT VECTOR, X DIRECTION             '
      TXTPLS(1,9)='MAGN. FIELD UNIT VECTOR, Y DIRECTION             '
      TXTPLS(1,10)='MAGN. FIELD UNIT VECTOR, Z DIRECTION             '
      TXTPLS(1,11)='MAGN. FIELD STRENGTH                             '
C     TXTPLS(1,12)='TO BE READ                                       '
      TXTPLS(1,13)='BULK ION KINETIC DRIFT ENERGY                    '
      TXTPLS(1,14)='ZONE VOLUMES                                     '
     
      TXTPLS(1,15)='SPACE-SPECIES WEIGHT FUNCTION                    '
      TXTPLS(1,16)='PERP. MAGN. FIELD VECTOR, X DIRECTION            '
      TXTPLS(1,17)='PERP. MAGN. FIELD VECTOR, Y DIRECTION            '
C
      DO J=1,NTALI
        IF (J.NE.12) THEN
          DO I=2,N1MX
            TEXT72=TXTPLS(1,J)
            TXTPLS(I,J)=TEXT72
          ENDDO
        ENDIF
      ENDDO
C
      TXTPUN(1,1)='EV                      '
      TXTPUN(1,2)='EV                      '
      TXTPUN(1,3)='CM**-3                  '
      TXTPUN(1,4)='CM**-3                  '
      TXTPUN(1,5)='CM/SEC                  '
      TXTPUN(1,6)='CM/SEC                  '
      TXTPUN(1,7)='CM/SEC                  '
      TXTPUN(1,8)=' ---                    '
      TXTPUN(1,9)=' ---                    '
      TXTPUN(1,10)=' ---                    '
      TXTPUN(1,11)='TESLA                   '
C     TXTPUN(1,12)='TO BE READ              '
      TXTPUN(1,13)='EV                      '
      TXTPUN(1,14)='CM**3                   '
      TXTPUN(1,15)=' ---                    '
      TXTPUN(1,16)=' ---                    '
      TXTPUN(1,17)=' ---                    '
C
      DO J=1,NTALI
        IF (J.NE.12) THEN
          DO I=2,N1MX
            TEXT24=TXTPUN(1,J)
            TXTPUN(I,J)=TEXT24
          ENDDO
        ENDIF
      ENDDO
      RETURN
C
      ENTRY EIRENE_STTXT1
C
C
      NFSTVI(1)=NATMI
      NFSTVI(2)=NMOLI
      NFSTVI(3)=NIONI
      NFSTVI(4)=NPHOTI
      NFSTVI(5)=NATMI
      NFSTVI(6)=NMOLI
      NFSTVI(7)=NIONI
      NFSTVI(8)=NPHOTI
      NFSTVI(9)=1
      NFSTVI(10)=NATMI
      NFSTVI(11)=NMOLI
      NFSTVI(12)=NIONI
      NFSTVI(13)=NPHOTI
      NFSTVI(14)=NPLSI
      NFSTVI(15)=1
      NFSTVI(16)=NATMI
      NFSTVI(17)=NMOLI
      NFSTVI(18)=NIONI
      NFSTVI(19)=NPHOTI
      NFSTVI(20)=NPLSI
      NFSTVI(21)=1
      NFSTVI(22)=NATMI
      NFSTVI(23)=NMOLI
      NFSTVI(24)=NIONI
      NFSTVI(25)=NPHOTI
      NFSTVI(26)=NPLSI
      NFSTVI(27)=1
      NFSTVI(28)=NATMI
      NFSTVI(29)=NMOLI
      NFSTVI(30)=NIONI
      NFSTVI(31)=NPHOTI
      NFSTVI(32)=NPLSI
      NFSTVI(33)=1
      NFSTVI(34)=1
      NFSTVI(35)=1
      NFSTVI(36)=1
      NFSTVI(37)=1
      NFSTVI(38)=1
      NFSTVI(39)=1
      NFSTVI(40)=1
      NFSTVI(41)=1
      NFSTVI(42)=1
      NFSTVI(43)=1
      NFSTVI(44)=1
      NFSTVI(45)=1
      NFSTVI(46)=1
      NFSTVI(47)=1
      NFSTVI(48)=1
      NFSTVI(49)=1
      NFSTVI(50)=1
      NFSTVI(51)=1
      NFSTVI(52)=1
      NFSTVI(53)=1
      NFSTVI(54)=1
      NFSTVI(55)=1
      NFSTVI(56)=1
      NFSTVI(NTALA)=NADVI
      NFSTVI(NTALC)=NCLVI
      NFSTVI(NTALT)=NSNVI
      NFSTVI(NTALM)=NCPVI
C     NFSTVI(NTALB) IS DEFINED IN SUBR. XSECT...
      NFSTVI(NTALB)=0
      NFSTVI(NTALR)=NALVI
      NFSTVI(75)=NATMI
      NFSTVI(76)=NMOLI
      NFSTVI(77)=NIONI
      NFSTVI(78)=NPHOTI
      NFSTVI(79)=NPLSI
      NFSTVI(80)=1
      NFSTVI(81)=1
      NFSTVI(82)=1
      NFSTVI(83)=1
      NFSTVI(84)=1
      NFSTVI(85)=NATMI
      NFSTVI(86)=NMOLI
      NFSTVI(87)=NIONI
      NFSTVI(88)=NPHOTI
      NFSTVI(89)=NATMI
      NFSTVI(90)=NMOLI
      NFSTVI(91)=NIONI
      NFSTVI(92)=NPHOTI
      NFSTVI(93)=NATMI
      NFSTVI(94)=NMOLI
      NFSTVI(95)=NIONI
      NFSTVI(96)=NPHOTI
      NFSTVI(97)=NPLSI
      NFSTVI(98)=NPLSI
      NFSTVI(99)=NPLSI
      NFSTVI(100)=NPLSI
 
C
C
      NFSTWI(1)=NATMI
      NFSTWI(2)=NATMI
      NFSTWI(3)=NATMI
      NFSTWI(4)=NATMI
      NFSTWI(5)=NATMI
      NFSTWI(6)=NATMI
      NFSTWI(7)=NMOLI
      NFSTWI(8)=NMOLI
      NFSTWI(9)=NMOLI
      NFSTWI(10)=NMOLI
      NFSTWI(11)=NMOLI
      NFSTWI(12)=NMOLI
      NFSTWI(13)=NIONI
      NFSTWI(14)=NIONI
      NFSTWI(15)=NIONI
      NFSTWI(16)=NIONI
      NFSTWI(17)=NIONI
      NFSTWI(18)=NIONI
      NFSTWI(19)=NPHOTI
      NFSTWI(20)=NPHOTI
      NFSTWI(21)=NPHOTI
      NFSTWI(22)=NPHOTI
      NFSTWI(23)=NPHOTI
      NFSTWI(24)=NPHOTI
      NFSTWI(25)=NPLSI
 
      NFSTWI(26)=NATMI
      NFSTWI(27)=NATMI
      NFSTWI(28)=NATMI
      NFSTWI(29)=NATMI
      NFSTWI(30)=NATMI
      NFSTWI(31)=NATMI
      NFSTWI(32)=NMOLI
      NFSTWI(33)=NMOLI
      NFSTWI(34)=NMOLI
      NFSTWI(35)=NMOLI
      NFSTWI(36)=NMOLI
      NFSTWI(37)=NMOLI
      NFSTWI(38)=NIONI
      NFSTWI(39)=NIONI
      NFSTWI(40)=NIONI
      NFSTWI(41)=NIONI
      NFSTWI(42)=NIONI
      NFSTWI(43)=NIONI
      NFSTWI(44)=NPHOTI
      NFSTWI(45)=NPHOTI
      NFSTWI(46)=NPHOTI
      NFSTWI(47)=NPHOTI
      NFSTWI(48)=NPHOTI
      NFSTWI(49)=NPHOTI
      NFSTWI(50)=NPLSI
 
      NFSTWI(51)=NATMI
      NFSTWI(52)=NMOLI
      NFSTWI(53)=NIONI
      NFSTWI(54)=NPHOTI
      NFSTWI(55)=NPLSI
      NFSTWI(56)=1
      NFSTWI(NTLSA)=NADSI
      NFSTWI(NTLSR)=NALSI
      NFSTWI(NTALS)=NSPTOT
C
C
      NFSTPI(1)=1
      NFSTPI(2)=NPLSI
      NFSTPI(3)=1
      NFSTPI(4)=NPLSI
      NFSTPI(5)=NPLSI
      NFSTPI(6)=NPLSI
      NFSTPI(7)=NPLSI
      NFSTPI(8)=1
      NFSTPI(9)=1
      NFSTPI(10)=1
      NFSTPI(11)=1
      NFSTPI(12)=NAINI
      NFSTPI(13)=NPLSI
      NFSTPI(14)=1
      NFSTPI(15)=NATMI+NMOLI+NIONI
      NFSTPI(16)=1
      NFSTPI(17)=1
C
C  INITIALISE SPECIES ARRAYS FOR VOLUME TALLIES
 
      N1=NPHOTI
      N2=N1+NATMI
      N3=N2+NMOLI
      N4=N3+NIONI
      N5=N4+NPLSI
      N6=N5+NADVI
      N7=N6+NALVI
      N8=N7+NCLVI
      N9=N8+NCPVI
      N10=N9+NBGVI
      N11=N10+NSNVI
 
      NSPAN(1)=N1+1
      NSPAN(2)=N2+1
      NSPAN(3)=N3+1
      NSPAN(4)=1
      NSPAN(5)=N1+1
      NSPAN(6)=N2+1
      NSPAN(7)=N3+1
      NSPAN(8)=1
      NSPAN(9)=0
      NSPAN(10)=N1+1
      NSPAN(11)=N2+1
      NSPAN(12)=N3+1
      NSPAN(13)=1
      NSPAN(14)=N4+1
      NSPAN(15)=0
      NSPAN(16)=N1+1
      NSPAN(17)=N2+1
      NSPAN(18)=N3+1
      NSPAN(19)=1
      NSPAN(20)=N4+1
      NSPAN(21)=0
      NSPAN(22)=N1+1
      NSPAN(23)=N2+1
      NSPAN(24)=N2+1
      NSPAN(25)=1
      NSPAN(26)=N4+1
      NSPAN(27)=0
      NSPAN(28)=N1+1
      NSPAN(29)=N2+1
      NSPAN(30)=N2+1
      NSPAN(31)=1
      NSPAN(32)=N4+1
      NSPAN(33)=0
      NSPAN(34)=0
      NSPAN(35)=0
      NSPAN(36)=0
      NSPAN(37)=0
      NSPAN(38)=0
      NSPAN(39)=0
      NSPAN(40)=0
      NSPAN(41)=0
      NSPAN(42)=0
      NSPAN(43)=0
      NSPAN(44)=0
      NSPAN(45)=0
      NSPAN(46)=0
      NSPAN(47)=0
      NSPAN(48)=0
      NSPAN(49)=0
      NSPAN(50)=0
      NSPAN(51)=0
      NSPAN(52)=0
      NSPAN(53)=0
      NSPAN(54)=0
      NSPAN(55)=0
      NSPAN(56)=0
      NSPAN(NTALA)=N5+1
      NSPAN(NTALC)=N7+1
      NSPAN(NTALT)=N10+1
      NSPAN(NTALM)=N8+1
      NSPAN(NTALB)=N9+1
      NSPAN(NTALR)=N6+1
C  GENERATION LIMIT TALLIES
      NSPAN(63)=N1+1
      NSPAN(64)=N2+1
      NSPAN(65)=N3+1
      NSPAN(66)=1
      NSPAN(67)=N1+1
      NSPAN(68)=N2+1
      NSPAN(69)=N3+1
      NSPAN(70)=1
      NSPAN(71)=N1+1
      NSPAN(72)=N2+1
      NSPAN(73)=N3+1
      NSPAN(74)=1
      NSPAN(75)=N1+1
      NSPAN(76)=N2+1
      NSPAN(77)=N3+1
      NSPAN(78)=1
      NSPAN(79)=N4+1
      NSPAN(80)=0
      NSPAN(81)=0
      NSPAN(82)=0
      NSPAN(83)=0
      NSPAN(84)=0
      NSPAN(85)=N1+1
      NSPAN(86)=N2+1
      NSPAN(87)=N3+1
      NSPAN(88)=1
      NSPAN(89)=N1+1
      NSPAN(90)=N2+1
      NSPAN(91)=N3+1
      NSPAN(92)=1
      NSPAN(93)=N1+1
      NSPAN(94)=N2+1
      NSPAN(95)=N3+1
      NSPAN(96)=1
      NSPAN(97)=N4+1
      NSPAN(98)=N4+1
      NSPAN(99)=N4+1
      NSPAN(100)=N4+1
 
      NSPEN(1)=N2
      NSPEN(2)=N3
      NSPEN(3)=N4
      NSPEN(4)=N1
      NSPEN(5)=N2
      NSPEN(6)=N3
      NSPEN(7)=N4
      NSPEN(8)=N1
      NSPEN(9)=0
      NSPEN(10)=N2
      NSPEN(11)=N3
      NSPEN(12)=N4
      NSPEN(13)=N1
      NSPEN(14)=N5
      NSPEN(15)=0
      NSPEN(16)=N2
      NSPEN(17)=N3
      NSPEN(18)=N4
      NSPEN(19)=N1
      NSPEN(20)=N5
      NSPEN(21)=0
      NSPEN(22)=N2
      NSPEN(23)=N3
      NSPEN(24)=N4
      NSPEN(25)=N1
      NSPEN(26)=N5
      NSPEN(27)=0
      NSPEN(28)=N2
      NSPEN(29)=N3
      NSPEN(30)=N4
      NSPEN(31)=N1
      NSPEN(32)=N5
      NSPEN(33)=0
      NSPEN(34)=0
      NSPEN(35)=0
      NSPEN(36)=0
      NSPEN(37)=0
      NSPEN(38)=0
      NSPEN(39)=0
      NSPEN(40)=0
      NSPEN(41)=0
      NSPEN(42)=0
      NSPEN(43)=0
      NSPEN(44)=0
      NSPEN(45)=0
      NSPEN(46)=0
      NSPEN(47)=0
      NSPEN(48)=0
      NSPEN(49)=0
      NSPEN(50)=0
      NSPEN(51)=0
      NSPEN(52)=0
      NSPEN(53)=0
      NSPEN(54)=0
      NSPEN(55)=0
      NSPEN(56)=0
      NSPEN(NTALA)=N6
      NSPEN(NTALC)=N8
      NSPEN(NTALT)=N11
      NSPEN(NTALM)=N9
      NSPEN(NTALB)=N10
      NSPEN(NTALR)=N7
C  GENERATION LIMIT TALLIES
      NSPEN(63)=N2
      NSPEN(64)=N3
      NSPEN(65)=N4
      NSPEN(66)=N1
      NSPEN(67)=N2
      NSPEN(68)=N3
      NSPEN(69)=N4
      NSPEN(70)=N1
      NSPEN(71)=N2
      NSPEN(72)=N3
      NSPEN(73)=N4
      NSPEN(74)=N1
      NSPEN(75)=N2
      NSPEN(76)=N3
      NSPEN(77)=N4
      NSPEN(78)=N1
      NSPEN(79)=N5
      NSPEN(80)=0
      NSPEN(81)=0
      NSPEN(82)=0
      NSPEN(83)=0
      NSPEN(84)=0
      NSPEN(85)=N2
      NSPEN(86)=N3
      NSPEN(87)=N4
      NSPEN(88)=N1
      NSPEN(89)=N2
      NSPEN(90)=N3
      NSPEN(91)=N4
      NSPEN(92)=N1
      NSPEN(93)=N2
      NSPEN(94)=N3
      NSPEN(95)=N4
      NSPEN(96)=N1
      NSPEN(97)=N5
      NSPEN(98)=N5
      NSPEN(99)=N5
      NSPEN(100)=N5
 
      DO IPHOT=1,NPHOTI
        ISPZ=IPHOT
        TXTSPC(IPHOT,4)=TEXTS(ISPZ)
        TXTSPC(IPHOT,8)=TEXTS(ISPZ)
        TXTSPC(IPHOT,13)=TEXTS(ISPZ)
        TXTSPC(IPHOT,19)=TEXTS(ISPZ)
        TXTSPC(IPHOT,25)=TEXTS(ISPZ)
        TXTSPC(IPHOT,31)=TEXTS(ISPZ)
        TXTSPC(IPHOT,66)=TEXTS(ISPZ)
        TXTSPC(IPHOT,70)=TEXTS(ISPZ)
        TXTSPC(IPHOT,74)=TEXTS(ISPZ)
        TXTSPC(IPHOT,78)=TEXTS(ISPZ)
        TXTSPC(IPHOT,88)=TEXTS(ISPZ)
        TXTSPC(IPHOT,92)=TEXTS(ISPZ)
        TXTSPC(IPHOT,96)=TEXTS(ISPZ)
        TXTSPC(IPHOT,100)=TEXTS(ISPZ)
      END DO
 
      DO 10 IATM=1,NATMI
        ISPZ=NSPH+IATM
        TXTSPC(IATM,1)=TEXTS(ISPZ)
        TXTSPC(IATM,5)=TEXTS(ISPZ)
        TXTSPC(IATM,10)=TEXTS(ISPZ)
        TXTSPC(IATM,16)=TEXTS(ISPZ)
        TXTSPC(IATM,22)=TEXTS(ISPZ)
        TXTSPC(IATM,28)=TEXTS(ISPZ)
        TXTSPC(IATM,63)=TEXTS(ISPZ)
        TXTSPC(IATM,67)=TEXTS(ISPZ)
        TXTSPC(IATM,71)=TEXTS(ISPZ)
        TXTSPC(IATM,75)=TEXTS(ISPZ)
        TXTSPC(IATM,85)=TEXTS(ISPZ)
        TXTSPC(IATM,89)=TEXTS(ISPZ)
        TXTSPC(IATM,93)=TEXTS(ISPZ)
        TXTSPC(IATM,97)=TEXTS(ISPZ)
10    CONTINUE
C
      DO 20 IMOL=1,NMOLI
        ISPZ=NSPA+IMOL
        TXTSPC(IMOL,2)=TEXTS(ISPZ)
        TXTSPC(IMOL,6)=TEXTS(ISPZ)
        TXTSPC(IMOL,11)=TEXTS(ISPZ)
        TXTSPC(IMOL,17)=TEXTS(ISPZ)
        TXTSPC(IMOL,23)=TEXTS(ISPZ)
        TXTSPC(IMOL,29)=TEXTS(ISPZ)
        TXTSPC(IMOL,64)=TEXTS(ISPZ)
        TXTSPC(IMOL,68)=TEXTS(ISPZ)
        TXTSPC(IMOL,72)=TEXTS(ISPZ)
        TXTSPC(IMOL,76)=TEXTS(ISPZ)
        TXTSPC(IMOL,86)=TEXTS(ISPZ)
        TXTSPC(IMOL,90)=TEXTS(ISPZ)
        TXTSPC(IMOL,94)=TEXTS(ISPZ)
        TXTSPC(IMOL,98)=TEXTS(ISPZ)
20    CONTINUE
C
      DO 30 IION=1,NIONI
        ISPZ=NSPAM+IION
        TXTSPC(IION,3)=TEXTS(ISPZ)
        TXTSPC(IION,7)=TEXTS(ISPZ)
        TXTSPC(IION,12)=TEXTS(ISPZ)
        TXTSPC(IION,18)=TEXTS(ISPZ)
        TXTSPC(IION,24)=TEXTS(ISPZ)
        TXTSPC(IION,30)=TEXTS(ISPZ)
        TXTSPC(IION,65)=TEXTS(ISPZ)
        TXTSPC(IION,69)=TEXTS(ISPZ)
        TXTSPC(IION,73)=TEXTS(ISPZ)
        TXTSPC(IION,77)=TEXTS(ISPZ)
        TXTSPC(IION,87)=TEXTS(ISPZ)
        TXTSPC(IION,91)=TEXTS(ISPZ)
        TXTSPC(IION,95)=TEXTS(ISPZ)
        TXTSPC(IION,99)=TEXTS(ISPZ)
30    CONTINUE
C
      DO 40 IPLS=1,NPLSI
        ISPZ=NSPAMI+IPLS
        TXTSPC(IPLS,14)=TEXTS(ISPZ)
        TXTSPC(IPLS,20)=TEXTS(ISPZ)
        TXTSPC(IPLS,26)=TEXTS(ISPZ)
        TXTSPC(IPLS,32)=TEXTS(ISPZ)
        TXTSPC(IPLS,79)=TEXTS(ISPZ)
40    CONTINUE
C
      TXTSPC(1,9)='ELECTRONS               '
      TXTSPC(1,15)='ELECTRONS               '
      TXTSPC(1,21)='ELECTRONS               '
      TXTSPC(1,27)='ELECTRONS               '
      TXTSPC(1,33)='ELECTRONS               '
      TXTSPC(1,39)='ELECTRONS               '
      TXTSPC(1,45)='ELECTRONS               '
      TXTSPC(1,51)='ELECTRONS               '
C
      TXTSPC(1,34)='ATOMS                   '
      TXTSPC(1,40)='ATOMS                   '
      TXTSPC(1,46)='ATOMS                   '
      TXTSPC(1,52)='ATOMS                   '
      TXTSPC(1,80)='ATOMS                   '
C
      TXTSPC(1,35)='MOLECULES               '
      TXTSPC(1,41)='MOLECULES               '
      TXTSPC(1,47)='MOLECULES               '
      TXTSPC(1,53)='MOLECULES               '
      TXTSPC(1,81)='MOLECULES               '
C
      TXTSPC(1,36)='TEST IONS               '
      TXTSPC(1,42)='TEST IONS               '
      TXTSPC(1,48)='TEST IONS               '
      TXTSPC(1,54)='TEST IONS               '
      TXTSPC(1,82)='TEST IONS               '
C
      TXTSPC(1,37)='PHOTONS                 '
      TXTSPC(1,43)='PHOTONS                 '
      TXTSPC(1,49)='PHOTONS                 '
      TXTSPC(1,55)='PHOTONS                 '
      TXTSPC(1,83)='PHOTONS                 '
C
      TXTSPC(1,38)='BULK IONS               '
      TXTSPC(1,44)='BULK IONS               '
      TXTSPC(1,50)='BULK IONS               '
      TXTSPC(1,56)='BULK IONS               '
      TXTSPC(1,84)='BULK IONS               '
C
 
      DO IPHOT=1,NPHOTI
        ISPZ=IPHOT
        TXTSPW(IPHOT,19)=TEXTS(ISPZ)
        TXTSPW(IPHOT,20)=TEXTS(ISPZ)
        TXTSPW(IPHOT,21)=TEXTS(ISPZ)
        TXTSPW(IPHOT,22)=TEXTS(ISPZ)
        TXTSPW(IPHOT,23)=TEXTS(ISPZ)
        TXTSPW(IPHOT,24)=TEXTS(ISPZ)
        TXTSPW(IPHOT,44)=TEXTS(ISPZ)
        TXTSPW(IPHOT,45)=TEXTS(ISPZ)
        TXTSPW(IPHOT,46)=TEXTS(ISPZ)
        TXTSPW(IPHOT,47)=TEXTS(ISPZ)
        TXTSPW(IPHOT,48)=TEXTS(ISPZ)
        TXTSPW(IPHOT,49)=TEXTS(ISPZ)
        TXTSPW(IPHOT,54)=TEXTS(ISPZ)
      END DO
 
      DO IATM=1,NATMI
        ISPZ=NSPH+IATM
        TXTSPW(IATM,1)=TEXTS(ISPZ)
        TXTSPW(IATM,2)=TEXTS(ISPZ)
        TXTSPW(IATM,3)=TEXTS(ISPZ)
        TXTSPW(IATM,4)=TEXTS(ISPZ)
        TXTSPW(IATM,5)=TEXTS(ISPZ)
        TXTSPW(IATM,6)=TEXTS(ISPZ)
        TXTSPW(IATM,26)=TEXTS(ISPZ)
        TXTSPW(IATM,27)=TEXTS(ISPZ)
        TXTSPW(IATM,28)=TEXTS(ISPZ)
        TXTSPW(IATM,29)=TEXTS(ISPZ)
        TXTSPW(IATM,30)=TEXTS(ISPZ)
        TXTSPW(IATM,31)=TEXTS(ISPZ)
        TXTSPW(IATM,51)=TEXTS(ISPZ)
      END DO
C
      DO IMOL=1,NMOLI
        ISPZ=NSPA+IMOL
        TXTSPW(IMOL,7)=TEXTS(ISPZ)
        TXTSPW(IMOL,8)=TEXTS(ISPZ)
        TXTSPW(IMOL,9)=TEXTS(ISPZ)
        TXTSPW(IMOL,10)=TEXTS(ISPZ)
        TXTSPW(IMOL,11)=TEXTS(ISPZ)
        TXTSPW(IMOL,12)=TEXTS(ISPZ)
        TXTSPW(IMOL,32)=TEXTS(ISPZ)
        TXTSPW(IMOL,33)=TEXTS(ISPZ)
        TXTSPW(IMOL,34)=TEXTS(ISPZ)
        TXTSPW(IMOL,35)=TEXTS(ISPZ)
        TXTSPW(IMOL,36)=TEXTS(ISPZ)
        TXTSPW(IMOL,37)=TEXTS(ISPZ)
        TXTSPW(IMOL,52)=TEXTS(ISPZ)
      END DO
C
      DO IION=1,NIONI
        ISPZ=NSPAM+IION
        TXTSPW(IION,13)=TEXTS(ISPZ)
        TXTSPW(IION,14)=TEXTS(ISPZ)
        TXTSPW(IION,15)=TEXTS(ISPZ)
        TXTSPW(IION,16)=TEXTS(ISPZ)
        TXTSPW(IION,17)=TEXTS(ISPZ)
        TXTSPW(IION,18)=TEXTS(ISPZ)
        TXTSPW(IION,38)=TEXTS(ISPZ)
        TXTSPW(IION,39)=TEXTS(ISPZ)
        TXTSPW(IION,40)=TEXTS(ISPZ)
        TXTSPW(IION,41)=TEXTS(ISPZ)
        TXTSPW(IION,42)=TEXTS(ISPZ)
        TXTSPW(IION,43)=TEXTS(ISPZ)
        TXTSPW(IION,53)=TEXTS(ISPZ)
      END DO
C
      DO IPLS=1,NPLSI
        ISPZ=NSPAMI+IPLS
        TXTSPW(IPLS,25)=TEXTS(ISPZ)
        TXTSPW(IPLS,50)=TEXTS(ISPZ)
        TXTSPW(IPLS,55)=TEXTS(ISPZ)
      END DO
C
      TXTSPW(1,42)='                        '
      TXTSPW(1,43)='                        '
      TXTSPW(1,44)='                        '
      TXTSPW(1,45)='                        '
C
      TXTPSP(1,1)='ELECTRONS               '
      TXTPSP(1,3)='ELECTRONS               '
      TXTPSP(1,8)=' ---                    '
      TXTPSP(1,9)=' ---                    '
      TXTPSP(1,10)=' ---                    '
      TXTPSP(1,11)=' ---                    '
      TXTPSP(1,14)=' ---                    '
      TXTPSP(1,16)=' ---                    '
      TXTPSP(1,17)=' ---                    '
C
C     TXTPSP(IAIN,12)='TO BE READ            '
C
      DO 50 ISPZ=1,NSPAMI
50      TXTPSP(ISPZ,15)=TEXTS(ISPZ)
C
      DO 80 IPLS=1,NPLSI
        ISPZ=NSPAMI+IPLS
        TXTPSP(IPLS,2)=TEXTS(ISPZ)
        TXTPSP(IPLS,4)=TEXTS(ISPZ)
        TXTPSP(IPLS,5)=TEXTS(ISPZ)
        TXTPSP(IPLS,6)=TEXTS(ISPZ)
        TXTPSP(IPLS,7)=TEXTS(ISPZ)
80      TXTPSP(IPLS,13)=TEXTS(ISPZ)
C
      RETURN
      END
 
 
 
 
 
 
 
