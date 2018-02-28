**==PPOUT
C
C=======================================================================
      SUBROUTINE PPOUT(IOPT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C***********************************************************************
C
C VERSION : V1.R1.M0
C
C ROUTINE : PPOUT
C
C PURPOSE : TO WRITE OUT DATA FOR POST-PROCESSOR
C
C INPUT   : (I*4) IOPT   - =1 WRITE DATA FOR CURRENT TIME SLICE
C                        - =2 WRITE TIME TRACES
C
C OUTPUT  : NONE
C
C HISTORY : V1.R1.M0 --- 12/03/93 --- CREATION
C
C***********************************************************************
C
      INCLUDE 'p01'
      INCLUDE 'p02'
      INCLUDE 'p09'
C
C../COMPPF/
      INCLUDE 'c00'
C../CINPUT/
      INCLUDE 'c02'
C../CPHYS/
      INCLUDE 'c03'
C../PARAM/
      INCLUDE 'c04'
C../PARAME/
      INCLUDE 'c05'
C../UNK/
      INCLUDE 'c06'
C../CADAS/
      INCLUDE 'c07'
C../CPOLYG/
      INCLUDE 'c09'
C../PSIDAT/
      INCLUDE 'c11'
C../C.FLUX/
      INCLUDE 'c32'
C
      REAL*8       BUFFER(MP),FLXC(MP,2),BUFF2D(5,MPG),DVOLL(MP,1)
      REAL*4       DEBWR(MDEBWR*MDEBPT*3)
      CHARACTER    CBUF6*6,CBUF8*8,CBUFF*80,GEONAM*176,SYSUID*8,PREFIX*7
      CHARACTER    ADDSP*8,ADDSPN*32
      CHARACTER*8  TCODE,TMACHID,TISHOT,TDATE
      CHARACTER    YEARH*2,YEARZ*2
      LOGICAL      LINIT,LPPFIN
C
      PARAMETER(NLOCDB=1)
      LOGICAL LINDEB,LOCDEB(NLOCDB),LHTRAN
C
      COMMON/CPPOUT/LPP,IFORM
C
      DATA LINDEB/.TRUE./,LOCDEB/NLOCDB*.FALSE./
      DATA LINIT/.TRUE./,LPPFIN/.FALSE./
      DATA IMAP/0/
      DATA SYSUID/'        '/
C
      SAVE LINDEB,LINIT,LPPFIN,LHTRAN
      SAVE LSUM
C
C-----------------------------------------------------------------------
C  SET UP DEBUG PRINT SWITCHES
C-----------------------------------------------------------------------
C
      IF(LINDEB.AND.LDBPRT)THEN
        LINDEB = .FALSE.
        CALL PPDEB(NLOCDB,LOCDEB)
        IF(LOCDEB(1)) IFORM=2
      ENDIF
C
C-----------------------------------------------------------------------
C  INITIALISATIONS
C-----------------------------------------------------------------------
C
      IF(LINIT)THEN
C
        LPP  = LTRANF
        LSUM = LSUMM
        IF(HWNAME.EQ.'IBM3090'   ) IFORM=0
        IF(HWNAME.EQ.'RS6000'    ) IFORM=1
        IF(HWNAME.EQ.'X86/LINUX' ) IFORM=1
        IF(IPPOUT.EQ.2           ) IFORM=0
C
C SET UP VOL ELEMENT ALLOWING PERIODIC POINTS
C
        DO 50 K=1,NP
          DVOLL(K,1) = DVOL(K,1)
          IF( ITAG(K,2).LT.0 )THEN
            NJL = ITAG(K,1)
            K1  = K - NJL + 1
            DVOLL(K,1) = DVOL(K1,1)
          ENDIF
 50    CONTINUE
C
       LHTRAN = LHVMAG_R(1)  .OR. LHVMAG_R(2)  .OR.
     &          LHVMAG_T(1)  .OR. LHVMAG_T(2)  .OR.
     &          LHVBDR_R(1)  .OR. LHVBDR_R(2)  .OR.
     &          LHVBDR_T(1)  .OR. LHVBDR_T(2)  .OR.
     &          LHVCENT_R(1) .OR. LHVCENT_R(2) .OR.
     &          LHVCENT_T(1) .OR. LHVCENT_T(2) .OR.
     &          LVEXB_R(1)   .OR. LVEXB_R(2)   .OR.
     &          LVEXB_T(1)   .OR. LVEXB_T(2)
C
      ENDIF
C
C-----------------------------------------------------------------------
C FIRST CALL ONLY
C-----------------------------------------------------------------------
C
      IF(LINIT)THEN
C
        LINIT = .FALSE.
C
C  WRITE CATALOGUE IDENTIFIER
C  --------------------------
C
        CALL NAMTRN(TCODE,TMACHID,TISHOT,TDATE)
        CALL PPCSTR( '*' )
        CALL PPCSTR( '*CATALOGUE IDENTIFIER' )
        CALL PPCSTR( ' ' // TCODE  // ' ' // TMACHID
     &                   // TISHOT // ' ' // TDATE   )
C
C  WRITE VERSION INFO
C  ------------------
C
        CALL PPCSTR( '*' )
        CALL PPCSTR( '*CODE VERSION' )
        CALL PPCSTR( ' VERSON :'//VERSON )
        CALL PPCSTR( ' SUPP   :'//SUPP )
        CALL PPCSTR( ' SUPP0  :'//SUPP0 )
        CALL PPCSTR( ' LINK   :'//LINK )
        CALL PPCSTR( ' NIMB   :'//NIMB )
        CALL PPCSTR( ' CORO   :'//CORO )
        CALL PPCSTR( ' CROS   :'//CROS )
        CALL PPCSTR( ' CRIBM  :'//CRIBM )
        CALL PPCSTR( ' MODEL  :'//MODEL )
C
C  WRITE USER INFO
C  ---------------
C
        CALL DMGUID(SYSUID,PREFIX)
        CALL PPCSTR( '*' )
        CALL PPCSTR( '*USER INFORMATION' )
        CALL PPCSTR( ' MACHINE:'//HWNAME )
        CALL PPCSTR( ' USER   :'//SYSUID )
C
C  PREPARE COMMENT SECTION
C  -----------------------
C
        CALL PPCSTR( '*' )
        CALL PPCSTR( '*RUN COMMENTS' )
C
C  WRITE CURRENT CONTROL DATA (FORMATTED)
C  --------------------------------------
C
        CALL PPCSTR( '*' )
        CALL PPCSTR( '*CONTROL DATA (FORMATTED)' )
        CALL PPCSTR( ' ADASIDH:'//USERID)
        WRITE(YEARH,'(I2)') IHYEAR
        CALL PPCSTR( ' YEARH  :'//YEARH )
        DO IZS=1,MAX(NZS,1)
          CBUF8 = ' '//ADDSP('ADASIDZ',IZS)
          CALL PPCSTR( CBUF8//':'//USERID)
          WRITE(YEARZ,'(I2)') IZYEAR(IZS)
          CBUF8 = ADDSP(' YEARZ',IZS)
          CALL PPCSTR( CBUF8//':'//YEARZ )
        ENDDO
C
C  WRITE GEOMETRY DATA
C  -------------------
C
        IF(NXPNT.NE.1)THEN
          MSG1 = 'UNEXPECTED GEOMETRY FOR POST-PROCESSOR'
          CALL ERRMSS(LOUT,'PPOUT',1,MSG1,' ',' ')
        ENDIF
C
        CALL PPCSTR( '*' )
        CALL PPCSTR( '*GEOMETRY' )
        CALL PPCSTR( ' '//FEQUIL )
C
        CALL SETGEO(NPL,IOPENL,NXWL,NCL,NROWL,JPRGTL,JPLFTL)
        CALL PPCSTR( 'NP,IOPEN,NXW,NC,NROW,JPRGT,JPLFT,IRFCT' )
        IF(IFORM.EQ.0)THEN
          WRITE(LPP)  NPL,IOPENL,NXWL,NCL,NROWL,JPRGTL,JPLFTL,IREFCT
        ELSE
          CBUFF = ' '
          WRITE(CBUFF,'(8I4)') NPL,IOPENL,NXWL,NCL,NROWL,JPRGTL,JPLFTL,
     &                         IREFCT
          CALL PPCSTR( CBUFF )
        ENDIF
C
        SF = 1.D-2
        O8 = 1.D0
        CALL PPRECW('RPX   ','X-POINT R-COORD ','M   ',RPX,SF,'R8',1,0)
        CALL PPRECW('ZPX   ','X-POINT Z-COORD ','M   ',ZPX,SF,'R8',1,0)
        CALL PPRECW('R0    ','AXIAL R-COORD   ','M   ',R0 ,SF,'R8',1,0)
        CALL PPRECW('Z0    ','AXIAL Z-COORD   ','M   ',Z0 ,SF,'R8',1,0)
        CALL PPRECW('B0    ','AXIAL B-TOROIDAL','T   ',B0 ,O8,'R8',1,0)
        CALL PPRECW('PSI   ','NORMALISED FLUX ',
     &              '    ',FLXPSI(1),O8,'R8',NP,0)
        CALL PPRECW('RSEPX ','R FOR SEPARATRIX',
     &              'M  ',RSEPX(1),SF,'R8',NSEPX,0)
        CALL PPRECW('ZSEPX ','Z FOR SEPARATRIX',
     &              'M  ',ZSEPX(1),SF,'R8',NSEPX,0)
        IF(NMAP2D.NE.1)THEN
          CALL PPRECW('ITAGDV',' ',' ',ITAGDV,1,'I4',NP,0)
        ENDIF
        CALL PPRECW('RMESH ','R-COORD   ','M  ',RMESH(1),SF,'R8',NP,0)
        CALL PPRECW('ZMESH ','Z-COORD   ','M  ',ZMESH(1),SF,'R8',NP,0)
        CALL PPRECW('RHO   ','RHO       ','   ',RHO(1)  ,O8,'R8',NP,0)
        CALL PPRECW('THETA ','THETA     ','   ',THETA(1),O8,'R8',NP,0)
        CALL PPRECW('HRHO  ','HRHO      ','   ',HRO(1)  ,O8,'R8',NP,0)
        CALL PPRECW('HTETA ','HTETA     ','   ',HTETA(1),O8,'R8',NP,0)
        CALL PPRECW('BFI   ','B-TOROIDAL','T  ',BFI(1),1.D-4,'R8',NP,0)
        CALL PPRECW('SH    ','BPOL/BTOT ','   ',SH(1)   ,O8,'R8',NP,0)
C
        CALL PPRECW('RVESM1','R-VESSEL 1','M  ',RVESM(1,1),1.E-2,'R4',
     &               NVESM,0)
        CALL PPRECW('ZVESM1','Z-VESSEL 1','M  ',ZVESM(1,1),1.E-2,'R4',
     &               NVESM,0)
        CALL PPRECW('RVESM2','R-VESSEL 2','M  ',RVESM(1,2),1.E-2,'R4',
     &               NVESM,0)
        CALL PPRECW('ZVESM2','Z-VESSEL 2','M  ',ZVESM(1,2),1.E-2,'R4',
     &               NVESM,0)
        CALL PPRECW('IPTVES','VESS STYLE','   ',IVESM(1,3),
     &               1,'I4',NVESM,0)
C
        IF( NVESP.GT.0 ) THEN
          CALL PPRECW('RPUMP1','R-PUMP   1','M  ',RVESM(NVESM+1,1),
     &                 1.E-2,'R4',NVESP,0)
          CALL PPRECW('ZPUMP1','Z-PUMP   1','M  ',ZVESM(NVESM+1,1),
     &                 1.E-2,'R4',NVESP,0)
          CALL PPRECW('RPUMP2','R-PUMP   2','M  ',RVESM(NVESM+1,2),
     &                 1.E-2,'R4',NVESP,0)
          CALL PPRECW('ZPUMP2','Z-PUMP   2','M  ',ZVESM(NVESM+1,2),
     &                 1.E-2,'R4',NVESP,0)
          CALL PPRECW('IPTPMP','PUMP STYLE','   ',IVESM(NVESM+1,3),
     &                 1,'I4',NVESP,0)
        END IF
C
        IF( NVESC.GT.0 ) THEN
          CALL PPRECW('XGAUGE','R-GAUGE   ','M  ',XVESC(1,1),
     &                1.E-2,'R4',NVESC,0)
          CALL PPRECW('YGAUGE','Z-GAUGE   ','M  ',XVESC(1,2),
     &                1.E-2,'R4',NVESC,0)
          CALL PPRECW('RGAUGE','RADIUS-GAUGE ','M  ',XVESC(1,3),
     &                1.E-2,'R4',NVESC,0)
        END IF
C
        CALL PPRECW('KORPG ',' ',' ',KORPG,1,'I4',NP,0)
        CALL PPRECW('NVERTP',' ',' ',NVERTP,1,'I4',NPOLYP,0)
C
C... SKIP TO AVOID RECORD LENGTH OVERFLOW ON THE MAINFRAME
        IF( TISHOT(1:1).NE.'I' ) THEN
         CALL MTRAN(MPG,5,RVERTP(1,1),BUFF2D(1,1))
         CALL PPRECW('RVERTP',' ',' ',BUFF2D(1,1),1.D-2,'R8',5*NPOLYP,0)
         CALL MTRAN(MPG,5,ZVERTP(1,1),BUFF2D(1,1))
         CALL PPRECW('ZVERTP',' ',' ',BUFF2D(1,1),1.D-2,'R8',5*NPOLYP,0)
        END IF
C
      ENDIF
C
C-----------------------------------------------------------------------
C RESULTS FOR CURRENT TIME STEP
C-----------------------------------------------------------------------
C
      IF(IOPT.EQ.1)THEN
C
C  WRITE HEADER FOR TIME STEP DATA
C  -------------------------------
C
        CALL PPCSTR( '*' )
        CALL PPCSTR( '*TIME STEP DATA' )
C
        CBUFF = ' '
        IF( LJETTO ) THEN
            WRITE(CBUFF,'(A,E11.5,A,I8,A,F8.4)')
     &                                    '-TIME  :',TIME
     &                                  , ' STEP  :',ISTEPA
     &                                  , ' JETTO TIME :',TIMJET(NTCOR)
        ELSE
            WRITE(CBUFF,'(A,E11.5,A,I8)')
     &                                    '-TIME  :',TIME
     &                                  , ' STEP  :',ISTEPA
        END IF
        CALL PPCSTR( CBUFF )
C
C BASIC QUANTITIES
C ----------------
C
        CALL PPRECW('HMASS ','ION MASS        ',
     &              'AMU   ',HMASS    ,1.D0,'R8',1,0)
        CALL PPRECW('HCH   ','ION ATOMIC  NUMBER',
     &              '      ',1        ,1,'I4',1,0)
        CALL PPRECW('NZ    ','NO. OF IMP. CHARGE STATES',
     &              '      ',NZC(1)   ,1,'I4',MAX(NZS,1),0)
        IF( NZ.GT.0 )THEN
          CALL PPRECW('ZMASS ','IMP.MASS   ',
     &                'AMU   ',ZMASS(1),1.D0,'R8',NZS,0)
          CALL PPRECW('ZCH   ','IMP.ATOMIC NUMBER ',
     &                '      ',IZ0(1) ,1,'I4',NZS,0)
        ENDIF
        CALL PPRECW('DEN   ','ION DENSITY     ',
     &              'M-3   ',DEN(1) ,1.D6,'R8',NP,0)
        CALL PPRECW('VPI   ','ION PAR.VEL     ',
     &              'M S-1 ',GAM(1)   ,1.D-2,'R8',NP,1)
        CALL PPRECW('PMACH ','MACH NUMBER (IONS)',
     &              '      ',PMACH    ,O8   ,'R8',NP,1)
        CALL PPRECW('TEV   ','ION TEMPERATURE ',
     &              'EV    ',TEV(1)   ,O8,'R8',NP,0)
        CALL PPRECW('DENEL ','EL. DENSITY     ',
     &              'M-3   ',DENEL(1) ,1.D6,'R8',NP,0)
        CALL PPRECW('TEVE  ','EL. TEMPERATURE ',
     &              'EV    ',TEVE(1)  ,O8,'R8',NP,0)
        CALL PPRECW('VRO   ','ION PERP.VEL    ',
     &              'M S-1 ',VRO(1)   ,1.D-2,'R8',NP,0)
        CALL VDIV(NP,URO0(1),ALFX5(1),BUFFER(1))
        CALL PPRECW('URO   ','ION CLASS.PERP.V',
     &              'M S-1 ',BUFFER(1),1.D-2,'R8',NP,0)
        CALL PPRECW('VPINCH','ION PINCH VELOCITY',
     &              'M S-1 ',VPINCH(1), 1.D-2,'R8',NP,0)
        IF( LHVMAG_R(1) .OR. LHVMAG_R(2) .OR.
     &      LHVMAG_T(1) .OR. LHVMAG_T(2)      )THEN
          CALL PPRECW('VROM  ','ION PERP.MAG.V',
     &                'M S-1 ',VMAGI(1,1),1.D-2,'R8',NP,0)
          CALL PPRECW('VPAM  ','ION PAR.MAG.V',
     &                'M S-1 ',VMAGI(1,2),1.D-2,'R8',NP,0)
          CALL PPRECW('VTRM  ','ION TRAN.MAG.V',
     &                'M S-1 ',VMAGI(1,3),1.D-2,'R8',NP,0)
        ENDIF
        IF( LHTRAN )
     &    CALL PPRECW('VTRI  ','ION TRAN.VEL',
     &                'M S-1 ',VTRANI(1),1.D-2,'R8',NP,0)
        CALL PPRECW('VPE   ','EL. PAR.VEL     ',
     &              'M S-1 ',VETE(1)  ,1.D-2,'R8',NP,0)
        CALL PPRECW('VROE  ','EL. PERP.VEL    ',
     &              'M S-1 ',VROE(1)  ,1.D-2,'R8',NP,0)
C
        IF( NZ.GT.0 )
     &    CALL PPRECW('ZEFF  ','Z EFFECTIVE',
     &              '      ',ZEFF(1)   ,O8,'R8',NP,0)
        CALL VDIV(NP,UROE0(1),ALFX5(1),BUFFER(1))
        CALL PPRECW('UROE  ','EL. CLASS.PERP.V',
     &              'M S-1 ',BUFFER(1) ,1.D-2,'R8',NP,0)
        IF( LHTRAN )
     &    CALL PPRECW('VTRE  ','EL. TRAN.VEL',
     &                'M S-1 ',VTRANE(1),1.D-2,'R8',NP,0)
        CALL PPRECW('JPAR  ','PARA.CURR.DENS. ',
     &              'A M-2 ',CURPAR(1)  ,1.D0/3.D5,'R8',NP,0)
        IF( LHTRAN )THEN
          CALL PPRECW('JRO   ','PERP.CURR.DENS.',
     &                'A M-2 ',CURDEN(1,1),1.D0/3.D5,'R8',NP,0)
          CALL PPRECW('JTR   ','TRAN.CURR.DENS.',
     &                'A M-2 ',CURDEN(1,3),1.D0/3.D5,'R8',NP,0)
        ENDIF
        CALL PPRECW('POT   ','PLASMA POTENTIAL',
     &              'V     ',POT(1)   ,3.D2,'R8',NP,0)
        CALL PPRECW('EPAR  ','PARA.ELEC.FIELD ',
     &              'V M-1 ',ELEFLD(1),3.D4,'R8',NP,0)
C
        CALL PPRECW('DA    ','NEUT. ATOM DENSITY',
     &              'M-3   ',DA(1)    ,1.D6,'R8',NP,0)
        CALL PPRECW('VA0    ','NEUT. ATOM VELOCITY    ',
     &              'M S-1 ',VNUTAP(1),1.D-2,'R8',NP,0)
        CALL PPRECW('VA0R   ','NEUT. ATOM VEL.R-COMP',
     &              'M S-1 ',VA(1,1),1.D-2,'R8',NP,0)
        CALL PPRECW('VA0Z   ','NEUT. ATOM VEL.Z-COMP',
     &              'M S-1 ',VA(1,2),1.D-2,'R8',NP,0)
        CALL PPRECW('VA0T   ','NEUT. ATOM VEL.TOR-COMP',
     &              'M S-1 ',VA(1,3),1.D-2,'R8',NP,0)
        CALL PPRECW('ENEUTA','NEUTRAL ATOM TEMP.',
     &              'EV    ',ENEUTA(1),0.6666D0/EV,'R8',NP,0)
        CALL PPRECW('DM    ','NEUT. MOL. DENSITY',
     &              'M-3   ',DM(1)    ,1.D6,'R8',NP,0)
        CALL PPRECW('ENEUTM','NEUTRAL MOL. TEMP.',
     &              'EV    ',ENEUTM(1),0.6666D0/EV,'R8',NP,0)
        CALL PPRECW('VM0    ','NEUT. MOL. VELOCITY',
     &              'M S-1 ',VNUTMP(1),1.D-2,'R8',NP,0)
        CALL PPRECW('VM0R   ','NEUT. MOL. VEL.R-COMP',
     &              'M S-1 ',VM(1,1),1.D-2,'R8',NP,0)
        CALL PPRECW('VM0Z   ','NEUT. MOL. VEL.Z-COMP',
     &              'M S-1 ',VM(1,2),1.D-2,'R8',NP,0)
        CALL PPRECW('VM0T   ','NEUT. MOL. VEL.TOR-COMP',
     &              'M S-1 ',VM(1,3),1.D-2,'R8',NP,0)
        IF( NZS.GT.0 )THEN
          DO IZS=1,NZS
            CALL PPRECW(ADDSP('DZ',IZS),
     &                  ADDSPN('NEUT. IMP. DENSITY',IZS),
     &              'M-3   ',DZ(1,IZS)    ,1.D6,'R8',NP,0)
          ENDDO
          DO IZS=1,NZS
            CALL PPRECW(ADDSP('ENEUTZ',IZS),
     &                  ADDSPN('NEUT. IMP. TEMP.',IZS),
     &              'EV    ',ENEUTZ(1,IZS),0.6666D0/EV,'R8',NP,0)
          ENDDO
        ENDIF
C
        CALL PPRECW('DPERP ','PART.DIFF.COEFF.',
     &              'M2 S-1',DPERP(1) ,1.D-4,'R8',NP,0)
        CALL PPRECW('CHII  ','ION THERMAL DIFFUSIVITY',
     &              'M2 S-1',CHIIL(1) ,1.D-4,'R8',NP,0)
        CALL PPRECW('CHIE  ','EL. THERMAL DIFFUSIVITY',
     &              'M2 S-1',CHIEL(1) ,1.D-4,'R8',NP,0)
        CALL PPRECW('ETAPAR ','PARA.VISCOSITY  ',
     &              'N S M-2',CETAH(1) ,1.D-1,'R8',NP,0)
        CALL PPRECW('ETAPER','PERP.VISCOSITY  ',
     &              'N S M-2',ETA(1)   ,1.D-1,'R8',NP,0)
        CALL PPRECW('DHA   ','H-ALPHA RAD. DENS.',
     &              'PH.M-3/S',DHA(1)    ,1.D6,'R8',NP,0)
        CALL PPRECW('PI    ','ION NON-HYD.STAT. PRESSURE ',
     &              'PA      ',PIH(1) ,1.D-1,'R8',NP,0)
C
C COMPUTED QUANTITIES
C -------------------
C
        FAC = 1.D6
        DO K=1,MP
          IF( K.LE.NP )THEN
            BUFFER(K) = DEN(K)*DA(K)*SVCX(K)
          ELSE
            BUFFER(K) = 0.0D0
          ENDIF
        ENDDO
        CALL PPRECW('RCX   ','CX RATE',
     &              'M-3 S-1',BUFFER(1),FAC,'R8',NP,0)
C
C ION PARTICLE BALANCE
C --------------------
C
        FAC = 1.D6/TWOPI
        CALL PPFLUX(FLXL,FLXR,FLXD,FLXU,FLXC)
        CALL PPRECW('PFLXL ','ION PART. FLUX - LEFT',
     &              'S-1   ',FLXC(1,1),1.0D0,'R8',NP,1)
        CALL PPRECW('PFLXD ','ION PART. FLUX - DOWN',
     &              'S-1   ',FLXC(1,2),1.0D0,'R8',NP,1)
        CALL VDIV(NP,PION(1),DVOLL(1,1),BUFFER(1))
        CALL PPRECW('SOUN  ','IONISATION SOURCE   ',
     &              'M-3 S-1',BUFFER(1),FAC,'R8',NP,0)
        CALL VDIV(NP,PREC(1),DVOLL(1,1),BUFFER(1))
        CALL PPRECW('SIREC ','RECOMBIN.  SOURCE   ',
     &              'M-3 S-1',BUFFER(1),FAC,'R8',NP,0)
        CALL VDIV(NP,PEXT(1),DVOLL(1,1),BUFFER(1))
        CALL PPRECW('SEXT  ','EXTERNAL ION SOURCE  ',
     &              'M-3 S-1',BUFFER(1),FAC,'R8',NP,0)
        CALL VDIV(NP,DERN(1),DVOLL(1,1),BUFFER(1))
        CALL PPRECW('SIDT  ','D/DT(ION DENS.)     ',
     &              'M-3 S-1',BUFFER(1),FAC,'R8',NP,0)
C
C IMPURITY PARTICLE BALANCE
C -------------------------
C
        IF( NZS.GT.0 )THEN
          DO IZS=1,NZS
            CALL PPFLUX(FLXZLS(1,IZS),FLXZRS(1,IZS),
     &                  FLXZDS(1,IZS),FLXZUS(1,IZS),FLXC)
            CALL PPRECW(ADDSP('PFLXZL',IZS),
     &                  ADDSPN('IMP.PART. FLUX - LEFT',IZS),
     &                'S-1   ',FLXC(1,1),1.0D0,'R8',NP,1)
            CALL PPRECW(ADDSP('PFLXZD',IZS),
     &                  ADDSPN('IMP.PART. FLUX - DOWN',IZS),
     &                'S-1   ',FLXC(1,2),1.0D0,'R8',NP,1)
            CALL VDIV(NP,PSRCZS(1,IZS),DVOLL(1,1),BUFFER(1))
            CALL PPRECW(ADDSP('SOUZ  ',IZS),
     &                  ADDSPN('IMP.IONISATION SOURCE',IZS),
     &                'M-3 S-1',BUFFER(1),FAC,'R8',NP,0)
            CALL VDIV(NP,DERNZS(1,IZS),DVOLL(1,1),BUFFER(1))
            CALL PPRECW(ADDSP('SZDT  ',IZS),
     &                  ADDSPN('D/DT(IMP DENS.)     ',IZS),
     &                'M-3 S-1',BUFFER(1),FAC,'R8',NP,0)
          ENDDO
        ENDIF
C
C ION MOMENTUM BALANCE
C --------------------
C
        FAC = 1.0D-5*HMASS*AMU
        CALL PPFLUX(GIL,GIR,GID,GIU,FLXC)
        CALL PPRECW('MFLXL ','ION.MOM. FLUX - LEFT',
     &              'N     ',FLXC(1,1),FAC,'R8',NP,0)
        CALL PPRECW('MFLXD ','ION.MOM. FLUX - DOWN',
     &              'N     ',FLXC(1,2),FAC,'R8',NP,0)
        FAC = FAC/TWOPI*1.D6
        CALL VDIV(NP,RF(1),DVOLL(1,1),BUFFER(1))
        CALL PPRECW('RIF  ','ION MOMENTUM SOURCE',
     &              'N M-3  ',BUFFER(1),FAC,'R8',NP,1)
        CALL VDIV(NP,RPG(1),DVOLL(1,1),BUFFER(1))
        CALL PPRECW('RIPG ','ION PRESSURE GRADIENT',
     &              'N M-3  ',BUFFER(1),FAC,'R8',NP,1)
        CALL VDIV(NP,REF(1),DVOLL(1,1),BUFFER(1))
        CALL PPRECW('RIEF ','E-FIELD FORCE ON IONS',
     &              'N M-3  ',BUFFER(1),FAC,'R8',NP,1)
        CALL VDIV(NP,RFR(1),DVOLL(1,1),BUFFER(1))
        CALL PPRECW('RIFR ','FRICTIONAL FORCE ON IONS',
     &              'N M-3  ',BUFFER(1),FAC,'R8',NP,1)
        CALL VDIV(NP,RTH(1),DVOLL(1,1),BUFFER(1))
        CALL PPRECW('RITH ','THERMOELEC. FORCE ON IONS',
     &              'N M-3  ',BUFFER(1),FAC,'R8',NP,1)
        CALL VDIV(NP,RSI(1),DVOLL(1,1),BUFFER(1))
        CALL PPRECW('RISI ','ION SEMI-IMPLICIT TERM',
     &              'N M-3  ',BUFFER(1),FAC,'R8',NP,1)
        CALL VDIV(NP,REXT(1),DVOLL(1,1),BUFFER(1))
        CALL PPRECW('RIEXT ','ION EXTERNAL MOM. SOURCE',
     &              'N M-3  ',BUFFER(1),FAC,'R8',NP,1)
        CALL VDIV(NP,DMIDT(1),DVOLL(1,1),BUFFER(1))
        CALL PPRECW('RIDT ','D/DT(ION MOMENTUM)',
     &              'N M-3  ',BUFFER(1),FAC,'R8',NP,1)
C
C IMPURITY MOMENTUM BALANCE
C -------------------------
C
        IF( NZS.GT.0 )THEN
          DO IZS=1,NZS
            FAC = 1.0D-5*ZMASS(IZS)*AMU
            CALL PPFLUX(GZL(1,IZS),GZR(1,IZS),
     &                  GZD(1,IZS),GZU(1,IZS),FLXC)
            CALL PPRECW(ADDSP('MFLXZL',IZS),
     &                  ADDSPN('IMP.MOM. FLUX - LEFT',IZS),
     &                  'N     ',FLXC(1,1),FAC,'R8',NP,0)
            CALL PPRECW(ADDSP('MFLXZD',IZS),
     &                  ADDSPN('IMP.MOM. FLUX - DOWN',IZS),
     &                  'N     ',FLXC(1,2),FAC,'R8',NP,0)
            FAC = FAC/TWOPI*1.D6
            CALL VDIV(NP,RFZ(1,IZS),DVOLL(1,1),BUFFER(1))
            CALL PPRECW(ADDSP('RZF  ',IZS),
     &                  ADDSPN('IMP.MOMENTUM SOURCE',IZS),
     &                  'N M-3  ',BUFFER(1),FAC,'R8',NP,1)
            CALL VDIV(NP,RPGZ(1,IZS),DVOLL(1,1),BUFFER(1))
            CALL PPRECW(ADDSP('RZPG ',IZS),
     &                  ADDSPN('IMP.PRESSURE GRADIENT',IZS),
     &                'N M-3  ',BUFFER(1),FAC,'R8',NP,1)
            CALL VDIV(NP,REFZ(1,IZS),DVOLL(1,1),BUFFER(1))
            CALL PPRECW(ADDSP('RZEF ',IZS),
     &                  ADDSPN('E-FIELD FORCE ON IMP.',IZS),
     &                  'N M-3  ',BUFFER(1),FAC,'R8',NP,1)
            CALL VDIV(NP,RFRZ(1,IZS),DVOLL(1,1),BUFFER(1))
            CALL PPRECW(ADDSP('RZFR ',IZS),
     &                  ADDSPN('FRICTIONAL FORCE ON IMP.',IZS),
     &                  'N M-3  ',BUFFER(1),FAC,'R8',NP,1)
            CALL VDIV(NP,RTHZ(1,IZS),DVOLL(1,1),BUFFER(1))
            CALL PPRECW(ADDSP('RZTH ',IZS),
     &                  ADDSPN('THERMOELEC. FORCE ON IMP.',IZS),
     &                  'N M-3  ',BUFFER(1),FAC,'R8',NP,1)
            CALL VDIV(NP,RSIZ(1,IZS),DVOLL(1,1),BUFFER(1))
            CALL PPRECW(ADDSP('RZSI ',IZS),
     &                  ADDSPN('IMP.SEMI-IMPLICIT TERM',IZS),
     &                'N M-3  ',BUFFER(1),FAC,'R8',NP,1)
            CALL VDIV(NP,RSPLTZ(1,IZS),DVOLL(1,1),BUFFER(1))
            CALL PPRECW(ADDSP('RZSPLT',IZS),
     &                  ADDSPN('IMP. SPLIT SOURCE TERM',IZS),
     &                  'N M-3  ',BUFFER(1),FAC,'R8',NP,1)
            CALL VDIV(NP,DMZDT(1,IZS),DVOLL(1,1),BUFFER(1))
            CALL PPRECW(ADDSP('RZDT ',IZS),
     &                  ADDSPN('D/DT(IMP.MOMENTUM)',IZS),
     &                  'N M-3  ',BUFFER(1),FAC,'R8',NP,1)
          ENDDO
        ENDIF
C
C ION ENERGY BALANCE
C ------------------
C
        FAC = 1.0D12/TWOPI
        CALL PPFLUX(FLXIL,FLXIR,FLXID,FLXIU,FLXC)
        CALL PPRECW('QIFLXL','ION.ENERGY FLUX - LEFT',
     &              'W     ',FLXC(1,1),1.0D6,'R8',NP,1)
        CALL PPRECW('QIFLXD','ION.ENERGY FLUX - DOWN',
     &              'W     ',FLXC(1,2),1.0D6,'R8',NP,1)
        CALL VDIV(NP,SQINT(1,3),DVOLL(1,1),BUFFER(1))
        CALL PPRECW('SQIA  ','ION ATOM. IONISATION',
     &              'W M-3  ',BUFFER(1),FAC,'R8',NP,0)
        CALL VDIV(NP,SQINT(1,4),DVOLL(1,1),BUFFER(1))
        CALL PPRECW('SQIM  ','ION MOL.IONIS/DISSOC',
     &              'W M-3  ',BUFFER(1),FAC,'R8',NP,0)
        CALL VDIV(NP,SQINT(1,5),DVOLL(1,1),BUFFER(1))
        CALL PPRECW('SQICX ','ION CHARGE EXCHANGE ',
     &              'W M-3  ',BUFFER(1),FAC,'R8',NP,0)
        CALL VDIV(NP,SQINT(1,6),DVOLL(1,1),BUFFER(1))
        CALL PPRECW('SQIELA','ATOMIC EL.SCATTERING',
     &              'W M-3  ',BUFFER(1),FAC,'R8',NP,0)
        CALL VDIV(NP,SQINT(1,7),DVOLL(1,1),BUFFER(1))
        CALL PPRECW('SQIELM','MOL.EL.SCATTERING',
     &              'W M-3  ',BUFFER(1),FAC,'R8',NP,0)
        CALL VDIV(NP,SQRCI(1)  ,DVOLL(1,1),BUFFER(1))
        CALL PPRECW('SQIRC ','ION VOL.RECOMBINATION',
     &              'W M-3  ',BUFFER(1),FAC,'R8',NP,0)
        CALL VDIV(NP,QEI(1)   ,DVOLL(1,1),BUFFER(1))
        CALL PPRECW('SQIEQ ','ION EQUIPARTITION    ',
     &              'W M-3  ',BUFFER(1),FAC,'R8',NP,0)
        CALL VADD(NP,QCMRXI(1),QCMRYI(1),BUFFER(1))
        CALL VDIV(NP,BUFFER(1),DVOLL(1,1),BUFFER(1))
        CALL PPRECW('SQICMR','ION COMPRESSION TERM ',
     &              'W M-3  ',BUFFER(1),FAC,'R8',NP,0)
        CALL VDIV(NP,QEXTI(1),DVOLL(1,1),BUFFER(1))
        CALL PPRECW('SQIEXT','ION EXTERNAL ENERGY SOURCE',
     &              'W M-3  ',BUFFER(1),FAC,'R8',NP,0)
        CALL VDIV(NP,DPIDT(1)  ,DVOLL(1,1),BUFFER(1))
        CALL PPRECW('SQIDT ','D/DT(ION ENERGY DENS.) ',
     &              'W M-3  ',BUFFER(1),FAC,'R8',NP,0)
C
C ELECTRON ENERGY BALANCE
C -----------------------
C
        CALL PPFLUX(FLXEL,FLXER,FLXED,FLXEU,FLXC)
        CALL PPRECW('QEFLXL','EL. ENERGY FLUX - LEFT',
     &              'W     ',FLXC(1,1),1.0D6,'R8',NP,1)
        CALL PPRECW('QEFLXD','EL. ENERGY FLUX - DOWN',
     &              'W     ',FLXC(1,2),1.0D6,'R8',NP,1)
        CALL VDIV(NP,SQENT(1,4),DVOLL(1,1),BUFFER(1))
        CALL PPRECW('SQEA  ','EL. ATOM. IONISATION',
     &              'W M-3  ',BUFFER(1),FAC,'R8',NP,0)
        CALL VDIV(NP,SQENT(1,6),DVOLL(1,1),BUFFER(1))
        CALL PPRECW('SQEM  ','EL. MOL.IONIS/DISSOC',
     &              'W M-3  ',BUFFER(1),FAC,'R8',NP,0)
        CALL VDIV(NP,SQENT(1,5),DVOLL(1,1),BUFFER(1))
        CALL PPRECW('SQEHRAD','HYDROG.RADIATION    ',
     &              'W M-3  ',BUFFER(1),FAC,'R8',NP,0)
        IF( NZS.LE.1 )THEN
          CALL VDIV(NP,SQENT(1,2),DVOLL(1,1),BUFFER(1))
          CALL PPRECW('SQEZRAD','IMPURITY RADIATION',
     &              'W M-3  ',BUFFER(1),FAC,'R8',NP,0)
        ELSE
          DO IZS=1,NZS
            CALL VDIV(NP,SQE2Z(1,IZS),DVOLL(1,1),BUFFER(1))
            CALL PPRECW(ADDSP('SQEZRAD',IZS),
     &                  ADDSPN('IMPURITY RADIATION',IZS),
     &                'W M-3  ',BUFFER(1),FAC,'R8',NP,0)
          ENDDO
        ENDIF
        CALL VDIV(NP,SQRCE(1)  ,DVOLL(1,1),BUFFER(1))
        CALL PPRECW('SQERC ','EL. VOL.RECOMBINATION',
     &              'W M-3  ',BUFFER(1),FAC,'R8',NP,0)
        CALL VDIV(NP,QIE(1)   ,DVOLL(1,1),BUFFER(1))
        CALL PPRECW('SQEEQ ','EL. EQUIPARTITION    ',
     &              'W M-3  ',BUFFER(1),FAC,'R8',NP,0)
        CALL VADD(NP,QCMRXE(1),QCMRYE(1),BUFFER(1))
        CALL VDIV(NP,BUFFER(1) ,DVOLL(1,1),BUFFER(1))
        CALL PPRECW('SQECMR','EL. COMPRESSION TERM',
     &              'W M-3  ',BUFFER(1),FAC,'R8',NP,0)
        CALL VDIV(NP,QEXTE(1),DVOLL(1,1),BUFFER(1))
        CALL PPRECW('SQEEXT','EL. EXTERNAL ENERGY SOURCE',
     &              'W M-3  ',BUFFER(1),FAC,'R8',NP,0)
        CALL VDIV(NP,DPEDT(1)  ,DVOLL(1,1),BUFFER(1))
        CALL PPRECW('SQEDT ','D/DT(EL. ENERGY DENS.) ',
     &              'W M-3  ',BUFFER(1),FAC,'R8',NP,0)
C
C IMPURITY STAGES
C ---------------
C
        DO 200 IZ=1,NZ
          CBUF6 = ' '
          WRITE(CBUF6,'(A,I2.2)') 'IMP ',IZ
          CBUF8 = ' '
          WRITE(CBUF8,'(A,I2.2)') 'DENZ',IZ
          CALL PPRECW(CBUF8,CBUF6//' DENSITY  ',
     &                'M-3' ,DENZ(1,IZ) ,1.D6,'R8',NP,0)
          CBUF8 = ' '
          WRITE(CBUF8,'(A,I2.2)') 'VPZ',IZ
          CALL PPRECW(CBUF8,CBUF6//' PAR.VEL  ',
     &                'M S-1' ,GAMZ(1,IZ) ,1.D-2,'R8',NP,1)
          CBUF8 = ' '
          WRITE(CBUF8,'(A,I2.2)') 'VROZ',IZ
          CALL PPRECW(CBUF8,CBUF6//' PERP.VEL  ',
     &                'M S-1 ',VROZ(1,IZ) ,1.D-2,'R8',NP,0)
          IF( LHTRAN )THEN
            CBUF8 = ' '
            WRITE(CBUF8,'(A,I2.2)') 'VTRZ',IZ
            CALL PPRECW(CBUF8,CBUF6//' TRAN.VEL  ',
     &                  'M S-1 ',VTRANZ(1,IZ),1.D-2,'R8',NP,0)
          ENDIF
          CBUF8 = ' '
          WRITE(CBUF8,'(A,I2.2)') 'ZI',IZ
          CALL PPRECW(CBUF8,CBUF6//' CHARGE   ',
     &                '    ',ZI(1,IZ) ,O8,'R8',NP,0)
          CBUF8 = ' '
          WRITE(CBUF8,'(A,I2.2)') 'ZISQ',IZ
          CALL PPRECW(CBUF8,CBUF6//' CHARGE**2',
     &                '    ',ZI2(1,IZ) ,O8,'R8',NP,0)
C
 200    CONTINUE
C
C  WRITE HEADER FOR NIMBUS WALL DATA
C  ---------------------------------
C
        CALL PPCSTR( '*' )
        CALL PPCSTR( '*NIMBUS WALL' )
C
        CBUFF = ' '
        WRITE(CBUFF,'(A,E11.5,A,I8)') '-TIME :',
     &               TIMMC(NTINE),' STEP :',ISTNEU(NTINE)
        CALL PPCSTR( CBUFF )
C
        CALL SETMAR(NVESM,GVESM(1)
     &             ,IOTS,IOTE,IODS,IODE,IIDS,IIDE,IITS,IITE
     &             ,IPRS,IPRE)
        CALL PPCSTR(
     &       '%10I IOTS,IOTE,IODS,IODE,IIDS,IIDE,IITS,IITE,IPRS,IPRE')
        IF(IFORM.EQ.0)THEN
          WRITE(LPP)  IOTS,IOTE,IODS,IODE,IIDS,IIDE,IITS,IITE,IPRS,IPRE
        ELSE
          CBUFF = ' '
          WRITE(CBUFF,'(10I4)') IOTS,IOTE,IODS,IODE,IIDS,IIDE,IITS,IITE
     &                        , IPRS,IPRE
          CALL PPCSTR( CBUFF )
        ENDIF
C
C QUANTITIES
C ----------
C
        CALL PPRECW('FLXHW1 ','H + H2 WALL FLUX',
     &              'M-2 S-1',FLUXHW(1,3),1.D4,'R8',NFLXHW,0)
        CALL PPRECW('FLXHW2 ','H + H+ WALL FLUX',
     &              'M-2 S-1',FLUXHW(1,4),1.D4,'R8',NFLXHW,0)
        CALL PPRECW('FLXHW3 ','Z-SPUTTERING WALL FLUX',
     &              'M-2 S-1',FLUXHW(1,5),1.D4,'R8',NFLXHW,0)
        CALL PPRECW('FLXHW4 ','Z-REDEPOSITION WALL FLUX',
     &              'M-2 S-1',FLUXHW(1,6),1.D4,'R8',NFLXHW,0)
        CALL PPRECW('FLXHW5 ','AVG. H-NEUTRAL ENERGY ON WALL',
     &              'EV     ',FLUXHW(1,7),1.D0,'R8',NFLXHW,0)
        CALL PPRECW('FLXHW6 ','H WALL FLUX',
     &              'M-2 S-1',FLUXHW(1,8),1.D4,'R8',NFLXHW,0)
        CALL PPRECW('FLXHW7 ','H IMPLANTATION FLUX',
     &              'M-2 S-1',FLUXHW(1,9),1.D4,'R8',NFLXHW,0)
C
C  WRITE HEADER FOR NIMBUS BOUNDARY
C  --------------------------------
C
        CALL PPCSTR( '*' )
        CALL PPCSTR( '*NIMBUS BOUNDARY' )
C
        CALL SETMAR(NPLASM,GPLASM(1)
     &             ,IOTS,IOTE,IODS,IODE,IIDS,IIDE,IITS,IITE
     &             ,IPRS,IPRE)
        CALL PPCSTR(
     &       '%10I IOTS,IOTE,IODS,IODE,IIDS,IIDE,IITS,IITE,IPRS,IPRE')
        IF(IFORM.EQ.0)THEN
          WRITE(LPP) IOTS,IOTE,IODS,IODE,IIDS,IIDE,IITS,IITE,IPRS,IPRE
        ELSE
          CBUFF = ' '
          WRITE(CBUFF,'(10I4)') IOTS,IOTE,IODS,IODE,IIDS,IIDE,IITS,IITE
     &                        , IPRS,IPRE
          CALL PPCSTR( CBUFF )
        ENDIF
C
C QUANTITIES
C ----------
C
        CALL PPRECW('RPLAM1','R-PLASMA BOUNDARY 1','M  ',RPLASM(1,1),
     &               1.E-2,'R4',NPLASM,0)
        CALL PPRECW('ZPLAM1','Z-PLASMA BOUNDARY 1','M  ',ZPLASM(1,1),
     &               1.E-2,'R4',NPLASM,0)
        CALL PPRECW('RPLAM2','R-PLASMA BOUNDARY 2','M  ',RPLASM(1,2),
     &               1.E-2,'R4',NPLASM,0)
        CALL PPRECW('ZPLAM2','Z-PLASMA BOUNDARY 2','M  ',ZPLASM(1,2),
     &               1.E-2,'R4',NPLASM,0)
C
        CALL PPRECW('FLXPB1 ','H + H2 PLASMA BOUNDARY FLUX',
     &              'M-2 S-1',FLUXPB(1,1),1.D4,'R8',NPLASM,0)
C
        IF( NJPOLY(1).GT.0 ) THEN
            DO J = 1 , NJPOLY(1)
               BUFFER(J) = RPOLY(1,J)
            ENDDO
            CALL PPRECW( 'RBCORE'  , 'R-CORE BOUNDARY' , 'M  '
     &                 , BUFFER(1) , 1.D-2 , 'R8'
     &                 , NJPOLY(1) , 0 )
            DO J = 1 , NJPOLY(1)
               BUFFER(J) = ZPOLY(1,J)
            ENDDO
            CALL PPRECW( 'ZBCORE'  , 'Z-CORE BOUNDARY' , 'M  '
     &                 , BUFFER(1) , 1.D-2 , 'R8'
     &                 , NJPOLY(1) , 0 )
        END IF
C
C  WRITE HEADER FOR NIMBUS OUTPUT
C  ------------------------------
C
        CALL PPCSTR( '*' )
        CALL PPCSTR( '*NIMBUS OUTPUT' )
C
C NEUTRAL TRAJECTORIES
C --------------------
C
        IF( IDEBWR.GT.0 ) THEN
C
            CALL PPRECW( 'NTHIST'
     &                 , 'HISTORY ASSOC. WITH NEUT. TRAJ.'
     &                 , '   ' , NDBHST(1) , 1 , 'I4' , IDEBWR , 0 )
            CALL PPRECW( 'NTTRACK'
     &                 , 'NO. OF TRACKS IN NEUT. TRAJ.'
     &                 , '   ' , KDEBWR(1) , 1 , 'I4' , IDEBWR , 0 )
C
C           PACK (X,Y,Z) COORDINATES INTO A 1-D ARRAY
            L = 0
            DO I = 1 , IDEBWR
               DO J = 1 , KDEBWR(I)
                  DO K = 1 , 3
                     L = L + 1
                     DEBWR(L) = XDEBWR(I,J,K)
                  ENDDO
               ENDDO
            ENDDO
C
            CALL PPRECW( 'NTCOORD'
     &                 , 'NEUT. TRAJ. (X,Y,Z) COORDS'
     &                 , 'M  ' , DEBWR(1) , 1.0E-02 , 'R4' , L , 0 )
C
        ENDIF
C
C-----------------------------------------------------------------------
C TIME TRACES
C-----------------------------------------------------------------------
C
      ELSE IF(IOPT.EQ.2)THEN
C
        IF(LPPFIN)THEN
          MSG1 = '2ND TIME TRACE DATA WRITE ATTEMPT'
          CALL ERRMSS(LOUT,'PPOUT',1,MSG1,' ',' ')
        ENDIF
        LPPFIN = .TRUE.
C
        CALL PPCSTR( '*' )
        CALL PPCSTR( '*TIME TRACES' )
        CALL PPRECW('ITSTEP','TRACE STEP NO.S ',
     &              '  ',ISTEPE(1),1,'I4',NEVOL,0)
        CALL PPRECW('TVEC  ','TRACE TIMES     ',
     &              'MS',TIMEE(1)  ,1.,'R4',NEVOL,0)
        CALL PPRECW('HCONTE','MAIN ION CONTENT',
     &              '  ',HCONTE(1),1.,'R4',NEVOL,0)
        IF( NZ.GT.0)THEN
          DO IZS=1,NZS
            CALL PPRECW(ADDSP('ZCONTE',IZS),
     &                  ADDSPN('IMPURITY CONTENT',IZS),
     &             '  ',ZCONTE(1,IZS),1.,'R4',NEVOL,0)
          ENDDO
        ENDIF
        CALL PPRECW('HMID  ','H MID-PLANE SEPX DENS',
     &              'M-3  ',HMID(1),1.D6,'R8',NEVOL,0)
        CALL PPRECW('DENSOT','ION DENSITY (OUTER TARGET)',
     &              'M-3   ',DENSOT(1),1.E6,'R4',NEVOL,0)
        CALL PPRECW('TISOT ','ION TEMPER. (OUTER TARGET)',
     &              'EV    ',TISOT(1),1.,'R4',NEVOL,0)
        CALL PPRECW('DENSIT','ION DENSITY (INNER TARGET)',
     &              'M-3   ',DENSIT(1),1.E6,'R4',NEVOL,0)
        CALL PPRECW('TISIT ','ION TEMPER. (INNER TARGET)',
     &              'EV    ',TISIT(1),1.,'R4',NEVOL,0)
        IF( NZ.GT.0 )THEN
          DO IZS=1,NZS
            CALL PPRECW(ADDSP('ZRAD',IZS),
     &                  ADDSPN('IMPURITY RADIATION',IZS),
     &              'MW    ',ZRAD(1,IZS),1.,'R4',NEVOL,0)
          ENDDO
        ENDIF
        CALL PPRECW('DEESOT','ELE DENSITY (OUTER TARGET)',
     &              'M-3   ',DEESOT(1),1.E6,'R4',NEVOL,0)
        CALL PPRECW('TESOT ','ELE TEMPER. (OUTER TARGET)',
     &              'EV    ',TESOT(1),1.,'R4',NEVOL,0)
        CALL PPRECW('DEESIT','ELE DENSITY (INNER TARGET)',
     &              'M-3   ',DEESIT(1),1.E6,'R4',NEVOL,0)
        CALL PPRECW('TESIT ','ELE TEMPER. (INNER TARGET)',
     &              'EV    ',TESIT(1),1.,'R4',NEVOL,0)
        CALL PPRECW('EMID  ','ELE MID-PLANE SEPX DENS',
     &              'M-3   ',EMID(1),1.E6,'R4',NEVOL,0)
        CALL PPRECW('TIMID ','ION MID-PLANE SEPX TEMP',
     &              'EV    ',TIMID(1),1.,'R4',NEVOL,0)
        CALL PPRECW('TEMID ','ELE MID-PLANE SEPX TEMP',
     &              'EV    ',TEMID(1),1.,'R4',NEVOL,0)
        CALL PPRECW('HCFACD','HCONT FRACTIONAL DIFF',
     &              '      ',HCFACD(1),1.,'R4',NEVOL,0)
        CALL PPRECW('PSFACD','PART. SRC. FACT. DIFF',
     &              '      ',PSFACD(1),1.,'R4',NEVOL,0)
        CALL PPRECW('QIFACD','ION ENERGY FACT. DIFF',
     &              '      ',QIFACD(1),1.,'R4',NEVOL,0)
        CALL PPRECW('QEFACD','ELE ENERGY FACT. DIFF',
     &              '      ',QEFACD(1),1.,'R4',NEVOL,0)
        CALL PPRECW('HBND  ','H AVERAGE CORE BOUND.DENS',
     &              'M-3  ',HBND(1),1.E6,'R4',NEVOL,0)
        CALL PPRECW('TIBND ','AVERAGE CORE BOUNDARY ION TEMP.',
     &              'EV   ',TIBND(1),1.,'R4',NEVOL,0)
        CALL PPRECW('TEBND  ','AVERAGE CORE BOUNDARY EL TEMP.',
     &              'EV   ',TEBND(1),1.,'R4',NEVOL,0)
        CALL PPRECW('FLUXIM ','AVERAGE CORE BOUNDARY ION FLUX',
     &              'S-1  ',FNBND(1),1.,'R4',NEVOL,0)
        CALL PPRECW('FIBND  ','AVERAGE CORE BOUNDARY ION POWER',
     &              'W    ',FIBND(1),1.E6,'R4',NEVOL,0)
        CALL PPRECW('FEBND  ','AVERAGE CORE BOUNDARY  EL.POWER',
     &              'W    ',FEBND(1),1.E6,'R4',NEVOL,0)
        CALL PPRECW('FLUX0M','NEUTRAL FLUX TO M.P.',
     &              'S-1   ',FLUX0M(1),1.,'R4',NEVOL,0)
        CALL PPRECW('FLUX0P','NEUTRAL FLUX TO PUMP',
     &              'S-1   ',FLUX0P(1),1.,'R4',NEVOL,0)
        CALL PPRECW('FLXWHP','WALL ION FLUX PUMPED',
     &              'S-1   ',FLXWHP(1),1.,'R4',NEVOL,0)
        CALL PPRECW('HEXTRA','EXTRA NEUTRAL FLUX',
     &              'S-1   ',HEXTRT(1),1.,'R4',NEVOL,0)
        IF( NZ.GT.0)THEN
          DO IZS=1,NZS
            CALL PPRECW(ADDSP('SPFAC ',IZS),
     &                  ADDSPN('SPUTTERING FACTOR',IZS),
     &              '      ',SPFAC(1,IZS),1.,'R4',NEVOL,0)
          ENDDO
          DO IZS=1,NZS
            CALL PPRECW(ADDSP('YLDT ',IZS),
     &                  ADDSPN('EFF. SPUTTERING YIELD',IZS),
     &              '      ',YLDT(1,IZS),1.,'R4',NEVOL,0)
          ENDDO
          DO IZS=1,NZS
            CALL PPRECW(ADDSP('FLUZ0M',IZS),
     &                  ADDSPN('IMP.NEUTRAL FLUX TO M.P.',IZS),
     &              '      ',FLUZ0M(1,IZS),1.,'R4',NEVOL,0)
          ENDDO
          DO IZS=1,NZS
            CALL PPRECW(ADDSP('FLUZ0P',IZS),
     &                  ADDSPN('IMP.NEUTRAL FLUX TO PUMP',IZS),
     &              '      ',FLUZ0P(1,IZS),1.,'R4',NEVOL,0)
          ENDDO
          DO IZS=1,NZS
            CALL PPRECW(ADDSP('GAMZT',IZS),
     &                  ADDSPN('IMP. ION FLUX OUT',IZS),
     &              '      ',GAMZT(1,IZS),1.,'R4',NEVOL,0)
          ENDDO
          DO IZS=1,NZS
            CALL PPRECW(ADDSP('ZEXTRT',IZS),
     &                  ADDSPN('EXTRA NEUTRAL IMPURITIES',IZS),
     &              '      ',ZEXTRT(1,IZS),1.,'R4',NEVOL,0)
          ENDDO
          DO IZS=1,NZS
            CALL PPRECW(ADDSP('ZBND',IZS),
     &                  ADDSPN('AVERAGE CORE BOUND.IMP.DENS',IZS),
     &              'M-3   ',ZBND(1,IZS),1.E6,'R4',NEVOL,0)
          ENDDO
          DO IZS=1,NZS
            CALL PPRECW(ADDSP('FZBND',IZS),
     &                  ADDSPN('CORE BOUNDARY IMP.FLUX',IZS),
     &              'S-1   ',FZBND(1,IZS),1.,'R4',NEVOL,0)
          ENDDO
        ENDIF
        CALL PPRECW('GAMIT','FLUX OUT',
     &              '      ',GAMIT(1),1.,'R4',NEVOL,0)
        IF( LJETTO )THEN
          CALL PPRECW('HTOTAL','CORE+EDGE PLASMA CONTENT',
     &                '  ',HTOTAL(1),1.,'R4',NEVOL,0)
          CALL PPRECW('HCONTJ','JETTO CORE PLASMA CONTENT',
     &                '  ',HCONTJ(1),1.,'R4',NEVOL,0)
          CALL PPRECW('HAVGJ ','JETTO CORE AVER. DENS.',
     &                'M-3   ',HAVGJ(1),1.E6,'R4',NEVOL,0)
          CALL PPRECW('HBNDJ ','JETTO BOUNDARY DENS.',
     &                'M-3   ',HBNDJ(1),1.E6,'R4',NEVOL,0)
          CALL PPRECW('TIBNDJ ','JETTO BOUNDARY ION TEMP.',
     &                'EV   ',TIBNDJ(1),1.,'R4',NEVOL,0)
          CALL PPRECW('TEBNDJ ','JETTO BOUNDARY EL TEMP.',
     &                'EV   ',TEBNDJ(1),1.,'R4',NEVOL,0)
          CALL PPRECW('FNBNDJ ','JETTO BOUNDARY ION FLUX',
     &                'S-1  ',FNBNDJ(1),1.,'R4',NEVOL,0)
          CALL PPRECW('FIBNDJ ','JETTO BOUNDARY ION POWER',
     &                'W    ',FIBNDJ(1),1.E6,'R4',NEVOL,0)
          CALL PPRECW('FEBNDJ ','JETTO BOUNDARY EL.POWER',
     &                'W    ',FEBNDJ(1),1.E6,'R4',NEVOL,0)
          CALL PPRECW('DPBNDJ ','JETTO BOUNDARY DPERP',
     &                'M2 S-1',DPBNDJ(1),1.E-4,'R4',NEVOL,0)
          CALL PPRECW('VPINJ ','JETTO BOUNDARY VPINCH',
     &                'M S-1',VPINJ(1),1.E-2,'R4',NEVOL,0)
          CALL PPRECW('XIBNDJ ','JETTO BOUNDARY CHII',
     &                'M2 S-1',XIBNDJ(1),1.E-4,'R4',NEVOL,0)
          CALL PPRECW('XEBNDJ ','JETTO BOUNDARY CHIE',
     &                'M2 S-1',XEBNDJ(1),1.E-4,'R4',NEVOL,0)
          CALL PPRECW('WITHJ ','JETTO CORE ION ENERGY',
     &                'J     ',WITHJ(1),1.,'R4',NEVOL,0)
          CALL PPRECW('WETHJ ','JETTO CORE EL. ENERGY',
     &                'J     ',WETHJ(1),1.,'R4',NEVOL,0)
          IF( NZ.GT.0 )THEN
            DO IZS=1,NZS
              CALL PPRECW(ADDSP('FZSN',IZS),
     &                    ADDSPN('SANCO BOUNDARY IMP.FLUX',IZS),
     &                'S-1   ',FZSN(1,IZS),1.,'R4',NEVOL,0)
            ENDDO
            DO IZS=1,NZS
              CALL PPRECW(ADDSP('DPSN',IZS),
     &                    ADDSPN('SANCO BOUNDARY IMP.DPERP',IZS),
     &                'M2 S-1',DPSN(1,IZS),1.E-4,'R4',NEVOL,0)
            ENDDO
            DO IZS=1,NZS
              CALL PPRECW(ADDSP('VPINSN',IZS),
     &                    ADDSPN('SANCO BOUNDARY IMP.VPINCH',IZS),
     &                'M S-1',VPINSN(1,IZS),1.E-2,'R4',NEVOL,0)
            ENDDO
            DO IZS=1,NZS
              CALL PPRECW(ADDSP('ZTOTAL',IZS),
     &                    ADDSPN('CORE+EDGE IMP.CONTENT',IZS),
     &                '     ',ZTOTAL(1,IZS),1.,'R4',NEVOL,0)
            ENDDO
            DO IZS=1,NZS
              CALL PPRECW(ADDSP('ZLEVSN',IZS),
     &                    ADDSPN('SANCO CORE IMP.CONTENT',IZS),
     &                '     ',ZLEVSN(1,IZS),1.,'R4',NEVOL,0)
            ENDDO
            DO IZS=1,NZS
              CALL PPRECW(ADDSP('ZRADSN',IZS),
     &                    ADDSPN('SANCO CORE IMP.RADIATION',IZS),
     &                'MW ',ZRADSN(1,IZS),1.E-6,'R4',NEVOL,0)
            ENDDO
          ENDIF
        ENDIF
        CALL PPRECW('FLXOT','ION FLUX TO OUTER TARGET',
     &              '      ',FLXOT(1),1.,'R4',NEVOL,0)
        CALL PPRECW('FLXIT','ION FLUX TO INNER TARGET',
     &              '      ',FLXIT(1),1.,'R4',NEVOL,0)
C
        CALL PPRECW('HBMID ','H MID-PLANE BND DENS.',
     &              'M-3  ',HBMID(1),1.E6,'R4',NEVOL,0)
        CALL PPRECW('HBTOP ','H TOP-PLANE BND DENS.',
     &              'M-3  ',HBTOP(1),1.E6,'R4',NEVOL,0)
C
        CALL PPCSTR( '*' )
        CALL PPCSTR( '*EOF' )
C
        LINIT = .TRUE.
        LPPFIN = .FALSE.
C
      ENDIF
C
C-----------------------------------------------------------------------
C
 9999 RETURN
C
      END
**++EOF
**==ADDSP
C
C=======================================================================
      FUNCTION ADDSP(VAR,IZS)
      IMPLICIT REAL*8(A-H,O-Z)
C
C***********************************************************************
C
C VERSION : V1.R1.M0
C
C ROUTINE : ADDSP
C
C PURPOSE : ADD IMPURITY SPECIES INDEX TO POST PROCESSOR VARIABLE NAME
C
C INPUT   : (C* ) VAR    - VARIABLE NAME
C           (I*4) IZS    - SPECIES INDEX
C
C OUTPUT  : (C*8) ADDSP  - CONSTRUCTED NAME
C
C HISTORY : V1.R1.M0 --- 30/10/96 --- CREATION
C
C***********************************************************************
C
      INCLUDE 'p01'
      INCLUDE 'p02'
C
      CHARACTER*(*) VAR
      CHARACTER*2   CBUF2
      CHARACTER*8   ADDSP
C
C../PARAM/
      INCLUDE 'c04'
C
      IF(LENSTR(VAR) .GT. 6 .AND.
     &       (VAR .NE. 'SQEZRAD' .AND. VAR .NE. 'ADASIDZ'))
     &       CALL ERRMSS(LOUT
     &           ,'PPOUT',1,'VARIABLE NAME ',VAR,' IS TOO LONG')
C
      IF( NZS.GT.1 )THEN
         LNAME = LENSTR(VAR)
         IF(LNAME .EQ. 7 ) LNAME = 5
         WRITE(CBUF2,'(A,I1)') '_',IZS
         ADDSP=VAR(1:LNAME)//CBUF2
      ELSE
         ADDSP=VAR
      ENDIF
C
C
      RETURN
      END
**++EOF
**==ADDSPN
C
C=======================================================================
      FUNCTION ADDSPN(DESC,IZS)
      IMPLICIT REAL*8(A-H,O-Z)
C
C***********************************************************************
C
C VERSION : V1.R1.M0
C
C ROUTINE : ADDSPN
C
C PURPOSE : ADD IMPURITY SPECIES NAME TO POST PROCESSOR DESCRIPTION
C
C INPUT   : (C* ) DESC   - DESCRIPTION
C           (I*4) IZS    - SPECIES INDEX
C
C OUTPUT  : (C*32) ADDSP  - CONSTRUCTED DESCRIPTION
C
C HISTORY : V1.R1.M0 --- 30/10/96 --- CREATION
C
C***********************************************************************
C
      INCLUDE 'p01'
      INCLUDE 'p02'
C
      INCLUDE 'c02'
C
      CHARACTER*(*) DESC
      CHARACTER*32  ADDSPN
      CHARACTER     CBUF*2,ELEMEN*15,ELETMP*15
C
C../PARAM/
      INCLUDE 'c04'
C
      ELETMP = ELEMEN(IZ0(IZS))
      ADDSPN = ELETMP(14:15)//'- '//DESC
C
C
      RETURN
      END
**++EOF
**==CONVER
C
C=======================================================================
      SUBROUTINE CONVER(IMAP,NDATA,RDATA,NDATN,RDATN)
      IMPLICIT REAL*8(A-H,O-Z)
C
C***********************************************************************
C
C VERSION : V1.R1.M0
C
C ROUTINE : CONVER
C
C PURPOSE : CONVERT DATA TO K SYSTEM
C
C INPUT   : (I*4) IMAP   - CONVERSION TYPE -   1 -> (II,JJ) TO K
C           (I*4) NDATA  - NUMBER OF DATA POINTS
C           ( - ) RDATA  - DATA (MAY BE I*4, R*4)
C           (I*4) NDATN  - DIMENSION OF RDATN
C
C OUTPUT  : (I*4) NDATN  - NUMBER OF DATA POINTS (NP)
C           ( - ) RDATN  - DATA (MAY BE I*4, R*4)
C
C NOTES   : ADDITIONAL IMAP'S MAY BE REQUIRED FOR OTHER CODES
C
C HISTORY : V1.R1.M0 --- 12/03/93 --- CREATION
C
C***********************************************************************
C
      INCLUDE 'p01'
      INCLUDE 'p02'
C
C../PARAM/
      INCLUDE 'c04'
C
      DIMENSION RDATA(MX,MY),RDATN(MP)
C
C
C CHECK INPUT
C -----------
C
      IF(NDATN.LT.NP)THEN
        CALL ERRMSS(LOUT,'CONVER',1,'RDATN ARRAY TOO SMALL',' ',' ')
      ENDIF
C
      IF(IMAP.NE.1)THEN
        CALL ERRMSS(LOUT,'CONVER',1,'UNKNOWN IMAP',' ',' ')
      ENDIF
C
      IF(IMAP.EQ.1 .AND. NDATA.NE.NI2D*NJ2D)THEN
        CALL ERRMSS(LOUT,'CONVER',1,
     &    'NDATA IS INCONSISTENT WITH IMAP',' ',' ')
      ENDIF
C
C CONVERT FROM (II,JJ)
C --------------------
C
      IF(IMAP.EQ.1)THEN
C
        NDATN = NP
        DO 100 K=1,NP
          II = IKOR(K)
          JJ = JKOR(K)
          IF(II.GT.0 .AND. JJ.GT.0)THEN
            RDATN(K) = RDATA(II,JJ)
          ENDIF
 100    CONTINUE
C
       ENDIF
C
C
      RETURN
      END
**++EOF
**==NAMTRN
C
C=======================================================================
      SUBROUTINE NAMTRN(TCODE,TMACHID,TISHOT,TDATE)
      IMPLICIT REAL*8(A-H,O-Z)
C
C***********************************************************************
C
C VERSION : V1.R4.M0
C
C ROUTINE : NAMTRN
C
C PURPOSE : TO CREATE LABELS USED FOR CATALOGUED TRANSFER FILE NAME
C
C INPUT   : NONE
C
C OUTPUT  : (C*8) TCODE   - CODE NAME LABEL
C           (C*8) TMACHID - MACHINE LABEL
C           (C*8) TISHOT  - SHOT LABEL
C           (C*8) TDATE   - DATE LABEL
C
C HISTORY : V1.R1.M0 --- 01/06/94 --- CREATION
C           V1.R2.M0 --- 26/09/94 --- APPEND EQUILIBRIUM VERSION TO
C                                     TISHOT AS TWO CHARACTER IDENTIFIER
C                                     (VX, WHERE X=1 TO 9 FOR VERS1->9
C                                                X=A TO Z FOR VERS10->35
C                                     ). IF THE STRING 'VERSXX'  DOES
C                                     NOT COMPLETE THE EQUILIBRIUM NAME
C                                     THEN NO VERSION WILL BE APPENDED.
C           V1.R3.M9 --- 07/12/94 --- ISHOT MAY BE DIFFERENT FROM NSHOT
C                                     WHERE, NSHOT IS READ FROM THE
C                                     GEOMETRY FILE
C           V1.R4.M0 --- 14/06/01 --- TMACHID
C
C***********************************************************************
C
      INCLUDE 'p01'
      INCLUDE 'p02'
C
      INCLUDE 'c00'
      INCLUDE 'c02'
      INCLUDE 'c04'
      INCLUDE 'c10'
C
C DUMMY ARGUMENTS
C
      CHARACTER*8  TCODE,TMACHID,TISHOT,TDATE
      CHARACTER*8  NDATE
      CHARACTER*80 DATIME
      CHARACTER*1  CITOA(35)
      CHARACTER*20 TEQUIL
      CHARACTER*2  CTEMP
C
      NAMELIST/PPROC/NDATE
C
      DATA CITOA/ '1' , '2' , '3' , '4' , '5', '6', '7', '8', '9'
     *          , 'a' , 'b' , 'c' , 'd' , 'e', 'f', 'g', 'h', 'i'
     *          , 'j' , 'k' , 'l' , 'm' , 'n', 'o', 'p', 'q', 'r'
     *          , 's' , 't' , 'u' , 'v' , 'w', 'x', 'y', 'z'/
c
c  DEFAULT ASSIGNMENT
C  ------------------
C
      TCODE = CODNAM
      CALL CHRUTL(TCODE,TCODE)
C
      IF( LENSTR(MACHID).LE.0 ) THEN
          IF( ISHOT.GE.90000. AND. ISHOT.LT.91000 ) THEN
              TMACHID = 'ITER'
          ELSE IF (ISHOT.GE.91000. AND. ISHOT.LT.92000) THEN
              TMACHID = 'ALCATOR'
          ELSE
              TMACHID = 'JET'
          ENDIF
      ELSE
          TMACHID = MACHID(1:8)
      ENDIF
      CALL CHRUTL(TMACHID,TMACHID)
C
      CALL CHRLTU(XEQUIL,TEQUIL)
C
      WRITE(TISHOT,'(I5)') ISHOT
C
C.. POSITION OF LAST CHARACTER IN TISHOT
      ITSHOT  = 0
      DO 100 I = 1 , 8
         IF( TISHOT(I:I).NE.' '               ) ITSHOT = I
  100 CONTINUE
C
C....................... GET VERSION OF GEOMETRY .......................
C          (APPEND VERSION NUMBER OF EQUILIBRIUM ONTO TISHOT)
C
      IF( ISHOT.EQ.NSHOT ) THEN
C
        IF( TEQUIL(JEQUIL+1:JEQUIL+4).EQ.'VERS' .AND. ITSHOT.LE.6 ) THEN
            IEQUIL = 0
            DO 110 I = 5 , 8
               IF( TEQUIL(JEQUIL+I:JEQUIL+I).GE.'0' .AND.
     &             TEQUIL(JEQUIL+I:JEQUIL+I).LE.'9'       ) THEN
                   IEQUIL = I
               ELSE IF( TEQUIL(JEQUIL+I:JEQUIL+I).NE.' ' ) THEN
                   GOTO 10
               END IF
  110       CONTINUE
            CTEMP = TEQUIL(JEQUIL+5:JEQUIL+IEQUIL)
            READ(CTEMP,'(I2)') I
            IF( I.GT.0 .AND. I.LE.35 )
     >          TISHOT = TISHOT(1:ITSHOT)//'v'//CITOA(I)
        END IF
C
      ELSE
C
        TISHOT = TISHOT(1:ITSHOT)//'X'
C
      END IF
C
C.......................................................................
C
  10  CONTINUE
C
      CALL TODAY(DATIME, 0)
      DO 120 I=17,25
        IF(DATIME(I:I).EQ.' ') DATIME(I:I) = '0'
 120  CONTINUE
      NDATE = DATIME(20:22)//DATIME(17:18)//DATIME(24:25)
C
C  OVERRIDE WITH NAMELIST
C  ----------------------
C
      CALL LOCSTR( LINP , '&PPROC' , .TRUE. , IRC )
C
      IF( IRC.GT.0 ) THEN
C
          IF( IRC.NE.2 ) CALL ERRMSS(LOUT,'PPOUT',1
     &      , '&PPROC NAMELIST LABEL MUST START ON 2ND COLUMN OF LINE'
     &      , ' ' , ' ' )
C
          READ(LINP,PPROC,IOSTAT=IERROR)
C
          IF( IERROR.GT.0 ) THEN
              CALL ERRMSS(LOUT
     &               ,'PPOUT',1,'NAMELIST /PPROC/ IS IN ERROR',' ',' ')
          ENDIF
      ENDIF
C
      TDATE = NDATE
      CALL CHRUTL( TDATE, TDATE )
C
 9999 RETURN
      END
**++EOF
**==PPFLUX
C
C=======================================================================
      SUBROUTINE PPFLUX(FLUXL,FLUXR,FLUXD,FLUXU,FLUXC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C***********************************************************************
C
C VERSION : V1.R1.M0
C
C ROUTINE : PPFLUX
C
C PURPOSE : TO COMBINE (L/R) AND (D/U) FLUX ARRAYS FOR POST PROCESSOR
C           IN COMBINED ARRAYS  A) POSITIVE => (L TO R) AND (D TO U)
C                               B) STAGGERED FLUXES REMAIN STAGGERED
C
C INPUT   : (R*8) FLUXL  - LEFT  BOUNDARY INWARD FLUX (K ORDERED)
C           (R*8) FLUXR  - RIGHT BOUNDARY INWARD FLUX
C           (R*8) FLUXD  - DOWN  BOUNDARY INWARD FLUX
C           (R*8) FLUXU  - UP    BOUNDARY INWARD FLUX
C
C OUTPUT  : (R*8) FLUXC  - (K,1) LEFT  BOUNDARY FLUX
C                          (K,2) DOWN  BOUNDARY FLUX
C
C HISTORY : V1.R1.M0 --- 02/08/94 --- CREATION
C
C***********************************************************************
C
      INCLUDE 'p01'
      INCLUDE 'p02'
      INCLUDE 'c04'
      INCLUDE 'c30'
C
      DIMENSION FLUXL(MP),FLUXR(MP),FLUXD(MP),FLUXU(MP),FLUXC(MP,2)
C
      DO 100 K=1,NP
        IF( IKOR(K).EQ.NI2D )THEN
          K1         = KORXYS(NI2D-1,JKOR(K))
          FLUXC(K,1) = -FLUXR(K1)
        ELSE
          FLUXC(K,1) = FLUXL(K)
        ENDIF
        IF( JKOR(K).EQ.NJ2D )THEN
          II         = IKOR(K)
          K1         = KORXYS(II,JM1P(II,JKOR(K)))
          FLUXC(K,2) = -FLUXU(K1)
        ELSE
          FLUXC(K,2) = FLUXD(K)
        ENDIF
 100  CONTINUE
C
 9999 RETURN
      END
**++EOF
**==PPDEB
C
C=======================================================================
      SUBROUTINE PPDEB(NFLG,LDBLST)
      IMPLICIT NONE
C
C***********************************************************************
C
C VERSION : V1.R1.M0
C
C ROUTINE : PPDEB
C
C PURPOSE : TO SET UP DEBUG CONTROL FLAGS FOR PPOUT ROUTINES - THIS
C           IS AN INTERFACE ROUTINE TO DEBFLG SYSTEM FOR EDGE2D/U.
C
C INPUT   : (I*4) NFLG  - NUMBER OF DEBUG SWITCHES
C
C OUTPUT  : (L)   LDBLST- DEBUG SWITCHES FOR PPOUT ROUTINES
C
C NOTES   : THIS WILL BE CODE DEPENDENT
C
C HISTORY : V1.R1.M0 --- 12/03/93 --- CREATION
C
C***********************************************************************
C
C DUMMY ARGUMENTS
C
      INTEGER NFLG
      LOGICAL LDBLST(NFLG)
C
C
      CALL DEBFLG('PPOUT   ',NFLG,LDBLST)
C
 9999 RETURN
      END
**++EOF
**==PPSUM
C
C=======================================================================
      SUBROUTINE PPSUM(LUN)
      IMPLICIT REAL*8(A-H,O-Z)
C
C***********************************************************************
C
C VERSION : V1.R1.M0
C
C ROUTINE : PPSUM
C
C PURPOSE : TO WRITE SUMMARY DATA SECTION TO POST PROCESSOR FILE
C           A TEMPORARY DATASET ON  LOGICAL UNIT LSUM IS USED
C           TO PERMIT FORMATTED TO BINARY CONVERSION
C
C INPUT   : (I*4) LUN    - LOGICAL UNIT NUMBER OF TEMPORARY FILE
C
C OUTPUT  : NONE
C
C HISTORY : V1.R1.M0 --- 02/08/94 --- CREATION
C
C***********************************************************************
C
      CHARACTER*133 BUFFER
C
      REWIND LUN
      CALL SUMPRT(LUN)
      REWIND LUN
C
 10   READ(LUN,'(A)',END=100) BUFFER
      CALL PPCSTR(BUFFER)
      GOTO 10
 100  CONTINUE
C
      RETURN
      END
**++EOF
**==SETGEO
C
C=======================================================================
      SUBROUTINE SETGEO(NPL,IOPENL,NXWL,NCL,NROW,JPRGTL,JPLFTL)
      IMPLICIT REAL*8(A-H,O-Z)
C
C***********************************************************************
C
C VERSION : V1.R1.M0
C
C ROUTINE : SETGEO
C
C PURPOSE : TO SET UP GEOMETRY GRID PARAMETERS FOR POST PROCESSOR FILE
C
C INPUT   : NONE
C
C OUTPUT  : (I*4) NPL    - TOTAL NUMBER OF K VALUES
C           (I*4) IOPENL - FIRST OPEN RING
C           (I*4) NXWL   - WALL RING
C           (I*4) NCL    - NUMBER OF RINGS
C           (I*4) NROW   - NUMBER OF ROWS
C           (I*4) JPRGTL - FIRST ROW INTERSECTING MAIN PLASMA
C           (I*4) JPLFTL - LAST  ROW INTERCEPTING MAIN PLASMA
C
C HISTORY : V1.R1.M0 --- 12/03/93 --- CREATION
C
C***********************************************************************
C
      INCLUDE 'p01'
      INCLUDE 'p02'
C
C../PARAM/
      INCLUDE 'c04'
C../PARAME/
      INCLUDE 'c05'
C
C
C DERIVE REQUIRED OUTPUT
C ----------------------
C
      NROW   = NR + 1
      IOPENL = IOPEN
      NXWL   = NXW
      NCL    = NC
      JPRGTL = JPRGT
      JPLFTL = JPLFT
C
C CHECK NP IS CONSISTENT
C ----------------------
C
      NMP  = (IOPEN-1)*(JPLFT-JPRGT+2)
      NPRV = (JPRGT-1+NROW-JPLFT)*(NC-NXW)
      NSOL = NROW*(NXW-IOPEN+1)
      NPL  = NMP+NPRV+NSOL
      IF( NPL.NE.NP )THEN
        MSG1 = 'UNEXPECTED NP'
        MSG2 =  '    NP,IOPEN,NXW,NC,NROW,JPRGT,JPLFT'
        WRITE(MSG3,*)  NPL,IOPENL,NXWL,NCL,NROW,JPRGTL,JPLFTL
        CALL ERRMSS(LOUT,'SETGEO',1,MSG1,MSG2,MSG3)
      ENDIF
C
C
      RETURN
      END
**++EOF
**==SETMAR
C
C=======================================================================
      SUBROUTINE SETMAR( NSEG,GSEG
     O                 , IOTS,IOTE,IODS,IODE,IIDS,IIDE,IITS,IITE
     O                 , IPRS,IPRE)
      IMPLICIT REAL*8(A-H,O-Z)
C
C***********************************************************************
C
C VERSION : V1.R2.M0
C
C ROUTINE : SETMAR
C
C PURPOSE : TO SET UP SEGMENT MARKING PARAMETERS FOR POST PROCESSOR FILE
C
C INPUT   : (I*4) NSEG   - NUMBER OF SEGMENTS
C           (C**) GSEG() - 'OT' ... OUTER TARGET
C                          'OC' ... OUTER CORNER
C                          'OD' ... OUTER DIVERTOR
C                          'MS' ... MAIN SOL
C                          'ID' ... INNER DIVERTOR
C                          'IC' ... INNER CORNER
C                          'IT' ... INNER TARGET
C                          'PV' ... PRIVATE VOID
C                          'BA' ... BAFFLE
C                          'CB' ... COMPOUND BAFFLE
C
C OUTPUT  : (I*4) IOTS   - FIRST SEGMENT ON OUTER TARGET
C           (I*4) IOTE   - LAST  SEGMENT ON OUTER TARGET
C           (I*4) IODS   - FIRST SEGMENT ON OUTER DIVERTOR
C           (I*4) IODE   - LAST  SEGMENT ON OUTER DIVERTOR
C           (I*4) IIDS   - FIRST SEGMENT ON INNER DIVERTOR
C           (I*4) IIDE   - LAST  SEGMENT ON INNER DIVERTOR
C           (I*4) IITS   - FIRST SEGMENT ON INNER TARGET
C           (I*4) IITE   - LAST  SEGMENT ON INNER TARGET
C           (I*4) IPRS   - FIRST SEGMENT ON PRIVATE REGION
C           (I*4) IPRE   - LAST  SEGMENT ON PRIVATE REGION
C
C NOTE    : ORDER IS CLOCK-WISE RANGING STARTING AT SEGMENT # 1
C
C HISTORY : V1.R1.M0 --- 06/04/94 --- CREATION
C           V1.R2.M0 --- 08/07/97 --- IPRS & IPRE INTO ARGUMENT LIST
C
C***********************************************************************
C
      CHARACTER GSEG(NSEG)*(*)
C
C
C INITIALISE QUANTITIES
C ---------------------
C
C
      IOTS = 1
      IOTE = 0
      IODS = 0
      IODE = 0
      IIDS = 0
      IIDE = 0
      IITS = 0
      IITE = 0
      IPRS = 0
      IPRE = 0
C
C
C DERIVE REQUIRED OUTPUT
C ----------------------
C
      DO 100 I = 1 , NSEG-1
         IF( GSEG(I).NE.'OT' .AND. GSEG(I+1).EQ.'OT' ) IOTS = I+1
         IF( GSEG(I).EQ.'OT' .AND. GSEG(I+1).NE.'OT' ) IOTE = I
         IF( GSEG(I).NE.'OD' .AND. GSEG(I+1).EQ.'OD' ) IODS = I+1
         IF( GSEG(I).EQ.'OD' .AND. GSEG(I+1).NE.'OD' ) IODE = I
         IF( GSEG(I).NE.'ID' .AND. GSEG(I+1).EQ.'ID' ) IIDS = I+1
         IF( GSEG(I).EQ.'ID' .AND. GSEG(I+1).NE.'ID' ) IIDE = I
         IF( GSEG(I).NE.'IT' .AND. GSEG(I+1).EQ.'IT' ) IITS = I+1
         IF( GSEG(I).EQ.'IT' .AND. GSEG(I+1).NE.'IT' ) IITE = I
         IF( GSEG(I).NE.'PV' .AND. GSEG(I+1).EQ.'PV' ) IPRS = I+1
         IF( GSEG(I).EQ.'PV' .AND. GSEG(I+1).NE.'PV' ) IPRE = I
  100 CONTINUE
C
      IF( IPRE.EQ.0 ) IPRE = NSEG
C
      RETURN
      END
**++EOF
