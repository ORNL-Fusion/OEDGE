c     -*Former Mode Specification*-
c
      subroutine DIVTRN(nizs,iter,niters,facta,factb,
     >                  title,job,equil,DESC,jfcb)
      IMPLICIT NONE
C
C  *********************************************************************
C  *                                                                   *
C  *  DIVTRN:  READS DIVIMP RESULTS AND REWRITES THEM INTO A TRANSFER  *
C  *           FILE COMPATIBLE WITH THE NEW POST-PROCESSOR.            *
C  *                                                                   *
C  *            LORNE HORTON  (JET)  FEBRUARY 1994                     *
C  *                                                                   *
C  *********************************************************************
C
C
C---- DIVIMP COMMON BLOCKS (LOADED BY CALL TO GET)
C
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'cneut'
      include 'dynam1'
      include 'dynam2'
      include 'dynam3'
      include 'dynam4'
      include 'pindata'
C
C---- INITIALISATION VARIABLES
C
      INTEGER   NIZS,ITER,NITERS
      REAL      TIME,TIME1,ZA02AS
      REAL      FACTA(-1:MAXIZS),FACTB(-1:MAXIZS)
      CHARACTER*(*) TITLE,JOB,EQUIL,DESC,jfcb
c
c     Local variables
c
      integer ik,ir,iz
C
C---- INITIALISATION
C
      TIME1 = ZA02AS (1)
      CALL XUFLOW (0)

c
c     Already loaded when called as a subroutine 
c
c     May lose support for iterated DIVIMP results when called
c     as a subroutine since only the current iteration results
c     will be written - need to check into this. This is not 
c     a commonly used DIVIMP feature though.
c
c
C---- LOAD DIVIMP RESULTS INTO COMMON BLOCKS
C
c      REWIND (8)
c   10 CONTINUE
c      CALL GET (TITLE,NIZS,JOB,EQUIL,FACTA,FACTB,ITER,NITERS)
c
c      Copy DDLIMS and DDTS to SDLIMS and SDTS for single precision write 
c
      call rzero(sdlims, MAXNKS*MAXNRS*(MAXIZS+2))
      call rzero(sdts  , MAXNKS*MAXNRS*(MAXIZS+2))
c
      do iz = -1,nizs
         do ir = 1,nrs
            do ik = 1,nks(ir)  
c
               sdlims(ik,ir,iz) = ddlims(ik,ir,iz)
               sdts(ik,ir,iz) = ddts(ik,ir,iz)
c
            end do
         end do
      end do
c
C
c      IF (ITER.EQ.1) THEN
c        RIZB = REAL (CIZB)
C
C---- WRITE TITLE TO PRINT FILE
C
c        WRITE ( 6,9010) VERSON,TITLE(1:58),JOB(1:36),JOB(37:72),
c     >                  EQUIL(1:54)
c        WRITE (52,9010) VERSON,TITLE(1:58),JOB(1:36),JOB(37:72),
c     >                  EQUIL(1:54)
c      ELSE
c        WRITE ( 6,9011) TITLE(61:80)
c        WRITE (52,9011) TITLE(61:80)
c      ENDIF
c
c     This routine is modeled on the EDGE2D one
C
C---- CALL PPOUT TO WRITE VALUES TO TRANSFER FILE FOR THIS ITERATION
C
      CALL PPOUT(NIZS,JOB,EQUIL,ITER,desc,jfcb)
C
c      IF (ITER.LT.NITERS) GOTO 10
C
      TIME = ZA02AS (1) - TIME1
      WRITE (6,9100) TIME
c
c      STOP
C
c 9010 FORMAT(/1X,62('*'),/1X,'*',60X,'*',
c     >  /1X,'*',16X,'RUN OF DIVTRN VERSION ',A5,17X,'*',
c     >  /1X,'*',18X,24('-'),18X,'*',/1X,'*',60X,'*',
c     >  /1X,'* ',A58,' *',/1X,'*',60X,'*',/1X,'* ',A36,22X,' *',
c     >  /1X,'*',60X,'*',/1X,'* ',A36,22X,' *',
c     >  /1X,'*',60X,'*',/1X,'* ',A54,4X,' *',
c     >  /1X,'*',60X,'*',/1X,62('*'),/)
c 9011 FORMAT(/1X,62('*'),/1X,'*',60X,'*',
c     >  /1X,'*',18X,A20,22X,'*',/1X,'*',60X,'*',/1X,62('*'),/)
 9100 FORMAT(/1X,'DIVTRN: TOTAL TIME USED = ',G11.4,' SEC',/)
C
      END
**==PPOUT
C
C=======================================================================
      SUBROUTINE PPOUT(NIZS,JOB,EQUIL,ITER,desc,jfcb)
      implicit none 
C
C***********************************************************************
C
C VERSION : V1.R1.M0
C
C ROUTINE : PPOUT
C
C PURPOSE : TO WRITE OUT DIVIMP DATA FOR POST-PROCESSOR
C
C INPUT   : (I*4)  NIZS  - MAXIMUM IMPURITY CHARGE STATE WHICH WAS
C                          FOLLOWED IN THIS RUN
C         : (C*72) JOB   - INFORMATION STRING WHICH SPECIFIES INPUT
C                          INPUT FILE NAME AND TIME OF DIVIMP RUN
C         : (C*72) EQUIL - STRING CONTAINING EQUILIBRIUM FILE NAME
C         : (I*4)  ITER  - ITERATION NUMBER TO BE WRITTEN
c         : (C*(*)) DESC - case description string 
C
C OUTPUT  : NONE
C
C HISTORY : V1.R1.M0 --- 22/02/94 --- CREATED FROM EDGE2D EQUIVALENT
C           V1.R1.M1 --- 31/01/97 --- TRANSLATE ADAS IDS TO JET3090
C
C***********************************************************************
C
      INTEGER       NIZS,ITER
      CHARACTER*(*) JOB,EQUIL,desc,jfcb
C
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'cneut'
      include 'dynam2'
      include 'dynam3'
      include 'dynam4'
      include 'pindata'
      include 'cadas'   
c slmod begin
      INCLUDE 'slcom'
c slmod end
C
      CHARACTER  CBUF6*6,CBUF8*8,CBUFF*80,SYSUID*7,PREFIX*7,TMPUID*7
      CHARACTER  TCODE*8,TSHOT*10,TDATE*8,tmachid*8
c      REAL*4     RSEPX(MAXNKS+1),ZSEPX(MAXNKS+1)
      REAL*4     TMPMAT(MAXNKS,MAXNRS),TMPMAT2(MAXNKS,MAXNRS)
      REAL*4     TMPVEC(MAXNDS)
      LOGICAL    LINIT
      INTEGER    L1, L2
c
c     Declare the variables passed to SETGEO
c
      integer npl, iopenl, nxwl,ncl,nrowl,jprgtl,jplftl
      integer nsepx
      real rsepx(maxnks+1),zsepx(maxnks+1)     
c
c     Declare the rest of the local variables
c
      real o8,time,sf
      integer ik,ir,iz,i
      integer iots,iote,iods,iode,iids,iide,iits,iite,iprs,ipre
c
c     integer     luid
c      DATA       LUID/10/
c
c     Local variables
c
      integer in,len,strlen,lenstr
      external lenstr 
C
      integer lpp,iform
      COMMON/CPPOUT/LPP,IFORM
C
      DATA LINIT/.TRUE./
C
      SAVE LINIT
c slmod begin
c...                  Need a tag for poloidal refinement
c                     only:

      IF (nbr.GT.0.OR.grdnmod.NE.0.OR.eirgrid.EQ.1) THEN
        IF (sloutput) WRITE(0,*) 'EXITING FROM PPOUT, FIX REQUIRED'
        RETURN
      ENDIF
c slmod end
C
C-----------------------------------------------------------------------
C FIRST CALL ONLY
C-----------------------------------------------------------------------
C
      IF(LINIT)THEN
C
        LPP   = TRANUNIT

c
C  SET IFORM = 0      FOR BINARY    OUTPUT IN TRANSFER FILE
C            = 1      FOR HEX       OUTPUT IN TRANSFER FILE
C            = 2      FOR FORMATTED OUTPUT IN TRANSFER FILE
C  --------------------------------------------------------
C
c
c  SET value for IFORM - BINARY is supposed to be the standard now  
c                        for all systems.
c
c        IFORM = 1
c
        IFORM = 0
c
        LINIT = .FALSE.
C
C  WRITE TRANSFER FILE NAME DEFINITION INFO
C  ----------------------------------------
C

        CALL NAMTRN(JOB,jfcb,TCODE,TSHOT,TDATE,tmachid,
     >              tmachine_opt,divshotid)
        CALL PPCSTR( '*' )
        CALL PPCSTR( '*CATALOGUE IDENTIFIER' )
        CALL PPCSTR( ' ' // TCODE  // ' ' // TMACHID
     &                   // TSHOT // ' ' // TDATE   )
c
        write(6,*) 'SAVING TRAN FILE:'
        write(6,*) 'TCODE  :',tcode
        write(6,*) 'TMACHID:',tmachid
        write(6,*) 'TSHOT  :',tshot
        write(6,*) 'TDATE  :',tdate
c
c        CALL PPCSTR( '*' )
c        CALL PPCSTR( '*CATALOGUE IDENTIFIER')
c        CALL PPCSTR( ' '//TCODE//' '//TSHOT//' '//TDATE )
c
C
C  WRITE VERSION INFO
C  ------------------
C
        CALL PPCSTR( '*' )
        CALL PPCSTR( '*CODE VERSION' )
        CALL PPCSTR( ' VERSON :'//VERSON )
C
C  WRITE USER INFO
C  ---------------
C
        CALL DMGUID(SYSUID,PREFIX)
c
C
        CALL PPCSTR( '*' )
        CALL PPCSTR( '*USER INFORMATION' )
        CALL PPCSTR( ' MACHINE:'//'JAC' )
        CALL PPCSTR( ' USER   :'//SYSUID )
C
C  PREPARE COMMENT SECTION
C  -----------------------
C
        CALL PPCSTR( '*' )
        CALL PPCSTR( '*RUN COMMENTS' )
        CALL PPCSTR( JOB )
c
c       Write case description into the TRAN file 
c       - maximum line length is 133
c
        len = lenstr(desc)
        do in = 1,len,133 
           strlen = min(len-in,132)
           call ppcstr(desc(in:in+strlen))
        end do

C
C  WRITE CURRENT CONTROL DATA
C  --------------------------
C
        CALL PPCSTR( '*' )
        CALL PPCSTR( '*CONTROL DATA' )
        REWIND(5)
 100    CBUFF = ' '
        READ(5,'(A80)',END=102) CBUFF
        CALL PPCSTR( CBUFF )
        GOTO 100
 102    REWIND(5)
        CALL PPCSTR( '*' )
        CALL PPCSTR( '*CONTROL DATA (FORMATTED)' )
C
        L1 = INDEX(USERIDH,'/')
        L2 = LENSTR(USERIDH)
        L1 = INDEX(USERIDH(L1+1:L2),'/') + L1
        L2 = INDEX(USERIDH(L1+1:L2),'/') + L1
        TMPUID = USERIDH(L1+1:L2-1)
c       CALL MAPUID(TMPUID,LUID)
        CALL PPCSTR( ' ADASIDH:'//TMPUID )
        CBUFF = ' YEARH  :XX'
        WRITE(CBUFF(10:11),'(I2.2)') IYEARH
        CALL PPCSTR( CBUFF )
C
        L1 = INDEX(USERIDZ,'/')
        L2 = LENSTR(USERIDZ)
        L1 = INDEX(USERIDZ(L1+1:L2),'/') + L1
        L2 = INDEX(USERIDZ(L1+1:L2),'/') + L1
        TMPUID = USERIDZ(L1+1:L2-1)
c       CALL MAPUID(TMPUID,LUID)
        CALL PPCSTR( ' ADASIDZ:'//TMPUID )
        CBUFF = ' YEARZ  :XX'
        WRITE(CBUFF(10:11),'(I2.2)') IYEARZ
        CALL PPCSTR( CBUFF )

c
C
C  WRITE GEOMETRY DATA
C  -------------------
C
        CALL PPCSTR( '*' )
        CALL PPCSTR( '*GEOMETRY' )
        CALL PPCSTR( ' '//EQUIL )
C
        CALL SETGEO(NPL,IOPENL,NXWL,NCL,NROWL,JPRGTL,JPLFTL,
     &              NSEPX,RSEPX,ZSEPX)
        CALL PPCSTR( 'NP,IOPEN,NXW,NC,NROW,JPRGT,JPLFT,IRFCT' )
        IF (IFORM.EQ.0) THEN
          WRITE(LPP) NPL,IOPENL,NXWL,NCL,NROWL,JPRGTL,JPLFTL,REFCT
        ELSE
          WRITE(CBUFF,'(8I4)') NPL,IOPENL,NXWL,NCL,NROWL,JPRGTL,JPLFTL,
     &                         REFCT
          CALL PPCSTR( CBUFF )
        ENDIF
C
        O8 = 1.D0
        CALL PPRECW('RPX   ','X-POINT R-COORD ','M ',RXP,O8,'R4',1,0)
        CALL PPRECW('ZPX   ','X-POINT Z-COORD ','M ',ZXP,O8,'R4',1,0)
        CALL PPRECW('R0    ','AXIAL R-COORD   ','M ',R0 ,O8,'R4',1,0)
        CALL PPRECW('Z0    ','AXIAL Z-COORD   ','M ',Z0 ,O8,'R4',1,0)
        CALL PPRECW('B0    ','AXIAL B-TOROIDAL','T ',CBPHI,O8,'R4',1,0)
        CALL PPRECW('RSEPX ','R FOR SEPARATRIX',
     &              'M  ',RSEPX(1),O8,'R4',NSEPX,0)
        CALL PPRECW('ZSEPX ','Z FOR SEPARATRIX',
     &              'M  ',ZSEPX(1),O8,'R4',NSEPX,0)
        CALL ICONV('ITAGDV','          ','     ',TAGDV(1,1)  , 1,NPL,0)
        CALL RTARG('RMESH ','R-COORD   ','M    ',
     &             RS(1,1)     ,RP(1)       ,O8,NPL,NDS,0)
        CALL RTARG('ZMESH ','Z-COORD   ','M    ',
     &             ZS(1,1)     ,ZP(1)       ,O8,NPL,NDS,0)
        CALL RCONV('RHO   ','RHO       ','     ',RHOG(1,1)   ,O8,NPL,0)
        CALL RCONV('THETA ','THETA     ','     ',THETAG(1,1) ,O8,NPL,0)
        CALL RCONV('HRHO  ','HRHO      ','     ',HRO(1,1)    ,O8,NPL,0)
        CALL RCONV('HTETA ','HTETA     ','     ',HTETA(1,1)  ,O8,NPL,0)
        CALL RCONV('BFI   ','B-TOROIDAL','T    ',BTS(1,1)    ,O8,NPL,0)
        CALL RCONV('PSI   ','PSI       ','     ',PSIFL(1,1)  ,O8,NPL,0)
C
        DO IR = 1,NRS
          DO IK = 1,NKS(IR)
            TMPMAT(IK,IR) = 1.0/KBFS(IK,IR)
          ENDDO
        ENDDO
        CALL RCONV('SH    ','BPOL/BTOT ','     ',TMPMAT(1,1) ,O8,NPL,0)
C
        CALL PPRECW('RVESM1','R-VESSEL 1','M  ',RVESM(1,1),O8,'R4',
     &               NVESM,0)
        CALL PPRECW('ZVESM1','Z-VESSEL 1','M  ',ZVESM(1,1),O8,'R4',
     &               NVESM,0)
        CALL PPRECW('RVESM2','R-VESSEL 2','M  ',RVESM(1,2),O8,'R4',
     &               NVESM,0)
        CALL PPRECW('ZVESM2','Z-VESSEL 2','M  ',ZVESM(1,2),O8,'R4',
     &               NVESM,0)
C
        IF( NVESP.GT.0 ) THEN
          CALL PPRECW('RPUMP1','R-PUMP   1','M  ',RVESM(NVESM+1,1),
     &                 O8,'R4',NVESP,0)
          CALL PPRECW('ZPUMP1','Z-PUMP   1','M  ',ZVESM(NVESM+1,1),
     &                 O8,'R4',NVESP,0)
          CALL PPRECW('RPUMP2','R-PUMP   2','M  ',RVESM(NVESM+1,2),
     &                 O8,'R4',NVESP,0)
          CALL PPRECW('ZPUMP2','Z-PUMP   2','M  ',ZVESM(NVESM+1,2),
     &                 O8,'R4',NVESP,0)
        END IF
C
        CALL ICONV ('KORPG ',' ',' ',KORPG , 1,NPL,0)
        CALL PPRECW('NVERTP',' ',' ',NVERTP, 1,'I4',  NPOLYP,0)
        CALL PPRECW('RVERTP',' ',' ',RVERTP,O8,'R4',5*NPOLYP,0)
        CALL PPRECW('ZVERTP',' ',' ',ZVERTP,O8,'R4',5*NPOLYP,0)
C
C  WRITE HEADER FOR TIME STEP DATA (ACTUALLY ITERATION DATA IN DIVIMP)
C  -------------------------------------------------------------------
C
        CALL PPCSTR( '*' )
        CALL PPCSTR( '*TIME STEP DATA' )
C
      ENDIF

C
C-----------------------------------------------------------------------
C RESULTS FOR CURRENT ITERATION
C-----------------------------------------------------------------------
C
      CBUFF = ' '
      TIME = 0.0

      WRITE(CBUFF,'(A,E11.5,A,I8)') '-TIME  :',TIME,' STEP  :',ITER
      CALL PPCSTR( CBUFF )
C
C MAIN PLASMA - BASIC QUANTITIES
C
      CALL PPRECW('HCH','ION ATOMIC NUMBER',' ',CIZB,1,'I4',1,0)
      CALL PPRECW('HMASS','ION MASS','AMU',CRMB,O8,'R4',1,0)
      CALL RTARG('DEN   ','ION DENSITY     ','M-3   ',
     &           KNBS(1,1)   ,KNDS(1)     ,O8,NPL,NDS,0)
      SF = RIZB
      CALL RTARG('DENEL ','EL. DENSITY','M-3   ',
     &           KNBS(1,1)   ,KNDS(1)     ,SF,NPL,NDS,0)
C
C LOAD SOUND SPEED INTO TARGET PARALLEL VELOCITY FOR CALCULATION
C OF ION SATURATION CURRENT
C
      DO I = 1,NDS
        TMPVEC(I) = 9.79E3 * SQRT (0.5*(KTEDS(I)+KTIDS(I))*
     &                       (1.0+RIZB)/CRMB)
        TMPVEC(I) = QTIM * TMPVEC(I)
      ENDDO
      SF = 1.0 / QTIM
      CALL RTARG('VTE   ','ION PAR.VEL     ','M/S   ',
     &           KVHS(1,1)   ,TMPVEC(1)   ,SF,NPL,NDS,0)
      CALL RTARG('TEV   ','ION TEMPERATURE ','EV    ',
     &           KTIBS(1,1)  ,KTIDS(1)    ,O8,NPL,NDS,0)
      CALL RTARG('TEVE  ','EL. TEMPERATURE ','EV    ',
     &           KTEBS(1,1)  ,KTEDS(1)    ,O8,NPL,NDS,0)
      SF = 1.0 / (QTIM*QTIM*EMI/CRMI)
      CALL RCONV('EPAR  ','PARA.ELEC.FIELD ',
     &           'V/M   ',KES(1,1)     ,SF,NPL,0)
C
C NEUTRAL HYDROGEN QUANTITIES
C
      CALL RCONV('DA','NEUT. ATOM DENSITY',
     &           'M-3',PINATOM(1,1),O8,NPL,0)
      CALL RCONV('DM','NEUT. MOL. DENSITY',
     &           'M-3',PINMOL(1,1),O8,NPL,0)
      CALL RCONV('DZ','NEUT. IMP. DENSITY',
     &           'M-3',PINZ0(1,1),O8,NPL,0)
      CALL RCONV('ENEUTA','NEUT. ATOM ENERGY',
     &           'EV',PINENA(1,1),O8,NPL,0)
      CALL RCONV('ENEUTM','NEUT. MOL. ENERGY',
     &           'EV',PINENM(1,1),O8,NPL,0)
      CALL RCONV('ENEUTZ','NEUT. IMP. ENERGY',
     &           'EV',PINENZ(1,1),O8,NPL,0)
      CALL RCONV('SOUN','IONISATION SOURCE',
     &           'M-3',PINION(1,1),O8,NPL,0)
      CALL RCONV('SOUZ','IMP. IONIS. SOURCE',
     &           'M-3',PINIONZ(1,1),O8,NPL,0)
      CALL RCONV('RIF','ION MOMENTUM SOURCE',
     &           'N M-3',PINMP(1,1),O8,NPL,0)
      CALL RCONV('SQI','ION ENERGY SOURCE',
     &           'W M-3',PINQI(1,1),O8,NPL,0)
      CALL RCONV('SQE','ELEC. ENERGY SOURCE',
     &           'W M-3',PINQE(1,1),O8,NPL,0)
C
C IMPURITIES
C
      CALL PPRECW('ZCH','IMPURITY ATOMIC NUMBER',' ',CION,1,'I4',1,0)
      CALL PPRECW('ZMASS','IMPURITY MASS','AMU',CRMI,O8,'R4',1,0)
      CALL PPRECW('NZ ','NUMBER OF CHARGE STATES FOLLOWED',
     &            ' ',NIZS,1,'I4',1,0)
C
C  FOR RECYCLING IMPURITIES ABSFAC=0 --> SET TO UNITY
C
      IF (ABSFAC.GT.0.0) THEN
        SF = ABSFAC
      ELSE
        SF = 1.0
      ENDIF
C
      DO 200 IZ=0,NIZS
        CBUF6 = ' '
        WRITE(CBUF6,'(A,I2.2)') 'IMP ',IZ
        CBUF8 = ' '
        WRITE(CBUF8,'(A,I2.2)') 'DENZ',IZ
        CALL RCONV(CBUF8,CBUF6//' DENSITY  ',
     &             'M-3' ,SDLIMS(1,1,IZ) ,SF,NPL,0)
 200  CONTINUE
C
      DO 210 IZ=0,NIZS
        CBUF6 = ' '
        WRITE(CBUF6,'(A,I2.2)') 'IMP ',IZ
        CBUF8 = ' '
        WRITE(CBUF8,'(A,I2.2)') 'TEVZ',IZ
        CALL RCONV(CBUF8,CBUF6//' TEMPERATURE',
     &             'EV' ,SDTS(1,1,IZ) ,O8,NPL,0)
 210  CONTINUE
C
      CALL RCONV('ZEFF','Z EFFECTIVE',' ',ZEFFS(1,1,3),O8,NPL,0)

c
C
C  WRITE HEADER FOR NIMBUS WALL DATA
C  ---------------------------------
C
        CALL PPCSTR( '*' )
        CALL PPCSTR( '*NIMBUS WALL' )
C
        CBUFF = ' '
        TIME = 0.0
        WRITE(CBUFF,'(A,E11.5,A,I8)') '-TIME :',TIME,' STEP :',ITER
        CALL PPCSTR( CBUFF )
C
        CALL SETMAR(NVESM,JVESM(1)
     &             ,IOTS,IOTE,IODS,IODE,IIDS,IIDE,IITS,IITE,iprs,ipre)
        CALL PPCSTR( '%10I IOTS,IOTE,IODS,IODE,IIDS,IIDE,'
     &                 //'IITS,IITE,IPRS,IPRE')
        IF(IFORM.EQ.0)THEN
          WRITE(LPP)  IOTS,IOTE,IODS,IODE,IIDS,IIDE,IITS,IITE,iprs,ipre
        ELSE
          CBUFF = ' '
          WRITE(CBUFF,'(10I4)')  IOTS,IOTE,IODS,IODE,IIDS,IIDE,
     >                    IITS,IITE,iprs,ipre
          CALL PPCSTR( CBUFF )
        ENDIF

C
C BASIC QUANTITIES
C ----------------
C
        CALL PPRECW('FLXHW1 ','NEUT. WALL FLUX',
     &              'M-2 S-1',FLUXHW,O8,'R4',NVESM,0)
        CALL PPRECW('FLXHW2 ','ION+ATOM WALL FLUX',
     &              'M-2 S-1',FLXHW2,O8,'R4',NVESM,0)
        CALL PPRECW('FLXHW3 ','IMP. SPUT. FLUX',
     &              'M-2 S-1',FLXHW3,O8,'R4',NVESM,0)
        CALL PPRECW('FLXHW4 ','IMP. REDEP FLUX',
     &              'M-2 S-1',FLXHW4,O8,'R4',NVESM,0)
        CALL PPRECW('FLXHW5 ','AVE. ATOM ENERGY',
     &              'EV',FLXHW5,O8,'R4',NVESM,0)
C
C-----------------------------------------------------------------------
C
      CALL PPCSTR( '*' )
      CALL PPCSTR( '*EOF' )
C
 9999 RETURN
C
      END
**++EOF
**==NAMTRN
C
C=======================================================================
      SUBROUTINE NAMTRN(JOB,jfcb,TCODE,TSHOT,TDATE,tmachid,tmachine_opt,
     >                  divshotid)
      IMPLICIT NONE
C
C***********************************************************************
C
C VERSION : V3.R1.M0
C
C ROUTINE : NAMTRN
C
C PURPOSE : PARSE THE PIECES NECESSARY TO CONSTRUCT THE ARCHIVE NAME
C           OF THE TRANSFER FILE
C
C INPUT   : (C *) JOB    - JOB INFORMATION FROM DIVIMP, INCLUDING
C                          INPUT FILE NAME AND DATE/TIME OF RUN
c           (I 4) TMACHINE_OPT - INTEGER VALUE DEFINING MACHINE TYPE
c           (C *) DIVSHOTID - SHOTID/# READ FROM THE INPUT FILE
C
C OUTPUT  : (C *) TCODE  - CODE FIELD OF ARCHIVED FILE NAME
C           (C *) TSHOT  - SHOT FIELD OF ARCHIVED FILE NAME
C           (C *) TDATE  - DATE FIELD OF ARCHIVED FILE NAME
c           (C *) TMACHID- ID STRING FOR MACHINE 
C
C HISTORY : V1.R1.M0 --- 11/04/94 --- CREATION
C         : V2.R1.M0 --- 08/11/95 --- VERSION FOR UNIX
C         : V3.R1.M0 --- 11/08/99 --- GET CASE NAME FROM FILE NAME
C                                     ATTACHED TO UNIT 5
C
C***********************************************************************
C
      CHARACTER*(*) JOB,jfcb,TCODE,TSHOT,TDATE,tmachid,divshotid
      integer tmachine_opt
C
      INTEGER       I,J,K,LJOB,LJFCB,LENSTR,IERR,len
      integer       startpos,endpos
c      CHARACTER*176 JFCB
      CHARACTER*1   POINT,OPEN,CLOSE,SLASH
      CHARACTER*4   EXT
      CHARACTER*3   MONTH(12)
      DATA          POINT/'.'/, OPEN/'('/, CLOSE/')'/, SLASH/'/'/
      DATA          EXT/'.d6i'/
      DATA          MONTH/'JAN','FEB','MAR','APR','MAY','JUN',
     >                    'JUL','AUG','SEP','OCT','NOV','DEC'/
c begin afmod
c for testing changes...
c      LOGICAL       ippchange
c      ippchange = .false.
c this section of code does not seem to accesss the ippchange common block value
c therefore have to define ippchange to test changes.

c end afmod      
C
C
C GET LENGTH OF JOB STRING
C ------------------------
C
      LJOB = LENSTR(JOB)
C
C HARDWIRE CODE TO BE DIVIMP
C --------------------------
C
      TCODE = 'DIVIMP'
C
C SET MACHINE TYPE
C --------------------------
C
      if (tmachine_opt.eq.0) then 
         tmachid = 'jet'
      elseif (tmachine_opt.eq.1) then 
         tmachid = 'diiid' 
      elseif (tmachine_opt.eq.2) then 
         tmachid = 'alcator'
      elseif (tmachine_opt.eq.3) then 
         tmachid = 'aug'
      elseif (tmachine_opt.eq.4) then 
         tmachid = 'iter'
      elseif (tmachine_opt.eq.5) then 
         tmachid = 'ignitor'
      elseif (tmachine_opt.eq.6) then 
         tmachid = 'fire'
      endif

C
c Done in RUNDIV and passed along
c
C GET NAME OF FILE ATTACHED TO UNIT 5
C -----------------------------------
C
c      WRITE(JFCB,'(176X)')
c      CALL IONAME(5,JFCB,99,IERR)
c      CALL IONAME(JFCB,IERR)
c
c      IF (IERR.NE.0) THEN
c        WRITE(6,*) ' NO FILE CONNECTED TO UNIT 5?!'
c        STOP
c      ENDIF
c
      LJFCB = LENSTR(JFCB)
C
C PARSE FILE NAME TO GET CASE NAME (LONG UNIX PATH NAMES
C HAVE BEEN RESULTING IN TRUNCATION IN JOB)
C ------------------------------------------------------
C
c      DO I = 1, LJFCB
c       IF (JFCB(I:I+3).EQ.EXT) GOTO 10
c      ENDDO
c
c
c     Use the case name to define the shotid if divshotid is blank
c
      if (divshotid.eq.' ') then
c
         endpos = 0 
c
c        Find end position of the case name
c
         DO I = 1, LJFCB
           IF (JFCB(I:I+3).EQ.EXT) then
              endpos = i-1
              exit
           endif 
         ENDDO
c
         if (endpos.gt.0) then 
c afmod begin
c        Krieger IPP/05 - startpos must be set to avoid crash in
c        absence of slash
c            IF (ippchange) THEN
c	       startpos=1
c            else 
c            endif   
c afmod end
c
c           jdemod - I don't see any reason to not initialize startpos
c     
            startpos = 1
c
            DO J = I-1, 1, -1
               IF (JFCB(J:J).EQ.SLASH) then 
                  startpos = j+1
                  exit
               endif
            ENDDO
c
c           Truncate SHOT identifier from case name to 8 characters
c     
c           Allow for up to 10 characters for TSHOT
c
            endpos = min(endpos,startpos+9)

c slmod begin
c
c     jdemod - this code should no longer be necessary with startpos 
c              initialized to 1 above
c
c... Intel compiler complaining about startpos=0, when string is declared
c    char(1:).  Not sure this is a good fix, but I am in a hurry, so
c    if there are problems check here first:
c            WRITE(0,*) 
c            WRITE(0,*) '---------------------------------------------- '
c            WRITE(0,*) ' SLOPPY FIX IN NAMTRN, CHECK HERE IF PROBLEMS '
c            WRITE(0,*) '---------------------------------------------- '
c            WRITE(0,*) 
c  
c           
c
c            startpos = MAX(1,startpos)
           
            TSHOT = jfcb(startpos:endpos)
c
c            TSHOT = jfcb(startpos:endpos)
c slmod end
c
c        Assign a default value if the case name is not available
c
         else

            TSHOT = 'NOSHOTNUM' 
  
         endif

         TSHOT = JFCB(J+1:I-1)

      else

         TSHOT = divshotid

      endif
c
c     Remove STOP error if name too long  
c
c   20 IF (I-J .GT. 10) GOTO 910
c

C
C SCAN THROUGH JOB TO FIND RUN DATE, WHICH IS DELIMITED BY SLASHES
C (THE FIRST 36 CHARACTERS OF JOB ARE RESERVED FOR THE INPUT FILE
C  NAME - START WITH CHARACTER 37!)
C ----------------------------------------------------------------
C
      TDATE = 'MMMDDYY'
      DO K = 37, LJOB
        IF (JOB(K:K).EQ.SLASH) GOTO 30
      ENDDO
      GOTO 920
C  NOTE THAT DATE HAS THE AMERICAN FORMAT!!!  (MM/DD/YY)
   30 TDATE(4:5) = JOB(K+1:K+2)
      TDATE(6:7) = JOB(K+4:K+5)
      READ(JOB(K-2:K-1),'(I2)') J
      TDATE(1:3) = MONTH(J)
C
c     Translate all output strings to lower case
c
      call chrutl(tcode,tcode)
      call chrutl(tmachid,tmachid)
      call chrutl(tshot,tshot)
      call chrutl(tdate,tdate)
c
      RETURN
C
  900 WRITE(6,*) 'ERROR FINDING CASE NAME IN JOB STRING!!!'
      STOP
  910 WRITE(6,*) 'CASE NAME TOO LONG (>8 CHARACTERS)!!!'
      STOP
  920 WRITE(6,*) 'ERROR FINDING DATA IN JOB STRING!!!'
      STOP
      END
**++EOF
**==RCONV
C
C=======================================================================
      SUBROUTINE RCONV(NAME,DESC,UNITS,RDATA,RSF,NDATA,ISTAG)
      IMPLICIT NONE
C
C***********************************************************************
C
C VERSION : V1.R2.M0
C
C ROUTINE : RCONV
C
C PURPOSE : CONVERT REAL DATA FROM DIVIMP MATRIX FORMAT TO K SYSTEM
C           AND WRITE TO TRANSFER FILE
C
C INPUT   : (C *) NAME   - DATA VARIABLE NAME
C           (C *) DESC   - DATA DESCRIPTION
C           (C *) UNITS  - DATA UNITS
C           (R*4) RDATA  - DATA
C           (R*4) RSF    - SCALE FACTOR
C           (I*4) NDATA  - NUMBER OF DATA VALUES
C           (I*4) ISTAG  - STAGGERED VARIABLE FLAG (0-UNSTAGGERED
C                                                   1-STAGGERED)
C
C OUTPUT  : NONE
C
C HISTORY : V1.R1.M0 --- 22/02/94 --- CREATION
C           V1.R2.M0 --- 31/01/95 --- ADD STAGGERED VARIABLE FLAG
C
C***********************************************************************
C
      include 'params'
      include 'cgeom'
      include 'comtor'
C
      CHARACTER*(*) NAME,DESC,UNITS
      REAL*4        RDATA(MAXNKS,MAXNRS),RSF
      INTEGER       NDATA,ISTAG
C
      REAL*4        RDATN(MAXNKS*MAXNRS)
      INTEGER       IK,IR,NPL
C
C
C CHECK INPUT
C -----------
C
      IF(NDATA.GT.MAXNKS*MAXNRS)THEN
        WRITE(6,*)
        WRITE(6,*) '*** ERROR(RCONV) : LOCAL ARRAY TOO SMALL'
        CALL EXIT(6)
      ENDIF
C
C CONVERT FROM (IK,IR)
C --------------------
C
c        write(6,'(a,50i5)') 'RCONV:',nrs,(nks(ir),ir=1,nrs)
c
        CALL RZERO(RDATN,MAXNKS*MAXNRS)
        NPL = 0
        DO 100 IR=1,NRS
          DO 110 IK=1,NKS(IR)
            if (kory(ir,ik).gt.0) then 
               RDATN(KORY(IR,IK)) = RDATA(IK,IR)
            endif
            NPL = NPL + 1
c
c            write(6,'(a,3i5,2(1x,e12.5))') name,kory(ir,ik),ir,ik,
c     >                 rdatn(kory(ir,ik)),rdata(ik,ir)
c
 110      CONTINUE
c
          IF (IR.GE.IRSEP .AND. .NOT.VIRTGRID) NPL = NPL + 2
c
c         Add code to fill the cells at the ends of the rings with 
c         non-zero values if that is necessary. 
c
c          IF (IR.GE.IRSEP .AND. .NOT.VIRTGRID) then 
c
c             NPL = NPL + 2
c
c            Duplicate ring end data
c
c             rdatn(kory(ir,1)-1) = rdatn(kory(ir,1))
c             rdatn(kory(ir,nks(ir))+1) = rdatn(kory(ir,nks(ir)))
c
c          endif 
c
c

 100    CONTINUE
C
C
C CONSISTENCY CHECK
C -----------------
C
      IF(NPL.NE.NDATA)THEN
        WRITE(6,*)
        WRITE(6,*) '*** ERROR(RCONV) : INCONSISTENT NUMBER OF POINTS'
        CALL EXIT(6)
      ENDIF
C
C
C PRINT RECORD
C ------------
C
      CALL PPRECW(NAME,DESC,UNITS,RDATN,RSF,'R4',NDATA,ISTAG)
C
C
      RETURN
      END
**++EOF
**==ICONV
C
C=======================================================================
      SUBROUTINE ICONV(NAME,DESC,UNITS,IDATA,ISF,NDATA,ISTAG)
      IMPLICIT NONE
C
C***********************************************************************
C
C VERSION : V1.R2.M0
C
C ROUTINE : ICONV
C
C PURPOSE : CONVERT INTEGER DATA FROM DIVIMP MATRIX FORMAT TO K SYSTEM
C           AND WRITE TO TRANSFER FILE
C
C INPUT   : (C *) NAME   - DATA VARIABLE NAME
C           (C *) DESC   - DATA DESCRIPTION
C           (C *) UNITS  - DATA UNITS
C           (I*4) IDATA  - DATA
C           (I*4) ISF    - SCALE FACTOR
C           (I*4) NDATA  - NUMBER OF DATA VALUES
C           (I*4) ISTAG  - STAGGERED VARIABLE FLAG (0-UNSTAGGERED
C                                                   1-STAGGERED)
C
C OUTPUT  : NONE
C
C HISTORY : V1.R1.M0 --- 22/02/94 --- CREATION
C           V1.R2.M0 --- 31/01/95 --- ADD STAGGERED VARIABLE FLAG
C
C***********************************************************************
C
      include 'params'
      include 'cgeom'
      include 'comtor'
C
      CHARACTER*(*) NAME,DESC,UNITS
      INTEGER       IDATA(MAXNKS,MAXNRS),ISF
      INTEGER       NDATA,ISTAG
C
      INTEGER       IDATN(MAXNKS*MAXNRS)
      INTEGER       IK,IR,NPL
C
C
C CHECK INPUT
C -----------
C
      IF(NDATA.GT.MAXNKS*MAXNRS)THEN
        WRITE(6,*)
        WRITE(6,*) '*** ERROR(ICONV) : LOCAL ARRAY TOO SMALL'
        CALL EXIT(6)
      ENDIF
C
C CONVERT FROM (IK,IR)
C --------------------
C
        CALL IZERO(IDATN,MAXNKS*MAXNRS)
        NPL = 0
        DO 100 IR=1,NRS
          DO 110 IK=1,NKS(IR)
            if (kory(ir,ik).gt.0) then  
               IDATN(KORY(IR,IK)) = IDATA(IK,IR)
            endif
            NPL = NPL + 1
 110      CONTINUE
c
          IF (IR.GE.IRSEP .AND. .NOT.VIRTGRID) NPL = NPL + 2
c
c         Add code to fill the cells at the ends of the rings with 
c         non-zero values if that is necessary. 
c
c          IF (IR.GE.IRSEP .AND. .NOT.VIRTGRID) then 
c
c             NPL = NPL + 2
c
c            Duplicate ring end data
c
c             idatn(kory(ir,1)-1) = idatn(kory(ir,1))
c             idatn(kory(ir,nks(ir))+1) = idatn(kory(ir,nks(ir)))
c
c          endif 
c
 100    CONTINUE
C
C
C CONSISTENCY CHECK
C -----------------
C
      IF(NPL.NE.NDATA)THEN
        WRITE(6,*)
        WRITE(6,*) '*** ERROR(ICONV) : INCONSISTENT NUMBER OF POINTS'
        CALL EXIT(6)
      ENDIF
C
C
C PRINT RECORD
C ------------
C
      CALL PPRECW(NAME,DESC,UNITS,IDATN,ISF,'I4',NDATA,ISTAG)
C
C
      RETURN
      END
**++EOF
**==RTARG
C
C=======================================================================
      SUBROUTINE RTARG(NAME,DESC,UNITS,RDATA,RTDATA,RSF,NDATA,NTDATA,
     >                 ISTAG)
      IMPLICIT NONE
C
C***********************************************************************
C
C VERSION : V1.R2.M0
C
C ROUTINE : RTARG
C
C PURPOSE : CONVERT REAL DATA FROM DIVIMP MATRIX FORMAT TO K SYSTEM,
C           APPEND TARGET DATA, AND WRITE TO TRANSFER FILE
C
C INPUT   : (C *) NAME   - DATA VARIABLE NAME
C           (C *) DESC   - DATA DESCRIPTION
C           (C *) UNITS  - DATA UNITS
C           (R*4) RDATA  - DATA
C           (R*4) RTDATA - TARGET DATA
C           (R*4) RSF    - SCALE FACTOR
C           (I*4) NDATA  - NUMBER OF DATA VALUES
C           (I*4) NTDATA - NUMBER OF TARGET DATA VALUES
C           (I*4) ISTAG  - STAGGERED VARIABLE FLAG (0-UNSTAGGERED
C                                                   1-STAGGERED)
C
C OUTPUT  : NONE
C
C HISTORY : V1.R1.M0 --- 25/02/94 --- CREATION
C           V1.R2.M0 --- 31/01/95 --- ADD STAGGERED VARIABLE FLAG
C
C***********************************************************************
C
      include 'params'
      include 'cgeom'
      include 'comtor'
C
      CHARACTER*(*) NAME,DESC,UNITS
      REAL*4        RDATA(MAXNKS,MAXNRS),RTDATA(MAXNDS),RSF
      INTEGER       NDATA,NTDATA,ISTAG
C
      REAL*4        RDATN(MAXNKS*MAXNRS+MAXNDS)
      INTEGER       IK,IR,ID,NPL
C
C
C CHECK INPUT
C -----------
C
      IF(NDATA+NTDATA.GT.MAXNKS*MAXNRS+MAXNDS)THEN
        WRITE(6,*)
        WRITE(6,*) '*** ERROR(RTARG) : LOCAL ARRAY TOO SMALL'
        CALL EXIT(6)
      ENDIF
C
C CONVERT FROM (IK,IR)
C --------------------
C
c        write(6,'(a,50i5)') 'RTARG:',nrs,(nks(ir),ir=1,nrs)
c
        CALL RZERO(RDATN,MAXNKS*MAXNRS+MAXNDS)
        NPL = 0
        DO 100 IR=1,NRS
          DO 110 IK=1,NKS(IR)
            if (kory(ir,ik).gt.0) then 
               RDATN(KORY(IR,IK)) = RDATA(IK,IR)
            endif
            NPL = NPL + 1
c
c            write(6,'(a,3i5,2(1x,e12.5))') name,kory(ir,ik),ir,ik,
c     >                 rdatn(kory(ir,ik)),rdata(ik,ir)
c

 110      CONTINUE
c
          IF (IR.GE.IRSEP .AND. .NOT.VIRTGRID) NPL = NPL + 2
c
c         Add code to fill the cells at the ends of the rings with 
c         non-zero values if that is necessary. 
c
c          IF (IR.GE.IRSEP .AND. .NOT.VIRTGRID) then 
c
c             NPL = NPL + 2
c
c            For Rtarg - write target data at the ends of the rings 
c
c             rdatn(kory(ir,1)-1) = rtdata(idds(ir,2))
c             rdatn(kory(ir,nks(ir))+1) = rtdata(idds(ir,1))
c
c          endif 
c
 100    CONTINUE
C
C APPEND TARGET DATA
C ------------------
C
        DO 120 ID=1,NDS
          NPL = NPL + 1
          RDATN(NPL) = RTDATA(ID)
 120    CONTINUE
C
C
C CONSISTENCY CHECK
C -----------------
C
      IF(NPL.NE.NDATA+NTDATA)THEN
        WRITE(6,*)
        WRITE(6,*) '*** ERROR(RTARG) : INCONSISTENT NUMBER OF POINTS'
        CALL EXIT(6)
      ENDIF
C
C
C PRINT RECORD
C ------------
C
      CALL PPRECW(NAME,DESC,UNITS,RDATN,RSF,'R4',NPL,ISTAG)
C
C
      RETURN
      END
**++EOF
**==SETGEO
C
C=======================================================================
      SUBROUTINE SETGEO(NPL,IOPENL,NXWL,NCL,NROW,JPRGTL,JPLFTL,
     &                  NSEPX,RSEPX,ZSEPX)
      implicit none
C
C***********************************************************************
C
C VERSION : V1.R1.M1
C
C ROUTINE : SETGEO
C
C PURPOSE : TO SET UP GEOMETRY GRID PARAMETERS FOR POST PROCESSOR FILE
C
C OUTPUT  : (I*4) NPL    - TOTAL NUMBER OF K VALUES
C           (I*4) IOPENL - FIRST OPEN RING
C           (I*4) NXWL   - WALL RING
C           (I*4) NCL    - NUMBER OF RINGS
C           (I*4) NROW   - NUMBER OF ROWS ON SOL RINGS (INCLUDING
C                          VIRTUAL POINTS)
C           (I*4) JPRGTL - FIRST ROW INTERSECTING MAIN PLASMA
C           (I*4) JPLFTL - LAST  ROW INTERCEPTING MAIN PLASMA
C           (I*4) NSEPX  - NUMBER OF POINTS ALONG SEPARATRIX
C           (R*4) RSEPX  - VECTOR OF R-VALUES OF SEPARATRIX POINTS
C           (R*4) ZSEPX  - VECTOR OF Z-VALUES OF SEPARATRIX POINTS
C
C HISTORY : V1.R1.M0 --- 22/02/94 --- CREATED FROM EDGE2D EQUIVALENT
C           V1.R1.M1 --- 25/02/94 --- ADDED SEPARATRIX CALCULATION
C
C***********************************************************************
C
      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
c
c     Declare the variables passed to SETGEO
c
      integer npl, iopenl, nxwl,ncl,nrow,jprgtl,jplftl
      integer nsepx
c
c      real rsepx(maxnks+1),zsepx(maxnks+1)     
C
      REAL*4    RSEPX(MAXNKS+1),ZSEPX(MAXNKS+1)
c
c     Local Variables
c
      integer nmp,nprv,nsol,np,ir,ik,k
C
C
C DERIVE REQUIRED OUTPUT
C ----------------------
C
      NROW   = NKS(IRSEP)
      IOPENL = IRSEP
      NXWL   = IRWALL
      NCL    = NRS
      JPRGTL = IKTO + 1
      JPLFTL = IKTI - 1
C
C
C DIVIMP IS OFTEN RUN WITHOUT VIRTUAL POINTS AT THE END OF OPEN RINGS
C (VIRTGRID=.FALSE.).  IF THIS IS THE CASE, PAD WITH ZEROS
C -------------------------------------------------------------------
C
      IF (.NOT.VIRTGRID) THEN
        NROW   = NROW + 2
        JPRGTL = JPRGTL + 1
        JPLFTL = JPLFTL + 1
      ENDIF
C
C CALCULATE THE TOTAL NUMBER OF POINTS, INCLUDING VIRTUAL POINTS
C --------------------------------------------------------------
C
      if (sonnet_grid_sub_type.eq.1) then 
         NMP  = (IOPENL-1)*(JPLFTL-JPRGTL+1)
      else
         NMP  = (IOPENL-1)*(JPLFTL-JPRGTL+2)
      endif
      NPRV = (JPRGTL-1+NROW-JPLFTL)*(NCL-NXWL)
      NSOL = NROW*(NXWL-IOPENL+1)
      NPL  = NMP+NPRV+NSOL
C
C CHECK THIS IS CONSISTENT WITH THE NUMBER ADDED RING BY RING
C -----------------------------------------------------------
C
      NP = 0
      DO IR = 1,NRS
        NP = NP + NKS(IR)
        IF (.NOT.VIRTGRID .AND. IR.GE.IOPENL) NP = NP + 2
      ENDDO
C
      IF( NPL.NE.NP )THEN
        WRITE(6,*)
        WRITE(6,*) '*** ERROR(SETGEO) : UNEXPECTED NP'
        WRITE(6,*) '    NP,IOPEN,NXW,NC,NROW,JPRGT,JPLFT'
        WRITE(6,'(10(1x,i8))') NP,IOPENL,NXWL,NCL,NROW,JPRGTL,JPLFTL
        WRITE(6,'(10(1x,i8))') NPL,NPRV,NSOL,NMP
c
        WRITE(0,*) '*** ERROR(SETGEO) : UNEXPECTED NP'
        WRITE(0,*) '    NP,IOPEN,NXW,NC,NROW,JPRGT,JPLFT'
        WRITE(0,'(10(1x,i8))') NP,IOPENL,NXWL,NCL,NROW,JPRGTL,JPLFTL
        WRITE(0,'(10(1x,i8))') NPL,NPRV,NSOL,NMP
        CALL EXIT(6)
      ENDIF
C
C BUILD A VECTOR OF SEPARATRIX POINTS FROM THE RELEVANT POLYGON VERTICES
C ----------------------------------------------------------------------
C
      IR = IRSEP
      NSEPX = 1
      DO IK = 1, NKS(IR)
        K = KORPG(IK,IR)
        IF (K.NE.0) THEN
          NSEPX = NSEPX + 1
          RSEPX(NSEPX) = RVERTP(1,K)
          ZSEPX(NSEPX) = ZVERTP(1,K)
        ENDIF
      ENDDO
      K = KORPG(NKS(IR),IR)
      IF (K.EQ.0) K = KORPG(NKS(IR)-1,IR)
      NSEPX = NSEPX + 1
      RSEPX(NSEPX) = RVERTP(4,K)
      ZSEPX(NSEPX) = ZVERTP(4,K)
C
C ADD EXTRAPOLATED END POINTS TO ENSURE INTERSECTION WITH TARGET ROWS
C -------------------------------------------------------------------
C
      RSEPX(1) = RSEPX(2) - (RSEPX(3)-RSEPX(2))
      ZSEPX(1) = ZSEPX(2) - (ZSEPX(3)-ZSEPX(2))
CPRNT WRITE(6,*) 'R 1,2,3'
C     WRITE(6,*) RSEPX(1),RSEPX(2),RSEPX(3)
C     WRITE(6,*) 'Z 1,2,3'
CPRNT WRITE(6,*) ZSEPX(1),ZSEPX(2),ZSEPX(3)
      NSEPX = NSEPX + 1
      RSEPX(NSEPX) = RSEPX(NSEPX-1) - (RSEPX(NSEPX-2)-RSEPX(NSEPX-1))
      ZSEPX(NSEPX) = ZSEPX(NSEPX-1) - (ZSEPX(NSEPX-2)-ZSEPX(NSEPX-1))
CPRNT WRITE(6,*) 'R N,N-1,N-2'
C     WRITE(6,*) RSEPX(NSEPX),RSEPX(NSEPX-1),RSEPX(NSEPX-2)
C     WRITE(6,*) 'Z N,N-1,N-2'
CPRNT WRITE(6,*) ZSEPX(NSEPX),ZSEPX(NSEPX-1),ZSEPX(NSEPX-2)
C
C
      RETURN
      END
**++EOF
**==RSUM
C
C=======================================================================
      SUBROUTINE RSUM(DATA,NIZS,TOTAL)
C
C***********************************************************************
C
C VERSION : V1.R1.M0
C
C ROUTINE : RSUM
C
C PURPOSE : TO SUM A VECTOR OF PROFILES OVER IONISATION STATE
C
C INPUT   : (R*4) DATA   - THE MATRIX OF PROFILES
C           (I*4) NIZS   - THE NUMBER OF AVAILABLE IONISATION STATES
C
C OUTPUT  : (R*4) TOTAL  - THE SUM OVER IZ=0,NIZS
C
C HISTORY : V1.R1.M0 --- 24/02/94 --- CREATED
C
C***********************************************************************
C
      include 'params'
      include 'comtor'
      include 'cgeom'
c
      real data,total 
      integer nizs,ik,ir,iz
C
      DIMENSION DATA(MAXNKS,MAXNRS,NIZS), TOTAL(MAXNKS,MAXNRS)
C
      CALL RZERO(TOTAL,MAXNKS*MAXNRS)
      DO IZ = 0, NIZS
        DO IR = 1,NRS
          DO IK = 1,NKS(IR)
            TOTAL(IK,IR) = TOTAL(IK,IR) + DATA(IK,IR,IZ)
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
**++EOF
**==SETMAR
C
C=======================================================================
      SUBROUTINE SETMAR( NSEG,JSEG
     >                 , IOTS,IOTE,IODS,IODE,IIDS,IIDE,IITS,IITE
     >                 , iprs,ipre )
      IMPLICIT REAL*8(A-H,O-Z)
C
C***********************************************************************
C
C VERSION : V2.R1.M0
C
C ROUTINE : SETMAR
C
C PURPOSE : TO SET UP SEGMENT MARKING PARAMETERS FOR POST PROCESSOR FILE
C
C INPUT   : (I*4) NSEG   - NUMBER OF SEGMENTS
C           (I*4) JSEG() - 1 = 'OT' ... OUTER TARGET
C                          2 = 'OC' ... OUTER CORNER
C                          3 = 'OD' ... OUTER DIVERTOR
C                          7 = 'MS' ... MAIN SOL
C                          6 = 'ID' ... INNER DIVERTOR
C                          5 = 'IC' ... INNER CORNER
C                          4 = 'IT' ... INNER TARGET
C                          8 = 'PV' ... PRIVATE VOID
C                          9 = 'BA' ... BAFFLE (NOT ALLOWED IN DIVIMP)
C                         10 = 'CB' ... COMPOUND BAFFLE (")
C
C OUTPUT  : (I*4) IOTS   - FIRST SEGMENT ON OUTER TARGET
C           (I*4) IOTE   - LAST  SEGMENT ON OUTER TARGET
C           (I*4) IODS   - FIRST SEGMENT ON OUTER DIVERTOR
C           (I*4) IODE   - LAST  SEGMENT ON OUTER DIVERTOR
C           (I*4) IIDS   - FIRST SEGMENT ON INNER DIVERTOR
C           (I*4) IIDE   - LAST  SEGMENT ON INNER DIVERTOR
C           (I*4) IITS   - FIRST SEGMENT ON INNER TARGET
C           (I*4) IITE   - LAST  SEGMENT ON INNER TARGET
C
C NOTE    : ORDER IS CLOCK-WISE RANGING STARTING AT SEGMENT # 1
C
C HISTORY : V1.R1.M0 --- 06/04/94 --- CREATION
C         : V2.R1.M0 --- 09/11/95 --- VERSION FOR DIVIMP WITH INTEGER
C                                     TAGS INSTEAD OF CHARACTER
C
C***********************************************************************
C
      integer nseg 
      INTEGER JSEG(NSEG)
      integer iots,iote,iods,iode,iids,iide,iits,iite,iprs,ipre,i
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
      iprs = 0
      ipre = 0
C
C
C DERIVE REQUIRED OUTPUT
C ----------------------
C
      DO 100 I = 1 , NSEG-1
         IF( JSEG(I).NE.1 .AND. JSEG(I+1).EQ.1 ) IOTS = I+1
         IF( JSEG(I).EQ.1 .AND. JSEG(I+1).NE.1) IOTE = I
         IF( JSEG(I).NE.3 .AND. JSEG(I+1).EQ.3 ) IODS = I+1
         IF( JSEG(I).EQ.3 .AND. JSEG(I+1).NE.3 ) IODE = I
         IF( JSEG(I).NE.6 .AND. JSEG(I+1).EQ.6 ) IIDS = I+1
         IF( JSEG(I).EQ.6 .AND. JSEG(I+1).NE.6 ) IIDE = I
         IF( JSEG(I).NE.4 .AND. JSEG(I+1).EQ.4 ) IITS = I+1
         IF( JSEG(I).EQ.4 .AND. JSEG(I+1).NE.4 ) IITE = I
         IF( JSEG(I).NE.8 .AND. JSEG(I+1).EQ.8 ) IPRS = I+1
         IF( JSEG(I).EQ.8 .AND. JSEG(I+1).NE.8 ) IPRE = I
  100 CONTINUE
C
      RETURN
      END
**++EOF
**==MAPUID
C
C=======================================================================
c      SUBROUTINE MAPUID(SYSUID,LUID)
C
C TRIES TO MAP USER ID FROM JAC TO JET3090.  IF NO MATCHING ID IS FOUND
C THEN THE UID IS LEFT UNCHANGED.
C
c      INTEGER       LUID
c      CHARACTER*(*) SYSUID
C
c      INTEGER       INIT, MAXUID, NUIDS
c      PARAMETER     (MAXUID=100)
c      CHARACTER     UNIXID(MAXUID)*8, IBMID(MAXUID)*6
c      DATA INIT/0/
C
c      IF (SYSUID(1:3).EQ.'JET') RETURN
C
c      IF (INIT.EQ.0) THEN
c        READ(LUID,*) NUIDS
c        IF (NUIDS.GT.MAXUID) THEN
c          WRITE(6,*) ' TOO MANY UIDS IN TRANSCRIPTION FILE!!!'
c          STOP
c        ENDIF
c        DO I = 1, NUIDS
c          READ(LUID,*) IBMID(I), UNIXID(I)
c        ENDDO
c        INIT = 1
c      ENDIF
C
c      DO I = 1, NUIDS
c        IF (SYSUID.EQ.UNIXID(I)) THEN
c          SYSUID = IBMID(I)
c          RETURN
c        ENDIF
c      ENDDO
C
c      RETURN
c      END
c**++EOF
