c     -*Fortran*-
c
      SUBROUTINE SOLEDGE(irlim1,irlim2,ikopt)
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_comhr
      use mod_slcom
      IMPLICIT  NONE
      integer   irlim1,irlim2,ikopt
C     INCLUDE   "PARAMS"
c     include 'params'
C     INCLUDE   "CGEOM"
c     include 'cgeom'
C     INCLUDE   "COMTOR"
c     include 'comtor'
c     include 'comhr'
c slmod begin
c     INCLUDE 'slcom'
c slmod end
C
C
C *****************************************************************
C *                                                               *
C * SOLEDGE:  THIS ROUTINE SETS THE TEMPERATURE AND DENSITY OF    *
C *           SOL PLASMA- IT ALSO ASSIGNS THE VELOCITY AND        *
C *           ELECTRIC FIELD IN THE SOL. THIS IS DONE USING A SET *
C *           OF SIMULTANEOUS EQUATIONS WHICH FORM A SOMEWHAT     *
C *           CONSISTENT MODEL OF THE PLASMA EDGE.                *
C *                                                               *
C * DAVID ELDER     OCT  1991                                     *
C *                                                               *
C *****************************************************************
C
      DOUBLE PRECISION ARG1,ARG2,ARG3,ARG4
      double precision n,v
      INTEGER IR,J,IK,NS,PLATEOPT,IRLIMIT
      INTEGER IKMID,ikstart,ikend,ikfirst,iklast
      DOUBLE PRECISION DS ,V0,V0I,PINF,PINFI
      REAL*8 ACT_PRESS
      DOUBLE PRECISION RCF,RCFI
      DOUBLE PRECISION CIS1,CIS2
      INTEGER RINGNO
      EXTERNAL CIS1,CIS2,RINGNO
      DOUBLE PRECISION DS1,DS2,DP1,DP2,DT1,DT2,NB1,NB2
      DOUBLE PRECISION S,SMAX,NBP,NBPI,TEBP,TEBPI,TIBP,TIBPI
      double precision spredi,spredo,sprev,predls,sinj
      DOUBLE PRECISION LPPA,LPPAI,LPPAEB,LPPAEI,LPPAIB,LPPAII
      DOUBLE PRECISION MASSI,TMP
      DOUBLE PRECISION DELTAS,pmax
      DOUBLE PRECISION FSRC,LMSRC,LNSRC
      DOUBLE PRECISION SRCION,SRCRAD, srcpei
      DOUBLE PRECISION LSSIZ,LPSIZ,MFACT,MFACT2
c
      double precision sl1,sl2,slv,v1,v1i,t1,t1i,t2,t2i,n1,n1i,
     >                 lprad,lpradi,news,ti1,ti1i
      double precision sl1i,sl2i,slvi
      double precision sl1a,sl1ai,sl1b,sl1bi  
c
      EXTERNAL SRCION,SRCRAD,srcpei
      CHARACTER*80 INFO
      REAL SOLVTEMP
      EXTERNAL SOLVTEMP
c
c     Variables for SOL 16+
c
      double precision rootn,gamman,sact,strt,stmp
      integer ikn,in
      double precision soltelast,soltecur,soltilast,solticur
      double precision solnelast,solnecur,solvellast,solvelcur
      double precision solcorlast,solcorcur
      double precision solprn,solprhn,solpcxn,solphn,solpein
      double precision helpi,soli
C
      double precision  l1r,l1ri,l2r,l2ri,lvr,lvri,ter,
     >                  teri,tir,tiri,nr,nri,vbm,vbmi,qr,qri,
     >                  n_exp,n_expi
c
c     Extra parameters - have no effect unless specified in 
c     for S21            additional per ring data sets.
c
c     l1a ... allow for more detail in the description of the 
c              density evolution in region A - three linear 
c              fitted regions. 
c
      double precision l1a,l1ai,l1b,l1bi,nr1a,nr1ai,nr1b,nr1bi,
     >                         ter1a,ter1ai,tir1a,tir1ai,
     >                         ter1b,ter1bi,tir1b,tir1bi
C
       DO 600 IR = IRlim1, IRLIM2
c
c       SET IK values for inner loops        
c
        call set_ikvals(ir,ikstart,ikend,ikopt)
c
        SMAX = KSMAXS(IR)
        pmax = kpmaxs(ir)
c
c     jdemod - do not read the starting target conditions from the grid
c            - they should be loaded from the specified target conditions
c              for each ring
c
c     
        tebp = kteds(idds(ir,2))
        tibp = ktids(idds(ir,2))
        nbp  = knds(idds(ir,2))
        v0   = -abs(kvds(idds(ir,2)))
c
        tebpi = kteds(idds(ir,1))
        tibpi = ktids(idds(ir,1))
        nbpi  = knds(idds(ir,1))
        v0i   = -abs(kvds(idds(ir,1)))
c
c        NBP = KNBS(1,IR)
c        NBPI = KNBS(NKS(IR),IR) 
C
C       SET UP ELECTRON AND ION PLATE TEMPERATURES
C
c        TEBP  = KTEBS(1,IR)
c        TEBPI = KTEBS(NKS(IR),IR)
c        TIBP  = KTIBS(1,IR)
c        TIBPI = KTIBS(NKS(IR),IR)
C
C       IF (TEBP.LT.1.0E+00) TEBP = 1.0E+00
C       IF (TIBP.LT.1.0E+00) TIBP = 1.0E+00
C       IF (TEBPI.LT.1.0E+00) TEBPI = 1.0E+00
C       IF (TIBPI.LT.1.0E+00) TIBPI = 1.0E+00
C
C
C       SAVE PLATE TEMPERATURES AND DENSITY FOR EACH LAUNCH POSITION
C
C       THE ARRAY IDDS GIVES THE NDS INDICES OF THE PLATE POINTS FOR
C       EACH RING. THE INNER PLATE IS 1 AND THE OUTER IS 2.
C
C
c        V0  = - SQRT(0.5*EMI*(TEBP+TIBP)*(1+RIZB)/CRMB)
c        V0I = - SQRT(0.5*EMI*(TEBPI+TIBPI)*(1+RIZB)/CRMB)
c
c        write (0,'(a,8g8.2)') 'SOLEDGE:',v0,v0i,tebpi,tibpi
C
c
c        
c     jdemod - these are set in initplasma and should not be       
c              overwritten here
c     
c     if (ikopt.eq.1.or.ikopt.eq.3) then 
c           KTEDS(IDDS(IR,2)) = TEBP
c           KTIDS(IDDS(IR,2)) = TIBP
c           KNDS(IDDS(IR,2)) = NBP
c           KVDS(idds(ir,2)) = V0
c        endif
c
c        if (ikopt.eq.2.or.ikopt.eq.3) then 
c           KTEDS(IDDS(IR,1)) = TEBPI
c           KTIDS(IDDS(IR,1)) = TIBPI
c           KNDS(IDDS(IR,1)) = NBPI
c           KVDS(idds(ir,1)) = -V0I
c        endif
C
C       Set power coefficients
C
        IF (CIOPTF.EQ.12.OR.CIOPTF.EQ.14.or.cioptf.eq.16
     >     .or.cioptf.eq.18.or.cioptf.eq.19) THEN
           LPPA = (2.0*TIBP+5.0*TEBP)*1.602192E-19*NBP*
     >             DABS(V0)
           LPPAI = (2.0*TIBPI+5.0*TEBPI)*1.602192E-19*NBPI*
     >             DABS(V0I)
        ENDIF
        IF (CIOPTF.EQ.13.OR.CIOPTF.EQ.15.or.cioptf.eq.17
     >      .or.cioptf.eq.20) THEN
           LPPAEB = 5.0*TEBP*1.602192E-19*NBP*
     >            DABS(V0)
           LPPAEI = 5.0*TEBPI*1.602192E-19*NBPI*
     >            DABS(V0I)
           LPPAIB = 2.0*TIBP*1.602192E-19*NBP*
     >            DABS(V0)
           LPPAII = 2.0*TIBPI*1.602192E-19*NBPI*
     >            DABS(V0I)
           LPPA = LPPAEB
           LPPAI= LPPAEI
        ENDIF
c
c       SOL 21 - load parameters
c
        if (cioptf.eq.21) then
c
           call load_s21params(ir,l1r,l1ri,l2r,l2ri,
     >                         lvr,lvri,ter,teri,tir,tiri,nr,nri,
     >                         vbm,vbmi,qr,qri,n_exp,n_expi,
     >                         l1a,nr1a,l1ai,nr1ai,
     >                         l1b,nr1b,l1bi,nr1bi,
     >                         ter1a,ter1ai,tir1a,tir1ai,
     >                         ter1b,ter1bi,tir1b,tir1bi)
c
c           write (6,'(a,i4,10(1x,g9.3))') 'S21:',ir,ter,teri,tir,tiri,
c     >                                    nr,nri,l1r,l1ri,l2r,l2ri
c
c
c          Length reference switch for SOL 21 - options
c
c             0 = relative to SMAX
c             1 = relative to PMAX
c             2 = absolute SMAX units  
c             3 = absolute PMAX units
c
           if (s21refsw.eq.0) then 
c
              sl1a = l1a * smax
              sl1b = l1b * smax
              sl1  = l1r * smax
              sl2  = l2r * smax
              slv  = lvr * smax
c
              sl1ai = l1ai * smax
              sl1bi = l1bi * smax
              sl1i  = l1ri * smax
              sl2i  = l2ri * smax
              slvi  = lvri * smax
c
           elseif (s21refsw.eq.1) then 
c
              sl1a = l1a * pmax
              sl1b = l1b * pmax
              sl1  = l1r * pmax
              sl2  = l2r * pmax
              slv  = lvr * pmax
c
              sl1ai = l1ai * pmax
              sl1bi = l1bi * pmax
              sl1i  = l1ri * pmax
              sl2i  = l2ri * pmax
              slvi  = lvri * pmax
c
           elseif (s21refsw.eq.2.or.s21refsw.eq.3) then 
c
              sl1a = l1a 
              sl1b = l1b 
              sl1  = l1r 
              sl2  = l2r 
              slv  = lvr 
c
              sl1ai = l1ai 
              sl1bi = l1bi 
              sl1i  = l1ri 
              sl2i  = l2ri 
              slvi  = lvri 
c
           endif    
c
c          Convert P-coordinates to S-coordinates since code operates in S
c
           if (s21refsw.eq.1.or.s21refsw.eq.3) then 
c
c             OUTER
c
              call cnvrtptos(sl1a,news,ir)
              sl1a = news 
c
              call cnvrtptos(sl1b,news,ir)
              sl1b = news 
c
              call cnvrtptos(sl1,news,ir)
              sl1 = news 
c
              call cnvrtptos(sl2,news,ir)
              sl2 = news 
c
              call cnvrtptos(slv,news,ir)
              slv = news 
c
c             INNER
c
              call cnvrtptos(sl1ai,news,ir)
              sl1ai = news 
c
              call cnvrtptos(sl1bi,news,ir)
              sl1bi = news 
c
              call cnvrtptos(sl1i,news,ir)
              sl1i = news 
c
              call cnvrtptos(sl2i,news,ir)
              sl2i = news 
c
              call cnvrtptos(slvi,news,ir)
              slvi = news 
c
           endif
c         
           lppa = (5.0+2.0*tibp/tebp+15.0/tebp)
     >             *nbp*dabs(v0)*tebp*1.602192e-19
           lppai= (5.0+2.0*tibpi/tebpi+15.0/tebpi)
     >             *nbpi*dabs(v0i)*tebpi*1.602192e-19
c
           lprad = qr * lppa / (sl2-sl1)
           lpradi= qri * lppai / (sl2i-sl1i)
c
           t1    = ter * tebp
           t1i   = teri * tebpi
c
           ti1   = tir * tibp 
           ti1i  = tiri * tibpi 
c
c          Pressure limited nratio option in SOL 21
c          is indicated by negative values of the 
c          stated ratio.
c
c          nrat applies to the first half of the ring and 
c          nrati to the second - as a result - they must check
c          for a defined upstream pressure on the opposite half
c          of the ring before pressure matching is invoked. 
c          If such a pressure is not defined the aboslute value
c          of nrat is used to calculate n1, n1i.
c
           ikmid = ikmids(ir)+1
c
           if (nr.lt.0.0) then 
c
c             Check for valid pressure
c
              if (kpress(ikmid,ir,2).gt.0.0) then 
c
c                Drop velocity term when estimating the 
c                required pressure multiplier? 
c
c                Do initially for simplicity
c
                 n1 = kpress(ikmid,ir,2)/(ech*(t1+ti1))
c
              else
c
c                 No pressure defined
c

                  n1 = abs(nr) * nbp

              endif
c
           else

               n1 = nr * nbp

           endif               
c
c          
c
           if (nri.lt.0.0) then 
c
c
c             Check for valid pressure
c
              if (kpress(ikmid-1,ir,2).gt.0.0) then 
c
c                Drop velocity term when estimating the 
c                required pressure multiplier? 
c
c                Do initially for simplicity
c
                 n1i = kpress(ikmid-1,ir,2)/(ech*(t1i+ti1i))
c
              else
c
c                 No pressure defined
c

                  n1i = abs(nri) * nbpi

              endif
c
           else

               n1i = nri * nbpi
c
           endif               
c
c          Store actual value of density rise ratio used
c
           nrat_used(ir,2) = n1/nbp
           nrat_used(ir,1) = n1i/nbpi
c
c
c           write (6,'(a,2i5,8(1x,g12.5))') 'KPRESS:',ir,ikmid,
c     >               kpress(ikmid,ir,2),kpress(ikmid-1,ir,2),
c     >               n1,nbp,nrat_used(ir,2),n1i,nbpi,nrat_used(ir,1)
c
c
           t2    = (t1**3.5+7.0/(2.0*ck0)*(lppa*(sl2-sl1)
     >                    +0.5*(sl2-sl1)**2*lprad))**(2.0/7.0)
           t2i   = (t1i**3.5+7.0/(2.0*ck0)*(lppai*(sl2i-sl1i)
     >                    +0.5*(sl2i-sl1i)**2*lpradi))**(2.0/7.0)
c
           v1    = (nbp * v0) / n1 * vbm
           v1i   = (nbpi * v0i) / n1i * vbmi
c
c           write (6,'(a,6(1x,g12.5))') 'SOL21o:',
c     >                              sl1,sl2,slv,smax,lppa,lprad
c           write (6,'(a,6(1x,g12.5))') 'SOL21i:',
c     >                           sl1i,sl2i,slvi,smax,lppai,lpradi
c           write (6,'(a,14(1x,g9.3))') 'SOL21b:',
c     >                         t1,t1i,n1,n1i,t2,t2i,v1,v1i,v0,v0i,
c     >                         vbm,vbmi,crmb
        endif
 
C
C
        ikmid = ikmids(ir)+1
C
C       Set up pressure value
C
        IF (CIOPTF.EQ.12.OR.CIOPTF.EQ.14.OR.CIOPTF.EQ.16
     >    .or.cioptf.eq.18.or.cioptf.eq.19) THEN
          PINF  = 4.0*NBP*ECH*TEBP
          PINFI = 4.0*NBPI*ECH*TEBPI
        ELSEIF (CIOPTF.EQ.13.OR.CIOPTF.EQ.15.or.cioptf.eq.17
     >      .or.cioptf.eq.20) THEN
          PINF  = NBP*ECH*(2.0*TEBP+2.0*TIBP)
          PINFI = NBPI*ECH*(2.0*TEBPI+2.0*TIBPI)
        ENDIF
C
        MASSI = CRMB * AMU
c
C     SET UP THE IONIZATION SOURCE ARRAY FIRST. THIS HAS MORE
C     POINTS THAN THE BACKGROUND N,V,T,E ARRAYS IN ORDER TO MAKE
C     THE NUMERICAL INTEGRATION FOR THE SOURCE TERMS MORE ACCURATE.
C
c       Execute the ionization source code only for
c       options that require it. i.e. NOT Sol 21
c
        if (cioptf.ne.21) then
c
C
C        WRITE(6,*) 'P,PI,MI:',PINF,PINFI,MASSI
C
C        WRITE(6,*) 'SMAX:',SMAX
C
        FSRC = 1.0
        IF (CSOPT.EQ.0) THEN
          LNSRC = CSOLLS
          IF (CSOLLS.GT.0.5) LNSRC = 0.5
          LMSRC = 0.0
        ELSEIF (CSOPT.EQ.1) THEN
          LNSRC = 0.5
          LMSRC = CSOLLS
        ELSEIF (CSOPT.EQ.4.OR.CSOPT.EQ.5) THEN
          LNSRC = CSOLLT
          LMSRC = CSOLLS
        ENDIF
c
c       Write out parameters for testing ...
c
c        write(6,*) smax,ir,nbp,tebp,v0,nbpi,tebpi,v0i
c
C
C       WILL REPLACE WITH SPECIFIED VALUES IF FOUND.
C
        IF (FLUXROPT.EQ.1)  CALL FLUXRLOOK(IR,FSRC,LMSRC,LNSRC)
C
C         WRITE(6,*) 'FLUX:',FSRC,LNSRC,LMSRC
C
C       CALL THE SETUP ROUTINE : IT ESTABLISHES ALL THE VALUES
C       REQUIRED FOR THE VARIOUS IONIZATION AND RADIATION OPTIONS
C       THAT UNDERLIE THE CALCULATION OF THE SOL CHARACTERISTICS.
C
        CALL SETUPVAL(SMAX,CSOPT,CPOPT,FSRC,LNSRC,LMSRC,
     >      NBP,V0,TEBP,NBPI,V0I,TEBPI,RCF,RCFI,LSSIZ,LPSIZ,IR,
     >      LPPA,LPPAI)
C
c
c     Set up predictor barrier limit
c
      spredo = 3.0 * csolls * smax
      if (spredo.gt.smax/2.0) spredo = smax/2.0
      spredi = spredo
c
c      write (6,*) 'sizes:',csolls,smax,spredi,spredo
c
c     Calculate an approximate average injection position
c     for use in calculating barrier value estimates.
c     The actual mean injection position is not known
c     until after the particles have been tracked.
c
      if ((ciopte.eq.2.or.ciopte.eq.3.or.ciopte.eq.5)
     >     .and.cneuta.eq.1) then
         sinj = (injf1+injf2)/2.0 * smax
      else
         sinj = 0.0
      endif
c
c     Endif to SOL opt NOT 21
c
      endif

C
C       OUTER PLATES ....
C

        if (ikopt.eq.1.or.ikopt.eq.3) then 
c
        sprev = 0.0
c
        DO 300 IK = ikstart, IKMID-1
C
          S  = KSS(IK,IR)
c
c     Calculate revised pressure
          if (sol13_pdist.gt.0.0) then 
             act_press = pinf *
     >         min((s/(sol13_pdist*smax)),1.0) *sol13_padd +
     >         pinf     
          else
             act_press = pinf * (1.0+sol13_padd)
          endif
             
          IF (CIOPTF.EQ.12) THEN
c
C           CALCULATE BACKGROUND TEMPERATURE
C
            KTEBS(IK,IR) = (TEBP**3.5 + 7.0/(2.0*CK0)*
     >       (LPPA*S+CIS2(S,SRCRAD,CPOPT,0)))**(2.0/7.0)
C
            KTIBS(IK,IR) = KTEBS(IK,IR)
C
C           CALCULATE DENSITY
C
            SOLI = CIS1(S,SRCION,CSOPT,0)
            GAMMAN = NBP*V0 + SOLI + RCF * S
C
            IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)
     >        GAMMAN = 0.0
C
            call calcnv(dble(ktebs(ik,ir)),dble(ktibs(ik,ir)),
     >            gamman,act_press,n,v)
c
            knbs(ik,ir) = n
            kvhs(ik,ir) = v
c
c            WRITE(6,*) 'NUMBERS:',IK,IR,KTEBS(IK,IR),KTIBS(IK,IR),
c     >             KNBS(IK,IR),KVHS(IK,IR)
c            WRITE(6,*) 'OTHERS1:',ROOTN,GAMMAN,S,act_press
c            WRITE(6,*) 'OTHERS2:',(PINF/(2*ECH*KTEBS(IK,IR)))**2,
c     >            -2.0*(MASSI*GAMMAN**2)/(ECH*KTEBS(IK,IR)),
c     >             NBP*V0,SOLI,RCF
C
c
c
c
          ELSEIF (CIOPTF.EQ.13.OR.CIOPTF.EQ.15) THEN
c
c
c
C
C           SET DIFFERENCE BETWEEN 13 AND 15
C
            IF (CIOPTF.EQ.13) THEN
               MFACT = 7.0/2.0
               MFACT2 = 1.0
            ELSEIF (CIOPTF.EQ.15) THEN
               MFACT = 7.0/4.0
               MFACT2 = S / (SMAX/2.0)
            ENDIF
C
C           CALCULATE ELECTRON TEMPERATURE
C
            TMP = CIS2(S,SRCRAD,CPOPT,0)
            KTEBS(IK,IR) = (TEBP**3.5 + (MFACT/CK0)*
     >       (LPPAEB*S*MFACT2+
     >       TMP))**(2.0/7.0)
C
c             WRITE (6,*) 'NUM1:',TMP,LPPAEB,TEBP,MFACT,
c     >                CK0,MFACT2
C
C           CALCULATE ION TEMPERATURE
C
            KTIBS(IK,IR) = (TIBP**3.5 + (MFACT/CK0I)*
     >       (LPPAIB*S*MFACT2))**(2.0/7.0)
C
C           CALCULATE DENSITY
C
            SOLI = CIS1(S,SRCION,CSOPT,0)
            GAMMAN = NBP*V0 + SOLI + RCF * S
 
            IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)
     >        GAMMAN = 0.0
c
            call calcnv(dble(ktebs(ik,ir)),dble(ktibs(ik,ir)),
     >            gamman,act_press,n,v)
c
            knbs(ik,ir) = n
            kvhs(ik,ir) = v
c
c            WRITE(6,*) 'NUMBERS:',IK,IR,KTEBS(IK,IR),KTIBS(IK,IR),
c     >             KNBS(IK,IR),KVHS(IK,IR),RCF
c            WRITE(6,*) 'OTHERS1:',ROOTN,GAMMAN,S
c            WRITE(6,*) 'OTHERS2:',(PINF/(ECH*(KTEBS(IK,IR)
c     >            +KTIBS(IK,IR)) ) ) **2,
c     >            -4.0*(MASSI*GAMMAN**2)/(ECH*(KTEBS(IK,IR)
c     >            +KTIBS(IK,IR))),
c     >             NBP*V0,SOLI,RCF*S
C
C
C
          ELSEIF (CIOPTF.EQ.14) THEN
C
C
c
            SOLI = CIS1(S,SRCION,CSOPT,0)
            GAMMAN = NBP*V0 + SOLI + RCF * S
C
            IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)
     >        GAMMAN = 0.0
C
C           CALCULATE THE BACKGROUND TEMPERATURE
C
            write(6,*) 'soli:',ir,ik,soli,gamman
            IF (IK.EQ.1.AND.S.LE.0.0) THEN
              KTEBS(IK,IR) = TEBP
            ELSEIF (IK.EQ.1.AND.S.GT.0.0) THEN
              ARG1 = KSS(IK,IR) * 5.0 * GAMMAN * ECH
              solprn  = CIS1(S,SRCRAD,CPOPT,0)
              ARG2 = KSS(IK,IR) * (LPPA+solprn)
              ARG3 = TEBP * CK0
              ARG4 = -CK0
C              write(6,*) 'args:',arg1,arg2,arg3,arg4,solprn,lppa
              KTEBS(IK,IR) = SOLVTEMP(ARG1,ARG2,ARG3,ARG4)
            ELSE
              DELTAS = KSS(IK,IR) - KSS(IK-1,IR)
              ARG1 = DELTAS * 5.0 * GAMMAN * ECH
              solprn = CIS1(S,SRCRAD,CPOPT,0)
              ARG2 = DELTAS * ( LPPA + solprn)
              ARG3 = KTEBS(IK-1,IR) * CK0
              ARG4 = -CK0
C              write(6,*) 'args:',arg1,arg2,arg3,arg4,solprn,lppa
              KTEBS(IK,IR) = SOLVTEMP(ARG1,ARG2,ARG3,ARG4)
            ENDIF
C
            WRITE (6,*) 'TEMPS:',IK,IR,KTEBS(IK,IR)
C
            KTIBS(IK,IR) = KTEBS(IK,IR)
C
C           CALCULATE DENSITY
C
            call calcnv(dble(ktebs(ik,ir)),dble(ktibs(ik,ir)),
     >            gamman,act_press,n,v)
c
            knbs(ik,ir) = n
            kvhs(ik,ir) = v
c
c
c
          elseif (cioptf.eq.16) then
c
c
c
c           This SOL option uses a slightly different calculational
c           method but should yield results comparable to SOL 14.
c
            if (ik.eq.1) then
              deltas = kss(ik,ir)/msolpt
              strt   = 0.0
            else
              deltas = (kss(ik,ir)-kss(ik-1,ir))/msolpt
              strt   = kss(ik-1,ir)
            endif
            write(6,*) 'deltas:',deltas,strt
            do 220 in = 1,msolpt
c
              stmp = strt + deltas * real(in)
              ikn = (ik-1)*msolpt + in
              solcorcur = stmp
c
              if (ikn.eq.1) then
                 soltelast = tebp
                 soltilast = tibp
                 solnelast = nbp
                 solvellast = v0
                 solcorlast = 0.0
                 if (ir.eq.cirhr) then
                   solte(0) = tebp
                   solti(0) = tibp
                   solne(0) = nbp
                   solvel(0) = v0
                   solcor(0) = 0.0
                 endif
              endif
c
              SOLI = CIS1(Stmp,SRCION,CSOPT,0)
              GAMMAN = NBP*V0 + SOLI + RCF * Stmp
C
              IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)
     >          GAMMAN = 0.0
C
C             CALCULATE THE BACKGROUND TEMPERATURE
C
              solprn =  cis1(stmp,srcrad,cpopt,0)
C
              soltecur =
     >            - ((lppa + solprn)
     >              *deltas
     >              + ck0*soltelast**3.5)
     >            / ( 5.0 * gamman * deltas * ech
     >              - ck0*soltelast**2.5)
C
              solticur = soltecur
C
C             CALCULATE DENSITY
C
              call calcnv(soltecur,solticur,
     >            gamman,act_press,solnecur,solvelcur)
c
c             Store hi-res background if on selected ring
c
              if (ir.eq.cirhr) then
                 solte(ikn) = soltecur
                 solti(ikn) = solticur
                 solne(ikn) = solnecur
                 solvel(ikn) = solvelcur
                 solcor(ikn) = solcorcur
              endif
c
c             Update old variables for next loop iteration
c
              soltelast = soltecur
              soltilast = solticur
              solnelast = solnecur
              solvellast = solvelcur
              solcorlast = solcorcur
c
c              write(6,*) 'sol16:t:',ikn,ir,gamman,rootn,
c     >            solte(ikn,ir),solti(ikn,ir),solne(ikn,ir)
c              write(6,*) 'cont:',solvel(ikn,ir),solcor(ikn,ir),
c     >                    tmp
c
c
 220        continue
c
c           Assign values to actual grid points
c
            ktebs(ik,ir) = soltecur
            ktibs(ik,ir) = solticur
            knbs(ik,ir) = solnecur
            kvhs(ik,ir) = solvelcur
c            write(6,*) 'sol16:',ik,ir,ktebs(ik,ir),ktibs(ik,ir),
c     >                 knbs(ik,ir),kvhs(ik,ir)
c
c           Calculate 1/2 up to mid-plane
c
            if (ik.eq.(ikmid-1)) then
              deltas = (kss(ik+1,ir)-kss(ik,ir))/msolpt
              strt   = kss(ik,ir)
              do 230 in = 1,msolpt/2
                stmp = strt + deltas * real(in)
                ikn = ik*msolpt + in
                solcorcur = stmp
c
                SOLI = CIS1(Stmp,SRCION,CSOPT,0)
                GAMMAN = NBP*V0 + SOLI + RCF * Stmp
C
                IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)
     >            GAMMAN = 0.0
C
C               CALCULATE THE BACKGROUND TEMPERATURE
C
                solprn = cis1(stmp,srcrad,cpopt,0)
                soltecur =
     >             - ((lppa + solprn)
     >               *deltas
     >               + ck0*soltelast**3.5)
     >             / ( 5.0 * gamman * deltas * ech
     >               - ck0*soltelast**2.5)
C
                solticur = soltecur
C
C               CALCULATE DENSITY
C
                call calcnv(soltecur,solticur,
     >              gamman,act_press,solnecur,solvelcur)
c
c               Store hi-res background if on selected ring
c
                if (ir.eq.cirhr) then
                   solte(ikn) = soltecur
                   solti(ikn) = solticur
                   solne(ikn) = solnecur
                   solvel(ikn) = solvelcur
                   solcor(ikn) = solcorcur
                endif
c
c               Update old variables for next loop iteration
c
                soltelast = soltecur
                soltilast = solticur
                solnelast = solnecur
                solvellast = solvelcur
                solcorlast = solcorcur
c
 230          continue
            endif
c
c
c
          elseif (cioptf.eq.17) then
c
c
c
c
c           This SOL option uses a slightly different calculational
c           method but should yield results comparable to SOL 14
c           except that Te is decoupled from Ti.
c
            if (ik.eq.1) then
              deltas = kss(ik,ir)/msolpt
              strt   = 0.0
            else
              deltas = (kss(ik,ir)-kss(ik-1,ir))/msolpt
              strt   = kss(ik-1,ir)
            endif
            write(6,*) 'deltas:',deltas,strt
            do 240 in = 1,msolpt
c
              stmp = strt + deltas * real(in)
              ikn = (ik-1)*msolpt + in
              solcorcur = stmp
c
              if (ikn.eq.1) then
                 soltelast = tebp
                 soltilast = tibp
                 solnelast = nbp
                 solvellast = v0
                 solcorlast = 0.0
                 if (ir.eq.cirhr) then
                   solte(0) = tebp
                   solti(0) = tibp
                   solne(0) = nbp
                   solvel(0) = v0
                   solcor(0) = 0.0
                 endif
              endif
c
c             Calculate Gamma
c
              SOLI = CIS1(Stmp,SRCION,CSOPT,0)
              GAMMAN = NBP*V0 + SOLI + RCF * Stmp
C
              IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)
     >          GAMMAN = 0.0
C
C             CALCULATE THE BACKGROUND TEMPERATURE
C
              solprn = cis1(stmp,srcrad,cpopt,0)
              soltecur =
     >            - ((lppaeb + solprn)
     >              *deltas
     >              + ck0*soltelast**3.5)
     >            / ( 2.5 * gamman * deltas * ech
     >              - ck0*soltelast**2.5)
C
              solticur =
     >            - (lppaib
     >              *deltas
     >              + ck0i*soltilast**3.5)
     >            / ( 2.5 * gamman * deltas * ech
     >              - ck0i*soltilast**2.5)
C
C             CALCULATE DENSITY
C
              call calcnv(soltecur,solticur,
     >              gamman,act_press,solnecur,solvelcur)
c
c             Store hi-res background if on selected ring
c
              if (ir.eq.cirhr) then
                 solte(ikn) = soltecur
                 solti(ikn) = solticur
                 solne(ikn) = solnecur
                 solvel(ikn) = solvelcur
                 solcor(ikn) = solcorcur
              endif
c
c             Update old variables for next loop iteration
c
              soltelast = soltecur
              soltilast = solticur
              solnelast = solnecur
              solvellast = solvelcur
              solcorlast = solcorcur
c
c              write(6,*) 'sol16:t:',ikn,ir,gamman,rootn,
c     >            solte(ikn,ir),solti(ikn,ir),solne(ikn,ir)
c              write(6,*) 'cont:',solvel(ikn,ir),solcor(ikn,ir),
c     >                    tmp
c
c
 240        continue
c
c           Assign values to actual grid points
c
            ktebs(ik,ir) = soltecur
            ktibs(ik,ir) = solticur
            knbs(ik,ir) = solnecur
            kvhs(ik,ir) = solvelcur
c            write(6,*) 'sol16:',ik,ir,ktebs(ik,ir),ktibs(ik,ir),
c     >                 knbs(ik,ir),kvhs(ik,ir)
c
c           Calculate 1/2 up to mid-plane
c
            if (ik.eq.(ikmid-1)) then
              deltas = (kss(ik+1,ir)-kss(ik,ir))/msolpt
              strt   = kss(ik,ir)
              do 245 in = 1,msolpt/2
                stmp = strt + deltas * real(in)
                ikn = ik*msolpt + in
                solcorcur = stmp
c
                SOLI = CIS1(Stmp,SRCION,CSOPT,0)
                GAMMAN = NBP*V0 + SOLI + RCF * Stmp
C
                IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)
     >            GAMMAN = 0.0
C
C               CALCULATE THE BACKGROUND TEMPERATURE
C
                solprn = cis1(stmp,srcrad,cpopt,0)
                soltecur =
     >             - ((lppaeb + solprn)
     >               *deltas
     >               + ck0*soltelast**3.5)
     >             / ( 2.5 * gamman * deltas * ech
     >               - ck0*soltelast**2.5)
C
              solticur =
     >            - (lppaib
     >              *deltas
     >              + ck0i*soltilast**3.5)
     >            / ( 2.5 * gamman * deltas * ech
     >              - ck0i*soltilast**2.5)
C
C               CALCULATE DENSITY
C
                call calcnv(soltecur,solticur,
     >             gamman,act_press,solnecur,solvelcur)
c
c               Store hi-res background if on selected ring
c
                if (ir.eq.cirhr) then
                   solte(ikn) = soltecur
                   solti(ikn) = solticur
                   solne(ikn) = solnecur
                   solvel(ikn) = solvelcur
                   solcor(ikn) = solcorcur
                endif
c
c               Update old variables for next loop iteration
c
                soltelast = soltecur
                soltilast = solticur
                solnelast = solnecur
                solvellast = solvelcur
                solcorlast = solcorcur
c
 245          continue
            endif
c
c
c
          elseif (cioptf.eq.18) then
c
c
c
c           Adds the kinetic term 1/2mv**2*gamma to the
c           solution of SOL16.
c
            if (ik.eq.1) then
              deltas = kss(ik,ir)/msolpt
              strt   = 0.0
            else
              deltas = (kss(ik,ir)-kss(ik-1,ir))/msolpt
              strt   = kss(ik-1,ir)
            endif
            write(6,*) 'deltas:',deltas,strt
            do 247 in = 1,msolpt
c
              stmp = strt + deltas * real(in)
              ikn = (ik-1)*msolpt + in
              solcorcur = stmp
c
              if (ikn.eq.1) then
                 soltelast = tebp
                 soltilast = tibp
                 solnelast = nbp
                 solvellast = v0
                 solcorlast = 0.0
                 if (ir.eq.cirhr) then
                   solte(0) = tebp
                   solti(0) = tibp
                   solne(0) = nbp
                   solvel(0) = v0
                   solcor(0) = 0.0
                 endif
              endif
c
              SOLI = CIS1(Stmp,SRCION,CSOPT,0)
              GAMMAN = NBP*V0 + SOLI + RCF * Stmp
C
              IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)
     >          GAMMAN = 0.0
C
C             CALCULATE THE BACKGROUND TEMPERATURE
C
              solprn = cis1(stmp,srcrad,cpopt,0)
              soltecur =
     >            - ((lppa + solprn
     >                + 0.5*massi*solvellast**2*gamman)
     >              *deltas
     >              + ck0*soltelast**3.5)
     >            / ( 5.0 * gamman * deltas * ech
     >              - ck0*soltelast**2.5)
C
              solticur = soltecur
C
C             CALCULATE DENSITY
C
              call calcnv(soltecur,solticur,
     >            gamman,act_press,solnecur,solvelcur)
c
c             Store hi-res background if on selected ring
c
              if (ir.eq.cirhr) then
                 solte(ikn) = soltecur
                 solti(ikn) = solticur
                 solne(ikn) = solnecur
                 solvel(ikn) = solvelcur
                 solcor(ikn) = solcorcur
              endif
c
c             Update old variables for next loop iteration
c
              soltelast = soltecur
              soltilast = solticur
              solnelast = solnecur
              solvellast = solvelcur
              solcorlast = solcorcur
c
c              write(6,*) 'sol18:t:',ikn,ir,gamman,rootn,
c     >            solte(ikn,ir),solti(ikn,ir),solne(ikn,ir)
c              write(6,*) 'cont:',solvel(ikn,ir),solcor(ikn,ir),
c     >                    tmp
c
c
 247        continue
c
c           Assign values to actual grid points
c
            ktebs(ik,ir) = soltecur
            ktibs(ik,ir) = solticur
            knbs(ik,ir) = solnecur
            kvhs(ik,ir) = solvelcur
c            write(6,*) 'sol16:',ik,ir,ktebs(ik,ir),ktibs(ik,ir),
c     >                 knbs(ik,ir),kvhs(ik,ir)
c
c           Calculate 1/2 up to mid-plane
c
            if (ik.eq.(ikmid-1)) then
              deltas = (kss(ik+1,ir)-kss(ik,ir))/msolpt
              strt   = kss(ik,ir)
              do 248 in = 1,msolpt/2
                stmp = strt + deltas * real(in)
                ikn = ik*msolpt + in
                solcorcur = stmp
c
                SOLI = CIS1(Stmp,SRCION,CSOPT,0)
                GAMMAN = NBP*V0 + SOLI + RCF * Stmp
C
                IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)
     >            GAMMAN = 0.0
C
C               CALCULATE THE BACKGROUND TEMPERATURE
C
                solprn = cis1(stmp,srcrad,cpopt,0)
                soltecur =
     >             - ((lppa + solprn
     >                + 0.5 *massi*solvellast**2*gamman)
     >               *deltas
     >               + ck0*soltelast**3.5)
     >             / ( 5.0 * gamman * deltas * ech
     >               - ck0*soltelast**2.5)
C
                solticur = soltecur
C
C               CALCULATE DENSITY
C
                call calcnv(soltecur,solticur,
     >              gamman,act_press,solnecur,solvelcur)
c
c               Store hi-res background if on selected ring
c
                if (ir.eq.cirhr) then
                   solte(ikn) = soltecur
                   solti(ikn) = solticur
                   solne(ikn) = solnecur
                   solvel(ikn) = solvelcur
                   solcor(ikn) = solcorcur
                endif
c
c               Update old variables for next loop iteration
c
                soltelast = soltecur
                soltilast = solticur
                solnelast = solnecur
                solvellast = solvelcur
                solcorlast = solcorcur
c
 248          continue
            endif
c
c
c
          elseif (cioptf.eq.19) then
c
c
c
c           Adds the kinetic term 1/2mv**2*gamma to the
c           solution of SOL16.
c
            if (ik.eq.1) then
              deltas = kss(ik,ir)/msolpt
              strt   = 0.0
            else
              deltas = (kss(ik,ir)-kss(ik-1,ir))/msolpt
              strt   = kss(ik-1,ir)
            endif
            write(6,*) 'deltas:',deltas,strt
            do 1000 in = 1,msolpt
c
              stmp = strt + deltas * real(in)
              ikn = (ik-1)*msolpt + in
              solcorcur = stmp
c
              if (ikn.eq.1) then
                 soltelast = tebp
                 soltilast = tibp
                 solnelast = nbp
                 solvellast = v0
                 solcorlast = 0.0
                 if (ir.eq.cirhr) then
                   solte(0) = tebp
                   solti(0) = tibp
                   solne(0) = nbp
                   solvel(0) = v0
                   solcor(0) = 0.0
                 endif
              endif
c
              SOLI = CIS1(Stmp,SRCION,CSOPT,0)
              GAMMAN = NBP*V0 + SOLI + RCF * Stmp
C
              IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)
     >          GAMMAN = 0.0
C
C             CALCULATE THE BACKGROUND TEMPERATURE
C
              helpi = 17.5 + (5.0+37.5/soltelast)
     >                * (1.0+0.25/soltelast)
     >                * log10(1.0d21/solnelast)
c
              solprn = cis1(stmp,srcrad,cpopt,0)
              soltecur =
     >            - ((lppa + solprn
     >                + 0.5*massi*solvellast**2*gamman
     >                + helpi * soli * 1.6d-19)
     >              *deltas
     >              + ck0*soltelast**3.5)
     >            / ( 5.0 * gamman * deltas * ech
     >              - ck0*soltelast**2.5)
C
              solticur = soltecur
C
C             CALCULATE DENSITY
C
              call calcnv(soltecur,solticur,
     >            gamman,act_press,solnecur,solvelcur)
c
c             Store hi-res background if on selected ring
c
              if (ir.eq.cirhr) then
                 solte(ikn) = soltecur
                 solti(ikn) = solticur
                 solne(ikn) = solnecur
                 solvel(ikn) = solvelcur
                 solcor(ikn) = solcorcur
              endif
c
c             Update old variables for next loop iteration
c
              soltelast = soltecur
              soltilast = solticur
              solnelast = solnecur
              solvellast = solvelcur
              solcorlast = solcorcur
c
c              write(6,*) 'sol18:t:',ikn,ir,gamman,rootn,
c     >            solte(ikn,ir),solti(ikn,ir),solne(ikn,ir)
c              write(6,*) 'cont:',solvel(ikn,ir),solcor(ikn,ir),
c     >                    tmp
c
c
 1000       continue
c
c           Assign values to actual grid points
c
            ktebs(ik,ir) = soltecur
            ktibs(ik,ir) = solticur
            knbs(ik,ir) = solnecur
            kvhs(ik,ir) = solvelcur
c            write(6,*) 'sol16:',ik,ir,ktebs(ik,ir),ktibs(ik,ir),
c     >                 knbs(ik,ir),kvhs(ik,ir)
c
c           Calculate 1/2 up to mid-plane
c
            if (ik.eq.(ikmid-1)) then
              deltas = (kss(ik+1,ir)-kss(ik,ir))/msolpt
              strt   = kss(ik,ir)
              do 1010 in = 1,msolpt/2
                stmp = strt + deltas * real(in)
                ikn = ik*msolpt + in
                solcorcur = stmp
c
                SOLI = CIS1(Stmp,SRCION,CSOPT,0)
                GAMMAN = NBP*V0 + SOLI + RCF * Stmp
C
                IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)
     >            GAMMAN = 0.0
C
C               CALCULATE THE BACKGROUND TEMPERATURE
C
                helpi = 17.5 + (5.0+37.5/soltelast)
     >                * (1.0+0.25/soltelast)
     >                * log10(1.0d21/solnelast)
c
                solprn = cis1(stmp,srcrad,cpopt,0)
                soltecur =
     >             - ((lppa + solprn
     >                + 0.5 *massi*solvellast**2*gamman
     >                + helpi * soli *1.6d-19)
     >               *deltas
     >               + ck0*soltelast**3.5)
     >             / ( 5.0 * gamman * deltas * ech
     >               - ck0*soltelast**2.5)
C
                solticur = soltecur
C
C               CALCULATE DENSITY
C
                call calcnv(soltecur,solticur,
     >              gamman,act_press,solnecur,solvelcur)
c
c               Store hi-res background if on selected ring
c
                if (ir.eq.cirhr) then
                   solte(ikn) = soltecur
                   solti(ikn) = solticur
                   solne(ikn) = solnecur
                   solvel(ikn) = solvelcur
                   solcor(ikn) = solcorcur
                endif
c
c               Update old variables for next loop iteration
c
                soltelast = soltecur
                soltilast = solticur
                solnelast = solnecur
                solvellast = solvelcur
                solcorlast = solcorcur
c
 1010         continue
            endif
c
c
c
          elseif (cioptf.eq.20) then
c
c
c
c           Adds the kinetic term 1/2mv**2*gamma to the
c           solution of SOL16.
c
            maxind = nks(ir)*msolpt+msolpt
            if (ik.eq.1) then
              deltas = kss(ik,ir)/msolpt
              strt   = 0.0
            else
              deltas = (kss(ik,ir)-kss(ik-1,ir))/msolpt
              strt   = kss(ik-1,ir)
            endif
            solpein = srcpei(0.0d0,0,0,deltas,1)
            write(6,*) 'deltas:',deltas,strt
            do 1020 in = 1,msolpt
c
              stmp = strt + deltas * real(in)
              ikn = (ik-1)*msolpt + in
              cind = ikn
              solcorcur = stmp
c
              if (ikn.eq.1) then
                 soltelast = tebp
                 soltilast = tibp
                 solnelast = nbp
                 solvellast = v0
                 solcorlast = 0.0
c
c                Hi-res background for calculating
c                Pei(s) and its integral.
c
                 soltehr(0) = tebp
                 soltihr(0) = tibp
                 solnehr(0) = nbp
                 solcorhr(0) = 0.0
c
                 if (ir.eq.cirhr) then
                   solte(0) = tebp
                   solti(0) = tibp
                   solne(0) = nbp
                   solvel(0) = v0
                   solcor(0) = 0.0
                 endif
              endif
c
              SOLI = CIS1(Stmp,SRCION,CSOPT,0)
              GAMMAN = NBP*V0 + SOLI + RCF * Stmp
C
              IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)
     >          GAMMAN = 0.0
C
C             CALCULATE THE BACKGROUND TEMPERATURE
C
              helpi = 17.5 + (5.0+37.5/soltelast)
     >                * (1.0+0.25/soltelast)
     >                * log10(1.0d21/solnelast)
c
              solprn = cis1(stmp,srcrad,cpopt,0)
              solpein = srcpei(stmp,ikn,0,deltas,0)
              soltecur =
     >            - ((lppaeb + solprn + solpein
     >                + 0.5*massi*solvellast**2*gamman
     >                + helpi * soli * 1.6d-19)
     >              *deltas
     >              + ck0*soltelast**3.5)
     >            / ( 2.5 * gamman * deltas * ech
     >              - ck0*soltelast**2.5)
              solticur =
     >            - ((lppaib -solpein)
     >              *deltas
     >              + ck0i*soltilast**3.5)
     >            / ( 2.5 * gamman * deltas * ech
     >              - ck0i*soltilast**2.5)
C
C             CALCULATE DENSITY
C
              call calcnv(soltecur,solticur,
     >            gamman,act_press,solnecur,solvelcur)
c
c             Store hi-res background if on selected ring
c
              soltehr(ikn) = soltecur
              soltihr(ikn) = solticur
              solnehr(ikn) = solnecur
              solcorhr(ikn) = solcorcur
              if (ir.eq.cirhr) then
                 solte(ikn) = soltecur
                 solti(ikn) = solticur
                 solne(ikn) = solnecur
                 solvel(ikn) = solvelcur
                 solcor(ikn) = solcorcur
              endif
c
c             Update old variables for next loop iteration
c
              soltelast = soltecur
              soltilast = solticur
              solnelast = solnecur
              solvellast = solvelcur
              solcorlast = solcorcur
c
c              write(6,*) 'sol18:t:',ikn,ir,gamman,rootn,
c     >            solte(ikn,ir),solti(ikn,ir),solne(ikn,ir)
c              write(6,*) 'cont:',solvel(ikn,ir),solcor(ikn,ir),
c     >                    tmp
c
c
 1020       continue
c
c           Assign values to actual grid points
c
            ktebs(ik,ir) = soltecur
            ktibs(ik,ir) = solticur
            knbs(ik,ir) = solnecur
            kvhs(ik,ir) = solvelcur
c            write(6,*) 'sol16:',ik,ir,ktebs(ik,ir),ktibs(ik,ir),
c     >                 knbs(ik,ir),kvhs(ik,ir)
c
c           Calculate 1/2 up to mid-plane
c
            if (ik.eq.(ikmid-1)) then
              deltas = (kss(ik+1,ir)-kss(ik,ir))/msolpt
              strt   = kss(ik,ir)
              do 1030 in = 1,msolpt/2
                stmp = strt + deltas * real(in)
                ikn = ik*msolpt + in
                cind = ikn
                solcorcur = stmp
c
                SOLI = CIS1(Stmp,SRCION,CSOPT,0)
                GAMMAN = NBP*V0 + SOLI + RCF * Stmp
C
                IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)
     >            GAMMAN = 0.0
C
C               CALCULATE THE BACKGROUND TEMPERATURE
C
                helpi = 17.5 + (5.0+37.5/soltelast)
     >                * (1.0+0.25/soltelast)
     >                * log10(1.0d21/solnelast)
c
                solprn = cis1(stmp,srcrad,cpopt,0)
                solpein = srcpei(stmp,ikn,0,deltas,0)
                soltecur =
     >            - ((lppaeb + solprn + solpein
     >                + 0.5*massi*solvellast**2*gamman
     >                + helpi * soli * 1.6d-19)
     >              *deltas
     >              + ck0*soltelast**3.5)
     >            / ( 2.5 * gamman * deltas * ech
     >              - ck0*soltelast**2.5)
                solticur =
     >            - ((lppaib -solpein)
     >              *deltas
     >              + ck0i*soltilast**3.5)
     >            / ( 2.5 * gamman * deltas * ech
     >              - ck0i*soltilast**2.5)
C
C               CALCULATE DENSITY
C
                call calcnv(soltecur,solticur,
     >              gamman,act_press,solnecur,solvelcur)
c
c               Store hi-res background if on selected ring
c
                soltehr(ikn) = soltecur
                soltihr(ikn) = solticur
                solnehr(ikn) = solnecur
                solcorhr(ikn) = solcorcur
                if (ir.eq.cirhr) then
                   solte(ikn) = soltecur
                   solti(ikn) = solticur
                   solne(ikn) = solnecur
                   solvel(ikn) = solvelcur
                   solcor(ikn) = solcorcur
                endif
c
c               Update old variables for next loop iteration
c
                soltelast = soltecur
                soltilast = solticur
                solnelast = solnecur
                solvellast = solvelcur
                solcorlast = solcorcur
c
 1030         continue
            endif
c
c         Sol option 21 ... model for a detached plasma.
c
          elseif (cioptf.eq.21) then
c
c           Test S-value and then calculate Te,Ti,N and V
c
            if (s.lt.sl1a) then 
c
              ktebs(ik,ir) = tebp + (tebp*ter1a-tebp) * (s/sl1a)
              ktibs(ik,ir) = tibp + (tibp*tir1a-tibp) * (s/sl1a)
c
              knbs(ik,ir)  = nbp + (nbp*nr1a-nbp) * (s/sl1a)**n_exp
c
              kvhs(ik,ir)  = (nbp * v0)/ knbs(ik,ir) *aux_vel21
c
            elseif (s.lt.sl1b) then 
c
              ktebs(ik,ir) = tebp*ter1a + (tebp*ter1b-tebp*ter1a)
     >                       * (s-sl1a)/(sl1b-sl1a)
              ktibs(ik,ir) = tibp*tir1a + (tibp*tir1b-tibp*tir1a)
     >                       * (s-sl1a)/(sl1b-sl1a)
c
              knbs(ik,ir)  = nbp*nr1a 
     >                     + (nbp*nr1b-nbp*nr1a) 
     >                       * ((s-sl1a)/(sl1b-sl1a))**n_exp
c
              kvhs(ik,ir)  = (nbp * v0)/ knbs(ik,ir) * aux_vel21
c
            elseif (s.le.sl1) then
c
              ktebs(ik,ir) = tebp*ter1b + (t1-tebp*ter1b) 
     >                       * (s-sl1b)/(sl1-sl1b)
              ktibs(ik,ir) = tibp*tir1b + (ti1-tibp*tir1b) 
     >                       * (s-sl1b)/(sl1-sl1b)
c
c              knbs(ik,ir)  = nbp + (n1-nbp) * (s/sl1)**n_exp
c
              knbs(ik,ir)  = nbp*nr1b 
     >                     + (n1-nbp*nr1b) 
     >                       * ((s-sl1b)/(sl1-sl1b))**n_exp
c
              kvhs(ik,ir)  = (nbp * v0)/ knbs(ik,ir)
c
            elseif (s.le.sl2) then
              ktebs(ik,ir) = (t1**3.5+7.0/(2.0*ck0)*(lppa*(s-sl1)
     >                    +0.5*(s-sl1)**2*lprad))**(2.0/7.0)
              ktibs(ik,ir) = ktebs(ik,ir)
              knbs(ik,ir)  = (n1*t1) / ktebs(ik,ir)
              if (s.le.slv) then
                kvhs(ik,ir)= v1 * (slv-s)/(slv-sl1)
              else
                kvhs(ik,ir)= 0.0
              endif
            else
              ktebs(ik,ir) = (t2**3.5+7.0/(2.0*ck0)*
     >                       ((1.0+qrat)*lppa*(s-sl2)))**(2.0/7.0)
              ktibs(ik,ir) = ktebs(ik,ir)
              knbs(ik,ir)  = (n1*t1) / ktebs(ik,ir)
              if (s.le.slv) then
                kvhs(ik,ir)= v1 * (slv-s)/(slv-sl1)
              else
                kvhs(ik,ir)= 0.0
              endif
            endif
          ENDIF
c
c
          IF (ABS(KVHS(IK,IR)).LT.1.0) KVHS(IK,IR) = 0.0
c
c         Temp at 3LS
c
          if (spredo.lt.s.and.spredo.ge.sprev) then
             if (spredo.lt.((s+sprev)/2.0)) then
                if (ik.eq.1) then
                   kti3ls(idds(ir,2)) = ktids(idds(ir,2))
                else
                   kti3ls(idds(ir,2)) = ktibs(ik-1,ir)
                endif
             else
                kti3ls(idds(ir,2)) = ktibs(ik,ir)
             endif
          endif
c
c         Temp at Sinj
c
          if (sinj.lt.s.and.sinj.ge.sprev) then
             if (sinj.lt.((s+sprev)/2.0)) then
                if (ik.eq.1) then
                   ktinj(idds(ir,2)) = ktids(idds(ir,2))
                else
                   ktinj(idds(ir,2)) = ktibs(ik-1,ir)
                endif
             else
                ktinj(idds(ir,2)) = ktibs(ik,ir)
             endif
          endif
c
          sprev = s
c
c         Calculate and store pressure for the cell
c

          kpress(ik,ir,2) = knbs(ik,ir) * ech * 
     >                     (ktebs(ik,ir) + ktibs(ik,ir))
     >            + crmb * amu * knbs(ik,ir) * kvhs(ik,ir)**2
c

 300    CONTINUE
c
c       ENDIF for IKOPT test
c
        endif 



C
C       INNER PLATES ...
C
        if (ikopt.eq.2.or.ikopt.eq.3) then 
c
        sprev = 0.0
C
        DO 400 IK = ikend, IKMID ,-1
C
          S  = SMAX - KSS(IK,IR)
C
c
C     Calculate revised PINF
c
c     Bug - inner plates should use pinfi not pinf
c
          if (sol13_pdist.gt.0.0) then 
             act_press = pinfi *
     >         min((s/(sol13_pdist*smax)),1.0) *sol13_padd +
     >         pinfi     
          else
             act_press = pinfi * (1.0+sol13_padd)
          endif
c     
c
          IF (CIOPTF.EQ.12) THEN
c
C
C           CALCULATE BACKGROUND TEMPERATURE
C
            KTEBS(IK,IR) = (TEBPI**3.5 + 7.0/(2.0*CK0)*
     >       (LPPAI*S+CIS2(S,SRCRAD,CPOPT,1)))**(2.0/7.0)
C
            KTIBS(IK,IR) = KTEBS(IK,IR)
C
C           CALCULATE DENSITY
C
            SOLI = CIS1(S,SRCION,CSOPT,1)
            GAMMAN = NBPI*V0I + SOLI + RCFI * S
C
            IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)
     >        GAMMAN = 0.0
C
C
            call calcnv(dble(ktebs(ik,ir)),dble(ktibs(ik,ir)),
     >            gamman,act_press,n,v)
c
            knbs(ik,ir) = n
            kvhs(ik,ir) = -v
c
C
C
          ELSEIF (CIOPTF.EQ.13.OR.CIOPTF.EQ.15) THEN
c
c
c
C
C           SET DIFFERENCE BETWEEN 13 AND 15
C
            IF (CIOPTF.EQ.13) THEN
               MFACT = 7.0/2.0
               MFACT2 = 1.0
            ELSEIF (CIOPTF.EQ.15) THEN
               MFACT = 7.0/4.0
               MFACT2 = S / (SMAX/2.0)
            ENDIF
C
C           CALCULATE ELECTRON TEMPERATURE
C
            KTEBS(IK,IR) = (TEBPI**3.5 + (MFACT/CK0)*
     >       (LPPAEI*S*MFACT2
     >       +CIS2(S,SRCRAD,CPOPT,1)))**(2.0/7.0)
C
C           CALCULATE ION TEMPERATURE
C
            KTIBS(IK,IR) = (TIBPI**3.5 + (MFACT/CK0I)*
     >       (LPPAII*S*MFACT2))**(2.0/7.0)
C
C           CALCULATE DENSITY
C
            SOLI = CIS1(S,SRCION,CSOPT,1)
            GAMMAN = NBPI*V0I + SOLI + RCFI * S
C
            IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)
     >        GAMMAN = 0.0
C
            call calcnv(dble(ktebs(ik,ir)),dble(ktibs(ik,ir)),
     >            gamman,act_press,n,v)
c
            knbs(ik,ir) = n
            kvhs(ik,ir) = -v
C
c
c
          ELSEIF (CIOPTF.EQ.14) THEN
C
C
c
            SOLI = CIS1(S,SRCION,CSOPT,1)
            GAMMAN = NBPI*V0I + SOLI + RCFI * S
C
            IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)
     >        GAMMAN = 0.0
C
C           CALCULATE THE BACKGROUND TEMPERATURE
C
            IF (IK.EQ.NKS(IR).AND.S.LE.0.0) THEN
              KTEBS(IK,IR) = TEBPI
            ELSEIF (IK.EQ.NKS(IR).AND.S.GT.0.0) THEN
              ARG1 = (SMAX-KSS(IK,IR)) * 5.0 * GAMMAN * ECH
              solprn = CIS1(S,SRCRAD,CPOPT,1)
              ARG2 = (SMAX-KSS(IK,IR)) * ( LPPAI + solprn)
              ARG3 = TEBPI * CK0
              ARG4 = -CK0
              KTEBS(IK,IR) = SOLVTEMP(ARG1,ARG2,ARG3,ARG4)
            ELSE
              DELTAS = KSS(IK+1,IR) - KSS(IK,IR)
              ARG1 = DELTAS * 5.0 * GAMMAN * ECH
              solprn = CIS1(S,SRCRAD,CPOPT,1)
              ARG2 = DELTAS * ( LPPAI + solprn)
              ARG3 = KTEBS(IK+1,IR) * CK0
              ARG4 = -CK0
              KTEBS(IK,IR) = SOLVTEMP(ARG1,ARG2,ARG3,ARG4)
            ENDIF
C
C            WRITE (6,*) 'TEMPS:',IK,IR,KTEBS(IK,IR)
C
            KTIBS(IK,IR) = KTEBS(IK,IR)
C
C           CALCULATE DENSITY
C
            call calcnv(dble(ktebs(ik,ir)),dble(ktibs(ik,ir)),
     >            gamman,act_press,n,v)
c
            knbs(ik,ir) = n
            kvhs(ik,ir) = -v
c
c
c
          elseif (cioptf.eq.16) then
c
c
c
c
c           This SOL option uses a slightly different calculational
c           method but should yield results comparable to SOL 12.
c
            if (ik.eq.nks(ir)) then
              deltas = (smax-kss(ik,ir))/msolpt
              strt   = kss(ik,ir)
            else
              deltas = (kss(ik+1,ir)-kss(ik,ir))/msolpt
              strt   = kss(ik,ir)
            endif
            do 350 in = msolpt,1,-1
c
              sact = strt + deltas * real(in-1)
              stmp = smax - sact
              ikn = ik*msolpt + (in-1)
              solcorcur = sact
c
              if (ikn.eq.nks(ir)*msolpt+msolpt-1) then
                 soltelast = tebpi
                 soltilast = tibpi
                 solnelast = nbpi
                 solvellast = v0i
                 solcorlast = smax
                 if (ir.eq.cirhr) then
                   solte(ikn+1) = tebpi
                   solti(ikn+1) = tibpi
                   solne(ikn+1) = nbpi
                   solvel(ikn+1) = v0i
                   solcor(ikn+1) = smax
                 endif
              endif
c
              SOLI = CIS1(Stmp,SRCION,CSOPT,1)
              GAMMAN = NBPi*V0i + SOLI + RCFi * Stmp
C
              IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)
     >          GAMMAN = 0.0
C
C             CALCULATE THE BACKGROUND TEMPERATURE
C
              solprn = cis1(stmp,srcrad,cpopt,1)
              soltecur =
     >            - ((lppai + solprn)
     >              *deltas
     >              + ck0*soltelast**3.5)
     >            / ( 5.0 * gamman * deltas * ech
     >              - ck0*soltelast**2.5)
C
              solticur = soltecur
C
C             CALCULATE DENSITY
C
              call calcnv(soltecur,solticur,
     >             gamman,act_press,solnecur,solvelcur)
              solvelcur = -solvelcur
c
c             Store hi-res background if on selected ring
c
              if (ir.eq.cirhr) then
                 solte(ikn) = soltecur
                 solti(ikn) = solticur
                 solne(ikn) = solnecur
                 solvel(ikn) = solvelcur
                 solcor(ikn) = solcorcur
              endif
c
c             Update old variables for next loop iteration
c
              soltelast = soltecur
              soltilast = solticur
              solnelast = solnecur
              solvellast = solvelcur
              solcorlast = solcorcur
c
 350        continue
c
c           Assign values to actual grid points
c
            ktebs(ik,ir) = soltecur
            ktibs(ik,ir) = solticur
            knbs(ik,ir) = solnecur
            kvhs(ik,ir) = solvelcur
c
c           Finish calculation to mid-plane
c
            if (ik.eq.ikmid) then
              strt = kss(ik-1,ir)
              deltas = (kss(ik,ir) - kss(ik-1,ir))/msolpt
              do 360 in = msolpt,msolpt/2+1,-1
                sact = strt + deltas * real(in-1)
                stmp = smax - sact
                ikn = (ik-1) * msolpt + in-1
                solcorcur = sact
c
                SOLI = CIS1(Stmp,SRCION,CSOPT,1)
                GAMMAN = NBPi*V0i + SOLI + RCFi * Stmp
C
                IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)
     >            GAMMAN = 0.0
C
C               CALCULATE THE BACKGROUND TEMPERATURE
C
                solprn = cis1(stmp,srcrad,cpopt,1)
                soltecur =
     >             - ((lppai + solprn)
     >               *deltas
     >               + ck0*soltelast**3.5)
     >             / ( 5.0 * gamman * deltas * ech
     >               - ck0*soltelast**2.5)
C
                solticur = soltecur
C
C               CALCULATE DENSITY
C
                call calcnv(soltecur,solticur,
     >             gamman,act_press,solnecur,solvelcur)
                solvelcur = -solvelcur
c
c               Store hi-res background if on selected ring
c
                if (ir.eq.cirhr) then
                   solte(ikn) = soltecur
                   solti(ikn) = solticur
                   solne(ikn) = solnecur
                   solvel(ikn) = solvelcur
                   solcor(ikn) = solcorcur
                endif
c
c               Update old variables for next loop iteration
c
                soltelast = soltecur
                soltilast = solticur
                solnelast = solnecur
                solvellast = solvelcur
                solcorlast = solcorcur
c
 360          continue
            endif
c
c
c
          elseif (cioptf.eq.17) then
c
c
c
c
c           This SOL option uses a slightly different calculational
c           method but should yield results comparable to SOL 12.
c
            if (ik.eq.nks(ir)) then
              deltas = (smax-kss(ik,ir))/msolpt
              strt   = kss(ik,ir)
            else
              deltas = (kss(ik+1,ir)-kss(ik,ir))/msolpt
              strt   = kss(ik,ir)
            endif
            do 370 in = msolpt,1,-1
c
              sact = strt + deltas * real(in-1)
              stmp = smax - sact
              ikn = ik*msolpt + (in-1)
              solcorcur = sact
c
              if (ikn.eq.nks(ir)*msolpt+msolpt-1) then
                 soltelast = tebpi
                 soltilast = tibpi
                 solnelast = nbpi
                 solvellast = v0i
                 solcorlast = smax
                 if (ir.eq.cirhr) then
                   solte(ikn+1) = tebpi
                   solti(ikn+1) = tibpi
                   solne(ikn+1) = nbpi
                   solvel(ikn+1) = v0i
                   solcor(ikn+1) = smax
                 endif
              endif
c
              SOLI = CIS1(Stmp,SRCION,CSOPT,1)
              GAMMAN = NBPi*V0i + SOLI + RCFi * Stmp
C
              IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)
     >          GAMMAN = 0.0
C
C             CALCULATE THE BACKGROUND TEMPERATURE
C
              solprn = cis1(stmp,srcrad,cpopt,1)
              soltecur =
     >            - ((lppaei + solprn)
     >              *deltas
     >              + ck0*soltelast**3.5)
     >            / ( 2.5 * gamman * deltas * ech
     >              - ck0*soltelast**2.5)
C
              solticur =
     >            - (lppaii
     >              *deltas
     >              + ck0i*soltilast**3.5)
     >            / ( 2.5 * gamman * deltas * ech
     >              - ck0i*soltilast**2.5)
C
C             CALCULATE DENSITY
C
              call calcnv(soltecur,solticur,
     >             gamman,act_press,solnecur,solvelcur)
              solvelcur = -solvelcur
c
c             Store hi-res background if on selected ring
c
              if (ir.eq.cirhr) then
                 solte(ikn) = soltecur
                 solti(ikn) = solticur
                 solne(ikn) = solnecur
                 solvel(ikn) = solvelcur
                 solcor(ikn) = solcorcur
              endif
c
c             Update old variables for next loop iteration
c
              soltelast = soltecur
              soltilast = solticur
              solnelast = solnecur
              solvellast = solvelcur
              solcorlast = solcorcur
c
 370        continue
c
c           Assign values to actual grid points
c
            ktebs(ik,ir) = soltecur
            ktibs(ik,ir) = solticur
            knbs(ik,ir) = solnecur
            kvhs(ik,ir) = solvelcur
c
c           Finish calculation to mid-plane
c
            if (ik.eq.ikmid) then
              strt = kss(ik-1,ir)
              deltas = (kss(ik,ir) - kss(ik-1,ir))/msolpt
              do 380 in = msolpt,msolpt/2+1,-1
                sact = strt + deltas * real(in-1)
                stmp = smax - sact
                ikn = (ik-1) * msolpt + in-1
                solcorcur = sact
c
                SOLI = CIS1(Stmp,SRCION,CSOPT,1)
                GAMMAN = NBPi*V0i + SOLI + RCFi * Stmp
C
                IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)
     >            GAMMAN = 0.0
C
C               CALCULATE THE BACKGROUND TEMPERATURE
C
                solprn = cis1(stmp,srcrad,cpopt,1)
                soltecur =
     >             - ((lppaei + solprn)
     >               *deltas
     >               + ck0*soltelast**3.5)
     >             / ( 2.5 * gamman * deltas * ech
     >               - ck0*soltelast**2.5)
C
                solticur =
     >            - (lppaii
     >              *deltas
     >              + ck0i*soltilast**3.5)
     >            / ( 2.5 * gamman * deltas * ech
     >              - ck0i*soltilast**2.5)
C
C               CALCULATE DENSITY
C
                call calcnv(soltecur,solticur,
     >             gamman,act_press,solnecur,solvelcur)
                solvelcur = -solvelcur
c
c               Store hi-res background if on selected ring
c
                if (ir.eq.cirhr) then
                   solte(ikn) = soltecur
                   solti(ikn) = solticur
                   solne(ikn) = solnecur
                   solvel(ikn) = solvelcur
                   solcor(ikn) = solcorcur
                endif
c
c               Update old variables for next loop iteration
c
                soltelast = soltecur
                soltilast = solticur
                solnelast = solnecur
                solvellast = solvelcur
                solcorlast = solcorcur
c
 380          continue
            endif
c
c
c
          elseif (cioptf.eq.18) then
c
c
c
c           Adds the kinetic term 1/2mv**2*gamma to the
c           solution of SOL16.
c
            if (ik.eq.nks(ir)) then
              deltas = (smax-kss(ik,ir))/msolpt
              strt   = kss(ik,ir)
            else
              deltas = (kss(ik+1,ir)-kss(ik,ir))/msolpt
              strt   = kss(ik,ir)
            endif
            do 385 in = msolpt,1,-1
c
              sact = strt + deltas * real(in-1)
              stmp = smax - sact
              ikn = ik*msolpt + (in-1)
              solcorcur = sact
c
              if (ikn.eq.nks(ir)*msolpt+msolpt-1) then
                 soltelast = tebpi
                 soltilast = tibpi
                 solnelast = nbpi
                 solvellast = v0i
                 solcorlast = smax
                 if (ir.eq.cirhr) then
                   solte(ikn+1) = tebpi
                   solti(ikn+1) = tibpi
                   solne(ikn+1) = nbpi
                   solvel(ikn+1) = v0i
                   solcor(ikn+1) = smax
                 endif
              endif
c
              SOLI = CIS1(Stmp,SRCION,CSOPT,1)
              GAMMAN = NBPi*V0i + SOLI + RCFi * Stmp
C
              IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)
     >          GAMMAN = 0.0
C
C             CALCULATE THE BACKGROUND TEMPERATURE
C
              solprn = cis1(stmp,srcrad,cpopt,1)
              soltecur =
     >            - ((lppai + solprn
     >               + 0.5*massi*solvellast**2*gamman)
     >              *deltas
     >              + ck0*soltelast**3.5)
     >            / ( 5.0 * gamman * deltas * ech
     >              - ck0*soltelast**2.5)
C
              solticur = soltecur
C
C             CALCULATE DENSITY
C
              call calcnv(soltecur,solticur,
     >             gamman,act_press,solnecur,solvelcur)
              solvelcur = -solvelcur
c
c             Store hi-res background if on selected ring
c
              if (ir.eq.cirhr) then
                 solte(ikn) = soltecur
                 solti(ikn) = solticur
                 solne(ikn) = solnecur
                 solvel(ikn) = solvelcur
                 solcor(ikn) = solcorcur
              endif
c
c             Update old variables for next loop iteration
c
              soltelast = soltecur
              soltilast = solticur
              solnelast = solnecur
              solvellast = solvelcur
              solcorlast = solcorcur
c
 385        continue
c
c           Assign values to actual grid points
c
            ktebs(ik,ir) = soltecur
            ktibs(ik,ir) = solticur
            knbs(ik,ir) = solnecur
            kvhs(ik,ir) = solvelcur
c
c           Finish calculation to mid-plane
c
            if (ik.eq.ikmid) then
              strt = kss(ik-1,ir)
              deltas = (kss(ik,ir) - kss(ik-1,ir))/msolpt
              do 390 in = msolpt,msolpt/2+1,-1
                sact = strt + deltas * real(in-1)
                stmp = smax - sact
                ikn = (ik-1) * msolpt + in-1
                solcorcur = sact
c
                SOLI = CIS1(Stmp,SRCION,CSOPT,1)
                GAMMAN = NBPi*V0i + SOLI + RCFi * Stmp
C
                IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)
     >            GAMMAN = 0.0
C
C               CALCULATE THE BACKGROUND TEMPERATURE
C
                solprn = cis1(stmp,srcrad,cpopt,1)
                soltecur =
     >             - ((lppai + solprn
     >                 + 0.5*massi*solvellast**2*gamman)
     >               *deltas
     >               + ck0*soltelast**3.5)
     >             / ( 5.0 * gamman * deltas * ech
     >               - ck0*soltelast**2.5)
C
                solticur = soltecur
C
C               CALCULATE DENSITY
C
                call calcnv(soltecur,solticur,
     >             gamman,act_press,solnecur,solvelcur)
                solvelcur = -solvelcur
c
c               Store hi-res background if on selected ring
c
                if (ir.eq.cirhr) then
                   solte(ikn) = soltecur
                   solti(ikn) = solticur
                   solne(ikn) = solnecur
                   solvel(ikn) = solvelcur
                   solcor(ikn) = solcorcur
                endif
c
c               Update old variables for next loop iteration
c
                soltelast = soltecur
                soltilast = solticur
                solnelast = solnecur
                solvellast = solvelcur
                solcorlast = solcorcur
c
 390          continue
            endif
c
c
c
          elseif (cioptf.eq.19) then
c
c
c
c           Adds the kinetic term 1/2mv**2*gamma to the
c           solution of SOL16.
c
            if (ik.eq.nks(ir)) then
              deltas = (smax-kss(ik,ir))/msolpt
              strt   = kss(ik,ir)
            else
              deltas = (kss(ik+1,ir)-kss(ik,ir))/msolpt
              strt   = kss(ik,ir)
            endif
            do 2000 in = msolpt,1,-1
c
              sact = strt + deltas * real(in-1)
              stmp = smax - sact
              ikn = ik*msolpt + (in-1)
              solcorcur = sact
c
              if (ikn.eq.nks(ir)*msolpt+msolpt-1) then
                 soltelast = tebpi
                 soltilast = tibpi
                 solnelast = nbpi
                 solvellast = v0i
                 solcorlast = smax
                 if (ir.eq.cirhr) then
                   solte(ikn+1) = tebpi
                   solti(ikn+1) = tibpi
                   solne(ikn+1) = nbpi
                   solvel(ikn+1) = v0i
                   solcor(ikn+1) = smax
                 endif
              endif
c
              SOLI = CIS1(Stmp,SRCION,CSOPT,1)
              GAMMAN = NBPi*V0i + SOLI + RCFi * Stmp
C
              IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)
     >          GAMMAN = 0.0
C
C             CALCULATE THE BACKGROUND TEMPERATURE
C
              helpi = 17.5 + (5.0+37.5/soltelast)
     >                * (1.0+0.25/soltelast)
     >                * log10(1.0d21/solnelast)
c
              solprn = cis1(stmp,srcrad,cpopt,1)
              soltecur =
     >            - ((lppai + solprn
     >               + 0.5*massi*solvellast**2*gamman
     >               + helpi * soli * 1.6d-19)
     >              *deltas
     >              + ck0*soltelast**3.5)
     >            / ( 5.0 * gamman * deltas * ech
     >              - ck0*soltelast**2.5)
C
              solticur = soltecur
C
C             CALCULATE DENSITY
C
              call calcnv(soltecur,solticur,
     >             gamman,act_press,solnecur,solvelcur)
              solvelcur = -solvelcur
c
c             Store hi-res background if on selected ring
c
              if (ir.eq.cirhr) then
                 solte(ikn) = soltecur
                 solti(ikn) = solticur
                 solne(ikn) = solnecur
                 solvel(ikn) = solvelcur
                 solcor(ikn) = solcorcur
              endif
c
c             Update old variables for next loop iteration
c
              soltelast = soltecur
              soltilast = solticur
              solnelast = solnecur
              solvellast = solvelcur
              solcorlast = solcorcur
c
 2000       continue
c
c           Assign values to actual grid points
c
            ktebs(ik,ir) = soltecur
            ktibs(ik,ir) = solticur
            knbs(ik,ir) = solnecur
            kvhs(ik,ir) = solvelcur
c
c           Finish calculation to mid-plane
c
            if (ik.eq.ikmid) then
              strt = kss(ik-1,ir)
              deltas = (kss(ik,ir) - kss(ik-1,ir))/msolpt
              do 2010 in = msolpt,msolpt/2+1,-1
                sact = strt + deltas * real(in-1)
                stmp = smax - sact
                ikn = (ik-1) * msolpt + in-1
                solcorcur = sact
c
                SOLI = CIS1(Stmp,SRCION,CSOPT,1)
                GAMMAN = NBPi*V0i + SOLI + RCFi * Stmp
C
                IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)
     >            GAMMAN = 0.0
C
C               CALCULATE THE BACKGROUND TEMPERATURE
C
                helpi = 17.5 + (5.0+37.5/soltelast)
     >                * (1.0+0.25/soltelast)
     >                * log10(1.0d21/solnelast)
c
                solprn = cis1(stmp,srcrad,cpopt,1)
                soltecur =
     >             - ((lppai + solprn
     >                 + 0.5*massi*solvellast**2*gamman
     >                 + helpi * soli *1.6d-19)
     >               *deltas
     >               + ck0*soltelast**3.5)
     >             / ( 5.0 * gamman * deltas * ech
     >               - ck0*soltelast**2.5)
C
                solticur = soltecur
C
C               CALCULATE DENSITY
C
                call calcnv(soltecur,solticur,
     >             gamman,act_press,solnecur,solvelcur)
                solvelcur = -solvelcur
c
c               Store hi-res background if on selected ring
c
                if (ir.eq.cirhr) then
                   solte(ikn) = soltecur
                   solti(ikn) = solticur
                   solne(ikn) = solnecur
                   solvel(ikn) = solvelcur
                   solcor(ikn) = solcorcur
                endif
c
c               Update old variables for next loop iteration
c
                soltelast = soltecur
                soltilast = solticur
                solnelast = solnecur
                solvellast = solvelcur
                solcorlast = solcorcur
c
 2010         continue
            endif
c
c
c
          elseif (cioptf.eq.20) then
c
c
c
c           Adds the kinetic term 1/2mv**2*gamma to the
c           solution of SOL16.
c
            maxind = nks(ir)*msolpt+msolpt
            if (ik.eq.nks(ir)) then
              deltas = (smax-kss(ik,ir))/msolpt
              strt   = kss(ik,ir)
            else
              deltas = (kss(ik+1,ir)-kss(ik,ir))/msolpt
              strt   = kss(ik,ir)
            endif
c
c           Initialize Pei source term
c
            solpein = srcpei(stmp,ikn,1,deltas,1)
c
            write(6,*) 'deltas:',deltas,strt
            do 2020 in = msolpt,1,-1
c
              sact = strt + deltas * real(in-1)
              stmp = smax - sact
              ikn = ik*msolpt + (in-1)
              cind = ikn
              solcorcur = sact
c
              if (ikn.eq.nks(ir)*msolpt+msolpt-1) then
                 soltelast = tebpi
                 soltilast = tibpi
                 solnelast = nbpi
                 solvellast = v0i
                 solcorlast = smax
c
c                Hi-res background for calculating
c                Pei(s) and its integral.
c
                 soltehr(ikn+1) = tebpi
                 soltihr(ikn+1) = tibpi
                 solnehr(ikn+1) = nbpi
                 solcorhr(ikn+1) = 0.0
c
                 if (ir.eq.cirhr) then
                   solte(ikn+1) = tebpi
                   solti(ikn+1) = tibpi
                   solne(ikn+1) = nbpi
                   solvel(ikn+1) = v0i
                   solcor(ikn+1) = smax
                 endif
              endif
c
              SOLI = CIS1(Stmp,SRCION,CSOPT,1)
              GAMMAN = NBPi*V0i + SOLI + RCFi * Stmp
C
              IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)
     >          GAMMAN = 0.0
C
C             CALCULATE THE BACKGROUND TEMPERATURE
C
              helpi = 17.5 + (5.0+37.5/soltelast)
     >                * (1.0+0.25/soltelast)
     >                * log10(1.0d21/solnelast)
c
              solprn = cis1(stmp,srcrad,cpopt,1)
              solpein = srcpei(stmp,ikn,1,deltas,0)
              soltecur =
     >            - ((lppaei + solprn + solpein
     >                + 0.5*massi*solvellast**2*gamman
     >                + helpi * soli * 1.6d-19)
     >              *deltas
     >              + ck0*soltelast**3.5)
     >            / ( 2.5 * gamman * deltas * ech
     >              - ck0*soltelast**2.5)
              solticur =
     >            - ((lppaii -solpein)
     >              *deltas
     >              + ck0i*soltilast**3.5)
     >            / ( 2.5 * gamman * deltas * ech
     >              - ck0i*soltilast**2.5)
C
C             CALCULATE DENSITY
C
              call calcnv(soltecur,solticur,
     >             gamman,act_press,solnecur,solvelcur)
              solvelcur = -solvelcur
c
c             Store hi-res background if on selected ring
c
              soltehr(ikn) = soltecur
              soltihr(ikn) = solticur
              solnehr(ikn) = solnecur
              solcorhr(ikn) = solcorcur
              if (ir.eq.cirhr) then
                 solte(ikn) = soltecur
                 solti(ikn) = solticur
                 solne(ikn) = solnecur
                 solvel(ikn) = solvelcur
                 solcor(ikn) = solcorcur
              endif
c
c             Update old variables for next loop iteration
c
              soltelast = soltecur
              soltilast = solticur
              solnelast = solnecur
              solvellast = solvelcur
              solcorlast = solcorcur
c
c              write(6,*) 'sol18:t:',ikn,ir,gamman,rootn,
c     >            solte(ikn,ir),solti(ikn,ir),solne(ikn,ir)
c              write(6,*) 'cont:',solvel(ikn,ir),solcor(ikn,ir),
c     >                    tmp
c
c
 2020       continue
c
c           Assign values to actual grid points
c
            ktebs(ik,ir) = soltecur
            ktibs(ik,ir) = solticur
            knbs(ik,ir) = solnecur
            kvhs(ik,ir) = solvelcur
c            write(6,*) 'sol16:',ik,ir,ktebs(ik,ir),ktibs(ik,ir),
c     >                 knbs(ik,ir),kvhs(ik,ir)
c
c           Finish calculation to mid-plane
c
            if (ik.eq.ikmid) then
              strt = kss(ik-1,ir)
              deltas = (kss(ik,ir) - kss(ik-1,ir))/msolpt
              do 2030 in = msolpt,msolpt/2+1,-1
                sact = strt + deltas * real(in-1)
                stmp = smax - sact
                ikn = (ik-1) * msolpt + in-1
                cind = ikn
                solcorcur = sact
c
                SOLI = CIS1(Stmp,SRCION,CSOPT,1)
                GAMMAN = NBPi*V0i + SOLI + RCFi * Stmp
C
                IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)
     >            GAMMAN = 0.0
C
C               CALCULATE THE BACKGROUND TEMPERATURE
C
                helpi = 17.5 + (5.0+37.5/soltelast)
     >                * (1.0+0.25/soltelast)
     >                * log10(1.0d21/solnelast)
c
              solprn = cis1(stmp,srcrad,cpopt,1)
              solpein = srcpei(stmp,ikn,1,deltas,0)
              soltecur =
     >            - ((lppaei + solprn + solpein
     >                + 0.5*massi*solvellast**2*gamman
     >                + helpi * soli * 1.6d-19)
     >              *deltas
     >              + ck0*soltelast**3.5)
     >            / ( 2.5 * gamman * deltas * ech
     >              - ck0*soltelast**2.5)
              solticur =
     >            - ((lppaii -solpein)
     >              *deltas
     >              + ck0i*soltilast**3.5)
     >            / ( 2.5 * gamman * deltas * ech
     >              - ck0i*soltilast**2.5)
C
C               CALCULATE DENSITY
C
                call calcnv(soltecur,solticur,
     >             gamman,act_press,solnecur,solvelcur)
                solvelcur = -solvelcur
c
c               Store hi-res background if on selected ring
c
                soltehr(ikn) = soltecur
                soltihr(ikn) = solticur
                solnehr(ikn) = solnecur
                solcorhr(ikn) = solcorcur
                if (ir.eq.cirhr) then
                   solte(ikn) = soltecur
                   solti(ikn) = solticur
                   solne(ikn) = solnecur
                   solvel(ikn) = solvelcur
                   solcor(ikn) = solcorcur
                endif
c
c               Update old variables for next loop iteration
c
                soltelast = soltecur
                soltilast = solticur
                solnelast = solnecur
                solvellast = solvelcur
                solcorlast = solcorcur
c
 2030         continue
            endif
c
c         Sol option 21 ... model for a detached plasma.
c
          elseif (cioptf.eq.21) then
c
c           Test S-value and then calculate Te,Ti,N and V
c
            if (s.lt.sl1ai) then
c
              ktebs(ik,ir) = tebpi + (tebpi*ter1ai-tebpi) 
     >                               * (s/sl1ai)
              ktibs(ik,ir) = tibpi + (tibpi*tir1ai-tibpi) 
     >                               * (s/sl1ai)
c
              knbs(ik,ir)  = nbpi + (nbpi*nr1ai-nbpi) 
     >                               * (s/sl1ai)**n_expi
c
              kvhs(ik,ir)  = - (nbpi * v0i)/ knbs(ik,ir) * aux_vel21
c
            elseif (s.lt.sl1bi) then
c
              ktebs(ik,ir) = tebpi*ter1ai 
     >                       + (tebpi*ter1bi-tebpi*ter1ai) 
     >                         * (s-sl1ai)/(sl1bi-sl1ai)
              ktibs(ik,ir) = tibpi*tir1ai 
     >                       + (tibpi*tir1bi-tibpi*tir1ai) 
     >                         * (s-sl1ai)/(sl1bi-sl1ai)
c
              knbs(ik,ir)  = nbpi*nr1ai 
     >                     + (nbpi*nr1bi-nbpi*nr1ai) 
     >                       * ((s-sl1ai)/(sl1bi-sl1ai))**n_expi
c
              kvhs(ik,ir)  = - (nbpi * v0i)/ knbs(ik,ir) * aux_vel21
c
            elseif (s.le.sl1i) then
c
              ktebs(ik,ir) = tebpi*ter1bi 
     >                       + (t1i -tebpi*ter1bi) 
     >                         * (s-sl1bi)/(sl1i-sl1bi)
              ktibs(ik,ir) = tibpi*tir1bi 
     >                       + (ti1i  -tibpi*tir1bi) 
     >                         * (s-sl1bi)/(sl1i-sl1bi)
c
c             knbs(ik,ir)  = nbpi + (n1i-nbpi) * (s/sl1i)**n_expi
c
              knbs(ik,ir)  = nbpi*nr1bi 
     >                     + (n1i-nbpi*nr1bi) 
     >                       * ((s-sl1bi)/(sl1i-sl1bi))**n_expi
c
              kvhs(ik,ir)  = - (nbpi * v0i)/ knbs(ik,ir)
c
            elseif (s.le.sl2i) then
              ktebs(ik,ir) = (t1i**3.5+7.0/(2.0*ck0)*(lppai*(s-sl1i)
     >                    +0.5*(s-sl1i)**2*lpradi))**(2.0/7.0)
              ktibs(ik,ir) = ktebs(ik,ir)
              knbs(ik,ir)  = (n1i*t1i) / ktebs(ik,ir)
              if (s.le.slvi) then
                kvhs(ik,ir)= - v1i * (slvi-s)/(slvi-sl1i)
              else
                kvhs(ik,ir)= 0.0
              endif
            else
              ktebs(ik,ir) = (t2i**3.5+7.0/(2.0*ck0)*
     >                       ((1.0+qrati)*lppai*(s-sl2i)))**(2.0/7.0)
              ktibs(ik,ir) = ktebs(ik,ir)
              knbs(ik,ir)  = (n1i*t1i) / ktebs(ik,ir)
              if (s.le.slvi) then
                kvhs(ik,ir)= - v1i * (slvi-s)/(slvi-sl1i)
              else
                kvhs(ik,ir)= 0.0
              endif
            endif
c
          ENDIF
c
          IF (ABS(KVHS(IK,IR)).LT.1.0) KVHS(IK,IR) = 0.0
c
c         Temp at 3LS
c
          if (spredi.lt.s.and.spredi.ge.sprev) then
             if (spredi.lt.((s+sprev)/2.0)) then
                if (ik.eq.nks(ir)) then
                   kti3ls(idds(ir,1)) = ktids(idds(ir,1))
                else
                   kti3ls(idds(ir,1)) = ktibs(ik+1,ir)
                endif
             else
                kti3ls(idds(ir,1)) = ktibs(ik,ir)
             endif
          endif
c
c         Temp at Sinj
c
          if (sinj.lt.s.and.sinj.ge.sprev) then
             if (sinj.lt.((s+sprev)/2.0)) then
                if (ik.eq.nks(ir)) then
                   ktinj(idds(ir,1)) = ktids(idds(ir,1))
                else
                   ktinj(idds(ir,1)) = ktibs(ik+1,ir)
                endif
             else
                ktinj(idds(ir,1)) = ktibs(ik,ir)
             endif
          endif
c
          sprev = s


c
c         Calculate and store pressure for the cell
c

          kpress(ik,ir,2) = knbs(ik,ir) * ech * 
     >                     (ktebs(ik,ir) + ktibs(ik,ir))
     >            + crmb * amu * knbs(ik,ir) * kvhs(ik,ir)**2
c


c
 400    CONTINUE
c
        endif
C
C       CALCULATE ELECTRIC FIELD
C
C
C       IN THE FOLLOWING EQUATIONS THE FACTOR E CANCELS WITH THE
C       SAME FACTOR USED IN CONVERTING T IN EV TO KT.
C
c       If OFIELD is turned ON then set the electric field to
c       zero for this case. Note: the electric field is
c       initialized to zero in the plasma.d3a module.
c
c       For OFIELD ... 0=off   1=on
c
        if (ofield.eq.0) then
c
c
          if (ikopt.eq.1.or.ikopt.eq.3) then  
c
          if (kss(1,ir).eq.0.0) then
            DS1 = KSS(2,IR) - KSS(1,IR)
            DP1 = (KNBS(2,IR)*KTEBS(2,IR)-KNBS(1,IR)*KTEBS(1,IR))
            DT1 = (KTEBS(2,IR)-KTEBS(1,IR))
            NB1 = 0.5*(KNBS(2,IR)+KNBS(1,IR))
C
            KES(1,IR) = -(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1
c
            KEDS(idds(ir,2))= kes(1,ir)
          else
            DS2 = KSS(2,IR) - KSS(1,IR)
            DP2 = (KNBS(2,IR)*KTEBS(2,IR)-KNBS(1,IR)*KTEBS(1,IR))
            DT2 = (KTEBS(2,IR)-KTEBS(1,IR))
            NB2 = 0.5*(KNBS(2,IR)+KNBS(1,IR))
c
            DS1 = KSS(1,IR)
            DP1 = KNBS(1,IR)*KTEBS(1,IR)-
     >              KNDS(idds(ir,2))*KTEDS(idds(ir,2))
            DT1 = KTEBS(1,IR)-KTEDS(idds(ir,2))
            NB1 = 0.5*(KNBS(1,IR)+KNDS(idds(ir,2)))
c
            KES(1,IR) = 0.5*((-(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1)
     >                    + (-(1/NB2)*DP2/DS2 - 0.71 * DT2/DS2))
c
            KEDS(idds(ir,2)) = -(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1
          endif
c
          endif
C
c
c
          if (ikopt.eq.2.or.ikopt.eq.3) then 
c
          if (kss(nks(ir),ir).eq.ksmaxs(ir)) then
            DS1 = KSS(NKS(IR),IR) - KSS(NKS(IR)-1,IR)
            DP1 = (KNBS(NKS(IR),IR)*KTEBS(NKS(IR),IR)
     >         -KNBS(NKS(IR)-1,IR)*KTEBS(NKS(IR)-1,IR))
            DT1 = (KTEBS(NKS(IR),IR)-KTEBS(NKS(IR)-1,IR))
            NB1 = 0.5*(KNBS(NKS(IR),IR)+KNBS(NKS(IR)-1,IR))
C
            KES(NKS(IR),IR) = -(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1
c
            KEDS(idds(ir,1))= kes(nks(ir),ir)
          else
c
            DS2 = KSS(NKS(IR),IR) - KSS(NKS(IR)-1,IR)
            DP2 = (KNBS(NKS(IR),IR)*KTEBS(NKS(IR),IR)
     >         -KNBS(NKS(IR)-1,IR)*KTEBS(NKS(IR)-1,IR))
            DT2 = (KTEBS(NKS(IR),IR)-KTEBS(NKS(IR)-1,IR))
            NB2 = 0.5*(KNBS(NKS(IR),IR)+KNBS(NKS(IR)-1,IR))
c
            DS1 = ksmaxs(ir) - KSS(nks(ir),IR)
            DP1 = KNDS(idds(ir,1))*KTEDS(idds(ir,1))-
     >             KNBS(nks(ir),IR)*KTEBS(nks(ir),IR)
            DT1 = KTEDS(idds(ir,1))-KTEBS(nks(ir),IR)
            NB1 = 0.5*(KNBS(nks(ir),IR)+KNDS(idds(ir,1)))
c
            KES(nks(ir),IR) = 0.5*((-(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1)
     >                    + (-(1/NB2)*DP2/DS2 - 0.71 * DT2/DS2))
c
            KEDS(idds(ir,1)) = -(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1
          endif
c
          endif

C
C        WRITE(6,*) 'KES:IR:',KES(1,IR),KES(NKS(IR),IR)
C
          if (ikopt.eq.1) then 
            ikfirst = ikstart+1
            iklast = ikmid -1
          elseif(ikopt.eq.2) then
            ikfirst = ikmid
            iklast = ikend-1
c slmod begin
c...BUG: Aug 8, 2000 -SL
          elseif(ikopt.eq.3) then
c
c          elseif(ikopt.eq.2) then
c slmod end
            ikfirst = ikstart
            iklast = ikend
          endif
c
          DO 500 IK = ikfirst,iklast
c slmod begin
c...BUG?
            IF (IK.EQ.1      ) THEN
              DS1 =KSS(IK,IR) - KSB(IK-1,IR)
              DP1 =KNBS(IK,IR)     *KTEBS(IK,IR)-
     .             KNDS(IDDS(IR,2))*KTEDS(IDDS(IR,2))
              DT1 =(KTEBS(IK,IR)-KTEDS(IDDS(IR,2)))
              NB1 =0.5*(KNBS(IK,IR)+KNDS(IDDS(IR,2)))
            ELSE
              DS1 =KSS(IK,IR) - KSS(IK-1,IR)
              DP1 =KNBS(IK,IR)*KTEBS(IK,IR)-KNBS(IK-1,IR)*KTEBS(IK-1,IR)
              DT1 =(KTEBS(IK,IR)-KTEBS(IK-1,IR))
              NB1 =0.5*(KNBS(IK,IR)+KNBS(IK-1,IR))
            ENDIF

            IF (IK.EQ.NKS(IR)) THEN
              DS2 =KSB(IK,IR) - KSS(IK,IR)
              DP2 =KNDS(IDDS(IR,1))*KTEDS(IDDS(IR,1))-
     .             KNBS(IK,IR)     *KTEBS(IK,IR)
              DT2 =KTEDS(IDDS(IR,1))-KTEBS(IK,IR)
              NB2 =0.5*(KNDS(IDDS(IR,1))+KNBS(IK,IR))
            ELSE
              DS2 =KSS(IK+1,IR) - KSS(IK,IR)
              DP2 =KNBS(IK+1,IR)*KTEBS(IK+1,IR)-KNBS(IK,IR)*KTEBS(IK,IR)
              DT2 =(KTEBS(IK+1,IR)-KTEBS(IK,IR))
              NB2 =0.5*(KNBS(IK+1,IR)+KNBS(IK,IR))
            ENDIF
c
c            DS1 = KSS(IK,IR) - KSS(IK-1,IR)
c            DP1 = KNBS(IK,IR)*KTEBS(IK,IR)-KNBS(IK-1,IR)*KTEBS(IK-1,IR)
c            DT1 = (KTEBS(IK,IR)-KTEBS(IK-1,IR))
c            NB1 = 0.5*(KNBS(IK,IR)+KNBS(IK-1,IR))
c            DS2 = KSS(IK+1,IR) - KSS(IK,IR)
c            DP2 = KNBS(IK+1,IR)*KTEBS(IK+1,IR)-KNBS(IK,IR)*KTEBS(IK,IR)
c            DT2 = (KTEBS(IK+1,IR)-KTEBS(IK,IR))
c            NB2 = 0.5*(KNBS(IK+1,IR)+KNBS(IK,IR))
c slmod end
            if (nb1.eq.0.0.or.nb2.eq.0.0.or.ds1.eq.0.0.or.ds2.eq.0.0) 
     >         then 
                 kes(ik,ir) = 0.0
            else 
               KES(IK,IR) = 0.5*((-(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1)
     >                    + (-(1/NB2)*DP2/DS2 - 0.71 * DT2/DS2))
            endif 
C
C            WRITE(6,*) 'KES:',IK,IR,KES(IK,IR)
C
 
 
 500      CONTINUE
c
c         End of test on OFIELD
c
        endif
c
c
c       SOL 21 does not involve an ionization source
c       Therefore - do not calculate these values.
c
        if (cioptf.ne.21) then
c
c       Calculate Near-plate predictor barrier values
c
        predls = csolls * smax
c
c       Barrier value array contents kpredbar(x,y,z)
c       x - target position index
c       y - 1,2,3 the three components of the barrier quantity
c       z - 1 at S=0, 2 at S=Sinj
c
c
c       Barrier at S=0
c
        kpredbar(idds(ir,2),1,1) =
     >    -2.35e-15*knds(idds(ir,2))*predls/ktids(idds(ir,2))**2
        kpredbar(idds(ir,2),2,1) =
     >     31.8 * log(kti3ls(idds(ir,2))/ktids(idds(ir,2)))
        kpredbar(idds(ir,2),3,1) = kpredbar(idds(ir,2),1,1)
     >                         + kpredbar(idds(ir,2),2,1)
c
        kpredbar(idds(ir,1),1,1) =
     >    -2.35e-15*knds(idds(ir,1))*predls/ktids(idds(ir,1))**2
        kpredbar(idds(ir,1),2,1) =
     >     31.8 * log(kti3ls(idds(ir,1))/ktids(idds(ir,1)))
        kpredbar(idds(ir,1),3,1) = kpredbar(idds(ir,1),1,1)
     >                         + kpredbar(idds(ir,1),2,1)
c
c       Barrier at S=Sinj
c
        kpredbar(idds(ir,2),1,2) =
     >    -2.35e-15*knds(idds(ir,2))*predls/ktids(idds(ir,2))**2
     >    * exp(-sinj/predls)
        kpredbar(idds(ir,2),2,2) =
     >     31.8 * log(kti3ls(idds(ir,2))/ktinj(idds(ir,2)))
        kpredbar(idds(ir,2),3,2) = kpredbar(idds(ir,2),1,2)
     >                         + kpredbar(idds(ir,2),2,2)
c
        kpredbar(idds(ir,1),1,2) =
     >    -2.35e-15*knds(idds(ir,1))*predls/ktids(idds(ir,1))**2
     >    * exp(-sinj/predls)
        kpredbar(idds(ir,1),2,2) =
     >     31.8 * log(kti3ls(idds(ir,1))/ktinj(idds(ir,1)))
        kpredbar(idds(ir,1),3,2) = kpredbar(idds(ir,1),1,2)
     >                         + kpredbar(idds(ir,1),2,2)
c
c
c        write(6,*) 'pred:',predls,spredo,spredi,kti3ls(idds(ir,2)),
c     >         kti3ls(idds(ir,1)),ir,idds(ir,2),idds(ir,1)
c        write(6,*) 'pred2:',kpredbar(idds(ir,2)),
c     >                    kpredbar(idds(ir,1))
c        write(6,*) 'O1:', 31.8 * log(kti3ls(idds(ir,2))
c     >                    /ktids(idds(ir,2))),ktids(idds(ir,2))
c        write(6,*) 'O2:',  -2.35e-15*knds(idds(ir,2))*predls
c     >           /ktids(idds(ir,2))**2 ,knds(idds(ir,2))
c        write(6,*) 'I1:', 31.8 * log(kti3ls(idds(ir,1))
c     >                    /ktids(idds(ir,1))),ktids(idds(ir,1))
c        write(6,*) 'I2:',  -2.35e-15*knds(idds(ir,1))*predls
c     >           /ktids(idds(ir,1))**2 ,knds(idds(ir,1))
c
c
c        Endif for SOL 21 test.
c
         endif
c
c
600   CONTINUE
c
c jdemod
c     MOVED to BGPLASMA module and controlled by bg_velocity_opt   
c
c slmod begin 
c     Adding some capacity to over-ride the SOL21 calcualted hydrogenic flow
c     velocity.  This should be moved in the future, preferentially to some
c     generic "over-ride" routine that is called after the solution has been
c     calculated, irrespective of the solver that is in use: 
c      IF (cioptf.EQ.21.AND.osmns28.GT.0) THEN
c        WRITE(0,*) 'CALLING vb OVER-RIDE CODE'
c
c        CALL PrescribeFlow 
c      ENDIF
c slmod end 
c jdemod
c
      RETURN
 9100 FORMAT('  DECAY VALUE LS FOR IONIZATION SOURCE ADJUSTED: ',/,
     >'   VALUE FOR RINGS ',I4,' AND GREATER NOW ',G16.7)
 
      END
C
C
      REAL FUNCTION SOLVTEMP (ARG1,ARG2,ARG3,ARG4)
      IMPLICIT NONE
      DOUBLE PRECISION ARG1,ARG2,ARG3,ARG4
C
C     THIS IS A FRONT-END TO THE EQUATION SOLVING ROUTINES.
C     IT'S SOLE PURPOSE IS TO SEGREGATE THE ROUTINES USED TO
C     SOLVE THE EQUATION FOR THE TEMPERATURE. THE ARGUMENTS ARE
C     THE COEFFICIENTS OF THE TEMPERATURE EQUATION WITH THE
C     TEMPERATURE IN EV.
C
      INTEGER MAXROOTS,MAXINTS
      REAL LOWTEMP,HIGHTEMP,XACC
      PARAMETER (MAXROOTS=7,MAXINTS=200,LOWTEMP=0.0,
     >           HIGHTEMP=200.0,XACC=1.0E-3)
      DOUBLE PRECISION FARG1,FARG2,FARG3,FARG4
      COMMON /FPARAMS/ FARG1,FARG2,FARG3,FARG4
      REAL XB1(MAXROOTS),XB2(MAXROOTS)
      INTEGER NB,N,I,J,K
      REAL SOL14FX,RTBIS
      EXTERNAL SOL14FX,RTBIS
C
      FARG1 = ARG1
      FARG2 = ARG2
      FARG3 = ARG3
      FARG4 = ARG4
C
      NB = MAXROOTS
      N = MAXINTS
      CALL ZBRAK(SOL14FX,LOWTEMP,HIGHTEMP,N,XB1,XB2,NB)
      IF  (NB.LE.0) THEN
        WRITE(6,*)
     > 'SOLVTEMP: ERROR IN EQUATION SOLVER - NO T ROOT FOUND IN RANGE'
         STOP
      ELSEIF (NB.GT.1) THEN
        WRITE(6,*) 'SOLVTEMP: ERROR - MULTIPLE ROOTS FOUND',NB
        DO 10 I = 1, NB
          WRITE(6,*) 'ROOT IN THE FOLLOWING BRACKETS:',XB1(I),XB2(I)
10      CONTINUE
C
C       SOLVE FOR LARGEST AVAILABLE ROOT
C
        SOLVTEMP = RTBIS(SOL14FX,XB1(NB),XB2(NB),XACC)
 
      ELSE
C
C       SOLVE FOR ROOT
C
        SOLVTEMP = RTBIS(SOL14FX,XB1(NB),XB2(NB),XACC)
      ENDIF
 
      RETURN
      END
C
C
      REAL FUNCTION SOL14FX(T2)
      IMPLICIT NONE
      REAL T2
C
C     THIS CALCULATES THE FUNCTIONAL VALUE OF F(T2) WHICH FOR
C     A ROOT WILL EQUAL ZERO.
C
      DOUBLE PRECISION FARG1,FARG2,FARG3,FARG4
      COMMON /FPARAMS/ FARG1,FARG2,FARG3,FARG4
C
      SOL14FX = FARG1*T2 + FARG4*T2**3.5 + FARG3*T2**2.5 + FARG2
C
      RETURN
      END
C
C
C
      SUBROUTINE ZBRAK(FX,X1,X2,N,XB1,XB2,NB)
      implicit none
C
C     GIVEN A FUNCTION FX DEFINED ON THE INTERVAL FROM X1 TO X2, SUBDIVI
C     THE INTERVAL INTO N EQUALLY SPACED SEGMENTS AND SEARCH FOR ZERO
C     CROSSINGS OF THE FUNCTION. NB IS INPUT AS THE MAXIMUM NUMBER OF RO
C     SOUGHT, AND IS RESET TO THE NUMBER OF BRACKETING PAIRS XB1, XB2 TH
C     ARE FOUND.
C
C     FROM "NUMERICAL RECIPES" - CHAPTER 9
C
      integer n,nb,i,nbb
      EXTERNAL FX
      REAL FX
      REAL XB1(NB),XB2(NB),x1,x2,dx,x,fp,fc
c
      NBB = NB
      NB = 0
      X = X1
      DX = (X2-X1)/N
      FP = FX(X)
      DO 11 I = 1,N
        X = X+DX
        FC = FX(X)
        IF (FC*FP.LT.0.0) THEN
          NB = NB + 1
          XB1(NB) = X-DX
          XB2(NB) = X
        ENDIF
        FP = FC
        IF (NBB.EQ.NB) RETURN
 11   CONTINUE
      RETURN
      END
C
C
C
      REAL FUNCTION RTBIS(FUNC,X1,X2,XACC)
      implicit none
C
C     USING BISECTION, FIND THE ROOT OF A FUNCTION FUNC KNOWN
C     TO LIE BETWEEN X1 AND X2. THE ROOT RETURNED AS RTBIS, WILL BE
C     REFINED UNTIL ITS ACCURACY IS +/- XACC
C
C     FROM "NUMERICAL RECIPES - CHAPTER 9"
C
C     LIMIT OF 40 BISECTIONS
C
      EXTERNAL FUNC
      REAL FUNC,x1,x2,xacc,dx,f,xmid,fmid
      integer jmax,j
      PARAMETER (JMAX=40)
c
c
      FMID = FUNC(X2)
      F = FUNC(X1)
      IF (F*FMID.GE.0.0) THEN
        WRITE(6,*) 'ROOT MUST BE BRACKETED FOR BISECTION'
        rtbis = X1
        RETURN
      ENDIF
      IF (F.LT.0.0) THEN
        RTBIS = X1
        DX = X2-X1
      ELSE
        RTBIS = X2
        DX = X1-X2
      ENDIF
      DO 11 J = 1,JMAX
        DX = DX * 0.5
        XMID = RTBIS+DX
        FMID = FUNC(XMID)
        IF (FMID.LE.0.0) RTBIS = XMID
        IF (ABS(DX).LT.XACC.OR.FMID.EQ.0.0) RETURN
 11   CONTINUE
      WRITE (6,*) 'TOO MANY BISECTIONS'
      RETURN
      END
C
C
C
      DOUBLE PRECISION FUNCTION SRCION(S,SOPT,PLATEOPT,SLIM,IND)
      use mod_params
      use mod_pindata
      use mod_comsol
      IMPLICIT NONE
      DOUBLE PRECISION S,SLIM
      INTEGER SOPT,PLATEOPT,IND
C
C     THIS FUNCTION RETURNS THE VALUE OF THE IONIZATION SOURCE
C     STRENGTH AT THE POSITION S FROM ONE OF TWO PLATES.
C     THIS IS NEEDED TO ENHANCE THE ACCURACY OF THE INTEGRATION
C     ROUTINES - WHICH ARE SOMEWHAT INACCURATE USING THE CURRENT
C     METHOD - HOWEVER - THIS ALTERNATE COULD BE COMPUTATIONALLY
C     EXPENSIVE AND SO THE ALTERNATE CODE WILL BE RETAINED.
C
C     DAVID ELDER    MAY 1, 1992
C
C     INCLUDE "PARAMS"
c     include 'params'
C     INCLUDE "PINDATA"
c     include 'pindata'
C     INCLUDE "COMSOL"
c     include 'comsol'
      INTEGER IPOS,IN
      EXTERNAL IPOS
      REAL*8 S0,S0A,S0B
      REAL STMP
c
c     The source function needs to maintain some data
c     so that the integrations can go more efficiently
c     Unfortunately ... since the integration routines
c     are generic, this data needs to be maintained in the
c     source functions themselves. This information
c     include the last integrated value. The inital
c     values for each of these would be set in the
c     setupval subroutine.
c
c     Functions:
c        ind = 1 : check for whole source integration
c        ind = 2 : check for S > last S and return partial
c                  integral
c        ind = 3 : update integrated data
c
c     This could be done utilizing the ENTRY statement to
c     allow multiple entry points to the same subroutine.
c     OR ... one could write separete routines for each
c     source function to do each item ... unfortunately this
c     would make using a generic intergration routine quite
c     difficult.
c
C
      IF (IND.EQ.1.AND.SLIM.Gt.LENSRC) THEN
         ind = -1
         if (plateopt.eq.0) then
            SRCION = IONINT
         elseif (plateopt.eq.1) then
            srcion = ioninti
         endif
         return
      endif
C
      if (ind.eq.2.and.slim.gt.sionl.and.
     >   plateopt.eq.pionl) then
         ind = -2
         srcion = iionl
         slim = sionl
         return
      endif
C
      if (ind.eq.3) then
         iionl = slim
         pionl = plateopt
         sionl = s
         srcion = iionl
         return
      endif
C
      IF (S.LE.LENSRC) THEN
        IF (SOPT.EQ.0.OR.SOPT.EQ.1) THEN
          IF (PLATEOPT.EQ.0 ) THEN
            S0 = S0OUT
          ELSE
            S0 = S0IN
          ENDIF
          IF (SOPT.EQ.0) THEN
            SRCION = S0
          ELSEIF(SOPT.EQ.1) THEN
            SRCION = S0 * EXP(-S/LAMSRC)
          ENDIF
        ELSEIF (SOPT.EQ.2.OR.SOPT.EQ.3) THEN
C
C         SOPT 2 AND 3 DIFFER IN THE VALUE OF FNORM
C
          IF (PLATEOPT.EQ.0) THEN
            STMP = S
            IN = IPOS(STMP,TSS,MAXIK)
            if (in.eq.1) then
               SRCION =  FNORM *  PINION(IN,IRN)
            else
               SRCION =  FNORM * ( PINION(IN-1,IRN) +
     >              (PINION(IN,IRN)-PINION(IN-1,IRN))/
     >              (TSS(IN)-TSS(IN-1))*
     >              (S-TSS(IN-1)))
            endif
C
C            WRITE(6,*) IN,IRN,S,TSS(IN),TSS(IN-1),FNORM,SRCION,
C     >                  PINION(IN,IRN),PINION(IN-1,IRN)
C
C
          ELSEIF (PLATEOPT.EQ.1) THEN
            STMP = TMAXS-S
            IN = IPOS(STMP,TSS,MAXIK)
            if (in.eq.maxik.and.stmp.gt.tss(maxik)) then
               SRCION =  FNORMI *  PINION(IN,IRN)
            else
               SRCION =  FNORMI * ( PINION(IN,IRN) +
     >              (PINION(IN-1,IRN)-PINION(IN,IRN))/
     >              (TSS(IN)-TSS(IN-1))*
     >              (TSS(IN)-STMP))
             endif
          ENDIF
        ELSEIF (SOPT.EQ.4.OR.SOPT.EQ.5) THEN
          IF (PLATEOPT.EQ.0 ) THEN
            S0A = S0AOUT
            S0B = S0BOUT
          ELSE
            S0A = S0AIN
            S0B = S0BIN
          ENDIF
          IF (SOPT.EQ.4) THEN
            IF (S.LE.LAMSRC) THEN
              SRCION = S0A
            ELSE
              SRCION = S0B
            ENDIF
          ELSEIF(SOPT.EQ.5) THEN
            SRCION = S0A * EXP(-S/LAMSRC) + S0B
          ENDIF
        ENDIF
      ELSE
        SRCION = 0.0D0
      ENDIF
C
      RETURN
      END
C
C
      DOUBLE PRECISION FUNCTION SRCRAD(S,POPT,PLATEOPT,SLIM,IND)
      use mod_params
      use mod_pindata
      use mod_comsol
      IMPLICIT NONE
      DOUBLE PRECISION S,SLIM
      INTEGER POPT,PLATEOPT,IND
C
C     THIS FUNCTION RETURNS THE VALUE OF THE RADIATION SOURCE
C     STRENGTH AT THE POSITION S FROM ONE OF TWO PLATES.
C     THIS IS NEEDED TO ENHANCE THE ACCURACY OF THE INTEGRATION
C     ROUTINES - WHICH ARE SOMEWHAT INACCURATE USING THE CURRENT
C     METHOD - HOWEVER - THIS ALTERNATE COULD BE COMPUTATIONALLY
C     EXPENSIVE AND SO THE ALTERNATE CODE WILL BE RETAINED.
C
C     DAVID ELDER    MAY 1, 1992
C
C     INCLUDE "PARAMS"
c     include 'params'
C     INCLUDE "PINDATA"
c     include 'pindata'
C     INCLUDE "COMSOL"
c     include 'comsol'
C
      REAL*8 P0
c     The source function needs to maintain some data
c     so that the integrations can go more efficiently
c     Unfortunately ... since the integration routines
c     are generic, this data needs to be maintained in the
c     source functions themselves. This information
c     include the last integrated value. The inital
c     values for each of these would be set in the
c     setupval subroutine.
c
c     Functions:
c        ind = 1 : check for whole source integration
c        ind = 2 : check for S > last S and return partial
c                  integral
c        ind = 3 : update integrated data
c
c     This could be done utilizing the ENTRY statement to
c     allow multiple entry points to the same subroutine.
c     OR ... one could write separete routines for each
c     source function to do each item ... unfortunately this
c     would make using a generic intergration routine quite
c     difficult.
c
C
      if (ind.eq.1.and.slim.gt.plensrc) then
         ind = -1
         if (plateopt.eq.0) then
            srcrad = radint
         elseif (plateopt.eq.1) then
            srcrad = radinti
         endif
         return
      endif
C
      if (ind.eq.2.and.slim.gt.sradl.and.
     >   plateopt.eq.pradl) then
         ind = -2
         srcrad = iradl
         slim = sradl
         return
      endif
C
      if (ind.eq.3) then
         iradl = slim
         srcrad = iradl
         pradl = plateopt
         sradl = s
         return
      endif
C
      IF (S.LE.PLENSRC) THEN
        IF (PLATEOPT.EQ.0 ) THEN
          P0 = P0OUT
        ELSE
          P0 = P0IN
        ENDIF
        IF (POPT.EQ.0.OR.POPT.EQ.2) THEN
          SRCRAD = P0
        ELSEIF (POPT.EQ.1.OR.POPT.EQ.3) THEN
          SRCRAD = P0 * EXP(-S/PLAMSRC)
        ENDIF
      ELSE
        SRCRAD = 0.0D0
      ENDIF
C
      RETURN
      END
C
C
C
C
C
      DOUBLE PRECISION FUNCTION SRCPEI(S,ikn,PLATEOPT,ds,IND)
      use mod_params
      use mod_pindata
      use mod_comtor
      use mod_comsol
      use mod_comhr
      IMPLICIT NONE
      DOUBLE PRECISION S,ds
      INTEGER PLATEOPT,IND,ikn
C
C     THIS FUNCTION RETURNS THE VALUE OF THE RADIATION SOURCE
C     STRENGTH AT THE POSITION S FROM ONE OF TWO PLATES.
C     THIS IS NEEDED TO ENHANCE THE ACCURACY OF THE INTEGRATION
C     ROUTINES - WHICH ARE SOMEWHAT INACCURATE USING THE CURRENT
C     METHOD - HOWEVER - THIS ALTERNATE COULD BE COMPUTATIONALLY
C     EXPENSIVE AND SO THE ALTERNATE CODE WILL BE RETAINED.
C
C     DAVID ELDER    MAY 1, 1992
C
C     INCLUDE "PARAMS"
c     include 'params'
C     INCLUDE "PINDATA"
c     include 'pindata'
c
c     include 'comtor'
c
C     INCLUDE "COMSOL"
c     include 'comsol'
c     include 'comhr'
C
      REAL*8 P0
c     The source function needs to maintain some data
c     so that the integrations can go more efficiently
c     Unfortunately ... since the integration routines
c     are generic, this data needs to be maintained in the
c     source functions themselves. This information
c     include the last integrated value. The inital
c     values for each of these would be set in the
c     setupval subroutine.
c
c     Functions:
c        ind = 1 : initialize to 0.0
c
c     This could be done utilizing the ENTRY statement to
c     allow multiple entry points to the same subroutine.
c     OR ... one could write separete routines for each
c     source function to do each item ... unfortunately this
c     would make using a generic intergration routine quite
c     difficult.
c
C
      double precision teb,tib,neb,srctmp
      integer ilow,imid,ihigh
c
c      write(6,*) 'srcpei1:',ind,s,popt,slim,plateopt,
c     >   cind,maxind
c
      if (ind.eq.1) then
         ind = -1
         ipeil = 0.0
         srcpei = ipeil
         ppeil = plateopt
         speil = s
         return
      endif
c
C
      if (plateopt.eq.0) then
         teb = soltehr(ikn-1)
         tib = soltihr(ikn-1)
         neb = solnehr(ikn-1)
         srctmp = 1.14d-32*neb**2*(teb-tib)/(crmb*teb**(1.5))
         ipeil = ipeil + srctmp * ds
         speil = s
         srcpei = ipeil
      elseif (plateopt.eq.1) then
         teb = soltehr(ikn+1)
         tib = soltihr(ikn+1)
         neb = solnehr(ikn+1)
         srctmp = 1.14d-32*neb**2*(teb-tib)/(crmb*teb**(1.5))
         ipeil = ipeil + srctmp * ds
         speil = s
         srcpei = ipeil
      endif
      RETURN
      END
C
C
C
      DOUBLE PRECISION FUNCTION CIS1(S,FVAL,SOPT,PLATEOPT)
      IMPLICIT NONE
      INTEGER PLATEOPT,SOPT
      DOUBLE PRECISION S,SMAX,FVAL
      EXTERNAL FVAL
C
C     CIS1 - PERFORMS ONE DIMENSIONAL NUMERICAL INTEGRATION ON FUNCTION
C            FVAL. THE RANGE OF INTEGRATION IS FROM 0 TO S.
C
C            THIS IS CALCULATED USING THE FORMULA:
C
C            SIGMA(I=1TON) SI * DELTAS
C
C            WHERE N IS THE NTH BIN AND CORRESPONDS TO POSITION S.
C
C            SOPT IS THE SOURCE OPTION THAT NEEDS TO BE
C            PASSED TO THE ROUTINE FVAL - IN ORDER TO EVALUATE THE
C            VALUE OF EITHER IONIZATION OR RADIATION CORRECTLY.
C
C
C     NSVAL -  SPECIFIES THE NUMBER OF SEGMENTS WHICH 0 TO S WILL BE
C     DIVIDED INTO FOR THE INTEGRATION. CAN BE MADE ARBITRARILY LARGE
C     FOR INCREASED ACCURACY AT THE COST OF PERFORMANCE.
C
C     THE SOURCE FUNCTION IS KNOWN TO BE IDENTICALLY ZERO BEYOND THE
C     POINT LSIZ.
C
      INTEGER NSVAL
      PARAMETER (NSVAL=1000)
      DOUBLE PRECISION DS,REM,INT1,STMP,F1,F2,RESULT,tmpval
      double precision strt,sret
      INTEGER NS,START,STOP,STEP,I,IND
      INTEGER INTTYPE
C
C
C
c      write(6,*) 'cis1:',s,sopt,plateopt
      IF (S.EQ.0.0) THEN
         CIS1 = 0.0D0
         RETURN
      ENDIF
C
C     This call tests if the integration range is greater
C     than the source length ... if it is a pre-calculated
c     value for the whole source is assigned and returned
c     Otherwise it falls through to the whole integration routine.
C
      IND = 1
      CIS1 = FVAL(S,SOPT,PLATEOPT,S,IND)
      IF (IND.EQ.-1) RETURN
C
C     This call retrieves the previous value of the integral
c     if it is appropriate.
C
      IND = 2
      tmpval = fval(s,sopt,plateopt,sret,ind)
      if (ind.eq.-2) then
         strt = sret
      else
         strt = 0.0d0
         tmpval = 0.0d0
      endif
C
c      write(6,*) 'before q:',strt,s,tmpval
      INTTYPE = 1
      CALL QSIMP(FVAL,strt,S,RESULT, SOPT,PLATEOPT,INTTYPE)
      CIS1 = RESULT + tmpval
c
c     Update the record in the source function with the
c     value of the integral up to this point.
c
      tmpval = result + tmpval
      ind = 3
      tmpval = FVAL(S,SOPT,PLATEOPT,tmpval,IND)
C
C      WRITE(6,*) 'COMPARISON:',RESULT,INT1
C
C      WRITE(6,*) 'LSIZ,STMP,DS:', LSIZ,STMP,DS,INT1
C
C
      RETURN
      END
C
C
C
      DOUBLE PRECISION FUNCTION CIS2(S,FVAL,SOPT,PLATEOPT)
      IMPLICIT NONE
      INTEGER SOPT,PLATEOPT
      DOUBLE PRECISION S,FVAL
      EXTERNAL FVAL
C
C
C     CIS2 - PERFORM A DOUBLE INTEGRATION OVER THE GIVEN SOURCE
C            FUNCTION.
C
C            SOPT IS THE SOURCE OPTION THAT NEEDS TO BE
C            PASSED TO THE ROUTINE FVAL - IN ORDER TO EVALUATE THE
C            VALUE OF EITHER IONIZATION OR RADIATION CORRECTLY.
C
C     NSVAL - SPECIFIES THE NUMBER OF PARTITIONS CONTRIBUTING TO THE
C     INTEGRATION.
C
      INTEGER NSVAL
      PARAMETER (NSVAL=1000)
      INTEGER I,NS,J,START,STOP,STEP
      DOUBLE PRECISION DS,INT2,INT1,STMP,DS2,F1,F2,RESULT
      INTEGER INTTYPE
C
C
      IF (S.EQ.0.0) THEN
         CIS2 = 0.0D0
         RETURN
      ENDIF
C
C
C     ADD CODE TO USE QUADRATURE FOR SECOND INTEGRAL AS WELL...
C
C     FOR N = 2    I = INT (X-T) F(T) DT  FOR T ON 0 TO X
C     SRCRAD RETURNS F(T) ... FUNCTION NEEDS TO BE (X-T)F(T)
C
C     NEED TO CREATE A FUNCTION WHICH RETURNS (X-T)F(T) ... OR ANOTHER
C     OPTION TO SRCRAD - INTTYPE
C
      INTTYPE = 2
      CALL QSIMP(FVAL,0.0D0,S,RESULT, SOPT,PLATEOPT,INTTYPE)
      CIS2 = RESULT
C
C      WRITE(6,*) 'COMPARISON:',RESULT
C
C
C
      RETURN
      END
C
C
C
      SUBROUTINE SETUPVAL(SMAX,SOPT,POPT,FSRC,LNSRC,LMSRC,
     >           N0,V0,T0,N0I,V0I,T0I,RCF,RCFI,LSSIZ,LPSIZ,IR,
     >           PAOUT,PAIN)
C
      use mod_params
      use mod_comtor
      use mod_cgeom
      use mod_comsol
      IMPLICIT NONE
      INTEGER SOPT,POPT,IR
      DOUBLE PRECISION FSRC,LNSRC,LMSRC,SMAX
      DOUBLE PRECISION N0,V0,T0
      DOUBLE PRECISION N0I,V0I,T0I
      DOUBLE PRECISION RCF,RCFI
      DOUBLE PRECISION LSSIZ,LPSIZ
      DOUBLE PRECISION PAOUT,PAIN
C
C     INCLUDE "PARAMS"
c     include 'params'
C     INCLUDE "COMTOR"
c     include 'comtor'
C     INCLUDE "CGEOM"
c     include 'cgeom'
C     INCLUDE "COMSOL"
c     include 'comsol'
C
C
C
C     SETUPVAL: THIS ROUTINE LOADS THE COMMON BLOCK VALUES NEEDED TO CAL
C               THE IONIZATION AND RADIATION SOURCE VALUES. THIS ALSO IN
C               THE EFFECT OF FLOW RECIRCULATION - IF ANY.
C
C     SOPT      - SOURCE IONIZATION OPTION
C     POPT      - SOURCE RADIATION OPTION
C     N0        - DENSITY AT PLATE (FOR THIS RING)
C     V0        - VELOCITY AT PLATE (FOR THIS RING)
C     T0        - TEMPERATURE AT PLATE (FOR THIS RING)
C     CSOLLS    - DECAY COEFFICIENT FOR SOURCE
C     CSOLPR    - INITIAL RADIATION COEFFICIENT
C     CSOLLR    - DECAY COEFFICIENT FOR RADIATION
C
C
C     DAVID ELDER   OCT  1991
C
C
C     OTHER DECLARATIONS
C
      INTEGER I,J,K,NS
      DOUBLE PRECISION TVAL,S0,STEMP,DS
      DOUBLE PRECISION CIS1,SRCION,srcrad
      double precision cisint
      EXTERNAL CIS1,SRCION,srcrad
C
C
C
c      write(6,*) 'setup:start:',SMAX,SOPT,POPT,FSRC,LNSRC,LMSRC,
c     >           N0,V0,T0,N0I,V0I,T0I,RCF,RCFI,LSSIZ,LPSIZ,IR,
c     >           PAOUT,PAIN
c
C
C     SET UP IONIZATION ARRAY FROM 0 TO SMAX FOR SPECIFIC RING
C     USING THE GIVEN OPTION FLAG
C
C
      IF (SOPT.EQ.0) THEN
         LENSRC = LNSRC * SMAX
         LAMSRC = LMSRC * SMAX
         S0OUT = -FSRC*N0*V0/LENSRC
         S0IN  = -FSRC*N0I*V0I/LENSRC
      ELSEIF (SOPT.EQ.1) THEN
         LENSRC = LNSRC * SMAX
         LAMSRC = LMSRC * SMAX
         S0OUT = -FSRC*N0*V0/(LAMSRC*
     >             (1-EXP(-LENSRC/LAMSRC)))
         S0IN  = -FSRC*N0I*V0I/(LAMSRC*
     >             (1-EXP(-LENSRC/LAMSRC)))
      ELSEIF (SOPT.EQ.2) THEN
         LENSRC = SMAX /2.0
         LAMSRC = SMAX /2.0
         IRN = IR
         MAXIK = NKS(IR)
         TMAXS = SMAX
         DO 500 I = 1,MAXIK
           TSS(I) = KSS(I,IR)
500      CONTINUE
         FNORM = 1.0
         FNORMI = 1.0
      ELSEIF (SOPT.EQ.3) THEN
         LENSRC = SMAX /2.0
         LAMSRC = SMAX /2.0
         IRN = IR
         MAXIK = NKS(IR)
         TMAXS = SMAX
         DO 600 I = 1,MAXIK
           TSS(I) = KSS(I,IR)
600      CONTINUE
C
C        FNORM AND FNORMI MUST BE SET TO 1 BEFORE THE CALL TO CIS1
C        BECAUSE THEY ARE USED WHEN THAT ROUTINE CALLS SRCION.
C
         FNORM = 1.0
         FNORMI = 1.0
C
C        PIN DOESN'T CALCULATE ON VIRTUAL RINGS, SO IF THEY ARE TO
C        USED IN DIVIMP, THE INTEGRAL OVER THE IONISATION SOURCE
C        MAY BE ZERO!
C
         CISINT = CIS1(SMAX/2.0,SRCION,SOPT,0)
         IF (CISINT.GT.0.0) THEN
            FNORM =  ABS(N0*V0 / CISINT)
         ELSE
            FNORM =  1.0
         ENDIF
         CISINT = CIS1(SMAX/2.0,SRCION,SOPT,1)
         IF (CISINT.GT.0.0) THEN
            FNORMI =  ABS(N0I*V0I / CISINT)
         ELSE
            FNORMI =  1.0
         ENDIF
      ELSEIF (SOPT.EQ.4) THEN
         LENSRC = SMAX* LNSRC
         LAMSRC = SMAX* LMSRC
         S0AOUT = -(1.0-CFIZ)*FSRC*N0*V0/LAMSRC
     >                 - FSRC*CFIZ*N0*V0/LENSRC
         S0BOUT = -FSRC*CFIZ*N0*V0/LENSRC
         S0AIN  = -(1.0-CFIZ)*FSRC*N0I*V0I/LAMSRC
     >                 - FSRC*CFIZ*N0I*V0I/LENSRC
         S0BIN = -FSRC*CFIZ*N0I*V0I/LENSRC
      ELSEIF (SOPT.EQ.5) THEN
         LENSRC = LNSRC * SMAX
         LAMSRC = LMSRC * SMAX
         S0AOUT = -FSRC*(1.0-CFIZ)*N0*V0/(LAMSRC*
     >             (1.0-EXP(-LENSRC/LAMSRC)))
         S0BOUT = -FSRC*CFIZ*N0*V0/LENSRC
         S0AIN  = -FSRC*(1.0-CFIZ)*N0I*V0I/(LAMSRC*
     >             (1.0-EXP(-LENSRC/LAMSRC)))
         S0BIN = -FSRC*CFIZ*N0I*V0I/LENSRC
      ENDIF
      LSSIZ = LENSRC
      IF (POPT.EQ.0) THEN
         PLENSRC = CSOLLR * SMAX
         PLAMSRC = 0.0
         P0OUT = CSOLPR
         P0IN  = CSOLPR
      ELSEIF (POPT.EQ.1) THEN
         PLENSRC = SMAX/2.0
         PLAMSRC = CSOLLR * SMAX
         P0OUT = CSOLPR
         P0IN  = CSOLPR
      ELSEIF (POPT.EQ.2) THEN
         PLENSRC = CSOLLR * SMAX
         PLAMSRC = 0.0
         P0OUT =  CSOLFR * PAOUT / PLENSRC
         P0IN = CSOLFR * PAIN / PLENSRC
      ELSEIF (POPT.EQ.3) THEN
         PLENSRC = SMAX/2.0
         PLAMSRC = CSOLLR * SMAX
         P0OUT = CSOLFR * PAOUT /
     >           (PLAMSRC*(1.0D0-DEXP(-PLENSRC/PLAMSRC)))
         P0IN  = CSOLFR * PAIN /
     >           (PLAMSRC*(1.0D0-DEXP(-PLENSRC/PLAMSRC)))
      ENDIF
      LPSIZ = PLENSRC
C
C
C      WRITE(6,*) 'POPT:',POPT,PLENSRC,PLAMSRC,PAOUT,PAIN,
C     >            P0OUT,P0IN,CSOLFR
C
C     CALCULATE CROSS FIELD TERMS IF FLUX RECIRCULATION HAS BEEN
C     SPECIFIED. DO NOT CALCULATE IT IF A NORMALIZED IONIZATION SOURCE
C     HAS BEEN SPECIFIED.
C
      IF (FLUXROPT.EQ.1.AND.FSRC.NE.1.0.AND.SOPT.NE.3) THEN
        RCF = -(N0*V0 + CIS1(SMAX/2.0,SRCION,SOPT,0) )
     >        /(SMAX/2.0)
        RCFI = -(N0I*V0I + CIS1(SMAX/2.0,SRCION,SOPT,1))
     >        /(SMAX/2.0)
      ELSE
        RCF = 0.0
        RCFI = 0.0
      ENDIF
C
C     Calculate full source integrals
C
      ionint  = cis1(lensrc,srcion,sopt,0)
      ioninti = cis1(lensrc,srcion,sopt,1)
      radint  = cis1(plensrc,srcrad,popt,0)
      radinti = cis1(plensrc,srcrad,popt,1)
c
c     Set-up the values needed for source function
c     integration.
c
c     Ion source
c
      sionl = 0.0d0
      pionl = 0
      iionl = 0.0d0
c
c     Rad source
c
      sradl = 0.0d0
      pradl = 0
      iradl = 0.0d0
C
C
C
c      write(6,*) 'setup:end:',SMAX,SOPT,POPT,FSRC,LNSRC,LMSRC,
c     >           N0,V0,T0,N0I,V0I,T0I,RCF,RCFI,LSSIZ,LPSIZ,IR,
c     >           PAOUT,PAIN
c      write(6,*) 'setup:end2:',p0out,p0in,plensrc,plamsrc
C
C
      RETURN
      END
C
C
      SUBROUTINE FLUXRLOOK(IR,FSRC,LMSRC,LNSRC)
      use mod_params
      use mod_comtor
      implicit none
      INTEGER IR
      DOUBLE PRECISION FSRC,LMSRC,LNSRC
C
C     INCLUDE "PARAMS"
c     include 'params'
C     INCLUDE "COMTOR"
c     include 'comtor'
C
C     IF THE FLUX RECIRCULATION SOURCE OPTION HAS BEEN SPECIFIED THIS
C     ROUTINE LOOKS TO SEE IF THE SOURCE STRENGTH MULTIPLIER, FSRC, THE
C     SOURCE DECAY LENGTH (FOR EXPONENTIAL), LMSRC, AND THE SOURCE
C     OVERALL LENGTH, LNSRC, HAVE BEEN SPECIFED FOR THIS RING.
C     IF THEY HAVE NOT BEEN THE VALUES PASSED IN REMAIN UNCHANGED -
C     OTHERWISE THEY ARE SET TO THE VALUES ENTERED IN THE INPUT FILE.
C
C     DAVID ELDER  MAY 1, 1992
C
      INTEGER I,J
C
C      WRITE(6,*) 'FLUXP:',FLUXPTS
C
C      DO 5 I = 1, FLUXPTS
C        WRITE(6,*) I,FLUXINFO(I,1),FLUXINFO(I,2),
C     >  FLUXINFO(I,3),FLUXINFO(I,4)
C 5    CONTINUE
C
      DO 10 I = 1, FLUXPTS
        IF (INT(FLUXINFO(I,1)) .EQ.IR) THEN
          FSRC = FLUXINFO(I,2)
          LNSRC = FLUXINFO(I,3)
          LMSRC = FLUXINFO(I,4)
          GOTO 20
        ENDIF
 10   CONTINUE
 
 20   CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE QSIMP(FUNC,A,B,S,OPT,PLATE,ITYP)
      use mod_params
      use mod_slcom
      implicit none
c slmod begin
c *TEMP*
c     INCLUDE 'params'
c     INCLUDE 'slcom'
c slmod end
      INTEGER JMAX,OPT,PLATE,ITYP
      DOUBLE PRECISION A,B,FUNC,S,EPS
      EXTERNAL FUNC
      PARAMETER (EPS=1.0D-5,JMAX=22)
C
C     QSIMP: THIS ROUTINE IS TAKEN FROM THE TEXT NUMERICAL RECIPES
C     IN FORTRAN. IT IMPLEMENTS SIMPSON'S METHOD OF CALCULATING
C     NUMERICAL QUADRATURE. IT CALLS THE ROUTINE TRAPZD TO REFINE
C     THE VALUE OF THE NUMERICAL INTEGRAL.
C
C     RETURNS AS S THE INTEGRAL OF THE FUNCTION FUNC FROM A TO B.
C     THE PARAMETER EPS CAN BE SET TO THE DESIRED FRACTIONAL ACCURACY
C     AND JMAX SO THAT 2 TO THE POWER JMAX-1 IS THE MAXIMUM ALLOWED NUMB
C     OF STEPS. INTEGRATION IS PERFORMED BY SIMPSONS RULE.
C
C     There are two types of integration performed by these
C     routines and the choice is controlled by the ityp
C     argument.
C
C     ITYP=1 : standard integration of f(x) from a to b
C     ITYP=2 : integration of (b-x)f(x) from a to b
C
      INTEGER J
      DOUBLE PRECISION OS,OST,ST
      OST = -1.0D30
      OS  = -1.0D30
      DO 10 J = 1,JMAX
        CALL TRAPZD (FUNC,A,B,ST,J,OPT,PLATE,ITYP)
        S = (4.0*ST-OST)/3.0
        IF (ABS(S-OS).LE.(EPS*ABS(OS))) RETURN
C
C
        OS = S
        OST = ST
C
C        WRITE(6,*) 'QSIMP:',OST,OS,ST
C
10    CONTINUE
C
C     ERROR CONDITION - INCOMPLETE CONVERGENCE
C
      WRITE(6,*) 'QSIMP: ERROR IN QUADRATURE',A,B,ST
      WRITE(6,*) 'QSIMP:',OST,OS,ST ,S
      write(6,*) 'ALSO :',A,B,J,OPT,PLATE,ITYP
c slmod begin
      IF (osmns28.NE.0.OR.iflexopt(8).EQ.11) THEN
c...                      980116023:
        S = 0.0D0
        WRITE(6,*) 'SOL28 SO NOT HALTING CODE'
        WRITE(0,*) 'QSIMP ERROR'
        RETURN
      ELSE
        STOP
      ENDIF
c
c      STOP
c slmod end
C      RETURN
      END
C
C
C
      SUBROUTINE TRAPZD (FUNC,A,B,S,N,OPT,PLATE,ITYP)
      implicit none
      INTEGER N,OPT,PLATE,ITYP
      DOUBLE PRECISION A,B,S,FUNC
      EXTERNAL FUNC
C
C     TRAPZD: THIS SUBROUTINE IS TAKEN FROM THE TEXT NUMERICAL RECIPES
C     IN FORTRAN. IT IS PART OF A SET OF ROUTINES FOR CALCULATING
C     NUMERICAL QUADRATURE OF A DISCRETE FUNCTION.
C
C     THIS ROUTINE COMPUTES THE NTH STAGE REFINEMENT OF AN EXTENDED
C     TRAPEZOIDAL RULE. FUNC IS INPUT AS THE NAME OF THE FUNCTION TO
C     BE INTEGRATED BETWEEN THE LIMITS OF A AND B, ALSO INPUT. WHEN
C     CALLED WITH N=1, THE ROUTINE RETURNS THE CRUDEST ESTIMATE OF
C     THE INTEGRAL. SUBSEQUENT CALLS WITH N=2,3... (IN THAT
C     SEQUENTIAL ORDER) WILL IMPROVE THE ACCURACY OF S BY ADDING
C     2**(N-2) ADDITIONAL INTERIOR POINTS. S SHOULD NOT BE MODIFIED
C     BETWEEN SEQUENTIAL CALLS.
C
      INTEGER IT,J,ind,dummy
      DOUBLE PRECISION DEL,SUM,TNM,X
c
c     Set ind = -1 so that the functions that return the
c     value of the source function will not return the
c     integrated value if S > Src length is passed. This is
c     necessary for the 2D integration through CIS2.
c
c     This function is also performed by setting slim = -1.0.
c
 
      ind = 0
      dummy = 0
      IF (N.EQ.1) THEN
         IF (ITYP.EQ.1) THEN
            S = 0.5*(B-A)*(FUNC(A,OPT,PLATE,-1.0d0,ind)+
     >           FUNC(B,OPT,PLATE,-1.0d0,ind))
         ELSEIF (ITYP.EQ.2) THEN
            S = 0.5*(B-A)*( (B-A)*FUNC(A,OPT,PLATE,-1.0d0,ind))
         ENDIF
C
C         WRITE(6,*) 'TRAPZD:',ITYP,A,B,
C     >                FUNC(A,OPT,PLATE,-1.0d0,ind),
C     >                FUNC(B,OPT,PLATE,-1.0d0,ind)
C
      ELSE
         IT = 2**(N-2)
         TNM = IT
         DEL = (B-A)/TNM
         X = A + 0.5*DEL
         SUM = 0.0
         DO 10 J = 1,IT
            IF (ITYP.EQ.1) THEN
               SUM = SUM + FUNC(X,OPT,PLATE,-1.0d0,ind)
            ELSEIF (ITYP.EQ.2) THEN
               SUM = SUM + (B-X) * FUNC(X,OPT,PLATE,-1.0d0,ind)
            ENDIF
            X = X + DEL
 10      CONTINUE
         S = 0.5*(S+(B-A)*SUM/TNM)
C
C         WRITE (6,*) 'TRAPZ:',S,SUM
C
      ENDIF
      RETURN
      END
c
c
      subroutine calcnv(te,ti,gamma,pinf,n,v)
      use mod_params
      use mod_comtor
      implicit none
      double precision te,ti,gamma,pinf,n,v
C     INCLUDE   "PARAMS"
c     include 'params'
C     INCLUDE   "CGEOM"
C     include 'cgeom'
C     INCLUDE   "COMTOR"
c     include 'comtor'
c
c     The purpose of this routine is to extract the
c     code that calculates the n,v values since it is
c     virtually identical for the SOL options avaliable.
c     This is due to the fact that n,v are calculated
c     from the pressure balance equation and the
c     flow equation after Te and Ti have been found.
c     Since the pressure balance and flow are unaffected
c     by additional terms in the equations for Te,Ti ...
c     It seems likely that this code can stay as it
c     is for some time.
c
c     Note: v is returned as -ve only and needs to be
c     adjusted depending on which side it is on.
c
c
      double precision rootn,massi
c
      MASSI = CRMB * AMU
C
C     CALCULATE DENSITY
C
      ROOTN = (PINF/(ECH*(te+ti)))**2
     >        - 4.0* (MASSI*GAMMA**2)
     >        / (ECH*(te+ti))
C
C
      IF (SROOTOPT.EQ.0.OR.
     >   (SROOTOPT.EQ.1.AND.ROOTN.GE.0.0)) THEN
C
         IF (ROOTN.LT.0.0) ROOTN = 0.0
C
C        Calculate density
C
         n = PINF/(2.0*ECH*(te+ti)) + 0.5 * SQRT(ROOTN)
C
C                 CALCULATE VELOCITY
C
         v = GAMMA/n
C
      ELSEIF (SROOTOPT.EQ.1) THEN
C
C        SET VELOCITY TO LOCAL SOUND SPEED.
C
         v = -SQRT(0.5*EMI*(te+ti)*(1+RIZB)/CRMB)
C
C        SET DENSITY
C
         n = ABS(GAMMA / v)
C
      ENDIF

c      write(6,'(a,10(1x,g12.5))') 'SOLEDGE:NV:',n,v,rootn,te,ti,
c     > gamma,pinf

      return
      end
c
c
c
      subroutine specplas(irlim1,irlim2,ikopt)
      use mod_params
      use mod_cgeom
      use mod_comtor
      implicit none
      integer irlim1,irlim2,ikopt
c
c     include 'params'
c     include 'cgeom'
c     include 'comtor'
c
c     SPECPLAS:
c
c     This subroutine applies a two-part linear
c     specified set of values for Te, Ti, n and v
c     for the background characteristics for the 
c     rings irlim1 to irlim2. 
c
c     For each of the quantities (Te, Ti ...) there
c     are four parameters specified.
c
c     e.g.
c
c     cteS1 - Te Length 1 - cteS1 * SMAX = first length
c     cteF1 - Te Fact   1 - cteF1 * Te(target) = Te at first length
c     cteS2 - Te Length 2 - cteS1 * SMAX = first length
c     cteF2 - Te FAct   2 - cteF1 * Te(target) = Te at second length
c 
c     The quantity is constant at Te(target) * cteM2 from teh second
c     point to the middle of the ring. It is linealy interpolated 
c     between each of the points - starting at Te(target) at S=0.
c
c     This routine may later allow for a different set of parameters
c     to be specified for both SOL and Private Plasma - though
c     it initially will use only one set of parameters. 
c    
c
      integer in,ir,ik,ikmid,ikstart,ikend
      real te0,ti0,ne0,vb0,smax
      real te1,te2,tes1,tes2
      real ti1,ti2,tis1,tis2
      real ne1,ne2,nes1,nes2
      real vb1,vb2,vbs1,vbs2
      real s 
c
      real lin_inter2
      external lin_inter2
c
c      write (6,*) 'SPEC PLAS:',irlim1,irlim2 
c
      do ir = irlim1,irlim2
c
         call set_ikvals(ir,ikstart,ikend,ikopt)
c
         smax = ksmaxs(ir)
c
C        Find ring midpoint
c
         ikmid = ikmids(ir) +1
c
c        Set up lengths - only needs to be done ONCE
c
         nes1 = cnes1 * smax
         nes2 = cnes2 * smax
         tes1 = ctes1 * smax
         tes2 = ctes2 * smax
         tis1 = ctis1 * smax
         tis2 = ctis2 * smax
         vbs1 = cvbs1 * smax
         vbs2 = cvbs2 * smax

c
         if (ikopt.eq.1.or.ikopt.eq.3) then  
c
c
c        Set up parameters for Outer half
c
         ne0 = KNBS(1,IR)
         te0 = KTEBS(1,IR)
         ti0 = KTIBS(1,IR)
         Vb0 = - SQRT(0.5*EMI*(TE0+TI0)*(1+RIZB)/CRMB)
c
         ne1 = cnef1 * ne0
         ne2 = cnef2 * ne0
         te1 = ctef1 * te0
         te2 = ctef2 * te0
         ti1 = ctif1 * ti0
         ti2 = ctif2 * ti0
         vb1 = cvbf1 * vb0
         vb2 = cvbf2 * vb0
c
c        Store TARGET values
c
         KNDS(IDDS(IR,2))  = Ne0
         KVDS(idds(ir,2))  = Vb0
         KTEDS(IDDS(IR,2)) = Te0
         KTIDS(IDDS(IR,2)) = ti0
c
c        Do OUTER half 
c
         do ik = 1,ikmid-1
c
            s = kss(ik,ir)
c
            ktebs(ik,ir) = lin_inter2(te0,te1,tes1,te2,tes2,s)            
c 
            ktibs(ik,ir) = lin_inter2(ti0,ti1,tis1,ti2,tis2,s)            
c 
            knbs(ik,ir) = lin_inter2(ne0,ne1,nes1,ne2,nes2,s)            
c 
            kvhs(ik,ir) = lin_inter2(vb0,vb1,vbs1,vb2,vbs2,s)            
c 
            kes(ik,ir) = 0.0
c
         end do
c
         endif
c
c
         if (ikopt.eq.2.or.ikopt.eq.3) then 
c
c        DO the INNER half now 
c
c        Set up parameters for INNER half
c
         ne0 = KNBS(nks(ir),IR)
         te0 = KTEBS(nks(ir),IR)
         ti0 = KTIBS(nks(ir),IR)
         Vb0 = SQRT(0.5*EMI*(TE0+TI0)*(1+RIZB)/CRMB)
c
         ne1 = cnef1 * ne0
         ne2 = cnef2 * ne0
         te1 = ctef1 * te0
         te2 = ctef2 * te0
         ti1 = ctif1 * ti0
         ti2 = ctif2 * ti0
         vb1 = cvbf1 * vb0
         vb2 = cvbf2 * vb0
c
c        Store the TARGET conditions  
c
         KNDS(IDDS(IR,1))  = Ne0
         KVDS(idds(ir,1))  = Vb0
         KTEDS(IDDS(IR,1)) = TE0
         KTIDS(IDDS(IR,1)) = ti0
c
c        Do INNER half 
c
         do ik = nks(ir),ikmid,-1
c
            s = smax - kss(ik,ir)
c
            ktebs(ik,ir) = lin_inter2(te0,te1,tes1,te2,tes2,s)            
c 
            ktibs(ik,ir) = lin_inter2(ti0,ti1,tis1,ti2,tis2,s)            
c 
            knbs(ik,ir) = lin_inter2(ne0,ne1,nes1,ne2,nes2,s)            
c 
            kvhs(ik,ir) = lin_inter2(vb0,vb1,vbs1,vb2,vbs2,s)            
c 
            kes(ik,ir) = 0.0
c
c
         end do
c
         endif
c
      end do  
c
      return
      end 
c
c
c
      real function lin_inter2(f0,f1,s1,f2,s2,s)
      implicit none
      real f0,f1,f2
      real s,s1,s2
c    
c     LIN_INTER2:
c    
c     This subroutine lineraly interpolates between three
c     points:
c     
c     S-value     Function-Value
c    
c       0.0           F0
c        S1           F1
c        S2           F2
c    
c        if (s>s2) then Freturned = F2 
c    
      real lin_inter1
      external lin_inter1 
c
      if (s.lt.s1) then 
c     
         lin_inter2 = lin_inter1(0.0,f0,s1,f1,s)
c     
      elseif (s.lt.s2) then
c     
         lin_inter2 = lin_inter1(s1,f1,s2,f2,s)
c     
      elseif (s.ge.s2) then 
c     
         lin_inter2 = f2
c      
      endif
c     
      return
      end
c
c
c
      real function lin_inter1(s1,f1,s2,f2,s)
      implicit none
      real f1,f2
      real s,s1,s2
c     
c     LIN_INTER1:
c     
c     This subroutine lineraly interpolates between two
c     points:
c     
c     S-value     Function-Value
c     
c        S1           F1
c        S2           F2
c     
      lin_inter1 = (s-s1)/(s2-s1) * (f2-f1) + f1  
c
      return
      end
c
c
c
      subroutine cnvrtptos(p,s,ir)
      use mod_params
      use mod_cgeom
      implicit none
      double precision p,s
      integer ir
c
c     Common blocks
c
c     include 'params'
c     include 'cgeom'
c
c     CNVRTPTOS: Finds the S value along the field line
c                that approximately corresponds to the given
c                P value. It interpolates between cell centers
c                after finding the correct one. 
c
      integer ik
c
      do ik = 1,nks(ir)-1
c
         if (ik.eq.1.and.p.lt.kps(ik,ir)) then 
            s = kss(ik,ir) * p / kps(ik,ir)
c            exit
            goto 100
         elseif (ik.eq.nks(ir)-1.and.p.gt.kps(ik+1,ir)) then             
            s = kss(ik+1,ir) + (ksmaxs(ir)-kss(ik+1,ir)) * 
     >             (p-kps(ik+1,ir))/(kpmaxs(ir)-kps(ik+1,ir)) 
c            exit 
            goto 100
         elseif (p.ge.kps(ik,ir).and.p.le.kps(ik+1,ir)) then  
            s = kss(ik,ir) + (kss(ik+1,ir)-kss(ik,ir)) * 
     >             (p-kps(ik,ir))/(kps(ik+1,ir)-kps(ik,ir)) 
c            exit 
            goto 100
         endif
c
      end do  
c
 100  continue      
c
c     Return 
c
      return
      end
c
c
c
      subroutine load_s21params(ir,l1r,l1ri,l2r,l2ri,
     >                         lvr,lvri,ter,teri,tir,tiri,nr,nri,
     >                         vbm,vbmi,qr,qri,n_exp,n_expi,
     >                         l1a,nr1a,l1ai,nr1ai,
     >                         l1b,nr1b,l1bi,nr1bi,
     >                         ter1a,ter1ai,tir1a,tir1ai,
     >                         ter1b,ter1bi,tir1b,tir1bi)
c
      use mod_params
      use mod_comtor
      implicit none
c     include 'params'
c     include 'comtor' 
c
      integer ir
      double precision l1r,l1ri,l2r,l2ri,lvr,lvri,
     >                 ter,teri,tir,tiri,nr,nri,vbm,vbmi,qr,qri,
     >                 n_exp,n_expi
c
c     Extra parameters - have no effect unless specified in 
c                        additional per ring data sets.
c
c     l1ra ... allow for more detail in the description of the 
c              density evolution in region A - three linear 
c              fitted regions. 
c
      double precision l1a,l1ai,l1b,l1bi,nr1a,nr1ai,nr1b,nr1bi,
     >                         ter1a,ter1ai,tir1a,tir1ai,
     >                         ter1b,ter1bi,tir1b,tir1bi
c
c     LOAD_S21PARAMS: This routine loads the inner and outer SOL 21
c     parameters for the specified ring.
c
c   
      integer ringno,ierr,in
      external ringno
c
c     Load the DEFAULT values to start and replace as necessary
c
c     Te ratio
c
      ter  = terat 
      teri = terati
c
c     Ti ratio
c
      tir  = tirat 
      tiri = tirati
c
c     Ne ratio at L1A
c
      nr1a  = 1.0
      nr1ai = 1.0
c
c     Te/Ti ratio at L1A
c
      ter1a  = 1.0
      ter1ai = 1.0
      tir1a  = 1.0
      tir1ai = 1.0
c
c     Ne ratio at L1B
c
      nr1b  = 1.0
      nr1bi = 1.0
c
c     Te/Ti ratio at L1B
c
      ter1b  = 1.0
      ter1bi = 1.0
      tir1b  = 1.0
      tir1bi = 1.0
c
c     Ne ratio at L1
c
      nr  = nrat 
      nri = nrati
c
c     Ne Exponent
c
      n_exp = nalph 
      n_expi= nalphi
c
c     Q ratio
c
      qr  = qrat 
      qri = qrati
c
c     VB multiplier
c
      vbm  = vbmult
      vbmi = vbmulti
c
c     L1A distance
c
      l1a  = 0.0
      l1ai = 0.0
c
c     L1B distance
c
      l1b  = 0.0
      l1bi = 0.0
c
c     L1 distance
c
      l1r  = l1rat
      l1ri = l1rati
c
c     L2 distance
c
      l2r  = l2rat
      l2ri = l2rati
c
c     LV distance
c
      lvr  = lvrat
      lvri = lvrati

c
c     Check for over-riding data on INNER
c
      if (ns21i.gt.0) then 
c
         in = RINGNO(ir,s21parmi,ns21i,maxnrs,9,ierr)
c
         if (ierr.eq.0) then
c
            teri  = s21parmi(in,2)
            tiri  = s21parmi(in,3)
            nri   = s21parmi(in,4)
            n_expi= s21parmi(in,5)
            qri   = s21parmi(in,6)
            l1ri  = s21parmi(in,7)
            l2ri  = s21parmi(in,8)
            lvri  = s21parmi(in,9)
            vbmi  = s21parmi(in,10)
c
         endif
c
      endif    
c
c     Check for over-riding data on OUTER
c
      if (ns21o.gt.0) then 
c
         in = RINGNO(ir,s21parmo,ns21o,maxnrs,9,ierr)
c
         if (ierr.eq.0) then
c
            ter  = s21parmo(in,2)
            tir  = s21parmo(in,3)
            nr   = s21parmo(in,4)
            n_exp= s21parmo(in,5)
            qr   = s21parmo(in,6)
            l1r  = s21parmo(in,7)
            l2r  = s21parmo(in,8)
            lvr  = s21parmo(in,9)
            vbm  = s21parmo(in,10)
c
         endif
c
      endif    
c
c     Check for over-riding extra/auxilliary input data on INNER
c
      if (aux_ns21i.gt.0) then 
c
         in = RINGNO(ir,aux_s21parmi,aux_ns21i,maxnrs,9,ierr)
c
         if (ierr.eq.0) then
c
            l1ai    = aux_s21parmi(in,2)
            l1bi    = aux_s21parmi(in,3)
            nr1ai   = aux_s21parmi(in,4)
            nr1bi   = aux_s21parmi(in,5)
            ter1ai  = aux_s21parmi(in,6)
            ter1bi  = aux_s21parmi(in,7)
            tir1ai  = aux_s21parmi(in,8)
            tir1bi  = aux_s21parmi(in,9)
c
         endif
c
      endif    
c
c     Check for over-riding extra/auxilliary input data on OUTER
c
      if (aux_ns21o.gt.0) then 
c
         in = RINGNO(ir,aux_s21parmo,aux_ns21o,maxnrs,9,ierr)
c
         if (ierr.eq.0) then
c
            l1a     = aux_s21parmo(in,2)
            l1b     = aux_s21parmo(in,3)
            nr1a    = aux_s21parmo(in,4)
            nr1b    = aux_s21parmo(in,5)
            ter1a   = aux_s21parmo(in,6)
            ter1b   = aux_s21parmo(in,7)
            tir1a   = aux_s21parmo(in,8)
            tir1b   = aux_s21parmo(in,9)
c
         endif
c
      endif    
c
c     END of LOAD_S21PARAMS 
c
      return
      end
