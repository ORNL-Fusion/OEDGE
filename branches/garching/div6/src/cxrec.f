c     -*-Fortran-*-
C@PROCESS OPT(1),VECTOR(LEV(0))
C
C
      SUBROUTINE CXREC (NIZS,CION,CIOPTI,RIZB,CRMB,CVCX,
     >                  CNHC,CNHO,CLAMHX,CLAMHY,cprint,cpinopt)
      IMPLICIT    NONE
      INTEGER     NIZS,CION,CIOPTI,cprint,cpinopt
      REAL        CRMB,CVCX,CNHC,CNHO,CLAMHX,CLAMHY,RIZB
C
C  *********************************************************************
C  *                                                                   *
C  *  CXREC: THIS ROUTINE DERIVED FROM NOTE 89                         *
C  *         VCX MODIFICATION GIVEN IN NOTE 173 APRIL 1988             *
C  *  DEALS WITH CHARGE EXCHANGE RECOMBINATION, IMPORTANT FOR          *
C  *  HYDROGENIC PLASMAS.                                              *
C  *                                                                   *
C  *  VECTORISATION SWITCHED OFF DUE TO POSSIBILITY OF DIVIDE BY ZERO  *
C  *  ERRORS WITHIN FINAL DO-LOOPS.  SHOULDN'T HAPPEN IN SCALAR SINCE  *
C  *  IF STATEMENTS TAKE CARE OF IT.  HOWEVER, THE IBM OPTIMISER       *
C  *  SHUFFLES THE LOOP AND THE ERRORS DO OCCUR.  HENCE THE SPECIAL    *
C  *  COMPILER OPTIONS ABOVE FOR THIS ROUTINE ONLY.                    *
C  *                                                                   *
C  *                                      C.M.FARRELL   JANUARY 1988   *
C  *                                                                   *
C  *********************************************************************
C
      include    'params'
      include    'cgeom'
      include    'cioniz'
      include    'pindata'
      include    'cedge2d'
      include    'cadas'
c
c     Temporary includes
c
      include 'temp'
C
      INTEGER IK,IR,IZ,in,ntemp
      REAL    VCX,V,Q(MAXIZS),SIGCX
      REAL    SMAX,S
c
      real sigvcx_maggi,sigvcx_maggi_mod
      external sigvcx_maggi,sigvcx_maggi_mod
c
c     Variables for helping with ADPAK/INEL CX data
c
      real    tmpte,tmpti,tmpne,rion,rrec,rcxr
c
c     neutral temperature added (may become a poloidal distribution
c     in future; Krieger, IPP/95)
      real    tbgr
c     functions to calculate boron IV CX; Krieger, IPP/95
      real    b4atotal
c
cfm
c      character*80  usercfm
cfm
C
C-----------------------------------------------------------------------
C     CALCULATE NEUTRAL HYDROGEN ATOM DENSITY  (FUNCTION OF X AND Y)
C     VARIOUS OPTIONS ALLOWED, KEYED WITH CIOPTI
C-----------------------------------------------------------------------
C
C  OPTION 0 - CXREC OFF - SET TO ZERO AND EXIT.
C  --------
      IF (CIOPTI.EQ.0) THEN
c
        CALL RZERO (KNHS, MAXNKS*MAXNRS)
c
c       Zero the CX recombination array and return to calling routine 
c
        call rzero (kfcxs,maxnks*maxnrs*(maxizs+1))
c
        return 

C
C  OPTION 1 AND 2
C  --------------
C
C  ALLOW BACKGROUND NEUTRAL DENSITY TO BE CONSTANT OVER ALL SPACE.
C
      ELSEIF (CIOPTI.EQ.1.OR.CIOPTI.EQ.2) THEN
        DO 110 IR = 1,NRS
          DO 110 IK = 1,NKS(IR)
            KNHS(IK,IR) = CNHO
 110    CONTINUE
C
C  OPTION 3
C  --------
C
C  BACKGROUND NEUTRAL DENSITY IS CONSTANT IN THE CORE AND DECAYS
C  EXPONENTIALLY FROM THE PLATES.
C
      ELSEIF (CIOPTI.EQ.3) THEN
C
C       CORE - CONSTANT
C
        DO 150 IR = 1,IRSEP-1
          DO 150 IK = 1,NKS(IR)
            KNHS(IK,IR) = CNHC
 150    CONTINUE
C
C       SOL AND TRAP - EXPONENTIAL DECAY FROM PLATES.
C
        DO 160 IR = IRSEP,NRS
          SMAX = KSMAXS(IR)
          DO 160 IK = 1,NKS(IR)
            IF (KSS(IK,IR).GT.0.5*SMAX) THEN
              S = SMAX - KSS(IK,IR)
            ELSE
              S = KSS(IK,IR)
            ENDIF
            KNHS(IK,IR) = CNHO * EXP(-S/CLAMHX)
 160    CONTINUE
C
C OPTION 4 5 6 7 AND 8 - setup
C --------------------------
C
C BACKGROUND NEUTRAL DENSITY IS READ FROM NIMBUS OR LOADED
C FROM A FLUID CODE RESULT.
C
      ELSE IF (CIOPTI.EQ.4.OR.CIOPTI.EQ.5.or.ciopti.eq.6.or.
c slmod begin
     >         ciopti.eq.7.or.ciopti.eq.8.or.ciopti.eq.9) THEN
c
c     >         ciopti.eq.7.or.ciopti.eq.8) THEN
c slmod end
c 
        if (cpinopt.eq.1) then 
           DO IR =1,NRS
              DO IK =1,NKS(IR)
                 KNHS(IK,IR)= PINATOM(IK,IR)
              end do
           end do  
        else
           DO IR =1,NRS
              DO IK =1,NKS(IR)
                 KNHS(IK,IR)= E2DATOM(IK,IR)
              end do
           end do  
        endif
c
c       added background temp., Krieger IPP/96
c
        tbgr = 5.0
c
      ENDIF
c
c     added background temp., Krieger IPP/96
c     (might become a 2D poloidal distribution)
c
c      DO IR =1,NRS                                                
c        DO IK =1,NKS(IR)                                           
c	  tbgr = 5.0
c        ENDDO
c      ENDDO
c
c
C-----------------------------------------------------------------------
c
c
c     ------------ CALCULATE CX TIMES ----------
c
c
C-----------------------------------------------------------------------
c     BASIC CX OPTIONS - 0 1 2 3 4 5 
C-----------------------------------------------------------------------
c
      if (ciopti.eq.0.or.ciopti.eq.1.or.ciopti.eq.2.or.
     >    ciopti.eq.3.or.ciopti.eq.4.or.ciopti.eq.5) then 
c
      DO 240 IR = 1, NRS
C
       DO 230 IK = 1, NKS(IR)
        IF (CIOPTI.EQ.0.OR.CIOPTI.EQ.1.OR.CIOPTI.EQ.3.OR.
     &      CIOPTI.EQ.4) THEN
          VCX = 1.56E4 * SQRT (KTIBS(IK,IR)/CRMB)
        ELSEIF (CIOPTI.EQ.2.OR.CIOPTI.EQ.5) THEN
          VCX = CVCX
        ENDIF
        V = 1.E-4 * VCX
C
c
C-----------------------------------------------------------------------
C       CALCULATE Q VALUES (FORMULAE DIFFERENT FOR EACH SPECIES)
c
C
        DO 200 IZ = 1, NIZS
          Q(IZ) = 0.0
  200   CONTINUE
C
C       HELIUM
C       ------
        IF     (CION.EQ.2) THEN
          Q(1) =-0.884 + 0.021*V
          Q(2) = 0.114 - 0.038*V + 0.001929*V*V
C
C       BERYLLIUM
C       ---------
        ELSEIF (CION.EQ.4) THEN
          Q(4) = 2.109 + 1.048*V - 0.00718*V*V

C                                                                       
C       BORON (added by Karl Krieger, IPP/96)                           
C       ----- there are only data available for BV->BIV CX rec.
C             since the CX-rate is a function of T_ion and T_neutral,
C             we include a "dummy" Q and do the real calculation
C             below
c
        ELSEIF (CION.EQ.5) THEN                                         
          Q(4) = 3.24-377.3*LOG(ABS(V))/V/V+276.6/V +0.32*V*EXP(-6.52/V)
          Q(5) = V * EXP (0.676 - 15.05/V)                              
C
C       CARBON
C       ------
        ELSEIF (CION.EQ.6) THEN
          Q(1) = 0.061 * V * (1.0 + V / 61.5)
          Q(2) = 9.16 * EXP (-18.68 / V)
          Q(3) = 21.93 - 1.93*V + 0.07*V*V - 0.000753*V*V*V
          Q(4) = 2.66 -1.76*LOG(ABS(V))/V + 17.14/V + 37.87*EXP(-8.51/V)
          Q(5) = 3.24-377.3*LOG(ABS(V))/V/V+276.6/V +0.32*V*EXP(-6.52/V)
          Q(6) = V * EXP (0.676 - 15.05/V)
C
C       OXYGEN
C       ------
        ELSEIF (CION.EQ.8) THEN
          Q(1) = 15.83  - 0.374 *V + 0.003913*V*V - 0.0000138*V*V*V
          Q(2) =  8.84  - 0.29  *V + 0.00283 *V*V
          Q(3) = 27.22  + 1.013 *V - 0.0203  *V*V
          Q(4) = 41.53  - 0.994 *V + 0.0147  *V*V
          Q(5) = 58.75  - 0.685 *V + 0.00927 *V*V
          Q(6) = 25.772 + 0.755 *V - 0.0143  *V*V + 0.0000794*V*V*V
          Q(7) = 31.25  - 0.0553*V
          Q(8) =  1.493 + 1.809 *V - 0.0116  *V*V
C
C       IRON   (NOTE 114:  SIGCX = 1.E-19 * Z, HENCE SET Q AS BELOW)
C       ----
        ELSEIF (CION.EQ.26) THEN
          DO 205 IZ = 1, NIZS
            Q(IZ) = 10.0 * REAL (IZ)
  205     CONTINUE
C
C       OTHER IMPURITIES
C       ----------------
        ELSEIF (IK.EQ.1.AND.IR.EQ.1.AND.CIOPTI.NE.0) THEN
          CALL PRC ('WARNING: NO CHARGE EXCHANGE RECOMBINATION DATA AVAI
     >LABLE FOR THIS IMPURITY')
        ENDIF
C
C       WRITE (6,'('' CXREC: IK,IR,V,NB,Q'',2I5,(1X,6G11.4))')
C    >    IK,IR,V,KNBS(IK,IR),(Q(IZ),IZ=1,NIZS)
C
C-----------------------------------------------------------------------
C     CALCULATE CX RECOMBINATION TIMES, STORE IN KFCXS ARRAY
C-----------------------------------------------------------------------
C
        DO 220 IZ = 1, NIZS
c
c         BV->BIV added by Krieger, IPP/96
c
	  if (cion.eq.5.and.iz.eq.4) then
       	    sigcx = b4atotal(tbgr,ktibs(ik,ir))
c
          else
c
c           Default
c
            SIGCX = Q(IZ) * 1.E-20 * VCX
	  endif
c
          IF (SIGCX.GT.0.0.and.knhs(ik,ir).gt.0.0) THEN
            KFCXS(IK,IR,IZ) = 1.0 / (KNHS(IK,IR) * SIGCX)
          ELSE
            KFCXS(IK,IR,IZ) = 0.0
          ENDIF
c
c
c          Old code - used to combine CX and EI - now done in TAU 
c
c          SIGCX = KNHS(IK,IR) * Q(IZ) * 1.E-20 * VCX / KNBS(IK,IR)
c          IF (KFRCS(IK,IR,IZ).GT.0.0) SIGCX = SIGCX +
c     >      1.0 / (RIZB * KNBS(IK,IR) * KFRCS(IK,IR,IZ))
c          IF (SIGCX.GT.0.0) THEN
c            KFCXS(IK,IR,IZ) = 1.0 / (RIZB * KNBS(IK,IR) * SIGCX)
c          ELSE
c            KFCXS(IK,IR,IZ) = 0.0
c          ENDIF
c
  220   CONTINUE
  230  CONTINUE
C
  240 CONTINUE
c
C-----------------------------------------------------------------------
c
c     CX OPTION 6 - ADAS  
c
C-----------------------------------------------------------------------
c
C     if ciopti = 6,
C     call ADAS for charge exchange recombination coefficient
C
      elseif (ciopti.eq.6) then 

        do 300 ir = 1,nrs
C
          do 310 ik = 1,nks(ir)
            pnesa(ik) = knbs(ik,ir) * rizb
            ptesa(ik) = ktebs(ik,ir)
  310     continue
C
C------ calculate charge exchange recombination rates...
C
cfm - this will disappear when we've got a full set of impurity data,
cfm   including CX, in one place
c          year = '96'
c          usercfm = '/u/cfm/adas'
c          call xxuid(usercfm)
cfm
c         SET base DIVIMP version to use default files 
c
          write(year,'(i2.2)') iyearz
          call xxuid(useridz)
          iclass = 3
C
          do 320 iz = 1,cion
            call adasrd(year,cion,iz,iclass,nks(ir),ptesa,pnesa,
     >                    pcoef(1,iz))
            do 330 ik = 1, nks(ir)
              if (knhs(ik,ir).gt.0.0.and.pcoef(ik,iz).gt.0.0) then
                kfcxs(ik,ir,iz) = 1.0 / (knhs(ik,ir) * pcoef(ik,iz))
                     WRITE(6,*) 'CX RATE:',ik,ir,kfcxs(ik,ir,iz)
              else
                kfcxs(ik,ir,iz) = 0.0
              endif
C
  330       continue
  320     continue
C
  300   continue
c
c       Print out cx cross-section data 
c

        if (cprint.eq.3.or.cprint.eq.9) then 
c
         ntemp = 30

         PTESA(1) = 0.1
         PTESA(2) = 0.25
         PTESA(3) = 0.50
         PTESA(4) = 0.75
         PTESA(5) = 1.0
         PTESA(6) = 1.25
         PTESA(7) = 1.50
         PTESA(8) = 1.75
         PTESA(9) = 2.0
         PTESA(10) = 2.25
         PTESA(11) = 2.5
         PTESA(12) = 2.75
         PTESA(13) = 3.0
         PTESA(14) = 4.0
         PTESA(15) = 5.0
         PTESA(16) = 6.0
         PTESA(17) = 7.0
         PTESA(18) = 8.0
         PTESA(19) = 9.0 
         PTESA(20) = 10.0
         PTESA(21) = 20.0
         PTESA(22) = 30.0
         PTESA(23) = 40.0
         PTESA(24) = 50.0
         PTESA(25) = 60.0
         PTESA(26) = 70.0
         PTESA(27) = 80.0
         PTESA(28) = 90.0
         PTESA(29) = 100.0
         PTESA(30) = 150.0

         do in = 1,4 

            tmpne = 1.0e17 * 10**(in) 

            do ik = 1,ntemp
 
               pnesa(ik) = tmpne

            end do 
C
C---- CALCULATE IONISATION AND RECOMBINATION RATES ...
C
            write(year,'(i2.2)') iyearz
            call xxuid(useridz)            
C
c           CX RECOMBINATION
c
            iclass = 3
c
            DO IZ = 1, CION
               CALL ADASRD(YEAR,CION,IZ,ICLASS,Ntemp,PTESA,PNESA,   
     >              PCOEF(1,IZ))                                      

            end do
c
            write (6,*) 'ADAS CX Recombination rate coefficients:'//
     >                  ' Density = ', tmpne
c
            write(6,1001)  (iz,iz=1,cion)
c
            do ik = 1,ntemp
               write (6,1000) ik,ptesa(ik),
     >                 (pcoef(ik,iz),iz=1,cion)
            end do
c
c        End density loop
c
         end do
c
1001  format(3x,'N',5x,'Tev',2x,12(4x,'IZ =',i2,4x)) 
1000  format(i5,1x,f7.2,12(1x,1p,g13.6))
c
        endif       
c
c
C-----------------------------------------------------------------------
c
c     CX option 7 - data is read from ADPAK tables - This option 
c     is only valied in combination with CDATOPT=2 or 3.
c
C-----------------------------------------------------------------------
c
      elseif (ciopti.eq.7) then 
c
c        ADPAK/STRAHL CX data 
c
         if (cdatopt.eq.2) then 
c
            do ir = 1,nrs
c    
               do ik = 1,nks(ir)
c       
                  do iz = 0,cion
c
c                    Scale the ion temperature by mass in AMU 
c
                     tmpne = knbs(ik,ir) * rizb                  
                     tmpti = ktibs(ik,ir) / crmb
                     tmpte = ktebs(ik,ir)  
c                  
                     call mcrates(tmpte,tmpti,tmpne,
     >                         iz,cion,cion,rion,rrec,rcxr,0)
c
                     IF (rcxr.GT.0.0) THEN
                        KFCXS(IK,IR,IZ)=1.0/(rcxr*knhs(ik,ir))
                     ELSE
                        KFCXS(IK,IR,IZ) = HI
                     ENDIF
c
c
c

               if (check_rates.and.iz.eq.1) then 

                  sigcx = rcxr

                  rate_calc(ik,ir,2) = 
     >                (sigcx * e2datom(ik,ir)) * e2dnzs(ik,ir,1)
     >                * karea2(ik,ir)
                  rate_calc(ik,ir,3) = rate_calc(ik,ir,2)
     >                               + rate_calc(ik,ir,1)

c
c                 Put output here for now
c
c                  if (rate_calc(ik,ir,3).gt.1.1*
c     >                    (e2drec(ik,ir)+e2dcxrec(ik,ir))
c     >                              .or.
c     >                rate_calc(ik,ir,3).lt.0.9*
c     >                    (e2drec(ik,ir)+e2dcxrec(ik,ir)))
c     >                write (6,*) '*** CALC NOTE ***'
                  
                  write (6,'(a,2i4,1p,9(1x,g12.5))') 'Rate_CALC A:',
     >                   ik,ir,
     >                   rate_calc(ik,ir,1), e2drec(ik,ir),
     >                   rate_calc(ik,ir,2), e2dcxrec(ik,ir),
     >                   rate_calc(ik,ir,3), 
     >                              e2drec(ik,ir)+e2dcxrec(ik,ir),
     >                   e2dnes(ik,ir),e2dnzs(ik,ir,1),e2datom(ik,ir) 
c
               endif


c
                  end do

               end do
       
            end do                   
c
c        INEL CX data 
c
         elseif (cdatopt.eq.3) then 
c
c           Loop through background plasma
c
            do ir = 1,nrs
c    
               do ik = 1,nks(ir)
c       
                  do iz = 0,cion
c
c                    Scale the ion temperature by mass in AMU 
c
                     tmpne = knbs(ik,ir) * rizb                  
                     tmpti = ktibs(ik,ir) / crmb
                     tmpte = ktebs(ik,ir)  
c                  
c                    Call for Te-based quantities 
c
                     call imprates(tmpte,iz,cion,rion,rrec,rcxr)
c

c
                     IF (rcxr.GT.0.0) THEN
                        KFCXS(IK,IR,IZ)=1.0/(rcxr*knhs(ik,ir))
                     ELSE
                        KFCXS(IK,IR,IZ) = HI
                     ENDIF


c
c
c

               if (check_rates.and.iz.eq.1) then 

                  sigcx = rcxr

                  rate_calc(ik,ir,2) = 
     >                (sigcx * e2datom(ik,ir)) * e2dnzs(ik,ir,1)
     >                * karea2(ik,ir)
                  rate_calc(ik,ir,3) = rate_calc(ik,ir,2)
     >                               + rate_calc(ik,ir,1)

c
c                 Put output here for now
c
c                  if (rate_calc(ik,ir,3).gt.1.1*
c     >                    (e2drec(ik,ir)+e2dcxrec(ik,ir))
c     >                              .or.
c     >                rate_calc(ik,ir,3).lt.0.9*
c     >                    (e2drec(ik,ir)+e2dcxrec(ik,ir)))
c     >                write (6,*) '*** CALC NOTE ***'
                  
                  write (6,'(a,2i4,1p,9(1x,g12.5))') 'Rate_CALC B:',
     >                   ik,ir,
     >                   rate_calc(ik,ir,1), e2drec(ik,ir),
     >                   rate_calc(ik,ir,2), e2dcxrec(ik,ir),
     >                   rate_calc(ik,ir,3), 
     >                              e2drec(ik,ir)+e2dcxrec(ik,ir),
     >                   e2dnes(ik,ir),e2dnzs(ik,ir,1),e2datom(ik,ir) 
c
               endif




c
                  end do
  
               end do
       
            end do                   
c
         endif
c
C-----------------------------------------------------------------------
c
c     CX OPTION 8 - CF MAGGI PHD THESIS <sigma V> reaction rates as
c                   interpolated/formularized by Tom Rognlien using 
c                   an exponential fit. These data are ONLY valid
c                   for Deuterium/Carbon and the data assumed
c                   Maxwellian velocity distributions for both species
c                   with Td = Tc.
c     CX OPTION 9 - Same as option 8 except that the fitting 
c                   parameters for CII (C+) recombination have been 
c                   modified by Tom Rognlien for a better low
c                   temperature approximation.
c
c
c   * NOTES: Other caveats - this is implemented only for comparison
c            to UEDGE. The CX rates should be using either deuterium
c            neutral temperatures or the carbon ion temperatures - the
c            rates have assumed these are equal and this is clearly 
c            not the case. As a result we use the deuterium ion  
c            temperature in the following calculations. (These values
c            come from a UEDGE background of the appropriate options 
c            have been specified.) 
c
C-----------------------------------------------------------------------
c
C     if ciopti = 8, ciopti = 9
c
C     use formula for charge exchange recombination coefficient
C
      elseif (ciopti.eq.8.or.ciopti.eq.9) then 

        do ir = 1,nrs
C
          do ik = 1,nks(ir)
c
            do iz = 1,cion

               if (ktibs(ik,ir).gt.0.0) then   
c
c                 Valid temperature
c
                  if (ciopti.eq.8) then 
                     sigcx = sigvcx_maggi(iz,ktibs(ik,ir)/crmb)
                  elseif (ciopti.eq.9) then
                     sigcx = sigvcx_maggi_mod(iz,ktibs(ik,ir)/crmb)
                  endif
c
                  if (knhs(ik,ir).gt.0.0.and.sigcx.gt.0.0) then   
c
c                    Check for valid sigmav and neutral density
c
                     kfcxs(ik,ir,iz) = 1.0 / (knhs(ik,ir) * sigcx)
c
                     WRITE(6,*) 'CX RATE 9:',ik,ir,kfcxs(ik,ir,iz)
                  else
                     kfcxs(ik,ir,iz) = 0.0
                  endif
c
               else

                  kfcxs(ik,ir,iz) = 0.0

               endif 
c
c
c              Check rates
c
c               if (check_rates) then 
c
c                     write(6,'(a,3i4,1p,10(1x,g10.4))') 'Rate2:',ik,ir,
c     >                 iz, e2dnzs(ik,ir,iz),e2datom(ik,ir),
c     >                 sigcx,
c     >                 (sigcx*e2datom(ik,ir))*e2dnzs(ik,ir,iz)c
c
c
c               endif
c

               if (check_rates.and.iz.eq.1) then 

                  if (ciopti.eq.8) then  
                     sigcx = sigvcx_maggi(iz,ktibs(ik,ir)/crmb)
                  elseif (ciopti.eq.9) then 
                     sigcx = sigvcx_maggi_mod(iz,ktibs(ik,ir)/crmb)
                  endif

                  rate_calc(ik,ir,2) = 
     >                (sigcx * e2datom(ik,ir)) * e2dnzs(ik,ir,1)
     >                * karea2(ik,ir)
                  rate_calc(ik,ir,3) = rate_calc(ik,ir,2)
     >                               + rate_calc(ik,ir,1)

c
c                 Put output here for now
c
c                  if (rate_calc(ik,ir,3).gt.1.1*
c     >                    (e2drec(ik,ir)+e2dcxrec(ik,ir))
c     >                              .or.
c     >                rate_calc(ik,ir,3).lt.0.9*
c     >                    (e2drec(ik,ir)+e2dcxrec(ik,ir)))
c     >                write (6,*) '*** CALC NOTE ***'
                  
                  write (6,'(a,2i4,1p,9(1x,g12.5))') 'Rate_CALC C:',
     >                   ik,ir,
     >                   rate_calc(ik,ir,1), e2drec(ik,ir),
     >                   rate_calc(ik,ir,2), e2dcxrec(ik,ir),
     >                   rate_calc(ik,ir,3), 
     >                              e2drec(ik,ir)+e2dcxrec(ik,ir),
     >                   e2dnes(ik,ir),e2dnzs(ik,ir,1),e2datom(ik,ir) 
c
               endif


C
            end do
C
          end do

        end do
c
c
c       Print out cx cross-section data 
c

        if (cprint.eq.3.or.cprint.eq.9) then 
c
         ntemp = 30
c
         PTESA(1) = 0.1
         PTESA(2) = 0.25
         PTESA(3) = 0.50
         PTESA(4) = 0.75
         PTESA(5) = 1.0
         PTESA(6) = 1.25
         PTESA(7) = 1.50
         PTESA(8) = 1.75
         PTESA(9) = 2.0
         PTESA(10) = 2.25
         PTESA(11) = 2.5
         PTESA(12) = 2.75
         PTESA(13) = 3.0
         PTESA(14) = 4.0
         PTESA(15) = 5.0
         PTESA(16) = 6.0
         PTESA(17) = 7.0
         PTESA(18) = 8.0
         PTESA(19) = 9.0 
         PTESA(20) = 10.0
         PTESA(21) = 20.0
         PTESA(22) = 30.0
         PTESA(23) = 40.0
         PTESA(24) = 50.0
         PTESA(25) = 60.0
         PTESA(26) = 70.0
         PTESA(27) = 80.0
         PTESA(28) = 90.0
         PTESA(29) = 100.0
         PTESA(30) = 150.0
c
         write (6,*) 'CFM CX Recombination rate coefficient '//
     >                                    'approximations:'
c
         write(6,1001)  (iz,iz=1,cion)
c
         do ik = 1,ntemp
c
            DO IZ = 1, CION
c
               if (ciopti.eq.8) then  
                  pcoef(ik,iz) = sigvcx_maggi(iz,ptesa(ik)/crmb)
               elseif (ciopti.eq.9) then 
                  pcoef(ik,iz) = sigvcx_maggi_mod(iz,ptesa(ik)/crmb)
               endif 
c
            end do
c
            write (6,1000) ik,ptesa(ik),
     >                     (pcoef(ik,iz),iz=1,cion)
c
         end do
c
        endif       
c
C
c     End of CIOPTI IF statement 
c
      endif 
C
      RETURN
      END
c
c
c
      real function b4atotal(t0arg, tzarg)

c     t0 = neutrals temperature (eV)
c     tz = impurity temperature (eV)
c     b4atotal = total charge exchange rate coefficient (m^3/s)

      integer dimt0, dimtz
      parameter (dimt0=7)
      parameter (dimtz=22)

      real t0_t(dimt0), tz_t(dimtz), atotal_t(dimt0,dimtz)
      save t0_t, tz_t, atotal_t

      real t0arg, tzarg, t0, tz, atotalbuf

      data t0_t /1.0e+00, 2.0e+00, 5.0e+00, 1.0e+01, 2.0e+01, 
     >           5.0e+01, 1.0e+02/
      data tz_t /1.0e+00, 2.0e+00, 3.0e+00, 4.0e+00, 5.0e+00, 
     >           6.0e+00, 7.0e+00, 8.0e+00, 9.0e+00, 1.0e+01, 
     >           2.0e+01, 3.0e+01, 4.0e+01, 5.0e+01, 6.0e+01, 
     >           7.0e+01, 8.0e+01, 9.0e+01, 1.0e+02, 2.0e+02, 
     >           5.0e+02, 1.0e+03/

c     atotal_t = charge exchange rate coefficient for process with
c                captured electron in n=3 state
      data atotal_t
     >     /8.85593e-15, 1.11251e-14, 1.45113e-14, 1.82145e-14,
     >      2.38727e-14, 3.65357e-14, 5.27652e-14,
     >      9.38637e-15, 1.14293e-14, 1.46674e-14, 1.83363e-14,
     >      2.39693e-14, 3.65478e-14, 5.29552e-14,
     >      9.85646e-15, 1.17091e-14, 1.48224e-14, 1.84578e-14,
     >      2.40617e-14, 3.66277e-14, 5.30235e-14,
     >      1.02740e-14, 1.19734e-14, 1.49738e-14, 1.85787e-14,
     >      2.41522e-14, 3.66968e-14, 5.30779e-14,
     >      1.06508e-14, 1.22207e-14, 1.51251e-14, 1.86994e-14,
     >      2.42389e-14, 3.67682e-14, 5.31124e-14,
     >      1.09921e-14, 1.24520e-14, 1.52725e-14, 1.88130e-14,
     >      2.43340e-14, 3.68297e-14, 5.31742e-14,
     >      1.13080e-14, 1.26724e-14, 1.54190e-14, 1.89337e-14,
     >      2.44243e-14, 3.68995e-14, 5.32277e-14,
     >      1.15993e-14, 1.28823e-14, 1.55596e-14, 1.90528e-14,
     >      2.45141e-14, 3.69721e-14, 5.32777e-14,
     >      1.18705e-14, 1.30831e-14, 1.57057e-14, 1.91691e-14,
     >      2.46042e-14, 3.70393e-14, 5.33378e-14,
     >      1.21233e-14, 1.32750e-14, 1.58488e-14, 1.92846e-14,
     >      2.46948e-14, 3.71033e-14, 5.33848e-14,
     >      1.40892e-14, 1.49434e-14, 1.72038e-14, 2.04012e-14,
     >      2.55748e-14, 3.77860e-14, 5.39276e-14,
     >      1.56218e-14, 1.63777e-14, 1.84555e-14, 2.14571e-14,
     >      2.64326e-14, 3.84636e-14, 5.44554e-14,
     >      1.69958e-14, 1.76925e-14, 1.96271e-14, 2.24553e-14,
     >      2.72667e-14, 3.91249e-14, 5.49838e-14,
     >      1.82634e-14, 1.89118e-14, 2.07241e-14, 2.34146e-14,
     >      2.80886e-14, 3.97801e-14, 5.55090e-14,
     >      1.94459e-14, 2.00526e-14, 2.17612e-14, 2.43307e-14,
     >      2.88892e-14, 4.04319e-14, 5.60337e-14,
     >      2.05535e-14, 2.11251e-14, 2.27477e-14, 2.52238e-14,
     >      2.96772e-14, 4.10812e-14, 5.65459e-14,
     >      2.15982e-14, 2.21407e-14, 2.36927e-14, 2.60897e-14,
     >      3.04537e-14, 4.17088e-14, 5.70648e-14,
     >      2.25918e-14, 2.31105e-14, 2.46041e-14, 2.69348e-14,
     >      3.12227e-14, 4.23487e-14, 5.75742e-14,
     >      2.35448e-14, 2.40437e-14, 2.54885e-14, 2.77618e-14,
     >      3.19692e-14, 4.29737e-14, 5.80862e-14,
     >      3.17798e-14, 3.21829e-14, 3.33746e-14, 3.53064e-14,
     >      3.89920e-14, 4.89469e-14, 6.29839e-14,
     >      5.09486e-14, 5.12485e-14, 5.21410e-14, 5.36049e-14,
     >      5.64491e-14, 6.44030e-14, 7.61847e-14,
     >      7.47587e-14, 7.49827e-14, 7.56512e-14, 7.67551e-14,
     >      7.89250e-14, 8.51633e-14, 9.48160e-14/

      t0 = min(max(t0arg, t0_t(1)), t0_t(dimt0))
      tz = min(max(tzarg, tz_t(1)), tz_t(dimtz))

      call linin2(t0_t,tz_t,atotal_t,dimt0,dimtz,t0,tz,atotalbuf)
      b4atotal=atotalbuf

      return

      end
c
c
c
      real function sigvcx_maggi_mod(iz,ti)
      implicit none
      integer iz
      real ti
c
c     Model parameter for Tom Rognlien's fit to CF Maggi's PhD thesis
c     CX recombination rate coefficients for D0 + C*
c
      real    M0(6),M1(6),M2(6)
c
c      Revised CX data for lower CX rate for CII
c
      data    M0 / -20.027  , -18.27 ,  -14.48, -14.85, -14.213,-17.576/
      data    M1 /   3.6433 , 2.3657 , 0.05715, 0.5219, 0.42193, 1.8758/
      data    M2 /  -0.59189,-0.29616,0.080947,0.048885,
     >                                             -0.033125, -0.095951/
c
c      data    M0 / -16.104  , -18.27 ,  -14.48, -14.85, -14.213,-17.576/
c      data    M1 /   0.5335 , 2.3657 , 0.05715, 0.5219, 0.42193, 1.8758/
c      data    M2 /-0.0009571,-0.29616,0.080947,0.048885,
c     >                                             -0.033125, -0.095951/
c
      sigvcx_maggi_mod = 10**(M0(iz)+M1(iz)*log10(ti) +
     >                          M2(iz)*log10(ti)**2 )
c
      return
      end   
c
c
c
      real function sigvcx_maggi(iz,ti)
      implicit none
      integer iz
      real ti
c
c     Model parameter for Tom Rognlien's fit to CF Maggi's PhD thesis
c     CX recombination rate coefficients for D0 + C*
c
      real    M0(6),M1(6),M2(6)
c
c      Revised CX data for lower CX rate -   
c
c      data    M0 / -20.027  , -18.27 ,  -14.48, -14.85, -14.213,-17.576/
c      data    M1 /   3.6433 , 2.3657 , 0.05715, 0.5219, 0.42193, 1.8758/
c      data    M2 /  -0.59189,-0.29616,0.080947,0.048885,
c
      data    M0 / -16.104  , -18.27 ,  -14.48, -14.85, -14.213,-17.576/
      data    M1 /   0.5335 , 2.3657 , 0.05715, 0.5219, 0.42193, 1.8758/
      data    M2 /-0.0009571,-0.29616,0.080947,0.048885,
     >                                             -0.033125, -0.095951/
c
      sigvcx_maggi = 10**(M0(iz)+M1(iz)*log10(ti) +
     >                          M2(iz)*log10(ti)**2 )
c
      return
      end   

