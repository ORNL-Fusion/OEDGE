c     -*-Fortran-*-
c
      SUBROUTINE NEUT (NATIZ,MATP,
     >                 MATT,NPROD,NPROD2,FTOT,FYTOT,neut2d_fytot,
     >                 RSTRUK,MTCSTRUK,
     >                 RMAIN,REXIT,
     >                 RATIZ,RNEUT,RWALLN,MTCWALLN,RCENT,RTMAX,SEED,
     >                 NRAND,
     >                 NEUTIM,RFAIL,NYMFS,STATUS)
      implicit none
      DOUBLE PRECISION SEED
      INTEGER   NRAND,NATIZ,MATT,MATP,NPROD,NYMFS,STATUS,NPROD2
      REAL      RFAIL,RSTRUK,NEUTIM,FTOT,FYTOT,RMAIN,REXIT
      real      neut2d_fytot
      REAL      RATIZ,RNEUT,RWALLN,RCENT,RTMAX
      REAL      MTCSTRUK,MTCWALLN
C
C  *********************************************************************
C  *                                                                   *
C  *  NEUT: CONTROL ROUTINE FOR SETTING UP PRIMARY NEUTRALS            *
C  *  -----------------------------------------------------            *
C  *                                                                   *
C  *    THIS ROUTINE CREATES A SET OF PRIMARY NEUTRALS DETAILS TO BE   *
C  *  PASSED TO LAUNCH ROUTINE.  DETAILS CREATED INCLUDE X POSITIONS,  *
C  *  Y POSITIONS, MAXIMUM RANDOM NUMBERS TO BE USED IN                *
C  *  VELOCITY CALCULATIONS AT LAUNCH POINTS, ETC.                     *
C  *                                                                   *
C  * INPUT ARGUMENTS :-                                                *
C  *  NPROD : NUMBER OF NEUTRALS TO LAUNCH                             *
C  *  NPROD2: NUMBER OF SUPPLEMENTARY NEUTRALS TO LAUNCH               *
C  *   SEED : RANDOM GENERATOR SEED  (ALSO PASSED BACK TO DIV)         *
C  *  NRAND : COUNTS TOTAL RANDOMS USED   (PASSED BACK TO DIV TOO)     *
C  *  NYMFS : NUMBER OF POINTS FOR INTERPOLATING YIELD MODIFIER FUNCTN *
C  *                                                                   *
C  * OUTPUTS :-                                                        *
C  * XATIZS : ARRAY WITH X COORDINATES OF INJECTED IONS - IN /CNEUT/   *
C  * YATIZS : ARRAY WITH Y COORDINATES OF INJECTED IONS - IN /CNEUT/   *
C  * PATIZS : ARRAY WITH P COORDINATES OF INJECTED IONS - IN /CNEUT/   *
C  *   VINS : ARRAY WITH VELOCITIES OF INJECTED IONS    - IN /CNEUT/   *
C  * SPUTYS : ARRAY WITH FRAGMENT SIZES OF NEUTRALS  (ALL 1'S) /CNEUT/ *
C  * FSRATE : TIMESTEP USED IN NEUT  - IN /COMTOR/                     *
C  *  NATIZ : NUMBER OF NEUTRALS WHICH IONISED                         *
C  *   MATT : MATERIAL OF TARGET REFERENCE FOR "YIELD" FUNCTION        *
C  *   FTOT : TOT INTEGRATED PRIMARY FLUX                              *
C  *  FYTOT : TOT INTEGRATED PRIMARY FLUX*YIELD                        *
C  *  RATIZ : TOTAL OF IONISED NEUTRAL FRAGMENTS                       *
C  *  RNEUT : TOTAL OF LAUNCHED NEUTRAL FRAGMENTS                      *
C  * RWALLN : TOTAL OF FRAGMENTS PLATING OUT ON WALLS                  *
C  *  RCENT : TOTAL OF FRAGMENTS REACHING CENTRAL MIRROR               *
C  *  RTMAX : TOTAL OF FRAGMENTS EXISTING AT TMAX                      *
C  * RSTRUK : TOTAL OF FRAGMENTS STRIKING TARGET                       *
C  *  RMAIN : TOTAL OF FRAGMENTS REACHING MAIN PLASMA                  *
C  *  REXIT : TOTAL OF FRAGMENTS EXITING MAIN PLASMA                   *
C  * NEUTIM : TIME SPENT TRACKING NEUTRALS SO FAR (CALC. IN LAUNCH)    *
C  *  RFAIL : NUMBER OF FAILED NEUTRAL LAUNCHES (V>VMAX TOO MANY TIMES)*
C  *                                                                   *
C  *                    CHRIS FARRELL (HUNTERSKIL) MARCH 1989          *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE "PARAMS"
      include    'params'
C     INCLUDE "DYNAM1"
      include    'dynam1'
C     INCLUDE "DYNAM3"
      include    'dynam3'
C     INCLUDE "DYNAM4"
      include    'dynam4'
C     INCLUDE "CYIELD"
      include    'cyield'
C     INCLUDE "COMTOR"
      include    'comtor'
C     INCLUDE "CGEOM"
      include    'cgeom'
C     INCLUDE "CIONIZ"
      include    'cioniz'
C     INCLUDE "CNEUT"
      include    'cneut'
      include    'cneut2'
c
c     Include line profile so that initialization may be called if the option is active
c
      include    'line_profile'
c
      include  'printopt'
c
c     Included for optional input variables
c
      include 'slcom' 
c
      EXTERNAL VLAN
C
      REAL      GAMMA,GAMBL,CS,YIELD,FTOT1,FYTOT1,FHTOT1,yldchem
      real      ftot2,fytot2,fhtot2,ftot3,fytot3,fhtot3,fytottemp
      REAL      FHTOT
      REAL      EMAX,RAN,X0,Y0,VLAN
      INTEGER   IPROD,IPOS,ID,IK,IR,IT,IX,IY
      INTEGER   TMPCNEUTB,TMPCNEUTC,tmpcneutc2
      INTEGER   NATIZ1,NATIZ2,natiz3
c
      integer   neuttype
c
      integer   nproda,nprod2a
      integer   nprod_neut2d
c
      real      fydata(maxpts,5),fyprob(maxpts)
      real      totfydata(3,5)
      integer   nfy,nfymap,pinsw,yieldsw,newcneutc,newcneutb
      integer   fymap(maxpts)
c
c      real      ntot 
c      real      wtp,wtc,wwp,wwc,wtot 
c
      integer   natiza, natiz2a 
      real      RATIZa,RNEUTa,RWALLNa,RCENTa,RTMAXa,RSTRUKa 
      real      RFAILa,RMAINa,REXITa  
      real      RATIZ2a,RNEUT2a,RWALLN2a,RCENT2a,RTMAX2a,RSTRUK2a 
      real      RFAIL2a,RMAIN2a,REXIT2a  
      real      ftota,ftot2a,fytota,fytot2a 
      real      ftottempa,ftottemp1,ftottemp2a,ftottemp2
      real      fytottempa,fytottemp1,fytottemp2a,fytottemp2
      real      ftottemp3,fytottemp3
c
      REAL      WLXIMAX
      REAL      RATIZ2,RNEUT2,RWALLN2,RCENT2,RTMAX2
      REAL      RSTRUK2,RFAIL2,RMAIN2,REXIT2
      REAL      RATIZ1,RNEUT1,RWALLN1,RCENT1,RTMAX1
      REAL      RSTRUK1,RFAIL1,RMAIN1,REXIT1
      REAL      RATIZ3,RNEUT3,RWALLN3,RCENT3,RTMAX3
      REAL      RSTRUK3,RFAIL3,RMAIN3,REXIT3
C
C     MOMENTUM TRANSFER COLLISION VARIABLES 
C      
      REAL      MTCSTRUKA,MTCSTRUK1,MTCSTRUK2,MTCSTRUK2A,
     >          MTCSTRUK3
      REAL      MTCWALLNA,MTCWALLN1,MTCWALLN2,MTCWALLN2A,
     >          MTCWALLN3
C
C-----------------------------------------------------------------------
C       SET UP SECTION
C-----------------------------------------------------------------------
C
      write(6,*) ' Entering NEUT ...'
c      write(50,*) ' Entering NEUT ...'
c      write(0,*) ' Entering NEUT ...'
c
c     Initialize
c
      call rzero(eprods,maximp)  
c
      nprod_neut2d = 0
      nproda  = 0
      nprod2a = 0
      natiza  = 0
      natiz2a = 0 
      natiz1  = 0
      natiz2  = 0
      natiz3  = 0 
c
c     Zero counters 
c
c     IPP/08 Krieger - variables were not all initialized
      RATIZ1=0.0
      RNEUT1=0.0
      RWALLN1=0.0
      RCENT1=0.0
      RTMAX1=0.0
      RSTRUK1=0.0
      RFAIL1=0.0
      RMAIN1=0.0
      REXIT1=0.0

      RATIZ2=0.0
      RNEUT2=0.0
      RWALLN2=0.0
      RCENT2=0.0
      RTMAX2=0.0
      RSTRUK2=0.0
      RFAIL2=0.0
      RMAIN2=0.0
      REXIT2=0.0
c     IPP/08 Krieger - end
c
c
      RATIZa=0.0
      RNEUTa=0.0
      RWALLNa=0.0
      RCENTa=0.0
      RTMAXa=0.0
      RSTRUKa=0.0
      RFAILa=0.0
      RMAINa=0.0
      REXITa=0.0
c
      RATIZ2a=0.0
      RNEUT2a=0.0
      RWALLN2a=0.0
      RCENT2a=0.0
      RTMAX2a=0.0
      RSTRUK2a=0.0
      RFAIL2a=0.0
      RMAIN2a=0.0
      REXIT2a=0.0
c
      RATIZ3=0.0
      RNEUT3=0.0
      RWALLN3=0.0
      RCENT3=0.0
      RTMAX3=0.0
      RSTRUK3=0.0
      RFAIL3=0.0
      RMAIN3=0.0
      REXIT3=0.0
c
      MTCSTRUKA=0.0
      MTCSTRUK1=0.0
      MTCSTRUK2=0.0
      MTCSTRUK2A=0.0
      MTCSTRUK3=0.0
c      
      MTCWALLNA=0.0
      MTCWALLN1=0.0
      MTCWALLN2=0.0
      MTCWALLN2A=0.0
      MTCWALLN3=0.0
c
      ftottempa = 0.0 
      ftottemp1 = 0.0 
      ftottemp2a = 0.0 
      ftottemp2 = 0.0 
      ftottemp3  = 0.0
c
      fytottemp = 0.0
      fytottempa = 0.0 
      fytottemp1 = 0.0 
      fytottemp2a = 0.0 
      fytottemp2 = 0.0 
      fytottemp3 = 0.0
c
c     Initialize line profile array
c
      if (line_profile_opt.ne.0) then 
         call init_line_profile_data
      endif

c
c     Set pin data availability switch 
c
      if (uedge_bg.eq.1) then
         pinsw = uedge_bg
      else
         pinsw = cpinopt 
      endif
c
c     Check for 2D NEUT options and prepare data
c
      if (neut2d_opt.eq.1.or.(cneutb.eq.5.and.nprod.gt.0)
     >                   .or.(cneuth.eq.5.and.nprod2.gt.0)) then  
c 
         call prep_neut2d
c
      endif
c
c
c     Copy over value of initial neutral V/A flag ... so that it can
c     be applied to the first group of neutrals launched.
c
c     Note that at the present time the initial neutral velocity/angle
c     flag is applied only to the first group of neutrals launched and
c     not to the group of supplementary neutrals (if any).
c
      TMPCNEUTC = CNEUTC
      CNEUTC = NVAOPT
C
C---- CHECK VALUE FOR FIMP IS OK FOR ALTERNATE SPUTTERING OPTIONS.
C
      NATIZ1 = 0
      NATIZ2 = 0
      CALL PRB
C
C---- CALCULATE GAMMA, FACTOR DETERMINING MAXIMUM ENERGY EXCHANGE.
C
      GAMMA  = 4.0 * CRMB * CRMI / ((CRMB+CRMI) * (CRMB+CRMI))
      GAMBL  = GAMMA * (1.0 - GAMMA)
c
c     CALCULATE values for nproda, nprod, nprod2a, nprod2 based 
c     on the selected Primary and supplementary launch options and 
c     any specified 2D launch option. 
c
      call redistribute_nprod(nproda,nprod,nprod2a,nprod2,
     >                        nprod_neut2d,pinsw,matt,matp)  
c
      if ((nproda+nprod+nprod2a+nprod2).eq.0.0) then 
        call prc ('ERROR: NO IMPURITIES TO FOLLOW!'//
     >            ' - ALL YIELDS ZERO')
        write (6,*) 'ERROR: NO IMPURITIES TO FOLLOW!'//
     >              ' - ALL YIELDS ZERO'
        write (0,*) 'ERROR: NO IMPURITIES TO FOLLOW!'//
     >              ' - ALL YIELDS ZERO'
        stop
      endif
c 
c
C-----------------------------------------------------------------------
c
c     FIRST LAUNCH (PRIMARY)
c
c-----------------------------------------------------------------------
c
c     Launch first groups 
c
      IF (sloutput) WRITE(0,*) 'DEBUG: NEUT -- FIRST LANUCH'
      IF (sloutput) WRITE(0,*)    nproda,nprod,cneutb
      if (nproda.gt.0.or.nprod.gt.0) then  
c
         if (cneutb.eq.0.or.cneutb.eq.3) then 
c
c           Launch first batch - if there is one.
c
            if (nproda.gt.0) then 
c
c              Load physical sputtering data on targets
c
               yieldsw = 0
c  
c              Deal with sputter option 7 - external target flux
c
               if (cneutd.eq.7) then  
                  IF (sloutput) WRITE(0,*) 'DEBUG: CALL TFY A.0'
                  call tfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >               3,yieldsw,matp,matt)              
               else
                  IF (sloutput) WRITE(0,*) 'DEBUG: CALL TFY A'
                  call tfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >               pinsw,yieldsw,matp,matt)               
               endif
c
               ftottempa = totfydata(3,1) 
               fytottempa = totfydata(3,5) 
c
c              Print
c 
               call printfy(fydata,fymap,fyprob,nfy,nfymap,totfydata)
c
               newcneutb = 0
               newcneutc = cneutc
c


               call neutbatch(newcneutc,newcneutb,yieldsw,pinsw,
     >                     0,nproda,0,natiza,gambl,
     >                     RSTRUKa,MTCSTRUKA,RMAINa,REXITa,RATIZa,
     >                     RNEUTa,RWALLNa,MTCWALLNA,RCENTa,RTMAXa,
     >                     rfaila,
     >                     SEED,NRAND,NEUTIM,STATUS,MATP,MATT,
     >                     fymap,fyprob,nfymap,fydata)

c
            endif
c
c           Launch second batch
c
            if (nprod.gt.0) then 
c
c              Physically sputtered particles
c
               if (cneutd.eq.0.or.cneutd.eq.1.or.
     >             cneutd.eq.3.or.cneutd.eq.4.or.
     >             cneutd.eq.7.or.cneutd.eq.8) then  
c
                  yieldsw = 0 
c
                  if (cneutd.eq.7) then  
                     IF (sloutput) WRITE(0,*) 'DEBUG: CALL TFY B.0'
                     call tfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  3,yieldsw,matp,matt)              
                  else
                     IF (sloutput) WRITE(0,*) 'DEBUG: CALL TFY B'
                     call tfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  pinsw,yieldsw,matp,matt)               
                  endif
c
                  newcneutc = cneutc
                  newcneutb = 0 
c
c              Chemically sputtered particles
c
               elseif (cneutd.eq.2.or.cneutd.eq.5.or.cneutd.eq.6) then 
c
                  yieldsw = 1
c
                  IF (sloutput) WRITE(0,*) 'DEBUG: CALL TFY C'
                  call tfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  pinsw,yieldsw,matp,matt)               
c
                  if (cneutd.eq.6) then  
                     newcneutc = 3
                     newcneutb = 0
                  else  
                     newcneutc = cneutc
                     newcneutb = 0
                  endif  
c
               endif
c
               ftottemp1 = totfydata(3,1) 
               fytottemp1 = totfydata(3,5) 
c
               call printfy(fydata,fymap,fyprob,nfy,nfymap,totfydata)
c

               call neutbatch(newcneutc,newcneutb,yieldsw,pinsw,
     >                     nproda,nprod,natiza,natiz1,gambl,
     >                     RSTRUK1,MTCSTRUK1,RMAIN1,REXIT1,RATIZ1,
     >                     RNEUT1,RWALLN1,MTCWALLN1,RCENT1,RTMAX1,
     >                     rfail1,
     >                     SEED,NRAND,NEUTIM,STATUS,MATP,MATT,
     >                     fymap,fyprob,nfymap,fydata)
c
             
            endif
c
c           Update variables
c
            natiz1 = natiz1 + natiza
c
c           NOTE: fluxes are not cumulative - it is the same flux
c           causing the sputtering of both nprod and nproda - however
c           there are times when nproda is non-zero while nprod is 
c           zero - this resulted in ftottemp not being assigned 
c           correctly.
c    
c
            if (ftottemp1.eq.0.0.and.ftottempa.ne.0.0) then 
               ftottemp1 = ftottempa 
            endif
c
c           ftottemp1 = ftottemp1 + ftottempa
c
            fytottemp1 = fytottemp1 + fytottempa
c
            RATIZ1  = RATIZ1  + RATIZA 
            RNEUT1  = RNEUT1  + RNEUTA 
            RWALLN1 = RWALLN1 + RWALLNA
            MTCWALLN1 = MTCWALLN1 + MTCWALLNA
            RCENT1  = RCENT1  + RCENTA 
            RTMAX1  = RTMAX1  + RTMAXA 
            RSTRUK1 = RSTRUK1 + RSTRUKA
            MTCSTRUK1 = MTCSTRUK1 + MTCSTRUKA
            RFAIL1  = RFAIL1  + RFAILA
            RMAIN1  = RMAIN1  + RMAINA 
            REXIT1  = REXIT1  + REXITA 


c
c        Free space point launch -  
c
         elseif (cneutb.eq.1.or.cneutb.eq.6.or.cneutb.eq.7) then

            call prb
            call prc('FREE-SPACE LAUNCH: FLUX AND YIELD DATA ARE')
            call prc('NOT PRINTED OR USED FOR THE LAUNCH.')
            call prb
c
            newcneutc = cneutc
            newcneutb = cneutb
c
            ftottemp1  = 0.0
            fytottemp1 = 0.0
c

            call neutbatch(newcneutc,newcneutb,yieldsw,pinsw,
     >                     0,nprod,0,natiz1,gambl,
     >                     RSTRUK1,MTCSTRUK1,RMAIN1,REXIT1,RATIZ1,
     >                     RNEUT1,RWALLN1,MTCWALLN1,RCENT1,RTMAX1,
     >                     rfail1,
     >                     SEED,NRAND,NEUTIM,STATUS,MATP,MATT,
     >                     fymap,fyprob,nfymap,fydata)
c
c         Wall launches
c
          elseif (cneutb.eq.2.or.cneutb.eq.4) then  
c
c            Launch first batch - if there is one.
c
             if (nproda.gt.0) then 
c
c               Load physical sputtering data on walls
c
                yieldsw = 0
c
                if (cneutd.eq.8) then 
                   call wfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  2,yieldsw,matp,matt)               
                else
                   IF (sloutput) WRITE(0,*) 'DEBUG: CALL WFY A'
                   call wfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  pinsw,yieldsw,matp,matt)               
                endif
c
                ftottempa = totfydata(3,1) 
                fytottempa = totfydata(3,5) 
c
c               Print
c
                call printfy(fydata,fymap,fyprob,nfy,nfymap,totfydata)
c
                newcneutb = 2
                newcneutc = cneutc
c

                call neutbatch(newcneutc,newcneutb,yieldsw,pinsw,
     >                     0,nproda,0,natiza,gambl,
     >                     RSTRUKa,MTCSTRUKA,RMAINa,REXITa,RATIZa,
     >                     RNEUTa,RWALLNa,MTCWALLNA,RCENTa,RTMAXa,
     >                     rfaila,
     >                     SEED,NRAND,NEUTIM,STATUS,MATP,MATT,
     >                     fymap,fyprob,nfymap,fydata)
c
             endif
c
c            Launch second batch
c
             if (nprod.gt.0) then 
c
c               Physical Sputtering 
c
                if (cneutd.eq.0.or.cneutd.eq.1.or.
     >              cneutd.eq.3.or.cneutd.eq.4.or.
     >              cneutd.eq.7.or.cneutd.eq.8) then  
c
                   yieldsw = 0 
c
                   if (cneutd.eq.8) then 
                      call wfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  2,yieldsw,matp,matt)               
                   else
                      IF (sloutput) WRITE(0,*) 'DEBUG: CALL WFY B'
                      call wfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  pinsw,yieldsw,matp,matt)               
                   endif 

c
                   newcneutc = cneutc
                   newcneutb = 2 
c
c               Chemical Sputtering 
c
                elseif (cneutd.eq.2.or.cneutd.eq.5.or.cneutd.eq.6) then 
c
                   yieldsw = 1
c
                   IF (sloutput) WRITE(0,*) 'DEBUG: CALL WFY C'
                   call wfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  pinsw,yieldsw,matp,matt)               
c
                   if (cneutd.eq.6) then  
                      newcneutc = 3
                      newcneutb = 2
                   else  
                      newcneutc = cneutc
                      newcneutb = 2
                   endif  
c
                endif
c
                ftottemp1 = totfydata(3,1) 
                fytottemp1 = totfydata(3,5) 
c
                call printfy(fydata,fymap,fyprob,nfy,nfymap,totfydata)
c

                call neutbatch(newcneutc,newcneutb,yieldsw,pinsw,
     >                     nproda,nprod,natiza,natiz1,gambl,
     >                     RSTRUK1,MTCSTRUK1,RMAIN1,REXIT1,RATIZ1,
     >                     RNEUT1,RWALLN1,MTCWALLN1,RCENT1,RTMAX1,
     >                     rfail1,
     >                     SEED,NRAND,NEUTIM,STATUS,MATP,MATT,
     >                     fymap,fyprob,nfymap,fydata)
     
c
             endif
c
c            Update variables
c
             natiz1 = natiz1 + natiza
c
c           NOTE: fluxes are not cumulative - it is the same flux
c           causing the sputtering of both nprod and nproda - however
c           there are times when nproda is non-zero while nprod is 
c           zero - this resulted in ftottemp not being assigned 
c           correctly.
c
            if (ftottemp1.eq.0.0.and.ftottempa.ne.0.0) then 
               ftottemp1 = ftottempa 
            endif
c
c            ftottemp1 = ftottemp1 + ftottempa
c
             fytottemp1 = fytottemp1 + fytottempa
c
             RATIZ1  = RATIZ1  + RATIZA 
             RNEUT1  = RNEUT1  + RNEUTA  
             RWALLN1 = RWALLN1 + RWALLNA
             MTCWALLN1 = MTCWALLN1 + MTCWALLNA
             RCENT1  = RCENT1  + RCENTA 
             RTMAX1  = RTMAX1  + RTMAXA 
             RSTRUK1 = RSTRUK1 + RSTRUKA
             MTCSTRUK1 = MTCSTRUK1 + MTCSTRUKA
             RFAIL1  = RFAIL1  + RFAILA
             RMAIN1  = RMAIN1  + RMAINA 
             REXIT1  = REXIT1  + REXITA 
c
c         2D free-space neutral launch 
c
          elseif (cneutb.eq.5) then 
c
            call prb
            call prc('2D FREE-SPACE LAUNCH: FLUX AND YIELD DATA ARE')
            call prc('NOT PRINTED OR USED FOR THE LAUNCH.')
            call prb
c
            newcneutc = cneutc
            newcneutb = cneutb
c
            ftottemp1  = 0.0
            fytottemp1 = 0.0
c

            call neutbatch(newcneutc,newcneutb,yieldsw,pinsw,
     >                     0,nprod,0,natiz1,gambl,
     >                     RSTRUK1,MTCSTRUK1,RMAIN1,REXIT1,RATIZ1,
     >                     RNEUT1,RWALLN1,MTCWALLN1,RCENT1,RTMAX1,
     >                     rfail1,
     >                     SEED,NRAND,NEUTIM,STATUS,MATP,MATT,
     >                     fymap,fyprob,nfymap,fydata)
c
            neut2d_fytot = neut2d_src
c
c
c         ENDIF - Corresponding to cneutb test above.
c
          endif
c
      endif 
C
C     RESTORE CNEUTC to its original value from the initial V/A flag.
C
      cneutc = tmpcneutc
c
c
c-----------------------------------------------------------------------
c
c     SECOND LAUNCH (SUPPLEMENTAL)
c
c-----------------------------------------------------------------------
c
c     Launch second batches of neutrals if there are any
c
      IF (sloutput) WRITE(0,*) 'DEBUG: NEUT -- SECOND LANUCH'
      IF (sloutput) WRITE(0,*)    nprod2a,nprod2
      if (nprod2a.gt.0.or.nprod2.gt.0) then 
c
         if (cneuth.eq.0.or.cneuth.eq.3) then 
c
c           Target launches  
c
c           Launch first batch - if there is one.
c
            if (nprod2a.gt.0) then 
c
               yieldsw = 0
c
c              Load physical sputtering data on targets
c
               if (cneutd2.eq.7) then  
                  IF (sloutput) WRITE(0,*) 'DEBUG: CALL TFY D.0'
                  call tfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >               3,yieldsw,matp,matt)              
               else
                  IF (sloutput) WRITE(0,*) 'DEBUG: CALL TFY D'
                  call tfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >               pinsw,yieldsw,matp,matt)               
               endif
c
               ftottemp2a = totfydata(3,1) 
               fytottemp2a = totfydata(3,5) 
c
c              Print
c
               call printfy(fydata,fymap,fyprob,nfy,nfymap,totfydata)
c
               newcneutb = 0
               newcneutc = cneuti
c

               call neutbatch(newcneutc,newcneutb,yieldsw,pinsw,
     >                     nproda+nprod,nprod2a,natiz1,natiz2a,gambl,
     >                     RSTRUK2a,MTCSTRUK2A,RMAIN2a,REXIT2a,RATIZ2a,
     >                     RNEUT2a,RWALLN2a,MTCWALLN2A,RCENT2a,RTMAX2a,
     >                     rfail2a,
     >                     SEED,NRAND,NEUTIM,STATUS,MATP,MATT,
     >                     fymap,fyprob,nfymap,fydata)
c
            endif
c
c           Launch second batch
c
            if (nprod2.gt.0) then 
c
               if (cneutd2.eq.0.or.cneutd2.eq.1.or.
     >             cneutd2.eq.3.or.cneutd2.eq.4.or.
     >             cneutd2.eq.7.or.cneutd2.eq.8) then  
c
                  yieldsw = 0
c
                  if (cneutd2.eq.7) then  
                     IF (sloutput) WRITE(0,*) 'DEBUG: CALL TFY E.0'
                     call tfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  3,yieldsw,matp,matt)              
                  else
                     IF (sloutput) WRITE(0,*) 'DEBUG: CALL TFY E'
                     call tfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  pinsw,yieldsw,matp,matt)               
                  endif
c
                  newcneutc = cneuti
                  newcneutb = 0 
c
               elseif (cneutd2.eq.2.or.cneutd2.eq.5.or.
     >                 cneutd2.eq.6) then 
c
                  yieldsw = 1
c
                  IF (sloutput) WRITE(0,*) 'DEBUG: CALL TFY F'
                  call tfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  pinsw,yieldsw,matp,matt)               
c
                  if (cneutd2.eq.6) then  
                     newcneutc = 3
                     newcneutb = 0
                  else  
                     newcneutc = cneuti
                     newcneutb = 0
                  endif  
c
               endif
c
               ftottemp2 = totfydata(3,1) 
               fytottemp2 = totfydata(3,5) 
c
               call printfy(fydata,fymap,fyprob,nfy,nfymap,totfydata)
c

               call neutbatch(newcneutc,newcneutb,yieldsw,pinsw,
     >                     nprod+nproda+nprod2a,nprod2,natiz1+natiz2a,
     >                     natiz2,
     >                     gambl,
     >                     RSTRUK2,MTCSTRUK2,RMAIN2,REXIT2,RATIZ2,
     >                     RNEUT2,RWALLN2,MTCWALLN2,RCENT2,RTMAX2,
     >                     rfail2,
     >                     SEED,NRAND,NEUTIM,STATUS,MATP,MATT,
     >                     fymap,fyprob,nfymap,fydata)
c
            endif
c
c           Update variables
c
            natiz2 = natiz2 + natiz2a
c
c           NOTE: fluxes are not cumulative - it is the same flux
c           causing the sputtering of both nprod2 and nprod2a - however
c           there are times when nproda is non-zero while nprod is 
c           zero - this resulted in ftottemp not being assigned 
c           correctly.
c
            if (ftottemp2.eq.0.0.and.ftottemp2a.ne.0.0) then 
               ftottemp2 = ftottemp2a 
            endif
c
c           ftottemp2 = ftottemp2 + ftottemp2a
c
            fytottemp2 = fytottemp2 + fytottemp2a
c
            RATIZ2  = RATIZ2  + RATIZ2A 
            RNEUT2  = RNEUT2  + RNEUT2A 
            RWALLN2 = RWALLN2 + RWALLN2A
            MTCWALLN2 = MTCWALLN2 + MTCWALLN2A
            RCENT2  = RCENT2  + RCENT2A 
            RTMAX2  = RTMAX2  + RTMAX2A 
            RSTRUK2 = RSTRUK2 + RSTRUK2A
            MTCSTRUK2 = MTCSTRUK2 + MTCSTRUK2A
            RFAIL2  = RFAIL2  + RFAIL2A
            RMAIN2  = RMAIN2  + RMAIN2A 
            REXIT2  = REXIT2  + REXIT2A 
c
         elseif (cneuth.eq.1) then 
c
c           Free space point launch  
c
            call prb
            call prc('FREE-SPACE LAUNCH: FLUX AND YIELD DATA ARE')
            call prc('NOT PRINTED OR USED FOR THE LAUNCH.')
            call prb
c
            newcneutc = cneuti
            newcneutb = cneuth
c
            ftottemp2  = 0.0
            fytottemp2 = 0.0
c

            call neutbatch(newcneutc,newcneutb,yieldsw,pinsw,
     >                     nprod+nproda,nprod2,natiz1,
     >                     natiz2,
     >                     gambl,
     >                     RSTRUK2,MTCSTRUK2,RMAIN2,REXIT2,RATIZ2,
     >                     RNEUT2,RWALLN2,MTCWALLN2,RCENT2,RTMAX2,
     >                     rfail2,
     >                     SEED,NRAND,NEUTIM,STATUS,MATP,MATT,
     >                     fymap,fyprob,nfymap,fydata)
c
c        WALL LAUNCHES
c
         elseif (cneuth.eq.2.or.cneuth.eq.4) then    
c
c           Launch first batch - if there is one.
c
            if (nprod2a.gt.0) then 
c
               yieldsw = 0
c
c              Load physical sputtering data on walls
c
               if (cneutd2.eq.8) then
                  call wfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >               2,yieldsw,matp,matt)               
               else 
                  IF (sloutput) WRITE(0,*) 'DEBUG: CALL WFY D'
                  call wfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >               pinsw,yieldsw,matp,matt)               
               endif
c  
               ftottemp2a = totfydata(3,1) 
               fytottemp2a = totfydata(3,5) 
c
c              Print
c
               call printfy(fydata,fymap,fyprob,nfy,nfymap,totfydata)
c
               newcneutb = 2
               newcneutc = cneuti
c
               call neutbatch(newcneutc,newcneutb,yieldsw,pinsw,
     >                     nproda+nprod,nprod2a,natiz1,natiz2a,gambl,
     >                     RSTRUK2a,MTCSTRUK2A,RMAIN2a,REXIT2a,RATIZ2a,
     >                     RNEUT2a,RWALLN2a,MTCWALLN2A,RCENT2a,RTMAX2a,
     >                     rfail2a,
     >                     SEED,NRAND,NEUTIM,STATUS,MATP,MATT,
     >                     fymap,fyprob,nfymap,fydata)
c
            endif
c
c           Launch second batch
c
            if (nprod2.gt.0) then 
c
               if (cneutd2.eq.0.or.cneutd2.eq.1.or.
     >             cneutd2.eq.3.or.cneutd2.eq.4.or.
     >             cneutd2.eq.7.or.cneutd2.eq.8) then  
c
                  yieldsw = 0
c
                  if (cneutd2.eq.8) then 
                     call wfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  2,yieldsw,matp,matt)               
                  else
                     IF (sloutput) WRITE(0,*) 'DEBUG: CALL WFY E'
                     call wfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  pinsw,yieldsw,matp,matt)               
                  endif
c
                  newcneutc = cneuti
                  newcneutb = 2 
c
               elseif (cneutd2.eq.2.or.cneutd2.eq.5.or.
     >                 cneutd2.eq.6) then 
c
                  yieldsw = 1 
c
                  IF (sloutput) WRITE(0,*) 'DEBUG: CALL WFY F'
                  call wfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  pinsw,yieldsw,matp,matt)               
c
                  if (cneutd2.eq.6) then  
                     newcneutc = 3
                     newcneutb = 2
                  else  
                     newcneutc = cneuti
                     newcneutb = 2
                  endif  
c
               endif
c
               ftottemp2 = totfydata(3,1) 
               fytottemp2 = totfydata(3,5) 
c
               call printfy(fydata,fymap,fyprob,nfy,nfymap,totfydata)
c
               call neutbatch(newcneutc,newcneutb,yieldsw,pinsw,
     >                     nproda+nprod+nprod2a,nprod2,natiz1+natiz2a,
     >                     natiz2,
     >                     gambl,
     >                     RSTRUK2,MTCSTRUK2,RMAIN2,REXIT2,RATIZ2,
     >                     RNEUT2,RWALLN2,MTCWALLN2,RCENT2,RTMAX2,
     >                     rfail2,
     >                     SEED,NRAND,NEUTIM,STATUS,MATP,MATT,
     >                     fymap,fyprob,nfymap,fydata)
c
            endif
c
c           Update variables
c
            natiz2 = natiz2 + natiz2a
c
c           NOTE: fluxes are not cumulative - it is the same flux
c           causing the sputtering of both nprod2 and nprod2a - however
c           there are times when nproda is non-zero while nprod is 
c           zero - this resulted in ftottemp not being assigned 
c           correctly.
c
            if (ftottemp2.eq.0.0.and.ftottemp2a.ne.0.0) then 
               ftottemp2 = ftottemp2a 
            endif
c
c            ftottemp2 = ftottemp2 + ftottemp2a
c
            fytottemp2 = fytottemp2 + fytottemp2a
c
            RATIZ2  = RATIZ2  + RATIZ2A 
            RNEUT2  = RNEUT2  + RNEUT2A 
            RWALLN2 = RWALLN2 + RWALLN2A
            MTCWALLN2 = MTCWALLN2 + MTCWALLN2A
            RCENT2  = RCENT2  + RCENT2A 
            RTMAX2  = RTMAX2  + RTMAX2A 
            RSTRUK2 = RSTRUK2 + RSTRUK2A
            MTCSTRUK2 = MTCSTRUK2 + MTCSTRUK2A
            RFAIL2  = RFAIL2  + RFAIL2A
            RMAIN2  = RMAIN2  + RMAIN2A 
            REXIT2  = REXIT2  + REXIT2A 
c
        elseif (cneuth.eq.5) then
c
c           2D Free space point launch  
c
            call prb
            call prc('2D NEUTRAL FREE-SPACE LAUNCH: FLUX'//
     >               ' AND YIELD DATA ARE')
            call prc('NOT PRINTED OR USED FOR THE LAUNCH.')
            call prb
c
            newcneutc = cneuti
            newcneutb = cneuth
c
            ftottemp2  = 0.0
            fytottemp2 = 0.0
c
            call neutbatch(newcneutc,newcneutb,yieldsw,pinsw,
     >                     nprod+nproda,nprod2,natiz1,
     >                     natiz2,
     >                     gambl,
     >                     RSTRUK2,MTCSTRUK2,RMAIN2,REXIT2,RATIZ2,
     >                     RNEUT2,RWALLN2,MTCWALLN2,RCENT2,RTMAX2,
     >                     rfail2,
     >                     SEED,NRAND,NEUTIM,STATUS,MATP,MATT,
     >                     fymap,fyprob,nfymap,fydata)
c
            neut2d_fytot = neut2d_src
        
c
c        ENDIF - Corresponding to cneuth test above.
c
         endif  
c
      endif 
c
c     NOW - launch any THIRD launch 2D neutrals that are called for ...
c
      if (neut2d_opt.eq.1) then 
c
c        2D Free space point launch  
c	 
         call prb
         call prc('-----------------------------------------------')
         call prb
         call prc('2D NEUTRAL FREE-SPACE LAUNCH: FLUX'//
     >            ' AND YIELD DATA ARE')
         call prc('NOT PRINTED OR USED FOR THE LAUNCH.')
         call prb
         call prc('-----------------------------------------------')
         call prb
c	 
         newcneutc = neut2d_vaopt
         newcneutb = 5
c	 
         ftottemp3  = 0.0
         fytottemp3 = 0.0
c
         call print_neut2d
c	 
         call neutbatch(newcneutc,newcneutb,yieldsw,pinsw,
     >                  nprod+nproda+nprod2a+nprod2,nprod_neut2d,
     >                  natiz1+natiz2,
     >                  natiz3,
     >                  gambl,
     >                  RSTRUK3,MTCSTRUK3,RMAIN3,REXIT3,RATIZ3,
     >                  RNEUT3,RWALLN3,MTCWALLN3,RCENT3,RTMAX3,
     >                  rfail3,
     >                  SEED,NRAND,NEUTIM,STATUS,MATP,MATT,
     >                  fymap,fyprob,nfymap,fydata)
c	 
         neut2d_fytot = neut2d_src
c
      endif 
c
c
c-----------------------------------------------------------------------
c
c     Common processing
c
c-----------------------------------------------------------------------
c      
c     Totals
c
c
c     Copy over fytottemp into fytot for passing back to DIV module
c
      fytot = fytottemp1 + fytottemp2 + fytottemp3
      ftot  = ftottemp1 + ftottemp2 + ftottemp3
c
      NATIZ = NATIZ1 + NATIZ2 + natiz3
      RATIZ = RATIZ1 + RATIZ2 + ratiz3
      RNEUT = RNEUT1 + RNEUT2 + rneut3
      RWALLN= RWALLN1+ RWALLN2 + rwalln3
      MTCWALLN= MTCWALLN1+ MTCWALLN2+ mtcwalln3
      RCENT = RCENT1 + RCENT2 + rcent3
      RTMAX = RTMAX1 + RTMAX2 + rtmax3
      RSTRUK= RSTRUK1+ RSTRUK2 +rstruk3
      MTCSTRUK= MTCSTRUK1+ MTCSTRUK2+ MTCSTRUK3 
      RFAIL = RFAIL1 + RFAIL2 + rfail3
      RMAIN = RMAIN1 + RMAIN2 + rmain3
      REXIT = REXIT1 + REXIT2 + rexit3
c
      write(6,'(a,7g12.5,2i6)') 'NEUT DATA:',rneut,ratiz,rstruk,rexit,
     >                rcent,rfail,rwalln,nprod,nprod2


c
c-----------------------------------------------------------------------
c
c     Load Physical sputtering on target Flux and Yield 
c     data into arrays for passing to the OUT program. 
c
      if (cneutd.eq.0.or.cneutd.eq.1.or.cneutd.eq.3.or.
     >    cneutd.eq.4.or.cneutd.eq.6.or.cneutd.eq.7.or.
     >    cneutd.eq.8) then 
c
         if (cneutd.eq.7) then  
            IF (sloutput) WRITE(0,*) 'DEBUG: CALL TFY G.0'
            call tfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  3,0,matp,matt)               
         else
            IF (sloutput) WRITE(0,*) 'DEBUG: CALL TFY G'
            call tfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  pinsw,0,matp,matt)               
         endif
c
      elseif (cneutd.eq.2.or.cneutd.eq.5) then 
c
         IF (sloutput) WRITE(0,*) 'DEBUG: CALL TFY H'
         call tfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  pinsw,1,matp,matt)               
      endif
c
      do id = 1,nds
         kflux(id) = fydata(id,1)
         kener(id) = fydata(id,2)
         kheat(id) = fydata(id,3)
         kyield(id)= fydata(id,4)
         kfy(id)   = fydata(id,5)
      end do                            
c
C
C-----------------------------------------------------------------------
C       STORE PRIMARY DETAILS IN "IONISATION -1" POSITION IN ARRAYS
C-----------------------------------------------------------------------
C
C---- WHEN PLOTTING GRAPHS, THE "-1" POSITION WILL BE CALLED THE
C---- "PRIMARY NEUTRAL" LINES, AND THE "0" POSITION WILL BE CALLED THE
C---- "TOTAL NEUTRAL" LINES.  HENCE NOW THAT ALL PRIMARIES HAVE BEEN
C---- LAUNCHED, BOTH POSITIONS SHOULD BE IDENTICAL.  FOR SUBSEQUENT
C---- LAUNCHES OF SECONDARY, TERTIARY, ..ETC NEUTRALS ONLY THE "0"
C---- POSITION WILL BE AFFECTED IN LAUNCH.
C
      DO 360 IR = 1, NRS
        DO 350 IK = 1, NKS(IR)
          DDLIMS(IK,IR,-1) = DDLIMS(IK,IR,0)
          DDTS(IK,IR,-1) = DDTS(IK,IR,0)
          TIZS  (IK,IR,-1) = TIZS  (IK,IR,0)
          DO 340 IT = 1, NTS
            LIMS(IK,IR,-1,IT) = LIMS(IK,IR,0,IT)
  340     CONTINUE
  350   CONTINUE
  360 CONTINUE
      DO 380 IK = 1, MAXNKS
        ELIMS(IK,1,-1) = ELIMS(IK,1,0)
        ELIMS(IK,2,-1) = ELIMS(IK,2,0)
        ELIMS(IK,3,-1) = ELIMS(IK,3,0)
  380 CONTINUE
      DO 390 ID = 1, NDS
        NEROS(ID,2) = NEROS(ID,3)
  390 CONTINUE


C
      RETURN
C
c     changed format, Krieger IPP, 6/95
c 9000 FORMAT(4X,'   R      Z   YMF   FLUX   ENERGY   YIELD    F*Y  ',
c     >  ' MAXRAN MAXVIN CUMPROB')                                       
c     changed format, Krieger IPP, 6/95
c 9002 FORMAT(2X,2(1x,F6.3),F5.2,1P,G9.2,0P,F7.2,1P,2G9.2,0P,F5.2,
c     >  1P,G8.1,0P,F6.3)
c
c
 9000 FORMAT(4X,'   R      Z   YMF   FLUX   ENERGY   YIELD    F*Y  ',
     >  '  Bt/Bth  LENGTH   NONORTH')
 9002 FORMAT(2X,2(1x,f6.3),F5.2,1P,G9.2,0P,F7.2,1P,2G9.2,0P,1x,
     >       F5.2,1P,
     >       1x,G9.2,0P,1x,G9.2)
 9008 FORMAT(5X,I4,5X,G12.5,4X,G12.5,6X,G14.8)
 9009 FORMAT('   INDEX  POSITION OF SEGMENT START (R,Z)',
     >       '          PROBABILITY')
      END
C
C
C
      SUBROUTINE LAUNCH (LPROD,NPROD,LATIZ,NATIZ,SSTRUK,MTCSTRUK,
     >                   SMAIN,SEXIT,
     >                   SATIZ,SNEUT,SWALLN,MTCWALLN,SCENT,STMAX,
     >                   SEED,NRAND,
     >                   NEUTIM,SFAIL,STATUS,MATP,MATT,neuttype)
      use mtc
      use velocity_dist
c slmod begin
      use mod_interface
      use mod_divimp
c slmod end
      IMPLICIT NONE
c
      INTEGER    LPROD,NPROD,LATIZ,NATIZ,NRAND,STATUS,MATP,MATT
      integer    neuttype
      REAL       SATIZ,SNEUT,SWALLN,SCENT,STMAX,SSTRUK,NEUTIM,SFAIL
      REAL       SMAIN,SEXIT
c slmod begin - temp
      INTEGER incell,target_loss,iteration_limit
      LOGICAL :: iteration_warning = .FALSE.
c slmod end
C
      REAL       MTCSTRUK,MTCWALLN
C
      DOUBLE PRECISION SEED
C
C  *********************************************************************
C  *                                                                   *
C  *  LAUNCH:  THIS ROUTINE CREATED TO IMPLEMENT NOTE 87 ON SELF -     *
C  *  SPUTTERING.  FOLLOWS NEUTRALS FROM LAUNCH POSITIONS INPUT        *
C  *  VIA THE ARGUMENT LIST UNTIL EVENTUAL IONISATION.  THE IONISATION *
C  *  POSITIONS ARE RETURNED, TOGETHER WITH LAUNCH VELOCITIES, ETC.    *
C  *    ALL NEUTRALS ARE FOLLOWED USING THE SAME TIMESTEP, FSRATE,     *
C  *  NO MATTER WHAT THE POSITION OF THE NEUTRAL HAPPENS TO BE.        *
C  *                                                                   *
C  *                                                                   *
C  *  INPUT:  XPRODS ) LAUNCH COORDINATES FROM /CNEUT/                 *
C  *          YPRODS )                                                 *
C  *          PPRODS )                                                 *
C  *          RMAXS    MAXIMUM RANDOM NUMBERS ALLOWED IN VIN CALCS.    *
C  *                   USED IN SELF-SPUTTER CASES TO HOLD EIMPS /CNEUT/*
C  *          SPUTYS   FRACTIONS OF ATOMS TO BE LAUNCHED - IN /CNEUT/  *
C  *          FSRATE   TIMESTEP TO BE USED IN FOLLOWING NEUTRALS       *
C  *          LPROD    FIRST NUMBER OF NEUTRAL FRAGMENTS TO START WITH *
C  *          NPROD    CONTINUE TO THIS POINT. TOTAL = NPROD-LPROD+1   *
C  *          SEED     RANDOM NUMBER GENERATOR SEED  (PASSED BACK TOO) *
C  *          NRAND    COUNTS TOTAL RANDOMS USED     (PASSED BACK TOO) *
C  *          NEUTIM   TIME USED TRACKING NEUTRALS   (PASSED BACK TOO) *
C  *          MATP     BOMBARDING ION TYPE, USED FOR YIELDS            *
C  *          MATT     TARGET  TYPE, USED FOR YIELDS                   *
C  *                                                                   *
C  *  OUTPUT: XATIZS ) FINAL POSITION COORDS AT IONISATION IN /CNEUT/  *
C  *          YATIZS )                                                 *
C  *          KATIZS )                                                 *
c  *          SATIZS )                                                 *
C  *          VINS     LAUNCH VELOCITY W.R.T Y DIRECTION - IN /CNEUT/  *
C  *          LATIZ    FIRST POSITION IN IONIZATION ARRAY TO FILL      *
C  *          NATIZ    NUMBER OF IONISED NEUTRAL FRAGMENTS ON RETURN   *
C  *          RSTRUK   TOTAL OF FRAGMENTS STRIKING TARGET              *
C  *          RATIZ    TOTAL OF IONISED NEUTRAL FRAGMENTS              *
C  *          RNEUT    TOTAL OF LAUNCHED NEUTRAL FRAGMENTS             *
C  *          RWALLN   TOTAL OF FRAGMENTS PLATING OUT ON WALLS         *
C  *          RCENT    TOTAL OF FRAGMENTS REACHING CENTRAL MIRROR      *
C  *          RTMAX    TOTAL OF FRAGMENTS EXISTING AT TMAX             *
C  *          RFAIL    TOTAL OF FRAGMENTS WITH LAUNCH FAILURES         *
C  *          RMAIN    TOTAL OF FRAGMENTS REACHING MAIN PLASMA         *
C  *                                                                   *
C  *                                      C.M.FARRELL   MARCH 1989     *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE    "PARAMS"
      include    'params'
C     INCLUDE    "DYNAM1"
      include    'dynam1'
C     INCLUDE    "DYNAM3"
      include    'dynam3'
C     INCLUDE    "DYNAM4"
      include    'dynam4'
C     INCLUDE    "CGEOM"
      include    'cgeom'
C     INCLUDE    "COMTOR"
      include    'comtor'
C     INCLUDE    "CIONIZ"
      include    'cioniz'
C     INCLUDE    "CNEUT"
      include    'cneut'
      include    'cneut2'
C     INCLUDE    "CYIELD"
      include    'cyield'
C     INCLUDE    "CRAND"
      include    'crand'
c
      include    'commv'
      include    'cadas'
      include    'printopt'
c
      include    'fperiph_com'
c
      include    'line_profile'
c slmod begin
      include    'slcom'
c slmod end
c
      include    'hc_global_opts'
C
      REAL    XTOT(2),YTOT(2),AATIZ(2),RATIZ(2)
      REAL    ATOT(2),VTOT(2),VTOTM(2),VTOTA(2),VTOTAM(2),RTMAX(2)
      REAL    VMULTT(2),VMULTM(2),ETOT(2),MFRAC(2),XATIZ(2),YATIZ(2)
      REAL    VATIZ(2),VATIZM(2),TATIZ(2),EATIZ(2),REJECT(2),KATIZ(2)
      REAL    RMAIN(2),REXIT(2),RSTRUK(2),RFAIL(2),RCENT(2),RWALLN(2)
      REAL    MTC_RSTRUK(2),MTC_RWALLN(2)
      REAL    KMATIZ(2),RNEUT(2),KSATIZ(2),SFRAC(2),XSATIZ(2),YSATIZ(2)
      REAL    SSATIZ(2)
      real    xatiz2(2),yatiz2(2),ratiz2(2)
      real    lionizdat(2,2,2,2,5)
      real    move_factor
c
c     These are used to accumulate the average temperature of the
c     ions over multiple launch cycles.
c
      real    eatiztot, satiztot,temav
      data    eatiztot /0.0/
      data    satiztot /0.0/
c
c
c     Moved to velocity_dist module
c
c     Add velocity and angular distribution diagnostics
c     - vin
c     - anglan
c     - tanlan
c     - anglan+tanlan
c     - psi, beta - where appropriate
c     
c      integer  maxangs, ia, iv, flag
c      real    vin_scale,eiv
c      real,allocatable :: vin_dist(:),psi_dist(:),beta_dist(:),
c     >            anglan_dist(:),tanlan_dist(:),angtot_dist(:)  
c
c
c      logical sect,cflrin,cflrex
c
      logical sect
      integer intersect_result
C
      REAL    R,Z,RANMAX,SPUTY
      REAL    ANGLE,BETA,PSI,VMULT,TANLAN,TANGNT,ANGLAN,VIN,YMF
      REAL    TEMN,RAN,STATIM,ZA02AS,TANTRU,RPROD,K
      REAL    CISTOT,CISMAX,RSTMAX,TSTEPN,XVELF,YVELF,S,SMAX
      INTEGER IPROD,IK,IR,IFATE,IPOS,IT,NREJEC,KK,KKLIM,M,JK,JR,JD,ID
      integer iklast, irlast
      INTEGER IX,IY,J,IS,in,isol,ifp,irflct,idstart,iwstart
      INTEGER TMPI,TMPJ
      CHARACTER FATE(6)*14
      DOUBLE PRECISION DSPUTY
c
      real*8 cist
      real*8 mtccist
      integer mtccnt
      real    tmpvelf,xvelftmp,yvelftmp,vintmp
      real    tmpsum
C
      REAL TDUM(MAXPTS),XDUM(MAXPTS),YDUM(MAXPTS),WORK(4*MAXPTS)
      REAL RESULT,ROLD,ZOLD

      INTEGER INDWORK(2,MAXPTS),NRFCNT,INDI,IND,TOTRF,MAXRF
      real nrfloss,refprob
      REAL RNEW,ZNEW,TNEW,tnorm,orgr,orgz,rtmp,ztmp
c
      real err_ok,err_out
c
      integer maxnrfcnt
      parameter (maxnrfcnt=500)
c
      logical griderr,reflect_neut
c
c     Variables related to FAR PERIPHERY NEUTRAL IONIZATION option 
c
      INTEGER   RES,FPERIPH
      EXTERNAL  FPERIPH
c
c     Thompson calculation
c
      real find_thompson_velocity
      external find_thompson_velocity  
c
      real*8 cistfp
      real fpdist,fplosstim,fptmax,diffr,xstart
      real fp_lamiz,fp_sigmav,fp_vin
      integer fp_ik,fp_ir,fp_in
      real fp_rc,fp_zc
c slmod begin
      integer i
c slmod ned
c
c     HC related variables
c
      integer hc_wbc_flag
C
      DATA FATE /'REACHED WALL',           'REACHED CENTRE',
     >           'TIME = TMAX',            'STRUCK TARGET',
     >           'IONISED TO 1',           'FAILED LAUNCH'/
C
c     INITIALIZATION
c

c slmod begin
      target_loss = 0

      IF (ALLOCATED(wall_flx)) wall_nlaunch = wall_nlaunch + 1  ! tmp -- cleaner if wall_imp inintialization moved from this routine...

      IF (sloutput) THEN
        WRITE(0,*)
        WRITE(0,*) '***** HERE IN LAUNCH ! *****',cneutb
        WRITE(0,*)
      ENDIF
c slmod end
      STATIM = ZA02AS (1)
c
      DEBUGN = .FALSE.
      IF (CSTEPN.GT.0.0) THEN
        WRITE (6,9004) NINT(CSTEPN),FSRATE
        DEBUGN = .TRUE.
      ENDIF
c
      CISTOT = 0.0
      CISMAX = 0.0
      TOTRF  = 0
      MAXRF  = 0
      NRFLOSS = 0.0 
c
      err_ok = 0.0
      err_out = 0.0 
C
C-----------------------------------------------------------------------
C        SET INITIAL VALUES FOR DIAGNOSTICS
C-----------------------------------------------------------------------
C
C---- XTOT  : TOTAL OF X POSITIONS AT PRODUCTION
C---- YTOT  : TOTAL OF Y VALUES AT PRODUCTION
C---- ATOT  : TOTAL OF ANGLES AT PRODUCTION, CALCULATED AS IF ALL
C----         NEUTRALS WERE LAUNCHED FROM +Y REGION.
C---- VTOTA : TOTAL OF VELOCITIES AT PRODUCTION WITHOUT VMULT FACTORS
C---- VTOTAM: MAXIMUM VELOCITY AT PRODUCTION WITHOUT VMULT FACTOR
C---- VTOT  : TOTAL OF VELOCITIES AT PRODUCTION
C---- VTOTM : MAXIMUM VELOCITY AT PRODUCTION
C---- VMULTT: TOTAL OF VELOCITY ANGLE MULTIPLIERS "VMULT"
C---- VMULTM: MAXIMUM OF VELOCITY ANGLE MULTIPLIERS "VMULT"
C---- ETOT  : TOTAL OF TEMPERATURES AT PRODUCTION
C---- RWALLN: NO OF NEUTRALS REACHING X=-AW
C---- RCENT : NO OF NEUTRALS REACHING X=A
C---- RTMAX : NO OF NEUTRALS EXISTING AT T=TMAX
C---- MFRAC : NO OF IONISATIONS THAT OCCUR INBOARD
C---- KMATIZ: TOTAL OF K POSITIONS FOR INBOARD IONISATIONS
C---- SFRAC : NO OF IONISATIONS THAT OCCUR IN SOL+TRAP
C---- KSATIZ: TOTAL OF K POSITIONS FOR SOL IONISATIONS
C---- KATIZ : TOTAL OF K POSITIONS FOR ALL IONISATIONS
C---- XATIZ : TOTAL OF X POSITIONS AT IONISATION
C---- YATIZ : TOTAL OF ABS(Y) AT IONISATION
C---- AATIZ : TOTAL OF ANGLES AT IONISATION
C---- VATIZ : TOTAL OF VELOCITIES AT IONISATION
C---- VATIZM: MAXIMUM      "      "       "
C---- EATIZ : TOTAL OF TEMPERATURES AT IONISATION
C---- TATIZ : TOTAL OF TIMES TO IONISATION
C---- REJECT: NUMBER OF VELOCITIES GREATER THAN MAX RANDOMS
C---- RNEUT : SUM OF FRAGMENTS LAUNCHED ) SET TO SMALL NO. TO PREVENT
C---- RATIZ : SUM OF FRAGMENTS IONISED  ) POSSIBLE DIVIDE BY ZERO LATER
C---- RSTRUK: NUMBER OF FRAGMENTS STRIKING TARGET
C---- RFAIL : NUMBER OF FAILED LAUNCHES  V > VMAX AT LEAST 1000
C----         TIMES, WHICH SUGGESTS THAT RMAX IS VERY NEAR ZERO!
C---- RMAIN : NUMBER ENTERING MAIN PLASMA
C
      call rzero(lionizdat,2*2*2*2*5) 
c
      DO 50 M = 1, 2
        XTOT  (M) = 0.0
        YTOT  (M) = 0.0
        ATOT  (M) = 0.0
        VTOT  (M) = 0.0
        VTOTM (M) = 0.0
        VTOTA (M) = 0.0
        VTOTAM(M) = 0.0
        VMULTT(M) = 0.0
        VMULTM(M) = 0.0
        ETOT  (M) = 0.0
        RWALLN(M) = 0.0
        MTC_RWALLN(M) = 0.0
        RCENT (M) = 0.0
        RTMAX (M) = 0.0
        MFRAC (M) = 0.0
        KMATIZ(M) = 0.0
        SFRAC (M) = 0.0
        KSATIZ(M) = 0.0
        XSATIZ(M) = 0.0
        YSATIZ(M) = 0.0
        SSATIZ(M) = 0.0
        KATIZ (M) = 0.0
        XATIZ (M) = 0.0
        YATIZ (M) = 0.0
        AATIZ (M) = 0.0
        VATIZ (M) = 0.0
        VATIZM(M) = 0.0
        TATIZ (M) = 0.0
        EATIZ (M) = 0.0
        REJECT(M) = 0.0
        RNEUT (M) = LO
        RATIZ (M) = LO
        RSTRUK(M) = 0.0
        MTC_RSTRUK(M) = 0.0
        RFAIL (M) = 0.0
        RMAIN (M) = 0.0
        REXIT (M) = 0.0
c
        XATIZ2(M) = 0.0
        YATIZ2(M) = 0.0
        RATIZ2(M) = LO
c
   50 CONTINUE
C
C---- SET  RSTMAX: MAX NUMBER OF ITERATIONS UP TO TIME = TMAX (0.1S)
C---- (DIFFERS FROM CSTMAX WHICH APPLIES TO IONS,  BY FSRATE/QTIM AND
C---- BECAUSE TMAX IS 0.1 HERE BUT 10 FOR IONS).
C---- SET NATIZ : NO OF NEUTRAL "FRAGMENTS" THAT GET IONISED
C---- SET RPROD : SUM OF FRAGMENTS PRODUCED
C
      RSTMAX = 0.1 / FSRATE
      NATIZ  = 0
      RPROD  = 0.0
C
C-----------------------------------------------------------------------
C      LOOP FOR EACH NEUTRAL FRAGMENT TO BE LAUNCHED .....
C-----------------------------------------------------------------------
C
C---- SET START COORDINATES, POINTERS, ETC...  UPDATE EROSION DISTRIB'N
C
      CALL SURAND (SEED, NPROD-LPROD+1, RANVA)
      CALL SURAND (SEED, NPROD-LPROD+1, RANVB)
      CALL SURAND (SEED, NPROD-LPROD+1, RANVC)
      NRAND = NRAND + 3 * (NPROD-LPROD+1)
      KK    = 1000 * ISECT
      KKLIM = KK - 10
C
c slmod begin
      if (sloutput) write(0,*) 'debug: NPROD-LPROD+1 = ',NPROD-LPROD+1
c slmod end
      DO 900 IPROD = 1, NPROD-LPROD+1
c
c     jdemod - temp debug
c     
c         if (iprod.eq.62) then 
c            cstepn = 1
c            debugn = .true.
c         else
c            cstepn = 0
c            debugn = .false.
c         endif
c
c

        IF (sloutput.AND.grdnmod.NE.0.AND.
     .      MOD(iprod,MAX(1,(NPROD-LPROD+1)/10)).EQ.0)   ! sltmp
     .    WRITE(0,*) 'debug: iprod',iprod,NPROD-LPROD+1

        TSTEPN= CSTEPN
        IF (DEBUGN) WRITE (6,9005)
        R     = XPRODS(IPROD+LPROD-1)
        Z     = YPRODS(IPROD+LPROD-1)

c slmod begin - tmp
        IF (sloutput) THEN
          IF (GRDNMOD.NE.0.AND.STOPOPT.LT.2000) THEN        ! sltmp
            IF (stopopt.LE.0) stopopt = 1
c           IF (STOPOPT.LT.MAXNWS) THEN
c            WRITE(0 ,*) 'DEBUG: WALKS',iprod,STOPOPT
c             WRITE(50,*) 'DEBUG: WALKS',iprod,STOPOPT
            WALKS(STOPOPT,1) = HI ! 100.0 * RMAX
            WALKS(STOPOPT,2) = HI ! 100.0 * ZMAX
            STOPOPT = STOPOPT + 1
            WALKS(STOPOPT,1) = R
            WALKS(STOPOPT,2) = Z
            WALKS(STOPOPT+1,1) = HI
            WALKS(STOPOPT+1,2) = HI
            STOPOPT = STOPOPT + 1
          ENDIF
        ENDIF
c slmod end   
c
        orgr = r 
        orgz = z
c
        IK    = 0
        IR    = 0
C
c       Zero Momentum Transfer Collision (MTC) Counter
c
        mtccnt = 0
c
c       IPP/08 Krieger - variable "is" is not initialized under
c       certain conditions -> run time error further down in
c       write statement for Intel compiler
        is = 0
c
c
c       Assign some initial values depending on launch option
c
        IF (CNEUTB.EQ.0.or.cneutb.eq.3) THEN
c
c         ik,ir - indices of nearest bin centre
c         id is the inex to the theta bin from which the launch
c         is occurring
c
          ID = IDPRODS(IPROD+LPROD-1)
          IK = IKDS(ID)
          IR = IRDS(ID)
c
c         Add to wall erosion
c
c          WRITE(0,*) 'what',id,wallindex(id)
 
          if (wallindex(id).ne.0.0) then  
             wallse(wallindex(id)) = wallse(wallindex(id)) 
     >                             +sputys(iprod+lprod-1)
c slmod begin
             if (.not.allocated(wall_flx)) then
               wall_n       = wallpts
               wall_nlaunch = 1
               allocate(wall_flx(wall_n))
               do i = 1, wall_n
                 wall_flx(i)%launch = 0.0
               enddo
             endif
             i = wallindex(id)
             wall_flx(i)%launch(wall_nlaunch) = 
     .         wall_flx(i)%launch(wall_nlaunch) + sputys(iprod+lprod-1)
c             WRITE(0,*) i,wall_nlaunch
c slmod end
          else 
             write(6,'(a,5i5,3(1x,g12.5))') 'Wallse:Target?:',
     >              id,ik,ir,
     >              iprod,lprod,
     >              wallindex(id)
             wallse(maxpts+1) = wallse(maxpts+1) 
     >                             +sputys(iprod+lprod-1)


          endif
c
          griderr = .false.
c
c          write (6,*) 'cneut:',iprod,r,z,id,ik,ir
Cw
       elseif (cneutb.eq.1.or.cneutb.eq.5.or.
     >         cneutb.eq.6.or.cneutb.eq.7) then
c
c         Since this is a launch in free space ... ik,ir are the
c         indices of the closest bin centre ... id contains the
c         value 0 and is not used.
c
          ID = IDPRODS(IPROD+LPROD-1)
          call gridpos(ik,ir,r,z,.true.,griderr)
c
c         Record free space launched neutrals in wallse(maxpts+1) so 
c         the total "erosion" source is available.
c
          wallse(maxpts+1) = wallse(maxpts+1) + sputys(iprod+lprod-1) 
c
c          IK = IKXYS(IX,IY)
c          IR = IRXYS(IX,IY)
c
       ELSEIF (CNEUTB.EQ.2.or.cneutb.eq.4) THEN
C
C         These quantities contain the IK,IR coordinates of the
C         grid point nearest the launch position.
C         ID contains the index of the wall point and IS is the
C         indicator of whether the launch is from the clockwise
C         or anti-clockwise wall segment coming from that wall point.
C
          ID = IDPRODS(IPROD+LPROD-1)
          IS = ISPRODS(IPROD+LPROD-1)
          call gridpos(ik,ir,r,z,.true.,griderr)
c
c         Record the particle in the Erosion array.
c
          if (id.lt.1.or.id.gt.wallpts) then  

             write(6,'(a,6i5)') 'Wallse:Target?:',id,
     >              ik,ir,
     >              iprod,lprod,
     >              wallindex(id)
             wallse(maxpts+1) = wallse(maxpts+1)+sputys(iprod+lprod-1)
          else
             wallse(id) = wallse(id) + sputys(iprod+lprod-1)
c slmod begin
             if (.not.allocated(wall_flx)) then
               wall_n       = wallpts
               wall_nlaunch = 1
               allocate(wall_flx(wall_n))
               do i = 1, wall_n
                 wall_flx(i)%launch = 0.0
               enddo
             endif
             i = id
             wall_flx(i)%launch(wall_nlaunch) = 
     .         wall_flx(i)%launch(wall_nlaunch) + sputys(iprod+lprod-1)
c             WRITE(0,*) i,wall_nlaunch
c slmod end
          endif

c
c          IK = IKXYS(IX,IY)
c          IR = IRXYS(IX,IY)
c
        ENDIF
c
c       Assign starting index - either target or wall coordinate depending on launch
c
        idstart = id
c
c       Assign starting wall index 
c
        if (neuttype.eq.1.or.neuttype.eq.2.or.neuttype.eq.3) then 
           iwstart = wallindex(idstart)
        elseif (neuttype.eq.4.or.neuttype.eq.5) then 
           iwstart = idstart
        else
           iwstart = -1
        endif

c
c       Use the nrs+1 element to record the total source from 
c       all elements. The portion that is ionized is recorded
c       in the rest of the array. 
c
c slmod begin
c        if (grdnmod.eq.0.or.neuttype.ne.0) then       
c
c jdemod - changed if clause
        if (neuttype.ne.0) then       
          if (idstart.ge.1.and.idstart.le.wallpts) then  
c
             wtsource(idstart,nrs+1,1,neuttype)=
     >                    wtsource(idstart,nrs+1,1,neuttype)
     >                     +SPUTYS(IPROD+LPROD-1)
c
          endif
        endif
c
cc        if (idstart.ge.1.and.idstart.le.wallpts) then  
cc
c           wtsource(idstart,nrs+1,1,neuttype)=
c     >                  wtsource(idstart,nrs+1,1,neuttype)
c     >                   +SPUTYS(IPROD+LPROD-1)
cc
c        endif
c slmod end
C
C         WRITE(6,*) 'NEUTRAL LAUNCH:',R,Z,IX,IY,ID
C         WRITE(6,*) 'INFO:',(WALLPT(ID,J),J=1,13)
C
c 
c        IF (griderr) THEN
c          WRITE(6,*) 'LAUNCH ERROR:',IX,IY,ID,R,Z,IK,IR
C
C         ALL CALCULATED LAUNCH POSITIONS ARE VALID - FOR LAUNCHES
C         FROM WALLS BASED ON THE OUTER RINGS - THIS IS NOT NECESSARILY
C         SO FOR "TRUE" WALL LAUNCHES.
C
C         IFXYS(IX,IY) = 1
c        ENDIF
C
c
c slmod begin
        if (cgridopt.eq.LINEAR_GRID.or.cgridopt.eq.RIBBON_GRID.or.
     >      nbr.gt.0) then
          if (ikds(id).eq.1) then
            M = 2                   ! I think this is right...
          else
            M = 1
          endif
        elseif (cgridopt.eq.0.or.cgridopt.eq.1.or.cgridopt.eq.3) then
c
c        if (cgridopt.eq.0.or.cgridopt.eq.1.or.cgridopt.eq.3) then
c slmod end
c
          M = 2
c
          IF ((CNEUTB.EQ.0.or.cneutb.eq.3).and.ID.LE.NDSIN) THEN
            M =1
          ELSEIF ((cneutb.eq.2.or.
     >             cneutb.eq.4.or.cneutb.eq.5)
     >           .AND.(KSS(IK,IR).GT.KSMAXS(IR)/2.0)) THEN
            M = 1
c
          elseif ((cneutb.eq.1.or.cneutb.eq.6.or.cneutb.eq.7).and.
     >             orgr.lt.r0) then 
c
            M = 1
c 
          ENDIF
c
        elseif (cgridopt.eq.2) then
c
          M = 2
          IF ((CNEUTB.EQ.0.or.cneutb.eq.3).and.
     >        (ID.gt.NDSIN.and.id.le.ndsin3)) THEN
            M =1
c
c         WARNING: NEED TO VERIFY THAT kss < 1/2 ksmaxs for the
c                  ir < irsep region corresponds to the inner
c                  region.
c
          ELSEIF ((cneutb.eq.2.or.
     >             cneutb.eq.4.or.cneutb.eq.5)
     >           .AND.((ir.ge.irsep.and.ir.le.irwall2).or.
     >                (((ir.ge.irtrap2.and.ir.le.nrs2).or.
     >                 (ir.ge.irtrap.and.ir.le.nrs) .or.
     >                 (ir.lt.irsep))
     >                .and.kss(ik,ir).lt.ksmaxs(ir)/2.0))) then
            M = 1
c
          elseif ((cneutb.eq.1.or.cneutb.eq.6.or.cneutb.eq.7).and.
     >             orgr.lt.r0) then 
c
            M = 1
c 
          ENDIF
c
        endif
c
        RANMAX= RMAXS(IPROD+LPROD-1)
        SPUTY = SPUTYS(IPROD+LPROD-1)
        DSPUTY= DBLE (SPUTY)
        RPROD = RPROD + SPUTY
        CIST  = 0.0
        MTCCIST = 0.0
c        IT    = IPOS (CIST, CTIMES(1,0), NTS+1)
        IT    = IPOS (real(CIST), CTIMES(1,0), NTS+1)

c
c Added recording of erosion in NEROS array for all launch options (?)
c
! ammod begin.
c        IF (CNEUTB.EQ.0.or.cneutb.eq.3) NEROS(ID,3)=NEROS(ID,3)+SPUTY
c
        IF (CNEUTB.EQ.0.or.cneutb.eq.3) Then
	   NEROS(ID,3)=NEROS(ID,3)+SPUTY
	ElseIf (CNeutB .eq. 2 .or. CNeutB .eq. 4) Then
	   ! Add removal from wall segment if it is also a target segment.
           ! Note: WALLPT is an array of type REAL.
	   If (INT (wallpt (ID,18)) .gt. 0) Then
	      ! Add to total removal from target.
	      NEROS (INT (wallpt (ID,18)),3) = 
     >          NEROS (INT (wallpt (ID,18)),3) + SPUTY
	   End If
	End If
! ammod end.

        CFLRIN= .TRUE.
        CFLREX= .TRUE.
        K     = KKS(IR)
        S     = KSS(IK,IR)
        SMAX  = KSMAXS(IR)


C
C-----------------------------------------------------------------------
C      SELECT ANGLE FOR LAUNCH ...
C      PETER'S NOTES 38,65,83,93,109,218
C-----------------------------------------------------------------------
C
c slmod begin
c...    I guess I am using an old vel/ang flag (3), since BETA and PSI 
c       are not assigned but are sent to RECORD_VDIST below:
        BETA = 0.0
        PSI  = 0.0
c slmod end

        IF     (CNEUTC.EQ.0) THEN
          ANGLAN = SIGN (ASIN (RANVA(IPROD)), RANVB(IPROD)-0.5)
        ELSEIF (CNEUTC.EQ.3.OR.CNEUTC.EQ.5.OR.CNEUTC.EQ.9) THEN
          ANGLAN = SIGN (ASIN (SQRT(RANVA(IPROD))), RANVB(IPROD)-0.5)
        ELSEIF (CNEUTC.EQ.1.OR.CNEUTC.EQ.2.OR.CNEUTC.EQ.4
     >          .or.cneutc.eq.14.or.cneutc.eq.16) THEN
          BETA   = ASIN (SQRT (RANVA(IPROD)))
          PSI    = 2.0 * PI * RANVB(IPROD)
          ANGLAN = ATAN (TAN (BETA) * COS (PSI))
        ELSEIF (CNEUTC.EQ.6.or.cneutc.eq.17.or.cneutc.eq.18) THEN
          ANGLAN = 0.0
        ELSEIF (CNEUTC.EQ.7.OR.CNEUTC.EQ.11) THEN
          ANGLAN = SIGN (ACOS ((1.0-RANVA(IPROD)) ** (1.0/3.0)),
     >                   RANVB(IPROD)-0.5)
        ELSEIF (CNEUTC.EQ.8.or.cneutc.eq.15) THEN
          ANGLAN = 2.0 * PI * RANVA(IPROD) - PI
        ELSEIF (CNEUTC.EQ.10) THEN
          BETA   = ACOS ((1.0-RANVA(IPROD)) ** (1.0/3.0))
          PSI    = 2.0 * PI * RANVB(IPROD)
          ANGLAN = BETA
          IF (PSI.LT.PI/2.0 .OR. PSI.GT.3.0*PI/2.0) ANGLAN = -BETA
        ELSEIF (CNEUTC.EQ.12.OR.CNEUTC.EQ.13) THEN
          BETA   = ACOS ((1.0-RANVA(IPROD))**(1./(CNIN+1.)))
          PSI    = 2.0 * PI *RANVB(IPROD)
          ANGLAN = ATAN (TAN (BETA) * COS (PSI))
c
c       3D isotropic
c
        elseif (cneutc.eq.19) then 
          BETA   = PI * RANVA(IPROD) - PI/2.0 
          PSI    = 2.0 * PI * RANVB(IPROD)
          ANGLAN = ATAN (TAN (BETA) * COS (PSI))
        ENDIF
C
        IF (CNEUTC.EQ.1.OR.CNEUTC.EQ.4.OR.CNEUTC.EQ.12.OR.
     >      CNEUTC.EQ.13.or.cneutc.eq.14.or.cneutc.eq.16.or.
     >      cneutc.eq.19) THEN
          VMULT = SQRT (ABS(COS(BETA)**2+(SIN(BETA)**2*COS(PSI)**2)))
        ELSE
          VMULT = 1.0
        ENDIF

C
C-----------------------------------------------------------------------
C     CALCULATE NORMAL TO TARGET AND ADD TO ANGLE
C     FOR SIDE PUFF GAS INJECTION, PUFF ALONG Y=0, IE RESET TANLAN TO 0
C     THIS MAY RESULT IN SOME NEUTRALS BEING LAUNCHED STRAIGHT AT THE
C     TARGET SURFACE: COUNT THESE BUT DO NOT FOLLOW THEM.
C     VARIABLES ANGLAN & TANLAN REFER TO +Y REGION ONLY.  TRUE VALUES
C     USED ARE ANGLE & TANGNT WHICH CORRECT FOR -Y OR +Y REGIONS.
C-----------------------------------------------------------------------
C
        if (cneutb.eq.1.or.cneutb.eq.5.or.
     >      cneutb.eq.6.or.cneutb.eq.7) then
c
          tantru = 0.0
c
        ELSEIF (CNEUTB.EQ.2.or.cneutb.eq.4) THEN
C
C         THE IS VARIABLE CONTAINS THE INFORMATION
C         ABOUT WHICH SIDE OF A WALL LAUNCH POINT THE
C         THE PARTICLE IS LAUNCHED FROM
C
           IF (IS.EQ.0) THEN
C
C             ANTI-CLOCKWISE
C
              TANTRU = WALLPT(ID,8) + 0.5 * PI
           ELSEIF (IS.EQ.1) THEN
C
C              CLOCKWISE
C
              TANTRU = WALLPT(ID,9) - 0.5 * PI
           ENDIF
        ELSEif (cneutb.eq.0.or.cneutb.eq.3) then 
          TANTRU = THETAS(ID)
        ENDIF
C
        IF (CNEUTE.EQ.0.or.(cneute.eq.1.and.status.gt.1)) THEN
          TANLAN = TANTRU
        else 
          TANLAN = CSNORM
        ENDIF
C
        ANGLE  = ANGLAN
        TANGNT = TANLAN
C
        XTOT (M) = XTOT(M) + R * SPUTY
        YTOT (M) = YTOT(M) + Z * SPUTY
        ATOT (M) = ATOT(M) + (ANGLAN+TANLAN) * SPUTY


        RNEUT(M) = RNEUT(M) + SPUTY

c
C
C-----------------------------------------------------------------------
C     CALCULATE LAUNCH VELOCITY.  FIRST SELECT RANDOM NUMBER IN THE
C     RANGE 0 < RAN < RANMAX.  IF NOT FOUND, REJECT NUMBER AND FIND
C     ANOTHER ONE ... TO PREVENT A POTENTIALLY INFINITE LOOP HERE, COUNT
C     NUMBER OF REJECTED VELOCITIES IN NREJEC AND LIMIT TO (SAY) 1000;
C     IF THIS LIMIT EXCEEDED THEN NEUTRAL IS CALLED A "FAILED LAUNCH"
C-----------------------------------------------------------------------
C
c       Set energy of sputtered particle based on input value 
c
        if (eprods(iprod).gt.0.0) then 

            VIN = 1.38E4 * SQRT (eprods(iprod)/CRMI)

        else
c
c         Set energy of sputtered particle based on specified Vel/Angle
c         flag. 
c

          RAN = RANVC(IPROD)
          NREJEC = 0
c slmod begin -25/09/2009
c...      Not looking in detail at why RAN must be less than RANMAX, or
c         how to avoid this in the code, but the 1000 iteration limit
c         is turning out to be a problem for Fe sputtering in the ITER
c         simluation, since RANMAX is on the order of 1E-03 and so a lot
c         of particles can be lost simply from this arbitrary 1000 iteration
c         limit.  So, increasing this based on the value of RANMAX (there's
c         some headroom here anyway, surely, with modern computers):
          ITERATION_LIMIT = MAX(1000,NINT(1.0/RANMAX)*1000)
          IF (.NOT.ITERATION_WARNING.AND.ITERATION_LIMIT.GT.1E6) THEN
            ITERATION_WARNING = .TRUE.
            WRITE(0,*) 
            WRITE(0,*) '*****************************************'
            WRITE(0,*) '* NEUT ITERATION LIMIT GREATER THAN 1E6 *'
            WRITE(0,*) '*****************************************'
            WRITE(0,*) 'INTERATION_LIMIT=',ITERATION_LIMIT
            WRITE(0,*) 'RANMAX          =',RANMAX
          ENDIF
c slmod end
  100     CONTINUE
          IF (RAN.GT.RANMAX) THEN
            REJECT(M) = REJECT(M) + SPUTY
            NREJEC = NREJEC + 1
            CALL SURAND2 (SEED, 1, RAN)
            NRAND = NRAND + 1  
c slmod begin
c...        See above:
            IF (NREJEC.LT.ITERATION_LIMIT) GOTO 100
c
c            IF (NREJEC.LT.1000) GOTO 100
c slmod end
            VIN   = 0.0
            TEMN  = 0.0
            RFAIL(M) = RFAIL(M) + SPUTY
            IFATE = 6
! ammod begin.       
	    ! WBC comparison addition for neut (charge state 0) failed launch.
            Call global_hc_wbc_comp(0,CRMI,VIN,TEMN,SPUTY)
! ammod end.	       


            GOTO 899
          ENDIF
C
          IF(CNEUTC.EQ.0.OR.CNEUTC.EQ.1.OR.CNEUTC.EQ.4.OR.
     >       CNEUTC.EQ.5.or.cneutc.eq.17)THEN
            VIN = 1.38E4 * SQRT (CEBD/(1.0/SQRT(RAN)-1.0)/CRMI) * VMULT
          ELSEIF (CNEUTC.EQ.2) THEN
            IF (RAN.GE.1.0) RAN = 0.999999
            VIN = 1.38E4 * SQRT (CTEM1*ABS(LOG(1.0-RAN))/CRMI)
          ELSEIF (CNEUTC.EQ.3.OR.CNEUTC.EQ.6.OR.CNEUTC.EQ.7.OR.
     >            CNEUTC.EQ.8.OR.CNEUTC.EQ.10.OR.CNEUTC.EQ.11) THEN
            VIN = 1.38E4 * SQRT (CTEM1/CRMI)
C
C------ NOTE 156 VEL/ANG9 FLAG.  NOTICE THAT LAUNCHES ALTERNATELY AT
C------ EIN1, EIN2,  BUT WE ALSO HAVE THE +/-Y THING.  HENCE HERE
C------ WE NEED TO LAUNCH 2 PARTICLES AT EIN1, FOLLOWED BY 2 AT EIN2,
C------ ETC TO ENSURE THAT SOME OF EACH ARE LAUNCHED ON EACH SIDE
C------ OF Y = 0.
C
          ELSEIF (CNEUTC.EQ.9) THEN
            IF (2*(IPROD/4).EQ.IPROD/2) THEN
              VIN = 1.38E4 * SQRT (CTEM1/CRMI)
            ELSE
              VIN = 1.38E4 * SQRT (CTEM2/CRMI)
            ENDIF
          ELSEIF (CNEUTC.EQ.12.or.cneutc.eq.19) THEN
            VIN = 1.38E4 * SQRT (CTEM1/CRMI) * VMULT
          ELSEIF (CNEUTC.EQ.13) THEN
            IF (RAN.GE.1.0) RAN = 0.999999
            VIN = 1.38E4 * SQRT (CTEM1*ABS(LOG(1.0-RAN))/CRMI) * VMULT
          elseif (cneutc.eq.14) then
            VIN = 1.38E4 * SQRT (KTIDS(ID)/CRMI) * VMULT * CVAMULT
          elseif (cneutc.eq.15) then
            VIN = 1.38E4 * SQRT (KTIBS(IK,IR)/CRMI) * VMULT * CVAMULT
          elseif (cneutc.eq.16) then
            VIN = find_thompson_velocity(neuttype,idstart,seed,nrand,
     >             1)
     >             * VMULT * CVAMULT
          elseif (cneutc.eq.18) then
            VIN = find_thompson_velocity(neuttype,idstart,seed,nrand,
     >             1)
     >             * VMULT * CVAMULT
          ENDIF
c
        endif 
c
        VTOTA (M) = VTOTA(M) + VIN/VMULT * SPUTY
        VTOTAM(M) = MAX (VTOTAM(M), VIN/VMULT)
        VMULTT(M) = VMULTT(M) + VMULT * SPUTY
        VMULTM(M) = MAX (VMULTM(M), VMULT)
        VTOT  (M) = VTOT(M) + VIN * SPUTY
        VTOTM (M) = MAX (VTOTM(M), VIN)
C
C------ CALCULATE X,Y COMPONENTS OF VELOCITY, PROB OF IONISATION, ETC
C
        XVELF  = VIN * COS(ANGLE+TANGNT) * FSRATE
        YVELF  = VIN * SIN(ANGLE+TANGNT) * FSRATE
        TEMN   = CRMI * (VIN/1.38E4) * (VIN/1.38E4)
        ETOT(M)= ETOT(M) + TEMN * SPUTY
c
c       Record detailed velocity debugging information if the 
c       option is active
c       This is now in the velocity_dist module to allow code 
c       sharing between regular launch and hc launch code
c
c       The setup/allocate and printing/deallocate calls are now
c       made in neutbatch so that it can support either regular or 
c       hc launches.
c
        call record_vdist(vin,beta,psi,tanlan,anglan,sputy)

c
c
c
C-----------------------------------------------------------------------
c
c       Check if the FAR PERIPHERY NEUTRAL IONIZATION OPTION is active.
c
C-----------------------------------------------------------------------
c
c       If it is - check the Lambda_iz value for the species 
c       given the local plasma conditions specified for the FP
c       or the plasma conditions taken from the nearest plasma
c       cell on the grid to the start point. This option can
c       only apply to actual wall elements - if it is a target 
c       element the launch continues normally. If LAMBDA_IZ is 
c       greater than the width of the FP the code will also continue
c       with a normal launch.
c
c       After return from this code - the result will be 
c       1) the particle returns to the wall
c       2) the particle is lost to the FP target
c       3) the particle transports into the main plasma 
c
c       In the case of (3) - an R,Z location representative of the 
c       edge of the plasma cell associated with the wall segment is 
c       assigned and the code continues as if ionization had occured
c       at this location. 
c
c        if (fp_neut_opt.gt.0) then 
c
c        jdemod - added test for iwstart.ne.-1 since fp neutral ionization 
c                 only makes sense at the moment for particles associated with
c                 a wall element. The code would need to be modified to
c                 accommodate free space launches.
c
        if (fp_neut_opt.gt.0.and.iwstart.ne.-1) then 
c
c          Split evaluation into 2 IF statements since iwstart is not 
c          always properly defined.       
c
           if(
     >     (wallpt(iwstart,16).ne.1.and.wallpt(iwstart,16).ne.4))
     >     then 
c
           fp_ik = wallpt(iwstart,26)
           fp_ir = wallpt(iwstart,27)
c
c          Calculate the LAMBDA_IZ assuming the particle 
c          starts with a velocity associated with an 
c          and intial energy of CTEM1.
c
           fp_VIN = 1.38E4 * SQRT (CTEM1/CRMI)
c
           if (fp_plasma_opt.eq.0) then 
              fp_lamiz =  fp_vin * kfizs(fp_ik,fp_ir,0) 
              write(6,'(a,3i4,4(1x,g12.5))') 'FP_NEUT:0:',
     >                 iwstart,fp_ik,fp_ir,fp_lamiz

           elseif (fp_plasma_opt.eq.1) then 
c
c             CALL ADASRD(YEAR,CION,IZ+1,ICLASS,NPTS,TE,NE,SIGMA_V)
c
              CALL ADASRD(YEAR,CION,1,2,1,fp_te,knbs(fp_ik,fp_ir),   
     >                    fp_sigmav)                                      
c
c             Check for valid values
c
              if (knbs(fp_ik,fp_ir)*fp_sigmav.gt.0.0) then 
c
                 fp_lamiz = fp_vin / (knbs(fp_ik,fp_ir)*fp_sigmav)
c
              endif
c
c             Print out for debugging purposes
c
              write(6,'(a,3i4,4(1x,g12.5))') 'FP_NEUT:1:',
     >                 iwstart,fp_ik,fp_ir,fp_lamiz,
     >                 fp_sigmav,fp_te,knbs(fp_ik,fp_ir)
c
           elseif (fp_plasma_opt.eq.2) then 
c
c             CALL ADASRD(YEAR,CION,IZ+1,ICLASS,NPTS,TE,NE,SIGMA_V)
c
              CALL ADASRD(YEAR,CION,1,2,1,fp_te,fp_ne,   
     >                    fp_sigmav)                                      
c
c             Check for valid values
c
              if (fp_ne*fp_sigmav.gt.0.0) then 
c
                 fp_lamiz = fp_vin / (fp_ne*fp_sigmav)
c
              endif
c
c             Print out for debugging purposes
c
              write(6,'(a,3i4,4(1x,g12.5))') 'FP_NEUT:1:',
     >                 iwstart,fp_ik,fp_ir,fp_lamiz,
     >                 fp_sigmav,fp_te,fp_ne
c
           endif
c
c          Check to see if fp_lamiz is less than the width
c          of the far periphery zone - this can be different 
c          for INNER and OUTER. Use the S value of the associated
c          plasma cell in order to remain consistent with the ion
c          routines in DIV.
c
c          USE the ion quantum time step for following the ion -
c          convert back to neutral time when recording data for
c          this routine.(?)
c             
c          Outer or Inner FP 
c	   
           if (kss(fp_ik,fp_ir).lt.(ksmaxs(ir)/2.0)) then 
c                         
c             OUTER (XPT UP) / INNER (XPT DOWN)
c	   
              fpdist = fpxmaxO
              fplosstim = qtim/fptimO 
c	   
           else
c	   
c             INNER          
c	   
              fpdist = fpxmaxI
              fplosstim = qtim/fptimI 
c	   
           endif
c
c          Invoke far periphery starting at (fpdist-fp_lamiz) into the 
c          far periphery
c
           if (fp_lamiz.le.fpdist) then 
c
              fp_neut_ent = fp_neut_ent + sputy
              fp_neut_data(iwstart,1) = fp_neut_data(iwstart,1) 
     >                                  + sputy 
c
              DIFFR = SQRT(2.0*CDPERPFP*QTIM)
c
              xstart = fpdist - fp_lamiz 
c
c             Allow maximum following time the same as neutrals (0.1 seconds) 
c             Time in FPERIPH uses the ion quantum time step.  
c
              FPTMAX = 0.1 / QTIM       
c
              RES = FPERIPH(0.0d0,cistfp,FPDIST,fplosstim,FPTMAX,
     >                      NRAND,DIFFR,SEED,xstart)

c
c             The return codes from FPERIPH mean the following:
C             1) DIFFUSES BACK INTO THE PLASMA
C             2) EXCEEDS ION TIME-LIMIT
C             3) HITS THE WALL AT EDGE OF PERIPHERY
C             4) IS ASSUMED TO HIT TARGET PLATE IN FP REGION
c
c             Particle enters plasma - calculate effective ionization
c             position in cell (fp_ik,fp_ir)
c
c             Adjust time for using QTIM instead of FSRATE.
c
              cist = cist + cistfp * qtim/fsrate

              if (res.eq.1) then 
c
c                Set R,Z and record ionized particle information
c
                 fp_neut_ioniz = fp_neut_ioniz + sputy
                 fp_neut_data(iwstart,2) = fp_neut_data(iwstart,2) 
     >                                   + sputy 
c
c                Calculate R,Z values for ionization position. 
c
                 fp_in = korpg(fp_ik,fp_ir)        
c
                 if (fp_in.eq.0.or.nvertp(fp_in).eq.0) then 
c
c                   No polygon data - this case should not occur 
c 
                    R = rs(fp_ik,fp_ir) 
                    Z = zs(fp_ik,fp_ir)
c
c                Polygon data exists
c                 
                 else           
c
c                   Find midpoint of "boundary side" of polygon
c
                    if (fp_ir.eq.irwall-1) then   
c
c                      Vertices 2,3 define the outer boundary
c                   
                       fp_rc = (rvertp(2,fp_in)+rvertp(3,fp_in))/2.0
                       fp_zc = (zvertp(2,fp_in)+zvertp(3,fp_in))/2.0
c
c                      Put particle just inside the plasma boundary
c
                       r = rs(fp_ik,fp_ir)
     >                     + 0.98 * (fp_rc-rs(fp_ik,fp_ir))
                       z = zs(fp_ik,fp_ir)
     >                     + 0.98 * (fp_zc-zs(fp_ik,fp_ir))
c                     
                    elseif (fp_ir.eq.irtrap+1) then 
c
c                      Vertices 1,4 define the outer boundary
c
                       fp_rc = (rvertp(1,fp_in)+rvertp(4,fp_in))/2.0
                       fp_zc = (zvertp(1,fp_in)+zvertp(4,fp_in))/2.0
c
c                      Put particle just inside the plasma boundary
c
                       r = rs(fp_ik,fp_ir) 
     >                     + 0.98 * (fp_rc-rs(fp_ik,fp_ir))
                       z = zs(fp_ik,fp_ir) 
     >                     + 0.98 * (fp_zc-zs(fp_ik,fp_ir))
c                     
                    else
c
c                      Not in a boundary ring - this should not occur
c
                       R = rs(fp_ik,fp_ir) 
                       Z = zs(fp_ik,fp_ir)
c
                       write(6,'(a,3i5)') 'WARNING: FP ION NOT'//
     >                                ' IN BOUNDARY RING',fp_ik,fp_ir
c 
                    endif 
c
                 endif
c
c                Set up data - assign ik,ir appropriately
c
                 ir = fp_ir
                 ik = fp_ik
c
                 TEMN = CTEM1
                 S = kss(ik,ir)
                 smax = ksmaxs(ir)  
                 k = kks(ir)
                 vin = fp_vin
c
c                Copy code from regular ionization routine            
c
c
                 TIZS(IK,IR,0) = TIZS(IK,IR,0) + SPUTY
c
                 NATIZ = NATIZ + 1
                 XATIZS(LATIZ+NATIZ-1) = R
                 YATIZS(LATIZ+NATIZ-1) = Z
                 KATIZS(LATIZ+NATIZ-1) = K
                 SATIZS(LATIZ+NATIZ-1) = MIN (S, SMAX-S)
                 TEMTIZS(LATIZ+NATIZ-1) = TEMN
                 snews(latiz+natiz-1)  = sputy
                 cistizs(latiz+natiz-1)  = cist * fsrate / qtim
c
c                Record starting wall or target element for neutral/ion.
c
                 IDATIZS(LATIZ+NATIZ-1,1) = idstart
                 IDATIZS(LATIZ+NATIZ-1,2) = neuttype
c
                 launchdat(LATIZ+NATIZ-1,3) = nrfcnt
                 launchdat(LATIZ+NATIZ-1,2) = launchdat(iprod,2)
c
                 IF (CNEUTG.EQ.0) THEN
                    VINS(LATIZ+NATIZ-1) = 0.0
                 ELSEIF (CNEUTG.EQ.1.OR.CNEUTG.EQ.3) THEN
                    CALL SURAND2 (SEED, 1, RAN)
                    NRAND = NRAND + 1
                    VINS(LATIZ+NATIZ-1) = SIGN (0.5*VIN, RAN-0.5)
                    IF (CNEUTG.EQ.3) THEN
                       CALL SURAND2 (SEED, 1, RAN)
                       NRAND = NRAND + 1
                       VINS(LATIZ+NATIZ-1)=
     >                          VINS(LATIZ+NATIZ-1)*2.0*SQRT(RAN)
                    ENDIF
                 ELSEIF (CNEUTG.EQ.2) THEN
                    VINS(LATIZ+NATIZ-1) = SIGN (VIN, 0.5*SMAX-S)
                 ENDIF
c
c                Record Ring and target index from which ion originated.
c                for target or wall launches. 
c
                 if (neuttype.ne.0) then 

                    wtsource(idstart,ir,1,neuttype)=
     >                        wtsource(idstart,ir,1,neuttype)+sputy
                    wtsource(idstart,ir,2,neuttype)=
     >                        wtsource(idstart,ir,2,neuttype)
     >                      + satizs(latiz+natiz-1)*sputy
c
c                   Accumulate data on the ionization of eroded 
c                   particles from each wall element
c
                    wallse_i(iwstart) =  wallse_i(iwstart) + sputy
c
                 endif  
c
c                NOTE:  David Elder,    1995  AUG 3
c	         
c                There were too many new variables being required 
c                to keep track of the various ionization data under
c                all the varied circumstances - so the information
c                was centralized into one array which contains
c                all of the information - there are some other
c                variables which accumulate the data for 
c                compatibility with the older code - but all of 
c                the new code uses the new setup.
c	         
c                Here is a summary of the array and its contents
c   
c                IONIZDAT(2    ,2  ,2    ,2      ,5)
c                         isol  m   ifp   irflct  quant
c	         
c                isol = 1 for SOL information 
c                     = 2 for MAIN plasma information
c  	         
c                m    = 1 for target 1   ik > nks(ir)/2
c                     = 2 for target 2   ik =< nks(i2)/2 
c
c                ifp  = 1 for a neutral resulting from a regular
c                         launch    
c                     = 2 for a neutral resulting from a Far Periphery
c                         relaunch
c
c                irflct = 1 for a neutral that has NOT been reflected
c                       = 2 for a neutral that has been reflected 
c	         
c                quant= 1 total weight of neutrals
c                     = 2 weighted R coordinate
c                     = 3 weighted Z coordinate 
c                     = 4 weighted K value
c                     = 5 weighted S value
c
c
                 if (ir.ge.irsep) then 
                    isol = 1
                 else
                    isol = 2
                 endif
c	         
                 if (fpropt.eq.1.and.launchdat(iprod,2).eq.1.0) then 
                    ifp = 2
                 else 
                    ifp = 1
                 endif   
c	         
                 if ((nrfopt.eq.1.or.nrfopt.eq.2).and.nrfcnt.gt.0) then 
                    irflct = 2
                 else 
                    irflct = 1
                 endif 
c	         
                 ionizdat(isol,m,ifp,irflct,1) = 
     >                ionizdat(isol,m,ifp,irflct,1) +  SPUTY
                 ionizdat(isol,m,ifp,irflct,2) = 
     >                ionizdat(isol,m,ifp,irflct,2) +  SPUTY * R
                 ionizdat(isol,m,ifp,irflct,3) = 
     >                ionizdat(isol,m,ifp,irflct,3) +  SPUTY * Z
                 ionizdat(isol,m,ifp,irflct,4) = 
     >                ionizdat(isol,m,ifp,irflct,4) +  SPUTY * K
                 ionizdat(isol,m,ifp,irflct,5) = 
     >                ionizdat(isol,m,ifp,irflct,5) 
     >                +  SPUTY * min(s,smax-s)
c
                 lionizdat(isol,m,ifp,irflct,1) = 
     >                lionizdat(isol,m,ifp,irflct,1) +  SPUTY
                 lionizdat(isol,m,ifp,irflct,2) = 
     >                lionizdat(isol,m,ifp,irflct,2) +  SPUTY * R
                 lionizdat(isol,m,ifp,irflct,3) = 
     >                lionizdat(isol,m,ifp,irflct,3) +  SPUTY * Z
                 lionizdat(isol,m,ifp,irflct,4) = 
     >                lionizdat(isol,m,ifp,irflct,4) +  SPUTY * K
                 lionizdat(isol,m,ifp,irflct,5) = 
     >                lionizdat(isol,m,ifp,irflct,5) 
     >                +  SPUTY * min(s,smax-s)
c	         
                 IF (IR.GE.IRSEP) THEN
                   SFRAC (M) = SFRAC (M) + SPUTY
                   KSATIZ(M) = KSATIZ(M) + SPUTY * K
                   XSATIZ(M) = XSATIZ(M) + SPUTY * R
                   YSATIZ(M) = YSATIZ(M) + SPUTY * Z
                   SSATIZ(M) = SSATIZ(M) + SPUTY * MIN (S, SMAX-S)
                 ELSE
                   MFRAC (M) = MFRAC (M) + SPUTY
                   KMATIZ(M) = KMATIZ(M) + SPUTY * K
                 ENDIF
c	         
                 KATIZ(M) = KATIZ(M) + K * SPUTY
                 TATIZ(M) = TATIZ(M) + CIST * FSRATE * SPUTY
                 XATIZ(M) = XATIZ(M) + R * SPUTY
                 YATIZ(M) = YATIZ(M) + Z * SPUTY
c	         
c                Record additional time information
c
                 cieizs(0) = cieizs(0) + cist * sputy
                 citizs(0) = citizs(0) + sputy
c	         
c slmod begin - ribbon dev
                 if (cgridopt.eq.0.or.cgridopt.eq.1.or.
     >               cgridopt.eq.3.or.
     >               cgridopt.eq.LINEAR_GRID.or.
     >               cgridopt.eq.RIBBON_GRID) then
c
c     >               cgridopt.eq.3) then
c slmod end
                   if (ik.gt.nks(ir)/2) then
                     xatiz2(1) = xatiz2(1) + R * sputy
                     yatiz2(1) = yatiz2(1) + Z * sputy
                     ratiz2(1) = ratiz2(1) +  sputy
                   else
                     xatiz2(2) = xatiz2(2) + R * sputy
                     yatiz2(2) = yatiz2(2) + Z * sputy
                     ratiz2(2) = ratiz2(2) +  sputy
                   endif
                 elseif (cgridopt.eq.2) then
                   if ((ir.ge.irsep2.and.ir.le.irwall).or.
     >                (((ir.ge.irtrap2.and.ir.le.nrs2).or.
     >                  (ir.ge.irtrap.and.ir.le.nrs).or.
     >                   ir.lt.irsep).and.
     >                  ik.gt.nks(ir)/2)) then
                     xatiz2(2) = xatiz2(2) + R * sputy
                     yatiz2(2) = yatiz2(2) + Z * sputy
                     ratiz2(2) = ratiz2(2) +  sputy
                   else
                     xatiz2(1) = xatiz2(1) + R * sputy
                     yatiz2(1) = yatiz2(1) + Z * sputy
                     ratiz2(1) = ratiz2(1) +  sputy
                   endif
                 endif
c
                 AATIZ(M) = AATIZ(M) + (ANGLAN+TANLAN) * SPUTY
                 VATIZ(M) = VATIZ(M) + VIN * SPUTY
                 EATIZ(M) = EATIZ(M) + TEMN * SPUTY
                 VATIZM(M)= MAX (VATIZM(M), VIN)
                 RATIZ(M) = RATIZ(M) + SPUTY
                 IFATE = 5
                 GOTO 899
c
c             Exceeded time limit 
c
              elseif (res.eq.2) then 
c
                 fp_neut_tmax = fp_neut_tmax + sputy
                 fp_neut_data(iwstart,3) = fp_neut_data(iwstart,3) 
     >                                   + sputy 

                 RTMAX(M) = RTMAX(M) + SPUTY
                 IFATE = 3
                 GOTO 899
c
c             Hits the wall 
c 
              elseif (res.eq.3) then 
c
                 fp_neut_wall = fp_neut_wall + sputy
                 fp_neut_data(iwstart,4) = fp_neut_data(iwstart,4) 
     >                                     + sputy 
c
                 RWALLN(M) = RWALLN(M) + SPUTY
C
c                Assign a value to ir that corresponds to the
c                index of the wall segment crossed.
c
                 if (fp_ir.eq.irwall-1) then
                    ir = irwall
                 elseif (fp_ir.eq.irtrap+1) then 
                    ir = irtrap
                 else
                    write (6,'(a,5i5,2(1x,f8.3))') 
     >                    'WARNING: FP Ion'//
     >                    ' not on wall boundary ',
     >                   iwstart,ir,ik,fp_ik,fp_ir,r,z
                 endif
c
                 WALLS(IK,IR,0) = WALLS(IK,IR,0) + SPUTY
                 WALLSN(iwstart) = WALLSN(iwstart) + SPUTY
                 wtdep(iwstart,iwstart,2) = 
     >                       wtdep(iwstart,iwstart,2) +sputy
                 IFATE = 1
                 WRITE(6,'(a,4(1x,i6),4(1x,g12.5))') 
     >                            'FP COLLISION WITH WALL  :',
     >                             iwstart,id,ik,ir,R,Z,
     >                      wallpt(iwstart,1),wallpt(iwstart,2)
                 GOTO 899
c
c             Hits target from FP region 
c
              elseif (res.eq.4) then 
c
                 fp_neut_targ = fp_neut_targ + sputy
                 fp_neut_data(iwstart,5) = fp_neut_data(iwstart,5) 
     >                                   + sputy 

c
c                Record NEUTRALS lost to FP target as 
c                being lost to the wall segment they 
c                were launched from for accounting 
c                purposes.  
c
                 RSTRUK(M) = RSTRUK(M) + SPUTY
c
                 WALLSN(iwstart) = WALLSN(iwstart) + SPUTY
                 wtdep(iwstart,iwstart,2) = 
     >                       wtdep(iwstart,iwstart,2) +sputy
c
                 IFATE = 4

                 WRITE(6,'(a,4(1x,i6),4(1x,g12.5))') 
     >            'FP COLLISION WITH TARGET:',
     >                      iwstart,id,IK,IR,R,Z,
     >                      wallpt(iwstart,1),wallpt(iwstart,2)
                 goto 899

              endif

           endif
c
c          Split of initial IF statement for this block of code           
c
           endif  
c
c       If the code falls through to here then the particle continues
c       as a normal neutral launch.    
c
        endif
c 

        IF (DEBUGN) THEN
c          WRITE (6,9003) IPROD,CIST,IK,IR,IX,IY,R,Z,K,
          WRITE (6,9003) IPROD,CIST,IK,IR,0,0,R,Z,K,
     >      VIN,TEMN,SPUTY,(ANGLE+TANGNT)*RADDEG,IT,'NEUTRAL LAUNCH'
        ENDIF
C
C------ CHECK IF NEUTRAL IS GOING TO STRIKE TARGET SURFACE
C
        IF ((cneutb.ne.1.and.cneutb.ne.5.and.
     >       cneutb.ne.6.and.cneutb.ne.7).and.
     >     ((ANGLAN+TANLAN .LE. TANTRU-PI/2.0) .OR.
     >      (ANGLAN+TANLAN .GE. TANTRU+PI/2.0))) THEN
          WRITE(6,*) 'OUT HERE',ANGLAN,TANLAN,
     .               TANTRU-PI/2.0,TANTRU+PI/2.0
          RSTRUK(M) = RSTRUK(M) + SPUTY
c
          if (mtccnt.gt.0) then
             MTC_RSTRUK(M) = MTC_RSTRUK(M) + SPUTY
          endif 
c
c slmod begin
c...      Already an IFATE = 4 above:
          IFATE = 8
c
c          IFATE = 4
c slmod end

! ammod begin.       
	  ! WBC comparison addition for neut struck target.
          Call global_hc_wbc_comp(0,CRMI,VIN,TEMN,SPUTY)
! ammod end.	       
          GOTO 899
        ENDIF
C
        if (iprod.eq.(iprod/100)*100.or.debugn) then 

           write (6,'(a,i6,3i4,5(1x,g13.5))')
     >         'LAUNCH:',iprod,id,ik,ir,r,z,vin,angle,tangnt

        endif 

c        write (0,'(a,i6,3i4,5(1x,g13.5))')
c     >         'LAUNCH:',iprod,id,ik,ir,r,z,vin,angle,tangnt


C
C       SET UP THE WALL DEFINITION IN THE ROUTINE THAT
C       DETERMINES WHETHER A POINT IS INSIDE OR OUTSIDE THE
C       BOUNDARY.
C
        NRFCNT = 0
C
C        WRITE (6,*) 'PCNT:',PCNT,MAXPTS
C        DO 199 J = 1,PCNT
C           WRITE(6,*) RW(J),ZW(J)
C199     CONTINUE
C
        CALL GA15A(PCNT,1,WORK,4*MAXPTS,INDWORK,MAXPTS,
     >             RW,ZW,TDUM,XDUM,YDUM,6)
C
C
C-----------------------------------------------------------------------
C  ITERATE AROUND MOTION FOR EACH TIMESTEP UNTIL IONISATION OCCURS.
C  GENERATE NEW SET OF RANDOMS WHEN OLD LOT ARE ALL USED UP
C-----------------------------------------------------------------------
C
c       Set initial cell occupied. 
c
        iklast = ik  
        irlast = ir
c
c       Check if initial position is inside wall 
c
        CALL GA15B(R,Z,RESULT,PCNT,1,WORK,4*MAXPTS,
     >             INDWORK,MAXPTS,RW,ZW,TDUM,XDUM,YDUM,6)
c
        if (result.lt.0.0) then 
           write(6,'(a,3(1x,g16.8))') 
     >            'SUBROUTINE LAUNCH: WARNING: NEUTRAL INITIAL'//
     >            ' POSITION IS APPARENTLY'//
     >              ' OUTSIDE WALL:',r,z,result
           write (6,'(a,5(1x,i6),5(1x,g13.5))')
     >         'DATA:',iprod,iwstart,id,ik,
     >                ir,r,z,vin,angle,tangnt

c           write (6,'(a,i6,3i4,5(1x,g13.5))')
c     >         'ERROR LAUNCH:',iprod,id,ik,ir,r,z,vmult,cvamult,
c     >                        ktids(id)

c
c          Move particle by one time step - it should only be outside for
c          numerical reasons. 
c
           R    = R + XVELF
           Z    = Z + YVELF
           cist = cist + 1.0

c slmod begin - tmp
          IF (sloutput) THEN
            IF (GRDNMOD.NE.0.AND.STOPOPT.LT.2000) THEN       ! sltmp
c            IF (STOPOPT.LT.MAXNWS) THEN
              WRITE(0,*) 'DEBUG: NEUT ERROR',STOPOPT
              WALKS(STOPOPT,1) = R
              WALKS(STOPOPT,2) = Z
              WALKS(STOPOPT+1,1) = HI
              WALKS(STOPOPT+1,2) = HI
              STOPOPT = STOPOPT + 1
            ENDIF
          ENDIF
c slmod end
c
c
c          Check again if initial position is inside wall 
c
           CALL GA15B(R,Z,RESULT,PCNT,1,WORK,4*MAXPTS,
     >             INDWORK,MAXPTS,RW,ZW,TDUM,XDUM,YDUM,6)
c
           if (result.ge.0.0) then 
               write(6,'(a,6(1x,g16.8))') 
     >            'SUBROUTINE LAUNCH: FOLLOWUP STATUS:'//
     >            ' OK : INSIDE AFTER ONE TIMESTEP:',r,z,result
               err_ok = err_ok + sputy   

           else
               write(6,'(a,6(1x,g16.8))') 
     >            'SUBROUTINE LAUNCH: FOLLOWUP STATUS:'//
     >            ' ERROR: STILL OUTSIDE AFTER ONE TIMESTEP:',r,z,result

               err_out = err_out + sputy   
           endif

c
        endif  
c
c       Write out diagnostic information   
c
        if ( status.le.10.and.
     >      (((nprod-lprod+1).lt.1000).or.
     >       ((nprod-lprod+1).lt.10000.and.(iprod/10)*10.0.eq.iprod).or.
     >       ((iprod/100)*100.0.eq.iprod))   ) then 
           write (6,'(a,2i6,4i4,7(1x,g12.5))') 
     >                   'NEUT-A:',iprod,lprod+iprod-1,
     >                       ik,ir,id,is,r,z,vin,xvelf,yvelf,
     >                       cist
c     >                       ,ZA02AS (1) - STATIM
        endif 
c
c        if ((lprod+iprod-1).eq.52) then 
c           debugn = .true.
c           cstepn = 1
c           tstepn = 1
c           write (6,*) 'DEBUGN Activated:'
c           write (0,*) 'DEBUGN Activated:'
c        endif        
c
c
c
c slmod begin - tmp
c
c       START OF NEUTRAL TRANSPORT LOOP (RELIC GOTO...)
c       ------------------------------------------------------------
c slmod end
  200   CONTINUE

        IF (KK.GT.KKLIM) THEN
          CALL SURAND (SEED, KK, RANV)
          NRAND = NRAND + KK
          KK = 0
        ENDIF
c
c       Adjust impurity velocity for change of cell if option is set. 
c
        if (cneutvel.eq.1.and.
     >      (ik.ne.iklast.or.ir.ne.irlast)) then 
c
c           Record old values
c 
            xvelftmp = xvelf
            yvelftmp = yvelf
            vintmp   = vin
c
c           Assign new velocity
c
            vin = 1.38e4 * sqrt(ktibs(ik,ir)/crmi)
c
c           Assign new temperature 
c
            TEMN   = CRMI * (VIN/1.38E4) * (VIN/1.38E4)
c
c           Reset velocity components
c
            xvelf = xvelftmp * vin/vintmp              
            yvelf = yvelftmp * vin/vintmp              
c
c            if (debugn) write(6,'(a,4(1x,g13.5),2i5,3(1x,g13.5))') 
c     >                    'debugv:',xvelf,yvelf,xvelftmp,yvelftmp,
c     >                    ik,ir,ktibs(ik,ir),vin,vintmp 
c
        endif
c
c       Record current cell
c
        iklast = ik
        irlast = ir 
C
C------UPDATE (X,Y) COORDINATES AND ARRAY INDICES
C
        ROLD = R
        ZOLD = Z
        R    = R + XVELF
        Z    = Z + YVELF
c slmod begin - tmp
        IF (sloutput) THEN
          IF (GRDNMOD.NE.0.AND.STOPOPT.LT.2000) THEN      
            WALKS(STOPOPT,1) = R
            WALKS(STOPOPT,2) = Z
            WALKS(STOPOPT+1,1) = HI
            WALKS(STOPOPT+1,2) = HI
            STOPOPT = STOPOPT + 1
          ENDIF
        ENDIF
c slmod end
! ammod begin.
!       Check if neutral is inside the WBC geometry boundary and count if not.
        if (global_hc_follow_option.gt.0) then 

           hc_wbc_flag = 0
c
           call global_hc_Check_WBC_Neut_Pos(R,Z,CRMI,VIN,
     >                 TEMN,SPUTY,
     >                 NPROD,LPROD,IPROD,real(CIST),hc_wbc_flag)
c 
           if (hc_wbc_flag.eq.1) then 
              IFATE = 7
              GOTO 899
           End If	
        endif 
! ammod end.	       



c
c        if (debugn) write (6,'(a,6(1x,g13.5))') 'DEBUG0:',
c     >                  cist,r,z,xvelf,yvelf,vin
c
        call gridpos(ik,ir,r,z,.false.,griderr)
C
        CALL GA15B(R,Z,RESULT,PCNT,1,WORK,4*MAXPTS,
     >             INDWORK,MAXPTS,RW,ZW,TDUM,XDUM,YDUM,6)
c
        IF (RESULT.GT.0.0) THEN
 
c          WRITE(0 ,*) 'RESULT:',ik,ir,result,griderr
c          WRITE(50,*) 'RESULT:',ik,ir,result,griderr

c          if (debugn) write (6,'(a,7(1x,g13.5),2i5,l5)') 
c     >                        'DEBUG1:',cist,r,z,xvelf,
c     >                        yvelf,vin,result,ik,ir,griderr
c
          CIST = CIST + 1.0
          K    = KKS(IR)
          S    = KSS(IK,IR)
          SMAX = KSMAXS(IR)
C
C
C         ALL OF THE FOLLOWING TESTS ARE USEFUL ONLY WHEN THE
C         PARTICLE IS INSIDE THE VESSEL WALL.
c
c         This has been changed to use the griderr value returned
c         by the gridpos subroutine.
c
C
          if (griderr) then  
c
c            Particle not in grid - record void region density data.
c            Estimate which region the neutral lies in - 
c
c
c            X-point up configuration
c
             if (zxp.gt.z0) then               
c
c               Divertor region.
c
                if (z.ge.zxp) then 
c
c                  Private Plasma void
c
                   if (r.gt.rp(ndsin).and.r.lt.rp(ndsin+1)) then 
c
                      ddvoid(2) = ddvoid(2) + dsputy   
c
c                  Other Divertor void regions. 
c
                   else
c
                      ddvoid(3) = ddvoid(3) + dsputy   
c
                   endif 
c
c               Main plasma void region
c
                elseif (z.lt.zxp) then  

                   ddvoid(1) = ddvoid(1) + dsputy   

                endif
c
c            X-point down configuration
c
             elseif (zxp.le.z0) then  
c
c               Divertor region.
c
                if (z.le.zxp) then 
c
c                  Private Plasma void
c
                   if (r.lt.rp(ndsin).and.r.gt.rp(ndsin+1)) then 
c
                      ddvoid(2) = ddvoid(2) + dsputy   
c
c                  Other Divertor void regions. 
c
                   else
c
                      ddvoid(3) = ddvoid(3) + dsputy   
c
                   endif 
c
c               Main plasma void region
c
                elseif (z.gt.zxp) then  

                   ddvoid(1) = ddvoid(1) + dsputy   

                endif
c              
             endif

c
          elseIF (.not.griderr) THEN
c
c           Particle is found on the grid.
c
 
C
C------ CHECK IF ENTERED MAIN PLASMA, REACHED CENTRE, OR SURVIVED TO
C------ TIME CUTOFF POINT ... OR HIT PLATES
C
            IF     (CFLRIN.AND.IR.LT.IRSEP) THEN
              CFLRIN = .FALSE.
              RMAIN(M) = RMAIN(M) + SPUTY
              ELIMS(IK,3,0) = ELIMS(IK,3,0) + SPUTY
            ELSEIF (CFLREX.AND.IR.GE.IRSEP) THEN
              IF (.NOT.CFLRIN) THEN
                CFLREX = .FALSE.
                REXIT(M) = REXIT(M) + SPUTY
                ELIMS(IK,1,0) = ELIMS(IK,1,0) + SPUTY
              ENDIF
            ENDIF
C
            IF (IR.EQ.1) THEN
              RCENT(M) = RCENT(M) + SPUTY
              IFATE = 2
              GOTO 899
            ELSEIF (CIST.GT.RSTMAX) THEN
              RTMAX(M) = RTMAX(M) + SPUTY
              IFATE = 3
              GOTO 899
            ENDIF
C
C------ SCORE PARTICLE IN DDLIMS ARRAY IN "IONISATION 0" POSITION.
C------ IF TIME POINT REACHED SCORE TIME POSITION ALSO, AND INCREMENT.
C
             DDLIMS(IK,IR,0) = DDLIMS(IK,IR,0) + DSPUTY
             DDTS(IK,IR,0) = DDTS(IK,IR,0) + DSPUTY * TEMN

c
c            At this time also score the velocity shifted contribution 
c            to a specified emission profile by calculating the 
c            component of the particle velocity in the direction
c            of the specified detector location. Get rate 
c            coefficient for each neutral line in use and add to 
c            the array containing the spread wavelength data for
c            this line. This option does not match the absolute 
c            magnitude of the emission - the objective here is to 
c            see what sort of shifts in the wavelength profile of 
c            the emission might be expected. 
c
             if (line_profile_opt.eq.1) then 
c                write(0,*) 'Updating:'  
                call update_line_profile(ik,ir,r,z,
     >                                   xvelf/fsrate,yvelf/fsrate,
     >                                   sputy,cion,rizb)
             endif
c
c
c            Record density separately for any chemically
c            sputtered component
c
             if (neuttype.eq.2.or.neuttype.eq.5) then 
c
                chemden(ik,ir) = chemden(ik,ir) + SPUTY
c
             endif
c
             if (nts.gt.0) then 
c
c               Determine the time bin for the particle 
c 
c                IT = IPOS (CIST, CTIMES(1,0), NTS+1)
                IT = IPOS (real(CIST), CTIMES(1,0), NTS+1)
c
c               IF (CIST.GE.CTIMES(IT,0)) THEN
c
                LIMS(IK,IR,0,IT) = LIMS(IK,IR,0,IT) + SPUTY
c
c               IT = IT + 1
c               ENDIF
c    
              endif
c
C
C------ SET NEW IONISATION PROBABILITY
C------ CHECK FOR IONISATION: IF IT HAS NOT OCCURED JUMP BACK FOR
C------ ANOTHER ITERATION.
C
            IF (DEBUGN) THEN
              IF (CIST.GE.TSTEPN) THEN
  495           TSTEPN = TSTEPN + CSTEPN
                IF (TSTEPN.LE.CIST) GOTO 495
c                WRITE (6,9003) IPROD,CIST,IK,IR,IX,IY,R,Z,K,
                WRITE (6,9003) IPROD,CIST,IK,IR,0,0,R,Z,K,
     >            VIN,TEMN,SPUTY,(ANGLE+TANGNT)*RADDEG,IT
c
c               jdemod - tmp debug
c
c                WRITE (0,9003) IPROD,CIST,IK,IR,0,0,R,Z,K,
c     >            VIN,TEMN,SPUTY,(ANGLE+TANGNT)*RADDEG,IT
              ENDIF
            ENDIF
C
            KK = KK + 1

            IF (RANV(KK).GE.KPCHS(IK,IR,0)) GOTO 200
c
c  State change event has occurred -
C  ------------------------------------------------------------------
C
C           1) ionization
c           2) charge exchange collision (?) - not energetically possible
c           3) momentum transfer collision
c
c
c           Check for and deal with momentum transfer collision first -
c           Since the particle will continue to be tracked as a neutral
c           after the direction of it's velocity is adjusted.
c
c           Note: the same random number is used since it still decides
c                 randomly between these probabilities there is no 
c                 need to draw another - simply check (for now) if
c                 it is greater than the ionization probability since
c                 the only additional state change currently implemented
c                 (11/14/97) is a momentum transfer collision.
c
            if (ranv(kk).gt.kpizs(ik,ir)) then 

c               call execute_mtc(1,mtccnt,mtccist,cist,mtcinf,
c     >               xvelf,yvelf,
c     >               sputy,vin,temn,cneutvel,fsrate,kk,crmi,ik,ir)

               call execute_mtc(1,mtccnt,mtccist,cist,mtcinf,
     >               xvelf,yvelf,
     >               sputy,vin,temn,cneutvel,fsrate,nrand,crmi,
     >               ktibs(ik,ir))
c
c              Momentum Transfer Collision event
c             
c              Record event statistics
c
c               mtccnt = mtccnt + 1 
c
c               mtcinf(1,1) = mtcinf(1,1) + sputy
c
c               if (mtccnt.eq.1) then 
c                  mtcinf(2,1) = mtcinf(2,1) + sputy * (cist-mtccist)
c                  mtcinf(4,1) = mtcinf(4,1) + sputy * 
c     >                               (cist-mtccist)*abs(vin)*fsrate
c               endif 
c
c               mtcinf(3,1) = mtcinf(3,1) + sputy * (cist-mtccist)
c               mtcinf(5,1) = mtcinf(5,1) + sputy *
c     >                            (cist-mtccist)*abs(vin)*fsrate
c               mtcinf(6,1) = mtcinf(6,1) + sputy * temn 
c               mtcinf(7,1) = mtcinf(7,1) + sputy * abs(vin)
c
cc              Update MTC time 
c
c               mtccist = cist
c
c
c              Calculate event (MTC)
c             
c              Choose randomly between +90 and -90 change in velocity
c              - even probability.
c
c               kk = kk + 1
c
c               if (ranv(kk).ge.0.5) then 
c
c                 counter-clockwise 90 degrees
c
c                  tmpvelf = xvelf
c                  xvelf   = -yvelf
c                  yvelf   = tmpvelf
cc
c               else  
c
c                 clockwise 90 degrees
c
c                  tmpvelf = xvelf
c                  xvelf   = yvelf
c                  yvelf   = -tmpvelf
c
c               endif
c
c              If Neutral Velocity Option 2 is in effect - change
c              the neutral velocity in magnitude based on local 
c              conditions.
c
c               if (cneutvel.eq.2) then 
c
c                 Record old values
c 
c                  xvelftmp = xvelf
c                  yvelftmp = yvelf
c                  vintmp   = vin
c
c                 Assign new velocity
c
c                  vin = 1.38e4 * sqrt(ktibs(ik,ir)/crmi)
c
c                 Assign new temperature
c
cc                  TEMN   = CRMI * (VIN/1.38E4) * (VIN/1.38E4)
c
c                 Reset velocity components
c
c                  xvelf = xvelftmp * vin/vintmp              
c                  yvelf = yvelftmp * vin/vintmp              
c
c               endif

c
c              Branch back to neutral following code
c 
               goto 200 
c
            endif

C
C  IONISATION HAS OCCURED : STORE PARTICLE DETAILS IN ARRAYS / TOTALS
C  ------------------------------------------------------------------
C
 250       continue          

c slmod begin - tmp
           IF (sloutput) THEN
             IF (GRDNMOD.NE.0.AND.STOPOPT.LT.2000) THEN      
c               WRITE(0 ,*) 'IK,IR:',ik,ir
c               WRITE(50,*) 'IK,IR:',ik,ir
               WALKS(STOPOPT,1) = R
               WALKS(STOPOPT,2) = Z
               WALKS(STOPOPT+1,1) = HI
               WALKS(STOPOPT+1,2) = HI
               STOPOPT = STOPOPT + 1
             ENDIF
           ENDIF
c slmod end

           TIZS(IK,IR,0) = TIZS(IK,IR,0) + SPUTY
c
           if (neuttype.eq.2.or.neuttype.eq.5) then 
c
              chemizs(ik,ir) = chemizs(ik,ir) + sputy 
c
           endif
c
           NATIZ = NATIZ + 1
           XATIZS(LATIZ+NATIZ-1) = R
           YATIZS(LATIZ+NATIZ-1) = Z
           KATIZS(LATIZ+NATIZ-1) = K
           SATIZS(LATIZ+NATIZ-1) = MIN (S, SMAX-S)
           TEMTIZS(LATIZ+NATIZ-1) = TEMN
           snews(latiz+natiz-1)  = sputy
           cistizs(latiz+natiz-1)  = cist * fsrate / qtim
c
c          Record starting wall or target element for neutral/ion.
c
           IDATIZS(LATIZ+NATIZ-1,1) = idstart
           IDATIZS(LATIZ+NATIZ-1,2) = neuttype
c
           launchdat(LATIZ+NATIZ-1,3) = nrfcnt
           launchdat(LATIZ+NATIZ-1,2) = launchdat(iprod,2)
c
           IF (CNEUTG.EQ.0) THEN
             VINS(LATIZ+NATIZ-1) = 0.0
           ELSEIF (CNEUTG.EQ.1.OR.CNEUTG.EQ.3) THEN
             CALL SURAND2 (SEED, 1, RAN)
             NRAND = NRAND + 1
             VINS(LATIZ+NATIZ-1) = SIGN (0.5*VIN, RAN-0.5)
             IF (CNEUTG.EQ.3) THEN
               CALL SURAND2 (SEED, 1, RAN)
               NRAND = NRAND + 1
               VINS(LATIZ+NATIZ-1)=VINS(LATIZ+NATIZ-1)*2.0*SQRT(RAN)
             ENDIF
           ELSEIF (CNEUTG.EQ.2) THEN
             VINS(LATIZ+NATIZ-1) = SIGN (VIN, 0.5*SMAX-S)
           ENDIF
c
c          Record Ring and target index from which ion originated.
c          for target or wall launches. 
c
           if (neuttype.ne.0) then 

              wtsource(idstart,ir,1,neuttype)=
     >                        wtsource(idstart,ir,1,neuttype)+sputy
              wtsource(idstart,ir,2,neuttype)=
     >                        wtsource(idstart,ir,2,neuttype)
     >                      + satizs(latiz+natiz-1)*sputy
c
c             Accumulate data on the ionization of eroded particles from each wall element
c
              wallse_i(iwstart) =  wallse_i(iwstart) + sputy
c
           endif  
c
c          NOTE:  David Elder,    1995  AUG 3
c
c          There were too many new variables being required 
c          to keep track of the various ionization data under
c          all the varied circumstances - so the information
c          was centralized into one array which contains
c          all of the information - there are some other
c          variables which accumulate the data for 
c          compatibility with the older code - but all of 
c          the new code uses the new setup.
c
c          Here is a summary of the array and its contents
c   
c          IONIZDAT(2    ,2  ,2    ,2      ,5)
c                   isol  m   ifp   irflct  quant
c
c          isol = 1 for SOL information 
c               = 2 for MAIN plasma information
c  
c          m    = 1 for target 1   ik > nks(ir)/2
c               = 2 for target 2   ik =< nks(i2)/2 
c
c          ifp  = 1 for a neutral resulting from a regular
c                   launch    
c               = 2 for a neutral resulting from a Far Periphery
c                   relaunch
c
c          irflct = 1 for a neutral that has NOT been reflected
c                 = 2 for a neutral that has been reflected 
c
c          quant= 1 total weight of neutrals
c               = 2 weighted R coordinate
c               = 3 weighted Z coordinate 
c               = 4 weighted K value
c               = 5 weighted S value
c
c
           if (ir.ge.irsep) then 
              isol = 1
           else
              isol = 2
           endif
c
           if (fpropt.eq.1.and.launchdat(iprod,2).eq.1.0) then 
              ifp = 2
           else 
              ifp = 1
           endif   
c
           if ((nrfopt.eq.1.or.nrfopt.eq.2).and.nrfcnt.gt.0) then 
              irflct = 2
           else 
              irflct = 1
           endif 
c
           ionizdat(isol,m,ifp,irflct,1) = 
     >          ionizdat(isol,m,ifp,irflct,1) +  SPUTY
           ionizdat(isol,m,ifp,irflct,2) = 
     >          ionizdat(isol,m,ifp,irflct,2) +  SPUTY * R
           ionizdat(isol,m,ifp,irflct,3) = 
     >          ionizdat(isol,m,ifp,irflct,3) +  SPUTY * Z
           ionizdat(isol,m,ifp,irflct,4) = 
     >          ionizdat(isol,m,ifp,irflct,4) +  SPUTY * K
           ionizdat(isol,m,ifp,irflct,5) = 
     >          ionizdat(isol,m,ifp,irflct,5) +  SPUTY * min(s,smax-s)
c
           lionizdat(isol,m,ifp,irflct,1) = 
     >          lionizdat(isol,m,ifp,irflct,1) +  SPUTY
           lionizdat(isol,m,ifp,irflct,2) = 
     >          lionizdat(isol,m,ifp,irflct,2) +  SPUTY * R
           lionizdat(isol,m,ifp,irflct,3) = 
     >          lionizdat(isol,m,ifp,irflct,3) +  SPUTY * Z
           lionizdat(isol,m,ifp,irflct,4) = 
     >          lionizdat(isol,m,ifp,irflct,4) +  SPUTY * K
           lionizdat(isol,m,ifp,irflct,5) = 
     >          lionizdat(isol,m,ifp,irflct,5) +  SPUTY * min(s,smax-s)
c
           IF (IR.GE.IRSEP) THEN
             SFRAC (M) = SFRAC (M) + SPUTY
             KSATIZ(M) = KSATIZ(M) + SPUTY * K
             XSATIZ(M) = XSATIZ(M) + SPUTY * R
             YSATIZ(M) = YSATIZ(M) + SPUTY * Z
             SSATIZ(M) = SSATIZ(M) + SPUTY * MIN (S, SMAX-S)
           ELSE
             MFRAC (M) = MFRAC (M) + SPUTY
             KMATIZ(M) = KMATIZ(M) + SPUTY * K
           ENDIF
c
           KATIZ(M) = KATIZ(M) + K * SPUTY
           TATIZ(M) = TATIZ(M) + CIST * FSRATE * SPUTY
           XATIZ(M) = XATIZ(M) + R * SPUTY
           YATIZ(M) = YATIZ(M) + Z * SPUTY
c
c          Record additional time information
c
           cieizs(0) = cieizs(0) + cist * sputy
           citizs(0) = citizs(0) + sputy
c
c slmod begin - ribbon grid
           if (cgridopt.eq.0.or.cgridopt.eq.1.or.cgridopt.eq.3.or.
     >         cgridopt.eq.LINEAR_GRID.or.cgridopt.EQ.RIBBON_GRID) then
c
c           if (cgridopt.eq.0.or.cgridopt.eq.1.or.cgridopt.eq.3) then
c slmod end
             if (ik.gt.nks(ir)/2) then
               xatiz2(1) = xatiz2(1) + R * sputy
               yatiz2(1) = yatiz2(1) + Z * sputy
               ratiz2(1) = ratiz2(1) +  sputy
             else
               xatiz2(2) = xatiz2(2) + R * sputy
               yatiz2(2) = yatiz2(2) + Z * sputy
               ratiz2(2) = ratiz2(2) +  sputy
             endif
           elseif (cgridopt.eq.2) then
             if ((ir.ge.irsep2.and.ir.le.irwall).or.
     >          (((ir.ge.irtrap2.and.ir.le.nrs2).or.
     >            (ir.ge.irtrap.and.ir.le.nrs).or.
     >             ir.lt.irsep).and.
     >            ik.gt.nks(ir)/2)) then
               xatiz2(2) = xatiz2(2) + R * sputy
               yatiz2(2) = yatiz2(2) + Z * sputy
               ratiz2(2) = ratiz2(2) +  sputy
             else
               xatiz2(1) = xatiz2(1) + R * sputy
               yatiz2(1) = yatiz2(1) + Z * sputy
               ratiz2(1) = ratiz2(1) +  sputy
             endif
           endif
c
           AATIZ(M) = AATIZ(M) + (ANGLAN+TANLAN) * SPUTY
           VATIZ(M) = VATIZ(M) + VIN * SPUTY
           EATIZ(M) = EATIZ(M) + TEMN * SPUTY
           VATIZM(M)= MAX (VATIZM(M), VIN)
           RATIZ(M) = RATIZ(M) + SPUTY
           IFATE = 5
           GOTO 899
         ENDIF
C
C
C
        ELSEIF (RESULT.LE.0.0) THEN
c
          if (debugn) write (6,'(a,7(1x,g13.5))') 
     >                       'DEBUG2:',cist,r,z,xvelf,
     >                        yvelf,vin,result
C
C       PARTICLE IS FOUND TO BE OUTSIDE BOUNDARY. NEED TO FIND
C       THE POSITION AND SEGMENT WHERE IT LEFT...
C
c         WRITE(6,'(a,2i4,1p,3g12.5,l4)') 'WALL COLL:',IK,IR,R,Z,
c     >                                 result,griderr
c         write(6,'(a,i4,1p,7g12.5)')  '         :',
c     >          iprod,rold,zold,
c     >          orgr,orgz,xvelf,yvelf,cist
C
c          WRITE(50,*) 'WHOA! PARTICLE OUTSIDE BOUNDARY!'  ! sltmp

c         BEST = HI
c         DSQ  = HI
c         IND = 1
c         DO 650 ID = 1,WALLPTS
c            DSQ = (WALLPT(ID,1)-R) ** 2 + (WALLPT(ID,2)-Z) ** 2
c            IF (DSQ.LT.BEST) THEN
c              BEST = DSQ
c              IND   = ID
c            ENDIF
c650      CONTINUE
C
c         WRITE(6,'(a,2i5,10(1x,g14.8))') 
c     >              'DSQ1:',ind,WALLPTS,DSQ,R,Z,ROLD,ZOLD
c
c         WRITE(6,'(a,i5,10(1x,g14.8))') 
c     >    'DSQ2:',ind,WALLPT(ind,20),wallpt(ind,21),wallpt(ind,1),
c     >                wallpt(ind,2),wallpt(ind,22),wallpt(ind,23)
C
c          CALL INTCALC(R,Z,ROLD,ZOLD,WALLPT(IND,1),WALLPT(IND,2),
c     >                 WALLPT(IND,8),WALLPT(IND,9),WALLPT(IND,5),
c     >                 WALLPT(IND,6),RNEW,ZNEW,TNEW,tnorm,
c     >                 SECT,nrfopt)
c          INDI = IND
C
c          WRITE(6,'(a,l4,8(1x,g13.6))') 'SECT:',SECT,
c     >              R,Z,ROLD,ZOLD,RNEW,ZNEW,orgr,orgz
C
c          IF (.NOT.SECT) THEN
c            DO 660 ID = 1,WALLPTS/2
c              INDI = IND + ID
C
c              WRITE(6,*) 'INDI1:' ,INDI,SECT
C
c              IF (INDI.GT.WALLPTS) INDI = INDI-WALLPTS
c              CALL INTCALC(R,Z,ROLD,ZOLD,WALLPT(INDI,1),WALLPT(INDI,2),
c     >                 WALLPT(INDI,8),WALLPT(INDI,9),WALLPT(INDI,5),
c     >                 WALLPT(INDI,6),RNEW,ZNEW,TNEW,tnorm,
c     >                 SECT,nrfopt)
c              IF (SECT) GOTO 670
c              INDI = IND - ID
C
c              WRITE(6,*) 'INDI2:' ,INDI,SECT
C
c jdemod - note - this is a bug - since INDI <= 0 - the quantity must be 
c                 ADDED to wallpts to get the proper index - not subtracted!
c
c              IF (INDI.LT.1) INDI = WALLPTS-INDI
c
c              IF (INDI.LT.1) INDI = WALLPTS+INDI
c              CALL INTCALC(R,Z,ROLD,ZOLD,WALLPT(INDI,1),WALLPT(INDI,2),
c     >                 WALLPT(INDI,8),WALLPT(INDI,9),WALLPT(INDI,5),
c     >                 WALLPT(INDI,6),RNEW,ZNEW,TNEW,tnorm,
c     >                 SECT,nrfopt)
c              IF (SECT) GOTO 670
c660         CONTINUE
c          ENDIF
c670       CONTINUE
c



c
c         Find location of wall intersection - use of indi for intersection index is 
c         for compatibility with existing code 
c

          call find_wall_intersection(r,z,rold,zold,rnew,znew,tnew,
     >                                tnorm,
     >                                nrfopt,indi,intersect_result,
c slmod begin
     >                                sect,cprint)
c
c     >                                sect)
c slmod end
c
c         Verify RNEW,ZNEW 
c
c          if ( (abs(r-rnew)+abs(rold-rnew)).ne.abs(r-rold).or.
c     >         (abs(z-znew)+abs(zold-znew)).ne.abs(z-zold)) then 
c
c             write (0,'(a,i5,l4,6(1x,g13.6))') 
c     >           'NEUT: WARNING: POSSIBLE RNEW,ZNEW ERROR:',indi,sect,
c     >                      r,z,rnew,
c     >                      znew,rold,zold
c             write (6,'(a,i5,l4,6(1x,g13.6))') 
c     >           'NEUT: WARNING: POSSIBLE RNEW,ZNEW ERROR:',indi,sect,
c     >                      r,z,rnew,
c     >                      znew,rold,zold
c
c             CALL GA15B(R,Z,RESULT,PCNT,1,WORK,4*MAXPTS,
c     >             INDWORK,MAXPTS,RW,ZW,TDUM,XDUM,YDUM,6)
c             write (6,'(a,6(1x,g13.6))') 
c     >           'DATA: R,Z,RESULT:',r,z,result
c
c             CALL GA15B(Rold,ZOLD,RESULT,PCNT,1,WORK,4*MAXPTS,
c     >             INDWORK,MAXPTS,RW,ZW,TDUM,XDUM,YDUM,6)
c             write (6,'(a,6(1x,g13.6))') 
c     >           'DATA: ROLD,ZOLD,RESULT:',r,z,result
c
c          endif
c



c
c         Dealing with possible reflections
c
c
c         Check for reflection - allowing for varying probabilities
c         for each segment - IF nrfopt is ON.
c
          reflect_neut = .false. 
c
c         Check if reflection is allowed
c
          if (nrfopt.ne.0) then 
c
c            Check for valid wall segment with non-zero reflection
c            probability
c
             if (indi.ge.1.and.indi.le.wallpts.and.
     >           abs(wallpt(indi,25)).gt.0.0) then 
c
                 refprob = min(abs(wallpt(indi,25)),1.0)
c         
                 NRAND = NRAND + 1
                 CALL SURAND2 (SEED, 1, RAN)
c
c                Draw random number and check for reflection
c
                 if (ran.le.refprob) reflect_neut = .true.
c
             endif                
c
          endif
c
c         Execute reflection or deposition
c
          IF (.not.reflect_neut) THEN
c
            IF (.NOT.SECT) THEN
              WRITE(6,'(a,2g12.5)') 'NEUT: ERROR: NEUTRAL'//
     >                              ' AT WALL BUT NO INTERSECTION'//
     >                              ' POINT FOUND',R,Z
              WRITE(0,'(a,2g12.5)') 'NEUT: ERROR: NEUTRAL'//
     >                              ' AT WALL BUT NO INTERSECTION'//
     >                              ' POINT FOUND',R,Z
            ENDIF
c
c           Target collision
c
c slmod begin
c
c     jdemod - using wlwall1, wlwall2 ... for generalized grids (or ribbon grid)
c              won't work ... so switch to using wallpt(indi,18) which is a pointer
c              to the matching target element or 0.0 otherwise
c
c              Looks like Steve partially implemented this for some cases
c
            if (wallpt(indi,18).ne.0.0) then 
c
c            IF ((GRDNMOD.NE.0.AND.WALLPT(INDI,18).NE.0.0).OR.
c     >          (GRDNMOD.EQ.0.AND.
c     >           ((INDI.GT.WLWALL2.AND.INDI.LT.WLTRAP1).OR.
c     >            (INDI.GT.WLTRAP2.AND.INDI.LE.WALLPTS)))) THEN
c
c            IF ((INDI.GT.WLWALL2.AND.INDI.LT.WLTRAP1).OR.
c     >           (INDI.GT.WLTRAP2.AND.INDI.LE.WALLPTS)) THEN
c slmod end
               RSTRUK(M) = RSTRUK(M) + SPUTY
c
               if (mtccnt.gt.0) then
                  MtC_RSTRUK(M) = MTC_RSTRUK(M) + SPUTY
               endif 
c
c              Record particles with invalid ID's in total 
c
               if (indi.lt.1.or.indi.gt.wallpts) then 
                  WALLSN(maxpts+1) = WALLSN(maxpts+1) + SPUTY
c
                  if (iwstart.ge.1.and.iwstart.le.wallpts) then  
                     wtdep(iwstart,maxpts+1,2) = 
     >                       wtdep(iwstart,maxpts+1,2) +sputy
                  endif
c
               else
                  WALLSN(INDI) = WALLSN(INDI) + SPUTY
c
                  if (iwstart.ge.1.and.iwstart.le.wallpts) then  
                     wtdep(iwstart,indi,2) = 
     >                       wtdep(iwstart,indi,2) +sputy
                  endif
c
               endif 

! ammod begin.
		  ! Add to deposition.
		  NEROS (INT (wallpt (INDI,18)),1) = 
     >              NEROS (INT (wallpt (INDI,18)),1) + SPUTY
! ammod end.

c
               IFATE = 4


! ammod begin.       
	       ! WBC comparison addition for neut struck target.
               ! Call Record_WBC_Ion_Event (0,CRMI,VIN,TEMN,SPUTY)
               Call global_hc_wbc_comp(0,CRMI,VIN,TEMN,SPUTY)
! ammod end.	       



               WRITE(6,'(a,5(1x,i6),4(1x,g12.5))') 
     >            'COLLISION WITH TARGET:',INDI,
     >                      iwstart,id,IK,IR,R,Z,
     >                      wallpt(indi,1),wallpt(indi,2)
               GOTO 899
c
c           Wall collision
c
            ELSE
c
               RWALLN(M) = RWALLN(M) + SPUTY
C
               IF (MTCCNT.GT.0) THEN 
                  MTC_RWALLN(M) = MTC_RWALLN(M) + SPUTY
               ENDIF    
C
c              Assign a value to ir that corresponds to the
c              index of the wall segment crossed.
c
c slmod begin
c
c              jdemod - this is more difficult for generalized grids that
c                       may have multiple nearest rings. 
c                       walls(ik,ir,iz) is outmoded now and may need replacing
c                       for more modern grids. 
c
c              This is known to be a wall collision so assign irwall unless 
c              it is in a range for irtrap. Note that wltrap1, wltrap2 are NOT
c              always well defined for all grids ... however, when they are the
c              following code should suffice.  
c
c               if (grdnmod.ne.0) then
c
                 if (indi.ge.wltrap1.and.indi.le.wltrap2) then 
c...                WLTRAP1,2 are still well defined for the generalized geometry, 
c                   and probably always will be (same login in NEUTONE.F):
                    ir = irtrap
                 else
c...                Everything else must be IRWALL (or not..?):
                    ir = irwall
                 endif
c
c
c               else
c                 if (indi.ge.wlwall1.and.indi.le.wlwall2) then
c                    ir = irwall
c                 elseif (indi.ge.wltrap1.and.indi.le.wltrap2) then 
c                    ir = irtrap
c                 else
c                    write (6,*) 'Neutral not on wall segment'//
c     >                      ' at collision',
c     >                     indi,ir,ik,r,z
c                 endif
c               endif
c
c               if (indi.ge.wlwall1.and.indi.le.wlwall2) then
c                  ir = irwall
c               elseif (indi.ge.wltrap1.and.indi.le.wltrap2) then 
c                  ir = irtrap
c               else
c                  write (6,*) 'Neutral not on wall segment'//
c     >                    ' at collision',
c     >                   indi,ir,ik,r,z
c               endif
c slmod end
c
               WALLS(IK,IR,0) = WALLS(IK,IR,0) + SPUTY
c
c
c              Record particles with invalid ID's in total 
c
               if (indi.lt.1.or.indi.gt.wallpts) then 
                  WALLSN(maxpts+1) = WALLSN(maxpts+1) + SPUTY
c
                  if (iwstart.ge.1.and.iwstart.le.wallpts) then  
                     wtdep(iwstart,maxpts+1,2) = 
     >                       wtdep(iwstart,maxpts+1,2) +sputy
                  endif
c
               else
                  WALLSN(INDI) = WALLSN(INDI) + SPUTY
c
                  if (iwstart.ge.1.and.iwstart.le.wallpts) then  
                     wtdep(iwstart,indi,2) = 
     >                       wtdep(iwstart,indi,2) +sputy
                  endif
c
               endif 
c
               IFATE = 1

! ammod begin.       
	       ! WBC comparison addition for neut struck vessel wall.
               ! Call Record_WBC_Ion_Event (0,CRMI,VIN,TEMN,SPUTY)
               Call global_hc_wbc_comp(0,CRMI,VIN,TEMN,SPUTY)
! ammod end.	       


               WRITE(6,'(a,5(1x,i6),4(1x,g12.5))') 
     >                            ' COLLISION WITH WALL  :',
     >                             INDI,iwstart,id,ik,ir,R,Z
               GOTO 899
            ENDIF
          ELSEIF (reflect_neut) THEN
            IF (.NOT.SECT) THEN
c
              WRITE(6,'(a,2g12.5)') 'NEUT: ERROR: NEUTRAL'//
     >                              ' AT WALL BUT NO INTERSECTION'//
     >                              ' POINT FOUND',R,Z
              WRITE(0,'(a,2g12.5)') 'NEUT: ERROR: NEUTRAL'//
     >                              ' AT WALL BUT NO INTERSECTION'//
     >                              ' POINT FOUND',R,Z
c
c             jdemod - change code for target impact detection
c
              if (wallpt(indi,18).ne.0.0) then
c
c              IF ((INDI.GT.WLWALL2.AND.INDI.LT.WLTRAP1).OR.a
c     >            (INDI.GT.WLTRAP2.AND.INDI.LE.WALLPTS)) THEN
c
                 RSTRUK(M) = RSTRUK(M) + SPUTY
c
                 if (mtccnt.gt.0) then
                    MTC_RSTRUK(M) = MTC_RSTRUK(M) + SPUTY
                 endif 
c
c                Record particles with invalid ID's in total 
c
                 if (indi.lt.1.or.indi.gt.wallpts) then 
                    WALLSN(maxpts+1) = WALLSN(maxpts+1) + SPUTY
c
                    if (iwstart.ge.1.and.iwstart.le.wallpts) then  
                       wtdep(iwstart,maxpts+1,2) = 
     >                       wtdep(iwstart,maxpts+1,2) +sputy
                    endif
c
                  else
                    WALLSN(INDI) = WALLSN(INDI) + SPUTY
c
                    if (iwstart.ge.1.and.iwstart.le.wallpts) then  
                       wtdep(iwstart,indi,2) = 
     >                       wtdep(iwstart,indi,2) +sputy
                    endif
c
                 endif 
c
! ammod begin.
		  ! Add to deposition.
                  if (wallpt(indi,18).ne.0) then 
                     NEROS (INT (wallpt (INDI,18)),1) = 
     >                  NEROS (INT (wallpt (INDI,18)),1) + SPUTY
                  endif
! ammod end.

                 IFATE = 4

! ammod begin.       
	         ! WBC comparison addition for neut struck target.
                 ! Call Record_WBC_Ion_Event (0,CRMI,VIN,TEMN,SPUTY)
                 Call global_hc_wbc_comp(0,CRMI,VIN,TEMN,SPUTY)
! ammod end.	       


                 GOTO 899
c
c             Impact is on a wall element
c
              ELSE
c
                 RWALLN(M) = RWALLN(M) + SPUTY
C
                 IF (MTCCNT.GT.0) THEN 
                    MTC_RWALLN(M) = MTC_RWALLN(M) + SPUTY
                 ENDIF    
C
                 WALLS(IK,IR,0) = WALLS(IK,IR,0) + SPUTY
c
c                Record particles with invalid ID's in total 
c
                 if (indi.lt.1.or.indi.gt.wallpts) then 
                    WALLSN(maxpts+1) = WALLSN(maxpts+1) + SPUTY
c 
                    if (iwstart.ge.1.and.iwstart.le.wallpts) then  
                       wtdep(iwstart,maxpts+1,2) = 
     >                       wtdep(iwstart,maxpts+1,2) +sputy
                    endif
c
                  else
                    WALLSN(INDI) = WALLSN(INDI) + SPUTY
c
                    if (iwstart.ge.1.and.iwstart.le.wallpts) then  
                       wtdep(iwstart,indi,2) = 
     >                       wtdep(iwstart,indi,2) +sputy
                    endif
c
                 endif 
c
                 IFATE = 1


! ammod begin.       
	         ! WBC comparison addition for neut struck vessel wall.
                 ! Call Record_WBC_Ion_Event (0,CRMI,VIN,TEMN,SPUTY)
                 Call global_hc_wbc_comp(0,CRMI,VIN,TEMN,SPUTY)
! ammod end.	       


                 GOTO 899
              ENDIF
c
c           Reflection from intersecting segment. 
c
            ELSE
c
c             Set various quantities related to reflection coefficients
c
c              if (wallpt(indi,25).eq.0.0) then 
c
c                Issue warning but reflect the particle anyway -
c                the code should NEVER execute this ...       
c
c                 write (6,*) 'WARNING: REFLECTION FROM SEGMENT'//
c     >                    ' WITH ZERO REFLECTION COEFFICIENT',indi
c
c             The implementation has changed - the probability
c             of reflection is used to determine if a particle is
c             reflected or not - it's weight does not change. This 
c             is an alternative implementation that may be used
c             later if the statictics are insufficient.        
c
c              elseif (wallpt(indi,25).gt.0.0) then 
c
c                Modify particle weight but leave all else unchanged
c 
c
c                 sputy = sputy * wallpt(indi,25)
c
c              elseif (wallpt(indi,25).lt.0.0) then 
c
              if (wallpt(indi,25).lt.0.0) then 
c
c                 Modify particle energy to
c                 specified thermal emission energy.
c
c                  sputy = sputy * abs(wallpt(indi,25))
c
                  VIN = 1.38E4 * SQRT (CTEM1/CRMI)
c  
c                 Assign new temperature as well 
c  
                  TEMN   = CRMI * (VIN/1.38E4) * (VIN/1.38E4)
c
              endif             
c 
              if (debugn) write (6,*) 'v1:',xvelf,yvelf,vin,tnew
c
              XVELF  = VIN * COS(TNEW) * FSRATE
              YVELF  = VIN * SIN(TNEW) * FSRATE
c
              if (debugn) write(6,*) 'v2:',xvelf,yvelf
c
              move_factor = 1.0
c
 5000         Rtmp = RNEW + move_factor * XVELF
              Ztmp = ZNEW + move_factor * YVELF
c
              CALL GA15B(Rtmp,Ztmp,RESULT,PCNT,1,WORK,4*MAXPTS,
     >             INDWORK,MAXPTS,RW,ZW,TDUM,XDUM,YDUM,6)
c
              if (result.lt.0) then 
c
                 move_factor = move_factor * 0.1

                 if (move_factor.lt.1.0e-8) then 

                    CALL GA15B(ROLD,ZOLD,RESULT,PCNT,1,WORK,4*MAXPTS,
     >                  INDWORK,MAXPTS,RW,ZW,TDUM,XDUM,YDUM,6)

                    write (6,'(a,8(1x,g13.6))') 
     >                 'WALL COLL: NEW R,Z OUTSIDE:',
     >                       RNEW,ZNEW,R,Z,XVELF,YVELF,TNEW
                    write (6,'(a,3i5,5(1x,g13.6))') 'MORE1:',ik,ir,
     >                           indi,rold,zold,result,
     >                            wallpt(indi,16),vin
                    write (6,'(a,8(1x,g13.6))') 'MORE2:',
     >                           wallpt(indi,1),wallpt(indi,2),
     >                           wallpt(indi,20),wallpt(indi,21),
     >                           wallpt(indi,22),wallpt(indi,23),
     >                           wallpt(indi,8),wallpt(indi,9)
                    rtmp = rnew + xvelf
                    ztmp = znew + yvelf

                    call gridpos(ik,ir,rold,zold,.false.,griderr)

                    CALL GA15B(Rtmp,Ztmp,RESULT,PCNT,1,WORK,4*MAXPTS,
     >                  INDWORK,MAXPTS,RW,ZW,TDUM,XDUM,YDUM,6)

c                    write (6,'(a,2i5,8(1x,g13.6),l4)') 'MORE3:',ik,ir,
c     >                           rtmp,ztmp,result,griderr

                    call gridpos(ik,ir,r,z,.false.,griderr)

                    CALL GA15B(R,Z,RESULT,PCNT,1,WORK,4*MAXPTS,
     >                  INDWORK,MAXPTS,RW,ZW,TDUM,XDUM,YDUM,6)

c                    write (6,'(a,2i5,8(1x,g13.6),l4)') 'MORE3:',ik,ir,
c     >                           r,z,result,griderr
c
                    RWALLN(M) = RWALLN(M) + SPUTY
C
                    IF (MTCCNT.GT.0) THEN 
                       MTC_RWALLN(M) = MTC_RWALLN(M) + SPUTY
                    ENDIF    
C
                    if (ir.ne.irtrap.and.ir.ne.irwall.and.
     >                  (ik.ne.1.or.ik.ne.nks(ir))) then
                       write (6,*)
     >                   'Neutral not in adjacent cell at collision',
     >                              ik,ir
c
c                       if (ir.gt.irtrap) then
c                          ir = irtrap
c                       else
c                          ir = irwall
c                       endif
c
                    endif

                    WALLS(IK,IR,0) = WALLS(IK,IR,0) + SPUTY
c
c                   Record particles with invalid ID's in total 
c
                    if (indi.lt.1.or.indi.gt.wallpts) then 
                       WALLSN(maxpts+1) = WALLSN(maxpts+1) + SPUTY
c
                       if (iwstart.ge.1.and.iwstart.le.wallpts) then  
                          wtdep(iwstart,maxpts+1,2) = 
     >                       wtdep(iwstart,maxpts+1,2) +sputy
                       endif
c
                    else
                       WALLSN(INDI) = WALLSN(INDI) + SPUTY
c
                       if (iwstart.ge.1.and.iwstart.le.wallpts) then  
                          wtdep(iwstart,indi,2) = 
     >                       wtdep(iwstart,indi,2) +sputy
                       endif
c
                    endif 
c
! ammod begin.       
	            ! WBC comparison addition for neut struck vessel wall.
                    ! Call Record_WBC_Ion_Event (0,CRMI,VIN,TEMN,SPUTY)
                    Call global_hc_wbc_comp(0,CRMI,VIN,TEMN,SPUTY)
! ammod end.	       


                    IFATE = 1
                    GOTO 899
                 endif     
c
                 goto 5000 
c
              endif
c
c              CIST = CIST + 1.0
c
              r = rtmp 
              z = ztmp
c slmod begin - tmp
              IF (sloutput) THEN
                IF (GRDNMOD.NE.0.AND.STOPOPT.LT.2000) THEN      
                  WALKS(STOPOPT,1) = R
                  WALKS(STOPOPT,2) = Z
                  WALKS(STOPOPT+1,1) = HI
                  WALKS(STOPOPT+1,2) = HI
                  STOPOPT = STOPOPT + 1
                ENDIF
              ENDIF
c slmod end

c
              NRFCNT = NRFCNT + 1
              TOTRF  = TOTRF +1
              MAXRF = MAX0(NRFCNT,MAXRF)
C
              WRITE(6,'(a,i5,1p,10g11.4)') 'NRF:',indi,
     >                       RNEW,ZNEW,ROLD,ZOLD,
     >                       tnew,R,Z,XVELF,YVELF,cist
C
c             Number of reflections as a neutral exceeds the maximum allowed
c
              IF (NRFCNT.GT.maxnrfcnt) THEN
c
                nrfloss = nrfloss + sputy 
c
c               Neutral at target segment
c
c               jdemod - change target detection
c
                if (wallpt(indi,18).ne.0.0) then 
c
c                IF ((INDI.GT.WLWALL2.AND.INDI.LT.WLTRAP1).OR.
c     >              (INDI.GT.WLTRAP2.AND.INDI.LE.WALLPTS)) THEN

                   write(0,'(a)') 'NEUT: ERROR: NEUTRAL AT WALL'//
     >                        ' REMOVED DUE TO TOO MANY REFLECTIONS'//
     >                        ': PARTICLE LIKELY TRAPPED'
                   write(6,'(a)') 'NEUT: ERROR: NEUTRAL AT WALL'//
     >                        ' REMOVED DUE TO TOO MANY REFLECTIONS'//
     >                        ': PARTICLE LIKELY TRAPPED'
                   write(6,'(a,3i6,4g12.5)') 'DATA:',ik,ir,indi,rnew,
     >                                     znew,xvelf,yvelf    

                   write(6,'(a,5i6)') 'DATA:',wlwall1,wlwall2,wltrap1,
     >                                       wltrap2,
     >                          wallpts
                   write(6,*) 'DATA:',
     >                         (id,':',wallpt(indi,id),',',id=1,13)
c
                   RSTRUK(M) = RSTRUK(M) + SPUTY
c
                   if (mtccnt.gt.0) then
                      MTC_RSTRUK(M) = MTC_RSTRUK(M) + SPUTY
                   endif 
c
c
! ammod begin.
		   ! Add to deposition.
                   if (wallpt(indi,18).ne.0) then
                      NEROS (INT (wallpt (INDI,18)),1) = 
     >                      NEROS (INT (wallpt (INDI,18)),1) + SPUTY
                   endif

! ammod end.

                   IFATE = 4


! ammod begin.       
	           ! WBC comparison addition for neut struck target.
                   ! Call Record_WBC_Ion_Event (0,CRMI,VIN,TEMN,SPUTY)
                   Call global_hc_wbc_comp(0,CRMI,VIN,TEMN,SPUTY)
! ammod end.	       



                   GOTO 899
c
c              Neutral on wall segment 
c 
               ELSE
                   write(0,'(a)') 'NEUT: ERROR: NEUTRAL AT WALL'//
     >                        ' REMOVED DUE TO TOO MANY REFLECTIONS'//
     >                        ': PARTICLE LIKELY TRAPPED'
                   write(6,'(a)') 'NEUT: ERROR: NEUTRAL AT WALL'//
     >                        ' REMOVED DUE TO TOO MANY REFLECTIONS'//
     >                        ': PARTICLE LIKELY TRAPPED'
                   write(6,'(a,3i6,4g12.5)') 'DATA:',ik,ir,indi,rnew,
     >                                     znew,xvelf,yvelf    

                   write(6,'(a,5i6)') 'DATA:',wlwall1,wlwall2,wltrap1,
     >                                       wltrap2,
     >                          wallpts
                   write(6,*) 'DATA:',(id,':',
     >                             wallpt(indi,id),',',id=1,13)
c
c
c
                   RWALLN(M) = RWALLN(M) + SPUTY
C
                   IF (MTCCNT.GT.0) THEN 
                      MTC_RWALLN(M) = MTC_RWALLN(M) + SPUTY
                   ENDIF    
C
                   if (ir.ne.irtrap.and.ir.ne.irwall) then
c
                      write (6,'(a,i6)')
     >                  'NEUT: WARNING: NEUTRAL NOT ASSOCIATED'//
     >                        ' WITH WALL RING AT WALL COLLISION',ir
c
                      if (ir.gt.irtrap) then
                         ir = irtrap
                      else
                         ir = irwall
                      endif
                   endif
c
                   WALLS(IK,IR,0) = WALLS(IK,IR,0) + SPUTY
c
c                  Record particles with invalid ID's in total 
c
                   if (indi.lt.1.or.indi.gt.wallpts) then 
                      WALLSN(maxpts+1) = WALLSN(maxpts+1) + SPUTY
c
                      if (iwstart.ge.1.and.iwstart.le.wallpts) then  
                         wtdep(iwstart,maxpts+1,2) = 
     >                       wtdep(iwstart,maxpts+1,2) +sputy
                      endif
c
                   else
                      WALLSN(INDI) = WALLSN(INDI) + SPUTY
c
                      if (iwstart.ge.1.and.iwstart.le.wallpts) then  
                         wtdep(iwstart,indi,2) = 
     >                       wtdep(iwstart,indi,2) +sputy
                      endif
c
                   endif 
c
                   IFATE = 1

! ammod begin.       
	           ! WBC comparison addition for neut struck target.
                   ! Call Record_WBC_Ion_Event (0,CRMI,VIN,TEMN,SPUTY)
                   Call global_hc_wbc_comp(0,CRMI,VIN,TEMN,SPUTY)
! ammod end.	       



                   GOTO 899
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
C
       GOTO 200
c slmod begin - tmp
c      ------------------------------------------------------------
c      END OF NEUTRAL TRANSPORT LOOP 
c
c slmod end
C
  899  CONTINUE

c
c     jdemod - temp debug
c
c       WRITE(0 ,*) 'IFATE:',ifate   ! sltmp
c       WRITE(0 ,*) '     :',r,z,ik,ir
       IF (ifate.EQ.8) target_loss = target_loss + 1
       WRITE(50,*) 'IFATE:',ifate   ! sltmp
       WRITE(50,*) '     :',r,z,ik,ir
c       call gridpos(ik,ir,r,z,.false.,griderr)
c       WRITE(0 ,*) '     :',r,z,ik,ir,griderr
c       WRITE(50,*) '     :',r,z,ik,ir,griderr
c       griderr = incell(ik,ir,r,z)
c       WRITE(0 ,*) '     :',r,z,ik,ir,griderr
c       WRITE(50,*) '     :',r,z,ik,ir,griderr
       WRITE(6,*) 'IFATE:',ifate   ! sltmp
c
c      Record MTC event frequency 
c
       if (mtccnt.gt.10) then 
          mtctotcnt(11,1) = mtctotcnt(11,1) + sputy
       else  
          mtctotcnt(mtccnt,1) = mtctotcnt(mtccnt,1) + sputy
       endif
c
       IF (DEBUGN.OR.(IFATE.EQ.6))
c     >   WRITE (6,9003) IPROD,CIST,IK,IR,IX,IY,R,Z,K,
     >   WRITE (6,9003) IPROD,CIST,IK,IR,0,0,R,Z,K,
     >     VIN,TEMN,SPUTY,(ANGLE+TANGNT)*RADDEG,IT,FATE(IFATE)
c
       CISTOT = CISTOT + CIST * SPUTY
c       CISMAX = MAX (CISMAX, CIST)
       CISMAX = MAX (CISMAX, real(CIST))
c
c      Main loop end 
c
  900  CONTINUE
C
C-----------------------------------------------------------------------
C      PRINT DIAGNOSTICS
C-----------------------------------------------------------------------
C
c slmod begin
      IF (sloutput) THEN
        WRITE(0 ,*) 'TARGET_LOSS:',target_loss
        WRITE(6 ,*) 'TARGET_LOSS:',target_loss
        WRITE(50,*) 'TARGET_LOSS:',target_loss
      ENDIF
c slmod end
      YMF = KMFPS(1)
c
      IF (STATUS.NE.1) YMF = KMFSS(1)
      IF (MATT.EQ.13.OR.MATT.EQ.14.OR.MATT.EQ.15
     >    .OR.MATT.EQ.16.OR.MATT.EQ.17.or.matt.eq.18
     >    .or.matt.eq.19) THEN
c
c      switch inner to outer for AUG; Krieger, 6/95
c
c        if (zxp.le.z0) then 
c           WRITE (7,9012) YMF                                              
c        else
c
           WRITE (7,9002) YMF,inner,outer
c
c        endif 
c
      ELSE
c
c      switch inner to outer for AUG; Krieger, 6/95
c
c        if (zxp.le.z0) then 
c           WRITE (7,9011) NINT(CETH(MATP,MATT)),
c     >            NINT(CETF(MATP,MATT)),CQ(MATP,MATT),YMF
c        else
c
           WRITE (7,9001) NINT(CETH(MATP,MATT)),
     >            NINT(CETF(MATP,MATT)),CQ(MATP,MATT),YMF,inner,outer
c
c        endif 
c
      ENDIF
C
      IF (CNEUTB.EQ.2.or.cneutb.eq.4) THEN
         CALL PRC('NEUTRALS LAUNCHED FROM WALLS:')
      ENDIF
C
      CALL PRI2 ('NUMBER OF NEUTRALS LAUNCHED        ',
     >                             NINT(RNEUT (1)),NINT(RNEUT (2)))
c
      tmpsum = 0.0
      do in = 1,wallpts
         tmpsum = tmpsum + wallse(in)
      end do
c
c      call prr  ('TEMP: TOTAL EROSION: ',tmpsum)
c
      IF (STATUS.LE.10) THEN
      CALL PRI2 ('NUMBER OF NEUTRALS STOPPED AT WALL ',
     >                             NINT(RWALLN(1)),NINT(RWALLN(2)))
c
      if (mtcopt.eq.1) then 
         CALL PRI2 ('NUMBER OF THESE NEUTRALS HAVING MTC',
     >                   NINT(MTC_RWALLN(1)),NINT(MTC_RWALLN(2)))
      endif 
c
      IF (NRFOPT.EQ.1.or.nrfopt.eq.2) THEN
      call prb
      call prc  ('REFLECTED NEUTRAL SUMMARY:')
      CALL PRI  ('TOTAL NUMBER OF WALL REFLECTIONS   ',TOTRF)
      CALL PRI  ('MAXIMUM REFLECTIONS/PARTICLE       ',MAXRF)
      call pri  ('MAXIMUM ALLOWED REFLECTIONS        ',maxnrfcnt)
      call prr  ('NUMBER OF NEUTRALS EXCEEDING MAX   ',nrfloss)
      call prb 
      ENDIF
c
c      CALL PRI2 ('NO OF NEUTS REACHING CENTRAL MIRROR',
c     >                             NINT(RCENT (1)),NINT(RCENT (2)))
c      CALL PRI2 ('NUMBER OF NEUTRALS EXISTING AT 0.1S',
c     >                             NINT(RTMAX (1)),NINT(RTMAX (2)))
c      CALL PRI2 ('NUMBER OF NEUTRALS STRIKING TARGET ',
c     >                             NINT(RSTRUK(1)),NINT(RSTRUK(2)))
c
      CALL PRR2 ('NO OF NEUTS REACHING CENTRAL MIRROR',
     >                             RCENT (1),RCENT (2))
      CALL PRR2 ('NUMBER OF NEUTRALS EXISTING AT 0.1S',
     >                             RTMAX (1),RTMAX (2))
      CALL PRR2 ('NUMBER OF NEUTRALS STOPPED AT TARGET',
     >                             RSTRUK(1),RSTRUK(2))
c
      if (mtcopt.eq.1) then  
         CALL PRR2 ('NUMBER OF THESE NEUTRALS HAVING MTC',
     >                             MTC_RSTRUK(1),MTC_RSTRUK(2))
      endif
c
c slmod begin
      CALL PRI2 ('NO FAILED LAUNCHES (100000 DISCARDS) ',
c
c      CALL PRI2 ('NO FAILED LAUNCHES (1000 DISCARDS) ',
c slmod end
     >                             NINT(RFAIL (1)),NINT(RFAIL (2)))
      CALL PRI2 ('NO OF VELOCITIES > VMAX (DISCARDED)',
     >                             NINT(REJECT(1)),NINT(REJECT(2)))
c
      CALL PRR2 ('AVERAGE R PRODUCTION POSITION (M)  ',
     >                            XTOT(1)/RNEUT(1),XTOT(2)/RNEUT(2))
      CALL PRR2 ('AVERAGE Z PRODUCTION POSITION (M)  ',
     >                            YTOT(1)/RNEUT(1),YTOT(2)/RNEUT(2))
      CALL PRR2 ('AVERAGE ANGLE AT PRODUCTION (DEGS) ',
     >                RADDEG*ATOT(1)/RNEUT(1),RADDEG*ATOT(2)/RNEUT(2))
      CALL PRR2 ('AVERAGE VELOCITY WITHOUT ANGLE FACT',
     >                           VTOTA(1)/RNEUT(1),VTOTA(2)/RNEUT(2))
      CALL PRR2 ('MAXIMUM VELOCITY WITHOUT ANGLE FACT',
     >                                   VTOTAM(1),VTOTAM(2))
      CALL PRR2 ('AVERAGE ANGULAR VEL CORRECTION FACT',
     >                           VMULTT(1)/RNEUT(1),VMULTT(2)/RNEUT(2))
      CALL PRR2 ('MAXIMUM ANGULAR VEL CORRECTION FACT',
     >                                   VMULTM(1),VMULTM(2))
      CALL PRR2 ('AVERAGE VELOCITY AT PRODUCTION(M/S)',
     >                           VTOT(1)/RNEUT(1),VTOT(2)/RNEUT(2))
      CALL PRR2 ('MAXIMUM VELOCITY AT PRODUCTION(M/S)',
     >                                   VTOTM(1),VTOTM(2))
      CALL PRR2 ('AVERAGE TEMP AT PRODUCTION (EV)    ',
     >                           ETOT(1)/RNEUT(1),ETOT(2)/RNEUT(2))
      CALL PRI2 ('NO OF NEUTRALS ENTERING MAIN PLASMA',
     >                           NINT(RMAIN(1)),NINT(RMAIN(2)))
      CALL PRI2 ('NO OF THESE NEUTS EXITING MAIN P.  ',
     >                           NINT(REXIT(1)),NINT(REXIT(2)))
      ENDIF
c
c     Print Far Periphery Neutral ionization statistics if the option is active 
c
      if (fp_neut_opt.gt.0) then 
c
         call prb
         call prc('FP NEUTRAL IONIZATION STATISTICS:')
         call prb
         call prr(' - TOTAL NUMBER OF NEUTRALS LAUNCHED: ',
     >             rneut(1)+rneut(2))
         call prr(' - TOTAL FP NEUTRAL IONIZATIONS     : ',fp_neut_ent)
         call prr('   - FP IONIZATIONS ENTERING PLASMA : ',
     >             fp_neut_ioniz)
         call prr('   - FP IONIZATIONS REACHING WALL   : ',fp_neut_wall)
         call prr('   - FP IONIZATIONS TO FP TARGET    : ',fp_neut_targ)
         call prr('   - FP IONIZATIONS EXISITNG TO TMAX: ',fp_neut_tmax)
c
      endif
C
      SNEUT  = RNEUT (1) + RNEUT (2)
      SWALLN = RWALLN(1) + RWALLN(2)
      MTCWALLN = MTC_RWALLN(1) + MTC_RWALLN(2)
      SCENT  = RCENT (1) + RCENT (2)
      STMAX  = RTMAX (1) + RTMAX (2)
      SSTRUK = RSTRUK(1) + RSTRUK(2)
      MTCSTRUK = MTC_RSTRUK(1) + MTC_RSTRUK(2)
      SFAIL  = RFAIL (1) + RFAIL (2)
      SATIZ  = RATIZ (1) + RATIZ (2)
      SMAIN  = RMAIN (1) + RMAIN (2)
      SEXIT  = REXIT (1) + REXIT (2)
c
c     Calculate average temperature of all ions at production.
c
      eatiztot = EATIZ (1) + EATIZ (2) + eatiztot
      satiztot = satiztot + satiz
      ctemav = eatiztot / satiztot
      temav = (eatiz(1)+eatiz(2))/satiz
c

c
      CALL PRB
      CALL PRI2 ('NUMBER OF NEUTRALS IONISED         ',
     >                           NINT(RATIZ(1)),NINT(RATIZ(2)))
C
      IF (SATIZ.GT.0.0 .AND. STATUS.LE.10) THEN
      CALL PRR2 ('AVG R IONIS POSN (TARGET SIDE) (M) ',
     >                          XATIZ(1)/RATIZ(1),XATIZ(2)/RATIZ(2))
      CALL PRR2 ('AVG Z IONIS POSN (TARGET SIDE) (M) ',
     >                          YATIZ(1)/RATIZ(1),YATIZ(2)/RATIZ(2))
      CALL PRR2 ('AVG R IONIS POSN (GRID SIDE)   (M) ',
     >                      XATIZ2(1)/RATIZ2(1),XATIZ2(2)/RATIZ2(2))
      CALL PRR2 ('AVG Z IONIS POSN (GRID SIDE)   (M) ',
     >                      YATIZ2(1)/RATIZ2(1),YATIZ2(2)/RATIZ2(2))
      CALL PRR2 ('AVERAGE ANGLE AT IONISATION (DEGS) ',
     >            RADDEG*AATIZ(1)/RATIZ(1),RADDEG*AATIZ(2)/RATIZ(2))
      CALL PRR2 ('AVERAGE VELOCITY AT IONISATION(M/S)',
     >                          VATIZ(1)/RATIZ(1),VATIZ(2)/RATIZ(2))
      CALL PRR2 ('MAXIMUM VELOCITY AT IONISATION(M/S)',
     >                                VATIZM(1),VATIZM(2))
      CALL PRR2 ('AVERAGE TEMP AT IONISATION (EV)    ',
     >                          EATIZ(1)/RATIZ(1),EATIZ(2)/RATIZ(2))
      call prr  ('AVERAGE TEMP AT IONIZATION IN THIS LAUNCH    (eV)'
     >,temav)
      call prr  ('AVERAGE TEMP OF ALL PARTICLES IONIZED SO FAR (eV)'
     >,ctemav)
      CALL PRR2 ('AVERAGE TIME TO IONISATION (S)     ',
     >                          TATIZ(1)/RATIZ(1),TATIZ(2)/RATIZ(2))
      WRITE (7,'(1X,A,2F11.6)') 'AVERAGE K FOR ALL IONISATIONS      ',
     >                          KATIZ(1)/RATIZ(1),KATIZ(2)/RATIZ(2)
      CALL PRR2 ('FRACTION OF IONISATIONS INBOARD    ',
     >                            MFRAC(1)/RATIZ(1),MFRAC(2)/RATIZ(2))
      IF (MFRAC(1)+MFRAC(2).GT.0.0)
     >WRITE (7,'(1X,A,2F11.6)') 'AVERAGE K FOR MAIN IONISATIONS     ',
     >           KMATIZ(1)/MAX(LO,MFRAC(1)),KMATIZ(2)/MAX(LO,MFRAC(2))
      IF (SFRAC(1)+SFRAC(2).GT.0.0) THEN
      WRITE (7,'(1X,A,2F11.6)') 'AVERAGE K FOR SOL IONISATIONS     ',
     >           KSATIZ(1)/MAX(LO,SFRAC(1)),KSATIZ(2)/MAX(LO,SFRAC(2))
      CALL PRR2 ('AVERAGE R FOR SOL IONISATIONS     ',
     >           XSATIZ(1)/MAX(LO,SFRAC(1)),XSATIZ(2)/MAX(LO,SFRAC(2)))
      CALL PRR2 ('AVERAGE Z FOR SOL IONISATIONS      ',
     >           YSATIZ(1)/MAX(LO,SFRAC(1)),YSATIZ(2)/MAX(LO,SFRAC(2)))
      CALL PRR2 ('AVERAGE S OR SMAX-S FOR SOL IONIZ  ',
     >           SSATIZ(1)/MAX(LO,SFRAC(1)),SSATIZ(2)/MAX(LO,SFRAC(2)))
c slmod begin
      CALL inOpenInterface('idl.divimp_launch',ITF_WRITE)
      CALL inPutData(SUM(ETOT (1:2))/SUM(RNEUT(1:2)),'TAVG_PROD','eV')
      CALL inPutData(SUM(EATIZ(1:2))/SUM(RATIZ(1:2)),'TAVG_IONI','eV')
      CALL inCloseInterface
c slmod end
c
c     Double check test
c    
      call prb
      call prc ('Double Check to test ionization accumulation')
      call prb
      call prioniz(1,3,3,lionizdat)
      call prioniz(2,3,3,lionizdat)
c
c     Print Out Ionization data for the Far Periphery launched   
c     Neutrals.  
c
      if (fpropt.eq.1) then  
c
      call prc  ('FAR PERIPHERY TARGET NEUTRAL LAUNCH SUMMARIES:') 
c
c        Print SOL  
c
         call prioniz(1,2,3,lionizdat)       
c
c        Print MAIN  
c
         call prioniz(2,2,3,lionizdat)       
c
c
      endif
c
      ENDIF
      ENDIF
c
C
C     Assign weights back from snews array to sputys array. This balance
c     the accounting for particles that do not become ionized.
C
      do 9000 in = latiz,natiz-1
         sputys(in) = snews(in)
 9000 continue
C

 
      WRITE (6,'(1X,A,I15)') 'AV. ITERS PER NEUT ',NINT(CISTOT/RPROD)
      WRITE (6,'(1X,A,I15)') 'MAX ITERS ANY NEUT ',NINT(CISMAX)
      WRITE (6,'(1X,A,I15)') 'TOTAL IONS CREATED ',NINT(SATIZ)
c
      write (6,*) 
c
      write (6,'(1x,a,1x,g12.2)') 
     >               'TOTAL PARTICLES INITIALLY OUTSIDE WALL:',
     >                             err_ok+err_out
      write (6,'(1x,a,1x,g12.2)') 
     >               'TOTAL PARTICLES OK AFTER ONE MOVEMENT :',
     >                             err_ok
      write (6,'(1x,a,1x,g12.2)') 
     >               'TOTAL PARTICLES STILL OUTSIDE WALL    :',
     >                             err_out
c
C
      NEUTIM = NEUTIM + ZA02AS (1) - STATIM
      RETURN
C
 9001 FORMAT(/1X,'ETH,ETF,Q,YMF',I4,I8,F8.4,F5.2,1X,a,4x,a)
c     >  'INNER    OUTER')
 9002 FORMAT(/1X,'ETH,ETF,Q,YMF',' N/A    N/A    N/A  ',F5.2,1X,
     >     a,4x,a) 

c     >  'INNER    OUTER')
c
c     Added 9011 and 9012 to support INNER/OUTER differences 
c     between geometries.
c
 9011 FORMAT(/1X,'ETH,ETF,Q,YMF',I4,I8,F8.4,F5.2,1X,                    
     >   a,4x,a)
c
c     >  'OUTER    INNER')                                               
c
 9012 FORMAT(/1X,'ETH,ETF,Q,YMF',' N/A    N/A    N/A  ',F5.2,1X,        
     >   a,4x,a)
c
c     >  'OUTER    INNER')                                               
c
 9003 FORMAT(1X,I5,F9.1,4I4,
     >  2F9.5,F9.4,F8.1,F7.2,F5.2,F8.3,I3,:,1X,A)
 9004 FORMAT(//1X,'NEUT DEBUG: DIAGNOSTICS TO BE PRINTED EVERY',I6,
     >  ' TIMESTEPS  (DELTA T =',1P,G10.3,' SECONDS).',//)
 9005 FORMAT(1X,'-NEUT-----TIME--IK--IR--IX--IY',
     >  '-----R--------Z--------K------VIN----TEMN--FRAC--ANGLE--IT',
     >  18('-'))
      END
C
C
C
      REAL FUNCTION VLAN (RAN)
      IMPLICIT NONE
      REAL    RAN
C
C  *********************************************************************
C  *                                                                   *
C  *  VLAN:  THIS FUNCTION RETURNS THE MAXIMUM POSSIBLE NEUTRAL LAUNCH *
C  *  VELOCITY FOR A GIVEN MAXIMUM RANDOM NUMBER "RAN" <= 1.0.         *
C  *                                                                   *
C  *                                      C.M.FARRELL   NOVEMBER 1987  *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE "PARAMS"
      include    'params'
C     INCLUDE "COMTOR"
      include    'comtor'
C
      RAN = MAX (0.0001, MIN (RAN, 0.9999))
C
      IF (CNEUTC.EQ.0.OR.CNEUTC.EQ.1.OR.CNEUTC.EQ.4.OR.CNEUTC.EQ.5) THEN
        VLAN = 1.38E4 * SQRT (CEBD / (1.0/SQRT(RAN)-1.0) / CRMI)
      ELSEIF (CNEUTC.EQ.2) THEN
        VLAN = 1.38E4 * SQRT (CTEM1 * ABS(LOG(1.0-RAN)) / CRMI)
      ELSEIF (CNEUTC.EQ.3.OR.CNEUTC.EQ.6.OR.CNEUTC.EQ.7.OR.
     >        CNEUTC.EQ.8.OR.CNEUTC.EQ.10.OR.CNEUTC.EQ.11) THEN
        VLAN = 1.38E4 * SQRT (CTEM1 / CRMI)
      ELSEIF (CNEUTC.EQ.9) THEN
        VLAN = 1.38E4 * SQRT (MAX (CTEM1,CTEM2) / CRMI)
      ENDIF
      RETURN
      END
c
c
c
      subroutine tfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >               pinsw,yieldsw,matp,matt)               
c slmod begin
      USE mod_divimp
c slmod ned
      implicit none
      include 'params'
c
      real fydata(maxpts,5),fyprob(maxpts)
      real totfydata(3,5)
      integer nfy,pinsw,yieldsw,nfymap,fymap(maxpts)
      integer matp,matt 
c
      include 'cgeom'
      include 'comtor'
      include 'pindata'
      include 'cneut2'
c
c
c     TFY  : Calculates Flux, Yield and launch data onto target
c            Segments for Physical or Chemical Sputtering. 
c
c            PINSW  = 0  No PIN/NIMBUS data available 
c                   = 1  Use PIN/NIMBUS data
c                   = 2  Use Wall plasma conditions (not implemented)
c                   = 3  Use external flux file 
c                   = 4  Use impurity influx data from PIN/NIMBUS - SL, 24/11/2009
c
c
c            YIELDSW= 0  Use Physiclal Sputtering Yields
c                   = 1  Use Chemical Sputtering Yields
c
c     Contents of FYDATA on Output:
c
c     FYDATA(ID,1)  = Flux
c               2   = Energy
c               3   = Heat
c               4   = Yield
c               5   = Flux * Yield
c
c     TOTFYDATA(N,1)= Total Flux
c                 2 =  
c                 3 = Total Heat
c                 4 = 
c                 5 = Total Flux * Yield 
c         
c                 N = 1 Inner 
c                     2 Outer
c                     3 Total 
c
c
c     Local Variables
c
      integer id,ik,ir,in,targ,ierr
      real cs,sheath_te,sheath_ti
      real ftot,fytot,fhtot,yldchem,yield
      external yldchem,yield
      integer ringno
      external ringno
c
c     Initialize 
c
      nfy = nds
      call rzero (fydata,maxpts*5)
c
c     Record Launch Type in unused space - for use in print routine
c
      totfydata(3,2) = 0.0
      totfydata(3,4) = yieldsw
c
c     Select work to be done based on pinsw 
c
      if (pinsw.eq.0.or.pinsw.eq.1.or.
     .    (pinsw.eq.4.and.yieldsw.eq.1)) then  
c
c     -------------------------------------------------------
C     
C---- CALCULATE FLUXES AND YIELDS
C     
      IF (sloutput) WRITE(0,*) 'CALL TFY: looping over targets...',
     .                         yieldsw,pinsw
      DO ID = 1,NDS
        IK = IKDS(ID)
        IR = IRDS(ID)
c
c       Set target identifier - TARG = 1 = OUTER for SONNET
c                                          INNER for GRID2D   
c                               TARG = 2 = INNER for SONNET
c                                          OUTER for GRID2D   
        if (ik.gt.nks(ir)/2) then 
           targ = 1
        else 
           targ = 2
        endif 
c
c       Assign Sheath quantities
c
        if ((targ.eq.1.and.nsheath_vali.eq.0).or.
     >      (targ.eq.2.and.nsheath_valo.eq.0)) then 
           sheath_te = kteds(id)
           sheath_ti = ktids(id)
        else
c sltmp
          IF (sloutput) STOP 'STOP: SHOULD NOT USE SHEATH VALUES'  
c
c          The code assumes for now that Te=Ti when the
c          sheath potential is directly specified. 
c
           if (targ.eq.1) then 

              in = RINGNO(ir,sheath_vali,nsheath_vali,
     >                    maxnrs,2,ierr)
              if (ierr.eq.0) then 
                 sheath_te = sheath_vali(in,2)
                 sheath_ti = sheath_vali(in,2)
              else
                 sheath_te = kteds(id)
                 sheath_ti = ktids(id)
              endif              

           elseif (targ.eq.2) then 

              in = RINGNO(ir,sheath_valo,nsheath_valo,
     >                    maxnrs,2,ierr)
              if (ierr.eq.0) then 
                 sheath_te = sheath_valo(in,2)
                 sheath_ti = sheath_valo(in,2)
              else
                 sheath_te = kteds(id)
                 sheath_ti = ktids(id)
              endif              
c
           endif

c           
        endif
c     
c       Calculate speed at target 
c     

c
c        if (cioptf.eq.22.or.cioptf.eq.23) then
c           if (ik.lt.nks(ir)/2) then
c              CS = 9.79E3 * SQRT (0.5*(KTEDS(ID)+KTIDS(ID))*
c     >                (1.0+RIZB)/CRMB) * cmachno(ir,2)
c           elseif  (ik.ge.nks(ir)/2) then
c              CS = 9.79E3 * SQRT (0.5*(KTEDS(ID)+KTIDS(ID))*
c     >                (1.0+RIZB)/CRMB) * cmachno(ir,1)
c           endif
c        else
c           CS = 9.79E3 * SQRT (0.5*(KTEDS(ID)+KTIDS(ID))*
c     >                (1.0+RIZB)/CRMB)
c        endif
c
        cs = abs(kvds(id))
c     
c       Calculate flux
c     
        if (northopt.eq.0) then
c     
          fydata(ID,1)   = KNDS(ID) * CS / KBFS(IK,IR)
c     
        elseif (northopt.eq.1.or.northopt.eq.2
     >          .or.northopt.eq.3) then
c     
          fydata(ID,1)   = KNDS(ID) * CS/KBFS(IK,IR) * COSTET(ID)  ! *** Assumption here that M=1 I think... ***
c     
        endif
c     
        IF (CNEUTD.EQ.1) THEN
c
c          fydata(ID,2) = 2.0*KTIDS(ID) + 3.0*REAL(CBOMBZ)*KTEDS(ID)
c          fydata(id,3) = fydata(id,2) + 2.0*REAL(CBOMBZ)*KTEDS(ID)
c
          fydata(ID,2) = 2.0*sheath_ti + 3.0*REAL(CBOMBZ)*sheath_te
          fydata(id,3) = fydata(id,2) + 2.0*REAL(CBOMBZ)*sheath_te
c
        ELSE
c
c          fydata(id,2) = 2.0*KTIDS(ID) + 3.0*RIZB*KTEDS(ID)
c          fydata(id,3) = fydata(id,2) + 2.0*RIZB*KTEDS(ID)
c
          fydata(id,2) = 2.0*sheath_ti + 3.0*RIZB*sheath_te
          fydata(id,3) = fydata(id,2) + 2.0*RIZB*sheath_te
c
        ENDIF
c
c
        if (fydata(id,1).gt.0.0.and.fydata(id,2).gt.0.0) then
c     
c          Physical Sputtering  
c     
           if (yieldsw.eq.0) then
c     
              fydata(id,4)  = YIELD (MATP,MATT,fydata(id,2),
     >                            sheath_te,sheath_ti) * KMFPS(ID)
c     
c          Chemical Sputtering 
c     
           elseif (yieldsw.eq.1) then
c     
              if (pinsw.eq.0) then  
c     
                 fydata(ID,4) = YLDCHEM (fydata(id,2),fydata(id,1),
     >                           MATP,MATT,tempds(id))*KMFCS(ID)
c     
              elseif (pinsw.eq.1.or.pinsw.eq.4) then 
c     
c                Test for valid index to NIMBUS data
c     
                 if (nimindex(id).eq.0) then 
c
c                    fydata(id,4) = 0.0
c
                    fydata(id,4) =  YLDCHEM (fydata(id,2),fydata(id,1),
     >                           MATP,MATT,tempds(id))*KMFCS(ID)
c
                 else
c     
c                   Use TOTAL hydrogenic flux for yield calculation 
c
                    fydata(ID,4) = YLDCHEM (fydata(id,2),
     >                           flxhw2(nimindex(id)),
     >                           MATP,MATT,tempds(id))*KMFCS(ID)
                 endif 
c     
              endif
      
           endif

        else
c
c          IF NO flux or No Energy - set yield = 0
c
           fydata(id,4) = 0.0
c
        endif    
c     
        fydata(id,5)  = fydata(id,1) * fydata(id,4)
c     
        if (cprint.eq.5.or.cprint.eq.9) then 
           write(6,'(a,3i8,10(1x,g18.8))') 'id:',id,ik,ir,
     >              kteds(id),ktids(id),knds(id),
     >              kbfs(ik,ir),kmfps(id),sheath_te,sheath_ti
        endif 
c
      end do 
c
c     -------------------------------------------------------
c
c     External fluxes specified and used  
c
      elseif (pinsw.eq.3) then 
c
c        Load external flux data array 
c
         call calc_extflx_yield(fydata,matt)
c
c slmod begin
c
c Put this back in when I start processing the impurity production in EIRENE
c that's coming from ions striking the target -- need to compare with what 
c DIVIMP is calculating.  Need to remove the pinsw.eq.4 references that occur
c in this subroutine, see above.
c
      elseif (pinsw.eq.4.and.yieldsw.EQ.0) then 
        IF (sloutput) WRITE(0,*) 'WHA-WHO! 2'
          IF (.NOT.ALLOCATED(wall_flx)) 
     .      CALL ER('TFY','WALL_FLX not allocated and PINSW=4')
c       FYDATA(ID,1)  = Flux
c                 2   = Energy
c                 3   = Heat
c                 4   = Yield
c                 5   = Flux * Yield
        DO id = 1, nds
          in = nimindex(id)
          IF (in.EQ.0) CYCLE  ! virtual rings 
          fydata(id,1) = wall_flx(in)%in_par_blk(1,0)    ! flxhw2(in)  ! FLUX OF HYDROGEN (ATOMS AND IONS) TO THE WALL
          fydata(id,2) = wall_flx(in)%in_ene_blk(1,0)    ! flxhw5(in)  ! AVERAGE ENERGY OF ATOMS HITTING THE WALL (EV)
          fydata(id,3) = 1.0 
          IF (wall_flx(in)%in_par_blk(1,0).NE.0.0)
     .      fydata(id,4) = wall_flx(in)%em_par_atm(2,1) / 
     .                     wall_flx(in)%in_par_blk(1,0)  ! flxhw3(in) / (flxhw2(in) + 1.0E-10)
          fydata(id,5) = wall_flx(in)%em_par_atm(2,1)    ! flxhw3(in)  ! Atoms (species=2) sputtering by bulk ions
          fydata(id,:) = fydata(id,:) * kmfps(id)

c          write(0,*) '  debug 1: fydata5',id,fydata(id,5),kmfps(id)
        enddo
c slmod end
      endif 
C     
C---- GENERATE CUMULATIVE DISTRIBUTION FUNCTIONS FOR NEUTRAL
C---- PRODUCTION
C     
c     
c     Initialize
c     
      FTOT  = 0.0
      FYTOT = 0.0
      FHTOT = 0.0
c     
      DO ID = 1, NDS
        FTOT  = FTOT  + fydata(id,1) * DDS(ID)
        FHTOT = FHTOT + fydata(id,1) * fydata(id,3) * DDS(ID)
        FYTOT = FYTOT + fydata(id,5) * DDS(ID)
c     
c       Inner totals
c     
        IF (ID.EQ.NDSIN) THEN
          TOTfydata(1,1) = FTOT
          totfydata(1,3) = FHTOT
          totfydata(1,5) = FYTOT
        ENDIF
c     
      end do 
c     
c     Grand totals 
c     
      TOTfydata(3,1) = FTOT
      totfydata(3,3) = FHTOT
      totfydata(3,5) = FYTOT
c     
c     Outer totals
c     
      TOTfydata(2,1) = TOTfydata(3,1)- TOTfydata(1,1)
      totfydata(2,3) = totfydata(3,3)- totfydata(1,3)
      totfydata(2,5) = totfydata(3,5)- totfydata(1,5)

      IF (sloutput) WRITE(0,*) 'CALL TFY: total yield =',totfydata(3,5)
c     
c     Calculate cumulative - normalized - non-zero launch 
c     probabilities. 
c     
c     
      nfymap = 0 
      call rzero(fyprob,maxpts)
c
      do id = 1,nds
c
c        If non-zero launch probability - add the point to map
c        and cumulative launch array - as well as accumulate probabilty
c        total.
c
         if ((fydata(id,5)*dds(id)).gt.0.0) then
            nfymap = nfymap+1
            fymap(nfymap) = id
            if (nfymap.eq.1) then   
               fyprob(nfymap) = (fydata(id,5) * dds(id))/totfydata(3,5) 
            else
               fyprob(nfymap) = fyprob(nfymap-1) 
     >                        + (fydata(id,5) * dds(id))/totfydata(3,5) 
            endif
         endif
      end do
c
c     Print out data for debugging
c      
      if (cprint.eq.3.or.cprint.eq.9) then

         write (6,*) 'FYDATA:',pinsw,yieldsw,nds
         write (6,'(4x,9x,''Flux'',7x,''Energy'',
     >           9x,''Heat'',8x,''Yield'',3x,''Flux*Yield'')') 
c
         do id = 1,nds
c
            write (6,'(i4,5g13.5)') id,(fydata(id,in),in=1,5) 
c   
         end do 
c
         write (6,*) 'Totals:'
         do id = 1,3
            write (6,'(i4,5g13.5)') id,(totfydata(id,in),in=1,5) 
         end do    
c
         write (6,*) 'MAP:',nfymap
         write (6,'(4x,''Index'',4x,''Cum. Prob'')')
c
         do id = 1,nfymap
c         
            write (6,'(2i4,1g13.5)') id,fymap(id),fyprob(id)
c
         end do 
c
      endif
c
c     Exit 
c
      return
      end
c
c
c
      subroutine wfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >               pinsw,yieldsw,matp,matt)               
c slmod begin
      use mod_divimp
c slmod end
      implicit none
      include 'params'
c
      real fydata(maxpts,5),fyprob(maxpts)
      real totfydata(3,5)
      integer nfy,pinsw,yieldsw,nfymap,fymap(maxpts)
      integer matp,matt 
c
      include 'cgeom'
      include 'comtor'
      include 'pindata'
      include 'cneut2'
      include 'slcom' 
c slmod begin - tmp
      LOGICAL warning, bug_message
      DATA    warning, bug_message /.FALSE., .TRUE./
      SAVE
c slmod end 
c
c
c     WFY  : Calculates Flux, Yield and launch data onto Wall
c            Segments for Physical or Chemical Sputtering. 
c
c            PINSW  = 0  No PIN/NIMBUS data available 
c                   = 1  Use PIN/NIMBUS data
c                   = 2  Use WALL plasma data
c                   = 3  Use external fluxes (not implemented)
c                   = 4  Use impurity influx data from PIN/NIMBUS - SL, 24/11/2009
c
c            YIELDSW= 0  Use Physical Sputtering Yields
c                   = 1  Use Chemical Sputtering Yields
c
c     Contents of FYDATA on Output:
c
c     FYDATA(ID,1)  = Flux
c               2   = Energy
c               3   = Heat
c               4   = Yield
c               5   = Flux * Yield
c
c     TOTFYDATA(N,1)= Total Flux
c                 2 =                - Indicates Target (1) or Wall (2)
c                 3 = Total Heat
c                 4 =                - Indicates Physical (0) or Chemical (1)
c                 5 = Total Flux * Yield 
c         
c                 N = 1 Inner 
c                     2 Outer
c                     3 Total 
c
c     Local Variables
c
      integer id,ik,ir,in
      real cs,ne,te,ti,ionflux,tsurf
      real ftot,fytot,fhtot,yldchem,yield
      external yldchem,yield
c
c     Initialize 
c
      write (6,*) 'WFY:',matt,pinsw,wallpts
c
c     Zero FLUX and YIELD array
c
      call rzero(fydata,maxpts*5)
      call rzero(totfydata,3*5)
c
      nfy = wallpts
c
c     Record Launch Type in unsed space - for use in print routine
c
      totfydata(3,2) = 1.0
      totfydata(3,4) = yieldsw
c
c     This routine requires PINSW to be on - since there is NO flux/yield
c     data without it - however - it will return the raw wall launch
c     information that is usually used for Launch Option 2 - this information
c     is copied out of the arrays fwlprob, wlind, and nwlind - which are 
c     calculated in WALLS.  
c
      if (pinsw.eq.0.or.nwlprob.gt.0) then 
c
c        Load basic Wall launch data.
c
c        Set non-zero totfy in order to allow a forced launch even though
c        no fy data is available. 
c
         totfydata(3,5) = 1.0
c
         nfymap = nwlind
c
         do id = 1,nwlind
            fymap(id) = wlind(id)            
            fyprob(id)= fwlprob(id) 
c slmod begin
            if (wallpt(wlind(id),7).eq.0.0) cycle

            if (id.eq.1.and.bug_message) then
 
              write(0,*) 
              write(0,*) '---------------------------------------------'
              write(0,*) '  WARNING! Bug correction related to assign'//
     .                   'ed wall launch probabilities, see neut.f'
              write(0,*) '---------------------------------------------'
              write(0,*) 

              write(0,*) 
              write(0,*) '---------------------------------------------'
              write(0,*) '  WARNING! Bug correction related to assign'//
     .                   'ed wall launch probabilities, see neut.f'
              write(0,*) '---------------------------------------------'
              write(0,*) 

              bug_message = .FALSE.
            endif

c...        The direct assignment of FWLPROB to FYDATA(,5) is not correct
c           since FWLPROB is an additive measure of wall launch probability,
c           with the last entry approaching 1.0, while FYDATA(,5) should contain
c           individual wall launch probabilities for each wall segment that are
c           independent of the lauch probabilities for the preceeding segments
c           in the list. - SL, 09/05/12

            if (id.eq.1) then
              fydata(wlind(id),5) = fwlprob(id) 
            else
              fydata(wlind(id),5) = fwlprob(id) - fwlprob(id-1)
            endif

            fydata(wlind(id),5) = fydata(wlind(id),5) * nabsfac / 
     .                            wallpt(wlind(id),7)

c
c     .        fydata(wlind(id),5) = fwlprob(id) * nabsfac / 
c     .                              wallpt(wlind(id),7)
c slmod end
         end do 
c
c     Do proper calculations if PIN/NIMBUS data is available
c
      elseif (pinsw.eq.1.or.pinsw.eq.2.or.pinsw.eq.4) then
c
         if (pinsw.eq.1) then 
C     
C----       CALCULATE FLUXES AND YIELDS - these are ONLY atomic fluxes
c           except for target wall elements where we need to extract the 
c           atomic component because it has a different energy spectrum
c           than the primary ions - these are taken care of in the target
c           yield routine. 
C     
            DO IN = 1,wallpts
c       
c             Calculate flux - must be taken from PIN/NIMBUS
c     
c
c             Target points
c
              if (wallpt(in,16).eq.1.or.wallpt(in,16).eq.4) then   
c
                 id = wallpt(in,18)
c
c                Modify the flux calculation removing the ion 
c                component. 
c
                 IK = IKDS(ID)
                 IR = IRDS(ID)
c     
c                Calculate speed at target 
c     
c
                 cs = abs(kvds(id))
c     
c                Calculate ion flux
c     
                 if (northopt.eq.0) then
c     
                    ionflux = KNDS(ID) * CS / KBFS(IK,IR)
c     
                 elseif (northopt.eq.1.or.northopt.eq.2
     >                   .or.northopt.eq.3) then
c     
                    ionflux = KNDS(ID) * CS/KBFS(IK,IR)*COSTET(ID)
c     
                 endif
c              
c                Calcuate NET flux
c
c                Minimum of 0.0
c
c                For Eirene99 use the atom flux directly reported by 
c                Eirene instead of deriving it.  
c  
c slmod begin
                 if (pincode.eq.2.or.pincode.eq.3.or.pincode.eq.4.or.
     >               pincode.eq.5) then 
c
c                 if (pincode.eq.2.or.pincode.eq.3) then 
c slmod end

                    fydata(in,1)=max(flxhw6(int(wallpt(in,17))),0.0)
c
                 else
c
                    fydata(in,1)=max(flxhw2(int(wallpt(in,17)))
     >                                    -ionflux,0.0)
                 endif 
c
c 
c             Regular Wall points
c
              else
c
                 fydata(In,1) = flxhw2(int(wallpt(in,17)))
c
              endif
c     
              fydata(IN,2) = flxhw5(int(wallpt(in,17)))
              fydata(in,3) = flxhw5(int(wallpt(in,17)))
c
c
              if (fydata(in,1).gt.0.0.and.fydata(in,2).gt.0.0) then 
c     
c                Physical Sputtering  
c      
                 if (yieldsw.eq.0) then
c     
c                   Supply plasma temperature at wall as zero for
c                   now - add check for zero temperature in yield so
c                   that regular yield code is used. 
c
                    fydata(in,4)  = YIELD (MATP,MATT,fydata(in,2),
     >                                  0.0,0.0) 
     >                               * kmfpws(in)
c     
c                Chemical Sputtering 
c     
                 elseif (yieldsw.eq.1) then
c
c                   Each segment now has a unique temperature 
c
                    tsurf = wallpt(in,19)
c
c     
c                   Use TOTAL hydrogenic flux for yield calculation 
c
                    fydata(In,4) = YLDCHEM(fydata(in,2),
     >                           flxhw2(int(wallpt(in,17))),
     >                           MATP,MATT,tsurf)* kmfcws(in)
c
                 endif
c
              else
c
c                IF NO flux or No Energy - set yield = 0
c
                 fydata(in,4) = 0.0
c
              endif
c     
              fydata(in,5)  = fydata(in,1) * fydata(in,4)
c     
            end do 
c
c        Deal with WALL PLASMA SPUTTERING - PHYSICAL
c 
c        - this code applies only to the option where a plasma has
c          been explicitly set for specific wall elements
c
         elseif (pinsw.eq.2) then 
c         
C     
C----       CALCULATE FLUXES AND YIELDS - these are based on the wall
c           plasma conditions specified by the wall_plasma_opt option
c           except for target wall elements where the actual target 
c           conditions are used.
C     
            DO IN = 1,wallpts
c       
c             Calculate flux - must be taken from PIN/NIMBUS
c
c             Target points
c
              if (wallpt(in,16).eq.1.or.wallpt(in,16).eq.4) then   
c
                 id = wallpt(in,18)
c
c                Modify the flux calculation removing the ion 
c                component. 
c
                 IK = IKDS(ID)
                 IR = IRDS(ID)
c              
c                Assign plasma conditions
c     
c                Calculate speed at target 
c
                 ne = knds(id)
                 te = kteds(id)
                 ti = ktids(id)
                 cs = abs(kvds(id))
c     
c                Calculate ion flux
c     
                 if (northopt.eq.0) then
c     
                    ionflux = KNDS(ID) * CS / KBFS(IK,IR)
c     
                 elseif (northopt.eq.1.or.northopt.eq.2
     >                   .or.northopt.eq.3) then
c     
                    ionflux = KNDS(ID) * CS/KBFS(IK,IR)*COSTET(ID)
c     
                 endif
c                 
              else
c              
c                Assign plasma conditions
c     
c                Calculate speed at target 
c
                 ne = wallpt(in,31)
                 te = wallpt(in,29)
                 ti = wallpt(in,30)
                 cs = 9.82e3*sqrt(2.0*te/crmi) 
                 ionflux = ne * cs
c
              endif
c
c             Calculate FY data     
c
              fydata(in,1) = ionflux 
              fydata(in,2) = 2.0*ti+3.0*te
              fydata(in,3) = fydata(in,2)+2.0*te
c
c             Calculate yield
c
              fydata(in,4)=yield(matp, matt, fydata(in,2),
     >                           ti, te) * kmfpws(in)
c
              fydata(in,5)=fydata(in,1)*fydata(in,4)
c
            end do
c
c slmod begin
         elseif (pinsw.eq.4.and.yieldsw.eq.0) then 
           IF (sloutput) WRITE(0,*) 'WHA-WHO! 1'
           IF (.NOT.ALLOCATED(wall_flx)) 
     .       CALL ER('WFY','WALL_FLX not allocated and PINSW=4')
c          FYDATA(ID,1)  = Flux
c                    2   = Energy
c                    3   = Heat
c                    4   = Yield
c                    5   = Flux * Yield
           do in = 1, wallpts
             id = NINT(wallpt(in,17))
             fydata(in,1) = flxhw2(id)  ! FLUX OF HYDROGEN (ATOMS AND IONS) TO THE WALL
             fydata(in,2) = flxhw5(id)  ! AVERAGE ENERGY OF ATOMS HITTING THE WALL (EV)
             fydata(in,3) = 1.0 
             IF (wall_flx(in)%in_par_atm(1,0).NE.0.0)
     .         fydata(in,4) = wall_flx(id)%em_par_atm(2,2) /
     .                        wall_flx(in)%in_par_atm(1,0)
             fydata(in,5) = wall_flx(id)%em_par_atm(2,2)  !  Atoms (species=2) sputtering by test atoms
c             fydata(in,4) = MAX(0.0,flxhw3(id) / (flxhw2(in) + 1.0E-10))  ! bug, FLXHW2(IN) should have been FLXHW6(ID) 
c             fydata(in,5) = MAX(0.0,flxhw3(id))  ! FLUX OF IMPURITIES SPUTTERED FROM THE WALL 
             fydata(id,:) = fydata(id,:) * kmfpws(id)

c          write(0,*) '  debug 2: fydata5',in,fydata(in,5)
c             IF (kmfpws(id).NE.0.0) THEN
c               WRITE(0,*) '  DEGUB: fydata4,5=',fydata(in,4:5),
c     .                     MAX(0.0,flxhw3(id) / (flxhw6(id) + 1.0E-10)),
c     .                     MAX(0.0,flxhw3(id))
c             ENDIF
           enddo
c slmod end
         endif 

C     
C----    GENERATE CUMULATIVE DISTRIBUTION FUNCTIONS FOR NEUTRAL
C----    PRODUCTION
C     
c     
c        Initialize
c     
         FTOT  = 0.0
         FYTOT = 0.0
         FHTOT = 0.0
c     
         DO IN = 1, WALLPTS
           FTOT  = FTOT  + fydata(in,1) * wallpt(in,13)
           FHTOT = FHTOT + fydata(in,1) * fydata(in,3) * wallpt(in,13)
           FYTOT = FYTOT + fydata(in,5) * wallpt(in,13)
c     
c          Inner totals - not meaningful for WALLS 
c     
           IF (IN.EQ.(WALLPTS/2)) THEN
              TOTfydata(1,1) = FTOT
              totfydata(1,3) = FHTOT
              totfydata(1,5) = FYTOT
           ENDIF
c     
         end do 
c     
c        Grand totals 
c     
         TOTfydata(3,1) = FTOT
         totfydata(3,3) = FHTOT
         totfydata(3,5) = FYTOT
c     
c        Outer totals
c     
         TOTfydata(2,1) = TOTfydata(3,1)- TOTfydata(1,1)
         totfydata(2,3) = totfydata(3,3)- totfydata(1,3)
         totfydata(2,5) = totfydata(3,5)- totfydata(1,5)
c     
c        Calculate cumulative - normalized - non-zero launch 
c        probabilities. 
c     
c     
         nfymap = 0 
         call rzero(fyprob,maxpts)
c
         do in = 1,wallpts
c
c           If non-zero launch probability - add the point to map
c           and cumulative launch array - as well as accumulate probabilty
c           total.
c
            if ((fydata(in,5)*wallpt(in,13)).gt.0.0) then
               nfymap = nfymap+1
               fymap(nfymap) = in
               if (nfymap.eq.1) then   
                  fyprob(nfymap) = (fydata(in,5)*wallpt(in,13))
     >                              /totfydata(3,5) 
               else
                  fyprob(nfymap) = fyprob(nfymap-1) 
     >                           + (fydata(in,5)*wallpt(in,13))
     >                             /totfydata(3,5) 
               endif
            endif
         end do
c
c     End of PINSW IF statement
c
      endif

      IF (sloutput) WRITE(0,*) 'CALL WFY: total yield =',totfydata(3,5)

c
c     Print out data for debugging
c      
      write (6,*) 'FYDATA:',pinsw,yieldsw,nfy
      write (6,'(4x,9x,''Flux'',7x,''Energy'',
     >           9x,''Heat'',8x,''Yield'',3x,''Flux*Yield'')') 
c
      do id = 1,nfy
c
         write (6,'(i4,5g13.5)') id,(fydata(id,in),in=1,5) 
c   
      end do 
c
      write (6,*) 'Totals:'
      do id = 1,3
         write (6,'(i4,5g13.5)') id,(fydata(id,in),in=1,5) 
      end do    
c
      write (6,*) 'MAP:',nfymap
      write (6,'(4x,''Index'',4x,''Cum. Prob'')')
c
      do id = 1,nfymap
c         
         write (6,'(2i4,1g13.5)') id,fymap(id),fyprob(id)
c
      end do 
c
c     Exit 
c
      return
      end
c
c
c
      subroutine neutbatch(newcneutc,newcneutb,yieldsw,pinsw,
     >                     nneut1,nneut2,nion1,nion2,gambl,
     >                     RSTRUKa,MTCSTRUKA,RMAINa,REXITa,RATIZa,
     >                     RNEUTa,RWALLNa,MTCWALLNA,RCENTa,RTMAXa,
     >                     rfaila,
     >                     SEED,NRAND,NEUTIM,STATUS,MATP,MATT,
     >                     fymap,fyprob,nfymap,fydata)
      use velocity_dist
      implicit none 
c
      include    'params'
      include    'hc_global_opts'
c
      integer  newcneutc,newcneutb,yieldsw,pinsw,	
     >         nneut1,nneut2,nion1,nion2,nrand,status,matp,matt,
     >         nfymap,fymap(maxpts)
      real     RSTRUKa,RMAINa,REXITa,RATIZa,RNEUTa,RWALLNa,
     >         RCENTa,RTMAXa,rfaila,fyprob(maxpts),fydata(maxpts,5),
     >         gambl 
      real     mtcstruka, mtcwallna  
      double precision seed 
c
c
c     NEUTBATCH: This routine calls the launch procedure with a 
c                Batch of neutrals from nneut1+1 to nneut1+nneut2 using
c                the velocity/angle flag and launch option that 
c                have been passed to it. 
c
c
c     Input:  yieldsw = yield switch : 0=physical 1=chemical
c             pinsw   = PIN switch   : 0=no PIN data 1=PIN data available
c
c
c
      include    'dynam1'
      include    'dynam3'
      include    'dynam4'
      include    'cyield'
      include    'comtor'
      include    'cgeom'
      include    'cioniz'
      include    'cneut'
      include    'cneut2'
c slmod begin
      include    'slcom'
c slmod end
c
c     Local variables
c
      integer ik,ir,itmp
      integer tmpcneutc,tmpcneutb,nprod,ipos,id,iprod,neuttype
      external ipos
      real neutim
      real emax,wlximax,ran
      real x0, y0
      logical nonzero_krmax
c
c     Set CNEUTC and CNEUTB to the values passed in.
c

c      write(0,*) 'Entering NEUTBATCH:'


      tmpcneutb = cneutb
      tmpcneutc = cneutc
c
      cneutb = newcneutb
      cneutc = newcneutc
c
C     RECALCULATE FOR THE NEW VALUE OF CNEUTC.
c     
C---- SET UP LIMITING RANDOM NUMBERS FOR CALCULATING LAUNCH VELOCITY.
c
c     For target launch - physical sputtering launch
C     
      nonzero_krmax = .false.
c
      if (cneutb.eq.0.or.cneutb.eq.3) then 
         IF (CNEUTC.EQ.1) CEMAXF = 1.0
C     
         DO ID = 1, NDS
            KRMAX(ID) = 1.0
c
c           Limits for cutting off the tail of physical sputtering energies 
c
            IF (CNEUTC.EQ.1 .OR. CNEUTC.EQ.4 .OR. CNEUTC.EQ.5) THEN
               IF (CNEUTD.EQ.1) THEN
                  EMAX = CEMAXF * fydata(id,2)
               ELSEif (cneutd.ne.1) then
                  if (northopt.eq.0.or.northopt.eq.2) then
                     EMAX = CEMAXF * (fydata(id,2) * GAMBL - CEBD)
                  elseif (northopt.eq.1.or.northopt.eq.3) then
                     if (matt.le.ntars) then 
                        EMAX = CEMAXF*CEBD * (fydata(id,2)
     >                              /CETH(MATP,MATT) - 1.0)
                     else
                        EMAX = CEMAXF * (fydata(id,2) * GAMBL - CEBD)
                     endif 
                  endif
               ENDIF
               IF (EMAX.GT.0.0) THEN
                  KRMAX(ID) = 1.0 / ((1.0+CEBD/EMAX) * (1.0+CEBD/EMAX))
               ELSE
                  KRMAX(ID) = 0.0
               ENDIF
            ENDIF
c
c           Check KRMAX to make sure that at least one element is non-zero 
c       
            if (krmax(id).gt.0.0) nonzero_krmax = .true. 
c
         end do 
c
      elseif (cneutb.eq.2.or.cneutb.eq.4) then 
c
         IF (CNEUTC.EQ.1) CEMAXF = 1.0
C     
         DO ID = 1, WALLPTS
            KRMAXW(ID) = 1.0
            IF (CNEUTC.EQ.1 .OR. CNEUTC.EQ.4 .OR. CNEUTC.EQ.5) THEN
c slmod begin
c...          Moved this up here and added the check for manual specification
c             of the neutral wall sputtering, as in WFY, so that the code 
c             doesn't complain if PIN was run but manual settings were
c             specified:  -SL, 09/09/20
c
c              PIN data not available 
c
               if (pinsw.eq.0.or.nwlprob.gt.0) then 
                  if (northopt.eq.0.or.northopt.eq.2) then
                     EMAX = CEMAXF * (CEIMP * GAMBL - CEBD)
                  elseif (northopt.eq.1.or.northopt.eq.3) then
                     if (matt.le.ntars) then 
                        EMAX = CEMAXF * CEBD * 
     >                         (CEIMP / CETH(MATP,MATT) - 1.0)
                     else
                        EMAX = CEMAXF * (CEIMP * GAMBL - CEBD)
                     endif
                  endif
c
c              Pin data available
c
               elseif (pinsw.eq.1.or.pinsw.eq.4) then 
c
                  IF (CNEUTD.EQ.1) THEN
                     EMAX = CEMAXF * fydata(id,2)
                  ELSEif (cneutd.ne.1) then
                    if (northopt.eq.0.or.northopt.eq.2) then
                        EMAX = CEMAXF * (fydata(id,2) * GAMBL - CEBD)
                     elseif (northopt.eq.1.or.northopt.eq.3) then
                        if (matt.le.ntars) then 
                           EMAX = CEMAXF*CEBD * (fydata(id,2)
     >                              /CETH(MATP,MATT) - 1.0)
                        else
                           EMAX = CEMAXF * (fydata(id,2) * GAMBL - CEBD)
                        endif 
                     endif
                  ENDIF
               endif
c
cc
cc              PIN data not available 
cc
c               elseif (pinsw.eq.0) then          
c                  if (northopt.eq.0.or.northopt.eq.2) then
c                     EMAX = CEMAXF * (CEIMP * GAMBL - CEBD)
c                  elseif (northopt.eq.1.or.northopt.eq.3) then
c                     if (matt.le.ntars) then 
c                        EMAX = CEMAXF * CEBD * 
c     >                         (CEIMP / CETH(MATP,MATT) - 1.0)
c                     else
c                        EMAX = CEMAXF * (CEIMP * GAMBL - CEBD)
c                     endif
c                  endif
c               endif
c slmod end

               IF (EMAX.GT.0.0) THEN
                  KRMAXW(ID) = 1.0 / ((1.0+CEBD/EMAX) * (1.0+CEBD/EMAX))
               ELSEif (cebd.eq.0.0) then 
                  krmaxw(id) = 1.0 
               else 
                  KRMAXW(ID) = 0.0
               ENDIF
            ENDIF
c
c           Check KRMAXW to make sure that at least one element is non-zero 
c       
            if (krmaxw(id).gt.0.0) nonzero_krmax = .true. 

         end do 
c
      ENDIF
c
c     IF there are no nonzero elements of KRMAX - i.e. they are all zero
c     then the code needs to issue an error message and exit.
c
c
      if (.not.nonzero_krmax.and.(cneutb.ne.1.and.
c slmod begin
     >         cneutb.ne.6.and.cneutb.ne.7).and.
     >         .not.((cneutb.eq.2.or.cneutb.eq.4).and.pinsw.eq.4)) then 
c
c     >         cneutb.ne.6.and.cneutb.ne.7)) then 
c slmod end
c
         write(6,'(a)') 'FATAL ERROR IN NEUTBATCH ROUTINE:'
         write(6,'(a)') '- LIMITING ENERGY VALUE FOR LAUNCH IS'//
     >                    ' ZERO FOR ALL LAUNCH ELEMENTS'
         write(6,'(a)') '- LAUNCH IS NOT POSSIBLE - CODE EXITING'
         write(6,'(a)') '- CHECK SELECTION OF VELOCITY/ANGLE FLAG'//
     >                  ' AND CORRESPONDING SPUTTER ENERGIES'
c
         write(0,'(a)') 'FATAL ERROR IN NEUTBATCH ROUTINE:'
         write(0,'(a)') '- LIMITING ENERGY VALUE FOR LAUNCH IS'//
     >                    ' ZERO FOR ALL LAUNCH ELEMENTS'
         write(0,'(a)') '- LAUNCH IS NOT POSSIBLE - CODE EXITING'
         write(0,'(a)') '- CHECK SELECTION OF VELOCITY/ANGLE FLAG'//
     >                  ' AND CORRESPONDING SPUTTER ENERGIES'
c sltmp
c         write(0,*) 'CNEUTB=',cneutb
c         write(0,*) 'PINSW =',pinsw
c         write(0,*) 'CEBD  =',cebd

         stop 1
c
      endif

c     
c     Calculate launch positions of particles.
c     
      CALL SURAND (SEED, Nneut2, RANVA)
      NRAND = NRAND + Nneut2

      CALL SURAND (SEED, Nneut2, RANVB)
      NRAND = NRAND + Nneut2

c      DO IPROd = 1,nneut2
c         write(6,'(a,i8,10(1x,g18.8))')
c     >      'ranv:',iprod,ranva(iprod),ranvb(iprod)
c         write(0,'(a,i8,10(1x,g18.8))')
c     >      'ranv:',iprod,ranva(iprod),ranvb(iprod)
c      end do


C     
      DO IPROD = 1,nneut2
         RAN    = RANVA (IPROD)
C
C------ DEPENDING ON OPTION CHOSEN,  SELECT LAUNCH POSITIONS
C
c         WRITE(0,*) 'ran:',iprod,ran

         IF (CNEUTB.EQ.0.or.cneutb.eq.3) THEN
c
  485       ID = IPOS (RAN, fyprob, nfymap)
            ID = fyMAP(ID)


c            WRITE(0,'(a,5i8,10(1x,g18.8))') 
c     >           'ran2:',cneutb,cneutf,iprod,id,nfymap,krmax(id),
c     >                    ran,rp(id),zp(id)
c            WRITE(6,'(a,5i8,10(1x,g18.8))') 
c     >           'ran2:',cneutb,cneutf,iprod,id,nfymap,krmax(id),
c     >                    ran,rp(id),zp(id)
c
            IF ( KRMAX(ID).LE.0.0.and.yieldsw.eq.0) THEN
              CALL SURAND2 (SEED, 1, RAN)
              NRAND = NRAND + 1
              GOTO 485
            ENDIF
c
c           Allow for neut spreading - if turned on can launch from any
c           position on the target element. 
c
            if (cneutf.eq.0) then 
               X0 = RP(ID)
               Y0 = ZP(ID)
            elseif (cneutf.eq.1) then 
c
c              Find launch position along target element. 
c               
               itmp = wallindex(id)
c
c               write(6,'(a,2i6,5(1x,g12.5))') 'TL1:',id,itmp,rp(id),
c     >                    wallpt(itmp,1),zp(id),wallpt(itmp,2)
c
               CALL SURAND2 (SEED, 1, RAN)
               NRAND = NRAND + 1
c
c              All wall segment halves are equal for wall/target elements 
c
               IF (RAN.LT.0.5) then 
                  CALL SURAND2 (SEED, 1, RAN)
                  NRAND = NRAND + 1
                  X0 = WALLPT(itmp,1) + RAN * WALLPT(itmp,5)
     >                * COS(wallpt(itmp,8))
                  Y0 = WALLPT(itmp,2) + RAN * WALLPT(itmp,5)
     >                * SIN(wallpt(itmp,8))
                  ISPRODS(IPROD+Nneut1) = 0
               ELSE
                  CALL SURAND2 (SEED, 1, RAN)
                  NRAND = NRAND + 1
                  X0 = WALLPT(itmp,1) + RAN * WALLPT(itmp,6)
     >                * COS(WALLPT(itmp,9))
                  Y0 = WALLPT(itmp,2) + RAN * WALLPT(itmp,6)
     >                * SIN(WALLPT(itmp,9))
                  ISPRODS(IPROD+Nneut1) = 1
               ENDIF

c              write(6,'(a,2i6,8(1x,g12.5))') 'TL2:',id,itmp,x0,y0,
c     >           wallpt(itmp,20),wallpt(itmp,21),wallpt(itmp,22),
c     >           wallpt(itmp,23)

           endif


c
            if (yieldsw.eq.0) then 
               RMAXS (IPROD+nneut1) = KRMAX(ID)
            elseif (yieldsw.eq.1) then 
               RMAXS (IPROD+nneut1) = 1.0
            endif
C
         ELSEIF (CNEUTB.EQ.1) THEN
c
            ID = 1
            X0 = CXSC
            Y0 = CYSC
            rmaxs(iprod+nneut1) = 1.0
C
         ELSEIF (CNEUTB.EQ.2.or.cneutb.eq.4) THEN
            CALL SURAND2 (SEED, 1, RAN)
            NRAND = NRAND + 1
            ID = IPOS(RAN,fyprob,nfymap)
            id = fymap(id)
c
c           Do not limit particle energies for chemical sputtering 
c
c slmod begin
            if (yieldsw.eq.0.and.pinsw.ne.4) then  
c
c            if (yieldsw.eq.0) then  
c slmod end
               RMAXS (IPROD+nneut1) = krmaxw(id)
            else
               RMAXS (IPROD+nneut1) = 1.0
            endif
c
            CALL SURAND2 (SEED, 1, RAN)
            NRAND = NRAND + 1
c
c           All wall segment halves are equal for NIMBUS walls. 
c
c           IF (RAN.LT.(WALLPT(id,10)/WALLPT(id,12))) then 
c
            IF (RAN.LT.0.5) then 
               CALL SURAND2 (SEED, 1, RAN)
               NRAND = NRAND + 1
               X0 = WALLPT(id,1) + RAN * WALLPT(id,5)
     >             * COS(wallpt(id,8))
               Y0 = WALLPT(id,2) + RAN * WALLPT(id,5)
     >             * SIN(wallpt(id,8))
               ISPRODS(IPROD+Nneut1) = 0
            ELSE
               CALL SURAND2 (SEED, 1, RAN)
               NRAND = NRAND + 1
               X0 = WALLPT(id,1) + RAN * WALLPT(id,6)
     >             * COS(WALLPT(id,9))
               Y0 = WALLPT(id,2) + RAN * WALLPT(id,6)
     >             * SIN(WALLPT(id,9))
               ISPRODS(IPROD+Nneut1) = 1
            ENDIF
C
         elseif (cneutb.eq.5) then 
c
            ID = IPOS (RAN, neut2d_prob, neut2d_num)
            IK = neut2d_index(id,1)
            IR = neut2d_index(id,2)
c
            X0 = RS(ik,ir)
            Y0 = ZS(ik,ir)
c
            RMAXS (IPROD+nneut1) = 1.0
            ID = IK
c
c        Line 
c
         ELSEIF (CNEUTB.EQ.6) THEN
c
            ID = 1
c
            NRAND = NRAND + 1
            CALL SURAND2 (SEED, 1, RAN)
            X0 = CXSCA + RAN * (CXSCB - CXSCA)
            Y0 = CYSCA + RAN * (CYSCB - CYSCA)
c
            rmaxs(iprod+nneut1) = 1.0
c
c        Box
c
         ELSEIF (CNEUTB.EQ.7) THEN
c
            ID = 1
c
            NRAND = NRAND + 1
            CALL SURAND2 (SEED, 1, RAN)
            X0 = CXSCA + RAN * (CXSCB - CXSCA)
c
            NRAND = NRAND + 1
            CALL SURAND2 (SEED, 1, RAN)
            Y0 = CYSCA + RAN * (CYSCB - CYSCA)
c
            rmaxs(iprod+nneut1) = 1.0
C
         ENDIF
C
C
C----- STORE LAUNCH DETAILS IN PRODUCTION ARRAYS
C
         XPRODS(IPROD+Nneut1) = X0
         YPRODS(IPROD+Nneut1) = Y0
         eprods(iprod+nneut1) = 0.0 
         SPUTYS(IPROD+Nneut1) = 1.0
         IDPRODS(IPROD+Nneut1) = ID


c         write(6,'(a,2i8,10(1x,g18.8))') 'Launch data:',iprod,nneut1,
c     >         iprod+nneut1,xprods(IPROD+Nneut1),yprods(IPROD+Nneut1),
c     >         eprods(iprod+nneut1),idprods(iprod+nneut1)


      end do
C
C-----------------------------------------------------------------------
C     LAUNCH AND FOLLOW PRIMARY NEUTRALS
C-----------------------------------------------------------------------
C
c     SET variable defining the type of the launch -
c
c     i.e. Wall or Target - Physical or Chemically sputtered.
c      
      neuttype =0 
c
c     Free Space= 0 - most data not recorded - cneutb=1
c     Targ+phys = 1
c     Targ+chem = 2
c     Targ-self = 3
c     Wall+phys = 4
c     Wall+chem = 5
c     2Dneutral = 6
c     Refl Ion  = 7
c
      if (cneutb.eq.0.or.cneutb.eq.3) then 
         if (yieldsw.eq.0) then 
            neuttype = 1  
         elseif (yieldsw.eq.1) then 
            neuttype = 2  
         endif
      elseif (cneutb.eq.2.or.cneutb.eq.4) then 
         if (yieldsw.eq.0) then 
            neuttype = 4  
         elseif (yieldsw.eq.1) then 
            neuttype = 5  
         endif
      elseif (cneutb.eq.5) then 
         neuttype = 6     
c slmod begin
      elseif (grdnmod.ne.0.and.cneutb.eq.1) then 
c...    Is this correct, or was CNEUTB=1 just overlooked? 
        WRITE(0,*) 'WARNING: NOT SURE THIS IS THE CORRECT '//
     .             'NEUTTYPE FOR CNEUTB=1'
        neuttype = 0
c slmod end
      endif   
c
      RATIZA  = 0.0
      RNEUTA  = 0.0
      RWALLNA = 0.0
      MTCWALLNA = 0.0
      RCENTA  = 0.0
      RTMAXA  = 0.0
      RSTRUKA = 0.0
      MTCSTRUKA = 0.0 
      RFAILA  = 0.0
      RMAINA  = 0.0
      REXITA  = 0.0
C
         !
         ! Initialize the velocity distribution debugging code - if required storage will be allocated
         !
      call init_velocity_dist(debug_neutv,debug_neutv_einmax,
     >                        debug_neutv_nbins,crmi)


! ammod begin.

      ! Insertion point for hydrocarbon following code.
      ! hc_follow_option is found in ComHC.
      ! Subroutine HC_Launch is found in hc_start.d6b. global_hc_launch is in hc_global_routines.d6a

      If ((global_hc_follow_option .eq. 1.and.
     >    (yieldsw .eq. 1.or.cneutb.eq.1))) Then ! Hydrocarbon following is activated.


	   ! Call HC to follow neutral HC fragments to eventual evolution to C+ ion or capture.
           ! Note there should be no differences from call to LAUNCH.
c
           Call global_HC_Launch (NNEUT1+1,NNEUT1+NNEUT2,NION1+1,NION2,
     >     RSTRUKA,MTCSTRUKA,RMAINA,REXITA,
     >     RATIZA,RNEUTA,RWALLNA,MTCWALLNA,RCENTA,RTMAXA,
     >     SEED,NRAND,NEUTIM,RFAILA,STATUS,MATP,MATT,
     >     neuttype,cneutb,cneutc)

!	   Write (0,*) "DIVIMP Launched chemsput:",
!     >     RSTRUKA,MTCSTRUKA,RMAINA,REXITA,RATIZA,
!     >     RNEUTA,RWALLNA,MTCWALLNA,RCENTA,RTMAXA  

      ! Additional options may be added here at a later date.   
      
      Else ! Physical sputtering if yieldsw = 0, or HC not activated.
           ! Call LAUNCH to follow neutrals from birth to eventual ionization or capture.

!           Write (0,*) "DIVIMP Physsput starting",NNEUT1+1,NNEUT1+
!     >       NNEUT2,NION1+1,NION2,"neutt",neuttype,"launch",cneutb,
!     >       "v/a",cneutc

	   CALL LAUNCH (NNEUT1+1,NNEUT1+NNEUT2,NION1+1,NION2,
     >     RSTRUKA,MTCSTRUKA,RMAINA,REXITA,
     >     RATIZA,RNEUTA,RWALLNA,MTCWALLNA,RCENTA,RTMAXA,
     >     SEED,NRAND,NEUTIM,RFAILA,STATUS,MATP,MATT,neuttype)

!	   Write (0,*) "DIVIMP Launched physsput:",
!     >     RSTRUKA,MTCSTRUKA,RMAINA,REXITA,RATIZA,
!     >     RNEUTA,RWALLNA,MTCWALLNA,RCENTA,RTMAXA  

      End If
! ammod end.



c
c     Print out Neutral velocity debug information and 
c     deallocate the arrays if the neutral velocity distribution debugging code
c     was active
c
c     Code has been moved to the velocity_dist module
c
      call print_vdist(crmi)


c
c Original call to launch
c
c      CALL LAUNCH (NNEUT1+1,NNEUT1+NNEUT2,NION1+1,NION2,RSTRUKA,
c     >             MTCSTRUKA,RMAINA,REXITA,
c     >             RATIZA,RNEUTA,RWALLNA,MTCWALLNA,RCENTA,RTMAXA,
c     >             SEED,NRAND,
c     >             NEUTIM,RFAILA,STATUS,MATP,MATT,neuttype)
c


c
c     Reset global variables 
c 
      CNEUTB = TMPCNEUTB
      CNEUTC = TMPCNEUTC
c
      return 
      end
c
c
c
      subroutine printfy(fydata,fymap,fyprob,nfy,nfymap,totfydata)
      implicit none
      include 'params'
c
      real fydata(maxpts,5),fyprob(maxpts)
      real totfydata(3,5)
      integer nfy,nfymap,fymap(maxpts)
c
      include 'cgeom'
      include 'comtor'
      include 'pindata'
      include 'cneut'
      include 'cneut2'
      include 'printopt' 
c
c     Local variables  
c
      integer id,ik,ir,in
      real    launchtype,sputtertype,ymf
c
c     Now in printopt common block - values set just after grids read 
c     in - in the tau module
c
c      character*5 inner,outer
c
      launchtype = totfydata(3,2)
      sputtertype = totfydata(3,4)
c
      call prb 
c
      CALL PRC ('SAMPLE PRIMARY FLUX AND YIELD DATA')
c
      if (cneutb.eq.1.or.cneutb.eq.6.or.cneutb.eq.7) then
       call prb
       call prc('FREE-SPACE LAUNCH: FLUX AND YIELD DATA ARE')
       call prc('PRINTED BUT NOT UTILIZED FOR THE LAUNCH')
       call prb
      endif
c
      if (sputtertype.eq.0.0) then 
         call prc('    - FOR A PHYSICALLY SPUTTERED') 
      elseif (sputtertype.eq.1.0) then 
         call prc('    - FOR A CHEMICALLY SPUTTERED') 
      endif      
c
      if (launchtype.eq.0.0) then 
         call prc('      TARGET LAUNCH (DUE TO ION FLUX)') 
      elseif (launchtype.eq.1.0) then 
         call prc('      WALL AND TARGET LAUNCH (DUE TO ATOM FLUX)') 
      endif      
c
c     Print target launch data 
c
      if (launchtype.eq.0.0) then 
         WRITE (7,9000)
         DO ID = 1, nfy
c
            if (sputtertype.eq.0) then 
               ymf = kmfps(id)
            elseif(sputtertype.eq.1) then 
               ymf = kmfcs(id)
            endif
c
            WRITE (7,9002) id,
     >         RP(ID),ZP(ID),ymf,fydata(id,1),fydata(id,2),
     >         fydata(id,4),fydata(id,5),kbfs(ikds(id),irds(id)),
     >         dds(id),costet(ID),tempds(id),
c slmod begin
     >         kteds(id),ktids(id),nimindex(id)
c slmod end
         end do
C
         call prb
c
         CALL PRR('PRIMARY INTEGRATED FLUX       '//INNER,
     >               totfydata(1,1))
         CALL PRR('PRIMARY INTEGRATED FLUX*YIELD '//INNER,
     >               totfydata(1,5))
         CALL PRR('THERMAL+KINETIC    HEAT FLUX  '//INNER,
     >               totfydata(1,3)*1.6E-19)
         CALL PRR('THERM+KINETIC+POT.EN HEAT FLX '//INNER,
     >        totfydata(1,3)*ech+totfydata(1,1)*15.8*ech)
c
         CALL PRR('PRIMARY INTEGRATED FLUX       '//OUTER,
     >               totfydata(2,1))
         CALL PRR('PRIMARY INTEGRATED FLUX*YIELD '//OUTER,
     >               totfydata(2,5))
         CALL PRR('THERMAL+KINETIC    HEAT FLUX  '//OUTER,
     >               totfydata(2,3)*ech)
         CALL PRR('THERM+KINETIC+POT.EN HEAT FLX '//OUTER,
     >        totfydata(2,3)*ech+totfydata(2,1)*15.8*ech)
c
         CALL PRR('TOTAL PRIMARY INTEGRATED FLUX      ',
     >               totfydata(3,1))
         CALL PRR('TOTAL PRIMARY INTEGRATED FLUX*YIELD',
     >               totfydata(3,5))
         CALL PRR('TOTAL THERMAL+KINETIC HEAT FLUX    ',
     >               totfydata(3,3)*ech)
         CALL PRR('TOTAL THERM+KINETIC+POT.EN HEAT FLX',
     >        totfydata(3,3)*ech+totfydata(3,1)*15.8*ech)
c
c
c         if (sputtertype.eq.0) then 
c            call prr('Yield Multiplication Factor Used:',
c     >               cymfs(1,2))
c         elseif (sputtertype.eq.1) then 
c            call prr('Yield Multiplication Factor Used:',
c     >               cymfs(1,4))
c         endif
c
      elseif (launchtype.eq.1.0) then 
c
         WRITE (7,9001)
         DO ID = 1, nfy
c
            if (sputtertype.eq.0) then 
               ymf = kmfpws(id)
            elseif(sputtertype.eq.1) then 
               ymf = kmfcws(id)
            endif
c
            WRITE (7,9003) id,
     >         wallpt(ID,1),wallpt(id,2),ymf,fydata(id,1),
     >         fydata(id,2),
     >         fydata(id,4),fydata(id,5),
     >         wallpt(id,7),wallpt(id,16),wallpt(id,19)
         end do
c
         call prb
c 
         CALL PRR('TOTAL PRIMARY INTEGRATED FLUX      ',
     >               totfydata(3,1))
         CALL PRR('TOTAL PRIMARY INTEGRATED FLUX*YIELD',
     >               totfydata(3,5))
         CALL PRR('TOTAL PRIMARY INTEGRATED HEAT FLUX ',
     >               totfydata(3,3)*1.6E-19)
c
c         if (sputtertype.eq.0) then 
c            call prc('Yield Multiplication Factors Used:')
c            call prr('   - for TARGET segments : ',cymfs(1,2))
c            call prr('   - for WALL segments   : ',cymfs(1,5))
c         elseif (sputtertype.eq.1) then 
c            call prc('Yield Multiplication Factors Used:')
c            call prr('   - for TARGET segments : ',cymfs(1,4))
c            call prr('   - for WALL segments   : ',cymfs(1,6))
c         endif
c
      endif 
c
      CALL PRB
c
      CALL PRC ('CUMULATIVE PROBABILITY DATA FOR THIS LAUNCH')
      WRITE (7,9009)
      DO ID = 1,nfymap
         if (id.eq.1) then
            WRITE(7,9008) id,fymap(id),fyprob(id),fyprob(id)
         else 
            WRITE(7,9008) id,fymap(id),fyprob(id)-fyprob(id-1),
     >                    fyprob(id)
         endif
      end do
c

 9000 FORMAT('IND   R    Z    YMF  FLUXDENS  ENERGY YIELD    F*Y  ',
     >  '  Bt/Bth  LENGTH   ORTH TEMP')
 9001 FORMAT('IND   R    Z    YMF  FLUXDENS  ENERGY YIELD    F*Y  ',
     >  '          LENGTH   TYPE TEMP')
 9002 FORMAT(I3,F5.2,1x,f5.2,F5.2,1P,1x,G9.2,0P,F7.2,1P,
     >      2G9.2,0P,1x,F5.2,1P,
c slmod begin
     >       G9.2,0P,1x,f5.2,f5.0,2x,2f5.1,i4)
c
c     >       G9.2,0P,1x,f5.2,f5.0)
c slmod end
 9003 FORMAT(I3,F5.2,1x,f5.2,F5.2,1P,1x,G9.2,0P,F7.2,1P,
     >      2G9.2,0P,1x,5x,1P,
     >       G9.2,0P,1x,f5.1,f5.0)
c
c 9002 FORMAT(I3,2F5.2,F5.2,1P,G9.2,0P,F7.2,1P,2G9.2,0P,1x,F5.2,1P,
c     >       1x,G9.2,0P,1x,f5.2,f5.0)
c 9003 FORMAT(I3,2F5.2,F5.2,1P,G9.2,0P,F7.2,1P,2G9.2,0P,1x,5x,1P,
c     >       1x,G9.2,0P,1x,f5.1,f5.0)
c
 9008 FORMAT(5X,I4,5X,I4,2x,G12.5,2X,G12.5)
 9009 FORMAT(4x,'INDEX',3x,'SEG.REF.',3x,'PROBABILITY',3x,'ACCUM.PROB.')
c
      return
      end
c
c
c
      subroutine prep_neut2d
      implicit none
c
      include    'params'
      include    'comtor'
      include    'cgeom'
      include    'cneut'
      include    'cneut2'
c
c
c     PREP_NEUT2D:
c
c     This subroutine prepares the arrays that are used for the
c     probability and mapping of the 2D neutral launch option. It 
c     uses the array neut2d_raw which is expected to contain 
c     neutral particles source rates in units of [#/m-toroidally/s].
c
c     It returns the total source strength and set up the ancilliary 
c     data arrays that are required to perform the launch. 
c
      integer in,ik,ir
      real    neut2d_tmp
c
      neut2d_src = 0.0
      neut2d_num = 0
c
c     Calculate TOTAL source strength - needed for probability scaling 
c
      do ir = 1,nrs
c     
         do ik = 1,nks(ir)
c
            if (neut2d_raw(ik,ir).gt.0.0) then 
c
               neut2d_src = neut2d_src + neut2d_raw(ik,ir)
c
            endif 
c  
         end do
c       
      end do 
c 
c     Build integrated probability array. 
c
      in = 0 
c
      do ir = 1,nrs
c
         do ik = 1,nks(ir)
c
            if (neut2d_raw(ik,ir).gt.0.0) then 
c
c              Increment counting index
c
               in = in + 1
c
c              Record cell indices
c
               neut2d_index(in,1) = ik 
               neut2d_index(in,2) = ir 
c
c              Assign cumulative probability 
c
               if (in.eq.1) then
c
                  neut2d_prob(in) = neut2d_raw(ik,ir)/neut2d_src
c
               else
c     
                  neut2d_prob(in) = neut2d_prob(in-1) + 
     >                             neut2d_raw(ik,ir)/neut2d_src
               endif 
c
            endif
c
         end do
c 
      end do
c 
      neut2d_num = in
c
c     Last entry should alread be 1.0 - but force it just in case of 
c     rounding errors.
c
      neut2d_prob(neut2d_num) = 1.0 
c
      return
      end 
c
c
c
      subroutine print_neut2d
      implicit none
c
      include    'params'
      include    'comtor'
      include    'cgeom'
      include    'cneut'
      include    'cneut2'
c 
c     PRINT_NEUT2D:
c
c     This routine prints out the 2D probability distribution and 
c     total influx for the 2D particle source.
c
      real prob_limit,cell_prob
      parameter (prob_limit=1.0e-5)
      integer in,prcount
      character*80 comment
c
      call prb  
      call prc('SUMMARY OF 2D NEUTRAL SOURCE LAUNCH PROBABILITY')
      call prr('- LISTING INCLUDES ALL CELLS WITH A PROBABILITY > ',
     >           prob_limit)
      call pri('- NUMBER OF CELLS WITH NON-ZERO LAUNCH PROBABILITY = ',
     >          neut2d_num)
      call prc('- LISTING IS SETS OF (IK,IR) PAIRS FOLLOWED'//
     >         ' BY PROBABILITY')
c
      call prb
      call prr('TOTAL STRENGTH OF 2D SOURCE #/m-tor/s = ',neut2d_src)  
      call prb
c
      prcount = 0    
C
      do in = 1,neut2d_num
c
         if (in.eq.1) then 
            cell_prob = neut2d_prob(in)
         else
            cell_prob = neut2d_prob(in) - neut2d_prob(in-1)
         endif
c
         if (cell_prob.gt.prob_limit) then 
c
c           Increment print count
c
            prcount = prcount + 1
c
c           write items to line
c
            write (comment((prcount-1)*20+1:prcount*20),100) 
     >           neut2d_index(in,1),neut2d_index(in,2),cell_prob
c
c           Print 4 entries / line
c
            if (prcount.eq.4) then 
c
               call prc(comment)
               prcount = 0
c
            endif
c
         endif 
c
      end do
c
      call prb
c

 100  format('(',i3,',',i3,')=',1p,g9.2,1x)
 
c
      return
      end 
c
c
c
      subroutine redistribute_nprod(nproda,nprod,nprod2a,nprod2,
     >                            nprod_neut2d,pinsw,matt,matp)
      implicit none
      integer nproda,nprod,nprod2a,nprod2,pinsw,matt,matp
      include    'params'
      include    'cyield'
      include    'comtor'
      include    'cgeom'
      include    'cneut'
      include    'cneut2'
c
c     REDISTRIBUTE_NPROD:
c
c     This routine analyses the requested neutral launch options and 
c     the corresponding sputter options and redistributes the numbers
c     of particles for each launch in an appropriate fashion. 
c
      integer   nprod_neut2d
      real      ntot 
      real      w1a,w1,w2a,w2,w3,wtot 
      real      fydata(maxpts,5),fyprob(maxpts)
      real      totfydata(3,5)
      integer   nfy,nfymap,yieldsw,newcneutc,newcneutb
      integer   fymap(maxpts)
c slmod begin
      integer   in
c slmod end
c      
      ntot = nprod + nprod2
c
      nproda = 0
      nprod2a = 0
      nprod_neut2d = 0      
c
      w1a = 0.0
      w1  = 0.0
      w2a = 0.0
      w2  = 0.0
      w3  = 0.0
      wtot= 0.0 
      IF (sloutput) WRITE(0,*) 'REDIST:',nprod,nprod2 
      IF (sloutput) WRITE(0,*) 'FIRST:',cneutb,cneutd
c
c     Calculate total and partial FY for first launch 
c
      if (nprod.gt.0) then  
c
         if (cneutb.eq.0.or.cneutb.eq.3) then 
c
            if (cneutd.eq.0.or.cneutd.eq.1.or.
     >          cneutd.eq.3.or.cneutd.eq.4.or.
     >          cneutd.eq.7.or.cneutd.eq.8) then 
c
c              Physical Sputtering On targets 
c     
               if (cneutd.eq.7) then 
                  call tfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  3,0,matp,matt)               
               else
                  IF (sloutput) WRITE(0,*) 'DEBUG: CALL TFY I'
                  call tfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  pinsw,0,matp,matt)               
               endif
c
               w1a = 0.0
               w1 = totfydata(3,5) 
c                              
            elseif (cneutd.eq.2.or.cneutd.eq.5) then 
c
c              Chemical Sputtering On targets 
c
               call tfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  pinsw,1,matp,matt)               
               w1a = 0.0
               w1 = totfydata(3,5) 
c
            elseif (cneutd.eq.6) then 
c
c              Physical Sputtering On targets 
c     
               call tfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  pinsw,0,matp,matt)               
               
               w1a = totfydata(3,5) 
c                              
c              Chemical Sputtering On targets 
c
               call tfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  pinsw,1,matp,matt)               
               w1 = totfydata(3,5) 
c
            endif
c
         elseif (cneutb.eq.1.or.cneutb.eq.2.or.
     >           cneutb.eq.6.or.cneutb.eq.7) then 
c
c           Use a weight value of -1 to indicate that NPROD should not be 
c           redistributed.
c
            w1a = -1.0
            w1  = -1.0
c
         elseif (cneutb.eq.4) then   
c
            if (cneutd.eq.0.or.cneutd.eq.1.or.
     >          cneutd.eq.3.or.cneutd.eq.4.or.cneutd.eq.8) then 
c
c              Physical Sputtering On Walls 
c
               if (cneutd.eq.8) then  
                  call wfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  2,0,matp,matt)               
               else
                  IF (sloutput) WRITE(0,*) 'DEBUG: CALL WFY G'
                  call wfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  pinsw,0,matp,matt)               
               endif
c 
               w1a = 0.0
               w1 = totfydata(3,5) 
c
            elseif (cneutd.eq.2.or.cneutd.eq.5) then 
c
c              Chemical Sputtering On Walls 
c
               IF (sloutput) WRITE(0,*) 'DEBUG: CALL WFY H'
               call wfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  pinsw,1,matp,matt)               
                w1a = 0.0
                w1 = totfydata(3,5) 
c
c
            elseif (cneutd.eq.6) then 
c
c               Physical Sputtering On Walls 
c
                IF (sloutput) WRITE(0,*) 'DEBUG: CALL WFY I',pinsw
                call wfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                   pinsw,0,matp,matt)               
                w1a = totfydata(3,5) 
c
c               Chemical Sputtering On Walls 
c
                IF (sloutput) WRITE(0,*) 'DEBUG: CALL WFY J',pinsw
                call wfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                   pinsw,1,matp,matt)               
                w1 = totfydata(3,5) 
c

            endif
c
         elseif (cneutb.eq.5) then 
c
            w1a = 0.0
            w1  = neut2d_src
c
c
         endif
c
      endif
c
c     Deal with supplementary launch
c
c
c     Calculate total and partial FY for second launch 
c
      IF (sloutput) WRITE(0,*) 'SECND:',cneuth,cneutd2
      if (nprod2.gt.0) then  
c
         if (cneuth.eq.0.or.cneuth.eq.3) then 
c
            if (cneutd2.eq.0.or.cneutd2.eq.1.or.
     >          cneutd2.eq.3.or.cneutd2.eq.4.or.
     >          cneutd2.eq.7.or.cneutd2.eq.8) then 
c
c              Physical Sputtering On targets 
c     
               if (cneutd2.eq.7) then 
                  call tfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  pinsw,0,matp,matt)               
               else
                  IF (sloutput) WRITE(0,*) 'DEBUG: CALL TFY J'
                  call tfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  pinsw,0,matp,matt)               
               endif
c 
               w2a = 0.0
               w2 = totfydata(3,5) 
c                              
            elseif (cneutd2.eq.2.or.cneutd2.eq.5) then 
c
c              Chemical Sputtering On targets 
c
               call tfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  pinsw,1,matp,matt)               
               w2a = 0.0
               w2 = totfydata(3,5) 
c
            elseif (cneutd2.eq.6) then 
c
c              Physical Sputtering On targets 
c     
               call tfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  pinsw,0,matp,matt)               
               
               w2a = totfydata(3,5) 
c                              
c              Chemical Sputtering On targets 
c
               call tfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  pinsw,1,matp,matt)               
               w2 = totfydata(3,5) 
c
            endif
c
         elseif (cneuth.eq.1.or.cneuth.eq.2) then 
c
c           Use a weight value of -1 to indicate that NPROD2 should not be 
c           redistributed.
c
            w2a = -1.0
            w2 = -1.0
c
         elseif (cneuth.eq.4) then   
c
c
            if (cneutd2.eq.0.or.cneutd2.eq.1.or.
     >          cneutd2.eq.3.or.cneutd2.eq.4.or.
     >          cneutd2.eq.8) then 
c
c              Physical Sputtering On Walls 
c
               if (cneutd2.eq.8) then 
                  call wfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  pinsw,0,matp,matt)               
               else
                  IF (sloutput) WRITE(0,*) 'DEBUG: CALL WFY K'
                  call wfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  pinsw,0,matp,matt)               
               endif

               w2a = 0.0
               w2 = totfydata(3,5) 
c
            elseif (cneutd2.eq.2.or.cneutd2.eq.5) then 
c
c              Chemical Sputtering On Walls 
c
               IF (sloutput) WRITE(0,*) 'DEBUG: CALL WFY L'
               call wfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                  pinsw,1,matp,matt)               
                w2a = 0.0
                w2 = totfydata(3,5) 
c
c
            elseif (cneutd2.eq.6) then 
c
c               Physical Sputtering On Walls 
c
                IF (sloutput) WRITE(0,*) 'DEBUG: CALL WFY M'
                call wfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                   pinsw,0,matp,matt)               
                w2a = totfydata(3,5) 
c
c               Chemical Sputtering On Walls 
c
                IF (sloutput) WRITE(0,*) 'DEBUG: CALL WFY N'
                call wfy(fydata,fymap,fyprob,nfy,nfymap,totfydata,
     >                   pinsw,1,matp,matt)               
                w2 = totfydata(3,5) 
c

            endif
c
         elseif (cneuth.eq.5) then 
c
            w2a = 0.0
            w2  = neut2d_src
c
         endif
c
      endif
c
c     Calculate the total FY for additional 2D launch option. 
c
      if (neut2d_opt.eq.0) then 
c  
         w3 = 0.0
c
      elseif (neut2d_opt.eq.1) then 
c
         w3 = neut2d_src
c
      endif  
c
c     FY weights for all allowed mechanisms have been assigned. 
c     Now redistribute the number of particles.       
c
c     Deal with the FY scalable sources first
c
      if (w1.ne.-1.0.and.w2.ne.-1.0) then  
c
         wtot = w1a + w1 + w2a + w2 + w3 
         ntot = nprod + nprod2      
c
         if (wtot.gt.0.0) then 

            nproda =  nint ( w1a/wtot * ntot)
            nprod  =  nint ( w1 /wtot * ntot)
            nprod2a=  nint ( w2a/wtot * ntot)
            nprod2 =  nint ( w2 /wtot * ntot)
c
            nprod_neut2d =  nint ( w3 /wtot * ntot)
c
         else
c
c           Issue error message - assign all particles to first 
c           primary launch mechanism 
c
            nproda = 0.0
            nprod  = 0.0
            nprod2a=0.0 
            nprod2 = 0.0
c
            write(6,*) 'ERROR: Total Weight for Neutral launch is zero'
            write(0,*) 'ERROR: Total Weight for Neutral launch is zero'
c
         endif 
c
c     Non-scalable first launch plus scalable second   
c
      elseif (w1.eq.-1.and.w2.ne.-1.0) then  
c
c        Nprod is left unchanged
c         
         wtot = w2a + w2 + w3 
         ntot = nprod2      
c
         if (ntot.gt.0.0.and.wtot.gt.0.0) then 
c
            nprod2a=  nint ( w2a/wtot * ntot)
            nprod2 =  nint ( w2 /wtot * ntot)
c
            nprod_neut2d =  nint ( w3 /wtot * ntot)
c
         else
c
c           Zero out nprod2,nprod2a since no particles launched 
c           from this method
c
            nprod2a      = 0
            nprod2       = 0
            nprod_neut2d = 0 
c
         endif
c
c
c     Non-scalable second launch plus scalable first  
c
      elseif (w1.ne.-1.and.w2.eq.-1.0) then  
c
c        Nprod2 is left unchanged
c         
         wtot = w1a + w1 + w3 
         ntot = nprod      
c
         if (ntot.gt.0.0.and.wtot.gt.0.0) then
c
            nproda =  nint ( w1a/wtot * ntot)
            nprod  =  nint ( w1 /wtot * ntot)
c
            nprod_neut2d =  nint ( w3 /wtot * ntot)
c
         else 
c
c           Zero out nprod,nproda since no particles launched 
c           from this method
c
            nproda       = 0
            nprod        = 0
            nprod_neut2d = 0
c
         endif  
c
c     Issue an error message if a 2D neutral launch has been
c     specified in combination with NON-scalable launch options.   
c
c     nprod values are left unchanged. 
c
      elseif (w1.eq.-1.and.w2.eq.-1.and.w3.ne.0.0) then 
c
         call prc('ERROR: 2D neutral launch called for BUT no')
         call prc('       other Flux/Yield data is available')
         call prc('       for scaling the particles to be launched.')
         call prc('       This is caused by choosing launch options')
         call prc('       which have NO flux/yield data available')
         call prc('       in combination with the 2D NUETRAL source.') 
 
 
      endif 
c
c     Print information
c
      write (6,'(a,6i8)') 'Nprods:',nproda,nprod,nprod2a,
     >                     nprod2,nprod_neut2d,nint(ntot)
      write (6,'(a,6g12.5)') 'Weights:',w1a,w1,w2a,w2,w3
c
c     Exit                 
c
      return
      end
c
c
c
      real function find_thompson_velocity(neuttype,id,seed,nrand,
     >                                     thom_opt)
      use error_handling
      implicit none 
      integer neuttype,id,nrand,thom_opt
      real*8 seed
c
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'pindata'
c
      common /thom_ye_params/  eimp,gamma,ebd
      real*8 eimp,gamma,ebd 

      common /thom_yv_params/  vimp,vgamma,vbd
      real*8 vimp,vgamma,vbd 

c
c     -Need a random number
c     -Need to calculate the integral over the distribution for this 
c      segment or have it stored in a pre-calculated array. 
c     -Option used for physical sputtering only.      
c     -What about the difference between atom and ion flux physical 
c      sputtering energies? These are part of Eimpact
c     -What is the value of gamma? gamma = 4 (mC*mD) / (mC+mD)**2
c     
c      intye_targ(maxnds)
c      intye_wall(maxpts)
c
c
c     Note: for now use the same formula for self sputtered particles as for
c           regular sputtering
c
c
c     Local variables  
c
      integer in,ik,ir,targ,ierr,iter_cnt
      real*8 sheath_te,sheath_ti
      real*8 emin,emax,e_result
      real*8 vmin,vmax,vtmp,v_result
c
      integer ringno
      real*8 thom_ye,thom_yv
      external thom_ye,thom_yv,ringno
c
c     Calculate Eimp and Vimp
c
c     Allow for self-sputtering to fall through here for now
c
      if ((neuttype.eq.1.or.neuttype.eq.3).and.
     >     id.ge.1.and.id.le.nds) then 

         IK = IKDS(ID)
         IR = IRDS(ID)
c        
c        Set target identifier - TARG = 1 = OUTER for SONNET
c                                          INNER for GRID2D   
c                               TARG = 2 = INNER for SONNET
c                                          OUTER for GRID2D   
         if (ik.gt.nks(ir)/2) then 
            targ = 1
         else 
            targ = 2
         endif 
c       
c        Assign Sheath quantities
c       
         if ((targ.eq.1.and.nsheath_vali.eq.0).or.
     >       (targ.eq.2.and.nsheath_valo.eq.0)) then 
            sheath_te = kteds(id)
            sheath_ti = ktids(id)
         else
c       
c           The code assumes for now that Te=Ti when the
c           sheath potential is directly specified. 
c       
            if (targ.eq.1) then 
        
               in = RINGNO(ir,sheath_vali,nsheath_vali,
     >                     maxnrs,2,ierr)
               if (ierr.eq.0) then 
                  sheath_te = sheath_vali(in,2)
                  sheath_ti = sheath_vali(in,2)
               else
                  sheath_te = kteds(id)
                  sheath_ti = ktids(id)
               endif              
        
            elseif (targ.eq.2) then 
        
               in = RINGNO(ir,sheath_valo,nsheath_valo,
     >                     maxnrs,2,ierr)
               if (ierr.eq.0) then 
                  sheath_te = sheath_valo(in,2)
                  sheath_ti = sheath_valo(in,2)
               else
                  sheath_te = kteds(id)
                  sheath_ti = ktids(id)
               endif              
c       
           endif
c           
         endif
c
         eimp = 2.0 * sheath_ti + 3.0 * RIZB * sheath_te
c
      elseif (neuttype.eq.4.and.id.ge.1.and.id.le.wallpts) then 
c
         eimp = flxhw5(int(wallpt(id,17)))          
c
      endif
c
      vimp = 1.38e4 * sqrt(eimp/crmi)    
c
c     Calculate Gamma and vgamma
c
      GAMMA  = 4.0 * CRMB * CRMI / ((CRMB+CRMI) * (CRMB+CRMI))
      vgamma = gamma
c
c     Calculate Ebd and vbd
c
      ebd = cebd
c
      vbd = 1.38e4 * sqrt(ebd/crmi)    
c
c     Calculate emin and emax limits for integration
c     Also vmin and vmax limits 
c     - note vmax requires a 
c       square root and instead of getting an error when negative
c       it would be imaginary - thus the sqrt is taken to get vmax
c       only after the quantity is verified. 
c
c
      emin = 0.0
      vmin = 0.0  
      emax = gamma*(1.0d0-gamma)*eimp-ebd
      vtmp = gamma*(1.0d0-gamma)*vimp**2-vbd**2
c
c     Check validity of Emax 
c
      if (thom_opt.eq.0.and.emax.lt.0.0) then 
c
         call errmsg('FIND_THOMSON_VELOCITY','EMAX<0')
c
         write(6,'(a,4(1x,g12.5))')
     >                  'ERROR : FIND_THOMPSON_VELOCITY : EMAX<0',
     >                  eimp,ebd,gamma,emax
c
c         stop 
c
c        Add a temporary correction where Emax is set to 
c        Ebd + 1% of Ebd - this will give a slow moving particle
c        and return some velocity for the particle.
c          
         eimp = 1.01  * ebd / ( gamma*(1.0d0-gamma))
         emax = gamma*(1.0d0-gamma)*eimp-ebd

      endif
c
c     Check validity of vmax
c
      if (thom_opt.eq.1.and.vtmp.lt.0.0) then 
c
         write(6,'(a,4(1x,g12.5))')
     >                  'ERROR : FIND_THOMPSON_VELOCITY :'//
     >                  ' VMAX IMAGINARY',
     >                  vimp,vbd,gamma,vtmp
         write(0,'(a,4(1x,g12.5))')
     >                  'ERROR : FIND_THOMPSON_VELOCITY :'//
     >                  ' VMAX IMAGINARY',
     >                  vimp,vbd,gamma,vtmp
c
c         stop 
c
c        Add a temporary correction where Vmax is set to 
c        Vbd + 1% of Vbd - this will give a slow moving particle
c        and return some velocity for the particle.
c          
         vimp = sqrt ( 1.01 * vbd**2 /  ( gamma*(1.0d0-gamma)) )
         vmax = sqrt ( gamma*(1.0d0-gamma)*vimp**2-vbd**2 )
c
      endif 
c
      vmax = sqrt(vtmp)
c
c     The next section of code is common two the various 
c     options available for the Thompson disribution. 
c
      if (thom_opt.eq.0) then  
         call evaluate_func(thom_ye,emin,emax,0.0d0,e_result,
     >                      seed,nrand,ierr)

      elseif (thom_opt.eq.1) then 
         call evaluate_func(thom_yv,vmin,vmax,0.0d0,v_result,
     >                      seed,nrand,ierr)
      endif
c
c     Energy desired is last_e - calculate the velocity
c
      if (thom_opt.eq.0) then 
         find_thompson_velocity = 1.38E4 * SQRT (e_result/CRMI) 
      elseif (thom_opt.eq.1) then   
         find_thompson_velocity = v_result
      endif
c
      return
      end
c
c
c
      subroutine evaluate_func(func,range_min,range_max,norm_val_in,
     >                         result_val,seed,nrand,ierr)
      implicit none
      integer nrand,ierr
      real*8 func,range_min,range_max,result_val,seed,norm_val_in
      external func
c
c     EVALUATE_FUNC: 
c
c     This routine takes a range minimum and maxium as well as a function.
c     It then selects a random number in the range [0,1] and finds the 
c     result_value that corresponds to the value of the range for which 
c     the integration of the function over the interval is equal to that 
c     fraction of the integration for the full range. 
c 
c     Local variables
c
      real*8 norm_val      
      real ran
      real*8 ran_frac,res_frac,last_r,rtest,rmin,rmax
      real*8 cumulative_result,current_result
      integer iter_cnt
c
      real*8 eps
      parameter(eps=0.00001d0)
c
c     Allow a value for norm_val to be passed in to save computation 
c     depending on usage. The input normalization value may be a 
c     constant expression. 
c
      norm_val = norm_val_in
c
      if (norm_val.le.0.0) then 
         call gen_qsimp(func,range_min,range_max,norm_val,0)
      endif 
c
c     Check validity of integration value
c
      if (norm_val.le.0.0) then 
         write(0,'(a,4(1x,g12.5))')
     >                  'ERROR : EVALUATE_FUNC : NORM <= 0',
     >                  range_min,range_max,norm_val
         write(6,'(a,4(1x,g12.5))')
     >                  'ERROR : EVALUATE_FUNC : NORM <= 0',
     >                  range_min,range_max,norm_val
         stop 
c
      endif 
c
c     Get random number for fraction of integration
c
      nrand = nrand + 1 
      CALL SURAND2 (SEED, 1, RAN)
      ran_frac = ran 
c
c     Perform binary search on integral until fraction is obtained
c     within error limit - start with rtest = 1/2 (range_max - range_min)
c
      rmin = range_min
      rmax = range_max 
      rtest = 0.5d0 * (rmax-rmin)
c
      res_frac = -1.0
c
c     Perform search loop - include maximum iteration test - set at 10,000 for now
c
      iter_cnt = 0
      cumulative_result  = 0.0d0
c
      do while ((abs(res_frac-ran_frac).gt.eps).and.(iter_cnt.lt.10000))
c
         iter_cnt = iter_cnt + 1
c
         call gen_qsimp(func,rmin,rtest,current_result,0)
c
         last_r   = rtest 
c
         res_frac = (cumulative_result+current_result) / norm_val 
c
c
c        Integral is less than the random fraction - need to increase the Energy
c        - add current integral to stored value
c        - set rmin to rtest
c        - set rtest to (rmin+rmax)/2
c        - integrate from rmin to rmax
c
         if (res_frac.lt.ran_frac) then 
c
            cumulative_result = cumulative_result + current_result 
            rmin = rtest
            rtest = 0.5d0*(rmin+rmax)
c
c        Integral is greater than the random fraction
c        - set rmax to rtest
c        - set rtest to (rmin + rmax)/2
c        - integrate from rmin to rmax
c
         else
c
            rmax = rtest
            rtest = 0.5d0*(rmin+rmax)
c
         endif
c
      end do
c
c     The value of the range for which the fractional condition is 
c     satisfied is stored in last_r
c
      result_val = last_r
c
c
c     return
      end
c
C
C
      SUBROUTINE GEN_QSIMP(FUNC,A,B,S,OPT)
      implicit none
      INTEGER JMAX,OPT,PLATE,ITYP
      DOUBLE PRECISION A,B,FUNC,S,EPS
      EXTERNAL FUNC
      PARAMETER (EPS=1.0D-6,JMAX=30)
C
C     GEN_QSIMP: THIS ROUTINE IS TAKEN FROM THE TEXT NUMERICAL RECIPES
C     IN FORTRAN. IT IMPLEMENTS SIMPSON'S METHOD OF CALCULATING
C     NUMERICAL QUADRATURE. IT CALLS THE ROUTINE TRAPZD TO REFINE
C     THE VALUE OF THE NUMERICAL INTEGRAL.
C
C     RETURNS AS S THE INTEGRAL OF THE FUNCTION FUNC FROM A TO B.
C     THE PARAMETER EPS CAN BE SET TO THE DESIRED FRACTIONAL ACCURACY
C     AND JMAX SO THAT 2 TO THE POWER JMAX-1 IS THE MAXIMUM ALLOWED NUMB
C     OF STEPS. INTEGRATION IS PERFORMED BY SIMPSONS RULE.
C
C     This is a generic version based on the code in soledge.d6a. The
c     main difference is that any specific function related parameters
c     are passed to the function using a common block so that the 
c     code itself can work for any function of the specified type. 
C
      INTEGER J
      DOUBLE PRECISION OS,OST,ST
c
c     Check for equal integration bounds 
c
      if (a.eq.b) then 
         s = 0.0
         return
      endif
c
      OST = -1.0D30
      OS  = -1.0D30
      DO 10 J = 1,JMAX
        CALL GEN_TRAPZD (FUNC,A,B,ST,J,OPT)

c
c       As a numerical check when dealing with small numbers
c       also check for s = os since if the numbers are small 
c       eps * abs(os) could underflow.  
c
        S = (4.0*ST-OST)/3.0
c        write(0,*) 'QSIMP:',s,os,s-os,eps*dabs(os)

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
      WRITE(6,'(a,5(1x,g12.5))') 
     >       'GEN_QSIMP: ERROR IN QUADRATURE',A,B,ST
      WRITE(6,'(a,5(1x,g12.5))') 
     >       'GEN_QSIMP:',OST,OS,ST ,S
      write(6,'(a,5(1x,g12.5))') 
     >        '    ALSO :',A,B,J,OPT
c
      WRITE(0,'(a,5(1x,g12.5))') 
     >       'GEN_QSIMP: ERROR IN QUADRATURE',A,B,ST
      WRITE(0,*) 
     >       'GEN_QSIMP: ERROR IN QUADRATURE',A,B,ST
c
      WRITE(0,'(a,5(1x,g12.5))') 
     >       'GEN_QSIMP:',OST,st,os ,S
      WRITE(0,*)
     >       'GEN_QSIMP:',OST,st,os ,S,s-os
c
      write(0,'(a,5(1x,g12.5))') 
     >        '    ALSO :',A,B,J,OPT
      STOP
C      RETURN
      END
C
C
C
      SUBROUTINE GEN_TRAPZD (FUNC,A,B,S,N,OPT)
      implicit none
      INTEGER N,OPT
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
      INTEGER IT,J
      DOUBLE PRECISION DEL,SUM,TNM,X,test
c
      IF (N.EQ.1) THEN
c
         S = 0.5*(B-A)*(FUNC(A)+FUNC(B))
C
      ELSE
         IT = 2**(N-2)
         TNM = IT
         DEL = (B-A)/TNM
         X = A + 0.5*DEL
         SUM = 0.0
         DO 10 J = 1,IT
            SUM = SUM + FUNC(X)
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
c
      double precision function thom_ye(e)
      implicit none
      double precision e
      common /thom_ye_params/  eimp,gamma,ebd
      real*8 eimp,gamma,ebd 
c
c     THOM_YE: 
c
c     This function returns the value of the Thomson 
c     energy yield function for a specific energy.
c
c
      thom_ye = e / (e+ebd)**3 * 
     >       (1.0d0 - sqrt((e+ebd)/(gamma*(1.0d0-gamma)*eimp)))
c
c      write(0,'(a,5(1x,g12.4))') 'THOM YE:',e,thom_ye,eimp,
c     >                  ebd,gamma
c
      return
      end 
c
c
c
      double precision function thom_yv(v)
      implicit none
      double precision v
      common /thom_yv_params/  vimp,vgamma,vbd
      real*8 vimp,vgamma,vbd 
c
c     THOM_YV: 
c
c     This function returns the value of the Thomson 
c     energy yield function for a specific velocity.
c
c
      thom_yv = v**3 / (v**2+vbd**2)**3 * 
     >       (1.0d0 - sqrt((v**2+vbd**2)/
     >            (vgamma*(1.0d0-vgamma)*vimp**2)))
c
c      write(0,'(a,5(1x,g12.4))') 'THOM YV:',v,thom_yv,vimp,
c     >                  vbd,vgamma
c
      return
      end 
c
c
c
      subroutine init_line_profile_data
      implicit none
      include 'params'
      include 'line_profile'
      
      !
      ! Initialize the line profile data array
      ! This is a separate routine in case future developments require more detailed initialization 
      ! or if the use in the HC module will require special initializations. 
      !

      line_profile = 0.0
     
      return 
      end
c
c
c
      subroutine update_line_profile(ik,ir,r,z,vr,vz,sputy,
     >                               cion,rizb)
      use velocity_dist
      implicit none
      integer ik,ir,cion
      real vr,vz,sputy,r,z,rizb
c
      include 'params'
      include 'cgeom'
      include 'line_profile' 
c
c     Data Required: - current particle position
c                    - current particle velocity 
c                    - current particle weight
c                    - observation position
c                    - number of lines for which a profile is required
c                    - ADAS data for each component of each line
c                    - local plasma conditions
c                    - viewing cone of the instrument
c                   
c                    
c     Calculate: - wavelength shift relative to the detector position
c                  based on the component of the particle velocity in the
c                  direction of the detector
c                - expected emission or emission weight in one time step
c                
c                - add the calculated amount of emission to the wavelength
c                  bin that corresponds to the shifted wavelength. 
c
c                - also calculate if the emission location is within the 
c                  viewing cone of the instrument
c 
c                - when the calculation is complete apply the instrument 
c                  measurement function to the final profile
c
c
c     Notes: - for now center the wavelength shift around zero - where zero is
c              is the wavelength of the unshifted profile - do this because it can
c              sometimes be difficult to properly specify the exact center 
c              wavelength - either this or use a nominal value for wavelength. 
c            - create some number of bins around the center wavelength for each
c              line and then save the data to the raw data file.            
c            - save value of option to the raw file and only load data if the
c              option was set.  
c
c
c     Local variables 
c
      real vobs,emission,wave,dlambda
      real ne,te,ti 
      logical inlos
      integer in
c
c     Write input:
c
c      write(6,'(a,3i5,2(1x,f8.5),5(1x,g12.5))') 
c     >             'UPDATE:',ik,ir,cion,r,z,vr,vz,rizb
c
c     Exit if option is not active      
c
      if (line_profile_opt.eq.0) return
c
c
c     1) Determine if particle location is within specified viewing cone from 
c        Robs, Zobs
c
      call check_particle_in_los(r,z,lp_robs,lp_zobs,
     >          lp_theta,lp_dtheta,inlos)
c
      if (.not.inlos) return
c
c     Particle is in LOS - get local ne and te
c
      ne = knbs(ik,ir)
      te = ktebs(ik,ir)
      ti = ktibs(ik,ir) 
c
c     Calculate velocity component in direction of observation position
c
      vobs = (vr * (lp_robs-r) + vz * (lp_zobs-z) )/
     >         sqrt((lp_robs-r)**2+(lp_zobs-z)**2)
c
c     If velocity debugging is on - record the velocity distribution contributing
c     to the line profile
c
c slmod begin
c     Krieger IPP/07 - SUN f90 refuses integers as logical expression
      if (debug_neutv_mod.ne.0) then
         call record_vlp_dist(vobs,sputy)
      endif
c
c      if (debug_neutv_mod) then
c         call record_vlp_dist(vobs,sputy)
c      endif
c slmod end
c
c     Call ADAS to obtain the sigmav value for the line using the data entered
c     in the input file. 
c
      call calc_adas_emission(ne,te,ti,sputy,emission,wave,
     >                        cion,rizb) 
c
c     Calculate shift in wavelength from data returned from ADAS - make the 
c     wavelength shift negative - i.e. shorter wavelengths for motion 
c     towards the detector.
c
      dlambda = -vobs/cspeed * wave 
c
c     Calulate the bin to score the emission   
c           
      if ( abs (dlambda) .le. (lp_bin_width/2.0) ) then 
         in = 0
      else  
         in =   min(  
     >           int(
     >             ((abs(dlambda)-lp_bin_width/2.0)/lp_bin_width)
     >              + 1.0),
     >               max_lp_bins)
     >            * int (sign(1.0,dlambda)) 
      endif
c
c      write (6,'(a,1x,i6,4(1x,g12.5))') 'LP:',in,wave,
c     >                  dlambda,line_profile(in),emission
c
c
c     Update appropriate bin in the line profile.  
c
      line_profile(in) = line_profile(in) + emission
c
      return
      end    
c
c
c
      subroutine calc_adas_emission(ne,te,ti,sputy,emission,wave,
     >                               cion,rizb)
      implicit none
      real ne,te,ti,sputy,emission,wave
      real rizb
      integer cion
c
      include 'params'
      include 'line_profile'
c
c     CALC_EMISSION:
c
c     This routine calls the ADAS SPEC routine to load 
c     PEC data for the spectral like specified in the 
c     input for this option. It then calculates the 
c     expected emission for the specific particle - this 
c     emission will then be scored in the appropriate 
c     wavelength shifted bin for the line profile. This 
c     routine also returns the wavelength of the line.
c
c     If the ne, te and ti are the same as the values on the last
c     call to this routine - stored values for the 
c     coefficients are used instead of calling the ADAS 
c     routines. This will save processing time when a 
c     particle remains in the same cell for some number of 
c     time steps
c
c     Local variables 
c
      real last_ne,last_te,last_pecae,last_wave
      data last_ne /0.0/ 
      data last_te /0.0/ 
      data last_pecae /0.0/ 
      data last_wave /0.0/ 
c
c      real last_pecar,last_pecax,last_ti
c      data last_pecar /0.0/ 
c      data last_pecax /0.0/ 
c      data last_ti /0.0/ 
c
c     ADAS variables
c
      real*8 te_adas,ne_adas,ti_adas
      real*8 wlngth,pecae,pecar,pecax
      logical*4 ltrng,ldrng  
      integer iz,cz,ircode,n_data_points
c
c      INTEGER   IR,IK,IADAS,NPAIRS,IKK
c      REAL*8    TADAS(20),DADAS(20)
c      REAL*8    WLNGTH,PECAE(20),PECAR(20),PECAX(20)
c      LOGICAL*4 LTRNG(20),LDRNG(20)
c
      CHARACTER ADASGR*8,ADASTY*80,PECTITLE*120
      CHARACTER XFESYM*2
      external xfesym
c
c     Set to appropriate species 
c     Set to neutral lines only for now
c
c
c     Initialization
c
c     Load for neutral impurity line only for now for a 
c     single data point. 
c
c      write(0,'(a,i5,6(1x,g12.5))') 'EMISS:',cion,te,ti,ne,
c     >      wave, rizb 
c
      cz = cion
      iz = 0
      n_data_points = 1
c
      WAVE = 0.0
      IRCODE = 0
      te_adas = te
      ti_adas = ti
      ne_adas = ne
      emission = 0.0
c
      CALL XXUID(LP_ADASID)
      IF (LP_ADASYR.GE.0) THEN
        ADASGR = 'pec??#'//XFESYM(CZ)
        WRITE(ADASGR(4:5),'(I2.2)') lp_ADASYR
      ELSE
        ADASGR = '*'
      ENDIF
      ADASTY = '*'
c
c      write(0,'(a,2i5,4(1x,a))') 
c     >        'ADAS:',lp_isele,len(lp_adasex),
c     >        adasgr,adasty,lp_adasex
c
      CALL XXSPEC(ADASGR,ADASTY,lp_ADASEX)
c
c     Load emission due to excitation. 
c
c      write(0,*) 'Loading PEC:'
c
      IF ((last_ne.ne.ne.or.last_te.ne.te).and.
     >     lp_ISELE.GT.0) THEN
         CALL SPEC(lp_ISELE,IZ,CZ,n_data_points,
     >             Te_adas,ne_ADAS, 
     >             WLNGTH,PECAE,LTRNG,LDRNG,PECTITLE,
     >             IRCODE)
         IF (IRCODE.NE.0) then 
            write(6,*) 
     >          'WARNING: NON-ZERO ADAS RETURN CODE=',ircode
            write(0,*) 
     >          'WARNING: NON-ZERO ADAS RETURN CODE=',ircode
            return
         endif  
c
c        Record values for reference on the next call
c
         last_ne = ne
         last_te = te 
         last_pecae = pecae
         last_wave = wlngth 
c
      else  
c
         pecae = last_pecae
         wlngth = last_wave
c
      endif 
c
c     Calculate emission and assign wavelength
c
      emission = 1.0e-6 * (pecae * rizb * ne_adas * sputy)  
c
c      write(0,*) 'Emission:',emission,pectitle
C
      WAVE = WLNGTH
      lp_wave = wlngth 
c
      return 
    
c
c     REFERENCE CODE:  
c
c     This code is included to show how the emission is 
c     calculated in OUT. 
c
c     RECOMBINATION AND CX COMPONENTS - IF EVER ADDED - MUST 
c     BE BASED ON THE ION VELOCITY AND WILL REQUIRE A NEW
c     ROUTINE.
c
c     The recombination and charge exchange components of 
c     the usual emission profile are not useful when calculating
c     the emission on a time step by time step basis since these 
c     events do not happen on a regular basis. 
c
c     
c
c      IF (ISELR.GT.0) THEN
c         CALL SPEC(ISELR,IZ,CZ,n_data_points,Te_adas,ne_adas,
c     >             WLNGTH,PECAR,LTRNG,LDRNG,PECTITLE,IRCODE)
c         IF (IRCODE.NE.0) RETURN
c      endif
c
c      IF (ISELX.GT.0) THEN
CLDH - USE ION TEMPERATURE FOR CX RATE
c         CALL SPEC(ISELX,IZ,CZ,n_data_points,Te_ADAS,ne_ADAS,
c     >             WLNGTH,PECAX,LTRNG,LDRNG,PECTITLE,IRCODE)
c         IF (IRCODE.NE.0) RETURN
c      endif
c 
c            IF (CZ.GT.1.0) THEN
c              CVALS(IKK,IR) = 1.D-6*
c     >        (PECAE(IADAS) * RIZB * KNBS(IKK,IR) * SDLIMS(IKK,IR,IZ)
c     >        +PECAR(IADAS) * RIZB * KNBS(IKK,IR) * SDLIMS(IKK,IR,IZ+1)
c     >        +PECAX(IADAS) * PINATOM(IKK,IR) * SDLIMS(IKK,IR,IZ+1))
c            endif  
c
      end    
c
c
c
      subroutine check_particle_in_los(r,z,robs,zobs,theta,dtheta,inlos)
      implicit none
      real r,z,robs,zobs,theta,dtheta
      logical inlos  
c
      include 'params' 
c
c
c     CHECK_PARTICLE_IN_LOS:
c 
c     This routine checks to see if the position R,Z is within a viewing region defined
c     by the observation position robs,zobs at an angle theta and width of dtheta.
c
      real theta_top,theta_bot,theta_part
      real atan2c
      external atan2c
c 
c      write(0,*) 'Checking LOS:'
c
c     Find theta_part and make sure it is in the range 0.0 to 2PI
c
      theta_part = atan2c(z-zobs,r-robs)
      if (theta_part.lt.0.0) theta_part = theta_part + 2.0 * PI
c
c     Check to see if theta_part lies between theta-dtheta/2 and theta+dtheta/2
c     Check to see if the observation cone cross the 0.0 degree line
c
      theta_top = theta + dtheta/2.0 
      theta_bot = theta - dtheta/2.0 
c
c
c     Check for wrap across 0.0 degrees
c
      if (theta_top.gt.2.0*PI) then 
c
         if (theta_part.le.(theta_top-2.0*PI).or.
     >       theta_part.ge.theta_bot) then 
            inlos=.true.
         else
            inlos=.false.
         endif
c
c     No wrap - easy check
c
      else
c
         if (abs(theta-theta_part).le.(dtheta/2.0)) then 
            inlos = .true.
         else
            inlos = .false.
         endif
c
      endif
c
c
      return
      end






