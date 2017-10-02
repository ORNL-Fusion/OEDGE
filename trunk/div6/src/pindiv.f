c     -*Fortran*-
c
      SUBROUTINE WRTPIN (title,equil,IOUNIT)
      use debug_options
      IMPLICIT none
      character*(*) title,equil                                 
      INTEGER IOUNIT                                                    
c
      include 'params'                                                  
      include 'comtor'                                                  
      include 'cgeom'                                                   
      include 'dynam4'                                                  
      include 'cadas'                                                  
c
      include 'pindata'
      include 'cedge2d' 
c
c  TEMPORARY LOCAL VARIABLES - LDH DIV4.07?
c
c      integer ihybrid, jhpuf1(2), jhpuf2(2)
c      real    ppcpuf, hpcpuf, tpufh
c
      logical lpvhpf
      integer recsigopt      
      real    target_flux
      integer calc_random_seed
      external calc_random_seed
c
C                                                                       
C     RECONSTRUCT THE FULL GEOMETRY (I.E. INCLUDING VIRTUAL POINTS)     
C     AND CONVERT PLASMA PROFILES TO CGS.                               
C     THIS WRITES OUT THE FILE PIN.DIV - NEEDED FOR THE PIN CODE        
C     IMPLEMENTATION.                                                   
C                                                                       
C     THE FOLLOWING ARE USED ONLY TO CREATE THE PIN FILE.               
C                                                                       
      INTEGER    MAXK ,IK,IR,K,NP,ID,IN,IFIRST,tmpcion,tmpsput
      PARAMETER (MAXK=4000)                                             
      REAL       DENEL(MAXK),PREEL(MAXK),DEN(MAXK),PRE(MAXK)              
      REAL       VTE(MAXK),VR0(MAXK),FLUXPX(MAXK),FLUXPY(MAXK)          
      REAL       CS, KKFLUX, FLXOUT                                  
      real       totrec
      DATA       IFIRST/0/  
C                                                                       
      call pr_trace('PINDIV','WRTPIN START')

      REWIND(IOUNIT)                                                    
c
c     SET UP PIN RANDOM NUMBER SEED - CURRENTLY LIMITED TO 6 CHARACTERS
c     IN NIMBUS.
c
c     IF PINISEED WAS INPUT AS 0 - A NEW SEED IS GENERATED EVERY TIME
c                             <0 - THE NIMBUS DEFAULT SEED OF 1 is USED
c                             >0 - THE SPECIFIED SEED IS USED - IT WILL
c                                  be TRUNCATED IN NIMBUS TO "I6"
c
      if (piniseed.eq.0) then  
         pinseed = calc_random_seed(6)
      elseif (piniseed.lt.0) then
         pinseed = 1  
      else
         pinseed = piniseed
      endif 

c
c     SET-UP PUFFING OPTION VALUES IF PIN PUFFING HAS BEEN SELECTED. 
c
c
c     Puffing-Off - set all puff related quantities to OFF
c
      if (pinpuff.eq.0) then  
         hpcpuf = 0.0
         ppcpuf = 0.0 
         tpufh =  0.1 
c
c         hextrl = 0.0
c         phxtra = 0.0
c
         lpvhpf = .false.
         jhpuf1(1) = -1000
         jhpuf1(2) = -1000
         jhpuf2(1) = -1
         jhpuf2(2) = -1
      elseif (pinpuff.eq.1) then
c
c        Zero Option 2 inputs
c
         ppcpuf = 0.0
c
c
c        Old implementation
c 
c         hextrl = 0.0
c         phxtra = 0.0
c
c        Set location switch
c
         if (swpvhpf.eq.0) then 
            lpvhpf = .false.
         elseif (swpvhpf.eq.1) then 
            lpvhpf = .true.
         endif
c
      elseif (pinpuff.eq.2) then
c
c        Zero Option 1 inputs
c
         hpcpuf = 0.0
c
c        Set location switch
c
         if (swpvhpf.eq.0) then 
            lpvhpf = .false.
         elseif (swpvhpf.eq.1) then 
            lpvhpf = .true.
         endif
      endif
c
c     Set the recombination option to match the required NIMBUS 
c     option. (Check ENTRY RECSIG in NIMBXS subroutine in NIMBUS.P4A
c     module for the NIMBUS options - at the present time they match
c     the mapping below.)
c
c     DIVIMP PATCH - ensures same recombination option in NIMBUS and
c                    DIVIMP
c                  - REQUIRED so that we can match exact conditions
c                    of Edge2D runs that may or may not use ADAS 
c                    based recombination
c
      if (crecopt.ge.0.and.crecopt.le.4) then 
c       
c        Other options
c
c        0 = OFF
c        1 = Gordeev
c        2 = Janev
c        3 = NRL
c        4 = ADAS 
c
         recsigopt = crecopt
c
      endif
c
c     Record values from LAST PIN iteration
c      
c
      hescpd_last = hescpd
      hescal_last = hescal   
      phfgal_last = phfgal 
      phfuga_last = phfuga
c 
C
C     Calculate the number of recombinations
C     Store the recombination array for later use.
c     This is needed for calculating the total source strength for 
c     the purposes of puffing.
C
      call calc_divrec(totrec)
C                                                                       
C                                                                       
C  USE KORY TO MAP THE DIVIMP MATRICES INTO NIMBUS VECTORS.  IF THE     
C  VIRTUAL POINTS HAVE BEEN REMOVED (VIRTGRID = .FALSE.) THEN NOTHING   
C  IS WRITTEN FOR THE VIRTUAL POINTS, WHICH AREN'T USED BY NIMBUS IN    
C  ANY CASE.                                                            
C                                                                       
      CALL RZERO(DENEL,MAXK)                                              
      CALL RZERO(PREEL,MAXK)                                              
      CALL RZERO(DEN,MAXK)                                             
      CALL RZERO(PRE,MAXK)                                             
      CALL RZERO(VTE,MAXK)                                              
      CALL RZERO(VR0,MAXK)                                              
      NP = 0                                                            

      call pr_trace('PINDIV','WRTPIN BEFORE PLASMA')

      DO 600 IR = 1,NRS                                                 
        NP = NP + NKS(IR)                                               
        DO 620 IK = 1,NKS(IR)                                           
c
c         changed by Krieger, 11/94 - check for 0'th index 
c
          if (kory(ir,ik).ne.0) then	
            DENEL(KORY(IR,IK))  = KNBS(IK,IR) * 1.0E-6                      
            PREEL(KORY(IR,IK))  = KNBS(IK,IR)*KTEBS(IK,IR)*ECH * 1.0E1
            DEN(KORY(IR,IK)) = DENEL(KORY(IR,IK))                          
            PRE(KORY(IR,IK)) = KNBS(IK,IR)*KTIBS(IK,IR)*ECH *  1.0E1
c            VTE(KORY(IR,IK))  = KVHS(IK,IR) / QTIM * 1.0E2     
            VTE(KORY(IR,IK))  = KVHS(IK,IR) * 1.0E2     
            VR0(KORY(IR,IK))  = 0.0
          endif				
c
 620    CONTINUE                                                        
 600  CONTINUE                                                          
      IF (.NOT.VIRTGRID) NP = NP + 2*((IRWALL-IRSEP+1)+(NRS-IRTRAP+1))  
C                                                                       
C  LINKPG/NIMBUS ARE NOW DESIGNED TO BE PASSED THE ION FLOW ACROSS      
C  THE PLASMA BOUNDARIES, NOT THE THE FLUX PER UNIT AREA.  THE          
C  PERPENDICULAR ION FLOW TO THE WALLS IS ASSUMED TO BE ZERO.  THE      
C  PARALLEL FLOW (PER UNIT LENGTH TOROIDALLY) IS CALCULATED USING THE   
C  TARGET SEGMENT LENGTHS FROM THE PLASMA POLYGON VERTEX POSITIONS.     
c  The new intfac (7/95) requires the total flux, i.e. not per unit
c  length toroidally.
C
C  10/10/95 - Modify to conform with TARGFLUX (in TAU) so that we have
C             only one definition of the flux to the target - ldh
C                                                                       
      CALL RZERO(FLUXPX,MAXK)                                           
      CALL RZERO(FLUXPY,MAXK)                                           
      flxout = 0.0

      call pr_trace('PINDIV','WRTPIN BEFORE FLUX')

C     
      DO 640 ID = 1,NDS
         IK = IKDS(ID)
         IR = IRDS(ID)
         in = 2
         if (ik.gt.nks(ir)/2) in = 1
c
c         CS = 9.79E3 * SQRT (0.5*(KTEDS(ID)+KTIDS(ID))*              
c     >                 (1.0+RIZB)/CRMB)
c
c - sol option 22 or 23 allows mach numbers at the target > 1
c
c         if (cioptf.eq.22.or.cioptf.eq.23) cs = cs * cmachno(ir,in)
c
         cs = abs(kvds(id))
c
c        Set flux - may be different due to E2D target condition option
c
         if (e2dtargopt.eq.3) then 
            target_flux = e2dtarg(ir,6,in) 
         else 
            target_flux = knds(id) * cs
         endif
c
         if (northopt.eq.0) then                                     
            KKFLUX = target_flux / KBFST(IR,IN)                     
         elseif (northopt.eq.1.or.northopt.eq.2) then
            KKFLUX = target_flux / KBFST(IR,IN) * COSTET(ID)       
         endif                                                       
c
c         KKFLUX = KNDS(ID) * CS / KBFST(IR,IN) * COSTET(ID)       
c                                                                       
c        changed by Krieger, 11/94 - added kory test
c
         if (kory(ir,ik).ne.0) then
           FLUXPY(KORY(IR,IK)) = KKFLUX * DDS2(ID) * 2.0*PI*RP(ID)
           flxout = flxout + fluxpy(kory(ir,ik))
         endif
c
 640  CONTINUE                                                          
c
c
c  flxout has a different normalisation in intfac!
c
      flxout = flxout / (2.*pi)
c
c  the flux of atoms+ions to the wall should be zeroed before
c  the first call to PIN.  After that, the value calculated
c  in the last call should be returned as input for the next.
c
      if (ifirst.eq.0) then
        call rzero(flxhw2,maxseg)
c
c       Set the initial pump/puff portion to zero. 
c
        phfgal = 0.0         
        phfuga = 0.0
c
        ifirst=1
c
c       On the first iteration HPCPUF is always effectively zero 
c       After that it uses the input value. 
c    
        acthpcpuf = 0.0
c
c       Initialize to zero on first iteration
c
        hescpd_last = 0.0
        hescal_last = 0.0
        phfgal_last = 0.0 
        phfuga_last = 0.0
c
      else
c
        acthpcpuf = hpcpuf
c
      endif
c
c
c     NIMBUS can not deal with recycling impurities being 
c     passed into it in CION. SO - the quantity passed to 
c     NIMBUS for CION can only be those which NIMBUS can 
c     deal with - at the moment the only one of these
c     (that I know of) is carbon - so tmpcion will be set 
c     to 6. 
c
c     Carbon
c
      if (cion.eq.6) then 
         tmpcion = cion
      else
         tmpcion = 6
      endif
c
c     Added modified sputtering data - but PIN/NIMBUS probably ONLY has the 
c     two - so make sure that only 1 or 2 gets passed. 
c
      if (csputopt.eq.2.or.csputopt.eq.3.or.csputopt.eq.5) then 
         tmpsput = 2
      else
         tmpsput = 1
      endif
c
c  LDH DIV4.07?
c  TEMPORARILY, hardwire these local variables and let David sort out 
c               getting them properly passed from the input file :-)
c
c     ihybrid = 0  no action, just use vessel wall in equil. file as before
c     ihybrid = 1  modify Mark I vessel wall to 'remote wall' for neutrals
c     ihybrid = 2  modify Mark IIA vessel wall to 'remote wall' for neutrals
c
c      ihybrid = 2
c
c     ppcpuf:  fraction of recycling ions which are not recycled but
c              rather reinserted as a puff - default 0.0
c      ppcpuf = 0.0
c
c     hpcpuf:  fraction of pumped (albedoed) neutrals which are
c              reinserted as a puff - default 0.0
c      hpcpuf = 1.0
c
c     lpvhpf:  logical to set H puff in private region (.TRUE.) or
c              onto SOL (.FALSE.) - default .FALSE.
c
c      lpvhpf = .TRUE. 
c
c     tpufh:   energy in eV of puffed gas particles - default 0.1 eV
c
c      tpufh = 0.1
c
c     jhpuf:   range of knots from which to puff.  Definition
c              depends on lpvhpf - default (puffing on SOL divertor
c              legs) jhpuf1(1)=-1000; jhpuf2(1)=-1
c                    jhpuf1(2)=-1000; jhpuf2(2)=-1
c
c      jhpuf1(1) = 1
c      jhpuf2(1) = 1
c      jhpuf1(2) = -1000
c      jhpuf2(2) = -1
c

      call pr_trace('PINDIV','WRTPIN BEFORE WRITE')

C                                                                       
C     WRITE OUT THE BACKGROUND PLASMA PROFILES FOR THE PIN CODE.        
C
      write(iounit,'(a)') title
      write(iounit,'(i6,a)') ishot,equil
      write(iounit,'(i6,a)') iyearh,useridh
      write(iounit,'(i6,a)') iyearz,useridz
      WRITE(IOUNIT,'(1P,2e12.4,2i3)') CRMB, flxout,TMPCION,tmpsput-1
      write(iounit,'(2(i2))') ihcorr,ihybrid
      write(iounit,'(3e12.4)') phfgal,phfuga,(totrec*2.0*PI*R0)
      write(iounit,'(l6)') lpvhpf
c
c     Modified for new NIMBUS version - Feb 25/98
c
      write(iounit,'(3e12.4)') ppcpuf,hpcpuf,tpufh
c
      write(iounit,'(2i6)') jhpuf1(1),jhpuf2(1)
      write(iounit,'(2i6)') jhpuf1(2),jhpuf2(2)
      write(iounit,'(2(i6))') recsigopt,pinseed
      WRITE(IOUNIT,'(1P,6E12.4)') (DENEL(K),K = 1,NP)                     
      WRITE(IOUNIT,'(1P,6E12.4)') (PREEL(K),K = 1,NP)                     
      WRITE(IOUNIT,'(1P,6E12.4)') (DEN(K),K = 1,NP)                    
      WRITE(IOUNIT,'(1P,6E12.4)') (PRE(K),K = 1,NP)                    
      WRITE(IOUNIT,'(1P,6E12.4)') (VTE(K),K = 1,NP)                     
      WRITE(IOUNIT,'(1P,6E12.4)') (VR0(K),K = 1,NP)                     
      WRITE(IOUNIT,'(1P,6E12.4)') (FLUXPX(K),K = 1,NP)                  
      WRITE(IOUNIT,'(1P,6E12.4)') (FLUXPY(K),K = 1,NP)                  
      WRITE(IOUNIT,'(1p,6E12.4)') (FLXHW2(K),K = 1,MAXSEG)
      CLOSE (IOUNIT)                                                    
C                                                                       
      RETURN                                                            
      END                                                               
C
C
      SUBROUTINE RINOUTU (OPT,RARRAY,N,IONUM)
C
C  *********************************************************************
C  *  RINOUTU: READS IN / WRITES OUT AN UNFORMATTED ARRAY OF REALS.    *
C  *  THE ARRAYS ARE READ/WRITTEN ON CHANNEL IONUM, TO A DATASET WITH  *
C  *  ATTRIBUTES BLKSIZE=6160, RECFM=VBS, LREC=6160, TRKS=(20,20)      *
C  *  OPT(1:1) SHOULD BE 'R' OR 'W', AND OPT(3:8) IS THE ARRAY NAME    *
C  *  (USED IN WRITE STATEMENT AT END OF ROUTINE).                     *
C  *                                                                   *
C  *          CHRIS FARRELL    MARCH 1989                              *
C  *          DAVID ELDER      JUNE  1992                              *
C  *********************************************************************
C
      INTEGER I,J,N,IBLOCK,IONUM
      CHARACTER OPT*8
      REAL RARRAY(N)
      DATA IBLOCK /1500/
C
      IF     (OPT(1:1).EQ.'R') THEN
        DO 100 I = 1, N, IBLOCK
          READ  (IONUM) (RARRAY(J),J=I,MIN(N,I+IBLOCK-1))
  100   CONTINUE
      ELSEIF (OPT(1:1).EQ.'W') THEN
        DO 200 I = 1, N, IBLOCK
          WRITE (IONUM) (RARRAY(J),J=I,MIN(N,I+IBLOCK-1))
  200   CONTINUE
      ENDIF
C      IF (4*N.GT.10000) WRITE (6,9001) OPT(3:8),REAL(4*N)
C 9001 FORMAT(1X,'RINOUT: SIZE OF ',A6,' =',-6P,F6.2,' MB')
C
      RETURN
      END
C
C
C
      SUBROUTINE PINPRN
C     INCLUDE "PARAMS"
      include 'params'
C     INCLUDE "PINDATA"
      include 'pindata'
C     INCLUDE "CGEOM"
      include 'cgeom'
      INTEGER I,J
 
      CALL PRRMATDIV(PINION,MAXNKS,nks(irsep),NRS,24,'IONIZATION')
      CALL PRRMATDIV(PINATOM,MAXNKS,nks(irsep),NRS,24,'NEUTRALS')
      CALL PRRMATDIV(PINALPHA,MAXNKS,nks(irsep),NRS,24,'HALPHA')
      CALL PRRMATDIV(PINMOL,MAXNKS,nks(irsep),NRS,24,'MOLECULES')
      CALL PRRMATDIV(PINZ0,MAXNKS,nks(irsep),NRS,24,'IMPURITY')
      CALL PRRMATDIV(PINIONZ,MAXNKS,nks(irsep),NRS,24,'IMP IONIZ')
      CALL PRRMATDIV(PINENZ,MAXNKS,nks(irsep),NRS,24, 'IMP ENERGY')
      call prrmatdiv(PINENA,maxnks,nks(irsep),nrs,24,'NEUT ENERGY')
C 
      RETURN
      END
C
C
      SUBROUTINE PINCHK
      IMPLICIT NONE
C     INCLUDE "PARAMS"
      include 'params'
C     INCLUDE "PINDATA"
      include 'pindata'
C     INCLUDE "CGEOM"
      include 'cgeom'
C     INCLUDE "CNEUT"
c      include 'cneut'
C     INCLUDE "COMTOR"
      include 'comtor'
      include 'printopt'
C
C     PINCHK: THIS ROUTINE TAKES THE IONIZATION DATA FROM PIN AND
C             PERFORMS VARIOUS SUMS OVER VOLUME AND THEN COMPARES THESE
C             TO CALCULATED PLATE FLUXES BASED ON SOUND SPEEDS AND PLATE
C             PLASMA DENSITIES.  IN ADDITION, CALCULATE THE TOTAL
C             VOLUME RECOMBINATION IN THE PLASMA.
C
C             MODIFIED 10/10/95 TO CONFORM TO THE DEFINITION OF THE 
C             TARGET FLUX IN TARGFLUX (IN TAU) - LDH
C
c
c     jdemod - 2012 - in the code kareas and karea2 are set to be the same 
c                     values which are calculated from the grid polygon areas
c                   - as a result - the code is being modified to output the 
c                     second set of data based on kvols rather than karea2
c
 
      INTEGER IR,IK,IN,IKMID,ID,IKTOP
      REAL TVR(MAXNRS,3),TVR2(MAXNRS,3),TIZR(MAXNRS,3),TIZR2(MAXNRS,3),
     >     TVS(3)       ,TVS2(3)       ,TIZS(3)       ,TIZS2(3)       ,
     >     TVC(3)       ,TVC2(3)       ,TIZC(3)       ,TIZC2(3)       ,
     >     TV(3)        ,TV2(3)        ,TIZ(3)        ,TIZ2(3)
      REAL TRCR(MAXNRS,3), TRCR2(MAXNRS,3), TRCS(3), TRCS2(3),
     >     TRCC(3)       , TRCC2(3)       , TRC(3) , TRC2(3)
      REAL TVD(3)       ,TVD2(3)       ,TIZD(3)       ,TIZD2(3),
     >     TRCD(3)      ,TRCD2(3)
      REAL FLUXR(MAXNRS,2),DISTR(MAXNRS,2),DISTR2(MAXNRS,2),
     >     SRCR(MAXNRS,3),SRCR2(MAXNRS,3),
     >     DIST(2),DIST2(2),SRC(3),SRC2(3)
      REAL CELVOL,CELVOL2,CS
 
C
C     INITIALIZATION
C
      CALL RZERO(TVR,MAXNRS*3)
      CALL RZERO(TVR2,MAXNRS*3)
      CALL RZERO(TIZR,MAXNRS*3)
      CALL RZERO(TIZR2,MAXNRS*3)
      CALL RZERO(TRCR,MAXNRS*3)
      CALL RZERO(TRCR2,MAXNRS*3)
      CALL RZERO(FLUXR,MAXNRS*2)
      CALL RZERO(DISTR,MAXNRS*2)
      CALL RZERO(DISTR2,MAXNRS*2)
      CALL RZERO(SRCR,MAXNRS*3)
      CALL RZERO(SRCR2,MAXNRS*3)
      DO 5 IN = 1,3
        TVS(IN) = 0.0
        TVS2(IN) = 0.0
        TIZS(IN) = 0.0
        TIZS2(IN) = 0.0
        TRCS(IN) = 0.0
        TRCS2(IN) = 0.0
        TVC(IN) = 0.0
        TVC2(IN) = 0.0
        TIZC(IN) = 0.0
        TIZC2(IN) = 0.0
        TRCC(IN) = 0.0
        TRCC2(IN) = 0.0
        TV(IN) = 0.0
        TV2(IN) = 0.0
        TIZ(IN) = 0.0
        TIZ2(IN) = 0.0
        TRC(IN) = 0.0
        TRC2(IN) = 0.0
        TVD(IN) = 0.0
        TVD2(IN) = 0.0
        TIZD(IN) = 0.0
        TIZD2(IN) = 0.0
        TRCD(IN) = 0.0
        TRCD2(IN) = 0.0
        SRC(IN) = 0.0
        SRC2(IN) = 0.0
 5    CONTINUE
      DO 6 IN = 1,2
        DIST(IN) = 0.0
        DIST2(IN) = 0.0
 6    CONTINUE
C
C     CALCULATE TOTAL IONIZATION AND IONIZATION/RING
C
      DO 10 IR = 1,NRS
C
C     OUTER HALF OF PLASMA
C
        IN = 2
        IKMID = NKS(IR) / 2
        IKTOP = NKS(IR)
        DO 20 IK = 1,IKMID
C
C         IONIZATION USING DIVIMP AREAS
C
          CELVOL = KAREAS(IK,IR)
          TIZR(IR,IN) = TIZR(IR,IN) + PINION(IK,IR)*CELVOL
          TIZ(IN)     = TIZ(IN)     + PINION(IK,IR)*CELVOL
          TVR(IR,IN)  = TVR(IR,IN)  + CELVOL
          TV(IN)      = TV(IN)      + CELVOL
C
C         RECOMBINATION USING DIVIMP AREAS
C
          TRCR(IR,IN) = TRCR(IR,IN) + PINREC(IK,IR)*CELVOL
          TRC(IN)     = TRC(IN)     + PINREC(IK,IR)*CELVOL
C
C         IONIZATION USING TRUE AREAS (AS USED IN NIMBUS)
C
c         jdemod - change output to use kvols 
c          CELVOL2 = KAREA2(IK,IR)
c
          CELVOL2 = KVOLS(IK,IR)
          TIZR2(IR,IN) = TIZR2(IR,IN) + PINION(IK,IR)*CELVOL2
          TIZ2(IN)     = TIZ2(IN)     + PINION(IK,IR)*CELVOL2
          TVR2(IR,IN)  = TVR2(IR,IN)  + CELVOL2
          TV2(IN)      = TV2(IN)      + CELVOL2
C
C         RECOMBINATION USING TRUE AREAS (AS USED IN NIMBUS)
C
          TRCR2(IR,IN) = TRCR2(IR,IN) + PINREC(IK,IR)*CELVOL2
          TRC2(IN)     = TRC2(IN)     + PINREC(IK,IR)*CELVOL2
C
C         IONIZATION GREATER AND LESS THAN ZXP IN SOL
C
          IF (IR.GE.IRSEP) THEN
            IF (ZS(IK,IR).LT.ZXP) THEN
              TIZD(1)  = TIZD(1)  + PINION(IK,IR)*CELVOL
              TRCD(1)  = TRCD(1)  + PINREC(IK,IR)*CELVOL
              TVD(1)   = TVD(1)   + CELVOL
              TIZD2(1) = TIZD2(1) + PINION(IK,IR)*CELVOL2
              TRCD2(1) = TRCD2(1) + PINREC(IK,IR)*CELVOL2
              TVD2(1)  = TVD2(1)  + CELVOL2
            ELSE
              TIZD(2)  = TIZD(2)  + PINION(IK,IR)*CELVOL
              TRCD(2)  = TRCD(2)  + PINREC(IK,IR)*CELVOL
              TVD(2)   = TVD(2)   + CELVOL
              TIZD2(2) = TIZD2(2) + PINION(IK,IR)*CELVOL2
              TRCD2(2) = TRCD2(2) + PINREC(IK,IR)*CELVOL2
              TVD2(2)  = TVD2(2)  + CELVOL2
            ENDIF
          ELSE
C
C         IONIZATION IN CORE
C
            TIZC(IN)  = TIZC(IN)  + PINION(IK,IR) * CELVOL
            TRCC(IN)  = TRCC(IN)  + PINREC(IK,IR) * CELVOL
            TVC(IN)   = TVC(IN)   + CELVOL
            TIZC2(IN) = TIZC2(IN) + PINION(IK,IR) * CELVOL2
            TRCC2(IN) = TRCC2(IN) + PINREC(IK,IR) * CELVOL2
            TVC2(IN)  = TVC2(IN)  + CELVOL2
          ENDIF
 
 20     CONTINUE
C
C       SUM TOTAL IONIZATION OUTER PLATE
C
        IF (IR.GE.IRSEP) THEN
          TIZS(IN)  = TIZS(IN)  + TIZR(IR,IN)
          TRCS(IN)  = TRCS(IN)  + TRCR(IR,IN)
          TVS(IN)   = TVS(IN)   + TVR(IR,IN)
          TIZS2(IN) = TIZS2(IN) + TIZR2(IR,IN)
          TRCS2(IN) = TRCS2(IN) + TRCR2(IR,IN)
          TVS2(IN)  = TVS2(IN)  + TVR2(IR,IN)
        ENDIF
C
C       INNER HALF OF PLASMA
C
        IN = 1
        DO 25 IK = IKMID+1,IKTOP
C
C         REMEMBER THAT FIRST AND LAST CELL ON CLOSED RINGS ARE THE SAME
C
          IF (IR.LT.IRSEP .AND. IK.EQ.NKS(IR)) GOTO 25
C
C         IONIZATION USING DIVIMP VOLUMES
C
          CELVOL = KAREAS(IK,IR)
          TIZR(IR,IN) = TIZR(IR,IN) + PINION(IK,IR)*CELVOL
          TIZ(IN)     = TIZ(IN)     + PINION(IK,IR)*CELVOL
          TVR(IR,IN)  = TVR(IR,IN)  + CELVOL
          TV(IN)      = TV(IN)      + CELVOL
C
C         RECOMBINATION USING DIVIMP VOLUMES
C
          TRCR(IR,IN) = TRCR(IR,IN) + PINREC(IK,IR)*CELVOL
          TRC(IN)     = TRC(IN)     + PINREC(IK,IR)*CELVOL
C
C         IONIZATION USING TRUE VOLUMES (AS USED IN NIMBUS)
C
c         jdemod - 
c          CELVOL2 = KAREA2(IK,IR)
          CELVOL2 = KVOLS(IK,IR)
c
          TIZR2(IR,IN) = TIZR2(IR,IN) + PINION(IK,IR)*CELVOL2
          TIZ2(IN)     = TIZ2(IN)     + PINION(IK,IR)*CELVOL2
          TVR2(IR,IN)  = TVR2(IR,IN)  + CELVOL2
          TV2(IN)      = TV2(IN)      + CELVOL2
C
C         RECOMBINATION USING TRUE VOLUMES (AS USED IN NIMBUS)
C
          TRCR2(IR,IN) = TRCR2(IR,IN) + PINREC(IK,IR)*CELVOL2
          TRC2(IN)     = TRC2(IN)     + PINREC(IK,IR)*CELVOL2
C
C         IONIZATION GREATER AND LESS THAN ZXP FOR SOL
C
          IF (IR.GE.IRSEP) THEN
            IF (ZS(IK,IR).LT.ZXP) THEN
              TIZD(1)  = TIZD(1)  + PINION(IK,IR) * CELVOL
              TRCD(1)  = TRCD(1)  + PINREC(IK,IR) * CELVOL
              TVD(1)   = TVD(1)   + CELVOL
              TIZD2(1) = TIZD2(1) + PINION(IK,IR) * CELVOL2
              TRCD2(1) = TRCD2(1) + PINREC(IK,IR) * CELVOL2
              TVD2(1)  = TVD2(1)  + CELVOL2
            ELSE
              TIZD(2)  = TIZD(2)  + PINION(IK,IR) * CELVOL
              TRCD(2)  = TRCD(2)  + PINREC(IK,IR) * CELVOL
              TVD(2)   = TVD(2)   + CELVOL
              TIZD2(2) = TIZD2(2) + PINION(IK,IR) * CELVOL2
              TRCD2(2) = TRCD2(2) + PINREC(IK,IR) * CELVOL2
              TVD2(2)  = TVD2(2)  + CELVOL2
            ENDIF
          ELSE
C
C         IONIZATION IN CORE
C
            TIZC(IN)  = TIZC(IN)  + PINION(IK,IR) * CELVOL
            TRCC(IN)  = TRCC(IN)  + PINREC(IK,IR) * CELVOL
            TVC(IN)   = TVC(IN)   + CELVOL
            TIZC2(IN) = TIZC2(IN) + PINION(IK,IR) * CELVOL2
            TRCC2(IN) = TRCC2(IN) + PINREC(IK,IR) * CELVOL2
            TVC2(IN)  = TVC2(IN)  + CELVOL2
          ENDIF
 25     CONTINUE
C
C       TOTAL IONIZATION AT INNER PLATE
C
        IF (IR.GE.IRSEP) THEN
          TIZS(IN)  = TIZS(IN)  + TIZR(IR,IN)
          TRCS(IN)  = TRCS(IN)  + TRCR(IR,IN)
          TVS(IN)   = TVS(IN)   + TVR(IR,IN)
          TIZS2(IN) = TIZS2(IN) + TIZR2(IR,IN)
          TRCS2(IN) = TRCS2(IN) + TRCR2(IR,IN)
          TVS2(IN)  = TVS2(IN)  + TVR2(IR,IN)
        ENDIF
C
C       TOTAL OF INNER AND OUTER HALVES
C
        TVR(IR,3)   = TVR(IR,1)   + TVR(IR,2)
        TVR2(IR,3)  = TVR2(IR,1)  + TVR2(IR,2)
        TIZR(IR,3)  = TIZR(IR,1)  + TIZR(IR,2)
        TIZR2(IR,3) = TIZR2(IR,1) + TIZR2(IR,2)
        TRCR(IR,3)  = TRCR(IR,1)  + TRCR(IR,2)
        TRCR2(IR,3) = TRCR2(IR,1) + TRCR2(IR,2)
        TV(3)   = TV(3)   + TVR(IR,3)
        TV2(3)  = TV2(3)  + TVR2(IR,3)
        TIZ(3)  = TIZ(3)  + TIZR(IR,3)
        TIZ2(3) = TIZ2(3) + TIZR2(IR,3)
        TRC(3)  = TRC(3)  + TRCR(IR,3)
        TRC2(3) = TRC2(3) + TRCR2(IR,3)
 10   CONTINUE
c
c     SOL
c
      TVS(3)   = TVS(1)   + TVS(2)
      TVS2(3)  = TVS2(1)  + TVS2(2)
      TIZS(3)  = TIZS(1)  + TIZS(2)
      TIZS2(3) = TIZS2(1) + TIZS2(2)
      TRCS(3)  = TRCS(1)  + TRCS(2)
      TRCS2(3) = TRCS2(1) + TRCS2(2)
c
c     Core
c
      TVC(3)   = TVC(1)   + TVC(2)
      TVC2(3)  = TVC2(1)  + TVC2(2)
      TIZC(3)  = TIZC(1)  + TIZC(2)
      TIZC2(3) = TIZC2(1) + TIZC2(2)  
      TRCC(3)  = TRCC(1)  + TRCC(2)
      TRCC2(3) = TRCC2(1) + TRCC2(2)
c
c     SOL ZXP Data          
c
      TVD(3)   = TVD(1)   + TVD(2)
      TVD2(3)  = TVD2(1)  + TVD2(2)
      TIZD(3)  = TIZD(1)  + TIZD(2)
      TIZD2(3) = TIZD2(1) + TIZD2(2)
      TRCD(3)  = TRCD(1)  + TRCD(2)
      TRCD2(3) = TRCD2(1) + TRCD2(2)
C
C     CALCULATE THE FLUXES TO THE PLATES FOR EACH RING USING
C     THE METHOD IN NEUT
C
      DO 30 ID = 1,NDS
        IK = IKDS(ID)
        IR = IRDS(ID)
        IN = 2
        IF (IK.GT.(NKS(IR)/2)) IN = 1
c        CS = 9.79E3 * SQRT (0.5*(KTEDS(ID)+KTIDS(ID))*
c     >                (1.0+RIZB)/CRMB)
c - sol option 22 or 23 allows mach numbers at the target > 1
c        if (cioptf.eq.22.or.cioptf.eq.23) cs = cs * cmachno(ir,in)
c
c
        cs = abs(kvds(id))
c
c       nonorth
c
        if (northopt.eq.0) then
           FLUXR(IR,IN)  = KNDS(ID) * CS / KBFST(IR,IN)
        elseif (northopt.eq.1.or.northopt.eq.2) then
           FLUXR(IR,IN)  = KNDS(ID) * CS / KBFST(IR,IN)*COSTET(ID)
        endif
c
c       nonorth
c
c        FLUXR(IR,IN)  = KNDS(ID) * CS / KBFST(IR,IN) * COSTET(ID)
c
        DISTR(IR,IN)  = DDS(ID)
        DISTR2(IR,IN) = DDS2(ID)
        SRCR(IR,IN)   = FLUXR(IR,IN) * DISTR(IR,IN)
        SRCR2(IR,IN)  = FLUXR(IR,IN) * DISTR2(IR,IN)
        DIST(IN)  = DIST(IN)  + DISTR(IR,IN)
        DIST2(IN) = DIST2(IN) + DISTR2(IR,IN)
        SRC(IN)   = SRC(IN)   + SRCR(IR,IN)
        SRC2(IN)  = SRC2(IN)  + SRCR2(IR,IN)
 30   CONTINUE
      DO 40 IR = 1,NRS
        SRCR(IR,3)  = SRCR(IR,1)  + SRCR(IR,2)
        SRCR2(IR,3) = SRCR2(IR,1) + SRCR2(IR,2)
 40   CONTINUE
      SRC(3)  = SRC(1)  + SRC(2)
      SRC2(3) = SRC2(1) + SRC2(2)
C
C     PRINT OUT A SUMMARY OF THE RESULTS FOUND
C
      WRITE(6,*) ' '
      WRITE(6,9000)
      WRITE(6,*) ' '
      WRITE(6,9005) INNER,OUTER
      WRITE(6,9010)
      DO 50 IR = 1,NRS
        WRITE(6,9015) IR,psitarg(ir,1),psitarg(ir,2),
     >    TVR(IR,1),TIZR(IR,1),TVR2(IR,1),TIZR2(IR,1),
     >    TVR(IR,2),TIZR(IR,2),TVR2(IR,2),TIZR2(IR,2),
     >    TVR(IR,3),TIZR(IR,3),TVR2(IR,3),TIZR2(IR,3)
 50   CONTINUE
      WRITE (6,*) ' '
      WRITE(6,9020) '   SOL',
     >  TVS(1),TIZS(1),TVS2(1),TIZS2(1),
     >  TVS(2),TIZS(2),TVS2(2),TIZS2(2),
     >  TVS(3),TIZS(3),TVS2(3),TIZS2(3)
      WRITE(6,9020) '  CORE',
     >  TVC(1),TIZC(1),TVC2(1),TIZC2(1),
     >  TVC(2),TIZC(2),TVC2(2),TIZC2(2),
     >  TVC(3),TIZC(3),TVC2(3),TIZC2(3)
      WRITE(6,9020) ' TOTAL',
     >  TV(1), TIZ(1), TV2(1), TIZ2(1),
     >  TV(2), TIZ(2), TV2(2), TIZ2(2),
     >  TV(3), TIZ(3), TV2(3), TIZ2(3)
      WRITE(6,*) ' '
      WRITE(6,9025)
      WRITE(6,9030)
      WRITE(6,9035)
      WRITE(6,9040)
     >  TVD(1), TIZD(1), TVD2(1), TIZD2(1),
     >  TVD(2), TIZD(2), TVD2(2), TIZD2(2),
     >  TVD(3), TIZD(3), TVD2(3), TIZD2(3)
      WRITE(6,*) ' '
      WRITE(6,9045)
      WRITE(6,9050) INNER,OUTER
      WRITE(6,9055)
      DO 60 IR = 1,NRS
        IF (IR.GE.IRSEP)
     >  WRITE(6,9060) IR,
     >    FLUXR(IR,1),DISTR(IR,1),DISTR2(IR,1),SRCR(IR,1),SRCR2(IR,1),
     >    FLUXR(IR,2),DISTR(IR,2),DISTR2(IR,2),SRCR(IR,2),SRCR2(IR,2),
     >    SRCR(IR,3),SRCR2(IR,3)
 60   CONTINUE
      WRITE (6,*) ' '
      WRITE(6,9065) ' TOTAL',
     >    DIST(1),DIST2(1),SRC(1),SRC2(1),
     >    DIST(2),DIST2(2),SRC(2),SRC2(2),
     >    SRC(3),SRC2(3)
      WRITE(6,*) ' '
      WRITE(6,9001)
      WRITE(6,*) ' '
      WRITE(6,9005) INNER,OUTER
      WRITE(6,9011)
      DO 70 IR = 1,NRS
c        WRITE(6,9015) IR,
        WRITE(6,9015) IR,psitarg(ir,1),psitarg(ir,2),
     >    TVR(IR,1),TRCR(IR,1),TVR2(IR,1),TRCR2(IR,1),
     >    TVR(IR,2),TRCR(IR,2),TVR2(IR,2),TRCR2(IR,2),
     >    TVR(IR,3),TRCR(IR,3),TVR2(IR,3),TRCR2(IR,3)
 70   CONTINUE
      WRITE (6,*) ' '
      WRITE(6,9020) '   SOL',
     >  TVS(1),TRCS(1),TVS2(1),TRCS2(1),
     >  TVS(2),TRCS(2),TVS2(2),TRCS2(2),
     >  TVS(3),TRCS(3),TVS2(3),TRCS2(3)
      WRITE(6,9020) '  CORE',
     >  TVC(1),TRCC(1),TVC2(1),TRCC2(1),
     >  TVC(2),TRCC(2),TVC2(2),TRCC2(2),
     >  TVC(3),TRCC(3),TVC2(3),TRCC2(3)
      WRITE(6,9020) ' TOTAL',
     >  TV(1), TRC(1), TV2(1), TRC2(1),
     >  TV(2), TRC(2), TV2(2), TRC2(2),
     >  TV(3), TRC(3), TV2(3), TRC2(3)
      WRITE(6,*) ' '
      WRITE(6,9026)
      WRITE(6,9030)
      WRITE(6,9036)
      WRITE(6,9040)
     >  TVD(1), TRCD(1), TVD2(1), TRCD2(1),
     >  TVD(2), TRCD(2), TVD2(2), TRCD2(2),
     >  TVD(3), TRCD(3), TVD2(3), TRCD2(3)
      WRITE(6,*) ' '
C
9000  FORMAT(' IONIZATION DATA ANALYSIS SUMMARY:'/
     >       '     AREAS, DISTANCES, AND INTEGRALS LABELED ''2'' ARE CAL
     >CULATED FROM THE CELL VOLUME')
9001  FORMAT(' RECOMBINATION DATA ANALYSIS SUMMARY:'/
     >       '     AREAS, DISTANCES, AND INTEGRALS LABELED ''2'' ARE CAL
     >CULATED FROM THE CELL VOLUME')
c9005  FORMAT('                        INNER
c     >      OUTER                                   TOTAL')
9005  FORMAT('                        ',a5,'
     >      ',a5,'                                   TOTAL')
9010  FORMAT('  RING   PSI1     PSI2     AREA     IONIZ      VOL      IO
     >NIZ2     AREA     IONIZ      VOL      IONIZ2     AREA     IONIZ   
     >   VOL      IONIZ2')
9011  FORMAT('  RING   PSI1     PSI2     AREA     RECOM      VOL      RE
     >COM2     AREA      RECOM      VOL      RECOM2     AREA     RECOM  
     >   VOL      RECOM2')
9015  FORMAT(I6,1P,14E10.3)
9020  FORMAT(A6,1P,12E10.3)
9025  FORMAT(' SOL IONIZATION:')
9026  FORMAT(' SOL RECOMBINATION:')
9030  FORMAT('                        < ZXP
     >      > ZXP                                   TOTAL')
9035  FORMAT('         AREA     IONIZ      VOL      IONIZ2     AREA
     >IONIZ      VOL      IONIZ2     AREA     IONIZ      VOL      IONIZ2
     > ')
9036  FORMAT('         AREA     RECOM      VOL      RECOM2     AREA
     >RECOM      VOL      RECOM2     AREA     RECOM      VOL      RECOM2
     > ')
9040  FORMAT(6X,1P,12E10.3)
9045  FORMAT(' SOURCE DATA ANALYSIS SUMMARY:')
c9050  FORMAT('                                 INNER
c     >                         OUTER                    TOTAL     TOTAL
c     >')
9050  FORMAT('                                 ',a5,'
     >                         ',a5,'                    TOTAL     TOTAL
     >')
9055  FORMAT('  RING   FLUX      DIST      DIST2    SOURCE    SOURCE2
     > FLUX      DIST      DIST2    SOURCE    SOURCE2   SOURCE    SOURCE
     >2')
9060  FORMAT(I6,1P,12E10.3)
9065  FORMAT(A6,1P,10X,4E10.3,10X,6E10.3)
C
      RETURN
      END
C                                                                       
C                                                                       
      SUBROUTINE READPIN                                                
      IMPLICIT NONE                                                     
C     INCLUDE "PARAMS"                                                  
      include 'params'                                                  
C     INCLUDE "PINDATA"                                                 
      include 'pindata'                                                 
C     INCLUDE "COMTOR"                                                  
      include 'comtor'                                                  
C     INCLUDE "CGEOM"                                                   
      include 'cgeom'                                                   
c slmod begin
      include 'slcom'
c slmod end
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C     READPIN: THIS SUBROUTINE IS DESIGNED TO READ IN THE               
C     DATA REQUIRED FROM SPECIFIC PIN OUTPUT FILES IN THE               
C     SPECIFIC DATA FORMAT USED BY THOSE FILES.                         
C                                                                       
C     IMPURITY WALL SOURCES DUE TO NEUTRAL HYDROGEN SPUTTERING          
C     ARE NOT PASSED BACK TO PIN (NOR IS THE WALL GEOMETRY, FOR         
C     THAT MATTER) SO THIS DATA MUST BE READ FROM THE NIMBUS            
C     PRINT FILE.  **NB** CARE MUST BE TAKEN WHEN IMPLEMENTING          
C     NEW VERSIONS OF NIMBUS TO GUARD AGAINST CHANGES IN THE FORMAT     
C     OF THE PRINT FILE.                                                
C                                                                       
C     DATA WHICH IS PASSED BACK TO PIN IS NOW WRITTEN IN A LABELED      
C     PASSING FILE BY THE GHOST SUBROUTINE IN PINPGX.  THIS MAKES       
C     ADDING ADDITIONAL INFORMATION MORE TRANSPARENT.  THIS DATA IS     
C     THEN LOADED INTO THE DIVIMP MATRICES USING THE VERSION OF KORY    
C     FOUND IN CGEOM.  IN THIS WAY, IF THE VIRTUAL POINTS HAVE BEEN     
C     REMOVED FROM THE OPEN RINGS IN TAU, APPROPRIATE MAPPINGS ARE      
C     MADE WITH THE MODIFIED KORY.                                      
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C     LOCAL VARIABLES                                                   
C                                                                       
      INTEGER I,J,K,ZONE,IND                                              
      CHARACTER*160 LINE                                                
      CHARACTER*4 INTER                                                 
      CHARACTER*10 TEMPSTR                                              
      REAL INFO(10)                                                     
      REAL TEMPDAT(MAXNKS*MAXNRS)                                       
      INTEGER IN,IR,IK,IND1,NP,IW
      INTEGER CNT,SPUTCNT,RES
      INTEGER NTRAJ, NSTEP(200)
      REAL CWSPUT(MAXPTS,2)
      INTEGER FINDPOL                                                   
      EXTERNAL FINDPOL                                                  
c
      integer lenstr,len
      external lenstr
C                                                                       
C     VARIABLES FOR ASSIGNING C-SPUTTERING                              
C                                                                       
      INTEGER NWLPROB2                                                  
      REAL    WLPROB2(MAXPTS,3)                                         
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C     READ IN CX-REC SPUTTERING DATA                                    
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      REWIND(INIMOUT)                                                   
 10   READ(INIMOUT,1000,end=40,err=40) LINE                             
C                                                                       
C      WRITE(6,*) 'CX:',LINE(1:27)                                      
C                                                                       
      IF (LINE(1:27).NE.' INTERACTIONS WITH SURFACES') GOTO 10          
C                                                                       
       WRITE(6,*) 'READING CX'                                          
C                                                                       
      READ(INIMOUT,1000) LINE                                           
C                                                                       
C     READ TABLE AND INDEX THE RESULTS BY POL NUMBER - THESE WILL       
C     BE USED TO ASSIGN THE RESULTS TO THE APPROPRIATE WALL SEGMENTS    
C                                                                       
 20   READ(INIMOUT,1000,end=40,err=40) LINE                             
      IF ((LINE(1:5).EQ.'     ').OR.(LINE(1:5).EQ.' ZONE')) GOTO 20     
C                                                                       
      CNT = 0                                                           
C                                                                       
C     READ FIRST LINE                                                   
C                                                                       
 30   READ(LINE,3000,END=40,ERR=40) ZONE,INTER,(INFO(I),I=1,7)          
      IF (ZONE.EQ.0) GOTO 40                                            
C                                                                       
C     READ THE REST UNTIL CARBON                                        
C                                                                       
 50   READ(INIMOUT,1000,END=40) LINE                                    
      READ(LINE,4000) INTER,(INFO(I),I=1,7)                             
      IF (INTER.EQ.'C   ') THEN                                         
         CNT = CNT +1                                                   
         CWSPUT(CNT,1) = ZONE                                           
         CWSPUT(CNT,2) = INFO(3)                                        
      ELSEIF (INTER.EQ.'    ') THEN                                     
         READ(INIMOUT,1000,END=40) LINE                                 
         GOTO 30                                                        
      ENDIF                                                             
      GOTO 50                                                           
C                                                                       
C     THE ZONE NUMBER SHOULD BE ZERO ON A BLANK LINE READ - OR          
C     INVALID DATA. THIS SHOULD CAUSE A BRANCH BEYOND THE SECTION       
C     READING IN THE TABLE                                              
C                                                                       
 40   CONTINUE                                                          
C                                                                       
       WRITE(6,*) 'CX CNT:',CNT                                         
C                                                                       
      SPUTCNT = CNT                                                     
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C     READ IN DATA FROM PIN PASSING FILE                                
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C     Zero the arrays                                                   
C                                                                       
      call rzero (hcorr,maxnks*maxnrs)                                
      call rzero (hval,maxnks*maxnrs)                                
      call rzero (pinatom,maxnks*maxnrs)                                
      call rzero (pinion,maxnks*maxnrs)                                 
      call rzero (pinalpha,maxnks*maxnrs)                               
      call rzero (pinmol,maxnks*maxnrs)                                 
      call rzero (pinz0,maxnks*maxnrs)                                  
      call rzero (pinionz,maxnks*maxnrs)                                
      call rzero (pinena,maxnks*maxnrs)                                 
      call rzero (pinenm,maxnks*maxnrs)                                 
      call rzero (pinenz,maxnks*maxnrs)                                 
      call rzero (pinqi,maxnks*maxnrs)                                  
      call rzero (pinqe,maxnks*maxnrs)                                  
      call rzero (pinmp,maxnks*maxnrs)                                  
      call rzero (pinvdist,3*14*maxnks*maxnrs)                          
      call rzero (rvesm,2*maxseg)
      call rzero (zvesm,2*maxseg)
      call rzero (fluxhw,maxseg)
      call rzero (hwalks,maxnws*2)
c
      h2ioniz = 0.0
C                                                                       
      REWIND(IPINOUT)                                                   
C                                                                       
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C     READ IN THE PIN DATA. USING THE KORY FROM CGEOM TO LOAD THEM      
C     INTO DIVIMP FORMAT ARRAYS.  NOTE ALSO THAT NIMBUS WORKS IN CGS    
C     SO CONVERSION TO MKS IS NECESSARY.                                
C                                                                       
C-----------------------------------------------------------------------
c                                                                       
 200  READ (IPINOUT,1010,END=400) LINE
      IF (LINE(3:9).EQ.'SRECYC:') THEN                               
C                                                                       
C     TOTAL RECYCLING NEUTRAL SOURCE (PER M TOROIDALLY, PER S)  
C
        READ(IPINOUT,9010) TEMPDAT(1)
        SRECYC = TEMPDAT(1) * 1.0E2
      ELSE IF (LINE(3:9).EQ.'SRECOM:') THEN                               
C                                                                       
C     TOTAL RECOMBINATION NEUTRAL SOURCE (PER M TOROIDALLY, PER S)  
C
        READ(IPINOUT,9010) TEMPDAT(1)
        SRECOM = TEMPDAT(1) * 1.0E2
      ELSE IF (LINE(3:9).EQ.'HESCPD:') THEN                               
C                                                                       
C     NUMBER OF NEUTRALS LOST FROM THE NEUTRAL CALCULATION EITHER TO 
C     THE CORE OR AN ALBEDO REGION
C
        READ(IPINOUT,9010) TEMPDAT(1)
        HESCPD = TEMPDAT(1) * 1.0E2
      ELSE IF (LINE(3:9).EQ.'HESCAL:') THEN                               
C                                                                       
C     NUMBER OF NEUTRALS LOST TO ALBEDO REGIONS (USUALLY THE PUMP)
C
        READ(IPINOUT,9010) TEMPDAT(1)
        HESCAL = TEMPDAT(1) * 1.0E2
      ELSE IF (LINE(3:9).EQ.'ZSPUT :') THEN                               
C                                                                       
C     NUMBER OF IMPURITIES SPUTTERED, LESS Z+0 REDEPOSITION
C
        READ(IPINOUT,9010) TEMPDAT(1)
        ZSPUT  = TEMPDAT(1) * 1.0E2
      ELSE IF (LINE(3:9).EQ.'ZSPUTN:') THEN                               
C                                                                       
C     NUMBER OF IMPURITIES SPUTTERED BY NEUTRALS,
C     LESS !!TOTAL!! Z+0 REDEPOSITION
C
        READ(IPINOUT,9010) TEMPDAT(1)
        ZSPUTN = TEMPDAT(1) * 1.0E2
      ELSE IF (LINE(3:9).EQ.'ZESCPD:') THEN                               
C                                                                       
C     NUMBER OF IMPURITIES LOST TO CENTRAL ESCAPE REGION
C
        READ(IPINOUT,9010) TEMPDAT(1)
        ZESCPD = TEMPDAT(1) * 1.0E2
      ELSE IF (LINE(3:9).EQ.'ZESCAL:') THEN                               
C                                                                       
C     NUMBER OF IMPURITIES LOST TO ALBEDO REGIONS (USUALLY THE PUMP)
C
        READ(IPINOUT,9010) TEMPDAT(1)
        ZESCAL = TEMPDAT(1) * 1.0E2
      ELSE IF (LINE(3:9).EQ.'HESCLK:') THEN                               
C                                                                       
C     NUMBER OF NEUTRALS LOST TO LEAK REGIONS - SUBSEQUENTLY REINJECTED
C     INSIDE NIMBUS FROM THE LEAK REINJECTION SEGMENTS.
C
        READ(IPINOUT,9010) TEMPDAT(1)
        HESCLK = TEMPDAT(1) * 1.0E2
      ELSE IF (LINE(3:9).EQ.'ZESCLK:') THEN                               
C                                                                       
C     NUMBER OF IMPURITIES LOST TO LEAK REGIONS
C
        READ(IPINOUT,9010) TEMPDAT(1)
        ZESCLK = TEMPDAT(1) * 1.0E2
      ELSE IF (LINE(3:9).EQ.'PHFGAL:') THEN                               
C                                                                       
C     FRACTION OF NEUTRALS LOST INTO ALBEDO REGIONS
C
        READ(IPINOUT,9010) TEMPDAT(1)
        PHFGAL = TEMPDAT(1)
      ELSE IF (LINE(3:9).EQ.'PHFUGA:') THEN                               
C                                                                       
C     FRACTION OF NEUTRALS LOST INTO REGIONS OTHER THAN ALBEDOS
C
        READ(IPINOUT,9010) TEMPDAT(1)
        PHFUGA = TEMPDAT(1)
      ELSE IF (LINE(4:9).EQ.'PROFA:') THEN                               
C                                                                       
C     NEUTRAL ATOM DENSITY                                              
C                                                                       
        READ(LINE,9000) TEMPSTR,NP                                      
        READ(IPINOUT,9010) (TEMPDAT(IN),IN=1,NP)                        
        DO 210 IR = 1,NRS                                               
          DO 210 IK = 1,NKS(IR)                                         
            PINATOM(IK,IR) = TEMPDAT(KORY(IR,IK)) * 1.0E6               
 210    CONTINUE                                                        
      ELSE IF (LINE(3:9).EQ.'PROFSN:') THEN                              
C                                                                       
C     IONIZATION SOURCE                                                 
C                                                                       
        READ(LINE,9000) TEMPSTR,NP                                      
        READ(IPINOUT,9010) (TEMPDAT(IN),IN=1,NP)                        
C
C     CALCULATE ALSO THE TOTAL IONISATION SOURCE
C
        HIONIZ = 0.0
        DO IR = 1,NRS                                               
          DO IK = 1,NKS(IR)                                         
            PINION(IK,IR) = TEMPDAT(KORY(IR,IK)) * 1.0E6                
c
c           Exclude the repeat point in the core from the
c           sum of the ionization. 
c
            if (.not.(ir.lt.irsep.and.ik.eq.nks(ir)))  
     >         HIONIZ = HIONIZ + PINION(IK,IR)*KAREA2(IK,IR)               
c
          end do
c
          if (ir.lt.irsep.and.
     >             pinion(1,ir).ne.(pinion(nks(ir),ir))) then 

             write(6,*) 'POSSIBLE ERROR: PINION - '//
     >                'FIRST AND LAST IN CORE NOT EQUAL',ir,
     >                  pinion(1,ir),pinion(nks(ir),ir)

          endif
c
        end do
c
      ELSE IF (LINE(3:9).EQ.'PROFHA:') THEN                             
C                                                                       
C     HALPHA EMISSION                                                   
C                                                                       
        READ(LINE,9000) TEMPSTR,NP                                      
        READ(IPINOUT,9010) (TEMPDAT(IN),IN=1,NP)                        
        DO 230 IR = 1,NRS                                               
          DO 230 IK = 1,NKS(IR)                                         
            PINALPHA(IK,IR) = TEMPDAT(KORY(IR,IK)) * 1.0E6              
 230    CONTINUE                                                        
      ELSE IF (LINE(4:9).EQ.'PROFM:') THEN                              
C                                                                       
C     MOLECULAR DENSITY                                                 
C                                                                       
        READ(LINE,9000) TEMPSTR,NP                                      
        READ(IPINOUT,9010) (TEMPDAT(IN),IN=1,NP)                        
        DO 240 IR = 1,NRS                                               
          DO 240 IK = 1,NKS(IR)                                         
            PINMOL(IK,IR) = TEMPDAT(KORY(IR,IK)) * 1.0E6                
 240    CONTINUE                                                        
      ELSE IF (LINE(4:9).EQ.'PROFZ:') THEN                              
C                                                                       
C     NEUTRAL IMPURITY DENSITY                                          
C                                                                       
        READ(LINE,9000) TEMPSTR,NP                                      
        READ(IPINOUT,9010) (TEMPDAT(IN),IN=1,NP)                        
        DO 250 IR = 1,NRS                                               
          DO 250 IK = 1,NKS(IR)                                         
            PINZ0(IK,IR) = TEMPDAT(KORY(IR,IK)) * 1.0E6                 
 250    CONTINUE                                                        
      ELSE IF (LINE(3:9).EQ.'PROFSZ:') THEN                             
C                                                                       
C     IMPURITY IONIZATION SOURCE                                        
C                                                                       
        READ(LINE,9000) TEMPSTR,NP                                      
        READ(IPINOUT,9010) (TEMPDAT(IN),IN=1,NP)                        
C
C     CALCULATE ALSO THE TOTAL IMPURITY IONISATION SOURCE
C
        ZIONIZ = 0.0
        DO 260 IR = 1,NRS                                               
          DO 260 IK = 1,NKS(IR)                                         
            PINIONZ(IK,IR) = TEMPDAT(KORY(IR,IK)) * 1.0E6
            ZIONIZ = ZIONIZ + PINIONZ(IK,IR)*KAREA2(IK,IR)               
 260    CONTINUE
        write (6,*) 'PIN:Zioniz',zioniz
      ELSE IF (LINE(3:9).EQ.'ENEUTA:') THEN                             
C                                                                       
C     AVERAGE ATOMIC HYDROGEN ENERGY (EV)                               
C                                                                       
        READ(LINE,9000) TEMPSTR,NP                                      
        READ(IPINOUT,9010) (TEMPDAT(IN),IN=1,NP)                        
        DO 270 IR = 1,NRS                                               
          DO 270 IK = 1,NKS(IR)                                         
            PINENA(IK,IR) = TEMPDAT(KORY(IR,IK)) / 1.6E-12             
 270    CONTINUE                                                        
      ELSE IF (LINE(3:9).EQ.'ENEUTM:') THEN                             
C                                                                       
C     AVERAGE MOLECULAR HYDROGEN ENERGY (EV)                            
C                                                                       
        READ(LINE,9000) TEMPSTR,NP                                      
        READ(IPINOUT,9010) (TEMPDAT(IN),IN=1,NP)                        
        DO 280 IR = 1,NRS                                               
          DO 280 IK = 1,NKS(IR)                                         
            PINENM(IK,IR) = TEMPDAT(KORY(IR,IK)) / 1.6E-12              
 280    CONTINUE                                                        
      ELSE IF (LINE(3:9).EQ.'ENEUTZ:') THEN                             
C                                                                       
C     AVERAGE NEUTRAL IMPURITY ENERGY (EV)                              
C                                                                       
        READ(LINE,9000) TEMPSTR,NP                                      
        READ(IPINOUT,9010) (TEMPDAT(IN),IN=1,NP)                        
        DO 290 IR = 1,NRS                                               
          DO 290 IK = 1,NKS(IR)                                         
            PINENZ(IK,IR) = TEMPDAT(KORY(IR,IK)) / 1.6E-12             
 290    CONTINUE                                                        
      ELSE IF (LINE(4:9).EQ.'PROFQ:') THEN                             
C                                                                       
C     ION ENERGY SOURCE TERM                                            
C                                                                       
        READ(LINE,9000) TEMPSTR,NP                                      
        READ(IPINOUT,9010) (TEMPDAT(IN),IN=1,NP)                        
        DO IR = 1,NRS                                                   
          DO IK = 1,NKS(IR)                                             
            PINQI(IK,IR) = TEMPDAT(KORY(IR,IK)) * 1.0E-1                
          ENDDO                                                         
        ENDDO                                                           
      ELSE IF (LINE(3:9).EQ.'PROFQE:') THEN                             
C                                                                       
C     ELECTRON ENERGY SOURCE TERM                                            
C                                                                       
        READ(LINE,9000) TEMPSTR,NP                                      
        READ(IPINOUT,9010) (TEMPDAT(IN),IN=1,NP)                        
        DO IR = 1,NRS                                                   
          DO IK = 1,NKS(IR)                                             
            PINQE(IK,IR) = TEMPDAT(KORY(IR,IK)) * 1.0E-1                
          ENDDO                                                         
        ENDDO                                                           
      ELSE IF (LINE(3:9).EQ.'PROFMP:') THEN                             
C                                                                       
C     ION MOMENTUM SOURCE TERM                                            
C                                                                       
        READ(LINE,9000) TEMPSTR,NP                                      
        READ(IPINOUT,9010) (TEMPDAT(IN),IN=1,NP)                        
        DO IR = 1,NRS                                                   
          DO IK = 1,NKS(IR)                                             
            PINMP(IK,IR) = TEMPDAT(KORY(IR,IK)) * 1.0E1                 
          ENDDO                                                         
        ENDDO                                                           
      ELSE IF (LINE(1:5).EQ.'VDIST') THEN                               
C                                                                       
C     NEUTRAL PARTICLE VELOCITY DISTRIBUTION INFORMATION FROM NIMBUS    
C                                                                       
        READ(LINE,9020) I,J,NP                                          
        READ(IPINOUT,9010) (TEMPDAT(IN),IN=1,NP)                        
        DO 300 IR = 1,NRS                                               
          DO 300 IK = 1,NKS(IR)                                         
            IF (J.EQ.1  .OR. J.EQ.2  .OR. J.EQ.3 .OR.                   
     >          J.EQ.13 .OR. J.EQ.14) THEN                              
              PINVDIST(I,J,IK,IR) = TEMPDAT(KORY(IR,IK)) * 1.0E6        
            ELSE IF (J.GE.4 .AND. J.LE.6) THEN                          
              PINVDIST(I,J,IK,IR) = TEMPDAT(KORY(IR,IK)) / 1.0E2        
            ELSE IF (J.GE.7 .AND. J.LE.12) THEN                         
              PINVDIST(I,J,IK,IR) = TEMPDAT(KORY(IR,IK)) / 1.0E4        
            ENDIF                                                       
        
            if (i.eq.1.and.j.eq.13.and.
     >         (.not.(ir.lt.irsep.and.ik.eq.nks(ir))))
     >         H2IONIZ = H2IONIZ+PINVDIST(i,j,ik,ir)*KAREA2(IK,IR)               


 300    CONTINUE                                                        
c
      ELSE IF (LINE(4:9).EQ.'NVESM:') THEN                             
C                                                                       
C     NUMBER OF SEGMENTS IN NIMBUS WALL  
C                                                                       
        READ(LINE,9000) TEMPSTR,NVESM
      ELSE IF (LINE(4:9).EQ.'NVESP:') THEN                             
C                                                                       
C     NUMBER OF SEGMENTS IN NIMBUS PUMP
C                                                                       
        READ(LINE,9000) TEMPSTR,NVESP                                 
      ELSE IF (LINE(3:9).EQ.'RVESM1:') THEN                             
C                                                                       
C     MAJOR RADIUS OF FIRST POINT IN WALL/PUMP SEGMENT
C                                                                       
        READ(LINE,9000) TEMPSTR,NP                                      
        READ(IPINOUT,9010) (TEMPDAT(IN),IN=1,NP)                        
        DO IW = 1,NVESM+NVESP
          RVESM(IW,1) = TEMPDAT(IW) * 0.01
        ENDDO                                                           
      ELSE IF (LINE(3:9).EQ.'ZVESM1:') THEN                             
C                                                                       
C     HEIGHT OF FIRST POINT IN WALL/PUMP SEGMENT
C                                                                       
        READ(LINE,9000) TEMPSTR,NP                                      
        READ(IPINOUT,9010) (TEMPDAT(IN),IN=1,NP)                        
        DO IW = 1,NVESM+NVESP
          ZVESM(IW,1) = TEMPDAT(IW) * 0.01
        ENDDO                                                           
      ELSE IF (LINE(3:9).EQ.'RVESM2:') THEN                             
C                                                                       
C     MAJOR RADIUS OF SECOND POINT IN WALL/PUMP SEGMENT
C                                                                       
        READ(LINE,9000) TEMPSTR,NP                                      
        READ(IPINOUT,9010) (TEMPDAT(IN),IN=1,NP)                        
        DO IW = 1,NVESM+NVESP
          RVESM(IW,2) = TEMPDAT(IW) * 0.01
        ENDDO                                                           
      ELSE IF (LINE(3:9).EQ.'ZVESM2:') THEN                             
C                                                                       
C     HEIGHT OF SECOND POINT IN WALL/PUMP SEGMENT
C                                                                       
        READ(LINE,9000) TEMPSTR,NP                                      
        READ(IPINOUT,9010) (TEMPDAT(IN),IN=1,NP)                        
        DO IW = 1,NVESM+NVESP
          ZVESM(IW,2) = TEMPDAT(IW) * 0.01
        ENDDO
      ELSE IF (LINE(4:9).EQ.'JVESM:') THEN                             
C                                                                       
C     LABEL FOR EACH WALL SEGMENT
C                                                                       
        READ(LINE,9000) TEMPSTR,NP                                      
        READ(IPINOUT,9030) (JVESM(IN),IN=1,NP)                        
      ELSE IF (LINE(3:9).EQ.'FLUXHW:') THEN                             
C                                                                       
C     FLUX OF HYDROGEN (ATOMS AND MOLECULES) TO THE WALL
C                                                                       
        READ(LINE,9000) TEMPSTR,NP
        READ(IPINOUT,9010) (TEMPDAT(IN),IN=1,NP)
        DO IW = 1,NVESM
          FLUXHW(IW) = TEMPDAT(IW) * 1.0E4
        ENDDO
      ELSE IF (LINE(3:9).EQ.'FLXHW2:') THEN                             
C                                                                       
C     FLUX OF HYDROGEN (ATOMS AND IONS) TO THE WALL
C                                                                       
        READ(LINE,9000) TEMPSTR,NP
        READ(IPINOUT,9010) (TEMPDAT(IN),IN=1,NP)
        DO IW = 1,NVESM
          FLXHW2(IW) = TEMPDAT(IW) * 1.0E4
        ENDDO
      ELSE IF (LINE(3:9).EQ.'FLXHW3:') THEN                             
C                                                                       
C     FLUX OF IMPURITIES SPUTTERED FROM THE WALL
C                                                                       
        READ(LINE,9000) TEMPSTR,NP
        READ(IPINOUT,9010) (TEMPDAT(IN),IN=1,NP)
        DO IW = 1,NVESM
          FLXHW3(IW) = TEMPDAT(IW) * 1.0E4
        ENDDO
      ELSE IF (LINE(3:9).EQ.'FLXHW4:') THEN                             
C                                                                       
C     FLUX OF IMPURITIES REDEPOSITED ONTO THE WALL
C                                                                       
        READ(LINE,9000) TEMPSTR,NP
        READ(IPINOUT,9010) (TEMPDAT(IN),IN=1,NP)
        DO IW = 1,NVESM
          FLXHW4(IW) = TEMPDAT(IW) * 1.0E4
        ENDDO
      ELSE IF (LINE(3:9).EQ.'FLXHW5:') THEN                             
C                                                                       
C     AVERAGE ENERGY OF ATOMS HITTING THE WALL (EV)
C                                                                       
        READ(LINE,9000) TEMPSTR,NP
        READ(IPINOUT,9010) (TEMPDAT(IN),IN=1,NP)
        DO IW = 1,NVESM
          FLXHW5(IW) = TEMPDAT(IW)
        ENDDO
      ELSE IF (LINE(3:9).EQ.'FLXHW6:') THEN                             
C                                                                       
C     FLUX OF HYDROGEN ATOMS TO THE WALL
C                                                                       
        READ(LINE,9000) TEMPSTR,NP
        READ(IPINOUT,9010) (TEMPDAT(IN),IN=1,NP)
        DO IW = 1,NVESM
          FLXHW6(IW) = TEMPDAT(IW) * 1.0E4
        ENDDO
      ELSE IF (LINE(3:9).EQ.'KDEBWR:') THEN                             
C                                                                       
C     NUMBER OF TRAJECTORIES SAVED IN NIMBUS
C      - HERE I'VE BROKEN MY RULE ON ORDERING -- THE
C        PARAMETERS NTRAJ AND NSTEP MUST BE PASSED
C        BEFORE THE COORDINATES OF THE TRAJECTORIES
C                                                                       
        READ(LINE,9000) TEMPSTR, NTRAJ
        READ(IPINOUT,9030) (NSTEP(IN),IN=1,NTRAJ)
      ELSE IF (LINE(2:9).EQ.'XDEBWR1:') THEN                             
C                                                                       
C     R COORDINATE OF TRAJECTORIES
C                                                                       
        READ(LINE,9000) TEMPSTR,NP
        READ(IPINOUT,9010) (TEMPDAT(IN),IN=1,NP)
        IW = 0
        IN = 0
        DO I = 1, NTRAJ
          DO K = 1, NSTEP(I)
            IF (IW.LT.MAXNWS-1) THEN
              IW = IW + 1
              IN = IN + 1
              HWALKS(IW,1) = TEMPDAT(IN) * 0.01
            ENDIF
          ENDDO
          IF (IW.LT.MAXNWS) THEN
            IW = IW + 1
            HWALKS(IW,1) = HI
          ENDIF
        ENDDO
      ELSE IF (LINE(2:9).EQ.'XDEBWR2:') THEN                             
C                                                                       
C     Z COORDINATE OF TRAJECTORIES
C                                                                       
        READ(LINE,9000) TEMPSTR,NP
        READ(IPINOUT,9010) (TEMPDAT(IN),IN=1,NP)
        IW = 0
        IN = 0
        DO I = 1, NTRAJ
          DO K = 1, NSTEP(I)
            IF (IW.LT.MAXNWS-1) THEN
              IW = IW + 1
              IN = IN + 1
              HWALKS(IW,2) = TEMPDAT(IN) * 0.01
            ENDIF
          ENDDO
          IF (IW.LT.MAXNWS) THEN
            IW = IW + 1
            HWALKS(IW,2) = HI
          ENDIF
        ENDDO
      ELSE IF (LINE(4:9).EQ.'HCORR:') THEN                              
C                                                                       
C     VOLUME CORRECTION FACTOR                                              
C                                                                       
        READ(LINE,9000) TEMPSTR,NP                                      
        READ(IPINOUT,9010) (TEMPDAT(IN),IN=1,NP)                        
        DO  IR = 1,NRS                                               
          DO IK = 1,NKS(IR)                                         
            HCORR(IK,IR) = TEMPDAT(KORY(IR,IK))
          end do
        end do
      ELSE IF (LINE(4:9).EQ.'    H:') THEN                              
C                                                                       
C     VOLUME CORRECTION FACTOR                                              
C                                                                       
        READ(LINE,9000) TEMPSTR,NP                                      
        READ(IPINOUT,9010) (TEMPDAT(IN),IN=1,NP)                        
        DO  IR = 1,NRS                                               
          DO  IK = 1,NKS(IR)                                         
            HVAL(IK,IR) = TEMPDAT(KORY(IR,IK))
          end do
        end do
      ENDIF
      GOTO 200                                     
c                                                                       
c                                                                       
 400  continue
c
c
c     End of Read of PIN file
c
c--------------------------------------------------------------
c
c     Read in GAUGE pressures from the .pinout file which is 
c     connected to unit IPINDAT - 15 at this time - need to 
c     check the value of LOUT in the PIN "plun" common block. 
c
      call rzero(gaugedat,5*2)
c     
      rewind(ipindat) 

 500  READ (IPINdat,1000,END=600) LINE


      IF (LINE(7:16).EQ.'GAUGE.PUMP') THEN                               
         len = lenstr(line) 
         read(line(17:len),*) (tempdat(in),in=1,10)
         gaugedat(1,1) = tempdat(1)
         gaugedat(2,1) = tempdat(2)
         gaugedat(3,1) = tempdat(7)
         gaugedat(4,1) = tempdat(8)
         write(6,*) 'GAUGE.PUMP:',(gaugedat(in,1),in=1,4)
      ELSEIF (LINE(7:15).EQ.'GAUGE.VES') THEN                               
         len = lenstr(line) 
         read(line(17:len),*) (tempdat(in),in=1,10)
         gaugedat(1,2) = tempdat(1)
         gaugedat(2,2) = tempdat(2)
         gaugedat(3,2) = tempdat(7)
         gaugedat(4,2) = tempdat(8)
         write(6,*) 'GAUGE.VES :',(gaugedat(in,1),in=1,4)
      endif

      goto 500

 600  continue

c
c     End of read on supplementary PIN file
c 
c--------------------------------------------------------------
c
c     Assign the calculated DIVIMP recombination to the PINREC quantity 
c     because NIMBUS does not return this. 
c
      do ir = 1,nrs
         do ik = 1,nks(ir)
            pinrec(ik,ir) = divrec(ik,ir)
         enddo
      enddo
c
c     Call targflux to calculate the correction factor, pincor,
c     necessary to conserve particles
c        
      call targflux
c
c     Apply the PIN correction factor necessary to guarantee 
c     particle conservation
c
        DO IR = 1,NRS 
          DO IK = 1,NKS(IR)                                         
            PINATOM(IK,IR) = PINATOM(IK,IR) * PINCOR
            PINION(IK,IR) = PINION(IK,IR) * PINCOR
            PINALPHA(IK,IR) = PINALPHA(IK,IR) * PINCOR
            PINMOL(IK,IR) = PINMOL(IK,IR) * PINCOR
            PINZ0(IK,IR) = PINZ0(IK,IR) * PINCOR
            PINIONZ(IK,IR) = PINIONZ(IK,IR) * PINCOR
            PINQI(IK,IR) = PINQI(IK,IR) * PINCOR
            PINQE(IK,IR) = PINQE(IK,IR) * PINCOR
            PINMP(IK,IR) = PINMP(IK,IR) * PINCOR
          ENDDO
        ENDDO
c
c       Debug - display PINQE/PINION ratio
c
c        do ir = 1,nrs
c           do ik = 1,nks(ir)
c              if (pinion(ik,ir).ne.0.0) then 
c                 write (6,*) ik,ir,pinqe(ik,ir),pinion(ik,ir),
c     >                    pinqe(ik,ir)/pinion(ik,ir),  
c     >                    pinqe(ik,ir)/pinion(ik,ir)/1.6e-19,
c     >                    ktebs(ik,ir),pinmol(ik,ir)
c              endif
c           end do
c        end do  
c
c
c       COPY all wall related quantities into arrays to store 
c       the original PIN data since the arrays that are originally
c       loaded may be changed or over-written if the wall is 
c       redefined. 
c
c       Correct wall fluxes - flxhw5 is not adjusted because it
c                             contains avarege temperatures and not
c                             fluxes.
c
        nvesm_pin = nvesm 
c
        do iw = 1,nvesm+nvesp
c
           fluxhw(iw) = fluxhw(iw) * pincor
           flxhw2(iw) = flxhw2(iw) * pincor
           flxhw3(iw) = flxhw3(iw) * pincor
           flxhw4(iw) = flxhw4(iw) * pincor
           flxhw6(iw) = flxhw6(iw) * pincor
c
c          Copy wall flux data
c
           fluxhw_pin(iw) = fluxhw(iw)
           flxhw2_pin(iw) = flxhw2(iw)
           flxhw3_pin(iw) = flxhw3(iw)
           flxhw4_pin(iw) = flxhw4(iw)
           flxhw5_pin(iw) = flxhw5(iw)
           flxhw6_pin(iw) = flxhw6(iw)
c
c          Copy vessel definitions  
c
           jvesm_pin(iw)  = jvesm(iw)
c
           do in = 1,2 
              rvesm_pin(iw,in)  = rvesm(iw,in)
              zvesm_pin(iw,in)  = zvesm(iw,in)
           end do

        end do 
c
        HIONIZ = HIONIZ * PINCOR
        H2ioniz = h2ioniz * pincor
        ZIONIZ = ZIONIZ * PINCOR
c
      write (6,*) 'srecO:',srecyc,hescpd,srecom 
c
c     Redo the walls - if the NIMBUS wall option is ON  
c
      if (cprint.eq.3.or.cprint.eq.9) then 
c
         write (6,*) 'Nves:',nvesm,nvesp,np
         do ik = 1,nvesm  
            write (6,'(4(1x,f18.10),i4)') rvesm(ik,1),
     >         zvesm(ik,1),
     >         rvesm(ik,2),zvesm(ik,2),jvesm(ik)
         end do 
c
         write (6,*) 'Nwall:',wallpts,nwall
         do ik = 1,wallpts
            write (6,'(4(1x,f18.10),i4)') 
     >        wallpt(ik,1) + wallpt(ik,5) * cos(wallpt(ik,8)),    
     >        wallpt(ik,2) + wallpt(ik,5) * sin(wallpt(ik,8)),    
     >        wallpt(ik,1) + wallpt(ik,6) * cos(wallpt(ik,9)),    
     >        wallpt(ik,2) + wallpt(ik,6) * sin(wallpt(ik,9)),
     >        ik 
         end do  
       
c
         write(6,*) 'NwallA:',wallpts     
         do ik = 1,wallpts
            write (6,'(2(1x,f18.10),i4)') 
     >      wallpt(ik,1),
     >      wallpt(ik,2),ik
         end do  
c
c slmod begin
         IF (grdnmod.NE.0) THEN
           WRITE(0,*) 'WARNING: STANDARD WALL INDICIES DO NOT APPLY'
         ENDIF
c slmod end
         write (6,*) 'Nwall:',wlwall2,wlwall1
         do ik = wlwall1,wlwall2
            write (6,'(4(1x,f18.10),i4)') 
     >        wallpt(ik,1) + wallpt(ik,5) * cos(wallpt(ik,8)),    
     >        wallpt(ik,2) + wallpt(ik,5) * sin(wallpt(ik,8)),    
     >        wallpt(ik,1) + wallpt(ik,6) * cos(wallpt(ik,9)),    
     >        wallpt(ik,2) + wallpt(ik,6) * sin(wallpt(ik,9)),
     >        ik
         end do  
c
         write (6,*) 'NInner:',wlwall2+1,wltrap1-1
         do ik = wlwall2+1,wltrap1-1
            write (6,'(4(1x,f18.10),i4)') 
     >        wallpt(ik,1) + wallpt(ik,5) * cos(wallpt(ik,8)),    
     >        wallpt(ik,2) + wallpt(ik,5) * sin(wallpt(ik,8)),    
     >        wallpt(ik,1) + wallpt(ik,6) * cos(wallpt(ik,9)),    
     >        wallpt(ik,2) + wallpt(ik,6) * sin(wallpt(ik,9)),
     >        ik    
        end do  
c
         write (6,*) 'NTrap:',wltrap1,wltrap2
         do ik = wltrap1,wltrap2
            write (6,'(4(1x,f18.10),i4)') 
     >        wallpt(ik,1) + wallpt(ik,5) * cos(wallpt(ik,8)),    
     >        wallpt(ik,2) + wallpt(ik,5) * sin(wallpt(ik,8)),    
     >        wallpt(ik,1) + wallpt(ik,6) * cos(wallpt(ik,9)),    
     >        wallpt(ik,2) + wallpt(ik,6) * sin(wallpt(ik,9)),
     >        ik    
         end do  
c
         write (6,*) 'NOuter:',wltrap2+1,wallpts
         do ik = wltrap2+1,wallpts
            write (6,'(4(1x,f18.10),i4)') 
     >        wallpt(ik,1) + wallpt(ik,5) * cos(wallpt(ik,8)),    
     >        wallpt(ik,2) + wallpt(ik,5) * sin(wallpt(ik,8)),    
     >        wallpt(ik,1) + wallpt(ik,6) * cos(wallpt(ik,9)),    
     >        wallpt(ik,2) + wallpt(ik,6) * sin(wallpt(ik,9)),
     >        ik    
         end do  
c
      endif 
c

C
C-----------------------------------------------------------------------
C
C     ASSIGN THE C-SPUTTERING DATA READ FROM PIN TO THE
C     WALL - LAUNCH PROBABILITY ARRAY. (IF REQUIRED)
C
C-----------------------------------------------------------------------
C
      IF ((WLPABS.EQ.2.OR.WLPABS.EQ.3).AND.CGEOOPT.NE.-1) THEN
C
C       THERE SHOULD BE A 1 TO 1 MAPPING OF THE WALL COORDINATES
C       AND CORRESPONDING POL VALUES IN THE WALLPOL ARRAY FOR THE
C       FIRST NWALL-1 SEGMENTS. THE NWALL POINT IS AN EXTRA
C       ARRAY ELEMENT CONTAINING THE END-POINT OF THE NWALL-1
C       WALL SEGMENT.
C
        NWLPROB2 = NWALL-1
        DO 2100 IK = 1, NWALL-1
           IF (WALLPOL(IK).LT.0) THEN
              WLPROB2(IK,1) = IK
              WLPROB2(IK,2) = IK
              WLPROB2(IK,3) = 0.0
           ELSE
              RES = FINDPOL(CWSPUT,SPUTCNT,WALLPOL(IK))
              IF (RES.EQ.0) THEN
                 WLPROB2(IK,1) = IK
                 WLPROB2(IK,2) = IK
                 WLPROB2(IK,3) = 0.0
              ELSE
                 WLPROB2(IK,1) = IK
                 WLPROB2(IK,2) = IK
                 WLPROB2(IK,3) = CWSPUT(RES,2)
              ENDIF
           ENDIF
 2100   CONTINUE
C
C       ADD C-SPUTTERING DATA FOR TRAP SEGMENT
C       AT THIS TIME, THIS ASSUMES THAT THE ONLY TRAP OPTION
C       USED IN COMBINATION WITH THE CX-SPUTTERING DATA
C       DATA FROM PIN IS THE ONE WHICH DRAWS THE TRAP
C       WALL FROM THE CORNER POINTS OF THE PLATES.
C
        NWLPROB2 = NWLPROB2 + 1
        WLPROB2(NWLPROB2,1) = NWLPROB2 + NPLAT -2
        WLPROB2(NWLPROB2,2) = NWLPROB2 + NPLAT -2
        WRITE(6,*) 'NWL2:',NWLPROB2,NPLAT,WLPROB2(NWLPROB2,1)
        RES = FINDPOL(CWSPUT,SPUTCNT,TRAPPOL)
        IF (RES.EQ.0) THEN
           WLPROB2(NWLPROB2,3) = 0.0
        ELSE
           WLPROB2(NWLPROB2,3) = CWSPUT(RES,2)
        ENDIF
C
        IF (NWLPROB.GT.0.AND.WLPABS.EQ.2) THEN
           DO 2200 IND = 1, NWLPROB
             DO 2200 IK = INT(WLPROB(IND,1)),INT(WLPROB(IND,2))
                WLPROB2(IK,3) = WLPROB2(IK,3) * WLPROB(IND,3)
 2200      CONTINUE
        ENDIF
C
        NWLPROB = NWLPROB2
        WRITE (6,*) 'NWL:',NWLPROB
        DO 2300 IK = 1,NWLPROB
           WLPROB(IK,1) = WLPROB2(IK,1)
           WLPROB(IK,2) = WLPROB2(IK,2)
           WLPROB(IK,3) = WLPROB2(IK,3)
           WRITE(6,*) IK,WLPROB(IK,1),WLPROB(IK,2),WLPROB(IK,3)
 2300   CONTINUE
C
      ENDIF
C
C-----------------------------------------------------------------------
C
c     Close the files that have been used for input
c
      close (inimout)
      close (ipinout)
c
 
 1000 FORMAT(A160)
 1010 FORMAT(A15)
 1020 FORMAT(A21)
 2000 FORMAT(I7,1P,9E13.5)
 2005 FORMAT(I7)
 3000 FORMAT(I4,2X,A4,1PE11.3,1X,6E11.3)
 4000 FORMAT(4X,2X,A4,1PE11.3,1X,6E11.3)
 9000 FORMAT(A9,2I6)
 9010 FORMAT(1P,6E12.4)
 9020 FORMAT(5X,I1,1X,I2,2I6)
 9030 FORMAT(12I6)
      RETURN
      END
C
C
      INTEGER FUNCTION FINDPOL(ARS,NARS,TNUM)
      IMPLICIT none
C     INCLUDE "PARAMS"
      include 'params'
      REAL ARS(MAXPTS,2)
      INTEGER NARS,TNUM
C
C-----------------------------------------------------------------------
C
C     THIS ROUTINE RETURNS THE INDEX INTO THE ARRAY ARS POINTING
C     TO THE LINE WHOSE CONTENTS RELATE TO THE QUANTITY TNUM.
C
C     THE FIRST ENTRY IN ARS IS THE NUMBER TNUM IN A REAL REPRESENTATION
C
C-----------------------------------------------------------------------
C
      INTEGER I
C
      FINDPOL = 0
C
      DO 100 I = 1,NARS
         IF (TNUM.EQ.INT(ARS(I,1))) THEN
            FINDPOL = I
            GOTO 200
         ENDIF
 100  CONTINUE
 200  CONTINUE
      RETURN
      END
c
c
c
      SUBROUTINE LOADGEO
      IMPLICIT NONE
C     INCLUDE "PARAMS"
      include 'params'
C     INCLUDE "COMTOR"
      include 'comtor'
C
C-----------------------------------------------------------------------
C
C     THIS SUBROUTINE CONTAINS THE PRE-ENTRED WALL,TARGET AND
C     TRAP SPECIFICATIONS FOR SPECIFIC GRID GEOMETRIES. THIS
C     SAVES BOTH TIME AND THE POSSIBILITY OF ERROR BY HAVING THE
C     DATA PRE-ENTERED IN A FORMAT THAT CAN BE CONVERTED FOR USE
C     DIRECTLY INSIDE DIVIMP.
C
C     THE DATA ARRAYS ARE TRANSFERRED TO THE EQUIVALENT ARRAYS THAT
C     WOULD BE USED FOR INPUT AND THEN PROCESSING PROCEEDS NORMALLY.
C     THIS SWITCH IS CONTROLLED BY THEN GEOOPT OPTION - THAT SELECTS
C     EITHER A PRE-DEFINED GEOMETRY OR ONE THAT IS ENTERED BY HAND.
C
C     THE DATA CONSISTS OF POL NUMBERS FOR EACH WALL-SEGMENT AS WELL
C     AS THE START CO-ORDINATE FOR EACH SEGMENT, THE END COORDINATE
C     IS THE START COORDINATE OF THE NEXT SEGMENT. THE POL NUMBERS ARE
C     USED TO ASSIGN THE C-SPUTTERING DATA FROM A PIN RUN TO THE
C     RESPECTIVE WALL SEGMENTS SO THAT A PROPERLY DISTRIBUTED WALL-LAUNC
C     CAN BE CALCULATED.
C
C-----------------------------------------------------------------------
C
C     LOCAL VARIABLES
C
      INTEGER I,J,K,IN
      INTEGER GEOOPT
      REAL    R1,Z1,R2,Z2,RC,ZC,RX,ZX
C
C     PARAMETER VALUES
C
      INTEGER G1WALLP,G1TARGP,G1TRAPP,G2WALLP,G2TARGP,G2TRAPP
      INTEGER G1TARGP2,G2TARGP2
      PARAMETER ( G1WALLP = 42, G1TARGP = 20, G1TRAPP = 1,
     >            G2WALLP = 44, G2TARGP = 18, G2TRAPP = 1)
      PARAMETER ( G1TARGP2 = G1TARGP/2, G2TARGP2= G2TARGP/2)
C
C     DATA VARIABLE DECLARATIONS
C
      INTEGER   MAXGEO
      PARAMETER (MAXGEO=2)
C
      INTEGER NWALLP(MAXGEO),NTARGP(MAXGEO,2),NTRAPP(MAXGEO)
      REAL    WALLDATA(MAXGEO,MAXPTS,3),TARGDATA(MAXGEO,MAXPTS,4)
      REAL    TRAPDATA(MAXGEO,MAXPTS,3)
C
C     DATA VALUES
C
      DATA NWALLP /G1WALLP,G2WALLP/
      DATA ((NTARGP(I,J),J=1,2),I=1,MAXGEO)
     >     /G1TARGP2,G1TARGP,
     >      G2TARGP2,G2TARGP/
      DATA NTRAPP /G1TRAPP,G2TRAPP/
C
C     DATA FOR GEOMETRY T56P302 SHOT 24719
C
C     WALLS
C
      DATA ((WALLDATA(1,J,K),K=1,3),J=1,G1WALLP+1)
     > /604, 2.806, 2.053,  603, 2.806, 2.051,  602, 2.807, 2.045,
     >  600, 2.823, 2.039,  598, 2.833, 2.036,  596, 2.843, 2.033,
     >  594, 2.910, 2.011,  592, 3.015, 1.963,  590, 3.140, 1.893,
     >  588, 3.324, 1.790,  586, 3.486, 1.666,  584, 3.663, 1.518,
     >  582, 3.818, 1.337,  580, 3.937, 1.113,  578, 4.040, 0.8398,
     >  576, 4.219, 0.5153, 574, 4.276, 0.1423, 572, 4.255,-0.3336,
     >  570, 4.068,-0.7856, 568, 3.872,-1.275,  566, 3.516,-1.641,
     >  564, 3.096,-1.918,  562, 2.624,-2.063,  560, 2.077,-1.803,
     >  558, 1.896,-1.344,  556, 1.830,-0.8434, 554, 1.752,-0.4143,
     >  552, 1.728,0.01131, 550, 1.747, 0.3694, 548, 1.795, 0.6890,
     >  546, 1.865, 0.9713, 544, 1.890, 1.216,  542, 1.912, 1.424,
     >  540, 1.959, 1.595,  538, 2.019, 1.723,  536, 2.078, 1.804,
     >  534, 2.111, 1.841,  532, 2.116, 1.846,  530, 2.119, 1.850,
     >  528, 2.129, 1.861,  526, 2.136, 1.869,  524, 2.138, 1.873,
     >  -2 , 2.137, 1.876/
C
C     TARGETS
C
      DATA ((TARGDATA(1,J,K),K=1,4),J=1,G1TARGP+2)
C
C     INNER
C
     > /614,2.137,1.876,14,  613,2.188,1.898,13,  612,2.232,1.920,12,
     >  611,2.279,1.946,11,  610,2.316,1.972,10,  609,2.342,1.998, 9,
     >  608,2.358,2.019, 8,  607,2.359,2.022,19,  606,2.361,2.025,18,
     >  605,2.363,2.028,17,  -2 ,2.366,2.034,-1,
C
C     OUTER
C
     >  624,2.542,2.064,17,  623,2.548,2.060,18,  622,2.551,2.058,19,
     >  621,2.554,2.056, 8,  620,2.557,2.054, 9,  619,2.584,2.043,10,
     >  618,2.622,2.037,11,  617,2.668,2.036,12,  616,2.715,2.040,13,
     >  615,2.759,2.046,14,  -2 ,2.806,2.053,-1/
C
C     TRAP
C
      DATA ((TRAPDATA(1,J,K),K=1,3),J=1,G1TRAPP+1)
     > /625, 2.366, 2.034,  -2 , 2.542, 2.064/
C
C     DATA FOR SH26308 GRID FILE GEOMETRY
C
C     OUTER WALLS
C
      DATA ((WALLDATA(2,J,K),K=1,3),J=1,G2WALLP+1)
     > / -1, 2.825, 1.973,  633, 2.840, 1.983,  631, 2.844, 1.981,
     >  629, 2.854, 1.979,  627, 2.873, 1.973,  625, 2.901, 1.963,
     >  623, 2.950, 1.944,  621, 3.015, 1.920,  619, 3.034, 1.913,
     >  617, 3.053, 1.905,  615, 3.344, 1.773,  613, 3.535, 1.631,
     >  611, 3.714, 1.468,  609, 3.842, 1.299,  607, 3.934, 1.092,
     >  605, 4.039, 0.8386, 603, 4.212, 0.5368, 601, 4.272, 0.1878,
     >  599, 4.260,-0.3011, 597, 4.079,-0.7705, 595, 3.866,-1.252,
     >  593, 3.519,-1.644,  591, 3.087,-1.870,  589, 2.610,-1.998,
     >  587, 2.071,-1.765,  585, 1.944,-1.243,  583, 1.821,-0.08044,
     >  581, 1.749,-0.03572,579, 1.732, 0.08305,577, 1.756, 0.4202,
     >  575, 1.803, 0.7196, 573, 1.869, 0.9838, 571, 1.943, 1.213,
     >  569, 1.946, 1.418,  567, 1.948, 1.587,  565, 1.987, 1.667,
     >  563, 1.991, 1.674,  561, 1.994, 1.681,  559, 2.017, 1.716,
     >  557, 2.042, 1.740,  555, 2.058, 1.754,  553, 2.067, 1.763,
     >  551, 2.072, 1.768,  -1 , 2.076, 1.771,  -2 , 2.152, 1.771/
C
C     TARGETS
C
      DATA ((TARGDATA(2,J,K),K=1,4),J=1,G2TARGP+2)
C
C     INNER
C
     >/642,2.152,1.771,16,  641,2.194,1.805,15,  640,2.228,1.838,14,
     > 639,2.252,1.868,13,  638,2.269,1.892,12,  637,2.278,1.906,11,
     > 636,2.279,1.909,21,  635,2.287,1.913,20,  634,2.284,1.919,19,
     > -2 ,2.289,1.930,-1,
C
C     OUTER
C
     > 651,2.644,2.026,19,  650,2.653,2.020,20,  649,2.659,2.016,21,
     > 648,2.662,2.015,11,  647,2.664,2.013,12,  646,2.678,2.006,13,
     > 645,2.703,1.996,14,  644,2.737,1.985,15,  643,2.778,1.977,16,
     > -2 ,2.825,1.973,-1/
C
C     TRAP
C
      DATA ((TRAPDATA(2,J,K),K=1,3),J=1,G2TRAPP+1)
     > /652, 2.289, 1.930,  -2 , 2.644, 2.026/
C
C     DEPENDING ON THE OPTIONS CHOSEN - ASSIGN THE LOADED DATA
C     TO THE EQUIVALENT POSITIONS IN THE INPUT ARRAYS - SO THAT THEY
C     CAN BE USED IN THE WALL AND PLATE ROUTINES DIRECTLY. ALSO, ASSIGN
C     THE POL VALUES TO A SECOND ARRAY SO THAT THE DATA READ FROM A PIN
C     RUN CAN BE APPROPRIATELY ASSIGNED BASED ON THE POL NUMBER.
C
      GEOOPT = CGEOOPT+1
C
C-----------------------------------------------------------------------
C
C     ASSIGN WALLS
C
C-----------------------------------------------------------------------
C
c
c
c
      if (cneur.eq.2.or.cneur.eq.3) then
c
c        REPLACED BY NEUTRAL WALL OPTION
c
c        Adjust the wall options to the corresponding values
c        as if the data had been read from the input file. The
c        only purpose of these additional wall-options is to
c        allow the specification of internally loaded data ...
c        otherwise they behave identically to their
c        counterparts.
c
c         wallswch = .true.
c         if (cionr.eq.4) then
c            cionr = 2
c         elseif (cionr.eq.5) then
c            cionr = 3
c         endif
c
         NWALL = NWALLP(GEOOPT) + 1
c
         WRITE(6,*)'GEO:',NWALLP(1),NWALLP(2) ,NWALL
         DO 100 I = 1,NWALL
            WALLCO(I,1) = WALLDATA(GEOOPT,I,2)
            WALLCO(I,2) = WALLDATA(GEOOPT,I,3)
            WALLPOL(I)  = INT(WALLDATA(GEOOPT,I,1))
 100     CONTINUE
         TRAPPOL = INT (TRAPDATA(GEOOPT,1,1))
      endif
C
C-----------------------------------------------------------------------
C
C     ASSIGN PLATES - INCLUDE RING NUMBER AND GENERATE AN ADDITIONAL
C     POINT FOR THE VIRTUAL RING. ALSO, NOTE THAT THE PLATES MUST BE
C     SPECIFIED BY BIN CENTRES - HOWEVER, THE DATA ENTERED HERE ARE
C     BIN EDGES AND SO SOME CONVERSION IS NECESSARY.
C
C-----------------------------------------------------------------------
C
C     INNER TARGETS
C
      if (ctargopt.eq.3) then
c
c        Ctargopt = 3 results in the pre-loaded values for the
c        plate coordinates being placed in the array in which
c        values normally read-in would be placed. The processing
c        then continues normally ... however, the ends of each
c        ring must be deleted.
c
         NPLAT = NTARGP(GEOOPT,1)+2
         DO 200 I = 1,NTARGP(GEOOPT,1)
C
C        ADD REGULAR POINTS
C
            R1 = TARGDATA(GEOOPT,I,2)
            Z1 = TARGDATA(GEOOPT,I,3)
            R2 = TARGDATA(GEOOPT,I+1,2)
            Z2 = TARGDATA(GEOOPT,I+1,3)
            RC = 0.5 * (R1+R2)
            ZC = 0.5 * (Z1+Z2)
            PLATCO(I+1,1) = TARGDATA(GEOOPT,I,4)
            PLATCO(I+1,2) = RC
            PLATCO(I+1,3) = ZC
C
C           ADD SPECIAL POINTS
C
            IF (I.EQ.1) THEN
C
C           ADD FIRST VIRTUAL POINT - INNER TARGET
C
               RX = 2.0 * R1 - RC
               ZX = 2.0 * Z1 - ZC
               PLATCO(I,1) = TARGDATA(GEOOPT,I,4) + 1.0
               PLATCO(I,2) = RX
               PLATCO(I,3) = ZX
            ELSEIF (I.EQ.NTARGP(GEOOPT,1)) THEN
C
C           ADD LAST VIRTUAL POINT - INNER TARGET
C
               RX = 2.0 * R2 - RC
               ZX = 2.0 * Z2 - ZC
               PLATCO(I+2,1) = TARGDATA(GEOOPT,I,4) - 1.0
               PLATCO(I+2,2) = RX
               PLATCO(I+2,3) = ZX
            ENDIF
 200     CONTINUE
C
C-----------------------------------------------------------------------
C
C        OUTER TARGETS
C
         DO 300 I = NTARGP(GEOOPT,2)+1,NTARGP(GEOOPT,1)+2,-1
C
C           ADD REGULAR POINTS
C
            IN =  NTARGP(GEOOPT,2) -I +1 +1
            R1 = TARGDATA(GEOOPT,I,2)
            Z1 = TARGDATA(GEOOPT,I,3)
            R2 = TARGDATA(GEOOPT,I+1,2)
            Z2 = TARGDATA(GEOOPT,I+1,3)
            RC = 0.5 * (R1+R2)
            ZC = 0.5 * (Z1+Z2)
            PLATCO(IN+1,4) = RC
            PLATCO(IN+1,5) = ZC
C
C           ADD SPECIAL POINTS
C
            IF (I.EQ.(NTARGP(GEOOPT,2)+1)) THEN
C
C           ADD FIRST VIRTUAL POINT - OUTER TARGET
C
               RX = 2.0 * R2 - RC
               ZX = 2.0 * Z2 - ZC
               PLATCO(IN,4) = RX
               PLATCO(IN,5) = ZX
            ELSEIF (I.EQ.(NTARGP(GEOOPT,1)+2)) THEN
C
C           ADD LAST VIRTUAL POINT - OUTER TARGET
C
               RX = 2.0 * R1 - RC
               ZX = 2.0 * Z1 - ZC
               PLATCO(IN+2,4) = RX
               PLATCO(IN+2,5) = ZX
            ENDIF
 300     CONTINUE
      endif
C
C-----------------------------------------------------------------------
C
      RETURN
      END
