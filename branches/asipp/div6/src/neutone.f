c     -*-Fortran-*-  
c
      SUBROUTINE LAUNCH_ONE (IMP,RST,ZST,RIZ,ZIZ,idstart,iwstart,
     >                   rc,temi,
     >                   cisti,sputy,
     >                   SSTRUK,MTCSTRUK,SMAIN,SEXIT,
     >                   SATIZ,SNEUT,SWALLN,MTCWALLN,SCENT,STMAX,
     >                   SEED,NRAND,
     >                   NEUTIM,SFAIL,neuttype,vout,vrec,RECLOSS)
      use mtc
      use error_handling
      IMPLICIT NONE
c      real rst,zst,riz,ziz,temi,cisti,sputy
      real rst,zst,riz,ziz,temi,sputy
      real*8  cisti
      integer rc,imp,idstart,iwstart
      real    SSTRUK,SMAIN,SEXIT,SATIZ,SNEUT,SWALLN,SCENT
      real    MTCSTRUK,MTCWALLN
      real    STMAX,sfail
      real*8  seed 
      real    neutim,vout,vrec,recloss
      integer nrand,neuttype
c
c     LAUNCH_ONE: This routine follows one recombined impurity atom 
c                 as a neutral from it's point of creation to 
c                 either re-ionization or other system removal.  
c
c                 1997, May 22
c
c     This routine uses either the 3-D free space launch Vel/Angle flag (8)
c     or the cosine distribution surface launch at a fixed energy (3) and
c     derives an estimate of the particle velocity from the input ion
c     temperature.   
c
c     NEUTTYPE defines the type of neutral being launched in case we
c     want to later gather additional statistics on these types of 
c     neutrals. Recombined neutrals are TYPE 6.
c
c     
c
c
      include    'params'
      include    'dynam1'
      include    'dynam3'
      include    'dynam4'
      include    'cgeom'
      include    'comtor'
      include    'cioniz'
      include    'cneut'
      include    'cneut2'
      include    'crand'
      include    'commv' 
c slmod begin
      include    'slcom'

      INTEGER in,min_in,count_relaunch
      REAL    dist,min_dist

      LOGICAL first_step,first_track,first_message,limit_message,
     .        result_zero,aggressive
      REAL    scale,last_vin
c slmod end
C
      CHARACTER FATE(6)*14
C
c     Variables related to wall reflection code
c
      REAL TDUM(MAXPTS),XDUM(MAXPTS),YDUM(MAXPTS),WORK(4*MAXPTS)
      REAL RESULT,ROLD,ZOLD
      INTEGER INDWORK(2,MAXPTS),NRFCNT,INDI,IND,TOTRF,MAXRF
      real refprob      
c
      integer maxnrfcnt
      parameter (maxnrfcnt=500)
c
      REAL RNEW,ZNEW,TNEW,Tnorm
      real statim
      integer iprod,ik,ir,mtccnt,mtcdat
c
      real*8 cistntot 
      integer id,ipos,it,ifate
      integer iklast,irlast
      logical sect
c
      integer intersect_result
c
      real    ranang,vin,yvelf,xvelf,ran,angle,temn
      real    ran1,ran2,anglan,tantru
      real    tmpvelf,xvelftmp,yvelftmp,vintmp
      real    rstmax,za02as,r,z,tstepn,k,s
      real*8  cist,mtccist
      external za02as,ipos
c
      integer :: iter_cnt

c
      logical griderr,reflect_neut
C
      DATA FATE /'REACHED WALL',           'REACHED CENTRE',
     >           'TIME = TMAX',            'STRUCK TARGET',
     >           'IONISED TO 1',           'FAILED LAUNCH'/
C
c     Total number of neutral timesteps that the particle has been in 
c     existence = cistntot
c
c     IPP/08 Krieger - initialized variable IT to avoid runtime error
c     in write statment further down
      it = 0
c slmod begin
c...  Initialization:
      aggressive = .TRUE.
      totrf = 0
      maxrf = 0
      first_step    = .TRUE.
      first_track   = .TRUE.
      IF (idstart.LE.0) first_track = .FALSE.  ! Only want the FIRST_TRACK code executing if the particle launched from a surface
      first_message = .TRUE.
      limit_message = .TRUE.
      scale = 1.0
      count_relaunch = 0
c slmod end
c
c     Set other initialization values - use iprod for local ion index 
c
      cistntot = cisti * qtim/fsrate  
      r = rst
      z = zst 
      iprod = imp
      mtccnt = 0
      mtccist = 0.0
      CIST  = 0.0
      iter_cnt = 0
c
c     USe neuttype to set up the routine for different behaviours ...
c
c     Value of MTCDAT is used as an index in the array that accumulates
c     MTC information - currently a maximum value of 3 is allowed.  
c
      if (neuttype.eq.7) then 
         mtcdat = 3
      else
         mtcdat = 2
      endif 
c
c     Random number limits - smaller than NEUT because only one at a time
c
c      KK    = 100 * ISECT
c      KKLIM = KK - 10
c
      statim =  ZA02AS (1) 
c
      DEBUGN = .FALSE.
c
      IF (CSTEPN.GT.0.0) THEN
        WRITE (6,9004) NINT(CSTEPN),FSRATE
        DEBUGN = .TRUE.
      ENDIF
c
C---- SET  RSTMAX: MAX NUMBER OF ITERATIONS UP TO TIME = TMAX (0.1S)
C---- (DIFFERS FROM CSTMAX WHICH APPLIES TO IONS,  BY FSRATE/QTIM AND
C---- BECAUSE TMAX IS 0.1 HERE BUT 10 FOR IONS).
C
      RSTMAX = 0.1 / FSRATE
      sneut = sneut + sputy
C
C---- SET START 
C

      TSTEPN= CSTEPN

      IF (DEBUGN) WRITE (6,9005)
c
      IK    = 0
      IR    = 0
C
      call gridpos(ik,ir,r,z,.true.,griderr)
C
      IF (griderr) THEN
c
c        Need to consider deposition of this particle
c
c
c        If the ion is not on the grid -exit with an appropriate error
c
         sfail = sfail + sputy
c
         call errmsg('NEUTONE:LAUNCH_ONE:',
     >               'INITIAL NEUTRAL PARTICLE POSITION IS OFF GRID')
c         WRITE(6,*) 'LAUNCH_ONE ERROR: R,Z,IK,IR<NEUTTYPE:',
c     >            R,Z,IK,IR,neuttype
c
         IFATE = 6
c
         goto 899
C
      ENDIF
c
c     SET status flags for entering and exiting core plasma 
c      
      CFLRIN= .TRUE.
      CFLREX= .TRUE.
c
c     Record K and S values for debug print-out and possible summaries
c
      K     = KKS(IR)
      S     = KSS(IK,IR)

c
C----------------------------------------------------------------------
c     NEUTTYPE = 6  (or other) - use VEL/ANG flag 8 - modified
c     NEUTTYPE = 7               use VEL/ANG flag 3   
C----------------------------------------------------------------------
c
c


c     
C     V/A flag 8 - isotropic - could be going any direction in 2-D 
c     
c     The proper modifier for the 3-D velocity is used later - this simply
c     decides which direction the particle will move. 
c     
      nrand = nrand + 1
      CALL SURAND2 (SEED, 1, RAN)
c     
c     Select angles for launch
c
      if (neuttype.eq.7) then 
         nrand = nrand + 1
         CALL SURAND2 (SEED, 1, RAN1)
c
         ANGLAN = SIGN (ASIN (SQRT(RAN)), RAN1-0.5)
         TANTRU = THETAS(IDSTART)
c
         ANGLE = ANGLAN + TANTRU
c
      else

         ANGLE = 2.0 * PI * RAN - PI

      endif
c
c      VMULT = 1.0
C     
C----------------------------------------------------------------------
C     CALCULATE LAUNCH VELOCITY - based on Ion temperature
C----------------------------------------------------------------------
C     

      if (cneutvel.eq.0) then 
         VIN = 1.38E4 * SQRT (TEMI/CRMI) * cvrmult
      elseif (cneutvel.eq.1.or.cneutvel.eq.2) then 
         VIN = 1.38E4 * SQRT (ktibs(ik,ir)/CRMI) * cvrmult
      endif 

c      write(6,*) 'debug: vin 1 ',vin,cneutvel
c     
c     Assign neutral temperature
c     
      TEMN   = CRMI * (VIN/1.38E4) * (VIN/1.38E4)
c     
c     Choose angle for 3rd dimension
c     
      if (neuttype.ne.7) then 
c
         nrand = nrand + 1
         CALL SURAND2 (SEED, 1, RAN)
         ranang = PI * ran - PI/2.0
c     
c        Adjust velocity for projection to 2-D poloidal plane
c
c slmod begin
c...
         do while (cos(ranang).lt.0.1) 
            nrand = nrand + 1
            CALL SURAND2 (SEED, 1, RAN)
            ranang = PI * ran - PI/2.0
            if (limit_message) then
               CALL WN('LAUNCH_ONE','Limiting 3D velocity correction')
               limit_message = .FALSE.
            endif
         enddo
c slmod end

         vin = vin * cos(ranang)
c
      endif
c
c      write(6,*) 'debug: vin 3 ',vin,neuttype

      vrec = vin
C     
C-----CALCULATE X,Y COMPONENTS OF VELOCITY, PROB OF IONISATION, ETC
C     
      XVELF  = VIN * COS(ANGLE) * FSRATE
      YVELF  = VIN * SIN(ANGLE) * FSRATE
c     
      IF (DEBUGN) THEN
         WRITE (6,9003) IPROD,CIST,IK,IR,0,0,R,Z,K,
     >     VIN,TEMN,SPUTY,ANGLE*RADDEG,IT,'RECOMBINED LAUNCH'
      ENDIF
c      WRITE (6,9003) IPROD,CIST,IK,IR,idstart,0,R,Z,K,
c     >  VIN,TEMN,SPUTY,ANGLE*RADDEG,IT,'RECOMBINED LAUNCH debug:'
C
C     SET UP THE WALL DEFINITION IN THE ROUTINE THAT
C     DETERMINES WHETHER A POINT IS INSIDE OR OUTSIDE THE
C     BOUNDARY.
C     
      NRFCNT = 0
C     
      CALL GA15A(PCNT,1,WORK,4*MAXPTS,INDWORK,MAXPTS,
     >           RW,ZW,TDUM,XDUM,YDUM,6)
C     
C     
C---------------------------------------------------------------------
C  ITERATE AROUND MOTION FOR EACH TIMESTEP UNTIL IONISATION OCCURS.
C  GENERATE NEW SET OF RANDOMS WHEN OLD LOT ARE ALL USED UP
C-----------------------------------------------------------------------
C
        iklast = ik
        irlast = ir
c
  200   CONTINUE
c
c       Increment iteration count
c        
        iter_cnt = iter_cnt + 1
c
c       Cutoff time
c
        IF (CIST.GT.RSTMAX) THEN
          STMAX = STMAX + SPUTY
          IFATE = 3
          GOTO 899
        ENDIF
c slmod begin
        last_vin = vin
c slmod end
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
c           Assign neutral temperature
c     
            TEMN   = CRMI * (VIN/1.38E4) * (VIN/1.38E4)
c     
c           Reset velocity components
c
            xvelf = xvelftmp * vin/vintmp              
            yvelf = yvelftmp * vin/vintmp              
c
        endif
c
c       Record current cell
c
        iklast = ik
        irlast = ir 
C
C------ UPDATE (X,Y) COORDINATES AND ARRAY INDICES
C
        ROLD = R
        ZOLD = Z
c slmod begin
        R    = R + XVELF * SCALE
        Z    = Z + YVELF * SCALE
c        R    = R + XVELF
c        Z    = Z + YVELF
c slmod end
c
        call gridpos(ik,ir,r,z,.false.,griderr)
C
        CALL GA15B(R,Z,RESULT,PCNT,1,WORK,4*MAXPTS,
     >             INDWORK,MAXPTS,RW,ZW,TDUM,XDUM,YDUM,6)

c slmod begin
c        write(6,*) '  debug: result',result,first_track

        IF (first_step.AND.result.LT.0.0) THEN
c...      The neutral particle has taken one step and is already outside the computational 
c         domain.  Try again, but with a reduced velocity and maybe (if this is a target
c         launch) also move the launch point closer to the centre of the target.  This is 
c         not a perfect solution since the particle accounting assumes a particular
c         time-step, which is no longer the case with this "fix", but so far it is a robust
c         way of preserving the neutral:

          IF (first_message) THEN
            CALL WN('LAUNCH_ONE','Neutral re-launch required')
            first_message = .FALSE.
          ENDIF

c          write(6,*) '  first_track,result:',first_track,result,scale

          IF (DEBUGN) 
     >       WRITE (6,9003) IPROD,CIST,IK,IR,0,0,R,Z,K,
     >        VIN,TEMN,SPUTY,ANGLE*RADDEG,IT,'RE-LAUNCHING'
c          WRITE (6,9003) IPROD,CIST,IK,IR,idstart,0,R,Z,K,
c     >      VIN,TEMN,SPUTY,ANGLE*RADDEG,IT,'RE-LAUNCHING debug:'

          ik = iklast  ! perhaps not required
          ir = irlast

          IF (first_track) THEN 
c           Move the launch point closer to the center of the target segment, and also
c           reduce the velocity (see the application of SCALE above):
            scale = scale * 0.9
            r = rst + (1.0 - scale) * (rp(idstart) - rst)
            z = zst + (1.0 - scale) * (zp(idstart) - zst)
          ELSE
c           Reduce the velocity, and therefore the distance traveled over the first
c           step, and reset the position to R,ZOLD from the same point:       
            scale = scale * 0.1
            r = rold
            z = zold
          ENDIF
            
c         Record an error and lose the particle of the above remedies are not working:
          IF ((     (first_track).AND.scale.LT.0.1    ).OR.
     >        (.NOT.(first_track).AND.scale.LT.0.00001)) THEN 
              write(0,'(a)') 'NEUTONE: ERROR: UNABLE TO '//
     >                       ' GENERATE A VALID LAUNCH'
              write(6,'(a)') 'NEUTONE: ERROR: UNABLE TO '//
     >                       ' GENERATE A VALID LAUNCH'
              write(6,'(a,3i6,4g12.5)') 
     >          'DATA:',ik,ir,indi,rnew,znew,xvelf,yvelf    
              write(6,'(a,5i6)') 
     >          'DATA:',wlwall1,wlwall2,wltrap1,wltrap2,wallpts
              write(6,*) 
     >          'DATA:',(id,':',wallpt(indi,id),',',id=1,13)

              if (aggressive) then

c               Find segment closest to current position and re-launch from   *** PUT IN A SEPARATE ROUTINE ***
c               the centre of the segment:                

                count_relaunch = count_relaunch + 1
                IF (count_relaunch.EQ.50) THEN
                  write(0,*) '  problem again 3'
                  write(6,*) '  problem again 3'
c                  stop
                ENDIF 

                min_dist = HI
                DO in = 1, wallpts
                  dist = (wallpt(in,1) - r)**2 + (wallpt(in,2) - z)**2
                  IF (dist.LT.min_dist) THEN
                    min_dist = dist
                    min_in   = in
                  ENDIF
                ENDDO               
                r     = wallpt(min_in,1)
                z     = wallpt(min_in,2)
                angle = wallpt(min_in,9) - 1.57  ! subtracting 90 degrees from the CW launch angle is required 
                vin   = last_vin
                xvelf = vin * COS(angle) * fsrate
                yvelf = vin * SIN(angle) * fsrate
                first_step = .TRUE.
                scale = 1.0
                CALL WN('LAUNCH_ONE','Attempting forced re-launch 3')
                WRITE(0,*) '  MIN_IN = ',min_in,vin,angle,
     .                     r,z,xvelf,yvelf
                WRITE(6,*) '  debug: MIN_IN = ',min_in,vin,angle,
     .                     r,z,xvelf,yvelf
                WRITE (6,9003) IPROD,CIST,IK,IR,idstart,0,R,Z,K,
     >            VIN,TEMN,SPUTY,ANGLE*RADDEG,IT,
     >            'FORCED RE-LAUNCH 3 debug:'



                GOTO 200

              endif

              SSTRUK = SSTRUK + SPUTY
	      ! Add to deposition.
              IF (WALLPT(INDI,18).EQ.0) THEN
                WRITE(0,*) 'ERROR: WALLPT(INDI,18).EQ.0 in '//
     >                     'NEUTONE.f when dropping trapped neuts'
                WRITE(6,*) 'ERROR: WALLPT(INDI,18).EQ.0 in '//
     >                     'NEUTONE.f when dropping trapped neuts'
              ELSE
                NEROS (INT (wallpt (INDI,18)),1) = 
     >                     NEROS (INT (wallpt (INDI,18)),1) + SPUTY
              ENDIF
              if (mtccnt.gt.0) then 
                 MTCSTRUK = MTCSTRUK + SPUTY
              endif
              IFATE = 4
c              STOP 'oops'
              GOTO 899            
           ENDIF

           GOTO 200  ! ...want some sauce with that spagetti?
        ENDIF
c...    Clear the flags and modifier once the particle is safely back inside the 
c       computational domain:
        first_step  = .FALSE.
        first_track = .FALSE.
        scale       = 1.0
c slmod end
c
c       If particle is within the defined outer wall
c
c
        IF (RESULT.GE.0.0) THEN
c slmod begin
          result_zero = .FALSE.
          IF (RESULT.EQ.0.0) result_zero = .TRUE.
c slmod end
          CIST = CIST + 1.0
          K    = KKS(IR)
          S    = KSS(IK,IR)
c
c         IF particle is NOT on the defined grid - but still inside wall.          
c
          if (griderr) then  
c
c            Particle not in grid - record void region density data.
c            Estimate which region the neutral lies in - 
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
                      ddvoid(2) = ddvoid(2) + sputy   
c
c                  Other Divertor void regions. 
c
                   else
c
                      ddvoid(3) = ddvoid(3) + sputy   
c
                   endif 
c
c               Main plasma void region
c
                elseif (z.lt.zxp) then  

                   ddvoid(1) = ddvoid(1) + sputy   

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
                      ddvoid(2) = ddvoid(2) + sputy   
c
c                  Other Divertor void regions. 
c
                   else
c
                      ddvoid(3) = ddvoid(3) + sputy   
c
                   endif 
c
c               Main plasma void region
c
                elseif (z.gt.zxp) then  

                   ddvoid(1) = ddvoid(1) + sputy   

                endif
c              
             endif

c
c         Particle is somewhere on the polygonal grid 
c
          elseIF (.not.griderr) THEN
C
C------   CHECK IF ENTERED MAIN PLASMA, REACHED CENTRE, OR SURVIVED TO
C------   TIME CUTOFF POINT ... OR HIT PLATES
C
c           Entered or exited MAIN
c
            IF     (CFLRIN.AND.IR.LT.IRSEP) THEN
              CFLRIN = .FALSE.
              SMAIN = SMAIN + SPUTY
              ELIMS(IK,3,0) = ELIMS(IK,3,0) + SPUTY
            ELSEIF (CFLREX.AND.IR.GE.IRSEP) THEN
              IF (.NOT.CFLRIN) THEN
                CFLREX = .FALSE.
                SEXIT = SEXIT + SPUTY
                ELIMS(IK,1,0) = ELIMS(IK,1,0) + SPUTY
              ENDIF
            ENDIF
C
c           Central Mirror
c
            IF (IR.EQ.1) THEN
              SCENT = SCENT + SPUTY
              IFATE = 2
              GOTO 899
            endif
c
c           Cutoff time
c
c            ELSEIF (CIST.GT.RSTMAX) THEN
c              STMAX = STMAX + SPUTY
c              IFATE = 3
c              GOTO 899
c            ENDIF
C
C------     SCORE PARTICLE IN DDLIMS ARRAY IN "IONISATION 0" POSITION.
C------     IF TIME POINT REACHED SCORE TIME POSITION ALSO, AND INCREMENT.
C
            DDLIMS(IK,IR,0) = DDLIMS(IK,IR,0) + SPUTY
            DDTS(IK,IR,0) = DDTS(IK,IR,0) + SPUTY * TEMN
c
c           Record time resolved data if any requested
c     
c           USE total time so far for time-resolved data
c     
            if (nts.gt.0) then
c               IT    = IPOS (SNGL(CISTNTOT+DBLE(CIST)),
               IT    = IPOS (SNGL(CISTNTOT+CIST),
     >                            CTIMES(1,0), NTS+1)
               LIMS(IK,IR,0,IT) = LIMS(IK,IR,0,IT) + SPUTY
            ENDIF
C
C------     SET NEW IONISATION PROBABILITY
C------     CHECK FOR IONISATION: IF IT HAS NOT OCCURED JUMP BACK FOR
c------     ANOTHER ITERATION.
C
            IF (DEBUGN) THEN
              IF (CIST.GE.TSTEPN) THEN
  495           TSTEPN = TSTEPN + CSTEPN
                IF (TSTEPN.LE.CIST) GOTO 495
                WRITE (6,9003) IPROD,CIST,IK,IR,0,0,R,Z,K,
     >            VIN,TEMN,SPUTY,ANGLE*RADDEG,IT
              ENDIF
            ENDIF
C
c          Check for ionization  
c
           nrand = nrand + 1
           CALL SURAND2 (SEED, 1, RAN)
c
           IF (RAN.GE.KPCHS(IK,IR,0)) GOTO 200
c
c  State change event has occurred -
C  ------------------------------------------------------------------
C
C           1) ionization
c           2) charge exchange collision
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
            if (ran.gt.kpizs(ik,ir)) then 

               call execute_mtc(mtcdat,mtccnt,mtccist,cist,mtcinf,
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
c               mtcinf(1,mtcdat) = mtcinf(1,mtcdat) + sputy
c
c               if (mtccnt.eq.1) then 
c                  mtcinf(2,mtcdat) = mtcinf(2,mtcdat) 
c     >                             + sputy * (cist-mtccist)
c                  mtcinf(4,mtcdat) = mtcinf(4,mtcdat) + sputy * 
c     >                               (cist-mtccist)*abs(vin)*fsrate
c               endif 
c
c               mtcinf(3,mtcdat) = mtcinf(3,mtcdat) 
c     >                          + sputy * (cist-mtccist)
c               mtcinf(5,mtcdat) = mtcinf(5,mtcdat) + sputy *
c     >                                (cist-mtccist)*abs(vin)*fsrate
c               mtcinf(6,mtcdat) = mtcinf(6,mtcdat) + sputy * temn 
c               mtcinf(7,mtcdat) = mtcinf(7,mtcdat) + sputy * abs(vin)
c
c              Update MTC time 
c
c               mtccist = cist
c
c              Calculate Event  
c
c              Choose randomly between +90 and -90 change in velocity
c              - even probability.
c
c               nrand = nrand + 1
c               CALL SURAND2 (SEED, 1, RAN)
c
c               if (ran.ge.0.5) then 
c
c                 counter-clockwise 90 degrees
c
c                  tmpvelf = xvelf
c                  xvelf   = -yvelf
c                  yvelf   = tmpvelf
c
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
c                 Assign neutral temperature
c     
c                  TEMN   = CRMI * (VIN/1.38E4) * (VIN/1.38E4)
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
           TIZS(IK,IR,0) = TIZS(IK,IR,0) + SPUTY
c
           RIZ = R
           ZIZ = Z
c
           NRAND = NRAND + 1
           CALL SURAND2 (SEED, 1, RAN)
c
c          This output velocity needs to be rescaled by the
c          ion timestep for use in the calling routines - in the case
c          of reionization of the neutral particle. The approximation
c          is 1/2 of the total particle velocity along the field line.  
c
c          This appropriate scaling is performed in the calling routines
c           
           vout =  SIGN (0.5*VIN, RAN-0.5) 
c
           sATIZ = sATIZ + SPUTY
           IFATE = 5
           GOTO 899
         ENDIF
C
C       Particle has stepped outside the wall boundary 
C
        ELSEIF (RESULT.LT.0.0) THEN
c slmod begin
          first_step = .TRUE.

c...      This check handles the (very) rare case that RESULT was exactly equal to 
c         0.0 on the previous cycle, in which case the wall intersection code
c         in find_wall_intersection (called below) will fail.  To account for this,
c         slide R,ZOLD back a little (0.1% of step length), to improve the chances 
c         that the intersection is well defined: 
          IF (result_zero) THEN
            CALL WN('LAUNCH_ONE','R,ZOLD slide in effect')
c            write(6,*) '  debug: sliding...'
            rold = r + 1.001 * (rold - r)
            zold = z + 1.001 * (zold - z)
            result_zero = .FALSE.
          ENDIF
c slmod end
c
c       WALL COLLISION AND BOUNDARY CODE 
c
C       PARTICLE IS FOUND TO BE OUTSIDE BOUNDARY. NEED TO FIND
C       THE POSITION AND SEGMENT WHERE IT LEFT...
C
c         WRITE(6,'(a,2i4,1p,3g12.5,l4)') 'WALL COLL:',IK,IR,
c     >                          R,Z,result,griderr
c         write(6,'(a,i6,1p,7g12.5)') '         :',iprod,rold,zold,
c     >          rst,zst,xvelf,yvelf,cist
C
c
c         Comment out old wall intersection calculation code
c
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
C         WRITE(6,*)  'DSQ:',DSQ,IND,WALLPTS
C
c          CALL INTCALC(R,Z,ROLD,ZOLD,WALLPT(IND,1),WALLPT(IND,2),
c     >                 WALLPT(IND,8),WALLPT(IND,9),WALLPT(IND,5),
c     >                 WALLPT(IND,6),RNEW,ZNEW,TNEW,tnorm,
c     >                 SECT,nrfopt)
c          INDI = IND
C
C          WRITE(6,*) 'SECT:',SECT
C
c          IF (.NOT.SECT) THEN
c            DO 660 ID = 1,WALLPTS/2
c              INDI = IND + ID
C
C              WRITE(6,*) 'INDI1:' ,INDI,SECT
C
c              IF (INDI.GT.WALLPTS) INDI = INDI-WALLPTS
c              CALL INTCALC(R,Z,ROLD,ZOLD,WALLPT(INDI,1),WALLPT(INDI,2),
c     >                 WALLPT(INDI,8),WALLPT(INDI,9),WALLPT(INDI,5),
c     >                 WALLPT(INDI,6),RNEW,ZNEW,TNEW,tnorm,
c     >                 SECT,nrfopt)
c              IF (SECT) GOTO 670
c              INDI = IND - ID
C
C              WRITE(6,*) 'INDI2:' ,INDI,SECT
C
c jdemod - note - this is a bug - since INDI <= 0 - the quantity must be 
c                 ADDED to wallpts to get the proper index - not subtracted!
c              IF (INDI.LT.1) INDI = WALLPTS-INDI
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
c
c         Find location of wall intersection - use of indi for intersection index is 
c         for compatibility with existing code 
c
c slmod begin          
          if (cprint.ge.1) 
     >      write(6,'(a,i10,2i6,5g18.10)')  
     >            'NEUTONE:FIND_WALL_INTERSESCTION-BEFORE:',
     >            iter_cnt,ik,ir,r,z,rold,zold,result
c
c          write(6,'(a,i10,2i6,5g18.10)')  
c     >          'NEUTONE:FIND_WALL_INTERSESCTION-BEFORE:',
c     >          iter_cnt,ik,ir,r,z,rold,zold,result
c slmod end

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
c slmod begin          
          if (cprint.ge.1) 
     >      write(6,'(a,i10,2i6,9g18.10,i10,l6)')  
     >            'NEUTONE:FIND_WALL_INTERSESCTION-AFTER:',
     >            iter_cnt,ik,ir,r,z,rold,zold,result,rnew,znew,
     >            tnew*raddeg,tnorm*raddeg,intersect_result,sect
c
c          write(6,'(a,i10,2i6,9g18.10,i10,l6)')  
c     >          'NEUTONE:FIND_WALL_INTERSESCTION-AFTER:',
c     >          iter_cnt,ik,ir,r,z,rold,zold,result,rnew,znew,
c     >          tnew*raddeg,tnorm*raddeg,intersect_result,sect
c slmod end
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
c             write(6,*) '  debug: index check',indi,wallpt(indi,25)
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
c                 write(6,*) '  debug: reflect_neut',ran,refprob
                 if (ran.le.refprob) reflect_neut = .true.
c                 write(6,*) '  debug: reflect_neut',reflect_neut
c
             endif                
c
          endif
c
c         Neutral is not reflected
c
          IF (.not.reflect_neut) THEN
            IF (.NOT.SECT) THEN
              WRITE(6,'(a,i10,2g12.5)') 'NEUTONE: ERROR: NEUTRAL'//
     >                              ' AT WALL BUT NO INTERSECTION'//
     >                              ' POINT FOUND',iter_cnt,R,Z
              WRITE(0,'(a,i10,2g12.5)') 'NEUTONE: ERROR: NEUTRAL'//
     >                              ' AT WALL BUT NO INTERSECTION'//
     >                              ' POINT FOUND',iter_cnt,R,Z
c slmod begin
c...          I don't quite understand this, since there's an error message
c             but the particle recording appears to be standard, i.e. a 
c             collision is registered eventhough no intersection is found:
c              IF (aggressive) THEN
c                stop 'aggressive 2'
c              ENDIF
c slmod end
            ENDIF
c
c           If the neutral has struck the target
c
c slmod begin
c
c
c            jdemod - modify target detection to use wallpt(indi,18)
c                     for all grids
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
c
               SSTRUK = SSTRUK + SPUTY
c

               ! Add to deposition.
               NEROS (INT (wallpt (INDI,18)),1) = 
     >              NEROS (INT (wallpt (INDI,18)),1) + SPUTY

               if (mtccnt.gt.0) then 
                  MTCSTRUK = MTCSTRUK + SPUTY
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
c
               IFATE = 4
               WRITE(6,*) ' COLLISION WITH TARGET',INDI,R,Z,IK,IR,
     >                      wallpt(indi,1),wallpt(indi,2)
               GOTO 899
c
c           Neutral struck a wall segment - that is not a target element
c
            ELSE
c
               SWALLN = SWALLN + SPUTY
c
               if (mtccnt.gt.0) then 
                  MTCWALLN = MTCWALLN + SPUTY
               endif
c
C
c              Assign a value to ir that corresponds to the
c              index of the wall segment crossed.
c
c slmod begin
c
c              jdemod - change to simple algorithm for all grids - if 
c                       impact not on trap ring it must be on wall
c                       However, there are now grid cases (ribbon) where
c                       wltrap1, wltrap2 are not well defined
c
c
c               if (grdnmod.ne.0) then
c
c
                 if (indi.ge.wltrap1.and.indi.le.wltrap2) then 
c...                WLTRAP1,2 are still well defined for the generalized geometry, 
c                   and probably always will be (same login in NEUT.F):
                    ir = irtrap
                 else
c...                Everything else must be IRWALL (or not..?):
                    ir = irwall
                 endif

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
               WRITE(6,*) ' COLLISION WITH WALL',INDI,R,Z,IK,IR
               GOTO 899
            ENDIF
c
c         Deal with neutral reflection 
c
          ELSEIF (reflect_neut) THEN
c
c           No intersection point found - record error and stop following 
c           particle
c
c            write(6,*) '  debug: reflect_neut happening'
            IF (.NOT.SECT) THEN
c
              WRITE(6,'(a,i10,2g12.5)') 'NEUTONE: ERROR: NEUTRAL'//
     >                              ' AT WALL BUT NO INTERSECTION'//
     >                              ' POINT FOUND',iter_cnt,R,Z
              WRITE(0,'(a,i10,2g12.5)') 'NEUTONE: ERROR: NEUTRAL'//
     >                              ' AT WALL BUT NO INTERSECTION'//
     >                              ' POINT FOUND',iter_cnt,R,Z
c slmod begin
              IF (aggressive) THEN
c               Find segment closest to current position and re-launch from   *** PUT IN A SEPARATE ROUTINE ***
c               the centre of the segment:                

                count_relaunch = count_relaunch + 1
                IF (count_relaunch.EQ.50) THEN
                  write(0,*) '  problem again'
                  write(6,*) '  problem again'
c                  stop
                ENDIF 

                min_dist = HI
                DO in = 1, wallpts
                  dist = (wallpt(in,1) - r)**2 + (wallpt(in,2) - z)**2
                  IF (dist.LT.min_dist) THEN
                    min_dist = dist
                    min_in   = in
                  ENDIF
                ENDDO               
                r     = wallpt(min_in,1)
                z     = wallpt(min_in,2)
                angle = wallpt(min_in,9) - 1.57  ! subtracting 90 degrees from the CW launch angle is required 
                vin   = last_vin
                xvelf = vin * COS(angle) * fsrate
                yvelf = vin * SIN(angle) * fsrate
                first_step = .TRUE.
                CALL WN('LAUNCH_ONE','Attempting forced re-launch')
                WRITE(0,*) '  MIN_IN = ',min_in,vin,angle,
     .                     r,z,xvelf,yvelf
                WRITE(6,*) '  debug: MIN_IN = ',min_in,vin,angle,
     .                     r,z,xvelf,yvelf
                WRITE (6,9003) IPROD,CIST,IK,IR,idstart,0,R,Z,K,
     >            VIN,TEMN,SPUTY,ANGLE*RADDEG,IT,
     >            'FORCED RE-LAUNCH debug:'
                GOTO 200
              ENDIF
c slmod end
c
c             jdemod - change to wallpt(indi,18)
c 
              if (wallpt(indi,18).ne.0.0) then 
c
c              IF ((INDI.GT.WLWALL2.AND.INDI.LT.WLTRAP1).OR.
c     >            (INDI.GT.WLTRAP2.AND.INDI.LE.WALLPTS)) THEN
c
                 SSTRUK = SSTRUK + SPUTY
c
  	         ! Add to deposition.
                  NEROS (INT (wallpt (INDI,18)),1) = 
     >               NEROS (INT (wallpt (INDI,18)),1) + SPUTY
c
                 if (mtccnt.gt.0) then 
                    MTCSTRUK = MTCSTRUK + SPUTY
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
                 IFATE = 4
                 GOTO 899
              ELSE
c
                 sWALLN = sWALLN + SPUTY
c
                 if (mtccnt.gt.0) then 
                    MTCWALLN = MTCWALLN + SPUTY
                 endif
c
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
                 GOTO 899
              ENDIF
c
c           Intersection point found
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
c     >                   ' WITH ZERO REFLECTION COEFFICIENT',indi
c
c              elseif (wallpt(indi,25).gt.0.0) then 
c
c                Modify particle weight but leave all else unchanged
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
c                 Assign neutral temperature
c     
                  TEMN   = CRMI * (VIN/1.38E4) * (VIN/1.38E4)
c
              endif             
c
              R = RNEW
              Z = ZNEW
c
c             This is handled after line 200 on the loop back IF we don't move the particle by 
c             one step. The movement of the particle by one step was an attempt to ensure that 
c             the particle was inside the boundary - however, there are geometric locations where
c             a single step will move the particle outside again. In order to deal with this 
c             problem the position_on_target code and find_wall_intersection have been revised
c             to return points garanteed to be inside the wall. 
c
c             Update Rold and Zold before moving the particle one step
c
c              rold = r
c              zold = z
c
c             Recalculate the particle trajectory  
c              
c slmod begin
c              WRITE(6,'(A,1P,8E14.6,0P)') 
c     .          'NRF 1:',RNEW,ZNEW,tnew,R,Z,XVELF,YVELF,fsrate
c slmod end
              XVELF  = VIN * COS(TNEW) * FSRATE
              YVELF  = VIN * SIN(TNEW) * FSRATE
c
c             Do NOT Update particle position by one step - the loop back to the top will take care of this
c
c              R = R + XVELF
c              Z = Z + YVELF
c
c              CIST = CIST + 1.0
c
              NRFCNT = NRFCNT + 1
              TOTRF  = TOTRF +1
              MAXRF = MAX0(NRFCNT,MAXRF)
C
c slmod begin
c              WRITE(6,'(A,1P,7E14.6,0P)') 
c     .          'NRF 2:',RNEW,ZNEW,tnew,R,Z,XVELF,YVELF
c              WRITE(6,*) 
c              WRITE(6,*) 'NRF:',RNEW,ZNEW,tnew,R,Z,XVELF,YVELF
c slmod end
C
              IF (NRFCNT.GT.maxnrfcnt) THEN
c slmod begin
                IF (aggressive) THEN
c                 Find segment closest to current position and re-launch from
c                 the centre of the segment:                
	       
                  count_relaunch = count_relaunch + 1
                  IF (count_relaunch.EQ.50) THEN
                    write(0,*) '  problem again 2'
                    write(6,*) '  problem again 2'
c                    stop
                  ENDIF 
	       
                  min_dist = HI
                  DO in = 1, wallpts
                    dist = (wallpt(in,1) - r)**2 + (wallpt(in,2) - z)**2
                    IF (dist.LT.min_dist) THEN
                      min_dist = dist
                      min_in   = in
                    ENDIF
                  ENDDO               
                  r     = wallpt(min_in,1)
                  z     = wallpt(min_in,2)
                  angle = wallpt(min_in,9) - 1.57  ! subtracting 90 degrees from the CW launch angle is required 
                  vin   = last_vin
                  xvelf = vin * COS(angle) * fsrate
                  yvelf = vin * SIN(angle) * fsrate
                  first_step = .TRUE.
                  CALL WN('LAUNCH_ONE','Attempting forced re-launch 2')
c                  WRITE(0,*) '  MIN_IN = ',min_in,vin,angle,
c     .                       r,z,xvelf,yvelf
                  WRITE(6,*) '  debug: MIN_IN = ',min_in,vin,angle,
     .                       r,z,xvelf,yvelf
                  WRITE (6,9003) IPROD,CIST,IK,IR,idstart,0,R,Z,K,
     >              VIN,TEMN,SPUTY,ANGLE*RADDEG,IT,
     >              'FORCED RE-LAUNCH 2 debug:'
                  GOTO 200
                ENDIF
c slmod end
c
                recloss = recloss + sputy
c
c             jdemod - change to wallpt(indi,18)
c
                if (wallpt(indi,18).ne.0.0) then 
c
c                IF ((INDI.GT.WLWALL2.AND.INDI.LT.WLTRAP1).OR.
c     >              (INDI.GT.WLTRAP2.AND.INDI.LE.WALLPTS)) THEN
c
                   write(0,'(a)') 'NEUTONE: ERROR: NEUTRAL AT WALL'//
     >                        ' REMOVED DUE TO TOO MANY REFLECTIONS'//
     >                        ': PARTICLE LIKELY TRAPPED'
                   write(6,'(a)') 'NEUTONE: ERROR: NEUTRAL AT WALL'//
     >                        ' REMOVED DUE TO TOO MANY REFLECTIONS'//
     >                        ': PARTICLE LIKELY TRAPPED'
                   write(6,'(a,3i6,4g12.5)') 'DATA:',ik,ir,indi,rnew,
     >                                     znew,xvelf,yvelf    

                   write(6,'(a,5i6)') 'DATA:',wlwall1,wlwall2,wltrap1,
     >                                       wltrap2,
     >                          wallpts
                   write(6,*) 'DATA:',
     >                      (id,':',wallpt(indi,id),',',id=1,13)

c
                   SSTRUK = SSTRUK + SPUTY
c
	           ! Add to deposition.
c slmod begin
c... This is causing problems on the extended ITER grid because WALLPT(INDI,18) is
c    zero: -SL, 01/07/2011
c    
                   IF (WALLPT(INDI,18).EQ.0) THEN
                     WRITE(0,*) 'ERROR: WALLPT(INDI,18).EQ.0 in '//
     >                          'NEUTONE.f when dropping trapped neuts'
                     WRITE(6,*) 'ERROR: WALLPT(INDI,18).EQ.0 in '//
     >                          'NEUTONE.f when dropping trapped neuts'
                   ELSE
                     NEROS (INT (wallpt (INDI,18)),1) = 
     >                          NEROS (INT (wallpt (INDI,18)),1) + SPUTY
                   ENDIF
c
c                   NEROS (INT (wallpt (INDI,18)),1) = 
c     >                        NEROS (INT (wallpt (INDI,18)),1) + SPUTY
c slmod end
c
                   if (mtccnt.gt.0) then 
                      MTCSTRUK = MTCSTRUK + SPUTY
                   endif
c
                   IFATE = 4
                   GOTO 899
                ELSE
                   write(0,'(a)') 'NEUTONE: ERROR: NEUTRAL AT WALL'//
     >                        ' REMOVED DUE TO TOO MANY REFLECTIONS'//
     >                        ': PARTICLE LIKELY TRAPPED'
                   write(6,'(a)') 'NEUTONE: ERROR: NEUTRAL AT WALL'//
     >                        ' REMOVED DUE TO TOO MANY REFLECTIONS'//
     >                        ': PARTICLE LIKELY TRAPPED'
                   write(6,'(a,3i6,4g12.5)') 'DATA:',ik,ir,indi,rnew,
     >                                     znew,xvelf,yvelf    

                   write(6,'(a,5i6)') 'DATA:',wlwall1,wlwall2,wltrap1,
     >                                       wltrap2,
     >                          wallpts
                   write(6,*) 'DATA:',
     >                      (id,':',wallpt(indi,id),',',id=1,13)
c
                   SWALLN = SWALLN + SPUTY
c
                   if (mtccnt.gt.0) then 
                      MTCWALLN = MTCWALLN + SPUTY
                   endif
c
                   if (ir.ne.irtrap.and.ir.ne.irwall) then
c
                      write (6,'(a,i6)')
     >                  'NEUTONE: WARNING: NEUTRAL NOT ASSOCIATED'//
     >                        ' WITH WALL RING AT WALL COLLISION',ir
c
                      if (ir.gt.irtrap) then
                         ir = irtrap
                      else
                         ir = irwall
                      endif
                   endif
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
                   GOTO 899
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
C
      GOTO 200
C
c     Particle has finished - one of 6 results - clean up and EXIT
c
  899 CONTINUE
c
c      IF (DEBUGN.OR.(IFATE.EQ.6))
c     >   WRITE (6,9003) iprod,CIST,IK,IR,0,0,R,Z,K,
c     >     VIN,TEMN,SPUTY,ANGLE*RADDEG,IT,FATE(IFATE)
c
c slmod begin
       IF (CPRINT.GT.1) THEN 
        WRITE (6,9003) iprod,CIST,IK,IR,0,0,R,Z,K,
     >     VIN,TEMN,SPUTY,ANGLE*RADDEG,IT,FATE(IFATE)
       ENDIF
c
c       WRITE (6,9003) iprod,CIST,IK,IR,0,0,R,Z,K,
c     >     VIN,TEMN,SPUTY,ANGLE*RADDEG,IT,FATE(IFATE)
c slmod end
c
c     Set return code to ifate
c
      rc = ifate 
c
c     Update lifetime of ion by adding time spent as a neutral.
c
c      cisti = (cistntot + dble(cist)) * fsrate/qtim 
      cisti = (cistntot + cist) * fsrate/qtim 
C
c     Update total time spent on neutrals 
c
      NEUTIM = NEUTIM + ZA02AS (1) - STATIM
c
c     Record MTC statistics 
c
      if (mtccnt.gt.10) then 
         mtctotcnt(11,mtcdat) = mtctotcnt(11,mtcdat) + sputy
      else  
         mtctotcnt(mtccnt,mtcdat) = mtctotcnt(mtccnt,mtcdat) + sputy
      endif
c
      RETURN
c
 9003 FORMAT('L-ONE:',1X,I5,F9.1,4I4,
     >  2F11.7,F9.4,F8.1,F7.2,F5.2,F8.3,I3,:,1X,A)
 9004 FORMAT(//1X,'NEUT DEBUG: DIAGNOSTICS TO BE PRINTED EVERY',I6,
     >  ' TIMESTEPS  (DELTA T =',1P,G10.3,' SECONDS).',//)
 9005 FORMAT(7X,'-NEUT-----TIME--IK--IR--IX--IY',
     >  '-------R----------Z--------K------VIN----TEMN--FRAC--ANGLE--IT',
     >  18('-'))
      END



