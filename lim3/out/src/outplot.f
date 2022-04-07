c
c     These are routines taken from the DIVIMP version of OUT to allow for
c     the implementation of a generalized LOS integration option in LIM.
c     These routines have been adapted to the LIM environment. 
c
c
      SUBROUTINE LOSINT (TVALS,TOUTS,TWIDS,NUMTHE,ROBS,ZOBS,AVPTS,VS,
     >                   MAXDIST,intopt)
      use mod_params
      use mod_comt2
      use mod_comxyt
      use mod_limpoly
      IMPLICIT NONE
C
C  *********************************************************************
C  *                                                                   *
C  *  LOSINT:  INTEGRATE THE VARIABLE VS ALONG A FAN OF SIGHT LINES    *
C  *           FROM A COMMON OBSERVATION POINT (ROBS,ZOBS).  SPATIAL   *
C  *           RESOLUTION IS SIMULATED BY AN AVPTS-POINT AVERAGE OF    *
C  *           A SET OF CHORDS SPANNING THE INTERVAL TWIDS.  AT THE    *
C  *           MOMENT THERE ARE NO WEIGHTS ON THIS AVERAGING WHICH     *
C  *           IMPLIES THAT WE ARE ASSUMING A RECTANGULAR, RATHER      *
C  *           THAN A CIRCULAR, VIEWING CONE.  THE INTEGRAL IS         *
C  *           PERFORMED BY FINDING THE PATH LENGTHS OF THE LOS IN     *
C  *           EACH PLASMA CELL.                                       *
c  *                                                                   *
c  *           INTOPT specifies a LOS integration that will be used    *
c  *                  for each component LOS that is used to calculate *
c  *                  the actual value for each LOS.                   *
c  *                  = 0 - weighted equally - equivalent to a         *
c  *                        rectangular view                           *
c  *                  = 1 - weighted by                                *
c  *                    sqrt(1-((theta-theta_center)/theta_width)**2)  *
c  *                    this should be the equivalent of a circular    *
c  *                    LOS.                                           *
C  *                                                                   *
C  *            CHRIS FARRELL  (HUNTERSKIL)  MARCH 1989                *
C  *            LORNE HORTON   (JET)         JULY  1993                *
C  *            DAVID ELDER    (TORONTO)     MAY   1998                *
C  *            - modifed to optionally return the MAX value in LOS    *
c  *                                                                   *
c  *                                                                   * 
c  *  NOTE: For use in LIM in this routine R=Y and Z=X which may       *
c  *        require transposing                                        *          
c  *        the LIM coordinate system for the code to work correctly   *
C  *        Rvertp = Yvertp   Zvertp = Xvertp                          *
c  *        In general LIM is toroidal while DIVIMP poloidal so the    *
c  *        coordinate systems do not directly overlap                 *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE "PARAMS"
c      include 'params'
c      include 'comxyt'
c      include 'comt2' 
c      include 'limpoly'
c
c      include 'cgeom'
c
      INTEGER NUMTHE,AVPTS,intopt
      REAL    TVALS(NUMTHE),TOUTS(NUMTHE),TWIDS(NUMTHE),
     >        ROBS,ZOBS,VS(MAXNXS,MAXNYS),maxdist
C
      INTEGER I,J,K,IX,IY,NINT,SIDE(2),in
      REAL    THETA,XB(2),WB(2),DIST(2),actdist
      real    weight_factor,total_weight,wfact
c
c     Finding location of maximum along LOS
c
      integer maxswitch,ixmax,iymax
      real maxval 

      integer :: pz=1
      
C     
c     If AVPTS is passed in as < 0 - this is used to indicate that
c     the LOS routine should return the MAXIMUM value of the function
c     along the LOS instead of the integral.
c
      maxswitch = 0
c
      if (avpts.lt.0.0) then
         avpts = abs(avpts)
         maxswitch = 1
      endif
c
c      write (6,*) 'debug:',avpts,maxswitch
c
      XB(1) = ROBS
      XB(2) = ZOBS
C
C  LOOP OVER SIGHT LINES
C
      DO 200 I = 1, NUMTHE
        TVALS(I) = 0.0
c
c       Zero out location of maximum value
c
        maxval= -HI
        ixmax = 0
        iymax = 0 
        total_weight = 0.0
c
C
C  LOOP OVER CHORDS FOR AVERAGING
C
        DO 100 J = 1, AVPTS
          IF (AVPTS.EQ.1) THEN
            THETA = TOUTS(I)
          ELSE
c
c           Modify averaging angles so that the "edge" of the 
c           outside of the averaging cone matches the edge of
c           the specified viewing cone width.
c
c            THETA = TOUTS(I) - 0.5*TWIDS(I) + (J-1)*TWIDS(I)/(AVPTS-1)
c
            THETA = TOUTS(I) - 0.5*TWIDS(I) 
     >                       + (real(J)-0.5)*TWIDS(I)/AVPTS
c
          ENDIF
c
c         Calculate weight_factor for sub-chord
c
          if (intopt.eq.0) then 
             weight_factor = 1.0
          elseif (intopt.eq.1) then   
c
             if (twids(i).gt.0.0) then 
                wfact = 1.0-((theta-touts(i))/(0.5*twids(i)))**2
             else
                wfact = 1.0
             endif
c
             if (wfact.lt.0.0) then 
                weight_factor = 0.0
             else
                weight_factor = sqrt(wfact)
             endif  

          else     
             weight_factor=1.0  
          endif 
c
c         NOTE: If weight factor is zero then do not sum into LOS - cycle
c               the loop at this point
c
          if (weight_factor.le.0.0) cycle

c
c         Sum to total weight
c
          total_weight = total_weight+weight_factor
c
          THETA = THETA*DEGRAD
          WB(1) = COS(THETA)
          WB(2) = SIN(THETA)
C
C  LOOP OVER PLASMA CELLS
C
          DO 20 IX = 1, NXS
c
            DO 10 IY = 1, NYS
              K = KORPG(IX,IY)
c
              if (k.eq.0.or.nvertp(k).eq.0) cycle
c
              CALL INTERS(NVERTP(K),YVERTP(1,K),XVERTP(1,K),
     >                    XB,WB,NINT,DIST,SIDE)
c

              IF (NINT.EQ.2 .AND. DIST(2).GT.0.0) THEN
c
c                write(6,'(a,3i5,12(1x,g12.5))') 'INTERS:',
c     >               ix,iy,k,
c     >               (xvertp(in,k),yvertp(in,k),in=1,4),
c     >               dist(1),dist(2),vs(ix,iy)
c
c
c               NOTE: This may not be the most efficient way to
c               implement the MAX on LOS functionality. This code
c               could be rewritten to return the series of cells
c               along the LOS for the one case and the path
c               length through the cells for the other. However,
c               the interpretation of "MAX" may change so this
c               seems like the best way to implement it for
c               now.
c
                if (maxdist.le.0.0.or.
     >             (maxdist.gt.0.0.and.dist(1).lt.maxdist)) then 
c
c                  Locate and record maximum along LOS.  
c
                   if (vs(ix,iy).gt.maxval) then 
                      maxval = vs(ix,iy)
                      ixmax = ix
                      iymax = iy 
                   endif 
c
c
                   if (maxswitch.eq.0) then
c
c
c                     Assign actual distance to be used for 
c                     integration  
c 
                      if (maxdist.gt.0.0.and.
     >                    dist(2).gt.maxdist) then 
                         actdist = maxdist
                      else
                         actdist = dist(2)
                      endif 
c
                      IF (DIST(1).LT.0.0) THEN
                        TVALS(I)=TVALS(I)+VS(IX,IY)*ACTDIST
     >                                    *weight_factor
                      ELSE
                        TVALS(I)=TVALS(I)+VS(IX,IY)*(ACTDIST-DIST(1))
     >                                    *weight_factor
                      ENDIF
c
                   else
c
c                     Set TVALS to maximum of function.
c
                      tvals(i) = max (tvals(i),abs(vs(ix,iy)))
c         
                   endif  
c
                endif
c
c                write(6,'(a,3i5,5(1x,g13.5))') 'los:',i,ix,iy,
c     >                     tvals(i),vs(ix,iy),
c     >                     dist(2),dist(1),maxdist
c
              ENDIF
   10       CONTINUE
   20     CONTINUE
  100   CONTINUE
c
c       Normalize for LOS integral averaging.
c
        if (maxswitch.eq.0) then
c
c           write(6,'(a,i5,3(1x,g12.5))') 'LOS:',i,tvals(i),total_weight,
c     >                              tvals(i)/total_weight
c
           TVALS(I) = TVALS(I) / total_weight
c
        endif
c
c       Write out information about maximum on LOS 
c
        if (maxswitch.eq.1.and.ixmax.eq.0.or.iymax.eq.0) then 

           write(6,'(a,2i4,1p,8(1x,g12.5))')
     >          'NO LOCATION OF MAXIMUM EMISSION LOS :',ixmax,iymax,
     >          robs,zobs,touts(i)

        elseif (maxswitch.eq.1) then 
           write(6,'(a,2i4,1p,8(1x,g12.5))')
     >          'LOCATION OF MAXIMUM EMISSION ON LOS :',ixmax,iymax,
     >          robs,zobs,touts(i),vs(ixmax,iymax),
     >          xouts(ixmax),youts(iymax),
     >          crnbs(ixmax,iymax,pz),ctembs(ixmax,iymax,pz)    
        endif
c
  200 CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE INTERS (NV,RV,ZV,XB,WB,NINT,DIST,SIDE)
      IMPLICIT NONE
C
C  *********************************************************************
C  *                                                                   *
C  *  INTERS:  FIND THE NUMBER OF INTERSECTIONS, IF ANY, BETWEEN A     *
C  *           LINE DEFINED BY THE POINT XB AND THE DIRECTION COSINES  *
C  *           WB WITH A POLYGON DEFINED BY ITS VERTICES               *
C  *           (RV(I),ZV(I)),I=1,NV.  THE ALGORITHM IS TAKEN FROM THE  *
C  *           NIMBUS ROUTINE G1.  POLYGONS ARE ASSUMED TO BE CONVEX   *
C  *           SO THAT THERE MUST BE EITHER 0 OR 2 INTERSECTIONS.  IN  *
C  *           ADDITION, IF INTERSECTIONS ARE FOUND, THE DISTANCE FROM *
C  *           XB TO EACH IS RETURNED, AS IS THE SIDE OF THE POLYGON   *
C  *           WHICH IS CROSSED IN EACH CASE.                          *
C  *                                                                   *
C  *            LORNE HORTON   (JET)         JULY  1993                *
C  *                                                                   *
C  *********************************************************************
C
      INTEGER NV,NINT,SIDE(2),LR(2)
      REAL    RV(NV),ZV(NV),XB(2),WB(2),DIST(2)
C
      INTEGER I,IP1
      REAL XI,YI,XIP1,YIP1,D
      REAL FNUM,DENO,T,XP,YP
      REAL FMIN,FMAX,TT(2),EPS
      DATA EPS/1.0E-5/
C
      NINT = 0
      TT(1) = 0.0
      DO 120 I = 1, NV
        IF (NINT.EQ.2) GOTO 120
        XI = RV(I)
        YI = ZV(I)
        IP1 = I + 1
        IF (IP1.GT.NV) IP1 = 1
        XIP1 = RV(IP1)
        YIP1 = ZV(IP1)
        D = YIP1 - YI
        IF (D) 90,80,90
   80   FNUM = YI - XB(2)
        DENO = WB(2)
        GOTO 100
   90   D = (XIP1-XI)/D
        FNUM = XI - XB(1) - D*(YI-XB(2))
        DENO = WB(1) - D*WB(2)
  100   IF (DENO) 110,120,110
  110   T = FNUM/DENO
        XP = XB(1) + WB(1)*T
        IF (XI.EQ.XIP1) XP = XI
        YP = XB(2) + WB(2)*T
        IF (YI.EQ.YIP1) YP = YI
        FMIN = AMIN1(XI,XIP1)
        FMAX = AMAX1(XI,XIP1)
        IF (XP.LT.FMIN .OR. XP.GT.FMAX) GO TO 120
        FMIN = AMIN1(YI,YIP1)
        FMAX = AMAX1(YI,YIP1)
        IF (YP.LT.FMIN .OR. YP.GT.FMAX) GO TO 120
        IF (NINT.EQ.1 .AND. ABS(T-TT(1)).LT.EPS) GO TO 120
        NINT = NINT + 1
        TT(NINT) = T
        LR(NINT) = I
  120 CONTINUE
      IF (NINT.EQ.0) GO TO 170
      DIST(1) = TT(1)
      DIST(2) = TT(2)
      SIDE(1) = LR(1)
      SIDE(2) = LR(2)
      IF (DIST(1).LT.DIST(2)) GO TO 170
      DIST(1) = TT(2)
      DIST(2) = TT(1)
      SIDE(1) = LR(2)
      SIDE(2) = LR(1)
  170 CONTINUE
C
      RETURN
      END
c
c
c
      subroutine load_limdata_array(tmpplot,iselect,istate,itype,
     >                           ylab,blab,ref,nizs,ierr)
      use mod_params
      use mod_comtor
      use mod_comt2
      use mod_dynam2
      use mod_dynam3
      use mod_adas_data_spec
      use mod_comxyt
      use mod_limpoly
      use mod_slcom
      implicit none
c      include 'params' 
c      include 'comxyt' 
c      include 'comt2'
c      include 'dynam2'
c      include 'dynam3'
c      include 'comtor'
c      include 'limpoly'
c
c      include 'slcom'
c      include 'adas_data_spec'
c
      real tmpplot(maxnxs,maxnys)
      integer iselect,istate,nizs,ierr,itype      
      character*(*) ylab,blab,ref 
c
c     LOAD_LIMDATA_ARRAY
c
c     This routine loads a 2D LIM array of size MAXNXS,MAXNYS with
c     a quantity specified by the values of iselect and istate. 
c     The allowed values of ISELECT are:
c
c     ITYPE specifies the type of plot - 0 = contour, 1 = integrated
c
c     This list applies to DIVIMP values - ones that are available in LIM
c     at the present time are marked with a "*".
c
c 
c     ISELECT = 1 = TOTAL H POWER LOSS  (W)
c *             2 = TOTAL IMPURITY POWER LOSS  (W)
c               3 = TOTAL POWER LOSS   (W)
c *             4 = SPECIFIED IMPURITY SPECTROSCOPIC LINE 
c                   - NEED TO READ ADAS DATA
c               5 = SPECIFIED HYDROGENIC SPECTROSCOPIC LINE 
c                   - NEED TO READ ADAS DATA
c               6 = PIN Halpha from PINALPHA array
c               7 = PIN HALPHA - By Component from Eirene - 6 for total
c                   - state specifies component
c                     1 - H ionisation
c                     2 - H+ recombination
c                     3 - H2 dissociation
c                     4 - H2+ dissociation
c                     5 - CX of H and H+
c                     6 - TOTAL 
c               8 = PIN HGAMMA - By component from Eirene - 6 for total
c                   - as above 
c               9 = Hydrogen Neutral Density 
c *            10 = Background Plasma Properties
c *                 1 = density
c *                 2 = electron temperature
c *                 3 = ion temperature
c                   4 = velocity
c                   5 = electric field
c *            11 = Impurity Species Density - specified by charge state
c *            12 = Impurity Species Temperature - specified by charge state
c              13 = Impurity Species Velocity - specified by charge state
c              14 = TOTAL H POWER LOSS (W/m3)
c *            15 = TOTAL IMPURITY POWER LOSS (W/m3)
c              16 = TOTAL POWER LOSS (W/m3)
c *            17 = Load PLRP (Particular Line Radiation Profile - see PLRP 
c                   module for istate values.
c              18 = Fluid code Background Plasma Properties
c                   1 = density
c                   2 = electron temperature
c                   3 = ion temperature
c                   4 = velocity
c                   5 = electric field
c              19 = Fluid code Impurity Species Density - specified by charge state
c              20 = Fluid code Impurity Species Temperature - specified by charge state
c              21 = Fluid code Impurity Species Velocity - specified by charge state
c *            22 = SPECIFIED IMPURITY SPECTROSCOPIC LINE AVERAGED TEMPERATURE
c                   - MAY NEED TO READ ADAS DATA
c
c
c     ADAS variables
c
      CHARACTER ADASID*80,graph3*80
      CHARACTER XFESYM*2
      character adasex*3
      integer   adasyr
      INTEGER ISELE,ISELR,ISELX,iseld,ircode
      integer line 
      REAL WLNGTH
c
c     Local variables
c
      real mfact
      integer ix,iy,iz,len,lenstr
      external lenstr

      integer :: pz=1
      
c     
c     Echo input
c
      write(6,'(a,3i5)') 'Loading LIMDATA:',iselect,istate,ierr
c      write(0,'(a,3i5)') 'Loading LIMDATA:',iselect,istate,ierr
      ierr = 0
c
c     Check for valid ISELECT as input
c
      if (iselect.lt.1.or.iselect.gt.22) then 
         write(6,'(a,i5)') 'LOAD_LIMDATA:INVALID SELECTOR:',
     >                       iselect
         ierr = 1
         return
      endif
c
c     Set the YLAB value
c
      call set_ylab(iselect,istate,itype,nizs,ylab) 
c
c     Set the BLAB value
c
      call set_blab(iselect,istate,itype,nizs,blab) 
c
c     Initialize scaling factor
c
      mfact = 1.0
c
c     Initialize data array
c
      call rzero(tmpplot,maxnxs*maxnys)
c
c----------------------------------------------------------
c
c     Hydrogenic power loss - LIM (N/A)
c
c----------------------------------------------------------
c       
      if (iselect.eq.1) then
c
c        Individual states
c
c         BLAB = 'BOLO H POWER LOSS (BOLO)'
c
c         if (istate.eq.0.or.istate.eq.1) then 
c
c            do iy = 1,nrs
c
c               do ix = 1, nks(iy)
c
c                  tmpplot(ix,iy) = tmpplot(ix,iy) 
c     >                              + hpowls(ix,iy,istate)
c     >                              * kareas(ix,iy)
c
c               end do
c
c            end do   
c
c        Total Hydrogenic 
c
c         else
c
c            do iy = 1,nrs
c
c               do ix = 1, nks(iy)
c
c                  do iz = 0,1
c
c                     tmpplot(ix,iy) = tmpplot(ix,iy) 
c     >                              + hpowls(ix,iy,iz)
c     >                              * kareas(ix,iy)
c
c                  end do
c  
c               end do
c
c            end do   
c
c         endif
c
c
c----------------------------------------------------------
c
c     Impurity power loss - LIM - POWLS array
c
c----------------------------------------------------------
c
      elseif (iselect.eq.2) then
c
c        Individual charge state 
c
c
c        Scale by MFACT if required
c
         IF (ABSFAC.GT.0.0) MFACT = MFACT * ABSFAC
c
c         BLAB = 'BOLO IMP POW LOSS'
c
         if (istate.ge.0.and.istate.le.nizs) then 
c
            do iy = 1,nys
c
               do ix = 1, nxs
c
                  tmpplot(ix,iy) = tmpplot(ix,iy) 
     >                              + powls(ix,iy,istate)*mfact
     >                              * kareas(ix,iy)
c
               end do
c
            end do   
c
c        Total Impurity
c
         else 
c
            do iy = 1,nys
c
               do ix = 1, nxs
c
                  do iz = 0,nizs
c
                     tmpplot(ix,iy) = tmpplot(ix,iy) 
     >                              + powls(ix,iy,iz)*mfact
     >                              * kareas(ix,iy)
c
                  end do
c  
               end do
c
            end do   


         endif
c
c----------------------------------------------------------
c
c     Total power loss - Hydrogenic + Impurity - LIM (N/A)
c
c----------------------------------------------------------
c
c      elseif (iselect.eq.3) then
c
c
c         BLAB = 'BOLO TOTAL POW LOSS'
c
c        Scale by MFACT if required
c
c         IF (ABSFAC.GT.0.0) MFACT = MFACT * ABSFAC
c        
c        Hydrogenic
c	 
c	 
c         do iy = 1,nys
c	 
c            do ix = 1, nxs
c	 
c               do iz = 0,1
c	 
c                  tmpplot(ix,iy) = tmpplot(ix,iy) 
c     >                           + hpowls(ix,iy,iz)
c     >                           * kareas(ix,iy)
c	 
c               end do
c  	 
c            end do
c	 
c         end do   
c	 
c        Impurity 
c	 
c         do iy = 1,nys
c	 
c            do ix = 1, nxs
c	 
c               do iz = 0,nizs
c	 
c                  tmpplot(ix,iy) = tmpplot(ix,iy) 
c     >                           + powls(ix,iy,iz)*mfact
c     >                           * kareas(ix,iy)
c	 
c               end do
c  	 
c            end do
c	 
c         end do   
c
c----------------------------------------------------------
c
c     ADAS based - Impurity spectral line 
c
c----------------------------------------------------------
c
      elseif (iselect.eq.4.or.iselect.eq.22) then
c
c         BLAB = 'CODE ADAS IMP PLRP'
c
c        Need to read in ADAS data spec to calculate radiation
c
         if (cadas_switch.eq.0) then  
        
            CALL RDG1 (GRAPH3,ADASID,adasyr,adasex,
     >              ISELE,ISELR,ISELX,ISELD,IERR)
c
c           Save the ADAS data read in into the common block for 
c           possible re-use. Do not set the cadas_switch.
c
            cadasid = adasid
            cadasyr = adasyr
            cadasex = adasex
            cisele  = isele
            ciselr  = iselr
            ciselx  = iselx
            ciseld  = iseld
c
c        Use ADAS data in common instead of reading from input 
c
         elseif (cadas_switch.eq.1) then 
c            
            adasid = cadasid
            adasyr = cadasyr
            adasex = cadasex
            isele  = cisele
            iselr  = ciselr
            iselx  = ciselx
            iseld  = ciseld
c
         endif
c
         if (ierr.ne.0) return 
c
         IF (ISTATE.GE.0.AND.ISTATE.LE.NIZS.AND.ISTATE.LT.CION)THEN
c
            call LDADAS(CION,ISTATE,ADASID,ADASYR,ADASEX,
     >                  ISELE,ISELR,ISELX,
     >                  tmpplot,Wlngth,IRCODE)
c
            IF (IRCODE.NE.0) THEN
               WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
               return   
            ENDIF
c
            if (iselect.eq.4) then   
               REF = 'ADAS PLRP XX XXXXX ('
            elseif (iselect.eq.22) then 
               REF = 'ADAS TEMP XX XXXXX ('
            endif             
c
            WRITE(REF(11:12),'(I2)') ISTATE
            WRITE(REF(14:18),'(I5)') NINT(WLNGTH)
            LEN = LENSTR(REF)
            IF (ISELE.GT.0) REF = REF(1:LEN) // 'E'
            LEN = LENSTR(REF)
            IF (ISELR.GT.0) REF = REF(1:LEN) // 'R'
            LEN = LENSTR(REF)
            IF (ISELX.GT.0) REF = REF(1:LEN) // 'C'
            LEN = LENSTR(REF)
            REF = REF(1:LEN) // ') '
            LEN = LENSTR(REF)
c
         endif 
c
c        Scale by MFACT
c
         IF (ABSFAC.GT.0.0) MFACT = MFACT * ABSFAC
c
         do iy = 1,nys
            do ix = 1,nxs
               tmpplot(ix,iy) = tmpplot(ix,iy) * mfact
            end do 
         end do
c
c        Scale by the impurity temperature for 
c        iselect option 22
c
         if (iselect.eq.22) then 
            do iy = 1,nys
               do ix = 1,nxs
                  tmpplot(ix,iy) = tmpplot(ix,iy) 
     >                           * sdts(ix,iy,istate)
c
               end do 
            end do
         endif
c
c----------------------------------------------------------
c
c     ADAS based - Hydrogenic spectral lines
c
c----------------------------------------------------------
c
      elseif (iselect.eq.5) then
c
c         BLAB = 'CODE ADAS H PLRP'
c
c        Need to read in ADAS data spec to calculate radiation
c
c
         if (cadas_switch.eq.0) then  
        
            CALL RDG1 (GRAPH3,ADASID,adasyr,adasex,
     >              ISELE,ISELR,ISELX,ISELD,IERR)
c
c        Use ADAS data in common instead of reading from input 
c
         elseif (cadas_switch.eq.1) then 
c            
            adasid = cadasid
            adasyr = cadasyr
            adasex = cadasex
            isele  = cisele
            iselr  = ciselr
            iselx  = ciselx
            iseld  = ciseld
c
         endif
c
         if (ierr.ne.0) return 
c
         IF (ISTATE.GE.0 .AND. ISTATE.LE.1) THEN
c
            call LDADAS(1,ISTATE,ADASID,ADASYR,ADASEX,
     >                  ISELE,ISELR,ISELX,
     >                  tmpplot,Wlngth,IRCODE)
c
            IF (IRCODE.NE.0) THEN
               WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
               return   
            ENDIF
c
            REF = 'ADAS H PLRP XX XXXXX ('
            WRITE(REF(13:14),'(I2)') IZ
            WRITE(REF(16:20),'(I5)') NINT(WLNGTH)
            LEN = LENSTR(REF)
            IF (ISELE.GT.0) REF = REF(1:LEN) // 'E'
            LEN = LENSTR(REF)
            IF (ISELR.GT.0) REF = REF(1:LEN) // 'R'
            LEN = LENSTR(REF)
            IF (ISELX.GT.0) REF = REF(1:LEN) // 'C'
            LEN = LENSTR(REF)
            REF = REF(1:LEN) // ') '
            LEN = LENSTR(REF)
c
         endif 
c
c        Scale by MFACT
c
c         IF (ABSFAC.GT.0.0) MFACT = MFACT * ABSFAC
c
         do iy = 1,nys
            do ix = 1,nxs
               tmpplot(ix,iy) = tmpplot(ix,iy) * mfact
            end do 
         end do
c
c----------------------------------------------------------
c
c     PIN Halpha and Hgamma - LIM (N/A)
c
c----------------------------------------------------------
c
c      elseif (iselect.eq.6) then
c
c
c        PIN Halpha - Total only 
c
c         BLAB = 'CODE CODE HALPHA'
c
c        Loop through array
c
c         do iy = 1,nys
c
c            do ix = 1, nxs
c
c               tmpplot(ix,iy) = pinalpha(ix,iy)
c
c            end do
c
c         end do   
c
c
c----------------------------------------------------------
c
c     PIN Halpha and Hgamma - By component from Eirene 
c     LIM (N/A) 
c
c----------------------------------------------------------
c
c      elseif (iselect.eq.7.or.iselect.eq.8) then
c
c        PIN Halpha  
c
c         if (iselect.eq.7) then 
c
c            BLAB = 'CODE CODE HALPHA'
c            line = H_BALPHA
c
c         elseif (iselect.eq.8) then 
c
c            BLAB = 'CODE CODE HGAMMA'
c            line = H_BGAMMA
c
c         endif
c
c        Loop through array
c
c         do iy = 1,nys
c
c            do ix = 1, nxs
c
c               tmpplot(ix,iy) = pinline(ix,iy,istate,line)
c
c            end do
c
c         end do   
c
c----------------------------------------------------------
c
c     PIN Hneutral Density - from Eirene - LIM (N/A)
c
c----------------------------------------------------------
c
c      elseif (iselect.eq.9) then  
c
c
c         do iy = 1,nys
c
c            do ix = 1, nxs
c
c               tmpplot(ix,iy) = pinatom(ix,iy)
c
c            end do
c
c         end do   
c
c----------------------------------------------------------
c
c     DIVIMP Background Plasma Properties - Ne, Te, Ti, Vb, E
c
c----------------------------------------------------------
c
      elseif (iselect.eq.10) then  
c
         do iy = 1,nys
c
            do ix = 1, nxs
c
               if (istate.eq.1) then 
                  tmpplot(ix,iy) = crnbs(ix,iy,pz)
               elseif (istate.eq.2) then 
                  tmpplot(ix,iy) = ctembs(ix,iy,pz)
               elseif (istate.eq.3) then 
                  tmpplot(ix,iy) = ctembsi(ix,iy,pz)
               elseif (istate.eq.4) then 
                  tmpplot(ix,iy) = velplasma(ix,iy,pz)
c               elseif (istate.eq.5) then 
c                  tmpplot(ix,iy) = kes(ix,iy)
               endif
c
            end do
c
         end do   

c
c----------------------------------------------------------
c
c     DIVIMP Impurity Species Densities
c
c----------------------------------------------------------
c


      elseif (iselect.eq.11) then  
c
c        Scaling factor 
c
         IF (ABSFAC.GT.0.0) MFACT = MFACT * ABSFAC
c
         do iy = 1,nys
c
            do ix = 1, nxs
c
               tmpplot(ix,iy) = sdlims(ix,iy,istate)*mfact
c
            end do
c
         end do   

c
c----------------------------------------------------------
c
c     DIVIMP Impurity Species Temperatures
c
c----------------------------------------------------------
c

      elseif (iselect.eq.12) then  
c
         do iy = 1,nys
c
            do ix = 1, nxs
c
               tmpplot(ix,iy) = sdts(ix,iy,istate)
c
            end do
c
         end do   
c
c
c----------------------------------------------------------
c
c     DIVIMP Impurity Species Velocities - LIM (N/A)
c
c----------------------------------------------------------
c
c      elseif (iselect.eq.13) then  
c
c         do iy = 1,nys
c
c            do ix = 1, nxs
c
c               tmpplot(ix,iy) = velavg(ix,iy,istate)
c
c            end do
c
c         end do   
c
c
c----------------------------------------------------------
c
c     Hydrogenic power loss - DIVIMP (W/m3) - LIM (N/A)
c
c----------------------------------------------------------
c       
c      elseif (iselect.eq.14) then
c
c        Individual states
c
c         if (istate.eq.0.or.istate.eq.1) then 
c
c            do iy = 1,nys
c
c               do ix = 1, nxs
c
c                  tmpplot(ix,iy) = tmpplot(ix,iy) 
c     >                              + hpowls(ix,iy,istate)
c
c               end do
c
c            end do   
c
c        Total Hydrogenic 
c
c         else
c
c            do iy = 1,nys
c
c               do ix = 1, nxs
c
c                  do iz = 0,1
c
c                     tmpplot(ix,iy) = tmpplot(ix,iy) 
c     >                              + hpowls(ix,iy,iz)
c
c                  end do
c  
c              end do
c
c            end do   
c
c         endif
c
c----------------------------------------------------------
c
c     Impurity power loss - DIVIMP (W/m3)
c
c----------------------------------------------------------
c
      elseif (iselect.eq.15) then
c
c        Individual charge state 
c
c
c        Scale by MFACT if required
c
         IF (ABSFAC.GT.0.0) MFACT = MFACT * ABSFAC
c
c         BLAB = 'BOLO IMP POW LOSS'
c
         if (istate.ge.0.and.istate.le.nizs) then 
c
            do iy = 1,nys
c
               do ix = 1, nxs
c
                  tmpplot(ix,iy) = tmpplot(ix,iy) 
     >                          + powls(ix,iy,istate)*mfact
c
               end do
c
            end do   
c
c        Total Impurity
c
         else 
c
            do iy = 1,nys
c
               do ix = 1, nxs
c
                  do iz = 0,nizs
c
                     tmpplot(ix,iy) = tmpplot(ix,iy) 
     >                              + powls(ix,iy,iz)*mfact
c
                  end do
c  
               end do
c
            end do   


         endif
c
c----------------------------------------------------------
c
c     Total power loss - Hydrogenic + Impurity - DIVIMP 
c           (W/m3)  - LIM (N/A)
c
c----------------------------------------------------------
c
c      elseif (iselect.eq.16) then
c
c         BLAB = 'BOLO TOTAL POW LOSS'
c
c        Scale by MFACT if required
c
c         IF (ABSFAC.GT.0.0) MFACT = MFACT * ABSFAC
c        
c        Hydrogenic
c	 
c	 
c         do iy = 1,nys
c	 
c            do ix = 1, nxs
c	 
c               do iz = 0,1
c	 
c                  tmpplot(ix,iy) = tmpplot(ix,iy) 
c     >                           + hpowls(ix,iy,iz)
c	 
c               end do
c  	 
c            end do
c	 
c         end do   
c	 
c        Impurity 
c	 
c         do iy = 1,nys
c	 
c            do ix = 1, nxs
c	 
c               do iz = 0,nizs
c	 
c                  tmpplot(ix,iy) = tmpplot(ix,iy) 
c     >                           + powls(ix,iy,iz)*mfact
c	 
c               end do
c  	 
c            end do
c	 
c         end do   
c
c----------------------------------------------------------
c
c     PLRP Calculated by the PLRP module
c
c----------------------------------------------------------
c
      elseif (iselect.eq.17) then  
c
         do iy = 1,nys
c
            do ix = 1, nxs
c
               tmpplot(ix,iy) = plrps(ix,iy,istate) * absfac
c
            end do
c
         end do   
c
c
c----------------------------------------------------------
c
c    FLUID CODE Solution Background Plasma Properties - Ne, Te, Ti, Vb, E
c
c    LIM (N/A)
c
c----------------------------------------------------------
c
c      elseif (iselect.eq.18) then  
c
c         do iy = 1,nys
c
c            do ix = 1, nxs
c
c               if (istate.eq.1) then 
c                  tmpplot(ix,iy) = e2dnbs(ix,iy)
c               elseif (istate.eq.2) then 
c                  tmpplot(ix,iy) = e2dtebs(ix,iy)
c               elseif (istate.eq.3) then 
c                  tmpplot(ix,iy) = e2dtibs(ix,iy)
c               elseif (istate.eq.4) then 
c                  tmpplot(ix,iy) = e2dvhs(ix,iy)
c               elseif (istate.eq.5) then 
c                  tmpplot(ix,iy) = e2des(ix,iy)
c               endif
c
c            end do
c
c         end do   
c
c----------------------------------------------------------
c
c     FLUID CODE Impurity Species Densities
c     LIM (N/A) 
c
c----------------------------------------------------------
c
c
c      elseif (iselect.eq.19) then  
c
c         do iy = 1,nys
c
c            do ix = 1, nxs
c
c               tmpplot(ix,iy) = e2dnzs(ix,iy,istate)
c
c            end do
c
c         end do   
c
c----------------------------------------------------------
c
c     FLUID CODE Impurity Species Temperatures
c     LIM (N/A)
c
c----------------------------------------------------------
c
c      elseif (iselect.eq.20) then  
c
c        Note: This assumes that impurity temperature is 
c              equal to the ion temperature for a fluid 
c              code. 
c
c         do iy = 1,nys
c
c            do ix = 1, nxs
c
c               tmpplot(ix,iy) = e2dtebs(ix,iy)
c
c            end do
c
c         end do   
c
c
c----------------------------------------------------------
c
c     FLUID CODE Impurity Species Velocities
c     LIM (N/A)
c
c----------------------------------------------------------
c
c      elseif (iselect.eq.21) then  
c
c         do iy = 1,nys
c
c            do ix = 1, nxs
c
c               tmpplot(ix,iy) = e2dvzs(ix,iy,istate)
c
c            end do
c
c         end do   
c
      endif
c
c
c
      return
      end
c
c
c
      subroutine set_ylab(iselect,istate,itype,nizs,ylab)
      implicit none
      integer iselect,istate,itype,nizs
      character*(*) ylab
c
c     SET_YLAB:
c
c      
c     This is a support routine to the 2D DIVIMP data loading and 
c     integration code. Depending on the values of iselect,istate and 
c     itype - this routine sets the y axis label to something 
c     reasonable.
c
c     Itype specifies the type of plot - 0 = contour, 1 = integrated
c
c  
      integer len,lenstr
      external lenstr
c
c----------------------------------------------------------
c     Hydrogen power loss
c----------------------------------------------------------
c       
      if (iselect.eq.1) then
c
         if (istate.eq.0) then
            YLAB = 'H-NEUTRAL POW LOSS (BOLO)'
         elseif (istate.eq.1) then 
            YLAB = 'H+ POW LOSS (BOLO)'
         else 
            YLAB = 'TOT H POW LOSS (BOLO)'
         endif 
c
         len = lenstr(ylab)
c
         if (itype.eq.0) then 
            ylab = ylab(1:len) // ' (W)'
         elseif (itype.eq.1) then 
            ylab = ylab(1:len) // ' (W*M)'
         endif
c
c----------------------------------------------------------
c     Impurity power loss 
c----------------------------------------------------------
c
      elseif (iselect.eq.2) then
c
c        Individual charge state 
c
         if (istate.ge.0.and.istate.le.nizs) then 
c
            write(ylab,'(''IMP POW LOSS IZ ='',i4)')
     >                                          istate
c
c        Total Impurity
c
         else 
c
            YLAB = 'TOT IMP POW LOSS (BOLO) '
c
         endif 
c
         len = lenstr(ylab)
c
         if (itype.eq.0) then 
            ylab = ylab(1:len) // ' (W)'
         elseif (itype.eq.1) then 
            ylab = ylab(1:len) // ' (W*M)'
         endif
c
c----------------------------------------------------------
c     Total power loss 
c----------------------------------------------------------
c
      elseif (iselect.eq.3) then
c
         YLAB = 'TOT H+IMP POW LOSS (BOLO)'
c
         len = lenstr(ylab)
c
         if (itype.eq.0) then 
            ylab = ylab(1:len) // ' (W)'
         elseif (itype.eq.1) then 
            ylab = ylab(1:len) // ' (W*M)'
         endif
c
c----------------------------------------------------------
c     Line Radiation:
c          iselect = 4 = ADAS impurity line
c          iselect = 5 = ADAS hydrogen line
c          iselect = 6 = PIN total Halpha
c          iselect = 7 = EIRENE Halpha  istate =6 = total
c          iselect = 8 = EIRENE Hgamma  istate =6 = total
c          iselect = 17= PLRP radiation profile
c----------------------------------------------------------
c
      elseif (iselect.eq.4.or.iselect.eq.5.or.
     >        iselect.eq.6.or.iselect.eq.7.or.
     >        iselect.eq.8.or.iselect.eq.17) then
c
         YLAB = 'LINE RADIATION (PHOTONS'
c
         len = lenstr(ylab)
c
         if (itype.eq.0) then 
            ylab = ylab(1:len) // '/M^3/S)'
         elseif (itype.eq.1) then 
            ylab = ylab(1:len) // '/M^2/S)'
         elseif (itype.eq.2) then 
            ylab = ylab(1:len) // '/M^2/S/SR)'
         endif
c
c----------------------------------------------------------
c
c     PIN - Hydrogen Neutral Density 
c
c----------------------------------------------------------
c
      elseif (iselect.eq.9) then   

         YLAB = 'DENSITY ('
c
         len = lenstr(ylab)
c
         if (itype.eq.0) then 
            ylab = ylab(1:len) // '/M^3)'
         elseif (itype.eq.1) then 
            ylab = ylab(1:len) // '/M^2)'
         endif
c
c----------------------------------------------------------
c
c     DIVIMP - Background Plasma Properties
c
c----------------------------------------------------------
c
      elseif (iselect.eq.10) then   

         if (istate.eq.1) then 
            YLAB = 'DENSITY (M^-3)'
         elseif(istate.eq.2) then 
            YLAB = 'ELEC TEMPERATURE (eV)'
         elseif(istate.eq.3) then 
            YLAB = 'ION TEMPERATURE (eV)'
         elseif(istate.eq.4) then 
            YLAB = 'VELOCITY (M/S)'
         elseif(istate.eq.5) then 
            YLAB = 'ELECTRIC FIELD (V/M(?))'
         endif
c
c----------------------------------------------------------
c
c     DIVIMP - Impurity Density
c
c----------------------------------------------------------
c
      elseif (iselect.eq.11) then   

         write(YLAB,'(''IMP DENSITY: STATE='',i4,
     >                ''(M^-3)'')') istate
c
c
c----------------------------------------------------------
c
c     DIVIMP - Impurity Temperature
c
c----------------------------------------------------------
c
      elseif (iselect.eq.12) then   

         write(YLAB,'(''IMP TEMPERATURE: STATE='',i4,
     >                ''(eV)'')') istate
c
c
c----------------------------------------------------------
c
c     DIVIMP - Impurity Velocity
c
c----------------------------------------------------------
c
      elseif (iselect.eq.13) then   

         write(YLAB,'(''IMP VELOCITY: STATE='',i4,
     >                ''(M/S)'')') istate

c
c----------------------------------------------------------
c     Hydrogen power loss  (W/m3)
c----------------------------------------------------------
c       
      elseif (iselect.eq.14) then
c
         if (istate.eq.0) then
            YLAB = 'H-NEUTRAL POW LOSS (BOLO)'
         elseif (istate.eq.1) then 
            YLAB = 'H+ POW LOSS (BOLO)'
         else 
            YLAB = 'TOTAL H POW LOSS (BOLO)'
         endif 
c
         len = lenstr(ylab)
c
         if (itype.eq.0) then 
            ylab = ylab(1:len) // ' (W/M3)'
         elseif (itype.eq.1) then 
            ylab = ylab(1:len) // ' (W/M2)'
         elseif (itype.eq.2) then 
            ylab = ylab(1:len) // ' (W)'
         endif
c
c----------------------------------------------------------
c     Impurity power loss (W/m3)
c----------------------------------------------------------
c
      elseif (iselect.eq.15) then
c
c        Individual charge state 
c
         if (istate.ge.0.and.istate.le.nizs) then 
c
            write(ylab,'(''IMP POW LOSS IZ ='',i4)')
     >                                          istate
c
c        Total Impurity
c
         else 
c
            YLAB = 'TOTAL IMP POW LOSS (BOLO) '
c
         endif 
c
         len = lenstr(ylab)
c
         if (itype.eq.0) then 
            ylab = ylab(1:len) // ' (W/M3)'
         elseif (itype.eq.1) then 
            ylab = ylab(1:len) // ' (W/M2)'
         elseif (itype.eq.2) then 
            ylab = ylab(1:len) // ' (W)'
         endif
c
c----------------------------------------------------------
c     Total power loss  (W/m3)
c----------------------------------------------------------
c
      elseif (iselect.eq.16) then
c
         YLAB = 'TOT H+IMP POW LOSS (BOLO)'
c
         len = lenstr(ylab)
c
         if (itype.eq.0) then 
            ylab = ylab(1:len) // ' (W/M3)'
         elseif (itype.eq.1) then 
            ylab = ylab(1:len) // ' (W/M2)'
         elseif (itype.eq.2) then 
            ylab = ylab(1:len) // ' (W)'
         endif


c
c----------------------------------------------------------
c
c     FLUID CODE - Background Plasma Properties
c
c----------------------------------------------------------
c
      elseif (iselect.eq.18) then   

         if (istate.eq.1) then 
            YLAB = 'FC DENSITY (M^-3)'
         elseif(istate.eq.2) then 
            YLAB = 'FC ELEC TEMPERATURE (eV)'
         elseif(istate.eq.3) then 
            YLAB = 'FC ION TEMPERATURE (eV)'
         elseif(istate.eq.4) then 
            YLAB = 'FC VELOCITY (M/S)'
         elseif(istate.eq.5) then 
            YLAB = 'FC ELECTRIC FIELD (V/M(?))'
         endif
c
c----------------------------------------------------------
c
c     FLUID CODE - Impurity Density
c
c----------------------------------------------------------
c
      elseif (iselect.eq.19) then   

         write(YLAB,'(''FC IMP DENSITY: STATE='',i4,
     >                ''(M^-3)'')') istate
c
c
c----------------------------------------------------------
c
c     FLUID CODE - Impurity Temperature
c
c----------------------------------------------------------
c
      elseif (iselect.eq.20) then   

         write(YLAB,'(''FC IMP TEMPERATURE: STATE='',i4,
     >                ''(eV)'')') istate
c
c
c----------------------------------------------------------
c
c     FLUID CODE - Impurity Velocity
c
c----------------------------------------------------------
c
      elseif (iselect.eq.21) then   

         write(YLAB,'(''FC IMP VELOCITY: STATE='',i4,
     >                ''(M/S)'')') istate

c
c----------------------------------------------------------
c
c     DIVIMP - Emission Weighted impurity temperature
c
c----------------------------------------------------------

      elseif (iselect.eq.22) then 
c
         YLAB = 'EMISSION WEIGHTED AV. ION TEMP (eV)'
c
      endif
c
c
c
      return
      end
c
c
c
      subroutine set_blab(iselect,istate,itype,nizs,blab)
      implicit none
      integer iselect,istate,itype,nizs
      character*(*) blab
c
c     SET_BLAB:
c
c      
c     This is a support routine to the 2D DIVIMP data loading and 
c     integration code. Depending on the values of iselect,istate and 
c     itype - this routine sets the plot label to something 
c     reasonable.
c
c     Itype specifies the type of plot - 0 = contour, 1 = integrated
c
c  
      integer len,lenstr
      external lenstr
c
c----------------------------------------------------------
c     Hydrogen power loss
c----------------------------------------------------------
c       
      if (iselect.eq.1) then
c
         if (itype.eq.0) then           
            BLAB = 'H POW LOSS (BOLO)'
         elseif (itype.eq.1) then 
            BLAB = 'CODE H POW LOSS (BOLO)'
         endif
c
c----------------------------------------------------------
c     Impurity power loss 
c----------------------------------------------------------
c
      elseif (iselect.eq.2) then 
c
         if (itype.eq.0) then           
            BLAB = 'IMP POW LOSS (BOLO)'
         elseif (itype.eq.1) then 
            BLAB = 'BOLO IMP POW LOSS'
         endif
c
c----------------------------------------------------------
c     Total power loss 
c----------------------------------------------------------
c
      elseif (iselect.eq.3) then 
c
         if (itype.eq.0) then           
            BLAB = 'TOTAL POW LOSS (BOLO)'
         elseif (itype.eq.1) then 
            BLAB = 'BOLO TOTAL POW LOSS'
         endif
c  
c----------------------------------------------------------
c     ADAS IMPURITY PLRP 
c----------------------------------------------------------
c
      elseif (iselect.eq.4) then  
c
         if (itype.eq.0) then           
            BLAB = 'ADAS IMP PLRP'
         elseif (itype.eq.1.or.itype.eq.2) then 
            BLAB = 'CODE ADAS IMP PLRP'
         endif
c  
c----------------------------------------------------------
c     ADAS HYDROGENIC PLRP 
c----------------------------------------------------------
c
      elseif (iselect.eq.5) then 
c   
         if (itype.eq.0) then           
            BLAB = 'ADAS HYDROGENIC PLRP'
         elseif (itype.eq.1.or.itype.eq.2) then 
            BLAB = 'CODE ADAS H PLRP'
         endif
c  
c----------------------------------------------------------
c     PIN TOTAL HALPHA
c----------------------------------------------------------
c
      elseif (iselect.eq.6) then 
c   
         if (itype.eq.0) then           
            BLAB = 'CODE HALPHA'
         elseif (itype.eq.1.or.itype.eq.2) then 
            BLAB = 'CODE CODE HALPHA'
         endif
c  
c----------------------------------------------------------
c     EIRENE HALPHA
c----------------------------------------------------------
c
      elseif (iselect.eq.7) then 
c   
         if (istate.eq.6) then 
            if (itype.eq.0) then           
               BLAB = 'EIRENE TOTAL HALPHA'
            elseif (itype.eq.1.or.itype.eq.2) then 
               BLAB = 'CODE EIRENE TOT HALPHA'
            endif
         else
            if (itype.eq.0) then           
               write(blab,'(a,i4)') 'EIRENE HALPHA COMP='
            elseif (itype.eq.1.or.itype.eq.2) then 
               write(blab,'(a,i4)') 'CODE EIRENE HALPHA COMP='
            endif
         endif
c  
c----------------------------------------------------------
c     EIRENE HGAMMA 
c----------------------------------------------------------
c
      elseif (iselect.eq.8) then 
c   
         if (istate.eq.6) then 
            if (itype.eq.0) then           
               BLAB = 'EIRENE TOTAL HGAMMA'
            elseif (itype.eq.1.or.itype.eq.2) then 
               BLAB = 'CODE EIRENE TOT HGAMMA'
            endif
         else
            if (itype.eq.0) then           
               write(blab,'(a,i4)') 'EIRENE HGAMMA COMP='
            elseif (itype.eq.1.or.itype.eq.2) then 
               write(blab,'(a,i4)') 'CODE EIRENE HGAMMA COMP='
            endif
         endif
c
c
c----------------------------------------------------------
c     PIN - Hydrogen Neutral Density 
c----------------------------------------------------------
c
      elseif (iselect.eq.9) then   

         BLAB = 'PIN NEUTRAL H DENSITY'
c
c
c----------------------------------------------------------
c
c     DIVIMP - Background Plasma Properties
c
c----------------------------------------------------------
c
      elseif (iselect.eq.10) then   

         if (istate.eq.1) then 
            BLAB = 'BG ION DENSITY'
         elseif(istate.eq.2) then 
            BLAB = 'BG ELECTRON TEMPERATURE'
         elseif(istate.eq.3) then 
            BLAB = 'BG ION TEMPERATURE'
         elseif(istate.eq.4) then 
            BLAB = 'BG VELOCITY'
         elseif(istate.eq.5) then 
            BLAB = 'BG ELECTRIC FIELD'
         endif
c
c----------------------------------------------------------
c
c     DIVIMP - Impurity Density
c
c----------------------------------------------------------
c
      elseif (iselect.eq.11) then   

         write(BLAB,'(''IMP DENSITY: STATE='',i4,
     >                ''(M^-3)'')') istate
c
c
c----------------------------------------------------------
c
c     DIVIMP - Impurity Temperature
c
c----------------------------------------------------------
c
      elseif (iselect.eq.12) then   

         write(BLAB,'(''IMP TEMPERATURE: STATE='',i4,
     >                ''(eV)'')') istate
c
c
c----------------------------------------------------------
c
c     DIVIMP - Impurity Velocity
c
c----------------------------------------------------------
c
      elseif (iselect.eq.13) then   

         write(BLAB,'(''IMP VELOCITY: STATE='',i4,
     >                ''(M/S)'')') istate

c
c----------------------------------------------------------
c     Hydrogen power loss (W/m3)
c----------------------------------------------------------
c       
      elseif (iselect.eq.14) then
c
         if (itype.eq.0) then           
            BLAB = 'H POW LOSS (BOLO)'
         elseif (itype.eq.1.or.itype.eq.2) then 
            BLAB = 'CODE H POW LOSS (BOLO)'
         endif
c
c----------------------------------------------------------
c     Impurity power loss 
c----------------------------------------------------------
c
      elseif (iselect.eq.15) then 
c
         if (itype.eq.0) then           
            BLAB = 'IMP POW LOSS (BOLO)'
         elseif (itype.eq.1.or.itype.eq.2) then 
            BLAB = 'CODE IMP POW LOSS'
         endif
c
c----------------------------------------------------------
c     Total power loss 
c----------------------------------------------------------
c
      elseif (iselect.eq.16) then 
c
         if (itype.eq.0) then           
            BLAB = 'TOTAL POW LOSS (BOLO)'
         elseif (itype.eq.1.or.itype.eq.2) then 
            BLAB = 'CODE TOTAL POW LOSS'
         endif
c
c
      elseif (iselect.eq.17) then  
c
         if (itype.eq.0) then           
            BLAB = 'CUSTOM IMP PLRP'
         elseif (itype.eq.1.or.itype.eq.2) then 
            BLAB = 'CODE CUSTOM IMP PLRP'
         endif
c
c----------------------------------------------------------
c
c     FLUID CODE - Background Plasma Properties
c
c----------------------------------------------------------
c
      elseif (iselect.eq.18) then   

         if (istate.eq.1) then 
            BLAB = 'FC ION DENSITY'
         elseif(istate.eq.2) then 
            BLAB = 'FC ELECTRON TEMPERATURE'
         elseif(istate.eq.3) then 
            BLAB = 'FC ION TEMPERATURE'
         elseif(istate.eq.4) then 
            BLAB = 'FC VELOCITY'
         elseif(istate.eq.5) then 
            BLAB = 'FC ELECTRIC FIELD'
         endif
c
c----------------------------------------------------------
c
c     FLUID CODE - Impurity Density
c
c----------------------------------------------------------
c
      elseif (iselect.eq.19) then   

         write(BLAB,'(''FC IMP DENSITY: STATE='',i4,
     >                ''(M^-3)'')') istate
c
c
c----------------------------------------------------------
c
c     FLUID CODE - Impurity Temperature
c
c----------------------------------------------------------
c
      elseif (iselect.eq.20) then   

         write(BLAB,'(''FC IMP TEMPERATURE: STATE='',i4,
     >                ''(eV)'')') istate
c
c
c----------------------------------------------------------
c
c     FLUID CODE - Impurity Velocity
c
c----------------------------------------------------------
c
      elseif (iselect.eq.21) then   

         write(BLAB,'(''FC IMP VELOCITY: STATE='',i4,
     >                ''(M/S)'')') istate


c  
c----------------------------------------------------------
c     DIVIMP - EMISSION WEIGHTED AVERAGE ION TEMPERATURE 
c----------------------------------------------------------
c
      elseif (iselect.eq.4) then  
c
         BLAB = 'ADAS-BASED AVERAGE ION TEMP'
c
      endif
c
c
c
      return
      end
c
c
c
      subroutine calc_lim_poly
      use mod_params
      use mod_comtor
      use mod_comxyt
      use mod_limpoly
      implicit none
c      include 'params'
c      include 'comxyt'
c      include 'comtor'
c      include 'limpoly'
c
c     CALC_LIM_POLY:
c
c     This routine calculates "polygons" for the LIM grid
c     by calculating the corner points for each bin that are
c     used for the density computation. This will simplify the
c     interface to all the DIVIMP based grid routines that 
c     may be added since the LIM "grid" will then be 
c     represented in exactly the same format even though it 
c     starts off as just a series of cell bounds in the X and 
c     Y directions. This change should make it possible to 
c     even include routines similar to gridpos if that would be
c     useful.
c
c
c     Other useful existing quantities:
c
c     XOUTS,YOUTS: cell center points
c     XWIDS,YWIDS: cell widths in X and Y
c     XS   ,YS   : cell boundaries 
c 
      integer ix,iy,ip
      real lxb,uxb,lyb,uyb
c
      npolyp = nxs*nys
c  
c     Lower bound of XS = CAW
c     Upper bound of XS = CA
c
c     Lower bound of YS = 0.0
c     Upper bound of YS = CL
c
c     The upper bounds are loaded into the arrays automatically
c     when LIM starts.  
c
      ip = 0
c
c      do ix = 1,nxs
c         write(6,'(a,i5,3(1x,g12.5))') 'XS:',ix,xs(ix),
c     >               xouts(ix),xwids(ix)
c      end do
c      do iy = 1,nys
c         write(6,'(a,i5,3(1x,g12.5))') 'YS:',iy,ys(iy),
c     >                       youts(iy),ywids(iy)
c      end do
c
      do ix = 1,nxs
c 
         if (ix.eq.1) then  
            lxb = caw
            uxb = xs(1)
         else
            lxb = xs(ix-1)
            uxb = xs(ix)
         endif
c
         do iy = 1,nys
c
            if (iy.eq.1) then  
c
               lyb = 0.0
               uyb = ys(1)
c
            else
c 
               lyb = ys(iy-1)
               uyb = ys(iy) 
c
            endif
c
c           Calculate polygon corners
c
            ip = ip+1  
            korpg(ix,iy) = ip
            nvertp(ip) = 4
c 
            xvertp(1,ip) = lxb
            xvertp(2,ip) = uxb
            xvertp(3,ip) = uxb
            xvertp(4,ip) = lxb
c
            yvertp(1,ip) = lyb
            yvertp(2,ip) = lyb
            yvertp(3,ip) = uyb
            yvertp(4,ip) = uyb
c
c           NOTE: The actual cell area depends on some LIM 
c                 options that may not be passed to OUT. 
c           1) Slab or Cylindrical geometry option - stored in XCYLS
c              - specified by IGEOM in the LIM input file              
c              - correct cell volumes for the different geometries
c           2) Plasma elongation 
c              - specified by CKO and CKI in the LIM input file
c              - correction factor stored in DELPS 
c
c           At the present time these are not included here 
c           since it will only affect cases which use these options           
c
            kareas(ix,iy) = xwids(ix) * ywids(iy)
c
         end do
c
      end do
c 
      return
      end
C
C
C
      SUBROUTINE LDADAS(CZ,IZ,ADASID,ADASYR,ADASEX,ISELE,ISELR,ISELX,
     >                  CVALS,WAVE,IRCODE)
      use mod_params
      use mod_comtor
      use mod_dynam2
      use mod_comt2
      use mod_comxyt
      use mod_pindata
      IMPLICIT NONE
C
C  *********************************************************************
C  *                                                                   *
C  *  LDADAS:  CODE TO LOAD THE REQUESTED LINE EMISSION PROFILE        *
C  *           INTO THE MATRIX CVALS.  THE LINE IS SPECIFIED BY THE    *
C  *           NUCLEAR CHARGE OF THE EMITTING ION, THE                 *
C  *           CHARGE OF THE EMITTING IONISATION STATE, AN ID FLAG     *
C  *           WHICH LOCATES THE INPUT FILE, AND THREE BLOCK SELECTOR  *
C  *           NUMBERS, ONE EACH FOR EMISSION BY ELECTRON EXCITATION,  *
C  *           RECOMBINATION FROM THE NEXT HIGHER IONISATION STATE,    *
C  *           AND CHARGE EXCHANGE FROM THE NEXT HIGHER IONISATION     *
C  *           STATE.  IN ADDITION TO THE EMISSION PROFILE, THE        *
C  *           ROUTINE RETURNS THE WAVELENGTH OF THE TRANSITION AND    *
C  *           AN ERROR CODE FROM THE ADAS EXTRACTION ROUTINE, SPEC.   *
C  *           NOTE THAT THERE IS NO CHECKING OF THE SELECTOR NUMBERS! *
C  *                                                                   *
C  *                                                                   *
C  *            LORNE HORTON   (JET)         SEPTEMBER 1993            *
C  *                                                                   *
C  *********************************************************************
C
c      include 'params'
c      include 'comxyt'
c      INCLUDE 'cgeom'
c      include 'pindata'
c      include 'dynam2'
c      include 'comtor'
c      include 'comt2'
C
      INTEGER   CZ,IZ,ADASYR,ISELE,ISELR,ISELX,IRCODE
      REAL      WAVE,CVALS(MAXNXS,MAXNYS)
      CHARACTER ADASID*(*),ADASEX*(*)
C
      INTEGER   IX,IY,IADAS,NPAIRS,IKK
      REAL*8    TADAS(20),DADAS(20)
      REAL*8    WLNGTH,PECAE(20),PECAR(20),PECAX(20)
      LOGICAL*4 LTRNG(20),LDRNG(20)
      CHARACTER ADASGR*8,ADASTY*80,PECTITLE*120
      CHARACTER XFESYM*2
c
c     Define RIZB
c
      real rizb

      integer :: pz = 1
C     
      call rzero (cvals,MAXNXS*MAXNYS)
c
      rizb = cizb
      WAVE = 0.0
      IRCODE = 0
      CALL XXUID(ADASID)
      IF (ADASYR.GE.0) THEN
        ADASGR = 'pec??#'//XFESYM(CZ)
        WRITE(ADASGR(4:5),'(I2.2)') ADASYR
      ELSE
        ADASGR = '*'
      ENDIF
      ADASTY = '*'
      CALL XXSPEC(ADASGR,ADASTY,ADASEX)
C
      DO IY = 1,NYS
        DO IX = 1,NXS,20
          NPAIRS = MIN0(20,NXS-(IX-1))
          DO IADAS = 1,NPAIRS
            TADAS(IADAS) = DBLE(CTEMBS(IX+(IADAS-1),IY,pz))
            DADAS(IADAS) = DBLE(1.E-6*RIZB*CRNBS(IX+(IADAS-1),IY,pz))
          ENDDO
C
          CALL DZERO(PECAE,NPAIRS)
          IF (ISELE.GT.0) THEN
            CALL SPEC(ISELE,IZ,CZ,NPAIRS,TADAS,DADAS,
     >                WLNGTH,PECAE,LTRNG,LDRNG,PECTITLE,IRCODE)
            IF (IRCODE.NE.0) RETURN
          ELSE IF (ISELE.EQ.-1) THEN
C
C  JUST LOAD EMISSION MEASURE FOR ISEL = -1
C    - SINCE THIS VALUE CAN EXCEED THE UNIX SINGLE PRECISION LIMIT,
C      WORK IN DENSITY UNITS OF 10**18
C
            CALL DINIT(PECAE,NPAIRS,1.D6*1.D-36)
            WLNGTH = 0.0
          ENDIF
C
          CALL DZERO(PECAR,NPAIRS)
          IF (ISELR.GT.0) THEN
            CALL SPEC(ISELR,IZ,CZ,NPAIRS,TADAS,DADAS,
     >                WLNGTH,PECAR,LTRNG,LDRNG,PECTITLE,IRCODE)
            IF (IRCODE.NE.0) RETURN
          ELSE IF (ISELR.EQ.-1) THEN
            CALL DINIT(PECAR,NPAIRS,1.D6*1.D-36)
            WLNGTH = 0.0
          ENDIF
C
C  FOR IMPURITIES USE THE THIRD SWITCH FOR CX, FOR HYDROGEN
C  ADD IN A MOLECULAR CONTRIBUTION INSTEAD
C
          CALL DZERO(PECAX,NPAIRS)
CLDH - USE ION TEMPERATURE FOR CX RATE
          IF (CZ.GT.1.0) THEN
            DO IADAS = 1,NPAIRS
              TADAS(IADAS) = DBLE(CTEMBSI(IX+(IADAS-1),IY,pz))
            ENDDO
          ENDIF
CLDH
          IF (ISELX.GT.0) THEN
            CALL SPEC(ISELX,IZ,CZ,NPAIRS,TADAS,DADAS,
     >                WLNGTH,PECAX,LTRNG,LDRNG,PECTITLE,IRCODE)
            IF (IRCODE.NE.0) RETURN
          ELSE IF (ISELX.EQ.-1) THEN
            CALL DINIT(PECAX,NPAIRS,1.D6*1.D-36)
            WLNGTH = 0.0
          ENDIF
C
          DO IADAS = 1,NPAIRS
            IKK = IX + (IADAS-1)
c
c            write(6,'(a,3i5,10(1x,g13.5))') 'LDADAS:', ikk,iy,iz,
c     >         rizb,crnbs(ikk,ir),ktebs(ikk,iy),sdlims(ikk,iy,iz),
c     >         sdlims(ikk,ir,iz+1),pinatom(ikk,iy),
c     >         pecae(iadas),pecar(iadas),pecax(iadas)
c
            IF (CZ.GT.1.0) THEN
              CVALS(IKK,IY) = 1.D-6*
     >        (PECAE(IADAS) * RIZB *CRNBS(IKK,IY,pz)*SDLIMS(IKK,IY,IZ)
     >        +PECAR(IADAS) * RIZB *CRNBS(IKK,IY,pz)*SDLIMS(IKK,IY,IZ+1)
     >        +PECAX(IADAS) * PINATOM(IKK,IY) * SDLIMS(IKK,IY,IZ+1))
c
c
c              if (cvals(ikk,iy).gt.0.0) then 
c                write(6,'(a,3i5,8(1x,g12.5))') 'LDADAS1:',ikk,iy,iadas,
c     >                pecae(iadas),pecar(iadas),pecax(iadas),
c     >                CRnbs(ikk,iy),pinatom(ikk,iy),sdlims(ikk,iy,iz),
c     >                sdlims(ikk,iy,iz+1),cvals(ikk,iy)
c              endif 
c 
            ELSE
C
C---- HYDROGEN DENSITIES ARE IN DIFFERENT ARRAYS
C
              CVALS(IKK,IY) = 1.D-6*
     >        (PECAE(IADAS) * RIZB * CRNBS(IKK,IY,pz) * PINATOM(IKK,IY)
     >        +PECAR(IADAS) * RIZB * CRNBS(IKK,IY,pz) * CRNBS(IKK,IY,pz)
     >        +PECAX(IADAS) * CRNBS(IKK,IY,pz)   * PINMOL(IKK,IY))

c
c              if (cvals(ikk,iy).gt.0.0) then 
c                write(6,'(a,3i5,8(1x,g12.5))') 'LDADAS2:',ikk,iy,iadas,
c     >                pecae(iadas),pecar(iadas),pecax(iadas),
c     >                CRnbs(ikk,iy),pinatom(ikk,iy),pinmol(ikk,iy),
c     >                cvals(ikk,iy)
c              endif 
c

            ENDIF
c
c            if (cvals(ikk,iy).gt.0.0)  then 
c              write(6,'(a,3i5,10(1x,g13.5))') 'LDADAS:', ikk,iy,iz,
c     >         rizb,crnbs(ikk,iy),ctembs(ikk,iy),sdlims(ikk,iy,iz),
c     >         sdlims(ikk,iy,iz+1),pinatom(ikk,iy),
c     >         pecae(iadas),pecar(iadas),pecax(iadas),cvals(ikk,iy)
c            endif
c
          ENDDO
        ENDDO
      ENDDO
C
      WAVE = WLNGTH
C
      RETURN
      END
c
C
C
      subroutine adjustout(touts,numthe,zadj,robs,zobs)
      use mod_params
      implicit none
c      include 'params'
      integer numthe
      real touts(numthe),zadj,robs,zobs
c
c     This subroutine maps the THETA values in the Touts array from a
c     specific observation position onto a projected R co-ordinate plane
c     at value of Z specified by Z-adjustment.
c
c     David Elder     1995, June 8.
c
      integer i,j
      real    angadj
c
      if (zobs.gt.zadj) then
         angadj = 270.0
      elseif (zobs.lt.zadj) then
         angadj = 90.0
      else
         write (6,*) 'ERROR in mapping THTEA to R plotting coordinates:'
         write (6,*) 'Observation Z =',zobs,
     >               ' is the same as adjustment plane = ',zadj
         return
      endif
c
      do i = 1,numthe
c
c        Adjust Touts
c
         touts(i) = robs + (zobs-zadj) * tan((touts(i)-angadj)*degrad)
c
      end do
c
      return
      end
c
c
c
      subroutine adjustoutz(touts,numthe,radj,robs,zobs)
      use mod_params
      implicit none
c      include 'params'
      integer numthe
      real touts(numthe),radj,robs,zobs
c
c     This subroutine maps the THETA values in the Touts array from a
c     specific observation position onto a projected Z co-ordinate plane
c     at value of R specified by R-adjustment.
c
c     David Elder     1995, June 8.
c
      integer i,j
      real    angadj
c
      if (robs.gt.radj) then
         angadj = 270.0
      elseif (robs.lt.radj) then
         angadj = 90.0
      else
         write (6,*) 'ERROR in mapping THTEA to Z plotting coordinates:'
         write (6,*) 'Observation R =',robs,
     >               ' is the same as adjustment plane = ',radj
         return
      endif
c
      do i = 1,numthe
c
c        Adjust Touts
c
         touts(i) = zobs + (robs-radj) * tan((touts(i)-angadj)*degrad)
c
      end do
c
      return
      end
c
c
c
      subroutine calc_expt(iseld,touts,tvals,maxnthe,numthe,maxdatx,
     >                     themin,themax,maxnngs,ngs,datatitle)
      implicit none
c
c      include 'params'
c
      integer iseld,ngs,numthe,maxnthe,maxnngs,maxdatx
      real themin,themax,touts(maxnthe),tvals(maxnthe,maxnngs)
      character*(*) datatitle
c
c     CALC_EXPT: This subroutine loads the selected
c                experimental data set from disk into
c                the array tvals at the position ngs. It
c                interpolates the experimental data to
c                obtain values at the same points as will
c                be plotted for the DIVIMP data - this
c                was done to enable proper scaling of all
c                the data when the plots are being drawn.
c
c                This data is then loaded into the array
c                TVALS for passing to DRAW. NGS is the
c                plot index to be used by the first set of
c                experimental data.
c
c
c                David Elder, November 24, 1997
c
c
c     Local variables
c
      integer dataunit,maxcols
      parameter (dataunit=13,maxcols=1)
c
c      integer maxdatx,dataunit,maxcols
c      parameter (maxdatx=1000,dataunit=13,maxcols=1)
c
      integer axis_type,in,cur_index,num_expt,ipos,ncols
      integer ind1,ind2
c
c
c     Experimental data
c
      real expt_axis(maxdatx),expt_data(maxdatx),theta
c
c     Write out input -
c
      write(6,*) 'Input:',iseld,maxnthe,numthe,maxnngs,ngs
c
      call load_expt_data(dataunit,iseld,expt_axis,axis_type,expt_data,
     >                    maxcols,maxdatx,num_expt,ncols,
     >                    datatitle)
c
      if (num_expt.le.0) then
c
         write(6,*) 'ERROR IN EXPERIMENTAL DATA: CAN NOT'//
     >                    ' LOAD DATASET # ',
     >                    iseld, ' - NO ELEMENTS FOUND'
         write(0,*) 'ERROR IN EXPERIMENTAL DATA: CAN NOT'//
     >                    ' LOAD DATASET # ',
     >                    iseld, ' - NO ELEMENTS FOUND'
c
         ngs = ngs - 1
c
         return
c
      endif
c
c     Check for single point LOS comparisons
c
      if (numthe.eq.1.and.num_expt.eq.1) then
c
         tvals(1,ngs) = expt_data(1)
         write (6,*) 'Calc_expt: SINGLE:',tvals(1,ngs)
c
      else
c
         do in = 1, numthe
c
            theta = touts(in)
c
            call arrpos(theta,expt_axis,num_expt,ind1,ind2)
c
            if (ind1.eq.-1) then
c
c              Index below position of first experimental data point
c              Set to zero ...
c
               tvals(in,ngs) = 0.0
c
            elseif (ind2.eq.-1) then
c
c              Out of bounds at top of range - again set to zero
c
               tvals(in,ngs) = 0.0
c
            else
c
c              Value in range - perform linear interpolation.
c
               tvals(in,ngs) = expt_data(ind1) +
     >                     (expt_data(ind2)-expt_data(ind1))
     >                    *(theta-expt_axis(ind1))
     >                    /(expt_axis(ind2)-expt_axis(ind1))
c
            endif
c
c           Write out debugging information to start ...
c
c slmod begin - new
             write (6,'(a,3i4,4(1x,g12.5),a,
c
c            write (6,'(a,2i4,4(1x,g12.5),a,
c slmod end
     >                    3(1x,g12.5))') 'Calc_expt:',in,ind1,ind2,
     >              expt_data(ind1),tvals(in,ngs),
     >              expt_data(ind2),tvals(in,1),'|',
     >              expt_axis(ind1), theta,expt_axis(ind2)
c
c
c
c            if (cur_index.eq.1) then
c
c              Index below position of first experimental data point
c              Set to zero ...
c
c               tvals(in,ngs) = 0.0
c
c            elseif (in.eq.num_expt.and.theta.gt.expt_axis(num_expt))
c     >               then
c
c              Out of bounds at top of range - again set to zero
c
c               tvals(in,ngs) = 0.0
c
c            else
c
c              Value in range - perform linear interpolation.
c
c               tvals(in,ngs) = expt_data(cur_index-1) +
c     >                     (expt_data(cur_index)-expt_data(cur_index-1))
c     >                    *(theta-expt_axis(cur_index-1))
c     >                    /(expt_axis(cur_index)-expt_axis(cur_index-1))
c
c            endif
c
c           Write out debugging information to start ...
c
c slmod begin - new
c            IF (cur_index.GT.1)
c     .        write (6,'(a,2i4,4(1x,g12.5),a,
c
c            write (6,'(a,2i4,4(1x,g12.5),a,
c slmod end
c     >                    3(1x,g12.5))') 'Calc_expt:',in,cur_index,
c     >              expt_data(cur_index-1),tvals(in,ngs),
c     >              expt_data(cur_index),tvals(in,1),'|',
c     >              expt_axis(cur_index-1), theta,expt_axis(cur_index)
c
         end do
c
      endif
c
c     Exit
c
      return
c
      end
c
c
c
      subroutine load_expt_data(dataunit,iseld,expt_axis,axis_type,
     >                    expt_data,maxcols,maxdatx,
     >                    num_expt,ncols,datatitle)
c slmod begin - new
c
c Input:
c ISELD         - data index number
c MAXCOLS       - maximum number of data columns
c MAXDATX       - maximum number of data items
c
c Output:
c EXPT_AXIS     - independent data
c AXIS_TYPE     - read from file 13 (optional listing), default is 1
c               = 1 
c               - 6 = PSIN data   
c  
c
c EXPT_DATA     - MAXDATX,MAXCOLS dependent data
c NUM_EXPT      - number of data items in the data block, default is 0
c NCOLS         - number of data columns in data block, default is 1
c DATATITLE     - title, default is 'NO TITLE'
c
c 
c
c
c slmod end
      implicit none
      integer dataunit,axis_type,num_expt,iseld,maxdatx,colindex
      integer maxcols,ncols
      character*(*) datatitle
      real expt_data(maxdatx,maxcols),expt_axis(maxdatx)
c
c     LOAD_EXPT_DATA: This routine loads the experimental data
c                     into the local arrays and passes
c                     back the relevant information.
c
c
c     Local variables
c
      character*200 buffer
      real xval,yval(maxcols),scalef
      real axis_shift,expt_norm  
      integer extstr,len,start,in,startn,endn,lenstr,colcnt
      external extstr,lenstr
      integer dataset_num,dataset_cnt,total_dataset_cnt
c
c     The following quantities were introduced to allow adjustments
c     to invlaid experimental values without actually changing the
c     contents of the file.
c
c     cutoff_val - this is for the elimination of spikes - any entries
c                  that are greater than this value will be set to zero.
c     offset_val - this is to graphically compensate for non-zero
c                  calibrartion offsets so plot comparison is easier.
c
      integer maxvals_cutoff
      parameter(maxvals_cutoff=25)
      real cutoff_val(maxvals_cutoff), offset_val,minexpt,maxexpt

      REAL    HI,LO
      PARAMETER( HI=1.E37 ,LO=1.E-37)
c
c
c     Initialize
c
      scalef = 1.0
      axis_shift = 0.0
      expt_norm = 0.0
      startn = 1
      endn = maxdatx
      offset_val = 0.0
      call rinit(cutoff_val,maxvals_cutoff,HI)
      maxexpt = -HI   
      minexpt = HI
      axis_type = 1
      num_expt = 0
      dataset_cnt = 0
      ncols = 1
      datatitle = 'NO TITLE'
c
      write (6,*) 'ISELD:',iseld
c
      rewind (dataunit)
c
c     Scan through data unit looking for INDEX keyword
c     with the appropriate data tag number.
c
 100  read(dataunit,'(a200)',end=500,err=500) buffer
c
c      len = lenstr(buffer)
c      write(6,*) 'BUFF:',buffer(1:len),':'
c
c     Ignore Empty Lines
c
c      if (buffer.eq.'') goto 100
c
c     Ignore comments
c
      if (buffer(1:1).eq.'$') goto 100
c
c     Check number of datasets in file
c
      if (buffer(1:10).eq.'FILECOUNT:') then
         read (buffer(11:),*) total_dataset_cnt
c
         write (6,*) 'DATASETS:',total_dataset_cnt
c
         if (iseld.gt.total_dataset_cnt) then
            write (6,*) 'REQUESTED EXPERIMENTAL DATA SET # ',
     >                  iseld,' DOES NOT EXIST'
            write (6,*) 'ONLY ', total_dataset_cnt,' DATA SETS ARE '//
     >                       'SPECIFIED BY FILECOUNT:'
            write (0,*) 'REQUESTED EXPERIMENTAL DATA SET # ',
     >                  iseld,' DOES NOT EXIST'
            write (0,*) 'ONLY ', total_dataset_cnt,' DATA SETS ARE '//
     >                       'SPECIFIED BY FILECOUNT:'
         endif
      endif
c
c     Look for INDEX and count for datasets
c
      if (buffer(1:6).eq.'INDEX:') then
         dataset_cnt = dataset_cnt+1
c
         read(buffer(7:),*) dataset_num
c
c         write (6,*) 'INDEX:',dataset_cnt,dataset_num,iseld
c
c        Indexing error - write out error message and continue
c
         if (dataset_cnt.ne.dataset_num) then
c
            write (6,*) 'ERROR IN EXPERIMENTAL DATA'//
     >               ' FILE: DATASET # ',dataset_num, ' IS IN'//
     >                  ' FILE POSITION ',DATASET_CNT
C
            write (0,*) 'ERROR IN EXPERIMENTAL DATA'//
     >               ' FILE: DATASET # ',dataset_num, ' IS IN'//
     >                  ' FILE POSITION ',DATASET_CNT
C
         endif
c
         if (dataset_num.eq.iseld.or.dataset_cnt.eq.iseld) goto 300
c
      endif
c
c     Loop back and read more data until exit
c
      goto 100
c
c     Index or count for dataset found - continue processing
c
 300  continue
c
c     Read in the rest of the HEADER block until data found - then
c     start processing data until EOF or next INDEX.
c
 350  read(dataunit,'(a200)',end=400,err=400) buffer
c
c      len = lenstr(buffer)
c      write(6,*) 'DATA:',len,':',buffer(1:len),':'
c
c     Ignore Empty lines
c
c      if (buffer.eq.'') goto 350
c
c     Ignore comments
c
      if (buffer(1:1).eq.'$') goto 350
c
c     Exit if Next INDEX is found
c
      if (buffer(1:6).eq.'INDEX:') goto 400
c
c     Extract a Title if one is specified
c
      if (buffer(1:6).eq.'TITLE:') then
         len = extstr(buffer(7:),start)
         datatitle = buffer(8:LEN_TRIM(buffer))
c         datatitle = buffer(6+start:len)
         write(6,*) 'TITLE:',datatitle,':'
      endif
c
c     Extract Scaling Factor if specified
c
      if (buffer(1:7).eq.'SCALEF:') read(buffer(8:),*) scalef
c
c     Extract Axis_type if given
c
      if (buffer(1:5).eq.'AXIS:') read(buffer(6:),*) axis_type
c
c     Extract Number of colums of data if given - one assumed
c
      if (buffer(1:6).eq.'NCOLS:') then
c
         read(buffer(7:),*) ncols
c
         write(6,*) 'NCOLS:',ncols
c
         if (ncols.gt.maxcols) then
c
c           Issue error message and only load the first column of
c           data
c
            write (6,*) 'REQUESTED EXPERIMENTAL DATA SET # ',
     >                  iseld,' HAS MORE COLUMNS THAN MAX =',maxcols
            write (0,*) 'REQUESTED EXPERIMENTAL DATA SET # ',
     >                  iseld,' HAS MORE COLUMNS THAN MAX =',maxcols
            ncols = 1
c
         endif
      endif
c
c     Extract Offset Value if specified
c
      if (buffer(1:7).eq.'OFFSET:') read(buffer(8:),*) offset_val
c
c     Extract Axis Shift if given
c
      if (buffer(1:6).eq.'SHIFT:') read(buffer(7:),*) axis_shift
c
c     Extract Normalization value if given
c
      if (buffer(1:5).eq.'NORM:') read(buffer(6:),*) expt_norm
c
c     Extract Cutoff Value if specified
c
      if (buffer(1:7).eq.'CUTOFF:') read(buffer(8:),*) 
     >                (cutoff_val(in),in=1,min(ncols,maxvals_cutoff))
c
c     Extract data limit counters if present
c
      if (buffer(1:6).eq.'COUNT:') then
         read(buffer(7:),*) startn,endn
         write (6,*) 'COUNT:',startn,endn
      endif
c
c     Other headers will be ignored - check for data and load it.
c
c     Data lines start with a blank
c
      if (buffer(1:1).eq.' ') then
c
         if(lenstr(buffer).gt.6) then
c
c
c            write (0,*) 'BUFFER:',buffer,':'
c
            read(buffer,*) in,xval,
     >                        (yval(colcnt),colcnt=1,ncols)
c
            if (in.ge.startn.and.in.le.endn) then
c
               expt_axis(in-startn+1) = xval
c
c               if (axis_type.eq.5) then
c
c                  expt_axis(in-startn+1) = expt_axis(in-startn+1)
c     >                                     /360.0 * 2.0 * PI
c               endif
c
               do colcnt = 1,ncols
                  expt_data(in-startn+1,colcnt) = yval(colcnt)*scalef
               end do
c
               num_expt = in-startn+1
c
c              Keep track of minimum value in first experimental data.
c
               minexpt = min(minexpt,expt_data(in-startn+1,1))
c
c              Keep track of maximum value in first experimental data.
c
               maxexpt = max(maxexpt,expt_data(in-startn+1,1))
c
c               write (6,'(A,2I6,5G12.4)') 'NUM:',
c     .           num_expt,in,xval,yval(1),yval(2),yval(3),scalef
c
            endif
c
         endif
c
      endif
c
      goto 350
c
c     Continue and wrap up file processing - experimental data
c     has been read.
c
 400  continue
c
c     Data cleanup processing ...
c
c     Check experimental data for cutoff and adjust by offset.
c
      if (offset_val.eq.-1.0) offset_val = minexpt
c
      write (6,*) 'OFFSET_VAL:',offset_val
c
      do in = 1,num_expt
c
         if (axis_shift.ne.0.0) expt_axis(in) = expt_axis(in) 
     >                          + axis_shift
c
         do colcnt = 1,ncols

            expt_data(in,colcnt) = expt_data(in,colcnt) - offset_val
c
            if (expt_norm.ne.0.0.and.maxexpt.ne.0.0) then 
               expt_data(in,colcnt) = 
     >                  expt_data(in,colcnt)/maxexpt * expt_norm
            endif
c
            if (expt_data(in,colcnt).gt.
     >          cutoff_val(min(colcnt,maxvals_cutoff)))
     >                                  expt_data(in,colcnt) = 0.0

         end do
c
      end do
c
      return
c
 500  continue
c
c     Error exit condition - DATA SET NOT FOUND
c
      write (6,*) 'ERROR IN EXPERIMENTAL DATA: EOF OR DATASET # ',
     >             iseld,' NOT FOUND'
      write (0,*) 'ERROR IN EXPERIMENTAL DATA: EOF OR DATASET # ',
     >             iseld,' NOT FOUND'
c
c     Exit
c
      return
      end
c
c
c
      subroutine arrpos(val,vals,nvals,ind1,ind2)
      implicit none
c
      integer nvals,ind1,ind2 
      real val,vals(nvals) 
c
c     ARRPOS: This routine determines the indices of the data in the
c             array vals which are the upper and lower bounds of the
c             input value. The vals array must be ordered - either
c             ascending or descending is fine.  
c
      integer in,ipos,jpos
      external ipos,jpos
c
c
c     Check boundaries as special cases
c
      if (val.eq.vals(1)) then 
         ind1 = 1
         ind2 = 2
         return
      endif
c
      if (val.eq.vals(nvals)) then 
         ind1 = nvals-1
         ind2 = nvals
         return
      endif
c
c     Find specific cell
c
      if (vals(1).lt.vals(nvals)) then 
         in = ipos(val,vals,nvals)   
         if (in.eq.1.and.val.lt.vals(1)) then 
            ind1 = -1
            ind2 =  1 
         elseif (in.eq.nvals.and.val.gt.vals(nvals)) then 
            ind1 = nvals
            ind2 = -1 
         else
            ind1 = in -1
            ind2 = in  
         endif 

      else
         in = jpos(val,vals,nvals)
         if (in.eq.1.and.val.gt.vals(1)) then 
            ind1 = -1
            ind2 =  1 
         elseif (in.eq.nvals.and.val.lt.vals(nvals)) then 
            ind1 = nvals
            ind2 = -1 
         else
            ind1 = in 
            ind2 = in +1
         endif 


      endif
c 
      return
c
      end   
