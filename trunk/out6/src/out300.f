c     -*-Fortran-*-
c
      subroutine out300(iref,graph,iopt,ierr)
      use mod_collector_probe
      implicit none
      integer iref,iopt,ierr
      character*(*) graph
c
      include 'params'
      include 'outcom'
c
c     Other common blocks
c
      include 'cgeom'
      include 'comtor'
c      include 'cneut2'
c      include 'dynam2'
      include 'dynam3'
c      include 'dynam4'
      include 'pindata'
c      include 'cadas'
c      include 'grbound'
c      include 'outxy'
      include 'cedge2d'
      include 'transcoef'
      include 'line_profile'
c
      include 'expt_data'
c
c      include 'cioniz'
c      include 'reiser' 
c      include 'printopt' 
c
c     Local Variables
c
c
      integer ik,ir,iz
      integer in,ii
c
c     Local Variables
c

c
c     For plot 326
c
      real minfrac,maxfrac
      integer axis_type,plot_type  

c
c     RCP/OSM probe variables
c
      integer osmvals, osmplots,rcpvals,rcpplots
      integer int_type,icnt
      real    r1p,z1p,r2p,z2p,rsect,zsect
      real    axmin,axmax
      integer exp_ds, exp_vcalcopt,exp_tcalcopt,exp_dataopt
      real    exp_param, exp_offset,exp_offsets(4)
c
c     collector probe
c
      real probe_diameter,probe_dperp
      integer axis_opt
c
c     For LP - line profiles plots - 345, 346
c
      real maxmod
      integer lp_plot_type, lp_plot_avg,lp_iexpt
      real axis_calibration,lp_axis_offset
      real expt_axis_offset
c
      integer ncell_add
      real low_bin_width,high_bin_width,lp_out,lp_cnt

c
c     Other variables
c
      character infile*100
c
      integer iz1



c
c     Input for generalized plots in 300 series
c
      integer iselect,istate,iexpt,iavg,iaxis
      integer ifact,navg,npts,minsteps
      integer iflag
      integer icntr_opt,ncntr_opt
      real    optval,stepsize
      logical plotr

C
C---- VARIABLES FOR 350+ SERIES PLOTS
C
      INTEGER MAV,IMAV,KM,I3
      REAL CX1,CY1,CZ1,CX2,CY2,CZ2,DOT,TH3,CTHC,CTHN
      REAL C1,C2,VNOR,THN,DETM,FM,FKM,PHI,TPHI
      REAL CMATR(2,3)
      REAL VN(3),VP(3),VT(3)
      REAL COUTS(MAXCH3,MAXCH3,3),WRES(MAXCH3,MAXCH3)

      REAL GFF



c
c     For multiple plots on the page - series 700
c
      real mvals(maxdatx,maxplts,maxngs)
      real mouts(maxdatx,maxplts,maxngs)
      real mwids(maxdatx,maxplts,maxngs)
      integer  ringnos(maxplts),ip,ig,ip2,pnks(maxplts,maxngs)
      integer  mdrawtype(maxplts,maxngs)
      integer  drawtype(maxngs)  
      integer  axistype
      integer  sctype,ngrm
      integer  pngs(maxplts)  
      character*36 mlabs(maxplts,maxngs)
      character*36 pltlabs(maxplts)
c
c     For plot 357
c
      integer plotid,nexpt,maxexpt
      parameter (maxexpt=10)
      integer expt_ds(maxexpt)
      integer expt_col(maxexpt)
c slmod begin
      REAL CalcPressure,GetCs

      INTEGER optflow,i1,i2,i4
      REAL    frac
      CHARACTER*128 cdum1,graph6
c slmod end
c      
      integer in1

      REAL TOUTS2(MAXTHE),TVALS2(MAXTHE,MAXNGS)
c
c     Arrays for averaging  
c
      real av_outs(maxdatx),av_vals(maxdatx,maxngs),
     >     av_cnt(maxdatx)
      integer num_av,tmp_av
      integer av_ref(maxdatx),id 
      logical exclude_pt,exp_av
c
      IF (IREF.LT.310) THEN
        ref = graph(5:41)
        CALL RDG2 (GRAPH2,ROBS,ZOBS,THEMIN,DTHE,THERES,NUMTHE,
     >              IZMIN,IZMAX,AVPTS,NUMSMOOTH,ATYPE,IERR)
c
c       Make adjustment for ASDEX UPGRADE plots ... their angular
c       coordinate system is 180 degrees out of phase with JET.
c       i.e. Angles measured from -ve R-axis ... so add 180 to THEMIN
c       for ASDEX plots.
c
c       Not all SONNET grids are Asdex Upgrade grids - remove this for
c       now - replace with input parameter that will allow for
c       an arbitrary offset. (implement when needed)
c
c        if (cgridopt.eq.3) themin = themin + 180.0
c
        IF (IERR.EQ.1) THEN
           WRITE(6,*) 'RDG2 ERROR READING 200+SERIES- GRAPH DETAILS'
           IERR = 0
           RETURN
        ENDIF
        IF (IOPT.EQ.0) RETURN
c
c     Put in new plots - read auxilliary data if required.  
      ELSEIF (IREF.LT.350) THEN
c
c         WRITE(0,*) 'IOPT=',iopt,ierr
c
c        Set ierr = 0 since not all plots in this section will call
c        a secondary input routine and should not inherit an old value
c        of IERR. 
c
         ierr = 0
c
         if (iopt.eq.0) return 
c
         if (iref.eq.311.or.iref.eq.313.or.iref.eq.315) then 
c
            call rdfn_bolo(graph2,infile,iseld,iflag,ierr)
c
c        Contour
c
         elseif (iref.eq.321) then 
c
            call rdg_contour(graph2,iselect,istate,
     >                   iexpt,optval,ierr) 
c
c        Along ring
c 
         elseif (iref.eq.326) then 
c
            call rdg_ring(graph2,iselect,istate,
     >                    iexpt,minfrac,maxfrac,axis_type,
     >                    plot_type,ierr) 
c
c        LOS
c
         elseif (iref.eq.331.or.iref.eq.333) then 
c
            call rdg_los(graph2,npts,navg,iselect,istate,
     >                   iexpt,iaxis,iavg,ifact,optval,ierr) 
c
         elseif (iref.eq.335) then 
c
            call rdg_los3D(graph2,iselect,istate,
     >                   iexpt,iaxis,minsteps,stepsize,
     >                   optval,ierr) 
c
         elseif (iref.eq.341) then 
c
            call rdg_xsection(graph2,r1p,z1p,r2p,z2p,npts,iselect,
     >                        istate,iexpt,iavg,ierr)
c
         endif 
c
         if (ierr.ne.0) return        
c
c
c     RCP probe plots
c
      elseif (iref.lt.380) then
c
        if (iopt.eq.0) return
c
        if (iref.eq.361) then 

           call rdgcol(graph3,r1p,z1p,r2p,z2p,probe_diameter,
     >                 probe_dperp,axis_opt,ierr)
           if (ierr.ne.0) then
               WRITE(6,*) 'RDGCOL ERROR READING COLLECTOR'//
     >                    ' PROBE GRAPH DETAILS'
               IERR = 0
               RETURN
           endif
           
        else
c
           call rdgrcp(graph3,r1p,z1p,r2p,z2p,int_type,exp_ds,
     >              exp_offsets,exp_dataopt,exp_vcalcopt,
     >              exp_tcalcopt,
     >              exp_param,ierr)
c
           exp_offset = exp_offsets(1) 
c
           if (ierr.eq.1) then
               WRITE(6,*) 'RDGRCP ERROR READING RCP GRAPH DETAILS'
               IERR = 0
               RETURN
           endif
        endif
c
      ELSEIF (IREF.LT.400) THEN
         REF  = GRAPH(5:41)
         CALL RDG3 (GRAPH3,ROBS,ZOBS,CX1,CY1,CZ1,CX2,CY2,CZ2,NUMTHE,
     >              THERES,IZMIN,IZMAX,AVPTS,NUMSMOOTH,ATYPE,IERR)
         IF (IERR.EQ.1) THEN
            WRITE(6,*) 'RDG3 ERROR READING 390+SERIES- GRAPH DETAILS'
            IERR = 0
            RETURN
         ENDIF
         IF (IOPT.EQ.0) RETURN

      endif  





      call init_plot(iref,graph,iopt)







c
c-----------------------------------------------------------------------
c

      IF (IREF.EQ.300) THEN
C
C     LINE-INTEGRATED BREMSSTRAHLUNG EMISSION
c     FREE-FREE ONLY
C
C     VERIFY PARAMETERS
C
         plotr = .false.
         if (numsmooth.lt.0) then
            plotr = .true.
            numsmooth = abs(numsmooth)
         endif
c
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
C           WRITE(6,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
C  IZMIN AND IZMAX IGNORED
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
C
C  CREATE THETA VECTOR
C
         DO II = 1, NUMTHE
           TOUTS(II) = THEMIN + (II-1)*DTHE
           TWIDS(II) = THERES
         ENDDO
         THEMAX = TOUTS(NUMTHE)
C
C  LABELS AND SCALE FACTORS
C
         YLAB = 'BREMS. (PH M-2 S-1 NM-1)'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)// '- NEW METHOD'
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(''RESOLUTION='',F12.6,'' DEGREES'')') THERES
         IF (ATYPE.EQ.0) THEN
           MFACT = 1.0
           WRITE(ANLY,'(''NO SCALE FACTOR APPLIED'')')
         ELSEIF (ATYPE.EQ.1) THEN
           MFACT = DTHE / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = '',G12.6,'' / (2*PI)'')')
     >                DTHE
         ELSEIF (ATYPE.EQ.2) THEN
           MFACT = 1.0 / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (2*PI)'') ')
         ELSEIF (ATYPE.EQ.3) THEN
           MFACT = 1.0 / ( 4.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (4*PI)'')')
           YLAB = 'BREMS. (PH M-2 S-1 NM-1 SR-1)'
         ENDIF
c
C  ASSUME BREMSSTRAHLUNG IS MEASURED AT 5235 A (523.5 nm) FOR NOW
c
         WLNGTH = 5235.0
C
C  LOAD INTEGRAND
C
         CALL LDBRFF(wlngth,plastmp,IRCODE,NIZS)
c
c delete
c
c
c         DO IR = 1,NRS
c           DO IK = 1,NKS(IR)
c             GFF = 1.456*(2.368/KTEBS(IK,IR))**(-0.149)
c             PLASTMP(IK,IR) = 9.55E-21*ZEFFS(IK,IR,3)*GFF
c     >                       *(RIZB*KNBS(IK,IR))**2
c     >                       /(SQRT(KTEBS(IK,IR))*WLNGTH)
c           ENDDO
c         ENDDO
c
         CALL RVALKR (CVALSA,PLASTMP,1,1,1,FT,FP,MFACT,
     >     XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
C
C  INTEGRATE
C
         CALL LOSINT(TVALS(1,1),TOUTS,TWIDS,NUMTHE,
     >               ROBS,ZOBS,AVPTS,CVALSA,0.0,0)
C
         WRITE (IPLOT,9012) NPLOTS,REF
         WRITE (IPLOT,*) NVIEW
         WRITE (IPLOT,*) PLANE
         WRITE (IPLOT,*) ANLY
c
c         Adjust TOUTS to against R instead of Theta if required.
c
         if (plotr) then
            call adjustout(touts,numthe,zadj,robs,zobs)
            themin = touts(1)
            themax = touts(numthe)
            write(XLAB,'(''RADIUS (M) AT Z='',f6.3)') zadj
         endif
C
C  DON'T PLOT SINGLE LINES OF SIGHT
C
         IF (NUMTHE.GT.1) THEN
c
            if (iseldef.ne.0) then

               ngs = 2

               call calc_expt(iseldef,touts,tvals,maxthe,numthe,
     >                   themin,themax,maxngs,ngs,datatitle)

               ELABS(1) = BLABS
               ELABS(2) = 'EXPT'//DATATITLE
c
               CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,ngs,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)
c
            else
c
               CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,1,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,BLABS,REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)
c
            endif
c
         ELSE
           WRITE(iplot,*) TOUTS(1),TVALS(1,1)
         ENDIF
C
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
C
c
c-----------------------------------------------------------------------
c
C
      ELSEIF (IREF.EQ.301) THEN
C
C     LINE-INTEGRATED BREMSSTRAHLUNG EMISSION - ASSUMING NO IMPURITIES
c     FREE-FREE ONLY
C
C     VERIFY PARAMETERS
C
         plotr = .false.
         if (numsmooth.lt.0) then
            plotr = .true.
            numsmooth = abs(numsmooth)
         endif
c
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
C           WRITE(6,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
C  IZMIN AND IZMAX IGNORED
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
C
C  CREATE THETA VECTOR
C
         DO II = 1, NUMTHE
           TOUTS(II) = THEMIN + (II-1)*DTHE
           TWIDS(II) = THERES
         ENDDO
         THEMAX = TOUTS(NUMTHE)
C
C  LABELS AND SCALE FACTORS
C
         YLAB = 'BREMS. (PH M-2 S-1 NM-1)'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)// '- NEW METHOD'
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(''RESOLUTION='',F12.6,'' DEGREES'')') THERES
         IF (ATYPE.EQ.0) THEN
           MFACT = 1.0
           WRITE(ANLY,'(''NO SCALE FACTOR APPLIED'')')
         ELSEIF (ATYPE.EQ.1) THEN
           MFACT = DTHE / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = '',G12.6,'' / (2*PI)'')')
     >                DTHE
         ELSEIF (ATYPE.EQ.2) THEN
           MFACT = 1.0 / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (2*PI)'') ')
         ELSEIF (ATYPE.EQ.3) THEN
           MFACT = 1.0 / ( 4.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (4*PI)'')')
           YLAB = 'BREMS. (PH M-2 S-1 NM-1 SR-1)'
         ENDIF
c
C  ASSUME BREMSSTRAHLUNG IS MEASURED AT 5235 A (523.5 nm) FOR NOW
c
         WLNGTH = 5235.0
C
C  LOAD INTEGRAND - NO IMPURITIES - USE NIZS=-1
C
         CALL LDBRFF(wlngth,plastmp,IRCODE,-1)
c
c delete
c
c
c         DO IR = 1,NRS
c           DO IK = 1,NKS(IR)
c             GFF = 1.456*(2.368/KTEBS(IK,IR))**(-0.149)
c             PLASTMP(IK,IR) = 9.55E-21*ZEFFS(IK,IR,3)*GFF
c     >                       *(RIZB*KNBS(IK,IR))**2
c     >                       /(SQRT(KTEBS(IK,IR))*WLNGTH)
c           ENDDO
c         ENDDO
c
         CALL RVALKR (CVALSA,PLASTMP,1,1,1,FT,FP,MFACT,
     >     XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
C
C  INTEGRATE
C
         CALL LOSINT(TVALS(1,1),TOUTS,TWIDS,NUMTHE,
     >               ROBS,ZOBS,AVPTS,CVALSA,0.0,0)
C
         WRITE (IPLOT,9012) NPLOTS,REF
         WRITE (IPLOT,*) NVIEW
         WRITE (IPLOT,*) PLANE
         WRITE (IPLOT,*) ANLY
c
c         Adjust TOUTS to against R instead of Theta if required.
c
         if (plotr) then
            call adjustout(touts,numthe,zadj,robs,zobs)
            themin = touts(1)
            themax = touts(numthe)
            write(XLAB,'(''RADIUS (M) AT Z='',f6.3)') zadj
         endif
C
C  DON'T PLOT SINGLE LINES OF SIGHT
C
         IF (NUMTHE.GT.1) THEN
c
            if (iseldef.ne.0) then

               ngs = 2

               call calc_expt(iseldef,touts,tvals,maxthe,numthe,
     >                   themin,themax,maxngs,ngs,datatitle)

               ELABS(1) = BLABS
               ELABS(2) = 'EXPT'//DATATITLE
c
               CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,ngs,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)
c
            else
c
               CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,1,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,BLABS,REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)
c
            endif
c
         ELSE
           WRITE(iplot,*) TOUTS(1),TVALS(1,1)
         ENDIF
C
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
c
c
c-----------------------------------------------------------------------
c
c     DIIID bolometry plots
c
c
c-----------------------------------------------------------------------
c
c     Hydrogen signals - HPOWLS
c
      elseif (iref.eq.311) then  
c
         call rzero(plastmp,maxnks*maxnrs) 
c
         do ir = 1,nrs
c
            do ik = 1, nks(ir)
c
               do iz = 0,1
c
                  plastmp(ik,ir) = plastmp(ik,ir) 
     >                   + 1.0e-6* hpowls(ik,ir,iz)
c
c     >                           * kareas(ik,ir) 
c
               end do
c
            end do
c
         end do   
c
         call  out_bolo_diiid(plastmp,infile,iseld,graph,iopt,iflag,
     >                        job,title,table,avs,navs,iplot,nplots)
c
c
c-----------------------------------------------------------------------
c
c     Impurity - POWLS
c
      elseif (iref.eq.313) then  
c
         call rzero(plastmp,maxnks*maxnrs) 
c
c        Calculate 2D sum of impurity power loss 
c
         do ir = 1,nrs
c
            do ik = 1, nks(ir)
c
               do iz = 0,nizs
c
                  plastmp(ik,ir) = plastmp(ik,ir) 
     >                           + 1.0e-6 * powls(ik,ir,iz) * absfac
c
c     >                           * kareas(ik,ir) 
c
               end do
c
            end do
c
         end do   
c
c        Plot power loss
c
         call  out_bolo_diiid(plastmp,infile,iseld,graph,iopt,iflag,
     >                        job,title,table,avs,navs,iplot,nplots)
c
c
c-----------------------------------------------------------------------
c
c     Total - HPOWLS + POWLS * absfac
c
      elseif (iref.eq.315) then  
c
         call rzero(plastmp,maxnks*maxnrs) 
c
         do ir = 1,nrs
c
            do ik = 1, nks(ir)
c
c              Impurity component 
c
               do iz = 0,nizs
c
                  plastmp(ik,ir) = plastmp(ik,ir) 
     >                           + 1.0e-6 * powls(ik,ir,iz) * absfac
c
c     >                           * kareas(ik,ir) 
c
               end do
c
c              Hydrogen component
c
               do iz = 0,1
c
                  plastmp(ik,ir) = plastmp(ik,ir) 
     >                           + 1.0e-6 * hpowls(ik,ir,iz)
c
c     >                           * kareas(ik,ir) 
c
               end do  
c
            end do
c
         end do   
c
c        Plot power loss
c
         call  out_bolo_diiid(plastmp,infile,iseld,graph,iopt,iflag,
     >                        job,title,table,avs,navs,iplot,nplots)
c
c
c
c-----------------------------------------------------------------------
c
c      Generalized Contour plot
c
c
       elseif (iref.eq.321) then
c
c         Set IOPT=2 to get units of /m3/s/sr for emission plots
c
          call plot_contour(iselect,istate,
     >                  iexpt,optval,
     >                  iopt,job,title,table,nplots,
     >                  iplot,nizs,
     >                  minscale,maxscale,
     >                  ierr) 
c
c         Restore colours after any changes in contour routine
c
          call setup_col(n_cols,col_opt)
c
c
c-----------------------------------------------------------------------
c
c      Generalized Along Ring plot
c
c
       elseif (iref.eq.326) then
c
          call plot_ring(iselect,istate,
     >                  iexpt,minfrac,maxfrac,
     >                  axis_type,plot_type,
     >                  iopt,job,title,table,nplots,
     >                  iplot,nizs,
     >                  ierr) 
c
c         Restore colours after any changes in along ring routine
c
          call setup_col(n_cols,col_opt)
c
c
c-----------------------------------------------------------------------
c
c      Generalized LOS plot
c
c
       elseif (iref.eq.331) then
c
          call plot_los(iselect,istate,npts,navg,
     >                  iexpt,iaxis,iavg,ifact,optval,graph,
     >                  iopt,job,title,table,avs,navs,nplots,
     >                  iplot,nizs,ierr) 
c
c
c-----------------------------------------------------------------------
c
c      Generalized 3D LOS-integration plot
c
c
       elseif (iref.eq.333) then
c
          call plot_3Dlos(iselect,istate,npts,navg,
     >                  iexpt,iaxis,iavg,ifact,optval,
     >                  iopt,job,title,table,avs,navs,nplots,
     >                  iplot,nizs,ierr) 
c
c         Restore colours after any changes in contour routine
c
          call setup_col(n_cols,col_opt)
c
c
c-----------------------------------------------------------------------
c
c      Generalized 3D LOS plot - camera images
c
c
       elseif (iref.eq.335) then
c
          call plot_3Dimage(iselect,istate,
     >                  iaxis,minsteps,stepsize,optval,
     >                  iopt,job,title,table,nplots,
     >                  iplot,nizs,ierr) 
c
c         Restore colours after any changes in contour routine
c
          call setup_col(n_cols,col_opt)

c
c-----------------------------------------------------------------------
c
c      Cross-section plots across 2D contour plots 
c

      elseif (iref.eq.341) then
c
c         Call routine to plot cross-section         
c
          call plot_xsection(r1p,z1p,r2p,z2p,npts,iselect,istate,iexpt,
     >                       iavg,iopt,
     >                       nizs,job,title,table,avs,navs,iplot,nplots)
c
c
c--------------------------------------------------------------------
c
c     Plot of LINE PROFILE Data
c

      ELSEIF (IREF.EQ.345) THEN
c
c       Plots of the code simulated Doppler shifted line profile 
c       with or without experimental data.  
c
c       This plot reads an optional input line which defines the 
c       experimental data set to be plotted and any modifications
c       that should be made to this data - like axis shifts or
c       other scaling effects. It also includes a switch for 
c       relative or absolute wavelength scaling and allows for a 
c       base calibration wavelength for the experimental data to
c       be specified. 
c
c '000 LP ADDITIONAL DATA'   IEXPT  LP_PLOT_TYPE  LP_PLOT_AVG  ... 
c       EXPT_AXIS_OFFSET  LP_AXIS_OFFSET  AXIS_CALIBRATION 
c
c       IEXPT - integer - index for experimental data set
c       (LP_,EXPT_)AXIS_OFFSET - real - amount to shift data axis (lp or expt)
c       AXIS_CALIBRATION - real - calibration wavelength for V=0 profile (lp and expt)
c       PLOT_TYPE - integer - selector for absolute wavelength scale or relative
c                             to the expected peak location. 
c                 = 0 = absolute
c                 = 1 = relative
c
c       PLOT_AVG - integer - option to turn on averaging for the calculated data
c                            so that the number of points and their location 
c                            matches the experimental data. 
c
c
c       IF IOPT is less than zero then the experimental dataset 
c       to be plotted is ABS(IOPT) 
c      
c       This option can still be used if the additional input line is not present.  
c       In which case default values are loaded in the following. 
c
        if (iopt.lt.0) then
           iexpt = abs(iopt)
        else
           iexpt = 0
        endif 
c
c       Set to no pre-loaded experimental data to start
c
        expt_data_available = 0
        lp_iexpt = 0
c
c       Averaging off
c
        lp_plot_avg = 0  
c 
c       Set to absolute scale as opposed to relative to unshifted peak
c
        lp_plot_type = 0
c
c       Axis calibration and offset to 0.0
c
        expt_axis_offset  = 0.0
        lp_axis_offset    = 0.0
        axis_calibration = 0.0
c
c       Check for an additional line of input data for the plot 
c        
        call rdg_lp_plotdata(graph,lp_iexpt,lp_plot_type,lp_plot_avg,
     >                    expt_axis_offset,
     >                    lp_axis_offset,axis_calibration,ierr)
c
c
c       Set axis_calibration - if not set in input then set it to
c       the lp_wave value returned by ADAS.
c
        if (axis_calibration.eq.0.0) then 
           axis_calibration = lp_wave 
        endif
c
c       Write NVIEW 
c
        write(NVIEW,'(a,f8.2)') 'BASE WAVELENGTH =',
     >                             axis_calibration
c
c       Load and modify the experimental data into the common block 
c
        if (lp_iexpt.ne.0) then 
c
c           Load experimental data
c
            call load_expt_data(expt_dataunit,lp_iexpt,
     >                 expt_data_axis,expt_data_axis_type,
     >                 expt_data_values,
     >                 expt_data_maxcols,expt_data_maxdatx,
     >                 expt_data_num,expt_data_ncols,
     >                 expt_data_title)
c     
            if (expt_data_num.le.0) then
c     
               write(6,*) 'ERROR IN EXPERIMENTAL DATA: CAN NOT'//
     >                    ' LOAD DATASET # ',
     >                 lp_iexpt, ' - NO ELEMENTS FOUND'
               write(0,*) 'ERROR IN EXPERIMENTAL DATA: CAN NOT'//
     >                    ' LOAD DATASET # ',
     >                 lp_iexpt, ' - NO ELEMENTS FOUND'
c      
c              In the case of an error - set the actual IEXPT to zero
c
               lp_iexpt = 0
c
c           Apply adjustments to experimental data axis
c
            else
c
               do in = 1,expt_data_num
c
c                 Apply expt_axis_offset to data 
c
                  expt_data_axis(in) = expt_data_axis(in)
     >                               + expt_axis_offset
c
c                 Remove unshifted peak location from axis leaving
c                 only the profile relative to the peak. 
c 
c                 The experimental data is exected to be in the datafile
c                 as a function of the absolute wavelength. 
c
                  If (lp_plot_type.eq.1) then 
c
                     expt_data_axis(in) = expt_data_axis(in)
     >                                  - axis_calibration    
c
                  endif
c
               end do

c
c              Set the flag indicating that the experimental data is available.
c
               iexpt = lp_iexpt
               expt_data_available =1 
c
            endif
c
        endif 

c
c       Set plot options and load the simulated results for plotting. 
c
c
c       Select normalized profiles no matter what value of IOPT was
c       input. 
c
        if (lp_plot_avg.eq.0) then
           iopt = 2 
        elseif (lp_plot_avg.eq.1) then 
           iopt = 10
        endif 
c
        ngs=2 
c
        ELABS(1) = 'RAW RAW '
        ELABS(2) = 'INSTINST'
c
        ELABS(3) = 'EXPTEXPT    '
c
        XLAB = '   DLAMBDA  (A) '
        YLAB = '   NORMALIZED LINE PROFILE'
        REF  = 'DOPPLER SHIFTED EMISSION PROFILE'
        anly = ' '
        plane= ' '
c
        WRITE (IPLOT,9012) NPLOTS,REF
c
        CALL RZERO (TVALS, MAXTHE*MAXNGS)
c
        maxmod = 0.0  
c
c
c       Plot average option 0 - no averaging to match experiment
c
        if (lp_plot_avg.eq.0) then 
c
           numthe = (max_lp_bins-1) - (-max_lp_bins+1) + 1
c
c          Assign simulation data to bins
c

           do in = -max_lp_bins+1,max_lp_bins-1        
c
              ip = in - (-max_lp_bins+1) + 1
c
              touts(ip) = in * lp_bin_width + lp_axis_offset
c
              if (lp_plot_type.eq.0) then 
c
                 touts(ip) = touts(ip) + axis_calibration
c
              endif
c
              twids(ip) = lp_bin_width
              tvals(ip,1) = line_profile(in) 
              tvals(ip,2) = modified_line_profile(in) 
c
              maxmod = max(maxmod,tvals(ip,2)) 
c
              write(6,'(a,2(1x,i5),3(1x,g12.3))') 'LP:',ip,numthe,
     >           touts(ip),tvals(ip,1),tvals(ip,2)

           end do
c
c
c       Calculate averaging to match the experimental data if that option is set
c
c
        elseif (lp_plot_avg.eq.1.and.expt_data_available.eq.1) then 
c
c          Need to bin and average the SIMULATED LP data to match the EXPT data.
c       
c          Take the span of the expt data and add 2 cells at either end. 
c
c          Take the size of the outermost bins at each end and replicate these
c          up and down. 
c
c          Scan through the LP arrays and create totals for the bins defined. 
c
c          Assume the bin boundaries are half-way between the bin centers. 
c
c
           ncell_add = 3
c
           low_bin_width = abs(expt_data_axis(2)
     >                   - expt_data_axis(1))
c               
           high_bin_width = abs(expt_data_axis(expt_data_num)
     >                    - expt_data_axis(expt_data_num-1))
c
           numthe = expt_data_num + 2 * ncell_add
c
           do ip = 1,numthe
c
              if (ip.le.ncell_add) then 

                 touts(ip) = expt_data_axis(1) 
     >                     - (ncell_add-ip+1) * low_bin_width
                 twids(ip) = low_bin_width   


              elseif (ip.ge.numthe-ncell_add+1) then 

                 touts(ip) = expt_data_axis(expt_data_num) 
     >                     + (ip-(numthe-ncell_add)) * high_bin_width
                 twids(ip) = high_bin_width   

              else
c
                 touts(ip) = expt_data_axis(ip-ncell_add)
c
                 if (ip.eq.ncell_add+1) then 

                    twids(ip) = low_bin_width/2.0 + 
     >                          abs(expt_data_axis(ip-ncell_add+1)
     >                             -expt_data_axis(ip-ncell_add)) /2.0

                 elseif (ip.eq.numthe-ncell_add) then 

                    twids(ip) = high_bin_width/2.0 + 
     >                          abs(expt_data_axis(ip-ncell_add)
     >                             -expt_data_axis(ip-ncell_add-1)) /2.0

                 else
c
                    twids(ip) = abs(expt_data_axis(ip-ncell_add+1)
     >                             -expt_data_axis(ip-ncell_add)) /2.0 +
     >                          abs(expt_data_axis(ip-ncell_add)
     >                             -expt_data_axis(ip-ncell_add-1)) /2.0
c
                 endif
c
              endif     
c
              write(6,'(a,2i5,5(1x,f18.6))') 'LPT:',ip,numthe,touts(ip),
     >             twids(ip)


           end do
c
c          Loop through the simulation LP and average the data for each bin 
c
           do ip = 1,numthe
c
              tvals(ip,1) = 0.0
              tvals(ip,2) = 0.0
              lp_cnt = 0.0
c	      
c             Loop through LP data
c     	      
              do in = -max_lp_bins+1,max_lp_bins-1        
c	      
c                Calculate the axis coordinate of the LP data
c	      
                 lp_out = in * lp_bin_width + lp_axis_offset  
c	      
                 if (lp_plot_type.eq.0) then 
                    lp_out = lp_out + axis_calibration
                 endif  
c	      
c                This assumes even spacing of the experimental data - which 
c                is true at the present time. 
c	      
                 if (lp_out.gt.touts(ip)-twids(ip)/2.0.and.
     >               lp_out.le.touts(ip)+twids(ip)/2.0) then 
                     tvals(ip,1) = tvals(ip,1) + line_profile(in)
                     tvals(ip,2) = tvals(ip,2) 
     >                             + modified_line_profile(in)
                     lp_cnt = lp_cnt + 1.0
c 
                     write(6,'(a,2i5,5(1x,f18.6))') 'LPA:',ip,in,
     >                     lp_out,touts(ip),twids(ip),
     >               tvals(ip,1),tvals(ip,2)
                 endif  
c	      
              end do
c	      
c             Calculate average 
c	      
              tvals(ip,1) = tvals(ip,1) / lp_cnt 
              tvals(ip,2) = tvals(ip,2) / lp_cnt 
c
              maxmod = max(maxmod,tvals(ip,2)) 
c
              write(6,'(a,2(1x,i5),3(1x,g12.3))') 'LP:',ip,numthe,
     >           touts(ip),tvals(ip,1),tvals(ip,2)
c
           end do
c
       endif  
c
c
c       Determine THEMIN - the modified line profile is always
c       wider than the regular one.
c
c       The instrument profile contains extremely small values 
c       which are not worth plotting - so the plotting range
c       is truncated to a fraction of the peak height. 
c
        do ip = 1,numthe        
c
           if (tvals(ip,2).gt.1.0e-5*maxmod) then 
              themin = touts(ip)
              exit
           endif
c
        end do   
c
c       Determine THEMAX  
c
        do ip = numthe,1,-1
c
           if (tvals(ip,2).gt.1.0e-5*maxmod) then 
              themax = touts(ip)
              exit
           endif
c
        end do   
c
c       Generate plot
c
        CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,ngs,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,iexpt)
c

c
c
c--------------------------------------------------------------------
c
c     Plot of LINE PROFILE Data
c

      ELSEIF (IREF.EQ.346) THEN
c
c        IF IOPT is less than zero then call with the experimental
c        data set specified - the negative value will result in a 
c        plot. 
c
        if (iopt.lt.0) then
           iexpt = iopt
        else
           iexpt = 0
        endif 
c
c       Select normalized profiles no matter what value of IOPT was
c       input. 
c
        iopt = 2 
        ngs=2 
c
        ELABS(1) = 'RAW RAW '
        ELABS(2) = 'INSTINST'
c
        ELABS(3) = 'EXPTEXPT    '
c
        XLAB = '   DLAMBDA  (A) '
        YLAB = '   NORMALIZED LINE PROFILE'
        REF  = 'DOPPLER SHIFTED EMISSION PROFILE'
        anly = ' '
        plane= ' '
c
        write(NVIEW,'(a,f8.2)') 'BASE WAVELENGTH =',lp_wave
c
        WRITE (IPLOT,9012) NPLOTS,REF
c
        CALL RZERO (TVALS, MAXTHE*MAXNGS)
c
        numthe = (max_lp_bins-1) - (-max_lp_bins+1) + 1
c
        maxmod = 0.0  
c
        do in = -max_lp_bins+1,max_lp_bins-1        
c
           ip = in - (-max_lp_bins+1) + 1
c
           touts(ip) = in * lp_bin_width + lp_wave
           twids(ip) = lp_bin_width
           tvals(ip,1) = line_profile(in) 
           tvals(ip,2) = modified_line_profile(in) 
c
           maxmod = max(maxmod,tvals(ip,2)) 
c
           write(6,'(a,i5,3(1x,g12.3))') 'LP:',ip,
     >           touts(ip),tvals(ip,1),tvals(ip,2)

        end do
c
c       Determine THEMIN - the modified line profile is always
c       wider than the regular one.
c
c       The instrument profile contains extremely small values 
c       which are not worth plotting - so the plotting range
c       is truncated to a fraction of the peak height. 
c
        do ip = 1,numthe        
c
           if (tvals(ip,2).gt.1.0e-5*maxmod) then 
              themin = touts(ip)
              exit
           endif
c
        end do   
c
c       Determine THEMAX  
c
        do ip = numthe,1,-1
c
           if (tvals(ip,2).gt.1.0e-5*maxmod) then 
              themax = touts(ip)
              exit
           endif
c
        end do   
c
c        write(0,*) 'M:',themin,themax
c
c       Generate plot
c
        CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,ngs,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,iexpt)
c

c
c-----------------------------------------------------------------------
c
c
c
c
c     RCP - reciprocating probe comparison plots
c
      elseif (iref.eq.351) then
c
c        The input data for the RCP probe contains the following
c        information in the following format:
c
c '000 OSM PROBE' R1P  Z1P  R2P  Z2P  INT_TYPE
c '000 RCP PROBE' EXP_DS EXP_DATAOPT EXP_VCALCOPT EXP_TCALCOPT  EXP_PARAM
c                        EXP_OFFSET
c
c
c     Tasks required:
c
c     1) Extract the ne,Te,Ti and Pressure from the grid along the
c        probe length. Map these as a function of the mid-plane
c        coordinate of the corresponding rings. At present, only
c        intersections with the main SOL are calculated - the
c        core and private plasma rings are not calculated.
c
c     2) Read in the experimental probe data - Jsat (1,2), Te and
c        the outer mid-plane coordinates for each set of data.
c        Using the exp_calcopt specified convert the data to
c        ne, Te, Ti and Pressure. Apply X-axis offset.
c
c     3) Plot each of the corresponding sets of data on the same
c        graph. This will require some work in specifying the bounds
c        of the plot region as well as the fact that the two plots
c        will not share a common set of X-axis values.
c
c        Extract OSM data
c
c        Use exp_dataopt to turn on and off averaging 
c
         if (exp_dataopt.lt.0) then 
            exp_av = .true.
            exp_dataopt = -exp_dataopt
         else
            exp_av = .false.
         endif
c 

         WRITE (REF,'(''OSM PROBE PLOT:'',4(1x,f6.3))')
     >                                r1p,z1p,r2p,z2p
c
         call osmprobe(lvals,louts,osmvals,osmplots,
     >                 r1p,z1p,r2p,z2p,int_type,crmb,qtim)
c
c        Extract RCP data - if dataset specified
c
         if (exp_ds.gt.0) then
            call rcpprobe(tvals,touts,rcpvals,rcpplots,
     >                 exp_ds,exp_offset,exp_dataopt,
     >                 exp_vcalcopt,exp_tcalcopt,exp_param,
     >                 lvals,louts,osmvals,osmplots,
     >                 int_type,rizb,crmb,datatitle,r1p,z1p,r2p,z2p)

            WRITE (REF,'(''351: EXPT/OSM PROBE PLOT: OFFSET ='',f9.5)')
     >                                exp_offset
c
         endif
c
c        Eliminate all RCP data not within the specified plotting
c        window for R,Z data
c
         axmin= -HI
         axmax = HI
c
         if (int_type.eq.2) then
c
            axmin = min(r1p,r2p)
            axmax = max(r1p,r2p)
c
         elseif (int_type.eq.3) then
c
            axmin = min(z1p,z2p)
            axmax = max(z1p,z2p)
c
c        Quick fix for PSI plots looking at SOL
c
         elseif (int_type.eq.6) then   
            axmin = 0.98
            axmax = 1.5
         endif
c
c        Remove RCP data not in plotting range
c
         if ((axmin.ne.-HI.or.axmax.ne.HI).and.
     >        rcpvals.gt.0) then
c
            icnt = 0
c
            do in = 1,rcpvals
c
c              check axis coordinate and remove if necessary
c
               if (touts(in).lt.axmin.or.touts(in).gt.axmax) then
c
                  icnt = icnt + 1
c
c              Copy to correct position if icnt non-zero
c
               elseif (icnt.ne.0) then

                   touts(in-icnt) = touts(in)
c
                   do ip = 1,rcpplots
c
                      tvals(in-icnt,ip) = tvals(in,ip)
c
                   end do
c
               endif
c
            end do
c
c           Fix rcpvals which defines number of data points
c
            rcpvals = rcpvals - icnt
c
            write (6,*) ' Modified EXPT probe results:',rcpvals
c
            do in = 1,rcpvals
c
               write(6,'(a,i4,6(1x,g15.6))') 'RCP:',in,touts(in),
     >           tvals(in,1),tvals(in,2),tvals(in,3),tvals(in,4),
     >           tvals(in,5)
c
            end do
c
c           If averaging is ON - loop through remaining RCP values and 
c           calculate average values for each specific Z coordinate.
c
c           Also - remove all datapoints with any value greater than 
c           2 x the calculated average. 
c
            if (exp_av) then 
c
               call izero(av_cnt,maxdatx)
               call izero(av_ref,maxdatx)
               call rzero(av_vals,maxdatx*maxngs)
               call rzero(av_outs,maxdatx)
c
               num_av = 0
               tmp_av = 0 
c
               do in = 1,rcpvals
c
                  if (num_av.ge.1) then 

                     do id = 1,tmp_av
c
c                        write(6,'(a,4i5,6(1x,g16.8))')
c     >                         'NUM_AV:',in,id,tmp_av,num_av,
c     >                          av_outs(id),touts(in)
c
                        if (av_outs(id).eq.touts(in)) then 
c
                           av_ref(in) = id
c
                           av_cnt(id) = av_cnt(id) + 1.0
c
                           do ip = 1,rcpplots
                              av_vals(id,ip) = av_vals(id,ip)
     >                                    + tvals(in,ip)
                           end do

                           exit

                        elseif (id.eq.tmp_av) then 

                           num_av = num_av + 1

                           av_ref(in) = num_av  
                           av_cnt(num_av) = 1.0  
                           av_outs(num_av) = touts(in)

                           do ip = 1,rcpplots
                              av_vals(num_av,ip) = tvals(in,ip)
                           end do

                          exit 

                        endif                  

                     end do 

                  else

                     num_av = num_av + 1
		  
                     av_ref(in) = num_av  
                     av_cnt(num_av) = 1.0
                     av_outs(num_av) = touts(in)
		  
                     do ip = 1,rcpplots
                        av_vals(num_av,ip) = tvals(in,ip)
                     end do

                  endif
c
                  tmp_av = num_av
c
               end do 
c  
c              Calculate averages for each distinct Axis-value  
c
               do id = 1,num_av
c        
                  do ip = 1,rcpplots 
c
                     if (av_cnt(id).gt.0) then 
                        av_vals(id,ip) = av_vals(id,ip)/av_cnt(id)
                     else  
                        av_vals(id,ip) = 0.0
                     endif 
c
                  end do 
c
               end do 
c
c              Write out the average data:
c
               write (6,*) ' Averaged EXPT probe results:',num_av
c
               do in = 1,num_av
c
                  write(6,'(a,i4,6(1x,g15.6))') 'RCP:',in,av_outs(in),
     >         av_vals(in,1),av_vals(in,2),av_vals(in,3),av_vals(in,4),
     >         av_vals(in,5)
c
               end do
c
c              Remove experimental data that exceeds twice the average at the
c              specific Z-position. 
c
               icnt = 0
c
               do in = 1,rcpvals
c
c                 check axis coordinate and remove if necessary
c
                  exclude_pt = .false.
                  id = av_ref(in)
c
                  do ip = 1,rcpplots
c
                     if (abs(tvals(in,ip)).gt.
     >                   abs(2.0*av_vals(id,ip))) then
                        exclude_pt = .true.
                        exit
                     endif
c
                  end do
c  
c                 Exclude data point
c
                  if (exclude_pt) then  
c
                     icnt = icnt + 1
c
c                 Copy to correct position if icnt non-zero
c
                  elseif (icnt.ne.0) then
c
                      touts(in-icnt) = touts(in)
c
                      do ip = 1,rcpplots
c
                         tvals(in-icnt,ip) = tvals(in,ip)
c
                      end do
c
                  endif
c
               end do
c
c              Fix rcpvals which defines number of data points
c
               rcpvals = rcpvals - icnt
c
               write (6,*) ' Average Modified EXPT probe results:',
     >                 rcpvals
c
               do in = 1,rcpvals
c
                  write(6,'(a,i4,6(1x,g15.6))') 'RCPA:',in,touts(in),
     >           tvals(in,1),tvals(in,2),tvals(in,3),tvals(in,4),
     >           tvals(in,5)
c
               end do
c
            endif  
c
         endif

c
c        Print out a data file for Kevin
c
         if (cgrprint.eq.1) then

            len = lenstr(title)
c
            write(56,*) title(1:len)
c
            in = 1
c
            do ir = irsep,irwall-1
c
               write(56,'(i4,14(1x,g13.6))')
     >           ir,middist(ir,2),
c
c                 Jsatt, Tet, Tit, Net
c
     >            knds(idds(ir,2)) *  (ech * (9.79E+03 *
     >            SQRT(0.5*(kteds(idds(ir,2))+KTIds(idds(ir,2)))
     >                 *(1.0+RIZB)/CRMB))),
     >            kteds(idds(ir,2)),ktids(idds(ir,2)),
     >            knds(idds(ir,2)),
c
c                 Jsatu, Teu, Tiu, Neu
c
     >            lvals(in,1) *  (ech * (9.79E+03 *
     >            SQRT(0.5*(lvals(in,2)+lvals(in,3))
     >                 *(1.0+RIZB)/CRMB))),
     >            lvals(in,2),lvals(in,3),lvals(in,1),
c
c                 Xperps - outer-e, outer-i, total e, total combined
c
     >            ochiperpe(ir),ochiperpi(ir),chiperpe(ir),
     >            xperpt(ir)
c
               in = in +1
c
            end do
c
         endif
c
c
c        Generate plots
c
c        Place the contents of corresponding plots on a common plot
c        i.e. plots of ne, Te, Ti and Pressure.
c
c
c
         ELABS(1) = 'OSM OSM '
         Len = lenstr(datatitle)
c         ELABS(2) = 'EXPTEXPT'//datatitle(1:len)
         ELABS(2) = 'EXPTEXPT'

         if (exp_av) ELABS(3) = 'EXAVEXAV'
c
         PLTLABS(1) =  'Ne         '
         PLTLABS(2) =  'Te         '
         PLTLABS(3) =  'Ti         '
         PLTLABS(4) =  'Pressure   '
c
         if (int_type.eq.1) then 
            XLAB = 'R (M) OUTER MIDPLANE'
         elseif(int_type.eq.4.or.int_type.eq.5) then
            XLAB = 'RHO (M)'
         elseif (int_type.eq.2) then
            XLAB = 'R (M) ALONG PROBE LINE'
         elseif (int_type.eq.3) then
            XLAB = 'Z (M) ALONG PROBE LINE'
         elseif (int_type.eq.6) then
            XLAB = 'PSIN' 
         endif
c
         YLAB   = 'Probe Results'
c
         NPLOTS = NPLOTS + 1
         WRITE (IPLOT,9012) NPLOTS,REF
c
         CALL rzero (mvals, maxdatx*maxngs*maxplts)
c
         sctype = iopt
         ngrm  = 4
         nplts = 4
c
c        Set to plot only OSM or OSM + RCP data
c
         if (exp_ds.gt.0) then
            ngs = 2
            if (exp_av) ngs = 3
         else
            ngs = 1
         endif
c
c        Loop through the data for the 4 plots
c
c        Loading one set of information for each plot from each
c        probe result - OSM and RCP
c
         do ip = 1,nplts
c
c           Load OSM first
c
            pnks(ip,1) = osmvals
c
            do ik = 1,osmvals
c
               mouts(ik,ip,1) = louts(ik)
c
               mvals(ik,ip,1) = lvals(ik,ip)
c
            end do
c
            if (ngs.gt.1) then
c
c              Load RCP next
c
               pnks(ip,2) = rcpvals
c
               do ik = 1,rcpvals
c
                  mouts(ik,ip,2) = touts(ik)
c
                  mvals(ik,ip,2) = tvals(ik,ip)
c
               end do
c
               if (ngs.gt.2) then  

                  pnks(ip,3) = num_av
c
                  do ik = 1,num_av
c
                     mouts(ik,ip,3) = av_outs(ik)
c
                     mvals(ik,ip,3) = av_vals(ik,ip)
c
                  end do
c
               endif 
c
            endif

c
         end do
c
c        Set drawing style for the different sets of data
c
         drawtype(1) = 1
         drawtype(2) = 2
         if (exp_av) drawtype(3) = 1
C
c slmod begin - new
         Len = lenstr(datatitle)
         do in = 1,min(len/20+1,8)
            extra_comments(in) = 
     >             datatitle((in-1)*20+1:min(in*20,len))         
         end do  
c
         in = min(len/20+1,8)
c
         CALL SetPlotComments(iref,job,extra_comments,in,0.0)
c slmod end
c
c        Set pltfact to a value that will result in a slight extension 
c        of the X-axis to make the plotting look neater.
c
         pltfact = 1.02
c
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
               mdrawtype(ip,ig) = drawtype(ig)
            end do  
            pngs(ip) = ngs
         end do
c
         CALL DRAWM (MOUTS,MWIDS,MVALS,Maxdatx,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,ngs,
     >              mdrawtype,1)
c
c
c-----------------------------------------------------------------------
c
      elseif (iref.eq.353) then
c
c        The input data for the RCP probe contains the following
c        information in the following format:
c
c '000 OSM PROBE' R1P  Z1P  R2P  Z2P  INT_TYPE
c '000 RCP PROBE' EXP_DS EXP_DATAOPT EXP_VCALCOPT EXP_TCALCOPT  EXP_PARAM
c                        EXP_OFFSETS (1 to 4)
c
c
c     Tasks required:
c
c     1) Extract the ne,Te,Ti and Pressure from the grid along the
c        probe length. Map these as a function of the mid-plane
c        coordinate of the corresponding rings. At present, only
c        intersections with the main SOL are calculated - the
c        core and private plasma rings are not calculated.
c
c     2) Read in the experimental probe data - Jsat (1,2), Te and
c        the outer mid-plane coordinates for each set of data.
c        Using the exp_calcopt specified convert the data to
c        ne, Te, Ti and Pressure. Apply X-axis offset.
c
c     3) Plot each of the corresponding sets of data on the same
c        graph. This will require some work in specifying the bounds
c        of the plot region as well as the fact that the two plots
c        will not share a common set of X-axis values.
c
c slmod begin
c...     Check if parallel plasma flow should be over-written:
         optflow = 0
         READ(5,'(A80)') graph1
         IF (graph1(8:11).EQ.'Flow'.OR.graph1(8:11).EQ.'FLOW'.OR.
     .       graph1(8:11).EQ.'flow') THEN
           READ(graph1,*) cdum1,optflow
           CALL CalcFlow(optflow)
           restoresolution = .TRUE.
         ELSE
           BACKSPACE 5
         ENDIF
c slmod end
c        Extract OSM data
c
         CALL rzero (mvals, maxdatx*maxngs*maxplts)
c
         REF = ' '
c
         call osmprobe(lvals,louts,osmvals,osmplots,
     >                 r1p,z1p,r2p,z2p,int_type,crmb,qtim)
c
c        Load probe data into plotting arrays for all 4 offsets
c
         do in = 1,4

            do ip2 = 1,4
c
c           Load OSM first
c
               ip = (in-1)*4+ip2
c
               pnks(ip,1) = osmvals
c
               do ik = 1,osmvals
c
                  mouts(ik,ip,1) = louts(ik)
c
                  mvals(ik,ip,1) = lvals(ik,ip2)
c
               end do
c
            end do  
c
         end do
c slmod begin
         IF (.FALSE.) THEN
c...       Plot the outer target data versus the scanning probe data
c          instead:
           do in = 1,4
             ip = (in-1)*4
             pnks(ip,1) = osmvals
             do ik = 1,osmvals
               id = idds(ik+irsep-1,1)
               mvals(ik,ip+1,1) = knds(id)
               mvals(ik,ip+2,1) = kteds(id)
               mvals(ik,ip+3,1) = ktids(id)
               mvals(ik,ip+4,1) = CalcPressure(knds(id),kteds(id),
     .                                         ktids(id),kvds(id)) * ECH
             end do
           end do
         ELSEIF (exp_vcalcopt.EQ.4) THEN
c...       Replace Ti data on the 3rd plot with the parallel flow velocity, and
c          calculate pressure without velocity:
           do in = 1,4
             ip = (in-1)*4
             do ik = 1,osmvals
               mvals(ik,ip+3,1) = lvals(ik,5) / GetCs(mvals(ik,ip+2,1),
     .                                                mvals(ik,ip+3,1))
               mvals(ik,ip+4,1) = CalcPressure
     .                              (lvals(ik,1),lvals(ik,2),
     .                               lvals(ik,3),lvals(ik,5)) * ECH
c               mvals(ik,ip+4,1) = CalcPressure
c     .                              (lvals(ik,1),lvals(ik,2),
c     .                               lvals(ik,3),0.0) * ECH
             end do
           end do
         ENDIF
c slmod end

c
c
c        SET Axis ranges
c
         axmin= -HI
         axmax = HI
c
         if (int_type.eq.2) then
c
            axmin = min(r1p,r2p)
            axmax = max(r1p,r2p)
c
         elseif (int_type.eq.3) then
c
            axmin = min(z1p,z2p)
            axmax = max(z1p,z2p)
c
         endif

c
c        Extract RCP data - if dataset specified
c
c
         if (exp_ds.gt.0) then
c
            do in1 = 1,4
c
               exp_offset = exp_offsets(in1) 
               write (6,*) 'EXP_OFFSET:',exp_offset 
c
               call rcpprobe(tvals,touts,rcpvals,rcpplots,
     >                 exp_ds,exp_offset,exp_dataopt,
     >                 exp_vcalcopt,exp_tcalcopt,exp_param,
     >                 lvals,louts,osmvals,osmplots,
     >                 int_type,rizb,crmb,datatitle,r1p,z1p,r2p,z2p)

               WRITE (REF,'(''EXPT/OSM PROBE PLOT:'')')

c
c              Remove RCP data not in plotting range
c
               if (axmin.ne.-HI.or.axmax.ne.HI) then
c
                  icnt = 0
c
                  do in = 1,rcpvals
c
c                   check axis coordinate and remove if necessary
c
                    if (touts(in).lt.axmin.or.touts(in).gt.axmax) then
c
                       icnt = icnt + 1
c
c                      Copy to correct position if icnt non-zero
c
                    elseif (icnt.ne.0) then
     
                       touts(in-icnt) = touts(in)
c
                       do ip = 1,rcpplots
c
                           tvals(in-icnt,ip) = tvals(in,ip)
c
                       end do
c
                    endif
c
                 end do
c
c                Fix rcpvals which defines number of data points
c
                 rcpvals = rcpvals - icnt
c

               endif

c
c              Store RCP data for each offset into plotting arrays
c

               do ip2 = 1,4
c
                  ip = (in1-1)*4 + ip2  
c
                  pnks(ip,2) = rcpvals
c
                  do ik = 1,rcpvals
c
                     mouts(ik,ip,2) = touts(ik)
c
                     mvals(ik,ip,2) = tvals(ik,ip2)
c
c slmod begin
                     IF (exp_vcalcopt.EQ.4.AND.ip2.EQ.3) 
     .                 mvals(ik,ip,2) = tvals(ik,5) / 
     .                                  GetCs(tvals(ik,2),tvals(ik,3))
c slmod end
                  end do
c
               end do
c
            end do
c
         endif
c
c
c        Generate plots
c
c        Place the contents of corresponding plots on a common plot
c        i.e. plots of ne, Te, Ti and Pressure.
c
c
c
         ELABS(1) = 'OSM OSM '
         Len = lenstr(datatitle)
         ELABS(2) = 'EXPTEXPT'//datatitle(1:Len)
c
         do in = 1,4
c
            ip = (in-1)*4
c
c           Initialize
c
            pltlabs(ip+1) = ' '
            pltlabs(ip+2) = ' '
            pltlabs(ip+3) = ' '
            pltlabs(ip+4) = ' '
c
c           Set plot labels
c
c            write(pltlabs(ip+1),'(a,f5.3)') 'Ne: SHIFT=',
c     >                               exp_offsets(in)
c            write(pltlabs(ip+2),'(a,f5.3)') 'Te: SHIFT=',
c     >                               exp_offsets(in)
c            write(pltlabs(ip+3),'(a,f5.3)') 'Ti: SHIFT=',
c     >                               exp_offsets(in)
c            write(pltlabs(ip+4),'(a,f5.3)') 'Pr: SHIFT=',
c     >                               exp_offsets(in)

            write(pltlabs(ip+1),'(a,f7.4)') 'Ne: SHIFT=',
     >                               exp_offsets(in)
            write(pltlabs(ip+2),'(a,f7.4)') 'Te: SHIFT=',
     >                               exp_offsets(in)
            write(pltlabs(ip+3),'(a,f7.4)') 'Ti: SHIFT=',
     >                               exp_offsets(in)
            write(pltlabs(ip+4),'(a,f7.4)') 'Pr: SHIFT=',
     >                               exp_offsets(in)
c
         end do
c
c         PLTLABS(1) =  'Ne         '
c         PLTLABS(2) =  'Te         '
c         PLTLABS(3) =  'Ti         '
c         PLTLABS(4) =  'Pressure   '
c
         if (int_type.eq.1) then 
            XLAB = 'R (M) OUTER MIDPLANE'
         elseif(int_type.eq.4.or.int_type.eq.5) then
            XLAB = 'RHO (M)'
         elseif (int_type.eq.2) then
            XLAB = 'R (M) ALONG PROBE LINE'
         elseif (int_type.eq.3) then
            XLAB = 'Z (M) ALONG PROBE LINE'
         endif
c
         YLAB   = 'Probe Results'
c
         NPLOTS = NPLOTS + 1
         WRITE (IPLOT,9012) NPLOTS,REF
c
c slmod begin - new
c...     Sets some formatting options:
         slopt3 = 3
c slmod end
         sctype = iopt
         ngrm  = 16
         nplts = 16
c
c       Set to plot only OSM or OSM + RCP data
c
        if (exp_ds.gt.0) then
           ngs = 2
        else
           ngs = 1
        endif
c
         write(6,*) 'NGS:',ngs,exp_ds
c
c        Set drawing style for the different sets of data
c
c         drawtype(1) = 1
         drawtype(1) = 3
         drawtype(2) = 2
c slmod begin - new
c...     Load fluid code data from E2Dxxx arrays:
c
c        Only load if data exists
c
         IF (cre2d.gt.0) THEN
           ngs = ngs + 1
           CALL RZero(touts2,MAXTHE)
           CALL RZero(tvals2,MAXTHE*MAXNGS)
           call fluidprobe(tvals2,touts2,osmvals,osmplots,
     >                     r1p,z1p,r2p,z2p,int_type,crmb,1.0)
           do in = 1,4
              do ip2 = 1,4
                 ip = (in-1)*4+ip2
                 pnks(ip,ngs) = osmvals
                 do ik = 1,osmvals
                    mouts(ik,ip,ngs) = touts2(ik)
                    mvals(ik,ip,ngs) = tvals2(ik,ip2)
                 end do
              end do  
           end do

           ELABS   (ngs) = 'UE  UEDGE'  
           drawtype(ngs) = 1
         ENDIF

         Len = lenstr(datatitle)
         extra_comments(1) = datatitle(1:len)         

         CALL SetPlotComments(iref,job,extra_comments,1,0.0)


c         WRITE(0,*) 'GRAPH:',graph
c         STOP 'rwewer'

         IF (.TRUE.) THEN
c...       Map independent variable to rho (one the grid only):




         ENDIF
c slmod end

c
c        Set pltfact to a value that will result in a slight extension 
c        of the X-axis to make the plotting look neater.
c
         pltfact = 1.02
c
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
               mdrawtype(ip,ig) = drawtype(ig)
            end do  
            pngs(ip) = ngs
         end do
c
c slmod begin
c...     Output data to be plotted:
         WRITE(6,*)
         WRITE(6,*) 'PLOT 353 DATA:'
         DO ip = 1, nplts
           WRITE(6,*) 
           WRITE(6,'(A,I6)') ' PLOT: ',ip
           DO i2 = 1, 2
             DO i1 = 1, pnks(ip,i2)
               WRITE(6,'(I6,A,1P,2E12.4,0P)') 
     .           i2,':',mouts(i1,ip,i2),mvals(i1,ip,i2)
             ENDDO
           ENDDO
         ENDDO

c...     Check if calculation of a norm is requested:
         nrmcalculate = .FALSE.
         READ(5,'(A128)',END=55) graph6
         IF (graph6(8:11).EQ.'Norm'.OR.graph6(8:11).EQ.'NORM'.OR.
     .       graph6(8:11).EQ.'norm') THEN
           READ(graph6,*) cdum1,nrmtype,nrmi1,nrmi2,nrmr1,nrmr2,
     .                    cdum1
           nrmcalculate = .TRUE.
         ELSE
           BACKSPACE 5
         ENDIF
55       CONTINUE
         
c '000    Norm'  TYPE IND1 IND2 XRANGE1 XRANGE2 COMMENT
         
         IF (nrmcalculate) THEN

c           DO ip = 0, 0
           DO ip = 0, nplts-1, 4

             DO in = 1, 4

c...           Update the norm index:
               nrmindex = nrmindex + 1
               CALL RSet(nrmdata,MAXTHE*2,LO)

               IF (in.EQ.3) THEN
c...             Norm is based on where the upstream flow direction
c                changes:

                 DO i4 = 1, 2
                   i3 = 0                      
                   DO i1 = 1, pnks(ip+in,i4)-1
                     IF ((mvals(i1  ,ip+in,i4).GT.0.0.AND.
     .                    mvals(i1+1,ip+in,i4).LT.0.0).OR.
     .                   (mvals(i1  ,ip+in,i4).LT.0.0.AND.
     .                    mvals(i1+1,ip+in,i4).GT.0.0)) THEN
                       IF (i3.EQ.0) THEN
                         i3 = 1
                         frac = -mvals(i1,ip+in,i4) /
     .                       (mvals(i1+1,ip+in,i4) - mvals(i1,ip+in,i4))
                         nrmdata(1,i4) = mouts(i1,ip+in,i4) + frac * 
     .                       (mouts(i1+1,ip+in,i4) - mouts(i1,ip+in,i4))
                       ELSE
c...                     Expecting the data to cross zero only once, so
c                        turn off this norm otherwise:
                         nrmdata(1,i4) = LO
                       ENDIF
                     ENDIF
                   ENDDO
                 ENDDO

                 WRITE(0,*) 'NORM 3:',nrmdata(1,1),nrmdata(1,2)

                 nrmdata(1,1) = nrmdata(1,1) - mouts(1,ip+in,1)
                 nrmdata(1,2) = nrmdata(1,2) - mouts(1,ip+in,1)

               ELSE

                 i3 = 0
                 DO i1 = 1, pnks(ip+in,2)

                   IF ((mouts(i1,ip+in,2).LT.nrmr1.OR.
     .                  mouts(i1,ip+in,2).GT.nrmr2).AND.
     .                 (nrmr1.NE.0.0.OR.nrmr2.NE.0.0).AND.
     .                 in.NE.4) CYCLE

c...               Assign data in the specified range (x-axis data):            
	           i3 = i3 + 1
                   nrmdata(i3,1) = mvals(i1,ip+in,2)

                   WRITE(0,*) 'Data:',in+ip,i3,
     .                   mouts(i1,ip+in,2)
	    
c...               Secondary data linearly interpolated:
                   nrmdata(i3,2) = LO
                   DO i4 = 1, pnks(ip+in,1)-1
                     IF ((mouts(i4  ,ip+in,1).LE.mouts(i1,ip+in,2).AND. 
     .                    mouts(i4+1,ip+in,1).GE.mouts(i1,ip+in,2)).OR.
     .                   (mouts(i4  ,ip+in,1).GE.mouts(i1,ip+in,2).AND. 
     .                    mouts(i4+1,ip+in,1).LE.mouts(i1,ip+in,2))) 
     .                 THEN 
                       frac = (mouts(i1  ,ip+in,2) - mouts(i4,ip+in,1))/
     .                        (mouts(i4+1,ip+in,1) - mouts(i4,ip+in,1))
                       nrmdata(i3,2) = mvals(i4,ip+in,1) +
     .                    frac * (mvals(i4+1,ip+in,1)-mvals(i4,ip+in,1))
                       EXIT
                     ENDIF
                   ENDDO

                   IF (nrmdata(i3,2).EQ.LO) THEN
                     IF (.TRUE..OR.in+ip.EQ.1) 
     .                 WRITE(0,*) 'No interplolated:',in+ip,
     .                   mouts(i1,ip+in,2)
                   ENDIF

                 ENDDO

               ENDIF

c...           Store the number of data sets recorded:
               nrmnum(nrmindex) = i3

               WRITE(nrmcomment(nrmindex),'(A,I6)') 
     .           cdum1(1:LEN_TRIM(cdum1)),ip+in
          
c...           Generate the norm:
               CALL CalculateNorm
             ENDDO
           ENDDO

         ENDIF
c slmod end
         CALL DRAWM (MOUTS,MWIDS,MVALS,Maxdatx,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,ngs,
     >              mdrawtype,1)
c
c
c-----------------------------------------------------------------------
c
      elseif (iref.eq.355) then
c
c        The input data for the RCP probe contains the following
c        information in the following format:
c
c '000 OSM PROBE' R1P  Z1P  R2P  Z2P  INT_TYPE
c '000 RCP PROBE' EXP_DS EXP_DATAOPT EXP_VCALCOPT EXP_TCALCOPT  EXP_PARAM
c                        EXP_OFFSETS (1 to 4)
c
c
c     Tasks required:
c
c     1) Extract the ne, Te and Pressure from the grid along the
c        probe length. Map these as a function of the mid-plane
c        coordinate of the corresponding rings. At present, only
c        intersections with the main SOL are calculated - the
c        core and private plasma rings are not calculated.
c
c     2) Read in the experimental probe data - Jsat (1,2), Te and
c        the outer mid-plane coordinates for each set of data.
c        Using the exp_calcopt specified convert the data to
c        ne, Te and Pressure. Apply X-axis offset.
c
c     3) Plot each of the corresponding sets of data on the same
c        graph. This will require some work in specifying the bounds
c        of the plot region as well as the fact that the two plots
c        will not share a common set of X-axis values.
c
c        Set gray scale colours
c
         call setup_col(n_cols,1)
c
c        Extract OSM data
c
         CALL rzero (mvals, maxdatx*maxngs*maxplts)
c
         REF = ' '
c
         call osmprobe(lvals,louts,osmvals,osmplots,
     >                 r1p,z1p,r2p,z2p,int_type,crmb,qtim)

c
c        Load probe data into plotting arrays for 2 offsets
c
         do in = 1,2

            do ip2 = 1,3
c
c           Load OSM first
c
c               ip = (in-1)*3+ip2
               ip = 2*ip2 + in - 2  
c
               pnks(ip,1) = osmvals
c
               do ik = 1,osmvals
c
                  mouts(ik,ip,1) = louts(ik)
c
c                 Load pressure instead of Ti for 3rd element
c
                  if (ip2.eq.3) then 
                     mvals(ik,ip,1) = lvals(ik,ip2+1)
                  else
                     mvals(ik,ip,1) = lvals(ik,ip2)
                  endif
c
               end do
c
            end do  
c
         end do
c
c
c        SET Axis ranges
c
         axmin= -HI
         axmax = HI
c
         if (int_type.eq.2) then
c
            axmin = min(r1p,r2p)
            axmax = max(r1p,r2p)
c
         elseif (int_type.eq.3) then
c
            axmin = min(z1p,z2p)
            axmax = max(z1p,z2p)
c
         endif

c
c        Extract RCP data - if dataset specified
c
c
         if (exp_ds.gt.0) then
c
            do in1 = 1,2
c
               exp_offset = exp_offsets(in1) 
               write (6,*) 'EXP_OFFSET:',exp_offset 
c
               call rcpprobe(tvals,touts,rcpvals,rcpplots,
     >                 exp_ds,exp_offset,exp_dataopt,
     >                 exp_vcalcopt,exp_tcalcopt,exp_param,
     >                 lvals,louts,osmvals,osmplots,
     >                 int_type,rizb,crmb,datatitle,r1p,z1p,r2p,z2p)

               WRITE (REF,'(''EXPT/OSM PROBE PLOT:'')')

c
c              Remove RCP data not in plotting range
c
               if (axmin.ne.-HI.or.axmax.ne.HI) then
c
                  icnt = 0
c
                  do in = 1,rcpvals
c
c                   check axis coordinate and remove if necessary
c
                    if (touts(in).lt.axmin.or.touts(in).gt.axmax) then
c
                       icnt = icnt + 1
c
c                      Copy to correct position if icnt non-zero
c
                    elseif (icnt.ne.0) then
     
                       touts(in-icnt) = touts(in)
c
                       do ip = 1,rcpplots
c
                           tvals(in-icnt,ip) = tvals(in,ip)
c
                       end do
c
                    endif
c
                 end do
c
c                Fix rcpvals which defines number of data points
c
                 rcpvals = rcpvals - icnt
c

               endif

c
c              Store RCP data for each offset into plotting arrays
c

               do ip2 = 1,3
c
c                  ip = (in1-1)*3 + ip2  
c
                  ip = 2*ip2 + in1 - 2  

                  pnks(ip,2) = rcpvals
c
                  do ik = 1,rcpvals
c
                     mouts(ik,ip,2) = touts(ik)
c
c                    Load pressure instead of Ti for 3rd element
c
                     if (ip2.eq.3) then   
                        mvals(ik,ip,2) = tvals(ik,ip2+1)
                     else
                        mvals(ik,ip,2) = tvals(ik,ip2)
                     endif 
c
                  end do
c
               end do
c
            end do
c
         endif
c
c
c        Generate plots
c
c        Place the contents of corresponding plots on a common plot
c        i.e. plots of ne, Te, Ti and Pressure.
c
c
c
         ELABS(1) = 'OSM OSM '
         Len = lenstr(datatitle)
         ELABS(2) = 'EXPTEXPT'//datatitle(1:Len)
c
         do in = 1,2
c
c            ip = (in-1)*3
c
            ip = in
c
c           Initialize
c
            pltlabs(ip) = ' '
            pltlabs(ip+2) = ' '
            pltlabs(ip+4) = ' '
c
c           Set plot labels
c
c            write(pltlabs(ip+1),'(a,f5.3)') 'Ne: SHIFT=',
c     >                               exp_offsets(in)
c            write(pltlabs(ip+2),'(a,f5.3)') 'Te: SHIFT=',
c     >                               exp_offsets(in)
c            write(pltlabs(ip+3),'(a,f5.3)') 'Ti: SHIFT=',
c     >                               exp_offsets(in)
c            write(pltlabs(ip+4),'(a,f5.3)') 'Pr: SHIFT=',
c     >                               exp_offsets(in)

            write(pltlabs(ip),'(a,f7.4)') 'n_e(m^-3): SHIFT=',
     >                               exp_offsets(in)
            write(pltlabs(ip+2),'(a,f7.4)') 'T_e(eV): SHIFT=',
     >                               exp_offsets(in)
c            write(pltlabs(ip+3),'(a,f7.4)') 'T_i: SHIFT=',
c     >                               exp_offsets(in)
            write(pltlabs(ip+4),'(a,f7.4)') 'p(Pa): SHIFT=',
     >                               exp_offsets(in)
c
         end do
c
c         PLTLABS(1) =  'Ne         '
c         PLTLABS(2) =  'Te         '
c         PLTLABS(3) =  'Ti         '
c         PLTLABS(4) =  'Pressure   '
c
         if (int_type.eq.1) then 
            XLAB = 'R (M) OUTER MIDPLANE'
         elseif(int_type.eq.4.or.int_type.eq.5) then
            XLAB = 'RHO (M)'
         elseif (int_type.eq.2) then
            XLAB = 'R (M) ALONG PROBE LINE'
         elseif (int_type.eq.3) then
            XLAB = 'Z (M) ALONG PROBE LINE'
         endif
c
         YLAB   = ' '
c
         NPLOTS = NPLOTS + 1
         WRITE (IPLOT,9012) NPLOTS,REF
c
c slmod begin - new
c...     Sets some formatting options:
         slopt3 = 3
c slmod end
         sctype = iopt
         ngrm  = 6
         nplts = 6
c
c       Set to plot only OSM or OSM + RCP data
c
        if (exp_ds.gt.0) then
           ngs = 2
        else
           ngs = 1
        endif
c
         write(6,*) 'NGS:',ngs,exp_ds
c
c        Set drawing style for the different sets of data
c
         drawtype(1) = 1
         drawtype(2) = 2
c slmod begin - new
c...     Load fluid code data from E2Dxxx arrays:
c
c        Only load if data exists
c
         IF (cre2d.gt.0) THEN
           ngs = ngs + 1
           CALL RZero(touts2,MAXTHE)
           CALL RZero(tvals2,MAXTHE*MAXNGS)
           call fluidprobe(tvals2,touts2,osmvals,osmplots,
     >                     r1p,z1p,r2p,z2p,int_type,crmb,1.0)
           do in = 1,2
              do ip2 = 1,3
c
c                ip = (in-1)*3+ip2
c
                 ip = 2*ip2 + in - 2  

                 pnks(ip,ngs) = osmvals
                 do ik = 1,osmvals
                    mouts(ik,ip,ngs) = touts2(ik)
c
c                   Load pressure instead of Ti
c
                    if (ip2.eq.3) then 
                       mvals(ik,ip,ngs) = tvals2(ik,ip2+1)
                    else 
                       mvals(ik,ip,ngs) = tvals2(ik,ip2)
                    endif  
                 end do
              end do  
           end do

           ELABS   (ngs) = 'UE  UEDGE'  
           drawtype(ngs) = 1
         ENDIF

         Len = lenstr(datatitle)
         extra_comments(1) = datatitle(1:len)         

         CALL SetPlotComments(iref,job,extra_comments,1,0.0)
c slmod end

c
c        Set pltfact to a value that will result in a slight extension 
c        of the X-axis to make the plotting look neater.
c
         pltfact = 1.02
c
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
               mdrawtype(ip,ig) = drawtype(ig)
            end do  
            pngs(ip) = ngs
         end do
c
         CALL DRAWM (MOUTS,MWIDS,MVALS,Maxdatx,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,ngs,
     >              mdrawtype,1)
c
c        Restore colours
c
         call setup_col(n_cols,col_opt)
c
c
c
c-----------------------------------------------------------------------
c
c     RCP - reciprocating probe comparison plots
c
      elseif (iref.eq.357) then
c
c        The input data for the RCP probe contains the following
c        information in the following format:
c
c '000 OSM PROBE' R1P  Z1P  R2P  Z2P  INT_TYPE
c '000 RCP PROBE' EXP_DS EXP_DATAOPT EXP_VCALCOPT EXP_TCALCOPT  EXP_PARAM
c                        EXP_OFFSET
c
c
c     Tasks required:
c
c     1) Extract the ne,Te,Ti and Pressure from the grid along the
c        probe length. Map these as a function of the mid-plane
c        coordinate of the corresponding rings. At present, only
c        intersections with the main SOL are calculated - the
c        core and private plasma rings are not calculated.
c
c     2) Read in the experimental probe data - Jsat (1,2), Te and
c        the outer mid-plane coordinates for each set of data.
c        Using the exp_calcopt specified convert the data to
c        ne, Te, Ti and Pressure. Apply X-axis offset.
c
c        - the experimental probe data is specified by a set of inputs
c          - one for each panel on the plot - so 4 sets altogether. 
c        - each line of experimental data sets may consists
c          of an integer specifying how many experimental data sets
c          will be plotted and then the index of each dataset in the 
c          experimental data file. 
c        - NOTE: Most of RCP probe data input is NOT used at this 
c          time for this plot.  
c
c
c     3) Plot each of the corresponding sets of data on the same
c        graph. This will require some work in specifying the bounds
c        of the plot region as well as the fact that the two plots
c        will not share a common set of X-axis values.
c
c        Extract OSM data
c 

         WRITE (REF,'(''OSM PROBE PLOT:'',4(1x,f6.3))')
     >                                r1p,z1p,r2p,z2p
c
         call osmprobe(lvals,louts,osmvals,osmplots,
     >                 r1p,z1p,r2p,z2p,int_type,crmb,qtim)
c
c        Do NOT load standard RCP data - this plot can NOT use it at this time 
c        This plot can only be generated vs. PSIN and the data must be in the 
c        experimental data file in that format. 
c
c        Extract RCP data - if dataset specified
c
c         if (exp_ds.gt.0) then
c            call rcpprobe(tvals,touts,rcpvals,rcpplots,
c     >                 exp_ds,exp_offset,exp_dataopt,
c     >                 exp_vcalcopt,exp_tcalcopt,exp_param,
c     >                 lvals,louts,osmvals,osmplots,
c     >                 int_type,rizb,crmb,datatitle,r1p,z1p,r2p,z2p)
c
c            WRITE (REF,'(''351: EXPT/OSM PROBE PLOT: OFFSET ='',f9.5)')
c     >                                exp_offset
c
c         endif
c
c
c        Generate plots
c
c        Place the contents of corresponding plots on a common plot
c        i.e. plots of ne, Te, Ti and Pressure.
c
c
c
c
         PLTLABS(1) =  'Ne         '
         PLTLABS(2) =  'Te         '
         PLTLABS(3) =  'Ti         '
         PLTLABS(4) =  'Pressure   '
c
         if (int_type.eq.1) then 
            XLAB = 'R (M) OUTER MIDPLANE'
         elseif(int_type.eq.4.or.int_type.eq.5) then
            XLAB = 'RHO (M)'
         elseif (int_type.eq.2) then
            XLAB = 'R (M) ALONG PROBE LINE'
         elseif (int_type.eq.3) then
            XLAB = 'Z (M) ALONG PROBE LINE'
         elseif (int_type.eq.6) then
            XLAB = 'PSIN'
         endif
c
         YLAB   = 'Probe Results'
c
         NPLOTS = NPLOTS + 1
         WRITE (IPLOT,9012) NPLOTS,REF
c
         CALL rzero (mvals, maxdatx*maxngs*maxplts)
         CALL rzero (mouts, maxdatx*maxngs*maxplts)
         CALL izero (pnks , maxngs*maxplts)
c
         sctype = iopt
         ngrm  = 3
         nplts = 3
c
c        Set OSM labels
c
         do ip = 1,nplts
c
            mlabs(ip,1) = 'OSM OSM '
            pngs(ip) = 1
c
         end do   
c
c        Load experimental data sets for plotting
c
         ierr = 0
c
c        Read in as many lines specifying experimental
c        data for each plot as are present. Reads until
c        start of next plot. (new IREF).
c
         do while (ierr.eq.0) 
c
            call rdg_expt(graph,plotid,nexpt,maxexpt,expt_ds,
     >                    expt_col,ierr)
c
            if (ierr.eq.0) then         
               do  in = 1,nexpt
c
                   pngs(plotid) = pngs(plotid) +1
c
                   call get_mexpt_data(mvals,mouts,mlabs,
     >                  pnks,plotid,pngs(plotid),expt_ds(in),
     >                  int_type,expt_col(in),0)
c
               end do 
            end if
         end do
c
c        Reset I/O error flag
c
         ierr = 0 
c
c        LOAD OSM data into mvals array as well. 
c
c        Loop through the data for the 4 plots
c
c        Loading one set of information for each plot for OSM
c        Experimental data is already loaded. 
c
         do ip = 1,nplts
c
c           Load OSM first
c
            pnks(ip,1) = osmvals
c
            do ik = 1,osmvals
c
               mouts(ik,ip,1) = louts(ik)
c
               mvals(ik,ip,1) = lvals(ik,ip)
c
            end do
c
         end do
c
c        Set drawing style for the different sets of data
c
c        Manually set the drawing styles for now. 
c
c        Ne 
c
         mdrawtype(1,1) = 1
         mdrawtype(1,2) = 1
         mdrawtype(1,3) = 1
c
c        Te
c
         mdrawtype(2,1) = 1
         mdrawtype(2,2) = 1
         mdrawtype(2,3) = 1
c
c        Ti 
c
         mdrawtype(3,1) = 1
         mdrawtype(3,2) = 2
c
c        Pressure
c
c         mdrawtype(4,1) = 3
c
C
c slmod begin - new
c         Len = lenstr(datatitle)
c         do in = 1,min(len/20+1,8)
c            extra_comments(in) = 
c     >             datatitle((in-1)*20+1:min(in*20,len))         
c         end do  
c
c         in = min(len/20+1,8)
c
c         CALL SetPlotComments(iref,job,extra_comments,in,0.0)
c slmod end
c
c        Set pltfact to a value that will result in a slight extension 
c        of the X-axis to make the plotting look neater.
c
         pltfact = 1.02
c
c
c        Set up data for modified call to DRAWM
c
c         do ip = 1,nplts
c            do ig = 1, pngs(ip)
c               write(0,'(''MLABS'',3i4,a)') ip,ig,
c     >                mdrawtype(ip,ig),mlabs(ip,ig) 
c            end do  
c         end do
c
         CALL DRAWM (MOUTS,MWIDS,MVALS,Maxdatx,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,ngs,
     >              mdrawtype,1)
c
c
c--------------------------------------------------------------------
c
c        Collector probe estimates
c
      elseif (iref.eq.361) then 
c
c        model collector probe
c

         call collector_probe(r1p,z1p,r2p,z2p,probe_diameter,
     >                        probe_dperp,axis_opt,iopt)


      ELSEIF (IREF.EQ.391) THEN
C
C     LINE-INTEGRATED PIN H-ALPHA EMISSION
C     3-d array
C     VERIFY PARAMETERS
C
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
C           WRITE(6,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
C  NEUTRAL!
         IF (IZMIN.NE.0) IZMIN = 0
         iz1=1
C  IZMAX IGNORED
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
C
C  NORMALISE DIRECTIONAL COSINE VECTORS V1 and V2 (extremes of fan)
C  IN CASE THE INPUT WAS NOT NORMALIZED
C
         c1=(cx1*cx1 + cy1*cy1 + cz1*cz1)**0.5
         if(c1.ne.0.)then
           cx1=cx1/c1
           cy1=cy1/c1
           cz1=cz1/c1
         endif
c
         c2=(cx2*cx2 + cy2*cy2 + cz2*cz2)**0.5
c
         if(c2.ne.0.)then
           cx2=cx2/c2
           cy2=cy2/c2
           cz2=cz2/c2
         endif
c
c  vp, cross product between v1 and v2
c
         vp(1)=cy1*cz2 - cz1*cy2
         vp(2)=cz1*cx2 - cx1*cz2
         vp(3)=cx1*cy2 - cy1*cx2
         vnor=(vp(1)*vp(1)+vp(2)*vp(2)+vp(3)*vp(3))**0.5
         vp(1)=vp(1)/vnor
         vp(2)=vp(2)/vnor
         vp(3)=vp(3)/vnor
C
C  CREATE DIRECTIONAL COSINE ARRAY
C
         DO II = 1, NUMTHE
C
c
c  dot product between v1 and v2
c
           dot=cx1*cx2 + cy1*cy2 + cz1*cz2
c
c  angle th3 between v1 and v2
c
           if(dot.gt.1.)dot=1.
           if(dot.lt.-1.)dot=-1.
           th3=acos(dot)
           write(6,*)' fan spans ',th3,'radians'
c
c  angular increment dthe
c
           dthe=th3/(float(numthe-1))
           touts(ii)=dthe*(ii-1)
           TWIDS(II) = THERES
           write (6,*) 'touts:',ii,touts(ii),twids(ii),dthe,th3
c
           if(numthe.eq.1)then
             vn(1)=cx1
             vn(2)=cy1
             vn(3)=cz1
           elseif(ii.eq.1)then
             vn(1)=cx1
             vn(2)=cy1
             vn(3)=cz1
           elseif(ii.eq.numthe)then
             vn(1)=cx2
             vn(2)=cy2
             vn(3)=cz2
           else
c
c  for each viewing line, ii, calculate the central
c  directional cosine vector, vn, from the three conditions:
c  (1) v1*vn = cos (thn)
c  (2) v2*vn = cos (th3-thn)
c  (3) vp*vn = 0
c
             cthn=cos(touts(ii))
c
             if(cthn.gt.1.)cthn=1.
             if(cthn.lt.-1.)cthn=-1.
c
             thn=touts(ii)
             cthc=cos(th3-thn)
c
             if(vp(3).ne.0.)then
               cmatr(1,1)=cx1*vp(3) - cz1*vp(1)
               cmatr(1,2)=cy1*vp(3) - cz1*vp(2)
               cmatr(1,3)=cthn*vp(3)
               cmatr(2,1)=cx2*vp(3) - cz2*vp(1)
               cmatr(2,2)=cy2*vp(3) - cz2*vp(2)
               cmatr(2,3)=cthc*vp(3)
               detm=cmatr(1,1)*cmatr(2,2)-cmatr(2,1)*cmatr(1,2)
               if(detm.eq.0.)then
               else
                 vn(1)=(cmatr(2,2)*cmatr(1,3)
     >                 -cmatr(1,2)*cmatr(2,3))/detm
                 vn(2)=(cmatr(1,1)*cmatr(2,3)
     >                 -cmatr(2,1)*cmatr(1,3))/detm
                 vn(3)=-1.*(vp(1)*vn(1) + vp(2)*vn(2))/vp(3)
               endif
c
             elseif(vp(2).ne.0.)then
c
               cmatr(1,1)=cx1*vp(2) - cy1*vp(1)
               cmatr(1,2)=cz1*vp(2) - cy1*vp(3)
               cmatr(1,3)=cthn*vp(2)
               cmatr(2,1)=cx2*vp(2) - cy2*vp(1)
               cmatr(2,2)=cz2*vp(2) - cy2*vp(3)
               cmatr(2,3)=cthc*vp(2)
               detm=cmatr(1,1)*cmatr(2,2)-cmatr(2,1)*cmatr(1,2)
c
               if(detm.eq.0.)then
               else
                 vn(1)=(cmatr(2,2)*cmatr(1,3)
     >                 -cmatr(1,2)*cmatr(2,3))/detm
                 vn(3)=(cmatr(1,1)*cmatr(2,3)
     >                 -cmatr(2,1)*cmatr(1,3))/detm
                 vn(2)=-1.*(vp(1)*vn(1) + vp(3)*Vn(3))/vp(2)
               endif
c
             else
               vn(1)=0.
               detm=cy1*cz2 - cz1*cy2
               if(detm.eq.0.)then
                 write(6,*)' 3-d array ill-defined'
               else
                 vn(2)=(cthn*cz2 - cthc*cz1)/detm
                 vn(3)=(cthc*cy1 - cthn*cy2)/detm
               endif
             endif
           endif
c
c  about each central vector define an array of vectors to
c  simulate the finite angular resolution, theres.
c
c
c  vt, cross product between vn and vp
c
           vt(1)=vn(2)*vp(3) - vn(3)*vp(2)
           vt(2)=vn(3)*vp(1) - vn(1)*vp(3)
           vt(3)=vn(1)*vp(2) - vn(2)*vp(1)
c
c  the array will span both directions orthogonal to the viewing
c  line, vp and vt, with M vectors in each direction ( 2M vectors
c  per viewing line) and weights WRES corresponding to the solid
c  angle associated with each vector.
c
           Mav=2*avpts - 1
           fm=float(mav)
           phi = THERES/(fm)
c
           do IMav = 1,Mav
             km = 2*(imav-avpts)
             fkm=abs(float(km))
             tphi = tan( km*phi )
             vnor = (1+tphi*tphi)**0.5
c
             do i3 =1,3
               couts(ii,imav,i3)  = (vn(i3)+tphi*vp(i3))/vnor
               couts(ii,imav+mav,i3) = (vn(i3)+tphi*vt(i3))/vnor
             enddo
c
             if(km.ne.0)then
               wres(ii,imav)=((fkm+1.)*(fkm+1.)
     >                       -(fkm-1.)*(fkm-1.))/(4.*fm*fm)
             else
               wres(ii,imav)=1./(2.*fm*fm)
             endif
             wres(ii,imav+mav)=wres(ii,imav)
           ENDDO
         ENDDO
c
C
C  LABELS AND SCALE FACTORS
C
         YLAB = 'PIN H-ALPHA (PH M-2 S-1)'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(''RESOLUTION='',F12.6,'' DEGREES'')') THERES
         IF (ATYPE.EQ.0) THEN
           MFACT = 1.0
           WRITE(ANLY,'(''NO SCALE FACTOR APPLIED'')')
         ELSEIF (ATYPE.EQ.1) THEN
           MFACT = DTHE / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = '',G12.6,'' / (2*PI)'')')
     >                DTHE
         ELSEIF (ATYPE.EQ.2) THEN
           MFACT = 1.0 / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (2*PI)'') ')
         ELSEIF (ATYPE.EQ.3) THEN
           MFACT = 1.0 / ( 4.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (4*PI)'')')
           YLAB = 'PIN H-ALPHA (PH M-2 S-1 SR-1)'
         ENDIF
C
C  LOAD INTEGRAND
C
         CALL RVALKR (CVALSA,PINALPHA,1,1,1,FT,FP,MFACT,
     >     XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
c
         write (6,*) 'pinalpha:',vmin,vmax

C
C  INTEGRATE
C
         write(6,*)' calling LOS3D'
         CALL LOS3DA(TVALS(1,IZ1),robs,0.0,zobs,couts,wres,cvalsa,
     >     numthe,avpts)
C
         WRITE (IPLOT,9012) NPLOTS,REF
         WRITE (IPLOT,*) NVIEW
         WRITE (IPLOT,*) PLANE
         WRITE (IPLOT,*) ANLY
c
         write (6,*) 'Plot:'
         do ii = 1,numthe
            write (6,*) touts(ii),twids(ii),touts(ii)*raddeg,
     >                  tvals(ii,iz1)
            touts(ii) = touts(ii) * raddeg
         end do
c
         themin=touts(1)
         themax=touts(numthe)
C
C  DON'T PLOT SINGLE LINES OF SIGHT
C
         IF (NUMTHE.GT.1) THEN
           CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,1,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,ZLABS(IZMIN),REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)
         ELSE
           WRITE(iplot,*) IZMIN,TOUTS(1),TVALS(1,IZ1)
         ENDIF
C
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
C
      endif 

      return
c
c     Fortmat statements
c

 9012 FORMAT(1X,'PLOT',I3,4X,A)
c
      end
c
c
c
      subroutine get_mexpt_data(mvals,mouts,mlabs,pnks,
     >                          ip,ig,expt_ds,axis_req,
     >                          data_req,proc_data)
      implicit none
      include 'params'
c
c     For multiple plots on the page - series 700
c
      real mvals(maxdatx,maxplts,maxngs)
      real mouts(maxdatx,maxplts,maxngs)
      character*36 mlabs(maxplts,maxngs)
      integer pnks(maxplts,maxngs)
c
      integer  ip,ig,expt_ds
      integer  axis_req,data_req,proc_data
       
c
c     This routine loads the experimental dataset specified by expt_ds
c     into the mvals,mouts arrays at the location ip,ig. It also sets the 
c     value of pnks(ip,ig) to correspond to the number of data points in the 
c     data read in.  
c
c     MVALS    : Data to be plotted
c     MOUTS    : X-axis data 
c     MLABS    : Labels/titles for each graph on each plot
c     PNKS     : Number of data points on each graph
c     IP       : Index of Plot to load
c     IG       : Index of graph on specified plot to load
c     EXPT_DS  : Index of experimental dataset to load into IP,IG
c     AXIS_REQ : Type of axis requested for plot - the experimental
c                data will need to be converted to the requested
c                axis type if it is not already in a compatible 
c                format.  
c     DATA_REQ : If multiple sets of experimental data are on each line
c                in the experimental data file - this specifies which
c                column to load. Default value is 1 being the first 
c                dataset after the index and axis data. 
c     PROC_DATA: Flag to indicate any extra processing that must 
c                be performed on the experimental data to make it
c                suitable for plotting.                 
c     
c
c     NOTE: Axis conversion is NOT yet supported - all data is expected
c           to already be in the required format. 
c
c
      integer dataunit,maxcols
      parameter (dataunit=13,maxcols=8)
c
      integer in   
      integer axis_type
      integer num_expt,ncols,ndata
      real expt_axis(maxdatx) 
      real expt_data(maxdatx,maxcols)
      character*60 datatitle
c
      integer start,end,extstr
      external extstr
c
      call load_expt_data(dataunit,expt_ds,expt_axis,axis_type,
     >                    expt_data,maxcols,maxdatx,
     >                    num_expt,ncols,datatitle)

c
c     Insert call to convert axes as necessary 
c
c
c      call convert_axis(expt_data,expt_axis,maxdatx,
c     >                  maxcols,num_expt,ncols,
c     >                  axis_type,axis_req)
c

c
c      Copy the loaded experimental data into the mvals arrays  
c
       end = extstr(datatitle,start)
       mlabs(ip,ig) =  datatitle(start:start+3)
     >               //datatitle(start:start+3)
       pnks(ip,ig)  = num_expt
c
       do in = 1,num_expt
c 
          mouts(in,ip,ig) = expt_axis(in)
          mvals(in,ip,ig) = expt_data(in,data_req)
c
       end do
c
c
c
       return
       end

