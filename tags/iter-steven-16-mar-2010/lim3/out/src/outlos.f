c     -*-Fortran-*-
c
      subroutine plot_los(iselect,istate,npts,nlines,
     >                  iexpt,iaxis,iavg,ifact,optval,graph_org,
     >                  iopt,job,title,table,avs,navs,nplots,
     >                  iplot,nizs,ierr) 
      implicit none
c
      include 'params'
c
c     The PSIN plotting feature from DIVIMP has been left out of
c     LIM since it isn't meaningful at the present time. 
c
c      include 'psin_data' 
c
      include 'adas_data_spec'
c
      integer iselect,istate,npts,nlines,iexpt,iaxis,ifact,ierr 
      integer iavg
      real optval 
      character*(*) job,title,table,graph_org
      integer navs,iopt,nplots,iplot,nizs
      real avs(0:navs)
c
c     PLOT_LOS: This routine plots a generalized LOS plot. Each LOS
c               plot may have an individual R,Z,Theta,dTheta and
c               length/distance over which the LOS will be integrated.
c
c               The code is based on the usual LOS integration routines 
c               except that only a single LOS value is calculated at a 
c               time. There are a number of axis options possible
c               including - R, Z, Theta, R-intersection and Z-intersection.
c
c     NOTE: IAVG is used to control the weighting applied to the sub-LOS that
c           are used to calculate the LOS for a given angle. A value of zero
c           weights all LOS equally while a value of 1 will result in a 
c           circular weighting.   
c
c
c
c     Local variables
c
      integer maxip
      parameter(maxip=100)
      integer nr,nz,nt,ndt,nd,nds,nca,ndg_d,ndg_w,ndg_l
      real rvals(maxip),zvals(maxip),theta(maxip),dtheta(maxip),
     >     dists(maxip),scalefs(maxip),custom_axis(maxip)
      real geofact_d(maxip),geofact_w(maxip),geofact_l(maxip)
C
      real tmp_store
c      real psin  
      integer targ,itarg,icnt
c
      real r,z,thet,dthet,dist,losval,losval2,losvalexpt,scalef
      real g_d,g_w,g_l
c
      real tmpplot(maxnxs,maxnys),tmpplot2(maxnxs,maxnys)
      real mfact
c
c     Plot labels   
c
      character*36 blabs(2),xlab,ylab,ylab2,blabs2(2)
      character*50 ref,plane,anly,nview,ref2
      character*32 datatitle
      character*36 custom_xaxis_label 
      character*256 tmpjob
c
      INTEGER IGNORS(MAXNGS),ngs
c
      REAL TOUTS(MAXTHE),TWIDS(MAXTHE),TVALS(MAXTHE,MAXNGS)
      real tmptvals(0:maxthe-1)
      REAL THEMIN,THEMAX
      integer itec,ismoth
c
      integer len,lenstr,ii,in
c
c     Variables for 2D experimental data - iexpt = -1
c
      character*256 graph
      character*256 filename
      integer  maxix,maxiy,nix,niy,ifnopt
      parameter(maxix=100,maxiy=100)
      real expt_array(maxix,maxiy),raxis(maxix),zaxis(maxiy)
      real lim2D_array(maxix,maxiy)
      real expt_rmin,expt_rmax,expt_zmin,expt_zmax,expt_dr,expt_dz
c
      integer iexpt_out
      integer itype
c
c     This flag is used to distinguish between different 
c     ways of interpreting the input data. Default value is zero - if
c     Dist is interpreted as DTHET then input_type is set to 1. 
c
      integer input_type 
c
c     Variables for surface distances
c
      integer max_ints
      parameter(max_ints=10)
      integer n_ints, n_int_closest
      real  x_ints(max_ints),y_ints(max_ints),dist_ints(max_ints)
      real zero_shift
c
c     -----------------------------------------------
c
c     End of declarations
c
c      write(0,*) 'Starting outlos:'
c      write(6,'(a,i6)') 'OUTLOS:IAXIS=',iaxis
c
c     Initialization
c
c      external_psin = .false.
c
      if (nplots.eq.1) then  
         call calc_lim_poly
      endif

c
      targ      = 0
      itarg     = 0    
      iexpt_out = 0
      input_type= 0
c
      do ii = 1,maxngs 
         ignors(ii) = 1
      end do
c
c     Number of plots initially set to 1.
c
      ngs =1 
      mfact = 1.0
c
c     Set NVIEW, PLANE ...
c
      NVIEW  = 'GENERALIZED LOS/INTEGRATION PLOT'
      plane = ' '
      REF = ' ' 
      YLAB = ' '
c
c     Based on the values of iselect and istate - load up the required 
c     2D array for the LOS plot 
c
c     Scaling by absfac (if required) is performed in the 
c     load_divdata_array routine. 
c
      call rzero(tmpplot,maxnxs*maxnys)
c
      itype = 1  
c
c     For Steradian weighted plots - load a different label - assuming LOS here. 
c
      if (ifact.eq.3) itype = 2
c
      call load_limdata_array(tmpplot,iselect,istate,itype,
     >                         ylab,blabs(1),ref,nizs,ierr)
c
      write(6,*) 'LOADING LIMDATA:' ,ierr
c
c     If ISELECT=22 then also load the emission data unmodified 
c     by the temperature for the denominator in a line-averaged
c     impurity ion temperature plot.
c
      if (iselect.eq.22) then 
c
         call rzero(tmpplot2,maxnxs*maxnys)
c
c        Use ADAS data loaded on the previous call  
c
         cadas_switch = 1 
c
         call load_limdata_array(tmpplot2,4,istate,itype,
     >                      ylab2,blabs2(1),ref2,nizs,ierr)
c
c        Reset the switch to load previously stored ADAS data
c
         cadas_switch = 0
c  
         write(6,*) 'LOADING LIMDATA2:' ,ierr
c
      endif

c
c     Check return code from data load and exit if non-zero   
c
      if (ierr.ne.0) return 
c
c     --------- OPTION REMOVED FOR LIM -------------------------
c
c     Check to see if a LOS integration of a 2D experimental data array is 
c     also requested.
c
c     A value of IEXPT of -2 means to also convert the DIVIMP data 
c     to the same format as the experimental data prior to the 
c     calculation of the values to be plotted.
c
c      if (iexpt.eq.-1.or.iexpt.eq.-2) then  
c
c        Set to 2 plots on page
c 
c         ngs =2  
c
c         blabs(2) = 'EXPT EXPT'
c
c         call rdfn(graph,filename,ifnopt,ierr)
c
c         if (ierr.ne.0) then  
c
c            len = lenstr(filename)
c            write(0,*) 'PLOT_LOS:ERROR READING 2D EXPT'
c     >                   //' DATAFILE INPUT LINE: '
c     >                   //filename(1:len),ifnopt,ierr
c            iexpt = 0
c            ngs = 1
c
c         else 
c
c            call load_2Dexpt_data(filename,datatitle,
c     >                  maxix,maxiy,nix,niy,expt_array,raxis,zaxis,
c     >                  expt_rmin,expt_rmax,expt_zmin,expt_zmax,
c     >                  expt_dr,expt_dz,ierr)
c
c            if (ierr.ne.0) then  
c
c               len = lenstr(filename)
c               write(0,*) 'PLOT_LOS:ERROR LOADING 2D EXPT DATAFILE: '
c     >                   //filename(1:len)
c               iexpt = 0
c               ngs = 1
c
c            else
c
c               len=lenstr(datatitle)
c               blabs(2) = 'EXPT '//datatitle(1:len)
c
c            endif 
c
c         endif
c
c        Convert the DIVIMP array to an experimental data format
c 
c         if (iexpt.eq.-2) then 
c
c            call create_2Ddata(tmpplot,maxnxs,maxnys,
c     >                  maxix,maxiy,nix,niy,
c     >                  lim2D_array,raxis,zaxis,0)
c
c         endif
c
c      endif
c
c     --------- OPTION REMOVED FOR LIM -------------------------
c
c     Now that the data has been successfully loaded - load the 
c     data describing the various lines of sight. 
c
c     If an error occurs in reading the input data then exit immediately.
c
c     Load R-values 
c
      call RDG_REAL_ARRAY(GRAPH,rvals,maxip,nr,ierr)
      if (ierr.ne.0) return
c
c     Load Z-values 
c
      call RDG_REAL_ARRAY(GRAPH,zvals,maxip,nz,ierr)
      if (ierr.ne.0) return
c
c     Load Theta-values 
c
      call RDG_REAL_ARRAY(GRAPH,theta,maxip,nt,ierr)
      if (ierr.ne.0) return
c
c     Load Delta-theta-values 
c
      call RDG_REAL_ARRAY(GRAPH,dtheta,maxip,ndt,ierr)
      if (ierr.ne.0) return
c
c     Load length-values 
c
      call RDG_REAL_ARRAY(GRAPH,dists,maxip,nd,ierr)
      if (ierr.ne.0) return
c
c     Load Additional Scaling factors
c
      call RDG_REAL_ARRAY(GRAPH,scalefs,maxip,nds,ierr)
      if (ierr.ne.0) return
c
c     Load customized axis coordinates if iaxis=9
c
      if (iaxis.eq.9) then 
         call RDG_REAL_ARRAY(GRAPH,custom_axis,maxip,nca,ierr)
         if (ierr.ne.0) return
         len = lenstr(graph)
         custom_xaxis_label = graph(4:len)
      endif 
c
c     Get filename and load data for an externally specified PSIN mapping
c
c      if (iaxis.eq.10) then 
c
c         call get_psin_filename(ierr)
c         if (ierr.ne.0) return 
c
c        Load PSIN reference data
c
c         call read_psin_data
c
c      endif
c
c
c     --------- OPTION REMOVED FOR LIM -------------------------
c
c     Load additional data required for geometry scaling option 5
c
c     These inputs describe the geometry of the DIIID bolometers and 
c     allow for the calculation of additional scaling factors. 
c
c     When these values are entered the array of dT values is overwritten
c     with values based on these quanitities.  
c
c      if (ifact.eq.5.or.ifact.eq.6) then  
c         call RDG_REAL_ARRAY(GRAPH,geofact_d,maxip,ndg_d,ierr)
c         if (ierr.ne.0) return
c         call RDG_REAL_ARRAY(GRAPH,geofact_w,maxip,ndg_w,ierr)
c         if (ierr.ne.0) return
c         call RDG_REAL_ARRAY(GRAPH,geofact_l,maxip,ndg_l,ierr)
c         if (ierr.ne.0) return
c
c        Recalculate the dtheta values specified based on the D,W,L values input
c
c         ndt = npts
c
c         write(6,'(a,i5)') 'OUTLOS: Recalculating DTHETA',ndt
c
c         do in = 1,npts
c
c            ii = min(in,ndg_d)
c            g_d = geofact_d(ii)
c
c            ii = min(in,ndg_w)
c            g_w = geofact_w(ii)
c
c            ii = min(in,ndg_l)
c            g_l = geofact_l(ii)
c             
c            dtheta(in) = 2.0 * atan ((g_d+g_w)/(2.0*g_l)) * raddeg
c
c            write(6,'(a,i5,5(1x,g12.5))') 'DT:',in,dtheta(in),
c     >                    g_d,g_w,g_l 
c
c         end do
c
c      endif
c
c     --------- OPTION REMOVED FOR LIM -------------------------
c
c     All data to calculate LOS plot is now available - loop through 
c     calculating values. 
c
      IF (ifact.LT.0 .OR. ifact.GT.10) ifact = 0
c
c     Set LOS scaling factor and note value for plot
c
      if (iavg.eq.0) then 
         write(anly,'(''RECT VIEW: '')')
      elseif (iavg.eq.1) then  
         write(anly,'(''CIRC VIEW: '')')
      endif
c
      IF (ifact.EQ.0) THEN
         mfact = 1.0
         WRITE(ANLY(12:),'(''NO GEO FACTOR APPLIED'')')
      ELSEIF (ifact.EQ.1) THEN
         WRITE(ANLY(12:),'(''GEO FACTOR = DTHETA/(2*PI)'')')
      ELSEIF (ifact.EQ.2) THEN
         MFACT = 1.0 / (2.0 * PI)
         WRITE(ANLY(12:),'(''GEO FACTOR = 1 / (2*PI)'') ')
      ELSEIF (ifact.EQ.3) THEN
         MFACT = 1.0 / ( 4.0 * PI)
         WRITE(ANLY(12:),'(''GEO FACTOR = 1 / (4*PI)'')')
      ELSEIF (ifact.EQ.4) THEN
         WRITE(ANLY(12:),'(''GEO FACTOR=TAN(DTHETA/2)^2/PI'')')
      ELSEIF (ifact.EQ.5) THEN
         WRITE(ANLY(12:),'(''2D GEO FACTORS FROM W,D,L'')')
      ELSEIF (ifact.Eq.6) THEN
         WRITE(ANLY(12:),'(''GEO FACTORS FROM W,D,L'')')
      ELSEIF (ifact.EQ.7) THEN
         WRITE(ANLY(12:),'(''MOD GEO FACTORS FROM W,D,L'')')
      ENDIF
c
c
C     X-axis labels
C
c     Determined by value of iaxis
c
c     1 = Plot versus theta
c     2 = Plot versus R-intersection
c     3 = plot versus Z-intersection
c     4 = plot versus R-observation
c     5 = plot versus Z-observation 
c     6 = Index number
c     7 = Plot versus Psin
c     8 = Plot versus channel number = Index number + opt_val
c     9 = Custom axis - user specified input
c    10 = Plot versus Psin calculated from externally loaded data    
c    11 = Plot as distance along LIMITER
c    12 = Plot as distance along LIMITER with center LOS at 0.0
c
      if (iaxis.eq.1) then 
c
         XLAB = 'THETA (DEGREES)'
c
      elseif (iaxis.eq.2) then 
c
         write(XLAB,
     >       '(''R-Interection (M) (with Z='',f8.4,'')'')') optval
c
      elseif (iaxis.eq.3) then 
c
         write(XLAB,
     >       '(''Z-Interection (M) (with R='',f8.4,'')'')') optval
c
      elseif (iaxis.eq.4) then 
c
         XLAB = 'R (M)'
c
      elseif (iaxis.eq.5) then 
c
         XLAB = 'Z (M)'
c
      elseif (iaxis.eq.6) then 
c
         XLAB = 'Index Number'
c
      elseif (iaxis.eq.7.or.iaxis.eq.10) then 
c
         XLAB = 'Psin'
c
      elseif (iaxis.eq.8) then 
c
         XLAB = 'Channel Number'
c
      elseif (iaxis.eq.9) then
c
         len = lenstr(custom_xaxis_label) 
         XLAB = custom_xaxis_label(1:len)
c
      elseif (iaxis.eq.11) then 
c
         XLAB = 'ALONG SURFACE'
c
      elseif (iaxis.eq.12) then 
c
         XLAB = 'ALONG SURFACE (CENTER=0)'
c
      endif     
c
c     If the scale factor has been set to anything other than 1.0 - then 
c     adjust the BLAB value. 
c
      if (nds.ne.1.or.(nds.eq.1.and.scalefs(1).ne.1.0)) then
         call set_blab(iselect,istate,2,nizs,blabs(1))
      endif
c
c
C     INTEGRATE
C
c     Perform a series of single LOS integrations  
c
c     Use the global npts that was specified in the first line - if
c     there are less that npts values specified for any line - the 
c     last value on each input line will be used for all remaining
c     lines of sight. 
c         
      do in = 1,npts 
c
c        Calculate the R, Z, Theta and Dtheta values for each 
c        LOS.
c                   
c        R-obs 
c
         ii = min(in,nr)
         r = rvals(ii)  
c
c        Z-obs
c
         ii = min(in,nz)
         z = zvals(ii)  
c
c        Theta
c
         ii = min(in,nt)
         thet = theta(ii)  
c
c        DTheta
c
         ii = min(in,ndt)
         dthet = dtheta(ii)  
c
c        Dist
c
         ii = min(in,nd)
         dist = dists(ii)  
c
c        Additional Scaling factor
c
         ii = min(in,nds)
         scalef = scalefs(ii)  
c
c        Additional Geomtery Factors
c
c slmod begin
c
c bugish - Getting an out-of-bounds error when running with IFACT=3, since NDG_D, NDG_W and NDG_L
c          are not assigned above. -SL (Jan 20, 2003)
c
c         IF (ndg_d.EQ.0) THEN
c           ii = in
c         ELSE
c           ii = min(in,ndg_d)
c         ENDIF
c         g_d = geofact_d(ii)  
c
c         IF (ndg_d.EQ.0) THEN
c           ii = in
c         ELSE
c           ii = min(in,ndg_w)
c         ENDIF
c         g_w = geofact_w(ii)  
c
c         IF (ndg_d.EQ.0) THEN
c           ii = in
c         ELSE
c           ii = min(in,ndg_l)
c         ENDIF
c         g_l = geofact_l(ii)  
c
c         ii = min(in,ndg_d)
c         g_d = geofact_d(ii)  
c
c         ii = min(in,ndg_w)
c         g_w = geofact_w(ii)  
c
c         ii = min(in,ndg_l)
c         g_l = geofact_l(ii)  
c slmod end
c
c
c        NOTE: IF the LOS data specified contains only ONE value for 
c              Theta, ONE value for DTheta, ONE value for Robs and ONE
c              value for Zobs and yet npts is greater than ONE - then 
c              the code assumes that instead of a series of degenerate/
c              identical views - what is actually wanted is a series of
c              views starting at THETA with a separation of Dtheta and 
c              a width of DIST.
c
         if (npts.gt.1.and.
     >       nr.eq.1.and.nz.eq.1.and.nt.eq.1.and.ndt.eq.1) then 
c
             thet  = thet + (in-1) * dthet
             dthet = dist
             dist  = 0.0
             input_type = 1 
c
         endif            
c
c        Set scale factor based on viewing width 
c
         if (ifact.eq.1) then 
c
c            Convert Dthet to radians for this scaling
c
             MFACT = (abs(DTHET)*degrad) / (2.0 * PI)
c
         elseif (ifact.eq.4) then 
c
c            Convert Dthet to radians for this scaling
c
             MFACT = tan((abs(DTHET)/2.0)*DEGRAD)**2 / PI
c
c        INCORRECT SCALING - JUST FOR REFERENCE
c
         elseif (ifact.eq.5) then 
c
c            This is the scaling factor for the second dimension of the cross-section
c
             MFACT = tan(abs(DTHET)*DEGRAD) / PI
c
         elseif (ifact.eq.6) then 
c
c            Convert Dthet to radians for this scaling
c
             MFACT = (abs(DTHET)*degrad) / (4.0 * PI)
c
         elseif (ifact.eq.7) then 
c
c            Convert Dthet to radians for this scaling
c
             MFACT = 1.0 / (4.0 * PI)
c
         endif  
c
c         write (6,'(a,i4,5g12.4)') 'LOS:',in,r,z,thet,dthet,dist
c         write (0,'(a,i4,5g12.4)') 'LOS:',in,r,z,thet,dthet,dist
c          
c        If using DIVIMP data converted to experimental data format
c
c         if (iexpt.eq.-2) then 
c
c            CALL LOSINTEXPT(losvalexpt,thet,dthet,1,
c     >               R,Z,nlines,dist,
c     >               maxix,maxiy,nix,niy,lim2D_array,raxis,zaxis)
c
c            tvals(in,1) = losvalexpt * mfact * scalef
c
c
c        DIVIMP data in regular or DIVIMP internal format 
c
c         else
c
c            write (6,'(a,i6)') 'LOS: ifact = ',ifact
c 
c            if (ifact.eq.5.or.ifact.eq.6.or.ifact.eq.7) then
c
c               CALL LOSINT_SCALE(losval,thet,dthet,1,
c     >                  R,Z,nlines,tmpplot,dist,g_d,g_w,g_l,ifact)
c
c            else
c 
               CALL LOSINT(losval,thet,dthet,1,
     >                  R,Z,nlines,tmpplot,dist,iavg)

c
c              Load denominator for line averaged quantities
c
               if (iselect.eq.22) then    

                  CALL LOSINT(losval2,thet,dthet,1,
     >                  R,Z,nlines,tmpplot2,dist,iavg)

c
c                 Scale losval2 appropriately
c
                  losval2 = losval2 * mfact * scalef  
c
               endif
c
c            endif
c
c           Assign values to plot array
c
            tvals(in,1) = losval * mfact * scalef
c
            if (iselect.eq.22.and.losval2.ne.0.0) then

               tvals(in,1) = tvals(in,1) / losval2 

            endif  
c
c         endif
c
c
c     --------- OPTION REMOVED FOR LIM -------------------------
c
c        If LOS over 2D experimental data is required
c
c         if (iexpt.eq.-1.or.iexpt.eq.-2) then 
c
c            CALL LOSINTEXPT(losvalexpt,thet,dthet,1,
c     >               R,Z,nlines,dist,
c     >               maxix,maxiy,nix,niy,expt_array,raxis,zaxis)
c
c            tvals(in,2) = losvalexpt * mfact * scalef
c
c         endif
c
c     --------- OPTION REMOVED FOR LIM -------------------------
c
c        Calculate axis values
c
c        Theta 
c
         if (iaxis.eq.1) then  
c
            touts(in) = thet
            twids(in) = abs(dthet)
c
c        R-intersection 
c
         elseif (iaxis.eq.2) then 
c   
            call adjustout(thet,1,optval,r,z)
c
            touts(in) = thet
            twids(in) = 1.0 
c
c        Z-intersection - needs testing
c
         elseif (iaxis.eq.3) then 
c   
            call adjustoutz(thet,1,optval,r,z)
c
            touts(in) = thet
            twids(in) = 1.0 
c
c        R-observation 
c
         elseif (iaxis.eq.4) then 
c   
            touts(in) = r
            twids(in) = 1.0 
c
c        Z-observation 
c
         elseif (iaxis.eq.5) then 
c   
            touts(in) = z
            twids(in) = 1.0 
c
c        Index number
c
         elseif (iaxis.eq.6) then 
c   
            touts(in) = in
            twids(in) = 1.0 
c
c         elseif (iaxis.eq.7.or.iaxis.eq.10) then 
c
c            call calcpsin(r,z,thet,psin,targ)
c
c            if (targ.ne.0.and.itarg.eq.0) then
c               itarg = targ
c            endif 
c
c           ONLY LOS that strike one target are desired - the code
c           assumes that the first LOS that intersects a target defines
c           the target of interest. All other points will be excluded. 
c           LOS that do strike the target are assigned a proper psin value
c           - all others will get a psin value of -100.0 - all these 
c           points will be filtered out before plotting.
c
c
c            if (itarg.ne.0.and.targ.eq.itarg) then  
c
c               touts(in) = psin
c               twids(in) = 1.0 
c
c            else
c
c               touts(in) = -100.0
c               twids(in) =  1.0 
c
c            endif      
c
c            write(6,'(a,2i6,f12.6)') 'ITARG:',targ,itarg,psin
c
c
c        Index number minus 1 - corresponds to channel number 0+
c
         elseif (iaxis.eq.8) then 
c   
            touts(in) = in + optval
            twids(in) = 1.0 
c
c        Custom/user specified axis - must be ordered - ascending or descending - not mixed
c 
         elseif (iaxis.eq.9) then  
c
c           There should be an axis value specified for every point - but protect
c           against an error just in case.  
c
            ii = min(in,nca)
            touts(in) = custom_axis(ii)
c
c        Find the distance along the LIMITER/TARGET surface for the chord
c
         elseif (iaxis.eq.11.or.iaxis.eq.12) then 
c
            call limiter_intersection(r,z,thet,x_ints,y_ints,dist_ints,
     >                                max_ints,n_ints,n_int_closest)
c
c           Assign touts to the closest dist_ints value to the viewpoint
c           if it exists - otherwise just assign theta.
c
            if (n_int_closest.gt.0) then 
               touts(in) = dist_ints(n_int_closest)
            else
               touts(in) = thet
            endif 

         endif         
c
c        End of LOS calculation loop 
c
      end do
c
c     If plotting versus the psin scale - filter out all points with a 
c     PSIN value that is set to -100.0 - these points either do not 
c     intersect a target or are not on the default target defined by the 
c     first LOS to intersect a target. 
c
c      if (iaxis.eq.7.or.iaxis.eq.10) then 
c
c         icnt = 0
c
c         do in = 1,npts
c
c            if (touts(in).le.-100.0) then
c
c               icnt = icnt + 1
c
c            else
c              
c               touts(in - icnt) = touts(in)
c               twids(in - icnt) = touts(in)
c               tvals(in - icnt,1) = tvals(in,1)
c
c            endif
c
c         end do 
c
c         write(6,'(a,3i6)') 'OUTLOS:PSI REMOVED:',npts,icnt,npts-icnt 
c
c         npts = npts - icnt
c
c
c      endif 
c
c     Reorder results so that axis is in ascending order
c
      if (touts(1).gt.touts(npts)) then
c
c        reverse order of touts and tvals entries.
c        
         do in = 1,npts/2
c
c           axis
c
            tmp_store = touts(in)
            touts(in) = touts(npts-in+1)
            touts(npts-in+1) = tmp_store
c
c           Divimp values
c
            tmp_store = tvals(in,1)
            tvals(in,1) = tvals(npts-in+1,1)
            tvals(npts-in+1,1) = tmp_store
c
c           Experimental values
c
            tmp_store = tvals(in,2)
            tvals(in,2) = tvals(npts-in+1,2)
            tvals(npts-in+1,2) = tmp_store

         end do   
c
      endif
c
c     For IAXIS = 12 - readjust the axis so that the center
c     LOS has a value of 0.0
c
      if (iaxis.eq.12) then 
c
        zero_shift = touts(npts/2)
c
        do in = 1,npts  
c
           touts(in) = touts(in) - zero_shift
c
        end do  
c
      endif
c
c     Set max and min of plotting range
c
      themin = touts(1)
      themax = touts(npts) 
c
c     Load value into PLANE if req'd 
c
      write(plane,'(i4,'' SUB-LOS '')') nlines
c
      if (nds.eq.1.and.scalefs(1).ne.1.0) then 
c
          write(PLANE(13:),'(a,g12.5)')
     >   ': SCALING FACTOR=',scalefs(1)
c
      elseif (nds.gt.1) then  
c
          write(PLANE(13:),'(a)')
     >   ': SCALING FACTORS APPLIED'
c
      elseif (input_type.eq.0.and.
     >        nd.eq.1.and.dists(1).gt.0.0) then
c
         write(PLANE(13:),'(a,f6.3,a)')
     >   'INTEG OVER ',dists(1),' (M) ONLY'
c
      endif   
c  
c     Document plot information
c
      len = lenstr(ref)
      WRITE (iplot,9012) NPLOTS, REF(1:len)
      len = lenstr(blabs(1)) 
      WRITE (iplot,'(a)') BLABS(1)(1:len)
      WRITE (iplot,*) NVIEW
      WRITE (iplot,*) PLANE
      WRITE (iplot,*) ANLY
c
c     Load experimental data to be plotted
c
c     A value greater than zero loads from the regular experimental 
c     data file.
c
c
      if (iexpt.gt.0.and.iexpt.le.100) then 
c
          iexpt_out = 0
c
          ngs = 2
c
          call calc_expt(iexpt,touts,tvals,maxthe,npts,maxdatx,
     >                   themin,themax,maxngs,ngs,datatitle)

          BLABS(2) = 'EXPT '//DATATITLE
c
      elseif (iexpt.gt.100) then  
c
c         Instructs OUT to load and plot the experimental data without
c         trying to interpolate or draw straight lines.  
c
          iexpt_out = iexpt-100      
c
      endif
c
C
C     DON'T PLOT SINGLE LINES OF SIGHT - Just write single numbers to file
C
      IF (Npts.GT.1) THEN
c
          do in = 1,npts
             write(6,*) 'TVALS:',in,touts(in),tvals(in,1),
     >                                  tvals(in,2)
          end do
c
          itec = 1
          ismoth = 99          
c
c         Save value of JOB and store plot reference in first half 
c
          tmpjob = job
          len = lenstr(graph_org)
          job(1:36) = graph_org(5:min(len,41))
c
          CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NPTS,ANLY,NGS,
     >              ISMOTH,THEMIN,THEMAX,-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >              JOB,TITLE,XLAB,YLAB,BLABS,REF,NVIEW,PLANE,
     >              TABLE,IOPT,2,1.0,iexpt_out)
c
c         Restore value of JOB
c
          job = tmpjob 
c
      ELSE
c
         len = lenstr(blabs(1))
c
         if (ngs.eq.1) then  

            WRITE(iplot,'(a,a,2(a,i4),a,f8.3,a,g12.5,a,g12.5)')
     >                'RESULT LOS: ',blabs(1)(1:len),
     >                'PLOT =',Iselect,
     >                ' STATE = ',istate,
     >                ' POS = ',TOUTS(1),
     >                ' VAL = ',TVALS(1,1)
         else

            WRITE(iplot,'(a,a,2(a,i4),a,f8.3,a,g12.5,a,g12.5)')
     >                'RESULT LOS: ',blabs(1)(1:len),
     >                'PLOT =',Iselect,
     >                ' STATE = ',istate,
     >                ' POS = ',TOUTS(1),
     >                ' VAL = ',TVALS(1,1),
     >                ' EXPT = ',tvals(1,2)
         endif
c
      ENDIF
c
c
c     Return from generalized LOS routine
c
      return
c
c     Format statements 
c

 9012 FORMAT(1X,'PLOT',I3,4X,A)
c
      end
c
c
c
      SUBROUTINE RDG_LOS(GRAPH,npts,nlines,iselect,
     >                   istate,iexpt,iaxis,iavg,ifact,optval,ierr)
      IMPLICIT  NONE
      INTEGER   Iselect,IERR,istate,npts,nlines,iaxis,iavg,ifact,iexpt
      real optval  
      CHARACTER GRAPH*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  RDG_LOS  : LOAD BASE PARAMETERS FOR GENERALIZED LOS PLOTS        *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
c
      ierr = 0
C
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDGLOS'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
c      IF (BUFFER(2:2).EQ.'#') THEN
c        CALL Read_AdditionalPlotData(BUFFER)
c        GOTO 100
c      ENDIF
c
c     If it is not a '000' designation 
c
      if (buffer(2:4).ne.'000') then
         ierr = 1
         rewind(5)
         return
      endif
C
      MESAGE = 'EXPECTING 1 CHARACTER STRING, 8 INTEGERS AND 1 REAL'
c
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,iselect,istate,
     >                         npts,nlines,iaxis,
     >                         iexpt,iavg,ifact,optval
      RETURN

 9998 ierr = 1
       
      return 
C
 9999 IERR = 1
      WRITE (7,'(1X,2A,3(/1X,A))')
     >  'RDG_LOS: ERROR READING ',GRAPH,MESAGE,
     >  'LAST LINE READ :-',BUFFER
      RETURN
      END
c
c
c
      SUBROUTINE RDG_REAL_ARRAY(GRAPH,data,maxpts,npts,ierr)
      IMPLICIT  NONE
      INTEGER   npts,ierr,maxpts
      real data(maxpts) 
      CHARACTER GRAPH*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  RDG_REAL_ARRAY  : READS A 1D ARRAY OF REAL DATA ITEMS            *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
      integer in
c
      ierr = 0
C
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDG_RA'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
c      IF (BUFFER(2:2).EQ.'#') THEN
c        CALL Read_AdditionalPlotData(BUFFER)
c        GOTO 100
c      ENDIF
c
c     If it is not a '000' designation 
c
      if (buffer(2:4).ne.'000') then
         ierr = 1
         rewind(5)
         return
      endif
C
      MESAGE = 'EXPECTING 1 CHARACTER STRING AND A REAL DATA ARRAY'
c
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,npts
c
      if (npts.le.maxpts) then 
c
         READ (BUFFER,*,ERR=9999,END=9999) GRAPH,npts,
     >                (data(in),in=1,npts)
c
      else
c 
         MESAGE = 'NUMBER OF DATA ITEMS TOO LARGE FOR ARRAY'
         goto 9999
c
      endif
c
      RETURN

 9998 ierr = 1
      WRITE (6,'(1X,2A,3(/1X,A))')
     >  'RDG_REAL_ARRAY: ERROR READING ',GRAPH,MESAGE,
     >  'LAST LINE READ :-',BUFFER
      WRITE (7,'(1X,2A,3(/1X,A))')
     >  'RDG_REAL_ARRAY: ERROR READING ',GRAPH,MESAGE,
     >  'LAST LINE READ :-',BUFFER
       
      return 
C
 9999 IERR = 1
      WRITE (7,'(1X,2A,3(/1X,A))')
     >  'RDG_REAL_ARRAY: ERROR READING ',GRAPH,MESAGE,
     >  'LAST LINE READ :-',BUFFER
      WRITE (6,'(1X,2A,3(/1X,A))')
     >  'RDG_REAL_ARRAY: ERROR READING ',GRAPH,MESAGE,
     >  'LAST LINE READ :-',BUFFER
      RETURN
      END
C
C
C
      SUBROUTINE RDG1 (GRAPH,ADASID,ADASYR,ADASEX,
     >                 ISELE,ISELR,ISELX,ISELD,IERR)
      implicit none
      INTEGER   ISELE,ISELR,ISELX,ISELD,IERR,ADASYR
      CHARACTER GRAPH*(*), ADASID*(*),ADASEX*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  RDG1 : READ IN SELECTOR SWITCHES FOR ADAS PLRP CALCULATIONS      *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
C
      IERR = 0
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDG1'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c     Feature Only useful in OUT
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
c      IF (BUFFER(2:2).EQ.'#') THEN
c        CALL Read_AdditionalPlotData(BUFFER)
c        GOTO 100
c      ENDIF
c
c      write(0,'(a,8i5)')
c     >  'RDG1:',len(adasid),len(adasex),adasyr,isele,iselr,iselx
C
      MESAGE = 'EXPECTING 2 CHAR, 1 INT, 1 CHAR  AND 4 INTEGERS'
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,ADASID,ADASYR,ADASEX,
     >                                  ISELE,ISELR,ISELX,ISELD
c
c      write(0,'(a,8i5)')
c     >  'RDG1:',len(adasid),len(adasex),adasyr,isele,iselr,iselx
c
c      write(0,'(3a)')
c     >  'RDG1:',buffer,':'
c      write(0,'(3a)')
c     >  'RDG1:',graph,':'
c      write(0,'(3a)')
c     >  'RDG1:',adasid,':'
c      write(0,'(3a)')
c     >  'RDG1:',adasex,':'
c

      RETURN
C
 9998 IERR = 1
      WRITE (6,'(1X,A,4(/1X,A))')
     >  'RDG1: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',BUFFER
      WRITE (7,'(1X,A,4(/1X,A))')
     >  'RDG1: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',BUFFER
      RETURN
C
 9999 IERR = 1
      WRITE (6,'(1X,A,4(/1X,A))')
     >  'RDG1: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',BUFFER
      WRITE (7,'(1X,A,4(/1X,A))')
     >  'RDG1: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',BUFFER
      RETURN
      END
c
c
c
      SUBROUTINE RDG_plot (GRAPH, IPLOT,name, ierr)
      IMPLICIT  none
      INTEGER   IPLOT,IERR                 
      CHARACTER NAME*(*),GRAPH*(*)                                              
C                                                                               
C***********************************************************************        
C                                                                               
C         THIS ROUTINE READS IN A LINE OF GRAPH DETAILS                         
C                                                                               
C     PARAMETERS :-                                                             
C     GRAPH  : TITLE OF GRAPH                                                   
C     IPLOT  : PLOTTING OPTIONS   (1 INTEGERS)                           
C     NAME   : FOR PRINTING IN ANY ERROR MESSAGES                               
C     IERR   : SET TO 1 IF AN ERROR FOUND                                       
C                                                                               
C       DON'T WANT EOF TO BE TREATED AS AN ERROR IN THIS CASE.  STILL           
C     SET IERR TO 1 BUT DON'T PRINT ANY MESSAGES (HENCE GOTO 9998)              
C                                                                               
C        CHRIS FARRELL    JAN 1988                                              
C
C       ADDED INTEGRATION RANGE DATA. VALUES OF 0.0 TO 0.0 ARE 
C     GENERALLY TREATED AS THE FULL RANGE.
C
C       DAVID ELDER,  OCT 29, 1990
C
C                                                                               
C***********************************************************************        
C                                                                               
      INCLUDE   'reader'                                                        
C     INCLUDE   (READER)                                                        
      CHARACTER COMENT*72,MESAGE*72                                             
C                                                                               
      MESAGE = 'END OF FILE ON UNIT 5'                                          
  100 IF (IBUF.EQ.0) READ (5,'(A72)',ERR=9998,END=9998) BUFFER                  
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDGPLT'                                   
      IF (BUFFER(1:1).EQ.'$') GOTO 100                                          
C                                                                               
      MESAGE = 'EXPECTING CHARACTER STRING, 1 INTEGER'             
C
C     REPLACE LIST-DIRECTED I/O ON INTERNAL FILES WITH A BACKSPACE
C     FOLLOWED BY A READ ON THE EXTERNAL FILE THAT IS THE SOURCE FOR 
C     THE BUFFER. THIS IS DONE BECAUSE AT THIS POINT IN TIME THE CFT 
C     AND CFT77 COMPILERS ON THE CRAY DO NOT SUPPORT THIS FEATURE.
C     
C                      DAVID ELDER NOV.15,1989
C
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH                                   
C
c      BACKSPACE(5)
c      READ (5,*,ERR=9999,END=9999) GRAPH                                   
c
      IF (graph(1:3).ne.'LOS'.and.graph(1:3).ne.'000') THEN                       
        IBUF = 1                                                                
      ELSE                                                                      
        IBUF = 0                                                                
C
c       jdemod - the code originally tried to read
c                graph and iplot on all lines - this
c                failed for certain lines - now the 
c                code checks for LOS before trying 
c                to read iplot
c
c        READ (BUFFER,*,ERR=9999,END=9999) GRAPH,IPLOT
c
        if (graph(1:3).eq.'000') then 
           iplot = 0  
        elseif (graph(1:3).eq.'LOS') then 
           READ (BUFFER,*,ERR=9999,END=9999) GRAPH,IPLOT
        endif
c
C
c        BACKSPACE(5)
c        READ (5,*,ERR=9999,END=9999) GRAPH,VMIN,VMAX,
c     >    GRIMIN,GRIMAX,IPLOT,JSMOTH,         
c     >    MAXIZ,IPLANE,IFOLD,IALLIZ,IVU                                         
c
      ENDIF                                                                     
      RETURN                                                                    
C                                                                               
 9998 IERR = 1                                                                  
c      write(0,*) 'RDG_PLOT: ERR 9998'
      RETURN                                                                    
C                                                                               
 9999 IERR = 1                                                                  
c      write(0,*) 'RDG_PLOT: ERR 9999'
      WRITE (7,'(1X,2A,3(/1X,A))')                                              
     >  'RDG_PLOT: ERROR READING ',
     >    NAME,MESAGE,'LAST LINE READ :-',BUFFER            
      RETURN                                                                    
      END                                                                       
c
c
c
      subroutine limiter_intersection(r,z,theta,x,y,dist,
     >                                max_ints,n_ints,n_closest)
      implicit none
      include 'params'
      include 'comxyt'
      include 'comt2'
      include 'comtor' 
c
      integer max_ints,n_ints,n_closest
      real r,z,theta
      real x(max_ints),y(max_ints),dist(max_ints)
c
c     limiter_intersection: This routine looks for the intersections 
c                           between the LOS given by r,z,theta and 
c                           the limiter edge defined by QEDGES(iqx,2) 
c                           where the second index defines the two sides
c                           of the limiter shape. 
c
c                           When it finds these intersections it returns
c                           the X-value, Y-value and QDIST (or distance 
c                           along the limiter from (0.0,0.0)). These
c                           values can be used as the ordinate in plotting
c                           routines. 
c
c                           This routine uses INTSECT2DP to calculate the 
c                           intersections.  
c
c     Y < 0 is in QEDGES(IQX,1)         
c     Y > 0 is in QEDGES(IQX,2)
c     
c     Dimensions:  QEDGES(-MAXQXS:0,2) 
c                  QDISTS(-MAXQXS:0,2)
c                  QXS(-MAXQXS:MAXQXS)
c     
c     NQXSO - number of outboard elements
c     IQX = 1-NQXSO to 0
c     
c     NOTE ON COORDINATE SYSTEMS:
c     
c     Usually a tokamak uses an R,Z coordinate system for radial 
c     and vertical coordinates. LIM uses an X,Y system. Depending on 
c     how you orient the LIM grid - the definition of X,Y relative to the 
c     R,Z system may change. In the LIM case, X is the perpendicular to 
c     field line coordinate while Y is the parallel to field line direction. 
c     When translating this to the R,Z system used by the code imported from
c     DIVIMP it has been chosen that the Y coordinate corresponds to R while
c     the X one corresponds to Z. This may cause some confusion.   
c     
c     
c     Variables for intectdp - in double precision  
c     
      real*8 ra,za,rb,zb,r1,z1,r2,z2,rint,zint
      integer sect
c     
      integer in
      real tmp_dist,min_dist,seglen,frac
      real max_dim
      real*8 length
c     
c     Initialize 
c     
      n_ints = 0
      n_closest = 0
      min_dist = hi
c     
c     MAX_SIZE defines the approximate scale size of the modeled region      
c     
      max_dim = max(2.0*cl,ca-caw)
c     
c     LENGTH defines the distance from the R,Z observation position to the 
c     calculated second endpoint of the line. It should be far enough 
c     away that is is garanteed to be on the other side of the region 
c     being modeled in LIM. 
c     
      length = sqrt(r**2+z**2) + 2.0 * max_dim
c     
      ra = r
      za = z
c     
      rb = length * cos(theta*degrad)  
      zb = length * sin(theta*degrad)
c     
c     Loop through limiter surfaces
c     
      do in = 1-nqxso,0 
c     
         r1 = -qedges(in,1)
         z1 = qxs(in) 
c     
c        Look at side 1 - Y < 0 side - first
c     
c        Z1 and Z2 retain the same values on both sides of the 
c        Limiter surface. 
c     
         if (in.eq.0) then 
c     
c           Limiter always ends at the tip at 0.0,0.0
c     
            r2 = 0.0
            z2 = 0.0  
c     
         else         
c     
            r2 = -qedges(in+1,1)
            z2 = qxs(in+1) 
c     
         endif
c     
c        Calculate length of element of limiter
c     
         seglen = (r2-r1)**2 + (z2-z1)**2
c     
c        Check for intersection on non-zero segment lengths    
c     
         if (seglen.gt.0.0) then 
c     
            call intsect2dp(ra,za,rb,zb,r1,z1,r2,z2,rint,zint,sect)
c     
c           Intersection found  
c     
            if (sect.eq.1) then 
c     
               if (n_ints.lt.max_ints) then  
                  n_ints = n_ints+1
c     
                  tmp_dist = sqrt((ra-rint)**2 + (za-zint)**2)
c     
                  if (tmp_dist.lt.min_dist) then 
                     n_closest = n_ints
                     min_dist = tmp_dist
                  endif  
c     
                  frac =  sqrt((rint-r1)**2+(zint-z1)**2)
c     
                  x(n_ints) = zint
                  y(n_ints) = rint
                  dist(n_ints) = qdists(in,1) - frac
c     
               else
                  write(6,*) 'ERROR:LIMITER_INTERSECTIONS: NUMBER OF'//
     >                    ' INTERSECTIONS EXCEEDS MAXIMUM = ',max_ints
c     
               endif 
            endif
         endif
c     
c        Now look at side 2 - Y > 0
c     
         r1 = qedges(in,2)           
c     
         if (in.eq.0) then 
c     
c           Limiter always ends at the tip at 0.0,0.0
c     
            r2 = 0.0
c     
         else         
c     
            r2 = qedges(in+1,2)
c     
         endif
c     
c        Calculate length of element of limiter
c     
         seglen = (r2-r1)**2 + (z2-z1)**2
c     
c        Check for intersection on non-zero segment lengths    
c     
         if (seglen.gt.0.0) then 
c     
            call intsect2dp(ra,za,rb,zb,r1,z1,r2,z2,rint,zint,sect)
c     
c           Intersection found  
c     
            if (sect.eq.1) then 
c     
               if (n_ints.lt.max_ints) then 
                  n_ints = n_ints+1
c     
                  tmp_dist = sqrt((ra-rint)**2 + (za-zint)**2)
c     
                  if (tmp_dist.lt.min_dist) then 
                     n_closest = n_ints
                     min_dist = tmp_dist
                  endif  
c     
                  frac =  sqrt((rint-r1)**2+(zint-z1)**2)
c     
                  x(n_ints) = zint
                  y(n_ints) = rint
                  dist(n_ints) = qdists(in,1) - frac
c     
               else
                  write(6,*) 'ERROR:LIMITER_INTERSECTIONS: NUMBER OF'//
     >                    ' INTERSECTIONS EXCEEDS MAXIMUM = ',max_ints
      
               endif
            endif 
         endif
c     
      end do
c      
c      write(6,*) 'DEBUG: LIMITER_INTERSECTIONS:',n_ints,n_closest
c      do in = 1,n_ints
c
c         write(6,'(a,i5,3(1x,g12.5))') 'INT:',in,x(in),y(in),dist(in)
c
c      end do 
c
      return
      end  


