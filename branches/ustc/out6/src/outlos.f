c     -*-Fortran-*-
c 
      subroutine plot_los(iselect,istate,npts,nlines,
     >                  iexpt,iaxis,iavg,ifact,optval,graph_org,
     >                  iopt,job,title,table,avs,navs,nplots,
     >                  iplot,nizs,ierr) 
      use subgrid_plots
      implicit none
c
      include 'params'
      include 'psin_data' 
      include 'adas_data_spec'
c
      integer iselect,istate,npts,nlines,iexpt,iaxis,ifact,ierr 
      integer iavg
      real optval 
      character*(*) job,title,table,graph_org
      integer navs,iopt,nplots,iplot,nizs
      real avs(0:navs)
c
c     Variables for subgrids
c     
      integer :: iflag
      real,allocatable :: subgrid_data(:,:)
      real,allocatable :: subgrid_raxis(:),subgrid_zaxis(:)

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
      real psin  
      integer targ,itarg,icnt
c
      real r,z,thet,dthet,dist,losval,losval2,losvalexpt,scalef
      real g_d,g_w,g_l
c
      real tmpplot(maxnks,maxnrs),tmpplot2(maxnks,maxnrs)
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
      real div2D_array(maxix,maxiy)
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
      write(6,'(a,i6)') 'OUTLOS:IAXIS=',iaxis

c
c     Initialization
c
      external_psin = .false.
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
      call rzero(tmpplot,maxnks*maxnrs)
c
      itype = 1  
c
c     For Steradian weighted plots - load a different label - assuming LOS here. 
c
      if (ifact.eq.3) itype = 2
c
c     Load Subgrid data and arrays separately from the regular DIVIMP arrays since the format is different
c     - use allocatable storage for this array
c

      if (iselect.eq.32.or.iselect.eq.33.or.iselect.eq.34.or.
     >    iselect.eq.35) then 
         !
         ! Allocate subgrid storage and load the array data - sg_rdim,sg_zdim are in subgrid_options 
         !
         allocate(subgrid_data(sg_rdim,sg_zdim),stat=iflag)
         allocate(subgrid_raxis(sg_rdim),stat=iflag)
         allocate(subgrid_zaxis(sg_zdim),stat=iflag)

         call load_subgrid_array(subgrid_data,iselect,istate,itype,
     >                         ylab,blabs(1),ref,nizs,ierr)

         call calc_subgrid_axis(subgrid_raxis,sg_rmin,sg_rmax,sg_rdim)
         call calc_subgrid_axis(subgrid_zaxis,sg_zmin,sg_zmax,sg_zdim)

      else

         call load_divdata_array(tmpplot,iselect,istate,itype,
     >                         ylab,blabs(1),ref,nizs,ierr)

      endif
c
      write(6,*) 'LOADING DIVDATA:' ,ierr
c
c     If ISELECT=22 then also load the emission data unmodified 
c     by the temperature for the denominator in a line-averaged
c     impurity ion temperature plot.
c
      if (iselect.eq.22) then 
c
         call rzero(tmpplot2,maxnks*maxnrs)
c
c        Use ADAS data loaded on the previous call  
c
         cadas_switch = 1 
c
         call load_divdata_array(tmpplot2,4,istate,itype,
     >                      ylab2,blabs2(1),ref2,nizs,ierr)
c
c        Reset the switch to load previously stored ADAS data
c
         cadas_switch = 0
c  
         write(6,*) 'LOADING DIVDATA2:' ,ierr
c
      endif

c
c     Check return code from data load and exit if non-zero   
c
      if (ierr.ne.0) goto 200
c
c     Check to see if a LOS integration of a 2D experimental data array is 
c     also requested.
c
c     A value of IEXPT of -2 means to also convert the DIVIMP data 
c     to the same format as the experimental data prior to the 
c     calculation of the values to be plotted.
c
      if (iexpt.eq.-1.or.iexpt.eq.-2) then  
c
c        Set to 2 plots on page
c 
         ngs =2  
c
         blabs(2) = 'EXPT EXPT'
c
         call rdfn(graph,filename,ifnopt,ierr)
c
         if (ierr.ne.0) then  
c
            len = lenstr(filename)
            write(0,*) 'PLOT_LOS:ERROR READING 2D EXPT'
     >                   //' DATAFILE INPUT LINE: '
     >                   //filename(1:len),ifnopt,ierr
            iexpt = 0
            ngs = 1
c
         else 
c
            call load_2Dexpt_data(filename,datatitle,
     >                  maxix,maxiy,nix,niy,expt_array,raxis,zaxis,
     >                  expt_rmin,expt_rmax,expt_zmin,expt_zmax,
     >                  expt_dr,expt_dz,ierr)
c
            if (ierr.ne.0) then  
c
               len = lenstr(filename)
               write(0,*) 'PLOT_LOS:ERROR LOADING 2D EXPT DATAFILE: '
     >                   //filename(1:len)
               iexpt = 0
               ngs = 1
c
            else

               len=lenstr(datatitle)
               blabs(2) = 'EXPT '//datatitle(1:len)
c
            endif 
c
         endif
c
c        Convert the DIVIMP array to an experimental data format
c 
         if (iexpt.eq.-2) then 
c
            call create_2Ddata(tmpplot,maxix,maxiy,nix,niy,
     >                  div2D_array,raxis,zaxis,0)
c
         endif
c
      endif
c
c     Now that the data has been successfully loaded - load the 
c     data describing the various lines of sight. 
c
c     If an error occurs in reading the input data then exit immediately.
c
c     Load R-values 
c
      call RDG_REAL_ARRAY(GRAPH,rvals,maxip,nr,ierr)
      if (ierr.ne.0) goto 200
c
c     Load Z-values 
c
      call RDG_REAL_ARRAY(GRAPH,zvals,maxip,nz,ierr)
      if (ierr.ne.0) goto 200
c
c     Load Theta-values 
c
      call RDG_REAL_ARRAY(GRAPH,theta,maxip,nt,ierr)
      if (ierr.ne.0) goto 200
c
c     Load Delta-theta-values 
c
      call RDG_REAL_ARRAY(GRAPH,dtheta,maxip,ndt,ierr)
      if (ierr.ne.0) goto 200
c
c     Load length-values 
c
      call RDG_REAL_ARRAY(GRAPH,dists,maxip,nd,ierr)
      if (ierr.ne.0) goto 200
c
c     Load Additional Scaling factors
c
      call RDG_REAL_ARRAY(GRAPH,scalefs,maxip,nds,ierr)
      if (ierr.ne.0) goto 200
c
c     Load customized axis coordinates if iaxis=9
c
      if (iaxis.eq.9) then 
         call RDG_REAL_ARRAY(GRAPH,custom_axis,maxip,nca,ierr)
         if (ierr.ne.0) goto 200
         len = lenstr(graph)
         custom_xaxis_label = graph(4:len)
      endif 
c
c     Get filename and load data for an externally specified PSIN mapping
c
      if (iaxis.eq.10) then 
c
         call get_psin_filename(ierr)
         if (ierr.ne.0) goto 200
c
c        Load PSIN reference data
c
         call read_psin_data
c
      endif

c
c     Load additional data required for geometry scaling option 5
c
c     These inputs describe the geometry of the DIIID bolometers and 
c     allow for the calculation of additional scaling factors. 
c
c     When these values are entered the array of dT values is overwritten
c     with values based on these quanitities.  
c
c     Initialization required
c      
      ndg_d = 0
      ndg_w = 0
      ndg_l = 0
c
      if (ifact.eq.5.or.ifact.eq.6) then  
         call RDG_REAL_ARRAY(GRAPH,geofact_d,maxip,ndg_d,ierr)
         if (ierr.ne.0) goto 200
         call RDG_REAL_ARRAY(GRAPH,geofact_w,maxip,ndg_w,ierr)
         if (ierr.ne.0) goto 200
         call RDG_REAL_ARRAY(GRAPH,geofact_l,maxip,ndg_l,ierr)
         if (ierr.ne.0) goto 200
c
c        Recalculate the dtheta values specified based on the D,W,L values input
c
         ndt = npts
c
         write(6,'(a,i5)') 'OUTLOS: Recalculating DTHETA',ndt
c
         do in = 1,npts
c
            ii = min(in,ndg_d)
            g_d = geofact_d(ii)
c
            ii = min(in,ndg_w)
            g_w = geofact_w(ii)
c
            ii = min(in,ndg_l)
            g_l = geofact_l(ii)
c             
            dtheta(in) = 2.0 * atan ((g_d+g_w)/(2.0*g_l)) * raddeg
c
            write(6,'(a,i5,5(1x,g12.5))') 'DT:',in,dtheta(in),
     >                    g_d,g_w,g_l 
c
         end do



      endif

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
c 
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
      endif     
c
c     If the scale factor has been set to anything other than 1.0 - then 
c     adjust the BLAB value. 
c
      if (nds.ne.1.or.(nds.eq.1.and.scalefs(1).ne.1.0)) then
         call set_blab(iselect,istate,2,nizs,blabs(1))
      endif
c
c     This is a utility program that reads a set of external 
c     theta values and maps them to psin based on current grid
c     and specified observation position
c
c      call map_file_data_to_psin(rvals(1),zvals(1))
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
         IF (ndg_d.le.0) THEN
           ii = min(in,maxip)
         ELSE
           ii = min(in,ndg_d)
         ENDIF

         g_d = geofact_d(ii)  

         IF (ndg_w.le.0) THEN
           ii = min(in,maxip)
         ELSE
           ii = min(in,ndg_w)
         ENDIF
         g_w = geofact_w(ii)  

         IF (ndg_l.le.0) THEN
           ii = min(in,maxip)
         ELSE
           ii = min(in,ndg_l)
         ENDIF

         g_l = geofact_l(ii)  
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
c
c            If DIST specifies 2 values and this is less than the number
c            of plot points - interpret the second value as an integration
c            limit for distance for all LOS.
c
             if (nd.eq.2.and.nd.ne.npts) then 
                dthet=dists(1)
                dist =dists(2)
             else
                dthet = dist
                dist  = 0.0
             endif

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
         if (iexpt.eq.-2) then 
c
            CALL LOSINTEXPT(losvalexpt,thet,dthet,1,
     >               R,Z,nlines,dist,
     >               maxix,maxiy,nix,niy,div2D_array,raxis,zaxis)
c
            tvals(in,1) = losvalexpt * mfact * scalef
c
c
c        DIVIMP data in regular or DIVIMP internal format 
c
         else

c            write (6,'(a,i6)') 'LOS: ifact = ',ifact
 
            if (ifact.eq.5.or.ifact.eq.6.or.ifact.eq.7) then

               CALL LOSINT_SCALE(losval,thet,dthet,1,
     >                  R,Z,nlines,tmpplot,dist,g_d,g_w,g_l,ifact)


            elseif (iselect.eq.32.or.iselect.eq.33.or.iselect.eq.34.or.
     >              iselect.eq.35)then
               !
               ! LOS Integration over the subgrid data
               !
               CALL LOSINTEXPT(losval,thet,dthet,1,
     >               R,Z,nlines,dist,
     >               sg_rdim,sg_zdim,sg_rdim,sg_zdim,subgrid_data,
     >               subgrid_raxis,subgrid_zaxis)
               

            else
 
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

            endif
c
c           Assign values to plot array
c
            tvals(in,1) = losval * mfact * scalef
c
            if (iselect.eq.22.and.losval2.ne.0.0) then

               tvals(in,1) = tvals(in,1) / losval2 

            endif  
c
         endif
c
c        If LOS over 2D experimental data is required
c
         if (iexpt.eq.-1.or.iexpt.eq.-2) then 
c
            CALL LOSINTEXPT(losvalexpt,thet,dthet,1,
     >               R,Z,nlines,dist,
     >               maxix,maxiy,nix,niy,expt_array,raxis,zaxis)
c
            tvals(in,2) = losvalexpt * mfact * scalef
c
         endif
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
         elseif (iaxis.eq.7.or.iaxis.eq.10) then 
c
            call calcpsin(r,z,thet,psin,targ)
c
            if (targ.ne.0.and.itarg.eq.0) then
               itarg = targ
            endif 
c
c           ONLY LOS that strike one target are desired - the code
c           assumes that the first LOS that intersects a target defines
c           the target of interest. All other points will be excluded. 
c           LOS that do strike the target are assigned a proper psin value
c           - all others will get a psin value of -100.0 - all these 
c           points will be filtered out before plotting.
c
c
            if (itarg.ne.0.and.targ.eq.itarg) then  
c
               touts(in) = psin
               twids(in) = 1.0 
c
            else
c
               touts(in) = -100.0
               twids(in) =  1.0 
c
            endif      
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
      if (iaxis.eq.7.or.iaxis.eq.10) then 
c
         icnt = 0
c
         do in = 1,npts
c
            if (touts(in).le.-100.0) then
c
               icnt = icnt + 1
c
            else
c              
               touts(in - icnt) = touts(in)
               twids(in - icnt) = touts(in)
               tvals(in - icnt,1) = tvals(in,1)
c
            endif
c
         end do 
c
         write(6,'(a,3i6)') 'OUTLOS:PSI REMOVED:',npts,icnt,npts-icnt 
c
         npts = npts - icnt
c
c
      endif 
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
          call calc_expt(iexpt,touts,tvals,maxthe,npts,
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
c     Exit and clean up code
c

 200  continue

c
c     Deallocate storage assigned for subgrid plots if it has been allocated
c     
      if (allocated(subgrid_data)) deallocate(subgrid_data)
      if (allocated(subgrid_raxis)) deallocate(subgrid_raxis)
      if (allocated(subgrid_zaxis)) deallocate(subgrid_zaxis)

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
      subroutine calcpsin(r,z,thet,psin,targ)
      implicit none
      integer targ 
      real r,z,thet,psin
c
      include 'params'
      include 'cgeom' 
      include 'comtor'
      include 'psin_data'
      include 'grbound'
c
c     CALCPSIN: This routine interpolates the target PSI 
c               values to return an estimated PSI value for 
c               a LOS starting at location R,Z with angle thet. 
c
c               This routine looks for target intersections only and
c               returns a success/fail flag in the quantity targ. 
c     
c               targ = 0 - no target intersection found
c               targ = 1 - target number 1
c               targ = 2 - target number 2
c
c               INPUT: R  - R coordinate of viewing location
c               INPUT: Z  - Z coordinate of viewing location
c               INPUT: THET - viewing angle
c               OUTPUT: PSIN - Psi value of target intersection
c               OUTPUT: TARG - LOS intersection flag
c
c     jdemod - this routine has an issue with extended grids that have target elements
c              all around the vessel. It may not find the desired psin intersection point. 
c              As a quick workaround ... the code has been modified to select the 
c              intersection farthest from the observation location. 
c
c
c     Local variables
c
      integer id,iw,ierr
      real rend,zend,r1t,z1t,r2t,z2t,rsect,zsect,scale_dist
c
      real*8 rdp,zdp,renddp,zenddp
      real resultvw
      integer req_int
c
      real targ_psival
      logical check_intsect
      external check_intsect,targ_psival
c
      integer nints,in_maxdist
      integer,parameter :: maxints = 10
      real maxdist
      real int_dist(maxints)
      integer int_targ(maxints)
      real int_psin(maxints)

c
c     Initialization
c      
      nints = 0
      targ =   0
      psin = -100.0
c
c     Calculate an endpoint to the segment that should be on the far
c     side of the targets. 
c
      scale_dist = max(rmax-rmin,zmax-zmin) 
      scale_dist = max(scale_dist,abs(r-r0))
      scale_dist = max(scale_dist,abs(z-z0))
c
      rend = r + 4.0 * scale_dist * cos(thet*DEGRAD)
      zend = z + 4.0 * scale_dist * sin(thet*DEGRAD)
c
c      write(6,'(a,8(1x,g12.5))') 'PSIN:',r,z,rend,zend,
c     >                             rmax,zmax,rmin,zmin
c
c     Use externally specified data to calculate PSIN - requires an
c     additional plot input which gives the required filename. 
c
      if (external_psin) then  
c
         rdp = r
         zdp = z
         renddp = rend
         zenddp = zend

         call calc_intersections(rdp,zdp,renddp,zenddp,.false.)
c
c        If r,z is outside the vessel wall then we want the second
c        intersection (they are sorted by distance) - otherwise we want 
c        the first. 
c         
         
         CALL GA15B(R,Z,RESULTVW,neutwpts,1,NWWORK,4*MAXPTS,
     >            NWINDW,MAXPTS,RNW,ZNW,NWTDUM,NWXDUM,NWYDUM,6)
c
c        Outside
c
         if (resultvw.lt.0.0) then
c
            req_int = 2  
c
c        Inside
c
         else
c
            req_int = 1
c
         endif 
c
c        Assign target and psi information - only working for X-point down grids. 
c
c        Inner target - Xpoint down
c
         if (r_int(req_int).lt.rpsimin) then          
c
            targ = 2
c
c        Outer target - Xpoint down
c
         else
c
            targ = 1
c
         endif
c
c        Assign PSIN
c
         psin = psi_int(req_int)
c
c         write(6,'(a,2i4,8(1x,g12.5))') 'TARG:',req_int,targ,
c     >                       r_int(req_int),z_int(req_int),
c     >                       psin   
c         
c
c     Using PSIN determined from the grid - when external PSIN data is not 
c     available or specified 
c
      else
c
c        The code assumes that a given LOS will only ever strike one target. 
c        NOTE: May not be a good general assumption - ok for now. 
c     
c        jdemod - need quick fix for this assumption. 
c        
c        Count from outside edges of target to try to minimize problems with this 
c        
c        Check target 1 first 
c        
         do id = 1,ndsin
c        
c           Only check for real target elements 
c        
            if (dds(id).gt.0.0) then 
c        
               iw = wallindex(id)            
c        
               if (check_intsect(r,z,rend,zend,
     >                           wallpt(iw,20),wallpt(iw,21),
     >                           wallpt(iw,22),wallpt(iw,23), 
     >                           rsect,zsect)) then 
                  nints = nints+1
                  
                  if (nints.le.maxints) then 
                     int_targ(nints) = 1
                     int_dist(nints) = (rsect-r)**2 + (zsect-z)**2
                     int_psin(nints) = targ_psival(id,rsect,zsect)
                  endif
c                  write(6,'(a,2i8,10(1x,g12.5))') 'PSI:',nints,
c     >                            int_targ(nints),
c     >                            thet,int_psin(nints),
c     >                            rsect,zsect
               endif
c        
            endif  
c        
         end do 
c        
c        Check second target 
c        
         do id = nds, ndsin+1, -1
c        
c           Only check for real target elements 
c        
            if (dds(id).gt.0.0) then 
c        
               iw = wallindex(id)            
c        
               if (check_intsect(r,z,rend,zend,
     >                           wallpt(iw,20),wallpt(iw,21),
     >                           wallpt(iw,22),wallpt(iw,23), 
     >                           rsect,zsect)) then 
                  nints = nints+1
                  if (nints.le.maxints) then 
                     int_targ(nints) = 2
                     int_dist(nints) = (rsect-r)**2 + (zsect-z)**2
                     int_psin(nints) = targ_psival(id,rsect,zsect)
c                  write(6,'(a,2i8,10(1x,g12.5))') 'PSI:',nints,
c     >                            int_targ(nints),
c     >                            thet,int_psin(nints),
c     >                            rsect,zsect
c
                  endif
               endif
c        
            endif  
         
         end do 
c
c        Figure out which intersection point to return 
c
         if (nints.gt.0) then

            in_maxdist = 0
            maxdist = 0

            do id = 1,nints
               if (int_dist(id).gt.maxdist) then 
                  maxdist = int_dist(id)
                  in_maxdist = id
               endif
            end do
            
            targ = int_targ(in_maxdist)
            psin = int_psin(in_maxdist)

         endif
c
      endif

c
c      write(6,'(a,i8,6(1x,g12.5))') 'CALCPSIN:',
c     >             targ,r,z,rend,zend,thet,psin
c
      return
      end
c
c
c
      subroutine calc_psibounds(id,psistart,psimid,psiend,ierr)
      implicit none
      integer id,ierr
      real psistart,psimid,psiend
c
      include 'params'
      include 'cgeom'
      include 'comtor'
c
c     CALC_PSIBOUNDS: This calculates the PSI values at the target element
c                     boundaries and at its centre.
c      
      integer ir,itarg,ida,ira,iw,iwa
      real psia,dista,distb
c
      iw = wallindex(id)
      ir = irds(id)
      itarg = min(int(id/(ndsin+1))+1,2)
c 
      psimid = psitarg(ir,itarg)
c
c     Initialize to mid value 
c
      psistart = psimid
      psiend = psimid
c
c     Check for virtual target element
c
      if (dds(id).eq.0.0) then
         ierr = 1
         return
      endif
c 
c     Find lower boundary psi value 
c         
      ida = id -1
c      
      if (dds(ida).gt.0.0) then 
c
         iwa = wallindex(ida)
         ira = irds(ida)
         psia = psitarg(ira,itarg)
c 
         dista = wallpt(iwa,6)
         distb = wallpt(iw,5)
c
         psistart = psia + (psimid-psia) * dista / (dista+distb)
c
      endif
c 
c     Find upper boundary psi value 
c         
      ida = id + 1
c      
      if (dds(ida).gt.0.0) then 
c
         iwa = wallindex(ida)
         ira = irds(ida)
         psia = psitarg(ira,itarg)
c 
         dista = wallpt(iwa,5)
         distb = wallpt(iw,6)
c
         psiend = psia + (psimid-psia) * dista / (dista+distb)
c
      endif
c
      if (dds(id-1).eq.0.0) then 
c
         psistart = psimid - (psiend-psimid) 

      endif
c
      if (dds(id+1).eq.0.0) then 
c
         psiend = psimid - (psistart-psimid) 

      endif
c      
c      write(6,'(a,i5,8(1x,g12.5))') 'PSI:',id,
c     >                dds(id-1),dds(id),dds(id+1), 
c     >                psistart,psimid,psiend
c
      return
      end
c
c
c
      real function targ_psival(id,rsect,zsect)
      implicit none
      integer id,sect
      real rsect,zsect
c
      include 'params'
      include 'cgeom'
      include 'comtor'
c           
c     TARG_PSIVAL: calculates an approximate value of PSI for 
c                  an intersection point on a given target element. 
c
c
      real psistart,psimid,psiend,rc,zc,re,ze,dista
      integer ierr,iw
c
      call calc_psibounds(id,psistart,psimid,psiend,ierr)
c
      if (ierr.ne.0) then 
         targ_psival = psimid
         return
      endif
c
      iw = wallindex(id) 
c 
      rc = wallpt(iw,1) 
      zc = wallpt(iw,2) 
c
c
c     Check first half  
c
      re = wallpt(iw,20)
      ze = wallpt(iw,21)  
c
      IF ((ABS(zsect-zc).LE.ABS(ze-zc)) .AND.
     >    (ABS(zsect-ze).LE.ABS(ze-zc)) .AND.
     >    (ABS(rsect-rc).LE.ABS(re-rc)) .AND.
     >    (ABS(rsect-re).LE.ABS(re-rc))) then 
c
         dista = sqrt((rsect-re)**2+(zsect-ze)**2) 

         targ_psival = psistart 
     >               + (psimid-psistart) * dista/wallpt(iw,5)
c
c         write(6,'(a,2i5,8(1x,g12.5))') 'TARGPSI:1:',id,iw,
c     >           targ_psival,psistart,psimid,dista,wallpt(iw,5) 
c         write(6,'(a,2i5,8(1x,g12.5))') 'TARGPSI:1:',id,iw,
c     >           rc,rsect,re,zc,zsect,ze 
c
         return
      endif
c
c
c     Check second half  
c
      re = wallpt(iw,22)
      ze = wallpt(iw,23)  
c
      IF ((ABS(zsect-zc).LE.ABS(ze-zc)) .AND.
     >    (ABS(zsect-ze).LE.ABS(ze-zc)) .AND.
     >    (ABS(rsect-rc).LE.ABS(re-rc)) .AND.
     >    (ABS(rsect-re).LE.ABS(re-rc))) then 
c
         dista = sqrt((rsect-re)**2+(zsect-ze)**2) 
c
         targ_psival = psiend 
     >               + (psimid-psiend) * dista/wallpt(iw,5)
c
c         write(6,'(a,2i5,8(1x,g12.5))') 'TARGPSI:2:',id,iw,
c     >           targ_psival,psistart,psimid,dista,wallpt(iw,5) 
c         write(6,'(a,2i5,8(1x,g12.5))') 'TARGPSI:2:',id,iw,
c     >           rc,rsect,re,zc,zsect,ze 
c
         return
c
      endif
c
c     Error condition - code shouldn't reach here
c
      targ_psival = psimid
c
      write(6,'(a,i5,8(1x,g12.5))') 'TARGPSI:ERROR:',id,
     >             rsect,zsect,targ_psival
c
      return
      end 
c
c
c
      subroutine calc_intersections(rstart,zstart,rend,zend,ignore_end)
      implicit none
      real*8 rstart,zstart,rend,zend
      logical ignore_end
c
      include 'params'
      include 'psin_data' 
c
c     CALC_INTERSECTIONS:
c
c     This routine loops through the entire array of surface data and reports all
c     of the intersections of the given line (rstart,zstart) (rend,zend) 
c     with the surface. The ignore_end flag instructs the code to not
c     calculate any intersections for the end-point in the case where it is known
c     to be on the wall. 
c
      integer in,id,sid,next_in
      real*8 rnew,znew
      real*8 tmpdist,tmppsi,tmpr,tmpz
      integer sect
c
      real eps
      parameter (eps=1.0e-6)
c
      n_int = 0
      call dzero(r_int,max_int)
      call dzero(z_int,max_int)
      call dzero(psi_int,max_int)
      call dzero(dist_int,max_int)
c
      do in = 1,n_elements
c
c         Initialize for call - most not needed
c
          rnew = 0.0
          znew = 0.0         
          sect = 0
c
          next_in = in+1
          if (next_in.gt.n_elements) next_in = 1 
c
          CALL INTSECT2DP(rend,zend,rstart,zstart,
     >                    surface_data(in,1),
     >                    surface_data(in,2),
     >                    surface_data(next_in,1),
     >                    surface_data(next_in,2),
     >                    rnew,znew,sect)
c
c         Record any intersections found
c       
          if (sect.eq.1) then 
c
c            If not (endpoint ignored and this is the endpoint)
c
             if (.not.(ignore_end.and.
     >             (abs(rnew-rend).lt.eps.and.
     >              abs(znew-zend).lt.eps))) then 
c
c               Record point
c
                if (n_int.lt.max_int) then 
                   n_int = n_int + 1
                   r_int(n_int) = rnew
                   z_int(n_int) = znew   
c
                   psi_int(n_int) = 
     >                (sqrt((rnew-surface_data(in,1))**2
     >                   +(znew-surface_data(in,2))**2)
     >                   /surface_data(in,4)
     >                 *( surface_data(next_in,3)
     >                   -surface_data(in,3)))
     >                + surface_data(in,3)
c
                   dist_int(n_int) = 
     >                sqrt((rnew-rstart)**2+(znew-zstart)**2)

c
c                   write(6,'(a,3i6,l4,4f12.6)') 'DBG:',n_int,in,
c     >                 next_in,sect,
c     >                 r_int(n_int),z_int(n_int),psi_int(n_int),
c     >                 dist_int(n_int)
c   
c     >                 surface_data(in,1),surface_data(in,2),
c     >                 surface_data(next_in,1),surface_data(next_in,2)
c 
c
                endif 
c
             endif
c
          endif  
c
      end do
c
c     Sort the intersections from closest to farthest away 
c
      if (n_int.gt.1) then  
c
         do in = 1,n_int-1
c      
            sid = in 
c
            do id = in+1,n_int
c      
c              Swap if distance is less to get smallest distance at top
c      
               if (dist_int(id).lt.dist_int(in)) then  
c      
                  sid = id
c      
               end if
c      
            end do
c
            if (sid.ne.in) then 
c      
c              Swap smallest to top
c      
               tmpdist = dist_int(sid) 
               tmppsi  = psi_int(sid) 
               tmpr    = r_int(sid) 
               tmpz    = z_int(sid) 
c      
               dist_int(sid) = dist_int(in)
               psi_int(sid)  = psi_int(in)
               r_int(sid)    = r_int(in)
               z_int(sid)    = z_int(in)
            
               dist_int(in) =  tmpdist
               psi_int(in)  =  tmppsi 
               r_int(in)    =  tmpr   
               z_int(in)    =  tmpz   
c
            endif
c      
         end do

      endif 
c
c     Print intersections
c
c      write(6,'(a,2(1x,f12.6),i6)') 'PSIN INTERSECTIONS:',
c     >        rstart,zstart,n_int
c      write(6,'(a,6(1x,f12.6))') 'R   :',(r_int(in),in=1,n_int)
c      write(6,'(a,6(1x,f12.6))') 'Z   :',(Z_int(in),in=1,n_int)
c      write(6,'(a,6(1x,f12.6))') 'PSI :',(psi_int(in),in=1,n_int)
c      write(6,'(a,6(1x,f12.6))') 'DIST:',(dist_int(in),in=1,n_int)
c  
      return
      end
c
c
c
      subroutine read_psin_data
      implicit none
      include 'params'
      include 'psin_data'
c
c     READ_PSIN_DATA: This routine reads an input file
c                        of column ordered data - R,Z,PSI
c                        representing the wall structure and 
c                        PSI data at each turning point. It 
c                        additionally calculates the length
c                        of each segment.
c
      integer infile
      parameter (infile = 98) 
c
      logical eof
      real*8 r,z,psi,psimin
      integer ierr,in,next_in,id,len,lenstr,wall_num
      external lenstr
c
c     Initialization
c
      external_psin = .false. 
      psimin = HI 
      rpsimin = 0.0 
c
      len = lenstr(psin_filename)
      if (len.le.0) then 
         write(6,'(a)') 'ERROR: EXTERNAL_PSIN: FILE NAME IS EMPTY' 
         write(0,'(a)') 'ERROR: EXTERNAL_PSIN: FILE NAME IS EMPTY' 
         return
      endif 
c
      open(UNIT=infile,FILE=psin_filename(1:len),STATUS='OLD',
     >     IOSTAT=ierr,ERR=100)
c
      eof = .false.
      n_elements = 0
c
      do while (.not.eof)
c
         read(infile,*,IOSTAT=ierr) wall_num,r,z,psi
c
         if (ierr.eq.0) then
c   
            n_elements= n_elements+1
            surface_data(n_elements,1) = r  
            surface_data(n_elements,2) = z  
            surface_data(n_elements,3) = psi  
c
         else
c
            eof = .true. 
c
         endif              

      end do 
c
c     Calculate lengths and find R value of Minimum PSIN
c
      if (n_elements.eq.0) return
c
      do in = 1,n_elements
c  
         if (surface_data(in,3).lt.psimin) then 
            psimin = surface_data(in,3)
            rpsimin= surface_data(in,1)
         endif
c
         if (in.eq.n_elements) then 
            next_in = 1
         else
            next_in = in+1
         endif
c
         surface_data(in,4) = sqrt(
     >      (surface_data(next_in,1)-surface_data(in,1))**2
     >     +(surface_data(next_in,2)-surface_data(in,2))**2)
c     
c
c         write(6,'(a,4(1x,f12.6))') 'TEST:',
c     >                  (surface_data(in,id),id=1,4)
c
      end do
c
      write(6,'(a,2g12.5)') 'PSIMIN:',psimin,rpsimin

c
      external_psin = .true.
c

      close(infile)
c
      return

100   continue

      len = lenstr(psin_filename)
      write(0,'(a)') 'ERROR: PSIN DATA FILE DOES NOT EXIST: '
     >                        //psin_filename(1:len)
      write(6,'(a)') 'ERROR: PSIN DATA FILE DOES NOT EXIST: '
     >                        //psin_filename(1:len)

      close(infile)

      return
      end 
c
c
c
      subroutine get_psin_filename(ierr)
      implicit none
      include 'params'
      include 'psin_data'
      integer ierr 

      call RDC(psin_filename,'PSIN FILENAME',ierr)
      if (ierr.ne.0) return
      
      return 
      end
c
c
c
      subroutine plot_los_multicase(iselect,istate,npts,nlines,
     >                  iexpt,iaxis,iavg,ifact,optval,
     >                  iopt,job,title,table,avs,navs,nplots,
     >                  iplot,nizs,ierr) 
      implicit none
c
      include 'params'
      include 'psin_data' 
      include 'adas_data_spec'
c
      integer iselect,istate,npts,nlines,iexpt,iaxis,ifact,ierr 
      integer iavg
      real optval 
      character*(*) job,title,table
      integer navs,iopt,nplots,iplot,nizs
      real avs(0:navs)
c
c     PLOT_LOS_MULTICASE: This routine plots a generalized LOS plot. Each LOS
c               plot may have an individual R,Z,Theta,dTheta and
c               length/distance over which the LOS will be integrated.
c
c               This code calculates the LOS integrals over multiple DIVIMP
c               solutions and then plots all of them on the final result. 
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
      real psin  
      integer targ,itarg,icnt
c
      real r,z,thet,dthet,dist,losval,losval2,losvalexpt,scalef
      real g_d,g_w,g_l
c
      real tmpplot(maxnks,maxnrs),tmpplot2(maxnks,maxnrs)
      real mfact
c
c     Plot labels   
c
      character*36 blabs(2),xlab,ylab,ylab2,blabs2(2)
      character*44 ref,plane,anly,nview,ref2
      character*32 datatitle
      character*36 custom_xaxis_label 
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
      real div2D_array(maxix,maxiy)
      real expt_rmin,expt_rmax,expt_zmin,expt_zmax,expt_dr,expt_dz
c
      integer iexpt_out
      integer itype

c
c     Variables for multicase output 
c
      
      character*256 basename,base_psin_filename,cmnd,
     >              name_format,psin_format,case_name
      integer start_index,name_increment,ncases,ic,case_index,nsize
      integer ntarg1,ntarg2
      real targ1outs(maxthe),targ1vals(maxthe,maxngs),targ1wids(maxthe)
      real targ1chan(maxthe)
      real targ2outs(maxthe),targ2vals(maxthe,maxngs),targ2wids(maxthe)
      real targ2chan(maxthe)
      integer ttarg(maxthe)
c
c     ADAS variables
c
      CHARACTER ADASID*80,graph3*80
      CHARACTER XFESYM*2
      character adasex*3
      integer   adasyr
      INTEGER ISELE,ISELR,ISELX,iseld,ircode
c
      write(6,'(a,i6)') 'OUTLOS:IAXIS=',iaxis
      write(6,'(a,2i6)') 'OUTLOS:STATE=',istate,iselect
c
c     Initialization
c
      external_psin = .false.
      targ      = 0
      itarg     = 0    
      iexpt_out = 0
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
c     Load ADAS data if it is required (iselect = 4,5)
c
      if (iselect.eq.4.or.iselect.eq.5.or.iselect.eq.22) then  
         CALL RDG1 (GRAPH3,ADASID,adasyr,adasex,
     >              ISELE,ISELR,ISELX,ISELD,IERR)
c
         if (ierr.ne.0) return 
c
c        Assign ADAS data to global variables for reuse
c
         cadasid = adasid
         cadasyr = adasyr
         cadasex = adasex
         cisele  = isele
         ciselr  = iselr
         ciselx  = iselx
         ciseld  = iseld
c
         write(6,'(a)') 'OUT_LOS_MC:ADAS data assigned' 
c
      endif
c
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
      if (iaxis.eq.10) then 
c
         call get_psin_filename(ierr)
         if (ierr.ne.0) return 
c
      endif
c
c     For iaxis = 11 - load multigrid data for both PSIN and grid names 
c
      if (iaxis.eq.11) then  
c
         call get_psin_filename(ierr) 
         if (ierr.ne.0) return
c
         call rdfn_multicase(graph,cmnd,basename,start_index,
     >                           name_increment,ncases,ierr)
         if (ierr.ne.0) return
c
      endif
c
c
c     Load additional data required for geometry scaling option 5
c
c     These inputs describe the geometry of the DIIID bolometers and 
c     allow for the calculation of additional scaling factors. 
c
c     When these values are entered the array of dT values is overwritten
c     with values based on these quanitities.  
c
      if (ifact.eq.5.or.ifact.eq.6) then  
         call RDG_REAL_ARRAY(GRAPH,geofact_d,maxip,ndg_d,ierr)
         if (ierr.ne.0) return
         call RDG_REAL_ARRAY(GRAPH,geofact_w,maxip,ndg_w,ierr)
         if (ierr.ne.0) return
         call RDG_REAL_ARRAY(GRAPH,geofact_l,maxip,ndg_l,ierr)
         if (ierr.ne.0) return
c
c        Recalculate the dtheta values specified based on the D,W,L values input
c
         ndt = npts
c
         write(6,'(a,i5)') 'OUTLOS: Recalculating DTHETA',ndt
c
         do in = 1,npts
c
            ii = min(in,ndg_d)
            g_d = geofact_d(ii)
c
            ii = min(in,ndg_w)
            g_w = geofact_w(ii)
c
            ii = min(in,ndg_l)
            g_l = geofact_l(ii)
c             
            dtheta(in) = 2.0 * atan ((g_d+g_w)/(2.0*g_l)) * raddeg
c
            write(6,'(a,i5,5(1x,g12.5))') 'DT:',in,dtheta(in),
     >                    g_d,g_w,g_l 
c
         end do
c
      endif

c
c     All data to calculate LOS plot is now available - loop through 
c     calculating values. 
c
      IF (ifact.LT.0 .OR. ifact.GT.10) ifact = 0
c
c     Set LOS scaling factor and note value for plot
c
      IF (ifact.EQ.0) THEN
         mfact = 1.0
         WRITE(ANLY,'(''NO GEOMETRY FACTOR APPLIED'')')
      ELSEIF (ifact.EQ.1) THEN
         WRITE(ANLY,'(''GEOMETRY FACTOR = DTHETA / (2*PI)'')')
      ELSEIF (ifact.EQ.2) THEN
         MFACT = 1.0 / (2.0 * PI)
         WRITE(ANLY,'(''GEOMETRY FACTOR = 1 / (2*PI)'') ')
      ELSEIF (ifact.EQ.3) THEN
         MFACT = 1.0 / ( 4.0 * PI)
         WRITE(ANLY,'(''GEOMETRY FACTOR = 1 / (4*PI)'')')
      ELSEIF (ifact.EQ.4) THEN
         WRITE(ANLY,'(''GEOMETRY FACTOR=TAN(DTHETA/2)^2 / PI'')')
      ELSEIF (ifact.EQ.5) THEN
         WRITE(ANLY,'(''2D GEOMETRY FACTORS FROM W,D,L'')')
      ELSEIF (ifact.EQ.6) THEN
         WRITE(ANLY,'(''GEOMETRY FACTORS FROM W,D,L'')')
      ELSEIF (ifact.EQ.7) THEN
         WRITE(ANLY,'(''MODIFIED GEOMETRY FACTORS FROM W,D,L'')')
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
c 
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
      endif     
c
c     If the scale factor has been set to anything other than 1.0 - then 
c     adjust the BLAB value. 
c
      if (nds.ne.1.or.(nds.eq.1.and.scalefs(1).ne.1.0)) then
         call set_blab(iselect,istate,2,nizs,blabs(1))
      endif
c
c     Setup for loop over cases
c
      ntarg1 = 0
      ntarg2 = 0
      base_psin_filename= psin_filename
c
c     Loop over the entire set of cases and save all results ...
c
      do ic = 1,ncases
c
c        Calculate case name and psin file name
c
         case_index = start_index + (ic-1) * name_increment

         if (case_index.gt.0) then  
            nsize = int(log10(real(case_index))) + 1
         endif
c
         write(name_format,'(''(A,I'',i1,'')'')')  nsize          
         write(psin_format,'(''(A,I'',i1,'',A)'')')  nsize          
c
         len = lenstr(basename)
         write(case_name,name_format) basename(1:len),case_index
     >               ,base_psin_filename
         len = lenstr(base_psin_filename)
         write(psin_filename,psin_format) base_psin_filename(1:len),
     >                                    case_index,'.dat'
c
c        Load case data into OUT arrays - base data is saved if this is first call
c
         len=lenstr(case_name)
         call loaddata(case_name(1:len),cmnd,-1,-1)
c
c        Load PSIN reference data
c
         if (iaxis.eq.11) then 
            call read_psin_data
         endif
c
c        Load data to integrate from the new solution 
c
c        Based on the values of iselect and istate - load up the required 
c        2D array for the LOS plot 
c
c        Scaling by absfac (if required) is performed in the 
c        load_divdata_array routine. 
c
         call rzero(tmpplot,maxnks*maxnrs)
c
         itype = 1  
c
c        For Steradian weighted plots - load a different label - assuming LOS here. 
c
         if (ifact.eq.3) itype = 2
c
c        Set ADAS data reuse flag  
c
         if (iselect.eq.4.or.iselect.eq.5.or.iselect.eq.22)
     >       cadas_switch = 1
c
c        Load array for integration 
c
         call load_divdata_array(tmpplot,iselect,istate,itype,
     >                         ylab,blabs(1),ref,nizs,ierr)
c
c        Load denomintor if averaged temperatures are being plotted.  
c
         if (iselect.eq.22) then 
c
            call rzero(tmpplot2,maxnks*maxnrs)
c
            call load_divdata_array(tmpplot2,4,istate,itype,
     >                      ylab2,blabs2(1),ref2,nizs,ierr)
c
         endif
c
c        Reset ADAS data reuse flag  
c
         if (iselect.eq.4.or.iselect.eq.5.or.iselect.eq.22) 
     >      cadas_switch = 0
c
c        Check return code from data load and exit if non-zero   
c
         
         if (ierr.ne.0) then 
             write(6,'(a,i6)') 'OUT_LOS_MC: Non-zero RC'//
     >                       ' from Load_divdata_array routine =',ierr
             return 
         endif
c
c
c        Check to see if a LOS integration of a 2D experimental data array is 
c        also requested.
c        
c        A value of IEXPT of -2 means to also convert the DIVIMP data 
c        to the same format as the experimental data prior to the 
c        calculation of the values to be plotted.
c        
         if (iexpt.eq.-1.or.iexpt.eq.-2) then  
c        
c           Set to 2 plots on page
c        
            ngs =2  
c        
            blabs(2) = 'EXPT EXPT'
c        
            call rdfn(graph,filename,ifnopt,ierr)
c        
            if (ierr.ne.0) then  
c        
               len = lenstr(filename)
               write(0,*) 'PLOT_LOS:ERROR READING 2D EXPT'
     >                      //' DATAFILE INPUT LINE: '
     >                      //filename(1:len),ifnopt,ierr
               iexpt = 0
               ngs = 1
c        
            else 
c        
               call load_2Dexpt_data(filename,datatitle,
     >                     maxix,maxiy,nix,niy,expt_array,raxis,zaxis,
     >                     expt_rmin,expt_rmax,expt_zmin,expt_zmax,
     >                     expt_dr,expt_dz,ierr)
c        
               if (ierr.ne.0) then  
c        
                  len = lenstr(filename)
                  write(0,*) 'PLOT_LOS:ERROR LOADING 2D EXPT DATAFILE: '
     >                      //filename(1:len)
                  iexpt = 0
                  ngs = 1
c        
               else
         
                  len=lenstr(datatitle)
                  blabs(2) = 'EXPT '//datatitle(1:len)
c        
               endif 
c        
            endif
c        
c           Convert the DIVIMP array to an experimental data format
c        
            if (iexpt.eq.-2) then 
c        
               call create_2Ddata(tmpplot,maxix,maxiy,nix,niy,
     >                     div2D_array,raxis,zaxis,0)
c        
            endif
c        
         endif
c        
c
C        INTEGRATE
C        
c        Perform a series of single LOS integrations  
c        
c        Use the global npts that was specified in the first line - if
c        there are less that npts values specified for any line - the 
c        last value on each input line will be used for all remaining
c        lines of sight. 
c            
         do in = 1,npts 
c        
c           Calculate the R, Z, Theta and Dtheta values for each 
c           LOS.
c                      
c           R-obs 
c        
            ii = min(in,nr)
            r = rvals(ii)  
c        
c           Z-obs
c        
            ii = min(in,nz)
            z = zvals(ii)  
c        
c           Theta
c        
            ii = min(in,nt)
            thet = theta(ii)  
c        
c           DTheta
c        
            ii = min(in,ndt)
            dthet = dtheta(ii)  
c        
c           Dist
c        
            ii = min(in,nd)
            dist = dists(ii)  
c        
c           Additional Scaling factor
c        
            ii = min(in,nds)
            scalef = scalefs(ii)  
c        
c           Additional Geomtery Factors
c        
            ii = min(in,ndg_d)
            g_d = geofact_d(ii)  
         
            ii = min(in,ndg_w)
            g_w = geofact_w(ii)  
         
            ii = min(in,ndg_l)
            g_l = geofact_l(ii)  
c        
c        
c           NOTE: IF the LOS data specified contains only ONE value for 
c                 Theta, ONE value for DTheta, ONE value for Robs and ONE
c                 value for Zobs and yet npts is greater than ONE - then 
c                 the code assumes that instead of a series of degenerate/
c                 identical views - what is actually wanted is a series of
c                 views starting at THETA with a separation of Dtheta and 
c                 a width of DIST.
c        
            if (npts.gt.1.and.
     >          nr.eq.1.and.nz.eq.1.and.nt.eq.1.and.ndt.eq.1) then 
c        
                thet  = thet + (npts-1) * dthet
                dthet = dist
                dist  = 0.0
c        
            endif            
c        
c           Set scale factor based on viewing width 
c        
            if (ifact.eq.1) then 
c        
c               Convert Dthet to radians for this scaling
c        
                MFACT = (abs(DTHET)*degrad) / (2.0 * PI)
c        
            elseif (ifact.eq.4) then 
c        
c               Convert Dthet to radians for this scaling
c        
                MFACT = tan((abs(DTHET)/2.0)*DEGRAD)**2 / PI
c        
c           INCORRECT SCALING - JUST FOR REFERENCE
c        
            elseif (ifact.eq.5) then 
c        
c               This is the scaling factor for the second dimension of the cross-section
c        
                MFACT = tan(abs(DTHET)*DEGRAD) / PI
c        
            elseif (ifact.eq.6) then 
c        
c               Convert Dthet to radians for this scaling
c        
                MFACT = (abs(DTHET)*degrad) / (4.0 * PI)
c        
            elseif (ifact.eq.7) then 
c        
c               Convert Dthet to radians for this scaling
c        
                MFACT = 1.0 / (4.0 * PI)
c        
            endif  
c        
c            write (6,'(a,i4,5g12.4)') 'LOS:',in,r,z,thet,dthet,dist
c            write (0,'(a,i4,5g12.4)') 'LOS:',in,r,z,thet,dthet,dist
c             
c           If using DIVIMP data converted to experimental data format
c        
            if (iexpt.eq.-2) then 
c        
               CALL LOSINTEXPT(losvalexpt,thet,dthet,1,
     >                  R,Z,nlines,dist,
     >                  maxix,maxiy,nix,niy,div2D_array,raxis,zaxis)
c        
               tvals(in,1) = losvalexpt * mfact * scalef
c        
c        
c           DIVIMP data in regular or DIVIMP internal format 
c        
            else
         
               write (6,'(a,i6)') 'LOS: ifact = ',ifact
         
               if (ifact.eq.5.or.ifact.eq.6.or.ifact.eq.7) then
         
                  CALL LOSINT_SCALE(losval,thet,dthet,1,
     >                     R,Z,nlines,tmpplot,dist,g_d,g_w,g_l,ifact)
         
               else
         
                  CALL LOSINT(losval,thet,dthet,1,
     >                     R,Z,nlines,tmpplot,dist,iavg)
         
c
c                 Load denominator for line averaged quantities
c
                  if (iselect.eq.22) then    

                     CALL LOSINT(losval2,thet,dthet,1,
     >                     R,Z,nlines,tmpplot2,dist,iavg)

c
c                    Scale losval2 appropriately
c
                     losval2 = losval2 * mfact * scalef  
c
                  endif

               endif
c
c              Assign values to plot array
c
               tvals(in,1) = losval * mfact * scalef
c
               if (iselect.eq.22.and.losval2.ne.0.0) then

                  tvals(in,1) = tvals(in,1) / losval2 

               endif  
c        
            endif
c        
c           If LOS over 2D experimental data is required
c        
            if (iexpt.eq.-1.or.iexpt.eq.-2) then 
c        
               CALL LOSINTEXPT(losvalexpt,thet,dthet,1,
     >                  R,Z,nlines,dist,
     >                  maxix,maxiy,nix,niy,expt_array,raxis,zaxis)
c        
               tvals(in,2) = losvalexpt * mfact * scalef
c        
            endif
c        
c           Calculate axis values
c        
c           Theta 
c        
            if (iaxis.eq.1) then  
c        
               touts(in) = thet
               twids(in) = abs(dthet)
c        
c           R-intersection 
c        
            elseif (iaxis.eq.2) then 
c        
               call adjustout(thet,1,optval,r,z)
c        
               touts(in) = thet
               twids(in) = 1.0 
c        
c           Z-intersection - needs testing
c        
            elseif (iaxis.eq.3) then 
c        
               call adjustoutz(thet,1,optval,r,z)
c        
               touts(in) = thet
               twids(in) = 1.0 
c        
c           R-observation 
c        
            elseif (iaxis.eq.4) then 
c        
               touts(in) = r
               twids(in) = 1.0 
c        
c           Z-observation 
c        
            elseif (iaxis.eq.5) then 
c        
               touts(in) = z
               twids(in) = 1.0 
c        
c           Index number
c        
            elseif (iaxis.eq.6) then 
c        
               touts(in) = in
               twids(in) = 1.0 
c        
            elseif (iaxis.eq.7.or.iaxis.eq.10) then 
c        
               call calcpsin(r,z,thet,psin,targ)
c        
               if (targ.ne.0.and.itarg.eq.0) then
                  itarg = targ
               endif 
c        
c              ONLY LOS that strike one target are desired - the code
c              assumes that the first LOS that intersects a target defines
c              the target of interest. All other points will be excluded. 
c              LOS that do strike the target are assigned a proper psin value
c              - all others will get a psin value of -100.0 - all these 
c              points will be filtered out before plotting.
c        
c        
               if (itarg.ne.0.and.targ.eq.itarg) then  
c        
                  touts(in) = psin
                  twids(in) = 1.0 
c        
               else
c        
                  touts(in) = -100.0
                  twids(in) =  1.0 
c        
               endif      
         
               write(6,'(a,2i6,f12.6)') 'ITARG:',targ,itarg,psin
         
c        
c           Index number minus 1 - corresponds to channel number 0+
c        
            elseif (iaxis.eq.8) then 
c        
               touts(in) = in + optval
               twids(in) = 1.0 
c        
c           Custom/user specified axis - must be ordered - ascending or descending - not mixed
c        
            elseif (iaxis.eq.9) then  
c        
c              There should be an axis value specified for every point - but protect
c              against an error just in case.  
c        
               ii = min(in,nca)
               touts(in) = custom_axis(ii)
c        
            elseif (iaxis.eq.11) then 
c        
               call calcpsin(r,z,thet,psin,targ)
c        
               touts(in) = psin
               twids(in) = 1.0
               ttarg(in) = targ
c        
            endif         
c        
c           End of LOS calculation loop 
c        
         end do
c
c        Save calculated LOS integration results including PSIN and target data        
c
         do in = 1,npts
c
c           Do not include points not assigned to a target
c           Accumulate data for each target separately
c
            if (ttarg(in).eq.1) then 
c
c              Target 1 - OUTER for Xpoint down
c
               ntarg1 = ntarg1+1
               targ1outs(ntarg1) = touts(in)   
               targ1vals(ntarg1,1) = tvals(in,1)
               targ1wids(ntarg1) = twids(in)
               targ1chan(ntarg1) = in
c
            elseif (ttarg(in).eq.2) then 
c
c              Target 2 - INNER for Xpoint down
c
               ntarg2 = ntarg2+1
               targ2outs(ntarg2) = touts(in)   
               targ2vals(ntarg2,1) = tvals(in,1)
               targ2wids(ntarg2) = twids(in)
               targ2chan(ntarg2) = in
c
            endif 
c
         end do

c
c     End of case iteration loop
c
      end do 

c
c     Reload base case data in case of other plots
c
      call loaddata(' ',' ',-1,-2)
c
c     If plotting versus the psin scale - filter out all points with a 
c     PSIN value that is set to -100.0 - these points either do not 
c     intersect a target or are not on the default target defined by the 
c     first LOS to intersect a target. 
c
      if (iaxis.eq.7.or.iaxis.eq.10) then 
c
         icnt = 0
c
         do in = 1,npts
c
            if (touts(in).le.-100.0) then
c
               icnt = icnt + 1
c
            else
c              
               touts(in - icnt) = touts(in)
               twids(in - icnt) = touts(in)
               tvals(in - icnt,1) = tvals(in,1)
c
            endif
c
         end do 
c
         write(6,'(a,3i6)') 'OUTLOS:PSI REMOVED:',npts,icnt,npts-icnt 
c
         npts = npts - icnt
c
c
      endif 
c
c     Reorder results so that axis is in ascending order
c
c      if (touts(1).gt.touts(npts)) then
c
c        reverse order of touts and tvals entries.
c        
c         do in = 1,npts/2
c
c           axis
c
c            tmp_store = touts(in)
c            touts(in) = touts(npts-in+1)
c            touts(npts-in+1) = tmp_store
c
c           Divimp values
c
c            tmp_store = tvals(in,1)
c            tvals(in,1) = tvals(npts-in+1,1)
c            tvals(npts-in+1,1) = tmp_store
c
c           Experimental values
c
c            tmp_store = tvals(in,2)
c            tvals(in,2) = tvals(npts-in+1,2)
c            tvals(npts-in+1,2) = tmp_store
c
c         end do   
c
c      endif


c
c     Load value into PLANE if req'd 
c
      if (nds.eq.1.and.scalefs(1).ne.1.0) then 
c
          write(PLANE,'(a,g12.5)')
     >   'SCALING FACTOR = ',scalefs(1)
c
      elseif (nds.gt.1) then  
c
          write(PLANE,'(a)')
     >   'ADDITIONAL SCALING FACTORS APPLIED'
c
      elseif (nd.eq.1.and.dists(1).gt.0.0) then
c
         write(PLANE,'(a,f6.3,a)')
     >   'INTEGRATION OVER ',dists(1),' (M) ONLY'
c
      endif   

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
          call calc_expt(iexpt,touts,tvals,maxthe,npts,
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
c     PLOT TARGET 1 DATA - OUTER TARGET 
c

      PLANE = 'OUTER TARGET'

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
c     Fiddle with expt_nsets data - first half to first target 
c
      call exchange_expt_datasets(1) 
c
c     Find max and min axis values
c
      themin = HI
      themax = -HI
c
      do in = 1,ntarg1
         themin = min(targ1outs(in),themin) 
         themax = max(targ1outs(in),themax) 
      end do
c
      if (ntarg1.gt.1) then 
c
c
          ngs = 1
c
          do in = 1,ntarg1
             write(6,'(a,i4,3(1x,g12.5))') 'TARG1VALS:OUTER',in,
     >                              targ1outs(in),
     >                              targ1vals(in,1),targ1chan(in)

          end do
c
          itec = 1
          ismoth = 99          
c
          CALL DRAW(TARG1OUTS,TARG1WIDS,TARG1VALS,MAXTHE,NTARG1,
     >              ANLY,NGS,
     >              ISMOTH,THEMIN,THEMAX,-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >              JOB,TITLE,XLAB,YLAB,BLABS,REF,NVIEW,PLANE,
     >              TABLE,IOPT,2,1.0,iexpt_out)
c
      endif 





c
c     PLOT TARGET 2 DATA - INNER TARGET 
c
      PLANE = 'INNER TARGET'
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
c
c     Fiddle with expt_nsets data - second half to second target 
c
      call exchange_expt_datasets(2) 
c
c     Find max and min axis values
c
      themin = HI
      themax = -HI
c
      do in = 1,ntarg2
         themin = min(targ2outs(in),themin) 
         themax = max(targ2outs(in),themax) 
      end do
c
      if (ntarg2.gt.1) then  

c
          ngs = 1
c
          do in = 1,ntarg2
             write(6,'(a,i4,3(1x,g12.5))') 'TARG2VALS:INNER',in,
     >                            targ2outs(in),
     >                            targ2vals(in,1),targ2chan(in)
          end do
c
          itec = 1
          ismoth = 99          
c
          CALL DRAW(TARG2OUTS,TARG2WIDS,TARG2VALS,MAXTHE,NTARG2,
     >              ANLY,NGS,
     >              ISMOTH,THEMIN,THEMAX,-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >              JOB,TITLE,XLAB,YLAB,BLABS,REF,NVIEW,PLANE,
     >              TABLE,IOPT,2,1.0,iexpt_out)
c
      endif 
c
c     Restore original experimental datasets array
c
      call exchange_expt_datasets(0) 
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
      subroutine exchange_expt_datasets(opt)
      implicit none
      integer opt
c
      include 'params'
      include 'outcom'
c
c     EXCHANGE_EXPT_DATASETS:
c
c     This routine saves the data in the experimantal dataset
c     inputs and then rearraneges it so that only some of the 
c     data is stored. In this way the input can be split 
c     to describe the experimental data for multiple plots within
c     a single plot designation. 
c
      logical saved
      data saved /.false./
      integer in
      integer tmp_expt_nsets
      integer tmp_expt_datasets(max_expt_datasets)
c
c     Restore data if SAVED is available 
c
      if (opt.eq.0) then  
         if (saved) then 
c
            expt_nsets = tmp_expt_nsets
c
            do in = 1,tmp_expt_nsets
               expt_datasets(in) = tmp_expt_datasets(in) 
            end do
c
            saved = .false. 
c
         endif

c 
c     Store First half of numbers - save original if not saved
c 
      elseif (opt.eq.1) then  
c
         if (.not.saved) then 
c
            tmp_expt_nsets = expt_nsets
c
            do in = 1,expt_nsets
               tmp_expt_datasets(in) = expt_datasets(in) 
            end do
c
            saved = .true.
c   
         endif
c
c        Copy first half
c
         expt_nsets = tmp_expt_nsets/2
c 
         do in = 1,tmp_expt_nsets/2
            expt_datasets(in) = tmp_expt_datasets(in)
         end do

c
c     Store Second half of numbers - save original if not saved
c
      elseif (opt.eq.2) then  
c
         if (.not.saved) then 
c
            tmp_expt_nsets = expt_nsets
c
            do in = 1,expt_nsets
               tmp_expt_datasets(in) = expt_datasets(in) 
            end do
c
            saved = .true.
c   
         endif
c
c        Copy second half
c
         expt_nsets = tmp_expt_nsets - tmp_expt_nsets/2
c 
         do in = tmp_expt_nsets/2+1, tmp_expt_nsets
            expt_datasets(in-tmp_expt_nsets/2) 
     >                       = tmp_expt_datasets(in)
         end do
c
      endif 

      return 
      end

      subroutine map_file_data_to_psin(r,z)
      implicit none
      include 'params'
      real r,z
      integer in,targ
      integer unit
      real thet,psin

      call find_free_unit_number(unit)

      write(9,'(a,2g12.5)') 'R,Z:',r,z

      open(unit,file="../results/input_data.txt")

 10   read(unit,*,end=20,err=20) thet

      call calcpsin(r,z,thet,psin,targ)

      write(9,'(a,2i6,4(1x,g20.12))') 
     >            'MAPPING THETA TO PSIN:',in,targ,
     >            r,z,thet,psin
      
      goto 10

 20   close(unit)

      return
      end
