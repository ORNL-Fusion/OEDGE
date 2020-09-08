      subroutine plot_xsection(r1p,z1p,r2p,z2p,npts,iselect,istate,
     >                       iexpt,iavg,iopt,
     >                       nizs,job,title,table,avs,navs,iplot,
     >                       nplots)
c
      use mod_params
      use mod_cgeom
      use mod_dynam3
      use mod_comtor
      use mod_pindata
      implicit none
c
c     include 'params' 
c     include 'cgeom'
c     include 'dynam3'
c     include 'comtor'
c     include 'pindata'
c
      real r1p,z1p,r2p,z2p 
      integer iopt,iplot,nplots,npts,iselect,istate,nizs,iexpt,iavg
      character*(*) job,title,table
      integer navs
      real avs(0:navs)
c
c     PLOT_XSECTION: This routine generates a cross-section across a 
c     specified 2D contour plot. The ISELECT input specifies which 
c     quantity is to be plotted. 
c
c
c     Local Variables 
c
      real tstart,tend,tstep
      integer ptype,ik,ir,iz
c
      real tmpplot(maxnks,maxnrs)
c
c     Plot labels   
c
      character*36 blabs(2),xlab,ylab
      character*44 ref,plane,anly,nview
      character*32 datatitle
c
      INTEGER IGNORS(MAXNGS),ngs
c
      REAL TOUTS(MAXTHE),TWIDS(MAXTHE),TVALS(MAXTHE,MAXNGS)
      real tmptvals(0:maxthe-1)
      REAL THEMIN,THEMAX
c
      real getvalrz,r,z
c
      external getvalrz
c
      integer len,lenstr,ii,ierr
c
c     Variables for 2D experimental data
c
      character*256 graph
      character*256 filename
      integer  maxix,maxiy,nix,niy,ifnopt
      parameter(maxix=100,maxiy=100)
      real expt_array(maxix,maxiy),raxis(maxix),zaxis(maxiy)
      real div2D_array(maxix,maxiy)
      real expt_rmin,expt_rmax,expt_zmin,expt_zmax,expt_dr,expt_dz
      real getvalexpt2D
      external getvalexpt2D
c
c     Variables for averaging
c
      real r2a,z2a,r1a,z1a,ra,za
      integer in,avgnum,icnt
c
c     Initialization
c
      do ii = 1,maxngs 
         ignors(ii) = 1
      end do
c
      ngs =1 
c
c     Set NVIEW
c
      NVIEW  = 'CROSS-SECTION OF 2D CONTOUR'
      REF = ' ' 
      YLAB = ' '
c
c     Determine plot axis
c     
c     If R1 = R2 - plot along Z
c     
      if (r1p.eq.r2p) then 
c     
         if (z1p.lt.z2p) then 
            tstart = z1p
            tend   = z2p
         else
            tstart = z2p
            tend   = z1p
         endif
         ptype = 1

         XLAB = 'Z (M)'
c     
c     If Z1 = Z2 - plot along R
c     
      elseif (z1p.eq.z2p) then            
c     
         if (r1p.lt.r2p) then 
            tstart = r1p
            tend   = r2p
         else
            tstart = r2p
            tend   = r1p
         endif
         ptype =2  
c
         XLAB = 'R (M)'
c     
c     Arbitrary cut   
c     
      else 
c     
         ptype = 3
c     
         tstart = 0.0
c     
         tend = sqrt((r2p-r1p)**2+(z2p-z1p)**2) 
c         
         XLAB = 'DIST (M)'
c
      endif 
c
c     Write the start and end points to reference strings 
c
      write(anly,'(a,''('',f8.5,'','',f8.5,'')'')') 'START: ',r1p,z1p
      write(plane,'(a,''('',f8.5,'','',f8.5,'')'')') 'END  : ',r2p,z2p
c
c
c     Fill the tmpplot array with the data specified by ISELECT
c
      call rzero(tmpplot,maxnks*maxnrs)
c
      call load_divdata_array(tmpplot,iselect,istate,0,
     >                           ylab,blabs(1),ref,nizs,ierr)
c
c     Check return code from data load and exit if non-zero   
c
      if (ierr.ne.0) return 
c
      CALL PRRMATDIV(tmpplot,MAXNKS,nks(irsep),NRS,6,
     >                                  '341 - TMPPLOT')
c
c     Load and calculate any 2D experimental data that is also needed
c     for the plot. IEXPT = -1 means to load a regular 2D array of
c     experimental data and take a cross-section at the same points as
c     was used for the DIVIMP data. 
c
c     A value of IEXPT of -2 means to also convert the DIVIMP data 
c     to the same scaling as the experimental data prior to the 
c     calculation of the values to be plotted.
c
      if (iexpt.eq.-1.or.iexpt.eq.-2) then  
c
c        Set to 2 plots on page
c 
         ngs =2  
         blabs(2) = 'EXPT EXPT DATA'
c
         call rdfn(graph,filename,ifnopt,ierr)
c
         if (ierr.ne.0) then  
c
            len = lenstr(filename)
            write(0,*) 'PLOT_XSECTION:ERROR READING 2D EXPT DATAFILE'
     >                   //' INPUT LINE: '
     >                   //filename(1:len)
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
               write(0,*) 'PLOT_XSECTION:ERROR LOADING'
     >                   //' 2D EXPT DATAFILE: '
     >                   //filename(1:len)
               iexpt = 0
               ngs = 1
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
c     Calculate the axis and the values to be plotted
c     
      tstep = (tend-tstart)/(npts-1) 
c
      do ii = 1,npts
c     
         touts(ii) = tstart + (ii-1)* tstep          
         TWIDS(II) = 1.0
c
         if (ptype.eq.1) then 
c
            r = r1p 
            z = touts(ii) 
c
         elseif (ptype.eq.2) then  
c 
            r = touts(ii)
            z = z1p
c
         elseif (ptype.eq.3) then 
c
            r = r1p + (ii-1) * (r2p-r1p)/(npts-1)
            z = z1p + (ii-1) * (z2p-z1p)/(npts-1)
c
         endif
c
c        Calculate DIVIMP value based on DIVIMP Data converted to the 
c        experimental format
c
         if (iexpt.eq.-2) then   
c
            tvals(ii,1) = getvalexpt2D(r,z,
     >                  maxix,maxiy,nix,niy,div2D_array,raxis,zaxis,
     >                  expt_rmin,expt_rmax,expt_zmin,expt_zmax,
     >                  expt_dr,expt_dz,ifnopt)
c
         else
c
c           Point is average of avgnum subpoints
c
            if (iavg.lt.0) then 
c
               avgnum = abs(iavg)
c
               r2a = r1p + (ii) * (r2p-r1p)/(npts-1)   
               z2a = z1p + (ii) * (z2p-z1p)/(npts-1)   
               r1a = r1p + (ii-2) * (r2p-r1p)/(npts-1)   
               z1a = z1p + (ii-2) * (z2p-z1p)/(npts-1)   
c
               tvals(ii,1) = 0.0
c
               do in = 1,avgnum
c
                  ra = r1a + in * (r2a-r1a)/(avgnum+1)
                  za = z1a + in * (z2a-z1a)/(avgnum+1)
c
                  tvals(ii,1) = tvals(ii,1) 
     >                      + getvalrz(tmpplot,ra,za)
c
               end do
c
c              Perform average
c
               tvals(ii,1) = tvals(ii,1) / avgnum 
c
c           No averaging 
c            
            else
c
               tvals(ii,1) = getvalrz(tmpplot,r,z)
c
            endif 
c
         endif
c
         if (iexpt.eq.-1.or.iexpt.eq.-2) then
c
            tvals(ii,2) = getvalexpt2D(r,z,
     >                  maxix,maxiy,nix,niy,expt_array,raxis,zaxis,
     >                  expt_rmin,expt_rmax,expt_zmin,expt_zmax,
     >                  expt_dr,expt_dz,ifnopt)
c
         endif
c
         write (6,'(a,i4,2f10.5,1x,g13.6)') 'PLOT:',ii,r,z,
     >                 tvals(ii,1),tvals(ii,2)
c
      end do
c
c     Perform global average if one is defined. Use Tvals(ii,maxngs) for
c     temporary storage. 
c
      if (iavg.gt.0) then 
c
         do ii = 1,npts
c
            icnt = 0 
c
            do in = -(iavg/2),(iavg/2)
c
               if ((ii+in).ge.1.and.(ii+in).le.npts) then
c
                  icnt = icnt + 1 
                  tvals(ii,maxngs) = tvals(ii,maxngs) + tvals(ii+in,1)
c
               endif
c
            end do  
c
c           Renormalize 
c
            tvals(ii,maxngs) = tvals(ii,maxngs)/icnt
c
         end do
c
c        Copy revised values into tvals(ii,1)
c
         do ii = 1,npts
c
            tvals(ii,1) = tvals(ii,maxngs)
c
         end do 
c
      endif
c
      NPLOTS = NPLOTS + 1
      WRITE (IPLOT,9012) NPLOTS,REF
c
c
      CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,npts,ANLY,NGS,
     >               99,Tstart,Tend,-HI,HI,IGNORS,0,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,BLABS,REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)
C
c
 9012 FORMAT(1X,'PLOT',I3,4X,A)

      return 
      end  
c
c
c
      subroutine load_2Dexpt_data(filename,datatitle,
     >                  maxix,maxiy,nix,niy,rzarrdata,raxis,zaxis,
     >                  rmin,rmax,zmin,zmax,dr,dz,ierr)
      use mod_params
      implicit none
      character*(*) filename,datatitle
      integer maxix,maxiy,nix,niy,ierr
c     include 'params'
      real rzarrdata(maxix,maxiy),raxis(maxix),zaxis(maxiy)
      real rmin,rmax,zmin,zmax,dr,dz
c
c     LOAD_2Dexpt_DATA:
c
c     This routine loads a set of R,Z,VALUE data into an array. It 
c     assumes that the R,Z are two sets of axes. The code finds the
c     minimum and maximum axis values and then checks to see if the 
c     spacing of the points is consistent with a regular grid. If
c     the grid is not regular dr and dz are set to zero. 
c 
      integer len,lenstr,ios,extstr
      external lenstr,extstr
      character*256 line 
      integer start 
      integer ix,iy,in,rnix,zniy,ir,iz
      real r,z,val,scalef,file_format
c
c     Cutoff options
c
      integer lcutopt,ucutopt
      real lcutval,ucutval
      real data_max,cutoff
c
c     Initialize
c 
      ierr =0  
      nix  =0
      niy  =0
      rnix = 0
      zniy = 0 
c
      rmin = HI
      rmax = -HI 
      zmin = HI
      zmax = -HI 
c
      scalef = 1.0 
      file_format=1
      datatitle = 'UNTITLED DATA'
c
c     Cutoff options
c     
      lcutopt = 0
      ucutopt = 0
      lcutval = 1.0
      ucutval = 1.0
c
      len = lenstr(filename) 
      open(tmpunit,FILE=filename(1:len),ERR=100,IOSTAT=ios)
      rewind(tmpunit)
c
      write(6,'(a,a)') 'OUTXSECTION:LOAD_2DEXPT_DATA:'//
     >                 'LOADING:',filename(1:len)
c
c     Read in header line and extract array size
c

      do while (nix.eq.0)
c
         read(tmpunit,'(a256)',END=100,ERR=200,IOSTAT=ios) line

c
c        Extract a Title if one is specified
c
         if (line(1:6).eq.'TITLE:') then
            len = extstr(line(7:),start)
            datatitle = line(6+start:len)
c            datatitle = line(8:LEN_TRIM(line))
            len = lenstr(datatitle)
            write(6,*) 'TITLE:',datatitle(1:len),':'
         endif
c          
c        Extract File Format if specified
c
         if (line(1:7).eq.'FORMAT:') then
            read(line(8:),*) file_format
         endif 
c
c        Extract Scaling Factor if specified
c
         if (line(1:7).eq.'SCALEF:') then 
           read(line(8:),*) scalef
         endif
c
c        Cutoffs for data manipulation
c
c        Extract Minimum Cutoff option
c
c        Options: 0 = off
c                 1 = percentage of maximum value
c                 2 = absolute cutoff
c
         if (line(1:8).eq.'LCUTOPT:') then 
           read(line(9:),*) lcutopt
         endif
c
c        Extract Minimum Cutoff value
c
         if (line(1:8).eq.'LCUTVAL:') then 
           read(line(9:),*) lcutval
         endif
c
c        Extract Maximum Cutoff option
c
c        Options: 0 = off
c                 1 = percentage of maximum value 
c                 2 = absolute cutoff
c
         if (line(1:8).eq.'UCUTOPT:') then 
           read(line(9:),*) ucutopt
         endif
c
c        Extract Maximum Cutoff value
c
         if (line(1:8).eq.'UCUTVAL:') then 
           read(line(9:),*) ucutval
         endif
c
C        Extract number of lines in the data for FORMAT 1
C
         if (file_format.eq.1.and.line(1:5).eq.'DATA:') then 

            read(line(6:),*) nix,niy
c
            if (nix.gt.maxix.or.niy.gt.maxiy) then 
c
               write(6,'(a)') 'LOAD_RZARRAY_DATA ERROR: '//
     >           'DATA STORAGE REQUIRED IS TOO LARGE'
               write(6,'(2(a,i5))') 'NIX = ',nix,' MAXIX = ',maxix
               write(6,'(2(a,i5))') 'NIY = ',niy,' MAXIY = ',maxiy
c
               ierr = -1
               return
            endif 
c
         endif 
c
c        Extract R-axis for format 2  
c
         if (file_format.eq.2.and.line(1:6).eq.'RAXIS:') then 

            read(line(7:),*) rnix

            if (rnix.gt.maxix) then 
c
               write(6,'(a)') 'LOAD_RZARRAY_DATA ERROR: '//
     >           'DATA STORAGE REQUIRED IS TOO LARGE'
               write(6,'(2(a,i5))') 'NIX = ',rnix,' MAXIX = ',maxix
c
               ierr = -1
               return
            endif 
c
            do in = 1,rnix                     
c               
               read(tmpunit,*,END=100,ERR=200,IOSTAT=ios) raxis(in) 
c
            end do                      
c
         endif   
c
c        Extract Z-axis for format 2  
c
         if (file_format.eq.2.and.line(1:6).eq.'ZAXIS:') then 

            read(line(7:),*) zniy

            if (zniy.gt.maxiy) then 
c
               write(6,'(a)') 'LOAD_RZARRAY_DATA ERROR: '//
     >           'DATA STORAGE REQUIRED IS TOO LARGE'
               write(6,'(2(a,i5))') 'NIY = ',zniy,' MAXIY = ',maxiy
c
               ierr = -1
               return
            endif 
c
            do in = 1,zniy                     
c               
               read(tmpunit,*,END=100,ERR=200,IOSTAT=ios) zaxis(in) 
c
            end do                      
c                                  
         endif   
c
         if (file_format.eq.2.and.line(1:5).eq.'DATA:') then 
            nix = rnix
            niy = zniy
         endif 
c  
      end do  
c
c     Read in array from datafile
c
      do ix = 1,nix
c
         do iy = 1,niy
c
            read(tmpunit,'(a256)',END=100,ERR=100,IOSTAT=ios) line
c
            if (lenstr(line).gt.1.and.line(1:1).ne.'$'.and.
     >          line(1:5).ne.'DATA:') then 
c
               if (file_format.eq.1) then 
                  read(line,*) r,z,val 
c
                  rmin = min(r,rmin)  
                  rmax = max(r,rmax)  
                  zmin = max(z,zmin)  
                  zmax = max(z,zmax)  
c
                  raxis(ix) = r 
                  zaxis(iy) = z         
                  rzarrdata(ix,iy) = val * scalef
c
               elseif (file_format.eq.2) then 

                  read(line,*) val,in,ir,iz
c
                  rzarrdata(ir+1,iz+1) = val * scalef
c                 
               endif
c
            endif 
c
         end do
 
      end do
c
c     Determine if the data is on a regularly spaced array  
c
      dr = (rmax-rmin)/nix
      dz = (zmax-zmin)/niy
c
c     Check for axis spaced evenly within 1% threshold
c       
c     Raxis 
c
      do ix = 2,nix
c
         if (abs((raxis(ix)-raxis(ix-1))-dr).gt.0.01*dr) then 
c
            dr = 0.0
            exit 
c
         endif  
c
      end do 
c
c     Zaxis
c
      do iy = 2,niy
c
         if (abs((zaxis(iy)-zaxis(iy-1))-dz).gt.0.01*dz) then 
c
            dz = 0.0
            exit 
c
         endif  
c
      end do 
c
      close(tmpunit)
c
c     Calculate and apply cutoffs if any
c
      if (lcutopt.ne.0.or.ucutopt.ne.0) then 

c
c        Process upper absolute cutoff first - affects maximum in data
c
         if (ucutopt.eq.2) then 
c
c           Determine cutoff value
c           
            cutoff = ucutval
c
c           Apply cutoff to data replacing values above the cutoff with 0.0
c
            do ix = 1,nix
               do iy = 1,niy
                  if (rzarrdata(ix,iy).gt.cutoff) then
                     rzarrdata(ix,iy) = 0.0
                  endif
               end do
            end do
c
         endif


c
c        If percentage cutoffs are specified then 
c        determine maximum value in the data
c         
         
         if (lcutopt.eq.1.or.ucutopt.eq.1) then 

            data_max = -HI
            
            do ix = 1,nix
               do iy = 1,niy
                  data_max = max(data_max,rzarrdata(ix,iy))
               end do
            end do

         endif

c
c        Process lower cutoff
c        
         if (lcutopt.ne.0) then 
c
c           Determine cutoff value
c           
            if (lcutopt.eq.1) then
               cutoff = lcutval * data_max
            elseif (lcutopt.eq.2) then 
               cutoff = lcutval
            endif
c
c           Apply cutoff to data replacing values below the cutoff with 0.0
c
            do ix = 1,nix
               do iy = 1,niy
                  if (rzarrdata(ix,iy).lt.cutoff) then
                     rzarrdata(ix,iy) = 0.0
                  endif
               end do
            end do
c
         endif
c
c        Process upper relative cutoff
c
         if (ucutopt.eq.1) then 
c
c           Determine cutoff value
c           
            cutoff = ucutval * data_max
c
c           Apply cutoff to data replacing values above the cutoff with 0.0
c
            do ix = 1,nix
               do iy = 1,niy
                  if (rzarrdata(ix,iy).gt.cutoff) then
                     rzarrdata(ix,iy) = 0.0
                  endif
               end do
            end do
c
         endif
c
      endif

c
c      write(6,*) 'DEBUG: Nix,niy:',nix,niy,rnix,zniy
c
c      do ix = 1,nix
c         do iy = 1,niy
c            write(6,'(a,2i4,3(1x,g12.5))') 'DATA:',ix,iy,raxis(ix),
c     >            zaxis(iy),rzarrdata(ix,iy)
c         end do
c      end do 
c
      return 
c
c     Deal with error opening input file
c
 100  write(6,'(a,i4)') 'LOAD_2DEXPT_DATA ERROR: EOF ON '
     >              //filename(1:len)
 200  write(6,'(a,i4)') 'LOAD_2DEXPT_DATA ERROR: FILE '
     >              //filename(1:len)//' HAS AN ERROR = ',ios
      if (nix.eq.0.or.niy.eq.0) then  
         write(6,'(a,i4)') 'LOAD_2DEXPT_DATA ERROR: DATA '
     >              //'INPUT LINE NOT FOUND'
      endif 

      ierr = ios 
      close(tmpunit)
      return 
      end  
c
c
c
      real function getvalexpt2D(r,z,
     >                  maxix,maxiy,nix,niy,expt_array,raxis,zaxis,
     >                  expt_rmin,expt_rmax,expt_zmin,expt_zmax,
     >                  expt_dr,expt_dz,zero_opt)
      implicit none
c
      real r,z 
      integer  maxix,maxiy,nix,niy,zero_opt
      real expt_array(maxix,maxiy),raxis(maxix),zaxis(maxiy)
      real expt_rmin,expt_rmax,expt_zmin,expt_zmax,expt_dr,expt_dz
c
c     GETVALEXPT2D:
C
C     This routine treats the 2D experimental data as being on a 
c     grid with axes defined by the raxis and zaxis arrays - which 
c     define the centre points of each bin. The bins are assumed to 
c     extend from the half-way points of the nearest neighbours on each
c     side. 
c
c     At present the routine does not perform any interpolation - it 
c     simply returns the value in the appropriate bin. 
c        
c     If expt_dr and expt_dz are both non-zero - the code assumes that
c     these values contain the spacings of a regular grid for the data
c     which starts at the min values and has the last data point at the 
c     max values - this significantly simplifies the calculation of 
c     which bin the point is in. 
c
      integer ix,iy,ipos,jpos
      external ipos,jpos
c
c     Initialize as zero return value
c
      getvalexpt2D = 0.0            
c
      if (expt_dr.gt.0.0.and.expt_dz.gt.0.0) then
c
c        Calculate R bin
c           
         ix = INT((r-(expt_rmin-expt_dr/2.0))/expt_dr) +1
c
c        Calculate Z bin 
c
         iy = INT((z-(expt_zmin-expt_dz/2.0))/expt_dz) +1
c
c
c     Non-uniform grid
c
      else
c
c        Determine R-index
c                   
         if (raxis(1).lt.raxis(nix)) then

            ix = ipos(r,raxis,nix)
c
            if (ix.gt.1) then 
c
               if ((r-raxis(ix-1)).lt.(raxis(ix)-r)) then 
c
                  ix = ix -1 
c
               endif
c
           endif                              
c
         else
c
            ix = jpos(r,raxis,nix)
c
            if (ix.lt.nix) then 
c
               if ((r-raxis(ix)).lt.(raxis(ix+1)-r)) then 
c
                  ix = ix +1
c
               endif
c
            endif                              
c
         endif 
c           
c        Determine Z-index           
c
         if (zaxis(1).lt.zaxis(niy)) then

            iy = ipos(z,zaxis,niy)
c
            if (iy.gt.1) then 
c
               if ((z-zaxis(iy-1)).lt.(zaxis(iy)-z)) then 
c
                  iy = iy -1 
c
               endif
c
           endif                              
c
         else
c
            iy = jpos(z,zaxis,niy)
c
            if (iy.lt.niy) then 
c
               if ((z-zaxis(iy)).lt.(zaxis(iy+1)-z)) then 
c
                  iy = iy +1
c
               endif
c
            endif                              
c
         endif 
c
      endif  
c
c     Extract value 
c
c     Out of range error first
c
      if (ix.lt.1.or.ix.gt.nix.or.iy.lt.1.or.iy.gt.niy) then
c
         getvalexpt2D = 0.0            
c
      else
c
         getvalexpt2D = expt_array(ix,iy)
c
      endif
c
      if (zero_opt.eq.1.and.getvalexpt2D.lt.0.0) getvalexpt2D=0.0 
c
      return
      end
c
c
c
      subroutine create_2Ddata(divdata,maxix,maxiy,nix,niy,
     >                  div2D_array,raxis,zaxis,flag)
      use mod_params
      implicit none
c     include 'params'
c
      integer maxix,maxiy,nix,niy,flag
      real raxis(nix),zaxis(niy)
      real divdata(maxnks,maxnrs)
      real div2D_array(maxix,maxiy) 
c
c     CREATE_2DDATA: This routine performs a Monte Carlo integration
c                    on a specified grid over an input DIVIMP density 
c                    array. The accuracy is limited by the number
c                    of Monte Carlo samples taken in each cell. The idea
c                    is to convert a DIVIMP array into one that 
c                    exactly matches the experimentally reported 
c                    results. 
c
c     FLAG=0 - perform Monte Carlo integration  
c     FLAG=1 - map cell centre values
c
c 
      integer maxsamples
      parameter(maxsamples=100) 
c
c     Local Variables
c
      integer in,ix,iy 
      double precision tmpsum 
      real r,z  
c
      real xrange_l,xrange_g,yrange_l,yrange_g
c slmod begin - new
      real getvalrz,get_randc
      external getvalrz,get_randc
c
c      real getvalrz
c      external getvalrz
c slmod end
c
      do ix = 1,nix
c
c        Define X/R range 
c
         if (flag.eq.0) then         
            if (ix.eq.1) then 
c
               xrange_g = (raxis(ix+1) - raxis(ix))/2.0
               xrange_l = xrange_g
c
            elseif (ix.eq.nix) then 
c
               xrange_l = (raxis(ix) - raxis(ix-1)) / 2.0
               xrange_g = xrange_l
c  
            else
c
               xrange_l = (raxis(ix) - raxis(ix-1))/2.0
               xrange_g = (raxis(ix+1) - raxis(ix))/2.0
c
            endif
         endif
c
         do iy = 1,niy
c
            if (flag.eq.0) then
c
c              Define Y/Z range  
c
               if (iy.eq.1) then 
c
                  yrange_g = (zaxis(iy+1) - zaxis(iy))/2.0
                  yrange_l = yrange_g
c
               elseif (iy.eq.niy) then 
c
                  yrange_l = (zaxis(iy) - zaxis(iy-1)) / 2.0
                  yrange_g = yrange_l
c
               else
c
                  yrange_l = (zaxis(iy) - zaxis(iy-1)) / 2.0
                  yrange_g = (zaxis(iy+1) - zaxis(iy))/2.0
c 
               endif
c  
            endif
c
c           Loop for cell to perform integral
c
            if (flag.eq.0) then 

               tmpsum = 0.0 
c
               do in = 1,maxsamples
c
                  r = get_randc(raxis(ix),xrange_l,xrange_g) 
                  z = get_randc(zaxis(iy),yrange_l,yrange_g) 
c
                  tmpsum = tmpsum + getvalrz(divdata,r,z)
c
c                  write(6,'(a,3i4,4(1x,g12.5))') 'TMPSUM:',
c     >                   ix,iy,in,r,z,tmpsum
c
               end do
c
               div2D_array(ix,iy) = tmpsum / maxsamples
c
            elseif (flag.eq.1) then 
c
               div2D_array(ix,iy) = 
     >                 getvalrz(divdata,raxis(ix),zaxis(iy))
            endif 

c
c            write(6,'(a,2i4,4(1x,g12.5))') 'Create2D:',ix,iy,
c     >       div2D_array(ix,iy),getvalrz(divdata,raxis(ix),zaxis(iy))
c
         end do
c
      end do
c
      return
      end
c
c
c
      real function get_randc(startval,extent1,extent2)
      implicit none  
      real startval,extent1,extent2
c
c     GET_RANDC: This routine returns a random coordinate 
c                in the range 
c                [startval-extent1,startval+extent2]
c
      double precision seed
      real   randval
c
      call surand2(seed,1,randval)
c
      get_randc = (startval-extent1) + randval * (extent2+extent1)
c
      return 
      end




