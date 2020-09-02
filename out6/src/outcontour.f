c     -*-Fortran-*-
c
c
      subroutine plot_contour(iselect,istate,
     >                  iexpt,optval,
     >                  iopt,job,title,table,nplots,
     >                  iplot,nizs,
     >                  minscale,maxscale,
     >                  ierr) 
      use error_handling
      use subgrid_plots
      use hc_get
      use mod_params
      implicit none
c
c     jdemod - use optval as a scaling factor if it is greater than zero
c
c
c     include 'params'
c
      integer iselect,istate,iexpt,ierr 
      real optval,minscale,maxscale
      character*(*) job,title,table
      integer iopt,nplots,iplot,nizs
c
c     PLOT_CONTOUR: This routine plots a generalized CONTOUR plot.
c
c
c     Local variables
c
      real tmpplot(maxnks,maxnrs)
      real mfact
c
c     Plot labels   
c
      character*36 blabs(2),xlab,ylab
      character*44 refd,refe,plane,anly,nview
      CHARACTER*72 SMOOTHd,smoothe
      character*72 datatitle
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
c     Contour plot options 
c
      real xcen,ycen,xnear,ynear
      real xmin,xmax,ymin,ymax 
      real uconts(maxpts)
      integer icntr,ncntr
c
c     Variables for subgrids
c     
      integer :: iflag
      real,allocatable :: subgrid_data(:,:)
      real,allocatable :: subgrid_raxis(:),subgrid_zaxis(:)
c
      real dr,dz,cellvol,tmp_abs,photons_per_part
c
c      real,external :: gabsfac
c
c     Dummy arguments
c
      real xouts(1),youts(1) 
c
c
c     
      integer itype
      real scalef
c     
c
c     Initialization
c
      ierr = 0
      icntr = 0
      ncntr = 0
c     
c     Default type is contour plot (option 0) with standard units
c
      itype = 0
      scalef = 1.0

      if (optval.ne.0.0) then 
         scalef = optval
      endif
c
c     If IOPT is 2 - then set up for giving units in /m3/s/sr
c     NOTE: should only use IOPT=2 with appropriate ISELECT values
c
      if (iopt.eq.2) then 
         itype=3
         scalef = scalef * 1.0 / (4.0 * PI) 
      elseif (iopt.eq.3) then 
         itype=4
         scalef = scalef * 1.0 / (4.0 * PI * 1.0e6) 
      endif


c
c     Load additional contour plot parameters 
c
c     write(0,*) 'Reading input:'
c

      call rdg_contopts(graph,icntr,ncntr,uconts,maxpts,
     >                  xcen,ycen,xnear,ynear,ierr)
c
      if (ierr.ne.0) return
c
c     Set up plotting ranges
c
      xmin = xcen - xnear
      xmax = xcen + xnear
      ymin = ycen - ynear
      ymax = ycen + ynear
c
c     Set up colours based on contour options  
c
c     Contour option 6 is used to select a colour scheme but the contouring 
c     option is the same as 3 so icntr is reset to 3.
c
      write(0,*) 'cntr:',icntr,ncntr
c
      if (icntr.eq.0.or.icntr.eq.2) then
          call setup_col(ncntr+1,2)
      elseif (icntr.eq.1) then 
          ncntr = 10
          call setup_col(ncntr+1,2)
      elseif (icntr.eq.3.or.icntr.eq.4) then
          !ncntr = ncntr
c         IPP/09 - Krieger - this should be 3 instead of 2?
          call setup_col(ncntr+1,3)
      elseif (icntr.eq.5) then
          ncntr = 41
          call setup_col(ncntr+1,4)
      elseif (icntr.eq.6) then 
          !ncntr = ncntr
          call setup_col(ncntr+1,6)
          icntr = 3
      elseif (icntr.eq.7) then 
          !ncntr = ncntr
          call setup_col(ncntr+1,7)
          icntr = -3
      endif
c
c      write(0,'(a,2i4,6(1x,g12.5))') 'PLOT_CONT:',icntr,ncntr,
c     >      xcen,ycen,xnear,ynear
c

c
c
c     Set NVIEW, PLANE ...
c
      NVIEW  = 'GENERALIZED CONTOUR PLOT'
      plane  = ' '
      REFD    = ' ' 
      REFE    = ' ' 
      XLAB   = '   R  (M)'
      YLAB   = '   Z  (M)'
      anly   = ' '
      smoothd = ' '
      smoothe = ' '
c
c     Based on the values of iselect and istate - load up the required 
c     2D array for the CONTOUR plot 
c
c     Scaling by absfac (if required) is performed in the 
c     load_divdata_array routine. 
c
      call rzero(tmpplot,maxnks*maxnrs)
c
c     Do not load DIVIMP data if ISELECT is zero - only plotting experimental 
c     data. 
c
c
c     Load subgrid based data
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
     >                         smoothd,refd,plane,nizs,ierr)

         call calc_subgrid_axis(subgrid_raxis,sg_rmin,sg_rmax,sg_rdim)
         call calc_subgrid_axis(subgrid_zaxis,sg_zmin,sg_zmax,sg_zdim)


         if (iselect.eq.34.or.iselect.eq.35) then 
            dr = (sg_rmax-sg_rmin)/real(sg_rdim)
            dz = (sg_zmax-sg_zmin)/real(sg_zdim)
            cellvol = dr*dz
            tmp_abs = gabsfac()


            if (tmp_abs.ne.0.0) then 
               photons_per_part = sum(subgrid_data(:,:))*cellvol/tmp_abs
            else
               photons_per_part = sum(subgrid_data(:,:))*cellvol 
            endif

            write(0,'(a,a,6(1x,g12.5))') 'Photons/part:', trim(plane),
     >                dr,dz,cellvol,tmp_abs,
     >                sum(subgrid_data(:,:)),photons_per_part
            write(6,'(a,a,6(1x,g12.5))') 'Photons/part:', trim(plane),
     >                dr,dz,cellvol,tmp_abs,
     >                sum(subgrid_data(:,:)),photons_per_part


            write(anly,'(a,1x,g12.5)') 'Photons/part:',photons_per_part
c            

         endif
c
c        Apply scaling factor if needed
c
         if (scalef.ne.1.0) then 
            subgrid_data = subgrid_data * scalef
         endif

      elseif (iselect.ne.0) then 
c
         call load_divdata_array(tmpplot,iselect,istate,itype,
     >                         smoothd,refd,plane,nizs,ierr)
c
c        Apply scaling factor if needed
c
         if (scalef.ne.1.0) then 
            tmpplot = tmpplot * scalef

         endif
c
      endif 
c
c     Check return code from data load and exit if non-zero   
c
      if (ierr.ne.0) then

c
c         Deallocate storage assigned for subgrid plots if it has been allocated
c      
         if (allocated(subgrid_data)) deallocate(subgrid_data)
         if (allocated(subgrid_raxis)) deallocate(subgrid_raxis)
         if (allocated(subgrid_zaxis)) deallocate(subgrid_zaxis)

c        IPP/09 - Krieger - should reset colors here
         call setup_col(ncntr,abs(icntr))
         return 
      endif


c
c     Check to see if a CONTOUR of a 2D experimental data array is 
c     also requested.
c
c     A value of IEXPT of -2 means to also convert the DIVIMP data 
c     to the same format as the experimental data prior to the 
c     calculation of the values to be plotted.
c
c
c     A value of IEXPT of -3 means to also convert the DIVIMP data 
c     to the same format as the experimental data prior to the 
c     calculation of the values to be plotted. However, the experimental
c     data itself will NOT be plotted.
c
c     A value of IEXPT = -4 means to load and plot the experimental data ONLY
c
c
      if (iexpt.eq.-1.or.iexpt.eq.-2.or.iexpt.eq.-3.or.iexpt.eq.-4) then  
c
         call rdfn(graph,filename,ifnopt,ierr)
c
         if (ierr.ne.0) then  
c
            len = lenstr(filename)
            call errmsg('OUTCONTOUR:PLOT_CONTOUR',
     >                   'ERROR LOADING 2D EXPT'
     >                   //' DATAFILE NAME: DATA='//
     >                  filename(1:len_trim(filename)))
            iexpt = 0
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
               call errmsg('OUTCONTOUR:PLOT_CONTOUR',
     >                   'ERROR READING 2D EXPT'
     >                   //' DATAFILE: INPUT FILE NAME:'
     >                   //filename(1:len))
               iexpt = 0
c
            else

               len=lenstr(datatitle)
               REFE = 'EXPT '//datatitle(1:len)
               SMOOTHE= 'LINE EMISSION (PHOTONS/M^3/S)'  
c
            endif 
c
         endif
c
c        Convert the DIVIMP array to an experimental data format
c 
         if (iexpt.eq.-2.or.iexpt.eq.-3) then 
c
            call create_2Ddata(tmpplot,maxix,maxiy,nix,niy,
     >                  div2D_array,raxis,zaxis,0)
c
         endif
c
      endif
c
c     write(0,*) 'CONTOUR DATA LOADED:'
c

c
c     Now that the data has been successfully loaded - plot the 
c     respective contour plots.
c
c     Plot the DIVIMP data contour 
c
      if (iselect.ne.0.and.iexpt.ne.-4) then 


c
c        Plot the subgrid contours
c
         if (iselect.eq.32.or.iselect.eq.33.or.iselect.eq.34.or.
     >       iselect.eq.35) then 
c
            WRITE (IPLOT,9012) NPLOTS,REFD
            CALL GRTSET (TITLE,REFD,NVIEW,PLANE,JOB,XMIN,XMAX,
     >         YMIN,YMAX,TABLE,XLAB,YLAB,2,SMOOTHD,1,ANLY,Ncntr)

            call CONTOURXY(1,Ncntr,subgrid_data,
     >                    xmin,xmax,ymin,ymax,ncntr,uconts,
     >                    icntr,minscale,maxscale,
     >                    sg_rdim,sg_zdim,sg_rdim,sg_zdim,
     >                    subgrid_raxis,subgrid_zaxis,1)
c
         elseif (iexpt.eq.-2.or.iexpt.eq.-3) then 
c
            WRITE (IPLOT,9012) NPLOTS,REFD
            CALL GRTSET (TITLE,REFD,NVIEW,PLANE,JOB,XMIN,XMAX,
     >         YMIN,YMAX,TABLE,XLAB,YLAB,2,SMOOTHD,1,ANLY,Ncntr)

            call CONTOURXY(1,Ncntr,div2D_array,
     >                    XMIN,XMAX,YMIN,YMAX,ncntr,uconts,
     >                    icntr,minscale,maxscale,
     >                    maxix,maxiy,nix,niy,raxis,zaxis,1)
c
         else
c
            WRITE (IPLOT,9012) NPLOTS,REFD
            CALL GRTSET (TITLE,REFD,NVIEW,PLANE,JOB,XMIN,XMAX,
     >         YMIN,YMAX,TABLE,XLAB,YLAB,2,SMOOTHD,1,ANLY,Ncntr)
            CALL CONTOUR (1,Ncntr,tmpplot,1,1,1,1.0,1.0,1.0,
     >                XOUTS,1,1,YOUTS,1,1,
     >                XMIN,XMAX,YMIN,YMAX,
     >                ncntr,uconts,icntr,minscale,maxscale)
c
         endif   



      endif
c
c     Plot the EXPERIMENTAL data contour 
c
      if (iexpt.eq.-1.or.iexpt.eq.-2.or.iexpt.eq.-4) then 
c
         plane = ' '
c
         WRITE (IPLOT,9012) NPLOTS,REFE
         CALL GRTSET (TITLE,REFE,NVIEW,PLANE,JOB,XMIN,XMAX,
     >      YMIN,YMAX,TABLE,XLAB,YLAB,2,SMOOTHE,1,ANLY,Ncntr)

         call CONTOURXY(1,Ncntr,expt_array,
     >                 XMIN,XMAX,YMIN,YMAX,ncntr,uconts,
     >                 icntr,minscale,maxscale,
     >                 maxix,maxiy,nix,niy,raxis,zaxis,1)
c

      endif 



c
c     Some custom print out for Adam - there must be a better place
c     for it - but include and comment it out for now. 
c
! ammod begin.
c      open (unit=95,file="temp_out.out")
c      ! Write data to integrate emission contour.
c      R_zero = 1.72
c      Z_zero = 0.0
c      Do Print_Cell=35,60,1 ! Going CCW.
c	 Int_Sum = 0.0
c	 Area_Sum = 0.0
c         write (95,*) "By Cell Column:",Print_Cell
c         write (95,9304) "Cell","Ring","Rcentre","Zcentre","Angle",
c     >     "Area","Intensity"
c     
c         Do Print_Ring = irsep,irwall-1 ! Going from out to in.
c
!           write (95,9305) Print_Cell,Print_Ring,
!     >       rs(print_cell,print_ring),zs(print_cell,print_ring),
!     >       ATAN((ABS(R_Zero-rs(print_cell,print_ring)))/
!     >       (ABS(Z_Zero-zs(print_cell,print_ring))))*180.0/3.14159,
!     >       kareas(print_cell,print_ring),tmpplot
!     >       (print_cell,print_ring)
c        Int_Sum = Int_Sum + tmpplot(print_cell,print_ring)*
c     >     kareas(print_cell,print_ring)
c        Area_Sum = Area_Sum + kareas(print_cell,print_ring)
c     
c 9304 Format (10(A10))
c 9305 Format (I10,I10,3F10.3,10E10.3)
!	   ,"cell",M,"counts",DDLIMS(M,N,1),
!     >     "length",ksb(M,N)-ksb(M-1,N),"Center S",kss(M,N),
!     >     "r centre",rs(M,N),"z centre",zs(M,N),
c         End Do
c	 write (95,*) "Column sum:   ",Int_Sum," Area Sum:  ",Area_Sum
c      End Do
c      close (95)
! ammod end.


c
c      Deallocate storage assigned for subgrid plots if it has been allocated
c      
      if (allocated(subgrid_data)) deallocate(subgrid_data)
      if (allocated(subgrid_raxis)) deallocate(subgrid_raxis)
      if (allocated(subgrid_zaxis)) deallocate(subgrid_zaxis)

c
c     Return from generalized CONTOUR routine
c
c     IPP/09 - Krieger - should reset colors here
      call setup_col(ncntr,abs(icntr))
      return
c
c     Format statements 
c

 9012 FORMAT(1X,'PLOT',I3,4X,A)
c
c
      end






