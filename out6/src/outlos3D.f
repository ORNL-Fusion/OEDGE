      subroutine plot_3Dimage(iselect,istate,iaxis,minsteps,
     >                  stepsize,optval,
     >                  iopt,job,title,table,nplots,
     >                  iplot,nizs,ierr) 
      implicit none
c
      include 'params'
      include 'slout'
c
      integer iselect,istate,iaxis,minsteps,ierr 
      real optval,stepsize 
      character*(*) job,title,table
      integer iopt,nplots,iplot,nizs
c
c     PLOT_3Dimage: This routine plots a generalized 3D LOS plot. 
c                 Given a camera observation position this routine
c                 will simulate a 1D/2D camera view of the selected 
c                 simulated signal from the specified position. 
c
c                 1D views will generate regular line plots - 2D views
c                 will generate pixellated contour plots with the 
c                 possibility of producing an actual camera image in 
c                 a standard graphics format. (TIFF likely).
c
c                 Selection between 1D and 2D views is made via the 
c                 iopt input value.
c                  
c
c     Local variables
c
c
c     Plot labels   
c
      character*256 graph,camera
      CHARACTER*72 SMOOTH
      character*36 blabs(2),xlab,ylab
      character*44 ref,plane,anly,nview
      character*32 datatitle
c 
c     Other
c
      integer ir,ic,slopt_tmp,reflect_opt 
      real*8 maxval,minval,view_height,view_width
      real*8 dr,dz
      integer nix,niy,maxix,maxiy
c
c     Contour plot options 
c
      real xcen,ycen,xnear,ynear
      real xmin,xmax,ymin,ymax 
      real uconts(maxpts)
      integer icntr,ncntr
      real minscale,maxscale
c
c      INTEGER IGNORS(MAXNGS),ngs
c
c      REAL TOUTS(MAXTHE),TWIDS(MAXTHE),TVALS(MAXTHE,MAXNGS)
c      real tmptvals(0:maxthe-1)
c      REAL THEMIN,THEMAX
c
      integer len,lenstr
c
c     Variables for 2D experimental data - iexpt = -1
c
c      character*256 graph
c      character*256 filename
c      integer  maxix,maxiy,nix,niy,ifnopt
c      parameter(maxix=100,maxiy=100)
c      real expt_array(maxix,maxiy),raxis(maxix),zaxis(maxiy)
c      real div2D_array(maxix,maxiy)
c      real expt_rmin,expt_rmax,expt_zmin,expt_zmax,expt_dr,expt_dz
c
c
c     Set up Allocatable array to hold the image
c
      real image(:,:),raxis(:),zaxis(:)  
      allocatable image,raxis,zaxis
      integer xres,yres
      real*8 position(3),direction(3),upvec(3),rightvec(3),lookat(3)
      logical mirror
c
      real tmpplot(maxnks,maxnrs)
c
c
      real*8 vect_mag
      external vect_mag
c
c---------------------------------------------------------------------- 
c
c     Initialization
c
c     Use IAXIS to carry reflection option information for this plot
c
      reflect_opt = iaxis
c
c---------------------------------------------------------------------- 
c
c     Load DIVIMP data to be imaged.
c
c     Based on the values of iselect and istate - load up the required 
c     2D array for the LOS plot 
c
c     Scaling by absfac (if required) is performed in the 
c     load_divdata_array routine. 
c
      call rzero(tmpplot,maxnks*maxnrs)
c
      call load_divdata_array(tmpplot,iselect,istate,1,
     >                         ylab,blabs(1),ref,nizs,ierr)
c
      write(6,'(a,3i4)') 'DIVIMP data loaded: ',iselect,istate,ierr
c
c     Check return code from data load and exit if non-zero   
c
      if (ierr.ne.0) then
         write(6,'(a,i5)') 'ERROR Loading Integration DATA for Image',
     >                  ierr
         return 
      endif
c
c
c---------------------------------------------------------------------- 
c
c     Load camera view data - 
c
c     - observation position 
c     - definition of camera viewing area 
c     - direction/location of center of image
c
c     - resolution X x Y  - if X or Y is 1 then 1D
c     - bit depth of image output
c 
c     Alternative method -
c    
c     - load individual LOS including plotting axis coordinate. 
c
      call rdg_camera(graph,xres,yres,position,direction,
     >                upvec,rightvec,
     >                lookat,camera,ierr)
c
c      write(0,'(a,2i4)') 'CAMERA data loaded: ',ierr
c
c NOTE: In the POVRAY data format used for the cameras -
c       X is North, Z is East and Y is Up. i.e. Right-Handed with
c       the vertical axis as Y not Z. 
c
c       The code in DIVIMP has been written in such a way that the 
c       Z axis is up corresponding to the typical tokamak 
c       nomenclature for 2D (R for horizontal and Z for Up). 
c
c       This requires that all of the input vectors be transformed 
c       from the right hand POVRAY system to the DIVIMP system.
c
c       This is accomplished by the following:
c
c       Xpov -> Ydivimp
c       Ypov -> Zdivimp
c       Zpov -> Xdivimp
c
c       All that is different is the labels for the different axes
c       - their locations and orientations are the same. By transforming
c       the vectors in this way the code below can remain the same and 
c       confusion about coordinate systems inside the code can be 
c       avoided. 
c
c       The alternative is to re-write all of the LOS code to use the
c       POVRAY coordinate system instead.
c
c
c 210 parallel:
c     location  <-203.21, 73.665, -113.035>
c     direction <168.306,31.659,-9.376>
c     up        <0,  72.611,   0>
c     right   <99.014, 0,  0>
c     look_at <-34.904,105.324,-122.411>
c
c
c     Check return code from camera load and exit if non-zero   
c
      if (ierr.ne.0) then
         write(6,'(a,i4)') 'ERROR Loading CAMERA Specifications: ',
     >                     ierr 
         return 
      endif
c
c     NOTE: If the component of the right vector is less than zero -
c           this implies that the orientation of the image should
c           by flipped right to left (mirror image). 
c
c
      if (rightvec(1).LT.0.0) then
         mirror = .false.
      else
         mirror = .true.
      endif
c
c     Transform camera coordinates: 
c
c       Xpov -> Ydivimp
c       Ypov -> Zdivimp
c       Zpov -> Xdivimp
c
      call shiftcoord(position) 
      call shiftcoord(direction) 
      call shiftcoord(upvec) 
      call shiftcoord(rightvec) 
      call shiftcoord(lookat) 
c
c
c      write(0,'(a,3f10.6)') 'Camera Data loaded: pos= ',
c     >                      position(1),position(2),position(3)
c      write(0,*) xres,yres 
c
c---------------------------------------------------------------------- 
c
c     Allocate space for image and check to see if space available 
c
c      write(0,'(a)') 'Beginning Image Allocation'
c
      allocate (image(xres,yres),STAT=ierr)
c
      if (ierr.ne.0) then 
         write(0,*) 'IMAGE array could not be allocated: ',xres,yres
         write(6,*) 'IMAGE array could not be allocated: ',xres,yres
         return
      endif        
c
c     Zero image 
c
      call rzero(image,xres*yres)
c
c      write(0,'(a)') 'Image Allocated'
c
c---------------------------------------------------------------------- 
c
c     Generate image
c
c
      call calc_image(image,xres,yres,position,direction,
     >                upvec,rightvec,lookat,tmpplot,maxval,minval,
     >                view_height,view_width,minsteps,stepsize,
     >                mirror,reflect_opt)
c     
c      write(0,'(a)') 'Image calculated'
c
c---------------------------------------------------------------------- 
c
c     Contour plot output 
c
      if (iopt.eq.1.or.iopt.eq.3) then   

c
c        Plot the image as a contour plot
c
         maxix = xres
         maxiy = yres
         nix   = xres
         niy   = yres       
c
c         write(0,'(a)') 'Plotting contours'
c
c        Read in contour options if contour plot was specified  
c
         call rdg_contopts(graph,icntr,ncntr,uconts,maxpts,
     >                  xcen,ycen,xnear,ynear,ierr)

         if (ierr.ne.0) goto 100

c
c        Set up colours based on contour options  
c
         if (icntr.eq.0.or.icntr.eq.2) then
             call setup_col(ncntr+1,2)
         elseif (icntr.eq.1) then 
             ncntr = 10
             call setup_col(ncntr+1,2)
         elseif (icntr.eq.3.or.icntr.eq.4) then
             ncntr = ncntr
             call setup_col(ncntr+1,2)
         endif
c
c         write(0,'(a,2i4,8(1x,g12.5))') 'PLOT_CONT:',icntr,ncntr,
c     >         xcen,ycen,xnear,ynear,maxval,minval
c
c        Set NVIEW, PLANE ...
c
         NVIEW  = 'IMAGE CONTOUR PLOT'
         write(plane,'(a,f5.2,a)') 'IMAGE PLANE AT ',
     >                  vect_mag(direction),' (M)'
         XLAB   = '   X  (M)'
         YLAB   = '   Y  (M)'

         len = lenstr(camera)         
         smooth = camera(1:len)
         anly   = ' '
c
c        Set min and max levels for contour plot
c
         maxscale = maxval
         minscale = 0.0
c
c
c        Set up axes for the contour plot call - use the 
c        base camera definition space at the given distance.
c
c
c        Allocate space for axes and check to see if space available
c
         allocate (raxis(xres),STAT=ierr)
c
         if (ierr.ne.0) then 
            write(0,*) 'RAXIS array could not be allocated: ',xres
            write(6,*) 'RAXIS array could not be allocated: ',xres
            goto 100
         endif        
c
         allocate (zaxis(yres),STAT=ierr)
c
         if (ierr.ne.0) then 
            write(0,*) 'ZAXIS array could not be allocated: ',yres
            write(6,*) 'ZAXIS array could not be allocated: ',yres
            goto 100
         endif        
c
c        Assign axis values
c
         dr = view_width/xres
c
c        Column center coordinates
c
         do ic = 1,xres
            raxis(ic) = -(view_width/2.0) + dr/2.0 + (ic-1)*dr 
c            write(6,'(a,i5,2(1x,g12.5))') 'RAXIS:',ic,raxis(ic),dr 
         end do 
c
c        Row center coordinates
c
         dz = view_height/yres
         do ir = 1,yres  
            zaxis(ir) = (view_height/2.0) 
     >                         - dz/2.0 - (ir-1)*dz 
c            write(6,'(a,i5,2(1x,g12.5))') 'ZAXIS:',ir,zaxis(ir),dz 
         end do 
c
c        Set up plotting ranges
c
         if (xnear.ne.0.0.and.ynear.ne.0.0) then 
c
            xmin = xcen - xnear
            xmax = xcen + xnear
            ymin = ycen - ynear
            ymax = ycen + ynear
c     
         else
c
            xmin = raxis(1) - dr/2.0
            xmax = raxis(xres) + dr/2.0
            ymin = zaxis(yres) - dz/2.0
            ymax = zaxis(1) + dz/2.0
c
         endif 
c
c        Draw contour plot 
c
         WRITE (IPLOT,9012) NPLOTS,REF
         CALL GRTSET (TITLE,REF,NVIEW,PLANE,blabs(1),XMIN,XMAX,
     >         YMIN,YMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,Ncntr)
 
c
c        Avoid drawing cell boundaries 
c
         slopt_tmp = slopt
         slopt = 42
c
         call CONTOURXY(1,Ncntr,image,
     >                    XMIN,XMAX,YMIN,YMAX,ncntr,uconts,
     >                    icntr,minscale,maxscale,
     >                    maxix,maxiy,nix,niy,raxis,zaxis,0)
         
c
         slopt = slopt_tmp
c
         deallocate(raxis)
         deallocate(zaxis)

      endif
 100  continue   

c
c---------------------------------------------------------------------- 
c
c
c     Image type output
c
      if (iopt.eq.2.or.iopt.eq.3) then 
c
         write(6,'(a)') 'Saving Image as JPEG'
c
c
c        Save the image to a file
c
         call save_image(image,xres,yres,optval,maxval,minval)
c 

      endif
c
c---------------------------------------------------------------------- 
c
c     Free space assigned to image array
c
      deallocate(image)
c
c---------------------------------------------------------------------- 
c
c     Return from 3D image routine
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
      subroutine shiftcoord(vec)
      implicit none
      real*8 vec(3)
      
c
c     Shift the coordinates of the vector so that: 
c
c       Xpov -> Ydivimp
c       Ypov -> Zdivimp
c       Zpov -> Xdivimp
c
      real*8 temp  
      integer in
c
      temp = vec(3) 
      do in = 2,1,-1
         vec(in+1) = vec(in)            
      end do
      vec(1) = temp
      return
      end
c
c
c
      real*8 function vect_mag(a)
      implicit none
      real*8 a(3)
c
c     Returns the magnitude of a 3 vector
c 
      vect_mag = 
     >   (a(1)**2.0d0+a(2)**2.0d0+a(3)**2.0d0)**0.5d0
c
c      write(6,'(a,4(1x,g12.5))') 'VECT_MAG:',
c     >                      vect_mag,a(1),a(2),a(3) 
c
      return
      end  

c
c     
c
      subroutine save_image(image,xres,yres,calibration,maxval,minval)
      implicit none
      integer xres,yres
      real image(xres,yres),calibration
      real*8 maxval,minval
c
c     Set up and write out a JPEG formatted image
c
c
      integer temp_val 
      integer,parameter::SHORT=1
c
      integer*1 image_buffer(:,:)  
      allocatable image_buffer 
      integer ierr,len,lenstr,ic,ir
      external lenstr    
      integer image_count 
      data image_count /0/
      character*(100) fname,filename
c
c---------------------------------------------------------------------- 
c
c     Allocate space for image and check to see if space available 
c
c      write(0,*) 'Allocating IMAGE_BUFFER:' 
c
      allocate (image_buffer(xres,yres),STAT=ierr)
c
      if (ierr.ne.0) then 
         write(0,*) 'IMAGE_BUFFER array could not be allocated: ',
     >               xres,yres
         write(6,*) 'IMAGE_BUFFER array could not be allocated: ',
     >               xres,yres
         return
      endif        
c
c      write(0,*) 'IMAGE_BUFFER Allocated' 
c
c     If a calibration is not given for the camera - scale the 
c     maximum signal to maximum intensity
c
c      calibration = 3.44E+22 / 255.0 / 4.0 / 3.1415
         WRITE(0,*) 'CALIBRATION 2:',calibration
      if (calibration.eq.0.0) then 
         WRITE(0,*) 'CALIBRATION:',maxval
         calibration = maxval/255.0
      endif
c
c     Loop through image to generate scaled image_buffer
c
      do ic = 1,xres
         do ir = 1,yres

            temp_val = int(image(ic,ir)/calibration)
            if (temp_val.gt.255) temp_val = 255
            if (temp_val.lt.0) temp_val = 0
            image_buffer(ic,ir) = int(temp_val,SHORT)

c            WRITE(0,*) 'IMAGE:',ic,ir,image(ic,ir)
     
         end do
      end do  
c
c     Determine file name for image
c
      image_count = image_count+1
c
c...  Load base case name from environment variable CASENAME:
      CALL CaseName(filename,ierr)
c        
      len = lenstr(filename)
c
      call create_file_name(fname,filename(1:len),image_count,1)
c
c      len = lenstr(fname)
c      write(6,*) 'File:',fname(1:len) 
c
c     Write out actual JPG file - using quality value of 100 for highest
c     quality. This routine is written in C and accesses the publicly
c     available jpegsrc library. The code is based on the example.c
c     file shipped with that library. 
c

      len = lenstr(fname)
c
c     For the SUN implementation - the null termination is added in the 
c     C code since the Fortran compiler seemed to have trouble figuring
c     out what was required.
c
      call write_JPEG_file (fname,%VAL(len),image_buffer, 
     >                      %VAL(yres),%VAL(xres), %VAL(100),%VAL(1))
c
c
c---------------------------------------------------------------------- 
c
c     Free space assigned to image_buffer array
c
      deallocate(image_buffer)
c
c---------------------------------------------------------------------- 

      return 
      end 
c
c
c
      subroutine calc_image(image,xres,yres,position,direction,
     >                      upvec,rightvec,lookat,tmpplot,
     >                      maxval,minval,view_height,view_width,
     >                      minsteps,stepsize,mirror,reflect_opt)
      implicit none
      include 'params'
      integer xres,yres,minsteps,reflect_opt
      logical mirror
      real image(xres,yres),stepsize
      real*8 position(3),direction(3),upvec(3),rightvec(3),lookat(3)
      real*8 maxval,minval,view_height,view_width
c
      real tmpplot(maxnks,maxnrs)
c
c
c     CALC_IMAGE:
c 
c     This routine takes the defined camera location and the DIVIMP
c     array and calculates the 2D image by performing a 3D integration
c     through the volume of rotation of the magnetic grid about the Z
c     or vertical axis. 
c
c     The R coordinate is defined as R = sqrt(x**2+y**2)
c
c     The view of the camera is defined from the given vectors:
c
c     position - viewing position of the camera
c     
c     direction - defines distance along the Z-axis to the 
c             viewing frame defined by UP and RIGHT.
c     upvec - defines the vertical extent of the camera image at a
c             diatance of |direction| from the camera postion along
c             the Z-axis. 
c             (In terms of the camera view - this is Y-axis-camera)
c     rightvec - defines the horizontal extent of the camera image at a
c             diatance of |direction| from the camera postion along
c             the Z-axis. 
c             (In terms of the camera view - this is X-axis-camera)
c     lookat - defines a location in space at which the camera is 
c             looking - it lies along the viewing axis of the 
c             camera.
c     xres,yres - the horizontal and vertical pixel resolution of the 
c             camera.         
c     tmpplot - the 2D DIVIMP data to be integrated over
c     image - array allocated to store the resulting image
c
c
c     Local variables:
c
      integer res,ierr
c
c     Corner vectors and other related variables
c
      real*8 pt_ul(3),pt_ur(3),pt_ll(3),pt_lr(3),pt_cp(3)  
      real*8 fdist,hdist,vdist,ldist
      real*8 test_vect(3),lookvec(3)
      real*8 angle_rot_rightaxis,angle_rot_upaxis
c
c     Transformation matrix
c
      real*8 m(3,3)
c
c     Vector angle and view related values
c
c     Vector calculation functions
c      
      real*8 dotprod,vect_mag
      integer xprod
      external dotprod,xprod,vect_mag
c
c     Looping over view space related quantities
c
      real*8 tval,vn(3)
      real*8 vxp1(3),vxp2(3)
      real*8 dotv1,dotv2,vtheta1,vtheta2,dthetav1,dthetav2
      real*8 doth1,doth2,htheta1,htheta2,dthetah1,dthetah2
      real*8 toutv1,toutv2
      real*8 hxp(3),vh1(3),vh2(3)
      real*8 htheta,dthetah,doth,touth
      real contrib(maxnks,maxnrs) 
      integer ic,ir,in
c 
      integer rightaxis,upaxis,viewaxis
c
c     LOS related quantities
c
      real*8 wres,step_size
      integer nchords
c
c     Initialization 
c      
      wres = 1.0
      nchords = 1
      step_size = stepsize
      call dzero(m,9) 
c
c     Extract view information from input camera data (vectors)
c
c     The UP and RIGHT vectors specify the full height and width of the
c     view 
c
      fdist = sqrt(direction(1)**2+direction(2)**2+direction(3)**2)
      vdist = sqrt(upvec(1)**2+upvec(2)+upvec(3)**2)/2.0
      hdist = sqrt(rightvec(1)**2+rightvec(2)**2+rightvec(3)**2)/2.0
c
c     Define the axes of upvec and rightvec - axis 3 or Z is assumed in 
c     the code to be the vertical axis (POVRAY uses Y).
c
c     It is expected that rightvec(2) will be non-zero and upvec(3) will be
c     non-zero - placing the image plane in the Y-Z plane.
c
c     The code checks for this.
c
c     Initialize axis vectors
c
      upaxis = 0
      rightaxis = 0
      viewaxis = 0
c
c
c     UP axis 
c
      if (upvec(3).ne.0.0) then
         upaxis = 3
      elseif (upvec(2).ne.0.0) then
         upaxis = 2
      elseif (upvec(1).ne.0.0) then 
         upaxis = 1
      else
         upaxis = 0
         ierr   = 1
         write (0,'(a,3(1x,g12.5))')
     >          'UP axis undefined:',(upvec(in),in=1,3)  
         return 
      endif
c
c     RIGHT axis
c
      if (rightvec(3).ne.0.0) then
         rightaxis = 3
      elseif (rightvec(2).ne.0.0) then
         rightaxis = 2
      elseif (rightvec(1).ne.0.0) then 
         rightaxis = 1
      else
         rightaxis = 0
         ierr   = 1
         write (0,'(a,3(1x,g12.5))')
     >          'RIGHT axis undefined:',(rightvec(in),in=1,3)  
         return
      endif
c
c     Check that up and right are NOT the same
c
      if (upaxis.eq.rightaxis) then 
         ierr   = 1
         write (0,'(a,6(1x,g12.5))')
     >          'RIGHT and UP axis the same:',
     >          ((upvec(in),rightvec(in)),in=1,3)  
         return 
      endif
c
c     The VIEW axis is the remaining unused one
c 
      viewaxis = 1
c
      if (rightaxis.eq.viewaxis.or.upaxis.eq.viewaxis) then 
         viewaxis = 2
         if (rightaxis.eq.viewaxis.or.upaxis.eq.viewaxis) then 
            viewaxis = 3
         endif        
      endif        
c
c
c     The rotations are performed about the RIGHT axis first and then 
c     about the UP axis so that the camera is aligned with the 
c     lookat position.
c     
c     The UP axis defines the vertical axis for the image - in this code
c     it should always be the Z-axis (Y-axis in POVRAY format).
c
      view_height = 2.0d0*vdist
      view_width  = 2.0d0*hdist
c
c     Define the corner points of the base view relative to 
c     the observation postion. The view is defined along the POSITIVE
c     viewaxis.
c 
c     Upper Left - 
c     Upper Right -
c     Lower Left -
c     Lower Right -
c     Center Point -   
c
      pt_ul(rightaxis) = +hdist
      pt_ul(upaxis)    = +vdist
      pt_ul(viewaxis)  = +fdist
c
      pt_ur(rightaxis) = -hdist
      pt_ur(upaxis)    = +vdist
      pt_ur(viewaxis)  = +fdist
c
      pt_ll(rightaxis) = +hdist
      pt_ll(upaxis)    = -vdist
      pt_ll(viewaxis)  = +fdist
c
      pt_lr(rightaxis) = -hdist
      pt_lr(upaxis)    = -vdist
      pt_lr(viewaxis)  = +fdist
c       
      pt_cp(rightaxis) = 0 
      pt_cp(upaxis)    = 0
      pt_cp(viewaxis)  = +fdist
c
c     Calculate the rotation angles required:
c
c     First about RIGHT-axis then about the UP-axis
c
c     This can be generalized - The image plane is defined 
c     by two axes XY, YZ, or XZ. Rotation about the vertical or UP
c     axis is done second - this will keep the top of the 
c     camera frame horizontal relative to that axis. 
c
      do in = 1,3
         lookvec(in) = lookat(in) - position(in) 
      end do
c
      call calc_angles(angle_rot_rightaxis,angle_rot_upaxis,
     >                 lookvec,rightaxis,upaxis,viewaxis)
c
c     The axes rotation angles are the negative of the angles 
c     required to map the vectors above.
c
      angle_rot_rightaxis = -angle_rot_rightaxis
      angle_rot_upaxis = -angle_rot_upaxis
c
c     Calculate distance from camera position to lookat point - this
c     will be used to test the transformation matrix. 
c 
      ldist = vect_mag(lookvec)
      test_vect(rightaxis) = 0.0
      test_vect(upaxis)    = 0.0
      test_vect(viewaxis)  = ldist 
c
c      call calc_transform(m,angle_rot_rightaxis,rightaxis,0)           
c
c      call transform_vect(m,test_vect)
c      write(0,'(a,8(1x,g12.5))') 'TEST1:',angle_rot_rightaxis*raddeg,
c     >        angle_rot_upaxis*raddeg,
c     >        lookvec(1),lookvec(2),lookvec(3),
c     >        test_vect(1),test_vect(2),test_vect(3)
c      call calc_transform(m,angle_rot_upaxis,upaxis,0)          
c      call transform_vect(m,test_vect)
c      write(0,'(a,8(1x,g12.5))') 'TEST2:',angle_rot_rightaxis*raddeg,
c     >        angle_rot_upaxis*raddeg,
c     >        lookvec(1),lookvec(2),lookvec(3),
c     >        test_vect(1),test_vect(2),test_vect(3)
c
c
c     Calculate the coordinate transformation matrix that 
c     will map the original corner point vector to the 
c     actual corner point vectors for the view centered on the 
c     lookat position.
c
c     MUST do the RIGHT-axis rotation first - 
c     followed by the UP-axis - otherwise
c     the axes are also rotated - the angles calculated are for 
c     rotations about the original axes.
c
      call calc_transform(m,angle_rot_rightaxis,rightaxis,0)           
      call calc_transform(m,angle_rot_upaxis,upaxis,1)           
c
c     Transform the corner vectors of the imaged region from the initial 
c     definition to the actual view based on the "look at" location.
c
      call transform_vect(m,pt_ul)
      call transform_vect(m,pt_ur)
      call transform_vect(m,pt_ll)
      call transform_vect(m,pt_lr)
      call transform_vect(m,pt_cp)
      call transform_vect(m,test_vect)
c
c      write(6,'(a,8(1x,g12.5))') 'TEST:',angle_rot_rightaxis*raddeg,
c     >        angle_rot_upaxis*raddeg,
c     >        lookvec(1),lookvec(2),lookvec(3),
c     >        test_vect(1),test_vect(2),test_vect(3)
c
c      write(0,'(a,8(1x,g12.5))') 'TEST:',angle_rot_rightaxis*raddeg,
c     >        angle_rot_upaxis*raddeg,
c     >        lookvec(1),lookvec(2),lookvec(3),
c     >        test_vect(1),test_vect(2),test_vect(3)
c
c     Set up nested loops to generate a series of vectors across the
c     image line by line from the ul to ur and then moving down row by
c     row. Then call the integration routine to generate a value for 
c     each of these LOS.  
c
c
c     Normalize the corner vectors 
c
      call norm_vect(pt_ul)
      call norm_vect(pt_ur)
      call norm_vect(pt_ll)
      call norm_vect(pt_lr)
      call norm_vect(pt_cp)
c
c     Set up to create two fans of side vectors and then scan across these
c
c     Calculating each individual line involves solving for the vector
c     with the following constraints.
c
c     - If v1 and v2 are the extremes of the given fan.
c     - vp is the cross-product between v1 and v2.
c     - Theta is the angle between v1 and v2.
c     - vn is the vector at the given angle theta1 between v1 and v2. 
c
c     For a given angle theta1 between v1 and v2 (0<=theta1<=theta)
c     then three conditions must be satisfied by the vector at this angle
c
c     (1) v1.vn = cos(theta1)
c     (2) v2.vn = cos(theta-theta1)
c     (3) vp.vn = 0
c
c     This becomes a set of 3 equations in 3 unknowns which can be represented
c     in matrix notation such that:
c
c     M * vn = R
c
c     where M = [ v1x v1y v1z ]  vn = [ vnx ]   R = [ cos(theta)         ]
c               [ v2x v2y v2z ]       [ vny ]       [ cos(theta1-theta)  ]
c               [ vpx vpy vpz ]       [ vnz ]       [  0.0               ]
c
c     The solution to this is:  vn = M**-1 * R
c
c     The subroutine solv_vect takes the vectors v1,v2,vp and 
c     theta and theta1 and solves
c     for vn. Solv_vect uses straight gaussian elimination to 
c     solve the equations
c     without proceeding to calculating an inverse matrix. 
c      
c
      res =  xprod(pt_ul,pt_ll,vxp1)
      call norm_vect(vxp1)  
c
      res = xprod(pt_ur,pt_lr,vxp2)
      call norm_vect(vxp2)  
c
c     vertical 
c
      dotv1= dotprod(pt_ul,pt_ll)
      dotv2= dotprod(pt_ur,pt_lr) 
c
      vtheta1= acos(dotv1) 
      vtheta2= acos(dotv2) 
c
      dthetav1 = vtheta1/dble(yres)
      dthetav2 = vtheta2/dble(yres)
c
c     horizontal
c
      doth1 = dotprod(pt_ul,pt_ur)
      doth2 = dotprod(pt_ll,pt_lr)
c
      htheta1= acos(doth1) 
      htheta2= acos(doth2) 
c
      dthetah1 = htheta1/dble(xres)
      dthetah2 = htheta2/dble(xres)
c
c
c      write(6,'(a,8(1x,g12.5))') 'ANGLESV:',vtheta1*raddeg,
c     >                          vtheta2*raddeg,
c     >                          dthetav1*raddeg,dthetav2*raddeg,
c     >                          dotv1,dotv2
c      write(6,'(a,8(1x,g12.5))') 'ANGLESH:',htheta1*raddeg,
c     >                          htheta2*raddeg,
c     >                          dthetah1*raddeg,dthetah2*raddeg,
c     >                          doth1,doth2
c
c     Initialize min and max
c
      maxval = -HI
      minval = HI

c
c     for each viewing line, ii, calculate the central
c     directional cosine vector, vn, from the three conditions:
c     (1) v1*vn = cos (thn)
c     (2) v2*vn = cos (th3-thn)
c     (3) vp*vn = 0
c

      do ir = 1,yres
c
c         toutv1 = dthetav1 * (dble(ir)-1.0d0)     
c         toutv2 = dthetav2 * (dble(ir)-1.0d0)     
         toutv1 = dthetav1 * (dble(ir)-0.5d0)     
         toutv2 = dthetav2 * (dble(ir)-0.5d0)     
c
         call solv_vect(pt_ul,pt_ll,vxp1,toutv1,
     >                  vtheta1-toutv1,vh1,ierr)
         call solv_vect(pt_ur,pt_lr,vxp2,toutv2,
     >                  vtheta2-toutv2,vh2,ierr)
c
c        Calculate characteristics of horizontal scan
c
         doth   = dotprod(vh1,vh2)
         htheta = acos(doth)
         dthetah= htheta/dble(xres)
         res    = xprod(vh1,vh2,hxp)
c
c         write(6,'(a,3(1x,g12.5))') 'VH1:',vh1(1),vh1(2),vh1(3) 
c         write(6,'(a,3(1x,g12.5))') 'VH2:',vh2(1),vh2(2),vh2(3) 
c         write(6,'(a,3(1x,g12.5))') 'HXP:',hxp(1),hxp(2),hxp(3) 
c         write(6,'(a,5(1x,g12.5))') 'DOT:',doth,htheta*raddeg,
c     >      dotprod(vh1,vh2),dotprod(vh1,hxp),dotprod(vh2,hxp)
c         


c
c         write(0,'(a,i5,2f10.4)') 'Working on row:',ir,
c     >              htheta*raddeg,dthetah*raddeg
c
c        Loop over each horizontal row between the two defining end vectors
c         
         do ic = 1, xres
c
c            touth = dthetah * (dble(ic)-1.0d0)
c
            if (mirror) then 
               touth = dthetah * (dble(xres-ic+1)-0.5d0)
            else  
               touth = dthetah * (dble(ic)-0.5d0)
            endif 
c            
            call solv_vect(vh1,vh2,hxp,touth,htheta-touth,vn,ierr)
c            
c           Call the routine to perform the LOS integration
c
c           Calculate the values along each horizontal scan line of the image.
c           Code in LOS3D will take care of doing each row if the end vectors
c           are properly defined.
c
c
            call los3dint(tval,position,vn,wres,
     >                 tmpplot,nchords,step_size,minsteps,contrib,
     >                 reflect_opt)
c
            image(ic,ir) = tval
c
c            write(6,'(a,2i5,1x,g12.5)') 'LOS:',ir,ic,tval
c
            maxval = max(tval,maxval)
            minval = min(tval,minval)
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
      subroutine calc_angles(angle_rot_rightaxis,angle_rot_upaxis,
     >                 lookvec,rightaxis,upaxis,viewaxis)
      implicit none
c
      include 'params'
c
      integer rightaxis,upaxis,viewaxis 
      real*8 angle_rot_rightaxis,angle_rot_upaxis,temp
      real*8 lookvec(3)
      real*8 datan2c,vect_mag
      external datan2c,vect_mag
c
c     Calculate first rotation angle about rightaxis
c
c     The is the angle above or below the rightaxis-viewaxis plane
c     required to point in the direction of the lookat position. 
c
      temp = lookvec(upaxis)/vect_mag(lookvec) 
c
      angle_rot_rightaxis = asin(temp)
c
c      write(0,'(a,3i4,8(1x,g12.5))') 'AngleR1:',
c     >           upaxis,rightaxis,viewaxis,temp
c      write(0,'(a,6(1x,g12.5))') 'AngleR2:',
c     >           angle_rot_rightaxis*raddeg,vect_mag(lookvec),
c     >           lookvec(1),lookvec(2),lookvec(3)
c
c
c     Calculate second rotation angle about upaxis 
c
      if (lookvec(rightaxis).eq.0.0.and.lookvec(viewaxis).eq.0.0) then 
         angle_rot_upaxis = 0.0
      else
c
c        View is always along the positive view axis to need to find
c        rotation angle from the positive view axis to the lookat
c        direction in the rightaxis-viewaxis plane - this gives the
c        rotation angle around the upaxis.
c

         angle_rot_upaxis = 
     >             datan2c(lookvec(rightaxis),lookvec(viewaxis))
c
c         write(0,'(a,3(1x,g12.5))') 'AngleU:',
c     >           angle_rot_upaxis*raddeg
c
      endif 
c
      return
      end
c
c
c
      subroutine solv_vect(v1,v2,vp,theta,theta1,vn,ierr)
      implicit none 
c
      include 'params'
c
      real*8 v1(3),v2(3),vp(3),theta,theta1,vn(3) 
      integer ierr
c
c     Calculates the solution to 3 equations in 3 unknowns using 
c     Gaussian elimination.
c
      real*8 m(3,4)
      real*8 fact,d1,d2,d3,dotprod
      external dotprod
      integer ir,ic 
c
      ierr = 0 
c
      do ic = 1,3
         m(1,ic) = v1(ic) 
         m(2,ic) = v2(ic) 
         m(3,ic) = vp(ic) 
      end do
c
      m(1,4) = cos(theta)
      m(2,4) = cos(theta1)
      m(3,4) = 0.0
c 
c     Put largest element in first place in first row.
c
      do ir = 2,3
c
         if (abs(m(ir,1)).gt.abs(m(1,1))) then
            call switchrow(1,ir,m,3,4)
         endif
c
      end do
c         
c     Now proceed with elimination and back substitution
c
      if (m(1,1).eq.0.0) then 
         ierr = 1
         write(0,*) 'SOLV_VECT: Matric solution problem'
         return
      endif
c
c     Eliminate first column
c
      do ir = 2,3
         fact = m(ir,1)/m(1,1)
         call subtractrow(1,ir,fact,m,3,4)
      end do
c 
c     Put largest element in remaining rows in second place
c
      if (abs(m(3,2)).gt.abs(m(2,2))) then 
         call switchrow(2,3,m,3,4)
      endif
c
c     Check for zero problem
c
      if (m(2,2).eq.0.0) then 
         ierr = 1
         write(0,*) 'SOLV_VECT: Matric solution problem'
         return
      endif
c
c     Eliminate second element of last row 
c
      fact = m(3,2)/m(2,2)
      call subtractrow(2,3,fact,m,3,4)
c
c     Subsitute to solve for vn
c
      if (m(3,3).eq.0.0) then 
         ierr = 1
         write(0,*) 'SOLV_VECT: Matric solution problem'
         return
      endif
c
c     Solve for components
c
      vn(3) = m(3,4) / m(3,3) 
      vn(2) = (m(2,4) - m(2,3) * vn(3)) / m(2,2)
      vn(1) = (m(1,4) - m(1,3) * vn(3) - m(1,2) * vn(2)) / m(1,1)
c
c      write(6,'(a,12(1x,g12.5))') 'M-D-1:',(m(1,ic),ic=1,4)
c      write(6,'(a,12(1x,g12.5))') 'M-D-2:',(m(2,ic),ic=1,4)
c      write(6,'(a,12(1x,g12.5))') 'M-D-3:',(m(3,ic),ic=1,4)
c      write(6,'(a,12(1x,g12.5))') 'M-DBG:',(vn(ic),ic=1,3)
c
c     Verify
c
c      d1 = dotprod(v1,vn)
c      d2 = dotprod(v2,vn)
c      d3 = dotprod(vp,vn)
c
c      write(6,'(a,10(1x,g12.5))') 'VERIFY:',
c     >         d1,cos(theta),d1-cos(theta),
c     >         d2,cos(theta1),d2-cos(theta1),
c     >         d3,theta*raddeg,theta1*raddeg
c
c
      return
      end
c
c      
c
      subroutine subtractrow(row1,row2,fact,m,dimrow,dimcol)
      implicit none
      integer row1,row2,dimrow,dimcol
      real*8  m(dimrow,dimcol),fact
c
c     Subtract row1 from row2 - element by element 
c  
      integer ic
c
c      write(6,'(a,2i4,4(1x,g12.5))') 
c     >              'SUR1B:',row1,dimcol,(m(row1,ic),ic=1,dimcol)
c      write(6,'(a,2i4,4(1x,g12.5))') 
c     >              'SUR2B:',row2,dimcol,(m(row2,ic),ic=1,dimcol)
c

      do ic = 1,dimcol
         m(row2,ic) = m(row2,ic) - fact * m(row1,ic)
      end do
c
c      write(6,'(a,2i4,4(1x,g12.5))') 
c     >              'SUR1A:',row1,dimcol,(m(row1,ic),ic=1,dimcol)
c      write(6,'(a,2i4,4(1x,g12.5))') 
c     >              'SUR2A:',row2,dimcol,(m(row2,ic),ic=1,dimcol)
c


c
      return 
      end 
c
c
c

      subroutine switchrow(row1,row2,m,dimrow,dimcol)
      implicit none
      integer row1,row2,dimrow,dimcol
      real*8  m(dimrow,dimcol)
c
c     Exchange two rows of the matrix
c       
      real*8 temp
      integer ic
c
c      write(6,'(a,2i4,4(1x,g12.5))') 
c     >              'SWR1B:',row1,dimcol,(m(row1,ic),ic=1,dimcol)
c      write(6,'(a,2i4,4(1x,g12.5))') 
c     >              'SWR2B:',row2,dimcol,(m(row2,ic),ic=1,dimcol)
c
      do ic = 1,dimcol
         temp       = m(row1,ic)
         m(row1,ic) = m(row2,ic)
         m(row2,ic) = temp
      end do
c
c      write(6,'(a,2i4,4(1x,g12.5))') 
c     >              'SWR1A:',row1,dimcol,(m(row1,ic),ic=1,dimcol)
c      write(6,'(a,2i4,4(1x,g12.5))') 
c     >              'SWR2A:',row2,dimcol,(m(row2,ic),ic=1,dimcol)
c
      return 
      end 
c
c
c
      subroutine norm_vect(vect)
      implicit none
      real*8 vect(3)
c
c     Normalize the vector components 
c
      real*8 norm      
c
      norm = sqrt(vect(1)**2+vect(2)**2+vect(3)**2)
      if (norm.ne.0.0) then 
         vect(1) = vect(1) / norm
         vect(2) = vect(2) / norm
         vect(3) = vect(3) / norm
      endif
c
      return
      end
c
c
c
      subroutine transform_vect(m,a)
      implicit none
      real*8 m(3,3),a(3)
c
c     Transform_VECT: rotate the vector to get its value in the new
c                     coordinate system - in which the rotation to map 
c                     the camera view will be entirely around the rotated
c                     X-axis.
c
      integer i,j
      real*8 tmpa(3) 
c
c      write(6,'(a,3(1x,g12.5))') 'START VECT:',(a(i),i=1,3)
c
      do i = 1,3 
c      
         tmpa(i) = m(i,1)*a(1) + m(i,2)*a(2)+m(i,3)*a(3)
c
      end do
c
c     Assign result to a
c
      do i = 1,3
c
         a(i) = tmpa(i)
c
      end do 
c
c      write(6,'(a,3(1x,g12.5))') 'END VECT  :',(a(i),i=1,3)
c
c      write(6,'(a,12(1x,g12.5))') 'M-TRANS:',
c     >                     ((m(i,j),j=1,3),i=1,3)
c
      return
      end
c
c
c
      real*8 function vect_angle(a,b)
      implicit none
      real*8 a(3),b(3)
      include 'params'
c
c     VECT_ANGLE: Calculates the angle between two three-vectors
c                 based on the following formula. 
c    
c                 The value returned is also relative to the rotation 
c                 axis if this axis corresponds to a coordinate axis for
c                 the points. 
c
c
c     A dot B = |A||B|cos(AB)
c
c     A X B = |A||B|sin(AB) - but also the projection along 
c                             the third axis is greater than zero 
c                             for a positive rotation and less than 
c                             zero for a negative one. 
c
c
c     Angle AB = cos -1 ( A dot B /  (|A||B|) )
c
      real*8 c(3)
      real*8 amag,bmag,vect_mag
      integer axis
      integer xprod
      external xprod,vect_mag
c
      amag = vect_mag(a)
      bmag = vect_mag(b)
c
      if (amag.ne.0.0.and.bmag.ne.0.0) then 
         vect_angle = acos( (a(1)*b(1) + a(2)*b(2) +a(3)*b(3))
     >                       /  (amag*bmag))
      else
         vect_angle = 0.0d0
      endif 
c
c     Calculate the components of the cross-product.
c  
      axis =  xprod(a,b,c)
c
      if (axis.gt.0) then 
         if (c(axis).lt.0.0) vect_angle = -vect_angle
      endif
c
c      write(6,'(a,i5,3(1x,g15.7))') 'VECTA:',axis,vect_angle*raddeg
c
      return
      end   
c
c
c
      real*8 function dotprod(a,b) 
      implicit none
      real*8 a(3),b(3)
c
c     Calculate dot product of two 3-vectors
c
      dotprod = a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
c      
      return
      end
c
c
c
      integer function xprod(a,b,c)
      implicit none
      real*8 a(3),b(3),c(3)
      real*8,parameter::eps=1.0d-10
      integer coord
c
      c(1) = a(2)*b(3)-a(3)*b(2)
      c(2) = a(3)*b(1)-a(1)*b(3)
      c(3) = a(1)*b(2)-a(2)*b(1)
c
      if (c(1).ne.0.0.and.abs(c(2)).lt.eps.and.abs(c(3)).lt.eps) then 
         xprod = 1
      elseif(c(2).ne.0.0.and.abs(c(1)).lt.eps.and.abs(c(3)).lt.eps)then 
         xprod = 2
      elseif(c(3).ne.0.0.and.abs(c(1)).lt.eps.and.abs(c(2)).lt.eps)then 
         xprod = 3 
      else
         xprod = 0 
      endif
c
c      write(6,'(a,i5,3(1x,g15.7))') 'A    :',xprod,a(1),a(2),a(3) 
c      write(6,'(a,i5,3(1x,g15.7))') 'B    :',xprod,b(1),b(2),b(3) 
c      write(6,'(a,i5,3(1x,g15.7))') 'XPROD:',xprod,c(1),c(2),c(3) 
c
      return 
      end 
c
c
c
      subroutine calc_transform(m,angle,axis,init)
      implicit none
c
      real*8 m(3,3)
      real*8 angle
      integer axis,init
c
c     CALC_TRANSFORM: calculate the transformation matrix resulting from a 
c                     rotation of magnitude angle about the given axis. 
c
c      
      real*8 trans(3,3)
      real*8 tempm(3,3)
      integer i,j
c
c     Initialize to identity matrix if initialization option on
c
      if (init.eq.0) then 
         call dzero(m,9)
         m(1,1) = 1.0 
         m(2,2) = 1.0 
         m(3,3) = 1.0 
      endif
c
      call dzero(trans,9)
      call dzero(tempm,9)
c
c     Calculate transformation matrix
c
c     Rotation about first-axis 
c
      if (axis.eq.1) then 
         trans(1,1) =  1.0
         trans(2,2) =  cos(angle)
         trans(3,3) =  cos(angle)
         trans(2,3) =  sin(angle) 
         trans(3,2) = -sin(angle) 
c
c     Rotation about second-axis 
c
      elseif (axis.eq.2) then 
         trans(1,1) =  cos(angle)
         trans(2,2) =  1.0
         trans(3,3) =  cos(angle)
         trans(1,3) =  sin(angle) 
         trans(3,1) = -sin(angle) 
c
c     Rotation about third-axis
c
      elseif (axis.eq.3) then 
         trans(1,1) =  cos(angle)
         trans(2,2) =  cos(angle)
         trans(3,3) =  1.0
         trans(1,2) =  sin(angle) 
         trans(2,1) = -sin(angle) 
      endif 
c
c      
c     Multiply transformation and matrix - pre-multiply 
c
      do i=1,3
         do j = 1,3
c                         
            tempm(i,j) =   trans(i,1) * m(1,j) 
     >                   + trans(i,2) * m(2,j) 
     >                   + trans(i,3) * m(3,j)
c
         end do 
      end do
c  
c     Assign resulting transformation to m
c      
c     m = tempm
c
      do i = 1,3
         do j = 1,3
            m(i,j) = tempm(i,j)
         end do
      end do 
c
c      write(0,'(a,i4,12(1x,g12.5))') 'M-CALC:',axis,
c     >                     ((m(i,j),j=1,3),i=1,3)
c
      return
      end
c
c
c
      subroutine plot_3Dlos(iselect,istate,npts,nlines,
     >                  iexpt,iaxis,iavg,ifact,optval,
     >                  iopt,job,title,table,avs,navs,nplots,
     >                  iplot,nizs,ierr) 
      implicit none
c
      include 'params'
      include 'cgeom'
c
      integer iselect,istate,npts,nlines,iexpt,iaxis,ifact,ierr 
      integer iavg
      real optval 
      character*(*) job,title,table
      integer navs,iopt,nplots,iplot,nizs
      real avs(0:navs)
c
c     PLOT_3DLOS:
c
c     This routine performs a generalized 2D LOS plot through 3D.
c     Each LOS is defined by a start point and either a direction vector
c     or a point along the line of sight. Each LOS also includes a defintion
c     of a viewing region - either a cone defined by dtheta or a 
c     square/rectangular viewing area. These different shapes are relevant when
c     calculating the experimental data in units of solid angle. 
c  
c     Local variables
c
      integer maxip
      parameter(maxip=100)
      integer nvr,nvz,nvp,nr,nz,np,ndt,nd,nrt
      real rvertex(maxip),zvertex(maxip),phivertex(maxip)
      real rvals(maxip),zvals(maxip),phivals(maxip),dtheta(maxip),
     >     dists(maxip),rtan(maxip)
C
      real tmp_store
c
c      real r,z,thet,dthet,dist
c
      real tmpplot(maxnks,maxnrs)
      real mfact
c
c     Plot labels   
c
      character*36 blabs(2),xlab,ylab
      character*44 ref,plane,anly,nview,smooth
      character*32 datatitle
c
      INTEGER IGNORS(MAXNGS),ngs
c
      REAL TOUTS(MAXTHE),TWIDS(MAXTHE),TVALS(MAXTHE,MAXNGS)
c      real tmptvals(0:maxthe-1)
      REAL THEMIN,THEMAX
      integer itec,ismoth
c
      integer len,lenstr,ii,ij,ik,in,ir
c
c     Variables for 2D experimental data - iexpt = -1
c
      character*256 graph
c      character*256 filename
c      integer  maxix,maxiy,nix,niy,ifnopt
c      parameter(maxix=100,maxiy=100)
c      real expt_array(maxix,maxiy),raxis(maxix),zaxis(maxiy)
c      real div2D_array(maxix,maxiy)
c      real expt_rmin,expt_rmax,expt_zmin,expt_zmax,expt_dr,expt_dz
c
c     LOS related quantities
c
      integer reflect_opt
      real*8 wres,step_size,dthet,tval
      real*8 obs(3),view(3),startvn(3),vn(3)
      real*8 dotprod
      real contrib(maxnks,maxnrs) 
      external dotprod 
      integer nchords
c
c     Contour plot options 
c
      real xcen,ycen,xnear,ynear
      real xmin,xmax,ymin,ymax 
      real uconts(maxpts)
      integer icntr,ncntr
      real minscale,maxscale,maxval
c
      real xouts(maxgxs),youts(maxgys)
      real sum_contrib 
c
c     Contour kludge 
c
c     Initialization
c
      do ii = 1,maxngs 
         ignors(ii) = 1
      end do
c
      step_size = optval 
      nchords = 1
      wres = 1.0
      xmin = rmin
      xmax = rmax 
      ymin = zmin
      ymax = zmax 
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
      smooth = ' '
c
c     Based on the values of iselect and istate - load up the required 
c     2D array for the LOS plot 
c
c     Scaling by absfac (if required) is performed in the 
c     load_divdata_array routine. 
c
      call rzero(tmpplot,maxnks*maxnrs)
c
      call load_divdata_array(tmpplot,iselect,istate,1,
     >                         ylab,blabs(1),ref,nizs,ierr)
c
c     Check return code from data load and exit if non-zero   
c
      if (ierr.ne.0) return 
c
c
c     Now that the data has been successfully loaded - load the 
c     data describing the various lines of sight. 
c
c     If an error occurs in reading the input data then exit immediately.
c
c     Load Rvertex-values 
c
      call RDG_REAL_ARRAY(GRAPH,rvertex,maxip,nvr,ierr)
      if (ierr.ne.0) return
c
c     Load Zvertex-values 
c
      call RDG_REAL_ARRAY(GRAPH,zvertex,maxip,nvz,ierr)
      if (ierr.ne.0) return
c
c     Load PHIvertex-values 
c
      call RDG_REAL_ARRAY(GRAPH,phivertex,maxip,nvp,ierr)
      if (ierr.ne.0) return
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
c     Load PHI-values 
c
      call RDG_REAL_ARRAY(GRAPH,phivals,maxip,np,ierr)
      if (ierr.ne.0) return
c
c     Load Delta-theta-values 
c
      call RDG_REAL_ARRAY(GRAPH,dtheta,maxip,ndt,ierr)
      if (ierr.ne.0) return
c
c     Load length-values 
c
c      call RDG_REAL_ARRAY(GRAPH,dists,maxip,nd,ierr)
c      if (ierr.ne.0) return
c
c     Read in R tangent values for plotting IF that option is
c     specified.
c
      if (iaxis.eq.2) then 
         call RDG_REAL_ARRAY(GRAPH,rtan,maxip,nrt,ierr)
         if (ierr.ne.0) return
      endif 
c
c     Read in contour options if contour plot of LOS contributions
c     was specified 
c
      if (iavg.gt.0.and.iavg.le.npts+1) then 
c
         call rdg_contopts(graph,icntr,ncntr,uconts,maxpts,
     >                  xcen,ycen,xnear,ynear,ierr)
c
c        Set up plotting ranges
c
         xmin = xcen - xnear
         xmax = xcen + xnear
         ymin = ycen - ynear
         ymax = ycen + ynear

c
c        Set up colours based on contour options  
c	 
         if (icntr.eq.0.or.icntr.eq.2) then
            call setup_col(ncntr+1,2)
         elseif (icntr.eq.1) then 
            ncntr = 10
            call setup_col(ncntr+1,2)
         elseif (icntr.eq.3.or.icntr.eq.4) then
            ncntr = ncntr
            call setup_col(ncntr+1,2)
         endif
c
      endif
c
c     Make sure scaling option is in range 
c
      IF (ifact.LT.0 .OR. ifact.GT.3) ifact = 0
c
c     All data to calculate LOS plot is now available - loop through 
c     calculating values. 
c
c
c
C     INTEGRATE
C
c     Perform a series of single 3D LOS integrations  
c
c     Use the global npts that was specified in the first line - if
c     there are less that npts values specified for any line - the 
c     last value on each input line will be used for all remaining
c     lines of sight. 
c         
c
      do in = 1,npts 
c
c        Calculate the X,Y,Z values for vertex.
c                   
         ii = min(in,nvr)
         ij = min(in,nvp)
         ik = min(in,nvz)
         obs(1) = rvertex(ii) * cos(phivertex(ij)*degrad)
         obs(2) = rvertex(ii) * sin(phivertex(ij)*degrad)
         obs(3) = zvertex(ik)
c
c        Calculate the X,Y,Z values for strike/viewing point
c
         ii = min(in,nr)
         ij = min(in,np)
         ik = min(in,nz)
         view(1) = rvals(ii) * cos(phivals(ij)*degrad)
         view(2) = rvals(ii) * sin(phivals(ij)*degrad)
         view(3) = zvals(ik)
c
c        Calculate viewing vector
c
         do ii = 1,3
            vn(ii) = view(ii) - obs(ii)
         end do
c
c        Normalize the view vector
c
         call norm_vect(vn)
c          
c        Save first vector for angle calculations in axis option 1
c
         if (in.eq.1) then 
            do ii = 1,3
               startvn(ii) = vn(ii)
            end do
         endif
c
c        DTheta
c
         ii = min(in,ndt)
         dthet = dtheta(ii)  
c
c        Dist
c
c         ii = min(in,nd)
c         dist = dists(ii)  
c
c        Set scale factor based on viewing width 
c
         if (ifact.eq.1) then 
c
             MFACT = DTHET / (2.0 * PI)
c
         endif  
c
c
c        Call the routine to perform the LOS integration
c
c        Calculate the values along each horizontal scan line of the image.
c        Code in LOS3D will take care of doing each row if the end vectors
c        are properly defined.
c
c        Set minrun to zero for now for these integrations ...
c
c        Turn off reflections for now with this diagnostic
c
         reflect_opt = 0
c
         call los3dint(tval,obs,vn,wres,
     >                 tmpplot,nchords,step_size,0,contrib,
     >                 reflect_opt)
c
         tvals(in,1) = tval
c
c        Produce contour plot of contributions if requested 
c 
c        Contour of contributions to specific line
c 
         if (in.eq.iavg.or.iavg.eq.(npts+1)) then  
c
c           Rescale contrib array to 0.0 to 1.0 
c
            sum_contrib = 0.0 
            do ir = 1,maxnrs
               do ik = 1,maxnks
                  sum_contrib = sum_contrib + contrib(ik,ir)
               enddo 
            enddo
c
c           Normalize contributions and plot
c
            if (sum_contrib.gt.0.0) then 
c
               maxval = 0.0
c
               do ir = 1,maxnrs
                  do ik = 1,maxnks
                     contrib(ik,ir) = contrib(ik,ir)/sum_contrib
                     maxval = max(maxval,contrib(ik,ir))
                  enddo 
               enddo

               minscale = 0.0
               maxscale = maxval 

               XLAB   = '   R  (M)'
               YLAB   = '   Z  (M)'
               NGS    = ncntr
               write(REF,'(''CONTRIB: LOS='',i4)') in
               write(PLANE,'(''RSTRIKE = '',f8.3)') rvals(in)
               WRITE (IPLOT,9012) NPLOTS,REF
c
               CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XMIN,XMAX,
     >              YMIN,YMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
               CALL CONTOUR (1,ncntr,contrib,1,1,1,1.0,0.0,1.0,
     >                XOUTS,1,maxgxs,YOUTS,1,maxgys,
     >                XMIN,XMAX,YMIN,YMAX,
     >                ncntr,uconts,icntr,minscale,maxscale)

            endif
c
         endif
c
c        Calculate axis values
c
c        Theta 
c
         if (iaxis.eq.1) then  
c
            touts(in) = acos(dotprod(startvn,vn))
            twids(in) = dthet
c
c        R-tangential - MUST be loaded as extra data
c
         elseif (iaxis.eq.2) then 
c   
            ii = min(in,nrt)
            touts(in) = rtan(ii)
            twids(in) = 1.0 
c
c        R-strike
c
         elseif (iaxis.eq.3) then 
c   
            ii = min(in,nr)
            touts(in) = rvals(ii)
            twids(in) = 1.0 
c
c        Z-strike
c
         elseif (iaxis.eq.4) then 
c   
            ii = min(in,nr)
            touts(in) = zvals(ii)
            twids(in) = 1.0 
c
c        Index number
c
         elseif (iaxis.eq.5) then 
c   
            touts(in) = in
            twids(in) = 1.0 
c
         endif         
c
c        End of LOS calculation loop 
c
      end do
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
c      if (nd.eq.1.and.dists(1).gt.0.0) then
c
c         write(PLANE,'(a,f6.3,a)')
c     >   'INTEGRATION OVER ',dists(1),' (M) ONLY'
c
c      endif   


c
c     Load plot labels
c

c
c     Set LOS scaling factor and note value for plot
c
      IF (ifact.EQ.0) THEN
         mfact = 1.0
         WRITE(ANLY,'(''NO SCALE FACTOR APPLIED'')')
      ELSEIF (ifact.EQ.1) THEN
         WRITE(ANLY,'(''SCALE FACTOR = DTHETA / (2*PI)'')')
      ELSEIF (ifact.EQ.2) THEN
         MFACT = 1.0 / (2.0 * PI)
         WRITE(ANLY,'(''SCALE FACTOR = 1 / (2*PI)'') ')
      ELSEIF (ifact.EQ.3) THEN
         MFACT = 1.0 / ( 4.0 * PI)
         WRITE(ANLY,'(''SCALE FACTOR = 1 / (4*PI)'')')
      ENDIF
c
C     X-axis labels
C
c     Determined by value of iaxis
c
c     1 = Plot versus angle - calculated from first LOS given
c     2 = Plot versus R-tangent - MUST be specified in input
c     3 = plot versus R-strike
c     4 = plot versus Z-strike
c     5 = Index number
c 
c
      if (iaxis.eq.1) then 
c
         XLAB = 'THETA FROM FIRST LOS (DEGREES)'
c
      elseif (iaxis.eq.2) then 
c
         write(XLAB,
     >       '(''R-Tangent (M)'')')
c
      elseif (iaxis.eq.3) then 
c
         write(XLAB,
     >       '(''R-Strike (M)'')')
c
      elseif (iaxis.eq.4) then 
c
         write(XLAB,
     >       '(''Z-Strike (M)'')')
c
      elseif (iaxis.eq.5) then 
c
         XLAB = 'Index Number'
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
      if (iexpt.gt.0) then 
c
          ngs = 2
c
          call calc_expt(iexpt,touts,tvals,maxthe,npts,
     >                   themin,themax,maxngs,ngs,datatitle)

          BLABS(2) = 'EXPT '//DATATITLE
c
      endif
c
C
C     DON'T PLOT SINGLE LINES OF SIGHT - Just write single numbers to file
C
      IF (Npts.GT.1) THEN
c
          do in = 1,npts
             write(6,*) 'TVALS3D:',in,touts(in),tvals(in,1),
     >                                  tvals(in,2)
          end do
c
          itec = 1
          ismoth = 99          
c
          CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NPTS,ANLY,NGS,
     >              ISMOTH,THEMIN,THEMAX,-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >              JOB,TITLE,XLAB,YLAB,BLABS,REF,NVIEW,PLANE,
     >              TABLE,IOPT,2,1.0,0)
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
c     Return from generalized 3D LOS routine
c
      return
c
c     Format statements 
c

 9012 FORMAT(1X,'PLOT',I3,4X,A)
c
      end


