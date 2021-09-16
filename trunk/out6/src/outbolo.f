      subroutine out_bolo_diiid(power_loss,infile,iseld,graph,iopt,
     >                         iflag,
     >                         job,title,table,avs,navs,iplot,nplots)
      use mod_params
      use mod_cgeom
      implicit none
c
c     include 'params' 
c     include 'cgeom'
c
      integer iopt,iseld,iplot,nplots,iflag
      real power_loss(maxnks,maxnrs)
      character*(*) infile,graph,job,title,table
      integer navs
      real avs(0:navs)
c
c     OUT_BOLO_DIIID: This routine generates a bolometry plot
c     based on the pre-calculated DIIID channel component 
c     matrices. These matrices contain the contributions to each 
c     bolometer channel from each cell of a square grid that 
c     overlays the vessel cross-section. This region covers the 
c     area  R,Z = [0.84:2.54,0.0:3.2] in the shifted DIIID
c     coordinate space. There are 165 elements in the R direction and 
c     325 in the Z. 
c
C     NOTE: The array is stored in a disk file - to save memory - we 
c           do not allocate an array with the required number of 
c           elements since this would be 325x165x71 or about 15 Mb for 
c           a real array - instead the data is simply processed as it is
c           read.
C
c
c     There are 71 channels in the file - we only calculate and plot
c     the first 48.
c       
c     Set plotting range - assume that if zmin is less than zero that 
c     the grid has been shifted to normal DIIID coordinates - i.e. down
c     1.6 m - otherwise use the numbers as is. 
c
c     According to Bill Meyer the channel centres start and end at the 
c     limits of the given range - as a result the end bins actually would 
c     go one half bin width beyond the defined region.
c
c
c     Local Variables 
c
      real rbnd1,rbnd2,zbnd1,zbnd2,z_shift
      integer nchan,nrbin,nzbin
c
      character*36 blabs(2),xlab,ylab
      character*44 ref,plane,anly,nview
      character*32 datatitle
c
      INTEGER IGNORS(MAXNGS),ngs

      REAL TOUTS(MAXTHE),TWIDS(MAXTHE),TVALS(MAXTHE,MAXNGS)
      real tmptvals(0:maxthe-1)
      REAL THEMIN,THEMAX
c
      parameter(rbnd1=0.84,rbnd2=2.54,zbnd1=0.0,zbnd2=3.2,z_shift=-1.6)
      parameter(nchan=48,nrbin=165,nzbin=325)
c 
      real r_min,r_max,z_min,z_max,drbin,dzbin
      real raxis(0:nrbin-1),zaxis(0:nzbin-1)
c
      real getvalrz,r,z,coeff
      external getvalrz
      integer ios,len,lenstr,ix,iy,ichan,ii
      character*200 inline
c
c     Local index arrays 
c
      real powvals(0:nrbin-1,0:nzbin-1)

c
c     Initialization
c
      do ii = 1,maxngs 
         ignors(ii) = 1
      end do
c
      ngs =1 
c
c     Set boundaries for matrix region 
c
      r_min = rbnd1
      r_max = rbnd2
c
      if (zmin.lt.0.0) then 
         z_min = zbnd1 + z_shift         
         z_max = zbnd2 + z_shift         
      else
         z_min = zbnd1         
         z_max = zbnd2 
      endif
c
      drbin = (r_max-r_min)/(2.0*(nrbin-1)) 
      dzbin = (z_max-z_min)/(2.0*(nzbin-1)) 
c
c     Pre-calculate the power array to speed processing
c
c
c     Calculate axis arrays
c     
      do ix = 0,nrbin-1
         raxis(ix) = ix * 2.0 * drbin + r_min
      end do
c     
      do iy = 0,nzbin-1
         zaxis(iy) = iy * 2.0 * dzbin + z_min
      end do  
c
      call create_2Ddata(power_loss,nrbin,nzbin,nrbin,nzbin,
     >                   powvals,raxis,zaxis,iflag)
c
c
c
C  CREATE THETA VECTOR
C
      DO II = 0, nchan-1
         TOUTS(II+1) = ii
         TWIDS(II+1) = 1.0
      ENDDO
c
      themin = touts(1)
      THEMAX = TOUTS(nchan-1)
C
C  LABELS AND SCALE FACTORS
C
      YLAB = 'PRAD (BOLO) (W)'
      XLAB = 'CHANNEL #'
      BLABS(1) = 'BOLO DIIID BOLOMETRY'
      REF  = GRAPH(5:41)
c
      NPLOTS = NPLOTS + 1
      WRITE (IPLOT,9012) NPLOTS,REF
c
      call rzero(tvals,maxthe*maxngs)
      call rzero(tmptvals,maxthe) 
c
c     Open file 
c
      open(tmpunit,FILE=infile,IOSTAT=ios,STATUS ='old')
c
      if (ios.ne.0) then 
c
         write(6,*) 'Error opening bolometry matrix file: ',ios,infile 
         write(0,*) 'Error opening bolometry matrix file: ',ios,infile 
         return
c
      endif
c
c     Read in first line containing number of lines in file
c
      read(tmpunit,'(a200)',IOSTAT=ios)  inline
c
      do while (ios.eq.0) 
c
c        read in each line
c 
         read(tmpunit,'(a200)',IOSTAT=ios,END=100)  inline
c
         read(inline,*,IOSTAT=ios) ix,iy,ichan,coeff
c
         if (ichan.lt.nchan.and.ios.eq.0) then
c
c           ix is in the range 0 to nrbin-1
c           iy is in the range 0 to nzbin-1
c           drbin,dzbin are bin half widths the extra drbin
c           shifts to bin centers
c
            tmptvals(ichan) = tmptvals(ichan)
     >                + coeff * powvals(ix,iy)
c
         endif
c      
      end do 
c
 100  write(6,*) 'End of matrix file:',ios
c
c     Close input file
c
      close(tmpunit)
c
c     Copy tmptvals
c
      do ii = 1,nchan
         tvals(ii,1) = tmptvals(ii-1)
      end do 
c
c     Load experimental data if any  
c
      if (iseld .gt.0) then 
c
         ngs = 2
c
         call calc_expt(iseld,touts,tvals,maxthe,nchan,
     >                 themin,themax,maxngs,ngs,datatitle)
c
         BLABS(2) = 'EXPT'//DATATITLE

      endif
c
      CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,nchan,ANLY,NGS,
     >               99,THEMIN,THEMAX,0.0,HI,IGNORS,0,AVS,NAVS,
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
      real function getvalrz(data_arr,r,z)
      use mod_params
      implicit none
c     include 'params'
      real r,z
      real data_arr(maxnks,maxnrs) 
c
c     Extract the value of array data_arr at the position r,z
c        
      integer ik,ir
      data ik,ir /0,0/
      logical newinj,outofgrid
c              
c     Initialization
c      
c      ik = 0
c      ir = 0
c
      getvalrz = 0.0
      newinj = .false.
      outofgrid = .false.
c 
      call gridpos(ik,ir,r,z,newinj,outofgrid)
c 
      if (outofgrid) then 
         getvalrz = 0.0
      else 
c
         getvalrz = data_arr(ik,ir)
c
      endif
c
      return
      end 




