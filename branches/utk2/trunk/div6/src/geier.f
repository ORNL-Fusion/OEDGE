      subroutine pointdist(x1, x2, y1, y2, xp, yp, dist, xi,yi)
      implicit none
      real x1,y1,x2,y2,xp,yp,xi,yi,dist
c      
c Calculates distance between point (xp, yp) and a straight line given by
c (x1, x2, y1, y2) as well as the intersection between the line and a 
c normal though the point 
c
c
c     Local variables
c
      real vAx,vAy,vBx,vBy,magA,factAB,dotAB
c
      vAx = x2 -x1
      vAy = y2 -y1
c
      vBx = xp-x1
      vBy = yp-y1
c
c     A dot B = |A||B| cos(theta)
c     The projection of B along A =  |B| cos(theta)
c
c     |B| cos(theta) = A dot B / |A|
c
c     The normalized fraction of the distance from x2 to x1
c     is then |B| cos(theta) / |A| - thus the factor in the 
c     following is  A dot B / |A|**2
c
c     xi = x1 + A dot B / |A|**2 (x2-x1)
c     yi = y1 + A dot B / |A|**2 (y2-y1)
c 
c     dist = sqrt((xp-xi)**2+(yp-yi)**2) 
c
      magA = sqrt(vAx**2+vAy**2)  
c
      dotAB = vAx*vBx + vAy * vBy
c 
      factAB = dotAB / magA**2
c
      xi = x1 + factAB * (x2 - x1)  
      yi = y1 + factAB * (y2 - y1)  
c
      dist = sqrt((xp-xi)**2+(yp-yi)**2)
c
      return
      end
c
c
c
      subroutine read_flx(extflx_file_opt,nds, flx,ierr)
      use mod_params
      implicit none
c     include 'params'
c
c reads the flux onto the target plates provided in an external file
c
c Geier IPP/02
c      
      integer nds, extflx_file_opt,ierr
      real flx(maxe2dizs,maxnds)
c
      integer cell,id
c
c
c     Read type 0 flux file: 
c
      if (extflx_file_opt.eq.0) then 
c
c        open(90, file='/u/geier/divimp/shots/ext_flux.dat',status='old')
c
         open(tmpunit, file='ext_flux.dat',
     >        status='old',err=100,iostat=ierr)
c
c
         do id=1, nds 
            read (tmpunit,'( i4, 21(e12.4E2))',
     >            err=110,end=110,iostat=ierr)   cell, flx(:,id)
c           write (0,1000) cell, flx(:,id)
c           read (tmpunit,*) flx(id)
c           write (0,*) flx(:,id)
         end do
         close(tmpunit)
c
      endif
c
c      	  
      return
c
c     Deal with error conditions
c
c     Open - likely file doesn't exist 
c
100   write(6,'(a,i5)') 
     >      'ERROR: Problem opening ext_flux.dat file: ERRNO=',ierr
      write(0,'(a,i5)') 
     >      'ERROR: Problem opening ext_flux.dat file: ERRNO=',ierr

      return 
c
c     Read errors
c
110   write(6,'(a,i5)') 
     >      'ERROR: Problem opening ext_flux.dat file: ERRNO=',ierr
      write(0,'(a,i5)') 
     >      'ERROR: Problem opening ext_flux.dat file: ERRNO=',ierr

      return 
      end
c
c
c
      subroutine calc_extflx_yield(fydata,matt)
      use mod_params
      use mod_cgeom
      use mod_cneut2
      implicit none
c     include 'params'
c     include 'cgeom' 
c     include 'cneut2'
      real fydata(maxpts,5)
      integer matt       
c
c     CALC_EXTFLX_YIELD:
c
c     This routine calculates the NET yield and F*Y for 
c     a set of external flux data read in from a file.
c
c     The external flux data can be modified to specify the 
c     bombarding species and could conceivably be used for 
c     simultaneous bombardment by several different species and 
c     charge states.
c
c
c     Load the flux data from a file - if an error occurs - issue
c     error message and load zero fluxes.
c 
c     Local variables 
c
      integer iopt,id,ch_st
      integer extflx_file_opt
      integer matext,ierr
      real extflx(maxe2dizs,maxnds)
      real fy_tmp(5)
c
      real yield
      external yield
c     
c
c     This variable can be moved to the input file as 
c     an optional input file to support a more general format
c     for the data file containing the external fluxes.  
c
      extflx_file_opt = 0
c
c     Read external fluxes
c
      ierr = 0
c
      call read_flx(extflx_file_opt,nds,extflx,ierr) 
c
c     Exit if there was a problem reading the external fluxes
c
      if (ierr.ne.0) return
c
      if (extflx_file_opt.eq.0) then 
c
c
c Geier IPP/02 include option for physical sputtering with fluxes from 
c external file
c
c        This option assumes sputtering by Argon and an input file with
c        a specific format.  
c
         matext = 9   

c
         do id = 1,nds
c 
c           physical sputtering for all fluids of interest in ext. file
c	   
c           first compute total flux of all ionization stages
c           used for normalization later:
c
            do ch_st=1,18
c
c              +3 because the 1st three fluids in the file are D+, He+ and He++      
c
	       fydata(id,1)=fydata(id,1)+extflx(ch_st+3, id)
c
            end do
c	       
c           File format assumes Argon
c
            do ch_st=1,18
c
               fy_tmp(1)=extflx(ch_st+3,id) 
               fy_tmp(2)=2.0*KTIDS(ID)+ 3.0*REAL(ch_st)*KTEDS(ID)
               fy_tmp(3)=fy_tmp(2) + 2.0*REAL(ch_st)*KTEDS(ID) 
c
c              Calculate yield by external plasma material
c
               fy_tmp(4) = YIELD (MATEXT,MATT,fy_tmp(2),
     >                         kteds(id),ktids(id)) * KMFPS(ID)
c
c              energies are normalized to the respective relative fluxes:
c
               fydata(id,2)=fydata(id,2)+fy_tmp(2)*
     >                     fy_tmp(1)/fydata(id,1)
               fydata(id,3)=fydata(id,3)+fy_tmp(3)*
     >                     fy_tmp(1)/fydata(id,1)
c
c              Sum up totals 
c
               fydata(id,4)=fydata(id,4)+fy_tmp(4)
               fydata(id,5)=fydata(id,5)+fy_tmp(1)*fy_tmp(4)
c 
	    end do
c
	    fydata(id,4)=fydata(id,5)/fydata(id,1)
c
          end do 
c
      endif 
c
      return 
      end 
c
c
c      	 
      real function ndrand(sigma, x0) result(ran_n) 
c
c calculates normally distributed random variable with
c mean x0 and standard deviation sigma
c based on numerical recipes
c
c Alexander Geier IPP/01, last revisited 4.12.2001
c
      implicit none
c
      integer :: nrand
      real :: ran1, ran2, ran_sqr, fac
      real :: sigma, x0
      real :: lb, ub
      real, dimension(2) :: ran
c
      real*8 :: seed
c 
      nrand=0
c
c      lb=4.656612873077392578e-10
c      ub=1.
c
      do
c
c         call r_addrans(ran,2,lb,ub)
c
         CALL SURAND (SEED, 2, RAN)
         nrand = nrand + 1
c
c         call random(ran,2,lb,ub)
c
         ran(1)=2.*ran(1)-1.
         ran(2)=2.*ran(2)-1.
	 ran_sqr=ran(1)**2+ran(2)**2
c
         if (ran_sqr.lt.1.and.ran_sqr.ne.0.and.
     >       nrand.lt.10)  exit
c
      end do
c
c     Check error condition
c
      if (nrand.ge.10) then 
c
         write(6,'(a)') 'ERROR: Problem in ndrand'//
     >              ' : Too many attempts to find random numbers'
         write(0,'(a)') 'ERROR: Problem in ndrand'//
     >              ' : Too many attempts to find random numbers'
c
         ran_n = x0
c
      else
c
         fac=sqrt(-2. *log(ran_sqr)/ran_sqr)
         ran_n=sigma*ran(1)*fac+x0  
c
      endif
c
      end function ndrand










