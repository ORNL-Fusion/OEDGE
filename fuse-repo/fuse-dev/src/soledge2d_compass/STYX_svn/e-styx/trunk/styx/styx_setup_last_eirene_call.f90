subroutine styx_setup_last_eirene_call
  use eirmod_precision
  use eirmod_parmmod
  use eirmod_cplot
  use eirmod_ctrcei
  use eirmod_ctext
  use eirmod_comusr
  use eirmod_csdvi
  use eirmod_cspei
  use eirmod_comprt, only : iunout
  use eirmod_cpes
  use styx2eirene
  implicit none
  integer :: ispec,i
  real(dp) :: CH2max
  character(3) :: splot
  

  if (my_pe == 0) then


! switch on plotting routines
  NLPL2D=.true.
  NLPL3D=.false.
  PLHST=.true.
  TRCHST=.true.
  PLCUT(3)=.true.
  PL1ST=.true.

! printout of EIRENE global balance

  TRCBLA=.true.
  TRCBLM=.true.
  TRCBLI=.true.
  TRCBLP=.true.
  TRCBLE=.true.

! set centre and size of the box

  CH2X0 = 0.5_dp*(Rmaxg+Rming)
  CH2Y0 = 0.5_dp*(Zmaxg+Zming)
  CH2Z0 = 0._dp
  CH2MX = 1.2_dp*0.5_dp*(Rmaxg-Rming)
  CH2MY = 1.2_dp*0.5_dp*(Zmaxg+abs(Zming))

! set the same size in X and Y to avoid distorsions
  CH2max=max(CH2MX,CH2MY)
  CH2MX=CH2max
  CH2MY=CH2max

! plot the particles for the stratum defined in the input file

  I1TRC = Npart_Eirene(1)*(stratum_plot-1)
  I2TRC = I1TRC + nhist_plot

! title for the geometry plot

  write(splot,'(i3)') stratum_plot

  if (stratum_plot > 1) then
      ! turn off flux weighted allocation of CPU for stratum
      ! otherwise number of particles per stratum not know before run
      ! avoid in production runs, may deteriorate statistics on last iteration
      	ALLOC=0._dp
  endif
     
  TXTRUN = fluid_code//'-Eirene'

  if (stratum_plot <= Nrecyc) then				
	   TXTRUN = trim(adjustl(TXTRUN))//': recycling source, stratum 1'
  elseif (stratum_plot > Nrecyc .and. stratum_plot <(Nrecyc+Nrecomb)) then
      	    TXTRUN = trim(adjustl(TXTRUN))//': recombination source, stratum 2'
  elseif(stratum_plot > Nrecyc+Nrecomb+1) then
      	   if (.not.timedep) then
      	     TXTRUN = trim(adjustl(TXTRUN))//': gas puff, stratum '//trim(adjustl(splot))
      			
       	   elseif (timedep .and. stratum_plot == NSTRATA) then
       	     TXTRUN = trim(adjustl(TXTRUN))//': source from previous time step, stratum '//trim(adjustl(splot))
       	   else
      	     TXTRUN = trim(adjustl(TXTRUN))//': gas puff, stratum '//trim(adjustl(splot))
       	   endif
  endif

  write(iunout,*) 'geometry plotting turned on, iteration # ',curiter

  call eirene_plt2D

  endif

!!!!!!!!!! set up data for noise estimation !!!!!!!!!!
! all the arrays needs to be redimentionned (done in initialization phase)

! species (add a loop if multispecies)
 ispec =1

! these as would be provided by find_param
 NCV = 9
 NSD = 28
 NSDW = 59

 
! otherwise reallocation not carried out
! deallocate(SDVI1)
! needs to be reallocated
! deallocate(SDVI2,SIGMAC,SGMCS,ISDVI)
! deallocate(IIHC,IGHC)

 
 call eirene_dealloc_csdvi
 call eirene_dealloc_cspei
 call eirene_alloc_csdvi(1)
 call eirene_alloc_cspei
 call eirene_alloc_csdvi(2)

! be careful these are pointers 
! volume tallies
 NSIGVI=28
! surface tallies
 NSIGSI=59
! covariances
 NSIGCI=9

! needed in mcarlo, otherwise =0 and nothing is done ...
 NSIGI=NSIGVI+NSIGSI+NSIGCI

 if (my_pe == 0) then

 allocate(italv(NSIGVI)) 
 italv = (/ 1, 2, 3, 4, 5, 6, 7, 9, 15, 21, 33, 38, 39, 44, 45, &
          50, 85, 86, 87, 89, 90, 91, 93, 94, 95, 97, 98, 99 /)
  
  do i=1,NSIGVI
  	IGH(i)=ispec
        IIH(i)=italv(i)
  enddo

! surface tallies
  
  do i=1,NSIGSI
  	IGHW(i)=ispec
        IIHW(i)=i
  enddo

! correlations between tallies (to estimate variance of the sum)
  
  
  IGHC(1,1)=1
  IIHC(1,1)=97 
  IGHC(2,1)=1 
  IIHC(2,1)=98
  
  IGHC(1,2)=1
  IIHC(1,2)=97
  IGHC(2,2)=1
  IIHC(2,2)=99

  IGHC(1,3)=1
  IIHC(1,3)=98
  IGHC(2,3)=1
  IIHC(2,3)=99

  IGHC(1,4)=1 
  IIHC(1,4)=33
  IGHC(2,4)=1
  IIHC(2,4)=39
 
  IGHC(1,5)=1
  IIHC(1,5)=33
  IGHC(2,5)=1
  IIHC(2,5)=45

  IGHC(1,6)=1
  IIHC(1,6)=39
  IGHC(2,6)=1
  IIHC(2,6)=45

  IGHC(1,7)=1
  IIHC(1,7)=38
  IGHC(2,7)=1
  IIHC(2,7)=44

  IGHC(1,8)=1
  IIHC(1,8)=38
  IGHC(2,8)=1
  IIHC(2,8)=50
 
  IGHC(1,9)=1
  IIHC(1,9)=44
  IGHC(2,9)=1
  IIHC(2,9)=50

  endif
 

end subroutine
