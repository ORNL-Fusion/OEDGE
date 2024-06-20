subroutine styx_init_surface_data(ical)
  use all_variables, only : interp_data2, global_parameters, reference_parameters
  use Mphysics
  use eirmod_precision
  use styx2eirene
  use eirmod_comprt, only : iunout
! these variables are initialized only for the second call (ical=2)
  use eirmod_parmmod, only : NSTS,NLIM
  use eirmod_comusr, only : BXIN,BYIN,BZIN,BFIN
  use eirmod_ctrig, only : PTRIX,PTRIY
  use eirmod_ccona, only : eps60
  implicit none

  integer, intent(in) :: ical

  integer :: ITRI,IK,ip1,ip2
  integer :: dummy,J,ISIDE
  integer :: IPRO,IS,isurf,irecsurf,ic
  integer :: icount,isou,iv,i
  integer :: Ncore
  integer :: isurff(2)

  real(dp) :: pi_
  real(dp) :: xG,yG,R1,Z1,R2,Z2
  real(dp) :: sumw
  integer  :: nextvert,nextvert0,nver

  integer, allocatable :: ixtri(:),iytri(:)
  integer :: ier
  real*8 :: dum

  real*8, allocatable :: triRG(:),triZG(:)

  real*8 :: Bn,BsX,BsY,BsZ,Bs,sina
  real*8 :: uperpX,uperpY,uperpZ
  real*8 :: n_X,n_Y

  real*8 :: tst,tend,omp_get_wtime
  integer :: ithr,omp_get_thread_num

  character(15) :: format,format2
  character(9) :: dums
  character(3) :: Nksis,Ntaus
  integer :: il,iu,im

  real*8, parameter :: epsilon0 = 8.85418782d-12

  integer*4 :: k
  real*8,allocatable :: dOmega_tab(:)
  logical :: isWall

  pi_=4._dp*atan(1._dp)

  if (ical == 1) then
    	
  ! now read matlab preprocessing output


  	open(unit=11,file=trim(adjustl(fluid_code))//'.npco_char')

  	read(11,*) NRKNOT_styx

  ! now get knots coordinates
  	allocate(xtriang(NRKNOT_styx),ytriang(NRKNOT_styx))
  	allocate(ixtri(NTRI_styx),iytri(NTRI_styx))

  	if (levgeo_styx == 1) then
        ! cylindrical geometry
  	! set LZ, unit=cm (periodicity still to be implemented ?) !
  		LZ=2._dp*pi_*reference_parameters%geometry%R0*100._dp
  	else
  	! toroidal geometry, LZ not needed
  		LZ=0._dp
  	endif

  ! main loop on cells
  	allocate(NVERT(3,NTRI_styx),NEIGH(3,NTRI_styx))
  	allocate(NSIDE(3,NTRI_styx),IPROP(3,NTRI_styx))

  	do IK=1,NRKNOT_styx
     		read(11,*) dummy,xtriang(IK),ytriang(IK)
  	enddo

  	close(11)


  	open(unit=22,file=trim(adjustl(fluid_code))//'.elemente')

  	read(22,*) NTRI_styx

  	do ITRI=1,NTRI_styx
     		read(22,*) dummy, NVERT(1,ITRI), NVERT(2,ITRI), NVERT(3,ITRI)
  	enddo

  	close(22)


  	open(unit=33,file=trim(adjustl(fluid_code))//'.neighbors')

  	read(33,*)

  	do ITRI=1,NTRI_styx
     	read(33,*) dummy, NEIGH(1,ITRI), NSIDE(1,ITRI), IPROP(1,ITRI), NEIGH(2,ITRI), NSIDE(2,ITRI), IPROP(2,ITRI), NEIGH(3,ITRI), NSIDE(3,ITRI), IPROP(3,ITRI), IXTRI(ITRI), IYTRI(ITRI)
  	enddo

  	close(33)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calculate weights for interpolations (distances)
  ! wv2c: vertices to triangular cell centers (styx2D->EIRENE) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  	allocate(wv2c(NTRI_styx,3))

  	wv2c=0._dp

  	do ITRI=1,NTRI_styx
     	! first center of mass coordinates
     		xG=1._dp/3._dp*(xtriang(NVERT(1,ITRI))+xtriang(NVERT(2,ITRI))+xtriang(NVERT(3,ITRI)))
     		yG=1._dp/3._dp*(ytriang(NVERT(1,ITRI))+ytriang(NVERT(2,ITRI))+ytriang(NVERT(3,ITRI)))

     		wv2c(ITRI,1:3)=1._dp/sqrt((xtriang(NVERT(1:3,ITRI))-xG)**2+(ytriang(NVERT(1:3,ITRI))-yG)**2)
     		sumw=sum(wv2c(ITRI,:),1)
     		wv2c(ITRI,1:3)=wv2c(ITRI,1:3)/sumw

  	enddo

  	! number of triangle sides which are NON ABSORBING material surfaces
  	! the frontier with the core should be absorbing
  	Nsou=0
  	do ITRI=1,NTRI_styx
  		do ISIDE=1,3
  			if (IPROP(ISIDE,ITRI) /= 0 .and. IPROP(ISIDE,ITRI) /=1) Nsou=Nsou+1
  		enddo
  	enddo

  	! number of absorbing surface on core edge boundary

  	Ncore=0
  	do ITRI=1,NTRI_styx
  		do ISIDE=1,3
  			if (IPROP(ISIDE,ITRI) == 1) Ncore=Ncore+1
  		enddo
  	enddo

	! total number of surfaces where fluxes have to be calculated
  	NSURF_TAL = Nsou+Ncore

  	! array that will contain cell volumes

    !*****************
    ! Modif by Patrick
   	! allocate(vol_tri_eirene(NTRI_styx))
   	allocate(vol_tri_eirene(Neir_cells))
   	!*****************

  ! now deallocate the stuff no longer needed (keep NVERT for checks in interpolation routine)

  	deallocate(ixtri,iytri)

  elseif (ical == 2) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! there are currently three surface labeling systems : 
!  1) EIRENE standard surface numbering (surface(isurf)%...)
!  2) Recycling surface numbering, i.e. absorbing surfaces (on core-edge interface) excluded (recsurf(irecsurf)%...)
!     This is to speed up exchanges of fluxes (just one loop), in eirene_get_fluxes
!       mappings are provided in the following arrays:
!       isurf    = kr2e(irecsurf)
!       irecsurf = ke2r(isurf)
!  3) Ordering following the wall, in trigo. direction starting from the hfs
!       ksurf provides the mapping between EIRENE numbering and ordered numbering
!        isurf    = ksurf(iord)
!        irecsurf = krecsurf(iord) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! getting the right surface numbers for flux data exchanges
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    NSURF_TAL=NSURF_TAL+NSTS+NLIM
    nsurf0=NSTS+NLIM+1

    isurf=0

    allocate(recsurf(1:Nsou))

    do ITRI=1,NTRI_styx
        do ISIDE=1,3
           if (IPROP(ISIDE,ITRI) /= 0 .and. IPROP(ISIDE,ITRI) /= 1) then
  		isurf=isurf+1

  		if (isurf > Nsou) then
   			write(*,*) ' isurf > Nsou in surface initialization, stop'
   			call eirene_exit_own(1)
  		endif

  		if (ISIDE == 1) then
                  IP1=NVERT(1,ITRI)
                  IP2=NVERT(2,ITRI)
                  R1=xtriang(NVERT(1,ITRI))
                  Z1=ytriang(NVERT(1,ITRI))
                  R2=xtriang(NVERT(2,ITRI))
                  Z2=ytriang(NVERT(2,ITRI))
                elseif (ISIDE == 2) then
                  IP1=NVERT(2,ITRI)
                  IP2=NVERT(3,ITRI)
                  R1=xtriang(NVERT(2,ITRI))
  		  Z1=ytriang(NVERT(2,ITRI))
  		  R2=xtriang(NVERT(3,ITRI))
  		  Z2=ytriang(NVERT(3,ITRI))
                else
                  IP1=NVERT(3,ITRI)
                  IP2=NVERT(1,ITRI)
                  R1=xtriang(NVERT(3,ITRI))
  		  Z1=ytriang(NVERT(3,ITRI))
  		  R2=xtriang(NVERT(1,ITRI))
  		  Z2=ytriang(NVERT(1,ITRI))
                endif

  		recsurf(isurf)%itri  = ITRI
   		recsurf(isurf)%iside = ISIDE
  		recsurf(isurf)%v1    = ip1
  		recsurf(isurf)%v2    = ip2
                recsurf(isurf)%R1    = R1
   		recsurf(isurf)%Z1    = Z1
  		recsurf(isurf)%R2    = R2
  		recsurf(isurf)%Z2    = Z2

                recsurf(isurf)%ds=sqrt((recsurf(isurf)%R1-recsurf(isurf)%R2)**2 + &
  						(recsurf(isurf)%Z1-recsurf(isurf)%Z2)**2)*1e-2_dp
   	   endif
  	enddo
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! inverse mapping irecsurf = f(itri,iside) ... wastefull in terms of storage !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(recsurfinv(3,Ntri_styx))

    do itri=1,Ntri_styx
      do iside=1,3
        do isurf=1,Nsou
          if (recsurf(isurf)%itri == itri .and. recsurf(isurf)%iside == iside) then
            recsurfinv(iside,itri)=isurf
          endif
        enddo
      enddo
    enddo
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! mapping between surface numbers in EIRENE and triangles number and sides
! as in read_triangles.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IPROPMIN=MINVAL(IPROP(1:3,1:NTRI_styx),MASK=(IPROP(1:3,1:NTRI_styx)/=0))
    IPROPMAX=MAXVAL(IPROP(1:3,1:NTRI_styx),MASK=(IPROP(1:3,1:NTRI_styx)/=0))

    allocate(surface(NSURF_TAL))
    allocate(Rmin(IPROPMIN:IPROPMAX),Rmax(IPROPMIN:IPROPMAX))
    allocate(Zmin(IPROPMIN:IPROPMAX),Zmax(IPROPMIN:IPROPMAX))
    allocate(smin(IPROPMIN:IPROPMAX),smax(IPROPMIN:IPROPMAX))

    isurf=NSTS+NLIM

    Rmin=1.e8_dp
    Rmax=0._dp 
    Zmin=1.e8_dp
    Zmax=0._dp

    do IPRO=IPROPMIN,IPROPMAX
      do IS=1,3
        do ITRI=1,NTRI_styx
          if (IPRO == IPROP(IS,ITRI)) then
            isurf=isurf+1
  	    surface(isurf)%itri=ITRI
  	    surface(isurf)%iside=IS
  	    surface(isurf)%iprop=IPRO
            ! map sides and vertex numbers
            if (IS == 1) then
              surface(isurf)%R1=xtriang(NVERT(1,ITRI))
              surface(isurf)%Z1=ytriang(NVERT(1,ITRI))
              surface(isurf)%R2=xtriang(NVERT(2,ITRI))
  	      surface(isurf)%Z2=ytriang(NVERT(2,ITRI))
  	      surface(isurf)%v1=NVERT(1,ITRI)
  	      surface(isurf)%v2=NVERT(2,ITRI)
  	    elseif (IS == 2) then
  	      surface(isurf)%R1=xtriang(NVERT(2,ITRI))
  	      surface(isurf)%Z1=ytriang(NVERT(2,ITRI))
  	      surface(isurf)%R2=xtriang(NVERT(3,ITRI))
  	      surface(isurf)%Z2=ytriang(NVERT(3,ITRI))
	      surface(isurf)%v1=NVERT(2,ITRI)
  	      surface(isurf)%v2=NVERT(3,ITRI)
  	    elseif (IS == 3) then
  	      surface(isurf)%R1=xtriang(NVERT(3,ITRI))
  	      surface(isurf)%Z1=ytriang(NVERT(3,ITRI))
  	      surface(isurf)%R2=xtriang(NVERT(1,ITRI))
  	      surface(isurf)%Z2=ytriang(NVERT(1,ITRI))
  	      surface(isurf)%v1=NVERT(3,ITRI)
  	      surface(isurf)%v2=NVERT(1,ITRI)
  	    endif
           ! length of the surface element in poloidal plane (m)
            surface(isurf)%ds=sqrt((surface(isurf)%R1-surface(isurf)%R2)**2 + &
  		(surface(isurf)%Z1-surface(isurf)%Z2)**2)*1e-2_dp

           ! now find origin for curvilinear coordinates along wall/cei : surface farthest on the low field side
! also for geometry plot, to set box size 
            if (surface(isurf)%R1 < Rmin(IPRO)) then
              Rmin(IPRO)=surface(isurf)%R1
              smin(IPRO)=isurf
            endif
            if (surface(isurf)%R1 > Rmax(IPRO)) then
              Rmax(IPRO)=surface(isurf)%R1
  	    endif
            if (surface(isurf)%Z1 < Zmin(IPRO)) then
              Zmin(IPRO)=surface(isurf)%Z1
  	    endif
            if (surface(isurf)%Z1 > Zmax(IPRO)) then
              Zmax(IPRO)=surface(isurf)%Z1
            endif
          endif
        enddo
      enddo
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now mapping between surface numbers in EIRENE and actual relative position of
! surface; ksurf(i) = row of the surface i, starting from the point closest from
! the low field side ; ssurf(i) = curvilinear distance from that origin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Rming=1.e8_dp
    Rmaxg=0._dp
    Zming=1.e8_dp
    Zmaxg=0._dp

    ! choose origin (-> isurf for which x, i.e. R, is smallest)
    do IPRO=IPROPMIN,IPROPMAX
      if (IPRO /= 1) then
        if (Rmin(IPRO) < Rming) then
          Rming=Rmin(IPRO)
          isurf=smin(IPRO)
        endif
        if (Rmax(IPRO) > Rmaxg) then
          Rmaxg=Rmax(IPRO)
        endif
        if (Zmin(IPRO) < Zming) then
  	  Zming=Zmin(IPRO)
  	endif
  	if (Zmax(IPRO) > Zmaxg) then
  	  Zmaxg=Zmax(IPRO)
  	endif
      endif
    enddo

    deallocate(Rmin,Rmax,Zmin,Zmax,smin,smax)
    allocate(ksurf(nsurf0:NSURF_TAL+1),ssurf(nsurf0:NSURF_TAL+1))

    ksurf=0
    ssurf=0._dp

    ksurf(nsurf0)=isurf
    ssurf(nsurf0)=surface(ksurf(nsurf0))%ds

    if (surface(isurf)%iprop == 1) then
      write(*,*) ' problem with cei/wall position ...'
      call eirene_exit_own(1)
    endif

    ! set positive direction along the wall
    if (surface(isurf)%Z1 > surface(isurf)%Z2) then
      nextvert=surface(isurf)%v2
    elseif (surface(isurf)%Z1 < surface(isurf)%Z2) then
      nextvert=surface(isurf)%v1
    else
    ! strict equality (horinzontal) unlikely ... do it in x
      if (surface(isurf)%R1 > surface(isurf)%R2) then
        nextvert=surface(isurf)%v2
      else
        nextvert=surface(isurf)%v1
      endif	
    endif

    nextvert0=nextvert

    open(unit=555,file=trim(adjustl(fluid_code))//'.wall_segments',status='replace')
   
    write(555,'(E14.7,1x,E14.7,1x,i3)') xtriang(nextvert),ytriang(nextvert),surface(isurf)%iprop

    ! follow the wall
    do i=2,NSou+1
      j=NSTS+NLIM+i
  			
      ! identify wall surfaces that contain vertex nextver
      ! each vertex belong to one or more triangles, and each of these triangles has up to two sides on the wall
      ! here there are only two hits possible because we check relation between vertex number and side
      icount=0
      do isou=nsurf0,NSURF_TAL
        if (surface(isou)%iprop /= 1) then
          ITRI=surface(isou)%itri
          ISIDE=surface(isou)%iside
          do iv=1,3
            if (NVERT(iv,itri) == nextvert) then
              if (ISIDE == 1) then
                if (iv == 1 .or. iv == 2) then
                  icount=icount+1
                  isurff(icount)=isou
  		endif
              elseif (ISIDE == 2) then
                if (iv == 2 .or. iv == 3) then
  		  icount=icount+1
                  isurff(icount)=isou
  	        endif
              elseif (ISIDE == 3) then
                if (iv == 1 .or. iv == 3) then
                  icount=icount+1
                  isurff(icount)=isou
                endif
              endif				
            endif
          enddo
        endif
      enddo
  			
      if (icount == 1) then
        ksurf(j)=isurff(1)
      elseif (icount == 2 ) then
      ! distinguish between previous and current cell, in case triangle sides where identical
        if (isurff(1) == isurf) then
          ksurf(j)=isurff(2)
        elseif (isurff(2) == isurf) then
          ksurf(j)=isurff(1)
        else
          write(*,*) 'problem, initial surface not found !'
          call eirene_exit_own(1)
        endif
      endif

      ! check that it is not absorbing otherwise problem ...
      if (surface(ksurf(j))%iprop == 1) then
        write(*,*) 'problem, next surface on CEI ... that should not be !'
  	call eirene_exit_own(1)
      endif			

      ! determine next vertex (one of the two making the new surface is the old one...)
      nver=surface(ksurf(j))%v1
      if (nver == nextvert) then
        nextvert=surface(ksurf(j))%v2
      else
        if (nextvert == surface(ksurf(j))%v2) then
  	  nextvert=surface(ksurf(j))%v1
        else
          write(*,*) 'problem with nextver, isurf = ',ksurf(j)
        endif
      endif
  				
      isurf=ksurf(j)
      ssurf(j)=ssurf(j-1)+surface(isurf)%ds
      ! nextver and isurf are now set

      write(555,'(E14.7,1x,E14.7,1x,i3)') xtriang(nextvert),ytriang(nextvert),surface(ksurf(j))%iprop

    enddo
  	
    ! check that the loop is closed !
    if (nextvert /= nextvert0) then
      write(*,*) ' problem, the wall surface does not form a closed contour !'
      write(*,*) ' nextverO = ',nextvert0
      write(*,*) ' nextver = ',nextvert
      call eirene_exit_own(1)
    endif
    close(555)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   Now mapping between full eirene numbering and recycling surface numbering            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(kr2e(Nsou),ke2r(Nsurf_tal))

kr2e=0
ke2r=0

ic=0

do isurf=nsurf0,Nsurf_tal
  do irecsurf=1,Nsou
    if (recsurf(irecsurf)%itri == surface(isurf)%itri .and. recsurf(irecsurf)%iside == surface(isurf)%iside) then
      kr2e(irecsurf)=isurf
      ke2r(isurf)=irecsurf
      ic=ic+1
    endif
  enddo
enddo

if (ic /= Nsou) then
  write(*,*) ' warning, mapping between surface elements incomplete … '
  write(*,*) ' This may lead to diagnostic errors '
endif

!!! direct correspondance between recycling surfaces and ordering on the wall

allocate(krecsurf(Nsou))

! ksurf(i+nsurf0) = isurf (index for surface(isurf)% ...
! krecsur(i) = irecsur

do i=1,Nsou
  krecsurf(i)=ke2r(ksurf(i+nsurf0))
enddo


    if (.not.is_3D) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!              Calculation of the radiative flux matrix on the wall                       !!
!!    Delta_Omega(i,itri) = angle under which wall element i is seen from triangle itri    !!
!!                   sum(Delta_omega,1) = 2pi for each triangle                            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(Delta_Omega(Nsurf_tal,Ntri_styx))

! test if the radiation matrix is there, otherwise recalculate it

    open(unit=666,file='radiation_matrix',status='old',iostat=ier)
    close(666)

    if (ier.ne.0) then
      rad_mat=.true.
      write(*,*) 'Warning: radiation matrix not found, recalculating it ...'
    endif

    if (rad_mat) then

       Delta_Omega=0.d0

       allocate(triRG(Ntri_styx),triZG(Ntri_styx))

      ! initialize
      Delta_Omega = 0.d0
      triRG=(xtriang(nvert(1,:))+xtriang(nvert(2,:))+xtriang(nvert(3,:)))/3.d0
      triZG=(ytriang(nvert(1,:))+ytriang(nvert(2,:))+ytriang(nvert(3,:)))/3.d0

      tst = omp_get_wtime()

      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ithr,isurf,dOmega_tab,k,isWall)

      allocate(dOmega_tab(Ntri_styx))
      !$OMP DO
      do isurf=nsurf0,nsurf_tal
         ithr = omp_get_thread_num()
         !         ithr=1
         !$OMP CRITICAL
         write(*,'(A9,i3,A25,i6)') 'thread # ',ithr,' processing wall element #',isurf
         !$OMP END CRITICAL
         call styx_radiation_weights(isurf,Ntri_styx,triRG,triZG,dOmega_tab,isWall)
         !$OMP CRITICAL
         if(isWall) then
            do k=1,Ntri_styx
               Delta_Omega(isurf,k)=dOmega_tab(k)
            end do
         end if
         !$OMP END CRITICAL
      enddo
      !$OMP END DO
       
      deallocate(dOmega_tab) 
      !$OMP END PARALLEL

      tend = omp_get_wtime()

      write(*,*) ' time for radiation weights calculation (s) = ',tend-tst

      do itri=1,Ntri_styx
        do i=nsurf0,nsurf_tal
          if (isnan(Delta_Omega(i,itri))) then
            write(*,*) 'nan catched, itri = ',itri,' ,isou = ',i
          endif
        enddo
      enddo

      open(unit=666,file='radiation_matrix',form='unformatted',status='replace')
      do itri=1,Ntri_styx
        write(666) (Delta_Omega(i,itri),i=1,Nsurf_tal)
      enddo
      close(666)

!     Omega=sum(Delta_omega,1)

!     open(unit=666,file='Omega',status='replace')
!     do itri=1,Ntri_styx
!       write(666,*) Omega(itri)
!     enddo
!     close(666)


      !!!!! yannick for nf paper

      open(unit=666,file='rad_weights_200',status='replace')
      do itri=1,Ntri_styx
        write(666,*) Delta_Omega(200,itri),Delta_Omega(200,itri),Delta_Omega(200,itri)
      enddo
      close(666)

      open(unit=666,file='rad_weights_300',status='replace')
      do itri=1,Ntri_styx
        write(666,*) Delta_Omega(300,itri),Delta_Omega(300,itri),Delta_Omega(300,itri)
      enddo
      close(666)

    else !rad_mat

  !!!!!!!!!!!! load radiative matrix !!!!!!!!!!!

      write(*,*) 'loading radiation matrix ...'

      open(unit=666,file='radiation_matrix',form='unformatted',status='old')
      read(666,iostat=ier) (Delta_Omega(i,1),i=1,Nsurf_tal)
      if (ier == 0) then
        do itri=2,Ntri_styx
          read(666) (Delta_Omega(i,itri),i=1,Nsurf_tal)
        enddo
        close(666)
   ! but need to test if there are still data (restarting stationnary from a time dep. run
      read(666,iostat=ier) dum
      if (ier == 0) then
        write(*,*) ' radiation matrix too big, exit ... '
        write(*,*) ' If restarting in stationary mode from a time dependent mode, recalculate matrix'
        call eirene_exit_own(1)
      endif
         

      elseif (ier < 0 .and. timedep) then
        write(*,*) '------------------------------------------------------------'
        write(*,*) ' loading of radiation matrix failed, end of file ...        '
        write(*,*) ' retry assuming previous run was not time dependent ...     '
        
        rewind(666)
        do itri=1,Ntri_styx
          read(666) (Delta_Omega(i,itri),i=1,Nsurf_tal-2)
        enddo
         
      endif
      write(*,*) 'ok ...'
    
 
    endif

    deallocate(NSIDE)

   endif !  not is_3D

   allocate(Radflx(nsurf_tal))


  endif !(ical)

end subroutine styx_init_surface_data
