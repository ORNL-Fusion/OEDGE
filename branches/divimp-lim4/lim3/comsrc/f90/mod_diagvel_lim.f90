module mod_diagvel
  use mod_params
  use mod_diagvel_unstruc
  use debug_options
  implicit none

  !
  !     this file contains the declarations for debugging the
  !     velocity distribution. when this option is not needed
  !     these variables can be minimized so that memory is not
  !     used.
  !
  !
  logical,public :: debugv = .false.
  !integer,public :: debug_v_opt = 0
  !integer,public :: ti_calc_opt = 0


  !data vtime /0.1, 0.2, 0.4, 0.6, 0.8, 1.0,1.5, 2.0, 3.0, 5.0/

  !
  ! jdemod - moved from mod_div6.f90 - only use these for LIM - not the rest
  !
  real*8,public,allocatable :: ddvs(:,:,:)
  real*8,public,allocatable :: ddvs2(:,:,:)

  real,public,allocatable :: sdtimp(:,:,:)
  real,public,allocatable :: vtig_array(:,:,:)



  public :: allocate_mod_diagvel,deallocate_mod_diagvel, init_diagvel, update_diagvel, finish_diagvel, print_diagvel,calc_vtig_array

contains


  subroutine init_diagvel
    implicit none

    ! return if velocity debugging is not on
    if (.not.debugv) return

    ddvs = 0.0
    ddvs2 = 0.0
    sdtimp = 0.0

  end subroutine init_diagvel



  subroutine update_diagvel(ix,iy,iz,sputy,fvel)
    implicit none
    integer :: ix,iy,iz
    real*8 :: fvel,sputy  
    ! note sputy includes timestep factor
    ! return if velocity debugging is not on
    if (.not.debugv) return

    ddvs(ix,iy,iz) = ddvs(ix,iy,iz) + sputy * fvel
    ddvs2(ix,iy,iz) = ddvs2(ix,iy,iz) + sputy * fvel**2.0

    !write(6,*) 'update_diagvel:',ix,iy,iz,sputy,fvel,ddvs(ix,iy,iz)

  end subroutine update_diagvel



  subroutine finish_diagvel(qtim,nizs)
    use mod_comxyt
    use mod_dynam1
    use mod_comtor
    implicit none
    integer :: nizs
    integer :: ix,iy,iz,iqx
    real :: qtim
    real :: vav1,vav2
    ! return if velocity debugging is not on
    if (.not.debugv) return

    ! combine -Y and +Y results - make them match
    ! do symmetric contributions for ddvs here - not inline in LIM
    ! normalization is ok here
    do ix = 1,nxs
       do iy = -nys,-1
          do iz = 1,nizs
             iqx = iqxs(ix)
             if (ddlims(ix,iy,iz).gt.0.0) then
                vav1 =  ddvs(ix,iy,iz) /ddlims(ix,iy,iz)/(qtim*qs(iqx))
             endif
             
             if (ddlims(ix,iy+nys+1,iz).gt.0.0) then
                vav2 =  ddvs(ix,iy+nys+1,iz) /ddlims(ix,iy+nys+1,iz)/(qtim*qs(iqx))
             endif
                
             !if (iz.eq.nizs) &
             !     write(6,'(a,4i8,20(1x,g12.5))') 'ddvs1:',ix,iy,iy+nys+1,iz,xouts(ix),youts(iy),ddvs(ix,iy,iz),ddvs(ix,iy+nys+1,iz),&
             !          ddlims(ix,iy,iz),ddlims(ix,iy+nys+1,iz),vav1,vav2
             ddvs(ix,iy+nys+1,iz) = ddvs(ix,iy+nys+1,iz) + ddvs(ix,iy,iz)
             ddvs(ix,iy,iz) =   ddvs(ix,iy+nys+1,iz)
             ddvs2(ix,iy+nys+1,iz) = ddvs2(ix,iy+nys+1,iz) + ddvs2(ix,iy,iz)
             ddvs2(ix,iy,iz) =   ddvs2(ix,iy+nys+1,iz)
          end do
       end do
    end do

    ! normalize and calculate temperature
    do ix = 1,nxs
       do iy = -nys,nys
          do iz = 1,nizs
             ! calculate mean ion velocity
             if (ddlims(ix,iy,iz).gt.0.0) then 
                iqx = iqxs(ix)
                !if (iz.eq.nizs) &
                !   write(6,'(a,4i8,20(1x,g12.5))') 'ddvs2:',ix,iy,iz,iqx,qs(iqx),xouts(ix),youts(iy),ddvs(ix,iy,iz),ddlims(ix,iy,iz),ddvs(ix,iy,iz)/ddlims(ix,iy,iz)/(qtim*qs(iqx))
                ddvs(ix,iy,iz) =   ddvs(ix,iy,iz) /ddlims(ix,iy,iz)/(qtim*qs(iqx))
                ddvs2(ix,iy,iz) =   ddvs2(ix,iy,iz) /ddlims(ix,iy,iz)/(qtim*qs(iqx))**2
                if (ti_calc_opt.eq.0) then 
                   sdtimp(ix,iy,iz)=(ddvs2(ix,iy,iz)-ddvs(ix,iy,iz)**2) / 9.639e7 * crmi
                elseif (ti_calc_opt.eq.1) then
                   sdtimp(ix,iy,iz)=(ddvs2(ix,iy,iz)-ddvs(ix,iy,iz)**2) / 2.4544e8 * crmi
                endif
             else
                iqx = iqxs(ix)
                !if (iz.eq.nizs) &
                !   write(6,'(a,4i8,20(1x,g12.5))') 'ddvs2:',ix,iy,iz,iqx,qs(iqx),xouts(ix),youts(iy),ddvs(ix,iy,iz),ddlims(ix,iy,iz),0.0
             endif
          end do
       end do
    end do


  end subroutine finish_diagvel

  subroutine print_diagvel(qtim,nizs)
    use mod_io_units
    use mod_vtig
    use mod_comtor
    use mod_comxyt
    use mod_comt2
    use mod_dynam1
    implicit none
    integer :: nizs
    real :: qtim
    
    real :: vb,vtig
    integer :: ix,iy,iz,pz,iqx
    integer :: outunit
    real :: tgscal,scal
    ! return if velocity debugging is not on
    if (.not.debugv) return

    ! assign poloidal zone 1 for now
    !pz =1 

    call find_free_unit_number(outunit)

    TGSCAL = (1.6E-19)/(CRMI*1.673E-27) * QTIM *QTIM 

    ! first output file
    open(outunit,file='velocity_diagnostic_data.out',form='formatted')

    do pz = 1,maxpzone

       write(outunit,'(100(a,1x))') 'DIAGVEL:','IX','IY','PZ','X','Y','vb','vtig','vtotal','IZ','vz','DENSITY','FLUX','dTi/ds','Ti'
    do ix = 1,nxs
       do iy = -nys,nys
          iqx = iqxs(ix)
          scal = tgscal * qs(iqx) * qs(iqx)
          ! print out

          vtig = calc_vtig(ix,iy,1,qtim)

          if (crnbs(ix,iy).gt.0.0) then 
             vtig = integration_const/crnbs(ix,iy,pz) * ctembsi(ix,iy,pz)**(1.5) * ctigs(ix,iy,pz) / scal
          else
             vtig = 0.0
          endif

          vb = velplasma(ix,iy,pz)
          !write(outunit,'(a,2i8,100(1x,g12.5))') 'VEL:',ix,iy,xouts(ix),youts(iy),vb,vtig,&
          !     (real(iz),ddvs(ix,iy,iz),ddlims(ix,iy,iz),ddlims(ix,iy,iz)*ddvs(ix,iy,iz),iz=1,nizs),ctigs(ix,iy),ctembsi(ix,iy),integration_const,scal,ctigs(ix,iy)/scal
          write(outunit,'(a,2i8,100(1x,g12.5))') 'VEL:',ix,iy,pz,xouts(ix),youts(iy),vb,vtig,vb+vtig,&
               (real(iz),ddvs(ix,iy,iz),ddlims(ix,iy,iz),ddlims(ix,iy,iz)*ddvs(ix,iy,iz),iz=nizs,nizs),ctigs(ix,iy,pz)/scal,ctembsi(ix,iy,pz)
       end do
    end do

    end do
    close(outunit)
    
    ! second output file - moved to OUT
    !open(outunit,file='density_velocity_flux.out',form='formatted')

    ! write out background velocities vb and vtig

    !  if (debugv) then
    !     IXOUT = IPOS (-1.E-10, XS, NXS-1)
    !     ix = (nxs+ixout)/2
    !     write(outunit,'(a)') 'VELOCITIES:'
    !     write(outunit,'(a,1000(1x,g12.5))') ' YOUTS: ',(youts(iy),iy=1,nys)
    !     write(outunit,'(a,1000(1x,g12.5))') ' VB: ',(velplasma(ix,iy,1),iy=1,nys)
    !     write(outunit,'(a,1000(1x,g12.5))') ' VTIG: ',(vtig_array(ix,iy,1),iy=1,nys)
    !     write(outunit,'(a,1000(1x,g12.5))') ' VTOT: ',(velplasma(ix,iy,1)+vtig_array(ix,iy,1),iy=1,nys)
    !  endif

    ! only outputing max charge state for now
    !do iz = nizs,nizs
    !     write(outunit,'(a)') ' '
    !     write(outunit,'(a,i8)') ' DENSITY IZ=',iz
    !     write(outunit,'(26x,1000(1x,g12.5))') (ywids(iy),iy=1,nys)
    !     write(outunit,'(26x,1000(1x,g12.5))') (youts(iy),iy=1,nys)
    !     do ix = 1,nxs
    !        write(outunit,'(1000(1x,g12.5))') xwids(ix),xouts(ix),&
    !                   (ddlims(ix,iy,iz),iy=1,nys)
    !     end do
    !     
    !     write(outunit,'(a)') ' '
    !     write(outunit,'(a,i8)') ' VELOCITY IZ=',iz
    !     write(outunit,'(26x,1000(1x,g12.5))') (ywids(iy),iy=1,nys)
    !     write(outunit,'(26x,1000(1x,g12.5))') (youts(iy),iy=1,nys)
    !     do ix = 1,nxs
    !        write(outunit,'(1000(1x,g12.5))') xwids(ix),xouts(ix),&
    !                   (ddvs(ix,iy,iz),iy=1,nys)
    !     end do
    !     
    !     write(outunit,'(a)') ' '
    !     write(outunit,'(a,i8)') ' FLUX IZ=',iz
    !     write(outunit,'(26x,1000(1x,g12.5))') (ywids(iy),iy=1,nys)
    !     write(outunit,'(26x,1000(1x,g12.5))') (youts(iy),iy=1,nys)
    !     do ix = 1,nxs
    !        write(outunit,'(1000(1x,g12.5))') xwids(ix),xouts(ix),&
    !                   (ddlims(ix,iy,iz)*ddvs(ix,iy,iz),iy=1,nys)
    !     end do
    !                  
    !end do
    !close(outunit)
    
  end subroutine print_diagvel

  real function calc_vtig(ix,iy,ip,qtim)
    use mod_io_units
    use mod_vtig
    use mod_comtor
    use mod_comxyt
    use mod_comt2
    use mod_dynam1
    use mod_lambda
    implicit none
    integer :: ix,iy,ip
    integer :: iqx
    real :: tgscal,scal,qtim
    real :: lambda
    integer :: pz
    pz = pzones(ip)
    
    TGSCAL = (1.6E-19)/(CRMI*1.673E-27) * QTIM *QTIM 

    iqx = iqxs(ix)
    scal = tgscal * qs(iqx) * qs(iqx)
    !
    ! jdemod - allow for spatial variation of lambda
    !
    if (lambda_vary_opt.eq.1) then
       lambda = coulomb_lambda(crnbs(ix,iy,pz),ctembsi(ix,iy,pz))
    else
       lambda = 1.0
    endif
       
    
    ! print out
    if (crnbs(ix,iy).gt.0.0) then 
       calc_vtig = integration_const/crnbs(ix,iy,pz)/lambda * ctembsi(ix,iy,pz)**(1.5) * ctigs(ix,iy,pz) / scal
    else
       calc_vtig = 0.0
    endif

  end function calc_vtig

  subroutine calc_vtig_array(qtim)
    use mod_comxyt
    use allocate_arrays
    implicit none
    real :: qtim
    integer :: ix,iy,ip,ierr

    ! allocate and calculate the vtig velocity array
    ! use temporary storage
    ierr = 0
    call allocate_array(vtig_array,1,maxnxs,-maxnys,maxnys,1,maxpzone,'vtig_array',ierr)
    
    ip = 1

    do ix = 1,nxs
       do iy = -maxnys,maxnys
          vtig_array(ix,iy,ip) = calc_vtig(ix,iy,ip,qtim)
       end do
    end do

  end subroutine calc_vtig_array


  subroutine allocate_mod_diagvel
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_diagvel','ALLOCATE')

    if (debug_v_opt.eq.1) then
       debugv=.true.
    endif

    write(0,*) 'debugv set:', debug_v_opt, debugv

    
    if (debugv) then 

       ! allocation includes initialization to 0.0
       call allocate_array(ddvs,1,maxnxs,-maxnys,maxnys,1,maxizs,'ddvs',ierr)
       call allocate_array(ddvs2,1,maxnxs,-maxnys,maxnys,1,maxizs,'ddvs2',ierr)
       call allocate_array(sdtimp,1,maxnxs,-maxnys,maxnys,1,maxizs,'sdtimp',ierr)

    endif

  end subroutine allocate_mod_diagvel



  subroutine deallocate_mod_diagvel
    implicit none

    call pr_trace('mod_diagvel','DEALLOCATE')

    if (allocated(ddvs)) deallocate(ddvs)
    if (allocated(ddvs2)) deallocate(ddvs2)
    if (allocated(sdtimp)) deallocate(sdtimp)
    if (allocated(vtig_array)) deallocate(vtig_array)
    
  end subroutine deallocate_mod_diagvel

end module mod_diagvel
