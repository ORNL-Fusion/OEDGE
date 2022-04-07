module mod_sol22_lim

  use mod_plasma_data

  implicit none

  private

    real*8,allocatable :: ne_local(:),te_local(:),ti_local(:),vb_local(:),spts_local(:)
    
    integer,parameter :: maxlen = 1000
    character*(maxlen) :: sol22_parameter_filename

    public :: sol22,load_sol22_parameter_file, set_sol22_parameter_file


contains




  
  subroutine sol22
    use mod_params
    !use mod_comt2
    !use mod_comtor
    !use allocate_arrays
    use mod_calcsol_interface
    !use yreflection
    use mod_solparams
    use mod_allocate_sol22_storage
    use mod_sol22_input_lim
    use mod_sol22_input
    use mod_sol22_utils
    !use mod_comt2
    !use mod_comxyt
    use mod_plasma_data

    integer :: in
    !character*100 :: sol22_fname

    integer :: ierr
    integer :: ncnt

    real*8 :: nb0,te0,ti0,ringlen,r8cizb,r8crmb

    ! The way SOL22 is setup it re-uses the storage for the two half-ring calculations.
    ! It also uses "+" ve target velocity which then needs to be intelligently assigned
    ! depending on the situation. 
    !real :: cs

    !integer :: pz

    !integer :: new_unit
    !integer :: iqx, ixout
    !integer,external :: ipos
    !real :: xbnd1,xbnd2
    !real :: pbnd1,pbnd2

    ! mxspts is assigned to the same value as maxn in set_plasma_data_axis - it is now a variable not a parameter
    mxspts = maxn

    call allocate_sol22_storage

    call allocate_locals(mxspts)

    ! the solver parameters can't be loaded until after the solver storage has been assigned. This arrangement is going to cause some churn
    ! with dynamically allocated memory as sol22 gets run for each row of the plasma background - but it is only done once for the run - have
    ! to assess if it is a significant performance impact.
    
    call load_sol22_parameter_file(ierr)

    !call allocate_array(nb,nys/2,'nb',ierr)
    !call allocate_array(te,nys/2,'te',ierr)
    !call allocate_array(ti,nys/2,'ti',ierr)
    !call allocate_array(vb,nys/2,'vb',ierr)
    !call allocate_array(spts,nys/2,'spts',ierr)

    !call allocate_array(nb,mxspts,'nb',ierr)
    !call allocate_array(te,mxspts,'te',ierr)
    !call allocate_array(ti,mxspts,'ti',ierr)
    !call allocate_array(vb,mxspts,'vb',ierr)
    !call allocate_array(spts,mxspts,'spts',ierr)

    !write(0,*) 'SOL22A: start',mxspts

    !if (nys/2.gt.mxspts) then
    !   call errmsg('MOD_SOL22_LIM:','MXSPTS IS SMALLER THAN NYS/2 - INCREASE MXSPTS IN MOD_SOLCOMMON.F90')
    !   write(0,*) 'MXSPTS=',mxspts,' NYS/2=',nys/2
    !   write(0,*) 'LIM exiting'
    !   stop 'MOD_SOL22_LIM: MXSPTS TOO SMALL'
    !endif

    !
    ! Find IX first outboard index

    !IXOUT = IPOS (-1.E-10, XS, NXS-1)

    ! background plasma and geometry information
    r8cizb = cizb_local
    r8crmb = crmb_local
    ringlen = ring_length

    ! initialize pzone to 1 for now
    !pz = 1 

    ! moved to tau.f - in the plasma_overlay routine
    !do ir = 1,nsol22_opt

    ! Initialize unstructured inputs
    !call sol22_initialize_unstructured_input

    ! pick a unit number - reassign stdin to the new unit number
    ! open the file with the SOL22 option specifications

    !call find_free_unit_number(new_unit)
    !call set_unit_numbers(in_stdin=new_unit)
    !open(file=trim(sol22_filenames(ir)),unit=new_unit,form='formatted',status='old')

    ! read input options for solver 
    !call readsol(ierr)

    !close(new_unit)
    ! reset the unit numbers
    !call reset_unit_numbers

    !xbnd1 = sol22_regions(ir,1)
    !xbnd2 = sol22_regions(ir,2)

    ! ignore pbnds for now
    ! create two zones to overlay SOL22
    ! pzone = 1 -> |P| =< CPC0 (width of probe)
    ! pzone = 2 -> |P| > CPC0
    ! CPC0 
    ! if maxpzones = 1 OR not 3D (MAXNPS =1) then only
    ! calculate for pzone = 1
    ! 
    ! if there are absorbing surfaces - use them
    !    - calculate plasma over [yabsorb1a,yabsorb2a]
    !    - NOTE: there will need to be duplication of
    !      all plasma arrays for 1:pzone areas
    !      e.g. ctembs(ix,iz,pzone)
    !           ctegs ... and forces by iz and pzone
    !
    ! if not ignore
    !    - calculate plasma over whole range and duplicate
    !    - plasma from limiter surface to limiter surface
    !    - SOL22 is NOT applied to pzone 2 since there are
    !      no surfaces there.
    !    - This means this is only applicable for non-3D cases 
    !          
    !
    !pbnd1 = sol22_regions(ir,3)
    !pbnd2 = sol22_regions(ir,4)
    ! SOL22 input specifies which poloidal zone to use the model 
    !pz = sol22_regions(ir,3)

    !do ix = 1,nxs

    ! calculate y locations where the plasma needs to be calculated and
    ! map them to 0 ... ymax (needs to be mapped back afterward)
    ! Also need starting nb, Te, Ti values extracted from existing plasma
    ! or input. Take limits from SOLEDGE code. 

    ! ctembs(ix,iy) ... iy = -nys,nys 
    ! Two sides 
    ! Y < 0 = nys ... nys/2+1 OR -1 ... -nys

    ! only replace specified regions with SOL22 plasma 

    ! take qedges into account for Y coordinates
    ! qedges(iqx,1) for Y< 0 limiter surface
    ! qedges(iqx,2) for Y> 0 limiter surface
    ! HOWEVER - values are +ve for both so edge for Y<0 is -qedges(iqx,1)
    !

    !if (xouts(ix).ge.xbnd1.and.xouts(ix).le.xbnd2) then 

    !   iqx = iqxs(ix)
    !   ncnt = 0

    !write(0,*) 'SOL22:',ix,iqx,ixout,xouts(ix),xbnd1,xbnd2


    ! only overlay SOL22 outboard when no absorbing surfaces are present
    !  if (yabsorb_opt.eq.0.and.ix.le.ixout) then 

    !     ncnt = 0
    !     do iy = nys,nys/2+1,-1
    !        ! values with spts < 0 will be assigned the value of spts = 0
    !        if ((2.0*cl - qedges(iqx,1) - youts(iy)) .gt. 0.0) then 
    !           ncnt = ncnt +1
    !           spts(ncnt) = 2.0 * cl - qedges(iqx,1) - youts(iy)
    !        endif
    !     end do

    !     te0 = qtembs(iqx,1)  ! Y < 0 
    !     ti0 = qtembsi(iqx,1) ! Y < 0 
    !     nb0 = qrnbs(iqx,1)   ! Y < 0 


    te0 = te_bnd_lower
    ti0 = ti_bnd_lower
    nb0 = n_bnd_lower
    !cs  = -getcs_sol22(sngl(te0),sngl(ti0))  ! negative flow toward first surface

    ncnt = mxspts/2  ! do first half of the ring 
    spts_local = spts

    ! calculate SOL

    call calcsol_interface(nb0,te0,ti0,ringlen,ncnt,spts,ne,te,ti,vb,r8crmb,r8cizb)

    ! assign values back to arrays - this is only the first half of the solution
    ! the second half will be copied cell by cell with a different spts_local assigned

    ne = ne_local
    te = te_local
    ti = ti_local
    vb = -vb_local ! negative flow towards first surface (lower bound surface)


    !do iy = maxn,maxn/2+1,-1
    !   if ((2.0*cl - qedges(iqx,1) - youts(iy)) .gt. 0.0) then 
    !      ncnt = ncnt +1
    !      crnbs(ix,iy,pz) = nb(ncnt)
    !      ctembs(ix,iy,pz) = te(ncnt)
    !      ctembsi(ix,iy,pz) = ti(ncnt)
    !      velplasma(ix,iy,pz) = -vb(ncnt)
    !   else
    !      crnbs(ix,iy,pz) = nb0
    !      ctembs(ix,iy,pz) = te0
    !      ctembsi(ix,iy,pz) = ti0
    !      velplasma(ix,iy,pz) = cs
    !   endif
    !end do


    ! Y > 0 is 1...nys/2 or -nys...-nys/2-1

    !ncnt = 0
    !do iy = 1,nys/2
    !   if ((youts(iy)-qedges(iqx,2)) .gt. 0.0) then 
    !      ncnt = ncnt + 1
    !      spts(ncnt) = youts(iy) - qedges(iqx,2)
    !   endif
    !end do


    !te0 = qtembs(iqx,2)  ! Y > 0 
    !ti0 = qtembsi(iqx,2) ! Y > 0 
    !nb0 = qrnbs(iqx,2)   ! Y > 0 

    te0 = te_bnd_upper
    ti0 = ti_bnd_upper
    nb0 = n_bnd_upper
    !cs  = getcs_sol22(sngl(te0),sngl(ti0))  ! positive flow towards second surface
    ! map revised spts values
    do in = mxspts,mxspts/2+1,-1
       spts_local(mxspts-in+1) = ring_length-spts(in)
    end do
    ncnt = mxspts-mxspts/2

    ! calculate SOL

    call calcsol_interface(nb0,te0,ti0,ringlen,ncnt,spts_local,ne_local,te_local,ti_local,vb_local,r8crmb,r8cizb)

    ! copy second call to SOL22 solver into the upper half of the requested plasma solution
    do in = mxspts,mxspts/2+1,-1
       ne(in) = ne_local(mxspts-in+1)
       te(in) = ti_local(mxspts-in+1)
       ti(in) = ti_local(mxspts-in+1)
       vb(in) = vb_local(mxspts-in+1)
    end do



    !ncnt = 0
    !do iy = 1,nys/2
    !   if ((youts(iy)-qedges(iqx,2)) .gt. 0.0) then 
    !      ncnt = ncnt +1
    !      
    !      crnbs(ix,iy,pz) = nb(ncnt)
    !      ctembs(ix,iy,pz) = te(ncnt)
    !      ctembsi(ix,iy,pz) = ti(ncnt)
    !      velplasma(ix,iy,pz) = vb(ncnt)
    !   else
    !      crnbs(ix,iy,pz) = nb0
    !      ctembs(ix,iy,pz) = te0
    !      ctembsi(ix,iy,pz) = ti0
    !      velplasma(ix,iy,pz) = cs
    !   endif

    !end do


    ! copy over modified ctembs etc to -nys:-1

    !do iy = 1,nys
    !   crnbs(ix,-nys+iy-1,pz) = crnbs(ix,iy,pz)
    !   ctembs(ix,-nys+iy-1,pz) = ctembs(ix,iy,pz)
    !   ctembsi(ix,-nys+iy-1,pz) =  ctembsi(ix,iy,pz)
    !   velplasma(ix,-nys+iy-1,pz) =  velplasma(ix,iy,pz)
    !end do

    ! copy velplasma from pz = 1 to pz = 2
    !do in = 2,maxpzone
    !   velplasma(:,:,in) = velplasma(:,:,1)
    !end do


    call deallocate_sol22_storage

    call deallocate_locals

    call calculate_tgrad_e
    
  end subroutine sol22

  subroutine set_sol22_parameter_file(filename)
    implicit none
    character*(*) filename
    integer :: tlen

    tlen=min(len(trim(filename)),maxlen)

    sol22_parameter_filename = filename(1:tlen)
    
  end subroutine set_sol22_parameter_file
  
  subroutine load_sol22_parameter_file(ierr)
  !subroutine load_sol22_parameter_file(filename,ierr)
    use mod_io_units
    use mod_sol22_input
    implicit none
    !character*(*) :: filename
    integer :: new_unit,ierr
    !     Initialize unstructured inputs
    call sol22_initialize_unstructured_input

    !     pick a unit number - reassign stdin to the new unit number
    !     open the file with the SOL22 option specifications

    call find_free_unit_number(new_unit)
    call set_unit_numbers(in_stdin=new_unit)
    open(file=trim(sol22_parameter_filename),unit=new_unit,form='formatted',status='old')

    !     read input options for solver 
    call readsol(ierr)

    close(new_unit)
    !     reset the unit numbers
    call reset_unit_numbers

  end subroutine load_sol22_parameter_file


  
  subroutine allocate_locals(n)
    use allocate_arrays
    implicit none
    integer :: n,ierr
      call allocate_array(spts_local,n,'Local sol22 s',ierr)
      call allocate_array(ne_local,n,'Local sol22 ne',ierr)
      call allocate_array(te_local,n,'Local sol22 te',ierr)
      call allocate_array(ti_local,n,'Local sol22 ti',ierr)
      call allocate_array(vb_local,n,'Local sol22 vb',ierr)
      return
  end subroutine allocate_locals

  subroutine deallocate_locals
    implicit none

    if (allocated(spts_local)) deallocate(spts_local)
    if (allocated(ne_local)) deallocate(ne_local)
    if (allocated(te_local)) deallocate(te_local)
    if (allocated(ti_local)) deallocate(ti_local)
    if (allocated(vb_local)) deallocate(vb_local)
    return

  end subroutine deallocate_locals

end module mod_sol22_lim
