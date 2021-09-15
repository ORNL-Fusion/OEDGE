module plasma_overlay
  use error_handling

  private

      real,allocatable :: rmap(:,:),zmap(:,:),dens(:,:),temp(:,:)
      real,allocatable :: raxis(:),zaxis(:)
  
      integer :: nr,nz


      real :: rmin,rmax,zmin,zmax

      public load_plasma_overlay,interpolate_overlay,get_overlay_limits



contains

  subroutine load_plasma_overlay(file,ierr)
    implicit none
    integer :: ierr
    character*(*) :: file

    integer :: iunit
    integer ir,iz,in,irn,izn,irt,izt,id
    !character*2048 :: line


    ! open and load the plasma file - the first line should contain the dimensions

    call find_free_unit_number(iunit)

    write(0,'(a)') 'EXTERNAL PLASMA OVERLAY: OPENING FILE = '//trim(file)

    open(iunit,file=trim(file),status='OLD',iostat=ierr)


    if (ierr.ne.0) then 
       call errmsg('ERROR:OVERLAY_PLASMA: PROBLEM OPENING FILE:',ierr)
       return
    endif

    read(iunit,*) nr,nz
    
    call allocate_storage(nr,nz,ierr)
    if (ierr.ne.0) return


    ! DTS data is expected on an evenly spaced regular mesh
    ! Both Raxis and zaxis should be in ascending order - these
    ! will be used to identify the cell for interpolation

    do ir = 1,nr
       do iz = 1,nz
          read(iunit,*) irn,izn,rmap(ir,iz),zmap(ir,iz),dens(ir,iz),temp(ir,iz)
          if (irn.ne.ir.or.izn.ne.iz) then 
             write(0,*) 'OVERLAY PLASMA: INDEX MISMATCH READING FILE:',irn,ir,izn,iz
          endif
          if (ir.eq.1) then 
             zaxis(iz) = zmap(ir,iz)
          endif
          !write(6,'(a,2i8,4(1x,g18.8))') 'OP:',ir,iz,rmap(ir,iz),zmap(ir,iz),dens(ir,iz),temp(ir,iz)

       end do
       raxis(ir) = rmap(ir,1)

    end do

    rmin = minval(raxis)
    rmax = maxval(raxis)

    zmin = minval(zaxis)
    zmax = maxval(zaxis)

    close(iunit)

  end subroutine load_plasma_overlay

  subroutine allocate_storage(nr,nz,ierr)
    implicit none
    integer :: nr,nz,ierr


    if (allocated(rmap)) deallocate(rmap)
    allocate(rmap(nr,nz),stat=ierr)

    if (allocated(zmap)) deallocate(zmap)
    allocate(zmap(nr,nz),stat=ierr)

    if (allocated(dens)) deallocate(dens)
    allocate(dens(nr,nz),stat=ierr)

    if (allocated(temp)) deallocate(temp)
    allocate(temp(nr,nz),stat=ierr)

    if (allocated(raxis)) deallocate(raxis)
    allocate(raxis(nr),stat=ierr)

    if (allocated(zaxis)) deallocate(zaxis)
    allocate(zaxis(nz),stat=ierr)

    if (ierr.ne.0) then 
       call errmsg('ERROR:OVERLAY_PLASMA: PROBLEM ALLOCATING STORAGE:',ierr)
       return
    endif


  end subroutine allocate_storage

  subroutine interpolate_overlay(r,z,te,ne)
    use mod_interpolate
    implicit none

    real :: r,z,te,ne
    real*8 :: rcell(4),zcell(4),dens_val(4),temp_val(4)
    real*8 :: interp_val
    integer :: irn,izn,irt,izt,in
    logical :: check
    integer,external :: ipos

    ! if point is outside box - loop
    if (r.lt.rmin.or.r.gt.rmax.or.z.lt.zmin.or.z.gt.zmax) then
       return
    endif

    ! find cell - gives nearest higher value

    irn = ipos(r,raxis,nr)
    izn = ipos(z,zaxis,nz)

    ! if either entry is one then it is an error
    if (irn.eq.1.or.izn.eq.1) then 
       return
    endif

    ! The desired cell is the four points
    ! (irn,izn-1),(irn-1,izn-1),(irn-1,izn),(irn,izn)

    ! check the values for the four vertices to make
    ! sure none are zero (no fitted value)

    check = .false.
    do irt = irn-1,irn
       do izt = izn-1,izn
          if (dens(irt,izt).le.0.0.or.temp(irt,izt).le.0.0) then 
             check = .true.
          endif
       end do
    end do

    if (check) then 
       return
    endif

    ! now need to interpolate ... convert to double precision to use
    ! existing interpolation routines

    rcell(1) = rmap(irn,izn-1)
    rcell(2) = rmap(irn-1,izn-1)
    rcell(3) = rmap(irn-1,izn)
    rcell(4) = rmap(irn,izn)

    zcell(1) = zmap(irn,izn-1)
    zcell(2) = zmap(irn-1,izn-1)
    zcell(3) = zmap(irn-1,izn)
    zcell(4) = zmap(irn,izn)

    dens_val(1) = dens(irn,izn-1)
    dens_val(2) = dens(irn-1,izn-1)
    dens_val(3) = dens(irn-1,izn)
    dens_val(4) = dens(irn,izn)

    temp_val(1) = temp(irn,izn-1)
    temp_val(2) = temp(irn-1,izn-1)
    temp_val(3) = temp(irn-1,izn)
    temp_val(4) = temp(irn,izn)

    !write(6,'(a,2(1x,i8))')    'index:',irn,izn
    !write(6,'(a,4(1x,g18.8))') 'rcell:',(rcell(in),in=1,4)
    !write(6,'(a,4(1x,g18.8))') 'zcell:',(zcell(in),in=1,4)
    !write(6,'(a,4(1x,g18.8))') 'dens :',(dens_val(in),in=1,4)
    !write(6,'(a,4(1x,g18.8))') 'temp :',(temp_val(in),in=1,4)


    call cell_interpolate(dble(r),dble(z),interp_val,rcell,zcell,dens_val)
    ! overwrite density
    ne = interp_val



    call cell_interpolate(dble(r),dble(z),interp_val,rcell,zcell,temp_val)
    ! overwrite temperature
    te = interp_val

  end subroutine interpolate_overlay

  subroutine get_overlay_limits(ormin,ormax,ozmin,ozmax)
    implicit none
    real :: ormin,ormax,ozmin,ozmax
    ormin = rmin
    ormax = rmax
    ozmin = zmin
    ozmax = zmax
    return
  end subroutine get_overlay_limits



end module plasma_overlay
