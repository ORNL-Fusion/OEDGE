module plasma_overlay
  use error_handling

  private

      real,allocatable :: rmap(:,:),zmap(:,:),data(:,:,:)
      real,allocatable :: raxis(:),zaxis(:)
  
      integer :: nr,nz


      real :: rmin,rmax,zmin,zmax

      public load_plasma_overlay,interpolate_overlay,interpolate_extdata,get_overlay_limits,load_extdata,po_deallocate_storage



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

    read(iunit,*,iostat=ierr) nr,nz
    if (ierr.ne.0) then 
       call errmsg('ERROR:OVERLAY_PLASMA: PROBLEM READING FILE:',file)
       return
    else
       ! Routine loads the data from the file - first 2 specifies two datasets, second is to read index numbers for checking
       call load_extdata(iunit,nr,nz,2,2,ierr)
    endif
    
    close(iunit)

  end subroutine load_plasma_overlay

  subroutine load_extdata(iunit,nr,nz,opt,form_opt,ierr)
    implicit none
    integer :: iunit,nr,nz,ierr
    integer :: opt
    integer :: form_opt  ! option to control format of data read - 1 = index numbers included, 2 = no index numbers 
    integer :: id
    integer :: irn,izn,ir,iz
    
    ierr = 0
    call po_allocate_storage(nr,nz,ierr,opt)
    if (ierr.ne.0) return

    ! Data is expected on an evenly spaced regular mesh
    ! Both Raxis and zaxis should be in ascending order - these
    ! will be used to identify the cell for interpolation
    ! Opt specifies the number of data points on a line
    ! 1 or 2 are currently supported

    do ir = 1,nr
       do iz = 1,nz
          if (form_opt.eq.1) then 
             ! no index numbers
             read(iunit,*,iostat=ierr) rmap(ir,iz),zmap(ir,iz),(data(ir,iz,id),id=1,opt)
             if (ierr.ne.0) then
                call errmsg('ERROR: PLASMA_OVERLAY.F90 : ERROR READING DATA LINE WITHOUT INDEXING: IR,IZ = ',[ir,iz])
                return
             endif
          elseif (form_opt.eq.2) then
             ! index numbers included
             read(iunit,*,iostat=ierr) irn,izn,rmap(ir,iz),zmap(ir,iz),(data(ir,iz,id),id=1,opt)
             if (ierr.ne.0) then
                call errmsg('ERROR: PLASMA_OVERLAY.F90 : ERROR READING DATA LINE WITH INDEXING: IR,IZ = ',[ir,iz])
                return
             endif
          endif

          if (form_opt.eq.2.and.(irn.ne.ir.or.izn.ne.iz)) then 
             call errmsg('ERROR: PLASMA_OVERLAY.F90 : LOAD_EXTDATA : INDEX NUMBER MISMATCH READING FILE')
             ierr = 1
             return
          endif

          if (ir.eq.1) then 
             zaxis(iz) = zmap(ir,iz)
          endif
          !write(6,'(a,2i8,4(1x,g18.8))') 'OP:',ir,iz,rmap(ir,iz),zmap(ir,iz),(data(ir,iz,id),in=1,opt)

       end do
       raxis(ir) = rmap(ir,1)

    end do

    rmin = minval(raxis)
    rmax = maxval(raxis)

    zmin = minval(zaxis)
    zmax = maxval(zaxis)

    !write(6,*) 'LOAD_EXTDATA:',rmin,rmax,zmin,zmax,opt
    !do ir = 1,nr
    !   do iz = 1,nz
    !      write(6,'(2i8,2(1x,f10.5),10(1x,g12.5))')  ir,iz,raxis(ir),zaxis(iz),(data(ir,iz,id),id=1,opt)
    !   end do
    !end do
    
  end subroutine load_extdata


  
  subroutine po_allocate_storage(nri,nzi,ierr,opt)
    use allocate_arrays
    implicit none
    integer :: nri,nzi,ierr
    integer :: opt

    ! set module values of nr and nz since these are used in the interpolation routines.
    nr = nri
    nz = nzi
    
    ! opt is used to specify the number of data sets on a line - 1 or 2 currently supported. 

    call allocate_array(rmap,nr,nz,'rmap',ierr)
    call allocate_array(zmap,nr,nz,'rmap',ierr)
    call allocate_array(data,nr,nz,opt,'data',ierr)

    call allocate_array(raxis,nr,'raxis data',ierr)
    call allocate_array(zaxis,nz,'zaxis data',ierr)

    if (ierr.ne.0) then 
       call errmsg('ERROR: PLASMA_OVERLAY.F90: PROBLEM ALLOCATING STORAGE:',ierr)
       return
    endif


  end subroutine po_allocate_storage

  subroutine interpolate_overlay(r,z,te,ne)
    implicit none

    real :: r,z,te,ne

    ! interpolate ne - data set 1
    call interpolate_val(r,z,ne,1,1)
    
    ! interpolate te - data set 2
    call interpolate_val(r,z,te,2,1)

  end subroutine interpolate_overlay


  subroutine interpolate_extdata(r,z,val,opt,check_opt)
    implicit none
    integer :: opt
    real   :: r,z,val
    integer,optional :: check_opt
    integer :: check
    
    if (.not.present(check_opt)) then
       check = 0
    else
       check = check_opt
    endif
    
    ! interpolate val - data set 1
    call interpolate_val(r,z,val,opt,check)
    
  end subroutine interpolate_extdata

  

  subroutine interpolate_val(r,z,val,opt,check_opt)
    use mod_interpolate
    implicit none
    real :: r,z,val
    integer :: check_opt   ! if the values at any of the corners are zero the routine exits without assigning a value
    integer :: opt         ! specifies which data set to interpolate

    real*8 :: rcell(4),zcell(4),data_val(4)
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

    if (check_opt.eq.1) then
       check = .false.
       do irt = irn-1,irn
          do izt = izn-1,izn
             if (data(irt,izt,opt).le.0.0) then 
                check = .true.
             endif
          end do
       end do

       if (check) then 
          return
       endif
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

    data_val(1) = data(irn,izn-1,opt)
    data_val(2) = data(irn-1,izn-1,opt)
    data_val(3) = data(irn-1,izn,opt)
    data_val(4) = data(irn,izn,opt)

    !write(6,'(a,2(1x,i8))')    'index:',irn,izn
    !write(6,'(a,4(1x,g18.8))') 'rcell:',(rcell(in),in=1,4)
    !write(6,'(a,4(1x,g18.8))') 'zcell:',(zcell(in),in=1,4)
    !write(6,'(a,4(1x,g18.8))') 'dens :',(dens_val(in),in=1,4)
    !write(6,'(a,4(1x,g18.8))') 'temp :',(temp_val(in),in=1,4)

    call cell_interpolate(dble(r),dble(z),interp_val,rcell,zcell,data_val)

    ! assign interpolated value
    val = interp_val

    !write(6,'(a,2i8,2(1x,f10.5),10(1x,g12.5))') 'IV:',irn,izn,r,z,val

    
  end subroutine interpolate_val


  subroutine get_overlay_limits(ormin,ormax,ozmin,ozmax)
    implicit none
    real :: ormin,ormax,ozmin,ozmax
    ormin = rmin
    ormax = rmax
    ozmin = zmin
    ozmax = zmax
    return
  end subroutine get_overlay_limits

  subroutine po_deallocate_storage
    implicit none
    if (allocated(rmap)) deallocate(rmap)
    if (allocated(zmap)) deallocate(zmap)
    if (allocated(data)) deallocate(data)
    if (allocated(raxis)) deallocate(raxis)
    if (allocated(zaxis)) deallocate(zaxis)

  end subroutine po_deallocate_storage

  


end module plasma_overlay
