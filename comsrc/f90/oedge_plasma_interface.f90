module oedge_plasma_interface
  use error_handling

implicit none

private

integer :: interpolate_opt

integer :: nrs,mks,irsep,nds,npolyp,nvert

integer,allocatable :: nks(:)

!-------------------------------------------------
!
! Geometry information
!
! rs,zs - cell centers
! rbnd,zbnd - r,z, coordinates of cell boundary mid-point along field lines
! 
!
real,allocatable :: rs(:,:),zs(:,:),rbnd(:,:),zbnd(:,:)

! Connection map for grid contains indices of adjacent cells
! ikins - knot index inward
! irins - ring index inward
! ikouts - knot index outward
! irouts - ring index outward
! idds(ir,1..2) - index to target elements for each ring - 1 = end of ring, 2=start of ring
integer,allocatable :: ikins(:,:),irins(:,:),ikouts(:,:),irouts(:,:),idds(:,:)

! Index from grid to polygon data
integer,allocatable :: korpg(:,:)

! polygon data 
! nvertp - vertex count for cell
! rvertp - r coordinates of cell vertices
! zvertp - z coordinates of cell vertices
integer,allocatable :: nvertp(:,:)
real,allocatable :: rvertp(:,:),zvertp(:,:)

! Magnetic field data
! btot - magntitude of magnetic field in each cell 
! br,bz,bt - unit vector components of the magnetic field direction in each cell
real,allocatable :: btot,br,bz,bt


!-------------------------------------------------
!
! Plasma information
!
! Note: plasma file contains some duplication of geometry information
!
! Plasma quantities
! - volume terms
! knbs - plasma density (m-3)
! ktebs - electron temperature (eV)
! ktibs - ion temperature (eV)
! kvhs - parallel velocity (m/s) - negative velocities are towards the low S or start of field line -
!                                - need to check REDEP sign convention
! kes - parallel electric field (V/m)
!
! - target quantities
! knds  - plasma density (m-3)
! kteds - electron temperature (eV)
! ktids - ion temperature (eV)
! kvds  - parallel flow velocity (m/s)
! keds  - electric field (V/m) 
!

real,allocatable :: knbs(:,:),ktebs(:,:),ktibs(:,:),kvhs(:,:),kes(:,:)
real,allocatable :: knds(:),kteds(:),ktids(:),kvds(:),keds(:)





public :: get_oedge_plasma,init_oedge_plasma




contains


subroutine init_oedge_plasma(filename,interpolate_option,ierr)
  implicit none
  integer :: interpolate_option
  character*256 :: gridfilename,plasmafilename
  integer :: igridfile,iplasmafile,ierr,ios

  ! Assign chosen interpolation option
  ! 0 - no interpolation - cell value is returned
  ! 1 - nearest neighbours interpolation
  
  interpolate_opt = interpolate_option

  ! open input files 
  ! read array sizes from the geometry file

  call find_free_unit_number(igridfile)
  call find_free_unit_number(iplasmafile)


  gridfilename = trim(filename)//'.grd'
  plasmafilename = trim(filename)//'.bgp'

!...  Open geometry file:
      OPEN(UNIT=igridfile,FILE=gridfilename,STATUS='OLD',IOSTAT=ios)
      
      if (ios.ne.0) then 

         call errmsg('ERROR OPENING OEDGE GRID FILE:'//trim(gridfilename),ios)

         stop 'ERROR OPENING OEDGE GRID FILE'

      endif

!...  Open plasma file:
      OPEN(UNIT=iplasmafile,FILE=plasmafilename,STATUS='OLD',IOSTAT=ios)
      
      if (ios.ne.0) then 

         call errmsg('ERROR OPENING OEDGE PLASMA FILE:'//trim(plasmafilename),ios)
         stop 'ERROR OPENING OEDGE PLASMA FILE'

      endif

  ! load geometry data and allocate storage

  call load_oedge_grid(igridfile,ierr)


  ! load plasma data

  call load_oedge_plasma(iplasmafile,ierr)


end subroutine init_oedge_plasma


subroutine close_oedge_plasma
  implicit none

  ! deallocate any allocated storage

  call deallocate_arrays


end subroutine close_oedge_plasma


subroutine get_oedge_plasma(r,z,ne,te,ti,vb,ef,btot,br,bz,bt)
  implicit none

  real :: r,z,ne,te,ti,vb,ef

  ! Get the OEDGE plasma conditions at the specified R,Z location
  ! Two options - value in cell and interpolated - value in cell is quicker






end subroutine get_oedge_plasma


subroutine load_oedge_grid(infile,ierr)
implicit none
integer :: infile,ierr

integer :: ik,ir,in
logical :: done

  !
  !     READ Title line
  !
  read (infile,10) buffer
  !
  if (buffer(1:6).ne.'DIVIMP') then
     call errmsg('NOT A DIVIMP GRID FILE',trim(buffer))
     stop 'Not a valid DIVIMP grid file'
  endif

  ! First line contains overal grid size information

  
  done = .false.

  ! initialization
  nrs   = 0
  irsep = 0
  nds   = 0
  npolyp= 0
  nvert = 0


  do while (.not.done)

     read(infile,10,iostat=ios,err=error_exit) buffer

     if (ios.eq.0) then 

        if (buffer(1:4).eq.'NRS:') then
           read (buffer,200) nrs,irsep,nds,npolyp,nvert

           ! Allocate the nks array
           if (allocated(nks)) deallocate(nks)
           if (nrs.gt.0) then 
              allocate(nks(nrs),stat=ierr)
              if (ierr.ne.0) then 
                 call errmsg('ERROR Allocating NKS array:',ierr)
                 stop 'Error allocating NKS array'
              endif
           else
              call errmsg('INVALID GRID FILE: NUMBER OF RINGS: NRS = ',nrs)
              stop 'INVALID GRID FILE'
           endif

        elseif(buffer(1:6).eq.'KNOTS:') then

           if (.not.allocated(nks)) then 
              call errmsg('ERROR NKS Array not yet allocated')
              stop 'NKS not allocated'
           endif

           read (infile,400)  (nks(ir),ir=1,nrs)

           mks = maxval(nks)

           ! Allocate the storage for all of the arrays 
           call allocate_storage

        elseif (buffer(1:3).eq.'RS:') then
           !
           !       R cell center coordinate
           !
           read (infile,500) ((rs(ik,ir),ik=1,nks(ir)),ir=1,nrs)

        elseif (buffer(1:3).eq.'ZS:') then
           !
           !       Z cell center coordinate
           !
           read (infile,500) ((zs(ik,ir),ik=1,nks(ir)),ir=1,nrs)

        elseif (buffer(1:5).eq.'RBND:') then
           !
           !       R cell bound coordinate
           !
           read (infile,500) ((rbnd(ik,ir),ik=0,nks(ir)),ir=1,nrs)

        elseif (buffer(1:5).eq.'ZBND:') then
           !
           !       Z cell bound coordinate
           !
           read (infile,500) ((zbnd(ik,ir),ik=0,nks(ir)),ir=1,nrs)

        elseif (buffer(1:6).eq.'IRINS:') then
           !
           !       IR Inward connection map
           !
           read (infile,400) ((irins(ik,ir),ik=1,nks(ir)),ir=1,nrs)

        elseif (buffer(1:6).eq.'IKINS:') then
           !
           !       IK Inward connection map
           !
           read (infile,400) ((ikins(ik,ir),ik=1,nks(ir)),ir=1,nrs)

        elseif (buffer(1:7).eq.'IROUTS:') then
           !
           !       IR outward connection map
           !
           read (infile,400) ((irouts(ik,ir),ik=1,nks(ir)),ir=1,nrs)

        elseif (buffer(1:7).eq.'IKOUTS:') then
           !
           !       IK outward connection map
           !
           read (infile,400) ((ikouts(ik,ir),ik=1,nks(ir)),ir=1,nrs)

        elseif (buffer(1:5).eq.'IDDS:') then
           !
           !       Target element connection map
           !
           read (infile,400) ((idds(ir,in),ir=1,nrs),in=1,2)

        elseif (buffer(1:6).eq.'KORPG:') then
           !
           !       Cell polygon data index
           !
           read (infile,400) ((korpg(ik,ir),ik=1,nks(ir)),ir=1,nrs)

        elseif (buffer(1:7).eq.'NVERTP:') then
           !
           !       Number of vertices in each polygon
           !
           read (infile,400) (nvertp(in=1,npolyp))

        elseif (buffer(1:7).eq.'RVERTP:') then
           !
           !    R coordinates of polygon vertices
           !
           read (infile,500) (rvertp(ik,in),ik=1,nvert,in=1,npolyp)

        elseif (buffer(1:7).eq.'ZVERTP:') then
           !
           !    Z coordinates of polygon vertices
           !
           read (infile,500) (zvertp(ik,in),ik=1,nvert,in=1,npolyp)

        elseif (buffer(1:5).eq.'BTOT:') then
           !
           !       B-field magnitude in cell
           !
           read (infile,500) ((btot(ik,ir),ik=1,nks(ir)),ir=1,nrs)

        elseif (buffer(1:3).eq.'BR:') then
           !
           !       R part of B-field unit vector
           !
           read (infile,500) ((br(ik,ir),ik=1,nks(ir)),ir=1,nrs)

        elseif (buffer(1:3).eq.'BZ:') then
           !
           !       Z part of B-field unit vector
           !
           read (infile,500) ((bz(ik,ir),ik=1,nks(ir)),ir=1,nrs)

        elseif (buffer(1:3).eq.'BT:') then
           !
           !       T part of B-field unit vector
           !
           read (infile,500) ((bt(ik,ir),ik=1,nks(ir)),ir=1,nrs)

        endif

     else
        done=.true.
     endif
     !
     !     Loop back for continued reading
     !
  end do


  CLOSE (infile)

  return

  error_exit:

  call errmsg('ERROR READING IN DIVIMP GRID:')
  stop 'ERROR reading in grid file'






!
!     Formatting
!
  10  format(a)
 100  format(a40)
 200  format('NRS:',i5,'IRSEP:',i5,'NDS:',i5,'NPOLYP:',i5,'NVERT:',i5)
 400  format(12i6)
 500  format(6e18.10)






end subroutine load_oedge_grid

subroutine allocate_storage
  use allocate_arrays
  implicit none
  integer :: ierr


  ! This routine uses the values set in nrs, mks, nds, npolyp and nvert to allocate the storage for all of the 
  ! grid and plasma array data.

  ! Allocate basic grid geometry arrays
  call allocate_array(rs,mks,nrs,'RS',ierr)
  call allocate_array(zs,mks,nrs,'ZS',ierr)
  call allocate_array(rbnd,0,mks,1,nrs,'RBND',ierr)
  call allocate_array(zbnd,0,mks,1,nrs,'ZBND',ierr)

  ! Allocate connection map arrays
  call allocate_array(irins,mks,nrs,'IRINS',ierr)
  call allocate_array(ikins,mks,nrs,'IKINS',ierr)
  call allocate_array(irouts,mks,nrs,'IROUTS',ierr)
  call allocate_array(ikouts,mks,nrs,'IKOUTS',ierr)
  call allocate_array(idds,nrs,2,'IDDS',ierr)

  ! Allocate polygon map and cell information
  call allocate_array(korpg,mks,nrs,'KORPG',ierr)
  call allocate_array(nvertp,npolyp,'NVERTP',ierr)
  call allocate_array(rvertp,nvert,npolyp,'RVERTP',ierr)
  call allocate_array(zvertp,nvert,npolyp,'ZVERTP',ierr)
  
  ! Allocate magnetic field arrays
  call allocate_array(btot,mks,nrs,'BTOT',ierr)
  call allocate_array(br,mks,nrs,'BR',ierr)
  call allocate_array(bz,mks,nrs,'BZ',ierr)
  call allocate_array(bt,mks,nrs,'BT',ierr)

  ! Allocate plasma arrays
  ! Volume plasma background
  call allocate_array(knbs,mks,nrs,'KNBS',ierr)
  call allocate_array(ktebs,mks,nrs,'KTEBS',ierr)
  call allocate_array(ktibs,mks,nrs,'KTIBS',ierr)
  call allocate_array(kvbs,mks,nrs,'KVBS',ierr)
  call allocate_array(kes,mks,nrs,'KES',ierr)

  ! Target plasma background
  call allocate_array(knds,nds,'KNDS',ierr)
  call allocate_array(kteds,nds,'KTEDS',ierr)
  call allocate_array(ktids,nds,'KTIDS',ierr)
  call allocate_array(kvds,nds,'KVDS',ierr)
  call allocate_array(keds,nds,'KEDS',ierr)


end subroutine allocate_storage

subroutine deallocate_storage
implicit none


! Deallocate storage used in the OEDGE plasma interface

  if (allocated(nks)) deallocate(nks)

  if (allocated(rs)) deallocate(rs)
  if (allocated(zs)) deallocate(zs)
  if (allocated(rbnd)) deallocate(rbnd)
  if (allocated(zbnd)) deallocate(zbnd)

  if (allocated(irins)) deallocate(irins)
  if (allocated(ikins)) deallocate(ikins)
  if (allocated(irouts)) deallocate(irouts)
  if (allocated(ikouts)) deallocate(ikouts)
  if (allocated(idds)) deallocate(idds)

  if (allocated(korpg)) deallocate(korpg)
  if (allocated(nvertp)) deallocate(nvertp)
  if (allocated(rvertp)) deallocate(rvertp)
  if (allocated(zvertp)) deallocate(zvertp)

  if (allocated(btot)) deallocate(btot)
  if (allocated(br)) deallocate(br)
  if (allocated(bz)) deallocate(bz)
  if (allocated(bt)) deallocate(bt)

  if (allocated(knbs)) deallocate(knbs)
  if (allocated(ktebs)) deallocate(ktebs)
  if (allocated(ktibs)) deallocate(ktibs)
  if (allocated(kvhs)) deallocate(kvhs)
  if (allocated(kes)) deallocate(kes)

  if (allocated(knds)) deallocate(knds)
  if (allocated(kteds)) deallocate(kteds)
  if (allocated(ktids)) deallocate(ktids)
  if (allocated(kvds)) deallocate(kvds)
  if (allocated(keds)) deallocate(keds)

end subroutine deallocate_storage




subroutine load_oedge_plasma(infile,ierr)
  implicit none
  integer :: infile,ierr

  !
  !     LOAD_OEDGE_PLASMA:The purpose of this routine is to read in the
  !               DIVIMP background plasma in a DIVIMP specific
  !               format.
  !


  character*256 :: buffer
  integer :: ik,ir,id
  integer :: tmpnrs,tmpnds,tmpirsep
  integer :: tmpnks(maxnrs)
  logical :: done
  !
  !     READ Title line
  !
  read (infile,10) buffer
  !
  if (buffer(1:6).ne.'DIVIMP') then
     call errmsg('NOT A DIVIMP PLASMA FILE',trim(buffer))
     stop 'Not a valid DIVIMP plasma file'
  endif

  done = .false.

  do while (.not.done)

     read(infile,10,iostat=ios,err=error_exit) buffer

     if (ios.eq.0) then 

        if (buffer(1:4).eq.'NRS:') then
           read (buffer,200) tmpnrs,tmpirsep,tmpnds
           !
           !        Check to see if this matches the current grid.
           !
           if (nrs.ne.tmpnrs.or.irsep.ne.tmpirsep.or.nds.ne.tmpnds) then
              !
              !           Grid characteristic mismatch - exit program.
              !
              call errmsg('DIVIMP PLASMA FILE DOES NOT MATCH GRID')
              write (0,*) 'NRS  :',nrs,tmpnrs
              write (0,*) 'IRSEP:',irsep,tmpirsep
              write (0,*) 'NDS  :',nds,tmpnds
              write (0,*) 'PROGRAM EXITING'
              stop 'DIVIMP PLASMA/GRID FILE MISMATCH'

           end if

        elseif(buffer(1:6).eq.'KNOTS:') then

           read (infile,400)  (tmpnks(ir),ir=1,nrs)
           !
           !        Check to see if knots match
           !
           do ir = 1, nrs

              if (nks(ir).ne.tmpnks(ir)) then

                 call errmsg('DIVIMP PLASMA FILE DOES NOT MATCH GRID')
                 write (0,*) 'IR     :',ir
                 write (0,*) 'NKS(IR):',nks(ir),tmpnks(ir)
                 write (0,*) 'PROGRAM EXITING'
                 stop 'DIVIMP PLASMA/GRID FILE MISMATCH'
              end if
           end do

        elseif (buffer(1:5).eq.'KNBS:') then
           !
           !        Density - volume
           !
           read (infile,500) ((knbs(ik,ir),ik=1,nks(ir)),ir=1,nrs)

        elseif (buffer(1:5).eq.'KNDS:') then
           !
           !        Density - target
           !
           read (infile,500) (knds(id),id=1,nds)

        elseif (buffer(1:6).eq.'KTEBS:') then
           !
           !        Te - volume
           !
           read (infile,500) ((ktebs(ik,ir),ik=1,nks(ir)),ir=1,nrs)

        elseif (buffer(1:6).eq.'KTEDS:') then
           !
           !        Te - target
           !
           read (infile,500) (kteds(id),id=1,nds)

        elseif (buffer(1:6).eq.'KTIBS:') then
           !
           !        Ti - volume
           !
           read (infile,500) ((ktibs(ik,ir),ik=1,nks(ir)),ir=1,nrs)

        elseif (buffer(1:6).eq.'KTIDS:') then
           !
           !        Ti - target
           !
           read (infile,500) (ktids(id),id=1,nds)

        elseif (buffer(1:5).eq.'KVHS:') then
           !
           !        Velocity - volume
           !
           read (infile,500) ((kvhs(ik,ir),ik=1,nks(ir)),ir=1,nrs)

        elseif (buffer(1:5).eq.'KVDS:') then
           !
           !        Velocity - target
           !
           read (infile,500) (kvds(id),id=1,nds)

        elseif (buffer(1:4).eq.'KES:') then
           !
           !        Electric Field - volume
           !
           read (infile,500) ((kes(ik,ir),ik=1,nks(ir)),ir=1,nrs)

        elseif (buffer(1:5).eq.'KEDS:') then
           !
           !        Electric Field - target
           !
           read (infile,500) (keds(id),id=1,nds)

        endif

     else
        done=.true.
     endif
     !
     !     Loop back for continued reading
     !
  end do


  CLOSE (infile)

  return

  error_exit:

  call errmsg('ERROR READING IN DIVIMP PLASMA BACKGROUND')
  stop 'ERROR reading in Plasma file'


  !
  !     Formatting
  !
10 format(a)
100 format(a40)
200 format('NRS:',i5,'IRSEP:',i5,'NDS:',i5)
400 format(12i6)
500 format(6e18.10)

end subroutine load_oedge_plasma




subroutine find_free_unit_number(unit)
  implicit none
  integer unit
  !
  !     FIND_FREE_UNIT_NUMBER:
  !
  !     This routine scans through unit numbers looking for one that
  !     is not currently in use. This number is returned. This code
  !     is based on the assumption that any unit numbers returned will
  !     be used before this routine is called again asking for another 
  !     number - otherwise it will likely return the previous value.
  !
  integer test_unit
  logical unit_open

  test_unit = 10
  unit_open = .false.

  ! Check for unit number assignment.  
  Do While (Unit_open)
     test_unit=test_unit + 1
     Inquire (Unit = test_unit, Opened = Unit_open)
  End Do

  unit = test_unit

  return
end subroutine find_free_unit_number



end module oedge_plasma_interface
