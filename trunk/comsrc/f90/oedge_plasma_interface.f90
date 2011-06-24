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
real,allocatable :: btotd(:),brd(:),bzd(:),btd(:)




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


subroutine get_oedge_plasma(r,z,ne,te,ti,vb,ef,btoto,bro,bzo,bto)
  implicit none

  real :: r,z,ne,te,ti,vb,ef,btot,br,bz,bt

  ! Get the OEDGE plasma conditions at the specified R,Z location
  ! Two options - value in cell and interpolated - value in cell is quicker


  ! Get cell and position in cell




  ! Interpolate result if required


  if (interpolate_opt.eq.0) then 
     ! no interpolation
     ne = knbs(ik,ir)
     te = ktebs(ik,ir)
     ti = ktibs(ik,ir)
     vb = kvhs(ik,ir)
     ef = kes(ik,ir)
     btoto = btot(ik,ir)
     bro = br(ik,ir)
     bzo = bz(ik,ir)
     bto = bt(iki,ir)

  elseif (interpolate_opt.eq.1) then
     ! interpolate using Jeff's algorithm
     ! Note: ne,te,ti,vb,ef have target values that are used for interpolation in the first half cell
     !       btot,br,bz,bt do not have target data .. as a result the target values are the same as the first cell center






  endif

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


  ! Assign target data for magnetic field until this is put in transfer file (if needed)
  do ir = irsep,nrs
     btotd(idds(ir,1)) = btot(nks(ir),ir)
     btotd(idds(ir,2)) = btot(1,ir)

     brd(idds(ir,1)) = br(nks(ir),ir)
     brd(idds(ir,2)) = br(1,ir)

     bzd(idds(ir,1)) = bz(nks(ir),ir)
     bzd(idds(ir,2)) = bz(1,ir)

     btd(idds(ir,1)) = bt(nks(ir),ir)
     btd(idds(ir,2)) = bt(1,ir)
  end do


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

  ! Target magnetic field 
  ! Set to values at first cell center
  call allocate_array(btotd,nds,'BTOTD',ierr)
  call allocate_array(brd,nds,'BTOTD',ierr)
  call allocate_array(bzd,nds,'BTOTD',ierr)
  call allocate_array(btd,nds,'BTOTD',ierr)


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
















c
      subroutine gridpos(ik,ir,r,z,newinj,outofgrid)
      implicit none
      integer ik,ir
      real    r,z
      logical newinj,outofgrid
!      include 'params'
!      include 'cgeom'
!      include 'comtor'
!      include 'grbound'
!
!
!     GRIDPOS: This subroutine determines which ik,ir bin
!              the position R,Z is in. It utilizes the
!              checking function INCELL which uses
!              vector cross-products to determine if
!              the test-point is in the given cell.
!
!              Alternatively - for grids for which the cell
!              vertices are not available - it will use the
!              IKXYS,IRXYS values from the older
!              implementation.
!
!              The implicit assumption is that the cell
!              vertices are ordered clockwise.
!
!              If IK and IR are non-zero ... this routine
!              will first check the IK,IR cell and then
!              all immediately neighbouring cells before
!              it starts scanning the entire grid for
!              the bin containing the R,Z position.
!
!              IF an IK,IR bin is not found an error
!              message is issued and the values IR = IRSEP
!              and IK = 1 are returned. The argument
!              OUTOFGRID is set to TRUE.
!              These are the default return values for a
!              target launch. In the case of a wall launch
!              the IK,IR indices returned will be for the
!              centre of the closest outermost ring segment.
!              If newinj = .false. for a wall injection
!              particle that does not fall on the grid then
!              the values of ir and ik passed to the will be
!              routine will be returned if they are valid.
!
!              Note: Under error conditions - when the
!              particle is outside the grid - this routine
!              will return the nearest IN grid bin centre.
!              This means that rings 2, irwall-1 and irtrap+1
!              become the real rings because the rest are
!              additional non-polygons that are usually
!              used by the fluid codes to impose boundary
!              conditions and do not represent valid grid
!              elements. The 1,nks(ir) elements are Ok because
!              the virtual points have usually already been
!              stripped.
!
!              Code which uses this routine MUST be able
!              to deal with outofgrid condition being set. This
!              routine attempts to make reasonable guesses
!              for out of bounds values but these can not
!              be treated as reasonable under all circumstances.
!
!              In addition - in an attempt to make the code
!              more efficient for out of grid cases - this routine
!              has been modified to use the GA15 routines for
!              determining if a point is inside a closed n-sided
!              polygon. There are two used - one defining the outer
!              edge of the grid and the other defining the edge
!              of the core plasma. If the particle is not a new
!              injection then the value of the lastin variable
!              is used to determine the region search strategy.
!              If it is a new wall injection - outside is checked
!              first - then a grid location is searched for.
!
!
!              David Elder, Dec 8 1993
!
!
      external incell
      logical  incell
      real     rsq,zsq,dsq,best,result
      integer  ix,iy,in,jk,jr,ikorg,irorg
      integer  iktmp,iktmp1,irtmp ,iklim,ikdiff
c
c     Initialization
c
      outofgrid = .false.
      irorg = ir
      ikorg = ik
c
c
c     Make initial assumption that newinj particles are IN
c     the grid somewhere.
c
      if (newinj) lastin = 0
c
c
c      write(6,*) 'gridpos:',ik,ir,r,z,newinj,outofgrid,
c     >           cgridopt,xygrid,outgrid,lastin
c
      if (xygrid.eq.0.and.(.not.outgrid)) then
c
c       Check if the particle is in the current cell
c
c        WRITE(50,*) 'GRIDPOS:',ik,ir,lastin,newinj
        if (ik.ne.0.and.ir.ne.0.and.lastin.eq.0) then
           if (incell(ik,ir,r,z)) return
c
c          If it is not in the current cell - check the
c          adjacent cells.
c
           irtmp = ir
           iktmp = ik
c
c          Check same ring - forward and back
c
c           WRITE(50,*) 'GRIDPOS: FOR + BAK'
           if (iktmp.ne.1) then
              ik = iktmp -1
              if (incell(ik,ir,r,z)) return
           endif
           if (iktmp.ne.nks(ir)) then
              ik = iktmp +1
              if (incell(ik,ir,r,z)) return
           endif
c
c          Check next inner ring - same ik - forward and back
c
c           WRITE(50,*) 'GRIDPOS: SAME IK, FOR + BAK'
           ir = irins(iktmp,irtmp)
           ik = ikins(iktmp,irtmp)
           iktmp1 = ik
           if (incell(ik,ir,r,z)) return
           if (iktmp1.ne.1) then
              ik = iktmp1 -1
              if (incell(ik,ir,r,z)) return
           endif
           if (iktmp1.ne.nks(ir)) then
              ik = iktmp1 +1
              if (incell(ik,ir,r,z)) return
           endif
c
c          Check next outer ring - same ik - forward and back
c
c           WRITE(50,*) 'GRIDPOS: NEXT OUTER RING, FOR + BAK'
           ir = irouts(iktmp,irtmp)
           ik = ikouts(iktmp,irtmp)
           iktmp1 = ik
           if (incell(ik,ir,r,z)) return
           if (iktmp1.ne.1) then
              ik = iktmp1 -1
              if (incell(ik,ir,r,z)) return
           endif
           if (iktmp1.ne.nks(ir)) then
              ik = iktmp1 +1
              if (incell(ik,ir,r,z)) return
           endif
c
c          Reset ik,ir to original values
c
           ik = iktmp
           ir = irtmp
         endif
c
c        Perform a search over the entire grid to find cell.
c        This should be optimized by checking the more
c        likely locations first - e.g. near walls and
c        plates. So it will search plates, irwall and irtrap
c        then search the inner plasma regions.
c        In fact - depending on the launch option employed
c        different search options may be chosen. One of
c        which could include checking the original IK,IR of the last
c        particle checked - first.
c
c
c        New particle launches
c
         if (newinj) then
c           WRITE(50,*) 'GRIDPOS: NEWINJ = .TRUE., ENTIRE GRID'
c
c           CNEUTB = 0 check near plates first
c
            if (cneutb.eq.0.and.cneuta.eq.0) then
               lastin = 0
               do 10 ir = irsep,nrs
                  if (float(nks(ir)/2) .eq.float(nks(ir))/2.0) then
                     iklim = nks(ir)/2 - 1
                  else
                     iklim = nks(ir)/2
                  endif
                  do 20 ikdiff = 0,iklim
                     ik = ikdiff + 1
                     if (incell(ik,ir,r,z)) return
                     ik = nks(ir) - ikdiff
                     if (incell(ik,ir,r,z)) return
 20               continue
 10            continue
c
c              Check core plasma
c
               do 30 ir = 1,irsep-1
c
c                 First cell is repeat of last - but incell handles this - relevant to 
c                 S values
c
c                  do 30 ik = 1,nks(ir) -1
c
                  do 30 ik = 1,nks(ir)
c
                     if (incell(ik,ir,r,z)) return
30             continue
c
c              Particle was NOT found in grid -
c              Take default action. For Plate cases
c              this means issuing an error message
c              and determining IK,IR.
c
c
c              Check if ouside plasma region.
c
               CALL GA15B(R,Z,RESULT,IONWPTS,1,iwWORK,4*MAXPTS,
     >              iwINDW,MAXPTS,RIW,ZIW,iwTDUM,iwXDUM,iwYDUM,6)
c
c               write(6,*) 'gp:ga15b1:',result,r,z
c
c              If < 0 particle is outside last ring.
c
               if (result.lt.0.0) then
                  lastin = -1
                  outofgrid = .true.
                  ik = 1
                  ir = irsep
c
c                  write(6,'(a,2i6,2(1x,g12.5)') 'GRIDPOS:'//
c     >                 ' PARTICLE OFF GRID (1):',ikorg,irorg,r,z
c
c              Check if in core
c
               else
c
                  CALL GA15B(R,Z,RESULT,IONcPTS,1,icWORK,4*MAXPTS,
     >              icINDW,MAXPTS,RCW,ZCW,icTDUM,icXDUM,icYDUM,6)
c
                  if (result.gt.0.0) then
c
c                    Particle in core
c
                     outofgrid = .true.
                     lastin =1
                     ik = 1
                     ir = 2
c                     write(6,'(a,2i6,2(1x,g12.5)') 'GRIDPOS:'//
c     >                 ' PARTICLE OFF GRID (2):',ikorg,irorg,r,z
                  else
c
c                    Particle is in a small region
c                    between polygon edges and the
c                    ionwall defined above - because
c                    they use different methods - find
c                    the closest ik,ir indices and
c                    return those - hopefully this
c                    will not be required often.
c
                     lastin = 0
                     outofgrid = .false.
c
                     call findwall(ik,ir,r,z)
c
c                     write(6,'(a,4i6,2(1x,g12.5)') 'GRIDPOS:'//
c     >                 ' PARTICLE OFF GRID (3):',ikorg,irorg,ik,ir,r,z
c
                  endif
               endif
               return
            elseif (cneutb.eq.2.and.cneuta.eq.0) then
c
c              New injection for wall launch - check perimeter
c              cells first. 
c
               lastin = 0
c
c              Check if ouside plasma region.
c
               CALL GA15B(R,Z,RESULT,IONWPTS,1,iwWORK,4*MAXPTS,
     >              iwINDW,MAXPTS,RIW,ZIW,iwTDUM,iwXDUM,iwYDUM,6)
c
c              If < 0 particle is outside last ring.
c
               if (result.lt.0.0) then
                  lastin = -1

c                 Particle not found in grid - since it is
c                 a new injection for a wall launch ... find
c                 the nearest wall bin to the injection
c                 position and set outofgrid.
c
                  outofgrid = .true.
c
                  call findwall(ik,ir,r,z)
c
c                  write(6,'(a,4i6,2(1x,g12.5)') 'GRIDPOS:'//
c     >             ' PARTICLE OFF GRID (4):',ikorg,irorg,ik,ir,r,z
c
                  return
               endif
c
c              Check rest of plasma for particle
c
               ir = irwall-1
               do 40 ik = 1,nks(ir)
                  if (incell(ik,ir,r,z)) return
 40            continue

c
c              jdemod - check for existence of PFZ before checking trap rings
c
               if (.not.nopriv) then
                  ir = irtrap+1
                  do 50 ik = 1,nks(ir)
                     if (incell(ik,ir,r,z)) return
 50               continue
               endif

c
c              Check rest of grid
c
               do 60 ir = irwall -2,irsep,-1
                  do 60 ik = 1,nks(ir)
                     if (incell(ik,ir,r,z)) return
 60            continue
c
c              jdemod - check for existence of PFZ before checking trap rings
c
               if (.not.nopriv) then 
                  do 70 ir = irtrap+2 ,nrs
                     do 70 ik = 1,nks(ir)
                        if (incell(ik,ir,r,z)) return
 70               continue
               endif
c
c              Check core
c
               do 80 ir = 1,irsep-1
c                  do 80 ik = 1,nks(ir)-1
                  do 80 ik = 1,nks(ir)
                     if (incell(ik,ir,r,z)) return
 80            continue
c
c              Particle not found in grid - since it is
c              a new injection for a wall launch ... find
c              if it is in the core OR find nearest
c              boundary bin.
c
c
c              Check if actually in core
c
               CALL GA15B(R,Z,RESULT,IONcPTS,1,icWORK,4*MAXPTS,
     >              icINDW,MAXPTS,RcW,ZcW,icTDUM,icXDUM,icYDUM,6)
c
               if (result.gt.0.0) then
c
c                 Particle in core
c
                  outofgrid = .true.
                  lastin =1
                  ik = 1
                  ir = 2
c                  write(6,'(a,2i6,2(1x,g12.5)') 'GRIDPOS:'//
c     >              ' PARTICLE OFF GRID (5):',ikorg,irorg,r,z
               else
c
c                 Particle is in a small region
c                 between polygon edges and the
c                 ionwall defined above - because
c                 they use different methods - find
c                 the closest ik,ir indices and
c                 return those - hopefully this
c                 will not be required often.
c
                  lastin = 0
                  outofgrid = .false.
c
                  call findwall(ik,ir,r,z)
c
                  write(6,'(a,4i6,2(1x,g12.5))') 'GRIDPOS:'//
     >                 ' PARTICLE OFF GRID (6):',
     >                  ikorg,irorg,ik,ir,r,z
c

               endif
               return
c
            else
c
c            elseif ((cneutb.eq.1.and.cneuta.eq.0).or.cneuta.eq.1) then
c
c              This was originally the intended strategy for the
c              conditions in the commented out elseif - but the
c              strategy seems to be a reasonable one for the
c              default case.
c
c              Check the last original start position first
c              just in case all the particles are originating
c              in the same bin.
c
               lastin = 0
c
               if ((irstold.ge.1.and.irstold.le.nrs).and.
c slmod begin
c...Problem for Intel compiler (non-starndard FORTRAN):
     >             (ikstold.ge.1.and.
     >              ikstold.le.nks(MAX(1,irstold)))) then
c
c     >             (ikstold.ge.1.and.ikstold.le.nks(irstold))) then
c slmod end
                  ik = ikstold
                  ir = irstold
                  if (incell(ik,ir,r,z)) return
               endif
c
c              Scan the full grid and set ikstold and irstold
c
c              Scan the SOL first - since injection is more likely
c              in this region.
c
               do 100 ir = irsep,nrs
                  do 100 ik = 1,nks(ir)
                     if (incell(ik,ir,r,z)) then
                        ikstold = ik
                        irstold = ir
                        return
                     endif
100            continue
c
c              Scan the core
c
               do 110 ir = 1,irsep-1
c
c                 First and last cell of core are the same except for S coordinate issue -
c                 Incell now checks which part of first or last cell to determine particle
c                 prescence.   
c
c                 do 110 ik = 1,nks(ir)-1
c
                  do 110 ik = 1,nks(ir)
c
                     if (incell(ik,ir,r,z)) then
                        ikstold = ik
                        irstold = ir
                        return
                     endif
110            continue
c
c              Error condition - particle may not be in grid
c
c
c              Check if ouside plasma region.
c
               CALL GA15B(R,Z,RESULT,IONWPTS,1,iwWORK,4*MAXPTS,
     >              iwINDW,MAXPTS,RIW,ZIW,iwTDUM,iwXDUM,iwYDUM,6)
c
c               write(6,*) 'gp:ga15b3:',ionwpts,result,r,z
c
c              If result < 0 particle is outside last ring.
c
               if (result.lt.0.0) then
                  lastin = -1
                  outofgrid = .true.
                  ir = irsep
                  ik = 1
c                  write(6,'(a,2i6,2(1x,g12.5)') 'GRIDPOS:'//
c     >               ' PARTICLE OFF GRID (7):',ikorg,irorg,r,z
c
c              Otherwise may be in core
c
               else

c
c                 Check if actually in core
c
                  CALL GA15B(R,Z,RESULT,IONcPTS,1,icWORK,4*MAXPTS,
     >               icINDW,MAXPTS,RcW,ZcW,icTDUM,icXDUM,icYDUM,6)
c
                  if (result.gt.0.0) then
c
c                    Particle in core
c
                     outofgrid = .true.
                     lastin =1
                     ik = 1
                     ir = 2
c                     write(6,'(a,2i6,2(1x,g12.5)') 'GRIDPOS:'//
c     >               ' PARTICLE OFF GRID (8):',ikorg,irorg,r,z
                  else
c
c                    Particle is in a small region
c                    between polygon edges and the
c                    ionwall defined above - because
c                    they use different methods - find
c                    the closest ik,ir indices and
c                    return those - hopefully this
c                    will not be required often.
c
                     lastin = 0
                     outofgrid = .false.
c
                     call findwall(ik,ir,r,z)
c                     write(6,'(a,4i6,2(1x,g12.5)') 'GRIDPOS:'//
c     >                  ' PARTICLE OFF GRID (9):',
c     >                     ikorg,irorg,ik,ir,r,z
c
                  endif
               endif
c
               return
c
            endif
         elseif (.not.newinj) then
c           WRITE(50,*) 'GRIDPOS: NEWINJ = .FALSE., ENTIRE GRID',lastin
c
c           Not a new injection - but it has fallen through to
c           this point - This means either it was not near the
c           last recorded bin or the last recorded bin was not set.
c           i.e. IK = IR = 0 at the start of the routine.
c
c           Since this routine does not scan for ions stepping along
c           the field lines - only initial injections - and neutrals
c           are not expected to take large arbitrary steps :-) - There
c           is no reason to start at any specific place. So scan the
c           SOL and CORE grid in that order and set any appropriate
c           error conditions.
c
c           Check the lastin region - if the particle is still
c           outside the grid region the one can return quickly.
c
c           One does not expect the particle to move from the
c           core to outside the plasma in one step ... so if
c           the particle has changed regions - check for it in the
c           grid.
c
            if (lastin.eq.-1) then
c
c              Check if particle is still outside.
c              return the values passed to the subroutine.
c              Otherwise scan the grid region.
c
               CALL GA15B(R,Z,RESULT,IONWPTS,1,iwWORK,4*MAXPTS,
     >              iwINDW,MAXPTS,RIW,ZIW,iwTDUM,iwXDUM,iwYDUM,6)
c
c               write(6,*) 'gp:ga15b4:',ionwpts,result,r,z
c
c              If < 0 particle is outside last ring.
c
               ik = ikorg
               ir = irorg
c               WRITE(50,*) 'GRIDPOS: POLYGON CHECK',result
c               WRITE(50,*) 'GRIDPOS:              ',ionwpts,maxpts
c               WRITE(50,*) 'GRIDPOS:              ',r,z
               if (result.lt.0.0) then
                  outofgrid = .true.
c                  write(6,'(a,2i6,2(1x,g12.5)') 'GRIDPOS:'//
c     >               ' PARTICLE OFF GRID (10):',ikorg,irorg,r,z
                  return
               endif
c
            elseif (lastin.eq.1) then
c
               CALL GA15B(R,Z,RESULT,IONcPTS,1,icWORK,4*MAXPTS,
     >              icINDW,MAXPTS,RcW,ZcW,icTDUM,icXDUM,icYDUM,6)
c
c              If > 0 particle is inside core ring.
c
c               write(6,*) 'gp:ga15b5:',ioncpts,result,r,z
c
               ik = ikorg
               ir = irorg
               if (result.gt.0.0) then
                  outofgrid = .true.
c                  write(6,'(a,2i6,2(1x,g12.5)') 'GRIDPOS:'//
c     >             ' PARTICLE OFF GRID (11):',ikorg,irorg,r,z
                  return
               endif
            endif
c
c
c           Assume particle is now in grid.
c
            lastin = 0
c
c           Scan the SOL first - since injection is more likely
c           in this region.
c
c            WRITE(50,*) 'GRIDPOS: SCAN SOL+PFZ'
            do 120 ir = irsep,nrs
               do 120 ik = 1,nks(ir)
                  if (incell(ik,ir,r,z)) return
120         continue
c
c           Scan the core
c
c            WRITE(50,*) 'GRIDPOS: SCAN CORE'
            do 130 ir = 1,irsep-1
c               do 130 ik = 1,nks(ir)-1
               do 130 ik = 1,nks(ir)
                  if (incell(ik,ir,r,z)) return
130         continue
c
c
c           Error condition - particle may not be in grid
c
c
c           Check if ouside plasma region.
c
c            WRITE(50,*) 'GRIDPOS: CHECK IF OUTSIDE PLASMA REGION'

            CALL GA15B(R,Z,RESULT,IONWPTS,1,iwWORK,4*MAXPTS,
     >              iwINDW,MAXPTS,RIW,ZIW,iwTDUM,iwXDUM,iwYDUM,6)
c
c            write(6,*) 'gp:ga15b6:',ionwpts,result,r,z
c
c           If result < 0 particle is outside last ring.
c
c            WRITE(50,*) 'GRIDPOS: CHECKING POLYGON AGAIN',result
c            WRITE(50,*) 'GRIDPOS:                       ',r,z
            if (result.lt.0.0) then
               lastin = -1
               outofgrid = .true.


               if (cneutb.eq.2.and.cneuta.eq.0) then
c
c                 If a wall launch that is not in the grid
c                 and not a new injection then return the
c                 ik and ir values that were passed to
c                 the routine - if they are valid.
c
                  if ((irorg.ge.1.and.irorg.le.nrs).and.
     >               (ikorg.ge.1.and.ikorg.le.nks(irorg))) then

                     ir = irorg
                     ik = ikorg
                  else
                     ir = irsep
                     ik = 1
                  endif
               else
c
c                 Default failure values
c
                  ir = irsep
                  ik = 1
               endif

c               write(6,'(a,2i6,2(1x,g12.5)') 'GRIDPOS:'//
c     >               ' PARTICLE OFF GRID (12):',
c     >                       ikorg,irorg,r,z


c
c           Otherwise may be in core
c
            else
c
c              Check if actually in core
c
               CALL GA15B(R,Z,RESULT,IONcPTS,1,icWORK,4*MAXPTS,
     >               icINDW,MAXPTS,RcW,ZcW,icTDUM,icXDUM,icYDUM,6)
c
c               WRITE(50,*) 'GRIDPOS: IN CORE?',result
c
c               write(6,*) 'gp:ga15b8:',ioncpts,result,r,z
c
               if (result.gt.0.0) then
c
c                 Particle in core
c
                  outofgrid = .true.
                  lastin =1
                  ik = 1
                  ir = 2
c
c                  write(6,'(a,2i6,2(1x,g12.5)') 'GRIDPOS:'//
c     >               ' PARTICLE OFF GRID (13):',
c     >                ikorg,irorg,r,z
c
               else
c
c                 Particle is in a small region
c                 between polygon edges and the
c                 ionwall defined above - because
c                 they may use different methods - find
c                 the closest ik,ir indices and
c                 return those - hopefully this
c                 will not be required often.
c
c                  WRITE(50,*) 'GRIDPOS: DISPIRATION SET IN...'
                  lastin = 0
                  outofgrid = .false.
c
                  call findwall(ik,ir,r,z)
c
c                  write(6,'(a,4i6,2(1x,g12.5)') 'GRIDPOS:'//
c     >               ' PARTICLE OFF GRID (14):',
c     >                ikorg,irorg,ik,ir,r,z

c
               endif
            endif
c
            return
c
         endif
      elseif (xygrid.eq.1.or.outgrid) then
        IX   = MAX (1, MIN (NXS, INT((R-RMIN)/DR)+1))
        IY   = MAX (1, MIN (NYS, INT((Z-ZMIN)/DZ)+1))
        call gridcoords(ix,iy,ik,ir,in)
c
c        IK   = IKXYS(IX,IY)
c        IR   = IRXYS(IX,IY)
c
c       Set outofgrid to true if particle is outside of the
c       defined grid region - even if one is using the
c       X,Y grid.
c
        if (in.ne.1) outofgrid = .true.
c
      endif

      return
      end
c



      logical function incell(ik,ir,r,z)
      implicit none
      integer ik,ir
      real r,z
!      include 'params'
!      include 'cgeom'
!
!     INCELL: This function returns a simple YES/NO decision
!             about whether the point R,Z is in the cell designated
!             by IK,IR with a set of vertices defined in an ordered
!             clockwise fashion. It takes the cross product
!             between the vector from the vertex to the test point
!             and the vector from the vertex to the next clockwise
!             vertex of the polygon. The cross-product must be
!             the same sign for all vertices - if the
!             point is outside the polygon it will fail this test
!             for at least one vertex. (i.e. the cross-product will
!             be less than zero.) (Suggested solution courtesy
!             of Ian Youle :-) )
!
!             David Elder, Dec 8, 1993
!
!             Note: the objectives of the solution method were
!             simplicity and reasonable computational cost.
!             This solution avoids the need for square roots
!             or trigonometric calculations.
!
!             Note: in the confined plasma the first and last cells
!                   are identical. However, S=0 and S=SMAX are at the 
!                   center of this cell. This causes some inconsistencies
!                   when calculating particle positions. In order
!                   to address this - this routine will consder a paricle
!                   in the first cell on a core ring when it is in the
!                   second half of the cell and in the last cell of a core
!                   ring when it is in the first half of the cell. 
!
      integer k,v,nextv,i,nv
      real vxr,vxz,vwr,vwz,cp,lastcp
!
      logical inpoly,res
      external inpoly
      real rc(4),zc(4)
c
      lastcp = 0.0
c
      incell = .false.
      k = korpg(ik,ir)
      if (k.eq.0) return
      nv = nvertp(k)
      if (nv.eq.0) return
      do 10 v = 1,nv
         if (v.eq.nv) then
            nextv = 1
         else
            nextv = v+1
         endif
c
c        Want the vector cross-product Rx X Rw
c
c         vxr = r - rvertp(v,k)
c         vxz = z - zvertp(v,k)
c         vwr = rvertp(nextv,k) - rvertp(v,k)
c         vwz = zvertp(nextv,k) - zvertp(v,k)
c
c         cp = vxr*vwz - vxz*vwr
c
c         if (cp.lt.0.0)  return
c
c          if (   (
c     >     ( (r-rvertp(v,k)) *
c     >       (zvertp(nextv,k)-zvertp(v,k)) )
c     >    -( (z-zvertp(v,k)) *
c     >       (rvertp(nextv,k)-rvertp(v,k)) )
c     >           )
c     >         .lt.0.0) return
c
          cp =    (
     >     ( (r-rvertp(v,k)) *
     >       (zvertp(nextv,k)-zvertp(v,k)) )
     >    -( (z-zvertp(v,k)) *
     >       (rvertp(nextv,k)-rvertp(v,k)) )
     >           )
c
c         There is a problem for points that should 
c         lie on the boundary of the cell - i.e. that 
c         are calculated based on the polygon corners and 
c         which are mathematically on the polygon surface. 
c         Numerically, these points can have a cross product
c         which is close to zero but can vary to either side. 
c         In order to consider these points in the cell - the 
c         cross products are set to zero for values less than
c         a specified limit. In this case the limit is set to 1.0e-7 
c
c         This value was determined by examining the range of cross 
c         product values generated when sampling 50,000 points 
c         calculated on a polygon with a scale size of 1.0m. 
c         The maximum error cross product in this case was 6e-8.
c
c         D. Elder, Dec 13, 2006
c
c         Upon consideration - it might be best to not allow these points
c         to be considered inside the cell since if they are detected they
c         can be moved slightly to an appropriate location. 
c
c          if (abs(cp).lt.1.0e-7) cp = 0.0 
c
          if (v.eq.1.and.cp.ne.0.0) lastcp = cp
c
          if ((lastcp * cp).lt.0.0) return
c
          if (cp.ne.0.0) lastcp = cp
c
10    continue
c
c     Particle has been found in cell
c
      incell = .true.
c
c     Check particles in the first or last cell of core rings
c     for more accurate assessement.
c
      if (ir.lt.irsep.and.
     >   (ik.eq.1.or.ik.eq.nks(ir))) then 
c
c        For a particle found to be in the first or last cell of 
c        a core ring - need to decide if it is in the first
c        half or second half and revise incell result accordingly. 
c
c        Only the second half of the first cell or the first half
c        of the last cell should return true
c
c
c        NOTE: Keep in mind that the first and last cells of core
c              rings are supposed to be identical - if this changes
c              then this code needs to be modified. 
c         
c        Check to see if particle is in the first half of the 
c        cell - set vertices for call to inpoly.-  
c         
c        k was set to cell geometry index at the beginning of this
c        routine.
c
         rc(1) = rvertp(1,k)
         rc(2) = rvertp(2,k)
         rc(3) = (rvertp(2,k) + rvertp(3,k)) /2.0
         rc(4) = (rvertp(1,k) + rvertp(4,k)) /2.0
c
         zc(1) = zvertp(1,k)
         zc(2) = zvertp(2,k)
         zc(3) = (zvertp(2,k) + zvertp(3,k)) /2.0
         zc(4) = (zvertp(1,k) + zvertp(4,k)) /2.0
c
c        Check to see of point is in first half of cell
c
         res = inpoly(r,z,4,rc,zc)
c
c        Change value of incell to false for either of the 
c        invalid cases
c        First half and in first cell
c        Last half and in last cell. 
c
c
         if ((ik.eq.1.and.res).or.
     >       (ik.eq.nks(ir).and.(.not.res))) then 
            incell = .false.
         endif

      endif


      return
      end


      logical function inpoly(r,z,nv,rvert,zvert)
      implicit none
      integer ik,ir
      integer nv,maxvert   
!      parameter (maxvert=8)
      real r,z,rvert(nv),zvert(nv)
!
!     INPOLY: This function returns a simple YES/NO decision
!             about whether the point R,Z is in the cell designated
!             by a set of vertices defined in an ordered fashion.
!             It takes the cross product
!             between the vector from the vertex to the test point
!             and the vector from the vertex to the next
!             vertex of the polygon. The cross-product must be
!             the same sign for all vertices - if the
!             point is outside the polygon it will fail this test
!             for at least one vertex. (i.e. the cross-product will
!             change sign) (Suggested solution courtesy
!             of Ian Youle :-) )
!
!             David Elder, Dec 8, 1993
!
!             Note: the objectives of the solution method were
!             simplicity and reasonable computational cost.
!             This solution avoids the need for square roots
!             or trigonometric calculations.
!
      integer v,nextv
      real cp,lastcp

      lastcp = 0.0 

      inpoly = .false.

      if (nv.eq.0) return  
!
!     Loop through vertices
!
      do v = 1,nv
!
         if (v.eq.nv) then
            nextv = 1
         else
            nextv = v+1
         endif
!
!        Want the vector cross-product Rx X Rw
!
!         vxr = r - rvert(v)
!         vxz = z - zvert(v)
!         vwr = rvert(nextv) - rvert(v)
!         vwz = zvert(nextv) - zvert(v)
!
!         cp = vxr*vwz - vxz*vwr
!
          cp =    (((r-rvert(v)) * (zvert(nextv)-zvert(v))) &
                 - ((z-zvert(v)) * (rvert(nextv)-rvert(v))))
!
!         There is a problem for points that should 
!         lie on the boundary of the cell - i.e. that 
!         are calculated based on the polygon corners and 
!         which are mathematically on the polygon surface. 
!         Numerically, these points can have a cross product
!         which is close to zero but can vary to either side. 
!         In order to consider these points in the cell - the 
!         cross products are set to zero for values less than
!         a specified limit. In this case the limit is set to 1.0e-7 
!
!         This value was determined by examining the range of cross 
!         product values generated when sampling 50,000 points 
!         calculated on a polygon with a scale size of 1.0m. 
!         The maximum error cross product in this case was 6e-8.
!
!         D. Elder, Dec 13, 2006
!
          if (abs(cp).lt.1.0e-7) cp = 0.0 

          if (v.eq.1.and.cp.ne.0.0) lastcp = cp
!
!         Look for change in sign of cp  
!
          if ((lastcp * cp).lt.0.0) return 
!
          if (cp.ne.0.0) lastcp = cp  
!
      end do
c
      inpoly = .true.
      return
    end function inpoly



      recursive subroutine position_in_poly(r,z,nvert,rvert,zvert,s_frac,cross_frac,base_frac,iter,maxiter,ierr)
      implicit none
      integer nvert,iter,maxiter,ierr

      real r,z,rvert(nvert),zvert(nvert)

      real s_frac,cross_frac,base_frac 
!
!
!
!     POSITION_IN_POLY: This routine is invoked recursively - the 
!                       intention is to determine roughly at what fraction
!                       of the way across the two axes of the cell the
!                       input point R,Z lies. It does this by recursively 
!                       dividing the cell into quarters and iterating 
!                       a specified number of times. The point R,Z is always
!                       within the polygon that is passed on to the 
!                       routine when it invokes itself. The s_frac and 
!                       cross_frac values record the cumulative relative
!                       position on the "s" and "cross" axes - (though 
!                       "cross" is not an accurate description of the 
!                       second axis for non-orthogonal cells.) This
!                       routine can be iterated as long as desired to 
!                       obtain any desired level of accuracy.
!
!
!    Local variables 
!
!
!
      integer nv,nvmax,iv,ivf,ivnext,ivlast,in
      parameter(nvmax=4)
      real rv(nvmax),zv(nvmax),rs(nvmax),zs(nvmax),rcp,zcp
      logical found,inpoly
      external inpoly  
!
!     Code is designed to work for 4-sided polygons only.
!  
      if (nvert.ne.4) return

!
!     s_frac and cross_frac both start at 0.0
! 
!     base_frac starts initially at 1.0 when passed in
!       
      base_frac = base_frac/2.0
!
!     Split the polygon into 4 pieces and determine which part of the
!     cell the point R,Z lies in - adjust 
!
      nv = nvert

      do iv = 1,nv

         ivnext = iv+1
         if (iv.eq.nv) ivnext = 1 

         rs(iv) = (rvert(iv)+rvert(ivnext))/2.0
         zs(iv) = (zvert(iv)+zvert(ivnext))/2.0

      end do

      rcp = (rs(1) + rs(3))/2.0
      zcp = (zs(1) + zs(3))/2.0

      ivf= 0
      iv = 1
      found = .false. 

      do while(iv.le.4.and.(.not.found)) 
         ivlast = iv-1
         if (iv.eq.1) ivlast = 4 
!
!        Determine corners of polygon to check 
!
         rv(mod(iv-1,4)+1)   = rvert(iv) 
         rv(mod(1+iv-1,4)+1) = rs(iv)
         rv(mod(2+iv-1,4)+1) = rcp
         rv(mod(3+iv-1,4)+1) = rs(ivlast)

         zv(mod(iv-1,4)+1)   = zvert(iv) 
         zv(mod(1+iv-1,4)+1) = zs(iv)
         zv(mod(2+iv-1,4)+1) = zcp
         zv(mod(3+iv-1,4)+1) = zs(ivlast)

         found = inpoly(r,z,nv,rv,zv) 

         if (found) ivf = iv

         iv = iv+1

      end do

!
!     Check to make sure location found
!
      if (ivf.ne.0) then 
!
!        The definition of the +/- S and CROSS axes are set here
!        to match the conventions used in DIVIMP relative to the sides of the
!        cell. 
!
!        Side 41 (between polygon corners 1 and 4) is considered INWARD and 
!                corresponds to a positive CROSS displacement.
!        Side 23 (between polygon corners 2 and 3) is considered OUTWARD and 
!                corresponds to a negative CROSS displacement.
!        Side 34 (between polygon corners 3 and 4) is considered UP (larger
!                value of S along the field line) and 
!                corresponds to a positive S displacement (relative to the 
!                cell center)
!        Side 12 (between polygon corners 1 and 2) is considered DOWN (smaller
!                value of S along the field line) and 
!                corresponds to a negative S displacement (relative to the 
!                cell center)
!
         if (ivf.eq.1.or.ivf.eq.2) then 
            s_frac = s_frac - base_frac
         else
            s_frac = s_frac + base_frac
         endif

         if (ivf.eq.1.or.ivf.eq.4) then 
            cross_frac = cross_frac + base_frac
         else
            cross_frac = cross_frac - base_frac
         endif

         if (iter.eq.maxiter) then 
            return
         else
!
!           Increment iteration
!
            iter = iter + 1
         
            call position_in_poly(r,z,nv,rv,zv,s_frac,cross_frac,base_frac,iter,maxiter,ierr)   

         endif

      else

!        Error    
 
         ierr = 1

         write(6,'(a,8(1x,g12.5))') 'ERROR in "position_in_poly": point not found in cell: ITER=',iter
         write(6,'(a,8(1x,g12.5))') 'Last poly:',((rv(iv),zv(iv)),iv=1,4)
         write(6,'(a,8(1x,g12.5))') 'Point R,Z,S,C:',r,z,s_frac,cross_frac

      endif

      return 
    end subroutine position_in_poly






end module oedge_plasma_interface
