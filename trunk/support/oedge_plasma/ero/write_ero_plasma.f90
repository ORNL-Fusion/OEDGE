program write_ero_plasma
  use error_handling
  use oedge_plasma_interface
  implicit none

  character*512 :: arg
  character*256 :: case_name
  integer :: interpolate_option,errmsg_unit
  integer :: ierr
  integer :: ir,iz,nr,nz
  integer :: ounit
  character*256 :: outfilename

  real*8 :: r_offset, z_offset

  real*8 :: rmin,rmax,zmin,zmax 
  real*8 :: dr,dz
  real*8 :: rt,zt,ne,te,ti,vb,ef,btot,br,bz,bt

  integer :: test_cnt,in
  real*8 :: test_r(18)
  real*8 :: test_z(18)

  integer :: narg
  integer,external :: iargc

  !
  ! This program was developed from the plasma_test sample code. This code writes out a plasma file for a grid
  ! covering a specified region. Jeff can read this file into his code. 
  !

  !
  !
  ! Three steps to using the interface
  ! 1) Initialization  -   init_oedge_plasma
  ! 2) Get plasma data -   get_oedge_plasma
  ! 3) Clean up        -   close_oedge_plasma 
  !
  ! Clean up is done when the GET routine will no longer be called - it frees up all allocated arrays
  !
  !

  ! Load command line input - case name, interpolation and error messaging options
  ! 


! read in the ero plasma specification file
! 
! it contains the following information
!
! divimp plasma data file name
! divimp geometry data file name
! number of ero plasas to write out
! ero plasma specification blocks
! file name for storing the ero solution
! 4 vertices expected to be a rectangle
! R,Z offset to map the divimp machine coordinates to the ero simulation coordinates
! grid resolution in the two directions ... directions are defined by vectors and should be perpendicular
! interpolation option to be used 
! error message unit number
! unit conversion flag (do conversion here or in ERO?)
!
! command line argument is the name of the data file to be read
!


! open ero plasma command file

  narg = iargc()

  if (narg.ne.1) then 
     write(0,'(a)') 'Incorrect number of arguments specified:'
     write(0,'(a)') 'USAGE: program_name <ero_plasma_sepcification_file [str]>'
     stop
  endif

  !
  ! Read and assign command line arguments
  !
  ! Command line is:  'file name'   interpolate_option   errmsg_unit
  !
  ! File name is the base case name without the .grd or .bgp file extensions
  !
  call getarg(1,arg)
  case_name = trim(arg)

  if (case_name.eq.'') then 
     ! no arguments specified
     write(0,'(a)') 'No arguments specified on call: This code requires the following input line'
     write(0,'(a)') '<program name>    ''input file name''  '
     stop
  endif


  ! Interpolate option = 0 = off (gives value in cell)
  !                    = 1 = on  (uses interpolation scheme)
  !call getarg(2,arg)
  !read(arg,*) interpolate_option

  !
  ! errmsg_unit specifies the unit number for error message output ... 
  ! A value less than 0 (e.g. -1 ) will suppress all error messages from the oedge_plasma_interface routines
  !
  !call getarg(3,arg)
  !read(arg,*) errmsg_unit

  write(0,'(a,a,a,i8,a,i8)') 'ERO PLASMA: Case = ',trim(case_name)



  integer :: specunit,ios

  call find_free_unit_number(specunit)

  open(specunit,file=trim(case_name),status='old',form='formatted',iostat=ios)

  if (ios.ne.0) then 
      write(0,'(a,a,a,i8,a,i8)') 'ERROR Opening Specification file: Case = ',trim(case_name)
      stop 'Spec file not found'
  endif


  ! read the name of the DIVIMP plasma 

10  read(specunit,'(a512)') line

  call get_plasma_file_name(specunit,plasma_file_name)




  subroutine get_plasma_file_name
    implicit none
    integer :: inunit
    



  ! 
  ! R_offset and Z_offset specify a coordinate system shift to be applied to coordinates passed into the OEDGE plasma routine. 
  ! e.g. If R_offset = 1.5 and Z_offset = -1.366 then a coordinate of 0,0 passed into the OEDGE plasma routine will be mapped to 1.5,-1.366 
  ! This makes it easy to deal with different coordinate systems with different origins. 
  !
  ! However, the OEDGE code works in MKS units ... thus distances and geometry are in meters. If the calling code uses a different coordinate scale (e.g. cm) then 
  ! appropriate conversions must be made before the call to get the plasma conditions. 
  !

  ! set offsets so that the center of DIMES will be at 0.0,0.0 for this simulation. 
  !
  ! The outer strike point on the DIMES Mo grid is at R=1.4864,Z=-1.250
  ! The ACTUAL strike point is at R=1.48043,Z=-1.250  (averaged over all exposure shots from 2500 to 4500ms
  ! The ACTUAL location of the DIMES center is at R=1.4859, Z=-1.25
  !
  ! The ACTUAL strike point for the Mo exposure on May 1st,2012 is R=1.47977m
  ! May 1, 2012 - Mo exposure - DIMES offset from strike point = 0.00613m
  ! Roffset = 1.4864 + 0.00613 = 1.49253 
  ! Difference in offset 0.00016
  !
  !  r_offset = 1.495835
  !  z_offset = -1.25
  ! In order to match simulation and experiment .. the grid strike point of R=1.4864 is matched to the actual average strike point of R=1.48043 
  ! The difference is 5.97e-3m
  ! The difference from the actual strike point to actual DIMES center is : 5.47e-3m
  ! Thus the actual DIMES center in the simulation is 5.47e-3m outboard of the simulation strike point
  ! This puts the DIMES strike point at R = 1.4864 + 5.47e-3 = 1.49187
  !
  ! r_offset = 0.0
  ! z_offset = 0.0
  !
  ! DIMES W Experiment
  !
  ! For this experiment 75% of the exposure was at one position and 25% at another ... so two sets of offsets are required.
  !
  ! Set the Offsets so that the DIMES center is at 0.0, 0.0 in the OEDGE coordinate frame for grid 148679 at 3000ms
  ! These offsets should put the DIMES center at 0,0 for this simulation 
  ! 75%   R_offset = 1.48943  Z_offset = 1.25000
  !
  r_offset =  1.48943
  z_offset = -1.25

  ! 25%   R_offset = 1.50483  Z_offset = 1.25000
  !
  !r_offset =  1.50483
  !z_offset = -1.25

  write(0,'(a,f15.6,a,f15.6)') 'OFFSETS: R_OFFSET=',r_offset, ' Z_OFFSET=',z_offset
  

  call init_oedge_plasma(case_name,r_offset,z_offset,interpolate_option,errmsg_unit,ierr)

  if (ierr.ne.0) then 
     call errmsg('Problem loading oedge plasma',ierr)
     stop 'plasma not loaded'
  endif


  ! Define sample space - call and report plasma conditions
  ! Values chosen here are arbitrary for testing purposes

  outfilename = trim(case_name)//'.plasma.dat'

  call find_free_unit_number(ounit)

  open(ounit,file=trim(outfilename),status='unknown',form='formatted')

  write(0,'(a)') 'Writing plasma to file:'//trim(outfilename)

  nr = 100
  nz = 50

  rmin = -0.05
  rmax = 0.05

  zmin = 0.0
  zmax = 0.05

  dr = (rmax-rmin)/nr
  dz = (zmax-zmin)/nz


  write(ounit,*) 
  write(ounit,'(2(1x,a7),20(7x,a10,2x))') 'R_INDEX','Z_INDEX','R(m)','Z(m)','ne(m-3)','Te(eV)','Ti(eV)','Vpara(m/s)','Epara(V/m))','Btot(T)','B_R','B_Z','B_T'
  write(ounit,*) 

  do ir = 0,nr

     rt = rmin + real(ir) * dr
     !write(0,'(a,i8,10(1x,g18.8))') 'R:',ir,rt,rmin,dr,rmin + ir*dr

     do iz = 0,nz

        zt = zmin + iz * dz

        !write(0,'(a,2i8,10(1x,g18.8))') 'S:',ir,iz,dr,dz,rt,zt,rmin,rmax,zmin,zmax
        !write(6,'(a,2i8,10(1x,g18.8))') 'S:',ir,iz,dr,dz,rt,zt,rmin,rmax,zmin,zmax

        call get_oedge_plasma(rt,zt,ne,te,ti,vb,ef,btot,br,bz,bt,ierr)

        if (ierr.eq.0) then 
           write(ounit,'(2i8,20(1x,g18.8))') ir,iz,rt,zt,ne,te,ti,vb,ef,btot,br,bz,bt
        else
           write(error_message_data,'(a,3i8,20(1x,g18.8))') 'ERROR GETTING PLASMA:',ierr,ir,iz,rt,zt,ne,te,ti,vb,ef,btot,br,bz,bt
           call errmsg('PLASMA TEST',error_message_data)
           !write(6,'(a)') error_message_data

        endif

     end do

  end do

!  write(ounit,*) 
!  write(ounit,'(a,i10)') 'Ordered by Z',interpolate_option
!  write(ounit,*) 


!  do iz = 0,nz
!
!     zt = zmin + iz * dz
!
!     do ir = 0,nr
!
!        rt = rmin + ir * dr
!
!        call get_oedge_plasma(rt,zt,ne,te,ti,vb,ef,btot,br,bz,bt,ierr)
!
!        if (ierr.eq.0) then 
!           write(ounit,'(2i8,20(1x,g18.8))') ir,iz,rt,zt,ne,te,ti,vb,ef,btot,br,bz,bt
!        else
!           write(error_message_data,'(a,3i8,20(1x,g18.8))') 'ERROR GETTING PLASMA:',ierr,ir,iz,rt,zt,ne,te,ti,vb,ef,btot,br,bz,bt
!           call errmsg('PLASMA TEST',error_message_data)
!        endif
!
!     end do
!     write(ounit,*) 
!
!  end do




  ! Clean up when done calling OEDGE plasma interface

  call close_oedge_plasma

  close(ounit)

end program write_dimes_plasma
