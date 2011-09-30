program plasma_test
  use error_handling
  use oedge_plasma_interface
  implicit none

  character*512 :: arg
  character*256 :: case_name
  integer :: interpolate_option,errmsg_unit
  integer :: ierr
  integer :: ir,iz,nr,nz
  real*8 :: r_offset, z_offset

  real*8 :: rmin,rmax,zmin,zmax 
  real*8 :: dr,dz
  real*8 :: rt,zt,ne,te,ti,vb,ef,btot,br,bz,bt

  integer :: test_cnt,in
  real*8 :: test_r(18)
  real*8 :: test_z(18)

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
     write(0,'(a)') '<program name>    ''input file case name''   interpolate_option    errmsg_unit'
     stop
  endif


  ! Interpolate option = 0 = off (gives value in cell)
  !                    = 1 = on  (uses interpolation scheme)
  call getarg(2,arg)
  read(arg,*) interpolate_option

  !
  ! errmsg_unit specifies the unit number for error message output ... 
  ! A value less than 0 (e.g. -1 ) will suppress all error messages from the oedge_plasma_interface routines
  !
  call getarg(3,arg)
  read(arg,*) errmsg_unit

  write(0,'(a,a,a,i8,a,i8)') 'OEDGE PLASMA: Case = ',trim(case_name),'  Interpolate Option =',interpolate_option,'    Error message output unit =',errmsg_unit


  ! 
  ! R_offset and Z_offset specify a coordinate system shift to be applied to coordinates passed into the OEDGE plasma routine. 
  ! e.g. If R_offset = 1.5 and Z_offset = -1.366 then a coordinate of 0,0 passed into the OEDGE plasma routine will be mapped to 1.5,-1.366 
  ! This makes it easy to deal with different coordinate systems with different origins. 
  !
  ! However, the OEDGE code works in MKS units ... thus distances and geometry are in meters. If the calling code uses a different coordinate scale (e.g. cm) then 
  ! appropriate conversions must be made before the call to get the plasma conditions. 
  !
  r_offset = 0.0
  z_offset = 0.0

  call init_oedge_plasma(case_name,r_offset,z_offset,interpolate_option,errmsg_unit,ierr)

  if (ierr.ne.0) then 
     call errmsg('Problem loading oedge plasma',ierr)
     stop 'plasma not loaded'
  endif


  ! Specific test data at known points (cell centers)
  ! for case d-105500-hrc-srcs-a1
  !  IK  IR    R           Z          BPHI     TEB     TIB     NB      E1      VB          S      BTOT/BTHETA    FEG1    FIG1
  !-----------------------------------------------------------------------------------------------------------------------------------



  !   76  27  1.53402436 -1.33341038  2.308    29.1    37.5 1.5E+19-2.17E+00 9.58E+03 66.088394       21.87   -8.67E-01-1.65E+01
  !   77  27  1.53788400 -1.34299994  2.302    28.8    35.0 1.5E+19-2.21E+00 1.08E+04 66.312065       21.40   -8.62E-01-1.88E+01
  !   78  27  1.54140437 -1.35168374  2.297    28.6    32.4 1.6E+19-9.82E-01 1.27E+04 66.510849       21.00   -8.43E-01-2.10E+01
  !   79  27  1.54426229 -1.35869575  2.292    28.4    30.1 1.6E+19 1.45E+01 1.57E+04 66.668907       20.71   -7.68E-01-2.11E+01
  !   80  27  1.54643154 -1.36398411  2.289    28.3    28.6 1.4E+19 2.88E+01 2.62E+04 66.786804       20.51   -7.09E-01-2.00E+01


  !   74  24  1.50918984 -1.32066262  2.346    35.7    50.8 2.3E+19-2.16E+00 9.39E+03 78.480873       23.34   -1.11E+00-1.72E+01
  !   75  24  1.51333034 -1.32989371  2.339    35.4    48.0 2.4E+19-2.72E+00 9.86E+03 78.714127       22.74   -1.13E+00-1.95E+01
  !   76  24  1.51705146 -1.33820403  2.334    35.0    45.3 2.4E+19-3.29E+00 1.05E+04 78.919022       22.24   -1.14E+00-2.22E+01
  !   77  24  1.52047384 -1.34587097  2.328    34.7    42.5 2.5E+19-3.70E+00 1.13E+04 79.104095       21.83   -1.14E+00-2.53E+01

  ! separatrix ring on this grid (outer strike point is next to bottom cell)
  !   76  22  1.51093364 -1.34098816  2.343    34.7    43.7 2.3E+19-3.28E+00 1.13E+04 103.29485       22.31   -1.08E+00-2.20E+01
  !   77  22  1.51422942 -1.34807062  2.338    34.4    41.1 2.4E+19-3.57E+00 1.22E+04 103.46770       21.92   -1.07E+00-2.49E+01
  !   78  22  1.51710665 -1.35428226  2.334    34.2    38.6 2.5E+19-2.67E+00 1.35E+04 103.61670       21.60   -1.06E+00-2.80E+01
  !   79  22  1.51959741 -1.35968626  2.330    34.0    36.2 2.5E+19 1.69E+01 1.58E+04 103.74448       21.34   -9.74E-01-2.83E+01
  !   80  22  1.52169108 -1.36424339  2.326    33.9    34.3 2.3E+19 3.55E+01 2.73E+04 103.85102       21.13   -9.06E-01-2.71E+01

  ! last ring in PFZ - separatrix is outboard of bottom cell
  !   37  73  1.51250100 -1.34896624  2.341    32.5    38.6 2.2E+19-3.17E+00 1.26E+04 41.143063       21.93   -1.05E+00-2.40E+01
  !   38  73  1.51530445 -1.35495400  2.336    32.3    36.3 2.3E+19-2.04E+00 1.39E+04 41.287086       21.61   -1.03E+00-2.69E+01
  !   39  73  1.51766908 -1.36003101  2.333    32.1    34.1 2.3E+19 1.76E+01 1.62E+04 41.407513       21.37   -9.53E-01-2.71E+01
  !   40  73  1.51966643 -1.36433446  2.330    32.0    32.4 2.1E+19 3.61E+01 2.72E+04 41.508469       21.17   -8.86E-01-2.59E+01


  ! Sample some of these points to test return values and make sure they match expected data

  test_cnt = 18
  test_r =  [ 1.53402436, 1.53788400, 1.54140437, 1.54426229, 1.54643154, 1.50918984, 1.51333034, 1.51705146, 1.52047384,&
       1.51093364, 1.51422942, 1.51710665, 1.51959741, 1.52169108, 1.51250100, 1.51530445, 1.51766908, 1.51966643]
  test_z =  [-1.33341038,-1.34299994,-1.35168374,-1.35869575,-1.36398411,-1.32066262,-1.32989371,-1.33820403,-1.34587097,&
       -1.34098816,-1.34807062,-1.35428226,-1.35968626,-1.36424339,-1.34896624,-1.35495400,-1.36003101,-1.36433446] 

  write (6,'(a)') 'Test Results:'

  do in = 1,test_cnt

     call get_oedge_plasma(test_r(in),test_z(in),ne,te,ti,vb,ef,btot,br,bz,bt,ierr)

     if (ierr.eq.0) then 
        write(6,'(1i8,20(1x,g18.8))') in,test_r(in),test_z(in),btot,te,ti,ne,ef,vb,br,bz,bt
     else
        write(error_message_data,'(a,3i8,20(1x,g18.8))') 'ERROR GETTING PLASMA:',ierr,in,test_r(in),test_z(in)
        call errmsg('PLASMA TEST',error_message_data)
     endif

  end do


  ! Define sample space - call and report plasma conditions
  ! Values chosen here are arbitrary for testing purposes

  nr = 100
  nz = 100

  rmin = 1.4
  rmax = 1.6

  zmin = -1.37
  zmax = -1.17

  dr = (rmax-rmin)/nr
  dz = (zmax-zmin)/nz

  write(6,*) 
  write(6,'(a,i10)') 'Ordered by R: Interpolation =',interpolate_option
  write(6,*) 

  do ir = 0,nr

     rt = rmin + ir * dr

     do iz = 0,nz

        zt = zmin + iz * dz

        call get_oedge_plasma(rt,zt,ne,te,ti,vb,ef,btot,br,bz,bt,ierr)

        if (ierr.eq.0) then 
           write(6,'(2i8,20(1x,g18.8))') ir,iz,rt,zt,ne,te,ti,vb,ef,btot,br,bz,bt
        else
           write(error_message_data,'(a,3i8,20(1x,g18.8))') 'ERROR GETTING PLASMA:',ierr,ir,iz,rt,zt,ne,te,ti,vb,ef,btot,br,bz,bt
           call errmsg('PLASMA TEST',error_message_data)
        endif

     end do
     write(6,*) 

  end do

  write(6,*) 
  write(6,'(a,i10)') 'Ordered by Z',interpolate_option
  write(6,*) 


  do iz = 0,nz

     zt = zmin + iz * dz

     do ir = 0,nr

        rt = rmin + ir * dr

        call get_oedge_plasma(rt,zt,ne,te,ti,vb,ef,btot,br,bz,bt,ierr)

        if (ierr.eq.0) then 
           write(6,'(2i8,20(1x,g18.8))') ir,iz,rt,zt,ne,te,ti,vb,ef,btot,br,bz,bt
        else
           write(error_message_data,'(a,3i8,20(1x,g18.8))') 'ERROR GETTING PLASMA:',ierr,ir,iz,rt,zt,ne,te,ti,vb,ef,btot,br,bz,bt
           call errmsg('PLASMA TEST',error_message_data)
        endif

     end do
     write(6,*) 

  end do




  ! Clean up when done calling OEDGE plasma interface

  call close_oedge_plasma



end program plasma_test
