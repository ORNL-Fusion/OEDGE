      program gfile_vessel
      implicit none

      integer i,j,nw,nh,nlim,nbbs
      real gr1
      real, dimension(1:1000):: xlim,ylim
      character(len=10) gc1,gc2,gc3,gc4,gc5,gc6

      character*256 :: input_file = 'gfile'
      character*256 :: output_file = 'vvfile.ogr'

      integer :: nargs,ierr

!
! Read and assign command line arguments
!
! Command line is:  'file name'   tmin    tmax   chi_lim -e 'elm file name'
!
      nargs = iargc()

      ! error - too many arguments
      if (nargs.gt.2) then 
         call printhelp(input_file,output_file)
         return
      elseif (nargs.eq.2) then
          ! input and output files specified
          call getarg(1,input_file)
          call getarg(2,output_file)
      elseif (nargs.eq.1) then
          !input file specified
          call getarg(1,input_file)
          output_file=trim(input_file)//'.out'
      endif


      
      open(unit=7,file=trim(input_file),status='old',iostat=ierr)
      if (ierr.ne.0) then 
         call printhelp(input_file,output_file)
         return
      endif

      open(unit=8,file=trim(output_file),status='unknown')

      read(7,2000) gc1,gc2,gc3,gc4,gc5,gc6,i,nw,nh
      read(7,2020) gr1
      read(7,2020) gr1
      read(7,2020) gr1
      read(7,2020) gr1
      read(7,2020) (gr1,i=1,nw)
      read(7,2020) (gr1,i=1,nw)
      read(7,2020) (gr1,i=1,nw)
      read(7,2020) (gr1,i=1,nw)
      read(7,2020) ((gr1,i=1,nw),j=1,nh)
      read(7,2020) (gr1,i=1,nw)
      read(7,2022) nbbs,nlim
      read(7,2020) (gr1,gr1,i=1,nbbs)
      read(7,2020) (xlim(i),ylim(i),i=1,nlim)
      do i=1,nlim
         write(8,*) xlim(i)*1000.,ylim(i)*1000.
      end do

 2000 format (6a8,3i4)
 2020 format (5e16.9)
 2022 format (2i5)
      end

      subroutine printhelp(input_file,output_file)
      implicit none
      character*(*) input_file,output_file
      
         write(0,'(a)') 'This routine extracts the wall definition'//
     >                  ' from a G EFIT file'
         write(0,'(a)') '<command> (<input file>) (<output file>)'
         write(0,'(a)') 'If output file is not specified'//
     >        ' <input file>.out is used'
         write(0,'(a)') 'If neither are specified then the following'//
     >        ' defaults are used'
         write(0,'(a,a)') '<input file>  = ',trim(input_file)
         write(0,'(a,a)') '<output file> = ',trim(output_file)

      return
      end
