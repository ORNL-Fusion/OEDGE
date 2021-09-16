program processfs
  use analyse_fs

  implicit none

  character*100 :: filename,filebase,filetype,fschan,filenameout
  integer:: nlines,in
  real :: tstart, tend

  integer ::  shot 

  real :: minval,maxval

  integer :: nbins,printopt
  integer :: ichan,itype


  character*2 :: filetypes(4)
  character*6 :: channel_list(3,4)

  channel_list(1,1) = 'fs03'
  channel_list(1,2) = 'fs03db'
  channel_list(1,3) = 'fs03c2'
  channel_list(1,4) = 'fs03c3'

  channel_list(2,1) = 'fs04'
  channel_list(2,2) = 'none'    ! fs04db
  channel_list(2,3) = 'none'    ! fs04c2
  channel_list(2,4) = 'fs04c3'

  channel_list(3,1) = 'fs05'
  channel_list(3,2) = 'fs05db'
  channel_list(3,3) = 'fs05c2'
  channel_list(3,4) = 'fs05c3'


  filetypes(1) = 'da'
  filetypes(2) = 'db'
  filetypes(3) = 'c2'
  filetypes(4) = 'c3'

  !
  ! Print opt = 0 ... gives baseline average evolution over time
  !
  ! Print opt = 1 ... gives all print data
  !
  !

  printopt = 0


  !
  ! Note: channel names for db and c3 include the type but da data just uses the base channel name
  !

  !filetype = 'da'
  !fschan = 'fs03'

  !minval = 8e17
  !maxval = 2e19

  nbins  = 200


  do ichan = 1,3

     do itype = 1,3

        filetype = filetypes(itype)
        fschan = channel_list(ichan,itype)

        ! exclude non-existant channels
        if (fschan.eq.'none') cycle

        
        write(0,'(a,2i4,2x,a,2x,a)') 'Analysing:',ichan,itype,trim(filetype),trim(fschan)

        if (filetype.eq.'da') then 
           minval = 2e19
           maxval = 6e20
        elseif (filetype.eq.'db') then 
           minval = 5e17
           maxval = 6e18
        elseif (filetype.eq.'c3') then 
           minval = 8e17
           maxval = 2e19
        endif


        do shot = 134580,134597

           write(filename,'(a,i6,a)') 'fs_'//trim(filetype)//'_data_',shot,'_1500.00_4500.00.dat'
           write(0,'(a,a)') '1 Filename:'//trim(filename),' Channel: '//trim(fschan)

           call setup_fs(minval,maxval,nbins,fschan)

           !write(0,'(a,a)') '2 Filename:'//trim(filename),' Channel: '//trim(fschan)

           call read_fs(filename,fschan,nlines)

           !write(0,'(a,a)') '3 Filename:'//trim(filename),' Channel: '//trim(fschan)

           call accumulate_data(nlines)

           write(filenameout,'(a,i6,a)') trim(fschan)//'_',shot,'.dat'

           !write(0,'(a,a)') '4 Filename:'//trim(filename),' Channel: '//trim(fschan)

           call analyse_print_data(filenameout,printopt)

           !write(0,'(a,a)') '5 Filename:'//trim(filename),' Channel: '//trim(fschan)

        end do

     end do

  end do


  ! Accumulate data over all shots


  do ichan = 1,3

     do itype = 1,3

        filetype = filetypes(itype)
        fschan = channel_list(ichan,itype)

        if (filetype.eq.'da') then 
           minval = 2e19
           maxval = 6e20
        elseif (filetype.eq.'db') then 
           minval = 5e17
           maxval = 6e18
        elseif (filetype.eq.'c3') then 
           minval = 8e17
           maxval = 2e19
        endif


        call setup_fs(minval,maxval,nbins,fschan)


        do shot = 134580,134597

           write(filename,'(a,i6,a)') 'fs_'//trim(filetype)//'_data_',shot,'_1500.00_4500.00.dat'
           write(0,'(a)') 'Filename:'//trim(filename)

           call read_fs(filename,fschan,nlines)

           call accumulate_data(nlines)


        end do

        write(filenameout,'(a)') trim(fschan)//'_summary.dat'
        call analyse_print_data(filenameout,printopt)



     end do

  end do






end program processfs
