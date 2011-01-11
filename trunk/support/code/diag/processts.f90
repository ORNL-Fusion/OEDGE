program processts
use analyse_ts 

implicit none


character*100 :: filename
integer:: nlines,in

filename = 'divts_'

call setup_bins

do in = 134580,134597

   write(filename,'(a,i6,a)') 'divts_',in,'.dat'
   write(0,*) 'Filename:',trim(filename)

   call read_ts(filename,nlines)

   call accumulate_data(nlines)

end do


call analyse_print_data






end program processts
