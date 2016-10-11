program read_irtv
use saf_files
implicit none

!
! This code reads in an IRTV video file and performs a number of actions depending on user input. 
! The UI is very primitive :)
!

character*256 filename


! Input IRTV file name

write(6,*) 'Please input the IRTV file name:'
read(5,*) filename

