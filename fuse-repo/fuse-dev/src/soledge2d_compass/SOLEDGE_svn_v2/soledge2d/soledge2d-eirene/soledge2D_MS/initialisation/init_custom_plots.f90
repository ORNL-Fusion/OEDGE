subroutine init_custom_plots()
  use Mplots
  use Mreaders
  use all_variables, only : custom_plots, N_custom_plots
  implicit none
  integer :: finput
  logical :: dir_e
  integer*4 :: k,n
  character(128) :: String
  integer*4,allocatable :: bufferI(:)
  inquire(File='custom_save.txt',exist=dir_e)
  if(dir_e) then
     finput=100
     open(unit=finput,file='custom_save.txt',status='unknown')
     read(finput,1) N_custom_plots
     allocate(custom_plots(1:N_custom_plots))
     do k=1,N_custom_plots
        call skip_line(finput,2)
        read(finput,2) custom_plots(k)%type
        read(finput,3) custom_plots(k)%nzones
        allocate(custom_plots(k)%zones(1:custom_plots(k)%nzones))
        allocate(bufferI(1:custom_plots(k)%nzones))
        read(finput,4) String
        call parse_line_integer(String,custom_plots(k)%nzones,bufferI)
        do n=1,custom_plots(k)%nzones
           custom_plots(k)%zones(n)=bufferI(n)
        end do
        read(finput,5) custom_plots(k)%coord
        if(custom_plots(k)%type.eq.3) then
           read(finput,6) custom_plots(k)%coord2
        end if
        deallocate(bufferI)
     end do
     close(finput)
  end if
1 format(8X,I2)
2 format(7X,I1)
3 format(9X,I2)
4 format(A128)  
5 format(8X,I3)
6 format(9X,I3)
end subroutine init_custom_plots
