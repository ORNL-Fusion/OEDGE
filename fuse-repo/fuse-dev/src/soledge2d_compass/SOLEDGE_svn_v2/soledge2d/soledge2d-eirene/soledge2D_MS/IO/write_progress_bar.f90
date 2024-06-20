subroutine write_progress_bar(n_ite,nchar,t_ite)
  use all_variables, only : global_parameters,global_variables, reference_parameters
  implicit none
  integer*4,intent(in) :: n_ite
  character(nchar) buffer
  character(8) buffer2
  character(13) buffer3 
  character(20) buffer4 
  character(20) format_empty
  integer*4,intent(in) :: nchar
  integer*4 :: frac
  integer*4 :: k
  double precision,intent(in) :: t_ite
  if(global_parameters%N_iterations>=1000) then
     if(modulo(dble(n_ite),dble(global_parameters%N_iterations)/1000.d0).eq.0) then    
        write(format_empty,40) nchar
        write(buffer,trim(format_empty))
        frac=floor(dble(n_ite)/dble(global_parameters%N_iterations)*dble(nchar-2))
        write(buffer(1:1),'(A1)') "["
        write(buffer(nchar:nchar),'(A1)') "]"
        do k=2,frac+1
           write(buffer(k:k),'(A1)') "#"
        end do
        write(buffer2(1:1),'(1X)')
        write(buffer2(2:6),'(F5.1)') dble(n_ite)/dble(global_parameters%N_iterations)*100.
        write(buffer2(7:7),'(A1)') "%"
        write(buffer2(8:8),'(1X)')
        write(buffer3(1:4),'(A4)') " dt="
        write(buffer3(5:13),'(ES9.2)') global_variables%dt*reference_parameters%fields%tau0
        write(buffer4(1:9),'(A9)') " (cptime="
        write(buffer4(10:18),'(ES9.2)') t_ite
        write(buffer4(19:20),'(A1)') ")"
        write(format_empty,41) nchar+8+13+20
        write(6,fmt=trim(format_empty)) '+',char(13),buffer//buffer2//buffer3//buffer4
     end if
  else
     write(format_empty,40) nchar
     write(buffer,trim(format_empty))
     frac=floor(dble(n_ite)/dble(global_parameters%N_iterations)*dble(nchar-2))
     write(buffer(1:1),'(A1)') "["
     write(buffer(nchar:nchar),'(A1)') "]"
     do k=2,frac+1
        write(buffer(k:k),'(A1)') "#"
     end do
     write(buffer2(1:1),'(1X)')
     write(buffer2(2:6),'(F5.1)') dble(n_ite)/dble(global_parameters%N_iterations)*100.
     write(buffer2(7:7),'(A1)') "%"
     write(buffer2(8:8),'(1X)')
     write(buffer3(1:4),'(A4)') " dt="
     write(buffer3(5:13),'(ES9.2)') global_variables%dt*reference_parameters%fields%tau0
     write(buffer4(1:9),'(A9)') " (cptime="
     write(buffer4(10:18),'(ES9.2)') t_ite
     write(buffer4(19:20),'(A1)') ")"
     write(format_empty,41) nchar+8+13+20
     write(6,fmt=trim(format_empty)) '+',char(13),buffer//buffer2//buffer3//buffer4
  end if
40 format('(',I2,'X)')
41 format('(a1,a1,a',I2,')')
end subroutine write_progress_bar
