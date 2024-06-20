subroutine check_fields(zone)
  use all_variables, only : global_parameters, flags
  use Mlog_message
  use MZone
  implicit none
  type(Tzone),intent(in) :: zone
  integer*4 :: n,Nx,Nz
  integer*4 :: i,j
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do n=1,global_parameters%N_ions
     do i=1,Nx
        do j=1,Nz
           if(isNaN(zone%species(n)%var(2)%density(i,j))) then
              call write_log_message(msg_simulation_crashed)
              write(*,*) 'N NaN ',zone%number,i,j,n
              write(*,*) zone%species(n)%var(1)%density(i,j)
              stop
           end if
           if(isNaN(zone%species(n)%var(2)%Gamma(i,j))) then
              call write_log_message(msg_simulation_crashed)
              write(*,*) 'G NaN ',zone%number,i,j,n
              write(*,*) zone%species(n)%var(1)%Gamma(i,j)
              stop
           end if
           if(isNaN(zone%species(n)%var(2)%temperature(i,j))) then
              call write_log_message(msg_simulation_crashed)
              write(*,*) 'T NaN ',zone%number,i,j,n
              write(*,*) zone%species(n)%var(1)%temperature(i,j)
              stop
           end if
           if(flags%turbulence_model.eq.1) then
              if(isNaN(zone%kepsilon(2)%k(i,j))) then
                 call write_log_message(msg_simulation_crashed)
                 write(*,*) 'k NaN ',zone%number,i,j
                 write(*,*) zone%kepsilon(1)%k(i,j)
                 stop
              end if
           end if
        end do
     end do
  end do
end subroutine check_fields
