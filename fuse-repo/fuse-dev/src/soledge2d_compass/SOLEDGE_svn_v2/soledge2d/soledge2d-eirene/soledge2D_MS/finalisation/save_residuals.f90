subroutine save_residuals(n_ite)
  use all_variables, only : global_variables, global_parameters
  implicit none
  integer*4,intent(in) :: n_ite
  integer*4 :: n
  if(modulo(dble(n_ite),dble(global_parameters%N_iterations)/1000.D0).eq.0) then
     open(unit=100,file='residuals',status='unknown',access='append')
     write(100,101) (global_variables%residuals(n)%resn,n=0,global_parameters%N_ions),&
         (global_variables%residuals(n)%resG,n=0,global_parameters%N_ions),&
         (global_variables%residuals(n)%resT,n=0,global_parameters%N_ions)
101 format(512es15.7) 
     close(100)
  end if
end subroutine save_residuals


subroutine save_residuals_force()
  use all_variables, only : global_variables, global_parameters
  implicit none
  integer*4 :: n
  open(unit=100,file='residuals',status='unknown',access='append')
  write(100,101) (global_variables%residuals(n)%resn,n=0,global_parameters%N_ions),&
       (global_variables%residuals(n)%resG,n=0,global_parameters%N_ions),&
       (global_variables%residuals(n)%resT,n=0,global_parameters%N_ions)
101 format(512es15.7) 
  close(100)
end subroutine save_residuals_force
