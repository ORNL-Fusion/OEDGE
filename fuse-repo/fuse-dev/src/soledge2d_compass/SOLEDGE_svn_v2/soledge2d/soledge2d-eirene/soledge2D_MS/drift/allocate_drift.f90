subroutine allocate_drift(zone)
use all_variables, only : global_parameters
  use MZone 
  implicit none
  type(TZone),intent(inout) :: zone
  integer*4 :: Nx,Nz
  integer*4 :: n
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do n=0,global_parameters%N_ions
     allocate(zone%species(n)%drifts%uEp(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%drifts%uEt(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%drifts%udp(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%drifts%udt(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%drifts%uBp(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%drifts%uBt(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%drifts%jdiam(1:Nx,1:Nz,1:4))
     allocate(zone%species(n)%drifts%jExB(1:Nx,1:Nz,1:4))
     allocate(zone%species(n)%drifts%jBxDB(1:Nx,1:Nz,1:4))
     allocate(zone%species(n)%drifts%uEps(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%drifts%uEts(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%drifts%uBps(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%drifts%uBts(0:Nx+1,0:Nz+1))
     zone%species(n)%drifts%uEp = 0.d0
     zone%species(n)%drifts%uEt = 0.d0
     zone%species(n)%drifts%udp = 0.d0
     zone%species(n)%drifts%udt = 0.d0
     zone%species(n)%drifts%uBp = 0.d0
     zone%species(n)%drifts%uBt = 0.d0
  end do
end subroutine allocate_drift
