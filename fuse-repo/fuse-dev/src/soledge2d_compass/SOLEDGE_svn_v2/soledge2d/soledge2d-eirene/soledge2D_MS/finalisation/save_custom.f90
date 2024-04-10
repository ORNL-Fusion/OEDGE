subroutine save_custom(n_ite)
  use all_variables, only : global_variables, global_parameters, zones, reference_parameters
  implicit none
  integer*4,intent(in) :: n_ite
  integer*4 :: n,i,j,k
!!$  if(modulo(dble(n_ite),dble(global_parameters%N_iterations)/1000.D0).eq.0) then
!!$     open(unit=100,file='custom',status='unknown')
!!$     do i=1,zones(2)%mesh%Nx
!!$        write(100,101) zones(2)%species(0)%var(1)%temperature(i,1), zones(2)%electric_fields(1)%phi(i,1), zones(2)%electric_fields(1)%vorticity(i,1)
!!$     end do
!!$     do i=1,zones(4)%mesh%Nx
!!$        write(100,101) zones(4)%species(0)%var(1)%temperature(i,1), zones(4)%electric_fields(1)%phi(i,1), zones(4)%electric_fields(1)%vorticity(i,1)
!!$     end do
!!$101  format(512es15.7) 
!!$     close(100)
!!$  end if

!!$  if(modulo(dble(n_ite),dble(global_parameters%N_iterations)/1000.D0).eq.0) then
!!$     open(unit=100,file='custom',status='unknown')
!!$     do k=1,4
!!$        do i=1,zones(k)%mesh%Nx
!!$           write(100,101) zones(k)%species(0)%var(1)%temperature(i,60)*(1.D0-zones(k)%masks%chi2(i,60)),&
!!$                zones(k)%electric_fields(1)%phi(i,60)*(1.D0-zones(k)%masks%chi2(i,60)),&
!!$                zones(k)%electric_fields(1)%vorticity(i,60)*(1.D0-zones(k)%masks%chi2(i,60)),&
!!$                zones(k)%kepsilon(1)%mu_t(i,60)*(1.D0-zones(k)%masks%chi2(i,60))*reference_parameters%fields%k0*reference_parameters%fields%tau0
!!$        end do
!!$     end do
!!$101  format(512es15.7) 
!!$     close(100)
!!$  end if
!!$  if(modulo(dble(n_ite),dble(global_parameters%N_iterations)/1000.D0).eq.0) then
!!$     open(unit=100,file='custom2',status='unknown')
!!$     do k=1,4
!!$        do i=1,zones(k)%mesh%Nx
!!$           write(100,101) zones(k)%species(0)%var(1)%temperature(i,1)*(1.D0-zones(k)%masks%chi2(i,1)),&
!!$                zones(k)%electric_fields(1)%phi(i,1)*(1.D0-zones(k)%masks%chi2(i,1)),&
!!$                zones(k)%electric_fields(1)%vorticity(i,1)*(1.D0-zones(k)%masks%chi2(i,1)),&
!!$                zones(k)%kepsilon(1)%mu_t(i,1)*(1.D0-zones(k)%masks%chi2(i,1))*reference_parameters%fields%k0*reference_parameters%fields%tau0
!!$        end do
!!$     end do
!!$     close(100)
!!$  end if


!!$  if(modulo(dble(n_ite),dble(global_parameters%N_iterations)/1000.D0).eq.0) then
!!$     open(unit=100,file='custom',status='unknown')
!!$     k=6
!!$     do i=1,zones(k)%mesh%Nx
!!$        write(100,101) zones(k)%species(0)%var(1)%temperature(i,40)*(1.D0-zones(k)%masks%chi2(i,40)),&
!!$             zones(k)%electric_fields(1)%phi(i,40)*(1.D0-zones(k)%masks%chi2(i,40)),&
!!$             zones(k)%electric_fields(1)%vorticity(i,40)*(1.D0-zones(k)%masks%chi2(i,40)),&
!!$             0.09*zones(k)%kepsilon(1)%k(i,40)**2/zones(k)%kepsilon(1)%epsilon(i,40)&
!!$             *(1.D0-zones(k)%masks%chi2(i,40))*reference_parameters%fields%k0*reference_parameters%fields%tau0
!!$     end do
!!$     k=21
!!$     do i=1,zones(k)%mesh%Nx
!!$        write(100,101) zones(k)%species(0)%var(1)%temperature(i,40)*(1.D0-zones(k)%masks%chi2(i,40)),&
!!$             zones(k)%electric_fields(1)%phi(i,40)*(1.D0-zones(k)%masks%chi2(i,40)),&
!!$             zones(k)%electric_fields(1)%vorticity(i,40)*(1.D0-zones(k)%masks%chi2(i,40)),&
!!$             0.09*zones(k)%kepsilon(1)%k(i,40)**2/zones(k)%kepsilon(1)%epsilon(i,40)&
!!$             *(1.D0-zones(k)%masks%chi2(i,40))*reference_parameters%fields%k0*reference_parameters%fields%tau0
!!$     end do
!!$     k=5
!!$     do i=1,zones(k)%mesh%Nx
!!$        write(100,101) zones(k)%species(0)%var(1)%temperature(i,40)*(1.D0-zones(k)%masks%chi2(i,40)),&
!!$             zones(k)%electric_fields(1)%phi(i,40)*(1.D0-zones(k)%masks%chi2(i,40)),&
!!$             zones(k)%electric_fields(1)%vorticity(i,40)*(1.D0-zones(k)%masks%chi2(i,40)),&
!!$             0.09*zones(k)%kepsilon(1)%k(i,40)**2/zones(k)%kepsilon(1)%epsilon(i,40)&
!!$             *(1.D0-zones(k)%masks%chi2(i,40))*reference_parameters%fields%k0*reference_parameters%fields%tau0
!!$     end do
!!$     k=4
!!$     do i=1,zones(k)%mesh%Nx
!!$        write(100,101) zones(k)%species(0)%var(1)%temperature(i,40)*(1.D0-zones(k)%masks%chi2(i,40)),&
!!$             zones(k)%electric_fields(1)%phi(i,40)*(1.D0-zones(k)%masks%chi2(i,40)),&
!!$             zones(k)%electric_fields(1)%vorticity(i,40)*(1.D0-zones(k)%masks%chi2(i,40)),&
!!$             0.09*zones(k)%kepsilon(1)%k(i,40)**2/zones(k)%kepsilon(1)%epsilon(i,40)&
!!$             *(1.D0-zones(k)%masks%chi2(i,40))*reference_parameters%fields%k0*reference_parameters%fields%tau0
!!$     end do
!!$101  format(512es15.7) 
!!$     close(100)
!!$  end if
!!$
!!$    if(modulo(dble(n_ite),dble(global_parameters%N_iterations)/1000.D0).eq.0) then
!!$     open(unit=100,file='custom2',status='unknown')
!!$     k=13
!!$     do i=1,zones(k)%mesh%Nx
!!$        write(100,101) zones(k)%species(0)%var(1)%temperature(i,40)*(1.D0-zones(k)%masks%chi2(i,40)),&
!!$             zones(k)%electric_fields(1)%phi(i,40)*(1.D0-zones(k)%masks%chi2(i,40)),&
!!$             zones(k)%electric_fields(1)%vorticity(i,40)*(1.D0-zones(k)%masks%chi2(i,40)),&
!!$             0.09*zones(k)%kepsilon(1)%k(i,40)**2/zones(k)%kepsilon(1)%epsilon(i,40)&
!!$             *(1.D0-zones(k)%masks%chi2(i,40))*reference_parameters%fields%k0*reference_parameters%fields%tau0
!!$     end do
!!$     k=24
!!$     do i=1,zones(k)%mesh%Nx
!!$        write(100,101) zones(k)%species(0)%var(1)%temperature(i,40)*(1.D0-zones(k)%masks%chi2(i,40)),&
!!$             zones(k)%electric_fields(1)%phi(i,40)*(1.D0-zones(k)%masks%chi2(i,40)),&
!!$             zones(k)%electric_fields(1)%vorticity(i,40)*(1.D0-zones(k)%masks%chi2(i,40)),&
!!$             0.09*zones(k)%kepsilon(1)%k(i,40)**2/zones(k)%kepsilon(1)%epsilon(i,40)&
!!$             *(1.D0-zones(k)%masks%chi2(i,40))*reference_parameters%fields%k0*reference_parameters%fields%tau0
!!$     end do
!!$     k=12
!!$     do i=1,zones(k)%mesh%Nx
!!$        write(100,101) zones(k)%species(0)%var(1)%temperature(i,40)*(1.D0-zones(k)%masks%chi2(i,40)),&
!!$             zones(k)%electric_fields(1)%phi(i,40)*(1.D0-zones(k)%masks%chi2(i,40)),&
!!$             zones(k)%electric_fields(1)%vorticity(i,40)*(1.D0-zones(k)%masks%chi2(i,40)),&
!!$             0.09*zones(k)%kepsilon(1)%k(i,40)**2/zones(k)%kepsilon(1)%epsilon(i,40)&
!!$             *(1.D0-zones(k)%masks%chi2(i,40))*reference_parameters%fields%k0*reference_parameters%fields%tau0
!!$     end do
!!$     k=11
!!$     do i=1,zones(k)%mesh%Nx
!!$        write(100,101) zones(k)%species(0)%var(1)%temperature(i,40)*(1.D0-zones(k)%masks%chi2(i,40)),&
!!$             zones(k)%electric_fields(1)%phi(i,40)*(1.D0-zones(k)%masks%chi2(i,40)),&
!!$             zones(k)%electric_fields(1)%vorticity(i,40)*(1.D0-zones(k)%masks%chi2(i,40)),&
!!$             0.09*zones(k)%kepsilon(1)%k(i,40)**2/zones(k)%kepsilon(1)%epsilon(i,40)&
!!$             *(1.D0-zones(k)%masks%chi2(i,40))*reference_parameters%fields%k0*reference_parameters%fields%tau0
!!$     end do
!!$     close(100)
!!$  end if

  !TS
  if(modulo(dble(n_ite),dble(global_parameters%N_iterations)/1000.D0).eq.0) then
     open(unit=100,file='custom',status='unknown')
     k=9
     j=28
     do i=0,zones(k)%mesh%Nx
        write(100,101) zones(k)%species(0)%var(1)%density(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
             zones(k)%electric_fields(1)%phi(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
             zones(k)%electric_fields(1)%vorticity(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
!             zones(k)%kepsilon(1)%mu_t(i,j)*(reference_parameters%geometry%rs0**2/reference_parameters%fields%tau0)&
!             *(1.D0-zones(k)%masks%chi2(i,j)),& 
 !            zones(k)%kepsilon(1)%k(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
             zones(k)%species(0)%var(1)%temperature(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
 !            zones(k)%kepsilon(1)%interchange(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
 !            zones(k)%kepsilon(1)%UE_shear(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
 !            zones(k)%species(1)%transport_perp%v_pinch(i,j),&
 !            zones(k)%kepsilon(1)%epsilon(i,j)*(1.D0-zones(k)%masks%chi2(i,j))
             zones(k)%species(1)%drifts%uEt(i,j)*(1.D0-zones(k)%masks%chi2(i,j))*zones(k)%metric_coefficients%c_tt(i,j)
     end do
     k=12
     do i=1,zones(k)%mesh%Nx
        write(100,101) zones(k)%species(0)%var(1)%density(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
             zones(k)%electric_fields(1)%phi(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
             zones(k)%electric_fields(1)%vorticity(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
!             zones(k)%kepsilon(1)%mu_t(i,j)*(reference_parameters%geometry%rs0**2/reference_parameters%fields%tau0)&
!             *(1.D0-zones(k)%masks%chi2(i,j)),& 
!             zones(k)%kepsilon(1)%k(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
             zones(k)%species(0)%var(1)%temperature(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
!             zones(k)%kepsilon(1)%interchange(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
!             zones(k)%kepsilon(1)%UE_shear(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
!             zones(k)%species(1)%transport_perp%v_pinch(i,j),&
!             zones(k)%kepsilon(1)%epsilon(i,j)*(1.D0-zones(k)%masks%chi2(i,j))
             zones(k)%species(1)%drifts%uEt(i,j)*(1.D0-zones(k)%masks%chi2(i,j))*zones(k)%metric_coefficients%c_tt(i,j)
     end do
     k=15
     do i=1,zones(k)%mesh%Nx
        write(100,101) zones(k)%species(0)%var(1)%density(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
             zones(k)%electric_fields(1)%phi(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
             zones(k)%electric_fields(1)%vorticity(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
!             zones(k)%kepsilon(1)%mu_t(i,j)*(reference_parameters%geometry%rs0**2/reference_parameters%fields%tau0)&
!             *(1.D0-zones(k)%masks%chi2(i,j)),& 
!             zones(k)%kepsilon(1)%k(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
             zones(k)%species(0)%var(1)%temperature(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
!             zones(k)%kepsilon(1)%interchange(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
!             zones(k)%kepsilon(1)%UE_shear(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
!             zones(k)%species(1)%transport_perp%v_pinch(i,j),&
!             zones(k)%kepsilon(1)%epsilon(i,j)*(1.D0-zones(k)%masks%chi2(i,j))
             zones(k)%species(1)%drifts%uEt(i,j)*(1.D0-zones(k)%masks%chi2(i,j))*zones(k)%metric_coefficients%c_tt(i,j)
     end do
     k=18
     do i=1,zones(k)%mesh%Nx
        write(100,101) zones(k)%species(0)%var(1)%density(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
             zones(k)%electric_fields(1)%phi(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
             zones(k)%electric_fields(1)%vorticity(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
 !            zones(k)%kepsilon(1)%mu_t(i,j)*(reference_parameters%geometry%rs0**2/reference_parameters%fields%tau0)&
 !            *(1.D0-zones(k)%masks%chi2(i,j)),& 
 !            zones(k)%kepsilon(1)%k(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
             zones(k)%species(0)%var(1)%temperature(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
!             zones(k)%kepsilon(1)%interchange(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
!             zones(k)%kepsilon(1)%UE_shear(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
!             zones(k)%species(1)%transport_perp%v_pinch(i,j),&
!             zones(k)%kepsilon(1)%epsilon(i,j)*(1.D0-zones(k)%masks%chi2(i,j))
             zones(k)%species(1)%drifts%uEt(i,j)*(1.D0-zones(k)%masks%chi2(i,j))*zones(k)%metric_coefficients%c_tt(i,j)
     end do
     k=21
     do i=1,zones(k)%mesh%Nx
        write(100,101) zones(k)%species(0)%var(1)%density(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
             zones(k)%electric_fields(1)%phi(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
             zones(k)%electric_fields(1)%vorticity(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
!             zones(k)%kepsilon(1)%mu_t(i,j)*(reference_parameters%geometry%rs0**2/reference_parameters%fields%tau0)&
!             *(1.D0-zones(k)%masks%chi2(i,j)),& 
!             zones(k)%kepsilon(1)%k(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
             zones(k)%species(0)%var(1)%temperature(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
!             zones(k)%kepsilon(1)%interchange(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
!             zones(k)%kepsilon(1)%UE_shear(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
!             zones(k)%species(1)%transport_perp%v_pinch(i,j),&
!             zones(k)%kepsilon(1)%epsilon(i,j)*(1.D0-zones(k)%masks%chi2(i,j))
             zones(k)%species(1)%drifts%uEt(i,j)*(1.D0-zones(k)%masks%chi2(i,j))*zones(k)%metric_coefficients%c_tt(i,j)
     end do
     k=24
     do i=1,zones(k)%mesh%Nx
        write(100,101) zones(k)%species(0)%var(1)%density(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
             zones(k)%electric_fields(1)%phi(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
             zones(k)%electric_fields(1)%vorticity(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
!             zones(k)%kepsilon(1)%mu_t(i,j)*(reference_parameters%geometry%rs0**2/reference_parameters%fields%tau0)&
!             *(1.D0-zones(k)%masks%chi2(i,j)),& 
!             zones(k)%kepsilon(1)%k(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
             zones(k)%species(0)%var(1)%temperature(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
!             zones(k)%kepsilon(1)%interchange(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
!             zones(k)%kepsilon(1)%UE_shear(i,j)*(1.D0-zones(k)%masks%chi2(i,j)),&
!             zones(k)%species(1)%transport_perp%v_pinch(i,j),&
!             zones(k)%kepsilon(1)%epsilon(i,j)*(1.D0-zones(k)%masks%chi2(i,j))
             zones(k)%species(1)%drifts%uEt(i,j)*(1.D0-zones(k)%masks%chi2(i,j))*zones(k)%metric_coefficients%c_tt(i,j)
     end do
101  format(512es15.7) 
     close(100)
  end if
!!$
!!$  if(modulo(dble(n_ite),dble(global_parameters%N_iterations)/1000.D0).eq.0) then
!!$     open(unit=100,file='custom2',status='unknown')
!!$     k=10
!!$     do i=1,zones(k)%mesh%Nx
!!$        write(100,101) zones(k)%species(0)%var(1)%density(i,50)*(1.D0-zones(k)%masks%chi2(i,50)),&
!!$             zones(k)%electric_fields(1)%phi(i,50)*(1.D0-zones(k)%masks%chi2(i,50)),&
!!$             zones(k)%electric_fields(1)%vorticity(i,50)*(1.D0-zones(k)%masks%chi2(i,50)),&
!!$             0.09*zones(k)%kepsilon(1)%k(i,50)**2/zones(k)%kepsilon(1)%epsilon(i,50)&
!!$             *(1.D0-zones(k)%masks%chi2(i,50))*reference_parameters%fields%k0*reference_parameters%fields%tau0,& 
!!$             zones(k)%species(0)%var(1)%temperature(i,50)*(1.D0-zones(k)%masks%chi2(i,50))
!!$     end do
!!$     k=13
!!$     do i=1,zones(k)%mesh%Nx
!!$        write(100,101) zones(k)%species(0)%var(1)%density(i,50)*(1.D0-zones(k)%masks%chi2(i,50)),&
!!$             zones(k)%electric_fields(1)%phi(i,50)*(1.D0-zones(k)%masks%chi2(i,50)),&
!!$             zones(k)%electric_fields(1)%vorticity(i,50)*(1.D0-zones(k)%masks%chi2(i,50)),&
!!$             0.09*zones(k)%kepsilon(1)%k(i,50)**2/zones(k)%kepsilon(1)%epsilon(i,50)&
!!$             *(1.D0-zones(k)%masks%chi2(i,50))*reference_parameters%fields%k0*reference_parameters%fields%tau0,& 
!!$             zones(k)%species(0)%var(1)%temperature(i,50)*(1.D0-zones(k)%masks%chi2(i,50))
!!$     end do
!!$     k=16
!!$     do i=1,zones(k)%mesh%Nx
!!$        write(100,101) zones(k)%species(0)%var(1)%density(i,50)*(1.D0-zones(k)%masks%chi2(i,50)),&
!!$             zones(k)%electric_fields(1)%phi(i,50)*(1.D0-zones(k)%masks%chi2(i,50)),&
!!$             zones(k)%electric_fields(1)%vorticity(i,50)*(1.D0-zones(k)%masks%chi2(i,50)),&
!!$             0.09*zones(k)%kepsilon(1)%k(i,50)**2/zones(k)%kepsilon(1)%epsilon(i,50)&
!!$             *(1.D0-zones(k)%masks%chi2(i,50))*reference_parameters%fields%k0*reference_parameters%fields%tau0,& 
!!$             zones(k)%species(0)%var(1)%temperature(i,50)*(1.D0-zones(k)%masks%chi2(i,50))
!!$     end do
!!$     k=19
!!$     do i=1,zones(k)%mesh%Nx
!!$        write(100,101) zones(k)%species(0)%var(1)%density(i,50)*(1.D0-zones(k)%masks%chi2(i,50)),&
!!$             zones(k)%electric_fields(1)%phi(i,50)*(1.D0-zones(k)%masks%chi2(i,50)),&
!!$             zones(k)%electric_fields(1)%vorticity(i,50)*(1.D0-zones(k)%masks%chi2(i,50)),&
!!$             0.09*zones(k)%kepsilon(1)%k(i,50)**2/zones(k)%kepsilon(1)%epsilon(i,50)&
!!$             *(1.D0-zones(k)%masks%chi2(i,50))*reference_parameters%fields%k0*reference_parameters%fields%tau0,& 
!!$             zones(k)%species(0)%var(1)%temperature(i,50)*(1.D0-zones(k)%masks%chi2(i,50))
!!$     end do
!!$     k=22
!!$     do i=1,zones(k)%mesh%Nx
!!$        write(100,101) zones(k)%species(0)%var(1)%density(i,50)*(1.D0-zones(k)%masks%chi2(i,50)),&
!!$             zones(k)%electric_fields(1)%phi(i,50)*(1.D0-zones(k)%masks%chi2(i,50)),&
!!$             zones(k)%electric_fields(1)%vorticity(i,50)*(1.D0-zones(k)%masks%chi2(i,50)),&
!!$             0.09*zones(k)%kepsilon(1)%k(i,50)**2/zones(k)%kepsilon(1)%epsilon(i,50)&
!!$             *(1.D0-zones(k)%masks%chi2(i,50))*reference_parameters%fields%k0*reference_parameters%fields%tau0,& 
!!$             zones(k)%species(0)%var(1)%temperature(i,50)*(1.D0-zones(k)%masks%chi2(i,50))
!!$     end do
!!$     k=25
!!$     do i=1,zones(k)%mesh%Nx
!!$        write(100,101) zones(k)%species(0)%var(1)%density(i,50)*(1.D0-zones(k)%masks%chi2(i,50)),&
!!$             zones(k)%electric_fields(1)%phi(i,50)*(1.D0-zones(k)%masks%chi2(i,50)),&
!!$             zones(k)%electric_fields(1)%vorticity(i,50)*(1.D0-zones(k)%masks%chi2(i,50)),&
!!$             0.09*zones(k)%kepsilon(1)%k(i,50)**2/zones(k)%kepsilon(1)%epsilon(i,50)&
!!$             *(1.D0-zones(k)%masks%chi2(i,50))*reference_parameters%fields%k0*reference_parameters%fields%tau0,& 
!!$             zones(k)%species(0)%var(1)%temperature(i,50)*(1.D0-zones(k)%masks%chi2(i,50))
!!$     end do
!!$     close(100)
!!$  end if
!!$
!!$  if(modulo(dble(n_ite),dble(global_parameters%N_iterations)/1000.D0).eq.0) then
!!$     open(unit=100,file='custom3',status='unknown')
!!$     k=2
!!$     do j=1,zones(k)%mesh%Nz
!!$        write(100,101) zones(k)%species(0)%var(1)%density(20,j)*(1.D0-zones(k)%masks%chi2(20,j)),&
!!$             zones(k)%electric_fields(1)%phi(20,j)*(1.D0-zones(k)%masks%chi2(20,j)),&
!!$             zones(k)%electric_fields(1)%vorticity(20,j)*(1.D0-zones(k)%masks%chi2(20,j)),&
!!$             0.09*zones(k)%kepsilon(1)%k(20,j)**2/zones(k)%kepsilon(1)%epsilon(20,j)&
!!$             *(1.D0-zones(k)%masks%chi2(20,j))*reference_parameters%fields%k0*reference_parameters%fields%tau0,& 
!!$             zones(k)%species(0)%var(1)%temperature(20,j)*(1.D0-zones(k)%masks%chi2(20,j))
!!$     end do
!!$     k=12
!!$     do j=1,zones(k)%mesh%Nz
!!$        write(100,101) zones(k)%species(0)%var(1)%density(20,j)*(1.D0-zones(k)%masks%chi2(20,j)),&
!!$             zones(k)%electric_fields(1)%phi(20,j)*(1.D0-zones(k)%masks%chi2(20,j)),&
!!$             zones(k)%electric_fields(1)%vorticity(20,j)*(1.D0-zones(k)%masks%chi2(20,j)),&
!!$             0.09*zones(k)%kepsilon(1)%k(20,j)**2/zones(k)%kepsilon(1)%epsilon(20,j)&
!!$             *(1.D0-zones(k)%masks%chi2(20,j))*reference_parameters%fields%k0*reference_parameters%fields%tau0,& 
!!$             zones(k)%species(0)%var(1)%temperature(20,j)*(1.D0-zones(k)%masks%chi2(20,j))
!!$     end do
!!$     k=13
!!$     do j=1,zones(k)%mesh%Nz
!!$        write(100,101) zones(k)%species(0)%var(1)%density(20,j)*(1.D0-zones(k)%masks%chi2(20,j)),&
!!$             zones(k)%electric_fields(1)%phi(20,j)*(1.D0-zones(k)%masks%chi2(20,j)),&
!!$             zones(k)%electric_fields(1)%vorticity(20,j)*(1.D0-zones(k)%masks%chi2(20,j)),&
!!$             0.09*zones(k)%kepsilon(1)%k(20,j)**2/zones(k)%kepsilon(1)%epsilon(20,j)&
!!$             *(1.D0-zones(k)%masks%chi2(20,j))*reference_parameters%fields%k0*reference_parameters%fields%tau0,& 
!!$             zones(k)%species(0)%var(1)%temperature(20,j)*(1.D0-zones(k)%masks%chi2(20,j))
!!$     end do
!!$     k=14
!!$     do j=1,zones(k)%mesh%Nz
!!$        write(100,101) zones(k)%species(0)%var(1)%density(20,j)*(1.D0-zones(k)%masks%chi2(20,j)),&
!!$             zones(k)%electric_fields(1)%phi(20,j)*(1.D0-zones(k)%masks%chi2(20,j)),&
!!$             zones(k)%electric_fields(1)%vorticity(20,j)*(1.D0-zones(k)%masks%chi2(20,j)),&
!!$             0.09*zones(k)%kepsilon(1)%k(20,j)**2/zones(k)%kepsilon(1)%epsilon(20,j)&
!!$             *(1.D0-zones(k)%masks%chi2(20,j))*reference_parameters%fields%k0*reference_parameters%fields%tau0,& 
!!$             zones(k)%species(0)%var(1)%temperature(20,j)*(1.D0-zones(k)%masks%chi2(20,j))
!!$     end do
!!$     close(100)
!!$  end if



!!$  if(modulo(dble(n_ite),dble(global_parameters%N_iterations)/1000.D0).eq.0) then
!!$     open(unit=100,file='custom',status='unknown')
!!$     k=13 
!!$     do i=1,zones(k)%mesh%Nx
!!$        write(100,101) zones(k)%species(0)%var(1)%density(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%electric_fields(1)%phi(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%electric_fields(1)%vorticity(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%kepsilon(1)%mu_t(i,24)*(reference_parameters%geometry%rs0*reference_parameters%fields%c0)&
!!$             *(1.D0-zones(k)%masks%chi2(i,24)),& 
!!$             zones(k)%kepsilon(1)%k(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%species(0)%var(1)%temperature(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%kepsilon(1)%interchange(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%kepsilon(1)%UE_shear(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%species(1)%transport_perp%v_pinch(i,24),&
!!$             zones(k)%kepsilon(1)%epsilon(i,24)*(1.D0-zones(k)%masks%chi2(i,24))
!!$     end do
!!$     k=24
!!$     do i=1,zones(k)%mesh%Nx
!!$        write(100,101) zones(k)%species(0)%var(1)%density(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%electric_fields(1)%phi(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%electric_fields(1)%vorticity(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%kepsilon(1)%mu_t(i,24)*(reference_parameters%geometry%rs0*reference_parameters%fields%c0)&
!!$             *(1.D0-zones(k)%masks%chi2(i,24)),& 
!!$             zones(k)%kepsilon(1)%k(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%species(0)%var(1)%temperature(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%kepsilon(1)%interchange(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%kepsilon(1)%UE_shear(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%species(1)%transport_perp%v_pinch(i,24),&
!!$             zones(k)%kepsilon(1)%epsilon(i,24)*(1.D0-zones(k)%masks%chi2(i,24))
!!$     end do
!!$     k=12
!!$     do i=1,zones(k)%mesh%Nx
!!$        write(100,101) zones(k)%species(0)%var(1)%density(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%electric_fields(1)%phi(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%electric_fields(1)%vorticity(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%kepsilon(1)%mu_t(i,24)*(reference_parameters%geometry%rs0*reference_parameters%fields%c0)&
!!$             *(1.D0-zones(k)%masks%chi2(i,24)),& 
!!$             zones(k)%kepsilon(1)%k(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%species(0)%var(1)%temperature(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%kepsilon(1)%interchange(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%kepsilon(1)%UE_shear(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%species(1)%transport_perp%v_pinch(i,24),&
!!$             zones(k)%kepsilon(1)%epsilon(i,24)*(1.D0-zones(k)%masks%chi2(i,24))
!!$     end do
!!$     k=11
!!$     do i=1,zones(k)%mesh%Nx
!!$        write(100,101) zones(k)%species(0)%var(1)%density(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%electric_fields(1)%phi(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%electric_fields(1)%vorticity(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%kepsilon(1)%mu_t(i,24)*(reference_parameters%geometry%rs0*reference_parameters%fields%c0)&
!!$             *(1.D0-zones(k)%masks%chi2(i,24)),& 
!!$             zones(k)%kepsilon(1)%k(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%species(0)%var(1)%temperature(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%kepsilon(1)%interchange(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%kepsilon(1)%UE_shear(i,24)*(1.D0-zones(k)%masks%chi2(i,24)),&
!!$             zones(k)%species(1)%transport_perp%v_pinch(i,24),&
!!$             zones(k)%kepsilon(1)%epsilon(i,24)*(1.D0-zones(k)%masks%chi2(i,24))
!!$     end do
!!$
!!$101  format(512es15.7) 
!!$     close(100)
!!$
!     open(unit=100,file='customp',status='unknown')
!     k=23
!     do j=1,zones(k)%mesh%Nz
!        write(100,101) zones(k)%species(0)%var(1)%density(20,j)*(1.D0-zones(k)%masks%chi2(20,j)),&
!             zones(k)%electric_fields(1)%phi(20,j),&
!             zones(k)%electric_fields(1)%vorticity(20,j)*(1.D0-zones(k)%masks%chi2(20,j))!,&
!!$             zones(k)%kepsilon(1)%mu_t(20,j)*(reference_parameters%geometry%rs0*reference_parameters%fields%c0)&
!!$             *(1.D0-zones(k)%masks%chi2(20,j)),& 
!!$             zones(k)%kepsilon(1)%k(20,j)*(1.D0-zones(k)%masks%chi2(20,j)),&
!!$             zones(k)%species(0)%var(1)%temperature(20,j)*(1.D0-zones(k)%masks%chi2(20,j)),&
!!$             zones(k)%kepsilon(1)%interchange(20,j)*(1.D0-zones(k)%masks%chi2(20,j)),&
!!$             zones(k)%kepsilon(1)%UE_shear(20,j)*(1.D0-zones(k)%masks%chi2(20,j))
!     end do
!     k=24
!     do j=1,zones(k)%mesh%Nz
!        write(100,101) zones(k)%species(0)%var(1)%density(20,j)*(1.D0-zones(k)%masks%chi2(20,j)),&
!             zones(k)%electric_fields(1)%phi(20,j),&
!             zones(k)%electric_fields(1)%vorticity(20,j)*(1.D0-zones(k)%masks%chi2(20,j))!,&
!!$             zones(k)%kepsilon(1)%mu_t(20,j)*(reference_parameters%geometry%rs0*reference_parameters%fields%c0)&
!!$             *(1.D0-zones(k)%masks%chi2(20,j)),& 
!!$             zones(k)%kepsilon(1)%k(20,j)*(1.D0-zones(k)%masks%chi2(20,j)),&
!!$             zones(k)%species(0)%var(1)%temperature(20,j)*(1.D0-zones(k)%masks%chi2(20,j)),&
!!$             zones(k)%kepsilon(1)%interchange(20,j)*(1.D0-zones(k)%masks%chi2(20,j)),&
!!$             zones(k)%kepsilon(1)%UE_shear(20,j)*(1.D0-zones(k)%masks%chi2(20,j))
!     end do
!     close(100)
!  end if
!!$

end subroutine save_custom
