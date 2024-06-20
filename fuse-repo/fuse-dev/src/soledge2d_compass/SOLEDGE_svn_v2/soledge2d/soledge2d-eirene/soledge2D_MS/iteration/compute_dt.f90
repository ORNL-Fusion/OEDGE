subroutine compute_dt(zone,dt)
  use all_variables, only : global_parameters,reference_parameters
  use Mzone
  use Mphysics
  implicit none
  Type(TZone),intent(in) :: zone
  real*8,intent(out) :: dt
  integer*4 :: Nx,Nz
  integer*4 :: i,j
  integer*4 :: n
  real*8 :: csad, vel
  real*8 :: dtz,dtx
  integer*4 :: dt_type, lim_i, lim_j
  real*8 :: R0,c0,eps,rs0
  dt=1.e6
  dt_type=0
  lim_i=0
  lim_j=0
  R0=reference_parameters%geometry%R0
  rs0=reference_parameters%geometry%rs0
  c0=reference_parameters%fields%c0
  eps=1.d-10
  ! advection terms
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do n=1,global_parameters%N_ions
     do i=1,Nx
        do j=1,Nz
           if(zone%masks%chi4(i,j).eq.0) then
              vel = zone%species(n)%var(1)%velocity(i,j)+(zone%species(n)%drifts%uEt(i,j)&
                   +zone%species(n)%drifts%uBt(i,j))*sqrt(zone%metric_coefficients%ctt(i,j))*(2.d0*pi*R0/rs0)/zone%metric_coefficients%G(i,j)
              csad  = sqrt((zone%species(n)%charge*zone%species(0)%var(1)%temperature(i,j)+zone%species(n)%var(1)%temperature(i,j))&
                   /zone%species(n)%element%mass)
              dtz=abs((zone%mesh%z_plus_1half(i,j)-zone%mesh%z_minus_1half(i,j))&
                   /(zone%metric_coefficients%G(i,j)*(csad+abs(vel))))
              if(dt>dtz) then
                 dt = dtz
                 dt_type=1
                 lim_i=i
                 lim_j=j
              end if
           end if
        end do
     end do
  end do
  ! Radial diffusion terms
  do i=1,Nx
     do j=1,Nz
        if(zone%masks%chi4(i,j).eq.0) then
           dtx=(0.5d0*(zone%mesh%x(i+1,j)-zone%mesh%x(i-1,j)))**(2.d0)&
                /(zone%metric_coefficients%cpp(i,j)&
                *zone%species(1)%transport_perp%D_p(i,j))
           if(dt>dtx) then
              dt = dtx
              dt_type=2
              lim_i=i
              lim_j=j
           end if
        end if
     end do
  end do
  ! Radial drifts
  do i=1,Nx
     do j=1,Nz
        if(zone%masks%chi2(i,j).eq.0) then
           dtx=0.5d0*(zone%mesh%x(i+1,j)-zone%mesh%x(i-1,j))&
                /(abs((zone%species(1)%drifts%uEp(i,j)&
                +zone%species(1)%drifts%uBp(i,j))*R0/c0)+eps)
           if(dt>dtx) then
              dt = dtx
              dt_type=3
              lim_i=i
              lim_j=j
           end if
        end if
     end do
  end do
!  write(*,100) zone%number, dt, dt_type, lim_i, lim_j
100 format('zone ',I3,es15.7,I3,I4,I4)
end subroutine compute_dt
