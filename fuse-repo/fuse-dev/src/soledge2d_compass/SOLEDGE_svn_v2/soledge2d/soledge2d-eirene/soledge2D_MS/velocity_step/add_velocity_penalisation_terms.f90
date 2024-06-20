subroutine add_velocity_penalisation_terms(zone)
  use all_variables, only : global_parameters, reference_parameters,&
       global_variables, penalisation_parameters
  use Mzone
  use Mphysics
  implicit none
  type(Tzone) :: zone
  integer*4 :: i,j,n
  integer*4 :: Nx,Nz
  real*8 :: dt
  real*8 :: alpham,alphap
  real*8 :: alpham1,alphap1
  real*8 :: error
  real*8 :: cs,m_i
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  dt=global_variables%dt
  do i=1,Nx
     do j=1,Nz
        do n=1,global_parameters%N_ions
           call add_velocity_penalisation_terms_point(zone,i,j,n)
        end do
     end do
  end do
end subroutine add_velocity_penalisation_terms


subroutine add_velocity_penalisation_terms_point(zone,i,j,n)
  use all_variables, only : global_parameters, reference_parameters,&
       global_variables, penalisation_parameters, flags
  use Mzone
  use Mphysics
  implicit none
  type(Tzone) :: zone
  integer*4,intent(in) :: i,j,n
  integer*4 :: Nx,Nz
  real*8 :: dt
  real*8 :: alpham,alphap
  real*8 :: alpham1,alphap1
  real*8 :: error
  real*8 :: cs,m_i,Mem, Mep, R0, c0, rs0
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  dt=global_variables%dt
  m_i=zone%species(n)%element%mass
  R0 = reference_parameters%geometry%R0
  rs0 = reference_parameters%geometry%rs0
  c0 = reference_parameters%fields%c0
 
  Mem = (zone%species(n)%drifts%uEt(i,j-1)+zone%species(n)%drifts%uBt(i,j-1))&
       *sqrt(zone%metric_coefficients%ctt(i,j-1))/(zone%metric_coefficients%G(i,j-1)&
       *sqrt((zone%species(n)%charge*zone%species(0)%var(1)%temperature(i,j-1)+&
       zone%species(0)%var(1)%temperature(i,j-1))/m_i))*(2.d0*pi*R0/rs0)
  Mep = (zone%species(n)%drifts%uEt(i,j+1)+zone%species(n)%drifts%uBt(i,j+1))&
       *sqrt(zone%metric_coefficients%ctt(i,j+1))/(zone%metric_coefficients%G(i,j+1)&
       *sqrt((zone%species(n)%charge*zone%species(0)%var(1)%temperature(i,j+1)+&
       zone%species(0)%var(1)%temperature(i,j+1))/m_i))*(2.d0*pi*R0/rs0)

  m_i=zone%species(n)%element%mass
  error=1.d0-(zone%species(n)%var(1)%Mach(i,j-1)+Mem)*global_variables%sign_metric
  alpham1=1.d0-Mem*global_variables%sign_metric + error*penalisation_parameters%gain
  zone%species(n)%penalisation_memories%alpham(i,j)=zone%species(n)%penalisation_memories%alpham(i,j)&
       *dble(penalisation_parameters%keep-1)/dble(penalisation_parameters%keep)+error
  alpham=alpham1+zone%species(n)%penalisation_memories%alpham(i,j)/penalisation_parameters%dump_pen*penalisation_parameters%gain
  alpham=max(alpham,1.D0-Mem*global_variables%sign_metric)
  alpham=min(alpham,2.d0)
  error=1.d0+(zone%species(n)%var(1)%mach(i,j+1)+Mep)*global_variables%sign_metric
  alphap1=1.d0 + Mep*global_variables%sign_metric + error*penalisation_parameters%gain
  zone%species(n)%penalisation_memories%alphap(i,j)=zone%species(n)%penalisation_memories%alphap(i,j)&
       *dble(penalisation_parameters%keep-1)/dble(penalisation_parameters%keep)+error
  alphap=alphap1+zone%species(n)%penalisation_memories%alphap(i,j)/penalisation_parameters%dump_pen*penalisation_parameters%gain
  alphap=max(alphap,1.D0+Mep*global_variables%sign_metric)
  alphap=min(alphap,2.d0)

  zone%species(n)%tridiag%S(i,j)=zone%species(n)%tridiag%S(i,j)&
                                !penalisation1
       +zone%masks%chi1(i,j)/penalisation_parameters%eta&
       *zone%species(n)%var(1)%density(i,j)&
       *sqrt((zone%species(n)%charge*zone%species(0)%var(1)%temperature(i,j)+zone%species(n)%var(1)%temperature(i,j))&
       /m_i)*global_variables%sign_metric*max(alpham,abs(zone%species(n)%var(1)%Mach(i,j-1)))&
                                !penalisation3
       -zone%masks%chi3(i,j)/penalisation_parameters%eta&
       *zone%species(n)%var(1)%density(i,j)&
       *sqrt((zone%species(n)%charge*zone%species(0)%var(1)%temperature(i,j)+zone%species(n)%var(1)%temperature(i,j))&
       /m_i)*global_variables%sign_metric*max(alphap,abs(zone%species(n)%var(1)%Mach(i,j+1)))
  zone%species(n)%tridiag%b(i,j)=zone%species(n)%tridiag%b(i,j)&
                                !penalisation
       +(zone%masks%chi4(i,j)/penalisation_parameters%eta2**2&
       +zone%masks%chi1(i,j)/penalisation_parameters%eta&
       +zone%masks%chi3(i,j)/penalisation_parameters%eta)
end subroutine add_velocity_penalisation_terms_point
