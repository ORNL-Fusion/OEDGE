subroutine set_gaussian_ballooning()
  use all_variables, only : zones, ballooning_parameters, global_parameters
  use Mphysics
  implicit none
  integer*4 :: k,n
  integer*4 :: Nx,Nz
  real*8, allocatable :: ball(:,:),zp(:,:)
  real*8 abal,bbal
  bbal=ballooning_parameters%minmaxbal/(1.-ballooning_parameters%minmaxbal)
  abal=1./(bbal+ballooning_parameters%sigmabal*sqrt(pi))
  do k=1,global_parameters%N_zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     allocate(ball(0:Nx+1,0:Nz+1))
     allocate(zp(0:Nx+1,0:Nz+1))
     zp=modulo(zones(k)%mesh%z+(0.5-ballooning_parameters%zbal),1.)-0.5
     ball=abal*(exp(-(zp/ballooning_parameters%sigmabal)**2.)+bbal)
     zones(k)%species(0)%transport_perp%chi_p=zones(k)%species(0)%transport_perp%chi_p*ball
     zones(k)%species(0)%transport_perp%chi_t=zones(k)%species(0)%transport_perp%chi_t*ball
     do n=1,global_parameters%N_ions
        zones(k)%species(n)%transport_perp%D_p=zones(k)%species(n)%transport_perp%D_p*ball
        zones(k)%species(n)%transport_perp%D_t=zones(k)%species(n)%transport_perp%D_t*ball
        zones(k)%species(n)%transport_perp%nu_p=zones(k)%species(n)%transport_perp%nu_p*ball
        zones(k)%species(n)%transport_perp%nu_t=zones(k)%species(n)%transport_perp%nu_t*ball
        zones(k)%species(n)%transport_perp%chi_p=zones(k)%species(n)%transport_perp%chi_p*ball
        zones(k)%species(n)%transport_perp%chi_t=zones(k)%species(n)%transport_perp%chi_t*ball
     end do
     deallocate(zp,ball)
  end do
end subroutine set_gaussian_ballooning
