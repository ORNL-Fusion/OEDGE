subroutine add_test_source(zone,FIELD)
  use test_var
  use all_variables, only : global_parameters
  use MZone
  use Mdefinitions
  implicit none
  type(Tzone),intent(inout) :: zone
  integer*4,intent(in) :: FIELD
  integer*4 :: i,j,n,k
  integer*4 :: Nx,Nz
  real*8 :: dt
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  k=zone%number
  do i=1,Nx
     do j=1,Nz
        do n=0,global_parameters%N_ions
           select case(FIELD)
           case(DENSITY_FIELD)
              zone%species(n)%tridiag%S(i,j)=zone%species(n)%tridiag%S(i,j)+&
                   test_sources(k,n)%Sn(i,j)
           case(VELOCITY_FIELD)
              zone%species(n)%tridiag%S(i,j)=zone%species(n)%tridiag%S(i,j)+&
                   test_sources(k,n)%SG(i,j)
           case(TEMPERATURE_FIELD)
              zone%species(n)%tridiag%S(i,j)=zone%species(n)%tridiag%S(i,j)+&
                   test_sources(k,n)%ST(i,j)
           end select
        end do
     end do
  end do
end subroutine add_test_source
