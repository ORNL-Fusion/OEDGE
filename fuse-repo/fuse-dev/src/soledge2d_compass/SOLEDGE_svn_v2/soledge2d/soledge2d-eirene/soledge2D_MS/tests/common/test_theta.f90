function test_theta(z,x)
  use test_var
  use Mphysics
  implicit none
  real*8,intent(in) :: z
  real*8,intent(in) :: x
  real*8 :: test_theta
  test_theta=test_reg_theta*z*2.D0*pi-(1.D0-test_reg_theta)*pi*(cos(z*pi)-1)+test_theta_shift*pi*x
!  test_theta=2.D0*pi*z+test_theta_shift*pi*x
end function test_theta
