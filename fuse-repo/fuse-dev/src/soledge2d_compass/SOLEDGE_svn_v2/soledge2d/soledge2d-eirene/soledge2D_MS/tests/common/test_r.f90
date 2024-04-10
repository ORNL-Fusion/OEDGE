function test_r(x)
  use test_var
  use Mphysics
  implicit none
  real*8,intent(in) :: x
  real*8 :: test_r
  test_r=test_a*(0.6D0+test_reg_r*0.8D0*x-(1.D0-test_reg_r)*0.8*(cos(x*pi)-1.D0)*0.5D0)
end function test_r
