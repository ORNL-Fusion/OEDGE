subroutine project(n,X1,X2,X3,inter)
  use all_variables, only : interp_data2
  implicit none
  real*8,intent(in) :: X1,X2,X3
  integer*4,intent(in) :: n
  real*8,intent(out) :: inter
  real*8 a,b,c
  a=interp_data2%knots_interp_points(n)%interp_matrix(1,1)*X1&
       +interp_data2%knots_interp_points(n)%interp_matrix(1,2)*X2&
       +interp_data2%knots_interp_points(n)%interp_matrix(1,3)*X3
  b=interp_data2%knots_interp_points(n)%interp_matrix(2,1)*X1&
       +interp_data2%knots_interp_points(n)%interp_matrix(2,2)*X2&
       +interp_data2%knots_interp_points(n)%interp_matrix(2,3)*X3
  c=interp_data2%knots_interp_points(n)%interp_matrix(3,1)*X1&
       +interp_data2%knots_interp_points(n)%interp_matrix(3,2)*X2&
       +interp_data2%knots_interp_points(n)%interp_matrix(3,3)*X3
  inter=a*interp_data2%knots_R(n)+b*interp_data2%knots_Z(n)+c
end subroutine project
