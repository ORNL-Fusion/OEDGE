subroutine compute_Zeff()
  use all_variables, only : interp_data2, global_parameters
  implicit none
  integer*4 :: n,k
  real*8 :: ne
  Interp_Data2%knots_Zeff=0.D0
  do k=1,Interp_Data2%N_knots
     ne=0.d0
     do n=1,global_parameters%N_ions
        Interp_Data2%knots_Zeff(k)=Interp_Data2%knots_Zeff(k)&
             +Interp_Data2%knots_density(k,n,1)*global_parameters%ions_list(n,2)**2
     end do
     do n=1,global_parameters%N_ions
        ne=ne&
             +Interp_Data2%knots_density(k,n,1)*global_parameters%ions_list(n,2)
     end do
     Interp_Data2%knots_Zeff(k)=Interp_Data2%knots_Zeff(k)/ne
  end do
end subroutine compute_Zeff
