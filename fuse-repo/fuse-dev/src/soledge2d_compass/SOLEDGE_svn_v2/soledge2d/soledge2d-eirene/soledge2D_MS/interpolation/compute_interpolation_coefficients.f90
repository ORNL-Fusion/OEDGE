subroutine compute_interpolation_coefficients()
  use all_variables, only : interp_data2, zones
  implicit none
  integer*4 :: n
  integer*4 :: i,j,k
  real*8 :: R1,Z1,R2,Z2,R3,Z3,det

    ! computation of interpolation coefs
  do n=1,Interp_data2%N_knots
     if(Interp_data2%knots_interp_points(n)%pass.ne.1) then
        if(Interp_data2%knots_interp_points(n)%pass.ne.5) then
           if(Interp_data2%knots_interp_points(n)%pass.eq.2) then
              i=Interp_data2%knots_interp_points(n)%sol(1,1)
              j=Interp_data2%knots_interp_points(n)%sol(1,2)
              k=Interp_data2%knots_interp_points(n)%sol(1,3)
              R1=zones(k)%mesh%Rgeom(i,j)*100.D0
              Z1=zones(k)%mesh%Zgeom(i,j)*100.D0
              i=Interp_data2%knots_interp_points(n)%sol(2,1)
              j=Interp_data2%knots_interp_points(n)%sol(2,2)
              k=Interp_data2%knots_interp_points(n)%sol(2,3)
              R2=zones(k)%mesh%Rgeom(i,j)*100.D0
              Z2=zones(k)%mesh%Zgeom(i,j)*100.D0
              i=Interp_data2%knots_interp_points(n)%sol(3,1)
              j=Interp_data2%knots_interp_points(n)%sol(3,2)
              k=Interp_data2%knots_interp_points(n)%sol(3,3)
              R3=zones(k)%mesh%Rgeom(i,j)*100.D0
              Z3=zones(k)%mesh%Zgeom(i,j)*100.D0
           else
              if(Interp_data2%knots_interp_points(n)%pass.ge.3) then
                 if(Interp_data2%knots_interp_points(n)%n_soledge.eq.1) then
                    i=Interp_data2%knots_interp_points(n)%sol(1,1)
                    j=Interp_data2%knots_interp_points(n)%sol(1,2)
                    k=Interp_data2%knots_interp_points(n)%sol(1,3)
                    R1=zones(k)%mesh%Rgeom(i,j)*100.D0
                    Z1=zones(k)%mesh%Zgeom(i,j)*100.D0
                    R2=Interp_data2%knots_R(Interp_data2%knots_interp_points(n)%eir(1))
                    Z2=Interp_data2%knots_Z(Interp_data2%knots_interp_points(n)%eir(1))
                    R3=Interp_data2%knots_R(Interp_data2%knots_interp_points(n)%eir(2))
                    Z3=Interp_data2%knots_Z(Interp_data2%knots_interp_points(n)%eir(2))
                 else
                    i=Interp_data2%knots_interp_points(n)%sol(1,1)
                    j=Interp_data2%knots_interp_points(n)%sol(1,2)
                    k=Interp_data2%knots_interp_points(n)%sol(1,3)
                    R1=zones(k)%mesh%Rgeom(i,j)*100.D0
                    Z1=zones(k)%mesh%Zgeom(i,j)*100.D0
                    i=Interp_data2%knots_interp_points(n)%sol(2,1)
                    j=Interp_data2%knots_interp_points(n)%sol(2,2)
                    k=Interp_data2%knots_interp_points(n)%sol(2,3)
                    R2=zones(k)%mesh%Rgeom(i,j)*100.D0
                    Z2=zones(k)%mesh%Zgeom(i,j)*100.D0
                    R3=Interp_data2%knots_R(Interp_data2%knots_interp_points(n)%eir(1))
                    Z3=Interp_data2%knots_Z(Interp_data2%knots_interp_points(n)%eir(1))
                 end if
              end if
           end if
           det=R1*Z2+Z1*R3+R2*Z3-R3*Z2-Z3*R1-R2*Z1
!           write(*,*) det, n, R1, Z1 , R2, Z2, R3, Z3
           Interp_data2%knots_interp_points(n)%interp_matrix(1,1)=1/det*(Z2-Z3)
           Interp_data2%knots_interp_points(n)%interp_matrix(1,2)=1/det*(Z3-Z1)
           Interp_data2%knots_interp_points(n)%interp_matrix(1,3)=1/det*(Z1-Z2)
           Interp_data2%knots_interp_points(n)%interp_matrix(2,1)=1/det*(R3-R2)
           Interp_data2%knots_interp_points(n)%interp_matrix(2,2)=1/det*(R1-R3)
           Interp_data2%knots_interp_points(n)%interp_matrix(2,3)=1/det*(R2-R1)
           Interp_data2%knots_interp_points(n)%interp_matrix(3,1)=1/det*(R2*Z3-R3*Z2)
           Interp_data2%knots_interp_points(n)%interp_matrix(3,2)=1/det*(R3*Z1-R1*Z3)
           Interp_data2%knots_interp_points(n)%interp_matrix(3,3)=1/det*(R1*Z2-R2*Z1)
        end if
     end if
  end do
end subroutine compute_interpolation_coefficients
