subroutine find_closest_point_ultimate_interp()
  use all_variables, only : interp_data2, zones
  implicit none
  integer*4 n,n2,i,j,k,npt
  real*8 d,dmin
  do n=1,Interp_data2%N_knots
     dmin=1.d15
     do n2=1,Interp_data2%knots_interp_points(n)%n_soledge
        i=Interp_data2%knots_interp_points(n)%sol(n2,1)
        j=Interp_data2%knots_interp_points(n)%sol(n2,2)
        k=Interp_data2%knots_interp_points(n)%sol(n2,3)
        d=sqrt((Interp_data2%knots_R(n)-zones(k)%mesh%Rgeom(i,j)*100.D0)**2.D0+&
             (Interp_data2%knots_Z(n)-zones(k)%mesh%Zgeom(i,j)*100.D0)**2.D0)
        if(d.lt.dmin) then
           Interp_data2%knots_interp_points(n)%type_ultimate_point=1 ! Soledge point as closest
           Interp_data2%knots_interp_points(n)%ultimate_coords(1)=i
           Interp_data2%knots_interp_points(n)%ultimate_coords(2)=j
           Interp_data2%knots_interp_points(n)%ultimate_coords(3)=k
        end if
     end do
     do n2=1,Interp_data2%knots_interp_points(n)%n_eirene
        npt=Interp_data2%knots_interp_points(n)%eir(n2)
        d=sqrt((Interp_data2%knots_R(n)-Interp_data2%knots_R(npt))**2.D0+&
             (Interp_data2%knots_Z(n)-Interp_data2%knots_Z(npt))**2.D0)
        if(d.lt.dmin) then
           if((Interp_data2%knots_interp_points(npt)%pass.eq.1).or.(Interp_data2%knots_interp_points(npt)%pass.eq.2)) then
              Interp_data2%knots_interp_points(n)%type_ultimate_point=2 ! Eirene point as closest
              Interp_data2%knots_interp_points(n)%ultimate_coords(1)=npt
           end if
        end if
     end do
  end do
end subroutine find_closest_point_ultimate_interp
