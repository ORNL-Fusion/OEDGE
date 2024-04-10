subroutine make_Xpoint_cells_orthogonal()
  use all_variables, only : global_parameters, zones
  implicit none
  integer*4 :: i,j,k
  integer*4 :: Nx,Nz
  integer*4 :: Z1,Z2
  integer*4 :: North, South, East, West
  do k=1,global_parameters%N_zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     !NE corner
     North=zones(k)%Neighbors(1)
     if(North.gt.0) then
        Z1=zones(North)%Neighbors(3)
     else
        Z1=-1
     end if
     East=zones(k)%Neighbors(3)
     if(East.gt.0) then
        Z2=zones(East)%Neighbors(1)
     else
        Z2=-1
     end if
     if((Z1.gt.0).and.(Z2.gt.0).and.(Z1.ne.Z2)) then
        !Xpoint detected - set cpt to zero
        zones(k)%metric_coefficients%cpt(Nx,Nz)=0.d0
        zones(k)%metric_coefficients%cpt(Nx+1,Nz)=0.d0
        zones(k)%metric_coefficients%cpt(Nx,Nz+1)=0.d0
        zones(k)%metric_coefficients%cpt(Nx+1,Nz+1)=0.d0
     end if
     !NW corner
     North=zones(k)%Neighbors(1)
     if(North.gt.0) then
        Z1=zones(North)%Neighbors(4)
     else
        Z1=-1
     end if
     West=zones(k)%Neighbors(4)
     if(West.gt.0) then
        Z2=zones(West)%Neighbors(1)
     else
        Z2=-1
     end if
     if((Z1.gt.0).and.(Z2.gt.0).and.(Z1.ne.Z2)) then
        !Xpoint detected - set cpt to zero
        zones(k)%metric_coefficients%cpt(Nx,1)=0.d0
        zones(k)%metric_coefficients%cpt(Nx+1,1)=0.d0
        zones(k)%metric_coefficients%cpt(Nx,0)=0.d0
        zones(k)%metric_coefficients%cpt(Nx+1,0)=0.d0
     end if
     !SE corner
     South=zones(k)%Neighbors(2)
     if(South.gt.0) then
        Z1=zones(South)%Neighbors(3)
     else
        Z1=-1
     end if
     East=zones(k)%Neighbors(3)
     if(East.gt.0) then
        Z2=zones(East)%Neighbors(2)
     else
        Z2=-1
     end if
     if((Z1.gt.0).and.(Z2.gt.0).and.(Z1.ne.Z2)) then
        !Xpoint detected - set cpt to zero
        zones(k)%metric_coefficients%cpt(1,Nz)=0.d0
        zones(k)%metric_coefficients%cpt(1,Nz+1)=0.d0
        zones(k)%metric_coefficients%cpt(0,Nz)=0.d0
        zones(k)%metric_coefficients%cpt(0,Nz+1)=0.d0
     end if
     !SW corner
     South=zones(k)%Neighbors(2)
     if(South.gt.0) then
        Z1=zones(South)%Neighbors(4)
     else
        Z1=-1
     end if
     West=zones(k)%Neighbors(4)
     if(West.gt.0) then
        Z2=zones(West)%Neighbors(2)
     else
        Z2=-1
     end if
     if((Z1.gt.0).and.(Z2.gt.0).and.(Z1.ne.Z2)) then
        !Xpoint detected - set cpt to zero
        zones(k)%metric_coefficients%cpt(1,1)=0.d0
        zones(k)%metric_coefficients%cpt(1,0)=0.d0
        zones(k)%metric_coefficients%cpt(0,1)=0.d0
        zones(k)%metric_coefficients%cpt(0,0)=0.d0
     end if
  end do
end subroutine make_Xpoint_cells_orthogonal
