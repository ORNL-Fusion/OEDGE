subroutine MD_broadcast_corners_B()
  use all_variables, only : global_parameters, zones
  implicit none
  integer*4 :: n
  integer*4 :: mNorth,mSouth,mEast,mWest
  integer*4 :: North,South,East,West
  integer*4 :: Nx_N,Nz_N
  integer*4 :: i,j,k,Nx,Nz
  do k=1,global_parameters%N_zones

     mNorth=zones(k)%MagNeighbors(1)
     mSouth=zones(k)%MagNeighbors(2)
     mEast=zones(k)%MagNeighbors(3)
     mWest=zones(k)%MagNeighbors(4)

     North=zones(k)%Neighbors(1)
     South=zones(k)%Neighbors(2)
     East=zones(k)%Neighbors(3)
     West=zones(k)%Neighbors(4)


     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz

     !treating corner NE
     if(mNorth.eq.0) then
        if(mEast.eq.0) then
           zones(k)%mesh%B(Nx+1,Nz+1)=0.5D0*(&
                zones(North)%mesh%B(1,Nz+1)&
                +zones(East)%mesh%B(Nx+1,1))
        else
           zones(k)%mesh%B(Nx+1,Nz+1)=(&
                zones(North)%mesh%B(1,Nz+1))
        end if
     else
        if(mEast.eq.0) then
           zones(k)%mesh%B(Nx+1,Nz+1)=(&
                zones(East)%mesh%B(Nx+1,1))
        else !real corner
           zones(k)%mesh%B(Nx+1,Nz+1)=0.5D0*(&
                zones(k)%mesh%B(Nx+1,Nz)&
                +zones(k)%mesh%B(Nx,Nz+1))
        end if
     end if

     !treating corner NW
     if(mNorth.eq.0) then
        if(mWest.eq.0) then
           zones(k)%mesh%B(Nx+1,0)=0.5D0*(&
                zones(North)%mesh%B(1,0)&
                +zones(West)%mesh%B(Nx+1,Zones(West)%mesh%Nz))
        else
           zones(k)%mesh%B(Nx+1,0)=(&
                zones(North)%mesh%B(1,0))
        end if
     else
        if(mWest.eq.0) then
           zones(k)%mesh%B(Nx+1,0)=(&
                zones(West)%mesh%B(Nx+1,Zones(West)%mesh%Nz))
        else !real corner
           zones(k)%mesh%B(Nx+1,0)=0.5D0*&
                (zones(k)%mesh%B(Nx+1,1)&
                +zones(k)%mesh%B(Nx,0))
        end if
     end if

     !treating corner SE
     if(mSouth.eq.0) then
        if(mEast.eq.0) then
           zones(k)%mesh%B(0,Nz+1)=0.5D0*(&
                zones(South)%mesh%B(Zones(South)%mesh%Nx,Nz+1)&
                +zones(East)%mesh%B(0,1))
        else
           zones(k)%mesh%B(0,Nz+1)=(&
                zones(South)%mesh%B(Zones(South)%mesh%Nx,Nz+1))
        end if
     else
        if(mEast.eq.0) then
           zones(k)%mesh%B(0,Nz+1)=(&
                zones(East)%mesh%B(0,1))
        else !real corner
           zones(k)%mesh%B(0,Nz+1)=0.5D0*&
                (zones(k)%mesh%B(0,Nz)&
                +zones(k)%mesh%B(1,Nz+1))
        end if
     end if

     !treating corner SW
     if(mSouth.eq.0) then
        if(mWest.eq.0) then
           zones(k)%mesh%B(0,0)=0.5D0*(&
                zones(South)%mesh%B(Zones(South)%mesh%Nx,0)&
                +zones(West)%mesh%B(0,Zones(West)%mesh%Nz))
        else
           zones(k)%mesh%B(0,0)=(&
                zones(South)%mesh%B(Zones(South)%mesh%Nx,0))
        end if
     else
        if(mWest.eq.0) then
           zones(k)%mesh%B(0,0)=(&
                zones(West)%mesh%B(0,Zones(West)%mesh%Nz))
        else !real corner
           zones(k)%mesh%B(0,0)=0.5d0*(&
                zones(k)%mesh%B(1,0)&
                +zones(k)%mesh%B(0,1))
        end if
     end if

!### Leybros modif ###!

     zones(k)%mesh%cornersB(2,1,1)=zones(k)%mesh%B(0,0)
     zones(k)%mesh%cornersB(1,1,1)=zones(k)%mesh%B(Nx+1,0)
     zones(k)%mesh%cornersB(1,2,1)=zones(k)%mesh%B(Nx+1,Nz+1)
     zones(k)%mesh%cornersB(2,2,1)=zones(k)%mesh%B(0,Nz+1)

     zones(k)%mesh%cornersB(2,1,2)=zones(k)%mesh%cornersB(2,1,1)
     zones(k)%mesh%cornersB(1,1,2)=zones(k)%mesh%cornersB(1,1,1)
     zones(k)%mesh%cornersB(1,2,2)=zones(k)%mesh%cornersB(1,2,1)
     zones(k)%mesh%cornersB(2,2,2)=zones(k)%mesh%cornersB(2,2,1)

  end do
end subroutine MD_broadcast_corners_B
