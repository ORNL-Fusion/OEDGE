subroutine Broadcast_drift3()
  Use all_variables, only : zones, reference_parameters, global_parameters
  use Mphysics
  implicit none
  integer*4 Nspec
  integer*4 North,South,East,West
  integer*4 mNorth,mSouth,mEast,mWest
  integer*4 k,Nx,Nz,j
  do Nspec=0,global_parameters%N_ions
     do k=1,global_parameters%N_zones
        Nx=zones(k)%mesh%Nx
        Nz=zones(k)%mesh%Nz
        ! on top
        North=zones(k)%Neighbors(1)
        mNorth=zones(k)%MagNeighbors(1)
        if(mNorth.eq.0) then
           ! if there is a neighboring domain, just broadcast
           zones(k)%species(Nspec)%drifts%uEp(Nx+1,1:Nz)=zones(North)%species(Nspec)%drifts%uEp(1,1:Nz)
           zones(k)%species(Nspec)%drifts%udp(Nx+1,1:Nz)=zones(North)%species(Nspec)%drifts%udp(1,1:Nz)
           zones(k)%species(Nspec)%drifts%uBp(Nx+1,1:Nz)=zones(North)%species(Nspec)%drifts%uBp(1,1:Nz)
           zones(k)%species(Nspec)%drifts%uEt(Nx+1,1:Nz)=zones(North)%species(Nspec)%drifts%uEt(1,1:Nz)
           zones(k)%species(Nspec)%drifts%udt(Nx+1,1:Nz)=zones(North)%species(Nspec)%drifts%udt(1,1:Nz)
           zones(k)%species(Nspec)%drifts%uBt(Nx+1,1:Nz)=zones(North)%species(Nspec)%drifts%uBt(1,1:Nz)
           zones(k)%species(Nspec)%drifts%uEps(Nx+1,1:Nz)=zones(North)%species(Nspec)%drifts%uEps(1,1:Nz)
           zones(k)%species(Nspec)%drifts%uBps(Nx+1,1:Nz)=zones(North)%species(Nspec)%drifts%uBps(1,1:Nz)
           zones(k)%species(Nspec)%drifts%uEts(Nx+1,1:Nz)=zones(North)%species(Nspec)%drifts%uEts(1,1:Nz)
           zones(k)%species(Nspec)%drifts%uBts(Nx+1,1:Nz)=zones(North)%species(Nspec)%drifts%uBts(1,1:Nz)
        else
           zones(k)%species(Nspec)%drifts%uEp(Nx+1,1:Nz)=zones(k)%species(Nspec)%drifts%uEp(Nx,1:Nz)
           zones(k)%species(Nspec)%drifts%udp(Nx+1,1:Nz)=zones(k)%species(Nspec)%drifts%udp(Nx,1:Nz)
           zones(k)%species(Nspec)%drifts%uBp(Nx+1,1:Nz)=zones(k)%species(Nspec)%drifts%uBp(Nx,1:Nz)
           zones(k)%species(Nspec)%drifts%uEt(Nx+1,1:Nz)=zones(k)%species(Nspec)%drifts%uEt(Nx,1:Nz)
           zones(k)%species(Nspec)%drifts%udt(Nx+1,1:Nz)=zones(k)%species(Nspec)%drifts%udt(Nx,1:Nz)
           zones(k)%species(Nspec)%drifts%uBt(Nx+1,1:Nz)=zones(k)%species(Nspec)%drifts%uBt(Nx,1:Nz)
           zones(k)%species(Nspec)%drifts%uEps(Nx+1,1:Nz)=zones(k)%species(Nspec)%drifts%uEps(Nx,1:Nz)
           zones(k)%species(Nspec)%drifts%uBps(Nx+1,1:Nz)=zones(k)%species(Nspec)%drifts%uBps(Nx,1:Nz)
           zones(k)%species(Nspec)%drifts%uEts(Nx+1,1:Nz)=zones(k)%species(Nspec)%drifts%uEts(Nx,1:Nz)
           zones(k)%species(Nspec)%drifts%uBts(Nx+1,1:Nz)=zones(k)%species(Nspec)%drifts%uBts(Nx,1:Nz)
        end if
        ! at the bottom
        South=Zones(k)%Neighbors(2)
        mSouth=Zones(k)%MagNeighbors(2)
        if(mSouth.eq.0) then
           ! if there is a neighboring domain, just broadcast
           zones(k)%species(Nspec)%drifts%uEp(0,1:Nz)=zones(South)%species(Nspec)%drifts%uEp(zones(South)%mesh%Nx,1:Nz)
           zones(k)%species(Nspec)%drifts%udp(0,1:Nz)=zones(South)%species(Nspec)%drifts%udp(zones(South)%mesh%Nx,1:Nz)
           zones(k)%species(Nspec)%drifts%uBp(0,1:Nz)=zones(South)%species(Nspec)%drifts%uBp(zones(South)%mesh%Nx,1:Nz)
           zones(k)%species(Nspec)%drifts%uEt(0,1:Nz)=zones(South)%species(Nspec)%drifts%uEt(zones(South)%mesh%Nx,1:Nz)
           zones(k)%species(Nspec)%drifts%udt(0,1:Nz)=zones(South)%species(Nspec)%drifts%udt(zones(South)%mesh%Nx,1:Nz)
           zones(k)%species(Nspec)%drifts%uBt(0,1:Nz)=zones(South)%species(Nspec)%drifts%uBt(zones(South)%mesh%Nx,1:Nz)
           zones(k)%species(Nspec)%drifts%uEps(0,1:Nz)=zones(South)%species(Nspec)%drifts%uEps(zones(South)%mesh%Nx,1:Nz)
           zones(k)%species(Nspec)%drifts%uBps(0,1:Nz)=zones(South)%species(Nspec)%drifts%uBps(zones(South)%mesh%Nx,1:Nz)
           zones(k)%species(Nspec)%drifts%uEts(0,1:Nz)=zones(South)%species(Nspec)%drifts%uEts(zones(South)%mesh%Nx,1:Nz)
           zones(k)%species(Nspec)%drifts%uBts(0,1:Nz)=zones(South)%species(Nspec)%drifts%uBts(zones(South)%mesh%Nx,1:Nz)
        else
           zones(k)%species(Nspec)%drifts%uEp(0,1:Nz)=zones(k)%species(Nspec)%drifts%uEp(1,1:Nz)
           zones(k)%species(Nspec)%drifts%udp(0,1:Nz)=zones(k)%species(Nspec)%drifts%udp(1,1:Nz)
           zones(k)%species(Nspec)%drifts%uBp(0,1:Nz)=zones(k)%species(Nspec)%drifts%uBp(1,1:Nz)
           zones(k)%species(Nspec)%drifts%uEt(0,1:Nz)=zones(k)%species(Nspec)%drifts%uEt(1,1:Nz)
           zones(k)%species(Nspec)%drifts%udt(0,1:Nz)=zones(k)%species(Nspec)%drifts%udt(1,1:Nz)
           zones(k)%species(Nspec)%drifts%uBt(0,1:Nz)=zones(k)%species(Nspec)%drifts%uBt(1,1:Nz)
           zones(k)%species(Nspec)%drifts%uEps(0,1:Nz)=zones(k)%species(Nspec)%drifts%uEps(1,1:Nz)
           zones(k)%species(Nspec)%drifts%uBps(0,1:Nz)=zones(k)%species(Nspec)%drifts%uBps(1,1:Nz)
           zones(k)%species(Nspec)%drifts%uEts(0,1:Nz)=zones(k)%species(Nspec)%drifts%uEts(1,1:Nz)
           zones(k)%species(Nspec)%drifts%uBts(0,1:Nz)=zones(k)%species(Nspec)%drifts%uBts(1,1:Nz)
        end if
        ! on the right
        East=Zones(k)%Neighbors(3)
        mEast=Zones(k)%MagNeighbors(3)
        if(mEast.eq.0) then
           ! if there is a neighboring domain, just broadcast
           zones(k)%species(Nspec)%drifts%uEp(1:Nx,Nz+1)=zones(East)%species(Nspec)%drifts%uEp(1:Nx,1)
           zones(k)%species(Nspec)%drifts%udp(1:Nx,Nz+1)=zones(East)%species(Nspec)%drifts%udp(1:Nx,1)
           zones(k)%species(Nspec)%drifts%uBp(1:Nx,Nz+1)=zones(East)%species(Nspec)%drifts%uBp(1:Nx,1)
           zones(k)%species(Nspec)%drifts%uEt(1:Nx,Nz+1)=zones(East)%species(Nspec)%drifts%uEt(1:Nx,1)
           zones(k)%species(Nspec)%drifts%udt(1:Nx,Nz+1)=zones(East)%species(Nspec)%drifts%udt(1:Nx,1)
           zones(k)%species(Nspec)%drifts%uBt(1:Nx,Nz+1)=zones(East)%species(Nspec)%drifts%uBt(1:Nx,1)
           zones(k)%species(Nspec)%drifts%uEps(1:Nx,Nz+1)=zones(East)%species(Nspec)%drifts%uEps(1:Nx,1)
           zones(k)%species(Nspec)%drifts%uBps(1:Nx,Nz+1)=zones(East)%species(Nspec)%drifts%uBps(1:Nx,1)
           zones(k)%species(Nspec)%drifts%uEts(1:Nx,Nz+1)=zones(East)%species(Nspec)%drifts%uEts(1:Nx,1)
           zones(k)%species(Nspec)%drifts%uBts(1:Nx,Nz+1)=zones(East)%species(Nspec)%drifts%uBts(1:Nx,1)
        else
           zones(k)%species(Nspec)%drifts%uEp(1:Nx,Nz+1)=zones(k)%species(Nspec)%drifts%uEp(1:Nx,Nz)
           zones(k)%species(Nspec)%drifts%udp(1:Nx,Nz+1)=zones(k)%species(Nspec)%drifts%udp(1:Nx,Nz)
           zones(k)%species(Nspec)%drifts%uBp(1:Nx,Nz+1)=zones(k)%species(Nspec)%drifts%uBp(1:Nx,Nz)
           zones(k)%species(Nspec)%drifts%uEt(1:Nx,Nz+1)=zones(k)%species(Nspec)%drifts%uEt(1:Nx,Nz)
           zones(k)%species(Nspec)%drifts%udt(1:Nx,Nz+1)=zones(k)%species(Nspec)%drifts%udt(1:Nx,Nz)
           zones(k)%species(Nspec)%drifts%uBt(1:Nx,Nz+1)=zones(k)%species(Nspec)%drifts%uBt(1:Nx,Nz)
           zones(k)%species(Nspec)%drifts%uEps(1:Nx,Nz+1)=zones(k)%species(Nspec)%drifts%uEps(1:Nx,Nz)
           zones(k)%species(Nspec)%drifts%uBps(1:Nx,Nz+1)=zones(k)%species(Nspec)%drifts%uBps(1:Nx,Nz)
           zones(k)%species(Nspec)%drifts%uEts(1:Nx,Nz+1)=zones(k)%species(Nspec)%drifts%uEts(1:Nx,Nz)
           zones(k)%species(Nspec)%drifts%uBts(1:Nx,Nz+1)=zones(k)%species(Nspec)%drifts%uBts(1:Nx,Nz)
        end if
        ! on the left
        West=Zones(k)%Neighbors(4)
        mWest=Zones(k)%MagNeighbors(4)
        if(mWest.eq.0) then
           ! if there is a neighboring domain, just broadcast
           zones(k)%species(Nspec)%drifts%uEp(1:Nx,0)=zones(West)%species(Nspec)%drifts%uEp(1:Nx,zones(West)%mesh%Nz)
           zones(k)%species(Nspec)%drifts%udp(1:Nx,0)=zones(West)%species(Nspec)%drifts%udp(1:Nx,zones(West)%mesh%Nz)
           zones(k)%species(Nspec)%drifts%uBp(1:Nx,0)=zones(West)%species(Nspec)%drifts%uBp(1:Nx,zones(West)%mesh%Nz)
           zones(k)%species(Nspec)%drifts%uEt(1:Nx,0)=zones(West)%species(Nspec)%drifts%uEt(1:Nx,zones(West)%mesh%Nz)
           zones(k)%species(Nspec)%drifts%udt(1:Nx,0)=zones(West)%species(Nspec)%drifts%udt(1:Nx,zones(West)%mesh%Nz)
           zones(k)%species(Nspec)%drifts%uBt(1:Nx,0)=zones(West)%species(Nspec)%drifts%uBt(1:Nx,zones(West)%mesh%Nz)
           zones(k)%species(Nspec)%drifts%uEps(1:Nx,0)=zones(West)%species(Nspec)%drifts%uEps(1:Nx,zones(West)%mesh%Nz)
           zones(k)%species(Nspec)%drifts%uBps(1:Nx,0)=zones(West)%species(Nspec)%drifts%uBps(1:Nx,zones(West)%mesh%Nz)
           zones(k)%species(Nspec)%drifts%uEts(1:Nx,0)=zones(West)%species(Nspec)%drifts%uEts(1:Nx,zones(West)%mesh%Nz)
           zones(k)%species(Nspec)%drifts%uBts(1:Nx,0)=zones(West)%species(Nspec)%drifts%uBts(1:Nx,zones(West)%mesh%Nz)
        else
           zones(k)%species(Nspec)%drifts%uEp(1:Nx,0)=zones(k)%species(Nspec)%drifts%uEp(1:Nx,1)
           zones(k)%species(Nspec)%drifts%udp(1:Nx,0)=zones(k)%species(Nspec)%drifts%udp(1:Nx,1)
           zones(k)%species(Nspec)%drifts%uBp(1:Nx,0)=zones(k)%species(Nspec)%drifts%uBp(1:Nx,1)
           zones(k)%species(Nspec)%drifts%uEt(1:Nx,0)=zones(k)%species(Nspec)%drifts%uEt(1:Nx,1)
           zones(k)%species(Nspec)%drifts%udt(1:Nx,0)=zones(k)%species(Nspec)%drifts%udt(1:Nx,1)
           zones(k)%species(Nspec)%drifts%uBt(1:Nx,0)=zones(k)%species(Nspec)%drifts%uBt(1:Nx,1)
           zones(k)%species(Nspec)%drifts%uEps(1:Nx,0)=zones(k)%species(Nspec)%drifts%uEps(1:Nx,1)
           zones(k)%species(Nspec)%drifts%uBps(1:Nx,0)=zones(k)%species(Nspec)%drifts%uBps(1:Nx,1)
           zones(k)%species(Nspec)%drifts%uEts(1:Nx,0)=zones(k)%species(Nspec)%drifts%uEts(1:Nx,1)
           zones(k)%species(Nspec)%drifts%uBts(1:Nx,0)=zones(k)%species(Nspec)%drifts%uBts(1:Nx,1)
        end if
     end do
  end do
end subroutine Broadcast_drift3
