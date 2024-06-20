subroutine compute_jacobian()
  use all_variables, only : zones, global_parameters
  use Mphysics
  implicit none
  integer*4 :: k
  integer*4 :: Nx,Nz
  real*8,pointer::T(:,:),P(:,:),R(:,:),z(:,:)
  do k=1,global_parameters%N_zones
     Nx=zones(k)%mesh%Nx
     Nz=Zones(k)%mesh%Nz
     allocate(T(0:Nx+1,0:Nz+1))
     allocate(P(0:Nx+1,0:Nz+1))
     allocate(R(0:Nx+1,0:Nz+1))
     allocate(z(0:Nx+1,0:Nz+1))
     T=zones(k)%mesh%z*2.*pi
     P=zones(k)%mesh%x
     R=zones(k)%mesh%Rgeom
     Z=zones(k)%mesh%Zgeom
     ! computation of psi and theta derivatives
     call dpsi_dtheta_centre(Nx,Nz,R,zones(k)%jacobian%dRdP,zones(k)%jacobian%dRdT,P,T)
     call dpsi_dtheta_centre(Nx,Nz,Z,zones(k)%jacobian%dzdP,zones(k)%jacobian%dzdT,P,T)
     deallocate(T,P,R,z)
  end do
  ! broadcast and computation of derivatives on the edges
  call load_Neigh_dpsi()
  call load_Neigh_dtheta()
  call compute_inverse_derivative()
  call MD_broadcast_jacobian()
end subroutine compute_jacobian


subroutine dpsi_dtheta_centre(Nx,Nz,X,dXdP,dXdT,P,T)
  use Mphysics
  implicit none
  integer*4 Nx
  integer*4 Nz
  !   real*8,dimension(0:Nx+1,0:Nz+1),intent(in) :: X
  real*8,dimension(0:Nx+1,0:Nz+1) :: X
  real*8,dimension(0:Nx+1,0:Nz+1),intent(out) :: dXdP
  real*8,dimension(0:Nx+1,0:Nz+1),intent(out) :: dXdT
  real*8,dimension(0:Nx+1,0:Nz+1) :: P
  real*8,dimension(0:Nx+1,0:Nz+1) :: T
  integer*4 i,j
  !second order derivative in the middle
  do i=1,Nx
     do j=1,Nz
        dXdP(i,j)=(X(i+1,j)*(P(i,j)-P(i-1,j))**2.d0&
             -X(i-1,j)*(P(i+1,j)-P(i,j))**2.d0&
             +X(i,j)*((P(i+1,j)-P(i,j))**2.d0-(P(i,j)-P(i-1,j))**2.d0))&
             /((P(i+1,j)-P(i,j))*(P(i,j)-P(i-1,j))*(P(i+1,j)-P(i-1,j)))
        dXdT(i,j)=(X(i,j+1)*(T(i,j)-T(i,j-1))**2.d0&
             -X(i,j-1)*(T(i,j+1)-T(i,j))**2.d0&
             +X(i,j)*((T(i,j+1)-T(i,j))**2.d0&
             -(T(i,j)-T(i,j-1))**2.d0))&
             /((T(i,j+1)-T(i,j))*(T(i,j)-T(i,j-1))&
             *(T(i,j+1)-T(i,j-1)))
!!$        dXdT(i,j)=(X(i,j+1)*modulo((T(i,j)-T(i,j-1)),2.d0*pi)**2.d0&
!!$             -X(i,j-1)*modulo((T(i,j+1)-T(i,j)),2.d0*pi)**2.d0&
!!$             +X(i,j)*(modulo((T(i,j+1)-T(i,j)),2.d0*pi)**2.d0&
!!$             -modulo((T(i,j)-T(i,j-1)),2.d0*pi)**2.d0))&
!!$             /(modulo((T(i,j+1)-T(i,j)),2.d0*pi)*modulo((T(i,j)-T(i,j-1)),2.d0*pi)&
!!$             *modulo((T(i,j+1)-T(i,j-1)),2.d0*pi))
    end do
  end do
end subroutine dpsi_dtheta_centre


subroutine load_Neigh_dpsi()
  use all_variables, only : zones, global_parameters
  implicit none
  integer*4 North,South,East,West
  integer*4 mNorth,mSouth,mEast,mWest
  integer*4 k,Nx,Nz
  integer*4 Vois_info(4)
  do k=1,global_parameters%N_zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     Vois_info(1)=zones(k)%Neighbors(1)
     Vois_info(2)=zones(k)%Neighbors(2)
     Vois_info(3)=zones(k)%Neighbors(3)
     Vois_info(4)=zones(k)%Neighbors(4)
     ! on top
     North=zones(k)%Neighbors(1)
     mNorth=zones(k)%MagNeighbors(1)
     if(mNorth.eq.0) then
        ! if there is a neighboring domain, just broadcast
        zones(k)%jacobian%dRdP(Nx+1,1:Nz)=zones(North)%jacobian%dRdP(1,1:Nz)
        zones(k)%jacobian%dzdP(Nx+1,1:Nz)=zones(North)%jacobian%dzdP(1,1:Nz)
     else
        call compute_Neigh_dpsi(Nx,Nz,zones(k)%mesh%Rgeom,zones(k)%mesh%x,zones(k)%jacobian%dRdP,1,Vois_info,zones(k)%mesh%cornerX,zones(k)%mesh%cornerRg)
        call compute_Neigh_dpsi(Nx,Nz,zones(k)%mesh%Zgeom,zones(k)%mesh%x,zones(k)%jacobian%dzdP,1,Vois_info,zones(k)%mesh%cornerX,zones(k)%mesh%cornerZg)
     end if
     ! at the bottom
     South=Zones(k)%Neighbors(2)
     mSouth=Zones(k)%MagNeighbors(2)
     if(mSouth.eq.0) then
        ! if there is a neighboring domain, just broadcast
        zones(k)%jacobian%dRdP(0,1:Nz)=zones(South)%jacobian%dRdP(zones(South)%mesh%Nx,1:Nz)
        zones(k)%jacobian%dzdP(0,1:Nz)=zones(South)%jacobian%dzdP(zones(South)%mesh%Nx,1:Nz)
     else
        call compute_Neigh_dpsi(Nx,Nz,zones(k)%mesh%Rgeom,zones(k)%mesh%x,zones(k)%jacobian%dRdP,2,Vois_info,zones(k)%mesh%cornerX,zones(k)%mesh%cornerRg)
        call compute_Neigh_dpsi(Nx,Nz,zones(k)%mesh%Zgeom,zones(k)%mesh%x,zones(k)%jacobian%dzdP,2,Vois_info,zones(k)%mesh%cornerX,zones(k)%mesh%cornerZg)
     end if
     ! on the rigth
     East=Zones(k)%Neighbors(3)
     mEast=Zones(k)%MagNeighbors(3)
     if(mEast.eq.0) then
        ! if there is a neighboring domain, just broadcast
        zones(k)%jacobian%dRdP(1:Nx,Nz+1)=zones(East)%jacobian%dRdP(1:Nx,1)
        zones(k)%jacobian%dzdP(1:Nx,Nz+1)=zones(East)%jacobian%dzdP(1:Nx,1)
     else
        call compute_Neigh_dpsi(Nx,Nz,zones(k)%mesh%Rgeom,zones(k)%mesh%x,zones(k)%jacobian%dRdP,3,Vois_info,zones(k)%mesh%cornerX,zones(k)%mesh%cornerRg)
        call compute_Neigh_dpsi(Nx,Nz,zones(k)%mesh%Zgeom,zones(k)%mesh%x,zones(k)%jacobian%dzdP,3,Vois_info,zones(k)%mesh%cornerX,zones(k)%mesh%cornerZg)
     end if
     ! on the left
     West=Zones(k)%Neighbors(4)
     mWest=Zones(k)%MagNeighbors(4)
     if(mWest.eq.0) then
        ! if there is a neighboring domain, just broadcast
        zones(k)%jacobian%dRdP(1:Nx,0)=zones(West)%jacobian%dRdP(1:Nx,Zones(West)%mesh%Nz)
        zones(k)%jacobian%dzdP(1:Nx,0)=zones(West)%jacobian%dzdP(1:Nx,Zones(West)%mesh%Nz)
     else
        call compute_Neigh_dpsi(Nx,Nz,zones(k)%mesh%Rgeom,zones(k)%mesh%x,zones(k)%jacobian%dRdP,4,Vois_info,zones(k)%mesh%cornerX,zones(k)%mesh%cornerRg)
        call compute_Neigh_dpsi(Nx,Nz,zones(k)%mesh%Zgeom,zones(k)%mesh%x,zones(k)%jacobian%dzdP,4,Vois_info,zones(k)%mesh%cornerX,zones(k)%mesh%cornerZg)
     end if
  end do
end subroutine load_Neigh_dpsi


! ######################################################################################### !
! #####   Special derivative in psi* on the edge (Boundary conditions)             ######## !
! #####   Second order decentered                                                  ######## !
! ######################################################################################### !
subroutine compute_Neigh_dpsi(Nx,Nz,X,P,dXdP,Neigh,Vois_info,cornerP,cornerX)
  implicit none
  integer*4 Nx,Nz
  real*8,dimension(0:Nx+1,0:Nz+1),intent(inout) :: X
  real*8,dimension(0:Nx+1,0:Nz+1),intent(inout) :: P
  real*8,dimension(0:Nx+1,0:Nz+1),intent(out) :: dXdP
  real*8,intent(in) :: cornerX(2,2,2)
  real*8,intent(in) :: cornerP(2,2,2)
  integer*4 Neigh !1N, 2S, 3E, 4W
  integer*4 Vois_info(4)
  select case(Neigh)
  case(1)        ! on top
     dXdP(Nx+1,1:Nz)=(X(Nx+1,1:Nz)-X(Nx,1:Nz))/(P(Nx+1,1:Nz)-P(Nx,1:Nz))
  case(2)       ! at the bottom
     dXdP(0,1:Nz)=(X(1,1:Nz)-X(0,1:Nz))/(P(1,1:Nz)-P(0,1:Nz))
  case(3)       ! on the right
     dXdP(2:Nx-1,Nz+1)=(X(3:Nx,Nz+1)*(P(2:Nx-1,Nz+1)-P(1:Nx-2,Nz+1))**2.d0&
          -X(1:Nx-2,Nz+1)*(P(3:Nx,Nz+1)-P(2:Nx-1,Nz+1))**2.d0&
          +X(2:Nx-1,Nz+1)*((P(3:Nx,Nz+1)-P(2:Nx-1,Nz+1))**2.d0-(P(2:Nx-1,Nz+1)-P(1:Nx-2,Nz+1))**2.d0))&
          /((P(3:Nx,Nz+1)-P(2:Nx-1,Nz+1))*(P(2:Nx-1,Nz+1)-P(1:Nx-2,Nz+1))*(P(3:Nx,Nz+1)-P(1:Nx-2,Nz+1)))
     if(Vois_info(2).gt.0) then !use neighbor info
        X(0,Nz+1)=cornerX(2,2,1) !SE
        P(0,Nz+1)=cornerP(2,2,1) !SE
        dXdP(1,Nz+1)=(X(2,Nz+1)*(P(1,Nz+1)-P(0,Nz+1))**2.d0&
             -X(0,Nz+1)*(P(2,Nz+1)-P(1,Nz+1))**2.d0&
             +X(1,Nz+1)*((P(2,Nz+1)-P(1,Nz+1))**2.d0-(P(1,Nz+1)-P(0,Nz+1))**2.d0))&
             /((P(2,Nz+1)-P(1,Nz+1))*(P(1,Nz+1)-P(0,Nz+1))*(P(2,Nz+1)-P(0,Nz+1)))
     else !corner decenter derivative
        dXdP(1,Nz+1)=(X(1,Nz+1)-X(2,Nz+1))/(P(1,Nz+1)-P(2,Nz+1))
     end if
     if(Vois_info(1).gt.0) then
        X(Nx+1,Nz+1)=cornerX(1,2,1) !NE
        P(Nx+1,Nz+1)=cornerP(1,2,1) !NE
        dXdP(Nx,Nz+1)=(X(Nx+1,Nz+1)*(P(Nx,Nz+1)-P(Nx-1,Nz+1))**2.d0&
             -X(Nx-1,Nz+1)*(P(Nx+1,Nz+1)-P(Nx,Nz+1))**2.d0&
             +X(Nx,Nz+1)*((P(Nx+1,Nz+1)-P(Nx,Nz+1))**2.d0-(P(Nx,Nz+1)-P(Nx-1,Nz+1))**2.d0))&
             /((P(Nx+1,Nz+1)-P(Nx,Nz+1))*(P(Nx,Nz+1)-P(Nx-1,Nz+1))*(P(Nx+1,Nz+1)-P(Nx-1,Nz+1)))
     else
        dXdP(Nx,Nz+1)=(X(Nx,Nz+1)-X(Nx-1,Nz+1))/(P(Nx,Nz+1)-P(Nx-1,Nz+1))
     end if
  case(4)      ! on the left
     dXdP(2:Nx-1,0)=(X(3:Nx,0)*(P(2:Nx-1,0)-P(1:Nx-2,0))**2.d0&
          -X(1:Nx-2,0)*(P(3:Nx,0)-P(2:Nx-1,0))**2.d0&
          +X(2:Nx-1,0)*((P(3:Nx,0)-P(2:Nx-1,0))**2.d0-(P(2:Nx-1,0)-P(1:Nx-2,0))**2.d0))&
          /((P(3:Nx,0)-P(2:Nx-1,0))*(P(2:Nx-1,0)-P(1:Nx-2,0))*(P(3:Nx,0)-P(1:Nx-2,0)))
     if(Vois_info(2).gt.0) then !use neighbor info
        X(0,0)=cornerX(2,1,1) !SW
        P(0,0)=cornerP(2,1,1) !SW
        dXdP(1,0)=(X(2,0)*(P(1,0)-P(0,0))**2.d0&
             -X(0,0)*(P(2,0)-P(1,0))**2.d0&
             +X(1,0)*((P(2,0)-P(1,0))**2.d0-(P(1,0)-P(0,0))**2.d0))&
             /((P(2,0)-P(1,0))*(P(1,0)-P(0,0))*(P(2,0)-P(0,0)))
     else !corner decenter derivative
        dXdP(1,0)=(X(1,0)-X(2,0))/(P(1,0)-P(2,0))
     end if
     if(Vois_info(1).gt.0) then
        X(Nx+1,0)=cornerX(1,1,1) !NW
        P(Nx+1,0)=cornerP(1,1,1) !NW
        dXdP(Nx,0)=(X(Nx+1,0)*(P(Nx,0)-P(Nx-1,0))**2.d0&
             -X(Nx-1,0)*(P(Nx+1,0)-P(Nx,0))**2.d0&
             +X(Nx,0)*((P(Nx+1,0)-P(Nx,0))**2.d0-(P(Nx,0)-P(Nx-1,0))**2.d0))&
             /((P(Nx+1,0)-P(Nx,0))*(P(Nx,0)-P(Nx-1,0))*(P(Nx+1,0)-P(Nx-1,0)))
     else
        dXdP(Nx,0)=(X(Nx,0)-X(Nx-1,0))/(P(Nx,0)-P(Nx-1,0))
     end if
  case default
  end select
end subroutine compute_Neigh_dpsi


subroutine load_Neigh_dtheta()
  use all_variables, only : zones, global_parameters
  use Mphysics
  implicit none
  integer*4 North,South,East,West
  integer*4 mNorth,mSouth,mEast,mWest
  integer*4 k,Nx,Nz
  integer*4 Vois_info(4)
  real*8,allocatable :: theta(:,:)
  ! for all domains
  do k=1,global_parameters%N_zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     allocate(theta(0:Nx+1,0:Nz+1))
     theta=zones(k)%mesh%z*2.*pi
     Vois_info(1)=Zones(k)%Neighbors(1)
     Vois_info(2)=Zones(k)%Neighbors(2)
     Vois_info(3)=Zones(k)%Neighbors(3)
     Vois_info(4)=Zones(k)%Neighbors(4)
     ! on top
     North=Zones(k)%Neighbors(1)
     mNorth=Zones(k)%MagNeighbors(1)
     if(mNorth.eq.0) then
        ! if there is a neighboring domain, just broadcast
        zones(k)%jacobian%dRdT(Nx+1,1:Nz)=zones(North)%jacobian%dRdT(1,1:Nz)
        zones(k)%jacobian%dzdT(Nx+1,1:Nz)=zones(North)%jacobian%dzdT(1,1:Nz)
     else
        call compute_Neigh_dtheta(Nx,Nz,zones(k)%mesh%Rgeom,theta,zones(k)%jacobian%dRdT,1,Vois_info,zones(k)%mesh%cornerZ*2._8*pi,zones(k)%mesh%cornerRg)
        call compute_Neigh_dtheta(Nx,Nz,zones(k)%mesh%Zgeom,theta,zones(k)%jacobian%dzdT,1,Vois_info,zones(k)%mesh%cornerZ*2._8*pi,zones(k)%mesh%cornerZg)
     end if
     ! at the bottom
     South=Zones(k)%Neighbors(2)
     mSouth=Zones(k)%MagNeighbors(2)
     if(mSouth.eq.0) then
        ! if there is a neighboring domain, just broadcast
        zones(k)%jacobian%dRdT(0,1:Nz)=zones(South)%jacobian%dRdT(Zones(South)%mesh%Nx,1:Nz)
        zones(k)%jacobian%dzdT(0,1:Nz)=zones(South)%jacobian%dzdT(Zones(South)%mesh%Nx,1:Nz)
     else
        call compute_Neigh_dtheta(Nx,Nz,zones(k)%mesh%Rgeom,theta,zones(k)%jacobian%dRdT,2,Vois_info,zones(k)%mesh%cornerZ*2._8*pi,zones(k)%mesh%cornerRg)
        call compute_Neigh_dtheta(Nx,Nz,zones(k)%mesh%Zgeom,theta,zones(k)%jacobian%dzdT,2,Vois_info,zones(k)%mesh%cornerZ*2._8*pi,zones(k)%mesh%cornerZg)
     end if
     ! on the right
     East=Zones(k)%Neighbors(3)
     mEast=Zones(k)%MagNeighbors(3)
     if(mEast.eq.0) then
        ! if there is a neighboring domain, just broadcast
        zones(k)%jacobian%dRdT(1:Nx,Nz+1)=zones(East)%jacobian%dRdT(1:Nx,1)
        zones(k)%jacobian%dzdT(1:Nx,Nz+1)=zones(East)%jacobian%dzdT(1:Nx,1)
     else
        call compute_Neigh_dtheta(Nx,Nz,zones(k)%mesh%Rgeom,theta,zones(k)%jacobian%dRdT,3,Vois_info,zones(k)%mesh%cornerZ*2._8*pi,zones(k)%mesh%cornerRg)
        call compute_Neigh_dtheta(Nx,Nz,zones(k)%mesh%Zgeom,theta,zones(k)%jacobian%dzdT,3,Vois_info,zones(k)%mesh%cornerZ*2._8*pi,zones(k)%mesh%cornerZg)
     end if
     ! on the left
     West=Zones(k)%Neighbors(4)
     mWest=Zones(k)%MagNeighbors(4)
     if(mWest.eq.0) then
        ! if there is a neighboring domain, just broadcast
        zones(k)%jacobian%dRdT(1:Nx,0)=zones(West)%jacobian%dRdT(1:Nx,Zones(West)%mesh%Nz)
        zones(k)%jacobian%dzdT(1:Nx,0)=zones(West)%jacobian%dzdT(1:Nx,Zones(West)%mesh%Nz)
     else
        call compute_Neigh_dtheta(Nx,Nz,zones(k)%mesh%Rgeom,theta,zones(k)%jacobian%dRdT,4,Vois_info,zones(k)%mesh%cornerZ*2._8*pi,zones(k)%mesh%cornerRg)
        call compute_Neigh_dtheta(Nx,Nz,zones(k)%mesh%Zgeom,theta,zones(k)%jacobian%dzdT,4,Vois_info,zones(k)%mesh%cornerZ*2._8*pi,zones(k)%mesh%cornerZg)
     end if
     deallocate(theta)
  end do
end subroutine load_Neigh_dtheta


  ! ######################################################################################### !
  ! #####   Special derivative in theta* on the edge (Boundary conditions)           ######## !
  ! #####   Second order decentered                                                  ######## !
  ! ######################################################################################### !
  subroutine compute_Neigh_dtheta(Nx,Nz,X,T,dXdT,Neigh,Vois_info,cornerT,cornerX)
    implicit none
    integer*4 Nx,Nz
    real*8,dimension(0:Nx+1,0:Nz+1),intent(inout) :: X
    real*8,dimension(0:Nx+1,0:Nz+1),intent(inout) :: T
    real*8,dimension(0:Nx+1,0:Nz+1),intent(out) :: dXdT
    real*8,intent(in) :: cornerT(2,2,2)
    real*8,intent(in) :: cornerX(2,2,2)
    integer*4 Neigh !1N, 2S, 3E, 4W
    integer*4 Vois_info(4)
    real*8 pi
    pi=4._8*atan(1._8)
    select case(Neigh)
    case(1) ! on top
       dXdT(Nx+1,2:Nz-1)=(X(Nx+1,3:Nz)*(T(Nx+1,2:Nz-1)-T(Nx+1,1:Nz-2))**2.d0&
            -X(Nx+1,1:Nz-2)*(T(Nx+1,3:Nz)-T(Nx+1,2:Nz-1))**2.d0&
            +X(Nx+1,2:Nz-1)*((T(Nx+1,3:Nz)-T(Nx+1,2:Nz-1))**2.d0&
            -(T(Nx+1,2:Nz-1)-T(Nx+1,1:Nz-2))**2.d0))&
            /((T(Nx+1,3:Nz)-T(Nx+1,2:Nz-1))*(T(Nx+1,2:Nz-1)-T(Nx+1,1:Nz-2))&
            *(T(Nx+1,3:Nz)-T(Nx+1,1:Nz-2)))
       if(Vois_info(4).gt.0) then
          X(Nx+1,0)=cornerX(1,1,2) !WN
          T(Nx+1,0)=cornerT(1,1,2) !WN
          dXdT(Nx+1,1)=(X(Nx+1,2)*(T(Nx+1,1)-T(Nx+1,0))**2.d0&
               -X(Nx+1,0)*(T(Nx+1,2)-T(Nx+1,1))**2.d0&
               +X(Nx+1,1)*((T(Nx+1,2)-T(Nx+1,1))**2.d0&
               -(T(Nx+1,1)-T(Nx+1,0))**2.d0))&
               /((T(Nx+1,2)-T(Nx+1,1))*(T(Nx+1,1)-T(Nx+1,0))&
               *(T(Nx+1,2)-T(Nx+1,0)))
       else
          dXdT(Nx+1,1)=(X(Nx+1,2)-X(Nx+1,1))/((T(Nx+1,2)-T(Nx+1,1)))
       end if
       if(Vois_info(3).gt.0) then
          X(Nx+1,Nz+1)=cornerX(1,2,2) !EN
          T(Nx+1,Nz+1)=cornerT(1,2,2) !EN
          dXdT(Nx+1,Nz)=(X(Nx+1,Nz+1)*(T(Nx+1,Nz)-T(Nx+1,Nz-1))**2.d0&
               -X(Nx+1,Nz-1)*(T(Nx+1,Nz+1)-T(Nx+1,Nz))**2.d0&
               +X(Nx+1,Nz)*((T(Nx+1,Nz+1)-T(Nx+1,Nz))**2.d0&
               -(T(Nx+1,Nz)-T(Nx+1,Nz-1))**2.d0))&
               /((T(Nx+1,Nz+1)-T(Nx+1,Nz))*(T(Nx+1,Nz)-T(Nx+1,Nz-1))&
               *(T(Nx+1,Nz+1)-T(Nx+1,Nz-1)))
       else
          dXdT(Nx+1,Nz)=(X(Nx+1,Nz)-X(Nx+1,Nz-1))/((T(Nx+1,Nz)-T(Nx+1,Nz-1)))
       end if
    case(2) ! at the bottom
       dXdT(0,2:Nz-1)=(X(0,3:Nz)*(T(0,2:Nz-1)-T(0,1:Nz-2))**2.d0&
            -X(0,1:Nz-2)*(T(0,3:Nz)-T(0,2:Nz-1))**2.d0&
            +X(0,2:Nz-1)*((T(0,3:Nz)-T(0,2:Nz-1))**2.d0&
            -(T(0,2:Nz-1)-T(0,1:Nz-2))**2.d0))&
            /((T(0,3:Nz)-T(0,2:Nz-1))*(T(0,2:Nz-1)-T(0,1:Nz-2))&
            *(T(0,3:Nz)-T(0,1:Nz-2)))
       if(Vois_info(4).gt.0) then
          X(0,0)=cornerX(2,1,2) !WS
          T(0,0)=cornerT(2,1,2) !WS
          dXdT(0,1)=(X(0,2)*(T(0,1)-T(0,0))**2.d0&
               -X(0,0)*(T(0,2)-T(0,1))**2.d0&
               +X(0,1)*((T(0,2)-T(0,1))**2.d0&
               -(T(0,1)-T(0,0))**2.d0))&
               /((T(0,2)-T(0,1))*(T(0,1)-T(0,0))&
               *(T(0,2)-T(0,0)))
       else
          dXdT(0,1)=(X(0,2)-X(0,1))/((T(0,2)-T(0,1)))
       end if
       if(Vois_info(3).gt.0) then
          X(0,Nz+1)=cornerX(2,2,2) !ES
          T(0,Nz+1)=cornerT(2,2,2) !ES
          dXdT(0,Nz)=(X(0,Nz+1)*(T(0,Nz)-T(0,Nz-1))**2.d0&
               -X(0,Nz-1)*(T(0,Nz+1)-T(0,Nz))**2.d0&
               +X(0,Nz)*((T(0,Nz+1)-T(0,Nz))**2.d0&
               -(T(0,Nz)-T(0,Nz-1))**2.d0))&
               /((T(0,Nz+1)-T(0,Nz))*(T(0,Nz)-T(0,Nz-1))&
               *(T(0,Nz+1)-T(0,Nz-1)))
       else
          dXdT(0,Nz)=(X(0,Nz)-X(0,Nz-1))/((T(0,Nz)-T(0,Nz-1)))
       end if
    case(3) ! on the right
       dXdT(1:Nx,Nz+1)=(X(1:Nx,Nz+1)-X(1:Nx,Nz))/(T(1:Nx,Nz+1)-T(1:Nx,Nz))
    case(4) ! on the left
       dXdT(1:Nx,0)=(X(1:Nx,1)-X(1:Nx,0))/(T(1:Nx,1)-T(1:Nx,0))
    case default
    end select
  end subroutine compute_Neigh_dtheta


  ! ######################################################################################### !
  ! ######     Computing the jacobian and inverse from derivatives             ############## !
  ! ######################################################################################### !
  subroutine compute_inverse_derivative()
    use all_variables, only : zones, global_parameters
    implicit none
    integer*4 i,j,k
    integer*4 Nx,Nz
    real*8 det
    do k=1,global_parameters%N_Zones
       Nx=zones(k)%mesh%Nx
       Nz=zones(k)%mesh%Nz
       do i=0,Nx+1
          do j=0,Nz+1
             det=zones(k)%jacobian%dRdT(i,j)*zones(k)%jacobian%dzdP(i,j)-zones(k)%jacobian%dzdT(i,j)*zones(k)%jacobian%dRdP(i,j)
             if(det.ne.0._8) then
                zones(k)%jacobian%dTdR(i,j)=zones(k)%jacobian%dzdP(i,j)/det
                zones(k)%jacobian%dTdz(i,j)=-zones(k)%jacobian%dRdP(i,j)/det
                zones(k)%jacobian%dPdR(i,j)=-zones(k)%jacobian%dzdT(i,j)/det
                zones(k)%jacobian%dPdz(i,j)=zones(k)%jacobian%dRdT(i,j)/det             
             else
                zones(k)%jacobian%dTdR(i,j)=0._8
                zones(k)%jacobian%dTdz(i,j)=0._8
                zones(k)%jacobian%dPdR(i,j)=0._8
                zones(k)%jacobian%dPdz(i,j)=0._8
             end if
          end do
       end do
    end do
  end subroutine compute_inverse_derivative
