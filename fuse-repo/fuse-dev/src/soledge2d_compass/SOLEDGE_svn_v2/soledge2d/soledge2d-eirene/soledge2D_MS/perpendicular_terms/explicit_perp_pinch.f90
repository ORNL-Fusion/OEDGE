module Mpinch_perp

contains

  function explicit_perp_pinch(Field,v_pinch,&
       zone,Nx,Nz) result(fluxes)
    use all_variables, only : reference_parameters
    use MZone
    use Mphysics
    implicit none
    integer*4 :: Nx,Nz
    real*8,intent(inout) :: Field(0:Nx+1,0:Nz+1)
    real*8,intent(in) :: v_pinch(0:Nx+1,0:Nz+1)
    type(TZone),intent(in) :: zone
    real*8 :: fluxes(1:Nx,1:Nz,1:4)
    integer*4 :: i,j
    real*8 :: v_p
    real*8 :: R0,rs0,c0
    fluxes=0. !init to zero for east and west
    R0=reference_parameters%geometry%R0
    rs0=reference_parameters%geometry%rs0
    c0=reference_parameters%fields%c0
    !North
    do i=1,Nx
       do j=1,Nz
          v_p=(v_pinch(i+1,j)+v_pinch(i,j))*0.5
          fluxes(i,j,1)=-(Field(i,j)+Field(i+1,j))*0.5d0*v_p&
               *2.d0*pi*R0/(rs0*c0)
       end do
    end do
    !South
    do i=1,Nx
       do j=1,Nz
          v_p=(v_pinch(i-1,j)+v_pinch(i,j))*0.5
          fluxes(i,j,2)=-(Field(i,j)+Field(i-1,j))*0.5d0*v_p&
               *2.d0*pi*R0/(rs0*c0)! pinch contribution
       end do
    end do
  end function explicit_perp_pinch
  
end module Mpinch_perp
