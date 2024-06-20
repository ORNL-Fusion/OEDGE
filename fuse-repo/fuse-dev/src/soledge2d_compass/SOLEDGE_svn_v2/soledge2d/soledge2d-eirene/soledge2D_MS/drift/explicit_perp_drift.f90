module MDrift_perp

contains

  function explicit_perp_drift(Field,Driftp,Driftt,&
       zone,Nx,Nz) result(fluxes)
    use all_variables, only : reference_parameters
    use MZone
    use Mphysics
    implicit none
    integer*4 :: Nx,Nz
    real*8,intent(inout) :: Field(0:Nx+1,0:Nz+1)
    real*8,intent(in) :: Driftp(0:Nx+1,0:Nz+1), Driftt(0:Nx+1,0:Nz+1)
    type(TZone),intent(in) :: zone
    real*8 :: fluxes(1:Nx,1:Nz,1:4)
    integer*4 :: i,j
    real*8 :: x(0:Nx+1,0:Nz+1),z(0:Nx+1,0:Nz+1)
    real*8 :: c_pp(0:Nx+1,0:Nz+1),c_tt(0:Nx+1,0:Nz+1)
    real*8 :: v_d
    real*8 :: R0,c0
    R0=reference_parameters%geometry%R0
    c0=reference_parameters%fields%c0
    x=zone%mesh%x
    z=zone%mesh%z
    c_pp=zone%metric_coefficients%c_pp
    c_tt=zone%metric_coefficients%c_tt
    !North
    do i=1,Nx
       do j=1,Nz
          v_d=(Driftp(i+1,j)+Driftp(i,j))*0.5d0
          fluxes(i,j,1)=-0.5d0*v_d*(Field(i,j)+Field(i+1,j))
       end do
    end do
    !South
    do i=1,Nx
       do j=1,Nz
          v_d=(Driftp(i-1,j)+Driftp(i,j))*0.5d0
          fluxes(i,j,2)=-0.5d0*v_d*(Field(i,j)+Field(i-1,j))
       end do
    end do
    !East
    do i=1,Nx
       do j=1,Nz
          v_d=0.5d0*(Driftt(i,j+1)+Driftt(i,j))
          fluxes(i,j,3)=-0.5d0*v_d*(Field(i,j+1)+Field(i,j))
       end do
    end do
    !West
    do i=1,Nx
       do j=1,Nz
          v_d=0.5d0*(Driftt(i,j-1)+Driftt(i,j))
          fluxes(i,j,4)=-0.5d0*v_d*(Field(i,j-1)+Field(i,j))
       end do
    end do
  end function explicit_perp_drift

end module MDrift_perp
