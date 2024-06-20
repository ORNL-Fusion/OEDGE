module MDiffusion_perp

contains

  function explicit_perp_diffusion(Field,Corners,Diffusivity_p,&
       zone,Nx,Nz) result(fluxes)
    use all_variables, only : reference_parameters
    use MZone
    implicit none
    integer*4 :: Nx,Nz
    real*8,intent(inout) :: Field(0:Nx+1,0:Nz+1)
    real*8,intent(in) :: Corners(2,2,2)
    real*8,intent(in) :: Diffusivity_p(0:Nx+1,0:Nz+1)
    type(TZone),intent(in) :: zone
    real*8 :: fluxes(1:Nx,1:Nz,1:4)
    integer*4 :: i,j
    real*8 :: x(0:Nx+1,0:Nz+1),z(0:Nx+1,0:Nz+1)
    real*8 :: cpp(0:Nx+1,0:Nz+1),cpt(0:Nx+1,0:Nz+1),ctt(0:Nx+1,0:Nz+1)
    real*8 :: D_p
    x=zone%mesh%x
    z=zone%mesh%z
    cpp=zone%metric_coefficients%cpp
    cpt=zone%metric_coefficients%cpt
    ctt=zone%metric_coefficients%ctt
    !North
    Field(Nx+1,0)=(Corners(1,1,1)+Corners(1,1,2))*0.5d0 !NW corner
    x(Nx+1,0)=zone%mesh%cornerx(1,1,1)
    z(Nx+1,0)=zone%mesh%cornerz(1,1,1)
    Field(Nx+1,Nz+1)=(Corners(1,2,1)+Corners(1,2,2))*0.5d0 !NE corner
    x(Nx+1,Nz+1)=zone%mesh%cornerx(1,2,1)
    z(Nx+1,Nz+1)=zone%mesh%cornerz(1,2,1)
    do i=1,Nx
       do j=1,Nz
          D_p=(Diffusivity_p(i+1,j)+Diffusivity_p(i,j))*0.5
          fluxes(i,j,1)=D_p&
               *((Field(i,j)-Field(i+1,j))/(x(i,j)-x(i+1,j))&
               *0.5d0*(cpp(i,j)+cpp(i+1,j))&
               +0.5d0*((Field(i,j+1)-Field(i,j-1))/(z(i,j+1)-z(i,j-1))&
               *cpt(i,j)&
               +(Field(i+1,j+1)-Field(i+1,j-1))/(z(i+1,j+1)-z(i+1,j-1))&
               *cpt(i+1,j))&
               )/(0.5D0*(sqrt(cpp(i,j))+sqrt(cpp(i+1,j))))
       end do
    end do
    !South
    Field(0,0)=(Corners(2,1,1)+Corners(2,1,2))*0.5d0 !SW corner
    x(0,0)=zone%mesh%cornerx(2,1,1)
    z(0,0)=zone%mesh%cornerz(2,1,1)
    Field(0,Nz+1)=(Corners(2,2,1)+Corners(2,2,2))*0.5d0 !SE corner
    x(0,Nz+1)=zone%mesh%cornerx(2,2,1)
    z(0,Nz+1)=zone%mesh%cornerz(2,2,1)
    do i=1,Nx
       do j=1,Nz
          D_p=(Diffusivity_p(i-1,j)+Diffusivity_p(i,j))*0.5
          fluxes(i,j,2)=D_p&
               *((Field(i,j)-Field(i-1,j))/(x(i,j)-x(i-1,j))&
               *0.5d0*(cpp(i,j)+cpp(i-1,j))&
               +0.5d0*((Field(i,j+1)-Field(i,j-1))/(z(i,j+1)-z(i,j-1))&
               *cpt(i,j)&
               +(Field(i-1,j+1)-Field(i-1,j-1))/(z(i-1,j+1)-z(i-1,j-1))&
               *cpt(i-1,j))&
               )/(0.5D0*(sqrt(cpp(i,j))+sqrt(cpp(i-1,j))))
       end do
    end do
    !East
    Field(Nx+1,Nz+1)=(Corners(1,2,2)+Corners(1,2,1))*0.5d0 !EN corner
    x(Nx+1,Nz+1)=zone%mesh%cornerx(1,2,2)
    z(Nx+1,Nz+1)=zone%mesh%cornerz(1,2,2)
    Field(0,Nz+1)=(Corners(2,2,2)+Corners(2,2,1))*0.5d0 !ES corner
    x(0,Nz+1)=zone%mesh%cornerx(2,2,2)
    z(0,Nz+1)=zone%mesh%cornerz(2,2,2)
    do i=1,Nx
       do j=1,Nz
          D_p=0.5d0*(Diffusivity_p(i,j+1)+Diffusivity_p(i,j))
          fluxes(i,j,3)=D_p&
               *(0.5d0*((Field(i+1,j)-Field(i-1,j))&
               /(x(i+1,j)-x(i-1,j))*cpt(i,j)&
               +(Field(i+1,j+1)-Field(i-1,j+1))&
               /(x(i+1,j+1)-x(i-1,j+1))*cpt(i,j+1)))&
               /(0.5d0*(sqrt(ctt(i,j))+sqrt(ctt(i,j+1))))
       end do
    end do
    !West
    Field(Nx+1,0)=(Corners(1,1,2)+Corners(1,1,1))*0.5d0 !WN corner
    x(Nx+1,0)=zone%mesh%cornerx(1,1,2)
    z(Nx+1,0)=zone%mesh%cornerz(1,1,2)
    Field(0,0)=(Corners(2,1,2)+Corners(2,1,1))*0.5d0 !WS corner
    x(0,0)=zone%mesh%cornerx(2,1,2)
    z(0,0)=zone%mesh%cornerz(2,1,2)
    do i=1,Nx
       do j=1,Nz
          D_p=0.5d0*(Diffusivity_p(i,j-1)+Diffusivity_p(i,j))
          fluxes(i,j,4)=D_p&
               *(0.5d0*((Field(i+1,j)-Field(i-1,j))&
               /(x(i+1,j)-x(i-1,j))*cpt(i,j)&
               +(Field(i+1,j-1)-Field(i-1,j-1))&
               /(x(i+1,j-1)-x(i-1,j-1))*cpt(i,j-1)))&
               /(0.5d0*(sqrt(ctt(i,j))+sqrt(ctt(i,j-1))))
       end do
    end do
  end function explicit_perp_diffusion


  function in_flux_surface_perp_diffusion(Field,Corners,D_p,D_t,&
       zone,Nx,Nz) result(fluxes)
    use all_variables, only : reference_parameters
    use MZone
    implicit none
    integer*4 :: Nx,Nz
    real*8,intent(inout) :: Field(0:Nx+1,0:Nz+1)
    real*8,intent(in) :: Corners(2,2,2)
    real*8,intent(in) :: D_p(0:Nx+1,0:Nz+1)
    real*8,intent(in) :: D_t(0:Nx+1,0:Nz+1)
    type(TZone),intent(in) :: zone
    real*8 :: fluxes(1:Nx,1:Nz,1:4)
    integer*4 :: i,j
    real*8 :: x(0:Nx+1,0:Nz+1),z(0:Nx+1,0:Nz+1)
    real*8 :: cpp(0:Nx+1,0:Nz+1),cpt(0:Nx+1,0:Nz+1),ctt(0:Nx+1,0:Nz+1),G(0:Nx+1,0:Nz+1)
    real*8 :: A
    real*8 :: coef_east, coef_west
    A=reference_parameters%geometry%A
    x=zone%mesh%x
    z=zone%mesh%z
    cpp=zone%metric_coefficients%cpp
    cpt=zone%metric_coefficients%cpt
    ctt=zone%metric_coefficients%ctt
    G=zone%metric_coefficients%G
    !North
    fluxes(:,:,1)=0.D0
    !South
    fluxes(:,:,2)=0.D0
    !East
    Field(Nx+1,Nz+1)=(Corners(1,2,2)+Corners(1,2,1))*0.5d0 !EN corner
    x(Nx+1,Nz+1)=zone%mesh%cornerx(1,2,2)
    z(Nx+1,Nz+1)=zone%mesh%cornerz(1,2,2)
    Field(0,Nz+1)=(Corners(2,2,2)+Corners(2,2,1))*0.5d0 !ES corner
    x(0,Nz+1)=zone%mesh%cornerx(2,2,2)
    z(0,Nz+1)=zone%mesh%cornerz(2,2,2)
    do i=1,Nx
       do j=1,Nz
          coef_east=((D_t(i,j)*(ctt(i,j)-1.d0/A**2.d0*G(i,j)**2.d0)&
               +(D_p(i,j)-D_t(i,j))*cpt(i,j)**2.d0/cpp(i,j))&
               +(D_t(i,j+1)*(ctt(i,j+1)-1.d0/A**2.d0*G(i,j+1)**2.d0)&
               +(D_p(i,j+1)-D_t(i,j+1))*cpt(i,j+1)**2.d0/cpp(i,j+1)))*0.5d0&
               /(0.5d0*(sqrt(ctt(i,j))+sqrt(ctt(i,j+1))))
          fluxes(i,j,3)=coef_east*(Field(i,j+1)-Field(i,j))&
               /(z(i,j+1)-z(i,j))
       end do
    end do
    !West
    Field(Nx+1,0)=(Corners(1,1,2)+Corners(1,1,1))*0.5d0 !WN corner
    x(Nx+1,0)=zone%mesh%cornerx(1,1,2)
    z(Nx+1,0)=zone%mesh%cornerz(1,1,2)
    Field(0,0)=(Corners(2,1,2)+Corners(2,1,1))*0.5d0 !WS corner
    x(0,0)=zone%mesh%cornerx(2,1,2)
    z(0,0)=zone%mesh%cornerz(2,1,2)
    do i=1,Nx
       do j=1,Nz
          coef_west = &
               ((D_t(i,j)*(ctt(i,j)-1.d0/A**2.d0*G(i,j)**2.d0)&
               +(D_p(i,j)-D_t(i,j))*cpt(i,j)**2.d0/cpp(i,j))&
               +(D_t(i,j-1)*(ctt(i,j-1)-1.d0/A**2.d0*G(i,j-1)**2.d0)&
               +(D_p(i,j-1)-D_t(i,j-1))*cpt(i,j-1)**2.d0/cpp(i,j-1)))*0.5d0&
               /(0.5d0*(sqrt(ctt(i,j))+sqrt(ctt(i,j-1))))
          fluxes(i,j,4)=coef_west*(Field(i,j)-Field(i,j-1))&
               /(z(i,j)-z(i,j-1))
       end do
    end do
  end function in_flux_surface_perp_diffusion


end module MDiffusion_perp
