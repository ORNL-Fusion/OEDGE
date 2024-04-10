!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutines which computes average incidence angle for each surface element  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine styx_incidence_diag(itri,iside)
  use eirmod_precision
  use eirmod_parmmod
  use eirmod_ctrig
  use eirmod_comprt
  use eirmod_cstep

  use styx2eirene
  implicit none
  integer, intent(in) :: itri,iside 
  integer :: isurf
  real(dp) :: VsX,VsY,VsZ,Vs,Vdn
  real(dp) :: cosa,cosb
  real(dp) :: VB,Vn,Vu
  integer :: opt

 
! same numerotation system as in recsurf

  isurf=recsurfinv(iside,itri)
  
  opt=1

  if (opt==0) then

!!!! first way of calculating angles !!!!!!!

    Vdn = VELX*PTRIX(iside,itri)+VELY*PTRIY(iside,itri)

    VsX = VELX -Vdn*PTRIX(iside,itri)
    VsY = VELY -Vdn*PTRIY(iside,itri)
    VsZ = VELZ

    Vs=sqrt(VsX*VsX+VsY*VsY+VsZ*VsZ)

    cosa = min((VsX*VELX+VsY*VELY+VsZ*VELz)/Vs,1._dp)
    cosa = max(cosa,-1._dp)

    cosb =min((VsX*sheath1D(isurf)%uparX+VsY*sheath1D(isurf)%uparY+VsZ*sheath1D(isurf)%uparZ),1._dp)
    cosb= max(cosb,-1._dp)

    sheath1D(isurf)%alpha_V = sheath1D(isurf)%alpha_V + Acos(cosa)
    sheath1D(isurf)%beta_V = sheath1D(isurf)%beta_V + Acos(cosb)
    sheath1D(isurf)%hit_V = sheath1D(isurf)%hit_V+1

!!!! second way: use change of coordinate matrix !!!!

  else

    sheath1D(isurf)%hit_V=sheath1D(isurf)%hit_V+1

    VB=xyz2Bnu(1,1,isurf)*VELX+xyz2Bnu(1,2,isurf)*VELY+xyz2Bnu(1,3,isurf)*VELZ
    Vn=xyz2Bnu(2,1,isurf)*VELX+xyz2Bnu(2,2,isurf)*VELY+xyz2Bnu(2,3,isurf)*VELZ
    Vu=xyz2Bnu(3,1,isurf)*VELX+xyz2Bnu(3,2,isurf)*VELY+xyz2Bnu(3,3,isurf)*VELZ

    ! take care of small numerical errors if normal angles     
    Vn=min(Vn,1._dp)
    Vn=max(Vn,-1._dp)

    sheath1D(isurf)%alpha_V = sheath1D(isurf)%alpha_V - Asin(Vn)

    if (abs(Vu) < 1e-10_dp) Vu=0._dp

    if (VB /= 0._dp) then
      sheath1D(isurf)%beta_V = sheath1D(isurf)%beta_V + Atan(Vu/VB)
    else
      sheath1D(isurf)%beta_V=0._dp
    endif

  endif


end subroutine
