subroutine styx_sample_sheath1D(itri,iside)
  use eirmod_precision
  use eirmod_parmmod
  use eirmod_comprt, only : VELX,VELY,VELZ,VEL,E0
  use eirmod_czt1, only : RSQDVP
  use eirmod_ctrig, only : PTRIX,PTRIY
  use eirmod_cstep
  use styx2eirene

  implicit none
  integer, intent(in) :: iside,itri
  real(dp) :: tau,ksi,alphaB
  real(dp) :: alpha,beta,V1,V2,V3,pi
  integer :: iu,il,im,itau,iksi,isurf
  integer :: opt

! isurf in recycling surfaces index

  isurf=recsurfinv(iside,itri)


! get sheath 1D dimensionless parameters
  tau    = sheath1D(isurf)%tau
  ksi    = sheath1D(isurf)%ksi

  pi=4._dp*atan(1._dp) 

  opt=2

  select case (opt)

    case(0)
! first try to test diagnostics: with these choices alpha=90°, beta=0° everywhere
   
      E0=2._dp*tistep(1,1,isurf)+3._dp*testep(1,isurf)
      VELX=ptrix(ipstep(1,isurf),irstep(1,isurf))
      VELY=ptriy(ipstep(1,isurf),irstep(1,isurf))
      VELZ=0._dp
      VEL=SQRT(E0)*RSQDVP(1)

    case(1)
 
! second test : same thing, but assume alpha=xx° and beta=yy°, then deduce VELX, VELY ..

      E0=2._dp*tistep(1,1,isurf)+3._dp*testep(1,isurf)
       
      alpha=45._dp*pi/180._dp
      beta =45._dp*pi/180._dp
 
      ! coordinates in the frame (upar,n,uperp) where uperp=uparxn, n normal to the triangle
      ! upar = B-(B.n)B

      V1 = cos(alpha)*cos(beta)
      V2 = -sin(alpha)
      V3 = cos(alpha)*sin(beta)

      ! change to coordinates frame (X,Y,Z)

      VELX = Bnu2xyz(1,1,isurf)*V1+Bnu2xyz(1,2,isurf)*V2+Bnu2xyz(1,3,isurf)*V3
      VELY = Bnu2xyz(2,1,isurf)*V1+Bnu2xyz(2,2,isurf)*V2+Bnu2xyz(2,3,isurf)*V3
      VELZ = Bnu2xyz(3,1,isurf)*V1+Bnu2xyz(3,2,isurf)*V2+Bnu2xyz(3,3,isurf)*V3

      VEL=SQRT(E0)*RSQDVP(1)

! now the real thing, but using average energy and angles
    case(2)

       ! first, find index of tau and ksi in the grids: binary search
      
      il=0
      iu=sheath1D_av%Ntau
     
      do while (iu-il > 1)
        im=0.5*(iu+il)
        if (tau >= sheath1D_av%tau(im)) then
          il=im
        else
          iu=im
        endif
       enddo

      itau=iu

      if (itau < 1) itau=1
      if (itau > sheath1D_av%Ntau) itau=sheath1D_av%Ntau

      il=0
      iu=sheath1D_av%Nksi
     
      do while (iu-il > 1)
        im=0.5*(iu+il)
        if (ksi >= sheath1D_av%ksi(im)) then
          il=im
        else
          iu=im
        endif
      enddo

      iksi=iu

      if (iksi < 1) iksi=1
      if (iksi > sheath1D_av%Nksi) iksi=sheath1D_av%Nksi

      ! then, get normalized energy and angles (rad)

      E0    = sheath1D(isurf)%E(itau,iksi)
      alpha = sheath1D(isurf)%alpha(itau,iksi)
      beta  = sheath1D(isurf)%beta(itau,iksi)
      
      ! denormalizing energy
      E0 = 0.5_dp*testep(1,isurf)/(sheath1D_av%ksi(iksi)*sheath1D_av%ksi(iksi))*E0

      ! coordinates in the frame (B-(Bn)n,n,uperp) where uperp=Bxn, n normal to the triangle

      V1 = cos(alpha)*cos(beta)
      V2 = -sin(alpha)
      V3 = cos(alpha)*sin(beta)

      ! change to coordinates frame (X,Y,Z)

      VELX = Bnu2xyz(1,1,isurf)*V1+Bnu2xyz(1,2,isurf)*V2+Bnu2xyz(1,3,isurf)*V3
      VELY = Bnu2xyz(2,1,isurf)*V1+Bnu2xyz(2,2,isurf)*V2+Bnu2xyz(2,3,isurf)*V3
      VELZ = Bnu2xyz(3,1,isurf)*V1+Bnu2xyz(3,2,isurf)*V2+Bnu2xyz(3,3,isurf)*V3

      VEL=SQRT(E0)*RSQDVP(1) 
 

  end select



end subroutine styx_sample_sheath1D
