subroutine styx_radiation_weights(isurf,Ntri_tot,triRG,triZG,dOmega_tab,isWall)
  use eirmod_precision
  use eirmod_parmmod
  use styx2eirene
  implicit none

  integer, intent(in) :: isurf,Ntri_tot
  real*8, intent(in) :: triRG(Ntri_tot),triZG(Ntri_tot)
  real*8,intent(out) :: dOmega_tab(Ntri_tot)
  logical,intent(out) :: isWall
  real*8 :: dOmega
  integer*4 :: k,itri
  real*8 :: d(Ntri_tot),PWR(Ntri_tot),PWZ(Ntri_tot)
  real*8 :: Omega(Ntri_tot)
  real*8 :: RA,ZA,RB,ZB,RAs,ZAs
  real*8 :: ds,nZ,nR,norm,sinte,prop,dm,Numtris,nnR,nnZ,P,sinthe,Rm,Zm
  real*8 :: Rp,Zp,PWRp,PWZp,dpp,dsp,t,costheta

  real*8 :: d1,d2,d3
  real*8, dimension(3):: RBi,ZBi

  integer :: Numtri,kside,atend,atwall,next_side,im
  integer :: side1,side2,side3,Nsidesa
  integer :: icp

  integer :: calc_ang,iside
  
  real*8, parameter :: epsi = 1d-8

! debug
  integer :: icount
  icount=0

  Omega  = 0.d0
  isWall=.false.
  dOmega_tab=0.d0

  if (surface(isurf)%iprop /= 1) then

          itri=surface(isurf)%itri
          iside=surface(isurf)%iside
          k=isurf
          isWall=.true.


    
! center of wall segment
          Rm=0.5d0*(surface(k)%R1+surface(k)%R2)
          Zm=0.5d0*(surface(k)%Z1+surface(k)%Z2)
        
          ds=surface(k)%ds*100._dp
    
!  normal to wall segment
          if (abs(surface(k)%Z2-surface(k)%Z1) < epsi) then
            nZ=1
            nR=0
          elseif (abs(surface(k)%R2-surface(k)%R1)< epsi) then
            nZ=0
            nR=1
          else
            nZ=-(surface(k)%R2-surface(k)%R1)/(surface(k)%Z2-surface(k)%Z1)
            nR=1
            ! normalisation
            norm=sqrt(nR**2+nZ**2)
            nZ=nZ/norm
            nR=nR/norm
          endif

          !   check orientation (not correct when nx=0)
          sinte=(surface(k)%R2-surface(k)%R1)*nZ-(surface(k)%Z2-surface(k)%Z1)*nR
          if (sinte <0) then
            nR=-nR
            nZ=-nZ
          endif
!   distance to all triangles for this wall element
!   points from the wall to the triangle
    
          PWR=triRG-Rm
          PWZ=triZG-Zm
          d=sqrt(PWR**2+PWZ**2)

          do itri=1,Ntri_styx
          ! do itri=15736,15736

             !if (isurf == 404) write(*,*) 'itri = ',itri 
             
             atend=0
             atwall=0
        
             nnR=PWR(itri)/d(itri)
             nnZ=PWZ(itri)/d(itri)
        
             ! check orientation (not correct when nx=0)
             sinte=(surface(k)%R2-surface(k)%R1)*nnZ-(surface(k)%Z2-surface(k)%Z1)*nnR
             if (sinte <0) then
               nnR=-nnR
               nnZ=-nnZ
             endif

             ! check angle to normal (shadowed if negative)
        
             !sinte=(nnR*nR+nnZ*nZ)
             !if (sinte < 0) then
             !  write(*,*) 'wrong face intercepted, itri = ',itri
             !  write(*,*) 'surface(k)%R1 = ',surface(k)%R1
             !  write(*,*) 'surface(k)%R2 = ',surface(k)%R2
             !  write(*,*) 'surface(k)%Z1 = ',surface(k)%Z1
             !  write(*,*) 'surface(k)%Z2 = ',surface(k)%Z2
             !  write(*,*) 'nR = ',nR
             !  write(*,*) 'nZ = ',nZ
             !  write(*,*) 'nnR = ',nnR
             !  write(*,*) 'nnZ = ',nnZ
             !  write(*,*) 'triRG(itri) = ',triRG(itri)
             !  write(*,*) 'triZG(itri) = ',triZG(itri)
             !  write(*,*) ' itri = ',itri          
             !endif

             ! now follow line from there
        
             Numtri=surface(k)%itri
             kside=surface(k)%iside
        
        !tri_list=Numtri
        
             prop=0
             RA=Rm
             ZA=Zm
       
             do while (prop == 0 .and. Numtri /= itri)
 
                !if (isurf == 404) then 
                !  write(*,*) 'Numtri = ',Numtri
                !  icount=icount+1
                !  if (icount>1000) stop
                !endif

            ! intersection with other face
            
                if (kside == 1) then
                  sinthe=(xtriang(nvert(3,Numtri))-RA)*nnZ-(ytriang(nvert(3,Numtri))-ZA)*nnR
                  if (sinthe > 0) then
                    ! intersection on face 3
                    if ((xtriang(nvert(1,Numtri))-xtriang(nvert(3,Numtri))) /= 0) then
                        P=(ytriang(nvert(1,Numtri))-ytriang(nvert(3,Numtri)))/(xtriang(nvert(1,Numtri))-xtriang(nvert(3,Numtri)))
                        if (nnR /= 0._dp) then
                          RB=(ytriang(nvert(3,Numtri))-ZA+nnZ/nnR*RA-P*xtriang(nvert(3,Numtri)))/(nnZ/nnR-P)
                          ZB=ZA+nnZ/nnR*(ytriang(nvert(3,Numtri))-ZA+P*(RA-xtriang(nvert(3,Numtri))))/(nnZ/nnR-P)
                        else
                          ! vertical ray
                          RB=RA
                          ZB=P*(RA-xtriang(nvert(3,Numtri)))+ytriang(nvert(3,Numtri))
                        endif
                    else
                        RB=xtriang(nvert(1,Numtri))
                        ZB=ZA+nnZ/nnR*(RB-RA)
                    endif
                    next_side=3
                    prop=iprop(3,Numtri)
                  else
                    ! intersection on face 2
                    if ((xtriang(nvert(3,Numtri))-xtriang(nvert(2,Numtri))) /= 0) then
                        P=(ytriang(nvert(3,Numtri))-ytriang(nvert(2,Numtri)))/(xtriang(nvert(3,Numtri))-xtriang(nvert(2,Numtri)))
                        if (nnR /= 0._dp) then
                          RB=(ytriang(nvert(2,Numtri))-ZA+nnZ/nnR*RA-P*xtriang(nvert(2,Numtri)))/(nnZ/nnR-P)
                          ZB=ZA+nnZ/nnR*(ytriang(nvert(2,Numtri))-ZA+P*(RA-xtriang(nvert(2,Numtri))))/(nnZ/nnR-P)
                        else
                          RB=RA
                          ZB=P*(RA-xtriang(nvert(2,Numtri)))+ytriang(nvert(2,Numtri))
                        endif
                    else
                        RB=xtriang(nvert(2,Numtri))
                        ZB=ZA+nnZ/nnR*(RB-RB)
                    endif
                    next_side=2
                    prop=iprop(2,Numtri)
                  endif
                elseif (kside == 2) then
                  sinthe=(xtriang(nvert(1,Numtri))-RA)*nnZ-(ytriang(nvert(1,Numtri))-ZA)*nnR
                  if (sinthe > 0) then
                    ! intersection on face 1
                    if ((xtriang(nvert(2,Numtri))-xtriang(nvert(1,Numtri))) /= 0) then
                        P=(ytriang(nvert(2,Numtri))-ytriang(nvert(1,Numtri)))/(xtriang(nvert(2,Numtri))-xtriang(nvert(1,Numtri)))
                        if (nnR /= 0._dp) then
                          RB=(ytriang(nvert(1,Numtri))-ZA+nnZ/nnR*RA-P*xtriang(nvert(1,Numtri)))/(nnZ/nnR-P)
                          ZB=ZA+nnZ/nnR*(ytriang(nvert(1,Numtri))-ZA+P*(RA-xtriang(nvert(1,Numtri))))/(nnZ/nnR-P)
                        else
                          RB=RA
                          ZB=P*(RA-xtriang(nvert(1,Numtri)))+ytriang(nvert(1,Numtri))
                        endif
                    else
                        RB=xtriang(nvert(1,Numtri))
                        ZB=ZA+nnZ/nnR*(RB-RA)
                    endif
                    next_side=1
                    prop=iprop(1,Numtri)
                  else
                    ! intersection on face 3
                    if ((xtriang(nvert(1,Numtri))-xtriang(nvert(3,Numtri))) /= 0) then
                        P=(ytriang(nvert(1,Numtri))-ytriang(nvert(3,Numtri)))/(xtriang(nvert(1,Numtri))-xtriang(nvert(3,Numtri)))
                        if (nnR /= 0._dp) then
                          RB=(ytriang(nvert(3,Numtri))-ZA+nnZ/nnR*RA-P*xtriang(nvert(3,Numtri)))/(nnZ/nnR-P)
                          ZB=ZA+nnZ/nnR*(ytriang(nvert(3,Numtri))-ZA+P*(RA-xtriang(nvert(3,Numtri))))/(nnZ/nnR-P)
                        else
                          RB=RA
                          ZB=P*(RA-xtriang(nvert(3,Numtri)))+ytriang(nvert(3,Numtri))
                        endif
                    else
                        RB=xtriang(nvert(1,Numtri))
                        ZB=ZA+nnZ/nnR*(RB-RA)
                    endif
                    next_side=3
                    prop=iprop(3,Numtri)
                  endif
                elseif (kside == 3) then
                  sinthe=(xtriang(nvert(2,Numtri))-RA)*nnZ-(ytriang(nvert(2,Numtri))-ZA)*nnR
                  if (sinthe > 0) then
                    ! intersection on face 2
                    if ((xtriang(nvert(3,Numtri))-xtriang(nvert(2,Numtri))) /= 0) then
                        P=(ytriang(nvert(3,Numtri))-ytriang(nvert(2,Numtri)))/(xtriang(nvert(3,Numtri))-xtriang(nvert(2,Numtri)))
                        if (nnR /= 0._dp) then
                          RB=(ytriang(nvert(2,Numtri))-ZA+nnZ/nnR*RA-P*xtriang(nvert(2,Numtri)))/(nnZ/nnR-P)
                          ZB=ZA+nnZ/nnR*(ytriang(nvert(2,Numtri))-ZA+P*(RA-xtriang(nvert(2,Numtri))))/(nnZ/nnR-P)
                        else
                          RB=RA
                          ZB=P*(RA-xtriang(nvert(2,Numtri)))+ytriang(nvert(2,Numtri))
                        endif
                    else
                        RB=xtriang(nvert(2,Numtri))
                        ZB=ZA+nnZ/nnR*(RB-RA)
                    endif
                    next_side=2
                    prop=iprop(2,Numtri)
                  else
                    ! intersection on face 1
                    if ((xtriang(nvert(2,Numtri))-xtriang(nvert(1,Numtri))) /= 0) then
                        P=(ytriang(nvert(2,Numtri))-ytriang(nvert(1,Numtri)))/(xtriang(nvert(2,Numtri))-xtriang(nvert(1,Numtri)))
                        if (nnR /= 0._dp) then
                          RB=(ytriang(nvert(1,Numtri))-ZA+nnZ/nnR*RA-P*xtriang(nvert(1,Numtri)))/(nnZ/nnR-P)
                          ZB=ZA+nnZ/nnR*(ytriang(nvert(1,Numtri))-ZA+P*(RA-xtriang(nvert(1,Numtri))))/(nnZ/nnR-P)
                        else
                          RB=RA
                          ZB=P*(RA-xtriang(nvert(1,Numtri)))+ytriang(nvert(1,Numtri))
                        endif
                    else
                        RB=xtriang(nvert(1,Numtri))
                        ZB=ZA+nnZ/nnR*(RB-RA)
                    endif
                    next_side=1
                    prop=iprop(1,Numtri)
                  endif
                
                endif


                !write(*,*) 'RB = ',RB
                !write(*,*) 'ZB = ',ZB
                        
!% if core/edge interface hit, use target triangle and explore both sides
!% first towards cei (direction -nnR,--nnZ), then towards wall (direction nnR,nnZ)
!% one needs to know which face is intersected
            
                if (prop == 1) then
                
                  RAs=RA
                  ZAs=ZA
                
                  Numtri=itri
                  ! determining side
                  side1=0
                  side2=0
                  side3=0
                
                  ! intersection with line making side 1
                  if ((xtriang(nvert(2,Numtri))-xtriang(nvert(1,Numtri))) /= 0) then
                    P=(ytriang(nvert(2,Numtri))-ytriang(nvert(1,Numtri)))/(xtriang(nvert(2,Numtri))-xtriang(nvert(1,Numtri)))
                    if (nnR /= 0._dp) then
                      RBi(1)=(ytriang(nvert(1,Numtri))-ZA+nnZ/nnR*RA-P*xtriang(nvert(1,Numtri)))/(nnZ/nnR-P)
                      ZBi(1)=ZA+nnZ/nnR*(ytriang(nvert(1,Numtri))-ZA+P*(RA-xtriang(nvert(1,Numtri))))/(nnZ/nnR-P)
                    else
                      RBi(1)=RA
                      ZBi(1)=P*(RA-xtriang(nvert(1,Numtri)))+ytriang(nvert(1,Numtri))
                    endif
                  else
                    RBi(1)=xtriang(nvert(1,Numtri))
                    ZBi(1)=ZA+nnZ/nnR*(RB-RA)
                  endif
                  ! intersection with line making side 2
                  if ((xtriang(nvert(3,Numtri))-xtriang(nvert(2,Numtri))) /= 0) then
                    P=(ytriang(nvert(3,Numtri))-ytriang(nvert(2,Numtri)))/(xtriang(nvert(3,Numtri))-xtriang(nvert(2,Numtri)))
                    if (nnR /= 0._dp) then
                      RBi(2)=(ytriang(nvert(2,Numtri))-ZA+nnZ/nnR*RA-P*xtriang(nvert(2,Numtri)))/(nnZ/nnR-P)
                      ZBi(2)=ZA+nnZ/nnR*(ytriang(nvert(2,Numtri))-ZA+P*(RA-xtriang(nvert(2,Numtri))))/(nnZ/nnR-P)
                    else
                      RBi(2)=RA
                      ZBi(2)=P*(RA-xtriang(nvert(2,Numtri)))+ytriang(nvert(2,Numtri))
                    endif
                  else
                    RBi(2)=xtriang(nvert(2,Numtri))
                    ZBi(2)=ZA+nnZ/nnR*(RB-RA)
                  endif
                  ! intersection with line making side 3
                  if ((xtriang(nvert(1,Numtri))-xtriang(nvert(3,Numtri))) /= 0) then
                    P=(ytriang(nvert(1,Numtri))-ytriang(nvert(3,Numtri)))/(xtriang(nvert(1,Numtri))-xtriang(nvert(3,Numtri)))
                    if (nnR /= 0._dp) then
                      RBi(3)=(ytriang(nvert(3,Numtri))-ZA+nnZ/nnR*RA-P*xtriang(nvert(3,Numtri)))/(nnZ/nnR-P)
                      ZBi(3)=ZA+nnZ/nnR*(ytriang(nvert(3,Numtri))-ZA+P*(RA-xtriang(nvert(3,Numtri))))/(nnZ/nnR-P)
                    else
                      RBi(3)=RA
                      ZBi(3)=P*(RA-xtriang(nvert(3,Numtri)))+ytriang(nvert(3,Numtri))
                    endif
                  else
                    RBi(3)=xtriang(nvert(1,Numtri))
                    ZBi(3)=ZA+nnZ/nnR*(RB-RA)
                  endif

                  if (RBi(1) >= min(xtriang(nvert(2,Numtri)),xtriang(nvert(1,Numtri))) .and. RBi(1) <= max(xtriang(nvert(2,Numtri)),xtriang(nvert(1,Numtri))) &
                        .and. ZBi(1) >= min(ytriang(nvert(2,Numtri)),ytriang(nvert(1,Numtri))) .and. ZBi(1) <= max(ytriang(nvert(2,Numtri)),ytriang(nvert(1,Numtri)))) then
                    side1=1
                    d1=sqrt((Rm-RBi(1))**2+(Zm-ZBi(1))**2)
                  endif
                
                  if (RBi(2) >= min(xtriang(nvert(2,Numtri)),xtriang(nvert(3,Numtri))) .and. RBi(2) <= max(xtriang(nvert(2,Numtri)),xtriang(nvert(3,Numtri))) &
                        .and. ZBi(2) >= min(ytriang(nvert(2,Numtri)),ytriang(nvert(3,Numtri))) .and. ZBi(2) <= max(ytriang(nvert(2,Numtri)),ytriang(nvert(3,Numtri)))) then
                    side2=1
                    d2=sqrt((Rm-RBi(2))**2+(Zm-ZBi(2))**2)
                  endif
                
                  if (RBi(3) >= min(xtriang(nvert(3,Numtri)),xtriang(nvert(1,Numtri))) .and. RBi(3) <= max(xtriang(nvert(3,Numtri)),xtriang(nvert(1,Numtri))) &
                        .and. ZBi(3) >= min(ytriang(nvert(3,Numtri)),ytriang(nvert(1,Numtri))) .and. ZBi(3) <= max(ytriang(nvert(3,Numtri)),ytriang(nvert(1,Numtri)))) then
                    side3=1
                    d3=sqrt((Rm-RBi(3))**2+(Zm-ZBi(3))**2)
                  endif
                                
                ! get side number (the first one encountered ...)
                  if (side1==0) then
                    dm=minval([1d8, d2, d3],1)
                    im=minloc([1d8, d2, d3],1)
                    kside=im
                    RA=RBi(im)
                    ZA=ZBi(im)
                  elseif (side2==0) then
                    dm=minval([d1, 1d8, d3],1)
                    im=minloc([d1, 1d8, d3],1)
                    kside=im
                    RA=RBi(im)
                    ZA=ZBi(im)
                  elseif (side3==0) then
                    dm=minval([d1, d2, 1d8],1)
                    im=minloc([d1, d2, 1d8],1)
                    kside=im
                    RA=RBi(im)
                    ZA=ZBi(im)
                  endif
                                
                ! start again to follow lines, in both directions ...
                
                                                              
                ! now on the CEI side
                  RA=RAs
                  ZA=ZAs
                
                  Numtri=itri
                
                ! look for neighbor in side (to start inside a triangle)
                ! the possibility that it is on the CEI has to be considered
                ! here
                
                  Nsidesa=kside
                
                  if (Nsidesa == 1) then
                    kside=nside(1,Numtri)
                    prop=iprop(1,Numtri)
                    Numtri=neigh(1,Numtri)
                  elseif (Nsidesa == 2) then
                    kside=nside(2,Numtri)
                    prop=iprop(2,Numtri)
                    Numtri=neigh(2,Numtri)
                  elseif (Nsidesa == 3) then
                    kside=nside(3,Numtri)
                    prop=iprop(3,Numtri)
                    Numtri=neigh(3,Numtri)
                  endif
                
                ! shearch in reverse direction
                  nnR=-nnR
                  nnZ=-nnZ
                  
                  do while (prop == 0)
                    
                    ! intersection with other face
                    
                    if (kside == 1) then
                        sinthe=(xtriang(nvert(3,Numtri))-RA)*nnZ-(ytriang(nvert(3,Numtri))-ZA)*nnR
                        if (sinthe > 0) then
                            ! intersection on face 3
                            if ((xtriang(nvert(1,Numtri))-xtriang(nvert(3,Numtri))) /= 0) then
                                P=(ytriang(nvert(1,Numtri))-ytriang(nvert(3,Numtri)))/(xtriang(nvert(1,Numtri))-xtriang(nvert(3,Numtri)))
                                if (nnR /= 0._dp) then
                                  RB=(ytriang(nvert(3,Numtri))-ZA+nnZ/nnR*RA-P*xtriang(nvert(3,Numtri)))/(nnZ/nnR-P)
                                  ZB=ZA+nnZ/nnR*(ytriang(nvert(3,Numtri))-ZA+P*(RA-xtriang(nvert(3,Numtri))))/(nnZ/nnR-P)
                                else
                                  RB=RA
                                  ZB=P*(RA-xtriang(nvert(3,Numtri)))+ytriang(nvert(3,Numtri))
                                endif
                            else
                                RB=xtriang(nvert(1,Numtri))
                                ZB=ZA+nnZ/nnR*(RB-RA)
                            endif
                            next_side=3
                            prop=iprop(3,Numtri)
                        else
                            ! intersection on face 2
                            if ((xtriang(nvert(3,Numtri))-xtriang(nvert(2,Numtri))) /= 0) then
                                P=(ytriang(nvert(3,Numtri))-ytriang(nvert(2,Numtri)))/(xtriang(nvert(3,Numtri))-xtriang(nvert(2,Numtri)))
                                if (nnR /= 0._dp) then
                                  RB=(ytriang(nvert(2,Numtri))-ZA+nnZ/nnR*RA-P*xtriang(nvert(2,Numtri)))/(nnZ/nnR-P)
                                  ZB=ZA+nnZ/nnR*(ytriang(nvert(2,Numtri))-ZA+P*(RA-xtriang(nvert(2,Numtri))))/(nnZ/nnR-P)
                                else
                                  RB=RA
                                  ZB=P*(RA-xtriang(nvert(3,Numtri)))+ytriang(nvert(2,Numtri))
                                endif
                            else
                                RB=xtriang(nvert(2,Numtri))
                                ZB=ZA+nnZ/nnR*(RB-RA)
                            endif
                            next_side=2
                            prop=iprop(2,Numtri)
                        endif
                    elseif (kside == 2) then
                        sinthe=(xtriang(nvert(1,Numtri))-RA)*nnZ-(ytriang(nvert(1,Numtri))-ZA)*nnR
                        if (sinthe > 0) then
                            ! intersection on face 1
                            if ((xtriang(nvert(2,Numtri))-xtriang(nvert(1,Numtri))) /= 0) then
                                P=(ytriang(nvert(2,Numtri))-ytriang(nvert(1,Numtri)))/(xtriang(nvert(2,Numtri))-xtriang(nvert(1,Numtri)))
                                if (nnR /= 0._dp) then
                                  RB=(ytriang(nvert(1,Numtri))-ZA+nnZ/nnR*RA-P*xtriang(nvert(1,Numtri)))/(nnZ/nnR-P)
                                  ZB=ZA+nnZ/nnR*(ytriang(nvert(1,Numtri))-ZA+P*(RA-xtriang(nvert(1,Numtri))))/(nnZ/nnR-P)
                                else
                                  RB=RA
                                  ZB=P*(RA-xtriang(nvert(1,Numtri)))+ytriang(nvert(1,Numtri))
                                endif
                            else
                                RB=xtriang(nvert(1,Numtri))
                                ZB=ZA+nnZ/nnR*(RB-RA)
                            endif
                            next_side=1
                            prop=iprop(1,Numtri)
                        else
                            ! intersection on face 3
                            if ((xtriang(nvert(1,Numtri))-xtriang(nvert(3,Numtri))) /= 0) then
                                P=(ytriang(nvert(1,Numtri))-ytriang(nvert(3,Numtri)))/(xtriang(nvert(1,Numtri))-xtriang(nvert(3,Numtri)))
                                if (nnR /= 0._dp) then
                                  RB=(ytriang(nvert(3,Numtri))-ZA+nnZ/nnR*RA-P*xtriang(nvert(3,Numtri)))/(nnZ/nnR-P)
                                  ZB=ZA+nnZ/nnR*(ytriang(nvert(3,Numtri))-ZA+P*(RA-xtriang(nvert(3,Numtri))))/(nnZ/nnR-P)
                                else
                                  RB=RA
                                  ZB=P*(RA-xtriang(nvert(3,Numtri)))+ytriang(nvert(3,Numtri))
                                endif
                            else
                                RB=xtriang(nvert(1,Numtri))
                                ZB=ZA+nnZ/nnR*(RB-RA)
                            endif
                            next_side=3
                            prop=iprop(3,Numtri)
                        endif
                    elseif (kside == 3) then
                        sinthe=(xtriang(nvert(2,Numtri))-RA)*nnZ-(ytriang(nvert(2,Numtri))-ZA)*nnR
                        if (sinthe > 0) then
                            ! intersection on face 2
                            if ((xtriang(nvert(3,Numtri))-xtriang(nvert(2,Numtri))) /= 0) then
                                P=(ytriang(nvert(3,Numtri))-ytriang(nvert(2,Numtri)))/(xtriang(nvert(3,Numtri))-xtriang(nvert(2,Numtri)))
                                if (nnR /= 0._dp) then
                                  RB=(ytriang(nvert(2,Numtri))-ZA+nnZ/nnR*RA-P*xtriang(nvert(2,Numtri)))/(nnZ/nnR-P)
                                  ZB=ZA+nnZ/nnR*(ytriang(nvert(2,Numtri))-ZA+P*(RA-xtriang(nvert(2,Numtri))))/(nnZ/nnR-P)
                                else
                                  RB=RA
                                  ZB=P*(RA-xtriang(nvert(2,Numtri)))+ytriang(nvert(2,Numtri))
                                endif
                            else
                                RB=xtriang(nvert(2,Numtri))
                                ZB=ZA+nnZ/nnR*(RB-RA)
                            endif
                            next_side=2
                            prop=iprop(2,Numtri)
                        else
                            ! intersection on face 1
                            if ((xtriang(nvert(2,Numtri))-xtriang(nvert(1,Numtri))) /= 0) then
                                P=(ytriang(nvert(2,Numtri))-ytriang(nvert(1,Numtri)))/(xtriang(nvert(2,Numtri))-xtriang(nvert(1,Numtri)))
                                if (nnR /= 0._dp) then
                                  RB=(ytriang(nvert(1,Numtri))-ZA+nnZ/nnR*RA-P*xtriang(nvert(1,Numtri)))/(nnZ/nnR-P)
                                  ZB=ZA+nnZ/nnR*(ytriang(nvert(1,Numtri))-ZA+P*(RA-xtriang(nvert(1,Numtri))))/(nnZ/nnR-P)
                                else
                                  RB=RA
                                  ZB=P*(RA-xtriang(nvert(1,Numtri)))+ytriang(nvert(1,Numtri))
                                endif
                            else
                                RB=xtriang(nvert(1,Numtri))
                                ZB=ZA+nnZ/nnR*(RB-RA)
                            endif
                            next_side=1
                            prop=iprop(1,Numtri)
                        endif
                        
                    endif
                                     
                    
                    if (prop ==0) then
                        if (next_side == 1) then
                            kside=nside(1,Numtri)
                            Numtri=neigh(1,Numtri)
                        elseif (next_side == 2) then
                            kside=nside(2,Numtri)
                            Numtri=neigh(2,Numtri)
                        elseif (next_side == 3) then
                            kside=nside(3,Numtri)
                            Numtri=neigh(3,Numtri)
                        endif
                    endif
                    
                    RA=RB
                    ZA=ZB
                    
                enddo
                
               ! at this point the search must be finished for the direction considered %%%%
               ! provided the line crosses the core
                if (prop==1)then
                    atend=1
                elseif (prop > 1) then
                    atwall=1
                endif    
                
            endif ! prop == 1, search on the side of the core
            
            if (prop > 1) then
                    atwall=1
            endif    
            
            ! now find which the neighboring triangle and the corresponding side
            
            if (Numtri /= itri .and. atend /= 1 .and. atwall == 0) then
                
                Numtris=Numtri
                                
                if (next_side == 1) then
                    kside=nside(1,Numtri)
                    Numtri=neigh(1,Numtri)
                elseif (next_side == 2) then
                    kside=nside(2,Numtri)
                    Numtri=neigh(2,Numtri)
                elseif (next_side == 3) then
                    kside=nside(3,Numtri)
                    Numtri=neigh(3,Numtri)
                endif
                
            endif
            ! update coordinates
            
            RA=RB
            ZA=ZB    
            
        enddo ! do while
        
        
        if (Numtri==itri .or. atend == 1) then

! if the triangle is too close from the wall, cut the wall element into pieces and sum the contributions
! otherwise, too large difference between arc length R*Deltatheta and actual length ds

            calc_ang =1

            if (calc_ang == 0) then

            if (d(itri) < 10.d0*ds) then
                dsp=ds/100.d0
                dOmega=0.d0
            	do icp=1,100
                        ! parametric representation of the wall element
                        ! Zp,Rp coordinates of the center of elementary wall elements
                        t=(float(icp-1)+0.5d0)*dsp/ds
                	Rp=(1.d0-t)*surface(k)%R1+t*surface(k)%R2
                        Zp=(1.d0-t)*surface(k)%Z1+t*surface(k)%Z2
                        PWRp=triRG(itri)-Rp
                        PWZp=triZG(itri)-Zp
                        dpp=sqrt(PWRp**2+PWZp**2)
                        dOmega = dOmega + abs(PWRp*nR+PWZp*nZ)*dsp/(dpp**2)
                enddo

            else

              dOmega=abs(PWR(itri)*nR+PWZ(itri)*nZ)*ds/(d(itri)**2)
            
            endif

            else
                 
                costheta=((surface(k)%R1-triRG(itri))*(surface(k)%R2-triRG(itri))+               & 
                             (surface(k)%Z1-triZG(itri))*(surface(k)%Z2-triZG(itri)))/              &
                             (sqrt((surface(k)%R1-triRG(itri))**2+(surface(k)%Z1-triZG(itri))**2)*  & 
                             sqrt((surface(k)%R2-triRG(itri))**2+(surface(k)%Z2-triZG(itri))**2))
     
                ! costheta can be slightly above 1 because of numerical precision   
                costheta=min(costheta,1._dp)
                !if (costheta<0._dp) then
                !   write(*,*) ' cos(theta) negative, cos(theta) =',costheta
                !   write(*,*) 'surface(k)%R1 = ',surface(k)%R1
                !   write(*,*) 'surface(k)%R2 = ',surface(k)%R2
                !   write(*,*) 'surface(k)%Z1 = ',surface(k)%Z1
                !   write(*,*) 'surface(k)%Z2 = ',surface(k)%Z2
                !   write(*,*) 'nR = ',nR
                !   write(*,*) 'nZ = ',nZ
                !   write(*,*) 'nnR = ',nnR
                !   write(*,*) 'nnZ = ',nnZ
                !   write(*,*) 'triRG(itri) = ',triRG(itri)
                !   write(*,*) 'triZG(itri) = ',triZG(itri)
                !   write(*,*) ' itri = ',itri
                !endif 
                dOmega=Acos(costheta)
                
                if (dOmega<0_dp) then
                  write(*,*) 'dOmega < 0, itri = ',itri
                endif

            endif
            
            dOmega_tab(itri)=dOmega
            
            if (isnan(dOmega)) then
                 write(*,*) 'itri = ',itri
                 write(*,*) 'k    = ',k
                 write(*,*) 'surface(k)%R1-triRG(itri) = ',surface(k)%R1-triRG(itri)
                 write(*,*) 'surface(k)%R2-triRG(itri) = ',surface(k)%R2-triRG(itri)
                 write(*,*) 'surface(k)%Z1-triZG(itri) = ',surface(k)%Z1-triZG(itri)
                 write(*,*) 'surface(k)%Z2-triZG(itri) = ',surface(k)%Z2-triZG(itri)
                 write(*,*) 'cos(theta) = ',costheta
            endif     
           
        endif
 


     enddo ! itri
    
    endif ! iprop==1
 

end subroutine
