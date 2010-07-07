      function flaplace (r,s,t,xvec,yvec,zvec,tevec)

      use precision
      use ccona
      
      implicit none

      real(dp), intent(in) :: r, s, t
      real(dp), intent(in), dimension(27) :: xvec, yvec, zvec, tevec

      real(dp) :: u1(3), u2(3), u3(3), eu1(3),eu2(3),eu3(3),guik(3,3),
     .            gradphi(3), eo1(3), eo2(3), eo3(3), gradphi2(3),
     .            goik(3,3), dfguik(3,3,3)
      real(dp) :: dxdr, dxds, dxdt, dydr, dyds, dydt, dzdr, dzds, dzdt,
     .            e123, gikd2te, sumj, dtedxj, flaplace
      integer :: ii, j, k, l

      REAL(DP), save :: DFVEC(3,27), DF2VEC(3,3,27)
      real(dp), save :: rold=1.E20_DP, sold=1.E20_DP, told=1.E20_DP

      
      if (abs(r-rold) + abs(s-sold) + abs(t-told) > eps30) then

! berechne 1. Ableitungen der Shapefunktionen
         CALL DFTRIQUA(DFVEC, r, s, t)

! berechne 2. Ableitungen der Shapefunktionen
         CALL DF2TRIQUA(DF2VEC, r, s, t)

         rold = r
         sold = s
         told = t

      end if

      dxdr=sum(xvec*dfvec(1,:))
      dxds=sum(xvec*dfvec(2,:))
      dxdt=sum(xvec*dfvec(3,:))
      dydr=sum(yvec*dfvec(1,:))
      dyds=sum(yvec*dfvec(2,:))
      dydt=sum(yvec*dfvec(3,:))
      dzdr=sum(zvec*dfvec(1,:))
      dzds=sum(zvec*dfvec(2,:))
      dzdt=sum(zvec*dfvec(3,:))
      
      eu1 = (/ dxdr, dydr, dzdr /)
      eu2 = (/ dxds, dyds, dzds /)
      eu3 = (/ dxdt, dydt, dzdt /)
      
      e123 = eu1(1)*eu2(2)*eu3(3) + eu2(1)*eu3(2)*eu1(3) +
     .       eu3(1)*eu1(2)*eu2(3) - eu1(3)*eu2(2)*eu3(1) -
     .       eu2(3)*eu3(2)*eu1(1) - eu3(3)*eu1(2)*eu2(1)
      
      call kreuzprod(eu2,eu3,e123,eo1)
      call kreuzprod(eu3,eu1,e123,eo2)
      call kreuzprod(eu1,eu2,e123,eo3)
      
      guik(1,1)=dxdr*dxdr + dydr*dydr + dzdr*dzdr
      guik(1,2)=dxdr*dxds + dydr*dyds + dzdr*dzds
      guik(1,3)=dxdr*dxdt + dydr*dydt + dzdr*dzdt
      
      guik(2,1)=dxds*dxdr + dyds*dydr + dzds*dzdr
      guik(2,2)=dxds*dxds + dyds*dyds + dzds*dzds
      guik(2,3)=dxds*dxdt + dyds*dydt + dzds*dzdt
      
      guik(3,1)=dxdt*dxdr + dydt*dydr + dzdt*dzdr
      guik(3,2)=dxdt*dxds + dydt*dyds + dzdt*dzds
      guik(3,3)=dxdt*dxdt + dydt*dydt + dzdt*dzdt
      
      goik(1,1) = sum(eo1*eo1)
      goik(1,2) = sum(eo1*eo2)
      goik(1,3) = sum(eo1*eo3)
      
      goik(2,1) = sum(eo2*eo1)
      goik(2,2) = sum(eo2*eo2)
      goik(2,3) = sum(eo2*eo3)
      
      goik(3,1) = sum(eo3*eo1)
      goik(3,2) = sum(eo3*eo2)
      goik(3,3) = sum(eo3*eo3)
      
      u1 = 1._dp/sqrt(guik(1,1))*eu1
      u2 = 1._dp/sqrt(guik(2,2))*eu2
      u3 = 1._dp/sqrt(guik(3,3))*eu3
      
!     berechne Ableitungen der guik
      do ii=1,3
         do k=1,3
            do l=1,3
!     dguik(i,k) / dx_l
               dfguik(ii,k,l)=sum(xvec*df2vec(ii,l,:))*
     .              sum(xvec*dfvec(k,:))
     .              + sum(xvec*dfvec(ii,:))*
     .              sum(xvec*df2vec(k,l,:))
     .              + sum(yvec*df2vec(ii,l,:))*
     .              sum(yvec*dfvec(k,:)) 
     .              + sum(yvec*dfvec(ii,:))*
     .              sum(yvec*df2vec(k,l,:))
     .              + sum(zvec*df2vec(ii,l,:))*
     .              sum(zvec*dfvec(k,:)) 
     .              + sum(zvec*dfvec(ii,:))*
     .              sum(zvec*df2vec(k,l,:))
            end do
         end do
      end do
      
      gikd2te=0._dp
      do ii=1,3
         do k=1,3
            gikd2te = gikd2te +
     .           goik(ii,k) * sum(tevec*df2vec(ii,k,:))
         end do
      end do
      
      sumj = 0._dp
      do j = 1, 3
         dtedxj = sum(tevec*dfvec(j,:))
         do l = 1, 3
            do ii = 1,3
               do k = 1, 3
                  sumj = sumj + goik(ii,k) * goik(j,l) * dtedxj *
     .                 (dfguik(ii,l,k) + dfguik(k,l,ii) -
     .                 dfguik(ii,k,l))
               end do
            end do
         end do
      end do
      
      flaplace = gikd2te - 0.5_dp * sumj
      
      return
      end function flaplace
