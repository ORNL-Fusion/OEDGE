      function eirene_femint (fecken, icell, x, y, z) result(res)
      
      use eirmod_precision
      use eirmod_parmmod
      use eirmod_clogau
      use eirmod_cgrid
      use eirmod_cgeom
      use eirmod_cpolyg
      use eirmod_ctrig
      use eirmod_ctetra

      implicit none

      real(dp), intent(in) :: fecken(*)
      real(dp), intent(in) :: x, y, z
      integer, intent(in) :: icell

      real(dp) :: x1, x2, x3, x4, y1, y2, y3, y4, f1, f2, f3, f4, 
     .            z1, z2, z3, z4, r, s, t, u, res
      integer :: ir, ip, it, ia, ib

      res = 0._dp

      if (((levgeo == 1) .and. nlrad .and. nlpol) .or.
     .    (levgeo == 2) .or. (levgeo == 3)) then

        call eirene_ncelln(icell,ir,ip,it,ia,ib,nr1st,np2nd,nt3rd,nbmlt,
     .              nlrad,nlpol,nltor)

        if (levgeo == 1) then
          x1=rsurf(ir)
          x2=rsurf(ir+1)
          x3=rsurf(ir+1)
          x4=rsurf(ir)
          y1=psurf(ip)
          y2=psurf(ip)
          y3=psurf(ip+1)
          y4=psurf(ip+1)
        else if ((levgeo == 2) .or. (levgeo == 3)) then
          x1=xpol(ir+1,ip)
          x2=xpol(ir,ip)
          x3=xpol(ir,ip+1)
          x4=xpol(ir+1,ip+1)
          y1=ypol(ir+1,ip)
          y2=ypol(ir,ip)
          y3=ypol(ir,ip+1)
          y4=ypol(ir+1,ip+1)
        end if

        f1=fecken(INDPOINT(IR+1,IP))
        f2=fecken(INDPOINT(IR,IP))
        f3=fecken(INDPOINT(IR,IP+1))
        f4=fecken(INDPOINT(IR+1,IP+1))
        
        call eirene_xyz_to_rst(icell, x1, y1, 0._dp, x2, y2, 0._dp,
     .                  x3, y3, 0._dp, x4, y4, 0._dp, 
     .                  x, y, z, r, s, t, u)

        res = f1 * 0.25_dp * (1._dp - r) * (1._dp - s)
     .      + f2 * 0.25_dp * (1._dp + r) * (1._dp - s)
     .      + f3 * 0.25_dp * (1._dp + r) * (1._dp + s)
     .      + f4 * 0.25_dp * (1._dp - r) * (1._dp + s)
        
      else if (levgeo == 4) then

        x1=xtrian(necke(1,icell)) 
        x2=xtrian(necke(2,icell)) 
        x3=xtrian(necke(3,icell)) 
        y1=ytrian(necke(1,icell)) 
        y2=ytrian(necke(2,icell)) 
        y3=ytrian(necke(3,icell))
        f1=fecken(necke(1,icell)) 
        f2=fecken(necke(2,icell)) 
        f3=fecken(necke(3,icell))

        call eirene_xyz_to_rst(icell, x1, y1, 0._dp, x2, y2, 0._dp, x3,
     .                  y3, 0._dp, x4, y4, 0._dp, 
     .                  x, y, z, r, s, t, u)

        res = f1*r + f2*s + f3*t

      else if (levgeo == 5) then

        x1=xtetra(nteck(1,icell))
        x2=xtetra(nteck(2,icell))
        x3=xtetra(nteck(3,icell))
        x4=xtetra(nteck(4,icell))
        y1=ytetra(nteck(1,icell))
        y2=ytetra(nteck(2,icell))
        y3=ytetra(nteck(3,icell))
        y4=ytetra(nteck(4,icell))
        z1=ztetra(nteck(1,icell))
        z2=ztetra(nteck(2,icell))
        z3=ztetra(nteck(3,icell))
        z4=ztetra(nteck(4,icell))
        f1=fecken(nteck(1,icell))
        f2=fecken(nteck(2,icell))
        f3=fecken(nteck(3,icell))
        f4=fecken(nteck(4,icell))
        
        call eirene_xyz_to_rst(icell, x1, y1, z1, x2, y2, z2,  
     .                  x3, y3, z3, x4, y4, z4, x, y, z, r, s, t, u)

        res = f1*r + f2*s + f3*t + f4*u

      end if

      return
      end function eirene_femint

