      subroutine splcfn(x,y,nx,c,ic,ier)
c
c  version : 17.03.97 14:48
c
c======================================================================
c*** Cubic spline interpolation
c***
c*** Input:
c***
c***  x    array containing the abscissae of the data points
c***       (x(i),y(i)) i=1,...,nx.
c***       Must be ordered so that  x(i) < x(i+1).
c***  y    array containing the function values
c***  nx   number of elements in x and y. Must be > 1.
c***  ic   row dimension of matrix c as
c***         specified in the calling program.
c***
c*** Output:
c***
c***  c    spline coefficients. c is an nx-1 by 3 matrix.
c***         the value of the spline approximation at t is
c***             s(t) = ((c(i,3)*d+c(i,2))*d+c(i,1))*d+y(i)
c***         where x(i) .le. t .lt. x(i+1) and d = t-x(i).
c***  ier  return code.
c***         ier = -1 : ic < nx-1
c***         ier = -2 : nx < 2
c***         ier = -3 : input abscissa are not ordered correctly
c======================================================================
      implicit none
      integer nx,ic,ier
      real x(nx),y(nx),c(ic,3)
      integer im1,i,jj,j,mm1,mp1,m,nm1,nm2
      real da,db,dt,g,ccn(3)
c======================================================================
c
      ier = -1
      nm1 = nx-1
      if(ic.lt.nm1) then
        write(*,*) 'splcfn: wrong dimensions of array c:',ic,' < ',nm1
        return
      end if
      ier = -2
      if(nx.lt.2) then
        write(*,*) 'splcfn: too few data points:',nx
        return
      end if
      ier = -3
      if (nx .eq. 2) then
c
c*** Linear interpolation
c
        if (x(1) .ge. x(2)) then
          write(*,*) 'splcfn: x(i) >= x(i+1) at i = 1'
          return
        end if
        ier = 0
        c(1,1) = (y(2)-y(1))/(x(2)-x(1))
        c(1,2) = 0.0
        c(1,3) = 0.0
        return
      else
c
c*** Real spline
c
        do m = 2,nm1
          mm1=m-1
          c(m,2) = x(m)-x(mm1)
          if (c(m,2).le.0.0) then
            write(*,*) 'splcfn: x(i) >= x(i+1) at i =',mm1
            return
          end if
          c(m,3) = (y(m)-y(mm1))/c(m,2)
        end do
        ccn(2) = x(nx)-x(nm1)
        if (ccn(2).le.0.0) then
          write(*,*) 'splcfn: x(i) >= x(i+1) at i =',nm1
          return
        end if
        ccn(3) = (y(nx)-y(nm1))/ccn(2)
        ier = 0
        nm2 = nx-2
        if (nx .le. 3) then
          c(1,3) = ccn(2)
          c(1,2) = c(2,2)+ccn(2)
          c(1,1) = ((c(2,2)+2.*c(1,2))*c(2,3)*ccn(2)+
     +                                         c(2,2)**2*ccn(3))/c(1,2)
        else
          c(1,3) = c(3,2)
          c(1,2) = c(2,2)+c(3,2)
          c(1,1) = ((c(2,2)+2.*c(1,2))*c(2,3)*c(3,2)+
     +                                         c(2,2)**2*c(3,3))/c(1,2)
          do m=2,nm2
            mp1=m+1
            mm1=m-1
            g = -c(mp1,2)/c(mm1,3)
            c(m,1) = g*c(mm1,1)+3.*c(m,2)*c(mp1,3)+3.*c(mp1,2)*c(m,3)
            c(m,3) = g*c(mm1,2)+2.*c(m,2)+2.*c(mp1,2)
          end do
        end if
        g = -ccn(2)/c(nm2,3)
        c(nm1,1) = g*c(nm2,1)+3.*c(nm1,2)*ccn(3)+3.*ccn(2)*c(nm1,3)
        c(nm1,3) = g*c(nm2,2)+2.*c(nm1,2)+2.*ccn(2)
        if (nx.le.3) then
          ccn(1)=2.*ccn(3)
          ccn(3)=1.
          g=-1./c(nm1,3)
        else
          g = c(nm1,2)+ccn(2)
          ccn(1) = ((ccn(2)+2.*g)*ccn(3)*c(nm1,2)+
     +                           ccn(2)**2*(y(nm1)-y(nx-2))/c(nm1,2))/g
          g = -g/c(nm1,3)
          ccn(3) = c(nm1,2)
        end if
        ccn(3) = g*c(nm1,2)+ccn(3)
        ccn(1) = (g*c(nm1,1)+ccn(1))/ccn(3)
        c(nm1,1) = (c(nm1,1)-c(nm1,2)*ccn(1))/c(nm1,3)
        do jj=1,nm2
          j = nm1-jj
          c(j,1) = (c(j,1)-c(j,2)*c(j+1,1))/c(j,3)
        end do
        do i=2,nm1
          im1 = i-1
          dt = c(i,2)
          da = (y(i)-y(im1))/dt
          db = c(im1,1)+c(i,1)-2.*da
          c(im1,2) = (da-c(im1,1)-db)/dt
          c(im1,3) = db/dt**2
        end do
        dt = ccn(2)
        da = (y(nx)-y(nm1))/dt
        db = c(nm1,1)+ccn(1)-2.*da
        c(nm1,2) = (da-c(nm1,1)-db)/dt
        c(nm1,3) = db/dt**2
      end if
c======================================================================
      end
