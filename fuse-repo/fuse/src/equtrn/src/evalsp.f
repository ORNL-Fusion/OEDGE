      subroutine evalsp(x,y,nx,c,ic,u,s,m,ier)
c
c  version : 11.03.97 21:14
c
c======================================================================
c*** Evaluation of a cubic spline from coefficients calculated
c*** in advance. The spline is represented in the simplest way,
c***    s=((c(k,3)*d+c(k,2))*d+c(k,1))*d+y(k),   where d=u-x(k).
c***
c*** Input
c***  x   - vector containing the abscissae of the data points;
c***        must be ordered so that x(i) < x(i+1)
c***  y   - vector containing the function values
c***  nx  - number of elements in x and y, must be > 1
c***  c   - spline coefficients. c is an nx-1 by 3 matrix.
c***  ic  - row dimension of matrix c as specified in the dimension
c***        statement in the calling program, must be > nx-2
c***  u   - vector containing the x values at which the cubic spline
c***        is to be evaluated
c***  m   - number of elements in u and s
c***
c*** Output
c***
c***  s   - interpolated values
c***        the value of the spline approximation at u(i) is
c***          s(i) = ((c(j,3)*d+c(j,2))*d+c(j,1))*d+y(j)
c***        where  x(j) .le. u(i) .lt. x(j+1)  and  d = u(i)-x(j)
c***  ier - error parameter, =1 if  u(i) < x(1)  or  u(i) > x(nx)
c***
c***  The abscissae x(i) of the data points must be in accending order.
c***  No check of this condition is made in the routine!
c======================================================================
      implicit none
      integer nx,ic,m,ier
      real x(nx),y(nx),c(ic,3),u(m),s(m)
      integer i,jer,ker,nx1,k
      real d,dd
      data i/1/
c======================================================================
c
      if (m .le. 0) return
      jer = 0
      ker = 0
      nx1 = nx-1
      if (i .gt. nx1) i = 1
c*** evaluate spline at m points
      do 40 k=1,m
c*** find the interval
         d = u(k)-x(i)
         if (d) 5,25,15
    5    if (i .eq. 1) go to 30
         i = i-1
         d = u(k)-x(i)
         if (d) 5,25,20
   10    i = i+1
         d = dd
   15    if (i .ge. nx) go to 35
         dd = u(k)-x(i+1)
         if (dd .ge. 0.) go to 10
         if (d .eq. 0.) go to 25
c*** perform evaluation
   20    s(k) = ((c(i,3)*d+c(i,2))*d+c(i,1))*d+y(i)
         go to 40
   25    s(k) = y(i)
         go to 40
c*** u(i) <  x(1)  -  extrapolation
   30    jer = 1
         go to 20
c*** u(i) >  x(nx)  -  extrapolation
   35    if (dd .gt. 0.) ker = 1
         d = u(k)-x(nx1)
         i = nx1
         go to 20
   40 continue
      ier = max0(jer,ker)
      if (jer .gt. 0) print *,'*** evalsp: u < x(1)'
      if (ker .gt. 0) print *,'*** evalsp: u > x(nx)'
c======================================================================
      end
