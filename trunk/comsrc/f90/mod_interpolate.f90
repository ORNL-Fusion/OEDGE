module mod_interpolate

  private

  ! set module debugging option
  logical :: debug_code = .false.

  public :: cell_interpolate,cell_area


contains



  subroutine distance( xstar,zstar,xalpha,zalpha,xbeta,zbeta, f )
    implicit none
    real*8 :: xstar,zstar,xalpha,zalpha,xbeta,zbeta,f
    real*8 :: m,x8,z8      

    real*8 :: nd,nx,nz,apx,apz,dn,dx,dz,dist

    !  to compute perp. distance from grid point (xstar,zstar) to line connecting
    !  points (xalpha,zalpha) and (xbeta,zbeta)

    !m=(zbeta-zalpha)/(xbeta-xalpha) +1.e-12   ! slope of line between alpha,beta points
    !m=(zbeta-zalpha)/(xbeta-xalpha +1.e-12)   ! slope of line between alpha,beta points
    !x8=(zstar +xstar/m  -zalpha +m*xalpha)/( m +1/m )
    !z8=zstar +(x8-xstar)/m
    !f=sqrt( (xstar-x8)**2 +(zstar-z8)**2 )


    !
    ! The following code is based on the vector solution to this problem rather than using slopes
    ! The issue with slopes is that they can be zero or infinity at two extremes which can cause issues
    ! in the interpolation code above. 
    ! The vector code produces exactly the same values as the above for non-edge cases but does produce
    ! better values in some of the edge cases where the 1e-12 comes into play. 
    !
    ! A = Alpha Point, P = Test Point)
    ! N = unit vector from Alpha to Beta
    ! Vector formulations dist = norm ( [A-P] - ([A-P]dot N) N ) 
    ! 

    
    nd = sqrt((xbeta-xalpha)**2+(zbeta-zalpha)**2)
    nx = (xbeta-xalpha)/nd
    nz = (zbeta-zalpha)/nd

    apx = xalpha-xstar
    apz = zalpha-zstar

    dn = apx*nx + apz*nz

    dx = apx - dn * nx
    dz = apz - dn * nz

    dist = sqrt (dx**2 + dz**2) 

    !if (debug_code) then 
    !   write(0,'(a,20(1x,g18.10))') 'DIST1:',xstar,zstar,xalpha,zalpha,xbeta,zbeta,m,x8,z8,f
    !   write(0,'(a,20(1x,g18.10))') 'DIST2:',nd,nx,nz,apx,apz,dn,dx,dz,f,dist
    !endif

    ! over-write f with the value in dist
    f = dist


    return
  end subroutine distance


  subroutine distance2(x0,y0,x1,y1,x2,y2,f)
    implicit none
    real*8 :: x0,y0,x1,y1,x2,y2,f
    real*8 :: dist

    ! formula for the distance from point (x0,y0) to the line defined by (x1,y1),(x2,y2)

    dist = sqrt((y2-y1)**2 + (x2-x1)**2)
    if (dist.ne.0.0) then 
       f = abs(((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1) / dist)
    else
       f = sqrt((x1-x0)**2 + (y1-y0)**2)
    endif

  end subroutine distance2




  subroutine cell_interpolate(rq,zq,quant,r,z,values)
    implicit none
    real*8 :: rq,zq,quant
    integer, parameter :: nvert = 4
    real*8 :: r(nvert), z(nvert)
    real*8 :: values(nvert)
    real*8 :: e1,e2,f1,f2,esum,fsum,v12,v34
    real*8 :: flow,fhigh,elow,ehigh
    real*8 :: flowt,fhight,elowt,ehight
    real*8, parameter :: eps = 1d-8
    integer :: in
    !
    ! This routine uses Jeff's interpolation algorithm
    !

    ! The r,z coordinates are the locations where value is specified
    ! The (1 -> 2) and (3 -> 4) sides should be across the field lines while sides (2 -> 3) and (4 -> 1) are parallel to the field lines
    ! The routine interpolates inside the polygon to obtain a value quant at the location rq,zq

    !    interpolation
    !    f = along field line - constant i moves along a r surface
    !    e = across the field line - constant j moves along a z surface


    ! Distance to side (1,2)
    call distance( rq,zq, r(1),z(1),r(2),z(2), flow )
    !call distance( r,z, r(i,j),zz(i,j),r(i,j+1),zz(i,j+1), flow )

    ! Distance to side (3,4)
    call distance( rq,zq, r(3),z(3),r(4),z(4), fhigh )
    !call distance( ru,zu, r(i+1,j),zz(i+1,j),r(i+1,j+1),zz(i+1,j+1), fhigh )

    ! Distance to side (2,3)
    call distance( rq,zq, r(2),z(2),r(3),z(3), elow )
    !call distance( ru,zu, r(i,j),zz(i,j),r(i+1,j),zz(i+1,j), elow )

    ! Distance to side (4,1)
    call distance( rq,zq, r(4),z(4),r(1),z(1), ehigh )
    !call distance( ru,zu, r(i,j+1),zz(i,j+1),r(i+1,j+1),zz(i+1,j+1), ehigh )

    esum = elow+ehigh
    fsum = flow+fhigh

    e1=ehigh/esum
    e2=elow/esum
    f1=fhigh/fsum
    f2=flow/fsum

    v12 = values(1) * e2 + values(2) * e1
    v34 = values(4) * e2 + values(3) * e1

    quant = v12 * f1 + v34 * f2

    if (debug_code) then 

       call distance2( rq,zq, r(1),z(1),r(2),z(2), flowt )
       if (abs(flow-flowt).gt.eps) then 
          write(0,*) 'Problem in distance: flow:',flow,flowt
       endif

       call distance2( rq,zq, r(3),z(3),r(4),z(4), fhight )
       if (abs(fhigh-fhight).gt.eps) then 
          write(0,*) 'Problem in distance: fhigh:',fhigh,fhight
       endif

       call distance2( rq,zq, r(2),z(2),r(3),z(3), elowt )

       if (abs(elow-elowt).gt.eps) then 
          write(0,*) 'Problem in distance: elow:',elow,elowt
       endif

       call distance2( rq,zq, r(4),z(4),r(1),z(1), ehight )
       if (abs(ehigh-ehight).gt.eps) then 
          write(0,*) 'Problem in distance: ehigh:',ehigh,ehight
       endif
       write(6,'(a,20(1x,g15.8))') 'CI:',elow,ehigh,esum,e1,e2,flow,fhigh,fsum,f1,f2,v12,v34,quant,(values(in),in=1,4)

    endif

    if (quant.lt.(minval(values)-eps).or.quant.gt.(maxval(values)+eps)) then 
       write(0,*) 'Interpolation ERROR: Quantity is outside range for cell:',quant,minval(values),maxval(values)
    endif

  end subroutine cell_interpolate


  real function cell_area(rv,zv,nv)
    implicit none
    integer :: nv
    real*8 :: rv(nv),zv(nv)
    integer :: iv,iv_next
    real*8 :: area
    
    area = 0.0
   
    do iv = 1,nv
       if (iv.eq.nv) then
          iv_next = 1
       else
          iv_next = iv+1
       endif
       area = area + (RV(iv_next)*ZV(iv)- RV(iv)*ZV(iv_next))
    end do

    cell_area = 0.5 * abs(area)

    !if (debug_code) then 
    !   write(0,'(a,i6,20(1x,g18.8))') 'Cell area:',nv,(rv(iv),zv(iv),iv=1,nv),cell_area
    !endif

    return

  end function cell_area



end module mod_interpolate
