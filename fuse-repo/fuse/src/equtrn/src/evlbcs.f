      subroutine evlbcs(xs,nxs,ys,nys,ss,nxsd,
     ,                  xm,nx,ym,ny,fm,nxd,cx,ncx,cy,ncy,cc,ncc,
     ,                                            wrk,nwr,wrky,nwy,ier)
c
c  version : 17.06.98 12:48
c
c======================================================================
c*** This routine calculates the values of bi-cubic spline at the nodes
c*** of a regular grid in the (x,y) space.
c*** The spline coefficients must be calculated in advance with
c*** the "splbcb" routine.
c***
c*** Input:
c***  xs, ys    : arrays containing the co-ordinates of the points
c***              at which the spline is to be evaluated
c***  nxs, nys  : dimensions of these arrays (must be >0)
c***  nxsd      : the first dimension of the array "ss" in the calling
c***              routine (must be not less than nxs)
c***  xm,...,ncc: the same as in "splbcb"
c***  nwr       : dimension of the working array "wrk".
c***              Must be not less than 4*nx
c***  nwy       : dimension of the working array "wrky".
c***              Must be not less than ny
c***
c*** Output:
c***  ss        : the spline values
c***  ier       : return code (0 means OK, 1 if extrapolation occured,
c***              -1 means error in dimensions, no data calculated)
c***
c*** This routine uses subroutine "evalsp" to evaluate the 1D splines
c======================================================================
      implicit none
      integer nx,ny,nxd,ncx,ncy,ncc,ier,nxs,nys,nxsd,nwr,nwy
      real xs(nxs),ys(nys),ss(nxsd,nys),wrk(nx,4),wrky(ny)
      real xm(nx),ym(ny),fm(nx,ny),cx(nx-1,3,ny),cy(ny-1,3,nx),
     ,                             cc(ny-1,3,nx-1,3)

c*** Local variables

      integer i, j, k, l
c======================================================================
c*** Check the dimensions

      if(nxs.lt.1) then
        write(*,*) 'Error in evlbcs: nxs < 1'
        ier=-1
        return
      end if

      if(nys.lt.1) then
        write(*,*) 'Error in evlbcs: nys < 1'
        ier=-1
        return
      end if

      if(nxsd.lt.nxs) then
        write(*,*) 'Error in evlbcs: nxsd < nxs'
        ier=-1
        return
      end if

      if(nx.le.1) then
        write(*,*) 'Error in evlbcs: nx < 2'
        ier=-1
        return
      end if

      if(ny.le.1) then
        write(*,*) 'Error in evlbcs: ny < 2'
        ier=-1
        return
      end if

      if(nx.gt.nxd) then
        write(*,*) 'Error in evlbcs: nx > nxd'
        ier=-1
        return
      end if

      if(ncx.lt.(nx-1)*ny*3) then
        write(*,*) 'Error in evlbcs: ncx < (nx-1)*ny*3'
        ier=-1
        return
      end if

      if(ncy.lt.(ny-1)*nx*3) then
        write(*,*) 'Error in evlbcs: ncy < (ny-1)*nx*3'
        ier=-1
        return
      end if

      if(ncc.lt.(nx-1)*(ny-1)*9) then
        write(*,*) 'Error in evlbcs: ncx < (nx-1)*(ny-1)*9'
        ier=-1
        return
      end if

      if(nwr.lt.4*nx) then
        write(*,*) 'Error in evlbcs: nwr < 4*nx'
        ier=-1
        return
      end if

      if(nwy.lt.ny) then
        write(*,*) 'Error in evlbcs: nwy < ny'
        ier=-1
        return
      end if

c----------------------------------------------------------------------
      ier=0

      do j=1,nys

c*** Calculate the spline coefficients for the output column

        do i=1,nx
          do k=1,ny
            wrky(k)=fm(i,k)
          end do
          call evalsp(ym,wrky,ny,cy(1,1,i),ny-1,ys(j),wrk(i,1),1,ier)
        end do
        do i=1,nx-1
          do l=1,3
            do k=1,ny
              wrky(k)=cx(i,l,k)
            end do
            call evalsp(ym,wrky,ny,cc(1,1,i,l),ny-1,ys(j),
     ,                                                wrk(i,l+1),1,ier)
          end do
        end do

c*** Calculate the data values for the output column

        call evalsp(xm,wrk,nx,wrk(1,2),nx,xs,ss(1,j),nxs,ier)

      end do
c======================================================================
      end
