      subroutine splbcb(xm,nx,ym,ny,fm,nxd,cx,ncx,cy,ncy,cc,ncc,
     ,                                                     wrk,nwr,ier)
c
c  version : 17.06.98 15:09
c
c======================================================================
c*** Calculation of coefficients of bi-cubic spline for data set on a
c*** regular grid, nx*ny points.
c*** The spline is a combination of 1D cubic splines, that is,
c*** the spline coefficients are first calculated for each column of the
c*** table and then a spline is created for each row of the matrix of
c*** the coefficients.
c***
c*** Input:
c***  nx, ny: actual dimensions of the input data grid
c***  nxd   : first dimension of the array "fm" in the calling routine
c***  xm, ym: data grid ordinates
c***  fm    : the data to be interpolated
c***  ncx   : dimension of the array "cx" in the calling routine,
c***          must be not less than (nx-1)*ny*3
c***  ncy   : dimension of the array "cy" in the calling routine,
c***          must be not less than (ny-1)*nx*3
c***  ncc   : dimension of the array "cc" in the calling routine,
c***          must be not less than (nx-1)*(ny-1)*9
c***  nwr   : dimension of the working array "wrk"
c***          must be not less than .......
c***
c*** Output:
c***  cx    : Spline coefficients for the columns of fm
c***  cy    : Spline coefficients for the rows of fm
c***  cc    : Spline coefficients for the rows of cx
c***  ier   : The return code (0 means OK)
c***
c*** Arrays xm, ym, fm, cx, and cc are used in the subroutine "evlbcs"
c*** to calculate the spline value for a particular point (x, y)
c***
c*** This routine uses subroutine "splcfn" to calculate the 1D splines
c======================================================================
      implicit none
      integer nx,ny,nxd,ncx,ncy,ncc,nwr,ier
      real xm(nx),ym(ny),fm(nx,ny),cx(nx-1,3,ny),cy(ny-1,3,nx),
     ,                                        cc(ny-1,3,nx-1,3),wrk(ny)

c*** Local variables

      integer i, j, k
c======================================================================
c*** Check the dimensions

      if(nx.le.1) then
        write(*,*) 'Error in splbcb: nx < 2'
        ier=1
        return
      end if

      if(ny.le.1) then
        write(*,*) 'Error in splbcb: ny < 2'
        ier=1
        return
      end if

      if(nx.gt.nxd) then
        write(*,*) 'Error in splbcb: nx > nxd'
        ier=1
        return
      end if

      if(ncx.lt.(nx-1)*ny*3) then
        write(*,*) 'Error in splbcb: ncx < (nx-1)*ny*3'
        ier=1
        return
      end if

      if(ncy.lt.(ny-1)*nx*3) then
        write(*,*) 'Error in splbcb: ncy < (ny-1)*nx*3'
        ier=1
        return
      end if

      if(ncc.lt.(nx-1)*(ny-1)*9) then
        write(*,*) 'Error in splbcb: ncx < (nx-1)*(ny-1)*9'
        ier=1
        return
      end if

      if(nwr.lt.ny) then
        write(*,*) 'Error in splbcb: nwr < ny'
        ier=1
        return
      end if

c----------------------------------------------------------------------
      ier=0

c*** Calculate spline coefficients for the data columns

      do j=1,ny
        call splcfn(xm,fm(1,j),nx,cx(1,1,j),nx-1,ier)
        if(ier.ne.0) then
          write(*,*) 'splbcb: error detected in splcfn during the ',
     ,                                            'first pass.  j = ',j
          return
        end if
      end do

c*** Calculate spline coefficients for the data rows

      do i=1,nx
        do j=1,ny
          wrk(j)=fm(i,j)
        end do
        call splcfn(ym,wrk,ny,cy(1,1,i),ny-1,ier)
        if(ier.ne.0) then
          write(*,*) 'splbcb: error detected in splcfn during the ',
     ,                                           'second pass.  i = ',i
          return
        end if
      end do

c*** Calculate spline coefficients for the coefficient rows

      do k=1,3
        do i=1,nx-1
          do j=1,ny
            wrk(j)=cx(i,k,j)
          end do
          call splcfn(ym,wrk,ny,cc(1,1,i,k),ny-1,ier)
          if(ier.ne.0) then
            write(*,*) 'splbcb: error detected in splcfn during the ',
     ,                                            'third pass.  i = ',i
            return
          end if
        end do
      end do
c======================================================================
      end
