      subroutine rdeqpb(lun,ngpr,ngpz,
     ,                             iret,nr,nz,rbtor,psilim,rgr,zgr,pfm)
c======================================================================
c*** Read the equilibrium data from TdeV
c***
c*** Input:
c***  lun     the logical unit number for the input
c***  ngpr    the maximum number of points in R direction
c***  ngpz    the maximum number of points in Z direction
c***
c*** Output:
c***  iret    return code (0 means OK)
c***  nr      the actual number of points in R direction
c***  nz      the actual number of points in Z direction
c***  rbtor   product of the toroidal magnetic field by R
c***  zgr     the Z values for the grid points
c***  pfm     the values of the (poloidal flux - separatrix flux)
c***
c*** If the input file contains no data on the toroidal field (old
c*** version), then the file "rbtor.dat" is looked for, and if there is
c*** no data there, then the user is asked to enter the value of R*Btor
c======================================================================
c
c  version : 27.06.98 20:47
c
      real*8 rgr(ngpr), zgr(ngpz), pfm(ngpr, ngpz)
c... toroidal field in tesla, radius in m
      real*8 rbtor,psilim
      logical errenc
      integer nnul
      parameter (nnul=72)
      character ul*(nnul)
c======================================================================
c*** Read the plasma equilibrium ...

      rbtor=0.
      psilim=0

 10   read(lun,'(a72)',end=980) ul
      if(ul.eq.' ') go to 10
      l=index(ul,'=')
      if(l.eq.0) go to 20
c----------------------------------------------------------------------
c*** Parse the header lines

      if(index(ul,'R*BTOR').gt.0) then
        read(ul(l+1:nnul),*,err=982) rbtor
      else if(index(ul,'PSILIM').gt.0) then
        read(ul(l+1:nnul),*,err=982) psilim
      else if(index(ul,'X1,2,N').gt.0) then
        read(ul(l+1:nnul),*,err=982) x1,x2,nr
      else if(index(ul,'Y1,2,N').gt.0) then
        read(ul(l+1:nnul),*,err=982) y1,y2,nz
      end if

      go to 10
c----------------------------------------------------------------------
c*** Check the dimensions and determine r,z co-ordinates of the grid

 20   errenc=.false.
      if(nr.gt.ngpr .or. nr.le.1) then
        errenc=.true.
        write(*,*) 'wrong value of nr: ',nr
      end if
      if(nz.gt.ngpz .or. nz.le.1) then
        errenc=.true.
        write(*,*) 'wrong value of nz: ',nz
      end if
      if(x1.ge.x2) then
        errenc=.true.
        write(*,*) 'wrong values of x1,x2: ',x1,x2
      end if
      if(y1.ge.y2) then
        errenc=.true.
        write(*,*) 'wrong values of y1,y2: ',y1,y2
      end if

      if(errenc) then
        write(*,*) 'rdeqpb: error(s) in the input file header detected'
        iret=8
        return
      end if

      d=(x2-x1)/(nr-1)
      x=x1
      do i=1,nr
        rgr(i)=x
        x=x+d
      end do

      d=(y2-y1)/(nz-1)
      y=y1
      do i=1,nz
        zgr(i)=y
        y=y+d
      end do
c----------------------------------------------------------------------
c*** Check the R*Btor. If not read, look for the file 'rbtor.dat'.
c*** If does not exist, ask user to enter

      if(rbtor.eq.0.) then
        inquire(file='rbtor.dat',exist=errenc)
        if(errenc) then
          open(lun+2,file='rbtor.dat')
          rewind(lun+2)
          read(lun+2,*,err=30) rbtor
          write(*,*) 'R*Btor value read from "rbtor.dat"'
          close(lun+2)
        end if
      end if

30    if(rbtor.eq.0.) then
        write(*,*) 'No R*Btor data found, neither in the input file ',
     ,                                             'nor in "rbtor.dat"'
        write(*,*) 'Please enter the value of R*Btor : \c'
        read(*,*) rbtor
      end if
      if(rbtor.eq.0.) then
        write(*,*) 'Zero value of R*Btor - program stops'
        iret=8
        return
      end if

c----------------------------------------------------------------------
c*** Read the psi values.

      read(lun,*,err=984,end=986) ((pfm(i,j),i=1,nr),j=1,nz)

c      do i=1,nr
c        do j=1,nz
c          pfm(i,j)=pfm(i,j)-psilim
c        end do
c      end do

      iret=0
      return
c======================================================================

 980  print *,'==== rdeqpb: no psi data found in the input file'
      iret=8
      return

 982  print *,'==== rdeqpb: error by parsing ',ul
      iret=8
      return

 984  print *,'==== rdeqpb: wrong data format in the input file'
      iret=8
      return
c
 986  print *,'==== rdeqpb: not enough data in the input file'
      iret=8
      return
c======================================================================
      end
