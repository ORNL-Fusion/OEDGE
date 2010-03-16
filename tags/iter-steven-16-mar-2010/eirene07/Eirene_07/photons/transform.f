      module transform
      use precision
      implicit none
      private

      public :: tr_init, tr_finalize, tr_make_ref, tr_show_ref,
     .   tr_make_tmatrix, tr_show_tmatrix, tr_atan
      public :: tr_cartsph, tr_rotate, tr_read, tr_write

! params
      integer, parameter :: dimens=3

! typedefs

! vardecls
      real, save :: pia
      integer, save :: nelems

      real, save, allocatable, dimension(:,:,:) :: refbasis, tmatrix

      contains

      subroutine tr_init(nnelems)
      implicit none
      integer, intent(in) :: nnelems
      write(*,*) 'tr_init:'
      if(dimens /= 3) then
         stop
      endif
      nelems=nnelems
      if(allocated(refbasis)) deallocate(refbasis)
      if(allocated(tmatrix)) deallocate(tmatrix)
      allocate(refbasis(nelems, dimens, dimens))
      allocate(tmatrix(nelems, dimens, dimens))
      refbasis=0.
      tmatrix=0.
  
      pia = 4.d0*atan(dble(1.d0))
      end subroutine tr_init

      subroutine tr_clear
      implicit none
       if(allocated(refbasis)) then
          refbasis=0.
          tmatrix=0.
       endif
      end subroutine tr_clear
   
      subroutine tr_finalize
       implicit none
       write(*,*) 'tr_finalize:'
       if(allocated(refbasis)) deallocate(refbasis)
       if(allocated(tmatrix)) deallocate(tmatrix)
      end subroutine tr_finalize

      subroutine tr_rotate(elem, v, mode)
       implicit none
       real(dp), intent(inout) :: v(:)
       integer, intent(in) :: mode,elem
       real, dimension(size(v)) :: vv
       real :: sum,t
       integer :: ir,ic
   
       do ir=1,dimens
          sum = 0.d0
          do ic=1,dimens
   
             if(mode >= 0) then
! transposed matrix = inverted matrix (cart.coords)
               t=tmatrix(elem,ic,ir)
             else
               t=tmatrix(elem,ir,ic)
             endif
   
             sum=sum+ t*v(ic)
   
          enddo
   
          vv(ir) = sum
       enddo
       v=vv
      end subroutine tr_rotate
   
      subroutine tr_cartsph(v, rho,theta,phi, mode)
       implicit none
       real(dp), intent(inout) :: v(:)
       real(dp), intent(inout) :: rho,theta,phi
       integer :: mode
       real :: r
       if(size(v) /= 3) then
          write(*,*) 'tr_cartsph: error, size(v) /= 3'
          stop
       endif
   
       if(mode >= 0) then
          rho   = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
          theta= tr_atan( dsqrt(v(1)*v(1) + v(2)*v(2))  , v(3))
          phi  = tr_atan(v(2), v(1))
!theta = tr_atan_dirk( v(3), sqrt(v(1)*v(1) + v(2)*v(2)) )
!phi   = tr_atan_dirk(v(1), v(2))
       else
          v(1) = rho*sin(theta)*cos(phi)
          v(2) = rho*sin(theta)*sin(phi)
          v(3) = rho*cos(theta)
       endif
      end subroutine tr_cartsph
   
      subroutine tr_show_tmatrix(elem)
       implicit none
       integer, intent(in) :: elem
       write(*,*) 'tmatrixelem=',elem
       write(*,'(3(2x,e14.6))') tmatrix(elem,1,:)
       write(*,'(3(2x,e14.6))') tmatrix(elem,2,:)
       write(*,'(3(2x,e14.6))') tmatrix(elem,3,:)
      end subroutine tr_show_tmatrix
   
      subroutine tr_make_tmatrix(elem)
       implicit none
       integer, intent(in) :: elem
       real, allocatable, dimension(:) :: v,w
       integer :: ir,ic
   
       allocate(v(dimens))
       allocate(w(dimens))
       tmatrix(elem,:,:) = 0.d0
   
! make transformation matrix from cartesian (system-)coordinates to new reference basis
       do ir = 1, dimens
          do ic = 1,dimens
   
             v = 0.d0
             v(ic) = 1.d0
   
             w = refbasis(elem,:,ir)
   
             tmatrix(elem, ir, ic) = tr_scalprod(w,v)
   
          enddo
       enddo
       
       deallocate(v)
       deallocate(w)
      end subroutine tr_make_tmatrix
   
      subroutine tr_show_ref(elem)
       implicit none
       integer, intent(in) :: elem
       write(*,*) 'refelem=',elem
       write(*,'(3(2x,e14.6))') refbasis(elem,:,1)
       write(*,'(3(2x,e14.6))') refbasis(elem,:,2)
       write(*,'(3(2x,e14.6))') refbasis(elem,:,3)
       write(*,*) 'a1 * a2 = ', 
     .  tr_scalprod(refbasis(elem,:,1), refbasis(elem,:,2))
       write(*,*) 'a2 * a3 = ', 
     .  tr_scalprod(refbasis(elem,:,2), refbasis(elem,:,3))
       write(*,*) 'a3 * a1 = ', 
     .  tr_scalprod(refbasis(elem,:,3), refbasis(elem,:,1))
       write(*,*) '(a1 curl a2) * a3 = ', 
     .  tr_spatprod( refbasis(elem,:,1), refbasis(elem,:,2), 
     .  refbasis(elem,:,3))
      end subroutine tr_show_ref
   
      subroutine tr_make_ref(elem, bb)
       implicit none
       integer, intent(in) :: elem
       real(dp), intent(in), dimension (:) :: bb
       real, dimension(size(bb)) :: b,a1,a2
       real :: bn,a1n,a2n
       real, parameter :: eps = 1.d-10
   
       if(elem < 0 .or. elem > nelems) return
       if(size(b) /= dimens) return
       b=bb
   
! assume b is new z-direction, z'
       bn = sqrt(b(1)**2 + b(2)**2 + b(3)**2)
       if(bn < 1.d-10) then
          write(*,*) 'tr_make_ref:'
          write(*,*) ' bn == 0, elem= ',elem
          !stop
          refbasis(elem,:,:)=0.d0
          return
       endif
       b=b/bn
      
! get two new right-handed orthogonal basis vectors
! first by hand, new x'
       if(b(1) > eps) then
          if(b(2) < eps .and. b(3) < eps) then
             a1(1) = 0.d0
             a1(2) = 1.d0
             a1(3) = 0.d0
          else
             a1(2) =  b(3)
             a1(3) = -b(2)
             a1(1) = -(a1(2)*b(2) + a1(3)*b(3)) /b(1)       
          endif
       elseif(b(2) > eps) then
          if(b(1) < eps .and. b(3) < eps) then
             a1(2) = 0.d0
             a1(3) = 1.d0
             a1(1) = 0.d0
          else
             a1(3) =  b(1)
             a1(1) = -b(3)
             a1(2) = -(a1(3)*b(3) + a1(1)*b(1)) /b(2)       
          endif
       elseif(b(3) > eps) then
          if(b(1) < eps .and. b(2) < eps) then
             a1(3) = 0.d0
             a1(1) = 1.d0
             a1(2) = 0.d0
          else
             a1(1) =  b(2)
             a1(2) = -b(1)
             a1(3) = -(a1(1)*b(1) + a1(2)*b(2)) /b(3)       
          endif
       endif
       a1n = sqrt(a1(1)**2 + a1(2)**2 + a1(3)**2)
       a1=a1/a1n
   
! second by curl, new y', right-handed
       a2(1) = -(a1(2)*b(3) - a1(3)*b(2))
       a2(2) = -(a1(3)*b(1) - a1(1)*b(3))
       a2(3) = -(a1(1)*b(2) - a1(2)*b(1))
       a2n = sqrt(a2(1)**2 + a2(2)**2 + a2(3)**2)
       a2=a2/a2n
       
! refbasis = 3 row vectors (a1 a2 b)
       refbasis(elem,:,1) = a1(:)
       refbasis(elem,:,2) = a2(:)
       refbasis(elem,:,3) =  b(:)   
      end subroutine tr_make_ref
   
      real function tr_scalprod(a,b) result(res)
       implicit none
       real, intent(in) :: a(:),b(:)
       integer :: i
       res=0.d0
       do i=1,size(a)
          res=res+ a(i) * b(i)
       enddo
      end function tr_scalprod
   
      real function tr_spatprod(a,b,c) result(res)
       implicit none
       real, intent(in) :: a(:), b(:), c(:)
       res=0.d0
       if(size(a) /= 3) return
   
       res=     a(1)*b(2)*c(3) + b(1)*c(2)*a(3) + c(1)*a(2)*b(3)
       res=res- a(1)*c(2)*b(3) - c(1)*b(2)*a(3) - b(1)*a(2)*c(3)
      end function tr_spatprod
   
      real function tr_atan_dirk(x,y) result(res)
       implicit none
       real, intent(in) :: x,y
       real :: z
       res=0.d0
       if (x.ge.0.d0.and.y.eq.0.d0) then
          res=0.d0
       else if (x.lt.0.d0.and.y.eq.0.d0) then
          res=pia
       else if (x.ge.0.d0.and.y.gt.0.d0) then
          z=x/y
          res=pia/2.d0-atan(z)
       else if (x.lt.0.d0.and.y.gt.0.d0) then
          z=x/y
          res=pia/2.d0-atan(z)
       else if (x.lt.0.d0.and.y.lt.0.d0) then
          z=x/y
          res=3.d0*pia/2.d0-atan(z)
       else if (x.ge.0.d0.and.y.lt.0.d0) then
          z=x/y
          res=3.d0*pia/2.d0-atan(z)
       endif
      end function tr_atan_dirk
   
      real(dp) function tr_atan(y,x) result(res)
       implicit none
       real(dp), intent(in) :: x,y
       real(dp) :: z
       res=0.d0
       if    (x == 0.d0 .and. y >0.d0) then
          res= pia/2.d0
       elseif(x == 0.d0 .and. y <0.d0) then
          res= 3.d0*pia/2.d0   
       
       elseif(y > 0.d0 .and. x>= 0.d0) then
          z=y/x
          res= atan(z)
       elseif(y < 0.d0 .and. x>= 0.d0) then
          z=y/x
          res= 2.d0*pia+atan(z)
   
       elseif(y > 0.d0 .and. x< 0.d0) then
          z=y/x
          res=  pia+atan(z)
       elseif(y < 0.d0 .and. x< 0.d0) then
          z=y/x
          res=  pia+atan(z)
   
       else
          !stop
          res=0.d0
       endif
      end function tr_atan
   
      subroutine tr_write(filename)
       implicit none
       character(len=*),intent(in) :: filename
       character(len=80) :: fname
       integer :: elem,rlen,rno,rnoi,rnoe,ir,n
       real :: x
       write(*,*) 'tr_write: ',filename
   
       fname = trim(filename)//trim('.transinfo')
       open(unit=77,file=fname)
       rewind(77)
       write(77,'(5(I6,1x))') nelems,dimens,kind(x)
       close(77)
   
       fname = trim(filename)//trim('.trans')
       rlen = 2*dimens**kind(x)
       open(unit=77,file=fname,
     .      access='DIRECT',form='UNFORMATTED',RECL=rlen)
       rno=1
       do elem=1,nelems
          rnoe = rno+ (elem-1) * dimens*dimens
          do ir=1,dimens
             rnoi = rnoe+ (ir-1)*dimens
             write(77, rec=rnoi)  
     .            (refbasis(elem,ir,n), n=1,dimens) , 
     .            (tmatrix(elem,ir,n), n=1,dimens)
          enddo
       enddo
       close(77)
      end subroutine tr_write
   
      subroutine tr_read(filename)
       implicit none
       character(len=*),intent(in) :: filename
       character(len=80) :: fname
       integer :: elem,rlen,rno,rnoi,rnoe,
     .            nelemsin,dimensin,kindin,ir,n
       real :: x
       logical :: ex
       write(*,*) 'tr_read: ', filename
   
       fname = trim(filename)//trim('.transinfo')
       inquire(file=fname, exist=ex)
       if(.not.ex) return
   
       open(unit=77,file=fname)
       read(77,'(5(I6,1x))') nelemsin,dimensin,kindin
       close(77)
   
       if(dimensin /= dimens ) then
          write(*,*) ' error, dimensin /= dimens'
          stop
       endif
       if(kindin /= kind(x) ) then
          write(*,*) ' error, kindin /= kind(x)'
          stop
       endif
   
       if(nelems /= nelemsin) then
          call tr_init(nelemsin)
       else
          call tr_clear
       endif
   
       fname = trim(filename)//trim('.trans')
       rlen = 2*dimens**kind(x)
       open(unit=77,file=fname,
     .              access='DIRECT',form='UNFORMATTED',RECL=rlen)
       rno=1
       do elem=1,nelems
          rnoe = rno+ (elem-1) * dimens*dimens
          do ir=1,dimens
             rnoi = rnoe+ (ir-1)*dimens
             read(77, rec=rnoi)  
     .            (refbasis(elem,ir,n), n=1,dimens) , 
     .            (tmatrix(elem,ir,n), n=1,dimens)
          enddo
       enddo
       close(77)
   
      end subroutine tr_read
      end module transform
