!***********************************************************************
!
!     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
!
!***********************************************************************
!----  4.4.96 Groten JFC=128 in JFC=256 geaendert
!                    (Produkt der dim darf max 256*256 sein)
!---- 04.03.98 Groten Schleife "Farbbereiche" geaendert + NEAREST
!              wegen nah benachbarter/gleicher Werte in Rasterviereck.
!              Fortran 90

      subroutine grhhfl (nx,ix,x,ny,iy,y,nw,w,jco,n1,f,intact)
!
      implicit none
! --- Arguments
      integer, intent(in)   :: nx,ix,ny,iy,nw,n1,intact
      real,dimension(n1,(ny-1)*iy+1),intent(in) :: f
      real,dimension(nw),intent(in)             :: w
      integer,dimension(nw+1),intent(in)        :: jco
      real,dimension((nx-1)*ix+1),intent(in)    :: x
      real,dimension((ny-1)*iy+1),intent(in)    :: y
! --- local Variables
      logical :: logi
      integer :: i,j,jj,ii,iq,k,ifa
      real :: pp5,pp6,pp7,pp8,xmi,xma,ymi,yma,xh,yh,zmi,zma,h,h1,delx
      real :: dely,z
! --- Common grdyn
! --- named Constants, same as in GRHHNL
      integer,parameter :: jfc = 256, ifc = 512
      integer,parameter :: maw=jfc, maxq=maw*maw
      integer,parameter :: id1=5*(ifc)**2/4-2*maxq, id=(ifc*6)*1.2
      integer,parameter :: id2=14*ifc+id*2-2*maw
      real,dimension(2, maxq)  :: q
      real,dimension(maw+2)    :: we
      integer,dimension(maw+1) :: ico
      real,dimension(id1)      :: dum1
      real,dimension(id2)      :: dum2 
      common / grdyn / q, we, ico, dum1, dum2
      save :: / grdyn /
!DEC$ PSECT /GRDYN/ NOSHR
! --- Common grpp
      real,dimension(15)       :: pp
      integer                  :: icolor
      real,dimension(2)        :: qq
      common / grpp / pp, icolor, qq
!DEC$ PSECT /GRPP/ NOSHR
      save ::  / grpp /

!     Executable

      if (nx*ny > maxq) then
         print *, "GRHHFL: Rasterung beim Aufruf zu dicht."
         stop 
      end if
      if (nw > maw) then
         print *, "GRHHFL: Zu viele Farbschichten."
         stop 
      end if
!
      pp5 = pp (5)
      pp6 = pp (6)
      pp7 = pp (7)
      pp8 = pp (8)
!
      call grhhtr (x(1), y(1), xh, yh)
      xmi = xh
      xma = xh
      ymi = yh
      yma = yh
!     Transformation, Extrema x,y
      do j = 1, ny
         jj = (j-1) * iy + 1
         do i = 1, nx
            ii = (i-1) * ix + 1
            call grhhtr (x(ii), y(jj), xh, yh)
            if (xh < xmi) then
               xmi = xh
            else if (xh > xma) then
               xma = xh
            end if
            if (yh < ymi) then
               ymi = yh
            else if (yh > yma) then
               yma = yh
            end if
            iq = nx * (j-1) + i
            q (1, iq) = xh
            q (2, iq) = yh
         enddo
      enddo
!     Monoton ?
      logi = .false.
      do i = 2, nw
         logi = logi .or. w (i) <= w (i-1)
      enddo
      if (logi) then
         print *, "GRHHFL: Die Werte sind nicht monoton steigend."
         stop 
      end if
!     Extrema z
      zmi = f (1, 1)
      zma = zmi
      do j = 1, ny
         do i = 1, nx
            z = f ((i-1)*ix+1, (j-1)*iy+1)
            if (z < zmi) then
               zmi = z
            else if (z > zma) then
               zma = z
            end if
         enddo
      enddo
! --- Farbbereiche
      we (1) = zmi
      ico (1) = jco (nw+1)
      k = 1
      do i = 1, nw
         if (w(i) > zmi) then
            ico (k) = jco (i)
            if (w(i) < zma) then
               k = k + 1
               we (k) = w (i)
               ico (k) = jco (i+1)
            else
               exit
            end if
         end if
      enddo
      we (k+1) = zma
! --- Bildausschnitt
      if (intact == 0) then
         call grsclv (xmi, ymi, xma, yma)
      else if (intact /= 1234567890) then
         dely = (pp(4)-pp(2)) * (xma-xmi) / (pp(3)-pp(1))
         if ((yma-ymi) < dely) then
            call grsclv(xmi,(yma+ymi-dely)*0.5,xma,(yma+ymi+dely)*0.5)
         else
            delx = (pp(3)-pp(1)) * (yma-ymi) / (pp(4)-pp(2))
            call grsclv((xma+xmi-delx)*0.5,ymi,(xma+xmi+delx)*0.5,yma)
         end if
      end if
! --- Einfaerben
      ifa = icolor
      do i = 1, k
         call grnwpn (ico(i))
         h = we(i)
         we(i) = nearest(we(i),-1.0)
         h1 = we(i+1)
         we(i+1) = nearest(we(i+1),1.0)
         call gr2isl (nx, ix, ny, iy, we(i), n1, f, q)
         we(i) = h
         we(i+1) = h1
      enddo
! --- alter Stand
      call grnwpn (ifa)
      if (intact .ne. 0) then
         call grsclv (pp5, pp6, pp7, pp8)
      end if
!
      end subroutine grhhfl
