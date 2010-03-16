      subroutine cell_to_corner (f, fcorner)      
      use precision
      use parmmod
      use ctrig
      use cgrid
      use cgeom
      use ctetra
      use ccona
      USE CPOLYG
      
      implicit none

      real(dp), intent(in) :: f(:)
      real(dp), intent(out) :: fcorner(:)

      real(dp), allocatable :: volsum(:)
      real(dp) :: dist, ages
      integer :: i, j, ir, ipart, ip, ic, in
      TYPE(CELL_ELEM), POINTER :: CUR


      if ((levgeo == 2) .or. (levgeo == 3)) then

         fcorner = 0.
         DO IR=1,NR1ST
           DO IPART=1,NPPLG
             DO IP=NPOINT(1,IPART),NPOINT(2,IPART)
               IC=INDPOINT(IR,IP)
               AGES = 0._DP
               FCORNER(IC) = 0.
               CUR => COORCELL(IC)%PCELL
               DO WHILE (ASSOCIATED(CUR))
                 IN=CUR%NOCELL
                 IF (NSTGRD(IN) == 0) THEN
                   DIST=1._DP/SQRT((XCOM(IN)-XPOL(IR,IP))**2+
     .                  (YCOM(IN)-YPOL(IR,IP))**2)
!pb                   AGES = AGES + DIST*XSTGRD(IN)
                   AGES = AGES + DIST
                   FCORNER(IC) = FCORNER(IC) + F(IN)*DIST
                 END IF
                 CUR => CUR%NEXT_CELL
               END DO
               FCORNER(IC) = FCORNER(IC) / (AGES+EPS60)
             END DO  ! ip
           END DO  ! ipart 
         END DO  ! ir
         
      elseif (levgeo == 4) then

         allocate(volsum(nrknot))
         volsum = eps60
         fcorner(1:nrknot) = 0.
         do i=1,ntrii
           do j=1,3
             DIST=1._DP/SQRT((XTRIAN(NECKE(J,I))-XCOM(I))**2+
     .                       (YTRIAN(NECKE(J,I))-YCOM(I))**2)  
             volsum(necke(j,i)) = volsum(necke(j,i)) + dist
             fcorner(necke(j,i)) = fcorner(necke(j,i)) + dist*f(i)
           end do
         end do
         fcorner(1:nrknot) = fcorner(1:nrknot)/volsum(1:nrknot)
         deallocate (volsum)

      elseif (levgeo.eq.5) then
         
         allocate(volsum(ncoord))
         volsum = eps60
         fcorner(1:ncoord) = 0.
         do i=1,ntet
           do j=1,4
             dist=1._dp/sqrt((xtetra(nteck(j,i))-xtcen(i))**2 +
     .                       (ytetra(nteck(j,i))-ytcen(i))**2 +
     .                       (ztetra(nteck(j,i))-ztcen(i))**2) 
             volsum(nteck(j,i))=volsum(nteck(j,i))+dist
             fcorner(nteck(j,i)) = fcorner(nteck(j,i)) + dist*f(i)
           enddo
         enddo
         fcorner(1:ncoor) = fcorner(1:ncoor)/volsum(1:ncoor)
         deallocate (volsum)

      else
         write (6,*) ' levgeo = ',levgeo,' to be written in',
     .               ' subroutine cell_to_corner '
         write (6,*) ' calculation abandonned '
         call exit_own(1)
      end if

      return
      end subroutine cell_to_corner

         
