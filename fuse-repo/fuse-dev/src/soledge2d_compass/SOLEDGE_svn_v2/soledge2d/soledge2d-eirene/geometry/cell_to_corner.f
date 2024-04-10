      subroutine eirene_cell_to_corner (f, fcorner)      
      use eirmod_precision
      use eirmod_parmmod
      use eirmod_ctrig
      use eirmod_cgrid
      use eirmod_cgeom
      use eirmod_ctetra
      use eirmod_ccona
      USE eirmod_CPOLYG
      USE eirmod_CLOGAU
      
      implicit none

      real(dp), intent(in) :: f(:)
      real(dp), intent(out) :: fcorner(:)

      real(dp), allocatable :: volsum(:)
      real(dp) :: dist, ages, xc, yc, dist1, dist2, dist3, dist4
      integer :: i, j, ir, ipart, ip, ic, in, nrk, ic1, ic2, ic3, ic4,
     .           it
      TYPE(CELL_ELEM), POINTER :: CUR

      
      if ((levgeo == 1) .and. nlrad .and. nlpol) then 

         nrk = indpoint(nr1st,np2nd)
         allocate(volsum(nrk))
         fcorner = 0._dp
         volsum = 0._dp

         IT = 1 
         DO IR=1,NR1STM
           XC = 0.5_DP * (RSURF(IR) + RSURF(IR+1)) 
           DO IP=1,NP2NDM
             YC = 0.5_DP * (PSURF(IP) + PSURF(IP+1)) 
             IN = IR + ((IP-1)+(IT-1)*NP2T3)*NR1P2

             IC1 = INDPOINT(IR,IP)
             dist1 = 1._DP/SQRT((XC-RSURF(IR))**2+
     .                          (YC-PSURF(IP))**2) 
             FCORNER(IC1) = FCORNER(IC1) + F(IN)*DIST1
             VOLSUM(IC1) = VOLSUM(IC1) + DIST1

             IC2 = INDPOINT(IR+1,IP)
             dist1 = 1._DP/SQRT((XC-RSURF(IR+1))**2+
     .                          (YC-PSURF(IP))**2) 
             FCORNER(IC2) = FCORNER(IC2) + F(IN)*DIST2
             VOLSUM(IC2) = VOLSUM(IC2) + DIST2

             IC3 = INDPOINT(IR+1,IP+1)
             dist3 = 1._DP/SQRT((XC-RSURF(IR+1))**2+
     .                          (YC-PSURF(IP+1))**2) 
             FCORNER(IC3) = FCORNER(IC3) + F(IN)*DIST3
             VOLSUM(IC3) = VOLSUM(IC3) + DIST3

             IC4 = INDPOINT(IR,IP+1)
             dist4 = 1._DP/SQRT((XC-RSURF(IR))**2+
     .                          (YC-PSURF(IP+1))**2) 
             FCORNER(IC4) = FCORNER(IC4) + F(IN)*DIST4
             VOLSUM(IC4) = VOLSUM(IC4) + DIST4
           END DO  ! ip
         END DO  ! ir

         fcorner(1:nrk) = fcorner(1:nrk)/(volsum(1:nrk)+eps60)
         deallocate (volsum)

      elseif ((levgeo == 2) .or. (levgeo == 3)) then

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
         call  eirene_exit_own(1)
      end if

      return
      end subroutine  eirene_cell_to_corner

         
