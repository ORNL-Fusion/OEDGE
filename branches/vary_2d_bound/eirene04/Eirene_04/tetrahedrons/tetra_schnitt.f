


      subroutine TETRA_schnitt
     .  (ebene, tetra, spanz, spunkt)

! Funktion zur Berechnung saemtlicher Schnittpunkte einer Ebene und
! eines Tetraeders (max.4).
! - ebene ist ein Vektor mit den Parametern (a,b,c,d) der Ebene 
!   (ax+by+cz+d=0)
! - tetra ist eine 4x3-Matrix mit den (3-dim.) Eckpunkten des Tetraeders
! - spanz gibt die Anzahl der gefundenen Schnittpunkte an
! - in spunkt sind die gefundenen Schnittpunkte gespeichert

        USE PRECISION
        USE CCONA
        implicit none

        real(dp), intent(in), dimension(4)       :: ebene
        real(dp), intent(in), dimension(4,3)     :: tetra
        integer, intent(out)                     :: spanz
        real(dp), intent(out), dimension(4,3)    :: spunkt
        real(dp),dimension(3)                    :: punkt, SPP
        real(dp),dimension(5,5)                  :: dmat
        REAL(DP)                                 :: DST, DSTMIN
        integer                                  :: i, j, k, l, ip
        INTEGER,dimension(2)                     :: ipdst   
        logical                                  :: exi_spunkt

        spanz=0; i=1; j=1
        DSTMIN=HUGE(1.D0)
        dmat=DSTMIN
    !Schleife ueber alle moeglichen Kombinationen der Tetra-Eckpunkte
    !als Punkte einer Geraden
        do i=1,3
          j=(i+1)
          do
            if (j==5) then
              exit
            end if
          
          !Berechnung eines Schnittpunktes der Ebene und einer 
          !Kante des Tetraeders
!           call schnitt_ger_eb(spunkt(spanz+1,:), exi_spunkt,
!    .                          ebene, tetra(i,:), tetra(j,:))
            call schnitt_ger_eb(spp, exi_spunkt,
     .                          ebene, tetra(i,:), tetra(j,:))

          
          !Abfrage, ob Schnittpunkt gefunden
            if (exi_spunkt) then
               DO IP=1,SPANZ
                 dmat(ip,spanz+1)=SQRT(SUM((SPUNKT(IP,:)-SPP)**2))
               END DO
               dstmin=minval(dmat)
               IF (spanz < 4) then
                 if (minval(dmat(1:spanz,spanz+1)) > eps5) THEN
                   spanz = spanz+1
                   SPUNKT(SPANZ,:) = SPP
                 END IF
               ELSE 
                 ipdst=minloc(dmat)
                 IF (IPDST(2) .NE. SPANZ+1) THEN
                   SPUNKT(IPDST(2),:) = SPP
                   dmat(1:ipdst(2)-1,ipdst(2))=
     .                  dmat(1:ipdst(2)-1,spanz+1)
                   dmat(ipdst(2),ipdst(2)+1:spanz)=
     .                  dmat(ipdst(2)+1:spanz,spanz+1)
                 END IF
               END IF
            end if
          
            j = j+1
          end do
        end do
    

    !bei vier Schnittpunkten
        if (spanz == 4) then
      
          call sort_ueberpruef(spunkt)

        end if

      end subroutine TETRA_schnitt



!----------------------------------------------------------------------

      subroutine schnitt_ger_eb(schnittpunkt,schn_ctrl,a,x0,x1)

! Funktion zur Berechnung eines Schnittpunktes von einer Ebene
! und einer Kante eines Tetraeders:
! Wenn ein Schnittpunkt vorhanden ist, ist schn_ctrl = .true.
! und der Schnittpunkt ist in schnittpunkt gespeichert.
! - a ist der Vektor mit den Koeffizienten der Ebene
! - x0 und x1 sind die Eckpunkte der Kanten des Tetraeders

        USE PRECISION
        USE CCONA
        implicit none

        logical, intent(out)               :: schn_ctrl
        real(dp), intent(in), dimension(4)   :: a
        real(dp), intent(in), dimension(3)   :: x0, x1
        real(dp), intent(out), dimension(3)   :: schnittpunkt
        real(dp),dimension(3)   :: x
        real(dp),dimension(3)   :: r, rnorm
        real(dp)                :: schnittfakt, hilfe
        real(dp)                :: betrag0, betrag1, betrag
    
        schn_ctrl = .true.
        
        !Berechnung des Richtungsvektors der Geraden,
        !die durch die Kante des Tetraeders beschrieben ist
        r = x1 - x0      
    
        !Skalarprodukt von Normalenvektor der Ebene und 
        !Richtungsvektor der Geraden
        hilfe = dot_product(a(2:4),r)

        !Abstand von x0 zur Ebene
        betrag0 = a(1) + dot_product(a(2:4),x0)
    
        !Abstand von x1 zur Ebene
        betrag1 = a(1) + dot_product(a(2:4),x1)
    
        !Abstand zwischen x0 und x1 
        betrag = sqrt (sum (r * r))
        
        !normierter Richtungsvektor der Geraden 
        rnorm = r / betrag

        schnittfakt = huge(1._dp)    

        !Kante und Ebene sind parallel -> kein Schnittpunkt
!        if (abs(hilfe/(betraga*betragr)) < 1.E-4_DP) then
!        if (abs(hilfe) < 1.E-4_DP) then
!         if ((abs(hilfe) < eps5).or.
!     .       (betrag0*betrag1 > eps10) .or.
        if ((abs(dot_product(a(2:4),rnorm)) < eps5).or.
     .       (abs(betrag0-betrag1)/betrag < eps5)) then
           schn_ctrl = .false.
        else
           !Berechnung des Parameters der Geraden, der den 
           !Schnittpunkt festlegt
           schnittfakt = -( (a(1) + dot_product(a(2:4),x0)) / hilfe )
           if (abs(schnittfakt) < eps5) schnittfakt=0._dp
        end if
    
    
        !Abfrage, ob der Schnittpunkt auch auf der Kante des Tetraeders liegt
        if ((schnittfakt >= 0._dp) .and. 
     .      (schnittfakt <= 1._dp+eps5)) then
           !Berechnung des Schnittpunktes der Geraden und der Ebene
           x = x0 + schnittfakt * r
           schnittpunkt = x
        else
           schn_ctrl = .false.
        end if
           
      end subroutine schnitt_ger_eb


      subroutine sort_ueberpruef(punkt)
  
         USE PRECISION
         USE CCONA
         implicit none

         real(dp), intent(inout), dimension(4,3) :: punkt
         real(dp), dimension(3,2)                :: a
         real(dp), dimension(3)                  :: b, t, cen
         real(dp)                                :: x1, x2
         real(dp)                                :: d1, d2, d3, d4, d
         logical                                 :: lsw1, lsw2
         
         cen(1)=sum(punkt(1:4,1))*0.25_dp
         cen(2)=sum(punkt(1:4,2))*0.25_dp
         cen(3)=sum(punkt(1:4,3))*0.25_dp

         d1 = sqrt(sum((punkt(1,:)-cen)**2))
         d2 = sqrt(sum((punkt(2,:)-cen)**2))
         d3 = sqrt(sum((punkt(3,:)-cen)**2))
         d4 = sqrt(sum((punkt(4,:)-cen)**2))
         d = 1._dp/max(min(d1,d2,d3,d4),eps10)

         lsw1=.true.
         lsw2=.true.
         do while (lsw1 .or. lsw2)

!  test p1-p2 and p3-p4
            a(1,1) = punkt(2,1) - punkt(1,1)
            a(2,1) = punkt(2,2) - punkt(1,2)
            a(3,1) = punkt(2,3) - punkt(1,3)
            a(1,2) = punkt(3,1) - punkt(4,1)
            a(2,2) = punkt(3,2) - punkt(4,2)
            a(3,2) = punkt(3,3) - punkt(4,3)

            b(1) = punkt(3,1) - punkt(1,1)
            b(2) = punkt(3,2) - punkt(1,2)
            b(3) = punkt(3,3) - punkt(1,3)
            
            a = a*d
            b = b*d
            
            call solve_lgs

            if ((x1 >= 0._dp) .and. (x1 <= 1._dp) .and.
     .          (x2 >= 0._dp) .and. (x2 <= 1._dp)) then
!  intersection found, switch points 2 and 3
               t(1:3) = punkt(2,1:3)
               punkt(2,1:3) = punkt(3,1:3)
               punkt(3,1:3) = t(1:3)
               lsw1 = .true.
            else
               lsw1 = .false.
            end if

!  test p1-p4 and p2-p3
            a(1,1) = punkt(4,1) - punkt(1,1)
            a(2,1) = punkt(4,2) - punkt(1,2)
            a(3,1) = punkt(4,3) - punkt(1,3)
            a(1,2) = punkt(2,1) - punkt(3,1)
            a(2,2) = punkt(2,2) - punkt(3,2)
            a(3,2) = punkt(2,3) - punkt(3,3)

            b(1) = punkt(2,1) - punkt(1,1)
            b(2) = punkt(2,2) - punkt(1,2)
            b(3) = punkt(2,3) - punkt(1,3)
            
            a = a*d
            b = b*d
            
            call solve_lgs

            if ((x1 >= 0._dp) .and. (x1 <= 1._dp) .and.
     .          (x2 >= 0._dp) .and. (x2 <= 1._dp)) then
!  intersection found, switch points 3 and 4
               t(1:3) = punkt(3,1:3)
               punkt(3,1:3) = punkt(4,1:3)
               punkt(4,1:3) = t(1:3)
               lsw2 = .true.
            else
               lsw2 = .false.
            end if

         end do                 ! do while

         contains

         subroutine solve_lgs
           real(dp) :: c(2,2), d(2)
           real(dp) :: detc, detc1, detc2
           
           c(1,1) = a(1,1)**2 + a(2,1)**2 + a(3,1)**2
           c(1,2) = a(1,1)*a(1,2) + a(2,1)*a(2,2) + a(3,1)*a(3,2)
           c(2,1) = c(1,2)
           c(2,2) = a(1,2)**2 + a(2,2)**2 + a(3,2)**2

           d(1) = a(1,1)*b(1) + a(2,1)*b(2) + a(3,1)*b(3)
           d(2) = a(1,2)*b(1) + a(2,2)*b(2) + a(3,2)*b(3)
           
           detc  = c(1,1)*c(2,2) - c(1,2)*c(2,1)
           detc1 = d(1)*c(2,2) - d(2)*c(1,2)
           detc2 = c(1,1)*d(2) - c(2,1)*d(1)

           if (abs(detc) < 1.e-6_dp) then
             x1 = 10._dp
             x2 = 10._dp
           else
             x1 = detc1/detc
             x2 = detc2/detc
           end if
            
       end subroutine solve_lgs
         
       end subroutine sort_ueberpruef


