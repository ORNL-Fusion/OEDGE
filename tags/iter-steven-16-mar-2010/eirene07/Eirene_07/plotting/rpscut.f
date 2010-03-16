      subroutine rpscut (aorig,values)

      use precision
      use parmmod
      use ctetra
      use ctrig
      use cgeom
      use ccona
      use cplot
      use module_avltree

      implicit none

      REAL(DP), INTENT(IN) :: AORIG(*)
      REAL(DP), INTENT(OUT) :: VALUES(*)

      integer, allocatable, save :: tetra_kanten(:,:), kanten(:,:), 
     .             iedge(:), icono(:), itetno(:)
      integer, save :: kanten_nummer(4,4) = reshape(
     .    (/ -1,  1,  2,  3,
     .        1, -1,  4,  5,
     .        2,  4, -1,  6,
     .        3,  5,  6, -1 /), (/ 4, 4 /) )
      integer, save :: angrenzende_seiten(2,6) = reshape (
     .     (/ 1, 2,
     .        1, 4,
     .        2, 4,
     .        1, 3,
     .        2, 3,
     .        3, 4 /), (/ 2, 6 /) )
      integer, save :: nkanten, ifirst=0
      real(dp) :: timi, timen, second_own
      real(dp), allocatable, save :: spar(:)
      integer :: i, itri
 
      
!pb   if (.not.allocated(tetra_kanten)) then
      if (ifirst == 0) then
        ifirst = 1
        allocate (tetra_kanten(6,ntet))
        allocate (kanten(2,3*ntet))
        tetra_kanten=0
        kanten=0
        nkanten=0
!        call suche_kanten
!      timen = second_own()
!      write (0,*) ' cpu time spend in suche_kanten ',timen-timi,' sec'
!        kanten=0
!        nkanten=0
        timi = second_own()
        call suche_kanten2
        timen = second_own()
        write (0,*) ' cpu time spend in suche_kanten2 ',
     .                timen-timi,' sec'

        call schneide_kanten 
        call berechne_koordinaten
        call bilde_dreiecke 

        deallocate (tetra_kanten)
        deallocate (kanten)
        deallocate (iedge)
        deallocate (icono)
        deallocate (spar)
      end if

      do i=1,ntrii
        values(i) = aorig(itetno(i))
      end do

      return


      contains


      subroutine suche_kanten

      implicit none
      type(tet_elem), pointer :: cur
      integer :: itet, j, nxt_side, akt_tet, akt_edge, isi, ip1
      integer :: ic, i, nump1, nump2, noedge, jc, no_side, nxt_tet

      nkanten = 0
! für jede Koordinate
      do ic=1,ncoord

! durchsuche die Kanten der angrenzenden Tetraeder
        cur => coortet(ic)%ptet
        do while (associated(cur))
          itet = cur%notet

! Koordinate ic ist Eckpunkt des Tetraeders itet, 
! bestimme die Nummer der Ecke
          do i=1,4
            if (nteck(i,itet) == ic) ip1 = i
        end do  

! bearbeite alle Kanten, die von Eckpunkt nump1 ausgehen
! d.h. die Kanten zu den anderen Eckpunkten 
        do i=1,4
          nump1 = ip1
          if (i == nump1) cycle
          nump2 = i
          noedge = kanten_nummer(nump1,nump2)
! die Kante ist schon bei einer anderen Koordinate gefunden worden 
            if (tetra_kanten(noedge,itet) /= 0) cycle
! die Kante ist neu
          jc = nteck(nump2,itet)
            nkanten = nkanten+1
          kanten(1,nkanten) = ic
          kanten(2,nkanten) = jc
          tetra_kanten(noedge,itet) = nkanten

! durchlaufe alle Nachbartetraeder und markiere die gemeinsame Kante
        
          akt_tet = itet
          akt_edge = noedge
            no_side = angrenzende_seiten(1,akt_edge)
          isi = 0
          do
            nxt_tet = ntbar(no_side,akt_tet)
          if (nxt_tet == 0) then
            isi = isi+1
! an beiden Seiten bis zum Rand gelaufen
            if (isi > 1) exit
                no_side = angrenzende_seiten(2,akt_edge)
            nxt_tet = ntbar(no_side,itet)
            if (nxt_tet == 0) exit
            akt_tet = itet
          end if
            nxt_side = ntseite(no_side,akt_tet)
          if (nxt_tet == itet) exit
! bestimme die Kantennummer auf dem neuen Tetraeder
          do j=1,4
            if (nteck(j,nxt_tet) == ic) nump1 = j
            if (nteck(j,nxt_tet) == jc) nump2 = j
          end do
          noedge = kanten_nummer(nump1,nump2)
            tetra_kanten(noedge,nxt_tet) = nkanten
          no_side = angrenzende_seiten(1,noedge) + 
     .              angrenzende_seiten(2,noedge) - nxt_side
          akt_tet = nxt_tet
          end do  
                
        end do  

          cur => cur%next_tet

        end do  
      
      end do

      write (0,*) 'min(tetra_kanten) ',minval(tetra_kanten(:,1:ntet))
      return
      end subroutine suche_kanten
      

      subroutine suche_kanten2

      implicit none
      integer :: itet, j, nxt_side, akt_tet, akt_edge, isi, ik
      integer :: ic, i, nump1, nump2, noedge, jc, no_side, nxt_tet

      integer :: kpunkte(2,6) = reshape (
     .           (/ 1,2,  1,3,  1,4,  2,3,  2,4,  3,4 /), (/ 2, 6 /) )

      nkanten = 0
! für jeden Tetraeder
      do itet=1,ntet

! durchsuche die Kanten des Tetraeders
        do ik=1,6
          
! die Kante ist schon bei einem anderen Tetraeder gefunden worden 
          if (tetra_kanten(ik,itet) /= 0) cycle

! die Kante ist neu
          nump1 = kpunkte(1,ik)
          nump2 = kpunkte(2,ik)
          ic = nteck(nump1,itet)
          jc = nteck(nump2,itet)

          nkanten = nkanten+1
          kanten(1,nkanten) = ic
        kanten(2,nkanten) = jc
        tetra_kanten(ik,itet) = nkanten
      
! durchlaufe alle Nachbartetraeder und markiere die gemeinsame Kante
        
          noedge = ik
        akt_tet = itet
        akt_edge = noedge
          no_side = angrenzende_seiten(1,akt_edge)
        isi = 0
        do
          nxt_tet = ntbar(no_side,akt_tet)
          if (nxt_tet == 0) then
            isi = isi+1
! an beiden Seiten bis zum Rand gelaufen
          if (isi > 1) exit
              no_side = angrenzende_seiten(2,akt_edge)
            nxt_tet = ntbar(no_side,itet)
            if (nxt_tet == 0) exit
            akt_tet = itet
          end if
          nxt_side = ntseite(no_side,akt_tet)
          if (nxt_tet == itet) exit
! bestimme die Kantennummer auf dem neuen Tetraeder
          do j=1,4
            if (nteck(j,nxt_tet) == ic) nump1 = j
            if (nteck(j,nxt_tet) == jc) nump2 = j
          end do
          noedge = kanten_nummer(nump1,nump2)
          tetra_kanten(noedge,nxt_tet) = nkanten
          no_side = angrenzende_seiten(1,noedge) + 
     .              angrenzende_seiten(2,noedge) - nxt_side
          akt_tet = nxt_tet
        end do  

        end do  ! ik
      
      end do  ! itet
      return
      end subroutine suche_kanten2


      subroutine schneide_kanten
      implicit none
      integer :: ik, k1, k2, ic
      real(dp) :: p1(3), p2(3), v(3), a(3)
      real(dp) :: d, xnen, t, x, y, z, dst
      logical :: inserted

      type(TAVLTree), pointer :: baum

      allocate(spar(nkanten))
      allocate(iedge(nkanten))
      allocate(icono(nkanten))

      spar = 1.E30_dp
      iedge = 0
      icono = 0
      nrknot = 0
      baum => NewTree()

      a(1:3) = cutplane(2:4)
      d = cutplane(1)

      do ik=1,nkanten

        k1 = kanten(1,ik)
        k2 = kanten(2,ik)
      
        p1(1) = xtetra(k1)
        p1(2) = ytetra(k1)
        p1(3) = ztetra(k1)

        p2(1) = xtetra(k2)
        p2(2) = ytetra(k2)
        p2(3) = ztetra(k2)

        v = p2 - p1 

        xnen = sum(a*v)
        if (abs(xnen) < eps10) cycle

        t = -(d + sum(a*p1)) / xnen
        
        if((t > -eps6) .and. (t < 1._dp+eps6)) then  

          x = xtetra(k1) + t*(xtetra(k2)-xtetra(k1)) 
          y = ytetra(k1) + t*(ytetra(k2)-ytetra(k1)) 
          z = ztetra(k1) + t*(ztetra(k2)-ztetra(k1)) 
          dst = sqrt( (xtetra(k2)-xtetra(k1))**2 +
     .                (ytetra(k2)-ytetra(k1))**2 +
     .                (ztetra(k2)-ztetra(k1))**2 )

          ic = nrknot + 1
          inserted=.false.
          call insert (baum, x, y, z, dst, ic, inserted)

          if (inserted) then
            nrknot = nrknot + 1
            iedge(nrknot) = ik
          end if

          spar(ik) = t
          icono(ik) = ic
        end if  

      end do

      call DestroyTree(baum)

      end subroutine schneide_kanten


      subroutine berechne_koordinaten

      implicit none
      real(dp) :: a(3), b(3), c(3), p(3), orig(3)
      real(dp) :: am(2,2), amm1(2,2), rhs(2), x(3)
      real(dp) :: spat, bnorm, cnorm, detam
      integer :: i, ik, k1, k2

! a ist die richtung der normalen zur schnittebene
      a(1:3) = cutplane(2:4)

! berechne vektor b, b senkrecht zu a
! d.h.  a1*b1 + a2*b2 + a3*b3 = 0

      if (abs(a(1)) > eps10) then
        b(1) = -(a(2)+a(3))/a(1)
        b(2) = 1._dp
        b(3) = 1._dp
      else if (abs(a(2)) > eps10) then
        b(1) = 1._dp
        b(2) = -(a(1)+a(3))/a(2)
        b(3) = 1._dp
      else if (abs(a(3)) > eps10) then
        b(1) = 1._dp
        b(2) = 1._dp
        b(3) = -(a(1)+a(2))/a(3)
      else 
        write (6,*) ' Problem in berechne_koordinaten '
        write (6,*) ' die Koeffizienten der Schnittebene sind 0 '
        write (6,*) cutplane
        call exit (1) 
      end if 

! berechne vektor c = a x b ==> c ist senkrecht zu a und b

      c(1) = a(2)*b(3) - b(2)*a(3)
      c(2) = a(3)*b(1) - b(3)*a(1)
      c(3) = a(1)*b(2) - b(1)*a(2)
      
! berechne spatprodukt, (a x b)c > 0 ==> rechtsystem, sonst c=-c

      spat =  a(1)*b(2)*b(3) + a(2)*b(3)*c(1) + a(3)*b(1)*c(2)
     .      - c(1)*b(2)*a(3) - c(2)*b(3)*a(1) - c(3)*b(1)*a(2)
      if (spat < 0._dp) c = -c
      
! normiere b und c
      bnorm = sqrt(sum(b*b))
      cnorm = sqrt(sum(c*c))  
      
      b = b / bnorm
      c = c / cnorm
      
! b und c sind eine orthonormalbasis der schnittebene 

              
      nknot = nrknot
      nknots = nknot
      ntri = ntet / 2
      ntris = ntri

      call dealloc_ctrig
      call alloc_ctrig

      ik = iedge(1)
      k1 = kanten(1,ik)
      k2 = kanten(2,ik)
      orig(1) = xtetra(k1) + spar(ik)*(xtetra(k2)-xtetra(k1))
      orig(2) = ytetra(k1) + spar(ik)*(ytetra(k2)-ytetra(k1))
      orig(3) = ztetra(k1) + spar(ik)*(ztetra(k2)-ztetra(k1))

      am(1,1) = sum(b*b)
      am(1,2) = sum(b*c)
      am(2,1) = sum(b*c)
      am(2,2) = sum(c*c)

      detam = am(1,1)*am(2,2) - am(1,2)*am(2,1)

      if (abs(detam) < eps10) then
        write (6,*) ' problem in berechne_koordinaten '
        write (6,*) ' gls zur koordinatentransformation nicht loesbar'
        call exit (1)
      end if

      amm1(1,1) = am(2,2) / detam
      amm1(1,2) = -am(1,2) / detam
      amm1(2,1) = -am(2,1) / detam
      amm1(2,2) = am(1,1) / detam

      do i=1,nrknot

        ik = iedge(i)
        k1 = kanten(1,ik)
        k2 = kanten(2,ik)

        x(1) = xtetra(k1) + spar(ik)*(xtetra(k2)-xtetra(k1)) 
        x(2) = ytetra(k1) + spar(ik)*(ytetra(k2)-ytetra(k1)) 
        x(3) = ztetra(k1) + spar(ik)*(ztetra(k2)-ztetra(k1)) 

        rhs(1) = sum(b*(x-orig))
        rhs(2) = sum(c*(x-orig))

        xtrian(i) = amm1(1,1)*rhs(1) + amm1(1,2)*rhs(2)
        ytrian(i) = amm1(2,1)*rhs(1) + amm1(2,2)*rhs(2)
      end do

      end subroutine berechne_koordinaten



      subroutine sort_ueberpruefen(ipunkt)
  
      implicit none

      integer, intent(inout), dimension(4) :: ipunkt
      real(dp), dimension(4,2)             :: punkt
      real(dp), dimension(2,2)             :: a
      real(dp), dimension(2)               :: b, t, cen
      real(dp)                             :: x1, x2
      real(dp)                             :: d1, d2, d3, d4, d, deta
      logical                              :: lsw1, lsw2
      integer                              :: ih, i
         
      punkt(1:4,1) = (/ (xtrian(ipunkt(i)), i=1,4) /)
      punkt(1:4,2) = (/ (ytrian(ipunkt(i)), i=1,4) /)

      cen(1)=sum(punkt(1:4,1))*0.25_dp
      cen(2)=sum(punkt(1:4,2))*0.25_dp

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
        a(1,2) = punkt(3,1) - punkt(4,1)
        a(2,2) = punkt(3,2) - punkt(4,2)

        b(1) = punkt(3,1) - punkt(1,1)
        b(2) = punkt(3,2) - punkt(1,2)
            
        a = a*d
        b = b*d
            
        deta = a(1,1)*a(2,2) - a(2,1)*a(1,2)
        if (abs(deta) < eps10) then
! p1-p2 und p3-p4 sind parallel
        lsw1 = .false.
        else
        x1 = (b(1)*a(2,2) - b(2)*a(1,2))/deta
        x2 = (b(2)*a(1,1) - b(1)*a(2,1))/deta
          if ((x1 >= 0._dp) .and. (x1 <= 1._dp) .and.
     .        (x2 >= 0._dp) .and. (x2 <= 1._dp)) then
!  intersection found, switch points 2 and 3
            t(1:2) = punkt(2,1:2)
            punkt(2,1:2) = punkt(3,1:2)
            punkt(3,1:2) = t(1:2)
          ih = ipunkt(2)
            ipunkt(2) = ipunkt(3)
            ipunkt(3) = ih
            lsw1 = .true.
          else
            lsw1 = .false.
          end if
        end if

!  test p1-p4 and p2-p3
        a(1,1) = punkt(4,1) - punkt(1,1)
        a(2,1) = punkt(4,2) - punkt(1,2)
        a(1,2) = punkt(2,1) - punkt(3,1)
        a(2,2) = punkt(2,2) - punkt(3,2)

        b(1) = punkt(2,1) - punkt(1,1)
        b(2) = punkt(2,2) - punkt(1,2)
            
        a = a*d
        b = b*d
            
        deta = a(1,1)*a(2,2) - a(2,1)*a(1,2)
        if (abs(deta) < eps10) then
! p1-p4 und p2-p3 sind parallel
        lsw2 = .false.
        else
        x1 = (b(1)*a(2,2) - b(2)*a(1,2))/deta
        x2 = (b(2)*a(1,1) - b(1)*a(2,1))/deta
      
          if ((x1 >= 0._dp) .and. (x1 <= 1._dp) .and.
     .        (x2 >= 0._dp) .and. (x2 <= 1._dp)) then
!  intersection found, switch points 3 and 4
            t(1:2) = punkt(3,1:2)
            punkt(3,1:2) = punkt(4,1:2)
            punkt(4,1:2) = t(1:2)
            ih = ipunkt(3)
            ipunkt(3) = ipunkt(4)
            ipunkt(4) = ih
            lsw2 = .true.
          else
            lsw2 = .false.
          end if
        end if  

      end do                 ! do while

      end subroutine sort_ueberpruefen



      subroutine bilde_dreiecke

      implicit none
      integer :: iknot(6), ipunkt(4)
      integer :: itet, icount, j, ih, i1, i2, i3, ico
      real(dp) :: ar, artri3

      if (.not.allocated(itetno)) allocate (itetno(ntri))
      itetno = 0

      ntrii = 0
      do itet=1,ntet
!        iknot(1:6) = (/ (icono(tetra_kanten(j,itet)),j=1,6) /)
!        icount = count(iknot > 0)
        icount=0
        iloop: do i=1,6
          ico = icono(tetra_kanten(i,itet))
          if (ico > 0) then
            do j=1,icount
              if (ico == iknot(j)) cycle iloop
            end do
            icount = icount + 1
            iknot(icount) = ico
          end if
        end do iloop

        if (icount == 3) then
! schnittgebilde ist ein Dreieck
          ntrii = ntrii + 1
          necke(1:3,ntrii) = iknot(1:3)
          i1=necke(1,ntrii)
          i2=necke(2,ntrii)
          i3=necke(3,ntrii)
!          isum=i1+i2+i3
!          imn=min(i1,i2,i3)
!          imx=max(i1,i2,i3)
!          imid=isum-imn-imx
!          if ((imn == imid) .or. (imx == imid)) then
!            ntrii = ntrii - 1
!            cycle
!          end if
          itetno(ntrii) = itet
          ar = artri3(xtrian(i1),ytrian(i1),0._dp,
     .                xtrian(i2),ytrian(i2),0._dp,
     .                xtrian(i3),ytrian(i3),0._dp)
          if (ar < 0._dp) then
! orientierung des dreiecks stimmt nicht, drehe die punkte 2 und 3
            ih = necke(2,ntrii)
            necke(2,ntrii) = necke(3,ntrii)
            necke(3,ntrii) = ih
          end if

        else if (icount ==4) then
! schnittgebilde ist ein Viereck
          ipunkt = iknot(1:4)

! sorge dafuer, dass die Linie p1,p2,p3,p4 ein Viereck bildet und
! sich nicht schneidet
          call sort_ueberpruefen (ipunkt)
! bilde erstes Dreieck p1, p2, p3
          ntrii = ntrii + 1
          necke(1:3,ntrii) = (/ ipunkt(1), ipunkt(2), ipunkt(3) /)   
          i1=necke(1,ntrii)
          i2=necke(2,ntrii)
          i3=necke(3,ntrii)
!          isum=i1+i2+i3
!          imn=min(i1,i2,i3)
!          imx=max(i1,i2,i3)
!          imid=isum-imn-imx
!          if ((imn < imid) .and. (imx > imid)) then
            itetno(ntrii) = itet
            ar = artri3(xtrian(i1),ytrian(i1),0._dp,
     .                  xtrian(i2),ytrian(i2),0._dp,
     .                  xtrian(i3),ytrian(i3),0._dp)
            if (ar < 0._dp) then
! orientierung des dreiecks stimmt nicht, drehe die punkte 2 und 3
              ih = necke(2,ntrii)
              necke(2,ntrii) = necke(3,ntrii)
              necke(3,ntrii) = ih
            end if
!          else
!            ntrii = ntrii-1
!          end if

! bilde zweites Dreieck p1, p3, p4
          ntrii = ntrii + 1
          necke(1:3,ntrii) = (/ ipunkt(1), ipunkt(3), ipunkt(4) /)   
          i1=necke(1,ntrii)
          i2=necke(2,ntrii)
          i3=necke(3,ntrii)
!          isum=i1+i2+i3
!          imn=min(i1,i2,i3)
!          imx=max(i1,i2,i3)
!          imid=isum-imn-imx
!          if ((imn == imid) .or. (imx == imid)) then
!            ntrii = ntrii - 1
!            cycle
!          end if
          itetno(ntrii) = itet
          ar = artri3(xtrian(i1),ytrian(i1),0._dp,
     .                xtrian(i2),ytrian(i2),0._dp,
     .                xtrian(i3),ytrian(i3),0._dp)
          if (ar < 0._dp) then
! orientierung des dreiecks stimmt nicht, drehe die punkte 2 und 3
            ih = necke(2,ntrii)
            necke(2,ntrii) = necke(3,ntrii)
            necke(3,ntrii) = ih
        end if
         
        end if       
      end do

      do itri=1,ntrii
        xcom(itri) = (xtrian(necke(1,itri)) + xtrian(necke(2,itri)) +
     .                xtrian(necke(3,itri))) / 3._dp
        ycom(itri) = (ytrian(necke(1,itri)) + ytrian(necke(2,itri)) +
     .                ytrian(necke(3,itri))) / 3._dp
      end do

      end subroutine bilde_dreiecke


      end subroutine rpscut
