      module EIRMOD_phot_reflec
      use EIRMOD_precision
      use EIRMOD_parmmod
      use EIRMOD_cinit
      use EIRMOD_comprt, only: iunout
      implicit none
 
      private
 
      private :: EIRENE_R, EIRENE_interpolate
      public :: EIRENE_init_refl_hlm, EIRENE_reflect_hollmann
 
      real(dp), save :: PI = 3.141592654_dp
      real(dp), save, dimension(9, 16) :: graphite
      real(dp), save, dimension(9, 16) :: mo
      real(dp), external :: ranf_eirene
 
      contains
 
      subroutine EIRENE_init_refl_hlm()
      integer :: i, ifile
! Fehler abfragen...
!pb      open(unit = 23, file = "graphite_ext.dat")
      DO IFILE=1, NDBNAMES
        IF (INDEX(DBHANDLE(IFILE),'gr_ext') /= 0) EXIT
      END DO
 
      IF (IFILE > NDBNAMES) THEN
        WRITE (IUNOUT,*) ' NO DATABASENAME FOR graphite.ext.dat DEFINED'
        WRITE (IUNOUT,*) ' CALCULATION ABANDONNED '
        CALL EIRENE_EXIT_OWN(1)
      END IF
 
      OPEN (UNIT=23,FILE=DBFNAME(IFILE))
 
      do i = 1, 9
        read(unit = 23, fmt = *)
     .   graphite(i, 1),graphite(i, 2),graphite(i, 3),graphite(i, 4),
     .   graphite(i, 5),graphite(i, 6),graphite(i, 7),graphite(i, 8),
     .   graphite(i, 9),graphite(i, 10),graphite(i, 11),graphite(i, 12),
     .   graphite(i, 13),graphite(i, 14),graphite(i, 15),graphite(i, 16)
      end do
      close(unit = 23)
 
!pb      open(unit = 23, file = "mo_ext.dat")
      DO IFILE=1, NDBNAMES
        IF (INDEX(DBHANDLE(IFILE),'mo_ext') /= 0) EXIT
      END DO
 
      IF (IFILE > NDBNAMES) THEN
        WRITE (IUNOUT,*) ' NO DATABASENAME FOR mo.ext.dat DEFINED'
        WRITE (IUNOUT,*) ' CALCULATION ABANDONNED '
        CALL EIRENE_EXIT_OWN(1)
      END IF
 
      OPEN (UNIT=23,FILE=DBFNAME(IFILE))
      do i = 1, 9
        read(unit = 23, fmt = *)
     .   mo(i, 1),mo(i, 2),mo(i, 3),mo(i, 4),mo(i, 5),mo(i, 6),
     .   mo(i, 7),mo(i, 8),mo(i, 9),mo(i, 10),mo(i, 11),mo(i, 12),
     .   mo(i, 13),mo(i, 14),mo(i, 15),mo(i, 16)
      end do
      close(unit = 23)
 
 
      end subroutine EIRENE_init_refl_hlm
 
      subroutine EIRENE_interpolate(theta_i, lambda_in, mat, theta_0,
     .  k, rho_0)
      real(dp), intent(in) :: theta_i, lambda_in
      integer, intent(in) :: mat
      real(dp), intent(out) :: theta_0, k, rho_0
      integer :: column, line = 1
      real(dp) :: k0, k1, t0, t1, r0, r1, theta, lambda, lam_min
      real(dp), dimension(9,16) :: tab
 
! welches Material ?
      if (mat == 1) then
         tab = graphite
      else if (mat == 2) then
         tab = mo
      end if
 
! zwischen welchen Spalten muss gesucht werden ?
      if (theta_i <= 30.0) then
         column = 2
         theta = 30.0
      else if (theta_i <= 45.0) then
         column = 5
         theta = 45.0
      else if (theta_i <= 60.0) then
         column = 8
         theta = 60.0
      else
         column = 11
         theta = 75.0
      end if
 
      lam_min = tab(1,1) + (tab(2,1)-tab(1,1))*1.e-3_dp
      lambda = max(lambda_in,lam_min)
      lambda = min(lambda,tab(9,1))
 
! Zeile suchen
! bzgl. lambda
      do while(lambda > tab(line, 1))
         line = line + 1
      end do
 
! Interpolieren...
      k0 = tab(line, column)-(tab(line,1)-lambda)*
     .     (tab(line,column)-tab(line-1,column))/
     .     (tab(line,1)-tab(line-1,1))
      k1 = tab(line, column+3)-(tab(line,1)-lambda)*
     .     (tab(line,column+3)-tab(line-1,column+3))/
     .     (tab(line,1)-tab(line-1,1))
 
      t0 = tab(line, column+1)-(tab(line,1)-lambda)*
     .     (tab(line,column+1)-tab(line-1,column+1))/
     .     (tab(line,1)-tab(line-1,1))
      t1 = tab(line, column+4)-(tab(line,1)-lambda)*
     .     (tab(line,column+4)-tab(line-1,column+4))/
     .     (tab(line,1)-tab(line-1,1))
 
      r0 = tab(line, column+2)-(tab(line,1)-lambda)*
     .     (tab(line,column+2)-tab(line-1,column+2))/
     .     (tab(line,1)-tab(line-1,1))
      r1 = tab(line, column+5)-(tab(line,1)-lambda)*
     .     (tab(line,column+5)-tab(line-1,column+5))/
     .     (tab(line,1)-tab(line-1,1))
 
      theta_0 = t1 - (theta-theta_i) * (t1-t0)/(15.0)
      k = k1 - (theta-theta_i) * (k1-k0)/(15.0)
      rho_0 = r1 - (theta-theta_i) * (r1-r0)/(15.0)
 
      end subroutine EIRENE_interpolate
 
! reflect
! Die Routine entscheidet, ob ein Photon mit Einfallwinkel theta_i und
! Wellenlaenge lambda_i reflektiert wird. Im Falle einer
! Reflektion ist flag 1, ansonsten 0.
      subroutine EIRENE_reflect_hollmann(theta_i, lambda_i, mat, flag,
     .                   theta_out, alpha_out, rprob)
      real(dp), intent(in) :: theta_i, lambda_i
      integer, intent(in) :: mat
      real(dp), intent(out) :: theta_out, alpha_out, rprob
      integer, intent(out) :: flag
 
      real(dp) :: theta_kor, R_kor, alpha
! Werte durch Interpolieren der Tabelle
      real(dp) :: k, theta_0, rho_0
! Radien innerhalb des Lobes
      real(dp) :: rnd, rad_1, rad_2, rho
 
! "runde" Einfallwinkel auf Winkel mit bekannten Parametern
      if (theta_i <= 22.5) then
         theta_kor = 15.0
      else if (theta_i <= 37.5) then
         theta_kor = 30.0
      else if (theta_i <= 52.5) then
         theta_kor = 45.0
      else if (theta_i <= 67.5) then
         theta_kor = 60.0
      else
         theta_kor = 75.0
      end if
 
      call EIRENE_interpolate(theta_kor, lambda_i, mat, theta_0, k,
     .  rho_0)
 
! berechne R des korrigierten Einfallwinkels
      call EIRENE_R(theta_kor, lambda_i, mat, R_kor)
      theta_kor = theta_kor * PI / 180.0_dp
      theta_0 = theta_0 * PI / 180.0_dp
      rprob = R_kor
 
 
!pb      call random_number(rnd)
      rnd=ranf_eirene()
! wird das Photon reflektiert ?
!pb      if (rnd < R_kor) then
!pb  berechne die Winkel immer, Entscheidung, ob reflektiert faellt in
!pb  reflect_photon
! ja:
! bestimme Reflektionswinkel, welche das Photon nicht ins Material reflektieren
         do
!pb            call random_number(rnd)
            rnd=ranf_eirene()
! theta_out = g^(-1) (rnd)
            theta_out = theta_0 + k*PI/2.0 - k/2.0* acos(2.0*rnd - 1.0)
!pb            call random_number(rnd)
            rnd=ranf_eirene()
! innerer Winkel
            alpha = rnd * 2.0 * PI
 
            rho = rho_0*cos((theta_out - theta_0)/k)
            rad_1 = sin(theta_out - theta_0) * rho
            rad_2 = cos(theta_out - theta_0) * rho
 
! beide Winkel ok, bzw. wird das Photon nicht ins Material reflektiert ?
! (z-Koordinate im kartesischem System positiv)
            if ((sin(theta_0)*sin(alpha)*rad_1+cos(theta_0)*rad_2)
     .           > 0.0) then
               exit
            end if
         end do
 
! Winkel alpha noch bestimmen bzgl x/y-Achse
! tan alpha = x-Koord./y-Koord.
         alpha_out = atan( (cos(alpha)*rad_1)/(cos(theta_0)*sin(alpha)*
     .               rad_1-sin(theta_0)*rad_2) )
         flag = 1
 
! x-Koord.
! cos(alpha_out)*rad_1
! y-Koord.
! cos(theta_0)*sin(alpha)*rad_1-sin(theta_0)*rad_2
! z-Koord.
! sin(theta_0)*sin(alpha)*rad_1+cos(theta_0)*rad_2
 
! nein:
!pb      else
!pb         flag = 0
!pb      end if
 
      end subroutine EIRENE_reflect_hollmann
 
! r
! berechnet Wahrscheinlichkeit einer Reflektion eines Photon mit
! Wellenlaenge  lambda Angstrom, welches im Winkel theta_i auftrifft
! (die Winkel im Winkelmass)
      subroutine EIRENE_R(theta_i, lambda, mat,  out)
      real(dp), intent(in) :: theta_i, lambda
      integer, intent(in) :: mat
      real(dp), intent(out) :: out
      real(dp) :: theta_0, k, rho_0, theta
      real(dp) :: n_0, n_1
 
      real(dp) :: rnd_alpha, rnd_theta, rnd, rho
      real(dp) :: rad_1, rad_2
      real(dp) :: area
! MonteCarlo Variablen
      integer :: i, N = 100, treffer
 
      out = 0._dp
 
      if (.not.(theta_i < 15.0 .or. theta_i > 75.0 .or.
     .     lambda < 3341.0 .or. lambda > 7729.0)) then
 
! fuer den einfachen Fall werden nur k und rho_0
! benoetigt
         call EIRENE_interpolate(theta_i, lambda, mat, theta_0, k,
     .  rho_0)
! umrechnen ins Bogenmass
         theta_0 = theta_0 * PI / 180.0;
         theta = theta_i * PI / 180.0;
! berechne linke und rechte Nullstelle von theta_0
         n_0 = -PI*k/2.0 + theta_0;
         n_1 = PI*k/2.0 + theta_0;
 
! Schnitt mit der Pi/2 - Achse ?
! nein:
         if (n_1 <= PI/2.0) then
            out = 2.0*PI*k*rho_0*((sin(k*Pi/2.0)-k)/(1.0-k**2.0))
! ja:
         else
            treffer = 0
 
! Oberflaeche berechnen
            area = 2.0*PI*k*rho_0*((sin(k*Pi/2.0)-k)/(1.0-k**2.0))
 
! n-mal wuerfeln
            do i=1, N
! Acceptance-Rejection
               do
!pb                  call random_number(rnd)
                  rnd = ranf_eirene()
                  rnd_theta = n_0 + rnd * (n_1-n_0)
                  rho = rho_0 * cos((rnd_theta - theta_0)/k)
!pb                  call random_number(rnd)
                  rnd = ranf_eirene()
                  if ((abs(sin(rnd_theta-theta_0)*rho)/rho_0) >= rnd)
     .                 then
                     exit
                  end if
               end do
 
! wuerfeln
!pb               call random_number(rnd)
               rnd = ranf_eirene()
               rnd_alpha = rnd * 2.0 * PI
 
               rad_1 = sin(rnd_theta - theta_0) * rho
               rad_2 = cos(rnd_theta - theta_0) * rho
 
               if ((sin(theta_0)*sin(rnd_alpha)*rad_1+
     .              cos(theta_0)*rad_2) > 0.0) then
                  treffer = treffer + 1
               end if
            end do
            out = area * treffer / N
         end if
      end if
 
      end subroutine EIRENE_R
 
      end module EIRMOD_phot_reflec
 
 
 
 
