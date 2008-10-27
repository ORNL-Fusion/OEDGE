C EIRENE04 COMPILATION
C ===== SOURCE: phot_reflec.f
      module phot_reflec
      use precision
      implicit none

      private

      private :: R, interpolate
      public :: init_refl_hlm, reflect_hollmann

      real(dp), save :: PI = 3.141592654_dp
      real(dp), save, dimension(9, 16) :: graphite
      real(dp), save, dimension(9, 16) :: mo  
      real(dp), external :: ranf_eirene

      contains
 
      subroutine init_refl_hlm()
      integer :: i
! Fehler abfragen...
      open(unit = 23, file = "graphite_ext.dat")
      do i = 1, 9
        read(unit = 23, fmt = *)  
     .   graphite(i, 1),graphite(i, 2),graphite(i, 3),graphite(i, 4),
     .   graphite(i, 5),graphite(i, 6),graphite(i, 7),graphite(i, 8),
     .   graphite(i, 9),graphite(i, 10),graphite(i, 11),graphite(i, 12),
     .   graphite(i, 13),graphite(i, 14),graphite(i, 15),graphite(i, 16)
      end do
      close(unit = 23)

      open(unit = 23, file = "mo_ext.dat")
      do i = 1, 9
        read(unit = 23, fmt = *)  
     .   mo(i, 1),mo(i, 2),mo(i, 3),mo(i, 4),mo(i, 5),mo(i, 6),
     .   mo(i, 7),mo(i, 8),mo(i, 9),mo(i, 10),mo(i, 11),mo(i, 12),
     .   mo(i, 13),mo(i, 14),mo(i, 15),mo(i, 16)
      end do
      close(unit = 23)


      end subroutine init_refl_hlm

      subroutine interpolate(theta_i, lambda_in, mat, theta_0, k, rho_0)
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
      
      end subroutine interpolate
  
! reflect
! Die Routine entscheidet, ob ein Photon mit Einfallwinkel theta_i und
! Wellenlaenge lambda_i reflektiert wird. Im Falle einer 
! Reflektion ist flag 1, ansonsten 0. 
      subroutine reflect_hollmann(theta_i, lambda_i, mat, flag, 
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
      
      call interpolate(theta_kor, lambda_i, mat, theta_0, k, rho_0)

! berechne R des korrigierten Einfallwinkels 
      call R(theta_kor, lambda_i, mat, R_kor)
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
  
      end subroutine reflect_hollmann

! r
! berechnet Wahrscheinlichkeit einer Reflektion eines Photon mit 
! Wellenlaenge  lambda Angstrom, welches im Winkel theta_i auftrifft
! (die Winkel im Winkelmass)
      subroutine R(theta_i, lambda, mat,  out)
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
         call interpolate(theta_i, lambda, mat, theta_0, k, rho_0)
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

      end subroutine R

      end module phot_reflec




C ===== SOURCE: reflec_photon.f
C
C
      SUBROUTINE REFLEC_photon
C
C  REFLECT ESCAPING PHOTONS
C  INPUT:
C       ILREF = 0  PERFECT ABSORPTION
C       ILREF = 1  ERIC HOLLMAN DATABASE
C
C       ITYP  = 0  INCIDENT PHOTON
C  OUTPUT:
C     LGPART= TRUE AND:
C       ITYP = 0  PHOTON IPHOT IS RETURNED TO CALLING PROGRAM
C     LGPART= FALSE  NO PARTICLE IS RETURNED (ABSORBTION)
C       ITYP = -1
C
      USE PRECISION
      USE PARMMOD
      USE PHOT_REFLEC
      USE CTRCEI
      USE CADGEO
      USE CLGIN
      USE COMUSR
      USE CLOGAU
      USE CESTIM
      USE COMPRT
      USE CRAND
      USE CREF
      USE CCONA
      USE PHOTON
      USE CPES

      IMPLICIT NONE
    
      REAL(DP), INTENT(IN) :: WMIN, XMP, XCP
      INTEGER, INTENT(IN) :: NPRIN
      INTEGER, INTENT(INOUT) :: IGASF, IGAST
      INTEGER, SAVE :: IFIRST=0, NPANOLD=0
      INTEGER :: ILIM, ISP, ICOUNT, ISTS, ISEE, MODREF, ISPZO,
     .           IMAT, IREFL, MSS
      REAL(DP) :: DUMMY, XMW, XCW, E0TERM, EBIND, PRFCF, PRFCT,
     .            EXPP, EXPE, EXPI, RINTG, EINTG, AINTG, COSIN,
     .            THETA, XLAMBDA, THETA_OUT, ALPHA_OUT, WFAC, 
     .            FR1, ZCPHI, ZSPHI, ZCTHET, ZSTHET, VX, VY, VZ,
     .            RPROB,WABS
      REAL(DP), EXTERNAL :: RANF_EIRENE, ranset_eirene
      INTEGER, EXTERNAL :: RANGET_EIRENE, IDEZ
C
C---------------------------------------------------------------------
C

C
C  INITIALIZE SURFACE REFLECTION MODELS FOR PHOTONS
C
      ENTRY REFLC0_PHOTON
C
      IF (IFIRST.EQ.1) RETURN
      IFIRST=1
C
      CALL INIT_REFL_HLM
C
      IF (TRCREF) THEN
        CALL LEER(2)
      ENDIF
C
C
C  PRINTOUT REFLECTION PROPERTIES OF SURFACES
C  ALREADY DONE FROM REFLC0
C
C
      RETURN
C
      ENTRY REFLC1_PHOTON (WMIN,XMP,XCP,NPRIN,IGASF,IGAST)
C
C  SYNCHRONIZE RANDOM NUMBERS
C
      IF (NLCRR.AND.(NPANU.NE.NPANOLD)) THEN
C  INITIALIZE RANDOM NUMBERS FOR EACH PARTICLE, TO GENERATE CORRELATION
C       Call RANSET_EIRENE(ISEED)
        dummy=ranset_eirene(iseedR)
        DUMMY=RANF_EIRENE( )
        ISEEDR=ranget_eirene(isee)
        ISEEDR=INTMAX-ISEEDR
        NPANOLD=NPANU
      END IF
C
C  SURFACE NUMBER  : MSURF (MSURF=0: DEFAULT MODEL)
C  SPECIES INDEX   : ISPZ
C
      MODREF=IDEZ(ILREF(MSURF),2,2)
      XMW=ZNML(MSURF)
      XCW=ZNCL(MSURF)
      E0TERM=EWALL(MSURF)
      EBIND=EWBIN(MSURF)
      PRFCF=RECYCF(ISPZ,MSURF)
      PRFCT=RECYCT(ISPZ,MSURF)
      EXPP=EXPPL(ISPZ,MSURF)
      EXPE=EXPEL(ISPZ,MSURF)
      EXPI=EXPIL(ISPZ,MSURF)
      RINTG=RINTEG(MSURF)
      EINTG=EINTEG(MSURF)
      AINTG=AINTEG(MSURF)
      ISPZO=ISPZ
C
C
C   TENTATIVELY ASSUME  REFLECTION
      LGPART=.TRUE.
C   COSINE OF ANGLE OF INCIDENCE
      COSIN=VELX*CRTX+VELY*CRTY+VELZ*CRTZ
      IF (COSIN.LT.0.D0) GOTO 993
C
C
C
C   MODREF=0: "PERFECTLY ABSORBING SURFACE, DEFAULT
C   MODREF=1: "DATABASE REFLECTION MODEL" 
C
      IF (MODREF.EQ.1.AND.PRFCF.GT.0.) THEN
        GOTO 100
      ELSE
C  ABSORB THIS PHOTON
        GOTO 700
      ENDIF
C
C  DATABASE REFLECTION MODEL STARTS HERE
C
100   CONTINUE
C
C   CHECK IF WALL REFLECTION DATA FOR IPHOT INCIDENT ON
C   XWALL/ZWALL ARE AVAILABLE
C
C  CARBON OR MOLYBDENUM ?

      IMAT = 2
      
C
      THETA = ACOS(COSIN)*RADDEG
      XLAMBDA = hpcl/E0*10._DP*1.E7_DP

      CALL REFLECT_HOLLMANN (THETA, XLAMBDA, IMAT, IREFL, THETA_OUT,
     .                       ALPHA_OUT, RPROB)
C
130   CONTINUE
C
C
C   DECIDE IF PARTICLE IS TO BE REFLECTED OR ABSORBED
C   (NO THERMAL RE-EMISSION MODEL FOR INCIDENT PHOTONS)
C
      IF (WEIGHT.GT.WMIN) THEN
C  WITH SUPPRESSION OF ABSORPTION
        WABS=WEIGHT*(1.D0-RPROB)
        IF (WABS.GT.0.D0) THEN
          IF (LSPUMP) SPUMP(ISPZ,MSURF)=SPUMP(ISPZ,MSURF)+WABS
        ENDIF  
        WEIGHT=WEIGHT-WABS
        IF (WEIGHT.LE.EPS30) GOTO 700
      ELSE
C  NO SUPPRESSION OF ABSORPTION
        FR1=RANF_EIRENE( )
        IF (FR1.GT.RPROB) GOTO 700
      ENDIF
C
C  SPECIES OF REFLECTED PARTICLE
      IF (IGASF.LT.1.OR.IGASF.GT.NPHOTI) GOTO 992
      IPHOT=IGASF
      ISPZ=IPHOT
      ITYP=0
C
C  ENERGIE (WAVELENGTH):  NOT MODIFIED
C
C
C  POLAR ANGLE OF REFLECTION  (THETA BEI OLIVER)
C
C
      ZCPHI=COS(THETA_OUT)
C  LIMIT COSINE OF POLAR ANGLE TO 85. DEGREES
C  (I.E., 5 DEGREES AGAINST SURFACE TANGENTIAL PLANE)
      ZCPHI=MIN(0.999999_DP,MAX(0.08716_DP,ZCPHI))
      ZSPHI=SQRT(1.-ZCPHI*ZCPHI)
C
C  AZIMUTAL ANGLE OF REFLECTION    (ALPHA BEI OLIVER)
C
      ZCTHET=COS(ALPHA_OUT)
      ZCTHET=MAX(-.999999_DP,MIN(0.999999_DP,ZCTHET))
      ZSTHET=SQRT(1.-ZCTHET*ZCTHET)
      ZSTHET=ZSTHET*SIGN(1._DP,(RANF_EIRENE( )-0.5_DP))
C
      VX=-ZCPHI
      VY=ZSPHI*ZSTHET
      VZ=ZSPHI*ZCTHET
      IF (COSIN.GT.0.999999) THEN
        CALL ROTATF (VELX,VELY,VELZ,VX,VY,VZ,CRTX,CRTY,CRTZ)
      ELSE
        CALL ROTATE (VELX,VELY,VELZ,VX,VY,VZ,CRTX,CRTY,CRTZ,COSIN)
      ENDIF
      RETURN
C
C  ABSORB PARTICLE AT THIS SURFACE
C
700   CONTINUE
      IF (LSPUMP.AND.(MSURF.GT.0)) 
     .   SPUMP(ISPZO,MSURF)=SPUMP(ISPZO,MSURF)+WEIGHT
      LGPART=.FALSE.
      WEIGHT=0.
      ITYP=-1
      RETURN
C
C  ERROR MESSAGES FROM SUBR. REFLEC_PHOTON
C
C
992   CONTINUE
      WRITE (6,*) 'ERROR IN SUBR. REFLEC_PHOTON '
      MSS=MSURF
      IF (MSS.GT.NLIM) MSS=-(MSURF-NLIM)
      WRITE (6,*) 'MSURF = ',MSS
      WRITE (6,*) 'IGASF, IGAST ?? '
      WRITE (6,*) 'STOP HISTORY NO. NPANU= ',NPANU
      GOTO 999
c
993   CONTINUE
      WRITE (6,*) 'ERROR IN SUBR. REFLEC_PHOTON '
      MSS=MSURF
      IF (MSS.GT.NLIM) MSS=-(MSURF-NLIM)
      WRITE (6,*) 'MSURF = ',MSS
      WRITE (6,*) 'COSIN.LT.0. ', COSIN
      WRITE (6,*) 'STOP HISTORY NO. NPANU= ',NPANU
      GOTO 999
C
999   IF (NLTRC)  CALL CHCTRC(X0,Y0,Z0,16,18)
      LGPART=.FALSE.
      WEIGHT=0.
      RETURN
C
      END
