C EIRENE07 COMPILATION
C ===== SOURCE: phot_reflec.f
      module phot_reflec
      use precision
      use parmmod
      use cinit
      use comprt, only: iunout
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
      integer :: i, ifile
! Fehler abfragen...
!pb      open(unit = 23, file = "graphite_ext.dat")
      DO IFILE=1, NDBNAMES
        IF (INDEX(DBHANDLE(IFILE),'gr_ext') /= 0) EXIT
      END DO

      IF (IFILE > NDBNAMES) THEN
        WRITE (IUNOUT,*) ' NO DATABASENAME FOR graphite.ext.dat DEFINED'
        WRITE (IUNOUT,*) ' CALCULATION ABANDONNED '
        CALL EXIT_OWN(1)
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
        CALL EXIT_OWN(1)
      END IF

      OPEN (UNIT=23,FILE=DBFNAME(IFILE))
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
     .            RPROB, WABS
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
      WRITE (iunout,*) 'ERROR IN SUBR. REFLEC_PHOTON '
      MSS=MSURF
      IF (MSS.GT.NLIM) MSS=-(MSURF-NLIM)
      WRITE (iunout,*) 'MSURF = ',MSS
      WRITE (iunout,*) 'IGASF, IGAST ?? '
      WRITE (iunout,*) 'STOP HISTORY NO. NPANU= ',NPANU
      GOTO 999
c
993   CONTINUE
      WRITE (iunout,*) 'ERROR IN SUBR. REFLEC_PHOTON '
      MSS=MSURF
      IF (MSS.GT.NLIM) MSS=-(MSURF-NLIM)
      WRITE (iunout,*) 'MSURF = ',MSS
      WRITE (iunout,*) 'COSIN.LT.0. ', COSIN
      WRITE (iunout,*) 'STOP HISTORY NO. NPANU= ',NPANU
      GOTO 999
C
999   IF (NLTRC)  CALL CHCTRC(X0,Y0,Z0,16,18)
      LGPART=.FALSE.
      WEIGHT=0.
      RETURN
C
      END
C ===== SOURCE: transform.f
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
