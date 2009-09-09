c
c
c ======================================================================
c for transfer to mod_filament.f 
c
      SUBROUTINE LoadFilamentData(fname,status)
      USE mod_filament
      USE mod_filament_params
      IMPLICIT none

      INTEGER  , INTENT(OUT) :: status
      CHARACTER, INTENT(IN)  :: fname*(*)

      INTEGER fp,i1,idum(5),ifilament,nfilament1
      LOGICAL problem
      REAL    version

      status = 0

      fp = 99
      OPEN(UNIT=fp,FILE=TRIM(fname),ACCESS='SEQUENTIAL',
     .     FORM='UNFORMATTED',STATUS='OLD',ERR=98)            
c...  Check file version number:
      READ(fp,ERR=98) version
      IF (version.NE.1.0)
     .  CALL ER('LoadFilamentData','Unsupported version',*99)
c...  Check parameters used  in definition of TYPE_FILAMENT:
      READ(fp,ERR=98) (idum(i1),i1=1,5)
      problem = .FALSE.
      IF (idum(1).NE.MAX_FIL_NVTX  ) problem = .TRUE.
      IF (idum(2).NE.MAX_FIL_NCELL ) problem = .TRUE.
      IF (idum(3).NE.MAX_FIL_VCELL ) problem = .TRUE.
      IF (idum(4).NE.MAX_FIL_PARAMS) problem = .TRUE.
      IF (idum(5).NE.MAX_IR        ) problem = .TRUE.
      IF (problem)
     .  CALL ER('LoadFilamentData','Problem with TYPE_FILAMENT '//
     .          'parameter size(s)',*99)
c...  Load data:
      IF (ALLOCATED(filament)) DEALLOCATE(filament)  ! *** HACK *** Move elsewhere...
      READ(fp,ERR=98) nfilament1
      ALLOCATE(filament(nfilament1))
      READ(fp,ERR=98) (filament(ifilament),ifilament=1,nfilament1)
      nfilament = nfilament1
      CLOSE (fp)
      
      RETURN
 98   CALL ER('LoadFilamentData','Problems reading data file',*99)
 99   status = -1
      WRITE(0,*) '  FILE NAME: "'//TRIM(fname)//'"'
      RETURN
      END 
c
c ======================================================================
c for transfer to mod_filament.f 
c
      SUBROUTINE SaveFilamentData(fname)
      USE mod_filament
      USE mod_filament_params
      IMPLICIT none

      CHARACTER*(*), INTENT(IN) :: fname

      INTEGER fp,ifilament

      fp = 99
      OPEN(UNIT=fp,FILE=TRIM(fname),ACCESS='SEQUENTIAL',
     .     FORM='UNFORMATTED',STATUS='REPLACE',ERR=98)            
      WRITE(fp,ERR=98) 1.0
      WRITE(fp,ERR=98) MAX_FIL_NVTX,
     .                 MAX_FIL_NCELL,
     .                 MAX_FIL_VCELL,
     .                 MAX_FIL_PARAMS,
     .                 MAX_IR         
      WRITE(fp,ERR=98) nfilament
      WRITE(fp,ERR=98) (filament(ifilament),ifilament=1,nfilament)
      CLOSE (fp)
      
      RETURN
 98   CALL ER('SaveFilamentData','Problems writing data file',*99)
 99   WRITE(0,*) '  FILENAME = "'//TRIM(fname)//'"'
      STOP
      END 
c
c ======================================================================
c
c subroutine: DefineFilaments
c
c
      SUBROUTINE DefineFilaments
      USE mod_filament
      USE mod_options
      IMPLICIT none

      INTEGER ifilament
      REAL*8  t_start(40)

      DATA t_start / 0.0D-06, 50.0D-06, 10.0D-06, 80.0D-06, 20.0D-06,
     .              10.0D-06, 30.0D-06, 70.0D-06, 20.0D-06, 10.0D-06,
     .              40.0D-06, 20.0D-06, 60.0D-06, 20.0D-06, 50.0D-06,
     .              40.0D-06, 10.0D-06, 40.0D-06, 30.0D-06, 10.0D-06,
     .              90.0D-06,120.0D-06,100.0D-06,180.0D-06,150.0D-06,
     .             100.0D-06,110.0D-06, 90.0D-06,170.0D-06,140.0D-06,
     .             110.0D-06,100.0D-06,150.0D-06,170.0D-06,110.0D-06,
     .             150.0D-06,120.0D-06,140.0D-06, 90.0D-06,160.0D-06/ 

      SELECTCASE (opt_fil%opt)

        CASE (1)

c          t_start = 0.0D0
          opt_fil%scale(1) = 0.05D0
          opt_fil%scale(2) = 0.02D0
          opt_fil%scale(3) = 0.01D0
          
          nfilament = 12 ! 20
          IF (ALLOCATED(filament)) DEALLOCATE(filament)  ! *** HACK *** Move elsewhere...
          ALLOCATE(filament(nfilament))

          DO ifilament = 1, nfilament
            filament(ifilament)%status = 0
            filament(ifilament)%t_duration = 100.0D-06
            filament(ifilament)%t_start = t_start(ifilament)  ! -1.0D0  ! Time when filament is/was initiated
            filament(ifilament)%t_end   = t_start(ifilament) +
     .                                    filament(ifilament)%t_duration
          
            filament(ifilament)%ne_opt = 2               ! 
            filament(ifilament)%ne_param(1) =   3.0E+00  ! Number of points
            filament(ifilament)%ne_param(2) =   0.0E-06  ! t1
            filament(ifilament)%ne_param(3) =  50.0E-06  ! t2
            filament(ifilament)%ne_param(4) = 100.0E-06  ! t3
            filament(ifilament)%ne_param(5) =   2.5E+18  ! ne1
            filament(ifilament)%ne_param(6) =   5.0E+18  ! ne2
            filament(ifilament)%ne_param(7) =   2.5E+18  ! ne3
          
            filament(ifilament)%te_opt = 2               ! 
            filament(ifilament)%te_param(1) =   3.0E+00  ! Number of points
            filament(ifilament)%te_param(2) =   0.0E-06  ! t1
            filament(ifilament)%te_param(3) =  50.0E-06  ! t2
            filament(ifilament)%te_param(4) = 100.0E-06  ! t3
            filament(ifilament)%te_param(5) =  30.0      ! te1
            filament(ifilament)%te_param(6) =  30.0      ! te2
            filament(ifilament)%te_param(7) =   7.5      ! te3
          
            filament(ifilament)%crs_opt = 2             ! Outer midplane cross-section option
            filament(ifilament)%crs_param(1) = 0.015D0  ! Radial extent, +/- (RAD) (m)    ...numbers here from Glenn...
            filament(ifilament)%crs_param(2) = 0.035D0  ! Toroidal extent +/- (TOR) (m)
            filament(ifilament)%crs_param(3) = 11.0D0   !1.0D0     ! NRAD 
            filament(ifilament)%crs_param(4) = 11.0D0   !1.0D0     ! NPHI
          
            filament(ifilament)%tor_opt = 1               ! Linear from t=p1 to p2, at v=p3
            filament(ifilament)%tor_param(1) =   0.0D-06  ! t1 (s)
            filament(ifilament)%tor_param(2) =  50.0D-06  ! t2 (s)
            filament(ifilament)%tor_param(3) =  5.50D+03  ! v  (m s-1)  
            filament(ifilament)%tor_position =  0.0D0 ! 90.0D0  
     .                                  + DBLE(ifilament-1) * 30.0D0 ! 18.0D0
          
            filament(ifilament)%rad_opt = 1               ! Linear from t=p1 to p2, at v=p3
            filament(ifilament)%rad_param(1) =  50.0D-06  ! t1 (s)
            filament(ifilament)%rad_param(2) = 100.0D-06  ! t2 (s)
            filament(ifilament)%rad_param(3) = 1.000D+03  ! v  (m s-1)  
            filament(ifilament)%rad_position = -0.01D0 ! 0.000001D0 ! -0.01D0    ! Radial location at t=t1 (rho)
          ENDDO

        CASE(2)

          opt_fil%scale(1) = 0.10D0
          opt_fil%scale(2) = 0.05D0
          opt_fil%scale(3) = 0.02D0

          nfilament = 1 ! 5 ! 15
          IF (ALLOCATED(filament)) DEALLOCATE(filament)  ! *** HACK *** Move elsewhere...
          ALLOCATE(filament(nfilament))
          
          DO ifilament = 1, nfilament
            filament(ifilament)%status = 1
            filament(ifilament)%t_duration =  2.0D0
            filament(ifilament)%t_start    = -1.0D0
            filament(ifilament)%t_end      =  1.0D0

            filament(ifilament)%ne_opt = 2               ! 
            filament(ifilament)%ne_param(1) =   3.0E+00  ! Number of points
            filament(ifilament)%ne_param(2) =   0.0E-06  ! t1
            filament(ifilament)%ne_param(3) =  50.0E-06  ! t2
            filament(ifilament)%ne_param(4) = 100.0E-06  ! t3
            filament(ifilament)%ne_param(5) =   1.0E+19  ! ne1
            filament(ifilament)%ne_param(6) =   1.0E+19  ! ne2
            filament(ifilament)%ne_param(7) =   1.0E+19  ! ne3
          
            filament(ifilament)%te_opt = 2               ! 
            filament(ifilament)%te_param(1) =   3.0E+00  ! Number of points
            filament(ifilament)%te_param(2) =   0.0E-06  ! t1
            filament(ifilament)%te_param(3) =  50.0E-06  ! t2
            filament(ifilament)%te_param(4) = 100.0E-06  ! t3
            filament(ifilament)%te_param(5) =  30.0      ! te1
            filament(ifilament)%te_param(6) =  30.0      ! te2
            filament(ifilament)%te_param(7) =  30.0      ! te3
          
            filament(ifilament)%crs_opt = 2             ! Outer midplane cross-section option
            filament(ifilament)%crs_param(1) = 0.010D0  ! Radial extent, +/- (RAD) (m)    ...numbers here from Glenn...
            filament(ifilament)%crs_param(2) = 0.020D0  ! Toroidal extent +/- (TOR) (m)
            filament(ifilament)%crs_param(3) = 11.0D0    ! NRAD 
            filament(ifilament)%crs_param(4) = 11.0D0    ! NPHI
          
            filament(ifilament)%tor_opt = 1               ! Linear from t=p1 to p2, at v=p3
            filament(ifilament)%tor_param(1) =   0.0D-06  ! t1 (s)
            filament(ifilament)%tor_param(2) =  50.0D-06  ! t2 (s)
            filament(ifilament)%tor_param(3) =  0.00D+03  ! v  (m s-1)  
            filament(ifilament)%tor_position = 45.0D0 
          
            filament(ifilament)%rad_opt = 1               ! Linear from t=p1 to p2, at v=p3
            filament(ifilament)%rad_param(1) =  50.0D-06  ! t1 (s)
            filament(ifilament)%rad_param(2) = 100.0D-06  ! t2 (s)
            filament(ifilament)%rad_param(3) = 0.000D+03  ! v  (m s-1)  
            filament(ifilament)%rad_position = 0.05001D0  ! Radial location at t=t1 (rho)
     .                                  + DBLE(ifilament-1) * 0.001D0
          ENDDO

        CASE (3)  ! Fun with the linear grid

          opt_fil%scale(1) = 0.010D0
          opt_fil%scale(2) = 0.005D0
          opt_fil%scale(3) = 0.005D0
          
          nfilament = 1
          IF (ALLOCATED(filament)) DEALLOCATE(filament)  ! *** HACK *** Move elsewhere...
          ALLOCATE(filament(nfilament))
          
          DO ifilament = 1, nfilament
            filament(ifilament)%status = 0
            filament(ifilament)%t_duration = 100.0D-06
            filament(ifilament)%t_start = -1.0D0
            filament(ifilament)%t_end   = -1.0D0
          
            filament(ifilament)%ne_opt = 2               ! 
            filament(ifilament)%ne_param(1) =   3.0E+00  ! Number of points
            filament(ifilament)%ne_param(2) =   0.0E-06  ! t1
            filament(ifilament)%ne_param(3) =  50.0E-06  ! t2
            filament(ifilament)%ne_param(4) = 100.0E-06  ! t3
            filament(ifilament)%ne_param(5) =   1.0E+20  ! ne1
            filament(ifilament)%ne_param(6) =   1.0E+20  ! ne2
            filament(ifilament)%ne_param(7) =   1.0E+20  ! ne3
          
            filament(ifilament)%te_opt = 2               ! 
            filament(ifilament)%te_param(1) =   3.0E+00  ! Number of points
            filament(ifilament)%te_param(2) =   0.0E-06  ! t1
            filament(ifilament)%te_param(3) =  50.0E-06  ! t2
            filament(ifilament)%te_param(4) = 100.0E-06  ! t3
            filament(ifilament)%te_param(5) =  20.0      ! te1
            filament(ifilament)%te_param(6) =  20.0      ! te2
            filament(ifilament)%te_param(7) =  20.0      ! te3
          
            filament(ifilament)%crs_opt = 2             ! Outer midplane cross-section option
            filament(ifilament)%crs_param(1) = 1.00D-03 ! Radial extent, +/- (RAD) (m)    ...numbers here from Glenn...
            filament(ifilament)%crs_param(2) = 1.00D-03 ! Toroidal extent +/- (TOR) (m)
            filament(ifilament)%crs_param(3) = 1.0D0    !1.0D0     ! NRAD 
            filament(ifilament)%crs_param(4) = 3.0D0    !1.0D0     ! NPHI
          
            filament(ifilament)%tor_opt = 1               ! Linear from t=p1 to p2, at v=p3
            filament(ifilament)%tor_param(1) =   0.0D-06  ! t1 (s)
            filament(ifilament)%tor_param(2) = 100.0D-06  ! t2 (s)
            filament(ifilament)%tor_param(3) =  0.00D+03  ! v  (m s-1)  
            filament(ifilament)%tor_position =  0.0D0
          
            filament(ifilament)%rad_opt = 1               ! Linear from t=p1 to p2, at v=p3
            filament(ifilament)%rad_param(1) =   0.0D-06  ! t1 (s)
            filament(ifilament)%rad_param(2) = 100.0D-06  ! t2 (s)
            filament(ifilament)%rad_param(3) = 0.000D+00  ! v  (m s-1)  
            filament(ifilament)%rad_position = 7.500D-03  ! 0.000001D0 ! -0.01D0    ! Radial location at t=t1 (rho)
            
          ENDDO

        CASE (4)  ! ELM

          opt_fil%scale(1) = 0.15D0 !0.05D0
          opt_fil%scale(2) = 0.05D0 !0.02D0
          opt_fil%scale(3) = 0.02D0 !0.01D0
c          t_start = 0.0D0

          nfilament = 10
          IF (ALLOCATED(filament)) DEALLOCATE(filament)  ! *** HACK *** Move elsewhere...
          ALLOCATE(filament(nfilament))
          
          DO ifilament = 1, nfilament
            filament(ifilament)%status = 0
            filament(ifilament)%t_duration = 100.0D-06
            filament(ifilament)%t_start = t_start(ifilament)  ! -1.0D0  ! Time when filament is/was initiated
            filament(ifilament)%t_end   = t_start(ifilament) +
     .                                    filament(ifilament)%t_duration
          
            filament(ifilament)%ne_opt = 2               ! 
            filament(ifilament)%ne_param(1) =   3.0E+00  ! Number of points
            filament(ifilament)%ne_param(2) =   0.0E-06  ! t1
            filament(ifilament)%ne_param(3) =  50.0E-06  ! t2
            filament(ifilament)%ne_param(4) = 100.0E-06  ! t3
            filament(ifilament)%ne_param(5) =   5.0E+18  ! ne1
            filament(ifilament)%ne_param(6) =   4.0E+19  ! ne2
            filament(ifilament)%ne_param(7) =   1.0E+19  ! ne3
          
            filament(ifilament)%te_opt = 2               ! 
            filament(ifilament)%te_param(1) =   3.0E+00  ! Number of points
            filament(ifilament)%te_param(2) =   0.0E-06  ! t1
            filament(ifilament)%te_param(3) =  50.0E-06  ! t2
            filament(ifilament)%te_param(4) = 100.0E-06  ! t3
            filament(ifilament)%te_param(5) = 100.0      ! te1
            filament(ifilament)%te_param(6) = 100.0      ! te2
            filament(ifilament)%te_param(7) =  10.0      ! te3
          
            filament(ifilament)%crs_opt = 2             ! Outer midplane cross-section option
            filament(ifilament)%crs_param(1) = 0.020D0  ! Radial extent, +/- (RAD) (m)    ...numbers here from Glenn...
            filament(ifilament)%crs_param(2) = 0.020D0  ! Toroidal extent +/- (TOR) (m)
            filament(ifilament)%crs_param(3) = 11.0D0  !11.0D0        ! NRAD 
            filament(ifilament)%crs_param(4) = 11.0D0  !11.0D0        ! NPHI
          
            filament(ifilament)%tor_opt = 1               ! Linear from t=p1 to p2, at v=p3
            filament(ifilament)%tor_param(1) =   0.0D-06  ! t1 (s)
            filament(ifilament)%tor_param(2) =  50.0D-06  ! t2 (s)
            filament(ifilament)%tor_param(3) = 10.00D+03  ! v  (m s-1)  
            filament(ifilament)%tor_position = 0.0D0 ! 90.0D0
     .                                  + DBLE(ifilament-1) * 
     .                                  360.0D0 / DBLE(nfilament)
          
            filament(ifilament)%rad_opt = 1               ! Linear from t=p1 to p2, at v=p3
            filament(ifilament)%rad_param(1) =  50.0D-06  ! t1 (s)
            filament(ifilament)%rad_param(2) = 100.0D-06  ! t2 (s)
            filament(ifilament)%rad_param(3) = 1.400D+03  ! v  (m s-1)  
            filament(ifilament)%rad_position = -0.01D0    ! Radial location at t=t1 (rho)
          ENDDO

        CASE (5)  ! Ohmic, upper divertor

          opt_fil%scale(1) = 0.15D0  !0.15D0
          opt_fil%scale(2) = 0.05D0  !0.05D0
          opt_fil%scale(3) = 0.01D0  !0.02D0
c          t_start = 0.0D0

          nfilament = 40   ! 20-50
          IF (ALLOCATED(filament)) DEALLOCATE(filament)  ! *** HACK *** Move elsewhere...
          ALLOCATE(filament(nfilament))
          
          DO ifilament = 1, nfilament
            filament(ifilament)%status = 0
            filament(ifilament)%t_duration = 100.0D-06
            filament(ifilament)%t_start = t_start(ifilament)  ! -1.0D0  ! Time when filament is/was initiated
            filament(ifilament)%t_end   = t_start(ifilament) +
     .                                    filament(ifilament)%t_duration
          
            filament(ifilament)%ne_opt = 2               ! 
            filament(ifilament)%ne_param(1) =   3.0E+00  ! Number of points
            filament(ifilament)%ne_param(2) =   0.0E-06  ! t1
            filament(ifilament)%ne_param(3) =  50.0E-06  ! t2
            filament(ifilament)%ne_param(4) = 100.0E-06  ! t3
            filament(ifilament)%ne_param(5) =   0.5E+18  ! ne1
            filament(ifilament)%ne_param(6) =   4.0E+18  ! ne2
            filament(ifilament)%ne_param(7) =   1.0E+18  ! ne3
          
            filament(ifilament)%te_opt = 2               ! 
            filament(ifilament)%te_param(1) =   3.0E+00  ! Number of points
            filament(ifilament)%te_param(2) =   0.0E-06  ! t1
            filament(ifilament)%te_param(3) =  50.0E-06  ! t2
            filament(ifilament)%te_param(4) = 100.0E-06  ! t3
            filament(ifilament)%te_param(5) =  30.0      ! te1
            filament(ifilament)%te_param(6) =  30.0      ! te2
            filament(ifilament)%te_param(7) =   5.0      ! te3
          
            filament(ifilament)%crs_opt = 2              ! Outer midplane cross-section option
            filament(ifilament)%crs_param(1) = 0.0075D0  ! Radial extent, +/- (RAD) (m)    ...numbers here from Glenn...
            filament(ifilament)%crs_param(2) = 0.0500D0  ! Toroidal extent +/- (TOR) (m)
            filament(ifilament)%crs_param(3) = 11.0D0    !11.0D0        ! NRAD 
            filament(ifilament)%crs_param(4) = 11.0D0    !11.0D0        ! NPHI
          
            filament(ifilament)%tor_opt = 1               ! Linear from t=p1 to p2, at v=p3
            filament(ifilament)%tor_param(1) =   0.0D-06  ! t1 (s)
            filament(ifilament)%tor_param(2) =  50.0D-06  ! t2 (s)
            filament(ifilament)%tor_param(3) =  4.00D+03  ! v  (m s-1)  
            filament(ifilament)%tor_position = 0.0D0 ! 90.0D0
     .                                  + DBLE(ifilament-1) * 
     .                                  360.0D0 / DBLE(nfilament)
          
            filament(ifilament)%rad_opt = 1               ! Linear from t=p1 to p2, at v=p3
            filament(ifilament)%rad_param(1) =  50.0D-06  ! t1 (s)
            filament(ifilament)%rad_param(2) = 100.0D-06  ! t2 (s)
            filament(ifilament)%rad_param(3) = 0.500D+03  ! v  (m s-1)  
            filament(ifilament)%rad_position = 0.00001D0 ! +0.02D0    ! Radial location at t=t1 (rho)
          ENDDO

        CASE DEFAULT
          STOP 'NO FILAMENT DEFINITIONS SELECTED'
      ENDSELECT          
          
      RETURN
99    STOP
      END      
c
c ======================================================================
c
c subroutine: SetupFilaments
c
c
      SUBROUTINE SetupFilaments(t_global)
      USE mod_filament_params
      USE mod_filament
      USE mod_options
      IMPLICIT none

      REAL*8, INTENT(IN) :: t_global

      REAL FindSeparatrixRadius

      INTEGER ifilament,i1,i2,nrad,nphi,n,maxnphi,ikm,id,maxrad,
     .        nrings,rings(MAX_IR),icell,chop
      LOGICAL initialise
      REAL    x,y,z
      REAL*8  rad,tor,frad,phi,fphi,rmid,xmin,frac,radpos,torpos,rad1,
     .        rval,pval,rsep,x1,z1,ang,t,t_local,velocity,radius,arc,
     .        angle,t_delta,displacement,arclen,
     .        len1,len2,v(3,10000)
            
! Need a sub cross section radius, i.e. half the distance between, or can it be larger than the 
! actual distance between the chords, in particular now since I'm just assigning a uniform plasma
! across the cross section... whilst scanning down each chord each tetrahedron needs to be assigned
! a distance from the midplane... hey, need to crop each chord at this distance so that the filament
! can grow...  

      rsep = DBLE(FindSeparatrixRadius(1))

      DO ifilament = 1, nfilament

        initialise = .FALSE.

        IF (filament(ifilament)%t_start.EQ.-1.0D0) THEN
          filament(ifilament)%t_start = t_global
          filament(ifilament)%t_end   = t_global + 
     .                                  filament(ifilament)%t_duration
          filament(ifilament)%t_last  = 0.0D0
          initialise = .TRUE.
          filament(ifilament)%status = 1
        ELSEIF (t_global+1.D-06.LT.filament(ifilament)%t_start.OR.
     .          t_global-1.D-06.GT.filament(ifilament)%t_end  ) THEN
          WRITE(0,'(A,I6,2X,I6,2X,2F10.4)') ' ==FIL OFF:',
     .      NINT(t_global/1.0D-06),ifilament,
     .      SNGL(filament(ifilament)%t_start)/1.0E-06,
     .      SNGL(filament(ifilament)%t_end/1.0E-06)
          filament(ifilament)%status = 0
          CYCLE
        ELSE
          IF (filament(ifilament)%status.EQ.0) THEN
            initialise = .TRUE.
            filament(ifilament)%status = 1
          ENDIF
        ENDIF


c   **** NOTE *** Time keeping is not perfect here, and there's also the problem
c   that partial delta_t's are not registered properly when moving the filament,
c   i.e. if the current moment in time is just after or part way through the previous interval then
c   the filament won't be moved, even though it should have been...
        t       = t_global - filament(ifilament)%t_start
        t_delta = t - filament(ifilament)%t_last


c            WRITE(0,*) '==TIME  :',t,t_global,
c     .        filament(ifilament)%tor_param(1)


c...    Set the size of the filament cross-section:
        rad = filament(ifilament)%crs_param(1)   ! *** HACK *** Leave this fixed for now... 
        tor = filament(ifilament)%crs_param(2)


c...    Calculate centre-of-mass location of the filament:
        radpos = filament(ifilament)%rad_position
        torpos = filament(ifilament)%tor_position


c        WRITE(0,'(A,I6,2X,I6,2X,I6,2F10.4)') ' ==FIL POS:',
c     .    NINT(t_global/1.0D-06),ifilament,
c     .    NINT(t/1.0D-06),SNGL(torpos),SNGL(radpos)


        angle = 0.0D0
        SELECTCASE (filament(ifilament)%tor_opt)
          CASE(0)
          CASE(1)
            IF ((t.LT.filament(ifilament)%tor_param(1)).OR.
     .          (t.GT.filament(ifilament)%tor_param(2))) THEN
              filament(ifilament)%tor_status = 0
            ELSE
              filament(ifilament)%tor_status = 1
              velocity = filament(ifilament)%tor_param(3)
              arc = velocity * t_delta 
              radius = radpos + rsep
              angle = arc / radius * 180.0D0 / DPI             
            ENDIF
          CASE DEFAULT
            CALL ER('SetupFilaments','Unknown TOR_OPT',*99)
        ENDSELECT
        torpos = torpos + angle
c        WRITE(0,*) '===============================TORPOS:',
c     .    SNGL(t_delta),SNGL(angle),torpos

c...    Update radial position:
        displacement = 0.0D0
        SELECTCASE (filament(ifilament)%rad_opt)
          CASE(0)
          CASE(1)
            IF ((t.LT.filament(ifilament)%rad_param(1)).OR.
     .          (t.GT.filament(ifilament)%rad_param(2))) THEN
              filament(ifilament)%rad_status = 0
            ELSE
              filament(ifilament)%rad_status = 1
              velocity = filament(ifilament)%rad_param(3)
              displacement = velocity * t_delta 
            ENDIF
          CASE DEFAULT
            CALL ER('SetupFilaments','Unknown TOR_OPT',*99)
        ENDSELECT
        radpos = radpos + displacement
c        WRITE(0,*) '==RADPOS:',
c     .    SNGL(t_delta),SNGL(displacement),SNGL(radpos)


        filament(ifilament)%rad_position = radpos
        filament(ifilament)%tor_position = torpos

c...    Plasma:

        WRITE(0,'(A,I6,2X,I6,2X,I6,2F10.4)') ' ==FIL POS:',
     .    NINT(t_global/1.0D-06),ifilament,
     .    NINT(t/1.0D-06),SNGL(torpos),SNGL(radpos)



c...    Setup cross-section vertex array:
        SELECTCASE(filament(ifilament)%crs_opt)   
          CASE(0)
c...        Debug:
            n = 4
            filament(ifilament)%vtx(1:3,1:11) = 0.0D0
c            filament(ifilament)%vtx(1,1)  = -0.005D0 !+ 0.02D0
c            filament(ifilament)%vtx(1,2)  = -0.004D0 !+ 0.02D0
c            filament(ifilament)%vtx(1,3)  = -0.003D0 !+ 0.02D0
c            filament(ifilament)%vtx(1,4)  = -0.002D0 !+ 0.02D0
            filament(ifilament)%vtx(1,1)  = -0.0001D0 + radpos
            filament(ifilament)%vtx(1,2)  =  0.0000D0 + radpos
            filament(ifilament)%vtx(1,3)  =  0.0001D0 + radpos
            filament(ifilament)%vtx(1,4)  =  0.0002D0 + radpos
c            filament(ifilament)%vtx(1,9)  =  0.003D0 !+ 0.02D0
c            filament(ifilament)%vtx(1,10) =  0.004D0 !+ 0.02D0
c            filament(ifilament)%vtx(1,11) =  0.005D0 !+ 0.02D0
          CASE(1)
c...        Concentric circles, but this isn't great because it gives even weighting to the radial
c           and toroidal directions, and the radial direction requires much higher spatial resolution
c           near the sepearatrix due to the strong shear:
c           *** NOT IN USE *** (may not work)
            STOP 'POSSIBLY DEFUNCT'
            nrad = 2 ! 4 ! 5
            maxnphi = 50 ! 10
            n = 0
            DO i1 = 1, nrad
c              IF (i1.EQ.1) nphi = 1
c              IF (i1.EQ.2) nphi = 3
c              IF (i1.GT.2) nphi = nphi * 2
              frad = DBLE(i1 - 1) / DBLE(nrad - 1)
              nphi = NINT(frad * (maxnphi - 1)) + 1
              DO i2 = 1, nphi
                fphi = DBLE(i2 - 1) / DBLE(nphi)
                phi = 2.0D0 * 3.14159265D0 * fphi
                n = n + 1
                filament(ifilament)%vtx(1,n) = frad * rad * DSIN(phi) + 
     .                                         radpos
                filament(ifilament)%vtx(2,n) = 0.0D0
                filament(ifilament)%vtx(3,n) = frad * tor * DCOS(phi)
c                WRITE(0,*) 'PAIN:',i1,i2,frad,fphi
c                WRITE(0,*) 'PAIN:',i1,i2,
c     .            SNGL(filament(ifilament)%vtx(1:3,n))
              ENDDO
            ENDDO

          CASE(2)
c...        Elliptical, but vertices are now set along radial segments, with
c           even spacing and much lower frequency in the toroidal direction:
            nrad = NINT(filament(ifilament)%crs_param(3))
            nphi = NINT(filament(ifilament)%crs_param(4))
c            WRITE(0,*) 'NRAD,NPHI,RSEP:',nrad,nphi,SNGL(rsep)
            n = 0
            DO i1 = 1, nphi
              IF (nphi.EQ.1) THEN
                fphi = 0.5D0
              ELSE
                fphi = DBLE(i1 - 1) / DBLE(nphi - 1)
              ENDIF
              frad = DCOS(DPI * (2.0D0 * fphi - 1.0D0) / 2.0D0)
c              frad = COS(3.14159265D0 * (2.0D0 * fphi - 1.0D0) / 2.0D0)
c              WRITE(0,*) 'FRAD:',frad,MAX(1,NINT(nrad*frad)),fphi
              maxrad = MAX(1,NINT(nrad*frad))
              DO i2 = 1, maxrad
                IF (maxrad.EQ.1) THEN
                  frac = 0.0D0
                ELSE
                  frac = 2.0D0 * DBLE(i2 - 1) / DBLE(maxrad - 1) - 1.0D0
                ENDIF
                n = n + 1

                radius = frac * frad * rad + radpos + rsep
                arclen = 2.0D0 * (0.5D0 - fphi) * tor

                IF (.TRUE.) THEN
                  x1 = radius
                  z1 = 0.0D0
                  ang = arclen / radius
                  filament(ifilament)%vtx(1,n) = DCOS(ang) * x1 -
     .                                           DSIN(ang) * z1
                  filament(ifilament)%vtx(3,n) = DSIN(ang) * x1 +
     .                                           DCOS(ang) * z1
c                  WRITE(0,*) 'RAD,ARC:',radius,arclen
                ELSE
                  filament(ifilament)%vtx(1,n) = radius
                  filament(ifilament)%vtx(2,n) = 0.0D0
                  filament(ifilament)%vtx(3,n) = arclen
                ENDIF

c                WRITE(0,*) 'FIL:',n,
c     .                      REAL(filament(ifilament)%vtx(1:1,n))-rsep,
c     .                      REAL(filament(ifilament)%vtx(2:3,n))
c...            Displace the filament toroidally:
                ang = torpos * DPI / 180.0D0 ! -1.0D0 * torpos
c        WRITE(0,*) '===============================ANG:',
c     .    ang,ang*180.0D0/DPI
                rval = DSQRT(filament(ifilament)%vtx(1,n)**2 + 
     .                       filament(ifilament)%vtx(3,n)**2)
c                WRITE(0,*) '   :',n,
c     .                      REAL(filament(ifilament)%vtx(1:3,n)),
c     .                      REAL(rval)-rsep
                x1 = filament(ifilament)%vtx(1,n)
                z1 = filament(ifilament)%vtx(3,n)
                filament(ifilament)%vtx(1,n) = DCOS(ang) * x1 -
     .                                         DSIN(ang) * z1
                filament(ifilament)%vtx(3,n) = DSIN(ang) * x1 +
     .                                         DCOS(ang) * z1
c Debug:
                rval = DSQRT(filament(ifilament)%vtx(1,n)**2 + 
     .                       filament(ifilament)%vtx(3,n)**2)
c                WRITE(0,*) '   :',n,
c     .                      REAL(filament(ifilament)%vtx(1:3,n)),
c     .                      REAL(rval)-rsep
              ENDDO
            ENDDO

          CASE DEFAULT
            CALL ER('SetupFilament','Unrecognised mode',*99)
        ENDSELECT

        filament(ifilament)%nvtx  = n
        filament(ifilament)%ncell = n  ! *** HACK *** This is convenient, for now...

c...    Project filament cross-section down onto the midplane:

c...    Set the length of an individual field line associated with 
c       a filament flux-tube:
        IF (initialise) THEN
          DO icell = 1, filament(ifilament)%ncell  ! *** HACK *** As above...
            x = SNGL(filament(ifilament)%vtx(1,icell))
            y = 0.0
            z = SNGL(filament(ifilament)%vtx(3,icell))
            chop = 1
            IF     (opt_fil%length1.EQ.-1.0.OR.
     .              opt_fil%length2.EQ.-1.0) THEN
              chop = 3
            ELSEIF (opt_fil%length1.GE.0.0.OR.
     .              opt_fil%length2.GE.0.0) THEN
              len1 = MAX(0.0D0,DBLE(opt_fil%length1))
              len2 = MAX(0.0D0,DBLE(opt_fil%length2))
              chop = 5
            ENDIF
            CALL TraceFieldLine_DIVIMP(x,y,z,2,chop,len1,len2,
     .                                 n,v,10000)  ! *** HACK *** (the 10000)
            filament(ifilament)%lcell(1,icell) = len1
            filament(ifilament)%lcell(2,icell) = len2
c            WRITE(0,*) 'LEN1,2:',len1,len2
          ENDDO
        ELSE
          IF (.TRUE.) THEN
            DO icell = 1, filament(ifilament)%ncell  ! *** HACK *** As above...
              len1 = filament(ifilament)%lcell(1,icell)
              len2 = filament(ifilament)%lcell(2,icell)          
              len1 = len1 + 0.10D0
              len2 = len2 + 0.10D0
              filament(ifilament)%lcell(1,icell) = len1
              filament(ifilament)%lcell(2,icell) = len2          
            ENDDO
          ENDIF
        ENDIF

c...    Identify which rings/tubes on the fluid grid are within striking
c       distance of the filament:
        DO i1 = 1, filament(ifilament)%nvtx
          rings = 0
          rad1 = DSQRT(filament(ifilament)%vtx(1,i1)**2 +  
     .                 filament(ifilament)%vtx(3,i1)**2) - rsep
          CALL SelectGridRegion_DIVIMP(SNGL(rad1),nrings,rings,MAX_IR)
          filament(ifilament)%ir_space(0:MAX_IR,i1) = 0
          filament(ifilament)%ir_space(0       ,i1) = nrings
          filament(ifilament)%ir_space(1:nrings,i1) = rings(1:nrings)
c          WRITE(0,*) 'RINGS:',i1,
c     .               filament(ifilament)%ir_space(0:nrings,i1)
        ENDDO

c...    Update time stamp:
        filament(ifilament)%t_last  = t

      ENDDO  ! Main filament loop
      
      RETURN
99    STOP
      END
c
c ======================================================================
c ======================================================================
c
c function: CalcPerp
c
c Calculate the perpendicular distance from a point to a line.
c
      REAL*8 FUNCTION CalcPerp(a,b,c,t) 

      IMPLICIT none

      REAL*8, INTENT(IN)  :: a(3),b(3),c(3)
      REAL*8, INTENT(OUT) :: t

      REAL*8 p(3),delta(3),dist

      REAL*8     DTOL
      PARAMETER (DTOL = 1.0D-7)


      CalcPerp = -1.0D0

      delta = b - a

      IF (DABS(delta(1)).GT.DTOL.OR.DABS(delta(2)).GT.DTOL.OR.
     .    DABS(delta(3)).GT.DTOL) THEN

        t = ((c(1) - a(1)) * delta(1) + (c(2) - a(2)) * delta(2) + 
     .       (c(3) - a(3)) * delta(3)) /
     .      (delta(1)**2 + delta(2)**2 + delta(3)**2)

        p(1:3) = a(1:3) + t * delta(1:3)

        IF ((t+DTOL).GE.0.0D0.AND.(t-DTOL).LE.1.0D0) THEN

          dist = DSQRT((p(1) - c(1))**2 + (p(2) - c(2))**2 +
     .                 (p(3) - c(3))**2)

c          WRITE(0,*) 'A=',a
c          WRITE(0,*) 'B=',b
c          WRITE(0,*) 'C=',c
c          WRITE(0,*) 'P=',p
c          WRITE(0,*) 'DELTA=',delta

          IF (dist.LT.DTOL*10.0D0) THEN
c           Point C is on the line AB:
            CalcPerp = 0.0D0
          ELSE
c           Point of perpendicular intersection is displaced from
c           the point C:
            CalcPerp = dist
          ENDIF
        ELSE
          CalcPerp = -1.0D0
        ENDIF
      ELSE
c       If the points are all identicle, then return a positive result,
c       otherwise indicate that the problem was ill-posed:
        IF (DABS(a(1) - c(1)).LT.DTOL.AND.DABS(a(2) - c(2)).LT.DTOL.AND.
     .      DABS(a(3) - c(3)).LT.DTOL) THEN
          CalcPerp =  0.0D0
        ELSE
          CalcPerp = -1.0D0
        ENDIF
      ENDIF

      RETURN
      END
c
c ======================================================================
c
c subroutine: gmCalcCrossProduct
c
      SUBROUTINE gmCalcCrossProduct(v1,v2,xproduct)
      IMPLICIT none

      REAL*8, INTENT(IN)  :: v1(3),v2(3)
      REAL*8, INTENT(OUT) :: xproduct(3)

      xproduct(1) = v1(2) * v2(3) - v1(3) * v2(2)
      xproduct(2) = v1(3) * v2(1) - v1(1) * v2(3)
      xproduct(3) = v1(1) * v2(2) - v1(2) * v2(1)

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: gmCalcDotProduct
c
      SUBROUTINE gmCalcDotProduct(v1,v2,dproduct)
      IMPLICIT none

      REAL*8, INTENT(IN)  :: v1(3),v2(3)
      REAL*8, INTENT(OUT) :: dproduct

      dproduct = v1(1) * v2(1) + v1(2) * v2(2) + v1(3) * v2(3)

      RETURN
99    STOP
      END
c
c ======================================================================
c
c function: CalcTriangleArea
c
c see http://geometryalgorithms.com/Archive/algorithm_0101/algorithm_0101.htm
c
      REAL*8 FUNCTION gmCalcTriangleArea(isrf)
      USE mod_geometry
      IMPLICIT none

      INTEGER, INTENT(IN) :: isrf 
      REAL*8 area,v1(3),v2(3),xprod(3)

      IF (srf(isrf)%nvtx.NE.3) 
     .  CALL ER('gmCalcTriangleArea','Surface not a triangle',*99)

c...  Calculate cross-product about vertex #2:
      v1(1:3) = vtx(1:3,srf(isrf)%ivtx(1)) -
     .          vtx(1:3,srf(isrf)%ivtx(2))
      v2(1:3) = vtx(1:3,srf(isrf)%ivtx(3)) -
     .          vtx(1:3,srf(isrf)%ivtx(2))

      CALL gmCalcCrossProduct(v1,v2,xprod)

c...  Quick check:
      IF (xprod(1).EQ.0.0D0.AND.xprod(2).EQ.0.0D0.AND.xprod(3).EQ.0.0D0) 
     .  CALL ER('gmCalcTriangleArea','Zero area detected',*99)

c...  Length of 0.5*|V1xV2| gives the area of a 3D triangle:
      area = 0.5D0 * DSQRT(xprod(1)**2 + xprod(2)**2 + xprod(3)**2)

      gmCalcTriangleArea = area

      RETURN
99    WRITE(0,*) '  ISRF,NVTX = ',isrf,srf(isrf)%nvtx
      STOP
      END
c
c ======================================================================
c
c function: gmCalcTetrahedronVolume
c
c see http://mathworld.wolfram.com/Tetrahedron.html
c
      REAL*8 FUNCTION gmCalcTetrahedronVolume(iobj)
      USE mod_geometry
      IMPLICIT none

      INTEGER, INTENT(IN) :: iobj

      INTEGER iside,isrf,ivtx,ipts(4)
      REAL*8  p(3),v1(3),v2(3),v3(3),xprod(3),dprod,volume

c...  Confirm that specified object is a tetrahedron:
      IF (obj(iobj)%nside.NE.4) 
     .  CALL ER('gmCalcTetrahedronVolume','Object not a '//
     .          'tetrahedron',*99)
      DO iside = 1, obj(iobj)%nside
        isrf = ABS(obj(iobj)%iside(iside))
        IF (srf(isrf)%nvtx.NE.3) 
     .    CALL ER('gmCalcTetrahedronVolume','Object not a '//
     .            'tetrahedron',*99)
      ENDDO

c...  Get vectors required for calculating the volume:
      isrf = ABS(obj(iobj)%iside(1))
      ipts(1:3) = srf(isrf)%ivtx(1:3)
      v1(1:3) = vtx(1:3,ipts(1)) - vtx(1:3,ipts(2))
      v2(1:3) = vtx(1:3,ipts(3)) - vtx(1:3,ipts(2))

c...  Forth point, or apex point of sorts, from any of the other sides:
      isrf = ABS(obj(iobj)%iside(2))
      DO ivtx = 1, 3
        IF (srf(isrf)%ivtx(ivtx).NE.ipts(1).AND.    ! *** REALLY NEED A FUNCTION FOR THIS ***
     .      srf(isrf)%ivtx(ivtx).NE.ipts(2).AND.
     .      srf(isrf)%ivtx(ivtx).NE.ipts(3)) EXIT
      ENDDO
      IF (ivtx.EQ.4) 
     .  CALL ER('gmCalcTetrahedronVolume','4th tetrahedron vertex '//
     .          'not found',*99)
      ipts(4) = srf(isrf)%ivtx(ivtx)
      v3(1:3) = vtx(1:3,ipts(4)) - vtx(1:3,ipts(2))

c...  Volume of a tetrahedron is given by (1/6)*|V1.(V2xV3)|:
      CALL gmCalcCrossProduct(v2,v3,xprod)
      CALL gmCalcDotProduct(v1,xprod,dprod)
      volume = (1.0D0 / 6.0D0) * DABS(dprod)

c...  Quick check:
      IF (volume.LT.1.0D-10) 
     .  CALL ER('gmCalcTetrahedronVolume','Zero volume detected',*99)

      gmCalcTetrahedronVolume = volume

      RETURN
99    WRITE(0,*) '  IOBJ,ISRF,NVTX = ',iobj,isrf,srf(isrf)%nvtx
      STOP
      END
c
c ======================================================================
c
c subroutine: gmCalcSurfaceArea
c
c
      REAL*8 FUNCTION gmCalcSurfaceArea(isrf)
      USE mod_geometry
      IMPLICIT none

      REAL*8 gmCalcTriangleArea  ! This will go in mod_geometry...
 
      INTEGER, INTENT(IN) :: isrf
      INTEGER ndim
      REAL*8  area

      ndim = srf(isrf)%nvtx

      SELECTCASE (ndim)
c        CASE (1)  ! Point, error!:
c        CASE (2)  ! Line segment, assume a toroidally symmetric surface:
        CASE (3)  ! Triangle:
          area = gmCalcTriangleArea(isrf)
        CASE DEFAULT
          CALL ER('gmCalcSurfaceArea','Invalid number of vertices',*99)
      ENDSELECT

      gmCalcSurfaceArea = area

      RETURN
 99   WRITE(0,*) '  ISRF,NVTX = ',isrf,ndim
      STOP
      END

c
c ======================================================================
c
c subroutine: CalcCentroid
c
c_
      SUBROUTINE CalcCentroid(index,mode,p)
      USE mod_eirene06_locals
      IMPLICIT none

      INTEGER, INTENT(IN)  :: index,mode
      REAL*8 , INTENT(OUT) :: p(3)
      
      INTEGER i1,iobj,iside,isrf,ivtx
      REAL*8  count

      SELECTCASE (mode)
        CASE (1)  ! Surface:
          IF (srf(index)%nvtx.EQ.0) 
     .      CALL ER('CalcCentroid','NVTX.EQ.0',*99)
          p = 0.0D0
          DO i1 = 1, srf(index)%nvtx
            ivtx = srf(index)%ivtx(i1)
            p(1:3) = p(1:3) + vtx(1:3,ivtx)
          ENDDO 
          p(1:3) = p(1:3) / DBLE(srf(index)%nvtx)
        CASE (2)  ! Object:
          iobj = index
          p = 0.0D0
          count = 0.0D0
          DO iside = 1, obj(iobj)%nside
            isrf = ABS(obj(iobj)%iside(iside))
            DO ivtx = 1, srf(isrf)%nvtx
              count = count + 1.0D0
              p(1:3) = p(1:3) + vtx(1:3,srf(isrf)%ivtx(ivtx))
            ENDDO
          ENDDO
          p(1:3) = p(1:3) / count
        CASE DEFAULT
          CALL ER('CalcCentroid','Unrecognised MODE',*99)
      ENDSELECT

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c
c
c 
      SUBROUTINE FindConnectedTetrahedron(found,omap,smap,iside2,
     .                                    iobj,iside,ivtx1,ivtx2)
      USE mod_geometry
      IMPLICIT none

      LOGICAL, INTENT(OUT) :: found
      INTEGER, INTENT(OUT) :: omap,smap,iside2
      INTEGER, INTENT(IN)  :: iobj,iside,ivtx1,ivtx2

      INTEGER i1,i2,i3,ivtx3,ivtx4

      found = .FALSE.
      omap = 0
      smap = 0
      iside2 = 0

      DO i1 = 1, 4  ! Loop over sides
        IF (i1.EQ.iside) CYCLE
        DO i2 = 1, 3  ! Loop of side line segments
          i3 = i2 + 1
          IF (i3.EQ.4) i3 = 1
          ivtx3 = srf(ABS(obj(iobj)%iside(i1)))%ivtx(i2)
          ivtx4 = srf(ABS(obj(iobj)%iside(i1)))%ivtx(i3)
          IF ((ivtx1.EQ.ivtx3.AND.ivtx2.EQ.ivtx4).OR.
     .        (ivtx1.EQ.ivtx4.AND.ivtx2.EQ.ivtx3)) THEN
            omap = obj(iobj)%omap(i1)
            smap = obj(iobj)%smap(i1)
            iside2 = i1
            found = .TRUE.
            EXIT
          ENDIF
        ENDDO
        IF (found) EXIT 
      ENDDO

      RETURN
 99   STOP
      END

c
c ======================================================================
c
c subroutine: CheckTetrahedronStructure
c
      SUBROUTINE CheckTetrahedronStructure
      USE mod_geometry
      USE mod_filament_params
      IMPLICIT none

      INTEGER fp,iobj,iside,isrf,ivtx,nmatch,i1,i2,s(12),count
      LOGICAL match
      REAL*8  cen_obj(3),cen(3),v1(3),v2(3),xprod(3),len_cen,len_xprod,
     .        cos_theta,theta

      fp = 0

      WRITE(fp,*) '=== CHECKING TETRAHEDRON STRUCTURE ==='

c...  Make sure surface orientation is correct:
      DO iobj = 1, nobj 
        IF (grp(obj(iobj)%group)%type.NE.GRP_TETRAHEDRON) CYCLE

c...    Get centroid of the object:
        CALL CalcCentroid(iobj,2,cen_obj) 
c...    Scan over sides:
        DO iside = 1, obj(iobj)%nside        
          isrf = obj(iobj)%iside(iside)
          IF (isrf.GT.0) THEN
c...        Calculate cross-product about vertex 2:
            v1(1:3) = vtx(1:3,srf(isrf)%ivtx(1)) -
     .                vtx(1:3,srf(isrf)%ivtx(2))
            v2(1:3) = vtx(1:3,srf(isrf)%ivtx(3)) -
     .                vtx(1:3,srf(isrf)%ivtx(2))

            CALL gmCalcCrossProduct(v1,v2,xprod)
c            xprod(1) = v1(2) * v2(3) - v1(3) * v2(2)
c            xprod(2) = v1(3) * v2(1) - v1(1) * v2(3)
c            xprod(3) = v1(1) * v2(2) - v1(2) * v2(1)

c...        Vector from object center to vertex 2:
            cen(1:3) = cen_obj(1:3) - vtx(1:3,srf(isrf)%ivtx(2))

            len_xprod = DSQRT(xprod(1)**2 + xprod(2)**2 + xprod(3)**2)
            len_cen   = DSQRT(cen  (1)**2 + cen  (2)**2 + cen  (3)**2)

c...        Use dot product to get the angle between the cross product
c           and the centroid, which should always be less than 90 degrees
c           if the surface is anticlockwise when viewed from outside
c           the tetrahedron:
            cos_theta = (xprod(1) * cen(1) + xprod(2) * cen(2) + 
     .                   xprod(3) * cen(3)) / len_xprod / len_cen
            theta = DACOS(cos_theta) * 180.0D0 / DPI

            IF (theta.GT.90.0D0) THEN
              WRITE(fp,*) '==SURFACE ORIENTATION VIOLATION============'
              WRITE(fp,*) '   IOBJ,ISRF,THETA:',iobj,isrf,theta
            ENDIF
          ENDIF
        ENDDO 
      ENDDO

c...  Check that there are only 4 vertices assigned to a given tetrahedron:
      DO iobj = 1, nobj
        IF (grp(obj(iobj)%group)%type.NE.GRP_TETRAHEDRON) CYCLE

        isrf = ABS(obj(iobj)%iside(1))
        s(1:3) = srf(isrf)%ivtx(1:3)
        nmatch = 3
        DO iside = 2, obj(iobj)%nside
          isrf = ABS(obj(iobj)%iside(iside))
          DO ivtx = 1, srf(isrf)%nvtx
            match = .FALSE.
            DO i1 = 1, nmatch
              IF (s(i1).EQ.srf(isrf)%ivtx(ivtx)) match = .TRUE.
            ENDDO
            IF (.NOT.match) THEN
              nmatch = nmatch + 1
              s(nmatch) = srf(isrf)%ivtx(ivtx)
              IF (nmatch.GT.4) THEN
                WRITE(fp,*) '==MORE THAN 4 VERTICES PER TETRAHEDRON==='
                WRITE(fp,*) '   IOBJ,ISRF,NMATCH:',iobj,isrf,nmatch
                WRITE(fp,*) '   S               :',s(1:nmatch)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO


c...  Check that the EIRENE assumptions are satisfied, i.e. surface #2 corresponds to 
c     side #1 on surface #1, etc.  Note however that the vertices are not 
c     guaranteed to match up exactly:
      DO iobj = 1, nobj
        IF (grp(obj(iobj)%group)%type.NE.GRP_TETRAHEDRON) CYCLE

        isrf = obj(iobj)%iside(1)
        IF (isrf.GT.0) THEN
          s(1:3) = srf(isrf)%ivtx(1:3)
        ELSE
          s(1) = srf(-isrf)%ivtx(3)  ! This must match the convention in WriteEireneObjects
          s(2) = srf(-isrf)%ivtx(2)
          s(3) = srf(-isrf)%ivtx(1)
        ENDIF

        DO iside = 2, obj(iobj)%nside
          isrf = ABS(obj(iobj)%iside(iside))
          count = 0
          DO ivtx = 1, srf(isrf)%nvtx
            i1 = iside - 1
            i2 = i1 + 1
            IF (i2.EQ.4) i2 = 1
            IF (srf(isrf)%ivtx(ivtx).EQ.s(i1).OR.
     .          srf(isrf)%ivtx(ivtx).EQ.s(i2)) 
     .        count = count + 1
          ENDDO
          IF (count.NE.2) THEN
            WRITE(0,*) 'PROBLEM-O!',iobj,iside,isrf
            obj(iobj)%segment(1) = 3
          ENDIF
        ENDDO
      ENDDO



c...  Make sure that only surface #1 is assocated with boundaries:




      WRITE(fp,*) '=== DONE ==='

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: DivideTetrahedron
c
c Some subtleties here if all of the tetrahedrons resulting from sub-
c division of the original tetrahedron are to have sides witht the
c correct counter-clockwise orienation when viewed from the outside.
c
c Only surface 1 of any tetrahedron should be associated with a wall
c surface.
c
c
c
c
c
c
      SUBROUTINE DivideTetrahedron(iobj)
      USE mod_eirene06_locals
      IMPLICIT none

      INTEGER, INTENT(IN) :: iobj

      TYPE(type_srf   ) newsrf
      TYPE(type_object) newobj

      INTEGER s(3,12)
      DATA    s / 1, 5, 2,    
     .            2, 5, 3,    
     .            3, 5, 1,    

     .            2, 5, 1, 
     .            1, 5, 4,  
     .            4, 5, 2,

     .            3, 5, 2, 
     .            2, 5, 4,  
     .            4, 5, 3,

     .            1, 5, 3, 
     .            3, 5, 4,  
     .            4, 5, 1 / 

c      INTEGER s(3,6)
c      DATA    s / 1, 5, 2,    
c     .            2, 5, 3,    
c     .            3, 5, 1,    
c     .            1, 5, 4, 
c     .            2, 5, 4,  
c     .            3, 5, 4 / 


c      INTEGER t(4,4)
c      DATA    t / 1,  5,   6,  7,
c     .            2,  8,  -9, -5,
c     .            3,  9, -10, -6,   
c     .            4, 10,  -8, -7 / 

      INTEGER u(3,12)
      DATA    u / 1, 2, 5,
     .            1, 4, 2,
     .            2, 4, 5,   
     .            5, 4, 1,  

     .            2, 3, 5,
     .            2, 4, 3,
     .            3, 4, 5,   
     .            5, 4, 2,  

     .            3, 1, 5,
     .            3, 4, 1,
     .            1, 4, 5,   
     .            5, 4, 3 /       

      INTEGER v(3,8)
      DATA    v / 1, 5, 3,
     .            1, 4, 5,
     .            5, 4, 3,   
     .            0, 0, 0,
    
     .            5, 2, 3,
     .            5, 4, 2, 
     .            0, 0, 0,
     .            3, 4, 5 /

      INTEGER fp,i1,i2,i3,i4,i5,isrf,ivtx,istart_srf,iside,omap,smap,
     .        ivtx1,ivtx2,ivtx3,ivtx4,iobj1,iside1,smap1,omap1,count,
     .        iside2,vertex_centroid,nmatch,idum1,obj_segment
      INTEGER vertex(5),surface(10),object(5),change,obj_iside(4)
      LOGICAL subdivide,found,finished,subdivision,output
      REAL*8  p(3)

   
      output = .FALSE.
      fp = 0


      IF (obj(iobj)%segment(1).LE.0)                   ! Need to refine this check probably...
     .  CALL ER('DivideTetrahedron','Invaid call',*99)
 

      subdivision = .FALSE.

      obj_segment = obj(iobj)%segment(1)

      SELECTCASE(0) 
        CASE (0)
c...      Divide into 4 polygons:

          object(1) = iobj
          obj(object(1))%segment(1) = -1  ! Tagged for deletion

c...      Collect list of vertices:
          isrf = obj(iobj)%iside(1)   ! Points 1 to 3
          IF (isrf.GT.0) THEN
            vertex(1:3) = srf(isrf)%ivtx(1:3)
          ELSE
            vertex(1) = srf(-isrf)%ivtx(3)
            vertex(2) = srf(-isrf)%ivtx(2)
            vertex(3) = srf(-isrf)%ivtx(1)
          ENDIF
c          isrf        = ABS(obj(iobj)%iside(2))   ! Point 4
c          vertex(4)   = srf(isrf)%ivtx(2)
          isrf = ABS(obj(iobj)%iside(2))                  ! Can't assume a particular orientation for any 
          DO ivtx = 1, srf(isrf)%nvtx                     ! surface with respect to the other surfaces in
            IF (srf(isrf)%ivtx(ivtx).NE.vertex(1).AND.    ! the tetrahedron so always need to search for
     .          srf(isrf)%ivtx(ivtx).NE.vertex(2).AND.    ! a point that doesn't belong on the base, which
     .          srf(isrf)%ivtx(ivtx).NE.vertex(3)) EXIT   ! should always be surface 1...
          ENDDO
          vertex(4) = srf(isrf)%ivtx(ivtx)
          p(1) = DBLE(obj(iobj)%x)    ! Point 5
          p(2) = DBLE(obj(iobj)%y)
          p(3) = DBLE(obj(iobj)%z)
          vertex(5) = AddVertex(p)
          IF (output) WRITE(fp,*) 'VERTEX:',vertex(1:5)

c...      Dividing the tetrahedron about the centroid:
          istart_srf = nsrf

c...      Divide tetrahedron in 4, leaving the surfaces of the original
c         tetrahedron in tact:
          DO iside = 1, 4
            IF (iside.EQ.1) THEN
              i1 = 1
            ELSE    
              i2 = iside - 1
              i3 = i2 + 1
              IF (i3.GT.3) i3 = 1
              i1 = 0
              DO i4 = 2, 4
                nmatch = 0
                isrf = ABS(obj(iobj)%iside(i4))
                DO ivtx = 1, srf(isrf)%nvtx
                  IF (vertex(i2).EQ.srf(isrf)%ivtx(ivtx).OR.
     .                vertex(i3).EQ.srf(isrf)%ivtx(ivtx))
     .               nmatch = nmatch + 1
c                  IF (output) WRITE(fp,*) '  --->',iside,isrf,
c     .              vertex(i2),vertex(i3),srf(isrf)%ivtx(ivtx),
c     .              nmatch
                ENDDO
                IF (nmatch.EQ.2) THEN
                  IF (i1.NE.0) 
     .              CALL ER('DAMN','SUPER DAMN',*99)
                  i1 = i4
                ENDIF
              ENDDO 
            ENDIF
            IF (output) WRITE(fp,*) '  ==SURFACE ORDER:',iside,i1
            surface(1) = obj(iobj)%iside(i1)
            obj_iside(iside) = i1
              IF (output) WRITE(fp,*) '   s:',0,0,0,
     .          srf( ABS(obj(iobj)%iside(i1)) )%ivtx(1:3)
            newsrf = srf(ABS(surface(1)))
            DO i2 = 1, 3
              newsrf%index(IND_SURFACE) = 0  ! Clear references for internal surfaces
              newsrf%index(IND_TARGET ) = 0
              newsrf%ivtx(1:3) = vertex(s(1:3,3*(iside-1)+i2))
              IF (output) WRITE(fp,*) '   s:',s(1:3,3*(iside-1)+i2),
     .                                vertex(s(1:3,3*(iside-1)+i2))
              surface(i2+1) = AddSurface(newsrf)
            ENDDO

c           Divide the tetrahedron about the centroid:
            newobj = obj(object(1))
c            IF (output) WRITE(fp,*) 'OBJ:',t(1:4,i1)
            newobj%segment(1) = obj_segment
            newobj%iside(1:4) = surface(1:4)
            object(i1+1) = AddObject(newobj)
            IF (output) WRITE(fp,*) '   ISIDE:',newobj%iside(1:4)
          ENDDO
          IF (output) WRITE(fp,*) 'OBJECTS:',object(1:5)

c           STOP 'KNARLY'

c...      Check if each tetrahedron can be further subdivided with the goal of
c         reducing the maximum length of the eventual sides. Poor explanation, sigh:
          vertex_centroid = vertex(5)
          DO iside = 1, 4 ! 4       ! 1, 4   
            omap = obj(iobj)%omap(obj_iside(iside))
            smap = obj(iobj)%smap(obj_iside(iside))
c            omap = obj(object(1))%omap(iside)
c            smap = obj(object(1))%smap(iside)
            
            IF (output) WRITE(fp,*) 'MAP:',iside,obj_iside(iside),
     .                              obj(iobj)%omap(1:4)

            IF (omap.EQ.0.OR.
     .          obj(MAX(1,omap))%segment(1).NE.0) THEN ! *** HACK: NEED TO DECOUPLE THIS FROM THE DELETE FLAG... ***
c     .          obj(MAX(1,omap))%segment(1).EQ.1) THEN 
c...          Subdivide:

              obj(object(iside+1))%segment(1) = -1  ! Tagged for deletion

              DO i1 = 1, 3  ! 3  ! Add three tetrahedons

                IF (output) WRITE(fp,*) '=SUBDIVIDE LEVEL 1:',
     .                                  iside,object(iside+1),i1,nobj

c...            Collect list of vertices - note the complete rebuild on each iteration:
c                isrf = ABS(obj(object(iside+1))%iside(1))        ! Points 1 to 3
c                vertex(1:3) = srf(isrf)%ivtx(1:3)
                isrf = obj(object(iside+1))%iside(1)        ! Points 1 to 3
                IF (isrf.GT.0) THEN
                  vertex(1:3) = srf(isrf)%ivtx(1:3)
                ELSE
                  vertex(1) = srf(-isrf)%ivtx(3)  
                  vertex(2) = srf(-isrf)%ivtx(2)  ! VTX 1 and 2 must always be on the 
                  vertex(3) = srf(-isrf)%ivtx(1)  ! along the line segments of the original 
                ENDIF                            ! tetrahedron (on the outside of the new ones...)
                vertex(4) = vertex_centroid
c                isrf      = ABS(obj(object(iside+1))%iside(2))   ! Point 4
c                vertex(4) = srf(isrf)%ivtx(2)
                isrf = ABS(obj(object(iside+1))%iside(1))
                CALL CalcCentroid(isrf,1,p)                        ! Point 5
                vertex(5) = AddVertex(p)

                IF (output) THEN
                  WRITE(fp,*) '  VERTEX:',vertex(1:5)
                  WRITE(fp,*) '  P:',p
                ENDIF

                newsrf = srf(ABS(obj(iobj)%iside(obj_iside(iside))))
c                newsrf = srf(ABS(obj(object(1))%iside(iside)))
                surface = 0
                DO i2 = 1, 4
                  newsrf%ivtx(1:3) = vertex(u(1:3,i2+4*(i1-1)))
                  IF (output) THEN 
                    WRITE(fp,*) '  U:',i2+4*(i1-1),u(1:3,i2+4*(i1-1))
                    WRITE(fp,*) '   ===>',newsrf%ivtx(1:newsrf%nvtx)
                  ENDIF
                  IF (i2.GT.1) THEN
                    newsrf%index(IND_SURFACE) = 0  ! Clear references for internal surfaces
                    newsrf%index(IND_TARGET ) = 0
                  ENDIF
                  IF (i2.EQ.2) THEN
c...                This surface already exists and will be lost when the original tetrahedron
c                   is deleted, so make an explicit assignment of this surface rather than
c                   trying to create another, which will fail because the new surface will have the 
c                   wrong orientation:
                    newsrf%svtx = SUM(newsrf%ivtx(1:newsrf%nvtx))
                    surface(i2) = 0
                    DO i3 = 1, obj(object(iside+1))%nside
                      isrf = obj(object(iside+1))%iside(i3)
                      IF (newsrf%svtx.EQ.srf(ABS(isrf))%svtx) 
     .                  surface(i2) = isrf
                      IF (output) WRITE(fp,*) '   ->',
     .                  newsrf%svtx,srf(ABS(isrf))%svtx
                    ENDDO
                    IF (surface(i2).EQ.0)  THEN
                      IF (output) 
     .                  WRITE(fp,*) '  OBJ:',iside+1,object(iside+1)
                      DO i3 = 1, obj(object(iside+1))%nside
                        isrf = obj(object(iside+1))%iside(i3)
                        IF (output) 
     .                    WRITE(fp,*) '  VTX:',srf(ABS(isrf))%ivtx(1:3)
                      ENDDO
                      DO i4 = 1, nsrf
                        IF (output) 
     .                    WRITE(fp,*) '   ISRF,SUM:',i4,srf(i4)%svtx
                      ENDDO
                      CALL ER('SHIT','DOUBLE SHIT',*99)
                    ENDIF
c                    surface(i2) = obj(object(iside+1))%iside(i1+1)
                  ELSE
                    surface(i2) = AddSurface(newsrf)
                  ENDIF
                ENDDO
                newobj = obj(object(1))
                newobj%segment(1) = obj_segment
                newobj%omap = 0                    ! Clear connection map settings
                newobj%smap = 0
                newobj%iside(1:4) = surface(1:4)
                idum1 = AddObject(newobj)
                IF (output) WRITE(fp,*) '  SURFACE:',idum1,surface(1:4)

c...            Check if this tetrahedron can be sub-divided:
                subdivide = .TRUE.

                iobj1  = iobj
                iside1 = iside
                omap1  = omap 
                smap1  = smap 

                count = 1
                ivtx1 = srf(ABS(obj(nobj)%iside(1)))%ivtx(1)
                ivtx2 = srf(ABS(obj(nobj)%iside(1)))%ivtx(2)
c                ivtx1 = srf(ABS(surface(1)))%ivtx(1)
c                ivtx2 = srf(ABS(surface(1)))%ivtx(2)
                change = 0
                finished = .FALSE.

c                IF ((iobj.EQ.1.AND.iside.NE.2).OR.
c     .              (iobj.EQ.2))
c     .              (iobj.EQ.2.AND.iside.NE.1).OR.
c     .              (i1.NE.1)) 
c     .            subdivide = .FALSE.  ! *** DEBUG ***
c                 subdivide = .FALSE. 


                IF (subdivide) THEN
                  IF (output) THEN
                    WRITE(fp,*) '=== SUBDIVIDE LEVEL 2:',
     .                          iside,object(iside+1)
                    WRITE(fp,*) '  START TETR:',iobj1,iside1
                    WRITE(fp,*) '            :',omap1,smap1
                    WRITE(fp,*) '            :',ivtx1,ivtx2
                  ENDIF
                ENDIF

                DO WHILE (.NOT.finished.AND.subdivide)  ! *** THIS IS WRONG... *** 

                  IF     (omap1.EQ.0) THEN
c                   Hit a wall, need to reverse direction.  If stuck between two 
c                   wall surfaces then just let things bounce back and forth:                    

                    CALL FindConnectedTetrahedron
     .                     (found,omap1,smap1,iside2,
     .                      iobj1,iside1,ivtx1,ivtx2)

                    IF (.NOT.found) 
     .                CALL ER('DivideTetrahedron','Unable to find '//
     .                        'corresponding side',*99)

                    IF (output) THEN
                      WRITE(fp,*) '  WALL FOUND:',iobj1,iside1
                      WRITE(fp,*) '            :',iside2
                      WRITE(fp,*) '            :',omap1,smap1,change
                    ENDIF

                    change = change + 1
                    IF (change.EQ.2) finished = .TRUE.

                  ELSEIF (obj(omap1)%segment(1).NE.0) THEN   ! *** HACK: NEED TO DECOUPLE THIS FROM THE DELETE FLAG... ***
c                  ELSEIF (obj(omap1)%segment(1).EQ.1) THEN   
c                   Identify common line segment between the current side
c                   and the neighbour:
                    iobj1  = omap1
                    iside1 = smap1

                    CALL FindConnectedTetrahedron
     .                     (found,omap1,smap1,iside2,
     .                      iobj1,iside1,ivtx1,ivtx2)

c...                Need this line because...
                    IF (omap1.EQ.0) iside1 = iside2

                    IF (.NOT.found) 
     .                CALL ER('DivideTetrahedron','Unable to find '//
     .                        'neighbouring object',*99)
                    IF (omap1.EQ.object(1).AND.change.EQ.0) 
     .                finished = .TRUE.

                    IF (output) THEN
                      WRITE(fp,*) '  NEXT TETRA:',iobj1,iside1
                      WRITE(fp,*) '            :',omap1,smap1
                    ENDIF
                  ELSE
                    IF (output) THEN
                      WRITE(fp,*) '  SORRY...  :',iobj1,iside1,
     .                            obj(omap1)%segment(1)
                      WRITE(fp,*) '            :',omap1,smap1
                    ENDIF
                    subdivide = .FALSE.
                  ENDIF
                  count = count + 1

                  IF (output) 
     .              WRITE(fp,*) '  SUBDIVIDE STATUS:',finished,subdivide
                ENDDO


                IF (subdivide) THEN

                  subdivision = .TRUE.

                  iobj1 = nobj

                  obj(iobj1)%segment(1) = -1  ! Tag for deletion...

c...              Collect verticies from current object (just added):
                  isrf = obj(iobj1)%iside(1)                
                  IF (isrf.GT.0) THEN
                    vertex(1:3) = srf(isrf)%ivtx(1:3)
                  ELSE
                    isrf = ABS(isrf)
                    vertex(1) = srf(isrf)%ivtx(2)
                    vertex(2) = srf(isrf)%ivtx(1)
                    vertex(3) = srf(isrf)%ivtx(3)
                  ENDIF
                  isrf = ABS(obj(iobj1)%iside(2))              
                  vertex(4) = srf(isrf)%ivtx(2)
                  p(1:3) = 0.5D0 * (vtx(1:3,vertex(1)) + 
     .                              vtx(1:3,vertex(2)))  
                  vertex(5) = AddVertex(p)

c...              Divide into 2 tetrahedrons:
                  DO i2 = 1, 2  
                    newsrf = srf(ABS(obj(iobj)%iside(obj_iside(iside))))
c                    newsrf = srf(ABS(obj(object(1))%iside(iside)))
                    surface = 0
                    DO i3 = 1, 4
                      IF (output) 
     .                  WRITE(fp,*) '  V:',i3+4*(i2-1),
     .                              v(1:3,i3+4*(i2-1))
                      IF     (i2.EQ.1.AND.i3.EQ.4) THEN
c...                    (See note above:)
                        surface(i3) = obj(iobj1)%iside(4)
                      ELSEIF (i2.EQ.2.AND.i3.EQ.3) THEN
c...                    (See note above:)
                        surface(i3) = obj(iobj1)%iside(3)
                      ELSE
                        IF (i3.GT.1) THEN
                          newsrf%index(IND_SURFACE) = 0  ! Clear references for new internal surfaces
                          newsrf%index(IND_TARGET ) = 0
                        ENDIF
                        newsrf%ivtx(1:3) = vertex(v(1:3,i3+4*(i2-1)))
                        surface(i3) = AddSurface(newsrf)
                      ENDIF
                    ENDDO
                    newobj = obj(object(1))
                    newobj%segment(1) = obj_segment
                    newobj%omap = 0                    
                    newobj%smap = 0
                    newobj%iside(1:4) = surface(1:4)
                    idum1 = AddObject(newobj)
                    IF (output) 
     .                WRITE(fp,*) '  SURFACE:',idum1,surface(1:4)
                  ENDDO  

                ENDIF

              ENDDO  ! i1 LOOP



            ENDIF
          ENDDO


c...     Need option now to futher subdivide the original surfaces into 3, and produce more tetrahedrons!

c  c       obj(iobj)% =     

      IF (subdivision) THEN
        IF (output) THEN
          WRITE(fp,*) '    ===================='
          WRITE(fp,*) '      YES, SUBDIVISION'
          WRITE(fp,*) '    ===================='
        ENDIF
      ENDIF


      ENDSELECT

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: BuildConnectionMap
c
c
      SUBROUTINE BuildConnectionMap_New
      USE mod_eirene06_locals
      IMPLICIT none

      INTEGER fp,iobj,iobj1,iside,iside1,isrf
      INTEGER i1,i2

      INTEGER                 nlist,ilist
      INTEGER, ALLOCATABLE :: list(:)

      INTEGER count_srf_obj,count_srf_side  ! *** TEMP ***
      count_srf_obj = 0
      count_srf_side = 0

      fp = 0

      WRITE(fp,*) '=== BUILDING CONNECTION MAP ==='

      CALL CalcDerivedQuantity(MODE_SRF_OBJ)
      CALL CalcDerivedQuantity(MODE_SRF_SIDE)

 
c      WRITE(0,*) 'OBJECT LISTING:'
c      DO iobj = 1, nobj
c        WRITE(0,*) iobj,obj(iobj)%iside(1:4),obj(iobj)%segment(1)
c      ENDDO
c      WRITE(0,*) 'DONE'
  
      ALLOCATE(list(nobj))  ! Not really needed now, but will be...
      nlist = 0
      DO iobj = 1, nobj
        IF (obj(iobj)%segment(1).EQ.-1) CYCLE
        nlist = nlist + 1
        list(nlist) = iobj 
      ENDDO
c      WRITE(0,*) 'LIST:',list(1:nlist)

      DO ilist = 1, nlist
        iobj = list(ilist)
        obj(iobj)%omap = 0
        obj(iobj)%smap = 0
      ENDDO

      DO ilist = 1, nlist
        iobj = list(ilist)
        DO iside = 1, obj(iobj)%nside         
          isrf = obj(iobj)%iside(iside)   
          IF (isrf.LT.0) THEN            
            iobj1  = srf_obj (-isrf)     
            iside1 = srf_side(-isrf)     
c            iobj1  = srf(-isrf)%obj     
c            iside1 = srf(-isrf)%side    
                                         
c            WRITE(0,*) '   -->',iobj,iside,isrf,iobj1,iside1

            IF (iobj1 .NE.srf(-isrf)%obj) THEN
              count_srf_obj = count_srf_obj + 1
c              WRITE(0,*) ' -->',-isrf,iobj1,srf(-isrf)%obj
c              STOP 'BAD SRF_OBJ'
            ENDIF
            IF (iside1.NE.srf(-isrf)%side) THEN
              count_srf_side = count_srf_side + 1
c              WRITE(0,*) ' SRF_SIDE ',-isrf,iside1,srf(-isrf)%side
c              WRITE(0,*) '          ',obj(iobj1)%iside
c              STOP 'BAD SRF_SIDE'
            ENDIF
 
            IF (iobj1.EQ.0.OR.iside1.EQ.0) THEN
              WRITE(0,*) '  TROUBLE IOBJ1:',iobj ,iside
              WRITE(0,*) '               :',iobj1,iside1
              WRITE(0,*) '               :',-isrf       
              DO i1 = 1, nobj                           
                DO i2 = 1, obj(i1)%nside
                  IF (obj(i1)%iside(i2).EQ. isrf) THEN
                    WRITE(0,*) '    -A->',i1,i2, isrf,obj(i1)%segment(1)
                  ENDIF
                  IF (obj(i1)%iside(i2).EQ.-isrf) THEN
                    WRITE(0,*) '    -B->',i1,i2,-isrf,obj(i1)%segment(1)
                  ENDIF
                ENDDO
              ENDDO 
              STOP 'HALTING PROGRAM'
            ENDIF


            IF (obj(iobj )%omap(iside ).NE.0.OR.
     .          obj(iobj1)%omap(iside1).NE.0) 
     .        CALL ER('BuildConnectionMap','Multiple assignments',*99)
            obj(iobj )%omap(iside ) = iobj1
            obj(iobj )%smap(iside ) = iside1
            obj(iobj1)%omap(iside1) = iobj 
            obj(iobj1)%smap(iside1) = iside
          ENDIF
        ENDDO
      ENDDO

      WRITE(fp,*) '  COUNT_SRF_OBJ   = ',count_srf_obj ,' (DEBUG)'
      WRITE(fp,*) '  COUNT_SRF_SIDE  = ',count_srf_side,' (DEBUG)'

      WRITE(fp,*) '  ORPHAN SURFACES = ',COUNT(srf_obj.EQ.0),' OF',nsrf

      WRITE(fp,*) '=== DONE ==='

      CALL ClearDerivedQuantity(MODE_SRF_OBJ)
      CALL ClearDerivedQuantity(MODE_SRF_SIDE)
      DEALLOCATE(list)

      RETURN
 99   WRITE(fp,*) '  IOBJ ,ISIDE =',iobj ,iside
      WRITE(fp,*) '  IOBJ1,ISIDE1=',iobj1,iside1
      STOP
      END
c
c ======================================================================
c
c subroutine: CalcDerivedQuantity
c
c
      SUBROUTINE CalcDerivedQuantity(mode)
      USE mod_geometry
      IMPLICIT none

      INTEGER, INTENT(IN) :: mode

      INTEGER iobj,iside,isrf

      CALL ClearDerivedQuantity(mode)

      SELECTCASE(mode)     

        CASE (MODE_OBJ_CENTRE)
          ALLOCATE(obj_centroid(3,nobj))
          DO iobj = 1, nobj
            CALL CalcCentroid(iobj,2,obj_centroid(1,iobj))
          ENDDO

        CASE (MODE_OBJ_VOLUME)
          ALLOCATE(obj_volume(nobj))
          obj_volume = 0.0D0
          DO iobj = 1, nobj
            IF (grp(obj(iobj)%group)%type.NE.GRP_TETRAHEDRON) CYCLE
            CALL gmCalcTetrahedronVolume(iobj)
          ENDDO

        CASE (MODE_OBJ_TUBE)
          ALLOCATE(obj_tube(nobj))
          obj_tube = 0
        CASE (MODE_OBJ_DISTANCE)
          ALLOCATE(obj_distance(2,nobj))
          obj_distance = 0.0D0

        CASE (MODE_SRF_OBJ)
          ALLOCATE(srf_obj(nsrf))
          srf_obj = 0
          DO iobj = 1, nobj
            IF (obj(iobj)%segment(1).EQ.-1) CYCLE  ! Object tagged for deletion...
            DO iside = 1, obj(iobj)%nside
              isrf = obj(iobj)%iside(iside)
              IF (isrf.GT.0) THEN
                IF (srf_obj(isrf).NE.0) THEN
                  WRITE(0,*) '  ISIDE,ISRF    = ',iside,isrf
                  WRITE(0,*) '  CURRENT VALUE = ',srf_obj(isrf)
                  WRITE(0,*) '  NEW     VALUE = ',iobj
                  CALL ER('CalcDerivedQuantity','Duplicate SRF_OBJ',*99)     
                ELSE
                  srf_obj(isrf) = iobj
                ENDIF
              ENDIF
            ENDDO
          ENDDO

        CASE (MODE_SRF_SIDE)
          ALLOCATE(srf_side(nsrf))
          srf_side = 0
          DO iobj = 1, nobj
            IF (obj(iobj)%segment(1).EQ.-1) CYCLE  ! Object tagged for deletion...
            DO iside = 1, obj(iobj)%nside  
              isrf = obj(iobj)%iside(iside)
              IF (isrf.GT.0) THEN
                IF (srf_side(isrf).NE.0) THEN
                  CALL ER('CalcDerivedQuantity','Duplicate '//
     .                    'SRF_SIDE',*99)     
                ELSE
                  srf_side(isrf) = iside
                ENDIF
              ENDIF
            ENDDO
          ENDDO

        CASE DEFAULT
          CALL ER('CalcDerivedQuantity','Unrecognised MODE',*99)
      ENDSELECT

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: CalcDerivedQuantity
c
c
      SUBROUTINE ClearDerivedQuantity(mode)
      USE mod_geometry
      IMPLICIT none

      INTEGER, INTENT(IN) :: mode

      INTEGER fp
     
      fp = 0

      SELECTCASE(mode)     
        CASE (MODE_SRF_OBJ)
          IF (ALLOCATED(srf_obj)) THEN
            WRITE(fp,*) 'DEALLOCATING SRF_OBJ'
            DEALLOCATE(srf_obj)
          ENDIF
        CASE (MODE_SRF_SIDE)
          IF (ALLOCATED(srf_side)) THEN
            WRITE(fp,*) 'DEALLOCATING SRF_SIDE'
            DEALLOCATE(srf_side)
          ENDIF
        CASE (MODE_OBJ_CENTRE)
          IF (ALLOCATED(obj_centroid)) THEN
            WRITE(fp,*) 'DEALLOCATING OBJ_CENTROID'
            DEALLOCATE(obj_centroid)
          ENDIF
        CASE (MODE_OBJ_VOLUME)
          IF (ALLOCATED(obj_volume)) THEN
            WRITE(fp,*) 'DEALLOCATING OBJ_VOLUME'
            DEALLOCATE(obj_volume)
          ENDIF
        CASE (MODE_OBJ_TUBE)
          IF (ALLOCATED(obj_tube)) DEALLOCATE(obj_tube)
        CASE (MODE_OBJ_DISTANCE)
          IF (ALLOCATED(obj_distance)) DEALLOCATE(obj_distance)
        CASE DEFAULT
          CALL ER('ClearDerivedQuantity','Unrecognised MODE',*99)
      ENDSELECT

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: SetupSimpleTetrahedrons
c
c For testing.
c
      SUBROUTINE SetupSimpleTetrahedrons
      USE mod_eirene06_locals
      IMPLICIT none

      TYPE(type_srf) newsrf

      nvtx = 4
      vtx(1,1) = 1.0D0
      vtx(2,1) = 0.0D0
      vtx(3,1) = 0.0D0

      vtx(1,2) = 0.3D0
      vtx(2,2) = 0.0D0
      vtx(3,2) = 0.5D0

      vtx(1,3) = 0.0D0
      vtx(2,3) = 0.0D0
      vtx(3,3) = 0.0D0

      vtx(1,4) = 0.3D0
      vtx(2,4) = 0.5D0
      vtx(3,4) = 0.25D0

      nobj = 1
      obj(1)%omap(1:4) = 0
      obj(1)%smap(1:4) = 0
      obj(1)%x = 0.4D0
      obj(1)%y = 0.25D0
      obj(1)%z = 0.25D0

      nsrf = 0
      newsrf%obj  = nobj
      newsrf%side = 1
      newsrf%nvtx = 3
      newsrf%ivtx(1) = 1
      newsrf%ivtx(2) = 2
      newsrf%ivtx(3) = 3
      obj(nobj)%iside(1) = AddSurface(newsrf)
      newsrf%side = 2
      newsrf%ivtx(1) = 1
      newsrf%ivtx(2) = 4
      newsrf%ivtx(3) = 2
      obj(nobj)%iside(2) = AddSurface(newsrf)
      newsrf%side = 3
      newsrf%ivtx(1) = 2
      newsrf%ivtx(2) = 4
      newsrf%ivtx(3) = 3
      obj(nobj)%iside(3) = AddSurface(newsrf)
      newsrf%side = 4
      newsrf%ivtx(1) = 3
      newsrf%ivtx(2) = 4
      newsrf%ivtx(3) = 1
      obj(nobj)%iside(4) = AddSurface(newsrf)

      WRITE(0,*) 'OBJ:',obj(1)%iside(1:4)

      IF (.TRUE.) THEN
        nvtx = 5
        vtx(1,5) = 0.8D0
        vtx(2,5) = 0.50D0
        vtx(3,5) = 0.75D0

        nobj = 2
        obj(nobj)%omap(1:4) = 0
        obj(nobj)%smap(1:4) = 0
c        obj(nobj)%x = 0.90
c        obj(nobj)%y = 0.35
c        obj(nobj)%z = 0.25
        obj(nobj)%x = 0.25D0 * (vtx(1,1) + vtx(1,2) +
     .                          vtx(1,4) + vtx(1,5))
        obj(nobj)%y = 0.25D0 * (vtx(2,1) + vtx(2,2) +
     .                          vtx(2,4) + vtx(2,5))
        obj(nobj)%z = 0.25D0 * (vtx(3,1) + vtx(3,2) +
     .                          vtx(3,4) + vtx(3,5))

        newsrf%obj  = nobj
        newsrf%nvtx = 3
        newsrf%side = 1
        newsrf%ivtx(1) = 1
        newsrf%ivtx(2) = 2
        newsrf%ivtx(3) = 4
        obj(nobj)%iside(1) = AddSurface(newsrf)
        newsrf%side = 2
        newsrf%ivtx(1) = 1
        newsrf%ivtx(2) = 5
        newsrf%ivtx(3) = 2
        obj(nobj)%iside(2) = AddSurface(newsrf)
        newsrf%side = 3
        newsrf%ivtx(1) = 2
        newsrf%ivtx(2) = 5
        newsrf%ivtx(3) = 4
        obj(nobj)%iside(3) = AddSurface(newsrf)
        newsrf%side = 4
        newsrf%ivtx(1) = 4
        newsrf%ivtx(2) = 5
        newsrf%ivtx(3) = 1
        obj(nobj)%iside(4) = AddSurface(newsrf)

        WRITE(0,*) 'OBJ:',obj(2)%iside(1:4)
      ENDIF

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c  subroutine: CleanSurfaceArray
c
c
      SUBROUTINE CleanSurfaceArray
      USE mod_geometry
      IMPLICIT none

      INTEGER fp,isrf,mval
c      INTEGER, ALLOCATABLE :: list_shift(:)

      fp = 0

      WRITE(fp,*) '=== REMOVING ORPHAN AND DUPLICATE SURFACES ==='
   
c...  Calculate vertex sum:
      DO isrf = 1, nsrf
        srf(isrf)%svtx = SUM(srf(isrf)%ivtx(1:srf(isrf)%nvtx))
      ENDDO

      mval = MAXVAL(srf(1:nsrf)%svtx)

      WRITE(0,*) 'MVAL:',mval
      STOP 'sdfsd'


      WRITE(fp,*) '=== DONE ==='

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c  subroutine: CleanObjectArray
c
c
      SUBROUTINE CleanObjectArray
      USE mod_geometry
      IMPLICIT none

      INTEGER fp,iobj,ndelete,iobj_shift,shift
cifirst,ilast,origin
      INTEGER, ALLOCATABLE :: list_shift(:)

      fp = 0

      WRITE(fp,*) '=== REMOVING OBJECTS MARKED FOR DELETION ==='
   
      ndelete = 0
      DO iobj = 1, nobj
        IF (obj(iobj)%segment(1).EQ.-1) ndelete = ndelete + 1
      ENDDO

      IF (ndelete.GT.0) THEN

        ALLOCATE(list_shift(nobj))

c...    Build list of objects to be deleted:
        shift = 0
        iobj = 0
        iobj_shift = 0
        DO WHILE(iobj_shift.LT.nobj)
          iobj_shift = iobj_shift + 1
          IF (obj(iobj_shift)%segment(1).EQ.-1) THEN 
            shift = shift + 1
          ELSE
            iobj = iobj + 1
            list_shift(iobj) = iobj + shift
          ENDIF
        ENDDO
 
        WRITE(0,*) 'SUMS:',nobj,ndelete,iobj,shift

c...    Delete objects:
        nobj = nobj - ndelete
        DO iobj = 1, nobj
c          WRITE(0,*) '  DELETING IOBJ,SHIFT:',iobj,list_shift(iobj)
          IF (iobj.NE.list_shift(iobj)) obj(iobj)=obj(list_shift(iobj))
        ENDDO

c...    Remap:
c        DO iobj = 1, nobj
c          DO iside = 1, obj(iobj)%nside
c            obj(iobj)%omap(iside) = list_shift(obj(iobj)%omap(iside))
c          ENDDO
c        ENDDO

        DEALLOCATE(list_shift)
      ENDIF

c...  Make sure all the plasma cells are at the beginning of the array:
c      ifirst = 0
c      ilast = 0
c      DO iobj = 1, nobj
c        origin = grp(obj(iobj)%group)%origin
c        IF (ifirst.EQ.0.AND.origin.NE.GRP_MAGNETIC_GRID) THEN  ! *** NEED A BETTER 
c          ifirst = iobj                                        ! DETECTOR SINCE IN FUTURE
c        ENDIF                                                  ! SOME NON MAG GRID CELS WILL
c        IF (origin.EQ.GRP_MAGNETIC_GRID) THEN                  ! HAVE PLASMA ***
c          ilast = iobj
c        ENDIF                                                  ! *** THIS IS AN ISSUE FOR OSM MORE THAN HERE, WHERE IN OSM
c      ENDDO                                                    ! IT IS OFTEN ASSUMED THAT THERE IS A 1:1 BETWEEN NCELL AND NOBJ 
c      WRITE(fp,*) '  FIRST VACUUM, LAST PLASMA CELL=',ifirst,ilast,nobj

      WRITE(fp,*) '=== DONE ==='

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c  subroutine: FixTetrahedrons
c
c  This has to do with the order that the sides appear in the OBJ%SIDE
c  listing, rather than the ordering of vertices.
c
c
      SUBROUTINE FixTetrahedrons(istart,iend)
      USE mod_eirene06_locals
      USE mod_eirene06
      USE mod_geometry
      IMPLICIT none

      INTEGER, INTENT(IN) :: istart,iend

      INTEGER iobj,iside,isrf,ivtx,iobj1,iside1,i1,i2,i3,count,
     .        ipts(3,4),imap(3)
      LOGICAL checking
      TYPE(type_object) tmpobj


      WRITE(eirfp,*) '    SORTING SIDES TO BASE'      
      checking = .FALSE.

 10   CONTINUE

      DO iobj = istart, iend
c...    Check if tetrahedron base is ... :
        isrf = obj(iobj)%iside(1)
        IF     (isrf.EQ.0) THEN
          CALL ER('BuildConnectionMap','Problem with tetrahedron',*99)
        ELSE
c        ELSEIF (isrf.LT.0) THEN
c...      Eirene assumes that the base of the tetrahedron is ordered
c         counter clock-wise when viewed from the outside and that
c         tetrahedron sides 2-4 correspond to triangle sides 1-3 on the
c         base (there is also an assumption about vertex referencing
c         for the sides, so that the tetrahedron can be described by
c         just 4 vertices, but that has to be handled on the fly 
c         when the data transfer files are created in WriteEireneObjects):

          IF (isrf.GT.0) THEN
            ipts(1:3,1) = srf(isrf)%ivtx(1:3)  
          ELSE
            ipts(1,1) = srf(-isrf)%ivtx(3)  ! This convention must be 
            ipts(2,1) = srf(-isrf)%ivtx(2)  ! consistent with what is done
            ipts(3,1) = srf(-isrf)%ivtx(1)  ! in WriteEireneObjects.
          ENDIF      
  
          DO iside = 2, obj(iobj)%nside         
            isrf = ABS(obj(iobj)%iside(iside)) 
            ipts(1:3,iside) = srf(isrf)%ivtx(1:3)
          ENDDO

          imap = 0
          DO i1 = 1, 3
            i2 = i1 + 1
            IF (i2.EQ.4) i2 = 1
c...        Map tetrahedron side to corresponding base triangle side:
            DO iside = 2, obj(iobj)%nside
              count = 0
              DO i3 = 1, 3
                IF (ipts(i1,1).EQ.ipts(i3,iside).OR.
     .              ipts(i2,1).EQ.ipts(i3,iside)) count = count + 1
              ENDDO
c              WRITE(eirfp,*) 'COUNTING:',iobj,i1,iside,count
              IF     (count.EQ.2) THEN
                imap(i1) = iside
              ELSEIF (count.EQ.0) THEN
                CALL ER('FixTetrahedrons','Something wicked',*99)
              ENDIF
            ENDDO
          ENDDO

          tmpobj = obj(iobj)

          DO i1 = 1, 3
            IF (imap(i1).NE.i1+1) THEN
c...          Need to move the tetrahedron side so that it matches up
c             with the proper base triangle side, and also update
c             the connection map:

              IF (checking) WRITE(eirfp,*) '* SWAPPING SHIT!:',
     .                                     i1+1,imap(i1)

              obj(iobj)%iside(i1+1) = tmpobj%iside(imap(i1))
              obj(iobj)%omap (i1+1) = tmpobj%omap (imap(i1))
              obj(iobj)%smap (i1+1) = tmpobj%smap (imap(i1))
c...          Update reverse map:
              iobj1  = obj(iobj)%omap(i1+1)
              iside1 = obj(iobj)%smap(i1+1)

              obj(iobj1)%omap(iside1) = iobj   
              obj(iobj1)%smap(iside1) = i1+1
            ENDIF
          ENDDO

        ENDIF
      ENDDO

      IF (.NOT.checking) THEN
        WRITE(eirfp,*) '    PASS 2...'
        checking = .TRUE.
        GOTO 10
      ENDIF

 20   CONTINUE

      WRITE(eirfp,*) '    CHECKING MAP'
c...  Check for problems:
      DO iobj = istart, iend
        DO iside = 1, obj(iobj)%nside         
          IF (iside.NE.1.AND.obj(iobj)%omap(iside).EQ.0) THEN
            WRITE(eirfp,*) 'PROBLEM WITH MAP: ISIDE.NE.1',iobj,iside
          ENDIF
        ENDDO
      ENDDO

      WRITE(eirfp,*) '  DONE'

      RETURN
 99   STOP
      END
