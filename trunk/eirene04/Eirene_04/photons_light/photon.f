      module photon
      USE PRECISION
      USE PARMMOD
      use COMXS
      use COMUSR
      use CGRID
      use COMPRT
      use CCONA
      use CZT1

      implicit none 

      private

      public :: ph_init, ph_plasma_profile, ph_energy, ph_getcoeff,
     .    ph_xsectph, ph_alloc_xsecta, ph_xsectp, 
     .    ph_post_energy, ph_alloc_xsectph, 
     .    ph_integrate, ph_lorvdw, ph_vdwprof, ph_b21

      integer, public, save :: actual_linenum
      integer, public, save  :: phv_nrota, phv_nrotph, phv_muldens
      integer, public, save, allocatable :: 
     .   PHV_LGAOT(:,:,:), PHV_LGPHOT(:,:,:), 
     .   PHV_NAOTI(:), PHV_NPHOTI(:), 
     .   PHV_N1STOTAT(:,:,:), PHV_N2NDOTAT(:,:,:),
     .   PHV_IESTOTAT(:,:,:), PHV_IESTOTPH(:,:,:),
     .   PHV_N1STOTPH(:,:,:), PHV_N2NDOTPH(:,:,:)

      real(dp), public, save :: CLIGHT,EV2HZ,HPLNK,hpcl
      real(dp), public, save :: hwvdw
    
csw external
      integer, external :: idez, learc1, learc2
      real(dp), external :: ranf_eirene

      contains



      REAL(dp) FUNCTION PH_ENERGY(il,icell,kk,ipl2,nldoppl) result(res)
      IMPLICIT NONE
cdr  sample frequency (here: energy) from emission profile
cdr  iprofiletype =0  only line centre (delta),  stationary atoms
cdr  iprofiletype =1  ditto plus doppler in calling program
cdr  iprofiletype =2  only lorentz (delta),  stationary atoms
cdr  iprofiletype =3  ditto plus doppler in calling program
cdr  iprofiletype =4  lorentz+vdw ,  stationary atoms
cdr  iprofiletype =5  ditto plus doppler in calling program (new, july 04)
cdr
cdr
cdr
cdr
cdr  july 04: nldoppl introduced
cdr             now ph_energy only samples in the rest frame of the
cdr             emitting atom.
cdr             it returns nldoppl=true, if doppler shift is to be
cdr                                      added in calling program
cdr             
      integer, intent(in) :: il,icell,kk,ipl2
      logical, intent(out) :: nldoppl 
      real(dp) :: e00,pnue0,dnd,gam,zep1

!  idreac: kk from last call to get_reaction
      if (idreac /= kk) call get_reaction(kk)
      e00=reaction%e0
c     pnue0=e00*ev2hz

c  lorentz plus doppler
cdr      call ph_voigtprof(ipl2,icell,dnd,gam)
cdr      zep1 = ph_sam_voigt(gam,dnd/sqrt(2.))
cdr      res=e00 + zep1
      call ph_homprof(icell, gam)
      res = ph_sam_lorentz(gam,e00)
      nldoppl=.true.
      weight = 1.

      return
      END FUNCTION PH_ENERGY

      
      SUBROUTINE PH_GETCOEFF(kkin,isp,ity,icell,iipl,idsc,fac,res) 
      IMPLICIT NONE
      integer, intent(in) :: kkin,isp,ity,icell,iipl, idsc
      real(dp), intent(out) :: fac,res
      integer :: iid,kk,ipl,ignd,iptype
      real(dp):: gam,e00,dnd,xx,yy,pnue,pnue0

      if (kkin /= idreac) call get_reaction(kkin)

      kk=kkin

      iid = reaction%ircart
      fac=0._dp
      hwvdw=0._dp

      ipl = reaction%ignd
      
c     P.2 PH_ABS OT
      if(lgvac(icell,ipl)) then
         res=0.
         return
      endif

      e00=reaction%e0
      pnue0 = e00*EV2HZ
      pnue  = e0 *EV2HZ

         
      iptype = reaction%iprofiletype
      call ph_voigtprof(iipl,icell,dnd,gam)
      xx=(e0-e00)/dnd
      yy=gam/dnd
      fac = DBLE(PH_FADDEEVA(xx,yy,dnd,icell,iipl))

      res=e00*reaction%b12
      res=res*fac*clight/(4._dp*PIA)
c result should have units cm**3/s --> multiply backgr. density!
      phv_muldens=1

      
      return
      END SUBROUTINE PH_GETCOEFF


      subroutine ph_voigtprof(ipl,icell,dnd,gam)
      implicit none
      integer, intent(in) :: ipl,icell
      real(dp), intent(out) :: dnd,gam

      call ph_dopplerprof(ipl,icell,dnd)
      call ph_homprof(icell,gam)
      return
      end subroutine ph_voigtprof


      subroutine ph_dopplerprof(iipl,icell,dnd)
      implicit none
      integer, intent(in) :: iipl,icell
      real(dp), intent(out) :: dnd
      real(dp) :: e00,t,e1,e2
      integer :: n1,n2,ipl

      ipl=iipl
      e00=reaction%e0
      e1=reaction%e1
      e2=e00+e1
      t=tiin(mplsti(ipl),icell)     
      dnd=e00*sqrt(t)*rsqdvp(ipl)/clight
      return
      end subroutine ph_dopplerprof


      subroutine PH_HOMPROF(icell, gam)
      IMPLICIT NONE
      integer, intent(in) :: icell
      real(dp), intent(out) :: gam
      real(dp) :: gam1,gam2

      gam=0.
      
c natural
      gam=gam+reaction%aik

      gam=gam/(2.*PIA)/ev2hz
      return
      END subroutine PH_HOMPROF


      real(dp) function ph_lorvdw(de,hw,shift,dvdw,
     .                            icell,ipl) result(res)
      implicit none
      real(dp), intent(in) :: de,hw,shift,dvdw
      integer, intent(in) :: icell,ipl
      res = 0._dp
      write (6,*) ' ph_lorvdw: no calculations for photons possible '
      return
      end function ph_lorvdw


      subroutine ph_plasma_profile(ind)
      implicit none
      integer, intent(in) :: ind
      write (6,*) ' ph_plasma_profile: no calculations for photons',
     .            ' possible '
      return
      end subroutine ph_plasma_profile 


      SUBROUTINE PH_INIT(ICAL)
      implicit none
      integer, intent(in) :: ical

      if (ical == 0) then
         CLIGHT = 2.99792458D10
         EV2HZ  = 2.41798834D14
         HPLNK  = 1./EV2HZ
c     h*c [eV*cm]
         HPCL   = CLIGHT*HPLNK
      elseif (ical == 1) then
         phv_nrota=0
         phv_nrotph=0
      end if

      return
      end subroutine PH_INIT


      subroutine ph_integrate(istr, scale, fpht,fatm,fmol,fion)
      implicit none
      real(dp), intent(in) :: scale,fpht,fatm,fmol,fion
      integer, intent(in) :: istr
      write (6,*) ' ph_integrate: no integration of background photons '
      call leer(1)
      return
      end subroutine ph_integrate


      SUBROUTINE PH_POST_ENERGY(icell,kk,iflg,il,
     .                          iold,itypold,vxo,vyo,vzo,vlo,e0o,itypnw)
      IMPLICIT NONE
      integer, intent(in) :: icell,kk,iflg,iold,il,itypold,itypnw
      real(dp),intent(in) :: vxo,vyo,vzo,vlo,e0o
      write (6,*) ' PH_POST_ENERGY: no calculations for photons',
     .            ' possible '
      return
      end subroutine PH_POST_ENERGY


      SUBROUTINE PH_VDWPROF(icell,hw,shift,dvdw,lscale) 
      implicit none
      integer, intent(in) :: icell
      real(dp),intent(out) :: hw,shift,dvdw
      logical, intent(in) :: lscale
      write (6,*) ' PH_VDWPROF: no calculations for photons possible '
      hw=0._dp
      shift=0._dp
      dvdw=0._dp
      return
      end subroutine PH_VDWPROF


      SUBROUTINE PH_ALLOC_XSECTA(nnrot)
      IMPLICIT NONE
      integer, intent(in) :: nnrot

      allocate(PHV_LGAOT(0:natm,0:nnrot,0:5))
      PHV_NROTA=nnrot
      phv_lgaot=0         
      allocate(PHV_NAOTI(natm))
      phv_naoti=0
      
      allocate(phv_iestotat(0:natm,nnrot,3))
      phv_iestotat=0
      
      allocate(phv_n1stotat(0:natm,nnrot,3))
      allocate(phv_n2ndotat(0:natm,nnrot,3))
      phv_n1stotat=0
      phv_n2ndotat=0
!pb      write (6,*) ' PH_ALLOC_XSECTA: no calculations for photons',
!pb     .            ' possible '
      return
      end subroutine PH_ALLOC_XSECTA


      SUBROUTINE PH_XSECTP(ipls,nrc,idsc,irrc)
      IMPLICIT NONE
      integer, intent(in) :: ipls,nrc,idsc,irrc
      write (6,*) ' PH_XSECTP: no calculations for photons possible '
      return
      end subroutine PH_XSECTP


      SUBROUTINE PH_ALLOC_XSECTPH(nnrot)
      IMPLICIT NONE
      integer, intent(in) :: nnrot

      allocate(PHV_LGPHOT(0:nphot,0:nnrot,0:5))
      phv_nrotph=nnrot
      phv_lgphot=0
      allocate(phv_nphoti(nphot))
      phv_nphoti=0
      allocate(phv_iestotph(0:nphot,nnrot,3))
      phv_iestotph=0
      allocate(phv_n1stotph(0:nphot,nnrot,3))
      allocate(phv_n2ndotph(0:nphot,nnrot,3))
      phv_n1stotph=0
      phv_n2ndotph=0
!pb      write (6,*) ' PH_ALLOC_XSECTPH: no calculations for photons',
!pb     .            ' possible '

      return
      END SUBROUTINE PH_ALLOC_XSECTPH



      SUBROUTINE PH_XSECTPH(ipht,nrc,idsc)
      IMPLICIT NONE
      integer, intent(in) :: ipht,nrc,idsc,ipl
      integer :: kk,ipl0,ipl1,ipl2,ityp0,ityp1,ityp2,il,n0,n1,n2,
     .    nh,nl,iid,ifnd,mode,updf, j, nseot4, ierr, ipl0ti
      real(dp) :: factkk, ebulk

      kk=ireacph(ipht,nrc)
      factkk=freacph(ipht,nrc)
      if(factkk == 0.) factkk=1.

      IPL0 =IDEZ(IBULKPH(ipht,nrc),3,3)
      IPL1 =IDEZ(ISCD1PH(ipht,nrc),3,3)
      IPL2 =IDEZ(ISCD2PH(ipht,nrc),3,3)
      ITYP0=IDEZ(IBULKPH(ipht,nrc),1,3)
      ITYP1=IDEZ(ISCD1PH(ipht,nrc),1,3)
      ITYP2=IDEZ(ISCD2PH(ipht,nrc),1,3)

c collect niveau information data, N0, N1, N2
      
      updf=1
      mode=0
      ifnd=999
      
      PHV_LGPHOT(ipht,idsc,0)=idsc
      PHV_LGPHOT(ipht,idsc,1)=ipl0
      PHV_LGPHOT(ipht,idsc,2)=ifnd
      PHV_LGPHOT(ipht,idsc,3)=kk
      PHV_LGPHOT(ipht,idsc,4)=updf
      PHV_LGPHOT(ipht,idsc,5)=mode

CDR  1ST SECONDARY
      PHV_N1STOTph(ipht,idsc,1) = ityp1
      PHV_N1STOTph(ipht,idsc,2) = ipl1
      PHV_N1STOTph(ipht,idsc,3) = 0
      IF (ityp1 < 4) 
     .  PHV_N1STOTph(ipht,idsc,3) = IDEZ(ISCD1PH(ipht,nrc),2,3)
CDR  2ND SECONDARY
      PHV_N2NDOTph(ipht,idsc,1) = ityp2
      PHV_N2NDOTph(ipht,idsc,2) = ipl2
      PHV_N2NDOTph(ipht,idsc,3) = 0
      IF (ityp2 < 4) 
     .  PHV_N2NDOTph(ipht,idsc,3) = IDEZ(ISCD2PH(ipht,nrc),2,3)
      

cdr  CROSS SECTION: HERE: BEAM-BEAM, NO DOPPLER FROM THERMAL MOTION
      MODCOL(7,1,IPHT,IPL0)=KK
cdr  COLLISION MODEL 4: BEAM-BEAM
      MODCOL(7,2,IPHT,IPL0)=4
c
c     DEFCX(IRCX)=LOG(CVELI2*PMASS)
c     EEFCX(IRCX)=LOG(CVELI2*TMASS)
C
C  3. BULK ION MOMENTUM LOSS RATE
C
C
C  4. BULK ION ENERGY LOSS RATE
C
cdr   NSEOT4=IDEZ(ISCDE,4,5)
cdr bulk energy loss not ready
      nseot4=0
      ebulk=0.
cdr
      IF (NSEOT4.EQ.0) THEN
C  4.A)  ENERGY LOSS RATE OF IMP. BULK ION = CONST.*RATECOEFF.
C        SAMPLE COLLIDING ION FROM DRIFTING MONOENERGETIC ISOTROPIC DISTRIBUTION
        write (6,*) ' in ph_xsectph, nseot4=0 '
        IF (EBULK.LE.0.D0) THEN
          write (6,*) ' in ph_xsectph, nseot4=0, ebulk <= 0 '
          IF (NSTORDR >= NRAD) THEN
            write (6,*) ' in ph_xsectph, nstordr>nrad ',idsc
            IPL0TI = MPLSTI(IPL0)
            DO J=1,NSBOX
              EPLOT3(Idsc,J,1)=1.5*TIIN(IPL0TI,J)+EDRIFT(IPL0,J)
            ENDDO
            NELROT(Idsc) = -3
          ELSE
            write (6,*) ' in ph_xsectph, nstordr<nrad '
            NELROT(Idsc) = -3
          END IF
        ELSE
          write (6,*) ' in ph_xsectph, nseot4=0, ebulk > 0 '
          IF (NSTORDR >= NRAD) THEN
            write (6,*) ' in ph_xsectph, nstordr>nrad '
            DO 251 J=1,NSBOX
              EPLOT3(Idsc,J,1)=EBULK+EDRIFT(IPL0,J)
251         CONTINUE
            NELROT(Idsc) = -2
          ELSE
            NELROT(Idsc) = -2
            EPLOT3(Idsc,1,1)=EBULK
            write (6,*) ' in ph_xsectph, nstordr<nrad '
          END IF
        ENDIF
        MODCOL(7,4,IPHT,IPL0)=3
        write (6,*) ' in ph_xsectph, Modcol(7,4,1,1) ',
     .                MODCOL(7,4,IPHT,IPL0)
      ELSEIF (NSEOT4.EQ.1) THEN
C  4.B) ENERGY LOSS RATE OF IMP. ION = 1.5*TI* RATECOEFF.
C       SAMPLE COLLIDING ION FROM DRIFTING MAXWELLIAN
        write (6,*) ' in ph_xsectph, nseot4=1 '
        IF (EBULK.LE.0.D0) THEN
          IF (NSTORDR >= NRAD) THEN
            IPL0TI = MPLSTI(IPL0)
            DO 252 J=1,NSBOX
              EPLOT3(Idsc,J,1)=1.5*TIIN(IPL0TI,J)+EDRIFT(IPL0,J)
  252       CONTINUE
            NELROT(Idsc) = -3
          ELSE
            NELROT(Idsc) = -3
          END IF
        ELSE
          WRITE (6,*) 'WARNING FROM SUBR. xsectph '
          WRITE (6,*) 'MODIFIED TREATMENT OF photon collision '
          WRITE (6,*) 'SAMPLE FROM MAXWELLIAN WITH T = ',EBULK/1.5
          WRITE (6,*) 'RATHER THEN WITH T = TIIN '
          CALL LEER(1)
          IF (NSTORDR >= NRAD) THEN
            DO 2511 J=1,NSBOX
              EPLOT3(Idsc,J,1)=EBULK+EDRIFT(IPL0,J)
2511        CONTINUE
            NELROT(Idsc) = -2
          ELSE
            NELROT(Idsc) = -2
            EPLOT3(Idsc,1,1)=EBULK
          END IF
        ENDIF
        MODCOL(7,4,IPHT,IPL0)=1
C     ELSEIF (NSECX4.EQ.2) THEN
C  use i-integral expressions. to be written
c     ELSEIF (NSECX4.EQ.3) THEN
C  4.B)  ENERGY LOSS RATE OF IMP. ION = EN.WEIGHTED RATE
C  4.C)  ENERGY LOSS RATE OF IMP. ION = EN.WEIGHTED RATE
      ELSE
        IERR=5
        GOTO 996
      ENDIF

C  ESTIMATOR FOR CONTRIBUTION TO COLLISION RATES FROM THIS REACTION

      PHV_IESTOTph(ipht,idsc,1) = IDEZ(IESTMPH(IPHT,idsc),1,3)
      PHV_IESTOTph(ipht,idsc,2) = IDEZ(IESTMPH(IPHT,idsc),2,3)
      PHV_IESTOTph(ipht,idsc,3) = IDEZ(IESTMPH(IPHT,idsc),3,3)
C
cdr  not ready
c
c     ITYP1=N1STX(IRCX,1)
c     ITYP2=N2NDX(IRCX,1)
c     IF (IESTCX(IRCX,1).NE.0.AND.(ITYP1.NE.1.OR.ITYP2.NE.4)) THEN
c       CALL LEER(1)
c       WRITE (6,*) 'WARNING: COLL.EST NOT AVAILABLE FOR PART.-BALANCE '
c       WRITE (6,*) 'IRCX = ',IRCX
c       WRITE (6,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
c       IESTCX(IRCX,1)=0
c     ENDIF
c     IF (IESTCX(IRCX,2).NE.0.AND.(ITYP1.NE.1.OR.ITYP2.NE.4)) THEN
c       CALL LEER(1)
c       WRITE (6,*) 'WARNING: COLL.EST NOT AVAILABLE FOR MOM.-BALANCE '
c       WRITE (6,*) 'IRCX = ',IRCX
c       WRITE (6,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
c       IESTCX(IRCX,2)=0
c     ENDIF
c     IF (IESTCX(IRCX,3).NE.0.AND.(ITYP1.NE.1.OR.ITYP2.NE.4)) THEN
c       CALL LEER(1)
c       WRITE (6,*) 'WARNING: COLL.EST NOT AVAILABLE FOR EN.-BALANCE '
c       WRITE (6,*) 'IRCX = ',IRCX
c       WRITE (6,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
c       IESTCX(IRCX,3)=0
c     ENDIF
      RETURN
C
      ENTRY XSTPH_2(Idsc,IPL)
C
c     CALL LEER(1)
c     WRITE (6,*) 'Photon REACTION NO. Idsc= ',Idsc
c     CALL LEER(1)
c     WRITE (6,*) 'Collision WITH BULK IONS IPLS:'
c     WRITE (6,*) '1ST AND 2ND NEXT GEN. SPECIES I2ND1, I2ND2:'
c     ITYP1=N1STX(IRCX,1)
c     ITYP2=N2NDX(IRCX,1)
c     ISPZ1=N1STX(IRCX,2)
c     ISPZ2=N2NDX(IRCX,2)
c     IF (ITYP1.EQ.1) TEXTS1=TEXTS(NSPH+ISPZ1)
c     IF (ITYP1.EQ.2) TEXTS1=TEXTS(NSPA+ISPZ1)
c     IF (ITYP1.EQ.3) TEXTS1=TEXTS(NSPAM+ISPZ1)
c     IF (ITYP1.EQ.4) TEXTS1=TEXTS(NSPAMI+ISPZ1)
c     IF (ITYP2.EQ.1) TEXTS2=TEXTS(NSPH+ISPZ2)
c     IF (ITYP2.EQ.2) TEXTS2=TEXTS(NSPA+ISPZ2)
c     IF (ITYP2.EQ.3) TEXTS2=TEXTS(NSPAM+ISPZ2)
c     IF (ITYP2.EQ.4) TEXTS2=TEXTS(NSPAMI+ISPZ2)
c     WRITE (6,*) 'IPLS= ',TEXTS(NSPAMI+IPL),'I2ND1= ',TEXTS1,
c    .                    'I2ND2= ',TEXTS2
c     CALL LEER(1)
      RETURN
C
990   CONTINUE
      WRITE (6,*) 'ERROR IN XSectph. EXIT CALLED  '
      WRITE (6,*) 'INVALID SPECIES INDEX FOR OT '
      CALL EXIT(1)
992   CONTINUE
      WRITE (6,*) 'ERROR IN XSectph: EXIT CALLED  '
      WRITE (6,*) 'MASS NUMBERS OF INTERACTING PARTICLES INCONSISTENT'
!pb   WRITE (6,*) 'KK,IATM,IPLS ',KK,IATM,IPLS
      CALL EXIT(1)
993   CONTINUE
      WRITE (6,*) 'ERROR IN XSectph: EXIT CALLED  '
      WRITE (6,*) 'EBULK_ION .LE.0, BUT MONOENERGETIC DISTRIBUTION?'
      WRITE (6,*) 'CHECK ENERGY FLAG ISCDEA'
!pb   WRITE (6,*) 'KK,IATM,IPLS,ISCDEA ',KK,IATM,IPLS,ISCDEA
      CALL EXIT(1)
996   CONTINUE
      WRITE (6,*) 'ERROR IN XSectph: EXIT CALLED  '
      WRITE (6,*) 'NO CROSS SECTION AVAILABLE FOR NON DEFAULT OT'
      WRITE (6,*) 'KK,IPHT,IPL0 ',KK,IPHT,IPL0
      WRITE (6,*) 'EITHER PROVIDE CROSS SECTION OR USE DIFFERENT '
      WRITE (6,*) 'POST COLLISION SAMPLING FLAG ISCDEA'
      CALL EXIT(1)
      return
      END SUBROUTINE PH_XSECTPH


      REAL(DP) FUNCTION PH_B21() result(res)
      IMPLICIT NONE
c calculates B21 Einstein coefficient in units: cm**2
      
      real(DP) :: e00,pnue0

      e00=reaction%e0         
      
      if(e00 > 0.0) THEN
         pnue0 = E00*EV2HZ
         res=reaction%aik
         res=res * (hplnk*clight)**2 / (2.*e00**3)
c convert [cm^2 / (eV*s) ] --> [cm^2]
         res=res*hplnk
      else
         write(6,*) 'PHOTON MODULE (PH_B21): e00 <= 0.'
         res=0.
      endif
      
      return
      END FUNCTION PH_B21



      REAL(dp) FUNCTION PH_SAM_LORENTZ(alph,totshift) result(res)
      IMPLICIT NONE
      real(dp), intent(in) :: alph, totshift
      real(dp) :: xx
      do
        xx=PIHA*(RANF_EIRENE()*2._DP-1._DP)
        res=alph*tan(xx)+totshift
        if (res > 0._dp) exit
      end do
      return
      END FUNCTION PH_SAM_LORENTZ


      complex(dp) function ph_faddeeva(x,y,dnd,icell,ipl) result(cres)
      implicit none
      real(dp), intent(in) :: x,y,dnd
      integer, intent(in) :: icell,ipl

      cres=ph_humlik(x,y)/(dnd*sqrt(pia))
      return
      end function ph_faddeeva
c
c
c
      complex(dp) FUNCTION PH_HUMLIK(x,y) result(cres)
      IMPLICIT NONE
c arguments
      real(dp), intent(in) :: x,y
c To calculate the FADDEEVA function with relative error less than 10^(-R).
c R0=1.51*EXP(1.144*R) and R1=1.60*EXP(0.554*R) can be set by the the user
c subject to the constraints 14.88<R0<460.4 and 4.85<R1<25.5
      REAL(dp) :: K,L
      real(dp), PARAMETER :: R0 = 146.7, R1 = 14.67 ! for R=4, region boundaries
            
c Constants
      real(dp), PARAMETER :: RRTPI = 0.56418958     ! 1/sqrt(pi)
      real(dp), PARAMETER :: Y0 = 1.5, Y0PY0 = Y0+Y0, Y0Q = Y0*Y0 ! for cpf12 algor.
      REAL(dp), save :: C(0:5), S(0:5), T(0:5)
c SAVE preserves values of C, S and T (static) arrays between procedure calls

      DATA C / 1.0117281,     -0.75197147,        0.012557727,
     .     0.010022008,   -0.00024206814,     0.00000050084806 /
      DATA S / 1.393237,       0.23115241,       -0.15535147,
     .     0.0062183662,   0.000091908299,   -0.00000062752596 /
      DATA T / 0.31424038,     0.94778839,        1.5976826,
     .     2.2795071,      3.0206370,         3.8897249 /
      
c Local variables
      INTEGER :: I, J              ! Loop variables
      INTEGER :: RG1, RG2, RG3     ! y polynomial flags
      REAL(dp) :: ABX, XQ, YQ, YRRTPI ! |x|, x^2, y^2, y/SQRT(pi)
      REAL(dp) :: XLIM0, XLIM1, XLIM2, XLIM3, XLIM4 ! |x| on region boundaries
      REAL(dp) :: A0, D0, D2, E0, E2, E4, H0, H2, H4, H6 ! W4 temporary variables
      REAL(dp) :: P0, P2, P4, P6, P8, Z0, Z2, Z4, Z6, Z8
      real(dp) :: b1,f1,f3,f5,q1,q3,q5,q7
      REAL(dp) :: XP(0:5), XM(0:5), YP(0:5), YM(0:5) ! CPF12 temporary values
      REAL(dp) :: MQ(0:5), PQ(0:5), MF(0:5), PF(0:5)
      REAL(dp) :: D, YF, YPY0, YPY0Q  
      
c**** Start of executable code *****************************************
      
      RG1 = 1                   ! Set flags
      RG2 = 1
      RG3 = 1
      YQ  = Y*Y                 ! y^2
      YRRTPI = Y*RRTPI          ! y/SQRT(pi)
      
c Region boundaries when both K and L are required or when R<>4 
      XLIM0 = R0 - Y 
      XLIM1 = R1 - Y
      XLIM3 = 3.097*Y - 0.45
      
      XLIM2 = 6.8 - Y
      XLIM4 = 18.1*Y + 1.65
      IF ( Y .LE. 0.000001 ) THEN ! When y<10^-6
         XLIM1 = XLIM0          ! avoid W4 algorithm
         XLIM2 = XLIM0
      ENDIF
c.....
      ABX = ABS ( X )           ! |x|
      XQ  = ABX*ABX             ! x^2
      IF ( ABX .GT. XLIM0 ) THEN ! Region 0 algorithm
         K = YRRTPI / (XQ + YQ)
         L = k*x/y        
      ELSEIF ( ABX .GT. XLIM1 ) THEN ! Humlicek W4 Region 1
         IF ( RG1 .NE. 0 ) THEN ! First point in Region 1
            RG1 = 0
            A0 = YQ + 0.5       ! Region 1 y-dependents
            D0 = A0*A0
            D2 = YQ + YQ - 1.0
            b1 = yq - 0.5
         ENDIF
         D = RRTPI / (D0 + XQ*(D2 + XQ))
         K = D*Y * (A0 + XQ)
         L = d*x * (b1 + xq)         
      ELSEIF ( ABX .GT. XLIM2 ) THEN ! Humlicek W4 Region 2 
         IF ( RG2 .NE. 0 ) THEN ! First point in Region 2
            RG2 = 0
            H0 =  0.5625 + YQ*(4.5 + YQ*(10.5 + YQ*(6.0 + YQ))) ! Region 2 y-dependents
            H2 = -4.5    + YQ*(9.0 + YQ*( 6.0 + YQ* 4.0))
            H4 = 10.5    - YQ*(6.0 - YQ*  6.0)
            H6 = -6.0    + YQ* 4.0
            E0 =  1.875  + YQ*(8.25 + YQ*(5.5 + YQ))
            E2 =  5.25   + YQ*(1.0  + YQ* 3.0)
            E4 =  0.75*H6
            f1 = -1.875 + yq*(5.25 + yq*(4.5 + yq))
            f3 =  8.25  - yq*(1.0  - yq* 3.0)
            f5 = -5.5   + yq* 3.0
         ENDIF
         D = RRTPI / (H0 + XQ*(H2 + XQ*(H4 + XQ*(H6 + XQ))))
         K = D*Y   *(E0 + XQ*(E2 + XQ*(E4 + XQ)))
         L = d*x   *(f1 + xq*(f3 + xq*(f5 + xq)))         
      ELSEIF ( ABX .LT. XLIM3 ) THEN ! Humlicek W4 Region 3
         IF ( RG3 .NE. 0 ) THEN ! First point in Region 3
            RG3 = 0
            Z0 = 272.1014     + Y*(1280.829 + Y*(2802.870 + Y*(3764.966 ! Region 3 y-dependents
     &         + Y*(3447.629 + Y*(2256.981 + Y*(1074.409 + Y*(369.1989
     &         + Y*(88.26741 + Y*(13.39880 + Y)))))))))
            Z2 = 211.678      + Y*(902.3066 + Y*(1758.336 + Y*(2037.310
     &           + Y*(1549.675 + Y*(793.4273 + Y*(266.2987
     &           + Y*(53.59518 + Y*5.0)))))))
            Z4 = 78.86585     + Y*(308.1852 + Y*(497.3014 + Y*(479.2576
     &           + Y*(269.2916 + Y*(80.39278 + Y*10.0)))))
            Z6 = 22.03523     + Y*(55.02933 + Y*(92.75679 + Y*(53.59518
     &           + Y*10.0)))
            Z8 = 1.496460     + Y*(13.39880 + Y*5.0)
            P0 = 153.5168     + Y*(549.3954 + Y*(919.4955 + Y*(946.8970
     &          + Y*(662.8097 + Y*(328.2151 + Y*(115.3772 + Y*(27.93941
     &          + Y*(4.264678 + Y*0.3183291))))))))
            P2 = -34.16955    + Y*(-1.322256+ Y*(124.5975 + Y*(189.7730
     &           + Y*(139.4665 + Y*(56.81652 + Y*(12.79458
     &           + Y*1.2733163))))))
            P4 = 2.584042     + Y*(10.46332 + Y*(24.01655 + Y*(29.81482
     &           + Y*(12.79568 + Y*1.9099744))))
            P6= -0.07272979  + Y*(0.9377051+ Y*(4.266322 + Y*1.273316))
            P8 = 0.0005480304 + Y*0.3183291
            q1 = 173.2355  + y*(508.2585 + y*(685.8378 + y*(557.5178
     .                     + y*(301.3208 + y*(111.0528 + y*(27.62940
     .                     + y*(4.264130 + y*0.3183291)))))))
            q3 = 18.97431  + y*(100.7375 + y*(160.4013 + y*(130.8905
     .                     + y*(55.88650 + y*(12.79239+y*1.273316)))))
            q5 = 7.985877  + y*(19.83766 + y*(28.88480 + y*(12.79239
     .                     + y*1.909974)))
            q7 = 0.6276985 + y*(4.264130 + y*1.273316)            
         ENDIF
         D =1.7724538 / (Z0 + XQ*(Z2 + XQ*(Z4 + XQ*(Z6 + XQ*(Z8+XQ)))))
         K = D*(P0 + XQ*(P2 + XQ*(P4 + XQ*(P6 + XQ*P8))))
         L = d*x*(q1+xq*(q3+xq*(q5+xq*(q7 + xq*0.3183291))))         
      ELSE                      ! Humlicek CPF12 algorithm
         YPY0 = Y + Y0
         YPY0Q = YPY0*YPY0
         K = 0.0
         L = 0.0
         DO J = 0, 5
            D = X - T(J)
            MQ(J) = D*D
            MF(J) = 1.0 / (MQ(J) + YPY0Q)
            XM(J) = MF(J)*D
            YM(J) = MF(J)*YPY0
            D = X + T(J)
            PQ(J) = D*D
            PF(J) = 1.0 / (PQ(J) + YPY0Q)
            XP(J) = PF(J)*D
            YP(J) = PF(J)*YPY0
            L=L+c(j)*(xm(j)+xp(j)) + s(j)*(ym(j)-yp(j))
         ENDDO
         
         IF ( ABX .LE. XLIM4 ) THEN ! Humlicek CPF12 Region I
            DO J = 0, 5
               K = K + C(J)*(YM(J)+YP(J)) - S(J)*(XM(J)-XP(J))
            ENDDO            
         ELSE                   ! Humlicek CPF12 Region II
            YF   = Y + Y0PY0
            DO J = 0, 5
               K = K
     &     + (C(J)*(MQ(J)*MF(J)-Y0*YM(J)) + S(J)*YF*XM(J)) / (MQ(J)+Y0Q)
     &     + (C(J)*(PQ(J)*PF(J)-Y0*YP(J)) - S(J)*YF*XP(J)) / (PQ(J)+Y0Q)
            ENDDO
            K = Y*K + EXP ( -XQ )
         ENDIF
      ENDIF
*.....
      cres = CMPLX(K,L)
      RETURN
      END FUNCTION PH_HUMLIK


      end module photon
