      module photon
      USE PRECISION
      USE PARMMOD
      use COMXS

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
     .   PHV_N1STOTPH(:,:,:), PHV_N2NDOTPH(:,:,:), PHV_IS_PLSPHOT(:)

      real(dp), public, save :: CLIGHT,EV2HZ,HPLNK,hpcl
      real(dp), public, save :: hwvdw
    

      contains



      REAL(dp) FUNCTION PH_ENERGY(il,icell,kk,iipl,nladd) result(res)

      implicit none 
      integer, intent(in) :: il,icell,kk,iipl
      logical, intent(out) :: nladd
      
      nladd = .false.
      res = 0._dp
      write (6,*) ' PH_ENERGY: no calculations for photons possible '
      return
      end function ph_energy

      
      SUBROUTINE PH_GETCOEFF(kkin,isp,ity,icell,iipl,idsc,fac,res) 
      IMPLICIT NONE
      integer, intent(in) :: kkin,isp,ity,icell,iipl,idsc
      REAL(DP), INTENT(OUT) :: FAC,RES
      FAC = 0._DP
      res = 0._dp
      write (6,*) ' PH_GETCOEFF: no calculations for photons possible '
      return
      end SUBROUTINE ph_getcoeff


      real(dp) function ph_lorvdw(de,hw,shift,dvdw,
     .                            icell,ipl) result(res)
      implicit none
      real(dp), intent(in) :: de,hw,shift,dvdw
      integer, intent(in) :: icell,ipl
      res = 0._dp
      write (6,*) ' PH_LORVDW: no calculations for photons possible '
      return
      end function ph_lorvdw


      subroutine ph_plasma_profile(ind)
      implicit none
      integer, intent(in) :: ind
      write (6,*) ' PH_PLASMA_PROFILE: no calculations for photons',
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
         
         allocate(phv_is_plsphot(npls))
         phv_is_plsphot=0
      end if
      write (6,*) ' PH_INIT: no calculations for photons possible '

      return
      end subroutine PH_INIT


      subroutine ph_integrate(istr, scale, fpht,fatm,fmol,fion)
      implicit none
      real(dp), intent(in) :: scale,fpht,fatm,fmol,fion
      integer, intent(in) :: istr
      write (6,*) ' PH_INTEGRATE: no calculations for photons possible '
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

      allocate(PHV_LGAOT(0:natm,0:0,0:5))
      PHV_NROTA=0
      phv_lgaot=0         
      allocate(PHV_NAOTI(natm))
      phv_naoti=0
      
      allocate(phv_iestotat(0:natm,nnrot,3))
      phv_iestotat=0
      
      allocate(phv_n1stotat(0:natm,nnrot,3))
      allocate(phv_n2ndotat(0:natm,nnrot,3))
      phv_n1stotat=0
      phv_n2ndotat=0
      write (6,*) ' PH_ALLOC_XSECTA: no calculations for photons',
     .            ' possible '
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

      allocate(PHV_LGPHOT(0:nphot,0:0,0:5))
      phv_nrotph=0
      phv_lgphot=0
      allocate(phv_nphoti(nphot))
      phv_nphoti=0
      allocate(phv_iestotph(0:nphot,nnrot,3))
      phv_iestotph=0
      allocate(phv_n1stotph(0:nphot,nnrot,3))
      allocate(phv_n2ndotph(0:nphot,nnrot,3))
      phv_n1stotph=0
      phv_n2ndotph=0
      write (6,*) ' PH_ALLOC_XSECTPH: no calculations for photons',
     .            ' possible '

      return
      END SUBROUTINE PH_ALLOC_XSECTPH



      SUBROUTINE PH_XSECTPH(ipht,nrc,idsc)
      IMPLICIT NONE
      integer, intent(in) :: ipht,nrc,idsc
      write (6,*) ' PH_XSECTPH: no calculations for photons possible '
      return
      END SUBROUTINE PH_XSECTPH


      REAL(DP) FUNCTION PH_B21() result(res)
      IMPLICIT NONE
c calculates B21 Einstein coefficient in units: cm**2
      
      write (6,*) ' PH_B21: no calculations for photons possible '
      res=0.
      return
      END FUNCTION PH_B21


      end module photon
