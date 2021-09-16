C  27.6.05: phv_nrota, phv_nrotph removed
      module photon
      USE PRECISION
      USE PARMMOD
      use COMXS
      use comprt
      use ccona

      implicit none

      private
      
      public :: ph_init, ph_plasma_profile, ph_energy, ph_getcoeff,
     .    ph_xsectph,
     .    ph_post_energy, ph_alloc_xsectph,
     .    ph_lorvdw, ph_vdwprof, ph_b21, line_cutoff

      integer, public, save :: actual_linenum
      integer, public, save  :: phv_muldens
      integer, public, save, allocatable ::
     .   PHV_LGAOT(:,:,:), PHV_LGPHOT(:,:,:),
     .   PHV_NAOTI(:), PHV_NPHOTI(:),
     .   PHV_N1STOTAT(:,:,:), PHV_N2NDOTAT(:,:,:),
     .   PHV_IESTOTAT(:,:,:), PHV_IESTOTPH(:,:,:),
     .   PHV_N1STOTPH(:,:,:), PHV_N2NDOTPH(:,:,:), PHV_IS_PLSPHOT(:)

!pb black body removal (core saturation)
      real(dp), allocatable, public, save  :: 
     .          xintleft(:,:), xintright(:,:), xint_inf(:,:)

      real(dp), public, save :: hwvdw


      contains



      REAL(dp) FUNCTION PH_ENERGY(icell,kk,iipl,vn,nladd,
     .                            e_min,e_max) result(res)
      implicit none
      integer, intent(in) :: icell,kk,iipl
      logical, intent(out) :: nladd
      real(dp), intent(in) :: vn
      real(dp) :: e_min,e_max
  
      nladd = .false.
      res = 0._dp
      write (iunout,*) 
     .  ' PH_ENERGY: no calculations for photons possible '
      return
      end function ph_energy


      SUBROUTINE PH_GETCOEFF(kkin,isp,ity,icell,iipl,fac,res)
      IMPLICIT NONE
      integer, intent(in) :: kkin,isp,ity,icell,iipl
      REAL(DP), INTENT(OUT) :: FAC,RES
      FAC = 0._DP
      res = 0._dp
      write (iunout,*) 
     .  ' PH_GETCOEFF: no calculations for photons possible '
      return
      end SUBROUTINE ph_getcoeff


      real(dp) function ph_lorvdw(de,hw,shift,dvdw,
     .                            icell,ipl) result(res)
      implicit none
      real(dp), intent(in) :: de,hw,shift,dvdw
      integer, intent(in) :: icell,ipl
      res = 0._dp
      write (iunout,*) 
     .  ' PH_LORVDW: no calculations for photons possible '
      return
      end function ph_lorvdw


      subroutine ph_plasma_profile(ind)
      implicit none
      integer, intent(in) :: ind
      write (iunout,*) 
     .   ' PH_PLASMA_PROFILE: no calculations for photons',
     .   ' possible '
      return
      end subroutine ph_plasma_profile


      SUBROUTINE PH_INIT(ICAL)
      implicit none
      integer, intent(in) :: ical

      if (ical == 1) then

         allocate(phv_is_plsphot(npls))
         phv_is_plsphot=0
      end if
      write (iunout,*) ' PH_INIT: no calculations for photons possible '

      return
      end subroutine PH_INIT


      SUBROUTINE PH_POST_ENERGY(icell,kk,iflg,il,
     .                          iold,itypold,vxo,vyo,vzo,vlo,e0o,itypnw)
      IMPLICIT NONE
      integer, intent(in) :: icell,kk,iflg,iold,il,itypold,itypnw
      real(dp),intent(in) :: vxo,vyo,vzo,vlo,e0o
      write (iunout,*) ' PH_POST_ENERGY: no calculations for photons',
     .            ' possible '
      return
      end subroutine PH_POST_ENERGY


      SUBROUTINE PH_VDWPROF(icell,hw,shift,dvdw,lscale)
      implicit none
      integer, intent(in) :: icell
      real(dp),intent(out) :: hw,shift,dvdw
      logical, intent(in) :: lscale
      write (iunout,*) 
     .  ' PH_VDWPROF: no calculations for photons possible '
      hw=0._dp
      shift=0._dp
      dvdw=0._dp
      return
      end subroutine PH_VDWPROF


      SUBROUTINE PH_ALLOC_XSECTPH(nnrot)
      IMPLICIT NONE
      integer, intent(in) :: nnrot

      allocate(PHV_LGPHOT(0:nphot,0:0,0:5))
      phv_lgphot=0
      allocate(phv_nphoti(nphot))
      phv_nphoti=0
      allocate(phv_iestotph(0:nphot,nnrot,3))
      phv_iestotph=0
      allocate(phv_n1stotph(0:nphot,nnrot,3))
      allocate(phv_n2ndotph(0:nphot,nnrot,3))
      phv_n1stotph=0
      phv_n2ndotph=0
      write (iunout,*) ' PH_ALLOC_XSECTPH: no calculations for photons',
     .            ' possible '

      return
      END SUBROUTINE PH_ALLOC_XSECTPH



      SUBROUTINE PH_XSECTPH(ipht,nrc,idsc)
      IMPLICIT NONE
      integer, intent(in) :: ipht,nrc,idsc
      write (iunout,*) 
     .  ' PH_XSECTPH: no calculations for photons possible '
      return
      END SUBROUTINE PH_XSECTPH


      REAL(DP) FUNCTION PH_B21() result(res)
      IMPLICIT NONE
c calculates B21 Einstein coefficient in units: cm**2

      write (iunout,*) ' PH_B21: no calculations for photons possible '
      res=0.
      return
      END FUNCTION PH_B21




      SUBROUTINE line_cutoff
      IMPLICIT NONE
      write (iunout,*) 
     .  ' line_cutoff: no calculations for photons possible '
      return
      END SUBROUTINE line_cutoff

      end module photon
