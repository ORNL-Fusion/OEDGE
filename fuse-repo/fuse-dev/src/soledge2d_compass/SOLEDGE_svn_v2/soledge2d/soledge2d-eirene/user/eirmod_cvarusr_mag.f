      module eirmod_cvarusr_mag


      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CGRID

      IMPLICIT NONE

      private

      public :: alloc_cvarusr
      public :: dealloc_cvarusr_mag

      real(dp), public, save :: radciel(200), drad(200)

      real(dp), public, allocatable, save :: flxciel(:,:), delph(:),
     .                                       area_ciel(:,:),area1d(:)
      real(dp), public, allocatable, save :: tciel(:,:),flxYchem(:,:)

      real(dp), public, allocatable, save :: vfunc(:)
      integer, public, allocatable, save :: ifunc(:)

      integer, public, save :: nrciel=0, nflxcl=0, nbincll=0

      real(dp), public, allocatable, save :: 
     .  pfl_in_at(:,:),pfl_out_at(:,:),enfl_in_at(:,:),enfl_out_at(:,:),
     .  pfl_in_ml(:,:),pfl_out_ml(:,:),enfl_in_ml(:,:),enfl_out_ml(:,:),
     .  pfl_in_io(:,:),pfl_out_io(:,:),enfl_in_io(:,:),enfl_out_io(:,:),
     .  pfl_in_pl(:,:),pfl_out_pl(:,:),enfl_in_pl(:,:),enfl_out_pl(:,:)

      real(dp), public, allocatable, save :: 
     .  tpfl_in_at(:),tpfl_out_at(:),tenfl_in_at(:),tenfl_out_at(:),
     .  tpfl_in_ml(:),tpfl_out_ml(:),tenfl_in_ml(:),tenfl_out_ml(:),
     .  tpfl_in_io(:),tpfl_out_io(:),tenfl_in_io(:),tenfl_out_io(:),
     .  tpfl_in_pl(:),tpfl_out_pl(:),tenfl_in_pl(:),tenfl_out_pl(:)

      real(dp), public, allocatable, save :: 
     .  cpfl_in_at(:,:),cpfl_out_at(:,:),
     .  cenfl_in_at(:,:),cenfl_out_at(:,:),
     .  cpfl_in_ml(:,:),cpfl_out_ml(:,:),
     .  cenfl_in_ml(:,:),cenfl_out_ml(:,:),
     .  cpfl_in_io(:,:),cpfl_out_io(:,:),
     .  cenfl_in_io(:,:),cenfl_out_io(:,:),
     .  cpfl_in_pl(:,:),cpfl_out_pl(:,:),
     .  cenfl_in_pl(:,:),cenfl_out_pl(:,:)

c     magnetic field data

      integer, public, allocatable, save :: celnumB(:)
      real(dp), public, allocatable, save :: BTSx(:),BTSy(:),
     .   BTSz(:),BTSf(:)
      integer, public ,save :: icel

      contains

      subroutine alloc_cvarusr (ical)

      implicit none
      integer, intent(in) :: ical

      if (ical == 1) then

        allocate (pfl_in_at(0:natm,nbincll))
        allocate (pfl_out_at(0:natm,nbincll))
        allocate (enfl_in_at(0:natm,nbincll))
        allocate (enfl_out_at(0:natm,nbincll))
        allocate (pfl_in_ml(0:nmol,nbincll))
        allocate (pfl_out_ml(0:nmol,nbincll))
        allocate (enfl_in_ml(0:nmol,nbincll))
        allocate (enfl_out_ml(0:nmol,nbincll))
        allocate (pfl_in_io(0:nion,nbincll))
        allocate (pfl_out_io(0:nion,nbincll))
        allocate (enfl_in_io(0:nion,nbincll))
        allocate (enfl_out_io(0:nion,nbincll))
        allocate (pfl_in_pl(0:npls,nbincll))
        allocate (pfl_out_pl(0:npls,nbincll))
        allocate (enfl_in_pl(0:npls,nbincll))
        allocate (enfl_out_pl(0:npls,nbincll))

        pfl_in_at = 0._dp
        pfl_out_at = 0._dp
        enfl_in_at = 0._dp
        enfl_out_at = 0._dp
        pfl_in_ml = 0._dp
        pfl_out_ml = 0._dp
        enfl_in_ml = 0._dp
        enfl_out_ml = 0._dp
        pfl_in_io = 0._dp
        pfl_out_io = 0._dp
        enfl_in_io = 0._dp
        enfl_out_io = 0._dp
        pfl_in_pl = 0._dp
        pfl_out_pl = 0._dp
        enfl_in_pl = 0._dp
        enfl_out_pl = 0._dp

        allocate (tpfl_in_at(0:natm))
        allocate (tpfl_out_at(0:natm))
        allocate (tenfl_in_at(0:natm))
        allocate (tenfl_out_at(0:natm))
        allocate (tpfl_in_ml(0:nmol))
        allocate (tpfl_out_ml(0:nmol))
        allocate (tenfl_in_ml(0:nmol))
        allocate (tenfl_out_ml(0:nmol))
        allocate (tpfl_in_io(0:nion))
        allocate (tpfl_out_io(0:nion))
        allocate (tenfl_in_io(0:nion))
        allocate (tenfl_out_io(0:nion))
        allocate (tpfl_in_pl(0:npls))
        allocate (tpfl_out_pl(0:npls))
        allocate (tenfl_in_pl(0:npls))
        allocate (tenfl_out_pl(0:npls))

        tpfl_in_at = 0._dp
        tpfl_out_at = 0._dp
        tenfl_in_at = 0._dp
        tenfl_out_at = 0._dp
        tpfl_in_ml = 0._dp
        tpfl_out_ml = 0._dp
        tenfl_in_ml = 0._dp
        tenfl_out_ml = 0._dp
        tpfl_in_io = 0._dp
        tpfl_out_io = 0._dp
        tenfl_in_io = 0._dp
        tenfl_out_io = 0._dp
        tpfl_in_pl = 0._dp
        tpfl_out_pl = 0._dp
        tenfl_in_pl = 0._dp
        tenfl_out_pl = 0._dp

      else if (ical == 2) then

        allocate (cpfl_in_at(0:natm,nbincll))
        allocate (cpfl_out_at(0:natm,nbincll))
        allocate (cenfl_in_at(0:natm,nbincll))
        allocate (cenfl_out_at(0:natm,nbincll))
        allocate (cpfl_in_ml(0:nmol,nbincll))
        allocate (cpfl_out_ml(0:nmol,nbincll))
        allocate (cenfl_in_ml(0:nmol,nbincll))
        allocate (cenfl_out_ml(0:nmol,nbincll))
        allocate (cpfl_in_io(0:nion,nbincll))
        allocate (cpfl_out_io(0:nion,nbincll))
        allocate (cenfl_in_io(0:nion,nbincll))
        allocate (cenfl_out_io(0:nion,nbincll))
        allocate (cpfl_in_pl(0:npls,nbincll))
        allocate (cpfl_out_pl(0:npls,nbincll))
        allocate (cenfl_in_pl(0:npls,nbincll))
        allocate (cenfl_out_pl(0:npls,nbincll))

        cpfl_in_at = 0._dp
        cpfl_out_at = 0._dp
        cenfl_in_at = 0._dp
        cenfl_out_at = 0._dp
        cpfl_in_ml = 0._dp
        cpfl_out_ml = 0._dp
        cenfl_in_ml = 0._dp
        cenfl_out_ml = 0._dp
        cpfl_in_io = 0._dp
        cpfl_out_io = 0._dp
        cenfl_in_io = 0._dp
        cenfl_out_io = 0._dp
        cpfl_in_pl = 0._dp
        cpfl_out_pl = 0._dp
        cenfl_in_pl = 0._dp
        cenfl_out_pl = 0._dp

      else if (ical == 3) then
      
      	allocate (celnumB(NSBOX))
        allocate (BTSx(NSBOX),BTSy(NSBOX),BTSz(NSBOX),BTSf(NSBOX))

        celnumB=0
      	BTSx=0._dp
        BTSy=0._dp
        BTSz=0._dp
        BTSf=0._dp

      end if

      return
      end subroutine alloc_cvarusr

      subroutine dealloc_cvarusr_mag
      implicit none

      deallocate(celnumB)
      deallocate(BTSx,BTSy,BTSz,BTSf)
      
      end subroutine dealloc_cvarusr_mag


      end module EIRMOD_cvarusr_mag
