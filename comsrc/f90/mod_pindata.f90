module mod_pindata
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  !
  !            original pin wall quantities
  !
  !
  !            copies of pin wall quantities that may be
  !            modified to match the divimp wall - if baffles
  !            have been included in a redefined wall.
  !
  ! slmod begin - new
  !...         added flxhw7, the average molecular hydrogen energy:
  !... jdemod  added flxhw8, the eirene reported wall ion flux
  !
  !     >       flxhw3,flxhw4,flxhw5,flxhw6,nvesm,
  ! slmod end
  !
  !
  !    >       ,hextrl,phxtra
  !
  ! common /pindata/ pinatom,pinion,pinalpha,pinmol,pinz0,pinionz,pinena,pinenm,pinenz,&
  !     pinqi,pinqe,pinmp,pinrec,divrec,hioniz,zioniz,h2ioniz,srecyc,srecom,hescpd,hescal,&
  !     zsput,zsputn,zescpd,zescal,hesclk,zesclk,hescpd_last,hescal_last,phfgal_last,&
  !     phfuga_last,pinvdist,hcorr,hval,hwalks,rvesm_pin,zvesm_pin,jvesm_pin,fluxhw_pin,&
  !     flxhw2_pin,flxhw3_pin,flxhw4_pin,flxhw5_pin,flxhw6_pin,nvesm_pin,rvesm,zvesm,&
  !     jvesm,fluxhw,flxhw2,flxhw3,flxhw4,flxhw5,flxhw6,flxhw7,flxhw8,nvesm,nvesp,pincor,&
  !     cnimbin,nlines,pinpuff,hpcpuf,phfgal,swpvhpf,jhpuf1,jhpuf2,tpufh,ihybrid,phfuga,&
  !     ppcpuf,acthpcpuf,gaugedat,ihcorr,iiterpin,piniz_info,pinprint,piniseed,pinseed
  !
  ! save /pindata/
  integer,public :: maxlines
  !
  ! slmod begin
  parameter (maxlines=50)
  !
  !      real pinatom(maxnks,maxnrs),pinion(maxnks,maxnrs),
  ! slmod end
  !
  !
  ! slmod begin - new
  !...       added flxhw7, the average molecular hydrogen energy:
  !... jdemod  added flxhw8, the eirene reported wall ion flux
  !
  !     >     flxhw5(maxseg),flxhw6(maxseg),
  ! slmod end
  !
  !
  !    >     ,hextrl,phxtra
  !
  real, target ,public :: hioniz,zioniz,h2ioniz,srecyc,srecom,hescpd,hescal,zsput,zsputn,&
       zescpd,zescal,hesclk,zesclk,hescpd_last,hescal_last,phfgal_last,phfuga_last,&
       pincor,hpcpuf,phfgal,tpufh,phfuga,ppcpuf,acthpcpuf
  real, target ,public,allocatable :: pinatom(:,:),pinion(:,:),pinalpha(:,:),pinmol(:,:),&
       pinz0(:,:),pinionz(:,:),pinena(:,:),pinenm(:,:),pinenz(:,:),pinqi(:,:),pinqe(:,:),&
       pinmp(:,:),pinrec(:,:),divrec(:,:),pinvdist(:,:,:,:),hwalks(:,:),rvesm_pin(:,:),&
       zvesm_pin(:,:),fluxhw_pin(:),flxhw2_pin(:),flxhw3_pin(:),flxhw4_pin(:),flxhw5_pin(:),&
       flxhw6_pin(:),rvesm(:,:),zvesm(:,:),fluxhw(:),flxhw2(:),flxhw3(:),flxhw4(:),&
       flxhw5(:),flxhw6(:),flxhw7(:),flxhw8(:),gaugedat(:,:),piniz_info(:,:),hcorr(:,:),&
       hval(:,:)

  ! jdemod - sol22 external power terms - used in the same place as PIN terms
  real,allocatable,public :: ext_epowsrc(:,:),ext_ipowsrc(:,:)              ! external electron and ion power terms
  real,allocatable,public :: div_tpowls(:,:),div_tcooliz(:,:),div_cool(:,:) ! impurity power terms
  
  !
  integer,public :: nvesm,nvesp,nlines,pinpuff,swpvhpf,ihybrid,ihcorr,iiterpin,pinprint,&
       piniseed,pinseed,nvesm_pin
  integer,public,allocatable :: jvesm(:),jhpuf1(:),jhpuf2(:),jvesm_pin(:)
  
  !
  character*80,public :: cnimbin(maxlines)

  public :: allocate_mod_pindata,deallocate_mod_pindata

contains

  subroutine allocate_mod_pindata
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_pindata','ALLOCATE')

    call allocate_array(pinatom,maxnks,maxnrs,'pinatom',ierr)
    call allocate_array(pinion,maxnks,maxnrs,'pinion',ierr)
    call allocate_array(pinalpha,maxnks,maxnrs,'pinalpha',ierr)
    call allocate_array(pinmol,maxnks,maxnrs,'pinmol',ierr)
    call allocate_array(pinz0,maxnks,maxnrs,'pinz0',ierr)
    call allocate_array(pinionz,maxnks,maxnrs,'pinionz',ierr)
    call allocate_array(pinena,maxnks,maxnrs,'pinena',ierr)
    call allocate_array(pinenm,maxnks,maxnrs,'pinenm',ierr)
    call allocate_array(pinenz,maxnks,maxnrs,'pinenz',ierr)
    call allocate_array(pinqi,maxnks,maxnrs,'pinqi',ierr)
    call allocate_array(pinqe,maxnks,maxnrs,'pinqe',ierr)
    call allocate_array(pinmp,maxnks,maxnrs,'pinmp',ierr)
    call allocate_array(pinrec,maxnks,maxnrs,'pinrec',ierr)
    call allocate_array(divrec,maxnks,maxnrs,'divrec',ierr)
    call allocate_array(pinvdist,1,3,1,14,1,maxnks,1,maxnrs,'pinvdist',ierr)
    !
    ! jdemod - add externally loaded quantities used in the plasma solver
    !          NOT necessarily sourced from PIN (EIRENE) but used in the
    !          same routines
    call allocate_array(ext_epowsrc,maxnks,maxnrs,'epowsrc',ierr)
    call allocate_array(ext_ipowsrc,maxnks,maxnrs,'ipowsrc',ierr)
    call allocate_array(div_tpowls,maxnks,maxnrs,'div_tpowls',ierr)
    call allocate_array(div_tcooliz,maxnks,maxnrs,'div_tcooliz',ierr)
    call allocate_array(div_cool,maxnks,maxnrs,'div_cool',ierr)
    !
    call allocate_array(hwalks,maxnws,2,'hwalks',ierr)
    call allocate_array(rvesm_pin,maxseg,2,'rvesm_pin',ierr)
    call allocate_array(zvesm_pin,maxseg,2,'zvesm_pin',ierr)
    call allocate_array(fluxhw_pin,maxseg,'fluxhw_pin',ierr)
    call allocate_array(flxhw2_pin,maxseg,'flxhw2_pin',ierr)
    call allocate_array(flxhw3_pin,maxseg,'flxhw3_pin',ierr)
    call allocate_array(flxhw4_pin,maxseg,'flxhw4_pin',ierr)
    call allocate_array(flxhw5_pin,maxseg,'flxhw5_pin',ierr)
    call allocate_array(flxhw6_pin,maxseg,'flxhw6_pin',ierr)
    call allocate_array(rvesm,maxseg,2,'rvesm',ierr)
    call allocate_array(zvesm,maxseg,2,'zvesm',ierr)
    call allocate_array(fluxhw,maxseg,'fluxhw',ierr)
    call allocate_array(flxhw2,maxseg,'flxhw2',ierr)
    call allocate_array(flxhw3,maxseg,'flxhw3',ierr)
    call allocate_array(flxhw4,maxseg,'flxhw4',ierr)
    call allocate_array(flxhw5,maxseg,'flxhw5',ierr)
    call allocate_array(flxhw6,maxseg,'flxhw6',ierr)
    call allocate_array(flxhw7,maxseg,'flxhw7',ierr)
    call allocate_array(flxhw8,maxseg,'flxhw8',ierr)
    call allocate_array(gaugedat,5,2,'gaugedat',ierr)
    call allocate_array(piniz_info,maxnrs,4,'piniz_info',ierr)
    call allocate_array(hcorr,maxnks,maxnrs,'hcorr',ierr)
    call allocate_array(hval,maxnks,maxnrs,'hval',ierr)
    call allocate_array(jvesm,maxseg,'jvesm',ierr)
    call allocate_array(jhpuf1,2,'jhpuf1',ierr)
    call allocate_array(jhpuf2,2,'jhpuf2',ierr)
    call allocate_array(jvesm_pin,maxseg,'jvesm_pin',ierr)

  end subroutine allocate_mod_pindata


  subroutine deallocate_mod_pindata
    implicit none

    call pr_trace('mod_pindata','DEALLOCATE')

    if (allocated(pinatom)) deallocate(pinatom)
    if (allocated(pinion)) deallocate(pinion)
    if (allocated(pinalpha)) deallocate(pinalpha)
    if (allocated(pinmol)) deallocate(pinmol)
    if (allocated(pinz0)) deallocate(pinz0)
    if (allocated(pinionz)) deallocate(pinionz)
    if (allocated(pinena)) deallocate(pinena)
    if (allocated(pinenm)) deallocate(pinenm)
    if (allocated(pinenz)) deallocate(pinenz)
    if (allocated(pinqi)) deallocate(pinqi)
    if (allocated(pinqe)) deallocate(pinqe)
    if (allocated(pinmp)) deallocate(pinmp)
    if (allocated(pinrec)) deallocate(pinrec)
    if (allocated(divrec)) deallocate(divrec)
    if (allocated(pinvdist)) deallocate(pinvdist)
    if (allocated(hwalks)) deallocate(hwalks)
    if (allocated(rvesm_pin)) deallocate(rvesm_pin)
    if (allocated(zvesm_pin)) deallocate(zvesm_pin)
    if (allocated(fluxhw_pin)) deallocate(fluxhw_pin)
    if (allocated(flxhw2_pin)) deallocate(flxhw2_pin)
    if (allocated(flxhw3_pin)) deallocate(flxhw3_pin)
    if (allocated(flxhw4_pin)) deallocate(flxhw4_pin)
    if (allocated(flxhw5_pin)) deallocate(flxhw5_pin)
    if (allocated(flxhw6_pin)) deallocate(flxhw6_pin)
    if (allocated(rvesm)) deallocate(rvesm)
    if (allocated(zvesm)) deallocate(zvesm)
    if (allocated(fluxhw)) deallocate(fluxhw)
    if (allocated(flxhw2)) deallocate(flxhw2)
    if (allocated(flxhw3)) deallocate(flxhw3)
    if (allocated(flxhw4)) deallocate(flxhw4)
    if (allocated(flxhw5)) deallocate(flxhw5)
    if (allocated(flxhw6)) deallocate(flxhw6)
    if (allocated(flxhw7)) deallocate(flxhw7)
    if (allocated(flxhw8)) deallocate(flxhw8)
    if (allocated(gaugedat)) deallocate(gaugedat)
    if (allocated(piniz_info)) deallocate(piniz_info)
    if (allocated(hcorr)) deallocate(hcorr)
    if (allocated(hval)) deallocate(hval)
    if (allocated(jvesm)) deallocate(jvesm)
    if (allocated(jhpuf1)) deallocate(jhpuf1)
    if (allocated(jhpuf2)) deallocate(jhpuf2)
    if (allocated(jvesm_pin)) deallocate(jvesm_pin)

  end subroutine deallocate_mod_pindata

end module mod_pindata
