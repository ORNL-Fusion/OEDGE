module mod_pindata

  use mod_params

  ! LIM does not currently support PIN/EIRENE ... however some of the plot code imported
  ! from DIVIMP contains these quantities so these are being added as potential future
  ! placeholders. 

  
!C                                                                       
!      COMMON /PINDATA/ PINATOM,PINION,PINALPHA,PINMOL
!c
!c     >       ,pinz0,pinionz,pinena,pinenm,pinenz,                        
!c     >       pinqi,pinqe,pinmp,pinrec,divrec,
!c     >       hioniz,zioniz,h2ioniz,                 
!c
!c     >       srecyc,srecom,hescpd,hescal,zsput,zsputn,
!c     >       zescpd,zescal,hesclk,zesclk,
!c     >       hescpd_last,hescal_last,phfgal_last,phfuga_last,
!c     >       pinvdist,hcorr,hval,hwalks,
!c
!c            Original PIN wall quantities
!c
!c     >       rvesm_pin,zvesm_pin,jvesm_pin,
!c     >       fluxhw_pin,flxhw2_pin,
!c     >       flxhw3_pin,flxhw4_pin,flxhw5_pin,
!c     >       flxhw6_pin,nvesm_pin,
!c
!c            Copies of PIN wall quantities that may be 
!c            modified to match the DIVIMP wall - if baffles
!c            have been included in a redefined wall.
!c
!c     >       rvesm,zvesm,jvesm,fluxhw,flxhw2,
!c slmod begin - new
!c...         Added FLXHW7, the average molecular hydrogen energy:
!c... jdemod  added flxhw8, the Eirene reported wall ion flux
!c     >       flxhw3,flxhw4,flxhw5,flxhw6,flxhw7,flxhw8,nvesm,
!c
!c     >       flxhw3,flxhw4,flxhw5,flxhw6,nvesm,
!c slmod end
!c
!c     >       nvesp,
!c     >       pincor,cnimbin,nlines,
!c     >       pinpuff,hpcpuf,phfgal,swpvhpf,jhpuf1,jhpuf2,
!c     >       tpufh,ihybrid,phfuga,ppcpuf,acthpcpuf,
!c     >       gaugedat,ihcorr,iiterpin,piniz_info,pinprint,
!c     >       piniseed,pinseed
!c
!c    >       ,hextrl,phxtra
!c
!c      integer maxlines
!c      parameter (maxlines=50)
!C                                                                       
!      REAL PINATOM(MAXNXS,MAXNYS),PINION(MAXNXS,MAXNYS),                
!     >     PINALPHA(MAXNXS,MAXNYS),PINMOL(MAXNXS,MAXNYS)
!c
!c     >     ,pinz0(maxnxs,maxnys),pinionz(maxnxs,maxnys),                 
!c     >     pinena(maxnxs,maxnys),pinenm(maxnxs,maxnys),                 
!c     >     pinenz(maxnxs,maxnys),pinqi(maxnxs,maxnys),                  
!c     >     pinqe(maxnxs,maxnys),pinmp(maxnxs,maxnys),
!c     >     pinrec(maxnxs,maxnys),divrec(maxnxs,maxnys),
!c     >     hioniz,zioniz,h2ioniz,
!c     >     srecyc,srecom,hescpd,hescal,zsput,zsputn,
!c     >     zescpd,zescal,hesclk,zesclk,
!c     >     hescpd_last,hescal_last,phfgal_last,phfuga_last,
!c     >     pinvdist(3,14,maxnxs,maxnys),hwalks(maxnws,2),pincor,
!c
!c     >     rvesm_pin(maxseg,2),zvesm_pin(maxseg,2),
!c     >     fluxhw_pin(maxseg),flxhw2_pin(maxseg),
!c     >     flxhw3_pin(maxseg),flxhw4_pin(maxseg),
!c     >     flxhw5_pin(maxseg),flxhw6_pin(maxseg),
!c
!c     >     rvesm(maxseg,2),zvesm(maxseg,2),
!c     >     fluxhw(maxseg),flxhw2(maxseg),
!c     >     flxhw3(maxseg),flxhw4(maxseg),
!c slmod begin - new
!c...       Added FLXHW7, the average molecular hydrogen energy:
!c... jdemod  added flxhw8, the Eirene reported wall ion flux
!c     >     flxhw5(maxseg),flxhw6(maxseg),flxhw7(maxseg),
!c     >     flxhw8(maxseg),
!c
!c     >     flxhw5(maxseg),flxhw6(maxseg),
!c slmod end
!c
!c     >     hpcpuf,phfgal,tpufh,phfuga,ppcpuf,acthpcpuf,
!c     >     gaugedat(5,2),piniz_info(maxnys,4),
!c     >     hcorr(maxnxs,maxnys),hval(maxnxs,maxnys)
!c
!c    >     ,hextrl,phxtra
!C                                                                       
!c      integer nvesm,nvesp,jvesm(maxseg),nlines,pinpuff,
!c     >        swpvhpf,jhpuf1(2),jhpuf2(2),ihybrid,ihcorr,
!c     >        iiterpin,pinprint,piniseed,pinseed,
!c     >        nvesm_pin,jvesm_pin(maxseg)
!c
!c
!c      character*80 cnimbin(maxlines)
!c


  implicit none
  private


      REAL,allocatable,public:: PINATOM(:,:),PINION(:,:),&
                                PINALPHA(:,:),PINMOL(:,:)

  
  public :: allocate_mod_pindata, deallocate_mod_pindata


contains

  subroutine allocate_mod_pindata
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(DTEV  ,maxnxs,'DTEV',ierr)
     call allocate_array(PINATOM ,MAXNXS,MAXNYS,'PINATOM ',ierr)
     call allocate_array(PINION  ,MAXNXS,MAXNYS,'PINION  ',ierr)
     call allocate_array(PINALPHA,MAXNXS,MAXNYS,'PINALPHA',ierr)
     call allocate_array(PINMOL  ,MAXNXS,MAXNYS,'PINMOL  ',ierr)


  end subroutine allocate_mod_pindata


  subroutine deallocate_mod_pindata
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()
    deallocate(PINATOM )
    deallocate(PINION  )
    deallocate(PINALPHA)
    deallocate(PINMOL  )

  end subroutine deallocate_mod_pindata



end module mod_pindata
