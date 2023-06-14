module mod_sl_input


  implicit none

  ! Moving some Steve Lisgo common block contents into a module. 

      !COMMON /OSMPMKCOM/ save_osmpmk,scaleplateau
      REAL*8,allocatable ::     save_osmpmk(:,:)
      LOGICAL,allocatable ::    scaleplateau(:)

      !COMMON /OSMOPTS/ osm_fixopt
      INTEGER,allocatable ::    osm_fixopt(:,:)


      !COMMON /STATCOM/ osm_ncnt ,osm_nerr,osm_temod,osm_relfr,osm_cadj,osm_tcon
      INTEGER,allocatable ::    osm_ncnt (:)    ,osm_nerr (:,:),&
                                osm_cadj (:,:,:),osm_tcon (:)
      REAL,allocatable ::       osm_temod(:)    ,osm_relfr(:)

      !COMMON /CFSCOM/ cfs_mage,cfs_magi,cfs_sume,cfs_sumi
      REAL,allocatable ::        cfs_mage(:),cfs_magi(:),&
                                 cfs_sume(:),cfs_sumi(:)


      !COMMON /PININIT/ pininit
      LOGICAL,allocatable ::     pininit(:)

      !COMMON /ERRNOTE/ con_nerr,con_errc,con_errir
      !INTEGER          con_nerr
      !integer,allocatable :: con_errc(:),con_errir(:)


     !COMMON /PEIMULCOM/ osm_peimul2,osmpeires,
     !.                   restrictosmpmk,
     !.                   restrictmachno
      LOGICAL::            restrictosmpmk,restrictmachno
      REAL*8,allocatable ::            osm_peimul2(:,:),osmpeires(:,:)

      
contains



  subroutine allocate_mod_sl_input
    use mod_params
    use allocate_arrays
    implicit none

    integer :: ierr
    call allocate_array(save_osmpmk ,0,maxnks,1,maxnrs,'save_osmpmk',ierr)
    call allocate_array(scaleplateau,maxnrs,'scaleplateau',ierr)

    call allocate_array(osm_fixopt,2,maxnrs,'osm_fixopt',ierr)

    call allocate_array(osm_ncnt ,maxnrs,'osm_ncnt',ierr)
    call allocate_array(osm_nerr ,2,maxnrs,'osm_nerr',ierr)
    call allocate_array(osm_cadj ,5,2,maxnrs,'osm_cadj',ierr)
    call allocate_array(osm_tcon ,maxnrs,'osm_tcon',ierr)
    call allocate_array(osm_temod,maxnrs,'osm_temod',ierr)
    call allocate_array(osm_relfr,maxnrs,'osm_relfr',ierr)

    call allocate_array(cfs_mage,maxnrs,'cfs_mage',ierr)
    call allocate_array(cfs_magi,maxnrs,'cfs_magi',ierr)
    call allocate_array(cfs_sume,maxnrs,'cfs_sume',ierr)
    call allocate_array(cfs_sumi,maxnrs,'cfs_sumi',ierr)
    
    call allocate_array(pininit,maxnrs,'cfs_sumi',ierr)

    !call allocate_array(con_errc ,maxnrs,'con_errc',ierr)
    !call allocate_array(con_errir,maxnrs,'con_errir',ierr)

    call allocate_array(osm_peimul2,2,maxnrs,'osm_peimul2',ierr)
    call allocate_array(osmpeires ,2,maxnrs,'osm_peires',ierr)

    
  end subroutine allocate_mod_sl_input


  subroutine deallocate_mod_sl_input
    implicit none

    if (allocated(save_osmpmk )) deallocate (save_osmpmk )
    if (allocated(scaleplateau)) deallocate (scaleplateau)

    if (allocated(osm_fixopt)) deallocate (osm_fixopt)

    if (allocated(osm_ncnt )) deallocate (osm_ncnt )
    if (allocated(osm_nerr )) deallocate (osm_nerr )
    if (allocated(osm_cadj )) deallocate (osm_cadj )
    if (allocated(osm_tcon )) deallocate (osm_tcon )
    if (allocated(osm_temod)) deallocate (osm_temod)
    if (allocated(osm_relfr)) deallocate (osm_relfr)

    if (allocated(cfs_mage)) deallocate (cfs_mage)
    if (allocated(cfs_magi)) deallocate (cfs_magi)
    if (allocated(cfs_sume)) deallocate (cfs_sume)
    if (allocated(cfs_sumi)) deallocate (cfs_sumi)

    if (allocated(pininit)) deallocate (pininit)

    !if (allocated(con_errc )) deallocate (con_errc )
    !if (allocated(con_errir)) deallocate (con_errir)

    if (allocated(osm_peimul2)) deallocate (osm_peimul2)
    if (allocated(osmpeires )) deallocate (osmpeires )
    

  end subroutine deallocate_mod_sl_input



end module mod_sl_input
