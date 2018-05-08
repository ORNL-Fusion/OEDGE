module mod_dynam3


  !                                                                               
  !      COMMON /DYNAM3/ TIZS,TIZ3,POWLS,LINES,PLRPS,ZEFFS,           
  !     >  PLRP3,LIM5,SAVES,SDTXS,SDTYS,SDYXS,SDYYS,
  !     >  SCTXS,SCTYS,
  !     >  POWL3,LINE3
  !      REAL TIZS  (MAXNXS,-MAXNYS:MAXNYS,-1:MAXIZS)                              
  !      REAL TIZ3  (MAXNXS,-MAXY3D:MAXY3D,-1:MAXIZS,-MAXNPS:MAXNPS)               
  !      REAL ZEFFS (MAXNXS,-MAXNYS:MAXNYS,6),  SAVES(MAXNXS,-2:MAXIZS+1)          
  !      REAL POWLS (MAXNXS,-MAXNYS:MAXNYS,-1:MAXIZS)                              
  !      REAL LINES (MAXNXS,-MAXNYS:MAXNYS,-1:MAXIZS)                              
  !      REAL PLRPS (MAXNXS,-MAXNYS:MAXNYS,   MAXNLS)                              
  !      REAL POWL3 (MAXNXS,-MAXY3D:MAXY3D,-1:MAXIZS,-MAXNPS:MAXNPS)               
  !      REAL LINE3 (MAXNXS,-MAXY3D:MAXY3D,-1:MAXIZS,-MAXNPS:MAXNPS)               
  !      REAL PLRP3 (MAXNXS,-MAXY3D:MAXY3D,   MAXNLS,-MAXNPS:MAXNPS)               
  !      REAL LIM5  (MAXNXS,-MAXY3D:MAXY3D,-1:MAXIZS,-MAXNPS:MAXNPS,MAXNTS)        
  !      REAL SDTXS (MAXNXS,MAXIZS),SDTYS(-MAXNYS:MAXNYS,MAXIZS)                   
  !      REAL SDYXS (MAXNXS,MAXIZS),SDYYS(-MAXNYS:MAXNYS,MAXIZS)                   
  !      REAL SCTXS (MAXNXS,MAXNLS),SCTYS(-MAXNYS:MAXNYS,MAXNLS)
  !



  implicit none

  private

  REAL,allocatable,public:: TIZS(:,:,:),TIZ3(:,:,:,:),ZEFFS (:,:,:),&
       SAVES(:,:),POWLS (:,:,:), LINES(:,:,:), PLRPS (:,:,:),&
       POWL3 (:,:,:,:),LINE3 (:,:,:,:),PLRP3 (:,:,:,:),&               
       LIM5  (:,:,:,:,:),SDTXS (:,:),SDTYS(:,:),SDYXS (:,:),SDYYS(:,:),SCTXS (:,:),SCTYS(:,:)

  public:: allocate_mod_dynam3,deallocate_mod_dynam3



contains


  subroutine allocate_mod_dynam3
    use mod_params
    use allocate_arrays

    implicit none
    integer :: ierr


    call allocate_array(TIZS ,1,maxnxs,-maxnys,maxnys,-1,maxizs,'TIZS ',ierr)
    call allocate_array(ZEFFS,1,maxnxs,-maxnys,maxnys,1,6,'ZEFFS',ierr)
    call allocate_array(SAVES,1,maxnxs,-2,maxizs+1,'SAVES',ierr)

    call allocate_array(TIZ3 ,1,maxnxs,-maxnys,maxnys,-1,maxizs,-maxnps,maxnps,'TIZ3 ',ierr)
    call allocate_array(POWLS,1,maxnxs,-maxnys,maxnys,-1,maxizs,'POWLS',ierr)
    call allocate_array(LINES,1,maxnxs,-maxnys,maxnys,-1,maxizs,'LINES',ierr)
    call allocate_array(PLRPS,1,maxnxs,-maxnys,maxnys, 1,maxnls,'PLRPS',ierr)

    call allocate_array(POWL3,1,maxnxs,-maxy3d,maxy3d,-1,maxizs,-maxnps,maxnps,'POWL3',ierr)
    call allocate_array(LINE3,1,maxnxs,-maxy3d,maxy3d,-1,maxizs,-maxnps,maxnps,'LINE3',ierr)
    call allocate_array(PLRP3,1,maxnxs,-maxy3d,maxy3d, 1,maxnls,-maxnps,maxnps,'PLRP3',ierr)

    call allocate_array(LIM5 ,1,maxnxs,-maxy3d,maxy3d,-1,maxizs,-maxnps,maxnps,1,maxnts,'LIM5 ',ierr)

    call allocate_array(SDTXS,1,maxnxs,1,maxizs,'SDTXS',ierr)
    call allocate_array(SDTYS,-maxnys,maxnys,1,maxizs,'SDTYS',ierr)
    call allocate_array(SDYXS,1,maxnxs,1,maxizs,'SDYXS',ierr)
    call allocate_array(SDYYS,-maxnys,maxnys,1,maxizs,'SDYYS',ierr)
    call allocate_array(SCTXS,1,maxnxs,1,maxnls,'SCTXS',ierr)
    call allocate_array(SCTYS,-maxnys,maxnys,1,maxnls,'SCTYS',ierr)


  end subroutine allocate_mod_dynam3


  subroutine deallocate_mod_dynam3
    use mod_params
    use allocate_arrays
    implicit none

    deallocate(TIZS )
    deallocate(TIZ3 )
    deallocate(ZEFFS)
    deallocate(SAVES)
    deallocate(POWLS)
    deallocate(LINES)
    deallocate(PLRPS)
    deallocate(POWL3)
    deallocate(LINE3)
    deallocate(PLRP3)
    deallocate(LIM5 )
    deallocate(SDTXS)
    deallocate(SDTYS)
    deallocate(SDYXS)
    deallocate(SDYYS)
    deallocate(SCTXS)
    deallocate(SCTYS)


  end subroutine deallocate_mod_dynam3


end module mod_dynam3
