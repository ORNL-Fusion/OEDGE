module mod_comxyt

  use mod_params


!C                                                                               
!      COMMON /COMXYT/ NXS,NYS,NQXSO,NQXSI,NQYS,NTS,XSCALO,XSCALI,YSCALE,        
!     >                NY3D,CSTMAX,XS,YS,QXS,QYS,IQXS,IQYS,XWIDS,YWIDS,          
!     >                PS,XCYLS,DWELTS,DWELFS,CTIMES,PWIDS,QS,XOUTS,             
!     >                YOUTS,DELPS,CTMAX,CXDEEP
!      INTEGER         NXS,NYS,NQXSO,NQXSI,NQYS,NTS,IQYS(MAXNYS),NY3D            
!      INTEGER         IQXS(MAXNXS)                                              
!      REAL            CTIMES(MAXNTS+1,0:MAXIZS),CSTMAX,XOUTS(MAXNXS)            
!      REAL            XS(MAXNXS),YS(MAXNYS),QYS(MAXQYS),XSCALI,YSCALE           
!      REAL            QXS(-MAXQXS:MAXQXS),XSCALO                                
!      REAL            XWIDS(MAXNXS),YWIDS(MAXNYS),PS(-MAXNPS:MAXNPS)            
!      REAL            XCYLS(MAXNXS),DWELTS(0:MAXIZS),DWELFS(MAXNTS)             
!      REAL            PWIDS(-MAXNPS:MAXNPS),QS(-MAXQXS:MAXQXS)                  
!      REAL            YOUTS(-MAXNYS:MAXNYS),DELPS(MAXNXS,MAXNYS)
!      REAL            CTMAX,CXDEEP      
!
  
  
  implicit none
  private


      INTEGER,public::         NXS,NYS,NQXSO,NQXSI,NQYS,NTS,IQYS(MAXNYS),NY3D            
      INTEGER,public::         IQXS(MAXNXS)                                              
      REAL,public::            CTIMES(MAXNTS+1,0:MAXIZS),CSTMAX,XOUTS(MAXNXS)            
      REAL,public::            XS(MAXNXS),YS(MAXNYS),QYS(MAXQYS),XSCALI,YSCALE           
      REAL,public::            QXS(-MAXQXS:MAXQXS),XSCALO                                
      REAL,public::            XWIDS(MAXNXS),YWIDS(MAXNYS)
      REAL,public::            XCYLS(MAXNXS),DWELTS(0:MAXIZS),DWELFS(MAXNTS)             
      REAL,public::            QS(-MAXQXS:MAXQXS)                  
      REAL,public::            YOUTS(-MAXNYS:MAXNYS),DELPS(MAXNXS,MAXNYS)
      REAL,public::            CTMAX,CXDEEP      

      !REAL,public::            PWIDS(-MAXNPS:MAXNPS),PS(-MAXNPS:MAXNPS)                  

      real, allocatable,public :: pwids(:),ps(:)

      integer, public:: nsurf, npbins
      real, allocatable,public:: pbin_bnds(:),surf_bnds(:,:)
      integer, allocatable,public:: pzone(:)
      
  public :: allocate_mod_comxyt, deallocate_mod_comxyt


contains

  subroutine allocate_mod_comxyt
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(DTEV  ,maxnxs,'DTEV',ierr)
    ! 1D arrays with both bounds specified have a different argumemnt order to differentiate signatures. 
    call allocate_array(ps,-maxnps,'P bin widths',maxnps,ierr)
    call allocate_array(pwids,-maxnps,'P bin widths',maxnps,ierr)
    call allocate_array(pzone,-maxnps,'P bin widths',maxnps,ierr)
    call allocate_array(pbin_bnds,2*maxnps+1,'P bin boundaries',ierr)  ! allow extra one to include 0.0 if desired
    call allocate_array(surf_bnds,max_nsurf,2,'Poilidal limiter surface bounds',ierr)  ! Surface is defined as P1 to P2 - sets should not overlap
    

  end subroutine allocate_mod_comxyt


  subroutine deallocate_mod_comxyt
    use mod_params
    use allocate_arrays
    implicit none

    deallocate(ps)
    deallocate(pwids)
    deallocate(pzone)
    deallocate(pbin_bnds)
    deallocate(surf_bnds)

  end subroutine deallocate_mod_comxyt



end module mod_comxyt
