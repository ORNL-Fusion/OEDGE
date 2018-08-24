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
      REAL,public::            XWIDS(MAXNXS),YWIDS(MAXNYS),PS(-MAXNPS:MAXNPS)            
      REAL,public::            XCYLS(MAXNXS),DWELTS(0:MAXIZS),DWELFS(MAXNTS)             
      REAL,public::            PWIDS(-MAXNPS:MAXNPS),QS(-MAXQXS:MAXQXS)                  
      REAL,public::            YOUTS(-MAXNYS:MAXNYS),DELPS(MAXNXS,MAXNYS)
      REAL,public::            CTMAX,CXDEEP      


  
  public :: allocate_mod_comxyt, deallocate_mod_comxyt


contains

  subroutine allocate_mod_comxyt
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(DTEV  ,maxnxs,'DTEV',ierr)


  end subroutine allocate_mod_comxyt


  subroutine deallocate_mod_comxyt
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()

  end subroutine deallocate_mod_comxyt



end module mod_comxyt
