module mod_comt2


  !C                                                                               
  !      COMMON /COMT2/  CEYS,CVHYS,CTEMBS,CRNBS,QEDGES,CYSCLS,CFVHXS,             
  !     >                QTANS,CFEXZS,CFIZS,CFPS,CFSS,CFTS,CFRCS,CPCHS,            
  !     >                CPRCS,CCCFPS,CXAFS,CXBFS,CFCXS,CNHS,CXCFS,                
  !     >                QTEMBS,QRNBS,CTOLDS,CYMFPS,CYMFSS,QDISTS,       
  !     >                CTEMBSI,QTEMBSI,CVS,CTIGS,CTEGS,CXDPS
  !      REAL            CEYS  (MAXQYS),CVHYS (MAXQYS),QDISTS(-MAXQXS:0,2)         
  !      REAL            QTEMBS(-MAXQXS:1,2),    QRNBS (-MAXQXS:1,2)               
  !      REAL            QTEMBSI(-MAXQXS:1,2) 
  !      REAL            QEDGES(-MAXQXS:0,2),    CYSCLS(-MAXQXS:0)                 
  !      REAL            QTANS (-MAXQXS:0,2),    CXAFS (-MAXQXS:MAXQXS,3)          
  !
  !      REAL            CTEMBS(MAXNXS,-MAXNYS:MAXNYS)                             
  !      REAL            CTEMBSI(MAXNXS,-MAXNYS:MAXNYS)
  !      REAL            CTEGS(MAXNXS,-MAXNYS:MAXNYS)
  !      REAL            CTIGS(MAXNXS,-MAXNYS:MAXNYS)
  !      REAL            CRNBS (MAXNXS,-MAXNYS:MAXNYS)                             
  !      REAL            CFVHXS(MAXNXS,-MAXNYS:MAXNYS)                             
  !      REAL            CFEXZS(MAXNXS,-MAXNYS:MAXNYS,MAXIZS)                      
  !      REAL            CFIZS (MAXNXS,-MAXNYS:MAXNYS,0:MAXIZS)                    
  !      REAL            CFPS  (MAXNXS,-MAXNYS:MAXNYS,MAXIZS)                      
  !      REAL            CFSS  (MAXNXS,-MAXNYS:MAXNYS,MAXIZS)                      
  !      REAL            CFTS  (MAXNXS,-MAXNYS:MAXNYS,MAXIZS)                      
  !      REAL            CFRCS (MAXNXS,-MAXNYS:MAXNYS,0:MAXIZS)                    
  !      REAL            CPCHS (MAXNXS,-MAXNYS:MAXNYS,0:MAXIZS)                    
  !      REAL            CPRCS (MAXNXS,-MAXNYS:MAXNYS,0:MAXIZS)                    
  !      REAL            CFCXS (MAXNXS,-MAXNYS:MAXNYS,0:MAXIZS)                    
  !      REAL            CCCFPS(MAXNXS,-MAXNYS:MAXNYS,MAXIZS)                      
  !      REAL            CNHS  (MAXNXS,-MAXNYS:MAXNYS)
  !
  !      REAL            CYMFPS(-MAXQXS:1,2)        
  !
  !      REAL            CXBFS (-MAXQXS:MAXQXS,3), CXCFS(-MAXQXS:MAXQXS,3)         
  !      REAL            CXDPS (-MAXQXS:MAXQXS,3)
  !      REAL            CTOLDS(MAXNXS,MAXIZS), CYMFSS(-MAXQXS:1,2)                
  !      REAL            CVS(-MAXQXS:MAXQXS,2)

  implicit none
  private

  

  REAL,allocatable,public :: CEYS(:),CVHYS(:),QDISTS (:,:),QTEMBS (:,:),&
       QRNBS  (:,:),QTEMBSI(:,:),QEDGES (:,:),CYSCLS (:),QTANS  (:,:),&    
       CXAFS  (:,:),CTEMBS (:,:),CTEMBSI(:,:),CTEGS  (:,:),CTIGS  (:,:),&
       CRNBS  (:,:),CFVHXS (:,:),CFEXZS (:,:,:),CFIZS  (:,:,:),CFPS   (:,:,:),&
       CFSS   (:,:,:),CFTS(:,:,:),CFRCS  (:,:,:),CPCHS  (:,:,:),CPRCS  (:,:,:),&
       CFCXS  (:,:,:),CCCFPS (:,:,:),CXBFS  (:,:),CXCFS  (:,:),CXDPS  (:,:),&
       CNHS   (:,:),CYMFPS (:,:),CTOLDS (:,:),CYMFSS (:,:),CVS(:,:)


  integer,public:: vel_efield_opt = 0
  real,allocatable,public:: efield(:,:,:),velplasma(:,:,:)

  real,public :: sf_tau = 1.0, sf_vdiff = 1.0
  
  public:: allocate_mod_comt2,deallocate_mod_comt2




contains

  subroutine allocate_mod_comt2
    use mod_params
    use allocate_arrays

    implicit none
    integer :: ierr


    call allocate_array(CEYS   ,maxqys,'CEYS   ',ierr)
    call allocate_array(CVHYS  ,maxqys,'CVHYS  ',ierr)

    call allocate_array(QDISTS ,-maxqxs,0,1,2,'QDISTS ',ierr)
    call allocate_array(QEDGES ,-maxqxs,0,1,2,'QEDGES ',ierr)
    call allocate_array(QTANS  ,-maxqxs,0,1,2,'QTANS  ',ierr)

    call allocate_array(QRNBS  ,-maxqxs,1,1,2,'QRNBS  ',ierr)
    call allocate_array(QTEMBS ,-maxqxs,1,1,2,'QTEMBS ',ierr)
    call allocate_array(QTEMBSI,-maxqxs,1,1,2,'QTEMBSI',ierr)

    call allocate_array(CYSCLS ,-maxqxs,'CYSCLS ',0,ierr)

    call allocate_array(CXAFS  ,-maxqxs,maxqxs,1,3,'CXAFS  ',ierr)


    ! adding poloidal zones to all background arrays
    call allocate_array(CTEMBS ,1,maxnxs,-maxnys,maxnys,'CTEMBS ',ierr)
    call allocate_array(CTEMBSI,1,maxnxs,-maxnys,maxnys,'CTEMBSI',ierr)
    call allocate_array(CTEGS  ,1,maxnxs,-maxnys,maxnys,'CTEGS  ',ierr)
    call allocate_array(CTIGS  ,1,maxnxs,-maxnys,maxnys,'CTIGS  ',ierr)
    call allocate_array(CRNBS  ,1,maxnxs,-maxnys,maxnys,'CRNBS  ',ierr)
    call allocate_array(CFVHXS ,1,maxnxs,-maxnys,maxnys,'CFVHXS ',ierr)

    call allocate_array(efield   ,1,maxnxs,-maxnys,maxnys,1,maxpzone,'EFIELD',ierr)
    call allocate_array(velplasma,1,maxnxs,-maxnys,maxnys,1,maxpzone,'VELPASMA',ierr)

    call allocate_array(CFEXZS ,1,maxnxs,-maxnys,maxnys, 1,maxizs,'CFEXZS ',ierr)
    call allocate_array(CFIZS  ,1,maxnxs,-maxnys,maxnys, 0,maxizs,'CFIZS  ',ierr)
    call allocate_array(CFPS   ,1,maxnxs,-maxnys,maxnys, 1,maxizs,'CFPS   ',ierr)
    call allocate_array(CFSS   ,1,maxnxs,-maxnys,maxnys, 1,maxizs,'CFSS   ',ierr)
    call allocate_array(CFTS   ,1,maxnxs,-maxnys,maxnys, 1,maxizs,'CFTS   ',ierr)
    call allocate_array(CFRCS  ,1,maxnxs,-maxnys,maxnys, 0,maxizs,'CFRCS  ',ierr)
    call allocate_array(CPCHS  ,1,maxnxs,-maxnys,maxnys, 0,maxizs,'CPCHS  ',ierr)
    call allocate_array(CPRCS  ,1,maxnxs,-maxnys,maxnys, 0,maxizs,'CPRCS  ',ierr)
    call allocate_array(CFCXS  ,1,maxnxs,-maxnys,maxnys, 0,maxizs,'CFCXS  ',ierr)
    call allocate_array(CCCFPS ,1,maxnxs,-maxnys,maxnys, 1,maxizs,'CCCFPS ',ierr)
    call allocate_array(CNHS   ,1,maxnxs,-maxnys,maxnys,'CNHS   ',ierr)


    call allocate_array(CTOLDS ,maxnxs,maxizs,'CTOLDS ',ierr)
    call allocate_array(CXBFS  ,-maxqxs,maxqxs,1,3,'CXBFS  ',ierr)
    call allocate_array(CXCFS  ,-maxqxs,maxqxs,1,3,'CXCFS  ',ierr)
    call allocate_array(CXDPS  ,-maxqxs,maxqxs,1,3,'CXDPS  ',ierr)

    call allocate_array(CYMFPS ,-maxqxs,1,1,2,'CYMFPS ',ierr)
    call allocate_array(CYMFSS ,-maxqxs,1,1,2,'CYMFSS ',ierr)
    call allocate_array(CVS    ,-maxqxs,maxqxs,1,2,'CVS    ',ierr)


  end subroutine allocate_mod_comt2


  subroutine deallocate_mod_comt2
    use mod_params
    use allocate_arrays
    implicit none

    deallocate(efield)
    deallocate(velplasma)

    deallocate(CEYS   )
    deallocate(CVHYS  )
    deallocate(QDISTS )
    deallocate(QTEMBS )
    deallocate(QRNBS  )
    deallocate(QTEMBSI)
    deallocate(QEDGES )
    deallocate(CYSCLS )
    deallocate(QTANS  )
    deallocate(CXAFS  )
    deallocate(CTEMBS )
    deallocate(CTEMBSI)
    deallocate(CTEGS  )
    deallocate(CTIGS  )
    deallocate(CRNBS  )
    deallocate(CFVHXS )
    deallocate(CFEXZS )
    deallocate(CFIZS  )
    deallocate(CFPS   )
    deallocate(CFSS   )
    deallocate(CFTS   )
    deallocate(CFRCS  )
    deallocate(CPCHS  )
    deallocate(CPRCS  )
    deallocate(CFCXS  )
    deallocate(CCCFPS )
    deallocate(CXBFS  )
    deallocate(CXCFS  )
    deallocate(CXDPS  )
    deallocate(CNHS   )
    deallocate(CYMFPS )
    deallocate(CTOLDS )
    deallocate(CYMFSS )
    deallocate(CVS    )


  end subroutine deallocate_mod_comt2



end module mod_comt2
