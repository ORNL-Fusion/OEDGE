module mod_comnet

  !c     -*-Fortran-*-
  !C                                                                               
  !      COMMON /COMNET/ OYS,ODS,OYOUTS,ODOUTS,OYMAX,ODMAX,                        
  !     >                OYWIDS,ODWIDS,NEROXS,NEROYS,NERODS,DEPS,WALLS,
  !     >                OYCOORD,OYIQX,CDFLUX,OYMAX2,nerods3
  !      REAL    OYS(MAXOS),ODS(MAXOS),OYOUTS(MAXOS),ODOUTS(MAXOS)                 
  !      REAL    OYMAX,ODMAX,NEROXS(MAXNXS,5,3),NEROYS(MAXOS,6)
  !      ! jdemod - code expects nerods3 to be indexed to 6 - not sure why it was only 5
  !      REAL    nerods3(maxos,-maxnps:maxnps,6)                    
  !      !REAL    nerods3(maxos,-maxnps:maxnps,5)                    
  !      REAL    OYWIDS(MAXOS),ODWIDS(MAXOS),NERODS(MAXOS,5)                       
  !      REAL    WALLS(-MAXNYS:MAXNYS,-2:MAXIZS+1),DEPS(MAXNXS,MAXIZS+1,3)           
  !      REAL    OYCOORD(MAXOS,3),CDFLUX(MAXOS,3),OYMAX2(2)
  !      INTEGER OYIQX(MAXOS)


  implicit none

  private


  REAL,allocatable,public ::  OYS(:),ODS(:),OYOUTS(:),ODOUTS(:)       

  real,allocatable,public ::      NEROXS(:,:,:),NEROYS(:,:),nerods3(:,:,:)    

  real,allocatable,public ::      OYWIDS(:),ODWIDS(:),NERODS(:,:),WALLS(:,:),DEPS(:,:,:)

  real,allocatable,public ::      OYCOORD(:,:),CDFLUX(:,:)


  REAL,public :: OYMAX,ODMAX,OYMAX2(2)

  INTEGER,allocatable,public :: OYIQX(:)

  public:: allocate_mod_comnet,deallocate_mod_comnet

contains


  subroutine allocate_mod_comnet
    use mod_params
    use allocate_arrays

    implicit none
    integer :: ierr

    call allocate_array(OYS    ,maxos,'OYS    ',ierr)
    call allocate_array(ODS    ,maxos,'ODS    ',ierr)
    call allocate_array(OYOUTS ,maxos,'OYOUTS ',ierr)
    call allocate_array(ODOUTS ,maxos,'ODOUTS ',ierr)
    call allocate_array(OYWIDS ,maxos,'OYWIDS ',ierr)
    call allocate_array(ODWIDS ,maxos,'ODWIDS ',ierr)

    call allocate_array(OYIQX  ,maxos,'OYIQX  ',ierr)
    call allocate_array(OYCOORD,maxos,3,'OYCOORD',ierr)
    call allocate_array(CDFLUX ,maxos,3,'CDFLUX ',ierr)

    call allocate_array(NEROXS ,maxnxs,5,3,'NEROXS ',ierr)
    call allocate_array(NEROYS ,maxos,6,'NEROYS ',ierr)

    call allocate_array(nerods3,1,maxos,-maxnps,maxnps,1,6,'nerods3',ierr)

    call allocate_array(NERODS ,maxos,5,'NERODS ',ierr)

    call allocate_array(WALLS  ,-maxnys,maxnys,-2,maxizs+1,'WALLS  ',ierr)
    call allocate_array(DEPS   ,maxnxs,maxizs+1,3,'DEPS   ',ierr)



  end subroutine allocate_mod_comnet


  subroutine deallocate_mod_comnet
    use mod_params
    use allocate_arrays
    implicit none

    deallocate(OYS    )
    deallocate(ODS    )
    deallocate(OYOUTS )
    deallocate(ODOUTS )
    deallocate(NEROXS )
    deallocate(NEROYS )
    deallocate(nerods3)
    deallocate(OYWIDS )
    deallocate(ODWIDS )
    deallocate(NERODS )
    deallocate(WALLS  )
    deallocate(DEPS   )
    deallocate(OYCOORD)
    deallocate(CDFLUX )
    deallocate(OYIQX  )



  end subroutine deallocate_mod_comnet



end module mod_comnet
