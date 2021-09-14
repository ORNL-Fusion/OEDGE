module mod_commv

  use mod_params

!C                                                                               
!      COMMON /COMMV/  CRTABS, CRTRCS, CRVABS, CRAVAV, CRXMIN, CTBS,             
!     >                CICIZS, CIFIZS, CILIZS, CISIZS, CICABS, CIFABS,           
!     >                CILABS, CISABS, CICUTS, CICCOL, CICSTX, CICSTZ,           
!     >                CICRXA, CIFRXA, CICY2L, CIFY2L, COUTS,  CICAY0,           
!     >                CICA2L, CICAAW, CTEXS,  CTTT,                             
!     >                CFLRXA, CFLY2L,                                           
!     >                CICRCS, CIFRCS, CILRCS, CISRCS, CISRXA, CITRXA,
!     >                CIY0IZ, CI2LIZ, ZIY0IZ, ZI2LIZ,            
!     >                CTY0IZ, CT2LIZ, ZTY0IZ, ZT2LIZ            
!      REAL            CRTABS(MAXIZS), CRTRCS(MAXIZS), CRVABS(MAXIZS)            
!      REAL            CRAVAV(MAXIZS), CRXMIN, CITRXA, CTBS(MAXIZS)              
!      REAL            CICIZS(MAXIZS), CIFIZS(MAXIZS), CILIZS(MAXIZS)            
!      REAL            CISIZS(MAXIZS), CICABS(MAXIZS), CIFABS(MAXIZS)            
!      REAL            CILABS(MAXIZS), CISABS(MAXIZS), CICUTS(MAXIZS)            
!      REAL            CICCOL, CICSTX, CICSTZ, CICRXA, CIFRXA, CICY2L            
!      REAL            CIFY2L, CICAY0, CICA2L, CICAAW, CISRXA, CTTT              
!      REAL            CICRCS(MAXIZS), CIFRCS(MAXIZS), CILRCS(MAXIZS)            
!      REAL            CISRCS(MAXIZS), COUTS(MAXIZS),  CTEXS(10)
!      REAL            CIY0IZ(MAXIZS), CI2LIZ(MAXIZS), ZIY0IZ(MAXIZS)
!      REAL            ZI2LIZ(MAXIZS)                 
!      REAL            CTY0IZ(MAXIZS), CT2LIZ(MAXIZS), ZTY0IZ(MAXIZS)
!      REAL            ZT2LIZ(MAXIZS)                 
!      LOGICAL         CFLRXA, CFLY2L                                            
!

  
  implicit none
  private

      REAL,ALLOCATABLE,PUBLIC:: CRTABS(:), CRTRCS(:), CRVABS(:)            
      REAL,ALLOCATABLE,PUBLIC:: CRAVAV(:), CTBS(:)                              
      REAL,ALLOCATABLE,PUBLIC:: CICIZS(:), CIFIZS(:), CILIZS(:)            
      REAL,ALLOCATABLE,PUBLIC:: CISIZS(:), CICABS(:), CIFABS(:)            
      REAL,ALLOCATABLE,PUBLIC:: CILABS(:), CISABS(:), CICUTS(:)            

      REAL,public::            CRXMIN, CITRXA
      REAL,public::            CICCOL, CICSTX, CICSTZ, CICRXA, CIFRXA, CICY2L            
      REAL,public::            CIFY2L, CICAY0, CICA2L, CICAAW, CISRXA, CTTT              

      REAL,ALLOCATABLE,PUBLIC:: CICRCS(:), CIFRCS(:), CILRCS(:)            
      REAL,ALLOCATABLE,PUBLIC:: CISRCS(:), COUTS(:),  CTEXS(:)       
      REAL,ALLOCATABLE,PUBLIC:: CIY0IZ(:), CI2LIZ(:), ZIY0IZ(:)  
      REAL,ALLOCATABLE,PUBLIC:: ZI2LIZ(:)                                  
      REAL,ALLOCATABLE,PUBLIC:: CTY0IZ(:), CT2LIZ(:), ZTY0IZ(:)  
      REAL,ALLOCATABLE,PUBLIC:: ZT2LIZ(:)                                  

      LOGICAL,public::         CFLRXA, CFLY2L                                            


  public :: allocate_mod_commv, deallocate_mod_commv


contains

  subroutine allocate_mod_commv
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(DTEV  ,maxnxs,'DTEV',ierr)
     call allocate_array(CRTABS,MAXIZS,'CRTABS',ierr)
     call allocate_array(CRTRCS,MAXIZS,'CRTRCS',ierr)
     call allocate_array(CRVABS,MAXIZS,'CRVABS',ierr)            
     call allocate_array(CRAVAV,MAXIZS,'CRAVAV',ierr)
     call allocate_array(CTBS  ,MAXIZS,'CTBS  ',ierr)                              
     call allocate_array(CICIZS,MAXIZS,'CICIZS',ierr)
     call allocate_array(CIFIZS,MAXIZS,'CIFIZS',ierr)
     call allocate_array(CILIZS,MAXIZS,'CILIZS',ierr)            
     call allocate_array(CISIZS,MAXIZS,'CISIZS',ierr)
     call allocate_array(CICABS,MAXIZS,'CICABS',ierr)
     call allocate_array(CIFABS,MAXIZS,'CIFABS',ierr)            
     call allocate_array(CILABS,MAXIZS,'CILABS',ierr)
     call allocate_array(CISABS,MAXIZS,'CISABS',ierr)
     call allocate_array(CICUTS,MAXIZS,'CICUTS',ierr)            
     call allocate_array(CICRCS,MAXIZS,'CICRCS',ierr)
     call allocate_array(CIFRCS,MAXIZS,'CIFRCS',ierr)
     call allocate_array(CILRCS,MAXIZS,'CILRCS',ierr)  
     call allocate_array(CISRCS,MAXIZS,'CISRCS',ierr)
     call allocate_array(COUTS ,MAXIZS,'COUTS ',ierr)
     call allocate_array(CTEXS ,10    ,'CTEXS ',ierr)       
     call allocate_array(CIY0IZ,MAXIZS,'CIY0IZ',ierr)
     call allocate_array(CI2LIZ,MAXIZS,'CI2LIZ',ierr)
     call allocate_array(ZIY0IZ,MAXIZS,'ZIY0IZ',ierr)  
     call allocate_array(ZI2LIZ,MAXIZS,'ZI2LIZ',ierr)                                  
     call allocate_array(CTY0IZ,MAXIZS,'CTY0IZ',ierr)
     call allocate_array(CT2LIZ,MAXIZS,'CT2LIZ',ierr)
     call allocate_array(ZTY0IZ,MAXIZS,'ZTY0IZ',ierr)  
     call allocate_array(ZT2LIZ,MAXIZS,'ZT2LIZ',ierr)                                  

  end subroutine allocate_mod_commv


  subroutine deallocate_mod_commv
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()
     deallocate(CRTABS)
     deallocate(CRTRCS)
     deallocate(CRVABS)            
     deallocate(CRAVAV)
     deallocate(CTBS  )                              
     deallocate(CICIZS)
     deallocate(CIFIZS)
     deallocate(CILIZS)            
     deallocate(CISIZS)
     deallocate(CICABS)
     deallocate(CIFABS)            
     deallocate(CILABS)
     deallocate(CISABS)
     deallocate(CICUTS)            
     deallocate(CICRCS)
     deallocate(CIFRCS)
     deallocate(CILRCS)  
     deallocate(CISRCS)
     deallocate(COUTS )
     deallocate(CTEXS )       
     deallocate(CIY0IZ)
     deallocate(CI2LIZ)
     deallocate(ZIY0IZ)  
     deallocate(ZI2LIZ)                                  
     deallocate(CTY0IZ)
     deallocate(CT2LIZ)
     deallocate(ZTY0IZ)  
     deallocate(ZT2LIZ)                                  

  end subroutine deallocate_mod_commv



end module mod_commv
