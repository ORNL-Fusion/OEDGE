module mod_cioniz
  use debug_options
  implicit none
  !include 'cioniz'

  private

  !     -*-Fortran-*-
  !                                                                      
  !      COMMON /CIONIZ/ C215A,C350A,C350B,KFIZS,KFRCS,KPCHS,KPRCS,KFCXS,  
  !     >  KNHS,KFPS,KFSS,KFTS,KKKFPS,KPLOS,KFSSMOD,kpmtc,kpizs,
  !     >  lmfp_mtc
  !      save /cioniz/

  !      REAL C215A,C350A,C350B,                                           
  !     >  KFIZS(MAXNKS,MAXNRS,0:MAXIZS),KFRCS(MAXNKS,MAXNRS,0:MAXIZS),    
  !     >  KPCHS(MAXNKS,MAXNRS,0:MAXIZS),KPRCS(MAXNKS,MAXNRS,0:MAXIZS),    
  !     >  KFCXS(MAXNKS,MAXNRS,0:MAXIZS),KNHS (MAXNKS,MAXNRS),             
  !     >  KFPS (MAXNKS,MAXNRS,MAXIZS)  ,KFSS (MAXNKS,MAXNRS,MAXIZS),      
  !     >  KFTS (MAXNKS,MAXNRS,MAXIZS) ,KKKFPS(MAXNKS,MAXNRS,MAXIZS),      
  !     >  KPLOS(MAXNKS,MAXNRS,MAXIZS) ,kfssmod(maxnks,maxnrs),
  !     >  kpmtc(maxnks,maxnrs),kpizs(maxnks,maxnrs),
  !     >  lmfp_mtc(maxnks,maxnrs)

  !      REAL C215A,C350A,C350B,                                           
  !     >  KFIZS(MAXNKS,MAXNRS,0:MAXIZS),KFRCS(MAXNKS,MAXNRS,0:MAXIZS),    
  !     >  KPCHS(MAXNKS,MAXNRS,0:MAXIZS),KPRCS(MAXNKS,MAXNRS,0:MAXIZS),    
  !     >  KFCXS(MAXNKS,MAXNRS,0:MAXIZS),KNHS (MAXNKS,MAXNRS),             
  !     >  KFPS (MAXNKS,MAXNRS,MAXIZS)  ,KFSS (MAXNKS,MAXNRS,MAXIZS),      
  !     >  KFTS (MAXNKS,MAXNRS,MAXIZS) ,KKKFPS(MAXNKS,MAXNRS,MAXIZS),      
  !     >  KPLOS(MAXNKS,MAXNRS,MAXIZS) ,kfssmod(maxnks,maxnrs),
  !     >  kpmtc(maxnks,maxnrs),kpizs(maxnks,maxnrs),
  !     >  lmfp_mtc(maxnks,maxnrs)


  REAL, public ::  C215A,C350A,C350B

  real, allocatable, public ::  KFIZS(:,:,:),KFRCS(:,:,:),&    
       KPCHS(:,:,:),KPRCS(:,:,:),KFCXS(:,:,:),KNHS (:,:),KFPS (:,:,:)  ,KFSS (:,:,:), &
       KFTS (:,:,:) ,KKKFPS(:,:,:),KPLOS(:,:,:) ,kfssmod(:,:),kpmtc(:,:),kpizs(:,:), lmfp_mtc(:,:)



  public :: allocate_cioniz, deallocate_cioniz



contains

  subroutine allocate_cioniz
    use global_parameters
    use allocate_arrays

    implicit none

    integer :: ierr

    call pr_trace('MOD_CIONIZ','ALLOCATE')

    call allocate_array(kfizs,1,maxnks,1,maxnrs,0,maxizs,"KFIZS",ierr)
    call allocate_array(kfrcs,1,maxnks,1,maxnrs,0,maxizs,"KFRCS",ierr)
    call allocate_array(kpchs,1,maxnks,1,maxnrs,0,maxizs,"KPCHS",ierr)
    call allocate_array(kprcs,1,maxnks,1,maxnrs,0,maxizs,"KPRCS",ierr)
    call allocate_array(kfcxs,1,maxnks,1,maxnrs,0,maxizs,"KFCXS",ierr)

    call allocate_array(knhs,maxnks,maxnrs,"KNHS",ierr)
    call allocate_array(kfps,maxnks,maxnrs,maxizs,"KFPS",ierr)
    call allocate_array(kfss,maxnks,maxnrs,maxizs,"KFSS",ierr)
    call allocate_array(kfts,maxnks,maxnrs,maxizs,"KFTS",ierr)
    call allocate_array(kkkfps,maxnks,maxnrs,maxizs,"KKKFPS",ierr)
    call allocate_array(kplos,maxnks,maxnrs,maxizs,"KPLOS",ierr)
    call allocate_array(kfssmod,maxnks,maxnrs,"KFSSMOD",ierr)
    call allocate_array(kpmtc,maxnks,maxnrs,"KPMTC",ierr)
    call allocate_array(kpizs,maxnks,maxnrs,"KPIZS",ierr)
    call allocate_array(lmfp_mtc,maxnks,maxnrs,"LMFP_MTC",ierr)


  end subroutine allocate_cioniz

  subroutine deallocate_cioniz

    implicit none

    call pr_trace('MOD_CIONIZ','DEALLOCATE')

    deallocate(kfizs)
    deallocate(kfrcs)
    deallocate(kpchs)
    deallocate(kprcs)
    deallocate(kfcxs)

    deallocate(knhs)
    deallocate(kfps)
    deallocate(kfss)
    deallocate(kfts)
    deallocate(kkkfps)
    deallocate(kplos)
    deallocate(kfssmod)
    deallocate(kpmtc)
    deallocate(kpizs)
    deallocate(lmfp_mtc)


  end subroutine deallocate_cioniz




end module mod_cioniz
