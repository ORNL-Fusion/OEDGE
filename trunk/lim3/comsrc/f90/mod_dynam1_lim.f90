module mod_dynam1


  !C                                                                               
  !      COMMON /DYNAM1/ DDLIMS,DDLIM3,DDTS,DDYS                                   
  !      DOUBLE PRECISION                                                          
  !     >     DDLIMS(MAXNXS,-MAXNYS:MAXNYS,-1:MAXIZS),                             
  !     >     DDLIM3(MAXNXS,-MAXY3D:MAXY3D,-1:MAXIZS,-MAXNPS:MAXNPS),              
  !     >     DDTS  (MAXNXS,-MAXNYS:MAXNYS,MAXIZS),                                
  !     >     DDYS  (MAXNXS,-MAXNYS:MAXNYS,MAXIZS)                                 



  implicit none

  private


  REAL*8, allocatable,public ::  DDLIMS(:,:,:),&                             
       DDLIM3(:,:,:,:),DDTS  (:,:,:), DDYS  (:,:,:)
!  REAL*8, allocatable,public ::  DDLIMS(:,:,:,:),&                             
!       DDLIM3(:,:,:,:),DDTS  (:,:,:,:), DDYS  (:,:,:)

  public:: allocate_mod_dynam1,deallocate_mod_dynam1


contains


  subroutine allocate_mod_dynam1
    use mod_params
    use allocate_arrays

    implicit none
    integer :: ierr


    call allocate_array(DDLIMS,1,maxnxs,-MAXNYS,MAXNYS,-1,maxizs,'DDLIMS',ierr)
    call allocate_array(DDLIM3,1,maxnxs,-MAXY3D,MAXY3D,-1,maxizs,-maxnps,maxnps,'DDLIM3',ierr)
    call allocate_array(DDTS  ,1,maxnxs,-MAXNYS,MAXNYS, 1,maxizs,'DDTS  ',ierr)
    call allocate_array(DDYS  ,1,maxnxs,-MAXNYS,MAXNYS, 1,maxizs,'DDYS  ',ierr)

    !call allocate_array(DDLIMS,1,maxnxs,-MAXNYS,MAXNYS,-1,maxizs,1,maxpzone,'DDLIMS',ierr)
    !call allocate_array(DDLIM3,1,maxnxs,-MAXY3D,MAXY3D,-1,maxizs,-maxnps,maxnps,'DDLIM3',ierr)
    !call allocate_array(DDTS  ,1,maxnxs,-MAXNYS,MAXNYS, 1,maxizs,1,maxpzone,'DDTS  ',ierr)
    !call allocate_array(DDYS  ,1,maxnxs,-MAXNYS,MAXNYS, 1,maxizs,'DDYS  ',ierr)

  end subroutine allocate_mod_dynam1


  subroutine deallocate_mod_dynam1
    use mod_params
    use allocate_arrays
    implicit none

    deallocate(DDLIMS)
    deallocate(DDLIM3)
    deallocate(DDTS  )
    deallocate(DDYS  )


  end subroutine deallocate_mod_dynam1

end module mod_dynam1
