module mod_cneut

  !                                                                               
  !      COMMON /CNEUT/  SPUTYS,XATIZS,YATIZS,PATIZS,VINS  ,SNEWS ,                
  !     >                XPRODS,YPRODS,PPRODS,RMAXS ,RANVA ,RANVB ,RANVC,          
  !     >                FLUX1 ,ENEGY1,YIELD1,FY1   ,FYCUM ,RMAX1 ,RANVD,          
  !     >                FLUX2 ,ENEGY2,YIELD2,FY2   ,FSPLIT,RMAX2                  
  !REAL SPUTYS(MAXIMP),XATIZS(MAXIMP),YATIZS(MAXIMP),PATIZS(MAXIMP)          
  !REAL VINS  (MAXIMP),XPRODS(MAXIMP),YPRODS(MAXIMP)                         
  !REAL PPRODS(MAXIMP),RMAXS (MAXIMP),RANVA (MAXIMP),RANVB (MAXIMP)          
  !REAL RANVC (MAXIMP),SNEWS (MAXIMP),RANVD (MAXIMP)                         
  !REAL FLUX1 (-MAXQXS:1,2), ENEGY1(-MAXQXS:1,2), YIELD1(-MAXQXS:1,2)        
  !REAL FY1   (-MAXQXS:1,2), FYCUM (-MAXQXS:1,2), RMAX1 (-MAXQXS:1,2)        
  !REAL FLUX2 (-MAXQXS:1,2), ENEGY2(-MAXQXS:1,2), YIELD2(-MAXQXS:1,2)        
  !REAL FY2   (-MAXQXS:1,2), FSPLIT(-MAXQXS:1,2), RMAX2 (-MAXQXS:1,2)        


  implicit none


  private

  REAL,allocatable,public ::  SPUTYS(:),XATIZS(:),YATIZS(:),PATIZS(:),&
       VINS  (:),XPRODS(:),YPRODS(:),&
       PPRODS(:),RMAXS (:),RANVA (:),RANVB (:),&
       RANVC (:),SNEWS (:),RANVD (:),&                         
       FLUX1 (:,:), ENEGY1(:,:), YIELD1(:,:),&        
       FY1   (:,:), FYCUM (:,:), RMAX1 (:,:),&        
       FLUX2 (:,:), ENEGY2(:,:), YIELD2(:,:),&        
       FY2   (:,:), FSPLIT(:,:), RMAX2 (:,:)        

  public:: allocate_mod_cneut,deallocate_mod_cneut


contains


  subroutine allocate_mod_cneut
    use mod_params
    use allocate_arrays

    implicit none
    integer :: ierr

    call allocate_array(sputys,maximp,'SPUTYS',ierr)
    call allocate_array(xatizs,maximp,'XATIZS',ierr)
    call allocate_array(yatizs,maximp,'YATIZS',ierr)
    call allocate_array(patizs,maximp,'PATIZS',ierr)
    call allocate_array(vins,maximp,'VINS',ierr)
    call allocate_array(xprods,maximp,'XPRODS',ierr)
    call allocate_array(yprods,maximp,'YPRODS',ierr)
    call allocate_array(pprods,maximp,'PPRODS',ierr)
    call allocate_array(rmaxs,maximp,'RMAXS',ierr)
    call allocate_array(ranva,maximp,'RANVA',ierr)
    call allocate_array(ranvb,maximp,'RANVB',ierr)
    call allocate_array(ranvc,maximp,'RANVC',ierr)
    call allocate_array(ranvd,maximp,'RANVD',ierr)
    call allocate_array(snews,maximp,'SNEWS',ierr)

    call allocate_array(flux1,-maxqxs,1,1,2,'FLUX1',ierr)
    call allocate_array(enegy1,-maxqxs,1,1,2,'ENEGY1',ierr)
    call allocate_array(yield1,-maxqxs,1,1,2,'YIELD1',ierr)
    call allocate_array(fy1,-maxqxs,1,1,2,'FY1',ierr)
    call allocate_array(fycum,-maxqxs,1,1,2,'FYCUM',ierr)
    call allocate_array(rmax1,-maxqxs,1,1,2,'RMAX1',ierr)

    call allocate_array(flux2,-maxqxs,1,1,2,'FLUX2',ierr)
    call allocate_array(enegy2,-maxqxs,1,1,2,'ENEGY2',ierr)
    call allocate_array(yield2,-maxqxs,1,1,2,'YIELD2',ierr)
    call allocate_array(fy2,-maxqxs,1,1,2,'FY2',ierr)
    call allocate_array(fsplit,-maxqxs,1,1,2,'FSPLIT',ierr)
    call allocate_array(rmax2,-maxqxs,1,1,2,'RMAX2',ierr)



  end subroutine allocate_mod_cneut


  subroutine deallocate_mod_cneut
    use mod_params
    use allocate_arrays
    implicit none

    deallocate(sputys)
    deallocate(xatizs)
    deallocate(yatizs)
    deallocate(patizs)
    deallocate(vins)
    deallocate(xprods)
    deallocate(yprods)
    deallocate(pprods)
    deallocate(rmaxs)
    deallocate(ranva)
    deallocate(ranvb)
    deallocate(ranvc)
    deallocate(ranvd)
    deallocate(snews)

    deallocate(flux1)
    deallocate(enegy1)
    deallocate(yield1)
    deallocate(fy1)
    deallocate(fycum)
    deallocate(rmax1)

    deallocate(flux2)
    deallocate(enegy2)
    deallocate(yield2)
    deallocate(fy2)
    deallocate(fsplit)
    deallocate(rmax2)


  end subroutine deallocate_mod_cneut


end module mod_cneut
