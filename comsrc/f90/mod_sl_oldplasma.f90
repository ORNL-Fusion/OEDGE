module mod_sl_oldplasma

  implicit none
  ! jdemod
  ! reorg of SL arrays dependent on parameters
  ! contents of oldplasma common block - including allocation

     ! COMMON /OLDPLASMA/ oldknbs ,oldktebs ,oldktibs ,oldkvhs ,
     !                   oldknbs2,oldktebs2,oldktibs2,oldkvhs2
  public
  
  REAL,allocatable :: &
           oldknbs  (:,:),oldktebs (:,:),&
           oldktibs (:,:),oldkvhs  (:,:),&
           oldktebs2(:,:),oldktibs2(:,:),&
           oldknbs2 (:,:),oldkvhs2 (:,:)

  !REAL    sval1,quant(MAXNKS,MAXNRS,4),tquant(MAXNDS,4)
  REAL,allocatable ::  quant(:,:,:),tquant(:,:)
  

contains
  
      SUBROUTINE MirrorOldPlasma(te,ti,ne,vb)
      use mod_params
      IMPLICIT none

      REAL te(MAXNKS,MAXNRS),ti(MAXNKS,MAXNRS),ne(MAXNKS,MAXNRS),&
           vb(MAXNKS,MAXNRS) 

!      COMMON /OLDPLASMA/ oldknbs ,oldktebs ,oldktibs ,oldkvhs ,
!     .                   oldknbs2,oldktebs2,oldktibs2,oldkvhs2
!      REAL oldktebs (MAXNKS,MAXNRS),oldktibs (MAXNKS,MAXNRS),
!     .     oldknbs  (MAXNKS,MAXNRS),oldkvhs  (MAXNKS,MAXNRS),
!     .     oldktebs2(MAXNKS,MAXNRS),oldktibs2(MAXNKS,MAXNRS),
!     .     oldknbs2 (MAXNKS,MAXNRS),oldkvhs2 (MAXNKS,MAXNRS)

      INTEGER ik,ir

      DO ir = 1, MAXNRS
        DO ik = 1, MAXNKS
          oldktebs(ik,ir) = te(ik,ir)
          oldktibs(ik,ir) = ti(ik,ir)
          oldknbs (ik,ir) = ne(ik,ir)
          oldkvhs (ik,ir) = vb(ik,ir)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE MirrorOldPlasma


  SUBROUTINE MirrorOldPlasma2(te,ti,ne,vb)
      use mod_params
      IMPLICIT none

      REAL te(MAXNKS,MAXNRS),ti(MAXNKS,MAXNRS),ne(MAXNKS,MAXNRS),&
           vb(MAXNKS,MAXNRS) 

!      COMMON /OLDPLASMA/ oldknbs ,oldktebs ,oldktibs ,oldkvhs ,
!     .                   oldknbs2,oldktebs2,oldktibs2,oldkvhs2
!      REAL oldktebs (MAXNKS,MAXNRS),oldktibs (MAXNKS,MAXNRS),
!     .     oldknbs  (MAXNKS,MAXNRS),oldkvhs  (MAXNKS,MAXNRS),
!     .     oldktebs2(MAXNKS,MAXNRS),oldktibs2(MAXNKS,MAXNRS),
!     .     oldknbs2 (MAXNKS,MAXNRS),oldkvhs2 (MAXNKS,MAXNRS)

      INTEGER ik,ir

      DO ir = 1, MAXNRS
        DO ik = 1, MAXNKS
          oldktebs2(ik,ir) = te(ik,ir)
          oldktibs2(ik,ir) = ti(ik,ir)
          oldknbs2 (ik,ir) = ne(ik,ir)
          oldkvhs2 (ik,ir) = vb(ik,ir)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE MirrorOldPlasma2


  subroutine allocate_mod_sl_oldplasma
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr
    call allocate_array(oldknbs,maxnks,maxnrs,'oldknbs',ierr)
    call allocate_array(oldktebs,maxnks,maxnrs,'oldknbs',ierr)
    call allocate_array(oldktibs,maxnks,maxnrs,'oldknbs',ierr)
    call allocate_array(oldkvhs,maxnks,maxnrs,'oldknbs',ierr)

    call allocate_array(oldknbs2,maxnks,maxnrs,'oldknbs2',ierr)
    call allocate_array(oldktebs2,maxnks,maxnrs,'oldknbs2',ierr)
    call allocate_array(oldktibs2,maxnks,maxnrs,'oldknbs2',ierr)
    call allocate_array(oldkvhs2,maxnks,maxnrs,'oldknbs2',ierr)

    call allocate_array(quant,maxnks,maxnrs,4,'quant',ierr)
    call allocate_array(tquant,maxnds,4,'tquant',ierr)
     

  end subroutine allocate_mod_sl_oldplasma
  

  subroutine deallocate_mod_sl_oldplasma
    implicit none
    if (allocated(oldknbs )) deallocate(oldknbs)
    if (allocated(oldktebs)) deallocate(oldktebs)
    if (allocated(oldktibs)) deallocate(oldktibs)
    if (allocated(oldkvhs )) deallocate(oldkvhs)

    if (allocated(oldknbs2 )) deallocate(oldknbs2)
    if (allocated(oldktebs2)) deallocate(oldktebs2)
    if (allocated(oldktibs2)) deallocate(oldktibs2)
    if (allocated(oldkvhs2 )) deallocate(oldkvhs2)

    if (allocated(quant )) deallocate(quant)
    if (allocated(tquant )) deallocate(tquant)

  end subroutine deallocate_mod_sl_oldplasma




end module mod_sl_oldplasma
