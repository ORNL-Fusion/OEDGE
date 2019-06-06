module mod_commv
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  ! common /commv/  crtabs, crtrcs, crvabs, cravav, crxmin, ctbs,cicizs, cifizs, cilizs,&
  !      cisizs, cicabs, cifabs, cifrin,cilabs, cisabs, cicuts, ciccol, cicrxa, cifrxa,&
  !      cflrin,cicrcs, cifrcs, cilrcs, cisrcs, cisrxa, citrxa, citrin,ctexs , cflrxa,&
  !      ckkmax, cicrin, cisrin, cxxx  , csss  ,cnnn  , ckkmin, cnnnx , clll  , clllx ,&
  !      cflrex, cnnns ,cllls , cnnnk , cnnnt , cnnnkt, cmmm  , cmmmx , cmmms ,cikrin, cktrin,&
  !      cicrno, cikrno, cktrno, cicrnj, cvvxc ,cvvzm , cvvxe , cvvxp , cvvxs , ciclos,&
  !      ciflos, cillos,cislos, cssss , cvvrm , cvvsm , cvvkm , cnorgs, cnorgr,cnorgz,&
  !      cirrno, cizrno, cisrno,cvvrefm,cvvnrfm,cvvfpref,cvvfpnrf,cieizs,citizs
  !
  ! save /commv/
  
  real,public :: crxmin,citrxa,cifrin,citrin,cisrin,cicrin,ckkmax,ciccol,cicrxa,cifrxa,&
       cisrxa,ckkmin,cikrin,cktrin,cicrno,cikrno,cktrno,cicrnj,cvvxc,cvvzm,cvvxe,cvvxp,&
       cvvxs,cvvrm,cvvkm,cvvsm,cirrno,cizrno,cisrno,cvvrefm,cvvnrfm,cvvfpref,cvvfpnrf
  real,public,allocatable :: crtabs(:),crtrcs(:),crvabs(:),ctbs(:),cravav(:),cicizs(:),&
       cifizs(:),cilizs(:),cisizs(:),cicabs(:),cifabs(:),ctexs(:),cilabs(:),cisabs(:),&
       cicuts(:),cicrcs(:),cifrcs(:),cilrcs(:),cisrcs(:),cxxx(:),csss(:),cnnn(:),cnnnx(:),&
       clll(:),clllx(:),cnnns(:),cllls(:),cnnnk(:),cnnnt(:),cnnnkt(:),cmmm(:),cmmmx(:),&
       cmmms(:),ciclos(:),ciflos(:),cillos(:),cislos(:),cssss(:),cnorgs(:),cnorgr(:),&
       cnorgz(:),cieizs(:),citizs(:)
  ! jdemod - moved from comtor
  real,public,allocatable :: cleakn(:,:)

  logical,public :: cflrxa,cflrin,cflrex

  public :: allocate_mod_commv,deallocate_mod_commv

contains

  subroutine allocate_mod_commv
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_commv','ALLOCATE')

    call allocate_array(crtabs,maxizs,'crtabs',ierr)
    call allocate_array(crtrcs,maxizs,'crtrcs',ierr)
    call allocate_array(crvabs,maxizs,'crvabs',ierr)
    call allocate_array(ctbs,maxizs,'ctbs',ierr)
    call allocate_array(cravav,maxizs,'cravav',ierr)
    call allocate_array(cicizs,maxizs,'cicizs',ierr)
    call allocate_array(cifizs,maxizs,'cifizs',ierr)
    call allocate_array(cilizs,maxizs,'cilizs',ierr)
    call allocate_array(cisizs,maxizs,'cisizs',ierr)
    call allocate_array(cicabs,maxizs,'cicabs',ierr)
    call allocate_array(cifabs,maxizs,'cifabs',ierr)
    call allocate_array(ctexs,10,'ctexs',ierr)
    call allocate_array(cilabs,maxizs,'cilabs',ierr)
    call allocate_array(cisabs,maxizs,'cisabs',ierr)
    call allocate_array(cicuts,maxizs,'cicuts',ierr)
    call allocate_array(cicrcs,maxizs,'cicrcs',ierr)
    call allocate_array(cifrcs,maxizs,'cifrcs',ierr)
    call allocate_array(cilrcs,maxizs,'cilrcs',ierr)
    call allocate_array(cisrcs,maxizs,'cisrcs',ierr)
    call allocate_array(cxxx,maxizs,'cxxx',ierr)
    call allocate_array(csss,maxizs,'csss',ierr)
    call allocate_array(cnnn,-1,'cnnn',maxizs,ierr)
    call allocate_array(cnnnx,maxizs,'cnnnx',ierr)
    call allocate_array(clll,-1,'clll',maxizs,ierr)
    call allocate_array(clllx,maxizs,'clllx',ierr)
    call allocate_array(cnnns,maxizs,'cnnns',ierr)
    call allocate_array(cllls,maxizs,'cllls',ierr)
    call allocate_array(cnnnk,maxizs,'cnnnk',ierr)
    call allocate_array(cnnnt,maxizs,'cnnnt',ierr)
    call allocate_array(cnnnkt,maxizs,'cnnnkt',ierr)
    call allocate_array(cmmm,-1,'cmmm',maxizs,ierr)
    call allocate_array(cmmmx,maxizs,'cmmmx',ierr)
    call allocate_array(cmmms,maxizs,'cmmms',ierr)
    call allocate_array(ciclos,maxizs,'ciclos',ierr)
    call allocate_array(ciflos,maxizs,'ciflos',ierr)
    call allocate_array(cillos,maxizs,'cillos',ierr)
    call allocate_array(cislos,maxizs,'cislos',ierr)
    call allocate_array(cssss,maxizs,'cssss',ierr)
    call allocate_array(cnorgs,maxizs,'cnorgs',ierr)
    call allocate_array(cnorgr,maxizs,'cnorgr',ierr)
    call allocate_array(cnorgz,maxizs,'cnorgz',ierr)
    call allocate_array(cieizs,0,'cieizs',maxizs,ierr)
    call allocate_array(citizs,0,'citizs',maxizs,ierr)
    call allocate_array(cleakn,maxpts,maxizs+1,'cleakn',ierr)

  end subroutine allocate_mod_commv


  subroutine deallocate_mod_commv
    implicit none

    call pr_trace('mod_commv','DEALLOCATE')

    if (allocated(crtabs)) deallocate(crtabs)
    if (allocated(crtrcs)) deallocate(crtrcs)
    if (allocated(crvabs)) deallocate(crvabs)
    if (allocated(ctbs)) deallocate(ctbs)
    if (allocated(cravav)) deallocate(cravav)
    if (allocated(cicizs)) deallocate(cicizs)
    if (allocated(cifizs)) deallocate(cifizs)
    if (allocated(cilizs)) deallocate(cilizs)
    if (allocated(cisizs)) deallocate(cisizs)
    if (allocated(cicabs)) deallocate(cicabs)
    if (allocated(cifabs)) deallocate(cifabs)
    if (allocated(ctexs)) deallocate(ctexs)
    if (allocated(cilabs)) deallocate(cilabs)
    if (allocated(cisabs)) deallocate(cisabs)
    if (allocated(cicuts)) deallocate(cicuts)
    if (allocated(cicrcs)) deallocate(cicrcs)
    if (allocated(cifrcs)) deallocate(cifrcs)
    if (allocated(cilrcs)) deallocate(cilrcs)
    if (allocated(cisrcs)) deallocate(cisrcs)
    if (allocated(cxxx)) deallocate(cxxx)
    if (allocated(csss)) deallocate(csss)
    if (allocated(cnnn)) deallocate(cnnn)
    if (allocated(cnnnx)) deallocate(cnnnx)
    if (allocated(clll)) deallocate(clll)
    if (allocated(clllx)) deallocate(clllx)
    if (allocated(cnnns)) deallocate(cnnns)
    if (allocated(cllls)) deallocate(cllls)
    if (allocated(cnnnk)) deallocate(cnnnk)
    if (allocated(cnnnt)) deallocate(cnnnt)
    if (allocated(cnnnkt)) deallocate(cnnnkt)
    if (allocated(cmmm)) deallocate(cmmm)
    if (allocated(cmmmx)) deallocate(cmmmx)
    if (allocated(cmmms)) deallocate(cmmms)
    if (allocated(ciclos)) deallocate(ciclos)
    if (allocated(ciflos)) deallocate(ciflos)
    if (allocated(cillos)) deallocate(cillos)
    if (allocated(cislos)) deallocate(cislos)
    if (allocated(cssss)) deallocate(cssss)
    if (allocated(cnorgs)) deallocate(cnorgs)
    if (allocated(cnorgr)) deallocate(cnorgr)
    if (allocated(cnorgz)) deallocate(cnorgz)
    if (allocated(cieizs)) deallocate(cieizs)
    if (allocated(citizs)) deallocate(citizs)
    if (allocated(cleakn)) deallocate(cleakn)

  end subroutine deallocate_mod_commv

end module mod_commv
