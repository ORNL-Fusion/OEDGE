module mod_slout

  use mod_params


!      COMMON /DRAWCOM/ plottype,plotnorm,slopt,char,
!     .                 map1x,map2x,map1y,map2y,
!     .                 opt_xscale,opt_yscale,restoresolution,
!     .                 slopt2,slopt3,slopt4,slopt5
!
!      INTEGER plottype(MAXNGS),plotnorm(MAXNGS),slopt,
!     .        opt_xscale,opt_yscale,
!     .        slopt2,slopt3,slopt4,slopt5
!      LOGICAL restoresolution
!      REAL    map1x,map2x,map1y,map2y
!
!      CHARACTER*100    char(30)
!
!      COMMON /GRMCOM/ grm_opt,ngs2,elabs2,grm_shade,grm_cell,plottype2,
!     .                ylab2
!      INTEGER         grm_opt,ngs2(MAXNGS),plottype2(8,MAXNGS)
!      REAL            grm_shade(2,MAXNGS),grm_cell(0:MAXNXS,MAXNGS)
!      CHARACTER*36    elabs2   (8,MAXNGS)
!      CHARACTER*128   ylab2    (-30:30)
!
!
!      INTEGER    MAXSHOW
!      PARAMETER (MAXSHOW=10000)
!
!      COMMON /LOSCOM/ WGHT0,LOSOPT,nshow,rshow,zshow,ashow
!      INTEGER         LOSOPT,nshow
!      REAL            WGHT0(MAXNXS,MAXNYS),
!     .                rshow(MAXSHOW),zshow(MAXSHOW),ashow(MAXSHOW)
!
!      COMMON /TIMECOM/ qt
!      REAL             qt
!
!      COMMON /LOADCOM/ loadstep,stepopt,nsteplist,steplist
!      INTEGER          loadstep,stepopt,nsteplist,steplist(100)
!      
!      COMMON /GENCOM/ sldata
!      REAL            sldata
!
!c      COMMON /MACHPLOT/ machdat
!c      REAL              machdat(2,MAXNYS)
!
!
!      COMMON /NORMCOM/ nrmindex,nrmcalculate,nrmdata,nrmvalue,
!     .                 nrmi1,nrmi2,nrmtype,nrmr1,nrmr2,nrmnum,
!     .                 nrmstep,nrmcomment
!      INTEGER        nrmindex,nrmtype,nrmi1,nrmi2,nrmnum(1024),
!     .               nrmstep(1024)
!      LOGICAL        nrmcalculate
!      REAL           nrmdata(MAXTHE,2),nrmvalue(1024),nrmr1,nrmr2
!      CHARACTER*1024 nrmcomment(1024)


  implicit none
  private


      INTEGER,public:: plottype(MAXNGS),plotnorm(MAXNGS),slopt,&
              opt_xscale,opt_yscale,&
              slopt2,slopt3,slopt4,slopt5
      LOGICAL,public:: restoresolution
      REAL,public::    map1x,map2x,map1y,map2y

      CHARACTER*100,public::    char(30)

      INTEGER,public::         grm_opt,ngs2(MAXNGS),plottype2(8,MAXNGS)
      REAL,public::            grm_shade(2,MAXNGS)
      !REAL,public::            grm_shade(2,MAXNGS),grm_cell(0:MAXNXS,MAXNGS)
      REAL,allocatable,public::            grm_cell(:,:)
      CHARACTER*36,public::    elabs2   (8,MAXNGS)
      CHARACTER*128,public::   ylab2    (-30:30)


      INTEGER,public::    MAXSHOW
      PARAMETER (MAXSHOW=10000)

      INTEGER,public::         LOSOPT,nshow
      REAL,public::            rshow(MAXSHOW),zshow(MAXSHOW),ashow(MAXSHOW)
      !REAL,public::            WGHT0(MAXNXS,MAXNYS),rshow(MAXSHOW),zshow(MAXSHOW),ashow(MAXSHOW)
      REAL,allocatable,public::            WGHT0(:,:)

      REAL,public::             qt

      INTEGER,public::          loadstep,stepopt,nsteplist,steplist(100)
      
      REAL,public::            sldata

      INTEGER,public::        nrmindex,nrmtype,nrmi1,nrmi2,nrmnum(1024),nrmstep(1024)
      LOGICAL,public::        nrmcalculate
      REAL,public::           nrmdata(MAXTHE,2),nrmvalue(1024),nrmr1,nrmr2
      CHARACTER*1024,public:: nrmcomment(1024)


  
 
  public :: allocate_mod_slout, deallocate_mod_slout


contains

  subroutine allocate_mod_slout
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(DTEV  ,maxnxs,'DTEV',ierr)

    call allocate_array(grm_cell,0,maxnxs,1,maxngs,'GRM_CELL',ierr)
    call allocate_array(wght0,maxnxs,maxnys,'WGHT0',ierr)
    
  end subroutine allocate_mod_slout


  subroutine deallocate_mod_slout
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()
    if (allocated(grm_cell)) deallocate(grm_cell)
    if (allocated(wght0)) deallocate(wght0)
    
  end subroutine deallocate_mod_slout



end module mod_slout
