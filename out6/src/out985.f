
c
c ======================================================================
c
c subroutine: Plot985
c
c 3D analysis
c
      SUBROUTINE Plot985(job,graph,ref,title,iopt,
     .                   xxmin,xxmax,yymin,yymax,ft,fp,zadj,
     .                   ismoth,ignors,itec,avs,navs)
c      USE mod_out985
c      USE mod_out985_variables
      use mod_params
      use mod_slcom
      IMPLICIT none
c     INCLUDE 'params'
c     INCLUDE 'slcom'

c...  Input:
      INTEGER IOPT,ISMOTH,IGNORS(MAXNGS),ITEC,NAVS
      REAL    XXMIN,XXMAX,YYMIN,YYMAX,FP,FT,zadj,AVS(0:100)
      CHARACTER TITLE*(*),JOB*72,GRAPH*80
      CHARACTER*36 REF

c      TYPE(type_options985) :: opt     

c      INTEGER nobj,MAX3D
c      TYPE(type_3D_object), ALLOCATABLE :: obj(:)

c      INTEGER npixel,MAXPIXEL
c      TYPE(type_view), ALLOCATABLE :: pixel(:)

c      REAL*8 image(1000,1000)

c      WRITE(0,*) 'HERE IN 985'
c      WRITE(0,*) '  ALLOCATING OBJECTS'
c      MAX3D = 200000
c      ALLOCATE(obj(MAX3D))
c      CALL ALLOC_SURFACE(-1,MP_INITIALIZE)
c      WRITE(0,*) '  ALLOCATING PIXELS'
c      MAXPIXEL=1000*1000
c      ALLOCATE(pixel(MAXPIXEL))
c      CALL ALLOC_CHORD(10000)  ! Just for viewing! (make smaller!)
c      opt%load = 1
c      grd_ntorseg = eirntorseg

      CALL Main985(iopt,title)
c      CALL Main985(iopt,opt,MAXPIXEL,npixel,pixel,MAX3D,nobj,obj,image)



c...  Put into a subroutine:
cc      CALL DEALLOC_OBJ
cc      CALL DEALLOC_PIXEL
c      IF (ALLOCATED(obj))   DEALLOCATE(obj)
c      IF (ALLOCATED(pixel)) DEALLOCATE(pixel)
c      CALL DEALLOC_CHORD
c      IF (ALLOCATED(vtx))   DEALLOCATE(vtx)
c      IF (ALLOCATED(srf))   DEALLOCATE(srf)

      RETURN
99    STOP
      END
