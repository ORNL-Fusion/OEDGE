c
c ======================================================================
c
      SUBROUTINE LoadOptions988(opt)
      USE MOD_OUT988
      IMPLICIT none

      TYPE(type_options988) :: opt

      INTEGER   fp,idum1
      REAL      rdum1
      CHARACTER buffer*1024,cdum1*1024

      fp = 5

      IF (.TRUE.) THEN
c...    Load plot options:
c        READ(fp,'(A1024)') buffer
c        WRITE(0,*) '>>'//buffer(1:100)//'<<'
c        IF (buffer(8:12).EQ.'Image'.OR.buffer(8:12).EQ.'IMAGE'.OR.
c     .      buffer(8:12).EQ.'image') THEN

c          READ(buffer,*) cdum1,opt%imagetype,idum1,rdum1,opt%fimage

c          WRITE(0,*) '>>'//opt%fimage(1:LEN_TRIM(opt%fimage))//'<<'

c        ELSE
c          WRITE(0,*) 'Inversion image file not specified'
c          BACKSPACE 5
c        ENDIF

c        READ(fp,'(A1024)') buffer
c        WRITE(0,*) '>>'//buffer(1:100)//'<<'
c        IF (buffer(8:10).EQ.'Map'.OR.buffer(8:10).EQ.'MAP'.OR.
c     .      buffer(8:10).EQ.'map') THEN

c          READ(buffer,*) cdum1,idum1,idum1,rdum1,opt%fmap

c          WRITE(0,*) '>>'//opt%fmap(1:LEN_TRIM(opt%fmap))//'<<'
c        ELSE
c          WRITE(0,*) 'Inversion map file not specified'
c          BACKSPACE 5
c        ENDIF
      ENDIF



      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE Output988(opt)
      USE MOD_OUT988
      IMPLICIT none

      TYPE(type_options988) :: opt

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'comgra'
      INCLUDE 'colours'
      include 'printopt'

      INCLUDE 'slcom'
      INCLUDE 'slout'

      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost


      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: Main988
c
c Line shapes.
c
c
c
      SUBROUTINE Main988(iopt)
      USE MOD_OUT985
      USE MOD_OUT988
      IMPLICIT none

c...  Input:
      INTEGER iopt


      TYPE(type_options988) :: opt


      TYPE(type_options985) :: opt985     
      INTEGER nobj,MAX3D
      TYPE(type_3D_object), ALLOCATABLE :: obj(:)
      INTEGER npixel,MAXPIXEL
      TYPE(type_view), ALLOCATABLE :: pixel(:)


c...  Load inversion options:
      CALL LoadOptions988(opt)


      IF (iopt.EQ.2) THEN
        WRITE(0,*) 'GENERATING MAP AND IMAGE FROM 985'
c...    
        MAX3D = 1000000
        MAXPIXEL=1000*1000
        ALLOCATE(pixel(MAXPIXEL))
        ALLOCATE(obj(MAX3D))
        CALL ALLOC_CHORD(10000)  ! *TEMP* Just for viewing!

        CALL Main985(1,opt985,MAXPIXEL,npixel,pixel,MAX3D,nobj,obj)

        WRITE(0,*) 'DATA:',opt985%nxbin,opt985%nybin

c...    Clear memory:
        IF (ALLOCATED(obj)) DEALLOCATE(obj)
        IF (ALLOCATED(pixel)) DEALLOCATE(pixel)
        CALL DEALLOC_CHORD  ! *TEMP* until this mess gets sorted out
      ENDIF


c...  Need to be able to handle large numbers of spectra, thinking of the
c     future -- solution, use disk storage (as for .map in 989): 




c...  Plots:
      CALL Output988(opt)



      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: Plot985
c
c 3D analysis
c
      SUBROUTINE Plot988(job,graph,ref,title,iopt,
     .                   xxmin,xxmax,yymin,yymax,ft,fp,zadj,
     .                   ismoth,ignors,itec,avs,navs)
      USE MOD_OUT988
      IMPLICIT none
      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'comgra'
      INCLUDE 'colours'
      include 'printopt'

      INCLUDE 'slcom'
      INCLUDE 'slout'

      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost

c...  Input:
      INTEGER IOPT,ISMOTH,IGNORS(MAXNGS),ITEC,NAVS
      REAL    XXMIN,XXMAX,YYMIN,YYMAX,FP,FT,zadj,AVS(0:100)
      CHARACTER TITLE*(*),JOB*72,GRAPH*80
      CHARACTER*36 REF

      CALL Main988(iopt)

      RETURN
 99   STOP
      END
