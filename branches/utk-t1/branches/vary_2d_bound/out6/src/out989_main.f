c
c ======================================================================
c
c subroutine: Swap2I
c
      SUBROUTINE Swap2I(A)
      IMPLICIT none
      INTEGER*2   A
      CHARACTER   Ktemp*1
      INTEGER     Itemp
      CHARACTER   Jtemp(2)*1
      EQUIVALENCE (Jtemp(1),Itemp)

      Itemp    = A
      Ktemp    = Jtemp(2)
      Jtemp(2) = Jtemp(1)
      Jtemp(1) = Ktemp
      A        = Itemp

      RETURN
      END
c
c ======================================================================
c
c subroutine: Swap4I
c
      SUBROUTINE  Swap4I(A)
      IMPLICIT none
      INTEGER*4   A
      CHARACTER   Ktemp*1
      INTEGER     Itemp
      CHARACTER   Jtemp(4)*1
      EQUIVALENCE (Jtemp(1),Itemp)

      Itemp  = A
      Ktemp  = Jtemp( 4 )
      Jtemp( 4 ) = Jtemp( 1 )
      Jtemp( 1 ) = Ktemp
      Ktemp  = Jtemp( 3 )
      Jtemp( 3 ) = Jtemp( 2 )
      Jtemp( 2 ) = Ktemp
      A      = Itemp

      RETURN
      END
c
c ======================================================================
c
c subroutine: Swap4R
c
      SUBROUTINE  Swap4R(A)
      IMPLICIT none
      REAL*4      A
      CHARACTER   Ktemp*1
      REAL        Itemp
      CHARACTER   Jtemp(4)*1
      EQUIVALENCE (Jtemp(1),Itemp)

      Itemp  = A
      Ktemp  = Jtemp( 4 )
      Jtemp( 4 ) = Jtemp( 1 )
      Jtemp( 1 ) = Ktemp
      Ktemp  = Jtemp( 3 )
      Jtemp( 3 ) = Jtemp( 2 )
      Jtemp( 2 ) = Ktemp
      A      = Itemp

      RETURN
      END
c
c ======================================================================
c
c subroutine: LoadImage
c
      SUBROUTINE LoadImage(opt,image,header,nxbin,nybin)
      USE MOD_OUT989
      IMPLICIT none

c...  Input:
      TYPE(type_options989) :: opt
      INTEGER nxbin,nybin
c...  Output:
      REAL*8 image(nxbin,nybin)
c...  Locals:
      INTEGER   fp,i1,iframe,ierr,ix,iy,ix1,ix2,ix3
      LOGICAL   output
      REAL      version
      CHARACTER cdum1*14,buffer*1024,file*1024
      INTEGER   iload(nxbin,nybin)
c      INTEGER   iload(1100,1100)

      
 

      TYPE(type_header) :: header

      output = .FALSE.

      fp = 98

      header%shot    = 0
      header%frame   = 0
      header%time    = 0
      header%channel = 0

      IF     (opt%imagetype.EQ.1) THEN   
c...    Simple ASCII image dump:
        fp = 99
        file = opt%fimage(1:LEN_TRIM(opt%fimage))//'.ray.img'
        OPEN(fp,FILE=file(1:LEN_TRIM(file)),
     .       FORM='FORMATTED',STATUS='OLD',ERR=98)
        READ(fp,*) version,opt%nxbin,opt%nybin
        IF (version.EQ.1.0) THEN
          DO iy = 1, opt%nybin
            DO ix = 1, opt%nxbin, 5
              READ(fp,'(1P,10E20.12)') 
     .          (image(ix1,iy),ix1=ix,MIN(ix+4,opt%nxbin))
            ENDDO
          ENDDO     
        ELSE
          CALL ER('LoadImage','Unsupported .img format',*99)
        ENDIF
        CLOSE(fp)

      ELSEIF (opt%imagetype.EQ.2) THEN   

        DO iframe = 1, 1

          IF (.TRUE.) THEN
c...        ASCII:
            file = opt%fimage(1:LEN_TRIM(opt%fimage))
            OPEN(fp,FILE=file(1:LEN_TRIM(file)),
     .           FORM='FORMATTED',STATUS='OLD',ERR=98)

            READ(fp,'(A)') buffer

            WRITE(0,*) 'buffer:',buffer(1:100)
            IF (buffer(1:8).EQ.'{version') THEN
              version = 1.1              
            ELSE
              version = 1.0
              BACKSPACE fp
            ENDIF

            READ(fp,'(A14,A8)')     cdum1,header%id         
            READ(fp,'(A14,I12)')    cdum1,header%size       
            READ(fp,'(A14,A8)')     cdum1,header%codec       
            READ(fp,'(A14,A20)')    cdum1,header%date_time  
            READ(fp,'(A14,I12)')    cdum1,header%shot       
            READ(fp,'(A14,I12)')    cdum1,header%frame       
            READ(fp,'(A14,F12.6)')  cdum1,header%time    
            READ(fp,'(A14,I12)')    cdum1,header%channel       
            READ(fp,'(A14,F12.8)')  cdum1,header%trigger    
            READ(fp,'(A14,A24)')    cdum1,header%lens       
            READ(fp,'(A14,A24)')    cdum1,header%filter     
            READ(fp,'(A14,A64)')    cdum1,header%view       
            READ(fp,'(A14,I12)')    cdum1,header%numframes  
            READ(fp,'(A14,A64)')    cdum1,header%camera     
            READ(fp,'(A14,I12)')    cdum1,header%width      
            READ(fp,'(A14,I12)')    cdum1,header%height     
            READ(fp,'(A14,I12)')    cdum1,header%depth      
            READ(fp,'(A14,I12)')    cdum1,header%orient     
            READ(fp,'(A14,I12)')    cdum1,header%taps       
            READ(fp,'(A14,I12)')    cdum1,header%color      
            READ(fp,'(A14,I12)')    cdum1,header%hBin       
            READ(fp,'(A14,I12)')    cdum1,header%left       
            READ(fp,'(A14,I12)')    cdum1,header%right      
            READ(fp,'(A14,I12)')    cdum1,header%vBin       
            READ(fp,'(A14,I12)')    cdum1,header%top        
            READ(fp,'(A14,I12)')    cdum1,header%bottom     
            READ(fp,'(A14,2I12)')   cdum1,(header%offset(i1),i1=1,2)  
            READ(fp,'(A14,2F12.8)') cdum1,(header%gain(i1),i1=1,2)    ! Gain not written properly from IDL..?
            READ(fp,'(A14,I12)')    cdum1,header%preExp     
            READ(fp,'(A14,I12)')    cdum1,header%shutter    
            READ(fp,'(A14,I12)')    cdum1,header%strobe     
            READ(fp,'(A14,F12.8)')  cdum1,header%temperature
c             Also want frame time, calibration reference, ...

c...        Some memory please:
c            IF (ALLOCATED(image)) DEALLOCATE(image)
c            ALLOCATE(image(header%width,header%height))
c...        Read image:
            IF (version.EQ.1.1) THEN
              DO iy = 1, header%height    ! Flip vertically
                DO ix = 1, header%width, 10
                  READ(fp,'(10F12.4)') 
     .              (image(ix1,iy),ix1=ix,MIN(ix+9,header%width))
                ENDDO
              ENDDO
            ELSE
              DO iy = 1, header%height    ! Flip vertically
                DO ix = 1, header%width, 10
                  READ(fp,'(10I6)') 
     .              (iload(ix1,iy),ix1=ix,MIN(ix+9,header%width))
                ENDDO
              ENDDO
            ENDIF

          ELSE
c...        Binary (disaster):
            file =  opt%fimage(1:LEN_TRIM(opt%fimage))
            OPEN(fp,FILE=file(1:LEN_TRIM(file))//'-header',
     .           FORM='UNFORMATTED',STATUS='OLD',ERR=98)

            READ(fp,ERR=97,END=96) header

            IF (.TRUE.) THEN
              CALL Swap4I(header%size)
              CALL Swap4I(header%shot)
              CALL Swap4R(header%trigger)
              CALL Swap4I(header%numFrames)
              CALL Swap2I(header%width)
              CALL Swap2I(header%height)
              CALL Swap2I(header%depth)
              CALL Swap4I(header%orient)
              CALL Swap2I(header%taps)
              CALL Swap2I(header%color)
              CALL Swap2I(header%hBin)
              CALL Swap2I(header%left)
              CALL Swap2I(header%right)
              CALL Swap2I(header%vBin)
              CALL Swap2I(header%top)
              CALL Swap2I(header%bottom)
              CALL Swap2I(header%offset(1))
              CALL Swap2I(header%offset(2))
              CALL Swap4R(header%gain(1))
              CALL Swap4R(header%gain(2))
              CALL Swap4I(header%preExp)
              CALL Swap4I(header%shutter)
              CALL Swap4I(header%strobe)
              CALL Swap4R(header%temperature)
            ENDIF
c...        Some memory please:
            STOP '989 : HERE A'
c            IF (ALLOCATED(image)) DEALLOCATE(image)
c            ALLOCATE(image(header%width,header%height))
c...        Read image:
            READ(fp,ERR=97,END=96,IOSTAT=ierr) 
     .        ((iload(ix,iy),ix=1,header%width),
     .                       iy=1,header%height)

          ENDIF

          opt%nxbin = header%width
          opt%nybin = header%height
          IF (version.EQ.1.0) 
     .      image(1:opt%nxbin,1:opt%nybin) = DFLOAT(iload(1:opt%nxbin,
     .                                                    1:opt%nybin))

        ENDDO

        IF (output) THEN
          WRITE(0,*) 'id         :',header%id
          WRITE(0,*) 'size       :',header%size      
          WRITE(0,*) 'codec      :',header%codec
          WRITE(0,*) 'date_time  :',header%date_time
          WRITE(0,*) 'shot       :',header%shot      
          WRITE(0,*) 'trigger    :',header%trigger   
          WRITE(0,*) 'lens       :',header%lens
          WRITE(0,*) 'filter     :',header%filter
          WRITE(0,*) 'view       :',header%view
          WRITE(0,*) 'numframes  :',header%numFrames 
          WRITE(0,*) 'camera     :',header%camera
          WRITE(0,*) 'width      :',header%width     
          WRITE(0,*) 'height     :',header%height    
          WRITE(0,*) 'depth      :',header%depth     
          WRITE(0,*) 'orient     :',header%orient
          WRITE(0,*) 'taps       :',header%taps      
          WRITE(0,*) 'color      :',header%color     
          WRITE(0,*) 'hBin       :',header%hBin      
          WRITE(0,*) 'left       :',header%left      
          WRITE(0,*) 'right      :',header%right      
          WRITE(0,*) 'vBin       :',header%vBin      
          WRITE(0,*) 'top        :',header%top       
          WRITE(0,*) 'bottom     :',header%bottom    
          WRITE(0,*) 'offset(1)  :',header%offset(1) 
          WRITE(0,*) 'offset(2)  :',header%offset(2) 
          WRITE(0,*) 'gain(1)    :',header%gain(1)   
          WRITE(0,*) 'gain(2)    :',header%gain(2)   
          WRITE(0,*) 'preexp     :',header%preExp  
          WRITE(0,*) 'shutter    :',header%shutter 
          WRITE(0,*) 'strobe     :',header%strobe    
          WRITE(0,*) 'tempreature:',header%temperature 
        ENDIF

      ELSE
        CALL ER('LoadImage','Invalid image type',*99)
      ENDIF

      CLOSE(fp)


      RETURN
 96   CALL ER('989','UNEXPECTED END OF IMAGE FILE',*99)
 97   CALL ER('989','ERROR READING IMAGE FILE',*99)
 98   CALL ER('989','IMAGE FILE NOT FOUND',*99)
 99   WRITE(0,*) 'FILE  : '//file(1:LEN_TRIM(file))
      WRITE(0,*) 'IOSTAT: ',ierr
      STOP
      END
c
c ======================================================================
c
      SUBROUTINE LoadOptions989(opt)
      USE MOD_OUT989
      IMPLICIT none

      TYPE(type_options989) :: opt

      INTEGER   fp1,idum1
      REAL      rdum1
      CHARACTER buffer*1024,cdum1*1024

      opt%out_suffix = 'none'

      IF (.TRUE.) THEN
c...    Load plot options:
        DO WHILE (.TRUE.)
          READ(5,'(A1024)') buffer
          IF (buffer(1:LEN_TRIM(buffer)).EQ.'$') EXIT

          WRITE(0,*) 'IMG>>'//buffer(1:30)//'<<'

          IF (buffer(8:12).EQ.'Image'.OR.buffer(8:12).EQ.'IMAGE'.OR.
     .        buffer(8:12).EQ.'image') THEN

            READ(buffer,*) cdum1,opt%imagetype,idum1,rdum1,opt%fimage

            WRITE(0,*) 'IMG>>'//opt%fimage(1:LEN_TRIM(opt%fimage))//'<<'

            IF (opt%imagetype.NE.0) THEN
              opt%nimg = opt%nimg + 1
              opt%img_type (opt%nimg) = opt%imagetype
              opt%img_scale(opt%nimg) = rdum1
              opt%img_fname(opt%nimg) = opt%fimage
            ENDIF

c          ELSE
c            BACKSPACE 5
c            EXIT
          ENDIF

c        ENDDO


c          READ(5,'(A1024)') buffer
c          WRITE(0,*) 'MAP>>'//buffer(1:30)//'<<'
          IF (buffer(8:10).EQ.'Map'.OR.buffer(8:10).EQ.'MAP'.OR.
     .        buffer(8:10).EQ.'map') THEN

            READ(buffer,*) cdum1,idum1,idum1,rdum1,opt%fmap

            WRITE(0,*) 'MAP>>'//opt%fmap(1:LEN_TRIM(opt%fmap))//'<<'

          ENDIF

          IF (buffer(8:13).EQ.'Output'.OR.buffer(8:13).EQ.'OUTPUT'.OR.
     .        buffer(8:13).EQ.'output') THEN

            READ(buffer,*) cdum1,opt%out_opt

            READ(5,*) cdum1,opt%out_suffix

            WRITE(0,*) 'OUT>>'//
     .        opt%out_suffix(1:LEN_TRIM(opt%out_suffix))//'<<'
          ENDIF

        ENDDO

      ENDIF

      IF (opt%nimg.EQ.0) THEN
        WRITE(0,*) 'Inversion image file not specified'
        BACKSPACE 5
c        CALL ER('989','Image file not specified',*99)
      ENDIF



      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: Main989
c
c Fun with Ax=b. 
c
c
c
      SUBROUTINE CallMain985
      USE mod_out985
      IMPLICIT none

      STOP 'Calling 985 from 989 not a go at the moment...'

c...  Input:

c...  Locals:

c      TYPE(type_options985) :: opt985     
c      INTEGER nobj,MAX3D
c      TYPE(type_3D_object), ALLOCATABLE :: obj(:)
c      INTEGER npixel,MAXPIXEL
c      TYPE(type_view), ALLOCATABLE :: pixel(:)
c      TYPE(type_options989) :: opt
c      REAL*8, ALLOCATABLE, TARGET :: image(:,:),image2(:,:,:)

c      TYPE(type_header) :: header
c      TYPE(type_header),ALLOCATABLE :: header2(:)

c      REAL qmax1, qmax2,qavg1,qavg2

c      INTEGER, ALLOCATABLE :: npts(:)
c      REAL,    ALLOCATABLE :: tpts(:),xpts(:,:),ypts(:,:),qpts(:,:)



c          WRITE(0,*) '  GENERATING MAP AND IMAGE FROM 985 - CAUTION!'

c...      Calculate on the fly:
c          MAX3D = 1100000
c          MAXPIXEL=1100*1100
c          ALLOCATE(pixel(MAXPIXEL))
c          ALLOCATE(obj(MAX3D))
c          CALL ALLOC_CHORD(11000)  ! *TEMP* Just for viewing!

c          IF (opt%imagetype.NE.0) THEN
c            opt985%img_opt = 1
c            opt985%img_nxibn = opt%nxbin
c            opt985%img_nybin = opt%nybin
c            opt985%img_image => image
c          ELSE
c            opt985%img_opt = 0
c          ENDIF

c          STOP 'ASFALKJFAS'  ! LEFT OFF -- PARHAPS NEED TO DECLARI OBJ() OUTSIDE Main985?

c          CALL Main985(1,opt985,MAXPIXEL,npixel,pixel,MAX3D,nobj,obj,
c     .                 image)

c          WRITE(0,*) '  985 NX,NY:',opt985%nxbin,opt985%nybin

c...      Lame plot setup, lots of assumptions:
c          ninv = opt985%n
c          WRITE(0,*) '  NINV:',ninv
c          ALLOCATE(npts(ninv))
c          ALLOCATE(xpts(4,ninv))
c          ALLOCATE(ypts(4,ninv))
c          ALLOCATE(tpts(ninv))
c          ALLOCATE(qpts(ninv,2))
c          DO i1 = 1, ninv
c            npts(i1) = 0
c            IF (obj(i1)%nsur.EQ.4) npts(i1) = 4  ! Toroidally continuous   ! ***LAME/WEAK***
c            IF (obj(i1)%nsur.EQ.5) npts(i1) = 3  ! Triangles
c            IF (obj(i1)%nsur.EQ.6) npts(i1) = 4  ! Quadrelateral           
c            IF (npts(i1).EQ.0) THEN
c              WRITE(0,*) '  989: PROBLEM OBJECT',i1,n
c              STOP 'HALTING CODE'
c            ENDIF
c            DO i2 = 1, npts(i1) 
c              xpts(i2,i1) = SNGL(obj(i1)%v(1,i2))
c              ypts(i2,i1) = SNGL(obj(i1)%v(2,i2))
c            ENDDO
c            tpts(i1) = obj(i1)%path
c            qpts(i1,1) = obj(i1)%quantity(1)
c          ENDDO 

c...      Clear memory (put into subroutine):
c          IF (ALLOCATED(obj))   DEALLOCATE(obj)
c          IF (ALLOCATED(pixel)) DEALLOCATE(pixel)
c          CALL DEALLOC_CHORD  ! *TEMP* until this mess gets sorted out
c          IF (ALLOCATED(vtx))   DEALLOCATE(vtx)
c          IF (ALLOCATED(srf))   DEALLOCATE(srf)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: Main989
c
c Fun with Ax=b. 
c
c
c
      SUBROUTINE Main989(iopt)
      USE mod_out985
      USE mod_out989
      IMPLICIT none

c...  Input:
      INTEGER iopt

      REAL Clock2

c...  Locals:
      INTEGER   i1,i2,fp,codec,ix,ix1,iy,m,n,ierr,ninv,volcnt,pixcnt,
     .          idum1,ios,count,ib,iimg
      LOGICAL   ascii,cont
      CHARACTER dummy*1024,file*1024,afile*1024
      REAL      rtime,ttime,rdum1,version
      REAL*8    colsum,rowsum
      REAL*8,   ALLOCATABLE :: A(:,:),x(:),b(:)

      TYPE(type_options985) :: opt985     
      INTEGER nobj,MAX3D
      TYPE(type_3D_object), ALLOCATABLE :: obj(:)
      INTEGER npixel,MAXPIXEL
      TYPE(type_view), ALLOCATABLE :: pixel(:)
      TYPE(type_options989) :: opt
      REAL*8, ALLOCATABLE, TARGET :: image(:,:),image2(:,:,:)

      TYPE(type_header) :: header
      TYPE(type_header),ALLOCATABLE :: header2(:)

      REAL qmax1, qmax2,qavg1,qavg2

      INTEGER, ALLOCATABLE :: npts(:)
      REAL,    ALLOCATABLE :: tpts(:),xpts(:,:),ypts(:,:),qpts(:,:)

      INTEGER toy,status
      EXTERNAL toy 
      REAL sigma



      IF (.FALSE.) THEN
        status = Toy(%VAL(1),%VAL(10),
     .    '/home/slisgo/maxent/files/m-lmr-0000x.channel_a.ray.A',
     .    %VAL(1),%VAL(1),%VAL(0),%VAL(0.01),%VAL(sigma))

c            status = Toy(1,10,
c     .        '/home/slisgo/maxent/files/m-lmr-0000x.channel_a.ray.A',
c     .        1,1,0,0.01,sigma)
        WRITE(0,*) 'STATUS=',status

        STOP 'SHIT'
      ENDIF



      ALLOCATE(image (1100,1100))
      ALLOCATE(image2(1100,1100,MAXNIMAGE))

      ALLOCATE(header2(MAXNIMAGE))


      opt985%load = 1   ! *** do a better job of this...

      ninv = 0

      WRITE(0,*) 'LOADING 989 OPTIONS'

c...  Load inversion options:
      opt%nimg = 0
      CALL LoadOptions989(opt)




      count = 0
      cont = .TRUE.
      DO WHILE (cont) 
        count = count + 1
        WRITE(0,*) 'LOOP PASS',count



        IF (iopt.EQ.2) THEN

          IF (opt%nimg.GT.1) STOP 'PROBLEMS WITH MULTIPLE IMAGES'

          CALL CallMain985

        ENDIF



c     iopt = 1  -normal
c     iopt = 2  -do ray trace from here
c     iopt = 3  -previous inversion
c



        IF (iopt.EQ.1.OR.iopt.EQ.2) THEN

c...      Analyse geometry map from detector onto toroidally symmetric mesh:
c         (the proverbial 'A'):
          fp = 99
          file = opt%fmap(1:LEN_TRIM(opt%fmap))//'.ray.map'
          afile = file


          ascii = .TRUE.
          OPEN(fp,FILE=file(1:LEN_TRIM(file)),
     .         FORM='FORMATTED',STATUS='OLD',ERR=95)
          READ(fp,'(A1024)') dummy
          WRITE(0,*) '  dummy:',dummy(1:10)
          IF (dummy(1:1).EQ.'*') THEN
            BACKSPACE fp
          ELSE
c...        Try binary format:
            CLOSE(fp)
            ascii = .FALSE.
            OPEN(fp,FILE=file(1:LEN_TRIM(file)),
     .           FORM='UNFORMATTED',STATUS='OLD',ERR=95)   
            READ(fp) idum1
            IF (idum1.NE.888888888) 
     .        CALL ER('Main989','Unable to open inversion map file',*99)
          ENDIF

          IF (ascii) THEN
            dummy(1:1) = '*'
            DO WHILE (dummy(1:1).EQ.'*')
              READ(fp,'(A1024)') dummy
            ENDDO
            BACKSPACE fp
            READ(fp,'(A14,F5.1)') dummy,version
            IF (version.EQ.1.0) THEN
             READ(fp,'(A14,I7)') dummy,codec
             READ(fp,'(A14,I7)') dummy,m
             READ(fp,'(A14,I7)') dummy,n
             READ(fp,'(A14,I7)') dummy,opt985%ndet
             READ(fp,'(A14,99I7)') dummy,opt985%det_nxbin(1:opt985%ndet)
             READ(fp,'(A14,99I7)') dummy,opt985%det_nybin(1:opt985%ndet)
            ELSE
              CALL ER('Main989','Unrecognized ASCII .map version',*99)
            ENDIF
          ELSE
c...        Binary format:
            READ(fp) version
            IF (version.EQ.1.0) THEN
              READ(fp) codec
              READ(fp) m,n,opt985%ndet,opt985%det_nxbin(1:opt985%ndet),
     .                                 opt985%det_nybin(1:opt985%ndet)
            ELSE
              CALL ER('Main989','Unrecognized binary .map version',*99)
            ENDIF
          ENDIF

          opt985%nxbin = opt985%det_nxbin(1)   ! For backwards compatibility... (for now)
          opt985%nybin = opt985%det_nybin(1)

          WRITE(0,*) '  MAP NX,NY:',opt985%det_nxbin(1:opt985%ndet)
          WRITE(0,*) '           :',opt985%det_nybin(1:opt985%ndet)

          IF (ascii.AND.codec.EQ.1) THEN
            ALLOCATE(A(m,n))
            DO ix = 1, m
              DO i1 = 1, n, 5  ! Better way to write this?
                READ(fp,'(10D20.12)') (A(ix,iy),iy=i1,MIN(i1+4,n))
              ENDDO
              IF (n.LE.5) 
     .          WRITE(0,'(A,I4,10D20.12)') '  A:',ix,(A(ix,iy),iy=1,n)
            ENDDO        
          ELSEIF (codec.EQ.2) THEN
            WRITE(0,*) '  CODEC.EQ.2 IN MAIN989'
          ELSE
            CALL ER('Main989','Bad codec in .map file',*99)
          ENDIF

c...      Get poloidal cross-section for inversion mesh from end of 
c         .map file:
          IF (iopt.EQ.1) THEN
            IF (ascii) THEN
              DO WHILE(.TRUE.) 
                READ(fp,'(A1024)',END=96) dummy
                IF (dummy(1:10).EQ.'* Inversio') EXIT
              ENDDO
              READ(fp,*) ninv
            ELSE
c...          Advance the file to the start of the 
              DO WHILE(.TRUE.) 
                READ(fp,END=96) idum1
                IF (idum1.EQ.999999999) EXIT
              ENDDO
              READ(fp) ninv
            ENDIF
            WRITE(0,*) '  NINV:',ninv
       
            ALLOCATE(npts(ninv))
            ALLOCATE(xpts(4,ninv))
            ALLOCATE(ypts(4,ninv))
            ALLOCATE(tpts(ninv))
            ALLOCATE(qpts(ninv,2))
            IF (ascii) THEN
              DO i1 = 1, n
                READ(fp,*) idum1,npts(i1),qpts(i1,1),tpts(i1),
     .                     (xpts(i2,i1),ypts(i2,i1),i2=1,npts(i1))
              ENDDO 
            ELSE
              DO i1 = 1, n
                READ(fp) idum1,npts(i1),qpts(i1,1),tpts(i1),
     .                   (xpts(i2,i1),ypts(i2,i1),i2=1,npts(i1))
              ENDDO 
            ENDIF

c            qpts(1,1) = -998.0             ! *** TURN OFF CHECK PLOTS -- LAME ***
  
          ENDIF
          CLOSE(fp)


          IF (codec.EQ.1) THEN
c...        Check if a row of 'A' is everywhere zero:
            pixcnt = m
            volcnt = n
            DO ix = 1, m
              rowsum = 0.0D0
              DO iy = 1, n
                rowsum = rowsum + A(ix,iy)
              ENDDO
              IF (rowsum.EQ.0.0D0) pixcnt = pixcnt - 1
            ENDDO
            DO iy = 1, n
              colsum = 0.0D0
              DO ix = 1, m
                colsum = colsum + A(ix,iy)
              ENDDO
              IF (colsum.EQ.0.0D0) volcnt = volcnt - 1
            ENDDO
            WRITE(0,'(2X,F8.1,A,I7,A,I7,A)') 
     .      REAL(pixcnt)/REAL(m)*100.0,'% OF DETECTOR VIEWS '//
     .      'SAMPLE THE INVERSION MESH  (',pixcnt,' of',m,')'
            WRITE(0,'(2X,F8.1,A,I7,A,I7,A)') 
     .        REAL(volcnt)/REAL(n)*100.0,'% OF INVERSION MESH '//
     .        'SAMPLED BY LOS INTEGRATION (',volcnt,' of',n,')'
c            DO ix = 1, m
c              A(ix,1:4) = A(ix,3:6)
c            ENDDO
c            n = n - 2  
          ENDIF



c...      Aquire image for inversion, either simulated for from a camera:
c         (the proverbial 'b'):


          IF (opt%nimg.GE.1) THEN

            npixel = 0

            DO iimg = 1, opt%nimg
              IF (opt%img_type(iimg).NE.0) THEN
  
                WRITE(0,*) '  LOADING IMAGE',iimg

                opt%imagetype = opt%img_type (iimg)   ! Backwards compatibilty... (for now)
                opt%fimage    = opt%img_fname(iimg)

                CALL LoadImage(opt,image2(1,1,iimg),header2(iimg),
     .                         1100,1100)

                header = header2(iimg)
                opt%img_nxbin(iimg) = opt%nxbin
                opt%img_nybin(iimg) = opt%nybin
                WRITE(0,*) 'BINS1:',iimg,
     .            opt%img_nxbin(iimg),opt%img_nybin(iimg)

c...            Scale image intensity:
                IF (opt%img_scale(iimg).EQ.0.0D0) THEN
                  WRITE(0,*) '***WARNING*** IMAGE SCALE = 0.0 DETECTED'
                  WRITE(0,*) '              SETTING TO 1.0'
                  opt%img_scale(iimg) = 1.0D0
                ENDIF
                DO ix = 1, opt%img_nxbin(iimg)
                  DO iy = 1, opt%img_nybin(iimg)
                    image2(ix,iy,iimg) = image2(ix,iy,iimg) * 
     .                                   opt%img_scale(iimg)
                  ENDDO
                ENDDO

c...            Binning:
                IF (opt%img_nxbin(iimg).NE.opt985%det_nxbin(iimg).OR.
     .              opt%img_nybin(iimg).NE.opt985%det_nybin(iimg)) THEN
                  WRITE(0,*) '    BINNING IMAGE IN 989'
                  CALL BinImage(opt%img_nxbin(iimg),
     .                          opt%img_nybin(iimg),image2(1,1,iimg),
     .                          opt985%det_nxbin(iimg),
     .                          opt985%det_nybin(iimg))     
                  WRITE(0,*) '    DONE'
                ENDIF

                npixel = npixel +opt%img_nxbin(iimg)*opt%img_nybin(iimg)

              ENDIF
            ENDDO
          ELSE
c...        Load from file:
            WRITE(0,*) '  LOADING IMAGE'
            IF (opt%imagetype.NE.0) CALL LoadImage(opt,image,header,
     .                                             1100,1100)
c...        Perform required binning, if any (might be complex binning eventually), and/or 
c           select pixels (viewing chords) to be used in the inversion:
            IF (opt%nxbin.NE.opt985%nxbin.OR.
     .          opt%nybin.NE.opt985%nybin) THEN
              WRITE(0,*) '  BINNING IN 989!'
              CALL BinImage(opt%nxbin,opt%nybin,image,
     .                      opt985%nxbin,opt985%nybin)     
              WRITE(0,*) '  DONE'
            ENDIF
            npixel = opt%nxbin * opt%nybin
          ENDIF


          ALLOCATE(b(npixel))

c...      Assign 'b' vector:
          IF (opt%nimg.GE.1) THEN
            ib = 0
            DO iimg = 1, opt%nimg
              DO iy = 1, opt%img_nybin(iimg)
                DO ix = 1, opt%img_nxbin(iimg)
                  ib = ib + 1
                  b(ib) = image2(ix,iy,iimg)
                ENDDO
              ENDDO
            ENDDO

c..         For plots (for now...)
            iimg = 1
            opt%nxbin =  opt%img_nxbin(iimg)
            opt%nybin =  opt%img_nybin(iimg)
            DO ix = 1, 1100
              DO iy = 1, 1100
                image(ix,iy) = image2(ix,iy,iimg)
              ENDDO
            ENDDO

          ELSE
            DO iy = 1, opt%nybin
              b((iy-1)*opt%nxbin+1:iy*opt%nxbin) = image(1:opt%nxbin,iy)
            ENDDO
          ENDIF
 

c...      Now solve for (yes, the proverbial) 'x':
          ALLOCATE(x(n))

          rtime = Clock2()

          IF (.TRUE.) THEN
            WRITE(0,*) '  M:',m,npixel
            WRITE(0,*) '  N:',n

            IF (m.NE.npixel) 
     .        CALL ER('Plot989','m.NE.npixel',*99)

            IF (n.GT.m) 
     .        CALL ER('Plot989','N>M no allowed, for some reason',*99)  ! This may be a bug, or numerical rule...

            CALL Invert_LSQR(m,n,afile,x,b,0,0.15D0)  ! 0.15)

            WRITE(0,*) '  BACK FROM INVERSION ROUTINE'

c            DO i1 = 1, n
c             IF (x(i1).LT.0.0D0) WRITE(6,*) ' -VE X:',i1,x(i1)
c            ENDDO

            IF (ALLOCATED(qpts)) qpts(1:n,2) = SNGL(x(1:n))
          ENDIF
 
          ttime = Clock2() - rtime

          WRITE(0,'(A,F10.2,A)') '  TOTAL INVERSION TIME:',ttime,' s'
          WRITE(0,'(A,F10.3,A)') '  AVERAGE TIME PER RANK:',
     .                           ttime/REAL(n*m),' s'
          WRITE(0,'(A,F10.1,A)') '  TIME FOR 1E12 RANK    :',
     .                           ttime/REAL(n*m)*1E+12/3600.0,' h'


          IF (ALLOCATED(A)) DEALLOCATE(A)
          IF (ALLOCATED(x)) DEALLOCATE(x)
          IF (ALLOCATED(b)) DEALLOCATE(b)


        ELSEIF (iopt.EQ.3) THEN

c...      Load a previous reconstruction:
          fp = 99
          file = './cal-rdb-17-0302-1000x1000-1.ray.cgm'
c          file = './cal-rda-13-0259-1000x1000-0.5.ray.cgm'
          OPEN(fp,FILE=file(1:LEN_TRIM(file)),
     .         FORM='FORMATTED',STATUS='OLD',ERR=95)   
          READ(fp,*) ninv

          ALLOCATE(npts(ninv))
          ALLOCATE(xpts(4,ninv))
          ALLOCATE(ypts(4,ninv))
          ALLOCATE(tpts(ninv))
          ALLOCATE(qpts(0:ninv,2))

          tpts(1:ninv) = 1.0

          DO i1 = 1, ninv
            READ(fp,'(I7,E22.12,I5,10(2F11.6))')  
     .        idum1,qpts(i1,2),
     .        npts(i1),(xpts(i2,i1),ypts(i2,i1),i2=1,npts(i1))
          ENDDO     
          CLOSE(fp)

          IF (.TRUE.) THEN
            file = './cal-rda-17-0302-1000x1000-1.ray.cgm'
            OPEN(fp,FILE=file(1:LEN_TRIM(file)),
     .         FORM='FORMATTED',STATUS='OLD',ERR=95)   
            READ(fp,*) ninv
            DO i1 = 1, ninv
              READ(fp,'(I7,E22.12)')
     .          idum1,qpts(i1,1)
            ENDDO     
            CLOSE(fp)
          ELSE
            qpts(1:ninv,1) = qpts(1:ninv,2)
          ENDIF

          IF (.TRUE.) THEN
            qavg1 = 0.0
            qavg2 = 0.0
            qmax1 = -1.0E+30
            DO i1 = 1, ninv
              qavg1 = qavg1 + qpts(i1,1)
              qavg2 = qavg2 + qpts(i1,2)
              qmax1 = MAX(qmax1,qpts(i1,1))
              qmax2 = MAX(qmax2,qpts(i2,1))
            ENDDO
            qavg1 = qavg1 / REAL(ninv)
            qavg2 = qavg2 / REAL(ninv)

            WRITE(0,*) 'AVG:',qavg1,qavg2

            IF (.FALSE.) THEN
              DO i1 = 1, ninv
                IF (qpts(i1,1).GT.500.0.AND.
     .              qpts(i1,2).GT.0.1*qmax2) THEN
c                IF (qpts(i1,1).GT.0.01*qmax1.AND.
c     .              qpts(i1,2).GT.0.10*qmax2.AND. 
c     .              qpts(i1,1).GT.1.00*qavg1.AND.
c     .              qpts(i1,2).GT.0.15*qavg2) THEN
                  qpts(i1,2) = qpts(i1,1) / (2.8 * qpts(i1,2))
                ELSE
                  qpts(i1,2) = 0.0D0
                ENDIF
              ENDDO
            ENDIF
          ENDIF

          IF (opt%nimg.GE.1) THEN

            DO iimg = 1, opt%nimg
              IF (opt%img_type(iimg).NE.0) THEN
  
                WRITE(0,*) '  LOADING IMAGE',iimg

                opt%imagetype = opt%img_type (iimg)   ! Backwards compatibilty... (for now)
                opt%fimage    = opt%img_fname(iimg)

                CALL LoadImage(opt,image2(1,1,iimg),header2(iimg),
     .                         1100,1100)

                header = header2(iimg)
                opt%img_nxbin(iimg) = opt%nxbin
                opt%img_nybin(iimg) = opt%nybin
                WRITE(0,*) 'BINS1:',iimg,
     .            opt%img_nxbin(iimg),opt%img_nybin(iimg)

c...            Scale image intensity:
                IF (opt%img_scale(iimg).EQ.0.0D0) THEN
                  WRITE(0,*) '***WARNING*** IMAGE SCALE = 0.0 DETECTED'
                  WRITE(0,*) '              SETTING TO 1.0'
                  opt%img_scale(iimg) = 1.0D0
                ENDIF
                DO ix = 1, opt%img_nxbin(iimg)
                  DO iy = 1, opt%img_nybin(iimg)
                    image2(ix,iy,iimg) = image2(ix,iy,iimg) * 
     .                                   opt%img_scale(iimg)
                  ENDDO
                ENDDO

                

c...            Binning:
                IF (.TRUE.) THEN
                  WRITE(0,*) '    BINNING IMAGE IN 989 PRES.'
                  CALL BinImage(opt%img_nxbin(iimg),
     .                          opt%img_nybin(iimg),image2(1,1,iimg),
     .                          100,
     .                          100)
                  WRITE(0,*) '    DONE'
                ENDIF

              ENDIF
            ENDDO

c..         For plots (for now...)
            iimg = 1
            opt%nxbin =  opt%img_nxbin(iimg)
            opt%nybin =  opt%img_nybin(iimg)
            DO ix = 1, 1100
              DO iy = 1, 1100
                image(ix,iy) = image2(ix,iy,iimg)
              ENDDO
            ENDDO

          ELSE
            opt%nxbin = 1
            opt%nybin = 1
          ENDIF
 
        ENDIF   ! End of IOPT selection

c...    Plots:
        CALL Output989(opt,1100,1100,image,image2,header,header2,
     .                 ninv,npts,xpts,ypts,tpts,qpts)



        IF (ALLOCATED(npts)) DEALLOCATE(npts)
        IF (ALLOCATED(xpts)) DEALLOCATE(xpts)
        IF (ALLOCATED(ypts)) DEALLOCATE(ypts)
        IF (ALLOCATED(tpts)) DEALLOCATE(tpts)
        IF (ALLOCATED(qpts)) DEALLOCATE(qpts)


c...    Loop exit conditions:
        IF     (count.EQ.1) THEN 
          opt985%load = 0
          cont = .FALSE.
        ELSEIF (.FALSE.) THEN
          cont = .FALSE.
        ELSE
          cont = .FALSE.
        ENDIF


      ENDDO


      IF (ALLOCATED(image) ) DEALLOCATE(image)
      IF (ALLOCATED(image2)) DEALLOCATE(image2)
      IF (ALLOCATED(header2)) DEALLOCATE(header2)




      WRITE(0,*) 'DONE INVERSION PROCESSING'



      RETURN
 95   CALL ER('Main989','Unable to open file',*99)
 96   CALL ER('Main989','Unexpected EOF when reading .map',*99)
 97   CALL ER('Main989','Unable to create image dump',*99)
 99   WRITE(0,*) 'FILE:',file(1:LEN_TRIM(file))
      STOP
      END
