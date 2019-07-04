c     -*-Fortran-*-
c
c ======================================================================
c
c subroutine: Output989
c
      SUBROUTINE Output989(opt,nxbin_space,nybin_space,image_space,
     .                     image2,header,header2,
     .                     ninv,npts,xpts,ypts,tpts,qpts)
      USE mod_out989
      USE mod_out985_plots
      USE mod_interface
      IMPLICIT none

      TYPE(type_options989) :: opt

      INTEGER nxbin_space,nybin_space
      REAL*8 image2(nxbin_space,nybin_space,MAXNIMAGE)
      REAL*8 image_space(nxbin_space,nybin_space)
      TYPE(type_header) :: header, header2(MAXNIMAGE)

      INTEGER ninv,npts(ninv),ncel,ndat
      REAL    xpts(4,ninv),ypts(4,ninv),tpts(ninv),qpts(ninv,2),
     .        xdat(ninv),ydat(ninv)

      LOGICAL, PARAMETER :: plot_option = .TRUE.

      REAL, PARAMETER :: HI = 1.0E+30
      INTEGER, PARAMETER :: ncols = 0


c      INCLUDE 'params'
c      INCLUDE 'cgeom'
c      INCLUDE 'comtor'
c      INCLUDE 'pindata'
c      INCLUDE 'comgra'
c      INCLUDE 'colours'
c      include 'printopt'
c      INCLUDE 'slcom'
c      INCLUDE 'slout'

      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost

      INTEGER   nxbin,nybin,nbin,ix,iy,i1,i2,i3,ik,ir,fp,iimg,n1,id,ike
      LOGICAL   outofgrid
      REAL      deltar,deltaz,qmin,qmax,qval,xcen,ycen,p1(4),p2(4),
     .          cxmin,cxmax,cymin,cymax
      REAL*8    maxpixel
      CHARACTER caption*1024,xlabel*36,ylabel*64,file*512,filename*512,
     .          tag_x*2,tag_y*2
      INTEGER, ALLOCATABLE :: nv(:)
      REAL   , ALLOCATABLE :: rv(:,:),zv(:,:),cq(:),image(:,:)


      CALL LINCOL(1)
      CALL FILCOL(1)

      nxbin = opt%nxbin
      nybin = opt%nybin

      ALLOCATE(image(nxbin,nybin))


c...  Copy image data to a properly sized array:
      image(1:nxbin,1:nybin) = SNGL(image_space(1:nxbin,1:nybin))  


c...  Activate MaxEnt + other output:
      IF (opt%nxbin.GT.1.AND.opt%nybin.GT.1) THEN    
        fp = 99
        IF (opt%out_suffix(1:4).NE.'none') THEN
C         Krieger IPP/07 - SUN compiler chokes on format
C         statements w/o explicit length - removed explicit format
C         should be fixed by DIVIMP guru 
c         -not sure what to do about this since specifying the length is a problem... -SL
          WRITE(filename,'(A,3(I0.8,A))')
     .      './images/',
     .      header%shot ,  '_',
     .      header%frame,  '_',
     .      header%channel,'_'//
     .      opt%out_suffix(1:LEN_TRIM(opt%out_suffix))//'.'
        ELSE
          WRITE(filename,'(A,3(I0.8,A))')
     .      './images/',
     .      header%shot ,  '_',
     .      header%frame,  '_',
     .      header%channel,'.'
        ENDIF
c...    Remove blanks:
        n1 = LEN_TRIM(filename)
        i1 = 1
        DO WHILE (i1.LT.n1)
          IF (filename(i1:i1).EQ.' ') THEN
            DO i2 = i1, n1-1
              filename(i2:i2) = filename(i2+1:i2+1)
            ENDDO
            filename(n1:n1) = ' '
            n1 = n1 - 1
          ELSE
            i1 = i1 + 1
          ENDIF
        ENDDO
      ELSE
        fp = 0
      ENDIF
      WRITE(0,*) 'MAXENT:',opt%nxbin,opt%nybin,fp

      IF (nxbin.GE.100) THEN
        nbin = nxbin / 50          ! Should also check opt%nybin
        nbin = 1
      ELSE
        nbin = 1
      ENDIF

      IF (.FALSE.) THEN
c...    Plot image:
c        nbin = 1 !6  ! Binning (NBIN=1 is no binning)
        IF (nbin.GT.1.AND.
     .      (nbin.GE.opt%nxbin.OR.
     .       nbin.GE.opt%nybin))
     .    CALL ER('989','Bin request greater than image resolution',*99)
        nxbin = opt%nxbin / nbin
        nybin = opt%nybin / nbin
        qmin =  0.0
        qmax = -HI
        DO ix = 1, nxbin
          DO iy = 1, nybin
            qval = 0.0
            DO i2 = 0, nbin-1
              DO i3 = 0, nbin-1
                IF (nbin*(ix-1)+1+i2.GT.opt%nxbin) WRITE(0,*) 'xBOUNDS!'
                IF (nbin*(iy-1)+1+i3.GT.opt%nybin) WRITE(0,*) 'yBOUNDS!'
                qval = qval + image(nbin*(ix-1)+1+i2,
     .                              nbin*(iy-1)+1+i3)
              ENDDO
            ENDDO           
            qmax = MAX(qmax,qval)  ! Need to average here, and below?  
          ENDDO
        ENDDO
        WRITE(0,*) 'QMAD :',qmin,qmax
        WRITE(0,*) 'NX,NY:',nxbin,nybin
        map1x = 0.05             ! Add sensitivity to nxbin, nybin...
        map2x = map1x + 0.95  ! 0.45
        map1y = 0.05          ! 0.52
        map2y = map1y + 0.95  ! 0.45
        CALL PSPACE(map1x,map2x,map1y,map2y)
        CALL MAP   (0.0,1.0,1.0,0.0)
        ALLOCATE(nv(1))
        ALLOCATE(rv(4,1))  
        ALLOCATE(zv(4,1))
        ALLOCATE(cq(1))
        deltar = 1.0 / REAL(nxbin)
        deltaz = 1.0 / REAL(nybin)
        DO ix = 1, nxbin
          DO iy = 1, nybin
            i1 = 1
            rv(1,i1) = REAL(ix - 1) * deltar 
            zv(1,i1) = REAL(iy - 1) * deltaz
            rv(2,i1) = REAL(ix - 1) * deltar 
            zv(2,i1) = REAL(iy    ) * deltaz
            rv(3,i1) = REAL(ix    ) * deltar 
            zv(3,i1) = REAL(iy    ) * deltaz
            rv(4,i1) = REAL(ix    ) * deltar 
            zv(4,i1) = REAL(iy - 1) * deltaz
            cq(i1) = 0.0
            DO i2 = 0, nbin-1
              DO i3 = 0, nbin-1
                cq(i1) = cq(i1) + image(nbin*(ix-1)+1+i2,
     .                                  nbin*(iy-1)+1+i3)
              ENDDO
            ENDDO
            CALL SetCol255_04(2,cq(i1),qmin,qmax)
            CALL FILCOL(255)
            CALL LINCOL(255) 
            CALL PTPLOT(rv(1,i1),zv(1,i1),1,4,1)
          ENDDO
        ENDDO
c...    Clear arrays:
        IF (ALLOCATED(nv)) DEALLOCATE(nv)
        IF (ALLOCATED(rv)) DEALLOCATE(rv)
        IF (ALLOCATED(zv)) DEALLOCATE(zv)
        IF (ALLOCATED(cq)) DEALLOCATE(cq)
c...    Annotate:
        CALL LinCol(ncols+1)
        CALL CTRMAG(12)
c        WRITE(caption,'(A,3I5)') 'NX,NY,BIN:',opt%nxbin,opt%nybin,nbin
        CALL PLOTST(0.02,0.03,caption(1:LEN_TRIM(caption)))
c        CALL DrawColourScale(1,2,qmin,qmax,'none')
c        CALL DrawFrame
        CALL Frame
      ENDIF


      IF (.NOT.plot_option) THEN
c        WRITE(0,*) 'INVERTED PROFILE A'
        nxbin = opt%nxbin
        nybin = opt%nybin
        qmin =  0.0
        qmax = -HI
        cxmin =  HI
        cxmax = -HI
        cymin =  HI
        cymax = -HI
        DO i1 = 1, ninv
          qmax = MAX(qmax,qpts(i1,2))
          DO i2 = 1, npts(i1)
            cxmin = MIN(cxmin,xpts(i2,i1))
            cxmax = MAX(cxmax,xpts(i2,i1))
            cymin = MIN(cymin,ypts(i2,i1))
            cymax = MAX(cymax,ypts(i2,i1))
          ENDDO
        ENDDO
        cxmin = cxmin - 0.1 * (cxmax - cxmin)
        cxmax = cxmax + 0.1 * (cxmax - cxmin)
        cymin = cymin - 0.1 * (cymax - cymin)
        cymax = cymax + 0.1 * (cymax - cymin)
        map1x = 0.55             
        map1y = 0.52
        IF (cxmax-cxmin.GT.cymax-cymin) THEN
          map2x = map1x + 0.45
          map2y = map1y + 0.45 * (cymax - cymin) / (cxmax - cxmin)
        ELSE
          map2x = map1x + 0.45 * (cxmax - cxmin) / (cymax - cymin)
          map2y = map1y + 0.45 
        ENDIF
        CALL PSPACE(map1x,map2x,map1y,map2y)
        CALL MAP   (cxmin,cxmax,cymin,cymax)
        DO i1 = 1, ninv
          CALL SetCol255_04(2,qpts(i1,2),qmin,qmax)
          CALL FILCOL(255)
          CALL LINCOL(255) 
          CALL PTPLOT(xpts(1,i1),ypts(1,i1),1,npts(i1),1)
        ENDDO
c...    Annotate:
c        CALL LinCol(ncols+1)
c        CALL CTRMAG(12)
c        WRITE(caption,'(A,3I4)') 'NX,NY,BIN:',opt%nxbin,opt%nybin,nbin
c        CALL PLOTST(0.02,0.02,caption(1:LEN_TRIM(caption)))
c...    Frame:
c        CALL Supimp('PARTIAL')
c        CALL DrawFrame
        CALL Frame
      ENDIF



c...  Plot inverted profile:
      IF (ninv.GT.0.AND.plot_option) THEN
c        WRITE(0,*) 'INVERTED PROFILE B'

        IF (.TRUE..AND.opt%nxbin.GT.1.AND.opt%nybin.GT.1) THEN
c...      Plot camera image:
c          nbin = 1 !6  ! Binning (NBIN=1 is no binning)

          IF (nbin.GT.1.AND.
     .        (nbin.GE.opt%nxbin.OR.
     .         nbin.GE.opt%nybin))
     .      CALL ER('989','Bin request greater than image res.',*99)
          nxbin = opt%nxbin / nbin
          nybin = opt%nybin / nbin


          IF (fp.NE.0) THEN
            IF (nbin.NE.1) CALL ER('Output989','Bin prob',*99)

            file = opt%fmap(1:LEN_TRIM(opt%fmap))//'.ray.b'
            WRITE(0,*) 'DUMP IMAGE:',file(1:LEN_TRIM(file))
            OPEN(fp,FILE=file(1:LEN_TRIM(file)),
     .           FORM='FORMATTED',STATUS='REPLACE',ERR=97)
            WRITE(fp,'(I7,2I12,F12.6,I6)') nxbin*nybin,
     .        header%shot,header%frame,header%time,header%channel

            file = opt%fmap(1:LEN_TRIM(opt%fmap))//'.idl.b'
            write(0,*) file(1:LEN_TRIM(file))
            CALL inOpenInterface(file(1:LEN_TRIM(file)),ITF_WRITE)

            CALL inPutData(header%shot   ,'shot'   ,'none')
            CALL inPutData(header%frame  ,'frame'  ,'none')
            CALL inPutData(header%time   ,'time'   ,'s')
            CALL inPutData(header%channel,'channel','none')

          ENDIF

          qmin =  0.0
          qmax = -HI
          DO ix = 1, nxbin
            DO iy = 1, nybin
              qval = 0.0
              DO i2 = 0, nbin-1
                DO i3 = 0, nbin-1
                  IF (nbin*(ix-1)+1+i2.GT.opt%nxbin) WRITE(0,*) 'xBNDS!'
                  IF (nbin*(iy-1)+1+i3.GT.opt%nybin) WRITE(0,*) 'yBNDS!'
                  qval = qval + image(nbin*(ix-1)+1+i2,
     .                                nbin*(iy-1)+1+i3)
                ENDDO
              ENDDO           
              qmax = MAX(qmax,qval)  ! Need to average here, and below?  
            ENDDO
          ENDDO
          WRITE(0,*) 'QMAD :',qmin,qmax
          WRITE(0,*) 'NX,NY:',nxbin,nybin

          map1x = 0.05            
          map2y = 0.70
          IF (nxbin.GE.nybin) THEN
            map2x = map1x + 0.60
            map1y = map2y - 0.60 * REAL(nybin) / REAL(nxbin)
          ELSE
            map2x = map1x + 0.60 * REAL(nxbin) / REAL(nybin)
            map2y = map1y - 0.60 
          ENDIF 
c          map1x = 0.90 ! 0.05             
c          map1y = 0.52 !
c          map2x = map1x + 0.45
c          map2y = map1y + 0.45
          CALL PSPACE(map1x,map2x,map1y,map2y)
          CALL MAP   (0.0,1.0,1.0,0.0)
          ALLOCATE(nv(1))
          ALLOCATE(rv(4,1))  
          ALLOCATE(zv(4,1))
          ALLOCATE(cq(1))
          deltar = 1.0 / REAL(nxbin)
          deltaz = 1.0 / REAL(nybin)
          DO iy = 1, nybin
            DO ix = 1, nxbin
              i1 = 1
              rv(1,i1) = REAL(ix - 1) * deltar 
              zv(1,i1) = REAL(iy - 1) * deltaz
              rv(2,i1) = REAL(ix - 1) * deltar 
              zv(2,i1) = REAL(iy    ) * deltaz
              rv(3,i1) = REAL(ix    ) * deltar 
              zv(3,i1) = REAL(iy    ) * deltaz
              rv(4,i1) = REAL(ix    ) * deltar 
              zv(4,i1) = REAL(iy - 1) * deltaz
              cq(i1) = 0.0
              DO i2 = 0, nbin-1
                DO i3 = 0, nbin-1
                  cq(i1) = cq(i1) + image(nbin*(ix-1)+1+i2,
     .                                    nbin*(iy-1)+1+i3)
                ENDDO
              ENDDO
              CALL SetCol255_04(2,cq(i1),qmin,qmax)

              IF (fp.GT.0) THEN
                WRITE(fp,'(I7,1P,E22.12)') ix+(iy-1)*nxbin,cq(i1)
                CALL inPutData(ix+(iy-1)*nxbin,'i','none')
                CALL inPutData(cq(i1),'data','ph m-2 s-1')                            
              ENDIF

              CALL FILCOL(255)
              CALL LINCOL(255) 
              CALL PTPLOT(rv(1,i1),zv(1,i1),1,4,1)
            ENDDO
          ENDDO
c...      Clear arrays:
          IF (ALLOCATED(nv)) DEALLOCATE(nv)
          IF (ALLOCATED(rv)) DEALLOCATE(rv)
          IF (ALLOCATED(zv)) DEALLOCATE(zv)
          IF (ALLOCATED(cq)) DEALLOCATE(cq)
c...      Annotate:
          CALL LinCol(ncols+1)
          CALL CTRMAG(12)
c          WRITE(caption,'(A,3I5)') 'NX,NY,BIN:',opt%nxbin,opt%nybin,nbin
          CALL PLOTST(0.02,0.03,caption(1:LEN_TRIM(caption)))
c          CALL DrawColourScale(2,2,qmin,qmax,'none')
c          CALL DrawFrame

          IF (fp.GT.0) THEN
            CLOSE(fp)
            CALL inCloseInterface
          ENDIF
        ENDIF

c        CALL LoadTriangles

        nxbin = opt%nxbin
        nybin = opt%nybin
        qmin =  0.0
        qmax = -HI
        cxmin =  HI
        cxmax = -HI
        cymin =  HI
        cymax = -HI
        DO i1 = 1, ninv
          qmax = MAX(qmax,qpts(i1,2))
          DO i2 = 1, npts(i1)
            cxmin = MIN(cxmin,xpts(i2,i1))
            cxmax = MAX(cxmax,xpts(i2,i1))
            cymin = MIN(cymin,ypts(i2,i1))
            cymax = MAX(cymax,ypts(i2,i1))
          ENDDO
        ENDDO
        WRITE(0,*) 'PROFILE QMIN,QMAX :',qmin,qmax

c...    HACK!
c        IF (cymin.LT.0.0.AND.cymax.GT.0.0) THEN 
c          WRITE(0,*) 'WARNING! SUPER-HACK ON PLOT RANGE'
c          cymax = cymin + 0.5 * (cymax - cymin)
c        ENDIF

        cxmin = cxmin - 0.05 * (cxmax - cxmin)
        cxmax = cxmax + 0.05 * (cxmax - cxmin)
        cymin = cymin - 0.05 * (cymax - cymin)
        cymax = cymax + 0.05 * (cymax - cymin)


        IF (.FALSE.) THEN
c...      Reference profile:
          map1x = 0.90 !05             
          map1y = 0.52 !02
          IF (cxmax-cxmin.GT.cymax-cymin) THEN
            map2x = map1x + 0.45
            map2y = map1y + 0.45 * (cymax - cymin) / (cxmax - cxmin)
          ELSE
            map2x = map1x + 0.45 * (cxmax - cxmin) / (cymax - cymin)
            map2y = map1y + 0.45 
          ENDIF
          CALL PSPACE(map1x,map2x,map1y,map2y)
          CALL MAP   (cxmin,cxmax,cymin,cymax)
          DO i1 = 1, ninv
            CALL SetCol255_04(2,qpts(i1,1),qmin,qmax)
            CALL FILCOL(255)
            CALL LINCOL(255) 
            CALL PTPLOT(xpts(1,i1),ypts(1,i1),1,npts(i1),1)
          ENDDO
c...      Annotate:
c          CALL LinCol(ncols+1)
c          CALL CTRMAG(12)
c          WRITE(caption,'(A,3I4)') 'NX,NY,BIN:',opt%nxbin,opt%nybin,nbin
c          CALL PLOTST(0.02,0.02,caption(1:LEN_TRIM(caption)))
c          CALL DrawGrid(95)
c          CALL DrawColourScale(2,2,qmin,qmax,'none')
c          CALL DrawFrame

          IF (fp.GT.0) THEN
            file = opt%fmap(1:LEN_TRIM(opt%fmap))//'.ray.x-ref'
            WRITE(0,*) 'DUMP IMAGE:',file(1:LEN_TRIM(file))
            OPEN(fp,FILE=file(1:LEN_TRIM(file)),
     .           FORM='FORMATTED',STATUS='REPLACE',ERR=97)
            WRITE(fp,'(I7)') ninv
            DO i1 = 1, ninv
              WRITE(fp,'(I7,1P,E22.12)') i1,qpts(i1,1)
            ENDDO     
            CLOSE(fp)
          ENDIF

        ENDIF

c...    Inverted profile:
        IF (fp.GT.0) THEN

          file = filename(1:LEN_TRIM(filename))//'cgm'
          WRITE(0,*) 'DUMP INVERSION:',file(1:LEN_TRIM(file))
          OPEN(fp,FILE=file(1:LEN_TRIM(file)),
     .         FORM='FORMATTED',STATUS='REPLACE',ERR=97)
          WRITE(fp,'(I7,2I12,F12.6,I6)') ninv,
     .          header%shot,header%frame,header%time,header%channel
          DO i1 = 1, ninv
            WRITE(fp,'(I7,1P,E22.12,0P,I5,10(2F11.6))')  
     .        i1,qpts(i1,2),
     .        npts(i1),(xpts(i2,i1),ypts(i2,i1),i2=1,npts(i1))
          ENDDO     
          CLOSE(fp)

          file = opt%fmap(1:LEN_TRIM(opt%fmap))//'.ray.x'
          OPEN(fp,FILE=file(1:LEN_TRIM(file)),
     .         FORM='FORMATTED',STATUS='REPLACE',ERR=97)
          WRITE(fp,'(I7,2I12,F12.6,I6)') ninv,
     .          header%shot,header%frame,header%time,header%channel
          WRITE(file,'(514X)')          
          file = opt%fmap(1:LEN_TRIM(opt%fmap))//'.idl.x'
          CALL inOpenInterface(file(1:LEN_TRIM(file)),ITF_WRITE)
          
          DO i1 = 1, ninv
            WRITE(fp,'(I7,1P,E22.12,0P,I5,10(2F11.6))')  
     .        i1,qpts(i1,2),
     .        npts(i1),(xpts(i2,i1),ypts(i2,i1),i2=1,npts(i1))
           CALL inPutData(i1,'i','none')
           CALL inPutData(npts(i1),'npts','none')            
           DO i2 = 1, npts(i1)
             WRITE(tag_x,'(A,I1)') 'x',i2
             WRITE(tag_y,'(A,I1)') 'y',i2
             CALL inPutData(xpts(i2,i1),tag_x,'m')            
             CALL inPutData(ypts(i2,i1),tag_y,'m')            
           ENDDO
           CALL inPutData(qpts(i1,2),'data','ph m-3 s-1')            
          ENDDO     
          CLOSE(fp)
          CALL inCloseInterface
  
        ENDIF

        map1x = 0.70 ! 0.90 !55
        map1y = 0.10 ! 0.02
        IF (cxmax-cxmin.GT.cymax-cymin) THEN
          map2x = map1x + 0.60
          map2y = map1y + 0.60 * (cymax - cymin) / (cxmax - cxmin)
        ELSE
          map2x = map1x + 0.60 * (cxmax - cxmin) / (cymax - cymin)
          map2y = map1y + 0.60 
        ENDIF

        CALL PSPACE(map1x,map2x,map1y,map2y)
        CALL MAP   (cxmin,cxmax,cymin,cymax)
        DO i1 = 1, ninv
          IF (qpts(i1,2).GT.0.0) THEN
            CALL SetCol255_04(2,qpts(i1,2),qmin,qmax)
          ELSE
            CALL SetCol255_04(2,0.0       ,qmin,qmax)
          ENDIF
          CALL FILCOL(255)
          CALL LINCOL(255) 
          CALL PTPLOT(xpts(1,i1),ypts(1,i1),1,npts(i1),1)
          IF (.FALSE.) THEN
            IF (tpts(i1).EQ.0.0) THEN
c...          Inversion mesh cell not sampled by LOS integration:
              IF (qpts(i1,2).EQ.0.0) THEN
                CALL LINCOL(ncols+4) 
              ELSE
                CALL LINCOL(ncols+1) 
              ENDIF
              DO i2 = 1, npts(i1)/2
                i3 = npts(i1)/2 + i2
c              IF (i3.EQ.npts(i1)+1) i3 = 1
                CALL POSITN(xpts(i2,i1),ypts(i2,i1))
                CALL JOIN  (xpts(i3,i1),ypts(i3,i1))
              ENDDO
            ELSEIF (qpts(i1,2).LT.0.0) THEN
c...          Negative emissivity:
              CALL LINCOL(ncols+5) 
              DO i2 = 1, npts(i1)/2
                i3 = npts(i1)/2 + i2
                CALL POSITN(xpts(i2,i1),ypts(i2,i1))
                CALL JOIN  (xpts(i3,i1),ypts(i3,i1))
              ENDDO
            ENDIF
          ENDIF
        ENDDO
c...    Annotate:
c        CALL LinCol(ncols+1)
c        CALL CTRMAG(12)
c        WRITE(caption,'(A,3I4)') 'NX,NY,BIN:',opt%nxbin,opt%nybin,nbin
c        CALL PLOTST(0.02,0.02,caption(1:LEN_TRIM(caption)))
c...    Frame:
c        CALL DrawGrid(-95)
c        CALL Supimp('PARTIAL')
c        Call DrawFrame

c        CALL Frame


c        ninv = ninv / 4


        CALL Frame
      ENDIF



c...  Presentation plots:
      IF (.FALSE..AND.ninv.GT.0) THEN
c        WRITE(0,*) 'INVERTED PROFILE B'


        IF (.TRUE..AND.opt%nxbin.GT.1.AND.opt%nybin.GT.1) THEN
c...      Camera image:

          IF (nbin.GT.1.AND.
     .        (nbin.GE.opt%nxbin.OR.
     .         nbin.GE.opt%nybin))
     .      CALL ER('989','Bin request greater than image res.',*99)
          nxbin = opt%nxbin / nbin
          nybin = opt%nybin / nbin

          qmin =  0.0
          qmax = -HI
          DO ix = 1, nxbin
            DO iy = 1, nybin
              qval = 0.0
              DO i2 = 0, nbin-1
                DO i3 = 0, nbin-1
                  IF (nbin*(ix-1)+1+i2.GT.opt%nxbin) WRITE(0,*) 'xBNDS!'
                  IF (nbin*(iy-1)+1+i3.GT.opt%nybin) WRITE(0,*) 'yBNDS!'
                  qval = qval + image(nbin*(ix-1)+1+i2,
     .                                nbin*(iy-1)+1+i3)
                ENDDO
              ENDDO           
              qmax = MAX(qmax,qval)  ! Need to average here, and below?  
            ENDDO
          ENDDO

c          qmax = qmax / 5.0
            
          WRITE(0,*) 'QMAD :',qmin,qmax
          WRITE(0,*) 'NX,NY:',nxbin,nybin
          map1x = 0.05             
          map1y = 0.05 
          map2x = map1x + 0.90
          map2y = map1y + 0.90
          CALL PSPACE(map1x,map2x,map1y,map2y)
          CALL MAP   (0.0,1.0,1.0,0.0)
          ALLOCATE(nv(1))
          ALLOCATE(rv(4,1))  
          ALLOCATE(zv(4,1))
          ALLOCATE(cq(1))
          deltar = 1.0 / REAL(nxbin)
          deltaz = 1.0 / REAL(nybin)
          DO iy = 1, nybin
            DO ix = 1, nxbin
              i1 = 1
              rv(1,i1) = REAL(ix - 1) * deltar 
              zv(1,i1) = REAL(iy - 1) * deltaz
              rv(2,i1) = REAL(ix - 1) * deltar 
              zv(2,i1) = REAL(iy    ) * deltaz
              rv(3,i1) = REAL(ix    ) * deltar 
              zv(3,i1) = REAL(iy    ) * deltaz
              rv(4,i1) = REAL(ix    ) * deltar 
              zv(4,i1) = REAL(iy - 1) * deltaz
              cq(i1) = 0.0
              DO i2 = 0, nbin-1
                DO i3 = 0, nbin-1
                  cq(i1) = cq(i1) + image(nbin*(ix-1)+1+i2,
     .                                    nbin*(iy-1)+1+i3)
                ENDDO
              ENDDO
              CALL SetCol255_04(2,MIN(qmax,cq(i1)),qmin,qmax)

              CALL FILCOL(255)
              CALL LINCOL(255) 
              CALL PTPLOT(rv(1,i1),zv(1,i1),1,4,1)
            ENDDO
          ENDDO
c...      Clear arrays:
          IF (ALLOCATED(nv)) DEALLOCATE(nv)
          IF (ALLOCATED(rv)) DEALLOCATE(rv)
          IF (ALLOCATED(zv)) DEALLOCATE(zv)
          IF (ALLOCATED(cq)) DEALLOCATE(cq)
c...      Annotate:
          CALL LinCol(ncols+1)
          CALL CTRMAG(12)
c          WRITE(caption,'(A,3I5)') 'NX,NY,BIN:',opt%nxbin,opt%nybin,nbin
          CALL PLOTST(0.02,0.03,caption(1:LEN_TRIM(caption)))
c          CALL DrawColourScale(2,2,qmin,qmax,'none')
c          CALL DrawFrame
        ENDIF


c        CALL LoadTriangles

        nxbin = opt%nxbin
        nybin = opt%nybin

        WRITE(0,*) 'C VALUES:',cxmin,cxmax,cymin,cymax

        IF (.TRUE.) THEN
c...      Reference profile:

          CALL Frame

          qmin =  0.0
          qmax = -HI
          cxmin =  HI
          cxmax = -HI
          cymin =  HI
          cymax = -HI
          DO i1 = 1, ninv
            qmax = MAX(qmax,qpts(i1,1))
            DO i2 = 1, npts(i1)
              cxmin = MIN(cxmin,xpts(i2,i1))
              cxmax = MAX(cxmax,xpts(i2,i1))
              cymin = MIN(cymin,ypts(i2,i1))
              cymax = MAX(cymax,ypts(i2,i1))
            ENDDO
          ENDDO
          WRITE(0,*) 'PROFILE QMIN,QMAX 1 :',qmin,qmax

          cxmin = cxmin - 0.05 * (cxmax - cxmin)
          cxmax = cxmax + 0.05 * (cxmax - cxmin)
          cymin = cymin - 0.05 * (cymax - cymin)
          cymax = cymax + 0.05 * (cymax - cymin)

          WRITE(0,*) 'PROFILE CDATA 1 :',cxmin,cxmax,cymin,cymax

          map1x = 0.05
          map1y = 0.05
          IF (cxmax-cxmin.GT.cymax-cymin) THEN
            map2x = map1x + 0.90
            map2y = map1y + 0.90 * (cymax - cymin) / (cxmax - cxmin)
          ELSE
            map2x = map1x + 0.90 * (cxmax - cxmin) / (cymax - cymin)
            map2y = map1y + 0.90 
          ENDIF
          CALL PSPACE(map1x,map2x,map1y,map2y)
          CALL MAP   (cxmin,cxmax,cymin,cymax)

c...      Black background for plot:
          ALLOCATE(rv(4,1))  
          ALLOCATE(zv(4,1))
          rv(1,1) = cxmax
          zv(1,1) = cymin
          rv(2,1) = cxmin
          zv(2,1) = cymin
          rv(3,1) = cxmin
          zv(3,1) = cymax
          rv(4,1) = cxmax
          zv(4,1) = cymax
          CALL SetCol255_04(2,0.0,qmin,qmax)       ! Black
          CALL PTPLOT(rv(1,1),zv(1,1),1,4,1)
          IF (ALLOCATED(rv)) DEALLOCATE(rv)
          IF (ALLOCATED(zv)) DEALLOCATE(zv)

          DO i1 = 1, ninv
            CALL SetCol255_04(2,qpts(i1,1),qmin,qmax)
            CALL FILCOL(255)
            CALL LINCOL(255) 
            CALL PTPLOT(xpts(1,i1),ypts(1,i1),1,npts(i1),1)
          ENDDO
c...      Annotate:
c          CALL LinCol(ncols+1)
c          CALL CTRMAG(12)
c          WRITE(caption,'(A,3I4)') 'NX,NY,BIN:',opt%nxbin,opt%nybin,nbin
c          CALL PLOTST(0.02,0.02,caption(1:LEN_TRIM(caption)))
c          CALL DrawGrid(95)
c          CALL DrawColourScale(2,2,qmin,qmax,'none')
c          CALL DrawFrame


        ENDIF


c...    Inverted profile:
        CALL Frame

        qmin =  0.0
        qmax = -HI
        cxmin =  HI
        cxmax = -HI
        cymin =  HI
        cymax = -HI
        DO i1 = 1, ninv
          qmax = MAX(qmax,qpts(i1,2))
          DO i2 = 1, npts(i1)
            cxmin = MIN(cxmin,xpts(i2,i1))
            cxmax = MAX(cxmax,xpts(i2,i1))
            cymin = MIN(cymin,ypts(i2,i1))
            cymax = MAX(cymax,ypts(i2,i1))
          ENDDO
        ENDDO
        WRITE(0,*) 'PROFILE QMIN,QMAX 2 :',qmin,qmax

        cxmin = cxmin - 0.05 * (cxmax - cxmin)
        cxmax = cxmax + 0.05 * (cxmax - cxmin)
        cymin = cymin - 0.05 * (cymax - cymin)
        cymax = cymax + 0.05 * (cymax - cymin)

        WRITE(0,*) 'PROFILE CDATA 2 :',cxmin,cxmax,cymin,cymax

        map1x = 0.05
        map1y = 0.05
        IF (cxmax-cxmin.GT.cymax-cymin) THEN
          map2x = map1x + 0.90
          map2y = map1y + 0.90 * (cymax - cymin) / (cxmax - cxmin)
        ELSE
          map2x = map1x + 0.90 * (cxmax - cxmin) / (cymax - cymin)
          map2y = map1y + 0.90 
        ENDIF
        CALL PSPACE(map1x,map2x,map1y,map2y)
        CALL MAP   (cxmin,cxmax,cymin,cymax)

c...    Black background for plot:
        ALLOCATE(rv(4,1))  
        ALLOCATE(zv(4,1))
        rv(1,1) = cxmax
        zv(1,1) = cymin
        rv(2,1) = cxmin
        zv(2,1) = cymin
        rv(3,1) = cxmin
        zv(3,1) = cymax
        rv(4,1) = cxmax
        zv(4,1) = cymax
        CALL SetCol255_04(2,0.0,qmin,qmax)       ! Black
        CALL PTPLOT(rv(1,1),zv(1,1),1,4,1)
        IF (ALLOCATED(rv)) DEALLOCATE(rv)
        IF (ALLOCATED(zv)) DEALLOCATE(zv)

        DO i1 = 1, ninv
          IF (qpts(i1,2).GT.0.0) THEN
            CALL SetCol255_04(2,qpts(i1,2),qmin,qmax)
          ELSE
            CALL SetCol255_04(2,0.0       ,qmin,qmax)
          ENDIF
          CALL FILCOL(255)
          CALL LINCOL(255) 
          CALL PTPLOT(xpts(1,i1),ypts(1,i1),1,npts(i1),1)
        ENDDO
c...    Annotate:
c        CALL LinCol(ncols+1)
c        CALL CTRMAG(12)
c        WRITE(caption,'(A,3I4)') 'NX,NY,BIN:',opt%nxbin,opt%nybin,nbin
c        CALL PLOTST(0.02,0.02,caption(1:LEN_TRIM(caption)))
c...    Frame:
c        CALL DrawGrid(95)
c        CALL DrawColourScale(2,2,qmin,qmax,'none')
c        Call DrawFrame

        IF (.FALSE.) THEN
          CALL Frame
          CALL PSPACE(map1x,map2x,map1y,map2y)
          CALL MAP   (cxmin,cxmax,cymin,cymax)
c...      Background:
          ALLOCATE(rv(4,1))  
          ALLOCATE(zv(4,1))
          rv(1,1) = cxmax
          zv(1,1) = cymin
          rv(2,1) = cxmin
          zv(2,1) = cymin
          rv(3,1) = cxmin
          zv(3,1) = cymax
          rv(4,1) = cxmax
          zv(4,1) = cymax
          CALL SetCol255_04(2,0.0,qmin,qmax)       ! Black
          CALL PTPLOT(rv(1,1),zv(1,1),1,4,1)
          IF (ALLOCATED(rv)) DEALLOCATE(rv)
          IF (ALLOCATED(zv)) DEALLOCATE(zv)
          DO i1 = 1, ninv
            IF (qpts(i1,2).GT.0.0) THEN
              CALL SetCol255_04(2,qpts(i1,2),qmin,qmax)
            ELSE
              CALL SetCol255_04(2,0.0       ,qmin,qmax)
            ENDIF
            CALL FILCOL(255)
            CALL LINCOL(255) 
            CALL PTPLOT(xpts(1,i1),ypts(1,i1),1,npts(i1),1)
          ENDDO
c...      Annotate:
c...      Frame:
c          CALL DrawGrid(-95)
c          CALL DrawColourScale(2,2,qmin,qmax,'none')
c          Call DrawFrame
        ENDIF

        CALL Frame

      ENDIF




c      IF (nbin.GT.1) THEN
c        nxbin = opt%nxbin
c        nybin = opt%nybin
c        maxpixel = DBLE(MAXVAL(image))
c        WRITE(0,*) 'SAVING JPEG'
c        CALL Save_Image(image,nxbin,nybin,0.0,maxpixel,0.0D0)  ! Generate a jpeg file 
c        WRITE(0,*) 'DONE'
c      ENDIF

      IF (ALLOCATED(image)) DEALLOCATE(image)

c...  Clear triangle arrays:
c      CALL DEALLOC_VERTEX  
c      CALL DEALLOC_SURFACE   
c      CALL DEALLOC_TRIANGLE

      WRITE(0,*) '  PLOTS DONE'



      RETURN
 97   CALL ER('Output989','File access problem',*99)
 98   CALL ER('Output989','Trouble opening .raw.idl',*99)
 99   WRITE(0,*) '  FILE=',file(1:LEN_TRIM(file))
      STOP
      END

c
c ======================================================================
c
c subroutine: SetCol255
c
c
      SUBROUTINE SetCol255_04(mode,qval,qmin,qmax)
      IMPLICIT none

      INTEGER mode 
      REAL    qval,qmin,qmax
    

c * OLD *
      REAL frac1,hue,pastel,midfrac,scale
      COMMON /PLOT982COM/ hardscale,colscale
      LOGICAL             hardscale,colscale
      LOGICAL grayscale

      INTEGER last_mode
      REAL bright
      REAL frac,frac5,fmod5,last_frac,last_bright   ! mode.eq.1

      DATA last_mode,last_frac,last_bright /-1, -1.0, -1.0/    

      SAVE


      IF     (mode.EQ.1) THEN

        frac = (qval - qmin) / (qmax - qmin + 1.0E-10)
        frac5 = 100.0*frac
        fmod5 = AMOD(frac5,2.0)
        frac = MIN(0.98,(frac5-fmod5)/100.0)

        bright = 1.0-(0.98-frac)**20

        IF (mode.NE.last_mode.OR.frac.NE.last_frac.OR.
     .      bright.NE.last_bright) 
     .    CALL ColSet(1.0-0.75*frac,1.0,bright,255)
c     .    CALL ColSet(0.75*frac+0.25,1.0,bright,255)

c      ELSEIF (mode.EQ.2) THEN

      ELSEIF (mode.EQ.2) THEN
c...    Grayscale:

        frac = (qval - qmin) / (qmax - qmin)

c        frac = frac * 0.9 + 0.1

        IF (mode.NE.last_mode.OR.frac.NE.last_frac) 
c     .    CALL ColSet(0.0,0.0,frac,255)  ! Default
c     .    CALL ColSet(0.50,1.0-frac,frac, ! Original DIVCAM
c     .                255)
     .    CALL ColSet(0.30+0.230*frac,
     .                MIN(1.0,4.0*(1.0-frac)    ),
     .                MIN(1.0,4.0*     frac     ),
     .                255)
c     .    CALL ColSet(0.25+1.000*frac,
c     .                MIN(1.0,4.0*(1.0-frac)    ),
c     .                MIN(1.0,4.0*     frac     ),
c     .                255)
c     .    CALL ColSet(0.50,1.0-2.0*(MAX(0.0,frac**0.5-0.5)),frac**0.5,
c     .                255)

c        CALL ColSet(0.0,0.0,frac,255)
c        CALL ColSet(0.0,0.0,1.0-frac,255)


      ELSE
        CALL ER('SetCol255_04','Invalid mode',*99)
      ENDIF


      last_mode = mode
      last_frac = frac
      last_bright = bright

      RETURN

99    STOP
      END


