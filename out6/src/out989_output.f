c
c
c
      SUBROUTINE DrawFrame
      IMPLICIT none
      INCLUDE 'params'
      INCLUDE 'slout'
c...  Finish plot:
      CALL LINCOL(1)
      CALL PSPACE(MAP1X,MAP2X,MAP1Y,MAP2Y)
      CALL MAP   (0.0,1.0,0.0,1.0)
      CALL POSITN(0.0,0.0)
      CALL JOIN  (0.0,1.0)              
      CALL POSITN(0.0,1.0)
      CALL JOIN  (1.0,1.0)              
      CALL POSITN(1.0,1.0)
      CALL JOIN  (1.0,0.0)              
      CALL POSITN(1.0,0.0)
      CALL JOIN  (0.0,0.0)                      
      RETURN
      STOP
      END
c
c
c
      SUBROUTINE Output989(opt,nxbin_space,nybin_space,image_space,
     .                     image2,header,header2,
     .                     ninv,npts,xpts,ypts,tpts,qpts)
      USE mod_out989
      USE mod_eirene04
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

      INTEGER   nxbin,nybin,nbin,ix,iy,i1,i2,i3,ik,ir,fp,iimg,n1,id,
     .          ike,n
      LOGICAL   outofgrid
      REAL      deltar,deltaz,qmin,qmax,qval,xcen,ycen,p1(4),p2(4),t,
     .          x1,x2,y1,y2
      REAL*8    maxpixel
      CHARACTER caption*1024,xlabel*36,ylabel*64,file*512,filename*512,
     .          tag_x*11,tag_y*11
      INTEGER, ALLOCATABLE :: nv(:)
      REAL   , ALLOCATABLE :: rv(:,:),zv(:,:),cq(:),image(:,:),
     .                        x(:),y(:),v(:),s(:)


      LOGICAL, PARAMETER :: plot_option = .TRUE.

      CALL LINCOL(1)
      CALL FILCOL(1)


      nxbin = opt%nxbin
      nybin = opt%nybin

      ALLOCATE(image(nxbin,nybin))


c...  Copy image data to a properly sized array:
      image(1:nxbin,1:nybin) = SNGL(image_space(1:nxbin,1:nybin))  


c...  Activate MaxEnt + other output:
      IF (opt%nxbin.GT.1.OR.opt%nybin.GT.1) THEN    
        fp = 99
        IF (opt%out_suffix(1:4).NE.'none') THEN
C         Krieger IPP/07 - SUN compiler chokes on format
C         statements w/o explicit length - removed explicit format
C         should be fixed by DIVIMP guru 
c         -not sure what to do about this since specifying the length is a problem... -SL
          WRITE(filename,'(A,3(I,A))')
     .      './images/',
     .      header%shot ,  '_',
     .      header%frame,  '_',
     .      header%channel,'_'//
     .      opt%out_suffix(1:LEN_TRIM(opt%out_suffix))//'.'
        ELSE
          WRITE(filename,'(A,3(I,A))')
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
        CALL DrawColourScale(1,2,qmin,qmax,'none')
        CALL DrawFrame
        CALL Frame
      ENDIF


      IF (qpts(1,1).EQ.-998.0) THEN
c...    Assign, in some way, the integration quantity in the poloidal plane
c       to the inversion mesh that has just been loaded:
        WRITE(0,*) '  PULLING DALPHA FROM EIRENE, CRUDE INTERPOLATION'
        ik = 1
        ir = irsep
        DO i1 = 1, ninv
          IF (.FALSE.) THEN
c...        
            qpts(i1,1) = 0.0
            xcen = 0.0
            ycen = 0.0
            ncel = 0
            DO i2 = 1, npts(i1)
              xcen = xcen + xpts(i2,i1)
              ycen = ycen + ypts(i2,i1)
              ik = 1
              ir = irsep-1
              CALL GridPos(ik,ir,xpts(i2,i1),ypts(i2,i1),
     .                     .TRUE.,outofgrid)
              IF (.NOT.outofgrid) THEN
                qpts(i1,1)=qpts(i1,1)+pinalpha(ik,ir)
                ncel = ncel + 1
              ENDIF
            ENDDO
            xcen = xcen / REAL(npts(i1))
            ycen = ycen / REAL(npts(i1))
            ik = 1
            ir = irsep-1
            CALL GridPos(ik,ir,xcen,ycen,.TRUE.,outofgrid)
            IF (.NOT.outofgrid) THEN
              qpts(i1,1) = qpts(i1,1)+pinalpha(ik,ir)
              ncel = ncel + 1
            ENDIF
            IF (ncel.GT.0) THEN
              qpts(i1,1) = qpts(i1,1) / REAL(ncel)
            ELSE
              qpts(i1,1) = 0.0
            ENDIF
          ELSEIF (.TRUE.) THEN
c...        Find cell geometry center (some assumptions):
            xcen = 0.0
            ycen = 0.0
            DO i2 = 1, npts(i1)
              xcen = xcen + xpts(i2,i1)
              ycen = ycen + ypts(i2,i1)
            ENDDO
            xcen = xcen / REAL(npts(i1))
            ycen = ycen / REAL(npts(i1))
c...        Find which magentic grid cell the inversion cell is in:
            CALL GridPos(ik,ir,xcen,ycen,.TRUE.,outofgrid)
            IF (.NOT.outofgrid) THEN
              qpts(i1,1) = pinalpha(ik,ir)
            ELSE
              qpts(i1,1) = 0.0
            ENDIF
          ENDIF

        ENDDO

        file = filename(1:LEN_TRIM(filename))//'osm_interpol'
        WRITE(0,*) 'DUMP INTERPOLATION:',file(1:LEN_TRIM(file))
        OPEN(fp,FILE=file(1:LEN_TRIM(file)),
     .       FORM='FORMATTED',STATUS='REPLACE',ERR=97)
        WRITE(fp,'(I7,2I12,F12.6,I6)') ninv,
     .        header%shot,header%frame,header%time,header%channel
        DO i1 = 1, ninv
          WRITE(fp,'(I7,1P,E22.12,0P,I5,10(2F11.6))')  
     .      i1,qpts(i1,1),
     .      npts(i1),(xpts(i2,i1),ypts(i2,i1),i2=1,npts(i1))
        ENDDO     
        CLOSE(fp)

        WRITE(0,*) '  DONE'
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
        CALL Supimp('PARTIAL')
        CALL DrawFrame
        CALL Frame
      ENDIF



c...  Plot inverted profile:
      IF (ninv.GT.0.AND.plot_option) THEN
c        WRITE(0,*) 'INVERTED PROFILE B'

        IF (.TRUE..AND.(opt%nxbin.GT.1.OR.opt%nybin.GT.1)) THEN
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
            CALL inOpenInterface(file,ITF_WRITE)
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
          map1x = 0.90 ! 0.05             
          map1y = 0.52 !
          map2x = map1x + 0.45
          map2y = map1y + 0.45
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
          CALL DrawColourScale(2,2,qmin,qmax,'none')
          CALL DrawFrame

          IF (fp.GT.0) THEN
            CLOSE(fp)
            CALL inCloseInterface
          ENDIF
        ENDIF


        CALL LoadTriangles

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
          CALL DrawGrid(95)
          CALL DrawColourScale(2,2,qmin,qmax,'none')
          CALL DrawFrame

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
          WRITE(fp,'(I7,2I12,F12.6,I6,5X,A)') ninv,
     .          header%shot,header%frame,header%time,header%channel,
     .          opt%out_suffix(1:LEN_TRIM(opt%out_suffix))
          WRITE(file,'(514X)')          
          file = opt%fmap(1:LEN_TRIM(opt%fmap))//'.idl.x'
          CALL inOpenInterface(file,ITF_WRITE)   ! TRIM(file) would not work, compiler bug...
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

        map1x = 0.90 !55
        map1y = 0.02
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
        CALL DrawGrid(-95)
c        CALL Supimp('PARTIAL')
        Call DrawFrame

c        CALL Frame


c        ninv = ninv / 4


c...    Overlay known emission profile and inversion:
        IF (.FALSE.) THEN

          slopt2 = 1
          iopt_ghost = 1 ! 2
          plottype(1) = 2
          plottype(2) = 3
          map1x = 0.05             
          map2x = map1x + 0.80 !95
          map1y = 0.07 !55
          map2y = map1y + 0.15
          cxmin = 1.0
          cxmax = REAL(ninv)
          cymin =  HI
          cymax = -HI
          DO i1 = 1, ninv
            cymin = MIN(cymin,qpts(i1,1))
            cymax = MAX(cymax,qpts(i1,1))
          ENDDO
          DO i1 = 1, ninv
            xdat(i1) = REAL(i1)
          ENDDO
          xlabel = 'pixel    '
          ylabel = 'signal   '
          CALL GRTSET_TRIM(' ',' ',' ',' ',' ',      ! TITLE,REF,nVIEW,PLANE,glabel,
     .                     cxmin,cxmax,cymin,cymax,
     .                     ' ',xlabel,ylabel,        ! TABLE,XLAB,YLAB,
     .                     0,' ',0,' ',1)            ! 0,smooth,0,ANLY,1)
          CALL GRTRAC(xdat,qpts(1,1),ninv,'ref ','LINE',1)
c          DO i1 = 1, ninv
c            WRITE(6,*) 'WHAT THE?',i1,qpts(i1,1)
c          ENDDO
          ydat(1:ninv) = qpts(1:ninv,2)
          DO i1 = 1, ninv
            IF (ydat(i1).EQ.0.0) ydat(i1) = LO
          ENDDO
          CALL GRTRAC(xdat,ydat,ninv,'cal ','LINE',1)
          CALL DrawFrame
        
c...      %diff:
          DO i1 = 1, ninv
            xdat(i1) = REAL(i1)
            IF (qpts(i1,1).NE.0.0.AND.tpts(i1).NE.0.0) THEN
              ydat(i1) = ABS(qpts(i1,2)-qpts(i1,1))/qpts(i1,1)*100.0
            ELSE
              ydat(i1) = 0.0
            ENDIF
          ENDDO
          xlabel = 'none     '
          ylabel = '%        '
          cxmin = 1.0
          cxmax = REAL(ninv)
          cymin = 0.0
          cymax = 100.0
          map1x = 0.05             ! Add sensitivity to nxbin, nybin...
          map2x = map1x + 0.80 !95
          map1y = map2y
          map2y = map1y + 0.15
          CALL GRTSET_TRIM(' ',' ',' ',' ',' ',      ! TITLE,REF,nVIEW,PLANE,glabel,
     .                     cxmin,cxmax,cymin,cymax,
     .                     ' ',xlabel,ylabel,
     .                     0,' ',0,' ',1)            ! 0,smooth,0,ANLY,1)
          iopt_ghost = 1 ! 2
          CALL GRTRAC(xdat,ydat,ninv,'%df ','LINE',1)
          CALL DrawFrame
c...      Just inversion:
          map1x = 0.05             ! Add sensitivity to nxbin, nybin...
          map2x = map1x + 0.80 !95
          map1y = map2y
          map2y = map1y + 0.15
          cxmin = 1.0
          cxmax = REAL(ninv)
          cymin =  HI
          cymax = -HI
          DO i1 = 1, ninv
            cymin = MIN(cymin,qpts(i1,2))
            cymax = MAX(cymax,qpts(i1,2))
          ENDDO
          xlabel = 'none     '
          ylabel = '???      '
          CALL GRTSET_TRIM(' ',' ',' ',' ',' ',      ! TITLE,REF,nVIEW,PLANE,glabel,
     .                     cxmin,cxmax,cymin,cymax,
     .                     ' ',xlabel,ylabel,        ! TABLE,XLAB,YLAB,
     .                     0,' ',0,' ',1)            ! 0,smooth,0,ANLY,1)
          CALL GRTRAC(xdat,qpts(1,2),ninv,'inv ','LINE',1)
          CALL DrawFrame
c...      Sample weight:
          map1x = 0.05             ! Add sensitivity to nxbin, nybin...
          map2x = map1x + 0.80 
          map1y = map2y
          map2y = map1y + 0.15
          cxmin = 1.0
          cxmax = REAL(ninv)
          cymin =  HI
          cymax = -HI
          DO i1 = 1, ninv
            cymin = MIN(cymin,tpts(i1))
            cymax = MAX(cymax,tpts(i1))
          ENDDO
          xlabel = 'none     '
          ylabel = '???      '
          CALL GRTSET_TRIM(' ',' ',' ',' ',' ',      ! TITLE,REF,nVIEW,PLANE,glabel,
     .                     cxmin,cxmax,cymin,cymax,
     .                     ' ',xlabel,ylabel,        ! TABLE,XLAB,YLAB,
     .                     0,' ',0,' ',1)            ! 0,smooth,0,ANLY,1)
          ydat(1:ninv) = tpts(1:ninv)
          DO i1 = 1, ninv
            IF (ydat(i1).EQ.0.0) ydat(i1) = LO
          ENDDO
          CALL GRTRAC(xdat,ydat,ninv,'inv ','LINE',1)
          CALL DrawFrame
        
          GOTO 10
c...      ???
          map1x = 0.05             ! Add sensitivity to nxbin, nybin...
          map2x = map1x + 0.80 
          map1y = map2y
          map2y = map1y + 0.15
          cxmin = 1.0
          cxmax = REAL(ninv)
          cymin =  HI
          cymax = -HI
          DO i1 = 1, ninv
            cymin = MIN(cymin,qpts(i1,2))
            cymax = MAX(cymax,qpts(i1,2))
          ENDDO
          xlabel = 'none     '
          ylabel = '???      '
          CALL GRTSET_TRIM(' ',' ',' ',' ',' ',      ! TITLE,REF,nVIEW,PLANE,glabel,
     .                     cxmin,cxmax,cymin,cymax,
     .                     ' ',xlabel,ylabel,        ! TABLE,XLAB,YLAB,
     .                     0,' ',0,' ',1)            ! 0,smooth,0,ANLY,1)
          CALL GRTRAC(xdat,qpts(1,2),ninv,'inv ','LINE',1)
          CALL DrawFrame
c...      ???
          map1x = 0.05             ! Add sensitivity to nxbin, nybin...
          map2x = map1x + 0.80 
          map1y = map2y
          map2y = map1y + 0.15
          cxmin = 1.0
          cxmax = REAL(ninv)
          cymin =  HI
          cymax = -HI
          DO i1 = 1, ninv
            cymin = MIN(cymin,qpts(i1,2))
            cymax = MAX(cymax,qpts(i1,2))
          ENDDO
          xlabel = 'none     '
          ylabel = '???      '
          CALL GRTSET_TRIM(' ',' ',' ',' ',' ',      ! TITLE,REF,nVIEW,PLANE,glabel,
     .                     cxmin,cxmax,cymin,cymax,
     .                     ' ',xlabel,ylabel,        ! TABLE,XLAB,YLAB,
     .                     0,' ',0,' ',1)            ! 0,smooth,0,ANLY,1)
          CALL GRTRAC(xdat,qpts(1,2),ninv,'inv ','LINE',1)
          CALL DrawFrame
 10       CONTINUE

        ENDIF

        CALL Frame
      ENDIF



c...  Presentation plots:
      IF (.TRUE..AND.ninv.GT.0.AND.plot_option) THEN
c        WRITE(0,*) 'INVERTED PROFILE B'


        IF (.TRUE..AND.(opt%nxbin.GT.1.OR.opt%nybin.GT.1)) THEN
c...      Plot camera image:

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
          CALL DrawColourScale(2,2,qmin,qmax,'none')
          CALL DrawFrame
        ENDIF


        CALL LoadTriangles

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

          IF (qmin.GT.0.0) THEN

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

c...        Background:
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
c...        Annotate:
c            CALL LinCol(ncols+1)
c            CALL CTRMAG(12)
c            WRITE(caption,'(A,3I4)') 'NX,NY,BIN:',opt%nxbin,opt%nybin,nbin
c            CALL PLOTST(0.02,0.02,caption(1:LEN_TRIM(caption)))
            CALL DrawGrid(95)
            CALL DrawColourScale(2,2,qmin,qmax,'none')
            CALL DrawFrame
          ENDIF
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

c...    Background:
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
c        DO i1 = 101, 200
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
        CALL DrawGrid(95)
        Call DrawFrame


        IF (.FALSE.) THEN
c...      Reference profile, plasma grid:

          CALL Frame

c...      Comment this out to have scale comes from inverted image:
          qmin =  0.0
          qmax = -HI
          DO ir = 2, nrs
            IF (idring(ir).EQ.BOUNDARY) CYCLE
            DO ik = 1, nks(ir)
              IF (rs(ik,ir).LT.cxmin.OR.rs(ik,ir).GT.cxmax.OR.
     .            zs(ik,ir).LT.cymin.OR.zs(ik,ir).GT.cymax) CYCLE
              qmax = MAX(qmax,pinalpha(ik,ir))
            ENDDO
          ENDDO

          WRITE(0,*) 'PROFILE QMIN,QMAX 1 :',qmin,qmax
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
          DO ir = 2, nrs
            IF (idring(ir).EQ.BOUNDARY) CYCLE
            ike = nks(ir)
            IF (ir.LT.irsep) ike = nks(ir) - 1
            DO ik = 1, ike
              id = korpg(ik,ir)
              p1(1:4) = rvertp(1:4,id)
              p2(1:4) = zvertp(1:4,id)          
              CALL SetCol255_04(2,pinalpha(ik,ir),qmin,qmax)
              CALL FILCOL(255)
              CALL LINCOL(255) 
              CALL PTPLOT(p1,p2,1,4,1)
            ENDDO
          ENDDO

c...      Annotate:
c          CALL LinCol(ncols+1)
c          CALL CTRMAG(12)
c          WRITE(caption,'(A,3I4)') 'NX,NY,BIN:',opt%nxbin,opt%nybin,nbin
c          CALL PLOTST(0.02,0.02,caption(1:LEN_TRIM(caption)))
          CALL DrawGrid(95)
          CALL DrawColourScale(2,2,qmin,qmax,'none')
          CALL DrawFrame

c...      Dump data for processing in IDL:
c          file = filename(1:LEN_TRIM(filename))//'osm'
c          WRITE(6,*) '989: Duming OSM data to interface file'
c          WRITE(6,*) '     FILE = >',file,'<'
c          CALL inOpenInterface(file,ITF_WRITE)
c         Moved to SLoutplot.f        

        ENDIF


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
          CALL DrawGrid(-95)
          Call DrawFrame
        ENDIF

        CALL Frame

c...    He ratio stuff I think:
        IF (.FALSE.) THEN
          map1x = 0.05             ! Add sensitivity to nxbin, nybin...
          map2x = map1x + 0.80 
          map1y = 0.05
          map2y = map1y + 0.15

          ndat = 0
          DO i1 = 1, ninv
            xcen = 0.0
            ycen = 0.0
            DO i2 = 1, npts(i1)
              xcen = xcen + xpts(i2,i1)
              ycen = ycen + ypts(i2,i1)
            ENDDO
            xcen = xcen / REAL(npts(i1))
            ycen = ycen / REAL(npts(i1))

            IF (ycen.GT.1.60.AND.ycen.LE.1.61) THEN
c            IF (ycen.GT.1.50.AND.ycen.LE.1.51) THEN
c            IF (ycen.GT.1.40.AND.ycen.LE.1.41) THEN
              ndat = ndat + 1
              xdat(ndat) = xcen
              ydat(ndat) = qpts(i1,1)
              WRITE(0,*) 'DATA:',i1,xcen,ycen,qpts(i1,2)
            ENDIF
          ENDDO

          cxmin = xdat(1)
          cxmax = xdat(ndat)
          cymin = 0.0
          cymax = 0.0
          DO i1 = 1, ndat
            cymin = MIN(cymin,ydat(i1))
            cymax = MAX(cymax,ydat(i1))
          ENDDO
          xlabel = 'pixel '
          ylabel = '???      '
          CALL GRTSET_TRIM(' ',' ',' ',' ',' ',      ! TITLE,REF,nVIEW,PLANE,glabel,
     .                     cxmin,cxmax,cymin,cymax,
     .                     ' ',xlabel,ylabel,        ! TABLE,XLAB,YLAB,
     .                     0,' ',0,' ',1)            ! 0,smooth,0,ANLY,1)
          CALL GRTRAC(xdat,ydat,ndat,'inv ','LINE',1)
          CALL DrawFrame


          map1x = 0.05             ! Add sensitivity to nxbin, nybin...
          map2x = map1x + 0.80 
          map1y = map1y + 0.15
          map2y = map1y + 0.15


          ndat = 0
          DO i1 = 1, ninv
            xcen = 0.0
            ycen = 0.0
            DO i2 = 1, npts(i1)
              xcen = xcen + xpts(i2,i1)
              ycen = ycen + ypts(i2,i1)
            ENDDO
            xcen = xcen / REAL(npts(i1))
            ycen = ycen / REAL(npts(i1))

            IF (ycen.GT.1.60.AND.ycen.LE.1.61) THEN
c            IF (ycen.GT.1.50.AND.ycen.LE.1.51) THEN
c            IF (ycen.GT.1.40.AND.ycen.LE.1.41) THEN
              ndat = ndat + 1
              xdat(ndat) = xcen
              ydat(ndat) = qpts(i1,2)
              WRITE(0,*) 'DATA:',i1,xcen,ycen,qpts(i1,2)
            ENDIF
          ENDDO

          cxmin = xdat(1)
          cxmax = xdat(ndat)
          cymin = 0.0
          cymax = 0.0
          DO i1 = 1, ndat
            cymin = MIN(cymin,ydat(i1))
            cymax = MAX(cymax,ydat(i1))
          ENDDO
          xlabel = 'pixel '
          ylabel = '???      '
          CALL GRTSET_TRIM(' ',' ',' ',' ',' ',      ! TITLE,REF,nVIEW,PLANE,glabel,
     .                     cxmin,cxmax,cymin,cymax,
     .                     ' ',xlabel,ylabel,        ! TABLE,XLAB,YLAB,
     .                     0,' ',0,' ',1)            ! 0,smooth,0,ANLY,1)
          CALL GRTRAC(xdat,ydat,ndat,'inv ','LINE',1)
          CALL DrawFrame

          CALL Frame
        ENDIF

      ENDIF




      IF (nbin.GT.1) THEN
        nxbin = opt%nxbin
        nybin = opt%nybin
        maxpixel = DBLE(MAXVAL(image))
        WRITE(0,*) 'SAVING JPEG'
        CALL Save_Image(image,nxbin,nybin,0.0,maxpixel,0.0D0)  ! Generate a jpeg file 
        WRITE(0,*) 'DONE'
      ENDIF

      IF (ALLOCATED(image)) DEALLOCATE(image)

c...  Clear triangle arrays:
      CALL DEALLOC_VERTEX  
      CALL DEALLOC_SURFACE   
      CALL DEALLOC_TRIANGLE

      WRITE(0,*) '  PLOTS DONE'




      RETURN
 97   CALL ER('Output989','File access problem',*99)
 98   CALL ER('Output989','Trouble opening .raw.idl',*99)
 99   WRITE(0,*) '  FILE=',file(1:LEN_TRIM(file))
      STOP
      END
