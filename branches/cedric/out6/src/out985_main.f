c     -*-Fortran-*-
c
c ====================================================================== 
c
c subroutine: InterpolateValue
c
c
      SUBROUTINE InterpolateValue(nbfield,bfield,type,p,shape,region,
     .                            qval,rval,zval)
      IMPLICIT none
 
      INTEGER nbfield,type,shape,region
      REAL*8  bfield(nbfield,6),p(6),rval,zval,qval,deltar,deltaz,qmin

      INTEGER, PARAMETER :: MAXPTS = 100

      INTEGER ndim,npts,ipts
      LOGICAL cont
      REAL*8  r(MAXPTS),z(MAXPTS),q(MAXPTS,2),rmin,rmax,zmin,zmax,
     .        r1,z1,r2,z2

c...  Shape = 1 - rectangle
c             2 - circle
c
c     Region = 1 - perimeter
c              2 - volume
c
      SELECTCASE (shape)
        CASE (1)  ! Rectangle
          r1 = p(1)
          z1 = p(2)
          deltar = p(3)
          deltaz = p(4)
          rmin = r1 - 0.5D0 * deltar
          rmax = r1 + 0.5D0 * deltar
          zmin = z1 - 0.5D0 * deltaz
          zmax = z1 + 0.5D0 * deltaz
          ndim = INT(SQRT(REAL(MAXPTS)))
        CASE (2)  ! Circle
        CASE DEFAULT
          CALL ER('InterpolateValue','Unknown shape',*99)
      ENDSELECT

      cont = .TRUE.

      DO WHILE (cont) 

c...    Setup list of R,Z points to interpolate to:
        SELECTCASE (shape)
c...      Rectangle:
          CASE (1)  
            rmin = MAX(rmin,r1-0.5D0*deltar)
            rmax = MIN(rmax,r1+0.5D0*deltar)
            zmin = MAX(zmin,z1-0.5D0*deltaz)
            zmax = MIN(zmax,z1+0.5D0*deltaz)

            SELECTCASE (region)
              CASE (1)  ! Volume
                npts = 0
                DO z2 = zmin, zmax, (zmax - zmin) / DBLE(REAL(ndim))
                  DO r2 = rmin, rmax, (rmax - rmin) / DBLE(REAL(ndim))
                    npts = npts + 1
                    WRITE(0,*) 'npts:',npts,MAXPTS
                    IF (npts.GT.MAXPTS) 
     .                CALL ER('InterpolateValue','MAXPTS exceeded',*99)
                    r(npts) = r2
                    z(npts) = z2
                  ENDDO
                ENDDO

              CASE (2)  ! Perimeter
              CASE DEFAULT
                CALL ER('InterpolateValue','Unknown region',*99)
            ENDSELECT
c...      Circle:
          CASE (2)  ! Circle
          CASE DEFAULT
            CALL ER('InterpolateValue','Unknown shape',*99)
        ENDSELECT

c...    Interpolate and derived quantity of interest (if necessary):
        SELECTCASE (type)
          CASE (1) ! Poloidal field
            CALL Interpolate
     .             (nbfield,bfield(1,1),bfield(1,2),bfield(1,4),
     .              npts,r,z,q(1,1),0)
            CALL Interpolate
     .             (nbfield,bfield(1,1),bfield(1,2),bfield(1,6),
     .              npts,r,z,q(1,2),0)
            DO ipts = 1, npts                                     ! Calcualte this in IDL instead..?
              q(ipts,1) = DSQRT(q(ipts,1)**2 + q(ipts,2)**2) 
            ENDDO

          CASE (2) ! PSIn
          CASE DEFAULT
            CALL ER('InterpolateValue','Unknown type',*99)
        ENDSELECT
 
c...    Adjust search parameters or exit loop:
        SELECTCASE (shape)
c...      Rectangle:
          CASE (1)  
            qmin = 1.0E+30
            DO ipts = 1, npts
              IF (DABS(q(ipts,1)-qval).LT.qmin) THEN
                qmin = DABS(q(ipts,1)-qval)
                r1   = r(ipts)
                z1   = z(ipts)
              ENDIF
            ENDDO
            deltar = 0.5D0 * deltar
            deltaz = 0.5D0 * deltaz         
c...      Circle:
          CASE (2)  
        ENDSELECT


      ENDDO

      RETURN
 99   STOP
      END
c
c ====================================================================== 
c
c subroutine: FindPoloidalNulls
c
c
      SUBROUTINE FindPoloidalNulls(nbfield,bfield,MAXXPTS,nxpts,xpts)
      IMPLICIT none
      
      INTEGER nbfield,MAXXPTS,nxpts
      REAL*8  bfield(nbfield,6),xpts(MAXXPTS,2)

      INTEGER i1,i2,ir,iz,ixpts(MAXXPTS,2)
      LOGICAL cont
      REAL*8  rmin,rmax,zmin,zmax,rmin1,rmax1,zmin1,zmin2,r,z,val,
     .        rxpts(MAXXPTS),zxpts(MAXXPTS),vxpts(MAXXPTS)

      nxpts = 0
      ixpts = 0
      vxpts = 1.0D+20 

c...  Find spatial extent of b-field data:
      rmin =  1.0D+20
      rmax = -1.0D+20
      zmin =  1.0D+20
      zmax = -1.0D+20
      DO i1 = 1, nbfield
        rmin = MIN(rmin,bfield(i1,1))
        rmax = MAX(rmax,bfield(i1,1))
        zmin = MIN(zmin,bfield(i1,2))
        zmax = MAX(zmax,bfield(i1,2))
      ENDDO

      IF (.TRUE.) THEN
        cont = .TRUE.
        DO WHILE (cont) 
c..
          DO iz = 1, 10
            DO ir = 1, 10
c              CALL FindMinPoloidalField
c     .               (nbfield,bfield,rmin1,rmax1,zmin1,zmax1,r,z,val)
c...
              DO i1 = 1, MAX(1,nxpts)
                IF (val.LT.vxpts(i1)) THEN
                  nxpts = MIN(nxpts+1,MAXXPTS)
c...              Make room:              
                  DO i2 = nxpts, i1+1
                    vxpts(i2) = vxpts(i2-1)

                  ENDDO
c...
                  vxpts(i1) = val
                  rxpts(i1) = r
                  zxpts(i1) = z
                  EXIT
                ENDIF
              ENDDO

            ENDDO
          ENDDO

          cont = .FALSE.
        ENDDO

      ENDIF

      RETURN
 99   STOP
      END

c
c
c ====================================================================== 
c
c subroutine: TraceMagneticFieldLines
c
c
      SUBROUTINE TraceMagneticFieldLines
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none


      INTEGER   fp,i1,i2,iline,nbfield
      REAL*8    r,z,psin
      CHARACTER file*512,dummy*1024

      REAL*8,ALLOCATABLE :: bfield(:,:)

c...  Load magnetic field data:
      fp = 99
      file = 'bfield.dat'
      OPEN(fp,FILE=file(1:LEN_TRIM(file)),
     .     FORM='FORMATTED',STATUS='OLD',ERR=98)
      DO WHILE(.TRUE.)
        READ(fp,'(A1024)',END=10) dummy
        SELECTCASE (dummy(1:6))
          CASE('[HEAD]')               ! Check version...  add EFIT run data here...
          CASE('[SHOT]')
          CASE('[TIME]')
          CASE('[NPTS]')
            READ(dummy(7:),*) nbfield            
            ALLOCATE(bfield(nbfield,6))
          CASE('[DATA]')
            IF (.NOT.ALLOCATED(bfield))
     .        CALL ER('TraceMag...','[NPTS] tag missing',*99)
            DO i1 = 1, nbfield
              READ(fp,*,END=97,ERR=96) (bfield(i1,i2),i2=1,6)  ! List what the components are...
            ENDDO
            EXIT  ! Perhaps in future cycle here so that b-field data 
                  ! can be loaded either side of the requested time slice and
                  ! interpolated...
          CASE DEFAULT
            CALL ER('TraceMagnetic...','Unknown bfield file tag',*99)
        ENDSELECT
      ENDDO
 10   CONTINUE
      CLOSE(fp)

c...  Find x-points:
c     (only needed really if specifying that the field line to be traced is
c      in the PFZ)


      
      IF (.TRUE.) THEN
        DO iline = 1, 1

c...      Choose starting point of field line to be traced:
          IF (.TRUE.) THEN
c...        R,Z specified:
            r = 1.365D0
            z = 0.000D0
c...        Find PSIn:
            CALL Interpolate
     .             (nbfield,bfield(1,1),bfield(1,2),bfield(1,3),
     .              1,r,z,psin,0)
            WRITE(0,*) 'R,Z,PSIN:',r,z,psin
            STOP 'DONE FOR NOW'
          ELSE
c...        PSIn,region specified:

          ENDIF



        
        ENDDO
      ENDIF
c...  Clear arrays:
      IF (ALLOCATED(bfield)) DEALLOCATE(bfield)

      RETURN
 96   CALL ER('TraceMag...','Problem reading b-field data file',*99)
 97   CALL ER('TraceMag...','Unexpected end of file',*99)
 98   CALL ER('TraceMag...','b-field data file not found',*99)
 99   WRITE(0,*) 'FILE:',file(1:LEN_TRIM(file))   
      STOP
      END
c
c ====================================================================== 
c
c subroutine: Interpolate
c
c
      SUBROUTINE Interpolate(n,x,y,f,n1,x1,y1,f1,mode)
      USE qshep2d
      IMPLICIT NONE

      INTEGER :: n,n1,mode
      REAL*8  :: x(n),y(n),f(n),x1(n1),y1(n1),f1(n1)

      INTEGER :: i, ier, j, k, lnext(n), nr
      INTEGER, ALLOCATABLE :: lcell(:,:)
      REAL*8  :: a(5,n), dx, dy, eps, eq, eqx, eqy,
     .           p(10), px, py, q, q1, 
     .           qx, qy, rmax, rq, rsq(n), xmin, yk, 
     .           ymin
      LOGICAL warning

c      REAL fx,fy,fq

! QSHEP2 PARAMETERS AND LOGICAL UNIT FOR OUTPUT

c      INTEGER, PARAMETER  :: lout = 6, n = 36, nq = 13, nr = 3, nw = 19
      INTEGER, PARAMETER  :: lout = 6, nq = 13, nw = 19

      warning = .FALSE.

      nr = INT(SQRT(REAL(n)/3.0))

      ALLOCATE(lcell(nr,nr))  ! Do this every call...?

! COMPUTE PARAMETERS DEFINING THE INTERPOLANT Q.
c      WRITE(0,*) 'INTERPOL:'
c      WRITE(0,*) n,nr
c      WRITE(0,*) 'x:',x
c      WRITE(0,*) 'y:',y
c      WRITE(0,*) 'f:',f

      CALL qshep2(n,x,y,f,nq,nw,nr,lcell,lnext,xmin,ymin,
     .            dx,dy,rmax,rsq,a,ier)
      IF (ier.NE.0) THEN
        CALL ER('985','Error in call to QSHEP2',*99)
      ELSE
! COMPUTE THE MACHINE PRECISION EPS.
        eps = EPSILON(1.0)

! COMPUTE INTERPOLATION ERRORS AND TEST FOR AGREEMENT IN THE
!   Q VALUES RETURNED BY QS2VAL AND QS2GRD.
        DO  i = 1, n1
          px = x1(i)
          py = y1(i)
          q1 = qs2val(px,py,n,x,y,f,nr,lcell,lnext,xmin,ymin,
     .                dx,dy,rmax,rsq,a)
          CALL qs2grd(px,py,n,x,y,f,nr,lcell,lnext,xmin,ymin,
     .                dx,dy,rmax,rsq,a,q,qx,qy,ier)
          IF     (ier.NE.0) THEN
            IF (ier.EQ.2) THEN
              WRITE(0,*) 'px,y:',px,py
              IF (.NOT.warning) THEN
                WRITE(0,*) 'WARNING Interpolate: Point(s) outside grid'
                warning = .TRUE.
              ENDIF
            ELSE
              CALL ER('985','Error detected in call to INTERPOL',*99)
            ENDIF
          ELSEIF (ABS(q1-q).GT.3.0*ABS(q)*eps) THEN
            CALL ER('985','Problem in INTERPOL',*99)
          ELSE
            f1(i) = q1
          ENDIF
        ENDDO
      ENDIF

      IF (ALLOCATED(lcell)) DEALLOCATE(lcell)

      RETURN
 99   WRITE(0,*) '  ERROR CODE, IERR=',ier
      STOP
      END
c
c ======================================================================
c
      SUBROUTINE BinImage(nxbin1,nybin1,image,nxbin2,nybin2)
      IMPLICIT none

      INTEGER nxbin1,nybin1,nxbin2,nybin2
      REAL*8  image(1100,1100)

      INTEGER ix,iy,ii
      LOGICAL blankpixel
      REAL    x,x1,step,y,y1,frac
      REAL*8  val


c...  Horizontal binning:
      IF (nxbin1.NE.nxbin2) THEN
        step = REAL(nxbin1) / REAL(nxbin2)
 
        WRITE(0,*) 'step?',step,nxbin1,nxbin2
        WRITE(0,*) '    ?',step,nybin1,nybin2

        DO iy = 1, nybin1
          x = 0.0
          DO ix = 1, nxbin2
            val = 0.0D0
            blankpixel = .FALSE.
            DO x1 = x, x+step-0.00001, 1.0 
              ii = INT(x1) + 1
              IF (x1.EQ.x) THEN
                frac = 1.0 - MOD(x,1.0)
              ELSEIF (x1+1.0.GT.x+step) THEN
                frac = MOD(x+step,1.0)
              ELSE
                frac = 1.0
              ENDIF
              IF (image(ii,iy).EQ.0.0D0) blankpixel = .TRUE.
              val = val + DBLE(frac) * image(ii,iy)
c              IF (frac.GT.0.0) val = val + DBLE(frac) * image(ii,iy)
c              WRITE(0,'(A,3I7,3F10.4)') 'BIN:',iy,ix,ii,frac,
c     .         MOD(x1,1.0),x1
            ENDDO
            IF (blankpixel) val = 0.0D0
            image(ix,iy) = val
            x = x + step
          ENDDO
        ENDDO
      ENDIF


c...  Vertical binning:
      IF (nybin1.NE.nybin2) THEN
        step = REAL(nybin1) / REAL(nybin2)
 
        WRITE(0,*) 'y step?',step,nxbin1,nxbin2
        WRITE(0,*) '      ?',step,nybin1,nybin2

        DO ix = 1, nxbin2
          y = 0.0
          DO iy = 1, nybin2
            val = 0.0D0
            blankpixel = .FALSE.
            DO y1 = y, y+step-0.00001, 1.0 
              ii = INT(y1) + 1
              IF (y1.EQ.y) THEN
                frac = 1.0 - MOD(y,1.0)
              ELSEIF (y1+1.0.GT.y+step) THEN
                frac = MOD(y+step,1.0)
              ELSE
                frac = 1.0
              ENDIF
              IF (image(ix,ii).EQ.0.0D0) blankpixel = .TRUE.
              val = val + DBLE(frac) * image(ix,ii)
c              image(ix,iy) = image(ix,iy) + DBLE(frac) * image(ix,ii)
c              WRITE(0,'(A,3I6,3F10.4)') 'BIN:',iy,ix,ii,frac,MOD(y1,1.0)
c     .                                  ,y1
            ENDDO
            IF (blankpixel) val = 0.0D0
            image(ix,iy) = val
            y = y + step
          ENDDO
        ENDDO
      ENDIF  


c...  Scale:
      val = DBLE(REAL(nxbin1) / REAL(nxbin2)) *
     .      DBLE(REAL(nybin1) / REAL(nybin2))    
      WRITE(0,*) 'BINNING IMAGE, CORRECTING:',val
      image = image / val

      nxbin1 = nxbin2
      nybin1 = nybin2

c      STOP 'CRAZY'

      RETURN
 99   STOP
      END
c
c
c ======================================================================
c
      SUBROUTINE GetMemoryI4(size,ptr)
c      SUBROUTINE GetMemoryI4(size,ptr,bytes)
      IMPLICIT none
c...  Input:
      INTEGER, POINTER :: ptr(:)
      INTEGER          :: size,bytes


      INTEGER, POINTER :: oldptr(:)
      INTEGER, TARGET, ALLOCATABLE :: targ(:)

      WRITE(0,*) 'HERE!'
c      WRITE(0,*) 'A:',ptr(1),ptr(100) 
c      oldptr => ptr
c      WRITE(0,*) 'A:',oldptr(1),oldptr(100) 
c      ALLOCATE(targ(bytes))
c      WRITE(0,*) 'A: done' 
c      targ(1:size) = ptr(1:size)
c      WRITE(0,*) 'A: done 2' 
c      ptr => targ
c      WRITE(0,*) 'A: done 3'        
c      DEALLOCATE(oldptr)
c      WRITE(0,*) 'A: done 4'        

      RETURN
 99   STOP
      END
c
c ====================================================================== 
c
      SUBROUTINE ProcessPixels(npixel,pixel)
      USE mod_out985
      USE mod_out985_variables
      USE mod_interface
      IMPLICIT none
c...  Input:
      INTEGER npixel
      TYPE(type_view) :: pixel(npixel)

      REAL Clock2

      INTEGER   i1,i2,i3,i4,status,fp,fp2,m,n,codec,iobj,ierr,volcnt,
     .          pixcnt,ndum1,npts,vcount
      LOGICAL   ascii,message
      REAL      rtime,ttime
      CHARACTER file*1024

      INTEGER, ALLOCATABLE :: idum1(:)

      REAL*8, ALLOCATABLE, TARGET :: ddum1(:)
      REAL*4, ALLOCATABLE, TARGET :: mapchk(:)
      REAL*4, ALLOCATABLE, TARGET :: rdum1(:)

      message = .TRUE.

      nchord = 0


      IF (.TRUE.) THEN
c...    Open a stream for storing the inversion matrix 'A':
        IF (.FALSE..OR.
     .      opt%det_nxbin(opt%ndet)*opt%det_nybin(opt%ndet).LT.20) THEN
c...      Nothing fancy for small maps:
          codec = 1
          ascii = .TRUE.
        ELSE
c...      Large array coming:
          codec = 2 
          ascii = .NOT..FALSE.
        ENDIF
        m = npixel   ! Row rank in A
        n = 0        ! Column rank in A  (???)
        DO iobj = 1, nobj   ! *** need to reorder objects to that first n are integration volumes? ***
          IF (obj(iobj)%type.EQ.OP_INTEGRATION_VOLUME) 
     .      n = MAX(n,obj(iobj)%ivolume)
        ENDDO
        opt%n = n
        WRITE(0,*) 'NUMBER OF INVERSION MESH PIXELS:',n

        fp = 98
        file = 'output.'//opt%fmap(1:LEN_TRIM(opt%fmap))//'.ray.map'
        WRITE(0,*) 'FMAP:',file(1:LEN_TRIM(file))

        IF (ascii) THEN 
          OPEN(fp,FILE=file(1:LEN_TRIM(file)),
     .         FORM='FORMATTED',STATUS='REPLACE',ERR=98)
          WRITE(fp,'(A   )') '* Header'
          WRITE(fp,'(A14,F5.1)') 'VERSION=',1.0
          WRITE(fp,'(A14,I7)'  ) 'CODEC=  ',codec
          WRITE(fp,'(A14,I7)'  ) 'M=      ',m
          WRITE(fp,'(A14,I7)'  ) 'N=      ',n
          WRITE(fp,'(A14,I7)'  ) 'NDET=   ',opt%ndet
          WRITE(fp,'(A14,99I7)') 'XBIN=   ',opt%det_nxbin(1:opt%ndet)
          WRITE(fp,'(A14,99I7)') 'YBIN=   ',opt%det_nybin(1:opt%ndet)
        ELSE
          OPEN(fp,FILE=file(1:LEN_TRIM(file)),
     .         FORM='UNFORMATTED',STATUS='REPLACE',ERR=98)
          WRITE(fp) 888888888
          WRITE(fp) 1.0
          WRITE(fp) codec
          WRITE(fp) m,n,opt%ndet,opt%det_nxbin(1:opt%ndet),
     .                           opt%det_nybin(1:opt%ndet)
        ENDIF

c...    Simple format output:
        fp2 = 97
        IF (fp2.NE.0) THEN
          file = 'output.'//opt%fmap(1:LEN_TRIM(opt%fmap))//'.ray.A'
          WRITE(0,*) 'FMAP:',file(1:LEN_TRIM(file))
          OPEN(fp2,FILE=TRIM(file),FORM='FORMATTED',STATUS='REPLACE',
     .         ERR=98)
          WRITE(fp2,'(2I7)') m,n

c...      Also generate mod_interface file:
c          WRITE(file,'(1024X)')          
c          file = 'output.'//TRIM(opt%fmap)//'.idl.A'
c          CALL inOpenInterface(TRIM(file))
c          CALL inPutData(m,'m','none')
c          CALL inPutData(n,'n','none')
        ENDIF

      ELSE
        m = npixel
        n = nobj
      ENDIF

c...
      ALLOCATE(mapchk(n))
      ALLOCATE(idum1(n))

      mapchk = 0.0

c...  Initialize arrays used to store view intersection data:  *** SHOULD BE MOVED TO INSIDE MPI LOOP? ***
      nvwinter = 0
      ngbinter = 0
      nobinter = 0
      ALLOCATE (vwinter(MAXINTER))
      ALLOCATE (gbinter(MAXINTER))
      ALLOCATE (obinter(MAXINTER))

c...  Assemble a list of wall surfaces to speed things up:  (30% faster for current setup...)
      WRITE(0,*) '  BUILDING WALL SURFACE LISTS'  
      nvwlist = 0
      ngblist = 0
      DO i1 = 1, nobj
        DO i2 = 1, MAX(obj(i1)%nsur,obj(i1)%nside)
c...      Check if surface is part of the vessel wall:
          IF (obj(i1)%tsur(i2).EQ.SP_VESSEL_WALL  ) nvwlist = nvwlist+1
          IF (obj(i1)%tsur(i2).EQ.SP_GRID_BOUNDARY) ngblist = ngblist+1
        ENDDO
      ENDDO
      ALLOCATE (vwlist(nvwlist,2))
      ALLOCATE (gblist(ngblist,2))
      i3 = 0
      i4 = 0
      DO i1 = 1, nobj
        DO i2 = 1, MAX(obj(i1)%nsur,obj(i1)%nside)
c...      Check if surface is part of the vessel wall:
          IF     (obj(i1)%tsur(i2).EQ.SP_VESSEL_WALL) THEN
            i3 = i3 + 1
            vwlist(i3,1) = i1
            vwlist(i3,2) = i2
          ELSEIF (obj(i1)%tsur(i2).EQ.SP_GRID_BOUNDARY) THEN
            i4 = i4 + 1
            gblist(i4,1) = i1
            gblist(i4,2) = i2
          ENDIF
        ENDDO
      ENDDO
      IF (i3.NE.nvwlist) 
     .  CALL ER('ProcessPixels','Something isn''t right',*99)
      IF (i4.NE.ngblist) 
     .  CALL ER('ProcessPixels','Something isn''t right',*99)
      WRITE(0,*) '  DONE',nvwlist,ngblist


      ALLOCATE(ddum1(n))  ! MPI problem?  nobj=m, should be # integration volumes
      ALLOCATE(rdum1(100))  ! MPI problem?  nobj=m, should be # integration volumes

      rtime = 0.0
      pixcnt = m



c...  MPI!  (for MPI to work the various lists -- vw,gb,ob -- need to be passed like the objects are passed, I think)



      vcount = 0

c...  Loop over pixels:
c      DO i1 = 200000, npixel
      DO i1 = 1, npixel

c        IF (MOD(i1,npixel/10).EQ.0) 
        IF (MOD(i1,1000).EQ.0) 
     .    WRITE(0,'(A,I7,A,I7)') 'PROCESSING PIXEL ',i1,' OF ',npixel

        status = 0

        ddum1 = 0.0D0
        pixel(i1)%track => ddum1  ! Dont' make part of view derived type definition, just send on its own?

        rdum1 = 0.0
c        pixel(i1)%spectrum => rdum1  ! also a problem...

        rtime = Clock2()

        IF (pixel(i1)%valid) 
     .    CALL IntegrateAlongChords(pixel(i1),status)   ! new name!
c     .    CALL IntegrateAlongChords(opt,pixel(i1),nobj,obj,status)   ! new name!

        ttime = ttime + (Clock2() - rtime)

        IF (pixel(i1)%valid) vcount = vcount + 1

        IF     (status.EQ.-1) THEN
          WRITE(0,*) '  STATUS.EQ.-1 FOR PIXEL, CONTINUING',i1,nchord
          status = 0
        ELSEIF (status.LT.-1) THEN
          WRITE(0,*) '  STATUS.LT.-1 FOR PIXEL, STOPPING',i1,nchord
          RETURN
        ENDIF

c...    Check if pixel sampled the inversion mesh: 
        IF (SUM(pixel(i1)%track).EQ.0.0D0) pixcnt = pixcnt - 1

c...    Check that all cells on the inverson mesh are sampled:
        mapchk(1:n) = mapchk(1:n) + SNGL(pixel(i1)%track(1:n))

c...    Store map data in inversion matrix A:   ! MPI PROBLEM -- THESE MAY COME BACK OUT OF ORDER 
        IF (.TRUE.) THEN
c...      Option to check: don't store rows for pixels that do not intersect the invesion mesh... 

          IF (fp2.NE.0) THEN
            DO i2 = 1, n
              IF (pixel(i1)%track(i2).NE.0.0D0) THEN
                WRITE(fp2,'(2I7,1P,E22.12)') i1,i2,pixel(i1)%track(i2)
c                CALL inPutData(i1,'i','none')
c                CALL inPutData(i2,'j','none')
c                CALL inPutData(pixel(i1)%track(i2),'data','m')
              ENDIF
            ENDDO
          ENDIF

          IF     (codec.EQ.1) THEN
            DO i2 = 1, n, 5  ! Better way to write this with a format statement..? 
              WRITE(fp,'(10D20.12)') 
     .          (pixel(i1)%track(i3),i3=i2,MIN(i2+4,n))
            ENDDO
          ELSEIF (codec.EQ.2) THEN
c...        Add some brains:            
            ndum1 = 0
            DO i2 = 1, n  
              IF (pixel(i1)%track(i2).NE.0.0D0) THEN 
                ndum1 = ndum1 + 1
                idum1(ndum1) = i2
              ENDIF  
            ENDDO
            IF (ndum1.GT.0) THEN
              IF (ascii) THEN
                WRITE(fp,'(A6,2I8)') 'pixel:',i1,ndum1
                WRITE(fp,'(4(I7,1P,D22.15))') 
     .            (idum1(i2),pixel(i1)%track(idum1(i2)),i2=1,ndum1)
              ELSE
                WRITE(fp) 777777777
                WRITE(fp) i1,ndum1,    
     .            (idum1(i2),pixel(i1)%track(idum1(i2)),i2=1,ndum1)
              ENDIF
            ENDIF  
          ELSE
            CALL ER('ProcessPixels','Unknown codec',*99)
          ENDIF
        ELSE
        ENDIF

      ENDDO  ! Pixel loop

c...  Confirm that all volumes on the inversion mesh were sampled by the
c     LOS integration:
      volcnt = n
      DO i1 = 1, n
        IF (mapchk(i1).EQ.0.0) volcnt = volcnt - 1
      ENDDO     
      WRITE(0,'(F8.1,A,I7,A,I7,A)') 
     .  REAL(pixcnt)/REAL(m)*100.0,'% OF DETECTOR VIEWS '//
     .  'SAMPLE THE INVERSION MESH  (',pixcnt,' of',m,')'
      WRITE(0,'(F8.1,A,I7,A,I7,A)') 
     .  REAL(volcnt)/REAL(n)*100.0,'% OF INVERSION MESH '//
     .  'SAMPLED BY LOS INTEGRATION (',volcnt,' of',n,')'
      WRITE(0,'(I7,A,I7,A)')
     .  vcount,' OF ',npixel,' PIXELS VALID (NOT MASKED)'

c...  MPI!
      WRITE(0,'(A,F10.2,A)') 'TOTAL TIME FOR THIS RUN:',ttime,' s'
      WRITE(0,'(A,F10.3,A)') 'AVERAGE TIME PER PIXEL :',
     .                       ttime/REAL(npixel),' s'
      WRITE(0,'(A,F10.1,A)') 'TIME FOR 1E6 PIXELS    :',
     .                       ttime/REAL(npixel)*1E+06/3600.0,' h'


      IF (fp2.NE.0) THEN
        CLOSE(fp2)
c        CALL inCloseInterface
      ENDIF

      IF (.TRUE.) THEN
c...    Write poloidal cross-section information for inversion grid at
c       the end of .map file:

        IF (ascii) THEN
          WRITE(fp,'(A)') '* Inversion grid poloidal cross-section'
          WRITE(fp,'(I10)') n
        ELSE
          WRITE(fp) 999999999  ! Marker for start of inversion mesh grid data
          WRITE(fp) n
        ENDIF
        DO iobj = 1, n
          obj(iobj)%path = mapchk(iobj)
          npts = 0
          IF (obj(iobj)%nsur.EQ.4) npts = 4  ! Toroidally continuous   ! ***LAME/WEAK***
          IF (obj(iobj)%nsur.EQ.5) npts = 3  ! Triangle             
          IF (obj(iobj)%nsur.EQ.6) npts = 4  ! Quadrelateral           
          IF (npts.EQ.0) THEN
            IF (message) THEN
              WRITE(0,*)
              WRITE(0,*) '----------------------------------------'
              WRITE(0,*) '     IMPROPER INVERSION CELL FOUND'
              WRITE(0,*) '----------------------------------------'
              WRITE(0,*)
              message = .FALSE.
            ENDIF
            CYCLE
c            WRITE(0,*) '  985: PROBLEM OBJECT',iobj,n
c            STOP 'HALTING CODE'
          ENDIF
          IF (ascii) THEN
            IF (obj(iobj)%nside.NE.0) THEN
c....         Make sure that there is not more than 1 surface per side here:
              STOP 'SHOULDNT BE HERE AAA'
            ELSE
              WRITE(fp,'(2I7,1P,2E12.4,0P,2X,10(2F11.7),2I7)') 
     .          iobj,npts,obj(iobj)%quantity(1),obj(iobj)%path,
     .          (SNGL(obj(iobj)%v(1,i2)),
     .           SNGL(obj(iobj)%v(2,i2)),i2=1,npts),
     .          (0.0,0.0,i2=1,6),
     .          obj(iobj)%ik,obj(iobj)%ir
            ENDIF
          ELSE
            IF (obj(iobj)%nside.NE.0) THEN
c....         Make sure that there is not more than 1 surface per side here:
              STOP 'SHOULDNT BE HERE BBB'
            ELSE
              WRITE(fp) 
     .          iobj,npts,obj(iobj)%quantity(1),obj(iobj)%path,
     .          (SNGL(obj(iobj)%v(1,i2)),
     .           SNGL(obj(iobj)%v(2,i2)),i2=1,npts)
            ENDIF
          ENDIF
        ENDDO 
c...    Close stream:
        IF (.NOT.ascii) WRITE(fp) 666666666
        CLOSE(fp)
      ENDIF
c...  Clear arrays:
      DEALLOCATE(idum1)
      DEALLOCATE(ddum1)
      DEALLOCATE(rdum1)
      DEALLOCATE(vwlist)
      DEALLOCATE(gblist)
      DEALLOCATE(vwinter)
      DEALLOCATE(gbinter)
      DEALLOCATE(obinter)
      DEALLOCATE(mapchk)

      RETURN
 98   CALL ER('ProcessPixels','Unable to open inversion matrix '//
     .        'data file',*99)
 99   STOP
      END
c
c ====================================================================== 
c
      SUBROUTINE BuildPixels(MAXNPIXEL,npixel,pixel)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

c...  Input:
c      TYPE(type_options985) :: opt     
      INTEGER MAXNPIXEL,npixel
      TYPE(type_view) :: pixel(MAXNPIXEL)

c...  Need to add some checks that these are not violated (incase code is built 
c     without -C):

      INTEGER ix,iy,nxbin,nybin,i1,i2,count,ipixel,rpixel,
     .        ix1,ix2,iy1,iy2
      INTEGER, ALLOCATABLE :: imap(:,:)
      LOGICAL debug
      REAL    t,tlast
      REAL*8  xcen,ycen,zcen,xwidth,ywidth,xangle,yangle,factx,facty,
     .        mat(3,3),angle,dist,f,r,v,
     .        thetax,thetay,x,y,deltax,deltay,boundx,boundy

      REAL*8     PI
      PARAMETER (PI=3.14159265359)
     

c...  Identify region where views originate.  Some geometric object such as a point,
c     circle, square, etc. and determine how the views relate to this object and how many
c     ... :


c...    Point:

c....   Rectangle (can then be tilted vertically and horizontally, and rolled):
c...    Allow view distortions: 
      IF (.TRUE.) THEN

c...    Detector parameters:
        debug = .FALSE.

        IF (debug) THEN
c...      Debug:
          xcen = 2.5D0
          ycen = 0.0D0
          zcen = 0.0D0
          xwidth = 0.20D-07
          ywidth = 0.20D-07
          xangle = -13.90099009900990D0
          yangle = -25.66336633663366D0

c...      Number of rows and columns of pixels:
          nxbin =  1
          nybin =  1
        ELSEIF (.NOT..TRUE.) THEN
          xcen = 2.5D0
          ycen = 0.0D0
          zcen = 0.0D0
c          xwidth = 0.20D-07
c          ywidth = 0.20D-07
          xwidth = 1.80D0
          ywidth = 1.80D0
c          xwidth = 0.20D0
c          ywidth = 0.20D0
c          xangle = 22.5D0
c          yangle = 22.5D0
c          xangle = 65.0D0
c          yangle = 65.0D0
          xangle = 0.0D0
          yangle = 0.0D0
c          xangle = 54.0D0
c          yangle = 54.0D0

c...      Number of rows and columns of pixels:
          nxbin =  1
          nybin =  1
        ELSEIF (.NOT..TRUE.) THEN
          xcen = 2.15D0
          ycen = 1.42D0
          zcen = 0.0D0
          xwidth = 0.20D-07
          ywidth = 0.20D-07
c          xwidth = 0.20D0
c          ywidth = 0.20D0
c          xangle = 22.5D0
c          yangle = 22.5D0
c          xangle = 54.0D0
c          yangle = 54.0D0
          xangle = 65.0D0
          yangle = 65.0D0

c...      Number of rows and columns of pixels:
          nxbin =  11
          nybin =  11
        ELSE
          xcen = opt%cen(1)
          ycen = opt%cen(2)
          zcen = opt%cen(3)
          xwidth = opt%width(1)
          ywidth = opt%width(2)
          xangle = opt%angle(1)
          yangle = opt%angle(2)
c...      Number of rows and columns of pixels:
          nxbin = opt%nxbin
          nybin = opt%nybin
        ENDIF

        WRITE(0,'(A,3F10.3)') 'DATA:',xcen,ycen,zcen
        WRITE(0,'(A,3F10.3)') '    :',opt%roll,opt%tilt,
     .                                opt%swing
        WRITE(0,'(A,2F10.3)') '    :',xwidth,ywidth
        WRITE(0,'(A,2F10.3)') '    :',xangle,yangle
        WRITE(0,'(A,2I10  )') '    :',nxbin,nybin
        WRITE(0,'(A,2I10  )') '    :',opt%sa_nxbin,opt%sa_nybin
        WRITE(0,'(A,2F10.3)') '    :',opt%focallength,opt%distortion


        IF (opt%distortion.NE.0.0D0.OR.opt%focallength.NE.0.0D0) THEN

          dist = opt%distortion
          f    = opt%focallength

          thetax = 0.5D0 * xangle * PI / 180.0D0
          thetay = 0.5D0 * yangle * PI / 180.0D0

          x = 0.0D0
          deltax = 0.01D0
          count = 0
          DO WHILE(.TRUE.) 
            count = count + 1
c            v = DATAN( MAX(0.5D0, 0.0D0 * x**2 + 1.0D0) * x / f )
c            v = DATAN( MAX(0.5D0, dist * x**2 + 1.0D0) * x / f )
            v = DATAN( x / MAX(0.5D0, dist * x**2 + 1.0D0) / f )
            IF     (DABS(v-thetax).LT.1.0D-08.OR.
     .              count.EQ.1000) THEN
c              WRITE(0,*) 'WORKS X!',x,thetax*180.0D0/PI
              boundx = x
              EXIT
            ELSEIF ((v.GT.thetax.AND.deltax.GT.0.0D0).OR.
     .              (v.LT.thetax.AND.deltax.LT.0.0D0)) THEN
              deltax = -0.3D0 * deltax
            ENDIF
            x = x + deltax
c           WRITE(0,*) 'V:',v,thetax
c           WRITE(0,*) 'X:',x,deltax
c           WRITE(0,*) 
          ENDDO

          y = 0.0D0
          deltay = 0.01D0
          count = 0
          DO WHILE(.TRUE.) 
            count = count + 1
c            v = DATAN( MAX(0.5D0, 0.0D0 * y**2 + 1.0D0) * y / f )
            v = DATAN( y / MAX(0.5D0, dist * y**2 + 1.0D0) / f )
c            v = DATAN( MAX(0.5D0, dist * y**2 + 1.0D0) * y / f )
            IF     (DABS(v-thetay).LT.1.0D-08.OR.
     .              count.EQ.1000) THEN
c              WRITE(0,*) 'WORKS Y!',y
              boundy = y
              EXIT
            ELSEIF ((v.GT.thetay.AND.deltay.GT.0.0D0).OR.
     .              (v.LT.thetay.AND.deltay.LT.0.0D0)) THEN
              deltay = -0.3D0 * deltay
            ENDIF
            y = y + deltay
c           WRITE(0,*) 'V:',v,thetay
c           WRITE(0,*) 'Y:',y,deltay
c           WRITE(0,*) 
          ENDDO

        ENDIF


c...    Store index locations of first pixel for this detector:
        opt%det_istart(opt%ndet) = npixel + 1
c        IF (ndetector.GT.1) THEN 
c          xstart = pixel(npixel)%xindex 
c          ystart = pixel(npixel)%yindex 
c        ENDIF     



        DO iy = 1, nybin
          DO ix = 1, nxbin
            npixel = npixel + 1

            pixel(npixel)%valid = .TRUE.

            pixel(npixel)%weight = 1.0D0

            pixel(npixel)%integral = 0.0D0

            pixel(npixel)%index  = opt%ndet
            pixel(npixel)%tmp    = npixel    !*** TEMPORARY ***
            pixel(npixel)%xindex = ix
            pixel(npixel)%yindex = iy

            pixel(npixel)%nxbin = opt%sa_nxbin
            pixel(npixel)%nybin = opt%sa_nybin

            factx = -1.0D0 * (0.5D0 * (1.0D0 / DBLE(nxbin) - 1.0D0) +
     .                        DBLE(ix - 1) / DBLE(nxbin))

            facty = -1.0D0 * (0.5D0 * (1.0D0 / DBLE(nybin) - 1.0D0) +
     .                        DBLE(iy - 1) / DBLE(nybin))

            pixel(npixel)%v1(1) = factx * xwidth            ! Spatial horizontal location 
            pixel(npixel)%v1(2) = facty * ywidth 
            pixel(npixel)%v1(3) = 0.0D0

            pixel(npixel)%xwidth  = xwidth / DBLE(nxbin)      ! Spatial horizontal width
            pixel(npixel)%ywidth  = ywidth / DBLE(nybin)            

            IF (opt%distortion.EQ.0.0D0) THEN
              pixel(npixel)%xangle  = factx * xangle            ! View angle
              pixel(npixel)%yangle  = facty * yangle
            ELSE

              x = 2.0D0 * factx * boundx
              y = 2.0D0 * facty * boundy

              r = DSQRT(x**2 + y**2)

              x = x / MAX(0.5D0, dist * r**2 + 1.0D0)
              y = y / MAX(0.5D0, dist * r**2 + 1.0D0)

c              x = MAX(0.5D0, dist * r**2 + 1.0D0) * x
c              y = MAX(0.5D0, dist * r**2 + 1.0D0) * y

              pixel(npixel)%xangle = DATAN( x / f ) * 180.0D0 / PI
              pixel(npixel)%yangle = DATAN( y / f ) * 180.0D0 / PI

c             WRITE(0,*) 'XANGLE,YANGLE:',npixel,pixel(npixel)%xangle
c             WRITE(0,*) '             :',npixel,pixel(npixel)%yangle

             
c              WRITE(0,'(A,I6,8F10.4)') 
c     .          'DATA:',ix,factx,MAX(0.5D0, dist * r**2 + 1.0D0),
c     .          r,v,boundx,x,
c     .          pixel(npixel)%xangle,
c     .          factx * xangle

c              STOP 'sdfsd'

c              pixel(npixel)%xangle  = factx * xangle            ! View angle
c              pixel(npixel)%yangle  = facty * yangle              

            ENDIF


            SELECTCASE (opt%sa_opt)
              CASE (1) 
              pixel(npixel)%dxangle=xangle/DBLE(nxbin)*DBLE(opt%sa_par1)     ! View width 
              pixel(npixel)%dyangle=yangle/DBLE(nybin)*DBLE(opt%sa_par2)         
c              pixel(npixel)%dxangle=xangle/DBLE(nxbin)*DBLE(opt%sa_par1)     ! View width 
c              pixel(npixel)%dyangle=yangle/DBLE(nybin)*DBLE(opt%sa_par2)         

              CASE DEFAULT
                CALL ER('BuildPixels','Unknown solid angle '//
     .                  'option',*99)
            ENDSELECT


            IF (debug) THEN
              pixel(npixel)%xangle = xangle            ! View angle            
              pixel(npixel)%yangle = yangle            ! View angle
            ENDIF


c...        These are only needed for plotting, not actually used during LOS integration:
            pixel(npixel)%v2(1) = DSIN(pixel(npixel)%xangle *D_DEGRAD) *    ! Correct?
     .                            DCOS(pixel(npixel)%yangle *D_DEGRAD) * 
     .                            7.0D0 + pixel(npixel)%v1(1)                ! 7 meters/units?
            pixel(npixel)%v2(2) = DCOS(pixel(npixel)%xangle *D_DEGRAD) * 
     .                            DSIN(pixel(npixel)%yangle *D_DEGRAD) * 
     .                            7.0D0 + pixel(npixel)%v1(2)
            pixel(npixel)%v2(3) = DCOS(pixel(npixel)%xangle *D_DEGRAD) * 
     .                            DCOS(pixel(npixel)%yangle *D_DEGRAD) * 
     .                            7.0D0 + pixel(npixel)%v1(3)
c            pixel(npixel)%v2(3) = -pixel(npixel)%v2(3)



c           *** better to store these in a detector array, rather than taking up space with each pixel, or are the savings nominal...
c           *** use opt%det_ array...
            pixel(npixel)%rot(1) = DBLE(opt%tilt *PI/180.0D0) ! x-axis (tilt) 
            pixel(npixel)%rot(2) = DBLE(opt%swing*PI/180.0D0) ! y-axis (swing)
            pixel(npixel)%rot(3) = DBLE(opt%roll *PI/180.0D0) ! z-axis (roll) ! Set detector orientation

            pixel(npixel)%trans(1) = xcen
            pixel(npixel)%trans(2) = ycen
            pixel(npixel)%trans(3) = zcen



            IF (.NOT..TRUE.) THEN
c...          Rotate about z-axis (roll):
              CALL Calc_Transform2(mat,0.0D0,1,0)
              angle = pixel(npixel)%rot(3)
              CALL Calc_Transform2(mat,angle,3,1)
              CALL Transform_Vect(mat,pixel(npixel)%v1)
              CALL Transform_Vect(mat,pixel(npixel)%v2)
c...          Rotate about x-axis (tilt):                       !... better to do swing before tilt? 
              CALL Calc_Transform2(mat,0.0D0,1,0)
              angle = pixel(npixel)%rot(1)
              CALL Calc_Transform2(mat,angle,1,1)
              CALL Transform_Vect(mat,pixel(npixel)%v1)
              CALL Transform_Vect(mat,pixel(npixel)%v2)
c...          Rotate about y-axis (swing):
              CALL Calc_Transform2(mat,0.0D0,1,0)
              angle = pixel(npixel)%rot(2)
              CALL Calc_Transform2(mat,angle,2,1)
              CALL Transform_Vect(mat,pixel(npixel)%v1)
              CALL Transform_Vect(mat,pixel(npixel)%v2)
c...          Translate:
              pixel(npixel)%v1(1:3) = pixel(npixel)%v1(1:3) + 
     .                                pixel(npixel)%trans(1:3)
              pixel(npixel)%v2(1:3) = pixel(npixel)%v2(1:3) + 
     .                                pixel(npixel)%trans(1:3)
            ENDIF
 

c            WRITE(0,*) 'X!;',npixel,pixel(npixel)%x2
          ENDDO
        ENDDO

      ENDIF


c...  Mask pixels:
      DO i1 = 1, opt%mask_num

        WRITE(0,*) 'MASK?', opt%mask_num,opt%mask_opt(i1)

        SELECTCASE (opt%mask_opt(i1))
          CASE (-1)  ! Image loaded but no masking...
          CASE (1)
c...        Scan loaded image file and blank any pixels corresponding to 
c           artifical black (image value is 0):
            IF (opt%nxbin.NE.opt%img_nxbin.OR.
     .          opt%nybin.NE.opt%img_nybin) 
     .        CALL ER('BuildPixels','Image and CDD resolutions do '//
     .                'not match',*99)

            DO ipixel = opt%det_istart(opt%ndet), npixel
              ix = pixel(ipixel)%xindex
              iy = pixel(ipixel)%yindex
              IF (opt%img_image(ix,iy).EQ.0.0D0) THEN
                pixel(ipixel)%valid=.FALSE.
c                WRITE(0,*) 'got one:',ipixel
              ENDIF
            ENDDO

          CASE (2) 
c...        Invert mask:
            DO ipixel = opt%det_istart(opt%ndet), npixel
              pixel(ipixel)%valid = .NOT.pixel(ipixel)%valid
            ENDDO
            
c ...        Circular:        
c            DO ipixel = opt%det_istart(opt%ndet), npixel
c              ix = (pixel(ipixel)%xindex - opt%nxbin / 2)
c              iy = (pixel(ipixel)%yindex - opt%nxbin / 2)
c              rpixel = NINT(SQRT(REAL(ix**2) + REAL(iy**2)))
c              IF (rpixel.GT.opt%nxbin/2) pixel(ipixel)%valid = .FALSE.
c            ENDDO

          CASE (3)
c...      
            DO ipixel = opt%det_istart(opt%ndet), npixel
              ix = pixel(ipixel)%xindex
              iy = pixel(ipixel)%yindex
c             IF (iy.LT.70) THEN
c             IF (iy.GT.40.AND.iy.LT.70) THEN
c             IF (iy.GT.85.AND.iy.LT.115) THEN
c             IF (iy.LT.115) THEN
c             IF (iy.LT.175) THEN
             IF (ix.LT.75.AND.iy.GT.75) THEN
                pixel(ipixel)%valid = .FALSE.
              ENDIF
             IF (ix.LT.40.AND.iy.GT.40) THEN
                pixel(ipixel)%valid = .FALSE.
              ENDIF
            ENDDO

          CASE (4)  ! Along a line...

            ALLOCATE(imap(opt%nxbin,opt%nybin))

            DO ipixel = opt%det_istart(opt%ndet), npixel
              imap(pixel(ipixel)%xindex,pixel(ipixel)%yindex) = ipixel
            ENDDO

            ix1 = MAX(1,opt%mask_vtx(i1,1) / opt%img_nxratio)
            iy1 = MAX(1,opt%mask_vtx(i1,2) / opt%img_nyratio)
            ix2 = MAX(1,opt%mask_vtx(i1,3) / opt%img_nxratio)
            iy2 = MAX(1,opt%mask_vtx(i1,4) / opt%img_nyratio)

            tlast = -1.0
            DO t = 0.0, 1.0, 0.0005
              IF (t.EQ.tlast) CYCLE
              tlast = t

              ix = MAX(1,NINT((1.0 - t) * REAL(ix1) + t * REAL(ix2)))
              iy = MAX(1,NINT((1.0 - t) * REAL(iy1) + t * REAL(iy2)))

c              WRITE(0,*) 'I:',t,ix,iy,imap(ix,iy)
              pixel(imap(ix,iy))%valid = .FALSE.
c              if (t.EQ.0.1) STOP
            ENDDO

            WRITE(0,*) 'I:',ix1,iy1,ix2,iy2


            DEALLOCATE(imap)
            WRITE(0,*) 'RATIOS:',opt%img_nxratio,opt%img_nyratio
c            STOP 'HERE'

          CASE DEFAULT
            WRITE(0,*) 'MASK TYPE:',i1,opt%mask_opt(i1)      
            CALL ER('BuildPixels','Unknown mask type',*99)
        ENDSELECT

      ENDDO


c.... Add index offset(s) for multiple detectors:
c      IF (ndetector.GT.1) THEN
c        DO ipixel = istart, npixel
c          pixel(ipixel)%yindex = pixel(ipixel)%yindex + ystart
c        ENDDO
c      ENDIF
      opt%det_iend(opt%ndet) = npixel


c...  Fill out pixel view as necessary, including focus options (need to think) -- add subroutine for creating the views:
c...  Rotate and translate the views :  
c...  Do ray trace and add reflection views to list:
c...  Do all this for one pixel at a time which is farmed out:
c...  As each pixel is returned, increment pixel value, which isn't necessarily the calling pixel if the 
c     focus is off, then clear the memory for the views associated with that pixel (or plot first, or store in file): 


c...  Allocate space for line shape storage for unmasked pixels:



      RETURN
 99   STOP
      END
c
c ====================================================================== 
c
      SUBROUTINE BuildObjects
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none


c      INCLUDE 'params'
c      INCLUDE 'cgeom'
c      INCLUDE 'comtor'
c      INCLUDE 'pindata'
c      INCLUDE 'slcom'


      



      INTEGER i1,i2,i3,ix,iy,in1,in2,ir,id,nsector,isector,
     .        ielement
      REAL    angle,dangle,ang,dangle2
      REAL*8  x1,z1,x2,z2,p1(3,8),p2(3,8)



c...  Need to have the target plates covered at some point... perhaps the 
c     grid polygons need to be included in the surface intersection check
c     with the surfaces along the targets flagged as solid surfaces...


      nobj = 0


      IF (.NOT..TRUE.) THEN
c...    Test surface:
        nobj = 1

        obj(nobj)%index       = nobj
        obj(nobj)%type        = OP_EMPTY
        obj(nobj)%mode        = 0      
        obj(nobj)%surface     = 1      ! SOLID
        obj(nobj)%wedge1      = 0
        obj(nobj)%wedge2      = 0
        obj(nobj)%colour      = 1
        obj(nobj)%orientation = 1      ! CW
        obj(nobj)%ik          = 0
        obj(nobj)%ir          = 0
        obj(nobj)%in          = i1
        obj(nobj)%ivolume     = 0
        obj(nobj)%nsur        = 1
        obj(nobj)%gsur(1:1)   = GT_TD
        obj(nobj)%nver        = 4

        obj(nobj)%nmap        = 0

        obj(nobj)%tsur(1) = SP_VESSEL_WALL
        obj(nobj)%npts(1) = 4
        obj(nobj)%ipts(1,1) = 1
        obj(nobj)%ipts(2,1) = 2
        obj(nobj)%ipts(3,1) = 3
        obj(nobj)%ipts(4,1) = 4

c        obj(nobj)%v(1,1) = 0.0 !-1.0
c        obj(nobj)%v(2,1) = 0.0
c        obj(nobj)%v(3,1) = 0.0
c        obj(nobj)%v(1,2) = 0.0 !-1.0
c        obj(nobj)%v(2,2) = 1.0
c        obj(nobj)%v(3,2) = 0.0
c        obj(nobj)%v(1,3) = 0.0
c        obj(nobj)%v(2,3) = 1.0
c        obj(nobj)%v(3,3) = 1.0
c        obj(nobj)%v(1,4) = 0.0
c        obj(nobj)%v(2,4) = 0.0
c        obj(nobj)%v(3,4) = 1.0

        obj(nobj)%v(1,4) = -1.0 
        obj(nobj)%v(2,4) = 0.0
        obj(nobj)%v(3,4) = 0.0
        obj(nobj)%v(1,3) = -1.0 
        obj(nobj)%v(2,3) = 1.0
        obj(nobj)%v(3,3) = 0.0
        obj(nobj)%v(1,2) = 0.0 !
        obj(nobj)%v(2,2) = 1.0
        obj(nobj)%v(3,2) = 1.0
        obj(nobj)%v(1,1) = 0.0 !
        obj(nobj)%v(2,1) = 0.0
        obj(nobj)%v(3,1) = 1.0

        RETURN
      ENDIF

c
c     *******************************************
c     VOLUME INTEGRATION OBJECTS MUST COME FIRST!
c     *******************************************
c

c...    Integration mesh:
c        IF (opt%ob_invgrd.GT.0) 
c     .    CALL BuildInversionMesh(0)


      IF (opt%obj_num.NE.0) THEN
c        STOP 'NOT READY YET'
        DO ielement = 1, opt%obj_num
          SELECTCASE (opt%obj_type(ielement))
            CASE (1)
              CALL ProcessMagneticGrid    (ielement)
            CASE (2)
              CALL ProcessTriangleGrid    (ielement)
            CASE (3)
              CALL BuildInversionMesh     (ielement)
            CASE (4)
              CALL LoadVesselStructures   (ielement)
            CASE (5)
              CALL LoadLineSegmentFile    (ielement)
            CASE (6)
              CALL ProcessTetrahedronGrid (ielement)
            CASE (7)
              CALL LoadImageReconstruction(ielement)
            CASE DEFAULT
              CALL User_CustomObjects     (ielement)
          ENDSELECT
        ENDDO



c        DO i1 = 1, nobj
c          IF (obj(i1)%ik.EQ.ikto.AND.obj(i1)%ir.EQ.irsep) THEN
cc             WRITE(0,*) 'r:',obj(i1)%v(1,1:obj(nobj)%nver)
c          ENDIF
c        ENDDO


      ENDIF

      IF (.TRUE.) THEN
c...    OLD METHOD OF SPECIFYING THE 3D OBJECTS



c        WRITE(0,*) 'MAG GRID OPT',opt%ob_stdgrd
cc...    Standard magnetic grid:
c        IF (opt%ob_stdgrd.GT.0) 
c     .    CALL ProcessMagneticGrid(0)

c...    Triangle grid:
c        IF (opt%ob_trigrd.GT.0) 
c     .    CALL ProcessTriangleGrid(0)

c...    Vessel structures from .raw CAD generated files:
c        IF (opt%ob_raw_num.GT.0) 
c     .    CALL LoadVesselStructures(0)

c...    Trace magnetic field lines:
c        IF (opt%ob_line.GT.0) 
c     .    CALL TraceMagneticFieldLines

      ENDIF
      

c      IF (opt%ob_wall.GT.0) THEN
cc...    Standard PIN wall (xVESM):
c
c        IF (opt%ob_wall.EQ.2) THEN
c
c          DO i1 = 1, nvesm+nvesp
c        
cc...        Don't draw anything above the x-point:
cc            IF (zvesm(i1,1).GT.zxp.OR.zvesm(i1,2).GT.zxp) CYCLE
cc            IF (zvesm(i1,1).GT.zxp.AND.rvesm(i1,1).GT.rxp) CYCLE
c       
cc...        Don't draw other things:
c            IF (jvesm(i1).NE.7.AND.jvesm(i1).NE.8.AND.
c     .          jvesm(i1).NE.0) CYCLE
c        
c            IF (nobj+1.GT.MAX3D) 
c     .        CALL ER('BuildObjects','Insufficient array bounds '//
c     .                'for all objects',*98)     
c
c            nobj = nobj + 1
c
c            obj(nobj)%index       = 1 ! nobj
c            obj(nobj)%type        = OP_EMPTY
c            obj(nobj)%mode        = 0      
c            obj(nobj)%surface     = 1      ! SOLID
c            obj(nobj)%wedge1      = 0
c            obj(nobj)%wedge2      = 0
c            obj(nobj)%colour      = 1
c            obj(nobj)%orientation = 1      ! CW
c            obj(nobj)%ik          = i1
c            obj(nobj)%ir          = 0
c            obj(nobj)%in          = i1
c            obj(nobj)%ivolume     = 0
c            obj(nobj)%nsur        = 1
c            obj(nobj)%gsur(1:1)   = GT_TC
c            obj(nobj)%nver        = 2
c            obj(nobj)%tsur(1)     = SP_VESSEL_WALL
c            obj(nobj)%reflec(1)   = opt%ob_wall_reflec
c            obj(nobj)%npts(1)     = 2
c            obj(nobj)%ipts(1,1)   = 1
c            obj(nobj)%ipts(2,1)   = 2
c            obj(nobj)%nmap(1) = 0
c
cc...        Vertices:
c            obj(nobj)%v(1,2) = DBLE(rvesm(i1,1))
c            obj(nobj)%v(2,2) = DBLE(zvesm(i1,1))
c            obj(nobj)%v(3,2) = 0.0D0
c            obj(nobj)%v(1,1) = DBLE(rvesm(i1,2))
c            obj(nobj)%v(2,1) = DBLE(zvesm(i1,2))
c            obj(nobj)%v(3,1) = 0.0D0
c
c          ENDDO
c
c        ELSEIF (opt%ob_wall.EQ.1.OR.opt%ob_wall.EQ.3) THEN
c
c          nsector = opt%ob_nsector
c          IF (nsector.EQ.-1) nsector = eirntorseg
c          dangle = 360.0 / REAL(nsector) / RADDEG
c          nsector = NINT(REAL(nsector) * 
c     .                   (opt%ob_angle_end - opt%ob_angle_start)/360.0)
c
c          dangle2 = dangle
c
c          DO i1 = 1, nvesm+nvesp
c        
cc...        Don't draw anything above the x-point:
cc            IF (zvesm(i1,1).GT.zxp.OR.zvesm(i1,2).GT.zxp) CYCLE
cc            IF (zvesm(i1,1).GT.zxp.AND.rvesm(i1,1).GT.rxp) CYCLE
c       
cc...        Don't draw other things:
c            IF (jvesm(i1).NE.7.AND.jvesm(i1).NE.8.AND.
c     .          jvesm(i1).NE.0) CYCLE
c        
cc            DO ang = 0.0, 2.0*PI-dangle, dangle   !DO ang = 0.0, 0.45*2.0*PI-dangle, dangle
c
c            DO isector = 1, nsector
c           
c              dangle = dangle2    ! For opt%ob_wall.EQ.3, temp...
c
c              ang = REAL(isector - 1) * dangle + 
c     .              opt%ob_angle_start / RADDEG
c
c              IF (opt%ob_wall.EQ.3) THEN
c                IF (MOD(isector,2).EQ.1) THEN
c                  dangle = dangle2 * 1.60
c                ELSE
c                  dangle = dangle2 * 0.40
c                ENDIF
c              ENDIF        
c
c              IF (nobj+1.GT.MAX3D) 
c     .          CALL ER('BuildObjects','Insufficient array bounds '//
c     .                  'for all objects',*98)     
c
c              nobj = nobj + 1
c        
c              obj(nobj)%index       = nobj
c              obj(nobj)%type        = OP_EMPTY
c              obj(nobj)%mode        = 0      
c              obj(nobj)%surface     = 1      ! SOLID
c              obj(nobj)%wedge1      = 0
c              obj(nobj)%wedge2      = 0
c              obj(nobj)%colour      = 1
c              obj(nobj)%orientation = 1      ! CW
c              obj(nobj)%ik          = 0
c              obj(nobj)%ir          = 0
c              obj(nobj)%in          = i1
c              obj(nobj)%ivolume     = 0
c              obj(nobj)%nsur        = 1
c              obj(nobj)%gsur(1:1)   = GT_TD
c              obj(nobj)%nver        = 4
c        
c              obj(nobj)%tsur(1) = SP_VESSEL_WALL
c              obj(nobj)%reflec(1) = opt%ob_wall_reflec
c              obj(nobj)%npts(1) = 4
c              obj(nobj)%nmap(1) = 0
c        
cc...          A little over done here, with all the p's (from old code), clean:
c              p1(1,1) = DBLE(rvesm(i1,1)) * DCOS(DBLE(-0.5 * dangle))
c              p1(2,1) = DBLE(zvesm(i1,1))
c              p1(3,1) = DBLE(rvesm(i1,1)) * DSIN(DBLE(-0.5 * dangle))
c              p2(1,1) = DBLE(rvesm(i1,2)) * DCOS(DBLE(-0.5 * dangle))
c              p2(2,1) = DBLE(zvesm(i1,2))
c              p2(3,1) = DBLE(rvesm(i1,2)) * DSIN(DBLE(-0.5 * dangle))
cc              p1(1,1) = DBLE(rvesm(i1,1))
cc              p1(2,1) = DBLE(zvesm(i1,1))
cc              p1(3,1) = DBLE(rvesm(i1,1)) * DTAN(DBLE(-0.5 * dangle))
cc              p2(1,1) = DBLE(rvesm(i1,2))
cc              p2(2,1) = DBLE(zvesm(i1,2))
cc              p2(3,1) = DBLE(rvesm(i1,2)) * DTAN(DBLE(-0.5 * dangle))
c	
c              p1(1,2) = p2(1,1)
c              p1(2,2) = p2(2,1)
c              p1(3,2) = p2(3,1)
c              p2(1,2) = DBLE(rvesm(i1,2)) * DCOS(DBLE(+0.5 * dangle))
c              p2(2,2) = DBLE(zvesm(i1,2))
c              p2(3,2) = DBLE(rvesm(i1,2)) * DSIN(DBLE(+0.5 * dangle))
cc              p2(1,2) = DBLE(rvesm(i1,2))
cc              p2(2,2) = DBLE(zvesm(i1,2))
cc              p2(3,2) = DBLE(rvesm(i1,2)) * DTAN(DBLE(+0.5 * dangle))
c
c        
c              p1(1,3) = p2(1,2)
c              p1(2,3) = p2(2,2)
c              p1(3,3) = p2(3,2)
c              p2(1,3) = DBLE(rvesm(i1,1)) * DCOS(DBLE(+0.5 * dangle))
c              p2(2,3) = DBLE(zvesm(i1,1))
c              p2(3,3) = DBLE(rvesm(i1,1)) * DSIN(DBLE(+0.5 * dangle))
cc              p2(1,3) = DBLE(rvesm(i1,1))
cc              p2(2,3) = DBLE(zvesm(i1,1))
cc              p2(3,3) = DBLE(rvesm(i1,1)) * DTAN(DBLE(+0.5 * dangle))
c	
c              p1(1,4) = p2(1,3)
c              p1(2,4) = p2(2,3)
c              p1(3,4) = p2(3,3)
c              p2(1,4) = p1(1,1)
c              p2(2,4) = p1(2,1)
c              p2(3,4) = p1(3,1)
c          
c              DO i2 = 1, 4
c                x1 = p1(1,i2)
c                z1 = p1(3,i2)
c                x2 = p2(1,i2)
c                z2 = p2(3,i2)
c                p1(1,i2) =  DCOS(DBLE(ang)) * x1 - DSIN(DBLE(ang)) *z1
c                p1(3,i2) = +DSIN(DBLE(ang)) * x1 + DCOS(DBLE(ang)) *z1
c                p2(1,i2) =  DCOS(DBLE(ang)) * x2 - DSIN(DBLE(ang)) *z2
c                p2(3,i2) = +DSIN(DBLE(ang)) * x2 + DCOS(DBLE(ang)) *z2
c         
c                obj(nobj)%ipts(i2,1) = i2
c                obj(nobj)%v(1:3,i2) = p1(1:3,i2)
c              ENDDO
c          
cc... Reverse:
c              p1(1:3,1) = obj(nobj)%v(1:3,2)
c              obj(nobj)%v(1:3,2) = obj(nobj)%v(1:3,4)  
c              obj(nobj)%v(1:3,4) = p1(1:3,1)
c  
c            ENDDO
c          ENDDO
c        ELSE
c          CALL ER('BuildObjects','Bad wall option',*99)
c        ENDIF
c      ENDIF
c
c      IF (opt%ob_targ.GT.0) THEN
cc...    Targets:
c        
c        IF (opt%ob_targ.EQ.2) THEN
c
c          DO ir = irsep, nrs
c            IF (idring(ir).EQ.BOUNDARY) CYCLE
c
c            DO i1 = 1, 2
c
c              IF (nobj+1.GT.MAX3D) 
c     .          CALL ER('BuildObjects','Insufficient array bounds '//
c     .                  'for all objects',*98)     
c
c              nobj = nobj + 1
c
c              obj(nobj)%index       = nobj
c              obj(nobj)%type        = OP_EMPTY
c              obj(nobj)%mode        = 0      
c              obj(nobj)%surface     = 1      ! SOLID
c              obj(nobj)%wedge1      = 0
c              obj(nobj)%wedge2      = 0
c              obj(nobj)%colour      = 4
c              obj(nobj)%orientation = 1      ! CW
c              obj(nobj)%ik          = 0
c              obj(nobj)%ir          = ir
c              obj(nobj)%in          = 0
c              obj(nobj)%ivolume     = 0
c              obj(nobj)%nsur        = 1
c              obj(nobj)%gsur(1)     = GT_TC
c              obj(nobj)%nver        = 2
c              obj(nobj)%tsur(1)     = SP_VESSEL_WALL
c              obj(nobj)%reflec(1)   = opt%ob_targ_reflec
c              obj(nobj)%npts(1)     = 2
c              obj(nobj)%ipts(1,1)   = 1
c              obj(nobj)%ipts(2,1)   = 2
c              obj(nobj)%nmap(1)     = 0 
c...          Vertices:
c              IF (i1.EQ.1) THEN
c...            Inner target:
c                id = korpg(1,ir)
c                in1 = 1
c                in2 = 2
c              ELSE
cc...            Outer: 
c                id = korpg(nks(ir),ir)
c                in1 = 3
c                in2 = 4
c              ENDIF
c              obj(nobj)%v(1,2) = DBLE(rvertp(in1,id))
c              obj(nobj)%v(2,2) = DBLE(zvertp(in1,id))
c              obj(nobj)%v(3,2) = 0.0D0
c              obj(nobj)%v(1,1) = DBLE(rvertp(in2,id))
c              obj(nobj)%v(2,1) = DBLE(zvertp(in2,id))
c              obj(nobj)%v(3,1) = 0.0D0
c
c            ENDDO
c
c          ENDDO
c       
c        ELSEIF (opt%ob_targ.EQ.1) THEN
c          nsector = opt%ob_nsector
c          IF (nsector.EQ.-1) nsector = eirntorseg
c
c          dangle = 360.0 / REAL(nsector) / RADDEG
c        
c          DO ir = irsep, nrs
c            IF (idring(ir).EQ.BOUNDARY) CYCLE
c        
c            DO ang = 0.0, 2.0*PI-dangle, dangle   !DO ang = 0.0, 0.45*2.0*PI-dangle, dangle
c        
c              IF (nobj+1.GT.MAX3D) 
c     .          CALL ER('BuildObjects','Insufficient array bounds '//
c     .                  'for all objects',*98)     
c
c              DO i1 = 1, 2
c
c                nobj = nobj + 1
c         
c                obj(nobj)%index       = nobj
c                obj(nobj)%type        = OP_EMPTY
c                obj(nobj)%mode        = 0      
c                obj(nobj)%surface     = 1      ! SOLID
c                obj(nobj)%wedge1      = 0
c                obj(nobj)%wedge2      = 0
c                obj(nobj)%colour      = 4
c                obj(nobj)%orientation = 1      ! CW
c                obj(nobj)%ik          = 0
c                obj(nobj)%ir          = ir
c                obj(nobj)%in          = 0
c                obj(nobj)%ivolume     = 0
c                obj(nobj)%nsur        = 1
c                obj(nobj)%gsur(1:1)   = GT_TD
c                obj(nobj)%nver        = 4
c                obj(nobj)%tsur(1) = SP_VESSEL_WALL
c                obj(nobj)%reflec(1) = opt%ob_targ_reflec
c                obj(nobj)%npts(1) = 4
c                obj(nobj)%nmap(1) = 0
c        
c                IF (i1.EQ.1) THEN
cc...              Inner target:
c                  id = korpg(1,ir)
c                  in1 = 1
c                  in2 = 2
c                ELSE
cc...              Outer: 
c                  id = korpg(nks(ir),ir)
c                  in1 = 3
c                  in2 = 4
c                ENDIF
c
cc...            A little over done here, with all the p's (from old code), clean:
c                p1(1,1) = DBLE(rvertp(in1,id))
c                p1(2,1) = DBLE(zvertp(in1,id))
c                p1(3,1) = DBLE(rvertp(in1,id)) * DTAN(DBLE(-0.5*dangle))
c                p2(1,1) = DBLE(rvertp(in2,id))
c                p2(2,1) = DBLE(zvertp(in2,id))
c                p2(3,1) = DBLE(rvertp(in2,id)) * DTAN(DBLE(-0.5*dangle))
c	
c                p1(1,2) = p2(1,1)
c                p1(2,2) = p2(2,1)
c                p1(3,2) = p2(3,1)
c                p2(1,2) = DBLE(rvertp(in2,id))
c                p2(2,2) = DBLE(zvertp(in2,id))
c                p2(3,2) = DBLE(rvertp(in2,id)) * DTAN(DBLE(+0.5*dangle))
c        
c                p1(1,3) = p2(1,2)
c                p1(2,3) = p2(2,2)
c                p1(3,3) = p2(3,2)
c                p2(1,3) = DBLE(rvertp(in1,id))
c                p2(2,3) = DBLE(zvertp(in1,id))
c                p2(3,3) = DBLE(rvertp(in1,id)) * DTAN(DBLE(+0.5*dangle))
c	
c                p1(1,4) = p2(1,3)
c                p1(2,4) = p2(2,3)
c                p1(3,4) = p2(3,3)
c                p2(1,4) = p1(1,1)
c                p2(2,4) = p1(2,1)
c                p2(3,4) = p1(3,1)
c
c                DO i2 = 1, 4
c                  x1 = p1(1,i2)
c                  z1 = p1(3,i2)
c                  x2 = p2(1,i2)
c                  z2 = p2(3,i2)
c                  p1(1,i2) =  DCOS(DBLE(ang)) * x1 - DSIN(DBLE(ang))*z1
c                  p1(3,i2) = +DSIN(DBLE(ang)) * x1 + DCOS(DBLE(ang))*z1
c                  p2(1,i2) =  DCOS(DBLE(ang)) * x2 - DSIN(DBLE(ang))*z2
c                  p2(3,i2) = +DSIN(DBLE(ang)) * x2 + DCOS(DBLE(ang))*z2
c         
c                  obj(nobj)%ipts(i2,1) = i2
c                  obj(nobj)%v(1:3,i2) = p1(1:3,i2)
c                ENDDO
cc...            Reverse:
c                p1(1:3,1) = obj(nobj)%v(1:3,2)
c                obj(nobj)%v(1:3,2) = obj(nobj)%v(1:3,4)  
c                obj(nobj)%v(1:3,4) = p1(1:3,1)
c              ENDDO  
c
c            ENDDO
c          ENDDO
c        ELSE
c          CALL ER('BuildObjects','Bad target option',*99)
c        ENDIF
c      ENDIF

 98   RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: DumpImage
c
      SUBROUTINE DumpImage(npixel,pixel)
c      SUBROUTINE DumpImage(opt,npixel,pixel)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

c      TYPE(type_options985) :: opt     
      INTEGER npixel
      TYPE(type_view)      :: pixel(npixel)

      INTEGER   fp,ix,ix1,iy,i1,ierr,nxbin,nybin,idet
      CHARACTER file*1024
      REAL*8, ALLOCATABLE :: image(:,:)


      DO idet = 1, opt%ndet

        nxbin = opt%det_nxbin(idet)
        nybin = opt%det_nybin(idet)

        ALLOCATE(image(nxbin,nybin))

        image = 0.0D0
        DO i1 = opt%det_istart(idet), opt%det_iend(idet)
          image(pixel(i1)%xindex,
     .          pixel(i1)%yindex) = pixel(i1)%integral(1)
        ENDDO

        IF (.TRUE.) THEN
          fp = 99
          WRITE(file,'(A,I1,A)')
          file = 'output.'//TRIM(opt%det_fname(idet))//'.ray.img'
          WRITE(0,*) 'DUMP IMAGE:',file(1:LEN_TRIM(file))
          OPEN(fp,FILE=file(1:LEN_TRIM(file)),
     .         FORM='FORMATTED',STATUS='REPLACE',ERR=97)
          WRITE(fp,'(F6.2,2I7)') 1.0,nxbin,nybin
          DO iy = 1, nybin
            DO ix = 1, nxbin, 5
              WRITE(fp,'(1P,10E20.12)')                  ! Better way to write this!
     .          (image(ix1,iy),ix1=ix,MIN(ix+4,nxbin))
            ENDDO
          ENDDO     
          CLOSE(fp)
          WRITE(0,*) 'DONE'

c          IF (.FALSE.) THEN
c            file =opt%det_fname(idet)(1:LEN_TRIM(opt%det_fname(idet)))//
c     .            '.b'
c            WRITE(0,*) 'DUMP IMAGE:',file(1:LEN_TRIM(file))
c            OPEN(fp,FILE=file(1:LEN_TRIM(file)),
c     .           FORM='FORMATTED',STATUS='REPLACE',ERR=97)
c            WRITE(fp,'(I7)') opt%det_iend(idet)
c            DO i1 = opt%det_istart(idet), opt%det_iend(idet)
c              WRITE(fp,'(I7,1P,E22.12)') i1,pixel(i1)%integral(1)
c            ENDDO     
c            CLOSE(fp)
c          ENDIF

        ELSE
          CALL ER('DumpImage','Unknown image format specified',*99)
        ENDIF

        DEALLOCATE(image)

      ENDDO

      RETURN
 97   CALL ER('Main989','Unable to create image file',*99)
 99   WRITE(0,*) 'FILE:',file(1:LEN_TRIM(file))
      STOP
      END
c
c ======================================================================
c
      SUBROUTINE LoadEmissionDatabase(dummy,ni,nr)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

 
      CHARACTER dummy*(*)
      INTEGER   ni,nr

      INTEGER   i1,idum1
      REAL      rdum1
      CHARACTER cdum1*256

      SELECTCASE (opt%int_database(opt%int_num))
        CASE (0)
          READ(dummy,*) cdum1,(idum1,i1=1,ni),
     .                        (rdum1,i1=1,nr),idum1,
     .                  opt%int_line(opt%int_num)
        CASE (1)
          READ(dummy,*) cdum1,(idum1,i1=1,ni),
     .                        (rdum1,i1=1,nr),idum1,
     .                  opt%int_adasid(opt%int_num),
     .                  opt%int_adasyr(opt%int_num),
     .                  opt%int_adasex(opt%int_num),
     .                  opt%int_isele (opt%int_num),
     .                  opt%int_iselr (opt%int_num),
     .                  opt%int_iselx (opt%int_num),
     .                  opt%int_iseld (opt%int_num)
        CASE DEFAULT
          STOP 'SORRY, NO USER OPTIONS YET FOR DATABASE'
      ENDSELECT

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: Main985
c
      SUBROUTINE Main985(iopt)
c      SUBROUTINE Main985(iopt,opt985,
c     .                   MAXPIXEL,npixel,pixel,
c     .                   MAX3D985,nobj985,obj985,
c     .                   image)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

      INTEGER                          iopt
c      TYPE(type_options985), TARGET :: opt985
c      INTEGER                          MAXPIXEL,npixel
c      TYPE(type_view)               :: pixel(MAXPIXEL)
c      INTEGER                          MAX3D985,nobj985
c      TYPE(type_3D_object), TARGET  :: obj985(MAX3D985)

      REAL*8, PARAMETER :: PI = 3.141592654

      INTEGER   fp,ndetector,status
      INTEGER   save_opt_ccd  ! Hack
      REAL*8    image(1100,1100)  ! Not sure this really needs to be passed around...
      CHARACTER dummy*1024

c... *TEMP* interpolation test variables, see bottom of routine:
      INTEGER i1,i2,ivtx,img_nxbin,img_nybin
      REAL x(100),y(100),f(100),x1(10),y1(10),f1(10)  
      REAL*8 mat(3,3) 

      INTEGER MAXNPIXEL,npixel
      TYPE(type_view), ALLOCATABLE :: pixel(:)

c      opt => opt985
c      obj => obj985
c      MAX3D = MAX3D985

      WRITE(0,*) 'HERE IN 985'

      WRITE(0,*) '  ALLOCATING OBJECTS'
      MAX3D = 400000 ! 200000
      ALLOCATE(obj(MAX3D))

      CALL ALLOC_SURFACE(-1,MP_INITIALIZE)

      WRITE(0,*) '  ALLOCATING PIXELS'
      MAXNPIXEL=480*640 ! 1000*1000
      ALLOCATE(pixel(MAXNPIXEL))

      CALL ALLOC_CHORD(10000)  ! Just for viewing! (make smaller!)

      opt%load = 1

      grd_ntorseg = 48

      WRITE(0,*) 'MAXNPIXEL:',maxnpixel

      fp = 5

      npixel = 0

      opt%img_nxratio = 1
      opt%img_nyratio = 1


      IF (opt%load.EQ.1) THEN
        opt%obj_num = 0
        opt%int_num = 0
        opt%ref_num = 0
        opt%ndet = 0
        CALL LoadOptions985_New(opt,ALL_OPTIONS,status)
      ENDIF

      IF (opt%load.EQ.1) THEN
        status = 0
        save_opt_ccd = 0
        DO WHILE(status.EQ.0)
          opt%ccd = 0  ! Hack
          opt%mask_num = 0

          CALL LoadOptions985_New(opt,DETECTOR_ONLY,status)
c          CALL LoadDetectorOptions985(opt)
          IF (opt%ccd.GT.0) THEN

            IF (save_opt_ccd.EQ.0) save_opt_ccd = opt%ccd

            opt%ndet = opt%ndet + 1
            opt%det_nxbin(opt%ndet) = opt%nxbin
            opt%det_nybin(opt%ndet) = opt%nybin
            opt%det_fname(opt%ndet) = opt%fmap

            IF(opt%img_opt.NE.0) THEN
              WRITE(0,*) '  BINNING IMAGE!'
              img_nxbin = opt%img_nxbin
              img_nybin = opt%img_nybin
              CALL BinImage(opt%img_nxbin,opt%img_nybin,opt%img_image,
     .                      opt%nxbin,opt%nybin)
              opt%img_nxratio = img_nxbin / opt%img_nxbin
              opt%img_nyratio = img_nybin / opt%img_nybin
              WRITE(0,*) '  DONE',opt%img_nxbin,opt%img_nybin
            ENDIF
            WRITE(0,*) '  BUILDING PIXELS',opt%ndet
            CALL BuildPixels(MAXNPIXEL,npixel,pixel)
            WRITE(0,*) '  DONE BUILDING PIXELS'
          ELSE
            IF (opt%ndet.GT.0) opt%fmap = opt%det_fname(1)
          ENDIF
c...      Check for another detector:
c          READ(fp,'(A256)') dummy
c          BACKSPACE 5
c          WRITE(0,*) 'DUMMY:',dummy(1:10)
c          IF   (dummy(8:10).EQ.'Ccd'.OR.dummy(8:10).EQ.'CCD'.OR.
c     .          dummy(8:10).EQ.'ccd') THEN
c          ELSE
c            EXIT
c          ENDIF

c... This is all a right mess -- need to store the detector data even though it is not needed
c    once the pixels are all setup (true? some savings in a reference from each pixel element to the
c    data that was used to define the view?).  Hack for now...

        ENDDO

        opt%ccd = save_opt_ccd   ! Hack

      ELSE
c...    Need to store mask images, or at least the filenames, in opt%det_ 
c       as well and then load/bin them again... + buildpixels ...
        STOP 'NOT SURE WHAT TO DO HERE'
      ENDIF

c      STOP 'DONE MAN'  ! LEFT OFF

c      opt%ob_nsector = -1
c      opt%ob_invgrd  =  0

c      IF (opt%ob_model.GT.0) THEN
        WRITE(0,*) '  BUILDING 3D OBJECTS'
c        WRITE(0,*) '    toroidal extent ',
c     .    opt%ob_nsector,opt%ob_angle_start,
c     .    opt%ob_angle_end,opt%ob_yrotation
c        WRITE(0,*) '    standard  grid  ',opt%ob_stdgrd
c        WRITE(0,*) '    triangle  grid  ',opt%ob_trigrd
c        WRITE(0,*) '    inversion grid  ',opt%ob_invgrd
c        WRITE(0,*) '    vessel wall     ',opt%ob_wall
c        WRITE(0,*) '    targets         ',opt%ob_targ
c        WRITE(0,*) '    flux tube(s)    ',opt%ob_tube
c        WRITE(0,*) '    field line(s)   ',opt%ob_line
c        WRITE(0,*) '    user defined    ',opt%ob_user(1:10)
 
        CALL BuildObjects

        opt%obj_angle_start = 0.0
        opt%obj_angle_end = 360.0
        opt%obj_yrotation = 0.0
        opt%obj_nsector = -1

c...    Rotate all object vertices about the y-axis, to simulate when
c       the camera is located somewhere other than straight in along
c       the X or Z axes: 
        CALL Calc_Transform2(mat,0.0D0,1,0)
        CALL Calc_Transform2(mat,DBLE(opt%obj_yrotation*PI/180.0),2,1)
        DO ivtx = 1, nvtx
          CALL Transform_Vect(mat,vtx(1,ivtx))
        ENDDO

        WRITE(0,*) '  DONE BUILDING 3D OBJECTS'
c      ENDIF


      IF (opt%int_num.GT.0) 
     .  CALL AssignEmissionData(MAXNPIXEL,npixel,pixel)


      WRITE(0,*) '  NCHORD',nchord
      WRITE(0,*) '  NPIXEL',npixel
      WRITE(0,*) '  NOBJ  ',nobj
      WRITE(0,*) '  NSRF  ',nsrf
      WRITE(0,*) '  NVTX  ',nvtx

      IF (opt%ndet.GT.0) THEN

        WRITE(0,*) '  PROCESSING PIXELS'
        CALL ProcessPixels(npixel,pixel)
c        CALL ProcessPixels(opt,npixel,pixel,nobj,obj)
c        CALL ProcessPixels
        WRITE(0,*) '  DONE PROCESSING PIXELS'

        CALL DumpImage(npixel,pixel)
      ENDIF


      IF (.FALSE.) THEN
        DO i1 = 1, 10
          DO i2 = 1, 10
            x(i2+(i1-1)*10) = REAL(i2)
            y(i2+(i1-1)*10) = REAL(i1)
            f(i2+(i1-1)*10) = REAL(i2)
          ENDDO
          x1(i1) = REAL(i1)+0.5
          y1(i1) = 5.0
        ENDDO
c        WRITE(0,*) 'x:',x
c        WRITE(0,*) 'y:',y
c        WRITE(0,*) 'f:',f
        CALL Interpolate(100,x,y,f,9,x1,y1,f1,0)
        WRITE(0,*) 'x1:',x1
        WRITE(0,*) 'y1:',y1
        WRITE(0,*) 'F1:',f1
      ENDIF

      CALL Output985(iopt,MAXNPIXEL,npixel,pixel,image)

c      nobj985 = nobj

      IF (ALLOCATED(spectrum)) DEALLOCATE(spectrum)  ! temp?

      IF (ALLOCATED(plasma)) DEALLOCATE(plasma)

c...  Put into a subroutine:
      IF (ALLOCATED(obj))   DEALLOCATE(obj)
      IF (ALLOCATED(pixel)) DEALLOCATE(pixel)
      CALL DEALLOC_CHORD
      IF (ALLOCATED(vtx))   DEALLOCATE(vtx)
      IF (ALLOCATED(srf))   DEALLOCATE(srf)

c ... dealloc srf and vtx here as well, once plots are moved to output985, same for out989 

      RETURN
 99   STOP
      END