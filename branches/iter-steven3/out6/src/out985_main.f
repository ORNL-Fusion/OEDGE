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

      f1 = -999.0D0

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
c              WRITE(0,*) 'px,y:',px,py
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
      SUBROUTINE ProcessPixels(npixel,pixel,title)
      USE mod_out985
      USE mod_out985_variables
      USE mod_interface
      IMPLICIT none
c...  Input:
      INTEGER npixel
      TYPE(type_view) :: pixel(npixel)
      CHARACTER, INTENT(IN) :: title*(*)

      REAL Clock2

      REAL, PARAMETER :: PI = 3.141593

      INTEGER   i1,i2,i3,i4,status,fp,fp2,m,n,codec,iobj,ierr,volcnt,
     .          pixcnt,ndum1,npts,vcount,ipixel,nybin,nxbin,idet,npro
      LOGICAL   ascii,message
      REAL      rtime,ttime,fact
      REAL*8    int_sum,solid_total,solid_angle
      CHARACTER file*1024,tag*9,dummy*256

      INTEGER, ALLOCATABLE :: idum1(:)

      REAL*8, ALLOCATABLE, TARGET :: ddum1(:),ddum2(:,:)
      REAL*4, ALLOCATABLE, TARGET :: mapchk(:)
      REAL*4, ALLOCATABLE, TARGET :: rdum1(:)

      message = .TRUE.

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
      ALLOCATE(vwinter(MAXINTER))
      ALLOCATE(gbinter(MAXINTER))
      ALLOCATE(obinter(MAXINTER))
c...  Assemble a list of wall surfaces to speed things up:  (30% faster for current setup...)
      WRITE(0,*) '  BUILDING WALL SURFACE LISTS'  
      nvwlist = 0
      ngblist = 0
      DO i1 = 1, nobj
        DO i2 = 1, MAX(obj(i1)%nsur,obj(i1)%nside)
c...      Check if surface is part of the vessel wall:
          IF (obj(i1)%tsur(i2).EQ.SP_VESSEL_WALL  ) nvwlist = nvwlist+1
          IF (obj(i1)%tsur(i2).EQ.SP_GRID_BOUNDARY.OR.
     .        (obj(i1)%tsur(i2).EQ.SP_VESSEL_WALL       .AND.
     .         obj(i1)%type    .EQ.OP_INTEGRATION_VOLUME))
     .      ngblist = ngblist+1
        ENDDO
      ENDDO
      ALLOCATE(vwlist(nvwlist,2))
      ALLOCATE(gblist(ngblist,2))
      i3 = 0
      i4 = 0
      DO i1 = 1, nobj
        DO i2 = 1, MAX(obj(i1)%nsur,obj(i1)%nside)
c...      Check if surface is part of the vessel wall:
          IF     (obj(i1)%tsur(i2).EQ.SP_VESSEL_WALL) THEN
            i3 = i3 + 1
            vwlist(i3,1) = i1
            vwlist(i3,2) = i2
          ENDIF
          IF (obj(i1)%tsur(i2).EQ.SP_GRID_BOUNDARY.OR.
     .        (obj(i1)%tsur(i2).EQ.SP_VESSEL_WALL       .AND.
     .         obj(i1)%type    .EQ.OP_INTEGRATION_VOLUME)) THEN
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


      ALLOCATE(ddum1(n                          ))  ! MPI problem?  nobj=m, should be # integration volumes
      ALLOCATE(ddum2(nobj,-12:MAX(1,opt%int_num)))  ! MPI problem?  nobj=m, should be # integration volumes
      ALLOCATE(rdum1(100                        ))  ! MPI problem?  nobj=m, should be # integration volumes

      rtime = 0.0
      pixcnt = m

c...  MPI!  (for MPI to work the various lists -- vw,gb,ob -- need to be passed like the objects are passed, I think)

c...  Collect data, time and case information:
      CALL ZA09AS(dummy(1:8))
      dummy(9:10) = dummy(1:2)  ! Switch to EU formay
      dummy(1:2 ) = dummy(4:5)
      dummy(4:5 ) = dummy(9:10)
      dummy(9:10) = '  '
      CALL ZA08AS(dummy(11:18))
      CALL CASENAME(dummy(21:),ierr)

      vcount = 0

c...  Loop over pixels:
c      DO i1 = 200000, npixel
c      DO i1 = 4325, 4325

      idet = 1
      opt%ccd = opt%det_ccd(idet)
      WRITE(0,*) '================= IDET=',idet,'================='

      DO i1 = 1, npixel

c        WRITE(0,*) '---- pixel%v1 ---->',SNGL(pixel(i1)%v1)
c        WRITE(0,*) '---- pixel%v2 ---->',SNGL(pixel(i1)%v2)

c...    Keep track of which detector the current pixel belongs to:
        IF (i1.GT.opt%det_iend(idet)) THEN
          idet = idet + 1
          opt%ccd = opt%det_ccd(idet)
          WRITE(0,*) '================= IDET=',idet,'================='
        ENDIF

        IF (MOD(i1,MAX(1,npixel/10)).EQ.0) 
c        IF (MOD(i1,1).EQ.0) 
c        IF (MOD(i1,1000).EQ.0) 
     .    WRITE(0,'(A,I7,A,I7)') 'PROCESSING PIXEL ',i1,' OF ',npixel

        status = 0

        ddum1 = 0.0D0
        ddum2 = 0.0D0
        pixel(i1)%track   => ddum1  ! Dont' make part of view derived type definition, just send on its own?
        pixel(i1)%nprofile = 0
        pixel(i1)%profile => ddum2  ! Dont' make part of view derived type definition, just send on its own?

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
c        WRITE(0,*) 'INTEGRFAL:',i1,pixel(i1)%integral(1)


        IF (.FALSE..AND.   ! *** PROFILE HACK ***
     .      (i1.EQ.opt%det_istart(idet).OR.
     .       i1.EQ.opt%det_iend  (idet).OR.
     .       MOD(i1-opt%det_istart(idet)+1,5).EQ.0).AND.
     .      (opt%det_nxbin(1).EQ.1.OR.opt%det_nybin(1).EQ.1)) THEN
          WRITE(file,'(A,I0.3,256X)')          
     .      'idl.profile_'//TRIM(opt%det_fname(idet))//'_',
     .      i1-opt%det_istart(idet)+1
c          WRITE(0,*) ' FILE:'//TRIM(file)
          CALL inOpenInterface(TRIM(file),ITF_WRITE)
          npro=pixel(i1)%nprofile
         CALL inPutData(INT(pixel(i1)%profile(1:npro,-12)),'CELL','N/A')
          CALL inPutData(pixel(i1)%profile(1:npro,-5 ),'PATH'  ,'m'  )
          CALL inPutData(pixel(i1)%profile(1:npro,-4 ),'DELTA' ,'m'  )
          CALL inPutData(pixel(i1)%profile(1:npro,-3 ),'WEIGHT','N/A')
          CALL inPutData(pixel(i1)%profile(1:npro,-2 ),'NE'    ,'m-3')
          CALL inPutData(pixel(i1)%profile(1:npro,-1 ),'TE'    ,'eV' )
          CALL inPutData(pixel(i1)%profile(1:npro, 0 ),'TI'    ,'eV' )
          CALL inPutData(pixel(i1)%profile(1:npro,-7 ),'N_D'   ,'m-3')
          CALL inPutData(pixel(i1)%profile(1:npro,-6 ),'N_D2'  ,'m-3')
c          write(0,*) 'opt%int_num=',opt%int_num
          DO i2 = 1, MAX(1,opt%int_num)
            WRITE(tag,'(A,I0.2,A)') 'SIGNAL_',i2
            CALL inPutData(pixel(i1)%profile(1:npro,i2),TRIM(tag),'N/A')
          ENDDO
          CALL inPutData(opt%int_z(1:opt%int_num),'Z','N/A' )
          CALL inPutData(opt%int_a(1:opt%int_num),'A','N/A' )
          CALL inPutData(opt%int_charge(1:opt%int_num),'CHARGE','N/A')
          CALL inPutData(opt%int_wlngth(1:opt%int_num),'WLNGTH','nm')
          CALL inPutData(pixel(i1)%integral(1:opt%int_num),'INT','N/A')
          CALL inPutData(pixel(i1)%global_v1,'V1','m')
          CALL inPutData(pixel(i1)%global_v2,'V2','m')
          CALL inCloseInterface
c...      Dump to a more accessible ASCII file:
          WRITE(file,'(A,I0.3,256X)')          
     .      'profile.'//TRIM(opt%det_fname(idet))//'_',
     .      i1-opt%det_istart(idet)+1
          fp = 99

          OPEN(fp,FILE=file(1:LEN_TRIM(file)),
     .         FORM='FORMATTED',STATUS='REPLACE',ERR=97)
          WRITE(fp,'(A)') '*'
          WRITE(fp,'(A)') '* ---------------------------------------'//
     .                      '-------------------------------'//
     .                      '-------------------------------'
          WRITE(fp,'(A)') '* CASE            '//TRIM(dummy(21:))
          WRITE(fp,'(A)') '* TITLE           '//TRIM(title)
          WRITE(fp,'(A)') '* DATE AND TIME   '//dummy(1:18)
          WRITE(fp,'(A)') '*'
          WRITE(fp,'(A)') '* ---------------------------------------'//
     .                      '-------------------------------'//
     .                      '-------------------------------'
          WRITE(fp,'(A)') '{DATA FILE VERSION}'
          WRITE(fp,'(A)') '     1.0'
          WRITE(fp,'(A)') '*'
          WRITE(fp,'(A)') '* ---------------------------------------'//
     .                      '-------------------------------'//
     .                      '-------------------------------'
          WRITE(fp,'(A)') '{VIEW IDENTIFIER AND CHORD INDEX}'
          WRITE(fp,*) '  '//TRIM(opt%det_fname(idet)),
     .                i1-opt%det_istart(idet)+1
          WRITE(fp,'(A)') '* ---------------------------------------'//
     .                      '-------------------------------'//
     .                      '-------------------------------'
          WRITE(fp,'(A)') '{R,Z START AND END POINTS OF THE CHORD (m)}'
          WRITE(fp,'(3F10.4)') pixel(i1)%global_v1(1:2)
          WRITE(fp,'(3F10.4)') pixel(i1)%global_v2(1:2)
          WRITE(fp,'(A)') '*'
          WRITE(fp,'(A)') '* ---------------------------------------'//
     .                      '-------------------------------'//
     .                      '-------------------------------'
          WRITE(fp,'(A)') '{EMISSION SIGNALS}'   
          WRITE(fp,*) opt%int_num
          WRITE(fp,'(A)') '*'
          WRITE(fp,'(A)') '* Z          - atomic number of particle'
          WRITE(fp,'(A)') '* A          - atomic weight'
          WRITE(fp,'(A)') '* CHARGE     - electric charge'
          WRITE(fp,'(A)') '* WAVELENGTH - of the line'
          WRITE(fp,'(A)') '* INTEGRAL   - signal strength integrated '//
     .      'along the chord, in units of [photons ster-1 '//
     .      'm-2 s-1]'
          WRITE(fp,'(A)') '*'
          WRITE(fp,'(A2,A6,3A8,2A12)') 
     .      '* ','SIGNAL','Z'   ,'A'   ,'CHARGE','WAVELENGTH','INTEGRAL'
          WRITE(fp,'(A2,A6,3A8,2A12)')
     .      '* ','     ','    ','     ','      ','(nm)'      ,'        '
          fact = 4.0 * PI
          DO i2 = 1, opt%int_num
            WRITE(fp,'(2X,I6,3I8,F12.2,1P,E12.2,0P)') i2,
     .        opt%int_z         (i2)     ,
     .        opt%int_a         (i2)     ,
     .        opt%int_charge    (i2)     ,
     .        opt%int_wlngth    (i2)     ,
     .        pixel(i1)%integral(i2)/fact
          ENDDO
          WRITE(fp,'(A)') '*'
          WRITE(fp,'(A)') '* ---------------------------------------'//
     .                      '-------------------------------'//
     .                      '-------------------------------'
          WRITE(fp,'(A)') '{INTEGRATION VOLUMES}'
          WRITE(fp,*) npro
          WRITE(fp,'(A)') '*'
          WRITE(fp,'(A)') '* DIST   - distance along the chord '//
     .      '(line-of-sight)'
          WRITE(fp,'(A)') '* DELTA  - length of the path through the '//
     .      'current integation volume (triangular "cell")'
          WRITE(fp,'(A)') '* WEIGHT - set to unity for now, but will '//
     .      'be used when reflections are included in the analysis'
          WRITE(fp,'(A)') '* n_e    - electron density'
          WRITE(fp,'(A)') '* T_e    - electron temperature'
          WRITE(fp,'(A)') '* T_i    - background hydrogenic ion '//
     .      'temperature'
          WRITE(fp,'(A)') '* n_D    - deuterium atom density'
          WRITE(fp,'(A)') '* n_D2   - deuterium molecule density'
          WRITE(fp,'(A)') '* n_i+X  - density of impurity with charge X'
          WRITE(fp,'(A)') '* SIGNAL - local emission in units of '//
     .      '[photons ster-1 m-3 s-1]'
          WRITE(fp,'(A)') '*'
          WRITE(fp,'(A2,A6,4A10,2A8,6A10,10(A9,I1))') 
     .      '* ','INDEX','DIST','DELTA','WEIGHT','n_e'  ,'T_e' ,'T_i',
     .      'n_D','n_D2',
     .      'n_i+0','n_i+1','n_i+2','n_i+3',
     .      ('SIGNAL_',i2,i2=1,opt%int_num)
          WRITE(fp,'(A2,A6,4A10,2A8,6A10)') 
     .      '* ',' '    ,'(m)' ,'(m)'  ,' '     ,'(m-3)','(eV)','(eV)',
     .      ('(m-3)',i2=1,6)
50 	  FORMAT(50X,10I10)
51 	  FORMAT(50X,1P,10F10.2,0P)
          DO i2 = 1, npro
            WRITE(fp,'(2X,I6,1P,4E10.2,0P,2F8.1,1P,16E10.2,0P)') i2,
     .        pixel(i1)%profile(i2,-5),  ! path
     .        pixel(i1)%profile(i2,-4),  ! delta
     .        pixel(i1)%profile(i2,-3),  ! weight
     .        pixel(i1)%profile(i2,-2),  ! ne
     .        pixel(i1)%profile(i2,-1),  ! te
     .        pixel(i1)%profile(i2, 0),  ! ti
     .        pixel(i1)%profile(i2,-7),  ! nD
     .        pixel(i1)%profile(i2,-6),  ! nD2
     .        pixel(i1)%profile(i2,-11:-8),  !  4 lowest impurity charge states
     .       (pixel(i1)%profile(i2,i3)/fact,i3=1,opt%int_num)
          ENDDO
          CLOSE(fp)
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

      IF (.FALSE.) THEN
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
          IF (obj(iobj)%nsur.EQ.3) npts = 3  ! Toroidally continuous, triangle
          IF (obj(iobj)%nsur.EQ.4) npts = 4  ! Toroidally continuous, quadrangle ***WEAK***
          IF (obj(iobj)%nsur.EQ.5) npts = 3  ! Triangle             
          IF (obj(iobj)%nsur.EQ.6) npts = 4  ! Quadrelateral           
          IF (npts.EQ.0) THEN
            IF (message) THEN
              WRITE(0,*)
              WRITE(0,*) '----------------------------------------'
              WRITE(0,*) '     IMPROPER INVERSION CELL FOUND'
              WRITE(0,*) '        i,nobj= ',iobj,n
              WRITE(0,*) '        nsur  = ',obj(iobj)%nsur
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

      IF (.TRUE.) THEN
c...    Dump pixel view trajectories:
        DO idet = 1, opt%ndet
          fp = 99
          file = 'output.'//TRIM(opt%det_fname(idet))//'.ray.pxv'
          OPEN(fp,FILE=file(1:LEN_TRIM(file)),
     .         FORM='FORMATTED',STATUS='REPLACE',ERR=97)
          nxbin = opt%det_nxbin(idet)
          nybin = opt%det_nybin(idet)
          WRITE(fp,'(A)') '* Detector pixel views'
          WRITE(fp,'(2I8)') nxbin,nybin
          DO ipixel = opt%det_istart(idet), opt%det_iend(idet)
            WRITE(fp,'(2I8,2X,2(3F10.6,2X))')
     .       pixel(ipixel)%xindex,pixel(ipixel)%yindex,
     .       pixel(ipixel)%global_v1(1:3),
     .       pixel(ipixel)%global_v2(1:3)
          ENDDO
          CLOSE(fp)
        ENDDO
      ENDIF

      IF (.TRUE.) THEN
c...    Total the signal for all pixels:
        solid_total = (opt%angle(1) * opt%angle(2)) / (180.D0 * 360.0D0)
        IF (solid_total.GT.1.0D0) 
     .    CALL ER('ProcessPixels','Solid angle greater than 4PI',*99)
        solid_angle = solid_total / DBLE(npixel)
        DO i1 = 1, MAX(1,opt%int_num)
          int_sum = 0.0
          DO i2 = 1, npixel
            int_sum = int_sum + pixel(i2)%integral(i1) * solid_angle /
     .                          6.0D0 *  ! 8.0D0 * ! 6.0D0
     .                   DCOS(pixel(i2)%xangle * 3.1415D0 / 180.0D0) *
     .                   DCOS(pixel(i2)%yangle * 3.1415D0 / 180.0D0)


c            pixel(i2)%integral(i1) =  pixel(i2)%integral(i1) *
c     .                   DCOS(pixel(i2)%xangle * 3.1415D0 / 180.0D0) *
c     .                   DCOS(pixel(i2)%yangle * 3.1415D0 / 180.0D0)
 
          ENDDO
c          WRITE(0,*)
c          WRITE(0,*) 'INTEGRAL SUMMATION:',i1,int_sum
c          WRITE(0,*)
        ENDDO
      ENDIF

c...  Clear arrays:
      DEALLOCATE(idum1)
      DEALLOCATE(ddum1)
      DEALLOCATE(ddum2)
      DEALLOCATE(rdum1)
      DEALLOCATE(vwlist)
      DEALLOCATE(gblist)
      DEALLOCATE(vwinter)
      DEALLOCATE(gbinter)
      DEALLOCATE(obinter)
      DEALLOCATE(mapchk)

      RETURN
 97   CALL ER('ProcessPixels','Unable to save pixel views',*99)
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

c...  Detector parameters:
      debug = .FALSE.

      IF (debug) THEN
c...    Debug:
        xcen = 2.5D0
        ycen = 0.0D0
        zcen = 0.0D0
        xwidth = 0.20D-07
        ywidth = 0.20D-07
        xangle = -13.90099009900990D0
        yangle = -25.66336633663366D0

c...    Number of rows and columns of pixels:
        nxbin =  1
        nybin =  1
      ELSEIF (.NOT..TRUE.) THEN
        xcen = 2.5D0
        ycen = 0.0D0
        zcen = 0.0D0
c        xwidth = 0.20D-07
c        ywidth = 0.20D-07
        xwidth = 1.80D0
        ywidth = 1.80D0
c        xwidth = 0.20D0
c        ywidth = 0.20D0
c        xangle = 22.5D0
c        yangle = 22.5D0
c        xangle = 65.0D0
c        yangle = 65.0D0
        xangle = 0.0D0
        yangle = 0.0D0
c        xangle = 54.0D0
c        yangle = 54.0D0

c...    Number of rows and columns of pixels:
        nxbin =  1
        nybin =  1
      ELSEIF (.NOT..TRUE.) THEN
        xcen = 2.15D0
        ycen = 1.42D0
        zcen = 0.0D0
        xwidth = 0.20D-07
        ywidth = 0.20D-07
c        xwidth = 0.20D0
c        ywidth = 0.20D0
c        xangle = 22.5D0
c        yangle = 22.5D0
c        xangle = 54.0D0
c        yangle = 54.0D0
        xangle = 65.0D0
        yangle = 65.0D0

c...    Number of rows and columns of pixels:
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
c...    Number of rows and columns of pixels:
        nxbin = opt%nxbin
        nybin = opt%nybin
      ENDIF

      WRITE(0,'(A,3F10.3)') 'DATA:',xcen,ycen,zcen
      WRITE(0,'(A,3F10.3)') '    :',opt%roll,opt%tilt,
     .                              opt%swing
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
c          v = DATAN( MAX(0.5D0, 0.0D0 * x**2 + 1.0D0) * x / f )
c          v = DATAN( MAX(0.5D0, dist * x**2 + 1.0D0) * x / f )
          v = DATAN( x / MAX(0.5D0, dist * x**2 + 1.0D0) / f )
          IF     (DABS(v-thetax).LT.1.0D-08.OR.count.EQ.1000) THEN
c            WRITE(0,*) 'WORKS X!',x,thetax*180.0D0/PI
            boundx = x
            EXIT
          ELSEIF ((v.GT.thetax.AND.deltax.GT.0.0D0).OR.
     .            (v.LT.thetax.AND.deltax.LT.0.0D0)) THEN
            deltax = -0.3D0 * deltax
          ENDIF
          x = x + deltax
        ENDDO

        y = 0.0D0
        deltay = 0.01D0
        count = 0
        DO WHILE(.TRUE.) 
          count = count + 1
c          v = DATAN( MAX(0.5D0, 0.0D0 * y**2 + 1.0D0) * y / f )
c          v = DATAN( MAX(0.5D0, dist * y**2 + 1.0D0) * y / f )
          v = DATAN( y / MAX(0.5D0, dist * y**2 + 1.0D0) / f )
          IF (DABS(v-thetay).LT.1.0D-08.OR.count.EQ.1000) THEN
c            WRITE(0,*) 'WORKS Y!',y
            boundy = y
            EXIT
          ELSEIF ((v.GT.thetay.AND.deltay.GT.0.0D0).OR.
     .            (v.LT.thetay.AND.deltay.LT.0.0D0)) THEN
            deltay = -0.3D0 * deltay
          ENDIF
          y = y + deltay
        ENDDO
      ENDIF

c...  Store index locations of first pixel for this detector:
      opt%det_istart(opt%ndet) = npixel + 1

      IF (opt%chord_n.GT.0) THEN 

        DO i1 = 1, opt%chord_n 
          npixel = npixel + 1
          pixel(npixel)%valid    = .TRUE.
          pixel(npixel)%weight   = 1.0D0
          pixel(npixel)%integral = 0.0D0
          pixel(npixel)%index    = opt%ndet
          pixel(npixel)%type     = 2
          pixel(npixel)%ccd      = opt%chord_opt(i1)
          pixel(npixel)%tmp      = npixel    !*** TEMPORARY ***
          pixel(npixel)%xindex   = i1
          pixel(npixel)%yindex   = 1
          pixel(npixel)%nxbin    = 1 ! opt%sa_nxbin
          pixel(npixel)%nybin    = 1 ! opt%sa_nybin
          pixel(npixel)%v1(1:3)  = opt%chord_v1(1:3,i1)
          pixel(npixel)%v2(1:3)  = opt%chord_v2(1:3,i1)

          pixel(npixel)%xwidth   = 0.0D0
          pixel(npixel)%ywidth   = 0.0D0
          pixel(npixel)%xangle   = 0.0D0 
          pixel(npixel)%yangle   = 0.0D0
          SELECTCASE (opt%sa_opt)
           CASE (1) 
            pixel(npixel)%dxangle=xangle/DBLE(nxbin)*DBLE(opt%sa_par1)     ! View width 
            pixel(npixel)%dyangle=yangle/DBLE(nybin)*DBLE(opt%sa_par2)         
           CASE DEFAULT
            CALL ER('BuildPixels','Unknown solid angle option (1)',*99)
          ENDSELECT
          pixel(npixel)%rot      = 0.0D0
          pixel(npixel)%trans    = 0.0D0
        ENDDO
      ELSE
        DO iy = 1, nybin
          DO ix = 1, nxbin
            npixel = npixel + 1
        
            pixel(npixel)%valid    = .TRUE.
            pixel(npixel)%weight   = 1.0D0
            pixel(npixel)%integral = 0.0D0
            pixel(npixel)%index    = opt%ndet
            pixel(npixel)%type     = 1
            pixel(npixel)%tmp      = npixel    !*** TEMPORARY ***
            pixel(npixel)%xindex   = ix
            pixel(npixel)%yindex   = iy
            pixel(npixel)%nxbin    = opt%sa_nxbin
            pixel(npixel)%nybin    = opt%sa_nybin
        
            factx = -1.0D0 * (0.5D0 * (1.0D0 / DBLE(nxbin) - 1.0D0) +
     .                        DBLE(ix - 1) / DBLE(nxbin))
            facty = -1.0D0 * (0.5D0 * (1.0D0 / DBLE(nybin) - 1.0D0) +
     .                        DBLE(iy - 1) / DBLE(nybin))
        
            pixel(npixel)%v1(1)  = factx * xwidth            ! Spatial horizontal location 
            pixel(npixel)%v1(2)  = facty * ywidth 
            pixel(npixel)%v1(3)  = 0.0D0
            pixel(npixel)%xwidth = xwidth / DBLE(nxbin)      ! Spatial horizontal width
            pixel(npixel)%ywidth = ywidth / DBLE(nybin)            
        
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
              CALL ER('BuildPixels','Unknown solid angle option (2',*99)
            ENDSELECT
        
            IF (debug) THEN
              pixel(npixel)%xangle = xangle            ! View angle            
              pixel(npixel)%yangle = yangle            ! View angle
            ENDIF
        
c...        These are only needed for plotting, not actually used during LOS integration:
            pixel(npixel)%v2(1) = DSIN(pixel(npixel)%xangle *D_DEGRAD) *    ! Correct?
     .                            DCOS(pixel(npixel)%yangle *D_DEGRAD) * 
     .                            50.0D0 + pixel(npixel)%v1(1)                ! 7 meters/units?
            pixel(npixel)%v2(2) = DCOS(pixel(npixel)%xangle *D_DEGRAD) * 
     .                            DSIN(pixel(npixel)%yangle *D_DEGRAD) * 
     .                            50.0D0 + pixel(npixel)%v1(2)
            pixel(npixel)%v2(3) = DCOS(pixel(npixel)%xangle *D_DEGRAD) * 
     .                            DCOS(pixel(npixel)%yangle *D_DEGRAD) * 
     .                            50.0D0 + pixel(npixel)%v1(3)
c            pixel(npixel)%v2(3) = -pixel(npixel)%v2(3)
        
c           *** better to store these in a detector array, rather than taking up space with each pixel, or are the savings nominal...
c           *** use opt%det_ array...
            pixel(npixel)%rot(1) = DBLE(opt%tilt *PI/180.0D0) ! x-axis (tilt) 
            pixel(npixel)%rot(2) = DBLE(opt%swing*PI/180.0D0) ! y-axis (swing)
            pixel(npixel)%rot(3) = DBLE(opt%roll *PI/180.0D0) ! z-axis (roll) ! Set detector orientation
        
            pixel(npixel)%trans(1) = xcen
            pixel(npixel)%trans(2) = ycen
            pixel(npixel)%trans(3) = zcen
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

      IF (opt%obj_num.NE.0) THEN

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
            CASE (8)
              CALL rayLoadITERFWP         (ielement)
            CASE DEFAULT
              CALL User_CustomObjects     (ielement)
          ENDSELECT
        ENDDO

      ENDIF



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
      USE mod_interface
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
c...      Also generate mod_interface file:
          WRITE(file,'(1024X)')          
          file = 'output.'//TRIM(opt%det_fname(idet))//'.ray.img'
          CALL inOpenInterface(TRIM(file),ITF_WRITE)
          DO iy = 1, nybin
            CALL inPutData(image(1:nxbin,iy),'data','unknown')
          ENDDO
          CALL inPutData(nxbin,'dx','none')
          CALL inPutData(nybin,'dy','none')
          CALL inCloseInterface
        ENDIF

        IF (.TRUE.) THEN
          fp = 99
          WRITE(file,'(A,I1,A)')
          file = 'output.'//TRIM(opt%det_fname(idet))//'.ray.img_old'
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
c subroutine: Main985
c
      SUBROUTINE Main985(iopt,title)
      USE mod_out985
      USE mod_out985_variables
      USE mod_out985_clean
      IMPLICIT none

      INTEGER  , INTENT(IN) :: iopt
      CHARACTER, INTENT(IN) :: title*(*)

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

c     Just in case there are previous allocations from OUT987 tetrahedron plots:
      CALL Wrapper_ClearObjects  
      CALL DEALLOC_ALL

      WRITE(0,*) '  ALLOCATING OBJECTS'
c      MAX3D = 10000000
c      MAX3D = 4000000 
c      MAX3D = 500000
       MAX3D = 10000
c      MAX3D = 1500000 
      ALLOCATE(obj(MAX3D))

      CALL ALLOC_SURFACE(-1,MP_INITIALIZE)

      MAXNPIXEL= 550*550 ! 480*640 ! 1000*1000
      WRITE(0,*) '  ALLOCATING PIXELS',MAXNPIXEL
      ALLOCATE(pixel(MAXNPIXEL))

c      CALL ALLOC_CHORD(MAXNPIXEL)  ! Just for viewing! (make smaller!)
      CALL ALLOC_CHORD(22500)  ! Just for viewing! (make smaller!)

      opt%load = 1

      grd_ntorseg = 48

      WRITE(0,*) 'MAXNPIXEL:',maxnpixel

      fp = 5

      npixel = 0
      nchord = 0
      n_pchord = 0

      opt%ref_opt = 0

      opt%img_nxratio = 1
      opt%img_nyratio = 1


      opt%nplots  = -1
      opt%chord_n =  0

      IF (opt%load.EQ.1) THEN
        opt%obj_num = 0
        opt%int_num = 0
        opt%ref_num = 0
        opt%ndet    = 0
        opt%rib_n   = 0

        CALL LoadOptions985_New(opt,ALL_OPTIONS,status)
      ENDIF

      IF (opt%load.EQ.1) THEN    ! *** HACK *** This is crap, i.e. the {DETECTOR} tags
        status = 0               ! have to come last in the input file.  Also, the 
        save_opt_ccd = 0         ! focal length, distortion, etc. are only recorded for 
        DO WHILE(status.EQ.0)    ! the last active detector.  Need to define TYPE_DETECTOR
          opt%ccd = 0  ! Hack    ! and use in in TYPE_OPTION985...
          opt%distortion  = 0.0      
          opt%focallength = 0.0
          opt%mask_num    = 0
          CALL LoadOptions985_New(opt,DETECTOR_ONLY,status)
          IF (opt%ccd.GT.0) THEN

            IF (save_opt_ccd.EQ.0) save_opt_ccd = opt%ccd

            opt%ndet = opt%ndet + 1
            opt%det_ccd  (opt%ndet) = opt%ccd  
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


      WRITE(0,*) '  BUILDING 3D OBJECTS'
c      WRITE(0,*) '    toroidal extent ',
c     .  opt%ob_nsector,opt%ob_angle_start,
c     .  opt%ob_angle_end,opt%ob_yrotation
c      WRITE(0,*) '    standard  grid  ',opt%ob_stdgrd
c      WRITE(0,*) '    triangle  grid  ',opt%ob_trigrd
c      WRITE(0,*) '    inversion grid  ',opt%ob_invgrd
c      WRITE(0,*) '    vessel wall     ',opt%ob_wall
c      WRITE(0,*) '    targets         ',opt%ob_targ
c      WRITE(0,*) '    flux tube(s)    ',opt%ob_tube
c      WRITE(0,*) '    field line(s)   ',opt%ob_line
c      WRITE(0,*) '    user defined    ',opt%ob_user(1:10)
 
      CALL BuildObjects

      opt%obj_angle_start = 0.0
      opt%obj_angle_end   = 360.0
      opt%obj_yrotation   = 0.0
      opt%obj_nsector     = -1
c...  Rotate all object vertices about the y-axis, to simulate when
c     the camera is located somewhere other than straight in along
c     the X or Z axes: 
      CALL Calc_Transform2(mat,0.0D0,1,0)
      CALL Calc_Transform2(mat,DBLE(opt%obj_yrotation*PI/180.0),2,1)
      DO ivtx = 1, nvtx
        CALL Transform_Vect(mat,vtx(1,ivtx))
      ENDDO

      WRITE(0,*) '  DONE BUILDING 3D OBJECTS'


      WRITE(0,*) '  NCHORD',nchord
      WRITE(0,*) '  NPIXEL',npixel
      WRITE(0,*) '  NOBJ  ',nobj
      WRITE(0,*) '  NSRF  ',nsrf
      WRITE(0,*) '  NVTX  ',nvtx

c...
      IF (opt%int_num.GT.0) 
     .  CALL AssignEmissionData(MAXNPIXEL,npixel,pixel)

c...
      IF (opt%ndet.GT.0) THEN
        WRITE(0,*) '  PROCESSING PIXELS'
        CALL ProcessPixels(npixel,pixel,title)
        WRITE(0,*) '  DONE PROCESSING PIXELS'
        CALL DumpImage(npixel,pixel)
      ENDIF

c...
      IF (opt%rib_n.GT.0) THEN
        WRITE(0,*) '  PROCESSING RIBBON'
        CALL rayGenerateRibbonGrid
        WRITE(0,*) '  DONE PROCESSING RIBBON'
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

      IF (opt%nplots.GT.0)
     .  CALL Output985(iopt,MAXNPIXEL,npixel,pixel,image)

      IF (.TRUE.)
     .  CALL DumpShinKajita(title)

c      nobj985 = nobj

      IF (ALLOCATED(spectrum)) DEALLOCATE(spectrum)  ! temp?

      IF (ALLOCATED(plasma)) DEALLOCATE(plasma)

c...  Put into a subroutine:
      nobj = 0
      IF (ALLOCATED(obj))   DEALLOCATE(obj)
      IF (ALLOCATED(pixel)) DEALLOCATE(pixel)
      CALL DEALLOC_CHORD
      IF (ALLOCATED(vtx))   DEALLOCATE(vtx)
      IF (ALLOCATED(srf))   DEALLOCATE(srf)

c ... dealloc srf and vtx here as well, once plots are moved to output985, same for out989 

      RETURN
 99   STOP
      END
