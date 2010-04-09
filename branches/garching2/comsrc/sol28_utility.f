c     -*-Fortran-*-
c
c ======================================================================
c
      SUBROUTINE SaveGrid(fname)
      USE mod_sol28_global
      IMPLICIT none

      CHARACTER*(*) fname
      INTEGER fp,i1,ion

      fp = 99
      OPEN(UNIT=fp,FILE=TRIM(fname),ACCESS='SEQUENTIAL',
     .     FORM='UNFORMATTED',STATUS='REPLACE',ERR=97)            
      WRITE(fp,ERR=98) 1.00
      WRITE(fp,ERR=98) grid
      WRITE(fp,ERR=98) nion,ntube,ncell,nfield,npin,nphoton,ndrift,
     .                 nkinetic,nfluid,nimpurity
      WRITE(fp,ERR=98) (tube (i1),i1=1,ntube)
      WRITE(fp,ERR=98) (cell (i1),i1=1,ncell)
      WRITE(fp,ERR=98) (field(i1),i1=1,nfield)
      DO ion = 1, nion
        WRITE(fp,ERR=98) (pin     (i1,ion),i1=1,npin)
        WRITE(fp,ERR=98) (photon  (i1,ion),i1=1,nphoton)
        WRITE(fp,ERR=98) (drift   (i1,ion),i1=1,ndrift)
        WRITE(fp,ERR=98) (kinetic (i1,ion),i1=1,nkinetic)
        WRITE(fp,ERR=98) (fluid   (i1,ion),i1=1,nfluid)
        WRITE(fp,ERR=98) (impurity(i1,ion),i1=1,nimpurity)
      ENDDO
      CLOSE (fp)

c      WRITE(0,*) 'REF:',fluid(1:100,1)%ne

      RETURN
 97   CALL ER('SaveGrid','Problem creating raw data file',*99)
 98   CALL ER('SaveGrid','Problem writing to raw data file',*99)
 99   WRITE(0,*) '  FILE NAME = ',TRIM(fname)
      STOP
      END
c
c ======================================================================
c
      SUBROUTINE LoadGrid(fname)
      USE mod_sol28_params
      USE mod_sol28_global
      IMPLICIT none

      INTEGER i1,fp,ion     ,itarget,itube
      REAL    vgrid  ! Version number
      CHARACTER*(*) fname

      fp = 99
      OPEN(UNIT=fp,FILE=fname(1:LEN_TRIM(fname)),ACCESS='SEQUENTIAL',
     .     FORM='UNFORMATTED',STATUS='OLD',ERR=98)            
      READ (fp,ERR=98) vgrid
      READ (fp,ERR=98) grid
c...  Check version numbers:
      READ (fp,ERR=98) nion,ntube,ncell,nfield,npin,nphoton,ndrift,
     .                 nkinetic,nfluid,nimpurity
      IF (ALLOCATED(tube    )) DEALLOCATE(tube    )
      IF (ALLOCATED(cell    )) DEALLOCATE(cell    )
      IF (ALLOCATED(field   )) DEALLOCATE(field   )
      IF (ALLOCATED(pin     )) DEALLOCATE(pin     )
      IF (ALLOCATED(photon  )) DEALLOCATE(photon  )
      IF (ALLOCATED(drift   )) DEALLOCATE(drift   )
      IF (ALLOCATED(kinetic )) DEALLOCATE(kinetic )
      IF (ALLOCATED(fluid   )) DEALLOCATE(fluid   )
      IF (ALLOCATED(impurity)) DEALLOCATE(impurity)
      ALLOCATE(tube    (ntube ))
      ALLOCATE(cell    (ncell ))
      ALLOCATE(field   (nfield))
      ALLOCATE(pin     (npin     ,nion))
      ALLOCATE(photon  (nphoton  ,nion))
      ALLOCATE(drift   (ndrift   ,nion))
      ALLOCATE(kinetic (nkinetic ,nion))
      ALLOCATE(fluid   (nfluid   ,nion))
      ALLOCATE(impurity(nimpurity,nion))
      READ (fp,ERR=98) (tube (i1),i1=1,ntube)
      READ (fp,ERR=98) (cell (i1),i1=1,ncell)
      READ (fp,ERR=98) (field(i1),i1=1,nfield)
      DO ion = 1, nion
        READ (fp,ERR=98) (pin     (i1,ion),i1=1,npin)
        READ (fp,ERR=98) (photon  (i1,ion),i1=1,nphoton)
        READ (fp,ERR=98) (drift   (i1,ion),i1=1,ndrift)
        READ (fp,ERR=98) (kinetic (i1,ion),i1=1,nkinetic)
        READ (fp,ERR=98) (fluid   (i1,ion),i1=1,nfluid)
        READ (fp,ERR=98) (impurity(i1,ion),i1=1,nimpurity)
      ENDDO

      CLOSE (fp)

c      WRITE(0,*) 'NTUBE:',ntube,nfluid
c      WRITE(0,*) 'REF LOAD:',fluid(1:100,1)%ne
c      STOP 'sdgfsdsdg'


      IF (.FALSE.) THEN
        ion = 1

        tube(1:ntube)%jsat  (LO,ion) = 0.0
        tube(1:ntube)%machno(LO)     = 0.0
        tube(1:ntube)%gamma (LO,ion) = 0.0
        tube(1:ntube)%ne    (LO)     = 0.0
        tube(1:ntube)%te    (LO)     = 0.0
        tube(1:ntube)%ni    (LO,ion) = 0.0
        tube(1:ntube)%vi    (LO,ion) = 0.0
        tube(1:ntube)%ti    (LO,ion) = 0.0

        tube(1:ntube)%jsat  (HI,ion) = 0.0
        tube(1:ntube)%machno(HI)     = 0.0
        tube(1:ntube)%gamma (HI,ion) = 0.0
        tube(1:ntube)%ne    (HI)     = 0.0
        tube(1:ntube)%te    (HI)     = 0.0
        tube(1:ntube)%ni    (HI,ion) = 0.0
        tube(1:ntube)%vi    (HI,ion) = 0.0
        tube(1:ntube)%ti    (HI,ion) = 0.0

        fluid(1:nfluid,ion)%ne = 0.0
        fluid(1:nfluid,ion)%te = 0.0
        fluid(1:nfluid,ion)%ni = 0.0
        fluid(1:nfluid,ion)%vi = 0.0
        fluid(1:nfluid,ion)%ti = 0.0
      ENDIF

c        WRITE(0,*) 'TARGET DATA 1:'
c        DO itarget = LO, HI
c          WRITE(0,*)
c          DO itube = 1, ntube
c            WRITE(0,'(I6,1P,E16.6,0P,2X,F6.2,2X,2F10.6,2X,F8.2)')
c     .        itube,
c     .        tube(itube)%jsat  (itarget,ion),
c     .        tube(itube)%machno(itarget),
c     .        tube(itube)%te    (itarget),
c     .        tube(itube)%ti    (itarget,ion),
c     .        tube(itube)%gamma (itarget,ion)
c          ENDDO
c        ENDDO

      RETURN
 98   CALL ER('LoadGrid','Problem loading OSM grid file, maybe '//
     .        'MOD_SOL28 has changed?',*99)
 99   WRITE(0,*) '  FILE NAME='//'"'//fname(1:LEN_TRIM(fname))//'"'
      STOP
      END
c
c ======================================================================
c
      SUBROUTINE osmClean    ! Move to mod_osm.f
      USE mod_sol28_global
      USE mod_sol28_wall
      IMPLICIT none

c...  Close streaming output files: 
c      IF (logfp.GT.0) CLOSE(logfp)

c...  Solution arrays:
      ntube     = 0
      ncell     = 0
      nfield    = 0
      npin      = 0
      nphoton   = 0
      ndrift    = 0
      nkinetic  = 0
      nfluid    = 0
      nimpurity = 0
      IF (ALLOCATED(tube    )) DEALLOCATE(tube    )
      IF (ALLOCATED(cell    )) DEALLOCATE(cell    )
      IF (ALLOCATED(field   )) DEALLOCATE(field   )
      IF (ALLOCATED(pin     )) DEALLOCATE(pin     )
      IF (ALLOCATED(photon  )) DEALLOCATE(photon  )
      IF (ALLOCATED(drift   )) DEALLOCATE(drift   )
      IF (ALLOCATED(kinetic )) DEALLOCATE(kinetic )
      IF (ALLOCATED(fluid   )) DEALLOCATE(fluid   )
      IF (ALLOCATED(impurity)) DEALLOCATE(impurity)

c...  Geometry:
      nwall = 0
      IF (ALLOCATED(wall)) DEALLOCATE(wall)

c...  Reference solution:
      ref_nion   = 0
      ref_ntube  = 0
      ref_ncell  = 0
      ref_nfluid = 0
      IF (ALLOCATED(ref_tube )) DEALLOCATE(ref_tube) 
      IF (ALLOCATED(ref_fluid)) DEALLOCATE(ref_fluid) 
      IF (ALLOCATED(ref_cell )) DEALLOCATE(ref_cell) 

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      REAL*8 FUNCTION GetNodeParticleFlux(inode,ion)
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none

      INTEGER, INTENT(IN) :: inode,ion

      REAL GetCs2

      INTEGER ic
      LOGICAL success
      REAL ne_node,vb_node,te_node,ti_node,flux

      flux = 0.0D0

      ic = node(inode)%icell

      success = .TRUE.

      IF     (inode.EQ.1) THEN
        flux = isat(ictarg(LO),ion)
      ELSEIF (inode.EQ.nnode) THEN
        flux = isat(ictarg(HI),ion)
      ELSE
        IF     (node(inode)%ne.NE.0.0) THEN
          ne_node = node(inode)%ne
        ELSEIF (ne(ic).NE.0.0D0) THEN
          ne_node = SNGL(ne(ic))
        ELSE
          success = .FALSE.          
        ENDIF
        IF     (node(inode)%te.NE.0.0) THEN
          te_node = node(inode)%te
        ELSEIF (te(ic).NE.0.0D0) THEN
          te_node = SNGL(te(ic))
        ELSE
          success = .FALSE.          
        ENDIF
        IF     (node(inode)%ti(ion).NE.0.0) THEN
          ti_node = node(inode)%ti(ion)
        ELSEIF (ti(ic,ion).NE.0.0D0) THEN
          ti_node = SNGL(ti(ic,ion))
        ELSE
          success = .FALSE.          
        ENDIF

        IF (success) THEN
          IF     (node(inode)%v     .NE.0.0) THEN
            vb_node = node(inode)%v
          ELSEIF (node(inode)%machno.NE.0.0) THEN
            vb_node = node(inode)%machno * GetCs2(te_node,ti_node)
          ELSE
            CALL ER('GetNodeParticleFlux','Anomalous particle '//
     .              'source set to match specified velocity, '//
     .              'but no velocity information found',*99)
          ENDIF
          flux = ne_node * vb_node 
        ENDIF
      ENDIF

      GetNodeParticleFlux = DBLE(flux)

      RETURN
 99   WRITE(0,*) '  NODE= ',inode
      STOP
      END
c
c ======================================================================
c
c function: GetCellPressure
c
c     Returns the total pressure of a grid cell in units of [Pa].
c
      REAL FUNCTION GetCellPressure(icell,mode)
      USE mod_sol28_global
      IMPLICIT none

      INTEGER, INTENT(IN) :: icell,mode

      REAL, PARAMETER :: ECH = 1.6022E-19, AMU = 1.67E-27

      INTEGER ion
      REAL    m,n,v,pe,pi,te,ti

      ion = 1

      m = 2.0 * AMU  ! *** MASS HARDCODED! ***

      n  = fluid(icell,ion)%ne
      v  = fluid(icell,ion)%vi
      te = fluid(icell,ion)%te 
      ti = fluid(icell,ion)%ti 

      pe = n *  te                     * ECH
      pi = n * (ti + (m / ECH) * v**2) * ECH

      SELECTCASE (mode)
        CASE (1)
          GetCellPressure = pe    
        CASE (2)
          GetCellPressure = pi    
        CASE (3)
          GetCellPressure = pe + pi    
        CASEDEFAULT
          CALL ER('GetCellPressure','Unrecognised MODE',*99)
      ENDSELECT

      RETURN
99    STOP
      END
c
c ======================================================================
c
      REAL*8 FUNCTION GetNodePressure(inode,ion)
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none

      INTEGER, INTENT(IN) :: inode,ion

      REAL GetCs2

      INTEGER target
      REAL    ne_node,vb_node,te_node,ti_node,pressure

      pressure = 0.0

      IF     (inode.EQ.1) THEN
        pressure = pe(ictarg(LO)) + pi(ictarg(LO),1)
      ELSEIF (inode.EQ.nnode) THEN
        pressure = pe(ictarg(HI)) + pi(ictarg(HI),1)
      ELSEIF (node(inode)%ne.NE.0.0.AND.
     .        node(inode)%pe.NE.0.0) THEN
        STOP 'NOT SURE WHAT TO DO'
      ELSE
        target = LO
        IF (inode.GT.mnode) target = HI
 
        IF ( opt%bc(target).EQ.3.0.OR.
     .      (opt%bc(target).EQ.1.0.AND.
     .       (node(inode)%te.LE.0.001.OR.node(inode)%ti(ion).LE.0.001)))
     .    THEN
          te_node = te(node(inode)%icell)
          ti_node = ti(node(inode)%icell,ion)
        ELSE
          te_node = node(inode)%te
          ti_node = node(inode)%ti(ion)
        ENDIF
        IF (te_node.LE.0.001.OR.ti_node.LE.0.001) 
     .    CALL ER('GetNodePressure','Temperature data not found',*99)

        IF     (node(inode)%v     .NE.0.0) THEN
          vb_node = node(inode)%v
        ELSEIF (node(inode)%machno.NE.0.0) THEN
          vb_node = node(inode)%machno * GetCs2(te_node,ti_node)
        ELSE
c...      Ignore flow for now:
          vb_node = 0.0
        ENDIF

c        ti_node = te_node

        IF (node(inode)%ne.GT.0.0) THEN
          SELECTCASE (0)
            CASE (0)
              ne_node = node(inode)%ne
              pressure = ne_node * ((te_node + ti_node) * SNGL(ECH) +  
     .                              SNGL(mi(ion)) * vb_node**2)
            CASEDEFAULT                                           
          ENDSELECT                                                         
        ELSEIF (node(inode)%pe.GT.0.0) THEN
          SELECTCASE (0)
            CASE (0)
              pressure = 2.0 * node(inode)%pe * SNGL(ECH)
            CASEDEFAULT 
          ENDSELECT
        ELSE
          CALL ER('GetNodePressure','Unknown situation',*99)
        ENDIF
      ENDIF

      GetNodePressure = DBLE(pressure)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      INTEGER FUNCTION GetNumberOfObjects(filename)
      IMPLICIT none

      CHARACTER, INTENT(IN) :: filename*(*)

      INTEGER   fp,idum1
      CHARACTER buffer*256,fname*512

      fname = TRIM(filename)
      IF (TRIM(filename).EQ.'default') fname = 'eirene.transfer'      

      fp = 99
      OPEN(UNIT=fp,FILE=TRIM(fname),ACCESS='SEQUENTIAL',
     .     STATUS='OLD',ERR=98)
      DO WHILE (.TRUE.)
        READ(fp,'(A256)',END=10) buffer
        IF     (buffer(1:16).EQ.'* BULK PARTICLES') THEN
          READ(fp,*,ERR=97) 
          READ(fp,*,ERR=97) GetNumberOfObjects
          CLOSE(fp)
          RETURN
        ENDIF
      ENDDO

 10   CALL ER('GetNumberOfObjects','EOF error',*99)

      RETURN
 97   CALL ER('GetNumberOfObjects','File error #1',*99)
 98   CALL ER('GetNumberOfObjects','File error #2',*99)
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE LoadTriangles_06
      USE mod_eirene06
      IMPLICIT none

      INTEGER fp,i1,i2
      REAL    version

      fp = 99
      OPEN(UNIT=fp,FILE='triangles.raw',ACCESS='SEQUENTIAL',
     .     FORM='UNFORMATTED',STATUS='OLD',ERR=98)            
      READ(fp,ERR=98) version,ntri,nver,nsurface

      IF (version.NE.1.0)
     .  CALL ER('LoadTriangles','Unsupporting version',*99)

      CALL ALLOC_VERTEX(nver)
      CALL ALLOC_SURFACE(nsurface)
      CALL ALLOC_TRIANGLE(ntri)
      READ(fp,ERR=98) (tri(i1),i1=1,ntri)
      READ(fp,ERR=98) ((ver(i1,i2),i2=1,3),i1=1,nver)
      READ(fp,ERR=98) (surface(i1),i1=1,nsurface)

c      WRITE(0,*) 'LOADING TRIANGLES:',ntri,nver,nsurface

c      READ(fp,ERR=98) tri,ver,add
      CLOSE (fp)
      
      RETURN
 98   CALL ER('LoadTriangles_06','Problems reading data file',*99)
 99   STOP
      END
c
c ======================================================================
c
      RECURSIVE SUBROUTINE LoadTriangleData(flag1,flag2,flag3,normalize,
     .                                      tdata,filename)
c      USE mod_eirene04
      IMPLICIT none

      INTEGER   flag1,flag2,flag3,normalize
      REAL      tdata(*)
      CHARACTER filename*(*)

      INTEGER GetNumberOfObjects

      INTEGER   fp,ntally,ndata,icount,i1,index(20),ntri,
     .          iblk,iatm,imol,iion,ipho,ilin
      LOGICAL   output
      REAL      rdum(30),volmin
      CHARACTER buffer*256,fname*512

      REAL, ALLOCATABLE :: tvol(:)      

      output = .FALSE. 

      IF (output) WRITE(0,*) 'LOADTRIANGLEDATA:',flag1,flag2,flag3

c      tdata = 0.0  ! Initialization... problem, size unknown...

      IF (normalize.EQ.1) THEN
c...    Load volumes:
        volmin = 1.0E+20
        ntri = GetNumberOfObjects(TRIM(filename))
        ALLOCATE(tvol(ntri))
        CALL LoadTriangleData(7,0,4,0,tvol,TRIM(filename))
      ENDIF

      fname = TRIM(filename)
      IF (TRIM(filename).EQ.'default') fname = 'eirene.transfer'

      fp = 99
      OPEN(UNIT=fp,FILE=TRIM(fname),ACCESS='SEQUENTIAL',
     .     STATUS='OLD',ERR=98)

      iblk = 0
      iatm = 0
      imol = 0
      ipho = 0
      ilin = 0
      DO WHILE (.TRUE.)
        READ(fp,'(A256)',END=10) buffer

        IF     (buffer(1:22).EQ.'* BULK PARTICLES - VOL') THEN
          IF (output) WRITE(0,*) 'FOUND: ',buffer(1:12)
          iblk = iblk + 1
          READ(fp,*,ERR=97) ntally
          READ(fp,*,ERR=97) ndata                         
          READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          icount = 0
          DO WHILE (icount.LT.ndata)
            CALL NextLine(fp,ntally,icount,rdum)
            IF (flag1.EQ.1.AND.flag2.EQ.iblk) tdata(icount)=rdum(flag3)
          ENDDO
        ELSEIF (buffer(1:22).EQ.'* BULK PARTICLES - SUR') THEN
        ELSEIF (buffer(1:18).EQ.'* TEST ATOMS - VOL') THEN
          IF (output) WRITE(0,*) 'FOUND: ',buffer(1:12)
          iatm = iatm + 1
          READ(fp,*,ERR=97) ntally
          READ(fp,*,ERR=97) ndata                         
          READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          icount = 0
          DO WHILE (icount.LT.ndata)
            CALL NextLine(fp,ntally,icount,rdum)
            IF (flag1.EQ.2.AND.flag2.EQ.iatm) tdata(icount)=rdum(flag3)
          ENDDO
        ELSEIF (buffer(1:18).EQ.'* TEST ATOMS - SUR') THEN
        ELSEIF (buffer(1:22).EQ.'* TEST MOLECULES - VOL') THEN
          IF (output) WRITE(0,*) 'FOUND: ',buffer(1:12)
          imol = imol + 1
          READ(fp,*,ERR=97) ntally
          READ(fp,*,ERR=97) ndata                         
          READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          icount = 0
          DO WHILE (icount.LT.ndata)
            CALL NextLine(fp,ntally,icount,rdum)
            IF (flag1.EQ.3.AND.flag2.EQ.imol) tdata(icount)=rdum(flag3)
          ENDDO
        ELSEIF (buffer(1:17).EQ.'* TEST IONS - VOL') THEN
        ELSEIF (buffer(1:17).EQ.'* TEST IONS - SUR') THEN
        ELSEIF (buffer(1:20).EQ.'* TEST PHOTONS - VOL') THEN
          IF (output) WRITE(0,*) 'FOUND: ',buffer(1:12)
          ipho = ipho + 1
          READ(fp,*,ERR=97) ntally
          READ(fp,*,ERR=97) ndata                         
          READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          icount = 0
          DO WHILE (icount.LT.ndata)
            CALL NextLine(fp,ntally,icount,rdum)
            IF (flag1.EQ.5.AND.flag2.EQ.ipho) tdata(icount)=rdum(flag3)
          ENDDO
        ELSEIF (buffer(1:14).EQ.'* LINE EMISSIO') THEN
          IF (output) WRITE(0,*) 'FOUND: ',buffer(1:12)
          ilin = ilin + 1
          READ(fp,*,ERR=97) ntally
          READ(fp,*,ERR=97) ndata   
          READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          icount = 0
          DO WHILE (icount.LT.ndata)
            CALL NextLine(fp,ntally,icount,rdum)
            IF (flag1.EQ.6.AND.flag2.EQ.ilin) tdata(icount)=rdum(flag3)
          ENDDO

        ELSEIF (buffer(1:6 ).EQ.'* MISC') THEN
c...      Check volumes:
          IF (output) WRITE(0,*) 'FOUND: ',buffer(1:12)
          READ(fp,*,ERR=97) ntally
          READ(fp,*,ERR=97) ndata   
          READ(fp,*,ERR=97) (index(i1),i1=1,ntally-1)     ! TEMP WITH THE -1     
          icount = 0
          DO WHILE (icount.LT.ndata)
            CALL NextLine(fp,ntally,icount,rdum)
            IF (flag1.EQ.7) tdata(icount) = rdum(flag3)
            IF (flag1.EQ.7.AND.flag3.EQ.4.AND.rdum(flag3).EQ.0.0) THEN
              WRITE(0,*) 'NULL TRIANGLE VOLUME:',icount
              STOP 'HALTING CODE'
c              tdata(icount)=1.0E+10
            ENDIF
          ENDDO
        ELSEIF (buffer(1:1 ).EQ.'*') THEN
        ELSE
c          CALL ER('LoadTriangleData','Problem with transfer file',*99)
        ENDIF
      ENDDO
 10   CONTINUE

      CLOSE (fp)

      IF (normalize.EQ.1) THEN
c...    Scale by triangle volume:
        DO i1 = 1, ntri
          volmin = MIN(volmin,tvol(i1))
          tdata(i1) = tdata(i1) / tvol(i1) * 1.0E+06  ! Conversion to m-3
        ENDDO
c        WRITE(0,*) 'VOLMIN:',volmin
        IF (ALLOCATED(tvol)) DEALLOCATE(tvol)
      ENDIF

      RETURN
 97   WRITE(0,*) 'ERROR: PROBLEM READING DATA TRANSFER FILE'
      STOP
 98   WRITE(0,*) 'WARNING: eirene.transfer DATA FILE NOT FOUND'
      RETURN
 99   WRITE(0,*) 'BUFFER: >'//TRIM(buffer)//'<'
      END
c
c ======================================================================
c Geometry functions:
c ======================================================================
c
c ======================================================================
c
c subroutine: PointOnLine
c
      LOGICAL FUNCTION PointOnLine(x,y,s,t,mode,output)
c      USE mod_eirene04
      IMPLICIT none

      REAL*8     DTOL
c      PARAMETER (DTOL=1.0D-05)

      INTEGER mode
      LOGICAL test,output
      REAL*8  x(0:2),y(0:2),s,t,length

      INTEGER fp

      fp = 88

      IF     (mode.EQ.1.OR.mode.EQ.4) THEN
c...    Include end points as .TRUE.:
        DTOL = +1.0D-05
      ELSEIF (mode.EQ.2.OR.mode.EQ.3) THEN
c...    Don't include end points as .TRUE.:
        DTOL = -1.0D-05
      ELSEIF (mode.EQ.5.OR.mode.EQ.7) THEN
c...    More strict than MODE.EQ.2 or MODE.EQ.3:
        DTOL = -1.0D-06
      ELSEIF (mode.EQ.6) THEN
c...    Less strict than MODE.EQ.1:
        DTOL = +1.0D-07
      ELSE
        CALL ER('PointOnLine','Invalid MODE',*99)
      ENDIF
 
      test = .TRUE.
      
      s = -999.0D0
      t = -999.0D0

      IF (DABS(x(0)-x(1)).LT.DABS(DTOL)) THEN
c...    This always fails for MODE.EQ.2:  ! WHY WAS I EVER HAPPY WITH THIS? 
        test = test.AND.DABS(x(0)-x(2)).LT.DABS(DTOL)
c        test = test.AND.DABS(x(0)-x(2)).LT.DTOL  ! BUG
      ELSE
        s = (x(2) - x(0)) / (x(1) - x(0))
        test = test.AND.s.GT.0.0D0-DTOL.AND.s.LT.1.0D0+DTOL
      ENDIF

      IF (output) WRITE(fp,'(A,2F12.6,1P,E14.7,L2)') 
     .  'S TEST:',s,0.0D0-DTOL,1.0D0+DTOL,DTOL

      IF (test) THEN
        IF (DABS(y(0)-y(1)).LT.DABS(DTOL)) THEN
          IF (output) THEN
            WRITE(fp,*) '  T1',y(0),DABS(DTOL)
            WRITE(fp,*) '  T1',y(1),DABS(y(0)-y(1))
            WRITE(fp,*) '  T1',y(2),DABS(y(0)-y(2))
          ENDIF
          test = test.AND.DABS(y(0)-y(2)).LT.DABS(DTOL)
c          test = test.AND.DABS(y(0)-y(2)).LT.DTOL ! BUG
        ELSE
          t = (y(2) - y(0)) / (y(1) - y(0))
          test = test.AND.t.GT.0.0D0-DTOL.AND.t.LT.1.0D0+DTOL
          IF (output) WRITE(fp,*) '  T2',t,test
        ENDIF
c...    This requirement is a little more relaxed due to small cells:   ! *PERHAPS SCALE ACCORDING TO SIDE LENGTH?*
        IF ((mode.EQ.1.OR.mode.EQ.2.OR.mode.EQ.5.OR.mode.EQ.6).AND. 
c        IF (mode.NE.3.AND.mode.NE.4.AND.
     .      s.NE.-999.0D0.AND.t.NE.-999.0D0)
     .    test=test.AND.DABS(s-t).LT.100.0D0*DABS(DTOL)
        IF (mode.EQ.7.AND.
     .      s.NE.-999.0D0.AND.t.NE.-999.0D0) THEN
          length = DSQRT((x(1)-x(0))**2+(y(1)-y(0))**2)
          test=test.AND.DABS(s-t)/t.LT.MAX(0.001,0.001*0.1/length)
c          test=test.AND.DABS(s-t).LT.1000.0D0*DABS(DTOL)
        ENDIF
      ENDIF

      IF (output) THEN
        length = DSQRT((x(1)-x(0))**2+(y(1)-y(0))**2)
        WRITE(fp,'(A,L2)') 'S&T:',test
        WRITE(fp,'(A,2F12.6,1P,E14.7)') 'S&T:',s,t,DABS(s-t)
        WRITE(fp,'(A,2F12.6)') 'S&T:',DABS(s-t)/t,
     .                                MAX(0.001,0.001*0.1/length)
        WRITE(fp,*) 'LENGTH:',DSQRT((x(1)-x(0))**2+(y(1)-y(0))**2)
      ENDIF

      IF (.NOT.test) s = -1.0D0

      PointOnLine = test

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      LOGICAL FUNCTION PointInPolygon(x,y,n,p)
      IMPLICIT none

      INTEGER, INTENT(IN) :: n
      REAL*8 , INTENT(IN) :: x,y,p(n,2)

      REAL*8 , PARAMETER :: DTOL=1.0D-07

      INTEGER iseg,ninter,i,j,fp
      LOGICAL debug
      REAL*8  x1,x2,x3,x4,y1,y2,y3,y4,s12,s34

      fp = 88
      debug = .FALSE.

      PointInPolygon = .FALSE.

      x1 = x
      y1 = y
      x2 = x1 + 5.0D0
      y2 = y1
      ninter = 0
      DO i = 1, n
        j = i + 1
        IF (j.EQ.n+1) j = 1
        x3 = p(i,1)
        y3 = p(i,2)
        x4 = p(j,1)
        y4 = p(j,2)
        CALL CalcInter(x1,y1,x2,y2,x3,y3,x4,y4,s12,s34) 
        IF (s12.GT.DTOL.AND.s34.GT.0.0D0.AND.s34.LT.1.0D0) 
     .    ninter = ninter + 1
        IF (debug) WRITE(fp,'(4X,A,2F14.7,I4,2F12.5)')
     .      '           :',s12,s34,ninter,x1,y1
      ENDDO  

      IF (ninter.GT.0.AND.MOD(ninter+1,2).EQ.0) PointInPolygon = .TRUE.

      RETURN
 99   STOP
      END
c
c function: CalcTriangleArea
c
      REAL*8 FUNCTION CalcTriangleArea(x1,y1,x2,y2,x3,y3)
      IMPLICIT none

      REAL*8, INTENT(IN) :: x1,y1,x2,y2,x3,y3

      INTEGER  i1,i2,i3
      REAL*8   area,vx(4),vy(4),t1,t2

      area = 0.0D0

      IF     (y1.EQ.y3) THEN
        area = 0.5D0 * DABS(x1 - x3) * DABS(y1 - y2)
      ELSEIF (y2.EQ.y3) THEN
        area = 0.5D0 * DABS(x2 - x3) * DABS(y2 - y1)
      ELSEIF (y1.EQ.y2) THEN
        area = 0.5D0 * DABS(x1 - x2) * DABS(y1 - y3)
      ELSE
        vx(1) = x1
        vy(1) = y1
        vx(2) = x2
        vy(2) = y2
        vx(3) = x3
        vy(3) = y3
        IF     (y1.LT.y2.AND.y1.GT.y3.OR.
     .          y1.LT.y3.AND.y1.GT.y2) THEN
          i1 = 2
          i2 = 1
          i3 = 3
        ELSEIF (y2.LT.y1.AND.y2.GT.y3.OR.
     .          y2.LT.y3.AND.y2.GT.y1) THEN
          i1 = 1
          i2 = 2
          i3 = 3
        ELSEIF (y3.LT.y2.AND.y3.GT.y1.OR.
     .          y3.LT.y1.AND.y3.GT.y2) THEN
          i1 = 1
          i2 = 3
          i3 = 2
        ELSE
          CALL ER('CalcTriangleArea','Invalid vertex data',*99)
        ENDIF

        CALL CalcInter(vx(i1),vy(i1),vx(i3),vy(i3),
     .                 -100.0D0,vy(i2),+100.0D0,vy(i2),t1,t2)

        IF (t1.GT.0.0D0.AND.t1.LT.1.0D0) THEN
          vx(4) = vx(i1) + t1 * (vx(i3) - vx(i1))
          vy(4) = vy(i1) + t1 * (vy(i3) - vy(i1))
          area = 0.5 * DABS(vx(i2) - vx(4)) * DABS(vy(i1) - vy(i2)) +
     .           0.5 * DABS(vx(i2) - vx(4)) * DABS(vy(i3) - vy(i2))
        ELSE
          CALL ER('CalcTriangleArea','Intersection not found',*99)
        ENDIF
      ENDIF

      CalcTriangleArea = area

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: CalcPolygonArea
c
c a better solution would be to scan each vertex to be the common vertex,
c instead of the "center" point, for intersections with cell sides when going
c from that vertex to all of the other verticies.  If there is an intersection
c then you don't use that vertex, otherwise simple convex polygons are okay.
c
      REAL*8 FUNCTION CalcPolygonArea(x,y,n)
      IMPLICIT none

      INTEGER, INTENT(IN) :: n
      REAL*8 , INTENT(IN) :: x(n),y(n)

      REAL*8 CalcTriangleArea

      REAL*8, PARAMETER :: DTOL = 1.0E-6

      INTEGER i,j
      REAL*8  area,xcen,ycen

      area = -1.0D0

c...  Check that the polygon is continuous (no breaks):
c      DO i = 1, n
c        j = i + 1
c        IF (i.EQ.n) j = 1
c        IF (DABS(x(j)-x(i)).LT.DTOL.OR.DABS(y(j)-y(i)).LT.DTOL) 
c     .    CALL ER('CalcPolygonArea','Polygon not continuous',*99)
c      ENDDO
c...  Approximate centroid of the polygon:
      xcen = SUM(x(1:n)) / DBLE(n)
      ycen = SUM(y(1:n)) / DBLE(n)

c      WRITE(0,*) 'CEN:',xcen,ycen

c...  Calculate area by dividing the polygon into triangles and adding
c     up their respective areas (can fail if the polygon is convex, but
c     no checking done at the moment):
      area = 0
      DO i = 1, n
        j = i + 1
        IF (i.EQ.n) j = 1
        area = area + CalcTriangleArea(x(i),y(i),x(j),y(j),xcen,ycen)
      ENDDO

      CalcPolygonArea = area

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: CalcPolygonToroidalVolume
c
c
c
      REAL*8 FUNCTION CalcPolygonToroidalVolume(x,y,n)
      IMPLICIT none

      INTEGER, INTENT(IN) :: n
      REAL*8 , INTENT(IN) :: x(n),y(n)

      REAL*8 CalcPolygonArea

      REAL*8, PARAMETER :: PI = 3.14159265358979323846D0  

      REAL*8  area,xcen

      area = CalcPolygonArea(x,y,n)

c...  Approximate radial centroid of the polygon:
      xcen = SUM(x(1:n)) / DBLE(n)

      CalcPolygonToroidalVolume = area * 2.0D0 * PI * xcen

      RETURN
 99   STOP
      END
c
c ========================================================================
c
c subroutine: CalcInter
c
c
      SUBROUTINE CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
      IMPLICIT none

      REAL*8, INTENT(IN)  :: a1,a2,b1,b2,c1,c2,d1,d2
      REAL*8, INTENT(OUT) :: tab,tcd

      REAL*8, PARAMETER :: DTOL = 1.0D-07, HI = 1.0D+37, LO = 1.0D-37

      REAL*8 t0,e1,e2,f1,f2,fact,crossprod

      IF (DABS(c1-d1).LT.DTOL.AND.DABS(c2-d2).LT.DTOL) THEN
c        WRITE(88,*) '   C = D'
        tab = HI
        tcd = HI
        RETURN
      ENDIF
c
c
c
      IF (((a1 - b1) * (c2 - d2)).EQ.0.0D0.AND.
     .    ((a2 - b2) * (c1 - d1)).EQ.0.0D0) THEN
c        WRITE(88,*) '   CP is zero'
        tab = HI
        tcd = HI
        RETURN
      ENDIF

c      lenab = SQRT((b1 - a1)**2.0 + (b2 - a2)**2.0)
c      lencd = SQRT((b1 - a1)**2.0 + (b2 - a2)**2.0)
c
c      dotprod = ((b1 - a1) * (d1 - c1) + (b2 - a2) * (d2 - c2)) /
c                (lenab * lencd)

c
c     Find projection of AB onto CD:
c
      t0 = ((a1 - c1) * (d1 - c1) + (a2 - c2) * (d2 - c2)) /
     .     ((c1 - d1)**2 + (c2 - d2)**2)

      e1 = c1 + t0 * (d1 - c1)
      e2 = c2 + t0 * (d2 - c2)
c
c     Calculate F, the intersection point between AB and CD:
c
c
c
c
      fact = (e1 - a1) * (b1 - a1) + (e2 - a2) * (b2 - a2)

      IF (fact.EQ.0.0D0) THEN
c        WRITE(SLOUT,'(A)') '    fact = 0'
c        WRITE(SLOUT,'(A,1P,E15.7)') '    cp   =  ',crossprod
c        WRITE(SLOUT,'(5(2F8.4,2X))') a1,a2,b1,b2,c1,c2,d1,d2,e1,e2
        tab = 0.0
      ELSE
        tab = ((a1 - e1)**2 + (a2 - e2)**2) / fact
      ENDIF
c
c     Determine the parametric location of F on CD:
c
      f1 = a1 + tab * (b1 - a1)
      f2 = a2 + tab * (b2 - a2)

      IF (DABS(d1-c1).GT.DABS(d2-c2)) THEN
        tcd = (f1 - c1) / (d1 - c1)
      ELSE
        tcd = (f2 - c2) / (d2 - c2)
      ENDIF

      RETURN
c
c     Error output:
c
99    CONTINUE
      WRITE(0,'(5X,A,2F13.10,A,2F13.10)')
     .  ' A1,A2 = ',a1,a2,' B1,B2 = ',b1,b2
      WRITE(0,'(5X,A,2F13.10,A,2F13.10)')
     .  ' C1,C2 = ',c1,c2,' D1,D2 = ',d1,d2
      WRITE(0,'(5X,A,2F13.10,A,2F13.10)')
     .  ' E1,E2 = ',e1,e2,' F1,F2 = ',f1,f2
      WRITE(0,'(5X,A,E10.3,A,3E10.3)')
     .  ' FACT  = ',fact,' T0,TAB,TCD = ',t0,tab,tcd
      STOP
      END
c
c ======================================================================
c

c
c ======================================================================
c
c
c
      REAL FUNCTION CalcFlux(itarget,itube)
      USE mod_sol28_params
      USE mod_sol28_global
      IMPLICIT none

      REAL GetCs

      INTEGER itarget,itube,ion
      REAL    area

      IF (itarget.NE.LO.AND.itarget.NE.HI) 
     .  CALL ER('CalcFlux','Invalid target region',*99)

      IF (itube.LT.grid%isep) THEN
        CalcFlux = 0.0
      ELSE
        ion = 1

        area = V_PI * (tube(itube)%rp (itarget) +
     .                 tube(itube)%dds(itarget) * 0.5)**2 - 
     .         V_PI * (tube(itube)%rp (itarget) -
     .                 tube(itube)%dds(itarget) * 0.5)**2

c        WRITE(88,*) 'AREA:',area

        area = 2.0 * V_PI * tube(itube)%rp (itarget) * 
     .                      tube(itube)%dds(itarget)

c        WRITE(88,*) 'AREA:',area

        CalcFlux = tube(itube)%ni    (itarget,ion) * 
     .             tube(itube)%vi    (itarget,ion) * 
     .             tube(itube)%bratio(itarget) * 
     .             tube(itube)%costet(itarget) *
     .             area
c        CalcFlux = knds(id) * kvds(id) * dds2(id) * brat *
c     .             2.0 * PI * rp(id) * costet(id) * eirsrcmul * 
c     .             eirtorfrac
      ENDIF

      RETURN
99    STOP
      END
c
c ======================================================================
c
      REAL FUNCTION GetCs2(te,ti)
      IMPLICIT none

      REAL, INTENT(IN) :: te,ti
      REAL a,z

c     Te,i in eV
c     a    in amu
c     result is m s-1

      z = 1.0
      a = 2.0

      GetCs2 = SQRT((te + ti) * 1.60E-19 / (a * 1.67E-27))  ! Needs improvement... 

c      GetCs2 = 9.78817E+03 * SQRT(0.5 * (1.0 + a) * (te + ti) / z)

      RETURN
99    STOP
      END
c
c ======================================================================
c
      REAL FUNCTION GetJsat2(te,ti,ne,v)
      IMPLICIT none

      REAL, INTENT(IN) :: te,ti,ne,v
      REAL                vb

      REAL GetCs2

      IF (v.EQ.1.0) THEN
        vb = GetCs2(te,ti)
      ELSE
        vb = v
      ENDIF

      GetJsat2 = ne * 1.6022E-19 * vb

      RETURN
99    STOP
      END

c ======================================================================
c
c subroutine: ProcessIterationBlocks
c
c Moved here for compatilibity with OUT.
c 
      SUBROUTINE ProcessIterationBlocks
      USE mod_sol28_params
      USE mod_sol28_global
      USE mod_legacy
      IMPLICIT none

      LOGICAL osmGetLine
      INTEGER, PARAMETER :: WITH_TAG = 1, NO_TAG = 2

      INTEGER   fp,i,idum1(1:5)
      CHARACTER buffer*1024,cdum1*512

      LOGICAL :: status = .TRUE., new_block = .FALSE., load_data

      opt%tube(1)      = 1
      opt%tube(2)      = 1E+8
      opt%iteration(1) = 1
      opt%iteration(2) = 1E+8
      nopt = 1
      opt_iteration(1) = opt
      IF (niteration.GT.0) THEN
c...    Load data to selectively modify input options for particular 
c       flux-tubes and iterations of the solver:
        fp = -1
        niteration = niteration + 1
        iteration_buffer(niteration) = '''{EXIT}'''
        DO WHILE(osmGetLine(fp,buffer,WITH_TAG))
c          WRITE(0,*) 'buffer >'//TRIM(buffer)//'<'                        
c...      Isolate tag string:
          DO i = 2, LEN_TRIM(buffer)
            IF (buffer(i:i).EQ.'}') EXIT
          ENDDO
c          WRITE(0,*) 'buffer >'//TRIM(buffer(2:i))//'<'                        
          SELECTCASE (buffer(3:i-1))
            CASE ('CON ITERATION DATA')
              status = .TRUE.
              opt_iteration(nopt) = opt
              nopt = nopt + 1
              READ(buffer,*) cdum1,idum1(1:5)
              SELECTCASE (idum1(1))
                CASE (0) ! Not active
                  nopt = nopt - 1              
                CASE (1) ! Iteration block always based on the initial block (master block)
                  opt = opt_iteration(1)
                CASE (2) ! Iteration block based on previous iteration block
                  opt = opt_iteration(nopt-1)
                CASE DEFAULT
                  CALL ER('ProcessIterationBlocks','Unknown option',*99)
              ENDSELECT
              IF (idum1(1).NE.0) THEN
                load_data = .TRUE.
                opt%iteration(1:2) = idum1(2:3)
                opt%tube     (1:2) = idum1(4:5)
              ELSE
                load_data = .FALSE.
              ENDIF
            CASE ('EXIT')
              status = .FALSE.
            CASEDEFAULT
              IF (load_data) CALL ProcessInputTag(fp,i,buffer,status)
          ENDSELECT
          IF (.NOT.status) EXIT
        ENDDO
        opt_iteration(nopt) = opt
        opt = opt_iteration(1)
      ENDIF

      WRITE(logfp,*)
      WRITE(logfp,*) 'NITERATION = ',niteration
      WRITE(logfp,*) 'NOPT       = ',nopt
      DO i = 1, nopt
        WRITE(logfp,*) '  RANGE ITERATIION,TUBE  = ',
     .    opt_iteration(i)%iteration,opt_iteration(i)%tube,
     .    opt_iteration(i)%p_ion
      ENDDO

      RETURN
 99   WRITE(0,*) ' IDUM1(1) = ',idum1(1)
      STOP
      END
c
c ======================================================================
c
      SUBROUTINE GetVertex(iobj,i,x,y)
      USE mod_geometry
      IMPLICIT none

      INTEGER, INTENT(IN ) :: iobj,i
      REAL*8 , INTENT(OUT) :: x,y

      INTEGER isrf,ivtx
 
c...  The first vertex of each side corresponds to i:      
      
      isrf = obj(iobj)%iside(i)

      IF (srf(ABS(isrf))%type.NE.SPR_LINE_SEGMENT) 
     .  CALL ER('GetVertex','Routine only handles line segments',*99)

      IF (isrf.LT.0) THEN
        isrf = -isrf
        ivtx = srf(isrf)%ivtx(srf(isrf)%nvtx)
        x = vtx(1,ivtx)
        y = vtx(2,ivtx)
      ELSE
        ivtx = srf(isrf)%ivtx(1)
        x = vtx(1,ivtx)
        y = vtx(2,ivtx)
      ENDIF

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c
c
      REAL FUNCTION osm_GetFlux(itarget,itube)
      IMPLICIT none

      REAL CalcFlux

      INTEGER, INTENT(IN) :: itarget,itube
      REAL    brat

      osm_GetFlux = CalcFlux(itarget,itube)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c
c
      REAL FUNCTION osm_GetGamma(itarget,itube)
      USE mod_sol28_params
      USE mod_sol28_global
      IMPLICIT none

      REAL, PARAMETER :: ECH = 1.6022E-19, AMU = 1.67E-27, PI = 3.141592

      INTEGER, INTENT(IN) :: itarget,itube

      INTEGER ion
      REAL    mi,t_ratio,delta_e,m_ratio,log_arguement

      ion = 1

      IF     (itarget.EQ.LO) THEN
      ELSEIF (itarget.EQ.HI) THEN
      ELSE
        CALL ER('osm_GetGamma','Invalid target specification',*99)
      ENDIF

      IF (itube.LT.grid%isep) THEN
        osm_GetGamma = 0.0
      ELSE
c...    Taken from Stangeby 1st edition, equation 25.46, pg 649.  Note that many effects
c       are missing, i.e. realistic secondary electron emission (0.0 here), e-i recombination
c       energy, atom-atom recombination energy, low collisionality effects, space charge
c       effects, etc. see the discussion by Stangeby pp 646-654. -SL, 29.03.2010
        mi = 2.0   ! *** MASS HARDCODED! ***
        t_ratio = tube(itube)%ti(itarget,ion) / tube(itube)%te(itarget)
        delta_e = 0.0
        m_ratio = 9.11E-31 / (mi * AMU)
        log_arguement = 2.0 * PI * m_ratio * (1.0 + t_ratio) * 
     .                 (1.0 - delta_e)**-2 
        osm_GetGamma = 2.5 * t_ratio + 2.0 / (1.0 - delta_e) - 
     .                 0.5 * LOG( log_arguement )
      ENDIF

      RETURN
99    STOP
      END
c
c ======================================================================
c
c
c
      REAL FUNCTION osm_GetHeatFlux(itarget,itube)
      USE mod_sol28_params
      USE mod_sol28_global
      IMPLICIT none

      REAL osm_GetFlux,osm_GetGamma

      REAL, PARAMETER :: ECH = 1.6022E-19

      INTEGER, INTENT(IN) :: itarget,itube
      INTEGER id,ir
      REAL    isat,gamma

      IF     (itarget.EQ.LO) THEN
      ELSEIF (itarget.EQ.HI) THEN
      ELSE
        CALL ER('osm_GetHeatFlux','Invalid target specification',*99)
      ENDIF

      IF (itube.LT.grid%isep) THEN
        osm_GetHeatFlux = 0.0
      ELSE
        isat  = osm_GetFlux (itarget,itube)
        gamma = osm_GetGamma(itarget,itube)
        osm_GetHeatFlux = gamma * ABS(isat) * ECH * 
     .                    tube(itube)%te(itarget)
      ENDIF

      RETURN
99    STOP
      END
