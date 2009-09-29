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
 98   CALL ER('LoadGrid','Problem loading OSM grid file',*99)
 99   WRITE(0,*) '  FILE NAME='//'"'//fname(1:LEN_TRIM(fname))//'"'
      STOP
      END
c
c ======================================================================
c
      SUBROUTINE CleanUp    ! Move to mod_osm.f
      USE mod_sol28_global
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
      REAL*8 FUNCTION GetNodePressure(inode,ion)
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none

      INTEGER, INTENT(IN) :: inode,ion

      INTEGER target
      REAL    te_node,ti_node
      REAL*8  pressure,p(2)

      pressure = 0.0D0

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

        IF (node(inode)%ne.GT.0.0) THEN
          SELECTCASE (0)
            CASE (0)
              pressure = ECH * DBLE(node(inode)%ne *              ! Inadequate...v.NE.0
     .                   (te_node + ti_node))  ! Account for kinetic ions 
c     .                   (node(inode)%te + node(inode)%ti(ion)))  ! Account for kinetic ions 
            CASEDEFAULT                                           ! NODE%TE may not be set... ion.NE.1            
          ENDSELECT                                                         
        ELSEIF (node(inode)%pe.GT.0.0) THEN
          SELECTCASE (0)
            CASE (0)
              pressure = 2.0D0 * DBLE(node(inode)%pe) * ECH     ! See above
            CASEDEFAULT 
          ENDSELECT
        ELSE
          CALL ER('GetNodePressure','Unknown situation',*99)
        ENDIF
      ENDIF

      GetNodePressure = pressure

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
        WRITE(0,*) 'VOLMIN:',volmin
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
c ========================================================================
c
c subroutine: CalcInter
c
c
      SUBROUTINE CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
      IMPLICIT none

c      INCLUDE 'params'
c      INCLUDE 'slcom'

c     Input:
      DOUBLE PRECISION a1,a2,b1,b2,c1,c2,d1,d2

c     Output:
      DOUBLE PRECISION tab,tcd

      DOUBLE PRECISION TOL
      PARAMETER (TOL = 1.0E-07)

      DOUBLE PRECISION t0,e1,e2,f1,f2,fact,crossprod

      REAL, PARAMETER :: HI=1.0E37, LO=1.0E-37
c
c
c
      IF (ABS(c1-d1).LT.TOL.AND.ABS(c2-d2).LT.TOL) THEN
        tab = HI
        tcd = HI
        RETURN
c        CALL ER('CalcInter','Invalid line segment',*99)
      ENDIF
c
c
c
      IF (((a1 - b1) * (c2 - d2)).EQ.0.0.AND.
     .    ((a2 - b2) * (c1 - d1)).EQ.0.0) THEN
c        WRITE(50,*) '   CP is zero'
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

      IF (fact.EQ.0.0) THEN
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

      IF (abs(d1-c1).GT.abs(d2-c2)) THEN
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
      REAL FUNCTION CalcFlux(target,itube)
      USE mod_sol28_params
      USE mod_sol28_global
      IMPLICIT none

      REAL GetCs

      INTEGER target,itube,ion
      REAL    area

      IF (target.NE.LO.AND.target.NE.HI) 
     .  CALL ER('CalcFlux','Invalid target region',*99)

      IF (itube.LT.grid%isep) THEN
        CalcFlux = 0.0
      ELSE
        ion = 1

        area = V_PI * (tube(itube)%rp (target) +
     .                 tube(itube)%dds(target) * 0.5)**2 - 
     .         V_PI * (tube(itube)%rp (target) -
     .                 tube(itube)%dds(target) * 0.5)**2

        WRITE(88,*) 'AREA:',area

        area = 2.0 * V_PI * tube(itube)%rp (target) * 
     .                      tube(itube)%dds(target)

        WRITE(88,*) 'AREA:',area

        CalcFlux = tube(itube)%ni    (target,ion) * 
     .             tube(itube)%vi    (target,ion) * 
     .             tube(itube)%bratio(target) * 
     .             tube(itube)%costet(target) *
     .             area
c        CalcFlux = knds(id) * kvds(id) * dds2(id) * brat *
c     .             2.0 * PI * rp(id) * costet(id) * eirsrcmul * 
c     .             eirtorfrac
      ENDIF




c      IF (supflx(region,ir).EQ.1) CalcFlux = CalcFlux * 1.0E-15


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

      LOGICAL GetLine
      INTEGER, PARAMETER :: WITH_TAG = 1, NO_TAG = 2

      INTEGER   fp,i,idum1(1:5)
      CHARACTER buffer*1024,cdum1*512

      LOGICAL :: status = .TRUE., new_block = .FALSE.

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
        iteration_buffer(niteration) = '''{999}'''
        DO WHILE(GetLine(fp,buffer,WITH_TAG))
c          WRITE(0,*) 'buffer >'//TRIM(buffer)//'<'                        
c...      Isolate tag string:
          DO i = 2, LEN_TRIM(buffer)
            IF (buffer(i:i).EQ.'}') EXIT
          ENDDO
c          WRITE(0,*) 'buffer >'//TRIM(buffer(2:i))//'<'                        
          IF (buffer(3:i-1).EQ.'CON ITERATION DATA') THEN
            status = .TRUE.
            opt_iteration(nopt) = opt
            nopt = nopt + 1
            READ(buffer,*) cdum1,idum1(1:5)
            IF (idum1(1).EQ.1) THEN
              opt = opt_iteration(nopt-1)
            ELSE
              opt = opt_iteration(1)
            ENDIF
            opt%iteration(1:2) = idum1(2:3)
            opt%tube     (1:2) = idum1(4:5)
          ELSE
            CALL ProcessInputTag(fp,i,buffer,status)
          ENDIF
          IF (.NOT.status) EXIT
        ENDDO
        opt_iteration(nopt) = opt
        opt = opt_iteration(1)
      ENDIF

      RETURN
 99   STOP
      END
c
c ======================================================================
c




