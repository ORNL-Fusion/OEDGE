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
      USE mod_sol28_global
      IMPLICIT none

      INTEGER i1,fp,ion
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

c      WRITE(0,*) 'NTUBE:',ntube
c      STOP 'sdgfsdsdg'

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
      ref_nfluid = 0
      IF (ALLOCATED(ref_tube )) DEALLOCATE(ref_tube) 
      IF (ALLOCATED(ref_fluid)) DEALLOCATE(ref_fluid) 

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

      INTEGER inode,ion
      REAL*8  pressure,p(2)

      pressure = 0.0D0

      IF     (inode.EQ.1) THEN
        pressure = pe(ictarg(LO)) + pi(ictarg(LO),1)
      ELSEIF (inode.EQ.nnode) THEN
        pressure = pe(ictarg(HI)) + pi(ictarg(HI),1)
      ELSEIF (node(inode)%ne.NE.0.0.AND.
     .        node(inode)%pe.NE.0.0) THEN
        STOP 'NOT SURE WHAT TO DO'
      ELSEIF (node(inode)%ne.GT.0.0) THEN
        SELECTCASE (0)
          CASE (0)
            pressure = ECH * DBLE(node(inode)%ne *              ! Inadequate...v.NE.0
     .                 (node(inode)%te + node(inode)%ti(ion)))  ! Account for kinetic ions 
          CASEDEFAULT                                           ! NODE%TE may not be set... ion.NE.1            
        ENDSELECT                                                         
      ELSEIF (node(inode)%pe.GT.0.0) THEN
        SELECTCASE (0)
          CASE (0)
            pressure = 2.0D0 * DBLE(node(inode)%pe) * ECH     ! See above
          CASEDEFAULT 
        ENDSELECT
      ENDIF

      GetNodePressure = pressure

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      INTEGER FUNCTION GetNumberOfObjects()
      IMPLICIT none

      INTEGER   fp,idum1
      CHARACTER buffer*256

      fp = 99
      OPEN(UNIT=fp,FILE='eirene.transfer',ACCESS='SEQUENTIAL',
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
      RECURSIVE SUBROUTINE LoadTriangleData(flag1,flag2,flag3,normalize,
     .                                      tdata)
c      USE mod_eirene04
      IMPLICIT none

      INTEGER flag1,flag2,flag3,normalize
      REAL    tdata(*)

      INTEGER GetNumberOfObjects

      INTEGER   fp,ntally,ndata,icount,i1,index(20),ntri,
     .          iblk,iatm,imol,iion,ipho,ilin
      REAL      rdum(30),volmin
      CHARACTER buffer*256

      REAL, ALLOCATABLE :: tvol(:)      

      WRITE(0,*) 'LOADTRIANGLEDATA:',flag1,flag2,flag3

c      tdata = 0.0  ! Initialization... problem, size unknown...

      IF (normalize.EQ.1) THEN
c...    Load volumes:
        volmin = 1.0E+20
        ntri = GetNumberOfObjects()
        ALLOCATE(tvol(ntri))
        CALL LoadTriangleData(7,0,4,0,tvol)
      ENDIF

      fp = 99
      OPEN(UNIT=fp,FILE='eirene.transfer',ACCESS='SEQUENTIAL',
     .     STATUS='OLD',ERR=98)

      iblk = 0
      iatm = 0
      imol = 0
      ipho = 0
      ilin = 0
      DO WHILE (.TRUE.)
        READ(fp,'(A256)',END=10) buffer

        IF     (buffer(1:16).EQ.'* BULK PARTICLES') THEN

          iblk = iblk + 1
          READ(fp,*,ERR=97) ntally
          READ(fp,*,ERR=97) ndata                         
          READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          icount = 0
          DO WHILE (icount.LT.ndata)
            CALL NextLine(fp,ntally,icount,rdum)
            IF (flag1.EQ.1.AND.flag2.EQ.iblk) tdata(icount)=rdum(flag3)
          ENDDO

        ELSEIF (buffer(1:12).EQ.'* TEST ATOMS') THEN

          iatm = iatm + 1
          READ(fp,*,ERR=97) ntally
          READ(fp,*,ERR=97) ndata                         
          READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          icount = 0
          DO WHILE (icount.LT.ndata)
            CALL NextLine(fp,ntally,icount,rdum)
            IF (flag1.EQ.2.AND.flag2.EQ.iatm) tdata(icount)=rdum(flag3)
          ENDDO

        ELSEIF (buffer(1:16).EQ.'* TEST MOLECULES') THEN

          imol = imol + 1
          READ(fp,*,ERR=97) ntally
          READ(fp,*,ERR=97) ndata                         
          READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          icount = 0
          DO WHILE (icount.LT.ndata)
            CALL NextLine(fp,ntally,icount,rdum)
            IF (flag1.EQ.3.AND.flag2.EQ.imol) tdata(icount)=rdum(flag3)
          ENDDO

        ELSEIF (buffer(1:11).EQ.'* TEST IONS') THEN
        ELSEIF (buffer(1:14).EQ.'* TEST PHOTONS') THEN

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
c TEMP          CALL ER('LoadTriangleData','Problem with transfer file',*99)
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

      IF (target.NE.LO.AND.target.NE.HI) 
     .  CALL ER('CalcFlux','Invalid target region',*99)

      IF (itube.LT.grid%isep) THEN
        CalcFlux = 0.0
      ELSE
        ion = 1
        CalcFlux = tube(itube)%ni    (target,ion) * 
     .             tube(itube)%vi    (target,ion) * 
     .             tube(itube)%dds   (target) * 
     .             tube(itube)%bratio(target) * 
     .             2.0 * V_PI * 
     .             tube(itube)%rp    (target) *
     .             tube(itube)%costet(target)
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

      GetCs2 = 9.78817E+03 * SQRT(0.5 * (1.0 + a) * (te + ti) / z)

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

      GetJsat2 = ne * 1.602E-19 * vb

      RETURN
99    STOP
      END





