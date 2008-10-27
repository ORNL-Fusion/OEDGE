c     -*-Fortran-*-
c
c ======================================================================
c
c DIVIMP/EIRENE interface:
c
c
c WrtEIRENE
c b2wrpl
c readeir
c
c
c ======================================================================
c
c  subroutine: WrtEIRENE
c
      subroutine wrteirene
      IMPLICIT none

      include 'params'
      include 'cgeom'
      include 'comtor'
c slmod begin
      INCLUDE 'slcom'

      INTEGER ik,ir,in1,in2,i1,id,ik1

c slmod end

c     the purpose of this routine is to bring first the arrays
c
c     bratio(ik,ir) : magnetic field ratio B_pol/B
c     knbs  (ik,ir) : density of background plasma
c     kvhs  (ik,ir) : parallel flow velocity background plasma
c     ktibs (ik,ir) : ion temperature of background plasma
c     ktebs (ik,ir) : electron temperature
c     rs    (ik,ir) : radius of bin center
c
c     to B2 standard form, i.e. to the DIVIMP arrays we add virtual
c     grid points at the end of the rings.
c     In the following call of b2wrpl, maptodiv converts the arrays
c     to the B2 ordering scheme; gfsub3w writes the arrays conforming
c     to the B2 file structure.
c
c
c     ALSO - the velocities and fluxes need to be specified as
c            CELL EDGE values - this means that extra values
c            for the ends - at the targets - must be passed to
c            the mapping routines.
c
c     the following variables are needed in addition:
c
c     from comtor:
c     cutring  = number of separatrix ring
c     cutpt1   = number of trap knots in first half
c     cutpt2   = number of trap knots in second half of trap rings
c     maxrings = total number of rings
c     maxkpts  = total number of knots on ring
c
c     from tauin1 -> raug:
c
c     nrs     = number of rings in divimp-array
c     nks(ir) = number of knots on ring ir
c     irsep   = innermost sol-ring (separatrix)
c
c     from tauin1 -> raug -> b2repl -> (gfsub3r & maptodiv):
c
c     bratio, knbs, kvhs, ktibs, ktebs
c slmod begin
c
c Expand the grid here because not all the variables included in the
c expansion are included in Karl's cell additions (below) so things
c could get messed up if I don't go first.
c

      IF     (eiropacity.GT.0) THEN
        CALL EstimateOpacityMultiplier
      ELSEIF (eiropacity.EQ.-4) THEN
        mulrec = 1.0
        mulion = 1.0
        CALL ReadTriangles
      ELSEIF (eiropacity.EQ.-5) THEN
c...    No opacity (but opaque rates used in SOL24/28(old)):
        mulrec = 1.0
        mulion = 1.0
      ELSE
        mulrec = 1.0
        mulion = 1.0
      ENDIF


      IF (grdnmod.NE.0) THEN

        CALL OutputData(85,'Before modifying grid in EIREDIV')


        CALL StructureGrid(COMPLETE)



        cutpt1   = ikto + 1

        cutpt2 = cutpt1 + nks(2) 
c        cutpt2   = ikti + 1
        cutring  = irsep - 1
        maxkpts  = nks(irsep) + 2
        maxrings = irwall

        CALL OutputGrid(86,'Immediately after structuring grid')
      ENDIF



      IF (grdnmod.EQ.0.AND.
     .    ((cgridopt.EQ.3.AND.cmodopt.EQ.1).OR.
     .     stagopt.GT.0.OR.stopopt.EQ.121.OR.eirgrid.EQ.1)) THEN

        CALL DB('Structuring grid')

        CALL BuildMap
        CALL OutputGrid(85,'Before modifying grid in EIREDIV')

c        STOP 'BEFORE MODIFYING GRID'

        CALL StructureGrid_Old(COMPLETE)

        cutpt1   = ikto + 1
        cutpt2   = ikti + 1
        cutring  = irsep - 1
        maxkpts  = nks(irsep) + 2
        maxrings = irwall

        CALL OutputGrid(86,'Immediately after structuring grid')

      ELSE
        CALL OutputGrid(86,'Before calling EIRENE')
      ENDIF

      IF (.FALSE..AND.grd_refine.EQ.11) THEN
        cutpt1   = ikto + 1
        cutpt2   = ikti + 1
        cutring  = irsep - 1
        maxkpts  = nks(irsep) + 2
        maxrings = irwall
      ENDIF

      IF (grdnmod.EQ.0.AND.cgridopt.EQ.0) THEN
c     Make sure these are not used in any routines other than in eirediv...

c...    CODE TAGGED FOR DELETION - OR NOT
c        STOP 'CODE TAGGED FOR DELETION 001'

        CALL DB('Setting cut points')

        IF (stopopt.EQ.104) THEN
          CALL OutputGrid(85,'Before messing with rings')

          CALL InsertRing(2       ,BEFORE,0)
          CALL InsertRing(irtrap+1,BEFORE,0)

          idring(2       ) = -2
          idring(irtrap+1) = -2

          vpolmin = (MAXNKS*MAXNRS - npolyp) / 2 + npolyp
          vpolyp  = vpolmin

          ir = 2

          DO ik = 1, nks(ir)
            vpolyp = vpolyp + 1
            in1    = vpolyp
            in2    = korpg(ik,ir+1)

            korpg(ik,ir) = in1

            rvertp(1,in1) = rvertp(1,in2)
            zvertp(1,in1) = zvertp(1,in2)
            rvertp(2,in1) = rvertp(1,in2)
            zvertp(2,in1) = zvertp(1,in2)

            rvertp(3,in1) = rvertp(4,in2)
            zvertp(3,in1) = zvertp(4,in2)
            rvertp(4,in1) = rvertp(4,in2)
            zvertp(4,in1) = zvertp(4,in2)

            rs(ik,ir) = 0.25 * (rvertp(1,in1) + rvertp(2,in1) +
     .                          rvertp(3,in1) + rvertp(4,in1))
            zs(ik,ir) = 0.25 * (zvertp(1,in1) + zvertp(2,in1) +
     .                          zvertp(3,in1) + zvertp(4,in1))
          ENDDO

          ir = irtrap + 1

          DO ik = 1, nks(ir)
            vpolyp = vpolyp + 1
            in1    = vpolyp
            in2    = korpg(ik,ir+1)

            korpg(ik,ir) = in1

            rvertp(1,in1) = rvertp(1,in2)
            zvertp(1,in1) = zvertp(1,in2)
            rvertp(2,in1) = rvertp(1,in2)
            zvertp(2,in1) = zvertp(1,in2)

            rvertp(3,in1) = rvertp(4,in2)
            zvertp(3,in1) = zvertp(4,in2)
            rvertp(4,in1) = rvertp(4,in2)
            zvertp(4,in1) = zvertp(4,in2)

            rs(ik,ir) = 0.25 * (rvertp(1,in1) + rvertp(2,in1) +
     .                          rvertp(3,in1) + rvertp(4,in1))
            zs(ik,ir) = 0.25 * (zvertp(1,in1) + zvertp(2,in1) +
     .                          zvertp(3,in1) + zvertp(4,in1))
          ENDDO

          CALL OutputGrid(86,'Rings in ')


c          STOP 'Ring tricks'
        ENDIF



        cutpt1   = ikto + 1
        cutpt2   = ikti + 1
        cutring  = irsep - 1
        maxkpts  = nks(irsep) + 2
        maxrings = irwall

      ENDIF



      WRITE(0,*) '??BALANCE??: ',irsep-2,nrs-irtrap


      IF ((irsep-2.NE.nrs-irtrap).and.(nrs-irtrap).gt.0) THEN
        WRITE(0,*)
        WRITE(0,*) '*** BALANCING GRID ***'
        WRITE(0,*) irsep-2,nrs-irtrap
        WRITE(0,*)

        CALL OutputGrid(86,'Immediately before balancing grid')

        CALL BalanceGrid
      ENDIF

c      WRITE(0,*) '??BALANCE??: ',irsep-2,nrs-irtrap

c      CALL OutputData(85,'dfsfs')
c      STOP 'fsdkj'

      CALL RZero(eirnlimi,5000)
      CALL WriteBlock03a(99,99,0)
      CALL WriteBlock03b(99,99,0)


      CALL WriteInputFile


c      STOP 'AFTER WRITING INPUT FILE'

      CALL WriteGeometryFile
      IF (eirbgk.GE.2) CALL WriteVacuumGrid

c      STOP 'JUST STRUCTURED GRID'

c...MARK
c      CALL OutputData(85,'dfsfs')
c      STOP 'fsdkj'


      IF (stopopt.EQ.8) STOP 'Just structured grid'
c slmod end
      do 392 ir = irsep,nrs
         do 394 ik = nks(ir),1,-1
            bratio(ik+1,ir) = bratio(ik,ir)
            knbs  (ik+1,ir) = knbs  (ik,ir)
            knes  (ik+1,ir) = knes  (ik,ir)
            kvhs  (ik+1,ir) = kvhs  (ik,ir)
            ktibs (ik+1,ir) = ktibs (ik,ir)
            ktebs (ik+1,ir) = ktebs (ik,ir)
            rs    (ik+1,ir) = rs    (ik,ir)
c slmod begin
            virtag(ik+1,ir) = virtag(ik,ir)
            mulrec(ik+1,ir) = mulrec(ik,ir)
            mulion(ik+1,ir) = mulion(ik,ir)
            mulqer(ik+1,ir) = mulqer(ik,ir)
            mulqei(ik+1,ir) = mulqei(ik,ir)
            DO i1 = 1, MAXBGK*eirnsdtor
              pinbgk(ik+1,ir,i1) = pinbgk(ik,ir,i1)
            ENDDO
c slmod end
c
c           Also adjust kss and ksb since they are needed
c           to generate the cell boundary quantities correctly.
c
            kss (ik+1,ir)  = kss(ik,ir)
            ksb (ik+1,ir)  = ksb(ik,ir)
c
  394    continue
c
c        Setting values in the virtual points
c
c        The velocity is treated differently because it will
c        be converted to a CELL EDGE quantity by averaging the
c        contents of adjacent cells - as a result - it is
c        necessary to load the first and last cells with
c        values that will be used for the actual target velocity.
c
c slmod begin
c...
         IF (osmbulkn.NE.0.0) THEN
           id = idds(ir,2)
           bratio(1,ir) = bratio(2,ir)
           knbs  (1,ir) = knds  (id)
           knes  (1,ir) = knes  (2,ir)
           kvhs  (1,ir) = kvds  (id)
           ktibs (1,ir) = ktids (id)
           ktebs (1,ir) = kteds (id)
           rs    (1,ir) = rs    (2,ir)
           id = idds(ir,1)
           bratio(nks(ir)+2,ir) = bratio(nks(ir)+1,ir)
           knbs  (nks(ir)+2,ir) = knds  (id) 
           knes  (nks(ir)+2,ir) = knes  (nks(ir)+1,ir)
           kvhs  (nks(ir)+2,ir) = kvds  (id) 
           ktibs (nks(ir)+2,ir) = ktids (id)
           ktebs (nks(ir)+2,ir) = kteds (id)
           rs    (nks(ir)+2,ir) = rs    (nks(ir)+1,ir)
         ELSE
           bratio(1,ir) = bratio(2,ir)
           knbs  (1,ir) = knbs  (2,ir)
           knes  (1,ir) = knes  (2,ir)
           kvhs  (1,ir) = kvds  (idds(ir,2))
           ktibs (1,ir) = ktibs (2,ir)
           ktebs (1,ir) = ktebs (2,ir)
           rs    (1,ir) = rs    (2,ir)
           bratio(nks(ir)+2,ir) = bratio(nks(ir)+1,ir)
           knbs  (nks(ir)+2,ir) = knbs  (nks(ir)+1,ir)
           knes  (nks(ir)+2,ir) = knes  (nks(ir)+1,ir)
           kvhs  (nks(ir)+2,ir) = kvds  (idds(ir,1))
           ktibs (nks(ir)+2,ir) = ktibs (nks(ir)+1,ir)
           ktebs (nks(ir)+2,ir) = ktebs (nks(ir)+1,ir)
           rs    (nks(ir)+2,ir) = rs    (nks(ir)+1,ir)
         ENDIF

         virtag(1,ir) = 1 
         virtag(nks(ir)+2,ir) = 1 
c
c         bratio(1,ir) = bratio(2,ir)
c         knbs  (1,ir) = knbs(2,ir)
cc         kvhs  (1,ir) = kvds(idds(ir,2))*qtim*csolvf
c         kvhs  (1,ir) = kvds(idds(ir,2))
c         ktibs (1,ir) = ktibs(2,ir)
c         ktebs (1,ir) = ktebs(2,ir)
c         rs    (1,ir) = rs(2,ir)
cc
c         bratio(nks(ir)+2,ir) = bratio(nks(ir)+1,ir)
c         knbs  (nks(ir)+2,ir) = knbs  (nks(ir)+1,ir)
cc         kvhs  (nks(ir)+2,ir) = kvds(idds(ir,1))*qtim*csolvf
c         kvhs  (nks(ir)+2,ir) = kvds(idds(ir,1))
c         ktibs (nks(ir)+2,ir) = ktibs (nks(ir)+1,ir)
c         ktebs (nks(ir)+2,ir) = ktebs (nks(ir)+1,ir)
c         rs    (nks(ir)+2,ir) = rs (nks(ir)+1,ir)
c slmod end
c
c        Put values in kss and ksb - even if they will not be used.
c
         kss(1,ir)         = kss(2,ir)
         kss(nks(ir)+2,ir) = kss(nks(ir)+1,ir)
c
c        ksb is a 0-index array
c
         ksb(1,ir)         = ksb(0,ir)
         ksb(nks(ir)+2,ir) = ksb(nks(ir)+1,ir)
c
         nks(ir)=nks(ir)+2           ! because we have 2 more points!
c
  392 continue


c...MARK
c      CALL OutputData(85,'dfsfs')
c      STOP 'fsdkj'

c slmod begin
c      IF (cgridst.EQ.1) CALL OutputGrid(86,'junk')
c slmod end
*     added northopt to argument list; Krieger, IPP 6/95
      call b2wrpl(maxrings,maxkpts,cutring,cutpt1,cutpt2,
     >            qtim,csolvf,rizb,crmb,northopt,cioptf,cioptg,cmachno)

      do 492 ir = irsep,nrs
         do 494 ik = 1,nks(ir)-1
            bratio(ik,ir) = bratio(ik+1,ir)
            knbs  (ik,ir) = knbs  (ik+1,ir)
            knes  (ik,ir) = knes  (ik+1,ir)
            kvhs  (ik,ir) = kvhs  (ik+1,ir)
            ktibs (ik,ir) = ktibs (ik+1,ir)
            ktebs (ik,ir) = ktebs (ik+1,ir)
            rs    (ik,ir) = rs    (ik+1,ir)
c slmod begin
            virtag(ik,ir) = virtag(ik+1,ir)
            mulrec(ik,ir) = mulrec(ik+1,ir)
            mulion(ik,ir) = mulion(ik+1,ir)
            mulqer(ik,ir) = mulqer(ik+1,ir)
            mulqei(ik,ir) = mulqei(ik+1,ir)
            DO i1 = 1, MAXBGK*eirnsdtor
              pinbgk(ik,ir,i1) = pinbgk(ik+1,ir,i1)
            ENDDO
c slmod end
c
c           Readjust kss and ksb - ksb(0,...) doesn't need fixing
c                                  because it wasn't changed.
c
            kss (ik,ir) = kss(ik+1,ir)
            ksb (ik,ir) = ksb(ik+1,ir)
c
  494    continue
         nks(ir)=nks(ir)-2           ! because we deleted 2 points
  492 continue

      return
      end
c
c ======================================================================
c
      subroutine b2wrpl(mrings,mkpts,cutring,cutpt1,cutpt2,
     >                    qtim,csolvf,rizb,crmb,northopt,
     >                    cioptf,cioptg,cmachno)
      implicit none
c
      include 'params'
      include 'cgeom'
c slmod begin
      INCLUDE 'slcom'


      INTEGER cioptg,ik1
      REAL    temax
c slmod end
c
      integer mrings,mkpts,cutring,cutpt1,cutpt2,northopt
      integer cioptf
      real qtim,csolvf,rizb,crmb
      real cmachno(maxnrs,2)
c
c     B2WRPL: (B2 WRite PLasma)
c
c     this subroutine is a "mirror" of b2repl to write the arrays
c     from wrteirene to a B2 conformant background plasma file
c
c     this can of course be adapted to other B2 conformant files
c     with DIVIMP array quantities
c
c
      integer maxix,maxiy,maxis
      integer i,j,k
      integer ix,iy,is,id,ir,ik,nplasf,nx,ny,nxd,nyd
c     parameter nfla must be initialized before initialization of maxis
c     changed by Krieger IPP/97
c slmod begin
      INTEGER GetModel,SymmetryPoint

      INTEGER i1,i2,i3,nradd,in,ind1
      REAL    t1,t2,t3,tn
c slbug
      parameter (maxix=MAXNKS,maxiy=100,maxis=MAXNFLA)
c      parameter (maxix=110,maxiy=50,maxis=maxnfla)


c
c      integer nfla
c      parameter (nfla=2)
c      parameter (maxix=100,maxiy=50,maxis=nfla)
c slmod end

      real cs
      real ndummy(0:maxix+1,0:maxiy+1,maxis)
      real tdummy(0:maxix+1,0:maxiy+1)
c slmod begin
      REAL vdummy(MAXNKS,MAXNRS)
c slmod end
      real uu(maxnks,maxnrs),fniy(maxnks,maxnrs)
c
      real sum1,sum2,sum
c
      parameter (nplasf=17)
c
c slmod begin - f90
      OPEN(UNIT=NPLASF,ACCESS='SEQUENTIAL',STATUS='REPLACE')
c
c      rewind(nplasf)
c slmod end
c
c slmod begin - f90
c      DO i1 = 0, MAXIX+1
c        DO i2 = 0, MAXIY+1
c          DO i3 = 1, MAXIS
c            ndummy(i1,i2,i3) = 0.0
c          ENDDO
c        ENDDO
c      ENDDO
      CALL RZero(ndummy(0,0,1),(MAXIX+2)*(MAXIY+2)*MAXIS)


      nfla = 2
      IF (stopopt.EQ.106.OR.stopopt.EQ.107.OR.stopopt.EQ.110.OR.
     .    stopopt.EQ.111.OR.stopopt.EQ.113.OR.stopopt.EQ.114.OR.
     .    eirbgk.GE.1) nfla = 6

c slmod end
      nx = mkpts-2
      ny = mrings-2
      nxd = maxix
      nyd = maxiy
c slmod begin
      IF (nx.GE.nxd.OR.ny.GE.nyd) THEN
        WRITE(0,*) 'NX,NY NXD,NYD =',nx,ny,nxd,nyd
        STOP 'ERROR: Insufficient array bounds in B2WRPL'
      ENDIF
      IF (eirdtimv.NE.0.0) THEN
        WRITE(0,*) 
        WRITE(0,*) '*** WHIPING PLASMA BEFORE CALL TO EIRENE ***'
        WRITE(0,*) 
        DO ir = 1, nrs
          DO ik = 1, nks(ir)
          ENDDO
        ENDDO
      ENDIF
c slmod end
c
c-----------------------------------------------------------------------
c
c     Plasma density knbs(ik,ir) (ni)
c     scaled by 1.0
c
      call maptob2(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >             nfla,ndummy(0,0,1),knbs,maxnks,maxnrs,1.0,0)
      WRITE(nplasf,*) '[plasma ion density]'
      CALL DB('[plasma ion density]')
      call gfsub3w(nplasf,nx,ny,nxd,nyd,nfla,ndummy(0,0,1))
c
c-----------------------------------------------------------------------
c
c     Plasma electron density knes(ik,ir) (ne)
c     scaled by 1.0
c
      call maptob2(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >             nfla,ndummy(0,0,1),knes,maxnks,maxnrs,1.0,0)
      WRITE(nplasf,*) '[plasma electron density]'
      CALL DB('[plasma electron density]')
      call gfsub3w(nplasf,nx,ny,nxd,nyd,nfla,ndummy(0,0,1))
c
c-----------------------------------------------------------------------
c
c     poloidal velocity (uu)
c     scaled by qtim*csolvf to obtain m/s
c
c      DO ir = irtrap+1, nrs
c        DO ik = 1, nks(ir)
c          kvhs(ik,ir) = 0.0
c        ENDDO
c      ENDDO

      do i=1,maxnrs
         do j=1,maxnks
            uu(j,i)=kvhs(j,i)*bratio(j,i)
            WRITE(SLOUT,'(A,2I6,1P,3E10.2,0P)') 
     .       'POLOIDAL VEL:',i,j,kvhs(j,i),bratio(j,i),uu(j,i) 
         enddo
      enddo

      CALL OutputData(85,'Poloidal Veloicty')
c      WRITE(0,*) 'CUTRING:',cutring
c

      call maptob2(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >                    nfla,ndummy(0,0,1),uu,maxnks,maxnrs,
     >                    1.0,1)

c slmod begin
      WRITE(nplasf,*) '[poloidal velocity]'
      CALL DB('[poloidal velocity]')
      WRITE(SLOUT ,*) '[poloidal velocity]'
c
      DO i = 0, nrs-1
        DO j = 0, nks(i+1)+1
          IF (.NOT.(ndummy(j,i,1).GE.0.0.OR.ndummy(j,i,1).LT.0.0))
     .      WRITE(0,'(5X,A,2I4,1P)')
     .        'NaNQ = ',j,i
        ENDDO
      ENDDO

c      STOP 'sdfsd'
c slmod end
      call gfsub3w(nplasf,nx,ny,nxd,nyd,nfla,ndummy(0,0,1))
c
c-----------------------------------------------------------------------
c
c     Radial velocity (vv)
c     set to zero
c
      call rzero(ndummy,(maxix+1)*(maxiy+1)*maxis)
c slmod begin
c...03 VVB
      WRITE(nplasf,*) '[radial velocity]'
      CALL DB('[radial velocity]')
c slmod end
      call gfsub3w(nplasf,nx,ny,nxd,nyd,nfla,ndummy(0,0,1))
c
c-----------------------------------------------------------------------
c
c     Electron temperature ktebs(ik,ir) (te)
c     scaled by 1.0/1.6E-19 (eV -> Joule)
c
      call maptob2(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >             1,tdummy(0,0),ktebs,maxnks,maxnrs,1.0/1.6e-19,0)
c slmod begin
c...04
      WRITE(nplasf,*) '[electron temperature]'
      CALL DB('[electron temperature]')
c slmod end
      call gfsub3w(nplasf,nx,ny,nxd,nyd,1,tdummy(0,0))
c
c
c-----------------------------------------------------------------------
c
c     ion temperature ktibs(ik,ir) (ti)
c     scaled by 1.0/1.6e-19 (eV -> Joule)
c
      call maptob2(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >             1,tdummy(0,0),ktibs,maxnks,maxnrs,1.0/1.6e-19,0)
c slmod begin
c...05
      WRITE(nplasf,*) '[ion temperature]'
      CALL DB('[ion temperature]')
c slmod end
      call gfsub3w(nplasf,nx,ny,nxd,nyd,1,tdummy(0,0))
c
c
c-----------------------------------------------------------------------
c
c     (pr)
c     set to zero
c
      call rzero(tdummy,(maxix+1)*(maxiy+1))
c slmod begin
c...06
      WRITE(nplasf,*) '[pr - set to zero]'
      CALL DB('[pr - set to zero]')
c slmod end
      call gfsub3w(nplasf,nx,ny,nxd,nyd,1,tdummy(0,0))
c
c
c-----------------------------------------------------------------------
c
c     parallel velocity kvhs(ik,ir) (up)
c     scaled by qtim*csolvf to obtain m/s
c
      call maptob2(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >                    nfla,ndummy(0,0,1),kvhs,maxnks,maxnrs,
     >                    1.0,1)
c slmod begin
c...07 UPB
      WRITE(nplasf,*) '[parallel velocity]'
      CALL DB('[parallel velocity]')
c slmod end
      call gfsub3w(nplasf,nx,ny,nxd,nyd,nfla,ndummy(0,0,1))
c
c
c-----------------------------------------------------------------------
c
c     B-field ratio bratio(ik,ir) (pit)
c     scaled by 1
c
      call maptob2(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >                    nfla,tdummy(0,0),bratio,maxnks,maxnrs,1.0,0)
c slmod begin
c...08
      WRITE(nplasf,*) '[b-field ratio]'
      CALL DB('[b-field ratio]')
c slmod end
      call gfsub3w(nplasf,nx,ny,nxd,nyd,1,tdummy(0,0))
c
c
c-----------------------------------------------------------------------
c
c     Particle current on target cells (fnix)
c     i.e. parallel flux * target segment length * cos(incl. angle) *
c                          2 * pi * major radius * B-field ratio
c
c     Naming inconsistency in use of fniy - one of (fnix,fniy)
c     represents the parallel flux/current direction - we believe
c     that it is fnix at the moment so the use of fniy below is
c     inconsistent - this conclusion may be backward :-)
c
      call rzero(ndummy,(maxix+1)*(maxiy+1)*maxis)
      call rzero(fniy,maxnrs*maxnks)
c
      sum1 = 0.0
      sum2 = 0.0
      sum = 0.0
c
c     These fluxes should only need to be defined on the actual
c     edges of the end cells - this would be index ik=1 and ik=nks(ir)-1
c     However, at the moment they are also being defined (as the same
c     values) for ik=2 and ik=nks(ir) - so that the value may be obtained
c     no matter how Eirene is interpreting  an array of cell edge data.
C     Note that the array is then mapped as a cell-centred one - since the
c     "Averaging" that is used to calculate the cell boundary data for
c     an entire ring are unneeded here.
c
c
      do 640 ir = irsep, nrs
c
c         First target
c
            ik = 2
            id = idds(ir,2)
c
c           Handle possibly supersonic flow from SOL 22
c
c slmod begin
c... also include SOL22p?
            if (cioptf.eq.22.OR.osm_model(IKLO,ir).EQ.22) then
c            if (cioptf.eq.22.OR.GetModel(IKLO,ir).EQ.22) then

c              WRITE(0,*) 'EIRDIV: SOL22 mode',IKLO,ir
c
c            if (cioptf.eq.22) then
c slmod end
                CS = 9.79E3 * SQRT (0.5*(KTEDS(ID)+KTIDS(ID))*
     >                (1.0+RIZB)/CRMB) * cmachno(ir,2)
            elseIF (cioptg.EQ.92.or.osm_model(IKLO,ir).EQ.28) THEN
c...          This is needed to get the flux correct for
c             supersonic targets:
              cs = ABS(kvds(idds(ir,2)))              
            else
              cs = 9.79e3 * sqrt (0.5*(kteds(id)+ktids(id))*
     >                  (1.0+rizb)/crmb)
            endif
c
c           write (6,997) 'testfl-i:',ir,ik,id,dds(id),costet(id),cs
*           added code for northopt=2; Krieger 6/95
c
            if (northopt.eq.0) then
c             minus sign f. B2
              fniy(ik,ir) =-knds(id)*cs*dds(id)*eirtorfrac*eirsrcmul*
     >                      bratio(ik,ir)*2.*pi*rp(id)
            else if (northopt.eq.1.or.northopt.eq.2.or.northopt.eq.3)
     >                          then
c             minus sign f. B2
c...fix? dds to dds2? kbfst instead of bratio
              fniy(ik,ir) =-knds(id)*cs*dds(id)*eirtorfrac*eirsrcmul*
     >                      bratio(ik,ir)*2.*pi*rp(id)*costet(id)
            endif
c
c slmod begin
c            IF (ir.NE.20) fniy(ik,ir) = -1.0E+20

            IF ((stopopt.EQ.90.OR.stopopt.EQ.100).AND.
c     .          ir.GE.irsep .AND.
c     .          ir.LT.irwall) THEN
c     .          ir.GE.irbreak .AND.
c     .          ir.LE.nrs) THEN
     .          ir.GE.irsep .AND.
     .          ir.LE.nrs) THEN


              WRITE(PINOUT,*) 'SUPPRESSING INNER TARGET FLUX ',ir

              fniy(ik,ir) = 1.0E-06 * fniy(ik,ir)

            ENDIF

            IF ((.FALSE..AND.stopopt .EQ.121).OR.
     .          (.FALSE..AND.eirnpuff.GT.0  ).OR.
     .          stopopt.EQ.130) THEN
              IF (ir.EQ.irsep) 
     .          WRITE(PINOUT,*) 'SUPPRESSING INNER TARGET FLUX '

              fniy(ik,ir) = 1.0E-06 * fniy(ik,ir)
            ENDIF               


c...        Legit development tool:
            IF ((iflexopt(7).EQ.5                 ).OR.
     .          (iflexopt(7).EQ.6.AND.ir.GE.irtrap).OR.
     .          (iflexopt(7).EQ.7.                ).OR.
     .          (iflexopt(7).EQ.8.AND.ir.LE.irwall).OR.
     .          (iflexopt(7).EQ.9.                ).OR.
     .          (supflx(IKLO,ir).EQ.1)) THEN
              WRITE(PINOUT,*) 'SUPPRESSING INNER TARGET FLUX',ir
              fniy(ik,ir) = 1.0E-15 * fniy(ik,ir)
            ENDIF

           IF (iflexopt(5).EQ.31.AND.ir.GT.irtrap) THEN
              WRITE(PINOUT,*) 'SUPPRESSING INNER PFZ TARGET FLUX',ir
             fniy(ik,ir) = 1.0E-15 * fniy(ik,ir)
           ENDIF
c slmod end
            fniy(ik-1,ir)=fniy(ik,ir)
c
c         Second target
c
            ik = nks(ir) -1
            id = idds(ir,1)
c
c           Handle possibly supersonic flow from SOL 22
c
c slmod begin
            if (cioptf.eq.22.OR.GetModel(IKHI,ir).EQ.22) then

c              WRITE(0,*) 'EIRDIV: SOL22 mode',IKHI,ir
c
c            if (cioptf.eq.22) then
c slmod end
               CS = 9.79E3 * SQRT (0.5*(KTEDS(ID)+KTIDS(ID))*
     >                (1.0+RIZB)/CRMB) * cmachno(ir,1)
            elseIF (cioptg.EQ.92.or.osm_model(IKLO,ir).EQ.28) THEN  
c... BUG! SOL28 super-sonics were not recognizd 25.5.2005
c...          This is needed to get the flux correct for
c             supersonic targets:
              cs = ABS(kvds(idds(ir,1)))              
            ELSE
              cs = 9.79e3 * sqrt (0.5*(kteds(id)+ktids(id))*
     >                   (1.0+rizb)/crmb)
            endif


c
c           write (6,997) 'testfl-o:',ir,ik,id,dds(id),costet(id),cs
997        FORMAT(A,3I4,1P,3E12.4,0P)
*           added code for northopt=2; Krieger 6/95
c
            if (northopt.eq.0) then
              fniy(ik,ir) = knds(id)*cs*dds(id)*eirtorfrac*eirsrcmul*
     >                      bratio(ik,ir)*2.*pi*rp(id)
            else if (northopt.eq.1.or.northopt.eq.2.or.northopt.eq.3)
     >                     then
c...fix? dds to dds2?  kbfst instead of bratio
              fniy(ik,ir) = knds(id)*cs*dds(id)*eirtorfrac*eirsrcmul*
     >                      bratio(ik,ir)*2.*pi*rp(id)*costet(id)

c slmod begin
              WRITE(SLOUT,'(A,I6,2F12.6)') 'OUTER TARG FACTOR: ',
     .          ir,rs(ik,ir),bratio(ik,ir)*costet(id)
     .          
c slmod end


c          INFLUX(IR) = KNDS (ID) * CS / KBFST (IR,1) *
c     >                 COSTET (ID) * DDS2 (ID) * cmachno(ir,1)

            endif
c slmod begin
c            IF (ir.NE.20) fniy(ik,ir) = 1.0E+20

            IF ((stopopt.EQ.90.OR.stopopt.EQ.100).AND.
c     .          ir.GE.irsep .AND.
c     .          ir.LT.irwall) THEN
c     .          ir.GE.irbreak .AND.
c     .          ir.LE.nrs) THEN
     .          ir.GE.irsep .AND.
     .          ir.LE.nrs) THEN

              WRITE(PINOUT,*) 'SUPPRESSING OUTER TARGET FLUX ',ir

              fniy(ik,ir) = 1.0E-06 * fniy(ik,ir)

            ENDIF

            IF ((.FALSE..AND.stopopt .EQ.121.AND.
     .           (ir.LT.20.OR.ir.GT.22)).OR.
     .          (.FALSE..AND.eirnpuff.GT.0  ).OR.
     .          stopopt.EQ.130) THEN
              IF (ir.EQ.irsep) 
     .          WRITE(PINOUT,*) 'SUPPRESSING OUTER TARGET FLUX '

              fniy(ik,ir) = 1.0E-06 * fniy(ik,ir)
            ENDIF               
            IF (.FALSE..AND.ir.GE.20.AND.ir.LE.22) THEN
              WRITE(PINOUT,*) 'SUPPRESSING OUTER TARGET FLUX ON ',ir
              fniy(ik,ir) = 1.0E-06 * fniy(ik,ir)
            ENDIF


c...        Legit development tool:

            IF ((iflexopt(7).EQ.4                 ).OR.
     .          (iflexopt(7).EQ.6                 ).OR.
     .          (iflexopt(7).EQ.7.AND.ir.GE.irtrap).OR.
     .          (iflexopt(7).EQ.8                 ).OR.
     .          (iflexopt(7).EQ.9.AND.ir.LE.irwall).OR.
     .          (supflx(IKHI,ir).EQ.1)) THEN
              WRITE(PINOUT,*) 'SUPPRESSING OUTER TARGET FLUX',ir
              fniy(ik,ir) = 1.0E-15 * fniy(ik,ir)
            ENDIF

           IF (iflexopt(5).EQ.31.AND.ir.GT.irtrap) THEN
             WRITE(PINOUT,*) 'SUPPRESSING OUTER PFZ TARGET FLUX',ir
             fniy(ik,ir) = 1.0E-15 * fniy(ik,ir)
           ENDIF
c slmod end
c
            fniy(ik+1,ir)=fniy(ik,ir)
c
            sum1 = sum1 + fniy(2,ir)
            sum2 = sum2 + fniy(nks(ir)-1,ir)
            sum = sum - fniy(2,ir) +  fniy(nks(ir)-1,ir)



            write(6,998) 'Target Flux: ',ir,idds(ir,2),fniy(2,ir),
     >                 idds(ir,1),fniy(nks(ir)-1,ir),sum1,sum2,sum
998         FORMAT(A,I4,2(I4,1P,E12.4,0P),2X,1P,3E12.4,0P)

c
  640 continue
c
c     Although fniy should be a cell-edge based quantity - it is being
c     mapped as a cell centered quantity because it is already being
c     calculated above as a cell-edge quantity and shouldn't need
c     any further mapping.
c
      call maptob2(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >                    nfla,ndummy(0,0,1),fniy,maxnks,maxnrs,1.0,0)
c slmod begin
      WRITE(nplasf,*) '[particle current on targets]'
      WRITE(SLOUT ,*) '[particle current on targets]'
      CALL DB('[particle current on targets]')
c
c      DO i = 0, nrs-1
c        DO j = 0, nks(i+1)+1
c          IF (.NOT.(ndummy(j,i,1).GE.0.0.OR.ndummy(j,i,1).LT.0.0))
c     .      WRITE(SLOUT,'(5X,A,2I4,1P)')
c     .        'NaNQ = ',j,i
c        ENDDO
c      ENDDO
c slmod end
      call gfsub3w(nplasf,nx,ny,nxd,nyd,nfla,ndummy(0,0,1))
c
c
c-----------------------------------------------------------------------
c
c     (fniy)
c     set to zero
c
      call rzero(ndummy,(maxix+1)*(maxiy+1)*maxis)
c slmod begin
c...10
      WRITE(nplasf,*) '[fnix - set to zero]'
      CALL DB('[fnix - set to zero]')
c slmod end
      call gfsub3w(nplasf,nx,ny,nxd,nyd,nfla,ndummy(0,0,1))
c
c
c-----------------------------------------------------------------------
c
c     (feix)
c     set to zero

      call rzero(tdummy,(maxix+1)*(maxiy+1))
c slmod begin
c...11
      WRITE(nplasf,*) '[feix - set to zero]'
      CALL DB('[feix - set to zero]')
c slmod end
      call gfsub3w(nplasf,nx,ny,nxd,nyd,1,tdummy(0,0))
c
c
c-----------------------------------------------------------------------
c
c     (feiy)
c     set to zero
c
      call rzero(tdummy,(maxix+1)*(maxiy+1))
c slmod begin
c...11
      WRITE(nplasf,*) '[feiy - set to zero]'
      CALL DB('[feiy - set to zero]')
c slmod end
      call gfsub3w(nplasf,nx,ny,nxd,nyd,1,tdummy(0,0))
c
c
c-----------------------------------------------------------------------
c
c     (feex)
c     set to zero
c
      call rzero(tdummy,(maxix+1)*(maxiy+1))
c slmod begin
c...12
      WRITE(nplasf,*) '[feex - set to zero]'
      CALL DB('[feex - set to zero]')
c slmod end
      call gfsub3w(nplasf,nx,ny,nxd,nyd,1,tdummy(0,0))

c
c-----------------------------------------------------------------------
c
c     (feey)
c     set to zero
c
      call rzero(tdummy,(maxix+1)*(maxiy+1))
c slmod begin
c...13
      WRITE(nplasf,*) '[feey - set to zero]'
      CALL DB('[feey - set to zero]')
c slmod end
      call gfsub3w(nplasf,nx,ny,nxd,nyd,1,tdummy(0,0))
c
c slmod begin
c-----------------------------------------------------------------------
c
c     (tabrcm)
c     set to 1.0
c

      IF (.NOT..TRUE.) THEN
        WRITE(0     ,*) 'WARNING: OVERRIDE OF LOCAL RECOMBINATION RATE'
        WRITE(PINOUT,*) 'WARNING: OVERRIDE OF LOCAL RECOMBINATION RATE'
        DO ir = irtrap+1, nrs
          ik1 = 0
          temax = 0.0
          DO ik = 1, nks(ir)
            IF (ktebs(ik,ir).GT.temax) THEN
              temax = ktebs(ik,ir)
              ik1 = ik
            ENDIF
          ENDDO
          DO ik = 1, ik1
c          DO ik = ik1, nks(ir)
            mulrec(ik,ir) = 0.0
c            mulrec(ik,ir) = mulrec(ik,ir) * 2.0
          ENDDO
        ENDDO
      ENDIF


c      eiropacity=0

      CALL RSet(ndummy,(maxix+1)*(maxiy+1)*maxis,1.0)

      IF     (iflexopt(5).EQ.31) THEN
c...    Note opacity multiplier gets passed even if eiropacity.EQ.0:

        IF (iflexopt(5).EQ.31) THEN  
c...      Turn off vol. rec. in the PFZ:
          WRITE(PINOUT,*) 'SUPPRESSING PFZ VOL. REC.',ir
          DO ir = irtrap+1, nrs
            DO ik = 1, nks(ir)
              mulrec(ik,ir) = 1.0E-08
            ENDDO  
          ENDDO
          WRITE(0,*) '***************************'
          WRITE(0,*) '* SUPPRESSING PFZ SOURCES *'
          WRITE(0,*) '***************************'
        ENDIF

        call maptob2(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >               nfla,ndummy(0,0,1),mulrec,maxnks,maxnrs,1.0,0)

      ELSEIF (stopopt2.EQ.3) THEN
        STOP 'STOPOPT2.EQ.3 NO LONGER SUPPORTED'
        DO ir = irsep, irwall-1
          DO ik = 1, SymmetryPoint(ir)
            mulrec(ik,ir) = 1.0
          ENDDO
        ENDDO
        call maptob2(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >               nfla,ndummy(0,0,1),mulrec,maxnks,maxnrs,1.0,0)
      ELSEIF (stopopt2.EQ.2.OR.stopopt2.EQ.4.OR.stopopt2.EQ.5.OR.
     .        eiropacity.EQ.1.OR.eiropacity.EQ.2.OR.
     .        eiropacity.EQ.3.OR.eiropacity.EQ.4.OR.
     .        eiropacity.EQ.5) THEN
        call maptob2(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >               nfla,ndummy(0,0,1),mulrec,maxnks,maxnrs,1.0,0)
      ELSEIF (eiropacity.EQ.-1) THEN
        STOP 'EIROPACITY.EQ.-1 NO LONGER SUPPORTED'
        CALL RSet(ndummy,(maxix+1)*(maxiy+1)*maxis,1.0E-15)
      ELSEIF (eiropacity.EQ.-2.OR.eiropacity.EQ.-3) THEN
c...    Non-standard use of the opacity multiplier.  It is used to turn off
c       recombination on specified rings -- or -- SOL28 recombination scaling:
        call maptob2(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >               nfla,ndummy(0,0,1),mulrec,maxnks,maxnrs,1.0,0)
      ELSEIF (stopopt2.EQ.1) THEN
        STOP 'STOPOPT2.EQ.1 NO LONGER SUPPORTED'
        CALL RSet(ndummy,(maxix+1)*(maxiy+1)*maxis,0.5)
      ENDIF
c...13
c      CALL RSet(ndummy,(maxix+1)*(maxiy+1)*maxis,1.0)

      WRITE(nplasf,*) '[mulrec]'
      CALL DB('[mulrec]')
      call gfsub3w(nplasf,nx,ny,nxd,nyd,nfla,ndummy(0,0,1))
c
c     ------------------------------------------------------------------
c
      CALL RSet(ndummy,(maxix+1)*(maxiy+1)*maxis,1.0)

      IF     (iflexopt(6).EQ.12) THEN
c...    Turn off ionisation in the PFZ:
        WRITE(0,*) 'TURNING OFF IONISATION IN THE PFZ'
        DO ir = irtrap+1, nrs
          DO ik = 1, nks(ir)
            mulion(ik,ir) = 1.0E-10
          ENDDO
        ENDDO
        call maptob2(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >               nfla,ndummy(0,0,1),mulion,maxnks,maxnrs,1.0,0)
      ELSEIF (iflexopt(6).EQ.13) THEN
c...    Turn off ionisation in the outer SOL:
        WRITE(0,*) 'TURNING OFF IONISATION IN THE OUTER SOL'
        DO ir = irsep, irwall-1
          DO ik = nks(ir)/2, nks(ir)
            mulion(ik,ir) = 1.0E-10
          ENDDO
        ENDDO
        call maptob2(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >               nfla,ndummy(0,0,1),mulion,maxnks,maxnrs,1.0,0)
      ELSEIF (iflexopt(6).EQ.14) THEN
c...    Turn off ionisation in the outer SOL and PFZ:
        WRITE(0,*) 'TURNING OFF IONISATION NEAR OUTER TARGET'
        DO ir = irsep, nrs
          DO ik = nks(ir)/2, nks(ir)
            IF (zs(ik,ir).LT.zxp) mulion(ik,ir) = 1.0E-10
          ENDDO
        ENDDO
        call maptob2(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >               nfla,ndummy(0,0,1),mulion,maxnks,maxnrs,1.0,0)
      ELSEIF (stopopt2.EQ.3) THEN
        DO ir = irsep, irwall-1
          DO ik = 1, SymmetryPoint(ir)
            mulion(ik,ir) = 1.0
          ENDDO
        ENDDO
c...    Turn off ionisation in the PFZ:
        IF (iflexopt(6).EQ.12) THEN
          DO ir = irtrap+1, nrs
            DO ik = 1, nks(ir)
              mulion(ik,ir) = 1.0E-10
            ENDDO
          ENDDO
        ENDIF
        call maptob2(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >               nfla,ndummy(0,0,1),mulion,maxnks,maxnrs,1.0,0)
      ELSEIF (stopopt2.EQ.2.OR.stopopt2.EQ.4.OR.eiropacity.EQ.2.OR.
     .        eiropacity.EQ.4.OR.eiropacity.EQ.6.OR.
     .        eiropacity.EQ.-4) THEN

        call maptob2(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >               nfla,ndummy(0,0,1),mulion,maxnks,maxnrs,1.0,0)

      ELSEIF (stopopt2.EQ.1) THEN
        CALL RSet(ndummy,(maxix+1)*(maxiy+1)*maxis,0.5)
      ENDIF
c...13
      WRITE(nplasf,*) '[mulion]'
      CALL DB('[mulion')
      call gfsub3w(nplasf,nx,ny,nxd,nyd,nfla,ndummy(0,0,1))
c
c     MULQER, the opacity correction multiplier for electron power loss
c     due to recombination:
c
      CALL RSet(ndummy,(MAXIX+1)*(MAXIY+1)*MAXIS,1.0)
      IF (stopopt2.EQ.4.OR.stopopt2.EQ.5.OR.
     .    eiropacity.EQ.1.OR.eiropacity.EQ.2.OR.eiropacity.EQ.3.OR.
     .    eiropacity.EQ.4) THEN
        call maptob2(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >               nfla,ndummy(0,0,1),mulqer,maxnks,maxnrs,1.0,0)
      ENDIF

      WRITE(nplasf,*) '[mulqer]'
      CALL DB('[mulqer]')
      call gfsub3w(nplasf,nx,ny,nxd,nyd,nfla,ndummy(0,0,1))
c
c     MULQEI, the opacity correction multiplier for electron power loss
c     due to ionisation:
c
      CALL RSet(ndummy,(MAXIX+1)*(MAXIY+1)*MAXIS,1.0)
      IF (stopopt2.EQ.4.OR.eiropacity.EQ.2.OR.eiropacity.EQ.4) THEN
        call maptob2(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >               nfla,ndummy(0,0,1),mulqei,maxnks,maxnrs,1.0,0)
      ENDIF

      WRITE(nplasf,*) '[mulqei]'
      CALL DB('[mulqei]')
      call gfsub3w(nplasf,nx,ny,nxd,nyd,nfla,ndummy(0,0,1))
c
c     ------------------------------------------------------------------
c
c     (diin,tiin,vxin,vyin,vzin from PINBGK)
c
      CALL RZero(ndummy,(MAXIX+1)*(MAXIY+1)*MAXIS)



c...HUGE BUG! Well, it seems as if the standard grid BGK data was not being
c   sent from EIRENE to DIVIMP.  This explains a lot of things, such as why the
c   BGK parameters did not seem to be converging.  I cannot recall when the nfla.EQ.6
c   flag was replaced with eirbgk.EQ.2, but it was quite a while ago.  Anyway,
c   I will have to re-run cases and see how the DIII-D cases are affected. 
c   SL- June 7, 2000
      IF ((eirbgk.EQ.1.OR.eirbgk.EQ.2.OR.eirbgk.EQ.3).AND.
     .    stopopt3.NE.14) THEN
c      IF (eirbgk.GE.1) THEN
c      IF (eirbgk.EQ.2) THEN

c...need a better flag for when including n-n collisions
c      IF (nfla.EQ.6) THEN
c        WRITE(0,*) 'WRITING BGK DATA TO EIRENE'
        DO in = 1, eirnsdtor
          DO i1 = 1, 5
            DO i2 = 1, nfla
              DO ir = 1, MAXNRS
                DO ik = 1, MAXNKS
                  vdummy(ik,ir)=pinbgk(ik,ir,(i2-1)*5+i1+(in-1)*MAXBGK)
                ENDDO
              ENDDO
              CALL MapToB2(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,nfla,
     .                     ndummy(0,0,i2),vdummy,MAXNKS,MAXNRS,1.0,0)
            ENDDO
            WRITE(nplasf,*) '[BGK DATA: ',in,i1,']'
            CALL Gfsub3w(nplasf,nx,ny,nxd,nyd,nfla,ndummy(0,0,1))
          ENDDO
        ENDDO
      ELSE
        IF (eirbgk.GT.0) THEN
          WRITE(0     ,*) 'BLANKING STANDARD GRID BGK DATA PASSED '//
     .                    'TO EIRENE'
          WRITE(PINOUT,*) 'BLANKING STANDARD GRID BGK DATA PASSED '//
     .                    'TO EIRENE'
        ENDIF
        DO in = 1, eirnsdtor
          DO i1 = 1, 5
            WRITE(nplasf,*) '[BGK DATA: ',i1,in,' (set to zero)]'
            CALL Gfsub3w(nplasf,nx,ny,nxd,nyd,nfla,ndummy(0,0,1))
          ENDDO
        ENDDO
      ENDIF
c
c-----------------------------------------------------------------------
c
      ind1 = 1

      IF (eirbgk.EQ.2.OR.eirbgk.EQ.3.OR.eirbgk.EQ.4) THEN
        IF (stopopt3.EQ.14.OR.stopopt.EQ.122) THEN
c        IF (iflexopt(8).EQ.1) THEN
c...BUG: This used to blank the additional cell data when a c08 grid
c        with extended field lines was being used.  I didn't realize this
c        when running some initial vacuum grid cases, so the results were
c        bogus.  Anyway, fixed now and the cases are being re-run. 
c        SL - June 23, 2000
          WRITE(0     ,*) 'BLANKING ADD CELL DATA PASSED TO EIRENE'
          WRITE(PINOUT,*) 'BLANKING ADD CELL DATA PASSED TO EIRENE'
          CALL RZero(pinasd(1,1,1,ind1),MAXASCDAT*MAXASD2*MAXASS)
        ELSE
          WRITE(PINOUT,*) 'PASSING ADDITIONAL CELL DATA TO EIRENE'
        ENDIF

c...NEW! TESTING!
        WRITE(PINOUT,*) 'BLANKING ADDITIONAL CELL BGK DATA FOR '//
     .                  'PRESSURE GAUGE(S)'
        DO i1 = 1, 10
          DO i2 = 1, 1+eirnpgdat
            pinasd(i2,1,i1,ind1) = 1.0E+02
            pinasd(i2,2,i1,ind1) = 1.0E-02
            pinasd(i2,3,i1,ind1) = 0.0
            pinasd(i2,4,i1,ind1) = 0.0
            pinasd(i2,5,i1,ind1) = 0.0
          ENDDO
        ENDDO

        CALL AssignAdditionalCellPlasma


c...TEMP
c        DO i1 = 7, 7
c          DO i2 = 1, 10
c            WRITE(0,*)  i1,i2,(pinasd(i2,i3,i1,ind1),i3=1,2)
c          ENDDO
c        ENDDO

        nradd = 1 + eirnpgdat + asc_ncell*ascncut

        WRITE(nplasf,*) 'ADDITIONAL CELL BGK DATA'
        WRITE(nplasf,*) nradd,indasd

        DO i1 = 5, 10
          WRITE(nplasf,'(1X,A5,5A15,I4)')
     .      'IN','DIIN','TIIN','VXIN','VYIN','VZIN',i1
          DO i2 = 1, nradd

c...!!! 1 should be 2 after relaxation implimented ??? Maybe not... !!!
            WRITE(nplasf,'(I6,1P,5E15.7)')
     .        i2+indasd,(pinasd(i2,i3,i1,ind1),i3=1,5)
          ENDDO
        ENDDO

c        WRITE(nplasf,*) 'ADDITIONAL CELL BGK DATA'
c        WRITE(nplasf,*) nradd,5146
c
c        DO i1 = 5, 10
c          WRITE(nplasf,'(A5,5A15,I4)')
c     .      'IN','DIIN','TIIN','VXIN','VYIN','VZIN',i1
c          DO i2 = 1, nradd
cc...!!! 1 should be 2 after relaxation implimented ??? Maybe not... !!!
c            WRITE(nplasf,'(I6,1P,5E15.7)')
c     .        1,(2.0,i3=2,6)
c          ENDDO
c        ENDDO

      ENDIF


c...  Write parameters for recording time-to-ionisation data:
      WRITE(nplasf,'(A)') '[TIME-TO-IONISTAION PARAMETERS]'
      WRITE(nplasf,'(A)') '* ...'
      WRITE(nplasf,'(I6)') eirniontime
      DO i1 = 1, eirniontime 
        WRITE(nplasf,'(I6,1P,4E16.6,0P,I6)') 
     .    i1,(eiriontime(i1,i2),i2=1,4),
     .    NINT(eiriontime(i1,5))
        IF (eiriontime(i1,5).GT.0.0) THEN
          IF     (eiriontime(i1,6).EQ.1.0) THEN
            t1 = eiriontime(i1,7)
            t2 = eiriontime(i1,8)
            tn = eiriontime(i1,5)
          ELSEIF (eiriontime(i1,6).EQ.2.0) THEN
c...        LOG10 mode:
            t1 = LOG10(eiriontime(i1,7))
            t2 = LOG10(eiriontime(i1,8))
            tn = eiriontime(i1,5)
          ELSE
            CALL ER('...','Invalid time-to-ionisation bin mode',*99)
          ENDIF

          WRITE(nplasf,'(6X,1P,E16.6)') 0.0

          DO i2 = 1, NINT(tn+1)
            IF     (eiriontime(i1,6).EQ.1.0) THEN
c...          Linear mode:
              t3 = t1 + (t2 - t1) * REAL(i2-1) / tn
              WRITE(nplasf,'(6X,1P,E16.6)') t3
            ELSEIF (eiriontime(i1,6).EQ.2.0) THEN
c...        LOG10 mode:
              t3 = 10.0**(t1 + (t2 - t1) * REAL(i2-1) / tn)
              WRITE(nplasf,'(6X,1P,E16.6)') t3
            ENDIF
          ENDDO

          WRITE(nplasf,'(6X,1P,E16.6)') HI

        ENDIF
      ENDDO
      WRITE(nplasf,*)


c slmod end
c
c-----------------------------------------------------------------------
c
c     we need to close the file because the buffers must be flushed
c
      close(nplasf)
c
c slmod begin

      IF (stopopt.EQ.22.AND.rel_step.GT.0)
     .  STOP 'After writing PIN files'
c slmod end
      return
c slmod begin
99    STOP
c slmod end
      end
c
c ======================================================================
c
c  subroutine: ReadEIRE
c
      subroutine readeire

      implicit none
      include 'params'
      include 'pindata'
      include 'comtor'
      include 'cgeom'
c slmod begin
      INCLUDE 'slcom'

      LOGICAL output
c slmod end
c
c     readeire: this subroutine is designed to read in the data
c     required from EIRENE output files in the specific data format
c     used by those files and to scale them to the units used in
c     DIVIMP
c
      integer ir,ik,topik
c slmod begin - david
      real sum,sum1,sum2,sum1a,sum2a,sum3,sum3a,totrec
c
c      real sum,sum1,sum2,sum1a,sum2a,sum3,sum3a
c slmod end
c
c     nfla = number of ion species in the Braams data file;
c     for the EIRENE file (which is in similar format, I guess
c     it has to be 1)

      integer maxix,maxiy,maxis
c slmod begin

c...bug
c
c   Karl changed the abvoe nfla to nfla2 but I missed this in the incorporation
c of his code changes (but I am not sure it is best to add it)... please review

c
c Need to make this larger!
c
c
      character coment*100

c      INCLUDE 'params'

c...trouble if maxiy is more than maxnrs-2, since in theory the
c   EIRENE data can span 50 rings, but then 2 more are added...something like
c   this anyway - Jun 9, 1999
      parameter (maxix=MAXNKS,maxiy=100,maxis=MAXNFLA)
c      parameter (maxix=110,maxiy=50,maxis=7)

      INTEGER SymmetryPoint,GetModel

      INTEGER i1,i2,i3,in,ikm,id
      CHARACTER*72 cdum1

      INTEGER    MAXSTRAT,MAXSDATA
      PARAMETER (MAXSTRAT=10,MAXSDATA=10)
      COMMON /EIRCOM/
     .        eir_07nsrc,eir_07opt ,eir_07stra,
     .        eir_07ind1,
     .        eir_07ind2,eir_07ind3,
     .        eir_07wght           ,eir_07data
      INTEGER eir_07nsrc,eir_07opt ,eir_07stra(MAXSTRAT),
     .        eir_07ind1(MAXSTRAT),
     .        eir_07ind2(MAXSTRAT),eir_07ind3(MAXSTRAT)
      REAL    eir_07wght           ,eir_07data(MAXSTRAT,MAXSDATA)

c
c      integer nfla
c      parameter (nfla=1)
c      parameter (maxix=100,maxiy=50,maxis=7)
c slmod end

      real ndummy(0:maxix+1,0:maxiy+1,maxis)
      real tdummy(0:maxix+1,0:maxiy+1)
      real signv
      integer nx,ny,nxd,nyd
c
c     needed for temporary vars
c
      real pinbuf(maxnks,maxnrs)
c
c     calls to read routine
c
c     BEWARE: the EIRENE pin file only writes out the fields without
c             the dead cells, i.e the number of rings and knots is
c             the numbers of the used grid -2 !
c
c     NOTE: -4 is used here because -
c              1) There are two dead cells on each ring that are
c                 not used in the Eirene version of the B2 file.
c                 This accounst for -2
c              2) The numbers maxkpts and maxrings are set
c                 in DIVIMP for use in loops from 1 to maxrings or
c                 maxkpts. Whereas the reading routine gfsub3r uses
c                 loops running from 0 to nx +1 - and from 0 to
c                 ny +1 - thus to get the correct number of both
c                 knots and rings in the read routine for the Eirene
c                 data - another -2 must be taken from the indices.
c
c
c slmod begin - Karl


c...TEMP
c      CALL DUMPGRID('EIRENE99 recovery')




      nfla = 2
      IF (stopopt.EQ.106.OR.stopopt.EQ.107.OR.stopopt.EQ.110.OR.
     .    stopopt.EQ.111.OR.stopopt.EQ.113.OR.stopopt.EQ.114.OR.
     .    eirbgk.GE.1) nfla = 6

c slmod end
      nx = maxkpts-4
      ny = maxrings-4
      nxd = maxix
      nyd = maxiy
c slmod begin
c      output = .TRUE.
      output = .FALSE.

      IF (nx.GE.nxd.OR.ny.GE.nyd) CALL ER('ReadEIRE','Insufficient '//
     .                                    'array bounds',*99)
c slmod end
c
c-----------------------------------------------------------------------
c
c     read in data from pin passing file
c
c-----------------------------------------------------------------------
c slmod begin
c
c July 16, 97 - There is a problem here, as ususal, in that the grid
c needs to be contracted before the plasma quantities are mapped
c to the DIVIMP arrays by Maptodiv... wait, maybe this isn't true.
c ... not true.
c
c      CALL Warn('readeire','Errors in boundary quantities.',0,0,0)
c slmod end
c
c     zero arrays (should be done below)
c
      call rzero (pinvdist,3*14*maxnks*maxnrs)
c
c slmod begin - f90
      OPEN(UNIT=IPINOUT,ACCESS='SEQUENTIAL',STATUS='OLD')
c
c      rewind(ipinout)
c slmod end
c
c-----------------------------------------------------------------------
c
c     read in the pin data and map them from B2 storing convention
c     into divimp format arrays. note also that eirene passes in
c     units per cell so conversion with cell volume is necessary
*     Krieger IPP/96: I think we made an error here:
*     1. ion source had wrong unit (see below, does not really matter
*                                   if flux tubes are renormalized)
*        from Ampere to particles/s -> divide by 1.602e-19
*                                      (done in EIRENE now)
*     2. as cell volumes, we have to use kvol2/kvols instead of
*        karea2/kareas as in our first version
*
*        this was corrected IPP/96 ...

c     Krieger, IPP/97

c     the structure of the pin data file written by EIRENE is now as
c     follows:
c     - H+ source                [particles/(cell x sec)]     PINION
c     - Impurity ion source      [particles/(cell x sec)]     PINIONZ
c     - Hydrogen momentum source [(kg x m)/(cell x sec^2)]    PINMP
c       This is a cell edge quantity and needs special treatment
c       for DIVIMP presumably.
c       Sign convention like for velocities.
c slmod begin
c       An email from Detlev Reiter (March 17, 1999) states that the momentum
c       source is a cell averaged quantity, not a cell edge quantity. -SL
c slmod end
c     - Impurity momentum source [(kg x m)/(cell x sec^2)]
c     - H neutral density        [particles/cell]             PINATOM
c     - Impurity neutral density [particles/cell]             PINZ0
c     - H2 density               [molecules/cell]             PINMOL
c...Changed on Dec 10, 1999 - SL
c     - H2+ density              [molecules/cell]             PINMOI
c OLD     - H2+ density              [molecules/cell]            +PINMOL
c     - Avg. neutr. H energy     [eV]                         PINENA
c     - Avg. neutr. imp. energy  [ev]                         PINENZ
c     - Avg. H2 energy           [eV]                         PINENM
c     - Avg. H2+ energy          [eV]  !!! (not yet added to DIVIMP) !!!
c     - Halpha emissivity        [photons/(cell x sec)]       PINALPHA
c       (5 blocks for different generation mechanisms)
c     - H rec. source            [particles/(cell x sec)]     PINREC
c     - H+ (+imp.) energy source [Watt/cell]                  PINQI
c     - Electron energy source   [Watt/cell]                  PINQE
c
c     Neutral particle velocity distribution information from nimbus
c     (pinvdist) has no EIRENE equivalent yet
c
c     NB: the EIRENE output comes already with the dead cells removed
c         at both targets. Therefore, if we have a target option with
c         dead cells included (like ctargopt=4), we need to add these!
c
c slmod begin
c-----------------------------------------------------------------------
c
c     Stratum data for ionisation, neutral atom density, neutral
c     molecule density, as well as momentum channel data:
c
      CALL RZero(pindata,MAXNKS*MAXNRS*MAXDATA)

      IF (.TRUE..OR.eir_07opt.EQ.2.OR.eir_07opt.EQ.5) THEN
        IF (output) WRITE(0,*) 'NEW PINDATA TRANSFER METHOD'
c        STOP 'BLAHJDFHKSDHF'

        DO i1 = 1, eir_07nsrc
          DO i2 = 1, 3
            WRITE(cdum1,'(72X)')
            READ (IPINOUT,*) cdum1
            IF (output) WRITE(0,*) cdum1

            call rzero  (pinbuf,maxnks*maxnrs)
            call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
            call deadadd(tdummy,nx,ny,nxd,nyd)
            call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     .                    1,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,
     .                    1.0,0)

            IF (.TRUE.) THEN
c...          Just record ionisation data for each stratum -- temporary:
              IF (i2.EQ.1) THEN
c                IF (output) WRITE(0,*) 'STRATUM NO.=',eir_07stra(i1)
                i3 = eir_07stra(i1)
                call sum_contrib(pindata(1,1,i3),pinbuf,nrs,nks,2)
              ENDIF
            ELSE
c...          Default method:
              i3 = eir_07stra(i1) + 4 * (i2 - 1)
              call sum_contrib(pindata(1,1,i3),pinbuf,nrs,nks,2)
            ENDIF

          ENDDO
        ENDDO

        DO i1 = 1, 9
c
c         normalize from "per cell" to "per volume"
c
          call eire_renorm(pindata(1,1,i1),kvols,nrs,nks,irsep,sum,0)
        ENDDO
c
c
c
        IF (output) WRITE(0,*) 'LOOKING FOR MOMENTUM'

        DO i1 = H_MP1, H_MP1+NMOMCHA-1
          READ (IPINOUT,*) cdum1
          IF (output)
     .    WRITE(0      ,*) 'MOMENTUM',i1,cdum1

          call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
          call deadadd(tdummy,nx,ny,nxd,nyd)
          call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     .                  1,maxnfla,tdummy(0,0),pindata(1,1,i1),
     .                  maxnks,maxnrs,1.0,0)
c
c         normalize from "per cell" to "per volume"
c
          call eire_renorm(pindata(1,1,i1),kvols,nrs,nks,irsep,sum,0)
        ENDDO

      ELSE
        IF (output) WRITE(0,*) 'STANDARD METHOD OF READING PINDATA'

        IF (.TRUE.) THEN
          DO i1 = 1, eirnstrata+eirnpuff
            DO i2 = 1, 3
              READ (IPINOUT,*) cdum1
              IF (output)
     .        WRITE(0      ,*) '-> STRATUM '//cdum1            

              call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),
     .                     0)
              call deadadd(tdummy,nx,ny,nxd,nyd)
              call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     .                      1,maxnfla,tdummy(0,0),pinstrata(1,1,i2,i1),
     .                      maxnks,maxnrs,1.0,0)
              call eire_renorm(pinstrata(1,1,i2,i1),kvols,nrs,nks,
     .                         irsep,sum,0)
            ENDDO
          ENDDO

c...      Copy over ionistaion data:
          DO ir = 1, nrs
            DO ik = 1, nks(ir)
              pindata(ik,ir,H_ION1) = pinstrata(ik,ir,1,1)
              pindata(ik,ir,H_ION2) = pinstrata(ik,ir,1,2)
              pindata(ik,ir,H_ION3) = pinstrata(ik,ir,1,3)
            ENDDO
          ENDDO

          DO i1 = 1, NMOMCHA

            READ (IPINOUT,*) cdum1
            IF (output)
     .      WRITE(0      ,*) '-> MOMENTUM '//cdum1            
	    
            call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
            call deadadd(tdummy,nx,ny,nxd,nyd)
            call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     .                    1,maxnfla,tdummy(0,0),pinploss(1,1,i1),
     .                    maxnks,maxnrs,1.0,0)
            call eire_renorm(pinploss(1,1,i1),kvols,nrs,nks,
     .                       irsep,sum,0)

          ENDDO

        ELSE

          DO i1 = 1, H_MP1+NMOMCHA-1-3

            READ (IPINOUT,*) cdum1
            IF (output)
     .      WRITE(0      ,*) '->  MOMENTUM '//cdum1
  
            i2 = i1

            IF (i1.EQ. 1) i2 = H_ION1
            IF (i1.EQ. 2) i2 = H_ATM1
            IF (i1.EQ. 3) i2 = H_MOL1
            IF (i1.EQ. 4) i2 = H_ION2
            IF (i1.EQ. 5) i2 = H_ATM2
            IF (i1.EQ. 6) i2 = H_MOL2
            IF (i1.EQ. 7) i2 = H_ION3
            IF (i1.EQ. 8) i2 = H_ATM3
            IF (i1.EQ. 9) i2 = H_MOL3
            IF (i1.GT. 9) i2 = H_MP1+(i1-9)-1

            call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
            call deadadd(tdummy,nx,ny,nxd,nyd)
            call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     .                    1,maxnfla,tdummy(0,0),pindata(1,1,i2),
     .                    maxnks,maxnrs,1.0,0)
c
c           normalize from "per cell" to "per volume"
c
            call eire_renorm(pindata(1,1,i2),kvols,nrs,nks,irsep,sum,0)
          ENDDO
        ENDIF
      ENDIF
c
c slmod end
c-----------------------------------------------------------------------
c
c     H+ source
c
c slmod begin
      READ (IPINOUT,*) cdum1
      IF (output)
     .WRITE(0      ,*) 'IONISATION '//cdum1

      WRITE(SLOUT,*) 'IONISATION '//cdum1
      call rzero (pinion,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,tdummy(0,0),pinion,maxnks,maxnrs,1.0,0)
c
c      call rzero (pinion,maxnks*maxnrs)
c      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla,maxnfla,tdummy(0,0),0)
c      call deadadd(tdummy,nx,ny,nxd,nyd)
c      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
c     >              nfla,maxnfla,tdummy(0,0),pinion,maxnks,maxnrs,1.0,0)
c slmod end
c
c     normalize from "per cell" to "per volume"
c
      call eire_renorm(pinion,kvols,nrs,nks,irsep,sum,0)
c
      write(6,*) 'Total ionizations per second (PIN): ',sum
c
c-----------------------------------------------------------------------
c
c     Impurity ion source
c
c slmod begin
      call rzero (pinionz,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla-1,maxnfla,ndummy(0,0,1),0)
      call deadadd(ndummy(0,0,1),nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,ndummy(0,0,1),pinionz,maxnks,maxnrs,1.0,0)
c
c      call rzero (pinionz,maxnks*maxnrs)
c      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla,maxnfla,tdummy(0,0),0)
c      call deadadd(tdummy,nx,ny,nxd,nyd)
c      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
c     >             nfla,maxnfla,tdummy(0,0),pinionz,maxnks,maxnrs,1.0,0)
c slmod end
c
c     normalize from "per cell" to "per volume"
c
      call eire_renorm(pinionz,kvols,nrs,nks,irsep,sum,0)
      write(6,*) 'Total imp. ionizations per second (PIN): ',sum
c
c-----------------------------------------------------------------------
c
c     Hydrogen momentum source
c
c slmod begin
      READ (IPINOUT,*) cdum1
      IF (output)
     .WRITE(0      ,*) cdum1

      call rzero (pinmp,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,tdummy(0,0),pinmp,maxnks,maxnrs,1.0,0)
c
c      call rzero (pinqi,maxnks*maxnrs)
c      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla,maxnfla,tdummy(0,0),0)
c      call deadadd(tdummy,nx,ny,nxd,nyd)
c      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
c     >              nfla,maxnfla,tdummy(0,0),pinmp,maxnks,maxnrs,1.0,0)
c slmod end
c
c     normalize from "per cell" to "per volume"
c
      call eire_renorm(pinmp,kvols,nrs,nks,irsep,sum,0)
c...temp
c      WRITE(0,*) 'BLANKING PINMP'
c      call rzero (pinmp,maxnks*maxnrs)
c
c-----------------------------------------------------------------------
c
c     Impurity momentum source  (NOT USED YET)
c
c slmod begin
      call rzero (pinbuf,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla-1,maxnfla,ndummy(0,0,1),0)
      call deadadd(ndummy(0,0,1),nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,ndummy(0,0,1),pinbuf,maxnks,maxnrs,1.0,0)
c
c      call rzero (pinbuf,maxnks*maxnrs)
c      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla,maxnfla,tdummy(0,0),0)
c      call deadadd(tdummy,nx,ny,nxd,nyd)
c      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
c     >              nfla,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
c slmod end
c
c     normalize from "per cell" to "per volume"
c
      call eire_renorm(pinbuf,kvols,nrs,nks,irsep,sum,1)
c
c
c-----------------------------------------------------------------------
c
c     H neutral density
c
c slmod begin
      READ (IPINOUT,*) cdum1
      IF (output)
     .WRITE(0      ,*) cdum1

      tagpinatom = .TRUE.

      call rzero (pinatom,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,tdummy(0,0),pinatom,maxnks,maxnrs,1.0,0)

      IF (eiropacity.GT.0) CALL EstimateOpacityMultiplier      

c
c      call rzero (pinatom,maxnks*maxnrs)
c      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla,maxnfla,tdummy(0,0),0)
c      call deadadd(tdummy,nx,ny,nxd,nyd)
c      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
c     >             nfla,maxnfla,tdummy(0,0),pinatom,maxnks,maxnrs,1.0,0)
c slmod end
c
c     normalize from "per cell" to "per volume"
c
      call eire_renorm(pinatom,kvols,nrs,nks,irsep,sum,0)
c
c-----------------------------------------------------------------------
c
c     Impurity neutral density (optional for certain launch options)
c
c slmod begin
      call rzero (pinz0,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,ndummy(0,0,1),0)
c      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla-1,maxnfla,ndummy(0,0,1),0)
      call deadadd(ndummy(0,0,1),nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,ndummy(0,0,1),pinz0,maxnks,maxnrs,1.0,0)
c
c      call rzero (pinz0,maxnks*maxnrs)
c      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla,maxnfla,tdummy(0,0),0)
c      call deadadd(tdummy,nx,ny,nxd,nyd)
c      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
c     >              nfla,maxnfla,tdummy(0,0),pinz0,maxnks,maxnrs,1.0,0)
c slmod end
c
c     normalize from "per cell" to "per volume"
c
      call eire_renorm(pinz0,kvols,nrs,nks,irsep,sum,0)
c
c-----------------------------------------------------------------------
c
c     H2 and H2+ molecular density, summed up in one array
c
c     first H2
c
c slmod begin
      READ (IPINOUT,*) cdum1
      IF (output)
     .WRITE(0      ,*) cdum1

      call rzero (pinmol,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,tdummy(0,0),pinmol,maxnks,maxnrs,1.0,0)

      call eire_renorm(pinmol,kvols,nrs,nks,irsep,sum,0)
c
c     then H2+
c
      READ (IPINOUT,*) cdum1
      IF (output)
     .WRITE(0      ,*) cdum1

      call rzero (pinmoi,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,tdummy(0,0),pinmoi,maxnks,maxnrs,1.0,0)
      call eire_renorm(pinmoi,kvols,nrs,nks,irsep,sum,0)

c
c      call rzero (pinmol,maxnks*maxnrs)
c      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla,maxnfla,tdummy(0,0),0)
c      call deadadd(tdummy,nx,ny,nxd,nyd)
c      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
c     >              nfla,maxnfla,tdummy(0,0),pinmol,maxnks,maxnrs,1.0,0)
c
c     then H2+
c
c      call rzero (pinbuf,maxnks*maxnrs)
c      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla,maxnfla,tdummy(0,0),0)
c      call deadadd(tdummy,nx,ny,nxd,nyd)
c      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
c     >              nfla,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
c     add up all contributions, then normalize ...
c      call sum_contrib(pinmol,pinbuf,nrs,nks,2)
c
c      call eire_renorm(pinmol,kvols,nrs,nks,irsep,sum,0)
c slmod end
c
c-----------------------------------------------------------------------
c
c     Average neutral H energy
c
c slmod begin
      READ (IPINOUT,*) cdum1
      IF (output)
     .WRITE(0      ,*) cdum1

      call rzero (pinena,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,tdummy(0,0),pinena,maxnks,maxnrs,1.0,0)
c
c      call rzero (pinena,maxnks*maxnrs)
c      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla,maxnfla,tdummy(0,0),0)
c      call deadadd(tdummy,nx,ny,nxd,nyd)
c      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
c     >              nfla,maxnfla,tdummy(0,0),pinena,maxnks,maxnrs,1.0,0)
c slmod end
c
c     normalize from "per cell" to "per volume"
c
      call eire_renorm(pinena,kvols,nrs,nks,irsep,sum,1)
c
c-----------------------------------------------------------------------
c
c     Average neutral imp. energy
c
c slmod begin
      call rzero (pinenz,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,ndummy(0,0,1),0)
c      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla-1,maxnfla,ndummy(0,0,1),0)
      call deadadd(ndummy(0,0,1),nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,ndummy(0,0,1),pinenz,maxnks,maxnrs,1.0,0)
c
c      call rzero (pinenz,maxnks*maxnrs)
c      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla,maxnfla,tdummy(0,0),0)
c      call deadadd(tdummy,nx,ny,nxd,nyd)
c      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
c     >              nfla,maxnfla,tdummy(0,0),pinenz,maxnks,maxnrs,1.0,0)
c slmod end
c
c     normalize from "per cell" to "per volume"
c
      call eire_renorm(pinenz,kvols,nrs,nks,irsep,sum,1)
c
c-----------------------------------------------------------------------
c
c     average H2 molecule energy
c
c slmod begin
c
c     From here on, changes made to the data input blocks are the same as
c     as for the above.
c
c slmod end
      call rzero (pinenm,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,tdummy(0,0),pinenm,maxnks,maxnrs,1.0,0)
c
c     normalize from "per cell" to "per volume"
c
      call eire_renorm(pinenm,kvols,nrs,nks,irsep,sum,1)
c
c-----------------------------------------------------------------------
c
c     average H2+ molecule energy      !!! (not yet added to DIVIMP) !!!
c
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
c
c     normalize from "per cell" to "per volume"
c
      call eire_renorm(pinbuf,kvols,nrs,nks,irsep,sum,1)
c
c-----------------------------------------------------------------------
c 
      CALL RZero (pinline(1,1,1,H_BALPHA),MAXNKS*MAXNRS*6)

      READ (IPINOUT,*) cdum1
      IF (output) WRITE(0,*) cdum1
c
c     Halpha emissivity (various contributions summed up)
c
c     halpha emission 1 (emitted by H up to ionization)
c
      call rzero (pinalpha,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,tdummy(0,0),pinalpha,maxnks,maxnrs,1.0,0)
      call sum_contrib(pinline(1,1,1,H_BALPHA),pinalpha,nrs,nks,2)
      call sum_contrib(pinline(1,1,6,H_BALPHA),pinalpha,nrs,nks,2)
c
c     halpha emission 2 (emitted in H+ recombination)
c
      call rzero (pinbuf,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
      call sum_contrib(pinalpha,pinbuf,nrs,nks,2)
      call sum_contrib(pinline(1,1,2,H_BALPHA),pinbuf,nrs,nks,2)
      call sum_contrib(pinline(1,1,6,H_BALPHA),pinbuf,nrs,nks,2)
c
c     halpha emission 3 (emitted by H2 up to dissociation)
c
      call rzero (pinbuf,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
      call sum_contrib(pinalpha,pinbuf,nrs,nks,2)
      call sum_contrib(pinline(1,1,3,H_BALPHA),pinbuf,nrs,nks,2)
      call sum_contrib(pinline(1,1,6,H_BALPHA),pinbuf,nrs,nks,2)
c
c     halpha emission 4 (emitted by H2+ up to dissociation)
c
      call rzero (pinbuf,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
      call sum_contrib(pinalpha,pinbuf,nrs,nks,2)
      call sum_contrib(pinline(1,1,4,H_BALPHA),pinbuf,nrs,nks,2)
      call sum_contrib(pinline(1,1,6,H_BALPHA),pinbuf,nrs,nks,2)
c
c     halpha emission 5 (emitted in CX of H and H+)
c
      call rzero (pinbuf,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
      call sum_contrib(pinalpha,pinbuf,nrs,nks,2)
      call sum_contrib(pinline(1,1,5,H_BALPHA),pinbuf,nrs,nks,2)
      call sum_contrib(pinline(1,1,6,H_BALPHA),pinbuf,nrs,nks,2)
c
c...cylinder
      call eire_renorm(pinalpha,kvols,nrs,nks,irsep,sum,0)
      DO i1 = 1, 6
        call eire_renorm(pinline(1,1,i1,H_BALPHA),kvols,nrs,nks,irsep,
     .                   sum,0)
      ENDDO
c
c-----------------------------------------------------------------------
      CALL RZero(pinline(1,1,1,H_BBETA),MAXNKS*MAXNRS*6)

      READ (IPINOUT,*) cdum1
      IF (output) WRITE(0,*) cdum1
c
c     Hgamma emissivity (various contributions summed up)
c
c     Hgamma emission 1 (emitted by H up to ionization)
c
      call rzero (pinbuf,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
      call sum_contrib(pinline(1,1,1,H_BBETA),pinbuf,nrs,nks,2)
      call sum_contrib(pinline(1,1,6,H_BBETA),pinbuf,nrs,nks,2)
c
c     Hgamma emission 2 (emitted in H+ recombination)
c
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
      call sum_contrib(pinline(1,1,2,H_BBETA),pinbuf,nrs,nks,2)
      call sum_contrib(pinline(1,1,6,H_BBETA),pinbuf,nrs,nks,2)
c
c     Hgamma emission 3 (emitted by H2 up to dissociation)
c
      call rzero (pinbuf,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
      call sum_contrib(pinline(1,1,3,H_BBETA),pinbuf,nrs,nks,2)
      call sum_contrib(pinline(1,1,6,H_BBETA),pinbuf,nrs,nks,2)
c
c     Hgamma emission 4 (emitted by H2+ up to dissociation)
c
      call rzero (pinbuf,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
      call sum_contrib(pinline(1,1,4,H_BBETA),pinbuf,nrs,nks,2)
      call sum_contrib(pinline(1,1,6,H_BBETA),pinbuf,nrs,nks,2)
c
c     Hgamma emission 5 (emitted in CX of H and H+)
c
      call rzero (pinbuf,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
      call sum_contrib(pinline(1,1,5,H_BBETA),pinbuf,nrs,nks,2)
      call sum_contrib(pinline(1,1,6,H_BBETA),pinbuf,nrs,nks,2)
c
      DO i1 = 1, 6
        call eire_renorm(pinline(1,1,i1,H_BBETA),kvols,nrs,nks,irsep,
     .                   sum,0)
      ENDDO
c
c-----------------------------------------------------------------------
c
      CALL RZero (pinline(1,1,1,H_BGAMMA),MAXNKS*MAXNRS*6)

      READ (IPINOUT,*) cdum1
      IF (output) WRITE(0,*) cdum1
c
c     Hgamma emissivity (various contributions summed up)
c
c     Hgamma emission 1 (emitted by H up to ionization)
c
      call rzero (pinbuf,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
      call sum_contrib(pinline(1,1,1,H_BGAMMA),pinbuf,nrs,nks,2)
      call sum_contrib(pinline(1,1,6,H_BGAMMA),pinbuf,nrs,nks,2)
c
c     Hgamma emission 2 (emitted in H+ recombination)
c
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
      call sum_contrib(pinline(1,1,2,H_BGAMMA),pinbuf,nrs,nks,2)
      call sum_contrib(pinline(1,1,6,H_BGAMMA),pinbuf,nrs,nks,2)
c
c     Hgamma emission 3 (emitted by H2 up to dissociation)
c
      call rzero (pinbuf,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
      call sum_contrib(pinline(1,1,3,H_BGAMMA),pinbuf,nrs,nks,2)
      call sum_contrib(pinline(1,1,6,H_BGAMMA),pinbuf,nrs,nks,2)
c
c     Hgamma emission 4 (emitted by H2+ up to dissociation)
c
      call rzero (pinbuf,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
      call sum_contrib(pinline(1,1,4,H_BGAMMA),pinbuf,nrs,nks,2)
      call sum_contrib(pinline(1,1,6,H_BGAMMA),pinbuf,nrs,nks,2)
c
c     Hgamma emission 5 (emitted in CX of H and H+)
c
      call rzero (pinbuf,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
      call sum_contrib(pinline(1,1,5,H_BGAMMA),pinbuf,nrs,nks,2)
      call sum_contrib(pinline(1,1,6,H_BGAMMA),pinbuf,nrs,nks,2)
c
      DO i1 = 1, 6
        call eire_renorm(pinline(1,1,i1,H_BGAMMA),kvols,nrs,nks,irsep,
     .                   sum,0)
      ENDDO
c
c slmod end
c-----------------------------------------------------------------------
c
c     H recombination source rate
c
c slmod begin
      READ (IPINOUT,*) cdum1
      IF (output)
     .WRITE(0      ,*) cdum1
c slmod end
      call rzero (pinrec,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,tdummy(0,0),pinrec,maxnks,maxnrs,1.0,0)
c
c     normalize from "per cell" to "per volume"
c
      call eire_renorm(pinrec,kvols,nrs,nks,irsep,sum,0)
c
      write(6,*) 'Total recombinations per second (PIN): ',sum
c
c-----------------------------------------------------------------------
c
c     H+ (and impurities) energy source
c
      call rzero (pinqi,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,tdummy(0,0),pinqi,maxnks,maxnrs,1.0,0)
c
c     normalize from "per cell" to "per volume"
c
      call eire_renorm(pinqi,kvols,nrs,nks,irsep,sum,0)
      write(6,*) 'Total ion energy source (PIN): ',sum
c
c-----------------------------------------------------------------------
c
c     Electron energy source
c
c slmod begin
c      READ (IPINOUT,*) cdum1
c      IF (output)
c     .WRITE(0      ,*) cdum1
c slmod end
      call rzero (pinqe,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              1,maxnfla,tdummy(0,0),pinqe,maxnks,maxnrs,1.0,0)
c
c     normalize from "per cell" to "per volume"
c
      call eire_renorm(pinqe,kvols,nrs,nks,irsep,sum,0)

      write(6,*) 'Total electron energy source (PIN): ',sum
c slmod begin
c
c-----------------------------------------------------------------------
c
c     Recombination particle source
c
      READ (IPINOUT,*) cdum1
      IF (output)
     .WRITE(0      ,*) cdum1

      call rzero   (pinior,maxnks*maxnrs)
      call gfsub3r (ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd (tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     .              1,maxnfla,tdummy(0,0),pinior,maxnks,maxnrs,1.0,0)
      call eire_renorm(pinior,kvols,nrs,nks,irsep,sum,0)

      READ (IPINOUT,*) cdum1
      IF (output)
     .WRITE(0      ,*) cdum1

      call rzero   (pinmpr,maxnks*maxnrs)
      call gfsub3r (ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd (tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     .              1,maxnfla,tdummy(0,0),pinmpr,maxnks,maxnrs,1.0,0)
      call eire_renorm(pinmpr,kvols,nrs,nks,irsep,sum,0)

      READ (IPINOUT,*) cdum1
      IF (output)
     .WRITE(0      ,*) cdum1

      call rzero   (pinqir,maxnks*maxnrs)
      call gfsub3r (ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd (tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     .              1,maxnfla,tdummy(0,0),pinqir,maxnks,maxnrs,1.0,0)
      call eire_renorm(pinqir,kvols,nrs,nks,irsep,sum,0)

      READ (IPINOUT,*) cdum1
      IF (output)
     .WRITE(0      ,*) cdum1

      call rzero   (pinqer,maxnks*maxnrs)
      call gfsub3r (ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
      call deadadd (tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     .              1,maxnfla,tdummy(0,0),pinqer,maxnks,maxnrs,1.0,0)
      call eire_renorm(pinqer,kvols,nrs,nks,irsep,sum,0)

      IF (stopopt3.EQ.8) THEN
        IF (output) WRITE(0,*) 'WARNING: ADDING RECOMINATION POWER TERM'
        DO ir = 1, MAXNRS
          DO ik = 1, MAXNKS
            pinqe(ik,ir) = pinqe(ik,ir) + pinqer(ik,ir)
          ENDDO
        ENDDO
      ENDIF

c...  Individual contributions to the hydrogenic ionisation source (PINION) 
c     for atoms (D), molecules (D2) and test ions (D2+):
      DO i1 = 1, 3
        READ(IPINOUT,'(A72)') cdum1
        IF (output) WRITE(0,*) cdum1(1:LEN_TRIM(cdum1))
        WRITE(PINOUT,*) cdum1(1:LEN_TRIM(cdum1))

        call rzero(pinioncomp(1,1,i1),maxnks*maxnrs)
        call gfsub3r(ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
        call deadadd(tdummy,nx,ny,nxd,nyd)
        call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >                1,maxnfla,tdummy(0,0),pinioncomp(1,1,i1),
     .                maxnks,maxnrs,1.0,0)
        call eire_renorm(pinioncomp(1,1,i1),kvols,nrs,nks,irsep,sum,0)
      ENDDO

c...  Algebraic volume tallies:
      call RZero(pinalgv,maxnks*maxnrs*MAXTALLY)
      DO i1 = 1, eirntally

        READ (IPINOUT,'(A32)') cdum1
        IF (output) WRITE(0,*) cdum1

        call gfsub3r (ipinout,nx,ny,nxd,nyd,1,maxnfla,tdummy(0,0),0)
        call deadadd (tdummy,nx,ny,nxd,nyd)
        call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     .                1,maxnfla,tdummy(0,0),pinalgv(1,1,i1),
     .                maxnks,maxnrs,1.0,0)

        IF (eirtally(i1,5)(1:3).EQ.'vol') THEN
          WRITE(0     ,*) 'RESCALING TALLY ',i1,' BY VOLUME'
          WRITE(PINOUT,*) 'RESCALING TALLY ',i1,' BY VOLUME'
          call eire_renorm(pinalgv(1,1,i1),kvols,nrs,nks,irsep,sum,0)
        ENDIF

c        ir = 30
c        DO ik = 1, nks(ir)
c           WRITE(0,*) '-->',ik,ir,pinion(ik,ir),pinalgv(ik,ir,i1)
c        ENDDO

      ENDDO



c slmod end
c
c-----------------------------------------------------------------------
c
c     Other quantities ... not read in yet.
c
c
c     Neutral particle velocity distribution information from nimbus
c     (pinvdist)
c
c     (not read in from EIRENE yet)
c
c-----------------------------------------------------------------------
c slmod begin
c
c     Wall fluxes:
c
      CALL ReadWallFlux
      CALL ReadPumpedFlux

c...EIRENE only writes this data if BGK collisions are being modeled, which
c   isn't really correct, but is the way things are at the moment.  It will be
c   fixed in the future at some point:
      CALL ReadAdditionalCellData


c...  Time-to ionisation data:
      IF (eirniontime.GT.0) CALL ReadIonisationTime



c
c     BGK related quantites:
c
c
      CALL RZero(pinbgk,MAXNKS*MAXNRS*MAXBGK*MAXTOR)

      READ (IPINOUT,'(A72)',END=10,ERR=98) cdum1
      IF (output) WRITE(0,'(A)') '"'//cdum1(1:10)//'"'

      BACKSPACE IPINOUT

      IF (cdum1(1:10).EQ.'[BGK: ION ') THEN

        IF (output) WRITE(0,*) 'READING BGK DATA'

        IF (nfla.NE.6)
     .    CALL ER('ReadEIRE','PINBGK array not sized correctly. NFLA='//
     .                       '6 required.',*99)

c...5 is the number of BGK related quantities
        DO in = 1, eirnsdtor
          DO i1 = 1, 5
            READ (IPINOUT,'(A72)') cdum1
            IF (output) WRITE(0,*) cdum1(1:LEN_TRIM(cdum1))

            call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla,maxnfla,
     .                   ndummy(0,0,1),0)
            DO i2 = 1, nfla
              call deadadd (ndummy(0,0,i2),nx,ny,nxd,nyd)
              call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     .                      nfla,maxnfla,ndummy(0,0,i2),
     .                      pinbgk(1,1,(i2-1)*5+i1+(in-1)*MAXBGK),
     .                      maxnks,maxnrs,1.0,0)
c              WRITE(0,*) '------------->',pinbgk(78,30,10+2)
c              call eire_renorm(pinbgk(1,1,(i2-1)*5+i1+(in-1)*MAXBGK),
c     .                         kvols,nrs,nks,irsep,sum,1)
c              WRITE(0,*) '------------->',pinbgk(78,30,10+2)
            ENDDO
          ENDDO
        ENDDO

      ENDIF
10    CONTINUE
      IF (output) WRITE(0,*) 'DONE READING BGK DATA'
c
c-----------------------------------------------------------------------
c slmod end
c
c     Close the files that have been used for input
c
c slmod begin



      READ (IPINOUT,'(A72)') cdum1
      IF (output)
     .WRITE(0      ,*) cdum1(1:LEN_TRIM(cdum1))

c... PROBLEMS:
c 	-could be due to algebric tallies being specified in a
c        custom (not written by DIVIMP) input file, and so
c        they are not expected by DIVIMP (see unstructured input
c        *E33)

      IF (cdum1(1:11).NE.'[LAST LINE]') THEN
        WRITE(0,*) 'ERROR: INCORRECT POSITION WHEN READING EIRENE'//
     .             ' TRANSFER FILE'
        WRITE(0,*) 'CDUM1:'//cdum1(1:LEN_TRIM(cdum1))
        STOP
      ENDIF


c...Old method of reading the end of the .eir EIRENE data file.
      IF (output) WRITE(0,*) 'READING ADDITIONAL SURFACE DATA'
      CALL ReadAdditionalSurfaceData
c...This is the new and improved additional cell (BGK related - typically)
c   data transfer routine.


c      IF (eirbgk.GT.0) CALL ReadAdditionalCellData
      IF (output) WRITE(0,*) 'READING PARTICLE TRACKS'

      IF (eirtrc1.GT.0.OR.eirtrc2.GT.0) CALL ReadParticleTracks
      
      IF (output) WRITE(0,*) 'DONE'

c slmod end
      close (ipinout)
c slmod begin
c
c I am contracting the grid so that the added rings/cells are
c removed before DIVIMP does any calculations with the PIN data.
c This will hopefully minimize the chance for errors and the
c number of arrays that have to be shifted by the cell and ring
c deletion routines.
c
c Jul 10, 97 - Ring deletion/addition is no longer being considered.
c
      CALL OutputEIRENE(65,'BACK FROM EIRENE, STRUCTURED GRID')

      CALL OutputData(87,'BACK FROM EIRENE, GRID STRUCTURED')

      IF ((cgridopt.EQ.3.AND.cmodopt.EQ.1).OR.stagopt.GT.1.OR.
     .    stopopt.EQ.121.OR.eirgrid.EQ.1) THEN


        CALL UnstructureGrid

        cutpt1   = ikto + 1
        IF (grdnmod.NE.0) THEN
          cutpt2 = cutpt1 + nks(2) 
        ELSE 
          cutpt2   = ikti + 1
        ENDIF
        cutring  = irsep - 1
        maxkpts  = nks(irsep)
        maxrings = irwall

        CALL OutputEIRENE(66,'GRID JUST UNSTRUCTURED')
c        CALL WriteRaw(91)
      ENDIF

      IF (stopopt.EQ.104) THEN
        CALL DeleteRing(2       )
        CALL DeleteRing(irtrap+1)

        CALL OutputGrid(87,'Rings out')

        vpolyp = vpolmin
      ENDIF


      DO ir = nrs, 1, -1
        IF (idring(ir).EQ.-2) CALL DeleteRing(ir)
      ENDDO
      vpolyp = vpolmin

c      CALL OutputGrid(87,'DONE MESSING WITH RINGS')
c      STOP 'FNDJK'

      CALL OutputStratumData
      CALL OutputLineRadiationData
      CALL OutputBGKData
      CALL OutputRecombinationData
      CALL OutputMomentumData

c...TMP: Calculate radial sources:
      IF     (s28cfpdrft.GT.0) THEN
        WRITE(0,*) 'CALCULATING RADIAL FLUX AFTER RETURNING FROM EIRENE'
        CALL CalcRadialDrift(-1)
      ELSEIF (s28cfpdrft.EQ.-1) THEN
c...    Turning it off:
        WRITE(0,*) 'TURING OFF RADIAL FLUX AFTER RETURNING FROM EIRENE'
        DO ir = 1, nrs
          DO ik = 1, nks(ir)
            osmcfpflx(ik,ir,2) = (1.0 - rel_frac) * osmcfpflx(ik,ir,2) 
          ENDDO
        ENDDO
      ENDIF

      DO ir = 1, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        DO ik = 1, nks(ir)
          IF (virtag(ik,ir).EQ.1) WRITE(SLOUT,*) 'DEBUG: VIRTUAL',ik,ir
          IF (zs(ik,ir).GT.1.5) WRITE(SLOUT,*) 'DEBUG: HIGH',ik,ir
          id = korpg(ik,ir)
          DO in = 1, 4
            IF (zvertp(in,id).GT.1.5) WRITE(SLOUT,*) 
     .        'DEBUG: HIGH P',ik,ir
          ENDDO
        ENDDO
      ENDDO

c      CALL OutputGrid(87,'Just after return from EIRENE')

      IF (stopopt.EQ.17) STOP 'Just after returning from EIRENE'
c slmod end
c
c     Calculate the DIVIMP values for the recombination.
c
c slmod begin - david
      call calc_divrec(totrec)
c
c      call calc_divrec
c slmod end
c
c     Calculate the target fluxes and compare particle balances
c
      call targflux
c
c     Print out some of the Eirene data in DIVIMP format
c
      call eireprn
c
      IF (output) WRITE(0,*) 'DONE READING EIRENE DATA FILE'

c      CALL DumpGrid('JUST BACK FROM EIRENE')

      return
c slmod begin
98    STOP 'FILE ACCESS ERROR IN READEIR'
99    STOP 'READEIR'
c slmod end
      end
c
c ======================================================================
c
c subroutine: ReadIonisationTime
c
c This routine reads the pumped flux data passed from EIRENE in the
c main transfer file and assigns HESCPD and HESCAL.
c
      SUBROUTINE ReadIonisationTime
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER         nradd
      COMMON /RADCOM/ nradd

      INTEGER    MAXNDATA
      PARAMETER (MAXNDATA=2000)

      INTEGER        fp,i1,i2,i3,idum1
      REAL           rdum1
      CHARACTER*1024 cdum1

      fp = IPINOUT

      READ(fp,'(A10)') cdum1
     
      IF (cdum1(1:10).NE.'[TIME-TO-I') 
     .  CALL ER('ReadIonisationTime','Wrong position in EIRENE '//
     .          'transfer file',*99)
 
      READ(fp,'(A72)') cdum1

      DO i1 = 1, eirniontime
        READ(fp,*) idum1,(rdum1,i2=1,4),eiriontime(i1,18),
     .                                  eiriontime(i1,19)
      ENDDO

      DO i1 = 1, eirniontime
        IF (eiriontime(i1,5).EQ.0.0) CYCLE
        READ(fp,'(A72)') cdum1
        DO i2 = 1, NINT(eiriontime(i1,5))+2
          READ(fp,*) idum1,idum1,(eiriontime(i1,20+3*(i2-1)+i3),i3=0,2)
        ENDDO
      ENDDO

      RETURN
99    WRITE(0,*) 'CDUM1:"'//cdum1(1:10)//'"'
      STOP
      END
c
c ======================================================================
c
c subroutine: ReadPumpedFlux
c
c This routine reads the pumped flux data passed from EIRENE in the
c main transfer file and assigns HESCPD and HESCAL.
c
      SUBROUTINE ReadPumpedFlux
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER         nradd
      COMMON /RADCOM/ nradd

      INTEGER    MAXNDATA
      PARAMETER (MAXNDATA=2000)

      INTEGER       fp,ndat,surface,i1
      REAL          area,flux,totflux
      CHARACTER*128 cdum1,particle

      fp = IPINOUT

      hescpd = 0.0
      hescal = 0.0

      WRITE(PINOUT,*)
      WRITE(PINOUT,*) 'EIRENE pumped flux data:'
      WRITE(PINOUT,90) 'NO.','SPECIES','SURFACE AREA (CM**2)',
     .                 'PUMPED FLUX (AMPS)','HESCPD','HESCAL'
90    FORMAT(1X,A6,2X,A10,2A24,2A12)

      READ(fp,'(A128)') cdum1

      IF (cdum1(1:LEN_TRIM(cdum1)).EQ.'[PUMPING SURFACES]') THEN
        READ(fp,*) ndat
        READ(fp,*)
        DO i1 = 1, ndat
          READ(fp,*) surface,particle,area,flux

          IF     (particle(1:8).EQ.'D       ') THEN
c...        The surface area is in cm**2 and the flux is in Amps:
            totflux = flux / ECH 
            hescpd  = hescpd + totflux
c...        The core ring surface number is -1:
            IF (surface.NE.-1) hescal = hescal + totflux
            WRITE(PINOUT,91) surface,particle(1:8),area,flux,hescpd,
     .                       hescal
91          FORMAT(1X,I6,2X,A10,1P,2E24.6,2E12.4,0P)
          ELSEIF (particle(1:8).EQ.'D2      ') THEN
            totflux = flux / ECH * 2.0
            hescpd  = hescpd + totflux
            IF (surface.NE.-1) hescal = hescal + totflux
            WRITE(PINOUT,91) surface,particle(1:8),area,flux,hescpd,
     .                       hescal
          ELSE
            CALL WN('ReadPumpedFluxData','Unknown particle type')
          ENDIF
        ENDDO
      ELSE
        CALL ER('ReadPumpedFluxData','Data not found in the EIRENE '//
     .                               'transfer file',*99)
      ENDIF

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: ReadAdditionalCellData
c
c
      SUBROUTINE ReadAdditionalCellData
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slcom'

      INTEGER         nradd
      COMMON /RADCOM/ nradd

      INTEGER       fp,i1,i2,i3,i4,i5,idum1,idum2,idum3,in,add,
     .              ii(MAXASCDAT),iv(MAXASCDAT),nstats,count,cindex,
     .              nstrai,nmoli,natmi
      LOGICAL       output,status,writeline
      REAL          dv(MAXASCDAT),rdum1
c     .              tmpasd(maxascdat,MAXASD2,MAXASS,1)
c     .              tmpasd(MAXASCDAT,MAXASD2,MAXASS,5)
      CHARACTER*128 cdum1,cdum2,cdum3,cdum4

      REAL, ALLOCATABLE :: tmpasd(:,:,:,:)


      INTEGER nread,init
      REAL    frac,var
     
      DATA   count,cindex,nread,init /0,0,0,0/

      SAVE

      IF (outmode.GE.3)  WRITE(0,*) 'READING ADDITIONAL CELL DATA'

      IF (init.EQ.0) THEN
c        WRITE(0,*) 'ALLOCATING TMPASD'
        ALLOCATE(tmpasd(1+eirnpgdat+asc_ncell*ascncut,MAXASD2,MAXASS,5))
        init = 1
      ENDIF


c      WRITE(0,*) '*** NOT READING ADDITIONAL CELL DATA ***'
c      RETURN

      output = .FALSE.
c      output = .TRUE.

      IF (output) WRITE(0,*) 'READING ADDITIONAL CELL DATA'

      IF (.NOT.thesis.AND.relreset.EQ.1.AND.rel_iter.EQ.1.AND.
     .    rel_opt.GT.0) THEN
        WRITE(0     ,*) 'RESETTING COUNT IN READADDITIONALCELLDATA'
        WRITE(PINOUT,*) 'RESETTING COUNT IN READADDITIONALCELLDATA'
        cindex = 0
        count  = 0
        CALL RZero(tmpasd,(asc_ncell*ascncut)*MAXASD2*MAXASS*5)
      ENDIF

      WRITE(PINOUT,*)
      WRITE(PINOUT,*) 'EIRENE additional cell data:'

c...  Get additional cell tallies from the main EIRENE-DIVIMP transfer
c     file:
      fp = IPINOUT
 
      READ(fp,'(A20)') cdum1 
      IF (output) WRITE(0,*) 'CDUM1:',cdum1(1:20)
      READ(fp,*) idum1

      IF (idum1.NE.1+eirnpgdat+asc_ncell*ascncut) 
     .  CALL ER('ReadAdditionalCellData','Number of additional cells '//
     .                                   'does not match',*99)

c...  Skip over unwanted data (main vacuum cell and pressure gauges):

c...  (TAKE THIS OUT -- BUT HAVE TO MAKE THE NECESSARY INDEXING ADJUSTMENTS IN OUT982 and OUT984)
      DO i1 = 1, 3+eirnpgdat
        READ(fp,*) 
      ENDDO

      DO i1 = 1, asc_ncell*ascncut
        READ(fp,*) idum1,rdum1,ascdata(i1,1),ascdata(i1,2)
      ENDDO




      fp = 98

      nradd = -1

      OPEN(UNIT=fp,FILE='addsur.dat',ACCESS='SEQUENTIAL',
     .     STATUS='OLD',ERR=96)

      status = .TRUE.

      DO WHILE(.TRUE.) 
c...    Read additional cell data for each EIRENE iteration (not
c       DIVIMP/EIRENE iteration -- although typically EIRENE will 
c       only be called once):

        status      = .TRUE.
        cdum1(1:10) = '          '

c...    Advance to the start of the additional cell data for the 
c       next EIRENE self-iteration:
        DO WHILE (cdum1(1:10).NE.'* DATA FOR') 
          READ(fp,'(A128)',END=20,ERR=98) cdum1        
          IF (output) WRITE(0,*) 'CDUM1: ',cdum1(1:10)
        ENDDO

        DO WHILE (status)
          READ(fp,'(A128)',END=10,ERR=98) cdum1

          IF (cdum1(1:10).EQ.'* DATA FOR') THEN
c...        Trigger exit from this loop:
            STATUS = .FALSE.
            BACKSPACE fp
            CYCLE
          ENDIF 

          IF     (cdum1(1:8).EQ.'* NRADD ') THEN
            IF (output) WRITE(0,*) 'FOUND NRADD'
            READ(fp,*) nradd
          ELSEIF (cdum1(1:8).EQ.'* SURFAC') THEN
            IF (output) WRITE(0,*) 'FOUND SURFACE DATA'
            i1 = 1
            READ(fp,*,END=97,ERR=98)
            DO i2 = 1, nradd
              READ(fp,*,END=97,ERR=98)
     .          idum2,idum3,pinasd(i2,5,i1,1)
              IF (i2.EQ.1) indasd = idum3
            ENDDO
          ELSEIF (cdum1(1:8).EQ.'* NATMI ') THEN
            IF (output) WRITE(0,*) 'FOUND NATMI'
            READ(fp,*) idum1
            IF (idum1.NE.2) CALL ER('ReadAdditionalCellData',
     .                              'Expecting 2 neutral atom species'//
     .                              ' (D and C)',*99)
            DO i1 = 1, 2
              READ(fp,*,END=97,ERR=98)
              DO i2 = 1, nradd
                READ(fp,*,END=97,ERR=98)
     .            idum2,idum3,(pinasd(i2,i3,i1,1),i3=1,2)
                IF (i2.EQ.1) indasd = idum3
              ENDDO
            ENDDO
          ELSEIF (cdum1(1:8).EQ.'* NMOLI ') THEN
            IF (output) WRITE(0,*) 'FOUND NMOLI'
            READ(fp,*) idum1
            IF (idum1.NE.1) CALL ER('ReadAdditionalCellData',
     .                              'Expecting 1 neutral molecule '//
     .                              'species (D2)',*99)
            i1 = 3
            READ(fp,*,END=97,ERR=98)
            DO i2 = 1, nradd
              READ(fp,*,END=97,ERR=98)
     .          idum2,idum3,(pinasd(i2,i3,i1,1),i3=1,2)
              IF (i2.EQ.1) indasd = idum3
            ENDDO
          ELSEIF (cdum1(1:8).EQ.'* NIONI ') THEN
            IF (output) WRITE(0,*) 'FOUND NIONI'
            READ(fp,*) idum1
            IF (idum1.NE.1) CALL ER('ReadAdditionalCellData',
     .                              'Expecting 1 test ion species '//
     .                              '(D2+)',*99)
            i1 = 4
            READ(fp,*,END=97,ERR=98)
            DO i2 = 1, nradd
              READ(fp,*,END=97,ERR=98)
     .          idum2,idum3,(pinasd(i2,i3,i1,1),i3=1,2)
              IF (i2.EQ.1) indasd = idum3
            ENDDO
          ELSEIF (cdum1(1:8).EQ.'* NPLSI ') THEN
            IF (output) WRITE(0,*) 'FOUND NPLSI (BGK)'
            READ(fp,*) idum1
            IF (idum1.NE.2.AND.idum1.NE.6)
     .        CALL ER('ReadAdditionalCellData','Expecting 2 or 6 '//
     .                'plasma species (D+,C+ or D+,C+,DD,D2D2,DD2'//
     .                ',D2D)',*99)
            DO i1 = 1+4, idum1+4
              READ(fp,*,END=97,ERR=98)
              DO i2 = 1, nradd
                READ(fp,*,END=97,ERR=98)     
     .            idum2,idum3,(pinasd(i2,i3,i1,1),i3=1,5)
                IF (i2.EQ.1) indasd = idum3
              ENDDO
            ENDDO
          ELSEIF (cdum1(1:8).EQ.'BGK RESI') THEN
c...        Storing MODBGK routine statistics:
            eirnres = eirnres + 1
            READ(fp,*) nstats
            IF (nstats.NE.6) 
     .        CALL ER('ReadAdditionalCellData','Expecting 6 n-n '//
     .                'species',*99)
            DO i1 = 1, nstats
              READ(fp,*) idum1,(eirres(i2,i1,eirnres),i2=1,7)
            ENDDO
            WRITE(PINOUT,*)
            WRITE(PINOUT,*) 'BGK STATS:'
            DO i1 = 1, nstats
              WRITE(PINOUT,'(2I4,1P,7E12.4,0P)') i1,eirnres,
     .                        (eirres(i2,i1,eirnres),i2=1,7)
            ENDDO

          ELSEIF (cdum1(1:8).EQ.'STRATUM ') THEN
c...        Storing stratum data for selected additional cells:
            READ(fp,*) eirnstrdat,eirnstrai,eirnatmi,eirnmoli
            DO i1 = 1, eirnstrdat
              READ(fp,*) eirstrlis(i1)
              DO i2 = 1, eirnstrai+1
                DO i3 = 1, eirnatmi
                  READ(fp,*) i4,i5,eirstrdat(i1,i2,2*i3-1+1),
     .                             eirstrdat(i1,i2,2*i3  +1)
                ENDDO
                DO i3 = 1, eirnmoli
                  READ(fp,*) i4,i5,eirstrdat(i1,i2,2*i3-1+2*eirnatmi+1),
     .                             eirstrdat(i1,i2,2*i3  +2*eirnatmi+1)
                ENDDO
              ENDDO
            ENDDO
            READ(fp,*) (eirfluxt(i1),i1=1,eirnstrai)
 

          ENDIF
        ENDDO

10      CONTINUE



c...    Blank data on tiny cells:
        DO i1 = 1, asc_ncell*ascncut
          IF (asc_area(i1).LT.1.0E-05) THEN
c            IF (i1.LE.asc_ncell) 
c     .        WRITE(0     ,*) 'TINY CELL, BLANKING DATA',i1
            WRITE(PINOUT,*) 'TINY CELL, BLANKING DATA',i1
            i2 = i1 + 1 + eirnpgdat
            DO i3 = 1, 10
              pinasd(i2,1,i3,1) = 1.0E+01
              pinasd(i2,2,i3,1) = 1.0E-02
              pinasd(i2,3,i3,1) = 0.0
              pinasd(i2,4,i3,1) = 0.0
              pinasd(i2,5,i3,1) = 0.0
            ENDDO
          ENDIF
        ENDDO


        IF (.TRUE.) THEN
c...      BARF! I am putting this here to accommodate the main version:        
          count  = count  + 1
          nread  = nread  + 1
          cindex = cindex + 1
          frac   = 1 / REAL(MIN(5,count))

c...      Copy data from PINASD(...,1) to TMPASD(...,count):
          IF (cindex.GT.5) cindex = 1
          DO i1 = 1, 1+eirnpgdat+asc_ncell*ascncut
            DO i2 = 1, MAXASD2
              DO i3 = 1, MAXASS         
                tmpasd(i1,i2,i3,cindex) = pinasd(i1,i2,i3,1)  
              ENDDO
            ENDDO
          ENDDO  

c...      Average:
          WRITE(PINOUT,*)
          WRITE(PINOUT,*) 'AVERAGING ACD (PINASD) ',cindex,count,frac
          DO i1 = 1, 1+eirnpgdat+asc_ncell*ascncut
            DO i2 = 1, MAXASD2
              DO i3 = 1, MAXASS
                pinasd(i1,i2,i3,2) = 0.0
                DO i4 = 1, MIN(count,5)
                  pinasd(i1,i2,i3,2) =      pinasd(i1,i2,i3,2 )+
     .                                 frac*tmpasd(i1,i2,i3,i4)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

        IF (nradd.EQ.-1) CALL ER('ReadAdditionalCellData','NRADD not '//
     .                           'set',*99)

        IF (output) WRITE(0,*) 'MARK: NRADD= ',nradd

        WRITE(79,'(A,1X,7I6)') '''ACD 01    1.01''',
     .    rel_step,rel_iter,rel_count,nread,indasd,nradd,MAXASS
        DO i1 = 1, MAXASS
          DO i2 = 1, nradd
            writeline = .FALSE.
            in        = i2 - 1 - eirnpgdat
            DO i4 = 1, eirnaout
              add = MAX(0,NINT(eiraout(i4,2)-1.0) * asc_ncell)
              IF ((eiraout(i4,1).EQ.1.AND.
     .             (eiraout(i4,3)+add.EQ.in.OR.
     .              eiraout(i4,4)+add.EQ.in.OR.
     .              eiraout(i4,5)+add.EQ.in.OR.
     .              eiraout(i4,6)+add.EQ.in)).OR.
     .            ((eiraout(i4,1).EQ.2.OR.eiraout(i4,1).EQ.3).AND.
     .             (eiraout(i4,3)+add.LE.in.AND.
     .              eiraout(i4,4)+add.GE.in.OR.
     .              eiraout(i4,5)+add.LE.in.AND.
     .              eiraout(i4,6)+add.GE.in)))
     .         writeline = .TRUE.
            ENDDO
            IF (writeline) 
     .        WRITE(79,90)  i1,i2,(pinasd(i2,i3,i1,1),i3=1,MAXASD2)
          ENDDO
        ENDDO

        WRITE(79,'(A,1X,7I6)') '''ACD 03    1.01''',
     .    rel_step,rel_iter,rel_count,nread,indasd,nradd,MAXASS
        DO i1 = 1, MAXASS
          DO i2 = 1, nradd
            writeline = .FALSE.
            in        = i2 - 1 - eirnpgdat
            DO i4 = 1, eirnaout
              add = MAX(0,NINT(eiraout(i4,2)-1.0) * asc_ncell)
              IF ((eiraout(i4,1).EQ.1.AND.
     .             (eiraout(i4,3)+add.EQ.in.OR.
     .              eiraout(i4,4)+add.EQ.in.OR.
     .              eiraout(i4,5)+add.EQ.in.OR.
     .              eiraout(i4,6)+add.EQ.in)).OR.
     .            (eiraout(i4,1).EQ.2.AND.
     .             (eiraout(i4,3)+add.LE.in.AND.
     .              eiraout(i4,4)+add.GE.in.OR.
     .              eiraout(i4,5)+add.LE.in.AND.
     .              eiraout(i4,6)+add.GE.in)))
     .         writeline = .TRUE.
            ENDDO
            IF (writeline) 
     .        WRITE(79,90)  i1,i2,(pinasd(i2,i3,i1,2),i3=1,MAXASD2)
          ENDDO
        ENDDO

90      FORMAT(2I6,2X,1P,99(E12.4:))

      ENDDO

20    CONTINUE

      CLOSE (fp)

      IF (.NOT.thesis.AND.rel_opt.GE.1.AND.rel_iter.EQ.rel_niter) THEN
c...    Generate some output for this iteration.  The conditions on 
c       this IF statement should be the same as for the block 
c       that writes the iter-vac.dat file in module output.d6a:
        fp = 97
        OPEN(UNIT=fp,FILE='iter-vac.dat',ACCESS='SEQUENTIAL',
     .       STATUS='UNKNOWN',POSITION='APPEND')
        CALL OutputVacuumCellTable(fp)
        CLOSE (fp)
      ENDIF


      IF (outmode.GE.3)  WRITE(0,*) 'DONE'

c...TEMP
c      DO i1 = 7, 7
c        DO i2 = 1, 10
c          WRITE(0,*)  i1,i2,(pinasd(i2,i3,i1,1),i3=1,2)
c        ENDDO
c      ENDDO


      RETURN
96    CALL ER('ReadAdditionalCellData','Unable to open file',*99)
97    CALL ER('ReadAdditionalCellData','EOF                ',*99)
98    CALL ER('ReadAdditionalCellData','Unable to read file',*99)
99    STOP
      END
c
c ======================================================================
c
c subroutine: ReadAdditionalSurfaceData
c
c
      SUBROUTINE ReadAdditionalSurfaceData
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slcom'

      INTEGER       fp,i1,i2,idum1,ii(MAXASCDAT),iv(MAXASCDAT)
      REAL          dv(MAXASCDAT)
      CHARACTER*128 cdum1,cdum2,cdum3,cdum4

      REAL   count,frac
      DATA   count /0.0/

      SAVE

c     ******************
c     TO BECOME OBSOLETE
c     ******************

c
c     No additional surfaces registered:
c
c      WRITE(0,*) 'NO READING ADDITIONAL SURFACE DATA'
c      RETURN

c      IF (eirnpgdat.EQ.0) RETURN

      IF (.NOT.thesis.AND.relreset.EQ.1.AND.rel_iter.EQ.1.AND.
     .    rel_opt.GT.0) THEN
        WRITE(0     ,*) 'RESETTING COUNT IN READADDITIONALSURFACEDATA'
        WRITE(PINOUT,*) 'RESETTING COUNT IN READADDITIONALSURFACEDATA'
        count = 0.0
      ENDIF

      IF (thesis) THEN
        frac = rel_frac
      ELSE
        count = count + 1.0
        frac  = 1.0 / count               
      ENDIF

      fp = PINOUT

      WRITE(fp,*)
      WRITE(fp,*) 'EIRENE pressure gauge data (EIRPGDAT):'

      DO WHILE (.TRUE.)
        READ(IPINOUT,*,END=10,ERR=98) cdum1,cdum2,cdum3,cdum4

c        WRITE(fp,*) '>',cdum1(1:LEN_TRIM(cdum1)),'<'
c        WRITE(fp,*) '>',cdum2(1:LEN_TRIM(cdum2)),'<'
c        WRITE(fp,*) '>',cdum3(1:LEN_TRIM(cdum3)),'<'
c        WRITE(fp,*) '>',cdum4(1:LEN_TRIM(cdum4)),'<'

        IF (cdum1(1:6).EQ.'TALLY:') THEN

          READ(IPINOUT,*,END=97,ERR=98) idum1
          IF (idum1.GT.MAXASCDAT)
     .      CALL ER('ReadAdditionalSurfaceData','Array bounds'//
     .              ' exceeded',*99)

          DO i1 = 1, idum1
            READ(IPINOUT,*,END=97,ERR=98) ii(i1),dv(i1)
          ENDDO

          IF     (cdum2(1:26).EQ.'ENERGY DENSITY (MOLECULES)'.AND.
     .            cdum3(1: 2).EQ.'D2'                        .AND.
     .            cdum4(1: 9).EQ.'EV*CM**-3') THEN

            WRITE(fp,'(3X,A)') cdum2(1:LEN_TRIM(cdum2))
            WRITE(fp,'(3X,A4,5A7,5A13)')
     .        'Code','x','y','T(deg)','r (m)','L (m)',
     .        'Pexp (mTorr)','Vol (m3)','P (eV cm-3)','P (mTorr)',
     .        'Prel (mTorr)'

            DO i1 = 1, eirnpgdat
              eirpgdat(i1, 9) = dv(i1+1)
              eirpgdat(i1,10) = dv(i1+1) * 1.602E-19 * 1.0E+06 * 0.67
              eirpgdat(i1,11) =        frac  * eirpgdat(i1,10) +
     .                          (1.0 - frac) * eirpgdat(i1,11)

              WRITE(fp,'(3X,I4,5F7.2,1P,5E13.3,0P)')
     .          INT(eirpgdat(i1,1)),(eirpgdat(i1,i2),i2=2,11)
            ENDDO
          ELSEIF (cdum2(1:26).EQ.'ENERGY DENSITY (ATOMS)'.AND.
     .            cdum3(1: 2).EQ.'D'                        .AND.
     .            cdum4(1: 9).EQ.'EV*CM**-3') THEN

            WRITE(fp,'(3X,A)') cdum2(1:LEN_TRIM(cdum2))
            WRITE(fp,'(3X,A4,5A7,5A13)')
     .        'Code','x','y','T(deg)','r (m)','L (m)',
     .        'Pexp (mTorr)','Vol (m3)','P (eV cm-3)','P (mTorr)',
     .        'Prel (mTorr)'

            DO i1 = 1, eirnpgdat
              eirpgdat(i1,12) = dv(i1+1)
              eirpgdat(i1,13) = dv(i1+1) * 1.602E-19 * 1.0E+06 * 0.67
              eirpgdat(i1,14) =        frac  * eirpgdat(i1,13) +
     .                          (1.0 - frac) * eirpgdat(i1,14)

              WRITE(fp,'(3X,I4,5F7.2,1P,5E13.3,0P)')
     .          INT(eirpgdat(i1,1)),(eirpgdat(i1,i2),i2= 2, 8),
     .                              (eirpgdat(i1,i2),i2=12,14)
            ENDDO

          ENDIF

        ELSEIF (cdum1(1:12).EQ.'STRATA DATA:') THEN
c...
          IF (cdum2(1:16).EQ.'NUMBER OF TRACKS') THEN
            READ(IPINOUT,*,END=97,ERR=98) idum1
            DO i1 = 1, idum1
              READ(IPINOUT,*,END=97,ERR=98) ii(i1),iv(i1)
            ENDDO

            IF (eirnstrata+eirnpuff.NE.idum1) 
     .        CALL ER('ReadAdditionalCellData','Number of strata '//
     .                'reported is not consistent',*99)

            DO i1 = 1, eirnstrata+eirnpuff
              eirntracks(i1) = iv(i1)
            ENDDO

            WRITE(fp,*)
            WRITE(fp,*) 'Data for strata:'
            WRITE(fp,*) '   '//cdum2(1:LEN_TRIM(cdum2))
            DO i1 = 1, idum1
              WRITE(fp,'(5X,2I7)') ii(i1),iv(i1)
            ENDDO
          ENDIF

        ENDIF
      ENDDO

10    CONTINUE

      DO i1 = 1, eirnpgdat
        WRITE(79,'(A,3I4,1X,I4,7E11.3)') '''PG-DATA   1.01''',
     .    rel_step,rel_iter,rel_count,i1,(eirpgdat(i1,i2),i2=1, 7)
        WRITE(79,'(A,12X,1X,4X,7E11.3)') '''              ''',
     .                                   (eirpgdat(i1,i2),i2=8,14)
      ENDDO


      RETURN
97    CALL ER('ReadAdditionalSurfaceData','EOF                ',*99)
98    CALL ER('ReadAdditionalSurfaceData','Unable to read file',*99)
99    STOP
      END
c
c ======================================================================
c
c subroutine: ReadParticleTracks
c
c
      SUBROUTINE ReadParticleTracks
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER fp,i1,idum1,idum2,count
      REAL    rdum1,rdum2,rdum3

      WRITE(PINOUT,*)
      WRITE(PINOUT,'(A)') ' Reading EIRENE particle tracks:'

      hwalks(1,1) = 999.0
      hwalks(1,2) = 999.0

      fp = 98

      OPEN(UNIT=fp,FILE='eirtrac.dat',ACCESS='SEQUENTIAL',
     .     STATUS='OLD',ERR=98)

      count = 0
      i1    = 0
      DO WHILE (.TRUE..AND.i1.LT.MAXNWS-2)
        READ(fp,*,END=10,ERR=97) idum1,rdum1,rdum2,rdum3,idum2

        IF (idum2.EQ.0.AND.i1.GT.0) THEN
          count = count + 1

          i1 = i1 + 1
          hwalks(i1,1) = HI
          hwalks(i1,2) = HI
        ENDIF
          
        i1 = i1 + 1
        hwalks(i1,1) = rdum1 * 0.01
        hwalks(i1,2) = rdum2 * 0.01 
      ENDDO
10    CONTINUE

      i1 = i1 + 1
      hwalks(i1,1) = 999.0
      hwalks(i1,2) = 999.0

      CLOSE(fp)

      WRITE(PINOUT,'(A,I8)') '   NO. OF TRACKS= ',count
c      WRITE(0     ,'(A,I8)') '   NO. OF TRACKS= ',count      

      RETURN
96    CALL ER('ReadParticleTracks','EOF'     ,*99)
97    CALL ER('ReadParticleTracks','Problems',*99)
c...  The data file could be missing because either the EIRENE run script
c     needs to be updated, or there were not enough particle tracks followed
c     in EIRENE to reach the track number range set for output in the 
c     EIRENE input file (this is set in block 11 in the template file at the
c     moment):
98    CALL WN('ReadParticleTracks','Unable to open data file')
      RETURN
99    STOP
      END
c
c ======================================================================
c
c DIVIMP writing EIRENE input file:
c
c
c ======================================================================
c
c
c
c
c ======================================================================
c
c subroutine: WriteInputFile
c
c
c BuildNeutralWall needs to be run before this routine can be called...
c since rves and zves have to be set up properly...
c
c
c
c     2
c*
c     2     2
cVL  2.0000E+00
c*
c     3     3
cVL  3.0000E+00
c
      SUBROUTINE WriteInputFile
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER    MAXSTRAT,MAXSDATA
      PARAMETER (MAXSTRAT=10,MAXSDATA=10)
      COMMON /EIRCOM/
     .        eir_07nsrc,eir_07opt ,eir_07stra,
     .        eir_07ind1,
     .        eir_07ind2,eir_07ind3,
     .        eir_07wght           ,eir_07data
      INTEGER eir_07nsrc,eir_07opt ,eir_07stra(MAXSTRAT),
     .        eir_07ind1(MAXSTRAT),
     .        eir_07ind2(MAXSTRAT),eir_07ind3(MAXSTRAT)
      REAL    eir_07wght           ,eir_07data(MAXSTRAT,MAXSDATA)

      INTEGER   ik,ik1,ik2,ir,i1,i2,i3,fp1,fp2,in,icnt,
     .          add,ilst(1024)
      LOGICAL   output,firstcall
      REAL      x0,y0,r,vcel(MAXASCDAT),zaa,roa,fact
      CHARACTER buffer*200,geostr*4

      DATA firstcall /.TRUE./
      SAVE
c
c     Check whether DIVIMP input option requests EIRENE data file:
c
      IF (eirdata.NE.1) RETURN

      output = .FALSE.

      IF (output) WRITE(0,*) 'WRITING EIRENE INPUT FILE'
c
c     Initialization:
c
      fp1   = EIRIN
      fp2   = EIROUT

      OPEN(UNIT=fp1,FORM='FORMATTED',ERR=95,STATUS='OLD')
      OPEN(UNIT=fp2,FORM='FORMATTED',ERR=95,STATUS='REPLACE')

c      fp2 = 0
c      REWIND(fp1)

      CALL MS('WriteInputFile','Using xVESM to store wall data')

c      IF (iflexopt(6).EQ.11) THEN
        eirtemp1 = -ctargt * 1.38E-23 / ECH
        eirtemp2 = -cwallt * 1.38E-23 / ECH
c      ELSE
c        eirtemp1 = ctargt * 1.38E-23 / ECH
c        eirtemp2 = cwallt * 1.38E-23 / ECH
c      ENDIF


c...  Correct simulated pressure gauge volumes for cyclindrical approximation:
      IF (eirnpgdat.GT.0) THEN
        IF     (eirzaa.EQ.-1.0) THEN
          WRITE(PINOUT,*) 'USING CYLINDIRCAL PRESSURE GAUGE VOLUME'
          DO i1 = 1, eirnpgdat
            eirpgdat(i1,8) = PI * eirpgdat(i1,5)**2.0 * 2.0 * PI * rxp * 
     .                       eirtorfrac
          ENDDO
        ELSEIF (eirzaa.LT.0.0.OR.iflexopt(7).EQ.21) THEN
          WRITE(PINOUT,*) 'USING TOROIDAL PRESSURE GAUGE VOLUME'
          DO i1 = 1, eirnpgdat
            eirpgdat(i1,8) = PI * eirpgdat(i1,5)**2.0 * 2.0 * PI *
     .                       eirpgdat(i1,2) * eirtorfrac
          ENDDO
        ELSE
          WRITE(PINOUT,*) 'USING *SPECIFIED* CYLINDIRCAL PRESSURE '//
     .                    'GAUGE VOLUME'
          DO i1 = 1, eirnpgdat
            eirpgdat(i1,8) = PI * eirpgdat(i1,5)**2.0 * 2.0 * PI * 
     .                       eirzaa * eirtorfrac
          ENDDO
        ENDIF
      ENDIF


c
c     This is a somewhat convoluted loop at the moment, which reads
c     through an existing EIRENE input data file that serves as a
c     template for the new data file being written.  Grid and neutral
c     wall specific data are substituted into the template, as well
c     as any EIRENE options/settings that are specifiable
c     from DIVIMP (such as EIRENE execution time):
c

      IF (.NOT.thesis) THEN
c...    Initial strata setup for DIII-D cases.  There are 3 strata:
c       inner target, outer target and volume recombination.
        IF (eirbgk.GE.1) THEN
          eir_07opt = 4
        ELSE
          eir_07opt = 0
        ENDIF

        eir_07wght = eiralloc

        IF (eirnstrata.EQ.0) THEN
          eirnstrata = 3

          eirstrata(1,1) =  1.0
          eirstrata(1,2) = -3.0
          eirstrata(1,3) = -3.0
          eirstrata(1,4) =  999999.0
          eirstrata(1,5) = -1.0

          eirstrata(2,1) =  2.0
          eirstrata(2,2) = -3.0
          eirstrata(2,3) = -3.0
          eirstrata(2,4) =  999999.0
          eirstrata(2,5) = -1.0

          eirstrata(3,1) =  3.0
          eirstrata(3,2) = -3.0
          eirstrata(3,3) = -3.0
          eirstrata(3,4) =  999999.0
          eirstrata(3,5) = -1.0
        ENDIF

        eir_07nsrc = eirnstrata
 
        DO i1 = 1, eirnstrata

          eir_07stra(i1) = i1

          IF     (eirstrata(i1,1).EQ.1.0.OR.
     .            eirstrata(i1,1).EQ.2.0) THEN

            IF (eirstrata(i1,1).EQ.1.0) eir_07ind1(i1) = 1
            IF (eirstrata(i1,1).EQ.2.0) eir_07ind1(i1) = nks(irsep) + 3

            IF     (eirstrata(i1,2).EQ.-1.0) THEN
              eir_07ind2(i1) = irsep  - 1
              eir_07ind3(i1) = irwall - 1
            ELSEIF (eirstrata(i1,2).EQ.-2.0) THEN
              eir_07ind2(i1) = 1
              eir_07ind3(i1) = irsep - 1
            ELSEIF (eirstrata(i1,2).EQ.-3.0) THEN
              eir_07ind2(i1) = 1
              eir_07ind3(i1) = irwall - 1
            ELSEIF (eirstrata(i1,2).GE.irsep.AND.
     .              eirstrata(i1,2).LE.nrs  .AND.
     .              eirstrata(i1,3).GE.irsep.AND.
     .              eirstrata(i1,3).LE.nrs  .AND.
     .              .TRUE.) THEN
c     .              eirstrata(i1,2).LE.eirstrata(i1,3)) THEN
              IF (NINT(eirstrata(i1,2)).GT.irwall-1) THEN
                eir_07ind2(i1) = NINT(eirstrata(i1,2)) - irtrap
              ELSE
                eir_07ind2(i1) = NINT(eirstrata(i1,2)) - 1
              ENDIF
              IF (NINT(eirstrata(i1,3)).GT.irwall-1) THEN
                eir_07ind3(i1) = NINT(eirstrata(i1,3)) - irtrap + 1
              ELSE
                eir_07ind3(i1) = NINT(eirstrata(i1,3)) 
              ENDIF
            ELSE
              CALL ER('WriteInputFile','Invalid stratum region',*99)
            ENDIF

          ELSEIF (eirstrata(i1,1).EQ.3.0) THEN
c...        Volume recombination source:
            IF     (eirstrata(i1,3).EQ.-3.0) THEN
              eir_07ind1(i1) = 1
              eir_07ind2(i1) = irwall - 2
              eir_07ind3(i1) = 0

            ELSEIF (eirstrata(i1,2).GE.irsep.AND.
     .              eirstrata(i1,2).LE.nrs  .AND.
     .              eirstrata(i1,3).GE.irsep.AND.
     .              eirstrata(i1,3).LE.nrs) THEN
c     .              eirstrata(i1,2).LE.eirstrata(i1,3)) THEN

              IF (NINT(eirstrata(i1,2)).GT.irwall-1) THEN
                eir_07ind1(i1) = NINT(eirstrata(i1,2)) - irtrap
              ELSE
                eir_07ind1(i1) = NINT(eirstrata(i1,2)) - 1
              ENDIF
              IF (NINT(eirstrata(i1,3)).GT.irwall-1) THEN
                eir_07ind2(i1) = NINT(eirstrata(i1,3)) - irtrap 
              ELSE
                eir_07ind2(i1) = NINT(eirstrata(i1,3)) 
              ENDIF

              eir_07ind3(i1) = 0



c              eir_07ind1(i1) = 1
c              eir_07ind2(i1) = irwall - 2
c              eir_07ind3(i1) = 0
cc...          Turn off recombination via the opacity multiplier:
c              DO ir = 1, nrs
c                IF (ir.GE.NINT(eirstrata(i1,2)).AND.
c     .              ir.LE.NINT(eirstrata(i1,3))) THEN
c                  DO ik = 1, nks(ir)
c                    IF (rs(ik,ir).GT.0.58) THEN
c                      mulrec(ik,ir) = 1.0
c                    ELSE
c                      mulrec(ik,ir) = 0.0
c                    ENDIF
c                  ENDDO
c                ELSE
c                  DO ik = 1, nks(ir)
c                    mulrec(ik,ir) = 0.0
c                  ENDDO
c                ENDIF
c              ENDDO

            ELSE
              CALL ER('WriteInputFile','Invalid stratum region',*99)
            ENDIF

          ELSE
            CALL ER('WriteInputFile','Invalid stratum type',*99)
          ENDIF

          eir_07data(i1,1) = eirstrata(i1,4)
        ENDDO




c        eir_07stra(1) = 1
c        eir_07ind1(1) = 1
c        eir_07ind2(1) = 1
c        eir_07ind3(1) = irwall - 1
c        eir_07stra(2) = 2
c        eir_07ind1(2) = nks(irsep) + 3
c        eir_07ind2(2) = 1
c        eir_07ind3(2) = irwall - 1
c        eir_07stra(3) = 3
c        eir_07ind1(3) = 1
c        eir_07ind2(3) = 1
c        eir_07ind3(3) = 1
c        eir_07ind1(4) = 1
c        eir_07ind2(4) = 1
c        eir_07ind3(4) = 1
c...    NPTS (particle track limit in EIRENE):
c        DO i1 = 1, eir_07nsrc
c          eir_07data(i1,1) = 999999.0
c        ENDDO






      ELSE
        IF (nrs.EQ.44) THEN
c          eir_07opt = 1
          eir_07opt = 2
          IF (stopopt.EQ.107.OR.stopopt.EQ.113.OR.stopopt.EQ.114.OR.
     .        eirbgk.GE.1)
     .      eir_07opt = 5

          IF (eir_07opt.EQ.1) THEN
            eir_07nsrc = 3
            eir_07wght = 1.0
            eir_07stra(1) = 1
            eir_07stra(2) = 2
            eir_07stra(3) = 3
            eir_07ind1(1) = 1
            eir_07ind2(1) = 1
            eir_07ind3(1) = irwall - 1
            eir_07ind1(2) = nks(irsep) + 3
            eir_07ind2(2) = 1
            eir_07ind3(2) = irwall - 1
            eir_07ind1(3) = 1
            eir_07ind2(3) = 1
            eir_07ind3(3) = 1

c...        NPTS   (particle track limit in EIRENE)
            eir_07data(1,1) = 100000.0
            eir_07data(2,1) = 100000.0
            eir_07data(3,1) = 100000.0
          ELSE
            eir_07nsrc    = 6
c...        ALLOC (weighting between NPTS*(1.0-ALLOC) and FLUX*ALLOC):
            eir_07wght    = 0.25

            eir_07stra(1) = 1
            eir_07stra(2) = 1
            eir_07stra(3) = 2
            eir_07stra(4) = 2
            eir_07stra(5) = 3
            eir_07stra(6) = 4

            eir_07ind1(1) = 1
            eir_07ind2(1) = 1
            eir_07ind3(1) = irbreak - 1
            eir_07ind1(2) = 1
            eir_07ind2(2) = irbreak - 1
            eir_07ind3(2) = irwall  - 1

            eir_07ind1(3) = nks(irsep) + 3
            eir_07ind2(3) = 1
            eir_07ind3(3) = 21
            eir_07ind1(4) = nks(irsep) + 3
            eir_07ind2(4) = 21
            eir_07ind3(4) = irwall - 1

            eir_07ind1(5) = 1
            eir_07ind2(5) = (ikti + ikto) / 2
            eir_07ind3(5) = 1

            eir_07ind1(6) = (ikti +ikto) / 2 + 1
            eir_07ind2(6) = nks(irsep) + 3
            eir_07ind3(6) = 1

c...        NPTS (particle track limit in EIRENE)
            eir_07data(1,1) = 500000.0
            eir_07data(2,1) = 010000.0
            eir_07data(3,1) = 500000.0
            eir_07data(4,1) = 005000.0
            eir_07data(5,1) = 900000.0
            eir_07data(6,1) = 900000.0
          ENDIF
        ELSE
          eir_07opt = 0
          IF (stopopt.EQ.106.OR.stopopt.EQ.110.OR.stopopt.EQ.111.OR.
     .        eirbgk.GE.1)
     .      eir_07opt = 4

          IF (eirpmode.GE.1) THEN
            eir_07nsrc  = 3 + eirnpuff
          ELSE
            eir_07nsrc  = 3
          ENDIF   

          eirnstrata = eir_07nsrc

        ENDIF
      ENDIF




10    CONTINUE

      CALL ReadLine(fp1,buffer,1,*50,*98)

20    CONTINUE

      IF (buffer(1:6).EQ.'*** 0.') THEN
c
c Need to remove the requirement that the template file have an intitial
c seciton labelled *** 0...
c
c       This section has been added to EIRENE and contains options
c       for the new EIRENE code that are related to the
c       generalization of the grid:
c
        WRITE(fp2,'(A)') '*** 0. DIVIMP RELATED SETUP DATA (DIVIMP)'
        WRITE(fp2,'(2A,I6)') '''Geometry option  (GEOMOPT)  ',
     .    '0-standard   1-from DIVIMP    ''',eirgeom
        WRITE(fp2,'(2A,I6)') '''Grid option      (GRIDOPT)  ',
     .    '0-structured 1-generalized    ''',eirgrid
c        WRITE(fp2,'(2A,I6)') '''AddUsr option    (ADDOPT)   ',
c     .    '0-execute    1-do not execute ''',eiradd
        WRITE(fp2,'(2A,I6)') '''Wall data option (NEUTOPT)  ',
     .    '0-standard   1-accurate       ''',eirneut
        WRITE(fp2,'(2A,I6)') '''Debug option     (DEBUGOPT) ',
     .    '0-off                         ''',eirdebug
        WRITE(fp2,'(2A,I6)') '''CX D2+ production(CXD2OPT)  ',
     .    '0-off 1-Dalpha only 2-full    ''',eircxd2
        WRITE(fp2,'(2A,I6)') '''User ID          (OPTUSER)  ',
     .    '                              ''',optuser
        WRITE(fp2,'( A,I6)') '''Transparent toroidal spans of non-'//
     .                      'standard surfaces       ''',eirntrans
        IF (eirntorseg.NE.0) THEN
          fact = -360.0 / eirzaa
        ELSE
          fact = 100.0
        ENDIF
        DO i1 = 1, eirntrans
          WRITE(fp2,'(F6.1,2F12.4)') 
     .      eirtrans(i1,1),(eirtrans(i1,i2)*fact,i2=2,3)
        ENDDO

        icnt = 0
        DO i1 = 1, eirnaout
          IF (eiraout(i1,1).EQ.1.0) THEN
            add = MAX(0,(eiraout(i1,2) - 1) * asc_ncell) 
            DO i2 = 3,6
              IF (eiraout(i1,i2).GT.0.0) THEN
                icnt = icnt + 1
                ilst(icnt) = eiraout(i1,i2) + add + 1 + eirnpgdat
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        WRITE(fp2,'( A,I6)') '''Additional cells for recodring '//
     .                       'stratum data''',icnt
        WRITE(fp2,'(6I6)') (ilst(i1),i1=1,icnt)        


        IF (eirzaa.LT.0.0) THEN
          WRITE(fp2,'(A,F9.3)') '''Global volume scaling''',eirtorfrac
        ELSE
          WRITE(fp2,'(A,F9.3)') '''Global volume scaling''',1.0
        ENDIF

        IF (grdnmod.NE.0) THEN
          WRITE(fp2,'(A,I9)') '''Radial connection map''',1
        ELSE
          WRITE(fp2,'(A,I9)') '''Radial connection map''',0
        ENDIF


c
c       Advance the template file until the next section tag
c       is found:
c
21      CALL ReadLine(fp1,buffer,1,*97,*98)
        IF (buffer(1:3).NE.'***') GOTO 21
        BACKSPACE fp1

      ELSEIF (buffer(1:6).EQ.'*** 1.') THEN
c
c       Set EIRENE execution time:
c
        CALL WriteBlock01(fp1,fp2)
c        WRITE(fp2,'(A)') '*** 1. DATA FOR OPERATING MODE (DIVIMP)'
c        CALL UpdateLine1I(fp1,fp2,buffer,3,eirtime)
c        CALL TransferLine(fp1,fp2,buffer,1)

      ELSEIF (buffer(1:6).EQ.'*** 2.') THEN
c
c       Specify grid size (number of rings and cells per ring):
c
        WRITE(fp2,'(A)') '*** 2. DATA FOR STANDARD MESH (DIVIMP)'

23      CALL TransferLine(fp1,fp2,buffer,1)
        IF (buffer(1:2).EQ.'* ') GOTO 23

        CALL TransferLine(fp1,fp2,buffer,2)
        CALL UpdateLine2I(fp1,fp2,buffer,1,3,irwall-1,nks(irsep)+3)
        CALL TransferLine(fp1,fp2,buffer,3)
        CALL UpdateLine1I(fp1,fp2,buffer,1,nks(irsep)+3)

        CALL TransferLine(fp1,fp2,buffer,2)
        CALL ReadLine(fp1,buffer,3,*50,*98)

        IF     (eirzaa.EQ.-1.0) THEN
          zaa = 2.0 * PI * rxp * 100.0 * eirtorfrac
          roa = 0.0
          geostr = 'TFFF'
          eirnttra = 0
        ELSEIF (eirzaa.LT.0.0) THEN
          zaa = 360.0
          roa = 1.0E-10
          geostr = 'FTFF'
          eirnttra = eirntorseg + 1
        ELSE
          zaa = eirzaa * 100.0
          roa = 0.0
          geostr = 'TFFF'
          eirnttra = 0
        ENDIF

        WRITE(fp2,'(A4)'       ) geostr
        WRITE(fp2,'(4I6)'      ) 0,0,eirnttra,0
        WRITE(fp2,'(1P,6E12.4)') 0.0,0.0,zaa,0.0,roa,eirtorfrac

c        IF (eirzaa.EQ.-1.0) THEN
c          WRITE(fp2,'(1P,6E12.4)') 0.0,0.0,2.0*PI*rxp*100.0*eirtorfrac,
c     .                             0.0,0.0,eirtorfrac
c        ELSE
c          WRITE(fp2,'(1P,6E12.4)') 0.0,0.0,eirzaa*100.0,0.0,0.0,
c     .                             eirtorfrac
c        ENDIF


        IF (eirnsdtor.GT.1) THEN
          CALL ReadLine(fp1,buffer,2,*50,*98)         
          WRITE(fp2,'(A )') 'T'
          WRITE(fp2,'(I6)') eirnsdtor
          i1 = 1
          DO WHILE (i1.LE.eirnsdtor)
            i2 = MAX(0,MIN(6,eirnsdtor-I1+1))
            WRITE(fp2,'(1P,6E12.5)') (eirsdvol(i3),i3=i1,i1+(i2-1))
            i1 = i1 + i2
          ENDDO
c          WRITE(fp2,'(10(E12.4:))') (eirsdvol(i1),i1=1,eirnsdtor) 
        ELSE
          CALL TransferLine(fp1,fp2,buffer,2)
        ENDIF

        CALL TransferLine(fp1,fp2,buffer,1)


        WRITE(fp2,'(I6)') 1+eirnpgdat+asc_ncell*ascncut

        IF (asc_ncell.GT.0) THEN
          CALL RSet(vcel,MAXASCDAT,1.0)

          DO i1 = 2, 1+eirnpgdat
            vcel(i1) = eirpgdat(i1-1,8) * 1.0E+06
          ENDDO
          DO i1 = 2+eirnpgdat, 2+eirnpgdat+asc_ncell*ascncut
            vcel(i1) = asc_vol(i1-(2+eirnpgdat)+1) * 1.0E+06
          ENDDO

c          WRITE(0,*) 'MARK: FANCY ADDITIONAL CELL VOLUME'
          i1 = 1
          DO WHILE (i1.LE.1+eirnpgdat+asc_ncell*ascncut)
            i2 = MAX(0,MIN(6,(1+eirnpgdat+asc_ncell*ascncut)-i1+1))
c            WRITE(0,*) i2
            WRITE(fp2,'(1P,6E12.5,0P)') (vcel(i3),i3=i1,i1+(i2-1))
            i1 = i1 + i2
          ENDDO

        ELSE

c          WRITE(0,*) 'MARK: STANDARD ADDITIONAL CELL VOLUME ',eirnpgdat

          IF (eirnpgdat.GT.0) THEN
            CALL RSet(vcel,MAXASCDAT,1.0)
            DO i1 = 2, 1+eirnpgdat
              vcel(i1) = eirpgdat(i1-1,8) * 1.0E+06
            ENDDO
            i1 = 1
            DO WHILE (i1.LE.1+eirnpgdat+asc_ncell*ascncut)
              i2 = MAX(0,MIN(6,(1+eirnpgdat+asc_ncell*ascncut)-i1+1))
              WRITE(fp2,'(1P,6E12.5,0P)') (vcel(i3),i3=i1,i1+(i2-1))
              i1 = i1 + i2
            ENDDO
          ELSE
            WRITE(fp2,'(1P,E12.5,0P)') 1.0
          ENDIF

cc... preset limit here is 6 (set by EIRENE code: in subroutine INPUT
cc where the format data for the READ statement is limited to 6, and
cc in the common block PARMUSR, where NADD (maximum number of additional
cc surfaces) is also set to 6)
c          WRITE(fp2,'(1P,6E12.5,0P)') 1.0,(eirpgdat(i1,6)*1.0E+06,
c     .                                     i1=1,eirnpgdat)
c          IF (1+eirnpgdat.GT.6) CALL ER('WriteInputFile','Too many'//
c     .                                   ' add surfaces',*99)
        ENDIF

        CALL ReadLine(fp1,buffer,2,*50,*98)


      ELSEIF (buffer(1:6).EQ.'*** 3a') THEN

        CALL WriteBlock03a(fp1,fp2,1)

      ELSEIF (buffer(1:6).EQ.'*** 3b'.OR.
     .        buffer(1:6).EQ.'*** 3B') THEN
c      ELSEIF ((buffer(1:6).EQ.'*** 3b'.OR.
c     .         buffer(1:6).EQ.'*** 3B').AND.eirneut.EQ.1) THEN

        CALL WriteBlock03b(fp1,fp2,1)

      ELSEIF (buffer(1:6).EQ.'*** 4.') THEN
        CALL WriteBlock04(fp1,fp2)

      ELSEIF (buffer(1:6).EQ.'*** 5.') THEN
        CALL WriteBlock05(fp1,fp2)

      ELSEIF (buffer(1:6).EQ.'*** 6.') THEN
        CALL WriteBlock06(fp1,fp2)

      ELSEIF (buffer(1:6).EQ.'*** 7.') THEN

        CALL WriteBlock07(fp1,fp2)

      ELSEIF (buffer(1:7).EQ.'*** 10.') THEN

        WRITE(fp2,35) '*** 10. DATA FOR ADDITIONAL TALLIES - DIVIMP'
        WRITE(fp2,36) 0,0,eirntally,0,0,0
        WRITE(fp2,35) '*** 10A.'
        WRITE(fp2,35) '*** 10B.'
        WRITE(fp2,35) '*** 10C.'
        DO i1 = 1, eirntally
          WRITE(fp2,37) eirtally(i1,1)(1:LEN_TRIM(eirtally(i1,1)))
          WRITE(fp2,37) eirtally(i1,2)(1:69)//eirtally(i1,5)(1:3)
          WRITE(fp2,37) eirtally(i1,3)(1:24),
     .                  eirtally(i1,4)(1:24)
        ENDDO
        WRITE(fp2,35) '*** 10D.'
        WRITE(fp2,35) '*** 10E.'        

35      FORMAT(A)
36      FORMAT(6(I6:))
37      FORMAT(6(A:))

c...    Read to end of input block in template file:
40      CALL ReadLine(fp1,buffer,1,*50,*98)
        IF (buffer(1:6).NE.'*** 11') GOTO 40
        BACKSPACE fp1

      ELSEIF (buffer(1:6).EQ.'*** 11') THEN
        WRITE(fp2,'(A)')
     .    '*** 11. DATA FOR NUMERICAL AND GRAPHICAL OUTPUT (DIVIMP)'

        CALL TransferLine(fp1,fp2,buffer,17)
        CALL UpdateLine1I(fp1,fp2,buffer,2,nks(irsep)+2)
        CALL TransferLine(fp1,fp2,buffer,11)
        CALL UpdateLine2I(fp1,fp2,buffer,1,2,eirtrc1,eirtrc2)
        CALL TransferLine(fp1,fp2,buffer,1 )

      ELSEIF (buffer(1:6).EQ.'*** 13') THEN

c        WRITE(0,*) '-->',eirdtimv

        IF (eirdtimv.NE.0.0) THEN
          WRITE(fp2,'(A)') '*** 13. DATA FOR NONLI. AND/OR TIME DEP.'
          WRITE(fp2,'( I6)') 199999
          WRITE(fp2,'(2I6)') 0,1
          WRITE(fp2,'(1P,2E12.4)') eirdtimv,0.0
          WRITE(fp2,'(A)') '** 13A. DATA FOR SNAPSHOT TALLIES'
          WRITE(fp2,'( I6)') 0

          CALL ReadLine(fp1,buffer,1,*50,*98)
        ELSE
          CALL WriteLine(fp2,buffer)
          CALL TransferLine(fp1,fp2,buffer,1)
        ENDIF

      ELSEIF (buffer(1:6).EQ.'*** 14') THEN

        CALL WriteBlock14(fp1,fp2)

        CALL ReadLine(fp1,buffer,1,*50,*98)
        CALL ER('WriteInputFile','Expected end of file',*99)

      ELSE
c
c       Input block identifier is not recognized, so
c       copy the entire section as is:
c
        CALL WriteLine(fp2,buffer)

        GOTO 10
      ENDIF

      CALL ReadLine(fp1,buffer,1,*97,*98)
      IF (buffer(1:3).NE.'***') THEN
        CALL ER('WriteInputFile','Invalid template format',*99)
      ENDIF

      GOTO 20

50    CONTINUE

      CLOSE (fp1)
      CLOSE (fp2)

      IF (.FALSE..AND.eirneut.EQ.0) THEN
        nvesm = 0
        write(0,*)
        write(0,*) ' TEMPORARY BLANKING OF NEUTRAL WALL DATA! '
        write(0,*)
      ENDIF

c      STOP 'WRITE INPUT FILE'
      RETURN
c
c     Error code:
c
95    WRITE(0,*) 'FILE ERROR'
      STOP
97    CALL ER('WriteInputFile','Unexpected end of file',*99)
98    CALL ER('WriteInputFile','Problems reading template file',*99)
99    WRITE(EROUT,*) '  Last line read: '
      WRITE(EROUT,*) '  "',buffer,'"'
      STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock01
c
c
c
c
      SUBROUTINE WriteBlock01(fp1,fp2)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER    MAXSTRAT,MAXSDATA
      PARAMETER (MAXSTRAT=10,MAXSDATA=10)
      COMMON /EIRCOM/
     .        eir_07nsrc,eir_07opt ,eir_07stra,
     .        eir_07ind1,
     .        eir_07ind2,eir_07ind3,
     .        eir_07wght           ,eir_07data
      INTEGER eir_07nsrc,eir_07opt ,eir_07stra(MAXSTRAT),
     .        eir_07ind1(MAXSTRAT),
     .        eir_07ind2(MAXSTRAT),eir_07ind3(MAXSTRAT)
      REAL    eir_07wght          ,eir_07data(MAXSTRAT,MAXSDATA)

      INTEGER   fp1,fp2,eirntime
      CHARACTER buffer*200

      WRITE(fp2,90) '*** 1. DATA FOR OPERATING MODE (DIVIMP)'

      eirntime = 0
      IF (eirdtimv.NE.0.0) eirntime = 1

      IF     (eir_07opt.EQ.4.OR.eir_07opt.EQ.5) THEN
c...    BGK:
        WRITE(fp2,91) 2,1,eirtime,00000,0,eirniter,0,eirntime
        WRITE(fp2,90) 'FFFTF'

      ELSEIF (eir_07opt.EQ.0.OR.eir_07opt.EQ.2) THEN
c...    Standard (no BGK):
        WRITE(fp2,91) 2,1,eirtime,00100,0,0,0,eirntime
        WRITE(fp2,90) 'FFFTF'

      ELSE
        CALL ER('WriteBlock01','Invalid EIR_07OPT',*99)
      ENDIF

40    CALL ReadLine(fp1,buffer,1,*97,*98)
      IF (buffer(1:3).NE.'***') GOTO 40
      BACKSPACE fp1

      RETURN
90    FORMAT(A)
91    FORMAT(3I6,I6.5,20(I6:))
97    CALL ER('WriteBlock01','Unexpected end of file'        ,*99)
98    CALL ER('WriteBlock01','Problems reading template file',*99)
99    STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock03b
c
c
c
c
      SUBROUTINE WriteBlock03b(fp1,fp2,mode)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER fp1,fp2,mode

      REAL       TOL
      PARAMETER (TOL=1.0E-07)

      INTEGER    i1,i2,i3,i4,i5,ir,iv,id,cvesm,ndivadsur,nvp,
     .           nvoidside,voidside(MAXPTS),iliin,iacell,defiliin,
     .           defilside,ilside,defilswch,ilswch,fp3,s1,v1,cell,
     .           ilcell,in,nsdtor,nlimi1,nlimi2,nlimi3,nlimi4,nlimi5,
     .           sbgki,cbgki,haddi,hstdi,nlimi,id0,id1,id2,
     .           addtorseg1,addtorseg2,numseg,iltor,addsection,
     .           numsec,nsec
      REAL       x0,y0,r,material(4),recyct,recycf,rtmp,ztmp,
     .           x(4),y(4),z(4),
     .           defrecyct,defrecycf,rlb,zval,ewall
      LOGICAL    status,done


      INTEGER isec,ishift
      REAL    zshift


      CHARACTER  buffer*200

      INTEGER    MAXSUR
      PARAMETER (MAXSUR=MAXASC)
      INTEGER nsur,csur(MAXSUR),nfunny,isur(MAXSUR)
      REAL    xsur(4,MAXSUR),ysur(4,MAXSUR),zsur(4,MAXSUR)

      COMMON  /EIRWALCOM/ walln,wallr,wallz,wallw,ebgki
      INTEGER walln,wallw(MAXPTS),ebgki
      REAL    wallr(MAXPTS,2),wallz(MAXPTS,2),wallt(MAXPTS)

      DATA material / 9642., 1206., 18474., 904./

      DATA addsection,addtorseg1,addtorseg2 /0,0,0/


      COMMON /SURFACEMAPCOM/ surfacemap,surfacesrc
      INTEGER surfacemap(5000),surfacesrc(5000)            

      SAVE

C...  Paramaterize and move to SLCOM:
      CALL IZero(surfacemap,5000)
      CALL IZero(surfacesrc,5000)
 

c...  Check if the neutral wall is being assigned from the radial boundary of the
c     magnetic grid:
c     jdemod - added ctrap=8 - to see if it will work for cases without a PFZ wall - i.e. limiters
c
      IF     (cneur.EQ.7.AND.(ctrap.EQ.7.or.ctrap.eq.8)) THEN
        WRITE(fp2,'(A)') '*** 3b. DATA FOR ADDITIONAL SURFACES (DIVIMP)'
        WRITE(fp2,'(I6)') 0

        IF (mode.EQ.1) THEN
30        CALL ReadLine(fp1,buffer,1,*97,*98)
          IF (buffer(1:3).NE.'***') GOTO 30
          BACKSPACE fp1
        ENDIF
 
        RETURN
      ELSEIF (cneur.EQ.7.OR. ctrap.EQ.7) THEN
        CALL ER('WriteBlock03b','Code development required if the '//
     .          'neutral wall is assigned to only part of the outer '//
     .          'grid boundary',*99) 
      ENDIF


c...  Assign defaults to surface properties:
      defiliin  = 1
      defilside = 2
      defilswch = 0
      defrecyct = 1.0
      defrecycf = 1.0

      ndivadsur = 0
      nfunny    = 0
      
c...  Search for any default over-rides in the surface properties array:
      DO i1 = 1, eirnspdat
        IF (eirspdat(i1,1).EQ.-2.0) THEN
          WRITE(0,*) 'ONLY ILIIN AND RECYCT DEFAULTS CAN BE SET'
          IF (eirspdat(i1,3).NE.99.0) defiliin  = NINT(eirspdat(i1,3))
          IF (eirspdat(i1,8).NE.99.0) defrecyct = eirspdat(i1,8)
          IF (eirspdat(i1,9).NE.99.0) defrecycf = eirspdat(i1,9)
        ENDIF
      ENDDO

c...  Assemble a list of sides in the WALLPT arrays that are not to
c     be passed to EIRENE.  This is necessary in the event that 2 sides
c     overlap in DIVIMP, and it isn't a problem, but it would be
c     in EIRENE (current special case of including the bypass in the
c     C-Mod plenum):
      nvoidside = 0
      DO i1 = 1, eirnspdat
        IF (eirspdat(i1,1).EQ.  2.0.AND.
     .      eirspdat(i1,3).EQ.999.0) THEN
          nvoidside = nvoidside + 1
          voidside(nvoidside) = NINT(eirspdat(i1,2))
        ENDIF
      ENDDO


      IF (.NOT.thesis) THEN

c...    For DIII-D case.  Wall data swapping code from SLTMP.D6A.  No call
c       to BuildNeutralWall necessary:

c...    Start the NIMBUS wall at IRTRAP on the "outer" target -- this
c       is how the wall is organized when calling NIMBUS:
        i1 = 0
        DO i2 = wallpts, 1, -1
          IF (wallpt(i2,16).EQ.1.0) i1 = i2
        ENDDO
        IF (i1.EQ.0) CALL ER('WriteInputFile','Bad wall data',*99)

        walln = 0
        DO i1 = wallpts, 1, -1
          IF (i1.EQ.wallpts) THEN
            i2 = 1
          ELSE
            i2 = i1 + 1
          ENDIF
          IF (wallpt(i1,18).EQ.0.0) THEN
            walln = walln + 1
            wallr(walln,1) =  wallpt(i2,20)
            wallz(walln,1) =  wallpt(i2,21)
            wallr(walln,2) =  wallpt(i1,20)
            wallz(walln,2) =  wallpt(i1,21)
            wallw(walln)   =  wallpt(i1,17)
            wallt(walln)   = -wallpt(i1,19) * 1.38E-23 / ECH

c...        Check if side should be omitted:
            DO i3 = 1, nvoidside
              IF (voidside(i3).EQ.NINT(wallpt(i1,17))) THEN
                walln        = walln - 1 
                voidside(i3) = 0
              ENDIF
            ENDDO

          ENDIF
        ENDDO

      ELSE

        STOP 'CODE NO LONGER IN USE'

c...    C-Mod modelling:

c.....  Make it work same as main version...
        CALL ER('...','Needs fixing',*99)
        walln = 0
        nvesm = wallpts
        DO i1 = wallpts, 1, -1
          rvesm(i1,1) = wallpt2(i1,1)
          zvesm(i1,1) = wallpt2(i1,2)
          IF (i1.EQ.wallpts) THEN
            i2 = 1
          ELSE
            i2 = i1 + 1
          ENDIF
          rvesm(i1,2) = wallpt2(i2,1)
          zvesm(i1,2) = wallpt2(i2,2)
          IF (i1.EQ.wallpts) THEN
            i2 = 1
          ELSE
            i2 = i1 + 1
          ENDIF
          IF (wallpt(i1,18).EQ.0) THEN
            walln = walln + 1
            wallr(walln,1) = wallpt2(i2,1)
            wallz(walln,1) = wallpt2(i2,2)
            wallr(walln,2) = wallpt2(i1,1)
            wallz(walln,2) = wallpt2(i1,2)
          ENDIF
        ENDDO
      ENDIF

c...  Surfaces specified in DIVIMP:
      IF (grdnmod.NE.0) THEN

      ELSEIF (stopopt.EQ.121) THEN

      ELSEIF (iflexopt(8).EQ.11) THEN
c...    980116023:

      ELSE
        DO i1 = 1, eirnasdat
c...      Only want to include 2D additional surfaces:
          IF (eirasdat(i1,1).NE.1.0) CYCLE

          ndivadsur = ndivadsur + 1

          IF     (eirasdat(i1,2).EQ.1.0) THEN
c...        Data point specified in DIVIMP input file:
            walln = walln + 1
            wallr(walln,2) = eirasdat(i1,3)
            wallz(walln,2) = eirasdat(i1,4)
          ELSEIF (eirasdat(i1,2).EQ.2.0) THEN
c...        Data point references grid:
            ir = INT(ABS(eirasdat(i1,3)))
            iv = INT(    eirasdat(i1,4) )
            IF   (eirasdat(i1,3).LT.0.0) THEN
              id = korpg(1      ,ir)
            ELSE
              id = korpg(nks(ir),ir)
            ENDIF
            walln = walln + 1
            wallr(walln,2) = rvertp(iv,id)
            wallz(walln,2) = zvertp(iv,id)
          ELSE
            CALL ER('WriteInputFile','Unsupported vertex code',*99)
          ENDIF

          IF     (eirasdat(i1,5).EQ.1.0) THEN
            wallr(walln,1) = eirasdat(i1,6)
            wallz(walln,1) = eirasdat(i1,7)
          ELSEIF (eirasdat(i1,5).EQ.2.0) THEN
            ir = INT(ABS(eirasdat(i1,6)))
            iv = INT(    eirasdat(i1,7) )
            IF   (eirasdat(i1,6).LT.0.0) THEN
              id = korpg(1      ,ir)
            ELSE
              id = korpg(nks(ir),ir)
            ENDIF
            wallr(walln,1) = rvertp(iv,id)
            wallz(walln,1) = zvertp(iv,id)
          ENDIF
c...      This works at the moment, but no
c         guarantees for the future:
          wallw(walln) = nvesm + i1
        ENDDO
      ENDIF



c...  Count funny additional surfaces:
      DO i1 = 1, eirnasdat
        IF (eirasdat(i1,1).EQ.2.0.OR.
     .      eirasdat(i1,1).EQ.7.0) nfunny = nfunny+1
      ENDDO
      nfunny = nfunny + eirnsdtor - 1

c...  Build list of additional EIRENE surfaces for vacuum grid:
      nsur = 0
      CALL IZero(csur,MAXSUR)
      CALL RZero(xsur,MAXSUR*4)
      CALL RZero(ysur,MAXSUR*4)
      CALL RZero(zsur,MAXSUR*4)

      IF (asc_3Dmode.EQ.1) THEN
c...    3D vacuum grid:

c...    First turn make sure that a given side is only written once
c       to the EIRENE input file:

        DO cell = 1, asc_ncell
          DO s1 = 1, asc_nvp3D(cell)
            IF (asc_link3D(s1,cell).LT.cell) asc_link3D(s1,cell) = -2 
          ENDDO
        ENDDO

        DO cell = 1, asc_ncell

          nvp = asc_nvp3D(cell)

          DO s1 = 1, nvp

            IF (asc_link3D(s1,cell).LE.-1) CYCLE

            nsur = nsur + 1
 
            DO v1 = 1, 4
              xsur(v1,nsur) = asc_xvp3D(s1,v1,cell)
              ysur(v1,nsur) = asc_yvp3D(s1,v1,cell)
              zsur(v1,nsur) = asc_zvp3D(s1,v1,cell)
            ENDDO

            csur(nsur) = asc_link3D(s1,cell) - asc_cell(cell)

          ENDDO
        ENDDO

        WRITE(SLOUT,*) 'TRYING...'
        DO i1 = 1, nsur
         WRITE(SLOUT,'(A,2I5,1X,4(3F10.5,1X),I5)') 'ADD SURFACES: ',
     .     i1,csur(i1),(xsur(i2,i1),ysur(i2,i1),zsur(i2,i1),i2=1,4)
        ENDDO

        fp3=98
        OPEN(UNIT=fp3,FILE='dump.dat',ACCESS='SEQUENTIAL',
     .       STATUS='REPLACE',ERR=96)
        DO i1 = 1, nsur
          DO i2 = 1, 4
            i3 = i2 + 1
            IF (i3.GT.4) i3 = 1
            WRITE(fp3,'(I6,2(3E14.6,2X))') 
     .        i1,xsur(i2,i1),ysur(i2,i1),zsur(i2,i1),
     .           xsur(i3,i1),ysur(i3,i1),zsur(i3,i1)
          ENDDO
        ENDDO
        CLOSE(fp3)

      ELSE

c...    Standard 2D vacuum grid:
        DO i2 = 1, asc_ncell
          DO i3 = 1, asc_nvp(i2)
c...TAG-A
c...BUG? Feb 14, 2001
c            IF (asc_link(i3,i2).LE.0) CYCLE
            IF (asc_link(i3,i2).LE.-1) CYCLE
            status = .TRUE.
            i5 = i3 + asc_nvp(i2)
            DO i4 = 1, nsur
              IF (ABS(asc_rvp(i3,i2)-xsur(1,i4)).LT.TOL.AND.
     .            ABS(asc_rvp(i5,i2)-xsur(2,i4)).LT.TOL.AND.
     .            ABS(asc_zvp(i3,i2)-ysur(1,i4)).LT.TOL.AND.
     .            ABS(asc_zvp(i5,i2)-ysur(2,i4)).LT.TOL)
     .          status = .FALSE.
            ENDDO
            IF (status) THEN
              nsur = nsur + 1
              isur(nsur) = i2
              xsur(1,nsur) = asc_rvp (i5,i2)
              ysur(1,nsur) = asc_zvp (i5,i2)
              xsur(2,nsur) = asc_rvp (i3,i2)
              ysur(2,nsur) = asc_zvp (i3,i2)
              csur(nsur)   = asc_link(i3,i2) - asc_cell(i2)
c...Check for "hanging" cells, which don't border a solid surface.  IF
c   one is found, invert the surface:
              IF (csur(nsur).LT.0) THEN
                WRITE(0,*) 'INVERTING BGK SIDE : ',nsur
                rtmp         = xsur(1,nsur) 
                ztmp         = ysur(1,nsur)
                xsur(1,nsur) = xsur(2,nsur)
                ysur(1,nsur) = ysur(2,nsur)
                xsur(2,nsur) = rtmp
                ysur(2,nsur) = ztmp
                csur(nsur)   = -csur(nsur)
              ENDIF
            ENDIF
          ENDDO
        ENDDO

        IF (asc_3dmode.EQ.2) THEN
          DO cell = 1, ascncut-1
            nsur = nsur + 1
 	    DO v1 = 1, 4
              xsur(v1,nsur) = asc_xvp3D(1,v1,cell)
              ysur(v1,nsur) = asc_yvp3D(1,v1,cell)
              zsur(v1,nsur) = asc_zvp3D(1,v1,cell)
            ENDDO
            csur(nsur) = asc_ncell
          ENDDO
        ENDIF

        DO i1 = 1, nsur
         WRITE(SLOUT,'(A,I5,1X,2(2F12.6,1X),I5)') 'ADD SURFACES: ',
     .     i1,xsur(1,i1),ysur(1,i1),xsur(2,i1),ysur(2,i1),csur(i1)
        ENDDO

      ENDIF

      haddi = 0
      nlimi = 0
c...  Figure out surface groups:
      DO i1 = 1, eirnasdat
        IF (eirasdat(i1,1).EQ.7) haddi = haddi + 1
        IF (eirasdat(i1,1).EQ.2) nlimi = nlimi + 1
      ENDDO

c...  Increase the number of 2-class additional surfaces to reflect the number
c     of toroidal sections requested:
      haddi2 = haddi      
c      WRITE(0,*) 'DFASD:',eirnsection,haddi
c      WRITE(0,*) 'DFASD:',haddi*(eirnsection-1),addsection
      haddi  = haddi + addsection
c      haddi  = haddi * eirnsection




      haddi = haddi + addtorseg1
      nlimi = nlimi + addtorseg2

      sbgki = walln + eirnpgdat + 1 
      ebgki = sbgki + nsur - 1
      cbgki = ebgki - ascncut + 1
      haddi = ebgki + haddi
      hstdi = haddi + eirnsdtor - 1
      nlimi = hstdi + nlimi

c      WRITE(0,*) '-->',sbgki,ebgki,cbgki,haddi,hstdi,nlimi

c
c     This section details the neutral wall segments which
c     are referred to as additional surfaces.  The geometry
c     data for these surfaces comes from the neutral wall
c     data in the DIVIMP input file or perhaps the grid file,
c     in the case of JET grids.
c
c     The ADDUSR routine in EIRENE clips surfaces around the targets,
c     but this has already been done in DIVIMP, so the ADDUSR
c     routine should not be executed (see EIREADD DIVIMP option):
c
      WRITE(fp2,'(A)') '*** 3b. DATA FOR ADDITIONAL SURFACES (DIVIMP)'
c...mod1 -- add manual additional surfaces
      WRITE(fp2,'(6I6)') nlimi,sbgki,cbgki,ebgki,haddi,hstdi
c      WRITE(fp2,'(4I6)') (walln+eirnpgdat+nsur+nfunny),
c     .                   walln+1,walln+nsur,walln+nsur-ascncut+1

      nlimi1 = 0

      DO i1 = 1, walln

        nlimi1 = nlimi1 + 1

        IF (i1.LE.walln-ndivadsur) THEN
          WRITE(fp2,'(A,I4,A)') '*',i1,' :'
        ELSE
          WRITE(fp2,'(A,I4,A)') '*',i1,' : DIVIMP specified '//
     .                          'additional surface '
        ENDIF
        
        done = .FALSE.

        DO i2 = 1, eirnspdat

          IF (     eirspdat(i2,1) .EQ.  2.0    .AND.
     .             eirspdat(i2,3) .NE.999.0    .AND.
     .        NINT(eirspdat(i2,2)).EQ.wallw(i1)) THEN

            IF (eirbgk.GT.0.AND.
     .          ((eirspdat(i2,3).EQ.1.0.AND.
     .           (eirspdat(i2,6).NE.0.0.OR.eirspdat(i2,7).NE.0.0)).OR.
     .            eirspdat(i2,3).LT.0.0)) THEN
c... This triggers the search routine in EIRENE that finds the additional
c    cell that the neutral transports into when passing through a 
c    transparent surface that is not part of the vacuum mesh:
              iacell = 999
            ELSE
              iacell = 0
            ENDIF

            WRITE(fp2,'(A)') ' 2.00000E+00 1.00000E+00 1.00000E-05'//
     .                       ' 1.00000E+05'
            WRITE(fp2,'(8I6)') (NINT(eirspdat(i2,i3)),i3=3,5),0,0,6,0,
     .                         iacell

            IF (eirneut.EQ.1) THEN
c...          Accurate wall data:
              WRITE(fp2,'(6E14.7)')
     .          wallr(i1,1)*100.0,wallz(i1,1)*100.0,-1.0E+20,
     .          wallr(i1,2)*100.0,wallz(i1,2)*100.0, 1.0E+20
            ELSE
c...          Low accuracy data:
              WRITE(fp2,'(1P,6E12.4,0P)')
     .          wallr(i1,1)*100.0,wallz(i1,1)*100.0,-1.0E+20,
     .          wallr(i1,2)*100.0,wallz(i1,2)*100.0, 1.0E+20
            ENDIF
            WRITE(fp2,'(A)') '     1     0     0     0'
            WRITE(fp2,'(1P,5E12.4)') material(eirmat2),wallt(i1),0.0,
     .                               eirspdat(i2,6),eirspdat(i2,7)
            WRITE(fp2,'(1P,6E12.4)') 1.0,eirspdat(i2,8),0.0,1.0,0.5,1.0

            done = .TRUE.
          ENDIF
        ENDDO


        IF (.NOT.done) THEN
          WRITE(fp2,'(A)') ' 2.00000E+00 1.00000E+00 1.00000E-05'//
     .                     ' 1.00000E+05'
c...TEMP ABSORBING WALL
          IF (eirdebug.LE.-20.AND.eirdebug.GT.-41) THEN
            WRITE(fp2,'(A)') '     2     2     0     0     0     1'//
     .                       '     0     0     0'
          ELSE
            WRITE(fp2,'(9I6)') defiliin,2,0,0,0,1,0,0,0
          ENDIF
          IF (eirneut.EQ.1) THEN
c...        Accurate wall data:
            WRITE(fp2,'(6E14.7)')
     .        wallr(i1,1)*100.0,wallz(i1,1)*100.0,-1.0E+20,
     .        wallr(i1,2)*100.0,wallz(i1,2)*100.0, 1.0E+20
          ELSE
c...        Low accuracy data:
            WRITE(fp2,'(1P,6E12.4,0P)')
     .        wallr(i1,1)*100.0,wallz(i1,1)*100.0,-1.0E+20,
     .        wallr(i1,2)*100.0,wallz(i1,2)*100.0, 1.0E+20
          ENDIF
c          WRITE(fp2,'(6E14.7)')
c     .      wallr(i1,1)*100.0,wallz(i1,1)*100.0,-1.0E+20,
c     .      wallr(i1,2)*100.0,wallz(i1,2)*100.0, 1.0E+20
          WRITE(fp2,'(A)') '     1     0     0     0'
c Replace with an array...
          WRITE(fp2,'(1P,2E12.4)') material(eirmat2),wallt(i1)
          WRITE(fp2,'(1P,6E12.4)') defrecycf,defrecyct,0.0,1.0,0.5,1.0
        ENDIF
      ENDDO



c
c
c
      DO i1 = 1, eirnpgdat
        nlimi1 = nlimi1 + 1
        x0 = eirpgdat(i1,2)
        y0 = eirpgdat(i1,3)
        r  = eirpgdat(i1,5)
        WRITE(fp2,'(A,I4,A)') '*',nlimi1,' : (pressure gauge)'
c        WRITE(fp2,'(A,I4,A)') '*',i1+walln,' : (pressure gauge)'
        WRITE(fp2,'(A)') ' 0.00000E+00 1.00000E+00 1.00000E-05'//
     .                   ' 1.00000E+05'
c...need to improve indexing here...since it fails
c   if the presure gauge is inside an additional surface
        WRITE(fp2,'(7I6,(3X,I3),I6)') -3,0,1000,0,0,2,0,i1,0
        WRITE(fp2,'(1P,6E12.5,0P)')
     .    (x0**2.0 + y0**2.0 - r**2.0)*10000.0,
     .    -2.0*x0*100.0,-2.0*y0*100.0,0.0,1.0,1.0
        WRITE(fp2,'(1P,6E12.5,0P)') 0.0,0.0,0.0,0.0,0.0,0.0
      ENDDO


      fp3=98
      nlimi3 = 0
      OPEN(UNIT=fp3,FILE='dump.dat',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)


c...  Vacuum grid additional surfaces:
      IF (asc_3Dmode.EQ.1) THEN

        STOP 'ASC_3DMODE=1 NO LONGER IN USE'

        DO i1 = 1, nsur

          nlimi1 = nlimi1 + 1

          DO i2 = 1, 4
            x(i2) = xsur(i2,i1)*100.0
            y(i2) = ysur(i2,i1)*100.0
            z(i2) = zsur(i2,i1)*100.0
          ENDDO

          WRITE(fp2,'(A,I4,A)') '*',nlimi1,': BGK mesh (3D)'
c          WRITE(fp2,'(A,I4,A)') '*',i1+eirnpgdat+walln,': BGK mesh (3D)'
          WRITE(fp2,'(A)') ' 4.00000E+00 1.00000E+00 1.00000E-05'//
     .                     ' 1.00000E+05'
          WRITE(fp2,'(7I6,(3X,I3),I6)') -1,0,1000,0,0,2,0,csur(i1),0

          WRITE(fp2,94) (x(i2),y(i2),z(i2),i2=1,2)
          WRITE(fp2,94) (x(i2),y(i2),z(i2),i2=4,3,-1)

          WRITE(fp2,'(1P,6E12.5,0P)') 0.0,0.0,0.0,0.0,0.0,0.0

        ENDDO

      ELSE 

        DO i1 = 1, nsur

          IF ( asc_3Dmode.EQ.0.OR.
     .        (asc_3dmode.EQ.2.AND.i1.LE.nsur-(ascncut-1))) THEN

            nlimi1 = nlimi1 + 1

            WRITE(fp2,'(A,I4,A,I5)') '*',nlimi1,': BGK mesh',isur(i1)
c            WRITE(fp2,'(A,I4,A)') '*',i1+eirnpgdat+walln,': BGK mesh'

            WRITE(fp2,'(A)') ' 2.00000E+00 1.00000E+00 1.00000E-05'//
     .                       ' 1.00000E+05'
            WRITE(fp2,'(7I6,(3X,I3),I6)') -1,0,2000,0,0,2,0,csur(i1),0
            WRITE(fp2,'(1P,6E14.7,0P)')
     .        xsur(1,i1)*100.0,ysur(1,i1)*100.0,-1.0E+20,
     .        xsur(2,i1)*100.0,ysur(2,i1)*100.0, 1.0E+20
            WRITE(fp2,'(1P,6E12.5,0P)') 0.0,0.0,0.0,0.0,0.0,0.0

            IF (xsur(1,i1).GT.xsur(2,i1)) THEN
              WRITE(PINOUT,*) 'BGK R-MESH BAD:',i1,xsur(1,i1),xsur(2,i1)
            ENDIF
            IF (ysur(1,i1).GT.ysur(2,i1)) THEN
              WRITE(PINOUT,*) 'BGK Z-MESH BAD:',i1,ysur(1,i1),ysur(2,i1)
            ENDIF

          ELSE

            iltor = 0

c...        Toroidal cut surface:
            DO i2 = 1, 4
              x(i2) = xsur(i2,i1)*100.0
              y(i2) = ysur(i2,i1)*100.0
              z(i2) = zsur(i2,i1)*100.0
            ENDDO

            IF (eirntorseg.NE.0) THEN
c...          Toroidal approximation used in EIRENE, so the toroidal location of this
c             surface needs to be determined, and the vertex coordinates rotated
c             appropriately.  The segment number for each surface also needs to be
c             identified and passed to EIRENE:

c...          Initialization call:
              numseg = -1
              CALL DivAddSurface(x,y,z,numseg,0,0,0)

	      CALL DivAddSurface(x,y,z,iltor,nlimi3,fp3,6)            
            ENDIF

            nlimi1 = nlimi1 + 1

            WRITE(fp2,'(A,I4,A)') '*',nlimi1,': BGK 3D TRIM'
            WRITE(fp2,'(A)') ' 4.00000E+00 1.00000E+00 1.00000E-05'//
     .                       ' 1.00000E+05'
            WRITE(fp2,'(7I6,(3X,I3),I6)') -1,0,1000,0,iltor,2,0,
     .                                    csur(i1),0
            WRITE(fp2,95) (x(i2),y(i2),z(i2),i2=1,2)
            WRITE(fp2,95) (x(i2),y(i2),z(i2),i2=4,3,-1)
c            WRITE(fp2,94) (x(i2),y(i2),z(i2),i2=1,2)
c            WRITE(fp2,94) (x(i2),y(i2),z(i2),i2=4,3,-1)
            WRITE(fp2,'(1P,6E12.5,0P)') 0.0,0.0,0.0,0.0,0.0,0.0

c...Dump:
            IF (eirntorseg.EQ.0) THEN
              nlimi3 = nlimi3 + 1
              DO i2 = 1, 4
                i3 = i2 + 1
                IF (i3.GT.4) i3 = 1
                WRITE(fp3,'(I6,2(3E14.6,2X),I6)') 
     .            nlimi3,x(i2)/100.0,y(i2)/100.0,z(i2)/100.0,
     .                   x(i3)/100.0,y(i3)/100.0,z(i3)/100.0,3
              ENDDO
            ENDIF

          ENDIF
        ENDDO

      ENDIF


c...  Non-standard additional surfaces that can ONLY be seen from
c     additional cells:

c Add a loop here to repeat the 'funny' additional surfaces a specified
c number of times, with a simple shift to z0 for each iteration.  There would have
c to be some accounting of this before hand of course, but really, this would
c the *only thing that need repeating.  Specifying surface properties gets a little more
c difficult of course, but that can't be helped -- maybe not much more difficult, if only
c the block or section that the surface properties refers to has to be given.
c
c More difficult to assign recycling sources from is perhaps 
c
c
c
c

c      IF (eirzaa.LE.0.0.AND.eirnsection.GT.1) 
c     .  CALL ER('WriteBlock03b','Automatied additional surface '//
c     .          'periodicity only available when cylindrical extent '//
c     .          '(EIRZAA) specified',*99)

      nlimi2 = 0
      nlimi4 = nlimi1

      DO isec = 1, eirnsection

        IF (eirntorseg.NE.0) THEN 
          zshift = SNGL(DBLE(-eirzaa)*DBLE(eirtorfrac)/
     .                  DFLOAT(eirnsection)*DFLOAT(isec-1))
        ELSE
          zshift = SNGL(DBLE(eirzaa)/DFLOAT(eirnsection)*
     .                  DFLOAT(isec-1))
        ENDIF

        ishift = haddi2 * (isec - 1)

        DO i1 = 1, eirnasdat
          IF (eirasdat(i1,1).NE.7.0) CYCLE

          rlb = 2.0
          x(1) = eirasdat(i1,3)*100.0
          y(1) = eirasdat(i1,4)*100.0
          z(1) = eirasdat(i1,8)*100.0 + zshift * 100.0
          x(2) = eirasdat(i1,6)*100.0
          y(2) = eirasdat(i1,7)*100.0
          z(2) = eirasdat(i1,9)*100.0 + zshift * 100.0

          IF (eirasdat(MIN(eirnasdat,i1+1),1).EQ.-1.0) THEN
            rlb = 4.0
            x(3) = eirasdat(i1+1,3)*100.0
            y(3) = eirasdat(i1+1,4)*100.0
            z(3) = eirasdat(i1+1,8)*100.0 + zshift * 100.0
            x(4) = eirasdat(i1+1,6)*100.0
            y(4) = eirasdat(i1+1,7)*100.0
            z(4) = eirasdat(i1+1,9)*100.0 + zshift * 100.0
          ENDIF        

c...      Close the gap over closed C-Mod ports. Nasty hardcoded business:
          IF (.TRUE..AND.
     .        ((eirnsection.EQ.4 .AND.isec.NE.1).OR.
     .         (eirnsection.EQ.2 .AND.isec.NE.1).OR.
     .         (eirnsection.EQ.3 .AND.isec.EQ.2).OR.
     .         (eirnsection.EQ.10.AND.
     .          (isec.EQ.1.OR.isec.EQ.2.OR.isec.EQ.3.OR.isec.EQ.6.OR.
     .           isec.EQ.9)))) THEN
            IF ((i1.GE.72 .AND.i1.LE.103).OR.
     .          (i1.GE.140.AND.i1.LE.151)) THEN
c...          Pass on the gussets and port boundaries:
            ELSE
              IF (eirasdat(i1  ,8).EQ.0.0908) z(1)=z(1)+3.72
              IF (eirasdat(i1  ,9).EQ.0.0908) z(2)=z(2)+3.72
              IF (eirasdat(i1+1,8).EQ.0.0908) z(3)=z(3)+3.72
              IF (eirasdat(i1+1,9).EQ.0.0908) z(4)=z(4)+3.72

              IF (eirasdat(i1  ,8).EQ.0.1692) z(1)=z(1)-3.82
              IF (eirasdat(i1  ,9).EQ.0.1692) z(2)=z(2)-3.82
              IF (eirasdat(i1+1,8).EQ.0.1692) z(3)=z(3)-3.82
              IF (eirasdat(i1+1,9).EQ.0.1692) z(4)=z(4)-3.82

            ENDIF
          ENDIF

          nlimi2 = nlimi2 + 1

c...      If a surface is infinite in toroidal extent, then do not duplicate
c         it in every toroidal section (this is only if more than one toroidal
c         section is defined by, i.e. EIRNSECTION > 1):
          IF (isec.GT.1.AND.
     .        (ABS(z(1)).GT.1.0E+10.OR.ABS(z(2)).GT.1.0E+10.OR.
     .         ABS(z(3)).GT.1.0E+10.OR.ABS(z(4)).GT.1.0E+10)) THEN
c            WRITE(0,*) 'CYCLING:',isec,
c     .                            nlimi2+nvesm+nvesp+eirnpgdat-ishift
            CYCLE
          ENDIF

          nlimi1 = nlimi1 + 1

c...      Tally the number of additional surfaces added if more than 1 toroidal
c         section is being specified:
          IF (isec.GT.1) addsection = addsection + 1

          iliin  = defiliin
          ilside = defilside
          ilswch = defilswch
          iltor  = 0
          ilcell = 0
          recyct = defrecyct
          recycf = defrecycf
          ewall  = eirtemp2
	  
c...      Search for non-default surface properties:
          DO i2 = 1, eirnspdat
            IF (eirspdat(i2,1).EQ.2.0.AND.
     .         ((eirspdat(i2,2 ).EQ.REAL(nvesm+nvesp+nlimi2-ishift).AND.
     .           eirspdat(i2,10).EQ.0.0).OR.
     .          (eirspdat(i2,2 ).EQ.REAL(nvesm+nvesp+nlimi2-ishift).AND.
     .           eirspdat(i2,10).EQ.REAL(isec)))) THEN
	  
c              WRITE(PINOUT,*) 'FOUND SPDATA FOR ADD SURF',nlimi1,nlimi2

c              IF (eirspdat(i2,10).EQ.REAL(isec)) THEN
c                WRITE(0,*) 'SUPER OVER-RIDE:',i2,nlimi1
c              ENDIF
            
              iliin  = NINT(eirspdat(i2,3))
              ilside = NINT(eirspdat(i2,4))
              ilswch = NINT(eirspdat(i2,5))
              recyct = eirspdat(i2,8)
              recycf = eirspdat(i2,9)
	  
c...          Additional cell index switching.  Use the serch routine
c             by default:
              IF (ilswch.EQ.1000) THEN
                IF (iliin .EQ.-9.0) THEN
c...              A transparent switching surface that can tolerate not finding
c                 and additional cell after calling CHKVAC in EIRENE, assuming
c                 that additional cell 1 has been entererd:
                  iliin  = -1.0
                  ilcell = 997
                ELSE
                  ilcell = 999       
                ENDIF
              ENDIF
              IF (ilswch.EQ.3000) THEN
                ilcell = 998998
                ilswch = 1000
              ENDIF
            ENDIF
          ENDDO
	  
c...      Search for non-default surface temperature:
          id0 = nvesm + nvesp + nlimi2 - ishift
          DO i2 = 1, nwltemp
            id1 = walltemp(i2,1)
            id2 = walltemp(i2,2)
            IF (id0.GE.id1.AND.id0.LE.id2) THEN
              ewall = -walltemp(i2,3) * 1.38E-23 / ECH
            ENDIF      
          ENDDO
	  
          IF (eirntorseg.NE.0) THEN
c...        Toroidal approximation used in EIRENE, so the toroidal location of this
c           surface needs to be determined, and the vertex coordinates rotated
c           appropriately.  The segment number for each surface also needs to be
c           identified and passed to EIRENE:

            numseg = -1
            CALL DivAddSurface(x,y,z,numseg,0,0,0)
	
            DO i4 = 1, numseg

              CALL DivAddSurface(x,y,z,iltor,nlimi3,fp3,7+MIN(0,iliin))            

              nlimi4 = nlimi4 + 1
 
              IF (i4.GT.1) THEN
                addtorseg1 = addtorseg1 + 1
              ELSE
                eirnlimi(i1+ishift) = nlimi4
              ENDIF

              surfacemap(nlimi4) = nlimi2 + nvesm + nvesp + eirnpgdat
              surfacesrc(nlimi4) = surfacemap(nlimi4) - ishift

              WRITE(fp2,'(A,I4,A,2I4,A)') 
     .          '*',nlimi4,' (',nlimi1,i4,') : funny'
              WRITE(fp2,93) rlb,1.0,1.0E-05,1.0E+05
              WRITE(fp2,92) iliin,ilside,ilswch,0,iltor,2,0,ilcell,0
              IF     (rlb.EQ.2.0) THEN
                WRITE(fp2,95) (x(i2),y(i2),z(i2),i2=1,2)
              ELSEIF (rlb.EQ.4.0) THEN
                WRITE(fp2,95) (x(i2),y(i2),z(i2),i2=1,2)
                WRITE(fp2,95) (x(i2),y(i2),z(i2),i2=4,3,-1)
c                WRITE(fp2,94) (x(i2),y(i2),z(i2),i2=1,2)
c                WRITE(fp2,94) (x(i2),y(i2),z(i2),i2=4,3,-1)
              ELSE
                CALL ER('Block3b','Unknown RLB value',*99)
              ENDIF


              WRITE(fp2,92) 1,0,0,0
              WRITE(fp2,94) material(eirmat2),ewall
              WRITE(fp2,94) recycf,recyct,0.0,1.0,0.5,1.0

            ENDDO


          ELSE
            nlimi4 = nlimi4 + 1

            surfacemap(nlimi4) = nlimi2 + nvesm + nvesp + eirnpgdat
            surfacesrc(nlimi4) = surfacemap(nlimi4)

            WRITE(fp2,'(A,I4,A)') '*',nlimi4,' : funny'
            WRITE(fp2,93) rlb,1.0,1.0E-05,1.0E+05
            WRITE(fp2,92) iliin,ilside,ilswch,0,0,2,0,ilcell,0
            IF     (rlb.EQ.2.0) THEN
              WRITE(fp2,95) (x(i2),y(i2),z(i2),i2=1,2)
            ELSEIF (rlb.EQ.4.0) THEN
              WRITE(fp2,95) (x(i2),y(i2),z(i2),i2=1,2)
              WRITE(fp2,95) (x(i2),y(i2),z(i2),i2=4,3,-1)
c              WRITE(fp2,94) (x(i2),y(i2),z(i2),i2=1,2)
c              WRITE(fp2,94) (x(i2),y(i2),z(i2),i2=4,3,-1)
            ELSE
              CALL ER('Block3b','Unknown RLB value',*99)
            ENDIF
            WRITE(fp2,92) 1,0,0,0                                    
            WRITE(fp2,94) material(eirmat2),ewall
            WRITE(fp2,94) recycf,recyct,0.0,1.0,0.5,1.0
	  
c...        Send surface geometry data to .dump file:
            nlimi3 = nlimi3 + 1
            DO i2 = 1, 4
              i3 = i2 + 1
              IF (i3.GT.4) i3 = 1
              WRITE(fp3,'(I6,2(3E14.6,2X),I6)') 
     .          nlimi3,x(i2)/100.0,y(i2)/100.0,z(i2)/100.0,
     .                 x(i3)/100.0,y(i3)/100.0,z(i3)/100.0,1
            ENDDO
          ENDIF

        ENDDO

      ENDDO

92    FORMAT(10(I6:))
93    FORMAT(1P,4(E12.5:),0P)       
94    FORMAT(1P,6(E12.4:),0P)       
95    FORMAT(1P,6(E14.7:),0P)       

c...  Standard grid NBLOCK switching additional surfaces:
      nsdtor = 0
      DO i1 = 1, eirnasdat
        IF (eirasdat(i1,1).NE.6.0) CYCLE

        nsdtor = 1

        IF (i1+2.LE.eirnasdat.AND.eirasdat(i1+2,1).EQ.-2.0) 
     .    nsdtor = 1 + eirasdat(i1+2,2)

        zval = eirasdat(i1,8)

        DO in = 1, nsdtor        

          IF (in.GT.1) zval = zval + eirasdat(i1+2,3)

          nlimi1 = nlimi1 + 1
          nlimi4 = nlimi4 + 1

          rlb = 2.0
          x(1) = eirasdat(i1,3)*100.0
          y(1) = eirasdat(i1,4)*100.0
          z(1) = zval          *100.0
          x(2) = eirasdat(i1,6)*100.0
          y(2) = eirasdat(i1,7)*100.0
          z(2) = zval          *100.0
          IF (eirasdat(MIN(eirnasdat,i1+1),1).EQ.-1.0) THEN
            rlb = 4.0
            x(3) = eirasdat(i1+1,3)*100.0
            y(3) = eirasdat(i1+1,4)*100.0
            z(3) = zval            *100.0
            x(4) = eirasdat(i1+1,6)*100.0
            y(4) = eirasdat(i1+1,7)*100.0
            z(4) = zval            *100.0
          ENDIF        

          iliin  = -1
          ilside =  0
          ilswch =  1000
          ilcell =  1000

          IF (eirntorseg.NE.0) THEN
c...        Toroidal approximation used in EIRENE, so the toroidal location of this
c           surface needs to be determined, and the vertex coordinates rotated
c           appropriately.  The segment number for each surface also needs to be
c           identified and passed to EIRENE:
	    
c...        Initialization call:
            numseg = -1
            CALL DivAddSurface(x,y,z,numseg,0,0,0)
	    
	    CALL DivAddSurface(x,y,z,iltor,nlimi3,fp3,8)            
c            STOP 'TORSEG work to do'
          ENDIF

          WRITE(fp2,'(A,I4,A)') '*',nlimi4,':NBLOCK'
	  WRITE(fp2,93) rlb,1.0,1.0E-05,1.0E+05
          WRITE(fp2,92) iliin,ilside,ilswch,0,iltor,2,0,ilcell,0
          IF     (rlb.EQ.2.0) THEN
            WRITE(fp2,95) (x(i2),y(i2),z(i2),i2=1,2)
          ELSEIF (rlb.EQ.4.0) THEN
            WRITE(fp2,95) (x(i2),y(i2),z(i2),i2=1,2)
            WRITE(fp2,95) (x(i2),y(i2),z(i2),i2=4,3,-1)
c            WRITE(fp2,94) (x(i2),y(i2),z(i2),i2=1,2)
c            WRITE(fp2,94) (x(i2),y(i2),z(i2),i2=4,3,-1)
          ELSE
            CALL ER('Block3B','Unknown RLB value',*99)
	  ENDIF

c...Dump:
          IF (eirntorseg.EQ.0) THEN
            nlimi3 = nlimi3 + 1
            DO i2 = 1, 4
              i3 = i2 + 1
              IF (i3.GT.4) i3 = 1
              WRITE(fp3,'(I6,2(3E14.6,2X),I6)') 
     .          nlimi3,x(i2)/100.0,y(i2)/100.0,z(i2)/100.0,
     .                 x(i3)/100.0,y(i3)/100.0,z(i3)/100.0,4
            ENDDO
          ENDIF

        ENDDO

      ENDDO

c...  Check that standard grid NBLOCK switching surface was found if
c     there is a gap in a non-standard default surface:
      IF (eirntrans.NE.0.AND.nsdtor.EQ.0) 
     .  CALL ER('Block07b','Transparent span in non-default '//
     .                     'standard surface and standard grid '//
     .                     'switching surface not found',*99)

      nlimi1 = nlimi4
      nlimi5 = nlimi2

c...  Non-standard additional surfaces that can be seen from everywhere:
      DO i1 = 1, eirnasdat
        IF (eirasdat(i1,1).NE.2.0) CYCLE

        nlimi1 = nlimi1 + 1        
        nlimi2 = nlimi2 + 1        

        rlb = 2.0
        x(1) = eirasdat(i1,3)*100.0
        y(1) = eirasdat(i1,4)*100.0
        z(1) = eirasdat(i1,8)*100.0
        x(2) = eirasdat(i1,6)*100.0
        y(2) = eirasdat(i1,7)*100.0
        z(2) = eirasdat(i1,9)*100.0
        IF (eirasdat(MIN(eirnasdat,i1+1),1).EQ.-1.0) THEN
          rlb = 4.0
          x(3) = eirasdat(i1+1,3)*100.0
          y(3) = eirasdat(i1+1,4)*100.0
          z(3) = eirasdat(i1+1,8)*100.0
          x(4) = eirasdat(i1+1,6)*100.0
          y(4) = eirasdat(i1+1,7)*100.0
          z(4) = eirasdat(i1+1,9)*100.0
        ENDIF        
        iliin  = defiliin
        ilside = defilside
        ilswch = defilswch
        iltor  = 0
        ilcell = 0
        recyct = defrecyct
        recycf = defrecycf
        ewall  = eirtemp2

c...    Search for non-default surface properties:
        IF (eirnsection.NE.1) 
     .    WRITE(0,*) 'FUNNY2 INDEX:',nlimi2,nvesm+nvesp+nlimi2

        DO i2 = 1, eirnspdat
          IF (eirspdat(i2,1).EQ.2.0                     .AND.
     .        eirspdat(i2,2).EQ.REAL(nvesm+nvesp+nlimi2)) THEN
c            WRITE(0     ,*) 'FOUND SPDATA FOR ADD SURF 2',nlimi1,nlimi2
            WRITE(PINOUT,*) 'FOUND SPDATA FOR ADD SURF 2',nlimi1,nlimi2
            iliin  = NINT(eirspdat(i2,3))
            ilside = NINT(eirspdat(i2,4))
            ilswch = NINT(eirspdat(i2,5))
            recyct = eirspdat(i2,8)
            recycf = eirspdat(i2,9)
c...        Additional cell index switching.  Use the serch routine
c           by default:
            IF (ilswch.EQ.1000) ilcell = 999
            IF (ilswch.EQ.3000) THEN
              ilcell = 998998
              ilswch = 1000
            ENDIF
          ENDIF
        ENDDO

c...    Search for non-default surface temperature:
        id0 = nvesm + nvesp + nlimi2
        DO i2 = 1, nwltemp
          id1 = walltemp(i2,1)
          id2 = walltemp(i2,2)
          IF (id0.GE.id1.AND.id0.LE.id2) 
     .      ewall = -walltemp(i2,3) * 1.38E-23 / ECH
        ENDDO

        IF (eirntorseg.NE.0) THEN
c...      Toroidal approximation used in EIRENE, so the toroidal location of this
c         surface needs to be determined, and the vertex coordinates rotated
c         appropriately.  The segment number for each surface also needs to be
c         identified and passed to EIRENE:

          numseg = -1
          CALL DivAddSurface(x,y,z,numseg,0,0,0)
          IF (numseg.GT.1) THEN
            numsec = eirnsection
            zshift = SNGL(DBLE(-eirzaa)*DBLE(eirtorfrac)/
     .                    DFLOAT(eirnsection))*100.0
          ELSE
            numsec = 1
            zshift = 0.0
          ENDIF
 
          DO nsec = 1, numsec
            IF (numsec.GT.1) THEN
c...          Reassign surface data and shift:
              x(1) = eirasdat(i1,3)*100.0
              y(1) = eirasdat(i1,4)*100.0
              z(1) = eirasdat(i1,8)*100.0 + zshift * REAL(nsec-1)
              x(2) = eirasdat(i1,6)*100.0
              y(2) = eirasdat(i1,7)*100.0
              z(2) = eirasdat(i1,9)*100.0 + zshift * REAL(nsec-1)
              IF (eirasdat(MIN(eirnasdat,i1+1),1).EQ.-1.0) THEN
                x(3) = eirasdat(i1+1,3)*100.0
                y(3) = eirasdat(i1+1,4)*100.0
                z(3) = eirasdat(i1+1,8)*100.0 + zshift * REAL(nsec-1)
                x(4) = eirasdat(i1+1,6)*100.0
                y(4) = eirasdat(i1+1,7)*100.0
                z(4) = eirasdat(i1+1,9)*100.0 + zshift * REAL(nsec-1)
              ENDIF        
c...          Reset number of toroidal segments the surface will be
c             divided into:
              numseg = -1
              CALL DivAddSurface(x,y,z,numseg,0,0,0)
              IF (nsec.GT.1) addtorseg2 = addtorseg2 + 1
            ENDIF

            nlimi5 = nlimi5 + 1
	
            DO i4 = 1, numseg

              nlimi4 = nlimi4 + 1

              CALL DivAddSurface(x,y,z,iltor,nlimi3,fp3,9)            
 
              IF (i4.GT.1) addtorseg2 = addtorseg2 + 1
  
              surfacemap(nlimi4) = nlimi5 + nvesm + nvesp + eirnpgdat
              surfacesrc(nlimi4) = nlimi2 + nvesm + nvesp + eirnpgdat

              WRITE(fp2,'(A,I4,A,2I4,A)') 
     .          '*',nlimi4,' (',nlimi1,i4,') : funny2'
              WRITE(fp2,93) rlb,1.0,1.0E-05,1.0E+05
              WRITE(fp2,92) iliin,ilside,ilswch,0,iltor,2,0,ilcell,0
              IF     (rlb.EQ.2.0) THEN
                WRITE(fp2,95) (x(i2),y(i2),z(i2),i2=1,2)
              ELSEIF (rlb.EQ.4.0) THEN
                WRITE(fp2,95) (x(i2),y(i2),z(i2),i2=1,2)
                WRITE(fp2,95) (x(i2),y(i2),z(i2),i2=4,3,-1)
c                WRITE(fp2,94) (x(i2),y(i2),z(i2),i2=1,2)
c                WRITE(fp2,94) (x(i2),y(i2),z(i2),i2=4,3,-1)
              ELSE
                CALL ER('Block3b','Unknown RLB value',*99)
              ENDIF
              WRITE(fp2,92) 1,0,0,0
              WRITE(fp2,94) material(eirmat2),ewall
              WRITE(fp2,94) recycf,recyct,0.0,1.0,0.5,1.0

            ENDDO

          ENDDO

        ELSE

          surfacemap(nlimi1) = nlimi2 + nvesm + nvesp + eirnpgdat
          surfacesrc(nlimi1) = surfacemap(nlimi1)

          WRITE(fp2,'(A,I4,A)') '*',nlimi1,' : funny2'
          WRITE(fp2,93) rlb,1.0,1.0E-05,1.0E+05
          WRITE(fp2,92) iliin,ilside,ilswch,0,0,2,0,ilcell,0
          IF     (rlb.EQ.2.0) THEN
            WRITE(fp2,95) (x(i2),y(i2),z(i2),i2=1,2)
          ELSEIF (rlb.EQ.4.0) THEN
            WRITE(fp2,95) (x(i2),y(i2),z(i2),i2=1,2)
            WRITE(fp2,95) (x(i2),y(i2),z(i2),i2=4,3,-1)
c            WRITE(fp2,94) (x(i2),y(i2),z(i2),i2=1,2)
c            WRITE(fp2,94) (x(i2),y(i2),z(i2),i2=4,3,-1)
          ELSE
            CALL ER('Block3b','Unknown RLB value',*99)
          ENDIF
          WRITE(fp2,92) 1,0,0,0
          WRITE(fp2,94) material(eirmat2),ewall
          WRITE(fp2,94) recycf,recyct,0.0,1.0,0.5,1.0

c...      Send geometry data to .dump file:
          nlimi3 = nlimi3 + 1
          DO i2 = 1, 4
            i3 = i2 + 1
            IF (i3.GT.4) i3 = 1
            WRITE(fp3,'(I6,2(3E14.6,2X),I6)') 
     .        nlimi3,x(i2)/100.0,y(i2)/100.0,z(i2)/100.0,
     .               x(i3)/100.0,y(i3)/100.0,z(i3)/100.0,6
          ENDDO

        ENDIF

      ENDDO












      CLOSE (fp3)

c
c
c


      IF (mode.EQ.1) THEN
40      CALL ReadLine(fp1,buffer,1,*97,*98)
        IF (buffer(1:3).NE.'***') GOTO 40
        BACKSPACE fp1

c...    Clear variables that store the number of additional surfaces
c       created when the toroidal geometry approximation is used:
        addtorseg1 = 0
        addtorseg2 = 0
        addsection = 0
      ENDIF




c      STOP 'skdfjldj'



      RETURN
96    CALL ER('WriteInputFile','Cannot create dump file',*99)
97    CALL ER('WriteInputFile','Unexpected end of file',*99)
98    CALL ER('WriteInputFile','Problems reading template file',*99)
99    WRITE(EROUT,*) '  Last line read: '
      WRITE(EROUT,*) '  "',buffer,'"'
      STOP
      END




c
c ======================================================================
c
c subroutine: WriteBlock03a
c
c
c
c
      SUBROUTINE WriteBlock03a(fp1,fp2,mode)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER    MAXSTRAT,MAXSDATA
      PARAMETER (MAXSTRAT=10,MAXSDATA=10)
      COMMON /EIRCOM/
     .        eir_07nsrc,eir_07opt ,eir_07stra,
     .        eir_07ind1,
     .        eir_07ind2,eir_07ind3,
     .        eir_07wght           ,eir_07data
      INTEGER eir_07nsrc,eir_07opt ,eir_07stra(MAXSTRAT),
     .        eir_07ind1(MAXSTRAT),
     .        eir_07ind2(MAXSTRAT),eir_07ind3(MAXSTRAT)
      REAL    eir_07wght           ,eir_07data(MAXSTRAT,MAXSDATA)


      INTEGER    MAXSUR,     MAXTAR    ,MAXREG
      PARAMETER (MAXSUR=MAXASC,MAXTAR=100,MAXREG=10)
      INTEGER nsur(MAXREG),ssur(MAXSUR,MAXREG),esur(MAXSUR,MAXREG),
     .                     csur(MAXSUR,MAXREG)

      INTEGER fp1,fp2,mode

      INTEGER     i1,i2,i3,ir,ir1,ik,ik1,ik2,id,ikm,iras,nstsi,
     .            itar,ilswch3,nlist,irlist(MAXNRS),ik1list(MAXNRS),
     .            ik2list(MAXNRS),
     .            ntar1,ttar1(MAXTAR),btar1(MAXTAR),etar1(MAXTAR),
     .            ntar2,ttar2(MAXTAR),btar2(MAXTAR),etar2(MAXTAR)
      INTEGER     iliin,ilside,ilcell,ring1,ring2,cell,cell2,nnss,
     .            ilblck,ilacll,ilblck2
      LOGICAL     output,load
      CHARACTER   buffer*200
      CHARACTER*6 ilswch(2),ilswch2
      REAL        material(4)

      INTEGER    defiliin,defilside,defilswch
      REAL       defrecyct,defrecycf
  
      DATA material / 9642., 1206., 18474., 904./

      SAVE


      output = .FALSE.
c      output = .TRUE.


c...  Assign defaults to surface properties:
      defiliin  = 1
      defilside = 2
      defilswch = 0
      defrecyct = 1.0
      defrecycf = 1.0

c...  Search for any default over-rides in the surface properties array:
      DO i1 = 1, eirnspdat
        IF (eirspdat(i1,1).EQ.-3.0) THEN
          WRITE(0,*) 'ONLY ILIIN, RECYCT AND RECYCF DEFAULTS CAN BE SET'
          IF (eirspdat(i1,3).NE.99.0) defiliin  = NINT(eirspdat(i1,3))
          IF (eirspdat(i1,8).NE.99.0) defrecyct = eirspdat(i1,8)
          IF (eirspdat(i1,9).NE.99.0) defrecycf = eirspdat(i1,9)
        ENDIF
      ENDDO

c...  Set up surface temperatures (this is not done in TAUIN1 until after
c     the plasma solution has been assigned, so a call to DOTEMP has
c     to be added here):
      CALL DoTemp

c...  Confirm that all target segment temperatures are set to the 
c     default target temperature (no WLTEMP array over-ride) since
c     EIRENE presently does a poor job of specifying individual 
c     target segment temperatures:
      DO i1 = 1, wallpts
        IF (wallpt(i1,18).NE.0.0.AND.wallpt(i1,19).NE.ctargt) 
     .    CALL ER('WriteBlock03a','Target temperature not equal to '//
     .            'default temperature',*99)
      ENDDO

      IF (.TRUE..OR.eirbgk.EQ.2.OR.eirbgk.EQ.3.OR.eirbgk.EQ.4) THEN
c      IF (eir_07opt.EQ.4.OR.eir_07opt.EQ.5) THEN

c...   Forcing reaction rates to *ALWAYS* be calculated
c      in the vacuum (temporary).

c       WRITE(0,*) 'ILSWCH OVERRIDE IN EIRENE.D6A'

c       ilswch(1) = '010200'
c       ilswch(2) = '020100'

       ilswch(1) = '010000'
       ilswch(2) = '020000'
      ELSE
        ilswch(1) = '010020'
        ilswch(2) = '020010'
      ENDIF

c...  Turn on NBLOCK search in EIRENE for non-default standard
c     grid surfaces:
      IF (eirnsdtor.GT.1) THEN
        ilblck = 999
      ELSE
        ilblck = 001
      ENDIF

      IF (stopopt3.EQ. 6.OR.stopopt3.EQ.7.OR.stopopt3.EQ.12.OR.
     .    stopopt3.EQ.13.OR.
     .    eirbgk.EQ.2) THEN

c...    Process vacuum grid:

        STOP 'STOP: VACUUM GRID METHOD NO LONGER SUPPORTED'
c        CALL MapToStandardGrid(nsur,ssur,esur,csur,MAXSUR,MAXREG)

c...MESS (LIKE EVERYTHING ELSE HERE):
      ELSEIF (eirbgk.EQ.3.OR.eirbgk.EQ.4) THEN
        nsur(1) = 2
      ELSE
c...    No vacuum grid:
        nsur(1) = 2
      ENDIF



c...  Specifying surface properties:
      IF (eirnspdat.GT.0) THEN
        ntar1 = 0
      ELSE
        ntar1 = 0
      ENDIF


      IF (eirnspdat.GT.0) THEN

        ntar2    =  0
        ttar2(1) = -1
        DO ir = 1, irwall-1

          IF (ntar2.GT.0.AND.ir.LE.etar2(ntar2)) CYCLE

          itar = 0
          DO i1 = 1, eirnspdat
            IF (    eirspdat(i1,1) .EQ.1.0.AND.
     .          INT(eirspdat(i1,2)).EQ.ir ) itar = i1
          ENDDO
c...      Default definition:
          IF (ttar2(1).EQ.-1) THEN
            ntar2 = ntar2 + 1
            ttar2(ntar2) = 0
            btar2(ntar2) = ir
            etar2(ntar2) = ir
          ELSE
            IF (itar.NE.0) THEN
c...          Input file definition:
              ntar2 = ntar2 + 1
              ttar2(ntar2) = itar
              btar2(ntar2) = ir
              IF (itar+1.LE.eirnspdat.AND.
     .            eirspdat(itar+1,1).EQ.-1.0) THEN
                etar2(ntar2) = INT(eirspdat(itar+1,2))
              ELSE
                etar2(ntar2) = ir
              ENDIF
            ELSE
              IF (ttar2(ntar2).GT.0) THEN
                ntar2 = ntar2 + 1
                ttar2(ntar2) = 0
                btar2(ntar2) = ir
                etar2(ntar2) = ir
              ELSE
                etar2(ntar2) = etar2(ntar2) + 1
              ENDIF
            ENDIF
          ENDIF
        ENDDO

      ELSE
        ntar2 = 0
      ENDIF


c...  Count the number of non-standard surfaces specified in the 
c     DIVIMP input file:
      nnss = 0
      DO i1 = 1, eirnasdat
        IF (eirasdat(i1,1).EQ.3.0.OR.eirasdat(i1,1).EQ.4.0) nnss=nnss+1
      ENDDO

c
c     Configure standard surfaces in the grid, such as targets,
c     cut points and core/wall/pfz boundaries.  The settings
c     for neutral production at the targets are currently
c     hard coded and so are not copied from the template file:
c
      WRITE(fp2,'(A)') '*** 3a. DATA FOR NON DEFAULT STANDARD '//
     .                 'SURFACES (DIVIMP)'

      IF (grdnmod.NE.0) THEN

        WRITE(fp2,'(I6)') nstsi
        nstsi = 0

      ELSEIF (nbr.GT.0) THEN
c       This will have to be improved for more general wall
c       contact which has more than one broken region...
        IF (stopopt.EQ.100.OR.stopopt.EQ.101) THEN
          WRITE(fp2,'(I6)') 6+nsur(1)+MAX(0,ntar1-1)+MAX(0,ntar2-1)

        ELSEIF (stopopt.EQ.121) THEN
          IF     (stopopt2.EQ. 9) THEN
c...        Plenum:
            STOP 'STOP: OLD'
          ELSEIF (stopopt2.EQ.10) THEN
c...        PFZ and plenum:
            STOP 'STOP: OLD'
          ELSE
c...        PFZ:
            WRITE(fp2,'(I6)') 5+nsur(1)+MAX(0,ntar1-1)+MAX(0,ntar2-1)+
     .                        nnss
          ENDIF
        ELSE
          WRITE(fp2,'(I6)') 5+nsur(1)+MAX(0,ntar1-1)+MAX(0,ntar2-1)+
     .                      nnss
        ENDIF

      ELSE
        IF (stopopt.EQ.100.OR.stopopt.EQ.101) THEN
          WRITE(fp2,'(I6)') 5+nsur(1)+MAX(0,ntar1-1)+MAX(0,ntar2-1)+
     .                      nnss
        ELSE
          WRITE(fp2,'(I6)') 4+nsur(1)+MAX(0,ntar1-1)+MAX(0,ntar2-1)+
     .                      nnss
        ENDIF
      ENDIF

c...add
      iras = 2
      DO WHILE (idring(iras).EQ.-2)
        iras = iras + 1
      ENDDO

c      WRITE(0,*) 'CORE AS = ',iras-1

      nstsi = 1
      WRITE(fp2,'(A,I2,A)') '* RADIAL FLUXSURFACE ',1,
     .                      ' ABSORBING (PLASMA CORE)'
      IF (grdnmod.NE.0) THEN
        WRITE(fp2,'(5I6)') nstsi,1,iras-1,cutpt1,cutpt2+1
      ELSE
        WRITE(fp2,'(5I6)') nstsi,1,iras-1,ikto+1,ikti+2
      ENDIF
c...rubber core ring
      IF (iflexopt(5).EQ.30) THEN
        WRITE(0,*) 'USING A RUBBER CORE'
        WRITE(fp2,'(8I6)') 3, 0,0,0,0,6,0,0
      ELSE
        WRITE(fp2,'(8I6)') 2,-2,0,0,0,6,0,0
      ENDIF

      IF (stopopt.EQ.100.OR.stopopt.EQ.101) THEN
        CALL WN('WriteInputFile','EIRENE needs further development'//
     .                           ' (see FOLNEUT)')
        nstsi = nstsi + 1
        WRITE(fp2,'(A,I2,A)') '* RADIAL FLUXSURFACE ',1,' ABSORBING'//
     .                        ' (SEPARATRIX)'
        WRITE(fp2,'(5I6)') nstsi,1,irsep-1,1,nks(irsep)
        WRITE(fp2,'(8I6)') 2,0,0,0,0,6,0,0
      ENDIF

      iras = irtrap + 1
      DO WHILE (idring(iras).EQ.-2)
        iras = iras + 1
      ENDDO

c      WRITE(0,*) 'PFZ AS = ',iras-irtrap

      IF ((stopopt3.NE. 6.AND.stopopt3.NE.7.AND.stopopt3.NE.12.AND.
     .     stopopt3.NE.13.AND.eirbgk.LT.2.AND.stopopt2.NE.10).OR.
     .     stopopt2.EQ.9.OR.
     .     vacregion(1).EQ.0) THEN


c...  Search surface properties array to see if the properties of the outer 
c     surface of the grid are being modified:

        IF (cneur.EQ.7.AND.ctrap.EQ.7) THEN
c...      Surface settings when the outer radial surface of the
c         grid is used as the neutral wall:
          iliin   =  1
          ilside  =  0
          ilswch2 = '     0'
          ilblck2 =  0
          ilacll  =  0
        ELSE
c...      Default transparent surface for switching regions as
c         a neutral enters and exits the standard (magnetic) grid:
          iliin   = -1
          ilside  =  0
          ilswch2 = ilswch(1)
          ilblck2 = ilblck
          ilacll  =  1
        ENDIF

        DO i1 = 1, eirnspdat
          IF (eirspdat(i1,1).EQ.3.0.AND.eirspdat(i1,2).EQ.1.0) THEN
            iliin   =  NINT(eirspdat(i1,3))
            ilside  =  NINT(eirspdat(i1,4))
            ilswch2 =  '     0'
            IF (iliin.GE.1) ilblck2 = 0
            ilacll  =  0
          ENDIF
        ENDDO
          
        nstsi = nstsi + 1
        WRITE(fp2,'(A,I2,A)') '* RADIAL FLUXSURFACE ',1,
     .                        ' TRANSPARENT (PFZ, LEFT TARGET)'
        WRITE(fp2,'(5I6)') nstsi,1,iras-irtrap,
     .                     virloc(irtrap+1,IKLO),ikto+1
        WRITE(fp2,'(2I6,A6,4I6,2I3.3)') iliin,ilside,ilswch2,0,0,6,0,
     .                                  ilblck2,ilacll

        IF (cneur.EQ.7.AND.ctrap.EQ.7) THEN
          WRITE(fp2,'(A)') '     1     0     0     0'
          WRITE(fp2,'(1P,2E12.4)') material(eirmat2),eirtemp2
          WRITE(fp2,'(1P,6E12.4)') defrecycf,defrecyct,0.0,1.0,0.5,1.0
        ENDIF

        nstsi = nstsi + 1
        WRITE(fp2,'(A,I2,A)') '* RADIAL FLUXSURFACE ',1,
     .                        ' TRANSPARENT (PFZ, RIGHT TARGET)'
        IF (grdnmod.NE.0) THEN
          WRITE(fp2,'(5I6)') nstsi,1,iras-irtrap, 
     .                     cutpt2+1,virloc(irtrap+1,IKHI)+(ikti-ikto)+2
        ELSE
          WRITE(fp2,'(5I6)') nstsi,1,iras-irtrap,
     .                       ikti+2,virloc(irtrap+1,IKHI)+(ikti-ikto)+2
        ENDIF
        WRITE(fp2,'(2I6,A6,4I6,2I3.3)') iliin,ilside,ilswch2,0,0,6,0,
     .                                  ilblck2,ilacll

        IF (cneur.EQ.7.AND.ctrap.EQ.7) THEN
          WRITE(fp2,'(A)') '     1     0     0     0'
          WRITE(fp2,'(1P,2E12.4)') material(eirmat2),eirtemp2
          WRITE(fp2,'(1P,6E12.4)') defrecycf,defrecyct,0.0,1.0,0.5,1.0
        ENDIF

      ELSEIF ((eirbgk.EQ.3.OR.eirbgk.EQ.4).AND.vacregion(1).EQ.1) THEN
        nstsi = nstsi + 1
        WRITE(fp2,'(A,I2,A)') '* RADIAL FLUXSURFACE ',1,
     .                        ' TRANSPARENT (PFZ, LEFT TARGET, BGK)'
        WRITE(fp2,'(5I6)') nstsi,1,iras-irtrap,
     .                     virloc(irtrap+1,IKLO),ikto+1
c...DEV:
        WRITE(fp2,'(A,I3,A)') '    -1     0'//ilswch(1)//'     0'//
     .                        '     0     6     0',ilblck,'998'

        nstsi = nstsi + 1

        WRITE(fp2,'(A,I2,A)') '* RADIAL FLUXSURFACE ',1,
     .                        ' TRANSPARENT (PFZ, RIGHT TARGET, BGK)'
        IF (grdnmod.NE.0) THEN
          WRITE(fp2,'(5I6)') nstsi,1,iras-irtrap, 
     .                     cutpt2+1,virloc(irtrap+1,IKHI)+(ikti-ikto)+2
        ELSE
          WRITE(fp2,'(5I6)') nstsi,1,iras-irtrap,
     .                       ikti+2,virloc(irtrap+1,IKHI)+(ikti-ikto)+2
        ENDIF
c        WRITE(fp2,'(5I6)') nstsi,1,iras-irtrap,
c     .                     ikti+2,virloc(irtrap+1,IKHI)+(ikti-ikto)+2
c...DEV:
        WRITE(fp2,'(A,I3,A)') '    -1     0'//ilswch(1)//'     0'//
     .                        '     0     6     0',ilblck,'998'

      ELSEIF (eirbgk.EQ.2) THEN 

c...    Write non default surfaces to the EIRENE input file:

        IF (stopopt2.NE.9) THEN
c...      PFZ:

          DO i1 = 1, nsur(1)
            nstsi = nstsi + 1
            WRITE(fp2,'(A,I2,A,I2,A)') '* RADIAL FLUXSURFACE ',1,
     .                                 ' TRANSPARENT (PFZ ',i1,')'
            IF (ssur(i1,1).LE.ikto) THEN
              WRITE(fp2,'(7I6)') nstsi,1,iras-irtrap,1,1,
     .                           ssur(i1,1),esur(i1,1)+1
            ELSE
              WRITE(fp2,'(7I6)') nstsi,1,iras-irtrap,1,1,
     .                           ssur(i1,1)+(ikti-ikto)+1,
     .                           esur(i1,1)+(ikti-ikto)+1+1
            ENDIF
            WRITE(fp2,'(A,I3,I3.3)') '    -1     0'//ilswch(1)//
     .        '     0     0     6     0',ilblck,
     .                            csur(i1,1)+eirnpgdat+1
          ENDDO

        ENDIF

      ENDIF


      IF ((eirbgk.EQ.3.OR.eirbgk.EQ.4).AND.vacregion(2).EQ.1) THEN

          IF (irbreak.LE.irwall) THEN
            ik1 = ikbreak(irbreak-1) + 3
            ik2 = virloc(irbreak-1,IKHI) + 3
            ir1 = irbreak-1
          ELSE
            ik1 = ikbreak(irbreak-1) + 2 + nks(1)
            ik2 = virloc(irbreak-1,IKHI) + 2 + nks(1)
            ir1 = irbreak - irtrap
          ENDIF
          
          nstsi = nstsi + 1
          WRITE(fp2,'(A,I2,A,I2,A)') '* RADIAL FLUXSURFACE ',1,
     .                               ' TRANSPARENT (PLENUM ',i1,')'
          WRITE(fp2,'(7I6)') nstsi,1,ir1,1,1,ik1,ik2
          WRITE(fp2,'(A,I3,A)') '    -1     0'//ilswch(1)//'     0'//
     .                          '     0     6     0',ilblck,'998'

      ELSEIF (eirbgk.EQ.2.AND.(stopopt2.EQ.9.OR.stopopt2.EQ.10)) THEN
c...    Plenum:
       
        IF (stopopt2.EQ.10) THEN
          i3 = 2
        ELSE
          i3 = 1
        ENDIF
       
        DO i1 = 1, nsur(i3)
          nstsi = nstsi + 1
          WRITE(fp2,'(A,I2,A,I2,A)') '* RADIAL FLUXSURFACE ',1,
     .                               ' TRANSPARENT (PLENUM ',i1,')'
          WRITE(fp2,'(7I6)') nstsi,1,irbreak-1,1,1,
     .                       ssur(i1,i3)+2,
     .                       esur(i1,i3)+2+1
          WRITE(fp2,'(A,I3,I3.3)') '    -1     0'//ilswch(1)//'     0'//
     .                             '     0     6     0',ilblck,
     .                             csur(i1,i3)+eirnpgdat+1
        ENDDO

      ENDIF


      IF (grdnmod.NE.0) THEN

c...    Scan the connection map for IR=IRWALL and add a radial EIRENE surface for
c       each broken grid boundary region:

c...    SHOULD BE ABLE TO MAKE THIS GENERIC!  NO NEED FOR DEDICATED NON-BROKEN CODE.

        nlist   = 0
        ik1list = 0
        ik2list = 0
        irlist  = 0

        ir = irwall

        DO ik = 1, nks(ir)

          ik1 = ikins(ik,ir)
          ir1 = irins(ik,ir)

c...      Ignore links to virtual cells:
          IF (virtag(ik1,ir1).EQ.1) CYCLE

c...      Check if this boundary region has already been identified,
c         and if not then add it to the list:
          load = .TRUE.
          DO i1 = 1, nlist
            IF (irlist(i1).EQ.ir1) load = .FALSE.
          ENDDO
          IF (load) THEN
            nlist = nlist + 1  
            irlist(nlist) = ir1
          ENDIF
 	  
c...      Update the IK indecies of the boundary region:
          DO i1 = 1, nlist
            IF (irlist(i1).EQ.ir1) THEN
              IF (ik1list(i1).EQ.0) THEN
                ik1list(i1) = ik1
              ELSE
                ik2list(i1) = ik1                
              ENDIF
            ENDIF
          ENDDO

        ENDDO

        ilswch2 = ilswch(2)
        ilswch2(2:2) = '3'

        DO i1 = 1, nlist
          WRITE(0,*) 'REGIONS:',irtrap,cutpt2-1
          WRITE(0,*) 'REGIONS:',i1,irlist(i1),ik1list(i1),ik2list(i1)

c...      Adjust indecies based on their relation to the grid cut points and the
c         separatrix:

          IF (irlist(i1).GT.irtrap) THEN
           IF (ik1list(i1).GT.ikto2(irlist(i1)))
     .       ik1list(i1)=ik1list(i1)+nks(2)+1
           IF (ik2list(i1).GT.ikto2(irlist(i1)))
     .       ik2list(i1)=ik2list(i1)+nks(2)+1
c           IF (ik1list(i1).GT.cutpt2-1) ik1list(i1)=ik1list(i1)+nks(2)+1
c           IF (ik2list(i1).GT.cutpt2-1) ik2list(i1)=ik2list(i1)+nks(2)+1
c            IF (ik1list(i1).GT.ikto) ik1list(i1) = ik1list(i1) +nks(1)+1
c            IF (ik2list(i1).GT.ikto) ik2list(i1) = ik2list(i1) +nks(1)+1
            irlist(i1) = irlist(i1) - irtrap + 1
          ELSE
            IF (ik1list(i1).GT.cutpt1-1) ik1list(i1) = ik1list(i1) + 1
            IF (ik1list(i1).GT.cutpt2-1) ik1list(i1) = ik1list(i1) + 1
            IF (ik2list(i1).GT.cutpt1-1) ik2list(i1) = ik2list(i1) + 1
            IF (ik2list(i1).GT.cutpt2-1) ik2list(i1) = ik2list(i1) + 1
c            IF (ik1list(i1).GT.ikto) ik1list(i1) = ik1list(i1) + 1
c            IF (ik1list(i1).GT.ikti) ik1list(i1) = ik1list(i1) + 1
c            IF (ik2list(i1).GT.ikto) ik2list(i1) = ik2list(i1) + 1
c            IF (ik2list(i1).GT.ikti) ik2list(i1) = ik2list(i1) + 1
          ENDIF

          WRITE(0,*) 'REGIONS:',i1,irlist(i1),ik1list(i1),ik2list(i1)

          nstsi = nstsi + 1
          WRITE(fp2,'(A,I2,A)') '* RADIAL FLUXSURFACE ',ir,
     .                          ' TRANSPARENT (WALL)'
          WRITE(fp2,1001) nstsi,1,irlist(i1),ik1list(i1),ik2list(i1)+1
          IF (eirbgk.eq.3.AND.i1.EQ.2) THEN
            WRITE(fp2,'(A,I3,A)') '    -1     0'//ilswch2//
     .                           '     0     0     6     0',ilblck,'998'
          ELSE
            WRITE(fp2,'(A,I3,A)') '    -1     0'//ilswch2//
     .                           '     0     0     6     0',ilblck,'001'
          ENDIF

1001      FORMAT(10(I6:))
        ENDDO

c        IF (mode.EQ.1) STOP 'sdfsd'

      ELSEIF (nbr.GT.0.AND.stopopt2.NE.9.AND.stopopt2.NE.10.AND.
     .    .NOT.((eirbgk.EQ.3.OR.eirbgk.EQ.4).AND.
     .          vacregion(2).EQ.1)) THEN
c     .         (vacregion.EQ.2.OR.vacregion.EQ.3))) THEN
c...much better way to do this: scan IRWALL (since it goes all the way around the
c   SOL and use IKINS to identify different
c   regions that should have independent surfaces -- also goes nicely with
c

c
c       Partial radial surfaces due to broken rings:
c
        WRITE(0,*) 'DETERMINING PARTIAL RADIAL SURFACE'

c PFZ BREAK
        DO ir = irsep, nrs
          IF (ir.EQ.irwall-1.OR.idring(ir).EQ.-1) CYCLE

          ik1 = 0
          ik2 = 0

          IF (stopopt.EQ.121) THEN

            DO ik = virloc(ir,IKLO), virloc(ir,IKHI)
              IF (irouts(ik,ir).NE.irwall) ik1 = ik
              IF (ik.EQ.virloc(ir,IKHI).AND.ik.NE.ik1) THEN
                ik2 = ik - 1

                IF (ir.GT.irtrap) THEN
                  IF (ik1.GT.ikto) ik1 = ik1 + nks(1) + 1
                  IF (ik2.GT.ikto) ik2 = ik2 + nks(1) + 1
                  ir1 = ir - irtrap + 1
                ELSE
                  IF (ik1.GT.ikto) ik1 = ik1 + 1
                  IF (ik1.GT.ikti) ik1 = ik1 + 1
                  IF (ik2.GT.ikto) ik2 = ik2 + 1
                  IF (ik2.GT.ikti) ik2 = ik2 + 1
                  ir1 = ir
                ENDIF
 
                nstsi = nstsi + 1

                WRITE(fp2,'(A,I2,A)') '* RADIAL FLUXSURFACE ',ir,
     .                                ' TRANSPARENT (WALL)'
                WRITE(fp2,'(5I6)') nstsi,1,ir1,ik1+1,ik2+2
                WRITE(fp2,'(A,I3,A)') '    -1     0'//ilswch(2)//
     .                     '     0     0     6     0',ilblck,'001'
                ik1 = 0
                ik2 = 0
              ENDIF
            ENDDO

          ELSE

            DO ik = virloc(ir,IKLO), virloc(ir,IKHI)
              IF (irouts(ik,ir).EQ.irwall.AND.ik1.EQ.0) ik1 = ik
              IF (irouts(ik,ir).NE.irwall.AND.ik1.NE.0) THEN
                ik2 = ik - 1
                IF (ik1.GT.ikto) ik1 = ik1 + 1
                IF (ik1.GT.ikti) ik1 = ik1 + 1
                IF (ik2.GT.ikto) ik2 = ik2 + 1
                IF (ik2.GT.ikti) ik2 = ik2 + 1
                nstsi = nstsi + 1
                WRITE(fp2,'(A,I2,A)') '* RADIAL FLUXSURFACE ',ir,
     .                                ' TRANSPARENT (WALL)'
                WRITE(fp2,'(5I6)') nstsi,1,ir,ik1,ik2+1
                WRITE(fp2,'(A,I3,A)') '    -1     0'//ilswch(2)//
     .                           '     0     0     6     0',ilblck,'001'
                ik1 = 0
                ik2 = 0
              ENDIF
            ENDDO

          ENDIF
        ENDDO
      ENDIF

c...  Search surface properties array to see if the properties of the outer 
c     surface of the grid are being modified:

c...  Assign defaults for outer SOL surface of grid:
      IF (cneur.EQ.7.AND.ctrap.EQ.7) THEN
        iliin   =  1
        ilside  =  0
        ilswch2 = '     0'
        ilblck2 =  0
        ilacll  =  0
      ELSE
        iliin   = -1
        ilside  =  0
        ilswch2 = ilswch(2)
        ilblck2 = ilblck
        ilacll  =  1
      ENDIF

      DO i1 = 1, eirnspdat
        IF (eirspdat(i1,1).EQ.3.0.AND.eirspdat(i1,2).EQ.1.0) THEN
          iliin   =  NINT(eirspdat(i1,3))
          ilside  =  NINT(eirspdat(i1,4))
          ilswch2 =  '     0'
          IF (iliin.GE.1) ilblck2 = 0
          ilacll  =  0
        ENDIF
      ENDDO

      IF ((eirbgk.EQ.3.OR.eirbgk.EQ.4).AND.vacregion(4).EQ.1) THEN
        nstsi = nstsi + 1 
        WRITE(fp2,'(A,I2,A)') '* RADIAL FLUXSURFACE ',irwall-1,
     .                        ', TRANSPARENT (OUTER WALL) SEARCHING'
        WRITE(fp2,'(5I6)') nstsi,1,irwall-1,virloc(irwall-1,IKLO),
     .                     virloc(irwall-1,IKHI)+3
        WRITE(fp2,'(2I6,A6,4I6,2I3)') iliin,ilside,ilswch2,0,0,6,0,
     .                                ilblck2,997

      ELSEIF (grdnmod.EQ.0) THEN
        nstsi = nstsi + 1
        WRITE(fp2,'(A,I2,A)') '* RADIAL FLUXSURFACE ',irwall-1,
     .                        ', TRANSPARENT (OUTER WALL)'
        WRITE(fp2,'(5I6)') nstsi,1,irwall-1,virloc(irwall-1,IKLO),
     .                     virloc(irwall-1,IKHI)+3
        WRITE(fp2,'(2I6,A6,4I6,2I3.3)') iliin,ilside,ilswch2,0,0,6,0,
     .                                  ilblck2,ilacll
      ENDIF

      IF (cneur.EQ.7.AND.ctrap.EQ.7) THEN
        WRITE(fp2,'(A)') '     1     0     0     0'
        WRITE(fp2,'(1P,2E12.4)') material(eirmat2),eirtemp2
        WRITE(fp2,'(1P,6E12.4)') defrecycf,defrecyct,0.0,1.0,0.5,1.0
      ENDIF

c...  TARGETS:

c...  Determine listing of "standard" poloidal surfaces, to account
c     for "specified" surface properties a la EIRSPDAT:

      nstsi = nstsi + 1
      WRITE(fp2,'(A,I2,A)') '* POLOIDAL FLUXSURFACE ',1,
     .                    ', SOURCE/REFLECTING (LEFT TARGET)'
      WRITE(fp2,'(5I6)') nstsi,2,1,1,irwall-1
      WRITE(fp2,'(6I6)') 1,0,0,0,0,6
      WRITE(fp2,'(2I6)') 1,2
      WRITE(fp2,'(1P,2E12.4)') material(eirmat1),eirtemp1
      WRITE(fp2,'(1P,6E12.4)') defrecycf,defrecyct,0.0,1.0,0.5,1.0


      IF (eirnspdat.GT.0) THEN

        DO i1 = 1, ntar2
          nstsi = nstsi + 1
          WRITE(fp2,'(2(A,I2))') '* POLOIDAL FLUXSURFACE ',nks(irsep)+3,
     .                           ', SOURCE (HIGH INDEX TARGET) ',i1
          WRITE(fp2,'(5I6)') nstsi,2,nks(irsep)+3,MAX(1,btar2(i1)-1),
     .                       etar2(i1)
          IF (ttar2(i1).GT.0) THEN
            i2 = ttar2(i1)
            IF (vacregion(8).EQ.1) THEN
c...This is completely lame, as is this entire vaccum grid region.  The assumption
c   is that if the hard-coded vacuum grid is being used, that it is 
c   behind the transparent targets, which will not always be the case.  It is also
c   impossible to have transparent targets that do not exit into the vacuum grid (for this
c   hard-coded option).
              IF (eirnsdtor.GT.1) THEN
               WRITE(fp2,'(8I6)') (INT(eirspdat(i2,i3)),i3=3,5),0,0,6,0,
     .                            999998
              ELSE
               WRITE(fp2,'(8I6)') (INT(eirspdat(i2,i3)),i3=3,5),0,0,6,0,
     .                            001998
              ENDIF
            ELSE
c...Sloppy Joe here -- specifying ILCELL for all of these surface types...
c...bug: of sorts I think.  This switching surface was bein activated ...
              IF (eirntrans.GT.0.AND.iflexopt(8).EQ.10) THEN
c              IF (eirnsdtor.GT.0.AND.iflexopt(8).EQ.10) THEN
c...LAME! set of conditions, but they will do for now:
               WRITE(fp2,'(8I6)') (INT(eirspdat(i2,i3)),i3=3,5),0,0,6,0,
     .                            999997
              ELSE
               WRITE(fp2,'(8I6)') (INT(eirspdat(i2,i3)),i3=3,5),0,0,6,0,
     .                            001001
              ENDIF
            ENDIF

            WRITE(fp2,'(2I6)') 1,2
            WRITE(fp2,'(1P,5E12.4)') material(eirmat1),eirtemp1,0.0,
     .                               eirspdat(i2,6),eirspdat(i2,7)
            WRITE(fp2,'(1P,6E12.4)') defrecycf,defrecyct,0.0,1.0,0.5,1.0

          ELSE

            WRITE(fp2,'(6I6)') 1,0,0,0,0,6
            WRITE(fp2,'(2I6)') 1,2
            WRITE(fp2,'(1P,2E12.4)') material(eirmat1),eirtemp1
            WRITE(fp2,'(1P,6E12.4)') defrecycf,defrecyct,0.0,1.0,0.5,1.0

          ENDIF

        ENDDO

      ELSE
        nstsi = nstsi + 1
        WRITE(fp2,'(A,I2,A)') '* POLOIDAL FLUXSURFACE ',nks(irsep)+3,
     .                        ', SOURCE/REFLECTING (RIGHT TARGET)'
        WRITE(fp2,'(5I6)') nstsi,2,nks(irsep)+3,1,irwall-1
        WRITE(fp2,'(6I6)') 1,0,0,0,0,6
        WRITE(fp2,'(2I6)') 1,2
        WRITE(fp2,'(1P,2E12.4)') material(eirmat1),eirtemp1
        WRITE(fp2,'(1P,6E12.4)') defrecycf,defrecyct,0.0,1.0,0.5,1.0
c        WRITE(fp2,'(2A)')
c     .    '  1.0000E 00  1.0000E 00  0.0000E 00  1.0000E 00',
c     .    '  0.5000E 00  1.0000E 00'
      ENDIF


c...  Radial surfaces specified in the DIVIMP input file:
      DO i1 = 1, eirnasdat
        IF (eirasdat(i1,1).NE.3.0) CYCLE

        WRITE(PINOUT,*) 'FOUND NON-STANDARD RADIAL SURFACE'

        IF (eirasdat(i1,2).EQ.5.0) THEN
c...      Radial surface on the separatrix:
          ring1 = irsep-1
          cell = virloc(irsep,IKLO)
          cell2 = virloc(irsep,IKHI) + 3
          IF (eirasdat(i1,3).NE.99.0) cell = NINT(eirasdat(i1,3))
          IF (eirasdat(i1,4).NE.99.0) cell2 = NINT(eirasdat(i1,4))
        ELSEIF (eirasdat(i1,2).GT.irsep.AND.eirasdat(i1,2).LE.nrs) THEN
c...THIS BLOWS!
          IF (NINT(eirasdat(i1,2)).GT.irwall) THEN
            ring1 = NINT(eirasdat(i1,2)) - irtrap 
          ELSE
            ring1 = NINT(eirasdat(i1,2)) - 1
          ENDIF  
          cell = virloc(NINT(eirasdat(i1,2)),IKLO)
          cell2 = virloc(NINT(eirasdat(i1,2)),IKHI) + 3
          IF (eirasdat(i1,3).NE.99.0) cell = NINT(eirasdat(i1,3))
          IF (eirasdat(i1,4).NE.99.0) cell2 = NINT(eirasdat(i1,4))
        ELSE
          CALL ER('Block03a','Unknown radial surface specification',*99)
        ENDIF 

c...    Absorbing surface by default:
        iliin   = 2
        IF (iflexopt(7).EQ.22) THEN
          iliin   = 3
          WRITE(0,*) 'ASSIGNING REFLECTING RADIAL SURFACE'
        ENDIF
        ilside  = 0
        ilswch3 = 0
        ilcell  = 0

c...    Seach for non-default surface properties:


c...    Write surface to EIRENE input file:
        nstsi = nstsi + 1 
        WRITE(fp2,'(A,I2,A)') '* RADIAL FLUXSURFACE (SEPARATRIX)'
        WRITE(fp2,90) nstsi,1     ,ring1  ,cell,cell2
        WRITE(fp2,90) iliin,ilside,ilswch3,0    ,0    ,6,0,ilcell

90      FORMAT(10(I6:))

      ENDDO


c...  Poloidal surfaces specified in the DIVIMP input file:
      DO i1 = 1, eirnasdat
        IF (eirasdat(i1,1).NE.4.0) CYCLE

        WRITE(PINOUT,*) 'FOUND NON-STANDARD POLOIDAL SURFACE'

        cell = NINT(eirasdat(i1,2))
        ring1 = 1
        ring2 = irwall-1
        IF (eirasdat(i1,3).NE.99.0) ring1 = NINT(eirasdat(i1,3))
        IF (eirasdat(i1,4).NE.99.0) ring2 = NINT(eirasdat(i1,4))

c...    Absorbing surface by default:
        iliin   = 2
        IF (iflexopt(7).EQ.22) THEN
          iliin   = 3
          WRITE(0,*) 'ASSIGNING REFLECTING POLOIDAL SURFACE'
        ENDIF
        ilside  = 0
        ilswch3 = 0
        ilcell  = 0

c...    Seach for non-default surface properties:


c...    Write surface to EIRENE input file:
        nstsi = nstsi + 1 
        WRITE(fp2,'(A,I2,A)') '* POLOIDAL FLUXSURFACE'
        WRITE(fp2,90) nstsi,2     ,cell  ,ring1,ring2
        WRITE(fp2,90) iliin,ilside,ilswch3,0    ,0    ,6,0,ilcell
      ENDDO

c...  Advance the file pointer in the template:
      IF (mode.EQ.1) THEN
30      CALL ReadLine(fp1,buffer,1,*97,*98)
        IF (buffer(1:3).NE.'***') GOTO 30
        BACKSPACE fp1
      ENDIF

      RETURN
97    CALL ER('WriteBlock03a','Unexpected end of file'        ,*99)
98    CALL ER('WriteBlock03a','Problems reading template file',*99)
99    STOP 'END OF WRITEINPUT03A'
      END
c
c
c ======================================================================
c
c subroutine: WriteBlock04
c
c
c
c
      SUBROUTINE WriteBlock04(fp1,fp2)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER    MAXSTRAT,MAXSDATA
      PARAMETER (MAXSTRAT=10,MAXSDATA=10)
      COMMON /EIRCOM/
     .        eir_07nsrc,eir_07opt ,eir_07stra,
     .        eir_07ind1,
     .        eir_07ind2,eir_07ind3,
     .        eir_07wght           ,eir_07data
      INTEGER eir_07nsrc,eir_07opt ,eir_07stra(MAXSTRAT),
     .        eir_07ind1(MAXSTRAT),
     .        eir_07ind2(MAXSTRAT),eir_07ind3(MAXSTRAT)
      REAL    eir_07wght          ,eir_07data(MAXSTRAT,MAXSDATA)

      INTEGER     fp1,fp2,nrcm,nrca,nrci,nmass
      LOGICAL     new
      CHARACTER   buffer*200
      CHARACTER*8 psym(2,3)

      nmass = NINT(crmb)

      IF (nmass.NE.1.AND.nmass.NE.2) 
     .  CALL ER('WriteBlock04','EIRENE can only run with H and D',*99)

      psym(1,1) = 'H'
      psym(2,1) = 'D'
      psym(1,2) = 'H2'
      psym(2,2) = 'D2'
      psym(1,3) = 'H2+'
      psym(2,3) = 'D2+'

c...  Set defaults:
      nrca = 1
      nrcm = 0
      nrci = 0

      WRITE(fp2,'(A)') '*** 4. DATA FOR SPECIES SPECIFICATION AND '//
     .                 'ATOMIC PHYSICS MODULE (DIVIMP)'

      IF (eiracx .GE.1) nrca = nrca + 1
      IF (eirbgk .GE.1) nrca = nrca + 2

      IF (eird2ion.EQ.1) nrcm = nrcm + 1
      IF (eird2dis.EQ.1) nrcm = nrcm + 2
      IF (eirph2  .EQ.1) nrcm = nrcm + 1
      IF (eircxd2 .EQ.2) nrcm = nrcm + 1
      IF (eirbgk  .GE.1) nrcm = nrcm + 2
      
      IF (eirfuji.EQ.1) nrci = nrci + 3

      WRITE(fp2,90) '* ATOMIC REACTION CARDS  NREACI='
      WRITE(fp2,90) '    22'
      WRITE(fp2,90) '  1 AMJUEL H.4 2.1.5     EI  0  1'
      WRITE(fp2,90) '  2 AMJUEL H.102.1.5     EI  0  1'
      WRITE(fp2,90) '  3 HYDHEL H.1 3.1.8     CX  1  1'
      WRITE(fp2,90) '  3 HYDHEL H.3 3.1.8     CX  1  1'
      WRITE(fp2,90) '  4 AMJUEL H.4 2.2.9     EI  0  2'
      WRITE(fp2,90) '  5 AMJUEL H.4 2.2.5     DS  0  2'
      WRITE(fp2,90) '  6 AMJUEL H.4 2.2.10    DS  0  2'
      WRITE(fp2,90) '  7 AMJUEL H.4 2.2.12    DS  0  2'
      WRITE(fp2,90) '  8 AMJUEL H.4 2.2.11    DS  0  2'
      WRITE(fp2,90) '  9 AMJUEL H.4 2.2.14    DS  0  2'
      WRITE(fp2,90) ' 10 AMJUEL H.8 2.2.14    DS  0  2'
      IF     (pincode.EQ.2) THEN
        WRITE(fp2,90) ' 11 AMJUEL H.1 0.2       EL  1  4'
        WRITE(fp2,90) ' 12 AMJUEL H.1 0.3       EL  1  2'
        WRITE(fp2,90) ' 12 AMJUEL H.3 0.3       EL  1  2'
      ELSEIF (pincode.EQ.3) THEN
        WRITE(fp2,90) ' 11 AMJUEL H.0 0.2T      EL  1  4'
        WRITE(fp2,90) ' 11 AMJUEL H.1 0.2T      EL  1  4'
        WRITE(fp2,90) ' 12 AMJUEL H.0 0.3T      EL  1  2'
        WRITE(fp2,90) ' 12 AMJUEL H.1 0.3T      EL  1  2'
        WRITE(fp2,90) ' 12 AMJUEL H.3 0.3T      EL  1  2'
      ELSE
        CALL ER('BLOCK 4:','INVALID PINCODE VALUE',*99)
      ENDIF
      WRITE(fp2,90) ' 13 HYDHEL H.2 2.3.9     EI  0  4'
      WRITE(fp2,90) ' 14 METHAN H.2 2.23      EI  0 12'
      WRITE(fp2,90) ' 15 AMJUEL H.4 2.1.8     RC  0  1'
      WRITE(fp2,90) ' 16 AMJUEL H.102.1.8     RC  0  1  1.3600E 01'
      WRITE(fp2,90) ' 17 CONST  H.2           EL  2  2'
      WRITE(fp2,90) ' -2.1091E 01  0.2500E 00  0.0000E 00  0.0000E '//
     .              '00  0.0000E 00  0.0000E 00'
      WRITE(fp2,90) '  0.0000E 00  0.0000E 00  0.0000E 00  0.0000E '//
     .              '00  0.0000E 00  0.0000E 00'
      WRITE(fp2,90) ' 18 CONST  H.2           EL  2  2'
      WRITE(fp2,90) ' -2.0589E 01  0.2500E 00  0.0000E 00  0.0000E '//
     .              '00  0.0000E 00  0.0000E 00'
      WRITE(fp2,90) '  0.0000E 00  0.0000E 00  0.0000E 00  0.0000E '//
     .              '00  0.0000E 00  0.0000E 00'
      WRITE(fp2,90) ' 19 CONST  H.2           EL  4  4'
      WRITE(fp2,90) ' -2.0357E 01  0.2500E 00  0.0000E 00  0.0000E '//
     .              '00  0.0000E 00  0.0000E 00'
      WRITE(fp2,90) '  0.0000E 00  0.0000E 00  0.0000E 00  0.0000E '//
     .              '00  0.0000E 00  0.0000E 00'
      WRITE(fp2,90) ' 20 AMJUEL H.3 3.2.3     CX  1  2'
      WRITE(fp2,90) ' 21 AMJUEL H.9 3.1.8     CX  1  1'
      WRITE(fp2,90) ' 22 AMJUEL H.2 3.1.8FJ   CX  1  1'
      WRITE(fp2,90) '*NEUTRAL ATOMS SPECIES CARDS: NATMI='
      WRITE(fp2,90) '     2'
      WRITE(fp2,93) 1,psym(nmass,1),nmass,1,1,0,1,-4,0,nrca
      WRITE(fp2,90) '     1   115   114     0 30000   000'
      WRITE(fp2,92) 2.0,0.0,0.0,0.0,eirscale(1)
      IF (eiracx.EQ.1) THEN
        WRITE(fp2,90) '     3   114   111   114 01001'
        WRITE(fp2,92) 0.0,0.0,0.0,0.0,eirscale(2)
      ENDIF
      IF (eirbgk.GE.1) THEN
        WRITE(fp2,90) '    17   314           0 01001   000   111'
        WRITE(fp2,92) 0.0,0.0,0.0,0.0,eirscale(3)
        WRITE(fp2,90) '    19   514           0 01001   000   112'
        WRITE(fp2,92) 0.0,0.0,0.0,0.0,eirscale(4)
      ENDIF
      WRITE(fp2,90) ' 2 C        12  6  1  0  2  2  0  1'
      WRITE(fp2,90) '    14   115   214     0 00000'
      WRITE(fp2,90) '  1.1300E 01  0.0000E 00'
      WRITE(fp2,90) '** 4b NEUTRAL MOLECULES SPECIES CARDS: NMOLI='
      WRITE(fp2,90) '     1'
      IF     (pincode.EQ.2) THEN
        WRITE(fp2,93) 1,psym(nmass,2),2*nmass,2,2,0,0,-1,0,nrcm
      ELSEIF (pincode.EQ.3) THEN
        WRITE(fp2,93) 1,psym(nmass,2),2*nmass,2,2,0,0, 2,0,nrcm
      ENDIF
      IF (eird2ion.EQ.1) THEN
        WRITE(fp2,90) '     4   115   113     0'
        WRITE(fp2,90) ' -1.5400E 01  0.0000E 00'
      ENDIF
      IF (eird2dis.EQ.1) THEN
        WRITE(fp2,90) '     5   115   121   000'
        WRITE(fp2,90) ' -1.0500E 01  0.0000E 00  3.0000E 00  3.0000E 00'
        WRITE(fp2,90) '     6   115   111   114'
        WRITE(fp2,90) ' -2.5000E 01  0.0000E 00  5.0000E 00  5.0000E 00'
      ENDIF
      IF (eirbgk.GE.1) THEN
        WRITE(fp2,90) '    18   414           0 01001   000   112'
        WRITE(fp2,92) 0.0,0.0,0.0,0.0,eirscale(9)
        WRITE(fp2,90) '    19   614           0 01001   000   111'
        WRITE(fp2,92) 0.0,0.0,0.0,0.0,eirscale(8)
      ENDIF
      IF (eircxd2.EQ.2) THEN
        WRITE(fp2,90) '    20   114   111   113 01001'
        WRITE(fp2,92) 0.0,0.0,0.0,0.0,eirscale(10)
      ENDIF
      IF (eirph2.EQ.1) THEN
        WRITE(fp2,90) '    12   114             01001'
        WRITE(fp2,92) 0.0,0.0,0.0,0.0,eirscale(12)
      ENDIF
      WRITE(fp2,90) '**4c TEST ION SPECIES CARDS:  NIONI ION SPECIE'//
     .              'S ARE CONSIDERED, NIONI='
      WRITE(fp2,90) '     1'
      WRITE(fp2,93) 1,psym(nmass,3),2*nmass,2,2,1,0,-4,0,nrci,-1
      IF (eirfuji.EQ.1) THEN
c...    Sawada and Fujimoto:
        WRITE(fp2,90) '     7   115   111   114                        '
        WRITE(fp2,90) ' -1.0400E 01  0.0000E 00  4.3000E 00  4.3000E 00'
        WRITE(fp2,90) '     8   115   124   000                        '
        WRITE(fp2,90) ' -1.5500E 01  0.0000E 00  0.2500E 00  0.2500E 00'
        WRITE(fp2,90) '     9   115   121   000 30002                  '
        WRITE(fp2,90) '  1.0000E 01  0.0000E 00  0.5000E 00  0.5000E 00'
      ELSE
c...    Janev (EIRENE default):
      ENDIF

c...  Advance to next block in template file:
40    CALL ReadLine(fp1,buffer,1,*97,*98)
      IF (buffer(1:3).NE.'***') GOTO 40
      BACKSPACE fp1

      RETURN
90    FORMAT(A)
91    FORMAT(A,2(I3:))
92    FORMAT(1P,6(E12.4:),0P)
93    FORMAT(I2,1X,A8,1X,12(I2,1X:))
97    CALL ER('WriteBlock04','Unexpected end of file'        ,*99)
98    CALL ER('WriteBlock04','Problems reading template file',*99)
99    STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock05
c
c
c
c
      SUBROUTINE WriteBlock05(fp1,fp2)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER    MAXSTRAT,MAXSDATA
      PARAMETER (MAXSTRAT=10,MAXSDATA=10)
      COMMON /EIRCOM/
     .        eir_07nsrc,eir_07opt ,eir_07stra,
     .        eir_07ind1,
     .        eir_07ind2,eir_07ind3,
     .        eir_07wght           ,eir_07data
      INTEGER eir_07nsrc,eir_07opt ,eir_07stra(MAXSTRAT),
     .        eir_07ind1(MAXSTRAT),
     .        eir_07ind2(MAXSTRAT),eir_07ind3(MAXSTRAT)
      REAL    eir_07wght          ,eir_07data(MAXSTRAT,MAXSDATA)

      INTEGER   fp1,fp2,nmass
      CHARACTER buffer*200
      CHARACTER*8 psym(2,5)

      nmass = NINT(crmb)

      IF (nmass.NE.1.AND.nmass.NE.2) 
     .  CALL ER('WriteBlock04','EIRENE can only run with H and D',*99)

      psym(1,1) = 'H+'
      psym(2,1) = 'D+'
      psym(1,2) = 'H(B)'
      psym(2,2) = 'D(B)'
      psym(1,3) = 'H2(B)'
      psym(2,3) = 'D2(B)'
      psym(1,4) = 'HH2(B)'
      psym(2,4) = 'DD2(B)'
      psym(1,5) = 'H2H(B)'
      psym(2,5) = 'D2D(B)'

      WRITE(fp2,90) '*** 5. DATA FOR PLASMA-BACKGROUND (DIVIMP)'

      IF     (eir_07opt.EQ.4.OR.eir_07opt.EQ.5) THEN
c...    BGK:
        WRITE(fp2,90) '*BULK ION SPECIES CARDS:  NPLSI ION SPECIES '//
     .                'ARE CONSIDERED, NPLSI='
        WRITE(fp2,90) '     6'
        WRITE(fp2,93) 1,psym(nmass,1),nmass,1,1,1,1,-1,2,1
        WRITE(fp2,90) '    15   115   111       30000'
        WRITE(fp2,91) 1.6E+01,0.0,0.0,0.0,eirsrcmul*eirscale(11)
        WRITE(fp2,90) ' 2 C+       12  6  1  1  1 -1  2  0'
        WRITE(fp2,93) 3,psym(nmass,2),  nmass,1,1,0,1,-1,0,0
        WRITE(fp2,93) 4,psym(nmass,3),2*nmass,2,2,0,1,-1,0,0
        WRITE(fp2,93) 5,psym(nmass,4),  nmass,1,1,0,1,-1,0,0
        WRITE(fp2,93) 6,psym(nmass,5),2*nmass,2,2,0,1,-1,0,0
        WRITE(fp2,90) '     7     7     7     7     6     6     0'//
     .                '     0     0     0     0     0'
      ELSEIF (eir_07opt.EQ.0.OR.eir_07opt.EQ.2) THEN
c...    Standard (no BGK):
        WRITE(fp2,90) '*BULK ION SPECIES CARDS:  NPLSI ION SPECIES '//
     .                'ARE CONSIDERED, NPLSI='
        WRITE(fp2,90) '     2'
        WRITE(fp2,93) 1,psym(nmass,1),nmass,1,1,1,1,-4,0,1
        WRITE(fp2,90) '    15   115   111       30000'
        WRITE(fp2,91) 1.6E+01,0.0,0.0,0.0,eirsrcmul*eirscale(11)
        WRITE(fp2,90) ' 1 C+       12  6  1  1  2  2  0  0'
        WRITE(fp2,90) '     6     6     6     6     6     6     0'//
     .                '     0     0     0     0     0'
      ELSE
        CALL ER('WriteBlock05','Invalid EIR_07OPT',*99)
      ENDIF

40    CALL ReadLine(fp1,buffer,1,*97,*98)
      IF (buffer(1:3).NE.'***') GOTO 40
      BACKSPACE fp1

      RETURN
90    FORMAT(A)
91    FORMAT(1P,5E12.4)
93    FORMAT(I2,1X,A8,1X,12(I2,1X:))
97    CALL ER('WriteBlock05','Unexpected end of file'        ,*99)
98    CALL ER('WriteBlock05','Problems reading template file',*99)
99    STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock06
c
c
c
c
      SUBROUTINE WriteBlock06(fp1,fp2)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER   fp1,fp2
      CHARACTER buffer*200

      WRITE(fp2,90) '*** 6. DATA FOR GENERAL REFLECTION MODEL (DIVIMP)'

      WRITE(fp2,90) 'TF'

      IF (eirtrim.EQ.1) THEN
c. *HARDCODED* Need to read the DIVIMP execution directory from an enviroment
c              variable:
        WRITE(fp2,90) 'PATH  ./TRIM/'
        WRITE(fp2,90) 'D_on_Mo'               
        WRITE(fp2,90) 'D_on_Fe'               
        WRITE(fp2,90) 'D_on_C'               
        WRITE(fp2,90) 'D_on_Be'               
      ENDIF

      WRITE(fp2,90) '  1.0000E 00  1.0000E 00'
      WRITE(fp2,90) '  0.2500E 00  0.2500E 00  0.5000E 00'
      WRITE(fp2,90) '  1.0000E 00  1.0000E 00'            
      WRITE(fp2,90) '  1.0000E 00  0.0000E 00  0.0000E 00'
      WRITE(fp2,91) eirermin,50.0,0.1,eirrinteg,eireinteg

c...  Advance input stream to the start of the next input block:
40    CALL ReadLine(fp1,buffer,1,*97,*98)
      IF (buffer(1:3).NE.'***') GOTO 40
      BACKSPACE fp1

      RETURN
90    FORMAT(A)
91    FORMAT(1P,5E12.4,0P)
97    CALL ER('WriteBlock06','Unexpected end of file'        ,*99)
98    CALL ER('WriteBlock06','Problems reading template file',*99)
99    STOP
      END
c
c
c ======================================================================
c
c subroutine: WriteBlock14
c
c
c
c
      SUBROUTINE WriteBlock14(fp1,fp2)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER fp1,fp2

      INTEGER    MAXSTRAT,MAXSDATA
      PARAMETER (MAXSTRAT=10,MAXSDATA=10)
      COMMON /EIRCOM/
     .        eir_07nsrc,eir_07opt ,eir_07stra,
     .        eir_07ind1,
     .        eir_07ind2,eir_07ind3,
     .        eir_07wght           ,eir_07data
      INTEGER eir_07nsrc,eir_07opt ,eir_07stra(MAXSTRAT),
     .        eir_07ind1(MAXSTRAT),
     .        eir_07ind2(MAXSTRAT),eir_07ind3(MAXSTRAT)
      REAL    eir_07wght           ,eir_07data(MAXSTRAT,MAXSDATA)




      INTEGER   i1,i2,nstrata
      CHARACTER buffer*200


      nstrata = 0
      DO i1 = 1, eirnstrata
        IF (eirstrata(i1,1).EQ.1.0.OR.eirstrata(i1,1).EQ.2.0) 
     .    nstrata = nstrata + 1
      ENDDO

      WRITE(fp2,'(A)') '*** 14. DATA FOR INTERFACING (DIVIMP)'

      IF     (eir_07opt.EQ.4) THEN

        WRITE(fp2,90) 'FTFFF'
        WRITE(fp2,90) '     6     0     1'
        WRITE(fp2,92) '     1     1  1.0000E 00',crmb
        WRITE(fp2,90) '     2     2  1.0000E 00 12.0000E 00'
        WRITE(fp2,90) '     3   -27'
        WRITE(fp2,90) '     4   -26'
        WRITE(fp2,90) '     5   -25'
        WRITE(fp2,90) '     6   -24'

        WRITE(fp2,91) nks(irsep)+2,irwall-2
        WRITE(fp2,91) nstrata
        WRITE(fp2,91) (1,i2=1,nstrata)

        DO i1 = 1, eirnstrata
          IF     (eirstrata(i1,1).EQ.1.0) THEN
            WRITE(fp2,'(6I6)') i1,eir_07ind1(i1)-1,-1,1,eir_07ind2(i1),
     .                                                  eir_07ind3(i1)
          ELSEIF (eirstrata(i1,1).EQ.2.0) THEN
            WRITE(fp2,'(6I6)') i1,eir_07ind1(i1)-1,+1,1,eir_07ind2(i1),
     .                                                  eir_07ind3(i1)
          ENDIF
        ENDDO

        WRITE(fp2,90) '  5.0000E 00  5.0000E 00  1.5000E 01'//
     .                '  1.0000e 01'
        WRITE(fp2,90) '     0'
        WRITE(fp2,90) '     0'
        WRITE(fp2,90) '*** link targets to surfaces'
        WRITE(fp2,90) '2   64'
        WRITE(fp2,90) '1   81'
        WRITE(fp2,90) '1   50'
        WRITE(fp2,90) '2   80'

        DO WHILE (.TRUE.)
          CALL ReadLine(fp1,buffer,1,*50,*98)
        ENDDO

      ELSEIF (eir_07opt.EQ.0) THEN

        WRITE(fp2,90) 'FTFFF'
        WRITE(fp2,90) '     2     0     1'
        WRITE(fp2,92) '     1     1  1.0000E 00',crmb
        WRITE(fp2,90) '     2     2  1.0000E 00 12.0000E 00'

        WRITE(fp2,91) nks(irsep)+2,irwall-2
        WRITE(fp2,91) nstrata
        WRITE(fp2,91) (1,i2=1,nstrata)

        DO i1 = 1, eirnstrata
          IF     (eirstrata(i1,1).EQ.1.0) THEN
            WRITE(fp2,'(6I6)') i1,eir_07ind1(i1)-1,-1,1,eir_07ind2(i1),
     .                                                  eir_07ind3(i1)
          ELSEIF (eirstrata(i1,1).EQ.2.0) THEN
            WRITE(fp2,'(6I6)') i1,eir_07ind1(i1)-1,+1,1,eir_07ind2(i1),
     .                                                  eir_07ind3(i1)
          ENDIF
        ENDDO

        WRITE(fp2,90) '  5.0000E 00  5.0000E 00  1.5000E 01'//
     .                '  1.0000e 01'
        WRITE(fp2,90) '     0'
        WRITE(fp2,90) '     0'
        WRITE(fp2,90) '*** link targets to surfaces'
        WRITE(fp2,90) '2   64'
        WRITE(fp2,90) '1   81'
        WRITE(fp2,90) '1   50'
        WRITE(fp2,90) '2   80'

        DO WHILE (.TRUE.)
          CALL ReadLine(fp1,buffer,1,*50,*98)
        ENDDO

      ELSEIF (eir_07opt.EQ.1.OR.eir_07opt.EQ.2.OR.eir_07opt.EQ.5) THEN

        STOP 'UNREASONABLE RESTRICTION ON VOLUME STRATA'

        IF (eir_07opt.EQ.5) THEN
c...      Passing 6 ions:
          CALL ReadLine(fp1,buffer,4,*97,*98)
          WRITE(fp2,90) 'FTFFF'
          WRITE(fp2,90) '     6     0     1'
          WRITE(fp2,92) '     1     1  1.0000E 00',crmb
          WRITE(fp2,90) '     2     2  1.0000E 00 12.0000E 00'
          WRITE(fp2,90) '     3   -27'
          WRITE(fp2,90) '     4   -26'
          WRITE(fp2,90) '     5   -25'
          WRITE(fp2,90) '     6   -24'
        ELSE
c...      Passing 2 ions (template set for D, not H):
          WRITE(0,*) 'WARNING: FLUID CODE PARICLE DATA PASSED FROM '//
     .               'EIRENE INPUT TEMPLATE FILE - MAY NOT BE '//
     .               'CONSISTENT WITH DIVIMP SETTINGS'
          CALL TransferLine(fp1,fp2,buffer,4)
        ENDIF

        WRITE(fp2,'(2I6)') nks(irsep)+2,irwall-2

        WRITE(fp2,'(I6     )') eir_07nsrc-2
        WRITE(fp2,'(20(I6:))') (1,i2=1,eir_07nsrc-2)

c ***********************************************************
c ARGH!  THIS BLOWS.  ONLY GOOD IF THERE ARE 2 VOLUME STRATA!
c ***********************************************************

        DO i1 = 1, eir_07nsrc-2
          IF     (eir_07stra(i1).EQ.1) THEN
            WRITE(fp2,'(6I6)') i1,eir_07ind1(i1)-1,-1,1,eir_07ind2(i1),
     .                                                  eir_07ind3(i1)
          ELSEIF (eir_07stra(i1).EQ.2) THEN
            WRITE(fp2,'(6I6)') i1,eir_07ind1(i1)-1,+1,1,eir_07ind2(i1),
     .                                                  eir_07ind3(i1)
          ENDIF
        ENDDO

        WRITE(fp2,'(A)') '  5.0000E 00  5.0000E 00  1.5000E 01'//
     .                   '  1.0000e 01'
        WRITE(fp2,'(A)') '     0'
        WRITE(fp2,'(A)') '     0'
        WRITE(fp2,'(A)') '*** link targets to surfaces'
        WRITE(fp2,'(A)') '2   64'
        WRITE(fp2,'(A)') '1   81'
        WRITE(fp2,'(A)') '1   50'
        WRITE(fp2,'(A)') '2   80'
c
c       Advance template file to EOF:
c
        DO WHILE (.TRUE.)
          CALL ReadLine(fp1,buffer,1,*50,*98)
        ENDDO

      ELSE
        CALL ER('WriteBlock07','Invalid option',*99)
      ENDIF

50    CONTINUE

      RETURN
90    FORMAT(A)
91    FORMAT(20(I6:))
92    FORMAT(A,1P,E12.4)
97    CALL ER('WriteBlock14','Unexpected end of file'        ,*99)
98    CALL ER('WriteBlock14','Problems reading template file',*99)
99    STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock07
c
c
c
c
      SUBROUTINE WriteBlock07(fp1,fp2)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER fp1,fp2

      REAL GetFlux,ATAN2C

      INTEGER ninitl

      INTEGER    MAXSTRAT,MAXSDATA
      PARAMETER (MAXSTRAT=10,MAXSDATA=10)
      COMMON /EIRCOM/
     .        eir_07nsrc,eir_07opt ,eir_07stra,
     .        eir_07ind1,
     .        eir_07ind2,eir_07ind3,
     .        eir_07wght           ,eir_07data
      INTEGER eir_07nsrc,eir_07opt ,eir_07stra(MAXSTRAT),
     .        eir_07ind1(MAXSTRAT),
     .        eir_07ind2(MAXSTRAT),eir_07ind3(MAXSTRAT)
      REAL    eir_07wght           ,eir_07data(MAXSTRAT,MAXSDATA)


c...  Move to SLCOM:
      COMMON  /EIRWALCOM/ walln,wallr,wallz,wallw,ebgki
      INTEGER walln,wallw(MAXPTS),ebgki
      REAL    wallr(MAXPTS,2),wallz(MAXPTS,2)

      INTEGER MAXNSR
c...  Has to be smaller than NSRFS in EIRENE (PARMUSR):
      PARAMETER (MAXNSR=200)

      INTEGER i1,i2,i3,ik,ir,in1,in2,in3,walli(MAXNSR),nsrfsi,
     .        nasor,npts,nemods,iktemp,irtemp,id,orientation,
     .        numseg(MAXNSR),ins,in4,in5,indseg(MAXNSR),
     .        cntseg(MAXNSR),idum1,count,iltor(MAXNSR),sftseg(MAXNSR),
     .        isec
      LOGICAL status
      REAL flux,x1(MAXNSR),x2(MAXNSR),y1(MAXNSR),y2(MAXNSR),
     .          z1(MAXNSR),z2(MAXNSR),totarea,brat,cost,
     .     zbnd1,zbnd2,sorlim,sorwgt(MAXNSR),sorifl,rad1,rad2,
     .     defzbnd1,defzbnd2,torlen,
     .     sorctx,sorcty,sorctz,len,dist,distmin,xcen,ycen,
     .     deltar,deltaz,alpha,beta,a1,a2,a3,b1,b2,b3,c1,c2,c3,
     .     te,ti,soreni,sorene,esheat,temp1,temp2,
     .     x(4),y(4),z(4),zshift,phiseg
      CHARACTER buffer*200


      IF     (piniseed.EQ.-1.OR. piniseed.EQ.0     ) THEN
        ninitl = -1
      ELSEIF (piniseed.GT. 0.AND.piniseed.LE.999999) THEN
        WRITE(PINOUT,*)
        WRITE(PINOUT,*) 'USING USER SPECIFIED RANDOM NUMBER SEED'
        WRITE(PINOUT,*)
        ninitl = piniseed
      ELSE
        CALL ER('WriteBlock07','Invalid PINISEED value',*99)
      ENDIF

      IF     (eirzaa.EQ.-1.0) THEN
        defzbnd1 = +0.1 + 0.0 
        defzbnd2 = -0.1 + 2.0 * PI * rxp * 100.0 * eirtorfrac
      ELSEIF (eirzaa.LT.0.0) THEN
        defzbnd1 = 0.0
        defzbnd2 = 0.0
      ELSE
        defzbnd1 = +0.1 + 0.0
        defzbnd2 = -0.1 + eirzaa * 100.0
      ENDIF

      zbnd1 = defzbnd1
      zbnd2 = defzbnd2

      WRITE(fp2,'(A)') '*** 7. PRIMARY NEUTRAL SOURCE (DIVIMP)'

      IF     (eir_07opt.EQ.0.OR.eir_07opt.EQ.4) THEN

        WRITE(fp2,84) eirnstrata+eirnpuff
        WRITE(fp2,84) 1,1,1
        WRITE(fp2,80) eiralloc

        DO i1 = 1, eirnstrata
          IF     (eirstrata(i1,1).EQ.1.0) THEN
            WRITE(fp2,81) '* DIVIMP SPECIFIED LOW INDEX TARGET'
            WRITE(fp2,81) 'FFFFF'
            WRITE(fp2,84)  INT(eir_07data(i1,1)),ninitl,3,1
            WRITE(fp2,80)  1.0
            WRITE(fp2,81) 'FFFTF'
            WRITE(fp2,84) -1
            WRITE(fp2,81) 'FFTFF'
            WRITE(fp2,84)  1
            WRITE(fp2,84)  1,4,0,0,0
            WRITE(fp2,80)  1.0,204.0,1.0,0.0,2000.0
            IF (eirnsdtor.GT.1) THEN
              WRITE(fp2,84) -1,0,0,999
            ELSE
              WRITE(fp2,84) -1,0,0,0
            ENDIF
            WRITE(fp2,80)  0.0,0.0,0.0,0.0,zbnd1,zbnd2
            WRITE(fp2,80)  3.0,0.5
            WRITE(fp2,80)  1.0,90.0

           ELSEIF (eirstrata(i1,1).EQ.2.0) THEN
            WRITE(fp2,81) '* DIVIMP SPECIFIED HIGH INDEX TARGET'
            WRITE(fp2,81) 'FFFFF'
            WRITE(fp2,84)  INT(eir_07data(i1,1)),ninitl,3,1
            WRITE(fp2,80)  1.0
            WRITE(fp2,81) 'FFFTF'
            WRITE(fp2,84) -1
            WRITE(fp2,81) 'FFTFF'
            WRITE(fp2,84)  1
            WRITE(fp2,84)  1,4,0,0,0
            WRITE(fp2,80)  1.0,204.0,1.0,0.0,1000.0
            IF (eirnsdtor.GT.1) THEN
              WRITE(fp2,84) -1,0,0,999
            ELSE
              WRITE(fp2,84) -1,0,0,0
            ENDIF
            WRITE(fp2,80)  0.0,0.0,0.0,0.0,zbnd1,zbnd2
            WRITE(fp2,80)  3.0,0.5
            WRITE(fp2,80)  1.0,90.0

          ELSEIF (eirstrata(i1,1).EQ.3.0) THEN
            WRITE(fp2,81) '* DIVIMP SPECIFIED VOLUME SOURCE'
            WRITE(fp2,81) 'FFFFF'
            WRITE(fp2,84)  INT(eir_07data(i1,1)),ninitl,3,1
            WRITE(fp2,80)  1.0
            WRITE(fp2,81) 'FFFTF'
            WRITE(fp2,84)  1
            WRITE(fp2,81) 'FFFTF'
            WRITE(fp2,84)  1
            WRITE(fp2,84)  1,4,0,eir_07ind1(i1),eir_07ind2(i1)+1,
     .                     1,nks(irsep)+3
            WRITE(fp2,80)  1.0,204.0,1.0,0.0,1000.0
            IF (eirnsdtor.GT.1) THEN
               WRITE(fp2,84) -1,0,0,999
            ELSE
              WRITE(fp2,84) -1,0,0,0
            ENDIF
            WRITE(fp2,80)  0.0,0.0,0.0,0.0,zbnd1,zbnd2
            WRITE(fp2,80)  3.0,0.5
            WRITE(fp2,80)  0.0,90.0
          ENDIF
        ENDDO




        IF (eirpmode.GT.0) THEN
c...DEV: Surface puffing
          IF (eirnpuff.GT.3) 
     .      CALL ER('WriteBlock07','Sorry, maximium of 3 puffing '//
     .                             'surfaces allowed at the moment',*99)

          DO i1 = 1, eirnpuff

            eir_07nsrc = eir_07nsrc + 1
c...        This isn't really correct.  If the stratum number is less than 
c           the number of stratum blocks in the EIRENE input file, then this
c           will not work (althoug the data will still be recored, it is just
c           that the data will not be stored sequentially):
            eir_07stra(eir_07nsrc) = eir_07nsrc

            flux   = 0.0
            nsrfsi = 1
            npts   = 999999
            nemods = 1
            IF     (eirpuff(i1,1).EQ.1.0) THEN
c...          Specified value:
              flux = eirpuff(i1,2) * ECH
            ELSEIF (eirpuff(i1,1).EQ.2.0) THEN
c...          Fraction of low index target flux:
              DO ir = irsep, nrs
                flux = flux - GetFlux(IKLO,ir)
              ENDDO
              flux = flux * eirpuff(i1,2) * ECH
            ELSEIF (eirpuff(i1,1).EQ.3.0) THEN
c...          Fraction of high index target flux:
              DO ir = irsep, nrs
                flux = flux + GetFlux(IKHI,ir)
              ENDDO
              flux = flux * eirpuff(i1,2) * ECH
            ELSEIF (eirpuff(i1,1).EQ.4.0) THEN
c...          Fraction of total target flux:
              DO ir = irsep, nrs
                flux = flux - GetFlux(IKLO,ir) + GetFlux(IKHI,ir)
              ENDDO
              flux = flux * eirpuff(i1,2) * ECH
            ELSEIF (eirpuff(i1,1).EQ.5.0) THEN

c..."realistic" energy distributions

c...          Range of "standard" EIRENE additional surfaces:
              IF (eirpuff(i1,3).GT.0.0.AND.eirpuff(i1,4).GT.0.0.AND.
     .            eirpuff(i1,3).LE.eirpuff(i1,4)) THEN

                IF (eirntorseg.NE.0) THEN

                  IF (NINT(eirpuff(i1,3)).LE.nvesm+nvesp) THEN
c...                Surface is a toroidally continuous additional surface:
                    nsrfsi = NINT(eirpuff(i1,4) - eirpuff(i1,3)) + 1
                  ELSE
c...                Surface is not toroidally continuous, and needs to
c                   be subdivided according to the number of toroidal segments
c                   in the toroidal geometry approximation:
c                    IF (eirnsection.GT.1) THEN
c                      WRITE(0,*) 'PUFFING FROM ADDITIONAL SURFACES '//
c     .                           'IN THE TOROIDAL APPROXIMATION '//
c     .                           'NEEDS WORK'
c                      STOP
c                    ELSE
c                      zshift = 0.0
c                    ENDIF

                    ins = nvesm + nvesp
                    nsrfsi = 0
                    DO in2 = eirpuff(i1,3), eirpuff(i1,4)
                      in3 = ins
c...                  Scan additional surfaces and find polygon data:
                      status = .FALSE.
                      DO in1 = 1, eirnasdat
                        IF (eirasdat(in1,1).EQ.7.0.OR.
     .                      eirasdat(in1,1).EQ.2.0) in3 = in3 + 1
                        IF (in3.EQ.in2) THEN
                          status = .TRUE.
                          EXIT
                        ENDIF
                      ENDDO
                      IF (.NOT.status) 
     .                  CALL ER('EIR07','Add surface not found',*99)                   

                      DO isec = 1, eirnsection                     

                        zshift = SNGL(DBLE(-eirzaa)*DBLE(eirtorfrac)/
     .                           DFLOAT(eirnsection)*DFLOAT(isec-1))

                        x(1) = eirasdat(in1,3)*100.0
                        y(1) = eirasdat(in1,4)*100.0
                        z(1) = eirasdat(in1,8)*100.0 + zshift*100.0
                        x(2) = eirasdat(in1,6)*100.0
                        y(2) = eirasdat(in1,7)*100.0
                        z(2) = eirasdat(in1,9)*100.0 + zshift*100.0
                        IF (eirasdat(MIN(eirnasdat,in1+1),1).EQ.-1.0) 
     .                    THEN
                          x(3) = eirasdat(in1+1,3)*100.0
                          y(3) = eirasdat(in1+1,4)*100.0
                          z(3) = eirasdat(in1+1,8)*100.0 + zshift*100.0
                          x(4) = eirasdat(in1+1,6)*100.0
                          y(4) = eirasdat(in1+1,7)*100.0
                          z(4) = eirasdat(in1+1,9)*100.0 + zshift*100.0
                        ELSE
                          CALL ER('Puff code','Only add surfaces',*99)
                        ENDIF


c...      Close the gap over closed C-Mod ports. Nasty hardcoded business:
          IF (.TRUE..AND.
     .        ((eirnsection.EQ.4 .AND.isec.NE.1).OR.
     .         (eirnsection.EQ.2 .AND.isec.NE.1).OR.
     .         (eirnsection.EQ.3 .AND.isec.EQ.2).OR.
     .         (eirnsection.EQ.10.AND.
     .          (isec.EQ.1.OR.isec.EQ.2.OR.isec.EQ.3.OR.isec.EQ.6.OR.
     .           isec.EQ.9)))) THEN
            IF ((in1.GE.72 .AND.in1.LE.103).OR.
     .          (in1.GE.140.AND.in1.LE.151)) THEN
c...          Pass on the gussets and port boundaries:
            ELSE
              IF (eirasdat(in1  ,8).EQ.0.0908) z(1)=z(1)+3.72
              IF (eirasdat(in1  ,9).EQ.0.0908) z(2)=z(2)+3.72
              IF (eirasdat(in1+1,8).EQ.0.0908) z(3)=z(3)+3.72
              IF (eirasdat(in1+1,9).EQ.0.0908) z(4)=z(4)+3.72

              IF (eirasdat(in1  ,8).EQ.0.1692) z(1)=z(1)-3.82
              IF (eirasdat(in1  ,9).EQ.0.1692) z(2)=z(2)-3.82
              IF (eirasdat(in1+1,8).EQ.0.1692) z(3)=z(3)-3.82
              IF (eirasdat(in1+1,9).EQ.0.1692) z(4)=z(4)-3.82
            ENDIF
          ENDIF


c...                    If a surface is infinite in toroidal extent, then do not duplicate
c                       it in every toroidal section (this is only if more than one toroidal
c                       section is defined by, i.e. EIRNSECTION > 1):
c                        WRITE(0,*) 'EIRNSECTION:',eirnsection,isec  
                        IF (isec.GT.1.AND.
     .                      (ABS(z(1)).GT.1.0E+10.OR.
     .                       ABS(z(2)).GT.1.0E+10.OR.
     .                       ABS(z(3)).GT.1.0E+10.OR.
     .                       ABS(z(4)).GT.1.0E+10)) THEN
                              WRITE(0,*) 'CYCLING PUFFING SURFACE:',isec,in2
                              CYCLE
                        ENDIF
                        count = -1
                        CALL DivAddSurface(x,y,z,count,0,99,0)
                        DO in4 = 1, count
                          nsrfsi = nsrfsi + 1
c...                      Store the location of the additional surface data, and the
c                         index of the toroidal subsegment so that the loop over
c                         NSRFSI below can function in a similar manner to the 
c                         simpler puffing geometries:
                          indseg(nsrfsi) = in1
                          cntseg(nsrfsi) = in4
                          sftseg(nsrfsi) = haddi2 * (isec - 1)
                        ENDDO
c                        WRITE(0,*) 'COUNT:',count,nsrfsi
                      ENDDO
                    ENDDO
                  ENDIF
                ELSE
                  nsrfsi = NINT(eirpuff(i1,4) - eirpuff(i1,3)) + 1
                ENDIF
              ELSE
                STOP 'problem'
              ENDIF

              totarea = 0.0
              npts = NINT(eirpuff(i1,10))
              nemods = 201
c              nemods = 217

c              WRITE(0,*) '*** MAKE SURE ADDITIONAL SURFACE SOURCES '//
c     .                   'ARE VALID ***'

              DO i2 = 1, nsrfsi

                iltor(i2) = 0
                orientation = 0

                in2 = NINT(eirpuff(i1,3)) + (i2 - 1)
	    
                IF (in2.LE.nvesm+nvesp) THEN
c...              "Standard" divimp wall surface:
	    
                  DO in1 = 1, walln
                    IF (wallr(in1,1).EQ.rvesm(in2,2).AND.
     .                  wallr(in1,2).EQ.rvesm(in2,1).AND.
     .                  wallz(in1,1).EQ.zvesm(in2,2).AND.
     .                  wallz(in1,2).EQ.zvesm(in2,1)) walli(i2) = in1 
                  ENDDO
                  IF (walli(i2).EQ.0) CALL ER('EIR07','EIRENE wall '//
     .                                     'index not found',*99)

c...              Assign coordinates of puffing region:
                  x1(i2) = wallr(walli(i2),1)
                  y1(i2) = wallz(walli(i2),1)
                  x2(i2) = wallr(walli(i2),2)
                  y2(i2) = wallz(walli(i2),2)
                  IF (eirntorseg.NE.0) THEN
                    zbnd1 = 0.0
                    zbnd2 = 2.0*PI * 0.5*(x1(i2) + x2(i2)) * eirtorfrac
                    z1(i2) = zbnd1
                    z2(i2) = zbnd2
                  ELSE
                    zbnd1 = defzbnd1
                    zbnd2 = defzbnd2
                    z1(i2) = zbnd1 / 100.0
                    z2(i2) = zbnd2 / 100.0
                  ENDIF

                  orientation = 1

                ELSE
	    
	    
c...              Make some assumption about orientation here.  The user
c                 will have to decide if it is correct:
                  IF (eirntorseg.NE.0) THEN
c...                Geometry data for the puffing surface is to be
c                   taken from a subsection of the additional surface (toroidal
c                   resolution outside the standard grid):

                    in1 = indseg(i2)

                    x(1) = eirasdat(in1,3)*100.0
                    y(1) = eirasdat(in1,4)*100.0
                    z(1) = eirasdat(in1,8)*100.0 + zshift*100.0
                    x(2) = eirasdat(in1,6)*100.0
                    y(2) = eirasdat(in1,7)*100.0
                    z(2) = eirasdat(in1,9)*100.0 + zshift*100.0

                    x(3) = eirasdat(in1+1,3)*100.0
                    y(3) = eirasdat(in1+1,4)*100.0
                    z(3) = eirasdat(in1+1,8)*100.0 + zshift*100.0
                    x(4) = eirasdat(in1+1,6)*100.0
                    y(4) = eirasdat(in1+1,7)*100.0
                    z(4) = eirasdat(in1+1,9)*100.0 + zshift*100.0

                    count = -1
                    CALL DivAddSurface(x,y,z,count,0,99,0)

                    DO in3 = 1, cntseg(i2)
                      CALL DivAddSurface(x,y,z,iltor(i2),0,99,0)
                    ENDDO

c...                Store vectors on surface for computing the normal:
                    a1 = x(1) - x(2)
                    a2 = y(1) - y(2)
                    a3 = z(1) - z(2)
                    b1 = x(3) - x(2)
                    b2 = y(3) - y(2)
                    b3 = z(3) - z(2)

                    x1(i2) = MIN(x(1),x(2),x(3),x(4))*0.01
                    x2(i2) = MAX(x(1),x(2),x(3),x(4))*0.01
                    y1(i2) = MIN(y(1),y(2),y(3),y(4))*0.01
                    y2(i2) = MAX(y(1),y(2),y(3),y(4))*0.01
                    z1(i2) = MIN(z(1),z(2),z(3),z(4))*0.01
                    z2(i2) = MAX(z(1),z(2),z(3),z(4))*0.01

c...                Index of surface in block 3b in the EIRENE input file:
                    walli(i2) = eirnlimi(in1+sftseg(i2))+cntseg(i2)-1

c                    WRITE(0,*) 'status:',in1,cntseg(i2),walli(i2)

c                    STOP 'DEV REQUIRED'

c... assign x and y
c... then find the cross-product
c... figure out walli

                  ELSE
c                    in2 = NINT(eirpuff(i1,3)) + (i2 - 1)
                    in3 = nvesm + nvesp
c...                Scan additional surfaces and find polygon data:
                    status = .FALSE.
                    DO in1 = 1, eirnasdat
                      IF (eirasdat(in1,1).EQ.7.0.OR.
     .                    eirasdat(in1,1).EQ.2.0) in3 = in3 + 1
                      IF (in3.EQ.in2) THEN
                        status = .TRUE.
                        EXIT
                      ENDIF
                    ENDDO
                    IF (.NOT.status) 
     .                CALL ER('Block07','Surface data not found',*99)	    

c....               Check that the surface is a rectangle - requried!:

c...                Store vectors on surface for computing the normal below:
                    a1 = eirasdat(in1  ,6) - eirasdat(in1,3)
                    a2 = eirasdat(in1  ,7) - eirasdat(in1,4)
                    a3 = eirasdat(in1  ,9) - eirasdat(in1,8)

                    b1 = eirasdat(in1+1,6) - eirasdat(in1,3)
                    b2 = eirasdat(in1+1,7) - eirasdat(in1,4)
                    b3 = eirasdat(in1+1,9) - eirasdat(in1,8)

                    x1(i2) = MIN(eirasdat(in1  ,3),eirasdat(in1  ,6),
     .                           eirasdat(in1+1,3),eirasdat(in1+1,6))
                    x2(i2) = MAX(eirasdat(in1  ,3),eirasdat(in1  ,6),
     .                           eirasdat(in1+1,3),eirasdat(in1+1,6))
                    y1(i2) = MIN(eirasdat(in1  ,4),eirasdat(in1  ,7),
     .                           eirasdat(in1+1,4),eirasdat(in1+1,7))
                    y2(i2) = MAX(eirasdat(in1  ,4),eirasdat(in1  ,7),
     .                           eirasdat(in1+1,4),eirasdat(in1+1,7))
                    z1(i2) = MIN(eirasdat(in1  ,8),eirasdat(in1  ,9),
     .                           eirasdat(in1+1,8),eirasdat(in1+1,9))
                    z2(i2) = MAX(eirasdat(in1  ,8),eirasdat(in1  ,9),
     .                           eirasdat(in1+1,8),eirasdat(in1+1,9))

c...                Index of surface in block 3b in the EIRENE input file:
                    walli(i2) = ebgki + NINT(eirpuff(i1,3)) + (i2 - 1) - 
     .                          (nvesm + nvesp) 
                  ENDIF

                  IF (eirpuff(i1,5).NE.0.0.AND.
     .                eirpuff(i1,6).NE.0.0) THEN
c...                Use a subsection of the surface:
                    temp1 = x1(i2) + eirpuff(i1,5) * (x2(i2) - x1(i2))
                    temp2 = x1(i2) + eirpuff(i1,6) * (x2(i2) - x1(i2))
                    x1(i2) = temp1
                    x2(i2) = temp2
                    temp1 = y1(i2) + eirpuff(i1,5) * (y2(i2) - y1(i2))
                    temp2 = y1(i2) + eirpuff(i1,6) * (y2(i2) - y1(i2))
                    y1(i2) = temp1
                    y2(i2) = temp2
                    temp1 = z1(i2) + eirpuff(i1,5) * (z2(i2) - z1(i2))
                    temp2 = z1(i2) + eirpuff(i1,6) * (z2(i2) - z1(i2))
                    z1(i2) = temp1
                    z2(i2) = temp2
                  ENDIF

	    
                ENDIF
	    
c...            Calculate the surface normal and make sure that it is
c               perpendicular to toroidal direction:
                IF (orientation.EQ.0) THEN
c...              Take the cross product to get the surface normal:
                  c1 = a2 * b3 - a3 * b2                
                  c2 = a3 * b1 - a1 * b3
                  c3 = a1 * b2 - a2 * b1                
c...              Examine the surface normal to decide if the surface is
c                 parallel or perpendicular to the toroidal axis:
                  IF     (ABS(c3).LT.1.0E-6) THEN
c...                Perpendicular:
                    orientation = 1
                  ELSEIF (ABS(c1).LT.1.0E-6.AND.ABS(c2).LT.1.0E-6) THEN
c...                Parallel:
                    orientation = 2
                  ELSE
                    CALL ER('Write07','Only surfaces parallel or '//
     .                      'perpendicular to the toroidal axis can '//
     .                      'be selected as puffing surfaces at '//
     .                      'this time',*99)
                  ENDIF
                ENDIF

                IF (.TRUE.) THEN
c...              Estimate the pitch of the magnetic field:
                  xcen = 0.5 * (x1(i2) + x2(i2))
                  ycen = 0.5 * (y1(i2) + y2(i2))
                  distmin = HI
                  DO ir = irsep, nrs
c...                Only want to check rings that are next to boundary rings:
                    IF (idring(ir).EQ.-1.OR.
     .                  (idring(    ir-1        ).NE.-1.AND.
     .                   idring(MIN(ir+1,MAXNRS)).NE.-1)) CYCLE
                    DO ik = 1, nks(ir)
c...                  Ignore virtual cells, if any:
                      IF (virtag(ik,ir).EQ.1) CYCLE
c...                  Calculate the distance to cell IK,IR and see if it is the 
c                     closest cell to the puffing surface in question (THIS IS
c                     A VERY CRUDE METHOD!  MAIN CHAMBER WALLS ARE POORLY HANDLED!):
                      dist = SQRT((xcen - rs(ik,ir))**2 +
     .                            (ycen - zs(ik,ir))**2)
                      IF (dist.LT.distmin) THEN
                        distmin = dist
                        brat = bratio(ik,ir)
                        id = korpg(ik,ir)
                        deltar = rvertp(1,id) - rvertp(4,id)
                        deltaz = zvertp(1,id) - zvertp(4,id)
                        iktemp = ik
                        irtemp = ir
                      ENDIF
                    ENDDO
                  ENDDO

c...              Estimate the poloidal inclination of the target:
                  cost = -1.0
                  DO WHILE (cost.LT.0.0)
                    alpha  = ATAN2C(deltar,deltaz)
                    deltar = x2(i2) - x1(i2)
                    deltaz = y2(i2) - y1(i2)
                    beta   = ATAN2C(deltar,deltaz) - alpha
                    cost   = COS(PI / 2.0 - beta)
                    IF (cost.LT.0.0) THEN
c                      WRITE(0,*) 'TRYING AGAIN'
                      id = korpg(iktemp,irtemp)
                      deltar = rvertp(4,id) - rvertp(1,id)
                      deltaz = zvertp(4,id) - zvertp(1,id)
                    ENDIF
                  ENDDO

c...              Calculate effective surface area:
                  IF     (orientation.EQ.1) THEN
                    sorwgt(i2) = SQRT((x2(i2) - x1(i2))**2 + 
     .                                (y2(i2) - y1(i2))**2) * 
     .                         (z2(i2) - z1(i2)) * brat * cost
                  ELSEIF (orientation.EQ.2) THEN
                    WRITE(0,*) 'WARNING: NOT SURE THIS WORKS PROPERLY'
                    brat = 1.0 - brat
                    sorwgt(i2) =(x2(i2) - x1(i2)) * 
     .                          (y2(i2) - y1(i2)) * brat * cost
                  ENDIF

                  IF (orientation.EQ.1.AND.eirntorseg.EQ.0) THEN
c...                Scale area according to major radius, but only if 
c                   the surface normal is perpendicular to the toroidal axis and
c                   the toroidal geometry approximation is not being used:

                    IF (eirzaa.LE.0.0) 
     .                CALL ER('Block07','Need to develop surface '//
     .                                  'area scaling',*99)

                    rad1 = eirzaa / (2.0 * PI * eirtorfrac)
                    rad2 = 0.5 * (x1(i2) + x2(i2))
                    sorwgt(i2) = sorwgt(i2) * (rad2 / rad1)
                  ENDIF

                  totarea = totarea + sorwgt(i2)

c...              Calcualte surface normal and normalize to unity:

c THIS WILL BE THE NORMAL FOR THE LAST SURFACE IN THE SET, SO BE CAREFUL
c IF BRAT IS EXTRAPOLATED IN THE FUTURE, OR THE POLOIDAL TRAJECTORY
c FOR THAT MATTER

                  sorctz = -COS(brat)
                  sorctx = -COS(PI/2.0 - alpha) * COS(PI/2.0 - brat)
                  sorcty = -COS(alpha)          * COS(PI/2.0 - brat)
                  len = SQRT(sorctx**2 + sorcty**2 + sorctz**2)
                  sorctx = sorctx / len 
                  sorcty = sorcty / len
                  sorctz = sorctz / len

c                  WRITE(0,*) 'WHERE :',iktemp,irtemp
c                  WRITE(0,'(1X,A,2I6,3F12.6,I6)') 
c     .              'ANGLES:',iktemp,irtemp,brat,cost,totarea,
c     .                        orientation
c                  WRITE(0,*) 'NORMAL:',sorctx,sorcty,sorctz


                ELSE
                  CALL ER('EIR07','Surface normal must be perpendi'//
     .                    'cular to the z-axis (in EIRENE)',*99)
                ENDIF


              ENDDO

              
c HOW IS THIS PASSED TO THE PARTICLE ACCOUNTING CODE IN TAU?
              IF (eirpuff(i1,2).EQ.-1.0) THEN
c...            Take flux from the outer target of the first PFZ ring:
                IF (lpdatsw.NE.1) 
     .            CALL ER('WriteBlock7','Target data must include '//
     .                                  'Isat',*99)
                flux = 0.0
                DO i3 = 1, nlpdati
                  IF (NINT(lpdati(i3,1)).EQ.irtrap+1) THEN                
                    flux = lpdati(i3,4) * totarea * eirsrcmul
                  ENDIF
                ENDDO
                IF (flux.EQ.0.0) 
     .            CALL ER('WriteBlock7','Isat data not found for '//
     .                                  'puff',*99)

              ELSEIF (eirpuff(i1,2).LT.-1.0) THEN
c...            Total flux specified in the input file:
                flux = -eirpuff(i1,2) * ECH * REAL(eirnsection) * 
     .                  eirsrcmul
              ELSE
                flux = eirpuff(i1,2) * totarea * eirsrcmul
              ENDIF

            ELSE
              CALL ER('WriteBlock07','Unrecognized flux option',*99)
            ENDIF

c...        This will only work for EIRPMODE=5 I think:
            eirpuff(i1,MAXPUFF-1) = eirpuff(i1,11)

            eirpuff(i1,MAXPUFF) = flux
            
            WRITE(PINOUT,*) 'PUFF FLUX = ',flux,eirpuff(i1,2)

            IF     (eirpmode.EQ.1) THEN

              WRITE(0,*) 'WARNING: USING OLD PUFFING MODE 1'

            ELSEIF (eirpmode.EQ.2) THEN
c...          Map xVESM index to EIRENE wall index:

              WRITE(0,*) 'WARNING: USING OLD PUFFING MODE 2'

            ELSEIF (eirpmode.EQ.3.OR.eirpmode.EQ.4) THEN

              sorifl = 1000.0

              WRITE(fp2,81) '* ADDITIONAL SURFACE SOURCE'
              WRITE(fp2,81) 'FFFFF'
              WRITE(fp2,84) npts,ninitl,nemods,1
              WRITE(fp2,80) flux
              IF     (eirpuff(i1,7).EQ.1.0) THEN
                WRITE(fp2,81) 'TFFFF'
              ELSEIF (eirpuff(i1,7).EQ.2.0) THEN
                WRITE(fp2,81) 'FTFFF'
              ELSEIF (eirpuff(i1,7).EQ.3.0) THEN
                WRITE(fp2,81) 'FFFTF'
              ELSE
                CALL ER('EIR07','Unknown puffing species',*99)
              ENDIF
              WRITE(fp2,81) '     1'
              WRITE(fp2,81) 'FFTFF'
              WRITE(fp2,84) nsrfsi 

              DO i2 = 1, nsrfsi
                IF (eirpuff(i1,1).EQ.5.0) THEN
                  IF (orientation.EQ.2) THEN                   
                    sorlim = -1.0
c                    sorlim = 0.2200E+02
                  ELSE
                    IF (x1(i2).EQ.x2(i2)) THEN
                      sorlim = 2.2000E+02
                    ELSE
                      sorlim = 2.0200E+02
                    ENDIF
                  ENDIF
                  sorifl = eirpuff(i1,9)
                  IF (eirpuff(i1,11).NE.0.0) THEN 
                    nasor = -NINT(eirpuff(i1,11)) - 1
                  ELSE
                    nasor = -1
                  ENDIF
                ELSE
c...              Map xVESM index to EIRENE wall index:
                  in2 = NINT(eirpuff(i1,3))
                  DO in1 = 1, walln
                    IF (wallr(in1,1).EQ.rvesm(in2,2).AND.
     .                  wallr(in1,2).EQ.rvesm(in2,1).AND.
     .                  wallz(in1,1).EQ.zvesm(in2,2).AND.
     .                  wallz(in1,2).EQ.zvesm(in2,1)) walli(i2) = in1 
                  ENDDO
                  IF (walli(i2).EQ.0) CALL ER('EIR07','EIRENE wall '//
     .                                        'index not found',*99)

c...              Find coordinates of puffing region:
                  x1(i2) = wallr(walli(i2),2) + eirpuff(i1,5) * 
     .                     (wallr(walli(i2),1) - wallr(walli(i2),2))
                  y1(i2) = wallz(walli(i2),2) + eirpuff(i1,5) * 
     .                     (wallz(walli(i2),1) - wallz(walli(i2),2))
                  x2(i2) = wallr(walli(i2),2) + eirpuff(i1,6) * 
     .                     (wallr(walli(i2),1) - wallr(walli(i2),2))
                  y2(i2) = wallz(walli(i2),2) + eirpuff(i1,6) * 
     .                     (wallz(walli(i2),1) - wallz(walli(i2),2))
                  IF (eirntorseg.NE.0) THEN
                    zbnd1 = 0.0
                    zbnd2 = 2.0*PI * 0.5*(x1(i2) + x2(i2)) * eirtorfrac
                    z1(i2) = zbnd1
                    z2(i2) = zbnd2
                  ELSE
                    zbnd1 = defzbnd1
                    zbnd2 = defzbnd2
                    z1(i2) = zbnd1 / 100.0
                    z2(i2) = zbnd2 / 100.0
                  ENDIF

                  sorwgt(i2) = 1.0

                  IF (wallr(walli(i2),1).EQ.wallr(walli(i2),2)) THEN
                    sorlim = 2.2000E+02
                  ELSE
                    sorlim = 2.0200E+02
                  ENDIF

                  IF     (eirpuff(i1,4).EQ.-2.0) THEN
c...                Use additional cell searching routine in EIRENE:
                    nasor = -1
                  ELSEIF (eirpuff(i1,4).EQ.-1.0) THEN
c...                Launch into EIRENE additional cell 1, not cell 1 on the
c                   vacuum grid:
                    nasor = 1
                  ELSEIF (NINT(eirpuff(i1,4)).GE.1        .AND.
     .                    NINT(eirpuff(i1,4)).LE.asc_ncell*ascncut) THEN
                    nasor = NINT(eirpuff(i1,4))+1+eirnpgdat
                  ELSE
                    CALL ER('WriteBlock07','Invalid injection cell',*99)
                  ENDIF

                  sorctx = 0.0
                  sorcty = 0.0
                  sorctz = 0.0

                ENDIF

                WRITE(fp2,84) 1,0,walli(i2),0,0,0,0
                WRITE(fp2,80) sorwgt(i2),sorlim,1.0,0.0,sorifl
c                WRITE(fp2,84) 0,0,iltor(i2),0,nasor
                WRITE(fp2,84) 0,0,0,0,nasor

                IF (eirntorseg.NE.0) THEN
c...              Convert the toroidal coordinate to angular coordinates for
c                 sampling in EIRENE (NLTRA=.TRUE, NLTOR=.FALSE. are valid):

c...              Toroidal circumference for this surface:
                  IF (iltor(i2).NE.0) THEN
                    torlen = 2.0 * PI * x2(i2)
                  ELSE
                    torlen = 2.0 * PI * 0.5 * (x1(i2) + x2(i2))
                  ENDIF

                  z1(i2) = z1(i2) / torlen * 360.0
                  z2(i2) = z2(i2) / torlen * 360.0

c                  WRITE(0,*) 'CHECKING:',torlen
c                  phiseg = 360.0 / REAL(eirntorseg) 
c                    WRITE(0,*) 'CHECKING:',i2,z1(i2),z2(i2),phiseg
c                  z1(i2) = z1(i2) + phiseg * REAL(iltor(i2)-1)
c                  z2(i2) = z2(i2) + phiseg * REAL(iltor(i2)-1)
c                    WRITE(0,*) 'CHECKING:',i2,z1(i2),z2(i2) 

                  z1(i2) = z1(i2) / 100.0
                  z2(i2) = z2(i2) / 100.0

                ENDIF

                WRITE(fp2,80) MIN(x1(i2)*100.0,x2(i2)*100.0),
     .                        MAX(x1(i2)*100.0,x2(i2)*100.0),
     .                        MIN(y1(i2)*100.0,y2(i2)*100.0),
     .                        MAX(y1(i2)*100.0,y2(i2)*100.0),
     .                        z1(i2)*100.0,z2(i2)*100.0
              ENDDO

c...          Velocity space distribution:
              IF (eirpuff(i1,8).EQ.-1.0) THEN
                id = idds(irtrap+1,1)
                ti = ktids(id)
                te = kteds(id)
c                WRITE(fp2,80) ktids(id),kteds(id),sorctx,
c     .                        sorcty,sorctz
              ELSE
                ti = eirpuff(i1,8)
                te = ti
c                WRITE(fp2,80) eirpuff(i1,8),eirpuff(i1,8),sorctx,
c     .                        sorcty,sorctz
              ENDIF

              IF (eirpuff(i1,7).EQ.3.0.AND.nemods.EQ.201) THEN
c...            Sheath energy gain from equation 38 in the EIRENE manual, but only
c               for "mixed" source (not strictly D or D2):
                esheat = 0.5 * LOG(te * crmb * 1835.0 / 
     .                   (2.0 * PI * (te + ti))) * cizb * te
                soreni = 3.0 * ti + 0.5 * te + esheat
                sorene = 0.0
              ELSE
                soreni = ti
                sorene = te
              ENDIF
              WRITE(fp2,80) soreni,sorene,sorctx,
     .                      sorcty,sorctz
              WRITE(fp2,81) '  1.0000E 00  9.0000E 01'

            ENDIF

          ENDDO

80        FORMAT(1P,6(E12.4:),0P)         
81        FORMAT(A)
82        FORMAT(6(E12.4:))
83        FORMAT(2E12.4)
84        FORMAT(7(I6:))
        ENDIF

      ELSEIF (eir_07opt.EQ.1.OR.eir_07opt.EQ.2.OR.eir_07opt.EQ.5) THEN

        WRITE(fp2,'(I6         )') eir_07nsrc
        WRITE(fp2,'(20(I6:)    )') (1,i2=1,eir_07nsrc)
        WRITE(fp2,'(1P,E12.4,0P)') eir_07wght

        DO i1 = 1, eir_07nsrc
          IF     (eir_07stra(i1).EQ.1) THEN
            WRITE(fp2,'(A,I2        )') '* INNER SURFACE (DIVIMP)',i1
            WRITE(fp2,'(A           )') 'FFFFF'
            WRITE(fp2,'(4I6         )') INT(eir_07data(i1,1)),ninitl,3,1
            WRITE(fp2,'(1P,E12.4,0P )') 1.0
            WRITE(fp2,'(A           )') 'FFFTF'
            WRITE(fp2,'(I6          )') -1
            WRITE(fp2,'(A           )') 'FFTFF'
            WRITE(fp2,'(I6          )')  1
            WRITE(fp2,'(5I6         )') 1,4,eir_07ind1(i1),
     .                                      eir_07ind2(i1),
     .                                      eir_07ind3(i1)
            WRITE(fp2,'(1P,5E12.4,0P)') 1.0,2.04E+02,1.0,0.0,2.0E+03
            IF (eirnsdtor.GT.1) THEN
              WRITE(fp2,84) -1,0,0,999
            ELSE
              WRITE(fp2,84) -1,0,0,0
            ENDIF
            WRITE(fp2,'(1P,6E12.4,0P)') 0.0,0.0,0.0,0.0,zbnd1,zbnd2
            WRITE(fp2,'(1P,2E12.4,0P)') 3.0,0.5
            WRITE(fp2,'(1P,2E12.4,0P)') 1.0,0.0
c            WRITE(fp2,'(1P,2E12.4,0P)') 1.0,90.0
          ELSEIF (eir_07stra(i1).EQ.2) THEN
            WRITE(fp2,'(A,I2        )') '* OUTER SURFACE (DIVIMP)',i1
            WRITE(fp2,'(A           )') 'FFFFF'
            WRITE(fp2,'(4I6         )') INT(eir_07data(i1,1)),ninitl,3,1
            WRITE(fp2,'(1P,E12.4,0P )') 1.0
            WRITE(fp2,'(A           )') 'FFFTF'
            WRITE(fp2,'(I6          )') -1
            WRITE(fp2,'(A           )') 'FFTFF'
            WRITE(fp2,'(I6          )')  1
            WRITE(fp2,'(5I6         )') 1,4,eir_07ind1(i1),
     .                                      eir_07ind2(i1),
     .                                      eir_07ind3(i1)
            WRITE(fp2,'(1P,5E12.4,0P)') 1.0,2.04E+02,2.0,0.0,1.0E+03
            IF (eirnsdtor.GT.1) THEN
              WRITE(fp2,84) -1,0,0,999
            ELSE
              WRITE(fp2,84) -1,0,0,0
            ENDIF
            WRITE(fp2,'(1P,6E12.4,0P)') 0.0,0.0,0.0,0.0,zbnd1,zbnd2
            WRITE(fp2,'(1P,2E12.4,0P)') 3.0,0.5
            WRITE(fp2,'(1P,2E12.4,0P)') 1.0,0.0
c            WRITE(fp2,'(1P,2E12.4,0P)') 1.0,90.0
          ELSEIF (eir_07stra(i1).EQ.3.OR.eir_07stra(i1).EQ.4) THEN
            WRITE(fp2,'(A,I2        )') '* VOLUME SOURCE (DIVIMP)'
            WRITE(fp2,'(A           )') 'FFFFF'
            WRITE(fp2,'(4I6         )') INT(eir_07data(i1,1)),ninitl,3,1
            WRITE(fp2,'(1P,E12.4,0P )') 1.0
            WRITE(fp2,'(A           )') 'FFFTF'
            WRITE(fp2,'(I6          )') 1
            WRITE(fp2,'(A           )') 'FFFTF'
            WRITE(fp2,'(I6          )')  1
            WRITE(fp2,'(7I6         )') 1,4,39,1,irwall-1,
     .                                  eir_07ind1(i1),eir_07ind2(i1)
c            WRITE(fp2,'(5I6         )') 1,4,0,0,0
            WRITE(fp2,'(1P,5E12.4,0P)') 1.0,2.04E+02,1.0,0.0,1.0E+03
            IF (eirnsdtor.GT.1) THEN
              WRITE(fp2,84) -1,0,0,999
            ELSE
              WRITE(fp2,84) -1,0,0,0
            ENDIF
c            WRITE(fp2,'(4I6         )') -1,0,0,0
            WRITE(fp2,'(1P,6E12.4,0P)') 0.0,0.0,0.0,0.0,zbnd1,zbnd2
            WRITE(fp2,'(1P,2E12.4,0P)') 3.0,0.5
            WRITE(fp2,'(1P,2E12.4,0P)') 0.0,90.0
          ELSE
            CALL ER('WriteBlock07','Unknown source specification',*99)
          ENDIF

        ENDDO

      ELSE
        CALL ER('WriteBlock07','Invalid option',*99)
      ENDIF
c
c     Advance template file to the beginning of the next input block:
c
      buffer(1:3) = '   '
      DO WHILE (buffer(1:3).NE.'***')
        CALL ReadLine(fp1,buffer,1,*97,*98)
      ENDDO
      BACKSPACE fp1


c      CLOSE (fp2)
c      STOP 'WRITE INPUT BLOCK 7'

      RETURN
97    CALL ER('WriteBlock07','Unexpected end of file'        ,*99)
98    CALL ER('WriteBlock07','Problems reading template file',*99)
99    STOP
      END
c
c ======================================================================
c
c Write EIRENE geometry file:
c
c ======================================================================
c
c subroutine: WriteGeometryFile
c
c
c
c
c
      SUBROUTINE WriteGeometryFile
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER nnks,dimxh,dimyh,nncut,nxcut1,nxcut2

      INTEGER           fp,ik,ik1,ir1,ir,id,ikreal,ir2(MAXNRS),i1,si,
     .                  eirik,eirir
      DOUBLE PRECISION  x1,x2,x3,x4,y1,y2,y3,y4
c
c     Check to see if DIVIMP should write the geometry
c     file for EIRENE:
c
      IF (eirgeom.EQ.0) RETURN

c...  Some tricky business to make sure that the surface normals
c     are calculated for ring IRBREAK in EIRENE:
      IF (.FALSE..AND.nbr.GT.0) THEN
        ir     = irbreak
        ikreal = 0
        IF (stopopt.EQ.121) THEN
          DO ik = nks(ir)/2, nks(ir)
            IF (virtag(ik,ir).EQ.0) THEN 
              ikreal = ik
            ELSE
              IF (ikreal.EQ.0) 
     .          CALL ER('WriteGeometryFile','Unable to approximate '//
     .                  'broken grid surfaces',*99)

              WRITE(0,*) 'FIXING NORMALS:',ik,ir,virtag(ik,ir),ikreal

              rvertp(1,korpg(ik,ir)) = rvertp(1,korpg(ikreal,ir))
              zvertp(1,korpg(ik,ir)) = zvertp(1,korpg(ikreal,ir))
              rvertp(4,korpg(ik,ir)) = rvertp(4,korpg(ikreal,ir))
              zvertp(4,korpg(ik,ir)) = zvertp(4,korpg(ikreal,ir))
            ENDIF
          ENDDO
        ELSE
          WRITE(0,*) 'CODE DEVELOPMENT REQUIRED: NEED TO MAKE THIS '//
     .               'WORK FOR GRIDS WITH BROKEN RINGS OUTSIDE THE '//
     .               'PFZ'
        ENDIF
      ENDIF


c      fp = 0
      fp = 52
c
c     Initialisation:
c
      nnks   = nks(irsep)
      dimxh  = nnks
      dimyh  = irwall-2
      nncut  = 2
      IF (grdnmod.NE.0) THEN
        nxcut1 = cutpt1 - 1
        nxcut2 = cutpt2 - 1
      ELSE
        nxcut1 = ikto
        nxcut2 = ikti
      ENDIF

      OPEN(UNIT=52,ACCESS='SEQUENTIAL',STATUS='REPLACE')
c
c     Write the header to the geometry file:
c
      WRITE(fp,'(A)') 'DIVIMP/EIRENE Geometry Data'
      WRITE(fp,*)
      WRITE(fp,'(A30,I6)') 'dimhx                         ',dimxh
      WRITE(fp,'(A30,I6)') 'dimyh                         ',dimyh
      WRITE(fp,'(A30,I6)') 'nncut                         ',nncut
      WRITE(fp,'(A30,I6)') 'nxcut1+1                      ',nxcut1+1
      WRITE(fp,'(A30,I6)') 'nxcut2                        ',nxcut2
      WRITE(fp,'(A30,I6)') '[dummy]                       ',0
      WRITE(fp,'(A30,I6)') '[dummy]                       ',0
      WRITE(fp,'(A30,I6)') '[dummy]                       ',0
      WRITE(fp,*)

      eirik = 0
      DO ik = 1, nnks
        eirik = eirik + 1
        eirir = 0
        DO ir = 2, irwall
          eirir = eirir + 1
c
c         Determine the cell index:
c
          IF (ir.LT.irsep) THEN
            IF (.FALSE..AND.grdnmod.NE.0) THEN
              IF (ik.LE.ikto2(irtrap+ir-1)) THEN
                ik1 = ik
                ir1 = irtrap + ir - 1
              ELSEIF (ik.GE.ikti2(irtrap+ir-1)+nks(ir)-1) THEN
                ik1 = ik - nks(ir) + 1 
                ir1 = irtrap + ir - 1
              ELSE
                ik1 = ik - ikto2(irtrap+ir-1)
                ir1 = ir
              ENDIF
              WRITE(SLOUT,*) 'WRITEGEOM:',ir,ik1,ir1
            ELSE
              IF (ik.LE.nxcut1) THEN
                ik1 = ik
                ir1 = irtrap+ir-1
              ELSEIF (ik.GE.nxcut2) THEN
                ik1 = ik-(nxcut2-nxcut1-1)
                ir1 = irtrap+ir-1
              ELSE
                ik1 = ik-nxcut1
                ir1 = ir
              ENDIF
            ENDIF
          ELSEIF (ir.EQ.irwall) THEN
            ik1 = ik
            ir1 = ir-1
          ELSE
            ik1 = ik
            ir1 = ir
          ENDIF

          id = korpg(ik1,ir1)
c
c         Write geometry data:
c
          x1 = rvertp(1,id)
          y1 = zvertp(1,id)
          x2 = rvertp(2,id)
          y2 = zvertp(2,id)

          IF (ir.EQ.irwall) THEN
            WRITE(fp,3334)
     .        x2,y2,0.0,0.0,ir1,ik1,' n',eirir,eirik
          ELSE
            WRITE(fp,3334)
     +        x1,y1,x2,y2,ir1,ik1,' n',eirir,eirik
          ENDIF
        ENDDO
c
c       Add in the extra cell information around cut points:
c
        IF (ik.EQ.nxcut1.OR.ik.EQ.nxcut2-1.OR.ik.EQ.nnks) THEN
          eirik = eirik + 1
          eirir = 0
          DO ir = 2, irwall
            eirir = eirir + 1
c
c           Determine the cell index:
c
            IF (ir.LT.irsep) THEN
              IF (.FALSE..AND.grdnmod.NE.0) THEN
                IF (ik.LE.ikto2(irtrap+ir-1)) THEN
                  ik1 = ik 
                  ir1 = irtrap + ir - 1
                ELSEIF (ik.LE.ikto2(irtrap+ir-1)+nks(ir)-1) THEN
                  ik1 = ik - ikto2(irtrap+ir-1)
                  ir1 = ir
c...   Check the .slm file and be sure that IK1,IR1 are correct.
                  STOP 'STOP: UNTESTED CODE'
                ELSE
                  ik1 = ik - nks(ir) + 1
                  ir1 = irtrap + ir - 1
                ENDIF
              ELSE
                IF (ik.EQ.nxcut1) THEN
                  ik1 = ik
                  ir1 = irtrap+ir-1
                ELSEIF (ik.EQ.nxcut2-1) THEN
                  ik1 = nks(ir)-1
                  ir1 = ir
                ELSEIF (ik.EQ.nnks) THEN
                  ik1 = nks(irtrap+ir-1)
                  ir1 = irtrap+ir-1
                ENDIF
              ENDIF
                WRITE(SLOUT,*) 'WRITEPAD :',ir,ik1,ir1

            ELSEIF (ir.EQ.irwall) THEN
              ik1 = ik
              ir1 = ir - 1
            ELSE
              ik1 = ik
              ir1 = ir
            ENDIF

            id = korpg(ik1,ir1)

            x3 = rvertp(3,id)
            y3 = zvertp(3,id)
            x4 = rvertp(4,id)
            y4 = zvertp(4,id)

            IF ((ir.EQ.irwall.AND.ik.EQ.nxcut1  ).OR.
     .          (ir.EQ.irwall.AND.ik.EQ.nxcut2-1).OR.
     .          (ir.EQ.irwall.AND.ik.EQ.nnks    )) THEN
              WRITE(fp,3334)
     .          x3,y3,0.0,0.0,ir1,ik1,' p',eirir,eirik
            ELSE
              WRITE(fp,3334)
     .          x4,y4,x3,y3,ir1,ik1,' p',eirir,eirik
            ENDIF
          ENDDO
        ENDIF

      ENDDO

c...  Also transfer connection map, with the appropriate adjustments
c     for the differences between OEDGE and EIRENE ring indexing:

      WRITE(52,*)
      WRITE(52,*) 'Connection map:'

      WRITE(52,*)
      WRITE(52,'(A4,2X,100(I4:))') 'IK',(i1-1,i1=2,irwall-1)
      WRITE(52,*)

      si = 0
      DO ik = 1, nks(irsep)
        ir2 = -999
        DO ir = 2, irwall-1 
          IF (ir.LT.irsep) THEN
            IF     (ik.LE.cutpt1-1) THEN
c...          Low IK index PFZ cell:     
              ir2(ir) = irins(ik,ir+irtrap-1) - irtrap 
              WRITE(SLOUT,*) 'CONMAP A:',ik,ir+irtrap-1,ir2(ir)
            ELSEIF (ik.GE.cutpt2-1) THEN
c...          High IK index PFZ cell:
              ir2(ir) = irins(ik-nks(2)+1,ir+irtrap-1) - irtrap
              WRITE(SLOUT,*) 'CONMAP B:',ik-nks(2)+1,ir+irtrap-1,ir2(ir)
           ELSE
c...          Core cell:
              ir2(ir) = irins(ik-cutpt1+1,ir) - 1
              WRITE(SLOUT,*) 'CONMAP C:',ik-cutpt1+1,ir,ir2(ir)
            ENDIF
          ELSE
            IF (ir.EQ.irsep.AND.(ik.LE.ikto.OR.ik.GE.ikti)) THEN
c            IF (ir.EQ.irsep.AND.(ik.LE.cutpt1-1.OR.ik.GE.cutpt2-1)) THEN
              ir2(ir) = irins(ik,ir) - irtrap
              WRITE(SLOUT,*) 'CONMAP D:',ik,ir,ir2(ir)
            ELSE
              ir2(ir) = irins(ik,ir) - 1
              WRITE(SLOUT,*) 'CONMAP E:',ik,ir,ir2(ir)
            ENDIF
          ENDIF
        ENDDO
 
        IF (ir2(ir).EQ.irwall-1) ir2(ir) = MIN(ir,ir2(ir))

        WRITE(52,3335) ik+si,(MAX(-1,ir2(i1)),i1=2,irwall-1)
c...
        IF (ik.EQ.cutpt1-1.OR.ik.EQ.cutpt2-2.OR.ik.EQ.nks(irsep)) THEN
          si = si + 1
              WRITE(SLOUT,*) 'CONMAP Z:',ik+si
          WRITE(52,3335) ik+si,(-1,i1=2,irwall-1)
        ENDIF
      ENDDO


      WRITE(52,*)
      WRITE(52,'(A4,2X,100(I4:))') 'IK',(i1-1,i1=2,irwall-1)
      WRITE(52,*)

      si = 0
      DO ik = 1, nks(irsep)
        ir2 = -999
        DO ir = 2, irwall-1 

          IF (ir.LT.irsep) THEN

            IF     (ik.LE.cutpt1-1) THEN
c...          Low IK index PFZ cell:     
              IF (ir.EQ.irsep-1) THEN
c...            This will fail if the target is discontinuous across the PFZ/SOL boundary:
                ir2(ir) = irouts(ik,ir+irtrap-1) - 1 
              ELSE
                ir2(ir) = irouts(ik,ir+irtrap-1) - irtrap 
              ENDIF
              WRITE(SLOUT,*) 'CONMAP A:',ik,ir+irtrap-1,ir2(ir)
            ELSEIF (ik.GE.cutpt2-1) THEN
c...          High IK index PFZ cell:
              IF (ir.EQ.irsep-1) THEN
c...            This will fail if the target is discontinuous across 
c               the PFZ/SOL boundary:
                ir2(ir) = irouts(ik-nks(2)+1,ir+irtrap-1) - 1
                WRITE(SLOUT,*) 'CONMAP B:',ik-nks(2)+1,ir+irtrap-1,
     .                         ir2(ir)
              ELSE
                ir2(ir) = irouts(ik-nks(2)+1,ir+irtrap-1) - irtrap
                WRITE(SLOUT,*) 'CONMAP C:',ik-nks(2)+1,ir+irtrap-1,
     .                         ir2(ir)
              ENDIF
            ELSE
c...          Core cell:
              ir2(ir) = irouts(ik-cutpt1+1,ir) - 1
              WRITE(SLOUT,*) 'CONMAP D:',ik-cutpt1+1,ir,ir2(ir)
            ENDIF

          ELSE

            ir2(ir) = irouts(ik,ir) - 1
            WRITE(SLOUT,*) 'CONMAP E:',ik,ir,ir2(ir)

          ENDIF

          IF (ir2(ir).EQ.irwall-1) ir2(ir) = MIN(ir,ir2(ir))
         
        ENDDO
        WRITE(52,3335) ik+si,(MAX(-1,ir2(i1)),i1=2,irwall-1)
c...
        IF (ik.EQ.cutpt1-1.OR.ik.EQ.cutpt2-2.OR.ik.EQ.nks(irsep)) THEN
          si = si + 1
          WRITE(52,3335) ik+si,(-1,i1=2,irwall-1)
          WRITE(SLOUT,*) 'CONMAP Z:',si
        ENDIF
      ENDDO

c...  With the advent of fully double-null generalized grids, I need to 
c     tell EIRENE which radial side of each ring is associated with 3A
c     sides:
      WRITE(52,*)
      WRITE(52,*) 'Non-default radial surface map:'
      WRITE(52,*)
      DO ir = 2, irwall-1
        IF (ir.EQ.2) THEN
          ir2(ir) = 0
        ELSE
          ir2(ir) = 1
c...      Identify secondary PFZ boundary ring(s):
c          DO ik = 1, nks(ir)
c            IF (irins(ik,ir).EQ.irwall) ir2(ir) = 0
c          ENDDO
c...*HACK*
          IF (ir.EQ.21) ir2(ir) = 0

        ENDIF
        WRITE(52,*) ir-1,ir2(ir)
      ENDDO
c...  One more for outer radial boundary of EIRENE grid:
      WRITE(52,*) ir-1,1


      CLOSE(52)
c      STOP 'sdfsd'


c      CALL OutputData(85,'Dumping geometry')
c      STOP 'sdfsdfd'

      RETURN
3334  FORMAT(4E15.7,5X,2I6,A2,2I6)
3335  FORMAT(I4,2X,100(I4:))
99    STOP
      END
c
c ======================================================================
c
c
c ======================================================================
c
c subroutine: EstimateOpacityMultipier
c
c
c
c
      SUBROUTINE EstimateOpacityMultiplier
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      REAL InterpolateValue,CalcWidth

      INTEGER    MAXIN1  ,MAXIN2
      PARAMETER (MAXIN1=6,MAXIN2=5)
      

      REAL rectmp(MAXIN1),recden(MAXIN2),
     .     rectra(MAXIN1,MAXIN2),iontra(MAXIN1,MAXIN2),
     .     recopa(MAXIN1,MAXIN2),ionopa(MAXIN1,MAXIN2)

c...OLD RATE DATA : OBSOLETE
      DATA rectmp /-0.69  , -0.22  , 0.176  , 0.653  , 1.000,  2.000/
c      DATA rectmp /0.2    , 0.6    , 1.5    , 4.5    , 10.0   /
      DATA recden /1.0E+16, 1.0E+18, 1.0E+20, 1.0E+21, 1.0E+22/

      DATA rectra /1.8E-12, 6.0E-13, 3.0E-13, 1.4E-13, 6.0E-14, 5.0E-15,
     .             5.0E-12, 1.0E-12, 4.0E-13, 1.5E-13, 7.0E-14, 5.0E-15,
     .             6.5E-11, 4.0E-12, 7.0E-13, 1.8E-13, 8.0E-14, 5.0E-15,
     .             5.0E-10, 1.6E-11, 1.4E-12, 2.5E-13, 9.0E-14, 5.0E-15,
     .             1.0E-08, 1.0E-10, 2.5E-12, 3.0E-13, 1.0E-13, 5.0E-15/
      DATA recopa /1.5E-12, 5.5E-13, 1.2E-13, 6.5E-14, 4.0E-14, 4.5E-15,
     .             4.0E-12, 8.0E-13, 1.3E-13, 6.5E-14, 4.0E-14, 4.5E-15,
     .             6.0E-11, 1.4E-12, 1.4E-13, 6.5E-14, 4.0E-14, 4.5E-15,
     .             5.0E-10, 5.0E-12, 1.5E-13, 6.5E-14, 4.0E-14, 4.5E-15,
     .             1.0E-08, 4.5E-11, 5.5E-13, 1.0E-13, 4.5E-14, 4.5E-15/

      DATA iontra /1.0E-15, 1.0E-15, 1.5E-12, 1.5E-12, 1.5E-12, 1.5E-12,
     .             1.0E-15, 1.0E-15, 2.0E-12, 2.0E-12, 2.0E-12, 2.0E-12,
     .             1.0E-15, 1.0E-15, 6.0E-12, 6.0E-12, 6.0E-12, 6.0E-12,
     .             1.0E-15, 1.0E-15, 2.0E-11, 2.0E-11, 2.0E-11, 2.0E-11,
     .             1.0E-15, 1.0E-15, 4.5E-11, 4.5E-11, 4.5E-11, 4.5E-11/
      DATA ionopa /2.5E-15, 2.5E-15, 5.0E-11, 5.0E-11, 5.0E-11, 5.0E-11,
     .             2.5E-15, 2.5E-15, 5.0E-11, 5.0E-11, 5.0E-11, 5.0E-11,
     .             2.5E-15, 2.5E-15, 5.0E-11, 5.0E-11, 5.0E-11, 5.0E-11,
     .             2.5E-15, 2.5E-15, 5.0E-11, 5.0E-11, 5.0E-11, 5.0E-11,
     .             2.5E-15, 2.5E-15, 5.0E-11, 5.0E-11, 5.0E-11, 5.0E-11/


      DOUBLE PRECISION RPS10        ,RPS30
      PARAMETER       (RPS10=1.0E-10,RPS30=1.0E-30)

      REAL GetEAD

      INTEGER ik,ir,id
      LOGICAL checkmfp,highrec(MAXNKS,MAXNRS)
      REAL    v1,v2,v3,v4,lambdaHLa,celldim,ratio,maxrec

c      v1 = LOG10(10.0)
c      v2 = 1.0E+15

      CALL RSet(mulrec,MAXNKS*MAXNRS,1.0)
      CALL RSet(mulion,MAXNKS*MAXNRS,1.0)
      CALL RSet(mulqer,MAXNKS*MAXNRS,1.0)
      CALL RSet(mulqei,MAXNKS*MAXNRS,1.0)

c      WRITE(0,*) GetEAD(1.0,1.0E+21,3,'H.4 '),
c     .           GetEAD(1.0,1.0E+21,4,'H.4 ')

      tagmulrec = .TRUE.

      checkmfp = .TRUE.


      IF (tagpinatom.AND.eiropacity.EQ.4) THEN
        WRITE(0,*) 'NEED TO CLARIFY EIROPACITY=4 HERE AND IN'//
     .             'CALCINITRECOM'
c        STOP 'HALTING PROGRAM'
      ENDIF


c      IF ((eiropacity.EQ.3.OR.eiropacity.EQ.4).AND..NOT.tagpinatom) THEN
      IF ((eiropacity.EQ.3.AND..NOT.tagpinatom).OR.
     .    eiropacity.EQ.4.OR.eiropacity.EQ.5.OR.eiropacity.EQ.6) THEN
        WRITE(0     ,*) '*** NOT CHECKING LYMAN PHOTON MFP ***'
        WRITE(PINOUT,*) '*** NOT CHECKING LYMAN PHOTON MFP ***'
        checkmfp = .FALSE.
      ENDIF

c...CHECK THAT EIRENE ATOMIC DATA IS LOADED

c...  Find maxium recombination:
      maxrec = -HI
      DO ir = 1, nrs
        IF (idring(ir).EQ.-1) CYCLE
        DO ik = 1, nks(ir)
          maxrec = MAX(maxrec,pinrec(ik,ir))
        ENDDO
      ENDDO
c...  Tag cells that have recombination greater than 1% of the maximum:
      DO ir = 1, nrs
        IF (idring(ir).EQ.-1) CYCLE
        highrec(ik,ir) = .FALSE.
        DO ik = 1, nks(ir)
          IF (pinrec(ik,ir).GT.0.01*maxrec) THEN
            highrec(ik,ir) = .TRUE.
c            WRITE(PINOUT,*) 'HIGH RECOMBINATION:',ik,ir
          ENDIF
        ENDDO
      ENDDO  


      DO ir = irsep, nrs
        IF (idring(ir).EQ.-1) CYCLE

        DO ik = 1, nks(ir)
c          v1= LOG10(ktebs(ik,ir))
c          v2= knbs(ik,ir)
c          v3= InterpolateValue(v1,v2,rectmp,recden,rectra,MAXIN1,MAXIN2)
c          v4= InterpolateValue(v1,v2,rectmp,recden,recopa,MAXIN1,MAXIN2)
c          mulrec(ik,ir) = v4 / (v3 + RPS30)


          IF     (eiropacity.EQ.6) THEN
          ELSEIF (eiropacity.EQ.5) THEN
            v3=GetEAD(ktebs(ik,ir),MAX(1.0E+13,knbs(ik,ir)), 3,'H.4 ')
            v4=GetEAD(ktebs(ik,ir),MAX(1.0E+13,knbs(ik,ir)),27,'H.4 ')
            v1= v4 / (v3 + RPS30)
            mulrec(ik,ir) = v1 
            mulqer(ik,ir) = v1
          ELSE
            v3=GetEAD(ktebs(ik,ir),MAX(1.0E+13,knbs(ik,ir)),3,'H.4 ')
            v4=GetEAD(ktebs(ik,ir),MAX(1.0E+13,knbs(ik,ir)),4,'H.4 ')
            v1= v4 / (v3 + RPS30)
            mulrec(ik,ir) = v1 * mulrecl(ir)
            mulqer(ik,ir) = v1
          ENDIF


c

          IF     (eiropacity.EQ.6) THEN

            v3=GetEAD(ktebs(ik,ir),knbs(ik,ir),1 ,'H.4 ')
            v4=GetEAD(ktebs(ik,ir),knbs(ik,ir),28,'H.4 ')
            v1= v4 / (v3 + RPS30)

            mulion(ik,ir) = v1

          ELSEIF (eiropacity.EQ.2.OR.eiropacity.EQ.6.OR.
     .        (eiropacity.EQ.4.AND.highrec(ik,ir))) THEN

            v3=GetEAD(ktebs(ik,ir),knbs(ik,ir),1,'H.4 ')
            v4=GetEAD(ktebs(ik,ir),knbs(ik,ir),2,'H.4 ')
            v1= v4 / (v3 + RPS30)

            mulion(ik,ir) = v1

            v3=GetEAD(ktebs(ik,ir),knbs(ik,ir),5,'H.10')
            v4=GetEAD(ktebs(ik,ir),knbs(ik,ir),6,'H.10')
            v1= v4 / (v3 + RPS30)

            mulqei(ik,ir) = v1

          ENDIF




c...bug: ratio was not being calculated and the multipliers were not
c        begin adjusted, so they were always on everywhere - Nov 4, 1999
          IF (checkmfp.OR.stopopt2.EQ.4.OR.stopopt2.EQ.5) THEN
c          IF (stopopt2.EQ.4.OR.stopopt2.EQ.5) THEN

            lambdaHLa = 1.6E-03 * SQRT(ktebs(ik,ir) / 1.0) *
     .                 (1.0E+20 / (pinatom(ik,ir) + RPS10))

            id = korpg(ik,ir)
c            celldim = CalcWidth(ik,ir,CENTER,TOTAL)
            celldim = SQRT((rvertp(1,id) - rvertp(3,id))**2.0 +
     .                     (zvertp(1,id) - zvertp(3,id))**2.0)

            ratio = MAX(1.0,MIN(100.0,lambdaHLa / celldim))

            mulrec(ik,ir) = mulrec(ik,ir) +
     .                     (ratio - 1.0) / 99.0 * (1.0 - mulrec(ik,ir))
            mulion(ik,ir) = mulion(ik,ir) +
     .                     (ratio - 1.0) / 99.0 * (1.0 - mulion(ik,ir))
            mulqer(ik,ir) = mulqer(ik,ir) +
     .                     (ratio - 1.0) / 99.0 * (1.0 - mulqer(ik,ir))
            mulqei(ik,ir) = mulqei(ik,ir) +
     .                     (ratio - 1.0) / 99.0 * (1.0 - mulqei(ik,ir))

c            ratio = MAX(10.0,MIN(100.0,ratio))
c            mulrec(ik,ir) = mulrec(ik,ir) +
c     .                     (ratio - 10.0) / 90.0 * (1.0 - mulrec(ik,ir))
c            mulion(ik,ir) = mulion(ik,ir) +
c     .                     (ratio - 10.0) / 90.0 * (1.0 - mulion(ik,ir))
c            WRITE(6,'(A,2I4,F6.2,1P,E10.2,0P,2X,2F12.6,2X,
c     .                F10.4,2F14.7)')
c     .        'TEST : ',ik,ir,10.0**v1,v2,lambdaHLa,celldim,
c     .        ratio,mulrec(ik,ir),mulion(ik,ir)
          ENDIF



c          v3=GetEAD(ktebs(ik,ir),knbs(ik,ir),7,'H.10')
c          v4=GetEAD(ktebs(ik,ir),knbs(ik,ir),8,'H.10')
c          v1= v4 / (v3 + RPS30)

c          WRITE(6,'(A,2I4,1P,2E12.4,0P,3F10.4,2X,F8.2,1P,E12.4,0P)')
c     .      'COOL= ',ik,ir,v3,v4,v1,ratio,mulqer(ik,ir),
c     .      ktebs(ik,ir),knbs(ik,ir)

c          v3=GetEAD(ktebs(ik,ir),knbs(ik,ir),5,'H.10')
c          v4=GetEAD(ktebs(ik,ir),knbs(ik,ir),6,'H.10')
c          v1= v4 / (v3 + RPS30)

c          WRITE(6,'(A,2I4,1P,2E12.4,0P,3F10.4,2X,F8.2,1P,E12.4,0P)')
c     .      'COOL= ',ik,ir,v3,v4,v1,ratio,mulqei(ik,ir),
c     .      ktebs(ik,ir),knbs(ik,ir)

c          v3=GetEAD(ktebs(ik,ir),knbs(ik,ir),3,'H.4 ')
c          v4=GetEAD(ktebs(ik,ir),knbs(ik,ir),4,'H.4 ')
c          v1= v4 / (v3 + RPS30)

          WRITE(6,'(A,2I4,1P,2E12.4,0P,3F10.4,2X,F8.2,1P,E12.4,0P)')
     .      'RECMUL= ',ik,ir,v3,v4,v1,ratio,mulrec(ik,ir),
     .      ktebs(ik,ir),knbs(ik,ir)


c          v3=GetEAD(ktebs(ik,ir),knbs(ik,ir),1,'H.4 ')
c          v4=GetEAD(ktebs(ik,ir),knbs(ik,ir),2,'H.4 ')
c          v1= v4 / (v3 + RPS30)

c          WRITE(6,'(A,2I4,1P,2E12.4,0P,3F10.4,2X,F8.2,1P,E12.4,0P)')
c     .      'IONMUL= ',ik,ir,v3,v4,v1,ratio,mulion(ik,ir),
c     .      ktebs(ik,ir),knbs(ik,ir)


        ENDDO
      ENDDO

c      STOP 'CALCULATING MULTIPLIER'

      RETURN
99    CONTINUE
      STOP
      END
c
c ======================================================================
c

c
c ======================================================================
c
c subroutine: BalanceGrid
c
      SUBROUTINE BalanceGrid
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER i1,ik,ir,in1,in2,irref,nring
      LOGICAL output

      output = .FALSE.

      IF (output) WRITE(0,*) 'BALANCING GRID'

      CALL OutputGrid(85,'Before balancing grid')
c
c
c
      IF (irsep-2.GT.nrs-irtrap) THEN
        irref = irtrap + 1
        nring = (irsep - 2) - (nrs - irtrap)
      ELSE
        irref = 2
        nring = (nrs - irtrap) - (irsep - 2)
      ENDIF
c
c
c
      DO i1 = 1, nring
        CALL InsertRing(irref,BEFORE,TEMPORARY)
        idring(irref) = -2
      ENDDO
c
c
c
      vpolmin = (MAXNKS*MAXNRS - npolyp) / 2 + npolyp
      vpolyp  = vpolmin
c
c
c
      DO ir = irref+nring-1, irref, -1

        IF (output) WRITE(0,*) 'ASSIGNING RING = ',ir

        DO ik = 1, nks(ir)
          IF (vpolyp.EQ.MAXNKS*MAXNRS)
     .      CALL ER('BalanceGrid','Out of bounds',*99)

          vpolyp = vpolyp + 1
          in1    = vpolyp
          in2    = korpg(ik,ir+1)

          korpg(ik,ir) = in1

          rvertp(1,in1) = rvertp(1,in2)
          zvertp(1,in1) = zvertp(1,in2)
          rvertp(2,in1) = rvertp(1,in2)
          zvertp(2,in1) = zvertp(1,in2)

          rvertp(3,in1) = rvertp(4,in2)
          zvertp(3,in1) = zvertp(4,in2)
          rvertp(4,in1) = rvertp(4,in2)
          zvertp(4,in1) = zvertp(4,in2)

          rs(ik,ir) = 0.25 * (rvertp(1,in1) + rvertp(2,in1) +
     .                        rvertp(3,in1) + rvertp(4,in1))
          zs(ik,ir) = 0.25 * (zvertp(1,in1) + zvertp(2,in1) +
     .                        zvertp(3,in1) + zvertp(4,in1))
        ENDDO

      ENDDO


      CALL OutputGrid(86,'Grid balanced')


      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: DivAddSurface
c
c Assumes as 4 sided surface.
c
c
c
c
c
      SUBROUTINE DivAddSurface(x,y,z,mode,nlimi3,fp3,colour)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'


      INTEGER mode,nlimi3,fp3
      REAL    x(4),y(4),z(4)


      INTEGER localmode,segmin,segmax,i1,i2,i3,orientation,index(4),
     .        istart,icount,status,itemp,iseg,colour
      REAL    x1(4),y1(4),z1(4),lencyl,phiseg,zmin,zmax,phimin,phimax,
     .        a1,a2,a3,b1,b2,b3,c1,c2,c3,zshift,lenseg,x0,phi,x2(4),
     .        y2(4),z2(4),xsave

      SAVE

c...  Default:
c      lencyl = 4.7

      IF (eirzaa.GE.0.0.OR.eirzaa.EQ.-1.0) 
     .  CALL ER('DivAddSurface','Invalid EIRZAA value',*99)

      lencyl = -eirzaa

      IF (mode.EQ.-1) THEN

c...    Determine the number of segments spanned by the additional surface:
    
        zmin =  HI
        zmax = -HI
        DO i1 = 1, 4
          zmin = MIN(zmin,z(i1) / 100.0)
          zmax = MAX(zmax,z(i1) / 100.0)
        ENDDO

        localmode = 0

        phiseg = 360.0 / REAL(eirntorseg)
        lenseg = 100.0 * lencyl / REAL(eirntorseg)

        IF (zmin.LT.-1.0E+03.OR.zmax.GT.1.0E+03) THEN
c...      Surface is assumed to be toroidally continuous:
          mode = 1
          localmode = -1
        ELSE
          phimin = zmin / lencyl * 360.0 + 0.5 * phiseg
          phimax = zmax / lencyl * 360.0 + 0.5 * phiseg

          segmin = INT(phimin/phiseg) + MAX(0,NINT(SIGN(1.0,phimin)))
          segmax = INT(phimax/phiseg) + MAX(0,NINT(SIGN(1.0,phimax)))

c          WRITE(0,*) 'phimin,phiseg',phimin,segmin,phiseg
c          WRITE(0,*) 'phimax       ',phimax,segmax

c...      Return the number of sections the surface will be divided into:
          mode = segmax - segmin + 1

c...      Store parent surface coorindates:
          DO i1 = 1, 4
            x1(i1) = x(i1)
            y1(i1) = y(i1)
            z1(i1) = z(i1)
          ENDDO 

          icount = 0

        ENDIF

        RETURN

      ELSEIF (localmode.EQ.-1) THEN
c...    The surface was determined to be toroidally continuous on the
c       previous call to this routine:
        mode = 0

      ELSE

c...    Figure out which segment the surface is being cropped to:
        iseg = segmin + icount

c...    Determine the orientation of the surface:

c...    Take the cross product to get the surface normal:
        a1 = x1(1) - x1(2)
        a2 = y1(1) - y1(2)
        a3 = z1(1) - z1(2)
	
        b1 = x1(3) - x1(2)
        b2 = y1(3) - y1(2)
        b3 = z1(3) - z1(2)

        c1 = a2 * b3 - a3 * b2                
        c2 = a3 * b1 - a1 * b3
        c3 = a1 * b2 - a2 * b1                

c...    Examine the surface normal to decide if the surface is
c       parallel or perpendicular to the toroidal axis:
        IF     (ABS(c3).LT.1.0E-6) THEN
c...      Perpendicular:
          IF     (c2.GT. 1.0E-6) THEN
c...        Points ordered clockwise:
            orientation = 1
          ELSEIF (c2.LT.-1.0E-6) THEN
c...        Points ordered counterclockwise:
            orientation = 2
          ELSEIF (c1.GT. 1.0E-6) THEN
            orientation = 3
          ELSEIF (c1.LT.-1.0E-6) THEN
            orientation = 4
          ELSE
            CALL ER('DivAddSurface','Unrecognized surface '//
     .              'orientation',*99)
          ENDIF
        ELSEIF (ABS(c1).LT.1.0E-6.AND.ABS(c2).LT.1.0E-6) THEN
c...      Parallel:
          orientation = 5
        ENDIF

c...    Standardize the order of the points, if necessary:
c
c           3  4
c
c           2  1
c
        IF (orientation.NE.5) THEN
c...      Surface normal is perpendicular to the cylindrical axis:

c...      Find starting point - the Xmax,Zmin point (or Ymax,Zmin if the 
c         surface normal is parallel to the x-axis):

c...      Order the points from lowest to highest z-value:             
          status = 0
          DO i1 = 1, 4
            index(i1) = i1
          ENDDO
          DO WHILE (status.EQ.0)
            status = 1
            DO i1 = 1, 3
              IF (z1(index(i1)).GT.z1(index(i1+1))) THEN
                itemp = index(i1)
                index(i1) = index(i1+1)
                index(i1+1) = itemp
                status = 0
              ENDIF
            ENDDO
          ENDDO

          IF     (orientation.EQ.1.OR.orientation.EQ.2) THEN
            IF (x1(index(1)).GT.x1(index(2))) THEN
              istart = index(1)
            ELSE
              istart = index(2)
            ENDIF
          ELSEIF (orientation.EQ.3.OR.orientation.EQ.4) THEN
            IF (y1(index(1)).GT.y1(index(2))) THEN
              istart = index(1)
            ELSE
              istart = index(2)
            ENDIF
          ENDIF

          IF (orientation.EQ.1.OR.orientation.EQ.3) THEN
            i3 = 0
            DO i1 = istart, istart + 3
              i2 = i1
              IF (i2.GT.4) i2 = i2 - 4
              i3 = i3 + 1
              x(i3) = x1(i2)
              y(i3) = y1(i2)
              z(i3) = z1(i2)
            ENDDO
          ELSE
c...        Reverse the order of the points:
            i3 = 0
            DO i1 = istart, istart-3, -1
              i2 = i1
              IF (i2.LT.1) i2 = i2 + 4
              i3 = i3 + 1
              x(i3) = x1(i2)
              y(i3) = y1(i2)
              z(i3) = z1(i2)
            ENDDO
          ENDIF
        ENDIF


c        WRITE(0,*) 'seg:',iseg,segmin,segmax
c        WRITE(0,*) 'orientation',orientation,istart

c        WRITE(0,*) '--- STAR ---'
c        WRITE(0,*) (x(i1),i1=1,4)
c        WRITE(0,*) (y(i1),i1=1,4)
c        WRITE(0,*) (z(i1),i1=1,4)

c...    Translate the surface so that the region to be cropped out
c       is located in the first segment, if necessary:
        IF (segmin.NE.1.OR.segmax.NE.1) THEN
          zshift = -100.0 * lencyl / REAL(eirntorseg) * REAL(iseg - 1)
          DO i1 = 1, 4
            z(i1) = z(i1) + zshift
          ENDDO
        ENDIF

c        WRITE(0,*) '--- TRAN ---',zshift
c        WRITE(0,*) (x(i1),i1=1,4)
c        WRITE(0,*) (y(i1),i1=1,4)
c        WRITE(0,*) (z(i1),i1=1,4)

c...    Crop the surface, if necessary:  
        IF (segmax.GT.segmin) THEN
          DO i1 = 1, 4
            z(i1) = MIN(MAX(z(i1),-0.5*lenseg),0.5*lenseg)
          ENDDO
        ENDIF

c        WRITE(0,*) '--- CROP ---',0.5*lenseg
c        WRITE(0,*) (x(i1),i1=1,4)
c        WRITE(0,*) (y(i1),i1=1,4)
c        WRITE(0,*) (z(i1),i1=1,4)

c...    Convert to toroidal coordinates:
        x0 = 0.5 * lenseg / TAN(0.5 * phiseg * DEGRAD) 
c        WRITE(0,*) 'x0:',x0
        DO i1 = 1, 4
          z(i1) = z(i1) * (x(i1) / x0)
        ENDDO


c...    Determine ILTOR:
        mode = iseg
        IF (mode.LT.1         ) mode = mode + eirntorseg
        IF (mode.GT.eirntorseg) mode = mode - eirntorseg

c...    Adjust the orientation of the surface so that it matches
c       the orientation of the parent surface, if necessary:
        IF (orientation.EQ.2.OR.orientation.EQ.4) THEN
c...      Reverse the order of the points:
          DO i1 = 1, 4
            x2(i1) = x(i1)
            y2(i1) = y(i1)
            z2(i1) = z(i1)
          ENDDO
          i3 = 0
          DO i1 = 4, 1, -1
            i3 = i3 + 1
            x(i3) = x2(i1)
            y(i3) = y2(i1)
            z(i3) = z2(i1)
          ENDDO
        ENDIF



c...    Output:

c...    Rotate through phi:
        phi = phiseg * REAL(mode - 1) * DEGRAD 
        DO i1 = 1, 4
          x2(i1) = COS(phi) * x(i1) - SIN(phi) * z(i1)
          y2(i1) = y(i1)
          z2(i1) = SIN(phi) * x(i1) + COS(phi) * z(i1)
        ENDDO
	
c...    Send geometry data to .dump file:
        IF (fp3.NE.99) THEN
          nlimi3 = nlimi3 + 1
          DO i2 = 1, 4
            i3 = i2 + 1
            IF (i3.GT.4) i3 = 1
            WRITE(fp3,'(I6,2(3E14.6,2X),I6)') 
     .        nlimi3,x2(i2)/100.0,y2(i2)/100.0,z2(i2)/100.0,
     .               x2(i3)/100.0,y2(i3)/100.0,z2(i3)/100.0,colour
          ENDDO
        ENDIF

c...    Update counter:
        icount = icount + 1
      ENDIF


      IF (localmode.EQ.-1) THEN

c...    Crop:
        DO i1 = 1, 4
          x2(i1) = x(i1)
          y2(i1) = y(i1)
          z2(i1) = MIN(MAX(z(i1),-0.5*lenseg),0.5*lenseg)
        ENDDO
	
c...    Convert to toroidal coordinates:
        x0 = 0.5 * lenseg / TAN(0.5 * phiseg * DEGRAD) 
        DO i1 = 1, 4
          z2(i1) = z2(i1) * (x2(i1) / x0)
        ENDDO

        DO iseg = 1, eirntorseg

          IF (iseg.EQ.1) THEN
            phi = 0.0
          ELSE
            phi = phiseg * DEGRAD 
          ENDIF
       
          DO i1 = 1, 4
            xsave = x2(i1)
            x2(i1) = COS(phi) * xsave - SIN(phi) * z2(i1)
            z2(i1) = SIN(phi) * xsave + COS(phi) * z2(i1)
          ENDDO

          IF (fp3.NE.99) THEN
            nlimi3 = nlimi3 + 1
            DO i2 = 1, 4
              i3 = i2 + 1
              IF (i3.GT.4) i3 = 1
              WRITE(fp3,'(I6,2(3E14.6,2X),I6)') 
     .          nlimi3,x2(i2)/100.0,y2(i2)/100.0,z2(i2)/100.0,
     .                 x2(i3)/100.0,y2(i3)/100.0,z2(i3)/100.0,-colour
            ENDDO
          ENDIF
        ENDDO

      ELSE

      ENDIF

      
      RETURN
99    WRITE(0,*) 'XYZ:',x,y,z
      STOP
      END
c
c ======================================================================
c ======================================================================
c
c Routines for writing plasma data in the standard B2-EIRENE format.
c
c
c
c ======================================================================
c
c
c
      subroutine maptob2(cutring,cutpt1,cutpt2,nx,ny,ndimx,ndimy,
     >                    ndims,dummy,divarr,dim1,dim2,scalef,valtype)
      implicit none
c
c     Include params and cgeom to get at geometry data needed for
c     mapping cell edge quantites.
c
      include 'params'
      include 'cgeom'
c slmod begin
      include 'slcom'
 
      INTEGER ik1
c slmod end
c
c     MAPTOB2:
c
c     transform divimp-array divarr to b2-format (dummy)
c     scalef is a scale factor
c     dummy = divarr / scalef
c
c     This routine will also convert a cell-centered DIVIMP quantity to
c     a cell-edge value by loading the average weighted value from two
c     adjacent DIVIMP cells.
c
c
      integer valtype
c
      integer cutring,cutpt1,cutpt2
      integer nx,ny,ndimx,ndimy,ndims,dim1,dim2
      integer ir,ik,mrings,mkpts,ikold,irold,ikact,iract
      real divarr(dim1,dim2),scalef
      real dummy(0:ndimx+1,0:ndimy+1,1:ndims)
c
      mkpts = nx + 2
      mrings= ny + 2
c
c
c     Handle cell-centered quantities
c
      if (valtype.eq.0) then
c
         do 100 ir = 1,mrings
            do 100 ik = 1,mkpts
               if (ir.le.cutring) then
                  if (ik.le.cutpt1) then
                     iract = mrings + ir
                     ikact = ik
                  elseif (ik.ge.cutpt2) then
                     iract = mrings + ir
                     ikact = ik - cutpt2 + cutpt1 +1
                  else
                     iract = ir
                     ikact = ik - cutpt1
                  endif
               else
                  iract = ir
                  ikact = ik
               endif
              dummy(ik-1,ir-1,1) = divarr(ikact,iract)/scalef
c               write (SLOUT,'(A,4I6,1P,2G10.2,0P,2F10.4)') 
c     .                  'maptob2 0:',ik,ir,ikact,iract,
c     >                  dummy(ik-1,ir-1,1),divarr(ikact,iract),
c     >                  kss(ikact+1,iract),kss(ikact,iract)
 100     continue
c
c     Handle cell-edge based quantities - by averaging cell-centre values.
c
c     This code assumes that additional values have been loaded into
c     the virtual cells - so that when averaged with the first real
c     cells they will recreate the actualt target value.
c
c
      elseif (valtype.eq.1) then
c
c        Locate matching DIVIMP cell
c
         do 200 ir = 1,mrings
            do 200 ik = 1,mkpts
c slmod begin
              IF (.FALSE..AND.grdnmod.NE.0) THEN
c...            Generalized grid mapping:
                if (ir.le.cutring) then
                   if (ik.le.ikto2(mrings+ir)+1) then
                      iract = mrings + ir
                      ikact = ik
                   elseif (ik.ge.ikti2(mrings+ir)+1+nks(ir)-1) then
                      iract = mrings + ir
                      ikact = ik - nks(ir) + 1
                   else
                      iract = ir
                      ikact = ik - ikto2(mrings+ir) - 1
                   endif
                else
                   iract = ir
                   ikact = ik
                endif
              ELSE
c...            Standard index mapping:
                if (ir.le.cutring) then
                   if (ik.le.cutpt1) then
                      iract = mrings + ir
                      ikact = ik
                   elseif (ik.ge.cutpt2) then
                      iract = mrings + ir
                      ikact = ik - cutpt2 + cutpt1 +1
                   else
                      iract = ir
                      ikact = ik - cutpt1
                   endif
                else
                   iract = ir
                   ikact = ik
                endif
              ENDIF
c
c               if (ir.le.cutring) then
c                  if (ik.le.cutpt1) then
c                     iract = mrings + ir
c                     ikact = ik
c                  elseif (ik.ge.cutpt2) then
c                     iract = mrings + ir
c                     ikact = ik - cutpt2 + cutpt1 +1
c                  else
c                     iract = ir
c                     ikact = ik - cutpt1
c                  endif
c               else
c                  iract = ir
c                  ikact = ik
c               endif
c slmod end
c
c              Test special cases and assign values
c
c              Maps to first - virtual cell - on DIVIMP ring
c
               if (ikact.eq.1) then
c
c                 Copy in value loaded into the virtual cell
c
                  dummy(ik-1,ir-1,1) = divarr(ikact,iract)/scalef

c
c              Maps to last - virtual cell - on DIVIMP ring
c
               elseif (ikact.eq.nks(iract)) then
c
c                 Copy in value loaded into the virtual cell
c
                  dummy(ik-1,ir-1,1) = divarr(ikact,iract)/scalef

c
c              Maps to second last - last real cell - on DIVIMP ring
c
               elseif (ikact.eq.nks(iract)-1) then
c
c                 Copy in value loaded into the virtual cell - since
c                 it has been already set to the boundary quantity.
c                 Do not do averaging.
c
                  dummy(ik-1,ir-1,1) = divarr(ikact+1,iract)/scalef
c
               else
c
c                 Perform boundary weighted average base on S.
c
c slmod begin
c
c                  IF (iract.EQ.71) THEN
c                    WRITE(0,*) 'IKACT,IRACT:',
c     .                ikact,iract,virtag(ikact,iract)
c                  ENDIF

                  IF (virtag(ikact,iract).EQ.1) THEN
                    IF (ikact.LE.ikto2(iract)+1) THEN
                      DO ik1 = ikto2(iract)+1, 2, -1
                        IF (virtag(ik1,iract).EQ.0) ikact = ik1
                      ENDDO
                    ELSE
                      DO ik1 = ikti2(iract)+1, nks(iract)-2
                        IF (virtag(ik1,iract).EQ.0) ikact = ik1 - 1
                      ENDDO
                    ENDIF
                  ENDIF
                  IF (iract.EQ.1.OR.iract.EQ.irwall.OR.
     .                iract.EQ.irtrap) THEN
                    dummy(ik-1,ir-1,1) = -1.0
                  ELSE
                    dummy(ik-1,ir-1,1) = (divarr(ikact,iract)
     >                    + ((ksb(ikact,iract)-kss(ikact,iract)) /
     >                      (kss(ikact+1,iract)-kss(ikact,iract))) *
     >                      (divarr(ikact+1,iract)-divarr(ikact,iract)))
     >                            / scalef
c                    WRITE(0,*) 'END'
                  ENDIF
c
c                  dummy(ik-1,ir-1,1) = (divarr(ikact,iract)
c     >                  + ((ksb(ikact,iract)-kss(ikact,iract)) /
c     >                    (kss(ikact+1,iract)-kss(ikact,iract))) *
c     >                   (divarr(ikact+1,iract)-divarr(ikact,iract)))
c     >                           / scalef
c
c slmod end
c
               endif

               write (SLOUT,'(A,4I6,1P,2F10.2,0P,2F10.4)') 
     .                  'maptob2 1:',ik,ir,ikact,iract,
     >                  dummy(ik-1,ir-1,1),divarr(ikact,iract),
     >                  kss(ikact+1,iract),kss(ikact,iract)

c
c               write (6,*) 'maptob2:',ik,ir,ikact,iract,
c     >                  dummy(ik-1,ir-1,1),divarr(ikact,iract),
c     >                  kss(ikact+1,iract),kss(ikact,iract)
c



 200     continue

      endif
c
      return
      end

