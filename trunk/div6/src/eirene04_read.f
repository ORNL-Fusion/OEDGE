c
c ======================================================================
c
      SUBROUTINE LoadEireneData_04
      USE mod_eirene04
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER fp,ntally,ndata,icount,index(30),ik,ir,i1,
     .        iblk,iatm,imol,iion,ipho,ilin,isur
      LOGICAL goodeof
      REAL    rdum(30),frac,
     .        sumion,amps,flux
      CHARACTER buffer*256,species*32

      REAL, ALLOCATABLE :: tdata(:,:,:)

      ALLOCATE(tdata(MAXNKS,MAXNRS,5))
      tdata = 0.0

      goodeof = .FALSE.

c      pinion = 0.0   ! TEMP... 
c      pinrec = 0.0

      pinalpha = 0.0
      pinline = 0.0
      pinatom = 0.0
      pinmol = 0.0

      hescpd = 0.0
      hescal = 0.0

      IF (rel_opt.EQ.1.OR.rel_opt.EQ.3) THEN
        frac = rel_frac
      ELSE
        frac = 1.0
      ENDIF
      IF (sloutput)
     .  WRITE(0     ,*) 'RELAXATION FRACTION FOR EIRENE04:',frac
      WRITE(PINOUT,*) 'RELAXATION FRACTION FOR EIRENE04:',frac

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
          READ(fp,*,ERR=97) ndata                         ! Check...
          READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          icount = 0
          DO WHILE (icount.LT.ndata)
            CALL NextLine(fp,ntally,icount,rdum)
            IF (tri(icount)%type.EQ.MAGNETIC_GRID) THEN  
              ik = tri(icount)%index(1)                   ! Should pull these from .transfer
              ir = tri(icount)%index(2)   
              IF (iblk.EQ.1) THEN                           ! Data for D+ only:
                tdata(ik,ir,1) = tdata(ik,ir,1) + rdum(5)   ! PINION
                tdata(ik,ir,2) = tdata(ik,ir,2) + rdum(6)   ! PINREC (relax?)
                tdata(ik,ir,3) = tdata(ik,ir,3) + rdum(7)   ! PINMP
                tdata(ik,ir,4) = tdata(ik,ir,4) + rdum(8)   ! PINQi
                tdata(ik,ir,5) = tdata(ik,ir,5) + rdum(9)   ! PINQe
              ENDIF
            ENDIF
          ENDDO

        ELSEIF (buffer(1:12).EQ.'* TEST ATOMS') THEN

          iatm = iatm + 1
          READ(fp,*,ERR=97) ntally
          READ(fp,*,ERR=97) ndata                         
          READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          icount = 0
          DO WHILE (icount.LT.ndata)
            CALL NextLine(fp,ntally,icount,rdum)
            IF (tri(icount)%type.EQ.MAGNETIC_GRID) THEN  
              ik = tri(icount)%index(1)                   
              ir = tri(icount)%index(2)   
              IF (iatm.EQ.1) THEN                            ! Only for D, presumably the 1st atom species, need check...
                pinatom(ik,ir) = pinatom(ik,ir) + rdum(1) 
              ENDIF
            ENDIF
          ENDDO

        ELSEIF (buffer(1:16).EQ.'* TEST MOLECULES') THEN

          imol = imol + 1
          READ(fp,*,ERR=97) ntally
          READ(fp,*,ERR=97) ndata                         
          READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          icount = 0
          DO WHILE (icount.LT.ndata)
            CALL NextLine(fp,ntally,icount,rdum)
            IF (tri(icount)%type.EQ.MAGNETIC_GRID) THEN  
              ik = tri(icount)%index(1)                   
              ir = tri(icount)%index(2)   
              IF (imol.EQ.1) THEN                            ! Need check...
                pinmol(ik,ir) = pinmol(ik,ir) + rdum(1)
              ENDIF
            ENDIF
          ENDDO

        ELSEIF (buffer(1:11).EQ.'* TEST IONS') THEN
        ELSEIF (buffer(1:14).EQ.'* TEST PHOTONS') THEN
        ELSEIF (buffer(1:14).EQ.'* LINE EMISSIO') THEN

          ilin = ilin + 1   
          READ(fp,*,ERR=97) ntally
          READ(fp,*,ERR=97) ndata   
          READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          icount = 0
          DO WHILE (icount.LT.ndata)
            CALL NextLine(fp,ntally,icount,rdum)
            IF (tri(icount)%type.EQ.MAGNETIC_GRID) THEN  
              ik = tri(icount)%index(1)                  
              ir = tri(icount)%index(2)   
              IF (ilin.EQ.1) THEN
                pinline(ik,ir,1:5,H_BALPHA)=pinline(ik,ir,1:5,H_BALPHA)+ 
     .                                      rdum(1:5)
                pinline(ik,ir,6  ,H_BALPHA)=pinline(ik,ir,6  ,H_BALPHA)+ 
     .                                      rdum(7)
              ELSEIF (ilin.EQ.2) THEN
                pinline(ik,ir,1:6,H_BGAMMA)=pinline(ik,ir,1:6,H_BGAMMA)+ 
     .                                      rdum(1:6)
              ENDIF  
            ENDIF
          ENDDO

        ELSEIF (buffer(1:6 ).EQ.'* MISC') THEN
c...      Check volumes:
        ELSEIF (buffer(1:13).EQ.'* PUMPED FLUX') THEN

          DO WHILE (.TRUE.) 
            READ(fp,'(A256)',END=97,ERR=97) buffer
            IF (buffer(1:6).EQ.'* DONE') THEN
              goodeof = .TRUE.
              EXIT
            ENDIF
            READ(buffer,*) isur,species,amps
            IF     (species(1:6).EQ.'D     '.OR.
     .              species(1:6).EQ.'D(N=1)') THEN
              flux = amps / ECH 
              hescpd  = hescpd + flux
              IF (isur.NE.-1) hescal = hescal + flux       ! Core ring assumption!!!!
            ELSEIF (species(1:6).EQ.'D2    ') THEN
              flux = amps / ECH * 2.0
              hescpd  = hescpd + flux
              IF (isur.NE.-1) hescal = hescal + flux
            ELSE
              CALL WN('LoadEireneData','Unknown particle type')
            ENDIF          
          ENDDO

        ELSEIF (buffer(1:1 ).EQ.'*') THEN
        ELSE
        ENDIF

      ENDDO
 10   CONTINUE

      CLOSE (fp)

c...  Relaxation into OSM arrays only?


c...  Normalize quantities (need to be careful vis-a-vis relaxation):
      sumion = 0.0
      DO ir = 1, nrs
c slmod begin
        IF (idring(ir).EQ.BOUNDARY) CYCLE
c slmod end        
        DO ik = 1, nks(ir)

c...      Volume normalization:
          tdata(ik,ir,1:5) = tdata(ik,ir,1:5) / kvols(ik,ir)

c...      Linear relaxation:
          pinion(ik,ir) = (1.0-frac)*pinion(ik,ir) + frac*tdata(ik,ir,1)
          pinrec(ik,ir) = (1.0-frac)*pinrec(ik,ir) + frac*tdata(ik,ir,2)
          pinmp (ik,ir) = (1.0-frac)*pinmp (ik,ir) + frac*tdata(ik,ir,3)
          pinqi (ik,ir) = (1.0-frac)*pinqi (ik,ir) + frac*tdata(ik,ir,4)
          pinqe (ik,ir) = (1.0-frac)*pinqe (ik,ir) + frac*tdata(ik,ir,5)

c          pinion (ik,ir) = (1.0 - frac) * pinion(ik,ir) + 
c     .                            frac  * tdata (ik,ir,1) / kvols(ik,ir)
c          pinrec (ik,ir) = (1.0 - frac) * pinrec(ik,ir) + 
c     .                                    tdata (ik,ir,2) / kvols(ik,ir)

c...      Not relaxed, volume normalize:
          pinatom(ik,ir) = pinatom(ik,ir) / kvols(ik,ir)
          pinmol (ik,ir) = pinmol (ik,ir) / kvols(ik,ir)
          pinline(ik,ir,1:6,H_BALPHA) = pinline(ik,ir,1:6,H_BALPHA) /
     .                                  kvols(ik,ir)
          pinline(ik,ir,1:6,H_BGAMMA) = pinline(ik,ir,1:6,H_BGAMMA) /
     .                                  kvols(ik,ir)

c          DO i1 = 1, 6
c            pinline(ik,ir,i1,H_BALPHA) = pinline(ik,ir,i1,H_BALPHA) /
c    .                                    kvols(ik,ir)
c            pinline(ik,ir,i1,H_BGAMMA) = pinline(ik,ir,i1,H_BGAMMA) /
c    .                                    kvols(ik,ir)
c          ENDDO
c...
          pinalpha(ik,ir) = pinline(ik,ir,6,H_BALPHA)

          sumion = sumion + pinion(ik,ir) * kvols(ik,ir)
        ENDDO
      ENDDO


      ir = irsep
      DO ik = 1, 10
        WRITE(0,*) 'CHECK:',
     .    pinline(ik,ir,6,H_BALPHA),
     .    pinline(ik,ir,6,H_BGAMMA),
     .    pinline(ik,ir,6,H_BGAMMA) /
     .    pinline(ik,ir,6,H_BALPHA)
      ENDDO


c      WRITE(0,*) 'SUMION:',sumion

      CALL OutputEIRENE(67,'WORKING WITH EIRENE04')

c...  Make sure that the whole EIRENE data file was there:
      IF (.NOT.goodeof)  CALL ER('LoadEireneData','Problem with '//
     .                           'eirene.transfer file',*99)

      DEALLOCATE(tdata)

      RETURN
 97   WRITE(0,*) 'ERROR: PROBLEM READING DATA TRANSFER FILE'
      STOP
 98   WRITE(0,*) 'WARNING: eirene.transfer DATA FILE NOT FOUND'
      RETURN
 99   STOP
      END

c
c ======================================================================
c
c subroutine: WriteEireneDatFileInfo
c
c
      SUBROUTINE WriteEireneDatFileInfo
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      REAL    GetFlux

      INTEGER fp,ik,ir
      LOGICAL first
      REAL    balance,sumion,sumrec,sumpuf,sumflx,lossc,lossp

      DATA first /.TRUE./
      SAVE

      fp = 98

      OPEN(UNIT=fp,FILE='eirene04.txt',ACCESS='SEQUENTIAL',
     .     STATUS='UNKNOWN',POSITION='APPEND',ERR=98)

c...  need to calculate the total flux to targets, recom, ionisation, balance, etc. ...


      sumion = 0.0  ! *** Also need sources/sinks that are outside the magnetic grid
      sumrec = 0.0
      sumflx = 0.0
      sumpuf = 0.0
      CALL VolInteg(pinion,3,1,nrs,sumion)
      CALL VolInteg(pinrec,3,1,nrs,sumrec)
      DO ir = irsep, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        sumflx = sumflx - GetFlux(IKLO,ir) + GetFlux(IKHI,ir)  
      ENDDO
c      sumion = sumion / eirtorfrac
c      sumrec = sumrec / eirtorfrac
c      sumflx = sumflx / eirtorfrac

      balance = (sumflx + sumrec + sumpuf) / (sumion + hescpd)

      lossc = (hescpd - hescal) / (sumflx + sumrec + sumpuf) * 100.0
      lossp = (         hescal) / (sumflx + sumrec + sumpuf) * 100.0

      IF (first) THEN 
        CALL HD(fp,'  PARTICLE BALANCE PER ITERATION','PARBAL-HD',5,85)
        WRITE(fp,*)
        WRITE(fp,'(4X,3A5,A7,4A10,2A6)')
     .    'ITER',' STEP','SUBI','BAL.','TARFLX','VOLREC','PUFF',
     .    'IONIZ','CORE','PUMP'
        WRITE(fp,'(4X,3A5,A7,4A10,2A6)')
     .    '     ','     ','     ','  ','(s-1)','(s-1)','(s-1)',
     .    '(s-1)','(%)','(%)'
        first = .FALSE.
      ENDIF

      WRITE(fp,'(4X,3I5,F7.3,1P,4E10.2,0P,2F6.1)')
     .  rel_count,rel_step,rel_iter,
     .  balance,sumflx,sumrec,sumpuf,sumion,lossc,lossp

      CLOSE (fp)

c...  Contribution to .rel.raw (awkward system):
      WRITE(79,'(A,4I6,F6.3,1P,4E10.2,0P,2F6.1)')
     .  '''PARBALAN  1.00''',rel_step,rel_iter,rel_count,7,    ! '7' indicates the number of quantites
     .  balance,sumflx,sumrec,sumpuf,sumion,lossc,lossp

      RETURN
 98   CALL ER('WriteEireneDatFileInfo','Unable to open .txt file',*99)
 99   STOP
      END
c
c ======================================================================
c
c subroutine: ReadParticleTracks
c
c
      SUBROUTINE ReadParticleTracks_04
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
c ======================================================================
c ======================================================================
c

c subroutine: DumpTriangles
c
      SUBROUTINE DumpTriangles
      USE mod_eirene04
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'


      REAL       TOL
      PARAMETER (TOL=1.0E-07)

      REAL GetMach,GetJsat,GetFlux 

      INTEGER fp,fp2,ik1,ik2,ir1,ir2,i1,i2,trin,id,nl,v1,v1n,v2,v2n,
     .        vern,in,tarside,region
      LOGICAL status
      REAL    x1,x2,y1,y2,deltax,deltay,Bfrac,fact,
     .        tarte,tarti,tarne,tarv,tarflux,tarisat,tarM
      REAL*8  Bx,By,Bz,beta,brat

      INTEGER, ALLOCATABLE :: nlink1(:,:)    ,nlink2(:,:),
     .                        llink1(:,:,:,:),llink2(:,:,:,:),
     .                        triik(:),triir(:),trisf(:),
     .                        trimap(:,:),triside(:,:),trivert(:,:),
     .                        trisurface(:,:)
      REAL, ALLOCATABLE :: trix(:,:),verx(:),verz(:),
     .                     triy(:,:),very(:)


c      RETURN

      WRITE(0,*) 'DUMPING TRIANGLES'

      ALLOCATE(nlink1(MAXNKS,MAXNRS))
      ALLOCATE(nlink2(MAXNKS,MAXNRS))
      ALLOCATE(llink1(MAXNKS,MAXNRS,100,2))
      ALLOCATE(llink2(MAXNKS,MAXNRS,100,2))

      ALLOCATE(triik(MAXNKS*MAXNRS))
      ALLOCATE(triir(MAXNKS*MAXNRS))
      ALLOCATE(trisf(MAXNKS*MAXNRS))
      ALLOCATE(trix(MAXNKS*MAXNRS,3))
      ALLOCATE(triy(MAXNKS*MAXNRS,3))
      ALLOCATE(trimap(MAXNKS*MAXNRS,3))
      ALLOCATE(triside(MAXNKS*MAXNRS,3))
      ALLOCATE(trisurface(MAXNKS*MAXNRS,3))
      ALLOCATE(trivert(MAXNKS*MAXNRS,3))

      ALLOCATE(verx(MAXNKS*MAXNRS))
      ALLOCATE(very(MAXNKS*MAXNRS))
      ALLOCATE(verz(MAXNKS*MAXNRS))

      vern = 0
      verx = 0.0
      very = 0.0
      verz = 0.0


      fp = 99
      fp2 = 98


c...  Build a list of reference quadrilaterals, and a list of connection
c     map references to each cell:    

      nlink1 = 0
      nlink2 = 0
      DO ir1 = 2, nrs
        IF (idring(ir1).EQ.-1) CYCLE
        DO ik1 = 1, nks(ir1)
          DO ir2 = 1, nrs
            IF (idring(ir2).EQ.BOUNDARY) CYCLE
            DO ik2 = 1, nks(ir2)
              IF (ir2.LT.irsep.AND.ik2.EQ.nks(ir2)) CYCLE
              IF     (irouts(ik2,ir2).EQ.ir1.AND.
     .                ikouts(ik2,ir2).EQ.ik1) THEN

                nlink1(ik1,ir1) = nlink1(ik1,ir1) + 1
                llink1(ik1,ir1,nlink1(ik1,ir1),1) = ik2
                llink1(ik1,ir1,nlink1(ik1,ir1),2) = ir2
              ELSEIF (irins (ik2,ir2).EQ.ir1.AND.
     .                ikins (ik2,ir2).EQ.ik1) THEN
                nlink2(ik1,ir1) = nlink2(ik1,ir1) + 1
                llink2(ik1,ir1,nlink2(ik1,ir1),1) = ik2
                llink2(ik1,ir1,nlink2(ik1,ir1),2) = ir2
              ENDIF
            ENDDO
          ENDDO
c...      
          IF (nlink1(ik1,ir1).LE.1) nlink1(ik1,ir1) = 0
          IF (nlink2(ik1,ir1).LE.1) nlink2(ik1,ir1) = 0
c...      Identify cells that adjacent to boundary rings:
          IF (idring(irins (ik1,ir1)).EQ.-1) nlink1(ik1,ir1) = -1
          IF (idring(irouts(ik1,ir1)).EQ.-1) nlink2(ik1,ir1) = -1
        ENDDO
c        DO ik1 = 1, nks(ir1)     
c         IF (ir1.EQ.irsep.AND.ik1.EQ.25) THEN
c           WRITE(0,*)  'DUMP:',ik1,nlink1(ik1,ir1),nlink2(ik1,ir1)
c         ENDIF
c        ENDDO
      ENDDO 

c...  Construct 2 triangles from each regular OEDGE cell:
c
c     counter-clockwise
c     Triangle 1 contains OEDGE cell verticies 1 and 2, and 
c     triangle 2 contains verticies 3 and 4.
c
c
      DO ir1 = 2, nrs
        IF (idring(ir1).EQ.-1) CYCLE

c        IF (ir1.NE.irtrap+1.AND.ir1.NE.irtrap+2) CYCLE
c        IF (ir1.NE.irsep.AND.ir1.NE.nrs.AND.ir1.NE.irsep+1.AND.
c     .      ir1.NE.irsep-1) CYCLE
c        IF (ir1.LT.irtrap) CYCLE

        DO ik1 = 1, nks(ir1)
          IF (ir1.LT.irsep.AND.ik1.EQ.nks(ir1)) CYCLE
 
          id = korpg(ik1,ir1)

          trin = trin + 1
          triik(trin) = ik1
          triir(trin) = ir1
          trisf(trin) = 23
          trix(trin,1) = rvertp(1,id)         
          triy(trin,1) = zvertp(1,id) 
          trix(trin,2) = rvertp(3,id) 
          triy(trin,2) = zvertp(3,id) 
          trix(trin,3) = rvertp(2,id) 
          triy(trin,3) = zvertp(2,id) 

          trin = trin + 1
          triik(trin) = ik1
          triir(trin) = ir1
          trisf(trin) = 14
          trix(trin,1) = rvertp(3,id) 
          triy(trin,1) = zvertp(3,id) 
          trix(trin,2) = rvertp(1,id) 
          triy(trin,2) = zvertp(1,id) 
          trix(trin,3) = rvertp(4,id) 
          triy(trin,3) = zvertp(4,id) 
        ENDDO
      ENDDO

c.... Slice triangles with more than one neighbour:
      i1 = 0
      status = .TRUE.
      DO WHILE (status)
        i1 = i1 + 1

        ik1 = triik(i1)
        ir1 = triir(i1)

        IF     (trisf(i1).EQ.23) THEN
          IF (nlink2(ik1,ir1).LE.0) CYCLE

          nl = nlink2(ik1,ir1)

c...      Make room for new triangles:
          DO i2 = trin, i1+1, -1
            triik(i2+nl-1)   = triik(i2)
            triir(i2+nl-1)   = triir(i2)
            trisf(i2+nl-1)   = trisf(i2)
            trix (i2+nl-1,1) = trix (i2,1)
            triy (i2+nl-1,1) = triy (i2,1)
            trix (i2+nl-1,2) = trix (i2,2)
            triy (i2+nl-1,2) = triy (i2,2)
            trix (i2+nl-1,3) = trix (i2,3)
            triy (i2+nl-1,3) = triy (i2,3)
          ENDDO
          trin = trin + nl - 1
  
c...      Duplicate triangle:
          DO i2 = i1+1, i1+nl-1
            triik(i2) =  triik(i1)
            triir(i2) =  triir(i1)
            trisf(i2) = -trisf(i1)
            trix(i2,1) = trix(i1,1)
            triy(i2,1) = triy(i1,1)
            trix(i2,2) = trix(i1,2)
            triy(i2,2) = triy(i1,2)
            trix(i2,3) = trix(i1,3)
            triy(i2,3) = triy(i1,3)
          ENDDO

c...      Update trailing vertex from connection map data:
          DO i2 = 2, nl
            ik2 = llink2(ik1,ir1,i2,1)
            ir2 = llink2(ik1,ir1,i2,2)
            id = korpg(ik2,ir2)
            trix(i1+i2-1,3) = rvertp(1,id) 
            triy(i1+i2-1,3) = zvertp(1,id) 
          ENDDO

c...      Make triangles in this series consistent with each subsequent
c         triangle in the series:
          DO i2 = 1, nl-1
            trix(i1+i2-1,2) = trix(i1+i2,3)  
            triy(i1+i2-1,2) = triy(i1+i2,3) 
          ENDDO

        ELSEIF (trisf(i1).EQ.14) THEN
          IF (nlink1(ik1,ir1).LE.0) CYCLE

          nl = nlink1(ik1,ir1)

c...      Make room for new triangles:
          DO i2 = trin, i1+1, -1
            triik(i2+nl-1)   = triik(i2)
            triir(i2+nl-1)   = triir(i2)
            trisf(i2+nl-1)   = trisf(i2)
            trix (i2+nl-1,1) = trix (i2,1)
            triy (i2+nl-1,1) = triy (i2,1)
            trix (i2+nl-1,2) = trix (i2,2)
            triy (i2+nl-1,2) = triy (i2,2)
            trix (i2+nl-1,3) = trix (i2,3)
            triy (i2+nl-1,3) = triy (i2,3)
          ENDDO
          trin = trin + nl - 1
  
c...      Duplicate triangle:
          DO i2 = i1+1, i1+nl-1
            triik(i2) =  triik(i1)
            triir(i2) =  triir(i1)
            trisf(i2) = -trisf(i1)
            trix(i2,1) = trix(i1,1)
            triy(i2,1) = triy(i1,1)
            trix(i2,2) = trix(i1,2)
            triy(i2,2) = triy(i1,2)
            trix(i2,3) = trix(i1,3)
            triy(i2,3) = triy(i1,3)
          ENDDO

c...      Update trailing vertex from connection map data:
          DO i2 = 2, nl
            ik2 = llink1(ik1,ir1,i2,1)
            ir2 = llink1(ik1,ir1,i2,2)
            id = korpg(ik2,ir2)
            trix(i1+i2-1,2) = rvertp(2,id) 
            triy(i1+i2-1,2) = zvertp(2,id) 
            IF (ik1.EQ.25.AND.
     .        ir1.EQ.irsep) 
     .        WRITE(0,*) 'ik,ir2=',ik2,ir2,rvertp(2,id),zvertp(2,id),nl
          ENDDO

c...      Make triangles in this series consistent with each subsequent
c         triangle in the series:
          DO i2 = 1, nl-1
            trix(i1+i2-1,3) = trix(i1+i2,2)  
            triy(i1+i2-1,3) = triy(i1+i2,2) 
          ENDDO

        ENDIF

        status = .FALSE.
        IF (i1.LT.trin) status = .TRUE.
      ENDDO


c         DO i2 = 1, trin
c           WRITE(0,'(A,I4,3F10.5)') '->',i2,trix(i2,1),triy(i2,1)
c           WRITE(0,'(A,I4,3F10.5)') '  ',i2,trix(i2,2),triy(i2,2)
c           WRITE(0,'(A,I4,3F10.5)') '  ',i2,trix(i2,3),triy(i2,3)
c         ENDDO
 


c...  Build connection map:

c...  THIS CAN BE MORE EFFICIENT:
      DO i1 = 1, trin
        DO v1 = 1, 3
          v1n = v1 + 1
          IF (v1n.EQ.4) v1n = 1
          DO i2 = 1, trin
            IF (i1.EQ.i2) CYCLE
            DO v2 = 1, 3    
              v2n = v2 - 1
              IF (v2n.EQ.0) v2n = 3
             
              IF (ABS(trix(i1,v1 )-trix(i2,v2 )).LT.TOL.AND.
     .            ABS(triy(i1,v1 )-triy(i2,v2 )).LT.TOL.AND.
     .            ABS(trix(i1,v1n)-trix(i2,v2n)).LT.TOL.AND.
     .            ABS(triy(i1,v1n)-triy(i2,v2n)).LT.TOL) THEN
                trimap (i1,v1) = i2
                triside(i1,v1) = v2n
              ENDIF

            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO i2 = 1, trin
        DO v1 = 1, 3
          IF (trimap(i2,v1).EQ.0) THEN
c...        Identify neighbourless surfaces and assign the appropriate mapping
c           to the appropriate non-default standard surface:
            ik1 = triik(i2)
            ir1 = triir(i2)
            IF     (v1.EQ.3.AND.ik1.EQ.1       .AND.ir1.GE.irsep) THEN
c...          Low IK index target segment, EIRENE poloidal suface 4:
              trisurface(i2,v1) = 2                            
c              trisurface(i2,v1) = 4                            
            ELSEIF (v1.EQ.3.AND.ik1.EQ.nks(ir1).AND.ir1.GE.irsep) THEN
c...          High IK index target segment, EIRENE poloidal surface 5:
              trisurface(i2,v1) = 3
c              trisurface(i2,v1) = 5                            
            ELSEIF (v1.EQ.2.AND.irins(ik1,ir1) .EQ.1) THEN
c...          Core boundary, EIRENE radial surface 1:
              trisurface(i2,v1) = 1                            
            ELSEIF ((v1.EQ.2.AND.irouts(ik1,ir1).EQ.irwall).OR.
c...          Added for MAST-DN secondary PFZ:
     .              (v1.EQ.2.AND.irins (ik1,ir1).EQ.irwall)) THEN
c...          SOL boundary, EIRENE radial surface 4:
              trisurface(i2,v1) = 4                            
c              trisurface(i2,v1) = 2                            
            ELSEIF (v1.EQ.2.AND.irins(ik1,ir1) .EQ.irtrap) THEN
c...          PFZ boundary, EIRENE radial surface 5:
              trisurface(i2,v1) = 5
c              trisurface(i2,v1) = 3                            
            ELSEIF (v1.EQ.1.OR.v1.EQ.3) THEN
              WRITE(0,*) 'PROBLEM 3:',i2,v1,triik(i2),triir(i2)
            ENDIF
             
          ELSEIF (v1.EQ.2.AND.ABS(trisf(i2)).EQ.23.AND.
     .            (triir(i2).EQ.irsep-1.OR.triir(i2).EQ.nrs)) THEN
c...        Label "right" side of separatrix as EIRENE non-standard surface 6:
            trisurface(i2,v1) = 6
          ELSEIF (v1.EQ.2.AND.ABS(trisf(i2)).EQ.14.AND.
     .            triir(i2).EQ.irsep) THEN
c...        Label "left" side of separatrix as EIRENE non-standard surface 7:
            trisurface(i2,v1) = 7

          ELSEIF (trimap(trimap(i2,v1),triside(i2,v1)).NE.i2) THEN
            WRITE(0,*) 'PROBLEM 4:',i2,v1,triik(i2),triir(i2)

          ELSEIF (triside(trimap(i2,v1),triside(i2,v1)).NE.v1) THEN
            WRITE(0,*) 'PROBLEM 5:',i2,v1,triik(i2),triir(i2)

          ENDIF
        ENDDO
      ENDDO

c...  All is well, so generate a master list of points with pointers to the respective
c     triangle verticies:
      DO i1 = 1, trin
        DO v1 = 1, 3

c...      Scan list of verticies and assign the relevant vertex
c         number to the triangle:
          status = .FALSE.
          DO i2 = 1, vern
            IF (ABS(verx(i2)-trix(i1,v1)).LT.TOL.AND.
     .          ABS(very(i2)-triy(i1,v1)).LT.TOL) THEN
              IF (trivert(i1,v1).EQ.0) THEN
                trivert(i1,v1) = i2
                status = .TRUE.
              ELSE
                CALL ER('DumpTriangles','Redundant vertex',*99)
              ENDIF
            ENDIF
          ENDDO
c...      Vertex was not found, so add it to the list:
          IF (.NOT.status) THEN
            vern = vern + 1
            verx(vern) = trix(i1,v1)
            very(vern) = triy(i1,v1)
            trivert(i1,v1) = vern
          ENDIF

        ENDDO
      ENDDO

c      WRITE(0,*) 'vern:',vern,trin

c...  Dump triangles:
      OPEN(UNIT=fp,FILE='triangles.dat',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
      WRITE(fp,*) trin
      DO i1 = 1, trin
        WRITE(fp,'(I6,3(2F10.6,2X))')
     .    i1,(verx(trivert(i1,i2)),very(trivert(i1,i2)),i2=1,3)
c        WRITE(fp,'(I6,3(2F10.6,2X))')
c     .    i1,(trix(i1,i2),triy(i1,i2),i2=1,3)
      ENDDO
      CLOSE(fp)      

c...  Dump triangles:
      OPEN(UNIT=fp,FILE='triangles.points',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
      WRITE(fp,*) vern
      DO i1 = 1, vern
        WRITE(fp,'(I6,3F12.6)') i1,verx(i1)*100.0,very(i1)*100.0,0.0
      ENDDO
      CLOSE(fp)      

c...  Dump triangles:
      OPEN(UNIT=fp,FILE='triangles.sides',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
      WRITE(fp,*) trin
      DO i1 = 1, trin
        WRITE(fp,'(I6,4X,3I6)') i1,(trivert(i1,v1),v1=1,3)
      ENDDO
      CLOSE(fp)      

c...  Dump triangles:
      OPEN(UNIT=fp,FILE='triangles.map',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
      WRITE(fp,*) trin
      DO i1 = 1, trin
        WRITE(fp,'(I6,4X,3(3I6,4X),2I6)') i1,
     .    (trimap(i1,v1),triside(i1,v1),trisurface(i1,v1),v1=1,3),
     .    triik(i1),triir(i1)
      ENDDO
      CLOSE(fp)      


      fact = qtim * qtim * emi / crmi

c...  Dump plasma data:
      OPEN(UNIT=fp ,FILE='triangles.plasma',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
      OPEN(UNIT=fp2,FILE='triangles.efield',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      

c...  Header:
      WRITE(fp,'(A)') '* VERSION 1.0 OSM PLASMA FILE FOR TRIANGULAR '//
     .                'GRID'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A)') '* BULK PLASMA DATA'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A7,2A8,2X,A10,2X,3A10,2X,3A12)') 
     .  '* Index','Te','Ti','ne','vx',
     .  'vy','vz','Bx','By','Bz'
      WRITE(fp,'(A7,2A8,2X,A10,2X,3A10,2X,3A12)') 
     .  '*      ','(eV)','(eV)','(cm-3)','(cm s-1)',
     .  '(cm s-1)','(cm s-1)','(Tesla)','(Tesla)','(Tesla)'


c...  Header:
      WRITE(fp2,'(A)') '* VERSION 1.0 OSM PLASMA FILE FOR TRIANGULAR '//
     .                'GRID'
      WRITE(fp2,'(A)') '*'
      WRITE(fp2,'(A)') '* BULK PLASMA DATA (PART DEUX)'
      WRITE(fp2,'(A)') '*'
      WRITE(fp2,'(A7,3A12)') 
     .  '* Index','Ex','Ey','Ez'
      WRITE(fp2,'(A7,3A12)')
     .  '*      ','V m-1','V m-1','V m-1'

      
      WRITE(fp,*) trin
      DO i1 = 1, trin
        ik1 = triik(i1)
        ir1 = triir(i1)
        id = korpg(ik1,ir1)
        x1 = 0.5 * (rvertp(1,id) + rvertp(2,id))
        y1 = 0.5 * (zvertp(1,id) + zvertp(2,id))
        x2 = 0.5 * (rvertp(3,id) + rvertp(4,id))
        y2 = 0.5 * (zvertp(3,id) + zvertp(4,id))
        deltax = (x2 - x1)
        deltay = (y2 - y1)
        brat = DBLE(bratio(ik1,ir1))
        beta = DBLE(deltax / deltay)
        Bz = DSQRT(1.0 - brat**2)
        By = brat * DSQRT(1.0D0 / (1.0D0 + beta**2)) * 
     .       DBLE(SIGN(1.0,deltay))
        Bx = beta * By
c...    CBPHI is the on-axis B-field value specified in the OSM input file:
        Bfrac = cbphi * r0 / rs(ik1,ir1) 

        WRITE(fp,'(I7,2F8.2,2X,1P,E10.2,2X,3E10.2,2X,
     .             3E12.4,0P,6X,2I4)') i1,
c     .             3E12.4,0P,4X,2I6,3F10.5)') i1,
c        WRITE(fp,'(I6,2F10.2,1P,E10.2,2X,5E10.2,0P,2I6,2F10.2))') i1,
     .         ktebs(ik1,ir1),ktibs(ik1,ir1),knbs(ik1,ir1)*1.0E-06,
     .         SNGL(Bx)*kvhs(ik1,ir1)*100.0,  !/qtim,
     .         SNGL(By)*kvhs(ik1,ir1)*100.0,  !/qtim,
     .         SNGL(Bz)*kvhs(ik1,ir1)*100.0,  !/qtim,
     .         SNGL(Bx)*Bfrac,SNGL(By)*Bfrac,SNGL(Bz)*Bfrac,
c...    Target quantities:
c     .         tarside,tarflux,tarte,tarti,tarne*1.0E-06,tarv*100.0,
c     .         tarM,tarisat*1.0E-04, 
c...    Grid cell where triangle originated
     .         ik1,ir1
c     .         SQRT(SNGL(Bx)**2+
c     .              SNGL(By)**2+
c     .              SNGL(Bz)**2),bratio(ik1,ir1),
c     .         SQRT(SNGL(Bx)**2+
c     .              SNGL(By)**2)
c     .         SQRT((SNGL(Bx)*kvhs(ik1,ir1))**2+
c     .              (SNGL(By)*kvhs(ik1,ir1))**2+
c     .              (SNGL(Bz)*kvhs(ik1,ir1))**2)*
c     .              sign(1.0,kvhs(ik1,ir1)),
c     .         kvhs(ik1,ir1),
c     .         ik1,ir1,
c     .         bratio(ik1,ir1),SNGL(DSQRT(Bx**2+By**2)/Bz)
cc     .         SNGL(Bx),SNGL(By),SNGL(Bz),
cc     .         DSQRT(Bx**2+By**2+Bz**2),ik1,ir1,
cc     .         bratio(ik1,ir1),SNGL(DSQRT(Bx**2+By**2)/Bz)


c...    Dump efield data:
        WRITE(fp2,'(I7,1P,3E12.4,10X,0P,2I4,1P,E12.4)') i1,
     .         SNGL(Bx)*kes(ik1,ir1)/fact,
     .         SNGL(By)*kes(ik1,ir1)/fact,
     .         SNGL(Bz)*kes(ik1,ir1)/fact,
     .         ik1,ir1,kes(ik1,ir1)/fact

      ENDDO

c      WRITE(fp,'(A)') '*'
c      WRITE(fp,'(A)') '* TARGET DATA'
      WRITE(fp,'(A)') '*'
c      WRITE(fp,'(A7,A6,A10,2A8,2A10,A6,A10') 
c     .  '* Index','Side','flux','Te','Ti','ne','v','Mach','Isat'
c      WRITE(fp,'(A7,A6,A10,2A8,2A10,A6,A10') 
c     .  '*      ','','(A)','(eV)','(eV)','(cm-3)','(cm s-1)','',
c     .  '(A cm-2)'
      
      WRITE(fp,*) 2*(nrs-irsep-1)

      DO i1 = 1, trin
        ik1 = triik(i1)
        ir1 = triir(i1)
        IF (ik1.NE.1.AND.ik1.NE.nks(ir1).OR.ir1.LT.irsep) CYCLE
        DO v1 = 1, 3
          IF (v1.NE.3.OR.trisurface(i1,v1).EQ.0) CYCLE
c...      Process target data:
          IF (ik1.EQ.1) THEN
            in = idds(ir1,2)
            region = IKLO
          ELSE
            in = idds(ir1,1)
            region = IKHI
          ENDIF
          tarside = v1
          tarte = kteds(in)
          tarti = ktids(in)
          tarne = knds(in)
          tarv = kvds(in)
          tarM = GetMach(tarv,tarte,tarti)
          tarisat = GetJsat(tarte,tarti,tarne,tarv)
          tarflux = GetFlux(region,ir1) / eirtorfrac / eirsrcmul * ECH

          WRITE(fp,'(I7,I6,1P,E10.2,0P,2F8.2,1P,2E10.2,0P,
     .               F6.2,1P,E10.2,0P,6X,2I4)') 
c...    Target quantities:
     .           i1,tarside,tarflux,tarte,tarti,tarne*1.0E-06,
     .           tarv*100.0,tarM,tarisat*1.0E-04, 
c...    Grid cell where triangle originated
     .           ik1,ir1
        ENDDO
      ENDDO

c...  All done:
      CLOSE(fp)      
      CLOSE(fp2)      

c...  Store triangle list in case opacity data is read from a separate
c     EIRENE run:
      IF (trin.GT.MAXNKS*MAXNRS) 
     .  CALL ER('DumpTrinagles','Arrays bounds error',*99)
      eirtrin = trin
      DO i1 = 1, trin
        eirtriik(i1) = triik(i1)
        eirtriir(i1) = triir(i1)
      ENDDO

c...  mod_eirene04:
      ntri = trin
      nver = vern

      DEALLOCATE(nlink1)
      DEALLOCATE(nlink2)
      DEALLOCATE(llink1)
      DEALLOCATE(llink2)

      DEALLOCATE(triik)
      DEALLOCATE(triir)
      DEALLOCATE(trix)
      DEALLOCATE(triy)
      DEALLOCATE(trimap)
      DEALLOCATE(triside)
      DEALLOCATE(trisurface)
      DEALLOCATE(trivert)

      DEALLOCATE(verx)
      DEALLOCATE(very)
      DEALLOCATE(verz)

c      STOP 'DUMPING'
      WRITE(0,*) 'DONE'

      RETURN
96    WRITE(0,*) 'DUMPTRIANGES: PROBLEMS WITH FILE ACCESS'
      STOP
99    STOP
      END
