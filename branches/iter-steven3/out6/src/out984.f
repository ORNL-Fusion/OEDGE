c     -*-Fortran-*-
c
c ======================================================================
c
c subroutine: Plot984
c
c 2D reaction rate plots
c
c
      SUBROUTINE Plot984(job,graph,ref,title,iopt,
     .                   xxmin,xxmax,yymin,yymax,ft,fp,zadj,
     .                   ismoth,ignors,itec,avs,navs,zval)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'dynam2'
      INCLUDE 'pindata'
      INCLUDE 'comgra'
      INCLUDE 'colours'
      INCLUDE 'reiser_com'
      INCLUDE 'slcom'
      INCLUDE 'slout'

      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost

      COMMON /NSMOOTH/ NUMSMOOTH,cgrprint
      INTEGER NUMSMOOTH,cgrprint

      COMMON /DUMCOM/ map1x_d,map2x_d,map1y_d,map2y_d,tag_d
      INTEGER         tag_d
      REAL            map1x_d,map2x_d,map1y_d,map2y_d

      COMMON /PLOT982COM/ hardscale,colscale
      LOGICAL             hardscale,colscale

c...  General:
      INTEGER       iopt,FP,FT,ISMOTH,J,id,i1,i2,i,ik,ir,nconts,
     .              IGNORS(MAXNGS),ITEC,NAVS
      LOGICAL       status,grayscale
      REAL          AVS(0:100),XXMIN,XXMAX,YYMIN,YYMAX,zadj
      CHARACTER     TITLE*(*),JOB*72,GRAPH*80,graph6*128
      CHARACTER*36  XLAB,REF,NVIEW,PLANE,ANLY,TABLE
      CHARACTER*72  YLAB,SMOOTH
      CHARACTER*128 cdum1


c...  984:
      REAL GetEAD

      INTEGER    MAXCELL      ,MAXPOINT
      PARAMETER (MAXCELL=20000,MAXPOINT=100)

      INTEGER   nc,nv(MAXCELL),shade,shift,set,avoid(MAXASCDAT),idum1,
     .          istart,iend,i3,cut1,size,line,component,iz
      LOGICAL   vacuum,standard,inside
      REAL      rv(MAXPOINT,MAXCELL),zv(MAXPOINT,MAXCELL),cq(MAXCELL),
     .          x1,y1,qmin,qmax,frac,spot,dspot,range1,range2,zval,xval,
     .          xmin,ymin,xmax,ymax,red,green,blue,frac5,fmod5,zlim,
     .          rlim,ebeam,vbeam,sigmav,sigmavb,lmin,lmax,vavg,tD2,tD,
     .          decade,decval,y2,val,ti,ne,t1a,t2a,xpos,ypos,bright,hue,
     .          nD,nD2,nD2p,sigmav1,sigmav2,MARsum,cellvol
      CHARACTER label*128


c...  Contour plot option:
      LOGICAL centerinside
      REAL*8 t1,t2
c
c     jdemod - added the allocatable attribute to the allocatable arrays below
c
      real,allocatable :: image(:,:),image1(:,:),raxis(:),zaxis(:)  
      integer,allocatable :: tpt(:,:),tpr(:,:)
c
      real uconts(maxpts),rdum1
      integer icntr,ncntr,xres,yres,maxix,maxiy,nix,niy,ic,
     .        icount,npt,iplot,rq(MAXCELL)
      real minscale,maxscale,dpt(100),vpt(100),tmpdpt,tmpvpt,dist,totdpt
      real xcen,ycen,xnear,ynear,rc(MAXCELL),zc(MAXCELL)
      character*36 blabs(2)
      CHARACTER*128 elabs(MAXNGS)  
      CHARACTER dummy*5000,caption*5000,graph3*80

c...  ADAS:
      CHARACTER ADASID*80,adasex*3
      integer   adasyr,ISELE,ISELR,ISELX,iseld,ierr,ircode,izmin
      REAL      wlngth,PLRPAD(MAXNKS,MAXNRS)


      call setup_col(10,5)

      WRITE(6,*) '984:'

c      CALL THICK2(4)

      CALL IZero(avoid,MAXASC)
      CALL RZero(cq   ,MAXCELL)
      CALL IZero(rq   ,MAXCELL)

      vacuum   = .FALSE.
      standard = .FALSE.

      MARsum = 0.0

c      slopt = 1

      grayscale = .FALSE.

c...  Individual plot setup:

      IF     (iopt.EQ.1) THEN
c...    p + H2 -> p+ H2 mean free path:
        standard = .TRUE.
        vacuum   = .TRUE.
        READ(5,'(A80)') dummy
        IF   (dummy(8:11).EQ.'Beam'.OR.dummy(8:11).EQ.'BEAM'.OR.
     .        dummy(8:11).EQ.'beam') THEN
          READ(dummy,*) cdum1,ebeam
        ELSE
          CALL ER('984','Expecting beam energy data line',*99)
        ENDIF
        IF (ebeam.EQ.99.0) THEN
          vbeam = 99.0
        ELSE
          vbeam = SQRT(2.0 * ebeam * ECH / (2.0 * crmb * AMU))
        ENDIF
        WRITE(char(28),'(A,F8.3,A)') 'Ebeam =',ebeam,' eV'
        WRITE(char(29),'(A       )') 'UNITS =   m'
        WRITE(char(30),'(A,I4    )') 'STEP   ',rel_step
        CALL LoadEIRENEAtomicData
      ELSEIF (iopt.EQ.2) THEN
c...    e + H2 -> ... mean free path:
        standard = .TRUE.
        vacuum   = .TRUE.

        READ(5,'(A80)') dummy
        IF   (dummy(8:11).EQ.'Beam'.OR.dummy(8:11).EQ.'BEAM'.OR.
     .        dummy(8:11).EQ.'beam') THEN
          READ(dummy,*) cdum1,ebeam
        ELSE
          BACKSPACE 5
          ebeam = 99.0
        ENDIF

        WRITE(char(28),'(A,F8.3,A)') 'Ebeam =',ebeam,' eV'
        WRITE(char(29),'(A   )') 'UNITS =   m'
        WRITE(char(30),'(A,I4)') 'STEP   ',rel_step
        CALL LoadEIRENEAtomicData
      ELSEIF (iopt.EQ.3) THEN
c...    e + H -> 2E + H+ mean free path:
        standard = .TRUE.
        vacuum   = .TRUE.
        WRITE(char(29),'(A   )') 'UNITS =   m'
        WRITE(char(30),'(A,I4)') 'STEP   ',rel_step
        CALL LoadEIRENEAtomicData

      ELSEIF (iopt.EQ.4) THEN
c...    p + H2(v) -> H + H2+ MFP:
        standard = .TRUE.
        vacuum   = .TRUE.
        READ(5,'(A80)') dummy
        IF   (dummy(8:11).EQ.'Beam'.OR.dummy(8:11).EQ.'BEAM'.OR.
     .        dummy(8:11).EQ.'beam') THEN
          READ(dummy,*) cdum1,ebeam
        ELSE
          CALL ER('984','Expecting beam energy data line',*99)
        ENDIF
        IF (ebeam.EQ.99.0) THEN
          vbeam = 99.0
        ELSE
          vbeam = SQRT(2.0 * ebeam * ECH / (2.0 * crmb * AMU))
        ENDIF
        WRITE(char(28),'(A,F8.3,A)') 'Ebeam =',ebeam,' eV'
        WRITE(char(29),'(A       )') 'UNITS =   m'
        WRITE(char(30),'(A,I4    )') 'STEP   ',rel_step
        CALL LoadEIRENEAtomicData

      ELSEIF (iopt.EQ.5) THEN
c...    e + H2+ -> H + H (MAR):
        standard = .TRUE.
        vacuum   = .TRUE.

        READ(5,'(A80)') dummy
        IF   (dummy(8:11).EQ.'Beam'.OR.dummy(8:11).EQ.'BEAM'.OR.
     .        dummy(8:11).EQ.'beam') THEN
          READ(dummy,*) cdum1,ebeam
        ELSE
          BACKSPACE 5
          ebeam = 99.0
        ENDIF

        WRITE(char(28),'(A,F8.3,A)') 'Ebeam =',ebeam,' eV'
        WRITE(char(29),'(A   )') 'UNITS =   m'
        WRITE(char(30),'(A,I4)') 'STEP   ',rel_step
        CALL LoadEIRENEAtomicData

      ELSEIF (iopt.EQ.6) THEN
c...    e + H2+ -> H + H+ + e (MAD):
        standard = .TRUE.
        vacuum   = .TRUE.

        READ(5,'(A80)') dummy
        IF   (dummy(8:11).EQ.'Beam'.OR.dummy(8:11).EQ.'BEAM'.OR.
     .        dummy(8:11).EQ.'beam') THEN
          READ(dummy,*) cdum1,ebeam
        ELSE
          BACKSPACE 5
          ebeam = 99.0
        ENDIF

        WRITE(char(28),'(A,F8.3,A)') 'Ebeam =',ebeam,' eV'
        WRITE(char(29),'(A   )') 'UNITS =   m'
        WRITE(char(30),'(A,I4)') 'STEP   ',rel_step
        CALL LoadEIRENEAtomicData

      ELSEIF (iopt.EQ.30) THEN
c...    D source from MAR:
        standard = .TRUE.
        vacuum   = .TRUE.
        WRITE(char(29),'(A   )') 'UNITS =   s m-3'
        WRITE(char(30),'(A,I4)') 'STEP   ',rel_step

        CALL LoadEIRENEAtomicData

      ELSEIF (iopt.EQ.40) THEN
c...    D2+ density on the standard grid:
        standard = .TRUE.
        vacuum   = .TRUE.
        WRITE(char(29),'(A   )') 'UNITS =   s'
        WRITE(char(30),'(A,I4)') 'STEP   ',rel_step

        CALL LoadEIRENEAtomicData

      ELSEIF (iopt.EQ.50) THEN
c...    plasma density everywhere:
        standard = .TRUE.
        vacuum   = .TRUE.
        WRITE(char(29),'(10X,A   )') 'UNITS =   m-3'
        WRITE(char(30),'(10X,A,I4)') 'STEP   ',rel_step

      ELSEIF (iopt.EQ.51) THEN
c...    D2+ density:
        standard = .TRUE.
        vacuum   = .TRUE.
        WRITE(char(29),'(A   )') 'UNITS =   m-3'
        WRITE(char(30),'(A,I4)') 'STEP   ',rel_step

        CALL LoadEIRENEAtomicData

      ELSEIF (iopt.EQ.52) THEN
c...    Recombination:
        standard = .TRUE.
        vacuum   = .FALSE.
        WRITE(char(29),'(A   )') 'UNITS =   m-3 s-1'
        WRITE(char(30),'(A,I4)') 'STEP   ',rel_step

      ELSEIF (iopt.EQ.53) THEN
c...    Ionisation:
        standard = .TRUE.
        vacuum   = .FALSE.
        WRITE(char(29),'(A   )') 'UNITS =   m-3 s-1'
        WRITE(char(30),'(A,I4)') 'STEP   ',rel_step

      ELSEIF (iopt.EQ.54) THEN
c...    Te:
        standard = .TRUE.
        vacuum   = .TRUE.
c        WRITE(char(29),'(A   )') 'UNITS =   eV'
        WRITE(char(30),'(A,I4)') 'STEP   ',rel_step

      ELSEIF (iopt.EQ.55) THEN
c...     D2+ density(?):

      ELSEIF (iopt.EQ.56) THEN
c...    Dalpha MFP (D density based):
        standard = .TRUE.
        vacuum   = .TRUE.

      ELSEIF (iopt.EQ.57) THEN
c...    Line radiation:
c       
c       Sub-option (LINE):
c         1 - Dalpha
c         2 - Dgamma     
c         3 - Dbeta
c
        READ(graph(14:15),*) line
        READ(graph(17:18),*) component
        WRITE(6,*) 'LINE:',line,component

        standard = .TRUE.
        vacuum   = .FALSE.

      ELSEIF (iopt.EQ.58) THEN
c...    Recombination + MAR (plot 30):
        standard = .TRUE.
        vacuum   = .TRUE.
        WRITE(char(29),'(A   )') 'UNITS =   m-3 s-1'
        WRITE(char(30),'(A,I4)') 'STEP   ',rel_step

      ELSEIF (iopt.EQ.59.OR.iopt.EQ.60) THEN
c...    Facelift: *TEMP*
        standard = .TRUE.
        vacuum   = .FALSE.
        WRITE(char(29),'(A   )') 'UNITS =  W m-3'
        WRITE(char(30),'(A,I4)') 'STEP   ',rel_step

      ELSEIF (iopt.GE.70.AND.iopt.LE.79) THEN
c...    Impurity density states 0 through 9:
        standard = .TRUE.
        vacuum   = .FALSE.
      ELSEIF (iopt.GE.80.AND.iopt.LE.89) THEN
c...    Impurity emission from ADAS:
        standard = .TRUE.
        vacuum   = .FALSE.
        READ(5,'(A1024)') dummy
        IF (dummy(8:11).EQ.'Adas'.OR.dummy(8:11).EQ.'ADAS'.OR.
     .      dummy(8:11).EQ.'adas') THEN
          BACKSPACE 5
          CALL RDG1 (GRAPH3,ADASID,adasyr,adasex,ISELE,ISELR,ISELX,
     .               ISELD,IERR)
          WRITE(0,*) 'ADAS SETTINGS:',adasid,adasyr,adasex
          WRITE(0,*) '             :',isele,iselr,iselx,iseld
          IF (IERR.NE.0) THEN
            WRITE(6,*) '984: ERROR READING ADAS DETAILS, IERR = ',IERR
            IERR = 0
            GOTO 99
          ENDIF
          izmin = iopt - 80
          WRITE(0,*) '             :',cion,izmin
          CALL LDADAS(cion,IZMIN,ADASID,ADASYR,ADASEX,ISELE,ISELR,ISELX,
     >                plrpad,Wlngth,IRCODE)
          WRITE(0,*) 'ADAS DATA:',izmin,wlngth,ircode
        ELSE
          CALL ER('984','ADAS data file not found',*99)
        ENDIF

c...    Impurity temperature to background ion temperature ratio:
c        standard = .TRUE.
c        vacuum   = .FALSE.
      ELSEIF (iopt.GE.90.AND.iopt.LE.99) THEN
c...    Impurity velocity to background ion parallel flow velocity ratio:
        standard = .TRUE.
        vacuum   = .FALSE.

      ELSEIF (iopt.GE.200) THEN
c DEFUNCT?
c...    Load external data array to be plotted:
        standard = .FALSE.
        vacuum   = .FALSE.


      ELSE
        CALL ER('982','Unknown plot option',*99)
      ENDIF

      XLAB   = '   R  (M)'
      YLAB   = '   Z  (M)'

      NVIEW  = '                                    '
      PLANE  = '                                    '
      SMOOTH = '                                    '//
     .         '                                    '
      ANLY   = '                                    '
      REF    = graph(14:LEN_TRIM(graph))



      nc = 0

c...  Check if data is to be plotted below a specified Z value only:
      zlim = -HI
c      zlim = HI
      READ(5,'(A80)',END=10) dummy
      IF   (dummy(8:11).EQ.'Zlim'.OR.dummy(8:11).EQ.'ZLIM'.OR.
     .      dummy(8:11).EQ.'zlim') THEN
        READ(dummy,*) cdum1,zlim
      ELSE
c        WRITE(0,*) 'dummy zlim:',TRIM(dummy)
        BACKSPACE 5
      ENDIF
10    CONTINUE
c...  Check if data is to be plotted below a specified R value only:
      rlim = 0.0
      READ(5,'(A80)',END=12) dummy
      IF   (dummy(8:11).EQ.'Rlim'.OR.dummy(8:11).EQ.'RLIM'.OR.
     .      dummy(8:11).EQ.'rlim') THEN
        READ(dummy,*) cdum1,rlim
      ELSE
c        WRITE(0,*) 'dummy rlim:',TRIM(dummy)
        BACKSPACE 5
      ENDIF
12    CONTINUE

      IF ((asc_3dmode.EQ.1.OR.asc_3dmode.EQ.2).AND.zval.EQ.-99.0) THEN
c...    Look for zval data in OUT input file:
        READ(5,'(A80)',END=15) dummy
        IF   (dummy(8:11).EQ.'Zval'.OR.dummy(8:11).EQ.'ZVAL'.OR.
     .        dummy(8:11).EQ.'zval') THEN
          READ(dummy,*) cdum1,zval
          IF (eirzaa.LT.0.0.AND.eirzaa.NE.-1.0) THEN
            zval = -zval / eirzaa * 360.0
            WRITE(char(29),'(A,F6.1,A)') 'Tval  = ',zval,' degrees'
          ELSE
            WRITE(char(29),'(A,F6.3,A)') 'Tval  = ',zval,' m'
          ENDIF
        ELSE
c          WRITE(0,*) 'dummy zval:',TRIM(dummy)
          BACKSPACE 5
        ENDIF
15      CONTINUE
      ENDIF

      xval = -99.0
      IF (asc_3dmode.EQ.2) THEN
c...    Look for xval data in OUT input file:
        READ(5,'(A80)',END=20) dummy
        IF   (dummy(8:11).EQ.'Xval'.OR.dummy(8:11).EQ.'XVAL'.OR.
     .        dummy(8:11).EQ.'xval') THEN
          READ(dummy,*) cdum1,xval
          WRITE(char(29),'(A,F6.3,A)') 'Xval  = ',xval,' m'
        ELSE
c          WRITE(0,*) 'dummy zval:',TRIM(dummy)
          BACKSPACE 5
        ENDIF
20      CONTINUE
      ENDIF

c...  Check for vacuum cells that are to be avoided:
25    READ(5,'(A80)',END=30) dummy
      IF   (dummy(8:11).EQ.'Kill'.OR.dummy(8:11).EQ.'KILL'.OR.
     .      dummy(8:11).EQ.'kill') THEN
        READ(dummy,*) cdum1,idum1
        avoid(idum1) = 1
        GOTO 25
      ELSE
c        WRITE(0,*) 'dummy kill:',TRIM(dummy)
        BACKSPACE 5
      ENDIF
30    CONTINUE
      IF (iopt.EQ.79) THEN
35      READ(5,'(A80)',END=40) dummy
        IF (dummy(8:9).EQ.'Iz'.OR.dummy(8:9).EQ.'IZ'.OR.
     .      dummy(8:9).EQ.'iz') THEN
          READ(dummy,*) cdum1,iz
        ELSE
          CALL ER('Plot984','79: ionisations state not '//
     .            'specified in the OUT input file',*99) 
        ENDIF
      ELSE
        iz = iopt-70
      ENDIF
 40   CONTINUE


      IF (standard) THEN
c...    Assign polygons from standard grid:
        IF (eirnsdtor.GT.1.AND.xval.NE.-99.0) THEN
          istart = 1
          iend   = eirnsdtor
        ELSE
          istart = 1
          iend   = 1
        ENDIF

        shift = 0

        IF (eirnsdtor.GT.1.AND.zval.NE.-99.0) THEN
          DO i1 = 1, eirnsdtor-1
            IF (zval.GE.eirsdtor(i1)) THEN
              shift = MAXBGK * (i1 - 1) 
              WRITE(6,*) '982: STD SHIFT BOOST ',zval,i1,shift
c              WRITE(0,*) '982: STD SHIFT BOOST ',zval,i1,shift,
c     .          eirsdtor(i1)
            ENDIF
          ENDDO
        ENDIF

        DO cut1 = istart, iend

          DO ir = 2, nrs
            DO ik = 1, nks(ir)

              id = korpg(ik,ir)
              status = .FALSE.

              IF (id.EQ.0.OR.nvertp(id).EQ.0) CYCLE

              IF (zs(ik,ir).LT.zlim) CYCLE
c              IF (zs(ik,ir).GT.zlim) CYCLE
              IF (rs(ik,ir).LT.rlim) CYCLE

              IF (xval.NE.-99.0) THEN
                xmin =  HI
                ymin =  HI
                xmax = -HI
                ymax = -HI
                DO i1 = 1, nvertp(id)
                  xmin = MIN(xmin,rvertp(i1,id))
                  ymin = MIN(ymin,zvertp(i1,id))
                  xmax = MAX(xmax,rvertp(i1,id))
                  ymax = MAX(ymax,zvertp(i1,id))
                ENDDO

c...            Decide if cell is in the plotting region:
                IF (xval.GE.xmin.AND.xval.LT.xmax) status = .TRUE.

                zmin = eirsdtor(cut1)
                zmax = eirsdtor(MIN(cut1+1,iend))
                IF (cut1.EQ.iend) zmax = eirzaa

                shift = MAXBGK * (cut1 - 1) 
              ELSE
                DO i1 = 1, nvertp(id)
                  x1 = rvertp(i1,id)
                  y1 = zvertp(i1,id)
c...              Decide if cell is in the plotting region:
                  IF (x1.GE.xxmin.AND.x1.LE.xxmax.AND.
     .                y1.GE.yymin.AND.y1.LE.yymax) status = .TRUE.
                ENDDO
              ENDIF

              IF (status) THEN
c...            If the cell is in the plotting area, then add
c               the appropriate data to the plot arrays:

                IF (xval.NE.-99.0) THEN
c...              Add cell geometry data:
                  nc = nc + 1
                  nv(nc) = 4
                  rv(1,nc) = zmax
                  zv(1,nc) = ymin
                  rv(2,nc) = zmin
                  zv(2,nc) = ymin
                  rv(3,nc) = zmin
                  zv(3,nc) = ymax
                  rv(4,nc) = zmax
                  zv(4,nc) = ymax
                ELSE
                  nc = nc + 1
                  nv(nc) = nvertp(id)
                  DO i1 = 1, nvertp(id)
                    rv(i1,nc) = rvertp(i1,id)
                    zv(i1,nc) = zvertp(i1,id)
                  ENDDO
                ENDIF

                ti = ktibs(ik,ir)
                ne = knbs (ik,ir)

                nD  = pinbgk(ik,ir,10+1+shift) * 1.0E+6
                nD2 = pinbgk(ik,ir,15+1+shift) * 1.0E+6

c...            Assign quantity to be plotted:
                IF     (iopt.EQ.1) THEN
c...              p + D2 -> p + D2
                  sigmavb = GetEAD(MAX(ti,tD2),0.0,12,'H.3 ')
                  IF (vbeam.EQ.99.0) THEN
                    tD2  = pinbgk(ik,ir,15+2+shift)
                    vavg = SQRT((3.0 * tD2 * ECH) / (2 * crmb * AMU))
                    cq(nc) = vavg / (ne * sigmavb * 1.0E-06)                 
                  ELSE
                    cq(nc) = vbeam / (ne * sigmavb * 1.0E-06) 
                  ENDIF
                ELSEIF (iopt.EQ.2) THEN
c...              Some sort of e + H2 reaction:

c                  IF (ir.LT.irsep) CYCLE
c                  IF (ir.LT.irwall.AND.ik.GT.nks(ir)/2) CYCLE

c...              e + H2 + e -> e + H2+
                  sigmav =          GetEAD(ti,ne,19,'H.4 ')
c...              e + H2 -> e + H + H+
                  sigmav = sigmav + GetEAD(ti,ne,13,'H.4 ')
c...              e + H2 -> e + H + H+
                  sigmav = sigmav + GetEAD(ti,ne,14,'H.4 ')

                  tD2  = pinbgk(ik,ir,15+2+shift)

                  IF (ebeam.NE.99.0) tD2 = ebeam

                  vavg = SQRT((3.0 * tD2 * ECH) / (2 * crmb * AMU))
                  cq(nc) = vavg / (ne * sigmav * 1.0E-06) 
                ELSEIF (iopt.EQ.3) THEN
c...              e + H -> 2e + H+:
                  sigmav = GetEAD(ti,ne,1,'H.4 ')

                  tD = pinbgk(ik,ir,10+2+shift)
                  vavg = SQRT((3.0 * tD * ECH) / (crmb * AMU))
                  cq(nc) = vavg / (ne * sigmav * 1.0E-06) 
                ELSEIF (iopt.EQ.4) THEN
c...              p + H2(v) -> H + H2+
                  sigmavb = GetEAD(MAX(ti,tD2),0.0,18,'H.3 ')
                  IF (vbeam.EQ.99.0) THEN
                    tD2  = pinbgk(ik,ir,15+2+shift)
                    vavg = SQRT((3.0 * tD2 * ECH) / (2 * crmb * AMU))
                    cq(nc) = vavg / (ne * sigmavb * 1.0E-06)                 
                  ELSE
                    cq(nc) = vbeam / (ne * sigmavb * 1.0E-06) 
                  ENDIF

                ELSEIF (iopt.EQ.5) THEN
c...              MAR:
c...              e + H2+ -> e + H + H
                  sigmav = GetEAD(ti,ne,17,'H.4 ')

                  tD2  = pinbgk(ik,ir,15+2+shift)

                  IF (ebeam.NE.99.0) tD2 = ebeam

                  vavg = SQRT((3.0 * tD2 * ECH) / (2 * crmb * AMU))
                  cq(nc) = vavg / (ne * sigmav * 1.0E-06) 

                ELSEIF (iopt.EQ.6) THEN
c...              MAD:
c...              e + H2+ -> e + H + H+
                  sigmav = GetEAD(ti,ne,15,'H.4 ')

                  tD2  = pinbgk(ik,ir,15+2+shift)

                  IF (ebeam.NE.99.0) tD2 = ebeam

                  vavg = SQRT((3.0 * tD2 * ECH) / (2 * crmb * AMU))
                  cq(nc) = vavg / (ne * sigmav * 1.0E-06) 

                ELSEIF (iopt.EQ.30) THEN
c...              D source rate from MAR:

                  IF (ne.LT.1.0E+15) CYCLE

c...              D2+ density:
                  sigmav = GetEAD(ti,ne,22,'H.11')
c                  sigmav = GetEAD(ti,ne,22,'H.12')
                  nD2p = sigmav * pinmol(ik,ir)

c                  WRITE(0,*) 'DATA:',ti,sigmav

c...              MAR: e + H2+ -> e + H + H
                  sigmav1 = GetEAD(ti,ne,17,'H.4 ')
c...              MAD: e + H2+ -> e + H + H+
                  sigmav2 = GetEAD(ti,ne,15,'H.4 ')

                  sigmav = sigmav1**2 / (sigmav1 + sigmav2)

                  cq(nc) = sigmav * 1.0E-06 * ne * nD2p * 2.0

                  MARsum = MARsum + cq(nc) * kvols(ik,ir)

                ELSEIF (iopt.EQ.40) THEN
c...              D2+ lifetime:

c...              MAR:
c...              e + H2+ -> e + H + H
                  sigmav =          GetEAD(ti,ne,17,'H.4 ')
c...              MAD:
c...              e + H2+ -> e + H + H+
                  sigmav = sigmav + GetEAD(ti,ne,15,'H.4 ')

                  IF (ne.LT.1.0E+15) CYCLE

                  cq(nc) = 1.0 / (ne * sigmav * 1.0E-06) 


                ELSEIF (iopt.EQ.50) THEN
c...              plasma density everywhere:
                  cq(nc) = ne

                ELSEIF (iopt.EQ.51) THEN
c...              D2+ density on standard grid:
c                  cq(nc) = pinmoi(ik,ir)

                  sigmav = GetEAD(ti,ne,22,'H.4 ')

                  cq(nc) = sigmav * pinmol(ik,ir)

                  WRITE(0,*) 'SIGMAV 1:',sigmav

                ELSEIF (iopt.EQ.52) THEN
c...              Recombination:
                  cq(nc) = pinrec(ik,ir)

                ELSEIF (iopt.EQ.53) THEN
c...              Ionisation:
                  cq(nc) = pinion(ik,ir)

                ELSEIF (iopt.EQ.54) THEN
c...              Ti:
c                  IF (ir.LT.irsep) CYCLE

                  cq(nc) = ti


                ELSEIF (iopt.EQ.55) THEN

                  sigmav = GetEAD(ti,ne,22,'H.4 ')

                  cq(nc) = sigmav * pinmol(ik,ir) / ne

                ELSEIF (iopt.EQ.56) THEN
c...              Dalpha MFP (D density based):

                  IF (Ti.EQ.0.0.OR.nD.LT.1.0E+16) CYCLE

c                  IF (.NOT.(ir.LT.irwall.AND.ik.GT.nks(ir)/2)) CYCLE


                  cq(nc) = 2.2E-4 * SQRT(Ti / 1.0) / (nD + nD2)*1.0E+20

c                  IF (ir.GT.19.AND.ik.GT.nks(ir)-10) 
c     .               WRITE(0,*) '--->',ik,ir,cq(nc)

                ELSEIF (iopt.EQ.57) THEN
c...              Line radiation:
                  cq(nc) = pinline(ik,ir,component,line)

                ELSEIF (iopt.EQ.58) THEN
c...              Recombination + MAR:
                  IF (ne.LT.1.0E+15) CYCLE

c...              D2+ equilibrium density:
                  sigmav = GetEAD(ti,ne,22,'H.11')
                  nD2p = sigmav * pinmol(ik,ir)

c...              MAR: e + H2+ -> e + H + H
                  sigmav1 = GetEAD(ti,ne,17,'H.4 ')
c...              MAD: e + H2+ -> e + H + H+
                  sigmav2 = GetEAD(ti,ne,15,'H.4 ')

                  sigmav = sigmav1**2 / (sigmav1 + sigmav2)

                  cq(nc) = sigmav * 1.0E-06 * ne * nD2p * 2.0

                  MARsum = MARsum + cq(nc) * kvols(ik,ir)

c...              Direct recombination:
                  cq(nc) = cq(nc) + pinrec(ik,ir)

                ELSEIF (iopt.EQ.59) THEN
c...              Line radiation:
c                  IF (eirpho2(ik,ir).GT.0.0) THEN
c                    cq(nc) = eirpho2(ik,ir)
                  IF (eirpho1(ik,ir).GT.0.0) THEN
                    cq(nc) = eirpho1(ik,ir)
                  ELSE
                    cq(nc) = 0.0
                  ENDIF

                ELSEIF (iopt.EQ.60) THEN
c...              Line radiation:
                  IF (eirpho2(ik,ir).LT.0.0) THEN
                    cq(nc) = -eirpho2(ik,ir)
                  ELSE
                    cq(nc) = 0.0
                  ENDIF

                ELSEIF (iopt.GE.70.AND.iopt.LE.79) THEN
c...              Impurity density for charge states 0 through 9:
                  IF (iz.EQ.-1) THEN  ! Add up all the ionisation states and plot
                    DO iz = 1, cion
                      cq(nc) = cq(nc) + MAX(0.0,sdlims(ik,ir,iz))
                    ENDDO
                  ELSE
                    cq(nc) = MAX(0.0,sdlims(ik,ir,iz))
c                    IF (iz.EQ.0.AND.ik.EQ.1.AND.ir.EQ.109) 
c     .                WRITE(0,*) sdlims(1:nks(109),109,iz)
                  ENDIF

                ELSEIF (iopt.GE.80.AND.iopt.LE.89) THEN
c...              Impurity emission from ADAS:
                  cq(nc) = plrpad(ik,ir)
c...              Impurity temperature ratio for states 0 through 9:
c                  cq(nc) = sdts(ik,ir,iopt-80) / (ktebs(ik,ir) + LO)

                ELSEIF (iopt.GE.90.AND.iopt.LE.99) THEN
c...              Impurity velocity ratio for states 0 through 9:
                  cq(nc) = velavg(ik,ir,iopt-90) / 
     .                     ((kvhs(ik,ir) + 1.0E-10) / qtim)

                ELSE

                  CALL ER('982','Unknown quantity to be plotted',*99)
                ENDIF
              ENDIF
            ENDDO    
          ENDDO

        ENDDO

      ENDIF

c...  Decide whether to load most recent iteration data (set=1) or average over
c     last 5 (currently) iterations (set=2):
      set = 2

      IF (vacuum) THEN
c...    Assign polygons from vacuum grid:

        shift = 1 + eirnpgdat

        IF (asc_3dmode.EQ.2.AND.zval.NE.-99.0) THEN
          DO i1 = 1, ascncut
            IF (zval.GE.asc_zmin3D(i1).AND.zval.LT.asc_zmax3D(i1)) THEN
              shift = shift + asc_ncell * (i1-1)
              WRITE(6,*) '982: VAC SHIFT BOOST ',i1-1
              EXIT
            ENDIF
          ENDDO
        ENDIF

        IF (xval.NE.-99.0) THEN
          istart = 1
          iend   = asc_ncell * ascncut
        ELSE
          istart = 1
          iend   = asc_ncell
        ENDIF

        DO i1 = istart, iend

c...      See if cell is black listed:
          IF (avoid(i1).EQ.1) CYCLE

c...      Find x,ymin and x,ymax for the cell:
          xmin =  HI
          ymin =  HI
          xmax = -HI
          ymax = -HI
          i3 = MOD(i1,asc_ncell)            
          IF (i3.EQ.0) i3 = asc_ncell
          DO i2 = 1, ascnvertex(i3)
            xmin = MIN(xmin,ascvertex(i2*2-1,i3))
            ymin = MIN(ymin,ascvertex(i2*2  ,i3))
            xmax = MAX(xmax,ascvertex(i2*2-1,i3))
            ymax = MAX(ymax,ascvertex(i2*2  ,i3))
          ENDDO

c...      If the cell is above the specified Z limit (DIVIMP
c         coordinate system now), then skip the cell:
c          IF (0.5*(ymin+ymax).LT.zlim) CYCLE
c          IF (0.5*(ymin+ymax).GT.zlim) CYCLE
c          IF (0.5*(xmin+xmax).LT.rlim) CYCLE

          IF (xval.NE.-99.0) THEN
c...        Toroidal plot specified:
            IF (.NOT.(xmin.LE.xval.AND.xmax.GT.xval)) CYCLE
c...        Find zmin and zmax for the cell:
            i3 = INT(REAL(i1) / (REAL(asc_ncell) + 0.5)) + 1
            zmin = asc_zmin3D(i3)
            zmax = asc_zmax3D(i3)
            nc = nc + 1
            nv(nc) = 4
            rv(1,nc) = zmax
            zv(1,nc) = ymin
            rv(2,nc) = zmin
            zv(2,nc) = ymin
            rv(3,nc) = zmin
            zv(3,nc) = ymax
            rv(4,nc) = zmax
            zv(4,nc) = ymax
          ELSE
            nc = nc + 1
            nv(nc) = ascnvertex(i1)
            DO i2 = 1, ascnvertex(i1)
              rv(i2,nc) = ascvertex(i2*2-1,i1)
              zv(i2,nc) = ascvertex(i2*2  ,i1)
            ENDDO
          ENDIF

c...      PINASD(index,quantity,species,relaxation):        

          ti = pinasd(i1+shift,2,5,1)
          ne = pinasd(i1+shift,1,5,1) * 1.0E+06

          nD  = pinasd(i1+shift,1,1,set) * 1.0E+06
          nD2 = pinasd(i1+shift,1,3,set) * 1.0E+06

          IF (ti.LE.0.1.OR.ne.LE.1.0E+12) CYCLE

          IF     (iopt.EQ.1) THEN
            sigmavb = GetEAD(MAX(ti,tD2),0.0,12,'H.3 ') * 1.0E-06
            IF (vbeam.EQ.99.0) THEN
              tD2  = pinasd(i1+shift,2,8,set)
              vavg = SQRT((3.0 * tD2 * ECH) / (2 * crmb * AMU))
              cq(nc) = vavg / (ne * sigmavb) 
            ELSE
              cq(nc) = vbeam / (ne * sigmavb) 
            ENDIF

          ELSEIF (iopt.EQ.2) THEN
c...        Some sort of e + H2 reaction:
            sigmav =          GetEAD(ti,ne,19,'H.3 ')
            sigmav = sigmav + GetEAD(ti,ne,13,'H.3 ')
            sigmav = sigmav + GetEAD(ti,ne,14,'H.3 ')
            tD2  = pinasd(i1+shift,2,8,set)
            IF (ebeam.NE.99.0) tD2 = ebeam
            vavg = SQRT((3.0 * tD2 * ECH) / (2 * crmb * AMU))
            cq(nc) = vavg / (ne * sigmav * 1.0E-06) 

          ELSEIF (iopt.EQ.3) THEN
c...        e + H -> 2E + H+:
            sigmav = GetEAD(ti,ne,1,'H.4 ')
            tD = pinasd(i1+shift,2,7,set)
            vavg = SQRT((3.0 * tD * ECH) / (crmb * AMU))
            cq(nc) = vavg / (ne * sigmav * 1.0E-06) 

          ELSEIF (iopt.EQ.4) THEN
c...        p + H2(v) -> H + H2+
            sigmavb = GetEAD(MAX(ti,tD2),0.0,18,'H.3 ') * 1.0E-06
            IF (vbeam.EQ.99.0) THEN
              tD2  = pinasd(i1+shift,2,8,set)
              vavg = SQRT((3.0 * tD2 * ECH) / (2 * crmb * AMU))
              cq(nc) = vavg / (ne * sigmavb) 
            ELSE
              cq(nc) = vbeam / (ne * sigmavb) 
            ENDIF

          ELSEIF (iopt.EQ.5) THEN
c...        MAR:
c...        e + H2+ -> e + H + H
            sigmav = GetEAD(ti,ne,17,'H.3 ')
            tD2  = pinasd(i1+shift,2,8,set)
            IF (ebeam.NE.99.0) tD2 = ebeam
            vavg = SQRT((3.0 * tD2 * ECH) / (2 * crmb * AMU))
            cq(nc) = vavg / (ne * sigmav * 1.0E-06) 

          ELSEIF (iopt.EQ.6) THEN
c...        MAD:
c...        e + H2+ -> e + H + H+
            sigmav = GetEAD(ti,ne,15,'H.3 ')
            tD2  = pinasd(i1+shift,2,8,set)
            IF (ebeam.NE.99.0) tD2 = ebeam
            vavg = SQRT((3.0 * tD2 * ECH) / (2 * crmb * AMU))
            cq(nc) = vavg / (ne * sigmav * 1.0E-06) 

          ELSEIF (iopt.EQ.30) THEN
c...        D source rate from MAR:
            IF (ne.LT.1.0E+15) CYCLE
	  
c...        D2+ density:
            sigmav = GetEAD(ti,ne,22,'H.3 ')
            nD2p = sigmav * pinasd(i1+shift,1,3,set) * 1.0E+06
	  
c...        MAR: e + H2+ -> e + H + H
            sigmav1 = GetEAD(ti,ne,17,'H.4 ')
c...        MAD: e + H2+ -> e + H + H+
            sigmav2 = GetEAD(ti,ne,15,'H.4 ')
            sigmav = sigmav1**2 / (sigmav1 + sigmav2)
            cq(nc) = sigmav * 1.0E-06 * ne * nD2p * 2.0
c...        * Estimated! *
            cellvol = (xmax - xmin) * (ymax - ymin) *
     .                2.0 * PI * 0.5 * (xmin + xmax) 
            WRITE(0,*) 'estimated cellvol:',cellvol
            MARsum = MARsum + cq(nc) * cellvol

          ELSEIF (iopt.EQ.40) THEN
c...        D2+ lifetime:
c...        MAR: e + H2+ -> e + H + H
            sigmav =          GetEAD(ti,ne,17,'H.4 ')
c...        MAD: e + H2+ -> e + H + H+
            sigmav = sigmav + GetEAD(ti,ne,15,'H.4 ')
            IF (ne.LT.1.0E+15) CYCLE
            cq(nc) = 1.0 / (ne * sigmav * 1.0E-06) 
            WRITE(0,*) 'SIGMA:',sigmav,ti,ne,cq(nc)

          ELSEIF (iopt.EQ.50) THEN
c...        plasma density everywhere:
            cq(nc) = ne

          ELSEIF (iopt.EQ.51) THEN
            sigmav = GetEAD(ti,ne,22,'H.4 ')
c            nD2 = pinasd(i1+shift,1,3,set) * 1.0E+06
            cq(nc) = sigmav * nD2
            WRITE(0,*) 'SIGMAV 2:',sigmav,nD2

          ELSEIF (iopt.EQ.54) THEN
c...        Ti:
            cq(nc) = ti

          ELSEIF (iopt.EQ.55) THEN
c...        D2+ density(?):
            sigmav = GetEAD(ti,ne,22,'H.4 ')
            nD2 = pinasd(i1+shift,1,3,set) * 1.0E+06
            cq(nc) = sigmav * nD2 / ne

          ELSEIF (iopt.EQ.56) THEN
c...        Dalpha MFP (D density based):
            IF (Ti.EQ.0.0.OR.nD.LT.1.0E+16) CYCLE
            CYCLE
            cq(nc) = 2.2E-4 * SQRT(Ti / 1.0) / nD * 1.0E+20

          ELSEIF (iopt.EQ.58) THEN
c...        Recombination + MAR:
            IF (ne.LT.1.0E+15) CYCLE

c...        D2+ equilibrium density:
            sigmav = GetEAD(ti,ne,22,'H.11')
            nD2p = sigmav * pinasd(i1+shift,1,3,set) * 1.0E+06

c...        MAR: e + H2+ -> e + H + H
            sigmav1 = GetEAD(ti,ne,17,'H.4 ')
c...        MAD: e + H2+ -> e + H + H+
            sigmav2 = GetEAD(ti,ne,15,'H.4 ')

            sigmav = sigmav1**2 / (sigmav1 + sigmav2)
            cq(nc) = sigmav * 1.0E-06 * ne * nD2p * 2.0

c...        * Estimated! *
            cellvol = (xmax - xmin) * (ymax - ymin) *
     .                2.0 * PI * 0.5 * (xmin + xmax) 
            WRITE(0,*) 'estimated cellvol:',cellvol
            MARsum = MARsum + cq(nc) * cellvol

c...         Direct recombination:
c             cq(nc) = cq(nc) + pinrec(ik,ir)


          ELSE
            CALL ER('982','Unknown quantity to be plotted on '//
     .                    'vacuum grid',*99)
          ENDIF

c...      Store the additional cell region of the cell for averaging limitations
c         in the contour plot:
          rq(nc) = asc_region(i1)


c...      Output data to .OUT file:
          WRITE(6,'(A,3I6,1P,E12.4,0P)') 
     .      '984 STD:',i1,i1+shift,nc,cq(nc)

        ENDDO
      ENDIF


c...  Graph outline:

c...  Get the size of the plot, if non-standard:
      READ(5,'(A256)') dummy
      IF   (dummy(8:11).EQ.'Size'.OR.dummy(8:11).EQ.'size'.OR.
     .      dummy(8:11).EQ.'SIZE') THEN
        READ(dummy,*) cdum1,map1x,map2x,map1y,map2y
      ELSE
        BACKSPACE 5
      ENDIF


      CALL CustomisePlot(title,xlab,ylab,elabs)


c...  For some reason the case name is not being passed from DIVIMP through GET?   
      CALL GetEnv('CASENAME',job(1:36))

      slopt4 = 1
      char(30) = ' '
      REF    = '                                    '
c      iopt_ghost = 1
      CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     .             YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NRS)


c...  Find plotting quantity minimum and maximum:
      nconts = 10
      qmax   = -HI
      qmin   =  HI
      DO i1 = 1, nc
c...    Decide if the cell is within the viewing range:
        inside = .FALSE.
        IF (cq(i1).EQ.0.0) CYCLE
        DO i2 = 1, nv(i1)
          IF (rv(i2,i1).GE.cxmin.AND.rv(i2,i1).LE.cxmax.AND.
     .        zv(i2,i1).GE.cymin.AND.zv(i2,i1).LE.cymax) inside = .TRUE.
        ENDDO
        IF (inside) THEN
          IF (cq(i1).NE.0.0.AND.qmin.GT.cq(i1)) qmin = cq(i1)
          IF (cq(i1).NE.0.0.AND.qmax.LT.cq(i1)) qmax = cq(i1)
        ENDIF
      ENDDO

c...  Decide the number of decades:
      IF     (iopt.GE.50.AND.iopt.LE.199) THEN

c * THIS REALLY BLOWS.  NEED TO MERGE PLOTS 982 and 984

        lmin = REAL(INT(LOG10(qmin)) - 1)
        lmax = REAL(INT(LOG10(qmax)) + 1)
 
        IF (iopt.EQ.50) THEN
          IF (grayscale) THEN
            lmin = MAX(lmin,19.0)
            lmax = MIN(lmax,21.0)
          ELSE
            lmin = MAX(lmin,18.0)
            lmax = MIN(lmax,21.0)
          ENDIF
        ENDIF

        IF (iopt.EQ.51) THEN
          lmax = 18.0
        ENDIF

        IF (iopt.EQ.52.OR.iopt.EQ.58) THEN
          lmin = 23
        ENDIF

        IF (iopt.EQ.53) THEN
          lmin = lmax - 1
        ENDIF

        IF (.not..false..and.iopt.EQ.54) THEN
          lmax = lmax - 1
        ENDIF

        IF (iopt.EQ.55) THEN
          lmin = -5
          lmax =  0
        ENDIF

        IF (iopt.EQ.56) THEN
          lmin = -5
          lmax =  1
        ENDIF


        WRITE(6,*) 'LMIN:',lmin,lmax,qmin,qmax



      ELSEIF (iopt.GE.40) THEN

        lmin = REAL(INT(LOG10(qmin)) - 1)
        lmax = REAL(INT(LOG10(qmax)) + 1)

        lmax = -5.0


      ELSEIF (iopt.GE.30) THEN

        lmin = REAL(INT(LOG10(qmin)) - 1)
        lmax = REAL(INT(LOG10(qmax)) + 1)

        lmin = 22.0
        lmax = 24.0


      ELSE
        lmin = MAX(REAL(INT(LOG10(qmin)) - 1),-4.0)
        lmax = MIN(REAL(INT(LOG10(qmax)) + 1), 0.0)

        lmin = -4.0
        lmax =  0.0

      ENDIF

c...  Check scale adjustment:
      READ(5,'(A80)') dummy
      IF   (dummy(8:11).EQ.'Scal'.OR.dummy(8:11).EQ.'SCAL'.OR.
     .      dummy(8:11).EQ.'scal') THEN
        READ (dummy,*) cdum1,lmin,lmax
        qmin = 10**lmin
        qmax = 10**lmax
      ELSE
        BACKSPACE 5
      ENDIF


c...  Check if a contour plot is to be plotted:
      READ(5,'(A80)') dummy
      IF   (dummy(8:11).EQ.'Cont'.OR.dummy(8:11).EQ.'CONT'.OR.
     .      dummy(8:11).EQ.'cont') THEN
        BACKSPACE 5
        GOTO 100
      ELSE
        BACKSPACE 5
      ENDIF





c...  Draw polygons:
      CALL PSPACE(MAP1X,MAP2X,MAP1Y,MAP2Y)
      CALL MAP   (CXMIN,CXMAX,CYMIN,CYMAX)
      CALL HSV
      DO i1 = 1, nc
        IF (cq(i1).GE.qmin) THEN
          frac = (LOG10(cq(i1)) - lmin) / (lmax - lmin)
          frac5 = 100.0 * 0.98 * frac
          fmod5 = AMOD(frac5,1.0)
          frac = (frac5 - fmod5) / 100.0
c...dev
          bright = 1.0-(1.0-frac)**10        
          frac = (1.0 - frac) * 0.90
          frac = frac + 0.34
          IF (frac.GT.1.0) frac = frac - 1.0

c...new
          val = cq(i1)
c         val = 1.5
          decade = REAL(INT(LOG10(val)))
          IF (val.LT.1.0) decade = decade - 1.0
          decval = 10.0**decade + LO
          val = (REAL(INT(val / decval)) + 0.5) * decval
          decval = 10.0**decade + LO
          IF (grayscale) THEN
            bright = 0.0
c...         Colour:
c            frac = 0.66 * (LOG10(val) - lmin) / (lmax - lmin) + 0.34
c            frac = 1.0 - frac + 0.34
            frac = (LOG10(val) - lmin) / (lmax - lmin)
            frac = 0.9 * frac + 0.1
          ELSEIF (lmax-lmin.LE.2) THEN
            bright = 1.0
            frac = REAL(decade - lmin) / REAL(lmax - 1 - lmin) * 0.40 
            frac = frac + (val - decval) / (8.5 * decval) * 0.40
          ELSE
            bright = (val - decval) / (9.0 * decval) * 0.5 + 0.5
c            bright = (val - decval) / (8.5 * decval) * 0.3 + 0.7
            frac = REAL(decade - lmin) / REAL(lmax - 1 - lmin) * 0.65
c            WRITE(0,*) 'FRAC:',bright,frac
          ENDIF
 
c           WRITE(0,*) '???:',cq(i1),val,decval,decade
c           STOP

          IF (grayscale) THEN
            CALL ColSet(0.0,0.0,1.0-frac,255)
          ELSEIF (frac.LT.0.0) THEN
            CALL RGB
            CALL ColSet(0.0,0.0,0.0,255)
            CALL HSV
c            STOP '1'
          ELSEIF (frac.GT.0.65) THEN
            CALL RGB
            CALL ColSet(1.0,0.0,0.0,255)
            CALL HSV
c           WRITE(0,*) '???:',cq(i1),val,decval,decade
c           WRITE(0,*) '???:',frac,lmin,lmax
c           STOP '2'
          ELSE
            CALL ColSet(frac,1.0,bright,255)
          ENDIF

          WRITE(6,'(A,I6,3F10.4,1P,E12.4)') 
     .      '984:',i1,rv(1,i1),zv(1,i1),frac,cq(i1)   
 
          CALL FILCOL(255)
          CALL LINCOL(255)
          CALL PTPLOT(rv(1,i1),zv(1,i1),1,nv(i1),1)
        ENDIF
      ENDDO

c...  Finish off the plot boarder:
      CALL FULL
      CALL LINCOL(defcol)
      CALL PSPACE(MAP1X,MAP2X,MAP1Y,MAP2Y)
      CALL MAP(0.0,1.0,0.0,1.0)
      CALL POSITN (0.0,1.0)
      CALL JOIN   (1.0,1.0)
      CALL JOIN   (1.0,0.0)

      GOTO 200

100   CONTINUE


         WRITE(0,*) 'contouring',iopt_ghost
c
c        Read in contour options if contour plot was specified  
c
         call rdg_contopts(graph,icntr,ncntr,uconts,maxpts,
     >                     xcen,ycen,xnear,ynear,ierr)

      

        colscale = .TRUE.
 
        xres = 99
        yres = 99
        READ(5,'(A80)') dummy
        IF   (dummy(8:11).EQ.'Cres'.OR.dummy(8:11).EQ.'CRES'.OR.
     .        dummy(8:11).EQ.'cres') THEN
          READ(dummy,*) cdum1,xres,yres
        ELSE
          BACKSPACE 5
        ENDIF


c
c        Plot the image as a contour plot
c
         maxix = xres
         maxiy = yres
         nix   = xres
         niy   = yres       

         if (ierr.ne.0) 
     .     CALL ER('982','Problem loading contour plot options',*99)

         IF (ncntr.EQ.1.AND.uconts(1).EQ.-1.0) THEN
           ncntr = 15
c           ncntr = 20
           ncntr = 30
c           ncntr = 19
           DO i1 = 1, ncntr
             frac = REAL(i1-1) / REAL(ncntr)
             uconts(i1) = lmin + frac * (lmax - lmin) + (1 - lmin) 
             IF (i1.GT.1) uconts(i1) = uconts(i1) - 0.0001
c             uconts(i1) = 10**(lmin + frac * (lmax - lmin))
             WRITE(0,*) 'uconts:',uconts(i1),lmin,lmax
             WRITE(6,*) 'uconts:',uconts(i1),lmin,lmax
           ENDDO
         ENDIF

c...     Take the logarithm of the data:
         DO i1 = 1, nc
           IF (cq(i1).EQ.0.0) THEN
             cq(i1) = LO
           ELSE
             cq(i1) = LOG10(cq(i1)) + (1 - lmin)
           ENDIF
         ENDDO

c         ncntr = ncntr - 1

c
c        Set up colours based on contour options  
c
         if (icntr.eq.0.or.icntr.eq.2) then
             call setup_col(ncntr+1,2)
         elseif (icntr.eq.1) then 
             ncntr = 10
             call setup_col(ncntr+1,2)
         elseif (icntr.eq.3.or.icntr.eq.4) then
             ncntr = ncntr
             call setup_col(ncntr+1,2)
         endif

         PLANE = '                                    '
         smooth = ' '
         anly   = ' '
c
c        Set min and max levels for contour plot
c
         maxscale = qmax
         minscale = qmin
c
c
c        Set up axes for the contour plot call - use the 
c        base camera definition space at the given distance.
c
c
c        Allocate space for axes and check to see if space available
c
         allocate (image(xres,yres),STAT=ierr)
         if (ierr.ne.0) then 
            write(0,*) 'IMAGE array could not be allocated: ',xres
            write(6,*) 'IMAGE array could not be allocated: ',xres
            goto 99 
         endif        
         allocate (image1(xres,yres),STAT=ierr)
         if (ierr.ne.0) then 
            write(0,*) 'IMAGE1 array could not be allocated: ',xres
            write(6,*) 'IMAGE1 array could not be allocated: ',xres
            goto 99 
         endif        
         allocate (tpt(xres,yres),STAT=ierr)
         if (ierr.ne.0) then 
            write(0,*) 'TPT array could not be allocated: ',xres
            write(6,*) 'TPT array could not be allocated: ',xres
            goto 99 
         endif        
         allocate (tpr(xres,yres),STAT=ierr)
         if (ierr.ne.0) then 
            write(0,*) 'TPR array could not be allocated: ',xres
            write(6,*) 'TPR array could not be allocated: ',xres
            goto 99 
         endif        
         allocate (raxis(xres),STAT=ierr)
         if (ierr.ne.0) then 
            write(0,*) 'RAXIS array could not be allocated: ',xres
            write(6,*) 'RAXIS array could not be allocated: ',xres
            goto 99 
         endif        
         allocate (zaxis(yres),STAT=ierr)
         if (ierr.ne.0) then 
            write(0,*) 'ZAXIS array could not be allocated: ',yres
            write(6,*) 'ZAXIS array could not be allocated: ',yres
            goto 99 
         endif        

         CALL IZero(tpt,xres*yres)
         CALL IZero(tpr,xres*yres)
c
c        Assign axis values
c
         WRITE(6,*) 'XXMIN,MAX:',xxmin,xxmax
         WRITE(6,*) 'YYMIN,MAX:',yymin,yymax

         dr = (xxmax - xxmin) / xres
c
c        Column center coordinates
c
         do ic = 1, xres
           raxis(ic) = xxmin + REAL(ic - 1) * dr + 0.5 * dr
            write(6,'(a,i5,2(1x,g12.5))') 'RAXIS:',ic,raxis(ic),dr 
         end do 
c
c        Row center coordinates
c
         dz = (yymax - yymin) / REAL(yres)
         do ir = 1,yres  
            zaxis(ir) = yymin + REAL(ir - 1) * dz + 0.5 * dz 
            write(6,'(a,i5,2(1x,g12.5))') 'ZAXIS:',ir,zaxis(ir),dz 
         end do 

c...     Assign data to image array from 2D cell plot:

c...     Method2: Some weighted average of nearby data points.  If the XY cell center is not 
c        inside the grid, then assign a null value:

         DO ir = 1, yres
           DO ic = 1, xres
 
             tpt(ic,ir) = 0
             
             xcen = raxis(ic)
             ycen = zaxis(ir)
             centerinside = .FALSE.
             i1 = 0
             DO WHILE(.NOT.centerinside.AND.i1.LT.nc)
               i1 = i1 + 1
               IF (cq(i1).EQ.LO) CYCLE
c...           Decide if cell center is inside or outside the parent cell polygon:
               icount = 0
               DO i2 = 1, nv(i1)
                 i3 = i2 + 1
                 IF (i2.EQ.nv(i1)) i3 = 1
                 CALL CalcInter
     .             (DBLE(rv(i2,i1)),DBLE(zv(i2,i1)),
     .              DBLE(rv(i3,i1)),DBLE(zv(i3,i1)),
     .              DBLE(xcen),DBLE(ycen),DBLE(xcen+1000.0),DBLE(ycen),
     .              t1,t2)
                 IF (t1.GT.0.0D0.AND.t1.LT.1.0D0.AND.
     .               t2.GT.0.0D0.AND.t2.LT.1.0D0) icount = icount + 1
               ENDDO
               IF (icount.EQ.0.OR.MOD(icount,2).EQ.0) THEN
               ELSE
                 centerinside = .TRUE.
                 tpr(ic,ir) = rq(i1)
               ENDIF
             ENDDO
             IF (centerinside) tpt(ic,ir) = i1
           ENDDO
         ENDDO

c...     Find the center point of each cell:
         DO i1 = 1, nc
           rc(i1) = 0.0
           zc(i1) = 0.0
           DO i2 = 1, nv(i1)
             rc(i1) = rc(i1) + 1.0 / REAL(nv(i1)) * rv(i2,i1)
             zc(i1) = zc(i1) + 1.0 / REAL(nv(i1)) * zv(i2,i1)
           ENDDO
         ENDDO           

c...     Find the 5 closest points with CQ greater than QMIN and 
c        average:

         DO ir = 1, yres
           DO ic = 1, xres

             image(ic,ir) = 0.0

             WRITE(6,'(A,4I6)') 'IMAGE:',ic,ir,tpt(ic,ir),tpr(ic,ir)

c *TEMP*   
             IF (tpt(ic,ir).EQ.0) CYCLE

             xcen = raxis(ic)
             ycen = zaxis(ir)

             npt = 1
             vpt(npt) = 0.0
             dpt(npt) = HI

             DO i1 = 1, nc

c *TEMP*
               IF (cq(i1).LT.lmin+(1-lmin)) CYCLE

c               IF (cq(i1).LT.qmin) CYCLE

c...           Don't average cells across different additional cell
c              regions:
c *TEMP*
               IF (rq(i1).NE.0.AND.tpr(ic,ir).NE.0.AND.
     .             rq(i1).NE.tpr(ic,ir)) CYCLE

               dist = SQRT((xcen - rc(i1))**2 + (ycen - zc(i1))**2)

               IF (dist.LT.dr+dz.OR.
     .             (rq(i1).NE.0.AND.dist.LT.0.01)) THEN
                 IF (npt.LT.100) npt = npt + 1
                 dpt(npt) = dist 
                 vpt(npt) = cq(i1) 
c...             Sort data points from closest to farthest:                 
                 i2 = 1
                 DO WHILE(i2.LT.npt)
                   i2 = i2 + 1
                   IF (dpt(i2-1).GT.dpt(i2)) THEN
                     tmpdpt = dpt(i2-1)
                     tmpvpt = vpt(i2-1)
                     dpt(i2-1) = dpt(i2)
                     vpt(i2-1) = vpt(i2)
                     dpt(i2) = tmpdpt
                     vpt(i2) = tmpvpt
                     i2 = 1   
                   ENDIF
                 ENDDO
               ENDIF
             ENDDO

c...         No value identified for cell:
             IF (xval.NE.-99.OR.(npt.EQ.1.AND.vpt(1).EQ.0.0)) THEN
               npt = 1
               vpt(1) = cq(tpt(ic,ir))
               dpt(1) = 1.0
             ENDIF

c...         Average the data points:
             totdpt = 0.0
             DO i1 = 1, npt
               dpt(i1) = 1.0 / (dpt(i1)**2)
               totdpt = totdpt + dpt(i1)
             ENDDO
             IF (totdpt.GT.0.0) THEN
               DO i1 = 1, npt
                 image(ic,ir) = image(ic,ir) + dpt(i1)/totdpt * vpt(i1)
               ENDDO
             ELSE
               image(ic,ir) = 0.0
             ENDIF
           ENDDO                          
         ENDDO

175      CONTINUE
c
c        Draw contour plot 
c

         CALL CTRMAG(10)
 
c...     Have GRTSET_TRIM called:

c         slopt4 = 1
c         char(30) = ' '
c         REF    = '                                    '

c
c         jdemod - commented out write statement since IPLOT has 
c                  not been assigned a value
c
c         WRITE (IPLOT,9012) NPLOTS,REF
c
c         CALL GRTSET (TITLE,REF,NVIEW,PLANE,blabs(1),
c     .         xXMIN,xXMAX,
c     >         yYMIN,yYMAX,
c     .         TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,Ncntr)
          CALL HSV

c...      Clean up IMAGE array:
          DO i1 = 1, xres
            DO i2 = 1, yres
              IF (image(i1,i2).GT.-1.0E-10.AND.
     .            image(i1,i2).LT. 1.0E-10) image(i1,i2) = 0.0
            ENDDO
          ENDDO

c...      Erase boarder:
          DO i1 = 1, xres
            image(i1,1) = 0.0
            image(i1,yres) = 0.0
          ENDDO
          DO i2 = 1, yres
            image(1,i2) = 0.0
            image(xres,i2) = 0.0
          ENDDO

c...      Sort contours in increasing absolute value:
          status = .TRUE.
          DO WHILE (status)
            status = .FALSE.
            DO i1 = 1, ncntr-1
              IF (ABS(uconts(i1)).GT.ABS(uconts(i1+1))) THEN
                rdum1 = uconts(i1+1)
                uconts(i1+1) = uconts(i1)
                uconts(i1) = rdum1
                status = .TRUE.
              ENDIF
            ENDDO
          ENDDO

          WRITE(0,*) 'NCTR:',ncntr

          DO i1 = 1, ncntr

c            frac = (uconts(i1) - qmin) / (qmax - qmin)
c            frac5 = 100.0*frac
c            fmod5 = AMOD(frac5,2.0)
c            frac = MIN(0.98,(frac5-fmod5)/100.0)
c            CALL SetCol255(frac,qmin,qmax)



            frac = MAX(0.0,(uconts(i1) - (1 - lmin)  - lmin) / 
     .                     (lmax - lmin))
            frac = 100.0 * 0.98 * frac
            fmod5 = AMOD(frac,1.0)
            frac = (frac - fmod5) / 100.0

c...dev
          bright = 1.0-(1.0-frac)**10
          frac = (1.0 - frac) * 0.90
          frac = frac + 0.34
          IF (frac.GT.1.0) frac = frac - 1.0


          val = 10.0**(uconts(i1) - (1 - lmin))
c          WRITE(0,*) 'FRAC:',val
          IF (val.LT.1.0) THEN
            decade = REAL(INT(LOG10(val)-0.99999)) 
          ELSE
            decade = REAL(INT(LOG10(val)))
          ENDIF

          decval = 10.0**decade + LO
c          IF (val.LT.1.0) decade = decade - 1.0
c          decval = 10.0**decade + LO
c          val = (REAL(INT(val / decval)) + 0.5) * decval
c          decval = 10.0**decade + LO
c          WRITE(0,*) 'LMAX:',lmax,lmin,lmax-lmin

          hue = 1.0
          IF (grayscale) THEN
            bright = 0.0
c...         Colour:
c            frac = 0.66 * (LOG10(val) - lmin) / (lmax - lmin) + 0.34
c            frac = 1.0 - frac + 0.34
            frac = (LOG10(val) - lmin) / (lmax - lmin)
            frac = 0.9 * frac + 0.1
            WRITE(0,*) 'FRAC C:',frac,decade,decval,val,LOG10(val)

          ELSEIF (lmax-lmin.LE.2) THEN
            bright = 1.0
            frac = REAL(decade - lmin) / REAL(lmax - 1 - lmin) * 0.40 
            frac = frac + (val - decval) / (8.5 * decval) * 0.40
            WRITE(0,*) 'FRAC B:',frac,decade,decval,val,LOG10(val)
          ELSE
c            hue = 0.5*(1.0 - (val - decval) / (8.5 * decval)) + 0.5
c            hue = SQRT(hue)
            bright = (val - decval) / (9.0 * decval) * 0.5 + 0.5
c            bright = (val - decval) / (8.5 * decval) * 0.5 + 0.5
c            bright = (val - decval) / (8.5 * decval) * 0.3 + 0.7
 
            frac = REAL(decade - lmin) / REAL(lmax - 1 - lmin) * 0.65

c            IF (iopt.EQ.54.AND.decade.EQ.0) frac = frac + 0.10

c            WRITE(0,*) 'FRAC A:',frac,bright,decval,val
          ENDIF

c...new
c          val = cq(i1)
c          decade = REAL(INT(LOG10(val)))
c          IF (decade.LT.0) decade = decade - 1.0
c          decval = 10.0**decade + LO
c          val = (REAL(INT(val / decval)) + 0.5) * decval
c          decval = 10.0**decade + LO
c          bright = (val - decval) / (8.5 * decval) * 0.3 + 0.7
c          frac = REAL(decade - lmin) / REAL(lmax - 1 - lmin) * 0.65

            IF (grayscale) THEN
              CALL ColSet(0.0,0.0,1.0-frac,255)
            ELSEIF (frac.LT.0.0) THEN
              CALL RGB
              CALL ColSet(0.0,1.0,0.0,255)
              CALL HSV
            ELSEIF (frac.GT.0.65) THEN
              CALL RGB
              CALL ColSet(1.0,0.0,0.0,255)
              CALL HSV
            ELSE
              CALL ColSet(frac,hue,bright,255) 
            ENDIF
c            CALL ColSet(frac,1.0,1.0-(1.0-frac)**20,255)
            CALL FILCOL(255)
            CALL LINCOL(255)


            map1x_d = map1x
            map2x_d = map2x
            map1y_d = map1y
            map2y_d = map2y

            tag_d = 1

c...        Copy IMAGE to IMAGE1 and account for 
c           negative contours:
            IF (uconts(i1).LT.0.0) THEN
              DO i2 = 1, xres
                DO i3 = 1, yres
                  IF (image(i2,i3).LT.-1.0E-10) THEN
                    image1(i2,i3) = -image(i2,i3)
                  ELSE
                    image1(i2,i3) = 0.0
                  ENDIF
                ENDDO
              ENDDO
            ELSE       
              DO i2 = 1, xres
                DO i3 = 1, yres
                  image1(i2,i3) = image(i2,i3)
                ENDDO
              ENDDO
            ENDIF

c...        Plot contour:
            CALL GRCONT (image1,1,xres,xres,1,yres,yres,
     >                   ABS(uconts(i1)),raxis,zaxis,'hello')

            tag_d = 0
          ENDDO

c...      Finish off the plot boarder:
          CALL FULL
          CALL LINCOL(defcol)
          CALL PSPACE(MAP1X,MAP2X,MAP1Y,MAP2Y)
          CALL MAP(0.0,1.0,0.0,1.0)
          CALL POSITN (0.0,1.0)
          CALL JOIN   (1.0,1.0)
          CALL JOIN   (1.0,0.0)
c
c
         deallocate(raxis)
         deallocate(zaxis)
         deallocate(image)
         deallocate(image1)
         deallocate(tpt)

200   CONTINUE




      READ(5,'(A256)') dummy
      IF   (dummy(8:14).EQ.'Noscale'.OR.dummy(8:14).EQ.'noscale'.OR.
     .      dummy(8:14).EQ.'NOSCALE') THEN
        GOTO 210
      ELSE
        BACKSPACE 5
      ENDIF

      IF (elabs(1)(1:4).EQ.'none') GOTO210




c...  Draw scale:
      CALL PSPACE (0.0, 1.35, 0.0,1.0)
      CALL CSPACE (0.0, 1.35, 0.0,1.0)
      CALL MAP    (0.0, 1.35, 0.0,1.0)
      CALL FULL
      CALL HSV

      DO decade = lmin, lmax-1
        decval = 10.0**decade + LO
        DO val = 1.5*decval, 9.6*decval, decval
          t1a = MAX(0.0,(LOG10(val - 0.5 * decval) - lmin) / 
     .                  (lmax - lmin))
          t2a =         (LOG10(val + 0.5 * decval) - lmin) / 
     .                  (lmax - lmin)
          y1 = map1y + t1a * (map2y - map1y)
          y2 = map1y + t2a * (map2y - map1y)
c          y1 = 0.42 + t1a * (0.82 - 0.42)
c          y2 = 0.42 + t2a * (0.82 - 0.42)
c...dev
c          WRITE(0,*) 'VAL 2:',val
          hue = 1.0
          IF (grayscale) THEN
            bright = 0.0
c...         Colour:
c            frac = 0.66 * (LOG10(val) - lmin) / (lmax - lmin) + 0.34
c            frac = 1.0 - frac + 0.34
            frac = (LOG10(val) - lmin) / (lmax - lmin)
            frac = 0.9 * frac + 0.1
            WRITE(0,'(A,5F10.3)') 
     .        'FRAC 2:',frac,decade,decval,val,LOG10(val)
          ELSEIF (lmax-lmin.LE.2) THEN
            bright = 1.0
            frac = REAL(decade - lmin) / REAL(lmax - 1 - lmin) * 0.40 
            frac = frac + (val - decval) / (8.5 * decval) * 0.40
            WRITE(0,*) 'FRAC 2:',frac,decade,decval,val
          ELSE
c            hue = 0.5*(1.0 - (val - decval) / (8.5 * decval)) + 0.5
c            hue = SQRT(hue)
            bright = (val - decval) / (8.5 * decval) * 0.5 + 0.5
c            bright = (val - decval) / (8.5 * decval) * 0.3 + 0.7
            frac = REAL(decade - lmin) / REAL(lmax - 1 - lmin) * 0.65 
          ENDIF
c          WRITE(0,*) 'SCALE:',frac,bright,val,hue
          IF (grayscale) THEN
            CALL ColSet(0.0,0.0,1.0-frac,255)
          ELSE
            CALL ColSet(frac,hue,bright,255)
          ENDIF
          CALL FILCOL(255)
          CALL LINCOL(255)
c...      Draw coloured box:
c          DSPOT = 0.016
c          CALL BOX (0.98-DSPOT,0.98+DSPOT,y1,y2)
          DSPOT = (MAP2Y - MAP1Y) * 0.02
c          DSPOT = (MAP2Y - MAP1Y) * 0.0373
          CALL BOX (MAP2X+0.03*(MAP2X-MAP1X),
     .              MAP2X+0.03*(MAP2X-MAP1X)+2.0*DSPOT,
     .              y1,y2)

c...      Put labels on scale:
          CALL ColSet(1.0,1.0,0.0,255)
          IF (val.LT.1.6*decval) THEN
            DIST = 0.045*(MAP2X-MAP1X)+2.0*DSPOT+0.010
            CALL PLOTST(MAP2X+DIST,y1,'10')
            CALL POSITN(MAP2X+DIST+0.015,y1+0.015)
            CALL TYPENI(INT(decade))
            CALL POSITN(MAP2X+0.03*(MAP2X-MAP1X)+2.0*DSPOT      ,y1)
            CALL JOIN  (MAP2X+0.03*(MAP2X-MAP1X)+2.0*DSPOT+0.005,y1)
          ENDIF
          IF (decade.EQ.lmax-1.AND.val.GT.9.4*decval) THEN
            DIST = 0.045*(MAP2X-MAP1X)+2.0*DSPOT+0.010
            CALL PLOTST(MAP2X+DIST,y2,'10')
            CALL POSITN(MAP2X+DIST+0.015,y2+0.015)
            CALL TYPENI(INT(decade+1.0))
            CALL POSITN(MAP2X+0.03*(MAP2X-MAP1X)+2.0*DSPOT      ,y2)
            CALL JOIN  (MAP2X+0.03*(MAP2X-MAP1X)+2.0*DSPOT+0.005,y2)
          ENDIF

        ENDDO
      ENDDO

c...  Draw box around scale:
      CALL POSITN(MAP2X+0.03*(MAP2X-MAP1X),map1y)
      CALL JOIN  (MAP2X+0.03*(MAP2X-MAP1X)+2.0*dspot,map1y)
      CALL JOIN  (MAP2X+0.03*(MAP2X-MAP1X)+2.0*dspot,map2y)
      CALL JOIN  (MAP2X+0.03*(MAP2X-MAP1X),map2y)
      CALL JOIN  (MAP2X+0.03*(MAP2X-MAP1X),map1y)

c...  Draw contour plot label:
      CALL CTRORI(90.0)
      CALL CTRMAG(14)      
      CALL PLOTST(MAP2X+0.03*(MAP2X-MAP1X)+2.0*DSPOT+0.075,map1y,
     .            elabs(1)(1:LEN_TRIM(elabs(1))))
      CALL CTRORI(0.0)


210   CONTINUE

c...  Print comments:
      CALL PSPACE (0.0, 1.35, 0.0,1.0)
      CALL CSPACE (0.0, 1.35, 0.0,1.0)
      CALL MAP    (0.0, 1.35, 0.0,1.0)
      CALL CTRMAG (12)
      DO i = 1, 10
        IF (char(i).NE.' ') CALL PLOTST(1.00,0.590+(i-1)*0.02,char(i))
      ENDDO
      DO i = 20, 30
        j = i - 19
        IF (char(i).NE.' ') CALL PLOTST(1.00,0.50-(j-1)*0.02,char(i))
      ENDDO

c...  Plot grid outline:
c      IF (xval.NE.-99.0) THEN
c      ELSE
c        CALL FULL
c        CALL LINCOL(1)
c        CALL SUPIMP('PARTIAL')      
c      ENDIF

c...  Add a comment to the plot:
79    READ(5,'(A5000)') dummy
      IF   (dummy(8:11).EQ.'Cmnt'.OR.dummy(8:11).EQ.'cmnt'.OR.
     .      dummy(8:11).EQ.'CMNT') THEN
        READ (dummy,*) cdum1,xpos,ypos,size,caption
        IF (TRIM(caption).EQ.'<charge state>') THEN
          WRITE(caption,'(A,I2,1000X)') '+',iz
        ENDIF
c...    Annotate graph:
        CALL PSPACE (map1x,map2x,map1y,map2y)
        CALL MAP    (0.0,1.0,0.0,1.0)
        CALL LinCol(1)
        CALL CTRMAG(size)
        CALL PLOTST(xpos,ypos,caption(1:LEN_TRIM(caption)))
c...    Another comment:        
        GOTO 79
      ELSE
        BACKSPACE 5
      ENDIF

c...  Add a caption to the plot:
      READ(5,'(A5000)') dummy
      IF   (dummy(8:11).EQ.'Note'.OR.dummy(8:11).EQ.'note'.OR.
     .      dummy(8:11).EQ.'NOTE') THEN
        READ(dummy,*) cdum1,xpos,ypos,size,caption
        CALL AddCaption(caption,xpos,ypos,size)
      ELSE
        BACKSPACE 5
      ENDIF

c...  Plot grid outline:
      IF (xval.NE.-99.0) THEN

      ELSE
        iplots = 0
        CALL FULL
        CALL LINCOL(55)
        slopt  = 1
        slopt2 = 0
        CALL SUPIMP('PARTIAL')      
      ENDIF


c...  Plot a line marking the injection location:
      IF (ciopte.EQ.9) THEN
        CALL PSPACE(map1x,map2x,map1y,map2y)
        CALL MAP   (cxmin,cxmax,cymin,cymax)
        CALL BROKEN(7,7,7,7)
        CALL POSITN(cxsca,cysca)
        CALL JOIN  (cxscb,cyscb)
        CALL FULL
      ENDIF



      READ(5,'(A256)') dummy
      IF   (dummy(8:14).EQ.'Noframe'.OR.dummy(8:14).EQ.'noframe'.OR.
     .      dummy(8:14).EQ.'NOFRAME') THEN
      ELSE
        CALL FRAME
        BACKSPACE 5
      ENDIF


c...  Write plot data to the OUT data file:
c      WRITE(6,*) '982: DATA'
c      DO i1 = 1, nc
c        WRITE(6,'(I6,1P,E12.4,0P)') i1,cq(i1)
c      ENDDO


      IF (zval.NE.-99.0) zval = -99.0

      IF (MARsum.NE.0.0) WRITE(0,*) 'MARSUM:',MARsum

      RETURN
 9012 FORMAT(1X,'PLOT',I3,4X,A)
 99   WRITE(0,*) 'DUMMY >'//TRIM(dummy)//'<'
      STOP
      END
