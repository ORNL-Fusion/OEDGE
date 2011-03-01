c     -*-Fortran-*-
c
c ======================================================================
c
c function PsinToR
c
      REAL FUNCTION PsinToR(psin,z)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      REAL psin,z

      INTEGER ir,ir0,ir2,ik,ir1,irshift,i1,id
      REAL    psi0,psi1,psi2,frac,minz,maxz,minr,maxr

c...  Add check for broken grid since the IR progression will fail...     

      PsinToR = 0.0


      ir = irtrap + 1
      IR_LOOP: DO WHILE (ir.NE.irwall)
        irshift = 0
c        IF (ir.GT.irwall.AND.z.GT.zxp) irshift = -irtrap + 1 
c      IF (irsep-2.NE.nrs-irtrap.AND.z.GT.zxp) 
c     .  CALL ER('PsinToR','This routine only works properly with '//    ! If trying to use the above condition, need to have
c     .          'balanced grids',*99)                                   ! this check


        ir1 = ir + irshift
        ir0 = irins (nks(ir1),ir1)
        ir2 = irouts(nks(ir1),ir1)

        psi0 = psitarg(ir0,1)
        psi1 = psitarg(ir1,1)
        psi2 = psitarg(ir2,1)

        IF (ir.EQ.irtrap+1) THEN
          psi0 = psi1                                          ! Inner half of ring gets missed... 
        ELSE
          frac = (rho(ir1,IN14 ) - rho(ir1,CELL1)) / 
     .           (rho(ir0,CELL1) - rho(ir1,CELL1))
          psi0 = frac * psi0 + (1.0 - frac) * psi1
        ENDIF

        IF (ir.EQ.irwall-1) THEN
          psi2 = psi1
        ELSE
          frac = (rho(ir1,OUT23) - rho(ir1,CELL1)) / 
     .           (rho(ir2,CELL1) - rho(ir1,CELL1))
          psi2 = frac * psi2 + (1.0 - frac) * psi1
        ENDIF

        IF (psin.GE.psi0.AND.psin.LE.psi2) THEN
c...      Find R:

          IF (irshift.NE.0) WRITE(0,*) 'FOUND CORE RING (CHECK)'

          IF (psi1.LT.1.0.AND.irshift.EQ.0) 
     .      WRITE(6,*) '   FOUND PFZ RING',ir1

          DO ik = nks(ir1), 2, -1
c          DO ik = nks(ir1), nks(ir1)/2, -1

c...        Find z limits:
            id = korpg(ik,ir1)
            minz =  HI
            maxz = -HI
            DO i1 = 1, nvertp(id)
              IF (zvertp(i1,id).LT.minz) THEN
                minz = zvertp(i1,id)
                minr = rvertp(i1,id)
              ENDIF
              IF (zvertp(i1,id).GT.maxz) THEN
                maxz = zvertp(i1,id)
                maxr = rvertp(i1,id)
              ENDIF
            ENDDO

            IF (z.GE.minz.AND.z.LT.maxz) THEN

              frac = (z - minz) / (maxz - minz)
              PsinToR = (1.0 - frac) * minr + frac * maxr

c            IF (z.GE.zs(ik,ir1).AND.z.LT.zs(ik-1,ir1)) THEN
c              frac = (z            - zs(ik,ir1)) /                   
c     .               (zs(ik-1,ir1) - zs(ik,ir1))                     
c              PsinToR = frac * rs(ik-1,ir1) + (1.0 - frac) * rs(ik,ir1) 
              EXIT IR_LOOP
            ENDIF
          ENDDO
c          WRITE(0,*) 'DTS MAP PROBLEM:',ir,psin,z
c          CALL ER('PsinToR','Unable to find Z',*99)
        ENDIF

        ir = irouts(nks(ir),ir)
      ENDDO IR_LOOP


      WRITE(6,*) 'PSIn IR  :',ir1,PsinToR,psin,z
          


c      PsinToR = 1.45 



      RETURN
 99   WRITE(0,*) 'IRSEP,IRTRAP,NRS=',irsep,irtrap,nrs
      STOP
      END
c
c
c ======================================================================
c
c
c
c ======================================================================
c
c subroutine LoadRCPData
c
c input:
c
c output:
c   raw - 1 - R
c         2 - Z
c         3 - cell index
c         4 - ring index
c         5 - denisty
c         6 - temperature
c
c
      SUBROUTINE LoadRCPData(graph,nraw,raw,MAXTDAT,MAXCOLS,
     .                       xshift1,yshift1)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      INTEGER    DATAUNIT   ,MAXBLOCK
      PARAMETER (DATAUNIT=13,MAXBLOCK=100)

c...  Input:
      INTEGER   MAXTDAT,MAXCOLS
      REAL      xshift1,yshift1
      CHARACTER graph*(*)
c...  Output:
      INTEGER nraw
      REAL    raw(MAXTDAT,MAXCOLS)

      INTEGER   etype,ndata,ncol,i1,i2,ik,ir,iblock,id,i,optmap
      LOGICAL   outofgrid,intersect
      REAL      xdata(MAXTDAT),edata(MAXTDAT,MAXCOLS),t,psin,
     .          count(MAXNKS,MAXNRS),f,val,xshift2,yshift2
      REAL*8    a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd      
      CHARACTER datatitle*128,cdum1*100

      CALL RZero(count,MAXNKS*MAXNRS)


      nraw     = 0


c...  Load parameters from OUT input file:
      READ(graph,*) cdum1,cdum1,xshift2,yshift2,optmap,iblock

      WRITE(6,*) 'LOAD RCP: ',xshift2,yshift2,optmap,iblock

c...  Read target data:
      IF (iblock.GT.0) THEN

        CALL Load_Expt_Data(DATAUNIT,iblock,xdata,etype,edata,
     .                      MAXCOLS,MAXTDAT,ndata,ncol,datatitle)

        WRITE(6,*) '============================================='
        WRITE(6,*) 'TITLE  = ',datatitle(1:LEN_TRIM(datatitle))
        WRITE(6,*) 'INDEX  = ',iblock
        WRITE(6,*) 'TYPE   = ',etype
        WRITE(6,*) 'NDATA  = ',ndata
        WRITE(6,*) 'NCOL   = ',ncol

      ELSE
        CALL ER('LoadRCPData','Invalid IBLOCK value',*99)
      ENDIF

      IF     (optmap.EQ.1) THEN
c...    DIII-D horizontal RCP probed at Z = -0.190 m, mapping from PSIn to R,Z:

c...    Assign shift to the location of the data:
        IF (xshift1.NE.0.0.OR.xshift2.NE.0.0.OR.
     .      yshift1.NE.0.0.OR.yshift2.NE.0.0) 
     .    CALL ER('LoadRCPData','Data shifting not ready',*99)

        DO ir = irsep, irwall-1
c...      Find intersection of the data with IR:
          a1 = DBLE(r0)
          a2 = -0.190
          b1 = DBLE(r0 + 100.0)
          b2 = -0.190
          intersect = .FALSE.
          DO ik = 1, nks(ir)
            id = korpg(ik,ir)
            c1 = DBLE(0.5 * (rvertp(1,id) + rvertp(2,id)))
            c2 = DBLE(0.5 * (zvertp(1,id) + zvertp(2,id)))
            d1 = DBLE(0.5 * (rvertp(3,id) + rvertp(4,id)))
            d2 = DBLE(0.5 * (zvertp(3,id) + zvertp(4,id)))
            CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
            IF (tab.GE.0.0.AND.tab.LE.1.0.AND.
     .          tcd.GE.0.0.AND.tcd.LE.1.0) THEN
              IF (.NOT.intersect) THEN
                nraw = nraw + 1
                raw(nraw,1) = SNGL(c1 + tcd * (d1 - c1))
                raw(nraw,2) = SNGL(c2 + tcd * (d2 - c2))
                raw(nraw,3) = REAL(ik)
                raw(nraw,4) = REAL(ir)
                intersect = .TRUE.
              ELSE
                CALL ER('LoadRCPData','Multiple intersections '//
     .                  'found',*99)
              ENDIF
            ENDIF
          ENDDO        
c...      Interpolate data based on PSIn from experimental data file and 
c         PSIn (outer target for DIII-D):
          IF (intersect) THEN

            psin = psitarg(ir,1)

            t = 0.0
            DO i1 = 1, ndata-1
              IF ((xdata(i1  ).LE.psin.AND.
     .             xdata(i1+1).GE.psin).OR.
     .            (xdata(i1  ).GE.psin.AND.
     .             xdata(i1+1).LE.psin)) THEN
                t = (psin - xdata(i1)) / (xdata(i1+1) - xdata(i1))
                raw(nraw,5) =edata(i1,1) + t*(edata(i1+1,1)-edata(i1,1))
                raw(nraw,6) =edata(i1,2) + t*(edata(i1+1,2)-edata(i1,2))
              ENDIF
            ENDDO
            IF (t.EQ.0.0) 
     .        CALL ER('LoadRCPData','No data interpolation '//
     .                'possible',*99)
          ENDIF
        ENDDO

        DO i1 = 1, nraw
          WRITE(0,'(A,2F10.4,2F5.0,1P,E10.2,0P,F10.2)')
     .      'RCP DATA:',(raw(i1,i2),i2=1,6)
        ENDDO

      ELSEIF (.FALSE.) THEN

        DO i1 = 1, ndata
          xdata(i1)   = xdata(i1)   + xshift1 + xshift2
          edata(i1,1) = edata(i1,1) + yshift1 + yshift2
          ik = 1
          ir = 1
          CALL GridPos(ik,ir,xdata(i1),edata(i1,1),.FALSE.,outofgrid)
          IF (.NOT.outofgrid) THEN
            IF (edata(i1,2).NE.0.0.AND.edata(i1,3).NE.0.0) THEN
              nraw = nraw + 1
              raw(nraw,1) = xdata(i1)
              raw(nraw,2) = edata(i1,1)
              raw(nraw,3) = REAL(ik)
              raw(nraw,4) = REAL(ir)
              raw(nraw,5) = edata(i1,2)
              raw(nraw,6) = edata(i1,3)
            ENDIF
          ENDIF
          WRITE(6,90) ' TD : ',i1,xdata(i1),edata(i1,1),ik,ir,
     .                outofgrid,rs(ik,ir),zs(ik,ir),
     .                xdata(i1),(edata(i1,i2),i2=1,ncol)
90        FORMAT(A,I4,2F8.4,2I4,L2,2F8.4,1P,4E10.2,0P,F6.3)
        ENDDO

      ELSE
        CALL ER('LoadRCPData','Invalid OPTMAP value',*99)
      ENDIF

c...  Return the total shift in XSHIFT1 for plot notes:
c      xshift1 = xshift1 + xshift2

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: LoadTargetData
c
c This routine loads the Langmuir probe target data that is stored
c in the .experiment file.
c
c
c
      SUBROUTINE LoadTargetData(graph,targdat,MAXTDAT,MAXCOLS)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'
      INCLUDE 'slout'

      INTEGER    DATAUNIT   ,MAXBLOCK
      PARAMETER (DATAUNIT=13,MAXBLOCK=100)

c...  Input:
      CHARACTER*80 graph
      INTEGER      MAXTDAT,MAXCOLS

c...  Output:
      REAL targdat(MAXNDS,5)

      REAL GetNe

      INTEGER   etype,ndata1,ndata2,ncol,i1,i2,ik,ir,lblock,hblock,in
      LOGICAL   swjsat
      REAL      xdata(MAXTDAT),edata(MAXTDAT,MAXCOLS),
     .          lshift,hshift,temp_lshift,temp_hshift,
     .          temp_lpdati(MAXINS,4),temp_lpdato(MAXINS,4)
      CHARACTER datatitle*128,dataline*128,cdum1*128

c...  Store current target data:
      temp_lshift = tarshift(IKLO)
      temp_hshift = tarshift(IKHI)
      DO i1 = 1, MAXINS
        DO i2 = 1, 4
          temp_lpdati(i1,i2) = lpdati(i1,i2)      
          temp_lpdato(i1,i2) = lpdato(i1,i2)      
        ENDDO
      ENDDO

c...  Load parameters from OUT input file:
      READ(graph,*) cdum1,lshift,hshift,lblock,hblock

      WRITE(6,*) 'LOAD TARG: ',lshift,hshift,lblock,hblock

c...  Use shifts specified in the OUT input file, instead of those
c     in the DIVIMP input file:
      IF (lshift.NE.99.0) tarshift(IKLO) = lshift
      IF (hshift.NE.99.0) tarshift(IKHI) = hshift

c...  Load low index target data:
      IF (lblock.GT.0) THEN
        ndata1 = 0
        CALL Load_Expt_Data(DATAUNIT,lblock,xdata,etype,edata,
     .                      MAXCOLS,MAXTDAT,ndata1,ncol,datatitle)

        WRITE(6,*) '============================================='
        WRITE(6,*) 'TITLE  = ',datatitle(1:LEN_TRIM(datatitle))
        WRITE(6,*) 'INDEX  = ',lblock
        WRITE(6,*) 'TYPE   = ',etype
        WRITE(6,*) 'NDATA  = ',ndata1
        WRITE(6,*) 'NCOL   = ',ncol

        DO i1 = 1, ndata1
          WRITE(6,'(A,I4,F8.1,2F10.4,1P,E12.4,0P)')
     .      ' LO -> ',i1,xdata(i1),edata(i1,1),edata(i1,2),edata(i1,3)
        ENDDO

        nlpdato = ndata1
        DO i1 = 1, ndata1
          lpdato(i1,1) = xdata(i1)
          lpdato(i1,2) = edata(i1,1)
          lpdato(i1,3) = edata(i1,2)
          lpdato(i1,4) = edata(i1,3)
        ENDDO

      ELSEIF (lblock.EQ.-1) THEN
c...    Use target data passed via the raw file:
        lshift = 0.0
        nlpdato = nrs - irsep + 1
        ndata1  = nlpdato
        DO i1 = 1, ndata1
          lpdato(i1,1) = REAL(i1+irsep-1)
          lpdato(i1,2) = kteds(idds(i1+irsep-1,2))
          lpdato(i1,3) = ktids(idds(i1+irsep-1,2))
          lpdato(i1,4) = knds (idds(i1+irsep-1,2))
        ENDDO
      ELSE
        CALL ER('LoadTargetData','Invalid LBLOCK',*99)
      ENDIF

      IF (hblock.GT.0) THEN
c...    Load high index target data:
        ndata2 = 0
        CALL Load_Expt_Data(DATAUNIT,hblock,xdata,etype,edata,
     .                      MAXCOLS,MAXTDAT,ndata2,ncol,datatitle)

        WRITE(6,*) '============================================='
        WRITE(6,*) 'TITLE  = ',datatitle(1:LEN_TRIM(datatitle))
        WRITE(6,*) 'INDEX  = ',lblock
        WRITE(6,*) 'TYPE   = ',etype
        WRITE(6,*) 'NDATA  = ',ndata2
        WRITE(6,*) 'NCOL   = ',ncol

        DO i1 = 1, ndata2
          WRITE(6,'(A,I4,F8.1,2F10.4,1P,E12.4,0P)')
     .      '  HI -> ',i1,xdata(i1),edata(i1,1),edata(i1,2),edata(i1,3)
        ENDDO

        nlpdati = ndata2
        DO i1 = 1, nlpdati
          lpdati(i1,1) = xdata(i1)
          lpdati(i1,2) = edata(i1,1)
          lpdati(i1,3) = edata(i1,2)
          lpdati(i1,4) = edata(i1,3)
        ENDDO

        WRITE(char(8),10) 'OUTER LANGMUIR DATA SHIFT =',
     .                    INT(tarshift(IKHI)*1000.0),' mm'
10      FORMAT(A,I4,A)
       
      ELSEIF (hblock.EQ.-1) THEN
c...    Use target data passed via the raw file:
        hshift = 0.0
        nlpdati = nrs - irsep + 1
        ndata2  = nlpdati
        DO i1 = 1, ndata2
          lpdati(i1,1) = REAL(i1+irsep-1)
          lpdati(i1,2) = kteds(idds(i1+irsep-1,1))
          lpdati(i1,3) = ktids(idds(i1+irsep-1,1))
          lpdati(i1,4) = knds (idds(i1+irsep-1,1))
        ENDDO
      ELSE
        CALL ER('LoadTargetData','Invalid HBLOCK',*99)
      ENDIF


      IF (ndata1.NE.ndata2) 
     .  CALL WN('LoadTargetData','Different number of data points for'//
     .          ' inner and outer targets.  Check the experimental'//
     .          ' data file.')

c...  Do the shifting:
      WRITE(6,*) 'SHIFT: ',tarshift(IKLO),tarshift(IKHI)
      CALL ShiftTargetData

c...  Decide if target data in the .experiment file uses Jsat or Ne:
      DO i1 = 1, ndata1
        IF (INT(lpdato(i1,1)).EQ.irsep) THEN
          IF (lpdato(i1,4).GT.1.0E+10) THEN
            swjsat = .FALSE. 
            WRITE(6,*) 'DENSITY DETECTED'
          ELSE
            swjsat = .TRUE.
            WRITE(6,*) 'JSAT DETECTED'
          ENDIF
        ENDIF
      ENDDO

c...  Assign the plotting array using the shifted target data:
      DO i1 = 1, nlpdato
        in = idds(NINT(lpdato(i1,1)),2)
        targdat(in,1) = lpdato(i1,2)
        targdat(in,2) = lpdato(i1,3)
        IF (swjsat) THEN
          targdat(in,3) = GetNe(lpdato(i1,2),lpdato(i1,3),
     .                          lpdato(i1,4),1.0)
        ELSE
          targdat(in,3) = lpdato(i1,4)
        ENDIF
      ENDDO

      DO i1 = 1, nlpdati
        in = idds(INT(lpdati(i1,1)),1)
        targdat(in,1) = lpdati(i1,2)
        targdat(in,2) = lpdati(i1,3)
        IF (swjsat) THEN
          targdat(in,3) = GetNe(lpdati(i1,2),lpdati(i1,3),
     .                          lpdati(i1,4),1.0)
        ELSE
          targdat(in,3) = lpdati(i1,4)
        ENDIF
      ENDDO

      DO i1 = 1, MAXNDS
        WRITE(6,'(A,I6,2I8,2F10.4,1P,E12.4,0P)')
     .    ' ..?',i1,irds(i1),ikds(i1),(targdat(i1,i2),i2=1,3)
      ENDDO

c...  Restore original target data:
      tarshift(IKLO) = temp_lshift
      tarshift(IKHI) = temp_hshift
      DO i1 = 1, MAXINS
        DO i2 = 1, 4
          lpdati(i1,i2) = temp_lpdati(i1,i2)      
          lpdato(i1,i2) = temp_lpdato(i1,i2)      
        ENDDO
      ENDDO

      RETURN
98    CALL ER('LoadTargetData','OUT input file access error',*99)
99    STOP
      END
c
c
c
c
c
      SUBROUTINE SetPlotComments(iref,job,extra_comments,ncomments,r1)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      include 'cedge2d'
      INCLUDE 'slcom'
      INCLUDE 'slout'

      INTEGER       iref,ncomments
      REAL          r1
      CHARACTER*(*) job
      character*(*) extra_comments(ncomments)
      
      INTEGER    idum1,idum2,idum3,i1,ncom
      CHARACTER filename*128

c      DO i1 = 1, 30
c        char(i1) = '          '
c      ENDDO

c...  Load case name from environment variable CASENAME:
      CALL IoName(filename,idum1)
c      CALL IoName(idum1,filename,idum2,idum3)


      IF (sldata.LT.1.6) THEN
        tarshift(IKHI) = -99.0

        te_mult_i = -99.0
        ti_mult_i = -99.0
        n_mult_i  = -99.0
      ENDIF

      IF (iref.EQ.296.OR.iref.EQ.980) THEN
c...    296: Dalpha LOS plot
c       980: Generic Dalpha LOS plots

        WRITE(char(20),10) 'OUTER TARGET DATA SHIFT =',
     .                    NINT(tarshift(IKHI)*1000.0),' mm'
        WRITE(char(21),10) 'VIEWING VERTEX R-SHIFT  =',
     .                     INT(tarshift(IKHI)*1000.0),' mm'
        WRITE(char(23),11) 'Te      MULTIPLIER =',te_mult_i
        WRITE(char(24),11) 'Ti      MULTIPLIER =',ti_mult_i
        WRITE(char(25),12) 'Ne/Jsat MULTIPLIER =',n_mult_i

        WRITE(char(27),13) filename(1:LEN_TRIM(filename))
        IF (rel_step.NE.0) 
     .  WRITE(char(28),14) 'STEP',rel_step

10      FORMAT(A,I4,A)
11      FORMAT(A,F7.3)
12      FORMAT(A,1P,E11.3,0P)
13      FORMAT(A)
14      FORMAT(A,I4)
15      FORMAT(A,I4,A,I2,A)
      ENDIF

      IF (iref.EQ.974) THEN
c...    974: Generic contour plot

        WRITE(char(23),10) 'OUTER TARGET DATA SHIFT =',
     .                    NINT(tarshift(IKHI)*1000.0),' mm'
        WRITE(char(25),11) 'Te      MULTIPLIER =',te_mult_i
        WRITE(char(26),11) 'Ti      MULTIPLIER =',ti_mult_i
        WRITE(char(27),12) 'Ne/Jsat MULTIPLIER =',n_mult_i

        WRITE(char(29),13) filename(1:LEN_TRIM(filename))
        IF (rel_step.NE.0) 
     .  WRITE(char(30),14) 'STEP',rel_step
      ENDIF

      IF (iref.EQ.966) THEN
c...    966: Divertor Th/OSM comparison

        WRITE(char(10),10) 'OUTER TARGET   DATA SHIFT =',
     .                    NINT(tarshift(IKHI)*1000.0),' mm'
        WRITE(char( 9),10) 'DTS            DATA SHIFT =',
     .                     INT(r1*1000.0),' mm'
        WRITE(char(6) ,11) 'Te      MULTIPLIER =',te_mult_i
        WRITE(char(5) ,11) 'Ti      MULTIPLIER =',ti_mult_i
        WRITE(char(4) ,12) 'Ne/Jsat MULTIPLIER =',n_mult_i

        WRITE(char(28),13) filename(1:LEN_TRIM(filename))
        IF (rel_step.NE.0) 
     .  WRITE(char(29),14) 'STEP',rel_step
        WRITE(char(30),13) job(38:72)
      ENDIF

      IF (iref.EQ.351.OR.iref.EQ.353.or.iref.eq.355.or.
     .    iref.EQ.978.OR.
     .    iref.EQ.970) THEN
c...    351: Core Th/OSM comparison
c       978: Combination plot
c       970: History plots

        WRITE(char(10),10) 'INNER TARGET DATA SHIFT =',
     .                    NINT(tarshift(IKLO)*1000.0),' mm'
        WRITE(char( 9),10) 'OUTER TARGET DATA SHIFT =',
     .                    NINT(tarshift(IKHI)*1000.0),' mm'
        WRITE(char(7) ,11) 'Te      MULTIPLIER =',te_mult_i
        WRITE(char(6) ,11) 'Ti      MULTIPLIER =',ti_mult_i
        WRITE(char(5) ,12) 'Ne/Jsat MULTIPLIER =',n_mult_i
         
        IF (iref.EQ.353.and.cre2d.gt.0) THEN
          WRITE(char(1),'(A)') '--- UEDGE'          
        ENDIF

        if (iref.eq.351.or.iref.eq.353.or.iref.eq.355) then 
c
c          Accept a maximum of 8 lines of comment for now
c
           ncom = min(ncomments,8)
           do i1 = 1,ncom
              char(19+i1) = extra_comments(i1)
           end do  
        endif

        WRITE(char(28),13) filename(1:LEN_TRIM(filename))
        IF (rel_step.NE.0) 
     .  WRITE(char(29),15) 'STEP',rel_step,' (',loadstep,')'
        WRITE(char(30),13) job(38:72)
      ENDIF

      RETURN
99    STOP
      END
c
c
c
c
c
      SUBROUTINE DrawAdditionalSurfaces(iopt)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comgra'
      INCLUDE 'pindata'
      INCLUDE 'colours'
      INCLUDE 'slcom'
      INCLUDE 'slout'

      INTEGER      i1,i2,i3,i4,i,j,in,iopt,cut1,nselect,select(100),l1,
     .             cell,cell2,fp1,index,ik1,ik2,ir,listi,id,ik,ir1,
     .             iks,ike
      LOGICAL      wrtnum,status,shiftthetav
      REAL         x,y,r,pvals(5,2),dr1,dz1,x1,y1,x2,y2,z1,z2,t1,t2,
     .             zval,thetav,listr(MAXNRS),listz(MAXNRS),fract,fracp,
     .             r1,r2,thetav1
C     IPP/01 - Krieger: fixed initialization in declaration statement
C     by adding appropiate data statement (picky SUN compiler)
      CHARACTER*36 NAME
      CHARACTER*99 graph6,cdum1
      CHARACTER*20 gaugename(110),celnum,segnum
      CHARACTER*3  cindex
c
      DATA NAME /' '/
c
c     INIT: 
c
      CALL LINCOL(ncols+1)
      CALL FILCOL(ncols+1)


      CALL PSPACE(map1x,map2x,map1y,map2y)
      CALL MAP   (cxmin,cxmax,cymin,cymax)


        call setup_col(ncols,5)

c...    I AM GOING TO MOVE THIS TO ANOTHER ROUTINE

        icol = 1

c        gaugename(001) = 'Midplane gauge'
c        gaugename(002) = 'Divertor gas box gauge'

        gaugename(101) = 'PBF1'
        gaugename(102) = 'PBF2'
        gaugename(103) = 'PBF3'
        gaugename(104) = 'PV1'
        gaugename(105) = 'PR2'
        gaugename(106) = 'VPLOWS'
        gaugename(107) = 'PCM105BAF'
        gaugename(108) = 'PCM240TOR'

        CALL FULL


c...  Set font size based on plot range:
      IF (iopt.NE.5.AND.iopt.NE.6) THEN
        dr1 = 0.0
        dz1 = 0.0
        wrtnum = .FALSE.    
      ELSEIF (cxmax-cxmin.LT.0.15) THEN
        CALL CTRMAG(8)
        dr1 = 0.014 * (cxmax-cxmin)
        dz1 = 0.010 * (cymax-cymin)
      ELSEIF (cxmax-cxmin.LT.0.65) THEN
        CALL CTRMAG(4)
        dr1 = 0.010 * (cxmax-cxmin)
        dz1 = 0.006 * (cymax-cymin)
      ELSE
        CALL CTRMAG(2)
        dr1 = 0.005 * (cxmax-cxmin)
        dz1 = 0.004 * (cymax-cymin)
      ENDIF

c...  Draw simulated pressure gauges:      
      IF (iopt.NE.7) THEN
        DO i = 1, eirnpgdat
          x = eirpgdat(i,2)
          y = eirpgdat(i,3)
          r = eirpgdat(i,5)
          WRITE(cindex,'(I3)') INT(eirpgdat(i,1))
c...      Lame way to draw a circle I know, but I can't get the CIRCLE routine
c         to work:
          DO J = 0, 355, 5
            x1 = x + r * COS(REAL(j  ) * PI / 180.0)
            y1 = y + r * SIN(REAL(j  ) * PI / 180.0)
            x2 = x + r * COS(REAL(j+5) * PI / 180.0)
            y2 = y + r * SIN(REAL(j+5) * PI / 180.0)
            CALL POSITN(x1,y1)
            CALL JOIN  (x2,y2)
c            CALL GRTRAC(pvals(1,1),pvals(1,2),2,name,'LINE',0)
          ENDDO
          CALL PLOTCS(x+r,y,' '//cindex//' '//
     .                    gaugename(INT(eirpgdat(i,1))))
        ENDDO
      ENDIF


c...  Check for 3D mode and a z-value specified in the OUT input file:
      zval = -99.0
      IF (asc_3dmode.EQ.1.OR.asc_3dmode.EQ.2) THEN
        READ(5,'(A80)',END=10) graph6
        IF   (graph6(8:11).EQ.'Zval'.OR.graph6(8:11).EQ.'ZVAL'.OR.
     .        graph6(8:11).EQ.'zval') THEN
          READ(graph6,*) cdum1,zval
          WRITE(char(28),'(A,F10.3,A)') 'Zval  = ',zval,' m'
        ELSE
          BACKSPACE 5
        ENDIF
10      CONTINUE
      ENDIF

c...  Select only certain cells to be plotted:
      nselect = 0
      READ(5,'(A80)',END=15) graph6
      IF   (graph6(8:11).EQ.'Cell'.OR.graph6(8:11).EQ.'CELL'.OR.
     .      graph6(8:11).EQ.'cell') THEN
        READ(graph6,*) cdum1,nselect,(select(i1),i1=1,nselect)
      ELSE
        BACKSPACE 5
      ENDIF
15    CONTINUE

c...  Draw vacuum grid:      
      DO i2 = 1, asc_ncell
        IF (dr1.NE.0.0) wrtnum = .TRUE.

c        CYCLE

c...    Allow for plotting of selected cells only (selection via input file):
        status = .TRUE.
        IF (nselect.GT.0) THEN
          status = .FALSE.
          wrtnum = .TRUE.
        ENDIF
        DO i1 = 1, nselect
          cell = MOD(select(i1),asc_ncell)            
          IF (cell.EQ.0) cell = asc_ncell
          IF (i2.EQ.cell) THEN
            status = .TRUE.
            cell2  = select(i1)
          ENDIF
        ENDDO
        IF (.NOT.status) CYCLE

c...    Don't plot cells that do not span the z-value specified
c       in the OUT input file is 3D mode is present:
        IF (asc_3Dmode.EQ.1.AND.zval.NE.-99.0.AND.
     .      .NOT.(asc_zmin3D(i2).LT.zval.AND.
     .            asc_zmax3D(i2).GT.zval)) CYCLE

        DO i3 = 1, asc_nvp(i2)
          IF (asc_link(i3,i2).GE.0) THEN
            x1 = asc_rvp(i3,i2)
            y1 = asc_zvp(i3,i2)
            x2 = asc_rvp(i3+asc_nvp(i2),i2)
            y2 = asc_zvp(i3+asc_nvp(i2),i2)
            CALL POSITN(x1,y1)
            CALL JOIN  (x2,y2)
          ENDIF
          IF (i3.EQ.1) THEN
            i4 = asc_nvp(i2)
          ELSE
            i4 = i3 - 1
          ENDIF
          IF (wrtnum.AND.asc_link(i3,i2).GE.0.AND.
     .                   asc_link(i4,i2).GE.0.AND.
     .        asc_rvp(i3,i2).GT.cxmin.AND.asc_rvp(i3,i2).LT.cxmax.AND.
     .        asc_zvp(i3,i2).GT.cymin.AND.asc_zvp(i3,i2).LT.cymax)
     .      THEN
            wrtnum = .FALSE.

c...        3D mode:
            cut1 = 1
            IF (asc_3dmode.EQ.2.AND.zval.NE.-99.0) THEN
              DO cut1 = 1, ascncut
                IF (asc_zmin3D(cut1).LE.zval.AND.
     .              asc_zmax3D(cut1).GT.zval) EXIT
              ENDDO
              IF (cut1.EQ.ascncut+1) THEN
                WRITE(0,*) 'DrawAdditionalSurfaces: Zval outside '//
     .                      'vacuum grid'
                cut1 = 1
              ENDIF
            ENDIF

            IF (nselect.GT.0) THEN
              CALL CTRMAG(10)
              WRITE(celnum,'(I6)') cell2
              l1 = -1
              DO i4 = 1, LEN_TRIM(celnum)
                IF (celnum(i4:i4).EQ.' ') l1 = i4
              ENDDO
              IF (l1.EQ.-1) l1=1
              CALL PLOTST(asc_rvp(i3,i2)+dr1,asc_zvp(i3,i2)+dz1,
     .                    celnum(l1:6))
            ELSE
              WRITE(celnum,'(I5,17X)') i2 + (cut1 - 1) * asc_ncell
              IF     (i3.EQ.1) THEN
                CALL PCSCEN(asc_rvp(i3,i2)-1.7*dr1,
     .                      asc_zvp(i3,i2)+    dz1,
     .                      celnum(1:5))
              ELSEIF (i3.EQ.2) THEN
                CALL PCSCEN(asc_rvp(i3,i2)+    dr1,
     .                      asc_zvp(i3,i2)+    dz1,
     .                      celnum(1:5))
              ELSEIF (i3.EQ.3) THEN
                CALL PCSCEN(asc_rvp(i3,i2)+    dr1,
     .                      asc_zvp(i3,i2)-1.8*dz1,
     .                    celnum(1:5))
              ELSEIF (i3.EQ.4) THEN
                CALL PCSCEN(asc_rvp(i3,i2)-1.7*dr1,
     .                      asc_zvp(i3,i2)-1.8*dz1,
     .                      celnum(1:5))
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDDO


c...  Draw additional surfaces:
      IF (iopt.EQ.6) THEN
        fp1=98
        OPEN(UNIT=fp1,FILE='dump.dat',ACCESS='SEQUENTIAL',
     .       STATUS='OLD',ERR=17)
        status = .TRUE.
        DO WHILE (status)
          READ(fp1,*,END=17) index,x1,y1,z1,x2,y2,z2
          IF (index.LT.ascncut) THEN
            CALL LINCOL(ncols+3)
          ELSE
            CALL LINCOL(ncols+4)
          ENDIF
          CALL POSITN (x1,y1)
          CALL JOIN   (x2,y2)
        ENDDO
17      CONTINUE
        CLOSE(fp1)
      ENDIF


      IF (iopt.EQ.2.OR.iopt.EQ.6.OR.iopt.EQ.8.OR.iopt.EQ.5) THEN
c...    SOL24 fronts:
        CALL LINCOL(ncols+3)
        DO i1 = 1, osm_nfnt-1
          CALL POSITN(osm_ionfnt(i1  ,2),osm_ionfnt(i1  ,3))
          CALL JOIN  (osm_ionfnt(i1+1,2),osm_ionfnt(i1+1,3))
        ENDDO
c...    SOL28 fronts:
        CALL LINCOL(ncols+1)

        DO i1 = 1, osmns28-1
          IF (osms28(i1,1).NE.osms28(i1+1,1)) CYCLE
          CALL POSITN(osms28(i1  ,5),osms28(i1  ,6))
          CALL JOIN  (osms28(i1+1,5),osms28(i1+1,6))
c          CALL POSITN(osms28(i1  ,3),osms28(i1  ,4))
c          CALL JOIN  (osms28(i1+1,3),osms28(i1+1,4))
        ENDDO
      ENDIF

c...  Draw reflections:
      CALL LINCOL(ncols+4)

      DO i1 = 1, nshow
c        IF (i1.GE.362.AND.i1.LE.1000) CYCLE
c        IF (i1.GT.1806) CYCLE
c        IF (i1.GE.22.AND.i1.LE.1000) CYCLE
 
        x1 = rshow(i1)
        y1 = zshow(i1)
        x2 = COS(ashow(i1)) * 10.0 + x1
        y2 = SIN(ashow(i1)) * 10.0 + y1

        WRITE(0,*) 'hsdf:',i1,x1,x2,y1,y2

        CALL POSITN(x1,y1)
        CALL JOIN  (x2,y2)
      ENDDO        
c...  Only plot views once:
      nshow = 0 

c...  Draw target probes:
      IF (iopt.GE.5) THEN
        CALL LINCOL(ncols+1)
        DO i1 = 1, prb_num(OFMP)
          x = prb_r(i1,OFMP)
          y = prb_z(i1,OFMP)
          r = 0.004
          CALL POSITN(x+r,y)
          DO J = 20, 340, 20
            x1 = x + r * COS(REAL(j+20) * PI / 180.0)
            y1 = y + r * SIN(REAL(j+20) * PI / 180.0)
            CALL JOIN(x1,y1)
          ENDDO
        ENDDO
        DO i1 = 1, prb_num(IFMP)
          x = prb_r(i1,IFMP)
          y = prb_z(i1,IFMP)
          r = 0.004
          CALL POSITN(x+r,y)
          DO J = 20, 340, 20
            x1 = x + r * COS(REAL(j+20) * PI / 180.0)
            y1 = y + r * SIN(REAL(j+20) * PI / 180.0)
            CALL JOIN(x1,y1)
          ENDDO
        ENDDO
      ENDIF

c...  Draw puffing surfaces:
      CALL LINCOL(ncols+1)
      CALL THICK (2)
      DO i1 = 1, 0  ! eirnpuff
        in = NINT(eirpuff(i1,3))
        t1 = eirpuff(i1,5)
        t2 = eirpuff(i1,6)
        x1 = rvesm(in,1) + t1 * (rvesm(in,2) - rvesm(in,1))
        y1 = zvesm(in,1) + t1 * (zvesm(in,2) - zvesm(in,1))
        x2 = rvesm(in,1) + t2 * (rvesm(in,2) - rvesm(in,1))
        y2 = zvesm(in,1) + t2 * (zvesm(in,2) - zvesm(in,1))
        CALL POSITN(x1,y1)
        CALL JOIN  (x2,y2)
      ENDDO

c...  Draw pumping wall surfaces:
      CALL LINCOL(ncols+3)
      CALL THICK (2)
      DO i1 = 1, eirnspdat
        IF (eirspdat(i1,1).EQ.2.0.AND.(eirspdat(i1,3).EQ.2.0.OR.
     .                                 eirspdat(i1,8).LT.1.0)) THEN
          in = NINT(eirspdat(i1,2))
          x1 = rvesm(in,1)
          y1 = zvesm(in,1)
          x2 = rvesm(in,2)
          y2 = zvesm(in,2)
          CALL POSITN(x1,y1)
          CALL JOIN  (x2,y2)        
        ENDIF
      ENDDO

c...  Draw indecies on wall surfaces:
      CALL CTRMAG(8)
      IF (iopt.EQ.5.OR.iopt.EQ.6) THEN
        CALL LINCOL(ncols+4)
        IF (nvesm.NE.0) THEN
          DO i1 = 1, nvesm+nvesp
            IF (jvesm(i1).NE.1.AND.jvesm(i1).NE.4) THEN
              x1 = 0.5 * (rvesm(i1,1) + rvesm(i1,2))
              y1 = 0.5 * (zvesm(i1,1) + zvesm(i1,2))
              IF (x1.GT.cxmin.AND.x1.LT.cxmax.AND.
     .            y1.GT.cymin.AND.y1.LT.cymax) THEN
                WRITE(segnum(1:3),'(I3)') i1
                CALL PCSCEN(x1,y1,segnum(1:3))
              ENDIF
            ENDIF
          ENDDO
        ELSEIF (nves.NE.0) THEN
          DO i1 = 1, nves-1
            x1 = 0.5 * (rves(i1) + rves(i1+1))
            y1 = 0.5 * (zves(i1) + zves(i1+1))
            IF (x1.GT.cxmin.AND.x1.LT.cxmax.AND.
     .          y1.GT.cymin.AND.y1.LT.cymax) THEN
              WRITE(segnum(1:3),'(I3)') i1
              CALL PCSCEN(x1,y1,segnum(1:3))
            ENDIF
          ENDDO
        ENDIF
      ENDIF

c...  Draw time-to-ionisation regions:
      IF (iopt.EQ.7) THEN
        CALL LINCOL(ncols+4)
        CALL THICK (1)
        CALL CTRMAG(10)
        DO i1 = 1, eirniontime
          x1 = eiriontime(i1,1)
          y1 = eiriontime(i1,3)
          x2 = eiriontime(i1,2)
          y2 = eiriontime(i1,4)
          CALL POSITN(x1,y1)
          CALL JOIN  (x2,y1)
          CALL POSITN(x2,y1)
          CALL JOIN  (x2,y2)
          CALL POSITN(x2,y2)
          CALL JOIN  (x1,y2)
          CALL POSITN(x1,y2)
          CALL JOIN  (x1,y1)
c...      Only plot average time to ionisation if the region
c         appears on the plot:
          IF (x1.GT.cxmin.AND.x1.LT.cxmax.AND.
     .        y1.GT.cymin.AND.y1.LT.cymax) THEN
            WRITE(cdum1,'(I2)') i1
c            WRITE(cdum1,'(A,1P,E9.2,0P,A)') 'Tavg =',eiriontime(i1,19),
c     .                                       ' s'
            CALL PLOTCS(0.9*x1+0.1*x2,0.5*y1+0.5*y2,
     .                  cdum1(1:LEN_TRIM(cdum1)))
          ENDIF
        ENDDO        
      ENDIF

c...  Draw THETAG location:
      IF (.FALSE.) THEN


c...    Draw poloidal surfaces:
        CALL LINCOL(ncols+2)
        CALL THICK(1)
        DO ir = 2, nrs
          IF (idring(ir).EQ.BOUNDARY) CYCLE
          DO ik = 1, nks(ir)
            id = korpg(ik,ir)
            CALL POSITN(rvertp(2,id),zvertp(2,id))
            CALL JOIN  (rvertp(3,id),zvertp(3,id))
          ENDDO
        ENDDO


c...    DTHETG is not passed in the main DIVIMP/OUT transfer (.RAW) file at
c       the moment, so loading the supplimental geometry data file is required:
        IF (dthetg.EQ.0.0) CALL ER('DrawAdditionalSurfaces',
     .                             'DTHETG=0.0',*99)

c...    Assemble list of points:

c...    Assign THETAV:
        ir1 = irtrap + 1
c        ir1 = 2
        iks = 1
        ike = nks(ir1)
        IF (ir1.LT.irsep) ike = ike - 1

        DO thetav1 = thetag(iks,ir1), 
     .               thetag(ike,ir1), 2.5
c     .               thetag(ike,ir1), 5.0

          ir = ir1 
          listi = 0
          thetav = thetav1
          ik2 = 0
          shiftthetav = .FALSE.

c          WRITE(0,*) 'THETAV',THETAV

          DO WHILE (ir.NE.irwall.AND.ir.GT.irtrap) 
c          DO WHILE (ir.NE.irwall) 

c...        Shift THETAV because the intersection is in the high index
c           PFZ:
            IF (shiftthetav) THEN
              thetav = thetav + dthetg
              shiftthetav = .FALSE.
            ENDIF

            ik2 = 0
            DO ik1 = 1, nks(ir)
              IF (thetav.GE.thetag(ik1,ir)) ik2 = ik1
            ENDDO

c...        Assume that THETAV is not valid for this ring:
            IF (.NOT.(ik2.EQ.0.OR.ik2.EQ.nks(ir))) THEN

c...          Find midpoint and adjust IK2 and necessary:
              fract = (thetav - thetag(ik2,ir)) / 
     .                (thetag(ik2+1,ir) - thetag(ik2,ir))
              fracp = (kpb(ik2  ,ir) - kps(ik2,ir)) /
     .                (kps(ik2+1,ir) - kps(ik2,ir))
              IF (fract.GT.fracp) ik2 = ik2 + 1

c...          Increase THETAV but DTHETG on the next iteration:
              IF (ir.EQ.nrs.AND.ik2.GE.ikti2(ir)) shiftthetav = .TRUE.

c              WRITE(0,*) '-->',ik2,ir,thetav   

              fract = (thetav - thetag(ik2,ir)) / 
     .                (thetag(ik2+1,ir) - thetag(ik2,ir))
              fracp = (kpb(ik2  ,ir) - kps(ik2,ir)) /
     .                (kps(ik2+1,ir) - kps(ik2,ir))

              IF (fract.LT.fracp) THEN
                fract = fract / fracp
                id = korpg(ik2,ir)
                r1 = rs(ik2,ir) 
                z1 = zs(ik2,ir) 
                r2 = 0.5 * (rvertp(3,id) + rvertp(4,id))
                z2 = 0.5 * (zvertp(3,id) + zvertp(4,id))
              ELSE
                fract = (fract - fracp) / (1.0 - fracp)
                id = korpg(ik2+1,ir)
                r1 = 0.5 * (rvertp(1,id) + rvertp(2,id))
                z1 = 0.5 * (zvertp(1,id) + zvertp(2,id))
                r2 = rs(ik2+1,ir)
                z2 = zs(ik2+1,ir)
              ENDIF

              listi = listi + 1
              listr(listi) = r1 + fract * (r2 - r1)
              listz(listi) = z1 + fract * (z2 - z1)
            ENDIF

            IF (ik2.EQ.0) ik2 = nks(ir) / 2
            ir = irouts(ik2,ir)

          ENDDO
      
          IF (listi.NE.0) THEN
            CALL LINCOL(ncols+1)
            DO i1 = 1, listi-1
c              CALL POSITN(listr(i1  ),listz(i1  ))
c              CALL JOIN  (listr(i1+1),listz(i1+1))
            ENDDO
          ENDIF



        ENDDO

      ENDIF


c...  Reset drawing settings:
c     IPP/09 - Krieger - also reset color settings
c     which is not possible right now because n_cols, col_opt not avail.
c     call setup_col(n_cols,col_opt)
      CALL FULL
      CALL LINCOL(1)
      CALL THICK(1)

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
        IF (char(i).NE.' ') CALL PLOTST(1.00,0.550-(j-1)*0.02,char(i))
      ENDDO

      RETURN
99    STOP
      END

c
c ======================================================================
c
c
c
c
c
      SUBROUTINE AsgnThomsonData(nraw,raw,ind,gxdata,gydata,gndata,
     .                           MAXNPS,MAXTYP,MAXTDAT,MAXCOLS)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'
      INCLUDE 'slout'

      INTEGER nraw,ind,MAXNPS,MAXTYP,MAXTDAT,MAXCOLS
      REAL    raw(MAXTDAT,MAXCOLS)

      INTEGER gndata(MAXNRS,MAXTYP)
      REAL    gydata(MAXNPS,MAXNRS,MAXTYP),gxdata(MAXNPS,MAXNRS,MAXTYP),
     .        gwdata(MAXNPS,MAXNRS)

      INTEGER CalcPoint

      INTEGER ik,ik1,ir,in,id,num,code,idum1
      LOGICAL incik,decik
      REAL    xdat(MAXNPS),ydat(MAXNPS)
      REAL*8  ar,az,br,bz,cr,cz,t
c
      logical print_warning
      data print_warning /.true./
      save print_warning
c
c     jdemod - variables to support new assignment option
c
      real s_approx,cross
      integer thom_asgn_opt
c
      thom_asgn_opt = 0
c
      if     (thom_asgn_opt.eq.0.and.print_warning) then 
         write(0,*) 'THOMSON ASSIGNMENT OPTION SET TO 0 -'//
     >              ' ORIGINAL MAPPING CODE'
         write(6,*) 'THOMSON ASSIGNMENT OPTION SET TO 0 -'//
     >              ' ORIGINAL MAPPING CODE'
         print_warning = .false.
      elseif (thom_asgn_opt.eq.1.and.print_warning) then 
         write(0,*) 'THOMSON ASSIGNMENT OPTION SET TO 1 -'//
     >              ' USING GETSCROSS_APPROX'
         write(6,*) 'THOMSON ASSIGNMENT OPTION SET TO 1 -'//
     >              ' USING GETSCROSS_APPROX'
         print_warning = .false.
      endif
c
      CALL IZero(gndata(1,ind),MAXNRS)

c
c     jdemod - note - the cell containing the data point is found in the load_thomson_data routine
c                     using gridpos and stored in ik=raw(in,3) and ir=raw(in,4)
c                   - if the intersection method originally implemented fails then I think that
c                     the improved getscross_approx routine should be used to estimate the S coordinate
c                     of the data point. 
c                   - the reason for this is that there are alot of cases where the code used below
c                     will fail - particularly with a thomson point near a convex joining of two cells
c                     - the perpendicular intersection point will lie on neither of the axes of the two
c                     adjacent cells. 
c                   - on the other hand - gridpos has already determined the actual cell occupied by the data
c                     point and an appropriate estimate of S is possible. 
c


      DO ir = 2, nrs

        num = 0
        ydat = 0.0
        xdat = 0.0


        DO in = 1, nraw

          incik = .FALSE.
          decik = .FALSE.

          IF (INT(raw(in,4)).EQ.ir) THEN

c
c          jdemod - thom_asgn_opt 0 - original 
c
           if (thom_asgn_opt.eq.0) then 


            ik = INT(raw(in,3))

            ik1 = ik

10          IF (incik.AND.decik) ik = ik1

            id = korpg(ik,ir)

c            IF (incik.OR.decik) WRITE(0,*) 'ATD: IK= ',ik,in,t,
c     .            incik,decik,ik1

            ar = DBLE(0.5 * (rvertp(1,id) + rvertp(2,id)))
            az = DBLE(0.5 * (zvertp(1,id) + zvertp(2,id)))
            br = DBLE(0.5 * (rvertp(3,id) + rvertp(4,id)))
            bz = DBLE(0.5 * (zvertp(3,id) + zvertp(4,id)))
            cr = DBLE(raw(in,1))
            cz = DBLE(raw(in,2))

            idum1 = CalcPoint(ar,az,br,bz,cr,cz,t)

            IF (ir.EQ.15) THEN
c              WRITE(6,*) 'MAP>',ik,ir,xdat(num),ydat(num)
c              WRITE(6,*) 'MAP>',ik,ir,t,raw(in,4+ind)
            ENDIF


c
c           jdemod - if perpendicular intersection of thomson point with 
c                    the cell axis lies on the line segment. 
c
            IF    ((t.GE.0.0.AND.t.LE.1.0).OR.(incik.AND.decik)) THEN
              num= num + 1
              xdat(num) = ksb(ik-1,ir)+REAL(t)*(ksb(ik,ir)-ksb(ik-1,ir)) 
c...          Nudge is required so that LoadArray creates an entry for points 
c             with the same (well, not any more) xdat value
              xdat(num) = xdat(num) * (1.0 + REAL(num)*1.0E-6) 
              ydat(num) = raw(in,4+ind)

C
C           jdemod - looks like debugging code ...
C
            IF (ir.EQ.15) THEN
              WRITE(6,*) 'MAP> 15 :',ik,ir,xdat(num),ydat(num)
c              WRITE(6,*) 'MAP>',ik,ir,t,raw(in,4+ind)
            ENDIF

c
c           jdemod - intersection is before start of line segment
c
            ELSEIF (t.LT.0.0) THEN
c
c...          Data point intersects a cell with lower IK:
              IF (ik.GT.1) THEN
                ik = ik - 1
                decik = .TRUE.
                GOTO 10
              ELSE
                IF (t.GT.-0.2) THEN
                  num= num + 1
                  xdat(num) = 0.0
                  ydat(num) = raw(in,4+ind)
                ELSE

                  write(6,*) 'AsgnThomsonData: Error mapping:',cr,cz
                  write(0,*) 'AsgnThomsonData: Error mapping:',cr,cz
                  CALL ER('AsgnThomsonData','Cannot position data',*99)
                ENDIF
              ENDIF
            ELSEIF (t.GT.1.0) THEN
c
c...          Data point intersects a cell with higher IK:
              IF (ik.LT.nks(ir)) THEN
                ik = ik + 1
                incik = .TRUE.
                GOTO 10
              ELSE
                IF (t.LT.1.2) THEN
                  num= num + 1
                  xdat(num) = ksmaxs(ir)
                  ydat(num) = raw(in,4+ind)
                ELSE
                  write(6,*) 'AsgnThomsonData: Error mapping:',cr,cz
                  CALL ER('AsgnThomsonData','Cannot position data',*99)
                ENDIF
              ENDIF
            ENDIF


c
c         jdemod - adding an alternate way to determine the S values
c                  mapped for thomson points. In general shouldn't be too 
c                  different
c
c          NOTE: offgrid thomson data is filtered out when it is loaded 
c
c
           elseif (thom_asgn_opt.eq.1) then 
c
c             Find approximate S,CROSS location
c
              call getscross_approx(raw(in,1),raw(in,2),s_approx,cross,
     >                          int(raw(in,3)),int(raw(in,4)))

c              write(6,'(a,2i6,2(1x,f6.1),4(1x,f12.5))') 'TH:',in,num,
c     >                raw(in,3),raw(in,4),raw(in,1),raw(in,2),
c     >                s_approx,cross
c
             num=num+1
c
c...          Nudge is required so that LoadArray creates an entry for points 
c             with the same (well, not any more) xdat value
c              xdat(num) = xdat(num) * (1.0 + REAL(num)*1.0E-6) 
c       
c            Getting messages from maparray about overwriting data - so nudge xdat in case
c            I am getting repeating values of s_approx - which can happen for points reasonably
c            close to each other. 
c
             xdat(num) = s_approx  * (1.0 + REAL(num)*1.0E-7)
c             xdat(num) = s_approx 
             ydat(num) = raw(in,4+ind)
c
c           jdemod - endif for thom_asgn_opt
c
           endif


          ENDIF


c          IF (num.GT.0) THEN
c            gndata(ir,ind) = 0

c            IF (ir.EQ.15) THEN
c              WRITE(6,*) 'MAP>',num,ydat(1:num)
c            ENDIF

c            CALL LoadArray(gxdata(1,ir,ind),gndata(ir,ind),xdat,1,num)
c            CALL MapArray (gxdata(1,ir,ind),gydata(1,ir,ind),
c     .                     1,gndata(ir,ind),
c     .                     xdat,ydat,1,num)
c          ENDIF

        ENDDO



          IF (num.GT.0) THEN
            gndata(ir,ind) = 0


            CALL LoadArray(gxdata(1,ir,ind),gndata(ir,ind),xdat,1,num)
            CALL MapArray (gxdata(1,ir,ind),gydata(1,ir,ind),
     .                     1,gndata(ir,ind),
     .                     xdat,ydat,1,num)

c            IF (ir.EQ.15) THEN
c              WRITE(6,'(A,2I6,10F10.2)') 
c     .          'MAP>',num,gndata(ir,ind),ydat(1:num)
c
c             DO ik = 1, gndata(ir,ind)  ! LEFT OFF
c              WRITE(6,'(A,3I6,2(1x,g14.4))')
c     .     ' TMAP 15 : ',ind,ik,ir,gxdata(ik,ir,ind),gydata(ik,ir,ind)
c              ENDDO
c
c            ENDIF

          ENDIF

      ENDDO



      DO ir = irsep, nrs
        DO ik = 1, gndata(ir,ind)
          WRITE(6,'(A,3I6,2(1x,G14.4))')
     .      ' TMAP    : ',ind,ik,ir,gxdata(ik,ir,ind),gydata(ik,ir,ind)
        ENDDO
      ENDDO


c      STOP

c...     PCS wants the scatter plotted -- best to generalize the x data and
c        then slap the scattered data wherever the cips fall -- then
c        slap on the average value as well  -- need new x-axis option that
c        allows xdata to be specified here:

c        CALL LoadArray(tdist,jd,rho    (1,CELL),irsep,irwall-1     )
c        CALL LoadArray(tdist,jd,prb_rho(1,FSP1),1    ,prb_num(FSP1))
c
c        CALL MapArray(tdist          ,dvals (1,1)   ,1,jd,
c     .                prb_rho(1,FSP1),prb_te(1,FSP1),1,prb_num(FSP1))
c        CALL MapArray(tdist          ,dvals (1,3)   ,1,jd,
c     .                prb_rho(1,FSP1),prb_ne(1,FSP1),1,prb_num(FSP1))


      RETURN
99    STOP
      END
c
c ======================================================================
c
      SUBROUTINE Plot966(nplts,ringnos,graph,nplots,ref,title,iopt,
     .                   ngrm,pltmins,pltmaxs,pltfact,iplot,job)
      IMPLICIT none


      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'cedge2d'
      INCLUDE 'colours'
      INCLUDE 'slcom'
      INCLUDE 'slout'

      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost

      REAL CalcPressure,GetEAD,GetCs

      INTEGER ik,ir

      integer iplot
      REAL    NP,NS,NT,FP,FT,TIME,TIME1,ZA02AS,ASPECT,POINT
      INTEGER IOPT,J,NIZS,ITER,NITERS,IZ,JZ,ISMOTH,JK,IN,inc
      integer  nplts,ringnos(maxplts),ip,pnks(maxplts)
      CHARACTER TITLE*(*),JOB*72,GRAPH*80,GRAPH1*80

      character*36 pltlabs(maxplts)

      CHARACTER*36 XLAB,YLAB,XPOINT
      CHARACTER*36 REF,NVIEW,PLANE,ANLY,TABLE,ZLABS(-2:MAXIZS+1)
      CHARACTER*36 NAME,ELABS(MAXNGS),PLABS(-2:MAXPLRP),KLAB

      CHARACTER*128 REF2

      INTEGER IX,IY,IG,NGS,IERR,IXMIN,IXMAX,NPLOTS,M,ID,JR

      INTEGER ii1,ii2,midnks,i1,i2,ret,plt_opt1,osm_cfs,npts,ntime,
     .        in1,in2,xtype,ytype,btype,id1,id2
      integer  sctype,ngrm
      real pltmins(maxplts),pltmaxs(maxplts)
      real pltmin,pltmax
      real pltfact

      REAL tauii,pii
      LOGICAL status,cont
      INTEGER nenum,tenum,opt_const,plot_mode(30),array,iter1,iter2,
     .        xaxis,ring,mode,inorm(MAXNGS)

      REAL          te1,ti1,ne1,te2,ti2,ne2,norm,
     .              rdum1,rdum2,rdum3,rdum4,rdum5,rdum6,rdum7,
     .              radum1(MAXNRS),radum2(MAXNRS),radum3(MAXNRS)
      REAL    nemin,nestep,temin,temax,neTe,frac1,
     .        max1,max2,ynorm(MAXNGS)


      REAL CalcDist
      INTEGER iks,ike,ik1
      REAL    deltas,dist(MAXNKS,MAXNRS),sum,magcfp,sumn,sump,factp,
     .        factn,pm,pl,xdata(0:MAXNKS,MAXNRS),ydata(0:MAXNKS,MAXNRS),
     .        xsh,ysh,targdat(MAXNDS,5)

      INTEGER    MAXTDAT     ,MAXCOLS   ,MAXNPS    ,MAXTYP
      PARAMETER (MAXTDAT=500,MAXCOLS=10,MAXNPS=1000,MAXTYP=10)
c      PARAMETER (MAXTDAT=4000,MAXCOLS=10,MAXNPS=1000,MAXTYP=10)

      INTEGER thnnraw,ind
      REAL    avgdat(MAXNKS,MAXNRS),thnraw(MAXTDAT,MAXCOLS),cnt

      INTEGER    MAXTHDAT
      PARAMETER (MAXTHDAT=15)
c      PARAMETER (MAXTHDAT=15)
c      INTEGER gndata(       MAXNRS,MAXTYP)
c      REAL    gydata(MAXNPS,MAXNRS,MAXTYP),
c     .        gxdata(MAXNPS,MAXNRS,MAXTYP)


c      INTEGER gndata(       MAXNRS,MAXTYP,MAXTHDAT),ngdata
c      REAL    gydata(MAXNPS,MAXNRS,MAXTYP,MAXTHDAT),
c     .        gxdata(MAXNPS,MAXNRS,MAXTYP,MAXTHDAT)


      INTEGER ngdata
      INTEGER, ALLOCATABLE :: gndata(:,:,:)
      REAL   , ALLOCATABLE :: gydata(:,:,:,:),gxdata(:,:,:,:)
     .        

      INTEGER idum1,idum2,idum3,ncase,casendata(MAXNRS,5)
      CHARACTER*128 filename

      real mvals(maxnps,maxplts,maxngs),hue,saturation,value
      real mouts (maxnps,maxplts),mwids (maxnps,maxplts)
c      real mvals(maxnks+2,maxplts,maxngs)
c      real mouts (maxnks+2,maxplts),mwids (maxnks+2,maxplts)

      LOGICAL langdat,plotuedge,convert
      CHARACTER cdum1*128
      INTEGER   step
      CHARACTER cname*256,cmnd*256,dummy*256

      REAL xrange1,xrange2,sbnd1,sbnd2

      INTEGER v1,cell,cell2,shift,cut1,i3,optscale,line,set
      REAL    xcen,ycen,xmin,xmax,ymin,ymax,fact,frac,old,dmax,gmax,
     .        xval,yval,zval

      REAL,ALLOCATABLE :: osmtmp(:,:),casedata(:,:,:,:),casetdata(:,:,:)
c
c jdemod
c
c Variables for averaged DTS plots - ytype options 40 and 42  
c
      integer ndatasets
      real avdata(maxnks,maxnrs),ncounts(maxnks,maxnrs)
      logical plotav
c 
c jdemod 
c
c     Array for extra comments 
c
      character*20 extra_comments(1)
c
      REAL      scale
      CHARACTER fname*1024

      REAL nD2p,ti,ne,sigmav1,sigmav2,sigmav

      langdat = .FALSE.
      plotuedge = .FALSE.
      ngdata = 0
      optscale = 0

      ALLOCATE(osmtmp(MAXNKS,MAXNRS))
      ALLOCATE(casedata(MAXNKS,MAXNRS,0:7,5))
      ALLOCATE(casetdata(MAXNDS,4,5))

      ALLOCATE (gndata(       MAXNRS,MAXTYP,MAXTHDAT))
      ALLOCATE (gydata(MAXNPS,MAXNRS,MAXTYP,MAXTHDAT))
      ALLOCATE (gxdata(MAXNPS,MAXNRS,MAXTYP,MAXTHDAT))

      CALL RSet(avgdat,MAXNKS*MAXNRS,LO)

      CALL CalculateRho 
c
c
c
       iopt_ghost = 1
c        grm_opt    = 1
       slopt      = 1

       CALL RZero(mvals,MAXNPS*MAXPLTS*MAXNGS)
       CALL RZero(mouts,MAXNPS*MAXPLTS)
       CALL RZero(mwids,MAXNPS*MAXPLTS)

       CALL RZero(xdata(0,1),(MAXNKS+1)*MAXNRS)
       CALL RZero(ydata(0,1),(MAXNKS+1)*MAXNRS)

       DO ip = 1, nplts
         WRITE (pltlabs(ip),'(i4)') ringnos(ip)
       ENDDO

       ref2 = ref

       xsh = tarshift(IKHI)

       READ(graph(14:15),*) xtype
       READ(graph(17:18),*) ytype

       WRITE(6,*) '966: xtype.ytype=',xtype,ytype

c ...change to an array...
       IF     (xtype.EQ.1.OR.xtype.EQ.6) THEN
c...     1 - standard
c        6 - MapArray called
c        9 - MapArray called, s inverted to outer target is s=0
         XLAB = '   s (m) '
       ELSEIF (xtype.EQ.9) THEN
c...     9 - MapArray called, s inverted to outer target is s=0
         XLAB = '   s (m) '
         DO ir = 1, nrs
           DO ik = 1, nks(ir)
             kss(ik,ir) = ksmaxs(ir) - kss(ik,ir)
           ENDDO
         ENDDO
       ELSEIF (xtype.EQ.2) THEN
         XLAB = '   theta  '
       ELSEIF (xtype.EQ.3.OR.xtype.EQ.8) THEN
         XLAB = '   s-pol (m)'
       ELSEIF (xtype.EQ.4) THEN
         XLAB = '   s / smax'
       ELSEIF (xtype.EQ.5) THEN
         XLAB = '   Te (eV)'
       ELSEIF (xtype.EQ.7) THEN
         XLAB = '   Z (m)'
       ELSE
         CALL ER('966','Unrecognized XTYPE specified',*99)
       ENDIF
c

c       PLOTTYPE(1) = 2
c       PLOTTYPE(2) = 3
c       PLOTTYPE(3) = 4
c       PLOTTYPE(4) = 5
c       PLOTTYPE(5) = 6
c       PLOTTYPE(6) = 7

c       ngs = 0
c...   Spagetti - for plotting results from more than one case:
c 5     CONTINUE

       IF     (ytype.EQ.1) THEN
         ylab     = 'plasma'
         elabs(1) = '    n'
         elabs(2) = '    Te'
         elabs(3) = '    Ti'
         elabs(4) = '    v'
         elabs(5) = '    ga'
         ngs      = 1
c         ngs      = 5
         btype    = 1
       ELSEIF (ytype.EQ.2) THEN
         ylab     = 'T (eV)'
         elabs(1) = '    Te'
         elabs(2) = '    Ti'
         ngs      = 2
         btype    = 1
       ELSEIF (ytype.EQ.3) THEN
         ylab     = ' '
         elabs(1) = '    Te'
         elabs(2) = '    Ti'
         elabs(3) = '    Ga'
         ngs      = 3
         elabs(1) = '    R'
         elabs(2) = '    Z'
         ngs      = 2
         btype    = 2
       ELSEIF (ytype.EQ.4) THEN
         ylab     = ' '
         elabs(1) = '    Te'
         elabs(2) = '    ne'
         elabs(3) = '    Ha'
         ngs      = 3
         btype    = 2
       ELSEIF (ytype.EQ.5) THEN
         ylab     = ' '
         elabs(1) = '    Te'
         elabs(2) = '    ne'
         elabs(3) = '    P'
         ngs      = 3
         btype    = 1
       ELSEIF (ytype.EQ.6) THEN
         ylab     = 'pin'
         elabs(1) = '    rec'
         elabs(2) = '    ion'
         ngs      = 2
         btype    = 2
       ELSEIF (ytype.EQ.7) THEN
c...     Total plasma pressure:
         mvals = LO
         ylab     = 'p (eV m-3)'
         elabs(1) = '    p'
         plottype(3) = -1
         plottype(4) = -2
         slopt2    =  2
         ngs      = 1
         btype    = 2
         DO ir = 1, nrs
           IF (idring(ir).EQ.BOUNDARY) CYCLE
           DO ik = 1, nks(ir)
             osmtmp(ik,ir) = CalcPressure(knbs(ik,ir),ktebs(ik,ir),
     .                                    ktibs(ik,ir),kvhs(ik,ir))
           ENDDO
         ENDDO
c         IF (pinqi (1,1).EQ.LO.AND.osmcfi(1,1).EQ.LO.AND.
c     .       osmpei(1,1).EQ.LO)
c     .     CALL ER('966','Null data for ion power plot',*9997)
c         ylab     = 'powi'
c         elabs(1) = '    qi'
c         elabs(2) = '    pei'
c         elabs(3) = '    cfi'
c         elabs(4) = '    net'
c         ngs      = 4
c         btype    = 2
       ELSEIF (ytype.EQ.8) THEN
         IF (pinqe (1,1).EQ.LO.AND.osmcfe(1,1).EQ.LO.AND.
     .       osmpei(1,1).EQ.LO)
     .     CALL ER('966','Null data for electron power plot',*9997)

         ylab     = 'powe'
         elabs(1) = '    qe (Phelpi)'
         elabs(2) = '    qi (CX)'
         elabs(3) = '    qe (rec)'
         elabs(4) = '    qi (rec)'
         elabs(5) = '    net'

c         elabs(3) = '    cfe'
c         elabs(2) = '    mock'
c         elabs(2) = '    pei'
         ngs      = 4
         btype    = 2
       ELSEIF (ytype.EQ.9) THEN
         WRITE(6,*) 'FLAGS: ',pinion(1,1),pinrec(1,1),osmcfp(1,1)

         IF (pinion(1,1).EQ.LO.AND.pinrec(1,1).EQ.LO.AND.
     .       osmcfp(1,1).EQ.LO)
     .     CALL ER('966','Null data for mass plot',*9997)

         ylab     = 'par'
         elabs(1) = '    ion'
         elabs(2) = '    rec'
         elabs(3) = '    cfp (modified)'
         ngs      = 3
c         elabs(4) = '    net'
c         ngs      = 4
         btype    = 2


C...is cfp CORRECT NORMALLY?

         DO ir = 9, 29
           iks = 1
           ike = nks(ir)
c           iks = ikbound(ir,1)
c           ike = ikbound(ir,2)
           sump = 0.0
           sumn = 0.0
           sum  = 0.0
           DO ik = iks, ike
             IF     (ik.GT.1      .AND.ik.EQ.iks) THEN
               deltas = ksb(ik,ir) - kss(ik  ,ir)
             ELSEIF (ik.LT.nks(ir).AND.ik.EQ.ike) THEN
               deltas = kss(ik,ir) - ksb(ik-1,ir)
             ELSE
               deltas = ksb(ik,ir) - ksb(ik-1,ir)
             ENDIF
C             dist(ik,ir) = deltas * CalcDist(kss(ik,ir),ir,12)

c             dist(ik,ir) = deltas * CalcDist(kss(ik,ir),ir,11)
             IF     (dist(ik,ir).GE.0.0) THEN
               sump = sump + dist(ik,ir)
             ELSEIF (dist(ik,ir).LT.0.0) THEN
               sumn = sumn + dist(ik,ir)
             ENDIF
             sum = sum + dist(ik,ir)
           ENDDO

           rdum1 = -(knds(idds(ir,1)) * ABS(kvds(idds(ir,1))) +
     .               knds(idds(ir,2)) * ABS(kvds(idds(ir,2))))
c           rdum1 = -(knbs(iks,ir) * ABS(kvhs(iks,ir) / qt) +
c     .               knbs(ike,ir) * ABS(kvhs(ike,ir) / qt))

           CALL CalcIntegral4(pinion,iks,ike,ir,rdum2,2)
           CALL CalcIntegral4(pinrec,iks,ike,ir,rdum3,2)

           magcfp = -rdum1 - (rdum2 - rdum3)

           CALL CalcIntegral3(osmcfp,iks,ike,ir,rdum4,2)

           factp = 1.0
           factn = 1.0
           IF     (magcfp.GT.0.0.AND.sum.LT.0.0) THEN
             factp = -1.0 - 2.0 * sumn / sump
             WRITE(0,*) 'SCALING POSITIVE DISTRIBUTION ',factp
           ELSEIF (magcfp.LT.0.0.AND.sum.GT.0.0) THEN
             factn = -1.0 - 2.0 * sump / sumn
             WRITE(0,*) 'SCALING NEGATIVE DISTRIBUTION ',factn
           ENDIF

           DO ik = iks, ike
             IF     (ik.GT.1      .AND.ik.EQ.iks) THEN
               deltas = ksb(ik,ir) - kss(ik  ,ir)
             ELSEIF (ik.LT.nks(ir).AND.ik.EQ.ike) THEN
               deltas = kss(ik,ir) - ksb(ik-1,ir)
             ELSE
               deltas = ksb(ik,ir) - ksb(ik-1,ir)
             ENDIF
             IF (dist(ik,ir).GE.0.0) THEN
               dist(ik,ir) = dist(ik,ir) / ABS(sum) / deltas * factp
             ELSE
               dist(ik,ir) = dist(ik,ir) / ABS(sum) / deltas * factn
             ENDIF
           ENDDO

           CALL CalcIntegral3(dist  ,iks,ike,ir,rdum6,2)

           DO ik = iks, ike
             osmcfp(ik,ir) = ABS(magcfp) * dist(ik,ir)
           ENDDO

           CALL CalcIntegral3(osmcfp,iks,ike,ir,rdum5,2)

           WRITE(0,10) '966: IR,MAG,INT= ',ir,magcfp,rdum4,rdum5,
     .                                     rdum6,sum,sumn,sump
 10        FORMAT(A,I4,1P,7E10.2,0P)
         ENDDO
       ELSEIF (ytype.EQ.10) THEN
         IF (pinmp(1,1).EQ.LO.AND.osmmp(1,1).EQ.LO)
     .     CALL ER('966','Null data for momentum source plot',*9997)

         ylab     = 'mom'
         elabs(1) = '    cx (Pa/m)'
         elabs(2) = '    rec (Pa/m)'
         elabs(3) = '    osm (Pa/m)'
         elabs(4) = '    p (Pa)'
         elabs(5) = '    p (PIN) (Pa)'
c         elabs(5) = '    p (10%)'
         elabs(6) = '    v (10%)'
c         elabs(6) = '    v ( 1%)'
         ngs      = 6
         btype    = 2

         slopt2 = 2     

         plottype(1) = 2
         plottype(2) = 3
         plottype(3) = 4
         plottype(4) = 5
         plottype(5) = 6
         plottype(6) = 7

         CALL SetPlotComments(966,job,extra_comments,0,xsh)

c...     "Extended" colour set:
         DO in = ncols+1, (ncols+1) + 6
           colour(in) = in
         ENDDO
         CALL RGB
         CALL ColSet(1.0,0.0,0.0,ncols+1)
         CALL ColSet(0.0,1.0,0.0,ncols+2)
         CALL ColSet(0.0,0.0,1.0,ncols+3)
         CALL ColSet(0.5,0.7,0.2,ncols+4)
         CALL ColSet(0.5,0.0,0.5,ncols+5)
         CALL ColSet(0.0,0.5,0.5,ncols+6)

       ELSEIF (ytype.EQ.11) THEN
         IF (pinmp(1,1).EQ.LO.AND.osmmp(1,1).EQ.LO)
     .     CALL ER('966','Null data for viscosity plot',*9997)

         ylab     = 'mom'
         elabs(1) = '    cx '
         elabs(2) = '    vis'
         ngs      = 2
         btype    = 2
       ELSEIF (ytype.EQ.12) THEN
         IF (pinmp(1,1).EQ.LO.OR.osmmp(1,1).EQ.LO)
     .     CALL ER('966','Null data for mom. channel plot',*9997)

         ylab     = 'mom'
         elabs(1) = '    2+3'
         elabs(2) = '    4  '
         elabs(3) = '    9  '
         elabs(4) = '    sum'
         elabs(5) = '    pin'
         ngs      = 5
         btype    = 2
       ELSEIF (ytype.EQ.13.OR.ytype.EQ.14) THEN
         PLOTTYPE(1) = 1
         ylab     = ' ion_te / ion_tot'
         elabs(1) = '    ion_te / ion_tot'
         ngs      = 1
         btype    = 2

c...     Need to do some calculations here in order to get the count, or
c        weighted by ionisation, temperatures, er... whatever:

         ringnos(1) = 2
         pltlabs(1) = 'DIVERTOR'
         nplts      = 1

         IF (ytype.EQ.13) THEN
           pltlabs(1) = 'INNER DIVERTOR (DETACHED RINGS)'
         ELSE
           pltlabs(1) = 'OUTER DIVERTOR (DETACHED RINGS)'
         ENDIF

         i2 = 2

         DO i1 = 0, 31
           xdata(i1,i2) = REAL(i1)
         ENDDO

         rdum1 = 0.0
         DO ir = irsep, irwall-1
           IF ((osm_model(IKLO,ir).NE.24.AND.ytype.EQ.13).OR.
     .         (osm_model(IKHI,ir).NE.24.AND.ytype.EQ.14)) CYCLE

           IF (ytype.EQ.13) THEN
             iks = 1
             ike = ikmids(ir)
           ELSE
             iks = ikmids(ir) + 1
             ike = nks(ir)
           ENDIF
           CALL CalcIntegral3(pinion,iks,ike,ir,rdum2,2)
           rdum1 = rdum1 + rdum2
         ENDDO

         DO ir = irsep, irwall-1
           IF ((osm_model(IKLO,ir).NE.24.AND.ytype.EQ.13).OR.
     .         (osm_model(IKHI,ir).NE.24.AND.ytype.EQ.14)) CYCLE

           IF (ytype.EQ.13) THEN
             iks = 1
             ike = ikmids(ir)
           ELSE
             iks = ikmids(ir) + 1
             ike = nks(ir)
           ENDIF
           DO i1 = 1, 30
             DO ik = iks, ike
               IF (ktebs(ik,ir).GE.0.5*(xdata(i1-1,i2)+xdata(i1,i2))
     .        .AND.ktebs(ik,ir).LE.0.5*(xdata(i1+1,i2)+xdata(i1,i2))
     .           ) THEN
                 rdum2 = pinion(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir))
                 ydata(i1,i2) = ydata(i1,i2) + rdum2 / (rdum1 + 1.0E-8)
               ENDIF
             ENDDO
           ENDDO
         ENDDO


         rdum4 = 0.0
         DO ir = irsep, irwall-1
           IF ((osm_model(IKLO,ir).NE.24.AND.ytype.EQ.13).OR.
     .         (osm_model(IKHI,ir).NE.24.AND.ytype.EQ.14)) CYCLE
           IF (ytype.EQ.13) THEN
             iks = 1
             ike = ikmids(ir)
           ELSE
             iks = ikmids(ir) + 1
             ike = nks(ir)
           ENDIF
           DO ik = iks, ike
             rdum2 = pinion(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir))
             rdum4 = rdum4 + ktebs(ik,ir) * rdum2 / (rdum1 + 1.0E-8)
           ENDDO
         ENDDO

         rdum1 = 0.0
         rdum3 = 0.0
         DO i1 = 1, 30
           rdum1 = rdum1 + ydata(i1,i2)
           rdum3 = rdum3 + xdata(i1,i2) * ydata(i1,i2)
         ENDDO
         WRITE(0,*) '966 : 13 : SUM, AVG Te= ',rdum1,rdum3,rdum4


c         DO ir = irsep, irwall-1
c           IF (ytype.EQ.13) THEN
c             iks = 1
c             ike = ikmids(ir)
c           ELSE
c             iks = ikmids(ir) + 1
c             ike = nks(ir)
c           ENDIF
c           CALL CalcIntegral3(pinion,iks,ike,ir,rdum1,2)
c           DO i1 = 0, 31
c             xdata(i1,ir) = REAL(i1)
c           ENDDO
c           DO i1 = 1, 30
c             DO ik = iks, ike
c               IF (ktebs(ik,ir).GE.0.5*(xdata(i1-1,ir)+xdata(i1,ir))
c     .        .AND.ktebs(ik,ir).LE.0.5*(xdata(i1+1,ir)+xdata(i1,ir))
c     .           ) THEN
c                 rdum2 = pinion(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir))
c                 ydata(i1,ir) = ydata(i1,ir) + rdum2 / (rdum1 + 1.0E-8)
c               ENDIF
c             ENDDO
c           ENDDO
c           rdum1 = 0
c           DO i1 = 1, 30
c             rdum1 = rdum1 + ydata(i1,ir)
c           ENDDO
c           WRITE(0,*) '966 : 13 : SUM= ',rdum1
c         ENDDO

       ELSEIF (ytype.EQ.15) THEN
c...     Plot of temperatures with Thomson data:
         CALL RSet(mvals,MAXNPS*MAXPLTS*MAXNGS,LO)

c         ref2 = '996:15 Divertor Thomson / OSM temperature comparison'

         ylab = 'T (eV)'

         ngs      = 2
         elabs(1) = '    Te'
         elabs(2) = '    Ti'

         plottype(3) = -1
         plottype(4) = -2

         btype    = 1
         ind      = 2

       ELSEIF (ytype.EQ.16) THEN
c...     Plot of density with Thomson data:
         CALL RSet(mvals,MAXNPS*MAXPLTS*MAXNGS,LO)

c         ref2 = '996:16 Divertor Thomson / OSM density comparison'

         ylab     = 'n (m-3)'
         elabs(1) = '    n'
         ngs      = 1
         btype    = 1
         ind      = 1

         plottype(2) = -1
         plottype(3) = -2

       ELSEIF (ytype.EQ.17) THEN
c...     Ionisation rate:
         CALL RSet(mvals,MAXNPS*MAXPLTS*MAXNGS,LO)

c...TEMP--MOVE
         CALL LoadEIRENEAtomicData

         ylab     = 'ioniz rate (s-1)'
         elabs(1) = '    Rion'
         ngs      = 1
         btype    = 2
       ELSEIF (ytype.EQ.18) THEN
c...     Mach number:
         CALL RSet(mvals,MAXNPS*MAXPLTS*MAXNGS,LO)

         ylab     = 'Mach no.'
         elabs(1) = '    M'
         ngs      = 1
         btype    = 1
       ELSEIF (ytype.EQ.19) THEN
c...     Mach number:
         CALL RSet(mvals,MAXNPS*MAXPLTS*MAXNGS,LO)

         ylab     = 'gamma (m-2 s-1)'
         elabs(1) = '    G'
         ngs      = 1
         btype    = 1

       ELSEIF (ytype.EQ.20.or.ytype.eq.46) THEN
c...     Plot of Te with Thomson data and UEDGE:
         CALL RSet(mvals,MAXNPS*MAXPLTS*MAXNGS,LO)
         ylab      = 'Te (eV)'
         ngs       = 1
c         ngs       =  ngs + 1
         elabs(ngs)  = '    OSM'
         plotuedge = .false.
         slopt2    =  2
         btype     =  1
         ind       =  2
       ELSEIF (ytype.EQ.21) THEN
c...     Plot of Ti with UEDGE data:
         CALL RSet(mvals,MAXNPS*MAXPLTS*MAXNGS,LO)
         ylab      = 'Ti (eV)'
         ngs       = 1
         elabs(1)  = '    OSM'
         plotuedge = .true.
         btype     = 1
         ind       = 2
       ELSEIF (ytype.EQ.22.or.ytype.eq.47) THEN
c...     Plot of density with Thomson data and UEDGE data:
         CALL RSet(mvals,MAXNPS*MAXPLTS*MAXNGS,LO)
         ylab        = 'n (m-3)'
         ngs         =  1
         elabs(1)    = '    OSM'
         plotuedge   = .FALSE.
         slopt2      =  2
         btype       =  1
         ind         =  1
       ELSEIF (ytype.EQ.23) THEN
c...     Plot of Te and another case:
         CALL RSet(mvals,MAXNPS*MAXPLTS*MAXNGS,LO)
         ylab        = 'Te (eV)'
         ngs         = 1
         elabs(1)    = '     5.0 eV'
         plottype(3) = -1
         plottype(4) = -2
         slopt2    =  2
         btype    = 1
         ind      = 2
       ELSEIF (ytype.EQ.24) THEN
c...     Plot of Ti and another case:
         CALL RSet(mvals,MAXNPS*MAXPLTS*MAXNGS,LO)
         ylab        = 'Ti (eV)'
         ngs         = 1
         elabs(1)    = '    2/3'
         plottype(3) = -1
         plottype(4) = -2
         slopt2    =  2
         btype    = 1
         ind      = 2
       ELSEIF (ytype.EQ.25) THEN
c...     Plot of n and another case:
         CALL RSet(mvals,MAXNPS*MAXPLTS*MAXNGS,LO)
         ylab        = 'n (eV)'
         ngs         = 1
         elabs(1)    = '    5.0 eV'
         plottype(3) = -1
         plottype(4) = -2
         slopt2    =  2
         btype    = 1
         ind      = 1
       ELSEIF (ytype.EQ.26) THEN
c...     Plot of momentum loss and another case:
         CALL RSet(mvals,MAXNPS*MAXPLTS*MAXNGS,LO)
         ylab     = 'Dalpha (...)'
         ngs      = 2
         elabs(1) = '    SOL22'
         btype    = 2
       ELSEIF (ytype.EQ.27) THEN
c...     Plot of D density and another case:
         CALL RSet(mvals,MAXNPS*MAXPLTS*MAXNGS,LO)
         ylab     = 'PINATOM (...)'
         ngs      = 2
         elabs(1) = '    SOL22'
         btype    = 2
       ELSEIF (ytype.EQ.28) THEN
c...     Plot of D2 density and another case:
         CALL RSet(mvals,MAXNPS*MAXPLTS*MAXNGS,LO)
         ylab     = 'PINMOL (...)'
         ngs      = 2
         elabs(1) = '    SOL22'
         btype    = 2
       ELSEIF (ytype.EQ.29) THEN
         ylab     = 'bratio'
         elabs(1) = '    brat'
         ngs      = 1
         btype    = 2
       ELSEIF (ytype.EQ.30) THEN
c...     Plotting the pressure along a toroidal chord:
         ylab     = 'p (mTorr)'
         elabs(1) = '    p'
         ngs      = 1
         btype    = 2

c...     Find the center point for each cell listed in the OUT data
c        line and find all the pressure data along that toroidal 
c        chord:
         DO ip = 1, nplts
           cell = ringnos(ip)       

           plottype(ip) = 1

           WRITE(6,*) 'TARGET CELL:',cell

           IF (asc_3dmode.EQ.1) THEN
             WRITE(0,*) '966: ASC3DMODE=1 NO LONGER SUPPORTED'
             STOP
           ELSEIF (asc_3dmode.EQ.2) THEN

c...         Loop over all additional cells and assemble data
c            for all cells that contain the center point:
             iks = 1
             ike = 0
             shift = 1 + eirnpgdat
             fact  = ECH * 1.0E+06 * 0.67 * 7.502
             DO cut1 = 1, ascncut
               ike   = ike + 1
               cell2 = cell+(cut1-1)*asc_ncell+shift
               xdata(ike,ip) = 0.5*(asc_zmin3D(cut1)+asc_zmax3D(cut1))
c...           D2:
c               ydata(ike,ip) = pinasd(cell2,2,3,1) * fact
               ydata(ike,ip) = pinasd(cell2,2,3,2) * fact
c...           D:
c               ydata(ike,ip) = pinasd(cell2,2,1,1) * fact
c               ydata(ike,ip) = pinasd(cell2,2,1,2) * fact
c               ydata(ike,ip) = ydata(ike,ip)+pinasd(cell2,2,1,2) * fact
             ENDDO

           ENDIF

           pnks(ip) = ike
 
         ENDDO
       ELSEIF (ytype.GE.31.AND.ytype.LE.35) THEN
c...     Plotting the pressure along a toroidal chord of the standard grid:
         ylab     = 'p (mTorr)'
         elabs(1) = '    p'
         ngs      = 1
         btype    = 2

c...     Find the center point for each cell listed in the OUT data
c        line and find all the pressure data along that toroidal 
c        chord:
         DO ip = 1, nplts
           cell = ringnos(ip)       

           plottype(ip) = 1

           WRITE(6,*) 'TARGET CELL:',cell

           ir = irsep+1
           ik = cell

           DO i1 = 1, eirnsdtor
             i2 = (i1-1)*MAXBGK+(ytype-31)*5
c             WRITE(0,'(A,3I3,F10.4,1P,5E12.4,0P)') 
c     .          '-->',ik,ir,(ytype-31)+1,
c     .          eirsdtor(i1),(pinbgk(ik,ir,i2+i3),i3=1,5)
             WRITE(6,'(A,3I3,F10.4,1P,5E12.4,0P)') 
     .          '-->',ik,ir,(ytype-31)+1,
     .          eirsdtor(i1),(pinbgk(ik,ir,i2+i3),i3=1,5)
           ENDDO
           WRITE(0,*) 
           WRITE(6,*) 
         ENDDO
 
c         CALL FRAME

         RETURN

       ELSEIF (ytype.GE.36.AND.ytype.LE.45) THEN
         i1 = ytype-36+1
         ref2 = ref(1:LEN_TRIM(ref))//' '//
     .          eirtally(i1,2)(1:LEN_TRIM(eirtally(i1,2)))
         ylab     = eirtally(i1,2)(1:LEN_TRIM(eirtally(i1,2)))//' ('//
     .              eirtally(i1,4)(1:LEN_TRIM(eirtally(i1,4)))//')'
         elabs(1) = '    '//eirtally(i1,3)(1:LEN_TRIM(eirtally(i1,3)))

         elabs(1) = '    MAR-EIRENE'
         elabs(2) = '    MAR-CALCULATED FROM nD2'


         ngs      = 1
         btype    = 2

         IF (ytype.EQ.39) ngs = 2

       ELSEIF (ytype.EQ.48) THEN
         ylab     = 'ionisation'
         elabs(1) = '    ATM'
         elabs(2) = '    MOL'
         elabs(3) = '    TES'
         elabs(4) = '    NET'
         elabs(5) = '    ION'
         ngs      = 4
         btype    = 2

         plottype(1) = 2
         plottype(2) = 3
         plottype(3) = 4
         plottype(4) = 5
         plottype(5) = 6

         slopt2 = 2

       ELSEIF (ytype.EQ.49) THEN
         ylab     = 'line (...)'
         elabs(1) = '    EIRENE'
         elabs(2) = '    camera'
         ngs      = 2
         btype    = 2

         plottype(1) = 2
         plottype(2) = 3
         plottype(3) = 4
         plottype(4) = 5
         plottype(5) = 6

         slopt2 = 2

c...     Load and map toroidal inversion data:
         scale = 0.0
         WRITE(fname,'(1024X)')
         READ(5,'(A256)') dummy
         IF   (dummy(8:11).EQ.'File'.OR.dummy(8:11).EQ.'file'.OR.
     .         dummy(8:11).EQ.'FILE') THEN
           READ(dummy,*) cdum1,line,optscale,scale,fname
           osmtmp = 0.0
           CALL LoadCameraData(osmtmp,fname,scale)
         ELSE
           CALL ER('966:49','Data file not specified',*99)
         ENDIF

       ELSEIF (ytype.EQ.50) THEN
c...     D MFP:
         CALL RSet(mvals,MAXNPS*MAXPLTS*MAXNGS,LO)
         ylab     = 'PINATOM (...)'
         ngs      = 1
         elabs(1) = '    ...'
         btype    = 2

         slopt2    =  2
         CALL LoadEIRENEAtomicData

       ELSEIF (ytype.EQ.51) THEN
c...     Plotting the pressure along a radial chord:
         ylab     = '...'
         elabs(1) = '    p'
         ngs      = 1
         btype    = 2
         set      = 2

c...     Find the center point for each cell listed in the OUT data
c        line and find all the pressure data along that toroidal 
c        chord:
         DO ip = 1, nplts
           cell = -1

           plottype(ip) = 1

           WRITE(6,*) 'TARGET CELL:',cell

           yval = -0.58
           zval =  0.13

           IF     (asc_3dmode.EQ.1) THEN
             WRITE(0,*) '966: ASC3DMODE=1 NO LONGER SUPPORTED'
             STOP
           ELSEIF (asc_3dmode.EQ.2) THEN

c...         Loop over all additional cells and assemble data
c            for all cells spanning YVAL:
             iks = 1
             ike = 0
             shift = 1 + eirnpgdat
             fact  = ECH * 1.0E+06 * 0.67 * 7.502

             DO cut1 = 1, ascncut
               IF (zval.GE.asc_zmin3D(cut1).AND.
     .             zval.LE.asc_zmax3D(cut1)) EXIT
             ENDDO
             IF (cut1.EQ.ascncut+1) 
     .         CALL ER('966:51','Toroidal section not found',*99)

             ike = 0
            
             DO cell = 1, asc_ncell
               cell2 = cell+(cut1-1)*asc_ncell+shift
c...           Find x,ymin and x,ymax for the cell:
               xmin =  HI
               ymin =  HI
               xmax = -HI
               ymax = -HI
               DO i2 = 1, ascnvertex(cell)
                 xmin = MIN(xmin,ascvertex(i2*2-1,cell))
                 ymin = MIN(ymin,ascvertex(i2*2  ,cell))
                 xmax = MAX(xmax,ascvertex(i2*2-1,cell))
                 ymax = MAX(ymax,ascvertex(i2*2  ,cell))
               ENDDO

               IF ( 0.5 * (xmin + xmax).LT.0.630) CYCLE
c               IF ( 0.5 * (xmin + xmax).LT.0.659) CYCLE

               IF (yval.GE.ymin.AND.yval.LT.ymax) THEN
                 ike = ike + 1
c                 WRITE(0,*) '966:51:',cell
                 xdata(ike,ip) = 0.5 * (xmin + xmax)
                 IF (ip.EQ.1) THEN
c...               D2 pressure:
c                   ydata(ike,ip) = pinasd(cell2,2,3,set) * fact

                   ydata(ike,ip) = pinasd(cell2,2,3,set) /
     .                             pinasd(cell2,1,3,set)

 
c..BGK y-vel pressure:
c                   ydata(ike,ip) = SQRT(pinasd(cell2,5,8,set)**2+
c     .                                  pinasd(cell2,5,8,set)**2 + 
c     .                                  pinasd(cell2,5,8,set)**2) * 
c     .                             1.0E+02
                 ELSEIF (ip.EQ.2) THEN
c...               D pressure:
c                   ydata(ike,ip) = pinasd(cell2,2,1,set) * fact

c..BGK x-vel pressure:
                   ydata(ike,ip) = pinasd(cell2,3,8,set) * 1.0E+02
c                   ydata(ike,ip) = pinasd(cell2,2,7,set) * fact
                 ELSEIF (ip.EQ.3) THEN
c...               D2 density:
                   ydata(ike,ip) = pinasd(cell2,1,3,set) * 1.0E+06
                 ELSEIF (ip.EQ.4) THEN
c...               D2 pressure:
                   ydata(ike,ip) = pinasd(cell2,2,3,set) * fact
c...               D density:
c                   ydata(ike,ip) = pinasd(cell2,1,1,set) * 1.0E+06
                 ENDIF
               ENDIF
             ENDDO

c...         Sort points from inside to outside:
             cont = .TRUE.
             DO WHILE (cont)
               cont = .FALSE.
               DO i1 = 1, ike-1
                 IF (xdata(i1,ip).GT.xdata(i1+1,ip)) THEN
                   xdata(0,ip) = xdata(i1,ip)
                   ydata(0,ip) = ydata(i1,ip)
                   xdata(i1,ip) = xdata(i1+1,ip)
                   ydata(i1,ip) = ydata(i1+1,ip)
                   xdata(i1+1,ip) = xdata(0,ip)
                   ydata(i1+1,ip) = ydata(0,ip)
                   cont = .TRUE.
                 ENDIF
               ENDDO
             ENDDO

           ENDIF

c           ydata(ike,4) = 10.0

           pnks(ip) = ike
 
         ENDDO


       ELSEIF (ytype.EQ.52) THEN
c...     Plot of pressure from Thomson data Te and ne:
         CALL RSet(mvals,MAXNPS*MAXPLTS*MAXNGS,LO)
         ylab      = 'p_e_st (mTorr)'
         ngs       = 1
         elabs(ngs)  = '    OSM'
         plotuedge = .FALSE.
         slopt2    =  2
         btype     =  1
         ind       =  3

       ELSEIF (ytype.EQ.53) THEN
c...     E-field:
         CALL RSet(mvals,MAXNPS*MAXPLTS*MAXNGS,LO)
         ylab      = 'E-field [V m-1]'
         ngs       = 1
         elabs(ngs)  = '    kes'
         plotuedge = .FALSE.
         slopt2    =  2
         btype     =  2

       ELSEIF (ytype.EQ.54) THEN
c...     Electric potential:
         CALL RSet(mvals,MAXNPS*MAXPLTS*MAXNGS,LO)
         ylab      = 'E-field [V m-1]'
         ngs       = 1
         elabs(ngs)  = '    kes'
         plotuedge = .FALSE.
         slopt2    =  2
         btype     =  2

       ELSE
         CALL ER('966','Invalid plot (YTYPE) specification',*99)
       ENDIF
c
c *** ENDIF ***
c


       CALL CustomisePlot(title,xlab,ylab,elabs)
c
c
c
       ncase = 0
       IF (ytype.EQ.7.OR.
     .     ytype.EQ.23.OR.ytype.EQ.24.OR.ytype.EQ.25.OR.
     .     ytype.EQ.26.OR.ytype.EQ.27.OR.ytype.EQ.28.OR.
     .     ytype.EQ.50) THEN

 

 12      READ(5,'(A80)',END=98) graph1
         WRITE(0,*) 'GRAPH1;',graph1
         IF (graph1(8:11).EQ.'Case'.OR.graph1(8:11).EQ.'CASE'.OR.
     .       graph1(8:11).EQ.'case') THEN
           READ(graph1,*) cdum1,elabs(ngs+1),cname,step,cmnd
           CALL LoadData(cname,cmnd,step,-1)
c...       Copy grid data (grids better be the same):
           ngs = ngs + 1
           plottype(ngs) = ngs
           ncase = ncase + 1 
           DO ir = 1, nrs
             casendata(ir,ncase) = nks(ir)
             DO ik = 1, nks(ir)
               casedata(ik,ir,0,ncase) = kss(ik,ir)
               casedata(ik,ir,1,ncase) = ktebs(ik,ir)
               casedata(ik,ir,2,ncase) = ktibs(ik,ir)
               casedata(ik,ir,3,ncase) = knbs (ik,ir)
               casedata(ik,ir,4,ncase) = pinalpha(ik,ir)
               casedata(ik,ir,5,ncase) = pinatom(ik,ir)
               casedata(ik,ir,6,ncase) = pinmol(ik,ir)
               casedata(ik,ir,7,ncase) = 
     .           CalcPressure(knbs(ik,ir),ktebs(ik,ir),
     .                        ktibs(ik,ir),kvhs(ik,ir))
             ENDDO
           ENDDO
c...       Copy target data:
           DO in = 1, nds
             casetdata(in,1,ncase) = kteds(in)
             casetdata(in,2,ncase) = ktids(in)
             casetdata(in,3,ncase) = knds (in)
           ENDDO
c...       Look for more data:        
           GOTO 12

         ELSE
c...       Restore original solution:
           IF (ncase.NE.0) THEN
             CALL LoadData(' ',' ',NULL,-2)                   
             WRITE(6,*) 'CASE DATA:'
             DO i1 = 1, ncase
               DO ir = irsep, nrs
                 IF (idring(ir).EQ.BOUNDARY) CYCLE
                 DO ik = 1, nks(ir)
                    WRITE(6,'(3I6,1P,10(E10.2))')
     .                i1,ik,ir,(casedata(ik,ir,i2,i1),i2=1,6)
                 ENDDO
               ENDDO
             ENDDO
           ENDIF

           BACKSPACE 5
         ENDIF

       ENDIF
c
c
c
c...   Check for UEDGE data:
       plotuedge = .false.
       IF (plotuedge) THEN
         ngs           = ngs + 1
         elabs(ngs)    = '    UEDGE'
         plottype(ngs) = 1
       ENDIF

c
c
c...   Check for AVERAGE data:

c       plotav = .TRUE.

       if (ytype.eq.46.or.ytype.eq.47) plotav = .true.

       IF (plotav) THEN
         ngs           = ngs + 1
         elabs(ngs)    = '    AVER.'
         plottype(ngs) = 1
       ENDIF
c
c
c
       IF (ytype.EQ.15.OR.ytype.EQ.16.OR.ytype.EQ.20.OR.
     .     ytype.EQ.22.OR.ytype.EQ.23.OR.ytype.EQ.24.OR.
     .     ytype.EQ.25.OR.ytype.EQ.52) THEN

c
c jdemod
c
c         Calculate averages for future use - even if not plotted
c
          ndatasets = 0
          call rzero(avdata,maxnks*maxnrs)
          call rzero(ncounts,maxnks*maxnrs)
c
c jdemod
c

c...     Look for Thomson data to plot:
15       READ(5,'(A80)') graph1
         IF (graph1(8:11).EQ.'Thom'.OR.graph1(8:11).EQ.'THOM'.OR.
     .       graph1(8:11).EQ.'thom') THEN
          ngdata        = ngdata + 1
          ngs           = ngs    + 1
          plottype(ngs) = -(ngdata+2)
c          plottype(ngs) = -55
          READ(graph1,*) cdum1,elabs(ngs)
          BACKSPACE 5
          xsh = 0.0
          ysh = 0.0
          CALL LoadThomsonData(thnnraw,thnraw,MAXTDAT,MAXCOLS,xsh,ysh,1)
c
          ndatasets = ndatasets + 1
c
c jdemod
c
c         Accumulate average data
c
          do i1 = 1,thnnraw
               ik = int(thnraw(i1,3))
               ir = int(thnraw(i1,4))
               avdata(ik,ir) = avdata(ik,ir) + thnraw(i1,ind+4)             
               ncounts(ik,ir) = ncounts(ik,ir) + 1.0
          end do
c
c jdemod
c

          CALL AsgnThomsonData(thnnraw,thnraw,1,
     .                        gxdata(1,1,1,ngdata),gydata(1,1,1,ngdata),
     .                        gndata(  1,1,ngdata),
     .                        MAXNPS,MAXTYP,MAXTDAT,MAXCOLS)

          CALL AsgnThomsonData(thnnraw,thnraw,2,
     .                        gxdata(1,1,1,ngdata),gydata(1,1,1,ngdata),
     .                        gndata(  1,1,ngdata),
     .                        MAXNPS,MAXTYP,MAXTDAT,MAXCOLS)

c
c          CALL AsgnThomsonData(thnnraw,thnraw,ind,
c     .                        gxdata(1,1,1,ngdata),gydata(1,1,1,ngdata),
c     .                        gndata(  1,1,ngdata),
c     .                        MAXNPS,MAXTYP,MAXTDAT,MAXCOLS)

c...      Calculate pressure:
          DO ir = 2, nrs
            IF (idring(ir).EQ.BOUNDARY) CYCLE
            gndata(ir,3,ngdata) = gndata(ir,1,ngdata)
            DO i1 = 1, gndata(ir,3,ngdata)
              gxdata(i1,ir,3,ngdata) = gxdata(i1,ir,1,ngdata)
              gydata(i1,ir,3,ngdata) = 
     .          CalcPressure(gydata(i1,ir,1,ngdata),
     .                       gydata(i1,ir,2,ngdata),
     .                       gydata(i1,ir,2,ngdata)*1.0,0.0)*ECH
            ENDDO
          ENDDO

          GOTO 15     
         ELSE
            BACKSPACE 5
c
c jdemod
c
c           Process Average data
c
            if (ndatasets.gt.0) then 
c
c              Set AVDATA with no average to LO for MVALS compatibility
c
               do ir = 1,nrs
                  do ik = 1,nks(ir)
                     if (ncounts(ik,ir).gt.0.0) then
                        avdata(ik,ir) = avdata(ik,ir) / ncounts(ik,ir)
                     else
                        avdata(ik,ir) = LO  
                     endif  
c
c                    Write the average data out to unit 6 
c
                     write(6,'(a,4i6,(1x,g12.5))') 'TS Averages:',
     >                        ytype,ind,ik,ir,avdata(ik,ir)
                  end do
               end do
               



            end if
c
c jdemod
c
         ENDIF

c...     Look for RCP data plot - from experimantal data file:
         READ(5,'(A80)',END=19) graph1
         IF (graph1(8:11).EQ.'RCP '.OR.graph1(8:11).EQ.'Rcp '.OR.
     .       graph1(8:11).EQ.'rcp ') THEN
           CALL LoadRCPData(graph1,thnnraw,thnraw,
     .                      MAXTDAT,MAXCOLS,xsh,ysh)
           ngdata        = ngdata + 1
           ngs           = ngs + 1
           plottype(ngs) = -(ngdata+2)
           READ(graph1,*) cdum1,elabs(ngs)
           CALL AsgnThomsonData(thnnraw,thnraw,ind,
     .                        gxdata(1,1,1,ngdata),gydata(1,1,1,ngdata),
     .                        gndata(  1,1,ngdata),
     .                        MAXNPS,MAXTYP,MAXTDAT,MAXCOLS)
         ELSE
           BACKSPACE 5
         ENDIF
19       CONTINUE


c...     Look for target data to plot - from experimantal data file:
         READ(5,'(A80)') graph1
         IF (graph1(8:11).EQ.'Targ'.OR.graph1(8:11).EQ.'TARG'.OR.
     .       graph1(8:11).EQ.'targ') THEN
           CALL LoadTargetData(graph1,targdat,MAXTDAT,MAXCOLS)
           ngs           = ngs + 1
           plottype(ngs) = -(ngs+1)
           elabs   (ngs) = '    LP'
c           plottype(ngs) = 2
           langdat = .TRUE.
         ELSE
           BACKSPACE 5
         ENDIF
20       CONTINUE
         CALL SetPlotComments(966,job,extra_comments,0,xsh)
       ENDIF
c
c jdemod
c
       IF (ytype.EQ.46.OR.ytype.EQ.47) THEN

          call rzero(avdata,maxnks*maxnrs)
          call rzero(ncounts,maxnks*maxnrs)

          ndatasets = 0

c...     Look for Thomson data to average:
17       READ(5,'(A80)',END=22) graph1
         IF (graph1(8:11).EQ.'Thom'.OR.graph1(8:11).EQ.'THOM'.OR.
     .       graph1(8:11).EQ.'thom') THEN
c
c            READ(graph1,*) cdum1,elabs(ngs)
c
            BACKSPACE 5
            xsh = 0.0
            ysh = 0.0
            CALL LoadThomsonData(thnnraw,thnraw,MAXTDAT,MAXCOLS,
     .                           xsh,ysh,1)

            ndatasets = ndatasets+1

            do i1 = 1,thnnraw

               ik = int(thnraw(i1,3))
               ir = int(thnraw(i1,4))
     
               avdata(ik,ir) = avdata(ik,ir) + thnraw(i1,ind+4)             
               ncounts(ik,ir) = ncounts(ik,ir) + 1.0

            end do

            GOTO 17
         ELSE
c
            BACKSPACE 5
c
c           Process averages 
c
            if (ndatasets.gt.0) then 
c
c              Set AVDATA with no average to LO for MVALS compatibility
c
               do ir = 1,nrs
                  do ik = 1,nks(ir)
                     if (ncounts(ik,ir).gt.0.0) then
                        avdata(ik,ir) = avdata(ik,ir) / ncounts(ik,ir)
                     else
                        avdata(ik,ir) = LO 
                     endif  
                  end do
               end do
            end if
         ENDIF
c
c jdemod   
c

c...     Look for target data to plot - from experimantal data file:
         READ(5,'(A80)',END=22) graph1
         IF (graph1(8:11).EQ.'Targ'.OR.graph1(8:11).EQ.'TARG'.OR.
     .       graph1(8:11).EQ.'targ') THEN
           CALL LoadTargetData(graph1,targdat,MAXTDAT,MAXCOLS)
           ngs           = ngs + 1
           plottype(ngs) = -(ngs+1)
           elabs   (ngs) = '    LP'
c           plottype(ngs) = 2
           langdat = .TRUE.
         ELSE
           BACKSPACE 5
         ENDIF
22       CONTINUE
         CALL SetPlotComments(966,job,extra_comments,0,xsh)
       ENDIF


c
c
c
       sctype = iopt
c
c       if (sctype.lt.1.or.sctype.gt.4) sctype = 1
c
c      IF SCTYPE = 7 - check for scaling information
c
       if (sctype.eq.7) then 

c...     Look for scaling data for plot 
         READ(5,'(A80)',END=30) graph1
         IF (graph1(8:12).EQ.'Scale'.OR.graph1(8:12).EQ.'SCALE'.OR.
     .       graph1(8:12).EQ.'scale') THEN

             read(graph1,*) cdum1,pltmin,pltmax
c
             call rinit(pltmins,maxplts,pltmin)
             call rinit(pltmaxs,maxplts,pltmax)
c
             write(6,'(a,2(1x,g12.5))') 'SCALE:',pltmin,pltmax 
c
         ELSE
           CALL ER('966','SCLTYPE=7 but scale data not found in'//
     .             'input file',*99)
c           BACKSPACE 5     
         ENDIF
c
       endif

30     continue

c
       NPLOTS = NPLOTS + 1
       WRITE (IPLOT,9012) NPLOTS,REF

c       CALL rzero (mvals,maxnks*maxngs*maxplts)
c
c
       DO ip = 1, nplts

c         WRITE(0,*) '966: IP= ',ip

         IF (ytype.EQ.30.OR.ytype.EQ.51) THEN
c...       No importance:
           ir = irsep
         ELSE
           ir = ringnos(ip)
         ENDIF

         IF (idring(ir).EQ.-1) CYCLE

         IF (ytype.EQ.49) THEN
           dmax = LO
           gmax = LO
           IF     (optscale.EQ.0) THEN
           ELSEIF (optscale.EQ.1) THEN
c...         Scale accoring to local peak on outside of ring:
             DO ik = nks(ir)/2, nks(ir)
               dmax = MAX(dmax,pinline(ik,ir,6,line))
               gmax = MAX(gmax,osmtmp(ik,ir))
             ENDDO
             IF (line.EQ.2) 
     .         dmax = dmax * (6.63E-34*3.0E+08) / (4340.0*1.0E-10)
           ELSEIF (optscale.EQ.2) THEN
c...         Scale accoring to local peak on outside separatrix:
             DO ik = nks(irsep)/2, nks(irsep)
               dmax = MAX(dmax,pinline(ik,irsep,6,line))
               gmax = MAX(gmax,osmtmp(ik,irsep))
             ENDDO
             IF (line.EQ.2) 
     .         dmax = dmax * (6.63E-34*3.0E+08) / (4340.0*1.0E-10)
             WRITE(0,*) '966 SCALING:',dmax/gmax
           ELSEIF (optscale.EQ.3) THEN
c *HARDCODED*
             dmax = 2.1E+13
             gmax = 1.0
           ELSE
             CALL ER('966:49','UNKNOWN CAMERA SCALING OPTION',*99)
           ENDIF
         ENDIF

         IF (ytype.NE.13.AND.ytype.NE.14) THEN
           IF (osm_model(IKLO,ir).EQ.24) THEN
             grm_shade(IKLO,ip) = kss(ikfluid(IKLO,ir),ir)
c             grm_shade(IKLO,ip) = kss(ikbound(ir,IKLO),ir)
             DO ik = 0, nks(ir)
               grm_cell(ik,ip) = ksb(ik,ir)
             ENDDO
           ENDIF
           IF (osm_model(IKHI,ir).EQ.24) THEN
             grm_shade(IKHI,ip) = kss(ikfluid(IKHI,ir),ir)
c             grm_shade(IKHI,ip) = kss(ikbound(ir,IKHI),ir)
             DO ik = 0, nks(ir)
               grm_cell(ik,ip) = ksb(ik,ir)
             ENDDO
           ENDIF
         ENDIF
c
c...     Set bounds of data index:
         IF (ytype.EQ.7.OR.
     .       ytype.EQ.15.OR.ytype.EQ.16.OR.ytype.EQ.20.OR.
     .       ytype.EQ.21.OR.ytype.EQ.22.OR.ytype.EQ.23.OR.
     .       ytype.EQ.24.OR.ytype.EQ.25.or.ytype.eq.46.or.
     .       ytype.eq.47.OR.ytype.EQ.50.OR.ytype.EQ.52) THEN
           iks = 1
           i1 = 0
           in = 0
           IF (btype.EQ.1) in = 1
           DO i2 = 1, ngdata
             CALL LoadArray(mouts(1+in,ip),i1,gxdata(1,ir,ind,i2),1,
     .                      gndata(ir,ind,i2))
           ENDDO
           DO i2 = 1, ncase
             CALL LoadArray(mouts(1+in,ip),i1,casedata(1,ir,0,i2),1,
     .                      casendata(ir,i2))
           ENDDO
           CALL LoadArray(mouts(1+in,ip),i1,kss(1,ir),1,nks(ir))
           ike = i1
c           IF (ir.EQ.irsep) THEN
c             DO I1 = 1, ike
c               WRITE(0,*) 'KIE:',i1,kss(MIN(i1,nks(ir)),ir),
c     .                     casedata(MIN(i1,casendata(ir,1)),ir),
c     .                     mouts(i1+in,ip)
c             ENDDO
c           ENDIF
c           WRITE(0,*) 'IKE:',ir,ike,nks(ir),ngdata

c           ike = nks(ir)
c           DO i2 = 1, ngdata
c             ike = ike + gndata(ir,ind,i2)
c           ENDDO
         ELSEIF (ytype.EQ.13.OR.ytype.EQ.14) THEN
           iks = 1
           ike = 30
         ELSEIF (ytype.EQ.30.OR.ytype.EQ.51) THEN
           iks = 1
           ike = pnks(ip)
         ELSE
           iks = 1
           ike = nks(ir)
         ENDIF

         if     (btype.EQ.2) THEN
            in  = 0
            inc = 0
         elseif (btype.EQ.1) THEN
           ir = ringnos(ip)
           in = 1
           inc = 2
           mouts(1,ip) = 0.0
           if     (xtype.EQ.1.OR.xtype.EQ.6.OR.xtype.EQ.8) then
             mouts(ike+inc,ip) = ksmaxs(ir)
           ELSEIF (xtype.EQ.9) THEN
             mouts(1      ,ip) = 0.0
             mouts(ike+inc,ip) = ksmaxs(ir)
           elseIF (xtype.EQ.2) THEN
             mouts(1,ip) = thetat(idds(ir,2))
             mWIDS(1,ip) = thetag(1,ir) - thetat(idds(ir,2))
             mouts(nks(ir)+inc,ip) = thetat(idds(ir,1))
             mwids(nks(ir)+inc,ip) = thetat(idds(ir,1)) -
     .                               thetag(nks(ir),ir)
           ELSEIF (xtype.EQ.3) THEN
             mWIDS(1,ip) = kps(1,ir)
             mouts(nks(ir)+inc,ip) = kpmaxs(ir)
             mwids(nks(ir)+inc,ip) = kpmaxs(ir) - kps(nks(ir),ir)
           ELSEIF (xtype.EQ.4) THEN
             mWIDS(1,ip) = kss(1,ir) / ksmaxs(ir)
             mouts(nks(ir)+inc,ip) = 1.0
             mwids(nks(ir)+inc,ip) = 1.0 - kss(nks(ir),ir) / ksmaxs(ir)
           ELSE
             CALL ER('964','Invalid x-axis data',*9997)
           endif

           IF (xtype.EQ.9) THEN
             id1 = idds(ir,1)
             id2 = idds(ir,2)
           ELSE
             id1 = idds(ir,2)
             id2 = idds(ir,1)
           ENDIF

           IF     (ytype.EQ.1) THEN
             mvals(1          ,ip,1) = knds(idds(ir,2))
             mvals(nks(ir)+inc,ip,1) = knds(idds(ir,1))
             mvals(1          ,ip,2) = kteds(idds(ir,2))
             mvals(nks(ir)+inc,ip,2) = kteds(idds(ir,1))
             mvals(1          ,ip,3) = ktids(idds(ir,2))
             mvals(nks(ir)+inc,ip,3) = ktids(idds(ir,1))
             mvals(1          ,ip,4) = kvds (id1)
             mvals(nks(ir)+inc,ip,4) = kvds (id2)
             mvals(1          ,ip,5) = kvds (id1) * knds(id1)
             mvals(nks(ir)+inc,ip,5) = kvds (id2) * knds(id2)
           ELSEIF (ytype.EQ.2) THEN
             mvals(1          ,ip,1) = kteds(idds(ir,2))
             mvals(nks(ir)+inc,ip,1) = kteds(idds(ir,1))
             mvals(1          ,ip,2) = ktids(idds(ir,2))
             mvals(nks(ir)+inc,ip,2) = ktids(idds(ir,1))
           ELSEIF (ytype.EQ.3) THEN
             mvals(1          ,ip,1) = kteds(id1)
             mvals(nks(ir)+inc,ip,1) = kteds(id2)
             mvals(1          ,ip,2) = ktids(id1)
             mvals(nks(ir)+inc,ip,2) = ktids(id2)
             mvals(1          ,ip,3) = kvds (id1) * knds(id1)
             mvals(nks(ir)+inc,ip,3) = kvds (id2) * knds(id2)
           ELSEIF (ytype.EQ.5) THEN
             fact = crmb * AMU / ECH

             mvals(1          ,ip,1) = kteds(id1)
             mvals(nks(ir)+inc,ip,1) = kteds(id2)
             mvals(1          ,ip,2) = knds (id1)
             mvals(nks(ir)+inc,ip,2) = knds (id2)
             mvals(1          ,ip,3) = knds(id1) *
     .         (kteds(id1) + ktids(id1) + fact * kvds(id1)**2.0)
             mvals(nks(ir)+inc,ip,3) = knds(id2) *
     .         (kteds(id2) + ktids(id2) + fact * kvds(id2)**2.0)
           ELSEIF (ytype.EQ.15) THEN
             mvals(1      ,ip,1) = kteds(id1)
             mvals(ike+inc,ip,1) = kteds(id2)
             mvals(1      ,ip,2) = ktids(id1)
             mvals(ike+inc,ip,2) = ktids(id2)
             IF (langdat) THEN
               mvals(1      ,ip,ngs) = targdat(id1,1)
               mvals(ike+inc,ip,ngs) = targdat(id2,1)
             ENDIF
           ELSEIF (ytype.EQ.16) THEN
             mvals(1      ,ip,1) = knds(idds(ir,2))
             mvals(ike+inc,ip,1) = knds(idds(ir,1))
             IF (langdat) THEN
               mvals(1      ,ip,ngs) = targdat(id1,3)
               mvals(ike+inc,ip,ngs) = targdat(id2,3)
             ENDIF
           ELSEIF (ytype.EQ.18) THEN
             mvals(1      ,ip,1) = kvds(id1) /
     .                        GetCs(kteds(idds(ir,2)),ktids(idds(ir,2)))
             mvals(ike+inc,ip,1) = kvds(id2) /
     .                        GetCs(kteds(idds(ir,1)),ktids(idds(ir,1)))
           ELSEIF (ytype.EQ.19) THEN
             mvals(1      ,ip,1) = kvds(id1) * knds(id1)
             mvals(ike+inc,ip,1) = kvds(id2) * knds(id2)
           ELSEIF (ytype.EQ.20.or.ytype.eq.46) THEN
             mvals(1      ,ip,1) = kteds(id1)
             mvals(ike+inc,ip,1) = kteds(id2)
             IF (langdat) THEN
               mvals(1      ,ip,ngs) = targdat(id1,1)
               mvals(ike+inc,ip,ngs) = targdat(id2,1)
             ENDIF
           ELSEIF (ytype.EQ.21) THEN
             mvals(1      ,ip,1) = ktids(id1)
             mvals(ike+inc,ip,1) = ktids(id2)
           ELSEIF (ytype.EQ.22.or.ytype.eq.47) THEN
             mvals(1      ,ip,1) = knds(id1)
             mvals(ike+inc,ip,1) = knds(id2)
             IF (langdat) THEN
               mvals(1      ,ip,ngs) = targdat(id1,3)
               mvals(ike+inc,ip,ngs) = targdat(id2,3)
             ENDIF
           ELSEIF (ytype.EQ.23) THEN
             mvals(1      ,ip,1) = kteds(id1)
             mvals(ike+inc,ip,1) = kteds(id2)
             DO i2 = 1, ncase
               mvals(1      ,ip,1+i2) = casetdata(id1,1,i2)
               mvals(ike+inc,ip,1+i2) = casetdata(id2,1,i2)
             ENDDO
             IF (langdat) THEN
               mvals(1      ,ip,ngs) = targdat(id1,1)
               mvals(ike+inc,ip,ngs) = targdat(id2,1)
             ENDIF
           ELSEIF (ytype.EQ.24) THEN
             mvals(1      ,ip,1) = ktids(id1)
             mvals(ike+inc,ip,1) = ktids(id2)
             DO i2 = 1, ncase
               mvals(1      ,ip,1+i2) = casetdata(id1,2,i2)
               mvals(ike+inc,ip,1+i2) = casetdata(id2,2,i2)
             ENDDO
             IF (langdat) THEN
               mvals(1      ,ip,ngs) = targdat(id1,2)
               mvals(ike+inc,ip,ngs) = targdat(id2,2)
             ENDIF
           ELSEIF (ytype.EQ.25) THEN
             mvals(1      ,ip,1) = knds(id1)
             mvals(ike+inc,ip,1) = knds(id2)
             DO i2 = 1, ncase
               mvals(1      ,ip,1+i2) = casetdata(id1,3,i2)
               mvals(ike+inc,ip,1+i2) = casetdata(id2,3,i2)
             ENDDO
             IF (langdat) THEN
               mvals(1      ,ip,ngs) = targdat(id1,3)
               mvals(ike+inc,ip,ngs) = targdat(id2,3)
             ENDIF
           ELSEIF (ytype.GE.36.AND.ytype.LE.45) THEN
             mvals(1      ,ip,1) = LO
             mvals(ike+inc,ip,1) = LO
             mvals(1      ,ip,2) = LO
             mvals(ike+inc,ip,2) = LO

           ELSEIF (ytype.EQ.52) THEN

             mvals(1      ,ip,1) = 
     .       CalcPressure(knds(id1),kteds(id1),ktids(id1),kvds(id1))*ECH
             mvals(ike+inc,ip,1) = 
     .       CalcPressure(knds(id2),kteds(id2),ktids(id2),kvds(id2))*ECH
             IF (.FALSE..AND.langdat) THEN
               mvals(1      ,ip,ngs) = targdat(id1,1)
               mvals(ike+inc,ip,ngs) = targdat(id2,1)
             ENDIF
           ENDIF
         ELSE
           CALL ER('964','Invalid boundary data type',*9997)
         endif

         pnks(ip) = ike + inc

c         WRITE(6,*) 'PNKS:',ip,pnks(ip)

         DO ik = iks, ike
c
c          Assign x-axis data:
c
           IF     (xtype.EQ.1) THEN
             mOUTS(IK+in,ip) = KSS(IK,IR)
             mWIDS(IK+in,ip) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
           ELSEIF (xtype.EQ.2) THEN
             mOUTS(IK+in,ip) = thetag(IK,IR)
             mWIDS(IK+in,ip) = 0.0
             IF (IK.GT.1) mWIDS(IK+in,ip)=0.5*(thetag(IK,IR)-
     .                                         thetag(IK-1,IR))
             IF (IK.LT.NKS(IR)) mWIDS(IK+in,ip) = mWIDS(IK+in,ip) +
     >                          0.5 * (thetag(IK+1,IR)-thetag(IK,IR))
           ELSEIF (xtype.EQ.3) THEN
             mOUTS(IK+in,ip) = KPS(IK,IR)
             mWIDS(IK+in,ip) = 0.0
             IF (IK.GT.1) mWIDS(IK+in,ip)=0.5*(KPS(IK,IR)-KPS(IK-1,IR))
             IF (IK.LT.NKS(IR)) mWIDS(IK+in,ip) = mWIDS(IK+in,ip) +
     .                              0.5 * (KPS(IK+1,IR)-KPS(IK,IR))
           ELSEIF (xtype.EQ.4) THEN
             mOUTS(IK+in,ip) = KSS(IK,IR) / ksmaxs(ir)
             mWIDS(IK+in,ip) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR)) /
     .                         ksmaxs(ir)
           ELSEIF (xtype.EQ.5) THEN
             mouts(ik+in,ip) = xdata(ik,ir)
             mwids(ik+in,ip) = 0.0
           ELSEIF (xtype.EQ.6.OR.xtype.EQ.8) THEN
c             IF (ik.EQ.iks) THEN
c               i1 = 0
c               DO i2 = 1, ngdata
c                 CALL LoadArray(mouts(1+in,ip),i1,gxdata(1,ir,ind,i2),1,
c     .                          gndata(ir,ind,i2))
c               ENDDO
c               DO i2 = 1, ncase
c                 CALL LoadArray(mouts(1+in,ip),i1,casedata(1,ir,0,i2),1,
c     .                          casendata(ir,i2))
c               ENDDO
c               CALL LoadArray(mouts(1+in,ip),i1,kss(1,ir),1,nks(ir))
c             ENDIF
           ELSEIF (xtype.EQ.9) THEN
             IF (ik.EQ.iks) THEN
               i1 = 0
               DO i2 = 1, ngdata
c...             Invert s-data:
                 DO i3 = 1, gndata(ir,ind,i2)
                   gxdata(i3,ir,ind,i2) =ksmaxs(ir)-gxdata(i3,ir,ind,i2)
                 ENDDO
                 CALL LoadArray(mouts(1+in,ip),i1,gxdata(1,ir,ind,i2),1,
     .                          gndata(ir,ind,i2))
               ENDDO
c...           KSS already inverted:
               CALL LoadArray(mouts(1+in,ip),i1,kss(1,ir),1,nks(ir))
             ENDIF
           ELSEIF (xtype.EQ.7) THEN
             mouts(ik+in,ip) = xdata(ik,ip)
             mwids(ik+in,ip) = 0.0
           ENDIF
c
c          Assign series data:
c
           IF     (ytype.EQ.1) THEN
             MVALS(IK+in,ip,1) = knbs (IK,IR)
             MVALS(IK+in,ip,2) = ktebs(IK,IR)
             MVALS(IK+in,ip,3) = ktibs(IK,IR)
             MVALS(ik+in,ip,4) = kvhs (ik,ir) / qt
             MVALS(ik+in,ip,5) = kvhs (ik,ir) / qt * knbs(ik,ir)
           ELSEIF (ytype.EQ.2) THEN
             MVALS(IK+in,ip,1) = ktebs(IK,IR)
             MVALS(IK+in,ip,2) = ktibs(IK,IR)
           ELSEIF (ytype.EQ.3) THEN
             MVALS(ik+in,ip,1) = ktebs(ik,ir)
             MVALS(ik+in,ip,2) = ktibs(ik,ir)
             MVALS(ik+in,ip,3) = kvhs (ik,ir) / qt * knbs(ik,ir)

             MVALS(ik+in,ip,1) = rs(ik,ir)
             MVALS(ik+in,ip,2) = zs(ik,ir)

           ELSEIF (ytype.EQ.4) THEN
             MVALS(ik+in,ip,1) = ktebs   (ik,ir)
             MVALS(ik+in,ip,2) = knbs    (ik,ir)
             MVALS(ik+in,ip,3) = pinalpha(ik,ir)
           ELSEIF (ytype.EQ.5) THEN
             MVALS(ik+in,ip,1) = ktebs(ik,ir)
             MVALS(ik+in,ip,2) = knbs (ik,ir)
             MVALS(ik+in,ip,3) = knbs (ik,ir) *
     .         (ktebs(ik,ir) + ktibs(ik,ir) + fact *
     .          ((kvhs(ik,ir) / qt)**2.0))
c             MVALS(ik+in,ip,4) = ktibs(ik,ir)
c             MVALS(ik+in,ip,4) = ktebs(ik,ir) + ktibs(ik,ir)
c             MVALS(ik+in,ip,4) = (kvhs(ik,ir) / qt)**2.0
           ELSEIF (ytype.EQ.6) THEN
             MVALS(ik+in,ip,1) = pinrec(ik,ir)
             MVALS(ik+in,ip,2) = pinion(ik,ir)
c             MVALS(ik+in,ip,1) = knbs  (ik,ir)
c             MVALS(ik+in,ip,2) = pinrec(ik,ir)
c             MVALS(ik+in,ip,3) = pinion(ik,ir)
           ELSEIF (ytype.EQ.7) THEN
             IF (ik.EQ.iks) THEN
               CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,1),1,ike,
     .                       kss(1,ir),osmtmp(1,ir),1,nks(ir))
               DO i2 = 1, ncase
                 CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,1+i2),1,ike,
     .                         kss(1,ir),casedata(1,ir,7,i2),1,nks(ir))
               ENDDO
             ENDIF
c             IF (pinqi (1,1).NE.LO) MVALS(ik+in,ip,1) = pinqi (ik,ir)
c             IF (osmpei(1,1).NE.LO) MVALS(ik+in,ip,2) = osmpei(ik,ir)
c             IF (osmcfi(1,1).NE.LO) MVALS(ik+in,ip,3) = osmcfi(ik,ir)
c             MVALS(ik+in,ip,4) = MVALS(ik+in,ip,1) +
c     .                           MVALS(ik+in,ip,2) +
c     .                           MVALS(ik+in,ip,3)
           ELSEIF (ytype.EQ.8) THEN
             IF (pinqe (1,1).NE.LO) MVALS(ik+in,ip,1) = pinqe (ik,ir)
             IF (pinqe (1,1).NE.LO) MVALS(ik+in,ip,2) = pinqi (ik,ir)
             IF (pinqe (1,1).NE.LO) MVALS(ik+in,ip,3) = pinqer(ik,ir)
             IF (pinqe (1,1).NE.LO) MVALS(ik+in,ip,4) = pinqir(ik,ir)

             MVALS(ik+in,ip,5) = MVALS(ik+in,ip,1) +
     .                           MVALS(ik+in,ip,2) +
     .                           MVALS(ik+in,ip,3) +
     .                           MVALS(ik+in,ip,4)

           ELSEIF (ytype.EQ.9) THEN
             IF (pinion(1,1).NE.LO) MVALS(ik+in,ip,1) = pinion(ik,ir)
             IF (pinrec(1,1).NE.LO) MVALS(ik+in,ip,2) =-pinrec(ik,ir)
             IF (osmcfp(1,1).NE.LO) MVALS(ik+in,ip,3) = osmcfp(ik,ir)

c             MVALS(ik+in,ip,4) = MVALS(ik+in,ip,1) +
c     .                           MVALS(ik+in,ip,2) +
c     .                           MVALS(ik+in,ip,3)

           ELSEIF (ytype.EQ.10) THEN
             mvals(ik+in,ip,1) = pinmp(ik,ir)
c...not really accurate...
             mvals(ik+in,ip,2) = -2.0 * AMU * kvhs(ik,ir) / qt *
     .                            pinrec(ik,ir)
c             mvals(ik+in,ip,2) = osmmp(ik,ir) - pinmp(ik,ir)
             mvals(ik+in,ip,3) = osmmp(ik,ir)
             mvals(ik+in,ip,4) = CalcPressure(knbs(ik,ir),ktebs(ik,ir),
     .                                 ktibs(ik,ir),kvhs(ik,ir)/qt)*ECH
c             mvals(ik+in,ip,4) = 0.1 *
c     .         CalcPressure(knbs (ik,ir),ktebs(ik,ir),
c     .                      ktibs(ik,ir),kvhs (ik,ir)/qt) * ECH


c...this should be iksyms, but it is not defined yet
             ik1 = ikmids(ir)
             pm = CalcPressure(knbs (ik1,ir),ktebs(ik1,ir),
     .                         ktibs(ik1,ir),kvhs (ik1,ir)/qt) * ECH

             IF (ik.LT.ik1) THEN
               CALL CalcIntegral3(osmmp,ik,ik1,ir,pl,2)

               mvals(ik+in,ip,5) = pm - pl
             ELSE
               CALL CalcIntegral3(osmmp,ik1,ik,ir,pl,2)

               mvals(ik+in,ip,5) = pm + pl
             ENDIF




             mvals(ik+in,ip,6) = 0.1 * kvhs(ik,ir) / qt
c             mvals(ik+in,ip,5) = 0.01 * kvhs(ik,ir) / qt


           ELSEIF (ytype.EQ.11) THEN
             mvals(ik+in,ip,1) = pinmp(ik,ir)

             IF     (ik.EQ.1) THEN
               rdum1 = kvds(idds(ir,2))
             ELSE
               rdum1 = (kvhs(ik-1,ir)/qt * (kss(ik  ,ir)-ksb(ik-1,ir))+
     .                  kvhs(ik  ,ir)/qt * (ksb(ik-1,ir)-kss(ik-1,ir)))/
     .                 (kss(ik,ir) - kss(ik-1,ir))
             ENDIF

             IF (ik.EQ.nks(ir)) THEN
               rdum2 = kvds(idds(ir,1))
             ELSE
               rdum2 = (kvhs(ik  ,ir)/qt * (kss(ik+1,ir)-ksb(ik,ir)) +
     .                  kvhs(ik+1,ir)/qt * (ksb(ik  ,ir)-kss(ik,ir))) /
     .                 (kss(ik+1,ir) - kss(ik,ir))
             ENDIF

             rdum3 = (rdum2 - rdum1) / (ksb(ik,ir) - ksb(ik-1,ir))

             WRITE(6,*) 'TESTING : ',ik,ir,rdum1,rdum2,rdum3

             tauii = 2.5 * 2.09E+13 * ktibs(ik,ir)**1.5 * SQRT(crmb) /
     .               (knbs(ik,ir) * 15.0)

             pii   = -4.0 / 9.0 * knbs(ik,ir) * ECH * ktibs(ik,ir) *
     .               tauii * rdum3

             mvals(ik+in,ip,2) = pii

           ELSEIF (ytype.EQ.12) THEN

             STOP 'STOP 966: PINDATA ARRAY NOT AVAILABLE'

             mvals(ik+in,ip,1) = pindata(ik,ir,H_MP2) +
     .                           pindata(ik,ir,H_MP3)
             mvals(ik+in,ip,2) = pindata(ik,ir,H_MP4)
             mvals(ik+in,ip,3) = pindata(ik,ir,H_MP9)
             mvals(ik+in,ip,4) = pindata(ik,ir,H_MP2) +
     .                           pindata(ik,ir,H_MP3) +
     .                           pindata(ik,ir,H_MP4)
             mvals(ik+in,ip,5) = pinmp (ik,ir)

           ELSEIF (ytype.EQ.13.OR.ytype.EQ.14) THEN
             mvals(ik+in,ip,1) = ydata(ik,ir)
           ELSEIF (ytype.EQ.15) THEN
             IF (ik.EQ.iks) THEN
               CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,1),1,ike,
     .                       kss(1,ir),ktebs (1,ir),1,nks(ir))
               CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,2),1,ike,
     .                       kss(1,ir),ktibs (1,ir),1,nks(ir))
               DO i2 = 1, ngdata
                 CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,2+i2),1,ike,
     .                         gxdata(1,ir,ind,i2),gydata(1,ir,ind,i2),
     .                         1,gndata(ir,ind,i2))
               ENDDO
             ENDIF
           ELSEIF (ytype.EQ.16) THEN
             IF (ik.EQ.iks) THEN
               CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,1),1,ike,
     .                       kss(1,ir),knbs (1,ir),1,nks(ir))
               DO i2 = 1, ngdata
                 CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,1+i2),1,ike,
     .                         gxdata(1,ir,ind,i2),gydata(1,ir,ind,i2),
     .                         1,gndata(ir,ind,i2))
               ENDDO
             ENDIF
           ELSEIF (ytype.EQ.17) THEN
             mvals(ik+in,ip,1) = knbs(ik,ir) * 1.0E-06 *
     .                         GetEAD(ktebs(ik,ir),knbs(ik,ir),1,'H.4 ')

           ELSEIF (ytype.EQ.18) THEN
             mvals(ik+in,ip,1) = kvhs(ik,ir) / qt /
     .                           GetCs(ktebs(ik,ir),ktibs(ik,ir))
           ELSEIF (ytype.EQ.19) THEN
             mvals(ik+in,ip,1) = kvhs(ik,ir) / qt * knbs(ik,ir)
           ELSEIF (ytype.EQ.20) THEN
             IF (ik.EQ.iks) THEN
               CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,1),1,ike,
     .                       kss(1,ir),ktebs  (1,ir),1,nks(ir))
               in2 = 1
c
               IF (plotuedge) THEN
                 in2 = in2+1
                 CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,in2),1,ike,
     .                         kss(1,ir),e2dtebs(1,ir),1,nks(ir))
               ENDIF
c
c jdemod       
c
               IF (plotav) THEN
                 in2 = in2+1
                 CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,in2),1,ike,
     .                         kss(1,ir),avdata(1,ir),1,nks(ir))
               ENDIF
c
c jdemod
c

               DO i2 = 1, ngdata
                 CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,in2+i2),
     .                         1,ike,
     .                         gxdata(1,ir,ind,i2),gydata(1,ir,ind,i2),
     .                         1,gndata(ir,ind,i2))
               ENDDO

             ENDIF
c
c jdemod
c
           ELSEIF (ytype.EQ.46) THEN
             IF (ik.EQ.iks) THEN
c 
c              OSM
c 
               CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,1),1,ike,
     .                       kss(1,ir),ktebs  (1,ir),1,nks(ir))
               in2 = 1
c
c              UEDGE 
c
               IF (plotuedge) THEN
                 in2 = in2+1
                 CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,in2),1,ike,
     .                         kss(1,ir),e2dtebs(1,ir),1,nks(ir))
               ENDIF
c
c              Thomson Average
c
               in2 = in2+1
               CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,in2),1,ike,
     .                         kss(1,ir),avdata(1,ir),1,nks(ir))
c
             endif
c
c jdemod
c
           ELSEIF (ytype.EQ.21) THEN
             IF (ik.EQ.iks) THEN
               CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,1),1,ike,
     .                       kss(1,ir),ktibs  (1,ir),1,nks(ir))
               CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,2),1,ike,
     .                       kss(1,ir),e2dtibs(1,ir),1,nks(ir))
             ENDIF
           ELSEIF (ytype.EQ.22) THEN
             IF (ik.EQ.iks) THEN
               CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,1),1,ike,
     .                       kss(1,ir),knbs  (1,ir),1,nks(ir))
               in2 = 1

               IF (plotuedge) THEN
                 in2 = in2+1
                 CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,in2),1,ike,
     .                         kss(1,ir),e2dnbs(1,ir),1,nks(ir))
               ENDIF
c
c jdemod       
c
               IF (plotav) THEN
                 in2 = in2+1
                 CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,in2),1,ike,
     .                         kss(1,ir),avdata(1,ir),1,nks(ir))
               ENDIF
c
c jdemod
c

               DO i2 = 1, ngdata
                 CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,in2+i2),
     .                         1,ike,
     .                         gxdata(1,ir,ind,i2),gydata(1,ir,ind,i2),
     .                         1,gndata(ir,ind,i2))
               ENDDO
             ENDIF

           ELSEIF (ytype.EQ.22) THEN
             IF (ik.EQ.iks) THEN
               CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,1),1,ike,
     .                       kss(1,ir),knbs  (1,ir),1,nks(ir))
               in2 = 1

               IF (plotuedge) THEN
                 in2 = in2+1
                 CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,in2),1,ike,
     .                         kss(1,ir),e2dnbs(1,ir),1,nks(ir))
               ENDIF
c
c jdemod       
c
               IF (plotav) THEN
                 in2 = in2+1
                 CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,in2),1,ike,
     .                         kss(1,ir),avdata(1,ir),1,nks(ir))
               ENDIF
c
c jdemod
c

               DO i2 = 1, ngdata
                 CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,in2+i2),
     .                         1,ike,
     .                         gxdata(1,ir,ind,i2),gydata(1,ir,ind,i2),
     .                         1,gndata(ir,ind,i2))
               ENDDO
             ENDIF

c
c jdemod
c
           ELSEIF (ytype.EQ.47) THEN
             IF (ik.EQ.iks) THEN
c
c              OSM
c
               CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,1),1,ike,
     .                       kss(1,ir),knbs  (1,ir),1,nks(ir))
               in2 = 1

c
c              UEDGE
c
               IF (plotuedge) THEN
                 in2 = in2+1
                 CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,in2),1,ike,
     .                         kss(1,ir),e2dnbs(1,ir),1,nks(ir))
               ENDIF
c
c              Average Thomson
c
               in2 = in2+1
               CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,in2),1,ike,
     .                         kss(1,ir),avdata(1,ir),1,nks(ir))
             ENDIF
c
c jdemod
c
           ELSEIF (ytype.EQ.23) THEN
             IF (ik.EQ.iks) THEN
               CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,1),1,ike,
     .                       kss(1,ir),ktebs(1,ir),1,nks(ir))
               DO i2 = 1, ncase
                CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,1+i2),1,ike,
     .                        casedata(1,ir,0,i2),casedata(1,ir,1,i2),
     .                        1,casendata(ir,i2))
               ENDDO
               i3 = 1 + ncase 
               DO i2 = 1, ngdata
                CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,i3+i2),1,ike,
     .                        gxdata(1,ir,ind,i2),gydata(1,ir,ind,i2),
     .                        1,gndata(ir,ind,i2))
               ENDDO
             ENDIF
           ELSEIF (ytype.EQ.24) THEN
             IF (ik.EQ.iks) THEN
               CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,1),1,ike,
     .                       kss(1,ir),ktibs(1,ir),1,nks(ir))
               DO i2 = 1, ncase
                 CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,1+i2),1,ike,
     .                         kss(1,ir),casedata(1,ir,2,i2),1,nks(ir))
               ENDDO
               i3 = 1 + ncase 
               DO i2 = 1, ngdata
                CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,i3+i2),1,ike,
     .                        gxdata(1,ir,ind,i2),gydata(1,ir,ind,i2),
     .                        1,gndata(ir,ind,i2))
               ENDDO
             ENDIF
           ELSEIF (ytype.EQ.25) THEN
             IF (ik.EQ.iks) THEN
               CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,1),1,ike,
     .                       kss(1,ir),knbs   (1,ir),1,nks(ir))
               DO i2 = 1, ncase
                 CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,1+i2),1,ike,
     .                         kss(1,ir),casedata(1,ir,3,i2),1,nks(ir))
               ENDDO
               i3 = 1 + ncase 
               DO i2 = 1, ngdata
                CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,i3+i2),1,ike,
     .                        gxdata(1,ir,ind,i2),gydata(1,ir,ind,i2),
     .                        1,gndata(ir,ind,i2))
               ENDDO
             ENDIF
           ELSEIF (ytype.EQ.26) THEN
             IF (ik.EQ.iks) THEN
               CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,1),1,ike,
     .                       kss(1,ir),pinalpha(1,ir),1,nks(ir))
               DO i2 = 1, ncase
                 CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,1+i2),1,ike,
     .                         kss(1,ir),casedata(1,ir,4,i2),1,nks(ir))
               ENDDO
             ENDIF
           ELSEIF (ytype.EQ.27) THEN
             IF (ik.EQ.iks) THEN
               CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,1),1,ike,
     .                       kss(1,ir),pinatom(1,ir),1,nks(ir))
               DO i2 = 1, ncase
                 CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,1+i2),1,ike,
     .                         kss(1,ir),casedata(1,ir,5,i2),1,nks(ir))
               ENDDO
             ENDIF
           ELSEIF (ytype.EQ.28) THEN
             IF (ik.EQ.iks) THEN
               CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,1),1,ike,
     .                       kss(1,ir),pinmol(1,ir),1,nks(ir))
               DO i2 = 1, ncase
                 CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,1+i2),1,ike,
     .                         kss(1,ir),casedata(1,ir,6,i2),1,nks(ir))
               ENDDO
             ENDIF
           ELSEIF (ytype.EQ.29) THEN
             MVALS(IK+in,ip,1) = bratio(IK,IR)
           ELSEIF (ytype.EQ.30.OR.ytype.EQ.51) THEN
             MVALS(IK+in,ip,1) = ydata(ik,ip)
           ELSEIF (ytype.GE.36.AND.ytype.LE.45) THEN
             MVALS(IK+in,ip,1) = pinalgv(ik,ir,ytype-36+1)
c *TEMP*
             IF (ytype.EQ.39) THEN
c...           D2+ density:
               ti = ktibs(ik,ir)
               ne = knbs (ik,ir)
               sigmav = GetEAD(ti,ne,22,'H.12')
               nD2p = sigmav * pinmol(ik,ir)

c...           MAR: e + H2+ -> e + H + H
               sigmav1 = GetEAD(ti,ne,17,'H.4 ')
c...           MAD: e + H2+ -> e + H + H+
               sigmav2 = GetEAD(ti,ne,15,'H.4 ')

               sigmav = sigmav1**2 / (sigmav1 + sigmav2)

               MVALS(IK+in,ip,2) = sigmav * 1.0E-06 * ne * nD2p * 2.0
             ENDIF

           ELSEIF (ytype.EQ.48) THEN
             mvals(ik+in,ip,1) = pinioncomp(ik,ir,1)
             mvals(ik+in,ip,2) = pinioncomp(ik,ir,2)
             mvals(ik+in,ip,3) = pinioncomp(ik,ir,3)
             mvals(ik+in,ip,4) = pinioncomp(ik,ir,1) + 
     .                           pinioncomp(ik,ir,2) +  
     .                           pinioncomp(ik,ir,3)
             mvals(ik+in,ip,5) = pinion(ik,ir)

           ELSEIF (ytype.EQ.49) THEN

             mvals(ik+in,ip,1) = pinline(ik,ir,6,line)
             IF (line.EQ.2) 
     .         mvals(ik+in,ip,1) = mvals(ik+in,ip,1)*
     .                    (6.63E-34*3.0E+08)/(4340.0*1.0E-10)

             mvals(ik+in,ip,2) = osmtmp(ik,ir) * dmax / gmax

           ELSEIF (ytype.EQ.50) THEN

             fact = GetEAD(ktebs(ik,ir),knbs(ik,ir),1,'H.4 ')
             mvals(ik+in,ip,1) = knbs(ik,ir) * fact

             DO i2 = 1, ncase
               plottype(1+i2) = i2 + 1
               fact = GetEAD(casedata(ik,ir,1,i2),
     .                       casedata(ik,ir,3,i2),1,'H.4 ')
               mvals(ik+in,ip,1+i2) = casedata(ik,ir,3,i2) * fact
c               WRITE(0,*) 'TE,NE:',casedata(ik,ir,1,i2),
c     .                             casedata(ik,ir,3,i2),fact
             ENDDO

           ELSEIF (ytype.EQ.52) THEN
             IF (ik.EQ.iks) THEN

               osmtmp = 0.0
               DO ik1 = 1, nks(ir)
                 osmtmp(ik1,ir) = 
     .             CalcPressure(knbs (ik1,ir),ktebs(ik1,ir),
     .                          ktebs(ik1,ir),0.0)*ECH
               ENDDO

               CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,1),1,ike,
     .                       kss(1,ir),osmtmp(1,ir),1,nks(ir))

               in2 = 1
               IF (plotuedge) THEN
                 in2 = in2+1
                 CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,in2),1,ike,
     .                         kss(1,ir),e2dtebs(1,ir),1,nks(ir))
               ENDIF
               DO i2 = 1, ngdata
                 CALL MapArray(mouts(1+in,ip),mvals(1+in,ip,in2+i2),
     .                         1,ike,
     .                         gxdata(1,ir,ind,i2),gydata(1,ir,ind,i2),
     .                         1,gndata(ir,ind,i2))
               ENDDO

             ENDIF

           ELSEIF (ytype.EQ.53) THEN
c...         Electric field:
             mvals(ik+in,ip,1) = kes(ik,ir)

           ELSEIF (ytype.EQ.54) THEN
c...         Electric potential:
             mvals(ik+in,ip,1) = osmpot(ik,ir)

           ENDIF

         ENDDO
c
c
c
         DO i1 = 1, ngs
           status = .FALSE.
           DO ik = 1, ike+inc
             IF (mvals(ik,ip,i1).NE.0.0) status = .TRUE.
           ENDDO
           IF (.NOT.status) THEN
             WRITE(6,*) '966: NULL DATA ARRAY DETECTED: IP,I1 = ',ip,i1
             CALL RSet(mvals(1,ip,i1),MAXNPS,LO)
           ENDIF
         ENDDO
c
c        Normalize plots:
c
         DO i1 = 1, ngs
           inorm(i1) = i1
         ENDDO
c
c        Special normalization for temperatures:
c
         IF ((.FALSE..AND.ytype.EQ.1).OR.
     .       (.FALSE..AND.ytype.EQ.2).OR.
     .       (.FALSE..AND.ytype.EQ.3).OR.
     .       ytype.EQ.6
     .       ) THEN
           IF     (ytype.EQ.1) THEN
             i1 = 2
             i2 = 3
           ELSEIF (ytype.EQ.2) THEN
             i1 = 1
             i2 = 2
           ELSEIF (ytype.EQ.3) THEN
             i1 = 1
             i2 = 2
           ELSEIF (ytype.EQ.6) THEN
             i1 = 2
             i2 = 3
           ENDIF

           max1 = LO
           max2 = LO

           DO ik = 1, ike+inc
             max1 = MAX(max1,ABS(mvals(ik,ip,i1)))
             max2 = MAX(max2,ABS(mvals(ik,ip,i2)))
           ENDDO

           IF (max1.GT.max2) THEN
             inorm(i1) = i1
             inorm(i2) = i1
           ELSE 
             inorm(i1) = i2
             inorm(i2) = i2
           ENDIF

c           WRITE(6,*) max1,max2,i1,i2,inorm(i1),inorm(i2)
         ENDIF

         IF (.FALSE..AND.
     .       (ytype.EQ.1.OR.ytype.EQ.3.OR.ytype.EQ.4.OR.
     .        ytype.EQ.5.OR.ytype.EQ.6)) THEN

           DO i1 = 1, ngs
             ynorm(i1) = LO
             DO ik = 1, ike+inc
               ynorm(i1) = MAX(ynorm(i1),ABS(mvals(ik,ip,inorm(i1))))
             ENDDO
           ENDDO

           DO i1 = 1, ngs
             DO ik = 1, ike+inc
               mvals(ik,ip,i1) = mvals(ik,ip,i1) / ynorm(i1)
             ENDDO
           ENDDO
         ENDIF
c
c
c
         write(6,'(a)') 'PLOT DATA:'
c

         do i1 = 1,ngs
c
            if (i1.eq.1.or.i1.eq.2.or.i1.eq.3.or.i1.eq.ngs) 
     >         write (6,'(a,i5)') 'DATASET:'//elabs(i1),i1
c
            do ik = iks, ike+inc
c               if (mvals(ik,ip,i1).gt.LO) then
                  write(6,'(2i5,2(1x,f12.5),1p,3(1x,e12.5))') ik,ir,
     >                  mouts(ik,ip),mouts(ike+inc,ip)-mouts(ik,ip),
     >                  mvals(ik,ip,i1)
c               endif
            end do
         end do




       ENDDO



c...   Output:
       DO ip = 1, nplts

         ir = ringnos(ip)
         IF (ytype.EQ.15.OR.ytype.EQ.16.OR.ytype.EQ.20.OR.
     .       ytype.EQ.21.OR.ytype.EQ.22.or.ytype.eq.46.or.
     .       ytype.eq.47) THEN
           iks = 1
           ike = nks(ir)
           DO i2 = 1, ngdata
             ike = ike + gndata(ir,ind,i2)
           ENDDO
         ELSEIF (ytype.EQ.13.OR.ytype.EQ.14) THEN
           iks = 1
           ike = 30
         ELSEIF (ytype.EQ.30.OR.ytype.EQ.51) THEN
           iks = 1
           ike = pnks(ip)
         ELSE
           iks = 1
           ike = nks(ir)
         ENDIF

c         DO ik = 1, nks(ir)
c           WRITE(0,*) '---->',ik,ir,kpb(ik,ir)
c         ENDDO 
c         STOP 'checking KPB'

         IF (xtype.EQ.8) THEN
c...       Convert from s-space to p-space:
           DO ik = iks, ike+inc
             convert = .FALSE.
             DO ik1 = 1, nks(ir)
               IF (mouts(ik,ip).GE.ksb(ik1-1,ir).AND.
     .             mouts(ik,ip).LE.ksb(ik1  ,ir)) THEN
                 frac = (mouts(ik,ip) - ksb(ik1-1,ir)) /
     .                  (ksb(ik1,ir)  - ksb(ik1-1,ir))
                 mouts(ik,ip)=(1.0-frac)*kpb(ik1-1,ir)+frac*kpb(ik1,ir)
                 convert = .TRUE.
                 EXIT
               ENDIF
             ENDDO
             IF (.NOT.convert) 
     .         CALL ER('966','Unable to convert from s to s-pol',*99)
           ENDDO
         ENDIF

       ENDDO
c
c...  Change plot labels:    
      CALL CustomisePlot(title,xlab,ylab,elabs)

      READ(5,'(A256)') dummy
      IF   (dummy(8:11).EQ.'Xran'.OR.dummy(8:11).EQ.'xran'.OR.
     .      dummy(8:11).EQ.'XRAN') THEN
        READ(dummy,*) cdum1,xrange1,xrange2
        DO ip = 1, nplts
c...      Shift data points to make room for x-axis range points:
          DO ik = pnks(ip)+1, 2, -1
            mouts(ik,ip) = mouts(ik-1,ip)
            DO i1 = 1, ngs
              mvals(ik,ip,i1) = mvals(ik-1,ip,i1)
            ENDDO
          ENDDO
          pnks(ip) = pnks(ip) + 2
c...      Add new data points that mark x-axis range:
          mouts(1       ,ip) = xrange1
          mouts(pnks(ip),ip) = xrange2
          DO i1 = 1, ngs
            mvals(1       ,ip,i1) = LO
            mvals(pnks(ip),ip,i1) = LO
          ENDDO          
        ENDDO
      ELSE
        BACKSPACE 5
      ENDIF



c...  Analysis of DTS data:

95    READ(5,'(A80)') graph1
      IF (graph1(8:11).EQ.'Wave'.OR.graph1(8:11).EQ.'WAVE'.OR.
     .    graph1(8:11).EQ.'wave') THEN

        WRITE(0,*) 'WAVING!'

        ind = 2

        READ(graph1,*) cdum1,sbnd1,sbnd2

       
c....   Assign ne data:         
c        CALL AsgnThomsonData(thnnraw,thnraw,1,
c     .                      gxdata(1,1,1,ngdata),gydata(1,1,1,ngdata),
c     .                      gndata(  1,1,ngdata),
c     .                      MAXNPS,MAXTYP,MAXTDAT,MAXCOLS)


        WRITE(0,*) pnks(ind),ip,ir,ngdata
        DO ir = irsep, irsep
          WRITE(0,*) 'GNL',gndata(ir,1,1),gndata(ir,1,2),gndata(ir,1,3)
          WRITE(0,*) '   ',gndata(ir,2,1),gndata(ir,2,2),gndata(ir,2,3)
          WRITE(0,*) '   ',gndata(ir,3,1),gndata(ir,3,2),gndata(ir,3,3)
          WRITE(0,*) '   ',gndata(ir,4,1),gndata(ir,4,2),gndata(ir,4,3)
        ENDDO


c...    Overwrite data for IND:
        pnks(ind) = 10        
        DO i1 = 1, 10
          mouts(i1,ind)   = REAL(i1)        
          mvals(i1,ind,1) = REAL(i1)        
          mvals(i1,ind,2) = REAL(i1) + 1.0        
        ENDDO

        ip = 1
        ir = ringnos(ip)
   
     

        DO i2 = 1, ngs
          DO i1 = 1, pnks(1)
            mvals(i1,ind,i2) = 0.0
          ENDDO
        ENDDO



        pnks(ind) = 0

        pnks(ind) = pnks(ind) + 1
        mouts(pnks(ind),ind)   = 0.0
        mvals(pnks(ind),ind,2) = 1.0
        mvals(pnks(ind),ind,1) = 0.0     

        DO i2 = 1, ngdata
          DO i1 = 1, gndata(ir,1,i2)
            pnks(ind) = pnks(ind) + 1
            mouts(pnks(ind),ind)   = gxdata(i1,ir,2,i2)
            mvals(pnks(ind),ind,2) = gydata(i1,ir,2,i2)
            mvals(pnks(ind),ind,1) = 0.0     
c            mvals(pnks(ind),ind,1) = REAL(i1) + 1.0        
          ENDDO
        ENDDO



c        STOP 'dfssd'

c...     Calculate pressure:
c         DO ir = 2, nrs
c           IF (idring(ir).EQ.BOUNDARY) CYCLE
c           gndata(ir,3,ngdata) = gndata(ir,1,ngdata)
c           DO i1 = 1, gndata(ir,3,ngdata)
c             gxdata(i1,ir,3,ngdata) = gxdata(i1,ir,1,ngdata)
c             gydata(i1,ir,3,ngdata) = 
c     .         CalcPressure(gydata(i1,ir,1,ngdata),
c     .                      gydata(i1,ir,2,ngdata),
c     .                      gydata(i1,ir,2,ngdata)*1.0,0.0)*ECH
c           ENDDO
c         ENDDO



         GOTO 95
       ELSE
         BACKSPACE 5
       ENDIF



         write(6,'(a)') 'PLOT DATA:'
c
        DO ip = 1, nplts
c
c jdemod - not sure why but the data set identifier is only being written for the first 3 and the last
c          some of the plots have 10 or more so the tags are missing and the data is grouped 
c          together in the print out
c
c            if (i1.eq.1.or.i1.eq.2.or.i1.eq.3.or.i1.eq.ngs) 
c     >         write (6,'(a,i5)') 'DATASET:'//elabs(i1),i1
c 
c
            write (6,'(a10,18X,20(A10:))') 
     .        'DATASET:',(elabs(i1)(1:10),i1=1,ngs)
c
c jdemod
c
            do ik = 1, pnks(ip)
c               if (mvals(ik,ip,i1).gt.LO) then
                  write(6,'(2i5,2(1x,f8.3),1p,20(E10.2:))') ik,ip,
     >                  mouts(ik,ip),mouts(ike+inc,ip)-mouts(ik,ip),
     >                  (mvals(ik,ip,i1),i1=1,ngs)
c               endif
            end do
       enddo
c


       CALL SLDRAWM (MOUTS,MWIDS,MVALS,MAXNPS,pnks,
     >             nplts,ngs,pltlabs,elabs,xlab,ylab,ref2,title,
     >             sctype,ngrm,pltmins,pltmaxs,pltfact)


      DEALLOCATE(osmtmp)
      DEALLOCATE(casedata)
      DEALLOCATE(casetdata)

      DEALLOCATE(gndata)
      DEALLOCATE(gydata)
      DEALLOCATE(gxdata)

      IF (xtype.EQ.9) THEN
c...    Restore KSS:
        DO ir = 1, nrs
          DO ik = 1, nks(ir)
            kss(ik,ir) = ksmaxs(ir) - kss(ik,ir)
          ENDDO
        ENDDO
      ENDIF

 9997 CONTINUE
      RETURN

 9012 FORMAT(1X,'PLOT',I3,4X,A)
98    CALL ER('966','End of file error',*99)
99    STOP
      END
