c     -*-Fortran-*-
c
c ======================================================================
c
c subroutine: OutputIonisationTimeData
c
      SUBROUTINE OutputIonisationTimeData(fp)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slcom'

      INTEGER fp,i1,i2,i3

      CALL HD(fp,'  TIME-TO-IONISTION STATISTICS','IONTIME-HD',5,72)      

      IF (eirniontime.GT.0) THEN

        WRITE(fp,*)
        WRITE(fp,'(3X,A8,4A8,2X,A10,A12)') 
     .    'REGION','R1(m)','R2(m)','Z1(m)','Z2(m)','COUNT','AVERAGE(s)'
        DO i1 = 1, eirniontime
          WRITE(fp,'(3X,I8,4F8.3,2X,F10.2,1P,E12.2,0P)')
     .      i1,(eiriontime(i1,i2),i2=1,4),(eiriontime(i1,i2),i2=18,19)
        ENDDO

        DO i1 = 1, eirniontime
          IF (eiriontime(i1,5).EQ.0.0) CYCLE

          WRITE(fp,*) 
          WRITE(fp,'(3X,2A8,2A10,2X,A14)')
     .      'REGION','BIN','T1(s)','T2(s)','FRACTION(%)'

          DO i2 = 1, NINT(eiriontime(i1,5))+2

            WRITE(fp,'(3X,2I8,1P,2E10.2,0P,2X,F14.2)')
     .        i1,i2,(eiriontime(i1,20+3*(i2-1)+i3),i3=0,2)
          ENDDO

        ENDDO
      ELSE
        WRITE(fp,*)
        WRITE(fp,*) '    No regions assigned.'
      ENDIF

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: OutputAnalysis
c
      SUBROUTINE OutputAnalysis
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER SymmetryPoint
      REAL    GetCs

      INTEGER ik,ir,ikm,id1,id2,i1
      REAL    rdum1,rdum2,rdum3,rdum4,rdum5,rdum6

      IF (osm_mode.GE.1) THEN

        CALL AnalyseSolution  (PINOUT)
        CALL AnalyseContinuity(PINOUT)

        CALL ShowStats

        DO ir = irsep, nrs
          IF (ikbound(ir,IKLO).EQ.1) THEN
            id1 = idds(ir,2)

            rdum1 = kteds(id1)
            rdum3 = ktids(id1)
            rdum5 = knds (id1)
          ELSE
            id1 = ikbound(ir,IKLO)

            rdum1 = ktebs(id1,ir)
            rdum3 = ktibs(id1,ir)
            rdum5 = LO
            DO ik = id1, SymmetryPoint(ir)
              rdum5 = MAX(rdum5,knbs(ik,ir))
            ENDDO
          ENDIF

          IF (ikbound(ir,IKLO).EQ.nks(ir)) THEN
            id2 = idds(ir,1)

            rdum2 = kteds(id2)
            rdum4 = ktids(id2)
            rdum6 = knds (id2)
          ELSE
            id2 = ikbound(ir,IKHI)

            rdum2 = ktebs(id2,ir)
            rdum4 = ktibs(id2,ir)
            rdum6 = LO
            DO ik = SymmetryPoint(ir)+1, id2
              rdum6 = MAX(rdum6,knbs(ik,ir))
            ENDDO
          ENDIF

          id1 = idds(ir,2)
          id2 = idds(ir,1)

          WRITE(79,'(A,1X,3I4,1X,I4,1P,8E12.4)') '''TARGET    1.00''',
     .      rel_step,rel_iter,rel_count,ir,
     .      rdum1,rdum2,rdum3,rdum4,rdum5,rdum6,
     .      cmachno(ir,2),cmachno(ir,1)
        ENDDO

        DO ir = irsep, nrs
          WRITE(79,'(A,1X,3I4,1X,I4,4F12.6)')    '''RELAX     1.00''',
     .      rel_step,rel_iter,rel_count,ir,
     .      rel_qemul (IKLO,ir),rel_qemul (IKHI,ir),
     .      osm_peimul(IKLO,ir),osm_peimul(IKHI,ir)
        ENDDO


        DO ir = irsep, nrs
          IF (idring(ir).EQ.-1) CYCLE
 
          ikm = SymmetryPoint(ir)

          CALL CalcIntegral3(osmmp,ikbound(ir,IKLO),ikm,ir,rdum1,2)
          CALL CalcIntegral3(osmmp,ikm,ikbound(ir,IKHI),ir,rdum2,2)
          CALL CalcIntegral3(osmmp,1  ,ikm    ,ir,rdum3,2)
          CALL CalcIntegral3(osmmp,ikm,nks(ir),ir,rdum4,2)
          CALL CalcIntegral3(pinrec,1  ,ikm    ,ir,rdum5,2)
          CALL CalcIntegral3(pinrec,ikm,nks(ir),ir,rdum6,2)

          WRITE(79,'(A,1X,3I4,1X,I4,1P,6E12.4)') '''MOM LOSS  1.00''',
     .      rel_step,rel_iter,rel_count,ir,
     .      -rdum1,rdum2,-rdum3,rdum4,rdum5,rdum6
        ENDDO


        DO ir = irsep, nrs
          IF (idring(ir).EQ.-1) CYCLE

          ikm = SymmetryPoint(ir)
          rdum1 = 0.0
          DO ik = ikbound(ir,IKLO), ikm
             rdum1 = MAX(rdum1,
     .                   kvhs(ik,ir)/
     .                   GetCs(ktebs(ik,ir),ktibs(ik,ir)))
          ENDDO
          rdum2 = 0.0
          DO ik = ikm+1, ikbound(ir,IKHI)
             rdum2 = MAX(rdum2,
     .                   kvhs(ik,ir)/
     .                   GetCs(ktebs(ik,ir),ktibs(ik,ir)))
          ENDDO
          WRITE(79,'(A,1X,3I4,1X,I4,2F12.6)') '''MAXMACH   1.00''',
     .      rel_step,rel_iter,rel_count,ir,
     .      rdum1,rdum2
        ENDDO

c...some other quantities for history plots
        WRITE(79,'(A,1X,3I4,1X,1P,8E12.4,0P)') '''REGINTEG  1.00''',
     .      rel_step,rel_iter,rel_count
        CALL OutputRegionIntegrals(pinion)
        CALL OutputRegionIntegrals(pinrec)
        CALL OutputRegionIntegrals(pinatom)
        CALL OutputRegionIntegrals(pinmol)
        CALL OutputRegionIntegrals(pinqe)
        CALL OutputRegionIntegrals(pinqi)
        CALL OutputRegionIntegrals(pinmp)
        CALL OutputRegionIntegrals(pinalpha)
        DO i1 = 1, 5
          CALL OutputRegionIntegrals(pinline(1,1,i1,H_BALPHA))
        ENDDO

c        CALL AnalyseContinuity(PINOUT)

      ENDIF

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: OutputBGKData
c
      SUBROUTINE OutputBGKData
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER fp,ik,ir,i1,i2,in

      fp = 98

      OPEN(UNIT=fp,FILE='bgk.dat',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE')

      WRITE(fp,*) 'BGK data:'
c
c     Ionisation:
c

      DO in = 1, eirnsdtor

        DO i2 = 1, 6
          WRITE(fp,'(2(A,I6),A)') 'Ion species ',i2,
     .                            ' toriodal segment ',in,':'

          DO ir = osm_watch1, osm_watch2
c           IPP/08 Krieger - ensure that idring(ir) is referenced inside
c           array bounds
c           IF (ir.LT.1.OR.ir.GT.nrs.OR.idring(ir).EQ.-1) CYCLE
            IF (ir.LT.1.OR.ir.GT.nrs) THEN
              CYCLE
            ELSEIF (idring(ir).EQ.-1) THEN
              CYCLE
            ENDIF

            WRITE(fp,10) 'ik','ir','KNBS','DIIN','TIIN',
     .                             'VXIN','VYIN','VZIN'
10          FORMAT(2A4,A12,2X,5A12)

            DO ik = 1, nks(ir)
              WRITE(fp,'(2I4,1P,E12.4,2X,5E12.4)')
     .          ik,ir,knbs(ik,ir),(pinbgk(ik,ir,(i2-1)*5+i1),i1=1,5)
            ENDDO

          ENDDO
        ENDDO
      ENDDO

      CLOSE(fp)

      RETURN
99    STOP
      END

c
c ======================================================================
c
c subroutine: OutputRecombinationData
c
      SUBROUTINE OutputRecombinationData
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER   fp,ik,ir,i1,i2
      CHARACTER note*20

      fp = 98

      OPEN(UNIT=fp,FILE='rec.dat',ACCESS='SEQUENTIAL',STATUS='REPLACE')

      WRITE(fp,*) 'Recombination data:'

      DO ir = osm_watch1, osm_watch2
c       IPP/08 Krieger - ensure that idring(ir) is referenced inside
c       array bounds
c       IF (ir.LT.1.OR.ir.GT.nrs.OR.idring(ir).EQ.-1) CYCLE
        IF (ir.LT.1.OR.ir.GT.nrs) THEN
          CYCLE
        ELSEIF (idring(ir).EQ.-1) THEN
          CYCLE
        ENDIF

        WRITE(fp,'(2A4,4A12)') 'ik','ir','pinior','pinmpr','pinqir',
     .                                   'pinqer'

        DO ik = 1, nks(ir)

          note(1:20) = '                   '
          IF (ik.EQ.ikto2 (ir)) note = note(1:LEN_TRIM(note))//' IKTO2'
          IF (ik.EQ.ikti2 (ir)) note = note(1:LEN_TRIM(note))//' IKTI2'
          IF (ik.EQ.ikmids(ir)) note = note(1:LEN_TRIM(note))//' IKMIDS'
          IF (ik.EQ.ikbound(ir,IKLO))
     .      note = note(1:LEN_TRIM(note))//' IK1'
          IF (ik.EQ.ikbound(ir,IKHI))
     .      note = note(1:LEN_TRIM(note))//' IK2'


          WRITE(fp,'(2I4,1P,4E12.4,0P,A)')
     .      ik,ir,pinior(ik,ir),pinmpr(ik,ir),pinqir(ik,ir),
     .            pinqer(ik,ir),
     .      note
        ENDDO
      ENDDO

      CLOSE(fp)

      RETURN
99    STOP
      END
c
c ======================================================================
c
c
c subroutine: OutputAdditionalCellData
c
      SUBROUTINE OutputAdditionalCellData
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER         nradd
      COMMON /RADCOM/ nradd

      INTEGER fp,i1,i2,i3


      fp = 98

      OPEN(UNIT=fp,FILE='acd.dat',ACCESS='SEQUENTIAL',
     .     STATUS='UNKNOWN',POSITION='APPEND')

      WRITE(fp,*)
      WRITE(fp,*) 'ITERATION= ',rel_count,':'
      WRITE(fp,*)

      DO i1 = 1, 10
        IF     (i1.GE.1.AND.i1.LE. 2) THEN
          WRITE(fp,'(A6,3A12,I3)') 'IN','PDENA','EDENA','VOLUME',i1
          DO i2 = 1, nradd
            WRITE(fp,'(I6,1P,3E12.4)')
     .        i2+indasd-1,(pinasd(i2,i3,i1,1),i3=1,2),pinasd(i2,5,i1,1)
          ENDDO
        ELSEIF (i1.EQ.3) THEN
          WRITE(fp,'(A6,2A12,I3)') 'IN','PDENM','EDENM',i1
          DO i2 = 1, nradd
            WRITE(fp,'(I6,1P,2E12.4)')
     .        i2+indasd-1,(pinasd(i2,i3,i1,1),i3=1,2)
          ENDDO
        ELSEIF (i1.EQ.4) THEN
          WRITE(fp,'(A6,2A12,I3)') 'IN','PDENI','EDENI',i1
          DO i2 = 1, nradd
            WRITE(fp,'(I6,1P,2E12.4)')
     .        i2+indasd-1,(pinasd(i2,i3,i1,1),i3=1,2)
          ENDDO
        ELSEIF (i1.GE.5.AND.i1.LE.10) THEN
          WRITE(fp,'(A6,5A12,I3)') 'IN','DIIN','TIIN','VXIN','VYIN',
     .                             'VZIN',i1
          DO i2 = 1, nradd
            WRITE(fp,'(I6,1P,5E12.4)')
     .        i2+indasd-1,(pinasd(i2,i3,i1,1),i3=1,5)
          ENDDO
        ENDIF
      ENDDO



      CLOSE(fp)


c...TEMP
c      DO i1 = 7, 7
c        DO i2 = 1, 10
c          WRITE(0,*)  i1,i2,(pinasd(i2,i3,i1,1),i3=1,2)
c        ENDDO
c      ENDDO


      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: OutputLineRadiationData
c
      SUBROUTINE OutputLineRadiationData
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER fp,ik,ir,i1,i2

      fp = 98

      OPEN(UNIT=fp,FILE='line.dat',ACCESS='SEQUENTIAL',STATUS='REPLACE')

      WRITE(fp,*) 'Line radiation data:'

      DO i1 = 1, 2
        WRITE(fp,*) ' '
        WRITE(fp,*) ' Line: ',i1

        DO ir = osm_watch1, osm_watch2
c         IPP/08 Krieger - ensure that idring(ir) is referenced inside
c         array bounds
C         IF (ir.LT.1.OR.ir.GT.nrs.OR.idring(ir).EQ.BOUNDARY) CYCLE
          IF (ir.LT.1.OR.ir.GT.nrs) THEN
            CYCLE
          ELSEIF (idring(ir).EQ.BOUNDARY) THEN
            CYCLE
          ENDIF

          WRITE(fp,'(2A4,6A12)') 'ik','ir','1','2','3','4','5','Sum'

          DO ik = 1, nks(ir)
            IF (i1.EQ.1) THEN
              WRITE(fp,'(2I4,1P,7E12.4)')
     .          ik,ir,(pinline(ik,ir,i2,i1),i2=1,6),pinalpha(ik,ir)
            ELSE
              WRITE(fp,'(2I4,1P,6E12.4)')
     .          ik,ir,(pinline(ik,ir,i2,i1),i2=1,6)
            ENDIF
          ENDDO
        ENDDO

      ENDDO

      CLOSE(fp)

      RETURN
99    STOP
      END
c
c
c ======================================================================
c
c subroutine: OutputStratumData
c
      SUBROUTINE OutputStratumData
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER fp,ik,ir,i1

      fp = 98

      OPEN(UNIT=fp,FILE='strata.dat',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE')

      WRITE(fp,*) 'Stratum data:'
c
c     Ionisation:
c
      WRITE(fp,*) ' '
      WRITE(fp,*) ' Ionisation: '

      DO ir = osm_watch1, osm_watch2
c       IPP/08 Krieger - ensure that idring(ir) is referenced inside
c       array bounds
c       IF (ir.LT.1.OR.ir.GT.nrs.OR.idring(ir).EQ.-1) CYCLE
        IF (ir.LT.1.OR.ir.GT.nrs) THEN
          CYCLE
        ELSEIF (idring(ir).EQ.-1) THEN
          CYCLE
        ENDIF

        WRITE(fp,'(2A4,6A12)') 'ik','ir','1','2','3','4','SUM','PINION'

        DO ik = 1, nks(ir)
          WRITE(fp,'(2I4,1P,6E12.4)')
     .      ik,ir,
     .     (pindata(ik,ir,i1),i1=H_ION1,H_ION4),
     .      pindata(ik,ir,H_ION1)+pindata(ik,ir,H_ION2)+
     .      pindata(ik,ir,H_ION3)+pindata(ik,ir,H_ION4),
     .      pinion(ik,ir)
        ENDDO
      ENDDO
c
c     H density:
c
      WRITE(fp,*) ' '
      WRITE(fp,*) 'Neutral hydrogen atom density: '

      DO ir = osm_watch1, osm_watch2
c       IPP/08 Krieger - ensure that idring(ir) is referenced inside
c       array bounds
c       IF (ir.LT.1.OR.ir.GT.nrs.OR.idring(ir).EQ.-1) CYCLE
        IF (ir.LT.1.OR.ir.GT.nrs) THEN
          CYCLE
        ELSEIF (idring(ir).EQ.-1) THEN
          CYCLE
        ENDIF

        WRITE(fp,'(2A4,6A12)') 'ik','ir','1','2','3','4','SUM','PINATOM'

        DO ik = 1, nks(ir)
          WRITE(fp,'(2I4,1P,6E12.4)')
     .      ik,ir,
     .     (pindata(ik,ir,i1),i1=H_ATM1,H_ATM4),
     .      pindata(ik,ir,H_ATM1)+pindata(ik,ir,H_ATM2)+
     .      pindata(ik,ir,H_ATM3)+pindata(ik,ir,H_ATM4),
     .      pinatom(ik,ir)
        ENDDO
      ENDDO
c
c     H2 density:
c
      WRITE(fp,*) ' '
      WRITE(fp,*) ' Neutral hydrogen molecule density: '

      DO ir = osm_watch1, osm_watch2
c       IPP/08 Krieger - ensure that idring(ir) is referenced inside
c       array bounds
c       IF (ir.LT.1.OR.ir.GT.nrs.OR.idring(ir).EQ.-1) CYCLE
        IF (ir.LT.1.OR.ir.GT.nrs) THEN
          CYCLE
        ELSEIF (idring(ir).EQ.-1) THEN
          CYCLE
        ENDIF

        WRITE(fp,'(2A4,6A12)') 'ik','ir','1','2','3','4','SUM','PINMOL'

        DO ik = 1, nks(ir)
          WRITE(fp,'(2I4,1P,6E12.4)')
     .      ik,ir,
     .     (pindata(ik,ir,i1),i1=H_MOL1,H_MOL3),
     .      pindata(ik,ir,H_MOL1)+pindata(ik,ir,H_MOL2)+
     .      pindata(ik,ir,H_MOL3)+pindata(ik,ir,H_MOL4),
     .      pinmol(ik,ir)
        ENDDO
      ENDDO

      CLOSE(fp)

      RETURN
99    STOP
      END

c
c...temporary substitution routine
c
      SUBROUTINE OutputGrid(fp,comment)
      IMPLICIT none

      INTEGER   fp
      CHARACTER comment*(*)

c     IPP/08 Krieger - changed output channel to 6 (lim file)
c     because of conflict with open statement in OutputData
c      write(fp,*) 'OutputGrid:',trim(comment)

      CALL OutputData(fp,comment)

      RETURN
      END
c
c     Write out just the grid data for comparisons - including polygons
c
      SUBROUTINE OutputGrid2(fp,comment)
      IMPLICIT none

      INTEGER   fp
      CHARACTER comment*(*)

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'grbound'
      INCLUDE 'slcom'

      INTEGER ik,ir,in,ii,i1,i2,id
      REAL minval(3,4),maxval(3,4),cs,GetCs

      CHARACTER*7   irtag(0:MAXNRS)
      CHARACTER*128 note

      write(6,*) 'OutputGrid2:',trim(comment)

      DO ir = 1, MAXNRS
        irtag(ir) = '      '
      ENDDO



c----------------------------------

      irtag(0) = 'PROBLEM'
      irtag(irsep)   = 'IRSEP  '
      irtag(irsep2)  = 'IRSEP2 '
      IF (nbr.GT.0) irtag(irbreak) = 'IRBREAK'
      irtag(irwall)  = 'IRWALL '
      irtag(irtrap)  = 'IRTRAP '
      irtag(nrs)     = 'NRS    '

      DO i1 = 1, 3
        DO i2 = 1, 4
          minval(i1,i2) =  HI
          maxval(i1,i2) = -HI
        ENDDO
      ENDDO

      OPEN(UNIT=fp,ACCESS='SEQUENTIAL',STATUS='REPLACE')

      CALL SetBounds

      WRITE(fp,'(2A)') 'DIVIMP data output: ',comment
      WRITE(fp,*)

      WRITE(fp,*)
      WRITE(fp,*) 'GRID DATA:'
      WRITE(fp,*)

      WRITE(fp,10) 'irsep    ',irsep   ,'irsep2   ',irsep2, 
     .             'irbreak  ',irbreak
      WRITE(fp,11) 'irwall   ',irwall  ,'irtrap   ',irtrap
      WRITE(fp,11) 'irwall2  ',irwall2 ,'irtrap2  ',irtrap2
      WRITE(fp,11) 'nbr      ',nbr     ,'nrs      ',nrs
      WRITE(fp,10) 'ikto     ',ikto    ,'ikti     ',ikti
      WRITE(fp,11) 'cutpt1   ',cutpt1  ,'cutpt2   ',cutpt2
      WRITE(fp,10) 'cutring  ',cutring ,'maxkpts  ',maxkpts
      WRITE(fp,12) 'maxrings ',maxrings
      WRITE(fp,10) 'npolyp   ',npolyp  ,'vpolyp   ',vpolyp
      WRITE(fp,12) 'vpolmin  ',vpolmin
      WRITE(fp,*)
      WRITE(fp,16) 'r0       ',r0      ,'z0       ',z0
      WRITE(fp,16) 'rxp      ',rxp     ,'zxp      ',zxp
      WRITE(fp,17) 'dthetg   ',dthetg

c f90
10    FORMAT(1X,3(A,I5,4X:))
11    FORMAT(1X,2(A,I5,4X))
12    FORMAT(1X,A,I5)
16    FORMAT(1X,2(A,F10.3,4X))
17    FORMAT(1X,A,F10.3)

      IF (outmode.EQ.3) WRITE(0,*) 'OUTPUTDATA: RING'

      WRITE(fp,*)

      WRITE(fp,*)
      WRITE(fp,*) 'RING DATA:'
      WRITE(fp,*)

      WRITE(fp,'(4A4,4A6,1X,A4,1X,2(A3,6X),3A10,A8,A5,1X,2A8,2A4)')
     .  'ir' ,'nks','id','org','idds2','idds1','ikto2','ikti2','ikm',
     .  'vl1','vl2','rhoIN (mm)','rho','rhoOUT','ksmaxs','ikb',
     .  'psin2','psin1','mp1','mp2'

      DO ir = 1, nrs
        WRITE(fp,'(4I4,4I6,1X,I4,1X,'//
     .           ' 2(I3,A,I2,A),3F10.5,1X,F7.3,I5,1X,2F8.4,
     .               2I4,1X,A)')
     .    ir,nks(ir),idring(ir),irorg2(ir),
     .    MAX(MIN(idds(ir,2),999),-999),
     .    MAX(MIN(idds(ir,1),999),-999),
     .    ikto2(ir),ikti2(ir),ikmids(ir),
     .    virloc(ir,IKLO),' (',virloc(ir,IKLO)-1,') ',
     .    virloc(ir,IKHI),' (',nks(ir)-virloc(ir,IKHI),') ',
     .    rho(ir,IN14)*1000.0,rho(ir,CELL1)*1000.0,
     .    rho(ir,OUT23)*1000.0,
     .    ksmaxs(ir),ikbreak(ir),psitarg(ir,2),psitarg(ir,1),
     .    ikmidplane(ir,IKLO),ikmidplane(ir,IKHI),irtag(ir)
      ENDDO

c
c      IF (outmode.EQ.3) WRITE(0,*) 'OUTPUTDATA: POLYGON'
c
      WRITE(fp,*)
      WRITE(fp,*) 'POLYGON DATA:'
c
      DO ir = 1, nrs
        WRITE(fp,*)
        WRITE(fp,'(2A3,A6,2A3,8A12,1X,A)')
     .    'ik','ir','in','vt','nv',
     .    'rvertp1','zvertp1','rvertp2','zvertp2',
     .    'rvertp3','zvertp3','rvertp4','zvertp4',irtag(ir)
c
        DO ik = 1, nks(ir)
          in = korpg(ik,ir)
c...temp: korpg=0
          IF (in.EQ.0) in = MAXNKS*MAXNRS
c
          note = ' '
          IF (ik.EQ.ikto2 (ir)) note = note(1:LEN_TRIM(note))//' IKTO2'
          IF (ik.EQ.ikti2 (ir)) note = note(1:LEN_TRIM(note))//' IKTI2'
          IF (ik.EQ.ikmids(ir)) note = note(1:LEN_TRIM(note))//' IKMIDS'
          IF (ik.EQ.ikbound(ir,IKLO))
     .      note = note(1:LEN_TRIM(note))//' IK1'
          IF (ik.EQ.ikbound(ir,IKHI))
     .      note = note(1:LEN_TRIM(note))//' IK2'
          if (in.eq. MAXNKS*MAXNRS)      
     >      note = trim(note)//' INVALID CELL'
       
c
          WRITE(fp,'(2I3,I6,2I3,8F12.8,A)')
     .      ik,ir,in,virtag(ik,ir),nvertp(in),
     .      rvertp(1,in),zvertp(1,in),rvertp(2,in),zvertp(2,in),
     .      rvertp(3,in),zvertp(3,in),rvertp(4,in),zvertp(4,in),
     .      note(1:LEN_TRIM(note))
        ENDDO
      ENDDO
c

      CLOSE(fp)
c-----------------------------------

      RETURN
      END

c
c ======================================================================
c
c subroutine: OuputData
c
c
      SUBROUTINE OutputData(fp,comment)
      IMPLICIT none

      INTEGER   fp
      CHARACTER comment*(*)

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'dynam1'
      INCLUDE 'pindata'
      INCLUDE 'grbound'
      INCLUDE 'slcom'

      INTEGER ik,ir,in,ii,i1,i2,id,iz,ike
      REAL minval(3,4),maxval(3,4),cs,GetCs

      CHARACTER*7   irtag(0:MAXNRS)
      CHARACTER*128 note
      CHARACTER     ring_tag(4)*4

      ring_tag(1) = 'SOL '
      ring_tag(2) = 'PFZ '
      ring_tag(3) = 'CORE'
      ring_tag(4) = 'N/A '

      write(6,*) 'OutputData:',fp,":",trim(comment)
c     IPP/09 Krieger - don't need this because we use unit 6
c     for write statement above
c     close(fp)

c f90 strange
      IF (stopopt.EQ.80) RETURN

      DO ir = 1, MAXNRS
        irtag(ir) = '      '
      ENDDO

      IF (outmode.EQ.3) WRITE(0,*) 'COMMENT:',comment

      irtag(0) = 'PROBLEM'
      irtag(irsep)   = 'IRSEP  '
      irtag(irsep2)  = 'IRSEP2 '
      IF (nbr.GT.0) irtag(irbreak) = 'IRBREAK'
      irtag(irwall)  = 'IRWALL '
      irtag(irtrap)  = 'IRTRAP '
      irtag(nrs)     = 'NRS    '

      DO i1 = 1, 3
        DO i2 = 1, 4
          minval(i1,i2) =  HI
          maxval(i1,i2) = -HI
        ENDDO
      ENDDO

      OPEN(UNIT=fp,ACCESS='SEQUENTIAL',STATUS='REPLACE')

      CALL SetBounds

      WRITE(fp,'(2A)') 'DIVIMP data output: ',comment
      WRITE(fp,*)

      IF (outmode.EQ.3) WRITE(0,*) 'OUTPUTDATA: HEADER'

      WRITE(fp,*)
      WRITE(fp,*) 'GRID DATA:'
      WRITE(fp,*)

      WRITE(fp,10) 'irsep    ',irsep   ,'irsep2   ',irsep2, 
     .             'irbreak  ',irbreak
      WRITE(fp,11) 'irwall   ',irwall  ,'irtrap   ',irtrap
      WRITE(fp,11) 'irwall2  ',irwall2 ,'irtrap2  ',irtrap2
      WRITE(fp,11) 'nbr      ',nbr     ,'nrs      ',nrs
      WRITE(fp,10) 'ikto     ',ikto    ,'ikti     ',ikti
      WRITE(fp,11) 'cutpt1   ',cutpt1  ,'cutpt2   ',cutpt2
      WRITE(fp,10) 'cutring  ',cutring ,'maxkpts  ',maxkpts
      WRITE(fp,12) 'maxrings ',maxrings
      WRITE(fp,10) 'npolyp   ',npolyp  ,'vpolyp   ',vpolyp
      WRITE(fp,12) 'vpolmin  ',vpolmin
      WRITE(fp,* ) 'vpolmin  ',vpolmin
      WRITE(fp,* ) 'nopriv   ',nopriv
      WRITE(fp,16) 'r0       ',r0      ,'z0       ',z0
      WRITE(fp,16) 'rxp      ',rxp     ,'zxp      ',zxp
      WRITE(fp,17) 'dthetg   ',dthetg

c f90
10    FORMAT(1X,3(A,I5,4X:))
11    FORMAT(1X,2(A,I5,4X))
12    FORMAT(1X,A,I5)
16    FORMAT(1X,2(A,F10.3,4X))
17    FORMAT(1X,A,F10.3)

      IF (outmode.EQ.3) WRITE(0,*) 'OUTPUTDATA: RING'

      WRITE(fp,*)

      WRITE(fp,*)
      WRITE(fp,*) 'RING DATA:'
      WRITE(fp,*)

      WRITE(fp,'(4A4,5A6,1X,A4,1X,2(A3,6X),3A10,A8,A5,1X,2A8,2A4)')
     .  'ir' ,'nks','id','org','type','idds2','idds1','ikto2','ikti2',
     .  'ikm','vl1','vl2','rhoIN (mm)','rho','rhoOUT','ksmaxs','ikb',
     .  'psin2','psin1','mp1','mp2'

      DO ir = 1, nrs
        i1 = ringtype(ir)
        IF (i1.LT.1.OR.i1.GT.3) i1 = 4
        WRITE(fp,'(4I4,2X,A4,4I6,1X,I4,1X,'//
     .           ' 2(I3,A,I2,A),3F10.5,1X,F7.3,I5,1X,2F8.4,
     .               2I4,1X,A)')
     .    ir,nks(ir),idring(ir),irorg2(ir),ring_tag(i1),
     .    MAX(MIN(idds(ir,2),999),-999),
     .    MAX(MIN(idds(ir,1),999),-999),
     .    ikto2(ir),ikti2(ir),ikmids(ir),
     .    virloc(ir,IKLO),' (',virloc(ir,IKLO)-1,') ',
     .    virloc(ir,IKHI),' (',nks(ir)-virloc(ir,IKHI),') ',
     .    rho(ir,IN14)*1000.0,rho(ir,CELL1)*1000.0,
     .    rho(ir,OUT23)*1000.0,
     .    ksmaxs(ir),ikbreak(ir),psitarg(ir,2),psitarg(ir,1),
     .    ikmidplane(ir,IKLO),ikmidplane(ir,IKHI),irtag(ir)
      ENDDO

      IF (outmode.EQ.3) WRITE(0,*) 'OUTPUTDATA: TARGET'

      WRITE(fp,*)
      WRITE(fp,*) 'TARGET DATA:'
      WRITE(fp,*)

      WRITE(fp,'(A10,I10)') ' ndsin    ',ndsin
      WRITE(fp,'(A10,I10)') ' ndsdiv   ',ndsdiv
      WRITE(fp,'(A10,I10)') ' nds      ',nds
      WRITE(fp,*)

      WRITE(fp,'(A4,2X,2A5,2X,2A8,2X,2A10,2X,4A12,2A4)')
     .  'in','ikds','irds','rp','zp','dds','dds2',
     .  'thetas','thetas2','costet','thetat','WI','NI'

c TEMP
c      DO in = nds, 1, -1
      DO in = 1, nds
        IF (nbr.GT.0.AND.in.EQ.ndsdiv+1) THEN
          note = 'WALL TARGETS'
        ELSE
          note = ' '
        ENDIF

        WRITE(fp,'(I4,2X,2I5,2X,2F8.4,2X,1P,2E10.2,0P,2X,4F12.6,2I4,
     .             1X,2A)')
     .    in,ikds(in),irds(in),rp(in),zp(in),dds(in),dds2(in),
     .    thetas(in)*180.0/PI,thetas2(in)*180.0/PI,costet(in),
     .    thetat(in),wallindex(in),nimindex(in),
     .    irtag(irds(in)),note(1:LEN_TRIM(note))

      ENDDO

      WRITE(fp,*)
      WRITE(fp,'(A4,2A8,5A12)')
     .  'in','kteds','ktids','knds','kvds','keds','kbfst1','kbfst2'

c TEMP
c      DO in = nds, 1, -1
      DO in = 1, nds
        IF (irds(in).EQ.0) CYCLE

        IF (nbr.GT.0.AND.in.EQ.ndsdiv+1) THEN
          note = 'WALL TARGETS'
        ELSE
          note = ' '
        ENDIF

        WRITE(fp,'(I4,2F8.3,1P,5E12.4,0P,1X,2A)')
     .    in,
     .    kteds(in),ktids(in),knds(in),kvds(in),keds(in),
     .    kbfst(irds(in),1),kbfst(irds(in),2),irtag(irds(in)),
     .    note(1:LEN_TRIM(note))
      ENDDO

      WRITE(fp,*)
      WRITE(fp,'(A4,4A14)')
     .  'ii','lpdati1 (ir)','lpdati2 (Te)','lpdati3 (Ti)','lpdati4 (n)'

      DO ii = 1, nlpdati
        WRITE(fp,'(I4,F14.1,2F14.5,1P,E14.6,0P,1X,A)')
     .    ii,lpdati(ii,1),lpdati(ii,2),lpdati(ii,3),lpdati(ii,4),
     .    irtag(INT(lpdati(ii,1)))
      ENDDO

      WRITE(fp,*)
      WRITE(fp,'(A4,4A14)')
     .  'ii','lpdato1 (ir)','lpdato2 (Te)','lpdato3 (Ti)','lpdato4 (n)'

      DO ii = 1, nlpdato
        WRITE(fp,'(I4,F14.1,2F14.5,1P,E14.6,0P,1X,A)')
     .    ii,lpdato(ii,1),lpdato(ii,2),lpdato(ii,3),lpdato(ii,4),
     .    irtag(INT(lpdato(ii,1)))
      ENDDO

      IF (outmode.EQ.3) WRITE(0,*) 'OUTPUTDATA: NEUTRAL WALL'

      WRITE(fp,*)
      WRITE(fp,*) 'NEUTRAL WALL DATA:'
      WRITE(fp,*)

      WRITE(fp,' (A,I10)  ') ' pcnt     ',pcnt
      WRITE(fp,' (A,I10)  ') ' ionwpts  ',ionwpts
      WRITE(fp,' (A,I10)  ') ' ioncpts  ',ioncpts
      WRITE(fp,'(2(A,I10))') ' nwall    ',nwall,' nwall2   ',nwall2
      WRITE(fp,'(2(A,I10))') ' nvesm    ',nvesm,' nvesp   ',nvesp
      WRITE(fp,'  (A,I10) ') ' nves     ',nves
      WRITE(fp,'  (A,I10) ') ' wallpts  ',wallpts
      WRITE(fp,' (A,2I10) ') ' wlwall1,2',wlwall1,wlwall2
      WRITE(fp,' (A,2I10) ') ' wltrap1,2',wltrap1,wltrap2

c      WRITE(fp,*)
c      WRITE(fp,'(A4,2A10,1X,2A10,1X,2A10,1X,4A10)')
c     .  'in','rw','zw',
c     .  'wallco1','wallco2','wallpt21','wallpt22',
c     .   'rvesm1','zvesm1','rvesm2','zvesm2'

      WRITE(fp,*)
      WRITE(fp,'(A4,1X,2A10)') 'in','rw','zw'
      DO in = 1, pcnt
        WRITE(fp,'(I4,1X,2F10.6)') in,rw(in),zw(in)
      ENDDO
      WRITE(fp,*)
      WRITE(fp,'(A4,1X,2A10)') 'in','riw','ziw'
      DO in = 1, ionwpts
        WRITE(fp,'(I4,1X,2F10.6)') in,riw(in),ziw(in)
      ENDDO
      WRITE(fp,*)
      WRITE(fp,'(A4,1X,2A10)') 'in','rcw','zcw'
      DO in = 1, ioncpts
        WRITE(fp,'(I4,1X,2F10.6)') in,riw(in),ziw(in)
      ENDDO

      WRITE(fp,*)
      WRITE(fp,'(A4,1X,2A10)') 'in','wallco1','wallco2'
      DO in = 1, nwall
        WRITE(fp,'(I4,1X,2F10.6)') in,wallco(in,1),wallco(in,2)
      ENDDO

      WRITE(fp,*)
      WRITE(fp,'(A4,1X,2A10)') 'in','wallco1_2','wallco2_2'
      DO in = 1, nwall2
        WRITE(fp,'(I4,1X,2F10.6)') in,wallco2(in,1),wallco2(in,2)
      ENDDO

      WRITE(fp,*)
      WRITE(fp,'(A4,1X,2A10)') 'in','wallpt21','wallpt22'
      DO in = 1, wallpts
        WRITE(fp,'(I4,1X,4F10.6)') in,wallpt2(in,1),wallpt2(in,2),
     .                             wallpt(in,16),wallpt(in,18)
      ENDDO

      WRITE(fp,*)
      WRITE(fp,'(A4,A6,12A10)') 'in','jvesm','rvesm1','zvesm1',
     .                          'rvesm2','zvesm2','fluxhw','flxhw2',
     .                          'flxhw3','flxhw4','flxhw5','flxhw6',
     .                          'flxhw7','flxhw8'
      DO in = 1, nvesm+nvesp
        WRITE(fp,'(I4,I6,4F10.6,1P,8E10.2,0P)')
     .    in,jvesm(in),rvesm(in,1),zvesm(in,1),rvesm(in,2),zvesm(in,2),
     .    fluxhw(in),flxhw2(in),flxhw3(in),flxhw4(in),flxhw5(in),
     .    flxhw6(in),flxhw7(in),flxhw8(in)
      ENDDO


      WRITE(fp,*)
      WRITE(fp,*) 'wallpt: '
      WRITE(fp,*)
      WRITE(fp,'(A4,9A10)') 'in','r','z',
     .  'weight_cc','weight_c','length_cc',
     .  'length_c' ,'length'  ,'angle_cc' ,'angle_c'
      WRITE(fp,'(A4,9A10)') '  ','1','2','3','4','5','6' ,'7','8' ,'9'

      DO in = 1, wallpts
        WRITE(fp,'(I4,9F10.5)') in,(wallpt(in,ii),ii=1,9)
      ENDDO

      WRITE(fp,*)
      WRITE(fp,'(A4,12A10)')
     .  'in'  ,'prob_cc','prob_c','prob','pprob_f',
     .  'type','nimbus' ,'target','temp','r1','z1','r2','z2'
      WRITE(fp,'(A4,12A10)')
     .  '  '  ,'10','11','12','13','16','17' ,'18','19','20','21','22',
     .         '23'

      DO in = 1, wallpts
        WRITE(fp,'(I4,12F10.5,1X,A)')
     .    in,(wallpt(in,ii),ii=10,13),(wallpt(in,ii),ii=16,23),    
     .    irtag(irds(MAX(1,NINT(wallpt(in,18)))))
      ENDDO

      WRITE(fp,*)
      WRITE(fp,'(A4,8A10)')
     .  'in'  ,'param','reflec.','IK','IR',
     .  'min dist.','Te' ,'Ti','ne'
      WRITE(fp,'(A4,12A10)')
     .  '  '  ,'24','25','26','27','28','29','30','31'

      DO in = 1, wallpts
        WRITE(fp,'(I4,8F10.5)')
     .    in,(wallpt(in,ii),ii=24,31)
      ENDDO

      IF (outmode.EQ.3) WRITE(0,*) 'OUTPUTDATA: PLASMA'

      WRITE(fp,*)
      WRITE(fp,*) 'PLASMA DATA:'

      DO ir = 1, nrs
        WRITE(fp,*) ' '
        WRITE(fp,'(2A3,8A11,3A10,1X,A7)')
     +    'ik'  ,'ir'  ,'kbfs' ,'bratio','kes',
     +    'kvhs','knbs','knes','natm','nmol',
     .    'Mach no.','ktibs','ktebs' ,irtag(ir)

        id = MAX(1,idds(ir,2))
        cs = GetCs(kteds(id),ktids(id)) + 1.0E-10

        IF (ir.GE.irsep)
     .    WRITE(fp,'(2I3,1P,E11.3,2(11X),2E11.3,3(11X),0P,
     .               3F10.4,A,I2,A)')
     .      idds(ir,2),ir,kbfst(ir,2),kvds(id),knds(id),
     .      kvds(id)/cs,ktids(id),
     .      kteds(id),' (',idds(ir,2),')'

        DO ik = 1, nks(ir)
          note = ' '
          IF (ik.EQ.ikto2 (ir)) note = note(1:LEN_TRIM(note))//' IKTO2'
          IF (ik.EQ.ikti2 (ir)) note = note(1:LEN_TRIM(note))//' IKTI2'
          IF (ik.EQ.ikmids(ir)) note = note(1:LEN_TRIM(note))//' IKMIDS'
          IF (ik.EQ.ikbound(ir,IKLO))
     .      note = note(1:LEN_TRIM(note))//' IK1'
          IF (ik.EQ.ikbound(ir,IKHI))
     .      note = note(1:LEN_TRIM(note))//' IK2'

          cs = GetCs(ktebs(ik,ir),ktibs(ik,ir)) + 1.0E-10

c          IF (kvhs(ik,ir).LT.1.0) cs = cs / qtim

          WRITE(fp,'(2I3,1P,3E11.4,5E11.3,0P,3F10.4,A)') ik,ir,
     .      kbfs (ik,ir),bratio(ik,ir),kes(ik,ir),
c     .      kvhs(ik,ir),knbs(ik,ir),knes(ik,ir),
     .      kvhs(ik,ir)/qtim,knbs(ik,ir),knes(ik,ir),
     .      pinatom(ik,ir),pinmol(ik,ir)*2,
c     .      kvhs(ik,ir)/cs,ktibs(ik,ir),ktebs (ik,ir),
     .      kvhs(ik,ir)/qtim/cs,ktibs(ik,ir),ktebs (ik,ir),
     .      note(1:LEN_TRIM(note))

        ENDDO

        id = MAX(1,idds(ir,1))
        cs = GetCs(kteds(id),ktids(id)) + 1.0E-10
        IF (ir.GE.irsep)
     .    WRITE(fp,'(2I3,1P,E11.3,2(11X),2E11.3,3(11X),0P,
     .               3F10.4,A,I2,A)') 
     .      idds(ir,1),ir,kbfst(ir,1),kvds(id),knds(id),
     .      kvds(id)/cs,
     .      ktids(id),kteds(id),' (',idds(ir,1),')'

      ENDDO

      IF (outmode.EQ.3) WRITE(0,*) 'OUTPUTDATA: CELL'

      WRITE(fp,*) ' '
      WRITE(fp,*) 'CELL DATA:'

      DO ir = 1, nrs
        WRITE(fp,*)
        WRITE(fp,'(2A4,5A4,4A10,3A11,1X,A7)')
     .    'ik','ir','vt','iki','iri','iko','iro',
     .    'kinds','koutds',
     .    'finds','foutds','kvols','kareas','thetag',irtag(ir)

        IF (ir.GE.irsep) THEN
          id = MAX(1,idds(ir,2))
          WRITE(fp,'(91X,F10.6,I4)') thetat(id),idds(ir,2)
        ENDIF

        DO ik = 1, nks(ir)
          note = ' '
          IF (ik.EQ.ikto2 (ir)) note = note(1:LEN_TRIM(note))//' IKTO2'
          IF (ik.EQ.ikti2 (ir)) note = note(1:LEN_TRIM(note))//' IKTI2'
          IF (ik.EQ.ikmids(ir)) note = note(1:LEN_TRIM(note))//' IKMIDS'
          IF (ik.EQ.ikbound(ir,IKLO))
     .      note = note(1:LEN_TRIM(note))//' IK1'
          IF (ik.EQ.ikbound(ir,IKHI))
     .      note = note(1:LEN_TRIM(note))//' IK2'

          WRITE(fp,'(2I4,5I4,4F10.6,1P,2E11.3,0P,F11.6,A)')
     .      ik,ir,virtag(ik,ir),
     .      ikins(ik,ir),irins (ik,ir),ikouts(ik,ir),irouts(ik,ir),
     .      kinds(ik,ir),koutds(ik,ir),finds (ik,ir),foutds(ik,ir),
     .      kvols(ik,ir),kareas(ik,ir),thetag(ik,ir),
c     .      kvols(ik,ir)*rxp/rs(ik,ir),kareas(ik,ir),thetag(ik,ir),
     .      note(1:LEN_TRIM(note))
        ENDDO

        IF (ir.GE.irsep) THEN
          id = MAX(1,idds(ir,1))
          WRITE(fp,'(91X,F10.6,I4)') thetat(id),idds(ir,1)
        ENDIF
      ENDDO

      IF (outmode.EQ.3) WRITE(0,*) 'OUTPUTDATA: GEOMETRY'

      WRITE(fp,*)
      WRITE(fp,*) 'GEOMETRY DATA:'

      DO ir = 1, nrs
        WRITE(fp,*)
        WRITE(fp,'(2A4,2A10,2A20,A20,A12,2A8,1X,A7)')
     .    'ik','ir',
     .    'rs','zs','kss (/ max)','ksb (% max)',
     .    'kps (/ max)','kpb','d_kpb','psi_n',
     .    irtag(ir)

        IF (ir.LT.irsep) THEN
          WRITE(fp,'(2I4,40X,F12.6,28X,F12.6)')
     .      0,ir,ksb(0,ir),kpb(0,ir)
        ELSE
          id = MAX(1,idds(ir,2))
          WRITE(fp,'(2I4,2F10.6,20X,F12.6,28X,F12.6,A,I2,A)')
     .      idds(ir,2),ir,rp(id),zp(id),ksb(0,ir),
     .      kpb(0,ir),' (',idds(ir,2),')'
        ENDIF

        DO ik = 1, nks(ir)
          note = ' '
          IF (ik.EQ.ikto2 (ir)) note = note(1:LEN_TRIM(note))//' IKTO2'
          IF (ik.EQ.ikti2 (ir)) note = note(1:LEN_TRIM(note))//' IKTI2'
          IF (ik.EQ.ikmids(ir)) note = note(1:LEN_TRIM(note))//' IKMIDS'
          IF (ik.EQ.ikbound(ir,IKLO))
     .      note = note(1:LEN_TRIM(note))//' IK1'
          IF (ik.EQ.ikbound(ir,IKHI))
     .      note = note(1:LEN_TRIM(note))//' IK2'

c          WRITE(fp,'(2I3,2F10.6,2(F12.6,F8.4),F12.6,F8.4,F12.6,A,F9.2)')
          WRITE(fp,'(2I4,2F10.6,2(F12.6,F8.4),F12.6,F8.4,F12.6,
     .               2F8.4,A)')
     .      ik,ir,
     .      rs (ik,ir),zs (ik,ir),
     .      kss(ik,ir),kss(ik,ir)/(ksmaxs(ir)+1.0E-10),
     .      ksb(ik,ir),(ksb(ik,ir)-ksb(ik-1,ir))/(ksmaxs(ir)+1.0E-10)*
     .                  100.0,
     .      kps(ik,ir),kps(ik,ir)/(kpmaxs(ir)+1.0E-10),
     .      kpb(ik,ir),kpb(ik,ir)-kpb(ik-1,ir),
     .      psifl(ik,ir),
     .      note(1:LEN_TRIM(note))
c     .      kss2(ik,ir)

        ENDDO

        IF (ir.GE.irsep) THEN
          id = MAX(1,idds(ir,1))
          WRITE(fp,'(2I4,2F10.6,2X,F10.6,60X,A,I2,A)')
     .      idds(ir,1),ir,rp(id),zp(id),ksmaxs(ir),
     .      ' (',idds(ir,1),')'
        ENDIF
      ENDDO

      IF (outmode.EQ.3) WRITE(0,*) 'OUTPUTDATA: POLYGON'

      WRITE(fp,*)
      WRITE(fp,*) 'POLYGON DATA:'

      DO ir = 1, nrs
        WRITE(fp,*)
        WRITE(fp,'(2A3,A6,2A3,8A12,1X,A)')
     .    'ik','ir','in','vt','nv',
     .    'rvertp1','zvertp1','rvertp2','zvertp2',
     .    'rvertp3','zvertp3','rvertp4','zvertp4',irtag(ir)

        DO ik = 1, nks(ir)
          in = korpg(ik,ir)

c...temp: korpg=0
          IF (in.EQ.0) in = MAXNKS*MAXNRS

          note = ' '
          IF (ik.EQ.ikto2 (ir)) note = note(1:LEN_TRIM(note))//' IKTO2'
          IF (ik.EQ.ikti2 (ir)) note = note(1:LEN_TRIM(note))//' IKTI2'
          IF (ik.EQ.ikmids(ir)) note = note(1:LEN_TRIM(note))//' IKMIDS'
          IF (ik.EQ.ikbound(ir,IKLO))
     .      note = note(1:LEN_TRIM(note))//' IK1'
          IF (ik.EQ.ikbound(ir,IKHI))
     .      note = note(1:LEN_TRIM(note))//' IK2'

          WRITE(fp,'(2I3,I6,2I3,8F12.7,A)')
     .      ik,ir,in,virtag(ik,ir),nvertp(in),
     .      rvertp(1,in),zvertp(1,in),rvertp(2,in),zvertp(2,in),
     .      rvertp(3,in),zvertp(3,in),rvertp(4,in),zvertp(4,in),
     .      note(1:LEN_TRIM(note))
        ENDDO
      ENDDO

      IF (ctestsol.GE.0.0) THEN
        IF (outmode.EQ.3) WRITE(0,*) 'OUTPUTDATA: IMPURITIES'
        WRITE(fp,*)
        WRITE(fp,*) 'DIVIMP IMPURITY DATA:'
        DO ir = 1, nrs
          WRITE(fp,*)
          WRITE(fp,'(2A4,2A10,9A10,1X,A7)')
     .      'ik','ir','r','z',
     .      '1','2','3','4','5','10','15','20','25',
     .      irtag(ir)
          ike = nks(ir) 
          IF (ir.LT.irsep) ike = nks(ir) - 1
          DO ik = 1, ike
            note = ' '
            IF (ik.EQ.ikto2 (ir)) note = TRIM(note)//' IKTO2'
            IF (ik.EQ.ikti2 (ir)) note = TRIM(note)//' IKTI2'
            IF (ik.EQ.ikmids(ir)) note = TRIM(note)//' IKMIDS'
            IF (ik.EQ.ikbound(ir,IKLO))
     .        note = note(1:LEN_TRIM(note))//' IK1'
            IF (ik.EQ.ikbound(ir,IKHI))
     .        note = note(1:LEN_TRIM(note))//' IK2'
            WRITE(fp,'(2I4,2F10.6,1P,9E10.2,0P)')
     .        ik,ir,
     .        rs (ik,ir),zs (ik,ir),
     .        (SNGL(ddlims(ik,ir,iz)),iz=1,5),
     .        (SNGL(ddlims(ik,ir,iz)),iz=10,25,5)
          ENDDO
        ENDDO
      ENDIF

      WRITE(fp,*)
      WRITE(fp,*) 'DALPHA DATA:'
c
      DO ir = 1, nrs
        WRITE(fp,*)
        WRITE(fp,'(2A3,7A12,1X,A)')
     .    'ik','ir',
     .    'pinalpha','pinline6','pinline1','pinline2',
     .    'pinline3','pinline4','pinline5',irtag(ir)
c
        DO ik = 1, nks(ir)
          note = ' '
          IF (ik.EQ.ikto2 (ir)) note = note(1:LEN_TRIM(note))//' IKTO2'
          IF (ik.EQ.ikti2 (ir)) note = note(1:LEN_TRIM(note))//' IKTI2'
          IF (ik.EQ.ikmids(ir)) note = note(1:LEN_TRIM(note))//' IKMIDS'
          IF (ik.EQ.ikbound(ir,IKLO))
     .      note = note(1:LEN_TRIM(note))//' IK1'
          IF (ik.EQ.ikbound(ir,IKHI))
     .      note = note(1:LEN_TRIM(note))//' IK2'
          if (in.eq. MAXNKS*MAXNRS)      
     >      note = trim(note)//' INVALID CELL'

          WRITE(fp,'(2I3,1P,7E12.4,0P,A)')
     .      ik,ir,
     .      pinalpha(ik,ir),
     .      pinline (ik,ir,6  ,H_BALPHA),
     .      pinline (ik,ir,1:5,H_BALPHA),
c     .      pinline (ik,ir,7  ,H_BALPHA),
     .      note(1:LEN_TRIM(note))
        ENDDO
      ENDDO

      CLOSE(fp)

      IF (outmode.EQ.3) WRITE(0,*) 'OUTPUTDATA: DONE'

      RETURN
      END

c ======================================================================
c
c subroutine: OutputEIRENE
c
c
c Volume quantities confirmed:
c
c
c Boundary quantities confirmed:
c
c
c
      SUBROUTINE OutputEIRENE(fp,note)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER      fp
      CHARACTER*(*) note
      INTEGER      ik,ir,in
      CHARACTER*20 tag

      CHARACTER*64 machine2

      CALL GetEnv('DIVNAME',machine2)

      IF (machine2(1:LEN_TRIM(machine2)).EQ.
     .    'juelich'.AND.nrs.GT.10) THEN
        WRITE(fp,*) 'CASE RUN IN GERMANY, NO EIRENE OUTPUT'
        RETURN
      ENDIF


      IF (osm_watch1.EQ.0.OR.osm_watch2.EQ.0) RETURN

      CALL MS('OutputEIRENE','Turned off -- sometimes')
      CALL DB('Writing PIN data')
      CALL SetBounds


      IF (sloutput.AND.citersol.NE.0) THEN
        WRITE(pinout,*) '***NOT WRITING PIN DATA FILE***'
        RETURN
      ENDIF

     
c      IF (cmodopt.EQ.1.AND.citersol.GT.0.AND.nitersol.GT.5.AND.
c     .    fp.NE.PINOUT.AND.
c     .    rel_count.NE.0.AND.rel_step.NE.rel_nstep) RETURN


      IF (sloutput) THEN
        WRITE(pinout,*) 'WRITING PIN DATA OUTPUT FILE'
        WRITE(pinout,*) '  '//note(1:LEN_TRIM(note))
      ENDIF

c      WRITE(0,*) 'SKIPPING TO AVOID BUS ERROR!?'
c      RETURN
c
c     Unit 67 is the PIN output before each SOL 22 iteration, so the file
c     needs to remain open:
c
      IF (fp.NE.67.AND.fp.NE.SLOUT.AND.fp.NE.PINOUT)
     .  OPEN(UNIT=fp,ACCESS='SEQUENTIAL',STATUS='REPLACE')

      WRITE(fp,*) 'EIRENE data '//note(1:LEN_TRIM(note))//':'

      DO ir = osm_watch1, osm_watch2

c       IPP/08 Krieger - ensure that idring(ir) is referenced inside
c       array bounds
c       IF (ir.LT.1.OR.ir.GT.nrs.OR.idring(ir).EQ.BOUNDARY) CYCLE
        IF (ir.LT.1.OR.ir.GT.nrs) THEN
          CYCLE
        ELSEIF (idring(ir).EQ.BOUNDARY) THEN
          CYCLE
        ENDIF

        WRITE(fp,*) ' '

        WRITE(fp,'(2A4,7A12,4A6)')
     +    'ik','ir','pinatom','pinion','osmion','pinrec','osmrec',
     .    'osmcfp','pinalpha','mulrec','mulion','mulqer','mulqei'

        DO ik = 1, nks(ir)
          tag = '                    '
          IF (ik.EQ.ikbound(ir,IKLO)) tag = tag(1:LEN_TRIM(tag))//'IK1'
          IF (ik.EQ.ikbound(ir,IKHI)) tag = tag(1:LEN_TRIM(tag))//'IK2'

          IF (ir.EQ.2.AND.ik.EQ.1) THEN
            WRITE(fp,*) 'RXP,RS',rxp,rs(ik,ir),pinion(ik,ir)
            WRITE(fp,*) 'VOLUME',kvols(ik,ir)*rxp/rs(ik,ir)
          ENDIF

          WRITE(fp,'(2I4,1P,7E12.4,0P,4F6.3,1X,A)')
     .      ik,ir,
     .      pinatom(ik,ir),pinion(ik,ir),
c     .      pinatom(ik,ir),pinion(ik,ir)*rs(ik,ir)/rxp*1.6E-19*1.0E-06,
     .      osmion  (ik,ir),pinrec(ik,ir),
     .      osmrec (ik,ir),osmcfp(ik,ir),pinalpha(ik,ir),
     .      mulrec (ik,ir),mulion(ik,ir),mulqer  (ik,ir),mulqei(ik,ir),
     .      tag(1:LEN_TRIM(tag))

        ENDDO

        WRITE(fp,*) ' '
        WRITE(fp,'(2A4,9A12)')
     +    'ik','ir',
     +    'pinionz','pinenz','pinmol','pinmoi','pinqi','pinqe','osmqe',
     .    'osmpei' ,'osmpmk'

        DO ik = 1, nks(ir)
          tag = '                    '
          IF (ik.EQ.ikbound(ir,IKLO)) tag = tag(1:LEN_TRIM(tag))//'IK1'
          IF (ik.EQ.ikbound(ir,IKHI)) tag = tag(1:LEN_TRIM(tag))//'IK2'

          WRITE(fp,'(2I4,1P,9E12.4,0P,2F8.2,1X,A)')
     +      ik,ir,
     +      pinionz(ik,ir),pinenz (ik,ir),pinmol (ik,ir),pinmoi(ik,ir),
     +      pinqi  (ik,ir),pinqe  (ik,ir),osmqe  (ik,ir),
     .      osmpei (ik,ir),osmpmk (ik,ir),
     .      osm_dp4(ik,ir),osm_dp6(ik,ir),
     .      tag(1:LEN_TRIM(tag))
        ENDDO

        WRITE(fp,*) ' '
        WRITE(fp,'(2A4,7A12)')
     +    'ik','ir',
     +    'pinena','osmcfe','osmcfi','pinmp','osmmp','osmmp2','pinqe2'


        DO ik = 1, nks(ir)
          tag = '                    '
          IF (ik.EQ.ikbound(ir,IKLO)) tag = tag(1:LEN_TRIM(tag))//'IK1'
          IF (ik.EQ.ikbound(ir,IKHI)) tag = tag(1:LEN_TRIM(tag))//'IK2'

          WRITE(fp,'(2I4,1P,7E12.4,0P,1X,A)')
     +      ik,ir,
     .      pinena(ik,ir),osmcfe(ik,ir),osmcfi(ik,ir),pinmp (ik,ir),
     .      osmmp (ik,ir),osmmp2(ik,ir),pinqe2(ik,ir),
     .      tag(1:LEN_TRIM(tag))
        ENDDO

      ENDDO

      IF (fp.NE.67.AND.fp.NE.SLOUT.AND.fp.NE.PINOUT) CLOSE(fp)

c      CALL DB('done')

      IF (sloutput) WRITE(pinout,*) 'DONE'

      RETURN
      END
c
c ======================================================================
c
c
c ======================================================================
c
c subroutine: WIO
c
      SUBROUTINE WIO(tag1,id1,optval,tag2,id2,maxlen)

      IMPLICIT none

      INTEGER   optval,id1,id2,maxlen
      CHARACTER tag1*(*),tag2*(*)

      INCLUDE 'params'
      INCLUDE 'slcom'

      CHARACTER tag3*100,spaces*40
      INTEGER   fpin,fpout,taglen

      fpin   = 98
      fpout  = 7
      spaces = '                                        '

      OPEN(UNIT=fpin,FILE='info.dat',ACCESS='SEQUENTIAL',
     .     STATUS='OLD',ERR=98)

      taglen = MIN(20,LEN_TRIM(tag1))

      WRITE(fpout,'(2A,I3,A)') spaces(1:id1)//tag1(1:taglen)//
     .                         ' OPTION: ',
     .                         spaces(1:20-taglen),optval,
     .                         spaces(1:20)//tag2

      WRITE(tag3,'(A,I3)') tag2,optval

      CALL GetInfo(fpin,fpout,tag3,id2,maxlen)

      CLOSE(fpin)

      RETURN
c
c     jdemod - changed to be non-critical error if the documentation
c              is not found
c

98    CALL WN('WIO','Unstructured Input Documentation'//
     >   ' Data File ''info.dat'' not found')
      WRITE(fpout,*) '   ERROR: INFO.DAT MESSAGE TEXT FILE NOT FOUND'
99    return
c
c     jdemod
c
      END
c
c ======================================================================
c
c subroutine: OutputOptionInformation
c
      SUBROUTINE OutputOptionInformation
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER ii,in,i1,i2,ir,fp


      fp = 7

      CALL PRC(' ')
      CALL PRC('SLMOD')
      CALL PRC(' ')
c
c     General options:
c
      CALL WIO('PROBE DATA',3,tarsource,'TARSOURCE',5,67)
      IF (tarsource.GT.0) THEN
        CALL PRI('     SHOT       = '  ,shot)
        CALL PRI('     SLICE TIME =   ',slice_time)

        CALL PRC(' ')
        WRITE(fp,'(5X,3A10)') 'Rho (m)','Te (eV)','Ne (m-3)'
        DO i1 = 1, NUMPRB
          DO i2 = 1, prb_num(i1)
            WRITE(fp,'(I5,F10.4,F10.2,1P,E10.2)')
     .        i1,prb_rho(i2,i1),prb_te(i2,i1),dip_ne(i2,i1)
          ENDDO
        ENDDO

        CALL PRC(' ')
        WRITE(fp,'(5X,A5,4A10)') 'Ring','S (m)','Te (eV)','Ti (eV)',
     .                           'Ne (m-3)'
        DO ir = irsep, irwall-1
          WRITE(fp,'(5X,I5,3F10.2,1P,E10.2)')
     .      ir,dip_s (ir,FSP1),dip_te(ir,FSP1),
     .         dip_ti(ir,FSP1),dip_ne(ir,FSP1)
        ENDDO
      ENDIF


      CALL WIO('PROBE ALIGNMENT',3,prb_align,'PRB_ALIGN',5,67)
      IF (prb_align.EQ.1) THEN
        CALL PRR('     FSP SHIFT = ',prb_shift)
      ENDIF
c
c     OSM output:
c
      IF (cioptf.EQ.22) THEN
        CALL PRC(' ')
        CALL PRR('   SOL 22 HARD ERROR REGION (*SMAX) = ',osm_range)
        CALL PRC(' ')
      ENDIF
      CALL PRC(' ')

c...  EIRENE options:
      CALL SLOPT02(fp)
c
c     Grid options:
c
      IF (nbr.GT.0) THEN
        CALL PRC('  BROKEN GRID DATA:')
        CALL PRI('     NUMBER OF BROKEN RINGS = ',nbr)
        CALL PRI('     FIRST BROKEN RING      = ',irbreak)
        CALL PRR('     INNER WALL SHIFT (mm)  = ',grd_shift*1000.0)
        CALL PRC(' ')
      ENDIF

      CALL WIO('TARGET REFINEMENT',3,grd_refine,'GRD_REFINE',5,67)
      IF (grd_refine.GT.0) THEN
        CALL PRR('     REFINEMENT DISTANCE = '  ,grd_range)
      ENDIF

      CALL WIO('GRID ADAPTATION',3,adp_opt,'ADP_OPT',5,67)
      IF (adp_opt.GT.0) THEN
        CALL PRI('     ADAPTATION REGION = '  ,adp_region)
        CALL PRR('     REFINEMENT LIMIT  = '  ,adp_upper)
        CALL PRR('     RECOVERY LIMIT    = '  ,adp_lower)
      ENDIF
c
c     Relaxation options:
c
      CALL WIO('RELAXATION',3,rel_opt,'REL_OPT',5,67)
      IF (rel_opt.EQ.1.OR.rel_opt.EQ.3.OR.rel_ndata.GT.0) THEN
        CALL PRR('     PIN SOURCE REL. FRACTION = '  ,rel_frac)
      ENDIF
      IF (rel_opt.EQ.2.OR.rel_opt.EQ.3.OR.rel_ndata.GT.0) THEN
        CALL PRI('     STEPS FOR BC RELAXATION   = '  ,rel_nstep)
        CALL PRI('     ITERATIONS PER STEP       = '  ,rel_niter)
        CALL PRR('     PACE FOR BC RELAXATION    = '  ,rel_pace)
      ENDIF
      IF (rel_opt.GT.0) THEN
        CALL PRI('     OUTPUT PIN SOURCES       = '  ,out_source)
        CALL PRI('     OUTPUT PLASMA DATA       = '  ,out_plasma)
        CALL PRI('     OUTPUT GEOMETRY DATA     = '  ,out_geom)
      ENDIF

c...  Setup option:
      CALL SLOPT01(fp)

      IF (nlpdati2.GT.0) THEN
        CALL PRC('  INNER (CMOD OUTER) BC RELAXATION DATA:')
        WRITE(fp,'(5X,A5,2X,2(3A10,2X))')
     .    'Ring','Te1 (eV)','Ti1 (eV)','Ne1/Isat1',
     .           'Te2 (eV)','Ti2 (eV)','Ne2/Isat2'
        DO in = 1, nlpdati2
          WRITE(fp,'(5X,F5.1,2X,2(2F10.2,1P,E10.2,0P,2X))')
     .      (lpdati2(in,ii),ii=1,7)
        ENDDO
      ENDIF
      IF (nlpdato2.GT.0) THEN
        CALL PRC('  OUTER (CMOD INNER) BC RELAXATION DATA:')
        WRITE(fp,'(5X,A5,2X,2(3A10,2X))')
     .    'Ring','Te1 (eV)','Ti1 (eV)','Ne1/Isat1',
     .           'Te2 (eV)','Ti2 (eV)','Ne2/Isat2'
        DO in = 1, nlpdato2
          WRITE(fp,'(8X,F5.1,2X,2(2F10.2,1P,E10.2,0P,2X))')
     .      (lpdato2(in,ii),ii=1,7)
        ENDDO
      ENDIF

      IF (rel_ndata.GT.0) THEN
        CALL PRC('  RELAXATION OVER-RIDE DATA:')
        WRITE(fp,'(5X,6A11)')
     .    '       Step','    # Iter.','  Rel. Opt.',' Rel. Frac.',
     .    'Adapt. Opt.',' EIRENE (s)'
        DO in = 1, rel_ndata
          WRITE(fp,'(5X,7F11.3)') (rel_data(in,ii),ii=1,7)
        ENDDO
      ENDIF
c
c     Output options:
c
      CALL PRC(' ')

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: SLOPT01
c
      SUBROUTINE SLOPT01(fp)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'printopt'
      INCLUDE 'slcom'

      INTEGER fp

c...  Target data shift:
      CALL WIO('TARGET DATA SHIFT OPTION',3,tarshiftopt,'TARSHIFTOPT',
     .         5,67)      

      CALL HD(fp,'  TARGET PROBE DATA SHIFT','TARSHIFT-HD',5,67)
      WRITE(fp,*)
      WRITE(coment,80) 'INNER JET; OUTER C-MOD,DIII-D SHIFT (m) = ',
     .                 tarshift(IKHI)
      CALL PRC(coment)
      WRITE(coment,80) 'OUTER JET; INNER C-MOD,DIII-D SHIFT (m) = ',
     .                 tarshift(IKLO)
      CALL PRC(coment)
80    FORMAT(4X,A,F10.3)

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: OutputPressureTable
c
c
      SUBROUTINE OutputPressureTable(fp)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'printopt'
      INCLUDE 'slcom'

      INTEGER      fp,i1,in
      REAL         p_exp,n_exp,t_exp
      CHARACTER*30 gaugename(110)


      do in = 1,110
         gaugename(in) = ' '
      end do


c      gaugename(001) = 'Midplane gauge'
c      gaugename(002) = 'Divertor gas box gauge'

      gaugename(101) = 'PBF1 pump chamber'
      gaugename(102) = 'PBF2'
      gaugename(103) = 'PBF3'
      gaugename(104) = 'PV1  private zone'
      gaugename(105) = 'PR2'
      gaugename(106) = 'VPLOWS'
      gaugename(107) = 'PCM105BAF'
      gaugename(108) = 'PCM240TOR'


      IF (eirnpgdat.GT.0) THEN
        WRITE(fp,*)
        WRITE(fp,30) 'CODE','X (m)','Y (m)','P_exp','P_D','P_D2'
        DO i1 = 1, eirnpgdat
c...      Multiply by 7.502 to convert from Pa to mTorr:
          WRITE(fp,31) INT(eirpgdat(i1,1)),
     .                 eirpgdat(i1, 2),eirpgdat(i1, 3),eirpgdat(i1,7),
     .                 pinasd(i1+1,2,1,1)*ECH*1.0E+06*0.67*7.502,
     .                 pinasd(i1+1,2,3,1)*ECH*1.0E+06*0.67*7.502,
     .                 gaugename(INT(eirpgdat(i1,1)))

c...      Multiply by 7.502 to convert from Pa to mTorr:
c          WRITE(fp,31) INT(eirpgdat(i1,1)),
c     .                 eirpgdat(i1, 2),eirpgdat(i1, 3),eirpgdat(i1,7),
c     .                 eirpgdat(i1,13)*7.502,
c     .                 eirpgdat(i1,10)*7.502,
c     .                 gaugename(INT(eirpgdat(i1,1)))
        ENDDO

        WRITE(fp,*)
        WRITE(fp,30) 'CODE','X (m)','Y (m)','T_exp','T_D','T_D2'
        DO i1 = 1, eirnpgdat
          t_exp = 300.0 * 1.38E-23 / ECH

          WRITE(fp,32) INT(eirpgdat(i1,1)),
     .          eirpgdat(i1, 2),eirpgdat(i1, 3),t_exp,
     .          pinasd(i1+1,2,1,1)/(pinasd(i1+1,1,1,1)+1.0E-10)*0.67,
     .          pinasd(i1+1,2,3,1)/(pinasd(i1+1,1,3,1)+1.0E-10)*0.67,' '
        ENDDO

        WRITE(fp,*)
        WRITE(fp,30) 'CODE','X (m)','Y (m)','n_exp','n_D','n_D2'
        DO i1 = 1, eirnpgdat
          p_exp = eirpgdat(i1,7) / 7.502
          n_exp = p_exp / (1.38E-23 * 300.0)

          WRITE(fp,31) INT(eirpgdat(i1,1)),
     .                 eirpgdat(i1, 2),eirpgdat(i1, 3),n_exp,
     .                 pinasd(i1+1,1,1,1)*1.0E+06,
     .                 pinasd(i1+1,1,3,1)*1.0E+06,' '
        ENDDO

        WRITE(fp,*)
        WRITE(fp,*) '    AVERAGES (over iteration):'

        WRITE(fp,*)
        WRITE(fp,30) 'CODE','X (m)','Y (m)','P_exp','P_D','P_D2'
        DO i1 = 1, eirnpgdat
c...      Multiply by 7.502 to convert from Pa to mTorr:
          WRITE(fp,31) INT(eirpgdat(i1,1)),
     .                 eirpgdat(i1, 2),eirpgdat(i1, 3),eirpgdat(i1,7),
     .                 pinasd(i1+1,2,1,2)*ECH*1.0E+06*0.67*7.502,
     .                 pinasd(i1+1,2,3,2)*ECH*1.0E+06*0.67*7.502,
     .                 gaugename(INT(eirpgdat(i1,1)))
        ENDDO

        WRITE(fp,*)
        WRITE(fp,30) 'CODE','X (m)','Y (m)','T_exp','T_D','T_D2'
        DO i1 = 1, eirnpgdat
          t_exp = 300.0 * 1.38E-23 / ECH

          WRITE(fp,32) INT(eirpgdat(i1,1)),
     .          eirpgdat(i1, 2),eirpgdat(i1, 3),t_exp,
     .          pinasd(i1+1,2,1,2)/(pinasd(i1+1,1,1,2)+1.0E-10)*0.67,
     .          pinasd(i1+1,2,3,2)/(pinasd(i1+1,1,3,2)+1.0E-10)*0.67,' '
        ENDDO

        WRITE(fp,*)
        WRITE(fp,30) 'CODE','X (m)','Y (m)','n_exp','n_D','n_D2'
        DO i1 = 1, eirnpgdat
          p_exp = eirpgdat(i1,7) / 7.502
          n_exp = p_exp / (1.38E-23 * 300.0)

          WRITE(fp,31) INT(eirpgdat(i1,1)),
     .                 eirpgdat(i1, 2),eirpgdat(i1, 3),n_exp,
     .                 pinasd(i1+1,1,1,2)*1.0E+06,
     .                 pinasd(i1+1,1,3,2)*1.0E+06,' '
        ENDDO

      ELSE
        WRITE(fp,*)
        WRITE(fp,*) '    No pressure gauges defined.'
      ENDIF
30    FORMAT(5X,A4,2A7     ,3A10)
31    FORMAT(5X,I4,2F7.3,1P,3E10.2,0P,2X,A)
32    FORMAT(5X,I4,2F7.3,   3F10.3   ,2X,A)



      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: OutputVacuumCellTable
c
c
      SUBROUTINE OutputVacuumCellTable(fp)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'printopt'
      INCLUDE 'slcom'

      INTEGER fp,i1,i2,i3,i4,i5,in,init,cell,cut1,add
      LOGICAL writeline
      REAL    teD,teD2,neD,neD2,pD,pD2,fact,x,y,z,normcont(MAXSTRATA),
     .        sum

      DATA init /0/
      SAVE
 
      IF (eirbgk.EQ.2.OR.eirbgk.EQ.3.OR.eirbgk.EQ.4) THEN
c...    EIRENE additional cell data (vacuum grid):

        IF (init.EQ.0) THEN
c...      Write header, but only for the first call to this routine:
          CALL HD(fp,'  NEUTRAL HYDROGEN VACUUM CELL DATA',
     .               'EIRADDCELL-HD',5,67)
          init = 1
        ENDIF

        IF (.NOT.thesis.AND.rel_opt.GE.1.AND.rel_iter.EQ.rel_niter) THEN
c...      The conditions on this IF statement should be the same as for
c         the block that writes the iter-vac.dat file in module
c         eirene.d6a:
          WRITE(fp,*)
          WRITE(fp,'(5X,A,I4)') 'STEP:',rel_step

          WRITE(fp,*)          
          WRITE(fp,5) 'TeMi',te_mult_o,'TiMi',ti_mult_o,'nMi',n_mult_o,
     .                'EIRTIME',eirtime
          WRITE(fp,5) 'TeMo',te_mult_i,'TiMo',ti_mult_i,'nMo',n_mult_i
5         FORMAT(2X,3(3X,A,'=',F7.2):3X,A,'=',I5)
        ENDIF

        WRITE(fp,*)
        WRITE(fp,30) 'Ind','x','y','z',' T_D','T_D2','  n_D',' n_D2',
     .               '    P_D','   P_D2'
        WRITE(fp,30) '   ','(m)','(m)','(m)','(eV)','(eV)',
     .               '(m-3)','(m-3)','(mTorr)','(mTorr)'
     .                               
        DO i1 = 1+(eirnpgdat+1), asc_ncell*ascncut+(eirnpgdat+1)

          writeline = .FALSE.
          in        = i1 - (1 + eirnpgdat)
          DO i4 = 1, eirnaout
            add = MAX(0,NINT(eiraout(i4,2)-1.0) * asc_ncell)
            IF ((eiraout(i4,1).EQ.1.AND.
     .           (eiraout(i4,3)+add.EQ.in.OR.
     .            eiraout(i4,4)+add.EQ.in.OR.
     .            eiraout(i4,5)+add.EQ.in.OR.
     .            eiraout(i4,6)+add.EQ.in)).OR.
     .          (eiraout(i4,1).EQ.3.AND.
     .           (eiraout(i4,3)+add.LE.in.AND.
     .            eiraout(i4,4)+add.GE.in.OR.
     .            eiraout(i4,5)+add.LE.in.AND.
     .            eiraout(i4,6)+add.GE.in)))
     .       writeline = .TRUE.
          ENDDO

          IF (.NOT.writeline) CYCLE

c...      Find coordinates of the center of the cell:
          x     = 0.0
          y     = 0.0
          z     = 0.0 
          cell = MOD(in,asc_ncell)            
          IF (cell.EQ.0) cell = asc_ncell
          DO i2 = 1, ascnvertex(cell)
            x = x + ascvertex(2*i2-1,cell) / REAL(ascnvertex(cell))
            y = y + ascvertex(2*i2  ,cell) / REAL(ascnvertex(cell))
          ENDDO          
          IF     (asc_3dmode.EQ.1) THEN
            z = 0.5 * (asc_zmin3D(cell) + asc_zmax3D(cell))
          ELSEIF (asc_3dmode.EQ.2) THEN
            cut1 = INT(REAL(in) / (REAL(asc_ncell) + 0.5)) + 1
            z = 0.5 * (asc_zmin3D(cut1) + asc_zmax3D(cut1))
          ENDIF

          teD  = pinasd(i1,2,1,2) / (pinasd(i1,1,1,2)+1.0E-10) * 0.667
          teD2 = pinasd(i1,2,3,2) / (pinasd(i1,1,3,2)+1.0E-10) * 0.667
          neD  = pinasd(i1,1,1,2) * 1.0E+06
          neD2 = pinasd(i1,1,3,2) * 1.0E+06

          fact = ECH * 1.0E+06 * 0.67 * 7.502
          pD   = pinasd(i1,2,1,2) * fact
          pD2  = pinasd(i1,2,3,2) * fact

          WRITE(fp,31) i1-(eirnpgdat+1),x,y,z,teD,teD2,neD,neD2,pD,pD2
        ENDDO


        IF (eirnstrdat.NE.0) THEN
c...      Data per stratum:
          WRITE(fp,*)
          WRITE(fp,*) '    DENSITY CONTRIBUTION PER STRATUM:'
          WRITE(fp,*)

c        DO i1 = 1, eirnstrdat
c          i3 = 1
c          i4 = eirnstrai + 1
c          WRITE(fp,'(I6,1P,10(E12.4:))')
c     .      eirstrlis(i1),(eirstrdat(i1,i2,2*i3-1+2*eirnatmi+1),
c     .      i2=1,eirnstrai),eirstrdat(i1,i4,2*i3-1+2*eirnatmi+1)
c        ENDDO


          i3 = 1
          i4 = eirnstrai + 1

          DO i1 = 1, eirnstrdat

            CALL RZero(normcont,MAXSTRATA)
            sum = 0.0
            DO i5 = 1, eirnstrai
              normcont(i5) = eirstrdat(i1,i5,2*i3-1+2*eirnatmi+1) /
     .                       eirfluxt(i5)
              sum = sum + normcont(i5)
            ENDDO
            DO i5 = 1, eirnstrai
              normcont(i5) = normcont(i5) / sum
            ENDDO
            sum = 0.0
            DO i5 = 1, eirnstrai
              sum = sum + eirfluxt(i5)
            ENDDO

            IF (eirstrdat(i1,i4,2*i3-1+2*eirnatmi+1).EQ.0.0) THEN
              WRITE(fp,'(4X,I6,2X,A)')
     .          eirstrlis(i1),' nil'
            ELSE
              WRITE(fp,'(4X,I6,10(F6.1,2(A,F4.2),A:))')
     .          eirstrlis(i1)-1-eirnpgdat,
     .          (eirstrdat(i1,i2,2*i3-1+2*eirnatmi+1)/
     .           eirstrdat(i1,i4,2*i3-1+2*eirnatmi+1)*100.0,
     .          '%(',normcont(i2),',',eirfluxt(i2)/sum,')',
     .          i2=1,eirnstrai)
            ENDIF
          ENDDO

        ENDIF


30      FORMAT(5X,A5,1X,3A7  ,1X,2A8     ,1X,2A10     ,1X,2A8  )
31      FORMAT(5X,I5,1X,3F7.3,1X,2F8.3,1P,1X,2E10.2,0P,1X,2F8.3)
      ENDIF

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: SLOPT02
c
c Output to the .dat file that is related to EIRENE.
c
c
c
      SUBROUTINE SLOPT02(fp)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'printopt'
      INCLUDE 'slcom'

      INTEGER fp

      INTEGER i1,i2,fp1,ios

      CHARACTER*30  surmat(20)
      CHARACTER*256 buffer

      surmat(1) = '    MOLYBDENUM'
      surmat(2) = '        CARBON'
      surmat(3) = '      TUNGSTEN'
      surmat(4) = '     BERYLLIUM'

c...  EIRENE options:
      CALL PRB
      CALL WIO('EIRENE GEOMETRY DATA',3,eirgeom,'EIRGEOM',5,67)
      WRITE(fp,*)

      CALL WIO('EIRENE INPUT FILE',3,eirdata,'EIRDATA',5,67)
      WRITE(fp,*)
      IF (eirdata.EQ.1) THEN
        CALL PRI('    RUN TIME (CPU seconds)           ',eirtime)
        CALL PRC('    SURFACE MATERIAL: TARGET  '//surmat(eirmat1))
        CALL PRC('                      WALL    '//surmat(eirmat2))
        CALL PRR('    SURFACE TEMPERATURE (eV): TARGET ',ABS(eirtemp1))
        CALL PRR('                              WALL   ',ABS(eirtemp2))
        CALL PRI('    GRID TYPE OPTION                 ',eirgrid)
c        CALL PRI('    WALL DATA PRECISION OPTION       ',eirneut)
c        CALL PRI('    DEBUG OPTION                     ',eirdebug)
        WRITE(fp,*)
      ENDIF

      CALL WIO('N-N COLLISION'       ,3,eirbgk    ,'EIRBGK'    ,5,67)
      WRITE(fp,*)

      CALL WIO('LYMAN ALPHA OPACITY' ,3,eiropacity,'EIROPACITY',5,67)
      WRITE(fp,*)

      CALL WIO('CX D2+ PRODUCTION'   ,3,eircxd2   ,'EIRCXD2'   ,5,67)
      WRITE(fp,*)

      CALL WIO('PROTON-D2 COLLISIONS',3,eirph2    ,'EIRPH2'    ,5,67)

c...  Time-to-ionsation data from EIRENE:
      CALL OutputIonisationTimeData(fp)

c...  Pressure gauge data:
      CALL HD(fp,'  PRESSURE GAUGE DATA','EIRPGAUGE-HD',5,67)
      CALL OutputPressureTable(fp)

c...  Iteration data for '04:
      fp1 = 97
      OPEN(UNIT=fp1,FILE='eirene04.txt',STATUS='OLD',IOSTAT=ios)
      IF (ios.EQ.0) THEN
c...    Iteration data file is there, so insert into .dat file:
        DO WHILE (.TRUE.)
          CALL ReadLine(fp1,buffer,1,*9 ,*97)
          WRITE(fp,'(A)') buffer(1:LEN_TRIM(buffer))
        ENDDO
 9      CLOSE (fp1)
      ENDIF

c...  Data for vacuum grid additional cells:
      fp1 = 97
      OPEN(UNIT=fp1,FILE='iter-vac.dat',STATUS='OLD',IOSTAT=ios)
      IF (ios.EQ.0) THEN
c...    Iteration data file is there, so insert into .dat file:
        DO WHILE (.TRUE.)
          CALL ReadLine(fp1,buffer,1,*10,*97)
          WRITE(fp,'(A)') buffer(1:LEN_TRIM(buffer))
        ENDDO
10      CLOSE (fp1)
      ELSE
c...    Otherwise, just output the data to the .dat file:
        CALL OutputVacuumCellTable(fp)
      ENDIF

c...  Fluxes to transparent and pumping surfaces.  Just transfer the
c     contents of tmp-001.dat to the .dat file:
      CALL HD(fp,'  TRANSPARENT AND ABSORBING WALL SURFACE FLUXES',
     .           'EIRNONSTSUR-HD',5,67)
      fp1 = 98
      OPEN(UNIT=fp1,FILE='tmp-001.dat',STATUS='OLD',IOSTAT=ios)
      IF (ios.EQ.0) THEN
c...    Data for each step is there, so insert into .dat file:
        DO WHILE (.TRUE.)
          CALL ReadLine(fp1,buffer,1,*20,*97)
          WRITE(fp,'(A)') buffer(1:LEN_TRIM(buffer))
        ENDDO
20      CLOSE (fp1)
      ELSE
c...    Otherwise:
        WRITE(fp,*)
        WRITE(fp,'(5X,A)') 'No data available.'
      ENDIF

c...  Fluxes to surfaces:
      CALL HD(fp,'  NEUTRAL HYDROGEN FLUX DATA','EIRNEUTFLUX-HD',5,67)
      WRITE(fp,*)

      WRITE(fp,40) 'Ind',' r ',' z ',' D_parflx',' D_avgeng',
     .                               'D2_parflx','D2_avgeng',
     .                               'Im_parflx','Im_avgeng',
     .                               'D+_parflx','D+_avgeng'
      WRITE(fp,40) '   ','(m)','(m)','(m-2 s-1)','   (eV)  ',   
     .                               '(m-2 s-1)','   (eV)  ',
     .                               '(m-2 s-1)','   (eV)  ',
     .                               '(m-2 s-1)','   (eV)  '
c...  Loop over NWALL for now, to avoid BGK grid:
      DO i1 = 1, nvesm+nvesp
c
c jdemod - adjusted to print table for entire wall - not just below Xpoint
c
c        IF (0.5*(zvesm(i1,1)+zvesm(i1,2)).LT.zxp) THEN
          WRITE(fp,41) i1,       
     .      0.5*(rvesm(i1,1)+rvesm(i1,2)),0.5*(zvesm(i1,1)+zvesm(i1,2)),
     .      flxhw6(i1),flxhw5(i1),fluxhw(i1)-flxhw6(i1),flxhw7(i1),
     .      flxhw3(i1),flxhw4(i1),flxhw8(i1)           ,-1.0
c
c        ENDIF
c

      ENDDO
40    FORMAT(5X,A3,2A7     ,8A11)
41    FORMAT(5X,I3,2F7.3,4(1P,E11.2,0P,F11.3))      



      RETURN
97    CALL ER('SLOPT02','Problem reading iteration data file',*99)
98    CALL ER('SLOPT02','Cannot find iteration data file'    ,*99)
99    STOP
      END
c
c
c
c
c ======================================================================
c
c
c
c
      SUBROUTINE OutputRegionIntegrals(source)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      REAL source(MAXNKS,MAXNRS)

      INTEGER SymmetryPoint

      INTEGER i1,ik,ir,ikm
      REAL    rdum(20)


      CALL RZero(rdum,20)

      DO ir = irsep, irwall-1
        ikm = SymmetryPoint(ir)
        DO ik = 1, ikm
          rdum(1) = rdum(1) + source(ik,ir) * kvols(ik,ir)
        ENDDO
        DO ik = ikm+1, nks(ir)
          rdum(2) = rdum(2) + source(ik,ir) * kvols(ik,ir)
        ENDDO
      ENDDO
      DO ir = irtrap+1, nrs
        DO ik = 1, nks(ir)
          rdum(3) = rdum(3) + source(ik,ir) * kvols(ik,ir)
        ENDDO
      ENDDO
      DO ir = 2, irsep-1
        DO ik = 1, nks(ir)-1
          rdum(4) = rdum(4) + source(ik,ir) * kvols(ik,ir)
        ENDDO
      ENDDO

      WRITE(79,'(5X,1P,8E12.4,0P)') (rdum(i1),i1=1,4)

      RETURN
99    STOP
      END


c
c ======================================================================
c
      SUBROUTINE ShowStats
      IMPLICIT   none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER ir,i1
      REAL    rms(4),peak(4)

      COMMON /OLDPLASMA/ oldknbs ,oldktebs ,oldktibs ,oldkvhs ,
     .                   oldknbs2,oldktebs2,oldktibs2,oldkvhs2
      REAL oldktebs (MAXNKS,MAXNRS),oldktibs (MAXNKS,MAXNRS),
     .     oldknbs  (MAXNKS,MAXNRS),oldkvhs  (MAXNKS,MAXNRS),
     .     oldktebs2(MAXNKS,MAXNRS),oldktibs2(MAXNKS,MAXNRS),
     .     oldknbs2 (MAXNKS,MAXNRS),oldkvhs2 (MAXNKS,MAXNRS)

      IF (rel_frac.EQ.0.0) RETURN

      DO i1 = 1, 4
        rms (i1) = 0.0
        peak(i1) = 0.0
      ENDDO

      WRITE(PINOUT,*)
      WRITE(PINOUT,'(A)') 'RMS/PEAK solution variation'

      WRITE(PINOUT,'(A4,3(A20,2X))') 'ir','te','ti','n'

      DO ir = irsep, irwall-1
        IF (idring(ir).EQ.-1) CYCLE

        CALL GetRelStat(ktebs,oldktebs2,ir,rms(1),peak(1))
        CALL GetRelStat(ktibs,oldktibs2,ir,rms(2),peak(2))
        CALL GetRelStat(knbs ,oldknbs2 ,ir,rms(3),peak(3))

        IF (nrs.EQ.6) THEN
          IF (rms(1)/rel_frac.LT.0.01) THEN
            WRITE(0     ,*) 'Steady Te profile'
            WRITE(PINOUT,*) 'Steady Te profile'
          ENDIF
          IF (rms(2)/rel_frac.LT.0.01) THEN
            WRITE(0     ,*) 'Steady Ti profile'
            WRITE(PINOUT,*) 'Steady Ti profile'
          ENDIF
          IF (rms(3)/rel_frac.LT.0.01) THEN
            WRITE(0     ,*) 'Steady ne profile'
            WRITE(PINOUT,*) 'Steady ne profile'
          ENDIF
        ENDIF

        WRITE(PINOUT,'(I4,3(2F10.4,2X),0P)')
     .    ir,(rms(i1)/rel_frac,peak(i1)/rel_frac,i1=1,3)

        WRITE(79,'(A,1X,4I4,1X,F10.4,2X,3(2F12.6,2X))')
     .    '''PLASMAVAR 1.01''',
     .    rel_step,rel_iter,rel_count,ir,
     .    rel_frac,(rms(i1)/rel_frac,peak(i1)/rel_frac,i1=1,3)

      ENDDO

c      WRITE(PINOUT,'(A4,4(A20,2X))') 'ir','ion','rec','qe','qi'
c
c      DO ir = irsep, irwall-1
c        IF (idring(ir).EQ.-1) CYCLE
c
c        CALL GetRelStat(pinion,pinion2,ir,rms(1),peak(1))
c        CALL GetRelStat(pinrec,pinrec2,ir,rms(2),peak(2))
c        CALL GetRelStat(pinqe ,pinqe2 ,ir,rms(3),peak(3))
c        CALL GetRelStat(pinqi ,pinqi2 ,ir,rms(4),peak(4))
c
c        WRITE(PINOUT,'(I4,1P,4(2F10.4,2X),0P)')
c     .    ir,(rms(i1),peak(i1),i1=1,4)
c
c      ENDDO




      RETURN
99    STOP
      END

c
c ======================================================================
c
      SUBROUTINE GetRelStat(array1,array2,ir,rms,peak)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      REAL    array1(MAXNKS,MAXNRS),array2(MAXNKS,MAXNRS),
     .        rms,peak,cnt,diff,squ
      INTEGER ik,ir

      squ  =  0.0
      rms  = -1.0
      peak = -1.0
      cnt  =  0.0

      DO ik = 1, nks(ir)
        IF (array1(ik,ir).NE.0.0.AND.array2(ik,ir).NE.0.0) THEN
          diff = ABS(array1(ik,ir) - array2(ik,ir)) / array2(ik,ir)
          peak = MAX(diff,peak)
          squ  = squ + diff**2.0
          cnt  = cnt + 1.0

c          WRITE(PINOUT,*) array1(ik,ir),array2(ik,ir)
        ENDIF
      ENDDO

      IF (cnt.GT.0.0) rms = SQRT(squ / cnt)

c      WRITE(PINOUT,*) rms,squ,cnt

      RETURN
99    STOP
      END


c
c ======================================================================
c
      SUBROUTINE MirrorOldPlasma(te,ti,ne,vb)
      IMPLICIT none

      INCLUDE 'params'

      REAL te(MAXNKS,MAXNRS),ti(MAXNKS,MAXNRS),ne(MAXNKS,MAXNRS),
     .     vb(MAXNKS,MAXNRS) 

      COMMON /OLDPLASMA/ oldknbs ,oldktebs ,oldktibs ,oldkvhs ,
     .                   oldknbs2,oldktebs2,oldktibs2,oldkvhs2
      REAL oldktebs (MAXNKS,MAXNRS),oldktibs (MAXNKS,MAXNRS),
     .     oldknbs  (MAXNKS,MAXNRS),oldkvhs  (MAXNKS,MAXNRS),
     .     oldktebs2(MAXNKS,MAXNRS),oldktibs2(MAXNKS,MAXNRS),
     .     oldknbs2 (MAXNKS,MAXNRS),oldkvhs2 (MAXNKS,MAXNRS)

      INTEGER ik,ir

      DO ir = 1, MAXNRS
        DO ik = 1, MAXNKS
          oldktebs2(ik,ir) = te(ik,ir)
          oldktibs2(ik,ir) = ti(ik,ir)
          oldknbs2 (ik,ir) = ne(ik,ir)
          oldkvhs2 (ik,ir) = vb(ik,ir)
        ENDDO
      ENDDO

      RETURN
99    STOP
      END
c
c ======================================================================
c
c
c
      SUBROUTINE AnalyseContinuity(fp)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INCLUDE 'solparams'
      INCLUDE 'solcommon'
      INCLUDE 'solswitch'

      INTEGER fp

      INTEGER GetModel,SymmetryPoint
      REAL    CalcPressure,GetJsat

      INTEGER          ik,ik1,ik2,ik3,ik4,ik5,ir,i1,id1,id2,
     .                 model1(MAXNRS),model2(MAXNRS)

      DOUBLE PRECISION intg(20)
      CHARACTER*7      irtag(0:MAXNRS)


      irtag(irsep)   = 'IRSEP  '
      irtag(irbreak) = 'IRBREAK'
      irtag(irwall)  = 'IRWALL '
      irtag(irtrap)  = 'IRTRAP '
      irtag(nrs)     = 'NRS    '

      WRITE(fp,*)
      WRITE(fp,*) 'Detailed continuity analysis:'


      DO i1 = 1, 20
        intg(i1) = 0.0
      ENDDO


      DO ir = irsep, nrs
        IF (idring(ir).EQ.-1) CYCLE

        model1(ir) = GetModel(IKLO,ir)
        model2(ir) = GetModel(IKHI,ir)
      ENDDO

      CALL SetBounds

c
c
c     Conservation:
c
c
c     Tally sources and fluxes:
c
      WRITE(fp,*)
      WRITE(fp,*) 'Fluxes:'
      WRITE(fp,'(3X,5A4,A5,6A12)')
     .  'ik1','ik2','ik3','ik4','ik5','ir',
     .  'flux1','flux2','flux3','flux4','flux5','total'

      DO ir = irsep, nrs
        IF (idring(ir).EQ.-1) CYCLE

        IF ((model1(ir).GE.22.AND.model1(ir).LE.24).AND.
     .      (model2(ir).GE.22.AND.model2(ir).LE.24)) THEN
          ik1 = 1
          ik2 = ikbound(ir,IKLO)
          ik3 = SymmetryPoint(ir)
          ik4 = ikbound(ir,IKHI)
          ik5 = nks(ir)

          intg(1) = -ABS(knds(idds(ir,2)) * kvds(idds(ir,2)))
          intg(2) =  gtarg(ir,2)
          intg(3) =  knbs(ik3,ir) * kvhs(ik3,ir)
          intg(4) = -gtarg(ir,1)
          intg(5) =  ABS(knds(idds(ir,1)) * kvds(idds(ir,1)))

          intg(6) = intg(1) - intg(5)

          WRITE(fp,'(3X,5I4,I5,1P,6E12.3,0P)')
     .      ik1,ik2,ik3,ik4,ik5,ir,(intg(i1),i1=1,6)

          WRITE(79,'(A,1X,3I4,1X,I4,5E12.4)') '''FLUXES    1.00''',
     .      rel_step,rel_iter,rel_count,ir,(intg(i1),i1=1,5)
        ENDIF
      ENDDO
c
c
c
c
c
      IF (switch(SWION  ).EQ.1.0.OR.
     .    switch(SWION  ).EQ.2.0) THEN

        WRITE(fp,*)
        WRITE(fp,*) 'Ionisation:'
        WRITE(fp,'(3X,5A4,A5,5A12)')
     .    'ik1','ik2','ik3','ik4','ik5','ir',
     .    'ion1','ion2','ion3','ion4','total'

        DO ir = irsep, nrs
          IF (idring(ir).EQ.-1) CYCLE

          IF ((model1(ir).GE.22.AND.model1(ir).LE.24).AND.
     .        (model2(ir).GE.22.AND.model2(ir).LE.24)) THEN
            ik1 = 1
            ik2 = ikbound(ir,IKLO)
            ik3 = SymmetryPoint(ir)
            ik4 = ikbound(ir,IKHI)
            ik5 = nks(ir)

            CALL CalcIntegral2(pinion(1,ir),ik1,ik2,ir,intg(1),2)
            CALL CalcIntegral2(pinion(1,ir),ik2,ik3,ir,intg(2),2)
            CALL CalcIntegral2(pinion(1,ir),ik3,ik4,ir,intg(3),2)
            CALL CalcIntegral2(pinion(1,ir),ik4,ik5,ir,intg(4),2)

            intg(5) = intg(1) + intg(2) + intg(3) + intg(4)

            WRITE(fp,'(3X,5I4,I5,1P,5E12.3,0P)')
     .        ik1,ik2,ik3,ik4,ik5,ir,(intg(i1),i1=1,5)

            WRITE(79,'(A,1X,3I4,1X,I4,5E12.4)') '''ION       1.00''',
     .        rel_step,rel_iter,rel_count,ir,(intg(i1),i1=1,5)
          ENDIF
        ENDDO
      ENDIF
c
c
c
c
c
      IF (switch(SWRECOM).EQ.1.0) THEN

        WRITE(fp,*)
        WRITE(fp,*) 'Recombination:'
        WRITE(fp,'(3X,5A4,A5,5A12)')
     .    'ik1','ik2','ik3','ik4','ik5','ir',
     .    'rec1','rec2','rec3','rec4','total'

        DO ir = irsep, nrs
          IF (idring(ir).EQ.-1) CYCLE

          IF ((model1(ir).GE.22.AND.model1(ir).LE.24).AND.
     .        (model2(ir).GE.22.AND.model2(ir).LE.24)) THEN
            ik1 = 1
            ik2 = ikbound(ir,IKLO)
            ik3 = SymmetryPoint(ir)
            ik4 = ikbound(ir,IKHI)
            ik5 = nks(ir)

            CALL CalcIntegral2(pinrec(1,ir),ik1,ik2,ir,intg(1),2)
            CALL CalcIntegral2(pinrec(1,ir),ik2,ik3,ir,intg(2),2)
            CALL CalcIntegral2(pinrec(1,ir),ik3,ik4,ir,intg(3),2)
            CALL CalcIntegral2(pinrec(1,ir),ik4,ik5,ir,intg(4),2)

            intg(5) = intg(1) + intg(2) + intg(3) + intg(4)

            WRITE(fp,'(3X,5I4,I5,1P,5E12.3,0P)')
     .        ik1,ik2,ik3,ik4,ik5,ir,(-intg(i1),i1=1,5)

            WRITE(79,'(A,1X,3I4,1X,I4,5E12.4)') '''REC       1.00''',
     .        rel_step,rel_iter,rel_count,ir,
     .        intg(1),intg(2),intg(3),intg(4),intg(5)
          ENDIF
        ENDDO
      ENDIF
c
c
c
c
c
      IF (switch(SWGPERP).EQ.2.0.AND.piniter) THEN

        WRITE(fp,*)
        WRITE(fp,*) 'Cross-field source:'
        WRITE(fp,'(3X,5A4,A5,5A12)')
     .    'ik1','ik2','ik3','ik4','ik5','ir',
     .    'cfp1','cfp2','cfp3','cfp4','total'

        DO ir = irsep, nrs
          IF (idring(ir).EQ.-1) CYCLE

c          IF ((model1(ir).EQ.22.OR.model1(ir).EQ.24).AND.
c     .        (model2(ir).EQ.22.OR.model2(ir).EQ.24)) THEN
            ik1 = 1
            ik2 = ikbound(ir,IKLO)
            ik3 = SymmetryPoint(ir)
            ik4 = ikbound(ir,IKHI)
            ik5 = nks(ir)

            CALL CalcIntegral2(osmcfp(1,ir),ik1,ik2,ir,intg(1),2)
            CALL CalcIntegral2(osmcfp(1,ir),ik2,ik3,ir,intg(2),2)
            CALL CalcIntegral2(osmcfp(1,ir),ik3,ik4,ir,intg(3),2)
            CALL CalcIntegral2(osmcfp(1,ir),ik4,ik5,ir,intg(4),2)

            intg(5) = intg(1) + intg(2) + intg(3) + intg(4)

            WRITE(fp,'(3X,5I4,I5,1P,5E12.3,0P)')
     .        ik1,ik2,ik3,ik4,ik5,ir,(intg(i1),i1=1,5)

            WRITE(79,'(A,1X,3I4,1X,I4,5E12.4)') '''CFP       1.00''',
     .        rel_step,rel_iter,rel_count,ir,(intg(i1),i1=1,5)
c          ENDIF
        ENDDO
      ENDIF



      IF (.TRUE.) THEN

        WRITE(fp,*)
        WRITE(fp,*) 'Hydrogen atom density:'
        WRITE(fp,'(3X,A5,5A12)') 'ir','atm1','atm2','atm3','atm4',
     .                           'total'

        DO ir = irsep, nrs
          IF (idring(ir).EQ.-1) CYCLE

          IF ((model1(ir).GE.22.AND.model1(ir).LE.24).AND.
     .        (model2(ir).GE.22.AND.model2(ir).LE.24)) THEN
            ik1 = 1
            ik2 = ikbound(ir,IKLO)
            ik3 = SymmetryPoint(ir)
            ik4 = ikbound(ir,IKHI)
            ik5 = nks(ir)

            CALL CalcIntegral2(pinatom(1,ir),ik1,ik2,ir,intg(1),2)
            CALL CalcIntegral2(pinatom(1,ir),ik2,ik3,ir,intg(2),2)
            CALL CalcIntegral2(pinatom(1,ir),ik3,ik4,ir,intg(3),2)
            CALL CalcIntegral2(pinatom(1,ir),ik4,ik5,ir,intg(4),2)

            intg(5) = intg(1) + intg(2) + intg(3) + intg(4)

            WRITE(fp,'(3X,I5,1P,5E12.3,0P)') ir,(intg(i1),i1=1,5)

            WRITE(79,'(A,1X,3I4,1X,I4,5E12.4)') '''NATOM     1.00''',
     .        rel_step,rel_iter,rel_count,ir,(intg(i1),i1=1,5)
          ENDIF
        ENDDO
      ENDIF

      IF (.TRUE.) THEN

        WRITE(fp,*)
        WRITE(fp,*) 'Hydrogen molecule (D2,D2+) density:'
        WRITE(fp,'(3X,A5,5A12)') 'IR','REGION1','REGION2','REGION3',
     .                           'REGION4','TOTAL'

        DO ir = irsep, nrs
          IF (idring(ir).EQ.-1) CYCLE

          IF ((model1(ir).GE.22.AND.model1(ir).LE.24).AND.
     .        (model2(ir).GE.22.AND.model2(ir).LE.24)) THEN
            ik1 = 1
            ik2 = ikbound(ir,IKLO)
            ik3 = SymmetryPoint(ir)
            ik4 = ikbound(ir,IKHI)
            ik5 = nks(ir)

            CALL CalcIntegral2(pinmol(1,ir),ik1,ik2,ir,intg(1),2)
            CALL CalcIntegral2(pinmol(1,ir),ik2,ik3,ir,intg(2),2)
            CALL CalcIntegral2(pinmol(1,ir),ik3,ik4,ir,intg(3),2)
            CALL CalcIntegral2(pinmol(1,ir),ik4,ik5,ir,intg(4),2)

            intg(5) = intg(1) + intg(2) + intg(3) + intg(4)

            WRITE(fp,'(3X,I5,1P,5E12.3,0P)') ir,(intg(i1),i1=1,5)

            WRITE(79,'(A,1X,3I4,1X,I4,5E12.4)') '''NMOL      1.00''',
     .        rel_step,rel_iter,rel_count,ir,(intg(i1),i1=1,5)
          ENDIF
        ENDDO
      ENDIF







      RETURN
99    STOP
      END













