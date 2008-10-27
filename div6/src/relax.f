c     -*-Fortran-*-
c
c ======================================================================
c
c subroutine: LoadPIN
c
      SUBROUTINE LoadPIN
      IMPLICIT   none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER SymmetryPoint,GetModel
      REAL    CalcPressure

      INTEGER   ihold(10),i1,i2,u1,u2,u3,hold_nrs,hold_irwall,ik,ik1,ir,
     .          hold_irbreak,hold_irsep,error1,error2,error3,id,fp
      LOGICAL   loadgeo
      REAL      rhold(10),rdum1,p,p1,pl,pl1,
     .          eirpgdat2(MAXNAS,MAXASD),eirasdat2(MAXNAS,MAXASD),
     .          eirspdat2(MAXNAS,MAXASD)
      CHARACTER fname*128

      COMMON /OLDPLASMA/ oldknbs ,oldktebs ,oldktibs ,oldkvhs ,
     .                   oldknbs2,oldktebs2,oldktibs2,oldkvhs2
      REAL oldktebs (MAXNKS,MAXNRS),oldktibs (MAXNKS,MAXNRS),
     .     oldknbs  (MAXNKS,MAXNRS),oldkvhs  (MAXNKS,MAXNRS),
     .     oldktebs2(MAXNKS,MAXNRS),oldktibs2(MAXNKS,MAXNRS),
     .     oldknbs2 (MAXNKS,MAXNRS),oldkvhs2 (MAXNKS,MAXNRS)

      COMMON /PININIT/ pininit
      LOGICAL          pininit(MAXNRS)

c
c     jdemod - moved variables to fix alignment issue 
c
      COMMON /PEIMULCOM/ osm_peimul2,osmpeires,
     .                   restrictosmpmk,
     .                   restrictmachno
      LOGICAL            restrictosmpmk,restrictmachno
      REAL*8             osm_peimul2(2,MAXNRS),osmpeires(2,MAXNRS)



      COMMON /OSMPMKCOM/ save_osmpmk,scaleplateau
      REAL*8             save_osmpmk(0:MAXNKS,MAXNRS)
      LOGICAL            scaleplateau(MAXNRS)

      CALL DB('Entering LoadPIN')

      fp = PINOUT

c      CALL OutputData(85,'Before loading')

c...Here?
      IF (osm_probe.GE.1) CALL InterpolateProbeData(osm_probe)

      ihold(1) = nitersol
      ihold(2) = rel_opt
      rhold(3) = rel_frac
      ihold(4) = rel_nstep
      ihold(5) = rel_niter
      ihold(6) = rel_step
      ihold(7) = rel_iter
      ihold(8) = rel_count
      ihold(9) = adp_opt

      hold_irsep   = irsep
      hold_irbreak = irbreak
      hold_irwall  = irwall
      hold_nrs     = nrs

      DO i1 = 1, MAXNAS
        DO i2 = 1, MAXASD
          eirpgdat2(i1,i2) = eirpgdat(i1,i2)
        ENDDO
      ENDDO

      error1 = 0
      error2 = 0
      error3 = 0

      u1 = 40
      u2 = 41
      u3 = 42

c      u1 = 150
c      u2 = 151
c      u3 = 152

      OPEN(u1,FILE='source.dat',FORM='UNFORMATTED',STATUS='OLD',ERR=97)
      OPEN(u2,FILE='plasma.dat',FORM='UNFORMATTED',STATUS='OLD',ERR=96)
      OPEN(u3,FILE='geomty.dat',FORM='UNFORMATTED',STATUS='OLD',ERR=98)

      CALL GetEnv('SAVENAME',fname)

      WRITE(fp,*)
      WRITE(fp,'(A   )') 'Loading background plasma:'
      WRITE(fp,'(A   )') '  FILE = '//fname(1:LEN_TRIM(fname))
      WRITE(fp,'(A,I3)') '  STEP = ',osm_store

c      IF (outmode.GE.2) WRITE(0,'(A,::)') 'Reading iteration '


c      CALL SetBounds
c      WRITE(0,*) 'MARK: IKBOUND12= ',ikbound(3,IKLO),ikbound(3,IKHI)


c      WRITE(0     ,*) 'NOT LOADING GEOMETRY DATA IN LOADPIN'
c      WRITE(PINOUT,*) 'NOT LOADING GEOMETRY DATA IN LOADPIN'

      IF (outmode.GE.2) WRITE(0,*) 'OSM_STORE=',osm_store,rel_step

      loadgeo = .TRUE.

      CALL OutputData(85,'Before loading solution')

      DO i1 = 0, osm_store
        IF (outmode.GE.2) WRITE(0,'(1X,2I2)') i1,osm_store

        CALL ReadSources(u1,error1)
        CALL ReadPlasma (u2,error2)

        IF (loadgeo.OR.adp_opt.GT.0) THEN
          loadgeo = .FALSE.
          WRITE(0,*) 'CALLING READGEOMETRY'
          CALL ReadGeometry(u3,error3)
        ENDIF

        IF (error1.NE.0.OR.error2.NE.0.OR.error3.NE.0)
     .    CALL ER('LoadPIN','Data unavailable',*99)
      ENDDO

      WRITE(0,'(A,I3,A)') ' Read step ',i1-1,' of '//
     .                    fname(1:LEN_TRIM(fname))
c      IF (outmode.GE.2) WRITE(0,'(A,I3,A)') ' Read step ',i1-1,' of '//
c     .                                      fname(1:LEN_TRIM(fname))


      CLOSE(u1)
      CLOSE(u2)
      CLOSE(u3)

      CALL OutputData(86,'After loading solution')

      tagpinatom = .TRUE.

c      WRITE(0,*) 'MARK: IKBOUND12= ',ikbound(3,IKLO),ikbound(3,IKHI)

c
c
c
c      WRITE(0,*)
c      WRITE(0,*)
c      WRITE(0,*)
c      WRITE(0,*) '***** DELETING OSMCFPFLX! *****'
c      WRITE(0,*)
c      WRITE(0,*)
c      WRITE(0,*)
c      osmcfpflx = 0.0  


c...check that target corner points match up...

      IF (hold_irsep .NE.irsep .OR.hold_irbreak.NE.irbreak.OR.
     .    hold_irwall.NE.irwall.OR.hold_nrs    .NE.nrs)
     .  CALL ER('LoadPIN','Stored grid is incompatable',*99)
c
c
c
      piniter = .TRUE.
c
c
c
c      CALL WN('LoadPIN','CALCULATING ALL OSM SOURCES!')
      CALL WN('LoadPIN','*NOT* CALCULATING OSM SOURCES')
c      CALL CalcSourceTerms
c
c
c
      nitersol  = ihold(1)
      rel_opt   = ihold(2)
      rel_frac  = rhold(3)
      rel_nstep = ihold(4)
      rel_niter = ihold(5)
      rel_step  = ihold(6)
      rel_iter  = ihold(7)
      rel_count = ihold(8)
      adp_opt   = ihold(9)

c...  Make sure pressure gauge specifications have not changed.  If they have
c     then discard the retrieved data:
      DO i1 = 1, MAXNAS
        IF (eirpgdat(i1,1).NE.eirpgdat2(i1,1).OR.
     .      eirpgdat(i1,2).NE.eirpgdat2(i1,2).OR.
     .      eirpgdat(i1,3).NE.eirpgdat2(i1,3)) THEN
          WRITE(0,*) 'DISCARDING PG DATA FOR GAUGE ',i1
          DO i2 = 1, MAXASD
            eirpgdat(i1,i2) = eirpgdat2(i1,i2)
          ENDDO
        ENDIF
      ENDDO
c
c
c
      DO ir = 1, nrs
        pininit(ir) = .TRUE.
      ENDDO
c
c
c
      rel_count    = rel_count + 1
      rel_viter(0) = 1


c      WRITE(0,*) 'MARK: IKBOUND12= ',ikbound(3,IKLO),ikbound(3,IKHI),
c     .   rel_count

c...Should make sure that the loaded and saved solutions have the same
c   version number, otherwise NULL arrays and such could be registered
c   as valid data by the following SaveSolution call:
      IF ((rel_opt.GT.1.AND.rel_count.EQ.0).OR.citersol.EQ.0) THEN
        WRITE(0     ,*) 'SAVING SOLUTION IN LOADPIN'
        WRITE(PINOUT,*) 'SAVING SOLUTION IN LOADPIN'
        CALL SaveSolution
        CALL DB('Done saving solution')
      ENDIF

      CALL OutputData(86,'After loading saved solution')

c      CALL OutputStratumData

c      IF (out_source.EQ.1) THEN
c        CALL SaveSources
c        CALL SavePlasma
c        CALL SaveGeometry
c      ENDIF


c      CALL WN('LoadPIN','Blanking momentum data')
c      CALL RZero(pinmp,MAXNKS*MAXNRS)
c      CALL RZero(osmmp,MAXNKS*MAXNRS)



      WRITE(0,*) 'MARK: IFLEXOPT4=',iflexopt(4)
      IF (iflexopt(4).EQ.2) THEN
        WRITE(0     ,*) 'REMOVING SOURCES ON OUTER HALF_RING 21'
        WRITE(PINOUT,*) 'REMOVING SOURCES ON OUTER HALF_RING 21'
        ir = 21
        DO ik = ikmids(ir)+1, nks(ir)
          osmmp (ik,ir) = 0.0
          pinrec(ik,ir) = 0.0
        ENDDO
        DO ik = 107, ikmids(ir)+1, -1
          pinion(ik+5,ir) = pinion(ik,ir)
          pinqe (ik+5,ir) = pinqe (ik,ir)
        ENDDO
      ENDIF
      IF (iflexopt(4).EQ.3) THEN
        WRITE(0     ,*) 'TRUNCATING CX MOMENTUM SOURCE'
        WRITE(PINOUT,*) 'TRUNCATING CX MOMENTUM SOURCE'
        DO ir = irsep, nrs
          DO ik = 1, ikmids(ir)
            osmmp(ik,ir) = MAX(osmmp(ik,ir),0.0)
            pinmp(ik,ir) = MAX(pinmp(ik,ir),0.0)
          ENDDO
          DO ik = ikmids(ir)+1, nks(ir)
            osmmp(ik,ir) = MIN(osmmp(ik,ir),0.0)
            pinmp(ik,ir) = MIN(pinmp(ik,ir),0.0)
          ENDDO
        ENDDO
      ENDIF
      IF (iflexopt(4).EQ.4) THEN
c...    Transition from "standard" method of precribing the shape of the
c       density peak in the prescription region, to the using the momentum
c       loss calculated by EIRENE.  Here, I need to overwrite the momentum
c       loss term in the prescription region so that, initially, the
c       shape of the density peak is preserved, and so the new shape can be
c       relaxed into:

        DO ir = irsep, nrs
          IF (GetModel(IKLO,ir).EQ.24) THEN
            ik1 = ikbound(ir,IKLO)

            DO ik = ik1-1, 1, -1
              p1 = CalcPressure(knbs (ik+1,ir),ktebs(ik+1,ir),
     .                          ktibs(ik+1,ir),kvhs (ik+1,ir))
              p  = CalcPressure(knbs (ik  ,ir),ktebs(ik  ,ir),
     .                          ktibs(ik  ,ir),kvhs (ik  ,ir))

              pl1 = osmmp(ik+1,ir) * (kss(ik+1,ir) - ksb(ik,ir)) / ECH

              IF (ik.EQ.ik1-1.AND.pl1.GT.0.5*(p1-p)) THEN
                WRITE(PINOUT,*) 'FIXING OSMMP'

                pl1 = 0.5 * (p1 - p)

                osmmp(ik+1,ir) = pl1 / (kss(ik+1,ir) - ksb(ik,ir)) * ECH
              ENDIF

              pl  = (p1 - p) - pl1

              osmmp(ik,ir) = pl / (ksb(ik,ir) - kss(ik,ir)) * ECH
              pinmp(ik,ir) = osmmp(ik,ir)

              WRITE(PINOUT,'(A,2I4,1P,6E10.2,0P)')
     .          'PLC : ',ik,ir,p1,p,(p1-p),pl1,pl,osmmp(ik,ir)
            ENDDO

            p1 = CalcPressure(knbs (ik1,ir),ktebs(ik1,ir),
     .                        ktibs(ik1,ir),kvhs (ik1,ir))

            p  = CalcPressure(knbs (1,ir),ktebs(1,ir),
     .                        ktibs(1,ir),kvhs (1,ir))

            CALL CalcIntegral3(osmmp,1,ik1,ir,rdum1,3)

            WRITE(PINOUT,'(A,1P,3E12.4,0P)') ' FINAL: ',
     .        p1,rdum1/ECH,p

            WRITE(0,*)

          ENDIF
        ENDDO
      ENDIF
      IF (iflexopt(4).EQ.5.AND.rflexopt(3).LT.1.0) THEN
        WRITE(0     ,*) 'SCALING RECOMBINATION SINK'
        WRITE(PINOUT,*) 'SCALING RECOMBINATION SINK'
        DO ir = irsep, nrs
          DO ik = 1, nks(ir)
            pinrec(ik,ir) = pinrec(ik,ir) * (rflexopt(3)**3.0)
          ENDDO
        ENDDO
      ENDIF
      IF (iflexopt(4).EQ.6) THEN
        WRITE(0     ,*) 'LIMITED HARDCODED SCALING RECOMBINATION SINK'
        WRITE(PINOUT,*) 'LIMITED HARDCODED SCALING RECOMBINATION SINK'
        ir = 3
        DO ik = 1, nks(ir)
          pinrec(ik,ir) = pinrec(ik,ir) * 0.50
        ENDDO
c        ir = 10
c        DO ik = 1, nks(ir)
c          pinrec(ik,ir) = pinrec(ik,ir) * 0.50
c        ENDDO
c        DO ir = irsep, irwall-1
c          DO ik = 1, nks(ir)
c            pinrec(ik,ir) = pinrec(ik,ir) * 0.75
c          ENDDO
c        ENDDO
      ENDIF
      IF (iflexopt(4).EQ.7.AND.rflexopt(4).GT.1.0) THEN
        WRITE(0     ,*) 'SCALING RECOMBINATION SINK'
        WRITE(PINOUT,*) 'SCALING RECOMBINATION SINK'
        DO ir = irsep, nrs
          DO ik = 1, nks(ir)
            pinrec(ik,ir) = pinrec(ik,ir) / rflexopt(4)
          ENDDO
        ENDDO
      ENDIF
      IF (iflexopt(4).EQ.8) THEN
        WRITE(0     ,*) 'CLEARING PINQI'
        WRITE(PINOUT,*) 'CLEARING PINQI'
        CALL RZero(pinqi,MAXNKS*MAXNRS)
      ENDIF
      IF (iflexopt(4).EQ.9) THEN
        WRITE(0,*) '>>>>>>>>>>>TEMP WHIPE<<<<<<<<<<<<'
        WRITE(0,*) '>>>>>>>>>>>TEMP WHIPE<<<<<<<<<<<<'
        WRITE(0,*) '>>>>>>>>>>>TEMP WHIPE<<<<<<<<<<<<'

c        pinion(ikbound(9,IKLO),ir) = pinion(ikbound(9,IKLO),ir) * 0.5
c        pinqe (ikbound(9,IKLO),ir) = pinqe (ikbound(9,IKLO),ir) * 0.5
c        pinion(ikbound(9,IKLO)+1,ir) = pinion(ikbound(9,IKLO)+1,ir)*0.5
c        pinqe (ikbound(9,IKLO)+1,ir) = pinqe (ikbound(9,IKLO)+1,ir)*0.5

c        DO ir = 1, nrs
        ir = 9
          DO ik = 1, ikmids(ir)
            IF (osm_dp6(ik,ir).NE.1.0)
     .            osm_dp6(ik,ir) = osm_dp6(ik,ir) * 0.60
          ENDDO
c          osm_peimul(IKLO,ir) = 0.1
c          osm_peimul(IKHI,ir) = 0.1
c        ENDDO
c        ir = 29
c          DO ik = 1, nks(ir)
c            osm_dp6(ik,ir) = 2.00
c          ENDDO
c          osm_peimul(IKLO,ir) = 0.1
c          osm_peimul(IKHI,ir) = 0.1
      ENDIF
      IF     (iflexopt(4).EQ.10) THEN
        CALL AnalyseDensityPeakWidth(PINOUT)
c        CALL FitDensityPeak
c        CALL CheckDensityLimit
      ELSEIF (iflexopt(4).EQ.11) THEN
        ir = 9

        osm_peimul(IKLO,ir) = 0.607878
        osm_peimul(IKHI,ir) = 1.410888

        DO ik = 1, nks(ir)
          osm_dp6(ik,ir) = 1.00
        ENDDO

        osm_dp6(ikbound(ik,IKLO)+1,ir) = -168.96
        osm_dp6(ikbound(ik,IKLO)+2,ir) = -0.33
        osm_dp6(ikbound(ik,IKLO)+3,ir) = -6.39

        osm_dp6(ikbound(ik,IKHI)-1,ir) = -975.56
        osm_dp6(ikbound(ik,IKHI)-2,ir) =  -12.21
        osm_dp6(ikbound(ik,IKHI)-3,ir) =   -6.73
        osm_dp6(ikbound(ik,IKHI)-4,ir) =   -0.81
      ELSEIF (iflexopt(4).EQ.12) THEN
        DO ir = irsep, nrs
          osm_peimul(IKLO,ir) = 0.5
          osm_peimul(IKHI,ir) = 0.5
        ENDDO

      ELSEIF (.FALSE..OR.iflexopt(4).EQ.13) THEN
        WRITE(0,*)
        WRITE(0,*)
        WRITE(0,*)
        WRITE(0,*) '*****************************************'
        WRITE(0,*) '* BLANKING PINASD AND PINBGK IN LOADPIN *'
        WRITE(0,*) '*****************************************'
        WRITE(0,*)
        WRITE(0,*)
        WRITE(0,*)
        indasd = 0
        CALL RZero(pinasd,MAXASCDAT*MAXASD2*MAXASS*2)
        CALL RZero(pinbgk,MAXNKS*MAXNRS*MAXBGK*MAXTOR)
      ENDIF
      WRITE(0,*) 'END OF IFLEXOPT(4)'



c      ir = 21
c      DO ik = 0.5*nks(ir) , nks(ir)
c        pinqe(ik,ir) = 0.83 * pinqe(ik,ir)
c      ENDDO


c
c     Need to load PINxxx2 arrays for subsequent source relaxation:
c
      IF (osm_powopt.EQ.1) CALL ReversePINQeMultiplier(-1,-1)

      CALL StoreSources(-1)

                  WRITE(0,*) 'pinqe: ',pinqe(1,8),pinqe2(1,8)

      IF (osm_powopt.EQ.1) CALL ApplyPINQeMultiplier(-1,-1)


      IF (rel_count.EQ.0) THEN
c        CALL OutputGrid(86,'After loading stored plasma')

        WRITE(67,'(A,I4,A)') 'ITERATION ',rel_count,' (LoadPIN)'
        WRITE(67,*)
        CALL OutputEIRENE(67,'LOADING SUPP. RAW FILES')
        WRITE(67,*)

        CALL DB('Writing line radiation data')
        CALL OutputLineRadiationData
      ENDIF
 

c      CALL WN('LoadPIN','Altering ionisation data')
c      CALL SetBounds
c      ir = 9
c      DO ik = nks(ir)/2, ikbound(ir,IKHI)
cc      DO ik = ikbound(ir,IKLO), nks(ir)/2
c        pinion(ik,ir) = 2.0 * pinion(ik,ir)
c        pinqe (ik,ir) = 2.0 * pinqe (ik,ir)
cc        pinion(ik,ir) = 0.5 * pinion(ik,ir)
cc        pinqe (ik,ir) = 0.5 * pinqe (ik,ir)
c      ENDDO

      IF (stopopt3.EQ.5) THEN
        CALL WN('LoadPIN','Altering momentum source data')
        CALL SetBounds
        DO ir = irsep, nrs
          IF (GetModel(IKLO,ir).EQ.22) THEN
            WRITE(0,*) 'MARK: BLANKING PINMP ON IKLO, RING= ',ir
            DO ik = 1, nks(ir)/2
              pinmp(ik,ir) = 0.0 * pinmp(ik,ir)
            ENDDO
          ENDIF
          IF (GetModel(IKHI,ir).EQ.22) THEN
            WRITE(0,*) 'MARK: BLANKING PINMP ON IKHI, RING= ',ir
            DO ik = nks(ir)/2+1, nks(ir)
              pinmp(ik,ir) = 0.0 * pinmp(ik,ir)
            ENDDO
          ENDIF
        ENDDO
        stopopt3 = 4
      ENDIF




c
c
c
c ... need to store osm_sympt data...
      IF     (osm_symopt.EQ.0) THEN
        DO ir = 1, nrs
          osm_sympt(ir) = ikmids(ir)
        ENDDO
      ELSEIF (osm_symopt.GT.2) THEN
        CALL WN('LoadPIN','Unsupported SYMPT option')
      ENDIF
c
c     For use in SOL 22:
c
      DO ir = 1, nrs
        DO ik = 1, nks(ir)
          oldknbs (ik,ir) = knbs (ik,ir)
          oldktebs(ik,ir) = ktebs(ik,ir)
          oldktibs(ik,ir) = ktibs(ik,ir)
          oldkvhs (ik,ir) = kvhs (ik,ir)
        ENDDO
      ENDDO

      CALL DB('Estimating opacity multiplier')
      IF (eiropacity.GT.0) CALL EstimateOpacityMultiplier      

c
c     Assign target quantities if requested:
c
      DO i1 = 1, nlpdato
        id = idds(INT(lpdato(i1,1)),2)

        IF (lpdato(i1,2).EQ.-1.0) lpdato(i1,2) = kteds(id)
        IF (lpdato(i1,3).EQ.-1.0) lpdato(i1,3) = ktids(id)
        IF (lpdato(i1,4).EQ.-1.0) THEN
          IF (lpdatsw.EQ.1) THEN
            lpdato(i1,4) = knds(id) * ECH * 9.79E+03 *
     .                     SQRT(0.5 * (kteds(id) + ktids(id)) *
     .                          (1.0 + rizb) / crmb)
          ELSE
            lpdato(i1,4) = knds(id)
          ENDIF
        ENDIF

      ENDDO

      DO i1 = 1, nlpdati
        id = idds(INT(lpdato(i1,1)),1)

        IF (lpdati(i1,2).EQ.-1.0) lpdati(i1,2) = kteds(id)
        IF (lpdati(i1,3).EQ.-1.0) lpdati(i1,3) = ktids(id)
        IF (lpdati(i1,4).EQ.-1.0) THEN
          IF (lpdatsw.EQ.1) THEN
            lpdati(i1,4) = knds(id) * ECH * 9.79E+03 *
     .                     SQRT(0.5 * (kteds(id) + ktids(id)) *
     .                          (1.0 + rizb) / crmb)
          ELSE
            lpdati(i1,4) = knds(id)
          ENDIF
        ENDIF
      ENDDO

      DO i1 = 1, nlpdato2
        id = idds(INT(lpdato2(i1,1)),2)

        IF (lpdato2(i1,2).LT.0.0)
     .    lpdato2(i1,2) = kteds(id) * ABS(lpdato2(i1,2))

        IF (lpdato2(i1,3).LT.0.0)
     .    lpdato2(i1,3) = ktids(id) * ABS(lpdato2(i1,3))

        IF (lpdato2(i1,4).LT.0.0) THEN
          lpdato2(i1,4) = knds (id) * ABS(lpdato2(i1,4))

          IF (lpdatsw.EQ.1)
     .      lpdato2(i1,4) = lpdato2(i1,4) * ECH * 9.79E+03 *
     .                      SQRT(0.5 * (lpdato2(i1,2) + lpdato2(i1,3)) *
     .                           (1.0 + rizb) / crmb)
        ENDIF

        IF (lpdato2(i1,5).LT.0.0)
     .    lpdato2(i1,5) = kteds(id) * ABS(lpdato2(i1,5))

        IF (lpdato2(i1,6).LT.0.0)
     .    lpdato2(i1,6) = ktids(id) * ABS(lpdato2(i1,6))

        IF (lpdato2(i1,7).LT.0.0) THEN
          lpdato2(i1,7) = knds (id) * ABS(lpdato2(i1,7))

          IF (lpdatsw.EQ.1)
     .      lpdato2(i1,7) = lpdato2(i1,7) * ECH * 9.79E+03 *
     .                      SQRT(0.5 * (lpdato2(i1,2) + lpdato2(i1,3)) *
     .                           (1.0 + rizb) / crmb)
        ENDIF

      ENDDO

      DO i1 = 1, nlpdati2
        id = idds(INT(lpdati2(i1,1)),1)

        IF (lpdati2(i1,2).LT.0.0)
     .    lpdati2(i1,2) = kteds(id) * ABS(lpdati2(i1,2))

        IF (lpdati2(i1,3).LT.0.0)
     .    lpdati2(i1,3) = ktids(id) * ABS(lpdati2(i1,3))

        IF (lpdati2(i1,4).LT.0.0) THEN
          lpdati2(i1,4) = knds (id) * ABS(lpdati2(i1,4))

          IF (lpdatsw.EQ.1)
     .      lpdati2(i1,4) = lpdati2(i1,4) * ECH * 9.79E+03 *
     .                      SQRT(0.5 * (lpdato2(i1,2) + lpdato2(i1,3)) *
     .                           (1.0 + rizb) / crmb)
        ENDIF


        IF (lpdati2(i1,5).LT.0.0)
     .    lpdati2(i1,5) = kteds(id) * ABS(lpdati2(i1,5))

        IF (lpdati2(i1,6).LT.0.0)
     .    lpdati2(i1,6) = ktids(id) * ABS(lpdati2(i1,6))

        IF (lpdati2(i1,7).LT.0.0) THEN
          lpdati2(i1,7) = knds (id) * ABS(lpdati2(i1,7))

          IF (lpdatsw.EQ.1)
     .      lpdati2(i1,7) = lpdati2(i1,7) * ECH * 9.79E+03 *
     .                      SQRT(0.5 * (lpdato2(i1,2) + lpdato2(i1,3)) *
     .                           (1.0 + rizb) / crmb)
        ENDIF

      ENDDO

      IF (1.EQ.2) THEN
        WRITE(0,*)
        DO i1 = 1, nlpdato
          WRITE(0,*) (lpdato(i1,i2),i2=1,4)
        ENDDO
        WRITE(0,*)
        DO i1 = 1, nlpdati
          WRITE(0,*) (lpdati(i1,i2),i2=1,4)
        ENDDO
        WRITE(0,*)
        DO i1 = 1, nlpdato2
          WRITE(0,*) (lpdato2(i1,i2),i2=1,7)
        ENDDO
        WRITE(0,*)
        DO i1 = 1, nlpdati2
          WRITE(0,*) (lpdati2(i1,i2),i2=1,7)
        ENDDO
      ENDIF


c      IF (nrs.EQ.6) THEN
c        ir = 3
c
c        WRITE(0,*) '>>>>>>>>>>>TEMP WHIPE<<<<<<<<<<<<'
c        WRITE(0,*) '>>>>>>>>>>>TEMP WHIPE<<<<<<<<<<<<'
c        WRITE(0,*) '>>>>>>>>>>>TEMP WHIPE<<<<<<<<<<<<'
c
c        rdum1 = 0.5
c
c        DO ik = 1, nks(ir)
c          osm_dp6(ik,ir) = SIGN(1.0,osm_dp6(ik,ir)) *
c     .        ((1.0-rdum1)*ABS(osm_dp6(ik,ir)) + rdum1*1.00)
c        ENDDO
c      ENDIF


c       DO ir = irsep, nrs
c         WRITE(0,*) 'KTEDS: ',ir,kteds(idds(ir,2)),kteds(idds(ir,1))
c       ENDDO

      DO ir = 1, MAXNRS
        osm_peimul2(IKLO,ir) = DBLE(osm_peimul(IKLO,ir))
        osm_peimul2(IKHI,ir) = DBLE(osm_peimul(IKHI,ir))
        DO ik = 0, MAXNKS
          osmpmk2    (ik,ir) = DBLE(osmpmk(ik,ir))
          save_osmpmk(ik,ir) = DBLE(osmpmk(ik,ir))
        ENDDO

        IF (stopopt3.EQ.20.AND.
     .      ir.GE.27.AND.ir.LE.38.AND.
     .      osm_model(IKLO,ir).EQ.24.AND.
     .      osm_model(IKHI,ir).EQ.22) THEN
          scaleplateau(ir) = .TRUE.

          WRITE(0,*) 'LOADPIN: SETTING SCALEPLATAEU',ir

        ENDIF

      ENDDO

c...  Silly, but clean up NVERTP on boundary rings (SL, 19.11.2004):
      DO ir = 1, nrs
        IF (idring(ir).NE.BOUNDARY) CYCLE
        DO ik = 1, nks(ir)
          nvertp(korpg(ik,ir)) = 0
        ENDDO
      ENDDO

c      CALL OutputData(86,'After loading')
c      STOP 'sdfsdfsd'

      RETURN
95    STOP 'End OF LoadPIN'
96    STOP 'ERROR LoadPIN: Cannot find plasma file'
97    STOP 'ERROR LoadPIN: Cannot find source file'
98    STOP 'ERROR LoadPIN: Cannot find geomty file'
99    WRITE(0,'(5X,A,3I4)') 'ERROR1-3 = ',
     .  error1,error2,error3
      STOP ' '
      END


c
c =SECTION==============================================================
c
c ======================================================================
c
c subroutine: SaveSolution
c
      SUBROUTINE SaveSolution
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INCLUDE 'solparams'
      INCLUDE 'solcommon'
      INCLUDE 'solswitch'

      osmcfp(1,1) = 0.0
      pinion(1,1) = 0.0
      pinrec(1,1) = 0.0
      osmcfe(1,1) = 0.0
      osmcfi(1,1) = 0.0
      osmpei(1,1) = 0.0
      pinqi (1,1) = 0.0
      pinqe (1,1) = 0.0
      osmqe (1,1) = 0.0

c...prad!
      IF (switch(SWGPERP).NE. 2.0.AND.
     .    switch(SWGPERP).NE. 7.0.AND.
     .    switch(SWGPERP).NE. 8.0) osmcfp(1,1) = LO
      IF (switch(SWION  ).EQ. 0.0) pinion(1,1) = LO
      IF (switch(SWRECOM).NE. 1.0) pinrec(1,1) = LO
      IF (switch(SWPOW  ).NE.13.0) osmcfe(1,1) = LO
      IF (switch(SWPOW  ).NE.13.0) osmcfi(1,1) = LO
      IF (switch(SWPEI  ).NE. 4.0.AND.
     .    switch(SWPEI  ).NE. 5.0) osmpei(1,1) = LO
      IF (osmmock.EQ.0)            osmpmk(1,1) = LO
      IF (switch(SWPCX  ).NE. 2.0) pinqi (1,1) = LO
      IF (switch(SWPHELP).NE. 2.0.AND.
     .    switch(SWPHELP).NE. 3.0) pinqe (1,1) = LO
      IF (switch(SWRECOM).NE. 1.0.OR.
     .    osm_recopt     .EQ. 0)   osmqe (1,1) = LO

      CALL SaveSources
      CALL SavePlasma
      CALL SaveGeometry

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: SaveVacuumGrid
c
      SUBROUTINE ProcessVacuumGrid(mode)
      IMPLICIT   none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      REAL      version,rdum1,rand
      INTEGER   fp,v1,s1,cell,idum1
      CHARACTER mode*4

      version = 1.0

      IF     (mode.EQ.'LOAD') THEN

        fp=98
        OPEN(UNIT=fp,FILE='vac-grid.dat',FORM='UNFORMATTED',
     .       STATUS='OLD',ERR=96)
 
        READ (fp) rdum1,asc_ncell,idum1,ascncut

        IF (asc_ncell.GT.MAXASC) 
     .    CALL ER('ProcessVacuumGrid','Unable to load grid.  Change'//
     .                                'MAXASC to:',*95)
        IF (asc_ncell*ascncut.GT.MAXASCDAT) 
     .    CALL ER('ProcessVacuumGrid','Unable to load grid.  Change'//
     .                                'MAXASCDAT to:',*93)
        IF (ascncut.GT.MAXASC3D) 
     .    CALL ER('ProcessVacuumGrid','Unable to load grid.  Change'//
     .                                'MAXASC3D to:',*92)
        IF (idum1.NE.MAXVACREGION)
     .    CALL ER('ProcessVacuumGrid','Unable to load grid.  Change '//
     .                                'MAXVACREGION to:',*94)

        READ (fp) 
     .    asc_nregion,asc_3Dmode,asccode,
     .    (asc_cell  (cell),asc_region(cell),asc_nvp(cell),
     .     ascnvertex(cell),
     .     (asc_link  (s1,cell),s1=1,4),
     .     (asc_grid  (s1,cell),s1=1,2),
     .     (ascvertex (s1,cell),s1=1,2*ascnvertex(cell)),
     .     (asc_rvp   (s1,cell),
     .      asc_zvp   (s1,cell),s1=1,8),
     .     cell=1,asc_ncell),
     .    (asc_vol(cell),
     .     cell=1,asc_ncell*ascncut),
     .    (asc_nvp3D (cell),asc_zmin3D(cell),
     .     asc_zmax3D(cell),
     .     (asc_link3D(s1,cell),s1=1,6),
     .     ((asc_xvp3D(s1,v1,cell),
     .       asc_yvp3D(s1,v1,cell),
     .       asc_zvp3D(s1,v1,cell),v1=1,8),s1=1,6),
     .     cell=1,ascncut),
     .    (asc_rstart(s1),asc_rend(s1),s1=1,10),
     .    (vacregion(s1),s1=1,MAXVACREGION)

        CLOSE(fp)

c        WRITE(0,*) 'LOADING VACUUM GRID, ASCCODE=',asccode

      ELSEIF (mode.EQ.'SAVE') THEN

c...    Generate code for vacuum grid, to make sure that remains synchronized with 
c       the connection map generated in EIRENE:
        CALL Surand2(0.0D0,0,rand)
        asccode = NINT(rand*1.0E+7)

c        WRITE(0,*) 'ASCCODE=',asccode

        fp=98
        OPEN(UNIT=fp,FILE='vac-grid.dat',FORM='UNFORMATTED',
     .       STATUS='REPLACE',ERR=96)
 
        WRITE(fp) version,asc_ncell,MAXVACREGION,ascncut

        WRITE(fp) 
     .    asc_nregion,asc_3Dmode,asccode,
     .    (asc_cell  (cell),asc_region(cell),asc_nvp(cell),
     .     ascnvertex(cell),
     .     (asc_link  (s1,cell),s1=1,4),
     .     (asc_grid  (s1,cell),s1=1,2),
     .     (ascvertex (s1,cell),s1=1,2*ascnvertex(cell)),
     .     (asc_rvp   (s1,cell),
     .      asc_zvp   (s1,cell),s1=1,8),
     .     cell=1,asc_ncell),
     .    (asc_vol(cell),
     .     cell=1,asc_ncell*ascncut),
     .    (asc_nvp3D (cell),asc_zmin3D(cell),
     .     asc_zmax3D(cell),
     .     (asc_link3D(s1,cell),s1=1,6),
     .     ((asc_xvp3D(s1,v1,cell),
     .       asc_yvp3D(s1,v1,cell),
     .       asc_zvp3D(s1,v1,cell),v1=1,8),s1=1,6),
     .     cell=1,ascncut),
     .    (asc_rstart(s1),asc_rend(s1),s1=1,10),
     .    (vacregion(s1),s1=1,MAXVACREGION)

        CLOSE(fp)

      ELSE
        CALL ER('SaveVacuumGrid','Unknown mode',*99)
      ENDIF


      RETURN
92    WRITE(0,*) ascncut
      STOP
93    WRITE(0,*) asc_ncell*ascncut
      STOP
94    WRITE(0,*) idum1
      STOP
95    WRITE(0,*) asc_ncell,ascncut
      STOP
96    CALL ER('ProcessVacuumGrid','Unable to open vacuum grid file',*99)
99    STOP
      END
c
c ======================================================================
c
c subroutine: SaveGeometry
c
      SUBROUTINE SaveGeometry
      IMPLICIT   none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      REAL    slver
      INTEGER fp,ik,ir,id,ii,i1,init

      DATA init /0/
      SAVE
c
c     Only store geometry data on the first iteration, and when
c     the grid adaptaion code is being executed:
c

      IF (adp_opt.EQ.0.AND.init.EQ.1) RETURN

      init = 1

      slver = 1.08
      fp    = 98
      OPEN(UNIT=fp,FILE='geomty.dat',FORM='UNFORMATTED',
     .     POSITION='APPEND')

      WRITE(fp) slver
      WRITE(fp)
     .  nitersol,rel_opt ,rel_frac ,rel_nstep,rel_niter,
     .  rel_step,rel_iter,rel_count,adp_opt

      WRITE(fp)
     .  nrs,nds,ndsin,irwall,irtrap,irbreak,nbr,ikto,ikti,npolyp,
     .  vpolyp,vpolmin,r0,z0,rxp,zxp,dthetg

      WRITE(fp) (
     .    nks(ir),ikmids(ir),ksmaxs(ir),idds(ir,1),idds(ir,2),
     .    idring(ir),ikto2(ir),ikti2(ir),
     .    rho(ir,IN14),rho(ir,CELL1),rho(ir,OUT23),
     .  ir=1,nrs),(
     .    thetat(id),ikds(id),irds(id),sepdist(id),sepdist2(id),
     .  id=1,nds),((
     .    rs    (ik,ir),zs    (ik,ir),kss   (ik,ir),kps  (ik,ir),
     .    kbacds(ik,ir),kfords(ik,ir),thetag(ik,ir),korpg(ik,ir),
     .    bratio(ik,ir),kbfs  (ik,ir),
     .    nvertp(MAX(1,korpg(ik,ir))),
     .    (rvertp(ii,MAX(1,korpg(ik,ir))),
     .     zvertp(ii,MAX(1,korpg(ik,ir))),
     .       ii=1,nvertp(MAX(1,korpg(ik,ir)))),
     .  ik=1,nks(ir)),ir=1,nrs),((
     .    ksb(ik,ir),kpb(ik,ir),
     .  ik=0,nks(ir)),ir=1,nrs),
     .  nvesm,nvesp,(rvesm(i1,1),rvesm(i1,2),zvesm(i1,1),zvesm(i1,2),
     .               jvesm(i1),i1=1,nvesm+nvesp)

      CLOSE(fp)

      RETURN
      END
c
c ======================================================================
c
c subroutine: SavePlasma
c
      SUBROUTINE SavePlasma
      IMPLICIT   none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      REAL    slver
      INTEGER fp,ik,ir,id,i1,count

      DATA count /-1/

      SAVE


      count = count + 1
      WRITE(PINOUT,'(1X,A,I3,A)') 'Saving plasma (',count,')'
      IF (sloutput) 
     .  WRITE(0     ,'(1X,A,I3,A)') 'Saving plasma (',count,')'
c      IF (rel_opt.GT.0)
c     .WRITE(0     ,'(1X,A,I3,A)') 'Saving plasma (',count,')'


      slver = 1.15
c      fp    = PINOUT3
      fp    = 98
      OPEN(UNIT=fp,FILE='plasma.dat',FORM='UNFORMATTED',
     .     POSITION='APPEND')

      WRITE(fp) slver
      WRITE(fp)
     .  nitersol,rel_opt ,rel_frac ,rel_nstep,rel_niter,
     .  rel_step,rel_iter,rel_count,adp_opt

      WRITE(fp) nds,(
     .    kteds(id),ktids(id),knds(id),kvds(id),keds(id),
     .  id=1,nds),
     .  nrs,(nks(ir),cmachno(ir,1),cmachno(ir,2),
     .      s28ionfrac(IKLO,ir),s28recfrac(IKLO,ir),s28momfrac(IKLO,ir),
     .      s28ionfrac(IKHI,ir),s28recfrac(IKHI,ir),s28momfrac(IKHI,ir),
     .    rel_hte  (ir),rel_hti  (ir),rel_hne  (ir)  ,rel_deltati(ir),
     .    rel_dirtg(ir),rel_symfr(ir),rel_prbfr(1,ir),rel_prbfr  (2,ir),
     .    osm_sympt(ir),osm_peimul(1,ir),osm_peimul(2,ir),
     .    (rel_hproe(i1,ir),rel_hproi(i1,ir),i1=1,3),
     .    rel_qemul(IKLO,ir),rel_qemul(IKHI,ir),
     .    osm_model(IKLO,ir),osm_model(IKHI,ir),
     .    osm_code (IKLO,ir),osm_code (IKHI,ir),     
     .    ikbound  (ir,IKLO),ikbound  (ir,IKHI),
     .   (kpress (ik,ir,1),kpress(ik,ir,2),
     .    ktebs  (ik,ir)  ,ktibs (ik,ir)  ,knbs(ik,ir),kvhs(ik,ir),
     .    kes    (ik,ir)  ,
     .    osm_dp6(ik,ir)  ,
     .  ik=1,nks(ir)),ir=1,nrs),
     .  0,s21_ndatai,(s21_datai(i1,7),i1=1,s21_ndatai),
     .    s21_ndatao,(s21_datao(i1,7),i1=1,s21_ndatao),
     .  tarshift(IKLO),tarshift(IKHI),
     .  te_mult_o,ti_mult_o,n_mult_o,te_mult_i,ti_mult_i,n_mult_i,
     .  s28ionset,s28recset,s28momset,
     .  s28ionsetpfz,s28recsetpfz,s28momsetpfz

      CLOSE(fp)

      RETURN
      END
c
c ======================================================================
c
c subroutine: SaveSources
c
      SUBROUTINE SaveSources
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      REAL    slver
      INTEGER fp,ik,ir,id,i1,i2,i3,i4,mode

      slver = 1.21
      fp    = PINOUT2
c      fp    = 98
      OPEN(UNIT=fp,FILE='source.dat',FORM='UNFORMATTED',
     .     POSITION='APPEND')

        WRITE(fp) slver

        WRITE(fp)
     .    nitersol,rel_opt ,rel_frac ,rel_nstep,rel_niter,
     .    rel_step,rel_iter,rel_count,adp_opt

        WRITE(fp) nrs,(nks(ir),(
     .      pinion (ik,ir),pinrec(ik,ir),pinqe   (ik,ir),pinqi (ik,ir),
     .      pinatom(ik,ir),pinmol(ik,ir),pinalpha(ik,ir),pinmp (ik,ir),
     .      pinena (ik,ir),osmpei(ik,ir),osmcfp  (ik,ir),osmcfi(ik,ir),
     .      osmcfe (ik,ir),osmmp (ik,ir),osmqe   (ik,ir),
     .      osmcve (ik,ir),osmcvi(ik,ir),osmcde  (ik,ir),osmcdi(ik,ir),
     .      osmion (ik,ir),osmrec(ik,ir),
     .      (osmcfpflx(ik,ir,i1),i1=1,5),
     .      pinqir (ik,ir),pinqer(ik,ir),
     .      6      ,(pinline(ik,ir,i1,H_BALPHA),
     .               pinline(ik,ir,i1,H_BGAMMA),i1=1,6      ),
     .      3,MAXSTRATA,
     .        ((pinstrata(ik,ir,i2,i1),i2=1,3),i1=1,MAXSTRATA),
     .      NMOMCHA         ,(pinploss(ik,ir,i1),i1=1,NMOMCHA),
     .      MAXBGK*eirnsdtor,(pinbgk  (ik,ir,i1),i1=1,MAXBGK*eirnsdtor),
     .    ik=1,nks(ir)),
     .                (nks(ir),(osmpmk(ik,ir)),
     .    ik=1,nks(ir)+1),
     .    ir=1,nrs),
     .    MAXNAS,MAXASD,((eirpgdat(i1,i2),i1=1,MAXNAS),i2=1,MAXASD),
     .    indasd,asc_ncell*ascncut,MAXASD2,MAXASS,2,
     .    ((((pinasd(i1,i2,i3,i4),i1=1,asc_ncell*ascncut),i2=1,MAXASD2),
     .                            i3=1,MAXASS           ),i4=1,2)

      CLOSE(fp)

      RETURN
      END
c
c =SECTION==============================================================
c
c ======================================================================
c
c subroutine: UpdateTargets
c
c ... make sure that nlpdati,o is not zero if rel_opt  = 2...
c
      SUBROUTINE UpdateTargets(iitersol)
      IMPLICIT   none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER iitersol

      INTEGER i1,i2,i3,i4,i5,sav_niter,ir,in
      REAL    frac,lpval1,lpval2,lpdat1,lpdat2,ratio,log1,log2,log3

      DATA sav_niter /-1/
      SAVE

      IF (iitersol.NE.-1.OR.
     .    (nlpdati2.NE.0.OR.nlpdato2.NE.0))
     .  WRITE(0,*) 'HERE IN UPDATE TARGETS',iitersol


      IF ((tarshift(IKLO).NE.0.0.AND.nlpdato2.LT.irwall-2).OR.
     .    (tarshift(IKHI).NE.0.0.AND.nlpdati2.LT.irwall-2)) 
     .  CALL ER('UpdateTargets','Target shift requested but the '//
     .          'reference target data is insufficient',*99)

      IF (sav_niter.EQ.-1) sav_niter = rel_niter

      IF (outmode.EQ.3) CALL MS('UpdateTargets','Relaxing BC''s...')

      WRITE(PINOUT,'(A,6I6)') 
     .  ' UPDATE TARGETS: ',iitersol,rel_iter,rel_niter,
     .                               rel_step,rel_nstep,citersol
c      WRITE(0     ,'(A,6I6)') 
c     .  ' UPDATE TARGETS: ',iitersol,rel_iter,rel_niter,
c     .                               rel_step,rel_nstep,citersol


      IF (.FALSE.) THEN
        IF (citersol.EQ.0) THEN
          frac = 0.0

          IF (rel_ndata.GT.0) CALL SetupRelaxation
        ELSE
          IF     (iitersol.EQ.-1) THEN
            frac      = 0.0
            rel_iter  = 0
            rel_step  = 0
          ELSEIF (iitersol.EQ.1) THEN
            frac      = 0.0
            rel_iter  = 1
            rel_step  = 0
            IF (rel_ndata.GT.0) CALL SetupRelaxation
          ELSEIF (rel_nstep.EQ.1) THEN
            IF (rel_step.EQ.0) rel_iter = 0
            frac     = 0.0
            rel_iter = rel_iter + 1
            rel_step = 1
            IF (rel_ndata.GT.0) CALL SetupRelaxation
          ELSE
            IF (rel_step.EQ.0.OR.rel_iter.GE.rel_niter) THEN
c            IF (rel_step.EQ.0.OR.rel_iter.GE.rel_niter) THEN
              rel_viter(rel_step) = rel_iter
              rel_iter            = 1
              rel_step            = rel_step + 1

              WRITE(PINOUT,'(2X,A,3I4)') 'NEW STEP = ',
     .          rel_step,rel_iter,rel_viter(rel_step-1)

              IF (rel_ndata.GT.0) THEN
                CALL SetupRelaxation
              ELSE
                rel_niter = sav_niter
              ENDIF
            ELSE
              rel_iter = rel_iter + 1
            ENDIF
            frac = REAL(rel_step - 1) / REAL(rel_nstep - 1)
            IF (rel_pace.GT.0.0) frac = frac**(1.0 / rel_pace)
          ENDIF
        ENDIF
      ELSEIF     (iitersol.EQ.-1) THEN
        frac      = 0.0
        rel_iter  = 0
        rel_step  = 0
      ELSE
c...   REL_ITER, REL_STEP, etc. updated in SetupIteration:
        frac = REAL(rel_step - 1) / REAL(rel_nstep - 1)
        IF (rel_pace.GT.0.0) frac = frac**(1.0 / rel_pace)
      ENDIF

      WRITE(PINOUT,'(A,5I5,F10.5)')
     .  ' UPDATE TARGETS: ',iitersol,rel_iter,rel_niter,
     .                               rel_step,rel_nstep,frac
c      IF (outmode.EQ.2)
c     .  WRITE(0     ,'(A,3I5,F10.5)')
c     .    'UPDATE TARGETS: ',iitersol,rel_step,rel_iter,frac
c
c     Inner (outer CMOD) target:
c

      DO i1 = 1, nlpdati2
        DO i2 = 1, nlpdati
          IF (lpdati(i2,1).EQ.lpdati2(i1,1)) THEN
            i4 = i2
            GOTO 10
          ENDIF
        ENDDO

        CALL MS('UpdateTargets','Adding inner (CMOD outer) data point')
        nlpdati = nlpdati + 1
        i4      = nlpdati

        lpdati(i4,1) = lpdati2(i1,1)

10      CONTINUE

c...FSP Te multiplier
        IF (stopopt2.NE.8) THEN
         lpval1 = lpdati2(i1,8)+rel_bound1*(lpdati2(i1,9)-lpdati2(i1,8))
         lpval2 = lpdati2(i1,8)+rel_bound2*(lpdati2(i1,9)-lpdati2(i1,8))
         rel_mfsp(IKHI,INT(lpdati2(i1,1))) = lpval1+frac*(lpval2-lpval1)
c        WRITE(0,*)'MFSP HI ',INT(lpdati2(i1,1)),
c     .            rel_mfsp(IKHI,INT(lpdati2(i1,1)))
        ENDIF

        DO i3 = 2, 4
          lpval1 = lpdati2(i1,i3) +
     .             rel_bound1 * (lpdati2(i1,i3+3) - lpdati2(i1,i3))
          lpval2 = lpdati2(i1,i3) +
     .             rel_bound2 * (lpdati2(i1,i3+3) - lpdati2(i1,i3))

          lpdat2        = lpdati(i4,i3)

          IF (rel_pace.EQ.-1.0) THEN
c...        Logarithmic source variation:
            log1 = LOG10(lpval1)
            log2 = LOG10(lpval2)
            log3 = log1 + frac * (log2 - log1)
            lpdati(i4,i3) = 10.0**log3
c            WRITE(0,'(A,I6,3F10.2,1P,E10.2,0P)') 
c     .        '-->',i3,log1,log2,log3,lpdati(i4,i3)
          ELSE
            lpdati(i4,i3) = lpval1 + frac * (lpval2 - lpval1)
          ENDIF

          IF (osm_matchs.EQ.1.AND.rel_step.GT.1) THEN
c
c
c ... not sure if I should be using a ratio or increment value...?
c
            ir = INT(lpdati(i4,1))

            DO in = 1, nbgplas
              IF (INT(bgplasopt(in,1)).LE.ir  .AND.
     .            INT(bgplasopt(in,2)).GE.ir  .AND.
     .                bgplasopt(in,3) .EQ.3.0 .AND.
     .                bgplasopt(in,5) .EQ.22.0) THEN

                ratio = lpdati(i4,i3) / lpdat2

                IF (rel_iter.EQ.1)
     .            WRITE(0,*) 'Updating inner ',ir,' by ',ratio

                DO i5 = 1, nlpdato
                  IF (INT(lpdato(i5,1)).EQ.ir)
     .              lpdato(i5,i3) = lpdato(i5,i3) * ratio
                ENDDO
              ENDIF
            ENDDO
          ENDIF


c          IF (outmode.EQ.2)
c     .      write(0,*) 'data -> ',
c     .        i4,i3,lpdati2(i1,i3), frac,lpdati2(i1,i3+3)
c          write(0     ,'(A,2I6,1P,4E12.4)') 'data -> ',
c     .        i4,i3,lpdati2(i1,i3), frac,lpdati2(i1,i3+3),lpdati(i4,i3)
c          write(PINOUT,'(A,2I6,1P,4E12.4)') 'data -> ',
c     .        i4,i3,lpdati2(i1,i3), frac,lpdati2(i1,i3+3),lpdati(i4,i3)

        ENDDO

        lpdati(i4,2) = lpdati(i4,2) * te_mult_i
        lpdati(i4,3) = lpdati(i4,3) * ti_mult_i
        lpdati(i4,4) = lpdati(i4,4) *  n_mult_i
      ENDDO

      IF (stopopt.EQ.23) STOP 'After updating outer CMOD target'
c
c     Outer (inner CMOD) target:
c
      DO i1 = 1, nlpdato2
        DO i2 = 1, nlpdato
          IF (lpdato(i2,1).EQ.lpdato2(i1,1)) THEN
            i4 = i2
            GOTO 20
          ENDIF
        ENDDO

        CALL MS('UpdateTargets','Adding outer (CMOD inner) data point')
        nlpdato = nlpdato + 1
        i4      = nlpdato

        lpdato(i4,1) = lpdato2(i1,1)

20      CONTINUE

c...FSP Te multiplier
        IF (stopopt2.NE.8) THEN
         lpval1 = lpdato2(i1,8)+rel_bound1*(lpdato2(i1,9)-lpdato2(i1,8))
         lpval2 = lpdato2(i1,8)+rel_bound2*(lpdato2(i1,9)-lpdato2(i1,8))
         rel_mfsp(IKLO,INT(lpdato2(i1,1))) = lpval1+frac*(lpval2-lpval1)
c          WRITE(0,*)'MFSP LO ',INT(lpdato2(i1,1)),
c     .      rel_mfsp(IKLO,INT(lpdato2(i1,1)))
        ENDIF

        DO i3 = 2, 4
          lpval1 =  lpdato2(i1,i3) +
     .              rel_bound1 * (lpdato2(i1,i3+3) - lpdato2(i1,i3))
          lpval2 =  lpdato2(i1,i3) +
     .              rel_bound2 * (lpdato2(i1,i3+3) - lpdato2(i1,i3))

          lpdat2        = lpdato(i4,i3)
          IF (rel_pace.EQ.-1.0) THEN
c...        Logarithmic source variation:
            log1 = LOG10(lpval1)
            log2 = LOG10(lpval2)
            log3 = log1 + frac * (log2 - log1)
            lpdato(i4,i3) = 10.0**log3
c            WRITE(0,'(A,I6,3F10.2,1P,E10.2,0P)') 
c     .        '-->',i3,log1,log2,log3,lpdati(i4,i3)
          ELSE
            lpdato(i4,i3) = lpval1 + frac * (lpval2 - lpval1)
          ENDIF

          IF (osm_matchs.EQ.2.AND.rel_step.GT.1) THEN
c
c
c ... not sure if I should be using a ratio or increment value...?
c
            ir = INT(lpdato(i4,1))

            DO in = 1, nbgplas
              IF (INT(bgplasopt(in,1)).LE.ir  .AND.
     .            INT(bgplasopt(in,2)).GE.ir  .AND.
     .                bgplasopt(in,3) .EQ.3.0 .AND.
     .                bgplasopt(in,5) .EQ.22.0) THEN

                ratio = lpdato(i4,i3) / lpdat2

                IF (rel_iter.EQ.1)
     .            WRITE(0,*) 'Updating outer ',ir,' by ',ratio

                DO i5 = 1, nlpdati
                  IF (INT(lpdati(i5,1)).EQ.ir)
     .              lpdati(i5,i3) = lpdati(i5,i3) * ratio
                ENDDO
              ENDIF
            ENDDO
          ENDIF
        ENDDO

        ir = NINT(lpdato(i2,1))

        IF (.FALSE..AND.
     .      (osm_model(IKLO,ir).EQ.24.OR.osm_model(IKLO,ir).EQ.28)) THEN
c...      This check needs to be added to wherever the target multipliers
c         are used:
          WRITE(0,*) 'NOT APPLYING TARGET DATA MULTIPLIERS TO INNER',ir
        ELSE
          lpdato(i4,2) = lpdato(i4,2) * te_mult_o
          lpdato(i4,3) = lpdato(i4,3) * ti_mult_o
          lpdato(i4,4) = lpdato(i4,4) *  n_mult_o
        ENDIF
      ENDDO

      WRITE(PINOUT,*)
      WRITE(PINOUT,*) 'Updating target data:'
        
      IF (.TRUE..OR.osm_mode.GE.2.OR.outmode.GE.2) THEN
        DO i1 = 1, nlpdato
          WRITE(PINOUT,'(5X,20(F5.1,2X,2F10.4,2X,1P,E15.7))')
     .      (lpdato(i1,i2),i2=1,4)
c          WRITE(0     ,'(5X,20(F5.1,2X,2F10.4,2X,1P,E15.7))')
c     .      (lpdato(i1,i2),i2=1,4)
         ENDDO
        WRITE(PINOUT,*)
        DO i1 = 1, nlpdati
          WRITE(PINOUT,'(5X,20(F5.1,2X,2F10.4,2X,1P,E15.7))')
     .      (lpdati(i1,i2),i2=1,4)
c          WRITE(0     ,'(5X,20(F5.1,2X,2F10.4,2X,1P,E15.7))')
c     .      (lpdati(i1,i2),i2=1,4)
        ENDDO
      ENDIF


c      IF (iitersol.NE.-1) STOP

      IF ((tarshift(IKLO).NE.0.0.OR.tarshift(IKHI).NE.0.0).AND. 
     .    iitersol.NE.-1) CALL ShiftTargetData


      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: SetupRelax
c
crelax
      SUBROUTINE SetupRelaxation
      IMPLICIT   none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER i1,ir,nvac,i2
      LOGICAL status,init

      status = .FALSE.

      nvac = 0

c      WRITE(0     ,*) '**** REL_STEP=',rel_step,rel_ndata
      WRITE(PINOUT,*) '**** REL_STEP=',rel_step,rel_ndata

      init = .TRUE.

      DO i1 = 1, rel_ndata

c...    RELMODE=20 option, where RELMODE can be changed on the fly:
        IF (rel_data(i1,1).EQ.-2.0) relmode = NINT(rel_data(i1,2))

        IF (NINT(rel_data(i1,1)).EQ.rel_step.OR.
     .      NINT(rel_data(i1,1)).EQ.-1) THEN

          status = .TRUE.

          IF (init.AND.relmode.EQ.2) THEN
            CALL ISet(supflx,2*MAXNRS,1)
            init = .FALSE.
          ENDIF

          IF     (relmode.EQ.0) THEN
            rel_niter = NINT(rel_data(i1,2))
            rel_opt   = NINT(rel_data(i1,3))
            rel_frac  =     rel_data(i1,4)
            adp_opt   = NINT(rel_data(i1,5))
            eirtime   = NINT(rel_data(i1,6))
            rel_tol   =     rel_data(i1,7)
          ELSEIF (relmode.EQ.1) THEN
            te_mult_o      = rel_data(i1,2)
            ti_mult_o      = rel_data(i1,3)
            n_mult_o       = rel_data(i1,4) 
            te_mult_i      = rel_data(i1,5)
            ti_mult_i      = rel_data(i1,6)
            n_mult_i       = rel_data(i1,7)
            tarshift(IKLO) = rel_data(i1,8)
            tarshift(IKHI) = rel_data(i1,9)
          ELSEIF (relmode.EQ.2) THEN
c...        Set targets whose flux will NOT be suppressed when calling
c           EIRENE (SUPFLX=0):
            DO ir = NINT(rel_data(i1,2)), NINT(rel_data(i1,3))
              IF (idring(ir).EQ.BOUNDARY) CYCLE
              IF (rel_data(i1,4).EQ.1.0.OR.
     .            rel_data(i1,4).EQ.3.0) supflx(IKLO,ir) = 0
              IF (rel_data(i1,4).EQ.2.0.OR.
     .            rel_data(i1,4).EQ.3.0) supflx(IKHI,ir) = 0
            ENDDO
          ELSEIF (relmode.EQ.3) THEN
            te_mult_o      = rel_data(i1,2)
            ti_mult_o      = rel_data(i1,3)
            n_mult_o       = rel_data(i1,4) 
            te_mult_i      = rel_data(i1,5)
            ti_mult_i      = rel_data(i1,6)
            n_mult_i       = rel_data(i1,7)
            CALL ISet(supflx,2*MAXNRS,0) 
            IF (rel_data(i1,8).NE.0.0) THEN
              CALL ISet(supflx,2*MAXNRS,1) 
              DO ir = 1, nrs
                IF (rel_data(i1,8).EQ.1.0.OR.
     .              rel_data(i1,8).EQ.3.0) supflx(IKLO,ir) = 0
                IF (rel_data(i1,8).EQ.2.0.OR.
     .              rel_data(i1,8).EQ.3.0) supflx(IKHI,ir) = 0
              ENDDO
            ENDIF
            eirtime = NINT(rel_data(i1,9))
          ELSEIF (relmode.EQ.4) THEN
c...        Over-ride additional cell plasma parameters:
            nvac = nvac + 1
            vacpla(nvac,6) = rel_data(i1,2)
            vacpla(nvac,7) = rel_data(i1,3)
          ELSEIF (relmode.EQ.5) THEN
c...        Over-ride additional cell plasma parameters and bulk plasma
c           parameters:
            nvac = nvac + 1
            vacpla(nvac,6) = rel_data(i1,2)
            vacpla(nvac,7) = rel_data(i1,3)
            osmbulkte = rel_data(i1,4)
            osmbulkti = rel_data(i1,5)
            osmbulkn  = rel_data(i1,6)
            osmbulkv  = rel_data(i1,7)
          ELSEIF (relmode.EQ.6) THEN
c...        Over-ride puff strength:
            nvac = nvac + 1
            eirpuff(nvac,2) = rel_data(i1,2)
          ELSEIF (relmode.EQ.7) THEN
c...        Over-ride bulk plasma values:
            osmbulkte = rel_data(i1,2)
            osmbulkti = rel_data(i1,3)
            osmbulkn  = rel_data(i1,4)
            osmbulkv  = rel_data(i1,5)
          ELSEIF (relmode.EQ.8) THEN
c...        Over-ride exponential decay parameters in InterpolateTargetData:
            DO ir = NINT(rel_data(i1,2)), MIN(nrs,NINT(rel_data(i1,3)))
              IF (idring(ir).EQ.BOUNDARY) CYCLE
              IF (rel_data(i1,4).EQ.1.0.OR.rel_data(i1,4).EQ.3.0)
     .          relexpdat(IKLO,1:3,ir) = rel_data(i1,5:7)
              IF (rel_data(i1,4).EQ.2.0.OR.rel_data(i1,4).EQ.3.0) 
     .          relexpdat(IKHI,1:3,ir) = rel_data(i1,5:7)  
            ENDDO
          ELSE
            CALL ER('SetupRelaxation','Invalid RELMODE value',*99)
          ENDIF
        ENDIF
      ENDDO

c...  Make sure that the additional cell plasma was fully specified:
      IF (((relmode.EQ.4.OR.relmode.EQ.5).AND.nvac.NE.vacnpla ).OR.
     .    ((relmode.EQ.6                ).AND.nvac.NE.eirnpuff)) 
     .  CALL ER('SetupRelaxation','Insufficient additional cell '//
     .          'plasma over-ride data',*99)

      IF (.NOT.status) CALL ER('SetupRelaxation','Cannot find step '//
     .                                           'data',*99)

      RETURN
99    WRITE(0,*) 'REL_STEP=',rel_step
      STOP
      END
c
c =SECTION==============================================================
c
c ======================================================================
c
c  subroutine: MatchProbe
c
      SUBROUTINE MatchProbe(region,ir,mode,status)
      IMPLICIT none

      INTEGER region,ir,mode,status

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INCLUDE 'solparams'
      INCLUDE 'solcommon'
      INCLUDE 'solswitch'

c
c     jdemod - moved variables to fix alignment issue 
c
      COMMON /PEIMULCOM/ osm_peimul2,osmpeires,
     .                   restrictosmpmk,
     .                   restrictmachno
      LOGICAL            restrictosmpmk,restrictmachno
      REAL*8             osm_peimul2(2,MAXNRS),osmpeires(2,MAXNRS)


      INTEGER GetModel     ,
     .        SymmetryPoint,ElectronPowerPoint  ,FlowReversalPoint
      REAL    IdentifyFlowReversal

      REAL       TADJUST      ,LOWFRAC
c      PARAMETER (TADJUST = 1.0,LOWFRAC = 1.0E-13)
      PARAMETER (TADJUST = 1.0,LOWFRAC = 1.0E-15)



      INTEGER   ik,ik1,ik2,i1,i2,flag,initialize,region2,mode_te,result,
     .          eflag,idum1,idum2,idum3,count,count_limit,countdp6

      REAL      norm,norm1,diffe,qem,signe,qem1,signe1,te2,
     .          hte1,hte2,hte3,frac,tmpprbfr,
     .          lgmul(2,MAXNRS),lgte(2,MAXNRS),temark,multi,oldhte1(2),
     .          peimul,hqem(2),dp6adj(MAXNKS)
      REAL*8 deltam,qem3,qem2 (2,MAXNRS),qem0,lgmul2(2,MAXNRS)

      INTEGER ii  (2,MAXNRS),iimark(2,MAXNRS),rstep(MAXNRS),
     .        code(2,MAXNRS),riter (MAXNRS),dmark(2,MAXNRS),
     .        emark(MAXNRS),fmark(MAXNRS),
     .        errorc(2),errorc2(2,MAXNRS),ecount(MAXNRS),flag1(2,MAXNRS)

      REAL    norm2(2,MAXNRS),rate  (2,MAXNRS),temod(MAXNRS),
     .        signe2(2,MAXNRS)


      REAL pmk,pmktot,dpcf,spcf,spcftot,qrat1,srat1,prat1,randv

      REAL oldpmk(MAXNKS,MAXNRS),oldcfe(MAXNKS,MAXNRS)


      INTEGER iks,iksft,isft
      REAL  rdum1,rdum2,rdum3

      REAL    osm_ionmul(2,MAXNRS)

      INTEGER count_faster(2)
      LOGICAL done_revert(2),allow_faster

      CHARACTER*2 tag(2)

      COMMON /OPTTEMP/ osm_matcht,forcet1
      INTEGER          osm_matcht,forcet1

      COMMON /OSMLOCK/ lockionptarg
      LOGICAL          lockionptarg
     


      DATA  initialize,mode_te
     .     /         0,      0/



      COMMON /CFSCOM/ cfs_mage,cfs_magi,cfs_sume,cfs_sumi
      REAL            cfs_mage(MAXNRS),cfs_magi(MAXNRS),
     .                cfs_sume(MAXNRS),cfs_sumi(MAXNRS)

      COMMON /STATCOM/ osm_ncnt ,osm_nerr,osm_temod,osm_relfr,osm_cadj,
     .                 osm_tcon
      INTEGER          osm_ncnt (MAXNRS)    ,osm_nerr (2,MAXNRS),
     .                 osm_cadj (5,2,MAXNRS),osm_tcon (MAXNRS)
      REAL             osm_temod(MAXNRS)    ,osm_relfr(MAXNRS)


      COMMON /OSMOPTS/ osm_fixopt
      INTEGER          osm_fixopt(2,MAXNRS)

      REAL Clock2
      COMMON /TIMESTATS/ timeint
      REAL               timeint(100)

      LOGICAL shapepeak(2)

      CHARACTER note*32

      SAVE

      CALL DB('Matching probe temperature')


c      IF (region.EQ.IKLO)
c     .       WRITE(0,'(A,I6,1P,3E12.4,0P)') 
c     .           'CHECK:',ikbound(ir,IKLO),
c     .           osm_peimul(region,ir),osm_dp6(ikbound(ir,IKLO),ir),
c     .           osmpmk(ikbound(ir,IKLO),ir)



      timeint(5) = Clock2()

      status = 1
      eflag  = 0

      peimul = osm_peimul(region,ir)
c
c     Checks:
c
      IF (idring(ir).EQ.-1)
     .  CALL ER('MatchProbe','Virtual rings are invalid',*99)

c...these are only need to be checked once (???TRUE???) so they could
c   be done elsewhere...
      IF (osm_matcht.GE.1.AND.osm_matcht.LE.3) THEN
        IF (osm_probe.NE.FSP1)
     .    CALL ER('MatchProbe','Invalid reference pressure'//
     .                         ' specification',*99)
        IF (prb_num(osm_probe).LE.0)
     .    CALL ER('MatchProbe','Probe data missing',*99)
        IF (osm_symopt.NE.1)
     .    CALL ER('MatchProbe','Invalid symmetry point option',*99)
      ENDIF
c
c
c
      IF (initialize.EQ.0) THEN
        tag(IKLO) = 'l '
        tag(IKHI) = 'h '

        DO i1 = 1, MAXNRS
          emark(i1) = 0
          fmark(i1) = 0

          DO i2 = 1, 2
            ii     (i2,i1) =   1
            iimark (i2,i1) =   1
            rstep  (   i1) = -99
            riter  (   i1) = -99
            code   (i2,i1) =   0
            rate   (i2,i1) = 1.0
            errorc2(i2,i1) =   0
            dmark  (i2,i1) =   0
            flag1  (i2,i1) =   0

            rel_prbfr(i2,i1) = MAX(0.05,rel_frac)
c            rel_prbfr(i2,i1) = MAX(0.1,rel_frac)
c            rel_prbfr(i2,i1) = MAX(0.01,rel_frac)
          ENDDO
        ENDDO


        initialize = 1
      ENDIF
c
c
c
      IF ((((mode.EQ.0.AND.region.EQ.IKLO).OR.
     .      mode.EQ.2.OR.mode.EQ.3).AND.
     .     (rel_step.NE.rstep(ir).OR.rel_iter.NE.riter(ir))).OR.
c...reset
     .     (rel_prbfr(IKLO,ir).EQ.999.0.AND.
     .      rel_prbfr(IKHI,ir).EQ.999.0)) THEN

        IF (rel_step.NE.rstep(ir)) THEN
          temod(ir) = MAX(0.0,temod(ir)-TADJUST)
          emark(ir) = 0
          fmark(ir) = 0
        ENDIF

        DO i1 = IKLO, IKHI
          ii        (i1,ir) = MIN(2,ii(i1,ir))
          iimark    (i1,ir) = ii(i1,ir)
          rate      (i1,ir) = 1.0
          code      (i1,ir) = 0
          dmark     (i1,ir) = 0
          flag1     (i1,ir) = 0
c          rel_prbfr (i1,ir) = MIN(1.0,MAX(LOWFRAC,rel_prbfr(i1,ir)))
          rel_prbfr (i1,ir) = MIN(0.1,MAX(LOWFRAC,rel_prbfr(i1,ir)))
          lgmul     (i1,ir) = 0.0
          lgmul2    (i1,ir) = 0.0D0
          lgte      (i1,ir) = 0.0
          osm_fixopt(i1,ir) = 0
          hqem      (i1)    = HI
        ENDDO


c...    Equalize relaxation rates for 24-24 rings:
        IF (osm_model(IKLO,ir).EQ.24.AND.osm_model(IKHI,ir).EQ.24) THEN
          rel_prbfr(IKLO,ir) = MIN(rel_prbfr(IKLO,ir),
     .                             rel_prbfr(IKHI,ir))
          rel_prbfr(IKHI,ir) = MIN(rel_prbfr(IKLO,ir),
     .                             rel_prbfr(IKHI,ir))
        ENDIF

c...GOG
        tmpprbfr = rel_prbfr(region,ir)

        ecount(ir) = 0

        rstep(ir) = rel_step
        riter(ir) = rel_iter

        done_revert(IKLO) = .FALSE.
        done_revert(IKHI) = .FALSE.

        count_faster(IKLO) = 0
        count_faster(IKHI) = 0

        count_limit = 0

        oldhte1(IKLO) = -1.0
        oldhte1(IKHI) = -1.0

        restrictosmpmk = .FALSE.
        restrictmachno = .FALSE.

        osmikshift(IKLO,ir) = 0
        osmikshift(IKHI,ir) = 0


        WRITE(0     ,*) 'RESETTING MATCHPROBE PARAMETERS FOR RING ',ir
c        WRITE(0     ,*) '  REL_FRAC = ',rel_frac
c        WRITE(PINOUT,*) 'Resetting relaxation parameters for ',ir

        WRITE(PINOUT,*)
        WRITE(PINOUT,'(3X,2A3,A5,A4,A12,2A8,1X,A3,A10,A6,2X,'//
     .               ' A12,1X,2A4)')
     .    'ir','ec','ii','ik1','norm','t1','t2','s','fr','temod',
     .    'qem1','ec1','ec2'


c...TEMP:
        WRITE(0     ,*) 'PEIMUL:',osm_peimul(IKLO,ir),
     .                            osm_peimul(IKHI,ir)
        WRITE(PINOUT,*) 'PEIMUL:',osm_peimul(IKLO,ir),
     .                            osm_peimul(IKHI,ir)


      ENDIF
c
c     Initialization:
c
      IF (region.EQ.IKLO) THEN
        region2 = IKHI
      ELSE
        region2 = IKLO
      ENDIF

      note = '                '
      norm = -1.0
      flag =  0

      errorc(IKLO) = err1
      errorc(IKHI) = err2



c...GOG
  
      IF (region.EQ.IKLO.AND.ir.EQ.23) THEN
        WRITE(PINOUT,*) 'LAST GO:',osm_peimul(IKLO,ir)
        WRITE(PINOUT,*) 'LAST GO:',cfs_sume(ir)
        WRITE(PINOUT,*) 'LAST GO:',cfs_mage(ir)
        WRITE(PINOUT,*) 'LAST GO:',errorc(region)
        WRITE(PINOUT,*) 'LAST GO:',ikerror(region,ir)
        WRITE(PINOUT,*) 'LAST GO:',ikfluid(region,ir)
      ENDIF

      IF (region.EQ.IKLO.AND.osm_peimul(IKLO,ir).LT.0.1.AND.
     .    .NOT.firstshape) THEN
      
C... DO THIS IN CONVERGE SOLUTION?


        countdp6 = 0
        DO WHILE((cfs_sume(ir).LE.0.0.OR.
     .            cfs_mage(ir).LE.0.0.OR.
     .            (errorc(region).GT.1.AND.
     .             ikerror(region,ir).LT.ikfluid(region,ir))).AND.
     .           countdp6.LT.3) 
      
c...      Knock out the multiplier:
      
          WRITE(0     ,*) '*********************************'
          WRITE(0     ,*) '* CRUSHING OSM_DP6 IN MATCHPROBE*'
          WRITE(0     ,*) '*********************************'
      
          WRITE(PINOUT,*) '**********************************'
          WRITE(PINOUT,*) '* CRUSHING OSM_DP6 IN MATCHPROBE *'
          WRITE(PINOUT,*) '**********************************'
      
          DO ik = 1, osm_sympt(ir)
            IF (osm_dp6(ik,ir).NE.1.0) 
     .        osm_dp6(ik,ir) = osm_dp6(ik,ir) * 0.10
          ENDDO
      

          CALL INITPLASMA(ir,ir,3)
          CALL SOL_PLASMA(ir,ir,3)
          CALL SOL       (ir,ir,3)

          countdp6 = countdp6 + 1
        ENDDO
      
        IF (countdp6.GE.3) THEN

          CALL SaveSolution

          osm_mode = 2
  
          CALL INITPLASMA(ir,ir,3)
          CALL SOL_PLASMA(ir,ir,3)
          CALL SOL       (ir,ir,3)

          CALL ER('MatchProbe','Trouble',*99)

        ENDIF
        

      ENDIF




c...  This stops the expert from falling down an abyss of reverts
c     if the error is not resolved after the first revert:
      IF (done_revert(region).AND.
     .    errorc(IKLO).LE.1.AND.errorc(IKHI).LE.1) 
     .  done_revert(region) = .FALSE.

c      WRITE(0     ,*) 'DONE_REVERT:',done_revert(region),
c     .                            count_faster(region),code(region,ir)
c      WRITE(PINOUT,*) 'DONE_REVERT:',done_revert(region),
c     .                            count_faster(region),code(region,ir)


c...  This isn't going to work well when a solution is started from a
c     previous iteration, since the restriction will all of a sudden
c     be gone, and so the solution will be way off:
c      IF ((stopopt3.EQ.19.OR.stopopt3.EQ.20).AND.
c     .    ir.GE.27.AND.ir.LE.30.AND.
c     .    osm_model(IKLO,ir).EQ.24.AND.osm_model(IKHI,ir).EQ.22.AND.
c     .    rel_prbfr(IKLO,ir).LT.1.0E-10.AND.
c     .    errorc(region).LE.1.AND.errorc(region2).LE.1) 
c     .  restrictosmpmk = .TRUE.

c      restrictosmpmk = .FALSE.

c      WRITE(0     ,*) 'RESTRICT OSMPMK:',restrictosmpmk,restrictmachno,
c     .                osmikshift(IKLO,ir),osmikshift(IKHI,ir)
c      WRITE(PINOUT,*) 'RESTRICT OSMPMK:',restrictosmpmk,restrictmachno,
c     .                osmikshift(IKLO,ir),osmikshift(IKHI,ir)

      ik1 = osm_sympt(ir)
      ik2 = ik1 + 1
c
c     Store last good fraction:
c
c
      IF (errorc(IKLO).LE.0.AND.errorc(IKHI).LE.0) THEN
        IF (restrictosmpmk) THEN
          lgmul2(region,ir) = osmpeires  (region,ir)
        ELSE
          lgmul (region,ir) = osm_peimul (region,ir)
          lgmul2(region,ir) = osm_peimul2(region,ir)
        ENDIF

        IF     (region.EQ.IKLO) THEN
          lgte(region,ir) = ktebs(ik1,ir)
        ELSEIF (region.EQ.IKHI) THEN
          lgte(region,ir) = ktebs(ik2,ir)
        ENDIF

        ecount(ir) = 0
c        ecount(ir) = MAX(ecount(ir)-1,0)

c      ELSEIF (errorc(region).GE.1.AND.errorc(region2,ir).GE.1) THEN
c        ecount(ir) = MAX(ecount(ir)+0,0)
c        ecount(ir) = MAX(ecount(ir)+2,0)

c      ELSEIF (errorc(region).EQ.1) THEN
c        ecount(ir) = MAX(ecount(ir)+1,0)

      ELSEIF (errorc(region).GT.1.AND.
     .        .NOT.(osm_code  (region,ir).EQ.2.AND.
     .              osm_fixopt(region,ir).EQ.0)) THEN
        ecount(ir) = MAX(ecount(ir)+2,0)

c      ELSEIF (region.EQ.IKHI.AND.
c     .        errorc(region).GT.1.AND.errorc2(region,ir).GT.1) THEN
c        ecount(ir) = MAX(ecount(ir)+2,0)
      ENDIF

c
c     Set upstream temperature to be matched:
c
      IF     (mode.EQ.0) THEN
        temark = prp_te(ir,osm_probe) * rel_mfsp(region,ir)

c     ELSEIF (mode.EQ.1) THEN
c...    Match to probe Ti.

      ELSEIF (mode.EQ.2) THEN
c...    Match inner and outer symmetry point Te values:
        IF (region.EQ.IKLO) THEN
c...      temp! until I can get the Te matching error condtion to work without
c         adjusting the mock power term for both half-rings...
          IF     (osm_matcht.EQ.5) THEN
            temark = ktebs(ik1,ir) * rel_mfsp(region,ir)
          ELSEIF (osm_matcht.EQ.4.AND.osm_model(IKHI,ir).EQ.24) THEN
            IF (tuprat.GT.0.0) THEN
              temark = tuprat 
            ELSE
              WRITE(0,*) 'STOP: NO INNER UPSTREAM TEMPERATURE FOUND'
              STOP
            ENDIF         
          ELSE
            temark = ktebs(ik2,ir) * rel_mfsp(region,ir)
          ENDIF
        ELSE 
          IF (osm_matcht.EQ.4) THEN
c            temark = ktebs(ik2+1,ir) + (ktebs(ik2,ir)-ktebs(ik2+1,ir))*
c     .               (kss(ik2-1,ir) - kss(ik2+1,ir)) / 
c     .               (kss(ik2  ,ir) - kss(ik2+1,ir))

c...LOGIC HERE IS SCREWED...REALLY NEED TO BE USING OSM_MATCHT=5,6,etc.

            IF     ((region.EQ.IKLO.AND.osm_model(IKHI,ir).EQ.22).OR.
     .              (region.EQ.IKHI.AND.osm_model(IKHI,ir).EQ.22)) THEN
 
              temark = ktebs(ik2,ir) * rel_mfsp(region,ir)

            ELSEIF ((region.EQ.IKHI.AND.osm_model(IKLO,ir).EQ.22).OR.
     .              (region.EQ.IKLO.AND.osm_model(IKLO,ir).EQ.22)) THEN
              temark = ktebs(ik1,ir) * rel_mfsp(region,ir)

            ELSE
c...tempA
              IF (tuprati.GT.0.0) THEN
                temark = tuprati
              ELSE
                WRITE(0,*) 'STOP: NO OUTER UPSTREAM TEMPERATURE FOUND'
                STOP
              ENDIF         
            ENDIF
          ELSE
            temark = ktebs(ik1,ir) * rel_mfsp(region,ir)
          ENDIF
        ENDIF

c      ELSEIF (mode.EQ.3) THEN
c...    Match inner and outer symmetry point Ti values:

      ELSE
        CALL ER('MatchProbe','Invalid MODE value',*99)
      ENDIF


      IF     (region.EQ.IKLO) THEN
c
c      jdemod - changed osmikshift(IKLO) to osmikshift(IKLO,ir)
c               in hte1 assignment
c
       hte1 = ktebs(ik1-osmikshift(IKLO,ir),ir)
c        hte1 = ktebs(ik1,ir)
        hte2 = temark
        hte3 = ktebs(ik2,ir)

        IF (errorc(region).GT.1) hte1 = 0.0

        diffe = (hte2 - hte1) / hte2

      ELSEIF (region.EQ.IKHI) THEN
        hte3 = ktebs(ik1,ir)
        hte2 = temark
        hte1 = ktebs(ik2,ir)

        IF (errorc(region).GT.1) hte1 = 0.0

        diffe = (hte2 - hte1) / hte2

      ENDIF


c...  Need a condition to prevent "FASTER" from occuring when 
c     TEMARK is a low temperature, and the algorithm is proceeding at a
c     reasonable pace.  Unfortunately, things are so sensitive at these
c     low temperatures that "FASTER" can cause the algorithm to skip
c     around the solution:
      allow_faster = .TRUE.
      IF (errorc(region).LE.1.AND.errorc(region2).LE.1.AND.
     .    temark.LT.3.0.AND.ABS(hte1-oldhte1(region)).GT.0.10.AND.
     .    hte1.LT.2*temark) THEN  
        allow_faster = .FALSE.
      ENDIF


c      WRITE(0 ,*) 'ALLOW_FASTER:',allow_faster,hte1,oldhte1(region)
c      WRITE(88,*) 'ALLOW_FASTER:',allow_faster,hte1,oldhte1(region)

      oldhte1(region) = hte1

 



      IF     (osm_powopt.EQ.1) THEN
        qem  =      rel_qemul(region,ir)
        qem0 = DBLE(rel_qemul(region,ir))
      ELSEIF (osm_powopt.EQ.2) THEN
        IF (restrictosmpmk) THEN
          qem  = SNGL(osmpeires(region,ir))
          qem0 =      osmpeires(region,ir)
        ELSE
          qem  = osm_peimul (region,ir)
          qem0 = osm_peimul2(region,ir)
        ENDIF
      ELSE
        CALL ER('MatchProbe','Invalid OSM_POWOPT value',*99)
      ENDIF


      signe = SIGN(1.0,diffe)
      diffe = ABS(diffe)
      norm  = diffe
c
c
c
c
c
c      WRITE(0,*) '  [',hte1,hte2,hte3,code(region2,ir),']'



c...  Scan temperatures in the peak shaping region to see if they are
c     becoming depressed:

c...peak
c      WRITE(0,*) 'FIRSTSHAPE=',firstshape

      shapepeak(IKLO) = .FALSE.
      shapepeak(IKHI) = .FALSE.

      IF ((osm_model(IKLO,ir).EQ.24.OR.osm_model(IKHI,ir).EQ.24).AND.
     .    osm_matcht.EQ.4.AND..NOT.firstshape.AND.
     .    errorc(IKLO).LE.1.AND.errorc(IKHI).LE.1.AND.
     .    region.EQ.IKLO) THEN

        IF (osmmock.NE.0) iksft = -1
        DO ik = ikbound(ir,IKLO)+1, ikfluid(IKLO,ir)-1
          IF (ktebs(ik,ir).LT.0.75*ktebs(ikbound(ir,IKLO),ir).AND.
     .        osm_dp6(ik+iksft,ir).NE.1.0) shapepeak(IKLO) = .TRUE.
        ENDDO

        iksft = 0
        DO ik = ikbound(ir,IKHI)-1, ikfluid(IKHI,ir)+1, -1
c          WRITE(0,*) '--->',ik,ir,ktebs(ik,ir),
c     .                      ktebs(ikbound(ir,IKLO),ir),
c     .                      osm_dp6(ik+iksft,ir)
          IF (ktebs(ik,ir).LT.0.75*ktebs(ikbound(ir,IKHI),ir).AND.
     .        osm_dp6(ik+iksft,ir).NE.1.0) THEN
c            WRITE(0,*) 'CHECKING: ',ik,ir,iksft
            shapepeak(IKHI) = .TRUE.
          ENDIF
        ENDDO
      ENDIF

c... 
      IF (.NOT.firstshape.AND.region.EQ.IKLO.AND.
     .    errorc(region).GT.1.AND.ikerror(region,ir).LT.
     .    ikfluid(region,ir)) THEN

        shapepeak(region) = .TRUE.

      ENDIF

      IF   (shapepeak(IKLO).OR.shapepeak(IKHI)) THEN     
c        CALL SaveSolution
c        STOP 'JSGsDFSD'

        WRITE(0     ,*) 'EVALUATE TERMS 0'
        WRITE(PINOUT,*) 'EVALUATE TERMS 0'

        count_limit = 0

        restrictosmpmk = .FALSE.
        restrictmachno = .FALSE.

c...    Do not let peak shaping region expand:
        idum1 = -1
        idum2 = -1

        IF (shapepeak(IKLO)) CALL EvaluateDensityPeak(IKLO,ir,idum1)
        IF (shapepeak(IKHI)) CALL EvaluateDensityPeak(IKHI,ir,idum2)
	
        note            = 'SHAPING PEAK'
        code(region,ir) =  0
        flag            =  0
        status          =  2

        IF (idum1.EQ.-1) THEN

c...      Fail:
          status = -3

        ENDIF

      ELSEIF ((ecount(ir).GT.200.AND.region.EQ.IKLO).OR.
     .        (code(region2,ir).EQ.12).OR.code(region2,ir).EQ.16) THEN

        note            = 'EXHAUSTED'
        code(region,ir) = 12
        status          = -3
c        status          = 1
        flag            = 0

        norm = 0.0

c        IF (region.EQ.IKHI) eflag = 1

      ELSEIF (code(region2,ir).EQ.18) THEN

        note            = 'BEST'
        code(region,ir) = 18
        flag            =  0
        status          =  0

      ELSEIF ((norm.LT.0.01*rel_tol.OR.
c..this condition ABS(norm).LT.0.01 is relative to the osm_peimul normalization!
c..I had forgotten about this other path to satisfying the upstream temperature match
c  requirement -- I assume it is here because the algorithm starts to choke if QEM is 
c  small:
     .        (norm.LT.0.10*rel_tol.AND.ABS(qem).LT.0.01).OR.
c...For small symmetry point temperature, be merciful:
     .        (ABS(hte1-hte2).LT.0.10.AND.temark.LT.3.0).OR.
     .        (ABS(hte1-hte2).LT.0.25.AND.temark.GE.3.0)).AND.
c     .        (ABS(hte1-hte2).LT.0.25)).AND.
     .     code(region2,ir).NE.14) THEN

c     .    ((ABS(ktebs(ik1,ir)-ktebs(ik2,ir)).LT.rel_tol.AND.
c     .          ktebs(ik1,ir)-ktebs(ik2,ir) .NE.0.0).AND.
c     .     ii(region,ir).EQ.1.AND.rel_count.EQ.1.AND.
c     .     osm_store.GT.-1)

        IF (fmark(ir).EQ.1.AND.region.EQ.IKHI.AND.
     .      code(IKLO,ir).EQ.6) THEN

          note            = 'DONE [GO]'

          code(region,ir) = 6
          status          = 1
          flag            = 0

          norm = 0.0

          eflag = 1
        ELSE
          result = 0

          IF (region.EQ.IKHI.AND.code(region2,ir).EQ.6) THEN

c            IF (.NOT.(code(IKHI,ir).EQ.21.AND.
c     .                osm_peimul(IKLO,ir).EQ.hqem(IKLO).AND.
c     .                osm_peimul(IKHI,ir).EQ.hqem(IKHI))) THEN

              IF (osm_mode.GT.0) THEN
                WRITE(0     ,*) 'EVALUATE TERMS 1'
                WRITE(PINOUT,*) 'EVALUATE TERMS 1'
              ENDIF

              count_limit = 0

              restrictosmpmk = .FALSE.
              restrictmachno = .FALSE.

              osm_fixopt(IKLO,ir) = 2
              osm_fixopt(IKHI,ir) = 2

              idum1 = 0
              idum2 = 0

              IF (GetModel(region2,ir).EQ.24) 
     .          CALL EvaluatePowerTerms(IKLO,ir,errorc(IKLO),idum1)
              IF (GetModel(region ,ir).EQ.24) 
     .          CALL EvaluatePowerTerms(IKHI,ir,errorc(IKHI),idum2)

              IF (GetModel(region2,ir).EQ.24) 
     .          CALL EvaluateDensityPeak(IKLO,ir,idum1)
              IF (GetModel(region ,ir).EQ.24) 
     .          CALL EvaluateDensityPeak(IKHI,ir,idum2)

              firstshape = .FALSE.

              IF (osm_mode.GT.0) WRITE(0,*) 'IDUM1,2=',idum1,idum2

              IF     (idum1.EQ.-1.OR.idum2.EQ.-1) THEN
                result = -1
              ELSEIF (idum1.EQ.-2.OR.idum2.EQ.-2) THEN
                result = -2
              ELSEIF (idum1.EQ.1.OR.idum2.EQ.1) THEN

                IF (osm_peimul(IKLO,ir).EQ.0.0.AND.idum1.EQ.1) THEN
                  WRITE(0     ,*) 'THIS MAKES FLATTENING PROFILE HARD!'
                  WRITE(PINOUT,*) 'THIS MAKES FLATTENING PROFILE HARD!'
                  osm_peimul(IKLO,ir) = 0.1
                ENDIF

                IF (osm_peimul(IKHI,ir).EQ.0.0.AND.idum2.EQ.1) THEN
                  WRITE(0     ,*) 'THIS MAKES FLATTENING PROFILE HARD!'
                  WRITE(PINOUT,*) 'THIS MAKES FLATTENING PROFILE HARD!'
                  osm_peimul(IKHI,ir) = 0.1
                ENDIF


                result =  1
              ENDIF

c            ELSE
c              WRITE(0     ,*) 'SOLUTION ALGORITHM FAILURE DETECTED'
c              WRITE(PINOUT,*) 'SOLUTION ALGORITHM FAILURE DETECTED'
c            ENDIF

c            cnt = cnt + 1
c            IF (cnt.LT.15) THEN
c              WRITE(0,*) '*** FORCING RESULT=1 ***'
c              result = 1
c            ENDIF

c              WRITE(0,*) '*** FORCING RESULT=0 ***'
c              result = 0

          ENDIF

          IF     (result.EQ.-1) THEN
c
c...        Reduce particle source:
            note            = 'CHECK SOURCE'
            code(region,ir) = 19
            status          = 0
            flag            = 2
            eflag           = 1
            emark(ir)       = 0
          ELSEIF (result.EQ.0) THEN
            note             = 'DONE'
            code (region,ir) = 6
            rate (region,ir) = 1.0
            flag             = 2
            status           = 0
            dmark(region,ir) = 1
          ELSE
            IF (.NOT.(code(IKHI,ir).EQ.21.AND.
     .                osm_peimul(IKLO,ir).EQ.hqem(IKLO).AND.
     .                osm_peimul(IKHI,ir).EQ.hqem(IKHI).AND.
     .                osmmock.EQ.0)) THEN
              note             = 'DONE [C]'
              code (region,ir) = 21
              hqem(IKLO) = osm_peimul(IKLO,ir)
              hqem(IKHI) = osm_peimul(IKHI,ir)
              flag             = 0
              status           = 3
c              code (region,ir) = 0
c              rel_prbfr(region ,ir) = rel_prbfr(region ,ir) * 0.1
c              rel_prbfr(region2,ir) = rel_prbfr(region2,ir) * 0.1
            ELSE
              WRITE(0     ,*) 'SOLUTION ALGORITHM FAILURE DETECTED'
              WRITE(PINOUT,*) 'SOLUTION ALGORITHM FAILURE DETECTED'
              note             = 'DONE [FAILURE]'
              code (region,ir) = 6
              rate (region,ir) = 1.0
              flag             = 2
              status           = 0
              dmark(region,ir) = 1
            ENDIF
          ENDIF
        ENDIF

      ELSEIF (1.EQ.2.AND.rel_prbfr(region,ir).EQ.LOWFRAC.AND.
     .        errorc(region).GT.0.AND.
     .        code(region2,ir).EQ.6) THEN

        result = 0

        STOP 'WOW, REALLY GET HERE'

        WRITE(0     ,*) 'EVALUATE TERMS 2'
        WRITE(PINOUT,*) 'EVALUATE TERMS 2'

        IF (switch(SWPOW3).EQ.4.0.OR.switch(SWPOW3).EQ.5.0.OR.
     .      switch(SWPOW3).EQ.1.0.OR.switch(SWPOW3).EQ.6.0.OR.
     .      switch(SWPOW3).EQ.7.0) THEN
          CALL EvaluatePowerTerms (region,ir,errorc(region),result)
          CALL EvaluateDensityPeak(region,ir,result)
        ENDIF

        IF (result.EQ.1) THEN

          note            = 'SPREAD'
          rel_prbfr(region,ir) = 0.1
          rel_prbfr(region,ir) = rel_prbfr(region,ir) * 0.2
          code(region,ir) = 15
          flag            = 0

        ELSEIF (result.EQ.-1) THEN

            WRITE(0,*) 'MARK: LGMUL1'

            IF (restrictosmpmk) THEN
              STOP 'CODE DEVELOPMENT REQUIRED A121'
            ELSE
              osm_peimul(IKLO,ir) = lgmul(IKLO,ir)
              osm_peimul(IKHI,ir) = lgmul(IKHI,ir)
            ENDIF

c            DO ik = 1, nks(ir)
c              osm_dp4(ik,ir) = 1.0
c            ENDDO

            DO ik = 1, nks(ir)
              IF (osm_dp4(ik,ir).LT.0.0) THEN
                osm_dp4(ik,ir) = 0.90 * osm_dp4(ik,ir)
              ENDIF
            ENDDO

            WRITE(0     ,*) 'EVALUATE TERMS X'
            WRITE(PINOUT,*) 'EVALUATE TERMS X'

            CALL EvaluatePowerTerms(IKLO,ir,errorc(IKLO),idum1)
            CALL EvaluatePowerTerms(IKHI,ir,errorc(IKHI),idum2)

            CALL EvaluateDensityPeak(IKLO,ir,idum1)
            CALL EvaluateDensityPeak(IKHI,ir,idum2)


            IF (idum1.EQ.-1.OR.idum2.EQ.-1) THEN
              CALL ER('MatchProbe','CRAP!',*99)
            ELSE
c...not sure this will do the trick, especially if there hasn't been a good
c solution, and if rel_prbfr needs to be adjusted before switching to
c symmetry point matching...

              note            = 'BEST 2'
              code(region,ir) = 18
              flag            =  0
              status          =  0
           ENDIF
        ELSE

          note            = 'CRAP'

          code(region,ir) = 14
          flag            = 0
          status          = 0

          norm = 0.0

          eflag = 1

        ENDIF


      ELSEIF (rel_prbfr(region ,ir).LT.LOWFRAC.OR.
     .        rel_prbfr(region2,ir).LT.LOWFRAC.OR.
     .        code(region2,ir).EQ.14) THEN

        note = 'LOWER LIMIT'


        IF (count_limit.EQ.5) THEN
c        IF (count_limit.EQ.1) THEN
c        IF (.TRUE..OR.count_limit.EQ.1) THEN
          IF (region.EQ.IKLO) THEN
            code(region,ir) = 14
          ELSE
            code(region ,ir) = 0
            code(region2,ir) = 0
          ENDIF
 
c...      GIVE UP ON THIS RING
          status = -3
          flag   = 0

          norm = 0.0
        ELSEIF (region.EQ.IKLO) THEN
          note = note(1:LEN_TRIM(note))//' [JOSTLE]'
          count_limit = count_limit + 1

c...      Adjust the REL_PRB scaling by some random amount:

c          CALL SURAND(relseed,1,randv)
c  
c          randv = randv + 0.01
c
c          rel_prbfr(IKLO,ir) = rel_prbfr(IKLO,ir) * randv * 1000.0  
c          rel_prbfr(IKHI,ir) = rel_prbfr(IKLO,ir) * randv * 1000.0 
c
c          WRITE(0     ,*) 'RANDV:',randv,count_limit
c          WRITE(PINOUT,*) 'RANDV:',randv,count_limit

          rel_prbfr(IKLO,ir) = rel_prbfr(IKLO,ir) * 100.0  
          rel_prbfr(IKHI,ir) = rel_prbfr(IKLO,ir) * 100.0 

c          restrictosmpmk = .TRUE.

c          osmpeires(IKLO,ir) = osm_peimul2(IKLO,ir)
c          osmpeires(IKHI,ir) = osm_peimul2(IKHI,ir)


c          qem  = 1.0
c          qem0 = 1.0D0
c          qem2(IKLO,ir) = 1.0D0
c          qem2(IKHI,ir) = 1.0D0

c          lgmul (IKLO,ir) = 1.0
c          lgmul (IKHI,ir) = 1.0
c          lgmul2(IKLO,ir) = 1.0D0
c          lgmul2(IKLO,ir) = 1.0D0

          norm1 = 0.0

          IF (osm_model(IKLO,ir).EQ.24.AND.
     .        osm_model(IKHI,ir).EQ.22) THEN
            osmikshift(region,ir) = osmikshift(region,ir) + 1
          ELSE
            count_limit = 5          
          ENDIF

        ENDIF


      ELSEIF (ii(region,ir).GT.1) THEN

        IF (code(region2,ir).EQ.1) THEN
c        IF (code(region2,ir).EQ.1.OR.errorc(region2).GT.1) THEN

            note            = 'WAITING'

            code(region,ir) = 9

            flag   = 1

        ELSEIF (1.EQ.1.AND.
c     .          norm.LT.norm2(region,ir).AND.
     .          signe.NE.signe2(region,ir).AND.
     .          ii(region,ir).NE.iimark(region,ir).AND.
     .          errorc(IKLO).LE.0.AND.errorc(IKHI).LE.0.AND.
     .          errorc2(IKLO,ir).LE.0.AND.errorc2(IKHI,ir).LE.0.AND.
     .          (code(region ,ir).EQ.0.OR.code(region ,ir).EQ. 6.OR.
     .           code(region ,ir).EQ.7.OR.code(region ,ir).EQ.21).AND.
     .          (code(region2,ir).EQ.0.OR.code(region2,ir).EQ. 6.OR.
     .           code(region2,ir).EQ.7.OR.code(region2,ir).EQ.21)
     .          ) THEN

          rel_prbfr(region,ir) = rel_prbfr(region,ir) * 0.2
          note            = 'FLIP'

c          IF (errorc(region).LE.0.AND.errorc(region2).LE.0)
c     .      rel_prbfr(region,ir) = MAX(rel_prbfr(region,ir),LOWFRAC)

          code(region,ir) = 1
          flag            = 1

        ELSEIF (norm.LT.rel_tol*10.0*LOWFRAC.AND.
     .          ii(region,ir).EQ.iimark(region,ir)) THEN

          IF (norm .LT.norm2 (region,ir).AND.
     .        signe.EQ.signe2(region,ir)) THEN
            note            = 'CONVERGE'
            code(region,ir) = 2
            rate(region,ir) = 1.0
            flag          = 2
            status        = 0
          ELSE
            note            = 'SCORE'
            code(region,ir) = 3
            rate(region,ir) = 1.0
            flag          = 2
            status       = 0
          ENDIF

c        ELSEIF (
c     .          (ii(region,ir).NE.iimark(region,ir).AND.
c     .           ((errorc (region    ).GT.1.AND.
c     .             errorc2(region ,ir).LE.1.AND.
c     .             errorc2(region2,ir).LE.1.AND.
c     .        .NOT.(signe2(region ,ir).GT.0.0.AND.
c     .              signe2(region2,ir).GT.0.0)).OR.
c     .            (errorc (region    ).EQ.1.AND.signe.EQ.+1)).AND.
c     .          (norm2(region2,ir).LT.10.0*rel_tol)
c     .           .AND.stopopt2.NE.7)
c...Removed May 25, 2001.  The solver behaves a little better now, so I am
c   relaxing these conditions a little bit.  If a sonic transition occurs, then
c   there is not much to do about it in most cases, since the sources 
c   are highly controlled at the moment.  Likely it is a problem with bounary conditions, or
c   how SOL24 is being setup. -SL
cc     .      .OR.(stopopt2.NE.7.AND.
cc     .           errorc (region    ).GT.0.AND.
cc     .           errorc (region2   ).GT.0.AND.
cc     .           errorc2(region,ir ).LE.0)
c     .      .OR.(stopopt2.EQ.7.AND.
c     .           errorc (region    ).GT.1.AND.
c     .           errorc (region2   ).GT.1.AND.
c     .           errorc2(region,ir ).LE.0)
c     .         ) THEN
c
c          rel_prbfr(region,ir) = rel_prbfr(region,ir) * 0.2
c
cc          IF ((errorc(region2).EQ.0.OR.errorc(region).LE.1).AND.
cc     .        code(region2,ir).NE.6)
cc     .      rel_prbfr(region,ir) = MAX(LOWFRAC,rel_prbfr(region,ir))
c
c          IF     (osm_fixopt          (region,ir).GE.  1.AND.
c     .            errorc              (region)   .GT.  0.AND.
c     .            GetModel            (region,ir).EQ. 24
cc     .            IdentifyFlowReversal(region,ir).NE.0.0.AND.
cc     .            rel_prbfr           (region,ir).LT.100.0*LOWFRAC
c     .            ) THEN
c
c            osm_fixopt(region,ir) = 2
c
c            WRITE(0     ,*) 'EVALUATE TERMS Y'
c            WRITE(PINOUT,*) 'EVALUATE TERMS Y'
c
c            CALL EvaluatePowerTerms (region,ir,errorc(region),idum1)
c            CALL EvaluateDensityPeak(region,ir,idum1)
c
c          ELSEIF (osm_fixopt          (region,ir).EQ.0 .AND.
c     .            errorc              (region)   .GT.0 .AND.
c     .            GetModel            (region,ir).EQ.24.AND.
c     .            IdentifyFlowReversal(region,ir).NE.0.0) THEN
c            osm_fixopt(region,ir) = 1
c          ENDIF
c
c          WRITE(0,*) 'FLOW REVERSAL ',
c     .      region,osm_fixopt(region,ir),
c     .      errorc(region),
c     .      GetModel(region,ir),IdentifyFlowReversal(region,ir)
c
c          note            = 'SUPPRESS'
c          code(region,ir) = 10
c          flag            = 0
c          status          = 1
c
        ELSEIF (norm.GT.2.0*norm2(region,ir).AND.
c...weak indication that the algorithm has been reset
     .          rel_prbfr(region,ir).NE.0.1.AND.
     .          ii(region,ir).NE.iimark(region,ir).AND.
     .          (code(region ,ir).EQ.0.OR.code(region ,ir).EQ. 6.OR.
     .           code(region ,ir).EQ.7.OR.code(region ,ir).EQ.21).AND.
     .          (code(region2,ir).EQ.0.OR.code(region2,ir).EQ. 6.OR.
     .           code(region2,ir).EQ.7.OR.code(region2,ir).EQ.21).AND.
     .          ((stopopt2.NE.7.AND.
     .            errorc(IKLO).LE.0.AND.errorc(IKHI).LE.0.AND.
     .            errorc2(IKLO,ir).LE.0.AND.errorc2(IKHI,ir).LE.0).OR.
     .           (stopopt2.EQ.7.AND.
     .            errorc(IKLO).LE.1.AND.errorc(IKHI).LE.1.AND.
c     .            errorc2(IKLO,ir).LE.1.AND.errorc2(IKHI,ir).LE.1))
		
     .            errorc2(IKLO,ir).LE.1.AND.errorc2(IKHI,ir).LE.1)).OR.

c...peak


c     .           (osmmock.NE.0.AND.region.EQ.IKLO.AND.
c     .            errorc(region).GT.1.AND.
c     .            ikerror(IKLO,ir).LT.ikfluid(IKLO,ir).AND.
c     .            osm_dp6(ikerror(IKLO,ir)-1,ir).NE.1.0.AND.
c     .            firstshape)
c...WHAT TO DO?
cc     .            .NOT.firstshape)

c...I put this in because a ring was getting stuck in an infinite loop when
c   when the solution was close enough to matching that 'FASTER' was not triggered, but
c   it kept overshooting and striking an error upstream.  The multipler would then recover,
c   but be right back where it started from.  Not sure the following will be a help
c   in all cases, but it works here:

     .             (errorc(region).GT.1.AND.errorc2(region,ir).LE.1.AND.
     .              code(region2,ir).EQ.6.AND.
     .              .NOT.done_revert(region)).OR.

c...This line addresses another rare problem, like the one above, so hopefully it 
c   won't cause an error more generally:
     .             (errorc(region).LE.1.AND.errorc2(region,ir).GT.1.AND.
     .              code(region,ir).EQ.0)

     .          ) THEN


c     .          errorc(IKLO).LE.0.AND.errorc(IKHI).LE.0.AND.
c     .          errorc2(IKLO,ir).LE.0.AND.errorc2(IKHI,ir).LE.0) THEN


          rel_prbfr(region,ir) = rel_prbfr(region,ir) * 0.2


          IF (errorc(region).LE.0.AND.errorc(region2).LE.0)
     .      rel_prbfr(region,ir) = MAX(rel_prbfr(region,ir),LOWFRAC)

          done_revert(region) = .TRUE.



          note            = 'REVERT'
          code(region,ir) = 5
          flag            = 1

c          IF (errorc2(IKLO,ir).LE.1.AND.errorc2(IKHI,ir))
c     .      flag = 3
          IF (errorc2(IKLO,ir).LE.1.AND.errorc2(IKHI,ir).LE.1.AND.
     .        lgmul(region,ir).NE.0.0) 
c... BUG: Mar 14, 2003
c     .        lgmul(region).NE.0.0) 
     .      flag = 3



        ELSEIF (allow_faster.AND.(
c     .          norm       .NE.norm2 (region,ir).AND.
     .          signe                 .EQ.signe2(region,ir).AND.
     .          (((norm2(region,ir)-norm).LT.0.5*norm.AND.
     .             (GetModel(region ,ir).EQ.24.AND.
     .              GetModel(region2,ir).EQ.24)).OR.
     .           ((norm2(region,ir)-norm).LT.0.5*norm.AND.
     .             (GetModel(region ,ir).EQ.22.OR.
     .              GetModel(region2,ir).EQ.22))).AND.
     .          ii(region,ir).NE.iimark(region,ir).AND.
     .          (code(region,ir).EQ. 0.OR.code(region,ir).EQ. 7.OR.
     .           code(region,ir).EQ.15.OR.code(region,ir).EQ.16.OR.
     .           code(region,ir).EQ.21)
     .          .AND.
     .          (errorc(region).GE.0.OR.errorc(region2).LE.1).AND.
     .          .NOT.(errorc    (region   ).GT.1.AND.
     .                osm_code  (region,ir).EQ.2.AND.
     .                osm_fixopt(region,ir).EQ.0))
     .          ) THEN


c...  Reset counter for how many 'FASTER' events there has been in a row:
          IF (code(region,ir).EQ.7) THEN
            count_faster(region) = count_faster(region) + 1
          ELSE
            count_faster(region) = 0
          ENDIF

          note            = 'FASTER'
          code(region,ir) = 7

          IF     (count_faster(region).GE.9.AND.
     .            DABS(osm_peimul2(region,ir)).GT.0.1D0.AND.
     .            norm.GT.0.05*rel_tol) THEN
            rel_prbfr(region,ir) = rel_prbfr(region,ir) *
     .                             (10.0 * rate(region,ir))
            count_faster(region) = 0
            note = note(1:LEN_TRIM(note))//' [JOLT]'
          ELSEIF (errorc(region).LE.1.AND.errorc(region2).LE.1.AND.
     .        hte1.LT.2.0*temark.AND.temark.LT.3.0) THEN
            rel_prbfr(region,ir) = rel_prbfr(region,ir) *
     .                             (1.02 * rate(region,ir))
            note = note(1:LEN_TRIM(note))//' [SUBDUED]'
          ELSEIF (errorc(region).LE.1) THEN
c...        Increase rate of OSM_PEIMUL adjustment:
            rel_prbfr(region,ir) = rel_prbfr(region,ir) *
     .                             (1.5 * rate(region,ir))
          ELSE
            rel_prbfr(region,ir) = rel_prbfr(region,ir) *
     .                             (1.1 * rate(region,ir))
          ENDIF

        ELSE
          rate(region,ir) = 1.0
          code(region,ir) = 0
        ENDIF

      ENDIF

      CALL DB('DONE ADJUSTING PARAMETERS')
c
c
c
c
c
      IF (((errorc(region2).LE.1.AND.errorc(region).GT.0.AND.
     .      norm2(region2,ir).LT.100.0*LOWFRAC*rel_tol).OR.
     .     (errorc(region2).LE.1.AND.errorc(region).GT.0.AND.
     .      rel_prbfr(region,ir).LT.0.01)).AND.
     .    ii(region,ir).NE.iimark(region,ir).AND.
     .    code(region,ir).NE.9.AND.osm_code(region,ir).EQ.2
     .   ) THEN

        WRITE(0     ,*) 'EVALUATE TERMS 3'
        WRITE(PINOUT,*) 'EVALUATE TERMS 3'

        result = 0

        IF (switch(SWPOW3).EQ.4.0.OR.switch(SWPOW3).EQ.5.0.OR.
     .      switch(SWPOW3).EQ.1.0.OR.switch(SWPOW3).EQ.6.0.OR.
     .      switch(SWPOW3).EQ.7.0) THEN
          CALL EvaluatePowerTerms (region,ir,errorc(region),result)
          CALL EvaluateDensityPeak(region,ir,result)
        ENDIF

        IF (result.EQ.-1) THEN
c
c           Reduce particle source:
c
            note            = 'CHECK SOURCE'
            code(region,ir) = 19
            status          = 0
            flag            = 2

            eflag           = 1
            emark(ir)       = 0

        ENDIF

      ENDIF
c
c
c     Check to see if the mock power term is getting rediculously large. At
c     the moment this occurs when there is flow reversal such that flow to
c     the offending target extends beyond the majority of the electron
c     power loss, so that the routine that ensures that the power flow reversal
c     point is outside the particle flow reversal point essentially stifles
c     the growth of the mock power term.  A more sophisticated solution would
c     be to check the integral of the electron power loss term in the flow-to-target
c     region to make sure that it does not contain the majority of the power loss
c     for the half-ring:
c
      IF (osm_powopt.EQ.2.AND.
     .    ((.NOT.restrictosmpmk.AND.osm_peimul(region,ir).GT.10.0  ).OR.
     .     (     restrictosmpmk.AND.osmpeires (region,ir).GT.10.0D0))
     .   ) THEN
        note            = 'OVERFLOW'
        code(region,ir) = 20
        status          = 0
        flag            = 2

        eflag           = 1
        emark(ir)       = 0
      ENDIF


c
c
c     Process error:
c
c
      IF (eflag.EQ.1) THEN

        WRITE(0,*) note(1:LEN_TRIM(note))

        IF (restrictosmpmk) STOP 'DEVELOPMENT NEEDED A364'

        IF (emark(ir).EQ.0) THEN
c
c          WRITE(0     ,*) 'Resetting mock power source'
c          WRITE(PINOUT,*) 'Resetting mock power source'
c


c          IF     (errorc(IKLO).GE.1.AND.errorc(IKHI).LE.0) THEN

           WRITE(0,*) 'FLOW REVERSAL : I,O =',
     .       IdentifyFlowReversal(IKLO,ir),
     .       IdentifyFlowReversal(IKHI,ir)
           WRITE(PINOUT,*) 'FLOW REVERSAL : I,O =',
     .       IdentifyFlowReversal(IKLO,ir),
     .       IdentifyFlowReversal(IKHI,ir)


           IF     (stopopt2.NE.8.AND.
     .                  (IdentifyFlowReversal(IKLO,ir).GE.
     .                          IdentifyFlowReversal(IKHI,ir).AND.
     .                          IdentifyFlowReversal(IKLO,ir).NE.0.0)
     .            ) THEN

              WRITE(0     ,*)
              WRITE(0     ,*) ' !!!Reducing lower ionisation source!!!'
              WRITE(0     ,*)
              WRITE(PINOUT,*) 'Reducing lower ionisation source'

              osm_cadj(1,region,ir) = 1

              osm_fixopt(IKLO,ir) = 0
              osm_fixopt(IKHI,ir) = 0

              WRITE(0,*)
              WRITE(0,*) '***TEMP WHIPE 01***'
              WRITE(PINOUT,*) '***TEMP WHIPE 01***'
              WRITE(0,*)

              DO ik = 1, nks(ir)
                osm_dp6(ik,ir) = 1.00
              ENDDO

              osm_peimul(IKLO,ir) = 0.1
              osm_peimul(IKHI,ir) = 0.1

              IF (restrictosmpmk) STOP 'DEVELOPMENT NEEDED A324'

              IF (switch(SWPOW3).EQ.4.0.OR.switch(SWPOW3).EQ.5.0.OR.
     .            switch(SWPOW3).EQ.1.0.OR.switch(SWPOW3).EQ.6.0.OR.
     .            switch(SWPOW3).EQ.7.0) THEN
                WRITE(0,*) 'ASSIGNING MOCK POWER EXTENT'

                CALL AssignMockPowerExtent(IKLO,ir)
                CALL AssignMockPowerExtent(IKHI,ir)
              ENDIF


              WRITE(79,'(A,1X,3I4,1X,I4,E12.4)') '''LOWER ION 1.00''',
     .          rel_step,rel_iter,rel_count,ir,0.9

              IF     (osm_code(IKLO,ir).EQ.1) THEN
                DO ik = ikbound(ir,IKLO), SymmetryPoint(ir)
                  pinion(ik,ir) = 0.9 * pinion(ik,ir)
                  pinqe (ik,ir) = 0.9 * pinqe (ik,ir)
                ENDDO
              ELSEIF (osm_code(IKLO,ir).EQ.2) THEN
                DO ik = SymmetryPoint(ir), ikbound(ir,IKLO), -1
                  IF (kss(ik,ir).GE.osm_dp5(IKLO,ir)) ik1 = ik
c                  WRITE(PINOUT,*) 'BUG? ',ik,ik1,ir,
c     .                       kss(ik,ir),osm_dp5(IKLO,ir)
                ENDDO


                ik = ik1
                DO WHILE (ABS(pinqe(ik,ir)).GT.0.20*ABS(pinqe(ik1,ir)))
c                DO WHILE (ABS(pinqe(ik,ir)).GT.0.05*ABS(pinqe(ik1,ir)))
c                  WRITE(PINOUT,*) 'BUG + ',ik,ik1,ir,
c     .                       ABS(pinqe(ik,ir)),ABS(pinqe(ik1,ir))
                  ik = ik - 1
                ENDDO

                ik1 = ik

                WRITE(0     ,*) 'REDUCING INNER IONISATION : IK1 = ',ik1
                WRITE(PINOUT,*) 'REDUCING INNER IONISATION : IK1 = ',ik1

                DO ik = ik1, SymmetryPoint(ir)
                  pinion(ik,ir) = 0.9 * pinion(ik,ir)
                  pinqe (ik,ir) = 0.9 * pinqe (ik,ir)
                ENDDO

c                DO ik = ikbound(ir,IKLO), SymmetryPoint(ir)
c                  pinion(ik,ir) = 0.9 * pinion(ik,ir)
c                  pinqe (ik,ir) = 0.9 * pinqe (ik,ir)
c                ENDDO
              ELSE
                CALL ER('WWW1','Unknown low index OSM_CODE',*99)
              ENDIF

c            ELSE
c              WRITE(0     ,*) 'FLOW REVERSAL NOT FOUND ON INNER LEG'
c              WRITE(PINOUT,*) 'FLOW REVERSAL NOT FOUND ON INNER LEG'
c
c              emark(ir) = 4
c            ENDIF






            ELSEIF (stopopt2.NE.8.AND.
     .                (IdentifyFlowReversal(IKHI,ir).GT.
     .                           IdentifyFlowReversal(IKLO,ir).AND.
     .                           IdentifyFlowReversal(IKHI,ir).NE.0.0)
     .             ) THEN

              WRITE(0     ,*)
              WRITE(0     ,*) ' !!!Reducing higher ionisation source!!!'
              WRITE(0     ,*)
              WRITE(PINOUT,*) 'Reducing higher ionisation source'

              osm_cadj(1,region,ir) = 1

              WRITE(0,*)
              WRITE(0,*) '***TEMP WHIPE 02***'
              WRITE(PINOUT,*) '***TEMP WHIPE 02***'
              WRITE(0,*)

              DO ik = 1, nks(ir)
                osm_dp6(ik,ir) = 1.00
              ENDDO

              osm_peimul(IKLO,ir) = 0.1
              osm_peimul(IKHI,ir) = 0.1

              IF (restrictosmpmk) STOP 'DEVELOPMENT NEEDED A643'

              IF (switch(SWPOW3).EQ.4.0.OR.switch(SWPOW3).EQ.5.0.OR.
     .            switch(SWPOW3).EQ.1.0.OR.switch(SWPOW3).EQ.6.0.OR.
     .            switch(SWPOW3).EQ.7.0) THEN
                WRITE(0,*) 'ASSIGNING MOCK POWER EXTENT'

                CALL AssignMockPowerExtent(IKLO,ir)
                CALL AssignMockPowerExtent(IKHI,ir)
              ENDIF

              osm_fixopt(IKLO,ir) = 0
              osm_fixopt(IKHI,ir) = 0




              WRITE(79,'(A,1X,3I4,1X,I4,E12.4)') '''HIGHER ION1.00''',
     .          rel_step,rel_iter,rel_count,ir,0.9

              IF     (osm_code(IKHI,ir).EQ.1) THEN
                DO ik = ikbound(ir,IKHI), SymmetryPoint(ir)+1, -1
                  pinion(ik,ir) = 0.9 * pinion(ik,ir)
                  pinqe (ik,ir) = 0.9 * pinqe (ik,ir)
                ENDDO
              ELSEIF (osm_code(IKHI,ir).EQ.2) THEN
                DO ik = SymmetryPoint(ir)+1, ikbound(ir,IKLO)
                  IF (kss(ik,ir).LT.osm_dp5(IKLO,ir)) ik1 = ik
                ENDDO

                ik = ik1
                DO WHILE (ABS(pinqe(ik,ir)).GT.0.05*ABS(pinqe(ik1,ir)))
                  ik = ik + 1
                ENDDO

                ik1 = ik

                WRITE(0     ,*) 'REDUCING OUTER IONISATION : IK1 = ',ik1
                WRITE(PINOUT,*) 'REDUCING OUTER IONISATION : IK1 = ',ik1

                DO ik = ik1, SymmetryPoint(ir)+1, -1
                  pinion(ik,ir) = 0.9 * pinion(ik,ir)
                  pinqe (ik,ir) = 0.9 * pinqe (ik,ir)
                ENDDO

              ELSE
                CALL ER('WWW1','Unknown high index OSM_CODE',*99)
              ENDIF

              emark(ir) =  0
              riter(ir) = -1

              CALL AssignMockPowerExtent(IKLO,ir)

            ELSEIF (stopopt2.EQ.8) THEN
c...          Increase the upstream temperature in order to require that more 
c             power enter the ring:

              rel_mfsp(IKLO,ir) = 1.05 * rel_mfsp(IKLO,ir)
              rel_mfsp(IKHI,ir) =        rel_mfsp(IKLO,ir)
              WRITE(0,*) 'INCREASE REL_MFSP BY 5%',
     .                   region,ir,rel_mfsp(IKLO,ir)

              emark(ir) =  0
              riter(ir) = -1

            ELSE
              WRITE(0     ,*) 'FLOW REVERSAL NOT FOUND'
              WRITE(PINOUT,*) 'FLOW REVERSAL NOT FOUND'

              emark(ir) = 4
            ENDIF

c          ELSE
c            emark(ir) = 4
c          ENDIF



          IF (emark(ir).EQ.0) THEN
c...temp1
c            osm_peimul(IKLO,ir) = 0.1
c            osm_peimul(IKHI,ir) = 0.1

            rel_prbfr(IKLO,ir) = 0.1
            rel_prbfr(IKHI,ir) = 0.1
c...eh?
            IF (GetModel(region,ir).EQ.22.0.AND.
     .          switch(SWMACH).NE.0.0) THEN
              IF (region.EQ.IKLO) THEN
                cmachno(ir,2) = initm0
              ELSE
                cmachno(ir,1) = initm0
              ENDIF
            ENDIF

            status = 3

          ELSE
            WRITE(0     ,*) 'Unable to find solution'
            WRITE(PINOUT,*) 'Unable to find solution'

            osm_nflag = -5
          ENDIF

        ELSEIF (emark(ir).EQ.1) THEN

          IF (fmark(ir).EQ.0) THEN
            WRITE(0     ,*) 'Suppressing midstream flow'
            WRITE(PINOUT,*) 'Suppressing midstream flow'

c...temp1
c            osm_peimul(IKLO,ir) = 0.1
c            osm_peimul(IKHI,ir) = 0.1

            temod(ir) = 0.0



c            multi = ABS(intrec1(nks(ir)+1,ir)) /
c     .              ABS(intion1(nks(ir)+1,ir))

c             multi = multi * 0.5
c             CALL ChopSource(-1,ir,pinion,multi)



            multi = (ABS(knds(idds(ir,2))*kvds(idds(ir,2))) +
     .               intrec1(SymmetryPoint(ir),ir)) /
     .               intion1(SymmetryPoint(ir),ir)



            WRITE(0     ,*) 'CHOPPING! multi = ',multi
c            WRITE(SLOUT ,*) 'CHOPPING! multi = ',multi
            WRITE(PINOUT,*) 'CHOPPING! multi = ',multi

            CALL ChopSource(IKLO,ir,pinion,multi)

            osm_ionmul(IKLO,ir) = multi


            multi = (ABS(knds(idds(ir,1))*kvds(idds(ir,1))) +
     .               (intrec1(nks(ir)+1        ,ir) -
     .                intrec1(SymmetryPoint(ir),ir))) /
     .               (intion1(nks(ir)+1        ,ir) -
     .                intion1(SymmetryPoint(ir),ir))

            WRITE(0     ,*) 'CHOPPING! multi = ',multi
c            WRITE(SLOUT ,*) 'CHOPPING! multi = ',multi
            WRITE(PINOUT,*) 'CHOPPING! multi = ',multi


            CALL ChopSource(IKHI,ir,pinion,multi)

            osm_ionmul(IKHI,ir) = multi

            fmark(ir) = 1
            status = 3

          ELSEIF (fmark(ir).EQ.1) THEN

            IF (errorc(region).EQ.0.AND.errorc(region2).EQ.0) THEN
              multi = (1.0 / osm_ionmul(IKLO,ir))**(1.0 / 10.0)

              WRITE(0     ,*) 'SPIFFING! multi = ',multi
c              WRITE(SLOUT ,*) 'SPIFFING! multi = ',multi
              WRITE(PINOUT,*) 'SPIFFING! multi = ',multi

              CALL ChopSource(IKLO,ir,pinion,multi)

              multi = (1.0 / osm_ionmul(IKHI,ir))**(1.0 / 10.0)

              WRITE(0     ,*) 'SPIFFING! multi = ',multi
c              WRITE(SLOUT ,*) 'SPIFFING! multi = ',multi
              WRITE(PINOUT,*) 'SPIFFING! multi = ',multi

              CALL ChopSource(IKHI,ir,pinion,multi)

              status = 3
            ELSE
              WRITE(0     ,*) 'Darn!'
              WRITE(PINOUT,*) 'Darn!'


              fmark(ir) = 0
              emark(ir) = 4

              osm_nflag = -5
            ENDIF

          ELSEIF (fmark(ir).EQ.2) THEN
            WRITE(0     ,*) 'ERROR RETURNED FROM RELAXATION ROUTINE'
            WRITE(PINOUT,*) 'ERROR RETURNED FROM RELAXATION ROUTINE'

            emark(ir) = 2
            osm_nflag = -5
          ELSE
            STOP 'NASTY ERROR 001'
          ENDIF


        ELSEIF (emark(ir).EQ.2) THEN

          WRITE(0     ,*) 'Reducing relaxation fraction'
          WRITE(PINOUT,*) 'Reducing relaxation fraction'

            WRITE(0,*) 'MARK: LGMUL2'

          osm_peimul(IKLO,ir) = lgmul(IKLO,ir)
          osm_peimul(IKHI,ir) = lgmul(IKHI,ir)

          temod(ir) = 0.0

          status    = -2
          emark(ir) =  3

        ELSEIF (emark(ir).EQ.3) THEN

          IF (lgmul(IKLO,ir).NE.0.0.AND.lgmul(IKHI,ir).NE.0.0) THEN

            WRITE(0     ,*) 'Increasing effective probe temperature'
            WRITE(PINOUT,*) 'Increasing effective probe temperature'

            WRITE(0,*) 'MARK: LGMUL3'

            osm_peimul(IKLO,ir) = lgmul(IKLO,ir)
            osm_peimul(IKHI,ir) = lgmul(IKHI,ir)

            temod(ir) = temod(ir) + MAX(TADJUST,
     .                                  MAX(lgte(IKLO,ir)-temark,
     .                                      lgte(IKHI,ir)-temark))

            status = 3

          ELSE

            WRITE(0     ,*) 'No base solution available'
            WRITE(PINOUT,*) 'No base solution available'

            osm_nflag = -5

          ENDIF

          emark(ir) = 4

        ELSE

          WRITE(0     ,*) 'Unable to find solution'
          WRITE(PINOUT,*) 'Unable to find solution'

          osm_nflag = -5

        ENDIF

        riter(ir) = -1

      ENDIF

c
c
c

      IF (rel_prbfr(region,ir).GT.20.0) THEN

         note = note(1:LEN_TRIM(note))//' [MAX]'

c...     Updated this for the attached (maybe) C-Mod case 990429019
c        but this may cause problems for 981116027 modelling:
         rel_prbfr(region,ir) = 50.0
c         rel_prbfr(region,ir) = 20.0

      ENDIF

      IF (      rel_prbfr(region ,ir).GT.
     .    100.0*rel_prbfr(region2,ir).AND.
     .    code(region2,ir).NE.6) THEN

         note = note(1:LEN_TRIM(note))//' [SLOW]'

         rel_prbfr(region,ir) = 10.0 * rel_prbfr(region2,ir)

      ENDIF


c
c
c

      IF (flag.EQ.1) THEN
        qem1   = SNGL(qem2  (region,ir))
        qem3   =      qem2  (region,ir)
        norm1  = norm2 (region,ir)
        signe1 = signe2(region,ir)
      ELSE
        qem1   = qem
        qem3   = qem0
        norm1  = norm
        signe1 = signe
      ENDIF

c      WRITE(0,*) 'HELP:',flag,norm1,norm1,qem1

c
c
c

        lockionptarg = .FALSE.

        IF (code(IKLO,ir).EQ.6.AND.osm_model(IKHI,ir).EQ.22) THEN
c          osm_mode = 2
          lockionptarg = .TRUE.
        ENDIF

c        WRITE(0     ,*) 'KTE1 C4=',osm_sympt(ir),ktebs(osm_sympt(ir),ir)
c     .    ,cmachno(ir,1)
c        WRITE(PINOUT,*) 'KTE1 C4=',osm_sympt(ir),ktebs(osm_sympt(ir),ir)
c
c      CALL INITPLASMA(ir,ir,3)
c      CALL SOL_PLASMA(ir,ir,3)
c      CALL SOL       (ir,ir,3)
c
c        WRITE(0     ,*) 'KTE1 D4=',osm_sympt(ir),ktebs(osm_sympt(ir),ir)
c     .   ,status
c     .    ,cmachno(ir,1),elecptarg(ir,3),ionptarg(ir,3)
c        WRITE(PINOUT,*) 'KTE1 D4=',osm_sympt(ir),ktebs(osm_sympt(ir),ir)
c     .   ,status


      IF (status.EQ.1.AND.
     .    .NOT.(osm_matcht.EQ.4.AND.region.EQ.IKHI.AND.
     .          osm_model(IKHI,ir).EQ.22)) THEN
c...little safety feature:
        IF (rel_prbfr(region,ir).LT.1.0E-10) THEN
          deltam = DBLE(rel_prbfr(region,ir))
        ELSE
          deltam = MIN(10.0D0,DBLE(norm1) * DBLE(rel_prbfr(region,ir)))
        ENDIF

c        WRITE(0,*) 'PRBFR  = ',rel_prbfr(region,ir)
c        WRITE(0,*) 'NORM1  = ',norm1
c        WRITE(0,*) 'SIGNE1 = ',signe1
c        WRITE(0,*) 'DELTAM = ',deltam
c        WRITE(0,*) 'QEM1   = ',qem1

        IF     (osm_powopt.EQ.1) THEN
          qem1 = SNGL(DBLE(qem1) + DBLE(signe1) * deltam)
          qem3 =           qem3  + DBLE(signe1) * deltam 
        ELSEIF (osm_powopt.EQ.2) THEN
          IF (flag.EQ.3) THEN
            qem1 = lgmul (region,ir)
            qem3 = lgmul2(region,ir)
          ELSE

c            WRITE(0,*) '******** deltam:',deltam

            qem1 = SNGL(DBLE(qem1) - DBLE(signe1) * deltam)
            qem3 =           qem3  - DBLE(signe1) * deltam 
          ENDIF
        ELSE
          CALL ER('MatchProbe','Invalid OSM_POWOPT value',*99)
        ENDIF


c        WRITE(0,*) 'QEM1   = ',qem1
      ENDIF
c
c
c
c      WRITE(0,*) 'HELP2:',qem1,osm_peimul(IKLO,ir),status

      IF (status.EQ.1) THEN

        IF (osm_powopt.EQ.1) THEN
           STOP 'REALLY OLD METHOD'
          IF     (flag.GE.0.AND.region.EQ.IKLO) THEN
            CALL ReversePINQeMultiplier(IKLO,ir)
            rel_qemul(IKLO,ir) = qem1
            CALL ApplyPINQeMultiplier(IKLO,ir)
          ELSEIF (flag.GE.0.AND.region.EQ.IKHI) THEN
            CALL ReversePINQeMultiplier(IKHI,ir)
            rel_qemul(IKHI,ir) = qem1
            CALL ApplyPINQeMultiplier(IKHI,ir)
          ENDIF
        ELSEIF (.TRUE.) THEN


c...temp1
c          IF (errorc(region).GT.1.AND.osm_code(region,ir).EQ.2.AND.
c     .        osm_fixopt(region,ir).EQ.0) qem1 = 0.0

          IF (restrictosmpmk) THEN
            osmpeires  (region,ir) = qem3
          ELSE

c...GOG
c            IF (ABS(0.9 * qem).GT.ABS(qem1)) THEN
c              WRITE(0     ,*) 'CLIPPERS:',qem,qem1
c              WRITE(PINOUT,*) 'CLIPPERS:',qem,qem1
c              qem1 = SIGN(1.0,qem1) * ABS(0.9 * qem)
c              qem3 = DBLE(qem1)
c              rel_prbfr(region,ir) = tmpprbfr
c              WRITE(0     ,*) 'CLIPPERS:',qem,qem1
c              WRITE(PINOUT,*) 'CLIPPERS:',qem,qem1
c            ENDIF
            IF (.FALSE..AND.
     .          ABS(qem).LT.0.1.AND.ABS(1.1 * qem).LT.ABS(qem1)) THEN
              WRITE(0     ,*) 'CLIPPERS:',qem,qem1
              WRITE(PINOUT,*) 'CLIPPERS:',qem,qem1
              qem1 = SIGN(1.0,qem1) * ABS(1.1 * qem)
              qem3 = DBLE(qem1)
              rel_prbfr(region,ir) = tmpprbfr
              WRITE(0     ,*) 'CLIPPERS:',qem,qem1
              WRITE(PINOUT,*) 'CLIPPERS:',qem,qem1
            ENDIF

            osm_peimul (region,ir) = qem1
            osm_peimul2(region,ir) = qem3
          ENDIF


c          WRITE(0,*) 'HELP3:',qem1,osm_peimul(IKLO,ir)
c          WRITE(0,*) 'MARK: QEM,QEM1= ',qem,qem1

c...Adjust osm_dp6, the local mock power multiplier:
c...Don't do this if there is an error, giving the algoritm a chance
c   to recover the solution:
c          IF ((SIGN(1.0,qem).EQ.SIGN(1.0,qem1).AND.
c     .        ABS(qem).LT.ABS(qem1).AND.qem.NE.0.0.AND.qem1.NE.0.0).AND.
c     .         errorc(region).LE.1) THEN
c          IF (SIGN(1.0,qem).EQ.SIGN(1.0,qem1).AND.
c     .        ABS(qem).LT.ABS(qem1).AND.qem.NE.0.0.AND.qem1.NE.0.0) THEN

c          WRITE(0,*) '--??',ir,region,ikerror(region,ir)



c          WRITE(0,*) 'SIGNS:',SIGN(1.0,qem),SIGN(1.0,qem1)


          IF (SIGN(1.0,qem).EQ.SIGN(1.0,qem1).AND.
     .        qem.NE.0.0.AND.qem1.NE.0.0.AND.
c...new
     .        .NOT.(errorc(region).GT.1.AND.
     .              osm_dp6(ikerror(region,ir),ir).NE.1.0).AND.

c... **** THIS IS A POOR CONDITION (OSM_MATCHT.EQ.4) : BASICALLY I JUST WANTED
c         TO GET THIS OUT OF HERE QUICK ****

     .         .NOT.(osm_matcht.EQ.4)) THEN
c     .         .NOT.(osm_matcht.EQ.4.AND..NOT.firstshape)) THEN
c     .        osm_dp6(ikerror(region,ir),ir).EQ.1.0) THEN


c            STOP 'STOP : SHOULD NOT BE HERE'

            IF     (region.EQ.IKLO) THEN
              ik1 = 1
              ik2 = SymmetryPoint(ir)
            ELSEIF (region.EQ.IKHI) THEN
              ik1 = SymmetryPoint(ir) + 1
              ik2 = nks(ir)
            ENDIF

            frac = qem1 / qem

            WRITE(0     ,*) 'TRYING TO SCALE OSM_DP6'
            WRITE(PINOUT,*) 'TRYING TO SCALE OSM_DP6'

            DO ik = ik1, ik2
              IF (osm_dp6(ik,ir).NE.1.0) THEN
                  osm_dp6(ik,ir) = osm_dp6(ik,ir) / frac
                  WRITE(PINOUT,*) 'MARK: SCALING OSM_DP6 MULTIPLIER ',
     .                            ik,ir,frac
c                  WRITE(0     ,*) 'MARK: SCALING OSM_DP6 MULTIPLIER ',
c     .                            ik,ir,frac
              ENDIF
            ENDDO

          ELSEIF (SIGN(1.0,qem).NE.SIGN(1.0,qem1).AND.
     .            qem.NE.0.0.AND.qem1.NE.0.0.AND.
     .            .NOT.restrictosmpmk) THEN
            IF     (region.EQ.IKLO) THEN
              ik1 = 1
              ik2 = osm_sympt(ir)
            ELSEIF (region.EQ.IKHI) THEN
              ik1 = osm_sympt(ir) + 1
              ik2 = nks(ir)
            ENDIF

            IF (ir.GT.irsep.AND.ir.LT.irwall) THEN
              frac = ABS(qem1 / qem)
            ELSE
              frac = 1.0
            ENDIF

            WRITE(0     ,*) 'MARK: TROUBLE DETECTED: ',qem1,qem,frac
            WRITE(PINOUT,*) 'MARK: TROUBLE DETECTED: ',qem1,qem,frac


c            WRITE(0,*) 'SETTING OSM_MODE TO 2'
c            osm_mode = 2
            

            DO ik = ik1, ik2
              IF (osm_dp6(ik,ir).NE.1.0) THEN
c                WRITE(0,*) '-->:',SIGN(1.0,osm_dp6(ik,ir)),
c     .                     SIGN(1.0,qem1)
                IF     (SIGN(1.0,osm_dp6(ik,ir)).EQ.SIGN(1.0,qem1)) THEN
                  osm_dp6(ik,ir) = -osm_dp6(ik,ir) / frac
                  WRITE(0     ,*) 'MARK: INVERTING OSM_DP6 MULTIPLIER ',
     .                            ik,ir
                  WRITE(PINOUT,*) 'MARK: INVERTING OSM_DP6 MULTIPLIER ',
     .                            ik,ir
                ELSEIF (SIGN(1.0,osm_dp6(ik,ir)).NE.SIGN(1.0,qem1)) THEN
                  osm_dp6(ik,ir) = -osm_dp6(ik,ir) / frac
                  WRITE(0     ,*) 'MARK: -INVERTING OSM_DP6 MULTIPLIER',
     .                            ik,ir
                  WRITE(PINOUT,*) 'MARK: -INVERTING OSM_DP6 MULTIPLIER',
     .                            ik,ir
                ENDIF
              ENDIF
            ENDDO

            WRITE(0,*) 'PRESERVING PEAK REGION'
            note = note(1:LEN_TRIM(note))//' [JOGGED]'
            osm_peimul(region,ir) = -qem           

          ELSE

          ENDIF


c...      The above scaling keeps the mock power term constant over the
c         range of OSM_DP6.NE.1.0, but it fails when the change in the 
c         cross-field power term is sufficient to generate negative
c         temperatures in this region.  Indeed, the above code stops
c         power loss increases via the mock power term from recovering
c         the solution.  So, one needs to account for the change in 
c         the cross-field power term and account for that by scaling
c         OSM_DP6 appropriately:

          IF (rel_prbfr(region,ir).LT.1.0E-10) THEN

            WRITE(0     ,*) 'REL_PRBFR SMALL, NOT SCALING PLATEAU '//
     .                      'REGION!'
            WRITE(PINOUT,*) 'REL_PRBFR SMALL, NOT SCALING PLATEAU '//
     .                      'REGION!'

          ELSEIF (osm_matcht.EQ.4.AND.osmmock.EQ.1.AND.
     .        osm_model(region,ir).EQ.24.AND.
     .        SIGN(1.0,qem).EQ.SIGN(1.0,qem1).AND.
     .        .NOT.(errorc(region).GT.1.AND.
     .              ((region.EQ.IKLO.AND.
     .                ikerror(IKLO,ir).LT.ikfluid(IKLO,ir)).OR.
     .               (region.EQ.IKHI.AND.
     .                ikerror(IKHI,ir).GT.ikfluid(IKHI,ir)))).AND.
     .        .NOT.firstshape.AND..NOT.restrictosmpmk) THEN
c     .        errorc(region).LE.1.AND..NOT.firstshape) THEN


             WRITE(PINOUT,*) 'FANCY OSM_DP6 ADJUSTING'

c          IF (region.EQ.IKLO.AND.osm_matcht.EQ.4.AND.osmmock.EQ.1.AND.
c     .        osm_model(IKLO,ir).EQ.24.AND.
c     .        SIGN(1.0,qem).EQ.SIGN(1.0,qem1))THEN

c...        Assume that the cross-field power term is distributed
c           uniformly.  Some tricky business because the true mock
c           power term is staggered between cell centers:

c...        Find bounds for OSM_DP6.NE.1.0 (not all that efficient 
c           at the moment -- should really be setting the bound 
c           elsewhere):

            ik1 = MAXNKS
            ik2 = 1

            IF     (region.EQ.IKLO) THEN

              DO ik = ikbound(ir,IKLO), osm_sympt(ir)               
                IF (osm_dp6(ik,ir).NE.1.0.AND.ik1.EQ.MAXNKS) ik1 = ik
                IF (osm_dp6(ik,ir).NE.1.0.AND.ik1.NE.MAXNKS) ik2 = ik
              ENDDO

              pmktot = intpmk1(osm_sympt(ir),ir) - intpmk1(ik2+1,ir)

            ELSE

              DO ik = osm_sympt(ir)+1, ikbound(ir,IKHI)-1
                IF (osm_dp6(ik,ir).NE.1.0.AND.ik1.EQ.MAXNKS) ik1 = ik
                IF (osm_dp6(ik,ir).NE.1.0.AND.ik1.NE.MAXNKS) ik2 = ik
              ENDDO

              pmktot = intpmk1(ik1,ir) - intpmk1(osm_sympt(ir)+1,ir) 
c              pmktot = intpmk1(ik2,ir) - intpmk1(osm_sympt(ir)+1,ir) 

            ENDIF

            IF    (osm_model(IKLO,ir).EQ.24.AND.
     .             osm_model(IKHI,ir).EQ.24) THEN
              spcftot = ksb(ikbound(ir,IKHI)  ,ir)-
     .                  ksb(ikbound(ir,IKLO)-1,ir)
            ELSEIF(osm_model(IKLO,ir).EQ.24.AND.
     .             osm_model(IKHI,ir).EQ.22) THEN
              spcftot = ksmaxs(ir) - ksb(ikbound(ir,IKLO)-1,ir)
            ELSEIF(osm_model(IKLO,ir).EQ.22.AND.
     .             osm_model(IKHI,ir).EQ.24) THEN
              spcftot = ksb(ikbound(ir,IKHI),ir)
            ELSE
              STOP 'STOP: UNSUPPORTED SOLVER CONFIGURATION'
            ENDIF
          
            DO ik = ikbound(ir,IKLO), ikbound(ir,IKHI)
c            DO ik = ik1, ik2
              IF (osm_dp6(ik,ir).EQ.1.0) CYCLE


c...Well, this doesn't quite work properly, since changin OSM_DP6 to account
c   for a change in the cross-field power term cannot be done here without
c   changing the cross-field power term again.  So, some iterative (or more clever
c   than I can think of at the moment) method is needed.  Anyway, it turns out
c   that if one is careful with the power distribution then the present approximation
c   is good enough.  It is also not much of a concern if the peak width is smaller
c   than the PINQE decay length:

c              spcf = kss(ik+1,ir) - kss(ik,ir)
c              dpcf = pmktot2 * (qem1 / qem - 1.0) * (spcf / spcftot)
c              pmk = (intpmk1(ik+1,ir) - intpmk1(ik,ir))
c              frac = (1.0 + dpcf / pmk) / (qem1 / qem)
c              osm_dp6(ik,ir) = osm_dp6(ik,ir) * frac
c              WRITE(0,*) 'DPCF:',dpcf,frac-1.0,pmk,dpcf/pmk
c
c !!!         jdemod - uncommenting the assignment of dpcf and pmk so 
c                      they don't generate error messages
c                      at compilation
c
              pmk  = 1.0
              dpcf = 0.0
c
c !!!         jdemod 
c
              spcf = kss(ik+1,ir) - kss(ik,ir)
              qrat1 = qem1 / qem
              srat1 = spcftot / spcf
              prat1 = pmktot / (intpmk1(ik+1,ir) - intpmk1(ik,ir)) 
              frac  = (((qrat1 - 1.0) / (srat1 - 1.0)) * prat1 + 1.0) /
     .                qrat1

c              osm_dp6(ik,ir) = osm_dp6(ik,ir) * frac 

              iksft = 0
              IF (ik.LT.osm_sympt(ir)) iksft = 1

              IF (flag.EQ.3) 
     .          WRITE(PINOUT,'(3X,A,2I6,2F8.3,1P,5E11.3,0P,F8.4)')
     .            'HAIL:',ik,ir,spcf,spcftot,dpcf,
     .            pmk,pmktot,frac-1.0,osm_dp6(ik,ir),
     .                ktebs(ik+iksft,ir)  

c              WRITE(0     ,'(3X,A,2I6,2F8.3,1P,5E11.3,0P,F8.4)')
c     .          'HAIL:',ik,ir,spcf,spcftot,dpcf,
c     .          pmk,pmktot,frac-1.0,osm_dp6(ik,ir),
c     .              ktebs(ik+iksft,ir)
            ENDDO


            WRITE(PINOUT,*) 'MAINTAINING OSM_DP6 IN MATCHPROBE'

            CALL CalcLocalMockAdjustment2(
     .           DBLE(qem1/qem-1.0),DBLE(pmktot),ir,dp6adj)
	    
            DO ik = ikbound(ir,IKLO), ikbound(ir,IKHI)
c            DO ik = ik1, ik2
              IF (osm_dp6(ik,ir).EQ.1.0) CYCLE
	    
c              WRITE(0,*) '-->',ik,dp6adj(ik)
c              WRITE(0,*) '-->',ik,((1.0+dp6adj(ik))/(qem1/qem))-1.0
	    
c              osm_dp6(ik,ir) = osm_dp6(ik,ir) * (1.0 + dp6adj(ik)) 
	    
              osm_dp6(ik,ir) = osm_dp6(ik,ir) * (1.0 + dp6adj(ik)) / 
     .                         (qem1 / qem)
	    
	    
            ENDDO



          ENDIF

        ENDIF
      ENDIF

c      WRITE(0,*) 'HELP3:',qem1,osm_peimul(IKLO,ir)

c
c
c
c
c
c
c
       IF (osm_mode.GE.1)
     .  WRITE(PINOUT,'(1X,A,2I4,I5,I4,F12.4,2F15.9,1X,I3,1P,E10.2,0P,'//
     .               ' F6.1,2X,2F23.15,1X,2I4,1X,A)')
     .    tag(region),ir,ecount(ir),ii(region,ir),ik1,
     .    norm*100.0,hte1,hte2,
     .    INT(signe),rel_prbfr(region,ir),temod(ir),
     .    osm_peimul2(region,ir),osmpeires(region,ir),
     .    errorc(IKLO),errorc(IKHI),note(1:LEN_TRIM(note))

       IF (osm_mode.GE.1)
     .  WRITE(0     ,'(1X,A,2I4,2F15.9,I4,F11.5,I3,1P,E10.2,0P,F6.1,'//
     .               ' 1X,2F23.15,1X,A)')
     .    tag(region),ir,ecount(ir),hte1,hte2,ii(region,ir),norm*100.0,
     .    INT(signe),rel_prbfr(region,ir),temod(ir),
     .    osm_peimul2(region,ir),osmpeires(region,ir),
     .    note(1:LEN_TRIM(note))


c        WRITE(0     ,*) 'KTE1 C3=',osm_sympt(ir),ktebs(osm_sympt(ir),ir)
c     .    ,cmachno(ir,1)
c        WRITE(PINOUT,*) 'KTE1 C3=',osm_sympt(ir),ktebs(osm_sympt(ir),ir)
c
c      CALL INITPLASMA(ir,ir,3)
c      CALL SOL_PLASMA(ir,ir,3)
c      CALL SOL       (ir,ir,3)
c
c        WRITE(0     ,*) 'KTE1 D3=',osm_sympt(ir),ktebs(osm_sympt(ir),ir)
c     .    ,cmachno(ir,1),elecptarg(ir,3),ionptarg(ir,3)
c        WRITE(PINOUT,*) 'KTE1 D3=',osm_sympt(ir),ktebs(osm_sympt(ir),ir)


c        WRITE(0,*) 'CODE:',code(region,ir)

c
c     Check for slam dunk:
c

c      IF (region.EQ.IKHI.AND.ii(IKHI,ir).EQ.1.AND.
c     .    (code(IKLO,ir).NE.6.OR.
c     .     code(IKHI,ir).NE.6).AND.
c     .    osm_fixopt(IKLO,ir).NE.2.AND.
c     .    osm_fixopt(IKHI,ir).NE.2.AND..FALSE.) THEN
c
c         WRITE(0,*) '***TEMP WHIPE 06***'
c         WRITE(PINOUT,*) '***TEMP WHIPE 06***'
c         WRITE(0,*)
c
c         DO ik = 1, nks(ir)
c           osm_dp6(ik,ir) = 1.00
c         ENDDO
c
c         rel_prbfr(IKLO,ir) = 1.0
c         rel_prbfr(IKHI,ir) = 1.0
c
c         osm_peimul(IKLO,ir) = 0.1
c         osm_peimul(IKHI,ir) = 0.1
c
c         IF (restrictosmpmk) STOP 'DEVELOPMENT A329'
c
c       ENDIF






c
c     Increment counter and store data (unless an error is registered):
c
      IF (flag.EQ.-1.OR.flag.EQ.0.OR.flag.EQ.2) THEN
        ii(region,ir) = ii(region,ir) + 1

        norm2 (region,ir)  = norm
        qem2  (region,ir)  = qem0
        signe2(region,ir)  = signe
        errorc2(region,ir) = errorc(region)
      ENDIF

      osm_temod(ir) = temod(ir)

      timeint(55) = timeint(55) + (Clock2() - timeint(5))

     
c      IF (rel_prbfr(region,ir).LT.1.0E-10) restrictmachno = .TRUE.



      RETURN
99    CONTINUE
      WRITE(0,*) 'OSM_POWOPT=',osm_powopt
      STOP
      END
c
c ======================================================================
c
c subroutine: EvaluatePowerTerms
c
c
c    Event                          PeiMul    dp1        dp2        dp3
c
c 1) Te becomes negative before     < 0.0     -          -          -
c    peak in mock power profile     > 0.0     -          -          decrease
c 2) Te becomes negative after      < 0.0     -          -          -
c    peak in mock power profile     > 0.0     decrease   -          -
c    but within scope of mock
c    power source
c 3) Sonic transition within        < 0.0     increase   -          -
c    scope of mock power profile    > 0.0     ?          -          -
c 4) Sonic transition outside       < 0.0     increase   -          -
c    scope of mock power profile    > 0.0     ?          -          -
c 5) Dip in temperature profile     < 0.0     ?          -          -
c    before peak in mock power      > 0.0     decrease   -          -
c    profile (converged solution)
c 6) Dip in temperature profile     < 0.0     ?          -          -
c    after peak in mock power       > 0.0     -          -          decrease
c    profile (converged solution)
c
c
      SUBROUTINE EvaluatePowerTerms(region,ir,error,result)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INCLUDE 'solparams'
      INCLUDE 'solcommon'
      INCLUDE 'solswitch'

      COMMON /OPTTEMP/ osm_matcht,forcet1
      INTEGER          osm_matcht,forcet1


      INTEGER region,ir,error,result



      INTEGER GetModel,SymmetryPoint
      REAL    GetL1


      INTEGER iklim1,iklim2,i1,action,ik,ik1,ik2,iks,ikm,ike,status
      REAL    rdum1,l1,se,sp,sm,s,qeint,qetot

      INTEGER    INCREASE_DP1  ,DECREASE_DP1  ,INCREASE_DP3   ,
     .           DECREASE_DP3  ,
     .           BEYOND_SCOPE  ,ACTION_UNSPECIFIED    ,NO_ACTION
      PARAMETER (INCREASE_DP1=1,DECREASE_DP1=2,INCREASE_DP3= 3,
     .           DECREASE_DP3=4,
     .           BEYOND_SCOPE=5,ACTION_UNSPECIFIED = 6,NO_ACTION=7)

      result = 1


c      result = 0
c      RETURN


c      WRITE(0,*) '*** NOT EVALUATING POWER TERMS FOR BOOSTING T ***'
c      RETURN

      IF (region.EQ.IKHI.AND.osm_matcht.EQ.4) THEN
        WRITE(PINOUT,*) 'AVOIDING TEMPERATURE PROFILE ANALYSIS ON '//
     .                  'SOL'
        result = 0
        RETURN
      ENDIF

      IF (osmmock.GT.0) THEN
        WRITE(PINOUT,*) 'TEMPERATURE PROFILE ANALYSIS NOT READY '//
     .                  'FOR TRUE MOCK'
        result = 0
        RETURN
      ENDIF


      l1 = GetL1(region,ir)

      rdum1 = (osm_dp2(region,ir) - l1) + 1.5 * osm_dp1(region,ir)

      WRITE(PINOUT,*) 'ADJUST s m = ',serr2(region),rdum1,
     .  '  ',l1,osm_dp2(region,ir)

      se = serr2(region)

      sp = osm_dp2(region,ir) - l1

c      WRITE(0,*) osm_dp2(region,ir),l1


      sm = sp + 5.0 * osm_dp1(region,ir)
c      sm = sp + 2.0 * osm_dp1(region,ir)



      IF     (region.EQ.IKLO) THEN
        sp = osm_dp5(IKLO,ir) / ksmaxs(ir) - l1

        IF (osm_code(IKLO,ir).EQ.1) THEN
          iks = ikbound      (ir,IKLO)
          ikm = SymmetryPoint(ir)
          ike = 0

          qeint = 0.0
          qetot = 0.0

          CALL CalcIntegral4(pinqe,iks,ikm,ir,qetot,2)

          ike = iks+1

          DO ik = iks+1, ikm
            IF (qeint.LT.0.90*qetot) THEN
c            IF (qeint.LT.0.95*qetot) THEN
c...optimise this -- poorly (inefficently) coded:
              CALL CalcIntegral4(pinqe,iks,ik,ir,qeint,2)
              ike = ik
c              WRITE(PINOUT,*) '   >> 95% Qe : ',ik,ike,ir,qeint,qetot
            ENDIF
          ENDDO

          sm = kss(ike,ir) / ksmaxs(ir) - l1

        ELSEIF (osm_code(IKLO,ir).EQ.2) THEN
          sm = sp

        ELSE
          CALL ER('EvaluatePowerTerms','Unknown low index OSM_CODE',*99)
        ENDIF

      ELSEIF (region.EQ.IKHI) THEN

        WRITE(PINOUT,*) 'CHECK = ',ksmaxs(ir),osm_dp5(IKHI,ir),l1,
     .    (ksmaxs(ir)-osm_dp5(IKHI,ir))/ksmaxs(ir)

        sp = (ksmaxs(ir) - osm_dp5(IKHI,ir)) / ksmaxs(ir) - l1

        IF (osm_code(IKHI,ir).EQ.1) THEN
          iks = ikbound      (ir,IKHI)
          ikm = SymmetryPoint(ir)+1
          ike = 0

          qeint = 0.0
          qetot = 0.0

          CALL CalcIntegral4(pinqe,ikm,iks,ir,qetot,2)

          ike = iks - 1

          DO ik = iks-1, ikm, -1
            IF (qeint.LT.0.90*qetot) THEN
c            IF (qeint.LT.0.95*qetot) THEN
              CALL CalcIntegral4(pinqe,ik,iks,ir,qeint,2)
              ike = ik
            ENDIF
          ENDDO

          sm = (ksmaxs(ir) - kss(ike,ir)) / ksmaxs(ir) - l1
        ELSEIF (osm_code(IKLO,ir).EQ.2) THEN
          sm = sp
        ELSE
          CALL ER('EvaluatePowerTerms','Unknown high ind. OSM_CODE',*99)
        ENDIF

      ELSE
        CALL ER('EvaluatePowerTerms','Invalid REGION',*99)
      ENDIF

c
c
c
c
c
      IF (error.EQ.0.OR.(stopopt2.EQ.7.and.error.EQ.1)) THEN

        se = 0.0

        IF     (region.EQ.IKLO) THEN
          ik1 = ikbound(ir,IKLO)
          ik2 = SymmetryPoint(ir)

c...temp
c          IF (status.NE.3) THEN
c            osm_mode = 2
c            CALL INITPLASMA(ir,ir,3)
c            CALL SOL_PLASMA(ir,ir,3)
c            CALL SOL       (ir,ir,3)
c            osm_mode = 1
c          ENDIF


          DO ik = ik2, ik1+1, -1
            IF     (ktebs(ik-1,ir).GT.ktebs(ik,ir)) THEN
              se = kss(ik,ir) / ksmaxs(ir) - l1
c              WRITE(0     ,*) 'TEMPERATURE DIP ',ik,ir,se,ik1,
c     .                                      ktebs(ik-1,ir),ktebs(ik,ir)
c              WRITE(PINOUT,*) 'TEMPERATURE DIP ',ik,ir,se,ik1,
c     .                                      ktebs(ik-1,ir),ktebs(ik,ir)
            ELSEIF (ik.EQ.2.AND.kteds(idds(ir,2)).GT.
     .                          ktebs(ik-1,ir   )) THEN
              se = kss(ik-1,ir) / ksmaxs(ir) - l1            
            ENDIF
          ENDDO

c...mark1
          IF (GetModel(region,ir).EQ.24.AND.
     .        osm_code(region,ir).EQ. 1.AND.
     .        ktebs(ik1+1,ir)-ktebs(ik1,ir).LT.osm_testep) THEN
            se = kss(ik1+1,ir) / ksmaxs(ir) - l1
c            WRITE(0     ,*) 'PROP UP INNER',se,ik1
c            WRITE(PINOUT,*) 'PROP UP INNER',se,ik1
          ENDIF

        ELSEIF (region.EQ.IKHI) THEN
          ik1 = SymmetryPoint(ir) + 1
          ik2 = ikbound(ir,IKHI)

          DO ik = ik1, ik2-1
            IF     (ktebs(ik2,ir).GT.ktebs(ik,ir)) THEN
              se =(ksmaxs(ir) - kss(ik  ,ir)) / ksmaxs(ir) - l1
            ELSEIF (ik.EQ.nks(ir)-1.AND.kteds(idds(ir,1)).GT.
     .                                  ktebs(ik+1,ir   )) THEN
              se =(ksmaxs(ir) - kss(ik+1,ir)) / ksmaxs(ir) - l1
            ENDIF
          ENDDO

c...mark1
          IF (GetModel(region,ir).EQ.24.AND.
     .        osm_code(region,ir).EQ. 1.AND.
c          IF (GetModel(region,ir).EQ.24.AND.
     .        ktebs(ik2-1,ir)-ktebs(ik2,ir).LT.osm_testep) THEN
            se = (ksmaxs(ir) - kss(ik2-1,ir)) / ksmaxs(ir) - l1
            WRITE(0,*) 'PROP UP OUTER'
          ENDIF
        ENDIF
      ENDIF
c
c
c
c
      ike = 0

      IF     (region.EQ.IKLO) THEN
        s = (se + l1) * ksmaxs(ir)

        DO ik = 1, nks(ir)
          IF (s.GE.ksb(ik-1,ir).AND.s.LE.ksb(ik,ir)) ike = ik
        ENDDO

      ELSEIF (region.EQ.IKHI) THEN
        s = ksmaxs(ir) - (se + l1) * ksmaxs(ir)

        DO ik = 1, nks(ir)
          IF (s.GE.ksb(ik-1,ir).AND.s.LE.ksb(ik,ir)) ike = ik
        ENDDO

      ENDIF
c
c
c
c
c
      IF     (error.EQ.4) THEN
c
c
c
        IF     (osm_peimul(region,ir).LT.0.0) THEN
c          action = ACTION_UNSPECIFIED
            action = DECREASE_DP3
        ELSEIF (osm_peimul(region,ir).GE.0.0) THEN
          IF     (se.LT.sp) THEN
            action = DECREASE_DP3
          ELSEIF (se.LT.sm) THEN
            action = DECREASE_DP3
c            action = INCREASE_DP1
c            action = DECREASE_DP1
          ELSE
            action = BEYOND_SCOPE
          ENDIF
        ENDIF

      ELSEIF (error.EQ.1.AND.stopopt2.NE.7) THEN
c
c
c
        IF     (osm_peimul(region,ir).LT.0.0) THEN
          IF     (se.LT.sm) THEN
            action = INCREASE_DP1
          ELSE
            action = BEYOND_SCOPE
          ENDIF
        ELSEIF (osm_peimul(region,ir).GE.0.0) THEN
          action = ACTION_UNSPECIFIED
        ENDIF

      ELSEIF (error.EQ.0.OR.(stopopt2.EQ.7.and.error.EQ.1)) THEN
c
c
c
c        IF     (osm_peimul(region,ir).LT.0.0) THEN
c          action = ACTION_UNSPECIFIED
c        ELSEIF (osm_peimul(region,ir).GE.0.0) THEN
          IF     (se.EQ.0.0) THEN
            action = NULL
          ELSEIF (se.LT.sp) THEN
            action = DECREASE_DP3
          ELSEIF (se.GT.0.0.AND.se.LT.sm) THEN
c...should really have tow options here... one if profiles shouldbe shrunk, the other if expanded
c (negative Te right near peak)...
            action = DECREASE_DP3
c            action = INCREASE_DP1
c            action = DECREASE_DP1
          ELSE
            action = BEYOND_SCOPE
          ENDIF
c        ENDIF

      ELSE
        action = NULL
      ENDIF
c
c
c
      WRITE(PINOUT,'(1X,2I4,I5,3G12.4,I4,A)')
     .  region,ir,error,sp,sm,se,action,' <NOTE>'
      WRITE(0     ,'(1X,2I4,I5,3G12.4,I4,A)')
     .  region,ir,error,sp,sm,se,action,' <NOTE>'
c
c
c
c
c
      WRITE(PINOUT,*) 'osm_dp1-3 = ',osm_dp1(region,ir),
     .                               osm_dp2(region,ir),
     .                               osm_dp3(region,ir)
      WRITE(0     ,*) 'osm_dp1-3 = ',osm_dp1(region,ir),
     .                               osm_dp2(region,ir),
     .                               osm_dp3(region,ir)

c      IF (region.EQ.IKLO) THEN
c        WRITE(0,*) 'BOGUS ACTION'
c        action = INCREASE_DP3
c      ENDIF

      IF     (action.EQ.INCREASE_DP1) THEN
        osm_dp1(region,ir) = osm_dp1(region,ir) * 1.10

      ELSEIF (action.EQ.DECREASE_DP1) THEN
        osm_dp1(region,ir) = osm_dp1(region,ir) * 0.90

      ELSEIF (action.EQ.INCREASE_DP3) THEN
c...dp3
        WRITE(PINOUT,*) 'Evaluating power terms:'
        WRITE(0     ,*) 'Evaluating power terms (NO IKOPT!):'

c...TEMP!
c        WRITE(0,*)
c        WRITE(0,*) '***TEMP WHIPE 03***'
c        WRITE(0,*)
c
c        IF (region.EQ.IKLO) THEN
c          DO ik = 1, SymmetryPoint(ir)
c            osm_dp6(ik,ir) = 1.00
c          ENDDO
c        ELSEIF (region.EQ.IKHI) THEN
c          DO ik = SymmetryPoint(ir)+1, nks(ir)
c            osm_dp6(ik,ir) = 1.00
c          ENDDO
c        ENDIF

        WRITE(PINOUT,*) 'CALLING TEMPERATURE LEVELLING ROUTINE'

        i1     = 0
        status = 1
        DO WHILE (status.GT.0.AND.i1.LT.200)
          i1 = i1 + 1

c          WRITE(0,*) i1,region

          IF (status.NE.3) THEN
            CALL INITPLASMA(ir,ir,3)
            CALL SOL_PLASMA(ir,ir,3)
            CALL SOL       (ir,ir,3)
          ENDIF

c          IF (stopopt3.EQ.9.AND.osm_code(region,ir).EQ.1.AND.
c     .        region.EQ.IKLO) THEN
c            CALL ShapeDensityPeak(region,ir,status)
c          ELSE
            CALL FlattenTeProfile(region,ir,status)
c          ENDIF

c...bad 3 instead of IKOPT!
        ENDDO

        WRITE(PINOUT,*) '  ITERATION COUNT, REGION = ',i1,region

        IF (status.EQ.-1) result = -1
c        osm_dp4(ike,ir) = osm_dp4(ike,ir) + 0.1
c        osm_dp3(region,ir) = osm_dp3(region,ir) + 0.1
      ELSEIF (action.EQ.DECREASE_DP3) THEN
c...dp3
        WRITE(PINOUT,*) 'Evaluating power terms:'
        WRITE(0     ,*) 'Evaluating power terms (NO IKOPT!):'

c...TEMP!
c        WRITE(0,*)
c        WRITE(0,*) '***TEMP WHIPE 04***'
c        WRITE(0,*)
c
c        IF (region.EQ.IKLO) THEN
c          DO ik = 1, SymmetryPoint(ir)
c            osm_dp6(ik,ir) = 1.00
c          ENDDO
c        ELSEIF (region.EQ.IKHI) THEN
c          DO ik = SymmetryPoint(ir)+1, nks(ir)
c            osm_dp6(ik,ir) = 1.00
c          ENDDO
c        ENDIF

        WRITE(PINOUT,*) 'CALLING TEMPERATURE LEVELLING ROUTINE'

        i1     = 0
        status = 1

        DO WHILE (status.GT.0.AND.i1.LT.200)
          i1 = i1 + 1

c          WRITE(0,*) '--)',i1,region

          IF (status.NE.3) THEN
            CALL INITPLASMA(ir,ir,3)
            CALL SOL_PLASMA(ir,ir,3)
            CALL SOL       (ir,ir,3)
          ENDIF

c          IF (stopopt3.EQ.9.AND.osm_code(region,ir).EQ.1.AND.
c     .        region.EQ.IKLO) THEN
c
c...Shouldn't be here... needs separate routine to check problems with density peak,
c   when that is okay, go here...
c
c            CALL ShapeDensityPeak(region,ir,status)
c          ELSE
            CALL FlattenTeProfile(region,ir,status)
c          ENDIF
        ENDDO

        WRITE(PINOUT,*) '  ITERATION COUNT, REGION = ',i1,region

        IF (status.EQ.-1) result = -1

c        osm_dp4(ike,ir) = osm_dp4(ike,ir) - 0.1
c        osm_dp3(region,ir) = osm_dp3(region,ir) - 0.1
      ELSE
        result = 0
      ENDIF



      WRITE(PINOUT,*) 'osm_dp1-3 = ',osm_dp1(region,ir),
     .                               osm_dp2(region,ir),
     .                               ike,osm_dp4(ike,ir)
      WRITE(0     ,*) 'osm_dp1-3 = ',osm_dp1(region,ir),
     .                               osm_dp2(region,ir),
     .                               ike,osm_dp4(ike,ir)



      RETURN
99    STOP
      END
c
c ======================================================================
c
c
c
      SUBROUTINE AssignMockPowerExtent(region,ir)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INCLUDE 'solparams'
      INCLUDE 'solcommon'
      INCLUDE 'solswitch'

      INTEGER region,ir

      INTEGER GetModel,SymmetryPoint
      REAL    FindPeak

      INTEGER ik,ik1,ik2,ik3,iks,ike,ikm,npeak,ikp,ihold1,iflag
      REAL    qeint,qetot,s1,s2,rdum1,rdum2
      LOGICAL peak

c
c
c
      CALL SetBounds



      

c
c
c
      IF (GetModel(region,ir).NE.22.AND.
     .    GetModel(region,ir).NE.24) RETURN


      DO ik = 1, nks(ir)
        osm_dp4(ik,ir) = 1.0
c        osm_dp6(ik,ir) = 1.0
      ENDDO

c      osm_code(region,ir) = 1




c      osm_code(IKLO,19) = 2
c      osm_code(IKLO,20) = 2
c      osm_code(IKLO,21) = 2
c      osm_code(IKLO,22) = 2
c      osm_code(IKLO,23) = 2
c
c     Assign OSM_DP5:
c

      ihold1 = osm_mode
      osm_mode = 1

c      WRITE(0,*) 'MARK: ASSIGN 1',osm_mode

c      IF (.NOT.thesis) THEN
c        WRITE(0,*) 'NOT SEARCHING FOR IONISATION PEAKS'
c        DO ir = irsep, nrs
c          IF (idring(ir).EQ.-1) CYCLE
c          osm_dp5(IKLO,ir) = 0.25 * ksmaxs(ir)
c          osm_dp5(IKHI,ir) = 0.75 * ksmaxs(ir)
c        ENDDO
c        GOTO 50
c      ENDIF

      IF (GetModel(region,ir).EQ.24) THEN
c      IF (ir.GE.irsep.AND.ir.LT.irwall.AND.
c     .    GetModel(region,ir).EQ.24) THEN
        CALL SetBounds

        iks = ikbound(ir,region)
        ikm = iks
        ike = SymmetryPoint(ir)

        IF     (region.EQ.IKLO) THEN
          npeak = 0
          DO ik = iks, ike
            IF (FindPeak(ik,ir,pinqe).GT.0.0) npeak = npeak + 1
c            WRITE(0,*) 'FINDPEAK2CHECK ',
c     .                ik,ir,pinqe(ik,ir),FindPeak(ik,ir,pinqe)
          ENDDO

          IF (npeak.LT.1.OR.npeak.GT.2) THEN
            WRITE(0,*) 'WARNING! Unexpected peak count'
          ENDIF

c... Added the reference to GetModel because this was halting the
c    code for 990429019 when n-n collisions were included in the
c    PFZ void. Jul 14, 2000
c          IF (npeak.EQ.0) 
          IF (npeak.EQ.0.AND.GetModel(IKLO,ir).EQ.24) 
     .      CALL ER('AssignMockPowerExtent','No peak found',*99)

          ik = iks
          DO WHILE (ik.LT.ike.AND.FindPeak(ik,ir,pinqe).EQ.0.0)

c            WRITE(0,*) 'FINDPEAK CHECK ',
c     .                 ik,ir,pinqe(ik,ir),FindPeak(ik,ir,pinqe)

            ik = ik + 1
          ENDDO

          ikp = ik




          osm_dp5(IKLO,ir) = ksb(ik-1,ir)

c...new method
c          osm_code(IKLO,ir) = 2
c          DO ik = iks, ike
c            osm_dp6(ik,ir) = 1.0
c          ENDDO

c      WRITE(0,*) 'MARK: ASSIGN 2',osm_mode

          IF (osm_mode.GE.1) THEN
c            WRITE(0     ,90) 'Lower peak found  ',iks,ik,ike,ir,
c     .                       osm_code(region,ir)
            WRITE(PINOUT,90) 'Lower peak found  ',iks,ik,ike,ir,
     .                       osm_code(region,ir)
90          FORMAT(A,4I6,2X,I6)
          ENDIF


c...      MANUAL OVERRIDE TO MAINATIN THE SOLUTION! - Sep 23. 1999
c         Fix this you ass!
c
c         SEE ADDITIONAL OVERRIDE CODE BELOW
c
c         The ionisation is so washed out that the first peak cannot
c         be found, and the midplane peak is instead...
c
          iflag = 0
          IF (zs(ikp,ir).GT.-0.1.AND.ir.LT.irbreak) THEN

            iflag = 1

            ikp   = iks + 1
            ik    = ikp
            npeak = 2

            WRITE(0,     90) '        override  ',iks,ik,ike,ir,
     .                       osm_code(region,ir)
            WRITE(PINOUT,90) '        override  ',iks,ik,ike,ir,
     .                       osm_code(region,ir)
          ENDIF

          osmikp(IKLO,ir,1) = ikp

c
c         Look for a second peak (likely due to recycling from the inner wall):
c

c...need parameter
          IF (npeak.GT.1) THEN
c...take last peak
            ik = ike
            DO WHILE (FindPeak(ik,ir,pinqe).LT.1.00)
c            DO WHILE (FindPeak(ik,ir,pinqe).LT.0.80)            
              ik = ik - 1
c            WRITE(0,*) 'FINDPEAK2CHECK ',
c     .                ik,ir,pinqe(ik,ir),FindPeak(ik,ir,pinqe)
            ENDDO

            IF (ik.NE.ikp) THEN
              IF (osm_code(IKLO,ir).EQ.2.OR.iflexopt(4).EQ.21) THEN
c...temp
                osm_code(IKLO,ir) = 2

                IF (osm_mode.GE.1) THEN
                  WRITE(0,     90) '        update    ',iks,ik,ike,ir,
     .                             osm_code(region,ir)
                  WRITE(PINOUT,90) '        update    ',iks,ik,ike,ir,
     .                             osm_code(region,ir)
                ENDIF

                osmikp(IKLO,ir,2) = ik

                osm_dp5(IKLO,ir) = ksb(ik-1,ir)
c...clear OSM_DP6 in case last OSM_CODE for this ring was 1, not 2
                DO ik = iks, ike
                  osm_dp6(ik,ir) = 1.0
                ENDDO
              ELSE
c...            Need to figure out how to make a smooth transition when 
c               OSM_CODE changes.  Perhaps the midplane peak should not be 
c               allowed to dominate until it is much larger that the
c               target peak?  Perhaps the ring needs to be whiped in some 
c               friendly manner?  It is likely, however, that the local
c               multiplier will be keeping the target peak alive -- so
c               doing away with 
                WRITE(0     ,*) 'OSM_CODE CHANGE DISSALLOWED 2'
                WRITE(PINOUT,*) 'OSM_CODE CHANGE DISSALLOWED 2'              
                DO ik = iks, ike
                  rdum1 = FindPeak(ik,ir,pinqe)
                  WRITE(PINOUT,'(10X,2I4,1P,E10.2,0P,F8.3)')
     .              ik,ir,pinqe(ik,ir),rdum1
                ENDDO
              ENDIF
            ELSE
              IF (osm_code(IKLO,ir).EQ.1) THEN
                osm_code(IKLO,ir) = 1            
              ELSE
                WRITE(0     ,*) 'OSM_CODE CHANGE DISSALLOWED 1'
                WRITE(PINOUT,*) 'OSM_CODE CHANGE DISSALLOWED 1'              
                DO ik = iks, ike
                  rdum1 = FindPeak(ik,ir,pinqe)
                  WRITE(PINOUT,'(10X,2I4,1P,E10.2,0P,F8.3)')  
     .              ik,ir,pinqe(ik,ir),rdum1
                ENDDO
              ENDIF
            ENDIF
          
          ELSE
c.. This is no longer set above becasue when deciding whether or not to
c   switch from OSM_CODE=1 to OSM_CODE=2, the last OSM_CODE value is used.
c   - Jan 30, 2000
            osm_code(IKLO,ir) = 1
          
          ENDIF


c...      MANUAL OVERRIDE!

          IF (iflag.EQ.1.AND.osm_code(IKLO,ir).NE.2) THEN
            CALL ER('AMPE','MANUAL OVERRIDE FAILED',*99)
          ENDIF



        ELSEIF (region.EQ.IKHI) THEN

          osm_code(IKHI,ir) = 1            

          ik = iks
          DO WHILE (ik.GT.ike.AND.FindPeak(ik,ir,pinqe).EQ.0.0)
            ik = ik - 1
          ENDDO

          IF (ik.EQ.ike.AND.GetModel(IKHI,ir).Eq.24) 
     .      CALL ER('AssignMockPowerExtent','No peak found',*99)

          osm_dp5(IKHI,ir) = ksb(ik,ir)

          IF (osm_mode.GE.1) THEN
c            WRITE(0     ,90) 'Higher peak found ',iks,ik,ike,ir,
c     .                       osm_code(region,ir)
            WRITE(PINOUT,90) 'Higher peak found ',iks,ik,ike,ir,
     .                       osm_code(region,ir)
          ENDIF


c...No check is presently done for a second peak on an outer leg



        ELSE
          CALL ER('AssignMockPowerExtent','Invalid region',*99)
        ENDIF
      ELSE
        IF (region.EQ.IKLO) THEN
          osm_dp5(IKLO,ir) = ksb(0,ir)
        ELSE
          osm_dp5(IKHI,ir) = ksmaxs(ir)
        ENDIF
      ENDIF

c      WRITE(PINOUT,*)
c     .  'REGION,IR OSM_CODE = ',region,ir,osm_code(region,ir)
c      WRITE(0     ,*)
c     .  'REGION,IR OSM_CODE = ',region,ir,osm_code(region,ir)



c      STOP

c
c     Remove cf power for <DP5 region on code=2:
c
c      IF (osm_code(IKLO,ir).EQ.2) THEN
c        iks = ikbound(ir,IKLO)
c        ike = SymmetryPoint(ir)
c
c        DO ik = iks, ike
c          IF (kss(ik,ir).LT.osm_dp5(IKLO,ir)) osm_dp4(ik,ir) = 0.0
c        ENDDO
c      ENDIF
c
c      IF (osm_code(IKHI,ir).EQ.2) THEN
c        iks = SymmetryPoint(ir) + 1
c        ike = ikbound(ir,IKHI)
c
c        DO ik = iks, ike
c          IF (kss(ik,ir).GT.osm_dp5(IKHI,ir)) osm_dp4(ik,ir) = 0.0
c        ENDDO
c      ENDIF


      IF (switch(SWPOW3).EQ.1.0.OR.switch(SWPOW3).EQ.6.0.OR.
     .    switch(SWPOW3).EQ.7.0) THEN

        osm_dp1(IKLO,ir) = 0.0
        osm_dp1(IKHI,ir) = 0.0

        CALL SetBounds

        ik1 = ikbound(ir,IKLO)
        ik2 = SymmetryPoint(ir)

        CALL CalcIntegral4(pinqe,ik1,ik2,ir,rdum1,2)

        ik3 = ik2
        DO ik = ik2, ik1, -1
          CALL CalcIntegral4(pinqe,ik1,ik,ir,rdum2,2)
          IF (rdum2.GT.0.95*rdum1) ik3 = ik
        ENDDO

        IF (ik3.EQ.1) THEN
          s1 = ksb(0  ,ir)
        ELSE
          s1 = kss(ik3,ir)
        ENDIF

        osm_dp2(IKLO,ir) = s1  / ksmaxs(ir) * 1.2

        ik1 = ikbound(ir,IKHI)
        ik2 = SymmetryPoint(ir) + 1

        CALL CalcIntegral4(pinqe,ik2,ik1,ir,rdum1,2)

        ik3 = ik2
        DO ik = ik2, ik1
          CALL CalcIntegral4(pinqe,ik,ik1,ir,rdum2,2)
          IF (rdum2.GT.0.95*rdum1) ik3 = ik
        ENDDO

        IF (ik3.EQ.nks(ir)) THEN
          s1 = ksb(nks(ir),ir)
        ELSE
          s1 = kss(ik3,ir)
        ENDIF

        osm_dp2(IKHI,ir) = (ksmaxs(ir) - s1) / ksmaxs(ir) * 1.2


c        WRITE(0,*) 'OSMDP2 ',osm_dp2(IKLO,ir),osm_dp2(IKHI,ir)


c        osm_dp2(IKLO,ir) = HI
c        osm_dp2(IKHI,ir) = HI



        osm_dp3(IKLO,ir) = 1.0
        osm_dp3(IKHI,ir) = 1.0

        RETURN
      ENDIF
c
c
c
c
c
      IF     (region.EQ.IKLO) THEN

c        WRITE(0     ,*) 'LENGTH = Inner'
c        WRITE(PINOUT,*) 'LENGTH = Inner'

        osm_dp3(IKLO,ir) = 1.0

        ik1 = ikbound(ir,IKLO)
        ik3 = SymmetryPoint(ir)
c
c       Only want to sample PINQE near target/prescription boundary:
c
        ik  = ik1
        ik2 = ik3
        DO WHILE (ik.LT.ik3.AND.ik2.EQ.ik3)
          ik = ik + 1

          IF (ABS(pinqe(ik,ir)).LT.0.05*ABS(pinqe(ik1,ir))) ik2 = ik
c          IF (pinqe(ik,ir).EQ.0.0) ik2 = ik
        ENDDO

        ik2 = MIN(MAX(ik2,ik1+5),ik3)

        CALL CalcIntegral4(pinqe,ik1,ik2,ir,rdum1,2)
        CALL CalcIntegral4(pinqe,ik1,ik3,ir,rdum2,2)

c        WRITE(0     ,*) 'LENGTH = ',ik1,ik2,ik3,rdum1,rdum2
c        WRITE(PINOUT,*) 'LENGTH = ',ik1,ik2,ik3,rdum1,rdum2

        IF (rdum1.LT.0.01*rdum2)
     .    CALL ER('AssignMockPowerExtent','Suspicious integral',*99)
c
c       Find PINQE peak:
c
        ik3 = ik2
        ik2 = ik1
        DO ik = ik1+1, ik3
          IF (ABS(pinqe(ik,ir)).GT.ABS(pinqe(ik2,ir))) ik2 = ik
        ENDDO

        ik1 = ik2
c
c
c
        CALL CalcIntegral4(pinqe,ik1,ik3,ir,qetot,2)

c        WRITE(PINOUT,*) 'LENGTH = ',qetot,rdum1,rdum2,ik1,ik3
c        WRITE(0     ,*) 'LENGTH = ',qetot,rdum1,rdum2,ik1,ik3

        ik2 = ik3

        DO ik = ik3, ik1, -1
          CALL CalcIntegral4(pinqe,ik1,ik,ir,qeint,2)

          IF (qeint/qetot.GE.0.90) ik2 = ik

c          WRITE(PINOUT,*) 'LENGTH = ',region,ir,qeint,qeint/qetot,ik2
c          WRITE(0     ,*) 'LENGTH = ',region,ir,qeint,qeint/qetot,ik2
        ENDDO

        IF (ik2.EQ.ik1) ik2 = ik1 + 3

        IF (ik2.EQ.ik1)
     .    CALL ER('AssignMockPowerExtent','No inner extent',*99)


        IF (ik1.EQ.1) THEN
          s1 = ksb(0  ,ir)
        ELSE
          s1 = kss(ik1,ir)
        ENDIF

        s2 = kss(ik2,ir)

        IF     (switch(SWPOW3).EQ.4.0) THEN
          osm_dp1(IKLO,ir) = (s2 - s1) / ksmaxs(ir)
          osm_dp2(IKLO,ir) =       s1  / ksmaxs(ir)

        ELSEIF (switch(SWPOW3).EQ.5.0) THEN
c          osm_dp1(IKLO,ir) = (s2 - s1) / ksmaxs(ir) * 0.5

          osm_dp1(IKLO,ir) = SQRT(-0.5 / LOG(0.1) *
     .                            ((s2 - s1) / ksmaxs(ir))**2.0)
          osm_dp2(IKLO,ir) = s1 / ksmaxs(ir)

        ELSE
          CALL ER('AssignMockPowerExtent','Invalid SWPOW3 option',*99)
        ENDIF

c        WRITE(PINOUT,*) 'LENGTH = ',s1,s2,ksmaxs(ir),
c     .     osm_dp1(IKLO,ir),osm_dp2(IKLO,ir),osm_dp3(IKLO,ir),
c     .     ik1,ik2,ik3,ikbound(ir,IKLO),SymmetryPoint(ir)
c        WRITE(0     ,*) 'LENGTH = ',s1,s2,ksmaxs(ir),
c     .     osm_dp1(IKLO,ir),osm_dp2(IKLO,ir),osm_dp3(IKLO,ir),
c     .     ik1,ik2,ik3,ikbound(ir,IKLO),SymmetryPoint(ir)

      ELSEIF (region.EQ.IKHI) THEN

c        WRITE(0     ,*) 'LENTH = Outer'
c        WRITE(PINOUT,*) 'LENTH = Outer'


        osm_dp3(IKHI,ir) = 1.0

        ik1 = ikbound(ir,IKHI)
        ik3 = SymmetryPoint(ir)
c
c       Only want to sample PINQE near target/prescription boundary:
c
        ik  = ik1
        ik2 = ik3
        DO WHILE (ik.GT.ik3.AND.ik2.EQ.ik3)
          ik = ik - 1

          IF (ABS(pinqe(ik,ir)).LT.0.05*ABS(pinqe(ik1,ir))) ik2 = ik
c          IF (pinqe(ik,ir).EQ.0.0) ik2 = ik
        ENDDO

        ik2 = MAX(MIN(ik2,ik1-5),ik3)

        CALL CalcIntegral4(pinqe,ik2,ik1,ir,rdum1,2)
        CALL CalcIntegral4(pinqe,ik3,ik1,ir,rdum2,2)

c        WRITE(0     ,*) 'LENGTH = ',ik1,ik2,ik3,rdum1,rdum2
c        WRITE(PINOUT,*) 'LENGTH = ',ik1,ik2,ik3,rdum1,rdum2

        IF (rdum1.LT.0.01*rdum2)
     .    CALL ER('AssignMockPowerExtent','Suspicious integral',*99)
c
c       Find PINQE peak:
c
        ik3 = ik2
        ik2 = ik1
        DO ik = ik1-1, ik3, -1
          IF (ABS(pinqe(ik,ir)).GT.ABS(pinqe(ik2,ir))) ik2 = ik
        ENDDO

        ik1 = ik2
c
c
c
        CALL CalcIntegral4(pinqe,ik3,ik1,ir,qetot,2)

c        WRITE(PINOUT,*) 'LENGTH = ',qetot,rdum1,rdum2,ik1,ik3
c        WRITE(0     ,*) 'LENGTH = ',qetot,rdum1,rdum2,ik1,ik3

        ik2 = ik3
        DO ik = ik3, ik1
          CALL CalcIntegral4(pinqe,ik,ik1,ir,qeint,2)

          IF (qeint/qetot.GE.0.90) ik2 = ik

c          WRITE(PINOUT,*) 'LENGTH = ',region,ir,qeint,qeint/qetot,ik2
c          WRITE(0     ,*) 'LENGTH = ',region,ir,qeint,qeint/qetot,ik2
        ENDDO

        IF (ik2.EQ.ik1) ik2 = ik1 - 3

        IF (ik2.EQ.ik1)
     .    CALL ER('AssignMockPowerExtent','No outer extent',*99)

        IF (ik1.EQ.nks(ir)) THEN
          s1 = ksb(nks(ir),ir)
        ELSE
          s1 = kss(ik1,ir)
        ENDIF

        s2 = kss(ik2,ir)

        IF     (switch(SWPOW3).EQ.4.0) THEN
          osm_dp1(IKHI,ir) = (s1         - s2) / ksmaxs(ir)
          osm_dp2(IKHI,ir) = (ksmaxs(ir) - s1) / ksmaxs(ir)

        ELSEIF (switch(SWPOW3).EQ.5.0) THEN
c          osm_dp1(IKHI,ir) = (s1 - s2) / ksmaxs(ir) * 0.5

          osm_dp1(IKHI,ir) = SQRT(-0.5 / LOG(0.1) *
     .                            ((s1 - s2) / ksmaxs(ir))**2.0)
          osm_dp2(IKHI,ir) = (ksmaxs(ir) - s1) / ksmaxs(ir)

        ELSE
          CALL ER('AssignMockPowerExtent','Invalid SWPOW3 option',*99)
        ENDIF


c        WRITE(PINOUT,*) 'LENGTH = ',s1,s2,ksmaxs(ir),
c     .     osm_dp1(IKHI,ir),osm_dp2(IKHI,ir),osm_dp3(IKHI,ir),
c     .     ik1,ik2,ik3,ikbound(ir,IKHI),SymmetryPoint(ir)
c        WRITE(0     ,*) 'LENGTH = ',s1,s2,ksmaxs(ir),
c     .     osm_dp1(IKHI,ir),osm_dp2(IKHI,ir),osm_dp3(IKHI,ir),
c     .     ik1,ik2,ik3,ikbound(ir,IKHI),SymmetryPoint(ir)
      ELSE
        CALL ER('AssignMockPowerExtent','Invalid REGION value',*99)
      ENDIF

      osm_mode = ihold1

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
      REAL FUNCTION IdentifyFlowReversal(region,ring)
      IMPLICIT none

      INTEGER region,ring

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER GetModel,SymmetryPoint

      INTEGER ik1,ik2,ir
      REAL    ionsrc,recsrc,netsrc,tarflx

      ir = ring

      IF     (region.EQ.IKLO) THEN
        ik1 = ikbound      (ir,IKLO)
        ik2 = SymmetryPoint(ir)

        CALL CalcIntegral3(pinion,ik1,ik2,ir,ionsrc,2)
        CALL CalcIntegral3(pinrec,ik1,ik2,ir,recsrc,2)

        netsrc = ionsrc - recsrc

        IF     (GetModel(IKLO,ir).EQ.22) THEN
          tarflx = ABS(knds(idds(ir,2)) * kvds(idds(ir,2)))
        ELSEIF (GetModel(IKLO,ir).EQ.24) THEN
          tarflx = ABS(knbs(ik1,ir) * kvhs(ik1,ir))
        ELSE
          CALL ER('IdentifyFlowReversal','Invalid SOL option',*99)
        ENDIF

      ELSEIF (region.EQ.IKHI) THEN
        ik1 = SymmetryPoint(ir)
        ik2 = ikbound      (ir,IKHI)

        CALL CalcIntegral3(pinion,ik1,ik2,ir,ionsrc,2)
        CALL CalcIntegral3(pinrec,ik1,ik2,ir,recsrc,2)

        netsrc = ionsrc - recsrc

        IF     (GetModel(IKHI,ir).EQ.22) THEN
          tarflx = ABS(knds(idds(ir,1)) * kvds(idds(ir,1)))
        ELSEIF (GetModel(IKHI,ir).EQ.24) THEN
          tarflx = ABS(knbs(ik2,ir) * kvhs(ik2,ir))
        ELSE
          CALL ER('IdentifyFlowReversal','Invalid SOL option',*99)
        ENDIF

      ELSE
        CALL ER('IdentifyFlowReversal','Invalid REGION',*99)
      ENDIF

c      WRITE(PINOUT,*) ' FLOW : ',region,ir,tarflx,netsrc,netsrc/tarflx
c      WRITE(0     ,*) ' FLOW : ',region,ir,tarflx,netsrc,netsrc/tarflx

c      IF (netsrc/tarflx.GT.1.0) THEN
c        IdentifyFlowReversal = 1
c      ELSE
c        IdentifyFlowReversal = 0
c      ENDIF

      IF (netsrc/tarflx.GT.1.1) THEN
        IdentifyFlowReversal = netsrc / tarflx
      ELSE
        IdentifyFlowReversal = 0.0
      ENDIF

      RETURN
99    WRITE(0,*) '  REGION RING SOL = ',region,ir,GetModel(region,ir)
      STOP
      END
c
c ======================================================================
c
c
c
      SUBROUTINE FlattenTeProfile(region,ring,status)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INCLUDE 'solparams'
      INCLUDE 'solswitch'


      COMMON /CFSCOM/ cfs_mage,cfs_magi,cfs_sume,cfs_sumi
      REAL            cfs_mage(MAXNRS),cfs_magi(MAXNRS),
     .                cfs_sume(MAXNRS),cfs_sumi(MAXNRS)

      LOGICAL           modify(MAXNKS)

      COMMON /STATCOM/ osm_ncnt ,osm_nerr,osm_temod,osm_relfr,osm_cadj,
     .                 osm_tcon
      INTEGER          osm_ncnt (MAXNRS)    ,osm_nerr (2,MAXNRS),
     .                 osm_cadj (5,2,MAXNRS),osm_tcon (MAXNRS)
      REAL             osm_temod(MAXNRS)    ,osm_relfr(MAXNRS)



      REAL rdum1,rdum2,rdum3,rdum4,rdum5,rdum6

      LOGICAL output

      REAL CalcDist

      INTEGER region,ring,status

      INTEGER SymmetryPoint

      INTEGER ir,ik,iks,ikb,ike,ikm,ikl,ikf,i1,i2,ikr,
     .        cmod2(MAXNKS),stepl,iterl,ringl,regionl,mode,model,modec,
     .        code,tmpcode
      REAL    qeint,qetot,pei1,pei2,pei3,pei4,
     .        te1,te2,teadj,temod
      REAL    cdetmp(0:MAXNKS)


      DATA stepl,iterl /-1, -1/

      SAVE


      output = .FALSE.
c      IF (stopopt3.EQ.9.OR.stopopt3.EQ.12.OR.stopopt3.EQ.13)
c     .  output = .TRUE.
c      IF (ring.EQ.24) output = .TRUE.
c      output = .TRUE.


      ir     = ring

      IF (rel_iter.NE.iterl.OR.rel_step.NE.stepl) THEN
        regionl = -1
        ringl   = -1

        iterl = rel_iter
        stepl = rel_step
      ENDIF

      IF (ring.NE.ringl) THEN
        CALL IZero(cmod2,MAXNKS)

        regionl = -1

        ringl = ring
      ENDIF

      IF (region.NE.regionl) THEN
c...add qualifiers!
        IF (osm_peimul(region,ir).EQ.0.0) osm_peimul(region,ir) = 0.1

        teadj = 0.0
        temod = 0.0
        ikr   =  0
        ikl   = -1

        model = -1
        modec =  0



c
c       Assigning CODE:
c
        IF     (osm_code(region,ring).EQ.1) THEN
          code = 1
        ELSEIF (osm_code(region,ring).EQ.2) THEN
          IF (region.EQ.IKLO) THEN
            ik  = 1
            ike = SymmetryPoint(ir)
            DO WHILE(ktebs(ik,ir).LE.ktebs(ik+1,ir).AND.ik.LT.ike)
              ik = ik + 1
            ENDDO
c...should be more qualifiers here
            IF (switch(SWPEI).EQ.5.0) THEN
              IF (osmpei(ik,ir).EQ.0.0) THEN
                code = 2
                WRITE(0,*) 'INCREASING IONISATION SOURCE  REG. =',region
              ELSE
                code = 1
                WRITE(0,*) 'STANDARD MOCK POWER ADJ.      REG. =',region
              ENDIF
              WRITE(PINOUT,*) 'OSM_CODE=2 SELECTION  CODE = ',code
            ELSEIF (osmmock.GT.0) THEN
              IF (osmpmk(ik,ir).EQ.0.0) THEN
                code = 2
                WRITE(0,*) 'INCREASING IONISATION SOURCE  REG. =',region
              ELSE
                code = 1
                WRITE(0,*) 'STANDARD MOCK POWER ADJ.      REG. =',region
              ENDIF
              WRITE(PINOUT,*) 'OSM_CODE=2 SELECTION  CODE = ',code
            ELSE
              CALL ER('FlattenTeProfile','Invalid mock power opt',*99)
            ENDIF
          ELSE
            code = 2
          ENDIF
        ELSE
          CALL ER('FlattenTeProfile','Invalid OSM_CODE value',*99)
        ENDIF



        regionl = region
      ENDIF

      status = 0


      osm_ncnt(ir)          = osm_ncnt(ir) + 1
      osm_cadj(2,region,ir) = code

c
c
c
c
c
      IF     (region.EQ.IKLO) THEN
c
c       Find index range:
c
        IF (ikl.EQ.-1) THEN
          iks = ikbound(ir,IKLO)
          ikm = SymmetryPoint(ir)
          ike = 0

          qeint = 0.0
          qetot = 0.0

          CALL CalcIntegral4(pinqe,iks,ikm,ir,qetot,2)

          ike = iks+1

          DO ik = iks, ikm
            modify(ik) = .FALSE.
          ENDDO

c...should be code=1, not osm_code=1? *****
          IF     (osm_code(region,ir).EQ.1) THEN
            DO ik = iks+1, ikm
              IF (qeint.LT.0.95*qetot) THEN
                CALL CalcIntegral4(pinqe,iks,ik,ir,qeint,2)
                ike = ik
              ENDIF
            ENDDO

          ELSEIF (osm_code(region,ir).EQ.2) THEN
            DO ik = iks+1, ikm
              IF (kss(ik,ir).LT.osm_dp5(region,ir)) ike = ik
            ENDDO

            CALL CalcIntegral4(pinion,iks,ike,ir,rdum5,2)

          ELSE
            CALL ER('FlattenTeProfile','Unknown low index OSM_CODE',*99)
          ENDIF

          ikl = iks + 1

          ikb = -1
        ENDIF

c
c
c
c...CHANGES
        DO ik = 1, ikm
          IF (ik.EQ.1) THEN
            rdum1 = ktebs(ik,ir)
            rdum2 = (ktebs(ik,ir) - kteds(idds(ir,2))) /
     .              (kss  (ik,ir) - ksb  (0,ir)      )

            cdetmp(ik) = -2000.0 * rdum1**(5.0/2.0) * rdum2
          ELSE
            rdum1 = ktebs(ik,ir)
            rdum2 = (ktebs(ik,ir) - ktebs(ik-1,ir)) /
     .              (kss  (ik,ir) - kss  (ik-1,ir))

            cdetmp(ik) = -2000.0 * rdum1**(5.0/2.0) * rdum2
          ENDIF
        ENDDO
c        DO ik = 2, ikm
c          rdum1 = ktebs(ik,ir)
c          rdum2 = (ktebs(ik,ir) - ktebs(ik-1,ir)) /
c     .            (kss  (ik,ir) - kss  (ik-1,ir))
c
c          cdetmp(ik) = -2000.0 * rdum1**(5.0/2.0) * rdum2
c        ENDDO

        IF (ikb.EQ.-1) THEN
          DO ik = iks, ike
            IF     (code.EQ.1) THEN
              IF (output)
     .          WRITE(PINOUT,'(10X,4I5,2F12.5,1P,5E10.2,0P)')
     .            ik,iks,ir,code,
     .            ktebs(ik,ir),osm_dp6(ik,ir),
     .            -osmpei(ik,ir),-osmpmk(ik,ir),
     .            pinqe(ik,ir),osmcfe(ik,ir),cdetmp(ik)

            ELSEIF (code.EQ.2) THEN
              IF (output)
     .          WRITE(PINOUT,'(10X,4I5,2F12.5,1P,3E10.2,0P,
     .                         F10.5,1P,2E10.2,0P)')
     .            ik,iks,ir,
     .            code,ktebs(ik,ir),osm_dp4(ik,ir),
     .            -osmpei(ik,ir),-osmpmk(ik,ir),
     .            pinqe(ik,ir),osm_dp6(ik,ir),osmcfe(ik,ir),cdetmp(ik)

            ENDIF
          ENDDO

          IF     (code.EQ.1) THEN
c...CHANGES!
            ikb = iks
c            ikb = iks + 1

            DO WHILE (ikb.LT.ike.AND.cdetmp(ikb).LT.0.0)
              ikb = ikb + 1
            ENDDO

c...mark1...should this be done for rings that are not detached?
            IF (ikb.GT.iks+1.AND.
     .          ktebs(iks+1,ir)-ktebs(iks,ir).LT.osm_testep)
     .        ikb = iks + 1

          ELSEIF (code.EQ.2) THEN
            ikb = iks + 1

            DO WHILE (ikb.LT.ike.AND.cdetmp(ikb).LT.0.0)
              ikb = ikb + 1
            ENDDO

            DO i1 = ikb+1, ikm
              modify(i1) = .FALSE.
            ENDDO

c            ikb = ike
          ENDIF

        ENDIF



        mode = 0

        IF     (code.EQ.1) THEN
          IF (cdetmp(ikb).GE.0.0) mode = 1
        ELSEIF (code.EQ.2) THEN
          DO ik = iks+1, ikb
            IF (cdetmp(ik).GE.0.0) mode = 1
          ENDDO
        ENDIF

        IF (cfs_sume(ir).LE.0.0.OR.cfs_mage(ir).LE.0.0) mode = 1
c...mark1
        IF (code.EQ.1.AND.ikb.EQ.iks+1.AND.
     .      ktebs(iks+1,ir)-ktebs(iks,ir).LT.osm_testep) mode = 1
c
c
c
        IF     (mode.EQ.1.AND.model.EQ.-1) THEN
          IF     (code.EQ.1) THEN
            teadj = -SIGN(1.0,osm_peimul(IKLO,ir))
            teadj = -1.0
          ELSEIF (code.EQ.2) THEN
            teadj =  0.1
          ENDIF
          modec = 0
          model = 1
        ELSEIF (mode.EQ.1.AND.model.EQ.2) THEN
          model = 1
          modec = modec + 1
          teadj = -0.5 * teadj
        ELSEIF (mode.EQ.0.AND.model.EQ.1.AND.modec.LT.3) THEN
          model = 2
          modec = modec + 1
          teadj = -0.5 * teadj
        ELSEIF (mode.EQ.0.AND.modec.GE.3) THEN
          model = -1
        ELSEIF (modec.EQ.0) THEN
          teadj = teadj * 2.0
          IF (output) WRITE(PINOUT,*) '    INCREASING TEADJ',teadj
        ENDIF



        ik = ikb

        IF (output)
     .    WRITE(PINOUT,*) 'l -> ',iks,ike,ikm,ikr,osm_peimul(IKLO,ir),
     .      cfs_mage(ir),ikl,temod,mode,model,modec,teadj,
     .      code,ikb

        IF (code.EQ.1) THEN
          IF (output.AND.ik.GT.1)
     .      WRITE(PINOUT,'(3I5,2F9.3,1P,7E12.4,0P,I4)')
     .        ik-1,ir,region,
     .        ktebs(ik-1,ir),osm_dp6(ik-1,ir),
     .        -osmpei(ik-1,ir) * (ksb(ik-1,ir) - ksb(ik-2,ir)),
     .        -osmpmk(ik-1,ir) * (ksb(ik-1,ir) - ksb(ik-2,ir)),
     .         osmcfe(ik-1,ir) * (ksb(ik-1,ir) - ksb(ik-2,ir)),
     .         pinqe (ik-1,ir) * (ksb(ik-1,ir) - ksb(ik-2,ir)),
     .         cdetmp(ik-1),-1.0,-1.0,0
        ELSEIF (code.EQ.2) THEN
          IF (output)
     .      WRITE(PINOUT,'(3I5,2F9.3,1P,7E12.4,0P,I4)')
     .        ik-1,ir,region,
     .        ktebs(ik-1,ir),osm_dp4(ik-1,ir),
     .        -osmpei(ik-1,ir) * (ksb(ik-1,ir) - ksb(ik-2,ir)),
     .        -osmpmk(ik-1,ir) * (ksb(ik-1,ir) - ksb(ik-2,ir)),
     .         osmcfe(ik-1,ir) * (ksb(ik-1,ir) - ksb(ik-2,ir)),
     .         pinqe (ik-1,ir) * (ksb(ik-1,ir) - ksb(ik-2,ir)),
     .         cdetmp(ik-1),-1.0,-1.0,0
        ENDIF

        DO WHILE (ik.LE.ikb.AND.status.EQ.0)

          IF     (1.EQ.1.AND.
     .            (cfs_sume(ir).LE.0.0.OR.cfs_mage(ir).LE.0.0)) THEN

            WRITE(0     ,*) 'BUMMER! LOW REGION FAILED'
            WRITE(PINOUT,*) 'BUMMER! LOW REGION FAILED'
            status = -1

          ELSEIF (model.GT.0) THEN

            IF (code.EQ.1) THEN
              ikr = ik

              DO i2 = ik-1, iks, -1
                IF (osm_dp6(i2,ir).EQ.osm_dp6(ik,ir)) ikr = i2
              ENDDO
c              WRITE(0,*) '--)',ik,ikr,iks,teadj,ikb,cdetmp(ikb),mode
            ELSE
              ikr = ik

c...new method
c              DO i2 = ik-1, iks, -1
c                IF (.NOT.modify(i2)) ikr = i2
c              ENDDO
c              WRITE(0,*) ik,ikr,iks,teadj,ikb,cdetmp(ikb),mode
              ikr = iks
            ENDIF

            DO i2 = ikr, ik

              IF     (code.EQ.1) THEN

                IF (osm_peimul(IKLO,ir).GT.0.0) THEN
                  osm_dp6(i2,ir) = osm_dp6(i2,ir) + teadj
                ELSE
                  osm_dp6(i2,ir) = osm_dp6(i2,ir) - teadj
                ENDIF

                IF (output)
     .            WRITE(PINOUT,'(A,3I4,3F10.4,1P,2E10.2,0P)')
     .              '             One ',
     .              ik,ikr,i2,ktebs(i2,ir),
     .              osm_dp6(i2,ir),teadj,-osmpei(i2,ir),-osmpmk(i2,ir)

              ELSEIF (code.EQ.2) THEN

                IF (i2.LT.ike) THEN
                  pinion(i2,ir) = pinion(i2,ir) * (1.0 + teadj)
                  pinqe (i2,ir) = pinqe (i2,ir) * (1.0 + teadj)
                ENDIF
c                osm_dp4(i2,ir) = osm_dp4(i2,ir) + teadj

                IF (output)
     .            WRITE(PINOUT,'(A,3I4,3F10.4,1P,3E10.2,0P,L2)')
     .              '             Two ',
     .              ik,ikr,i2,ktebs(i2,ir),
     .              osm_dp4(i2,ir),teadj,osmcfe(i2,ir),
     .              pinion(i2,ir),pinqe(i2,ir),modify(i2)

              ENDIF

            ENDDO

            status = 1

          ELSEIF (ik.GT.iks+1.AND.osmcde(ik,ir).EQ.osmcde(ik-1,ir)) THEN

            status = 2

          ELSEIF (ik.GT.ikl) THEN

            status = 3

          ENDIF

c...new
          IF     (code.EQ.1) THEN
            IF (output)
     .        WRITE(PINOUT,'(3I5,2F9.3,1P,7E12.4,0P,I4)')
     .          ik,ir,region,ktebs(ik,ir),osm_dp6(ik,ir),
     .         -osmpei(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .         -osmpmk(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .          osmcfe(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .          pinqe (ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .          cdetmp(ik),-1.0,-1.0,cmod2(ik)

          ELSEIF (code.EQ.2) THEN
            IF (output)
     .        WRITE(PINOUT,'(3I5,2F9.3,1P,7E12.4,0P,I4)')
     .          ik,ir,region,ktebs(ik,ir),osm_dp4(ik,ir),
     .         -osmpei(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .         -osmpmk(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .          osmcfe(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .          pinqe (ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .          cdetmp(ik),-1.0,-1.0,cmod2(ik)

          ENDIF

          ik = ik + 1
        ENDDO
c
c
c

        IF (status.EQ.0.AND.ikb.LT.ike) status = 3

        IF     (status.EQ.0.OR.status.EQ.-1) THEN

c          tmpcode = tmpcode + 1
c          WRITE(0,*) 'INCREASING TMPCODE ',tmpcode

          ikl = -1
          ikr = 0
          teadj = 0.0
          temod = 0.0

          DO ik = iks, ike
c...new
            IF     (code.EQ.1) THEN
              IF (output)
     .          WRITE(PINOUT,'(10X,2I5,2F12.5)')
     .            ik,ir,ktebs(ik,ir),osm_dp6(ik,ir)

            ELSEIF (code.EQ.2) THEN
              IF (output)
     .          WRITE(PINOUT,'(10X,2I5,2F12.5)')
     .            ik,ir,ktebs(ik,ir),osm_dp4(ik,ir)

            ENDIF
          ENDDO

          IF (code.EQ.2) THEN
            CALL CalcIntegral4(pinion,iks,ike,ir,rdum6,2)

            WRITE(PINOUT,*) 'DELTA ION: ',
     .                      ir,rdum5,rdum6,rdum6/(rdum5+1.0)
          ENDIF

c          WRITE(PINOUT,*) 'END OF FLATTENING: '
c          osm_mode = 2
c          CALL INITPLASMA(ir,ir,3)
c          CALL SOL_PLASMA(ir,ir,3)
c          CALL SOL       (ir,ir,3)
c          osm_mode = 1
c
c          CALL SaveSolution

          model = -1
        ELSEIF (status.EQ.1) THEN
          ikl = ik - 1
          cmod2(ikl) = 1
        ELSEIF (status.EQ.3) THEN
          ikl = ik + 1
          ikr = 0

          DO i1 = iks, ikb
            modify(i1) = .TRUE.
          ENDDO

          ikb = -1
          teadj = 0.0
          temod = 0.0
        ENDIF

c...new
        IF     (code.EQ.1) THEN
          IF (output)
     .      WRITE(PINOUT,'(3I5,2F9.3,1P,7E12.4,0P,2I4)')
     .        ik,ir,region,ktebs(ik,ir),osm_dp6(ik,ir),
     .       -osmpei(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .       -osmpmk(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .        osmcfe(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .        pinqe (ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .        cdetmp(ik),-1.0,-1.0,status,ikl

        ELSEIF (code.EQ.2) THEN
          IF (output)
     .      WRITE(PINOUT,'(3I5,2F9.3,1P,6E12.4,0P,2I4)')
     .        ik,ir,region,ktebs(ik,ir),osm_dp4(ik,ir),
     .       -osmpei(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .       -osmpmk(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .        osmcfe(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .        pinqe (ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .        cdetmp(ik),-1.0,-1.0,status,ikl

        ENDIF

      ELSEIF (region.EQ.IKHI) THEN

c        CALL ER('FlattenTeProfile','Sorry, no go for region=IKHI',*99)








c
c       Find index range:
c
        IF (ikl.EQ.-1) THEN
          iks = ikbound(ir,IKHI)
          ikm = SymmetryPoint(ir) + 1
          ike = 0

          qeint = 0.0
          qetot = 0.0

          CALL CalcIntegral4(pinqe,ikm,iks,ir,qetot,2)

          ike = iks - 1

          IF     (osm_code(region,ir).EQ.1) THEN
            DO ik = iks-1, ikm, -1
              IF (qeint.LT.0.95*qetot) THEN
                CALL CalcIntegral4(pinqe,ik,iks,ir,qeint,2)
                ike = ik
              ENDIF
            ENDDO

          ELSEIF (osm_code(region,ir).EQ.2) THEN
            DO ik = iks-1, ikm, -1
              IF (kss(ik,ir).GT.osm_dp5(region,ir)) ike = ik
            ENDDO

          ELSE
            CALL ER('FlattenTeProfile','Unknown high index OSM_CODE',
     .              *99)
          ENDIF

          ikl = iks - 1

          ikb = -1
        ENDIF

c
c
c
        DO ik = ikm, nks(ir)-1
          rdum1 = ktebs(ik,ir)
          rdum2 = (ktebs(ik  ,ir) - ktebs(ik+1,ir)) /
     .            (kss  (ik+1,ir) - kss  (ik  ,ir))

          cdetmp(ik) = -2000.0 * rdum1**(5.0/2.0) * rdum2
        ENDDO


        IF (ikb.EQ.-1) THEN
          DO ik = iks, ike, -1
            IF     (code.EQ.1) THEN
              IF (output)
     .          WRITE(PINOUT,'(10X,3I5,2F12.5,1P,5E10.2,0P)')
     .            ik,iks,ir,ktebs(ik,ir),osm_dp6(ik,ir),
     .            -osmpei(ik,ir),-osmpmk(ik,ir),
     .            pinqe(ik,ir),osmcfe(ik,ir),cdetmp(ik)

            ELSEIF (code.EQ.2) THEN
              IF (output)
     .          WRITE(PINOUT,'(10X,3I5,2F12.5,1P,3E10.2,0P,
     .                         F10.5,1P,2E10.2,0P)')
     .            ik,iks,ir,ktebs(ik,ir),osm_dp4(ik,ir),
     .            -osmpei(ik,ir),-osmpmk(ik,ir),
     .            pinqe(ik,ir),osm_dp6(ik,ir),osmcfe(ik,ir),cdetmp(ik)

            ENDIF
          ENDDO

          IF     (code.EQ.1) THEN
            ikb = iks - 1

            DO WHILE (ikb.GT.ike.AND.cdetmp(ikb).LT.0.0)
              ikb = ikb - 1
            ENDDO

c...mark1
            IF (ikb.LT.iks-1.AND.
     .          ktebs(iks-1,ir)-ktebs(iks,ir).LT.osm_testep)
     .        ikb = iks - 1
          ELSEIF (code.EQ.2) THEN
            ikb = ike
          ENDIF

        ENDIF

        mode = 0
        DO ik = iks-1, ikb, -1
          IF (cdetmp(ik).GE.0.0) mode = 1
        ENDDO

        IF (cfs_sume(ir).LE.0.0.OR.cfs_mage(ir).LE.0.0) mode = 1
c...mark1
        IF (ktebs(iks-1,ir)-ktebs(iks,ir).LT.osm_testep) mode = 1

        IF     (mode.EQ.1.AND.model.EQ.-1) THEN
          IF     (code.EQ.1) THEN
            teadj = -1.0
          ELSEIF (code.EQ.2) THEN
            teadj =  2.0
          ENDIF
          modec = 0
          model = 1
        ELSEIF (mode.EQ.1.AND.model.EQ.2) THEN
          model = 1
          modec = modec + 1
          teadj = -0.5 * teadj
        ELSEIF (mode.EQ.0.AND.model.EQ.1.AND.modec.LT.3) THEN
          model = 2
          modec = modec + 1
          teadj = -0.5 * teadj
        ELSEIF (modec.GE.3) THEN
          model = -1
        ELSEIF (modec.EQ.0) THEN
          teadj = teadj * 2.0
          WRITE(PINOUT,*) 'INCREASING TEADJ',teadj
        ENDIF

        ik = ikb

        IF (output)
     .    WRITE(PINOUT,*) 'h -> ',iks,ike,ikm,ikr,osm_peimul(IKLO,ir),
     .      cfs_mage(ir),ikl,temod,mode,model,modec,teadj,
     .      code,ikb

        IF (code.EQ.1) THEN
          IF (output)
     .      WRITE(PINOUT,'(3I5,2F9.3,1P,7E12.4,0P,I4)')
     .        ik-1,ir,region,
     .        ktebs(ik-1,ir),osm_dp6(ik-1,ir),
     .        -osmpei(ik-1,ir) * (ksb(ik-1,ir) - ksb(ik-2,ir)),
     .        -osmpmk(ik-1,ir) * (ksb(ik-1,ir) - ksb(ik-2,ir)),
     .         osmcfe(ik-1,ir) * (ksb(ik-1,ir) - ksb(ik-2,ir)),
     .         pinqe (ik-1,ir) * (ksb(ik-1,ir) - ksb(ik-2,ir)),
     .         cdetmp(ik-1),-1.0,-1.0,0
        ELSEIF (code.EQ.2) THEN
          IF (output)
     .      WRITE(PINOUT,'(3I5,2F9.3,1P,7E12.4,0P,I4)')
     .        ik-1,ir,region,
     .        ktebs(ik-1,ir),osm_dp4(ik-1,ir),
     .        -osmpei(ik-1,ir) * (ksb(ik-1,ir) - ksb(ik-2,ir)),
     .        -osmpmk(ik-1,ir) * (ksb(ik-1,ir) - ksb(ik-2,ir)),
     .         osmcfe(ik-1,ir) * (ksb(ik-1,ir) - ksb(ik-2,ir)),
     .         pinqe (ik-1,ir) * (ksb(ik-1,ir) - ksb(ik-2,ir)),
     .         cdetmp(ik-1),-1.0,-1.0,0
        ENDIF

        DO WHILE (ik.GE.ikb.AND.status.EQ.0)

          IF     (cfs_sume(ir).LE.0.0.OR.cfs_mage(ir).LE.0.0) THEN

            WRITE(0     ,*) 'BUMMER! HIGH REGION FAILED'
            WRITE(PINOUT,*) 'BUMMER! HIGH REGION FAILED'
            status = -1

          ELSEIF (model.GT.0) THEN

            ikr = iks

            DO i2 = ikr, ik, -1

              IF     (code.EQ.1) THEN

                IF (osm_peimul(IKHI,ir).GT.0.0) THEN
                  osm_dp6(i2,ir) = osm_dp6(i2,ir) + teadj
                ELSE
                  osm_dp6(i2,ir) = osm_dp6(i2,ir) - teadj
                ENDIF

                IF (output)
     .            WRITE(PINOUT,'(A,3I4,3F10.4,1P,2E10.2,0P)')
     .              '                 ',
     .              ik,ikr,i2,ktebs(i2,ir),
     .              osm_dp6(i2,ir),teadj,-osmpei(i2,ir),-osmpmk(ik,ir)

              ELSEIF (code.EQ.2) THEN

                pinion(i2,ir) = pinion(i2,ir) * (1.0 + teadj)
                pinqe (i2,ir) = pinqe (i2,ir) + (1.0 + teadj)
c                osm_dp4(i2,ir) = osm_dp4(i2,ir) + teadj

                IF (output)
     .            WRITE(PINOUT,'(A,3I4,3F10.4,1P,3E10.2,0P)')
     .              '                 ',
     .              ik,ikr,i2,ktebs(i2,ir),
     .              osm_dp4(i2,ir),teadj,osmcfe(i2,ir),
     .              pinion(i2,ir),pinqe(i2,ir)

              ENDIF

            ENDDO

            status = 1

          ELSEIF (ik.LT.iks-1.AND.osmcde(ik,ir).EQ.osmcde(ik+1,ir)) THEN

            status = 2

          ELSEIF (ik.GT.ikl) THEN

            status = 3

          ENDIF

c...new
          IF     (code.EQ.1) THEN
            IF (output)
     .        WRITE(PINOUT,'(3I5,2F9.3,1P,7E12.4,0P,I4)')
     .          ik,ir,region,ktebs(ik,ir),osm_dp6(ik,ir),
     .         -osmpei(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .         -osmpmk(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .          osmcfe(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .          pinqe (ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .          cdetmp(ik),-1.0,-1.0,cmod2(ik)

          ELSEIF (code.EQ.2) THEN
            IF (output)
     .        WRITE(PINOUT,'(3I5,2F9.3,1P,7E12.4,0P,I4)')
     .          ik,ir,region,ktebs(ik,ir),osm_dp4(ik,ir),
     .         -osmpei(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .         -osmpmk(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .          osmcfe(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .          pinqe (ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .          cdetmp(ik),-1.0,-1.0,cmod2(ik)

          ENDIF

          ik = ik - 1
        ENDDO
c
c
c

        IF (status.EQ.0.AND.ikb.GT.ike) status = 3

        IF     (status.EQ.0.OR.status.EQ.-1) THEN
          tmpcode = tmpcode + 1

          ikl = -1
          ikr = 0
          teadj = 0.0
          temod = 0.0

          DO ik = iks, ike, -1
c...new
            IF     (code.EQ.1) THEN
              IF (output)
     .          WRITE(PINOUT,'(10X,2I5,2F12.5)')
     .            ik,ir,ktebs(ik,ir),osm_dp6(ik,ir)

            ELSEIF (code.EQ.2) THEN
              IF (output)
     .          WRITE(PINOUT,'(10X,2I5,2F12.5)')
     .            ik,ir,ktebs(ik,ir),osm_dp4(ik,ir)

            ENDIF
          ENDDO

          model = -1
        ELSEIF (status.EQ.1) THEN
          ikl = ik - 1
          cmod2(ikl) = 1
        ELSEIF (status.EQ.3) THEN
          ikl = ik + 1
          ikr = 0
          ikb = -1
          teadj = 0.0
          temod = 0.0
        ENDIF

c...new
        IF     (code.EQ.1) THEN
          IF (output)
     .      WRITE(PINOUT,'(3I5,2F9.3,1P,7E12.4,0P,2I4)')
     .        ik,ir,region,ktebs(ik,ir),osm_dp6(ik,ir),
     .       -osmpei(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .       -osmpmk(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .        osmcfe(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .        pinqe (ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .        cdetmp(ik),-1.0,-1.0,status,ikl

        ELSEIF (code.EQ.2) THEN
          IF (output)
     .      WRITE(PINOUT,'(3I5,2F9.3,1P,7E12.4,0P,2I4)')
     .        ik,ir,region,ktebs(ik,ir),osm_dp4(ik,ir),
     .       -osmpei(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .       -osmpmk(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .        osmcfe(ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .        pinqe (ik,ir) * (ksb(ik,ir) - ksb(ik-1,ir)),
     .        cdetmp(ik),-1.0,-1.0,status,ikl

        ENDIF







      ELSE
        CALL ER('FlattenTeProfile','Invalid REGION',*99)
      ENDIF



c      STOP 'xxx1'

      RETURN
99    STOP
      END



























c
c =SECTION==============================================================
c
c ======================================================================
c
c subroutine: StoreSources
c
      SUBROUTINE StoreSources(ir)
      IMPLICIT none

      INTEGER ir

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      REAL    storefrac

      CALL DB('Storing sources')

      storefrac = rel_frac
      rel_frac  = -1.0

c      WRITE(0,*) 'MARK: UPDATESOURCES B'
      CALL UpdateSources(ir)

      rel_frac = storefrac

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: UpdateSources
c
      SUBROUTINE UpdateSources(ir)
      IMPLICIT none

      INTEGER ir

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      COMMON /OPTTEMP/ osm_matcht,forcet1
      INTEGER          osm_matcht,forcet1


      INTEGER i1,i2,i3
      REAL    hold_rel_frac


      CALL DB('Updating sources')

c...  This is here because of PINASD concerns:
c      IF (.NOT.thesis) THEN
c        WRITE(0,*) 'SORRY, SOURCE RELAXATION NOT AVAILABLE IN THE '//
c     .             'STANDARD VERSION AT THE MOMENT' 
c        RETURN
c      ENDIF



      WRITE(PINOUT,*)
      WRITE(PINOUT,'(A     )') 'Relaxing sources:'
      WRITE(PINOUT,'(A,F8.3)') '  REL_FRAC =',rel_frac

      CALL RelaxSource(ir,pinion  ,pinion2  ,'PINION')
      IF (stopopt2.EQ.6.AND.rel_frac.NE.-1.0) THEN
        WRITE(0     ,*)
        WRITE(0     ,*) '*** NOT RELAXING PINREC ***'
        WRITE(0     ,*)
        WRITE(PINOUT,*)
        WRITE(PINOUT,*) '*** NOT RELAXING PINREC ***'
        WRITE(PINOUT,*)
      ELSE
        CALL RelaxSource(ir,pinrec  ,pinrec2  ,'PINREC')
      ENDIF
      CALL RelaxSource(ir,pinqe   ,pinqe2   ,'PINQE')
      CALL RelaxSource(ir,pinqi   ,pinqi2   ,'PINQI')
      CALL RelaxSource(ir,pinmol  ,pinmol2  ,'PINMOL')
      CALL RelaxSource(ir,pinena  ,pinena2  ,'PINENA')
      CALL RelaxSource(ir,pinenz  ,pinenz2  ,'PINENZ')
      CALL RelaxSource(ir,pinenm  ,pinenm2  ,'PINENM')
c      CALL RelaxSource(ir,pinvdist,pinvdist2,'PINVDIST')

      WRITE(0,*) 'NOT RELAXING STRATUM OR PLOSS DATA'
c      DO i1 = 1, MAXDATA
c        CALL RelaxSource(ir,pindata(1,1,i1),pindata2(1,1,i1),'PINDATA')
c      ENDDO

      WRITE(0,*) 'NOT RELAXING STANDARD GRID BGK DATA'
c      DO i1 = 1, MAXBGK
c        CALL RelaxSource(ir,pinbgk(1,1,i1),pinbgk2(1,1,i1),'PINBGK')
c      ENDDO

      WRITE(0,*) 'NOT RELAXING PINATOM'
c      CALL RelaxSource(ir,pinatom ,pinatom2 ,'PINATOM')

      IF (eiropacity.GT.0) CALL EstimateOpacityMultiplier      

c...  PINHALPHA is not relaxed because it is an observable, as
c     opposed to a source term used in the solution.  It doesn't matter
c     as long as the solution is converged or the PIN run time is
c     sufficient for good statistics:
c      CALL RelaxSource(ir,pinalpha,pinalpha2,'PINALPHA')
c      DO i1 = 1, 6
c        CALL RelaxSource(ir,pinline (1,1,i1,H_BALPHA),
c     .                      pinline2(1,1,i1,H_BALPHA),'PINLINE1')
c        CALL RelaxSource(ir,pinline (1,1,i1,H_BGAMMA),
c     .                      pinline2(1,1,i1,H_BGAMMA),'PINLINE2')
c      ENDDO

      CALL RelaxSource(ir,pinionz ,pinionz2 ,'PINIONZ')
      CALL RelaxSource(ir,pinz0   ,pinz02   ,'PINZ02')

      IF (.TRUE.) THEN
c        CALL RelaxSource(ir,osmpei  ,osmpei2  ,'OSMPEI')

        hold_rel_frac = rel_frac

        IF (rel_frac.NE.-1.0) THEN
          IF     (iflexopt(7).EQ. 2) THEN
            rel_frac = rel_frac * 0.5
          ELSEIF (iflexopt(7).EQ. 3) THEN
            rel_frac = rel_frac * 0.1
          ELSEIF ((iflexopt(7).EQ.10.OR.iflexopt(7).EQ.13).AND.
     .            rel_count.NE.0) THEN
c            WRITE(0,*) 'ACCELERATING MOMENTUM SOURCE RELAXATION (BUT '//
c     .                 'NOT OSMMP)'
            WRITE(0,*) 'ACCELERATING OSMMP RELAXATION'
            rel_frac = rel_frac * 3.0
c            rel_frac = rel_frac * 10.0
          ENDIF
        ENDIF

        IF (thesis) THEN
          WRITE(0,*) 'NOT RELAXING PINMP'
          WRITE(6,*) 'NOT RELAXING PINMP'
        ELSE
          CALL RelaxSource(ir,pinmp   ,pinmp2   ,'PINMP')
        ENDIF

c...    Rescaling OSMMP in CheckDensityPeak causes
c       problems when relaxing OSMMP, since OSMMP
c       can vary significantly between iterations.  I am not
c       really sure why at the moment, but likely it has
c       much to do with the volatility of PINMP in general.
c       At any rate, it would likely be better for a momentum
c       multiplier to be assigned in CheckDensityPeak, and then
c       apply it in SOL22, rather than directly rescaling OSMMP, 
c       although that leads to complications when plotting:
c        WRITE(PINOUT,*) '*** NOT RELAXING OSMMP,Qe ***'
c        WRITE(0     ,*) '*** NOT RELAXING OSMMP,Qe ***'
        CALL RelaxSource(ir,osmmp   ,osmmp2   ,'OSMMP')
        CALL RelaxSource(ir,osmqe   ,osmqe2   ,'OSMQE')

        rel_frac = hold_rel_frac

        IF (rel_frac.NE.-1.0) WRITE(0,*) 'RELAXING OSM SOURCES'
      ELSE
        WRITE(0,*) 'NOT RELAXING OSM SOURCES'
      ENDIF

c...Have to do PINASC manually since it does not have the standard
c   IK,IR format:

      WRITE(0,*) 'NOT RELAXING PINASD IN UPDATESOURCES'
c      DO i1 = 1, MAXNAS2
c        DO i2 = 1, MAXASD2
c          DO i3 = 1, MAXASS
c            pinasd(i1,i2,i3,3) = pinasd(i1,i2,i3,1)
c            pinasd(i1,i2,i3,2) =        rel_frac  * pinasd(i1,i2,i3,1) +
c     .                           (1.0 - rel_frac) * pinasd(i1,i2,i3,2)
c            pinasd(i1,i2,i3,1) = pinasd(i1,i2,i3,2)
c          ENDDO
c        ENDDO
c      ENDDO

      CALL RelaxSource(ir,pinior,pinior2,'PINIOR')
      CALL RelaxSource(ir,pinmpr,pinmpr2,'PINMPR')
      CALL RelaxSource(ir,pinqir,pinqir2,'PINQIR')
      CALL RelaxSource(ir,pinqer,pinqer2,'PINQER')

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: RelaxSource
c
      SUBROUTINE RelaxSource(mode,source,source2,tag)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER SymmetryPoint,GetModel

c     Input:
      INTEGER   mode
      REAL      source(MAXNKS,MAXNRS),source2(MAXNKS,MAXNRS),frac,
     .          integ,integ2,rdum1,rdum2,src1,src2,val_rel_frac
      CHARACTER tag*(*)

      INTEGER ik,ik1,ik2,ir,ir1,ir2,iks,ike



      INTEGER rel_method


      IF (rel_opt.NE.1.AND.rel_opt.NE.3) RETURN


      IF (iflexopt(7).EQ.1.OR.iflexopt(7).EQ.2 .OR.
     .    iflexopt(7).EQ.3.OR.iflexopt(7).EQ.10) THEN
        rel_method = 2
      ELSEIF (iflexopt(7).EQ.11) THEN
        rel_method = 3
      ELSEIF (iflexopt(7).EQ.12) THEN
        rel_method = 4
      ELSE
        rel_method = 1
      ENDIF

      IF (mode.EQ.-1) THEN
        ir1 = 1
        ir2 = nrs
      ELSE
        ir1 = mode
        ir2 = mode
      ENDIF

      CALL SetBounds

      IF     (mode.EQ.-1.AND.rel_frac.NE.-1.0.AND.stopopt.NE.78) THEN
        IF (tag.EQ.'PINION'.OR.tag.EQ.'PINREC'.OR.tag.EQ.'PINQE' .OR.
     .      tag.EQ.'PINMP' .OR.tag.EQ.'OSMMP' .OR.tag.EQ.'PINATOM') THEN
          WRITE(PINOUT,*)
          WRITE(PINOUT,'(1X,3A4,1X,A4,1X,2A12,1X,2A10,I4,1X,A)')
     .      'ik1','ik2','nks','ir','integ','integ2','frac','frac2',
     .      rel_method,tag
        ENDIF
      ELSEIF (stopopt.EQ.78.AND.tag.EQ.'PINION') THEN
        WRITE(PINOUT,*)
        WRITE(PINOUT,*) '(restoring sources)'
      ENDIF

      DO ir = ir1, ir2
        IF (idring(ir).EQ.-1) CYCLE

c...    Do not relax OSMMP for SOL28 rings:
        IF (GetModel(IKLO,ir).EQ.28.AND.tag.EQ.'OSMMP') CYCLE

        IF (rel_frac.NE.-1.0.AND.stopopt.NE.78) THEN
c          IF (tag.EQ.'PINION'.OR.tag.EQ.'PINREC') THEN
            ik1 = 1
            ik2 = nks(ir)
c          ELSE
c            ik1 = ikbound(ir,IKLO)
c            ik2 = ikbound(ir,IKHI)
c          ENDIF
          CALL CalcIntegral4(source ,ik1,ik2,ir,integ ,2)
          CALL CalcIntegral4(source2,ik1,ik2,ir,integ2,2)
          IF     (rel_method.EQ.1) THEN
            rdum1 =        rel_frac  * integ
            rdum2 = (1.0 - rel_frac) * integ2
            IF (rdum1.EQ.0.0.AND.rdum2.EQ.0.0) THEN
              frac = 0.0
            ELSE
              frac = rdum1 / (rdum1 + rdum2)
            ENDIF
          ELSEIF (rel_method.EQ.2) THEN
            IF (integ.EQ.0.0.OR.integ2.EQ.0.0) THEN
              rdum1 =        rel_frac
              rdum2 = (1.0 - rel_frac)
            ELSE
              rdum1 =        rel_frac  * (integ + integ2) / integ
              rdum2 = (1.0 - rel_frac) * (integ + integ2) / integ2
            ENDIF
            frac  = rdum1 / (rdum1 + rdum2)
          ELSEIF (rel_method.EQ.3.OR.rel_method.EQ.4) THEN
          ELSE
            CALL ER('RelaxSource','Invalid REL_METHOD option',*99)
          ENDIF

          IF (rel_method.NE.3.AND.rel_method.NE.4.AND.
     .        (tag.EQ.'PINION'.OR.tag.EQ.'PINREC'.OR.tag.EQ.'PINQE'.OR.
     .         tag.EQ.'PINMP' .OR.tag.EQ.'OSMMP' .OR.
     .         tag.EQ.'PINATOM')) THEN
            IF (mode.EQ.-1) THEN
              WRITE(PINOUT,90)
     .          ik1,ik2,nks(ir),ir,integ,integ2,frac,1.0-frac,
     .          frac*integ/(frac*integ+(1.0-frac)*integ2),' '
            ELSE
              WRITE(PINOUT,90)
     .          ik1,ik2,nks(ir),ir,integ,integ2,frac,1.0-frac,
     .          frac*integ/(frac*integ+(1.0-frac)*integ2),tag
            ENDIF
90          FORMAT(1X,3I4,1X,I4,1X,1P,2E12.4,0P,1X,3F10.5,1X,A)
          ENDIF
        ENDIF




c        IF (.FALSE..AND.stopopt.EQ.90.AND.nbr.GT.0.AND.
c     .      ir.GE.irbreak.AND.ir.LT.irwall
c     .     ) THEN
c          WRITE(0     ,*) 'NO RELAXATION ON INNER BROKEN RING ',ir
c          WRITE(PINOUT,*) 'NO RELAXATION ON INNER BROKEN RING ',ir
c
c          iks = SymmetryPoint(ir)
c        ENDIF



        DO ik = 1, nks(ir)

          IF (rel_frac.NE.-1.0) THEN

            val_rel_frac = rel_frac

            IF (s21_mode.EQ.5.AND.osm_model(IKLO,ir).EQ.24.AND.
     .          ik.LE.osm_sympt(ir).AND.citersol.GT.0.AND.
     .          .NOT.(relreset.GE.1.AND.rel_niter.EQ.1).AND.
     .          (tag.EQ.'OSMMP' .OR.tag.EQ.'PINQE'.OR.tag.EQ.'PINQI'.OR.
     .           tag.EQ.'PINION'.OR.tag.EQ.'PINMP')) THEN
c              WRITE(PINOUT,*) 'PROTECTING A:',ik,ir,tag
              rel_frac = 0.0
            ENDIF

            IF (iflexopt(6).EQ.20.AND.ikfluid(IKHI,ir).NE.nks(ir).AND.
     .          osm_model(IKHI,ir).EQ.22.AND.
c            IF (iflexopt(6).EQ.20.AND.ir.GE.14.AND.ir.LE.16.AND.
     .          ik.GT.osm_sympt(ir).AND.citersol.GT.0.AND.
     .          .NOT.(relreset.GE.1.AND.rel_niter.EQ.1).AND.
     .          (tag.EQ.'OSMMP' .OR.tag.EQ.'PINQE'.OR.tag.EQ.'PINQI'.OR.
     .           tag.EQ.'PINION'.OR.tag.EQ.'PINMP')) THEN
c              WRITE(PINOUT,*) 'PROTECTING B:',ik,ir,tag
              rel_frac = 0.0
            ENDIF

            IF (iflexopt(6).EQ.24.AND.ir.GE.14.AND.ir.LE.16.AND.
     .          ik.GT.osm_sympt(ir).AND.citersol.GT.0.AND.
     .          .NOT.(relreset.GE.1.AND.rel_niter.EQ.1).AND.
     .          (tag.EQ.'OSMMP' .OR.tag.EQ.'PINQE'.OR.tag.EQ.'PINQI'.OR.
     .           tag.EQ.'PINION')) THEN
c              WRITE(PINOUT,*) 'PROTECTING C:',ik,ir,tag
              rel_frac = 0.0
            ENDIF

          ENDIF 



          IF     (rel_frac.EQ.-1.0) THEN
            source2(ik,ir) = source(ik,ir)
          ELSEIF (stopopt.EQ.78) THEN
            source(ik,ir) =
     .        (source(ik,ir) - (1.0 - rel_frac) * source2(ik,ir)) /
     .         rel_frac
          ELSE
            IF     (rel_method.EQ.1) THEN
                source(ik,ir) =        rel_frac  * source (ik,ir) +
     .                          (1.0 - rel_frac) * source2(ik,ir)
            ELSEIF (rel_method.EQ.2) THEN
              source(ik,ir) =        frac  * source (ik,ir) +
     .                        (1.0 - frac) * source2(ik,ir)
            ELSEIF (rel_method.EQ.3) THEN
c...          Relax each cell individually:
              src1 = ABS(source (ik,ir))
              src2 = ABS(source2(ik,ir))

              IF (src1.EQ.0.0.OR.src2.EQ.0.0) THEN
                rdum1 =        rel_frac
                rdum2 = (1.0 - rel_frac)
              ELSE
                rdum1 =        rel_frac  * (src1 + src2) / src1
                rdum2 = (1.0 - rel_frac) * (src1 + src2) / src2
              ENDIF

              frac  = rdum1 / (rdum1 + rdum2)             

              rdum1 = frac

              IF (frac.GT.rel_frac) frac = rel_frac

              IF (tag.EQ.'OSMMP'.OR.tag.EQ.'PINMP') THEN
                WRITE(PINOUT,'(A,2I6,1P,2E10.2,0P,2F10.4)') 
     .            '   DETAIL: ',ik,ir,source(ik,ir),
     .                      source2(ik,ir),frac,rdum1
              ENDIF

              source(ik,ir) =        frac  * source (ik,ir) +
     .                        (1.0 - frac) * source2(ik,ir)

            ELSEIF (rel_method.EQ.4) THEN
c...          Relax each cell individually, mindless:

              frac = rel_frac

c              IF (tag.EQ.'OSMMP' .OR.tag.EQ.'PINMP'.OR.
c     .            tag.EQ.'PINION'.OR.tag.EQ.'PINREC') THEN
c                WRITE(PINOUT,'(A,2I6,1P,2E10.2,0P,2F10.4)') 
c     .            '      DETAILS: ',ik,ir,source(ik,ir),
c     .                             source2(ik,ir),frac
c              ENDIF

              source(ik,ir) =        frac  * source (ik,ir) +
     .                        (1.0 - frac) * source2(ik,ir)
            ELSE
              CALL ER('RelaxSource','Invalid REL_METHOD value',*99)
            ENDIF

          ENDIF

          IF (rel_frac.NE.-1.0) THEN

            rel_frac = val_rel_frac

          ENDIF
          
        ENDDO

      ENDDO


c...  Output:
      IF     (mode.EQ.-1.AND.rel_frac.NE.-1.0.AND.stopopt.NE.78.AND.
     .        (rel_method.EQ.3.OR.rel_method.EQ.4).AND.
     .        (tag.EQ.'PINION'.OR.tag.EQ.'PINREC'.OR.tag.EQ.'PINQE' .OR.
     .         tag.EQ.'PINMP' .OR.tag.EQ.'OSMMP' .OR.tag.EQ.'PINATOM')
     .       ) THEN

        DO ir = irsep, nrs
          IF (idring(ir).EQ.-1) CYCLE
	
          CALL CalcIntegral4(source ,ik1,ik2,ir,integ ,2)
          CALL CalcIntegral4(source2,ik1,ik2,ir,integ2,2)
	  
          IF (integ.EQ.0.0.OR.integ2.EQ.0.0) THEN
            rdum1 =        rel_frac
            rdum2 = (1.0 - rel_frac)
          ELSE
            rdum1 =        rel_frac  * (integ + integ2) / integ
            rdum2 = (1.0 - rel_frac) * (integ + integ2) / integ2
          ENDIF
          frac  = rdum1 / (rdum1 + rdum2)
	  
          IF (mode.EQ.-1) THEN
            WRITE(PINOUT,90)
     .        ik1,ik2,nks(ir),ir,integ,integ2,frac,1.0-frac,
     .        frac*integ/(frac*integ+(1.0-frac)*integ2),' '
          ELSE
            WRITE(PINOUT,90)
     .        ik1,ik2,nks(ir),ir,integ,integ2,frac,1.0-frac,
     .        frac*integ/(frac*integ+(1.0-frac)*integ2),tag
          ENDIF

91          FORMAT(1X,3I4,1X,I4,1X,1P,2E12.4,0P,1X,3F10.5,1X,A)
  	
        ENDDO

      ENDIF

      RETURN
99    STOP
      END


