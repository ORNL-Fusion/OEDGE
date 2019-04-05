c
c ======================================================================
c
c subroutine: AnalyseSolution
c
c ... assumes SetBounds has been called
c
      SUBROUTINE AnalyseSolution(fp)
      use debug_options
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_pindata
      use mod_slcom
      use mod_slout
      use mod_solparams
      use mod_solcommon
      use mod_solswitch
      IMPLICIT none

c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'comtor'
c     INCLUDE 'pindata'
c     INCLUDE 'slcom'
c     INCLUDE 'slout'

c     INCLUDE 'solparams'
c     INCLUDE 'solcommon'
c     INCLUDE 'solswitch'

      INTEGER    IN_ION    ,IN_REC    ,IN_CFP,
     .           IN_PQE    ,IN_PEI    ,IN_CFE,
     .           IN_PQI    ,IN_RAD    ,IN_CFI
      PARAMETER (IN_ION = 1,IN_REC = 2,IN_CFP = 3,
     .           IN_PQE = 4,IN_PEI = 5,IN_CFE = 6,
     .           IN_PQI = 7,IN_RAD = 8,IN_CFI = 9)

      INTEGER GetModel,SymmetryPoint
      LOGICAL ChkInt
      REAL    CalcPressure,GetJsat,GetFlux,GetCs,ATan2C,GetHeatFlux,
     .        GetGamma

      INTEGER ik,ik1,ik2,ir,i1,i2,model1(MAXNRS),model2(MAXNRS),
     .        id,maxik1,maxik2,fp,ikm,model,ring,target
      REAL    parflx,eleflx,ionflx,p,p1,p2,jsat,area,
     .        partot,eletot,iontot,maxs,rdum(20),mach,beta,alpha,
     .        radpow(MAXNKS,MAXNRS),fact,cost,flux(2),deltar,deltaz,
     .        fluxmax(2),isat,qpara,gamma,sum_target,sum_total
      REAL*8  intg(20),ddum1,ddum2,a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,
     .        f1,f2,tab,tcd
      CHARACTER*7 irtag(0:MAXNRS)
      CHARACTER*1 note1,note2

      COMMON /OPTTEMP/ osm_matcht,forcet1
      INTEGER          osm_matcht,forcet1

      DOUBLE PRECISION EPS10
      PARAMETER       (EPS10=1.0D-10)

c..temp
      REAL    delta,r1,area2

      call pr_trace('ANALYSESOLUTION','START')

c      IF (cioptg.NE.91) THEN
c        CALL WN('AnalyseSolution','Avoiding')
c        RETURN
c      ENDIF

c      WRITE(0,*) 'ANALYSING THE SOLUTION'

      CALL DB('Here in AnalyseSolution')

c
c jdemod - initialization
c
      qt = qtim 
c
c jdemod
c
c      fp = PINOUT

      DO ir = 1, MAXNRS
        irtag(ir) = '      '
      ENDDO

      irtag(irsep)   = 'IRSEP  '
      irtag(irbreak) = 'IRBREAK'
      irtag(irwall)  = 'IRWALL '
      irtag(irtrap)  = 'IRTRAP '
      irtag(nrs)     = 'NRS    '


      CALL HD(fp,'Solution analysis','SOLANAL-TITLE',5,67)
c      WRITE(fp,*)
c      WRITE(fp,*) 'Solution analysis:'

      DO i1 = 1, 20
        intg(i1) = 0.0
      ENDDO

      DO ir = irsep, nrs
        model1(ir) = osm_model(IKLO,ir)
        model2(ir) = osm_model(IKHI,ir)
      ENDDO
c
c
c     Target parameters:
c
c
      CALL HD(fp,'Target conditions','SOLANAL-TARGETS',5,67)
c      WRITE(fp,*)
c      WRITE(fp,*) 'Target conditions:'
      DO i1 = 2, 1, -1
        WRITE(fp,*)
        WRITE(fp,'(1X,A4,A5,2A9,1X,2A12,2A10,1X,A10,A6,A12)')
     .    'ir','sol','psin','rho','jsat','ne','Te','Ti','v','M','p'
c        WRITE(fp,'(1X,A4,A5,2A9,1X,2A10,2A12,1X,A10,A6,A12)')
c     .    'ir','sol','psin','rho','Te','Ti','Ne','Jsat','v','M','p'
c        DO ir = irsep, nrs
        DO ring = 1, nrs-irsep+1
          IF (ring.LT.nrs-irtrap+1) THEN
            ir = irtrap + ring
          ELSE
            ir = ring - (nrs - irtrap + 1) + irsep 
          ENDIF
          IF (idring(ir).EQ.-1) CYCLE
          IF (i1.EQ.1) THEN
            model = model2(ir)
          ELSE
            model = model1(ir)
          ENDIF
          id   = idds(ir,i1)
          p    = CalcPressure(knds(id),kteds(id),ktids(id),kvds(id))
          jsat = GetJsat     (kteds(id),ktids(id),knds(id),kvds(id))
          

          !
          ! jdemod - trying to track down a divide by zero error
          !
          if (GetCs(kteds(id),ktids(id)).ne.0.0) then
             mach = kvds(id) / GetCs(kteds(id),ktids(id))
          else
             write(0,*) 'WARNING: ANALYSESOLUTION: GetCs=0.0'
             mach = 0.0
          endif


          WRITE(fp,'(1X,I4,I5,2F9.5,1X,1P,2E12.4,0P,2F10.4,1X,1P,'//
     .             'E10.2,0P,F6.2,1P,E12.4,0P,1X,A,2F10.5,2X)')
     .      ir,model,psitarg(ir,1),rho(ir,CELL1),
     .      jsat,knds(id),kteds(id),ktids(id),kvds(id),mach,p,
     .      irtag(ir),rp(id),zp(id)
c          WRITE(fp,'(1X,I4,I5,2F9.5,1X,2F10.4,1P,2E12.4,1X,E10.2,0P'//
c     .             ' ,F6.2,1P,E12.4,0P,1X,A,2F10.5)')
c     .      ir,model,psitarg(ir,1),rho(ir,CELL1),
c     .      kteds(id),ktids(id),knds(id),jsat,kvds(id),mach,p,
c     .      irtag(ir),rp(id),zp(id)
        ENDDO
      ENDDO
c
c
c
      CALL HD(fp,'Target heat flux','SOLANAL-TARGETS-HEAT',5,67)
      sum_total = 0.0
      DO i1 = 2, 1, -1
        IF (i1.EQ.2) target = IKLO
        IF (i1.EQ.1) target = IKHI
        WRITE(fp,*)
        WRITE(fp,'(1X,A4,A5,2A9,2X,1A10,2A7,2X,A7,2A10)')
     .    'ir','sol','psin','rho','jsat','Te','Ti','gamma',
     .    'P Flux','H Flux'
        sum_target = 0.0
        DO ring = 1, nrs-irsep+1
          IF (ring.LT.nrs-irtrap+1) THEN
            ir = irtrap + ring
          ELSE
            ir = ring - (nrs - irtrap + 1) + irsep 
          ENDIF
          IF (idring(ir).EQ.-1) CYCLE
          IF (i1.EQ.1) THEN
            model = model2(ir)
          ELSE
            model = model1(ir)
          ENDIF
          id = idds(ir,i1)
          jsat  = ABS(GetJsat(kteds(id),ktids(id),knds(id),kvds(id)))
          isat  = ABS(GetFlux(target,ir))
          gamma = GetGamma   (target,ir)
          qpara = GetHeatFlux(target,ir)

          if (rp(id).ne.0.0.and.dds2(id).ne.0.0.and.
     >       (rho(ir,OUT23)-rho(ir,IN14)).ne.0.0.and.costet(id).ne.0.0)
     >       then
            WRITE(fp,'(1X,I4,I5,2F9.5,2X,1P,E10.2,0P,2F7.2,2X,'//
     .             'F7.2,1P,2E10.2,0P,F7.2,2X,F10.5,2X,3F6.2,2X,A)')
     .      ir,model,psitarg(ir,1),rho(ir,CELL1),
     .      jsat,kteds(id),ktids(id),gamma,isat,qpara,
     .      qpara/(2.0*PI*rp(id)*dds2(id)*1.0E+6),dds2(id),
     .      dds2(id)             /(rho(ir,OUT23)-rho(ir,IN14)),
     .      (dds2(id)*costet(id))/(rho(ir,OUT23)-rho(ir,IN14)),
     .      1.0/costet(id),
     .      irtag(ir)
          else 
            WRITE(fp,'(1X,I4,I5,2F9.5,2X,1P,E10.2,0P,2F7.2,2X,'//
     .             'F7.2,1P,2E10.2,0P,F7.2,2X,F10.5,2X,3F6.2,2X,A)')
     .      ir,model,psitarg(ir,1),rho(ir,CELL1),
     .      jsat,kteds(id),ktids(id),gamma,isat,qpara,
     .      rp(id),dds2(id),
     .      0.0,
     .      (rho(ir,OUT23)-rho(ir,IN14)),
     .      costet(id),
     .      irtag(ir)//' WARNING'

          endif


          sum_target = sum_target + qpara
        ENDDO
        sum_total = sum_total + sum_target
        WRITE(fp,*) 
        WRITE(fp,*) 'TARGET TOTAL = ',sum_target/1.0E+6,' MW'
      ENDDO
      WRITE(fp,*) 
      WRITE(fp,*) 'TOTAL = ',sum_total/1.0E+6,' MW'
c
c     High index "symmetry" point parameters:
c
c
      CALL HD(fp,'High index symmetry point','SOLANAL-HISYM',5,67)
      DO i1 = 2, 1, -1
        WRITE(fp,*)
        WRITE(fp,'(1X,A4,A5,1X,A10,2A10,A12,1X,2A12,1X,A10)')
     .    'ir','sol','s','Te','Ti','ne','p (with v)','p (no v)','Mach'
        DO ir = irsep, nrs
          IF (idring(ir).EQ.-1) CYCLE

          IF (i1.EQ.1) THEN
            model = model2(ir)
          ELSE
            model = model1(ir)
          ENDIF

          ik   = ikmids(ir) + 1
          p1   = CalcPressure(knbs (ik,ir),ktebs(ik,ir),
     .                        ktibs(ik,ir),kvhs (ik,ir)/qt)
          p2   = CalcPressure(knbs (ik,ir),ktebs(ik,ir),
     .                        ktibs(ik,ir),0.0         )

          if ( GetCs(ktebs(ik,ir),ktibs(ik,ir)) .ne.0.0) then 
             mach = kvhs(ik,ir) / qt / 
     .           GetCs(ktebs(ik,ir),ktibs(ik,ir))
          else
             write(0,*) 'WARNING: ANALYSESOLUTION: GetCs=0.0'
             mach = 0.0
          endif

          WRITE(fp,'(1X,I4,I5,1X,F10.4,2F10.4,1P,E12.4,1X,2E12.4,0P,
     .               1X,F10.4,1X,A)')
     .      ir,model,kss(ik,ir),ktebs(ik,ir),ktibs(ik,ir),knbs(ik,ir),
     .      p1,p2,mach,irtag(ir)
        ENDDO
      ENDDO

c
c     FSP pressure error (by excluding the plasma velocity from the
c     calculation):
c...note: Presently dip_v is set to 0.0, so this doesn't show much.
c
      IF (cmodopt.EQ.1.AND.osm_matcht.NE.0) THEN
        WRITE(fp,*)
        WRITE(fp,'(1X,A4,2A12,2X,A12)') 'ir','n2T','mnv2','mnv2 / n2T'
        DO ir = irsep, irwall-1
          IF (idring(ir).EQ.-1) CYCLE

          p1 = CalcPressure(prp_ne(ir,osm_probe),prp_te(ir,osm_probe),
     .                      prp_ti(ir,osm_probe),0.0                 )
          p2 = CalcPressure(prp_ne(ir,osm_probe),0.0                 ,
     .                      0.0                 ,dip_v (ir,osm_probe))

          ! jdemod
          if (p1.ne.0.0) then 
             WRITE(fp,'(1X,I4,1P,2E12.4,0P,1X,F12.2,A)')
     .      ir,p1,p2,p2/p1*100.0,'% '//irtag(ir)
             
          else
             WRITE(fp,'(1X,I4,1P,2E12.4,0P,1X,F12.2,A)')
     .      ir,p1,p2,0.0,'% '//irtag(ir)//' WARNING P1=0'
          endif
        ENDDO
      ENDIF
c
c
c
c
c
      CALL HD(fp,'Density peaks','SOLANAL-DENPEAK',5,67)
c      WRITE(fp,*)
c      WRITE(fp,*) 'Density peaks:'
      WRITE(fp,'(2X,A2,2(1X,2(A4,A6,A10),2X,2A6,2X))')
     .  'ir','ik1','s1','ne1','ikm','sm','nem','te1','te1+1',
     .       'ik2','s2','ne2','ikm','sm','nem','te2','te2-1'
      DO ir = irsep, irwall-1
        IF (idring(ir).EQ.-1.OR.ikbound(ir,IKLO).EQ.0.OR.
     .                          ikbound(ir,IKHI).EQ.0) CYCLE

        ik1 = 1
        ik2 = ikmids(ir)
        maxik1 = 1
        DO ik = ik1, ik2
          IF (knbs(ik,ir).GT.knbs(maxik1,ir)) maxik1 = ik
        ENDDO

        ik1 = ik2 + 1
        ik2 = nks(ir)
        maxik2 = nks(ir)
        DO ik = ik1, ik2
          IF (knbs(ik,ir).GT.knbs(maxik2,ir)) maxik2 = ik
        ENDDO

        ik1 = ikbound(ir,IKLO)
        ik2 = ikbound(ir,IKHI)

        note1 = ' '
        note2 = ' '
        IF (ik1.NE.maxik1) note1 = '*'
        IF (ik2.NE.maxik2) note2 = '*'

        maxs = ksmaxs(ir)

        if (maxs.ne.0.0) then 
          WRITE(fp,'(2X,I2,2(1X,2(I4,F6.3,1P,E10.2,0P),'//
     .           ' A,1X,2F6.2,2X))')
     .    ir,   ik1,      kss(   ik1,ir) /maxs,knbs(   ik1,ir),
     .       maxik1,      kss(maxik1,ir) /maxs,knbs(maxik1,ir),note1,
     .       ktebs(ik1,ir),ktebs(ik1+1,ir),
     .          ik2,(maxs-kss(   ik2,ir))/maxs,knbs(   ik2,ir),
     .       maxik2,(maxs-kss(maxik2,ir))/maxs,knbs(maxik2,ir),note2,
     .       ktebs(ik2,ir),ktebs(ik2-1,ir)
        else
           write(0,*) 'WARNING:ANALYSESOLUTION: KSMAXS(IR)=0',maxs
        endif

      ENDDO
c
      call pr_trace('ANALYSESOLUTION','BEFORE CONSERVATION')
c
c     Conservation:
c
c
c     Tally sources and fluxes:
c
      CALL HD(fp,'Continuity','SOLANAL-CONTINUITY',5,67)
c      WRITE(fp,*)
c      WRITE(fp,*) 'Continuity:'
      WRITE(fp,'(1X,2A4,A5,4A12,2X,A15)')
     .  'ik1','ik2','ir','targ','ion','rec','cfp','net'

      DO ir = irsep, nrs
        IF (idring(ir).EQ.-1.OR.
     .      .NOT.(model1(ir).EQ.22.OR.model1(ir).EQ.24).OR.
     .      .NOT.(model2(ir).EQ.22.OR.model2(ir).EQ.24)) CYCLE

        ik1 = ikbound(ir,IKLO)
        ik2 = ikbound(ir,IKHI)

        CALL CalcIntegral3(pinion,ik1,ik2,ir,rdum(1),2)
        CALL CalcIntegral3(pinrec,ik1,ik2,ir,rdum(2),2)
        CALL CalcIntegral3(osmcfp,ik1,ik2,ir,rdum(3),2)

        parflx = 0.0
        IF (ik1.EQ.1) THEN
          parflx = parflx + knds(idds(ir,2)) * kvds(idds(ir,2))
        ELSE
          parflx = parflx + knbs(ik1,ir)     * kvhs(ik1,ir)
        ENDIF
        IF (ik2.EQ.nks(ir)) THEN
          parflx = parflx - knds(idds(ir,1)) * kvds(idds(ir,1))
        ELSE
          parflx = parflx - knbs(ik2,ir)     * kvhs(ik2,ir)
        ENDIF

        partot = parflx + rdum(1) - rdum(2) + rdum(3)

        !
        ! jdemod
        !
        if (parflx.ne.0.0) then 
          WRITE(fp,'(1X,2I4,I5,1P,4E12.3,2X,E15.5,0P)')
     .    ik1,ik2,ir,parflx,rdum(1),-rdum(2),rdum(3),partot/parflx
        else
          WRITE(fp,'(1X,2I4,I5,1P,4E12.3,2X,E15.5,0P)')
     .    ik1,ik2,ir,parflx,rdum(1),-rdum(2),rdum(3),0.0

        endif
          
      ENDDO
c
c     Source accounting:
c
c
      call pr_trace('ANALYSESOLUTION','BEFORE ACCOUNTING')
c
c
      CALL RZero(rdum,20)

c     Target flux in the SOL:
      DO ir = irsep, irwall-1
        rdum(1) = rdum(1) - GetFlux(IKLO,ir)
        rdum(2) = rdum(2) + GetFlux(IKHI,ir)
      ENDDO

c     Target flux in the PFZ:
      DO ir = irtrap+1, nrs
        rdum(3) = rdum(3) - GetFlux(IKLO,ir) + GetFlux(IKHI,ir)
      ENDDO

c     Target flux to the inner wall:
      IF (nbr.GT.0) THEN
        DO ir = irbreak, irwall-1
          rdum(5) = rdum(5) - GetFlux(IKLO,ir)
        ENDDO
      ENDIF
c
c     Recombination and ionisation in the SOL:
      CALL VolInteg(pinion,IKLO,irsep,irwall-1,rdum(6))
      CALL VolInteg(pinrec,IKLO,irsep,irwall-1,rdum(7))
      CALL VolInteg(pinion,IKHI,irsep,irwall-1,rdum(8))
      CALL VolInteg(pinrec,IKHI,irsep,irwall-1,rdum(9))
c
c     Recombination and ionisation in the PFZ:
      CALL VolInteg(pinion,3,irtrap+1,nrs,rdum(10))
      CALL VolInteg(pinrec,3,irtrap+1,nrs,rdum(11))
c
c     Recombination and ionisation in the core:
      CALL VolInteg(pinion,3,2,irsep-1,rdum(12))
      CALL VolInteg(pinrec,3,2,irsep-1,rdum(13))

c...Total ionisation:
      rdum(17) = rdum( 6) + rdum( 8) + rdum(10) + rdum(12)
c...Total recombination:
      rdum(14) = rdum( 7) + rdum( 9) + rdum(11) + rdum(13)
      rdum(15) = rdum( 1) + rdum( 2) + rdum( 3)
      rdum(16) = rdum(15) + rdum(14)



      ir    = 2
      area  = 0.0
      area2 = 0.0      
      DO ik = 1, nks(ir)-1
        id = korpg(ik,ir)
        r1 = 0.5 * (rvertp(1,id) + rvertp(4,id))
c
C       IPP/01 - Krieger: incredibly enough, the FUJI f90 compiler
C       chokes over v**2.0 if v is negative :-(
C       fixed by v**2.0 -> v**2
c
        delta = SQRT((rvertp(1,id) - rvertp(4,id))**2 +
     .               (zvertp(1,id) - zvertp(4,id))**2)      
c
c        delta = SQRT((rvertp(1,id) - rvertp(4,id))**2.0 +
c     .               (zvertp(1,id) - zvertp(4,id))**2.0)      
c
        area  = area  + 2.0 * PI * r1  * delta
        area2 = area2 + 2.0 * PI * rxp * delta
      ENDDO 
      WRITE(fp,*)
      WRITE(fp,*) 'AREA OF CORE RING (m2) = ',area
      WRITE(fp,*) 'AREA OF CORE RING (m2) = ',area2
      
c...balance
      CALL HD(fp,'Full plasma analysis','SOLANAL-FULLPLASMA',5,67)
c      WRITE(fp,*)
c      WRITE(fp,*) 'Full plasma analysis:'
      WRITE(fp, 9) 'GLOBAL','LOW SOL','HIGH SOL','CORE','PFZ','IW'
9     FORMAT(3X,30X,6A11)
      WRITE(fp,10) 'Continuity      ion/(flx+rec)',
     .  rdum(17)/(rdum(16)+EPS10),
     .  rdum( 6)/(rdum( 1)+rdum( 7)+EPS10),
     .  rdum( 8)/(rdum( 2)+rdum( 9)+EPS10)
      WRITE(fp,14) 'Ionisation           ion/tot ',
     .  rdum( 6)/(rdum(17)+EPS10),rdum( 8)/(rdum(17)+EPS10),
     .  rdum(12)/(rdum(17)+EPS10),rdum(10)/(rdum(17)+EPS10)
      WRITE(fp,10) 'Recombination   rec/(flx+rec)',
     .  rdum(14)/(rdum(16)+EPS10),rdum( 7)/(rdum(16)+EPS10),
     .  rdum( 9)/(rdum(16)+EPS10),rdum(13)/(rdum(16)+EPS10),
     .  rdum(11)/(rdum(16)+EPS10)
      WRITE(fp,16) 'Target flux     flx/(flx+rec)',
     .  rdum(15)/(rdum(16)+EPS10),rdum( 1)/(rdum(16)+EPS10),
     .  rdum( 2)/(rdum(16)+EPS10),rdum( 3)/(rdum(16)+EPS10),
     .  rdum( 5)/(rdum(16)+EPS10)
      WRITE(fp,10) 'Source ratio         rec/flx ',
     .  rdum(14)/(rdum(15)+EPS10)
10    FORMAT(3X,A30,5(F11.3:))
14    FORMAT(3X,A30,11X,5(F11.3:))
16    FORMAT(3X,A30,3F11.3,11X,2(F11.3:))

      WRITE(fp,*)
      WRITE(fp,11) 'Ionisation                   ',
     .  rdum(17),rdum( 6),rdum( 8),rdum(12),rdum(10)
      WRITE(fp,11) 'Recombination                ',
     .  rdum(14),rdum( 7),rdum( 9),rdum(13),rdum(11)
      WRITE(fp,15) 'Target flux                  ',
     .  rdum(15),rdum( 1),rdum( 2),rdum( 3),rdum( 5)
11    FORMAT(3X,A30,1P,5(E11.3:))
15    FORMAT(3X,A30,1P,3E11.3,11X,2E11.3)
c
c
c
c
c

      CALL RZero(rdum,20)
c
c     Target flux:
      DO ir = irsep, irwall-1
        IF (model1(ir).EQ.22) rdum(1) = rdum(1) - GetFlux(IKLO,ir)
        IF (model2(ir).EQ.22) rdum(2) = rdum(2) + GetFlux(IKHI,ir)
      ENDDO
c
c     Recombination and ionisation in the SOL:
      DO ir = irsep, irwall-1
        ikm = ikmids(ir)
        IF (model1(ir).EQ.22) THEN
          DO ik = 1, ikm
            rdum(3) = rdum(3) + pinion(ik,ir) * kvols(ik,ir)
            rdum(4) = rdum(4) + pinrec(ik,ir) * kvols(ik,ir)
          ENDDO
        ENDIF
        IF (model2(ir).EQ.22) THEN
          DO ik = ikm+1, nks(ir)
            rdum(5) = rdum(5) + pinion(ik,ir) * kvols(ik,ir)
            rdum(6) = rdum(6) + pinrec(ik,ir) * kvols(ik,ir)
          ENDDO
        ENDIF
      ENDDO

      rdum(7) = rdum(1) + rdum(4)
      rdum(8) = rdum(2) + rdum(6)

      CALL HD(fp,'Attached SOL plasma analysis','SOLANAL-ATTPLASMA',
     .        5,67)
c      WRITE(fp,*)
c      WRITE(fp,*) 'Attached SOL plasma analysis:'
      WRITE(fp, 9) 'LOW SOL','HIGH SOL'
      WRITE(fp,10) 'Continuity      ion/(flx+rec)',
     .  rdum(3)/(rdum(7)+EPS10),rdum(5)/(rdum(8)+EPS10)
      WRITE(fp,10) 'Source ratio         rec/flx ',
     .  rdum(4)/(rdum(1)+EPS10),rdum(6)/(rdum(2)+EPS10)
      WRITE(fp,11) 'Ionistaion                   ',rdum(3),rdum(5)
      WRITE(fp,11) 'Recombination                ',rdum(4),rdum(6)
      WRITE(fp,11) 'Target flux                  ',rdum(1),rdum(2)

c...temp
      IF (fp.EQ.PINOUT)
     .  WRITE(79,'(A,1X,3I4,1X,I4,1P,8E12.4,0P)') '''SRCANAATT 1.00''',
     .    rel_step,rel_iter,rel_count,ir,(rdum(i1),i1=1,8)
c
c
c
c
c
      CALL RZero(rdum,20)
c
c     Target flux:
      DO ir = irsep, irwall-1
        IF (model1(ir).EQ.24) rdum(1) = rdum(1) - GetFlux(IKLO,ir)
        IF (model2(ir).EQ.24) rdum(2) = rdum(2) + GetFlux(IKHI,ir)
      ENDDO

      call pr_trace('ANALYSESOLUTION','BEFORE SOL REC/ION')

c
c     Recombination and ionisation in the SOL:
      DO ir = irsep, irwall-1
        ikm = ikmids(ir)
        IF (model1(ir).EQ.24) THEN
          DO ik = 1, ikm
            rdum(3) = rdum(3) + pinion(ik,ir) * kvols(ik,ir)
            rdum(4) = rdum(4) + pinrec(ik,ir) * kvols(ik,ir)
          ENDDO
        ENDIF
        IF (model2(ir).EQ.24) THEN
          DO ik = ikm+1, nks(ir)
            rdum(5) = rdum(5) + pinion(ik,ir) * kvols(ik,ir)
            rdum(6) = rdum(6) + pinrec(ik,ir) * kvols(ik,ir)
          ENDDO
        ENDIF
      ENDDO

      rdum(7) = rdum(1) + rdum(4)
      rdum(8) = rdum(2) + rdum(6)

      CALL HD(fp,'Detached (prescription) SOL plasma analysis',
     .           'SOLANAL-DETPLASMA',5,67)
c      WRITE(fp,*)
c      WRITE(fp,*) 'Detached (prescription) SOL plasma analysis:'
      WRITE(fp, 9) 'LOW SOL','HIGH SOL'
      WRITE(fp,10) 'Continuity      ion/(flx+rec)',
     .  rdum(3)/(rdum(7)+EPS10),rdum(5)/(rdum(8)+EPS10)
      WRITE(fp,10) 'Source ratio         rec/flx ',
     .  rdum(4)/(rdum(1)+EPS10),rdum(6)/(rdum(2)+EPS10)
      WRITE(fp,11) 'Ionistaion                   ',rdum(3),rdum(5)
      WRITE(fp,11) 'Recombination                ',rdum(4),rdum(6)
      WRITE(fp,11) 'Target flux                  ',rdum(1),rdum(2)

c...temp
      IF (fp.EQ.PINOUT)
     .  WRITE(79,'(A,1X,3I4,1X,I4,1P,8E12.4,0P)') '''SRCANADET 1.00''',
     .    rel_step,rel_iter,rel_count,ir,(rdum(i1),i1=1,8)
c
c
c
c
c
c
      IF (switch(SWION  ).EQ.1.0.OR.
     .    switch(SWION  ).EQ.2.0) THEN


        CALL HD(fp,'Low index ionisation distribution by stratum',
     .             'SOLANAL-LIION',5,67)
c        WRITE(fp,*)
c        WRITE(fp,*) 'Low index ionisation distribution by stratum:'
        WRITE(fp,'(3X,A4,4A10)') 'ir','1','2','3','4'

        DO ir = irsep, nrs
          IF (model1(ir).NE.22.AND.model1(ir).NE.24) CYCLE

          ik1 = 1
          ikm = ikmids (ir)

          CALL CalcIntegral3(pindata(1,1,H_ION1),ik1,ikm,ir,rdum(1),2)
          CALL CalcIntegral3(pindata(1,1,H_ION2),ik1,ikm,ir,rdum(2),2)
          CALL CalcIntegral3(pindata(1,1,H_ION3),ik1,ikm,ir,rdum(3),2)
          CALL CalcIntegral3(pindata(1,1,H_ION4),ik1,ikm,ir,rdum(4),2)
          CALL CalcIntegral3(pinion             ,ik1,ikm,ir,rdum(5),2)

c...temp
          CALL CalcIntegral3(pinion,ikbound(ir,IKLO),ikm,ir,rdum(6),2)
          CALL CalcIntegral3(pinrec,ikbound(ir,IKLO),ikm,ir,rdum(7),2)
          rdum(8) = knbs(ikbound(ir,IKLO),ir) *
     .              kvhs(ikbound(ir,IKLO),ir)
          rdum(9) = rdum(6) / (rdum(7) - rdum(8)+EPS10)

          WRITE(fp,'(3X,I4,4F10.4,2X,1P,5E10.2,0P,F10.2)')
     .      ir,rdum(1)/(rdum(5)+EPS10),rdum(2)/(rdum(5)+EPS10),
     .         rdum(3)/(rdum(5)+EPS10),rdum(4)/(rdum(5)+EPS10),
     .        (rdum(i1),i1=1,5),rdum(9)
        ENDDO

        CALL HD(fp,'High index ionisation distribution by stratum',
     .             'SOLANAL-HIION',5,67)
c        WRITE(fp,*)
c        WRITE(fp,*) 'High index ionisation distribution by stratum:'
        WRITE(fp,'(3X,A4,4A10)') 'ir','1','2','3','4'
        DO ir = irsep, nrs
          IF (model1(ir).NE.22.AND.model1(ir).NE.24) CYCLE

          ik1 = nks   (ir)
          ikm = ikmids(ir) + 1

          CALL CalcIntegral3(pindata(1,1,H_ION1),ikm,ik1,ir,rdum(1),2)
          CALL CalcIntegral3(pindata(1,1,H_ION2),ikm,ik1,ir,rdum(2),2)
          CALL CalcIntegral3(pindata(1,1,H_ION3),ikm,ik1,ir,rdum(3),2)
          CALL CalcIntegral3(pindata(1,1,H_ION4),ikm,ik1,ir,rdum(4),2)
          CALL CalcIntegral3(pinion             ,ikm,ik1,ir,rdum(5),2)

c...temp
          CALL CalcIntegral3(pinion,ikm,ikbound(ir,IKHI),ir,rdum(6),2)
          CALL CalcIntegral3(pinrec,ikm,ikbound(ir,IKHI),ir,rdum(7),2)
          rdum(8) = knbs(ikbound(ir,IKHI),ir) *
     .              kvhs(ikbound(ir,IKHI),ir)
          rdum(9) = rdum(6) / (rdum(7) + rdum(8)+EPS10)

          WRITE(fp,'(3X,I4,4F10.4,2X,1P,5E10.2,0P,F10.2)')
     .      ir,rdum(1)/(rdum(5)+EPS10),rdum(2)/(rdum(5)+EPS10),
     .         rdum(3)/(rdum(5)+EPS10),rdum(4)/(rdum(5)+EPS10),
     .        (rdum(i1),i1=1,5),rdum(9)
        ENDDO
      ENDIF
c
c
      call pr_trace('ANALYSESOLUTION','BEFORE ELECTRON POWER')
c
c
c
      CALL HD(fp,'Electron power','SOLANAL-ELECPOWER',5,67)
c      WRITE(fp,*)
c      WRITE(fp,*) 'Electron power:'
      WRITE(fp,'(1X,2A4,A5,4A12,2X,A15)')
     .  'ik1','ik2','ir','targ','Qe','-Pei','cfe','net'

      DO ir = irsep, nrs
        IF (idring(ir).EQ.-1) CYCLE

        IF ((model1(ir).EQ.22.OR.model1(ir).EQ.24).AND.
     .      (model2(ir).EQ.22.OR.model2(ir).EQ.24)) THEN

          ik1 = ikbound(ir,IKLO)
          ik2 = ikbound(ir,IKHI)

c...prad!
          IF (switch(SWPHELP).GE.2.0) THEN
            CALL CalcIntegral2(pinqe(1,ir),ik1,ik2,ir,ddum1,2)
            CALL CalcIntegral2(osmqe(1,ir),ik1,ik2,ir,ddum2,2)
            intg(IN_PQE) = ddum1 + ddum2
          ENDIF
          IF (switch(SWPEI  ).EQ.4.0.OR.
     .        switch(SWPEI  ).EQ.5.0)
     .      CALL CalcIntegral2(osmpei(1,ir),ik1,ik2,ir,intg(IN_PEI),2)
          IF (switch(SWPOW  ).EQ.12.0.OR.
     .        switch(SWPOW  ).EQ.13.0)
     .       CALL CalcIntegral2(osmcfe(1,ir),ik1,ik2,ir,intg(IN_CFE),2)

          eleflx = elecptarg(ir,3)
          eletot = eleflx + intg(IN_PQE) - intg(IN_PEI) + intg(IN_CFE)

          WRITE(fp,'(1X,2I4,I5,1P,4E12.3,2X,E15.5,0P)')
     .      ik1,ik2,ir,eleflx,intg(IN_PQE),
     .     -intg(IN_PEI),intg(IN_CFE),eletot/(eleflx+EPS10)
        ENDIF
      ENDDO

      CALL HD(fp,'Ion power','SOLANAL-IONPOWER',5,67)
c      WRITE(fp,*)
c      WRITE(fp,*) 'Ion power:'
      WRITE(fp,'(1X,2A4,A5,4A12,2X,A15)')
     .  'ik1','ik2','ir','targ','Qi','Pei','cfi','net'

      DO ir = irsep, nrs
        IF (idring(ir).EQ.-1) CYCLE

        IF ((model1(ir).EQ.22.OR.model1(ir).EQ.24).AND.
     .      (model2(ir).EQ.22.OR.model2(ir).EQ.24)) THEN

          ik1 = ikbound(ir,IKLO)
          ik2 = ikbound(ir,IKHI)


          IF (switch(SWPCX  ).GE.2.0)
     .      CALL CalcIntegral2(pinqi(1,ir),ik1,ik2,ir,intg(IN_PQI),2)
          IF (switch(SWPEI  ).EQ.4.0.OR.
     .        switch(SWPEI  ).EQ.5.0)
     .      CALL CalcIntegral2(osmpei(1,ir),ik1,ik2,ir,intg(IN_PEI),2)
          IF (switch(SWPOW  ).EQ.12.0.OR.
     .        switch(SWPOW  ).EQ.13.0)
     .      CALL CalcIntegral2(osmcfi(1,ir),ik1,ik2,ir,intg(IN_CFI),2)


          ionflx = ionptarg(ir,3)
          iontot = ionflx + intg(IN_PQI) + intg(IN_PEI) + intg(IN_CFI)

          WRITE(fp,'(1X,2I4,I5,1P,4E12.3,2X,E15.5,0P)')
     .      ik1,ik2,ir,ionflx,intg(IN_PQI),
     .      intg(IN_PEI),intg(IN_CFI),iontot/(ionflx+EPS10)
        ENDIF
      ENDDO


      CALL CalcRadiatedPower(radpow,1)

c REGION:
c
c 1 - IKLO
c 2 - IKHI
c 3 - ENTIRE RING
c 4 - INNER HALF-RING, BELOW X-POINT
c 5 - INNER HALF-RING, ABOVE X-POINT
c 6 - OUTER HALF-RING, BELOW X-POINT
c 7 - OUTER HALF-RING, ABOVE X-POINT
c
      call pr_trace('ANALYSESOLUTION','BEFORE VOLUME POWER')

      CALL HD(fp,'EIRENE P_rad_h VOLUME INTEGRATED POWER',
     .           'SOLANA-HRADPOW-HD',5,67)

      CALL RZero(rdum,20)

c...  SOL:
      CALL VolInteg(radpow,4,irsep,irwall-1,rdum(1))
      CALL VolInteg(radpow,5,irsep,irwall-1,rdum(2))
      CALL VolInteg(radpow,6,irsep,irwall-1,rdum(3))
      CALL VolInteg(radpow,7,irsep,irwall-1,rdum(4))

c...  PFZ:
      CALL VolInteg(radpow,IKLO,irtrap+1,nrs,rdum(5))
      CALL VolInteg(radpow,IKHI,irtrap+1,nrs,rdum(6))

c...  Core:
      CALL VolInteg(radpow,IKLO,2,irsep-1,rdum(7))
      CALL VolInteg(radpow,IKHI,2,irsep-1,rdum(8))

      CALL VolInteg(radpow,3,2,nrs,rdum(9))

c...  Convert to mega-Watts:
      DO i1 = 1, 9
        rdum(i1) = rdum(i1) * 1.0E-06
      ENDDO

      fact = 2.0 * PI * rxp

      WRITE(fp,*) 
      WRITE(fp,'(4X,A,F8.3,A)') 'TOROIDAL CIR.= ',fact,' m'

      WRITE(fp,*) 
      WRITE(fp,20) 'REGION','INTEGRAL/m','INTEGRAL'
      WRITE(fp,20) ' ','(MW/m)','(MW)'

      WRITE(fp,*) 
      WRITE(fp,21) 'Low  index SOL below X-point',rdum(1),rdum(1)*fact
      WRITE(fp,21) 'Low  index SOL above X-point',rdum(2),rdum(2)*fact
      WRITE(fp,21) 'High index SOL below X-point',rdum(3),rdum(3)*fact
      WRITE(fp,21) 'High index SOL above X-point',rdum(4),rdum(4)*fact

      WRITE(fp,*) 
      WRITE(fp,21) 'Low  index PFZ',rdum(5),rdum(5)*fact
      WRITE(fp,21) 'High index PFZ',rdum(6),rdum(6)*fact

      WRITE(fp,*) 
      WRITE(fp,21) 'Low  index core',rdum(7),rdum(7)*fact
      WRITE(fp,21) 'High index core',rdum(8),rdum(8)*fact

      WRITE(fp,*) 
      WRITE(fp,21) 'Low  index SOL', rdum(1)+rdum(2),
     .                              (rdum(1)+rdum(2))*fact
      WRITE(fp,21) 'High index SOL', rdum(3)+rdum(4),
     .                              (rdum(3)+rdum(4))*fact
      WRITE(fp,21) 'SOL', rdum(1)+rdum(2)+rdum(3)+rdum(4),
     .                   (rdum(1)+rdum(2)+rdum(3)+rdum(4))*fact
      WRITE(fp,21) 'PFZ', rdum(5)+rdum(6),
     .                   (rdum(5)+rdum(6))*fact
      WRITE(fp,21) 'Core', rdum(7)+rdum(8),
     .                    (rdum(7)+rdum(8))*fact
      WRITE(fp,21) 'Entire plasma',rdum(9),rdum(9)*fact


20    FORMAT(4X,A29,2A12  )
21    FORMAT(4X,A29,2F12.3)







      call pr_trace('ANALYSESOLUTION','BEFORE MOMENTUM')

c
c     Some momentum loss analysis:
c
      CALL HD(fp,'Momentum loss (SOL22,24 only)','SOLANAL-MOMLOSS',5,67)
c      WRITE(fp,*)
c      WRITE(fp,*) 'Momentum loss (SOL22,24 only):'
      WRITE(fp,'(1X,2A4,A5,3A12,2X,3A12)')
     .  'ik1','ik2','ir','total rec','total mloss','model mloss',
     .  'total/rec','model/rec','model/total'

      DO ir = irsep, nrs
        IF (idring(ir).EQ.-1) CYCLE

        ik1 = ikbound(ir,IKLO)
        ik2 = ikmids(ir)

        IF (model1(ir).EQ.22.OR.model1(ir).EQ.24) THEN
          CALL CalcIntegral3(pinrec,1  ,ik2,ir,rdum(1),2)
          CALL CalcIntegral3(osmmp ,1  ,ik2,ir,rdum(2),2)
          CALL CalcIntegral3(osmmp ,ik1,ik2,ir,rdum(3),2)
          
          WRITE(fp,'(1X,2I4,I5,1P,3E12.4,2X,3E12.4,0P)')
     .      ik1,ik2,ir,rdum(1),rdum(2),rdum(3),
     .      rdum(2)/(rdum(1)+EPS10),
     .      rdum(3)/(rdum(1)+EPS10),rdum(3)/(rdum(2)+EPS10)
        ENDIF
      ENDDO

      WRITE(fp,*)
      DO ir = irsep, nrs
        IF (idring(ir).EQ.-1) CYCLE

        ik1 = ikmids(ir) + 1
        ik2 = ikbound(ir,IKHI)

        IF (model2(ir).EQ.22.OR.model2(ir).EQ.24) THEN
          CALL CalcIntegral3(pinrec,1  ,ik2,ir,rdum(1),2)
          CALL CalcIntegral3(osmmp ,1  ,ik2,ir,rdum(2),2)
          CALL CalcIntegral3(osmmp ,ik1,ik2,ir,rdum(3),2)

          WRITE(fp,'(1X,2I4,I5,1P,3E12.4,2X,3E12.4,0P)')
     .      ik1,ik2,ir,rdum(1),rdum(2),rdum(3),
     .      rdum(2)/(rdum(1)+EPS10),
     .      rdum(3)/(rdum(1)+EPS10),rdum(3)/(rdum(2)+EPS10)
        ENDIF
      ENDDO


c      STOP 'Temp in AS'




c
c
c     Density peak flow and recombination:
c
c
      DO ir = 1, nrs
        IF (idring(ir).EQ.-1.OR.model1(ir).NE.24) CYCLE




      ENDDO
      DO ir = 1, nrs
        IF (idring(ir).EQ.-1.OR.model2(ir).NE.24) CYCLE



      ENDDO


c
c     Symmetry point continuity:
c

      IF (.FALSE.) THEN

c...    Figure out the maximum ion flux into the divertor.  This is 
c       specific to grid SL2 for 990429019 at the moment:

c...    Assign surface across the divertor throat:

c       Ring 19 is the outermost ring for the inner SOL that is still in
c       the divertor:
        ir = 19
        ik = 1
        a1 = DBLE(rvertp(2,korpg(ik,ir)))
        a2 = DBLE(zvertp(2,korpg(ik,ir)))

c       Ring 19 is the outermost ring for the outer SOL that is still in
c       the divertor:
        ir = 19
        ik = nks(ir)
        b1 = DBLE(rvertp(3,korpg(ik,ir)))
        b2 = DBLE(zvertp(3,korpg(ik,ir)))

        flux   (IKLO) = 0.0
        flux   (IKHI) = 0.0
        fluxmax(IKLO) = 0.0
        fluxmax(IKHI) = 0.0
        DO ir = irsep, irwall-1
          DO ik = 1, osm_sympt(ir)
            id = korpg(ik,ir)
            c1 = DBLE(rvertp(1,id))
            c2 = DBLE(zvertp(1,id))
            d1 = DBLE(rvertp(4,id))
            d2 = DBLE(zvertp(4,id))
            IF (ChkInt(a1,a2,b1,b2,c1,c2,d1,d2)) THEN
              deltar = rvertp(4,id) - rvertp(1,id)
              deltaz = zvertp(4,id) - zvertp(1,id)
	      alpha  = ATAN2C(deltar,deltaz)
              deltar = SNGL(b1 - a1)
              deltaz = SNGL(b2 - a2)
              beta   = ATAN2C(deltar,deltaz) - alpha
              cost   = COS(PI / 2.0 - beta)
c...          Find segment of divertor boundary in cell IK,IR:
              CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
              e1 = c1 + tcd * (d1 - c1)
              e2 = c2 + tcd * (d2 - c2)
              c1 = DBLE(rvertp(2,id))
              c2 = DBLE(zvertp(2,id))
              d1 = DBLE(rvertp(3,id))
              d2 = DBLE(zvertp(3,id))
              CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
              f1 = c1 + tcd * (d1 - c1)
              f2 = c2 + tcd * (d2 - c2)
              area = SQRT(SNGL((f1 - e1)**2 + (f2 - e2)**2)) * 
     .               2.0 * PI * rs(ik,ir) * eirtorfrac
              flux(IKLO) = flux(IKLO) + knbs(ik,ir) * kvhs(ik,ir) * 
     .                                  cost * bratio(ik,ir) * area
              fluxmax(IKLO) = fluxmax(IKLO) + knbs(ik,ir) * 
     .           GetCs(ktibs(ik,ir),ktebs(ik,ir)) * cost * 
     .           bratio(ik,ir) * area
c
c             jdemod - check for division by zero cases
c              
              if (GetCs(ktibs(ik,ir),ktebs(ik,ir)).ne.0.0) then 
                 WRITE(0,'(A,2I6,2F12.6,1P,2E10.2,0P,F7.3)') 
     .          '-->',ik,ir,cost,area,flux(IKLO),fluxmax(IKLO),
     .          kvhs(ik,ir)/GetCs(ktibs(ik,ir),ktebs(ik,ir))
              else
                 WRITE(0,'(A,2I6,2F12.6,1P,2E10.2,0P,F7.3)') 
     .          '-->',ik,ir,cost,area,flux(IKLO),fluxmax(IKLO),
     .           0.0
              endif


            ENDIF
          ENDDO
c...      High index contribution:
          DO ik = osm_sympt(ir)+1, nks(ir)
            id = korpg(ik,ir)
            c1 = DBLE(rvertp(1,id))
            c2 = DBLE(zvertp(1,id))
            d1 = DBLE(rvertp(4,id))
            d2 = DBLE(zvertp(4,id))
            IF (ChkInt(a1,a2,b1,b2,c1,c2,d1,d2)) THEN
              deltar = rvertp(1,id) - rvertp(4,id)
              deltaz = zvertp(1,id) - zvertp(4,id)
	      alpha  = ATAN2C(deltar,deltaz)
              deltar = SNGL(b1 - a1)
              deltaz = SNGL(b2 - a2)
              beta   = ATAN2C(deltar,deltaz) - alpha
              cost   = COS(PI / 2.0 - beta)
c...          Find segment of divertor boundary in cell IK,IR:
              CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
              e1 = c1 + tcd * (d1 - c1)
              e2 = c2 + tcd * (d2 - c2)
              c1 = DBLE(rvertp(2,id))
              c2 = DBLE(zvertp(2,id))
              d1 = DBLE(rvertp(3,id))
              d2 = DBLE(zvertp(3,id))
              CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
              f1 = c1 + tcd * (d1 - c1)
              f2 = c2 + tcd * (d2 - c2)
              area = SQRT(SNGL((f1 - e1)**2 + (f2 - e2)**2)) * 
     .               2.0 * PI * rs(ik,ir) * eirtorfrac
              flux(IKHI) = flux(IKHI) + knbs(ik,ir) * kvhs(ik,ir) * 
     .                                  cost * bratio(ik,ir) * area
              fluxmax(IKHI) = fluxmax(IKHI) + knbs(ik,ir) * 
     .           GetCs(ktibs(ik,ir),ktebs(ik,ir)) * cost * 
     .           bratio(ik,ir) * area
c
c             jdemod - check for division by zero cases
c              
              if (GetCs(ktibs(ik,ir),ktebs(ik,ir)).ne.0.0) then 
                 WRITE(0,'(A,2I6,2F12.6,1P,2E10.2,0P,F7.3)') 
     .          '-->',ik,ir,cost,area,flux(IKHI),fluxmax(IKHI),
     .          kvhs(ik,ir)/GetCs(ktibs(ik,ir),ktebs(ik,ir))
              else
                 WRITE(0,'(A,2I6,2F12.6,1P,2E10.2,0P,F7.3)') 
     .          '-->',ik,ir,cost,area,flux(IKHI),fluxmax(IKHI),
     .           0.0
              endif

            ENDIF
          ENDDO



        ENDDO

      ENDIF

      call pr_trace('ANALYSESOLUTION','BEFORE TARGET FLUX')

c     Target flux by target region:
      WRITE(fp,*)
      WRITE(fp,*) 'Target flux by region:'

      rdum = 0.0
 
      DO i1 = 1, grdntreg(IKLO)
        rdum(1) = 0.0
        DO i2 = 1, grdntseg(i1,IKLO)
          ir = grdtseg(i2,i1,IKLO)
          IF (idring(ir).EQ.BOUNDARY) CYCLE
          rdum(1) = rdum(1) - GetFlux(IKLO,ir)      
        ENDDO
        rdum(2) = rdum(2) + rdum(1)
        WRITE(fp,'(A,I4,1P,E11.3,0P,F11.1)') 
     .    '  IKLO:',i1,rdum(1),rdum(1)*ECH
      ENDDO

      DO i1 = 1, grdntreg(IKHI)
        rdum(1) = 0.0
        DO i2 = 1, grdntseg(i1,IKHI)
          ir = grdtseg(i2,i1,IKHI)
          IF (idring(ir).EQ.BOUNDARY) CYCLE
          rdum(1) = rdum(1) + GetFlux(IKHI,ir)      
        ENDDO
        rdum(2) = rdum(2) + rdum(1)
        WRITE(fp,'(A,I4,1P,E11.3,0P,F11.1)') 
     .    '  IKHI:',i1,rdum(1),rdum(1)*ECH
      ENDDO
      WRITE(fp,'(A,4X,1P,E11.3,0P,F11.1)') 
     .  ' TOTAL:',rdum(2),rdum(2)*ECH

      CALL DB('Done')
      call pr_trace('ANALYSESOLUTION','END')


      RETURN
99    STOP 'AnalyseSolution'
      END


