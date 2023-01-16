module mod_sol22_support

  !use mod_sol22_sources
  

  implicit none


contains







      SUBROUTINE SOL22Status(region,ir,deltat,serr,spts,npts,errcode)
! ======================================================================


! ======================================================================

! subroutine: SOL22Status

      use mod_params
      use mod_cgeom
      use mod_slcom
      use mod_solparams
      use mod_solcommon
      use mod_solswitch
!     INCLUDE 'params'
!     INCLUDE 'cgeom'
!     INCLUDE 'slcom'
!     INCLUDE 'solparams'
!     INCLUDE 'solcommon'
!     INCLUDE 'solswitch'
      IMPLICIT none
      REAL Clock2
      INTEGER       region,ir,npts,errcode
      REAL          deltat,per
      REAL*8  spts(mxspts)
      REAL*8        serr
      CHARACTER*128 cdum2
      CHARACTER*20  errstr(0:9)
      errstr(0) = '                 '
      errstr(1) = 'Imaginary N FOUND'
      errstr(2) = 'INVALID          '
      errstr(3) = 'Negative T  ERROR'
      errstr(4) = 'Step-size   ERROR'
      errstr(5) = 'Excessive  T Drop'
      errstr(6) = 'Negative N  ERROR'
      errstr(7) = 'NaNQ-Solver ERROR'
      errstr(8) = 'SPECIFIED ERR OPT'
            
      IF (region.EQ.IKLO) THEN
        WRITE(cdum2,'(128X)')
        deltat = Clock2() - deltat
        tim1 = INT(deltat)
        serr2(IKLO) = 0.0
        err1 = errcode
        IF (errcode.EQ.1) THEN
          serr2  (IKLO)    = simag1 /  ringlen
          ikerror(IKLO,ir) = ierror
        ELSE
          serr2  (IKLO)    = (MAX(soffset,serr) - soffset) / ringlen
          ikerror(IKLO,ir) = ierror
        ENDIF
        IF (osm_mode.GE.1) THEN

!           jdemod - changed halflen to halfringlen

          IF (ABS(errcode).EQ.1) THEN
             WRITE(0,'(A,I3,A,F5.1,A,I3,1X,A,1X,F6.2,A,F6.2,A)')'SOL22: Low  ',ir,' (',deltat,' s)  E: ',&
                  err1,errstr(INT(errcode)),simag1 / (halfringlen - soffset) * 100.0,&
                  ' - ',simag2 / (halfringlen - soffset) * 100.0,'%'
          ELSE
            WRITE(0,'(A,I3,A,F5.1,A,I3,1X,A)')'SOL22: Low  ',ir,' (',deltat,' s)  E: ',err1,errstr(INT(errcode))
          ENDIF
        ENDIF
        WRITE(cdum2(1:128),'(A,I2,1X,2F7.3,1P,E11.3,0P,A,F5.1,A,I3)')'SOL22 L: ',ir,te0,ti0,n0,' (',deltat,' s) E:',errcode
        IF (ABS(errcode).EQ.1) THEN

!           jdemod - changed halflen to halfringlen

          per = (serr - soffset) / (spts(npts) - soffset) * 100.0
          WRITE(cdum2(LEN_TRIM(cdum2):128),'(A,2G8.1,A,I3,A,I3,A,F4.1)')'  '//errstr(ABS(errcode)),simag1,simag2,&
               ' (',INT(simag1 / (halfringlen - soffset) * 100.0),'-',INT(simag2 / (halfringlen - soffset) * 100.0),'%)  O: ',&
               actswerror
        ELSEIF (ABS(errcode).GT.0) THEN
          per = (serr - soffset) / (spts(npts) - soffset) * 100.0
          WRITE(cdum2(LEN_TRIM(cdum2):128),'(A,G10.3,A,I3,A,F4.1)')'  '//errstr(ABS(errcode)),serr,' (',INT(per),'%)  O: ',&
               actswerror
        ENDIF
        IF (osm_mode.EQ.2) THEN
          WRITE(PINOUT,'(A)') cdum2(1:LEN_TRIM(cdum2))
          CALL CloseStorageFiles
          IF (errcode.LE.1.OR.switch(swerror).EQ.0.0) THEN
            CALL InsertFile('osm1tmp.dat',PINOUT)
            CALL InsertFile('osm6tmp.dat',PINOUT)
            CALL InsertFile('osm2tmp.dat',PINOUT)
            CALL InsertFile('osm3tmp.dat',PINOUT)
            CALL InsertFile('osm5tmp.dat',PINOUT)
            CALL InsertFile('osm4tmp.dat',PINOUT)
            WRITE(PINOUT,*)
          ENDIF
        ENDIF
      ELSEIF (region.EQ.IKHI) THEN
        WRITE(cdum2,'(128X)')
        deltat = Clock2() - deltat
        tim2 = INT(deltat)
        serr2(IKHI) = 0.0
        err2 = errcode
        IF (errcode.EQ.1) THEN
          serr2  (IKHI)    = simag1 / ringlen
          ikerror(IKHI,ir) = nks(ir) - ierror + 1
        ELSE
          serr2  (IKHI)    = (MAX(soffset,serr) - soffset) / ringlen
          ikerror(IKHI,ir) = nks(ir) - ierror + 1
        ENDIF
        IF (osm_mode.GE.1) THEN

!           jdemod - changed halflen to halfringlen

          IF (ABS(errcode).EQ.1) THEN
             WRITE(0,'(A,I3,A,F5.1,A,I3,1X,A,1X,F6.2,A,F6.2,A)')'SOL22: High ',ir,' (',deltat,' s)  E: ',err2,&
                  errstr(INT(errcode)),simag1 / (halfringlen - soffset) * 100.0,' - ',simag2 / (halfringlen - soffset) * 100.0,'%'
          ELSE
            WRITE(0,'(A,I3,A,F5.1,A,I3,1X,A)')'SOL22: High ',ir,' (',deltat,' s)  E: ',err2,errstr(INT(errcode))
          ENDIF
        ENDIF
        WRITE(cdum2(1:128),'(A,I2,1X,2F7.3,1P,E11.3,0P,A,F5.1,A,I3)')'SOL22 H: ',ir,te0,ti0,n0,' (',deltat,' s) E:',errcode
        IF (ABS(errcode).EQ.1) THEN

!           jdemod - changed halflen to halfringlen

          per = (serr - soffset) / (spts(npts) - soffset) * 100.0
          WRITE(cdum2(LEN_TRIM(cdum2):128),'(A,2G8.1,A,I3,A,I3,A,F4.1)')'  '//errstr(ABS(errcode)),simag1,simag2,&
               ' (',INT(simag1 / (halfringlen - soffset) * 100.0),'-',INT(simag2 / (halfringlen - soffset) * 100.0),'%)  O: ',&
               actswerror
        ELSEIF (ABS(errcode).GT.0) THEN
          per = (serr - soffset) / (spts(npts) - soffset) * 100.0
          WRITE(cdum2(LEN_TRIM(cdum2):128),'(A,G10.3,A,I3,A,F4.1)')'  '//errstr(ABS(errcode)),serr,' (',INT(per),'%)  O: ',&
               actswerror
        ENDIF
        IF (osm_mode.EQ.2) THEN
          WRITE(PINOUT,'(A)') cdum2(1:LEN_TRIM(cdum2))
          CALL CloseStorageFiles
          IF (errcode.LE.1.OR.switch(swerror).EQ.0.0) THEN
            CALL InsertFile('osm1tmp.dat',PINOUT)
            CALL InsertFile('osm6tmp.dat',PINOUT)
            CALL InsertFile('osm2tmp.dat',PINOUT)
            CALL InsertFile('osm3tmp.dat',PINOUT)
            CALL InsertFile('osm5tmp.dat',PINOUT)
            CALL InsertFile('osm4tmp.dat',PINOUT)
            WRITE(PINOUT,*)
          ENDIF
!          WRITE(74,*)
!          WRITE(74,'(A)') cdum2(1:LEN_TRIM(cdum2))

!          CALL CloseStorageFiles

!          CALL InsertFile('osm1tmp.dat',PINOUT)
!          CALL InsertFile('osm6tmp.dat',PINOUT)
!          CALL InsertFile('osm2tmp.dat',PINOUT)
!          CALL InsertFile('osm3tmp.dat',PINOUT)
!          CALL InsertFile('osm5tmp.dat',PINOUT)
!          CALL InsertFile('osm4tmp.dat',PINOUT)

!          WRITE(PINOUT,*)
!        ENDIF

! jdemod - Commented these out - the code should not exit here - error STOPs
!          are handled in the calling routine - in addition - the STOP
!          comment and the DB one are inconsistent

!        IF (err2.GE.6) STOP 'Fatal error: outer'

!        CALL DB('Done finding solution for inner half-ring')

! jdemod

        ENDIF
      ENDIF
      RETURN
99    STOP


    END SUBROUTINE SOL22Status




      SUBROUTINE CloseStorageFiles

      
! ======================================================================

! subroutine: CloseStorageFiles

      IMPLICIT   none
      CALL DB('Closing storage files')
      CLOSE(70)
      CLOSE(71)
      CLOSE(72)
      CLOSE(73)
      CLOSE(74)
      CLOSE(75)
      RETURN


    END SUBROUTINE CloseStorageFiles

      SUBROUTINE OpenStorageFiles(ir,target,ta)
! ======================================================================

! subroutine: OpenStorageFiles
      use mod_params
      use mod_slcom
!     INCLUDE 'params'
!     INCLUDE 'slcom'
      IMPLICIT   none
      INTEGER   ir,target
!      INTEGER rir,rtarg,rstep,riter
!      DATA    rir,rtarg,rstep,riter /-1, -1, -1, -1/
!      SAVE
      CHARACTER ta*(*)
!      IF (rir  .NE.ir      .OR.rtarg.NE.target  .OR.
!     .    rstep.NE.rel_step.OR.riter.NE.rel_iter) THEN
!        OPEN(UNIT=74,FILE='osm1'//ta,ACCESS='SEQUENTIAL',
!     .       STATUS='REPLACE')
!      ELSE
!        OPEN(UNIT=74,FILE='osm1'//ta,ACCESS='SEQUENTIAL',
!     .       STATUS='OLD',POSITION='APPEND')
!      ENDIF
      CALL DB('Opening storage files')
      OPEN(UNIT=74,FILE='osm1'//ta,ACCESS='SEQUENTIAL',STATUS='REPLACE')
      OPEN(UNIT=70,FILE='osm2'//ta,ACCESS='SEQUENTIAL',STATUS='REPLACE')
      OPEN(UNIT=71,FILE='osm3'//ta,ACCESS='SEQUENTIAL',STATUS='REPLACE')
      OPEN(UNIT=72,FILE='osm4'//ta,ACCESS='SEQUENTIAL',STATUS='REPLACE')
      OPEN(UNIT=73,FILE='osm5'//ta,ACCESS='SEQUENTIAL',STATUS='REPLACE')
!      rir   = ir
!      rtarg = target
!      rstep = rel_step
!      riter = rel_iter
      OPEN(UNIT=75,FILE='osm6'//ta,ACCESS='SEQUENTIAL',STATUS='REPLACE')
      RETURN

    END SUBROUTINE OpenStorageFiles



      SUBROUTINE InsertFile(fname,fp2)
        ! ======================================================================

        ! subroutine: InsertFile

        IMPLICIT   none
        INTEGER   fp2
        CHARACTER fname*(*),buffer*256
        CALL DB('Inserting file')
        OPEN(UNIT=70,FILE=fname,STATUS='OLD',ACCESS='SEQUENTIAL',ERR=15)
10      CALL ReadLine(70,buffer,1,*15,*15)
        WRITE(fp2,'(A)') buffer(1:LEN_TRIM(buffer))
        GOTO 10
15      CONTINUE
        CLOSE(70)
        RETURN


      END SUBROUTINE InsertFile


      SUBROUTINE CalcInitSrc(region,ir)

        ! ======================================================================

        ! subroutine: CalcInitSrc
        use mod_solparams
        use mod_solswitch
        use mod_solcommon
        use mod_params
        use mod_cgeom
        use mod_pindata
        use mod_slcom
        !     INCLUDE 'solparams'
        !     INCLUDE 'solswitch'
        !     INCLUDE 'solcommon'
        !     INCLUDE 'params'
        !     INCLUDE 'cgeom'
        !     INCLUDE 'pindata'
        !     INCLUDE 'slcom'
        IMPLICIT none
        INTEGER GetModel,SymmetryPoint
        INTEGER ir,region
        INTEGER ik,iks,ike
        REAL    s,helpi,alpha,parflux
        COMMON /QECOM/ qemul
        REAL*8         qemul
        COMMON /PININIT/ pininit
        LOGICAL          pininit(MAXNRS)
        COMMON /OPTTEMP/ osm_matcht,forcet1
        INTEGER          osm_matcht,forcet1
        !      WRITE(0,*) '*** NOT SETTING PININIT = .TRUE. ***'
        !      IF (region.EQ.IKHI) pininit(ir) = .TRUE.
        !      WRITE(PINOUT,*)
        !      WRITE(PINOUT,*) 'Generating initial sources:'
        CALL DB('Calculating initial PIN source estimates')
        IF (region.EQ.IKLO) THEN
           iks     =  ikbound(ir,IKLO)
           ike     =  SymmetryPoint(ir)
           parflux = -gtarg(ir,2)
        ELSE
           iks     =  SymmetryPoint(ir) + 1
           ike     =  ikbound(ir,IKHI)
           parflux = -gtarg(ir,1)
           !... REMEMBER TO UNDO ADJUSTMENT IN SOLASCV1 ALSO!
           !      IF (osm_powopt.EQ.2) WRITE(0,*) 'INITIAL PINQE AND PINQI AT 10%'
        ENDIF
        DO ik = iks, ike
           IF (region.EQ.IKLO) THEN
              s = kss2(ik,ir)
           ELSE
              s = ksmaxs2(ir) - kss2(ik,ir)
              !...    Generate momentum source:
           ENDIF
           IF (actswnmom.EQ.11.0.OR.actswnmom.EQ.12.0) THEN

              IF (s.LT.actlenmom*ringlen) THEN
                 IF (region.EQ.IKLO.AND.osm_model(IKLO,ir).EQ.24.AND.(iflexopt(6).EQ.20.OR.iflexopt(6).EQ.21.OR.&
                      iflexopt(6).EQ.22.OR.iflexopt(6).EQ.24)) THEN
                    pinmp(ik,ir) = smom0 * EXP(-(s - soffset) /(3.0*lammom * ringlen))
                 ELSE
                    pinmp(ik,ir) = smom0 * EXP(-(s - soffset) /(lammom * ringlen))
                    !            WRITE(PINOUT,* ) '--pinmp--?',ik,smom0,
                    !     .               EXP(-(s - soffset) / (lammom * ringlen)),
                    !     .               pinmp(ik,ir)
                 ENDIF
              ELSE
                 pinmp(ik,ir) = 0.0
              ENDIF
              IF (region.EQ.IKHI) pinmp(ik,ir) = -pinmp(ik,ir)
              osmmp(ik,ir) = pinmp(ik,ir)

              !       Generate ionisation source:

           ENDIF
           IF     (actswioni.EQ.0.0) THEN
              IF (osm_matcht.GT.0.AND.osm_powopt.EQ.1) THEN
                 pinion(ik,ir) = 0.0
                 IF (region.EQ.IKLO.AND.ik.EQ.iks.OR.region.EQ.IKHI.AND.ik.EQ.ike) THEN
                    pinion(ik,ir) = 0.0
                 ELSE
                    IF (region.EQ.IKLO.AND.ik.EQ.ike.OR.region.EQ.IKHI.AND.ik.EQ.iks) THEN
                       pinion(ik,ir) = 0.0
                    ELSE
                       pinion(ik,ir) = MAX(SNGL(0.15*ringlen-(s-soffset)),0.0)
                    ENDIF
                 ENDIF
              ELSE
                 pinion(ik,ir) = s0 * EXP(-(s - soffset) / ssrcdecay)
              ENDIF
           ELSEIF (actswioni.EQ.4.0) THEN
              IF (s.LT.ssrcfi) pinion(ik,ir) = s0
           ELSEIF (actswioni.EQ.6.0) THEN
              alpha         = 2.5 / (s5gausslen**2.0)
              pinion(ik,ir) = s0 * (s**5.0) * EXP(-alpha * (s**2.0))
           ELSEIF (actswioni.EQ.9.0) THEN
              pinion(ik,ir) = s0 * (s**5.0) *EXP(-s5alph * ((s + 0.5 * s5gausslen)**2.0))
           ELSE
              CALL ER('CalcInitSrc','ACTSWIONI value not supported',*99)

              !       Generate Phelpi:

           ENDIF
           IF (actswphelp.EQ.1.0.OR.actswphelp.EQ.2.0) THEN
              !...fix
              IF (knbs(ik,ir).GT.1.0E+20) THEN
                 helpi = 17.5 + (5.0 + 37.50 / ktebs(ik,ir)) *(1.0 +  0.25 / ktebs(ik,ir))
              ELSE
                 helpi = 17.5 + (5.0 + 37.50 / ktebs(ik,ir)) *(1.0 +  0.25 / ktebs(ik,ir)) *LOG10(1.0E+21 / knbs(ik,ir))
                 !          WRITE(PINOUT,*) '   ',ik,ir,ktebs(ik,ir),ktebs(ik,ir),
                 !     .                         knbs(ik,ir),LOG10(1.0E+21 / knbs(ik,ir))
              ENDIF
              IF     (osm_matcht.EQ.0) THEN
                 pinqe(ik,ir) = -helpi * pinion(ik,ir) * ECH
              ELSEIF (osm_powopt.EQ.2) THEN
                 !            pinqe(ik,ir) = 0.1 * pinqe(ik,ir)
                 pinqe(ik,ir) = -helpi * pinion(ik,ir) * ECH
                 pinqe(ik,ir) = MIN(-1.0,pinqe(ik,ir))
              ELSE
                 pinqe(ik,ir) = 0.0
                 IF (region.EQ.IKLO.AND.ik.EQ.iks.OR.region.EQ.IKHI.AND.ik.EQ.ike) THEN
                    pinqe(ik,ir) = 0.0
                 ELSEIF (region.EQ.IKLO.AND.ik.EQ.ike.OR.region.EQ.IKHI.AND.ik.EQ.iks) THEN
                    pinqe(ik,ir) = 0.0

                    !             jdemod - this assignment to pinqe would appear to make no
                    !                      sense since it is assigning a length value to pinqe
                    !                      presumably the length should be multiplied by some
                    !                      sort of source term

                 ELSE
                    !              pinqe(ik,ir) = -MAX(SNGL(0.95*(halflen-soffset)-(s-soffset))
                    !     .                            ,0.0)
                    !              pinqe(ik,ir) = -ABS(SNGL((s-soffset)-0.5*(halflen-soffset))
                    pinqe(ik,ir) = -MAX(SNGL(0.15*ringlen-(s-soffset)),0.0)
                 ENDIF
              ENDIF
           ELSEIF (actswphelp.EQ.0.0) THEN
              pinqe(ik,ir) = 0.0
           ELSE
              CALL ER('CalcInitSrc','ACTSWPHELP value not supported',*99)

              !       Generate Pcx:

           ENDIF
           IF     (actswpcx.EQ.1.0.OR.actswpcx.EQ.2.0.OR.actswpcx.EQ.5.0) THEN
              !          pinqi(ik,ir) = pinqi(ik,ir) * 0.1
              pinqi(ik,ir) = -1.5 * ktibs(ik,ir) * CEICF * pinion(ik,ir) *ECH
           ELSEIF (actswpcx.EQ.0.0) THEN
              pinqi(ik,ir) = 0.0
           ELSE
              CALL ER('CalcInitSrc','ACTSWPCX value not supported',*99)
              !...prad!
           ENDIF
           IF     (actswprad.EQ.0.0) THEN
              osmrad(ik,ir) = 0.0
           ELSEIF (actswprad.EQ.1.0) THEN
              IF (s.LT.lenr) THEN
                 osmrad(ik,ir) = lamr * prad0 * (1 - exp(-s / lamr))
              ELSE
                 osmrad(ik,ir) = 0.0
              ENDIF
              !...    NOT VALID AT THE MOMENT SINCE THIS RADIATION OPTION IS
              !       NOT APPLIED WITH THE INITIAL SOLUTION:
              !       osmrad(ik,ir) = radsrc_mult * pinqe(ik,ir)
              !      jdemod - added rectangular radiation source
           ELSEIF (actswprad.EQ.3.0) THEN
           elseif (actswprad.eq.6.0) then 
              if (s.lt.lamr) then 
                 osmrad(ik,ir) = 0.0
              elseif (s.lt.lenr) then 
                 ! jdemod
                 ! if it needs the radiation in the cell .. then I am not
                 ! sure what to put in
                 ! Radiation from the entire region is prad0 ... can't put in
                 ! the integral used elsewhere. 
                 ! Closest easy answer is total radiation divided by the distance
                 ! over which it is emitted. 
                 osmrad(ik,ir) =  prad0 /(lenr-lamr)
              else
                 osmrad(ik,ir) = 0.0
              endif
           ELSE
              IF (rel_opt.NE.0.AND.rel_frac.NE.1.0)CALL WN('CalcInitSrc','ACTSWPRAD value not supported')
           ENDIF
           !...bug! (IF statement wasn't there and SOL22 rings didn't have any
           !         ionisation in the first cell, which exaggerated imaginary
           !         solutions.)
        ENDDO
        IF (GetModel(region,ir).EQ.24) THEN
           IF (region.EQ.IKLO) THEN
              pinion(iks,ir) = 0.0
           ELSE
              pinion(ike,ir) = 0.0
           ENDIF

           !     Recombination source for SOL22p:

        ENDIF
        IF     (actswrecom.EQ.1.0) THEN
           IF (GetModel(region,ir).EQ.24) Call CalcInitRecom(region,ir)
           !...  Do nothing.
        ELSEIF (actswrecom.EQ.0.0) THEN
        ELSE
           CALL ER('CalcInitSrc','ACTSWRECOM value not supported',*99)

           !     Calculate momentum loss profile:

           !...  Assume half the pressure is lost:
           !      IF (iflexopt(4).EQ.20) THEN
           !...Need flag set so that I know data is available...
           !        pu = CalcPressure(prp_ne(ir,i1),prp_te(ir,i1),
           !     .                    prp_ti(ir,i1),dip_v (ir,i1))

           !        pl = pu * 0.5

           !        CALL CalcIntegral4(pinion,iks,ike,ir,integ,2)

           !        fact = pl / integ

           !        DO ik = iks, ike
           !          osmmp(ik,ir) = fact * ABS(pinion(ik,ir))
           !        ENDDO
           !      ENDIF
           !...  Whipe pinqi for now, since it is causing trouble on the 24-24 rings.  What happens
           !     is the global mock power multiplier ends up changing sign when pinqi terms grows,
           !     and this wreaks havok on the local mock multipliers:
        ENDIF
        IF ((stopopt3.EQ.17.OR.stopopt3.EQ.18.OR.stopopt3.EQ.19.OR.stopopt3.EQ.20).AND.osm_model(IKLO,ir).EQ.24.AND.&
             osm_model(IKHI,ir).EQ.24) THEN
           WRITE(0,*) 'BLANKING PINQI ON RING',ir
           CALL RZero(pinqi(1,ir),MAXNKS)
        ENDIF
        !....   Get rid of PINQI for now:
        IF ((iflexopt(6).EQ.20.OR.iflexopt(6).EQ.22.OR.iflexopt(6).EQ.24).AND.osm_model(IKLO,ir).EQ.24) THEN
           WRITE(PINOUT,*) 
           WRITE(PINOUT,*) '***********************'
           WRITE(PINOUT,*) '  WHIPING INNER Qi ON ',ir
           WRITE(PINOUT,*) '***********************'
           WRITE(PINOUT,*) 
           DO ik = 1, osm_sympt(ir)
              pinqi(ik,ir) = 0.0
           ENDDO
        ENDIF
        !....   Get rid of PINQI for now:
        IF ((iflexopt(6).EQ.22).AND.osm_model(IKHI,ir).EQ.22) THEN
           WRITE(PINOUT,*) 
           WRITE(PINOUT,*) '***********************'
           WRITE(PINOUT,*) '  WHIPING OUTER Qi ON ',ir
           WRITE(PINOUT,*) '***********************'
           WRITE(PINOUT,*) 
           DO ik = osm_sympt(ir)+1, nks(ir)
              pinqi(ik,ir) = 0.0
           ENDDO
        ENDIF
        RETURN
99      STOP

      END SUBROUTINE CalcInitSrc

      SUBROUTINE CalcInitRecom(region,ir)
        ! ======================================================================

        ! subroutine: CalcInitRecom

        use mod_params
        use mod_cgeom
        use mod_comtor
        use mod_pindata
        use mod_cadas
        use mod_slcom
        IMPLICIT none
        !     INCLUDE 'params'
        !     INCLUDE 'cgeom'
        !     INCLUDE 'comtor'
        !     INCLUDE 'pindata'
        !     INCLUDE 'cadas'
        !     INCLUDE 'slcom'
        INTEGER region,ir
        INTEGER GetModel
        REAL    GetEAD
        INTEGER in,i1,ik,iks,ike
        REAL    maxcoef,rdum1,temin,nemax,te,ne
        temin = 0.0
        nemax = HI
        IF (osm_preopt.GT.0.AND.(GetModel(region,ir).EQ.22.OR.GetModel(region,ir).EQ.24)) THEN
           IF     (region.EQ.IKLO) THEN
              iks = 1
              !          ike = ikbound(ir,IKLO) - 1
              ike = osm_sympt(ir)
              temin = ktebs(MAX(1,ikbound(ir,IKLO)-1),ir)
              nemax = knbs (MAX(1,ikbound(ir,IKLO)-1),ir)
           ELSEIF (region.EQ.IKHI) THEN
              !          iks = ikbound(ir,IKHI) + 1
              iks = osm_sympt(ir) + 1
              ike = nks(ir)
              temin = kteds(idds(ir,1))
              nemax = knds (idds(ir,1))
           ELSE
              CALL ER('CalcInitRecom','Invalid region',*99)
           ENDIF
        ELSEIF (GetModel(region,ir).EQ.28) THEN
           IF     (region.EQ.IKLO) THEN
              iks = 1
              ike = osm_sympt(ir)
           ELSEIF (region.EQ.IKHI) THEN
              iks = osm_sympt(ir) + 1
              ike = nks(ir)
           ELSE
              CALL ER('CalcInitRecom','Invalid region',*99)
           ENDIF
        ELSE
           CALL ER('CalcInitRecom','Error',*99)
           !      WRITE(PINOUT,*) 'RECOM IKS,IKE = ',region,iks,ike
        ENDIF
        !...    Use ADAS recombination data:
        IF (.FALSE.) THEN
           in = 0
           DO ik = iks, ike
              in = in + 1
              ptesa(in) = ktebs(ik,ir)
              pnesa(in) = knbs (ik,ir) * rizb
              pnbs (in) = knbs (ik,ir)
           ENDDO
           WRITE(year,'(i2.2)') iyearh
           CALL xxuid(useridh)
           iclass = 1
           IF (crecopt.EQ.4) THEN
              CALL ADASRD(year,1,1,iclass,in,ptesa,pnesa,pcoef)
           ELSE
              CALL OtherRec(in,ptesa,pnesa,pcoef,crecopt)
           ENDIF
           maxcoef = 0.0
           DO i1 = 1, in
              maxcoef = MAX(maxcoef,pcoef(i1,1))
           ENDDO
           DO i1 = 1, in
              IF (pcoef(i1,1).EQ.0.0) pcoef(i1,1) = maxcoef
           ENDDO
           in = 0
           DO ik = iks, ike
              in = in + 1
              pinrec(ik,ir) = 1.0 * (pnesa(in) * pcoef(in,1)) * pnbs(in)
           ENDDO
           !...    Use Lyman transparent AMJUEL (EIRENE database) recombination data:
        ELSEIF (eiropacity.EQ.0.OR.eiropacity.EQ.1.OR.eiropacity.EQ.2.OR.eiropacity.EQ.6.OR.eiropacity.EQ.-3) THEN
           DO ik = iks, ike
              te = MAX(temin,ktebs(ik,ir))
              ne = MIN(nemax,knbs (ik,ir))
              !          pinrec(ik,ir) = GetEAD(te,ne,3,'H.4 ') *
              !     .                    1.0E-06 * knbs(ik,ir) * knbs(ik,ir) *
              !     .                    eirscale(11)
              pinrec(ik,ir) = GetEAD(te,ne,3,'H.4 ') *1.0E-06 * knbs(ik,ir) * knbs(ik,ir)
           ENDDO
           !...    Use non-Lyman opaque data with local multiplier:
           !...    Not sure this belongs here.  We are really looking for an initial guess here, and
           !       it may not be appropriate to be using PINATOM data, even if it is available...
        ELSEIF ((eiropacity.EQ.3.OR.eiropacity.EQ.4).AND.tagpinatom.AND.tagmulrec) THEN
           STOP 'CALCINITRECOM THIS SHOULD NOT BE CALLED'
           DO ik = iks, ike
              te = MAX(temin,ktebs(ik,ir))
              ne = MIN(nemax,knbs (ik,ir))
              !          pinrec(ik,ir) = GetEAD(te,ne,3,'H.4 ') *
              !     .                    1.0E-06 * knbs(ik,ir) * knbs(ik,ir) *
              !     .                    eirscale(11) * mulrec(ik,ir)
              pinrec(ik,ir) = GetEAD(te,ne,3,'H.4 ') *1.0E-06 * knbs(ik,ir) * knbs(ik,ir) *mulrec(ik,ir)
           ENDDO
           ! *NOT QUITE RIGHT*
           !...    Use Lyman opaque AMJUEL (EIRENE database) recombination data:
        ELSEIF (eiropacity.EQ.3.OR.eiropacity.EQ.4.OR.(eiropacity.EQ.5.AND.s28mode.LT.2.0).OR.eiropacity.EQ.-4.OR.&
             eiropacity.EQ.-5) THEN
           DO ik = iks, ike
              te = MAX(temin,ktebs(ik,ir))
              ne = MIN(nemax,knbs (ik,ir))
              pinrec(ik,ir) = GetEAD(te,ne,4,'H.4 ') *1.0E-06 * knbs(ik,ir) * knbs(ik,ir)
           ENDDO

           ! ***NOT IN USE BECAUSE THE OLD SOL24/28 CRASHES***

           !...    Use Lyman alpha opaque AMJUEL (EIRENE database) recombination data:
        ELSEIF (eiropacity.EQ.5) THEN
           DO ik = iks, ike
              te = MAX(temin,ktebs(ik,ir))
              ne = MIN(nemax,knbs (ik,ir))
              pinrec(ik,ir) = GetEAD(te,ne,27,'H.4 ') *1.0E-06 * knbs(ik,ir) * knbs(ik,ir)
           ENDDO
        ELSE
           CALL ER('CalcInitRecom','Unrecognized option',*99)
        ENDIF
        RETURN
99      WRITE(0,*) 'IR,GETMODEL=',ir,region,GetModel(region,ir)
        STOP

        ! ======================================================================

      END SUBROUTINE CalcInitRecom





      subroutine echosol(s1,s2,sp,coment,outer,inner)
      use mod_params
      use mod_solparams
      use mod_solswitch
      use mod_solcommon
      !use mod_printopt
      use mod_cgeom
      use mod_comtor
      use mod_slcom
      use mod_sol22_output
      
!     This subroutine echoes the input values to standard out - it also
!     prints the additional calculated values - and is followed by the
!     tabular output of s,te,ti,n and v at each point S ... which would
!     be suitable for plotting on a spreadsheet or may be plotted by
!     calling GHOST routines if the graph option is set.

!     include 'params'
!     include 'solparams'
!     include 'solswitch'
!     include 'solcommon'

!     include 'printopt'

!     include 'cgeom'
!     include 'comtor'
! slmod begin - new
!     INCLUDE 'slcom'
      implicit none
! slmod end

      character*(*) s1,s2,sp,coment,outer,inner
      
      INTEGER fp,i1,i2
      integer i,irlim,ir,cnt,in
      real*8 spts(mxspts)

      character*20 errstr(9)

!      character*24 sp
!      character  s1*6,s2*12


!      CHARACTER  COMENT*80

!     Initialization

!      sp = '                        '
!      s1 = '      '
!      s2 = '            '

      real swtmp
      errstr(1) = 'Imaginary N FOUND'
      errstr(2) = 'INVALID          '
      errstr(3) = 'Negative T  ERROR'
      errstr(4) = 'Step-size   ERROR'
      errstr(5) = 'Excessive  T Drop'
      errstr(6) = 'Negative N  ERROR'
      errstr(7) = 'NaNQ-Solver ERROR'

      errstr(8) = 'SPECIFIED ERR OPT'
      IF (CIOPTO.EQ.0.or.ciopto.eq.2.or.ciopto.eq.3.or.ciopto.eq.4) THEN
         IRLIM = IRWALL
      ELSEIF (CIOPTO.EQ.1) THEN
         IRLIM = NRS

      ENDIF
      call prb
      CALL PRC ('  SOL OPTION 22:  SUB-OPTIONS AND RESULTS SUMMARY')
      call prb
      call prc (s1//'SOL22: SUMMARY OF OPTIONS AND INPUT VALUES')

!     Options

      call prb
      CALL PRQ (S1//'GAMMA CORRECTION FACTOR IN GAMMAI  ', GAMCOR)

      CALL PRQ (S1//'GAMMA CORRECTION FACTOR IN GAMMAE  ', GAMECOR)

      call prb
      if (recfrac.eq.1.0) then
         CALL PRC(S1//'RECYCLING SOURCE FRACTION IS OFF')
         CALL PRC(S1//'- RECYCLING SOURCE FRACTION = 1.0')
      elseif (recfrac.ne.1.0) then
         CALL PRC(S1//'RECYCLING SOURCE FRACTION IS ON')
         CALL PRQ(S1//'- RECYCLING SOURCE FRACTION = ',RECFRAC)

      endif
      CALL PRB
      CALL PRI (S1//'INITIAL NUMBER OF RUNGE-KUTTA STEPS BETWEEN'//' EACH GRID KNOT IN SOLVER:', NDIV)


!     Indicate Forced Te=Ti or NOT

      CALL PRB
      if (forcet.eq.0) then
         CALL PRC (S1//'T FORCE OPTION 0: TE AND TI ARE FOLLOWED'//' INDEPENDENTLY.')
      ELSEIF (FORCET.EQ.1) THEN
         CALL PRC (S1//'T FORCE OPTION 1: TE AND TI ARE LOCKED'//' TOGETHER.')
         CALL PRC (SP//'BASED ON A COMBINED'//' ENERGY EQUATION.')

      endif

!     Velocity Error Correction Option:

      call prb
      if (velsw.eq.0) then
         CALL PRC (S1//'VEL/COR OPT   0 : VELOCITY SET TO CS WHEN'//' IMAGINARY RESULT FOUND.')
      elseif (velsw.eq.1) then
         CALL PRC (S1//'VEL/COR OPT   1 : VELOCITY HELD CONSTANT WHEN'//' IMAGINARY RESULT FOUND.')
      elseif (velsw.eq.2) then
         CALL PRC (S1//'VEL/COR OPT   2 : PRESSURE SET TO MINIMUM'//' NECESSARY TO AVOID IMAGINARY.')
         CALL PRC (SP//'DENSITY IS SET FOR THIS'//' PRESSURE VALUE.')
         CALL PRC (SP//'VELOCITY IS SET USING FLUX'//' CONSERVATION.')
         CALL PRC (SP//'ADDITIONAL PRESSURE REQUI'//'RED IS CARRIED FORWARD.')
      elseif (velsw.eq.3) then
         CALL PRC (S1//'VEL/COR OPT   3 : PRESSURE SET TO MINIMUM'//' NECESSARY TO AVOID IMAGINARY.')
         CALL PRC (SP//'DENSITY IS SET FOR THIS'//' PRESSURE VALUE.')
         CALL PRC (SP//'VELOCITY IS SET USING FLUX'//' CONSERVATION.')
         CALL PRC (SP//'ADDITIONAL PRESSURE REQUI'//'RED IS NOT CARRIED FORWARD.')

      endif
      if (lensind.eq.0) then
         CALL PRC (S1//'LENGTH OPTION 0 : IONIZATION LENGTHS ARE IN'//' ABSOLUTE UNITS (M)'//' UNLESS')
         call prc (sp//'INDICATED OTHERWISE.')
      elseif (lensind.eq.1) then
         CALL PRC (S1//'LENGTH OPTION 1 : IONIZATION LENGTHS ARE IN'//' RELATIVE UNITS * SMAX'//' UNLESS')
         call prc (sp//'INDICATED OTHERWISE.')


!     Main Ionization options

      endif
      call prb

      CALL PRC(S1//'MAIN IONIZATION OPTION: ')

      call prb

      if (switch(swion).eq.0.0) then
        CALL PRB
        CALL PRC (S1//'IONIZATION OPT 0: EXPONENTIAL'//' DECAY IONIZATION SOURCE')

        CALL PRB
        CALL PRQ (SP//'LENGTH OF IONIZATION SOURCE          ',LENSFI)

        CALL PRQ (SP//'DECAY LENGTH OF IONIZATION SOURCE    ',LAMS)

      elseif (switch(swion).eq.1.0) then
        CALL PRC (S1//'IONIZATION OPT 1: PIN DATA READ FOR'//' IONIZATION SOURCE')
        CALL PRC (SP//'DATA IS NORMALIZED TO NOVO')
        CALL PRB
        CALL PRC (SP//'DATA FOR IONIZATION SOURCE IS RETURNED')
        CALL PRC (SP//'FROM A PIN/NIMBUS RUN AND IS LINEARLY')
        CALL PRC (SP//'INTERPOLATED FOR VALUES BETWEEN GRID')
        CALL PRC (SP//'POINTS. THE IONIZATION SOURCE IS')
        CALL PRC (SP//'INTEGRATED OVER THE HALF-RING AND SET')
        CALL PRC (SP//'EQUAL TO THE FLOW TO THE TARGET NOVO')

        CALL PRC (SP//'AS A NORMALIZATION FACTOR.')

      elseif (switch(swion).eq.2.0) then
        CALL PRC (S1//'IONIZATION OPT 2: PIN DATA READ FOR'//' IONIZATION SOURCE')
        CALL PRC (SP//'DATA IS UNNORMALIZED.')
        CALL PRB
        CALL PRC (SP//'DATA FOR IONIZATION SOURCE IS RETURNED')
        CALL PRC (SP//'FROM A PIN/NIMBUS RUN AND IS LINEARLY')
        CALL PRC (SP//'INTERPOLATED FOR VALUES BETWEEN GRID')
        CALL PRC (SP//'POINTS. THE DATA IS USED AS-IS AND IS')
        CALL PRC (SP//'NOT NORMALIZED TO THE TARGET FLUX FOR')

        CALL PRC (SP//'THE RING.')

      elseif (switch(swion).eq.8.0) then
        CALL PRC (S1//'IONIZATION OPT 8: PIN DATA READ FOR'//' IONIZATION SOURCE')
        CALL PRC (SP//'DATA IS UNNORMALIZED')
        CALL PRB
        CALL PRC (SP//'THE INTEGRATED STRENGTH OF THE PIN'//' IONIZATION')
        CALL PRC (SP//'IS USED TO NORMALIZE THE ANALYTIC OPTION')
        CALL PRC (SP//'SPECIFIED FOR THE INITIAL PLASMA ON'//' SUBSEQUENT')

        CALL PRC (SP//'ITERATIONS.')

      elseif (switch(swion).eq.3.0) then
        CALL PRC(S1//'IONIZATION OPT 3: IMPOSED TRIANGULAR'//' IONIZATION SOURCE')
        CALL PRQ(SP//'EXTENDING FROM   :',LENSST)
        CALL PRQ(SP//'          TO     :',LENSFI)
        CALL PRC(SP//'INTEGRAL OF SOURCE NORMALIZED'//' TO RING TARGET FLUX.')

        call prb

      elseif (switch(swion).eq.4.0) then
        CALL PRC(S1//'IONIZATION OPT 4: IMPOSED RECTANGULAR'//' IONIZATION SOURCE')
        CALL PRQ(SP//'EXTENDING FROM   :',LENSST)
        CALL PRQ(SP//'          TO     :',LENSFI)
        CALL PRC(SP//'INTEGRAL OF SOURCE NORMALIZED'//' TO RING TARGET FLUX.')

        call prb

      elseif (switch(swion).eq.5.0) then
        CALL PRC(S1//'IONIZATION OPT 5: ALGORITHMIC RECT/TRI'//' IONIZATION SOURCE')
        CALL PRC(SP//'IF NT > 1.0E19  - TRIANGULAR SOURCE')
        CALL PRC(SP//'   FROM L1= (13 - 10*TET) M (TET <  1.3 EV)')
        CALL PRC(SP//'     OR L1= 0.0           M (TET >= 1.3 EV)')
        CALL PRC(SP//'     TO L2=L1+2 (M)')
        CALL PRC(SP//'IF NT <= 1.0E19 - RECTANGULAR SOURCE')
        CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 13-TET (M)'//' (TET<10EV)')
        CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 2 (M)'//' (TET>10EV)')

        call prb

      elseif (switch(swion).eq.6.0) then
        CALL PRC(S1//'IONIZATION OPT 6: IMPOSED S**5 GAUSSIAN'//' IONIZATION SOURCE')
        CALL PRC(SP//'OF FORM:  A * S**5 * EXP(-ALPHA * S**2)')
        CALL PRC(SP//'EXTENDING FROM :  0.0 (M)')
        CALL PRQ(SP//'CUT OFF AT     :',LENSFI)
        CALL PRC(SP//'WHERE          : ALPHA = 2.5 / WF**2')
        CALL PRQ(SP//'WITH WIDTH FACTOR (WF):',LAMS)
        CALL PRC(SP//'INTEGRAL OF SOURCE NORMALIZED'//' to ring target flux.')

        call prb

      elseif (switch(swion).eq.7.0) then
        CALL PRC(S1//'IONIZATION OPT 7: ALGORITHMIC'//' RECT/S**5GAUSS IONIZATION'//' SOURCE')
        CALL PRC(SP//'IF NT > 1.0E19  - S5 GAUSSIAN SOURCE')
        CALL PRC(SP//'   WITH WIDTH FACTOR WF = (14-10TET)'//' (M) (TET <  1.3EV)')
        CALL PRC(SP//'   WITH WIDTH FACTOR WF = 1.0       '//' (M) (TET >= 1.3EV)')
        CALL PRC(SP//'   FROM L1= 0.0 TO L2 = 1/2 RING LENGTH')
        CALL PRC(SP//'IF NT <= 1.0E19 - RECTANGULAR SOURCE')
        CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 13-TET (M)'//' (TET<10EV)')
        CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 2 (M)'//' (TET>10EV)')

        call prb

      elseif (switch(swion).eq.9.0) then
        CALL PRC(S1//'IONIZATION OPT 9: IMPOSED OFFSET S**5'//' GAUSSIAN IONIZATION SOURCE')
        CALL PRC(SP//'OF FORM:  A * (S+L)**5 *'//' EXP(-ALPHA * (S+L)**2)')
        CALL PRC(SP//'EXTENDING FROM :  0.0 (M)')
        CALL PRQ(SP//'   CUT OFF AT  :',LENSFI)
        CALL PRC(SP//'WHERE          : ALPHA = 2.5 / WF**2')
        CALL PRQ(SP//'WITH WIDTH FACTOR (WF):',LAMS)
        CALL PRC(SP//'AND L = WF/2.0')
        CALL PRC(SP//'INTEGRAL OF SOURCE NORMALIZED'//' TO RING TARGET FLUX.')

        CALL PRB

      elseif (switch(swion).eq.10.0) then
        CALL PRC(S1//'IONIZ. OPTION 10: ALGORITHMIC RECT/'//'OFFSET S**5GAUSS IONIZATION SOURCE')
        CALL PRC(SP//'IF NT > 1.0E19  - OFFSET S5'//' GAUSSIAN SOURCE')
        CALL PRC(SP//'   WITH WIDTH FACTOR WF = (28-20TET)'//' (M) (TET <  1.3EV)')
        CALL PRC(SP//'   WITH WIDTH FACTOR WF = 2.0       '//' (M) (TET >= 1.3EV)')
        CALL PRC(SP//'   FROM L1= 0.0 TO L2 = 1/2 RING LENGTH')
        CALL PRC(SP//'IF NT <= 1.0E19 - RECTANGULAR SOURCE')
        CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 13-TET (M)'//' (TET<10EV)')
        CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 2 (M)'//' (TET>10EV)')

        call prb

!     Initial Ionization Option


      endif

!        call prb
!        call prr (s1//'PIN ITERATED IONIZATION SOURCE OPTION'//
!     >                ' USED:',
!     >           switch(swion))

      if (switch(swion).eq.1.or.switch(swion).eq.2.or.switch(swion).eq.8) then
        call prb
        call prc (s1//'INITIAL SEED PLASMA OPTION:')

        call prb

        if (e2dstart.eq.1) then
          call prc(s1//'Seed plasma is read in from corresponding'//' Edge2D case')

          call prb

!         regular seed plasma options

        elseif (e2dstart.eq.0) then

          if (switch(swioni).eq.0.0) then
            CALL PRC(S1//'INIT IONIZ OPT 0: EXPONENTIAL DECAY'//' IONIZATION SOURCE')
            CALL PRQ (SP//'LENGTH OF IONIZATION SOURCE          ',LENSFI)

            CALL PRQ (SP//'DECAY LENGTH OF IONIZATION SOURCE    ',LAMS)

          elseif (switch(swioni).eq.3.0) then
            CALL PRC(S1//'INIT IONIZ OPT 3: IMPOSED TRIANGULAR'//' IONIZATION SOURCE')
            CALL PRQ(SP//'EXTENDING FROM   :',LENSST)
            CALL PRQ(SP//'          TO     :',LENSFI)

            CALL PRC(SP//'INTEGRAL OF SOURCE NORMALIZED'//' TO RING TARGET FLUX.')

          elseif (switch(swioni).eq.4.0) then
            CALL PRC(S1//'INIT IONIZ OPT 4: IMPOSED'//' RECTANGULAR IONIZATION SOURCE')
            CALL PRQ(SP//'EXTENDING FROM   :',LENSST)
            CALL PRQ(SP//'          TO     :',LENSFI)

            CALL PRC(SP//'INTEGRAL OF SOURCE NORMALIZED'//' TO RING TARGET FLUX.')

          elseif (switch(swioni).eq.5.0) then
            CALL PRC(S1//'INIT IONIZ OPT 5: ALGORITHMIC'//' RECT/TRI IONIZATION SOURCE')
            CALL PRC(SP//'IF NT > 1.0E19  - TRIANGULAR SOURCE')
            CALL PRC(SP//'   FROM L1= (13 - 10*TET) M (TET <'//'  1.3 EV)')
            CALL PRC(SP//'     OR L1= 0.0           M (TET >'//'= 1.3 EV)')
            CALL PRC(SP//'   TO L2=L1+2 (M)')
            CALL PRC(SP//'IF NT <= 1.0E19 - RECTANGULAR SOURCE')
            CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 13-TET (M)'//' (TET<10EV)')

            CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 2 (M)'//' (TET>10EV)')

          elseif (switch(swioni).eq.6.0) then
            CALL PRC(S1//'INIT IONIZ OPT 6: IMPOSED S**5'//' GAUSSIAN IONIZATION SOURCE')
            CALL PRC(SP//'OF FORM:  A * S**5 * EXP(-ALPHA'//' * S**2)')
            CALL PRC(SP//'EXTENDING FROM :  0.0 ')
            CALL PRQ(SP//'    CUT OFF AT :',LENSFI)
            CALL PRC(SP//'WHERE          : ALPHA = 2.5 / WF**2')
            CALL PRQ(SP//'WITH WIDTH FACTOR (WF):',LAMS)

            CALL PRC(SP//'INTEGRAL OF SOURCE NORMALIZED'//' TO RING TARGET FLUX.')

          elseif (switch(swioni).eq.7.0) then
            CALL PRC(S1//'INIT IONIZ OPT 7: ALGORITHMIC'//' RECT/S**5GAUSS IONIZATION')
            CALL PRC(SP//'IF NT > 1.0E19  - S5 GAUSSIAN SOURCE')
            CALL PRC(SP//'   WITH WIDTH FACTOR WF = (14-10TET)'//' (M) (TET <  1.3EV)')
            CALL PRC(SP//'   WITH WIDTH FACTOR WF = 1.0       '//' (M) (TET >= 1.3EV)')
            CALL PRC(SP//'   FROM L1= 0.0 TO L2 = 1/2 RING LENGTH')
            CALL PRC(SP//'IF NT <= 1.0E19 - RECTANGULAR SOURCE')
            CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 13-TET (M)'//' (TET<10EV)')

            CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 2 (M)'//' (TET>10EV)')

          elseif (switch(swioni).eq.9.0) then
            CALL PRC(S1//'INIT IONIZ OPT 9: IMPOSED'//' OFFSET S**5 GAUSSIAN')
            CALL PRC(SP//'OF FORM:  A * (S+L)**5 *'//' EXP(-ALPHA * (S+L)**2)')
            CALL PRC(SP//'EXTENDING FROM : 0.0 (M)')
            CALL PRQ(SP//'   CUT OFF AT  :',LENSFI)
            CALL PRC(SP//'WHERE          : ALPHA = 2.5 / WF**2')
            CALL PRQ(SP//'WITH WIDTH FACTOR (WF):',LAMS)
            CALL PRC(SP//'AND L = WF/2.0')

            CALL PRC(SP//'INTEGRAL OF SOURCE NORMALIZED'//' TO RING TARGET FLUX.')

          elseif (switch(swioni).eq.10.0) then
            CALL PRC(S1//'INIT IONIZ OP 10: ALGORITHMIC'//' RECT/OFFSET S**5GAUSS'//' IONIZATION SOURCE')
            CALL PRC(SP//'IF NT > 1.0E19  - OFFSET S5'//' GAUSSIAN SOURCE')
            CALL PRC(SP//'   WITH WIDTH FACTOR WF = (28-20TET)'//' (M) (TET <  1.3EV)')
            CALL PRC(SP//'   WITH WIDTH FACTOR WF = 2.0       '//' (M) (TET >= 1.3EV)')
            CALL PRC(SP//'   FROM L1= 0.0 TO L2 = 1/2 RING LENGTH')
            CALL PRC(SP//'IF NT <= 1.0E19 - RECTANGULAR SOURCE')
            CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 13-TET (M)'//' (TET<10EV)')

            CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 2 (M)'//' (TET>10EV)')

          elseif (switch(swioni).eq.11.0) then

            CALL PRC(s1//'INIT IONIZ OP 11: IONIZATION'//' SOURCE DATA READ FROM'//' EDGE2D INPUT FOR CASE.')

          elseif (switch(swioni).eq.12.0) then
            CALL PRC(S1//'INIT IONIZ OP 12: PIN IS RUN'//' WITH AN EDGE2D BACKGROUND')
            call prc(sp//'IN THE SOL ONCE BEFORE SOL 22 IS INVOKED.')

            CALL PRC(SP//'THIS GENERATES POWER TERMS BUT WILL'//' NOT INCLUDE PUFFING.')

          elseif (switch(swioni).eq.13.0) then
            CALL PRC(S1//'INIT IONIZ OP 13: PIN IS RUN'//' WITH AN EDGE2D BACKGROUND')
            call prc(sp//'IN THE SOL TWICE BEFORE SOL 22 IS INVOKED.')

            CALL PRC(SP//'THIS GENERATES POWER TERMS AND'//' PROPER PUFFING.')

          elseif (switch(swioni).eq.14.0) then
            CALL PRC(S1//'INIT IONIZ OPT 14: PIN IS RUN'//' WITH AN EDGE2D BACKGROUND')
            call prc(sp//'EVERYWHERE TWICE BEFORE SOL 22 IS INVOKED.')
            CALL PRC(SP//'THIS GENERATES POWER TERMS AND'//' PROPER PUFFING.')

            CALL PRC(SP//'THIS SHOULD BE USED'//' IN CONJUNCTION'//' WITH CORE OPTION -1.')

          elseif (switch(swioni).eq.15.0) then
            CALL PRC(s1//'INIT IONIZ OP 15: IONIZATION'//' SOURCE DATA READ FROM'//' EDGE2D INPUT FOR CASE.')

            call prc(sp//'EDGE2D PLASMA ASSIGNED AS OLD'//' FOR ENTIRE BACKGROUND')

          elseif (switch(swioni).eq.16.0) then
            CALL PRC(S1//'INIT IONIZ OP 16: PIN IS RUN'//' WITH A PREVIOUSLY CALCULATED DIVIMP BACKGROUND')
            call prc(sp//'EVERYWHERE - BEFORE SOL 22 IS INVOKED.')

            CALL PRC(SP//'THIS GENERATES POWER TERMS BUT WILL'//' NOT INCLUDE PUFFING.')

          elseif (switch(swioni).eq.17.0) then
            CALL PRC(S1//'INIT IONIZ OP 17: PIN IS RUN'//' WITH A PREVIOUSLY CALCULATED DIVIMP BACKGROUND')
            call prc(sp//'EVERYWHERE - BEFORE SOL 22 IS INVOKED.')
            CALL PRC(SP//'THIS GENERATES POWER TERMS BUT WILL'//' NOT INCLUDE PUFFING.')
            call prc(sp//'EACH SUBSEQUENT ITERATION WILL RE-LOAD THE')
            call prc(sp//'CORE PLASMA SOLUTION (AND THE PRIVATE'//' PLASMA SOLUTION')
            call prc(sp//'IF SPECIFIED) FROM THE ORIGINAL INPUT. ONLY')

!         End of Seed plasma options

            call prc(sp//'THE SOL SOLUTION WILL BE ITERATED.')

!       End of E2dstart

          endif

!     End of Ionization Option 1,2,8

        endif

!     Private plasma ionization option - if applicable


      endif

!     Only if the TGRAD option is turned ON.

      if (ciopto.eq.1) then
        call prb
        CALL PRC(S1//'PRIVATE PLASMA IONIZATION OPTION:')
        call prc(s1//'- NOTE.1: If the main Ionization Option has been'//' specified as')
        call prc(s1//'  PIN iterated (e.g. Options 1,2,8). Then the'//' value')
        call prc(s1//'  specified here is ONLY used for the seed'//' plasma solution.')
        call prc(s1//'  The PIN iterated ionization option'//' will be used')
        call prc(s1//'  the rest of the time - except in the case'//' of option -2')
        call prc(s1//'  which will use the prescribed private plasma.')
        call prc(s1//'- NOTE.2: The Edge2D seed plasma option takes')


! slmod begin - new
        call prc(s1//'  prcedence over option -2.')
        if     (switch(swionp).eq.-6.0) then
          CALL PRC(S1//'PP IONIZ OPT  -6: A UNIFORM PLASMA IS ASSIGNED')
          CALL PRC(SP//'TO EACH RING IN THE PRIVATE FLUX ZONE FROM A')
          CALL PRC(SP//'LISTING OF TEMPERATURE AND DENSITY IN THE')
          CALL PRC(SP//'DIVIMP INPUT FILE.  THE TARGET VALUES ARE SET')
          CALL PRC(SP//'FROM THE BULK PLASMA VALUES.')
          fp = 7
          WRITE(fp,1158) 'Ring','Region','Te (eV)','Ti (eV)','n (m-3)','v (m s-1)'
          DO i1 = 1, osmnppv
            WRITE(fp,1159) (INT(osmppv(i1,i2)),i2=1,2),(osmppv(i1,i2) ,i2=3,6)
          ENDDO
        elseif (switch(swionp).eq.-5.0) then
          CALL PRC(S1//'PP IONIZ OPT  -5: A UNIFORM PLASMA IS ASSIGNED')
          CALL PRC(SP//'TO EACH RING IN THE PRIVATE FLUX ZONE FROM A')
          CALL PRC(SP//'LISTING OF TEMPERATURE AND DENSITY IN THE')
          CALL PRC(SP//'DIVIMP INPUT FILE.  THE TARGET VALUES ARE NOT')
          CALL PRC(SP//'MODIFIED.')
          fp = 7
          WRITE(fp,1158) 'Ring','Region','Te (eV)','Ti (eV)','n (m-3)','v (m s-1)'
          DO i1 = 1, osmnppv
            WRITE(fp,1159) (INT(osmppv(i1,i2)),i2=1,2),(osmppv(i1,i2) ,i2=3,6)
          ENDDO
1158      FORMAT(5X,2A8,2A10,2A12)
1159      FORMAT(5X,2I8,2F10.2,1P,2E12.2,0P)
        elseif (switch(swionp).eq.-4.0) then
          CALL PRC(S1//'PP IONIZ OPT  -4: THOMSON DATA APPLIED TO PP')
          CALL PRC(SP//'AVERAGE VALUE OF TH DATA ON A RING IS ASSIGNED')
          CALL PRC(SP//'TO EVERY CELL ON THE RING.  RINGS WITHOUT ')
          CALL PRC(SP//'DATA ARE INTERPOLATED.  THE TARGET FLUX IS')
          CALL PRC(SP//'ASSIGNED USING THE THOMSON DATA.')
        elseif (switch(swionp).eq.-3.0) then
          CALL PRC(S1//'PP IONIZ OPT  -3: THOMSON DATA APPLIED TO PP')
          CALL PRC(SP//'AVERAGE VALUE OF TH DATA ON A RING IS ASSIGNED')
          CALL PRC(SP//'TO EVERY CELL ON THE RING.  RINGS WITHOUT ')
          CALL PRC(SP//'DATA ARE INTERPOLATED.  THE TARGET FLUX IS')
          CALL PRC(SP//'SPECIFIED IN THE DIVIMP INPUT FILE.')

!        if (switch(swionp).eq.-2) then
! slmod end

        elseif (switch(swionp).eq.-2.0) then
          CALL PRC(S1//'PP IONIZ OPT  -2: ALL PP CONDITIONS'//' ARE PRESCRIBED')

          CALL PRC(SP//'BY THE FOLLOWING FACTORS')

          CALL PRC(SP//'ALL DISTANCES ARE PROPORTIONS OF SMAX'//' FOR THE RING')
          CALL PRC(SP//'QUANTITIES ARE LINEARLY INTERPOLATED'//' BETWEEN POINTS')

          CALL PRC(SP//'AND CONSTANT OUTSIDE RANGE AT END-POINT VALUES')
          write(coment,1160)  ctes1,ctef1,ctes2,ctef2
          call prc(coment)
          write(coment,1170)  ctis1,ctif1,ctis2,ctif2
          call prc(coment)
          write(coment,1180)  cnes1,cnef1,cnes2,cnef2
          call prc(coment)
          write(coment,1190)  cvbs1,cvbf1,cvbs2,cvbf2

          call prc(coment)
 1160     format(6x,'  (S,Te) : (0.0,Tet) -> (',f6.3,',',f6.3,'*Tet) -> (',f6.3,',',f6.3,'*Tet)')
 1170     format(6x,'  (S,Ti) : (0.0,Tit) -> (',f6.3,',',f6.3,'*Tit) -> (',f6.3,',',f6.3,'*Tit)')
 1180     format(6x,'  (S,Ne) : (0.0,Net) -> (',f6.3,',',f6.3,'*Net) -> (',f6.3,',',f6.3,'*Net)')

 1190     format(6x,'  (S,Vb) : (0.0,Vbt) -> (',f6.3,',',f6.3,'*Vbt) -> (',f6.3,',',f6.3,'*Vbt)')

        elseif (switch(swionp).eq.-1) then
          CALL PRC(S1//'PP IONIZ OPT  -1: PP IONIZATION OPTIONS'//' ARE THE SAME')

          CALL PRC(SP//'AS FOR THE MAIN SOL.')

        elseif (switch(swionp).eq.0.0) then
          call prb
          CALL PRC(S1//'PP IONIZ OPT   0: EXPONENTIAL'//' DECAY IONIZATION SOURCE')

          CALL PRB
          CALL PRQ (SP//'LENGTH OF IONIZATION SOURCE          ',LENSFI)

          CALL PRQ (SP//'DECAY LENGTH OF IONIZATION SOURCE    ',LAMS)

        elseif (switch(swionp).eq.3.0) then
          CALL PRC(S1//'PP IONIZ OPT   3: IMPOSED'//' TRIANGULAR IONIZATION SOURCE')
          CALL PRQ(SP//'EXTENDING FROM   :',LENSST)
          CALL PRQ(SP//'          TO     :',LENSFI)
          CALL PRC(SP//'INTEGRAL OF SOURCE NORMALIZED'//' TO RING TARGET FLUX.')

          CALL PRB

        elseif (switch(swionp).eq.4.0) then
          CALL PRC(S1//'PP IONIZ OPT   4: IMPOSED'//' RECTANGULAR IONIZATION SOURCE')
          CALL PRQ(SP//'EXTENDING FROM   :',LENSST)
          CALL PRQ(SP//'          TO     :',LENSFI)
          CALL PRC(SP//'INTEGRAL OF SOURCE NORMALIZED'//' TO RING TARGET FLUX.')

          CALL PRB

        elseif (switch(swionp).eq.5.0) then
          CALL PRC(S1//'PP IONIZ OPT   5: ALGORITHMIC'//' RECT/TRI IONIZATION'//' SOURCE')
          CALL PRC(SP//'IF NT > 1.0E19  - TRIANGULAR SOURCE')
          CALL PRC(SP//'   FROM L1= (13 - 10*TET) M (TET < '//' 1.3 EV)')
          CALL PRC(SP//'     OR L1= 0.0           M (TET >='//' 1.3 EV)')
          CALL PRC(SP//'     TO L2=L1+2 (M)')
          CALL PRC(SP//'IF NT <= 1.0E19 - RECTANGULAR SOURCE')
          CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 13-TET (M)'//' (TET<10EV)')
          CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 2 (M)'//' (TET>10EV)')

          call prb

        elseif (switch(swionp).eq.6.0) then
          CALL PRC(S1//'PP IONIZ OPT   6: IMPOSED'//' S**5 GAUSSIAN IONIZATION'//' SOURCE')
          CALL PRC(SP//'OF FORM:  A * S**5 * EXP(-ALPHA'//' * S**2)')
          CALL PRC(SP//'EXTENDING FROM : 0.0 (M)')
          CALL PRQ(SP//'   CUT OFF AT  :',LENSFI)
          CALL PRC(SP//'WHERE          : ALPHA = 2.5 / WF**2')
          CALL PRQ(SP//'WITH WIDTH FACTOR (WF):',LAMS)
          CALL PRC(SP//'INTEGRAL OF SOURCE NORMALIZED'//' TO RING TARGET FLUX.')

          call prb

        elseif (switch(swionp).eq.7.0) then
          CALL PRC(S1//'PP IONIZ OPT   7: ALGORITHMIC'//' RECT/S**5GAUSS IONIZATION'//' SOURCE - ON')
          CALL PRC(SP//'IF NT > 1.0E19  - S5 GAUSSIAN SOURCE')
          CALL PRC(SP//'   WITH WIDTH FACTOR WF = (14-10TET)'//' (M) (TET <  1.3EV)')
          CALL PRC(SP//'   WITH WIDTH FACTOR WF = 1.0       '//' (M) (TET >= 1.3EV)')
          CALL PRC(SP//'   FROM L1= 0.0 TO L2 = 1/2 RING LENGTH')
          CALL PRC(SP//'IF NT <= 1.0E19 - RECTANGULAR SOURCE')
          CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 13-TET (M)'//' (TET<10EV)')
          CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 2 (M)'//' (TET>10EV)')

          call prb

        elseif (switch(swionp).eq.9.0) then
          CALL PRC(S1//'PP IONIZ OPT   9: IMPOSED'//' OFFSET S**5 GAUSSIAN'//' IONIZATION SOURCE')
          CALL PRC(SP//'OF FORM:  A * (S+L)**5 *'//' EXP(-ALPHA * (S+L)**2)')
          CALL PRC(SP//'EXTENDING FROM : 0.0 (M)')
          CALL PRQ(SP//'   CUT OFF AT  :',LENSFI)
          CALL PRC(SP//'WHERE          : ALPHA = 2.5 / WF**2')
          CALL PRQ(SP//'WITH WIDTH FACTOR (WF):',LAMS)
          CALL PRC(SP//'AND L = WF/2.0')
          CALL PRC(SP//'INTEGRAL OF SOURCE NORMALIZED'//' TO RING TARGET FLUX.')

          call prb

        elseif (switch(swionp).eq.10.0) then
          CALL PRC(S1//'PP IONIZ OPT  10: ALGORITHMIC'//' RECT/OFFSET S**5GAUSS'//' IONIZATION SOURCE')
          CALL PRC(SP//'IF NT > 1.0E19  - OFFSET S5'//' GAUSSIAN SOURCE')
          CALL PRC(SP//'   WITH WIDTH FACTOR WF = (28-20TET)'//' (M) (TET <  1.3EV)')
          CALL PRC(SP//'   WITH WIDTH FACTOR WF = 2.0       '//' (M) (TET >= 1.3EV)')
          CALL PRC(SP//'   FROM L1= 0.0 TO L2 = 1/2 RING LENGTH')
          CALL PRC(SP//'IF NT <= 1.0E19 - RECTANGULAR SOURCE')
          CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 13-TET (M)'//' (TET<10EV)')
          CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 2 (M)'//' (TET>10EV)')

          call prb

        elseif (switch(swionp).eq.11.0) then

          CALL PRC(s1//'PP IONIZ OPT  11: IONIZATION'//' SOURCE DATA READ FROM'//' EDGE2D INPUT FOR CASE.')

        elseif (switch(swionp).eq.12.0) then
          CALL PRC(S1//'PP IONIZ OPT  12: PIN IS RUN'//' WITH AN EDGE2D BACKGROUND')
          call prc(sp//'IN THE SOL ONCE BEFORE SOL 22 IS INVOKED.')

          CALL PRC(SP//'THIS GENERATES POWER TERMS BUT WILL'//' NOT INCLUDE PUFFING.')

        elseif (switch(swionp).eq.13.0) then
          CALL PRC(S1//'PP IONIZ OPT  13: PIN IS RUN'//' WITH AN EDGE2D BACKGROUND')
          call prc(sp//'IN THE SOL ONCE BEFORE SOL 22 IS INVOKED.')

          CALL PRC(SP//'THIS GENERATES POWER TERMS AND'//' PROPER PUFFING.')

        elseif (switch(swionp).eq.14.0) then
          CALL PRC(S1//'PP IONIZ OPT  14: PIN IS RUN'//' WITH AN EDGE2D BACKGROUND')
          call prc(sp//'IN THE SOL ONCE BEFORE SOL 22 IS INVOKED.')
          CALL PRC(SP//'THIS GENERATES POWER TERMS AND'//' PROPER PUFFING.')

          CALL PRC(SP//'THIS SHOULD BE USED'//' IN CONJUNCTION'//' WITH CORE OPTION -1.')

        elseif (switch(swionp).eq.15.0) then
          CALL PRC(s1//'PP IONIZ OPT  15: IONIZATION'//' SOURCE DATA READ FROM'//' EDGE2D INPUT FOR CASE.')

          call prc(sp//'EDGE2D PLASMA ASSIGNED AS OLD'//' FOR ENTIRE BACKGROUND')

!     CIOPTO = 0

        endif

      elseif (ciopto.eq.0.or.ciopto.eq.2.or.ciopto.eq.3.or.ciopto.eq.4) then
         CALL PRC(S1//'PRIVATE PLASMA IS NOT CALCULATED USING SOL 22.')
         CALL PRC(S1//'PRIVATE PLASMA WILL BE CONSTANT'//' AT TARGET VALUES')

!     ENDIF for CIOPTO

         CALL PRC(S1//'IF NO OTHER PRIVATE PLASMA OPTION HAS BEEN'//' SPECIFIED.')

      endif

!     Convection Term number 1

      call prb
      if (switch(swcond).eq.0.0) then
         CALL PRC(S1//'5/2 NV*kT OPT 0 : THIS TERM IF OFF')
      else
         CALL PRC(S1//'5/2 NV*kT OPT 1 : THIS TERM IF ON')

      endif

!     Convection Term number 2

      call prb
      if (switch(swconv).eq.0.0) then
         CALL PRC(S1//'1/2 mV^3N OPT 0 : THIS TERM IS OFF')
      else
         CALL PRC(S1//'1/2 mV^3N OPT 1 : THIS TERM IS ON')

!     Prad term

      endif

      call prb
      if (switch(swprad).eq.0.0) then
         CALL PRC(S1//'PRAD OPTION   0 : PRAD TERM IS OFF')
      elseif (switch(swprad).eq.1.0) then
         CALL PRC(S1//'PRAD OPTION   1 : EXPONENTIAL'//' DECAY RADIATION SOURCE')
         CALL PRQ(SP//'LENGTH OF RADIATION SOURCE         ', LENR)
         CALL PRQ(SP//'DECAY LENGTH OF RADIATION  SOURCE  ', LAMR)
         CALL PRQ(SP//'SOURCE STRENGTH FRACTION (FRR)     ', FRR)
      elseif (switch(swprad).eq.2.0) then
         CALL PRC(S1//'PRAD OPTION   2 : ANALYTIC'//' (GARCHING) RADIATION SOURCE')
         CALL PRQ(SP//'ALPHA = NIMP/NE      ',ALFIMP)
         CALL PRQ(SP//'TBASE FACTOR         ',TALIMP)
         CALL PRQ(SP//'FIRST EXPONENT       ',EX1IMP)
         CALL PRQ(SP//'SECOND EXPONENT      ',EX2IMP)
      elseif (switch(swprad).eq.3.0) then
         CALL PRC(S1//'PRAD OPTION   3 : PROPORTIONAL TO PINQE')
         CALL PRC(SP//'IMPURITY RADIATION IS A GIVEN FACTOR')
         CALL PRC(SP//'TIMES THE CALCULATED PINQE.')
         CALL PRQ(SP//'MULTIPLICATION FACTOR = ',radsrc_mult)
      elseif (switch(swprad).eq.4.0) then
         CALL PRC(S1//'PRAD OPTION   4 : EXTERNAL RADIATION SOURCE')
         CALL PRC(SP//'IMPURITY RADIATION IS READ IN FROM A FILE')
         CALL PRC(SP//'USUALLY RESULTING FROM A PREVIOUS DIVIMP')
         CALL PRC(SP//'CASE.')
         CALL PRQ(SP//'MULTIPLICATION FACTOR = ',radsrc_mult)
      elseif (switch(swprad).eq.5.0) then 
         CALL PRC(S1//'PRAD OPTION   5 : EXTERNAL RADIATION SOURCE')
         call prc(s1//' **** WARNING ****')
         call prc(s1//' IN DEVELOPMENT - NOT YET IMPLEMENTED')
         call prc(s1//' **** WARNING ****')
         CALL PRC(SP//'TOTAL RADIATION TO BE APPLIED IN SOL22 IS')
         CALL PRC(SP//'SPECIFIED IN THE INPUT ON A PER REGION BASIS')
         CALL PRC(SP//'BASED ON EXPERIMENTAL BOLOMETRY DATA')
         CALL PRC(SP//'THIS IS THEN BE DISTRIBUTED TO EACH RING AND')
!         call prq(sp//'REGION N ', prad_for_region)
         CALL PRC(SP//'CELL USING VARIOUS MECHANISMS')
      elseif (switch(swprad).eq.6.0) then
         CALL PRC(S1//'PRAD OPTION   6 : RECTANGULAR'//' RADIATION SOURCE')
         CALL PRQ(SP//'END LENGTH OF RADIATION SOURCE     ', LENR)
         CALL PRQ(SP//'START LENGTH OF RADIATION  SOURCE  ', LAMR)
         CALL PRQ(SP//'SOURCE STRENGTH FRACTION (FRR)     ', FRR)

      endif

!     Phelpi term

      call prb
      if (switch(swphelp).eq.0.0) then
         CALL PRC(S1//'PHELPI OPTION 0 : PHELPI TERM  IS OFF')
      elseif (switch(swphelp).eq.1.0) then
         CALL PRC(S1//'PHELPI OPTION 1 : INTERNAL'//' PHELPI TERM  IS ON')
         CALL PRC(SP//'SEE MANUAL FOR DETAILS.')
      elseif (switch(swphelp).eq.2.0) then
         CALL PRC(S1//'PHELPI OPTION 2 : PINQE TERM  IS ON')
         CALL PRC(SP//'INTERNAL TERM IS ON FOR SEED PLASMA')
         CALL PRC(SP//'SEE MANUAL FOR DETAILS.')
      elseif (switch(swphelp).eq.3.0) then
         CALL PRC(S1//'PHELPI OPTION 3 : PINQE TERM  IS ON')
         CALL PRC(SP//'INTERNAL PHELPI TERM IS OFF FOR SEED PLASMA')

      endif
      if (tcutqe.gt.0.0) then
         CALL PRB
         CALL PRQ(S1//'QE TEMPERATURE CUTOFF (EV) = ',TCUTQE)
         CALL PRC(S1//'- ALL CONTRIBUTIONS TO PHELPI TERM')
         call prc(s1//'  FROM REGIONS WITH TE<TCUT HAVE BEEN SET')
         call prc(s1//'  TO ZERO.')

      endif

      !     Pei Term

      call prb
      if (switch(swpei).eq.0.0) then
         CALL PRC(S1//'PEI OPTION    0 : PEI TERM  IS OFF')
      elseif (switch(swpei).eq.1.0) then
         CALL PRC(S1//'PEI OPTION    1 : PEI TERM  IS ON')
         CALL PRQ(SP//'PEI CORRECTION FACTOR =',PEICF)
         CALL PRC(SP//'SEE MANUAL FOR DETAILS.')
      elseif (switch(swpei).eq.3.0) then
         CALL PRC(S1//'PEI OPTION    3 : PEI TERM  IS'//' CALCULATED')
         CALL PRC(SP//'BUT NOT USED.')

      endif

!     PCX term

      call prb
      if (switch(swpcx).eq.0.0) then
         CALL PRC(S1//'PCX OPTION    0 : PCX TERM  IS OFF')
      elseif (switch(swpcx).eq.1.0) then
         CALL PRC(S1//'PCX OPTION    1 : INTERNAL PCX TERM  IS ON')
         CALL PRQ(SP//'CX POWER COEFF (CEICF) =',CEICF)
         CALL PRC(SP//'SEE MANUAL FOR DETAILS.')
      elseif (switch(swpcx).eq.2.0) then
         CALL PRC(S1//'PCX OPTION    2 : PINQI TERM IS ON')
         CALL PRC(SP//'INTERNAL TERM IS ON FOR SEED PLASMA')
         CALL PRQ(SP//'CEICF FOR SEED PLASMA',CEICF)
         CALL PRC(SP//'SEE MANUAL FOR DETAILS.')
      elseif (switch(swpcx).eq.3.0) then
         CALL PRC(S1//'PCX OPTION    3 : PINQI TERM IS ON')
         CALL PRC(SP//'INTERNAL PCX TERM IS OFF FOR SEED PLASMA')
      elseif (switch(swpcx).eq.4.0) then
         CALL PRC(S1//'PCX OPTION    4 : DIVIMP QI TERM  IS ON')
         CALL PRC(SP//'INTERNAL TERM IS OFF FOR SEED PLASMA')
         call prb
         CALL PRC(S2//'THE FOUR COMPONENTS OF DIVIMP QI'//' ARE CALCULATED USING')
         CALL PRC(S2//'THE FOLLOWING SUB-OPTIONS')

!        Add switch print-outs for swqidatiz,swqidmliz,swqidcx,swqidrec

         call prb

         CALL PRC(S2//'PINQID SUB-OPTION LISTING:')

         call prb
         if (switch(swqidatiz).eq.0.0) then
            CALL PRC(S2//'ATOMIC IONIZATION OPT 0 : OFF')
         elseif (switch(swqidatiz).eq.1.0) then
            CALL PRC(S2//'ATOMIC IONIZATION OPT 1 : ON')
            CALL PRC(S2//' - 3/2 K TATOM * (NH/(NH+NH2)) * SIONIZ')
            CALL PRQ(S2//' - CUTOFF TEMPERATURE (EV) = ',TCUTATIZ)
            CALL PRC(S2//' - CONTRIBUTIONS = 0 IF TI < TCUT')
         elseif (switch(swqidatiz).eq.2.0) then
            CALL PRC(S2//'ATOMIC IONIZATION OPT 2 : ON')
            CALL PRC(S2//' - 3/2 K TATOM * NE * NH * <SV>ATIZ (ADAS)')
            CALL PRQ(S2//' - CUTOFF TEMPERATURE (EV) = ',TCUTATIZ)
            CALL PRC(S2//' - CONTRIBUTIONS = 0 IF TI < TCUT')

         endif
         if (switch(swqidmliz).eq.0.0) then
            CALL PRC(S2//'MOLECULAR IONIZATION OPT 0 : OFF')
         elseif (switch(swqidmliz).eq.1.0) then
            CALL PRC(S2//'MOLECULAR IONIZATION OPT 1 : ON')
            CALL PRC(S2//' - 3.0 * (NH2/(NH+NH2)) * SIONIZ')
            CALL PRQ(S2//' - CUTOFF TEMPERATURE (EV) = ',TCUTMLIZ)
            CALL PRC(S2//' - CONTRIBUTIONS = 0 IF TI < TCUT')
         elseif (switch(swqidmliz).eq.2.0) then
            CALL PRC(S2//'MOLECULAR IONIZATION OPT 2 : ON')
            CALL PRC(S2//' - 3.0 * NE * NH2 * <SV>ATIZ (ADAS)')
            CALL PRC(S2//' - USES ATOMIC IONIZATION AS AN')
            CALL PRC(S2//'   APPROXIMATION TO MOLECULAR IONZIATION.')
            CALL PRQ(S2//' - CUTOFF TEMPERATURE (EV) = ',TCUTMLIZ)
            CALL PRC(S2//' - CONTRIBUTIONS = 0 IF TI < TCUT')

         endif
         if (switch(swqidrec).eq.0.0) then
            CALL PRC(S2//'RECOMBINATION OPTION 0 : OFF')
         elseif (switch(swqidrec).eq.1.0) then
            CALL PRC(S2//'RECOMBINATION OPTION 1 : ON')
            CALL PRC(S2//' - -3/2 * K TI * NI * NE * <SV>REC (ADAS)')
            CALL PRQ(S2//' - CUTOFF TEMPERATURE (EV) = ',TCUTREC)
            CALL PRC(S2//' - CONTRIBUTIONS = 0 IF TI < TCUT')

         endif
         if (switch(swqidcx).eq.0.0) then
            CALL PRC(S2//'CHARGE EXCHANGE OPTION 0 : OFF')
         elseif (switch(swqidcx).eq.1.0) then
            CALL PRC(S2//'CHARGE EXCHANGE OPTION 1 : ON')
            CALL PRC(S2//' - 3/2 K TREF * RCXMULT(TI) * CEICF *'//' SIONIZ')
            CALL PRQ(S2//' - TREF (EV) = ', TREFCX)
            CALL PRQ(S2//' - CEICF     = ',CEICF)
            CALL PRQ(S2//' - CUTOFF TEMPERATURE (EV) = ',TCUTCX)
            CALL PRC(S2//' - CONTRIBUTIONS = 0 IF TI < TCUT')
         elseif (switch(swqidcx).eq.2.0) then
            CALL PRC(S2//'CHARGE EXCHANGE OPTION 2 : ON')
            CALL PRC(S2//' - 3/2 * (TATOM - TI) * NE * NH * <SV>CX')
            CALL PRC(S2//' - DEFINE EAV = 3/2 K (TATOM+ TI) / 2.0')
            CALL PRC(S2//' - <SV>CX = 10E-13           '//'  FOR EAV  > 1000 EV')
            CALL PRC(S2//' - <SV>CX = 10E-14 EAV**(1/3)'//'  FOR EAV =< 1000 EV')
            CALL PRQ(S2//' - CUTOFF TEMPERATURE (EV) = ',TCUTCX)
            CALL PRC(S2//' - CONTRIBUTIONS = 0 IF TI < TCUT')
         elseif (switch(swqidcx).eq.3.0) then
            CALL PRC(S2//'CHARGE EXCHANGE OPTION 3 : ON')
            CALL PRC(S2//' - 3/2 * (TATOM - TI) * NE * NH *'//' <SV>CX (ADAS)')
            CALL PRC(S2//' - ADAS CX DATA IS CURRENTLY UNRELIABLE')
            CALL PRC(S2//'   CHECK BEFORE USING THIS OPTION.')
            CALL PRQ(S2//' - CUTOFF TEMPERATURE (EV) = ',TCUTCX)
            CALL PRC(S2//' - CONTRIBUTIONS = 0 IF TI < TCUT')
         elseif (switch(swqidcx).eq.4.0) then
            CALL PRC(S2//'CHARGE EXCHANGE OPTION 4 : ON')
            CALL PRC(S2//' - 3/2 * (TATOM - TI) * NE * NH * <SV>CX')
            CALL PRC(S2//' - DEFINE EAV = 3/2 K (TATOM+ TI) / 2.0')
            CALL PRC(S2//' - <SV>CX = 10E-13           '//'  FOR EAV  > 1000 EV')
            CALL PRC(S2//' - <SV>CX = 10E-14 EAV**(1/3)'//'  FOR EAV =< 1000 EV')
            CALL PRQ(S2//' - CUTOFF TEMPERATURE (EV) = ',TCUTCX)
            CALL PRC(S2//' - CONTRIBUTIONS = 0 IF TI < TCUT')
            CALL PRC(S2//' - TERM IS LIMITED TO ION COOLING ONLY')

         endif
      elseif (switch(swpcx).eq.5.0) then
         CALL PRC(S1//'PCX OPTION    5 : PINQI TERM IS ON')
         call prc(sp//'PINQI TERM IS CLIPPED - IT IS NOT ALLOWED')
         call prc(sp//'TO HEAT THE PLASMA.')
         CALL PRC(SP//'INTERNAL TERM IS ON FOR SEED PLASMA')
         CALL PRQ(SP//'CEICF FOR SEED PLASMA',CEICF)
         CALL PRC(SP//'SEE MANUAL FOR DETAILS.')

      endif
      if (switch(swpcx).ne.4.0.and.tcutcx.gt.0.0) then
         call prb
         CALL PRQ(S1//'QI TEMPERATURE CUTOFF (EV) = ',TCUTCX)
         CALL PRC(S1//'- ALL CONTRIBUTIONS FROM REGIONS WITH TI<TCUT')
         CALL PRC(S1//'  HAVE BEEN SET TO ZERO.')

      endif

!     External Epower term

      if (switch(swepow).eq.0.0) then
         CALL PRC(S1//'EXT EPOW OPTION   0 : EXTERNAL Epower TERM  IS OFF')
      elseif (switch(swepow).eq.1.0) then
         CALL PRC(S1//'EXT EPOW OPTION   1 : EXTERNAL Epower TERM  IS ON')
         CALL PRC(SP//'DATA IS READ FROM THE DIV AUX INPUT FILE WITH TAG "EXTEPOW:"')
      elseif (switch(swepow).eq.2.0) then
         CALL PRC(S1//'EXT EPOW OPTION   2 : EXTERNAL Epower TERM  IS ON')
         CALL PRC(SP//'DATA IS INTERPOLATED FROM EXTERNAL R,Z DATA PROVIDED IN FILE: '//trim(ext_epow_fn))
      endif
      

!     External Ipower term

      if (switch(swipow).eq.0.0) then
         CALL PRC(S1//'EXT IPOW OPTION   0 : EXTERNAL Ipower TERM  IS OFF')
      elseif (switch(swipow).eq.1.0) then
         CALL PRC(S1//'EXT IPOW OPTION   1 : EXTERNAL Ipower TERM  IS ON')
         CALL PRC(SP//'DATA IS READ FROM THE DIV AUX INPUT FILE WITH TAG "EXTIPOW:"')
      elseif (switch(swipow).eq.2.0) then
         CALL PRC(S1//'EXT IPOW OPTION   2 : EXTERNAL Ipower TERM  IS ON')
         CALL PRC(SP//'DATA IS INTERPOLATED FROM EXTERNAL R,Z DATA PROVIDED IN FILE: '//trim(ext_ipow_fn))
      endif

!     Private plasma electron TARGET power LOSS compensation term
      
      call prb
      if (switch(swppelec).eq.0.0) then
         CALL PRC(S1//'PP ELEC POW OPTION 0: PP ELECTRON POWER'//' LOSS COMPENSATION TERM IS OFF')
      elseif (switch(swppelec).eq.1.0) then
         CALL PRC(S1//'PP ELEC POW OPTION 1: PP ELECTRON POWER'//' LOSS COMPENSATION TERM IS ON')
         call prc(sp//'ELECTRON POWER LOST TO EACH ELEMENT OF PP')
         call prc(sp//'TARGET IS REMOVED FROM MAIN SOL RINGS IN')
         call prc(sp//'THE CORRESPONDING POSITION RELATIVE TO THE')
         CALL PRC(SP//'SEPARATRIX. THIS POWER LOSS')
         call prc(sp//'IS DISTRTIBUTED EVENLY TO THE XPOINT ON')
         call prc(sp//'THE MAIN SOL RINGS')
      elseif (switch(swppelec).eq.2.0) then
         CALL PRC(S1//'PP ELEC POW OPTION 2: PP ELECTRON POWER'//' LOSS COMPENSATION TERM IS ON')
         call prc(sp//'ELECTRON POWER LOST TO EACH ELEMENT OF PP')
         call prc(sp//'TARGET IS REMOVED FROM MAIN SOL RINGS IN')
         call prc(sp//'THE CORRESPONDING POSITION RELATIVE TO THE')
         CALL PRC(SP//'SEPARATRIX. THIS POWER LOSS')
         call prq(sp//'IS DISTRTIBUTED EVENLY OVER SMAX *',pp_pow_dist)
         call prc(sp//'FROM EACH TARGET ON THE MAIN SOL RINGS.')
      elseif (switch(swppelec).eq.3.0) then
         CALL PRC(S1//'PP ELEC POW OPTION 3: PP ELECTRON POWER'//' LOSS COMPENSATION TERM IS ON')
         call prc(sp//'ION POWER LOST TO LOCAL PFZ TARGET')
         call prc(sp//'IS DISTRIBUTED TO MAIN SOL RINGS USING')
         call prc(sp//'ONE OF 5 POSSIBLE RING DISTRIBUTION OPTIONS')
         call prc(sp//'IT IS DISTRTIBUTED EVENLY TO THE XPOINT ON')
         call prc(sp//'THE MAIN SOL RINGS')

!     Private plasma ION TARGET power LOSS compensation term

      ENDIF

      call prb
      if (switch(swppion).eq.0.0) then
         CALL PRC(S1//'PP ION POW OPTION 0 : PP ION POWER'//' LOSS COMPENSATION TERM IS OFF')
      elseif (switch(swppion).eq.1.0) then
         CALL PRC(S1//'PP ION POW OPTION 1 : PP ION POWER'//' LOSS COMPENSATION TERM IS ON')
         call prc(sp//'ION POWER LOST TO EACH ELEMENT OF PP')
         call prc(sp//'TARGET IS REMOVED FROM MAIN SOL RINGS IN')
         call prc(sp//'THE CORRESPONDING POSITION RELATIVE TO THE')
         CALL PRC(SP//'SEPARATRIX. THIS POWER LOSS')
         call prc(sp//'IS DISTRTIBUTED EVENLY TO THE XPOINT ON')
         call prc(sp//'THE MAIN SOL RINGS')
      elseif (switch(swppion).eq.2.0) then
         CALL PRC(S1//'PP ION POW OPTION 2 : PP ION POWER'//' LOSS COMPENSATION TERM IS ON')
         call prc(sp//'ION POWER LOST TO EACH ELEMENT OF PP')
         call prc(sp//'TARGET IS REMOVED FROM MAIN SOL RINGS IN')
         call prc(sp//'THE CORRESPONDING POSITION RELATIVE TO THE')
         CALL PRC(SP//'SEPARATRIX. THIS POWER LOSS')
         call prq(sp//'IS DISTRTIBUTED EVENLY OVER SMAX *',pp_pow_dist)
         call prc(sp//'FROM EACH TARGET ON THE MAIN SOL RINGS.')
      elseif (switch(swppion).eq.3.0) then
         CALL PRC(S1//'PP ION POW OPTION 3 : PP ION POWER'//' LOSS COMPENSATION TERM IS ON')
         call prc(sp//'ION POWER LOST TO LOCAL PFZ TARGET')
         call prc(sp//'IS DISTRIBUTED TO MAIN SOL RINGS USING')
         call prc(sp//'ONE OF 5 POSSIBLE RING DISTRIBUTION OPTIONS')
         call prc(sp//'IT IS DISTRTIBUTED EVENLY TO THE XPOINT ON')
         call prc(sp//'THE MAIN SOL RINGS')

!     Private plasma PRESSURE LOSS compensation term

      ENDIF

      call prb
      if (switch(swppress).eq.0.0) then
         CALL PRC(S1//'PP PRESSURE OPTION 0 : PP PRESSURE'//' LOSS COMPENSATION TERM IS OFF')
      elseif (switch(swppress).eq.1.0) then
         CALL PRC(S1//'PP PRESSURE OPTION 1 : PP PRESSURE'//' LOSS COMPENSATION TERM IS ON')
         call prc(sp//'PRESSURE LOST TO EACH ELEMENT OF PP')
         call prc(sp//'TARGET IS ADDED TO MAIN SOL RINGS IN')
         call prc(sp//'THE CORRESPONDING POSITION RELATIVE TO THE')
         CALL PRC(SP//'SEPARATRIX. THIS PRESSURE')
         call prc(sp//'IS DISTRTIBUTED EVENLY TO THE XPOINT ON')
         call prc(sp//'THE MAIN SOL RINGS')
      elseif (switch(swppress).eq.2.0) then
         CALL PRC(S1//'PP PRESSURE OPTION 2 : PP PRESSURE'//' LOSS COMPENSATION TERM IS ON')
         call prc(sp//'PRESSURE LOST TO EACH ELEMENT OF PP')
         call prc(sp//'TARGET IS ADDED TO MAIN SOL RINGS IN')
         call prc(sp//'THE CORRESPONDING POSITION RELATIVE TO THE')
         CALL PRC(SP//'SEPARATRIX. THIS PRESSURE')
         call prq(sp//'IS DISTRTIBUTED EVENLY OVER SMAX *',pp_pow_dist)
         call prc(sp//'FROM EACH TARGET ON THE MAIN SOL RINGS.')
      elseif (switch(swppress).eq.3.0) then
         CALL PRC(S1//'PP PRESSURE OPTION 3 : PP PRESSURE'//' LOSS COMPENSATION TERM IS ON')
         call prc(sp//'PRESSURE LOST TO LOCAL PFZ TARGET')
         call prc(sp//'IS DISTRIBUTED TO MAIN SOL RINGS USING')
         call prc(sp//'ONE OF 5 POSSIBLE RING DISTRIBUTION OPTIONS')
         call prc(sp//'IS DISTRTIBUTED EVENLY TO THE XPOINT ON')
         call prc(sp//'THE MAIN SOL RINGS')

!     Viscosity term

      ENDIF
      if (switch(swvisc1).ne.0.0) then
         call prb
         CALL PRC(S1//'DENSITY VISCOSITY CORRECTION TERM'//'  IS NOT IMPLEMENTED')

!     Momentum Loss term

      endif

      call prb
      if (switch(swnmom).eq.0.0) then
         CALL PRC(S1//'MOMENTUM OPT  0 : MOMENTUM LOSS TERM IS OFF')
      elseif (switch(swnmom).eq.1.0) then
         CALL PRC(S1//'MOMENTUM OPT  1 : MOMENTUM LOSS TERM IS ON')
         CALL PRC(SP//'SMOM = SMOM0    S < L * SMAX        ')
         CALL PRC(SP//'     = 0        S > L * SMAX        ')
         CALL PRC(SP//'SMOM0= PT/(L * SMAX) * (1/FFRIC -1)')
         if (n_extffric.eq.0) then 
            CALL PRQ(SP//'FFRIC= ',FFRIC)
            CALL PRQ(SP//'L    = ',LENMOM)
         else
            CALL PRQ(SP//'DEFAULT FFRIC= ',FFRIC)
            CALL PRQ(SP//'DEFAULT L    = ',LENMOM)
            call prc(sp//'OVER-RIDDEN BY THE'//' FOLLOWING VALUES ON SPECIFIED RINGS:')
            call prc(sp//'         FFRIC  LENMOM    FFRIC LENMOM')
            call prc(sp//'RING    '//OUTER//'      '//INNER)
            do in = 1, n_extffric
               write(coment,'(i4,2x,4(1x,f8.3))')int(extffric(in,1)),extffric(in,2),extffric(in,3),extffric(in,4),extffric(in,5)
               call prc(sp//coment)
            end do 
            call prc(sp//'A VALUE <= 0.0 = DEFAULT')
         end if
      elseif (switch(swnmom).eq.2.0) then
         CALL PRC(S1//'MOMENTUM OPT  2 : MOMENTUM LOSS TERM IS ON')
         CALL PRC(SP//'SMOM = SMOM0 * EXP ( -S / (LM * SMAX))  FOR'//' S < L * SMAX')
         CALL PRC(SP//'     = 0        S > L * SMAX        ')
         CALL PRC(SP//'SMOM0= PT/(LM * SMAX)*(1/FFRIC -1) / EXPF')
         CALL PRC(SP//'EXPF = (1 - EXP(-(L * SMAX)/(LM * SMAX)))')
         CALL PRQ(SP//'LM   = ',LAMMOM)
         CALL PRQ(SP//'L    = ',LENMOM)
         if (n_extffric.eq.0) then 
            CALL PRQ(SP//'FFRIC= ',FFRIC)
         else
            CALL PRQ(SP//'DEFAULT FFRIC= ',FFRIC)
            call prc(sp//'OVER-RIDDEN BY THE'//' FOLLOWING VALUES ON SPECIFIED RINGS:')
            call prc(sp//'         FFRIC  LAMMON   FFRIC   LAMMOM')
            call prc(sp//'RING    '//OUTER//'              '//INNER)
            do in = 1, n_extffric
               write(coment,'(i4,2x,8(1x,f8.3))')int(extffric(in,1)),extffric(in,2),extffric(in,3),extffric(in,4),&
                                 extffric(in,5)    !,extffric(in,6),extffric(in,7)
               !write(coment,'(i4,2x,2(1x,f8.3))')int(extffric(in,1)),extffric(in,2),extffric(in,3)
               call prc(sp//coment)
            end do 
         end if
      elseif (switch(swnmom).eq.3.0) then
         CALL PRC(S1//'MOMENTUM OPT  3 : MOMENTUM LOSS TERM IS ON')
         CALL PRC(SP//'SMOM = FACT * SIZ(S) / INT (0 TO L) '//'[ SIZ(S'') DS'' ] ')
         CALL PRC(SP//'FACT = PT *  (1/FFRIC -1) ')
         CALL PRC(SP//'L    = SMAX / 2.0')
         if (n_extffric.eq.0) then 
            CALL PRQ(SP//'FFRIC= ',FFRIC)
         else
            CALL PRQ(SP//'DEFAULT FFRIC= ',FFRIC)
            call prc(sp//'OVER-RIDDEN BY THE'//' FOLLOWING VALUES ON SPECIFIED RINGS:')
            call prc(sp//'RING    '//OUTER//'   '//INNER)
            do in = 1, n_extffric
               write(coment,'(i4,2x,2(1x,f8.3))')int(extffric(in,1)),extffric(in,2),extffric(in,3)
               call prc(sp//coment)
            end do 
         end if
      elseif (switch(swnmom).eq.4.0) then
         CALL PRC(S1//'MOMENTUM OPT  4 : MOMENTUM LOSS TERM IS ON')
         CALL PRC(SP//'SMOM = -MB * VB(S) * RCXMULT * RCXIZ'//' *  SIZ(S) ')
         CALL PRQ(SP//'RCXIZ= ',RCXMOM)
         CALL PRC(SP//'RCXMULT = A * EXP(-B * T) ')
         CALL PRC(SP//'B = 6.907/(TCX - 1) ')
         CALL PRC(SP//'A = 1000 * EXP (B)  ')
         CALL PRQ(SP//'TCX = ',TCXMOM)
         CALL PRC(SP//'WHERE: 1.0 < RCXMULT < 1500.0 IS FORCED')
      elseif (switch(swnmom).eq.5.0) then
         CALL PRC(S1//'MOMENTUM OPT  5 : MOMENTUM LOSS TERM IS ON')
         CALL PRC(SP//'SMOM IS READ DIRECTLY FROM NIMBUS')
      elseif (switch(swnmom).eq.6.0) then
         CALL PRC(S1//'MOMENTUM OPT  6 : MOMENTUM LOSS TERM IS ON')
         CALL PRC(SP//'NEUTRAL BG MOMENTUN LOSS TERM IS TAKEN')
         CALL PRC(SP//'FROM PIN RESULTS EXCEPT FOR THE FIRST')
         CALL PRC(SP//'ITERATION, IN WHICH CASE: ')
         CALL PRC(SP//'MOMENTUM LOSS TERM IS OFF')
      elseif (switch(swnmom).eq.7.0) then
         CALL PRC(S1//'MOMENTUM OPT  7 : MOMENTUM LOSS TERM IS ON')
         CALL PRC(SP//'NEUTRAL BG MOMENTUN LOSS TERM IS TAKEN')
         CALL PRC(SP//'FROM PIN RESULTS EXCEPT FOR THE FIRST')
         CALL PRC(SP//'ITERATION, IN WHICH CASE: ')
         CALL PRC(SP//'INITIAL NEUTRAL BG MOMENTUM LOSS TERM IS:')
         CALL PRC(SP//'SMOM = SMOM0    S < L * SMAX        ')
         CALL PRC(SP//'     = 0        S > L * SMAX        ')
         CALL PRC(SP//'SMOM0= PT/(L * SMAX) * (1/FFRIC -1)')
         CALL PRQ(SP//'FFRIC= ',FFRIC)
         CALL PRQ(SP//'L    = ',LENMOM)
      elseif (switch(swnmom).eq.8.0) then
         CALL PRC(S1//'MOMENTUM OPT  8 : MOMENTUM LOSS TERM IS ON')
         CALL PRC(SP//'NEUTRAL BG MOMENTUN LOSS TERM IS TAKEN')
         CALL PRC(SP//'FROM PIN RESULTS EXCEPT FOR THE FIRST')
         CALL PRC(SP//'ITERATION, IN WHICH CASE: ')
         CALL PRC(SP//'INITIAL NEUTRAL BG MOMENTUM LOSS TERM IS:')
         CALL PRC(SP//'SMOM = SMOM0 * EXP ( -S / (LM * SMAX))  FOR'//' S < L * SMAX')
         CALL PRC(SP//'     = 0        S > L * SMAX        ')
         CALL PRC(SP//'SMOM0= PT/(LM * SMAX)*(1/FFRIC -1) / EXPF')
         CALL PRC(SP//'EXPF = (1 - EXP(-(L * SMAX)/(LM * SMAX)))')
         CALL PRQ(SP//'FFRIC= ',FFRIC)
         CALL PRQ(SP//'LM   = ',LAMMOM)
         CALL PRQ(SP//'L    = ',LENMOM)
      elseif (switch(swnmom).eq.9.0) then
         CALL PRC(S1//'MOMENTUM OPT  9 : MOMENTUM LOSS TERM IS ON')
         CALL PRC(SP//'BASED ON CX EXCHANGE MOMENTUM CROSS-'//'SECTIONS')
         CALL PRC(SP//'TAKEN FROM THE EDGE2D/NIMBUS'//' IMPLEMENTATION.')
         CALL PRC(SP//'REF: WOJCIECH')
         CALL PRC(SP//'MOMENTUM LOSS MULTIPLIED BY RCXIZ')
         CALL PRQ(SP//'RCXIZ= ',RCXMOM)
      elseif (switch(swnmom).eq.10.0) then
         CALL PRC(S1//'MOMENTUM OPT 10 : MOMENTUM LOSS TERM IS ON')
         CALL PRC(SP//'BASED ON CX EXCHANGE MOMENTUM CROSS-'//'SECTIONS TAKEN FROM')
         CALL PRC(SP//'EDGE2D/NIMBUS IMPLEMENTATION.')
         CALL PRC(SP//'REF: WOJCIECH')
         CALL PRC(SP//'     SOME SORT OF ALTERNATIVE TO'//'OPTION 9')
         CALL PRC(SP//'     FACTORING IN H0 VELOCITY'//' DISTRIBUTIONS')
         CALL PRC(SP//'MOMENTUM LOSS MULTIPLIED BY RCXIZ')
         CALL PRQ(SP//'RCXIZ= ',RCXMOM)

!     Mulitplication factor for Momentum loss option

      endif
      if (switch(swnmom).ne.0.0) then
         CALL PRQ(SP//'MOMENTUM SOURCE MULT= ',smom_mult)

!     Iterative Mach number

      endif

      call prb
      if (switch(swmach).eq.0.0) then
         CALL PRC (S1//'MACH OPTION   0 : ITERATIVE MACH'//' SOLUTION IS OFF')
         CALL PRQ (SP//'IMPOSED MACH NUMBER AT THE TARGET ', M0)
      elseif (switch(swmach).eq.1.0) then
         CALL PRC (S1//'MACH OPTION   1 : ITERATIVE MACH'//' SOLUTION IS ON')
         CALL PRC (SP//'TARGET DENSITY IS FIXED WITH MACH NUMBER')
         CALL PRC (SP//'AS A RESULT, THE TARGET FLUX WILL CHANGE.')
      elseif (switch(swmach).eq.2.0) then
         CALL PRC (S1//'MACH OPTION   2 : ITERATIVE MACH'//' SOLUTION IS ON')
         CALL PRC (SP//'TARGET DENSITY VARIES WITH MACH NUMBER TO'//' MAINTAIN A CONSTANT TARGET FLUX.')
      elseif (switch(swmach).eq.3.0) then
         CALL PRC (S1//'MACH OPTION   3 : FIXED MACH NUMBER'//' SOLUTION IS ON')
         CALL PRC (SP//'MACH NUMBER IS SET BASED ON INPUT TARGET'//' FLOW VELOCITY RELATIVE TO SOUND SPEED.')

!     Edge2D compatibility

      endif

      call prb
      if (e2dstart.eq.1) then
         CALL PRC (S1//'WARNING: ')
         CALL PRB
         CALL PRC (S1//'EDGE2D STARTER PLASMA OPTION IS ON')
         CALL PRC (S1//'THE INITIAL PLASMA SOLUTION FOR SOL OPTION 22')
         CALL PRC (S1//'HAS BEEN READ FROM AN EDGE2D INPUT FILE.')
         CALL PRC (S1//'NONE OF THE SPECIFIED SWITCHES APPLY TO THE')
         CALL PRC (S1//'STARTER PLASMA.')
         CALL PRC (S1//'NOTE: THIS IS SET BY USING -OPT FOR THE EDGE2D')
         CALL PRC (S1//'         COMPATIBILITY OPTION.')
         CALL PRB

      endif
      if (switch(swe2d).eq.0.0) then
         CALL PRC(S1//'EDGE2D BG OPT 0 : COMPATIBILITY IS OFF')
         CALL PRC(SP//'SOLVER RUNS FROM THE ACTUAL TARGET')
         call prc(sp//'IF E2D TARGET OPTION 5 IS SELECTED THEN THE')
         CALL PRC(SP//'EDGE2D DOWN POWER FLUXES WILL ALSO BE USED IN')
         CALL PRC(SP//'THE SOLVER')
      elseif (switch(swe2d).eq.1.0) then
         CALL PRC(S1//'EDGE2D BG OPT 1 : COMPATIBILITY IS ON')
         CALL PRC(SP//'SOLVER RUNS FROM THE MIDDLE OF THE FIRST CELL')
         CALL PRC(SP//'EDGE2D DATA IS READ FOR FIRST CELL')
         CALL PRC(SP//'EDGE2D PRESSURE IS MATCHED AT FIRST')
         CALL PRC(SP//'CELL CENTRE')
      elseif (switch(swe2d).eq.2.0) then
         CALL PRC(S1//'EDGE2D BG OPT 2 : COMPATIBILITY IS ON')
         CALL PRC(SP//'SOLVER RUNS FROM THE MIDDLE OF THE FIRST CELL')
         CALL PRC(SP//'EDGE2D DATA IS READ FOR FIRST CELL')
         CALL PRC(SP//'VELOCITY IS READ FROM THE EDGE2D OUTPUT')
         CALL PRC(SP//'AT THE FIRST CELL CENTRE')
      elseif (switch(swe2d).eq.3.0) then
         CALL PRC(S1//'EDGE2D BG OPT 3 : COMPATIBILITY IS ON')
         CALL PRC(SP//'SOLVER RUNS FROM THE MIDDLE OF THE FIRST CELL')
         CALL PRC(SP//'EDGE2D DATA IS READ FOR FIRST CELL')
         CALL PRC(SP//'VELOCITY IS SET BY EXTRACTING THE EDGE2D')
         CALL PRC(SP//'FLUXES ENTERING AND LEAVING THE FIRST CELL')
         CALL PRC(SP//'AND TAKING THEIR AVERAGE AND DIVIDING BY')
         CALL PRC(SP//'EDGE2D CELL DENSITY.')
         call prc(sp//'** NOTE: SOL 22 IS NOT RUN. THE EDGE2D')
         CALL PRC(SP//'         SOLUTION IS USED FOR THE BACKGROUND')
         CALL PRC(SP//'         PLASMA. THIS OPTION IS USED ONLY TO')
         CALL PRC(SP//'         PRODUCE DETAILED FLUX ANALYSES FOR')
         CALL PRC(SP//'         DEBUGGING.')
      elseif (switch(swe2d).eq.4.0) then
         CALL PRC(S1//'EDGE2D BG OPT 4 : COMPATIBILITY IS ON')
         CALL PRC(SP//'EDGE2D DATA IS READ FOR FIRST CELL')
         CALL PRC(SP//'VELOCITY IS READ FROM THE EDGE2D OUTPUT')
         CALL PRC(SP//'AT THE FIRST CELL CENTRE')
         call prc(sp//'** NOTE: SOL 22 IS NOT RUN. THE EDGE2D')
         CALL PRC(SP//'         SOLUTION IS USED FOR THE BACKGROUND')
         CALL PRC(SP//'         PLASMA. THIS OPTION IS USED ONLY TO')
         CALL PRC(SP//'         PRODUCE DETAILED FLUX ANALYSES FOR')
         CALL PRC(SP//'         DEBUGGING.')
      elseif (switch(swe2d).eq.5.0) then
         CALL PRC(S1//'EDGE2D BG OPT 5 : COMPATIBILITY IS OFF')
         call prc(sp//'FLUXES ARE CALCULATED FROM THE TARGET')
         call prc(sp//'** NOTE: SOL 22 IS NOT RUN. THE EDGE2D')
         CALL PRC(SP//'         SOLUTION IS USED FOR THE BACKGROUND')
         CALL PRC(SP//'         PLASMA. THIS OPTION IS USED ONLY TO')
         CALL PRC(SP//'         PRODUCE DETAILED FLUX ANALYSES FOR')
         CALL PRC(SP//'         DEBUGGING.')
      elseif (switch(swe2d).eq.6.0) then
         CALL PRC(S1//'EDGE2D BG OPT 6 : COMPATIBILITY IS ON')
         CALL PRC(SP//'EDGE2D DATA IS READ FOR FIRST CELL')
         CALL PRC(SP//'EDGE2D PRESSURE IS MATCHED AT FIRST')
         CALL PRC(SP//'CELL CENTRE')
         call prc(sp//'** NOTE: SOL 22 IS NOT RUN. THE EDGE2D')
         CALL PRC(SP//'         SOLUTION IS USED FOR THE BACKGROUND')
         CALL PRC(SP//'         PLASMA. THIS OPTION IS USED ONLY TO')
         CALL PRC(SP//'         PRODUCE DETAILED FLUX ANALYSES FOR')
         CALL PRC(SP//'         DEBUGGING.')
      elseif (switch(swe2d).eq.7.0) then
         CALL PRC(S1//'EDGE2D BG OPT 7 : COMPATIBILITY IS ON')
         CALL PRC(SP//'EDGE2D DATA IS READ FOR FIRST CELL')
         CALL PRC(SP//'VELOCITY IS SET BY EXTRACTING THE EDGE2D')
         CALL PRC(SP//'FLUXES ENTERING AND LEAVING THE FIRST CELL')
         CALL PRC(SP//'AND TAKING THEIR AVERAGE AND DIVIDING BY')
         CALL PRC(SP//'EDGE2D CELL DENSITY.')
         call prc(sp//'** NOTE: SOL 22 IS NOT RUN. THE EDGE2D')
         CALL PRC(SP//'         SOLUTION IS USED FOR THE BACKGROUND')
         CALL PRC(SP//'         PLASMA. THIS OPTION IS USED ONLY TO')
         CALL PRC(SP//'         PRODUCE DETAILED FLUX ANALYSES FOR')
         CALL PRC(SP//'         DEBUGGING.')
      elseif (switch(swe2d).eq.8.0) then
         CALL PRC(S1//'EDGE2D BG OPT 8 : COMPATIBILITY IS ON')
         CALL PRC(SP//'SOLVER RUNS FROM THE MIDDLE OF THE FIRST CELL')
         CALL PRC(SP//'EDGE2D DATA IS REQUIRED FOR THE ENTIRE RING.')
         CALL PRC(SP//'VELOCITY IS SET BY EXTRACTING THE EDGE2D')
         CALL PRC(SP//'FLUXES ENTERING AND LEAVING THE FIRST CELL')
         CALL PRC(SP//'AND TAKING THEIR AVERAGE AND DIVIDING BY')
         CALL PRC(SP//'EDGE2D CELL DENSITY.')
         call prc(sp//'THESE FLUXES ARE EXTRACTED FROM THE EDGE2D')
         call prc(sp//'DOWN FLUX LISTING.')
         call prc(sp//'IF E2D TARGET OPTION 5 IS SELECTED THEN THE')
         CALL PRC(SP//'EDGE2D DOWN POWER FLUXES AT THE CELL CENTRE')
         CALL PRC(SP//'WILL ALSO BE USED IN THE SOLVER')
      elseif (switch(swe2d).eq.9.0) then
         CALL PRC(S1//'EDGE2D BG OPT 9 : COMPATIBILITY IS ON')
         CALL PRC(SP//'SOLVER RUNS FROM THE MIDDLE OF A SPECIFIED CELL')
         CALL PRI(SP//'STARTING KNOT INDEX = ',IKE2D)
         CALL PRC(SP//'EDGE2D DATA IS REQUIRED FOR THE ENTIRE RING.')
         CALL PRC(SP//'VELOCITY IS SET BY EXTRACTING THE EDGE2D')
         CALL PRC(SP//'FLUXES ENTERING AND LEAVING THE CELL')
         CALL PRC(SP//'AND TAKING THEIR AVERAGE AND DIVIDING BY')
         CALL PRC(SP//'EDGE2D CELL DENSITY.')
         call prc(sp//'THESE FLUXES ARE EXTRACTED FROM THE EDGE2D')
         call prc(sp//'DOWN FLUX LISTING.')
         call prc(sp//'IF E2D TARGET OPTION 5 IS SELECTED THEN THE')
         CALL PRC(SP//'EDGE2D DOWN POWER FLUXES AT THE CELL CENTRE')
         CALL PRC(SP//'WILL ALSO BE USED IN THE SOLVER')

!     Power Distribution Option


      endif

      call prb
      if (switch(swpow).eq.0.0) then
         CALL PRC(S1//'POWER DIST OPT 0: DISTRIBUTED INFLUX IS OFF')
      elseif (switch(swpow).eq.1.0) then
         CALL PRC(S1//'POWER DIST OPT 1: DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'ALONG THE HALF RING TO THE MIDPOINT')
      elseif (switch(swpow).eq.2.0) then
         CALL PRC(S1//'POWER DIST OPT 2: DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'FROM THE X-POINT REGION TO THE MIDPOINT')
      elseif (switch(swpow).eq.3.0) then
         CALL PRC(S1//'POWER DIST OPT 3: DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'ALONG THE HALF RING TO THE MIDPOINT.')
         CALL PRC(SP//'CORRECTED FOR MAJOR RADIUS EFFECTS. THIS')
         CALL PRC(SP//'OPTION USEFUL ONLY WITH MAJOR RADIUS'//' OPTION 4')
      elseif (switch(swpow).eq.4.0) then
         CALL PRC(S1//'POWER DIST OPT 4: DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS SUMMED FOR BOTH')
         CALL PRC(SP//'TARGETS AND DISTRIBUTED LINEARLY')
         CALL PRC(SP//'OVER THE ENTIRE RING.')
      elseif (switch(swpow).eq.5.0) then
         CALL PRC(S1//'POWER DIST OPT 5: DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'ALONG THE HALF RING TO THE MIDPOINT')
         CALL PRC(SP//'IF PIN BASED VOLUME POWER SOURCES ARE ON')
         CALL PRC(SP//'THEN THE INTEGRAL OF THESE IS ALSO'//' DISTRIBUTED.')
      elseif (switch(swpow).eq.6.0) then
         CALL PRC(S1//'POWER DIST OPT 6: DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS SUMMED FOR BOTH')
         CALL PRC(SP//'TARGETS AND DISTRIBUTED LINEARLY')
         CALL PRC(SP//'OVER THE ENTIRE RING.')
         CALL PRC(SP//'IF PIN BASED VOLUME POWER SOURCES ARE ON')
         CALL PRC(SP//'THEN THE INTEGRAL OF THESE IS ALSO'//' DISTRIBUTED.')
      elseif (switch(swpow).eq.7.0) then
         CALL PRC(S1//'POWER DIST OPT 7: DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'FROM A SPECIFIED POINT TO THE MIDPOINT')
         CALL PRQ(SP//'SPECIFIED POINT = SMAX * ',SPOWBEG)
      elseif (switch(swpow).eq.8.0) then
         CALL PRC(S1//'POWER DIST OPT 8: DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'FROM A SPECIFIED POINT TO THE MIDPOINT')
         CALL PRC(SP//'IF PIN BASED VOLUME POWER SOURCES ARE ON')
         CALL PRC(SP//'THEN THE INTEGRAL OF THESE IS ALSO'//' DISTRIBUTED.')
         CALL PRQ(SP//'SPECIFIED POINT = SMAX * ',SPOWBEG)
      elseif (switch(swpow).eq.9.0) then
         CALL PRC(S1//'POWER DIST OPT 9: DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'OVER A SPECIFIED RANGE.')
         CALL PRQ(SP//'STARTING POINT = SMAX * ',SPOWBEG)
         CALL PRQ(SP//'ENDING   POINT = SMAX * ',SPOWLEN)
      elseif (switch(swpow).eq.10.0) then
         CALL PRC(S1//'POWER DIST OP 10: DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'OVER A SPECIFIED RANGE.')
         CALL PRC(SP//'IF PIN BASED VOLUME POWER SOURCES ARE ON')
         CALL PRC(SP//'THEN THE INTEGRAL OF THESE IS ALSO'//' DISTRIBUTED.')
         CALL PRQ(SP//'STARTING POINT = SMAX * ',SPOWBEG)
         CALL PRQ(SP//'ENDING   POINT = SMAX * ',SPOWLEN)
      elseif (switch(swpow).eq.11.0) then
         CALL PRC(S1//'POWER DIST OP 11: DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS SUMMED FOR BOTH')
         CALL PRC(SP//'TARGETS AND DISTRIBUTED LINEARLY')
         CALL PRC(SP//'OVER THE RING FROM :')
         CALL PRQ(SP//'STARTING POINT = SMAX * ',SPOWBEG)
         CALL PRC(SP//'(FROM EACH TARGET) ')

      endif

!     Private Plasma Power Distribution

      call prb
      if (switch(swpowp).eq.0.0) then
         CALL PRC(S1//'PP POWER OPT  0 : DISTRIBUTED INFLUX IS OFF')
      elseif (switch(swpowp).eq.1.0) then
         CALL PRC(S1//'PP POWER OPT  1 : DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'ALONG THE HALF RING TO THE MIDPOINT')
      elseif (switch(swpowp).eq.2.0) then
         CALL PRC(S1//'PP POWER OPT  2 : DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'FROM THE X-POINT REGION TO THE MIDPOINT')
      elseif (switch(swpowp).eq.3.0) then
         CALL PRC(S1//'PP POWER OPT  3 : DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'ALONG THE HALF RING TO THE MIDPOINT.')
         CALL PRC(SP//'CORRECTED FOR MAJOR RADIUS EFFECTS. THIS')
         CALL PRC(SP//'OPTION USEFUL ONLY WITH MAJOR RADIUS'//' OPTION 4')
      elseif (switch(swpowp).eq.4.0) then
         CALL PRC(S1//'PP POWER OPT  4 : DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS SUMMED FOR BOTH')
         CALL PRC(SP//'TARGETS AND DISTRIBUTED LINEARLY')
         CALL PRC(SP//'OVER THE ENTIRE RING.')
      elseif (switch(swpowp).eq.5.0) then
         CALL PRC(S1//'PP POWER OPT  5 : DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'ALONG THE HALF RING TO THE MIDPOINT')
         CALL PRC(SP//'IF PIN BASED VOLUME POWER SOURCES ARE ON')
         CALL PRC(SP//'THEN THE INTEGRAL OF THESE IS ALSO'//' DISTRIBUTED.')
      elseif (switch(swpowp).eq.6.0) then
         CALL PRC(S1//'PP POWER OPT  6 : DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS SUMMED FOR BOTH')
         CALL PRC(SP//'TARGETS AND DISTRIBUTED LINEARLY')
         CALL PRC(SP//'OVER THE ENTIRE RING.')
         CALL PRC(SP//'IF PIN BASED VOLUME POWER SOURCES ARE ON')
         CALL PRC(SP//'THEN THE INTEGRAL OF THESE IS ALSO'//' DISTRIBUTED.')
      elseif (switch(swpowp).eq.7.0) then
         CALL PRC(S1//'PP POWER OPT  7 : DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'FROM A SPECIFIED POINT TO THE MIDPOINT.')
         CALL PRQ(SP//'SPECIFIED POINT = SMAX * ',SPOWBEG)
      elseif (switch(swpowp).eq.8.0) then
         CALL PRC(S1//'PP POWER OPT  8 : DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'FROM A SPECIFIED POINT TO THE MIDPOINT')
         CALL PRC(SP//'IF PIN BASED VOLUME POWER SOURCES ARE ON')
         CALL PRC(SP//'THEN THE INTEGRAL OF THESE IS ALSO'//' DISTRIBUTED.')
         CALL PRQ(SP//'SPECIFIED POINT = SMAX * ',SPOWBEG)
      elseif (switch(swpowp).eq.9.0) then
         CALL PRC(S1//'PP POWER OPT  9 : DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'OVER A SPECIFIED RANGE.')
         CALL PRQ(SP//'STARTING POINT = SMAX * ',SPOWBEG)
         CALL PRQ(SP//'ENDING   POINT = SMAX * ',SPOWLEN)
      elseif (switch(swpowp).eq.10.0) then
         CALL PRC(S1//'PP POWER OPT  10: DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'OVER A SPECIFIED RANGE.')
         CALL PRC(SP//'IF PIN BASED VOLUME POWER SOURCES ARE ON')
         CALL PRC(SP//'THEN THE INTEGRAL OF THESE IS ALSO'//' DISTRIBUTED.')
         CALL PRQ(SP//'STARTING POINT = SMAX * ',SPOWBEG)
         CALL PRQ(SP//'ENDING   POINT = SMAX * ',SPOWLEN)
      elseif (switch(swpow).eq.11.0) then
         CALL PRC(S1//'PP POWER OPT  11: DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS SUMMED FOR BOTH')
         CALL PRC(SP//'TARGETS AND DISTRIBUTED LINEARLY')
         CALL PRC(SP//'OVER THE RING FROM :')
         CALL PRQ(SP//'STARTING POINT = SMAX * ',SPOWBEG)
         CALL PRC(SP//'(FROM EACH TARGET) ')

      endif

!     Perpendicular flux option

      call prb
      if (switch(swgperp).eq.0.0) then
         CALL PRC(S1//'PERP FLUX OPT  0: TERM IS OFF')
      elseif (switch(swgperp).eq.1.0) then
         CALL PRC(S1//'PERP FLUX OPT  1: TERM IS ON')
         CALL PRC(SP//'NET FLUX ALONG FIELD LINE ->0 AT THE'//' MIDPOINT')
         CALL PRC(SP//'USING AN EVENLY DISTRIBUTED'//' PERPENDICULAR FLUX.')
      elseif (switch(swgperp).eq.2.0) then
         CALL PRC(S1//'PERP FLUX OPT  2: TERM IS ON')
         CALL PRC(SP//'NET FLUX OVER ENTIRE FIELD LINE ->0')
         CALL PRC(SP//'USING AN EVENLY DISTRIBUTED'//' PERPENDICULAR FLUX.')
      elseif (switch(swgperp).eq.3.0) then
         CALL PRC(S1//'PERP FLUX OPT  3: TERM IS ON')
         CALL PRC(SP//'NET FLUX OVER ENTIRE FIELD LINE ->0')
         CALL PRC(SP//'USING AN EVENLY DISTRIBUTED'//' PERPENDICULAR FLUX.')
         CALL PRC(SP//'FOR UNDER-IONIZED CASES AND A SOURCE')
         CALL PRC(SP//'DISTRIBUTED PROPORTIONAL TO NE FROM THE')
         CALL PRC(SP//'PREVIOUS ITERATION FOR AN OVER-IONIZED')
         CALL PRC(SP//'CASE.')
      elseif (switch(swgperp).eq.4.0) then
         CALL PRC(S1//'PERP FLUX OPT  4: TERM IS ON')
         CALL PRC(SP//'NET FLUX OVER THE HALF FIELD LINE ->0')
         CALL PRC(SP//'AT THE MIDPOINT -  ')
         CALL PRC(SP//'USING AN EVENLY DISTRIBUTED'//' PERPENDICULAR FLUX.')
         CALL PRC(SP//'FOR UNDER-IONIZED CASES AND A SOURCE')
         CALL PRC(SP//'DISTRIBUTED PROPORTIONAL TO NE FROM THE')
         CALL PRC(SP//'PREVIOUS ITERATION FOR AN OVER-IONIZED')
         CALL PRC(SP//'CASE.')
      elseif (switch(swgperp).eq.5.0) then
         CALL PRC(S1//'PERP FLUX OPT  5: TERM IS ON')
         CALL PRC(SP//'NET FLUX ALONG FIELD LINE ->0 AT THE'//' MIDPOINT')
         CALL PRC(SP//'USING AN EVENLY DISTRIBUTED'//' PERPENDICULAR FLUX OVER THE 1/2 FIELD LINE')
         CALL PRC(SP//'AND A RECTANGULAR SOURCE BETWEEN THE'//' FOLLOWING POINTS FROM EACH TARGET')
         CALL PRQ(SP//'START = SMAX * ',GPERPBEGF)
         CALL PRQ(SP//'END   = SMAX * ',GPERPENDF)
         CALL PRQ(SP//'FRACTION IN RECTANGULAR = ',GPERPFRAC)
         CALL PRQ(SP//'FRACTION UNIFORM        = ',1.0D0-GPERPFRAC)
      elseif (switch(swgperp).eq.6.0) then
         CALL PRC(S1//'PERP FLUX OPT  6: TERM IS ON')
         CALL PRC(SP//'NET FLUX OVER THE ENTIRE FILED LINE ->0')
         CALL PRC(SP//'USING AN EVENLY DISTRIBUTED'//' PERPENDICULAR FLUX OVER THE WHOLE FIELD LINE')
         CALL PRC(SP//'AND A RECTANGULAR SOURCE BETWEEN THE'//' FOLLOWING POINTS FROM EACH TARGET')
         CALL PRQ(SP//'START = SMAX * ',GPERPBEGF)
         CALL PRQ(SP//'END   = SMAX * ',GPERPENDF)
         CALL PRQ(SP//'FRACTION IN RECTANGULAR = ',GPERPFRAC)
         CALL PRQ(SP//'FRACTION UNIFORM        = ',1.0D0-GPERPFRAC)
      elseif (switch(swgperp).eq.7.0) then
         CALL PRC(S1//'PERP FLUX OPT  7: TERM IS ON')
         CALL PRC(SP//'NET FLUX OVER ENTIRE FIELD LINE ->0')
         CALL PRC(SP//'USING A PERPENDICULAR FLUX PROPORTIONAL')
         CALL PRC(SP//'TO THE SECOND GRADIENT IN THE DENSITY')
         CALL PRC(SP//'FROM AN EDGE2D SOLUTION OR A PREVIOUS')
         CALL PRC(SP//'SOL22 SOLUTION. THIS OPTION IS CHANGED')
         CALL PRC(SP//'TO A UNIFORM DISTRIBUTION FOR ANY RINGS')
         CALL PRC(SP//'WHERE THE RATIO OF EITHER THE POSITIVE')
         CALL PRC(SP//'OR NEGATIVE CONTRIBUTION TO THE INTEGRATED')
         CALL PRC(SP//'VALUE FOR TEH RING EXCEEDS 5.0' )
      elseif (switch(swgperp).eq.8.0) then
         CALL PRC(S1//'PERP FLUX OPT  8: TERM IS ON')
         CALL PRC(SP//'NET FLUX OVER ENTIRE FIELD LINE ->0')
         CALL PRC(SP//'USING A PERPENDICULAR FLUX CALCULATED')
         CALL PRC(SP//'FROM THE SECOND GRADIENT IN THE DENSITY')
         call PRC(SP//'AND AN IMPOSED VALUE FOR THE DIFFUSION')
         CALL PRR(SP//'RATE:  DPERP [M2/S]  = ',CDPERP)
         CALL PRC(SP//'ANY REMAINING EXECESS OR SHORTFALL OF')
         CALL PRC(SP//'FLUX IS COMPENSATED FOR BY THE APPLICATION')
         CALL PRC(SP//'OF AN ADDITIONAL UNIFORM SOURCE OVER THE')
         CALL PRC(SP//'ENTIRE FLUX TUBE.')

      endif

      call prb
      if (switch(swgperpp).eq.0.0) then
         CALL PRC(S1//'PP PERPF  OPT  0: TERM IS OFF')
      elseif (switch(swgperpp).eq.1.0) then
         CALL PRC(S1//'PP PERPF  OPT  1: TERM IS ON')
         CALL PRC(SP//'NET FLUX ALONG FIELD LINE ->0 AT THE'//' MIDPOINT')
         CALL PRC(SP//'USING AN EVENLY DISTRIBUTED'//' PERPENDICULAR FLUX.')
      elseif (switch(swgperpp).eq.2.0) then
         CALL PRC(S1//'PP PERPF  OPT  2: TERM IS ON')
         CALL PRC(SP//'NET FLUX OVER ENTIRE FIELD LINE ->0')
         CALL PRC(SP//'USING AN EVENLY DISTRIBUTED'//' PERPENDICULAR FLUX.')
      elseif (switch(swgperpp).eq.3.0) then
         CALL PRC(S1//'PP PERPF  OPT  3: TERM IS ON')
         CALL PRC(SP//'NET FLUX OVER ENTIRE FIELD LINE ->0')
         CALL PRC(SP//'USING AN EVENLY DISTRIBUTED'//' PERPENDICULAR FLUX.')
         CALL PRC(SP//'FOR UNDER-IONIZED CASES AND A SOURCE')
         CALL PRC(SP//'DISTRIBUTED PROPORTIONAL TO NE FROM THE')
         CALL PRC(SP//'PREVIOUS ITERATION FOR AN OVER-IONIZED')
         CALL PRC(SP//'CASE.')
      elseif (switch(swgperpp).eq.4.0) then
         CALL PRC(S1//'PP PERPF  OPT  4: TERM IS ON')
         CALL PRC(SP//'NET FLUX OVER THE HALF FIELD LINE ->0')
         CALL PRC(SP//'AT THE MIDPOINT -  ')
         CALL PRC(SP//'USING AN EVENLY DISTRIBUTED'//' PERPENDICULAR FLUX.')
         CALL PRC(SP//'FOR UNDER-IONIZED CASES AND A SOURCE')
         CALL PRC(SP//'DISTRIBUTED PROPORTIONAL TO NE FROM THE')
         CALL PRC(SP//'PREVIOUS ITERATION FOR AN OVER-IONIZED')
         CALL PRC(SP//'CASE.')
      elseif (switch(swgperpp).eq.5.0) then
         CALL PRC(S1//'PP PERPF  OPT  5: TERM IS ON')
         CALL PRC(SP//'NET FLUX ALONG FIELD LINE ->0 AT THE'//' MIDPOINT')
         CALL PRC(SP//'USING AN EVENLY DISTRIBUTED PERPENDICULAR')
         CALL PRC(SP//'FLUX OVER THE 1/2 FIELD LINE')
         CALL PRC(SP//'AND A RECTANGULAR SOURCE BETWEEN THE')
         call prc(sp//'FOLLOWING POINTS FROM EACH TARGET')
         CALL PRQ(SP//'START = SMAX * ',GPERPBEGF)
         CALL PRQ(SP//'END   = SMAX * ',GPERPENDF)
         CALL PRQ(SP//'FRACTION IN RECTANGULAR = ',GPERPFRAC)
         CALL PRQ(SP//'FRACTION UNIFORM        = ',1.0D0-GPERPFRAC)
      elseif (switch(swgperpp).eq.6.0) then
         CALL PRC(S1//'PP PERPF  OPT  6: TERM IS ON')
         CALL PRC(SP//'NET FLUX OVER THE ENTIRE FILED LINE ->0')
         CALL PRC(SP//'USING AN EVENLY DISTRIBUTED PERPENDICULAR')
         CALL PRC(SP//'FLUX OVER THE WHOLE FIELD LINE')
         CALL PRC(SP//'AND A RECTANGULAR SOURCE BETWEEN THE')
         CALL PRC(SP//'FOLLOWING POINTS FROM EACH TARGET')
         CALL PRQ(SP//'START = SMAX * ',GPERPBEGF)
         CALL PRQ(SP//'END   = SMAX * ',GPERPENDF)
         CALL PRQ(SP//'FRACTION IN RECTANGULAR = ',GPERPFRAC)
         CALL PRQ(SP//'FRACTION UNIFORM        = ',1.0D0-GPERPFRAC)
      elseif (switch(swgperpP).eq.7.0) then
         CALL PRC(S1//'PERP FLUX OPT  7: TERM IS ON')
         CALL PRC(SP//'NET FLUX OVER ENTIRE FIELD LINE ->0')
         CALL PRC(SP//'USING A PERPENDICULAR FLUX PROPORTIONAL')
         CALL PRC(SP//'TO THE SECOND GRADIENT IN THE DENSITY')
         CALL PRC(SP//'FROM AN EDGE2D SOLUTION OR A PREVIOUS')
         CALL PRC(SP//'SOL22 SOLUTION. THIS OPTION IS CHANGED')
         CALL PRC(SP//'TO A UNIFORM DISTRIBUTION FOR ANY RINGS')
         CALL PRC(SP//'WHERE THE RATIO OF EITHER THE POSITIVE')
         CALL PRC(SP//'OR NEGATIVE CONTRIBUTION TO THE INTEGRATED')
         CALL PRC(SP//'VALUE FOR TEH RING EXCEEDS 5.0' )
      elseif (switch(swgperpP).eq.8.0) then
         CALL PRC(S1//'PERP FLUX OPT  8: TERM IS ON')
         CALL PRC(SP//'NET FLUX OVER ENTIRE FIELD LINE ->0')
         CALL PRC(SP//'USING A PERPENDICULAR FLUX CALCULATED')
         CALL PRC(SP//'FROM THE SECOND GRADIENT IN THE DENSITY')
         call PRC(SP//'AND AN IMPOSED VALUE FOR THE DIFFUSION')
         CALL PRR(SP//'RATE:  DPERP [M2/S]  = ',CDPERP)
         CALL PRC(SP//'ANY REMAINING EXECESS OR SHORTFALL OF')
         CALL PRC(SP//'FLUX IS COMPENSATED FOR BY THE APPLICATION')
         CALL PRC(SP//'OF AN ADDITIONAL UNIFORM SOURCE OVER THE')
         CALL PRC(SP//'ENTIRE FLUX TUBE.')

!     Extra Gperp source/sink

      endif
      if (switch(swextra).eq.0.0) then
         CALL PRC(S1//'EXTRA PERP FLUX 0 : TERM IS OFF')
      elseif (switch(swextra).eq.1.0) then
         CALL PRC(S1//'EXTRA PERP FLUX 1 : TERM IS ON')
         CALL PRC(SP//'AN EXTRA SOURCE AND SINK ARE SUPERIMPOSED ON')
         CALL PRC(SP//'THE FLOW PATTERN FOR THE FLUX TUBE. THIS')
         CALL PRC(SP//'SOURCE AND SINK EXACTLY CANCEL BUT WILL')
         CALL PRC(SP//'AFFECT THE FLOW PATTERN ON THE FLUX TUBE')
         call prq(sp//'SOURCE SIZE = (TOTAL TARGET FLUX ON RING) * ',gextra_mult)
         call prc(sp//'THE SOURCE IS IMPOSED OVER THE REGION:')
         call prq2(sp//'   - SMAX  * [F1,F2]  = ',gextra_src_start,gextra_src_stop)
         call prc(sp//'THE SINK IS IMPOSED OVER THE REGION:')
         call prq2(sp//'   - SMAX  * [F1,F2]  = ',gextra_sink_start,gextra_sink_stop)

      endif

!     Major Radius Option


      call prb

      call prb
      if (switch(swmajr).eq.0.0) then
         CALL PRC(S1//'MAJOR RADIUS   0: MAJOR RADIUS FACTOR IS OFF')
      elseif (switch(swmajr).eq.1.0) then
         CALL PRC(S1//'MAJOR RADIUS   1: MAJOR RADIUS FACTOR IS ON')
         CALL PRC(SP//'ALL TARGET FLUX ADJUSTED BY RTARG/R0')
      elseif (switch(swmajr).eq.2.0) then
         CALL PRC(S1//'MAJOR RADIUS   2: MAJOR RADIUS FACTOR IS ON')
         CALL PRC(SP//'ALL IONIZATION ADJUSTED BY RCELL/R0')
      elseif (switch(swmajr).eq.3.0) then
         CALL PRC(S1//'MAJOR RADIUS   3: MAJOR RADIUS FACTOR IS ON')
         CALL PRC(SP//'ALL IONIZATION ADJUSTED BY R0/RCELL')
      elseif (switch(swmajr).eq.4.0) then
         CALL PRC(S1//'MAJOR RADIUS   4: MAJOR RADIUS FACTOR IS ON')
         CALL PRC(SP//'GENERALIZED R-CORRECTION APPLIED')

!     Core Flux option

      endif
      if (switch(swcore).ne.0.0) then
         call prb
         CALL PRC(S1//'CORE IONIZATION CORRECTION TERM'//' IS NOT IMPLEMENTED')

!      if (switch(swcore).eq.0.0) then
!         call prc('       Cross-field Core Ionization Source is OFF')
!      elseif (switch(swcore).eq.1.0) then
!         call prc('       Cross-field Core Ionization Source is ON')
!         call prq('       - Core source fraction : ',corefrac)
!         call prc('       - Core source flux is distribited linearly')
!         call prc('         along the entire ring to the midpoint')
!      elseif (switch(swcore).eq.2.0) then
!         call prc('       Cross-field Core Ionization Source is ON')
!         call prq('       - Core source fraction : ',corefrac)
!         call prc('       - Core source flux is distribited linearly')
!         call prc('         from the X-point region to the midpoint')
!      endif


!     Recombination term

      endif

      call prb
      if (switch(swrecom).eq.0.0) then
         CALL PRC(S1//'RECOMB OPT    0 : TERM IS OFF')
         CALL PRC(SP//'NO  RECOMBINATION PARTICLE SOURCE')
      elseif (switch(swrecom).eq.1.0) then
         CALL PRC(S1//'RECOMB OPT    1 : TERM IS ON')
         CALL PRC(SP//'RECOMBINATION PARTICLE SOURCE ADDED TO FLUX')
         CALL PRR(SP//'TE CUTOFF FOR RECOMBINATION = ',TRECCUT)
         CALL PRC(SP//'RECOMBINATION RATES ARE LIMITED BY THE'//' T LISTED.')
      elseif (switch(swrecom).eq.2.0) then
         CALL PRC(S1//'RECOMB OPT    2 : TERM IS ON')
         CALL PRC(SP//'RECOMBINATION PARTICLE SOURCE ADDED TO FLUX')
         CALL PRC(SP//'RECOMBINATION CALCULATED FROM EDGE2D BG'//' INPUT')

      endif

!     Smoothing Option

      call prb
      if (switch(swsmooth).eq.0.0) then
         CALL PRC(S1//'SMOOTHING OPT 0 : SMOOTHING IS'//' TURNED OFF')
         CALL PRC(sp//'VALUES WILL NOT MATCH AT THE'//' MID-POINT')
      elseif (switch(swsmooth).eq.1.0) then
         CALL PRC(S1//'SMOOTHING OPT 1 : SMOOTHING IS'//' TURNED ON')
         CALL PRC(SP//'BG QUANTITIES ARE ADJUSTED'//' TO MATCH AT THE MID-POINT')
      endif

!     Detached Solution Option

      call prb
      if (switch(swdetach).eq.0.0) then
         call prc(s1//'DETACHED  OPT 0 : DETACHED PLASMA'//' OVERRIDE OPTION IS OFF')
      elseif (switch(swdetach).eq.1.0) then
         call prc(s1//'DETACHED  OPT 1 : DETACHED PLASMA MODEL'//' OVERRIDES')

         call prc(sp//'SOL22 FOR '//outer//' HALF RING (IK=1)')
       if (s21refsw.eq.0) then
        call prc(sp//'NOTE: REF OPTION 0 - ')
        call prc(sp//'All distances are specified')
        call prc(sp//'in units of SMAX. (MF=SMAX)')
       elseif (s21refsw.eq.1) then
        call prc(sp//'NOTE: REF OPTION 1 - ')
        call prc(sp//'All distances are specified')
        call prc(sp//'in units of PMAX. (MF=PMAX)')
        call prc(sp//'Converted to S for each'//' ring.')
       elseif (s21refsw.eq.2) then
        call prc(sp//'NOTE: REF OPTION 2 - ')
        call prc(sp//'All distances are specified')
        call prc(sp//'in meters along S.(MF=1.0)')
       elseif (s21refsw.eq.3) then
        call prc(sp//'NOTE: REF OPTION 3 - ')
        call prc(sp//'All distances are specified')
        call prc(sp//'in meters along P.(MF=1.0)')
        call prc(sp//'Converted to S for each'//' ring.')

!      OUTER

       endif

       call prc ('                       '//outer//' HALF OF SOL:')
       call prc ('                       Three regions :   ')
       call prr ('                       A :  0 < MF * ',l1rat)
       write (coment,'(''                       B :'',g8.3,'' * MF < '',g8.3, '' * MF '')') l1rat,l2rat
       call prc(coment)
       call prr ('                       C : S  >  MF * ',l2rat)
       call prr ('                       Te at end of A: Te0 * ',terat)
       call prr ('                       Ti at end of A: Ti0 * ',tirat)
       call prr ('                       Ne at end of A: Ne0 * ',nrat)
       call prc ('                       Ne(s) = Ne0 + (Ne1-Ne0)'//' * (s/L)**ALPHA')
       call prr ('                       ALPHA = ',nalph)
       call prr ('                       Radiated power in B (Qrad): D0* ',qrat)
       call prc ('                       Te increases linearly in A')
       call prc ('                       Ti increases lineraly in A')
       call prc ('                       Ne increases linearly in A')
       call prc ('                       Velocity in A: v(s)=N0V0/n(s)')
       CALL PRC ('                       B: T=(T1**3.5+7/2K0*(D0(s-L1)+')
       call prc ('                         1/2(s-L1)**2*(Qrad/Lrad)))**(2/7)')
       CALL PRC ('                       C: T=(T2**3.5+7/2K0*(Qtot(s-L2)))**(2/7)')
       call prc ('                          Qtot = D0 + Qrad')
       call prc ('                       B,C: N(s) = N1 * (T(s)/T1)**(-1)')

       call prr ('                       Veocity in B,C linearly->0 at MF * ',lvrat)

      elseif (switch(swdetach).eq.2.0) then
         call prc(s1//'DETACHED  OPT 2 : DETACHED PLASMA MODEL'//' OVERRIDES')

         call prc(sp//'SOL22 FOR'//INNER//' HALF RING (IK=NKS(IR))')
       if (s21refsw.eq.0) then
        call prc(sp//'NOTE: REF OPTION 0 - ')
        call prc(sp//'All distances are specified')
        call prc(sp//'in units of SMAX. (MF=SMAX)')
       elseif (s21refsw.eq.1) then
        call prc(sp//'NOTE: REF OPTION 1 - ')
        call prc(sp//'All distances are specified')
        call prc(sp//'in units of PMAX. (MF=PMAX)')
        call prc(sp//'Converted to S for each'//' ring.')
       elseif (s21refsw.eq.2) then
        call prc(sp//'NOTE: REF OPTION 2 - ')
        call prc(sp//'All distances are specified')
        call prc(sp//'in meters along S.(MF=1.0)')
       elseif (s21refsw.eq.3) then
        call prc(sp//'NOTE: REF OPTION 3 - ')
        call prc(sp//'All distances are specified')
        call prc(sp//'in meters along P.(MF=1.0)')
        call prc(sp//'Converted to S for each'//' ring.')


!      INNER

       endif
       call prb

       call prc ('                       '//INNER//' HALF OF SOL:')
       call prc ('                       Three regions :   ')
       call prr ('                       A :  0 < MF * ',l1rati)
       write (coment,'(''                       B :'',g8.3,'' * MF < '',g8.3, '' * MF '')') l1rati,l2rati
       call prc(coment)
       call prr ('                       C : S  >  MF * ',l2rati)
       call prr ('                       Te at end of A: Te0 * ',terati)
       call prr ('                       Ti at end of A: Ti0 * ',tirati)
       call prr ('                       Ne at end of A: Ne0 * ',nrati)
       call prc ('                       Ne(s) = Ne0 + (Ne1-Ne0)'//' * (s/L)**ALPHA')
       call prr ('                       ALPHA = ',nalphi)
       call prr ('                       Radiated power in B (Qrad): Q0* ',qrati)
       call prc ('                       Te increases linearly in A')
       call prc ('                       Ti increases linearly in A')
       call prc ('                       Ne increases linearly in A')
       call prc ('                       Velocity in A: v(s)=N0V0/n(s)')
       CALL PRC ('                       B: T=(T1**3.5+7/2K0*(Q0(s-L1)+')
       call prc ('                         1/2(s-L1)**2*(Qrad/Lrad)))**(2/7)')
       CALL PRC ('                       C: T=(T2**3.5+7/2K0*(Qtot(s-L2)))**(2/7)')
       call prc ('                          Qtot = Q0 + Qrad')
       call prc ('                       B,C: N(s) = N1 * (T(s)/T1)**(-1)')

       call prr ('                       Veocity in B,C linearly->0 at MF * ',lvrati)

      endif

!     Error Correction

      call prb
      if (switch(swerror).eq.0.0) then
         CALL PRC(S1//'ERROR OPTION  0 : AUTOMATIC ERROR'//' CORRECTION IS TURNED OFF')
      elseif (switch(swerror).gt.0.0) then
         CALL PRC(S1//'ERROR OPTION  N : AUTOMATIC ERROR'//' CORRECTION IS TURNED ON')
         CALL PRR(SP//'ERROR CORRECTION STARTS AT LEVEL :',switch(swerror))

!     End of OPTIONS section

      endif
      call prb
      call prc(s1//'SOL 22: REPORT OF RESULTS')

!     Summary of over/under ionization

      call prb
      if (switch(swgperp).eq.2.or.switch(swgperp).eq.3.or.switch(swgperp).eq.7.or.switch(swgperp).eq.8.or.switch(swgperpp).eq.2.or.switch(swgperpp).eq.3.or.switch(swgperpp).eq.7.or.switch(swgperpp).eq.8) then
         call prc(s1//'SUMMARY OF OVER/UNDER IONIZATION BY RING:')
         call prb
         if (switch(swgperp).eq.7.0.or.switch(swgperpp).eq.7.0) then
            call prc('     IR        RCONST    '//'  Gperp RATIO    Gperp OPTION USED')

            do ir = irsep,nrs
               if (ir.le.irwall) then
                  swtmp = switch(swgperp)
               else
                  swtmp = switch(swgperpp)

               endif

               if (swtmp.eq.7.0.and.gperprat(ir).gt.5.0) swtmp = 2.0
               write(coment,'(3x,i5,4x,1p,g13.6,0p,2(4x,f7.2))')ir,rconst(ir),gperprat(ir),swtmp
               call prc(coment)

            end do

         else

            call prc('     IR        RCONST     Gperp OPTION USED')

            do ir = irsep,nrs
               if (ir.le.irwall) then
                  swtmp = switch(swgperp)
               else
                  swtmp = switch(swgperpp)

               endif
               write(coment,'(3x,i5,4x,1p,g13.6,0p,4x,f7.2)')ir,rconst(ir),swtmp
               call prc(coment)

            end do

         endif

         call prb

!     Rings which were fixed by default conditions.

      endif

      if (switch(swerror).eq.0.0) then
         cnt = 0
         do ir = irsep,irlim
            if (cerr(ir,1).ne.0.or.cerr(ir,2).ne.0) then
               cnt = cnt + 1
               if (cnt.eq.1) then
                  CALL PRC(S1//'SUMMARY OF ERRORS:')
                  CALL PRC(S1//'   - SOLVER HAD PROBLEMS WITH THESE RINGS:')
                  CALL PRC(S1//'RING        CODE   DESCRIPTION'//'       POSITION     ERROR OPTION')

               endif
               if (cerr(ir,2).ne.0) then
                  write(coment,1060) ir,outer,cerr(ir,2),errstr(cerr(ir,2)),cserr(ir,2),cdeferropt(ir,2)
                  call prc(coment)

               endif
               if (cerr(ir,1).ne.0) then
                  write(coment,1060) ir,inner,cerr(ir,1),errstr(cerr(ir,1)),cserr(ir,1),cdeferropt(ir,1)
                  call prc(coment)


               endif
            endif

         end do

         call prb
 1060    format(6x,i4,1x,a5,':',i4,3x,a,1x,g13.6,1x,f7.1)

      elseif (switch(swerror).ne.0.0) then
         CALL PRB
         CALL PRC (S1//'ERROR CORRECTION OPTION ACTIVATED:')

!        Summarize corrected rings.

         CALL PRB
         cnt = 0
         do ir = irsep,irlim
            if (cdeferr(ir,1).ne.0.or.cdeferr(ir,2).ne.0) then
               cnt = cnt + 1

!                 Print Error Correction description

               if (cnt.eq.1) then

                  call prerrdesc(switch(swerror))

                  call PRC(S1//'     ERROR SOLVER HAD A PROBLEM'//' WITH THESE RINGS:')

                  CALL PRC(S1//'RING        CODE   DESCRIPTION'//'       POSITION     ERROR OPTION')

               endif
               if (cdeferr(ir,2).ne.0) then
                  write(coment,1060) ir,outer,cdeferr(ir,2),errstr(cdeferr(ir,2)),cdefserr(ir,2),cdeferropt(ir,2)
                  call prc(coment)

               endif
               if (cdeferr(ir,1).ne.0) then
                  write(coment,1060) ir,inner,cdeferr(ir,1),errstr(cdeferr(ir,1)),cdefserr(ir,1),cdeferropt(ir,1)
                  call prc(coment)

               endif
            endif

         end do

         if (cnt.eq.0) then

            CALL PRC (S1//'NO RINGS REQUIRED ERROR'//' CORRECTION METHODS')

         endif

         call prb

!     Summary of Mach number results

      end if

      if (switch(swmach).eq.1.0.or.switch(swmach).eq.2.0) then
         call prc(s1//'MACH SOLVER ACTIVATED')

         CALL PRC(s1//'SUMMARY OF TARGET MACH NUMBERS:')
         do ir = irsep, irlim
            write(coment,1000) ir,inner,cmachno(ir,1),outer,cmachno(ir,2)
            call prc(coment)

         end do

         call prb
         CALL PRQ (S1//'DELTA M0 FOR FIRST MACH ITERATION  : ',DELTAM0)

!        Format for Mach number table

         CALL PRQ (S1//'ULTIMATE RESOLUTION OF MACH NUMBER : ',M0RES)

1000     format('       Ring : ',i4,'  ',a5,' M#: ',f12.6,'  ',a5,' M#: ',f12.6)

      elseif (switch(swmach).eq.3.0) then
         call prc(s1//'FIXED TARGET MACH NUMBERS:')

         CALL PRC(s1//'SUMMARY OF TARGET MACH NUMBERS:')
         do ir = irsep, irlim
            write(coment,1000) ir,inner,cmachno(ir,1),outer,cmachno(ir,2)
            call prc(coment)

         end do

      endif
      call prb
      call prc(s1//'SUMMARY OF INTEGRATED CONDUCTIVE POWER RATIOS:')

      call prb

      CALL PRC(S1//'RING        CONDe/Qetot  CONDi/Qitot'//'  COND/Qtot')

!        Target 2 - outer for JET geometries

      do ir = irsep,irlim
         write(coment,1061) ir,outer,sol22_power_ratio(ir,2,1),sol22_power_ratio(ir,2,2),sol22_power_ratio(ir,2,3)

!        Target 1 - inner for JET geometries

         call prc(coment)
         write(coment,1061) ir,inner,sol22_power_ratio(ir,1,1),sol22_power_ratio(ir,1,2),sol22_power_ratio(ir,1,3)

         call prc(coment)
      end do
 1061    format(6x,i4,' ',a5,':',3(1x,f11.7))



!     Output controls

!      call prb
!      if (graph.eq.0) then
!         call prc ('Plotting is turned OFF')
!      elseif (graph.eq.1) then
!         call prc ('Plotting is turned ON')
!      endif

!      if (graphaux.eq.0) then
!         call prc ('Auxilliary plots and tables are turned OFF')
!      elseif (graphaux.eq.1) then
!         call prc ('Auxilliary plots and tables are turned ON')
!      endif

!      if (graphvel.eq.0) then
!         call prc ('Velocity plots and tables are turned OFF')
!      elseif (graphvel.eq.1) then
!         call prc ('Velocity plots and tables are turned ON')
!      endif

!      call prb
!      if (graphran.eq.0.0) then
!         call prc ('Close up plots are turned OFF')
!      elseif (graphvel.eq.1) then
!         call prc ('Close up plots are turned ON')
!         call prr ('Range for close up plots is 0.0  to ',graphran)
!      endif


      call prb
      return

    end subroutine echosol



    end module mod_sol22_support
