c     -*-Fortran-*-
c
c ======================================================================
c
c
c ======================================================================
c
c subroutine: SOL22Status
c
c
c
      SUBROUTINE SOL22Status(region,ir,deltat,serr,spts,npts,errcode)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'
      INCLUDE 'solparams'
      INCLUDE 'solcommon'
      INCLUDE 'solswitch'

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
          IF (ABS(errcode).EQ.1) THEN
c
c           jdemod - changed halflen to halfringlen
c
            WRITE(0,'(A,I2,A,F5.1,A,I3,1X,A,1X,F6.2,A,F6.2,A)')
     .        'SOL22: Low  ',ir,' (',deltat,' s)  E: ',err1,
     .        errstr(INT(errcode)),
     .        simag1 / (halfringlen - soffset) * 100.0,' - ',
     .        simag2 / (halfringlen - soffset) * 100.0,'%'
          ELSE
            WRITE(0,'(A,I2,A,F5.1,A,I3,1X,A)')
     .        'SOL22: Low  ',ir,
     .        ' (',deltat,' s)  E: ',err1,errstr(INT(errcode))
          ENDIF
        ENDIF

        WRITE(cdum2(1:128),'(A,I2,1X,2F7.3,1P,E11.3,0P,A,
     .                       F5.1,A,I3)')
     .    'SOL22 L: ',ir,te0,ti0,n0,' (',deltat,' s) E:',errcode

        IF (ABS(errcode).EQ.1) THEN
          per = (serr - soffset) / (spts(npts) - soffset) * 100.0
c
c           jdemod - changed halflen to halfringlen
c
          WRITE(cdum2(LEN_TRIM(cdum2):128),
     .          '(A,2G8.1,A,I3,A,I3,A,F4.1)')
     .      '  '//errstr(ABS(errcode)),simag1,simag2,' (',
     .        INT(simag1 / (halfringlen - soffset) * 100.0),'-',
     .        INT(simag2 / (halfringlen - soffset) * 100.0),
     .      '%)  O: ',actswerror
        ELSEIF (ABS(errcode).GT.0) THEN
          per = (serr - soffset) / (spts(npts) - soffset) * 100.0
          WRITE(cdum2(LEN_TRIM(cdum2):128),'(A,G10.3,A,I3,A,F4.1)')
     .      '  '//errstr(ABS(errcode)),serr,' (',INT(per),
     .      '%)  O: ',actswerror
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
          IF (ABS(errcode).EQ.1) THEN
c
c           jdemod - changed halflen to halfringlen
c
            WRITE(0,'(A,I2,A,F5.1,A,I3,1X,A,1X,F6.2,A,F6.2,A)')
     .        'SOL22: High ',ir,' (',deltat,' s)  E: ',err2,
     .        errstr(INT(errcode)),
     .        simag1 / (halfringlen - soffset) * 100.0,' - ',
     .        simag2 / (halfringlen - soffset) * 100.0,'%'
          ELSE
            WRITE(0,'(A,I2,A,F5.1,A,I3,1X,A)')
     .        'SOL22: High ',ir,' (',deltat,' s)  E: ',err2,
     .        errstr(INT(errcode))
          ENDIF
        ENDIF

        WRITE(cdum2(1:128),'(A,I2,1X,2F7.3,1P,E11.3,0P,A,
     .                       F5.1,A,I3)')
     .    'SOL22 H: ',ir,te0,ti0,n0,' (',deltat,' s) E:',errcode

        IF (ABS(errcode).EQ.1) THEN
          per = (serr - soffset) / (spts(npts) - soffset) * 100.0
c
c           jdemod - changed halflen to halfringlen
c
          WRITE(cdum2(LEN_TRIM(cdum2):128),'(A,2G8.1,A,I3,A,I3,A,
     .                                       F4.1)')
     .      '  '//errstr(ABS(errcode)),simag1,simag2,' (',
     .        INT(simag1 / (halfringlen - soffset) * 100.0),'-',
     .        INT(simag2 / (halfringlen - soffset) * 100.0),
     .      '%)  O: ',actswerror
        ELSEIF (ABS(errcode).GT.0) THEN
          per = (serr - soffset) / (spts(npts) - soffset) * 100.0
          WRITE(cdum2(LEN_TRIM(cdum2):128),'(A,G10.3,A,I3,A,F4.1)')
     .      '  '//errstr(ABS(errcode)),serr,' (',INT(per),
     .      '%)  O: ',actswerror
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

c          WRITE(74,*)
c          WRITE(74,'(A)') cdum2(1:LEN_TRIM(cdum2))
c
c          CALL CloseStorageFiles
c
c          CALL InsertFile('osm1tmp.dat',PINOUT)
c          CALL InsertFile('osm6tmp.dat',PINOUT)
c          CALL InsertFile('osm2tmp.dat',PINOUT)
c          CALL InsertFile('osm3tmp.dat',PINOUT)
c          CALL InsertFile('osm5tmp.dat',PINOUT)
c          CALL InsertFile('osm4tmp.dat',PINOUT)
c
c          WRITE(PINOUT,*)
c        ENDIF
c
c jdemod - Commented these out - the code should not exit here - error STOPs
c          are handled in the calling routine - in addition - the STOP 
c          comment and the DB one are inconsistent
c
c        IF (err2.GE.6) STOP 'Fatal error: outer'
c
c        CALL DB('Done finding solution for inner half-ring')
c
c jdemod
c
      ENDIF

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: SOL22Headers
c
c
c
      SUBROUTINE SOL22Headers
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slcom'
      INCLUDE 'solparams'
      INCLUDE 'solcommon'


      IF (osm_mode.LE.1) RETURN


      IF (miter.EQ.1) THEN
        WRITE(75,*)
        WRITE(75,'(A,A3,A8,1X,2A7,A10,1X,A10,A5,2(1X,A10))')
     .    '`','in','s','Ti','Te','ne','Vb','M','Ga','P'

        IF (outmode.GE.3) THEN
          WRITE(71,*)
          WRITE(71,'(A,A3,4(1X,A10,10X))')
     .      '`','in','pais','paes','peis','srcf'
        ENDIF

        IF (forcet.EQ.0.OR.forcet.EQ.2.OR.forcet.EQ.3) THEN
          WRITE(72,*)
          WRITE(72,'(A,A3,1X,A7,3A10,10X,A9,2X,3A10)')
     .      '`','in','Ti','Pcf','Pcx','Pei','Pu/Pt',
     .      'Conv','Cond','Total'
        ENDIF

        WRITE(73,*)
        WRITE(73,'(A,A3,1X,A7,4A10,A9,2X,3A10)')
     .    '`','in','Te','Pcf','PHi','Pei','Prad','Pu/Pt',
     .    'Conv','Cond','Total'
c        WRITE(73,*)
c        WRITE(73,'(A,A3,1X,A7,3A10,A9,2X,3A10)')
c     .    '`','in','Te','Pcf','PHi','Pei','Pu/Pt',
c     .    'Conv','Cond','Total'

      ENDIF

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: NoName01
c
c
c
      SUBROUTINE NoName01(imflag, i,  loopstart,spts,npts,serr,
     .                    exitcond,  note,te,ti,ne,ne2,ga,vb,vb2,
     .                    pir,pii,vgradn,act_press,exp_press,
     .                    vsound,vsubs,vsupers,peiv,peiv2,pradv,pcxv,
     .                    phelpiv,scxv,convi,conve,condi,conde,errcode,
     .                    *)

      IMPLICIT none

      INCLUDE 'solparams'
      INCLUDE 'solcommon'

      INCLUDE 'params'
      INCLUDE 'slcom'

      integer i,i1,k,flag,imflag,lastflag,vcount
      integer loopstart
      integer errcode,npts,ir,irsep
      integer exitcond,negerrflag

      real*8 serr
      REAL*8  spts(mxspts)
      real*8  te(mxspts),ti(mxspts),ne(mxspts),vb(mxspts),
     >        exp_press(mxspts),act_press(mxspts),
     >        prad(mxspts)
      real*8    ga(mxspts),vb2(mxspts)
      real*8    ne2(mxspts),vsupers(mxspts),vsubs(mxspts)
      real*8    vsound(mxspts),scxv(mxspts)
      real*8    peiv(mxspts),pcxv(mxspts),phelpiv(mxspts)
      real*8    peiv2(mxspts),pradv(mxspts)
      real*8    pir(mxspts),pii(mxspts),vgradn(mxspts)
      real*8    condi(mxspts),conde(mxspts)
      real*8    convi(mxspts),conve(mxspts)

      CHARACTER*2 note(MXSPTS)


      IF (imflag.EQ.1) THEN
        IF (simag1.EQ.HI) THEN
          simag1 = simag - soffset
          simag2 = simag - soffset
        ELSE
          simag2 = simag - soffset
        ENDIF
      ENDIF

      IF ((i.GT.loopstart.AND.
     .    ((imflag.EQ.1.AND.spts(i).GT.osm_range*spts(npts)).OR.
     .     exitcond.GE.3)).OR.exitcond.EQ.99) THEN

c        WRITE(PINOUT,*)
c        WRITE(PINOUT,'(I4,3F10.4,2I4)')
c          i,spts(i),osm_range,osm_range*spts(npts),imflag,exitcond

        ierror = i

        DO i1 = i, npts
          IF (osm_mode.EQ.2) note(i1) = ' e'

          te (i1) = te (i-1)
          ti (i1) = ti (i-1)
          ne (i1) = ne (i-1)
          ne2(i1) = ne2(i-1)
          ga (i1) = ga (i-1)
          vb (i1) = vb (i-1)
          vb2(i1) = vb2(i-1)

          pir   (i1) = pir   (i-1)
          pii   (i1) = pii   (i-1)
          vgradn(i1) = vgradn(i-1)

          act_press(i1) = act_press(i-1)
          exp_press(i1) = exp_press(i-1)

          vsound (i1) = vsound (i-1)
          vsubs  (i1) = vsubs  (i-1)
          vsupers(i1) = vsupers(i-1)

          peiv   (i1) = peiv   (i-1)
          peiv2  (i1) = peiv2  (i-1)
          prad   (i1) = prad   (i-1)
          pradv  (i1) = pradv  (i-1)
          pcxv   (i1) = pcxv   (i-1)
          phelpiv(i1) = phelpiv(i-1)
          scxv   (i1) = scxv   (i-1)

          convi(i1) = convi(i-1)
          conve(i1) = conve(i-1)
          condi(i1) = condi(i-1)
          conde(i1) = conde(i-1)
        ENDDO

        IF (spts(i).GT.osm_range*spts(npts)) THEN
          IF ((exitcond.GE.3.AND.exitcond.LT.99).AND.
     .       snegerr.LT.simag) THEN
            errcode = exitcond
            serr    = snegerr
          ELSEIF (imflag.GT.0) THEN
            errcode = imflag
            serr    = simag
          ENDIF
          errcode = -errcode
          RETURN 1
c          GOTO 2000
        ENDIF
      ENDIF

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: SOL22Output
c
c
c
      SUBROUTINE SOL22Output(loopstart,spts,npts,conde,condi,conve,
     .                       convi,pcxv,peiv,phelpiv,pradv,te,ti,ne,
     .                       vb,ga,act_press,pmloss,note)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'slcom'
      INCLUDE 'solparams'
      INCLUDE 'solcommon'

      COMMON /OUTPUTJUNK/ dp4        ,dp6
      REAL                dp4(MXSPTS),dp6(MXSPTS)

      REAL     GetCs
      REAL*8   cond,conv,paes,pais,pmomloss,gperpf
      EXTERNAL cond,conv,paes,pais,pmomloss,gperpf

      INTEGER i,loopstart,npts
      REAL*8  srcf,powi,powe,mach,
     .        te(MXSPTS),spts (MXSPTS),pmloss(MXSPTS),exp_press(MXSPTS),
     .        ne(MXSPTS),prad (MXSPTS),pcxv  (MXSPTS),act_press(MXSPTS),
     .        ga(MXSPTS),peiv (MXSPTS),pradv (MXSPTS),phelpiv  (MXSPTS),
     .        ti(MXSPTS),condi(MXSPTS),conde (MXSPTS),
     .        vb(MXSPTS),convi(MXSPTS),conve (MXSPTS)
      CHARACTER*2 note(MXSPTS)


      IF (osm_mode.LE.1) RETURN


c     Pcx > 0 => cooling
c     PHi > 0 => cooling
c     Pei > 0 => electron cooling and ion heating

      DO i = loopstart, npts
        powi = (pai - pais(spts(i))) - pcxv   (i) + peiv(i)
        powe = (pae - paes(spts(i))) - phelpiv(i) - peiv(i)
        mach = DBLE(GetCs(SNGL(te(i)),SNGL(ti(i))))

        WRITE(70,'(1X,I3,F8.3,1X,2F7.2,1P,E10.2,
     .             1X,E10.2,0P,F5.2,1P,2(1X,E10.2),1X,2E10.2,0P,
     .             F10.4,A)')
     .    i,spts(i),
     .    ti(i),te(i),ne(i),vb(i),DABS(vb(i)/mach),ga(i),
     .    act_press(i)/ECONV,
     .    ionsrc(i),gperpf(spts(i)),
c     .      srcf(spts(i)),gperpf(spts(i)),
     .    pmloss(i),
     .    note(i)

        IF (forcet.EQ.0.OR.forcet.EQ.2.OR.forcet.EQ.3) THEN
          WRITE(72,'(1X,I3,1X,F7.2,1P,3E10.2,0P,10X,F9.3,
     .               2X,1P,3E10.2,0P,A)')
     .      i,ti(i),-(pais(spts(i))-pai),-pcxv (i), peiv(i),powi/pai,
     .      convi(i),condi(i),convi(i)+condi(i),note(i)
        ENDIF                            

c...Prad!                            
        WRITE(73,'(1X,I3,1X,F7.2,1P,4E10.2,0P,F9.3,'//
     .           ' 2X,1P,3E10.2,0P,1X,2F8.2,A)')
     .    i,te(i),-(paes(spts(i))-pae),-phelpiv(i),-peiv(i),-pradv(i),
     .    powe/pae,conve(i),conde(i),conve(i)+conde(i),dp4(i),dp6(i),
     .    note(i)

c        WRITE(73,'(1X,I3,1X,F7.2,1P,3E10.2,0P,F9.3,
c     .             2X,1P,3E10.2,0P,1X,2F8.2,A)')
c     .      i,te(i),-(paes(spts(i))-pae),-phelpiv(i),-peiv(i),powe/pae,
c     .      conve(i),conde(i),conve(i)+conde(i),dp4(i),dp6(i),note(i)

        IF (outmode.GE.3)
     .    WRITE(71,'(1X,I3,1P,4(1X,2E10.2),3X,E10.2,A)')
     .      i,
     .      pais(spts(i)),dumpai1,
     .      paes(spts(i)),dumpae1,
     .      peiv(i      ),dumpei1,
     .      srcf(spts(i)),dumsrc1,
     .      pmomloss(spts(i),1,vb(i),te(i),ti(i)),note(i)
      ENDDO


      RETURN
99    STOP
      END
c
c ======================================================================
c
c
c ======================================================================
c
c subroutine: AssignPP
c
c Assign a uniform plasma to the private flux zone from data listed in
c DIVIMP input file.
c
c
      SUBROUTINE AssignPP(irlim1,irlim2,ikopt,targ_con_opt)
      IMPLICIT none
c
c     irlim1 - start ring
c     irlim2 - end ring
c     ikopt  - apply to whole or partial ring
c     targ_con_opt   =5 - do not apply plasma values to target
c                    =6 - apply plasma values to target
c
      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER SymmetryPoint
      REAL    GetCs

      INTEGER irlim1,irlim2,ikopt,targ_con_opt

      INTEGER ik,ikm,ir,i1
      LOGICAL output
      REAL    te(2,MAXNRS),ti(2,MAXNRS),ne(2,MAXNRS),vb(2,MAXNRS)

      output = .TRUE.

      IF (osmnppv.EQ.0) CALL ER('AssignPP','No plasma data from input'//
     .                          'file found',*99)

      CALL RZero(te,MAXNRS*2)
      CALL RZero(ti,MAXNRS*2)
      CALL RZero(ne,MAXNRS*2)
      CALL RZero(vb,MAXNRS*2)

c...  Search input file listing for ring plasma data:
      DO ir = irlim1, irlim2
        DO i1 = 1, osmnppv
          IF (INT(osmppv(i1,1)).EQ.ir) THEN
            IF (osmppv(i1,2).EQ.1.0.OR.osmppv(i1,2).EQ.3.0) THEN
              te(IKLO,ir) =  osmppv(i1,3)
              ti(IKLO,ir) =  osmppv(i1,4)
              ne(IKLO,ir) =  osmppv(i1,5)
              vb(IKLO,ir) = -osmppv(i1,6)
            ENDIF
            IF (osmppv(i1,2).EQ.2.0.OR.osmppv(i1,2).EQ.3.0) THEN
              te(IKHI,ir) =  osmppv(i1,3)
              ti(IKHI,ir) =  osmppv(i1,4)
              ne(IKHI,ir) =  osmppv(i1,5)
              vb(IKHI,ir) =  osmppv(i1,6)
            ENDIF
          ENDIF
        ENDDO
      ENDDO

c...  Assign plasma to rings with Thomson data:
      DO ir = irlim1, irlim2
        ikm = SymmetryPoint(ir)

        IF (ikopt.EQ.1.OR.ikopt.EQ.3) THEN
          IF (te(IKLO,ir).EQ.0.0)
     .      CALL ER('AssignPP','Required low index plasma data not '//
     .                         'found',*99)
          DO ik = 1, ikm
            ktebs(ik,ir) = te(IKLO,ir)
            ktibs(ik,ir) = te(IKLO,ir)
            knbs (ik,ir) = ne(IKLO,ir)
            kvhs (ik,ir) = vb(IKLO,ir)
          ENDDO
        ENDIF
        IF (ikopt.EQ.2.OR.ikopt.EQ.3) THEN
          IF (te(IKHI,ir).EQ.0.0)
     .      CALL ER('AssignPP','Required high index plasma data not '//
     .                         'found',*99)
          DO ik = ikm+1, nks(ir)
            ktebs(ik,ir) = te(IKHI,ir)
            ktibs(ik,ir) = te(IKHI,ir)
            knbs (ik,ir) = ne(IKHI,ir)
            kvhs (ik,ir) = vb(IKHI,ir)
          ENDDO
        ENDIF

c...    Assign target values using Thomson data:
        IF (targ_con_opt.EQ.6) THEN
          IF (ikopt.EQ.2.OR.ikopt.EQ.3) THEN
            kteds(idds(ir,1)) = ktebs(nks(ir),ir) * te_mult_i
            ktids(idds(ir,1)) = ktibs(nks(ir),ir) * ti_mult_i
            IF (n_mult_i.GT.10) THEN
c...          Sometimes Isat or Ne data in the input file uses a multiplier
c             so that the numbers listed in the target data arrays do not
c             have to include the exponents (e.g. 1.0E+20 can be entered as
c             1.0 if a 1.0E+20 multiplier is being used).  However,
c             such vile conduct will break this code, which is assigning
c             the target data from the bulk plasma prescription but also applies
c             the target multiplier:
              CALL ER('AssignPP','Suspicious high index target data '//
     .                           'multiplier detected',*99)
            ELSE
              knds(idds(ir,1)) = knbs(nks(ir),ir) * n_mult_i
            ENDIF
            kvds (idds(ir,1)) =
     .              GetCs(kteds(idds(ir,1)),ktids(idds(ir,1)))
          ENDIF
          IF (ikopt.EQ.1.OR.ikopt.EQ.3) THEN
            kteds(idds(ir,2)) = ktebs(1,ir) * te_mult_o
            ktids(idds(ir,2)) = ktibs(1,ir) * ti_mult_o
            IF (n_mult_o.GT.10) THEN
c...          See the above comment for the high index target:
              CALL ER('AssignPP','Suspicious low index target data '//
     .                           'multiplier detected',*99)
            ELSE
              knds(idds(ir,2)) = knbs(1,ir) * n_mult_o
            ENDIF
            kvds (idds(ir,2)) =
     .             -GetCs(kteds(idds(ir,2)),ktids(idds(ir,2)))
          ENDIF
        ENDIF
      ENDDO

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: CloseStorageFiles
c
      SUBROUTINE CloseStorageFiles
      IMPLICIT   none

      CALL DB('Closing storage files')

      CLOSE(70)
      CLOSE(71)
      CLOSE(72)
      CLOSE(73)
      CLOSE(74)
      CLOSE(75)

      RETURN
      END
c
c ======================================================================
c
c subroutine: OpenStorageFiles
c
      SUBROUTINE OpenStorageFiles(ir,target,ta)
      IMPLICIT   none

      INCLUDE 'params'
      INCLUDE 'slcom'

      INTEGER   ir,target
      CHARACTER ta*(*)

c      INTEGER rir,rtarg,rstep,riter
c      DATA    rir,rtarg,rstep,riter /-1, -1, -1, -1/
c      SAVE

      CALL DB('Opening storage files')

c      IF (rir  .NE.ir      .OR.rtarg.NE.target  .OR.
c     .    rstep.NE.rel_step.OR.riter.NE.rel_iter) THEN
c        OPEN(UNIT=74,FILE='osm1'//ta,ACCESS='SEQUENTIAL',
c     .       STATUS='REPLACE')
c      ELSE
c        OPEN(UNIT=74,FILE='osm1'//ta,ACCESS='SEQUENTIAL',
c     .       STATUS='OLD',POSITION='APPEND')
c      ENDIF

      OPEN(UNIT=74,FILE='osm1'//ta,ACCESS='SEQUENTIAL',STATUS='REPLACE')
      OPEN(UNIT=70,FILE='osm2'//ta,ACCESS='SEQUENTIAL',STATUS='REPLACE')
      OPEN(UNIT=71,FILE='osm3'//ta,ACCESS='SEQUENTIAL',STATUS='REPLACE')
      OPEN(UNIT=72,FILE='osm4'//ta,ACCESS='SEQUENTIAL',STATUS='REPLACE')
      OPEN(UNIT=73,FILE='osm5'//ta,ACCESS='SEQUENTIAL',STATUS='REPLACE')
      OPEN(UNIT=75,FILE='osm6'//ta,ACCESS='SEQUENTIAL',STATUS='REPLACE')

c      rir   = ir
c      rtarg = target
c      rstep = rel_step
c      riter = rel_iter

      RETURN
      END
c
c ======================================================================
c
c subroutine: InsertFile
c
      SUBROUTINE InsertFile(fname,fp2)
      IMPLICIT   none

      INTEGER   fp2
      CHARACTER fname*(*),buffer*256

      CALL DB('Inserting file')

      OPEN(UNIT=70,FILE=fname,STATUS='OLD',ACCESS='SEQUENTIAL',ERR=15)

10    CALL ReadLine(70,buffer,1,*15,*15)
      WRITE(fp2,'(A)') buffer(1:LEN_TRIM(buffer))
      GOTO 10
15    CONTINUE

      CLOSE(70)

      RETURN
      END






c
c ======================================================================
c
c
c ======================================================================
c
c subroutine: CalcInitSrc
c
      SUBROUTINE CalcInitSrc(region,ir)
      IMPLICIT none

      INCLUDE 'solparams'
      INCLUDE 'solswitch'
      INCLUDE 'solcommon'
      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

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


      CALL DB('Calculating initial PIN source estimates')

c      WRITE(0,*) '*** NOT SETTING PININIT = .TRUE. ***'
c      IF (region.EQ.IKHI) pininit(ir) = .TRUE.

c      WRITE(PINOUT,*)
c      WRITE(PINOUT,*) 'Generating initial sources:'

      IF (region.EQ.IKLO) THEN
        iks     =  ikbound(ir,IKLO)
        ike     =  SymmetryPoint(ir)
        parflux = -gtarg(ir,2)
      ELSE
        iks     =  SymmetryPoint(ir) + 1
        ike     =  ikbound(ir,IKHI)
        parflux = -gtarg(ir,1)
      ENDIF


c... REMEMBER TO UNDO ADJUSTMENT IN SOLASCV1 ALSO!
c      IF (osm_powopt.EQ.2) WRITE(0,*) 'INITIAL PINQE AND PINQI AT 10%'



      DO ik = iks, ike
        IF (region.EQ.IKLO) THEN
          s = kss2(ik,ir)
        ELSE
          s = ksmaxs2(ir) - kss2(ik,ir)
        ENDIF

c...    Generate momentum source:
        IF (actswnmom.EQ.11.0.OR.actswnmom.EQ.12.0) THEN
        


          IF (s.LT.lenmom*ringlen) THEN
            IF (region.EQ.IKLO.AND.osm_model(IKLO,ir).EQ.24.AND.
     .          (iflexopt(6).EQ.20.OR.iflexopt(6).EQ.21.OR.
     .           iflexopt(6).EQ.22.OR.iflexopt(6).EQ.24)) THEN  
              pinmp(ik,ir) = smom0 * EXP(-(s - soffset) / 
     .                                     (3.0*lammom * ringlen))
            ELSE
              pinmp(ik,ir) = smom0 * EXP(-(s - soffset) / 
     .                                     (lammom * ringlen))
            ENDIF

c            WRITE(PINOUT,* ) '--pinmp--?',ik,smom0,
c     .               EXP(-(s - soffset) / (lammom * ringlen)),
c     .               pinmp(ik,ir)        
          ELSE
            pinmp(ik,ir) = 0.0
          ENDIF

          IF (region.EQ.IKHI) pinmp(ik,ir) = -pinmp(ik,ir)

          osmmp(ik,ir) = pinmp(ik,ir)


        ENDIF
c
c       Generate ionisation source:
c
        IF     (actswioni.EQ.0.0) THEN
          IF (osm_matcht.GT.0.AND.osm_powopt.EQ.1) THEN
            pinion(ik,ir) = 0.0
            IF (region.EQ.IKLO.AND.ik.EQ.iks.OR.
     .          region.EQ.IKHI.AND.ik.EQ.ike) THEN
              pinion(ik,ir) = 0.0
            ELSE
              IF (region.EQ.IKLO.AND.ik.EQ.ike.OR.
     .            region.EQ.IKHI.AND.ik.EQ.iks) THEN
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
          pinion(ik,ir) = s0 * (s**5.0) *
     .                    EXP(-s5alph * ((s + 0.5 * s5gausslen)**2.0))
        ELSE
          CALL ER('CalcInitSrc','ACTSWIONI value not supported',*99)
        ENDIF
c
c       Generate Phelpi:
c
        IF (actswphelp.EQ.1.0.OR.actswphelp.EQ.2.0) THEN
          IF (knbs(ik,ir).GT.1.0E+20) THEN
c...fix
            helpi = 17.5 + (5.0 + 37.50 / ktebs(ik,ir)) *
     .                     (1.0 +  0.25 / ktebs(ik,ir))
          ELSE
            helpi = 17.5 + (5.0 + 37.50 / ktebs(ik,ir)) *
     .                     (1.0 +  0.25 / ktebs(ik,ir)) *
     .              LOG10(1.0E+21 / knbs(ik,ir))
          ENDIF

c          WRITE(PINOUT,*) '   ',ik,ir,ktebs(ik,ir),ktebs(ik,ir),
c     .                         knbs(ik,ir),LOG10(1.0E+21 / knbs(ik,ir))

          IF     (osm_matcht.EQ.0) THEN
            pinqe(ik,ir) = -helpi * pinion(ik,ir) * ECH
          ELSEIF (osm_powopt.EQ.2) THEN
            pinqe(ik,ir) = -helpi * pinion(ik,ir) * ECH
c            pinqe(ik,ir) = 0.1 * pinqe(ik,ir)
            pinqe(ik,ir) = MIN(-1.0,pinqe(ik,ir))
          ELSE
            pinqe(ik,ir) = 0.0
            IF (region.EQ.IKLO.AND.ik.EQ.iks.OR.
     .          region.EQ.IKHI.AND.ik.EQ.ike) THEN
               pinqe(ik,ir) = 0.0
            ELSEIF (region.EQ.IKLO.AND.ik.EQ.ike.OR.
     .              region.EQ.IKHI.AND.ik.EQ.iks) THEN
               pinqe(ik,ir) = 0.0
            ELSE
c
c             jdemod - this assignment to pinqe would appear to make no
c                      sense since it is assigning a length value to pinqe
c                      presumably the length should be multiplied by some
c                      sort of source term
c
              pinqe(ik,ir) = -MAX(SNGL(0.15*ringlen-(s-soffset))
     .                            ,0.0)
c              pinqe(ik,ir) = -MAX(SNGL(0.95*(halflen-soffset)-(s-soffset))
c     .                            ,0.0)
c              pinqe(ik,ir) = -ABS(SNGL((s-soffset)-0.5*(halflen-soffset))
            ENDIF
          ENDIF
        ELSEIF (actswphelp.EQ.0.0) THEN
          pinqe(ik,ir) = 0.0
        ELSE
          CALL ER('CalcInitSrc','ACTSWPHELP value not supported',*99)
        ENDIF
c
c       Generate Pcx:
c
        IF     (actswpcx.EQ.1.0.OR.actswpcx.EQ.2.0.OR.
     .          actswpcx.EQ.5.0) THEN
          pinqi(ik,ir) = -1.5 * ktibs(ik,ir) * CEICF * pinion(ik,ir) *
     .                    ECH
c          pinqi(ik,ir) = pinqi(ik,ir) * 0.1
        ELSEIF (actswpcx.EQ.0.0) THEN
          pinqi(ik,ir) = 0.0
        ELSE
          CALL ER('CalcInitSrc','ACTSWPCX value not supported',*99)
        ENDIF

c...prad!
      IF     (actswprad.EQ.0.0) THEN
        osmrad(ik,ir) = 0.0
      ELSEIF (actswprad.EQ.1.0) THEN
        IF (s.LT.lenr) THEN
          osmrad(ik,ir) = lamr * prad0 * (1 - exp(-s / lamr))
        ELSE
          osmrad(ik,ir) = 0.0
        ENDIF
      ELSEIF (actswprad.EQ.3.0) THEN
c...    NOT VALID AT THE MOMENT SINCE THIS RADIATION OPTION IS 
c       NOT APPLIED WITH THE INITIAL SOLUTION:
c       osmrad(ik,ir) = radsrc_mult * pinqe(ik,ir)
      ELSE
        IF (rel_opt.NE.0.AND.rel_frac.NE.1.0) 
     .    CALL WN('CalcInitSrc','ACTSWPRAD value not supported')
      ENDIF

      ENDDO

c...bug! (IF statement wasn't there and SOL22 rings didn't have any
c         ionisation in the first cell, which exaggerated imaginary
c         solutions.)
      IF (GetModel(region,ir).EQ.24) THEN
        IF (region.EQ.IKLO) THEN
          pinion(iks,ir) = 0.0
        ELSE
          pinion(ike,ir) = 0.0
        ENDIF
      ENDIF
c
c     Recombination source for SOL22p:
c
      IF     (actswrecom.EQ.1.0) THEN
        IF (GetModel(region,ir).EQ.24) Call CalcInitRecom(region,ir)
      ELSEIF (actswrecom.EQ.0.0) THEN
c...  Do nothing.
      ELSE
        CALL ER('CalcInitSrc','ACTSWRECOM value not supported',*99)
      ENDIF
c
c     Calculate momentum loss profile:
c
c...  Assume half the pressure is lost:
c      IF (iflexopt(4).EQ.20) THEN
c...Need flag set so that I know data is available...
c        pu = CalcPressure(prp_ne(ir,i1),prp_te(ir,i1),
c     .                    prp_ti(ir,i1),dip_v (ir,i1))
c
c        pl = pu * 0.5
c
c        CALL CalcIntegral4(pinion,iks,ike,ir,integ,2)
c
c        fact = pl / integ
c
c        DO ik = iks, ike
c          osmmp(ik,ir) = fact * ABS(pinion(ik,ir))
c        ENDDO
c      ENDIF


c...  Whipe pinqi for now, since it is causing trouble on the 24-24 rings.  What happens
c     is the global mock power multiplier ends up changing sign when pinqi terms grows,
c     and this wreaks havok on the local mock multipliers:
      IF ((stopopt3.EQ.17.OR.stopopt3.EQ.18.OR.stopopt3.EQ.19.OR.
     .     stopopt3.EQ.20).AND.
     .    osm_model(IKLO,ir).EQ.24.AND.osm_model(IKHI,ir).EQ.24) THEN
        WRITE(0,*) 'BLANKING PINQI ON RING',ir
        CALL RZero(pinqi(1,ir),MAXNKS)
      ENDIF


      IF ((iflexopt(6).EQ.20.OR.iflexopt(6).EQ.22.OR.
     .     iflexopt(6).EQ.24).AND.
     .    osm_model(IKLO,ir).EQ.24) THEN
c....   Get rid of PINQI for now:
        WRITE(PINOUT,*) 
        WRITE(PINOUT,*) '***********************'
        WRITE(PINOUT,*) '  WHIPING INNER Qi ON ',ir
        WRITE(PINOUT,*) '***********************'
        WRITE(PINOUT,*) 
        DO ik = 1, osm_sympt(ir)
          pinqi(ik,ir) = 0.0
        ENDDO
      ENDIF

      IF ((iflexopt(6).EQ.22).AND.
     .    osm_model(IKHI,ir).EQ.22) THEN
c....   Get rid of PINQI for now:
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
99    STOP
      END
c
c ======================================================================
c
c subroutine: CalcInitRecom
c
      SUBROUTINE CalcInitRecom(region,ir)
      IMPLICIT none

      INTEGER region,ir

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'cadas'
      INCLUDE 'slcom'

      INTEGER GetModel
      REAL    GetEAD

      INTEGER in,i1,ik,iks,ike
      REAL    maxcoef,rdum1,temin,nemax,te,ne

      temin = 0.0
      nemax = HI

      IF (osm_preopt.GT.0.AND.
     .    (GetModel(region,ir).EQ.22.OR.
     .     GetModel(region,ir).EQ.24)) THEN

        IF     (region.EQ.IKLO) THEN
          iks = 1
          ike = osm_sympt(ir)
c          ike = ikbound(ir,IKLO) - 1

          temin = ktebs(MAX(1,ikbound(ir,IKLO)-1),ir)
          nemax = knbs (MAX(1,ikbound(ir,IKLO)-1),ir)
        ELSEIF (region.EQ.IKHI) THEN
          iks = osm_sympt(ir) + 1
c          iks = ikbound(ir,IKHI) + 1
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
      ENDIF

c      WRITE(PINOUT,*) 'RECOM IKS,IKE = ',region,iks,ike

      IF (.FALSE.) THEN
c...    Use ADAS recombination data:

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

      ELSEIF (eiropacity.EQ.0.OR.eiropacity.EQ.1.OR.
     .        eiropacity.EQ.2.OR.eiropacity.EQ.6.OR.
     .        eiropacity.EQ.-3) THEN
c...    Use Lyman transparent AMJUEL (EIRENE database) recombination data:
        DO ik = iks, ike
          te = MAX(temin,ktebs(ik,ir))
          ne = MIN(nemax,knbs (ik,ir))
          pinrec(ik,ir) = GetEAD(te,ne,3,'H.4 ') * 
     .                    1.0E-06 * knbs(ik,ir) * knbs(ik,ir) 
c          pinrec(ik,ir) = GetEAD(te,ne,3,'H.4 ') * 
c     .                    1.0E-06 * knbs(ik,ir) * knbs(ik,ir) * 
c     .                    eirscale(11)
        ENDDO

      ELSEIF ((eiropacity.EQ.3.OR.eiropacity.EQ.4).AND.
     .        tagpinatom.AND.tagmulrec) THEN
c...    Use non-Lyman opaque data with local multiplier:

c...    Not sure this belongs here.  We are really looking for an initial guess here, and
c       it may not be appropriate to be using PINATOM data, even if it is available...

        STOP 'CALCINITRECOM THIS SHOULD NOT BE CALLED'

        DO ik = iks, ike
          te = MAX(temin,ktebs(ik,ir))
          ne = MIN(nemax,knbs (ik,ir))
          pinrec(ik,ir) = GetEAD(te,ne,3,'H.4 ') * 
     .                    1.0E-06 * knbs(ik,ir) * knbs(ik,ir) * 
     .                    mulrec(ik,ir)
c          pinrec(ik,ir) = GetEAD(te,ne,3,'H.4 ') * 
c     .                    1.0E-06 * knbs(ik,ir) * knbs(ik,ir) * 
c     .                    eirscale(11) * mulrec(ik,ir)
        ENDDO

      ELSEIF (eiropacity.EQ.3.OR.eiropacity.EQ.4.OR.
     .        (eiropacity.EQ.5.AND.s28mode.LT.2.0).OR.
c *NOT QUITE RIGHT*
     .        eiropacity.EQ.-4.OR.eiropacity.EQ.-5) THEN
c...    Use Lyman opaque AMJUEL (EIRENE database) recombination data:
        DO ik = iks, ike
          te = MAX(temin,ktebs(ik,ir))
          ne = MIN(nemax,knbs (ik,ir))
          pinrec(ik,ir) = GetEAD(te,ne,4,'H.4 ') * 
     .                    1.0E-06 * knbs(ik,ir) * knbs(ik,ir) 
        ENDDO

      ELSEIF (eiropacity.EQ.5) THEN
c
c ***NOT IN USE BECAUSE THE OLD SOL24/28 CRASHES***
c

c...    Use Lyman alpha opaque AMJUEL (EIRENE database) recombination data:
        DO ik = iks, ike
          te = MAX(temin,ktebs(ik,ir))
          ne = MIN(nemax,knbs (ik,ir))
          pinrec(ik,ir) = GetEAD(te,ne,27,'H.4 ') * 
     .                    1.0E-06 * knbs(ik,ir) * knbs(ik,ir) 
        ENDDO

      ELSE
        CALL ER('CalcInitRecom','Unrecognized option',*99)
      ENDIF


      RETURN
99    WRITE(0,*) 'IR,GETMODEL=',ir,region,GetModel(region,ir)
      STOP
      END
c
c ======================================================================
c







