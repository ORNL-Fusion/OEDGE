c     -*-Fortran-*- 
c
c ======================================================================
c
c subroutine: SetupToroidalSurfaces
c
c
      SUBROUTINE ReadSOLPSDataFile(filename,array)
      USE mod_grid_divimp
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      CHARACTER filename*(*)
      REAL      array(MAXNKS,MAXNRS)

      INTEGER    MAXIK    ,MAXIR
      PARAMETER (MAXIK=200,MAXIR=200)

      INTEGER i1,ik,ik1,ir,ir1,fp,ikindex(MAXIK),irindex(MAXIR),
     .        numik,numir
      CHARACTER buffer*4096

      REAL, ALLOCATABLE :: b2dat(:,:)

      ALLOCATE(b2dat(MAXIK,MAXIR))

      fp = 98

      WRITE(0,*) 'LOADING SOLPS DATA '//filename(1:LEN_TRIM(filename))
c...  
      OPEN(UNIT=fp,FILE=filename(1:LEN_TRIM(filename)),
     .     ACCESS='SEQUENTIAL',STATUS='OLD',ERR=90)
 
      READ(fp,'(A4096)') buffer
      numik = 0
      DO WHILE(.TRUE.) 
        numik = numik + 1
        READ(buffer,'(200(I14))',ERR=10) (ikindex(i1),i1=1,numik)
      ENDDO
 10   CONTINUE
      numik = numik - 1

c      WRITE(0,*) ikindex(1:numik)

      numir = 0
      DO WHILE(.TRUE.) 
        numir = numir + 1
        READ(fp,'(I4,1X,200(E14.6))',ERR=20,END=20) 
     .    irindex(numir),(b2dat(i1,numir),i1=1,numik)
      ENDDO
 20   CONTINUE
      numir = numir - 1

c      WRITE(0,*) 'IR:',irindex(1:numir)

c...  Shift to match SONNET grid index:
      ikindex(1:numik) = ikindex(1:numik) + 1
      irindex(1:numir) = irindex(1:numir) + 1

c...  Assign B2 data:

      array = 0.0

      DO ir = 1, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE

        array(1:MAXNKS,ir) = -999.0

        DO ik = 1, nks(ir)     


          IR_LOOP: DO ir1 = 1, numir
            DO ik1 = 1, numik
              IF (divimp_ik(ik,ir).EQ.ikindex(ik1).AND.
     .            divimp_ir(ik,ir).EQ.irindex(ir1)) THEN    
                array(ik,ir) = b2dat(ik1,ir1)
                EXIT IR_LOOP
              ENDIF
            ENDDO
          ENDDO IR_LOOP

          IF (array(ik,ir).EQ.-999.0) THEN
            WRITE(0,*) 'PROBS:',ik,ir,irwall
          ENDIF

        ENDDO
      ENDDO

      CLOSE (fp)

      DEALLOCATE (b2dat)

      DO ir = 1, nrs
        DO ik = 1, nks(ir)
          WRITE(SLOUT,*) 'GRIDMAP:',ik,ir,' * ',
     .                   divimp_ik(ik,ir),divimp_ir(ik,ir)
        ENDDO
      ENDDO

      RETURN
 90   CALL ER('LoadSOLPSData','Sna_ion1 file not found',*99)
 99   STOP
      END

c ======================================================================
c
c subroutine: LoadSOLPSData
c
c
      SUBROUTINE LoadSOLPSData
      USE mod_grid

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'cedge2d'
      INCLUDE 'slcom'


      REAL, ALLOCATABLE :: fluvol(:,:)  ! Necessary, can I just do FLUVOL(MAXNKS,MAXNRS) ?

      ALLOCATE(fluvol(MAXNKS,MAXNRS))

      CALL LoadSOLPSData_OSM
      RETURN

c...  Files from Sergey (SOLPS5.0, 007666 at 220 ms, drifts and no drifts):
c	Sna_ion0  ionization source for neutrals
c	Sna_ion1  ionization source for ions
c	Sna_rec0  recombination source for neutrals
c	Sna_rec1  recombination source for ions
c	Smq_fr0   plasma momentum source from neutral-plasma friction for neutrals
c	Smq_fr1   plasma momentum source from neutral-plasma friction for ions
c	Shi_ion   ion energy source due to ionization
c	Shi_rec   ion energy source due to recombination
c	She_loss  source term for electron heat loss
c	fmax0     parallel momentum flux of neutrals between (ix-1,iy) and (ix,iy) cells.
c	fmax1     parallel momentum flux of ions???  between (ix-1,iy) and (ix,iy) cells.
c	fmay0     parallel momentum flux of neutrals between  (ix,iy-1) and (ix,iy) cells.
c	fmay1     parallel momentum flux of ions???  between  (ix,iy-1) and (ix,iy) cells.
c	fnax0     flux of neutrals between the (ix-1,iy) and (ix,iy) cells.
c	fnax1     flux of ions between the (ix-1,iy) and (ix,iy) cells.                      Pretty sure this is parallel flux
c	fnay0     flux of neutrals between the (ix,iy-1) and (ix,iy) cells.
c	fnay1     flux of ions between the (ix,iy-1) and (ix,iy) cells.
c	na0       density of neutrals on the (ix,iy) cell
c	na1       density of ions on the (ix,iy) cell
c	te        electron  temperature on the (ix,iy) cell
c	ti        ion  temperature on the (ix,iy) cell
c	uacc0     parallel velocity of neutrals on the (ix,iy) cell
c	uacc1     parallel velocity of ions in the (ix,iy) cell 
c	fhex      electron heat flux between the (ix-1,iy) and (ix,iy)  cells
c	fhey      electron heat flux between the (ix,iy-1) and (ix,iy)  cells
c	fhix      ion  heat flux between the (ix-1,iy) and (ix,iy)  cells
c	fhiy      ion heat flux between the (ix,iy-1) and (ix,iy)  cells

c...  
      CALL ReadSOLPSDataFile('vol.dat',fluvol)

      CALL ReadSOLPSDataFile('Sna_rec0.dat',e2drec)
      CALL ReadSOLPSDataFile('Sna_ion1.dat',e2dion)

      CALL ReadSOLPSDataFile('Smq_fr1.dat',e2dpepara)

      e2drec    = e2drec    / fluvol 
      e2dion    = e2dion    / fluvol 
      e2dpepara = e2dpepara / fluvol 

      CALL ReadSOLPSDataFile('na1.dat'  ,e2dnbs)
      CALL ReadSOLPSDataFile('uacc1.dat',e2dvhs)
      CALL ReadSOLPSDataFile('ti.dat'   ,e2dtibs)
      CALL ReadSOLPSDataFile('te.dat'   ,e2dtebs)

      e2dtibs = e2dtibs / ECH
      e2dtebs = e2dtebs / ECH

c...  Testing:
      CALL ReadSOLPSDataFile('fnax1.dat',e2dcxrec)
      e2dcxrec = e2dcxrec 

c      WRITE(0,*) 'E2DION:',e2dion

      DEALLOCATE(fluvol)

      WRITE(0,*) 'DONE'

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: SetupToroidalSurfaces
c
c
      SUBROUTINE SetupToroidalSurfaces
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      INTEGER i1,i2
      REAL    z1,z2,zlen

      DO i1 = 1, eirnasdat
        IF (eirasdat(i1,1).NE.6.0) CYCLE

        IF (i1.EQ.eirnasdat.OR.eirasdat(i1+1,1).NE.-1.0) 
     .    CALL ER('SetupToroidalSurfaces','Standard grid switching '//
     .            'surfaces requires continuation line',*99)

        IF (eirasdat(i1,8).EQ.eirasdat(i1  ,9).AND.
     .      eirasdat(i1,8).EQ.eirasdat(i1+1,8).AND.       
     .      eirasdat(i1,8).EQ.eirasdat(i1+1,9)) THEN

          eirnsdtor = eirnsdtor + 1
          eirsdtor(eirnsdtor) = eirasdat(i1,8)

          IF (i1+2.LE.eirnasdat.AND.eirasdat(i1+2,1).EQ.-2.0) THEN
c...        Add regularly spaced toroidal switching surfaces:
            DO i2 = 1, NINT(eirasdat(i1+2,2))
              eirnsdtor = eirnsdtor + 1
              eirsdtor(eirnsdtor) = eirsdtor(eirnsdtor-1) + 
     .                              eirasdat(i1+2,3)              
            ENDDO
          ENDIF

          IF (eirntorseg.NE.0) THEN
c...        Convert to angular coordinates for the toroidal 
c           approximation:
            DO i2 = 1, eirnsdtor
              eirsdtor(i2) = -eirsdtor(i2) / eirzaa * 360.0
            ENDDO
          ENDIF
        ELSE
          CALL ER('SetupToroidalSurfaces','All polygon verticies must'//
     .            ' have the same toroidal coordinate',*99)
        ENDIF
        
      ENDDO


      IF (eirnsdtor.GT.1) THEN

        IF     (eirzaa.EQ.-1.0) THEN
c...      Cylindrical approximation with the "toroidal circumference" based
c         on the x-point r-coordinate:
          zlen = 2.0 * PI * rxp * eirtorfrac
        ELSEIF (eirzaa.LT.0.0) THEN
c          CALL ER('SetupToroidalSurfaces','Length in z-direction not '//
c     .            'compatible with 3D standard grid',*99)
          zlen = 360.0 * eirtorfrac
        ELSE
          zlen = eirzaa
        ENDIF

      
        DO i1 = 1, eirnsdtor
          z1 = eirsdtor(i1)
          IF (i1.LT.eirnsdtor) THEN 
            z2 = eirsdtor(i1+1)
          ELSE
            z2 = zlen
          ENDIF

          eirsdvol(i1) = (z2 - z1) / zlen

          WRITE(SLOUT,'(A,I6,F12.5)') 'TOROIDAL SURFACES: ',
     .      i1,eirsdtor(i1),eirsdvol(i1)

        ENDDO

      ENDIF



      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: AssignAdditionalCellPlasma
c
c
      SUBROUTINE AssignAdditionalCellPlasma
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      INTEGER i1,i2,i3,ncut,species,cell,cut1
      REAL    xcen,ycen,xsidemin,xsidemax,ysidemin,ysidemax,frac,te1,
     .        te2,ne1,ne2,xp1,xp2,neval,teval

      DO i1 = 1, vacnpla

        IF     (vacpla(i1,1).EQ.-4.0) THEN
c...      Plasma data interpolated in the R-Z plane:

c...      Determine if the supplimental data lines are there:
          IF (vacnpla.LT.i1+3.OR.vacpla(i1+1,1).NE.0.0.OR.
     .                           vacpla(i1+2,1).NE.0.0.OR.
     .                           vacpla(i1+3,1).NE.0.0) 
     .      CALL ER('AssignAdditionalCellPlasma','Supplimental data '//
     .              'lines not found',*99)

c...      Check if a cell is inside the bounding box:
c          WRITE(0,*) 'ASC_NCELL=',asc_ncell

          DO cell = 1, asc_ncell

c...        Only include cells in the specified vaccum grid region:
            IF (NINT(vacpla(i1,8)).NE.asc_region(cell)) CYCLE

c...        Find approximate cell center:
            xsidemin =  HI
            xsidemax = -HI
            ysidemin =  HI
            ysidemax = -HI
            DO i2 = 1, ascnvertex(cell)
              xsidemin = MIN(xsidemin,ascvertex(2*i2-1,cell))
              xsidemax = MAX(xsidemax,ascvertex(2*i2-1,cell))
              ysidemin = MIN(ysidemin,ascvertex(2*i2  ,cell))
              ysidemax = MAX(ysidemax,ascvertex(2*i2  ,cell))
            ENDDO
            xcen = 0.5 * (xsidemin + xsidemax)
            ycen = 0.5 * (ysidemin + ysidemax)

            IF (xcen.GT.vacpla(i1  ,2).AND.xcen.LT.vacpla(i1+1,2).AND.
     .          ycen.GT.vacpla(i1+2,4).AND.ycen.LT.vacpla(i1  ,4)) THEN

c...          Calculate plasma quantities by interpolating between
c             the plasma data points:

c...          Vertical interpolation:
              frac = (ycen         - vacpla(i1+2,5)) / 
     .               (vacpla(i1,5) - vacpla(i1+2,5))
              te1 = (1.0-frac) * vacpla(i1+2,6) + frac * vacpla(i1,6)
              ne1 = (1.0-frac) * vacpla(i1+2,7) + frac * vacpla(i1,7)
              xp1 = (1.0-frac) * vacpla(i1+2,3) + frac * vacpla(i1,3)

              frac = (ycen           - vacpla(i1+3,5)) / 
     .               (vacpla(i1+1,5) - vacpla(i1+3,5))
              te2 = (1.0-frac) * vacpla(i1+3,6) + frac * vacpla(i1+1,6)
              ne2 = (1.0-frac) * vacpla(i1+3,7) + frac * vacpla(i1+1,7)
              xp2 = (1.0-frac) * vacpla(i1+3,3) + frac * vacpla(i1+1,3)

c...          Horizontal interpolation:
              frac = (xcen - xp1) / (xp2 - xp1)

c...          Limit FRAC for horizontal interpolation for now because it is generating
c             large T's when extrapolating to the wall:
              frac = MAX(0.0,MIN(1.0,frac))

              teval = (1.0 - frac) * te1 + frac * te2
              neval = (1.0 - frac) * ne1 + frac * ne2

              IF (cell.EQ.780.OR.cell.EQ.781.OR.cell.EQ.792.OR.
     .            cell.EQ.793) 
     .          WRITE(0,*) '2222:',cell,frac,teval

              DO cut1 = 1, ascncut

                i2 = cell + asc_ncell * (cut1 - 1) + eirnpgdat + 1

c...            The plasma quantities for the first ion species
c               are specified here, which has traditionally been
c               hydrogen ions.  If this changes in the future, 
c               then this code may need to be revised:
                species = 1 + 4

                pinasd(i2,1,species,1) = neval * 1.0E-6
                pinasd(i2,2,species,1) = teval

                pinasd(i2,3,species,1) = 0.0
                pinasd(i2,4,species,1) = 0.0
                pinasd(i2,5,species,1) = 0.0

                IF (cut1.EQ.1) 
     .            WRITE(PINOUT,'(A,2I6,3F10.2,1P,E10.2)')
     .              'ASSIGNING:',cell,cut1,xcen,ycen,teval,neval
              ENDDO

            ENDIF
          ENDDO

        ELSEIF (vacpla(i1,1).EQ.-2.0.OR.vacpla(i1,1).EQ.-3.0) THEN
c...      Assign vaccum grid plasma to a specified x-y region:
          DO cell = 1, asc_ncell

c...        Only include cells in the specified vaccum grid region:
            IF (vacpla(i1,1).EQ.-3.0.AND.
     .          NINT(vacpla(i1,8)).NE.asc_region(cell)) CYCLE

c...        Find approximate cell center:
            xsidemin =  HI
            xsidemax = -HI
            ysidemin =  HI
            ysidemax = -HI
            DO i2 = 1, ascnvertex(cell)
              xsidemin = MIN(xsidemin,ascvertex(2*i2-1,cell))
              xsidemax = MAX(xsidemax,ascvertex(2*i2-1,cell))
              ysidemin = MIN(ysidemin,ascvertex(2*i2  ,cell))
              ysidemax = MAX(ysidemax,ascvertex(2*i2  ,cell))
            ENDDO
            xcen = 0.5 * (xsidemin + xsidemax)
            ycen = 0.5 * (ysidemin + ysidemax)

            IF (xcen.GT.vacpla(i1,2).AND.xcen.LT.vacpla(i1,3).AND.
     .          ycen.GT.vacpla(i1,4).AND.ycen.LT.vacpla(i1,5)) THEN

              DO cut1 = 1, ascncut

                i2 = cell + asc_ncell * (cut1 - 1) + eirnpgdat + 1

c...            The plasma quantities for the first ion species
c               are specified here, which has traditionally been
c               hydrogen ions.  If this changes in the future, 
c               then this code may need to be revised:
                species = 1 + 4

                IF     (vacpla(i1,7).EQ.-1.0) THEN
                  pinasd(i2,1,species,1) = osmbulkn*1.E-6
                ELSEIF (vacpla(i1,7).EQ.-2.0) THEN
                  pinasd(i2,1,species,1) = knbs(nks(irtrap+1),irtrap+1)*
     .                                     1.E-6*2.0
                ELSE
                  pinasd(i2,1,species,1) = vacpla(i1,7)*1.E-6
                ENDIF
                IF     (vacpla(i1,6).EQ.-1.0) THEN
                  pinasd(i2,2,species,1) = osmbulkte
                ELSEIF (vacpla(i1,6).EQ.-2.0) THEN
                  pinasd(i2,2,species,1) = ktebs(nks(irtrap+1),irtrap+1)
                ELSE
                  pinasd(i2,2,species,1) = vacpla(i1,6)
                ENDIF
                pinasd(i2,3,species,1) = 0.0
                pinasd(i2,4,species,1) = 0.0
                pinasd(i2,5,species,1) = 0.0

                IF (cut1.EQ.1) 
     .            WRITE(PINOUT,'(A,2I6,F10.2,1P,E10.2)')
     .              'ASSIGNING:',cell,cut1,vacpla(i1,6),vacpla(i1,7)
              ENDDO

            ENDIF
          ENDDO

        ELSEIF (vacpla(i1,1).EQ.-1.0) THEN
c...      Assign vacuum grid plasma everywhere (including pressure gauges):

          i3 = NINT(vacpla(i1,3))+4

          DO i2 = 2, 2+eirnpgdat+asc_ncell*ascncut
           IF (vacpla(i1,5).NE.99.) pinasd(i2,1,i3,1)=vacpla(i1,5)*1.E-6
           IF (vacpla(i1,4).NE.99.) pinasd(i2,2,i3,1)=vacpla(i1,4)
           IF (vacpla(i1,6).NE.99.) pinasd(i2,3,i3,1)=vacpla(i1,6)*1.E-2
           IF (vacpla(i1,7).NE.99.) pinasd(i2,4,i3,1)=vacpla(i1,7)*1.E-2
           IF (vacpla(i1,8).NE.99.) pinasd(i2,5,i3,1)=vacpla(i1,8)*1.E-2
          ENDDO

        ELSEIF (vacpla(i1,1).EQ. 0.0) THEN
c...      Supplimental data line:

        ELSEIF (NINT(vacpla(i1,1)).GE.1        .AND.
     .          NINT(vacpla(i1,2)).GE.1        .AND.
     .          NINT(vacpla(i1,1)).LE.asc_ncell.AND.
     .          NINT(vacpla(i1,2)).LE.asc_ncell.AND.
     .          vacpla(i1,1).LE.vacpla(i1,2)) THEN
c...      Assign plasma to a specific range of cells:

          DO ncut = 1, ascncut
           DO i2 = NINT(vacpla(i1,1))+1+eirnpgdat+(ncut-1)*asc_ncell,
     .             NINT(vacpla(i1,2))+1+eirnpgdat+(ncut-1)*asc_ncell
           i3 = NINT(vacpla(i1,3)) + 4
           IF (vacpla(i1,5).NE.99.) pinasd(i2,1,i3,1)=vacpla(i1,5)*1.E-6
           IF (vacpla(i1,4).NE.99.) pinasd(i2,2,i3,1)=vacpla(i1,4)
           IF (vacpla(i1,6).NE.99.) pinasd(i2,3,i3,1)=vacpla(i1,6)*1.E-2
           IF (vacpla(i1,7).NE.99.) pinasd(i2,4,i3,1)=vacpla(i1,7)*1.E-2
           IF (vacpla(i1,8).NE.99.) pinasd(i2,5,i3,1)=vacpla(i1,8)*1.E-2
           ENDDO
          ENDDO
        ELSE
          CALL ER('AssignAdditionalCellPlasma','Unidentified cell',*99)
        ENDIF
      ENDDO

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine:  InitializeRelaxation
c
c
      SUBROUTINE InitializeRelaxation
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER i1,i2,status

c...DEV1
      IF (citersol.GT.0.AND.
     .    (rel_opt.EQ.1.OR.rel_opt.EQ.2.OR.rel_opt.EQ.3)) THEN
c      IF (citersol.GT.0.AND.(rel_opt.EQ.2.OR.rel_opt.EQ.3)) THEN
        CALL MS('ReadIn','Overwriting NITERSOL in IODIV')

c...NEW
        nitersol = rel_nstep * rel_niter 
c        nitersol = rel_nstep * rel_niter + 1

        WRITE(PINOUT,*) 'NITERSOL=',nitersol

        IF (nitersol.GT.MAXPINITER) THEN
          CALL ER('ReadIn','Relaxation exceeds maximum allowed number'//
     .                     ' of iterations',*99)
        ENDIF
      ENDIF

      IF (rel_ndata.GT.0) THEN
        DO i1 = 0, rel_nstep
          status = 0
          DO i2 = 1, rel_ndata
            IF (NINT(rel_data(i2,1)).EQ.i1.OR.
     .          NINT(rel_data(i2,1)).EQ.-1.OR.
     .          NINT(rel_data(i2,1)).EQ.-2) status = 1
          ENDDO
          IF (status.EQ.0) THEN
            IF     (relmode.EQ.0) THEN
              rel_ndata = rel_ndata + 1
              rel_data(rel_ndata,1) = REAL(i1)
              rel_data(rel_ndata,2) = REAL(rel_niter)
              rel_data(rel_ndata,3) = REAL(rel_opt)
              rel_data(rel_ndata,4) = rel_frac
              rel_data(rel_ndata,5) = REAL(adp_opt)
              rel_data(rel_ndata,6) = REAL(eirtime)
              rel_data(rel_ndata,7) = rel_tol
            ELSEIF (relmode.EQ.1) THEN
              rel_ndata = rel_ndata + 1
              rel_data(rel_ndata,1) = REAL(i1)
              rel_data(rel_ndata,2) = te_mult_o     
              rel_data(rel_ndata,3) = ti_mult_o     
              rel_data(rel_ndata,4) = n_mult_o      
              rel_data(rel_ndata,5) = te_mult_i     
              rel_data(rel_ndata,6) = ti_mult_i     
              rel_data(rel_ndata,7) = n_mult_i      
              rel_data(rel_ndata,8) = tarshift(IKLO)
              rel_data(rel_ndata,9) = tarshift(IKHI)
            ELSEIF (relmode.EQ.2) THEN
              rel_ndata = rel_ndata + 1
              rel_data(rel_ndata,1) = REAL(i1)
              rel_data(rel_ndata,2) = irsep
              rel_data(rel_ndata,3) = nrs
              rel_data(rel_ndata,4) = 3.0
            ELSEIF (relmode.EQ.3) THEN
              rel_ndata = rel_ndata + 1
              rel_data(rel_ndata,1) = REAL(i1)
              rel_data(rel_ndata,2) = te_mult_o     
              rel_data(rel_ndata,3) = ti_mult_o     
              rel_data(rel_ndata,4) = n_mult_o      
              rel_data(rel_ndata,5) = te_mult_i     
              rel_data(rel_ndata,6) = ti_mult_i     
              rel_data(rel_ndata,7) = n_mult_i      
              rel_data(rel_ndata,8) =  0.0
              rel_data(rel_ndata,9) = -1.0
            ELSEIF (relmode.EQ.4.OR.relmode.EQ.5.OR.relmode.EQ.7) THEN
              IF (.NOT.(i1.EQ.0.AND.
     .                  (osm_store.NE.-1.OR.relreset.NE.0))) THEN
                 WRITE(0,*) 'I1=',i1
                 CALL ER('InitializeRelaxation',
     .                   'Plasma over-ride data not found',*99)
              ENDIF
            ELSEIF (relmode.EQ.6) THEN
              DO i2 = 1, eirnpuff
                IF (eirpuff(i2,1).NE.5.0) 
     .            CALL ER('InitializeRelaxation','080 1.6 over-ride '//
     .                    'only works with puff mode 5',*99)
              ENDDO
            ELSE
              CALL ER('InitializeRelaxation','Invalid RELMODE',*99)
            ENDIF
          ENDIF
        ENDDO

        IF (relmode.EQ.3) THEN
          DO i1 = 1, rel_ndata
            IF (rel_data(i1,9).EQ.-1.0) rel_data(i1,9) = REAL(eirtime)
          ENDDO
        ENDIF


c...??? What was this for ??? Aug 9, 2000 -SL
c        DO i1 = 1, rel_ndata
c          IF (rel_data(i1,1).EQ.0.0) rel_data(i1,2) = 1.0
c        ENDDO

        WRITE(PINOUT,*) 'Step parameters:'
        DO i1 = 1, rel_ndata
          WRITE(PINOUT,'(1P,9E10.2)') (rel_data(i1,i2),i2=1,9)
c          IF (outmode.GE.2) WRITE(0,'(7F10.3)') (rel_data(i1,i2),i2=1,7)          
        ENDDO

        IF (relmode.EQ.0) THEN
          nitersol = 0
          DO i1 = 1, rel_ndata                              
            IF (rel_data(i1,1).NE.0.0)
     .        nitersol = nitersol + INT(rel_data(i1,2))
          ENDDO

          nitersol = nitersol + 1

          IF (nitersol.GT.MAXPINITER)
     .      CALL ER('ReadIn','Relaxation exceeds maximum allowed '//
     .                       'number of iterations',*99)
        ENDIF

      ENDIF

c...  Count the loaded solution as the initial DIVIMP
c     solution:
c      IF (osm_store.NE.-1) nitersol = nitersol - 1

c      IF (osm_mode.GE.2) THEN
c        WRITE(0     ,*) 'NITERSOL= ',nitersol
c        WRITE(PINOUT,*) 'NITERSOL= ',nitersol
c      ENDIF


      IF (s21_ndatai.GT.0) THEN
        s21_datai(s21_ndatai+1,2) = terati
        s21_datai(s21_ndatai+1,3) = tirati
        s21_datai(s21_ndatai+1,4) = nrati
        s21_datai(s21_ndatai+1,5) = qrati
        s21_datai(s21_ndatai+1,6) = l1rati
        s21_datai(s21_ndatai+1,7) = l2rati
        s21_datai(s21_ndatai+1,8) = lvrati
      ENDIF

      IF (s21_ndatao.GT.0) THEN
        s21_datao(s21_ndatao+1,2) = terat
        s21_datao(s21_ndatao+1,3) = tirat
        s21_datao(s21_ndatao+1,4) = nrat
        s21_datao(s21_ndatao+1,5) = qrat
        s21_datao(s21_ndatao+1,6) = l1rat
        s21_datao(s21_ndatao+1,7) = l2rat
        s21_datao(s21_ndatao+1,8) = lvrat
      ENDIF

c      DO i1 = 1, s21_ndatai+1
c        WRITE(0,'(A,8F7.2)') 'SOL21 I ',(s21_datai(i1,i2),i2=1,8)
c      ENDDO

c      DO i1 = 1, s21_ndatao+1
c        WRITE(0,'(A,8F7.2)') 'SOL21 O ',(s21_datao(i1,i2),i2=1,8)
c      ENDDO

       IF (osm_matchs.EQ.1.AND.nlpdato2.NE.0)
     .   CALL ER('ReadIn','Can''t relax and float boundary values',*99)

       IF (osm_matchs.EQ.2.AND.nlpdati2.NE.0)
     .   CALL ER('ReadIn','Can''t relax and float boundary values',*99)

c      IF (inputflag(56).NE.1) STOP 'Missing input flag 56'
c      IF (inputflag(57).NE.1) STOP 'Missing input flag 57'

c...  Initialize call to target data relaxation routine:
      CALL UpdateTargets(-1)
      
      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: ReadWallFlux
c
c
c
      SUBROUTINE ReadWallFlux
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER    MAXATM  ,MAXMOL  ,MAXPLS  ,MAXDAT
      PARAMETER (MAXATM=6,MAXMOL=1,MAXPLS=6,MAXDAT=20000)

      REAL       TOL
      PARAMETER (TOL=1.0E-06)

      INTEGER   nseg,ntar,natm,nmol,npls,i1,i2,i3,i4,i5,idum1,idum2,
     .          cvesm(MAXDAT),iw,in,iliin(MAXDAT),ndone,idone(1000)
      LOGICAL   done
      REAL      sdat(MAXDAT,7),len,cir,
     .          adatp(MAXDAT,MAXATM,3),adate(MAXDAT,MAXATM,3),
     .          mdatp(MAXDAT,MAXMOL,3),mdate(MAXDAT,MAXMOL,3),
     .          pdatp(MAXDAT,MAXPLS,3),pdate(MAXDAT,MAXPLS,3),
     .          rdum(6)
      CHARACTER buffer*72

      INTEGER fp

      COMMON /SURFACEMAPCOM/ surfacemap,surfacesrc
      INTEGER surfacemap(5000),surfacesrc(5000)            

      fp = PINOUT


c ALLOCATE PDATP, etc. DYNAMICALLY!


      READ(IPINOUT,'(A72)',ERR=97,END=98) buffer
c      WRITE(0,*) buffer
      IF (buffer(2:6).NE.'[SURF') CALL ER('ReadWallFLux','Data not '//
     .                                    'found',*99)

c...  Read wall flux data and store in local arrays:

      CALL RZero(sdat ,MAXDAT*7)
      CALL RZero(adatp,MAXDAT*MAXATM*3)
      CALL RZero(adate,MAXDAT*MAXATM*3)
      CALL RZero(mdatp,MAXDAT*MAXMOL*3)
      CALL RZero(mdate,MAXDAT*MAXMOL*3)
      CALL RZero(pdatp,MAXDAT*MAXPLS*3)
      CALL RZero(pdate,MAXDAT*MAXPLS*3)
      CALL IZero(cvesm,MAXDAT)

c...  Vessel wall geometry:
      READ(IPINOUT,*)
      READ(IPINOUT,*) nseg

      IF (nseg.GT.MAXDAT) 
     .  CALL ER('ReadWallFlux','Array bounds exceeded, increase '//
     .                         'MAXDAT',*99)

      DO i1 = 1, nseg
        READ(IPINOUT,*) 
     .    idum1,iliin(i1),
     .    sdat(i1,4),sdat(i1,5),sdat(i1,6),
     .    sdat(i1,1),sdat(i1,2),sdat(i1,3)

        IF ((sdat(i1,6).NE.-1.0E+18.OR.sdat(i1,3).NE.+1.0E+18).OR.
     .      (sdat(i1,1).EQ.sdat(i1,4).AND.
     .       sdat(i1,2).EQ.sdat(i1,5)))
     .    sdat(i1,7) = -1.0
      ENDDO

c...  Particle and energy flux of neutral atoms:
      READ(IPINOUT,*)
      READ(IPINOUT,*) natm
      IF (natm.GT.MAXATM) CALL ER('ReadFluxData','MAXATM exceeded',*99)
      DO i1 = 1, nseg
        READ(IPINOUT,*) idum1,idum2,(adatp(i1,i2,1),i2=1,natm),
     .                              (adate(i1,i2,1),i2=1,natm),
     .                              (adatp(i1,i2,2),i2=1,natm),
     .                              (adate(i1,i2,2),i2=1,natm)
      ENDDO

c...  Particle and energy flux of neutral molecules:
      READ(IPINOUT,*)
      READ(IPINOUT,*) nmol
      IF (nmol.GT.MAXMOL) CALL ER('ReadFluxData','MAXMOL exceeded',*99)
      DO i1 = 1, nseg
        READ(IPINOUT,*) idum1,idum2,(mdatp(i1,i2,1),i2=1,nmol),
     .                              (mdate(i1,i2,1),i2=1,nmol),
     .                              (mdatp(i1,i2,2),i2=1,nmol),
     .                              (mdate(i1,i2,2),i2=1,nmol)
      ENDDO

c...  Particle and energy flux of bulk ions:
      READ(IPINOUT,*)
      READ(IPINOUT,*) npls
      IF (npls.GT.MAXPLS) CALL ER('ReadFluxData','MAXPLS exceeded',*99)
      DO i1 = 1, nseg
        READ(IPINOUT,*) idum1,idum2,(pdatp(i1,i2,1),i2=1,npls),
     .                              (pdate(i1,i2,1),i2=1,npls)
      ENDDO

c...  Calculate net flux for atoms and molecules:
      DO i1 = 1, nseg
        DO i2 = 1, natm
          adatp(i1,i2,3) = adatp(i1,i2,1) + adatp(i1,i2,2)
          adate(i1,i2,3) = adate(i1,i2,1) + adate(i1,i2,2)
        ENDDO
        DO i2 = 1, nmol
          mdatp(i1,i2,3) = mdatp(i1,i2,1) + mdatp(i1,i2,2)
          mdate(i1,i2,3) = mdate(i1,i2,1) + mdate(i1,i2,2)
        ENDDO
      ENDDO    

      DO i4 = 1, 2
c...    Dump data to PINOUT for other non-standard (DIVIMP) additional
c       surfaces:
        IF (i4.EQ.1) THEN
          fp = PINOUT
        ELSE
          fp = 98
          OPEN(UNIT=fp,FILE='tmp-001.dat',ACCESS='SEQUENTIAL',
     .         STATUS='UNKNOWN',POSITION='APPEND')
        ENDIF

        WRITE(fp,*)
        WRITE(fp,'(5X,A,I4)') 'STEP:',rel_step
        WRITE(fp,*)
        WRITE(fp,90) 'NON-STANDARD DIVIMP WALL SURFACES (ATOM FLUXES):'
        WRITE(fp,'(5X,12X,2(2A29   ))')
     .    '  ---------- ATOMS ----------',
     .    '  -------- MOLECULES --------'
        WRITE(fp,'(5X,2A6,2(2A10,A9))')
     .    'INDEX','ILIIN','Pflux','Eflux','Eavg','Pflux','Eflux','Eavg'

        DO i3 = 1, 3
          WRITE(fp,*)
          IF (i3.EQ.1) WRITE(fp,90) 'POSITIVE FLUX'
          IF (i3.EQ.2) WRITE(fp,90) 'NEGATIVE FLUX'
          IF (i3.EQ.3) WRITE(fp,90) 'NET      FLUX'
          ndone = 0

          DO i1 = 1, nseg
            IF (iliin(i1).EQ.-2.OR.iliin(i1).EQ.2) THEN

c...          Should really map these surfaces to DIVIMP indicies -- of some
c             sort anyway:

c...          Check if this surface has been done already:
              done = .FALSE.
              DO i2 = 1, ndone              
                IF (idone(i2).EQ.surfacemap(i1)) done = .TRUE.
              ENDDO

              IF (.NOT.done) THEN
c...            Register surface as finished:
                ndone = ndone + 1
                idone(ndone) = surfacemap(i1)

                CALL RZero(rdum,6)
                DO i5 = 1, nseg
c...              Scan through all surfaces and add them up:
                  i2 = 1
                  IF (surfacemap(i1).EQ.surfacemap(i5).AND.
     .                iliin(i1).EQ.iliin(i5)) THEN
                    rdum(1) = rdum(1) + adatp(i5,i2,i3)
                    rdum(2) = rdum(2) + adate(i5,i2,i3)
                    rdum(3) = rdum(3) + adate(i5,i2,i3)
                    rdum(4) = rdum(4) + mdatp(i5,i2,i3) * 2.0
                    rdum(5) = rdum(5) + mdate(i5,i2,i3)
                    rdum(6) = rdum(6) + mdate(i5,i2,i3)
                  ENDIF
                ENDDO

c...            Get averages:
                rdum(3) = rdum(3) / (rdum(1) + 1.0E-10)
c...BUG!  This caused me a lot of gried (June 11, 2003):
c                rdum(6) = rdum(6) / (rdum(4) + 1.0E-10) * 2.0
                rdum(6) = rdum(6) / (rdum(4) + 1.0E-10)

                WRITE(fp,'(5X,2I6,2(1P,2E10.2,0P,F9.3),2X,2I6)')
     .            i1,iliin(i1),(rdum(i5),i5=1,6),surfacesrc(i1),
     .            surfacemap(i1)
              ENDIF


c...          Hydrogen only:
c              i2 = 1
c              WRITE(fp,'(5X,2I6,2(1P,2E10.2,0P,F9.3))')
c     .          i1,iliin(i1),
c     .          adatp(i1,i2,i3),adate(i1,i2,i3),
c     .          adate(i1,i2,i3)/(adatp(i1,i2,i3)+1.0E-10),
c     .          2.0*mdatp(i1,i2,i3),mdate(i1,i2,i3),
c     .          mdate(i1,i2,i3)/(mdatp(i1,i2,i3)+1.0E-10)
            ENDIF
          ENDDO



c...      Old output:
c          WRITE(fp,*)
c          DO i1 = 1, nseg
c            IF (iliin(i1).EQ.-2.OR.iliin(i1).EQ.2) THEN
cc...          Hydrogen only:
c              i2 = 1
c              WRITE(fp,'(5X,2I6,2(1P,2E10.2,0P,F9.3),2X,I6)')
c     .          i1,iliin(i1),
c     .          adatp(i1,i2,i3),adate(i1,i2,i3),
c     .          adate(i1,i2,i3)/(adatp(i1,i2,i3)+1.0E-10),
c     .          2.0*mdatp(i1,i2,i3),mdate(i1,i2,i3),
c     .          mdate(i1,i2,i3)/(mdatp(i1,i2,i3)+1.0E-10),
c     .          surfacemap(i1)
c            ENDIF
c          ENDDO

        ENDDO
      ENDDO


      CLOSE(fp)


      DO i1 = 1, nseg
        IF (iliin(i1).EQ.-2.OR.iliin(i1).EQ.2) THEN
c...      Hydrogen only:
          i2 = 1
          WRITE(79,'(A,6I6,4(1P,2E10.2,0P,F9.3))')
     .      '''PINFLUX   1.01''',
     .      rel_step,rel_iter,rel_count,i1,i2,iliin(i1),
     .      (adatp(i1,i2,i3),adate(i1,i2,i3),
     .       adate(i1,i2,i3)/(adatp(i1,i2,i3)+1.0E-10),
     .       2.0*mdatp(i1,i2,i3),mdate(i1,i2,i3),
     .       mdate(i1,i2,i3)/(mdatp(i1,i2,i3)+1.0E-10),i3=1,2)
        ENDIF
      ENDDO

      WRITE(PINOUT,*)
      DO i1 = nvesm+nvesp+eirnpgdat, nseg
        DO i2 = 1, natm
          IF (adatp(i1,i2,1).NE.0.0.OR.adate(i1,i2,1).NE.0.0) 
     .      WRITE(PINOUT,'(5X,A,3I6,1P,2E12.4,0P)')
     .        'ATM FLUXES: ',rel_count,i1,i2,
     .            adatp(i1,i2,1),adate(i1,i2,1)
        ENDDO
        DO i2 = 1, nmol
          IF (mdatp(i1,i2,1).NE.0.0.OR.mdate(i1,i2,1).NE.0.0) 
     .      WRITE(PINOUT,'(5X,A,3I6,1P,2E12.4,0P)')
     .        'MOL FLUXES: ',rel_count,i1,i2,
     .        2.0*mdatp(i1,i2,1),mdate(i1,i2,1)
        ENDDO
      ENDDO    

c...  Assign FLXHWx arrays:
c
c     FLUXHW - FLUX OF HYDROGEN (ATOMS AND MOLECULES) TO THE WALL
c     FLXHW2 - FLUX OF HYDROGEN (ATOMS AND IONS) TO THE WALL
c     FLXHW3 - FLUX OF IMPURITIES SPUTTERED FROM THE WALL (N/A)
c     FLXHW4 - FLUX OF IMPURITIES REDEPOSITED ONTO THE WALL (N/A)
c     FLXHW5 - AVERAGE ENERGY OF ATOMS HITTING THE WALL (EV)
c     FLXHW6 - FLUX OF HYDROGEN ATOMS TO THE WALL
c     FLXHW7 - AVERAGE ENERGY OF MOLECULES HITTING THE WALL (eV)
c     FLXHW8 - EIRENE REPORTED HYDROGEN ION FLUXES TO THE WALL 
c
      CALL RZero(fluxhw,MAXSEG)
      CALL RZero(flxhw2,MAXSEG)
      CALL RZero(flxhw3,MAXSEG)
      CALL RZero(flxhw4,MAXSEG)
      CALL RZero(flxhw5,MAXSEG)
      CALL RZero(flxhw6,MAXSEG)
      CALL RZero(flxhw7,MAXSEG)
      CALL RZero(flxhw8,MAXSEG)

      DO i1 = 1, nvesm+nvesp
        DO i2 = 1, nseg
          IF (sdat(i2,7).EQ.-1.0) CYCLE

c...      Scan over local wall/target segment data and find corresponding
c         xVESM index by matching segment verticies:
          IF ((ABS(sdat(i2,1)-rvesm(i1,1)).LT.TOL.AND.
     .         ABS(sdat(i2,2)-zvesm(i1,1)).LT.TOL.AND.
     .         ABS(sdat(i2,4)-rvesm(i1,2)).LT.TOL.AND.
     .         ABS(sdat(i2,5)-zvesm(i1,2)).LT.TOL).OR.
     .        (ABS(sdat(i2,1)-rvesm(i1,2)).LT.TOL.AND.
     .         ABS(sdat(i2,2)-zvesm(i1,2)).LT.TOL.AND.
     .         ABS(sdat(i2,4)-rvesm(i1,1)).LT.TOL.AND.
     .         ABS(sdat(i2,5)-zvesm(i1,1)).LT.TOL)) THEN

            cir = 2 * PI * 0.5 * (rvesm(i1,1) + rvesm(i1,2))

            len = SQRT((rvesm(i1,1) - rvesm(i1,2))**2.0 +
     .                 (zvesm(i1,1) - zvesm(i1,2))**2.0)

            fluxhw(i1) = (adatp(i2,1,1) + mdatp(i2,1,1)) / len / cir
            flxhw2(i1) = (adatp(i2,1,1) + pdatp(i2,1,1)) / len / cir
c            flxhw3(i1) = pdatp(i2,1,1)
c            flxhw4(i1) = adatp(i2,1,1)
            flxhw3(i1) = 0.0
            flxhw4(i1) = 0.0
            flxhw5(i1) = adate(i2,1,1) / (adatp(i2,1,1) + 1.0E-10)
            flxhw6(i1) = adatp(i2,1,1) / len / cir
            flxhw7(i1) = mdate(i2,1,1) / (mdatp(i2,1,1) + 1.0E-10)
c
c           jdemod - added flxhw8 - Eirene ion flux to walls
c
c...        If this is fixed on the EIRENE side so that it is no longer
c           statistical, then the use of FLXHW8 in the
c           CALC_TARGFLUXDATA routine can be removed:
            flxhw8(i1) = pdatp(i2,1,1) / len / cir 
c
            cvesm(i1) = 1
          ENDIF
        ENDDO
      ENDDO

c...  Check that flux data was assigned to every wall/target segment:
      DO i1 = 1, nvesm+nvesp
        IF (cvesm(i1).EQ.0) WRITE(0,*) 'WARNING: NO FLUX DATA, I1= ',i1
      ENDDO

c...  Not sure what PINCOR is about at the moment, so issue a warning
c     if it is not unity:
      IF (pincor.NE.1.0) THEN
        WRITE(PINOUT,*)
        WRITE(PINOUT,*) 'WARNING: PINCOR not unity, setting to 1.0'      
        WRITE(PINOUT,*)
        pincor = 1.0
      ENDIF

c...  Not sure why at the moment, but this is done in READPIN and the _PIN
c     arrays are referenced in TAU:
      DO iw = 1, nvesm+nvesp
        fluxhw(iw) = fluxhw(iw) * pincor
        flxhw2(iw) = flxhw2(iw) * pincor
        flxhw3(iw) = flxhw3(iw) * pincor
        flxhw4(iw) = flxhw4(iw) * pincor
        flxhw6(iw) = flxhw6(iw) * pincor

c...    Copy wall flux data:
        fluxhw_pin(iw) = fluxhw(iw)
        flxhw2_pin(iw) = flxhw2(iw)
        flxhw3_pin(iw) = flxhw3(iw)
        flxhw4_pin(iw) = flxhw4(iw)
        flxhw5_pin(iw) = flxhw5(iw)
        flxhw6_pin(iw) = flxhw6(iw)

c...    Copy vessel definitions:
        jvesm_pin(iw)  = jvesm(iw)
        DO in = 1,2
          rvesm_pin(iw,in)  = rvesm(iw,in)
          zvesm_pin(iw,in)  = zvesm(iw,in)
        ENDDO
      ENDDO

c...  From READPIN:
      IF ((WLPABS.EQ.2.OR.WLPABS.EQ.3).AND.CGEOOPT.NE.-1) THEN
        CALL ER('ReadFluxData','Probability data not passed from '//
     .                         'EIRENE',*99)
      ENDIF
c
c     HCORR   - volume correction factor
c     HVAL    - volume correction factor
c
      CALL RSet(hcorr,MAXNKS*MAXNRS,1.0)
      CALL RSet(hval ,MAXNKS*MAXNRS,1.0)

      RETURN
90    FORMAT(5X,A)
97    CALL ER('ReadWallFlux','Problem reading EIRENE transfer file',*99)
98    CALL ER('ReadWallFlux','End of file',*99)
99    STOP
      END
c
c ======================================================================
c
c subroutine: SetupIteration
c
      SUBROUTINE SetupIteration(iitersol)
      IMPLICIT none

      INTEGER region,ir

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      COMMON /OPTTEMP/ osm_matcht,forcet1
      INTEGER          osm_matcht,forcet1

      INTEGER iitersol,i1,i2,sav_niter

      DATA sav_niter /-1/  ! No idea what this is for, copied from UpdateTargets
      SAVE

c
c Target condition relaxation:
c

      rel_count = rel_count + 1

c      WRITE(0,*) 'SETUPITERATION',iitersol

      IF (sav_niter.EQ.-1) sav_niter = rel_niter

      IF (citersol.EQ.0.OR.rel_opt.EQ.0) THEN
        rel_iter  = 0
        rel_step  = 0
      ELSE
        IF     (iitersol.EQ.-1) THEN
          rel_iter  = 0
          rel_step  = 0
        ELSEIF (iitersol.EQ.1) THEN
          rel_iter  = 1
          rel_step  = 0
        ELSEIF (rel_nstep.EQ.1) THEN
          IF (rel_step.EQ.0) rel_iter = 0
          rel_iter = rel_iter + 1
          rel_step = 1
        ELSE
          IF (rel_step.EQ.0.OR.rel_iter.GE.rel_niter) THEN
            rel_viter(rel_step) = rel_iter
            rel_iter            = 1
            rel_step            = rel_step + 1

            WRITE(PINOUT,'(2X,A,3I4)') 'NEW STEP = ',
     .        rel_step,rel_iter,rel_viter(rel_step-1)

          ELSE
            rel_iter = rel_iter + 1
          ENDIF
        ENDIF
      ENDIF


      IF (rel_ndata.GT.0) CALL SetupRelaxation



      IF (rel_opt.EQ.2.OR.rel_opt.EQ.3) THEN

c...    Not sure if I really want this to be the standard, but
c       if the sources are being reset at the beginning of each
c       step, then all steps should have the same number of
c       iterations:
        IF (rel_count.EQ.0.AND.relreset.GE.1) THEN
          osm_mode  = 1
          rel_count = rel_count + 1
          rel_iter  = rel_niter
          iitersol  = iitersol + 1
 
          IF (osm_matcht.EQ.0) CALL SaveSolution
        ENDIF
      ENDIF





      IF (tarninter(IKLO).NE.NULL.OR.tarninter(IKHI).NE.NULL) THEN
        CALL InterpolateTargetData

      ELSEIF ((rel_opt.EQ.2.OR.rel_opt.EQ.3).AND.
     .    tarninter(IKLO).EQ.NULL.AND.tarninter(IKHI).EQ.NULL) THEN
        STOP 'WHAT THE?'
        CALL UpdateTargets(iitersol)   ! *** not sure this still works! ***

      ELSEIF (rel_opt.EQ.2.OR.rel_opt.EQ.3) THEN
        CALL ER('SetupIteration','Need to check if target data '//
     .          'interpolation will work if only 1 target '//
     .          'specified',*99)
      ENDIF





      IF (rel_opt.EQ.0) THEN
        rel_step            = rel_count
        rel_iter            = 1
        rel_viter(rel_step) = rel_iter
      ENDIF



      IF (iflexopt(6).EQ.10.AND.rel_opt.EQ.0) THEN
c...    Read in the background plasma and sources:
        CALL LoadPIN

c...    Adjust REL_COUNT, since it is modified in LOADPIN:
        rel_count = rel_count - 1

      ELSEIF (relreset.GE.1.AND.rel_count.GT.0.AND.
     .        rel_iter.GE.1.AND.rel_opt.GE.1) THEN

        IF (relreset.EQ.1) THEN
c...      Clear n-n collision data on the vacuum grid:
          indasd = 0
          CALL RZero(pinasd,MAXASCDAT*MAXASD2*MAXASS*2)
c...      Set additional cell "plasma" values to EIRENE vacuum 
c         defaults:
          DO i1 = 1, MAXASS
            DO i2 = 1, MAXASCDAT
              pinasd(i2,1,i1,1) = 1.0E+01
              pinasd(i2,2,i1,1) = 1.0E-02
              pinasd(i2,3,i1,1) = 0.0
              pinasd(i2,4,i1,1) = 0.0
              pinasd(i2,5,i1,1) = 0.0
            ENDDO
          ENDDO
c...      Clear n-n collision data on the standard grid:
          CALL RZero(pinbgk,MAXNKS*MAXNRS*MAXBGK*MAXTOR)
        ENDIF
       
c...    Reset PINATOM for opacity calculations:
        tagpinatom = .FALSE.
        CALL RZero(pinatom ,MAXNKS*MAXNRS)

        CALL RZero(pinion,MAXNKS*MAXNRS)
        CALL RZero(pinrec,MAXNKS*MAXNRS)
        CALL RZero(pinmp,MAXNKS*MAXNRS)
        CALL RZero(pinqe,MAXNKS*MAXNRS)
        CALL RZero(pinqi,MAXNKS*MAXNRS)
        CALL RZero(osmion,MAXNKS*MAXNRS)
        CALL RZero(osmrec,MAXNKS*MAXNRS)
        CALL RZero(osmmp,MAXNKS*MAXNRS)
        CALL RZero(osmqe,MAXNKS*MAXNRS)
        CALL RZero(osmqi,MAXNKS*MAXNRS)
        CALL RZero(osmqi,MAXNKS*MAXNRS)
        CALL RZero(osmpei,MAXNKS*MAXNRS)
        CALL RZero(osmpmk(0,1),(MAXNKS+1)*MAXNRS)
        CALL DZero(osmpmk2(0,1),(MAXNKS+1)*MAXNRS)


        IF (osm_store.EQ.-2) THEN
c...      Set REL_STORE to REL_STEP if REL_STORE is equal to -2:
          WRITE(0,*)
          WRITE(0,*) '*** LOADING BACKGROUND PLASMA ***',rel_step
          WRITE(0,*)

c...MAKE SURE THIS ROUTINE IS NOT CALLED MORE THAN ONCE PER STEP...
c...CHECK THAT OSM_STORE LESS THAN ZERO IS OKAY ELSEWHERE...
c...SET OSM_STORE BACK TO 

          osm_store = rel_step

c...      Read in the background plasma and sources:
          CALL LoadPIN

          osm_store = -2

c...      Adjust REL_COUNT, since it is modified in LOADPIN:
          rel_count = rel_count - 1

        ELSEIF (iflexopt(6).EQ.10) THEN
          WRITE(0,*)
          WRITE(0,*) '*** LOADING BACKGROUND PLASMA ***'
          WRITE(0,*)

c...      Read in the background plasma and sources:
          CALL LoadPIN

c...      Adjust REL_COUNT, since it is modified in LOADPIN:
          rel_count = rel_count - 1
        ELSE
c...Hack job here: Blanking PIN sources at the start of each step! - OPTION FIXED
          WRITE(0,*)
          WRITE(0,*) '*** RESETTING PINITER, PINATOM AND BGK DATA '//
     .               'ARRAYS ***'
          WRITE(0,*)

          piniter = .FALSE.
        ENDIF        

      ENDIF

      IF (osm_probe.GE.1) CALL InterpolateProbeData(osm_probe)

      IF (tarsource.NE.0.AND.osm_matcht.NE.0) CALL ProbePath

      WRITE(PINOUT,*)
      IF (rel_opt.EQ.2.OR.rel_opt.EQ.3) THEN
        WRITE(PINOUT,'(4(A,I3),A,L4)')
     .    'ITERATION ',rel_count,'  STEP ',rel_step,
     .    '  SUB-ITER. ',rel_iter,' (',rel_count,')',piniter
      ELSE
        WRITE(PINOUT,'(A,I3)') 'ITERATION ',rel_count
      ENDIF
      WRITE(PINOUT,*)


c...  This call copies the plasma arrays into OLDx2 arrays 
c     in the OLDPLASMA common block associated with the SOL24 
c     code.  There has been some duplication of arrays (OLDx in the 
c     PIN_CFD common block) because the code segments
c     were developed in parallel, and this should be 
c     sorted out at some point:

      CALL MirrorOldPlasma(ktebs,ktibs,knbs,kvhs)



      IF (rel_opt.NE.2.AND.rel_opt.NE.3.AND.    
     .    (tarshift(IKLO).NE.0.0.OR.tarshift(IKHI).NE.0.0)) 
     .  CALL ShiftTargetData



      RETURN
99    STOP
      END
c
c ======================================================================
c
c
c function: GetRelaxationFraction
c
      REAL FUNCTION GetRelaxationFraction()
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slcom' 

      REAL frac

      IF (rel_opt.EQ.2.OR.rel_opt.EQ.3) THEN
        IF     (rel_step.EQ.0.OR.rel_nstep.EQ.1) THEN
          frac = 0.0
        ELSEIF (rel_step.EQ.rel_nstep+1) THEN
          frac = 1.0
        ELSEIF (rel_step.GT.rel_nstep+1) THEN
          CALL ER('GetRelaxationFraction','Invalid REL_STEP value',*99)
        ELSE
          frac = REAL(rel_step - 1) / REAL(rel_nstep - 1)
        ENDIF
        IF (rel_pace.GT.0.0) frac = frac**(1.0 / rel_pace)
        frac = rel_bound1 + frac * (rel_bound2 - rel_bound1) 
c        WRITE(0,*) 'FRAC:',frac,rel_step,rel_nstep
      ELSE
        frac = 0.0
      ENDIF

      GetRelaxationFraction = frac

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c
c subroutine: InterpolateTargetData
c
      SUBROUTINE InterpolateTargetData
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      REAL GetRelaxationFraction

      INTEGER iitersol

      INTEGER ir,ii,ind1,ind2,i1,i2,i3,region,method,ring,target,
     .        no_data_warning
      LOGICAL specific,apply,repeat
      REAL    dum1,dum2,dum3,dum4,exp2,exp3,exp4,dpsin,psin0,frac,
     .        tedat(MAXNRS),tidat(MAXNRS),nedat(MAXNRS)
      DATA no_data_warning / 0 / 
      SAVE

      IF (connected) THEN
        WRITE(0,*) 'CHEATING ON PSITARG'
        DO ir = irsep, nrs
          psitarg(ir,1) = REAL(ir-irsep)/100.0 + 1.0
          psitarg(ir,2) = REAL(ir-irsep)/100.0 + 1.0
        ENDDO
      ENDIF

      CALL OutputData(85,'What the fll')

      WRITE(PINOUT,*)
      WRITE(PINOUT,*) 'INTERPOLATING TARGET DATA'

c...  Need to check that PSITARG is assigned:
c      IF (psitarg(irsep,2).EQ.0.0) 
c     .  CALL ER('InterpolateTargetData','PSITARG does not appear to '//
c     .                                  'be assigned',*99)

      repeat = .TRUE.


      frac = GetRelaxationFraction()

c...  Overrides any explicit assignments in the input file.  Not strictly 
c     required, but leaving this in until it is a problem:
      nlpdato = 0
      nlpdati = 0

c...  Assign target data (LPDATO, LPDATI) from the TARINTER
c     data listed in the DIVIMP/OEDGE input file:
 10     WRITE(PINOUT,*) 'REPEAT',repeat
        DO region = IKLO, IKHI
        ii = 0
        DO ir = irsep, nrs
          ii = ii + 1

          IF (idring(ir).EQ.BOUNDARY) THEN
c...        Assign default values to the virtual rings:
            lpdati(ii,1) = REAL(ir)
            lpdati(ii,2) = 1.0
            lpdati(ii,3) = 1.0
            lpdato(ii,1) = REAL(ir)
            lpdato(ii,2) = 1.0
            lpdato(ii,3) = 1.0
            IF (lpdatsw.EQ.0) THEN
              lpdati(ii,4) = 1.0E+12
              lpdato(ii,4) = 1.0E+12
            ELSE
              lpdati(ii,4) = 1.0E+00
              lpdato(ii,4) = 1.0E+00
            ENDIF
            CYCLE
          ENDIF

c...      Check if target data has been assigned specifically for this ring: 
          specific = .FALSE.
          method = 1
          ind1 = 0
          ind2 = 0

          DO i1 = 1, tarninter(region)
            apply = .FALSE.
            IF (tarinter(i1,1,region).EQ.-1.0) THEN
              specific = .TRUE.

              IF     (ir.GE.NINT(tarinter(i1,2,region)).AND.
     .                ir.LE.NINT(tarinter(i1,3,region))) THEN
c...            Apply target data to ring:
                apply = .TRUE.
              ELSEIF (tarinter(i1,2,region).LT.0.0.AND.
     .                tarinter(i1,3,region).LT.0.0) THEN
c...            Make sure that ring group exists:
                IF (-NINT(tarinter(i1,3,region)).GT.grdntreg(region))
     .            CALL ER('InterpolateTargetData','Specified '//
     .                    'region index is not valid',*99)
c...            Check if ring IR is a member of the specified group(s)
c               of rings:
                DO i2 = -NINT(tarinter(i1,2,region)), 
     .                  -NINT(tarinter(i1,3,region))
                  DO i3 = 1, grdntseg(i2,region)
                    IF (ir.EQ.grdtseg(i3,i2,region)) apply = .TRUE. 

c                    IF (region.EQ.1.AND.ir.EQ.38) THEN
c                      write(0,*) '1:',i1,NINT(tarinter(i1,2,region)),
c     .                                   NINT(tarinter(i1,3,region))
c                      write(0,*) ' :',i2,ir,grdtseg(i3,i2,region),apply
c                    ENDIF

                  ENDDO
                ENDDO                           
              ENDIF
              IF (apply) THEN
c...            Interpolation scheme:
                IF (tarinter(i1,4,region).EQ.1.0) THEN              
c...              Versus PSIn, exponential decay:
                  method = 2
                ELSE
c...              Standard:
                  method = 1
                ENDIF
                ind1 = 0
                ind2 = tarninter(region)
                DO i2 = i1+1, tarninter(region)-1
                  IF (ind1.EQ.0.AND.
     .                tarinter(i2  ,1,region).NE.-1.0) ind1 = i2
                  IF (ind1.NE.0.AND.
     .                tarinter(i2+1,1,region).EQ.-1.0) THEN
                    ind2 = i2 
                    EXIT
                  ENDIF
                ENDDO
                WRITE(PINOUT,*) 'APPLYING:',ir,ind1,ind2,region
                IF (ind2-ind1+1.LT.2) 
     .            CALL ER('InterpolateTargetData','Insufficient '//
     .                    'target data for interpolation',*99)
              ENDIF

            ENDIF
          ENDDO

          IF (.NOT.specific) THEN  
c...        Target data was not specified for specific rings, so treat
c           the entire TARINTER array as target data to be interpolated:
            ind1 = 1
            ind2 = tarninter(region)
          ELSEIF (ind1.EQ.0.OR.ind2.EQ.0) THEN
c...        Data was not found for this ring:
            IF (no_data_warning.EQ.0) no_data_warning = 1
            WRITE(PINOUT,*) 'WARNING: TARGET DATA NOT FOUND FOR',
     .                      ir,region
            CYCLE
          ENDIF

          IF     (method.EQ.1) THEN
c...        Linearly interpolate target data from TARINTER arrays:
            IF (region.EQ.IKLO) THEN
              dum1 = psitarg(ir,2)
             ELSE
              dum1 = psitarg(ir,1)
            ENDIF

c...        Make a list:            
            i3 = 0
            DO i2 = ind1, ind2
              i3 = i3 + 1
              tedat(i3) = (1.0 - frac) * tarinter(i2,2,region) + 
     .                           frac  * tarinter(i2,5,region)
              tidat(i3) = (1.0 - frac) * tarinter(i2,3,region) + 
     .                           frac  * tarinter(i2,6,region)
              nedat(i3) = (1.0 - frac) * tarinter(i2,4,region) + 
     .                           frac  * tarinter(i2,7,region)
            ENDDO
c            IF (ir.EQ.8) THEN
c              WRITE(0,*) 'TE:',tedat(1:2),frac
c              WRITE(0,*) 'Ti:',tidat(1:2),rel_step
c              WRITE(0,*) 'Ne:',nedat(1:2)
c            ENDIF
            CALL Fitter(ind2-ind1+1,tarinter(ind1,1,region),tedat,
     .                  1,dum1,dum2,'LINEAR')
            CALL Fitter(ind2-ind1+1,tarinter(ind1,1,region),tidat,
     .                  1,dum1,dum3,'LINEAR')
            CALL Fitter(ind2-ind1+1,tarinter(ind1,1,region),nedat,
     .                  1,dum1,dum4,'LINEAR')
c            CALL Fitter(ind2-ind1+1,tarinter(ind1,1,region),
c     .                  tarinter(ind1,2,region),1,dum1,dum2,'LINEAR')
c            CALL Fitter(ind2-ind1+1,tarinter(ind1,1,region),
c     .                  tarinter(ind1,3,region),1,dum1,dum3,'LINEAR')
c            CALL Fitter(ind2-ind1+1,tarinter(ind1,1,region),
c     .                  tarinter(ind1,4,region),1,dum1,dum4,'LINEAR')
          ELSEIF (method.EQ.2) THEN
c...        Exponential decay along the target versus PSIn:
            target = 0
            IF (tarinter(ind1,1,region).LT.0.0) THEN
c...          Decay is scaled by the values for a specified target segement:
              ring   = -NINT(tarinter(ind1,1,region))
              target =  NINT(tarinter(ind1,2,region))
              IF     (target.EQ.1) THEN
                DO i2 = 1, nlpdato
                  IF (NINT(lpdato(i2,1)).EQ.ring) THEN
c                  IF (NINT(lpdati(i2,1)).EQ.ring) THEN  ! bug? 22.05.06 -SL
                    tarinter(ind1,2,region) = lpdato(i2,2) / te_mult_o
                    tarinter(ind1,3,region) = lpdato(i2,3) / ti_mult_o
                    tarinter(ind1,4,region) = lpdato(i2,4) /  n_mult_o
                    tarinter(ind1,5,region) = tarinter(ind1,2,region)
                    tarinter(ind1,6,region) = tarinter(ind1,3,region)
                    tarinter(ind1,7,region) = tarinter(ind1,4,region)
                  ENDIF
                ENDDO
                psin0 = psitarg(ring,2)
              ELSEIF (target.EQ.2) THEN
                DO i2 = 1, nlpdati
                  IF (NINT(lpdati(i2,1)).EQ.ring) THEN
                    tarinter(ind1,2,region) = lpdati(i2,2) / te_mult_i
                    tarinter(ind1,3,region) = lpdati(i2,3) / ti_mult_i
                    tarinter(ind1,4,region) = lpdati(i2,4) /  n_mult_i
                    tarinter(ind1,5,region) = tarinter(ind1,2,region)
                    tarinter(ind1,6,region) = tarinter(ind1,3,region)
                    tarinter(ind1,7,region) = tarinter(ind1,4,region)
                  ENDIF
                ENDDO
                psin0 = psitarg(ring,1)
              ELSE
                CALL ER('InterpolateTargetData','Invalid target',*99)
              ENDIF
            ELSE
              psin0 = tarinter(ind1,1,region)
            ENDIF
            IF (region.EQ.IKLO) THEN
              dpsin = ABS(psin0 - psitarg(ir,2))
            ELSE
              dpsin = ABS(psin0 - psitarg(ir,1))
            ENDIF
c...        Interpolate target profile parameters (in case target boundary
c           relaxation is being used):
            DO i2 = 1, 3
              tedat(i2) = (1.0 - frac) * tarinter(ind1+i2-1,2,region) +
     .                           frac  * tarinter(ind1+i2-1,5,region)
              tidat(i2) = (1.0 - frac) * tarinter(ind1+i2-1,3,region) +
     .                           frac  * tarinter(ind1+i2-1,6,region)
              nedat(i2) = (1.0 - frac) * tarinter(ind1+i2-1,4,region) +
     .                           frac  * tarinter(ind1+i2-1,7,region)
            ENDDO
c...        Allow over-ride (RELMODE=8) of profile parameters at each step,
c           when boundary condition relaxation is being used, REL_OPT=2,3: 
            IF (tedat(1).EQ.-1.0) THEN
              IF (relmode.EQ.8.AND.relexpdat(region,1,ir).NE.0.0) THEN
                tedat(1) = relexpdat(region,1,ir)
              ELSEIF (rel_step.EQ.0) THEN
c...            The 0th EIRENE iteration is (currently) being avoided when target relaxation
c               is being used, so just assign some arbitrary nedat(1)!=-1.0.  The value from 
c               RELEXPDAT will be assigned on a subsequent call to this routine before the 
c               plasma is calculated (when REL_STEP=1, see in BgPlasma routine):
                tedat(1) = 1.0E+01
              ELSE
                CALL ER('InterpolateTargetData','Bad Te over-ride',*99)
              ENDIF
            ENDIF
            IF (tidat(1).EQ.-1.0) THEN
              IF (relmode.EQ.8.AND.relexpdat(region,2,ir).NE.0.0) THEN
                tidat(1) = relexpdat(region,2,ir)
              ELSEIF (rel_step.EQ.0) THEN
                tidat(1) = 1.0E+01
              ELSE
                CALL ER('InterpolateTargetData','Bad Ti over-ride',*99)
              ENDIF
            ENDIF
            IF (nedat(1).EQ.-1.0) THEN
              IF (relmode.EQ.8.AND.relexpdat(region,3,ir).NE.0.0) THEN
                nedat(1) = relexpdat(region,3,ir)
              ELSEIF (rel_step.EQ.0) THEN
                nedat(1) = 1.0E+01
              ELSE
                CALL ER('InterpolateTargetData','Bad ne over-ride',*99)
              ENDIF
            ENDIF

c...        Set exponents:
            exp2 = EXP(-dpsin / tedat(3))
            exp3 = EXP(-dpsin / tidat(3))
            exp4 = EXP(-dpsin / nedat(3))
c...        Calculate target boundary conditions:
            dum2 = (tedat(1) - tedat(2)) * exp2 + tedat(2)
            dum3 = (tidat(1) - tidat(2)) * exp3 + tidat(2)
            dum4 = (nedat(1) - nedat(2)) * exp4 + nedat(2)
c...        (Need to replace the target data specifier which was overwritten 
c           above):
            IF (target.NE.0) tarinter(ind1,2,region) = REAL(target)
          ENDIF
c...      Assign LPDATx arrays (target data used to assign target data
c         KxDS arrays):
          IF (region.EQ.IKLO) THEN
            nlpdato = MAX(nlpdato,ii)
            lpdato(ii,1) = REAL(ir)
            lpdato(ii,2) = dum2 * te_mult_o
            lpdato(ii,3) = dum3 * ti_mult_o
            lpdato(ii,4) = dum4 *  n_mult_o
          ELSE
            nlpdati = MAX(nlpdati,ii)
            lpdati(ii,1) = REAL(ir)
            lpdati(ii,2) = dum2 * te_mult_i
            lpdati(ii,3) = dum3 * ti_mult_i
            lpdati(ii,4) = dum4 *  n_mult_i
          ENDIF
c...      End of IR loop:
        ENDDO
c...    End of REGION loop:
      ENDDO

      IF (repeat) THEN
c...    Need to do this twice:
        repeat = .FALSE.
        GOTO 10
      ENDIF

c      CALL Outputdata(85,'sdfds')
c      STOP 'sdgsdg'
	
      IF (no_data_warning.EQ.1) THEN
c        WRITE(0,*)
c        WRITE(0,*) '****************************************'
c        WRITE(0,*) '* TARGET DATA NOT FOUND FOR SOME RINGS *'
c        WRITE(0,*) '****************************************'
c        WRITE(0,*)
        CALL WN('InterpolateTargetData','Data not found for some rings')
        no_data_warning = 2
      ENDIF

      RETURN
99    CONTINUE
      STOP
      END



c
c ======================================================================
c
c subroutine: ApplyProbeData
c
c...THIS ROUTINE IS TO BE RETIRED.
c
c Need to check LPDATSW, and adjust density accordingly...
c
c
      SUBROUTINE ApplyProbeData

      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER CalcPoint

      INTEGER ir,ii,i2
      REAL    dum1,dum2,dum3,dum4,cs

      INTEGER result,i1,idum1,id
      REAL*8  r1,z1,r2,z2,r3,z3,t


      WRITE(SLOUT,*)
      WRITE(SLOUT,*) 'APPLYING PROBE DATA'

      WRITE(0,*) 'WARNING: IT IS ASSUMED THAT THE DENSITY '//
     .           'DATA IN THE PROBE DATA FILE IS n_infinity'

      ii = 0
      DO ir = irsep, nrs
        dum1 = rho(ir,CELL1)

        ii = ii + 1

        lpdati(ii,1) = REAL(ir)
        lpdato(ii,1) = REAL(ir)

        IF (ir.EQ.irwall.OR.ir.EQ.irtrap) THEN
          lpdati(ii,2) = 1.0
          lpdati(ii,3) = 1.0
          lpdati(ii,4) = 1.0
          lpdato(ii,2) = 1.0
          lpdato(ii,3) = 1.0
          lpdato(ii,4) = 1.0
        ELSE

c...      High index target:
          IF (idring(ir).EQ.TARTOTAR.OR.
     .        idring(ir).EQ.WALTOTAR.OR.
     .        stopopt.EQ.121) THEN
            CALL Fitter(prb_num(OFMP),prb_rho(1,OFMP),prb_te(1,OFMP),
     .                  1,dum1,dum2,'LINEAR')
            CALL Fitter(prb_num(OFMP),prb_rho(1,OFMP),prb_ti(1,OFMP),
     .                  1,dum1,dum3,'LINEAR')
            CALL Fitter(prb_num(OFMP),prb_rho(1,OFMP),prb_ne(1,OFMP),
     .                  1,dum1,dum4,'LINEAR')
          ELSE
            CALL ER('ApplyProbeData','Development needed A',*99)
          ENDIF
	
          lpdati(ii,2) = dum2
          lpdati(ii,3) = dum3
c...      The 0.5 converts from n_infinity to n_target:
          lpdati(ii,4) = dum4 * 0.5
	
c...      Low index target:
          IF (((tarsource.EQ.3.OR.tarsource.EQ.7).AND.
     .         dum1.GE.prb_rho(1            ,IFMP).AND.
     .         dum1.LE.prb_rho(prb_num(IFMP),IFMP)).OR.
     .        tarsource.EQ.4.OR.tarsource.EQ.8) THEN
	
            IF (idring(ir).EQ.TARTOTAR.OR.
     .          stopopt.EQ.121) THEN
             CALL Fitter(prb_num(IFMP),prb_rho(1,IFMP),prb_te(1,IFMP),
     .                   1,dum1,dum2,'LINEAR')
             CALL Fitter(prb_num(IFMP),prb_rho(1,IFMP),prb_ti(1,IFMP),
     .                   1,dum1,dum3,'LINEAR')
             CALL Fitter(prb_num(IFMP),prb_rho(1,IFMP),prb_ne(1,IFMP),
     .                   1,dum1,dum4,'LINEAR')
            ELSEIF (idring(ir).NE.WALTOTAR) THEN
              CALL ER('ApplyProbeData','Development needed B',*99)
            ENDIF
          ENDIF
	
          lpdato(ii,2) = dum2
          lpdato(ii,3) = dum3
c...      The 0.5 converts from n_infinity to n_target:
          lpdato(ii,4) = dum4 * 0.5

        ENDIF

      ENDDO




      nlpdati = ii
      nlpdato = ii
c
c     Convert density data to ion saturation current data
c     depending on LPDATSW:
c
      WRITE(SLOUT,*)
      WRITE(SLOUT,'(A)') 'Outer (CMOD) target data'
      DO ii = 1, nlpdati
        IF (lpdatsw.EQ.1) THEN
          cs = SQRT(0.5 * 2.0 * lpdati(ii,2) * (1.0 + RIZB) / CRMB)
          lpdati(ii,4) = lpdati(ii,4) * ECH * 9.79E+03 * cs
        ENDIF

        WRITE(SLOUT,90) (lpdati(ii,i2),i2=1,4)
90      FORMAT(F10.1,2X,2F10.2,2X,1P,E10.2)

        lpdati(ii,2) = lpdati(ii,2) * te_mult_i
        lpdati(ii,3) = lpdati(ii,3) * ti_mult_i
        lpdati(ii,4) = lpdati(ii,4) *  n_mult_i
      ENDDO

      WRITE(SLOUT,*)
      WRITE(SLOUT,'(A)') 'Inner (CMOD) target data'
      DO ii = 1, nlpdato
        IF (lpdatsw.EQ.1) THEN
          cs = SQRT(0.5 * 2.0 * lpdato(ii,2) * (1.0 + RIZB) / CRMB)
          lpdato(ii,4) = lpdato(ii,4) * ECH * 9.79E+03 * cs
        ENDIF

        WRITE(SLOUT,90) (lpdato(ii,i2),i2=1,4)

        lpdato(ii,2) = lpdato(ii,2) * te_mult_o
        lpdato(ii,3) = lpdato(ii,3) * ti_mult_o
        lpdato(ii,4) = lpdato(ii,4) *  n_mult_o
      ENDDO

      RETURN
c
c     Error code:
c
99    CONTINUE
      STOP
      END
c
c ======================================================================
c


c
c =========
c
c subroutine: SetupGrid
c
c
      SUBROUTINE SetupGrid

      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      INTEGER GetModel
      LOGICAL OutsideBreak

      INTEGER ik,ir,ir1,ir2,iki,iko,id1,id2,id3,ii,id,in,midnks,i1,
     .        ikto3,ikti3,ik1,ik2,ik3,ir3,count
      LOGICAL recalculate,cheat1,status
      REAL rhozero,ikintersec(MAXNRS,2)

      REAL*8 a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd
c 

c     jdemod
      real*8 rtemp,ztemp
c     jdemod
c
      

      REAL       TOL
      PARAMETER (TOL=1.0E-06)

      DATA cheat1 /.TRUE./
      SAVE


      CALL SetBounds

c...should be elsewhere
      IF (cioptg.EQ.4.AND.nbgplas.NE.0)
     .  CALL WN('SetupGrid','CIOPTG=4 but NBGPLAS not equal to 0')

      DO ir = irsep, nrs
        osm_model(IKLO,ir) = GetModel(IKLO,ir)
        osm_model(IKHI,ir) = GetModel(IKHI,ir)

        WRITE(PINOUT,*) 'MODEL : ',osm_model(IKLO,ir),osm_model(IKHI,ir)
      ENDDO



c
c     Initialization:
c

      CALL IZero(ikto2,MAXNRS)
      CALL IZero(ikti2,MAXNRS)
      CALL IZero(virloc,MAXNRS*2)


c...  RINGTYPE:
      DO ir = 2, irsep-1
        ringtype(ir) = CORE
      ENDDO
      DO ir = irsep, irwall-1
        ringtype(ir) = SOL1
      ENDDO
      DO ir = irtrap+1, nrs
        ringtype(ir) = PFZ
      ENDDO
c...  Check for secondary PFZ by scanning inward along the connection
c     map to see if the core is encoutered -- if not, then set as a 
c     PFZ ring:
      IF (cgridopt.NE.LINEAR_GRID.AND.ikouts(1,irsep).NE.0) THEN  ! Check (lame) if the connection map is defined
        DO ir = irsep, irwall-1
          status = .TRUE.
          DO ik = 1, nks(ir)
            ik1 = ik
            ir1 = ir
            count = 0
            DO WHILE (idring(ir1).NE.BOUNDARY.AND.count.LE.nrs)
              ik2 = ikins(ik1,ir1)
              ir2 = irins(ik1,ir1)
              WRITE(88,'(A,6I6,L2)') ' PFZ-:',ik,ir,ik2,ir2,
     .                               nks(ir),irsep,status
              IF (ir1.LE.irsep) THEN               ! Changed 02/09/2010
c              IF (ir1.LT.irsep) THEN              ! This scheme is poor!
                WRITE(88,*) 'BOUNCE!'              ! Need to use the OSM method of line intersections!
                status = .FALSE.
                EXIT
              ENDIF
              ik1 = ik2
              ir1 = ir2
              count = count + 1
            ENDDO
            IF (count.EQ.nrs+1) THEN
              CALL WN('SetupGrid','Problem with connection '//
     .                'map when searching for private flux regions')
              WRITE(0,*) '  IK,IR= ',ik,ir
              EXIT
            ENDIF
            IF (.NOT.status) EXIT
          ENDDO
          IF (status) ringtype(ir) = PFZ
        ENDDO
      ENDIF


      IF (grdnmod.EQ.0) THEN
c
c       IDRING:
c

c       REMOVE FROM MAIN VERSION.  THIS IS NOW EXCLUSIVELY ASSIGNED
c       IN THE GRID TAILORING ROUTINES.

        DO ir = 1, nrs
          IF (nbr.GT.0) THEN
c PFZ BREAK
c...grd
            IF (OutsideBreak(ir)) THEN
              IF (stopopt.EQ.121) THEN
                idring(ir) = TARTOWAL
              ELSE
                idring(ir) = WALTOTAR
              ENDIF
            ELSE                                   
              idring(ir) = TARTOTAR                  
            ENDIF                                   
          ELSE                      
            idring(ir) = TARTOTAR        
          ENDIF
        ENDDO                          
                           
      ENDIF

      idring(1)      = BOUNDARY
      idring(irwall) = BOUNDARY
      idring(irtrap) = BOUNDARY  
c
c     VIRLOC:
c
      DO ir = irsep, nrs
        virloc(ir,IKHI) = nks(ir)

        DO ik = 1, nks(ir)
          IF     (virloc(ir,IKLO).EQ.0.AND.virtag(ik,ir).EQ.0) THEN
            virloc(ir,IKLO) = ik
          ELSEIF (virloc(ir,IKLO).NE.0.AND.virtag(ik,ir).EQ.1.AND.
     .            virloc(ir,IKHI).EQ.nks(ir)) THEN
            virloc(ir,IKHI) = ik - 1
          ENDIF
        ENDDO
      ENDDO

c...  Tighten up the grid :
c..THIS MAY MAKE THE CODE BELOW REDUNDANT

c
c     IKTO2 and IKTI2:
c
c...  Check first that IKTO and IKTI are correct:
      
      IF     (connected) THEN
        rxp = rvertp(4,korpg(ikto,irsep)) 
        zxp = zvertp(4,korpg(ikto,irsep)) 
      ELSEIF (nopriv) THEN
        rxp = rvertp(1,korpg(1,irsep)) 
        zxp = zvertp(1,korpg(1,irsep)) 

        ikto = 0
        ikti = nks(irsep) + 1
        ikto2 = ikto
        ikti2 = ikti
        GOTO 50  ! *** SPAGETTI FOR NOW ***
      ELSE
        id1 = korpg(ikto,irsep)      
        id2 = korpg(ikti,irsep)
        IF (ABS(rvertp(4,id1)-rvertp(1,id2)).GT.TOL.OR.
     .      ABS(zvertp(4,id1)-zvertp(1,id2)).GT.TOL) THEN
           IF (sloutput) 
     .       WRITE(pinout,*) 'WARNING: HEALING IKTO AND IKTI'
c...      Find IKTO and IKTI from the separatrix ring:
          ir = irsep
          DO ik1 = 1, nks(ir)-2
            id1 = korpg(ik1,ir)
            DO ik2 = ik1+2, nks(ir)
              id2 = korpg(ik2,ir)
              IF (ABS(rvertp(4,id1)-rvertp(1,id2)).LT.TOL.AND.
     .            ABS(zvertp(4,id1)-zvertp(1,id2)).LT.TOL) THEN
                ikto = ik1
                ikti = ik2
              ENDIF
            ENDDO
          ENDDO        
c         CALL ER('GridSpace','IKTO and IKTI not consistent',*99)
        ENDIF
        rxp = rvertp(4,korpg(ikto,irsep)) 
        zxp = zvertp(4,korpg(ikto,irsep)) 
      ENDIF

      ikto2(irsep) = ikto
      ikti2(irsep) = ikti

c ...unnecessary...
      DO ir = irsep+1, irwall-1
        ikto2(ir) = -1
        ikti2(ir) = -1
      ENDDO
      midnks = 0
      DO ir = irsep, irwall-1
        IF (ikto2(ir).NE.-1.AND.ikti2(ir).NE.-1)
     .    midnks = MAX(ikti2(ir) - ikto2(ir) - 1,midnks)
c Bug:
c     .    midnks = MAX(ikti2(ir) - ikto2(ir) + 1,midnks)
      ENDDO


      midnks = ikti - ikto - 1
      DO ir = irsep+1, irwall-1
c        ikto3 = (virloc(ir,IKHI) - midnks) / 2 + virloc(ir,IKLO)
        IF (connected.AND.ir.EQ.irsep2) THEN
c          WRITE(0,*) irsep,irsep2,irwall,connected
c          STOP 'sdfsdf'      
          id1 = korpg(ikto,irsep)
          id2 = korpg(ikti,irsep)
          DO ik = 1, nks(ir)
            id = korpg(ik,ir)
            IF (ABS(rvertp(1,id)-rvertp(4,id1)).LT.TOL.AND.
     .          ABS(zvertp(1,id)-zvertp(4,id1)).LT.TOL) ikti2(ir) = ik
            IF (ABS(rvertp(4,id)-rvertp(1,id2)).LT.TOL.AND.
     .          ABS(zvertp(4,id)-zvertp(1,id2)).LT.TOL) ikto2(ir) = ik
          ENDDO
          IF (ikto2(ir).EQ.-1.OR.ikti2(ir).EQ.-1) 
     .      CALL ER('SetupGrid','IKTI,O2 problem for IRSEP2',*99)
        ELSE
          ikto3 = virloc(ir,IKLO) - 1 +
     .            (virloc(ir,IKHI) - virloc(ir,IKLO) + 1 - midnks) / 2
          ikti3 = ikto3 + midnks + 1
          ikto2(ir) = MAX(virloc(ir,IKLO),ikto3)
          ikti2(ir) = MIN(virloc(ir,IKHI),ikti3)
        ENDIF  
      ENDDO

      DO ir = nrs, irtrap+1, -1
        IF (idring(ir).EQ.-99) CYCLE
        WRITE(SLOUT,*) 'CUTPOINTS IN PFZ: IR= ',ir
     
        ikto2(ir) = -1
        ikti2(ir) = -1

        IF (ir.EQ.nrs) THEN
          IF (connected) THEN
            ir1 = irsep
            ir2 = irsep2
          ELSE
            ir1 = irsep
            ir2 = ir1
          ENDIF
        ELSE
          ir1 = ir + 1
          ir2 = ir1
        ENDIF

        iko = ikto2(ir1)
        iki = ikti2(ir2)
        id1 = korpg(MAX(1,iko),ir1)
        id2 = korpg(MAX(1,iki),ir2)

        write (slout,'(a,10i6,6(1x,g12.5))') 'PFZ CUTPOINTS:',
     >       ir,ir1,ir2,
     >       ikto2(ir),ikti2(ir),iko,iki,id1,id2

        DO ik = 1, nks(ir)
          id3 = korpg(ik,ir)
          IF (ABS(rvertp(4,id1)-rvertp(3,id3)).LT.TOL.AND.iko.NE.-1.AND.
     .        ABS(zvertp(4,id1)-zvertp(3,id3)).LT.TOL) ikto2(ir) = ik
          IF (ABS(rvertp(1,id2)-rvertp(2,id3)).LT.TOL.AND.iki.NE.-1.AND.
     .        ABS(zvertp(1,id2)-zvertp(2,id3)).LT.TOL) ikti2(ir) = ik

           write (slout,'(a,5i6,6(1x,g12.5))') '   RING:',
     >       ir,ik,id1,id2,id3,
     >       ikto2(ir),ikti2(ir),
     >       rvertp(4,id1)-rvertp(3,id3),zvertp(4,id1)-zvertp(3,id3),
     >       rvertp(1,id2)-rvertp(2,id3),zvertp(1,id2)-zvertp(2,id3)

        ENDDO

        IF (ikti2(ir).EQ.-1.OR.ikto2(ir).EQ.-1)
     .    CALL WN('SetupGrid','Cannot find cut points in PFZ')
c     .    CALL ER('SetupGrid','Cannot find cut points in PFZ',*99)
      ENDDO
c
c     Set IKTO,I2 for double-null grids:
c     ------------------------------------------------------------------    
      DO ir = irsep+1, irwall-1          ! *** REPLACE CHECK WITH FUNCTION OR DOUBLE_NULL FLAG ***
        IF (ringtype(ir).EQ.PFZ) EXIT
      ENDDO
      IF     (ir.NE.irwall.AND.irsep.EQ.irsep2) THEN
        WRITE(0,*) 
        WRITE(0,*) 'DOUBLE NULL GRID DETECTED, BUT IRSEP=IRSEP2.' ! This is showing up on the big JET grid for j-bfg-0004c.
        WRITE(0,*) 'NOT EXECUTING DOUBLE-NULL CODE.'              ! It may just be a 'transient' problem, that goes away once the
        WRITE(0,*)                                                ! rings are properly ordered, but clearly needs some sorting out.
      ELSEIF (ir.NE.irwall) THEN                                  ! (A proper double-null grid detector...) -SL, 09/04/2010
        ir1 = irouts(1,irsep2)
        DO ik1 = 1, nks(ir1)
          IF (irins(ik1,ir1).NE.irsep2) EXIT
        ENDDO
        IF (ik1.EQ.nks(ir1)+1) THEN
          CALL DumpGrid('Identify problems')
          STOP 'DAMNA'
        ELSE
c...      Outer SOL - IKTO:
          ik = ikouts(ikto2(irsep2),irsep2)
          ir = ir1
          DO WHILE (idring(ir).NE.BOUNDARY)
            ikto2(ir) = ik
            ik3 = ik
            ir3 = ir
            ik = ikouts(ik3,ir3)
            ir = irouts(ik3,ir3)
          ENDDO
c...      Outer SOL - IKTI:
          ik = ik1
          ir = ir1
          DO WHILE (idring(ir).NE.BOUNDARY)
            ikti2(ir) = ik
            ik3 = ik
            ir3 = ir
            ik = ikouts(ik3,ir3)
            ir = irouts(ik3,ir3)
          ENDDO
c...      Secondary PFZ:
          ik = ikins(ik1,ir1)
          ir = irins(ik1,ir1)
          DO WHILE (idring(ir).NE.BOUNDARY)
            ikti2(ir) = ik
            ik3 = ik
            ir3 = ir
            ik = ikins(ik3,ir3)
            ir = irins(ik3,ir3)
          ENDDO
        ENDIF
c
        ir1 = irouts(nks(irsep2),irsep2)
        DO ik1 = nks(ir1), 1, -1
          IF (irins(ik1,ir1).NE.irsep2) EXIT
        ENDDO
        IF (ik1.EQ.0) THEN
          STOP 'DAMNB'
        ELSE
c...      Outer SOL - IKTI:
          ik = ikouts(ikti2(irsep2),irsep2)
          ir = ir1
c            write(0,*) 'go man go - start',ik,ir  
          DO WHILE (idring(ir).NE.BOUNDARY)
            ikti2(ir) = ik
c            write(0,*) 'go man go',ir,ikti2(ir)
            ik3 = ik
            ir3 = ir
            ik = ikouts(ik3,ir3)
            ir = irouts(ik3,ir3)
          ENDDO
c...      Outer SOL - IKTI:
          ik = ik1
          ir = ir1
          DO WHILE (idring(ir).NE.BOUNDARY)
            ikto2(ir) = ik
            ik3 = ik
            ir3 = ir
            ik = ikouts(ik3,ir3)
            ir = irouts(ik3,ir3)
          ENDDO
c...      Secondary PFZ:
          ik = ikins(ik1,ir1)
          ir = irins(ik1,ir1)
          DO WHILE (idring(ir).NE.BOUNDARY)
            ikto2(ir) = ik
            ik3 = ik
            ir3 = ir
            ik = ikins(ik3,ir3)
            ir = irins(ik3,ir3)
          ENDDO
        ENDIF

      ENDIF

 50   CONTINUE ! For NOPRIV...
c
c
c     IKMIDPLANE:
c
c
      ikmidplane = 0
      DO ir = irsep, irwall-1
        DO ik = 1, nks(ir)
          id = korpg(ik,ir)
          c1 = 0.5D0 * DBLE(rvertp(1,id) + rvertp(2,id))
          d1 = 0.5D0 * DBLE(rvertp(3,id) + rvertp(4,id))
          c2 = 0.5D0 * DBLE(zvertp(1,id) + zvertp(2,id))
          d2 = 0.5D0 * DBLE(zvertp(3,id) + zvertp(4,id))

          a1 = DBLE(r0)
          b1 = DBLE(r0) - 100.0D0
          a2 = DBLE(z0)
          b2 = DBLE(z0)
c          a2 = 0.0D0
c          b2 = 0.0D0
          CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
          IF (tab.GE.0.0.AND.tab.LE.1.0.AND.
     .        tcd.GE.0.0.AND.tcd.LE.1.0) THEN
            ikmidplane(ir,IKLO) = ik
            ikintersec(ir,IKLO) = r0 - SNGL(tab) * 100.0
          ENDIF

          a1 = DBLE(r0)
          b1 = DBLE(r0) + 100.0D0
          a2 = DBLE(z0)
          b2 = DBLE(z0)
c          a2 = 0.0D0
c          b2 = 0.0D0
          CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
          IF (tab.GE.0.0.AND.tab.LE.1.0.AND.
     .        tcd.GE.0.0.AND.tcd.LE.1.0) THEN
            ikmidplane(ir,IKHI) = ik
            ikintersec(ir,IKHI) = r0 + SNGL(tab) * 100.0
          ENDIF
        ENDDO
      ENDDO
c
c     RHO:
c
c      NEEDS REWORKING!  BIG TIME!
      rho = 0.0
      WRITE(SLOUT,*) 'CALCULATING RHO', r0,z0

      a1 = r0
      a2 = z0   ! This was not a good idea... removed - SL, 19/10/2009  yes it was... 10/12/2009
c      a2 = 0.0D0           
      b1 = r0 + 100.0D0
      b2 = z0
c      b2 = 0.0D0

      DO ir = 2, irwall-1
        DO ik = 1, nks(ir)
          id = korpg(ik,ir)

          c1 = rvertp(1,id)
          c2 = zvertp(1,id)
          d1 = rvertp(4,id)
          d2 = zvertp(4,id)
          CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
          IF (tab.GE.0.0D0.AND.tab.LE.1.0D0.AND.
     .        tcd.GE.0.0D0.AND.tcd.LE.1.0D0)
     .      rho(ir,IN14) = r0 + SNGL(tab) * 100.0D0

          c1 = 0.5 * (rvertp(1,id) + rvertp(2,id))
          c2 = 0.5 * (zvertp(1,id) + zvertp(2,id))
          d1 = 0.5 * (rvertp(3,id) + rvertp(4,id))
          d2 = 0.5 * (zvertp(3,id) + zvertp(4,id))
          CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
          IF (tab.GE.0.0D0.AND.tab.LE.1.0D0.AND.
     .        tcd.GE.0.0D0.AND.tcd.LE.1.0D0)
     .      rho(ir,CELL1) = r0 + SNGL(tab) * 100.0

          c1 = rvertp(2,id)
          c2 = zvertp(2,id)
          d1 = rvertp(3,id)
          d2 = zvertp(3,id)
          CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
          IF (tab.GE.0.0D0.AND.tab.LE.1.0D0.AND.
     .        tcd.GE.0.0D0.AND.tcd.LE.1.0D0)
     .      rho(ir,OUT23) = r0 + SNGL(tab)* 100.0
        ENDDO
      ENDDO

      IF (connected) THEN
        rhozero = rho(irsep2,IN14) 
      ELSE
        rhozero = rho(irsep ,IN14)
      ENDIF

      DO ir = 2, irwall-1
        IF (rho(ir,CELL1).NE.0.0) THEN
          rho(ir,IN14 ) = rho(ir,IN14)  - rhozero
          rho(ir,CELL1) = rho(ir,CELL1)  - rhozero
          rho(ir,OUT23) = rho(ir,OUT23) - rhozero
        ENDIF
      ENDDO

c This sucks... it is good to get rho from grid.. what if not a CMOD
c grid being used...?
      rho(nrs,IN14 ) = 2 * rho(irsep,IN14) - rho(irsep,OUT23)
      rho(nrs,OUT23) = rho(irsep,IN14)
      rho(nrs,CELL1) = 0.5 * (rho(nrs,IN14) + rho(nrs,OUT23))

      DO ir = nrs-1, irtrap+1, -1
        rho(ir,OUT23) = rho(ir+1,IN14)
        rho(ir,IN14 ) = rho(ir+1,IN14) - 
     .                  (rho(ir+1,OUT23) - rho(ir+1,IN14))
c        rho(ir,IN14)  = 2 * rho(ir+1,IN14) - rho(ir,OUT23)
        rho(ir,CELL1)  = 0.5 * (rho(ir,IN14) + rho(ir,OUT23))
      ENDDO
c...  Find inner midplane "rho", for rings that do not intesect
c     the outer midplane:
      a1 = DBLE(r0)
      b1 = DBLE(r0) - 100.0D0
      a2 = 0.0D0
      b2 = 0.0D0
c...  Find RHOZERO:
      IF (connected) THEN 
        ir = irsep2
      ELSE
        ir = irsep
      ENDIF
      rhozero = 0.0
      DO ik = 1, nks(ir)
        id = korpg(ik,ir)
        c1 = DBLE(rvertp(1,id))
        c2 = DBLE(zvertp(1,id))
        d1 = DBLE(rvertp(4,id))
        d2 = DBLE(zvertp(4,id))
        CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
        IF (tab.GE.0.0D0.AND.tab.LE.1.0D0.AND.
     .      tcd.GE.0.0D0.AND.tcd.LE.1.0D0)
     .    rhozero = r0 - SNGL(tab) * 100.0 
      ENDDO

      DO ir = irsep, irwall-1      
        IF (ikmidplane(ir,IKLO).NE.0.AND.
     .      ikmidplane(ir,IKHI).EQ.0) 
     .    rho(ir,CELL1) = -(ikintersec(ir,IKLO) - rhozero)
      ENDDO

      DO ir = 1, nrs
        WRITE(SLOUT,'(A,I4,1P,3E15.7)')
     .    'IR RHO = ',ir,(rho(ir,in),in=1,3)
      ENDDO
c
c
c
      CALL GridSpace2
c
c     OSM_DP5:
c
      DO ir = 1, nrs
        osm_dp5(IKLO,ir) = ksb(0,ir)
        osm_dp5(IKHI,ir) = ksmaxs(ir)
      ENDDO
c
c     Deal with boundary rings:
c
      DO ir = 1, nrs
        IF (idring(ir).EQ.-1) ksmaxs(ir) = MAX(ksmaxs(ir),1.0)
      ENDDO
c
c jdemod
c
c Changes were less than 1e-6m - not sure why this code is required
c
c...  Recalculate cell centers:
      recalculate = .FALSE.

      DO ir = 1, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        DO ik = 1, nks(ir)
          id = korpg(ik,ir)
          rtemp = 0.0
          ztemp = 0.0
          if (id.ne.0) then  

            DO in = 1, nvertp(id)
              rtemp = rtemp + rvertp(in,id)
              ztemp = ztemp + zvertp(in,id)
            ENDDO
            rtemp = rtemp/DBLE(nvertp(id))
            ztemp = ztemp/DBLE(nvertp(id))

c           Check for changes and record them

            if ( (abs(rtemp-rs(ik,ir)).gt.1.0e-4).or.
     >           (abs(ztemp-zs(ik,ir)).gt.1.0e-4)) then  

c                write(dbgunit,'(a,2i6,4(1x,g18.10))')
c     >              'CHANGE: NEW RS,ZS:',
c     >                ik,ir,rs(ik,ir),rtemp,zs(ik,ir),ztemp

              
              IF (.NOT.recalculate.AND.sloutput) 
     .          WRITE(0,*) 'NOTE: CELL CENTER RECALCULATION '//
     .                     'APPEARS NECESSARY'
              recalculate = .TRUE.
              
c              WRITE(0,*) 'HALTING CODE'
c              STOP
              WRITE(6,'(A,4I6,4F10.6)') 
     .  'CENCALC:',ik,ir,irsep,irtrap,rs(ik,ir),rtemp,zs(ik,ir),ztemp
             

              rs(ik,ir)=rtemp 
              zs(ik,ir)=ztemp


            endif 


          endif

        ENDDO
      ENDDO
c
c...  Recalculate cell centers:
c      DO ir = 1, nrs
c        DO ik = 1, nks(ir)
c          id = korpg(ik,ir)
c          rs(ik,ir) = 0.0
c          zs(ik,ir) = 0.0
c          DO in = 1, nvertp(id)
c            rs(ik,ir) = rs(ik,ir) + 1.0/REAL(nvertp(id))*rvertp(in,id)
c            zs(ik,ir) = zs(ik,ir) + 1.0/REAL(nvertp(id))*zvertp(in,id)
c          ENDDO
c        ENDDO
c      ENDDO
c jdemod
c


      RETURN
99    WRITE(EROUT,'(A,2I4)') 'IR   IR1   = ',ir,ir1
      WRITE(EROUT,'(A,2I4)') 'IKO  IKI   = ',iko,iki
      WRITE(EROUT,'(A,2I4)') 'IKTO IKTI  = ',ikto,ikti

      zvertp(4,id1) = 0.0
      zvertp(4,id2) = 0.0

      WRITE(0,*) 'IRs:',irsep,irsep2,ir
      WRITE(0,*) 'IKs:',ikto,ikti
      WRITE(0,*) 'IKs:',ikto2(irsep2),ikti2(irsep2)
      CALL DumpGrid('Connected bastard')
      STOP 
      END
c
c ======================================================================
c
      SUBROUTINE InitializeVariables
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slcom'
      INCLUDE 'cgeom'
      include 'cadas'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'cedge2d'  

      INCLUDE 'solparams'
      INCLUDE 'solswitch'
      INCLUDE 'solcommon'

      include 'reiser_com'  
c
c      include 'line_profile' 
c
      include 'driftvel'
      include 'fperiph_com'
      include 'dperpz'
      include 'slcom_sol28'
c

      INTEGER i1,i2,ir

      COMMON /OLDPLASMA/ oldknbs ,oldktebs ,oldktibs ,oldkvhs ,
     .                   oldknbs2,oldktebs2,oldktibs2,oldkvhs2
      REAL oldktebs (MAXNKS,MAXNRS),oldktibs (MAXNKS,MAXNRS),
     .     oldknbs  (MAXNKS,MAXNRS),oldkvhs  (MAXNKS,MAXNRS),
     .     oldktebs2(MAXNKS,MAXNRS),oldktibs2(MAXNKS,MAXNRS),
     .     oldknbs2 (MAXNKS,MAXNRS),oldkvhs2 (MAXNKS,MAXNRS)


      COMMON /OPTTEMP/ osm_matcht,forcet1
      INTEGER          osm_matcht,forcet1

      COMMON /PININIT/ pininit
      LOGICAL          pininit(MAXNRS)

      COMMON /ERRNOTE/ con_nerr,con_errc,con_errir
      INTEGER          con_nerr,con_errc(MAXNRS),con_errir(MAXNRS)

      COMMON /CFSCOM/ cfs_mage,cfs_magi,cfs_sume,cfs_sumi
      REAL            cfs_mage(MAXNRS),cfs_magi(MAXNRS),
     .                cfs_sume(MAXNRS),cfs_sumi(MAXNRS)

      COMMON /INPUTCHK/ inputflag
      INTEGER           inputflag(100)

      COMMON /STATCOM/ osm_ncnt ,osm_nerr,osm_temod,osm_relfr,osm_cadj,
     .                 osm_tcon
      INTEGER          osm_ncnt (MAXNRS)    ,osm_nerr (2,MAXNRS),
     .                 osm_cadj (5,2,MAXNRS),osm_tcon (MAXNRS)
      REAL             osm_temod(MAXNRS)    ,osm_relfr(MAXNRS)

      COMMON /PEAKSHAPECOM/ ipeakshape
      INTEGER               ipeakshape(2)

      COMMON /OSMSHAPE/ osmsplim   ,osmpksym
      REAL              osmsplim(2),osmpksym(2)

      COMMON /OSMOPTS/ osm_fixopt
      INTEGER          osm_fixopt(2,MAXNRS)

      COMMON /MACHCOM/ machine2
      CHARACTER*64     machine2

      COMMON /TIMESTATS/ timeint
      REAL               timeint(100)

      COMMON /OSMPMKCOM/ save_osmpmk,scaleplateau
      REAL*8             save_osmpmk(0:MAXNKS,MAXNRS)
      LOGICAL            scaleplateau(MAXNRS)

      INTEGER    MAXTAG
      PARAMETER (MAXTAG=1000)
      COMMON /INPCHK/ ntaglist,taglist
      INTEGER     ntaglist
      CHARACTER*3 taglist(MAXTAG)

      ntaglist = 0
      DO i1 = 1, MAXTAG
        taglist(i1) = '   '
      ENDDO

      WRITE(machine2,'(64X)')

      CALL LoadEIRENEAtomicData

      CALL ISet(nvertp,MAXNKS*MAXNRS,0)


      CALL RZero(flxhw7,MAXSEG)

      CALL RZero(timeint,100)

      CALL IZero(virtag,MAXNKS*MAXNRS)

      sldebug = 0
      connected = .FALSE.

      rvertp = 0.0
      zvertp = 0.0

      idring = TARTOTAR
      irorg2 = 0
      ikmidplane = 0
      psindat = 0
      psidat = 0.0
      pincode = 0

      CALL ISet(disindex,MAXNRS,-1)
      CALL RZero(dislev,MAXNRS)

      indasd = 0
      CALL RZero(pinasd,MAXASCDAT*MAXASD2*MAXASS*2)
c...  Set additional cell "plasma" values to EIRENE vacuum defaults:
      DO i1 = 1, MAXASS
        DO i2 = 1, MAXASCDAT
          pinasd(i2,1,i1,1) = 1.0E+01
          pinasd(i2,2,i1,1) = 1.0E-02
          pinasd(i2,3,i1,1) = 0.0
          pinasd(i2,4,i1,1) = 0.0
          pinasd(i2,5,i1,1) = 0.0
        ENDDO
      ENDDO

      asc_ncell = 0
      CALL IZero(asc_cell  ,MAXASC)
      CALL IZero(asc_region,MAXASC)
      CALL IZero(asc_link  ,MAXASC*4)
      CALL IZero(asc_grid  ,MAXASC*2)
      CALL IZero(asc_nvp   ,MAXASC)
      CALL RZero(asc_rvp   ,MAXASC*8)
      CALL RZero(asc_zvp   ,MAXASC*8)
      CALL RZero(ascdata   ,MAXASCDAT*5)

      CALL ISet(idring,MAXNRS,TARTOTAR)

      CALL IZero(supflx,MAXNRS*2)

      CALL RSet(mulrec,MAXNKS*MAXNRS,1.0)
      CALL RSet(mulrecl,MAXNRS,1.0)
      CALL RSet(mulion,MAXNKS*MAXNRS,1.0)
      CALL RSet(mulqer,MAXNKS*MAXNRS,1.0)
      CALL RSet(mulqei,MAXNKS*MAXNRS,1.0)


      CALL IZero(osm_ncnt ,MAXNRS)
      CALL IZero(osm_nerr ,MAXNRS*2)
      CALL IZero(osm_cadj ,MAXNRS*2*5)
      CALL IZero(osm_tcon ,MAXNRS)
      CALL RZero(osm_temod,MAXNRS)
      CALL RZero(osm_relfr,MAXNRS)


      CALL IZero(ikbreak,MAXNRS)

      DO i1 = 1, MAXNRS
        pininit(i1) = .FALSE.
      ENDDO

      WRITE(50,*) 'Initialising unstructured input options: '
      WRITE(50,*) ' '


      cgridst  = 0

      CALL IZero(eiraout,MAXASD*10) 
 
      CALL RSet(eirscale,100,1.0)

      setmach = 0.0
      muldensity = 1.0
      addtemp = 99.0

      eirntrans = 0
      eirniontime = 0

c...  Remove from common block.  Replaced with EIRNTORSEG.
      eirnttra = 0

      eirnaout = 1
      eiraout(1,1) =  2
      eiraout(1,2) = -MAXASCDAT
      eiraout(1,3) =  MAXASCDAT
      eiraout(1,4) =  0
      eiraout(1,5) =  0

      eirntally  =  0
      DO i1 = 1, MAXTALLY
        DO i2 = 1, 10
          WRITE(eirtally(i1,i2),'(72X)')
        ENDDO
      ENDDO
      eirphoton  =  0
      eirtime    =  0
      eirtrim    =  0
      eirdata    = -1
      eirgeom    = -1
      eirgrid    = -1
      eiradd     = -1
      eirneut    = -1
      eirniter   =  1
      eirdebug   =  0
      eirmat1    =  0
      eirmat2    =  0
      eirtemp1   =  0.0
      eirtemp2   =  0.0
      eirnpgdat  =  0
      eirnasdat  =  0
      eirnspdat  =  0
      eirspdatmode = 0
      eirnpuff   =  0
      eirbgk     =  0
      eiropacity =  0
      eircxd2    =  0
      eiracx     =  1
      eirpmode   =  0
c... Should really be 1?
      eird2ion = 1
      eird2dis = 1
      eirph2 = 0
      eiralloc   = 1.0
      eirnstrata = 0
      eirnsection = 1
      CALL IZero(eirnlimi,5000)
      CALL RZero(eirstrata,MAXSTRATA*10)
      eirntorseg = 0
      eirermin   =  1.0
      eirrinteg  =  0.0
      eireinteg  =  0.0
      eirtrc1    =  0
      eirtrc2    =  0
      eird2fu    =  0
      eirzaa     = -1.0
      eirtorfrac =  1.0
      eirsrcmul  =  1.0
      eirfuji    =  0
      eirnsdtor  = 1
      CALL RZero(eirsdtor,MAXTOR)

c...  Time-dependent mode parameters:
      eirdtimv = 0.0

      vacnseg   = 0
      CALL IZero(vacregion,MAXVACREGION)

      CALL RZero(eirpgdat,MAXNAS*MAXASD)
      CALL RZero(eirasdat,MAXNAS*MAXASD)
      CALL RZero(eirspdat,MAXNAS*MAXASD)

      eirnres = 0
      CALL RZero(eirres,6*7*MAXPINITER)
c      CALL RZero(eirres,20*MAXPINITER)

      cmodopt  = 0
      stagopt  = 0

      thesis = .FALSE.

      stopopt  = 0
      stopopt2 = 0
      stopopt3 = 0

      CALL IZero(iflexopt,MAXFLEX)
      CALL RZero(rflexopt,MAXFLEX)
c...Presently, allowed over/underionisation in AdjustParticleSource,
c   i.e. +-10%:
      rflexopt(1) = -99.0
      rflexopt(2) = 0.00
c      rflexopt(2) = 0.50
      rflexopt(3) = 1.00
      rflexopt(4) = 1.00
      rflexopt(5) = 1.00


      outwall   = 1
      outplasma = 1
      outcell   = 1
      outpoly   = 1
      outgeom   = 1
      outtarget = 1
      outpin    = 1
      outmode   = 0

      out_source = 0
      out_plasma = 0
      out_geom   = 0

      tarninter = NULL
      tarinter = 0.0
      tarintermode = 0.0

      psitarg = 0.0

      tarsource  = 0
      tarshift(IKLO) = 0.0
      tarshift(IKHI) = 0.0

      haldata    = 0

      grd_minpl = 1.0E-4
      grd_refine = 0
      grd_range = 1.0
      grd_shift = 0.0

      firstshape = .TRUE.
      ipeakshape(IKLO) = MAXNKS
      ipeakshape(IKHI) = MAXNKS

      osmbulkn = 0.0
      osmbulkte = 0.0
      osmbulkti = 0.0
      osmbulkv = 0.0


      osmsplim(IKLO) = 0.0
      osmsplim(IKHI) = 0.0
      osmpksym(IKLO) = 0.0
      osmpksym(IKHI) = 0.0

      osm_range  = 1.0
      osm_store  = -1
      osm_watch1 = 0
      osm_watch2 = 0
      osm_probe  = 0
      osm_nflag  = 0
      osm_mode   = 0
      osm_matchs = 0
      osm_matchp = 0
      osm_preopt = 0                     
      osm_recopt = 0

      CALL IZero(osmikshift,2*MAXNRS)

      osm_testep = 0.0

      s21_mode   = 0

      s28mode = 0.0
      s28ion = 0
      s28rec = 1
      s28cfp = 0
      s28cfpnt = 0
      s28cfpdrft = 0
      s28ionpfz = 0
      s28recpfz = 1
      s28cfppfz = 1
      s28cfppfznt = 0
      s28cfppfzdrft = 0
      s28mom = 0
      s28cfm = 0
      s28mompfz = 0
      s28cfmpfz = 0
      s28te = 0
      s28ti = 0
      s28tipfz = 0
      s28tiratio = 1.0
      s28fit = 0
      s28superdet = 1
      s28superdetpfz = 1
      s28probe = 0
      s28probepfz = 0
      s28nemode = 0
      s28nemodepfz = 0
      s28momfr = 0.3
      s22pfzion = 0
      s28b34 = 2.0 / 4.0
      s28b34det = 3.0 / 4.0
      s28momTe = 10.0
      s28sym = 1

      s28recfrac = 1.0
      s28ionfrac = 1.0
      s28momfrac = 1.0

      s28ionset = .FALSE.
      s28recset = .FALSE.
      s28momset = .FALSE.
      s28ionsetpfz = .FALSE.
      s28recsetpfz = .FALSE.
      s28momsetpfz = .FALSE.

c...  slcom_sol28:
      new = .FALSE.

      rel_iter   =  1
      rel_count  = -1
      rel_opt    =  0
      relmode    =  0
      rel_frac   =  1.0
      rel_nstep  =  0
      rel_step   =  0
      rel_niter  =  0
      rel_pace   =  1.0
      rel_ndata  =  0
      rel_bound1 =  0.0
      rel_bound2 =  1.0
      relexpdat  =  0.0
      osm_symopt =  0
      rel_tol    =  1.0
      osmmock    =  0

      prb_align = 0
      prb_shift = 0.0

      DO i1 = 1, MAXNRS

        rel_deltati   (i1) = 0.0
        rel_hte       (i1) = 0.0
        rel_hti       (i1) = 0.0
        rel_hne       (i1) = 0.0
        rel_symfr     (i1) = 1.0
        rel_prbfr(IKLO,i1) = 1.0
        rel_prbfr(IKHI,i1) = 1.0

        osm_fixopt(IKLO,i1) = 0
        osm_fixopt(IKHI,i1) = 0

        DO i2 = 1, 3
          rel_hproe(i2,i1) = HI
          rel_hproi(i2,i1) = HI
        ENDDO

        osm_model(IKLO,i1) = 0
        osm_model(IKHI,i1) = 0

        osm_dp1(IKLO,i1) = 0.15
        osm_dp1(IKHI,i1) = 0.15

        osm_dp2(IKLO,i1) = 0.0
        osm_dp2(IKHI,i1) = 0.0

        osm_dp3(IKLO,i1) = 1.0
        osm_dp3(IKHI,i1) = 1.0

        osm_dp5(IKLO,i1) = 0.0
        osm_dp5(IKHI,i1) = 0.0

        osm_code(IKLO,i1) = 1
        osm_code(IKHI,i1) = 1

        DO i2 = 1, MAXNKS
          osm_dp4(i2,i1) = 1.0
          osm_dp6(i2,i1) = 1.0
        ENDDO




        rel_qemul (IKLO,i1) = 1.0
        rel_qemul (IKHI,i1) = 1.0
c...temp1
        osm_peimul(IKLO,i1) = -0.1
        osm_peimul(IKHI,i1) = -0.1

        cfs_mage(i1) = 0.0
        cfs_magi(i1) = 0.0
        cfs_sume(i1) = 0.0
        cfs_sumi(i1) = 0.0
      ENDDO

      nlpdati2  =  0
      nlpdato2  =  0

      adp_opt    = 0
      adp_region = 3
      adp_upper  = HI
      adp_lower  = 1.0

      con_nerr = 0
      DO i1 = 1, MAXNRS
        cmachno(i1,1) = 1.0
        cmachno(i1,2) = 1.0

        con_errc (i1) = 0
        con_errir(i1) = 0

        idds(i1,1) = MAXNDS
        idds(i1,2) = MAXNDS

        ikerror(IKLO,i1) = MAXNKS
        ikerror(IKHI,i1) = MAXNKS
      ENDDO

      DO i1 = 1, MAXNRS
        DO i2 = 1, MAXNKS
          korpg(i2,i1) = MAXNKS*MAXNRS
        ENDDO
      ENDDO

      CALL RZero(lpdati2  ,MAXINS*9)
      CALL RZero(lpdato2  ,MAXINS*9)

      DO i1 = 1, MAXINS
        lpdati2(i1,8) = 1.0
        lpdati2(i1,9) = 1.0
        lpdato2(i1,8) = 1.0
        lpdato2(i1,9) = 1.0

        DO i2 = IKLO, IKHI
          rel_mfsp(i2,i1) = 1.0
        ENDDO
      ENDDO

      CALL RSet (rel_data ,MAXINS*15,1.0)
      CALL IZero(osm_sympt,MAXNRS)
      CALL IZero(rel_viter,(1+MAXSTEP))

      CALL RZero(oldktebs,MAXNKS*MAXNRS)
      CALL RZero(oldktibs,MAXNKS*MAXNRS)
      CALL RZero(oldknbs ,MAXNKS*MAXNRS)
      CALL RZero(oldkvhs ,MAXNKS*MAXNRS)

      s21_ndatai = 0
      s21_ndatao = 0

      CALL RZero(s21_datai,MAXNRS*9)
      CALL RZero(s21_datao,MAXNRS*9)

      CALL RZero(bgplasopt,2*MAXNRS*9)

c      switch(SWPOW2) = 1.0
c      switch(SWPOW3) = 4.0
c      switch(SWION2) = 1.0

      CALL RZero(prp_te,MAXNRS*NUMPRB)
      CALL RZero(prp_ti,MAXNRS*NUMPRB)
      CALL RZero(prp_ne,MAXNRS*NUMPRB)

      CALL RZero(dip_te,MAXPRB*NUMPRB)
      CALL RZero(dip_ti,MAXPRB*NUMPRB)
      CALL RZero(dip_ne,MAXPRB*NUMPRB)
      CALL RZero(dip_s ,MAXPRB*NUMPRB)
      CALL RZero(dip_v ,MAXPRB*NUMPRB)

      CALL IZero(prb_num,NUMPRB)

      CALL RZero(prb_te ,MAXPRB*NUMPRB)
      CALL RZero(prb_ti ,MAXPRB*NUMPRB)
      CALL RZero(prb_ne ,MAXPRB*NUMPRB)
      CALL RZero(prb_rho,MAXPRB*NUMPRB)
      CALL RZero(prb_r  ,MAXPRB*NUMPRB)
      CALL RZero(prb_z  ,MAXPRB*NUMPRB)


      CALL RZero(pinior ,MAXNKS*MAXNRS)
      CALL RZero(pinior2,MAXNKS*MAXNRS)
      CALL RZero(pinmpr ,MAXNKS*MAXNRS)
      CALL RZero(pinmpr2,MAXNKS*MAXNRS)
      CALL RZero(pinqir ,MAXNKS*MAXNRS)
      CALL RZero(pinqir2,MAXNKS*MAXNRS)
      CALL RZero(pinqer ,MAXNKS*MAXNRS)
      CALL RZero(pinqer2,MAXNKS*MAXNRS)

      CALL RZero(pinion  ,MAXNKS*MAXNRS)
      CALL RZero(pinrec  ,MAXNKS*MAXNRS)
      CALL RZero(pinqi   ,MAXNKS*MAXNRS)
      CALL RZero(pinqe   ,MAXNKS*MAXNRS)
      CALL RZero(pinmol  ,MAXNKS*MAXNRS)
      CALL RZero(pinena  ,MAXNKS*MAXNRS)
      CALL RZero(pinenz  ,MAXNKS*MAXNRS)
      CALL RZero(pinenm  ,MAXNKS*MAXNRS)
c      CALL RZero(pinvdist,...)
      CALL RZero(pinmp   ,MAXNKS*MAXNRS)
      CALL RZero(pinatom ,MAXNKS*MAXNRS)
      CALL RZero(pinalpha,MAXNKS*MAXNRS)
      CALL RZero(pinline(1,1,1,H_BALPHA),MAXNKS*MAXNRS*6)
      CALL RZero(pinline(1,1,1,H_BGAMMA),MAXNKS*MAXNRS*6)
      CALL RZero(pinionz ,MAXNKS*MAXNRS)
      CALL RZero(pinz0   ,MAXNKS*MAXNRS)
      CALL RZero(pindata ,MAXNKS*MAXNRS*MAXDATA)
      CALL RZero(pindata2,MAXNKS*MAXNRS*MAXDATA)

      CALL RZero(pinbgk  ,MAXNKS*MAXNRS*MAXBGK*MAXTOR)
      CALL RZero(pinbgk2 ,MAXNKS*MAXNRS*MAXBGK)



      CALL RZero(osmpei     ,MAXNKS*MAXNRS)
      CALL RZero(osmcde     ,MAXNKS*MAXNRS)
      CALL RZero(osmcdi     ,MAXNKS*MAXNRS)
      CALL RZero(osmcve     ,MAXNKS*MAXNRS)
      CALL RZero(osmcvi     ,MAXNKS*MAXNRS)
      CALL RZero(osmpmk (0,1),(MAXNKS+1)*MAXNRS)
      CALL DZero(osmpmk2(0,1),(MAXNKS+1)*MAXNRS)
      CALL DZero(save_osmpmk(0,1),(MAXNKS+1)*MAXNRS)
      DO ir = 1, MAXNRS
        scaleplateau(ir) = .FALSE.
      ENDDO

      CALL RZero(pinion2  ,MAXNKS*MAXNRS)
      CALL RZero(pinrec2  ,MAXNKS*MAXNRS)
      CALL RZero(pinqi2   ,MAXNKS*MAXNRS)
      CALL RZero(pinqe2   ,MAXNKS*MAXNRS)
      CALL RZero(pinmol2  ,MAXNKS*MAXNRS)
      CALL RZero(pinena2  ,MAXNKS*MAXNRS)
      CALL RZero(pinenz2  ,MAXNKS*MAXNRS)
      CALL RZero(pinenm2  ,MAXNKS*MAXNRS)
c      CALL RZero(pinvdist2,...)
      CALL RZero(pinmp2   ,MAXNKS*MAXNRS)
      CALL RZero(pinatom2 ,MAXNKS*MAXNRS)
      CALL RZero(pinalpha2,MAXNKS*MAXNRS)
      CALL RZero(pinionz2 ,MAXNKS*MAXNRS)
      CALL RZero(pinz02   ,MAXNKS*MAXNRS)


      CALL RZero(osmpei2,MAXNKS*MAXNRS)

      CALL RZero(osmcfp ,MAXNKS*MAXNRS)
      osmcfpflx = 0.0
      CALL RZero(osmcfi ,MAXNKS*MAXNRS)
      CALL RZero(osmcfe ,MAXNKS*MAXNRS)

      CALL RZero(osmmp  ,MAXNKS*MAXNRS)
      CALL RZero(osmmp2 ,MAXNKS*MAXNRS)
      CALL RZero(osmqe  ,MAXNKS*MAXNRS)
      CALL RZero(osmqe2 ,MAXNKS*MAXNRS)

      DO i1 = 1, MXSPTS
        DO i2 = 1, 3
          gtarg    (i1,i2) = 1.0
          elecptarg(i1,i2) = 1.0
          ionptarg (i1,i2) = 1.0
        ENDDO
      ENDDO

      DO i1 = 1, 100
        inputflag(i1) = 0
      ENDDO

c...  Tags:
      tagmulrec = .FALSE.
      tagpinatom = .FALSE.

c
c     Call routine to initialize the unstructured input data 
c     This is a separate routine since unlike normal variable 
c     initialization - the values set in this routine can have
c     a significant impact on the simulation and thus the default
c     values need to well documented.
c
      call InitializeUnstructuredInput

      CALL osm_InitializeOptions
c
      return
      END





