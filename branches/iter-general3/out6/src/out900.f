c     -*-Fortran-*-
c
      subroutine out900(iref,graph,iopt,ierr)
      implicit none
      integer iref,iopt,ierr
      character*(*) graph
c
      include 'params'
      include 'outcom'
c
c     Other common blocks
c
      include 'cgeom'
      include 'comtor'
c      include 'cneut2'
c      include 'dynam2'
c      include 'dynam3'
c      include 'dynam4'
c      include 'pindata'
      include 'cadas'
c      include 'grbound'
c      include 'outxy'
c      include 'cedge2d'
c      include 'transcoef'
c      include 'cioniz'
c      include 'reiser' 
c      include 'printopt'  
c
c     Local Variables
c
      integer ik,ir,iz
      integer in,ii
c      integer ind
c
c     Local Variables
c
c
c
c     Thompson plots
c
      integer interp_opt,nsets,nthompson
      integer max_tcols
      parameter (max_tcols=5)
      real    plotdata(maxpts,max_tcols),
     >        thompson_pts(maxpts,max_tcols),
     >        axis_offset_r,axis_offset_z
      external grid_interpolate
      real grid_interpolate

 
      character*80 label 
c slmod begin
      CHARACTER dummy*5000
c slmod end

c
c     For multiple plots on the page - series 700
c
      real mvals(maxdatx,maxplts,maxngs)
      real mouts(maxdatx,maxplts,maxngs)
      real mwids(maxdatx,maxplts,maxngs)
      integer  ringnos(maxplts),ip,ig,ip2,pnks(maxplts,maxngs)
      integer  mdrawtype(maxplts,maxngs)
      integer  drawtype(maxngs)
      integer  axistype
      integer  sctype,ngrm
      integer  pngs(maxplts)  
      character*36 mlabs(maxplts,maxngs)
      character*36 pltlabs(maxplts)
c

      integer i1,i2
      INTEGER ntime,stime,dtime


c slmod begin
c
c     Possible bug - this array is used in call to ADASRD to store the 
c     coefficients returned - but if it isn't big enough - will get 
c     problems.  
c
      real  dum1(MAXTHE)

      CHARACTER*50 cdum1,tag
      CHARACTER*174 title2

      INTEGER nenum,tenum,opt_const,plot_mode(30)
      REAL    nemin,nestep,temin,temax,neTe,frac1
      REAL plotte(MAXGXS,MAXNGS),plotne(MAXGXS,MAXNGS),
     .     plot_y (MAXGXS,MAXNGS),fac32,fac42,fac52,rec,pop,
     .     fac21,fac31,fac41,fac51

      REAL GetEAD
c slmod end

      include 'colours'
c
c
c     Thompson comparison plots
c
      if (iref.lt.920) then
c
        if (iopt.eq.0) return
c
c       Load data for Thompson plots
c
        call rdg6(GRAPH,interp_opt,iseld,
     >            axis_offset_r,axis_offset_z,
     >            nsets,
     >            plotdata,maxpts,max_tcols,ierr)
c
c
c slmod begin
      ELSEIF (IREF.LE.999) THEN

        IF (IOPT.EQ.0) RETURN
        ref = graph(19:80)
c
        IF     (IREF.EQ.950) THEN
          READ(5,*) cdum1,iclass,cion,nizs,opt_xscale,opt_yscale,
     .              tenum,temin,temax,opt_const,nenum,nemin,nestep
        ELSEIF (IREF.EQ.956) THEN
          READ(5,*) cdum1,tag,stime,dtime,ntime
        ELSEIF (IREF.EQ.966) THEN
          call rdg4(graph4,ngrm,nplts,ringnos,maxplts,pltfact,ierr)

          do ip = 1,nplts
            write (pltlabs(ip),'(''Ring:'',i2)') ringnos(ip)
          end do
        ELSEIF (IREF.EQ.988) THEN
          CALL RDG2 (GRAPH2,ROBS,ZOBS,THEMIN,DTHE,THERES,NUMTHE,
     >                IZMIN,IZMAX,AVPTS,NUMSMOOTH,ATYPE,IERR)
c
c          if (cgridopt.eq.3) themin = themin + 180.0
c
          IF (IERR.EQ.1) THEN
             WRITE(6,*) 'RDG2 ERROR READING 980/988 SERIES- GRAPH '//
     .                  'DETAILS'
             IERR = 0
             RETURN
          ENDIF
        ENDIF
      endif
c


      call init_plot(iref,graph,iopt) 




c
C-----------------------------------------------------------------------
c
c     Thompson plots - plots of quantities at specified R,Z positions
c                      interpolated as necessary
c
C-----------------------------------------------------------------------
c
c
c     901 - plot of background density
c
      if (iref.eq.901) then
c
c        RDG6 has been called to load additional information.
c
c
c        Calculate axis for RZ data points
c
         if (nsets.gt.0) then
c
            plotdata(1,3) = 0.0
c
            do in = 2,nsets
c
               plotdata(in,3) = plotdata(in-1,3) +
     >                    sqrt((plotdata(in,2)-plotdata(in-1,2))**2+
     >                         (plotdata(in,1)-plotdata(in-1,1))**2)
c
            end do
c
         endif

c
c        Load experimental data set
c
         if (iseld.gt.0) then
c
            call load_rzdata(iseld,nthompson,thompson_pts,
     >                      maxpts,max_tcols,axis_offset_r,
     >                      axis_offset_z,datatitle)
c
c
c           Copy Thompson RZ data over to sample set if nsets is zero
c
            if (nsets.le.0) then
c
               nsets = nthompson
c
               do in = 1,nthompson

                  plotdata(in,1) = thompson_pts(in,1)
                  plotdata(in,2) = thompson_pts(in,2)
                  plotdata(in,3) = thompson_pts(in,3)

               enddo

            endif

         endif
c
c        Check for error condition - only proceed if nsets is greater than
c        zero
c
         if (nsets.le.0) then
c
            write(6,*) 'ERROR IN PLOT 901 - NO DATA TO PLOT - '
     >                //'ABORTING PLOT'
            return
c
         endif
c
c
c        Loop through points and interpolate appropriate values
c        on grid.
c
         do in = 1,nsets
c
            plotdata(in,4) = grid_interpolate(
     >                              plotdata(in,1),plotdata(in,2),
     >                              interp_opt,knbs,1,1,knds,1)
c
            plotdata(in,5) = grid_interpolate(
     >                              plotdata(in,1),plotdata(in,2),
     >                              interp_opt,ktebs,1,1,kteds,1)
c
         end do
c
         write(6,*) 'DATA:',nsets
         do in = 1,nsets
            write(6,'(5(1x,g12.5))') (plotdata(in,ik),ik=1,5)
         end do
c
c        Load data arrays for plotting - set number of curves to plot,
c        calculate axes.
c
c
c        Generate plots
c
c        Place the contents of corresponding plots on a common plot
c        i.e. plots of ne and Te
c
c
         ELABS(1) = 'DIV DIV '
         Len = lenstr(datatitle)
         ELABS(2) = 'THOMTHOM '//datatitle(1:Len)
c
         PLTLABS(1) =  'Ne         '
         PLTLABS(2) =  'Te         '
c
         XLAB = 'DIST (M) ALONG PROBE LINE'
         YLAB   = 'DIVIMP and Probe Results'
c
         NPLOTS = NPLOTS + 1
         WRITE (IPLOT,9012) NPLOTS,REF
c
         CALL rzero (mvals, maxnks*maxngs*maxplts)
c
         sctype = iopt
         ngrm  = 2
         nplts = 2
c
c        Set to plot only DIVIMP or DIVIMP + THOMSON data
c
         if (iseld.gt.0) then
            ngs = 2
         else
            ngs = 1
         endif
c
c        Loop through the data for the 2 plots
c
c        Loading one set of information for each plot from each
c        source - DIVIMP and EXPERIMENT
c
         do ip = 1,nplts
c
c           Load DIVIMP first
c
            pnks(ip,1) = nsets
c
            do ik = 1,nsets
c
               mouts(ik,ip,1) = plotdata(ik,3)

               mvals(ik,ip,1) = plotdata(ik,3+ip)
c
            end do
c
            if (ngs.gt.1) then
c
c              Load Experimental data next
c
               pnks(ip,2) = nthompson
c
               do ik = 1,nthompson
c
                  mouts(ik,ip,2) = thompson_pts(ik,3)
c
                  mvals(ik,ip,2) = thompson_pts(ik,3+ip)
c
               end do
c
           endif
c
         end do
c
c        Set drawing style for the different sets of data
c
         drawtype(1) = 1
         drawtype(2) = 1
C
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
               mdrawtype(ip,ig) = drawtype(ig)
            end do  
            pngs(ip) = ngs
         end do
c
         CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,ngs,
     >              mdrawtype,1)
c
      endif
c




c slmod begin
      IF (iref.EQ.950) THEN
c
c       cion  = atomic charge
c       nizs  = ionisation state
c       temin = lower limit of temperature range
c       temax = upper limit of temperature range
c       nemin = lower limit of density rang
c       nemax = upper limit of density range
c       nenum = number of density values to plot
c
C        ICLASS = 1: RECOMBINATION RATE COEFFICIENT
C                 2: IONISATION RATE COEFFICIENT
C                 3: CHARGE EXCHANGE RECOMBINATION COEFFICIENT
C                 4: POWER COEF. FOR RECOMBINATION AND BREMSSTRAHLU
C                 5: POWER COEFFICIENT FOR LINE RADIATION
C                 6: POWER COEFFICIENT FOR CHARGE EXCHANGE
C                 added by Krieger, IPP 5/95
C                 7: PHOTON EMISSIVITY FOR DIAGNOSTIC LINES
c
c Amjuel database support:
c        ICLASS = -1:
c                 -2:
c                 -3:
c
c
c
c
c
        IF (tenum.GT.MAXGXS) CALL ER('950','MAXGXS exceeded',*9999)

        CALL RZero(plot_y ,MAXGXS*MAXNGS)
        CALL RZero(dwids  ,MAXNDS+2)

        slopt = 1

        CALL WN('950','No other plots allowed for this OUT execution')

        ylab2(1)  = '<sigma-v>rec (m3 s-1)   '
        ylab2(2)  = '<sigma-v>ion (m3 s-1)   '
        ylab2(3)  = '<sigma-v>cex (m3 s-1)   '
        ylab2(4)  = 'Power coef. (...)       '
        ylab2(5)  = 'Power coef. (J s-1 m-3) '
        ylab2(6)  = 'Power coef. (...)       '
        ylab2(7)  = 'Emissivity (J m-2 s-1)  '
        ylab2(15) = 'SXB (ionisations/photon)'
        ylab2(16) = 'SXB (recomb./photon)    '
        ylab2(17) = 'Helpi (eV) (Harrison)   '
        ylab2(18) = 'Energy / ionisation (eV)'
        ylab2(19) = 'Emissivity (ph m3 s-1)  '
        ylab2(20) = 'Emissivity (ph m3 s-1)  '
        ylab2(21) = 'Emissivity (ph m3 s-1)  '

c...    AMJUEL data:
        ylab2(-1) = ' ...                    '
        ylab2(-2) = ' ...                    '
        ylab2(-3) = ' ...                    '
        ylab2(-4) = ' ...                    '

        plot_mode(15) = 1
        plot_mode(19) = 1
        plot_mode(20) = 2
        plot_mode(21) = 3

        xlab  = 'Te (eV)'
        ylab  = ylab2(iclass)

        title2 = '                                                  '


c...    Read plot options:

c...    Get the size of the plot, if non-standard:
        READ(5,'(A256)') dummy
        IF   (dummy(8:11).EQ.'Size'.OR.dummy(8:11).EQ.'size'.OR.
     .        dummy(8:11).EQ.'SIZE') THEN
          READ(dummy,*) cdum1,map1x,map2x,map1y,map2y
        ELSE
          BACKSPACE 5
        ENDIF


        IF (iclass.GE.1) THEN
c...      Initialize ADAS data:
          write(year,'(i2.2)') adasyr
          call xxuid(adasid)
          ii = LEN_TRIM(adasid)

          title2 = 'Atomic Data and Analysis Structure'

          WRITE(char(10),'(A,I4)') 'NUCLEAR CHARGE   = ',cion
          WRITE(char(9) ,'(A,I4)') 'IONISATION STATE = ',nizs
          WRITE(char(7) ,'(A)')    'ADAS ID   = "'//adasid(1:ii)//'"'
          WRITE(char(6) ,'(A,I4)') 'ADAS YEAR = ',adasyr
          WRITE(char(5) ,'(A)')    'ADAS EXT  = "'//adasex//'"'
          WRITE(char(4) ,'(A,I4)') 'ISELE     = ',isele
          WRITE(char(3) ,'(A,I4)') 'ISELR     = ',iselr
          WRITE(char(2) ,'(A,I4)') 'ISELX     = ',iselx
          WRITE(char(1) ,'(A,I4)') 'ISELD     = ',iseld
        ELSE
          title2 = 'AMJUEL EIRENE Database'

          CALL LoadEIRENEAtomicData

          slopt2 = 2

          call setup_col(ncols,5)
        ENDIF
c 
c
        IF     (opt_const.EQ.1) THEN
          DO i1 = 1, nenum
            DO i2 = 1, tenum
              frac1 =  REAL(i2 - 1) / REAL(tenum - 1)
              IF (opt_xscale.EQ.2) frac1 = frac1**2.0

              plotte(i2,i1) = temin + frac1 * (temax - temin)
              plotne(i2,i1) = nemin + nestep * REAL(i1-1)
            ENDDO
            WRITE(elabs(i1),'(4X,1P,E10.2,A)') plotne(1,i1),' m-3'
          ENDDO
        ELSEIF (opt_const.EQ.2) THEN
          DO i1 = 1, nenum
            neTe = nemin * nestep**REAL(i1-1)
            DO i2 = 1, tenum
              frac1 =  REAL(i2 - 1) / REAL(tenum - 1)
              IF (opt_xscale.EQ.2) frac1 = frac1**2.0

              plotte(i2,i1) = temin + frac1 * (temax - temin)
              plotne(i2,i1) = neTe / plotte(i2,i1)
            ENDDO
            WRITE(elabs(i1),'(4X,1P,E10.2,A)') neTe,' m-3 eV'
          ENDDO
        ELSEIF (opt_const.EQ.4) THEN
          DO i1 = 1, nenum
            DO i2 = 1, tenum
              frac1 =  REAL(i2 - 1) / REAL(tenum - 1)
              IF (opt_xscale.EQ.2) frac1 = frac1**2.0

              plotte(i2,i1) = temin + frac1 * (temax - temin)
              
              plotne(i2,i1) = nemin * (nestep**REAL(i1-1))
            ENDDO
            WRITE(elabs(i1),'(4X,1P,E10.2,A)') plotne(1,i1),' m-3'
          ENDDO

        ELSEIF (opt_const.EQ.3) THEN
c...      Density hard coded:
          nenum = 4
          plotne(1,1) = 1.0E+20
          plotne(1,2) = 5.0E+20
          plotne(1,3) = 1.0E+21
          plotne(1,4) = 1.5E+21
          DO i1 = 1, nenum
            DO i2 = 1, tenum
              frac1 =  REAL(i2 - 1) / REAL(tenum - 1)
              IF (opt_xscale.EQ.2) frac1 = frac1**2.0

              plotte(i2,i1) = temin + frac1 * (temax - temin)
              plotne(i2,i1) = plotne(MAX(1,i2-1),i1)
            ENDDO
            WRITE(elabs(i1),'(4X,1P,E10.2,A)') plotne(1,i1),' m-3'
          ENDDO

        ELSE
          DO i1 = 1, nenum
            DO i2 = 1, tenum
              frac1 =  REAL(i2 - 1) / REAL(tenum - 1)
              IF (opt_xscale.EQ.2) frac1 = frac1**2.0

              plotte(i2,i1) = temin + frac1 * (temax - temin)
              plotne(i2,i1) = nemin * nestep**REAL(i1-1)
            ENDDO
            WRITE(elabs(i1),'(4X,1P,E10.2,A)') plotne(1,i1),' m-3'
          ENDDO
        ENDIF
c
c
c
c...    Temp:
        iopt_ghost = 1

        IF     (iclass.EQ.-1) THEN
c...      AMJUEL: 
          DO i1 = 1, nenum
            plottype(i1) = i1 + 1
c            plottype(i1+nenum) = i1 + 1
            DO i2 = 1, tenum
              plotne(i2,i1) = MAX(1.0E+13,plotne(i2,i1))

              plot_y(i2,i1) = GetEAD(plotte(i2,i1),
     .                               plotne(i2,i1),3,'H.4 ') /
     .                        GetEAD(plotte(i2,i1),
     .                               plotne(i2,i1),4,'H.4 ')


c              plot_y(i2,i1) = GetEAD(plotte(i2,i1),
c     .                               plotne(i2,i1),24,'H.4 ') /
c     .                        GetEAD(plotte(i2,i1),
c     .                               plotne(i2,i1),23,'H.4 ')

c              plot_y(i2,i1) = GetEAD(plotte(i2,i1),
c     .                               plotne(i2,i1),24,'H.4 ')*1.0E-6 
c              plot_y(i2,i1+nenum) = GetEAD(plotte(i2,i1),
c     .                               plotne(i2,i1),23,'H.4 ')*1.0E-6


c              plot_y(i2,i1) = GetEAD(plotte(i2,i1),
c     .                               plotne(i2,i1),4,'H.4 ')*1.0E-6

c              plot_y(i2,i1) = GetEAD(plotte(i2,i1),
c     .                               plotne(i2,i1),3,'H.4 ')*1.0E-6
            ENDDO
          ENDDO

c          nenum = nenum * 2

        ELSEIF (iclass.EQ.-2) THEN
c...      AMJUEL: e + H2 reactions:
          ylab = "e + H2 -> ..."//elabs(i1)(1:LEN_TRIM(elabs(i1)))
          DO i1 = 1, nenum
            DO i2 = 1, tenum
             plotne(i2,        i1) = MAX(1.0E+13,plotne(i2,i1))
             plot_y(i2,        i1) = GetEAD(plotte(i2,i1),
     .                                      plotne(i2,i1),19,'H.4 ')
             plotne(i2,  nenum+i1) = plotne(i2,i1)
             plotte(i2,  nenum+i1) = plotte(i2,i1)
             plot_y(i2,  nenum+i1) = GetEAD(plotte(i2,i1),
     .                                      plotne(i2,i1),13,'H.4 ')
             plotne(i2,2*nenum+i1) = plotne(i2,i1)
             plotte(i2,2*nenum+i1) = plotte(i2,i1)
             plot_y(i2,2*nenum+i1) = GetEAD(plotte(i2,i1),
     .                                      plotne(i2,i1),14,'H.4 ')
            ENDDO
          ENDDO
          DO i1 = 1, nenum
            WRITE(elabs(        i1),'(4X,A)') 'e + H2 -> e + H2+'
            WRITE(elabs(  nenum+i1),'(4X,A)') 'e + H2 -> e + H + H'
            WRITE(elabs(2*nenum+i1),'(4X,A)') 'e + H2 -> e + H + H+'
          ENDDO
          nenum = 3 * nenum
          DO i1 = 1, nenum
            plottype(i1) = i1 + 1
          ENDDO
        ELSEIF (iclass.EQ.-3) THEN
c...      AMJUEL: p + H2 -> p + H2
          DO i1 = 1, nenum
            DO i2 = 1, tenum
              plot_y(i2,        i1) = GetEAD(plotte(i2,i1),
     .                                       plotne(i2,i1),11,'H.1 ') * 
     .                                1.0E+04
              plotne(i2,  nenum+i1) = plotne(i2,i1)
              plotte(i2,  nenum+i1) = plotte(i2,i1)
              plot_y(i2,  nenum+i1) = GetEAD(plotte(i2,i1),
     .                                       plotne(i2,i1),12,'H.3 ')
              plotne(i2,2*nenum+i1) = plotne(i2,i1)
              plotte(i2,2*nenum+i1) = plotte(i2,i1)
              plot_y(i2,2*nenum+i1) = GetEAD(plotte(i2,i1),
     .                                       plotne(i2,i1),18,'H.3 ')
            ENDDO
          ENDDO
          DO i1 = 1, nenum
            WRITE(elabs(        i1),'(4X,A)') 'p + H2 -> p + H2  H.1'
            WRITE(elabs(  nenum+i1),'(4X,A)') 'p + H2 -> p + H2  H.3'
            WRITE(elabs(2*nenum+i1),'(4X,A)') 'p + H2 -> H + H2+'
          ENDDO
          nenum = 3 * nenum
          DO i1 = 1, nenum
            plottype(i1) = i1 + 1
          ENDDO

        ELSEIF (iclass.EQ.-4) THEN
c...      AMJUEL: e + H2+ reactions:
          ylab = "e + H2+ -> ..."//elabs(i1)(1:LEN_TRIM(elabs(i1)))
          DO i1 = 1, nenum
            DO i2 = 1, tenum
             plotne(i2,        i1) = MAX(1.0E+13,plotne(i2,i1))
             plot_y(i2,        i1) = GetEAD(plotte(i2,i1),
     .                                      plotne(i2,i1),15,'H.4 ')

             plotne(i2,  nenum+i1) = plotne(i2,i1)
             plotte(i2,  nenum+i1) = plotte(i2,i1)
             plot_y(i2,  nenum+i1) = GetEAD(plotte(i2,i1),
     .                                      plotne(i2,i1),17,'H.4 ')

             plotne(i2,2*nenum+i1) = plotne(i2,i1)
             plotte(i2,2*nenum+i1) = plotte(i2,i1)
             plot_y(i2,2*nenum+i1) = GetEAD(plotte(i2,i1),
     .                                      plotne(i2,i1),16,'H.4 ')

            ENDDO
          ENDDO
          DO i1 = 1, nenum
c            WRITE(elabs(        i1),'(4X,A)') 'e + H2+ -> e + H + H+'
c            WRITE(elabs(  nenum+i1),'(4X,A)') 'e + H2+ -> H + H'
c            WRITE(elabs(2*nenum+i1),'(4X,A)') 'e + H2+ -> 2e + H+ + H+'
          ENDDO
          nenum = 2 * nenum
c          nenum = 3 * nenum
          DO i1 = 1, nenum
            plottype(i1) = i1 + 1
          ENDDO

        ELSEIF (iclass.EQ.-5.OR.iclass.EQ.-6.OR.iclass.EQ.-7.OR.
     .          iclass.EQ.-8) THEN
c...      AMJUEL: ...
          ylab = '...'//elabs(i1)(1:LEN_TRIM(elabs(i1)))

c...      Balmer Gamma transition probability (from EIRENE source HALFA.F):          
         FAC32=2.530E+06
c...      Balmer Alpha transition probability (from EIRENE source HALFA.F):          
c          FAC32=4.410E7

c...      Load specific atomic data into entry NREAC:
c          CALL SLREAC(50,'AMJUEL  ','H.12','2.1.8d   ','OT ')          
c          CALL SLREAC(50,'AMJUEL  ','H.12','2.1.8a   ','OT ')

          DO i1 = 1, nenum
            DO i2 = 1, tenum
             plotne(i2,i1) = MAX(1.0E+13,plotne(i2,i1))

             pop = GetEAD(plotte(i2,i1),plotne(i2,i1),20,'H.12')
c             pop = GetEAD(plotte(i2,i1),plotne(i2,i1),21,'H.12')

             IF     (iclass.EQ.-5) THEN
               plot_y(i2,i1) = pop
             ELSEIF (iclass.EQ.-6) THEN
c...           Optically thin recombination rate:
               ylab = 'SxBrec (thin)'//elabs(i1)(1:LEN_TRIM(elabs(i1)))

               rec = GetEAD(plotte(i2,i1),plotne(i2,i1),3 ,'H.4 ')

               plot_y(i2,i1) = (rec * 1.0E-06 * plotne(i2,i1)) *
     .                                          plotne(i2,i1)
c               plot_y(i2,i1) = rec * 1.0E-06 * plotne(i2,i1) / 
c     .                         (pop * FAC32)
             ELSEIF (iclass.EQ.-7) THEN
c...           Optically thick recombination rate:
               ylab = 'SxBrec (thick)'//elabs(i1)(1:LEN_TRIM(elabs(i1)))

               rec = GetEAD(plotte(i2,i1),plotne(i2,i1),4 ,'H.4 ')

               plot_y(i2,i1) = (rec * 1.0E-6 * plotne(i2,i1)) *
     .                                         plotne(i2,i1)

c               plot_y(i2,i1) = rec * 1.0E-06 * plotne(i2,i1) / 
c     .                         (pop * FAC32)
             ELSEIF (iclass.EQ.-8) THEN
               ylab = 'Dgamma emission, rec (kW m-3)'

               plot_y(i2,i1) = pop * plotne(i2,i1) * FAC32 * 
     .                         (6.63E-34*3.0E+08)/(4340.0*1.0E-10) * 
     .                         1.0E-03
             ENDIF

            ENDDO
          ENDDO
c          DO i1 = 1, nenum
c            WRITE(elabs(        i1),'(4X,A)') '...'
c          ENDDO
          DO i1 = 1, nenum
            plottype(i1) = i1 + 1
          ENDDO

        ELSEIF (iclass.EQ.-9) THEN
c...      Examining (Balmer) line ratios for recombination signatures:

c...      AMJUEL: ...
          ylab = '...'//elabs(i1)(1:LEN_TRIM(elabs(i1)))

c...      Einstein coefficients for Balmer lines (from HALPHA.F in EIRENE):
          FAC32=4.410E+07
          FAC42=8.419E+06
          FAC52=2.530E+06

          FAC21=4.699E+08
          FAC31=5.575E+07
          FAC41=1.278E+07
          FAC51=4.125E+06

          IF (iopt.EQ.1) title2 = 'Ba_gamma:Ba_alpha'
          IF (iopt.EQ.2) title2 = 'Ba_gamma:Ba_beta '
          IF (iopt.EQ.3) title2 = 'Ba_gamma:B_beta:Ba_alpha'
          IF (iopt.EQ.4) title2 = 'Ba_beta:Ba_alpha'
          IF (iopt.EQ.5) title2 = 'Ly_alpha:Ba_alpha'
          IF (iopt.EQ.6) title2 = 'Ba_gamma PEC'

c...      Load specific atomic data into entry NREAC:
c          CALL SLREAC(50,'AMJUEL  ','H.12','2.1.8d   ','OT ')          

          DO i1 = 1, nenum
            DO i2 = 1, tenum
              plotne(i2,i1) = MAX(1.0E+13,plotne(i2,i1))

              IF     (iopt.EQ.1) THEN
                pop = (GetEAD(plotte(i2,i1),plotne(i2,i1),20,'H.12') /
     .                 GetEAD(plotte(i2,i1),plotne(i2,i1),21,'H.12')) *
     .                (FAC52/FAC32)
              ELSEIF (iopt.EQ.2) THEN
                pop = (GetEAD(plotte(i2,i1),plotne(i2,i1),20,'H.12') /
     .                 GetEAD(plotte(i2,i1),plotne(i2,i1),26,'H.12')) *
     .                (FAC52/FAC42)
c                pop = GetEAD(plotte(i2,i1),plotne(i2,i1),26,'H.12')
              ELSEIF (iopt.EQ.3) THEN
                pop = (GetEAD(plotte(i2,i1),plotne(i2,i1),20,'H.12') /
     .                 GetEAD(plotte(i2,i1),plotne(i2,i1),21,'H.12') /
     .                 GetEAD(plotte(i2,i1),plotne(i2,i1),26,'H.12')) *
     .                (FAC52/FAC32/FAC42)
              ELSEIF (iopt.EQ.4) THEN
                pop = (GetEAD(plotte(i2,i1),plotne(i2,i1),26,'H.12') /
     .                 GetEAD(plotte(i2,i1),plotne(i2,i1),21,'H.12')) *
     .                (FAC42/FAC32)
              ELSEIF (iopt.EQ.5) THEN
                pop = (GetEAD(plotte(i2,i1),plotne(i2,i1),25,'H.12') /
     .                 GetEAD(plotte(i2,i1),plotne(i2,i1),21,'H.12')) *
     .                (FAC21/FAC32)
              ELSEIF (iopt.EQ.6) THEN
                pop = GetEAD(plotte(i2,i1),plotne(i2,i1),20,'H.12')

              ELSE
                WRITE(0,*) 'ERROR 950: UNKNOWN IOPT VALUE',iopt
              ENDIF

              plot_y(i2,i1) = pop
            ENDDO
          ENDDO
c...      Set coloured plot lines:
          DO i1 = 1, nenum
            plottype(i1) = i1 + 1
          ENDDO
c...      Over-write IOPT:
          iopt = 1

        ELSEIF (iclass.eq.1.or.iclass.eq.3.or.
     .          iclass.eq.4.or.iclass.eq.6) THEN

          DO i1 = 1, nenum
            CALL ADASRD(YEAR,cion,nizs,iclass,tenum,
     .                  plotte(1,i1),plotne(1,i1),plot_y(1,i1))
          ENDDO
        ELSEIF (iclass.eq.2.or.iclass.eq.5.or.iclass.eq.7) THEN

          DO i1 = 1, nenum
            CALL ADASRD(YEAR,cion,nizs+1,iclass,tenum,
     .                  plotte(1,i1),plotne(1,i1),plot_y(1,i1))
          ENDDO

        ELSEIF (iclass.EQ.15.OR.iclass.EQ.19.OR.
     .          iclass.EQ.20.OR.iclass.EQ.21) THEN
          call LDADAS2(CION,NIZS,ADASID,ADASYR,ADASEX,ISELE,ISELR,ISELX,
     .                 tenum,plotte,nenum,plotne,plot_mode(iclass),
     .                 plot_y,Wlngth,IRCODE)

          IF (iclass.EQ.15) THEN
            DO i1 = 1, nenum
              CALL ADASRD(YEAR,cion,nizs+1,2,tenum,
     .                    plotte(1,i1),plotne(1,i1),dum1)
              DO i2 = 1, tenum
                plot_y(i2,i1) = dum1(i2) / plot_y(i2,i1)
              ENDDO
            ENDDO
          ENDIF

        ELSEIF (iclass.EQ.16) THEN
          call LDADAS2(CION,NIZS,ADASID,ADASYR,ADASEX,ISELE,ISELR,ISELX,
     .                 tenum,plotte,nenum,plotne,2,
     .                 plot_y,Wlngth,IRCODE)

          DO i1 = 1, nenum
            CALL ADASRD(YEAR,cion,nizs+1,1,tenum,
     .                  plotte(1,i1),plotne(1,i1),dum1)
            DO i2 = 1, tenum
              plot_y(i2,i1) = dum1(i2) / plot_y(i2,i1)
            ENDDO
          ENDDO

        ELSEIF (iclass.EQ.17) THEN
          DO i1 = 1, nenum
            DO i2 = 1, tenum
              plot_y(i2,i1) = 17.5 + (5.0 + 37.5 / plotte(i2,i1)) *
     .                        (1.0 + 0.25 / plotte(i2,i1)) *
     .                        LOG10(1.0E+21 / plotne(i2,i1))
            ENDDO
          ENDDO

        ELSEIF (iclass.EQ.18) THEN
          DO i1 = 1, nenum
            CALL ADASRD(YEAR,cion,nizs+1,5,tenum,
     .                  plotte(1,i1),plotne(1,i1),plot_y(1,i1))
          ENDDO
          DO i1 = 1, nenum
            CALL ADASRD(YEAR,cion,nizs+1,2,tenum,
     .                  plotte(1,i1),plotne(1,i1),dum1)
            DO i2 = 1, tenum
              plot_y(i2,i1) = plot_y(i2,i1)/ dum1(i2) / ECH + 13.6
            ENDDO
          ENDDO

        ELSE
          STOP 'ERROR 950: Invalid ICLASS value'
        ENDIF
c
c
c
        DO i1 = 1, nenum
          WRITE(6,*)
          DO i2 = 1, tenum
            WRITE(6,'(1P,3E15.7)')
     .        plotte(i2,i1),plotne(i2,i1),plot_y(i2,i1)
          ENDDO
        ENDDO
c
c        if (iclass.GE.1.AND.opt_yscale.eq.2) then
c
c          pltmin = HI
c
c          DO i1 = 1, nenum
c            DO i2 = 1, tenum
c              if (plot_y(i2,i1).gt.0.0) then
c                 plot_y(i2,i1) = log10(plot_y(i2,i1))
c                 pltmin = min(pltmin,plot_y(i2,i1))
c              else
c                 plot_y(i2,i1) = -9999.0
c              endif
c            ENDDO
c          enddo
c
c          DO i1 = 1, nenum
c            DO i2 = 1, tenum
c              if (plot_y(i2,i1).eq.-9999.0) then
c                 plot_y(i2,i1) = pltmin -2.0
c              endif
c            ENDDO
c          enddo
c
c          len = lenstr(ylab)  
c          ylab = 'LOG10 ( '//ylab(1:len)//' )'
c          write(nview,'(''MIN ACTVAL= '',g9.3)') pltmin
c
c        endif
c
c        DO i1 = 1, nenum
c          WRITE(6,*)
c          DO i2 = 1, tenum
c            WRITE(6,'(1P,3E15.7)')
c     .        plotte(i2,i1),plotne(i2,i1),plot_y(i2,i1)
c          ENDDO
c        ENDDO
c
c
        CALL DRAW(plotte(1,1),DWIDS,plot_y,MAXGXS,tenum,ANLY,
c     >    nenum,99,0.0,temax,0.0,150.0,IGNORS,ITEC,AVS,NAVS,
     >    nenum,99,temin,temax,-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE2,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,IOPT,2,1.0,0)

c         IPP/09 Krieger - I think colors should be reset here
          call setup_col(ncols,col_opt)
      endif 

c
c ----------------------------------------------------------------------
c
c ----------------------------------------------------------------------
c
c slmod begin - new


      IF (iref.EQ.966) THEN
c...    Multi plots:
        call setup_col(n_cols,5)
        CALL Plot966(nplts,ringnos,graph,nplots,ref,title,iopt,
     >               ngrm,pltmins,pltmaxs,pltfact,iplot,job)
        call setup_col(n_cols,col_opt)

      ELSEIF (IREF.EQ.970) THEN
c...    Iteration history plots:
        call setup_col(n_cols,5)
        CALL Plot970(job,graph,ref,title,iopt)
        call setup_col(n_cols,col_opt)

      ELSEIF (iref.EQ.972) THEN
c...    2D contour plots:
        call setup_col(n_cols,5)
        CALL Plot972(cngs,job,graph,nplots,ref,title,iopt,iplot,
     .               xxmin,xxmax,yymin,yymax,icntr,ft,fp,
     .               xouts,youts,nconts,conts,cntropt,n_cols,col_opt)
        call setup_col(n_cols,col_opt)

      ELSEIF (iref.EQ.974) THEN
c...    2D contour plots:
        call setup_col(n_cols,5)
        CALL Plot974(cngs,job,graph,nplots,ref,title,iopt,iplot,
     .               xxmin,xxmax,yymin,yymax,icntr,ft,fp,
     .               xouts,youts,nconts,conts,cntropt,n_cols,col_opt)
        call setup_col(n_cols,col_opt)

      ELSEIF (iref.EQ.978) THEN
c...    Combo plot:
        call setup_col(n_cols,5)
        CALL Plot978(graph,nplots,ref,title,iopt,pltmins,pltmaxs,iplot,
     .               job)
        call setup_col(n_cols,col_opt)

      ELSEIF (iref.EQ.980) THEN
c...    LOS plots:
c
c       Set up colours for plot980 
c
        call setup_col(n_cols,5)
c
        CALL Plot980(job,graph,ref,title,iopt,
     .               xxmin,xxmax,yymin,yymax,ft,fp,zadj,
     .               ismoth,ignors,itec,avs,navs)
c...    This is here so that 980 is called from 981, and 981
c       may need to draw on the plot after 980 is called:
        CALL FRAME
c
        call setup_col(n_cols,col_opt)
c
      ELSEIF (iref.EQ.981) THEN
c...    Fancy LOS plots:
        call setup_col(n_cols,5)
        CALL Plot981(job,graph,ref,title,iopt,
     .               xxmin,xxmax,yymin,yymax,ft,fp,zadj,
     .               ismoth,ignors,itec,avs,navs)
        call setup_col(n_cols,col_opt)
      ELSEIF (iref.EQ.982) THEN
c...    2D neutral pressure plot:
        call setup_col(n_cols,5)
        CALL Plot982(job,graph,ref,title,iopt,
     .               xxmin,xxmax,yymin,yymax,ft,fp,zadj,
     .               ismoth,ignors,itec,avs,navs,-99.0)
        call setup_col(n_cols,col_opt)
      ELSEIF (iref.EQ.983) THEN
c...    3D geometry plot:
        call setup_col(n_cols,5)
        CALL Plot983(job,graph,ref,title,iopt,
     .               xxmin,xxmax,yymin,yymax,ft,fp,zadj,
     .               ismoth,ignors,itec,avs,navs)
        call setup_col(n_cols,col_opt)
      ELSEIF (iref.EQ.984) THEN
c...    2D MFP plots:
        call setup_col(n_cols,5)
        CALL Plot984(job,graph,ref,title,iopt,
     .               xxmin,xxmax,yymin,yymax,ft,fp,zadj,
     .               ismoth,ignors,itec,avs,navs,-99.0)
        call setup_col(n_cols,col_opt)
      ELSEIF (iref.EQ.985) THEN
c...    New 3D analysis routines:
        call setup_col(n_cols,5)
        CALL Plot985(job,graph,ref,title,iopt,
     .               xxmin,xxmax,yymin,yymax,ft,fp,zadj,
     .               ismoth,ignors,itec,avs,navs)
        call setup_col(n_cols,col_opt)
      ELSEIF (iref.EQ.986) THEN
c...    Multi-case bar and line charts:
        call setup_col(n_cols,5)
        CALL Plot986(graph,nplots,ref,title,iopt,
     >               pltmins,pltmaxs,iplot,job)
        call setup_col(n_cols,col_opt)
      ELSEIF (iref.EQ.987) THEN
c...    New 2D plots (including EIRENE triangle plots):
        call setup_col(n_cols,5)
        CALL Plot987(job,graph,ref,title,iopt,
     .               xxmin,xxmax,yymin,yymax,ft,fp,zadj,
     .               ismoth,ignors,itec,avs,navs,nizs)
        call setup_col(n_cols,col_opt)
      ELSEIF (iref.EQ.988) THEN
c...    Line shapes:
        call setup_col(n_cols,5)
c        CALL Plot988(job,graph,ref,title,iopt,
c     .               xxmin,xxmax,yymin,yymax,ft,fp,zadj,
c     .               ismoth,ignors,itec,avs,navs)
        call setup_col(n_cols,col_opt)
      ELSEIF (iref.EQ.989) THEN
c...    Tangentially viewing camera image processing:
        call setup_col(n_cols,5)
        CALL Plot989(job,graph,ref,title,iopt,
     .               xxmin,xxmax,yymin,yymax,ft,fp,zadj,
     .               ismoth,ignors,itec,avs,navs)
        call setup_col(n_cols,col_opt)
      ELSEIF (iref.EQ.990) THEN
c...    Thomson shift analysis plot:
        call setup_col(n_cols,5)
        CALL Plot990(nplts,ringnos,graph,nplots,ref,title,iopt,
     .               ngrm,pltmins,pltmaxs,pltfact,iplot)

        call setup_col(n_cols,col_opt)

      ELSEIF (iref.EQ.998) THEN
c...    Dump EIRENE reaction data to file:
        CALL DumpRates(iopt)

      ELSEIF (iref.EQ.999) THEN
        call setup_col(n_cols,5)
        CALL Development(iopt,nizs,cizsc,crmi,cion,absfac,title)
        call setup_col(n_cols,col_opt)
c slmod end

      ENDIF

      return 

 9999 CONTINUE
      WRITE(0,*) 'MAXGXS=',MAXGXS     
      ierr= -99 

      return

c
c     Format statements    
c
 9012 FORMAT(1X,'PLOT',I3,4X,A)
 

      end

