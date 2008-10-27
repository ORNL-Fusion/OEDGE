      subroutine out700(iref,graph,iopt,ierr)
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
      include 'pindata'
c      include 'cadas'
c      include 'grbound'
c      include 'outxy'
      include 'cedge2d'
      include 'transcoef'
c      include 'cioniz'
c      include 'reiser' 
c      include 'printopt' 
c
c     Local Variables
c
      integer ik,ir,iz,it
      integer in
c
c     Local Variables
c


      integer inc,irlim

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

      integer ikstart,ikend,nstep
      real dist


c
        IF (IOPT.EQ.0) RETURN
c
        if (iref.ge.701.and.iref.le.750) then

c         Series 700 Plots - for the multiple plots on one page
c         read in a data line specifying:
c
c         Number of plots/page
c         Total number of plots
c         List of rings to be plotted - Listed are a number of
c         rings equal to the total number to be plotted.
c
c
           call rdg4(graph4,ngrm,nplts,ringnos,maxplts,
     >               pltfact,ierr)
c
           IF (IERR.EQ.1) THEN
              WRITE(6,*) 'RDG4 ERROR READING 700 SERIES- GRAPH DETAILS'
              IERR = 0
              RETURN
           ENDIF
c
c          Make up plot labels for each ring
c
           do ip = 1,nplts
             write (pltlabs(ip),'(''Ring:'',i2)') ringnos(ip)
           end do
c
        endif


      call init_plot(iref,graph,iopt) 


c
C-----------------------------------------------------------------------
C
C     700 Series Plots - Plots with multiple graphs on each page.
c     701 703 705 709 711 713 - Ne Te Ti Vb Gamma and Ionization
c                               cross-field gradn, gradTe, gradTi
c                               for DIVIMP and EDGE2D
c
c     721 - Pressure - actual and expected
C
C-----------------------------------------------------------------------
C
c
c     701 - Density
c
      if (iref.eq.701) then
c
c        Plots Ne for E2D and DIVIMP for 8 rings on one
c        sheet of paper for easier comparison. It is set up
c        so that the plots can be modified for either vs. S or vs. P
c
        REF = 'COMPARATIVE DENSITY PLOTS: FC vs DIV'
c
        ELABS(1) = 'NBFCNBFC    '
        ELABS(2) = 'NBD NBD    '
c
c       Set up the ring numbers for the plots
c
c       Set up for 2 sets of data on each of 8 plots
c
        ngs = 2
        sctype = iopt
c
        IF (IREF.EQ.701) THEN
          XLAB = '   S  (M)'
          axistype = 1
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
          axistype = 2
        ENDIF
c
        YLAB   = 'Density'
c
        NPLOTS = NPLOTS + 1
c        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)

        WRITE (IPLOT,9012) NPLOTS,REF

        CALL rzero (mvals,maxnks*maxngs*maxplts)
c
c
        do ip = 1, nplts
c
c         Access ring number and store number of knots info.
c
          ir = ringnos(ip)
c
c         Modify the following to include values for S=0 ... even if this
c         is not part of the grid. Assume that if one end is not a grid
c         point then neither is the other.
c
          if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
             in = 0
             inc = 0
          else
             in = 1
             inc = 2
          endif
c
c         Set number of points on plots
c
          pnks(ip,1) = nks(ir)+inc
c
c         Load axis
c
          call loadm_axis(mouts,mwids,ir,ip,axistype,in)
c
c         Load data
c
c         Load target data for rings that need it.
c
          if (in.ne.0) then
c
            mVALS(1,ip,1) = e2dnbs(1,ir)
            mVALS(1,ip,2) = KNdS (idds(ir,2))

            mVALS(nks(ir)+2,ip,1) = e2dnbs(nks(ir),ir)
            mVALS(nks(ir)+2,ip,2) = KNdS (idds(ir,1))
c
          endif
c
          DO IK = 1, NKS(IR)
c
            MVALS(IK+in,ip,1) = e2dnbs(ik,ir)
            MVALS(IK+in,ip,2) = KNBS (IK,IR)
c
          enddo
C
        enddo
c
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
      endif
c
c     703 - Te
c
      if (iref.eq.703) then
c
c        Plots Te for E2D and DIVIMP for 8 rings on one
c        sheet of paper for easier comparison. It is set up
c        so that the plots can be modified for either vs. S or vs. P
c
        REF = 'COMPARATIVE Te PLOTS: FC vs. DIV'
c
        ELABS(1) = 'TeFCTeFC   '
        ELABS(2) = 'TeD TeD    '
c
c       Set up the ring numbers for the plots
c
c       Set up for 2 sets of data on each of 8 plots
c
        ngs = 2
        sctype = iopt
c
        IF (IREF.EQ.703) THEN
          XLAB = '   S  (M)'
          axistype = 1
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
          axistype = 2
        ENDIF
c
        YLAB   = 'Te'
c
        NPLOTS = NPLOTS + 1
c        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)

        WRITE (IPLOT,9012) NPLOTS,REF

        CALL rzero (mvals,maxnks*maxngs*maxplts)
c
c
        do ip = 1, nplts
c
c         Access ring number and store number of knots info.
c
          ir = ringnos(ip)
c
c         Modify the following to include values for S=0 ... even if this
c         is not part of the grid. Assume that if one end is not a grid
c         point then neither is the other.
c
          if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
             in = 0
             inc = 0
          else
             in = 1
             inc = 2
          endif
c
c         Set number of points on plots
c
          pnks(ip,1) = nks(ir)+inc
c
c         Load axis
c
          call loadm_axis(mouts,mwids,ir,ip,axistype,in)
c
c         Load data
c
c         Load target data for rings that need it.
c
          if (in.ne.0) then
c
             mVALS(1,ip,1) = e2dtarg(ir,2,2)
             mVALS(1,ip,2) = KtedS (idds(ir,2))
c
             mVALS(nks(ir)+2,ip,1) = e2dtarg(ir,2,1)
             mVALS(nks(ir)+2,ip,2) = KtedS (idds(ir,1))
c
          endif
c
c         Loop through ring and load the rest of the data
c
          DO IK = 1, NKS(IR)
c
            MVALS(IK+in,ip,1) = e2dtebs(ik,ir)
            MVALS(IK+in,ip,2) = KteBS (IK,IR)
c
          enddo
C
        enddo
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
      endif
c
c     705 - Ti
c
      if (iref.eq.705) then
c
c        Plots Ti for E2D and DIVIMP for 8 rings on one
c        sheet of paper for easier comparison. It is set up
c        so that the plots can be modified for either vs. S or vs. P
c
        REF = 'COMPARATIVE Ti PLOTS: FC vs. DIV'
c
        ELABS(1) = 'TiFCTiFC   '
        ELABS(2) = 'TiD TiD    '
c
c       Set up the ring numbers for the plots
c
c       Set up for 2 sets of data on each of 8 plots
c
        ngs = 2
        sctype = iopt
c
        IF (IREF.EQ.705) THEN
          XLAB = '   S  (M)'
          axistype = 1
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
          axistype = 2
        ENDIF
c
        YLAB   = 'Ti'
c
        NPLOTS = NPLOTS + 1
c        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)

        WRITE (IPLOT,9012) NPLOTS,REF

        CALL rzero (mvals,maxnks*maxngs*maxplts)
c
c
        do ip = 1, nplts
c
c         Access ring number and store number of knots info.
c
          ir = ringnos(ip)
c
c         Modify the following to include values for S=0 ... even if this
c         is not part of the grid. Assume that if one end is not a grid
c         point then neither is the other.
c
          if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
             in = 0
             inc = 0
          else
             in = 1
             inc = 2
          endif
c
c         Set number of points on plots
c
          pnks(ip,1) = nks(ir)+inc
c
c         Load axis
c
          call loadm_axis(mouts,mwids,ir,ip,axistype,in)
c
c         Load data
c
c         Load target data for rings that need it.
c
          if (in.ne.0) then
c
            mVALS(1,ip,1) = e2dtarg(ir,3,2)
            mVALS(1,ip,2) = KtidS (idds(ir,2))

            mVALS(nks(ir)+2,ip,1) = e2dtarg(ir,3,1)
            mVALS(nks(ir)+2,ip,2) = KtidS (idds(ir,1))
c
          endif
c
          DO IK = 1, NKS(IR)
c
            MVALS(IK+in,ip,1) = e2dtibs(ik,ir)
            MVALS(IK+in,ip,2) = KtiBS (IK,IR)
c
          enddo
C
        enddo
c
c        enldist = mgst*kouts(nks(ir)+inc)
c
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
      endif
c
c     707 - Background Velocity
c
      if (iref.eq.707) then
c
c        Plots Vb for E2D and DIVIMP for 8 rings on one
c        sheet of paper for easier comparison. It is set up
c        so that the plots can be modified for either vs. S or vs. P
c
        REF = 'COMPARATIVE VELOCITY PLOTS: FC vs. DIV'
c
        ELABS(1) = 'VbFCVbFC    '
        ELABS(2) = 'VbD VbD    '
c
c       Set up the ring numbers for the plots
c
c       Set up for 2 sets of data on each of 8 plots
c
        ngs = 2
        sctype = iopt
c
        IF (IREF.EQ.707) THEN
          XLAB = '   S  (M)'
          axistype = 1
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
          axistype = 2
        ENDIF
c
        YLAB   = 'Velocity'
c
        NPLOTS = NPLOTS + 1
c        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)

        WRITE (IPLOT,9012) NPLOTS,REF

        CALL rzero (mvals,maxnks*maxngs*maxplts)
c
c
        do ip = 1, nplts
c
c         Access ring number and store number of knots info.
c
          ir = ringnos(ip)
c
c         Modify the following to include values for S=0 ... even if this
c         is not part of the grid. Assume that if one end is not a grid
c         point then neither is the other.
c
          if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
             in = 0
             inc = 0
          else
             in = 1
             inc = 2
          endif
c
c         Set number of points on plots
c
          pnks(ip,1) = nks(ir)+inc
c
c         Load axis
c
          call loadm_axis(mouts,mwids,ir,ip,axistype,in)
c
c         Load data
c
c         Load target data for rings that need it.
c
          if (in.ne.0) then
c
            mVALS(1,ip,1) = e2dtarg(ir,4,2)
            mVALS(1,ip,2) = KvdS (idds(ir,2))

            mVALS(nks(ir)+2,ip,1) = e2dtarg(ir,4,1)
            mVALS(nks(ir)+2,ip,2) = KvdS (idds(ir,1))
c
          endif
c
c          write (6,*) 'Vels1:', kvds(idds(ir,2)),kvds(idds(ir,1)),
c     >                             idds(ir,2),idds(ir,1)
c
          DO IK = 1, NKS(IR)
c
            MVALS(IK+in,ip,1) = e2dvhs(ik,ir)
            MVALS(IK+in,ip,2) = KvhS (IK,IR) / qtim
c
c            write (6,*) 'vels:',ik,ir,qtim,kvhs(ik,ir),kvhs(ik,ir)/qtim
c
          enddo
C
        enddo
c
c        enldist = mgst*kouts(nks(ir)+inc)
c
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
      endif
c
c     709 - Gamma
c
      if (iref.eq.709) then
c
c        Plots Ne for E2D and DIVIMP for 8 rings on one
c        sheet of paper for easier comparison. It is set up
c        so that the plots can be modified for either vs. S or vs. P
c
        REF = 'COMPARATIVE Gamma PLOTS: FC vs. DIV'
c
        ELABS(1) = 'GaFCGaFC   '
        ELABS(2) = 'GaD GaD    '
        ELABS(3) = 'GaEDGaED   '
c
c       Set up the ring numbers for the plots
c
c       Set up for 2 sets of data on each of 8 plots
c
        ngs = 3
        sctype = iopt
c
        IF (IREF.EQ.709) THEN
          XLAB = '   S  (M)'
          axistype = 1
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
          axistype = 2
        ENDIF
c
        YLAB   = 'Gamma (Flux)'
c
        NPLOTS = NPLOTS + 1
c        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)

        WRITE (IPLOT,9012) NPLOTS,REF

        CALL rzero (mvals,maxnks*maxngs*maxplts)
c
c
        do ip = 1, nplts
c
c         Access ring number and store number of knots info.
c
          ir = ringnos(ip)
c
c         Modify the following to include values for S=0 ... even if this
c         is not part of the grid. Assume that if one end is not a grid
c         point then neither is the other.
c
          if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
             in = 0
             inc = 0
          else
             in = 1
             inc = 2
          endif
c
c         Set number of points on plots
c
          pnks(ip,1) = nks(ir)+inc
c
c         Load axis
c
          call loadm_axis(mouts,mwids,ir,ip,axistype,in)
c
c         Load data
c
c         Load target data for rings that need it.
c
          if (in.ne.0) then
c
            mVALS(1,ip,1) = e2dtarg(ir,5,2)
            mVALS(1,ip,2) = KNdS(idds(ir,2))*KvdS(idds(ir,2))
            mVALS(1,ip,3) = e2dgpara(1,ir)

            mVALS(nks(ir)+2,ip,1) = e2dtarg(ir,5,1)
            mVALS(nks(ir)+2,ip,2) = KNdS(idds(ir,1))*KvdS(idds(ir,1))
            mVALS(nks(ir)+2,ip,3) = e2dgpara(nks(ir)+1,ir)
c
          endif
c
          DO IK = 1, NKS(IR)
c
            MVALS(IK+in,ip,1) = e2dnbs(ik,ir) * e2dvhs(ik,ir)
            MVALS(IK+in,ip,2) = KNBS (IK,IR)* kvhs(ik,ir)/qtim
            MVALS(IK+in,ip,3) = (e2dgpara(ik,ir)+e2dgpara(ik+1,ir))/2.0
c
          enddo
C
        enddo
c
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
      endif
c
c     711 - Ionization
c
      if (iref.eq.711) then
c
c        Plots Ne for E2D and DIVIMP for 8 rings on one
c        sheet of paper for easier comparison. It is set up
c        so that the plots can be modified for either vs. S or vs. P
c
        REF = 'COMPARATIVE IONIZ PLOTS: FC vs. DIV'
c
        ELABS(1) = 'HzFZHzFC    '
        ELABS(2) = 'HizDHizD    '
c
c       Set up the ring numbers for the plots
c
c       Set up for 2 sets of data on each of 8 plots
c
        ngs = 2
        sctype = iopt
c
        IF (IREF.EQ.711) THEN
          XLAB = '   S  (M)'
          axistype = 1
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
          axistype = 2
        ENDIF
c
        YLAB   = 'H Ionization'
c
        NPLOTS = NPLOTS + 1
c        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)

        WRITE (IPLOT,9012) NPLOTS,REF

        CALL rzero (mvals,maxnks*maxngs*maxplts)
c
c
        do ip = 1, nplts
c
c         Access ring number and store number of knots info.
c
          ir = ringnos(ip)
c
c         Set number of points on plots
c
          pnks(ip,1) = nks(ir)
c
c         Load axis
c
          call loadm_axis(mouts,mwids,ir,ip,axistype,0)
c
c         Load data
c
          DO IK = 1, NKS(IR)
c
            MVALS(IK+in,ip,1) = e2dion(ik,ir)
            MVALS(IK+in,ip,2) = pinion(IK,IR)
c
          enddo
C
        enddo
c
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
      endif
c
c     713 - Gradn
c
      if (iref.eq.713) then
c
c        Plots gradNe for E2D and DIVIMP for 8 rings on one
c        sheet of paper for easier comparison. It is set up
c        so that the plots can be modified for either vs. S or vs. P
c
        REF = 'COMPARATIVE GRADN PLOTS: FC vs. DIV'
c
        ELABS(1) = 'GRnEGRnE    '
        ELABS(2) = 'GRnDGRnD    '
c
c       Set up the ring numbers for the plots
c
c       Set up for 2 sets of data on each of 8 plots
c
        ngs = 2
        sctype = iopt
c
        IF (IREF.EQ.713) THEN
          XLAB = '   S  (M)'
          axistype = 1
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
          axistype = 2
        ENDIF
c
        YLAB   = 'Ne Gradients'
c
        NPLOTS = NPLOTS + 1
c        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)

        WRITE (IPLOT,9012) NPLOTS,REF

        CALL rzero (mvals,maxnks*maxngs*maxplts)
c
c
        do ip = 1, nplts
c
c         Access ring number and store number of knots info.
c
          ir = ringnos(ip)
c
c         Set number of points on plots
c
          pnks(ip,1) = nks(ir)
c
c         Load axis
c
          call loadm_axis(mouts,mwids,ir,ip,axistype,0)
c
c         Load data
c
          DO IK = 1, NKS(IR)
c
            MVALS(IK+in,ip,1) = e2dgradn(ik,ir)
            MVALS(IK+in,ip,2) = gradn(IK,IR)
c
          enddo
C
        enddo
c
c        enldist = mgst*kouts(nks(ir)+inc)
c
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
      endif
c
c     715 - Gradte
c
      if (iref.eq.715) then
c
c        Plots Gradte for E2D and DIVIMP for 8 rings on one
c        sheet of paper for easier comparison. It is set up
c        so that the plots can be modified for either vs. S or vs. P
c
        REF = 'COMPARATIVE GRAD-Te PLOTS: FC vs. DIV'
c
        ELABS(1) = 'GTeEGTeE    '
        ELABS(2) = 'GTeDGTeD    '
c
c       Set up the ring numbers for the plots
c
c       Set up for 2 sets of data on each of 8 plots
c
        ngs = 2
        sctype = iopt
c
        IF (IREF.EQ.715) THEN
          XLAB = '   S  (M)'
          axistype = 1
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
          axistype = 2
        ENDIF
c
        YLAB   = 'Te Gradients'
c
        NPLOTS = NPLOTS + 1
c        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)

        WRITE (IPLOT,9012) NPLOTS,REF

        CALL rzero (mvals,maxnks*maxngs*maxplts)
c
c
        do ip = 1, nplts
c
c         Access ring number and store number of knots info.
c
          ir = ringnos(ip)
c
c         Set number of points on plots
c
          pnks(ip,1) = nks(ir)
c
c         Load axis
c
          call loadm_axis(mouts,mwids,ir,ip,axistype,0)
c
c         Load data
c
          DO IK = 1, NKS(IR)
c
            MVALS(IK+in,ip,1) = e2dgradte(ik,ir)
            MVALS(IK+in,ip,2) = gradte(IK,IR)
c
          enddo
C
        enddo
c
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
      endif
c
c     717 - GradTi
c
      if (iref.eq.717) then
c
c        Plots GradTi for E2D and DIVIMP for 8 rings on one
c        sheet of paper for easier comparison. It is set up
c        so that the plots can be modified for either vs. S or vs. P
c
        REF = 'COMPARATIVE GRAD-Ti PLOTS: FC vs. DIV'
c
        ELABS(1) = 'GTiEGTiE    '
        ELABS(2) = 'GTiDGTiD    '
c
c       Set up the ring numbers for the plots
c
c       Set up for 2 sets of data on each of 8 plots
c
        ngs = 2
        sctype = iopt
c
        IF (IREF.EQ.717) THEN
          XLAB = '   S  (M)'
          axistype = 1
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
          axistype = 2
        ENDIF
c
        YLAB   = 'Ti Gradients'
c
        NPLOTS = NPLOTS + 1
c        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)

        WRITE (IPLOT,9012) NPLOTS,REF

        CALL rzero (mvals,maxnks*maxngs*maxplts)
c
c
        do ip = 1, nplts
c
c         Access ring number and store number of knots info.
c
          ir = ringnos(ip)
c
c         Set number of points on plots
c
          pnks(ip,1) = nks(ir)
c
c         Load axis
c
          call loadm_axis(mouts,mwids,ir,ip,axistype,0)
c
c         Load data
c
          DO IK = 1, NKS(IR)
c
            MVALS(IK+in,ip,1) = e2dgradti(ik,ir)
            MVALS(IK+in,ip,2) = gradti(IK,IR)
c
          enddo
C
        enddo
c
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
      endif
c
c     Comparative Gamma plots
c
c
c     719 - Gamma Plots - EDGE2D comparisons
c
      if (iref.eq.719) then
c
c       Comparisons of Egde2D Fluxes calculated in different ways
c
c
        XLAB = '   S  (M)'
        YLAB   = 'Gamma (Flux)'
c
c
        ELABS(1) = 'GaE GaE    '
        ELABS(2) = 'GaC GaC    '
        ELABS(3) = 'GaCRGaCRcor'
c
c       Set up the ring numbers for the plots
c
c       Set up for 2 sets of data on each of 8 plots
c
        ngs = 3
        sctype = iopt
c
c       Plot along S
c
        axistype = 1
c
c       ***********  EDGE2D 0  *************************
c
c       The following code is repeated for each of the 6
c       types of Gamma plots.
c
        REF = 'Fluxes : E2D vs. CALC:EDGE2D F'
c
        NPLOTS = NPLOTS + 1
c
        WRITE (IPLOT,9012) NPLOTS,REF
c
        CALL rzero (mvals,maxnks*maxngs*maxplts)
c
c
        do ip = 1, nplts
c
c         Access ring number and store number of knots info.
c
          ir = ringnos(ip)
c
c         Modify the following to include values for S=0 ... even if this
c         is not part of the grid. Assume that if one end is not a grid
c         point then neither is the other.
c
          if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
             in = 0
             inc = 0
          else
             in = 1
             inc = 2
          endif
c
c         Set number of points on plots
c
          pnks(ip,1) = nks(ir)+inc
c
c         Load axis
c
          call loadm_axis(mouts,mwids,ir,ip,axistype,in)
c
c         Load data
c
c         Load target data for rings that need it.
c
          if (in.ne.0) then
c
            mVALS(1,ip,1) = e2dtarg(ir,5,2)
            mVALS(1,ip,2) = e2dtarg(ir,5,2)

            mVALS(nks(ir)+2,ip,1) = e2dtarg(ir,5,1)
            mVALS(nks(ir)+2,ip,2) = e2dtarg(ir,5,1)
c
          endif
c
          DO IK = 1, NKS(IR)
c
            MVALS(IK+in,ip,1) = e2dnbs(ik,ir) * e2dvhs(ik,ir)
            MVALS(IK+in,ip,2) = fluxes(ik,ir,9)
            MVALS(IK+in,ip,3) = fluxes(ik,ir,10)
c
          enddo
C
        enddo
c
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
c       ***********  EDGE2D 1/2  *************************
c
c       The following code is repeated for each of the 6
c       types of Gamma plots.
c
        REF = 'Fluxes : E2D vs. CALC:EDGE2D 1/2'
c
        NPLOTS = NPLOTS + 1
c
        WRITE (IPLOT,9012) NPLOTS,REF
c
        CALL rzero (mvals,maxnks*maxngs*maxplts)
c
c
        do ip = 1, nplts
c
c         Access ring number and store number of knots info.
c
          ir = ringnos(ip)
c
c         Modify the following to include values for S=0 ... even if this
c         is not part of the grid. Assume that if one end is not a grid
c         point then neither is the other.
c
          if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
             in = 0
             inc = 0
          else
             in = 1
             inc = 2
          endif
c
c         Set number of points on plots
c
          pnks(ip,1) = nks(ir)+inc
c
c         Load axis
c
          call loadm_axis(mouts,mwids,ir,ip,axistype,in)
c
c         Load data
c
c         Load target data for rings that need it.
c
          if (in.ne.0) then
c
            mVALS(1,ip,1) = e2dtarg(ir,5,2)
            mVALS(1,ip,2) = e2dtarg(ir,5,2)

            mVALS(nks(ir)+2,ip,1) = e2dtarg(ir,5,1)
            mVALS(nks(ir)+2,ip,2) = e2dtarg(ir,5,1)
c
          endif
c
          DO IK = 1, NKS(IR)
c
            MVALS(IK+in,ip,1) = e2dnbs(ik,ir) * e2dvhs(ik,ir)
            MVALS(IK+in,ip,2) = fluxes(ik,ir,11)
            MVALS(IK+in,ip,3) = fluxes(ik,ir,12)
c
          enddo
C
        enddo
c
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
c       ***********  EDGE2D POL F *************************
c
c       The following code is repeated for each of the 6
c       types of Gamma plots.
c
        REF = 'Fluxes : E2D vs. CALC:EDGE2D POL F'
c
        NPLOTS = NPLOTS + 1
c
        WRITE (IPLOT,9012) NPLOTS,REF
c
        CALL rzero (mvals,maxnks*maxngs*maxplts)
c
c
        do ip = 1, nplts
c
c         Access ring number and store number of knots info.
c
          ir = ringnos(ip)
c
c         Modify the following to include values for S=0 ... even if this
c         is not part of the grid. Assume that if one end is not a grid
c         point then neither is the other.
c
          if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
             in = 0
             inc = 0
          else
             in = 1
             inc = 2
          endif
c
c         Set number of points on plots
c
          pnks(ip,1) = nks(ir)+inc
c
c         Load axis
c
          call loadm_axis(mouts,mwids,ir,ip,axistype,in)
c
c         Load data
c
c         Load target data for rings that need it.
c
          if (in.ne.0) then
c
            mVALS(1,ip,1) = e2dtarg(ir,5,2)
            mVALS(1,ip,2) = e2dtarg(ir,5,2)

            mVALS(nks(ir)+2,ip,1) = e2dtarg(ir,5,1)
            mVALS(nks(ir)+2,ip,2) = e2dtarg(ir,5,1)
c
          endif
c
          DO IK = 1, NKS(IR)
c
            MVALS(IK+in,ip,1) = e2dnbs(ik,ir) * e2dvhs(ik,ir)
            MVALS(IK+in,ip,2) = fluxes(ik,ir,1)
            MVALS(IK+in,ip,3) = fluxes(ik,ir,2)
c
          enddo
C
        enddo
c
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
c       ***********  EDGE2D POL 1/2  *************************
c
c       The following code is repeated for each of the 6
c       types of Gamma plots.
c
        REF = 'Fluxes : E2D vs. CALC:EDGE2D POL 1/2'
c
        NPLOTS = NPLOTS + 1
c
        WRITE (IPLOT,9012) NPLOTS,REF
c
        CALL rzero (mvals,maxnks*maxngs*maxplts)
c
c
        do ip = 1, nplts
c
c         Access ring number and store number of knots info.
c
          ir = ringnos(ip)
c
c         Modify the following to include values for S=0 ... even if this
c         is not part of the grid. Assume that if one end is not a grid
c         point then neither is the other.
c
          if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
             in = 0
             inc = 0
          else
             in = 1
             inc = 2
          endif
c
c         Set number of points on plots
c
          pnks(ip,1) = nks(ir)+inc
c
c         Load axis
c
          call loadm_axis(mouts,mwids,ir,ip,axistype,in)
c
c         Load data
c
c         Load target data for rings that need it.
c
          if (in.ne.0) then
c
            mVALS(1,ip,1) = e2dtarg(ir,5,2)
            mVALS(1,ip,2) = e2dtarg(ir,5,2)

            mVALS(nks(ir)+2,ip,1) = e2dtarg(ir,5,1)
            mVALS(nks(ir)+2,ip,2) = e2dtarg(ir,5,1)
c
          endif
c
          DO IK = 1, NKS(IR)
c
            MVALS(IK+in,ip,1) = e2dnbs(ik,ir) * e2dvhs(ik,ir)
            MVALS(IK+in,ip,2) = fluxes(ik,ir,3)
            MVALS(IK+in,ip,3) = fluxes(ik,ir,4)
c
          enddo
C
        enddo
c
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
c       ***********  DIVIMP F  *************************
c
c       The following code is repeated for each of the 6
c       types of Gamma plots.
c
        REF = 'Fluxes : E2D vs. CALC:DIV/PIN F'
c
        NPLOTS = NPLOTS + 1
c
        WRITE (IPLOT,9012) NPLOTS,REF
c
        CALL rzero (mvals,maxnks*maxngs*maxplts)
c
c
        do ip = 1, nplts
c
c         Access ring number and store number of knots info.
c
          ir = ringnos(ip)
c
c         Modify the following to include values for S=0 ... even if this
c         is not part of the grid. Assume that if one end is not a grid
c         point then neither is the other.
c
          if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
             in = 0
             inc = 0
          else
             in = 1
             inc = 2
          endif
c
c         Set number of points on plots
c
          pnks(ip,1) = nks(ir)+inc
c
c         Load axis
c
          call loadm_axis(mouts,mwids,ir,ip,axistype,in)
c
c         Load data
c
c         Load target data for rings that need it.
c
          if (in.ne.0) then
c
            mVALS(1,ip,1) = e2dtarg(ir,5,2)
            mVALS(1,ip,2) = KNdS(idds(ir,1))*KvdS(idds(ir,1))

            mVALS(nks(ir)+2,ip,1) = e2dtarg(ir,5,1)
            mVALS(nks(ir)+2,ip,2) = KNdS(idds(ir,1))*KvdS(idds(ir,1))
c
          endif
c
          DO IK = 1, NKS(IR)
c
            MVALS(IK+in,ip,1) = e2dnbs(ik,ir) * e2dvhs(ik,ir)
            MVALS(IK+in,ip,2) = fluxes(ik,ir,13)
            MVALS(IK+in,ip,3) = fluxes(ik,ir,14)
c
          enddo
C
        enddo
c
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
c       ***********  DIV/PIN 1/2  *************************
c
c       The following code is repeated for each of the 6
c       types of Gamma plots.
c
        REF = 'Fluxes : E2D vs. CALC:DIV/PIN 1/2'
c
        NPLOTS = NPLOTS + 1
c
        WRITE (IPLOT,9012) NPLOTS,REF
c
        CALL rzero (mvals,maxnks*maxngs*maxplts)
c
c
        do ip = 1, nplts
c
c         Access ring number and store number of knots info.
c
          ir = ringnos(ip)
c
c         Modify the following to include values for S=0 ... even if this
c         is not part of the grid. Assume that if one end is not a grid
c         point then neither is the other.
c
          if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
             in = 0
             inc = 0
          else
             in = 1
             inc = 2
          endif
c
c         Set number of points on plots
c
          pnks(ip,1) = nks(ir)+inc
c
c         Load axis
c
          call loadm_axis(mouts,mwids,ir,ip,axistype,in)
c
c         Load data
c
c         Load target data for rings that need it.
c
          if (in.ne.0) then
c
            mVALS(1,ip,1) = e2dtarg(ir,5,2)
            mVALS(1,ip,2) = KNdS(idds(ir,1))*KvdS(idds(ir,1))

            mVALS(nks(ir)+2,ip,1) = e2dtarg(ir,5,1)
            mVALS(nks(ir)+2,ip,2) = KNdS(idds(ir,1))*KvdS(idds(ir,1))
c
          endif
c
          DO IK = 1, NKS(IR)
c
            MVALS(IK+in,ip,1) = e2dnbs(ik,ir) * e2dvhs(ik,ir)
            MVALS(IK+in,ip,2) = fluxes(ik,ir,15)
            MVALS(IK+in,ip,3) = fluxes(ik,ir,16)
c
          enddo
C
        enddo
c
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
      endif
c
c
c     721 - Pressure - actual and expected from DIVIMP
c
c
      if (iref.eq.721) then
c
c        Plots Ne for E2D and DIVIMP for 8 rings on one
c        sheet of paper for easier comparison. It is set up
c        so that the plots can be modified for either vs. S or vs. P
c
        REF = 'COMPARATIVE PRESSURE PLOTS: FC vs. DIV'
c
        ELABS(1) = 'PrFCPrFC   '
        ELABS(2) = 'PACTPACT   '
        ELABS(3) = 'PEXPPEXP   '
c
c       Set up the ring numbers for the plots
c
c       Set up for 2 sets of data on each of 8 plots
c
        ngs = 2
        sctype = iopt
c
        IF (IREF.EQ.721) THEN
          XLAB = '   S  (M)'
          axistype = 1
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
          axistype = 2
        ENDIF
c
        YLAB   = 'Pressure'
c
        NPLOTS = NPLOTS + 1
c        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)

        WRITE (IPLOT,9012) NPLOTS,REF

        CALL rzero (mvals,maxnks*maxngs*maxplts)
c
c
        do ip = 1, nplts
c
c         Access ring number and store number of knots info.
c
          ir = ringnos(ip)
c
c         Modify the following to include values for S=0 ... even if this
c         is not part of the grid. Assume that if one end is not a grid
c         point then neither is the other.
c
          if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
             in = 0
             inc = 0
          else
             in = 1
             inc = 2
          endif
c
c         Set number of points on plots
c
          pnks(ip,1) = nks(ir)+inc
c
c         Load axis
c
          call loadm_axis(mouts,mwids,ir,ip,axistype,in)
c
c         Load data
c
c         Load target data for rings that need it.
c
          if (in.ne.0) then
c
c           Load actual target pressures.
c
            mVALS(1,ip,1) = knds(idds(ir,2))* ((kteds(idds(ir,2))
     >                             +ktids(idds(ir,2))) * ech +
     >             crmb * amu * kvds(idds(ir,2))**2)
            mVALS(1,ip,2) = mvals(1,ip,1)
c
c            mVALS(1,ip,3) = mvals(1,ip,1)
c
            mVALS(nks(ir)+2,ip,1) = knds(idds(ir,1))
     >                          *(( kteds(idds(ir,1))
     >                             +ktids(idds(ir,1))) * ech +
     >             crmb * amu * kvds(idds(ir,1))**2 )

            mVALS(nks(ir)+2,ip,2) =  mVALS(nks(ir)+2,ip,1)
c
c            mVALS(nks(ir)+2,ip,3) =  mVALS(nks(ir)+2,ip,1)
c
          endif
c
          DO IK = 1, NKS(IR)
c
c           E2D
c
            MVALS(ik+in,ip,1) = e2dnbs(ik,ir)*( ((e2dtebs(ik,ir)
     >                             +e2dtibs(ik,ir))) * ech +
     >             crmb * amu * e2dvhs(ik,ir)**2)
c
c           Actually used
c
            MVALS(IK+in,ip,2) = knbs(ik,ir)*( ((ktebs(ik,ir)
     >                             +ktibs(ik,ir))) * ech +
     >             crmb * amu * (kvhs(ik,ir)/qtim)**2)
c
c           Expected
c
c            MVALS(IK+in,ip,1) = kpress(ik,ir,1)
c

c
            write(6,'(a,2i4,8(1x,g12.5))')  'COMP:',ik,ir,
     >                  knbs(ik,ir),e2dnbs(ik,ir),
     >                  kvhs(ik,ir),e2dvhs(ik,ir),
     >                  ktibs(ik,ir),e2dtibs(ik,ir),
     >                  ktebs(ik,ir),e2dtebs(ik,ir)
c
c
          enddo
C
        enddo
c
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
      endif
c
c
c     723 - Ionization and Recombination
c
      if (iref.eq.723) then
c
c       Plots of ionization and recombination
c
c
        REF = 'IONIZATION/RECOMBINATION PLOTS - PIN'
c
        ELABS(1) = 'Hiz Hiz      '
        ELABS(2) = 'HrecHrec     '
c
c       Set up the ring numbers for the plots
c
c       Set up for 2 sets of data on each of 8 plots
c
        ngs = 2
        sctype = iopt
        if (sctype.lt.1.or.sctype.gt.6) sctype = 1
c
        IF (IREF.EQ.723) THEN
          XLAB = '   S  (M)'
          axistype = 1
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
          axistype = 2
        ENDIF
c
        YLAB   = 'H Ion/Rec - PIN'
c
        NPLOTS = NPLOTS + 1
c        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)

        WRITE (IPLOT,9012) NPLOTS,REF

        CALL rzero (mvals,maxnks*maxngs*maxplts)
c
c
        write (6,*) 'Ionization/Recombination'
c
        do ip = 1, nplts
c
c         Access ring number and store number of knots info.
c
          ir = ringnos(ip)
c
c         Set number of points on plots
c
          pnks(ip,1) = nks(ir)
c
c         Load axis
c
          call loadm_axis(mouts,mwids,ir,ip,axistype,0)
c
c         Load data
c
          DO IK = 1, NKS(IR)
c
            MVALS(IK,ip,1) = pinion(ik,ir)
            MVALS(IK,ip,2) = pinrec(IK,IR)
c
c            write (6,'(2i4,3g13.6)') ir,ik,mouts(ik,ip),pinion(ik,ir),
c     >                       pinrec(ik,ir)
c
          enddo
C
        enddo
c
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
      endif
c
c
c     725 - E2D Pressures - cell centres -
c
c
      if (iref.eq.725) then
c
c        Plots Ne for E2D and DIVIMP for 8 rings on one
c        sheet of paper for easier comparison. It is set up
c        so that the plots can be modified for either vs. S or vs. P
c
        REF = 'COMPARATIVE PRESSURE PLOTS: FC vs. DIV'
c
        ELABS(1) = 'PrFCPrFC   '
        ELABS(2) = 'PACTPACT   '
        ELABS(3) = 'PEXPPEXP   '
c
c       Set up the ring numbers for the plots
c
c       Set up for 2 sets of data on each of 8 plots
c
        ngs = 2
        sctype = iopt
c
        IF (IREF.EQ.725) THEN
          XLAB = '   S  (M)'
          axistype = 1
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
          axistype = 2
        ENDIF
c
        YLAB   = 'Pressure'
c
        NPLOTS = NPLOTS + 1
c        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)

        WRITE (IPLOT,9012) NPLOTS,REF

        CALL rzero (mvals,maxnks*maxngs*maxplts)
c
c
        do ip = 1, nplts
c
c         Access ring number and store number of knots info.
c
          ir = ringnos(ip)
c
c         Modify the following to include values for S=0 ... even if this
c         is not part of the grid. Assume that if one end is not a grid
c         point then neither is the other.
c
          if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
             in = 0
             inc = 0
          else
             in = 1
             inc = 2
          endif
c
c         Set number of points on plots
c
          pnks(ip,1) = nks(ir)+inc
c
c         Load axis
c
          call loadm_axis(mouts,mwids,ir,ip,axistype,in)
c
c         Load data
c
c         Load target data for rings that need it.
c
          if (in.ne.0) then
c
c           Load actual target pressures.
c
            mVALS(1,ip,1) = knds(idds(ir,2))* ((kteds(idds(ir,2))
     >                             +ktids(idds(ir,2))) * ech +
     >             crmb * amu * kvds(idds(ir,2))**2)
            mVALS(1,ip,2) = mvals(1,ip,1)
c
c            mVALS(1,ip,3) = mvals(1,ip,1)
c
            mVALS(nks(ir)+2,ip,1) = knds(idds(ir,1))
     >                          *(( kteds(idds(ir,1))
     >                             +ktids(idds(ir,1))) * ech +
     >             crmb * amu * kvds(idds(ir,1))**2 )

            mVALS(nks(ir)+2,ip,2) =  mVALS(nks(ir)+2,ip,1)
c
c            mVALS(nks(ir)+2,ip,3) =  mVALS(nks(ir)+2,ip,1)
c
          endif
c
          DO IK = 1, NKS(IR)
c
c           E2D
c
            MVALS(ik+in,ip,1) = e2dnbs(ik,ir)*( ((e2dtebs(ik,ir)
     >                             +e2dtibs(ik,ir))) * ech +
     >             crmb * amu * e2dvhs(ik,ir)**2)
c
c           Actually used
c
            MVALS(IK+in,ip,2) = knbs(ik,ir)*( ((ktebs(ik,ir)
     >                             +ktibs(ik,ir))) * ech +
     >             crmb * amu * (kvhs(ik,ir)/qtim)**2)
c
c           Expected
c
c            MVALS(IK+in,ip,1) = kpress(ik,ir,1)
c

c
c            write(6,'(a,2i4,8(1x,g12.5))')  'COMP:',ik,ir,
c     >                  knbs(ik,ir),e2dnbs(ik,ir),
c     >                  kvhs(ik,ir),e2dvhs(ik,ir),
c     >                  ktibs(ik,ir),e2dtibs(ik,ir),
c     >                  ktebs(ik,ir),e2dtebs(ik,ir)
c
c
          enddo
C
        enddo
c
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
      endif
c
c
c
c
c********************************************************************
c
c
c     PLOT 751 - Plot of Cross-field transport coefficients
c                vs. R at the Outer Midplane
c
c     ONE PLOT - WHOLE RING
c
c
      IF (IREF.EQ.751) THEN
        ELABS(1) =  'DP  DPERP      '
        ELABS(2) =  'XPE XPERP E    '
        ELABS(3) =  'XPI XPERP I    '
c
        IF (IREF.EQ.751) THEN
          XLAB = 'R (M) OUTER MIDPLANE'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
c
        YLAB   = 'Transport Coeff.'
        NPLOTS = NPLOTS + 1
        WRITE (REF,'(''PLOT OF TRANSPORT COEFF.'')')
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (kvals, maxnks*maxngs)
c
        DO Ir = irsep, irwall-2
c
          in = ir - irsep + 1

          if (xpoint_up) then 

             KOUTS(in) = rcouter(ir)

          else

             KOUTS(in) = rcinner(ir)

          endif
c
c         Approximate width
c
          kwids(in) = (kinds(oumid,ir)+koutds(oumid,ir) )/2.0
c
          KVALS(In,1) = dperp(ir)
          KVALS(In,2) = chiperpe(ir)
          KVALS(In,3) = chiperpi(ir)
c
          write(6,'(a4,2i4,6(3x,f13.5))') 'TC:',in,ir,dperp(ir),
     >           chiperpe(ir),chiperpi(ir),rcouter(ir),rcinner(ir)

c
        enddo
C
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,(irwall-2)-irsep+1,ANLY,
     >    3,99,KOUTS(1),KOUTS((irwall-2)-irsep+1),-HI,HI,IGNORS,
     >    ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
c

      ENDIF
c
c     PLOT 753 - Plot of Cross-field transport coefficients
c                vs. R at the Outer Midplane (OUTER COEFFICIENTS)
c
c     ONE PLOT - OUTER
c
c
      IF (IREF.EQ.753) THEN
c
        ELABS(1) =  'DPO DPERP   OUT'
        ELABS(2) =  'XPEOXPERP E OUT'
        ELABS(3) =  'XPIOXPERP I OUT'
        ELABS(4) =  'XPTOXPERP T OUT'
c
        IF (IREF.EQ.753) THEN
          XLAB = 'R (M) OUTER MIDPLANE'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
c
        YLAB   = 'OUTER Transport'
        NPLOTS = NPLOTS + 1
        WRITE (REF,'(''PLOT OF OUTER TRANSPORT COEFF.'')')
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (kvals, maxnks*maxngs)
c
c        write (6,*) 'Plot 753:',iref,iopt,ir
c
c
        DO Ir = irsep, irwall-2
c
          in = ir - irsep + 1

          if (xpoint_up) then 

             KOUTS(in) = rcouter(ir)

          else

             KOUTS(in) = rcinner(ir)

          endif 

c
c         Approximate width
c
          kwids(in) = (kinds(oumid,ir)+koutds(oumid,ir) )/2.0
c
          if (xpoint_up) then 

             KVALS(In,1) = odperp(ir)
             KVALS(In,2) = ochiperpe(ir)
             KVALS(In,3) = ochiperpi(ir)
             KVALS(In,4) = oxperpt(ir)

          else

             KVALS(In,1) = idperp(ir)
             KVALS(In,2) = ichiperpe(ir)
             KVALS(In,3) = ichiperpi(ir)
             KVALS(In,4) = ixperpt(ir)

          endif
c
          if (xpoint_up) then 
c
             write(6,'(a4,2i4,5(3x,f13.5))') 'TC:',in,ir,odperp(ir),
     >          ochiperpe(ir),ochiperpi(ir),oxperpt(ir),rcouter(ir)

          else
             write(6,'(a4,2i4,5(3x,f13.5))') 'TC:',in,ir,idperp(ir),
     >          ichiperpe(ir),ichiperpi(ir),ixperpt(ir),rcinner(ir)
          endif  
c
        enddo
C
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,(irwall-2)-irsep+1,ANLY,
     >    4,99,KOUTS(1),KOUTS((irwall-2)-irsep+1),-HI,HI,IGNORS,
     >    ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
c
      ENDIF
c
c     PLOT 755 - Plot of Cross-field transport coefficients
c                vs. R at the Outer Midplane (INNER COEFFICIENTS)
c
c     ONE PLOT - INNER
c
      IF (IREF.EQ.755) THEN
        ELABS(1) =  'DPI DPERP   IN'
        ELABS(2) =  'XPEIXPERP E IN'
        ELABS(3) =  'XPIIXPERP I IN'
        ELABS(4) =  'XPTIXPERP T IN'
c
        IF (IREF.EQ.755) THEN
          XLAB = 'R (M) OUTER MIDPLANE'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
c
        YLAB   = 'INNER Transport'
        NPLOTS = NPLOTS + 1
        WRITE (REF,'(''PLOT OF INNER TRANSPORT COEFF.'')')
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (kvals, maxnks*maxngs)
c
c        write (6,*) 'Plot 755:',iref,iopt,ir
c
c
        DO Ir = irsep, irwall-2
c
          in = ir - irsep + 1

          if (xpoint_up) then 

             KOUTS(in) = rcouter(ir)
      
          else

             KOUTS(in) = rcinner(ir)

          endif 

c
c         Approximate width
c
          kwids(in) = (kinds(oumid,ir)+koutds(oumid,ir) )/2.0
c
          if (xpoint_up) then 

             KVALS(In,1) = idperp(ir)
             KVALS(In,2) = ichiperpe(ir)
             KVALS(In,3) = ichiperpi(ir)
             KVALS(In,4) = ixperpt(ir)

          else

             KVALS(In,1) = odperp(ir)
             KVALS(In,2) = ochiperpe(ir)
             KVALS(In,3) = ochiperpi(ir)
             KVALS(In,4) = oxperpt(ir)

          endif          
c
          if (xpoint_up) then 
c
             write(6,'(a4,2i4,5(3x,f13.5))') 'TC:',in,ir,idperp(ir),
     >          ichiperpe(ir),ichiperpi(ir),ixperpt(ir),rcouter(ir)

          else

             write(6,'(a4,2i4,5(3x,f13.5))') 'TC:',in,ir,odperp(ir),
     >          ochiperpe(ir),ochiperpi(ir),oxperpt(ir),rcinner(ir)

          endif  
c
        enddo
C
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,(irwall-2)-irsep+1,ANLY,
     >    4,99,KOUTS(1),KOUTS((irwall-2)-irsep+1),-HI,HI,IGNORS,
     >    ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
c
      ENDIF
c
c     PLOT 757 - Plot of Cross-field transport coefficients
c                vs. R at the Outer Midplane
c
c     4 PLOTS/PAGE - WHOLE RING
c
c
      IF (IREF.EQ.757) THEN
c
        ELABS(1) = '    WHOLE RING'
c
        PLTLABS(1) =  'DPERP      '
        PLTLABS(2) =  'XPERP E    '
        PLTLABS(3) =  'XPERP I    '
        PLTLABS(4) =  'XPERP TOTAL'
c
        IF (IREF.EQ.757) THEN
          XLAB = 'R (M) OUTER MIDPLANE'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
c
        YLAB   = 'Transport Coeff.'
        NPLOTS = NPLOTS + 1
        WRITE (REF,'(''PLOT OF TRANSPORT COEFF.'')')
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (mvals, maxnks*maxngs*maxplts)
c
        sctype = iopt
        ngrm = 4
        nplts = 4
        ngs = 1
c
        irlim = irwall-2
c
        DO Ir = irsep, irlim
c
          in = ir - irsep + 1
c
          if (xpoint_up) then 
             
             mOUTS(in,1,1) = rcouter(ir)
             mOUTS(in,2,1) = rcouter(ir)
             mOUTS(in,3,1) = rcouter(ir)
             mOUTS(in,4,1) = rcouter(ir)

          else

             mOUTS(in,1,1) = rcinner(ir)
             mOUTS(in,2,1) = rcinner(ir)
             mOUTS(in,3,1) = rcinner(ir)
             mOUTS(in,4,1) = rcinner(ir)

          endif 

c
c         Approximate width
c
          mwids(in,1,1) = (kinds(oumid,ir)+koutds(oumid,ir) )/2.0
          mwids(in,2,1) = (kinds(oumid,ir)+koutds(oumid,ir) )/2.0
          mwids(in,3,1) = (kinds(oumid,ir)+koutds(oumid,ir) )/2.0
          mwids(in,4,1) = (kinds(oumid,ir)+koutds(oumid,ir) )/2.0
c
          mVALS(In,1,1) = dperp(ir)
          mVALS(In,2,1) = chiperpe(ir)
          mVALS(In,3,1) = chiperpi(ir)
          mVALS(In,4,1) = xperpt(ir)
c
          write(6,'(a4,2i4,4(3x,f13.5))') 'TC757:',in,ir,dperp(ir),
     >               chiperpe(ir),chiperpi(ir),rcouter(ir)
c
        enddo
c
        pnks(1,1) = irlim -irsep + 1
        pnks(2,1) = irlim -irsep + 1
        pnks(3,1) = irlim -irsep + 1
        pnks(4,1) = irlim -irsep + 1
c
        pltmins(1) = 0.0
        pltmaxs(1) = 0.3
c
        pltmins(2) = 0.0
        pltmaxs(2) = 3.0
c
        pltmins(3) = 0.0
        pltmaxs(3) = 3.0
c
        pltmins(4) = 0.0
        pltmaxs(4) = 3.0
C
C
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
      ENDIF
c
c     PLOT 759 - Plot of Cross-field transport coefficients
c                vs. R at the Outer Midplane
c
c     4 PLOTS - OUTER
c
c
      IF (IREF.EQ.759) THEN
c
        ELABS(1) = '    OUTER'
c
        PLTLABS(1) =  'DPERP   OUT'
        PLTLABS(2) =  'XPERP E OUT'
        PLTLABS(3) =  'XPERP I OUT'
        PLTLABS(4) =  'XPERP T OUT'
c
        IF (IREF.EQ.759) THEN
          XLAB = 'R (M) OUTER MIDPLANE'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
c
        YLAB   = 'Transport Coeff.'
        NPLOTS = NPLOTS + 1
        WRITE (REF,'(''PLOT OF TRANSPORT COEFF.'')')
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (mvals, maxnks*maxngs*maxplts)
c
        sctype = iopt
        ngrm = 4
        nplts = 4
        ngs = 1
c
        DO Ir = irsep, irwall-2
c
          in = ir - irsep + 1
c
          if (xpoint_up) then 
             
             mOUTS(in,1,1) = rcouter(ir)
             mOUTS(in,2,1) = rcouter(ir)
             mOUTS(in,3,1) = rcouter(ir)
             mOUTS(in,4,1) = rcouter(ir)

             mVALS(In,1,1) = odperp(ir)
             mVALS(In,2,1) = ochiperpe(ir)
             mVALS(In,3,1) = ochiperpi(ir)
             mVALS(In,4,1) = oxperpt(ir)

             write(6,'(a4,2i4,5(3x,f13.5))') 'TC:',in,ir,odperp(ir),
     >          ochiperpe(ir),ochiperpi(ir),oxperpt(ir),rcouter(ir)
          else

             mOUTS(in,1,1) = rcinner(ir)
             mOUTS(in,2,1) = rcinner(ir)
             mOUTS(in,3,1) = rcinner(ir)
             mOUTS(in,4,1) = rcinner(ir)

             mVALS(In,1,1) = idperp(ir)
             mVALS(In,2,1) = ichiperpe(ir)
             mVALS(In,3,1) = ichiperpi(ir)
             mVALS(In,4,1) = ixperpt(ir)

             write(6,'(a4,2i4,5(3x,f13.5))') 'TC:',in,ir,idperp(ir),
     >          ichiperpe(ir),ichiperpi(ir),ixperpt(ir),rcinner(ir)

          endif 
c
c         Approximate width
c
          mwids(in,1,1) = (kinds(oumid,ir)+koutds(oumid,ir) )/2.0
          mwids(in,2,1) = (kinds(oumid,ir)+koutds(oumid,ir) )/2.0
          mwids(in,3,1) = (kinds(oumid,ir)+koutds(oumid,ir) )/2.0
          mwids(in,4,1) = (kinds(oumid,ir)+koutds(oumid,ir) )/2.0
c
        enddo
c
        pnks(1,1) = (irwall-2) -irsep + 1
        pnks(2,1) = (irwall-2) -irsep + 1
        pnks(3,1) = (irwall-2) -irsep + 1
        pnks(4,1) = (irwall-2) -irsep + 1
c
        pltmins(1) = 0.0
        pltmaxs(1) = 0.2
c
        pltmins(2) = 0.0
        pltmaxs(2) = 2.0
c
        pltmins(3) = 0.0
        pltmaxs(3) = 2.0
C
        pltmins(4) = 0.0
        pltmaxs(4) = 2.0
C
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
      ENDIF
c
c
c     PLOT 761 - Plot of Cross-field transport coefficients
c                vs. R at the Outer Midplane
c
c     4 PLOTS - INNER
c
      IF (IREF.EQ.761) THEN
c
        ELABS(1) = '    INNER'
c
        PLTLABS(1) =  'DPERP   INN'
        PLTLABS(2) =  'XPERP E INN'
        PLTLABS(3) =  'XPERP I INN'
        PLTLABS(4) =  'XPERP T INN'
c
        IF (IREF.EQ.761) THEN
          XLAB = 'R (M) OUTER MIDPLANE'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
c
        YLAB   = 'Transport Coeff.'
        NPLOTS = NPLOTS + 1
        WRITE (REF,'(''PLOT OF TRANSPORT COEFF.'')')
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (mvals, maxnks*maxngs*maxplts)
c
        sctype = iopt
        ngrm = 4
        nplts = 4
        ngs = 1
c
        DO Ir = irsep, irwall-2
c
          in = ir - irsep + 1
c
          if (xpoint_up) then 
             
             mOUTS(in,1,1) = rcouter(ir)
             mOUTS(in,2,1) = rcouter(ir)
             mOUTS(in,3,1) = rcouter(ir)
             mOUTS(in,4,1) = rcouter(ir)

             mVALS(In,1,1) = idperp(ir)
             mVALS(In,2,1) = ichiperpe(ir)
             mVALS(In,3,1) = ichiperpi(ir)
             mVALS(In,4,1) = ixperpt(ir)
c
             write(6,'(a4,2i4,4(3x,f13.5))') 'TC:',in,ir,idperp(ir),
     >           ichiperpe(ir),ichiperpi(ir),ixperpt(ir),rcouter(ir)
c
          else

             mOUTS(in,1,1) = rcinner(ir)
             mOUTS(in,2,1) = rcinner(ir)
             mOUTS(in,3,1) = rcinner(ir)
             mOUTS(in,4,1) = rcinner(ir)

             mVALS(In,1,1) = odperp(ir)
             mVALS(In,2,1) = ochiperpe(ir)
             mVALS(In,3,1) = ochiperpi(ir)
             mVALS(In,4,1) = oxperpt(ir)
c
             write(6,'(a4,2i4,4(3x,f13.5))') 'TC:',in,ir,odperp(ir),
     >           ochiperpe(ir),ochiperpi(ir),oxperpt(ir),rcinner(ir)
c

          endif 

c
c         Approximate width
c
          mwids(in,1,1) = (kinds(oumid,ir)+koutds(oumid,ir) )/2.0
          mwids(in,2,1) = (kinds(oumid,ir)+koutds(oumid,ir) )/2.0
          mwids(in,3,1) = (kinds(oumid,ir)+koutds(oumid,ir) )/2.0
          mwids(in,4,1) = (kinds(oumid,ir)+koutds(oumid,ir) )/2.0
c
        enddo
c
        pnks(1,1) = (irwall-2) -irsep + 1
        pnks(2,1) = (irwall-2) -irsep + 1
        pnks(3,1) = (irwall-2) -irsep + 1
        pnks(4,1) = (irwall-2) -irsep + 1
c
        pltmins(1) = 0.0
        pltmaxs(1) = 0.2
c
        pltmins(2) = 0.0
        pltmaxs(2) = 2.0
c
        pltmins(3) = 0.0
        pltmaxs(3) = 2.0
c
        pltmins(4) = 0.0
        pltmaxs(4) = 2.0
C
C
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
      ENDIF
c
c
c     PLOT 763 - Plot of Cross-field transport coefficients
c                vs. R at the Outer Midplane
c
c     2 PLOTS - JUST DPERP and XPERP TOTALS - WHOLE RING
c
      IF (IREF.EQ.763) THEN
c
        ELABS(1) = '    WHOLE RING'
c
        PLTLABS(1) =  'DPERP   '
        PLTLABS(2) =  'XPERP   '
c
        IF (IREF.EQ.763) THEN
          XLAB = 'R (M) OUTER MIDPLANE'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
c
        YLAB   = 'Transport Coeff.'
        NPLOTS = NPLOTS + 1
        WRITE (REF,'(''PLOT OF TRANSPORT COEFF.'')')
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (mvals, maxnks*maxngs*maxplts)
c
        sctype = iopt
        ngrm = 2
        nplts = 2
        ngs = 1
c
        DO Ir = irsep, irwall-2
c
          in = ir - irsep + 1
c
          if (xpoint_up) then 

             mOUTS(in,1,1) = rcouter(ir)
             mOUTS(in,2,1) = rcouter(ir)

          else

             mOUTS(in,1,1) = rcouter(ir)
             mOUTS(in,2,1) = rcouter(ir)

          endif
c
c         Approximate width
c
          mwids(in,1,1) = (kinds(oumid,ir)+koutds(oumid,ir) )/2.0
          mwids(in,2,1) = (kinds(oumid,ir)+koutds(oumid,ir) )/2.0
c
          mVALS(In,1,1) = dperp(ir)
          mVALS(In,2,1) = xperpt(ir)
c
          write(6,'(a4,2i4,4(3x,f13.5))') 'TC763:',in,ir,dperp(ir),
     >           chiperpe(ir),chiperpi(ir),rcouter(ir),rcinner(ir)
c
        enddo
c
        pnks(1,1) = (irwall-2) -irsep + 1
        pnks(2,1) = (irwall-2) -irsep + 1
c
        pltmins(1) = 0.0
        pltmaxs(1) = 0.2
c
        pltmins(2) = 0.0
        pltmaxs(2) = 2.0
c
        pltmins(3) = 0.0
        pltmaxs(3) = 2.0
C
C
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
      ENDIF
c
c
c
c     PLOT 781 - Cross-field plots of density and interpolation
c                vs. D.
c
c
c
      IF (IREF.EQ.781) THEN
c
        ELABS(1) = 'NE2DNE2D'
        ELABS(2) = 'NDIVNDIV'
c
        CALL rzero (mvals, maxnks*maxngs*maxplts)
c
        sctype = iopt
        ngrm = 8
        nplts = 8
        ngs = 2
        nstep = nplts
c
        do ikstart = 1,nks(irsep),nstep
c
          ikend = min(ikstart + nplts -1,nks(irsep))
c
          nplts = min(nplts,ikend-ikstart+1)
c
          do ik = ikstart,ikend
c
          ip = ik - ikstart + 1
c
          write(PLTLABS(ip),'(a,i4)') 'KNOT=',ik
c
          IF (IREF.EQ.781) THEN
            XLAB = 'DIST (M)'
          ELSE
            XLAB = '   POLOIDAL DIST (M)'
          ENDIF
c
          YLAB   = 'Denisty'
          NPLOTS = NPLOTS + 1
          WRITE (REF,'(''PLOT OF DENSITY Cross-field'')')
          WRITE (IPLOT,9012) NPLOTS,REF
c
c         Set up
c
          dist = 0.0
c
          do ir = irsep,irwall-1
c
             in = ir - irsep +1
c
             mOUTS(in,ip,1) = dist
c
             mwids(in,ip,1) = (kinds(ik,ir)+koutds(ik,ir) )/2.0
c
             mVALS(In,ip,1) = e2dnbs(ik,ir)
             mVALS(In,ip,2) = knbs(ik,ir)
c
             dist = dist + koutds(ik,ir)
c
          end do
c
          pnks(ip,1) = (irwall-1) -irsep + 1
c
c
        enddo
c
C
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
      end do

c
      ENDIF
c
c
c
c     PLOT 783 - Cross-field plots of Te
c                vs. D.
c
c
c
      IF (IREF.EQ.783) THEN
c
        ELABS(1) = 'TEE TE-E2D'
        ELABS(2) = 'TED TE-DIV'
c
        CALL rzero (mvals, maxnks*maxngs*maxplts)
c
        sctype = iopt
        ngrm = 8
        nplts = 8
        ngs = 2
        nstep = nplts
c
        do ikstart = 1,nks(irsep),nstep
c
          ikend = min(ikstart + nplts -1,nks(irsep))
c
          nplts = min(nplts,ikend-ikstart+1)
c
          do ik = ikstart,ikend
c
          ip = ik - ikstart + 1
c
          write(PLTLABS(ip),'(a,i4)') 'KNOT=',ik
c
          IF (IREF.EQ.783) THEN
            XLAB = 'DIST (M)'
          ELSE
            XLAB = '   POLOIDAL DIST (M)'
          ENDIF
c
          YLAB   = 'Te'
          NPLOTS = NPLOTS + 1
          WRITE (REF,'(''PLOT OF Te Cross-field'')')
          WRITE (IPLOT,9012) NPLOTS,REF
c
c         Set up
c
          dist = 0.0
c
          do ir = irsep,irwall-1
c
             in = ir - irsep +1
c
             mOUTS(in,ip,1) = dist
c
             mwids(in,ip,1) = (kinds(ik,ir)+koutds(ik,ir) )/2.0
c
             mVALS(In,ip,1) = e2dtebs(ik,ir)
             mVALS(In,ip,2) = ktebs(ik,ir)
c
             dist = dist + koutds(ik,ir)
c
          end do
c
          pnks(ip,1) = (irwall-1) -irsep + 1
c
c
        enddo
c
C
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
      end do

c
      ENDIF
c
c
c
c     PLOT 785 - Cross-field plots Ti
c                vs. D.
c
c
c
      IF (IREF.EQ.785) THEN
c
        ELABS(1) = 'TIE TI-E2D'
        ELABS(2) = 'TID TI-DIV'
c
        CALL rzero (mvals, maxnks*maxngs*maxplts)
c
        sctype = iopt
        ngrm = 8
        nplts = 8
        ngs = 2
        nstep = nplts
c
        do ikstart = 1,nks(irsep),nstep
c
          ikend = min(ikstart + nplts -1,nks(irsep))
c
          nplts = min(nplts,ikend-ikstart+1)
c
          do ik = ikstart,ikend
c
          ip = ik - ikstart + 1
c
          write(PLTLABS(ip),'(a,i4)') 'KNOT=',ik
c
          IF (IREF.EQ.785) THEN
            XLAB = 'DIST (M)'
          ELSE
            XLAB = '   POLOIDAL DIST (M)'
          ENDIF
c
          YLAB   = 'Ti'
          NPLOTS = NPLOTS + 1
          WRITE (REF,'(''PLOT OF Ti Cross-field'')')
          WRITE (IPLOT,9012) NPLOTS,REF
c
c         Set up
c
          dist = 0.0
c
          do ir = irsep,irwall-1
c
             in = ir - irsep +1
c
             mOUTS(in,ip,1) = dist
c
             mwids(in,ip,1) = (kinds(ik,ir)+koutds(ik,ir) )/2.0
c
             mVALS(In,ip,1) = e2dtibs(ik,ir)
             mVALS(In,ip,2) = ktibs(ik,ir)
c
             dist = dist + koutds(ik,ir)
c
          end do
c
          pnks(ip,1) = (irwall-1) -irsep + 1
c
c
        enddo
c
C
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
      end do

c
      ENDIF
c
c
c
c     PLOT 787 - Cross-field plots of density gradient
c                vs. D.
c
c
c
      IF (IREF.EQ.787) THEN
c
        ELABS(1) = 'GNE GN-E2D'
        ELABS(2) = 'GND GN-DIV'
c
        CALL rzero (mvals, maxnks*maxngs*maxplts)
c
        sctype = iopt
        ngrm = 8
        nplts = 8
        ngs = 2
        nstep = nplts
c
        do ikstart = 1,nks(irsep),nstep
c
          ikend = min(ikstart + nplts -1,nks(irsep))
c
          nplts = min(nplts,ikend-ikstart+1)
c
          do ik = ikstart,ikend
c
          ip = ik - ikstart + 1
c
          write(PLTLABS(ip),'(a,i4)') 'KNOT=',ik
c
          IF (IREF.EQ.787) THEN
            XLAB = 'DIST (M)'
          ELSE
            XLAB = '   POLOIDAL DIST (M)'
          ENDIF
c
          YLAB   = 'Denisty Gradient'
          NPLOTS = NPLOTS + 1
          WRITE (REF,'(''PLOT OF DENSITY GRADIENT'')')
          WRITE (IPLOT,9012) NPLOTS,REF
c
c         Set up
c
          dist = 0.0
c
          do ir = irsep,irwall-1
c
             in = ir - irsep +1
c
             mOUTS(in,ip,1) = dist
c
             mwids(in,ip,1) = (kinds(ik,ir)+koutds(ik,ir) )/2.0
c
             mVALS(In,ip,1) = e2dgradn(ik,ir)
             mVALS(In,ip,2) = gradn(ik,ir)
c
             dist = dist + koutds(ik,ir)
c
         end do
c
         pnks(ip,1) = (irwall-1) -irsep + 1
c
c
       enddo
c
C
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
      end do

c
      ENDIF
c
c
c     PLOT 789 - Cross-field plots of Te Gradient
c                vs. D.
c
c
c
      IF (IREF.EQ.789) THEN
c
        ELABS(1) = 'GTeEGTeE2D'
        ELABS(2) = 'GTeDGTeDIV'
c
        CALL rzero (mvals, maxnks*maxngs*maxplts)
c
        sctype = iopt
        ngrm = 8
        nplts = 8
        ngs = 2
        nstep = nplts
c
        do ikstart = 1,nks(irsep),nstep
c
          ikend = min(ikstart + nplts -1,nks(irsep))
c
          nplts = min(nplts,ikend-ikstart+1)
c
          do ik = ikstart,ikend
c
          ip = ik - ikstart + 1
c
          write(PLTLABS(ip),'(a,i4)') 'KNOT=',ik
c
          IF (IREF.EQ.789) THEN
            XLAB = 'DIST (M)'
          ELSE
            XLAB = '   POLOIDAL DIST (M)'
          ENDIF
c
          YLAB   = 'Te Gradient'
          NPLOTS = NPLOTS + 1
          WRITE (REF,'(''PLOT OF Te Gradient CF'')')
          WRITE (IPLOT,9012) NPLOTS,REF
c
c         Set up
c
          dist = 0.0
c
          do ir = irsep,irwall-1
c
             in = ir - irsep +1
c
             mOUTS(in,ip,1) = dist
c
             mwids(in,ip,1) = (kinds(ik,ir)+koutds(ik,ir) )/2.0
c
             mVALS(In,ip,1) = e2dgradte(ik,ir)
             mVALS(In,ip,2) = gradte(ik,ir)
c
             dist = dist + koutds(ik,ir)
c
         end do
c
         pnks(ip,1) = (irwall-1) -irsep + 1
c
c
       enddo
c
C
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
      end do

c
      ENDIF
c
c
c     PLOT 791 - Cross-field plots of Ti Gradient
c                vs. D.
c
c
c
      IF (IREF.EQ.791) THEN
c
        ELABS(1) = 'GTiEGTiE2D'
        ELABS(2) = 'GTiDGTiDIV'
c
        CALL rzero (mvals, maxnks*maxngs*maxplts)
c
        sctype = iopt
        ngrm = 8
        nplts = 8
        ngs = 2
        nstep = nplts
c
        do ikstart = 1,nks(irsep),nstep
c
          ikend = min(ikstart + nplts -1,nks(irsep))
c
          nplts = min(nplts,ikend-ikstart+1)
c
          do ik = ikstart,ikend
c
          ip = ik - ikstart + 1
c
          write(PLTLABS(ip),'(a,i4)') 'KNOT=',ik
c
          IF (IREF.EQ.791) THEN
            XLAB = 'DIST (M)'
          ELSE
            XLAB = '   POLOIDAL DIST (M)'
          ENDIF
c
          YLAB   = 'Ti Gradient'
          NPLOTS = NPLOTS + 1
          WRITE (REF,'(''PLOT OF Ti Gradient CF'')')
          WRITE (IPLOT,9012) NPLOTS,REF
c
c         Set up
c
          dist = 0.0
c
          do ir = irsep,irwall-1
c
             in = ir - irsep +1
c
             mOUTS(in,ip,1) = dist
c
             mwids(in,ip,1) = (kinds(ik,ir)+koutds(ik,ir) )/2.0
c
             mVALS(In,ip,1) = e2dgradti(ik,ir)
             mVALS(In,ip,2) = gradti(ik,ir)
c
             dist = dist + koutds(ik,ir)
c
         end do
c
         pnks(ip,1) = (irwall-1) -irsep + 1
c
c
       enddo
c
C
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
      end do

c
      ENDIF
c
c
c
c     PLOT 793- Plot of Cross-field transport coefficients
c                vs. R at the Outer Midplane
c
c     4 PLOTS/PAGE - WHOLE RING
c
c
      IF (IREF.EQ.793) THEN
c
        ELABS(1) = '    Transport'
c
        PLTLABS(1) =  'Dperp'
        PLTLABS(2) =  'Xperpe'
        PLTLABS(3) =  'Xperpi'
        PLTLABS(4) =  'Xperpt'
c
        IF (IREF.EQ.793) THEN
          XLAB = 'Mid-plane distance'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
c
        YLAB   = ' '
        NPLOTS = NPLOTS + 1
        WRITE (REF,'('' '')')
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (mvals, maxnks*maxngs*maxplts)
c
        sctype = iopt
        ngrm = 6
        nplts = 4
        ngs = 1
c
        irlim = irwall-1
c
        DO Ir = irsep, irlim
c
          in = ir - irsep + 1
c
          mOUTS(in,1,1) = rcouter(ir)
          mOUTS(in,2,1) = rcouter(ir)
          mOUTS(in,3,1) = rcouter(ir)
          mOUTS(in,4,1) = rcouter(ir)
c
c         Approximate width
c
          mwids(in,1,1) = (kinds(oumid,ir)+koutds(oumid,ir) )/2.0
          mwids(in,2,1) = (kinds(oumid,ir)+koutds(oumid,ir) )/2.0
          mwids(in,3,1) = (kinds(oumid,ir)+koutds(oumid,ir) )/2.0
          mwids(in,4,1) = (kinds(oumid,ir)+koutds(oumid,ir) )/2.0
c
          mVALS(In,1,1) = dperp(ir)
          mVALS(In,2,1) = chiperpe(ir)
          mVALS(In,3,1) = chiperpi(ir)
          mVALS(In,4,1) = xperpt(ir)
c
c          write(6,'(a4,2i4,4(3x,f13.5))') 'TC:',in,ir,dperp(ir),
c     >               chiperpe(ir),chiperpi(ir),rcouter(ir)
c
        enddo
c
        pnks(1,1) = irlim -irsep + 1
        pnks(2,1) = irlim -irsep + 1
        pnks(3,1) = irlim -irsep + 1
        pnks(4,1) = irlim -irsep + 1
c
        pltmins(1) = 0.0
        pltmaxs(1) = 0.3
c
        pltmins(2) = 0.0
        pltmaxs(2) = 3.0
c
        pltmins(3) = 0.0
        pltmaxs(3) = 3.0
c
        pltmins(4) = 0.0
        pltmaxs(4) = 3.0
C
C
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
      ENDIF
c
c
c     PLOT 797 - Plot of Background Quantities on a specified RING
c                vs. S - Four Plots - Te Ti Ne Ga
c
c     4 PLOTS/PAGE
c
c
      IF (IREF.EQ.797) THEN
c
        ELABS(1) = 'FC  FC'
        ELABS(2) = 'OSM OSM'

c
        PLTLABS(1) =  'TE         '
        PLTLABS(2) =  'TI         '
        PLTLABS(3) =  'NE         '
        PLTLABS(4) =  'GAMMA      '
c
        YLAB   = 'Plasma Quantity'
        NPLOTS = NPLOTS + 1
        WRITE (REF,'(''PLOT OF BACKGROUND'')')
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (mvals, maxnks*maxngs*maxplts)
c
        sctype = 4
        ngrm = 6
        nplts = 4
        ngs = 2
c
        ir = iopt
c
        IF (IREF.EQ.797) THEN
          XLAB = '   S  (M)'
          axistype = 1
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
          axistype = 2
        ENDIF
c
        CALL rzero (mvals,maxnks*maxngs*maxplts)
c
        do ip = 1, nplts
c
c         Access ring number and store number of knots info.
c
c
c         Modify the following to include values for S=0 ... even if this
c         is not part of the grid. Assume that if one end is not a grid
c         point then neither is the other.
c
          if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
             in = 0
             inc = 0
          else
             in = 1
             inc = 2
          endif
c
c         Set number of points on plots
c
          pnks(ip,1) = nks(ir)+inc
c
c         Load axis
c
          call loadm_axis(mouts,mwids,ir,ip,axistype,in)
c
c         Load data
c
c         Load target data for rings that need it.
c
          if (in.ne.0) then
c
            if (ip.eq.1) then
c
               mVALS(1,ip,1) = e2dtebs(1,ir)
               mVALS(1,ip,2) = KtedS (idds(ir,2))

               mVALS(nks(ir)+2,ip,1) = e2dtebs(nks(ir),ir)
               mVALS(nks(ir)+2,ip,2) = KtedS (idds(ir,1))

            elseif (ip.eq.2) then
c
               mVALS(1,ip,1) = e2dtibs(1,ir)
               mVALS(1,ip,2) = KtidS (idds(ir,2))

               mVALS(nks(ir)+2,ip,1) = e2dtibs(nks(ir),ir)
               mVALS(nks(ir)+2,ip,2) = KtidS (idds(ir,1))

            elseif (ip.eq.3) then
c
               mVALS(1,ip,1) = e2dnbs(1,ir)
               mVALS(1,ip,2) = KNdS (idds(ir,2))

               mVALS(nks(ir)+2,ip,1) = e2dnbs(nks(ir),ir)
               mVALS(nks(ir)+2,ip,2) = KNdS (idds(ir,1))

            elseif (ip.eq.4) then
c
               mVALS(1,ip,1) = e2dtarg(ir,5,2)
               mVALS(1,ip,2) = KNdS(idds(ir,2))*KvdS(idds(ir,2))

               mVALS(nks(ir)+2,ip,1) = e2dtarg(ir,5,1)
               mVALS(nks(ir)+2,ip,2) = KNdS(idds(ir,1))*KvdS(idds(ir,1))
c
            endif
c
          endif
c
          pnks(ip,1) = nks(ir)+inc
c
          DO IK = 1, NKS(IR)
c
            if (ip.eq.1) then

               MVALS(IK+in,ip,1) = e2dtebs(ik,ir)
               MVALS(IK+in,ip,2) = KteBS (IK,IR)

            elseif (ip.eq.2) then

               MVALS(IK+in,ip,1) = e2dtibs(ik,ir)
               MVALS(IK+in,ip,2) = KtiBS (IK,IR)

            elseif (ip.eq.3) then

               MVALS(IK+in,ip,1) = e2dnbs(ik,ir)
               MVALS(IK+in,ip,2) = KNBS (IK,IR)

            elseif (ip.eq.4) then

               MVALS(IK+in,ip,1) = e2dnbs(ik,ir) * e2dvhs(ik,ir)
               MVALS(IK+in,ip,2) = KNBS (IK,IR)* kvhs(ik,ir)/qtim

            endif
c
          enddo
C
        enddo
c
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
      endif
c
c     PLOT 798 - Plot of Background Quantities on a specified RING
c                vs. S - Four Plots - Te Ti Ne Ga
c
c     4 PLOTS/PAGE
c
c
      IF (IREF.EQ.798) THEN
c
c        ELABS(1) = 'FC  FC'
        ELABS(1) = 'OSM OSM'

c
        PLTLABS(1) =  'T_e        '
        PLTLABS(2) =  'T_i        '
        PLTLABS(3) =  'n_e        '
        PLTLABS(4) =  'Parallel Flux'
c
        YLAB   = 'Plasma Quantity'
        NPLOTS = NPLOTS + 1
        WRITE (REF,'(''PLOT OF BACKGROUND ON RING: '',i4)') iopt
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (mvals, maxnks*maxngs*maxplts)
c
        sctype = 4
        ngrm = 6
        nplts = 4
        ngs = 1
c
        ir = iopt
c
        IF (IREF.EQ.798) THEN
          XLAB = '   S  (M)'
          axistype = 1
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
          axistype = 2
        ENDIF
c
        CALL rzero (mvals,maxnks*maxngs*maxplts)
c
        do ip = 1, nplts
c
c         Access ring number and store number of knots info.
c
c
c         Modify the following to include values for S=0 ... even if this
c         is not part of the grid. Assume that if one end is not a grid
c         point then neither is the other.
c
          if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
             in = 0
             inc = 0
          else
             in = 1
             inc = 2
          endif
c
c         Set number of points on plots
c
          pnks(ip,1) = nks(ir)+inc
c
c         Load axis
c
          call loadm_axis(mouts,mwids,ir,ip,axistype,in)
c
c         Load data
c
c         Load target data for rings that need it.
c
          if (in.ne.0) then
c
            if (ip.eq.1) then
c
c               mVALS(1,ip,1) = e2dtebs(1,ir)
               mVALS(1,ip,1) = KtedS (idds(ir,2))

c               mVALS(nks(ir)+2,ip,1) = e2dtebs(nks(ir),ir)
               mVALS(nks(ir)+2,ip,1) = KtedS (idds(ir,1))

            elseif (ip.eq.2) then
c
c               mVALS(1,ip,1) = e2dtibs(1,ir)
               mVALS(1,ip,1) = KtidS (idds(ir,2))

c               mVALS(nks(ir)+2,ip,1) = e2dtibs(nks(ir),ir)
               mVALS(nks(ir)+2,ip,1) = KtidS (idds(ir,1))

            elseif (ip.eq.3) then
c
c               mVALS(1,ip,1) = e2dnbs(1,ir)
               mVALS(1,ip,1) = KNdS (idds(ir,2))

c               mVALS(nks(ir)+2,ip,1) = e2dnbs(nks(ir),ir)
               mVALS(nks(ir)+2,ip,1) = KNdS (idds(ir,1))

            elseif (ip.eq.4) then
c
c               mVALS(1,ip,1) = e2dtarg(ir,5,2)
               mVALS(1,ip,1) = KNdS(idds(ir,2))*KvdS(idds(ir,2))

c               mVALS(nks(ir)+2,ip,1) = e2dtarg(ir,5,1)
               mVALS(nks(ir)+2,ip,1) = KNdS(idds(ir,1))*KvdS(idds(ir,1))
c
            endif
c
          endif
c
          pnks(ip,1) = nks(ir)+inc
c
          DO IK = 1, NKS(IR)
c
            if (ip.eq.1) then

c               MVALS(IK+in,ip,1) = e2dtebs(ik,ir)
               MVALS(IK+in,ip,1) = KteBS (IK,IR)

            elseif (ip.eq.2) then

c               MVALS(IK+in,ip,1) = e2dtibs(ik,ir)
               MVALS(IK+in,ip,1) = KtiBS (IK,IR)

            elseif (ip.eq.3) then

c               MVALS(IK+in,ip,1) = e2dnbs(ik,ir)
               MVALS(IK+in,ip,1) = KNBS (IK,IR)

            elseif (ip.eq.4) then

c               MVALS(IK+in,ip,1) = e2dnbs(ik,ir) * e2dvhs(ik,ir)
               MVALS(IK+in,ip,1) = KNBS (IK,IR)* kvhs(ik,ir)/qtim

            endif
c
          enddo
C
        enddo
c
c
c        Set up data for modified call to DRAWM
c
         do ip = 1,nplts
            do ig = 1, ngs            
               mlabs(ip,ig) = elabs(ig)
            end do  
            pngs(ip) = ngs
         end do
c
        CALL DRAWM (MOUTS,MWIDS,MVALS,MAXDATX,maxplts,maxngs,pnks,
     >              nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >              sctype,ngrm,pltmins,pltmaxs,pltfact,1,
     >              mdrawtype,0)
c
      endif

      return

c
c     Format statements    
c
 9012 FORMAT(1X,'PLOT',I3,4X,A)

      end


