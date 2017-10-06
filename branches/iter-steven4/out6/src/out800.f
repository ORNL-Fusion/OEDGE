c     -*-Fortran-*-
c
      subroutine out800(iref,graph,iopt,ierr)
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
      include 'dynam3'
c      include 'dynam4'
c      include 'pindata'
c      include 'cadas'
c      include 'grbound'
c      include 'outxy'
c      include 'cedge2d'
c      include 'transcoef'
c      include 'cioniz'
c      include 'reiser' 
      include 'printopt'  
c
c     Local Variables
c
      integer ik,ir,iz
      integer in,id
      integer ind
c
c     Local Variables
c

      real tmpsum


c
c     Variables for the 800 series bar charts
c
      integer maxsets
      parameter (maxsets=3)
c
      real valsts(maxpts,maxsets),ymin,ymax,totsrc,totleak
      real totdep
      real wlzmax
      integer iw,iwstart,iwend
      integer irstart,irend
      integer istart,istop, novals,nosets,is,wlmax
      character*15 pnames1(maxpts),cnames(maxsets)
      character*15 pnames2(maxpts)
      character*80 grtitle
c
c     Variables used in plot 821
c
      integer targid, wallind
      real geofact 
c
c     Variables for plot 831
c
      integer iselect, iscale, ival, iflag
c
      character*100 cmd,fn
c
      real tot1,tot2
      real,allocatable :: depdata(:,:,:),absfacs(:,:)
c
      real scalef_821

c
      IF (IOPT.EQ.0) RETURN



      call init_plot(iref,graph,iopt) 






c
C-----------------------------------------------------------------------
c
c     800 series plots - starting off with leakage barcharts!
c
C-----------------------------------------------------------------------
c
      if (iref.eq.801) then
c
c        Plots a sectional bar chart showing source
c        from inner and outer targets, pp wall and main wall due
c        to ion and atom fluxes.
c
c
c        Set up the framework
c
         REF='Summary of Source and Leakage'
         PLANE=' '
c
         CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >      YYMIN,YYMAX,TABLE,XLAB,YLAB,3,SMOOTH,1,ANLY,3)

c
c        This is the Source and Leakage summary bar chart
c
c        Plot Source First
c
c
c        Set title
c
         ylab = 'Percent %'
c
c        Set up labels
c
         pnames1(1) = INNER//' Target'
         pnames1(2) = OUTER//' Target'
         pnames1(3) = INNER//' Target'
         pnames1(4) = OUTER//' Target'
         pnames1(5) = ' '
         pnames1(6) = 'PP Wall'
         pnames1(7) = 'Main Wall'
c
         pnames2(1) = 'Ion Impact'
         pnames2(2) = 'Ion Impact'
         pnames2(3) = 'Atom Impact'
         pnames2(4) = 'Atom Impact'
         pnames2(5) = ' '
         pnames2(6) = 'Atom Impact'
         pnames2(7) = 'Atom Impact'
c
         cnames(1)  = 'Physical'
         cnames(2)  = 'Chemical'
         cnames(3)  = 'Self'
c
c        Set parameters
c
         nosets = maxsets
         novals = maxpts
c
         istart= 1
         istop = 7
c
c        Assign values to the valsts array and
c        calculate the min and max values for the plot
c
c        Copy over target sources
c
         totsrc = (targsrc(3,4) + wallsrc(5,3))/100.0
c
c        Check for valid totsrc before plotting
c
         if (totsrc.le.0.0) then 
            write(6,'(a,i6,g12.5)') 
     >             'ERROR:TOTSRC IS INVALID IN PLOT:',iref,totsrc
            write(0,'(a,i6,g12.5)') 
     >             'ERROR:TOTSRC IS INVALID IN PLOT:',iref,totsrc
            return
         endif 
c
         write(grtitle,'(a,1x,f10.3,1x,a)') 'Total Source:',
     >           totsrc*100.0,'Particles'
c
         call rzero (valsts,maxpts*maxsets)
c
         ymin = 0.0
c
         do in = 1,3
            valsts(1,in) = targsrc(1,in)/totsrc
            valsts(2,in) = targsrc(2,in)/totsrc
         end do
c
         ymax = max(targsrc(1,4)/totsrc,targsrc(2,4)/totsrc)
c
c        Copy over wallsources
c
         do in = 1,2
            valsts(in+2,1) = wallsrc(in,1)/totsrc
            valsts(in+2,2) = wallsrc(in,2)/totsrc
            valsts(in+2,3) = 0.0
            ymax = max (ymax,wallsrc(in,3)/totsrc)
         end do
c
c        Insert empty entry
c
         valsts(5,1) = 0.0
         valsts(5,2) = 0.0
         valsts(5,3) = 0.0
c
         do in = 3,4
            valsts(in+3,1) = wallsrc(in,1)/totsrc
            valsts(in+3,2) = wallsrc(in,2)/totsrc
            valsts(in+3,3) = 0.0
            ymax = max (ymax,wallsrc(in,3)/totsrc)
         end do
c
c        Plot bar chart
c
         call grbar(valsts,istart,istop,novals,nosets,
     >              ymin,ymax,1,pnames1,pnames2,grtitle,cnames,
     >              ylab)
c
c
c        Plot Leakage Next
c
c        Keep labels from previous half
c
c        Set title
c
c         ylab = 'Percent %'
c
c        Set up labels
c
c         pnames1(1) = 'Inner Target'
c         pnames1(2) = 'Outer Target'
c         pnames1(3) = 'Inner Target'
c         pnames1(4) = 'Outer Target'
c         pnames1(5) = ' '
c         pnames1(6) = 'PP Wall'
c         pnames1(7) = 'Main Wall'
c
c         pnames2(1) = 'Ion Impact'
c         pnames2(2) = 'Ion Impact'
c         pnames2(3) = 'Atom Impact'
c         pnames2(4) = 'Atom Impact'
c         pnames2(5) = ' '
c         pnames2(6) = 'Atom Impact'
c         pnames2(7) = 'Atom Impact'
c
c        Zero out cnames so that teh symbol table
c        entries are not done again
c
         cnames(1)  = ' '
         cnames(2)  = ' '
         cnames(3)  = ' '
c
c        Set parameters
c
         nosets = maxsets
         novals = maxpts
c
         istart= 1
         istop = 7
c
c        Assign values to the valsts array and
c        calculate the min and max values for the plot
c
c        Copy over target sources
c
         totleak = (targleak(3,4) + wallleak(5,3)) /100.0
c
c        Check for valid totleak before plotting
c
         if (totleak.le.0.0) then 
c
c           Set totleak = 1.0 to allow plot to finish
c             
            totleak = 1.0
c
         endif 
c
         write(grtitle,'(a,1x,f10.3,1x,a,1x,f7.2,a)') 'Total Leakage:',
     >       totleak*100.0,'Particles or',totleak/totsrc*100.0,'%'
c
         call rzero (valsts,maxpts*maxsets)
c
         ymin = 0.0
c
         do in = 1,3
            valsts(1,in) = targleak(1,in)/totleak
            valsts(2,in) = targleak(2,in)/totleak
         end do
c
         ymax = max(targleak(1,4)/totleak,targleak(2,4)/totleak)
c
c        Copy over wallsources
c
         do in = 1,2
            valsts(in+2,1) = wallleak(in,1)/totleak
            valsts(in+2,2) = wallleak(in,2)/totleak
            valsts(in+2,3) = 0.0
            ymax = max (ymax,wallleak(in,3)/totleak)
         end do
c
c        Insert empty entry
c
         valsts(5,1) = 0.0
         valsts(5,2) = 0.0
         valsts(5,3) = 0.0
c
         do in = 3,4
            valsts(in+3,1) = wallleak(in,1)/totleak
            valsts(in+3,2) = wallleak(in,2)/totleak
            valsts(in+3,3) = 0.0
            ymax = max (ymax,wallleak(in,3)/totleak)
         end do


c
c        Plot bar chart
c
         call grbar(valsts,istart,istop,novals,nosets,
     >              ymin,ymax,2,pnames1,pnames2,grtitle,cnames,
     >              ylab)
c
c        Clean up - finish plot
c
         call frame
c
      endif
c
C-----------------------------------------------------------------------
c
c     803 - plot of all sources - wall and target
c
      if (iref.eq.803) then
c
c        This plots a sectional bar chart showing combined
c        sources of particles from every wall element.
c
c        Set up the framework
c
         REF='Detailed Source Plot'
         PLANE=' '
c
         CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >      YYMIN,YYMAX,TABLE,XLAB,YLAB,3,SMOOTH,1,ANLY,3)
c
c        Calculate title and totsrc
c
         totsrc = (targsrc(3,4) + wallsrc(5,3))/100.0
c
c        Check for valid totsrc before plotting
c
         if (totsrc.le.0.0) then 
            write(6,'(a,i6,g12.5)') 
     >             'ERROR:TOTSRC IS INVALID IN PLOT:',iref,totsrc
            write(0,'(a,i6,g12.5)') 
     >             'ERROR:TOTSRC IS INVALID IN PLOT:',iref,totsrc
            return
         endif 
c
         write(grtitle,'(a,1x,f10.3,1x,a)') 'Total Source:',
     >           totsrc*100.0,'Particles'
c
         call rzero(valsts,maxpts*maxsets)
c
         novals = maxpts
         nosets = maxsets
         istart = 1
         istop  = wallpts
c
         wlmax = 1
         wlzmax = wallpt(1,2)
c
         do in = 1,wallpts
c
            pnames1(in) = ' '
            pnames2(in) = ' '
c
            if (wallpt(in,2).gt.wlzmax) then
               wlzmax = wallpt(in,2)
               wlmax = in
            endif
c
            do ir = 1,nrs
c
c              Physically sputtered sources - atom
c
               valsts(in,1) =  valsts(in,1)
     >                       + wtsource(in,ir,1,4)
c
c              Chemically sputtered sources - atom
c
               valsts(in,2) =  valsts(in,2)
     >                       + wtsource(in,ir,1,5)
c
c              Access corresponding target index for data
c              if there is one
c
               if (wallpt(in,18).gt.0) then
c
                  ind = wallpt(in,18)
c
c                 Physically sputtered sources - ion
c
                  valsts(in,1) =  valsts(in,1)
     >                       + wtsource(ind,ir,1,1)
c
c                 Chemical sputtered sources - ion
c
                  valsts(in,2) =  valsts(in,2)
     >                       + wtsource(ind,ir,1,2)
c
c                 Self Sputtered Sources - ion
c
                  valsts(in,3) =  valsts(in,3)
     >                       + wtsource(ind,ir,1,3)
c
               endif
c
            end do
c
         end do
c
c        Renormalize to percentage basis and find ymax
c
         ymin = 0.0
         ymax = lo
c
         do in = 1,wallpts
            tmpsum = 0.0
            do is = 1,nosets
               valsts(in,is) = valsts(in,is) / totsrc
               tmpsum = tmpsum + valsts(in,is)
            end do
            ymax = max(ymax,tmpsum)
         end do

c
         call setup_pnames(pnames1,pnames2,wlmax,0)

c
         cnames(1)  = 'Physical'
         cnames(2)  = 'Chemical'
         cnames(3)  = 'Self'
c
         ylab = 'Percent %'
c
c
c        Plot bar chart
c
         call grbar(valsts,istart,istop,novals,nosets,
     >              ymin,ymax,3,pnames1,pnames2,grtitle,cnames,
     >              ylab)
c
c
c        Clean up - finish plot
c
         call frame
c
      endif

c
C-----------------------------------------------------------------------
c
c     805 - plot of all leakage - wall and target
c
      if (iref.eq.805) then
c
c        NOTE: Core leakage may NOT include leakage by neutrals which 
c              are not ionized in the core. 
c
c
c        This plots a sectional bar chart showing combined
c        leakage of particles from every wall element.
c
c        Set up the framework
c
         if (iopt.eq.1) then 
            REF='Detailed Leakage Plot'
         elseif (iopt.eq.2) then 
            REF='Detailed Leakage Probabiltiy Plot'
         endif

         PLANE=' '
c
         CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >      YYMIN,YYMAX,TABLE,XLAB,YLAB,3,SMOOTH,1,ANLY,3)
c
c        Calculate title and totleak
c
         totsrc = (targsrc(3,4) + wallsrc(5,3))/100.0
c
c        Check for valid totsrc before plotting
c
         if (totsrc.le.0.0) then 
            write(6,'(a,i6,g12.5)') 
     >             'ERROR:TOTSRC IS INVALID IN PLOT:',iref,totsrc
            write(0,'(a,i6,g12.5)') 
     >             'ERROR:TOTSRC IS INVALID IN PLOT:',iref,totsrc
            return
         endif 
c
         totleak = (targleak(3,4) + wallleak(5,3)) /100.0
c
c        Check for valid totleak before plotting
c
         if (totleak.le.0.0) then 
            write(6,'(a,i6,g12.5)') 
     >             'ERROR:TOTLEAK IS INVALID IN PLOT:',iref,totleak
            write(0,'(a,i6,g12.5)') 
     >             'ERROR:TOTLEAK IS INVALID IN PLOT:',iref,totleak
            return
         endif 
c
c
         write(grtitle,'(a,1x,f10.3,1x,a,1x,f7.2,a)') 'Total Leakage:',
     >       totleak*100.0,'Particles or',totleak/totsrc*100.0,'%'
c
         call rzero(valsts,maxpts*maxsets)
c
         novals = maxpts
         nosets = maxsets
         istart = 1
         istop  = wallpts
c
         wlmax = 1
         wlzmax = wallpt(1,2)
c
         do in = 1,wallpts
c
            pnames1(in) = ' '
            pnames2(in) = ' '
c
            if (wallpt(in,2).gt.wlzmax) then
               wlzmax = wallpt(in,2)
               wlmax = in
            endif
c
            do ir = 1,nrs
c
c              Physically sputtered leakage - atom
c
               valsts(in,1) =  valsts(in,1)
     >                       + wtsource(in,ir,3,4)
c
c              Chemically sputtered leakage - atom
c
               valsts(in,2) =  valsts(in,2)
     >                       + wtsource(in,ir,3,5)
c
c              Access corresponding target index for data
c              if there is one
c
               if (wallpt(in,18).gt.0) then
c
                  ind = wallpt(in,18)
c
c                 Physically sputtered leakage - ion
c
                  valsts(in,1) =  valsts(in,1)
     >                       + wtsource(ind,ir,3,1)
c
c                 Chemical sputtered leakage - ion
c
                  valsts(in,2) =  valsts(in,2)
     >                       + wtsource(ind,ir,3,2)
c
c                 Self Sputtered leakage - ion
c
                  valsts(in,3) =  valsts(in,3)
     >                       + wtsource(ind,ir,3,3)
c
               endif
c
            end do
c
         end do
c
c        Renormalize to percentage basis and find ymax
c
         if (iopt.eq.1) then  
            ymin = 0.0
            ymax = lo
c
            do in = 1,wallpts
               tmpsum = 0.0
               do is = 1,nosets
                  valsts(in,is) = valsts(in,is) / totleak
                  tmpsum = tmpsum + valsts(in,is)
               end do
               ymax = max(ymax,tmpsum)
            end do
c
c        Calcualte leakage probability percent for each wall segment.
c
         elseif (iopt.eq.2) then  

            if (cgrprint.ne.0) then 
c
               write(6,'(a)') 'PLOT 805:'//
     >                 ' SUMMARY OF LEAKAGE PROBABILITY DATA'
c
            endif
c
            ymin = 0.0
            ymax = 100.0
c
            do is = 1,nosets
               do in = 1,wallpts
c
                  if (cgrprint.ne.0) then 
c
                     write(6,'(2i5,3(1x,f12.3),1x,g16.8)')
     >                      in,is,valsts(in,is),
     >                      wallse(in),wallse_i(in),
     >                      valsts(in,is)/wallse(in)
c
                  endif
c
                  if (wallse(in).gt.0.0) then 
                     valsts(in,is) = valsts(in,is)/wallse(in)
     >                            *100.0
                  else 
                     valsts(in,is) = 0.0   
                  endif   



               end do
            end do
         endif


c
         call setup_pnames(pnames1,pnames2,wlmax,0)


c-------------- jdemod - remove later
c
c        Set up Pnames and cnames labels
c
c         pnames1(wltrap1) = '|'
c         pnames1(wltrap2) = '|'
c         pnames1(wlwall1) = '|'
c         pnames1(wlwall2) = '|'
c
c         pnames1((wltrap1+wltrap2)/2) = 'PP'
c         pnames2((wltrap1+wltrap2)/2) = 'Wall'
c
c
c         pnames1(wlwall1+ INT(0.1*(wlwall2-wlwall1))) = Outer
c         pnames1(wlwall1+ INT(0.9*(wlwall2-wlwall1))) = Inner
c         pnames1((wlwall1+wlwall2)/2) = 'Top'
c         pnames2((wlwall1+wlwall2)/2) = 'Z (m)'
c
c        Mark Z-distances for some wall sections
c
c
c        Max of wall
c
c         pnames1(wlmax) = 'Top'
c
c         pnames2(wlmax) = 'Main Wall'
c
c         write(pnames2(wlmax),'(f5.1)') wallpt(wlmax,2)
c
c        0.02 - Outside
c
c         pnames1(wlwall1+ INT(0.02*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.02*(wlwall2-wlwall1))),'(a,f4.1,
c     >     a)') 'Z=',wallpt(wlwall1+ INT(0.02*(wlwall2-wlwall1)),2),'m'
c
c        0.2
c
c         pnames1(wlwall1+ INT(0.2*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.2*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.2*(wlwall2-wlwall1)),2)
c
c        0.4
c
c         pnames1(wlwall1+ INT(0.4*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.4*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.4*(wlwall2-wlwall1)),2)
c
c        0.6
c
c         pnames1(wlwall1+ INT(0.6*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.6*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.6*(wlwall2-wlwall1)),2)
c
c        0.8
c
c         pnames1(wlwall1+ INT(0.8*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.8*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.8*(wlwall2-wlwall1)),2)
c
c        0.98 - Inside
c
c         pnames1(wlwall1+ INT(0.98*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.98*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.98*(wlwall2-wlwall1)),2)
c
c         pnames2((wlwall1+wlwall2)/2) = 'Main Wall'
c
c         pnames1((wlwall2+wltrap1)/2) = INNER
c         pnames2((wlwall2+wltrap1)/2) = 'Target'
c         pnames1((wltrap2+wallpts)/2) = OUTER
c         pnames2((wltrap2+wallpts)/2) = 'Target'
c
c----------------


         cnames(1)  = 'Physical'
         cnames(2)  = 'Chemical'
         cnames(3)  = 'Self'
c
         if (iopt.eq.1) then  
            ylab = 'Leakage Fraction Percent %'
         elseif (iopt.eq.2) then 
            ylab = 'Leakage Probability Percent %'
         endif
c
c        Plot bar chart
c
c
         call grbar(valsts,istart,istop,novals,nosets,
     >              ymin,ymax,3,pnames1,pnames2,grtitle,cnames,
     >              ylab)
c
c
c        Clean up - finish plot
c
         call frame
c
      endif
c
C-----------------------------------------------------------------------
c
c     807 - combined wall source and leakage by element on one plot.
c
      if (iref.eq.807) then
c
c        This plots a sectional bar chart showing combined
c        sources of particles from every wall element.
c
c        IOPT = 1 plots the total source of particles that are ionized
c        IOPT = 2 plots the total source of all particles 
c
c        Set up the framework
c
         REF='Detailed Source and Leakage'
         PLANE=' '
c
         CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >      YYMIN,YYMAX,TABLE,XLAB,YLAB,3,SMOOTH,1,ANLY,3)
c
c        Calculate title and totsrc
c
         if (iopt.eq.1) then 

            totsrc = (targsrc(3,4) + wallsrc(5,3))/100.0

         elseif (iopt.eq.2) then

            totsrc = 0.0
c
c           Target source contributons
c 
            do in = 1,nds
               totsrc = totsrc + wtsource(in,nrs+1,1,1)  
               totsrc = totsrc + wtsource(in,nrs+1,1,2)  
               totsrc = totsrc + wtsource(in,nrs+1,1,3)  
            end do    
c
c           Wall source contributions
c
            do in = 1,wallpts
               totsrc = totsrc + wtsource(in,nrs+1,1,4)  
               totsrc = totsrc + wtsource(in,nrs+1,1,5)  
            end do    

            totsrc = totsrc/100.0 
  
         endif
c
c        Check for valid totsrc before plotting
c
         if (totsrc.le.0.0) then 
            write(6,'(a,i6,g12.5)') 
     >             'ERROR:TOTSRC IS INVALID IN PLOT:',iref,totsrc
            write(0,'(a,i6,g12.5)') 
     >             'ERROR:TOTSRC IS INVALID IN PLOT:',iref,totsrc
            return
         endif 
c
         if (iopt.eq.1) then 

            write(grtitle,'(a,1x,f10.3,1x,a)') 'Total Ion Source:',
     >           totsrc*100.0,'Particles'

         else if (iopt.eq.2) then   

            write(grtitle,'(a,1x,f10.3,1x,a)') 'Total Source:',
     >           totsrc*100.0,'Particles'

         endif 
c
         call rzero(valsts,maxpts*maxsets)
c
         novals = maxpts
         nosets = maxsets
         istart = 1
         istop  = wallpts
c
         wlmax = 1
         wlzmax = wallpt(1,2)
c
         do in = 1,wallpts
c
            pnames1(in) = ' '
            pnames2(in) = ' '
c
            if (wallpt(in,2).gt.wlzmax) then
               wlzmax = wallpt(in,2)
               wlmax = in
            endif
c
            if (iopt.eq.1) then 
               irstart = 1
               irend   = nrs
            elseif (iopt.eq.2) then
               irstart = nrs+1
               irend   = nrs+1
            endif  
c
            do ir = irstart,irend
c
c              Physically sputtered sources - atom
c	       
               valsts(in,1) =  valsts(in,1)
     >                    + wtsource(in,ir,1,4)
c	       
c              Chemically sputtered sources - atom
c	       
               valsts(in,2) =  valsts(in,2)
     >                    + wtsource(in,ir,1,5)
c	       
c              Access corresponding target index for data
c              if there is one
c	       
               if (wallpt(in,18).gt.0) then
c	       
                  ind = wallpt(in,18)
c	       
c                 Physically sputtered sources - ion
c	       
                  valsts(in,1) =  valsts(in,1)
     >                       + wtsource(ind,ir,1,1)
c
c                 Chemical sputtered sources - ion
c
                  valsts(in,2) =  valsts(in,2)
     >                       + wtsource(ind,ir,1,2)
c
c                 Self Sputtered Sources - ion
c
                  valsts(in,3) =  valsts(in,3)
     >                       + wtsource(ind,ir,1,3)
c
               endif
c
            end do
c
         end do
c
c        Renormalize to percentage basis and find ymax
c
         ymin = 0.0
         ymax = lo
c
         do in = 1,wallpts
            tmpsum = 0.0
            do is = 1,nosets
               valsts(in,is) = valsts(in,is) / totsrc
               tmpsum = tmpsum + valsts(in,is)
            end do
            ymax = max(ymax,tmpsum)
         end do

c
         call setup_pnames(pnames1,pnames2,wlmax,0)

c--------- jdemod - remove later
c
c        Set up Pnames and cnames labels
c
c         pnames1(wltrap1) = '|'
c         pnames1(wltrap2) = '|'
c         pnames1(wlwall1) = '|'
c         pnames1(wlwall2) = '|'
c
c         pnames1((wltrap1+wltrap2)/2) = 'PP'
c         pnames2((wltrap1+wltrap2)/2) = 'Wall'
c
c
c         pnames1(wlwall1+ INT(0.1*(wlwall2-wlwall1))) = Outer
c         pnames1(wlwall1+ INT(0.9*(wlwall2-wlwall1))) = Inner
c
c         pnames1((wlwall1+wlwall2)/2) = 'Top'
c         pnames2((wlwall1+wlwall2)/2) = 'Z (m)'
c
c        Mark Z-distances for some wall sections
c
c
c        Max of wall
c
c         pnames1(wlmax) = 'Top'
c
c         pnames2(wlmax) = 'Main Wall'
c
c         write(pnames2(wlmax),'(f5.1)') wallpt(wlmax,2)
c
c        0.02 - Outside
c
c         pnames1(wlwall1+ INT(0.02*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.02*(wlwall2-wlwall1))),'(a,f4.1,
c     >     a)') 'Z=',wallpt(wlwall1+ INT(0.02*(wlwall2-wlwall1)),2),'m'
c
c        0.2
c
c         pnames1(wlwall1+ INT(0.2*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.2*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.2*(wlwall2-wlwall1)),2)
c
c        0.4
c
c         pnames1(wlwall1+ INT(0.4*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.4*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.4*(wlwall2-wlwall1)),2)
c
c        0.6
c
c         pnames1(wlwall1+ INT(0.6*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.6*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.6*(wlwall2-wlwall1)),2)
c
c        0.8
c
c         pnames1(wlwall1+ INT(0.8*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.8*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.8*(wlwall2-wlwall1)),2)
c
c        0.98 - Inside
c
c         pnames1(wlwall1+ INT(0.98*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.98*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.98*(wlwall2-wlwall1)),2)
c
c         pnames2((wlwall1+wlwall2)/2) = 'Main Wall'
c
c         pnames1((wlwall2+wltrap1)/2) = INNER
c         pnames2((wlwall2+wltrap1)/2) = 'Target'
c         pnames1((wltrap2+wallpts)/2) = OUTER
c         pnames2((wltrap2+wallpts)/2) = 'Target'
c
c-------------------

         cnames(1)  = 'Physical'
         cnames(2)  = 'Chemical'
         cnames(3)  = 'Self'
c
         ylab = 'Percent %'
c
c
c        Plot bar chart
c
         call grbar(valsts,istart,istop,novals,nosets,
     >              ymin,ymax,1,pnames1,pnames2,grtitle,cnames,
     >              ylab)
c
c
c        This plots a sectional bar chart showing combined
c        leakage of particles from every wall element.
c
c        Calculate title and totleak
c
         totleak = (targleak(3,4) + wallleak(5,3)) /100.0
c
         write(grtitle,'(a,1x,f10.3,1x,a,1x,f7.2,a)') 'Total Leakage:',
     >       totleak*100.0,'Particles or',totleak/totsrc*100.0,'%'
c
c
c        Check for valid totleak before plotting
c
         if (totleak.le.0.0) then 
c
c           Set totleak = 1.0 to allow plot to finish
c             
            totleak = 1.0
c
         endif 
c
         call rzero(valsts,maxpts*maxsets)
c
         novals = maxpts
         nosets = maxsets
         istart = 1
         istop  = wallpts
c
         do in = 1,wallpts
c
c            Leave labels as in the first half
c
c            pnames1(in) = ' '
c            pnames2(in) = ' '
c
            do ir = 1,nrs
c
c              Physically sputtered leakage - atom
c
               valsts(in,1) =  valsts(in,1)
     >                       + wtsource(in,ir,3,4)
c
c              Chemically sputtered leakage - atom
c
               valsts(in,2) =  valsts(in,2)
     >                       + wtsource(in,ir,3,5)
c
c              Access corresponding target index for data
c              if there is one
c
               if (wallpt(in,18).gt.0) then
c
                  ind = wallpt(in,18)
c
c                 Physically sputtered leakage - ion
c
                  valsts(in,1) =  valsts(in,1)
     >                       + wtsource(ind,ir,3,1)
c
c                 Chemical sputtered leakage - ion
c
                  valsts(in,2) =  valsts(in,2)
     >                       + wtsource(ind,ir,3,2)
c
c                 Self Sputtered leakage - ion
c
                  valsts(in,3) =  valsts(in,3)
     >                       + wtsource(ind,ir,3,3)
c
               endif
c
            end do
c
         end do
c
c        Renormalize to percentage basis and find ymax
c
         ymin = 0.0
         ymax = lo
c
         do in = 1,wallpts
            tmpsum = 0.0
            do is = 1,nosets
               valsts(in,is) = valsts(in,is) / totleak
               tmpsum = tmpsum + valsts(in,is)
            end do
            ymax = max(ymax,tmpsum)
         end do
c
c        Labels are the same as in the first half
c
c        Change symbol table names so that they are not plotted again.
c
         cnames(1)  = ' '
         cnames(2)  = ' '
         cnames(3)  = ' '
c
         ylab = 'Percent %'
c
c
c        Plot bar chart
c
         call grbar(valsts,istart,istop,novals,nosets,
     >              ymin,ymax,2,pnames1,pnames2,grtitle,cnames,
     >              ylab)
c
c
c        Clean up - finish plot
c
         call frame
c
      endif
c
c
c
c
C-----------------------------------------------------------------------
c
c     809 - combined wall source and leakage by element on one plot.
c
c         - horizontal scale contains wall index numbers.
c
      if (iref.eq.809) then
c
c        This plots a sectional bar chart showing combined
c        sources of particles from every wall element.
c
c        Set up the framework
c
         REF='Detailed Source and Leakage'
         PLANE=' '
c
         CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >      YYMIN,YYMAX,TABLE,XLAB,YLAB,3,SMOOTH,1,ANLY,3)
c
c        Calculate title and totsrc
c
         totsrc = (targsrc(3,4) + wallsrc(5,3))/100.0
c
c
c        Check for valid totsrc before plotting
c
         if (totsrc.le.0.0) then 
            write(6,'(a,i6,g12.5)') 
     >             'ERROR:TOTSRC IS INVALID IN PLOT:',iref,totsrc
            write(0,'(a,i6,g12.5)') 
     >             'ERROR:TOTSRC IS INVALID IN PLOT:',iref,totsrc
            return
         endif 
c
         write(grtitle,'(a,1x,f10.3,1x,a)') 'Total Source:',
     >           totsrc*100.0,'Particles'
c
         call rzero(valsts,maxpts*maxsets)
c
         novals = maxpts
         nosets = maxsets
         istart = 1
         istop  = wallpts
c
         wlmax = 1
         wlzmax = wallpt(1,2)
c
         do in = 1,wallpts
c
            pnames1(in) = ' '
            pnames2(in) = ' '
c
            if (in.lt.wlwall2.and.((real(in)/5.0).eq.(in/5))) then
               write(pnames1(in),'(i3)') in
            endif
c
            if (wallpt(in,2).gt.wlzmax) then
               wlzmax = wallpt(in,2)
               wlmax = in
            endif
c
            do ir = 1,nrs
c
c              Physically sputtered sources - atom
c
               valsts(in,1) =  valsts(in,1)
     >                       + wtsource(in,ir,1,4)
c
c              Chemically sputtered sources - atom
c
               valsts(in,2) =  valsts(in,2)
     >                       + wtsource(in,ir,1,5)
c
c              Access corresponding target index for data
c              if there is one
c
               if (wallpt(in,18).gt.0) then
c
                  ind = wallpt(in,18)
c
c                 Physically sputtered sources - ion
c
                  valsts(in,1) =  valsts(in,1)
     >                       + wtsource(ind,ir,1,1)
c
c                 Chemical sputtered sources - ion
c
                  valsts(in,2) =  valsts(in,2)
     >                       + wtsource(ind,ir,1,2)
c
c                 Self Sputtered Sources - ion
c
                  valsts(in,3) =  valsts(in,3)
     >                       + wtsource(ind,ir,1,3)
c
               endif
c
            end do
c
         end do
c
c        Renormalize to percentage basis and find ymax
c
         ymin = 0.0
         ymax = lo
c
         do in = 1,wallpts
            tmpsum = 0.0
            do is = 1,nosets
               valsts(in,is) = valsts(in,is) / totsrc
               tmpsum = tmpsum + valsts(in,is)
            end do
            ymax = max(ymax,tmpsum)
         end do
c
         call setup_pnames(pnames1,pnames2,wlmax,0)

c--------- jdemod - remove later
c
c        Set up Pnames and cnames labels
c
c         pnames1(wltrap1) = '|'
c         pnames1(wltrap2) = '|'
c         pnames1(wlwall1) = '|'
c         pnames1(wlwall2) = '|'
c
c         pnames1((wltrap1+wltrap2)/2) = 'PP'
c         pnames2((wltrap1+wltrap2)/2) = 'Wall'
c
c
c         pnames2(wlwall1+ INT(0.1*(wlwall2-wlwall1))) = Outer
c         pnames2(wlwall1+ INT(0.9*(wlwall2-wlwall1))) = Inner
c
c        Max of wall
c
c         pnames2(wlmax) = 'Top'
c
c
c         pnames1((wlwall2+wltrap1)/2) = INNER
c         pnames2((wlwall2+wltrap1)/2) = 'Target'
c         pnames1((wltrap2+wallpts)/2) = OUTER
c         pnames2((wltrap2+wallpts)/2) = 'Target'
c
c
c---------------

         cnames(1)  = 'Physical'
         cnames(2)  = 'Chemical'
         cnames(3)  = 'Self'
c
         ylab = 'Percent %'
c
c
c        Plot bar chart
c
         call grbar(valsts,istart,istop,novals,nosets,
     >              ymin,ymax,1,pnames1,pnames2,grtitle,cnames,
     >              ylab)
c
c
c        This plots a sectional bar chart showing combined
c        leakage of particles from every wall element.
c
c        Calculate title and totleak
c
         totleak = (targleak(3,4) + wallleak(5,3)) /100.0
c
         write(grtitle,'(a,1x,f10.3,1x,a,1x,f7.2,a)') 'Total Leakage:',
     >       totleak*100.0,'Particles or',totleak/totsrc*100.0,'%'
c
c
c        Check for valid totleak before plotting
c
         if (totleak.le.0.0) then 
c
c           Set totleak = 1.0 to allow plot to finish
c             
            totleak = 1.0
c
         endif 
c
         call rzero(valsts,maxpts*maxsets)
c
         novals = maxpts
         nosets = maxsets
         istart = 1
         istop  = wallpts
c
         do in = 1,wallpts
c
c            Leave labels as in the first half
c
c            pnames1(in) = ' '
c            pnames2(in) = ' '
c
            do ir = 1,nrs
c
c              Physically sputtered leakage - atom
c
               valsts(in,1) =  valsts(in,1)
     >                       + wtsource(in,ir,3,4)
c
c              Chemically sputtered leakage - atom
c
               valsts(in,2) =  valsts(in,2)
     >                       + wtsource(in,ir,3,5)
c
c              Access corresponding target index for data
c              if there is one
c
               if (wallpt(in,18).gt.0) then
c
                  ind = wallpt(in,18)
c
c                 Physically sputtered leakage - ion
c
                  valsts(in,1) =  valsts(in,1)
     >                       + wtsource(ind,ir,3,1)
c
c                 Chemical sputtered leakage - ion
c
                  valsts(in,2) =  valsts(in,2)
     >                       + wtsource(ind,ir,3,2)
c
c                 Self Sputtered leakage - ion
c
                  valsts(in,3) =  valsts(in,3)
     >                       + wtsource(ind,ir,3,3)
c
               endif
c
            end do
c
         end do
c
c        Renormalize to percentage basis and find ymax
c
         ymin = 0.0
         ymax = lo
c
         do in = 1,wallpts
            tmpsum = 0.0
            do is = 1,nosets
               valsts(in,is) = valsts(in,is) / totleak
               tmpsum = tmpsum + valsts(in,is)
            end do
            ymax = max(ymax,tmpsum)
         end do
c
c        Labels are the same as in the first half
c
c        Change symbol table names so that they are not plotted again.
c
         cnames(1)  = ' '
         cnames(2)  = ' '
         cnames(3)  = ' '
c
         ylab = 'Percent %'
c
c
c        Plot bar chart
c
         call grbar(valsts,istart,istop,novals,nosets,
     >              ymin,ymax,2,pnames1,pnames2,grtitle,cnames,
     >              ylab)
c
c
c        Clean up - finish plot
c
         call frame
c
      endif
c
C-----------------------------------------------------------------------
c
c     811 - plot of all deposition - ion and neutral
c
      if (iref.eq.811) then
c
c        This plots a sectional bar chart showing combined
c        ion and neutral deposition of particles on every wall element.
c
c        IOPT = 1 = combined Ion + Neutral Deposition
c             = 2 = ION Deposition
c             = 3 = NEUTRAL Deposition
c             = 4 = DIVERTOR LEAKED ION Deposition
c
c         call print_deposition(57,cgridopt)
c
         if (iopt.lt.1.or.iopt.gt.4) then 
            write(6,'(a,i4)') 
     >           'Plot 811: IOPT value is invalid (not 1,2,3):',iopt
            write(0,'(a,i4)') 
     >           'Plot 811: IOPT value is invalid (not 1,2,3):',iopt
            return
         endif  
c
c        Set up the framework
c
         if (iopt.eq.1) then 
            REF='Detailed Deposition Plot'
         elseif (iopt.eq.2) then
            REF='Detailed ION Deposition Plot'
         elseif (iopt.eq.3) then 
            REF='Detailed NEUTRAL Deposition Plot'
         elseif (iopt.eq.3) then 
            REF='Detailed ION LEAKAGE Deposition Plot'
         endif

         PLANE=' '
c
         CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >      YYMIN,YYMAX,TABLE,XLAB,YLAB,3,SMOOTH,1,ANLY,3)
c
c        Calculate title and plot scale factor
c
         if (iopt.eq.1) then 
            totsrc = wallsn(maxpts+1) + wallsi(maxpts+1)
         elseif (iopt.eq.2) then 
            totsrc = wallsi(maxpts+1)
         elseif (iopt.eq.3) then 
            totsrc = wallsn(maxpts+1)
         elseif (iopt.eq.2) then 
            totsrc = wallsil(maxpts+1)
         endif  
c
c        Check for valid totsrc before plotting
c
         if (totsrc.le.0.0) then 
            write(6,'(a,i6,g12.5)') 
     >             'ERROR:TOTSRC IS INVALID IN PLOT:',iref,totsrc
            write(0,'(a,i6,g12.5)') 
     >             'ERROR:TOTSRC IS INVALID IN PLOT:',iref,totsrc
            return
         endif 
c
c
         if (iopt.eq.1) then 
            write(grtitle,'(a,1x,f10.3,1x,a,1x,f7.2,a)')
     >          'Total Deposition:',
     >          totsrc,'Particles'
         elseif (iopt.eq.2) then 
            write(grtitle,'(a,1x,f10.3,1x,a,1x,f7.2,a)')
     >          'Total ION Deposition:',
     >          totsrc,'Particles'
         elseif (iopt.eq.3) then 
            write(grtitle,'(a,1x,f10.3,1x,a,1x,f7.2,a)')
     >          'Total NEUT Deposition:',
     >          totsrc,'Particles'
         elseif (iopt.eq.2) then 
            write(grtitle,'(a,1x,f10.3,1x,a,1x,f7.2,a)')
     >          'Total LEAKED ION Dep:',
     >          totsrc,'Particles'
         endif

c
         call rzero(valsts,maxpts*maxsets)
c
         novals = maxpts
c
         if (iopt.eq.1) then 
            nosets = 2
         else
            nosets = 1
         endif 
c
         istart = 1
         istop  = wallpts
c
         wlmax = 1
         wlzmax = wallpt(1,2)
c
         write (6,*) 'Wall deposition:'
         do in = 1,wallpts
            write (6,'(i5,3(1x,g12.5))') in,wallsn(in),
     >                                wallsi(in),wallsil(in)
         end do
         write (6,*) 'Wall refs:',wltrap1,wltrap2,
     >               wlwall1,wlwall2
c
         do in = 1,wallpts
c
            pnames1(in) = ' '
            pnames2(in) = ' '
c
            if (wallpt(in,2).gt.wlzmax) then
               wlzmax = wallpt(in,2)
               wlmax = in
            endif

            if (iopt.eq.1) then 
c
c              Ion deposition
c
               valsts(in,1) =  valsts(in,1)
     >                       + wallsi(in)
c
c              Neutral depostition
c
               valsts(in,2) =  valsts(in,2)
     >                       + wallsn(in)
c
            elseif (iopt.eq.2) then 
c
c              Ion deposition
c
               valsts(in,1) =  valsts(in,1)
     >                       + wallsi(in)
c
            elseif (iopt.eq.3) then 
c
c              Neutral depostition
c
               valsts(in,1) =  valsts(in,1)
     >                       + wallsn(in)
c
            elseif (iopt.eq.4) then 
c
c              LEAKED Ion deposition
c
               valsts(in,1) =  valsts(in,1)
     >                       + wallsil(in)
c
            endif


         end do
c
c        Renormalize to percentage basis and find ymax
c
         ymin = 0.0
         ymax = lo
c
         do in = 1,wallpts
            tmpsum = 0.0
            do is = 1,nosets
               valsts(in,is) = valsts(in,is) / totsrc
               tmpsum = tmpsum + valsts(in,is)
            end do
            ymax = max(ymax,tmpsum)
         end do
c
c        Set up Pnames and cnames labels
c
         call setup_pnames(pnames1,pnames2,wlmax,0)

c-----  jdemod - remove later
c     
c
c         if (cgridopt.ne.RIBBON_GRID) then 
c            pnames1(wltrap1) = '|'
c            pnames1(wltrap2) = '|'
c            pnames1(wlwall1) = '|'
c            pnames1(wlwall2) = '|'
c
c            pnames1((wltrap1+wltrap2)/2) = 'PP'
c            pnames2((wltrap1+wltrap2)/2) = 'Wall'
c         endif
c
c
c         pnames1(wlwall1+ INT(0.1*(wlwall2-wlwall1))) = Outer
c         pnames1(wlwall1+ INT(0.9*(wlwall2-wlwall1))) = Inner
c         pnames1((wlwall1+wlwall2)/2) = 'Top'
c         pnames2((wlwall1+wlwall2)/2) = 'Z (m)'
c
c        Mark Z-distances for some wall sections
c
c
c        Max of wall
c
c         pnames1(wlmax) = 'Top'
c
c         pnames2(wlmax) = 'Main Wall'
c
c         write(pnames2(wlmax),'(f5.1)') wallpt(wlmax,2)
c
c        0.02 - Outside
c
c         pnames1(wlwall1+ INT(0.02*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.02*(wlwall2-wlwall1))),'(a,f4.1,
c     >     a)') 'Z=',wallpt(wlwall1+ INT(0.02*(wlwall2-wlwall1)),2),'m'
c
c        0.2
c
c         pnames1(wlwall1+ INT(0.2*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.2*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.2*(wlwall2-wlwall1)),2)
c
c        0.4
c
c         pnames1(wlwall1+ INT(0.4*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.4*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.4*(wlwall2-wlwall1)),2)
c
c        0.6
c
c         pnames1(wlwall1+ INT(0.6*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.6*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.6*(wlwall2-wlwall1)),2)
c
c        0.8
c
c         pnames1(wlwall1+ INT(0.8*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.8*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.8*(wlwall2-wlwall1)),2)
c
c        0.98 - Inside
c
c         pnames1(wlwall1+ INT(0.98*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.98*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.98*(wlwall2-wlwall1)),2)
c
c         pnames2((wlwall1+wlwall2)/2) = 'Main Wall'
c
c         pnames1((wlwall2+wltrap1)/2) = INNER
c         pnames2((wlwall2+wltrap1)/2) = 'Target'
c         pnames1((wltrap2+wallpts)/2) = OUTER
c         pnames2((wltrap2+wallpts)/2) = 'Target'
c
c-------------------------

         if (iopt.eq.1) then 
            cnames(1)  = 'Ion  Dep'
            cnames(2)  = 'Neut Dep'
         elseif (iopt.eq.2) then 
            cnames(1)  = 'Ion  Dep'
         elseif (iopt.eq.3) then
            cnames(1)  = 'Neut Dep'
         elseif (iopt.eq.4) then 
            cnames(1)  = 'L-IonDep'
         endif
c
         ylab = 'Fraction of Total'
c
c
c        Plot bar chart
c
         call grbar(valsts,istart,istop,novals,nosets,
     >              ymin,ymax,3,pnames1,pnames2,grtitle,cnames,
     >              ylab)
c
c
c        Clean up - finish plot
c
         call frame
c
      endif


c
C-----------------------------------------------------------------------
c
c     813 - plot of all erosion
c
      if (iref.eq.813) then
c
c        This plots a sectional bar chart showing 
c        erosion of particles from every wall element.
c
c        Set up the framework
c
         REF='Detailed Erosion Plot'
         PLANE=' '
c
         CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >      YYMIN,YYMAX,TABLE,XLAB,YLAB,3,SMOOTH,1,ANLY,3)
c
c        Calculate title and totleak
c
         totsrc = wallse(maxpts+1)
c
c        Check for valid totsrc before plotting
c
         if (totsrc.le.0.0) then 
            write(6,'(a,i6,g12.5)') 
     >             'ERROR:TOTSRC IS INVALID IN PLOT:',iref,totsrc
            write(0,'(a,i6,g12.5)') 
     >             'ERROR:TOTSRC IS INVALID IN PLOT:',iref,totsrc
            return
         endif 
c
c
         write(grtitle,'(a,1x,f10.3,1x,a,1x,f7.2,a)') 'Total Erosion:',
     >       totsrc,'Particles'
c
         call rzero(valsts,maxpts*maxsets)
c
         novals = maxpts
         nosets = 1
         istart = 1
         istop  = wallpts
c
         wlmax = 1
         wlzmax = wallpt(1,2)
c
         do in = 1,wallpts
c
            pnames1(in) = ' '
            pnames2(in) = ' '
c
            if (wallpt(in,2).gt.wlzmax) then
               wlzmax = wallpt(in,2)
               wlmax = in
            endif
c
c           Erosion
c
            valsts(in,1) =  valsts(in,1)
     >                    + wallse(in)
c
         end do
c
c        Renormalize to percentage basis and find ymax
c
         ymin = 0.0
         ymax = lo
c
         do in = 1,wallpts
            tmpsum = 0.0
            do is = 1,nosets
               valsts(in,is) = valsts(in,is) / totsrc
               tmpsum = tmpsum + valsts(in,is)
            end do
            ymax = max(ymax,tmpsum)
         end do
c
c        Set up Pnames and cnames labels
c

         call setup_pnames(pnames1,pnames2,wlmax,0)

c------  jdemod - remove later
c
c         if (cgridopt.ne.RIBBON_GRID) then 
c            pnames1(wltrap1) = '|'
c            pnames1(wltrap2) = '|'
c            pnames1(wlwall1) = '|'
c            pnames1(wlwall2) = '|'
c
c            pnames1((wltrap1+wltrap2)/2) = 'PP'
c            pnames2((wltrap1+wltrap2)/2) = 'Wall'
c
c         endif
c
c
c         pnames1(wlwall1+ INT(0.1*(wlwall2-wlwall1))) = Outer
c         pnames1(wlwall1+ INT(0.9*(wlwall2-wlwall1))) = Inner
c         pnames1((wlwall1+wlwall2)/2) = 'Top'
c         pnames2((wlwall1+wlwall2)/2) = 'Z (m)'
c
c        Mark Z-distances for some wall sections
c
c
c        Max of wall
c
c         pnames1(wlmax) = 'Top'
c
c         pnames2(wlmax) = 'Main Wall'
c
c         write(pnames2(wlmax),'(f5.1)') wallpt(wlmax,2)
c
c        0.02 - Outside
c
c         pnames1(wlwall1+ INT(0.02*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.02*(wlwall2-wlwall1))),'(a,f4.1,
c     >     a)') 'Z=',wallpt(wlwall1+ INT(0.02*(wlwall2-wlwall1)),2),'m'
c
c        0.2
c
c         pnames1(wlwall1+ INT(0.2*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.2*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.2*(wlwall2-wlwall1)),2)
c
c        0.4
c
c         pnames1(wlwall1+ INT(0.4*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.4*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.4*(wlwall2-wlwall1)),2)
c
c        0.6
c
c         pnames1(wlwall1+ INT(0.6*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.6*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.6*(wlwall2-wlwall1)),2)
c
c        0.8
c
c         pnames1(wlwall1+ INT(0.8*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.8*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.8*(wlwall2-wlwall1)),2)
c
c        0.98 - Inside
c
c         pnames1(wlwall1+ INT(0.98*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.98*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.98*(wlwall2-wlwall1)),2)
c
c         pnames2((wlwall1+wlwall2)/2) = 'Main Wall'
c
c         pnames1((wlwall2+wltrap1)/2) = INNER
c         pnames2((wlwall2+wltrap1)/2) = 'Target'
c         pnames1((wltrap2+wallpts)/2) = OUTER
c         pnames2((wltrap2+wallpts)/2) = 'Target'
c
c------------
c

         cnames(1)  = 'Erosion'
c
         ylab = 'Fraction of Total'
c
c
c        Plot bar chart
c
         call grbar(valsts,istart,istop,novals,nosets,
     >              ymin,ymax,3,pnames1,pnames2,grtitle,cnames,
     >              ylab)
c
c
c        Clean up - finish plot
c
         call frame
c
      endif

c
C-----------------------------------------------------------------------
c
c     815 - plot of net deposition - ion and neutral
c
      if (iref.eq.815) then
c
c        This plots a sectional bar chart showing combined
c        net deposition of particles from every wall element.
c
c        Set up the framework
c
         REF='Detailed Net-Deposition Plot'
         PLANE=' '
c
         CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >      YYMIN,YYMAX,TABLE,XLAB,YLAB,3,SMOOTH,1,ANLY,3)
c
c        Calculate title and totleak
c
         totsrc = wallse(maxpts+1)
c
         totdep = wallsi(maxpts+1) + wallsn(maxpts+1)
c
         write(grtitle,'(a,1x,f10.3,1x,a,1x,f10.3,a,f7.2)')
     >         'Total Erosion:',totsrc,' Total Dep:',totdep
c
         call rzero(valsts,maxpts*maxsets)
c
         novals = maxpts
         nosets = 1
         istart = 1
         istop  = wallpts
c
         wlmax = 1
         wlzmax = wallpt(1,2)
c
         do in = 1,wallpts
c
            pnames1(in) = ' '
            pnames2(in) = ' '
c
            if (wallpt(in,2).gt.wlzmax) then
               wlzmax = wallpt(in,2)
               wlmax = in
            endif
c
c           Net-Erosion
c
            valsts(in,1) =  wallsi(in)+wallsn(in)-wallse(in)
c
         end do
c
c        Renormalize to percentage basis and find ymax
c
         ymin =  hi
         ymax = -hi
c
         do in = 1,wallpts
            ymax = max(ymax,valsts(in,1))
            ymin = min(ymin,valsts(in,1))
         end do
c
c        Set up Pnames and cnames labels
c

         call setup_pnames(pnames1,pnames2,wlmax,0)

c------------- jdemod - remove later
c
c         if (cgridopt.ne.RIBBON_GRID) then 
c            pnames1(wltrap1) = '|'
c            pnames1(wltrap2) = '|'
c            pnames1(wlwall1) = '|'
c            pnames1(wlwall2) = '|'
c
c            pnames1((wltrap1+wltrap2)/2) = 'PP'
c            pnames2((wltrap1+wltrap2)/2) = 'Wall'
c         endif
c
c
c
c         pnames1(wlwall1+ INT(0.1*(wlwall2-wlwall1))) = Outer
c         pnames1(wlwall1+ INT(0.9*(wlwall2-wlwall1))) = Inner
c         pnames1((wlwall1+wlwall2)/2) = 'Top'
c         pnames2((wlwall1+wlwall2)/2) = 'Z (m)'
c
c        Mark Z-distances for some wall sections
c
c
c        Max of wall
c
c         pnames1(wlmax) = 'Top'
c
c         pnames2(wlmax) = 'Main Wall'
c
c         write(pnames2(wlmax),'(f5.1)') wallpt(wlmax,2)
c
c        0.02 - Outside
c
c         pnames1(wlwall1+ INT(0.02*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.02*(wlwall2-wlwall1))),'(a,f4.1,
c     >     a)') 'Z=',wallpt(wlwall1+ INT(0.02*(wlwall2-wlwall1)),2),'m'
c
c        0.2
c
c         pnames1(wlwall1+ INT(0.2*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.2*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.2*(wlwall2-wlwall1)),2)
c
c        0.4
c
c         pnames1(wlwall1+ INT(0.4*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.4*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.4*(wlwall2-wlwall1)),2)
c
c        0.6
c
c         pnames1(wlwall1+ INT(0.6*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.6*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.6*(wlwall2-wlwall1)),2)
c
c        0.8
c
c         pnames1(wlwall1+ INT(0.8*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.8*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.8*(wlwall2-wlwall1)),2)
c
c        0.98 - Inside
c
c         pnames1(wlwall1+ INT(0.98*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.98*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.98*(wlwall2-wlwall1)),2)
c
c         pnames2((wlwall1+wlwall2)/2) = 'Main Wall'
c
c         pnames1((wlwall2+wltrap1)/2) = INNER
c         pnames2((wlwall2+wltrap1)/2) = 'Target'
c         pnames1((wltrap2+wallpts)/2) = OUTER
c         pnames2((wltrap2+wallpts)/2) = 'Target'
c
c--------------

         cnames(1)  = 'Net Deposition'
c
         ylab = 'Net Number Deposited'
c
c
c        Plot bar chart
c
         call grbar(valsts,istart,istop,novals,nosets,
     >              ymin,ymax,4,pnames1,pnames2,grtitle,cnames,
     >              ylab)
c
c
c        Clean up - finish plot
c
         call frame
c
      endif


c
C-----------------------------------------------------------------------
c
c     817 - plot of particle endpoint distribution probability
c
      if (iref.eq.817) then
c
c        This plots a sectional bar chart showing the probability
c        of particles from each wall element ending in the specified region.
c
c        The value of IOPT defines the region 
c
c        IOPT = 1 = outer target     (total source)
c               2 = inner target     (total source) 
c               3 = main vessel wall (total source)  
c               4 = PFZ wall         (total source)
c        IOPT = 5 = outer target     (ion source)
c               6 = inner target     (ion source) 
c               7 = main vessel wall (ion source)  
c               8 = PFZ wall         (ion source)
c
c
c        Set labels and element integration ranges
c
         if (iopt.eq.1.or.iopt.eq.5) then 
            plane = 'OUTER TARGET' 
c
c           If the xpoint is up then the second target is the outer
c
            if (xpoint_up) then 
               iwstart = wltrap2+1
               iwend   = wallpts 
            else
               iwstart = wlwall2+1
               iwend   = wltrap1-1 
            endif 
c
         elseif (iopt.eq.2.or.iopt.eq.6) then 
            plane = 'INNER TARGET' 
c
c           If Xpoint is up then the first target is the inner
c
            if (xpoint_up) then 
               iwstart = wlwall2+1
               iwend   = wltrap1-1 
            else
               iwstart = wltrap2+1
               iwend   = wallpts 
            endif 
c
         elseif (iopt.eq.3.or.iopt.eq.7) then 
            plane = 'MAIN WALL' 
            iwstart = wlwall1
            iwend   = wlwall2
         elseif (iopt.eq.4.or.iopt.eq.8) then 
            plane = 'PFZ WALL' 
            iwstart = wltrap1
            iwend   = wltrap2  
         else
            return
         endif
c
c        Set up the framework
c
         if (iopt.ge.1.and.iopt.le.4) then 
            REF='TOTAL DEPOSITION PROBABILITY'
         elseif (iopt.ge.5.and.iopt.le.8) then 
            REF='ION DEPOSITION PROBABILITY'
         endif 
c
         CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >      YYMIN,YYMAX,TABLE,XLAB,YLAB,3,SMOOTH,1,ANLY,3)
c
         call rzero(valsts,maxpts*maxsets)
c
         novals = maxpts
         nosets = 1
         istart = 1
         istop  = wallpts
c
         wlmax = 1
         wlzmax = wallpt(1,2)
c
         do in = 1,wallpts
c
            pnames1(in) = ' '
            pnames2(in) = ' '
c
            if (wallpt(in,2).gt.wlzmax) then
               wlzmax = wallpt(in,2)
               wlmax = in
            endif
c
c           Caclulate number of particles from each wall element that end in the 
c           desired region
c
            do iw = iwstart,iwend

               if (iopt.ge.1.and.iopt.le.4) then 
                  valsts(in,1) =  valsts(in,1) + wtdep(in,iw,3)
               elseif (iopt.ge.5.and.iopt.le.8) then 
                  valsts(in,1) =  valsts(in,1) + wtdep(in,iw,1)
               endif                

            end do
c
         end do
c
c        Renormalize to percentage basis and find ymax
c
c        For Deposition probability plots - set the minimum and maximum to 0.0,1.0  
c         
         ymin = 0.0
         ymax = 1.0 
c
c         ymin =  hi
c         ymax = -hi
c
         totdep = 0.0 
c
         do in = 1,wallpts
c
c           Sum up total number of particles deposited in region. 
c
            totdep = totdep + valsts(in,1) 
c
            if (iopt.ge.1.and.iopt.le.4) then 
c
               if (wallse(in).gt.0) then   
                  valsts(in,1) = valsts(in,1) / wallse(in)
               else
                  valsts(in,1) = 0.0
               endif
c
            elseif (iopt.ge.5.and.iopt.le.8) then 
c
               if (wallse_i(in).gt.0) then   
                  valsts(in,1) = valsts(in,1) / wallse_i(in)
               else
                  valsts(in,1) = 0.0
               endif
c      
            endif
c
c            ymax = max(ymax,valsts(in,1))
c            ymin = min(ymin,valsts(in,1))
c
         end do
c
c        Calculate title 
c
         if (iopt.ge.1.and.iopt.le.4) then 
         
            totsrc = wallse(maxpts+1)
            write(grtitle,'(a,1x,f10.3,1x,a,1x,f10.3,a,f7.2)')
     >         'Total Erosion:',totsrc,
     >         ' Total Region Dep:',totdep

         elseif (iopt.ge.5.and.iopt.le.8) then

            totsrc = wallse_i(maxpts+1)
            write(grtitle,'(a,1x,f10.3,1x,a,1x,f10.3,a,f7.2)')
     >         'Total Ionized Erosion:',totsrc,
     >         ' Total Region Dep:',totdep

         endif
c
c        Set up Pnames and cnames labels
c
         call setup_pnames(pnames1,pnames2,wlmax,0)

c--------- jdemod - remove later
c
c         if (cgridopt.ne.RIBBON_GRID) then 
c            pnames1(wltrap1) = '|'
c            pnames1(wltrap2) = '|'
c            pnames1(wlwall1) = '|'
c            pnames1(wlwall2) = '|'
c
c            pnames1((wltrap1+wltrap2)/2) = 'PP'
c            pnames2((wltrap1+wltrap2)/2) = 'Wall'
c         endif
c
c
c         pnames1(wlwall1+ INT(0.1*(wlwall2-wlwall1))) = Outer
c         pnames1(wlwall1+ INT(0.9*(wlwall2-wlwall1))) = Inner
c
c        Mark Z-distances for some wall sections
c
c
c        Max of wall
c
c         pnames1(wlmax) = 'Top'
c
c         pnames2(wlmax) = 'Main Wall'
c
c         write(pnames2(wlmax),'(f5.1)') wallpt(wlmax,2)
c
c        0.02 - Outside
c
c         pnames1(wlwall1+ INT(0.02*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.02*(wlwall2-wlwall1))),'(a,f4.1,
c     >     a)') 'Z=',wallpt(wlwall1+ INT(0.02*(wlwall2-wlwall1)),2),'m'
c
c        0.2
c
c         pnames1(wlwall1+ INT(0.2*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.2*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.2*(wlwall2-wlwall1)),2)
c
c        0.4
c
c         pnames1(wlwall1+ INT(0.4*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.4*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.4*(wlwall2-wlwall1)),2)
c
c        0.6
c
c         pnames1(wlwall1+ INT(0.6*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.6*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.6*(wlwall2-wlwall1)),2)
c
c        0.8
c
c         pnames1(wlwall1+ INT(0.8*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.8*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.8*(wlwall2-wlwall1)),2)
c
c        0.98 - Inside
c
c         pnames1(wlwall1+ INT(0.98*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.98*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.98*(wlwall2-wlwall1)),2)
c
c         pnames2((wlwall1+wlwall2)/2) = 'Main Wall'
c
c         pnames1((wlwall2+wltrap1)/2) = INNER
c         pnames2((wlwall2+wltrap1)/2) = 'Target'
c         pnames1((wltrap2+wallpts)/2) = OUTER
c         pnames2((wltrap2+wallpts)/2) = 'Target'
c
c---------------

         cnames(1)  = 'Probablity'
c
         if (iopt.ge.1.and.iopt.le.4) then
            ylab = 'Region Deposition Probability'
         elseif (iopt.ge.5.and.iopt.le.8) then
            ylab = 'Region Ion Deposition Probability'
         endif 
c
c
c        Plot bar chart
c
         call grbar(valsts,istart,istop,novals,nosets,
     >              ymin,ymax,4,pnames1,pnames2,grtitle,cnames,
     >              ylab)
c
c
c        Clean up - finish plot
c
         call frame
c
      endif

c
C-----------------------------------------------------------------------
c
c     819 - plot of all deposition - ion and neutral
c           as a function of the distance along the wall.
c         - wall distance is initially calculated from the inner
c           midplane counter clockwise
c         - distances are "adjusted" to match experimental data (or
c           experimental data coordinates are adjusted to match the 
c           grid
c         - supports plotting of only a section of the wall based on
c           an input min and max value
c         - The deposition can be scaled by an input scaling factor for
c           cases where ABSFAC is not valid
c
      if (iref.eq.819) then

         call plot_deposition(iopt,cgridopt)


      endif




C
C-----------------------------------------------------------------------
c
c     821 - Plot of source across targets as a function of PSIN 
c
c
      IF (IREF.EQ.821) THEN
c
        ELABS(1) = 'A-P ATOM-PHYS'
        ELABS(2) = 'A-C ATOM-CHEM'
        ELABS(3) = 'I-P ION -PHYS'
        ELABS(4) = 'I-C ION -CHEM'
        ELABS(5) = 'SELFSELF '
        ELABS(6) = 'TOT TOTAL'
        XLAB = '   PSIN   '
        YLAB = '   FLUX (PART/M^2/S)   '

c
c
c       Scaling factor
c
c
c       Total weight of particles launched
c
        totsrc = (targsrc(3,4) + wallsrc(5,3))
c
        if (absfac.gt.0.0.and.totsrc.gt.0.0) then
c
c          Calculate scaling factor to convert actual particles
c          launched to a scaled fraction of total amount of 
c          impurity influx.
c          
c          Do not include major radius adjustments since actual
c          DIVIMP source is not corrected for this. 
c  
c           scalef_821 = absfac * 2.0 * PI * R0 / totsrc
c 
           scalef_821 = absfac / totsrc
c
        elseif (totsrc.gt.0.0) then 

           scalef_821 = 1.0/totsrc  
c
        else 
c
           scalef_821 = 1.0
c
        endif 

c
c       Data types in wtsource array (last index)
c       Note: Targ refers to ion sputtering which occurs only
c             on target elements while Wall refers to atom 
c             sputtering which occurs on the entire wall including
c             the target elements. The nomenclature can be 
c             confusing.
c       NOTE!: Types 1 to 3 are indexed in wtsource from 1 to nds
c              while types 4,5 are based on 1,wallpts. This is due to 
c              the different types of launches required for the 
c              different particles in neut. Could use some cleaning up.
c
c
c       Targ+phys = 1
c       Targ+chem = 2
c       Targ-self = 3
c       Wall+phys = 4
c       Wall+chem = 5
c       2Dneutral = 6
c       Refl Ion  = 7
c
c
c       First "INNER" target  
c
        REF  = INNER//' TARGET IMPURITY SOURCES'
c
        WRITE (IPLOT,9012) NPLOTS,REF
c
        call rzero(douts,maxnds+2)
        call rzero(dvals,(maxnds+2)*maxngs)
c
        do id = ndsin-1,2,-1
c
c          First target 
c
           ind = (ndsin-1) - id + 1 
c
           targid  = 1
           wallind = wallindex(id) 
c
           douts(ind) = psitarg(irds(id),targid)
           dwids(ind) = 1.0
c
           if (wallind.ne.0) then 
c
              do ir = 1,nrs
c
c                Physically sputtered sources - atom
c
                 dvals(ind,1) = dvals(ind,1) +
     >                       wtsource(wallind,ir,1,4)   
c
c                Chemically sputtered sources - atom
c
                 dvals(ind,2) = dvals(ind,2) +
     >                       wtsource(wallind,ir,1,5)   
c
c                Physically sputtered sources - ion
c
                 dvals(ind,3) = dvals(ind,3) +
     >                       wtsource(id,ir,1,1)   
c
c                Chemical sputtered sources - ion
c
                 dvals(ind,4) = dvals(ind,4) +
     >                       wtsource(id,ir,1,2)   
c
c                Self Sputtered Sources - ion
c
                 dvals(ind,5) = dvals(ind,5) +
     >                       wtsource(id,ir,1,3)   
c
              end do
c
c             Total
c
              do in = 1,5
c 
                 dvals(ind,6) = dvals(ind,6) + dvals(ind,in)
c
              end do
c
c             Scaling to absolute flux values - based on 
c             per meter toroidally per second values
c
c              geofact = dds(id) * 2.0 * PI * rp(id)
c
              geofact = dds(id) 
c
              dvals(ind,1) = dvals(ind,1) * scalef_821 / geofact
              dvals(ind,2) = dvals(ind,2) * scalef_821 / geofact
              dvals(ind,3) = dvals(ind,3) * scalef_821 / geofact
              dvals(ind,4) = dvals(ind,4) * scalef_821 / geofact
              dvals(ind,5) = dvals(ind,5) * scalef_821 / geofact
              dvals(ind,6) = dvals(ind,6) * scalef_821 / geofact
c
              write(6,'(a,3i4,7(1x,g12.5))') 'DVAL:',
     >            ind,id,wallind,douts(ind),
     >            (dvals(ind,in),in=1,6)

c
           endif  

        end do
c
c       Generate plot 
c
        CALL DRAW (DOUTS,DWIDS,DVALS,MAXNDS+2,NDSIN-2,ANLY,
     >    6,99,DOUTS(1),DOUTS(NDSIN-2),0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,2,1.0,0)
c
c
c
c       Second "OUTER" target  
c
c
c
        REF  = OUTER//' TARGET IMPURITY SOURCES'
c
        WRITE (IPLOT,9012) NPLOTS,REF
c
        call rzero(douts,maxnds+2)
        call rzero(dvals,(maxnds+2)*maxngs)
c
        do id = ndsin+2,nds-1
c
c          Second target 
c
           ind = id - (ndsin+1)
c
           targid  = 2
           wallind = wallindex(id) 
c
           douts(ind) = psitarg(irds(id),targid)
           dwids(ind) = 1.0
c
           if (wallind.ne.0) then 
c
              do ir = 1,nrs
c
c                Physically sputtered sources - atom
c
                 dvals(ind,1) = dvals(ind,1) +
     >                       wtsource(wallind,ir,1,4)   
c
c                Chemically sputtered sources - atom
c
                 dvals(ind,2) = dvals(ind,2) +
     >                       wtsource(wallind,ir,1,5)   
c
c                Physically sputtered sources - ion
c
                 dvals(ind,3) = dvals(ind,3) +
     >                       wtsource(id,ir,1,1)   
c
c                Chemical sputtered sources - ion
c
                 dvals(ind,4) = dvals(ind,4) +
     >                       wtsource(id,ir,1,2)   
c
c                Self Sputtered Sources - ion
c
                 dvals(ind,5) = dvals(ind,5) +
     >                       wtsource(id,ir,1,3)   
c
              end do
c
c             Total
c
              do in = 1,5
c 
                 dvals(ind,6) = dvals(ind,6) + dvals(ind,in)
c
              end do
c
c             Scaling to absolute flux values - based on 
c             per meter toroidally per second values
c
c              geofact = dds(id) * 2.0 * PI * rp(id)
c
              geofact = dds(id)
c
              dvals(ind,1) = dvals(ind,1) * scalef_821 / geofact
              dvals(ind,2) = dvals(ind,2) * scalef_821 / geofact
              dvals(ind,3) = dvals(ind,3) * scalef_821 / geofact
              dvals(ind,4) = dvals(ind,4) * scalef_821 / geofact
              dvals(ind,5) = dvals(ind,5) * scalef_821 / geofact
              dvals(ind,6) = dvals(ind,6) * scalef_821 / geofact
c
              write(6,'(a,3i4,7(1x,g12.5))') 'DVAL:',
     >            ind,id,wallind,douts(ind),
     >            (dvals(ind,in),in=1,6)
c
           endif  

        end do
c
c       Generate plot 
c
        CALL DRAW (DOUTS,DWIDS,DVALS,MAXNDS+2,NDS-NDSIN-2,ANLY,
     >    6,99,DOUTS(1),DOUTS(NDS-ndsin-2),0.0,HI,IGNORS,
     >    ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,2,1.0,0)
c
      ENDIF
c
C-----------------------------------------------------------------------
c
c     831 - multicase plot of net deposition - ion and neutral
c
c     NOTE: This plot is currently customized for the C12/C13 
c           study - labels and titles would need changing for 
c           other cases however the plot should work correctly
c           regardless. 
c
c
      if (iref.eq.831) then
c
c        This plots a sectional bar chart showing combined
c        net deposition of particles from every wall element.
c
c
c        Load secondary input parameters
c
c        ISELECT = 1 - combined net erosion/deposition
c                  2 - combined erosion
c                  3 - combined deposition
c
c        ISCALE  = 1 - ABSFAC_neut from both
c                  2 - specified external scale factor (scalef)
c                  3 - auto scale unsscaled source based on
c                      inner target deposition/erosion
c
         call rdg_3i1r(graph,iselect,iscale,ival,scalef,ierr)
c
c        Get file name and execution command of second case to load
c
         call rdfn_cmd(graph,cmd,fn,iflag,ierr)
c
c        Allocate variables to store the results - or used fixed
c        size locals and check to make sure the size of the 
c        second case parameters isn't too large.
c      
         allocate (absfacs(2,2),stat=ierr) 
c
         if (ierr.ne.0) then 
            write(0,*) 'ERROR PLOT 831: Allocate for ABSFACS'//
     >                   ' data storage failed: RETURNING',ierr
            write(6,*) 'ERROR PLOT 831: Allocate for ABSFAS'//
     >                   ' data storage failed: RETURNING',ierr
            return
         endif 
c
         allocate (depdata(maxpts+1,4,2),stat=ierr)

c
         if (ierr.ne.0) then 
            write(0,*) 'ERROR PLOT 831: Allocate for DEPDATA'//
     >                   ' data storage failed: RETURNING',ierr
            write(6,*) 'ERROR PLOT 831: Allocate for DEPDATA'//
     >                   ' data storage failed: RETURNING',ierr
            return
         endif 
c
c        Assume that all cases have exactly the same wall 
c        geometry - comparion isn't meaningful otherwise - 
c        However - I can have some checks for this. 
c
         call load_erodep_data(absfacs,depdata,cmd,fn,1)
c
c        Calculate scaling factor for second case if that plot option is 
c        selected. (As much as possible make net ero/dep = 0.0 for worst inner
c        target main SOL element.
c
         plane = ' '
c
         if (iscale.eq.1) then 
            if (absfacs(1,1).eq.0.0.or.absfacs(1,2).eq.0.0) then 
               write(0,'(a,2(1x,g12.5))') 
     >                'ERROR: PLOT 831: Scale option 1: '//
     >                'ABSOLUTE factor for one or both cases is'//
     >                ' undefined:',absfacs(1,1),absfacs(1,2)          
               write(6,'(a,2(1x,g12.5))') 
     >                'ERROR: PLOT 831: Scale option 1: '//
     >                'ABSOLUTE factor for one or both cases is'//
     >                ' undefined:',absfacs(1,1),absfacs(1,2)          
               return
            endif 
c
            plane = 'Scale factors from case data'
c
         elseif (iscale.eq.2) then 
c
c           Load imposed scaling factor for all cases where absfac is undefined 
c
            if (absfacs(1,1).eq.0.0) then   
               absfacs(1,1) = scalef
            endif 
c
            if (absfacs(1,2).eq.0.0) then   
               absfacs(1,2) = scalef
            endif 
c
            plane = 'Scale factors assigned if undefined'
c
         elseif (iscale.eq.3) then 
c
c           Calculate scale factor for second (or undefined) dataset 
c           and assign it to the appropriate ABSFACS entry.
c
            call calc_erodep_scalef(absfacs,depdata,scalef,ierr)
c
            if (ierr.eq.0) then 
               plane = 'Second scale factor calculated'
            elseif (ierr.eq.1) then 
               plane = 'Second scale factor input'
            elseif (ierr.eq.2) then 
               plane = 'Second scale factor equals first'
            endif
c
         endif
c
c        Combine the data including the requisite scaling factors
c        - write both scaling factors to the graph
c
c
c        Set up the framework
c
         if (iselect.eq.1) then 
            REF='Detailed Combined Net-Deposition'
         elseif (iselect.eq.2) then 
            REF='Detailed Combined Erosion'
         elseif (iselect.eq.3) then 
            REF='Detailed Combined Deposition'
         endif  
c
c        Calculate title and totleak
c
         tot1 = absfacs(1,1)
         tot2 = absfacs(1,2)
c
         write(grtitle,'(a,1x,1p,e10.3,1x,a,1x,e10.3)')
     >         'Original C13 Puff:',tot2,
     >         ' Local C13 Sputtering:',tot1
c
c        Zero data array 
c
         call rzero(valsts,maxpts*maxsets)
c
         novals = maxpts
c
c        Set number of data sets depending on type of data
c        to be plotted. 
c
         if (iselect.eq.1) then 
            nosets = 1
         elseif (iselect.eq.2.or.iselect.eq.3) then
            nosets = 2
         endif
c
         istart = 1
         istop  = wallpts
c
         wlmax = 1
         wlzmax = wallpt(1,2)
c
         write(6,'(a,i5)') 'PLOT 831: Plot data: ISELECT =',iselect
c
         do in = 1,wallpts
c
            pnames1(in) = ' '
            pnames2(in) = ' '
c
            if (wallpt(in,2).gt.wlzmax) then
               wlzmax = wallpt(in,2)
               wlmax = in
            endif
c
c           Net-Erosion
c
            if (iselect.eq.1) then 
c
c              Combine the net erosion from the two cases 
c
               valsts(in,1) =  depdata(in,1,1) *absfacs(1,1)
     >                       + depdata(in,1,2) *absfacs(1,2)
c
            elseif (iselect.eq.2) then 
c
c              Erosion from the two cases
c
               valsts(in,1) =  depdata(in,2,1) *absfacs(1,1)
               valsts(in,2) =  depdata(in,2,2) *absfacs(1,2)
c
            elseif (iselect.eq.3) then 
c
c              Deposition from the two cases
c
               valsts(in,1) =  depdata(in,3,1) *absfacs(1,1)
               valsts(in,2) =  depdata(in,3,2) *absfacs(1,2)
c
            endif
c
c           Change units to particles/m2/s if option set
c           Divide by length of wall segment
c       
            if (ival.eq.1.and.wallpt(in,7).gt.0.0) then 
               valsts(in,1) = valsts(in,1) / wallpt(in,7) 
               valsts(in,2) = valsts(in,2) / wallpt(in,7) 
            endif
c
            write(6,'(i5,2(1x,g12.5))') in,valsts(in,1),
     >                 valsts(in,2)
c
         end do
c
c        Find Ymin and Ymax for plotting
c
         ymin =  hi
         ymax = -hi
c
         do in = 1,wallpts
            if (iselect.eq.1) then 
               ymax = max(ymax,valsts(in,1))
               ymin = min(ymin,valsts(in,1))
            elseif (iselect.eq.2.or.iselect.eq.3) then 
               tmpsum = 0.0
               do is = 1,nosets
                  tmpsum = tmpsum + valsts(in,is)
               end do
               ymax = max(ymax,tmpsum)
               ymin = min(ymin,tmpsum)
            endif 
         end do
c
c        Draw labels on the chart
c
         write(nview,'(a,1x,1p,g12.5)') 'Max. Flux=',ymax
c        
         CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >      YYMIN,YYMAX,TABLE,XLAB,YLAB,3,SMOOTH,1,ANLY,3)
c
c        Define labels
c
         if (iselect.eq.1) then 
            cnames(1)  = 'Combined Net'
         elseif (iselect.eq.2) then 
c
            cnames(1)  = 'Sputtered C13'
            cnames(2)  = 'Original C13'
c
         elseif (iselect.eq.3) then 
c
            cnames(1)  = 'Sputtered C13'
            cnames(2)  = 'Original C13'
c
         endif
c
c        Set up Pnames labels
c
         call setup_pnames(pnames1,pnames2,wlmax,0)

c----- jdemod - remove later
c
c         if (cgridopt.ne.RIBBON_GRID) then 
c            pnames1(wltrap1) = '|'
c            pnames1(wltrap2) = '|'
c            pnames1(wlwall1) = '|'
c            pnames1(wlwall2) = '|'
c
c            pnames1((wltrap1+wltrap2)/2) = 'PP'
c            pnames2((wltrap1+wltrap2)/2) = 'Wall'
c         endif
c
c
c         pnames1(wlwall1+ INT(0.1*(wlwall2-wlwall1))) = Outer
c         pnames1(wlwall1+ INT(0.9*(wlwall2-wlwall1))) = Inner
c         pnames1((wlwall1+wlwall2)/2) = 'Top'
c         pnames2((wlwall1+wlwall2)/2) = 'Z (m)'
c
c        Mark Z-distances for some wall sections
c
c
c        Max of wall
c
c         pnames1(wlmax) = 'Top'
c
c         pnames2(wlmax) = 'Main Wall'
c
c         write(pnames2(wlmax),'(f5.1)') wallpt(wlmax,2)
c
c        0.02 - Outside
c
c         pnames1(wlwall1+ INT(0.02*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.02*(wlwall2-wlwall1))),'(a,f4.1,
c     >     a)') 'Z=',wallpt(wlwall1+ INT(0.02*(wlwall2-wlwall1)),2),'m'
c
c        0.2
c
c         pnames1(wlwall1+ INT(0.2*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.2*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.2*(wlwall2-wlwall1)),2)
c
c        0.4
c
c         pnames1(wlwall1+ INT(0.4*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.4*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.4*(wlwall2-wlwall1)),2)
c
c        0.6
c
c         pnames1(wlwall1+ INT(0.6*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.6*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.6*(wlwall2-wlwall1)),2)
c
c        0.8
c
c         pnames1(wlwall1+ INT(0.8*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.8*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.8*(wlwall2-wlwall1)),2)
c
c        0.98 - Inside
c
c         pnames1(wlwall1+ INT(0.98*(wlwall2-wlwall1))) = '^'
c         write(pnames2(wlwall1+ INT(0.98*(wlwall2-wlwall1))),'(f4.1)')
c     >              wallpt(wlwall1+ INT(0.98*(wlwall2-wlwall1)),2)
c
c         pnames2((wlwall1+wlwall2)/2) = 'Main Wall'
c
c
c        Modify these labels to mark separatrix element
c
c         pnames1((wlwall2+wltrap1)/2) = INNER
c         pnames2((wlwall2+wltrap1)/2) = 'Target'
c         pnames1((wltrap2+wallpts)/2) = OUTER
c         pnames2((wltrap2+wallpts)/2) = 'Target'
c
c         pnames1(wallindex(idds(irsep,1))) = '^'
c         pnames2((wlwall2+wltrap1)/2) = INNER//' Target'
c         pnames1(wallindex(idds(irsep,2))) = '^'
c         pnames2((wltrap2+wallpts)/2) = OUTER//' Target'
c
c-----------


         if (ival.eq.1) then          
            ylab = 'Particles/m2/s'
         else   
            ylab = 'Particles/m-tor/s'
         endif
c
c
c        Plot bar chart
c
         call grbar(valsts,istart,istop,novals,nosets,
     >              ymin,ymax,4,pnames1,pnames2,grtitle,cnames,
     >              ylab)
c
c
c        Clean up - finish plot
c
         call frame
c
c
c        Cleanup storage allocation
c
         deallocate(depdata)
         deallocate(absfacs)

      endif



c
c     End of 800 series plots
c
      return 
c
c     Format statements    
c
 9012 FORMAT(1X,'PLOT',I3,4X,A)
c
      end
c
c
c
      subroutine load_erodep_data(absfacs,depdata,cmd,fn,nsets)
      implicit none
      include 'params'
      include 'dynam3'
      include 'comtor'
c
c     Arguments
c
      integer nsets
      real absfacs(2,nsets+1)
      real depdata(maxpts+1,4,nsets+1)
      character*(*) cmd(nsets)
      character*(*) fn(nsets)
c     
c     LOAD_ERODEP_DATA:
c
c     Loads the erosion/deposition data from multiple cases and 
c     assigns the data to the depdata array.
c
c     depdata(wall_element,info_type,case)
c
c     info_type = 1 = net deposition = deposition - erosion
c                 2 = erosion
c                 3 = deposition (ion+neutral)
c                 4 = deposition (ion)
c
c     
c     Local variables
c     
      integer in, is, it
      real normf 
c     
c     Copy needed values for current case.
c     Load current raw data
c     
c     - absfac and absfac_neut
c     - wallse, wallsi, wallsn
c     
      absfacs(1,1) = absfac_neut
      absfacs(2,1) = absfac
c
c     absfac_neut was not saved to the raw file until this 
c     plot was created.   
c
      if (absfacs(1,1).eq.0.0) then 
         absfacs(1,1) = absfacs(2,1)
      endif
c     
      do in = 1,maxpts+1
c     
c        Save data for each of the cases 
c     
         depdata(in,1,1) = wallsi(in)+wallsn(in)-wallse(in)
         depdata(in,2,1) = wallse(in)
         depdata(in,3,1) = wallsi(in)+wallsn(in)
         depdata(in,4,1) = wallsi(in)
      end do  
c
c     Load additional cases and store data to depdata
c     
      do is = 1,nsets
c 
         call loaddata(fn(is)(1:len_trim(fn(is))),
     >                 cmd(is)(1:len_trim(cmd(is))),-1,-1)
c
c     
c        Copy depdata and absfacs 
c     
         absfacs(1,is+1) = absfac_neut
         absfacs(2,is+1) = absfac
c
         if (absfacs(1,is+1).eq.0.0) then 
            absfacs(1,is+1) = absfacs(2,is+1)
         endif
c
         do in = 1,maxpts+1
c     
c           Save data for each of the cases 
c     
            depdata(in,1,is+1) = wallsi(in)+wallsn(in)-wallse(in)
            depdata(in,2,is+1) = wallse(in)
            depdata(in,3,is+1) = wallsi(in)+wallsn(in)
            depdata(in,4,is+1) = wallsn(in)
         end do  
c
c        If the total erosion is tabulated as zero then this
c        was probably from a neutral puff (prior to fix) or
c        an ion injection case. For these situations - assign 
c        the total deposition to the total erosion value. 
c
         if (depdata(maxpts+1,2,is+1).eq.0.0) then 
            depdata(maxpts+1,2,is+1) = depdata(maxpts+1,3,is+1)
         endif
c     
      end do
c     
c     Reload base case data
c     - option strings are not relevant - just the -2 option in the last position
c     
c     
      call loaddata(fn(1)(1:len_trim(fn(is))),
     >              cmd(1)(1:len_trim(cmd(is))),-1,-2)
c     
c     Normalize all the sets of erosion/deposition data to a rate of 1 particle/s/m-tor
c     to make comparisons between the cases easier. 
c
      do is = 1,nsets+1
c 
         normf = depdata(maxpts+1,2,is) 
c
         if (normf.ne.0) then 
c
            do in = 1,maxpts+1
c
               do it = 1,4
c
                  depdata(in,it,is) = depdata(in,it,is) /normf
c
               end do 
c
            end do
c
         endif  
c
      end do
c
      return 
      end
c
c
c
      subroutine calc_erodep_scalef(absfacs,depdata,scalef,ierr) 
      implicit none
      include 'params' 
      include 'comtor'
      include 'cgeom'
      include 'printopt'
c     
      real scalef
      real absfacs(2,2)      
      real depdata(maxpts+1,4,2)      
      integer ierr
c     
c     Total source data is in the maxpts+1 elements 
c     - the total erosion is the equivalent of the absfac_neut scaling
c       factor - each particle count has a weight of absfac/n_eroded 
c     
c     
c     This scaling factor is difficult to calculate and is based on an 
c     analysis of the inner target main SOL only. 
c     
c     1) Check to see if there is any net erosion on the inner target main SOL
c        - if there are none assign the scalef value or assign equal if scalef=0
c     
c     2) Look for inner target elements which have deposition from the undefined
c        scale factor case and net erosion. Find the multiplier which will 
c        give for the worst net erosion case (where there is simultaneous 
c        deposition) a result of zero net erosion.  
c     
c     3) If there are no elements with secondary deposition which also 
c        have net erosion then set the scale factor as in 1.  
c        - if there are none assign the scalef value or assign equal if scalef=0
c     
c     Local variables
c     
      integer iw, is
      integer inner_start,inner_end,istep
      real max_erodep_ratio
c     
c     Initialization
c     
      ierr = 0
c     
      if (inner.eq.'INNER') then 
      
         inner_start = wallindex(idds(irsep,1))
         inner_end   = wallindex(idds(irwall-1,1))
         istep = -1
      
      else
      
         inner_start =  wallindex(idds(irsep,2))
         inner_end   =  wallindex(idds(irwall-1,2))
         istep = 1
      
      endif
c     
      write(6,*) 'DEBUG: PLOT 831 TARGET INDICES:',
     >            inner_start,inner_end,istep
c      
c     Are there any net erosion elements on the inner target from 
c     the base case? Do they also have deposition from the second case?
c     What is the required ratio to get zero net erosion?
c     
      max_erodep_ratio = 0.0 
c     
      write(6,'(a)') 'Calculating MAX ERO/DEP RATIO'//
     >               ' for zero net erosion:'
c
      do iw = inner_start,inner_end,istep 
c     
c        If base case element is net erosion ...
c     
         if (depdata(iw,1,1).lt.0.0) then 
c     
c           Is the second case net deposition?
c     
            if (depdata(iw,1,2).gt.0) then 
c     
c              What is the maximum ratio for elements satisfying both
c              these conditions?
c     
               max_erodep_ratio = max(max_erodep_ratio,
     >                          abs(depdata(iw,1,1)/depdata(iw,1,2)))
c
               write(6,'(a,i5,3(1x,g12.4))') 'Max Ratio:',iw,
     >             max_erodep_ratio,
     >             depdata(iw,1,1),depdata(iw,1,2)
c     
            endif
c     
         endif
c     
      end do         
c     
c     If there is no scaling ratio set  
c     
      if (max_erodep_ratio.eq.0.0) then 
c     
c        Set to specified scale factor
c     
         if (scalef.gt.0.0) then  
c     
            absfacs(1,2) = scalef
            absfacs(2,2) = scalef   
c     
            ierr = 1 
c     
c        Set equal in to each other in second case
c     
         else   
c     
            absfacs(1,2) = absfacs(1,1)
            absfacs(2,2) = absfacs(2,1)
c     
            ierr = 2  
c     
c     
         endif 
c     
c     Set scale factor based on maximum ratio 
c     
      else
c     
         absfacs(1,2) = max_erodep_ratio * absfacs(1,1)
         absfacs(2,2) = max_erodep_ratio * absfacs(2,1)
c     
      endif
c
      write(6,'(a,3(1x,g12.5))') 'PLOT 831:'//
     >     ' Scale calculated: ABS FACTORS:',
     >       absfacs(1,1),
     >       absfacs(1,2),max_erodep_ratio
c     
      return 
      end
c
c
c
      subroutine print_deposition(gridopt)
      implicit none
      include 'params'
      include 'cgeom'
      include 'dynam3'
      include 'comtor'
c      include 'walls_com'
      include 'printopt'
c
      integer ounit,gridopt
c
c     PRINT_DEPOSITION: 
c
c     This routine prints out the deposition data in wallsi and wallsn 
c     across the targets including R,Z coordinate data, cell length, distance
c     along the target from the separatrix and other geometric data. 
c
c     This information is printed to the specified unit
c
c     This routine may be later expanded to include other wall elements as 
c     well as erosion information 
c
      integer id,iw,icnt,ocnt
      real den_m2_tor
      character*1024 :: filename
c
c
      call find_free_unit_number(ounit)

      filename = 'impurity_wall_deposition.dat'

      open(ounit,file=trim(filename),form='formatted')

c
c     Make sure the wall length coordinate calculation has been done - these need to be moved to a central location for efficiency and the code 
c     should set a flag indicating it has been run. 
c
      if (gridopt.eq.RIBBON_GRID) then 
         call calc_wall_length_coordinate(2)
      else
         call calc_wall_length_coordinate(1)
      endif
c
c     Calculate number of data elements that will be listed. 
c
      icnt = 0
      ocnt = 0

      do id = 1,nds

         if (wallindex(id).ne.0) then 
 
            if (id.le.ndsin) then 
               icnt = icnt + 1
            else
               ocnt = ocnt + 1
            endif
         endif
      enddo
c
c
      write(ounit,'(a)')
c
c
      write(ounit,'(a,g18.8)') ' ABSFAC_NEUT:',absfac_neut
      write(ounit,'(a,g18.8)') ' ABSFAC_ION:',absfac
      write(ounit,'(a,g18.8)') ' TOT_NEUT_DEP:',wallsn(maxpts+1)
      write(ounit,'(a,g18.8)') ' TOT_ION_DEP:',wallsi(maxpts+1)
c 
c     First target - "INNER" for X-point up 
c
      write(ounit,'(a,3i6)') 'DEPOSITION DATA: NUMBER OF DATA POINTS:',
     >                icnt,ocnt,wallpts
      write(ounit,'(a)')
      write(ounit,'(a)') 'DEPOSITION DATA for '//INNER//' target'
c
      write(ounit,'(a)')
c
      write(ounit,100) 
     >       'ID','IW','R','Z','DDS','SEPDIST','WALLDIST',
     >       'ION-DEP','NEUT-DEP','DEN/M2-TOR','EROSION'
      
c
      do id = 1, ndsin
c
c
         if (dds(id).gt.0.0) then
            den_m2_tor = (wallsi(wallindex(id))+wallsn(wallindex(id)))
     >          /(wallsi(maxpts+1)+wallsn(maxpts+1))/dds(id)*absfac_neut
         else
            den_m2_tor = 0.0
         endif
c
         if (id.gt.0.and.wallindex(id).gt.0) then
c
c           Use negative values of sepdist to indicate locations
c           in the PFZ
c
            if (irds(id).le.irwall) then 

               write(ounit,200)
     >           id, wallindex(id),
     >           rp(id),zp(id),dds(id),sepdist2(id),
     >           wallpt(wallindex(id),32),
     >           wallsi(wallindex(id)),wallsn(wallindex(id)),
     >           den_m2_tor,wallse(wallindex(id))

             else

               write(ounit,200)
     >           id, wallindex(id),
     >           rp(id),zp(id),dds(id),-sepdist2(id),
     >           wallpt(wallindex(id),32),
     >           wallsi(wallindex(id)),wallsn(wallindex(id)),
     >           den_m2_tor,wallse(wallindex(id))

             endif
c
         endif
c
      end do
c
      write(ounit,'(a)')
c 
c     Second target - "OUTER" for X-point up 
c
      write(ounit,'(a)') 'DEPOSITION DATA for '//OUTER//' target'
c
      write(ounit,'(a)')
c
      write(ounit,100) 
     >       'ID','IW','R','Z','DDS','SEPDIST','WALLDIST',
     >       'ION-DEP','NEUT-DEP','DEN/M2-TOR'
c
      do id = nds,ndsin+1,-1
c
         if (dds(id).gt.0.0) then
            den_m2_tor = (wallsi(wallindex(id))+wallsn(wallindex(id)))
     >       /(wallsi(maxpts+1)+wallsn(maxpts+1))/dds(id)*absfac_neut
         else
            den_m2_tor = 0.0
         endif
c
         if (id.gt.0.and.wallindex(id).gt.0) then
c
c           Use negative values of sepdist to indicate locations
c           in the PFZ
c
            if (irds(id).le.irwall) then 

                write(ounit,200)
     >             id,wallindex(id),
     >             rp(id),zp(id),dds(id),sepdist2(id),
     >             wallpt(wallindex(id),32),
     >             wallsi(wallindex(id)),wallsn(wallindex(id)),
     >             den_m2_tor,wallse(wallindex(id))

            else

                write(ounit,200)
     >             id,wallindex(id),
     >             rp(id),zp(id),dds(id),-sepdist2(id),
     >             wallpt(wallindex(id),32),
     >             wallsi(wallindex(id)),wallsn(wallindex(id)),
     >             den_m2_tor,wallse(wallindex(id))

            endif

         endif
c
      end do
c
      write(ounit,'(a)')
c
      write(ounit,'(a,1x,g18.8)') 'TOTAL_ION_DEPOSITION:', 
     >        wallsi(maxpts+1)
      write(ounit,'(a,1x,g18.8)') 'TOTAL_NEUTRAL_DEPOSITION:', 
     >        wallsn(maxpts+1)
c
      write(ounit,'(a)')
c

c 
c     Entire Wall Deposition 
c
      write(ounit,'(a)') 'DEPOSITION DATA for Entire Wall'
c
      write(ounit,'(a)')
c
      write(ounit,100) 
     >       'IW','ID','R','Z','LEN','       ','WALLDIST',
     >       'ION-DEP','NEUT-DEP','DEN/M2-TOR'

c
      do iw = 1, wallpts
c
         if (wallpt(iw,7).gt.0.0) then
            den_m2_tor = (wallsi(iw)+wallsn(iw))
     >       /(wallsi(maxpts+1)+wallsn(maxpts+1))/wallpt(iw,7)
     >       * absfac_neut
         else
            den_m2_tor = 0.0
         endif
c
         write(ounit,300)
     >       iw,int(wallpt(iw,18)),
     >       wallpt(iw,1),wallpt(iw,2),wallpt(iw,7),wallpt(iw,32),
     >       wallsi(iw),wallsn(iw),den_m2_tor,wallse(iw)
c
      end do
c
c     Quick code to calculate the distance from the separatrix 
c     at the injection position for each ring in the SOL. 
c
      call print_specific_sepdist(45)


      close(ounit)
c
c     Printing formats
c
      ! Labels
      !       'ID','IW','R','Z','DDS','SEPDIST','WALLDIST',
      !       'ION-DEP','NEUT-DEP','DEN/M2-TOR','EROSION'
 100  format(5x,a,5x,a,10x,a,10x,a,8x,a,4x,a,3x,a,8x,a,11x,a,
     >       10x,a,10x,a)
      ! Target Data
 200  format(2(1x,i6),5(1x,f10.6),4(1x,g18.8))
      ! Wall data
 300  format(2(1x,i6),3(1x,f10.6),11x,(1x,f10.6),4(1x,g18.8))
c
      return
      end
c
c
c
      subroutine print_specific_sepdist(ik)
      implicit none
      integer ik
      include 'params'
      include 'cgeom'
c
c     Ik is on separatrix ring
c
      real rsep, zsep, rsep1,zsep1,sdist
      real region_sign
      integer in,ikn,ir
c
c     Calculate rsep and zsep
c
c     A) from grid polygons
c        - separatrix is between corners 1 and 4
c       
c
      in = korpg(ik,irsep)
c
      rsep=(rvertp(1,in) + rvertp(4,in))/2.0
      zsep=(zvertp(1,in) + zvertp(4,in))/2.0
c
      write(6,'(a,i6)') 'Separatrix Distance (cm) for Knot IK = ',ik
c
      do ir = 1,irwall-1
        if (ir.lt.irsep) then 
           ikn=ikins(ik,irsep) 
           region_sign = -1.0
        else
           ikn=ik
           region_sign = 1.0
        endif

        sdist=sqrt( (rs(ikn,ir)-rsep)**2 +
     >              (zs(ikn,ir)-zsep)**2)
     >        * region_sign * 100.0
c

        write(6,'(2(1x,i5),f9.4)') ir,ikn,sdist
c
      enddo 

      return 
      end
c
c
c
      subroutine plot_deposition(iopt,gridopt)
      implicit none
      integer iopt,gridopt

      include 'params'
      include 'walls_com'
      include 'dynam3'
      include 'cgeom'
      include 'outcom'
      include 'printopt'
     
      character*100 graph_819
      real mindist,maxdist,shift_dist
      real scale_min,scale_max,scale_factor
      integer ndatasets,in
      integer iexpt,ierr
      real wall_coord(maxpts)
      real wall_wids(maxpts)
      real wall_dep(maxpts,3)
      real dist
c
      real min_wall_len
      integer start_wall_len,cnt
c
c     Load up case name from environment to include on plot
c
c      CHARACTER casename*128
c



c
c     Read an additional line of plot details required for plot 819.
c
      call rdg_819(graph_819,mindist,maxdist,shift_dist,
     >             scale_min,scale_max,scale_factor,iexpt,ierr)
c
      if (ierr.ne.0) then 
         write(6,*) 'ERROR Reading details for plot 819'
         write(0,*) 'ERROR Reading details for plot 819'
         return
      endif

c
c     Make sure the wall length coordinate is loaded. 
c
      if (gridopt.eq.RIBBON_GRID) then 
         call calc_wall_length_coordinate(2)
      else
         call calc_wall_length_coordinate(1)
      endif

c
c     Load 3 datasets and plot ones selected by iopt
c

      ndatasets = 3
c
c     Load data to suitable arrays for call to DRAW
c
      min_wall_len = HI
c
c     Find coordinate starting point
c
      do in = 1,wallpts
         if (wallpt(in,32).lt.min_wall_len) then 
            start_wall_len = in
            min_wall_len= wallpt(in,32)
         endif
      enddo


      cnt = 1
c
c     Data needs to be ordered properly for plotting - reorganize it 
c     starting at start_wall_len and going backwards/counter clockwise
c

      do

         in = start_wall_len - cnt +1

         if (in.lt.1) in = in + wallpts

         wall_coord(cnt) = wallpt(in,32) + shift_dist  ! Coordinate of wall segment
         wall_wids(cnt)  = wallpt(in,7)                ! length of wall segment 
         ! wall deposition density in terms of 1 particle entering /s
         ! den/m2/m-tor
         if (wallpt(in,7).gt.0.0) then
            wall_dep(cnt,1) = (wallsi(in)+wallsn(in))
     >             /(wallsi(maxpts+1)+wallsn(maxpts+1))/wallpt(in,7)
     >             * scale_factor
            wall_dep(cnt,2) = (wallsi(in))
     >             /(wallsi(maxpts+1)+wallsn(maxpts+1))/wallpt(in,7)
     >             * scale_factor
            wall_dep(cnt,3) = (wallsn(in))
     >             /(wallsi(maxpts+1)+wallsn(maxpts+1))/wallpt(in,7)
     >             * scale_factor
         else
            wall_dep(in,1:3) = 0.0
         endif
c
         write(6,'(a,3i6,16(1x,g12.5))') 'DEP 819:',cnt,in,wallpts,
     >             wallpt(in,1),wallpt(in,2),
     >             wall_coord(cnt),
     >             wall_dep(cnt,1),wall_dep(cnt,2),
     >             wall_dep(cnt,3),
     >             wall_wids(cnt),
     >             wallsi(in),wallsn(in),
     >             scale_factor,shift_dist
c
c
c        Increment the wall element count
c
         cnt = cnt+1
c
c        Loop exit condition
c
         if (cnt.ge.wallpts+1) exit
c
      enddo 

      if (mindist.lt.wall_coord(1)) mindist = wall_coord(1)
      if (maxdist.gt.wall_coord(wallpts)) maxdist= wall_coord(wallpts)



c
c     Set up scaling
c
      if (scale_min.eq.0.0)  scale_min = -HI
      if (scale_max.eq.0.0)  scale_max =  HI

c
c     Do not overwrite Title or job
c
      anly = ' '
      ref  = ' '
      plane= ' '
      nview= ' '
c
      table= 'DEPOSITION'
      xlab = 'Distance Along Wall (m)' 
      ylab = 'Deposition'
      Elabs(1) = 'I+N  Ion+Neut Deposition'
      Elabs(2) = 'I    Ion      Deposition'
      Elabs(3) = 'N    Neut     Deposition'
c
c     Load strike point locations into 2 data fields
c
c
c     INNER for JET - PLANE
c      
      in = wallindex(idds(irsep,1))
      dist = wallpt(in,32) - wallpt(in,7)/2.0
      write(plane,'(a,f8.2)') INNER//' Strike Point =', dist+shift_dist
c
c     OUTER for JET - NVIEW
c
      in = wallindex(idds(irsep,2))
      dist = wallpt(in,32) + wallpt(in,7)/2.0
      write(nview,'(a,f8.2)') OUTER//' Strike Point =', dist+shift_dist
c
c     Describe the axis used - ANLY and REF
c
      write(ref(1:36),'(a)') 'Axis is CCW from inner midplane'
      write(anly,'(a,f7.3)') 'Axis shift =',shift_dist

c...  Load case name from environment variable CASENAME: - REF
c      CALL IoName(casename,ierr)
c
c      write(ref,'(a)') casename(1:len_trim(casename))
c
c
      ignors = 0
c
      if (iopt.ge.4) then 
         ignors = 1
      else
         ignors(iopt) =1 
      endif
c
c     Plot the deposition - include experimental data set if one is specified
c
      CALL DRAW (wall_coord,wall_wids,wall_dep,MAXPTS,
     >        wallpts,ANLY,
     >        ndatasets,99,mindist,maxdist,scale_min,scale_max,
     >        IGNORS,ITEC,AVS,NAVS,
     >        JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,
     >        PLANE,TABLE,1,2,1.0,iexpt)

      ignors = 1


      return
      end



      subroutine setup_pnames(pnames1,pnames2,wlmax,opt)
      implicit none
      include 'params'
      include 'walls_com'
      include 'printopt'
      character*15 pnames1(maxpts),pnames2(maxpts)
      integer wlmax
      integer opt

c
c        Set up Pnames and cnames labels
c
         if ((wltrap1.gt.0.and.wltrap1.le.maxpts).and.
     >      (wltrap2.gt.0.and.wltrap2.lt.maxpts)) then 

           pnames1(wltrap1) = '|'
           pnames1(wltrap2) = '|'
           pnames1(int((wltrap1+wltrap2)/2.0)) = 'PP'
           pnames2(int((wltrap1+wltrap2)/2.0)) = 'Wall'
           pnames1((wltrap2+wallpts)/2) = OUTER
           pnames2((wltrap2+wallpts)/2) = 'Target'

         endif

         if ((wlwall1.gt.0.and.wlwall1.lt.maxpts).and.
     >       (wlwall2.gt.0.and.wlwall2.lt.maxpts)) then 

            pnames1(wlwall1) = '|'
            pnames1(wlwall2) = '|'
            pnames1(wlwall1+ INT(0.1*(wlwall2-wlwall1))) = Outer
            pnames1(wlwall1+ INT(0.9*(wlwall2-wlwall1))) = Inner
c           pnames1((wlwall1+wlwall2)/2) = 'Top'
c           pnames2((wlwall1+wlwall2)/2) = 'Z (m)'
c
c        0.02 - Outside
c
            pnames1(wlwall1+ INT(0.02*(wlwall2-wlwall1))) = '^'
         write(pnames2(wlwall1+ INT(0.02*(wlwall2-wlwall1))),'(a,f4.1,
     >     a)') 'Z=',wallpt(wlwall1+ INT(0.02*(wlwall2-wlwall1)),2),'m'
c
c        0.2
c
            pnames1(wlwall1+ INT(0.2*(wlwall2-wlwall1))) = '^'
         write(pnames2(wlwall1+ INT(0.2*(wlwall2-wlwall1))),'(f4.1)')
     >              wallpt(wlwall1+ INT(0.2*(wlwall2-wlwall1)),2)
c
c        0.4
c
            pnames1(wlwall1+ INT(0.4*(wlwall2-wlwall1))) = '^'
         write(pnames2(wlwall1+ INT(0.4*(wlwall2-wlwall1))),'(f4.1)')
     >              wallpt(wlwall1+ INT(0.4*(wlwall2-wlwall1)),2)
c
c        0.6
c
            pnames1(wlwall1+ INT(0.6*(wlwall2-wlwall1))) = '^'
         write(pnames2(wlwall1+ INT(0.6*(wlwall2-wlwall1))),'(f4.1)')
     >              wallpt(wlwall1+ INT(0.6*(wlwall2-wlwall1)),2)
c
c        0.8
c
            pnames1(wlwall1+ INT(0.8*(wlwall2-wlwall1))) = '^'
         write(pnames2(wlwall1+ INT(0.8*(wlwall2-wlwall1))),'(f4.1)')
     >              wallpt(wlwall1+ INT(0.8*(wlwall2-wlwall1)),2)
c
c        0.98 - Inside
c
            pnames1(wlwall1+ INT(0.98*(wlwall2-wlwall1))) = '^'
         write(pnames2(wlwall1+ INT(0.98*(wlwall2-wlwall1))),'(f4.1)')
     >              wallpt(wlwall1+ INT(0.98*(wlwall2-wlwall1)),2)
c

            pnames1((wlwall2+wltrap1)/2) = INNER
            pnames2((wlwall2+wltrap1)/2) = 'Target'


         endif

c
c        Mark Z-distances for some wall sections
c
c        Max of wall
c
         pnames1(wlmax) = 'Top'
c
c         pnames2(wlmax) = 'Main Wall'
c
c         write(pnames2(wlmax),'(f5.1)') wallpt(wlmax,2)
c         pnames2((wlwall1+wlwall2)/2) = 'Main Wall'
c
       return
       end


