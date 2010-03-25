c     -*-Fortran-*-
c
c ======================================================================
c
c subroutine: XLOGSCALE
c
      SUBROUTINE XLogScale
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slout'
      include 'comgra'
c
c     jdemod - use the include for comgra so if changes are made they 
c              are reflected everywhere
c 
c      COMMON /COMGRA/ CXMIN,CXMAX,CYMIN,CYMAX,IPLOTS,ICOL,NPLOTS,ISPOT
c      REAL            CXMIN,CXMAX,CYMIN,CYMAX
c      INTEGER         IPLOTS,ICOL,NPLOTS,ISPOT
c
      REAL decade,decval,t,val,x1

      cxmin = LOG10(cxmin)
      cxmax = LOG10(cxmax)

      CALL PSPACE (0.0, 1.35, 0.0, 1.0)
      CALL MAP    (0.0, 1.35, 0.0, 1.0)

      DO decade = REAL(INT(cxmin)-1), REAL(INT(cxmax)+2)
        decval = 10.0**decade

        DO val = decval, decval*10.0, decval
          t = (LOG10(val) - cxmin) / (cxmax - cxmin)

          IF (t.LT.0.0.OR.t.GT.1.0) CYCLE

          x1 = map1x + t * (map2x - map1x)
c          x1 = 0.1 + t * (0.9 - 0.1)

          IF (val.EQ.decval) THEN
            CALL POSITN(x1,map1y)
            CALL JOIN  (x1,map1y+0.01)
            CALL POSITN(x1,map2y)
            CALL JOIN  (x1,map2y-0.01)

            CALL PCSCEN(x1,map1y-0.025,'10')
            CALL POSITN(x1+0.006,map1y-0.025+0.009)
            CALL TYPENI(INT(decade))
          ELSE
            CALL POSITN(x1,map1y)
            CALL JOIN  (x1,map1y+0.005)
            CALL POSITN(x1,map2y)
            CALL JOIN  (x1,map2y-0.005)
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END
c
c ======================================================================
c
c subroutine: YLOGSCALE
c
      SUBROUTINE YLogScale
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slout'
      include 'comgra'
c
c     jdemod - use the include for comgra so if changes are made they 
c              are reflected everywhere
c 
c      COMMON /COMGRA/ CXMIN,CXMAX,CYMIN,CYMAX,IPLOTS,ICOL,NPLOTS,ISPOT
c      REAL            CXMIN,CXMAX,CYMIN,CYMAX
c      INTEGER         IPLOTS,ICOL,NPLOTS,ISPOT
c
      REAL decade,decval,t,val,y1

      cymin = LOG10(cymin+LO)
      cymax = LOG10(cymax+LO)

      CALL PSPACE (0.0, 1.35, 0.0, 1.0)
      CALL MAP    (0.0, 1.35, 0.0, 1.0)

      DO decade =  REAL(INT(cymin)-1), REAL(INT(cymax)+2)

        decval = 10.0**decade + LO

        DO val = decval, decval*10.0, decval

          t = (LOG10(val) - cymin) / (cymax - cymin)

          IF (t.LT.0.0.OR.t.GT.1.0) CYCLE

          y1 = map1y + t * (map2y - map1y)
c          y1 = 0.11 + t * (0.89 - 0.11)

          IF (val.EQ.decval) THEN
            CALL POSITN(map1x      ,y1)
            CALL JOIN  (map1x+0.010,y1)
            CALL POSITN(map2x      ,y1)
            CALL JOIN  (map2x-0.010,y1)

            CALL PCSCEN(map1x-0.037,y1,'10')
            CALL POSITN(map1x-0.037+0.008,y1+0.009)
            CALL TYPENI(INT(decade))
          ELSE
            CALL POSITN(map1x      ,y1)
            CALL JOIN  (map1x+0.005,y1)

            CALL POSITN(map2x      ,y1)
            CALL JOIN  (map2x-0.005,y1)
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END




c slmod begin test
      subroutine SLDRAWM (MOUTS,MWIDS,MVALS,MAXNPS,pnks,
     >                  nplts,ngs,pltlabs,elabs,xlab,ylab,ref,title,
     >                  sctype,numplots,pltmins,pltmaxs,pltfact)
c     >                  sctype,numplots,pltmins,pltmaxs,pltfact,
c     .                  grm_opt,ngs2,grm_shade,elabs2)

      implicit none

      INCLUDE 'params'
      INCLUDE 'slout'

      integer MAXNPS,numplots
      real mouts(maxnps,maxplts)
      real mvals(maxnps,maxplts,maxngs)
      real mwids(maxnps,maxplts)
      real pltmins(maxplts),pltmaxs(maxplts)
      real pltfact
      integer pnks(maxplts)
      integer nplts,ngs,sctype
      character*(*) pltlabs(maxplts)
      character*(*) elabs(maxngs)
c
c      subroutine DRAWM (MOUTS,MWIDS,MVALS,MAXNMS,maxplts,maxngs,pnks,
c     >                  nplts,ngs,pltlabs,elabs,xlab,ylab,ref,title,
c     >                  sctype,numplots,pltmins,pltmaxs,pltfact)
c      implicit none
c      integer maxnms,maxplts,maxngs,numplots
c      real mouts(maxnms,maxplts)
c      real mvals(maxnms,maxplts,maxngs)
c      real mwids(maxnms,maxplts)
c      real pltmins(maxplts),pltmaxs(maxplts)
c      real pltfact
c      integer pnks(maxplts)
c      integer nplts,ngs,sctype
c      character*(*) pltlabs(maxplts)
c      character*(*) elabs(maxngs)
c slmod end
      character*(*) xlab,ylab,ref,title
c
      include 'grminfo'
c
c slmod begin
c      INTEGER         grm_opt,ngs2(MAXNGS),slopt,plottype(MAXNGS)
c      REAL            grm_shade(2,MAXNGS)
c      CHARACTER*36    elabs2   (8,MAXNGS)
c
c      real hi
c      parameter(hi=1.0e37)
      CHARACTER cdum1*256,dummy*256
c slmod end
C
C  *********************************************************************
C  *                                                                   *
C  * DRAWM:   DRAW A GRAPH OF VALUES IN MVALS AGAINST VALUES IN MOUTS  *
C  *          OVER THE RANGE OF VALUES IN MOUTS. PUTTING UP TO 8 PLOTS *
C  *          ON ONE PAGE                                              *
C  *                                                                   *
C  * ARGUMENTS :-                                                      *
C  *                                                                   *
C  * MOUTS  - VALUES TO GO ON X AXIS OF GRAPHS (FOR EACH DATA SET)     *
C  * MWIDS  - BIN WIDTHS  (FOR SMOOTHING PURPOSES)                     *
C  * MVALS  - VALUES TO GO ON Y AXIS OF GRAPHS (FOR EACH DATA SET)     *
C  * MAXNMS - FIRST DIMENSION OF M ARRAYS (IK INDEX)                   *
C  * MAXPLTS- SECOND DIMENSION OF M ARRAYS - MAX NUM OF PLOTS          *
C  * MAXNGS - MAX NUMBER OF DATA SETS FOR EACH PLOT                *
C  * PNKS   - NUMBER OF IK VALUES IN EACH SET OF DATA              *
C  * NPLTS  - ACTUAL NUMBER OF PLOTS TO BE MADE                        *
C  * NGS    - ACTUAL NUMBER OF DATA SETS ON EACH PLOT                  *
C  * PLTLABS- LABELS SPECIFIC TO EACH PLOT                             *
C  * ELABS  - LABELS SPECIFIC TO EACH DATA SET                 *
C  * XLAB   - X-AXIS LABEL                         *
C  * YLAB   - Y-AXIS LABEL                         *
C  * REF    - REFERENCE LABEL FOR ALL PLOTS                *
C  * TITLE  - TITLE OF CASE                            *
c  * PLTMINS- EXTERNAL PLOT MINIMUMS                                   *
c  * PLTMAXS- EXTERNAL PLOT MAXIMUMS                                   *
C  * SCTYPE - TYPE OF PLOT SCALING TO BE EMPLOYED              *
C  *          1 - ALL PLOTS ON PAGE SHARE THE SAME SCALING - MIN AND   *
C  *              MAX CALCULATED FOR ALL DATA                          *
C  *          2 - MINIMUM AND MAXIMUM CALCULATED FOR EACH PLOT         *
C  *              SEPARATELY ON THE PAGE.                              *
C  *          3 - MAX CALCULATED FOR ALL PLOTS ON PAGE - MIN IS ZERO   *
C  *              UNLESS PLOT CROSSES X-AXIS.                  *
C  *          4 - MAXIMUM CALCULATED FOR EACH PLOT SEPARATELY - MIN IS *
C  *              SET TO 0.0 UNLESS PLOT CROSSES X-AXIS            *
C  *          5 - Logarithmic Scaling                                  *
C  *              MAX CALCULATED FOR ALL PLOTS ON PAGE - MIN IS ZERO   *
C  *              UNLESS PLOT CROSSES X-AXIS.                  *
C  *          6 - Logarithmic Scaling                                  *
C  *              MAXIMUM CALCULATED FOR EACH PLOT SEPARATELY - MIN IS *
C  *              SET TO 0.0 UNLESS PLOT CROSSES X-AXIS            *
c  *          7 - Externally provided plot minimums and maximums.      *
C  *                                           *
C  *                                           *
C  *          David Elder, NOV 21, 1995                    *
C  *                                                                   *
C  *********************************************************************
C
c
c     Local Variables
c
      integer in,ip,ik,ir,ipmax,id,posin
      real pltmax,pltmin,axmax,axmin
c slmod begin
      INTEGER i1
c slmod end
c
c     Set up number of plots / page
c
      pageplots = numplots
c
c     Add overall Page labels and references
c
c slmod begin
      CALL THICK2(4)

      call SLgrmtitle(title,ref,nplts)
c
c      call SLgrmtitle(title,ref)
c slmod end
c
c     Loop through all the plots on the page
c
c
c     Find the Max and Min for the plots - the axis scale
c     Max and Min are always calculated separately for each plot.
c     The axis max and min assume that the array of axis coordinates
c     is in ascending order.
c
c     Global scaling
c
      write (6,*) 'sctype:',sctype,':',xlab,':',ylab,':'
c
c     Write out Mvals:
c
c      write(6,*) 'nplts:',nplts,ngs,numplots
c      do ip = 1,nplts
c        do ik = 1,pnks(ip)
c          write(6,*) 'info2:',ik,ip,pnks(ip),mouts(ik,ip)
c          write(6,'(i4,'':'',3(2x,e13.4))') ngs,
c     >                     (mvals(ik,ip,in),in=1,ngs)
c        end do
c      end do
c
c
c     Convert all values to logarithmic for sctype 5 and 6.
c
      if (sctype.eq.5.or.sctype.eq.6) then
        do ip = 1,nplts
          do ik = 1,pnks(ip)
            do in = 1,ngs
              if (mvals(ik,ip,in).ne.0.0) then
                mvals(ik,ip,in) = log10(mvals(ik,ip,in))
              else
                mvals(ik,ip,in) = 0.0
              endif
            end do
          end do
        end do
      endif
c
      pltmax = -hi
      pltmin =  hi
      if (sctype.eq.1.or.sctype.eq.3.or.sctype.eq.5) then
        do ip = 1,nplts
          write (6,*) 'ip:',ip,pnks(ip),ngs
          do ik = 1,pnks(ip)
            do in = 1,ngs
              pltmax = max(pltmax,mvals(ik,ip,in))
              pltmin = min(pltmin,mvals(ik,ip,in))
c              write (6,*) 'plt:',ik,ip,in,pltmax,pltmin,
c     >                     mvals(ik,ip,in)
            end do
          end do
        end do
        if (sctype.eq.3) then
          if (pltmax.ge.0.0.and.pltmin.ge.0.0) then
            pltmin = 0.0
          elseif (pltmax.le.0.0.and.pltmin.le.0.0) then
            pltmax = 0.0
          endif
        endif
c
c       Rescale plot max and min - just so the plots don't hit the edges
c
        if (pltmax.ge.0.0) then
           pltmax = pltmax * 1.05
        else
           pltmax = pltmax * 0.95
        endif
c
        if (pltmin.ge.0.0) then
           pltmin = pltmin * 0.95
        else
           pltmin = pltmin * 1.05
        endif
c
      endif
c
      do id = 1,nplts,pageplots
c
        ipmax = min(nplts,id+pageplots-1)
c
        write (6,*)'id:',id,nplts,pageplots,ipmax,grm_opt,ngs
c
        do ip = id,ipmax
c slmod begin
          IF (pltlabs(ip)(1:3).EQ.'-99') CYCLE

          WRITE(6,*) 'pnks:',ip,id,ipmax,pnks(ip)

          IF (pnks(ip).EQ.0) CYCLE

          IF (grm_opt.EQ.1) THEN
            ngs = ngs2(ip)

            DO i1 = 1, ngs
              ylab         = ylab2    (ip)
              elabs   (i1) = elabs2   (ip,i1)
              plottype(i1) = plottype2(ip,i1)
            ENDDO
          ENDIF
c slmod end
c
c         Set externally imposed scales for scale option 7
c
          if (sctype.eq.7) then
             pltmin = pltmins(ip)
             pltmax = pltmaxs(ip)
          endif
c
c         SET UP BOX FOR EACH PLOT - BASED ON POSITION/PLOT NUMBER
c
c slmod begin
c          call SLgrmbox(ip)
c
c         Add labels for plot (individual plot titles and annotation.
c
c          call SLgrmlabels(ip,pltlabs(ip),elabs,ngs)
c slmod end
c
c         Label and plot the axes for the given plot position
c
          axmin = mouts(1,ip)
          axmax = mouts(pnks(ip),ip)
c slmod begin
          write (6,*) 'range 1:',axmin,axmax,pltfact,
     .      abs(pltfact),abs(pltfact) * (axmax-axmin)
c slmod end
c
c         Allow for closeup plots at either end of the ring by passing
c         another parameter - if 0.0 it plots the whole scale.
c         if < 0.0 it plots from axmin to abs(pltfact)*(axmax-axmin)
c         and if greater than 0.0 it plots from
c         (1.0-pltfact)*(axmax-axmin) to axmax.
c
          if (pltfact.lt.0.0) then
             if (abs(pltfact).gt.1.0) pltfact = 0.1
c slmod begin
             axmax = abs(pltfact) * (axmax-axmin) + axmin
c
c             axmax = abs(pltfact) * (axmax-axmin)
c slmod end
          elseif (pltfact.gt.0.0) then
             if (pltfact.gt.1.0) pltfact = 0.1
c slmod begin
             axmin = (1.0-pltfact) * (axmax-axmin) + axmin
c
c             axmin = (1.0-pltfact) * (axmax-axmin)
c slmod end
          endif
c
c slmod begin
          write (6,*) 'range 2:',axmin,axmax,pltfact
c slmod end

          if (sctype.eq.2.or.sctype.eq.4.or.sctype.eq.6) then
            pltmax = -hi
            pltmin =  hi
            do ik = 1, pnks(ip)
              do in = 1,ngs
c slmod begin
                IF (mouts(ik,ip).GE.axmin.AND.
     .              mouts(ik,ip).LE.axmax) THEN
                  pltmax = max(pltmax,mvals(ik,ip,in))
                  pltmin = min(pltmin,mvals(ik,ip,in))
                ENDIF
c
c                pltmax = max(pltmax,mvals(ik,ip,in))
c                pltmin = min(pltmin,mvals(ik,ip,in))
c slmod end
              end do
            end do
            if (sctype.eq.4) then
              if (pltmax.ge.0.0.and.pltmin.ge.0.0) then
                pltmin = 0.0
              elseif (pltmax.le.0.0.and.pltmin.le.0.0) then
                pltmax = 0.0
              endif
            endif
c
c           Rescale plot max and min - just so the plots
c                                      don't hit the edges
c
c slmod begin
            IF (grm_opt.EQ.1) THEN
              if (pltmax.ge.0.0) then
                 pltmax = pltmax + 0.15 * ABS(pltmax - pltmin)
              else
                 pltmax = pltmax + 0.15 * ABS(pltmax - pltmin)
              endif

              if (pltmin.ge.0.0) then
                 pltmin = pltmin - 0.05 * ABS(pltmax - pltmin)
              else
                 pltmin = pltmin - 0.05 * ABS(pltmax - pltmin)
              endif
            ELSE
              if (pltmax.ge.0.0) then
                 pltmax = pltmax * 1.05
              else
                 pltmax = pltmax * 0.95
              endif

              if (pltmin.ge.0.0) then
                 pltmin = pltmin * 0.95
              else
                 pltmin = pltmin * 1.05
              endif
            ENDIF
c
c            if (pltmin.ge.0.0) then
c               pltmin = pltmin * 0.95
c            else
c               pltmin = pltmin * 1.05
c            endif
c
c            if (pltmax.ge.0.0) then
c               pltmax = pltmax * 1.05
c            else
c               pltmax = pltmax * 0.95
c            endif
c
c            if (pltmin.ge.0.0) then
c               pltmin = pltmin * 0.95
c            else
c               pltmin = pltmin * 1.05
c            endif
c slmod end
c
          endif
c
          write (6,*) 'range 3:',axmin,axmax,pltmin,pltmax
c slmod begin
          call SLgrmlabels(ip,pltlabs(ip),elabs,ngs)
c slmod end
c
          call SLgrmaxes(ip,pltmin,pltmax,axmin,axmax,
     >                 xlab,ylab,sctype)
c
c         Plot data for this plot
c
c slmod begin
          call SLgrmdata(mvals,mouts,pnks,maxnps,maxplts,
     >               maxngs,ip,ngs,
     >               pltmin,pltmax,axmin,axmax,elabs,sctype,
     .               slopt,slopt2,slopt3,plottype)

c... Symbol overlay:
          IF (slopt4.GT.0) THEN
            plottype(1) = slopt4 + 1

            call SLgrmlabels(ip,pltlabs(ip),elabs,ngs)

            call SLgrmdata(mvals,mouts,pnks,maxnps,maxplts,
     >                 maxngs,ip,ngs,
     >                 pltmin,pltmax,axmin,axmax,elabs,sctype,
     .                 slopt,slopt2,slopt3,plottype)

          ENDIF

          call SLgrmbox(ip,axmin,axmax,pltmin,pltmax,MAXNPS)

c
c          call SLgrmdata(mvals,mouts,pnks,maxnms,maxplts,
c     >               maxngs,ip,ngs,
c     >               pltmin,pltmax,axmin,axmax,elabs,sctype)
c slmod end
c
        end do
c
c slmod end
      READ(5,'(A256)') dummy
      IF   (dummy(8:14).EQ.'Noframe'.OR.dummy(8:14).EQ.'noframe'.OR.
     .      dummy(8:14).EQ.'NOFRAME') THEN
      ELSE
        CALL FRAME
        CALL THICK2(4)
        BACKSPACE 5
      ENDIF
c
c        call frame
c slmod end
c
      end do
c
      return
      end
c
c
c
c slmod begin
      subroutine SLgrmtitle(title,ref,nplts)
      implicit none
      character*(*) title,ref
      INTEGER nplts

      INCLUDE 'params'
      INCLUDE 'slout'
c
c      subroutine grmtitle(title,ref)
c      implicit none
c      character*(*) title,ref

c slmod end
c
c     GRMTITLE: This function places the titles on the page.
c
      integer l,lenstr,i,j
      external lenstr
c
      CALL PSPACE (0.0, 1.0, 0.0, 1.0)
      CALL MAP    (0.0, 1.0, 0.0, 1.0)
      CALL CTRMAG (10)
c
c slmod begin
      IF (nplts.EQ.1) THEN
        IF (slopt2.EQ.1) THEN
c         IPP/08 Krieger - revert to thin lines for the moment
c         CALL THICK (2)
          CALL THICK (1)
          CALL CTRMAG(30)
          CALL PLOTCS(0.03,0.94,TITLE(1:MIN(45,LEN_TRIM(TITLE))))
          L = LENSTR(REF)
          CALL THICK (1)
          CALL CTRMAG(8)
          CALL PCSCEN (1.12,0.005, ref(:l))
        ELSE
          L = LENSTR(TITLE)
          CALL PCSCEN (0.25+0.05,0.99, title(:l))
          L = LENSTR(REF)
          CALL PCSCEN (0.75+0.05,0.99, ref(:l))
        ENDIF
      ELSE
        L = LENSTR(TITLE)
        CALL PCSCEN (0.5+0.05,0.99, title(:l))

        L = LENSTR(REF)
        CALL PCSCEN (0.5+0.05,0.975, ref(:l))
      ENDIF

      CALL CTRMAG (10)
      DO l = 1, 10
        IF (char(l).NE.' ') CALL PLOTST (1.04,0.737+(l-1)*0.02,char(l))
      ENDDO

      DO i = 11, 19
        j = i - 10
        IF (char(i).NE.' ') CALL PLOTST (1.04,0.500-(j-1)*0.02,char(i))
      ENDDO

      DO i = 20, 30
        j = i - 19
        IF (char(i).NE.' ') CALL PLOTST (1.04,0.250-(j-1)*0.02,char(i))
      ENDDO

c
c      L = LENSTR(TITLE)
c      CALL PCSCEN (0.5,0.99, title(:l))
c
c      L = LENSTR(REF)
c      CALL PCSCEN (0.5,0.975, ref(:l))
c slmod end
c
      return
      end
c
c
c
c slmod begin
      subroutine SLgrmbox(boxindex,axmin,axmax,pltmin,pltmax,MAXNPS)
c
c      subroutine grmbox(boxindex)
c slmod end
      implicit none
      integer MAXNPS,boxindex
c slmod begin
      INCLUDE 'params'
      INCLUDE 'slcom'
      INCLUDE 'slout'
      INCLUDE 'colours'

      INTEGER i1,i2
      REAL    rv(4),zv(4)
c slmod end
c
c     GRMBOX: DRAWS THE OUTSIDE BOX FOR THE PLOT INDICATED BY
c             THE BOXINDEX.
c
c
      real xpt,ypt,xwid,ywid,xsep,ysep
c slmod begin
      real pltmin,pltmax,axmin,axmax
c slmod end
c
      call SLgrmfindxy(boxindex,xpt,ypt,xwid,ywid,xsep,ysep)
c slmod begin
      IF (.TRUE..AND.grm_shade(IKLO,boxindex).NE.0.0) THEN
        CALL PSPACE(xpt  ,xpt+xwid,ypt   , ypt+ywid)
        CALL MAP   (axmin,axmax   ,pltmin,pltmax   )

c        CALL FILCOL(ncols)
c        CALL LINCOL(ncols)
        CALL FILCOL(ncols-1)
        CALL LINCOL(ncols-1)

        rv(1) = axmin
        zv(1) = pltmin
        rv(2) = axmin
        zv(2) = pltmax
        rv(3) = GRM_shade(IKLO,boxindex)
        zv(3) = pltmax
        rv(4) = GRM_shade(IKLO,boxindex)
        zv(4) = pltmin

         CALL POSITN (rv(3),zv(3))
         CALL JOIN   (rv(4),zv(4))

c        CALL PTPLOT(rv(1),zv(1),1,4,1)


c        WRITE(0,*) 'MARK: NCOLS= ',ncols
        CALL FILCOL(8)
        CALL LINCOL(8)
        i2 = 1
        DO i1 = 1, MAXNKS-1
c        DO i1 = 1, MAXNPS-1
          IF (grm_cell(i1-1,boxindex).LT.axmax.AND.
     .        grm_cell(i1  ,boxindex).LT.axmax.AND.
     .        .NOT.(grm_cell(i1-1,boxindex).EQ.0.0.AND.
     .              grm_cell(i1  ,boxindex).EQ.0.0)) THEN

            IF (i2.EQ.1) THEN
c              WRITE(0,*) 'MARK: DRAW I1= ',i1,grm_cell(i1-1,boxindex),
c     .                                        grm_cell(i1  ,boxindex)
              rv(1) = grm_cell(i1-1,boxindex)
              zv(1) = pltmin * 0.05
              rv(2) = grm_cell(i1-1,boxindex)
              zv(2) = pltmax * 0.05
              rv(3) = grm_cell(i1  ,boxindex)
              zv(3) = pltmax * 0.05
              rv(4) = grm_cell(i1  ,boxindex)
              zv(4) = pltmin * 0.05

c              CALL PTPLOT(rv(1),zv(1),1,4,1)

              i2 = 0
            ELSE
              i2 = 1
            ENDIF
          ENDIF
        ENDDO

        CALL LINCOL(1)

        grm_shade(IKLO,boxindex) = 0.0

      ENDIF

      IF (.FALSE..AND.grm_shade(IKHI,boxindex).NE.0.0) THEN
        CALL PSPACE(xpt  ,xpt+xwid,ypt   ,ypt+ywid)
        CALL MAP   (axmin,axmax   ,pltmin,pltmax  )

c        CALL FILCOL(ncols)
c        CALL LINCOL(ncols)
        CALL FILCOL(ncols-1)
        CALL LINCOL(ncols-1)

        rv(1) = grm_shade(IKHI,boxindex)
        zv(1) = pltmin
        rv(2) = grm_shade(IKHI,boxindex)
        zv(2) = pltmax
        rv(3) = axmax
        zv(3) = pltmax
        rv(4) = axmax
        zv(4) = pltmin

        CALL PTPLOT(rv(1),zv(1),1,4,1)

        CALL LINCOL(1)

        grm_shade(IKHI,boxindex) = 0.0

      ENDIF

      CALL PSPACE(0.0,1.0,0.0,1.0)
      call MAP   (0.0,1.0,0.0,1.0)
c slmod end
c
c     Set up vector/ND space mapping
c
      call pspace(0.0,1.0,0.0,1.0)
      call map (0.0,1.0,0.0,1.0)
      call full
      CALL LINCOL (1)
      CALL THICK  (1)
c
c     Start drawing box
c
      call positn(xpt,ypt)
c
c     Draw box
c
      call join(xpt+xwid,ypt)
      call join(xpt+xwid,ypt+ywid)
      call join(xpt,ypt+ywid)
      call join(xpt,ypt)
     
c      WRITE(0,*) '-->',xpt,xpt+xwid,ypt,ypt+ywid

c
      return
      end
c
c
c
      subroutine SLgrmlabels(ip,pltlabs,elabs,ngs)
      implicit none
      integer ip,ngs
      character*(*) pltlabs
      character*(*) elabs(ngs)
c
c     GRMLABELS : This routine places the labels on
c                 each individual plot.
c
c
      integer in,l,lenstr
      external lenstr
      real xpt,ypt,xwid,ywid,xsep,ysep,xp,yp
c
      include 'colours'
c slmod begin
      INCLUDE 'params'
      INCLUDE 'slout'

      INTEGER j,in1,in2,in2a
      CHARACTER*1 ch
c slmod end
c
c      INTEGER COLOUR(8)
c
c      integer icol
c
c      DATA COLOUR /1,2,3,4,9,7,6,8/
c
      call SLgrmfindxy(ip,xpt,ypt,xwid,ywid,xsep,ysep)
c
      CALL PSPACE (0.0, 1.35, 0.0, 1.0)
      CALL MAP    (0.0, 1.35, 0.0, 1.0)
      CALL CTRMAG (10)
      call full
      call thick(1)
c
      L = LENSTR(PLTLABS)
c slmod begin
      IF (grm_opt.EQ.1.OR..TRUE.) THEN
        IF (slopt2.EQ.1) THEN
          CALL CTRMAG(20)
          CALL PCSEND(xpt+0.90*xwid,ypt+0.90*ywid,pltlabs(:L))
        ELSE
          CALL PCSEND(xpt+xwid-0.15*xsep,ypt+ywid-0.35*ysep,pltlabs(:L))
c          CALL PCSEND(xpt+xwid-0.10*xsep,ypt+ywid-0.35*ysep,pltlabs(:L))
        ENDIF
      ELSE
        CALL PCSEND (xpt+xwid/2.0,ypt+ywid+0.60*ysep, pltlabs(:L))
      ENDIF
c
c      CALL PCSCEN (xpt+xwid/2.0,ypt+ywid+0.75*ysep, pltlabs(:L))
c slmod end
c
      icol = 0
c
      DO IN = 1,NGS
c
c       Set line colour
c
        icol = icol + 1
        IF (ICOL.GT.NCOLS) ICOL = 1
        CALL LINCOL (COLOUR(ICOL))
c
c       Set line type
c
        IF (in.LE.1) THEN
           CALL FULL
        ELSEIF (in.LE.5) THEN
           CALL BROKEN (3*in,2*in,3*in,2*in)
        ELSE
           CALL BROKEN (2*in,1*in,2*in,1*in)
        ENDIF
c slmod begin
        IF (slopt.GT.0.OR.slopt2.GT.0) THEN
          IF (plottype(icol).GT.1) THEN
c            CALL FULL
            CALL LinCol(colour(plottype(icol)+ncols))
          ENDIF
        ENDIF
c slmod end
c
c       Draw line segment.
c
c slmod begin
        IF (slopt2.EQ.1) THEN
          xp = xpt + 1.02*xwid + 0.25 * xsep
          yp = ypt + ywid - (in - 1) * (ywid / ngs) * 0.4
          call full
c         IPP/08 Krieger - revert to thin lines for the moment
c         CALL THICK(2)
          CALL THICK(1)
          call positn(xp,yp)
          call join  (xp+0.05,yp)
          CALL THICK(1)
        ELSE
c          xp = xpt + (in-1) * 0.95*xwid/ngs  + 0.025*xwid/ngs
c          xp = xpt + in * (xwid/(ngs+1)) - 0.05
          IF (grm_opt.EQ.1) THEN
            in1 = MOD(in,5)
            IF (in1.EQ.0) in1 = 5
            xp = xpt+(in1-1)*0.95*xwid/MIN(5,ngs)+0.025*xwid/MIN(5,ngs)
            yp = ypt+ ywid - 0.25 * ysep
            IF (ngs.GT.5.AND.in.GT.5) yp = yp - 0.25*ysep
          ELSE
            in1 = MOD(in,4)
            IF (in1.EQ.0) in1 = 4
            xp = xpt+(in1-1)*0.95*xwid/MIN(4,ngs)+0.025*xwid/MIN(4,ngs)

c
c jdemod - adjust for more than 2 lines of text
c
c           determine total rows and which row for specific label
c
            in2  = int((in-1)/4) + 1
            in2a = int((ngs-1)/4) + 1
             
            yp = ypt + ywid + (in2a-in2+1) * 0.22 * ysep - 0.06 * ysep
c
c            IF (ngs.GT.4.AND.in.LE.4) yp = yp + 0.25*ysep
c
c
c jdemod
c

          ENDIF

          IF (slopt.GT.0.AND.((slopt2.EQ.0.AND.plottype(in).GT. 1).OR.
     .                        (slopt2.EQ.2.AND.plottype(in).LT.-2)))
     .      THEN
            CALL HATOPT(1)
            CALL FILCOL(ABS(plottype(in))+ncols)
            CALL LINCOL(1)
            CALL FULL
            CALL BOX (xp+0.015-0.004,xp+0.015+0.004,
     .                yp+0.002-0.004,yp+0.002+0.004)
          ELSEIF (slopt.GT.0.AND.(slopt2.EQ.0.OR.slopt2.EQ.2).AND.
     .            plottype(in).LT.0) THEN
            IF     (plottype(in).EQ.-1) THEN
c...          Plots a '+':
              CALL LINCOL (1)
              CALL FULL
              CALL POSITN(xp+0.025      ,yp+0.002-0.004)
              CALL JOIN  (xp+0.025      ,yp+0.002+0.004)
              CALL POSITN(xp+0.025-0.004,yp+0.002)
              CALL JOIN  (xp+0.025+0.004,yp+0.002)
            ELSEIF (plottype(in).EQ.-2) THEN
c...          Plots a 'x':
              CALL LINCOL (1)
              CALL FULL
              CALL POSITN(xp+0.025-0.004,yp+0.002-0.004)
              CALL JOIN  (xp+0.025+0.004,yp+0.002+0.004)
              CALL POSITN(xp+0.025-0.004,yp+0.002+0.004)
              CALL JOIN  (xp+0.025+0.004,yp+0.002-0.004)
            ELSE
              CH = '+'

              CALL PCSCEN(xp+0.025,yp,ch)
            ENDIF
          ELSE
            call positn(xp     ,yp+0.002)
            call join  (xp+0.03,yp+0.002)
          ENDIF
        ENDIF


        CALL LinCol(1)

        IF (slopt2.EQ.1) THEN
          CALL CTRMAG(20)
          L = LENSTR(ELABS(IN))
          CALL PLOTCS(xp+0.07,yp,elabs(in)(5:L))
          CALL PLOTCS(xp+0.07,yp,elabs(in)(5:L))

c...      Print notes:
          CALL CTRMAG(15)
          DO j = 1, 10
            IF (char(j).NE.' ')
     .        CALL PLOTST(xp,0.40-(j-1)*0.03,char(j))
          ENDDO
        ELSE
          L = LENSTR(ELABS(IN))
          CALL PLOTCS(xp+0.04,yp,elabs(in)(5:L))
        ENDIF
c
c        xp = xpt + in * (xwid/(ngs+1)) - 0.05
c        yp = ypt + ywid + 0.25 * ysep
c        call positn(xp,yp)
c        call join(xp+0.03,yp)
c        call full
c        L = LENSTR(ELABS(IN))
c        CALL PLOTCS(xp+0.04,yp,elabs(in)(5:L))
c slmod end
c
      end do
c
      return
      end
c
c
c
      subroutine SLgrmaxes(ip,pltmin,pltmax,axmin,axmax,xlab,ylab,
     >                   sctype)
      implicit none
      integer ip,iten,sctype
      real pltmin,pltmax,axmin,axmax
      character*(*) xlab,ylab
c
c     GRMAXES: This routine plots the axes on the
c              particular box and sets the mapping between
c              the ND-space and vector-space in preparation
c              for the plotting.
c
c
      real xpt,ypt,xwid,ywid,xsep,ysep
      real tmin,tmax,power
      integer in,ik,l,lenstr,iexp
      external lenstr,iexp
c slmod begin
      INCLUDE 'params'
      INCLUDE 'slout'
      INCLUDE 'grminfo'
c slmod end
c
      call SLgrmfindxy(ip,xpt,ypt,xwid,ywid,xsep,ysep)
c
c     Draw Axes scales as appropriate.
c
      CALL LINCOL (1)
c slmod begin
      CALL THICK  (1)
c
c      CALL THICK  (2)
c slmod end
      call full
c
c     Xscale
c
      ITEN = IEXP(AXMIN,AXMAX)
c
c     Xaxis scale
c
      POWER = 10.0**(-ITEN)
      TMIN = AXMIN * POWER
      TMAX = AXMAX * POWER
      CALL PSPACE (xpt, xpt+xwid,ypt, ypt+ywid)
      write(6,*) 'xscale:',axmin,axmax,iten,power,tmin,tmax
      CALL MAP    (tmin, tmax, ypt, ypt+ywid)
      CALL CTRMAG (10)
c slmod begin
c      WRITE(0,*) ip,MOD(ip,4)

      IF (slopt2.EQ.1) CALL CTRMAG(13)

      IF (grm_opt.EQ.1.AND.MOD(ip,4).NE.0.AND.pageplots.EQ.8) GOTO 10

      IF (xlab(1:LEN_TRIM(xlab)).EQ.'none') GOTO 10
c slmod end
      CALL XSCALE
c
c     Xaxis labels
c
      CALL PSPACE (xpt, xpt+xwid,ypt, ypt+ywid)
      call map    (xpt, xpt+xwid,ypt, ypt+ywid)
c
c slmod begin
      IF (slopt2.EQ.1) THEN
        IF (ITEN .NE. 0) THEN
          CALL CTRMAG (15)
          CALL PLOTCS (xpt+0.85*xwid,ypt-0.6*ysep,'x10')
          CALL CTRMAG (15)
          CALL SUPFIX
          CALL TYPENI (ITEN)
          CALL NORMAL
        ENDIF
c       IPP/08 Krieger - revert to thin lines for the moment
c       CALL THICK(2)
        CALL THICK(1)
        CALL CTRMAG (20)
        L=LENSTR(xlab)
        CALL PCSCEN (xpt+0.5*xwid,ypt-0.9*ysep,xlab(1:L))
        CALL THICK(1)
      ELSE
        IF (ITEN .NE. 0) THEN
          CALL CTRMAG (8)
          CALL PLOTCS (xpt+0.85*xwid,ypt-0.9*ysep,'X10')
          CALL CTRMAG (10)
          CALL TYPENI (ITEN)
        ENDIF
        CALL CTRMAG (10)
        L=LENSTR(xlab)
        CALL PCSCEN (xpt+0.5*xwid,ypt-0.9*ysep,xlab(1:L))
      ENDIF
10    CONTINUE

      IF (ylab(1:LEN_TRIM(ylab)).EQ.'none') RETURN
c
c      IF (ITEN .NE. 0) THEN
c        CALL CTRMAG (8)
c        CALL PLOTCS (xpt+0.85*xwid,ypt-0.6*ysep,'X10')
c        CALL CTRMAG (10)
c        CALL TYPENI (ITEN)
c      ENDIF
c      CALL CTRMAG (10)
c      L=LENSTR(xlab)
c      CALL PCSCEN (xpt+0.5*xwid,ypt-0.6*ysep,xlab(1:L))
c10    CONTINUE
c slmod end
C
c     yscales
c
      if (sctype.eq.1.or.sctype.eq.2
     >    .or.sctype.eq.3.or.sctype.eq.4.or.sctype.eq.7) then
        ITEN = IEXP(PLTMIN, PLTMAX)
c
        POWER = 10.0**(-ITEN)
        TMIN = pltMIN * POWER
        TMAX = pltMAX * POWER
        CALL PSPACE (xpt, xpt+xwid,ypt, ypt+ywid)
        CALL MAP    (xpt, xpt+xwid,tmin,tmax)
        write(6,*) 'yscale:',pltmin,pltmax,iten,power,tmin,tmax
c slmod begin
        IF (slopt2.EQ.1) THEN
          CALL CTRMAG(13)
        ELSE
          CALL CTRMAG (10)
        ENDIF
c
c        CALL CTRMAG (10)
c slmod end
        CALL YSCALE
c
c     Y = 0 line if the scale crosses zero
c
        if (pltmin*pltmax.lt.0.0) then
           call full
           call thick(1)
           CALL PSPACE (xpt, xpt+xwid,ypt, ypt+ywid)
           CALL MAP    (axmin,axmax,pltmin,pltmax)
           call positn(axmin,0.0)
           call join(axmax,0.0)
        endif
      elseif (sctype.eq.5.or.sctype.eq.6) then
c
c     Logarithmic Y-axis
c
c        ITEN = IEXP(PLTMIN, PLTMAX)
c
c        POWER = 10.0**(-ITEN)
c        TMIN = pltMIN * POWER
c        TMAX = pltMAX * POWER
c        CALL PSPACE (xpt, xpt+xwid,ypt, ypt+ywid)
c        CALL MAP    (xpt, xpt+xwid,tmin,tmax)
c        write(6,*) 'yscall:',pltmin,pltmax,iten,power,tmin,tmax
c        CALL CTRMAG (10)
c        CALL YSCALL
c
        ITEN = IEXP(PLTMIN, PLTMAX)
c
        POWER = 10.0**(-ITEN)
        TMIN = pltMAX - 4.0
        TMAX = pltMAX
        CALL PSPACE (xpt, xpt+xwid,ypt, ypt+ywid)
        CALL MAP    (xpt, xpt+xwid,tmin,tmax)
        write(6,*) 'yscall:',pltmin,pltmax,iten,power,tmin,tmax
        CALL CTRMAG (10)
        CALL YSCALE
      endif
c
c     Y axis labels
c
      CALL PSPACE (xpt, xpt+xwid,ypt, ypt+ywid)
      CALL MAP    (xpt, xpt+xwid,ypt, ypt+ywid)
c
c     Rotate text 90 degrees
c
      CALL CTRORI (90.0)
c slmod begin
      IF (slopt2.EQ.1) THEN
        IF (ITEN .NE. 0) THEN
          CALL CTRMAG (20)
          CALL PLOTCS (xpt-1.0*xsep,ypt+0.90*ywid,'X10')
          CALL SUPFIX
          CALL TYPENI (ITEN)
        ENDIF
        CALL NORMAL
        CALL CTRMAG (20)
c       IPP/08 Krieger - revert to thin lines for the moment
c       CALL THICK(2)
        CALL THICK(1)
        L = LENSTR (YLAB)
        CALL PCSCEN (xpt-1.0*xsep,ypt+0.5*ywid,YLAB(:L))
        CALL THICK(1)
      ELSE
        IF (grm_opt.EQ.1.AND.pageplots.EQ.4) THEN
          IF (ITEN .NE. 0) THEN
            CALL CTRMAG (8)
            CALL PLOTCS (xpt-1.1*xsep,ypt+0.9*ywid,'X10')
            CALL CTRMAG (10)
            CALL TYPENI (ITEN)
          ENDIF
          CALL CTRMAG (10)
          L = LENSTR (YLAB2(ip))
          CALL PCSCEN (xpt-1.1*xsep,ypt+0.5*ywid,YLAB2(ip)(:L))
        ELSE
          IF (ITEN .NE. 0) THEN
            CALL CTRMAG (8)
            CALL PLOTCS (xpt-1.20*xsep,ypt+0.82*ywid,'X10')
            CALL CTRMAG (10)
            CALL TYPENI (ITEN)
          ENDIF
          CALL CTRMAG (10)
          L = LENSTR (YLAB)
          CALL PCSCEN (xpt-1.2*xsep,ypt+0.5*ywid,YLAB(:L))
        ENDIF
      ENDIF
c
c      IF (ITEN .NE. 0) THEN
c         CALL PLOTCS (xpt-0.6*xsep,ypt+0.90*ywid,'X10')
c         CALL CTRMAG (10)
c         CALL TYPENI (ITEN)
c      ENDIF
c      CALL CTRMAG (10)
c      L = LENSTR (YLAB)
c      CALL PCSCEN (xpt-0.6*xsep,ypt+0.5*ywid,YLAB(:L))
c slmod end
c
      if (sctype.eq.5.or.sctype.eq.6) then
         CALL PCSCEN (xpt-0.6*xsep,ypt+0.1*ywid,'(LOG)')
      endif
c
      CALL CTRORI (0.0)
c
      return
      end
c
c
c
      subroutine SLgrmdata (mvals,mouts,pnks,maxnms,
     >                    maxplts,maxngs,ip,ngs,
     >                    pltmin,pltmax,axmin,axmax,elabs,
c slmod begin
     >                    sctype,slopt,slopt2,slopt3,plottype)
c
c     >                    sctype)
c slmod end
      implicit none
      integer maxnms,maxplts,maxngs,ngs,ip,sctype
      integer pnks(maxplts)
      real pltmin,pltmax,axmin,axmax
      real mvals(maxnms,maxplts,maxngs)
      real mouts(maxnms,maxplts)
      character*(*) elabs(maxngs)
c slmod begin
      COMMON /RSEPCOM/ rsep
      REAL             rsep

      INTEGER slopt,slopt2,slopt3,plottype(MAXNGS)
c slmod end
c
c     GRMDATA: This routine plots the sepcified data
c              on one of the small plots on the page
c              specified by the ip value.
c
      real xpt,ypt,xwid,ywid,xsep,ysep
      real tmin,tmax
      integer in,ik
c
      include 'colours'
c      integer icol
c slmod begin
      REAL       LO
      PARAMETER (LO=1.0E-37)

      INTEGER i1
      REAL    x1,y1,awid,pwid
      CHARACTER*1 ch
c slmod end
c
c      INTEGER COLOUR(8)
c
c      DATA COLOUR /1,2,3,4,9,7,6,8/
c
      call SLgrmfindxy(ip,xpt,ypt,xwid,ywid,xsep,ysep)
c
      CALL PSPACE (xpt, xpt+xwid,ypt, ypt+ywid)
c
      if (sctype.eq.1.or.sctype.eq.2
     >    .or.sctype.eq.3.or.sctype.eq.4.or.sctype.eq.7) then
         CALL MAP  (axmin,axmax,pltmin,pltmax)
      elseif (sctype.eq.5.or.sctype.eq.6) then
         CALL MAPYL(axmin,axmax,pltmax-6.0,pltmax)
      endif
c
c slmod begin
c...  Ad hoc way of drawing the separatrix on the target
c     data plots:
      IF     (slopt3.EQ.4) THEN
        call thick(1)
        CALL BROKEN (3,3,3,3)
        call positn(1.0,pltmax)      
        call join  (1.0,pltmin)
      ELSEIF (slopt3.EQ.5) THEN
        call thick(1)
        CALL BROKEN (3,3,3,3)
        call positn(0.0,pltmax)      
        call join  (0.0,pltmin)
      ELSEIF (slopt3.EQ.6) THEN
        call thick(1)
        CALL BROKEN (3,3,3,3)
        call positn(rsep,pltmax)
        call join  (rsep,pltmin)
      ENDIF

c      CALL RGB
c      CALL COLSET(0.7,1.0,0.,8)
c      CALL COLSET(1.0,0.7,0.,9)
c slmod end
C
C
      call ctrmag(8)
      call thick(1)
c
      icol = 0
c
      do in = 1, ngs
c slmod begin
c        IF (pnks(in).EQ.0) CYCLE
        WRITE(6,*) 'PLOTTING: ',in,ngs
c slmod end
c
c       Set line colour
c
        icol = icol + 1
        IF (ICOL.GT.NCOLS) ICOL = 1
        CALL LINCOL (COLOUR(ICOL))
c
c       Set line type
c
        IF (in.LE.1) THEN
           CALL FULL
        ELSEIF (in.LE.5) THEN
           CALL BROKEN (3*in,2*in,3*in,2*in)
        ELSE
           CALL BROKEN (2*in,1*in,2*in,1*in)
        ENDIF
c slmod begin

        IF (slopt2.GT.0) THEN
          CALL FULL
          
          IF (plottype(icol).GT.1) THEN
c            CALL FULL
            CALL LinCol(colour(plottype(icol)+ncols))
          ENDIF
c         IPP/08 Krieger - revert to thin lines for the moment
c         IF (slopt2.EQ.1) CALL THICK(2)
          IF (slopt2.EQ.1) CALL THICK(1)
        ENDIF
c slmod end
c
c       Draw each plot
c
c slmod begin
        i1 = 1
        DO WHILE (mvals(i1,ip,in).EQ.LO.AND.i1.LT.pnks(ip))
          i1 = i1 + 1
        ENDDO
        call positn(mouts(i1,ip),mvals(i1,ip,in))
c
c        call positn(mouts(1,ip),mvals(1,ip,in))
c slmod end

c
c slmod begin
c        WRITE(0,*) '???',in,plottype(in)

        DO ik = i1, pnks(ip)

          IF (mvals(ik,ip,in).EQ.LO) CYCLE

          IF (slopt.GT.0.AND.((slopt2.EQ.0.AND.plottype(in).GT. 1).OR.
     .                        (slopt2.EQ.2.AND.plottype(in).LT.-2)))
     .      THEN
            CALL HATOPT (1)
            CALL FILCOL (ABS(plottype(in))+ncols)
            CALL LINCOL (1)
            CALL FULL

            x1 = mouts(ik,ip)
            y1 = mvals(ik,ip,in)

            IF (X1.GE.AXMIN.AND.X1.LE.AXMAX.AND.Y1.GE.pltMIN.AND.
     >          Y1.LE.PLTMAX.AND.Y1.NE.LO)
c     .        CALL BOX (x1-0.004*(AXMAX-AXMIN),
c     .                  x1+0.004*(AXMAX-AXMIN),
c     .                  y1-0.004*(PLTMAX-PLTMIN)*xwid/ywid,
c     .                  y1+0.004*(PLTMAX-PLTMIN)*xwid/ywid)
     .        CALL BOX (x1-0.008*(AXMAX-AXMIN),
     .                  x1+0.008*(AXMAX-AXMIN),
     .                  y1-0.008*(PLTMAX-PLTMIN)*xwid/ywid,
     .                  y1+0.008*(PLTMAX-PLTMIN)*xwid/ywid)
c              CALL LINCOL (1)
c              CALL FULL
c              CALL POSITN(x1-0.008*(AXMAX-AXMIN),y1)
c              CALL JOIN  (x1+0.008*(AXMAX-AXMIN),y1)
c              CALL POSITN(x1,y1-0.008*(PLTMAX-PLTMIN)*xwid/ywid)
c              CALL JOIN  (x1,y1+0.008*(PLTMAX-PLTMIN)*xwid/ywid)


c            write(6,90) 'plt 1:',in,ip,ik,
c     >               mouts(ik,ip),mvals(ik,ip,in)

          ELSEIF (slopt.GT.0.AND.(slopt2.EQ.0.OR.slopt2.EQ.2).AND.
     .            plottype(in).LT.0) THEN

            x1 = mouts(ik,ip)
            y1 = mvals(ik,ip,in)

            IF (plottype(in).EQ.-1) THEN
c...          Plots a '+':
              CALL LINCOL (1)
              CALL FULL
              CALL POSITN(x1-0.008*(AXMAX-AXMIN),y1)
              CALL JOIN  (x1+0.008*(AXMAX-AXMIN),y1)
              CALL POSITN(x1,y1-0.008*(PLTMAX-PLTMIN)*xwid/ywid)
              CALL JOIN  (x1,y1+0.008*(PLTMAX-PLTMIN)*xwid/ywid)
            ELSEIF (plottype(in).EQ.-2) THEN
c...          Plots a 'x':
              CALL LINCOL (1)
              CALL FULL
              awid = AXMAX  - AXMIN
              pwid = PLTMAX - PLTMIN
              CALL POSITN(x1-0.008*awid,y1-0.008*pwid*xwid/ywid)
              CALL JOIN  (x1+0.008*awid,y1+0.008*pwid*xwid/ywid)
              CALL POSITN(x1-0.008*awid,y1+0.008*pwid*xwid/ywid)
              CALL JOIN  (x1+0.008*awid,y1-0.008*pwid*xwid/ywid)
            ELSE
c...          Slight error where this is plotted (a tiny bit too "high" on plots):
              CH = '+'

              CALL PLOTST (mouts(ik,ip),mvals(ik,ip,in),ch)
 
c            write(6,90) 'plt 2:',in,ip,ik,
c     >               mouts(ik,ip),mvals(ik,ip,in)
            ENDIF

          ELSE
            call join(mouts(ik,ip),mvals(ik,ip,in))
            IF (IK.EQ.3.OR.IK.EQ.pnks(ip)-5)
     >         CALL PLOTST (mouts(ik,ip),mvals(ik,ip,in),
     >                      elabs(in)(1:4))
c            write(6,90) 'plt 3:',in,ip,ik,
c     >               mouts(ik,ip),mvals(ik,ip,in)
90          FORMAT(A,3I6,1P,2E10.2,0P)
          ENDIF
        end do
c
      end do
c
      call full
c
      return
      end
c
c
c
      subroutine SLgrmfindxy (boxindex,xpt,ypt,xwid,ywid,
     >                      xsep,ysep)
      implicit none
      integer boxindex
      real xpt,ypt,xwid,ywid,xsep,ysep
c
      include 'grminfo'
c slmod begin
      INCLUDE 'params'
      INCLUDE 'slout'

      INTEGER overridecount
      LOGICAL override
      CHARACTER*256 dummy,cdum1

      DATA override,overridecount /.FALSE.,0/
c slmod end
c
c     GRMFINDXY: This routine finds the X,Y co-ordinates and
c                size of the box to be used for the particular
c                plot.
c
      integer ix,iy,n,xplts,yplts,plotindex
      real    xpts(2),ypts(4),xwid0,ywid0,xsep0,ysep0
c
c     The following describes the box location in ND space
c     X,Y are the location of the lower left corner of
c     the box.
c     Xwid and Ywid are the X and Y extent of the box
c     Xsep and Ysep are the anount of X and Y border
c     around each box.
c
      integer in
c
c slmod begin
      IF (override) THEN
        WRITE(0,*) 'overdire:',map1x
        xpt = map1x
        ypt = map1y
        xwid = map2x - map1x
        ywid = map2y - map1y
        xsep = 0.05
        ysep = 0.05

c... HAHAHAHAH, WHAT CRAP! GETTING DESPIRATE!
        overridecount = overridecount + 1
        IF (overridecount.EQ.3) THEN
          override = .FALSE.
          overridecount = 0
        ENDIF
        RETURN
      ELSE
        READ(5,'(A256)') dummy
        IF   (dummy(8:11).EQ.'Size'.OR.dummy(8:11).EQ.'size'.OR.
     .        dummy(8:11).EQ.'SIZE') THEN
          READ(dummy,*) cdum1,map1x,map2x,map1y,map2y
          override = .TRUE.
        ELSE
          BACKSPACE 5
        ENDIF
      ENDIF

      xsep0 = 0.05
      ysep0 = 0.05
c
c      data    xsep0 /0.05/
c      data    ysep0 /0.05/
c slmod end
c
c     Draw Box - Boxindex is 1 to pageplots
c     Modify so that this is so.
c
      n = (boxindex -1)/ pageplots
      plotindex = boxindex - n * pageplots
c slmod begin
      IF (pageplots.EQ.1) THEN
        xplts = 1
        IF (slopt2.EQ.1) THEN
          xsep0 = 0.10
          ysep0 = 0.10
        ENDIF
      ELSE
        xplts = 2
      ENDIF
c
c      xplts = 2
c slmod end
      yplts = (pageplots-1) /2 + 1
c
      xwid0 = (1.0 - xplts * 2.0 * xsep0) / xplts
c slmod begin
      IF (grm_opt.EQ.1) THEN
        IF     (pageplots.EQ.8) THEN
          ywid0 = (1.0 - 2.0 * ysep0 - 0.01 * yplts) / yplts
        ELSEIF (pageplots.EQ.4) THEN
          ywid0 = (1.0 - 2.0 * ysep0 - 0.04 * yplts) / yplts
        ELSEIF (pageplots.EQ.2) THEN
          ywid0 = (1.0 - 2.0 * ysep0 - 0.04 * yplts) / yplts
        ENDIF
      ELSE
        ywid0 = (1.0 - yplts * 2.0 * ysep0 - 0.5 * ysep0) / yplts
      ENDIF
c
c      ywid0 = (1.0 - yplts * 2.0 * ysep0) / yplts
c slmod end
c
      do in = 1,xplts
c slmod begin
         xpts(in) = 2.0 * xsep0 + (in-1) * (xwid0+2.0*xsep0)
c
c         xpts(in) = xsep0 + (in-1) * (xwid0+2.0*xsep0)
c slmod end
      end do
c
      do in = yplts,1,-1
c slmod begin
         IF (grm_opt.EQ.1) THEN
           IF     (pageplots.EQ.8) THEN
             ypts(in) = ysep0 + (yplts-in) * (ywid0 + 0.01)
           ELSEIF (pageplots.EQ.4) THEN
             ypts(in) = 1.2 * ysep0 + (yplts-in) * (ywid0 + 1.5 * ysep0)
           ELSEIF (pageplots.EQ.2) THEN
             ypts(in) = 1.2 * ysep0 + (yplts-in) * (ywid0 + 1.5 * ysep0)
           ENDIF
         ELSE
           ypts(in) = 1.2 * ysep0 + (yplts-in)*(ywid0+2.0*ysep0)
         ENDIF
c
c         ypts(in) = ysep0+(yplts-in)*(ywid0+2.0*ysep0)
c slmod end
      end do
c
c     iy is 1 to pageplots/2
c     ix is 1 to 2
c
c slmod begin - a
      IF (grm_opt.EQ.1) THEN
        IF     (pageplots.EQ.8) THEN
          ix = ((plotindex-1) / 4) + 1
          iy = plotindex - ((plotindex - 1) / 4) * 4
        ELSEIF (pageplots.EQ.4) THEN
          ix = ((plotindex-1) / 2) + 1
          iy = plotindex - ((plotindex - 1) / 2) * 2
        ELSEIF (pageplots.EQ.2) THEN
          ix = ((plotindex-1) / 1) + 1
          iy = plotindex - ((plotindex - 1) / 1) * 1
        ELSE
          STOP 'WHOA! NOT READY FOR THIS FROM 978'
        ENDIF

c        WRITE(0,*) plotindex,ix,iy
      ELSE
        ix = 1
        if (float(plotindex/2).eq.float(plotindex)/2.0) ix = 2

        iy =  (plotindex-1)/2 + 1
      ENDIF
c
c      ix = 1
c      if (float(plotindex/2).eq.float(plotindex)/2.0) ix = 2
c
c      iy =  (plotindex-1)/2 + 1
c slmod end
c
c slmod begin
      IF (slopt2.EQ.1) THEN
        xpt = xpts(ix) - 0.06
        ypt = ypts(iy)
      ELSE
        xpt = xpts(ix)
        ypt = ypts(iy)
      ENDIF
c
c      xpt = xpts(ix)
c      ypt = ypts(iy)
c slmod end
c
      xwid = xwid0
      ywid = ywid0
c
      xsep = xsep0
      ysep = ysep0
c
      write (6,*) 'GRMFind:',plotindex,boxindex,n,iy,ix,
     >                   xpt,ypt,xwid,
     >                   ywid,xsep,ysep,pageplots
c
      return
      end
c
c ======================================================================
c
c
c ======================================================================
c
c subroutine: SLSET
c
      SUBROUTINE SLSET(x1,x2,y1,y2,min1,max1,min2,max2)

      IMPLICIT none

      REAL x1,x2,y1,y2,min1,max1,min2,max2

      INCLUDE 'params'
      INCLUDE 'slcom'
      INCLUDE 'slout'

      include 'comgra'    
      include 'colours'

c      COMMON /COMGRA/ CXMIN,CXMAX,CYMIN,CYMAX,IPLOTS,ICOL,NPLOTS,ISPOT
c      REAL            CXMIN,CXMAX,CYMIN,CYMAX
c      INTEGER         IPLOTS,ICOL,NPLOTS,ISPOT


      map1x = x1
      map2x = x2
      map1y = y1
      map2y = y2

      cxmin = min1
      cxmax = max1
      cymin = min2
      cymax = max2

      iplots = 1
      icol   = 1
      ispot  = 12

      plottype(1) = 1

      RETURN
      END













