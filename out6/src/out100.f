c     -*-Fortran-*-
c
      subroutine out100(iref,graph,iopt,ierr)
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
      include 'dynam2'
      include 'dynam3'
      include 'dynam4'
c      include 'pindata'
c      include 'cadas'
c      include 'grbound'
      include 'outxy'
      include 'cedge2d'
c slmod being
      include 'colours'
c slmod end
c      include 'transcoef'
c      include 'cioniz'
c      include 'reiser' 
c      include 'printopt' 
c
c     Local Variables
c

      integer startin,endin,stepin
      integer ik,ir,iz
      integer jk,jr,jz
      integer id,in,ii,jd
      real r,z
      integer ig

      real dist
c
c     PEC Atomic data plots  
c 
      logical sxb
      real coef,coef2(20),nad(20),tad(20)
      character*2 yeart 
      integer iclass 
c

      integer lw,uw,iw
      integer lt,ut,it
      integer cnt
      integer tiz 

      real VAL,VALTOT



      character refinf*9

c
      real    kvalmin,minkval
c
c     Temp array for actual powls
c
      real tmppowls(maxnks,maxnrs,-1:maxizs)


c
c     Array for storing calculated TAUP's
c
      real    taup(maxnks,maxnrs),nregion(maxnks,maxnrs)
      real    ntotal,norm_fact


      integer iz1,iz2

      real   tbrem(20),brempec(20)



      integer   cz
      real      den
      REAL PLRPAD(MAXNKS,MAXNRS)
      REAL PLRPAD2(MAXNKS,MAXNRS)

      REAL ROUTS(MAXNRS),RWIDS(MAXNRS),RVALS(MAXNRS,MAXNGS)
      REAL LEVEL
      REAL KTOTA1


      IF (IREF.LE.170.or.(iref.ge.190.and.iref.lt.200)) THEN
        IF (IOPT.LT.-2.OR.IOPT.GT.NIZS+1) return

      ELSEIF (iref.eq.174.or.IREF.eq.176.or.iref.eq.177) THEN
         IF (IOPT.LT.0.OR.IOPT.GE.80) return

      else
        if (iopt.eq.0) return

      endif

      call init_plot(iref,graph,iopt) 
c     IPP/08 Krieger - added this because original color opt might 
c     get lost by other plots
      call setup_col(n_cols,col_opt)


C
C-----------------------------------------------------------------------
C     CONTOURS OF STEADY-STATE DENSITY
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.101.OR.IREF.EQ.102) THEN
        XLAB   = '   R  (M)'
        YLAB   = '   Z  (M)'
        NGS    = cngs
        IZ     = IOPT
        REF    = 'DENSITY ' // ZLABS(IZ)(5:11) // XPOINT
        write (6,*) 'range:',xxmin,xxmax,yymin,yymax,zminp,zmaxp
        WRITE (IPLOT,9012) NPLOTS,REF
c
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,SDLIMS,IZ+2,NIZS+2,MAXIZS+2,FT,FP,
     >                scalef,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
C
C-----------------------------------------------------------------------
C     CONTOURS OF TIME-DEPENDENT DENSITY
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.103.OR.IREF.EQ.104) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NGS    = CNGS
        IZ     = IOPT
        DO 1040 IT = 1, NTS
          WRITE (NAME,'(1P,E8.1,A)') DWELTS(IZ)*DWELFS(IT),'S'
          len = lenstr(name) 
          REF  = 'DENSITY ' // ZLABS(IZ)(5:11) // NAME(1:len) // XPOINT
          WRITE (IPLOT,9012) NPLOTS,REF
          CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >      YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
          CALL CONTOUR (ICNTR,NGS,LIMS(1,1,-1,IT),IZ+2,NIZS+2,MAXIZS+2,
     >                  FT,FP,scalef,
     >                  XOUTS,1,NXS,YOUTS,1,NYS,
     >                  XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
          IF (IT.NE.NTS) NPLOTS = NPLOTS + 1
 1040   CONTINUE
      ENDIF
C
C-----------------------------------------------------------------------
C     CONTOURS OF STEADY-STATE IMPURITY DENSITY / NE  = Fz
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.105.OR.IREF.EQ.106) THEN
        XLAB   = '   R  (M)'
        YLAB   = '   Z  (M)'
        NGS    = cngs
        IZ     = IOPT
        REF    = 'FZ ' // ZLABS(IZ)(5:11) // XPOINT
        write (6,*) 'range:',xxmin,xxmax,yymin,yymax,zminp,zmaxp
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
c
c       Calculate Fz for requested charge state.
c
        call rzero(plastmp,maxnks*maxnrs)
c
        if (iz.gt.nizs) then
           iz1 = 0
           iz2 = nizs
        else
           iz1 = iz
           iz2 = iz
        endif
c
        do iz = iz1,iz2
           do ir = 1,nrs
              do ik = 1,nks(ir)
                 if (knbs(ik,ir).ne.0.0) then
                    plastmp(ik,ir) = plastmp(ik,ir)
     >                    + sdlims(ik,ir,iz)/knbs(ik,ir)
                 endif
              end do
           end do
        end do
c
        CALL CONTOUR (ICNTR,NGS,PLASTMP,1,1,1,FT,FP,scalef,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
C
C-----------------------------------------------------------------------
C     CONTOURS OF STEADY-STATE DENSITY - FLUID CODE RESULTS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.107.OR.IREF.EQ.108) THEN
        XLAB   = '   R  (M)'
        YLAB   = '   Z  (M)'
        NGS    = cngs
        IZ     = IOPT
        REF    = 'FLUID CODE - DENS ' // ZLABS(IZ)(5:11) // XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,e2dnzs,IZ+1,cre2dizs+1,MAXe2dIZS+1,
     >                FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
c
c kkmod
c
C
C-----------------------------------------------------------------------
C     CONTOURS OF STEADY-STATE IMPURITY TEMPERATURE, Krieger IPP/97
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.109.OR.IREF.EQ.110) THEN
        XLAB   = '   R  (M)'
        YLAB   = '   Z  (M)'
        IZ     = IOPT
        REF    = 'TEMPERATURE ' // ZLABS(IZ)(5:11) // XPOINT
        write (6,*) 'range:',xxmin,xxmax,yymin,yymax,zminp,zmaxp
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,SDTS,IZ+2,NIZS+2,MAXIZS+2,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
c
c kkmod
c
C
C-----------------------------------------------------------------------
C     TIME-DEP & STEADY STATE DENSITY AGAINST S FOR GIVEN CONTOUR
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.111.OR.IREF.EQ.112.OR.IREF.EQ.113.OR.IREF.EQ.114) THEN
        NPLOTS = NPLOTS - 1
        XLAB   = '   S  (M)  '
        YLAB   = '           '
        IZ     = IOPT
        READ (GRAPH(38:44),'(I4)') IR
        ISMOTH = 99
        IF (IREF.EQ.111.OR.IREF.EQ.112) THEN
          LT = 1
          UT = NTS
          IF (IREF.EQ.112) ISMOTH = 0
        ELSE
          LT = NTS+1
          UT = NTS+1
          IF (IREF.EQ.114) ISMOTH = 0
        ENDIF
        DO 1140 IT = LT, UT
          CALL rzero (kvals, maxnks*maxngs)
          IF (IT.LE.NTS) THEN
            WRITE (NAME,'(A,1P,E8.1,A)') ' AT',DWELTS(IZ)*DWELFS(IT),'S'
            FACT = FACTA(IZ)
          ELSE
            NAME = ' INTEGRATED'
            FACT = FACTB(IZ)
            write(6,*) 'fact:',iz,factb(iz),facta(iz)
          ENDIF
          WRITE (KLAB,'(''ALONG RING'',I3,'', K='',F10.6)') IR,KKS(IR)
          DO 1130 IK = 1, NKS(IR)
            KOUTS(IK) = KSS(IK,IR)
            KWIDS(IK) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
            IF (IT.LE.NTS) THEN
              IF (IZ.EQ.NIZS+1) THEN
                DO 1110 JZ = 1, NIZS
                  KVALS(IK,1) = KVALS(IK,1) + LIMS(IK,IR,JZ,IT)
 1110           CONTINUE
              ELSEIF (IZ.EQ.-2) THEN
                KVALS(IK,1) = FT*LIMS(IK,IR,0,IT) - FP*LIMS(IK,IR,-1,IT)
              ELSE
                KVALS(IK,1) = LIMS(IK,IR,IZ,IT)
              ENDIF
            ELSE
              IF (IZ.EQ.NIZS+1) THEN
                DO 1120 JZ = 1, NIZS
                  KVALS(IK,1) = KVALS(IK,1) + SDLIMS(IK,IR,JZ)
 1120           CONTINUE
              ELSEIF (IZ.EQ.-2) THEN
                KVALS(IK,1) = FT*SDLIMS(IK,IR,0) - FP*SDLIMS(IK,IR,-1)
              ELSE
                KVALS(IK,1) = SDLIMS(IK,IR,IZ)
              ENDIF
            ENDIF
c
*           normalization below seems strange, deleted, Krieger IPP/98
*           KVALS(IK,1) = FACT / KWIDS(IK) * KVALS(IK,1)
c
            KVALS(IK,1) = KVALS(IK,1) * absfac
            write(iplot,"(1x,f7.4,3x,1p,e9.3)")
     >         kouts(ik),kvals(ik,1)
 1130     CONTINUE
          REF = 'DENSITY ' // ZLABS(IZ)(5:11) // NAME
          NPLOTS = NPLOTS + 1
          WRITE (IPLOT,9012) NPLOTS,REF
          ELABS(1) = ' '
          CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >      1,ISMOTH,KOUTS(1),KOUTS(NKS(IR)),0.0,HI,IGNORS,ITEC,AVS,
     >      NAVS,JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,KLAB,TABLE,
     >      2,2,1.0,0)
c
          if (clsup.eq.1) then
c
          CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >      1,ISMOTH,KOUTS(1),mgst*KOUTS(NKS(IR)),0.0,HI,
     >      IGNORS,ITEC,AVS,
     >      NAVS,JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,KLAB,TABLE,
     >      2,2,1.0,0)

c         IPP/01 Krieger - missing plot added
          CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >      1,ISMOTH,mgnd*KOUTS(NKS(IR)),KOUTS(NKS(IR)),0.0,HI,
     >      IGNORS,ITEC,AVS,
     >      NAVS,JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,KLAB,TABLE,2,2,1.0)
c
          endif
c
 1140   CONTINUE
      ENDIF
C
C-----------------------------------------------------------------------
C     PLRP CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.121.OR.IREF.EQ.122) THEN
       XLAB   = '   X/A  (M)'
       YLAB   = '   Y/A  (M)'
       NGS    = CNGS
       IZ     = IOPT
       TIZ    = IOPT
C
C      NEED TO TREAT A REQUEST FOR SECONDARY
C      NEUTRALS DIFFERENTLY BECAUSE THE PLRPS NEUTRALS NO LONGER
C      NECESSARILY CORRESPOND TO THE FIRST TWO ARRAY ENTRIES.
C
       IF (IZ.EQ.-2) THEN
         PSHIFT = PNCNT
         IF (PLAMS(-1).GT.0.0) THEN
           REF = 'PLRP ' // PLABS(-2)(5:11) // XPOINT
           WRITE (IPLOT,9012) NPLOTS,REF
c
c           CALL RVALUE (CVALS,PLRPS,IZ+2,PLRPCNT+2,MAXPLRP+2,FT,FP,
c     >       IXMIN,IXMAX,IYMIN,IYMAX,VMIN,VMAX)
c
           CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >       YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
c
c           CALL SUPIMP ('SELECT')
c           DO 1205 IG = 1, NGS-1
c             LEVEL = VMIN+REAL(IG)/REAL(NGS)/REAL(NGS+1-IG) *
c     >          (VMAX-VMIN)
c             if (ig.eq.1.and.level.gt.(vmin+0.01*(vmax-vmin)))
c     >           level = vmin + 0.01 * (vmax-vmin)
c             WRITE (NAME,'(4X,1P,E8.1)') LEVEL
c             IF (IG.EQ.NGS) NAME(13:18) = ' (MAX)'
c             CALL GRCONT (CVALS,IXMIN,IXMAX,MAXGXS,IYMIN,IYMAX,
c     >                  MAXGYS,LEVEL,XOUTS,YOUTS,NAME)
c 1205      CONTINUE
c           CALL FRAME
c
c kkmod
c
c          changed to updated contour routine; Krieger IPP/97
c
           DO IR = 1,NRS
             DO IK = 1,NKS(IR)
              PLASTMP(IK,IR)=FT*plrps(IK,IR,pncnt-2)-FP*plrps(IK,IR,-1)
             enddo
           enddo
           CALL CONTOUR (ICNTR,NGS,PLASTMP,1,1,1,FT,FP,1.0,
     >                   XOUTS,1,NXS,YOUTS,1,NYS,
     >                   XXMIN,XXMAX,YYMIN,YYMAX,
     >                   nconts,conts,cntropt,minscale,maxscale)
c
c kkmod
c
         ENDIF
         PSHIFT = 1
       ELSE
         DO 1220 IN = -1,PLRPCNT
           IF (PIZS(IN).EQ.TIZ) THEN
C
C          IF THE LINE MATCHES THE REQUESTED IONIZATION STATE THEN
C          PLOT IT.
C
C          SET IZ TO THE PROPER INDEX FOR SPECIFIC PLRP
C
           IZ = IN
C
           IF (PLAMS(IZ).GT.0.0) THEN
             REF = 'PLRP ' // PLABS(IZ)(5:11) // XPOINT
             WRITE (IPLOT,9012) NPLOTS,REF
c
c             CALL RVALUE (CVALS,PLRPS,IZ+2,PLRPCNT+2,MAXPLRP+2,
c     >         FT,FP,
c     >         IXMIN,IXMAX,IYMIN,IYMAX,VMIN,VMAX)

             CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >         YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)

c
c             CALL SUPIMP ('SELECT')
c             DO 1210 IG = 1, NGS-1
c               LEVEL = VMIN+REAL(IG)/REAL(NGS)/REAL(NGS+1-IG) *
c     >                 (VMAX-VMIN)
c               if (ig.eq.1.and.level.gt.(vmin+0.01*(vmax-vmin)))
c     >             level = vmin + 0.01 * (vmax-vmin)
c               WRITE (NAME,'(4X,1P,E8.1)') LEVEL
c               IF (IG.EQ.NGS) NAME(13:18) = ' (MAX)'
c               CALL GRCONT (CVALS,IXMIN,IXMAX,MAXGXS,IYMIN,IYMAX,
c     >                      MAXGYS,LEVEL,XOUTS,YOUTS,NAME)
c 1210        CONTINUE
c             CALL FRAME
c
c
c kkmod
c
c            changed to updated contour routine; Krieger IPP/97
c
             DO IR = 1,NRS
               DO IK = 1,NKS(IR)
                 PLASTMP(IK,IR)=plrps(IK,IR,iz)
               enddo
             enddo
             CALL CONTOUR (ICNTR,NGS,PLASTMP,1,1,1,FT,FP,1.0,
     >                     XOUTS,1,NXS,YOUTS,1,NYS,
     >                     XXMIN,XXMAX,YYMIN,YYMAX,
     >                     nconts,conts,cntropt,minscale,maxscale)
c
c kkmod
c
            ENDIF
           ENDIF
 1220    CONTINUE
        ENDIF
      ENDIF
C
C-----------------------------------------------------------------------
C     PLRP RATIO CONTOURS USING ADAS ATOMIC PHYSICS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.123.OR.IREF.EQ.124) THEN
        XLAB   = '     R  (M)'
        YLAB   = '     Z  (M)'
        PLANE    = 'RATIO PLOT'
        NGS    = CNGS

c
c       Load data for the first line
c
c       Verify that it is an acceptable line
c
        if ((z_atom.eq.cion.and.iz_state.ge.0.and.iz_state.le.nizs).or.
     >      (z_atom.eq.cizb.and.iz_state.ge.0.and.iz_state.le.cizb))
     >       then
c
          call LDADAS(z_atom,IZ_state,ADASID,ADASYR,ADASEX,
     >                  ISELE,ISELR,ISELX,
     >                  plrpad,Wlngth,IRCODE)
c
          IF (IRCODE.NE.0) THEN
             WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
             return
          ENDIF
c
          REF = 'ADAS Z=XX XX XXXXX ('
          WRITE(ref(8:9),'(I2)') Z_atom
          WRITE(REF(11:12),'(I2)') IZ_state
          WRITE(REF(14:18),'(I5)') NINT(WLNGTH)
          LEN = LENSTR(REF)
          IF (ISELE.GT.0) REF = REF(1:LEN) // 'E'
          LEN = LENSTR(REF)
          IF (ISELR.GT.0) REF = REF(1:LEN) // 'R'
          LEN = LENSTR(REF)
          IF (ISELX.GT.0) REF = REF(1:LEN) // 'C'
          LEN = LENSTR(REF)
          REF = REF(1:LEN) // ') '
          LEN = LENSTR(REF)
          REF = REF(1:LEN) // XPOINT
          ref = 'NUMER:'//ref
          WRITE (IPLOT,9012) NPLOTS,REF
c
        else
          WRITE(6,*) 'PLOT 123, 124 - INVALID ION AND CHARGE STATE'//
     >               ' SPECIFIED WITH FIRST ADAS INPUT LINE'
          return
        endif
c
c
c       Load data for the second line
c
c       Verify that it is an acceptable line
c
        if ((z_atom2.eq.cion.and.iz_state2.ge.0.and.iz_state2.le.nizs)
     >            .or.
     >      (z_atom2.eq.cizb.and.iz_state2.ge.0.and.iz_state2.le.cizb))
     >       then
c
          call LDADAS(z_atom2,IZ_state2,ADASID2,ADASYR2,ADASEX2,
     >                  ISELE2,ISELR2,ISELX2,
     >                  plrpad2,Wlngth,IRCODE)
c
          IF (IRCODE.NE.0) THEN
             WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
             return
          ENDIF
c
          anly = 'ADAS Z=XX XX XXXXX ('
          WRITE(anly(8:9),'(I2)') Z_atom2
          WRITE(anly(11:12),'(I2)') IZ_state2
          WRITE(anly(14:18),'(I5)') NINT(WLNGTH)
          LEN = LENSTR(anly)
          IF (ISELE2.GT.0) anly = anly(1:LEN) // 'E'
          LEN = LENSTR(anly)
          IF (ISELR2.GT.0) anly = anly(1:LEN) // 'R'
          LEN = LENSTR(anly)
          IF (ISELX2.GT.0) anly = anly(1:LEN) // 'C'
          LEN = LENSTR(anly)
          anly = anly(1:len) // ') '
          LEN = LENSTR(anly)
          anly = anly(1:LEN) // XPOINT
          anly = 'DENOM:'//anly
          WRITE (IPLOT,9012) NPLOTS,anly
c
        else
          WRITE(6,*) 'PLOT 123, 124 - INVALID ION AND CHARGE STATE'//
     >               ' SPECIFIED WITH SECOND ADAS INPUT LINE'
          return
        endif
c
c       Calculate RATIO and call contour plotting routine to plot it -
c       if the denominator contribution is 0.0 then the ratio is not calculated
c       for that cell.
c
        do ir = 1,nrs
           do ik = 1,nks(ir)
              if (plrpad2(ik,ir).gt.0.0) then
                 plrpad(ik,ir) = plrpad(ik,ir) / plrpad2(ik,ir)
              else
                 plrpad(ik,ir) = 0.0
              endif
           end do
        end do
C
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >      YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,PLRPAD,1,1,1,FT,FP,1.0,
     >                  XOUTS,1,NXS,YOUTS,1,NYS,
     >                  XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
C
C-----------------------------------------------------------------------
C     PLRP CONTOURS USING ADAS ATOMIC PHYSICS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.125.OR.IREF.EQ.126) THEN
        XLAB   = '     R  (M)'
        YLAB   = '     Z  (M)'
        PLANE    = 'PHOTONS/M^3/S'
        NGS    = CNGS
        IZ     = IOPT
        IF (IZ.GE.0 .AND. IZ.LE.NIZS .AND. IZ.LT.CION) THEN
c
          call LDADAS(CION,IZ,ADASID,ADASYR,ADASEX,ISELE,ISELR,ISELX,
     >                  plrpad,Wlngth,IRCODE)
c
          IF (IRCODE.NE.0) THEN
             WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
             RETURN
          ENDIF
C
          REF = 'ADAS PLRP XX XXXXX ('
          WRITE(REF(11:12),'(I2)') IZ
          WRITE(REF(14:18),'(I5)') NINT(WLNGTH)
          LEN = LENSTR(REF)
          IF (ISELE.GT.0) REF = REF(1:LEN) // 'E'
          LEN = LENSTR(REF)
          IF (ISELR.GT.0) REF = REF(1:LEN) // 'R'
          LEN = LENSTR(REF)
          IF (ISELX.GT.0) REF = REF(1:LEN) // 'C'
          LEN = LENSTR(REF)
          REF = REF(1:LEN) // ') '
          LEN = LENSTR(REF)
          REF = REF(1:LEN) // XPOINT
          WRITE (IPLOT,9012) NPLOTS,REF
C
          CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >      YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
          CALL CONTOUR (ICNTR,NGS,PLRPAD,1,1,1,FT,FP,1.0,
     >                  XOUTS,1,NXS,YOUTS,1,NYS,
     >                  XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
        ENDIF
      ENDIF
C
C-----------------------------------------------------------------------
C     H PLRP CONTOURS USING ADAS ATOMIC PHYSICS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.127.OR.IREF.EQ.128) THEN
        XLAB   = '     R  (M)'
        YLAB   = '     Z  (M)'
        PLANE    = 'PHOTONS/M^3/S'
        NGS    = CNGS
        IZ     = IOPT
        IF (IZ.EQ.0) THEN
c
          call LDADAS(1,IZ,ADASID,ADASYR,ADASEX,ISELE,ISELR,ISELX,
     >                  plrpad,Wlngth,IRCODE)
c
          IF (IRCODE.NE.0) THEN
             WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
             RETURN
          ENDIF
C
          REF = 'ADAS H PLRP XX XXXXX ('
          WRITE(REF(13:14),'(I2)') IZ
          WRITE(REF(16:20),'(I5)') NINT(WLNGTH)
          LEN = LENSTR(REF)
          IF (ISELE.GT.0) REF = REF(1:LEN) // 'E'
          LEN = LENSTR(REF)
          IF (ISELR.GT.0) REF = REF(1:LEN) // 'R'
          LEN = LENSTR(REF)
          IF (ISELX.GT.0) REF = REF(1:LEN) // 'C'
          LEN = LENSTR(REF)
          REF = REF(1:LEN) // ') '
          LEN = LENSTR(REF)
          REF = REF(1:LEN) // XPOINT
          WRITE (IPLOT,9012) NPLOTS,REF
C
          CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >      YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
          CALL CONTOUR (ICNTR,NGS,PLRPAD,1,1,1,FT,FP,1.0,
     >                  XOUTS,1,NXS,YOUTS,1,NYS,
     >                  XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
        ENDIF
      ENDIF
C
C-----------------------------------------------------------------------
C     IONISATION CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.131.OR.IREF.EQ.132) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NGS    = cngs
        IZ     = IOPT
        REF = 'IONISATION ' // ZLABS(IZ)(5:11) // XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,TIZS,IZ+2,NIZS+2,MAXIZS+2,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
C
C-----------------------------------------------------------------------
C     PLRP CONTOURS USING ADAS ATOMIC PHYSICS - Photons/S
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.135.OR.IREF.EQ.136) THEN
        XLAB   = '     R  (M)'
        YLAB   = '     Z  (M)'
        PLANE    = 'PHOTONS/S'
        NGS    = CNGS
        IZ     = IOPT
        IF (IZ.GE.0 .AND. IZ.LE.NIZS .AND. IZ.LT.CION) THEN
c
          call LDADAS(CION,IZ,ADASID,ADASYR,ADASEX,ISELE,ISELR,ISELX,
     >                  plrpad,Wlngth,IRCODE)
c
          IF (IRCODE.NE.0) THEN
             WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
             RETURN
          ENDIF
c
c         Rescale data to represent absolute number instead of density
c
          do ir = 1,nrs
             do ik = 1,nks(ir)
                plrpad(ik,ir) = plrpad(ik,ir) * kareas(ik,ir)
             end do
          end do
c
          REF = 'ADAS PLRP XX XXXXX ('
          WRITE(REF(11:12),'(I2)') IZ
          WRITE(REF(14:18),'(I5)') NINT(WLNGTH)
          LEN = LENSTR(REF)
          IF (ISELE.GT.0) REF = REF(1:LEN) // 'E'
          LEN = LENSTR(REF)
          IF (ISELR.GT.0) REF = REF(1:LEN) // 'R'
          LEN = LENSTR(REF)
          IF (ISELX.GT.0) REF = REF(1:LEN) // 'C'
          LEN = LENSTR(REF)
          REF = REF(1:LEN) // ') '
          LEN = LENSTR(REF)
          REF = REF(1:LEN) // XPOINT
          WRITE (IPLOT,9012) NPLOTS,REF
c
          if (cgrprint.eq.1) then
             CALL PRRMATDIV(plrpad,MAXNKS,nks(irsep),NRS,49,
     >                      REF)
          endif
C
          CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >      YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
          CALL CONTOUR (ICNTR,NGS,PLRPAD,1,1,1,FT,FP,1.0,
     >                  XOUTS,1,NXS,YOUTS,1,NYS,
     >                  XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
        ENDIF
      ENDIF
C
C-----------------------------------------------------------------------
C     H PLRP CONTOURS USING ADAS ATOMIC PHYSICS - Photons/S
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.137.OR.IREF.EQ.138) THEN
        XLAB   = '     R  (M)'
        YLAB   = '     Z  (M)'
        PLANE    = 'PHOTONS/S'
        NGS    = CNGS
        IZ     = IOPT
        IF (IZ.EQ.0) THEN
c
          call LDADAS(1,IZ,ADASID,ADASYR,ADASEX,ISELE,ISELR,ISELX,
     >                  plrpad,Wlngth,IRCODE)
c
          IF (IRCODE.NE.0) THEN
             WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
             RETURN
          ENDIF
c
c
c         Rescale data to represent absolute number instead of density
c
          do ir = 1,nrs
             do ik = 1,nks(ir)
                plrpad(ik,ir) = plrpad(ik,ir) * kareas(ik,ir)
             end do
          end do
c
          REF = 'ADAS H PLRP XX XXXXX ('
          WRITE(REF(13:14),'(I2)') IZ
          WRITE(REF(16:20),'(I5)') NINT(WLNGTH)
          LEN = LENSTR(REF)
          IF (ISELE.GT.0) REF = REF(1:LEN) // 'E'
          LEN = LENSTR(REF)
          IF (ISELR.GT.0) REF = REF(1:LEN) // 'R'
          LEN = LENSTR(REF)
          IF (ISELX.GT.0) REF = REF(1:LEN) // 'C'
          LEN = LENSTR(REF)
          REF = REF(1:LEN) // ') '
          LEN = LENSTR(REF)
          REF = REF(1:LEN) // XPOINT
          WRITE (IPLOT,9012) NPLOTS,REF
c
          if (cgrprint.eq.1) then
             CALL PRRMATDIV(plrpad,MAXNKS,nks(irsep),NRS,49,
     >                      REF)
          endif

C
          CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >      YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
          CALL CONTOUR (ICNTR,NGS,PLRPAD,1,1,1,FT,FP,1.0,
     >                  XOUTS,1,NXS,YOUTS,1,NYS,
     >                  XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
        ENDIF
      ENDIF
C
C-----------------------------------------------------------------------
C     POWER LOSS CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.141.OR.IREF.EQ.142) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NGS    = CNGS
        IZ     = IOPT
        REF = 'POWER LOSS ' // ZLABS(IZ)(5:11) // XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,POWLS,IZ+2,NIZS+2,MAXIZS+2,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
C
C-----------------------------------------------------------------------
C     ABSOLUTE POWER CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.143.OR.IREF.EQ.144) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NGS    = CNGS
        IZ     = IOPT
        write (6,*) '144:',absfac
        REF = 'ABS. POWER ' // ZLABS(IZ)(5:11) // XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,POWLS,IZ+2,NIZS+2,MAXIZS+2,FT,FP,ABSFAC,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
C
C-----------------------------------------------------------------------
C     ABSOLUTE H POWER CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.145.OR.IREF.EQ.146) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NGS    = CNGS
        IZ     = IOPT
        REF = 'ABS. H POWER ' // HLABS(IZ)(5:11) // XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL HCONTOUR (ICNTR,NGS,HPOWLS,IZ+1,2,2,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                conts,nconts,cntropt)
      ENDIF

C
C-----------------------------------------------------------------------
C     ABSOLUTE ACTUAL POWER - NOT DENSITY
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.147.OR.IREF.EQ.148) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NGS    = CNGS
        IZ     = IOPT
        REF = 'POWER (W)' // ZLABS(IZ)(5:11) // XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
c
c       Convert POWLS
c
        do iz = -1,nizs
           do ir = 1,nrs
              do ik = 1,nks(ir)
                 tmppowls(ik,ir,iz) = powls(ik,ir,iz) * kareas(ik,ir)
              end do
           end do
        end do
c
c       Set up contour level array.
c
        nclev = 2
        if (iz.ge.0) then
           clev(1) = pradclev(iz)
        else
           clev(1) = pradclev(0)
        endif
c
        clev(2) = 1.0
c
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
c
c       This plot is set to always be in false colour.
c
        CALL CONTOUR (1,NGS,tmpPOWLS,IZ+2,NIZS+2,MAXIZS+2,FT,FP,
     >                ABSFAC,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nclev,clev,cntropt,minscale,maxscale)
      ENDIF
C
C-----------------------------------------------------------------------
C     POWER LOSS CONTOURS - CALCULATED FROM EDGE2D IMPURITY DATA
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.149.OR.IREF.EQ.150) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NGS    = CNGS
        IZ     = IOPT
        REF = 'E2D-ABS POWER ' // ZLABS(IZ)(5:11) // XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,e2dPOWLS,IZ+1,cre2dizs+1,MAXe2dIZS+2,
     >                FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
C
C-----------------------------------------------------------------------
C     LINE RADIATION CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.151.OR.IREF.EQ.152) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NGS    = CNGS
        IZ     = IOPT
        REF    = 'LINE RADN ' // ZLABS(IZ)(5:11) // XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,LINES,IZ+2,NIZS+2,MAXIZS+2,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
C
C-----------------------------------------------------------------------
C     H LINE RADIATION CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.153.OR.IREF.EQ.154) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NGS    = CNGS
        IZ     = IOPT
        REF    = 'H LINE RADN ' // HLABS(IZ)(5:11) // XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL HCONTOUR (ICNTR,NGS,HLINES,IZ+1,2,2,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                conts,nconts,cntropt)
      ENDIF
C
C-----------------------------------------------------------------------
C     FLUID CODE LINE RADIATION CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.155.OR.IREF.EQ.156) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NGS    = CNGS
        IZ     = IOPT
        REF    = 'E2D-LINE RAD ' // ZLABS(IZ)(5:11) // XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,e2dLINES,IZ+1,cre2dizs+1,MAXe2dIZS+2,
     >                FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
C
C-----------------------------------------------------------------------
C     BREMSSTRAHLUNG CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.157.OR.IREF.EQ.158) THEN
c
C       ASSUME BREMSSTRAHLUNG IS MEASURED AT 5235 A (523.5 nm) FOR NOW
c
        WLNGTH = 5235.0
C
C       LOAD INTEGRAND
C
        CALL LDBREM(wlngth,plastmp,IRCODE,NIZS)
c
        XLAB   = '   R  (M)'
        YLAB   = '   Z  (M)'
        NGS    = CNGS
        REF    = 'BREMSSTRAHLUNG '// XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,plastmp,1,1,1,
     >                FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
c
c       Print out for now
c
        if (cgrprint.eq.1) then
c
           CALL PRRMATDIV(plastmp,MAXNKS,nks(irsep),NRS,49,
     >                    REF)
        endif
c
      ENDIF
C
C-----------------------------------------------------------------------
C     BREMSSTRAHLUNG CONTOURS - NO IMPURITY
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.159.OR.IREF.EQ.160) THEN
c
C       ASSUME BREMSSTRAHLUNG IS MEASURED AT 5235 A (523.5 nm) FOR NOW
c
        WLNGTH = 5235.0
C
C       LOAD INTEGRAND
C
        CALL LDBREM(wlngth,plastmp,IRCODE,-1)
c
        XLAB   = '   R  (M)'
        YLAB   = '   Z  (M)'
        NGS    = CNGS
        REF    = 'BREMSSTRAHLUNG-NO Z '// XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,plastmp,1,1,1,
     >                FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
c
c       Print out for now
c
        if (cgrprint.eq.1) then
c
           CALL PRRMATDIV(plastmp,MAXNKS,nks(irsep),NRS,49,
     >                    REF)
        endif
c
      ENDIF
C
C-----------------------------------------------------------------------
C     BREMSSTRAHLUNG CONTOURS FREE-FREE - NO IMPURITY
C-----------------------------------------------------------------------
C
c     WARNING: NUMBERING SYSTEM INCONSISTENCY
c
c
      IF (IREF.EQ.139.OR.IREF.EQ.140) THEN
c
C       ASSUME BREMSSTRAHLUNG IS MEASURED AT 5235 A (523.5 nm) FOR NOW
c
        WLNGTH = 5235.0
C
C       LOAD INTEGRAND
C
        CALL LDBRFF(wlngth,plastmp,IRCODE,-1)
c
        XLAB   = '   R  (M)'
        YLAB   = '   Z  (M)'
        NGS    = CNGS
        REF    = 'BREMS. FREE-FREE - NO Z '// XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,plastmp,1,1,1,
     >                FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
c
c       Print out for now
c
        if (cgrprint.eq.1) then
c
           CALL PRRMATDIV(plastmp,MAXNKS,nks(irsep),NRS,49,
     >                    REF)
        endif
c
      ENDIF
c
C-----------------------------------------------------------------------
c     Contour plots of TAUP and Total leakage contributions
c     for different regions.
C-----------------------------------------------------------------------
c
      if (iref.ge.161.and.iref.le.170) then
c
c       Initialize
c
        call rzero (nregion,maxnks*maxnrs)
        call rzero (taup,maxnks*maxnrs)
c
        if (iref.eq.161.or.iref.eq.162) then
           norm_fact = core_content
           REFINF = 'CORE     '
        elseif (iref.eq.163.or.iref.eq.164) then
           norm_fact = edge_content
           REFINF = 'EDGE     '
        elseif (iref.eq.165.or.iref.eq.166) then
           norm_fact = pp_content
           REFINF = 'TRAP     '
        elseif (iref.eq.167.or.iref.eq.168) then
           norm_fact = div_content
           REFINF = 'DIVERTOR '
        elseif (iref.eq.169.or.iref.eq.170) then
           norm_fact = main_content
           REFINF = 'MAIN SOL '
        endif
c
c       Loop over the entire grid
c
        ntotal = 0.0
c
        do ir = 1,nrs
           do ik = 1,nks(ir)
c
              if (iref.eq.161.or.iref.eq.162) then
                 nregion(ik,ir) = ncore (ik,ir)
              elseif (iref.eq.163.or.iref.eq.164) then
                 nregion(ik,ir) = nedge (ik,ir)
              elseif (iref.eq.165.or.iref.eq.166) then
                 nregion(ik,ir) = ntrap (ik,ir)
              elseif (iref.eq.167.or.iref.eq.168) then
                 nregion(ik,ir) = ndivert (ik,ir)
              elseif (iref.eq.169.or.iref.eq.170) then
                 nregion(ik,ir) = nmsol (ik,ir)
              endif
c
              ntotal = ntotal + nregion(ik,ir)
c
           end do
        end do
c
        do ir = 1,nrs
           do ik = 1,nks(ir)
c
              if (ntotal.gt.0.0) then
                 nregion(ik,ir) = nregion(ik,ir)/ ntotal * norm_fact
              endif
c
              if (tizs(ik,ir,0).ne.0.0) then
                 taup(ik,ir) = nregion(ik,ir)/
     >                         (tizs(ik,ir,0)*kareas(ik,ir))
              endif
c
           end do
        end do
c
        XLAB   = '   R  (M)'
        YLAB   = '   Z  (M)'
        NGS    = cngs
c
        REF = REFINF // 'CONTRIBUTIONS '// XPOINT
        PLANE = 'PER NEUTRAL/S/M-TOR'
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,nregion,1,1,1,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
c
        REF = REFINF // 'LEAKAGE TAUP '// XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,taup,1,1,1,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
c
      endif
C
C-----------------------------------------------------------------------
C     TOTAL DENSITIES, AVERAGED OVER K, PLOTTED VS REF, TIME DEPENDENT
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.171) THEN
        XLAB = '   Z  (M) '
        YLAB = '          '
        IZ   = IOPT
        JR   = 0
        JK   = ikref
        CALL RZERO (RVALS, MAXNRS*MAXNGS)
c
        if (zs(ikref,irsep).lt.zs(ikref,irwall)) then
           startin = 1
           endin   = irwall
           stepin  = 1
        else
           startin = irwall
           endin   = 1
           stepin  = -1
        endif
c
        DO 1666 IR = startin, endin, stepin
          JR = JR + 1
c
          jk = nks(ir)/2 +1
c
c          IF (IR.EQ.IRSEP-1) JK = JK - IKTO
c
          ROUTS(JR) = ZS(JK,IR)
          RWIDS(JR) = 0.5 * (KINDS(JK,IR) + KOUTDS(JK,IR))
          DO 1664 IT = 1, NTS
            WRITE (ELABS(IT),'(3X,1P,E8.1,A)') DWELTS(IZ)*DWELFS(IT),'S'
            VALTOT = 0.0
            DO 1662 IK = 1, NKS(IR)
              IF (IZ.EQ.NIZS+1) THEN
                VAL = 0.0
                DO 1660 JZ = 1, NIZS
                  VAL = VAL + LIMS(IK,IR,JZ,IT)
 1660           CONTINUE
              ELSEIF (IZ.EQ.-2) THEN
                VAL = FT*LIMS(IK,IR,0,IT) - FP*LIMS(IK,IR,-1,IT)
              ELSE
                VAL = LIMS(IK,IR,IZ,IT)
              ENDIF
              VALTOT = VALTOT + VAL * KAREAS(IK,IR)
 1662       CONTINUE
            RVALS(JR,IT) = VALTOT / KTOTAS(IR)
 1664     CONTINUE
 1666   CONTINUE
        REF = 'AV. DENSITY ON REF. ' // ZLABS(IZ)(5:11)
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL DRAW (ROUTS,RWIDS,RVALS,MAXNRS,IRWALL,ANLY,NTS,99,
     >    ROUTS(1),ROUTS(IRWALL),0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,2,1.0,0)
      ENDIF
C
C-----------------------------------------------------------------------
C     TOTAL DENSITIES, AVERAGED OVER K, PLOTTED VS REFERENCE LINE
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.172) THEN
        XLAB   = '   Z  (M) '
        YLAB   = '          '
        JR = 0
        JK = IKREF
c        jk = nks(irsep)/2 + 1
        CALL RZERO (RVALS, MAXNRS*MAXNGS)
c
c        write (6,*) '172a:',ikref,irsep,irwall,
c     >                   zs(ikref,irsep),zs(ikref,irwall),
c     >                   rs(ikref,irsep),rs(ikref,irwall)
c
        if (zs(ikref,irsep).lt.zs(ikref,irwall-1)) then
           startin = 2
           endin   = irwall-1
           stepin  = 1
        else
           startin = irwall-1
           endin   = 2
           stepin  = -1
        endif

c
        DO 1726 IR = startin, endin, stepin
          JR = JR + 1
c
c         This works for regular grids
c
          if (ir.lt.irsep) then 
             jk = ikins(ikref,irsep)
          else
             jk = ikref
          endif 
c
c          IF (IR.EQ.IRSEP-1) JK = jk - ikto
c
          ROUTS(JR) = ZS(JK,IR)
          RWIDS(JR) = 0.5 * (KINDS(JK,IR) + KOUTDS(JK,IR))
c
          DO 1724 IZ = -2, NIZS+1
            VALTOT = 0.0
            DO 1722 IK = 1, NKS(IR)
              IF (IZ.EQ.NIZS+1) THEN
                VAL = 0.0
                DO 1720 JZ = 1, NIZS
                  VAL = VAL + SDLIMS(IK,IR,JZ)
 1720           CONTINUE
              ELSEIF (IZ.EQ.-2) THEN
                VAL = FT*SDLIMS(IK,IR,0) - FP*SDLIMS(IK,IR,-1)
              ELSE
                VAL = SDLIMS(IK,IR,IZ)
              ENDIF
              VALTOT = VALTOT + VAL * KAREAS(IK,IR)
 1722       CONTINUE
            RVALS(JR,IZ+3) = VALTOT / KTOTAS(IR)
 1724     CONTINUE

c
          write(6,'(a,4(1x,i6),30(1x,g12.5))') '172:',ir,jk,
     >                 nks(ir)/2+1,jr,routs(jr),
     >            rwids(jr),zs(jk,ir),(rvals(jr,in),in=1,nizs+4)
c

 1726   CONTINUE
        REF = 'AVERAGE DENSITY ON REF'
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL DRAW (ROUTS,RWIDS,RVALS,MAXNRS,JR,ANLY,NIZS+4,99,
     >    ROUTS(1),ROUTS(JR),0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ZLABS,REF,NVIEW,PLANE,TABLE,2,2,1.0,0)
      ENDIF
C
C-----------------------------------------------------------------------
C     DENSITIES FOR Z < ZD, AVERAGED OVER K, PLOTTED VS REFERENCE LINE
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.173) THEN
        XLAB   = '   Z  (M) '
        YLAB   = '          '
        JR = 0
        JK = IKREF
        CALL RZERO (RVALS, MAXNRS*MAXNGS)
c
        if (zs(ikref,irsep).lt.zs(ikref,irwall)) then
           startin = 2
           endin   = irwall
           stepin  = 1
        else
           startin = irwall
           endin   = 2
           stepin  = -1
        endif
c
        DO 1736 IR = startin,endin,stepin
          JR = JR + 1
          jk = nks(ir)/2 +1
c          IF (IR.EQ.IRSEP-1) JK = JK - IKTO
          ROUTS(JR) = ZS(JK,IR)
          RWIDS(JR) = 0.5 * (KINDS(JK,IR) + KOUTDS(JK,IR))
          DO 1734 IZ = -2, NIZS+1
            VALTOT = 0.0
            KTOTA1 = 0.0
            DO 1732 IK = 1, NKS(IR)
c
c             All grids are now X-point down in OUT.
c
              IF (ZS(IK,IR).GE.CZD) then
c
c              IF ( (ZS(IK,IR).LE.CZD.and.cgridopt.eq.0
c     >                              .and.refct.eq.0).or.
c     >           (ZS(IK,IR).ge.CZD.and.
c     >           (cgridopt.eq.3.or.cgridopt.eq.2.or.cgridopt.eq.1
c     >           (cgridopt.eq.0.and.refct.eq.1))
c     >           )) THEN
c
                IF (IZ.EQ.NIZS+1) THEN
                  VAL = 0.0
                  DO 1730 JZ = 1, NIZS
                    VAL = VAL + SDLIMS(IK,IR,JZ)
 1730             CONTINUE
                ELSEIF (IZ.EQ.-2) THEN
                  VAL = FT*SDLIMS(IK,IR,0) - FP*SDLIMS(IK,IR,-1)
                ELSE
                  VAL = SDLIMS(IK,IR,IZ)
                ENDIF
                VALTOT = VALTOT + VAL * KAREAS(IK,IR)
                KTOTA1 = KTOTA1 + KAREAS(IK,IR)
              ENDIF
 1732       CONTINUE
            IF (KTOTA1.LE.0.0) THEN
               RVALS(JR,IZ+3) = 0.0
            ELSE
               RVALS(JR,IZ+3) = VALTOT / KTOTA1
            ENDIF
 1734     CONTINUE
 1736   CONTINUE
        if (cgridopt.eq.0) then
          WRITE(REF,'(''AV. DENSITY ON REF : Z<'',G11.5)') CZD
        elseif (cgridopt.eq.1.or.cgridopt.eq.2.or.cgridopt.eq.3) then
          WRITE(REF,'(''AV. DENSITY ON REF : Z>'',G11.5)') CZD
        endif
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL DRAW (ROUTS,RWIDS,RVALS,MAXNRS,JR,ANLY,NIZS+4,99,
     >    ROUTS(1),ROUTS(JR),0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ZLABS,REF,NVIEW,PLANE,TABLE,2,2,1.0,0)
      ENDIF
c
c-----------------------------------------------------------------------
c     174 - Added plot 174 - similar to 176 and 177 except it plots the 
c     ratio of emission rate to ionization rate
c-----------------------------------------------------------------------
c
c
C-----------------------------------------------------------------------
c     PLOTS OF SELECTED ADAS DATA SETS - 174 < 10eV  iopt < 50
c                                            <200eV  iopt >=50 iz=iopt-50 
C-----------------------------------------------------------------------
c
      if (iref.eq.174) then
c
c        Load up the ADAS based PEC data for the selected lines and
c        specified densities and then plot it.
c
c
         sxb =.true.
c
         if (iopt.ge.50) then 
            iz = iopt-50
         else
            IZ = IOPT
         endif
c
         read (graph(38:44),'(i4)') cz
c
         XLAB   = 'T (eV)'
c
         YLAB   = '(Emiss rate/Ionization rate)'
c        YLAB   = 'log10(Emiss rate/Ionization rate))'
c
         write (ref,'(''-'',i2,'' PEC data :'',i4)') iz,adasyr
         REF    = XFESYM(CZ) // ref
         nview   = adasid(1:36)
c
c        Setup data required for ADAS call
c
         IRCODE = 0
         CALL XXUID(ADASID)
         IF (ADASYR.GE.0) THEN
            ADASGR = 'pec??#'//XFESYM(CZ)
            WRITE(ADASGR(4:5),'(I2.2)') ADASYR
         ELSE
            ADASGR = '*'
         ENDIF
c
         ADASTY = '*'
         CALL XXSPEC(ADASGR,ADASTY,ADASEX)
C
c        For each of the requested data sets - load at 3 densities and
c        over a range of temperatures
c
         npairs = 20
c
         if (iopt.ge.50) then

            TADAS(1) = 1.0
            TADAS(2) = 2.0
            TADAS(3) = 3.0
            TADAS(4) = 4.0
            TADAS(5) = 5.0
            TADAS(6) = 7.5
            TADAS(7) = 10.0
            TADAS(8) = 20.0
            TADAS(9) = 30.0
            TADAS(10) = 40.0
            TADAS(11) = 50.0
            TADAS(12) = 60.0
            TADAS(13) = 70.0
            TADAS(14) = 80.0
            TADAS(15) = 90.0
            TADAS(16) = 100.0
            TADAS(17) = 125.0
            TADAS(18) = 150.0
            TADAS(19) = 175.0
            TADAS(20) = 200.0

         else

            TADAS(1) =   5.0
            TADAS(2) =   6.0
            TADAS(3) =   7.0
            TADAS(4) =   8.0
            TADAS(5) =   9.0
            TADAS(6) =  10.0
            TADAS(7) =  11.0
            TADAS(8) =  12.0
            TADAS(9) =  13.0
            TADAS(10) = 14.0
            TADAS(11) = 15.0
            TADAS(12) = 16.0
            TADAS(13) = 17.0
            TADAS(14) = 18.0
            TADAS(15) = 20.0
            TADAS(16) = 22.0
            TADAS(17) = 24.0
            TADAS(18) = 26.0
            TADAS(19) = 28.0
            TADAS(20) = 30.0

c            TADAS(1) = 0.5
c            TADAS(2) = 0.75
c            TADAS(3) = 1.0
c            TADAS(4) = 1.25
c            TADAS(5) = 1.50
c            TADAS(6) = 1.75
c            TADAS(7) = 2.0
c            TADAS(8) = 2.25
c            TADAS(9) = 2.5
c            TADAS(10) = 2.75
c            TADAS(11) = 3.0
c            TADAS(12) = 4.0
c            TADAS(13) = 5.0
c            TADAS(14) = 5.5
c            TADAS(15) = 6.0
c            TADAS(16) = 6.5
c            TADAS(17) = 7.0
c            TADAS(18) = 8.0
c            TADAS(19) = 9.0
c            TADAS(20) = 10.0

         endif
c
c
c        TWIDS not set to correct  values - should not be relevant
c
         do iadas = 1,npairs
            touts(iadas) = tadas(iadas)
            twids(iadas) = 1.0
         end do
c
         nplts = 0
c
c        This plot loads IONIZATION EMISSION ONLY
c
c        Do first PEC load and copy to TVALS
c
c        Load ISELE component
c
         if (isele.ne.0) then
c
            do in = 1,3

               den = 1.0e18 * 10**(in-1)
c
               pectitle = ' '
               call ldpec (den,tadas,dadas,npairs,cz,iz,isele,pecae,
     >                  ircode,wlngth,pectitle)
c
               if (ircode.eq.0) then
                  nplts = nplts + 1
                  write (6,*) '>',pectitle,'<'
                  elabs(nplts)= ' '
                  write(elabs(nplts),1761) 'E',nplts,'E',nint(wlngth),
     >                              den, isele
c
                  do iadas = 1,npairs
                     tvals(iadas,nplts) = pecae(iadas) * 1.0e-6
                     if (sxb) then
c                         
c                       Load ionization rate data
c 
                       CALL ADASRD('96',CZ,IZ+1,2,1,real(tadas(iadas)),
     >                              den,COEF)   
c
                        if (coef.ne.0.0) then
                           tvals(iadas,nplts) = tvals(iadas,nplts)/coef
                        else 
                           tvals(iadas,nplts) = 0.0
                        endif
                     endif
                  end do
c
               endif
c
            end do
c
         end if
c
c 1761    format(a1,i2,1x,1a,1x,i5,1x,'N=',1p,g8.1,1x,'D=',i4)
c
c
c        Now that all of the desired PEC's at four denisties have been
c        loaded - plot them
c
         write (6,*) '174:',nplts,iz,cz,isele,iselr,iselx,ircode
c
         do in = 1,nplts
            write (6,*) 'Plot:',in,elabs(in)
c
            pltmax = -hi
            pltmin = hi
c
            do iadas = 1,npairs
c
               if (tvals(iadas,in).gt.0.0) then
                  pltmin = min(pltmin,tvals(iadas,in))
                  pltmax = max(pltmax,tvals(iadas,in))
               endif
c
            end do

            do iadas = 1,npairs

               if (tvals(iadas,in).le.0.0) then
                  tvals(iadas,in) = pltmin
               endif

               write (6,*) touts(iadas),tvals(iadas,in)

            end do


c            do iadas = 1,npairs
c
c               if (tvals(iadas,in).gt.0.0) then
c                  pltmin = min(pltmin,log10(tvals(iadas,in)))
c                  pltmax = max(pltmax,log10(tvals(iadas,in)))
c               endif
c
c            end do
c
c            do iadas = 1,npairs
c
c               if (tvals(iadas,in).gt.0.0) then
c                  tvals(iadas,in) = log10(tvals(iadas,in))
c               else
c                  tvals(iadas,in) = pltmin
c               endif
c
c               write (6,*) touts(iadas),tvals(iadas,in),
c     >                     10.0**tvals(iadas,in)
c
c            end do
c
         end do
c
         CALL DRAW (TOUTS,TWIDS,TVALS,MAXTHE,NPAIRS,ANLY,
     >    nplts,99,0.0,touts(npairs),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,2,1.0,0)
c
      endif

c slmod begin
c
c-----------------------------------------------------------------------
c Total density for cells along a refernce IK value (1-D)
c-----------------------------------------------------------------------
c
      IF (IREF.EQ.175) THEN

        XLAB   = 'DIST (M) FROM OUTER'
        YLAB   = 'DENSITY'

        read (graph(38:44),'(i4)') cz

        if (cz.lt.0.or.cz.gt.nizs) then 
           startin = 1
           endin =nizs
        else
           startin = cz
           endin = cz
        endif

        CALL RZERO (RVALS, MAXNRS*MAXNGS)
c
        ir = irwall
        ik = iopt
        dist = 0.0
        cnt= 0
c
 1750   if (irins(ik,ir).ne.ir
     >     .and.irins(ik,ir).ne.1.and.irins(ik,ir).ne.irtrap) then
c
           jr = irins(ik,ir)
           jk = ikins(ik,ir)
c
           cnt = cnt + 1
c
           dist = dist + sqrt((zs(jk,jr)-zs(ik,ir))**2+
     >                        (rs(jk,jr)-rs(ik,ir))**2)
c
c          Set Width to 1.0 for now - only useful when averaging.
c
           routs(cnt) = dist
           rwids(cnt) = 1.0
c
           valtot = 0.0
c
           DO JZ = startin, endin
              VALTOT = VALTOT + SDLIMS(JK,JR,JZ)
           ENDDO
c
           rvals(cnt,1) = valtot
c
           write (6,*) '175:',jk,jr,cnt,routs(cnt),rvals(cnt,1)
c
           ik = jk
           ir = jr
c
c          Loop back adding points until all the way into core OR
c          to the wall in the private plasma
c
           goto 1750
c
        endif
c
        write(REF,'(''1D DENSITY FOR IK (WALL) ='',i4)') iopt
        write(NVIEW,'(''OVER RINGS: '',i4,'' TO '',i4)') irwall-1,ir
c
        write (elabs(1),'(a,i4,a,i4)') 'I   Ion Density IZ= ',startin,
     >             ' TO ', endin
c        ELABS(1) = 'I    Ion Density (Total)'
c
c        WRITE (IPLOT,9012) NPLOTS,REF
c
        CALL DRAW (ROUTS,RWIDS,RVALS,MAXNRS,CNT,ANLY,1,99,
     >    0.0,ROUTS(CNT),0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,4,2,1.0,0)
      ENDIF
c
C-----------------------------------------------------------------------
c     PLOTS OF SELECTED ADAS DATA SETS - 176 < 200eV - 177 < 10eV
C-----------------------------------------------------------------------
c
      if (iref.eq.176.or.iref.eq.177) then
c
c        Load up the ADAS based PEC data for the selected lines and
c        specified densities and then plot it.
c
c
         if (iopt.ge.50) then 
            sxb=.true.
            iz = iopt-50
         else
            sxb=.false. 
            IZ = IOPT
         endif

         read (graph(38:44),'(i4)') cz
c
         XLAB   = 'T (eV)'
c
         if (sxb) then 
            YLAB   = 'log10(Emiss rate/Ionization rate))'
         else
            YLAB   = 'log10(Emiss.[photon-m3/s])'
         endif 
c
         write (ref,'(''-'',i2,'' PEC data :'',i4)') iz,adasyr
         REF    = XFESYM(CZ) // ref
         nview   = adasid(1:36)
c
c        Setup data required for ADAS call
c
         IRCODE = 0
         CALL XXUID(ADASID)
         IF (ADASYR.GE.0) THEN
            ADASGR = 'pec??#'//XFESYM(CZ)
            WRITE(ADASGR(4:5),'(I2.2)') ADASYR
         ELSE
            ADASGR = '*'
         ENDIF
c
         ADASTY = '*'
         CALL XXSPEC(ADASGR,ADASTY,ADASEX)
C
c        For each of the requested data sets - load at 3 densities and
c        over a range of temperatures
c
         npairs = 20
c
         if (iref.eq.176) then

            TADAS(1) = 0.25
            TADAS(2) = 0.5
            TADAS(3) = 0.6
            TADAS(4) = 0.7
            TADAS(5) = 0.8
            TADAS(6) = 0.9
            TADAS(7) = 1.0
            TADAS(8) = 1.1
            TADAS(9) = 1.2
            TADAS(10) = 1.3
            TADAS(11) = 1.4
            TADAS(12) = 1.5
            TADAS(13) = 1.75
            TADAS(14) = 2.0
            TADAS(15) = 2.5
            TADAS(16) = 3.0
            TADAS(17) = 3.5
            TADAS(18) = 4.0
            TADAS(19) = 4.5
            TADAS(20) = 5.0

         elseif (iref.eq.177) then


            TADAS(1) =   5.0
            TADAS(2) =   6.0
            TADAS(3) =   7.0
            TADAS(4) =   8.0
            TADAS(5) =   9.0
            TADAS(6) =  10.0
            TADAS(7) =  11.0
            TADAS(8) =  12.0
            TADAS(9) =  13.0
            TADAS(10) = 14.0
            TADAS(11) = 15.0
            TADAS(12) = 16.0
            TADAS(13) = 17.0
            TADAS(14) = 18.0
            TADAS(15) = 20.0
c            TADAS(16) = 22.0
c            TADAS(17) = 24.0
c            TADAS(18) = 26.0
c            TADAS(19) = 28.0
c            TADAS(20) = 30.0

            TADAS(16) = 30.0
            TADAS(17) = 40.0
            TADAS(18) = 50.0
            TADAS(19) = 70.0
            TADAS(20) =100.0

c            TADAS(1) = 0.5
c            TADAS(2) = 0.75
c            TADAS(3) = 1.0
c            TADAS(4) = 1.25
c            TADAS(5) = 1.50
c            TADAS(6) = 1.75
c            TADAS(7) = 2.0
c            TADAS(8) = 2.25
c            TADAS(9) = 2.5
c            TADAS(10) = 2.75
c            TADAS(11) = 3.0
c            TADAS(12) = 4.0
c            TADAS(13) = 5.0
c            TADAS(14) = 5.5
c            TADAS(15) = 6.0
c            TADAS(16) = 6.5
c            TADAS(17) = 7.0
c            TADAS(18) = 8.0
c            TADAS(19) = 9.0
c            TADAS(20) = 10.0

         endif
c
c
c        TWIDS not set to correct  values - should not be relevant
c
         do iadas = 1,npairs
            touts(iadas) = tadas(iadas)
            twids(iadas) = 1.0
         end do
c
         nplts = 0
c
c        Do first PEC load and copy to TVALS
c
c        Load ISELE component
c
         if (isele.ne.0) then
c
            do in = 1,4

               den = 1.0e17 * 10**(in-1)
c
               pectitle = ' '
               call ldpec (den,tadas,dadas,npairs,cz,iz,isele,pecae,
     >                  ircode,wlngth,pectitle)
c
               if (ircode.eq.0) then
                  nplts = nplts + 1
                  write (6,*) '>',pectitle,'<'
                  elabs(nplts)= ' '
                  write(elabs(nplts),1761) 'E',nplts,'E',nint(wlngth),
     >                              den, isele
c
                  do iadas = 1,npairs
                     tvals(iadas,nplts) = pecae(iadas) * 1.0e-6
                     if (sxb) then
c                         
c                       Load ionization rate data
c 
                       CALL ADASRD('96',CZ,IZ+1,2,1,real(tadas(iadas)),
     >                              den,COEF)   
c
                        if (coef.ne.0.0) then
                           tvals(iadas,nplts) = tvals(iadas,nplts)/coef
                        else 
                           tvals(iadas,nplts) = 0.0
                        endif
                     endif
                  end do
c
               endif
c
            end do
c
         end if
c
c        Load ISELR pieces
c
         if (iselr.ne.0) then
c
            do in = 1,4

               den = 1.0e17 * 10**(in-1)
c
               pectitle = ' '
               call ldpec (den,tadas,dadas,npairs,cz,iz,iselr,pecar,
     >                  ircode,wlngth,pectitle)
c
               if (ircode.eq.0) then
                  nplts = nplts + 1
                  write (6,*) '>',pectitle,'<'
                  elabs(nplts)= ' '
                  write(elabs(nplts),1761) 'R',nplts,'R',nint(wlngth),
     >                               den,iselr
c
                  do iadas = 1,npairs
                     tvals(iadas,nplts) = pecar(iadas) * 1.0e-6
                     if (sxb) then
c                         
c                       Load ionization rate data
c 
                        CALL ADASRD('96',CZ,IZ+1,2,1,real(tadas(iadas)),
     >                              den,COEF)   
                        if (coef.ne.0.0) then
                           tvals(iadas,nplts) = tvals(iadas,nplts)/coef
                        else 
                           tvals(iadas,nplts) = 0.0
                        endif
                     endif 
                  end do
c
               endif
c
            end do
c
         end if
c
c        Load ISELX pieces
c
         if (iselx.ne.0) then
c
            do in = 1,4

               den = 1.0e17 * 10**(in-1)
c
               pectitle = ' '
               call ldpec (den,tadas,dadas,npairs,cz,iz,iselx,pecax,
     >                  ircode,wlngth,pectitle)
c
               if (ircode.eq.0) then
                  nplts = nplts + 1
                  write (6,*) '>',pectitle,'<'
                  elabs(nplts)= ' '
                  write(elabs(nplts),1761) 'X',nplts,'X',nint(wlngth),
     >                               den,iselx
c
                  do iadas = 1,npairs
                     tvals(iadas,nplts) = pecax(iadas) * 1.0e-6
                     if (sxb) then
c                         
c                       Load ionization rate data
c 
                        CALL ADASRD('96',CZ,IZ+1,2,1,real(tadas(iadas)),
     >                              den,COEF)   
                        if (coef.ne.0.0) then
                           tvals(iadas,nplts) = tvals(iadas,nplts)/coef
                        else 
                           tvals(iadas,nplts) = 0.0
                        endif
                     endif 
                  end do
c
               endif
c
            end do
c
         end if
c
 1761    format(a1,i2,1x,1a,1x,i5,1x,'N=',1p,g8.1,1x,'D=',i4)

c
c        Now that all of the desired PEC's at four denisties have been
c        loaded - plot them
c
         write (6,*) '176:',nplts,iz,cz,isele,iselr,iselx,ircode
c
         do in = 1,nplts
            write (6,*) 'Plot:',in,elabs(in)
c
            pltmax = -hi
            pltmin = hi
c
            do iadas = 1,npairs
c
               if (tvals(iadas,in).gt.0.0) then
                  pltmin = min(pltmin,log10(tvals(iadas,in)))
                  pltmax = max(pltmax,log10(tvals(iadas,in)))
               endif
c
            end do
c
            do iadas = 1,npairs
c
               if (tvals(iadas,in).gt.0.0) then
                  tvals(iadas,in) = log10(tvals(iadas,in))
               else
                  tvals(iadas,in) = pltmin
               endif
c
               write (6,*) touts(iadas),tvals(iadas,in),
     >                     10.0**tvals(iadas,in)
c
            end do
c
         end do
c
         CALL DRAW (TOUTS,TWIDS,TVALS,MAXTHE,NPAIRS,ANLY,
     >    nplts,99,0.0,touts(npairs),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,2,1.0,0)
c
      endif


c
C-----------------------------------------------------------------------
c     PLOTS OF SELECTED BREMSTRAHLUNG DATA - 178 < 200eV - 179 < 10eV
C-----------------------------------------------------------------------
c
      if (iref.eq.178.or.iref.eq.179) then
c
c        Load the Bremsshrahlung coefficients
c
         read (graph(38:44),'(i4)') cz
c
         XLAB   = 'T (eV)'
         YLAB   = 'log10(Emiss.[photon-m3/s/nm])'
         REF    = 'Bremstrahlung - No impurity'
c
c        For each of the requested data sets - load at 3 densities and
c        over a range of temperatures
c
         npairs = 20
c
         if (iref.eq.178) then

            TBREM(1) = 1.0
            TBREM(2) = 2.0
            TBREM(3) = 3.0
            TBREM(4) = 4.0
            TBREM(5) = 5.0
            TBREM(6) = 7.5
            TBREM(7) = 10.0
            TBREM(8) = 20.0
            TBREM(9) = 30.0
            TBREM(10) = 40.0
            TBREM(11) = 50.0
            TBREM(12) = 60.0
            TBREM(13) = 70.0
            TBREM(14) = 80.0
            TBREM(15) = 90.0
            TBREM(16) = 100.0
            TBREM(17) = 125.0
            TBREM(18) = 150.0
            TBREM(19) = 175.0
            TBREM(20) = 200.0

         elseif (iref.eq.179) then

            TBREM(1) = 1.0
            TBREM(2) = 1.25
            TBREM(3) = 1.50
            TBREM(4) = 1.75
            TBREM(5) = 2.0
            TBREM(6) = 2.25
            TBREM(7) = 2.5
            TBREM(8) = 2.75
            TBREM(9) = 3.0
            TBREM(10) = 3.5
            TBREM(11) = 4.0
            TBREM(12) = 4.5
            TBREM(13) = 5.0
            TBREM(14) = 5.5
            TBREM(15) = 6.0
            TBREM(16) = 6.5
            TBREM(17) = 7.0
            TBREM(18) = 8.0
            TBREM(19) = 9.0
            TBREM(20) = 10.0

         endif
c
c        TWIDS not set to correct  values - should not be relevant
c
         do iadas = 1,npairs
            touts(iadas) = tbrem(iadas)
            twids(iadas) = 1.0
         end do
c
         nplts = 0
c
c        Do BREM load and copy to TVALS
c
         do in = 1,4
c
            den = 1.0e17 * 10**(in-1)
c
            wlngth = 5235.0
c
            call ldbrem_spec(wlngth,npairs,den,tbrem,brempec,
     >                       ircode)
c
            if (ircode.eq.0) then
               nplts = nplts + 1
               write (6,*) '>BREM<'
               elabs(nplts)= ' '
               write(elabs(nplts),1761) 'E',nplts,'E',nint(wlngth),
     >                           den
c
               do iadas = 1,npairs
                  write (6,'(a,2i4,4g13.6)') 'BREM:',iadas,nplts,den,
     >                   tbrem(iadas),brempec(iadas),
     >                   (brempec(iadas)/den)/den
                  tvals(iadas,nplts) = (brempec(iadas)/den)/den
               end do
c
            endif
c
         end do
c
c        Now that all of the BREM's at four denisties have been
c        loaded - plot them
c
         do in = 1,nplts
            write (6,*) 'Plot:',in,elabs(in)
c
            pltmax = -hi
            pltmin = hi
c
            do iadas = 1,npairs
c
               if (tvals(iadas,in).gt.0.0) then
                  pltmin = min(pltmin,log10(tvals(iadas,in)))
                  pltmax = max(pltmax,log10(tvals(iadas,in)))
               endif
c
            end do
c
            do iadas = 1,npairs
c
               if (tvals(iadas,in).gt.0.0) then
                  tvals(iadas,in) = log10(tvals(iadas,in))
               else
                  tvals(iadas,in) = pltmin
               endif
c
               write (6,*) touts(iadas),tvals(iadas,in),
     >                     10.0**tvals(iadas,in)
c
            end do
c
         end do
c
         CALL DRAW (TOUTS,TWIDS,TVALS,MAXTHE,NPAIRS,ANLY,
     >    nplts,99,0.0,touts(npairs),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,2,1.0,0)
c
      endif
C
C-----------------------------------------------------------------------
C     TRACKS OF SELECTED IONS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.181) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NPLOTS = NPLOTS - 1
        LW     = 1
        UW     = 0
        IG     = 1
        NAME = '     TRACK'
C
 1810   CONTINUE
        DO 1812 IW = LW, MAXNWS
          IF (WALKS(IW,1).GE.10.0*RMAX) GOTO 1815
          UW = IW
 1812   CONTINUE
 1815   CONTINUE
        WRITE (REF,'(''ION TRACK'',I3,''  ('',I5,'' PTS)'')') IG,UW-LW+1
        IF (LW.LT.MAXNWS.AND.UW.GT.LW.AND.IG.LE.IOPT) THEN
          NPLOTS = NPLOTS + 1
          WRITE (IPLOT,9012) NPLOTS,REF
          CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >      YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
          CALL SUPIMP ('SELECT')
          CALL GRTRAC (WALKS(LW,1),WALKS(LW,2),UW-LW+1,NAME,'LINE',1)
          CALL FRAME
          IG = IG + 1
        ENDIF
        LW = UW + 2
        IF (LW.LT.MAXNWS) GOTO 1810
      ENDIF
c
c slmod begin
c
c ----------------------------------------------------------------------
c     Near X-point ion tracks
c ----------------------------------------------------------------------
c
      IF (IREF.EQ.182) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NPLOTS = NPLOTS - 1
        LW     = 1
        UW     = 0
        IG     = 1
        NAME = '     TRACK'
C
 1816   CONTINUE
        DO 1817 IW = LW, MAXNWS
          IF (WALKS(IW,1).GE.10.0*RMAX) GOTO 1818
          UW = IW
 1817   CONTINUE
 1818   CONTINUE
        WRITE (REF,'(''ION TRACK'',I3,''  ('',I5,'' PTS)'')') IG,UW-LW+1
        IF (LW.LT.MAXNWS.AND.UW.GT.LW.AND.IG.LE.IOPT) THEN
          NPLOTS = NPLOTS + 1
          WRITE (6,9012) NPLOTS,REF

          CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,
     +                 XXMIN,XXMAX,YYMIN,YYMAX,
     +                 TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)

c          CALL SUPIMP ('FULL')
          CALL SUPIMP ('PARIAL')
c          CALL THICK (2)
          DEFCOL = NCOLS+2
          WRITE(0,*) 'TRACKS:',lw,UW
          CALL GRTRAC (WALKS(LW,1),WALKS(LW,2),UW-LW+1,NAME,'LINE',-1)
c          CALL THICK (1)
          DEFCOL = 1
          CALL FRAME
          IG = IG + 1
        ENDIF
        LW = UW + 2
        IF (LW.LT.MAXNWS) GOTO 1816
      ENDIF
c slmod end



c
c jdemod
c
c Add HC particle track plots - most of the HC plots will be implemented
c using the generalized plotting code and adding the appropriate loading
c code in the load_divdata_array routine in outplot.o6a
c
c jdemod
c
! ammod start
C
C-----------------------------------------------------------------------
C     Hydrocarbon Particle Tracks
C-----------------------------------------------------------------------
C     183 TRACKS OF SELECTED Hydrocarbons (taken from option 181)
c     184 Magnified hydrocarbon track (taken from option 182)
C     185 PIN - TRACKS OF ALL RECORDED hydrocarbons (taken from option 651)
C
      IF (IREF.EQ.183.or.iref.eq.184.or.iref.eq.185) THEN
C
        call hc_plot_particle_track(iref,iopt)
c
      ENDIF
c
! ammod end.

C
C-----------------------------------------------------------------------
C     TOTAL RADIATED POWER CONTOURS FROM IMPURITIES + HYDROGEN
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.188.OR.IREF.EQ.189) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NGS    = CNGS
        DO IR = 1,NRS
          DO IK = 1,NKS(IR)
            PLASTMP(IK,IR) = 0.0
            DO IZ = 0,NIZS
              PLASTMP(IK,IR) = PLASTMP(IK,IR)
     >                       + ABSFAC*POWLS(IK,IR,IZ)
            ENDDO
            DO IZ = 0,1
              PLASTMP(IK,IR) = PLASTMP(IK,IR)
     >                       + HPOWLS(IK,IR,IZ)
            ENDDO
          ENDDO
        ENDDO
        REF = 'TOTAL Z + H POWER ' // XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,PLASTMP,1,1,1,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
c
C-----------------------------------------------------------------------
C
c     Plots along specific field lines (rings) of impurity
c     density and other quantities - as required
C
C-----------------------------------------------------------------------
c
c     IMPURITY DENSITY
c
      IF (IREF.EQ.191.OR.IREF.EQ.192.or.iref.eq.195.or.iref.eq.196
     >    .or.iref.eq.199) THEN
        iz = iopt
        write (elabs(1),'(i3,1x,i3,1x,''ion'')') iz,iz
        IF (IREF.EQ.191.or.iref.eq.195.or.iref.eq.199) THEN
          XLAB = '   S  (M)'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
        if (iref.eq.191.or.iref.eq.192) then
           YLAB   = 'Impurity Density'
        elseif (iref.eq.195.or.iref.eq.196) then
           YLAB   = 'LOG10(Impurity Density)'
        elseif (iref.eq.199) then
           YLAB   = 'LOG10(Impurity/deltaS)'
        endif
        READ (GRAPH(38:44),'(I4)') IR
        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)
        WRITE (IPLOT,9012) NPLOTS,REF
c
        CALL rzero (kvals, maxnks*maxngs)
        kvalmin = hi
        minkval = log10(lo)
c
        DO 1910 IK = 1, NKS(IR)
          IF (IREF.EQ.191.or.iref.eq.195.or.iref.eq.199) THEN
            KOUTS(IK) = KSS(IK,IR)
            KWIDS(IK) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
          ELSE
            KOUTS(IK) = KPS(IK,IR)
            KWIDS(IK) = 0.0
            IF (IK.GT.1) KWIDS(IK) = 0.5 * (KPS(IK,IR)-KPS(IK-1,IR))
            IF (IK.LT.NKS(IR)) KWIDS(IK) = KWIDS(IK) +
     >                               0.5 * (KPS(IK+1,IR)-KPS(IK,IR))
          ENDIF
c
          if (alphae.eq.0.0) then
           if (iz.eq.nizs+1) then
             do 1920 in = 0,nizs
                kvals(ik,1) = kvals(ik,1) + sdlims(ik,ir,in)
 1920        continue
           ELSEIF (IZ.EQ.-2) THEN
             KVALS(IK,1) = FT*SDLIMS(IK,IR,0) - FP*SDLIMS(IK,IR,-1)
           else
             KVALS(IK,1) = sdlims(ik,ir,iz)
           endif
c
c          Plot 199 plots the density as a function of the field line
c          length in the cell and not as a function of the cell
c          volume/area.
c
           if (iref.eq.199) then
              kvals(ik,1) = kvals(ik,1) * kareas(ik,ir) / kwids(ik)
           endif
c
          else
c
             call setkval(kvals(ik,1),alphae,ik,ir)
c
          endif
c
c
          if (iref.eq.195.or.iref.eq.196.or.iref.eq.199) then
             if (kvals(ik,1).le.0.0) then
c
c               Fix the lowest value possible ... if this
c               unreasonable ... it can be changed later.
c               For now ... approximate 0.0 as log10(LO).
c
                kvals(ik,1) = minkval
             else
                kvals(ik,1) = log10(kvals(ik,1))
                kvalmin = amin1(kvalmin,kvals(ik,1))
             endif
          endif
 1910   CONTINUE
c

        if (iref.eq.195.or.iref.eq.196.or.iref.eq.199) then
          do 1915 ik = 1,nks(ir)
             if (kvals(ik,1).eq.minkval) then
                kvals(ik,1) = kvalmin - 3.0
             endif
 1915     continue
        endif
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    1,99,KOUTS(1),KOUTS(NKS(IR)),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
c
c        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
c     >    1,99,KOUTS(1),mgst*KOUTS(NKS(IR)),-HI,HI,IGNORS,ITEC,AVS,NAVS,
c     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
c
        if (iref.eq.195.or.iref.eq.196.or.iref.eq.199) then
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    1,99,KOUTS(1),0.02*KOUTS(NKS(IR)),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >   1,99,KOUTS(1),0.005*KOUTS(NKS(IR)),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
c
c        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
c     >    1,99,mgnd*KOUTS(nks(ir)),KOUTS(NKS(IR)),-HI,HI,IGNORS,
c     >    ITEC,AVS,NAVS,
c     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    1,99,0.98*KOUTS(nks(ir)),KOUTS(NKS(IR)),-HI,HI,IGNORS,
     >    ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    1,99,0.995*KOUTS(nks(ir)),KOUTS(NKS(IR)),-HI,HI,IGNORS,
     >    ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
       endif
      ENDIF
c
C-----------------------------------------------------------------------
c
c     Impurity radiation - POWLS
c
C-----------------------------------------------------------------------
c
      IF (IREF.EQ.193.OR.IREF.EQ.194.or.iref.eq.197.or.iref.eq.198) THEN
        iz = iopt
        write (elabs(1),'(i3,1x,i3,1x,''ion'')') iz,iz
        IF (IREF.EQ.193.or.iref.eq.197) THEN
          XLAB = '   S  (M)'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
        if (iref.eq.193.or.iref.eq.194) then
           YLAB   = 'Impurity Power Loss'
        elseif (iref.eq.197.or.iref.eq.198) then
           YLAB   = 'LOG10(Imp Power Loss)'
        endif
        READ (GRAPH(38:44),'(I4)') IR
        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)
        WRITE (IPLOT,9012) NPLOTS,REF
c
        CALL rzero (kvals, maxnks*maxngs)
        kvalmin = hi
        minkval = log10(lo)
c
        DO 1930 IK = 1, NKS(IR)
          IF (IREF.EQ.193.or.iref.eq.197) THEN
            KOUTS(IK) = KSS(IK,IR)
            KWIDS(IK) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
          ELSE
            KOUTS(IK) = KPS(IK,IR)
            KWIDS(IK) = 0.0
            IF (IK.GT.1) KWIDS(IK) = 0.5 * (KPS(IK,IR)-KPS(IK-1,IR))
            IF (IK.LT.NKS(IR)) KWIDS(IK) = KWIDS(IK) +
     >                               0.5 * (KPS(IK+1,IR)-KPS(IK,IR))
          ENDIF
          if (iz.eq.nizs+1) then
             do 1940 in = 0,nizs
                kvals(ik,1) = kvals(ik,1) + powls(ik,ir,in)
 1940        continue
          ELSEIF (IZ.EQ.-2) THEN
             KVALS(IK,1) = FT*powls(IK,IR,0) - FP*powls(IK,IR,-1)
          else
             KVALS(IK,1) = powls(ik,ir,iz)
          endif
          if (iref.eq.197.or.iref.eq.198) then
             if (kvals(ik,1).le.0.0) then
c
c               Fix the lowest value possible ... if this
c               unreasonable ... it can be changed later.
c               For now ... approximate 0.0 as log10(LO).
c
                kvals(ik,1) = minkval
             else
                kvals(ik,1) = log10(kvals(ik,1))
                kvalmin = amin1(kvalmin,kvals(ik,1))
             endif
           endif
 1930   CONTINUE
c
        do 1935 ik = 1,nks(ir)
           if (kvals(ik,1).eq.minkval) then
              kvals(ik,1) = kvalmin - 3.0
           endif
 1935   continue
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    1,99,KOUTS(1),KOUTS(NKS(IR)),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)

        if (clsup.eq.1)
     >   CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    1,99,KOUTS(1),mgst*KOUTS(NKS(IR)),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
      ENDIF


      return

c
c     Fortmat statements
c

 9012 FORMAT(1X,'PLOT',I3,4X,A)

      end 

