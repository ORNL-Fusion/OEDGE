      subroutine out500(iref,graph,iopt,ierr)
      use mod_params
      use mod_outcom
      use mod_cgeom
      use mod_comtor
      use mod_dynam2
      use mod_dynam3
      use mod_reiser_com
      implicit none
      integer iref,iopt,ierr
      character*(*) graph

c
c     include 'params'
c     include 'outcom'
c
c     Other common blocks
c
c     include 'cgeom'
c     include 'comtor'
c      include 'cneut2'
c     include 'dynam2'
c     include 'dynam3'
c      include 'dynam4'
c      include 'pindata'
c      include 'cadas'
c      include 'grbound'
c      include 'outxy'
c      include 'cedge2d'
c      include 'transcoef'
c      include 'cioniz'
c     include 'reiser_com' 
c      include 'printopt' 
c
c     Local Variables
c
c
      integer ik,ir,iz
      integer in
c
c     Local Variables
c




      integer inc





c
c psmod
c
c     Variables used for average force plots 
c
      REAL KKpar,DDparpar,Kaa,Kab,Kaba,Kabb,Kac,Daa,Dab,Dac
      real LAMBDA2
      REAL TB,NB,TBG,VBG,KBetaf,BISCT,Vimp,DDIV,TAUpara
      INTEGER ZB,ZI,MI,MB,XPER,kounter,jcounter,negvalue,posvalue
      INTEGER FORCE,FILED
      REAL XVALUE(100),WIDTHS(100),YVALUE(100,MAXNGS),ROOTS(5)
      REAL Fsum(MAXNKS,MAXNRS,MAXIZS)
c
      REAL TARGETVB,Frictionf
      REAL CHIpara,CHIperp,CHI
      integer icounter,aa
      REAL Phi00,Phi01,Phi10,Phi11,Phi20,Phi21,Psi00,Psi01
      REAL Psi02,Psi10,Psi11,Psi12,Psi20,Psi21,Psi22
      REAL Cplus,Cminus

c
c psmod
c



      IF (IOPT.EQ.0) return



c
c psmod
c
      KKpar    = 0.0
      KBetaf   = 0.0
      DDparpar = 0.0
      CHIpara  = 0.0
      CHIperp  = 0.0
      Kaa      = 0.0
      Kab      = 0.0
      Kaba     = 0.0
      Kabb     = 0.0
      Kac      = 0.0
      Daa      = 0.0
      Dab      = 0.0
      Dac      = 0.0
      LAMBDA2  = 0.0
      VBG      = 0.0
      TBG      = 0.0
      TB       = 0.0
      NB       = 0.0
      DDIV     = 0.0
      TAUpara  = 0.0
      FORCE    = 0
      ZI       = 0
      ZB       = 0
      MB       = 0
      MI       = 0
      Xper     = 0
      FILED    = 0
      Frictionf= 0.0
c
c psmod
c

c
c psmod
c
c    DATA FOR REISER/DIVIMP COMPARISON OF FF+FTH VS CHI
c
      IF (IREF.EQ.579.OR.IREF.EQ.580.OR.IREF.EQ.581.OR.IREF.EQ.584)
     > THEN
        CALL RDG579(GRAPH,ZB,ZI,MB,MI,NB,TBG,TB,VBG,XPER,FILED,IERR)
        IF (IERR.NE.0) THEN
          WRITE(6,*) 'ERROR READING REISER/DIVIMP DETAILS, IREF = ',
     >                IREF
          IERR = 0
          return
        ENDIF
      ENDIF
c
c psmod
c



      call init_plot(iref,graph,iopt)



C
C-----------------------------------------------------------------------
C     CONTOURS OF ELECTRON TEMPERATURE
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.501.OR.IREF.EQ.502) THEN
        XLAB   = '   R  (M)'
        YLAB   = '   Z  (M)'
        NGS    = CNGS
        REF    = 'ELECTRON TEMPERATURE ' // XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,KTEBS,1,1,1,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
C
C-----------------------------------------------------------------------
C     CONTOURS OF ION TEMPERATURE
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.503.OR.IREF.EQ.504) THEN
        XLAB   = '   R  (M)'
        YLAB   = '   Z  (M)'
        NGS    = CNGS
        REF    = 'ION TEMPERATURE ' // XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,KTIBS,1,1,1,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
C
C-----------------------------------------------------------------------
C     CONTOURS OF ELECTRON DENSITY
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.505.OR.IREF.EQ.506) THEN
        XLAB   = '   R  (M)'
        YLAB   = '   Z  (M)'
        NGS    = CNGS
        REF    = 'ELECTRON DENSITY ' // XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,KNBS,1,1,1,FT,FP,RIZB,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
C
C-----------------------------------------------------------------------
C     CONTOURS OF PLASMA FLOW VELOCITY
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.507.OR.IREF.EQ.508) THEN
        XLAB   = '   R  (M)'
        YLAB   = '   Z  (M)'
        NGS    = CNGS
        REF    = 'PLASMA FLOW VELOCITY ' // XPOINT
c
c       Rescale VElocity by dividing by qtim
c
c        do ir = 1,nrs
c           do ik = 1,nks(ir)
c              cvalsa(ik,ir) = kvhs(ik,ir)/qtim
c           end do
c        end do

c       compute mach number

        do ir = 1,nrs
           do ik = 1,nks(ir)
              cvalsa(ik,ir) = kvhs(ik,ir)/9.79e3/
     >                        sqrt((ktebs(ik,ir)+ktibs(ik,ir))/crmb)
           end do
        end do

        WRITE (IPLOT,9012) NPLOTS,REF
c
c        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
c     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
c        CALL CONTOUR (ICNTR,NGS,cvalsa,1,1,1,FT,FP,1.0,
c     >                XOUTS,1,NXS,YOUTS,1,NYS,
c     >                XXMIN,XXMAX,YYMIN,YYMAX,
c     >                nconts,conts,cntropt,minscale,maxscale)
c
c       normalization added by Krieger IPP/97
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
c
c       if filled plot selected, use special color table for
c       directed quantities  Krieger, IPP/97
c
        if (icntr.eq.0) then
          CALL CONTOUR (ICNTR,NGS,KVHS,1,1,1,FT,FP,1.0/qtim,
     >                  XOUTS,1,NXS,YOUTS,1,NYS,
     >                  XXMIN,XXMAX,YYMIN,YYMAX,
     >                  nconts,conts,cntropt,minscale,maxscale)
        else
c
c         Setup 40 colour plots
c         
          call setup_col(23,4)

c         CALL CONTOUR (ICNTR,NGS,KVHS,1,1,1,FT,FP,1.0/qtim,
          CALL CONTOUR (ICNTR,NGS,cvalsa,1,1,1,FT,FP,1.0/qtim,
     >                  XOUTS,1,NXS,YOUTS,1,NYS,
     >                  XXMIN,XXMAX,YYMIN,YYMAX,
     >                  nconts,conts,5,minscale,maxscale)

c
c         restore colours to original
c
          call setup_col(n_cols,col_opt)

        endif
c
      ENDIF
C
C-----------------------------------------------------------------------
C     CONTOURS OF ELECTRIC FIELD
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.509.OR.IREF.EQ.510) THEN
        XLAB   = '   R  (M)'
        YLAB   = '   Z  (M)'
        NGS    = CNGS
        REF    = 'ELECTRIC FIELD ' // XPOINT
c
c       Rescale Electirc field by dividing by qtim
c
C        FACT  = QTIM * QTIM * EMI / CRMI
c
C        do ir = 1,nrs
C           do ik = 1,nks(ir)
C              cvalsa(ik,ir) = kes(ik,ir) / fact
C           end do
C        end do
c
        WRITE (IPLOT,9012) NPLOTS,REF
C
c       normalization added by Krieger IPP/97
C
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
C
c       if filled plot selected, use special color table for
c       directed quantities  Krieger, IPP/97
C
        if (icntr.eq.0) then
          CALL CONTOUR (ICNTR,NGS,KES,1,1,1,FT,FP,
     >                  1.0/(qtim*qtim*emi/crmi),
     >                  XOUTS,1,NXS,YOUTS,1,NYS,
     >                  XXMIN,XXMAX,YYMIN,YYMAX,
     >                  nconts,conts,cntropt,minscale,maxscale)
        else
c
c         Setup 40 colour plots
c         
          call setup_col(41,4)
c
          CALL CONTOUR (ICNTR,NGS,KES,1,1,1,FT,FP,
     >                  1.0/(qtim*qtim*emi/crmi),
     >                  XOUTS,1,NXS,YOUTS,1,NYS,
     >                  XXMIN,XXMAX,YYMIN,YYMAX,
     >                  nconts,conts,5,minscale,maxscale)

c
c         restore colours to original
c
          call setup_col(n_cols,col_opt)
c
        endif
c
      ENDIF
c
c kkmod
c
C
C-----------------------------------------------------------------------
C     CONTOURS OF MEAN CHARGE STATE  added Krieger IPP/98
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.511.OR.IREF.EQ.512) THEN
        XLAB   = '   R  (M)'
        YLAB   = '   Z  (M)'
c
        REF    = 'Zmean ' // XPOINT
c
        write (6,*) 'range:',xxmin,xxmax,yymin,yymax,zminp,zmaxp
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
c
c       Calculate imp. density for requested charge state.
c
        call rzero(plastmp,maxnks*maxnrs)
c
c       we want here the "real" density, not per injected
c       particle -> multiply by absfac
c
        do iz = 0,nizs
          do ir = 1,nrs
            do ik = 1,nks(ir)
              plastmp(ik,ir) = plastmp(ik,ir) + sdlims(ik,ir,iz)
            end do 
          end do
        end do 
        do ir = 1,nrs
          do ik = 1,nks(ir)
            plastmp(ik,ir) = zeffs(ik,ir,1)/(plastmp(ik,ir)*absfac)
          end do
        end do 
c
        CALL CONTOUR (ICNTR,NGS,PLASTMP,1,1,1,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,1,minscale,maxscale)
      ENDIF
c
c kkmod 





















c
c psmod
c
C
C-----------------------------------------------------------------------
C  CONTOURS OF FRICTIONAL FORCE 
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.551.OR.IREF.EQ.552) THEN
        XLAB   = '   R  (M)'
        YLAB   = '   Z  (M)'
        NGS    = CNGS
        READ (GRAPH(39:43),'(I5)') IZ
        IF(IREF.EQ.551) THEN
          WRITE (REF,'(''FRICTIONAL FORCE FOR IZ='',I3)')IZ
        ELSE
          WRITE (REF,'(''FRICTION FORCE NEAR X PT: IZ='',I3)')IZ
        ENDIF 
c
        WRITE (IPLOT,9012) NPLOTS,REF
c
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
c
C		IPP/09 Krieger - 40 Color plot option added 
C     
        if (icntr.eq.0) then
        CALL CONTOUR (ICNTR,NGS,Ffi,IZ,IZ,MAXIZS,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
        else
c
c         Setup 40 colour plots
         WRITE(6,*) 'Out551 using 40 color plots'
          call setup_col(41,4)
c		  		    
          CALL CONTOUR (ICNTR,NGS,Ffi,IZ,IZ,MAXIZS,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,5,minscale,maxscale) 
c
c         restore colours to original
c
          call setup_col(n_cols,col_opt)
     
        endif   
      ENDIF
C
C-----------------------------------------------------------------------
C  CONTOURS OF THERMAL FORCE (VIA REISER'S DRIFT-KINETIC FORMULATION OR
C  VIA THE FLUID APPROX. OF NEUHAUSER AS SELECTED IN THE INPUT DATA FILE)
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.553.OR.IREF.EQ.554) THEN
        XLAB   = '   R  (M)'
        YLAB   = '   Z  (M)'
        NGS    = CNGS
        READ (GRAPH(39:43),'(I5)') IZ
        IF(IREF.EQ.553) THEN
          WRITE (REF,'(''THERMALAL FORCE FOR IZ='',I3)')IZ
        ELSE
          WRITE (REF,'(''THERMAL FORCE NEAR X PT: IZ='',I3)')IZ
        ENDIF 
c
        WRITE (IPLOT,9012) NPLOTS,REF
c
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
c
C		IPP/09 - Krieger 40 Color plot option added 
C
        if (icntr.eq.0) then
          CALL CONTOUR (ICNTR,NGS,Fthi,IZ,IZ,MAXIZS,FT,FP,1.0,
     >                  XOUTS,1,NXS,YOUTS,1,NYS,
     >                  XXMIN,XXMAX,YYMIN,YYMAX,
     >                  nconts,conts,cntropt,minscale,maxscale)
        else
c
c         Setup 40 colour plots
          WRITE(6,*) 'Out551 using 40 color plots'
          call setup_col(41,4)
c		  		             
          CALL CONTOUR (ICNTR,NGS,Fthi,IZ,IZ,MAXIZS,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,5,minscale,maxscale)
c
c         restore colours to original
c
          call setup_col(n_cols,col_opt)
        endif 
      ENDIF
C
C-----------------------------------------------------------------------
C     CONTOURS OF VELOCITY-GRADIENT FORCE
C     (VIA REISER'S DRIFT-KINETIC FORMULATION)
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.555.OR.IREF.EQ.556) THEN
        XLAB   = '   R  (M)'
        YLAB   = '   Z  (M)'
        NGS    = CNGS
        READ (GRAPH(39:43),'(I5)') IZ
        IF(IREF.EQ.555) THEN
          WRITE (REF,'(''Vb-GRAD FORCE FOR IZ='',I3)')IZ
        ELSE
          WRITE (REF,'(''Vb-GRAD FORCE NEAR X PT: IZ='',I3)')IZ
        ENDIF 
c
        WRITE (IPLOT,9012) NPLOTS,REF
c
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
c
C		IPP/09 - Krieger 40 Color plot option added 
C
        if (icntr.eq.0) then
          CALL CONTOUR (ICNTR,NGS,Fvbg,IZ,IZ,MAXIZS,FT,FP,1.0,
     >                  XOUTS,1,NXS,YOUTS,1,NYS,
     >                  XXMIN,XXMAX,YYMIN,YYMAX,
     >                  nconts,conts,cntropt,minscale,maxscale)
        else
c
c         Setup 40 colour plots
         WRITE(6,*) 'Out551 using 40 color plots'
          call setup_col(41,4)
C
          CALL CONTOUR (ICNTR,NGS,Fvbg,IZ,IZ,MAXIZS,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,5,minscale,maxscale)
c
c         restore colours to original
c
          call setup_col(n_cols,col_opt)
        endif      
c
      ENDIF
C
C-----------------------------------------------------------------------
C     CONTOURS OF COMBINED FRICTIONAL AND THERMAL FORCE
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.557.OR.IREF.EQ.558) THEN
        XLAB   = '   R  (M)'
        YLAB   = '   Z  (M)'
        NGS    = CNGS
        READ (GRAPH(39:43),'(I5)') IZ
        IF(IREF.EQ.557) THEN
          WRITE (REF,'(''FRICT AND THERM FORCE FOR IZ='',I3)')IZ
        ELSE
          WRITE (REF,'(''FRICT + THERM FORCE NEAR X PT:IZ='',
     >                   I3)')IZ
        ENDIF 
c
        WRITE (IPLOT,9012) NPLOTS,REF
c
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
c
        CALL RZERO(Fsum,MAXNKS*MAXNRS*MAXIZS)
c
        DO ir = 1,nrs
          DO ik = 1,nks(ir)
            Fsum(ik,ir,iz) = Ffi(ik,ir,IZ)+ Fthi(ik,ir,IZ) 
          END DO
        END DO
c
c        WRITE(0,*)'Fsum',Fsum(25,8,iz),' ',ffi(25,8,iz),' ',
c     >                                      fthi(25,8,iz)
c
C		IPP/09 - Krieger 40 Color plot option added 
C
        if (icntr.eq.0) then
          CALL CONTOUR (ICNTR,NGS,Fsum,IZ,IZ,MAXIZS,FT,FP,1.0,
     >                  XOUTS,1,NXS,YOUTS,1,NYS,
     >                  XXMIN,XXMAX,YYMIN,YYMAX,
     >                  nconts,conts,cntropt,minscale,maxscale)
        else
c
c         Setup 40 colour plots
         WRITE(6,*) 'Out557 using 40 color plots'
          call setup_col(41,4)
C
          CALL CONTOUR (ICNTR,NGS,Fsum,IZ,IZ,MAXIZS,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,5,minscale,maxscale)
c
c         restore colours to original
c
          call setup_col(n_cols,col_opt)
        endif      

      ENDIF
C
C-----------------------------------------------------------------------
C     CONTOURS OF THE SUM OF ALL THE FORCES ACTING UPON 
C     AN IMPURITY ION AS SELECTED IN THE INPUT DATA FILE
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.559.OR.IREF.EQ.560) THEN
        XLAB   = '   R  (M)'
        YLAB   = '   Z  (M)'
        NGS    = CNGS
        READ (GRAPH(39:43),'(I5)') IZ
        IF(IREF.EQ.559) THEN
          WRITE (REF,'(''TOTAL FORCES FOR IZ='',I3)')IZ
        ELSE
          WRITE (REF,'(''TOTAL FORCES NEAR X PT: IZ='',I3)')IZ
        ENDIF 
c
        WRITE (IPLOT,9012) NPLOTS,REF
c
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
c
c        WRITE(0,*)'Fcell',Fcell(25,8,iz),' ',ffi(25,8,iz),' ',
c     >                                      fthi(25,8,iz)
c
C		IPP/09 - Krieger 40 Color plot option added 
C
        if (icntr.eq.0) then
          CALL CONTOUR (ICNTR,NGS,Fcell,IZ,IZ,MAXIZS,FT,FP,1.0,
     >                  XOUTS,1,NXS,YOUTS,1,NYS,
     >                  XXMIN,XXMAX,YYMIN,YYMAX,
     >                  nconts,conts,cntropt,minscale,maxscale)
        else
c
c         Setup 40 colour plots
         WRITE(6,*) 'Out557 using 40 color plots'
          call setup_col(41,4)
C
          CALL CONTOUR (ICNTR,NGS,Fcell,IZ,IZ,MAXIZS,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,5,minscale,maxscale)
c
c         restore colours to original
c
          call setup_col(n_cols,col_opt)

        endif      
      ENDIF
C
C-----------------------------------------------------------------------
C     CONTOURS OF THE RESULTANT VELOCITY DISTRIBUTION OF 
C     AN IMPURITY ION AS SELECTED IN THE INPUT DATA FILE
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.561.OR.IREF.EQ.562) THEN
        XLAB   = '   R  (M)'
        YLAB   = '   Z  (M)'
        NGS    = CNGS
        READ (GRAPH(39:43),'(I5)') IZ
        IF(IREF.EQ.561) THEN
          WRITE (REF,'(''AVE. VELOCITY FOR IZ='',I3)')IZ
        ELSE
          WRITE (REF,'(''AVE. VELOCITY NEAR X PT: IZ='',I3)')IZ
        ENDIF 
c
        WRITE (IPLOT,9012) NPLOTS,REF
c
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
c
c        WRITE(0,*)'Fcell',Fcell(25,8,iz),' ',ffi(25,8,iz),' ',
c     >                                      fthi(25,8,iz)
c
        CALL CONTOUR (ICNTR,NGS,VELavg,IZ,IZ,MAXIZS,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)

      ENDIF
C
C-----------------------------------------------------------------------
C     AVERAGED FORCES (Frictional,Thermal,Vbgrad) ALONG THE FIELD LINE
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.575) THEN
c
        ELABS(1) = 'FF  FF     '
        ELABS(2) = 'FIG FIG    '
        ELABS(3) = 'FvbgFvbg   '
        ELABS(4) = 'DIFFDIFF   '
*        ELABS(4) = 'FnetFnet   '
        ELABS(5) = 'FNETFNET(D)'

        XLAB = '   S  (M)'
        YLAB = ' FORCE(N)  '

        NVIEW = 'AVERAGED FORCES ALONG B'        

*        FACT = QTIM**2 * EMI / CRMI
c
c Add ionisation state info here
c
        READ (GRAPH(28:31),'(I4)') IZ
        READ (GRAPH(36:39),'(I4)') IR
        READ (GRAPH(43:44),'(I2)') FILED
c
c        IF (FILED.EQ.1) THEN
c          WRITE(121,*)'575-REISER FORCES'
c        ENDIF
c 
        NPLOTS = NPLOTS + 1
c
        WRITE (REF,'(''Ring: '',I3,'', K: '',F7.4,'' ,State: '',I3)') 
     +        IR,KKS(IR),IZ
c
        WRITE (IPLOT,9012) NPLOTS,REF
c
        write(iplot,'(2x,''IK'',2x,''IR'',5x,''FF'',5x,5x'//
     >                   ',''FIG'',4x,5x,''FEG'',4x,5x'//
     >                    ',''FE'',5x,5x,''NET'',4x,6x,''S'',5x'//
     >                    ',6x,''R'',5x,6x,''Z'')')
c
c
        CALL rzero (kvals, maxnks*maxngs)
c
c       Modify the following to include values for S=0 ... even if this
c       is not part of the grid. Assume that if one end is not a grid
c       point then neither is the other.
c
c
        if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
           in  = 0
           inc = 0
        else
           in       = 1
           inc      = 2
           kouts(1) = 0.0
c
c Fix constant for density units
c

           KWIDS(1)         = kss(1,ir)
           kouts(nks(ir)+2) = ksmaxs(ir)
           kwids(nks(ir)+2) = ksmaxs(ir) - kss(nks(ir),ir)

c For S=0 meters (at target plate), the Y-axis values will be taken 
c as those calculated for the first grid cell: 
                
          KVALS(1,1) = Ffi(1,IR,IZ)
          KVALS(1,2) = Fthi(1,IR,IZ)
          KVALS(1,3) = Fvbg(1,IR,IZ)
          KVALS(1,4) = DIFF(1,IR,IZ) 
*          KVALS(1,4) = KVALS(1,1) + KVALS(1,2) + KVALS(1,3)
          KVALS(1,5) = Fcell(1,IR,IZ)
           

c For S=Smax (at the other target plate), the Y-axis values will be
c taken as those calculated for the last grid cell: 

          KVALS(nks(ir)+2,1) = Ffi(nks(ir),ir,iz)
          KVALS(nks(ir)+2,2) = Fthi(nks(ir),ir,iz)
          KVALS(nks(ir)+2,3) = Fvbg(nks(ir),ir,iz)
          KVALS(nks(ir)+2,4) = DIFF(nks(ir),ir,iz)
*          KVALS(nks(ir)+2,4) = KVALS(nks(ir)+2,1)+KVALS(nks(ir)+2,2)
*     +                       + KVALS(nks(ir)+2,3)
          KVALS(nks(ir)+2,5) = Fcell(nks(ir),ir,iz)  
         
c
        endif
c
*        DO IK = 1, NKS(IR)
*           WRITE(0,*) Fthi(IK,IR,IZ)
*        END DO
c
c    !!!
c       Due to the fact that arrays Ffi, Fthi, Fvbg, DIFF, and Fcell are 
c       initially evaluated in units of distance and then converted back 
c       into units of force for the plotting routine, the values of the 
c       arrays become prone to round of error because they are REAL valued 
c       instead of DOUBLE PRECISION.  As a result, when many impurity 
c       particles are released (1000) with a small time step, QTIM, (5e-8 s)
c       the average value of the array for those grid cells where many 
c       iterations occur will be substantially different than the true value 
c       and will produce a misleading profile.  Note that the actual ion 
c       motions are correctly calculated, it is only the generated profiles
c       that are affected.  This is most evident when using the FIG of the 
c       fluid approximation model, since its value is always the same
c       for any grid cell, yet the force profile will be different for 
c       different QTIM's and total impurity ion numbers.
c    !!!  

        DO 1270 IK = 1, NKS(IR)

          KOUTS(IK+in) = KSS(IK,IR)
          KWIDS(IK+in) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))

          KVALS(IK+in,1) = Ffi(IK,IR,IZ)
          KVALS(IK+in,2) = Fthi(IK,IR,IZ)
          KVALS(IK+in,3) = Fvbg(IK,IR,IZ)
          KVALS(IK+in,4) = DIFF(IK,IR,IZ)
*          KVALS(IK+in,4) = KVALS(IK+IN,1) + KVALS(IK+IN,2) 
*     +                   + KVALS(IK+IN,3)
          KVALS(IK+in,5) = Fcell(IK,IR,IZ)
c
          WRITE(iplot,'(2I4,8G12.4)') IK,IR,KVALS(IK+IN,1),
     +      KVALS(IK+IN,2),
     +      KVALS(IK+IN,3),KVALS(IK+IN,4),KVALS(IK+IN,5),KSS(IK,IR),
     +      RS(IK,IR),ZS(IK,IR)
c
c          IF(FILED.EQ.1)THEN
c            WRITE(121,8)KOUTS(IK+in),KVALS(IK+in,1),
c     >                   KVALS(IK+in,2),KVALS(IK+in,3),
c     >                   KVALS(IK+in,4),KVALS(IK+in,5)
c    8       FORMAT(6(1X,G11.5))
c          ENDIF
c
 1270   CONTINUE

C
        enldist = mgst*kouts(nks(ir)+inc)
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >    5,99,KOUTS(1),KOUTS(NKS(IR)+inc),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0)
c
        IF (clsup.eq.1) THEN
c
          CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >      5,99,KOUTS(1),enldist   ,-HI,HI,IGNORS,ITEC,AVS,
     >      NAVS,
     >      JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0)

          CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >      5,99,KOUTS(NKS(IR)+inc)-enldist,KOUTS(NKS(IR)+inc),
     >      -HI,HI,IGNORS,ITEC,AVS,NAVS,
     >      JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0)
        ENDIF
c
      ENDIF
C
c-----------------------------------------------------------------------
c  Force Balance Plots for Reiser Formulation when Impurity Velocity = 0
c-----------------------------------------------------------------------
c
      IF (IREF.EQ.576) THEN

        ELABS(1) = 'FF  FF     '
        ELABS(2) = 'FIG FIG    '
        ELABS(3) = 'FvbgFvbg   '
        ELABS(4) = 'FnetFnet   '

        XLAB = '   S  (M)'
        YLAB = ' FORCE(N)  '

        NVIEW = 'REISER: FORCE BALANCE WHEN VEL = 0'
        
        READ (GRAPH(32:35),'(I4)') IZ
        READ (GRAPH(40:44),'(I4)') IR
        READ (GRAPH(48:49),'(I2)') FILED
c
c        IF (FILED.EQ.1) THEN
c          WRITE(121,*)'576-FORCE BALANCE, VZ = 0'
c        ENDIF
c
        NPLOTS = NPLOTS + 1

        WRITE (REF,'(''Ring: '',I3,'', K: '',F7.4,'' ,State: '',I3)') 
     +        IR,KKS(IR),IZ

        WRITE (IPLOT,9012) NPLOTS,REF

        write(iplot,'(2x,''IK'',2x,''IR'',5x,''FF'',5x,5x'//
     >                   ',''FIG'',4x,5x,''FEG'',4x,5x'//
     >                    ',''FE'',5x,5x,''NET'',4x,6x,''S'',5x'//
     >                    ',6x,''R'',5x,6x,''Z'')')

        CALL rzero (kvals, maxnks*maxngs)

c       Modify the following to include values for S=0 ... even if this
c       is not part of the grid. Assume that if one end is not a grid
c       point then neither is the other.

        if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
           in  = 0
           inc = 0
        else
           in       = 1
           inc      = 2
           kouts(1) = 0.0

           KWIDS(1)         = kss(1,ir)
           kouts(nks(ir)+2) = ksmaxs(ir)
           kwids(nks(ir)+2) = ksmaxs(ir) - kss(nks(ir),ir)
        endif

        CALL Vbgrad
        CALL Coeff(NIZS)


c For s = 0 meters (at outer target plate)

      FACT = QTIM * QTIM * EMI / CRMI

      TARGETVB = -SQRT((KTEDS(IDDS(IR,2))+KTIDS(IDDS(IR,2)))
     >           *1.602E-19/(CRMB*1.67E-27))*QTIM

      KVHGS(1,IR) = (KVHS(1,IR)-TARGETVB)/KSS(1,IR)

      LAMBDA1(IZ) = 0.23995*Coulomb_log*REAL(IZ*IZ)*
     >              REAL(RIZB*RIZB)/(CRMI*CRMI)
 
      ALPHAI(1,IR) = SQRT(CRMB/(2*KTIDS(IDDS(IR,2))*EMI))
          
      CC1(1,IR) = 1.8024009E17*KTIDS(IDDS(IR,2))*KFIDS(IDDS(IR,2))/
     >   (Coulomb_log*REAL(RIZB*RIZB*RIZB*RIZB)*KNDS(IDDS(IR,2))*FACT)

      CC2(1,IR)=-7.361216E12*KVHGS(1,IR)*SQRT(CRMB*
     >      KTIDS(IDDS(IR,2))*KTIDS(IDDS(IR,2))*KTIDS(IDDS(IR,2)))/
     >      (Coulomb_log*RIZB*RIZB*RIZB*(1.2*RIZB + 0.8485)*
     >      KNDS(IDDS(IR,2))*QTIM) 

      CHIpara = ALPHAI(1,IR)*(-KVDS(IDDS(IR,2)))
      LAMBDA2 = LAMBDA1(IZ)*KNDS(IDDS(IR,2))

      CALL Coulomb_Coll(KKpar,DDparpar,CHIpara,LAMBDA2,
     >                    1,IR,Kaa,Kab,Kac,Daa,Dab,Dac)

      KVALS(1,1) = Kaa * CRMI * amu 
      KVALS(1,2) = Kab * CRMI * amu
      KVALS(1,3) = Kac * CRMI * amu
      KVALS(1,4) = KVALS(1,1) + KVALS(1,2) + KVALS(1,3) 

*      WRITE(0,*)'576:'
*      WRITE(0,*)KVALS(1,1)
*      WRITE(0,*)KVALS(1,2)
*      WRITE(0,*)KVALS(1,3)
*      WRITE(0,*)KVALS(1,4)

c For s = smax meters (at inner target plate)
     
      TARGETVB = SQRT((KTEDS(IDDS(IR,1))+KTIDS(IDDS(IR,1)))
     >           *1.602E-19/(CRMB*1.67E-27))*QTIM

      KVHGS(NKS(IR),IR) = 
     >          (TARGETVB-KVHS(NKS(IR),IR))/
     >          (KSMAXS(IR)-KSS(NKS(IR),IR))

      LAMBDA1(IZ) = 0.23995*Coulomb_log*REAL(IZ*IZ)*
     >              REAL(RIZB*RIZB)/(CRMI*CRMI)
 
      ALPHAI(NKS(IR),IR) = SQRT(CRMB/(2*KTIDS(IDDS(IR,1))*EMI))
          
      CC1(NKS(IR),IR) = 1.8024009E17*KTIDS(IDDS(IR,1))*
     >     KFIDS(IDDS(IR,1))/(Coulomb_log*REAL(RIZB*RIZB*RIZB*RIZB)*
     >     KNDS(IDDS(IR,1))*FACT)

      CC2(NKS(IR),IR)=-7.361216E12*KVHGS(NKS(IR),IR)*SQRT(CRMB*
     >     KTIDS(IDDS(IR,1))*KTIDS(IDDS(IR,1))*KTIDS(IDDS(IR,1)))/
     >     (Coulomb_log*RIZB*RIZB*RIZB*(1.2*RIZB + 0.8485)*
     >     KNDS(IDDS(IR,1))*QTIM)   

      CHIpara = ALPHAI(NKS(IR),IR)*(-KVDS(IDDS(IR,1)))
      LAMBDA2 = LAMBDA1(IZ)*KNDS(IDDS(IR,1))

      CALL Coulomb_Coll(KKpar,DDparpar,CHIpara,LAMBDA2,
     >                  NKS(IR),IR,Kaa,Kab,Kac,Daa,Dab,Dac)

      KVALS(nks(IR)+inc,1) = Kaa * CRMI * amu 
      KVALS(nks(IR)+inc,2) = Kab * CRMI * amu
      KVALS(nks(IR)+inc,3) = Kac * CRMI * amu
      KVALS(nks(IR)+inc,4) = KVALS(nks(IR)+inc,1) + 
     >           KVALS(nks(IR)+inc,2) + KVALS(nks(IR)+inc,3) 

*      WRITE(0,*)KVALS(nks(IR)+inc,1)
*      WRITE(0,*)KVALS(nks(IR)+inc,2)
*      WRITE(0,*)KVALS(nks(IR)+inc,3)
*      WRITE(0,*)KVALS(nks(IR)+inc,4)
c
        CALL Vbgrad
        CALL Coeff(NIZS)
c        
        DO IK = 1,NKS(IR)
         
         KOUTS(IK+in) = KSS(IK,IR)
         KWIDS(IK+in) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
c
         CHIpara = ALPHAI(IK,IR)*(-KVHS(IK,IR))
         LAMBDA2 = LAMBDA1(IZ)*KNBS(IK,IR)
c         CALL Coulomb_Coll2(KKpar,DDparpar,CHIpara,LAMBDA2,
c     >                    IK,IR,Kaa,Kab,Kac,Daa,Dab,Dac)
         CALL Coulomb_Coll(KKpar,DDparpar,CHIpara,LAMBDA2,
     >                    IK,IR,Kaa,Kab,Kac,Daa,Dab,Dac)
         KVALS(IK+in,1) = Kaa * CRMI * amu
         KVALS(IK+in,2) = Kab * CRMI * amu
         KVALS(IK+in,3) = Kac * CRMI * amu
         KVALS(IK+in,4) = KVALS(IK+in,1)+KVALS(IK+in,2)+KVALS(IK+in,3) 

*         WRITE(0,*)IK+in,KVALS(IK+in,1),KVALS(IK+in,2),KVALS(IK+in,3),
*     >                   KVALS(IK+in,4),CHIpara,LAMBDA2

         WRITE(iplot,'(2I4,8G12.4)') IK,IR,KVALS(IK+IN,1),
     +      KVALS(IK+IN,2),
     +      KVALS(IK+IN,3),KVALS(IK+IN,4),KVALS(IK+IN,5),KSS(IK,IR),
     +      RS(IK,IR),ZS(IK,IR) 
c
c          IF(FILED.EQ.1)THEN
c            WRITE(121,9)KOUTS(IK+in),KVALS(IK+in,1),
c     >                  KVALS(IK+in,2),KVALS(IK+in,3),
c     >                  KVALS(IK+in,4)
c    9       FORMAT(5(1X,G11.5))
c          ENDIF
c 
        END DO

c For S=0 meters (at target plate), the Y-axis values will be taken 
c as those calculated for the first grid cell: 
                
*          KVALS(1,1) = KVALS(1+in,1)
*          KVALS(1,2) = KVALS(1+in,2)
*          KVALS(1,3) = KVALS(1+in,3) 
*          KVALS(1,4) = KVALS(1,1) + KVALS(1,2) + KVALS(1,3) 

c For S=Smax (at the other target plate), the Y-axis values will be
c taken as those calculated for the last grid cell: 

*          KVALS(nks(ir)+inc,1) = KVALS(nks(ir)+in,1)
*          KVALS(nks(ir)+inc,2) = KVALS(nks(ir)+in,2)
*          KVALS(nks(ir)+inc,3) = KVALS(nks(ir)+in,3)
*          KVALS(nks(ir)+inc,4) = KVALS(nks(ir)+inc,1) +
*     +                KVALS(nks(ir)+inc,2) + KVALS(nks(ir)+inc,3)  

         CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >    4,99,KOUTS(1),KOUTS(NKS(IR)+inc),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0)
c
        enldist = mgst*kouts(nks(ir)+inc)
c
        IF (clsup.eq.1) THEN
          CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >      4,99,KOUTS(1),enldist   ,-HI,HI,IGNORS,ITEC,AVS,
     >      NAVS,
     >      JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0)

          CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >      4,99,KOUTS(NKS(IR)+inc)-enldist,KOUTS(NKS(IR)+inc),
     >      -HI,HI,IGNORS,ITEC,AVS,NAVS,
     >      JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0)
        ENDIF

      ENDIF
c
C
c-----------------------------------------------------------------------
c  Force Balance Plots for Reiser Formulation when Relative Velocity = 0
c-----------------------------------------------------------------------
c
      IF (IREF.EQ.577) THEN

        ELABS(1) = 'FF  FF     '
        ELABS(2) = 'FIG FIG    '
        ELABS(3) = 'FvbgFvbg   '
        ELABS(4) = 'FnetFnet   '

        XLAB = '   S  (M)'
        YLAB = ' FORCE(N)  '

        NVIEW = 'REISER: FORCE BALANCE WHEN CHI = 0'
        
        READ (GRAPH(32:35),'(I4)') IZ
        READ (GRAPH(40:44),'(I4)') IR
        READ (GRAPH(48:49),'(I2)') FILED
c
c        IF (FILED.EQ.1) THEN
c          WRITE(121,*)'577-FORCE BALANCE, CHI = 0'
c        ENDIF
c
        NPLOTS = NPLOTS + 1

        WRITE (REF,'(''Ring: '',I3,'', K: '',F7.4,'' ,State: '',I3)') 
     +        IR,KKS(IR),IZ

        WRITE (IPLOT,9012) NPLOTS,REF

        write(iplot,'(2x,''IK'',2x,''IR'',5x,''FF'',5x,5x'//
     >                   ',''FIG'',4x,5x,''FEG'',4x,5x'//
     >                    ',''FE'',5x,5x,''NET'',4x,6x,''S'',5x'//
     >                    ',6x,''R'',5x,6x,''Z'')')

        CALL rzero (kvals, maxnks*maxngs)

c       Modify the following to include values for S=0 ... even if this
c       is not part of the grid. Assume that if one end is not a grid
c       point then neither is the other.

        if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
           in  = 0
           inc = 0
        else
           in       = 1
           inc      = 2
           kouts(1) = 0.0

           KWIDS(1)         = kss(1,ir)
           kouts(nks(ir)+2) = ksmaxs(ir)
           kwids(nks(ir)+2) = ksmaxs(ir) - kss(nks(ir),ir)
        endif

        CALL Vbgrad
        CALL Coeff(NIZS)   

c For s = 0 meters (at outer target plate)

      FACT = QTIM * QTIM * EMI / CRMI
 
      TARGETVB = -SQRT((KTEDS(IDDS(IR,2))+KTIDS(IDDS(IR,2)))
     >           *1.602E-19/(CRMB*1.67E-27))*QTIM

      KVHGS(1,IR) = (KVHS(1,IR)-TARGETVB)/KSS(1,IR)

      LAMBDA1(IZ) = 0.23995*Coulomb_log*REAL(IZ*IZ)*
     >              REAL(RIZB*RIZB)/(CRMI*CRMI)
 
      ALPHAI(1,IR) = SQRT(CRMB/(2*KTIDS(IDDS(IR,2))*EMI))
          
      CC1(1,IR) = 1.8024009E17*KTIDS(IDDS(IR,2))*KFIDS(IDDS(IR,2))/
     >   (Coulomb_log*REAL(RIZB*RIZB*RIZB*RIZB)*KNDS(IDDS(IR,2))*FACT)

      CC2(1,IR)=-7.361216E12*KVHGS(1,IR)*SQRT(CRMB*
     >      KTIDS(IDDS(IR,2))*KTIDS(IDDS(IR,2))*KTIDS(IDDS(IR,2)))/
     >      (Coulomb_log*RIZB*RIZB*RIZB*(1.2*RIZB + 0.8485)*
     >      KNDS(IDDS(IR,2))*QTIM) 


      CHIpara = 0.0
      LAMBDA2 = LAMBDA1(IZ)*KNDS(IDDS(IR,2))

      CALL Coulomb_Coll(KKpar,DDparpar,CHIpara,LAMBDA2,
     >                   1,IR,Kaa,Kab,Kac,Daa,Dab,Dac)

      KVALS(1,1) = Kaa * CRMI * amu 
      KVALS(1,2) = Kab * CRMI * amu
      KVALS(1,3) = Kac * CRMI * amu
      KVALS(1,4) = KVALS(1,1) + KVALS(1,2) + KVALS(1,3)

*      WRITE(0,*)'577:'
*      WRITE(0,*)KVALS(1,1)
*      WRITE(0,*)KVALS(1,2)
*      WRITE(0,*)KVALS(1,3)
*      WRITE(0,*)KVALS(1,4) 

c For s = smax meters (at inner target plate)
     
      TARGETVB = SQRT((KTEDS(IDDS(IR,1))+KTIDS(IDDS(IR,1)))
     >           *1.602E-19/(CRMB*1.67E-27))*QTIM

      KVHGS(NKS(IR),IR) = 
     >          (TARGETVB-KVHS(NKS(IR),IR))/
     >          (KSMAXS(IR)-KSS(NKS(IR),IR))

      LAMBDA1(IZ) = 0.23995*Coulomb_log*REAL(IZ*IZ)*
     >              REAL(RIZB*RIZB)/(CRMI*CRMI)
 
      ALPHAI(NKS(IR),IR) = SQRT(CRMB/(2*KTIDS(IDDS(IR,1))*EMI))
          
      CC1(NKS(IR),IR) = 1.8024009E17*KTIDS(IDDS(IR,1))*
     >     KFIDS(IDDS(IR,1))/(Coulomb_log*REAL(RIZB*RIZB*RIZB*RIZB)*
     >     KNDS(IDDS(IR,1))*FACT)

      CC2(NKS(IR),IR)=-7.361216E12*KVHGS(NKS(IR),IR)*SQRT(CRMB*
     >     KTIDS(IDDS(IR,1))*KTIDS(IDDS(IR,1))*KTIDS(IDDS(IR,1)))/
     >     (Coulomb_log*RIZB*RIZB*RIZB*(1.2*RIZB + 0.8485)*
     >     KNDS(IDDS(IR,1))*QTIM)    

      CHIpara = 0.0
      LAMBDA2 = LAMBDA1(IZ)*KNDS(IDDS(IR,1))

      CALL Coulomb_Coll(KKpar,DDparpar,CHIpara,LAMBDA2,
     >                  NKS(IR),IR,Kaa,Kab,Kac,Daa,Dab,Dac)

      KVALS(nks(IR)+inc,1) = Kaa * CRMI * amu 
      KVALS(nks(IR)+inc,2) = Kab * CRMI * amu
      KVALS(nks(IR)+inc,3) = Kac * CRMI * amu
      KVALS(nks(IR)+inc,4) = KVALS(nks(IR)+inc,1) + 
     >           KVALS(nks(IR)+inc,2) + KVALS(nks(IR)+inc,3) 

*      WRITE(0,*)KVALS(nks(IR)+inc,1)
*      WRITE(0,*)KVALS(nks(IR)+inc,2)
*      WRITE(0,*)KVALS(nks(IR)+inc,3)
*      WRITE(0,*)KVALS(nks(IR)+inc,4)

        CALL Vbgrad
        CALL Coeff(NIZS)     
        
        DO IK = 1,NKS(IR)
         
         KOUTS(IK+in) = KSS(IK,IR)
         KWIDS(IK+in) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))

         CHIpara = 0.0
         LAMBDA2 = LAMBDA1(IZ)*KNBS(IK,IR)
c         CALL Coulomb_Coll2(KKpar,DDparpar,CHIpara,LAMBDA2,
c     >                    IK,IR,Kaa,Kab,Kac,Daa,Dab,Dac)
         CALL Coulomb_Coll(KKpar,DDparpar,CHIpara,LAMBDA2,
     >                    IK,IR,Kaa,Kab,Kac,Daa,Dab,Dac)
         KVALS(IK+in,1) = Kaa * CRMI * amu
         KVALS(IK+in,2) = Kab * CRMI * amu
         KVALS(IK+in,3) = Kac * CRMI * amu
         KVALS(IK+in,4) = KVALS(IK+in,1)+KVALS(IK+in,2)+KVALS(IK+in,3) 


         WRITE(iplot,'(2I4,8G12.4)') IK,IR,KVALS(IK+IN,1),
     +      KVALS(IK+IN,2),
     +      KVALS(IK+IN,3),KVALS(IK+IN,4),KVALS(IK+IN,5),KSS(IK,IR),
     +      RS(IK,IR),ZS(IK,IR)
c
c          IF(FILED.EQ.1)THEN
c            WRITE(121,25)KOUTS(IK+in),KVALS(IK+in,1),
c     >                   KVALS(IK+in,2),KVALS(IK+in,3),
c     >                   KVALS(IK+in,4)
c   25       FORMAT(5(1X,G11.5))
c          ENDIF
c
        END DO

c For S=0 meters (at target plate), the Y-axis values will be taken 
c as those calculated for the first grid cell: 
                
*          KVALS(1,1) = KVALS(1+in,1)
*          KVALS(1,2) = KVALS(1+in,2)
*          KVALS(1,3) = KVALS(1+in,3) 
*          KVALS(1,4) = KVALS(1,1) + KVALS(1,2) + KVALS(1,3) 

c For S=Smax (at the other target plate), the Y-axis values will be
c taken as those calculated for the last grid cell: 

*          KVALS(nks(ir)+inc,1) = KVALS(nks(ir)+in,1)
*          KVALS(nks(ir)+inc,2) = KVALS(nks(ir)+in,2)
*          KVALS(nks(ir)+inc,3) = KVALS(nks(ir)+in,3)
*          KVALS(nks(ir)+inc,4) = KVALS(nks(ir)+inc,1) +
*     +                 KVALS(nks(ir)+inc,2) + KVALS(nks(ir)+inc,3)  

         CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >    4,99,KOUTS(1),KOUTS(NKS(IR)+inc),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0)
c
        enldist = mgst*kouts(nks(ir)+inc)
c
        IF (clsup.eq.1) THEN
          CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >      4,99,KOUTS(1),enldist   ,-HI,HI,IGNORS,ITEC,AVS,
     >      NAVS,
     >      JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0)

          CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >      4,99,KOUTS(NKS(IR)+inc)-enldist,KOUTS(NKS(IR)+inc),
     >      -HI,HI,IGNORS,ITEC,AVS,NAVS,
     >      JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0)
        ENDIF

      ENDIF
c
c-----------------------------------------------------------------------
c  Force Balance Plots for Reiser Formulation with Selected Impurity
c  Velocity, Vz in (m/s)
c-----------------------------------------------------------------------
c
      IF (IREF.EQ.578) THEN

        ELABS(1) = 'FF  FF     '
        ELABS(2) = 'FIG FIG    '
        ELABS(3) = 'FvbgFvbg   '
        ELABS(4) = 'FnetFnet   '

        XLAB = '   S  (M)'
        YLAB = ' FORCE(N)  '

        NVIEW = 'REISER: FORCE BALANCE: Vz specified'
        
        READ (GRAPH(20:26),'(F7.1)') Vimp
        READ (GRAPH(32:35),'(I4)') IZ
        READ (GRAPH(40:44),'(I4)') IR
        READ (GRAPH(48:49),'(I2)') FILED
c
c        IF (FILED.EQ.1) THEN
c          WRITE(121,*)'578-FORCE BALANCE, VZ SELECTED'
c        ENDIF
c
        NPLOTS = NPLOTS + 1

*        write(0,*)'vimp',vimp*QTIM,'vb',KVHS(1,8)

        WRITE (REF,'(''Ring: '',I3,'', Vz: '',E7.1,'' ,State: '',I3)') 
     +        IR,Vimp,IZ

        WRITE (IPLOT,9012) NPLOTS,REF

        write(iplot,'(2x,''IK'',2x,''IR'',5x,''FF'',5x,5x'//
     >                   ',''FIG'',4x,5x,''FEG'',4x,5x'//
     >                    ',''FE'',5x,5x,''NET'',4x,6x,''S'',5x'//
     >                    ',6x,''R'',5x,6x,''Z'')')

        CALL rzero (kvals, maxnks*maxngs)

c       Modify the following to include values for S=0 ... even if this
c       is not part of the grid. Assume that if one end is not a grid
c       point then neither is the other.
c
c       !!!
c          Pete's comment: These are the target surface values and those 
c          of the first or last grid cells along a ring.  This allows for 
c          a continous plot from the first and last grid cell to the solid 
c          target surfaces.  Note: the velocity values calculated for S=0 
c          and S=Smax are determined from the target temperatures using the 
c          formula for the plasma sound speed, cs.  This is fine for SOL 
c          options 12 and 13, however, this may not give a smooth velocity 
c          profile leading up to the targets when the target velocity can 
c          be selected independenlty of the target temperatures as in the 
c          case of SOL option 7.
c       !!!      

        if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
           in  = 0
           inc = 0
        else
           in       = 1
           inc      = 2
           kouts(1) = 0.0

           KWIDS(1)         = kss(1,ir)
           kouts(nks(ir)+2) = ksmaxs(ir)
           kwids(nks(ir)+2) = ksmaxs(ir) - kss(nks(ir),ir)
        endif

c For s = 0 meters (at outer target plate)

      FACT = QTIM * QTIM * EMI / CRMI

      TARGETVB = -SQRT((KTEDS(IDDS(IR,2))+KTIDS(IDDS(IR,2)))
     >           *1.602E-19/(CRMB*1.67E-27))*QTIM

      KVHGS(1,IR) = (KVHS(1,IR)-TARGETVB)/KSS(1,IR)

      LAMBDA1(IZ) = 0.23995*Coulomb_log*REAL(IZ*IZ)*
     >              REAL(RIZB*RIZB)/(CRMI*CRMI)
 
      ALPHAI(1,IR) = SQRT(CRMB/(2*KTIDS(IDDS(IR,2))*EMI))
          
*      CC1(1,IR) = 1.20246E16*KTIDS(IDDS(IR,2))*KFIDS(IDDS(IR,2))/
*     >          (REAL(RIZB*RIZB*RIZB*RIZB)*KNDS(IDDS(IR,2))*FACT)
*
*      CC2(1,IR) = -4.9113E11*KVHGS(1,IR)*SQRT(CRMB*
*     >    KTIDS(IDDS(IR,2))*KTIDS(IDDS(IR,2))*KTIDS(IDDS(IR,2)))/
*     > (RIZB*RIZB*RIZB*(1.2*RIZB + 0.8485)*KNDS(IDDS(IR,2))*QTIM) 

      CC1(1,IR) = 1.8024009E17*KTIDS(IDDS(IR,2))*KFIDS(IDDS(IR,2))/
     >   (Coulomb_log*REAL(RIZB*RIZB*RIZB*RIZB)*KNDS(IDDS(IR,2))*FACT)

      CC2(1,IR)=-7.361216E12*KVHGS(1,IR)*SQRT(CRMB*
     >      KTIDS(IDDS(IR,2))*KTIDS(IDDS(IR,2))*KTIDS(IDDS(IR,2)))/
     >      (Coulomb_log*RIZB*RIZB*RIZB*(1.2*RIZB + 0.8485)*
     >      KNDS(IDDS(IR,2))*QTIM) 

*      WRITE(121,*)CC1(1,IR),CC2(1,IR),KVHGS(1,IR),KTEDS(IDDS(IR,2)),
*     >            KTIDS(IDDS(IR,2)),TARGETVB/QTIM

      CHIpara = ALPHAI(1,IR)*(Vimp-KVDS(IDDS(IR,2)))
      LAMBDA2 = LAMBDA1(IZ)*KNDS(IDDS(IR,2))

      CALL Coulomb_Coll(KKpar,DDparpar,CHIpara,LAMBDA2,
     >                   1,IR,Kaa,Kab,Kac,Daa,Dab,Dac)

      KVALS(1,1) = Kaa * CRMI * amu 
      KVALS(1,2) = Kab * CRMI * amu
      KVALS(1,3) = Kac * CRMI * amu
      KVALS(1,4) = KVALS(1,1) + KVALS(1,2) + KVALS(1,3)

*      WRITE(0,*)'578:'
*      WRITE(0,*)KVALS(1,1)
*      WRITE(0,*)KVALS(1,2)
*      WRITE(0,*)KVALS(1,3)
*      WRITE(0,*)KVALS(1,4)

c For s = smax meters (at inner target plate)
     
      TARGETVB = SQRT((KTEDS(IDDS(IR,1))+KTIDS(IDDS(IR,1)))
     >           *1.602E-19/(CRMB*1.67E-27))*QTIM

      KVHGS(NKS(IR),IR) = 
     >          (TARGETVB-KVHS(NKS(IR),IR))/
     >          (KSMAXS(IR)-KSS(NKS(IR),IR))

      LAMBDA1(IZ) = 0.23995*Coulomb_log*REAL(IZ*IZ)*
     >              REAL(RIZB*RIZB)/(CRMI*CRMI)
 
      ALPHAI(NKS(IR),IR) = SQRT(CRMB/(2*KTIDS(IDDS(IR,1))*EMI))
          
*      CC1(NKS(IR),IR) = 
*     >          1.20246E16*KTIDS(IDDS(IR,1))*KFIDS(IDDS(IR,1))/
*     >          (REAL(RIZB*RIZB*RIZB*RIZB)*KNDS(IDDS(IR,1))*FACT)
*
*      CC2(NKS(IR),IR) = -4.9113E11*KVHGS(1,IR)*SQRT(CRMB*
*     >    KTIDS(IDDS(IR,1))*KTIDS(IDDS(IR,1))*KTIDS(IDDS(IR,1)))/
*     > (RIZB*RIZB*RIZB*(1.2*RIZB + 0.8485)*KNDS(IDDS(IR,1))*QTIM)

      CC1(NKS(IR),IR) = 1.8024009E17*KTIDS(IDDS(IR,1))*
     >     KFIDS(IDDS(IR,1))/(Coulomb_log*REAL(RIZB*RIZB*RIZB*RIZB)*
     >     KNDS(IDDS(IR,1))*FACT)

      CC2(NKS(IR),IR)=-7.361216E12*KVHGS(NKS(IR),IR)*SQRT(CRMB*
     >     KTIDS(IDDS(IR,1))*KTIDS(IDDS(IR,1))*KTIDS(IDDS(IR,1)))/
     >     (Coulomb_log*RIZB*RIZB*RIZB*(1.2*RIZB + 0.8485)*
     >     KNDS(IDDS(IR,1))*QTIM)  

*      WRITE(121,*)CC1(NKS(IR),IR),CC2(NKS(IR),IR),KVHGS(NKS(IR),IR),
*     >            KTEDS(IDDS(IR,1)),KTIDS(IDDS(IR,1)),TARGETVB/QTIM

      CHIpara = ALPHAI(NKS(IR),IR)*(Vimp-KVDS(IDDS(IR,1)))
      LAMBDA2 = LAMBDA1(IZ)*KNDS(IDDS(IR,1))

      CALL Coulomb_Coll(KKpar,DDparpar,CHIpara,LAMBDA2,
     >                  NKS(IR),IR,Kaa,Kab,Kac,Daa,Dab,Dac)

      KVALS(nks(IR)+inc,1) = Kaa * CRMI * amu 
      KVALS(nks(IR)+inc,2) = Kab * CRMI * amu
      KVALS(nks(IR)+inc,3) = Kac * CRMI * amu
      KVALS(nks(IR)+inc,4) = KVALS(nks(IR)+inc,1) + 
     >           KVALS(nks(IR)+inc,2) + KVALS(nks(IR)+inc,3) 

*      WRITE(0,*)KVALS(nks(IR)+inc,1)
*      WRITE(0,*)KVALS(nks(IR)+inc,2)
*      WRITE(0,*)KVALS(nks(IR)+inc,3)
*      WRITE(0,*)KVALS(nks(IR)+inc,4)

      CALL Vbgrad
      CALL Coeff(NIZS)

*      DO IK = 1,NKS(IR)
*        WRITE(121,*)IK,IR,CC1(IK,IR),CC2(IK,IR),KVHGS(IK,IR),
*     >              KVHS(IK,IR)/QTIM
*      END DO
        
        DO IK = 1,NKS(IR)
         
         KOUTS(IK+in) = KSS(IK,IR)
         KWIDS(IK+in) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))

         CHIpara = ALPHAI(IK,IR)*((Vimp*QTIM)-KVHS(IK,IR))
         LAMBDA2 = LAMBDA1(IZ)*KNBS(IK,IR)

         CALL Coulomb_Coll(KKpar,DDparpar,CHIpara,LAMBDA2,
     >                    IK,IR,Kaa,Kab,Kac,Daa,Dab,Dac)

         KVALS(IK+in,1) = Kaa * CRMI * amu
         KVALS(IK+in,2) = Kab * CRMI * amu
         KVALS(IK+in,3) = Kac * CRMI * amu
         KVALS(IK+in,4) = KVALS(IK+in,1)+KVALS(IK+in,2)+KVALS(IK+in,3) 


         WRITE(iplot,'(2I4,8G12.4)') IK,IR,KVALS(IK+IN,1),
     +      KVALS(IK+IN,2),
     +      KVALS(IK+IN,3),KVALS(IK+IN,4),KVALS(IK+IN,5),KSS(IK,IR),
     +      RS(IK,IR),ZS(IK,IR)
c
c          IF(FILED.EQ.1)THEN
c            WRITE(121,11)KOUTS(IK+in),KVALS(IK+in,1),
c     >                   KVALS(IK+in,2),KVALS(IK+in,3),
c     >                   KVALS(IK+in,4)
c   11       FORMAT(5(1X,G11.5))
c          ENDIF
c 
        END DO

c For S=0 meters (at target plate), the Y-axis values will be taken 
c as those calculated for the first grid cell: 
                
*          KVALS(1,1) = KVALS(1+in,1)
*          KVALS(1,2) = KVALS(1+in,2)
*          KVALS(1,3) = KVALS(1+in,3) 
*          KVALS(1,4) = KVALS(1,1) + KVALS(1,2) + KVALS(1,3)

c For S=Smax (at the other target plate), the Y-axis values will be
c taken as those calculated for the last grid cell: 

*          KVALS(nks(ir)+inc,1) = KVALS(nks(ir)+in,1)
*          KVALS(nks(ir)+inc,2) = KVALS(nks(ir)+in,2)
*          KVALS(nks(ir)+inc,3) = KVALS(nks(ir)+in,3)
*          KVALS(nks(ir)+inc,4) = KVALS(nks(ir)+inc,1) +
*     +                 KVALS(nks(ir)+inc,2) + KVALS(nks(ir)+inc,3)  

         CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >    4,99,KOUTS(1),KOUTS(NKS(IR)+inc),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0)
c
        enldist = mgst*kouts(nks(ir)+inc)
c
        IF (clsup.eq.1) THEN
          CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >      4,99,KOUTS(1),enldist   ,-HI,HI,IGNORS,ITEC,AVS,
     >      NAVS,
     >      JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0)

          CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >      4,99,KOUTS(NKS(IR)+inc)-enldist,KOUTS(NKS(IR)+inc),
     >      -HI,HI,IGNORS,ITEC,AVS,NAVS,
     >      JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0)
        ENDIF

      ENDIF
c
c-----------------------------------------------------------------------
c     COMBINED FRICTIONAL AND THERMAL FORCE PLOTS 
c-----------------------------------------------------------------------
c
      IF (IREF.EQ.579) THEN

        ELABS(1) = 'FkinFfriction + Fthermal(KINETIC)'
        ELABS(2) = 'FfluFfriction + Fthermal(FLUID)'
*         ELABS(1) = 'SIN SIN   '
*         ELABS(2) = 'COS COS   '
 
        XLAB = '   CHI   '
        YLAB = ' FORCE(N)'

        NVIEW = 'FRICTION + THERMAL(KIN./FLU.) FORCES'
*        PLANE = 'Plane123456789 123456789 123456789 123456789 123456789'
     
*        WRITE(0,*)GRAPH,ZB,ZI,MB,MI,NB,TB,TBG,IERR

        NPLOTS = NPLOTS + 1

        WRITE (REF,'(''Tbgrd '',F6.2,'',TB '',F6.2,'' ,State: '',I3)') 
     +        TBG,TB,ZI


        WRITE (IPLOT,9012) NPLOTS,REF

        CALL RZERO(YVALUE,100*MAXNGS)
        CALL RZERO(XVALUE,100)
        CALL RZERO(WIDTHS,100)
c
*        CALL RZERO(ROOTS,5)
        DO kounter = 1,5
         ROOTS(kounter) = 99.9999
        END DO

        CHIpara  = 0.0
        KKpar    = 0.0
        kounter  = 0
        jcounter = 0
        negvalue = 0
        posvalue = 0
        FORCE    = 4 
 
        DO CHIpara = -3.0,3.1,0.1
*         DO CHIpara = -4.0,4.1,0.1

           kounter = kounter + 1
           IF(negvalue.eq.-1.and.posvalue.eq.1)THEN
             jcounter = jcounter + 1
             ROOTS(jcounter)=BISCT(CHIpara-0.2,CHIpara-0.1,
     >                       VBG,TBG,TB,NB,MB,MI,ZI,ZB,XPER,
     >                       FORCE)
             negvalue=0
             posvalue=0
           ENDIF

           CALL Coulomb_Coll2(Kaa,Kab,Kaba,Kabb,Kac,KBetaf,
     >       Daa,Dab,Dac,Coulomb_log,Frictionf,CHIpara,Xper,
     >       VBG,TBG,TB,NB,MB,MI,ZI,ZB)
c
c           IF(CHIpara.EQ.0.0)THEN
c              WRITE(0,*)'k=',kounter
c           ENDIF
c           
           XVALUE(kounter)   = CHIpara
           WIDTHS(kounter)   = 0.1
*           YVALUE(kounter,1) = COS(CHIpara)
*           YVALUE(kounter,2) = COS(CHIpara)
*           YVALUE(kounter,1) = (Kaa+Kab+Kac)*REAL(MI)*amu
           YVALUE(kounter,1) = (Kaa+Kab)*REAL(MI)*amu
           YVALUE(kounter,2) = (Kaa+KBetaf)*REAL(MI)*amu

          IF(YVALUE(kounter,1).gt.0.0)THEN
            posvalue = 1
          ELSE
            negvalue = -1
          ENDIF  
           
        END DO
        CHIpara=0.00000

       WRITE (ANLY,'(''ROOTS(1-3)'',f8.4,f8.4,f8.4)') 
     +        ROOTS(1), ROOTS(2), ROOTS(3)
       WRITE (PLANE,'(''ROOTS(4-5)'',f8.4,f8.4)') 
     +        ROOTS(4), ROOTS(5)


       CALL DRAW(XVALUE,WIDTHS,YVALUE,100,61,anly,2,99,-3.0,3.0,
     >           -HI,HI,IGNORS,0,AVS,NAVS,JOB,TITLE,XLAB,YLAB,ELABS,
     >           REF,NVIEW,PLANE,TABLE,1,6,1.0)

      ENDIF
C
c------------------------------------------------------------------------
c FRICTIONAL, THERMAL, AND VISCOUS FORCE PLOTS (SHOWING COMP. SEPERATELY) 
c------------------------------------------------------------------------
c
      IF (IREF.EQ.580) THEN

        ELABS(1) = 'FFk FRICTION(KINETIC)  FFk'
        ELABS(2) = 'FFf FRICTION(FLUID)    FFf'
        ELABS(3) = 'FIGkTHERMAL (KINETIC) FIGk'
        ELABS(4) = 'FIGfTHERMAL (FLUID)   FIGf'
        ELABS(5) = 'FvbgVELOCITY GRADIENT Fvbg'
        ELABS(6) = 'Fnk Fnk=FFk+FIGk+Fvbg  Fnk'
        ELABS(7) = 'Fnf Fnf=FFf+FIGf       Fnf'
 
        XLAB = '   CHI   '
        YLAB = ' FORCE(N)'

        NVIEW = 'COMPARING COMPONENT FORCES'
     
        NPLOTS = NPLOTS + 1

        WRITE (REF,'(''Tbgrd '',F6.2,'',TB '',F6.2,'' ,State: '',I3)') 
     +        TBG,TB,ZI


        WRITE (IPLOT,9012) NPLOTS,REF

        CHIpara = 0.0
        KKpar   = 0.0
        kounter = 0
        jcounter = 0
        negvalue = 0
        posvalue = 0
        FORCE    = 3 

        CALL RZERO(YVALUE,100*MAXNGS)
        CALL RZERO(XVALUE,100)
        CALL RZERO(WIDTHS,100)

        DO kounter = 1,5
         ROOTS(kounter) = 99.9999
        END DO

        kounter = 0
c
c       IF (FILED.EQ.1) THEN
c         WRITE(121,*)'580-DRIFT COMPONENTS'
c       ENDIF
c 
        DO CHIpara = -3.0,3.1,0.1
           kounter = kounter + 1

            IF(negvalue.eq.-1.and.posvalue.eq.1)THEN
             jcounter = jcounter + 1
             ROOTS(jcounter)=BISCT(CHIpara-0.2,CHIpara-0.1,
     >                       VBG,TBG,TB,NB,MB,MI,ZI,ZB,XPER,
     >                       FORCE)
             negvalue=0
             posvalue=0
            ENDIF

           CALL Coulomb_Coll2(Kaa,Kab,Kaba,Kabb,Kac,KBetaf,
     >       Daa,Dab,Dac,Coulomb_log,Frictionf,CHIpara,Xper,
     >       VBG,TBG,TB,NB,MB,MI,ZI,ZB) 
           
           XVALUE(kounter)   = CHIpara
           WIDTHS(kounter)   = 0.1
           YVALUE(kounter,1) = Kaa * REAL(MI)*amu
           YVALUE(kounter,2) = Frictionf * REAL(MI)*amu
           YVALUE(kounter,3) = Kab * REAL(MI)*amu
           YVALUE(kounter,4) = Kbetaf * REAL(MI)*amu
           YVALUE(kounter,5) = Kac * REAL(MI)*amu
           YVALUE(kounter,6) = (Kaa + Kab + Kac) * REAL(MI)
     >                                           * amu
           YVALUE(kounter,7) = YVALUE(kounter,2) +
     >                         YVALUE(kounter,4)
c           IF(FILED.EQ.1)THEN
c            WRITE(121,12)XVALUE(kounter),YVALUE(kounter,1),
c     >                 YVALUE(kounter,2),YVALUE(kounter,3),
c     >                 YVALUE(kounter,4),YVALUE(kounter,5),
c     >                 YVALUE(kounter,6),YVALUE(kounter,7)
c   12       FORMAT(8(1X,G11.4))
c           ENDIF

          IF(YVALUE(kounter,5).gt.0.0)THEN
            posvalue = 1
          ELSE
            negvalue = -1
          ENDIF  
           
        END DO

       WRITE (ANLY,'(''ROOTS(1-3)'',f8.4,f8.4,f8.4)') 
     +        ROOTS(1), ROOTS(2), ROOTS(3)
       WRITE (PLANE,'(''ROOTS(4-5)'',f8.4,f8.4)') 
     +        ROOTS(4), ROOTS(5)

       CALL DRAW(XVALUE,WIDTHS,YVALUE,100,61,anly,7,99,-3.0,3.0,
     >           -HI,HI,IGNORS,0,AVS,NAVS,JOB,TITLE,XLAB,YLAB,ELABS,
     >           REF,NVIEW,PLANE,TABLE,1,6,1.0)

      ENDIF
c
c-----------------------------------------------------------------------
c     DIFFUSIVE FORCE PLOTS (SHOWING COMP. SEPERATELY) 
c-----------------------------------------------------------------------
c
      IF (IREF.EQ.581) THEN

        ELABS(1) = 'Dnb DIFFUSIVE TERMS: Dnb'
        ELABS(2) = 'DTbg                 DTbg'
        ELABS(3) = 'Dvbg                 Dvbg'
        ELABS(4) = 'Dnet                 Dnet'
        ELABS(5) = 'DDIV                 DDIV'
 
        XLAB = '    CHI    '
        YLAB = ' (m^2/s^3) '

        NVIEW = 'COMPARING DIFFUSIVE COMPONENTS'
     
        NPLOTS = NPLOTS + 1

        WRITE (REF,'(''Tbgrd '',F6.2,'',TB '',F6.2,'' ,State: '',I3)') 
     +        TBG,TB,ZI


        WRITE (IPLOT,9012) NPLOTS,REF

        CHIpara = 0.0
        KKpar   = 0.0
        kounter = 0
        jcounter = 0
        negvalue = 0
        posvalue = 0
        FORCE    = 3 

        TAUpara = 1.47e13*MI*SQRT(TB/MB)/(NB*ZI*ZI*ZB*ZB*15.0)
        DDIV = EMI*2.0/(MI*TAUpara)

        CALL RZERO(YVALUE,100*MAXNGS)
        CALL RZERO(XVALUE,100)
        CALL RZERO(WIDTHS,100)

        DO kounter = 1,5
         ROOTS(kounter) = 99.9999
        END DO

        kounter = 0

c        IF (FILED.EQ.1) THEN
c          WRITE(121,*)'581-DIFFUSION COMPONENTS'
c        ENDIF
 
        DO CHIpara = -3.0,3.1,0.1
           kounter = kounter + 1

*            IF(negvalue.eq.-1.and.posvalue.eq.1)THEN
*             jcounter = jcounter + 1
*             ROOTS(jcounter)=BISCT(CHIpara-0.2,CHIpara-0.1,
*     >                       VBG,TBG,TB,NB,MB,MI,ZI,ZB,XPER,
*     >                       FORCE)
*             negvalue=0
*             posvalue=0
*            ENDIF

           CALL Coulomb_Coll2(Kaa,Kab,Kaba,Kabb,Kac,KBetaf,
     >       Daa,Dab,Dac,Coulomb_log,Frictionf,CHIpara,Xper,
     >       VBG,TBG,TB,NB,MB,MI,ZI,ZB) 
      
           XVALUE(kounter)   = CHIpara
           WIDTHS(kounter)   = 0.1
           YVALUE(kounter,1) = Daa
           YVALUE(kounter,2) = Dab
           YVALUE(kounter,3) = Dac
           YVALUE(kounter,4) = Daa+Dab+Dac
           YVALUE(kounter,5) = DDIV

c           IF(FILED.EQ.1)THEN
c            WRITE(121,13)XVALUE(kounter),YVALUE(kounter,1),
c     >                 YVALUE(kounter,2),YVALUE(kounter,3),
c     >                 YVALUE(kounter,4),YVALUE(kounter,5)
c   13       FORMAT(6(1X,G11.4))
c           ENDIF

*          IF(YVALUE(kounter,4).gt.0.0)THEN
*            posvalue = 1
*          ELSE
*            negvalue = -1
*          ENDIF  
           
        END DO

*       WRITE (ANLY,'(''ROOTS(1-3)'',f8.4,f8.4,f8.4)') 
*     +        ROOTS(1), ROOTS(2), ROOTS(3)
*       WRITE (PLANE,'(''ROOTS(4-5)'',f8.4,f8.4)') 
*     +        ROOTS(4), ROOTS(5)

       CALL DRAW(XVALUE,WIDTHS,YVALUE,100,61,anly,5,99,-3.0,3.0,
     >           -HI,HI,IGNORS,0,AVS,NAVS,JOB,TITLE,XLAB,YLAB,ELABS,
     >           REF,NVIEW,PLANE,TABLE,1,6,1.0)

      ENDIF
C
C-----------------------------------------------------------------------
C     PLOT OF 1/C1 VERSUS CHI
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.582)THEN

        ELABS(1)  = 'C1p 1/C1plus '
        ELABS(2)  = 'C1m 1/C1minus'
        ELABS(3)  = 'C1p 1/C1plus '
        ELABS(4)  = 'C1m 1/C1minus'
        ELABS(5)  = 'C1p 1/C1plus '
        ELABS(6)  = 'C1m 1/C1minus'
        ELABS(7)  = 'C1p 1/C1plus '
        ELABS(8)  = 'C1m 1/C1minus'
        ELABS(9)  = 'C1p 1/C1plus '
        ELABS(10) = 'C1m 1/C1minus'

        XLAB = '    CHI    '
        YLAB = '   1/C1    '

        NVIEW = 'DETERMINING C1 FOR VARIOUS CHI'

        READ (GRAPH(45:46),'(I2)') FILED
     
        NPLOTS = NPLOTS + 1

        WRITE (REF,'(''Tbgrd '',F6.2,'',TB '',F6.2,'' ,State: '',I3)') 
     +        TBG,TB,ZI


        WRITE (IPLOT,9012) NPLOTS,REF

        CALL RZERO(YVALUE,100*MAXNGS)
        CALL RZERO(XVALUE,100)
        CALL RZERO(WIDTHS,100)

c        IF (FILED.EQ.1) THEN
c          WRITE(121,*)'582-C1 PROFILE'
c        ENDIF

        aa       = 1        
        icounter = 0

        DO icounter = 0,4
         jcounter = 1
         DO CHI = 0.0,10.0,0.1
          CHIpara = CHI*COS(icounter*3.1415927/4.0)
          CHIperp = CHI*SIN(icounter*3.1415927/4.0)
          Call Potentials(CHIpara,CHIperp,CHI,Phi00,Phi01,Phi10,
     >    Phi11,Phi20,Phi21,Psi00,Psi01,Psi02,Psi10,Psi11,Psi12,
     >    Psi20,Psi21,Psi22)
          Call COEFF_C1(Cplus,Cminus,CHIpara,CHIperp,CHI,
     >                  Psi01,Psi02,Psi10,Psi11,Psi12)
          XVALUE(jcounter)   = CHI
          WIDTHS(jcounter)   = 0.1
          YVALUE(jcounter,aa)= 1.0/Cplus
          YVALUE(jcounter,aa+1)= 1.0/Cminus
          jcounter=jcounter + 1
         END DO
         aa = aa + 2 
        END DO

c        IF(FILED.EQ.1)THEN
c         DO jcounter = 1,100
c          WRITE(121,14)XVALUE(jcounter),YVALUE(jcounter,1),
c     >     YVALUE(jcounter,2),YVALUE(jcounter,9),YVALUE(jcounter,10)
c   14     FORMAT(5(1X,G11.4))
c         END DO
c        ENDIF
    
       CALL DRAW(XVALUE,WIDTHS,YVALUE,100,100,anly,10,99,0.0,10.0,
     >           -HI,HI,IGNORS,0,AVS,NAVS,JOB,TITLE,XLAB,YLAB,ELABS,
     >           REF,NVIEW,PLANE,TABLE,1,6,1.0)

      ENDIF
C
C-----------------------------------------------------------------------
C     PLOT OF 1/C2 VERSUS CHI
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.583)THEN

        ELABS(1)  = 'C2p 1/C2plus '
        ELABS(2)  = 'C2m 1/C2minus'
        ELABS(3)  = 'C2p 1/C2plus '
        ELABS(4)  = 'C2m 1/C2minus'
        ELABS(5)  = 'C2p 1/C2plus '
        ELABS(6)  = 'C2m 1/C2minus'
        ELABS(7)  = 'C2p 1/C2plus '
        ELABS(8)  = 'C2m 1/C2minus'
        ELABS(9)  = 'C2p 1/C2plus '
        ELABS(10) = 'C2m 1/C2minus'

        XLAB = '    CHI    '
        YLAB = '   1/C2    '

        NVIEW = 'DETERMINING C2 FOR VARIOUS CHI'

        READ (GRAPH(45:46),'(I2)') FILED
     
        NPLOTS = NPLOTS + 1

        WRITE (REF,'(''Tbgrd '',F6.2,'',TB '',F6.2,'' ,State: '',I3)') 
     +        TBG,TB,ZI


        WRITE (IPLOT,9012) NPLOTS,REF

        CALL RZERO(YVALUE,100*MAXNGS)
        CALL RZERO(XVALUE,100)
        CALL RZERO(WIDTHS,100)

c        IF (FILED.EQ.1) THEN
c          WRITE(121,*)'583-C2 PROFILE'
c        ENDIF

        aa       = 1        
        icounter = 0

        DO icounter = 0,4
         jcounter = 1
         DO CHI = 0.0,10.0,0.1
          CHIpara = CHI*COS(icounter*3.1415927/4.0)
          CHIperp = CHI*SIN(icounter*3.1415927/4.0)
          Call Potentials(CHIpara,CHIperp,CHI,Phi00,Phi01,Phi10,
     >    Phi11,Phi20,Phi21,Psi00,Psi01,Psi02,Psi10,Psi11,Psi12,
     >    Psi20,Psi21,Psi22)
          Call COEFF_C2(Cplus,Cminus,CHIpara,CHIperp,CHI,
     >                  Psi01,Psi02,Psi20,Psi21,Psi22)
          XVALUE(jcounter)   = CHI
          WIDTHS(jcounter)   = 0.1
          YVALUE(jcounter,aa)= 1.0/Cplus
          YVALUE(jcounter,aa+1)= 1.0/Cminus
          jcounter=jcounter + 1
         END DO
         aa = aa + 2 
        END DO

c        IF(FILED.EQ.1)THEN
c         DO jcounter = 1,300
c          WRITE(121,16)XVALUE(jcounter),YVALUE(jcounter,1),
c     >     YVALUE(jcounter,2),YVALUE(jcounter,9),YVALUE(jcounter,10)
c   16     FORMAT(5(1X,G11.4))
c         END DO
c        ENDIF
    
       CALL DRAW(XVALUE,WIDTHS,YVALUE,100,100,anly,10,99,0.0,10.0,
     >           -HI,HI,IGNORS,0,AVS,NAVS,JOB,TITLE,XLAB,YLAB,ELABS,
     >           REF,NVIEW,PLANE,TABLE,1,6,1.0)

      ENDIF
C
C-----------------------------------------------------------------------
C     COMPONENTS OF FIG(R)
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.584) THEN

        ELABS(1) = 'FIG1COMPONENT1        FIG1'
        ELABS(2) = 'FIG2COMPONENT2        FIG2'
        ELABS(3) = 'FIGRNET THERMAL FORCE FIGR'

        XLAB = '   CHI   '
        YLAB = ' FORCE(N)'

        NVIEW = 'COMPONENTS OF FIG(R)'
     
        NPLOTS = NPLOTS + 1

        WRITE (REF,'(''Tbgrd '',F6.2,'',TB '',F6.2,'' ,State: '',I3)') 
     +        TBG,TB,ZI


        WRITE (IPLOT,9012) NPLOTS,REF

        CHIpara = 0.0
        KKpar   = 0.0
        kounter = 0
        jcounter = 0
        negvalue = 0
        posvalue = 0
        FORCE    = 3 

        CALL RZERO(YVALUE,100*MAXNGS)
        CALL RZERO(XVALUE,100)
        CALL RZERO(WIDTHS,100)

        DO kounter = 1,5
         ROOTS(kounter) = 99.9999
        END DO

        kounter = 0

c       IF (FILED.EQ.1) THEN
c         WRITE(121,*)'584-FIG(R) COMPONENTS'
c       ENDIF
 
        DO CHIpara = -3.0,3.1,0.1
           kounter = kounter + 1

            IF(negvalue.eq.-1.and.posvalue.eq.1)THEN
             jcounter = jcounter + 1
             ROOTS(jcounter)=BISCT(CHIpara-0.2,CHIpara-0.1,
     >                       VBG,TBG,TB,NB,MB,MI,ZI,ZB,XPER,
     >                       FORCE)
             negvalue=0
             posvalue=0
            ENDIF

           CALL Coulomb_Coll2(Kaa,Kab,Kaba,Kabb,Kac,KBetaf,
     >       Daa,Dab,Dac,Coulomb_log,Frictionf,CHIpara,Xper,
     >       VBG,TBG,TB,NB,MB,MI,ZI,ZB) 
        
           XVALUE(kounter)   = CHIpara
           WIDTHS(kounter)   = 0.1
           YVALUE(kounter,1) = Kaba * REAL(MI)*amu
           YVALUE(kounter,2) = Kabb * REAL(MI)*amu
           YVALUE(kounter,3) = Kab  * REAL(MI)*amu
        
c           IF(FILED.EQ.1)THEN
c            WRITE(121,21)XVALUE(kounter),YVALUE(kounter,1),
c     >                 YVALUE(kounter,2),YVALUE(kounter,3)
c   21       FORMAT(4(1X,G11.4))
c           ENDIF

          IF(YVALUE(kounter,5).gt.0.0)THEN
            posvalue = 1
          ELSE
            negvalue = -1
          ENDIF  
           
        END DO

       WRITE (ANLY,'(''ROOTS(1-3)'',f8.4,f8.4,f8.4)') 
     +        ROOTS(1), ROOTS(2), ROOTS(3)
       WRITE (PLANE,'(''ROOTS(4-5)'',f8.4,f8.4)') 
     +        ROOTS(4), ROOTS(5)

       CALL DRAW(XVALUE,WIDTHS,YVALUE,100,61,anly,3,99,-3.0,3.0,
     >           -HI,HI,IGNORS,0,AVS,NAVS,JOB,TITLE,XLAB,YLAB,ELABS,
     >           REF,NVIEW,PLANE,TABLE,1,6,1.0)

      ENDIF      
C
C-----------------------------------------------------------------------
C     AVERAGED ION VELOCITY PLOT ALONG THE FIELD LINE
C-----------------------------------------------------------------------
C 
      IF (IREF.EQ.585)THEN

       ELABS(1) = 'Vel VEL    '
       ELABS(2) = 'Alp ALPHA  '
       ELABS(3) = ' ChiCHIpara' 
    
       XLAB = ' S (M)'
       YLAB = ' VELOCITY (M/S)'

       NVIEW = 'AVERAGED ION VELOCITY ALONG B'

       READ (GRAPH(30:33),'(I4)') IZ
       READ (GRAPH(38:41),'(I4)') IR
       READ (GRAPH(45:46),'(I2)') FILED

c       IF (FILED.EQ.1) THEN
c         WRITE(121,*)'585-ION VELOCITY'
c       ENDIF

       NPLOTS = NPLOTS + 1

       WRITE (REF,'(''Ring: '',I3,'', K: '',F7.4,'' ,State: '',I3)') 
     +        IR,KKS(IR),IZ

       WRITE(IPLOT,9012) NPLOTS,REF

*      WRITE(iplot,)

       CALL RZERO (KVALS,MAXNKS*MAXNGS)

c       Modify the following to include values for S=0 ... even if this
c       is not part of the grid. Assume that if one end is not a grid
c       point then neither is the other.

        if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
           in  = 0
           inc = 0
        else
           in       = 1
           inc      = 2
           kouts(1) = 0.0

           KWIDS(1)         = kss(1,ir)
           kouts(nks(ir)+2) = ksmaxs(ir)
           kwids(nks(ir)+2) = ksmaxs(ir) - kss(nks(ir),ir)

        endif

c    !!!
c       Due to the fact that array Velavg is initially evaluated in units 
c       of distance and then converted back into units of velocity for the 
c       plotting routine, the values of the array become prone to round of 
c       error because they are REAL valued instead of DOUBLE PRECISION.  
c       As a result, when many impurity particles are released (1000) with 
c       a small time step, QTIM, (5e-8 s) the average value of the array for 
c       those grid cells where many iterations occur will be substantially 
c       different than the true value and will produce a misleading profile.  
c       Note that the actual ion motions are correctly calculated, it is only 
c       the generated profile that is affected.  This is most evident when 
c       using the FIG of the fluid approximation model, since its value is 
c       always the same for any grid cell, yet the force profile will be 
c       different for different QTIM's and total impurity ion numbers.  
c    !!!  

       DO 1271 IK = 1, NKS(IR)

         KOUTS(IK+in) = KSS(IK,IR)
         KWIDS(IK+in) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
          
         KVALS(IK+in,1) = VELavg(IK,IR,IZ)
         KVALS(IK+in,2) = SQRT(CRMB/(2*KTIBS(IK,IR)*EMI))
         KVALS(IK+in,3) = (KVALS(IK+in,1)-(KVHS(IK,IR)/QTIM))
     +                    * KVALS(IK+in,2)

*         WRITE(0,*)(KVHS(IK,IR)/QTIM),VELavg(IK,IR,IZ),KVALS(IK+in,3)
*         WRITE(0,*)'VELavg',VELavg(IK,IR,IZ)
*         WRITE(0,*)'CHI', KVALS(IK+in,3)

         WRITE(iplot,'(2I4,8G12.4)') IK,IR,KVALS(IK+IN,1),
     +      KVALS(IK+IN,2),
     +      KVALS(IK+IN,3),KVALS(IK+IN,4),KVALS(IK+IN,5),KSS(IK,IR),
     +      RS(IK,IR),ZS(IK,IR)

c          IF(FILED.EQ.1)THEN
c            WRITE(121,17)KOUTS(IK+in),KVALS(IK+in,1),
c     >                   KVALS(IK+in,2),KVALS(IK+in,3)
c   17       FORMAT(4(1X,G11.5))
c          ENDIF

 1271  CONTINUE

c For S=0 meters (at target plate), the Y-axis values will be taken 
c as those calculated for the first grid cell: 
                
          KVALS(1,1) = KVALS(1+in,1)
          KVALS(1,2) = KVALS(1+in,2)
          KVALS(1,3) = KVALS(1+in,3) 

c For S=Smax (at the other target plate), the Y-axis values will be
c taken as those calculated for the last grid cell: 

          KVALS(nks(ir)+inc,1) = KVALS(nks(ir)+in,1)
          KVALS(nks(ir)+inc,2) = KVALS(nks(ir)+in,2)
          KVALS(nks(ir)+inc,3) = KVALS(nks(ir)+in,3)

       CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >   3,99,KOUTS(1),KOUTS(NKS(IR)+inc),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >   JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0)
c
        enldist = mgst*kouts(nks(ir)+inc)
c
        IF (clsup.eq.1) THEN
          CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >      3,99,KOUTS(1),enldist   ,-HI,HI,IGNORS,ITEC,AVS,
     >      NAVS,
     >      JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0)

          CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >      3,99,KOUTS(NKS(IR)+inc)-enldist,KOUTS(NKS(IR)+inc),
     >      -HI,HI,IGNORS,ITEC,AVS,NAVS,
     >      JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0)
        ENDIF
  
      ENDIF

c
c psmod
c


c
c     Fortmat statements
c

 9012 FORMAT(1X,'PLOT',I3,4X,A)

      return
      end

