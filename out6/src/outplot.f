c     -*-Fortran-*-
c
      SUBROUTINE RVALUE (CVALS,VS,II,NIIS,MAXIIS,FT,FP,
     >                   IXMIN,IXMAX,IYMIN,IYMAX,VMIN,VMAX)
      use mod_params
      use mod_cgeom
      implicit none
c      IMPLICIT LOGICAL (A-Z)
C     INCLUDE "PARAMS"
c     include 'params'
C     INCLUDE "CGEOM"
c     include 'cgeom'
      INTEGER II,NIIS,MAXIIS,IXMIN,IXMAX,IYMIN,IYMAX,IX,IY,JJ,IK,IR
      REAL VS(MAXNKS,MAXNRS,MAXIIS)
      REAL CVALS(MAXGXS,MAXGYS),FT,FP,VMIN,VMAX
C
C  *********************************************************************
C  *                                                                   *
C  *  RVALUE:  MULTIPURPOSE ROUTINE TO EXTRACT CVALS FOR PLOTTING.     *
C  *                                                                   *
C  *            CHRIS FARRELL  (HUNTERSKIL)  FEBRUARY 1989             *
C  *                                                                   *
C  *********************************************************************
C
      logical griderr
      real r,z
c
      VMIN = HI
      VMAX =-HI
      CALL RZERO (CVALS, MAXGXS*MAXGYS)
C
      DO 130 IX = IXMIN, IXMAX
        R = (RMAX-RMIN) * REAL(IX)/REAL(NXS) + RMIN - 0.5 * DR
        DO 120 IY = IYMIN, IYMAX
          Z = (ZMAX-ZMIN) * REAL(IY)/REAL(NYS) + ZMIN - 0.5 * DZ
          call gridpos(ik,ir,r,z,.false.,griderr)
          if (griderr) goto 120
c
c          IK = IKXYS(IX,IY)
c          IR = IRXYS(IX,IY)
c
c          IF (IFXYS(IX,IY).EQ.0) GOTO 120
c

          IF (II.EQ.NIIS+1) THEN
            DO 110 JJ = 3, NIIS
              CVALS(IX,IY) = CVALS(IX,IY) + VS(IK,IR,JJ)
  110       CONTINUE
          ELSEIF (II.EQ.0) THEN
            CVALS(IX,IY) = FT * VS(IK,IR,1+PSHIFT) - FP * VS(IK,IR,1)
          ELSE
            CVALS(IX,IY) = VS(IK,IR,II)
          ENDIF
          VMIN = MIN (VMIN, CVALS(IX,IY))
          VMAX = MAX (VMAX, CVALS(IX,IY))
  120   CONTINUE
  130 CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE SUPIMP (OPTION)
      use mod_params
      use mod_cgeom
      use mod_comtor
      IMPLICIT none
      CHARACTER*(*) OPTION
C
C  *********************************************************************
C  *                                                                   *
C  *  SUPIMP:  ROUTINE TO SUPERIMPOSE A CONTOUR RING ON THE PLOT.      *
C  *                                                                   *
C  *            CHRIS FARRELL  (HUNTERSKIL)  MARCH 1989                *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE "PARAMS"
c     include 'params'
C     INCLUDE "CGEOM"
c     include 'cgeom'
C     INCLUDE "COMTOR"
c     include 'comtor'
C     IPP/01 - Krieger: fixed initialization in declaration statement
C     by adding appropiate data statement (picky SUN compiler)
      CHARACTER*36 NAME
      REAL   DVALS(MAXNDS,MAXNGS)
      INTEGER IR,ID
      DATA NAME /' '/
C
C
c     If a JET grid is in use - invoke the superposition routine
c     for JET grids - placing this decision here is more efficient
c     than testing at every place that this routine can be called
c     from.
c
c slmod begin
      if (cgridopt.eq.0.or.cgridopt.eq.3.or.cgridopt.eq.6.or.
     .    cgridopt.eq.8) then
c
c      if (cgridopt.eq.0.or.cgridopt.eq.3) then
c slmod end
         call supimp2(option)
         return
      endif
c
      IF (OPTION.EQ.'FULL') THEN
       DO 100 IR = 1, NRS
         CALL GRTRAC (RS(1,IR),ZS(1,IR),NKS(IR),NAME,'LINE',-1)
         CALL GRTRAC (RS(1,IR),ZS(1,IR),NKS(IR),NAME,'POINT',-1)
  100  CONTINUE
      ELSE
       CALL GRTRAC(RS(1,1     ),ZS(1,1     ),NKS(1     ),NAME,'LINE',-1)
       CALL GRTRAC(RS(1,IRSEP ),ZS(1,IRSEP ),NKS(IRSEP ),NAME,'LINE',-1)
       CALL GRTRAC(RS(1,IRWALL),ZS(1,IRWALL),NKS(IRWALL),NAME,'LINE',-1)
       CALL GRTRAC(RS(1,IRTRAP),ZS(1,IRTRAP),NKS(IRTRAP),NAME,'LINE',-1)
       if (cgridopt.eq.2) then
          CALL GRTRAC (RS(1,IRsep2),ZS(1,IRsep2),NKS(IRsep2),
     >                 NAME,'LINE',-1)
          CALL GRTRAC (RS(1,IRWALL2),ZS(1,IRWALL2),NKS(IRWALL2),
     >                 NAME,'LINE',-1)
          CALL GRTRAC (RS(1,IRTRAP2),ZS(1,IRTRAP2),NKS(IRTRAP2),
     >                 NAME,'LINE',-1)
          write (6,*) 'ir:',irsep2,irwall2,irtrap2,nks(irtrap2)
       endif
      ENDIF
C
      CALL GRTRAC (R0 ,Z0 ,1,NAME,'POINT',-1)
      CALL GRTRAC (RXP,ZXP,1,NAME,'POINT',-1)
C
      DO 120 ID = 1, NDS
        DVALS(ID,1) = RS(IKDS(ID),IRDS(ID))
        DVALS(ID,2) = ZS(IKDS(ID),IRDS(ID))
  120 CONTINUE
      CALL GRTRAC (DVALS(1,1),DVALS(1,2),NDSIN,NAME,'LINE',-1)
      if (cgridopt.eq.2) then
         CALL GRTRAC (DVALS(NDSIN+1,1),DVALS(NDSIN+1,2),
     >                            NDSin2-NDSIN,NAME,'LINE',-1)
         CALL GRTRAC (DVALS(NDSIN2+1,1),DVALS(NDSIN2+1,2),
     >                           NDSin3-NDSIN2,NAME,'LINE',-1)
         CALL GRTRAC (DVALS(NDSIN3+1,1),DVALS(NDSIN3+1,2),
     >                              NDS-NDSIN3,NAME,'LINE',-1)
      else
         CALL GRTRAC (DVALS(NDSIN+1,1),DVALS(NDSIN+1,2),
     >                               NDS-NDSIN,NAME,'LINE',-1)
      endif
C
      RETURN
      END
C
C
C
      SUBROUTINE SUPIMPOLD (OPTION)
      use mod_params
      use mod_cgeom
      use mod_comtor
      IMPLICIT none
      CHARACTER*(*) OPTION
C
C  *********************************************************************
C  *                                                                   *
C  *  SUPIMP:  ROUTINE TO SUPERIMPOSE A CONTOUR RING ON THE PLOT.      *
C  *                                                                   *
C  *            CHRIS FARRELL  (HUNTERSKIL)  MARCH 1989                *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE "PARAMS"
c     include 'params'
C     INCLUDE "CGEOM"
c     include 'cgeom'
C     INCLUDE "COMTOR"
c     include 'comtor'
C     IPP/01 - Krieger: fixed initialization in declaration statement
C     by adding appropiate data statement (picky SUN compiler)
      CHARACTER*36 NAME
      REAL   DVALS(MAXNDS,MAXNGS)
      INTEGER IR,ID
c
      DATA NAME /' '/
C
      IF (OPTION.EQ.'FULL') THEN
       DO 100 IR = 1, NRS
         CALL GRTRAC (RS(1,IR),ZS(1,IR),NKS(IR),NAME,'LINE',-1)
         CALL GRTRAC (RS(1,IR),ZS(1,IR),NKS(IR),NAME,'POINT',-1)
  100  CONTINUE
      ELSE
       CALL GRTRAC(RS(1,1     ),ZS(1,1     ),NKS(1     ),NAME,'LINE',-1)
       CALL GRTRAC(RS(1,IRSEP ),ZS(1,IRSEP ),NKS(IRSEP ),NAME,'LINE',-1)
       CALL GRTRAC(RS(1,IRWALL),ZS(1,IRWALL),NKS(IRWALL),NAME,'LINE',-1)
       CALL GRTRAC(RS(1,IRTRAP),ZS(1,IRTRAP),NKS(IRTRAP),NAME,'LINE',-1)
       if (cgridopt.eq.2) then
          CALL GRTRAC (RS(1,IRsep2),ZS(1,IRsep2),NKS(IRsep2),
     >                 NAME,'LINE',-1)
          CALL GRTRAC (RS(1,IRWALL2),ZS(1,IRWALL2),NKS(IRWALL2),
     >                 NAME,'LINE',-1)
          CALL GRTRAC (RS(1,IRTRAP2),ZS(1,IRTRAP2),NKS(IRTRAP2),
     >                 NAME,'LINE',-1)
          write (6,*) 'ir:',irsep2,irwall2,irtrap2,nks(irtrap2)
       endif
      ENDIF
C
      CALL GRTRAC (R0 ,Z0 ,1,NAME,'POINT',-1)
      CALL GRTRAC (RXP,ZXP,1,NAME,'POINT',-1)
C
      DO 120 ID = 1, NDS
        DVALS(ID,1) = RS(IKDS(ID),IRDS(ID))
        DVALS(ID,2) = ZS(IKDS(ID),IRDS(ID))
  120 CONTINUE
      CALL GRTRAC (DVALS(1,1),DVALS(1,2),NDSIN,NAME,'LINE',-1)
      if (cgridopt.eq.2) then
         CALL GRTRAC (DVALS(NDSIN+1,1),DVALS(NDSIN+1,2),
     >                            NDSin2-NDSIN,NAME,'LINE',-1)
         CALL GRTRAC (DVALS(NDSIN2+1,1),DVALS(NDSIN2+1,2),
     >                           NDSin3-NDSIN2,NAME,'LINE',-1)
         CALL GRTRAC (DVALS(NDSIN3+1,1),DVALS(NDSIN3+1,2),
     >                              NDS-NDSIN3,NAME,'LINE',-1)
      else
         CALL GRTRAC (DVALS(NDSIN+1,1),DVALS(NDSIN+1,2),
     >                               NDS-NDSIN,NAME,'LINE',-1)
      endif
C
      RETURN
      END
C
C
C
      SUBROUTINE INTLOS (TVALS,TOUTS,TWIDS,ROBS,ZOBS,DRAD,NUMTHE,
     >           AVPTS,ATYPE,THEMIN,THEMAX,DTHE,IIMAX,IZMIN,IZMAX,
     >           VS,MAXIIS,
     >           WS,WS2,WOPT,ANLY,PSWITCH,PIZS,
     >           FT,FP)
C
      use mod_params
      use mod_cgeom
      implicit none
C     INCLUDE "PARAMS"
c     include 'params'
C     INCLUDE "CGEOM"
c     include 'cgeom'
      LOGICAL PSWITCH
      INTEGER MAXIIS,IZMIN,IZMAX,NUMTHE,AVPTS,IIMAX,WOPT,ATYPE
      INTEGER PIZS(MAXIIS)
      REAL DTHE,THEMIN,DRAD,THEMAX
      REAL VS(MAXNKS,MAXNRS,MAXIIS)
      REAL WS(MAXNKS,MAXNRS,MAXIIS)
      REAL WS2(MAXNKS,MAXNRS)
      REAL TVALS(MAXTHE,MAXIIS)
      REAL TOUTS(MAXTHE)
      REAL TWIDS(MAXTHE)
      REAL ROBS,ZOBS,FT,FP
      CHARACTER*36 ANLY
C
C  *********************************************************************
C  *                                                                   *
C  *  INTLOS:  THIS ROUTINE INTEGRATES ALONG LINES OF SIGHT FROM A     *
C  *           GIVEN OBSERVATION POSITION. IT DIVIDES THE THEMIN/THEMAX*
C  *           ARC FACING THE PLASMA INTO A NUMBER OF EQUALLY SPACED   *
C  *           ARCS AND THEN INTEGRATES CONTRIBUTIONS ALONG EACH OF    *
C  *           THE LINES OF SIGHT.                                     *
C  *           THE CODE IS IMPLEMENTED SO THAT THEMIN AND DTHETA       *
C  *           ARE SPECIFIED IN DEGREES COUNTER-CLOCKWISE FROM THE     *
C  *           POSITIVE R-AXIS.                                        *
C  *                                                                   *
C  *           DAVID ELDER   (TORONTO)       AUGUST 1991               *
C  *                                                                   *
C  *********************************************************************
C
      REAL RBASE,ZBASE,DIST,THETA,DTHE2,MFACT
      REAL SINTHE,COSTHE,D,R,Z,TMPTOT,TMPTOT2
      REAL WVALS(MAXTHE,MAXNGS)
      INTEGER I,IZ,IT,IK,IR,COUNT,IX,IY,J,IA,IN,IZ1
      logical griderr
C
      ik = 0
      ir = 0
c
      IF ((NUMTHE*AVPTS).GT.MAXTHE) THEN
         AVPTS = 1
         WRITE(6,*) 'NUMBER OF AVERAGING BINS EXCEEDS ALLOWABLE'
      ENDIF
C
      CALL RZERO(TVALS,MAXTHE*MAXIIS)
      CALL RZERO(WVALS,MAXTHE*MAXNGS)
      CALL RZERO(TOUTS,MAXTHE)
      CALL RZERO(TWIDS,MAXTHE)
C
      RBASE = R0-ROBS
      ZBASE = Z0-ZOBS
C
      THEMIN = THEMIN * DEGRAD
      DTHE   = DTHE * DEGRAD
      THEMAX = THEMIN + (NUMTHE-1) * DTHE
      DTHE2  = DTHE/AVPTS
      DIST   = SQRT(RBASE**2+ZBASE**2)
C
      DO 400 IT = 1,NUMTHE
       DO 400 IA = 1,AVPTS
         THETA = THEMIN + (IT-1)*DTHE
     >           - 0.5*DTHE +0.5 * DTHE2 + (IA-1)*DTHE2
         IN = (IT-1) * AVPTS + IA
         SINTHE = SIN(THETA)
         COSTHE = COS(THETA)
C
C         SINTHE = SIN(THETA+PSI0)
C         COSTHE = COS(THETA+PSI0)
C
         DO 400 IZ = IZMIN,IZMAX
            IF (PSWITCH) THEN
              IZ1 = PIZS(IZ) + 2
            ELSE
              IZ1 = IZ
            ENDIF
C
C            WRITE(6,*) 'PLRP TEST:',PSWITCH,IZ1,PIZS(IZ),IZ,
C     >                  PSHIFT,PNCNT
C

            COUNT = 0
100         CONTINUE
            COUNT = COUNT +1
            D = COUNT * DRAD
            R = ROBS + D * COSTHE
            Z = ZOBS + D * SINTHE
C
C           THIS CUT OFF FOR D IS ARBITRARY BUT SHOULD BE
C           BEYOND THE PLASMA BOUNDARIES SINCE THE BASELINE DISTANCE
C           IS DETERMINED BY THE DISTANCE BETWEEN THE OBSERVATION
C           POSITION AND THE PLASMA CENTRE AND ONLY AFTER COUNTING THIS
C           DISTANCE CAN A LINE BE DISCARDED.
C
            IF (D.GT.DIST.AND.
     >         (Z.GT.ZMAX.OR.Z.LT.ZMIN.OR.R.GT.RMAX.OR.R.LT.RMIN))
     >         GOTO 300
c
             call gridpos(ik,ir,r,z,.false.,griderr)
             if (griderr) goto 100
c
c            IX = MAX(1,MIN(NXS,INT((R-RMIN)/DR)+1))
c            IY = MAX(1,MIN(NYS,INT((Z-ZMIN)/DZ)+1))
c            IF (IFXYS(IX,IY).EQ.0) GOTO 100
C
c            IK = IKXYS(IX,IY)
c            IR = IRXYS(IX,IY)
c
            IF (IZ.EQ.IIMAX) THEN
               DO 200 I = 3,IIMAX-1
                  TVALS(IN,IZ-IZMIN+1) = TVALS(IN,IZ-IZMIN+1)
     >                          + VS(IK,IR,I)
                  IF (WOPT.EQ.1) THEN
                     WVALS(IN,IZ-IZMIN+1) = WVALS(IN,IZ-IZMIN+1)
     >                        + WS(IK,IR,I)*VS(IK,IR,I)
                  ELSEIF (WOPT.EQ.2) THEN
                     WVALS(IN,IZ-IZMIN+1) = WVALS(IN,IZ-IZMIN+1)
     >                        + WS2(IK,IR)*VS(IK,IR,I)
                  ENDIF
200            CONTINUE
            ELSEIF (IZ.EQ.0) THEN
               TVALS(IN,IZ-IZMIN+1) = TVALS(IN,IZ-IZMIN+1)
     >                + FT * VS(IK,IR,2) - FP * VS(IK,IR,1)
                  IF (WOPT.EQ.1) THEN
                    WVALS(IN,IZ-IZMIN+1) = WVALS(IN,IZ-IZMIN+1)
     >                      + FT*WS(IK,IR,2)*VS(IK,IR,1+PSHIFT)
     >                      - FP*WS(IK,IR,1)*VS(IK,IR,1)
                  ELSEIF (WOPT.EQ.2) THEN
                    WVALS(IN,IZ-IZMIN+1) = WVALS(IN,IZ-IZMIN+1)
     >                      + FT*WS2(IK,IR)*VS(IK,IR,1+PSHIFT)
     >                      - FP*WS2(IK,IR)*VS(IK,IR,1)
                  ENDIF
            ELSE
               TVALS(IN,IZ-IZMIN+1) = TVALS(IN,IZ-IZMIN+1)
     >                          + VS(IK,IR,IZ)
               IF (WOPT.EQ.1) THEN
                 WVALS(IN,IZ-IZMIN+1) = WVALS(IN,IZ-IZMIN+1)
     >                     + WS(IK,IR,IZ1)*VS(IK,IR,IZ)
               ELSEIF (WOPT.EQ.2) THEN
                 WVALS(IN,IZ-IZMIN+1) = WVALS(IN,IZ-IZMIN+1)
     >                     + WS2(IK,IR)*VS(IK,IR,IZ)
               ENDIF
            ENDIF
            GOTO 100
 300     CONTINUE
C
         TVALS(IN,IZ-IZMIN+1)=TVALS(IN,IZ-IZMIN+1) * DRAD
         IF (WOPT.EQ.1.OR.WOPT.EQ.2)
     >       WVALS(IN,IZ-IZMIN+1)=WVALS(IN,IZ-IZMIN+1)*DRAD
C

 400  CONTINUE
C
C
      DO 500 IT = 1, NUMTHE
        THETA = THEMIN + (IT-1)*DTHE
        TOUTS(IT) = THETA * RADDEG
C
C       LEAVE TWIDS IN RADIANS FOR NOW
C
        TWIDS(IT) = DTHE
        DO 550 IZ = IZMIN,IZMAX
          TMPTOT = 0.0
          TMPTOT2 = 0.0
          DO 600 IA = 1,AVPTS
            TMPTOT = TMPTOT + TVALS((IT-1)*AVPTS+IA,IZ-IZMIN+1)
            IF (WOPT.EQ.1.OR.WOPT.EQ.2)
     >         TMPTOT2 = TMPTOT2 + WVALS((IT-1)*AVPTS+IA,IZ-IZMIN+1)
600       CONTINUE
          IF (WOPT.EQ.1.OR.WOPT.EQ.2) THEN
            IF (TMPTOT.LE.0.0) THEN
              TVALS(IT,IZ-IZMIN+1) = 0.0
            ELSE
              TVALS(IT,IZ-IZMIN+1) = TMPTOT2/TMPTOT
            ENDIF
          ELSE
            TVALS(IT,IZ-IZMIN+1) = TMPTOT/AVPTS
          ENDIF
550     CONTINUE
500   CONTINUE

C
C      APPLY THE MULTIPLICATIVE NORMALIZATION SPECIFIED BY ATYPE
C
C      0  -  NONE
C      1  -  DTHE / ( 2 * PI )
C      2  -  1 / ( 2 * PI )
C      3  -  1 / (4 * PI )
C
C
      IF (ATYPE.EQ.0) THEN
        WRITE(ANLY,'(''NO SCALE FACTOR APPLIED'')')
      ELSEIF (ATYPE.NE.0) THEN
        IF (ATYPE.EQ.1) THEN
          MFACT = DTHE / (2.0 * PI)
          WRITE(ANLY,'(''SCALE FACTOR = '',G12.6,'' / (2*PI)'')')
     >               DTHE
        ELSEIF (ATYPE.EQ.2) THEN
          MFACT = 1.0 / (2.0 * PI)
          WRITE(ANLY,'(''SCALE FACTOR = 1 / (2*PI)'') ')
        ELSEIF (ATYPE.EQ.3) THEN
          MFACT = 1.0 / ( 4.0 * PI)
          WRITE(ANLY,'(''SCALE FACTOR = 1 / (4*PI)'')')
        ENDIF
C
        DO 700 IT = 1 , NUMTHE
          DO 700 IZ = IZMIN,IZMAX
            TVALS(IT,IZ-IZMIN+1) = TVALS(IT,IZ-IZMIN+1) * MFACT
 700    CONTINUE
C
      ENDIF
C
      THEMIN = (THEMIN-DTHE/2.0) * RADDEG
      THEMAX = (THEMAX+DTHE/2.0) * RADDEG
      DTHE   = DTHE * RADDEG

C      WRITE(6,*) 'AVPTS',AVPTS
C
C      WRITE(6,*) 'CENTRE:',R0,Z0
C      WRITE(6,*) 'NUMTHE:',NUMTHE
C      WRITE(6,*) 'THEMAX/MIN:', THEMAX, THEMIN,RADDEG,
C     >           THEMAX*DEGRAD,
C     >           THEMIN*DEGRAD,DTHE,DTHE*DEGRAD
C      WRITE(6,*) 'DRAD,IIMAX,IZMIN,IZMAX:',DRAD,IIMAX,IZMIN,
C     >           IZMAX
C      WRITE(6,*) 'THETAS:',RBASE,ZBASE,DIST,PSI0
C      WRITE(6,*) 'DRAD:', DRAD
C      WRITE(6,*) 'MAX/MIN:',RMIN,RMAX,ZMIN,ZMAX
C      WRITE(6,*) 'TOUTS:',(TOUTS(I),I=1,NUMTHE)
C      WRITE(6,*) 'TWIDS:',(TWIDS(I),I=1,NUMTHE)
C      DO 1000 J = 1,IZMAX-IZMIN+1
C        WRITE(6,*) 'TVALS:',J,':',(TVALS(I,J),I=1,NUMTHE)
C1000  CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE CONTOUR (ICNTR,NGS,VS,II,NIIS,MAXIIS,FT,FP,MFACT,
     >                    XOUTS,IXMIN,IXMAX,YOUTS,IYMIN,IYMAX,
     >                    XXMIN,XXMAX,YYMIN,YYMAX,nconts,conts,
     >                    cntropt,minscale,maxscale)
      use mod_params
      use mod_cgeom
      use mod_comgra
      use mod_plot_switches
      use mod_slout
      IMPLICIT NONE
C     INCLUDE "PARAMS"
c     include 'params'
C     INCLUDE "CGEOM"
c     include 'cgeom'
c
c     need this common block here to set iplots to 0 for
c     superposition of separatrix (moving this plot routine
c     to the end of subr. contour produces dashed lines
c     because iplots is incremented by the contour/false color
c     routines). Krieger IPP/97
c
c     include 'comgra' 
c     include 'plot_switches'
c
c slmod begin - new
c     INCLUDE 'slout'

      INTEGER i,j
c slmod end
c
      INTEGER ICNTR,NGS,II,NIIS,MAXIIS,IXMIN,IXMAX,IYMIN,IYMAX
      integer nconts,cntropt
      real conts(maxpts)
      REAL VS(MAXNKS,MAXNRS,MAXIIS),FT,FP,MFACT
      REAL XOUTS(MAXGXS),YOUTS(MAXGYS)
      REAL XXMIN,XXMAX,YYMIN,YYMAX
      real minscale,maxscale
C
C  *********************************************************************
C  *                                                                   *
C  *  CONTOUR: IF ICNTR = 0, GENERATES A CONTOUR PLOT BY MAPPING       *
C  *           THE INPUT ARRAY VS ONTO THE XY GRID.  IF ICNTR = 1      *
C  *           A FALSE COLOUR PLOT OF THE SAME DATA IS GENERATED.      *
C  *                                                                   *
C  *            LORNE HORTON   (JET)         JULY     1993             *
C  *                                                                   *
C  *********************************************************************
C
      INTEGER IG,tmpngs,increment,npos,nneg
      REAL    VMIN,VMAX,VLO,VHI
      REAL    CVALXY(MAXGXS,MAXGYS), CVALKR(MAXNKS,MAXNRS)
      real    tmpconts(maxpts)
      CHARACTER*36 NAME
      integer cntr_order,start_ngs,end_ngs,step_ngs
c
      real logtable(21)
      data logtable /0.0e+0, 0.5e-6, 1.0e-6,
     >               0.2e-5, 0.5e-5, 1.0e-5,
     >               0.2e-4, 0.5e-4, 1.0e-4,
     >               0.2e-3, 0.5e-3, 1.0e-3,
     >               0.2e-2, 0.5e-2, 1.0e-2,
     >               0.2e-1, 0.5e-1, 1.0e-1,
     >               0.2e+0, 0.5e+0, 1.0e+0/
C
      IF (ICNTR.EQ.0) THEN
        CALL RVALXY(CVALXY,VS,II,NIIS,MAXIIS,FT,FP,MFACT,
     >              IXMIN,IXMAX,IYMIN,IYMAX,VMIN,VMAX)
      ELSE
        CALL RVALKR(CVALKR,VS,II,NIIS,MAXIIS,FT,FP,MFACT,
     >              XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
      ENDIF

      cntr_order = 1
      ! reverse contour order if less than zero
      if (cntropt.lt.0) then 
         cntropt = abs(cntropt)
         cntr_order= -1
      endif

c      write(0,'(a,8(1x,g18.8))') 'Contour:',xxmin,xxmax,yymin,yymax,
c     >              vmin,vmax,minscale,maxscale

c
c slmod begin
      if (minscale.ne.maxscale) then
        vmin=minscale
        vmax=maxscale
      elseif (minscale.lt.0.0.and.maxscale.lt.0.0) then
        vmin = -minscale * vmax
      endif
c
c      if (minscale.ne.maxscale) then
c        vmin=minscale
c        vmax=maxscale
c      endif
c slmod end
c
c     SUPIMP moved to end; in case of false color plots, the separatrix
c     was hidden by the coloring procedures - Krieger IPP/97
c
c     CALL SUPIMP('SELECT')
c
c     Set up the contour information.
c
c
c     Original default option - quadratic
c
      if (cntropt.eq.0) then
         tmpngs = ngs
         do ig = 1,tmpngs
            tmpconts(ig) = REAL(IG)/REAL(NGS)/REAL(NGS+1-IG)
         end do
c
c     10% levels - note 0.0 is excluded since in most cases this will 
c                  result in the rest of the entire grid being coloured.
c
      elseif (cntropt.eq.1) then
         tmpngs = 9
         do ig = 1,tmpngs
            tmpconts(ig) = ig * 0.1
         end do
c
c     Logarithmic
c
      elseif (cntropt.eq.2) then
         tmpngs = ngs
         do ig = 1,tmpngs
            tmpconts(ig) = real(2**ig)/real(2**ngs)
         end do
c
c     User specified
c
      elseif (cntropt.eq.3) then
c
         tmpngs = min(ngs, nconts)
         increment = nconts-tmpngs
         do ig = 1,tmpngs
            tmpconts(ig) = conts(increment+ig)
         end do  
c
c         tmpngs = nconts
c         do ig = 1,tmpngs
c            tmpconts(ig) = conts(ig)
c         end do
c
c     User specified - Absolute values
c
      elseif (cntropt.eq.4) then
c
         tmpngs = nconts
         do ig = 1,tmpngs
            tmpconts(ig) = conts(ig)
         end do
c
c     B2/Eirene Velocity plots 
c
      elseif (cntropt.eq.5) then 
c
         vmax=max(abs(vmin),abs(vmax))
         vmin=-vmax
c
         nneg=20
         npos=20
c
         do ig=1,nneg
           tmpconts(ig)=vmin*logtable(22-ig)
         enddo
         do ig=1,npos+1
           tmpconts(nneg+ig)=vmax*logtable(ig)
         enddo
c
         tmpngs=nneg+npos
c
c
c     10% levels - including 0.0
c                  
c
      elseif (cntropt.eq.6) then
         tmpngs = 10
         do ig = 1,tmpngs
            tmpconts(ig) = (ig-1) * 0.1
         end do
      endif
c
c     Debug
c
      write (6,*) 'Contour:',vmin,vmax,tmpngs,ngs,nconts
c
c     Set maximum value
c
      if (cntropt.eq.4.or.cntropt.eq.5) then
         tmpconts(tmpngs+1) = vmax
      else
         tmpconts(tmpngs+1) = 1.0
      endif
c
c...dev
c jdemod - put contour labels with highest values at top
c
c      DO IG = 1,tmpNGS
c
      if (cntr_order.eq.1) then 
         start_ngs = 1
         end_ngs = tmpngs
         step_ngs = 1
      elseif (cntr_order.eq.-1) then 
         start_ngs = tmpngs
         end_ngs = 1
         step_ngs = -1
      endif

      DO IG = start_ngs,end_ngs,step_ngs
c
        if (cntropt.eq.4.or.cntropt.eq.5) then
c
           VLO = tmpconts(ig)
           VHI = tmpconts(ig+1)
c
        else

           VLO = VMIN + tmpconts(ig)  *(VMAX-VMIN)
           VHI = VMIN + tmpconts(ig+1)*(VMAX-VMIN)
c
        endif
c
c       Plots with more than 20 contour levels
c
        if (tmpngs.gt.20) then 
c
           WRITE (NAME,'(4X,1P,E8.1)') VLO
           IF (IG.EQ.tmpNGS) THEN
              NAME(13:27) = '(MX= 0.0e-00)'
              write(name(17:24),'(1p,e8.1)') vhi
           endif
c
        else 
c
           WRITE (NAME,'(4X,1P,E8.1,1x,''TO'',1x,E8.1)') VLO,VHI
           IF (IG.EQ.TMPNGS) NAME(25:36) = ' (MAX)'
c
        endif
c
c       For first contour drawing - we want to INCLUDE the lower boundary
c       in the plot. Otherwise it should be exlcuded.
c 
        if (ig.eq.1) then  
           first_contour = .true.
        else
           first_contour = .false.
        endif
c
        IF (ICNTR.EQ.0) THEN
          CALL GRCONT (CVALXY,IXMIN,IXMAX,MAXGXS,IYMIN,IYMAX,MAXGYS,
     >                 VLO,XOUTS,YOUTS,NAME)
        ELSE
          CALL GRCOLR (CVALKR,VLO,VHI,NAME)
        ENDIF
      ENDDO
c
c     Reset the first_contour switch - so it won't affect any other routines
c
      first_contour = .false. 
c
c slmod begin - new
c...  Print comments:
      CALL PSPACE (0.0, 1.35, 0.0,1.0)
      CALL CSPACE (0.0, 1.35, 0.0,1.0)
      CALL MAP    (0.0, 1.35, 0.0,1.0)
      CALL CTRMAG (12)
      DO i = 1, 10
        IF (char(i).NE.' ') CALL PLOTST(1.00,0.590+(i-1)*0.02,char(i))
      ENDDO
      DO i = 20, 30
        j = i - 19
        IF (char(i).NE.' ') CALL PLOTST(1.00,0.550-(j-1)*0.02,char(i))
      ENDDO

c...  For toroidal camera plot 972:
      IF (slopt.EQ.3) THEN
        iplots = 0
        RETURN
      ENDIF
c slmod end

c
c`    Reset cntropt sign to value on entry
c
      cntropt = cntropt *  cntr_order
c
c     moved to end; in case of false color plots, the separatrix
c     was hidden by the coloring procedures - Krieger IPP/97
c
      iplots = 0
c slmod begin
      CALL FULL
c slmod end
      CALL SUPIMP('SELECT')
c
      CALL FRAME
C
      RETURN
      END
C
C
C
      SUBROUTINE CONTOURXY(ICNTR,NGS,VS,
     >                    XXMIN,XXMAX,YYMIN,YYMAX,nconts,conts,
     >                    cntropt,minscale,maxscale,
     >                    maxix,maxiy,nix,niy,raxis,zaxis,
     >                    overlay_grid)
      use mod_params
      use mod_comgra
      use mod_slout
      IMPLICIT NONE
c     include 'params'
c     include 'comgra' 

      integer maxix,maxiy,nix,niy,overlay_grid
      real vs(maxix,maxiy),raxis(maxix),zaxis(maxiy)
      integer ix,iy 
c
c slmod begin - new
c     INCLUDE 'slout'

      INTEGER i,j
c slmod end
c
      INTEGER ICNTR,NGS
      integer nconts,cntropt
      real conts(maxpts)
      REAL XXMIN,XXMAX,YYMIN,YYMAX
      real minscale,maxscale
C
C  *********************************************************************
C  *                                                                   *
C  *  CONTOURXY: IF ICNTR = 1                                          *
C  *           A FALSE COLOUR PLOT OF THE DATA IS GENERATED.           *
c
c
c     overlay_grid = 0 = off - grid is not overlaid
c                  = 1 = on  - call to supimp('SELECT') is made  
C  *                                                                   *
C  *********************************************************************
C
      INTEGER IG,tmpngs,increment,npos,nneg
      REAL    VMIN,VMAX,VLO,VHI
      real    tmpconts(maxpts)
      CHARACTER*36 NAME
c
      real logtable(21)
      data logtable /0.0e+0, 0.5e-6, 1.0e-6,
     >               0.2e-5, 0.5e-5, 1.0e-5,
     >               0.2e-4, 0.5e-4, 1.0e-4,
     >               0.2e-3, 0.5e-3, 1.0e-3,
     >               0.2e-2, 0.5e-2, 1.0e-2,
     >               0.2e-1, 0.5e-1, 1.0e-1,
     >               0.2e+0, 0.5e+0, 1.0e+0/
C
      if (minscale.ne.maxscale) then
        vmin=minscale
        vmax=maxscale
      else
        call find_minmax(vs,maxix,maxiy,nix,niy,vmin,vmax,
     >                   xxmin,xxmax,yymin,yymax,raxis,zaxis)
      endif
c
c     Original default option - quadratic
c
      if (cntropt.eq.0) then
         tmpngs = ngs
         do ig = 1,tmpngs
            tmpconts(ig) = REAL(IG)/REAL(NGS)/REAL(NGS+1-IG)
         end do
c
c     10% levels
c
      elseif (cntropt.eq.1) then
         tmpngs = 10
         do ig = 1,tmpngs
            tmpconts(ig) = ig * 0.1
         end do
c
c     Logarithmic
c
      elseif (cntropt.eq.2) then
         tmpngs = ngs
         do ig = 1,tmpngs
            tmpconts(ig) = real(2**ig)/real(2**ngs)
         end do
c
c     User specified
c
      elseif (cntropt.eq.3) then
c
         tmpngs = min(ngs, nconts)
         increment = nconts-tmpngs
         do ig = 1,tmpngs
            tmpconts(ig) = conts(increment+ig)
         end do  
c
c         tmpngs = nconts
c         do ig = 1,tmpngs
c            tmpconts(ig) = conts(ig)
c         end do
c
c     User specified - Absolute values
c
      elseif (cntropt.eq.4) then
c
         tmpngs = nconts
         do ig = 1,tmpngs
            tmpconts(ig) = conts(ig)
         end do
c
c     B2/Eirene Velocity plots 
c
      elseif (cntropt.eq.5) then 
c
         vmax=max(abs(vmin),abs(vmax))
         vmin=-vmax
c
         nneg=20
         npos=20
c
         do ig=1,nneg
           tmpconts(ig)=vmin*logtable(22-ig)
         enddo
         do ig=1,npos+1
           tmpconts(nneg+ig)=vmax*logtable(ig)
         enddo
c
         tmpngs=nneg+npos
c
      endif
c
c     Debug
c
      write (6,*) 'Contourxy:',vmin,vmax,tmpngs,ngs,nconts
c
c     Set maximum value
c
      if (cntropt.eq.4.or.cntropt.eq.5) then
         tmpconts(tmpngs+1) = vmax
      else
         tmpconts(tmpngs+1) = 1.0
      endif
c
c     jdemod - reorder contour plotting from lowest to highest 
c              to match color scheme (blue low and red high)
c              and to match what is done in the contour subroutine.
c
      DO IG = 1,tmpNGS
c      DO IG = tmpNGS,1,-1
c
        if (cntropt.eq.4.or.cntropt.eq.5) then
c
           VLO = tmpconts(ig)
           VHI = tmpconts(ig+1)
c
        else

           VLO = VMIN + tmpconts(ig)  *(VMAX-VMIN)
           VHI = VMIN + tmpconts(ig+1)*(VMAX-VMIN)
c
        endif
c
c       Plots with more than 20 contour levels
c
        if (tmpngs.gt.20) then 
c
           WRITE (NAME,'(4X,1P,E8.1)') VLO
           IF (IG.EQ.tmpNGS) THEN
              NAME(13:27) = '(MX= 0.0e-00)'
              write(name(17:24),'(1p,e8.1)') vhi
           endif
c
        else 
c
           WRITE (NAME,'(4X,1P,E8.1,1x,''TO'',1x,E8.1)') VLO,VHI
           IF (IG.EQ.TMPNGS) NAME(25:36) = ' (MAX)'
c
        endif
c
c slmod begin
        CALL GRCOLRXY (VS,maxix,maxiy,nix,niy,raxis,zaxis,
     >                   VLO,VHI,vmin,vmax,NAME)
c
c        CALL GRCOLRXY (VS,maxix,maxiy,nix,niy,raxis,zaxis,
c     >                   VLO,VHI,NAME)
c slmod end
c
      ENDDO

c
c slmod begin - new
c...  Print comments:
      IF (slopt.EQ.4) THEN
c...    For plot 982:
      ELSE
        CALL PSPACE (0.0, 1.35, 0.0,1.0)
        CALL CSPACE (0.0, 1.35, 0.0,1.0)
        CALL MAP    (0.0, 1.35, 0.0,1.0)
        CALL CTRMAG (12)
        DO i = 1, 10
          IF (char(i).NE.' ') CALL PLOTST(1.00,0.590+(i-1)*0.02,char(i))
        ENDDO
        DO i = 20, 30
          j = i - 19
          IF (char(i).NE.' ') CALL PLOTST(1.00,0.550-(j-1)*0.02,char(i))
        ENDDO

c...    For toroidal camera plot 972:
        IF (slopt.EQ.3) THEN
          iplots = 0
          RETURN
        ENDIF
      ENDIF
c slmod end

c
c     moved to end; in case of false color plots, the separatrix
c     was hidden by the coloring procedures - Krieger IPP/97
c
      iplots = 0
c
      if (overlay_grid.eq.1)  CALL SUPIMP('SELECT')
c
c slmod begin
c...  For plot 982:
      IF (slopt.EQ.4) RETURN
c slmod end
      CALL FRAME
C
      RETURN
      END
c
c
c
      subroutine find_minmax(vs,maxix,maxiy,nix,niy,
     >                       minval,maxval,
     >                       xmin,xmax,ymin,ymax,xaxis,yaxis) 
      use mod_params
      implicit none
c     include 'params'
      integer maxix,maxiy,nix,niy
      real vs(maxix,maxiy)
      real minval,maxval
      real xmin,xmax,ymin,ymax
      real xaxis(maxix),yaxis(maxiy)
c
c     FIND_MINMAX: Finds the minimum and maximum values in a 2D array
c
      integer ix,iy
c
      minval = HI
      maxval =-HI
c
      do ix = 1,nix
         do iy = 1,niy
            if ((xaxis(ix).ge.xmin.and.xaxis(ix).le.xmax).and.
     >          (yaxis(iy).ge.ymin.and.yaxis(iy).le.ymax)) then     
               if (vs(ix,iy).lt.minval) minval = vs(ix,iy)
               if (vs(ix,iy).gt.maxval) maxval = vs(ix,iy)
            endif
         end do          
      end do 
c
      return 
      end
C
C
C
      SUBROUTINE RVALKR (CVALS,VS,II,NIIS,MAXIIS,FT,FP,MFACT,
     >                   XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
      use mod_params
      use mod_cgeom
      IMPLICIT NONE
C     INCLUDE "PARAMS"
c     include 'params'
C     INCLUDE "CGEOM"
c     include 'cgeom'
      INTEGER II,NIIS,MAXIIS
      REAL CVALS(MAXNKS,MAXNRS),VS(MAXNKS,MAXNRS,MAXIIS),FT,FP,MFACT
      REAL XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX
C
C  *********************************************************************
C  *                                                                   *
C  *  RVALKR:  MULTIPURPOSE ROUTINE TO EXTRACT CVALS IN RING/KNOT      *
C  *           SPACE FOR PLOTTING.                                     *
C  *                                                                   *
C  *            CHRIS FARRELL  (HUNTERSKIL)  FEBRUARY 1989             *
C  *            LORNE HORTON   (JET)         JULY     1993             *
C  *                                                                   *
C  *********************************************************************
C
      INTEGER IR,IK,JJ
C
      VMIN = HI
      VMAX =-HI
      CALL RZERO (CVALS, MAXNKS*MAXNRS)
C
      DO 130 IR = 1, NRS
        DO 120 IK = 1, NKS(IR)
          IF (II.EQ.NIIS+1) THEN
C  SUM OVER IONISATION STATES
            DO 110 JJ = 2, NIIS
              CVALS(IK,IR) = CVALS(IK,IR) + VS(IK,IR,JJ)
  110       CONTINUE
          ELSEIF (II.EQ.0) THEN
C  SECONDARY NEUTRALS
            CVALS(IK,IR) = FT * VS(IK,IR,2) - FP * VS(IK,IR,1)
          ELSE
            CVALS(IK,IR) = VS(IK,IR,II)
          ENDIF
C  OVERALL SCALE FACTOR
          CVALS(IK,IR) = MFACT * CVALS(IK,IR)
C  FIND MINIMUM AND MAXIMUM ONLY FOR CELLS IN PLOT WINDOW
          IF ((RS(IK,IR).GE.XXMIN .AND. RS(IK,IR).LE.XXMAX) .AND.
     >        (ZS(IK,IR).GE.YYMIN .AND. ZS(IK,IR).LE.YYMAX)) THEN
            VMIN = MIN (VMIN, CVALS(IK,IR))
            VMAX = MAX (VMAX, CVALS(IK,IR))
          ENDIF
  120   CONTINUE
  130 CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE RVALXY (CVALS,VS,II,NIIS,MAXIIS,FT,FP,MFACT,
     >                   IXMIN,IXMAX,IYMIN,IYMAX,VMIN,VMAX)
      use mod_params
      use mod_cgeom
      implicit none
c      IMPLICIT LOGICAL (A-Z)
C     INCLUDE "PARAMS"
c     include 'params'
C     INCLUDE "CGEOM"
c     include 'cgeom'
      INTEGER II,NIIS,MAXIIS,IXMIN,IXMAX,IYMIN,IYMAX,IX,IY,JJ,IK,IR
      REAL VS(MAXNKS,MAXNRS,MAXIIS)
      REAL CVALS(MAXGXS,MAXGYS),FT,FP,MFACT,VMIN,VMAX
C
C  *********************************************************************
C  *                                                                   *
C  *  RVALXY:  MULTIPURPOSE ROUTINE TO EXTRACT CVALS FOR PLOTTING.     *
C  *                                                                   *
C  *            CHRIS FARRELL  (HUNTERSKIL)  FEBRUARY 1989             *
C  *            LORNE HORTON   (JET)         JULY     1993             *
C  *                                                                   *
C  *********************************************************************
C
      logical griderr
      real r,z
c
      ik = 0
      ir = 0
c
      VMIN = HI
      VMAX =-HI
      CALL RZERO (CVALS, MAXGXS*MAXGYS)
C
c
      DO 130 IX = IXMIN, IXMAX
        R = (RMAX-RMIN) * REAL(IX)/REAL(NXS) + RMIN - 0.5 * DR
        DO 120 IY = IYMIN, IYMAX
          Z = (ZMAX-ZMIN) * REAL(IY)/REAL(NYS) + ZMIN - 0.5 * DZ
c
          call gridpos(ik,ir,r,z,.false.,griderr)
c
          if (griderr) goto 120
c
c          IK = IKXYS(IX,IY)
c          IR = IRXYS(IX,IY)
c          IF (IFXYS(IX,IY).EQ.0) GOTO 120
c
          IF (II.EQ.NIIS+1) THEN
            DO 110 JJ = 2, NIIS
              CVALS(IX,IY) = CVALS(IX,IY) + VS(IK,IR,JJ)
  110       CONTINUE
          ELSEIF (II.EQ.0) THEN
            CVALS(IX,IY) = FT * VS(IK,IR,2) - FP * VS(IK,IR,1)
          ELSE
            CVALS(IX,IY) = VS(IK,IR,II)
          ENDIF
C  OVERALL SCALE FACTOR
          CVALS(IX,IY) = MFACT * CVALS(IX,IY)
C
          VMIN = MIN (VMIN, CVALS(IX,IY))
          VMAX = MAX (VMAX, CVALS(IX,IY))
c
  120   CONTINUE
  130 CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE HCONTOUR (ICNTR,NGS,VS,II,NIIS,MAXIIS,FT,FP,MFACT,
     >                     XOUTS,IXMIN,IXMAX,YOUTS,IYMIN,IYMAX,
     >                     XXMIN,XXMAX,YYMIN,YYMAX,
     >                     conts,nconts,cntropt)
      use mod_params
      use mod_cgeom
      IMPLICIT NONE
C     INCLUDE "PARAMS"
c     include 'params'
C     INCLUDE "CGEOM"
c     include 'cgeom'
c
      INTEGER ICNTR,NGS,II,NIIS,MAXIIS,IXMIN,IXMAX,IYMIN,IYMAX
      integer nconts,cntropt
      real conts(maxpts)
      REAL VS(MAXNKS,MAXNRS,MAXIIS),FT,FP,MFACT
      REAL XOUTS(MAXGXS),YOUTS(MAXGYS)
      REAL XXMIN,XXMAX,YYMIN,YYMAX
C
C  *********************************************************************
C  *                                                                   *
C  * HCONTOUR: IF ICNTR = 0, GENERATES A CONTOUR PLOT BY MAPPING       *
C  *           THE INPUT ARRAY VS ONTO THE XY GRID.  IF ICNTR = 1      *
C  *           A FALSE COLOUR PLOT OF THE SAME DATA IS GENERATED.      *
C  *                                                                   *
C  *            LORNE HORTON   (JET)         JULY     1993             *
C  *                                                                   *
C  *********************************************************************
C
      INTEGER IG,tmpngs
      REAL    VMIN,VMAX,VLO,VHI
      REAL    CVALXY(MAXGXS,MAXGYS), CVALKR(MAXNKS,MAXNRS)
      real    tmpconts(maxpts)
      CHARACTER*36 NAME
C
      IF (ICNTR.EQ.0) THEN
        CALL HVALXY(CVALXY,VS,II,NIIS,MAXIIS,FT,FP,MFACT,
     >              IXMIN,IXMAX,IYMIN,IYMAX,VMIN,VMAX)
      ELSE
        CALL HVALKR(CVALKR,VS,II,NIIS,MAXIIS,FT,FP,MFACT,
     >              XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
      ENDIF
      CALL SUPIMP('SELECT')
c
c     Set up the contour information.
c
c
c     Original default option - quadratic
c
      if (cntropt.eq.0) then
         tmpngs = ngs
         do ig = 1,ngs
            tmpconts(ig) = REAL(IG)/REAL(NGS)/REAL(NGS+1-IG)
         end do
c
c     10% levels
c
      elseif (cntropt.eq.1) then
         tmpngs = 10
         do ig = 1,ngs
            tmpconts(ig) = ig * 0.1
         end do
c
c     Logarithmic
c
      elseif (cntropt.eq.2) then
         tmpngs = ngs
         do ig = 1,ngs
            tmpconts(ig) = real(2**ig)/real(2**ngs)
         end do
c
c     User specified
c
      elseif (cntropt.eq.3.or.cntropt.eq.4) then
         tmpngs = nconts
         do ig = 1,ngs
            tmpconts(ig) = conts(ig)
         end do
      endif
c
c     Set maximum value
c
      if (cntropt.eq.4) then 
         tmpconts(tmpngs+1) = vmax
      else
         tmpconts(tmpngs+1) = 1.0
      endif
c
c
c      DO IG = 1,tmpngs
      DO IG = tmpngs,1,-1
c
c       let ngs=10 be a special case for 10% levels
c
        if (cntropt.eq.4) then 
           VLO = tmpconts(ig) 
           VHI = tmpconts(ig+1)
        else 
           VLO = VMIN + tmpconts(ig)  *(VMAX-VMIN)
           VHI = VMIN + tmpconts(ig+1)*(VMAX-VMIN)
        endif 
c
        WRITE (NAME,'(4X,1P,E8.1)') VLO
        IF (IG.EQ.NGS) NAME(13:18) = ' (MAX)'
c
        IF (ICNTR.EQ.0) THEN
          CALL GRCONT (CVALXY,IXMIN,IXMAX,MAXGXS,IYMIN,IYMAX,MAXGYS,
     >                 VLO,XOUTS,YOUTS,NAME)
        ELSE
          CALL GRCOLR (CVALKR,VLO,VHI,NAME)
        ENDIF
      ENDDO
      CALL FRAME
C
      RETURN
      END
C
C
C
      SUBROUTINE HVALKR (CVALS,VS,II,NIIS,MAXIIS,FT,FP,MFACT,
     >                   XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
      use mod_params
      use mod_cgeom
      IMPLICIT NONE
C     INCLUDE "PARAMS"
c     include 'params'
C     INCLUDE "CGEOM"
c     include 'cgeom'
      INTEGER II,NIIS,MAXIIS
      REAL CVALS(MAXNKS,MAXNRS),VS(MAXNKS,MAXNRS,MAXIIS),FT,FP,MFACT
      REAL XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX
C
C  *********************************************************************
C  *                                                                   *
C  *  HVALKR:  MULTIPURPOSE ROUTINE TO EXTRACT CVALS IN RING/KNOT      *
C  *           SPACE FOR PLOTTING.  VERSION FOR HYDROGEN DATA.         *
C  *                                                                   *
C  *            CHRIS FARRELL  (HUNTERSKIL)  FEBRUARY 1989             *
C  *            LORNE HORTON   (JET)         JULY     1993             *
C  *                                                                   *
C  *********************************************************************
C
      INTEGER IR,IK,JJ
C
      VMIN = HI
      VMAX =-HI
      CALL RZERO (CVALS, MAXNKS*MAXNRS)
C
      DO 130 IR = 1, NRS
        DO 120 IK = 1, NKS(IR)
          IF (II.EQ.NIIS+1) THEN
C  SUM OVER IONISATION STATES
            DO 110 JJ = 1, NIIS
              CVALS(IK,IR) = CVALS(IK,IR) + VS(IK,IR,JJ)
  110       CONTINUE
          ELSE
            CVALS(IK,IR) = VS(IK,IR,II)
          ENDIF
C  OVERALL SCALE FACTOR
          CVALS(IK,IR) = MFACT * CVALS(IK,IR)
C  FIND MINIMUM AND MAXIMUM ONLY FOR CELLS IN PLOT WINDOW
          IF ((RS(IK,IR).GE.XXMIN .AND. RS(IK,IR).LE.XXMAX) .AND.
     >        (ZS(IK,IR).GE.YYMIN .AND. ZS(IK,IR).LE.YYMAX)) THEN
            VMIN = MIN (VMIN, CVALS(IK,IR))
            VMAX = MAX (VMAX, CVALS(IK,IR))
          ENDIF
  120   CONTINUE
  130 CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE HVALXY (CVALS,VS,II,NIIS,MAXIIS,FT,FP,MFACT,
     >                   IXMIN,IXMAX,IYMIN,IYMAX,VMIN,VMAX)
      use mod_params
      use mod_cgeom
      implicit none
c      IMPLICIT LOGICAL (A-Z)
C     INCLUDE "PARAMS"
c     include 'params'
C     INCLUDE "CGEOM"
c     include 'cgeom'
      INTEGER II,NIIS,MAXIIS,IXMIN,IXMAX,IYMIN,IYMAX,IX,IY,JJ,IK,IR
      REAL VS(MAXNKS,MAXNRS,MAXIIS)
      REAL CVALS(MAXGXS,MAXGYS),FT,FP,MFACT,VMIN,VMAX
C
C  *********************************************************************
C  *                                                                   *
C  *  HVALXY:  MULTIPURPOSE ROUTINE TO EXTRACT CVALS FOR PLOTTING.     *
C  *           VERSION FOR HYDROGEN DATA.                              *
C  *                                                                   *
C  *            CHRIS FARRELL  (HUNTERSKIL)  FEBRUARY 1989             *
C  *            LORNE HORTON   (JET)         JULY     1993             *
C  *                                                                   *
C  *********************************************************************
C
      logical griderr
      real r,z
c
      ik = 0
      ir = 0
c
      VMIN = HI
      VMAX =-HI
      CALL RZERO (CVALS, MAXGXS*MAXGYS)
C
      DO 130 IX = IXMIN, IXMAX
        R = (RMAX-RMIN) * REAL(IX)/REAL(NXS) + RMIN - 0.5 * DR
        DO 120 IY = IYMIN, IYMAX
          Z = (ZMAX-ZMIN) * REAL(IY)/REAL(NYS) + ZMIN - 0.5 * DZ
          call gridpos(ik,ir,r,z,.false.,griderr)
          if (griderr) goto 120
c
c          IK = IKXYS(IX,IY)
c          IR = IRXYS(IX,IY)
c          IF (IFXYS(IX,IY).EQ.0) GOTO 120
c
          IF (II.EQ.NIIS+1) THEN
            DO 110 JJ = 1, NIIS
              CVALS(IX,IY) = CVALS(IX,IY) + VS(IK,IR,JJ)
  110       CONTINUE
          ELSE
            CVALS(IX,IY) = VS(IK,IR,II)
          ENDIF
C  OVERALL SCALE FACTOR
          CVALS(IX,IY) = MFACT * CVALS(IX,IY)
C
          VMIN = MIN (VMIN, CVALS(IX,IY))
          VMAX = MAX (VMAX, CVALS(IX,IY))
  120   CONTINUE
  130 CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE SUPIMP2 (OPTION)
      use mod_params
      use mod_cgeom
      use mod_pindata
      use mod_comtor
      use mod_colours
      use mod_slcom
      IMPLICIT NONE
C     INCLUDE "PARAMS"
c     include 'params'
      CHARACTER*(*) OPTION
c
c     NOTE: This superposition routine is dependent on data that is
c           available only for JET grids at this time. Furthermore
c           since the alternative SUPIMP routine that originally
c           was in use can fulfill this function for both JET grids
c           as well ITER and ASDEX grids without losing its
c           generality ... it was decided to use the original in place
c           of this routine for all grid superposition drawing.
c           This routine is still available if it proves useful
c           in the future to shift to this method of plotting the
c           grid lines. The benefit of this routine is that it uses
c           the actual cell boundaries for the grid overlays.
c
c           D. Elder  /  Oct 21 1993
c
C
C  *********************************************************************
C  *                                                                   *
C  *  SUPIMP2: ROUTINE TO SUPERIMPOSE THE EQUILIBRIUM GRID ON A        *
C  *           CONTOUR PLOT.                                           *
C  *                                                                   *
C  *            CHRIS FARRELL  (HUNTERSKIL)  MARCH 1989                *
C  *            LORNE HORTON   (JET)         JULY  1993                *
C  *                                                                   *
C  *********************************************************************
C
c     include 'cgeom'
c     include 'pindata'
c slmod begin
c     include 'comtor'
c     include 'colours'
c     include 'slcom'

      INTEGER   rind,rcol,rmode
      CHARACTER dummy*24,cdum1*1024,cdum2*32
      LOGICAL   reset_colour,tindex,barebones
      REAL      xcen,ycen
c slmod end
c
C     IPP/01 - Krieger: fixed initialization in declaration statement
C     by adding appropiate data statement (picky SUN compiler)
      CHARACTER*36 NAME
      INTEGER I,J,IK,K,IR,ID,fp
      REAL KVALS(MAXNKS+1,2),RVALS(MAXNRS+1,2),PVALS(5,2),
     >     wvals(maxseg+1,2)
c
      integer cnt,bafflecnt
c
      DATA NAME /' '/

c
c      REAL KVALS(MAXNKS+1,2),RVALS(MAXNRS+1,2),PVALS(5,2)
C
      IF (OPTION.EQ.'FULL') THEN
        write (6,*) 'FULL:'

        CALL CTRMAG(8) 
        call setup_col(ncols,5)
c slmod begin
        tindex = .FALSE.
 24     READ(5,'(A512)') cdum1
        IF   (cdum1(8:11).EQ.'Tind'.OR.cdum1(8:11).EQ.'tind'.OR.
     .        cdum1(8:11).EQ.'TIND') THEN
          tindex = .TRUE.
        ELSE
          BACKSPACE 5
        ENDIF          

c slmod begin
        DO ir = 2, nrs
c
c        DO ir = 1, nrs
c slmod end

           IF (ir.EQ.irwall) CYCLE
           IF (idring(ir).LT.0) CYCLE
c
c        DO ir = 2, nrs
c           IF (ir.NE.49) CYCLE
c slmod end
           DO ik = 1, nks(ir)

              DEFCOL = NCOLS+1

              i = korpg(ik,ir)
c
c             Check for valid polygons and numbers of sides
c             before plotting
c
              if (i.gt.0) then 
c
                 if (nvertp(i).gt.0) then  

c slmod begin
                    IF (nbr.GT.0) THEN
                      IF (ik.EQ.1) THEN
                        PVALS(1,1) = RVERTP(1,I)
                        PVALS(1,2) = ZVERTP(1,I)
                        PVALS(2,1) = RVERTP(2,I)
                        PVALS(2,2) = ZVERTP(2,I)
                        CALL GRTRAC (PVALS(1,1),PVALS(1,2),
     .                               2,NAME,'LINE',-1)        
                      ENDIF
                      PVALS(1,1) = RVERTP(3,I)
                      PVALS(1,2) = ZVERTP(3,I)
                      PVALS(2,1) = RVERTP(4,I)
                      PVALS(2,2) = ZVERTP(4,I)
                      CALL GRTRAC (PVALS(1,1),PVALS(1,2),
     .                             2,NAME,'LINE',-1)        
                      IF (ir.EQ.2.OR.ir.EQ.irtrap+1) THEN
                        PVALS(1,1) = RVERTP(1,I)
                        PVALS(1,2) = ZVERTP(1,I)
                        PVALS(2,1) = RVERTP(4,I)
                        PVALS(2,2) = ZVERTP(4,I)
                        CALL GRTRAC (PVALS(1,1),PVALS(1,2),
     .                               2,NAME,'LINE',-1)        
                      ENDIF
                      PVALS(1,1) = RVERTP(2,I)
                      PVALS(1,2) = ZVERTP(2,I)
                      PVALS(2,1) = RVERTP(3,I)
                      PVALS(2,2) = ZVERTP(3,I)
                      CALL GRTRAC (PVALS(1,1),PVALS(1,2),
     .                             2,NAME,'LINE',-1)        
                    ELSE
                      DO J = 1, NVERTP(I)
                         PVALS(J,1) = RVERTP(J,I)
                         PVALS(J,2) = ZVERTP(J,I)
                      ENDDO
                      PVALS(NVERTP(I)+1,1) = PVALS(1,1)
                      PVALS(NVERTP(I)+1,2) = PVALS(1,2)
                      CALL GRTRAC (PVALS(1,1),PVALS(1,2),
     .                          NVERTP(I)+1,NAME,'LINE',-1)
                    ENDIF
c
c                     DO J = 1, NVERTP(I)
c                       PVALS(J,1) = RVERTP(J,I)
c                       PVALS(J,2) = ZVERTP(J,I)
c                    ENDDO
c                    PVALS(NVERTP(I)+1,1) = PVALS(1,1)
c                    PVALS(NVERTP(I)+1,2) = PVALS(1,2)
c                    CALL GRTRAC (PVALS(1,1),PVALS(1,2),
c     .                        NVERTP(I)+1,NAME,'LINE',-1)
c slmod end
                 endif 

              endif 


c             IF (.TRUE..AND.ik.EQ.111.AND.ir.EQ.18) THEN
c                CALL CTRMAG(10)
c                WRITE(dummy,'(A1)') '+'
c                CALL PCSCEN(rs(ik,ir),zs(ik,ir),
c     .                      dummy(1:LEN_TRIM(dummy)))
c             ENDIF



           ENDDO

c...       Show ring index at target:
           IF (tindex) THEN
              call lincol(ncols+3)
              WRITE(dummy,'(I3)') ir
              CALL PCSCEN(rs(1      ,ir),zs(1      ,ir),
     .                    dummy(1:LEN_TRIM(dummy)))
              CALL PCSCEN(rs(nks(ir),ir),zs(nks(ir),ir),
     .                    dummy(1:LEN_TRIM(dummy)))
           ENDIF

        ENDDO

c...    Highlight a ring:
 25     READ(5,'(A512)') cdum1
        IF   (cdum1(8:11).EQ.'Ring'.OR.cdum1(8:11).EQ.'ring'.OR.
     .        cdum1(8:11).EQ.'RING') THEN
          READ(cdum1,*) cdum2,rind,rcol,rmode

          IF     (rind.EQ.-1 ) THEN
            ir = irsep
          ELSEIF (rind.EQ.-2 ) THEN
            ir = irsep2
          ELSEIF (rind.LT.-10) THEN
            ir = ir - (rind + 10)
          ELSE
            ir = rind
          ENDIF
          DEFCOL = NCOLS + rcol
          DO ik = 1, nks(ir)
             IF (rmode.EQ.1.AND.ik.NE.ikto.AND.ik.NE.ikti) CYCLE
             IF (idring(ir).EQ.BOUNDARY) THEN
               IF (irins(ik,ir).EQ.ir) THEN
                 i = korpg(ikouts(ik,ir),irouts(ik,ir))
               ELSE
                 i = korpg(ikins(ik,ir),irins(ik,ir))
               ENDIF
             ELSE
               i = korpg(ik,ir)
             ENDIF
             if (i.gt.0) then 
                if (nvertp(i).gt.0) then  
                  IF (.FALSE.) THEN
                    PVALS(1,1) = RVERTP(1,I)
                    PVALS(1,2) = ZVERTP(1,I)
                    PVALS(2,1) = RVERTP(4,I)
                    PVALS(2,2) = ZVERTP(4,I)
                    CALL GRTRAC (PVALS(1,1),PVALS(1,2),
     .                           2,NAME,'LINE',-1)        
                  ELSEIF (idring(ir).EQ.BOUNDARY) THEN
                    IF (irins(ik,ir).EQ.ir) THEN
                      PVALS(1,1) = RVERTP(1,I)
                      PVALS(1,2) = ZVERTP(1,I)
                      PVALS(2,1) = RVERTP(4,I)
                      PVALS(2,2) = ZVERTP(4,I)
                    ELSE
                      PVALS(1,1) = RVERTP(2,I)
                      PVALS(1,2) = ZVERTP(2,I)
                      PVALS(2,1) = RVERTP(3,I)
                      PVALS(2,2) = ZVERTP(3,I)
                    ENDIF
                    CALL GRTRAC (PVALS(1,1),PVALS(1,2),
     .                           2,NAME,'LINE',-1)
                  ELSE
                    DO J = 1, NVERTP(I)
                      PVALS(J,1) = RVERTP(J,I)
                      PVALS(J,2) = ZVERTP(J,I)
                    ENDDO
                    PVALS(NVERTP(I)+1,1) = PVALS(1,1)
                    PVALS(NVERTP(I)+1,2) = PVALS(1,2)
                    CALL GRTRAC (PVALS(1,1),PVALS(1,2),
     .                           NVERTP(I)+1,NAME,'LINE',-1)
                  ENDIF

                endif 
             endif 
          ENDDO

          GOTO 25
        ELSE
          DEFCOL = 0
          BACKSPACE 5
        ENDIF
c
c        DO I = 1, NPOLYP
c          DO J = 1, NVERTP(I)
c            PVALS(J,1) = RVERTP(J,I)
c            PVALS(J,2) = ZVERTP(J,I)
c          ENDDO
c          PVALS(NVERTP(I)+1,1) = PVALS(1,1)
c          PVALS(NVERTP(I)+1,2) = PVALS(1,2)
c          CALL GRTRAC (PVALS(1,1),PVALS(1,2),NVERTP(I)+1,NAME,'LINE',-1)
c          write (6,*) npolyp,i,nvertp(i),pvals(1,1),pvals(1,2)
c
c        ENDDO
c slmod end
c
c also plot the nimbus wall and pump (if any)
c
c       IF available
c

        if (nvesm.ne.0) then
c
c         Need to modify this to plot the baffles independently if they
c         are present.
c
c         Plot each segment separately - as is done with the pump structure
c
          DEFCOL = NCOLS + 4
c          DEFCOL = 1  ! BLACK

          do i = 1,nvesm
c             if (nbr.GT.0.AND.jvesm(i).ne.0) CYCLE
             if (jvesm(i).ne.9) then
                wvals(1,1) = rvesm(i,1)
                wvals(1,2) = zvesm(i,1)
                wvals(2,1) = rvesm(i,2)
                wvals(2,2) = zvesm(i,2)
                call grtrac (wvals(1,1),wvals(1,2),2,
     >                       name,'LINE',-1)
             endif
          enddo

          DEFCOL = 1
c
c
c       OLD method
c
c        if (nvesm.ne.0) then
c          do i = 1,nvesm
c            wvals(i,1) = rvesm(i,1)
c            wvals(i,2) = zvesm(i,1)
c          enddo
c          wvals(nvesm+1,1) = rvesm(nvesm,2)
c          wvals(nvesm+1,2) = zvesm(nvesm,2)
c          call grtrac (wvals(1,1),wvals(1,2),nvesm+1,name,'LINE',-1)
c
c
c
c the pump is not a continuous nose-to-tail structure
c
          do i = 1,nvesp
            wvals(1,1) = rvesm(nvesm+i,1)
            wvals(1,2) = zvesm(nvesm+i,1)
            wvals(2,1) = rvesm(nvesm+i,2)
            wvals(2,2) = zvesm(nvesm+i,2)
            call grtrac (wvals(1,1),wvals(1,2),2,name,'LINE',-1)
          enddo
c
c       Plot the wall read from the grid file if NIMBUS wall is
c       not available
c
        elseif (nves.ne.0) then

          write (6,*) 'Nves:',nves

          do i = 1,nves
            wvals(i,1) = rves(i)
            wvals(i,2) = zves(i)
          enddo
          wvals(nves+1,1) = rves(1)
          wvals(nves+1,2) = zves(1)
          call grtrac (wvals(1,1),wvals(1,2),nves+1,name,'LINE',-1)
        endif
c

c slmod begin
      ELSEIF (OPTION.EQ.'TRIANGLES') THEN
        write (6,*) 'TRIANGLES:'

        fp = 99
        OPEN(UNIT=fp,FILE='triangles.dat',ACCESS='SEQUENTIAL',
     .       STATUS='OLD',ERR=99)      

        reset_colour = .TRUE.
        READ(fp,*) i
c        WRITE(0,*) 'i:',i
        DO j = 1, i
          READ(fp,*) id,(pvals(k,1),pvals(k,2),k=1,3)
c          WRITE(0,*) id,(pvals(k,1),pvals(k,2),k=1,3)

c          IF (id.NE.5342) CYCLE

c          CYCLE

c          IF     (id.EQ.4019) THEN
c            DEFCOL = NCOLS+1
c            reset_colour = .TRUE.
c          ELSEIF (reset_colour) THEN
c            DEFCOL = 0
c            reset_colour = .FALSE.
c          ENDIF

c          READ(fp,*) id,(pvals(k,1),pvals(k,2),k=1,3),xcen,ycen
          PVALS(4,1) = PVALS(1,1)
          PVALS(4,2) = PVALS(1,2)
          CALL GRTRAC (PVALS(1,1),PVALS(1,2),4,NAME,'LINE',-1)
          PVALS(1,1) = xcen
          PVALS(1,2) = ycen
          PVALS(2,1) = xcen * 1.001
          PVALS(2,2) = ycen 
          CALL GRTRAC (PVALS(1,1),PVALS(1,2),2,NAME,'LINE',-1)
        ENDDO
        CLOSE(fp)              

        CALL RGB
        CALL COLSET(1.0,0.0,0.0,NCOLS+1)


 21     READ(5,'(A512)') cdum1
        IF   (cdum1(8:11).EQ.'Show'.OR.cdum1(8:11).EQ.'show'.OR.
     .        cdum1(8:11).EQ.'SHOW') THEN
          READ(cdum1,*) cdum2,rind
          WRITE(0,*) 'SHOWING TRIANGLE: ',rind
          DEFCOL = NCOLS+1
          OPEN(UNIT=fp,FILE='triangles.dat',ACCESS='SEQUENTIAL',
     .         STATUS='OLD',ERR=99)      
          READ(fp,*) i
          DO j = 1, i
            READ(fp,*) id,(pvals(k,1),pvals(k,2),k=1,3)
            IF (id.NE.rind) CYCLE
            PVALS(4,1) = PVALS(1,1)
            PVALS(4,2) = PVALS(1,2)
            CALL GRTRAC (PVALS(1,1),PVALS(1,2),4,NAME,'LINE',-1)
            PVALS(1,1) = xcen
            PVALS(1,2) = ycen
            PVALS(2,1) = xcen * 1.001
            PVALS(2,2) = ycen 
            CALL GRTRAC (PVALS(1,1),PVALS(1,2),2,NAME,'LINE',-1)
          ENDDO
          CLOSE(fp)              
          GOTO 21
        ELSE
          BACKSPACE 5
        ENDIF          
          

        DEFCOL = 0


c
c also plot the nimbus wall and pump (if any)
c
c       IF available
c

c       DEFCOL = NCOLS+3
        IF (.FALSE.) THEN

        if (nvesm.ne.0) then
c
c         Need to modify this to plot the baffles independently if they
c         are present.
c
c         Plot each segment separately - as is done with the pump structure
c
          do i = 1,nvesm
c             if (nbr.GT.0.AND.jvesm(i).NE.0) CYCLE
             if (jvesm(i).ne.9) then
                wvals(1,1) = rvesm(i,1)
                wvals(1,2) = zvesm(i,1)
                wvals(2,1) = rvesm(i,2)
                wvals(2,2) = zvesm(i,2)
                call grtrac (wvals(1,1),wvals(1,2),2,
     >                       name,'LINE',-1)
             endif
          enddo
c
c the pump is not a continuous nose-to-tail structure
c
          do i = 1,nvesp
            wvals(1,1) = rvesm(nvesm+i,1)
            wvals(1,2) = zvesm(nvesm+i,1)
            wvals(2,1) = rvesm(nvesm+i,2)
            wvals(2,2) = zvesm(nvesm+i,2)
            call grtrac (wvals(1,1),wvals(1,2),2,name,'LINE',-1)
          enddo
c
c       Plot the wall read from the grid file if NIMBUS wall is
c       not available
c
        elseif (nves.ne.0) then

          write (6,*) 'Nves:',nves

          do i = 1,nves
            wvals(i,1) = rves(i)
            wvals(i,2) = zves(i)
          enddo
          wvals(nves+1,1) = rves(1)
          wvals(nves+1,2) = zves(1)
          call grtrac (wvals(1,1),wvals(1,2),nves+1,name,'LINE',-1)
        endif

        ENDIF
c        DEFCOL = 0

c slmod end
      ELSE

c
c      jdemod - nbr does not seem to be properly set automatically
c               for broken/extended grids. Need to fix this since
c               it affects some plot options among other things.
c
c       write(0,*) 'SUPIMP2:Select',nbr,cgridopt
c       nbr =1 
c
        IR = IRWALL - 1
c slmod begin
c...    Changes made to accommodate broken grids (look
c       out, here they come): 
        barebones = .FALSE.
        IF     (barebones) THEN
        ELSEIF (NBR.GT.0.OR.CGRIDOPT.EQ.LINEAR_GRID.OR.
     .                      CGRIDOPT.EQ.RIBBON_GRID) THEN
          DEFCOL = NCOLS + 1
	  IR = IRWALL
          DO IK = 1, NKS(IR)
            K = KORPG(IKINS(IK,IR),IRINS(IK,IR))
            IF     (IROUTS(IKINS(IK,IR),IRINS(IK,IR)).EQ.IRWALL) THEN
              KVALS(1,1) = RVERTP(2,K)
              KVALS(1,2) = ZVERTP(2,K)
              KVALS(2,1) = RVERTP(3,K)
              KVALS(2,2) = ZVERTP(3,K)
            ELSEIF (IRINS(IKINS(IK,IR),IRINS(IK,IR)).EQ.IRWALL) THEN
              KVALS(1,1) = RVERTP(1,K)
              KVALS(1,2) = ZVERTP(1,K)
              KVALS(2,1) = RVERTP(4,K)
              KVALS(2,2) = ZVERTP(4,K)
            ELSE
              WRITE(0,*) 'ERROR: PROBLEM WITH CONNECTION MAP IN SUPIMP2'
              STOP
            ENDIF
            CALL GRTRAC(KVALS(1,1),KVALS(1,2),2,NAME,'LINE',-1)
          ENDDO
          DEFCOL = 0
          DEFCOL = NCOLS+2
        ELSE
          DO IK = 1, NKS(IR)
            K = KORPG(IK,IR)
            KVALS(IK,1) = RVERTP(2,K)
            KVALS(IK,2) = ZVERTP(2,K)
          ENDDO
          K = KORPG(NKS(IR),IR)
          KVALS(NKS(IR)+1,1) = RVERTP(3,K)
          KVALS(NKS(IR)+1,2) = ZVERTP(3,K)
          CALL GRTRAC(KVALS(1,1),KVALS(1,2),NKS(IR)+1,NAME,'LINE',-1)
        ENDIF

        IF (.NOT.barebones.AND.IRSEP.NE.IRSEP2.AND.IRSEP2.GT.0) THEN
          IR = IROUTS(1,IRSEP2)
          DO IK = 1, NKS(IR)
            K = KORPG(IK,IR)
            KVALS(IK,1) = RVERTP(1,K)
            KVALS(IK,2) = ZVERTP(1,K)
          ENDDO
          K = KORPG(NKS(IR),IR)
          KVALS(NKS(IR)+1,1) = RVERTP(4,K)
          KVALS(NKS(IR)+1,2) = ZVERTP(4,K)
          CALL GRTRAC(KVALS(1,1),KVALS(1,2),NKS(IR)+1,NAME,'LINE',-1)
          IR = IROUTS(NKS(IRSEP2),IRSEP2)
          DO IK = 1, NKS(IR)
            K = KORPG(IK,IR)
            KVALS(IK,1) = RVERTP(1,K)
            KVALS(IK,2) = ZVERTP(1,K)
          ENDDO
          K = KORPG(NKS(IR),IR)
          KVALS(NKS(IR)+1,1) = RVERTP(4,K)
          KVALS(NKS(IR)+1,2) = ZVERTP(4,K)
          CALL GRTRAC(KVALS(1,1),KVALS(1,2),NKS(IR)+1,NAME,'LINE',-1)
        ENDIF
c
c        DO IK = 1, NKS(IR)
c          K = KORPG(IK,IR)
c          KVALS(IK,1) = RVERTP(2,K)
c          KVALS(IK,2) = ZVERTP(2,K)
c        ENDDO
c        K = KORPG(NKS(IR),IR)
c        KVALS(NKS(IR)+1,1) = RVERTP(3,K)
c        KVALS(NKS(IR)+1,2) = ZVERTP(3,K)
c        CALL GRTRAC(KVALS(1,1),KVALS(1,2),NKS(IR)+1,NAME,'LINE',-1)
c slmod end
C
        IF (.NOT.barebones) THEN
          IR = 2
          DO IK = 1, NKS(IR)
            K = KORPG(IK,IR)
            KVALS(IK,1) = RVERTP(4,K)
            KVALS(IK,2) = ZVERTP(4,K)
          ENDDO
          CALL GRTRAC(KVALS(1,1),KVALS(1,2),NKS(IR),NAME,'LINE',-1)
        ENDIF
C
        IF (CGRIDOPT.NE.LINEAR_GRID.AND..NOT.NOPRIV.AND.
     .      CGRIDOPT.NE.RIBBON_GRID) THEN

          IF (.NOT.barebones) THEN 
            IR = IRTRAP + 1
            DO IK = 1, NKS(IR)
              K = KORPG(IK,IR)
              KVALS(IK,1) = RVERTP(1,K)
              KVALS(IK,2) = ZVERTP(1,K)
            ENDDO
            K = KORPG(NKS(IR),IR)
            KVALS(NKS(IR)+1,1) = RVERTP(4,K)
            KVALS(NKS(IR)+1,2) = ZVERTP(4,K)
            CALL GRTRAC(KVALS(1,1),KVALS(1,2),NKS(IR)+1,NAME,'LINE',-1)
          ENDIF
C
          IR = IRSEP
          DEFCOL = NCOLS + 3
          DO IK = 1, NKS(IR)
            K = KORPG(IK,IR)
            KVALS(IK,1) = RVERTP(1,K)
            KVALS(IK,2) = ZVERTP(1,K)
          ENDDO
          K = KORPG(NKS(IR),IR)
          KVALS(NKS(IR)+1,1) = RVERTP(4,K)
          KVALS(NKS(IR)+1,2) = ZVERTP(4,K)
          CALL GRTRAC(KVALS(1,1),KVALS(1,2),NKS(IR)+1,NAME,'LINE',-1)
          DEFCOL = NCOLS + 1
        ENDIF
C
c slmod begin
        IF     (barebones) THEN
C
C       jdemod - the following code should work for all grids both broken and not ... so might as well use it since it
c                will work for unbroken double null grids (which have issues with the base code below).
C
        ELSEIF (.true.) THEN
c
c        ELSEIF (NBR.GT.0) THEN
c...      Draw targets:
c        
c
          DEFCOL = NCOLS + 2
          DO IR = IRSEP, NRS
            IF (IDRING(IR).EQ.-1) CYCLE
            K = KORPG(1,IR)
            RVALS(1,1) = RVERTP(1,K)
            RVALS(1,2) = ZVERTP(1,K)
            RVALS(2,1) = RVERTP(2,K)
            RVALS(2,2) = ZVERTP(2,K)
            CALL GRTRAC(RVALS(1,1),RVALS(1,2),2,NAME,'LINE',-1)
            K = KORPG(NKS(IR),IR)
            RVALS(1,1) = RVERTP(3,K)
            RVALS(1,2) = ZVERTP(3,K)
            RVALS(2,1) = RVERTP(4,K)
            RVALS(2,2) = ZVERTP(4,K)
            CALL GRTRAC(RVALS(1,1),RVALS(1,2),2,NAME,'LINE',-1)
          ENDDO
c...      Restore black ink:
          DEFCOL = 1
        ELSE

            ID = 0
            DO IR = IRTRAP+1,NRS
              ID = ID + 1
              K = KORPG(1,IR)
              RVALS(ID,1) = RVERTP(1,K)
              RVALS(ID,2) = ZVERTP(1,K)
            ENDDO
            DO IR = IRSEP,IRWALL-1
              ID = ID + 1
              K = KORPG(1,IR)
              RVALS(ID,1) = RVERTP(1,K)
              RVALS(ID,2) = ZVERTP(1,K)
            ENDDO
            ID = ID + 1
            K = KORPG(1,IRWALL-1)
            RVALS(ID,1) = RVERTP(2,K)
            RVALS(ID,2) = ZVERTP(2,K)
            CALL GRTRAC(RVALS(1,1),RVALS(1,2),ID,NAME,'LINE',-1)
C
            ID = 0
            DO IR = IRSEP,IRWALL-1
              ID = ID + 1
              K = KORPG(NKS(IR),IR)
              RVALS(ID,1) = RVERTP(4,K)
              RVALS(ID,2) = ZVERTP(4,K)
            ENDDO
            ID = ID + 1
            K = KORPG(NKS(IRWALL-1),IRWALL-1)
            RVALS(ID,1) = RVERTP(3,K)
            RVALS(ID,2) = ZVERTP(3,K)
            CALL GRTRAC(RVALS(1,1),RVALS(1,2),ID,NAME,'LINE',-1)

        ENDIF

        DEFCOL = 0
c
c        ID = 0
c        DO IR = IRTRAP+1,NRS
c          ID = ID + 1
c          K = KORPG(1,IR)
c          RVALS(ID,1) = RVERTP(1,K)
c          RVALS(ID,2) = ZVERTP(1,K)
c        ENDDO
c        DO IR = IRSEP,IRWALL-1
c          ID = ID + 1
c          K = KORPG(1,IR)
c          RVALS(ID,1) = RVERTP(1,K)
c          RVALS(ID,2) = ZVERTP(1,K)
c        ENDDO
c        ID = ID + 1
c        K = KORPG(1,IRWALL-1)
c        RVALS(ID,1) = RVERTP(2,K)
c        RVALS(ID,2) = ZVERTP(2,K)
c        CALL GRTRAC(RVALS(1,1),RVALS(1,2),ID,NAME,'LINE',-1)
cC
c        ID = 0
c        DO IR = IRTRAP+1,NRS
c          ID = ID + 1
c          K = KORPG(NKS(IR),IR)
c          RVALS(ID,1) = RVERTP(4,K)
c          RVALS(ID,2) = ZVERTP(4,K)
c        ENDDO
c        DO IR = IRSEP,IRWALL-1
c          ID = ID + 1
c          K = KORPG(NKS(IR),IR)
c          RVALS(ID,1) = RVERTP(4,K)
c          RVALS(ID,2) = ZVERTP(4,K)
c        ENDDO
c        ID = ID + 1
c        K = KORPG(NKS(IRWALL-1),IRWALL-1)
c        RVALS(ID,1) = RVERTP(3,K)
c        RVALS(ID,2) = ZVERTP(3,K)
c        CALL GRTRAC(RVALS(1,1),RVALS(1,2),ID,NAME,'LINE',-1)
c slmod end
c
c also plot the nimbus wall and pump (if any)
c
c slmod begin
        if     (barebones) then
        elseif (nvesm.ne.0) then
c
c        if (nvesm.ne.0) then
c slmod end
c
          do i = 1,nvesm
c
c             if (nbr.GT.0.AND.jvesm(i).NE.0) CYCLE
             if (jvesm(i).ne.9) then
                wvals(1,1) = rvesm(i,1)
                wvals(1,2) = zvesm(i,1)
                wvals(2,1) = rvesm(i,2)
                wvals(2,2) = zvesm(i,2)
                call grtrac (wvals(1,1),wvals(1,2),2,
     >                       name,'LINE',-1)
             endif
          enddo
c
c          do i = 1,nvesm
c            wvals(i,1) = rvesm(i,1)
c            wvals(i,2) = zvesm(i,1)
c          enddo
c          wvals(nvesm+1,1) = rvesm(nvesm,2)
c          wvals(nvesm+1,2) = zvesm(nvesm,2)
c          call grtrac (wvals(1,1),wvals(1,2),nvesm+1,name,'LINE',-1)
c
c the pump is not a continuous nose-to-tail structure
c
          do i = 1,nvesp
            wvals(1,1) = rvesm(nvesm+i,1)
            wvals(1,2) = zvesm(nvesm+i,1)
            wvals(2,1) = rvesm(nvesm+i,2)
            wvals(2,2) = zvesm(nvesm+i,2)
            call grtrac (wvals(1,1),wvals(1,2),2,name,'LINE',-1)
          enddo
c
c       Plot the wall read from the grid file if NIMBUS wall is
c       not available
c
        elseif (nves.ne.0) then
          do i = 1,nves
            wvals(i,1) = rves(i)
            wvals(i,2) = zves(i)
          enddo
          wvals(nves+1,1) = rves(1)
          wvals(nves+1,2) = zves(1)
          call grtrac (wvals(1,1),wvals(1,2),nves+1,name,'LINE',-1)
        endif
      ENDIF
C
      RETURN
c slmod begin
 99   WRITE(0,*) 'SUPIMP2 ERROR: UNABLE TO OPEN TRIANGLES FILE'
      RETURN
c      STOP
c slmod end
      END
C
C
C
      SUBROUTINE LOSINT (TVALS,TOUTS,TWIDS,NUMTHE,ROBS,ZOBS,AVPTS,VS,
     >                   MAXDIST,intopt)
      use mod_params
      use mod_cgeom
      IMPLICIT NONE
C
C  *********************************************************************
C  *                                                                   *
C  *  LOSINT:  INTEGRATE THE VARIABLE VS ALONG A FAN OF SIGHT LINES    *
C  *           FROM A COMMON OBSERVATION POINT (ROBS,ZOBS).  SPATIAL   *
C  *           RESOLUTION IS SIMULATED BY AN AVPTS-POINT AVERAGE OF    *
C  *           A SET OF CHORDS SPANNING THE INTERVAL TWIDS.  AT THE    *
C  *           MOMENT THERE ARE NO WEIGHTS ON THIS AVERAGING WHICH     *
C  *           IMPLIES THAT WE ARE ASSUMING A RECTANGULAR, RATHER      *
C  *           THAN A CIRCULAR, VIEWING CONE.  THE INTEGRAL IS         *
C  *           PERFORMED BY FINDING THE PATH LENGTHS OF THE LOS IN     *
C  *           EACH PLASMA CELL.                                       *
c  *                                                                   *
c  *           INTOPT specifies a LOS integration that will be used    *
c  *                  for each component LOS that is used to calculate *
c  *                  the actual value for each LOS.                   *
c  *                  = 0 - weighted equally - equivalent to a         *
c  *                        rectangular view                           *
c  *                  = 1 - weighted by                                *
c  *                    sqrt(1-((theta-theta_center)/theta_width)**2)  *
c  *                    this should be the equivalent of a circular    *
c  *                    LOS.                                           *
C  *                                                                   *
C  *            CHRIS FARRELL  (HUNTERSKIL)  MARCH 1989                *
C  *            LORNE HORTON   (JET)         JULY  1993                *
C  *            DAVID ELDER    (TORONTO)     MAY   1998                *
C  *            - modifed to optionally return the MAX value in LOS    *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE "PARAMS"
c     include 'params'
C     INCLUDE "CGEOM"
c     include 'cgeom'
      INTEGER NUMTHE,AVPTS,intopt
      REAL    TVALS(NUMTHE),TOUTS(NUMTHE),TWIDS(NUMTHE),
     >        ROBS,ZOBS,VS(MAXNKS,MAXNRS),maxdist
C
      INTEGER I,J,K,IK,IR,NINT,SIDE(2)
      REAL    THETA,XB(2),WB(2),DIST(2),actdist
      real    weight_factor,total_weight,wfact
c
c     Finding location of maximum along LOS
c
      integer maxswitch,ikmax,irmax
      real maxval 
C
c     If AVPTS is passed in as < 0 - this is used to indicate that
c     the LOS routine should return the MAXIMUM value of the function
c     along the LOS instead of the integral.
c
      maxswitch = 0
c
      if (avpts.lt.0.0) then
         avpts = abs(avpts)
         maxswitch = 1
      endif
c
c      write (6,*) 'debug:',avpts,maxswitch
c
      XB(1) = ROBS
      XB(2) = ZOBS
C
C  LOOP OVER SIGHT LINES
C
      DO 200 I = 1, NUMTHE
        TVALS(I) = 0.0
c
c       Zero out location of maximum value
c
        maxval= -HI
        ikmax = 0
        irmax = 0 
        total_weight = 0.0
c
C
C  LOOP OVER CHORDS FOR AVERAGING
C
        DO 100 J = 1, AVPTS
          IF (AVPTS.EQ.1) THEN
            THETA = TOUTS(I)
          ELSE
c
c           Modify averaging angles so that the "edge" of the 
c           outside of the averaging cone matches the edge of
c           the specified viewing cone width.
c
c            THETA = TOUTS(I) - 0.5*TWIDS(I) + (J-1)*TWIDS(I)/(AVPTS-1)
c
            THETA = TOUTS(I) - 0.5*TWIDS(I) 
     >                       + (real(J)-0.5)*TWIDS(I)/AVPTS
c
          ENDIF
c
c         Calculate weight_factor for sub-chord
c
          if (intopt.eq.0) then 
             weight_factor = 1.0
          elseif (intopt.eq.1) then   
c
             if (twids(i).gt.0.0) then 
                wfact = 1.0-((theta-touts(i))/(0.5*twids(i)))**2
             else
                wfact = 1.0
             endif
c
             if (wfact.lt.0.0) then 
                weight_factor = 0.0
             else
                weight_factor = sqrt(wfact)
             endif  

          else     
             weight_factor=1.0  
          endif 
c
c         NOTE: If weight factor is zero then do not sum into LOS - cycle
c               the loop at this point
c
          if (weight_factor.le.0.0) cycle

c
c         Sum to total weight
c
          total_weight = total_weight+weight_factor
c
          THETA = THETA*DEGRAD
          WB(1) = COS(THETA)
          WB(2) = SIN(THETA)
C
C  LOOP OVER PLASMA CELLS
C
          DO 20 IR = 1, NRS

            if (ir.eq.1.or.ir.eq.irwall.or.ir.eq.irtrap) cycle

            DO 10 IK = 1, NKS(IR)
              K = KORPG(IK,IR)
c
              if (k.eq.0.or.nvertp(k).eq.0) cycle
c
              CALL INTERS(NVERTP(K),RVERTP(1,K),ZVERTP(1,K),
     >                    XB,WB,NINT,DIST,SIDE)
c
              IF (NINT.EQ.2 .AND. DIST(2).GT.0.0) THEN
c
c               NOTE: This may not be the most efficient way to
c               implement the MAX on LOS functionality. This code
c               could be rewritten to return the series of cells
c               along the LOS for the one case and the path
c               length through the cells for the other. However,
c               the interpretation of "MAX" may change so this
c               seems like the best way to implement it for
c               now.
c
                if (maxdist.le.0.0.or.
     >             (maxdist.gt.0.0.and.dist(1).lt.maxdist)) then 
c
c                  Locate and record maximum along LOS.  
c
                   if (vs(ik,ir).gt.maxval) then 
                      maxval = vs(ik,ir)
                      ikmax = ik
                      irmax = ir 
                   endif 
c
c
                   if (maxswitch.eq.0) then
c
c
c                     Assign actual distance to be used for 
c                     integration  
c 
                      if (maxdist.gt.0.0.and.
     >                    dist(2).gt.maxdist) then 
                         actdist = maxdist
                      else
                         actdist = dist(2)
                      endif 
c
                      IF (DIST(1).LT.0.0) THEN
                        TVALS(I)=TVALS(I)+VS(IK,IR)*ACTDIST
     >                                    *weight_factor
                      ELSE
                        TVALS(I)=TVALS(I)+VS(IK,IR)*(ACTDIST-DIST(1))
     >                                    *weight_factor
                      ENDIF
c
                   else
c
c                     Set TVALS to maximum of function.
c
                      tvals(i) = max (tvals(i),abs(vs(ik,ir)))
c         
                   endif  
c
                endif
c
c                write(6,'(a,3i5,5(1x,g13.5))') 'los:',i,ik,ir,
c     >                     tvals(i),vs(ik,ir),
c     >                     dist(2),dist(1),maxdist
c
              ENDIF
   10       CONTINUE
   20     CONTINUE
  100   CONTINUE
c
c       Normalize for LOS integral averaging.
c
        if (maxswitch.eq.0) then
c
c           write(6,'(a,i5,3(1x,g12.5))') 'LOS:',i,tvals(i),total_weight,
c     >                              tvals(i)/total_weight
c
           TVALS(I) = TVALS(I) / total_weight
c
        endif
c
c       Write out information about maximum on LOS 
c
        if (maxswitch.eq.1.and.ikmax.eq.0.or.irmax.eq.0) then 

           write(6,'(a,2i4,1p,8(1x,g12.5))')
     >          'NO LOCATION OF MAXIMUM EMISSION LOS :',ikmax,irmax,
     >          robs,zobs,touts(i)

        elseif (maxswitch.eq.1) then 
           write(6,'(a,2i4,1p,8(1x,g12.5))')
     >          'LOCATION OF MAXIMUM EMISSION ON LOS :',ikmax,irmax,
     >          robs,zobs,touts(i),vs(ikmax,irmax),
     >          rs(ikmax,irmax),zs(ikmax,irmax),knbs(ikmax,irmax),
     >          ktebs(ikmax,irmax)    
        endif
c
  200 CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE LOSINT_SCALE (TVALS,TOUTS,TWIDS,NUMTHE,ROBS,
     >                   ZOBS,AVPTS,VS,
     >                   MAXDIST,g_d,g_w,g_l,ifact)
      use mod_params
      use mod_cgeom
      IMPLICIT NONE
C
C  *********************************************************************
C  *                                                                   *
C  *  LOSINT_SCALE: INTEGRATE THE VARIABLE VS ALONG A FAN OF SIGHT LINES*
C  *           FROM A COMMON OBSERVATION POINT (ROBS,ZOBS).  SPATIAL   *
C  *           RESOLUTION IS SIMULATED BY AN AVPTS-POINT AVERAGE OF    *
C  *           A SET OF CHORDS SPANNING THE INTERVAL TWIDS.  AT THE    *
C  *           MOMENT THERE ARE NO WEIGHTS ON THIS AVERAGING WHICH     *
C  *           IMPLIES THAT WE ARE ASSUMING A RECTANGULAR, RATHER      *
C  *           THAN A CIRCULAR, VIEWING CONE.  THE INTEGRAL IS         *
C  *           PERFORMED BY FINDING THE PATH LENGTHS OF THE LOS IN     *
C  *           EACH PLASMA CELL.                                       *
c  * 
c  * This version of the routine will apply a scaling factor to each 
c  * LOS that contributes to the total and will then SUM them instead
c  * of averaging them. 
c  *        
c  *
C  *                                                                   *
C  *            CHRIS FARRELL  (HUNTERSKIL)  MARCH 1989                *
C  *            LORNE HORTON   (JET)         JULY  1993                *
C  *            DAVID ELDER    (TORONTO)     MAY   1998                *
C  *            - modifed to optionally return the MAX value in LOS    *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE "PARAMS"
c     include 'params'
C     INCLUDE "CGEOM"
c     include 'cgeom'
      INTEGER NUMTHE,AVPTS,ifact
      REAL    TVALS(NUMTHE),TOUTS(NUMTHE),TWIDS(NUMTHE),
     >        ROBS,ZOBS,VS(MAXNKS,MAXNRS),maxdist,
     >        g_d,g_w,g_l
c
      INTEGER I,J,K,IK,IR,NINT,SIDE(2)
      REAL    THETA,XB(2),WB(2),DIST(2),actdist
c
      real :: thetfact, dtheta, phi,scalef, response_function
c     
c     Finding location of maximum along LOS
c
      integer maxswitch,ikmax,irmax,in
      real maxval 
C
c     If AVPTS is passed in as < 0 - this is used to indicate that
c     the LOS routine should return the MAXIMUM value of the function
c     along the LOS instead of the integral.
c
      maxswitch = 0
c
      if (avpts.lt.0.0) then
         avpts = abs(avpts)
         maxswitch = 1
      endif
c
c     Find theta response function scaling factor 
c
      if (avpts.le.1) then 
         thetfact = twids(1)
      elseif (ifact.eq.7) then 

         thetfact = 0.0 

         do j = 1,avpts
c
c           Modify averaging angles so that the "edge" of the 
c           outside of the averaging cone matches the edge of
c           the specified viewing cone width.
c
c            THETA = TOUTS(1) - 0.5*TWIDS(1) 
c     >                       + (J-1)*TWIDS(1)/(AVPTS-1)
c
            THETA = TOUTS(1) - 0.5*TWIDS(1) 
     >                       + (real(J)-0.5)*TWIDS(1)/AVPTS
c
            dtheta = twids(1)/avpts * degrad
            phi = abs( theta - touts(1) ) * degrad
            response_function = 
     >           (g_w +g_d - 2.0 * g_l * tan(phi)) / (2.0 * g_d)
c
            thetfact = thetfact * response_function * dtheta
c
         end do 
     
         write(6,'(a,5(1x,g12.5))') 
     >          'THETFACT:',thetfact,thetfact/twids(1)  

       endif

c
c      write (6,*) 'debug:',avpts,maxswitch
c
      XB(1) = ROBS
      XB(2) = ZOBS
C
C  LOOP OVER SIGHT LINES
C
      DO 200 I = 1, NUMTHE
        TVALS(I) = 0.0
c
c       Zero out location of maximum value
c
        maxval= -HI
        ikmax = 0
        irmax = 0 
C
C  LOOP OVER CHORDS FOR SCALING AND SUMMING UP
C
        DO 100 J = 1, AVPTS
          IF (AVPTS.EQ.1) THEN
            THETA = TOUTS(I)
          ELSE
c
c           Modify averaging angles so that the "edge" of the 
c           outside of the averaging cone matches the edge of
c           the specified viewing cone width.
c
c            THETA = TOUTS(1) - 0.5*TWIDS(1) 
c     >                       + (J-1)*TWIDS(1)/(AVPTS-1)
c
            THETA = TOUTS(I) - 0.5*TWIDS(I) 
     >                       + (real(J)-0.5)*TWIDS(I)/AVPTS
c
          ENDIF
c
          dtheta = twids(i)/avpts * degrad
          phi = abs( theta - touts(i) ) * degrad
          response_function = 
     >           (g_w +g_d - 2.0 * g_l * tan(phi)) / (2.0 * g_d)
     
c
c         add geometry factor to functional response
c
          if (ifact.eq.5) then 
             scalef = response_function * tan(dtheta/2.0)
          elseif (ifact.eq.6) then 
             scalef = response_function * dtheta *
     >                sin( (90.0 - phi * raddeg)*degrad)
          elseif (ifact.eq.7) then 
             scalef = response_function * dtheta *
     >                sin( (90.0 - phi * raddeg)*degrad) *
     >                thetfact
          else
             scalef = 1.0
          endif 
c
          write(6,'(a,g12.5)') 'SCALEF:',scalef
c
          THETA = THETA*DEGRAD
          WB(1) = COS(THETA)
          WB(2) = SIN(THETA)
C
C  LOOP OVER PLASMA CELLS
C
          DO 20 IR = 1, NRS

            if (ir.eq.1.or.ir.eq.irwall.or.ir.eq.irtrap) cycle

            DO 10 IK = 1, NKS(IR)
              K = KORPG(IK,IR)
c
              if (k.eq.0.or.nvertp(k).eq.0) cycle
c
              CALL INTERS(NVERTP(K),RVERTP(1,K),ZVERTP(1,K),
     >                    XB,WB,NINT,DIST,SIDE)
c
              IF (NINT.EQ.2 .AND. DIST(2).GT.0.0) THEN
c
c               NOTE: This may not be the most efficient way to
c               implement the MAX on LOS functionality. This code
c               could be rewritten to return the series of cells
c               along the LOS for the one case and the path
c               length through the cells for the other. However,
c               the interpretation of "MAX" may change so this
c               seems like the best way to implement it for
c               now.
c
                if (maxdist.le.0.0.or.
     >             (maxdist.gt.0.0.and.dist(1).lt.maxdist)) then 
c
c                  Locate and record maximum along LOS.  
c
                   if (vs(ik,ir).gt.maxval) then 
                      maxval = vs(ik,ir)
                      ikmax = ik
                      irmax = ir 
                   endif 
c
c
                   if (maxswitch.eq.0) then
c
c
c                     Assign actual distance to be used for 
c                     integration  
c 
                      if (maxdist.gt.0.0.and.
     >                    dist(2).gt.maxdist) then 
                         actdist = maxdist
                      else
                         actdist = dist(2)
                      endif 
c
                      IF (DIST(1).LT.0.0) THEN
                        TVALS(I)=TVALS(I)
     >                        +VS(IK,IR)*ACTDIST*scalef
                      ELSE
                        TVALS(I)=TVALS(I)
     >                        +VS(IK,IR)*(ACTDIST-DIST(1))*scalef
                      ENDIF
c
                   else
c
c                     Set TVALS to maximum of function.
c
                      tvals(i) = max (tvals(i),abs(vs(ik,ir)))
c         
                   endif  
c
                endif
c
c                write(6,'(a,3i5,5(1x,g13.5))') 'los:',i,ik,ir,
c     >                     tvals(i),vs(ik,ir),
c     >                     dist(2),dist(1),maxdist
c
              ENDIF
   10       CONTINUE
   20     CONTINUE
  100   CONTINUE
c
c       Normalize for LOS integral averaging.
c
        if (maxswitch.eq.0.and.ifact.ne.5) then
c
           TVALS(I) = TVALS(I) / AVPTS
c
        endif
c
c       Write out information about maximum on LOS 
c
        if (maxswitch.eq.1.and.ikmax.eq.0.or.irmax.eq.0) then 

           write(6,'(a,2i4,1p,8(1x,g12.5))')
     >          'NO LOCATION OF MAXIMUM EMISSION LOS :',ikmax,irmax,
     >          robs,zobs,touts(i)

        elseif (maxswitch.eq.1) then 
           write(6,'(a,2i4,1p,8(1x,g12.5))')
     >          'LOCATION OF MAXIMUM EMISSION ON LOS :',ikmax,irmax,
     >          robs,zobs,touts(i),vs(ikmax,irmax),
     >          rs(ikmax,irmax),zs(ikmax,irmax),knbs(ikmax,irmax),
     >          ktebs(ikmax,irmax)    
        endif
c
  200 CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE LOSINTEXPT (TVALS,TOUTS,TWIDS,NUMTHE,ROBS,ZOBS,AVPTS,
     >                   MAXDIST,
     >                   maxix,maxiy,nix,niy,expt_array,raxis,zaxis)
      use mod_params
      IMPLICIT NONE
C
C  *********************************************************************
C  *                                                                   *
C  *  LOSINTRZ:INTEGRATE THE 2D EXPERIMENTAL DATA                      *
C  *           ALONG A FAN OF SIGHT LINES                              *
C  *           FROM A COMMON OBSERVATION POINT (ROBS,ZOBS).  SPATIAL   *
C  *           RESOLUTION IS SIMULATED BY AN AVPTS-POINT AVERAGE OF    *
C  *           A SET OF CHORDS SPANNING THE INTERVAL TWIDS.  AT THE    *
C  *           MOMENT THERE ARE NO WEIGHTS ON THIS AVERAGING WHICH     *
C  *           IMPLIES THAT WE ARE ASSUMING A RECTANGULAR, RATHER      *
C  *           THAN A CIRCULAR, VIEWING CONE.  THE INTEGRAL IS         *
C  *           PERFORMED BY FINDING THE PATH LENGTHS OF THE LOS IN     *
C  *           EACH PLASMA CELL.                                       *
C  *                                                                   *
C  *            CHRIS FARRELL  (HUNTERSKIL)  MARCH 1989                *
C  *            LORNE HORTON   (JET)         JULY  1993                *
C  *            DAVID ELDER    (TORONTO)     MAY   1998                *
C  *            - modifed to optionally return the MAX value in LOS    *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE "PARAMS"
c     include 'params'
      INTEGER NUMTHE,AVPTS
      REAL    TVALS(NUMTHE),TOUTS(NUMTHE),TWIDS(NUMTHE),
     >        ROBS,ZOBS,VS(MAXNKS,MAXNRS),maxdist
C
      INTEGER I,J,K,IX,IY,NINT,SIDE(2)
      REAL    THETA,XB(2),WB(2),DIST(2),actdist
c
c     2D experimental data arrays and characteristics
c
      integer  maxix,maxiy,nix,niy,ifnopt
      real expt_array(maxix,maxiy),raxis(nix),zaxis(niy)
c
c     Define cell corners for INTERS routine
c
      real r1,r2,z1,z2 
      real rvert(4),zvert(4)
      integer nvert
c
c     Maximum along LOS functions:
c
      real    maxval  
      integer maxswitch,ixmax,iymax
C
c     If AVPTS is passed in as < 0 - this is used to indicate that
c     the LOS routine should return the MAXIMUM value of the function
c     along the LOS instead of the integral.
c
      maxswitch = 0
c
      if (avpts.lt.0.0) then
         avpts = abs(avpts)
         maxswitch = 1
      endif
c
c      write (6,*) 'debug:',avpts,maxswitch,nix,niy
c
      XB(1) = ROBS
      XB(2) = ZOBS
C
C  LOOP OVER SIGHT LINES
C
      DO 200 I = 1, NUMTHE
        TVALS(I) = 0.0
c
c       Zero out location of maximum value
c
        maxval= -HI
        ixmax = 0
        iymax = 0 
C
C  LOOP OVER CHORDS FOR AVERAGING
C
        DO 100 J = 1, AVPTS
          IF (AVPTS.EQ.1) THEN
            THETA = TOUTS(I)
          ELSE
            THETA = TOUTS(I) - 0.5*TWIDS(I) + (J-1)*TWIDS(I)/(AVPTS-1)
          ENDIF
c
          THETA = THETA*DEGRAD
          WB(1) = COS(THETA)
          WB(2) = SIN(THETA)
C
C  LOOP OVER EXPERIMENTAL CELLS
C
          DO 20 IX = 1, NIX

            DO 10 IY = 1, NIY
c
c             Set up cell for INTERS routine 
c
              nvert = 4
c
c             Calculate R,Z vertices around the centre of the ix,iy cell.
c
c             R
c
              if (ix.eq.1) then 
c                
                 r2 = (raxis(ix+1)+raxis(ix))/2.0
                 r1 = raxis(ix) + (raxis(ix)-r2)
c
              elseif (ix.eq.nix) then 
c
                 r1 = (raxis(ix)+raxis(ix-1))/2.0
                 r2 = raxis(ix) + (raxis(ix)-r1)
c
              else
c
                 r1 = (raxis(ix)+raxis(ix-1))/2.0
                 r2 = (raxis(ix+1)+raxis(ix))/2.0
c
              endif 
c
c             Z
c
              if (iy.eq.1) then 
c                
                 z2 = (zaxis(iy+1)+zaxis(iy))/2.0
                 z1 = zaxis(iy) + (zaxis(iy)-z2) 
c
              elseif (iy.eq.niy) then 
c
                 z1 = (zaxis(iy)+zaxis(iy-1))/2.0
                 z2 = zaxis(iy) + (zaxis(iy)-z1)
c
              else
c
                 z1 = (zaxis(iy)+zaxis(iy-1))/2.0
                 z2 = (zaxis(iy+1)+zaxis(iy))/2.0
c
              endif 
c
c             Combine these values in the correct order to define the
c             vertices of the polygon.
c
              rvert(1) = r1
              zvert(1) = z1 
c
              rvert(2) = r1 
              zvert(2) = z2
c
              rvert(3) = r2 
              zvert(3) = z2
c
              rvert(4) = r2 
              zvert(4) = z1
c
c             calculate cell intersection
c
              CALL INTERS(NVERT,RVERT,ZVERT,
     >                    XB,WB,NINT,DIST,SIDE)
c
              IF (NINT.EQ.2 .AND. DIST(2).GT.0.0) THEN
c
c               NOTE: This may not be the most efficient way to
c               implement the MAX on LOS functionality. This code
c               could be rewritten to return the series of cells
c               along the LOS for the one case and the path
c               length through the cells for the other. However,
c               the interpretation of "MAX" may change so this
c               seems like the best way to implement it for
c               now.
c


                if (maxdist.le.0.0.or.
     >             (maxdist.gt.0.0.and.dist(1).lt.maxdist)) then 
c
c                  Locate and record maximum along LOS.  
c
                   if (expt_array(ix,iy).gt.maxval) then 
                      maxval = expt_array(ix,iy)
                      ixmax = ix
                      iymax = iy 
                   endif 
c
                   if (maxswitch.eq.0) then
c
c                     Assign actual distance to be used for 
c                     integration  
c 
                      if (maxdist.gt.0.0.and.
     >                    dist(2).gt.maxdist) then 
                         actdist = maxdist
                      else
                         actdist = dist(2)
                      endif 
c
                      IF (DIST(1).LT.0.0) THEN
                        TVALS(I)=TVALS(I)+expt_array(ix,iy)*ACTDIST
                      ELSE
                        TVALS(I)=TVALS(I)
     >                          +expt_array(ix,iy)*(ACTDIST-DIST(1))
                      ENDIF
c
                   elseif (maxswitch.eq.1) then 
c
c                     Set TVALS to maximum of function.
c
                      tvals(i) = max (tvals(i),expt_array(ix,iy))
c         
                   endif  
c
                endif
c
c                write(6,'(a,4i5,8(1x,g13.5))') 'rzl:',i,j,ix,iy,
c     >                     tvals(i),expt_array(ix,iy),
c     >                     dist(2),dist(1),maxdist
c
              ENDIF
   10       CONTINUE
   20     CONTINUE
  100   CONTINUE
c
c       Normalize for LOS integral averaging.
c
        if (maxswitch.eq.0) then
c
           TVALS(I) = TVALS(I) / AVPTS
c
        endif
c
c       Write out information about maximum on LOS 
c
        if (ixmax.gt.0.and.iymax.gt.0) then  
           write(6,'(a,2i4,1p,8(1x,g12.5))')
     >          'LOCATION OF MAXIMUM EMISSION OF EXPT:',ixmax,
     >          iymax,
     >          robs,zobs,touts(i),expt_array(ixmax,iymax),
     >          raxis(ixmax),zaxis(iymax)    
        else
           write(6,'(a,2i4,1p,8(1x,g12.5))')
     >          'WARNING: IXMAX,IYMAX = 0 ON LOS :',ixmax,iymax,
     >          robs,zobs,touts(i)

        endif 

  200 CONTINUE
C
      RETURN
      END
c
c
c
      SUBROUTINE LOS3DINT(TVAL,init_position,dircos,wres,vs,nchords,
     >              step_size,minrun_invessel,contrib,reflect_opt)
      use mod_params
      use mod_cgeom
      use mod_grbound
      IMPLICIT NONE
c     include 'params'
c     include 'cgeom'
c     include 'grbound'
c
      integer nchords,minrun_invessel,reflect_opt
      real*8 tval,init_position(3)
      real*8 wres(nchords),dircos(nchords,3)
      real*8 step_size
      real vs(maxnks,maxnrs)
      real contrib(maxnks,maxnrs)
C
C  *********************************************************************
C  *                                                                   *
C  *  LOS3DINT :  INTEGRATE THE VARIABLE VS ALONG THE GIVEN SIGHT LINE    *
C  *           FROM A COMMON OBSERVATION POINT (X,Y,Z).  SPATIAL       *
C  *           RESOLUTION IS SIMULATED BY NCHORDS INDIVIDUAL SIGHT     *
C  *           LINES FOR THE GIVEN VIEW - EACH WITH A WEIGHTING WRES   *
C  *                                                                   *
C  *********************************************************************
C
      integer i,ik,ir,invessel,num_reflections
      real resultvw,resultc,resultg

      real maxrun,maxrun_length
      parameter (maxrun_length=100.0)

C
      real*8 mult_fact,couts(3),step,position(3)
c
      real run,x,y,z,r,run_entered,run_reflected
      real ang,angt,angm,rv1,rv2,zv1,zv2
      real ctot
      real lastr,lastz 
      logical r_increasing,z_increasing
      integer ierr
c
      logical newinj,outofgrid,on_grid
C
c     Initialization
c
c
      newinj = .false.
      outofgrid = .false.
      call rzero(contrib,maxnks*maxnrs)
c
c     Set maximum number of steps along LOS
c
      if (step_size.gt.0.0) then 
         maxrun = maxrun_length/step_size
      else
c
c        If step_size is less than or equal to zero - issue error message and 
c        return
c
         write(6,'(a,1x,g16.8)') 'LOS integration step_size'//
     >       ' is less than or equal to zero in LOS3DINT:',step_size  
         write(0,'(a,1x,g16.8)') 'LOS integration step_size'//
     >       ' is less than or equal to zero in LOS3DINT:',step_size  
         return
      endif  
C
C  LOOP OVER INDIVIDUAL VIEWING CHORDS
C
      TVAL = 0.0
c
      DO I = 1, NCHORDS
c
c        Set number of reflections to zero  
c
         num_reflections = 0
         mult_fact = 1.0d0
c
         run=0.0
         run_entered = 0.0
         run_reflected = 0.0
c
         invessel = 0 
c
         couts(1) = dircos(i,1)
         couts(2) = dircos(i,2)
         couts(3) = dircos(i,3)
c
         position(1) = init_position(1) 
         position(2) = init_position(2) 
         position(3) = init_position(3) 
c
         step = step_size
c 
c        Normalize COUTS - just in case
c
         ctot = sqrt(couts(1)**2+couts(2)**2+couts(3)**2)
c
         on_grid = .false. 
c
         if (ctot.ne.0.0) then 
            couts(1) = couts(1)/ctot
            couts(2) = couts(2)/ctot
            couts(3) = couts(3)/ctot
         else
            exit
         endif  
c
         lastr = sqrt(position(1)**2+position(2)**2)
         lastz = position(3)
c
c
c  Begin moving along sight line...
c
c  Add test so that if R > RMAX or Z < ZMIN or Z > ZMAX and moving away - stop LOS
c 
c
         do while (run.lt.maxrun.and.invessel.lt.2) 

            run=run+1.0
            x=position(1)+((run-run_reflected)-0.5)*step*couts(1)
            y=position(2)+((run-run_reflected)-0.5)*step*couts(2)
            z=position(3)+((run-run_reflected)-0.5)*step*couts(3)
            r= sqrt(x**2 + y**2)
c
            if ((z.gt.zmax.and.z.gt.lastz).or.
     >          (z.lt.zmin.and.z.lt.lastz).or.
     >          (r.gt.rmax.and.r.gt.lastr)) then
c
c               LOS is moving away from the vessel - stop following it.
c
c                write(6,'(a,1x,f12.1))') 'STOPPED FOLLOWING LOS:',run
c                write(6,'(a,l4,4(1x,g12.5))') 'ZMAX:',
c     >                  z.gt.zmax.and.z.gt.lastz,z,zmax,lastz    
c                write(6,'(a,l4,4(1x,g12.5))') 'ZMIN:',
c     >                 z.lt.zmin.and.z.lt.lastz,z,zmin,lastz    
c                write(6,'(a,l4,4(1x,g12.5))') 'RMAX:',
c     >                r.gt.rmax.and.r.gt.lastr,r,rmax,lastr    
c
                exit
c            
            endif
c
c
c           If LOS is on_grid then skip boundary checks. 
c
C
            if (on_grid) then 

               outofgrid = .false.
c
c              Find grid cell and add contribution
c
               call gridpos(ik,ir,r,z,newinj,outofgrid)

               if (.not.outofgrid) then
                  TVAL = TVAL + VS(IK,IR)*step*wres(i)*mult_fact
                  contrib(ik,ir) = contrib(ik,ir) + 1.0
               else 
c
                  on_grid = .false.
c
c                 Check to see if LOS is outside vessel wall and
c                 do reflection if necessary - can occur at targets. 
c
                  CALL GA15B(R,Z,RESULTVW,neutwpts,1,
     >                       NWWORK,4*MAXPTS,
     >                       NWINDW,MAXPTS,RNW,ZNW,NWTDUM,
     >                       NWXDUM,NWYDUM,6)
c
                  if (resultvw.lt.0.0.and.invessel.eq.1) then 
c
                     if (reflect_opt.ne.0) then 
c
                        num_reflections = num_reflections+1 
c 
c                       EXIT vessel if more than 1 reflection
c
c                       Note: minrun option is not compatible with
c                       reflection and that input is ignored.
c
C
                        if (num_reflections.gt.1) then 

                           invessel = 2
c
c                       Calculate reflection of vector and continue
c                       summimg up. 
c
                        else 
c
                           ierr = 0
                           run_reflected = run
c
c                           write(6,'(a,8(1x,f15.5))') 'REFLECT:',
c     >                              run,run_entered,x,y,z,r
c                           write(6,'(a,8(1x,f15.5))') 'COUTSIN:',
c     >                             couts(1),couts(2),couts(3)
c                       
c
                           call calc_reflection(x,y,z,r,couts,
     >                              step,mult_fact,reflect_opt,
     >                              ierr)
c
c                          Assign reflection point as new observation
c                          position for this LOS. 
c
                           position(1) = x
                           position(2) = y
                           position(3) = z
c
c
c                          Error occurred in reflection - stop LOS and issue
c                          error message.
c
                           if (ierr.ne.0) then 
                              invessel = 2
                              write(6,
     >            '(''ERROR in 3D LOS reflection:'',i5)')  
     >                          ierr                             
                              ierr = 0 
                           endif

                        endif
c
                     else 
c
                        if ((run-run_entered).gt.minrun_invessel)
     >                       then 
                           invessel = 2
c
c                        write(6,'(a,3(1x,g12.5))') 'EXITING VESSEL:',
c     >                         run,r,z 
c
                        endif 
c
                     endif
c
                  endif  
c
               endif
c          
            else
c
c              LOS position was outside some boundary - check
c              if it has entered the vessel. 
c
               CALL GA15B(R,Z,RESULTVW,neutwpts,1,NWWORK,4*MAXPTS,
     >            NWINDW,MAXPTS,RNW,ZNW,NWTDUM,NWXDUM,NWYDUM,6)
c
               IF (RESULTVW.ge.0.0) THEN
c
                  if (invessel.eq.0) then 
                     run_entered = run
                     invessel = 1
                  endif
c
c                 Check if inside grid
c
                  CALL GA15B(R,Z,RESULTG,
     >               IONWPTS,1,IWWORK,4*MAXPTS,
     >               IWINDW,MAXPTS,RIW,ZIW,
     >               IWTDUM,IWXDUM,IWYDUM,6)
c
                  if (resultg.ge.0.0) then  
c
c                    Check if outside core
c
                     CALL GA15B(R,Z,RESULTC,
     >                  IONCPTS,1,ICWORK,4*MAXPTS,
     >                  ICINDW,MAXPTS,RCW,ZCW,
     >                  ICTDUM,ICXDUM,ICYDUM,6)
c
                     if (resultc.lt.0.0) then  

                        outofgrid = .false.
c
c                       Find grid cell and add contribution
c
                        call gridpos(ik,ir,r,z,newinj,outofgrid)

                        if (.not.outofgrid) then
                           on_grid = .true.
                           TVAL = TVAL + VS(IK,IR)*step*wres(i)
     >                                   *mult_fact
                           contrib(ik,ir) = contrib(ik,ir)+1.0
                        else
                           on_grid = .false. 
                        endif
c
                     else
                        on_grid = .false.          
                     endif
c
                  else
                     on_grid = .false.    
                  endif
c
               elseif (resultvw.lt.0.0.and.invessel.eq.1) then 
c
                  on_grid = .false. 
c
                  if (reflect_opt.ne.0) then 
c
                     num_reflections = num_reflections+1 
c 
c                    EXIT vessel if more than 1 reflection
c
c                    Note: minrun option is not compatible with
c                    reflection and that input is ignored.
c
C
                     if (num_reflections.gt.1) then 

                        invessel = 2
c
c                    Calculate reflection of vector and continue summimg up. 
c
                     else 
c
                        ierr = 0
                        run_reflected = run
c
c                        write(6,'(a,8(1x,f15.5))') 'REFLECT:',
c     >                              run,run_entered,x,y,z,r
c                        write(6,'(a,8(1x,f15.5))') 'COUTSIN:',
c     >                             couts(1),couts(2),couts(3)
c                       
c
                        call calc_reflection(x,y,z,r,couts,
     >                              step,mult_fact,reflect_opt,
     >                              ierr)
c
c                       Assign reflection point as new observation
c                       position for this LOS. 
c
                        position(1) = x
                        position(2) = y
                        position(3) = z
c
c
c                       Error occurred in reflection - stop LOS and issue
c                       error message.
c
                        if (ierr.ne.0) then 
                           invessel = 2
                           write(6,
     >                      '(''ERROR in 3D LOS reflection:'',i5)')  
     >                       ierr                             
                           ierr = 0 
                        endif

                     endif
c
                  else 
c
                     if ((run-run_entered).gt.minrun_invessel) then 
                        invessel = 2
c
c                     write(6,'(a,3(1x,g12.5))') 'EXITING VESSEL:',
c     >                         run,r,z 
c
                     endif 
c
                  endif
c
               endif
c
            endif
c
c            write(6,'(a,l4,3i4,9(1x,g12.5))')  'LOS3DINT',
c     >           on_grid,ik,ir,invessel,run,resultvw,resultg,resultc,
c     >           x,y,z,r,tval   
c
         end do 
c
      end do
c
c     Renormalize result for integration over N chords 
c
      tval = tval/nchords
c
c      write(6,'(a,3(1x,f16.2))') 'RUN:',run,run_entered,run_reflected 
c
      RETURN
      END
c
c
c
      subroutine calc_reflection(x,y,z,r,couts,
     >                        step,mult_fact,reflect_opt,ierr)
      use mod_params
      use mod_cgeom
      use mod_grbound
      implicit none
      integer reflect_opt,ierr
      real x,y,z,r
      real*8 couts(3),mult_fact,step
c
c     include 'params'
c     include 'cgeom'
c     include 'grbound'
c
c     CALC_REFLECTION: This routine calculates the reflection of the
c                      vector defined by couts when it strikes the 
c                      vessel wall. 
c
c
c     X,Y,Z,R - location that is just outside the boundary. Previous step
c               was inside. 
c
c     COUTS - incoming vector     
c     
c     Need surface normal at point of impact - use this to calculate
c     outgoing vector assuming specular reflection. Call a routine -
c     given angle of incidence that will return a reflection coefficient. 
c
c     Local variables
c
      real xnew,ynew,znew,resultw
      real*8 normvect(3),negvect(3),refvect(3),xvect(3),
     >       doti,in_angle,ref_angle
      integer res,xprod,in
      real*8 dotprod,reflection_coefficient
      external dotprod,xprod,reflection_coefficient
c
c     Initialize
c
      ierr = 0
c
c     REFINE Intersection point and find normal vector to wall
c
      call find_reflection(x,y,z,step,couts,
     >                     xnew,ynew,znew,normvect,ierr)
c
c      write(6,'(a,8(1x,f15.5))') 'REFPOINT:',xnew,ynew,znew
c      write(6,'(a,8(1x,f15.5))') 'COUTS   :',
c     >                             couts(1),couts(2),couts(3)
c
      if (ierr.ne.0) return 
c
c     Calculate angle of incidence - this is the angle between the 
c     surface normal and the negative of the incoming vector. 
c
      do in = 1,3
         negvect(in) = -couts(in)
      end do 
c
      doti = dotprod(normvect,negvect)
      in_angle = acos(doti)
      ref_angle = 2.0d0 * in_angle      
c
      res = xprod(negvect,normvect,xvect) 
c
c      write(6,'(a,8(1x,f15.5))') 'ANGLES :',
c     >            in_angle*raddeg,ref_angle*raddeg
c
c     Solve for reflection vector  
c
      call solv_vect(normvect,negvect,xvect,in_angle,
     >               ref_angle,refvect,ierr)
c
c     Exit in case of error 
c
      if (ierr.ne.0) return
c
c     Print tests of reflection vector
c
c      write(6,'(a,8(1x,g12.5))') 'REF_VECT:',refvect(1),
c     >     refvect(2),refvect(3),
c     >         acos(dotprod(refvect,normvect))*raddeg,
c     >         acos(dotprod(refvect,negvect))*raddeg,
c     >         acos(dotprod(refvect,xvect))*raddeg
c
c     Assign reflected vector to couts 
c
      do in = 1,3
         couts(in) = refvect(in)
      end do
      
c
c     Calculate reflection coefficient 
c
      mult_fact = reflection_coefficient(in_angle,reflect_opt) 
c
c     Set current location to point of reflection.
c
      x = xnew
      y = ynew
      z = znew
      r = sqrt(x**2+y**2) 
c
c     X,Y,Z should be inside the vessel wall at this point - unless the
c     reflection + partial step caries the LOS outside the vessel again -
c     in which case the LOS ends here. 
c
c      CALL GA15B(R,Z,RESULTW,neutwpts,1,NWWORK,4*MAXPTS,
c     >            NWINDW,MAXPTS,RNW,ZNW,NWTDUM,NWXDUM,NWYDUM,6)
c
c     Verify position
c
c      if (resultw.lt.0.0) then 
c
c         write(6,'(a,4(x,g16.8))') 
c     >         'ERROR: Reflection point outside wall:',x,y,z,r
c
c         ierr = 1
c
c      endif    
c
c     Return 
c
      return 
      end 
c
c
c
      real*8 function reflection_coefficient(in_angle,reflect_opt)
      implicit none
      real*8 in_angle
      integer reflect_opt
c
c     REFLECTION COEFFICIENT - assign a fixed value for now 
c
c      
      reflection_coefficient = 0.0d0
c
      if (reflect_opt.eq.1) then 
         reflection_coefficient = 0.1d0
      elseif (reflect_opt.eq.2) then 
         reflection_coefficient = 1.0d0
      endif
c
      return
      end
c
c
c
      subroutine find_reflection(x,y,z,step,couts,
     >                           xnew,ynew,znew,normvect,ierr)
      use mod_params
      use mod_comtor
      use mod_grbound
      implicit none
      real x,y,z
      real xnew,ynew,znew
      integer ierr
      real*8 step,couts(3),normvect(3)
c     include 'params'
c     include 'comtor'
c     include 'grbound'
c
c
c     Find intersection with wall for LOS.
c
c     
c     Local variables
c
      real xs,ys,zs,rs,xe,ye,ze,re
      real resultvw
      integer in
c
      real atan2c
      external atan2c 
c
      real xvect,yvect,zvect,rvect,normangle,projangle
c
c 
c     Variables required for neut wall intersection code
c
      real best,dsq,resulta
c slmod begin
c...NRFOPT is already declared in COMTOR (reported by AIX compiler).  NRFOPT
c   changed to NRFOPT2 in the rest of the routine as well.
      integer ind,indi,id,nrfopt2
c
c      integer ind,indi,id,nrfopt
c slmod end
      logical sect  
      real rnewtmp,znewtmp,tnew,tnorm
c 
c     The intersection point lies between location x,y,z and the previous
c     point along the LOS - one step backward along the LOS.
c       
      xs = x - step * couts(1)
      ys = y - step * couts(2)
      zs = z - step * couts(3)
      rs = sqrt(xs**2+ys**2)
c
c     Make sure that rs,zs is inside the vessel wall. 
c
      CALL GA15B(Rs,Zs,RESULTVW,neutwpts,1,NWWORK,4*MAXPTS,
     >            NWINDW,MAXPTS,RNW,ZNW,NWTDUM,NWXDUM,NWYDUM,6)
c
c     Move one more step back - steps back along the LOS should be OK.
c
      if (resultvw.lt.0.0) then   
         xs = xs - step * couts(1)
         ys = ys - step * couts(2)
         zs = zs - step * couts(3)
         rs = sqrt(xs**2+ys**2)
         write(6,'(a,6(1x,g16.8))') 
     >               'RS,ZS MOVED IN:',rs,zs,resultvw
         CALL GA15B(Rs,Zs,RESULTVW,neutwpts,1,NWWORK,4*MAXPTS,
     >            NWINDW,MAXPTS,RNW,ZNW,NWTDUM,NWXDUM,NWYDUM,6)
         write (6,'(a,6(1x,g16.8))') 
     >          'CHECK IN:',rs,zs,resultvw
      endif  
c
      xe = x 
      ye = y 
      ze = z 
      re = sqrt(xe**2+ye**2)
c
c     Make sure that re,ze is outside the vessel wall. 
c

      CALL GA15B(Re,Ze,RESULTVW,neutwpts,1,NWWORK,4*MAXPTS,
     >            NWINDW,MAXPTS,RNW,ZNW,NWTDUM,NWXDUM,NWYDUM,6)
c
c     Move a fraction of a step out
c
      if (resultvw.ge.0.0) then   
         xe = x + 0.1 * step * couts(1)
         ye = y + 0.1 * step * couts(2)
         ze = z + 0.1 * step * couts(3)
         re = sqrt(xe**2+ye**2)
         write(6,'(a,6(1x,g16.8))') 
     >               'RE,ZE MOVED OUT:',re,ze,resultvw

         CALL GA15B(Re,Ze,RESULTVW,neutwpts,1,NWWORK,4*MAXPTS,
     >            NWINDW,MAXPTS,RNW,ZNW,NWTDUM,NWXDUM,NWYDUM,6)
         write (6,'(a,6(1x,g16.8))') 'CHECK OUT:',re,ze,resultvw
c
      endif  
c
c
c     Find likely closest element of the wall 
c
c     The following code is taken from the wall intersection code in neut.d6a
c

c
c     Set the reflection type to specular for now
c
      nrfopt2 = 1
c
      BEST = HI
      DSQ  = HI
      IND = 1
      DO 650 ID = 1,WALLPTS
         DSQ = (WALLPT(ID,1)-RE) ** 2 + (WALLPT(ID,2)-ZE) ** 2
         IF (DSQ.LT.BEST) THEN
           BEST = DSQ
           IND   = ID
         ENDIF
650   CONTINUE
C
c         WRITE(6,'(a,2i5,10(1x,g14.8))') 
c     >              'DSQ1:',ind,WALLPTS,DSQ,R,Z,ROLD,ZOLD
c
c         WRITE(6,'(a,i5,10(1x,g14.8))') 
c     >    'DSQ2:',ind,WALLPT(ind,20),wallpt(ind,21),wallpt(ind,1),
c     >                wallpt(ind,2),wallpt(ind,22),wallpt(ind,23)
C
      CALL INTCALC(RE,ZE,RS,ZS,WALLPT(IND,1),WALLPT(IND,2),
     >             WALLPT(IND,8),WALLPT(IND,9),WALLPT(IND,5),
     >             WALLPT(IND,6),RNEWtmp,ZNEWtmp,TNEW,TNORM,
     >             SECT,nrfopt2)
      INDI = IND
C
c          WRITE(6,'(a,l4,8(1x,g13.6))') 'SECT:',SECT,
c     >              R,Z,ROLD,ZOLD,RNEW,ZNEW,orgr,orgz
C
      IF (.NOT.SECT) THEN
         DO 660 ID = 1,WALLPTS/2
           INDI = IND + ID
C      
c           WRITE(6,*) 'INDI1:' ,INDI,SECT
C      
           IF (INDI.GT.WALLPTS) INDI = INDI-WALLPTS
c
           CALL INTCALC(RE,ZE,RS,ZS,WALLPT(INDI,1),WALLPT(INDI,2),
     >             WALLPT(INDI,8),WALLPT(INDI,9),WALLPT(INDI,5),
     >             WALLPT(INDI,6),RNEWtmp,ZNEWtmp,TNEW,TNORM,
     >             SECT,nrfopt2)
c

           IF (SECT) GOTO 670
           INDI = IND - ID
C      
c           WRITE(6,*) 'INDI2:' ,INDI,SECT
C      
           IF (INDI.LT.1) INDI = WALLPTS-INDI
c
           CALL INTCALC(RE,ZE,RS,ZS,WALLPT(INDI,1),WALLPT(INDI,2),
     >             WALLPT(INDI,8),WALLPT(INDI,9),WALLPT(INDI,5),
     >             WALLPT(INDI,6),RNEWtmp,ZNEWtmp,TNEW,TNORM,
     >             SECT,nrfopt2)
           IF (SECT) GOTO 670
660      CONTINUE
      ENDIF
670   CONTINUE

c
c     Check for lack of wall intersection
c
      if (.not.sect) then 
c     
         ierr = 1  
         write (6,
     >      '(''ERROR in 3D LOS Reflection: Wall element not found'',
     >       2g16.8)') rnewtmp,znewtmp 
         CALL GA15B(Rs,Zs,RESULTVW,neutwpts,1,NWWORK,4*MAXPTS,
     >            NWINDW,MAXPTS,RNW,ZNW,NWTDUM,NWXDUM,NWYDUM,6)
         write (6,'(a,6(1x,g13.6))') 
     >          'RS,ZS,RESULT:',rs,zs,resultvw
         CALL GA15B(Re,Ze,RESULTVW,neutwpts,1,NWWORK,4*MAXPTS,
     >            NWINDW,MAXPTS,RNW,ZNW,NWTDUM,NWXDUM,NWYDUM,6)
         write (6,'(a,6(1x,g13.6))') 
     >          'RE,ZE,RESULT:',re,ze,resultvw

         return
c
      endif

c
c     Verify RNEW,ZNEW 
c
      if ( (abs(re-rnewtmp)+abs(rs-rnewtmp)).ne.abs(re-rs).or.
     >     (abs(ze-znewtmp)+abs(zs-znewtmp)).ne.abs(ze-zs)) then 
      
         write (6,'(a,i5,l4,6(1x,g13.6))') 
     >       'POSSIBLE RNEW,ZNEW ERROR:',indi,sect,re,ze,rnewtmp,
     >                  znewtmp,rs,zs

         CALL GA15B(Rs,Zs,RESULTVW,neutwpts,1,NWWORK,4*MAXPTS,
     >            NWINDW,MAXPTS,RNW,ZNW,NWTDUM,NWXDUM,NWYDUM,6)
         write (6,'(a,6(1x,g13.6))') 
     >       'RS,ZS,RESULT:',rs,zs,resultvw
         CALL GA15B(Re,Ze,RESULTVW,neutwpts,1,NWWORK,4*MAXPTS,
     >            NWINDW,MAXPTS,RNW,ZNW,NWTDUM,NWXDUM,NWYDUM,6)
         write (6,'(a,6(1x,g13.6))') 
     >       'RE,ZE,RESULT:',re,ze,resultvw
         CALL GA15B(Rnewtmp,Znewtmp,RESULTVW,neutwpts,1,NWWORK,
     >            4*MAXPTS,
     >            NWINDW,MAXPTS,RNW,ZNW,NWTDUM,NWXDUM,NWYDUM,6)
         write (6,'(a,6(1x,g13.6))') 
     >       'RN,ZN,RESULT:',rnewtmp,znewtmp,resultvw
         ierr = 2

         return      
c 
      endif
c
c
c     Wall intersection found at rnewtmp,znewtmp - calculate 3D intersection
c     point as well as the surface normal vector 
c
c     Intersection point
c
      xnew = (rnewtmp-rs)/(re-rs) * (xe-xs) + xs 
      ynew = (rnewtmp-rs)/(re-rs) * (ye-ys) + ys 
      znew = znewtmp
c
c     3-space normal vector
c
      normangle = tnorm
c
      rvect = cos(normangle)
      zvect = sin(normangle)
c
      projangle = atan2c(ynew,xnew)
c
      xvect = rvect * cos(projangle)
      yvect = rvect * sin(projangle)
c
      normvect(1) = xvect
      normvect(2) = yvect
      normvect(3) = zvect
c
      call norm_vect(normvect)
c
c      write(6,'(a,i4,8(1x,f15.5))') 'NORMVECT:',
c     >         indi,(normvect(in),in=1,3),
c     >         normangle*raddeg,projangle*raddeg
c
      return 
      end  
c
c
c
      logical function inplasma(r,z) 
      use mod_params
      use mod_cgeom
      implicit none
c
c
c     include 'params'
C     INCLUDE "CGEOM"
c     include 'cgeom'
c      include 'PPPARA'
c      include 'PPUNIT'
c      include 'PPGEOM'
C
      INTEGER IK,IR,IERR
      REAL    R,Z
C
      INTEGER I,J,K
      CHARACTER LEQUIL*72
      integer lshot
      LOGICAL INCELLA,SETUP
      DATA    SETUP/.TRUE./,LEQUIL/' '/,lshot/0/
C
      INTEGER   IW,LND
      PARAMETER (LND = 128, IW = 4*LND)
      INTEGER   NP,NE,KIND,INDP(2,LND),INDE(2,LND),LP
      REAL      WORKP(IW),WORKE(IW),
     >          XP(2*(MAXNKS+MAXNRS)+1),YP(2*(MAXNKS+MAXNRS)+1),
     >          XE(MAXNKS),YE(MAXNKS),
     >          T,XD,YD,RESULT
      DATA      KIND/1/,LP/-1/
c slmod begin
c...AIX compiler complained that INPLASMA was not assigned
c   in the routine.
      inplasma=.FALSE.
c slmod end
C
C  WHEN TREATING MULTIPLE FILES FROM THE D.A.M.U. POST PROCESSOR
C  THE EQUILIBRIUM MAY CHANGE.  CHECK FOR THIS AND REINITIALISE, IF
C  NECESSARY.
C
c      IF (lequil.NE.EQUIL) THEN
c        SETUP = .TRUE.
c        LEQUIL = EQUIL
c      ENDIF
      IF (lshot.NE.ishot) THEN
        SETUP = .TRUE.
        Lshot = ishot
      ENDIF
C
C  FIRST TIME IN, SETUP POLYGONS DEFINING PLASMA BOUNDARY AND
C  ESCAPE REGION
C
      IF (SETUP) THEN
C
C  PLASMA BOUNDARY
C
        NP = 0
        J = IRWALL - 1
        DO I = 1,NKS(J)
          K = KORPG(I,J)
C  CHECK TO SEE IF VIRTUAL POINTS EXIST AT ENDS OF RINGS!
          IF (K.NE.0) THEN
            NP = NP + 1
            XP(NP) = RVERTP(2,K)
            YP(NP) = ZVERTP(2,K)
          ENDIF
        ENDDO
        DO J = IRWALL-1,IRSEP,-1
          NP = NP + 1
          I = NKS(J)
          K = KORPG(I,J)
          IF (K.EQ.0) K = KORPG(I-1,J)
          XP(NP) = RVERTP(3,K)
          YP(NP) = ZVERTP(3,K)
        ENDDO
        DO J = NRS,IRTRAP+1,-1
          NP = NP + 1
          I = NKS(J)
          K = KORPG(I,J)
          IF (K.EQ.0) K = KORPG(I-1,J)
          XP(NP) = RVERTP(3,K)
          YP(NP) = ZVERTP(3,K)
        ENDDO
        J = IRTRAP + 1
        DO I = NKS(J),1,-1
          K = KORPG(I,J)
          IF (K.NE.0) THEN
            NP = NP + 1
            XP(NP) = RVERTP(4,K)
            YP(NP) = ZVERTP(4,K)
          ENDIF
        ENDDO
        DO J = IRTRAP+1,NRS
          NP = NP + 1
          I = 1
          K = KORPG(I,J)
          IF (K.EQ.0) K = KORPG(I+1,J)
          XP(NP) = RVERTP(1,K)
          YP(NP) = ZVERTP(1,K)
        ENDDO
        DO J = IRSEP,IRWALL-1
          NP = NP + 1
          I = 1
          K = KORPG(I,J)
          IF (K.EQ.0) K = KORPG(I+1,J)
          XP(NP) = RVERTP(1,K)
          YP(NP) = ZVERTP(1,K)
        ENDDO
        NP = NP + 1
        XP(NP) = RVERTP(2,K)
        YP(NP) = ZVERTP(2,K)
C
C  ESCAPE REGION (NOTE THAT LAST CELL REPEATS FIRST)
C
        NE = 0
        J = 2
        DO I = 1,NKS(J)
          NE = NE + 1
          K = KORPG(I,J)
          XE(NE) = RVERTP(4,K)
          YE(NE) = ZVERTP(4,K)
        ENDDO
C
C  INITIALISE REGIONS
C
        CALL GA15A(NP,KIND,WORKP,IW,INDP,LND,XP,YP,T,XD,YD,LP)
        IF (INDP(1,1).GT.0) THEN
          WRITE(6,*) ' GA15A ERROR ',INDP(1,1),' FOR PLASMA BOUNDARY'
          STOP ' '
        ENDIF
        CALL GA15A(NE,KIND,WORKE,IW,INDE,LND,XE,YE,T,XD,YD,LP)
        IF (INDE(1,1).GT.0) THEN
          WRITE(6,*) ' GA15A ERROR ',INDE(1,1),' FOR ESCAPE REGION'
          STOP ' '
        ENDIF
C
        SETUP = .FALSE.
      ENDIF
C
      IERR = 0
C
C  IF THERE IS A CURRENT CELL, CHECK IT AND ITS NEIGHBOURS
C
      IF (IK.GT.0) THEN
C
C  START WITH LAST KNOWN CELL
C
        IF (INCELLA(R,Z,IK,IR)) RETURN
C
C  OTHERWISE, CHECK NEAREST NEIGHBOURS
C
        IF (INCELLA(R,Z,IKINS(IK,IR),IRINS(IK,IR))) THEN
          I  = IKINS(IK,IR)
          IR = IRINS(IK,IR)
          IK = I
          RETURN
        ENDIF
        IF (INCELLA(R,Z,IKOUTS(IK,IR),IROUTS(IK,IR))) THEN
          I  = IKOUTS(IK,IR)
          IR = IROUTS(IK,IR)
          IK = I
          RETURN
        ENDIF
        IF (IK.GT.1) THEN
          IF (INCELLA(R,Z,IK-1,IR)) THEN
            IK = IK - 1
            RETURN
          ENDIF
        ENDIF
        IF (IK.LT.NKS(IR)) THEN
          IF (INCELLA(R,Z,IK+1,IR)) THEN
            IK = IK + 1
            RETURN
          ENDIF
        ENDIF
      ENDIF
C
C  CHECK IF POINT IS IN PLASMA
C
      CALL GA15B(R,Z,RESULT,NP,KIND,WORKP,IW,INDP,LND,XP,YP,T,XD,YD,LP)
      IF (RESULT.LE.0) THEN
        IK = -1
        IR = -1
        RETURN
      ENDIF
C
C  CHECK IF POINT IS IN ESCAPE REGION
C
      CALL GA15B(R,Z,RESULT,NE,KIND,WORKE,IW,INDE,LND,XE,YE,T,XD,YD,LP)
      IF (RESULT.GE.0) THEN
        IK = -2
        IR = -2
        RETURN
      ENDIF
C
C  FINALLY, LOOP OVER PLASMA CELLS
C
      DO 20 J = 1, NRS
        DO 10 I = 1, NKS(J)
          IF (INCELLA(R,Z,I,J)) THEN
            IK = I
            IR = J
            RETURN
          ENDIF
   10   CONTINUE
   20 CONTINUE
C
C  POINT STILL NOT FOUND
C
      IERR = 1
      IK = 0
      IR = 0
      RETURN
      END
c
c
c
      SUBROUTINE LOS3d(TVALS,robs,zobs,couts,wres,vs,numthe,avpts)
      use mod_params
      use mod_cgeom
      IMPLICIT NONE
C
C  *********************************************************************
C  *                                                                   *
C  *  LOS3D :  INTEGRATE THE VARIABLE VS ALONG numthe SIGHT LINES      *
C  *           FROM A COMMON OBSERVATION POINT (ROBS,ZOBS).  SPATIAL   *
C  *           RESOLUTION IS SIMULATED BY 2*(2*avpts - 1) CHORDS       *
C  *           WITHIN A CONE CENTERED ON EACH SIGHT LINE, AND ADDED    *
C  *           WITH A WEIGHTING wres.                                  *
C  *                                                                   *
C  *            CHRIS FARRELL  (HUNTERSKIL)  MARCH 1989                *
C  *            LORNE HORTON   (JET)         JULY  1993                *
C  *            JOHN O'ROURKE  (JET)         AUG   1993                *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE "PARAMS"
c     include 'params'
C     INCLUDE "CGEOM"
c     include 'cgeom'
      INTEGER NUMTHE,AVPTS
      integer i,j,k,m,mc,ifnd,ir,ik,iv,iv1
      REAL    TVALS(MAXTHE),TOUTS(MAXTHE),
     >        ROBS,ZOBS,VS(MAXNKS,MAXNRS)
      REAL COUTS(MAXCH3,MAXCH3,3),wres(maxch3,maxch3)
C
      real del,run,x,y,z,rmaj
      real ang,angt,angm,rv1,rv2,zv1,zv2
c
      logical newinj,outofgrid
C
c     Initialization
c
      newinj = .true.
      outofgrid = .false.
c
      m=2*(2*avpts - 1)
      mc=3*avpts - 1
      del=0.01
C
C  LOOP OVER SIGHT LINES
C
      DO 200 I = 1, NUMTHE
      write(6,*)' chord number',i,numthe,m,mc
        TVALS(I) = 0.0
C
C  LOOP OVER CHORDS FOR AVERAGING
C
        DO 100 J = 1, M
c
c  Because the sight cone is scanned twice (in perpendicular directions)
c  the central sight line is counted twice -- no need to recalculate it!
c
          if(j.eq.mc)then
          tvals(j)=tvals(avpts)
          goto 100
          endif
c
c  Begin moving along sight line...
c
          run=0.
300       run=run+1.
          x=robs+ (run-0.5)*del*couts(i,j,1)
          y=      (run-0.5)*del*couts(i,j,2)
          z=zobs+ (run-0.5)*del*couts(i,j,3)
          rmaj= (x*x + y*y)**0.5
c
c
c  ... until you are outside the divertor
c
c          if(z.lt.-1.8.or.z.gt.-1.2)goto 301
c          write (6,*) 'Info:',run,x,y,z
c
c          if(z.lt.1.8.and.z.gt.1.2) then
c          write (6,*) 'Info2:',rmaj
c          if(rmaj.lt.2.2.or.rmaj.gt.3.2)goto 301
C
C  LOOP OVER PLASMA CELLS
C
          if (z.gt.zmax.or.z.lt.zmin.or.
     >        rmaj.gt.rmax.or.rmaj.lt.rmin) goto 301

          outofgrid = .false.
          call gridpos(ik,ir,rmaj,z,newinj,outofgrid)
          if (.not.outofgrid) then
             TVALS(I) = TVALS(I) + VS(IK,IR)*del*wres(i,j)
          endif
c
          goto 300

c          ifnd=0
c          DO 20 IR = 1, NRS
c            DO 10 IK = 1, NKS(IR)
c
c  Once you have found that you are in a cell, there is no need to
c  check the rest.  This also avoids double counting due to round-off
c  errors.
c
c          if(ifnd.ne.0)goto 10
c
c  Calculate the angles subtended at (rmaj,z)
c  by adjoining pairs of vertices of a given polygon.
c  If (rmaj,z) is within the polygon, the sum of the three
c  smallest angles will be greater than 180 degrees.
c
c              K = KORPG(IK,IR)
c              angt=0.
c              angm=0.
c              do 30 iv=1,4
c              if(iv.le.3)then
c              iv1=iv+1
c              else
c              iv1=1
c              endif
c              rv1=rvertp(iv,k)
c              zv1=zvertp(iv,k)
c              rv2=rvertp(iv1,k)
c              zv2=zvertp(iv1,k)
c              call angl(rmaj,z,rv1,zv1,rv2,zv2,ang)
c              angt=angt+ang
c              if(ang.gt.angm)angm=ang
c   30         continue
c              angt=angt-angm
c              if(angt.ge.pi)then
c              TVALS(I) = TVALS(I) + VS(IK,IR)*del*wres(i,j)
c              ifnd=1
c
c             if(i.eq.1.and.j.eq.1)then
c
c              write(6,401)rmaj,z,k,vs(ik,ir)
c  401         format(' r=',f6.2,' z=',f6.2,' polygon #',i4,' e=',e12.3)
c              write (6,*) i,j,ik,ir,tvals(i),del,wres(i,j)
c
c             endif
c
c              endif
c   10       CONTINUE
c   20     CONTINUE
c

c
  301   continue
  100   CONTINUE
  200 CONTINUE
C
      RETURN
      END

C
C
C
      SUBROUTINE angl(r0,z0,r1,z1,r2,z2,ang)
      IMPLICIT NONE
c
c     Called from LOS3d - written by John O'Rourke at JET
c     It seems to return a direction cosine for 2 vectors
c
      REAL R0,z0,r1,z1,r2,z2,ang,v1n,v2n,dot
      real v1(2),v2(2)
c
      v1(1)=r1-r0
      v1(2)=z1-z0
      v2(1)=r2-r0
      v2(2)=z2-z0
      v1n=( v1(1)*v1(1) + v1(2)*v1(2) )**0.5
      v2n=( v2(1)*v2(1) + v2(2)*v2(2) )**0.5
c
      if(v1n.eq.0..or.v2n.eq.0.)then
         ang=0.
      else
         dot=( v1(1)*v2(1) + v1(2)*v2(2) ) / (v1n*v2n)
         if(dot.gt.1.)dot=1.
         if(dot.lt.-1.)dot=-1.
         ang=acos(dot)
      endif
c
      return
      end
C
C
C
      SUBROUTINE INTERS (NV,RV,ZV,XB,WB,NINT,DIST,SIDE)
      IMPLICIT NONE
C
C  *********************************************************************
C  *                                                                   *
C  *  INTERS:  FIND THE NUMBER OF INTERSECTIONS, IF ANY, BETWEEN A     *
C  *           LINE DEFINED BY THE POINT XB AND THE DIRECTION COSINES  *
C  *           WB WITH A POLYGON DEFINED BY ITS VERTICES               *
C  *           (RV(I),ZV(I)),I=1,NV.  THE ALGORITHM IS TAKEN FROM THE  *
C  *           NIMBUS ROUTINE G1.  POLYGONS ARE ASSUMED TO BE CONVEX   *
C  *           SO THAT THERE MUST BE EITHER 0 OR 2 INTERSECTIONS.  IN  *
C  *           ADDITION, IF INTERSECTIONS ARE FOUND, THE DISTANCE FROM *
C  *           XB TO EACH IS RETURNED, AS IS THE SIDE OF THE POLYGON   *
C  *           WHICH IS CROSSED IN EACH CASE.                          *
C  *                                                                   *
C  *            LORNE HORTON   (JET)         JULY  1993                *
C  *                                                                   *
C  *********************************************************************
C
      INTEGER NV,NINT,SIDE(2),LR(2)
      REAL    RV(NV),ZV(NV),XB(2),WB(2),DIST(2)
C
      INTEGER I,IP1
      REAL XI,YI,XIP1,YIP1,D
      REAL FNUM,DENO,T,XP,YP
      REAL FMIN,FMAX,TT(2),EPS
      DATA EPS/1.0E-5/
C
      NINT = 0
      TT(1) = 0.0
      DO 120 I = 1, NV
        IF (NINT.EQ.2) GOTO 120
        XI = RV(I)
        YI = ZV(I)
        IP1 = I + 1
        IF (IP1.GT.NV) IP1 = 1
        XIP1 = RV(IP1)
        YIP1 = ZV(IP1)
        D = YIP1 - YI
        IF (D) 90,80,90
   80   FNUM = YI - XB(2)
        DENO = WB(2)
        GOTO 100
   90   D = (XIP1-XI)/D
        FNUM = XI - XB(1) - D*(YI-XB(2))
        DENO = WB(1) - D*WB(2)
  100   IF (DENO) 110,120,110
  110   T = FNUM/DENO
        XP = XB(1) + WB(1)*T
        IF (XI.EQ.XIP1) XP = XI
        YP = XB(2) + WB(2)*T
        IF (YI.EQ.YIP1) YP = YI
        FMIN = AMIN1(XI,XIP1)
        FMAX = AMAX1(XI,XIP1)
        IF (XP.LT.FMIN .OR. XP.GT.FMAX) GO TO 120
        FMIN = AMIN1(YI,YIP1)
        FMAX = AMAX1(YI,YIP1)
        IF (YP.LT.FMIN .OR. YP.GT.FMAX) GO TO 120
        IF (NINT.EQ.1 .AND. ABS(T-TT(1)).LT.EPS) GO TO 120
        NINT = NINT + 1
        TT(NINT) = T
        LR(NINT) = I
  120 CONTINUE
      IF (NINT.EQ.0) GO TO 170
      DIST(1) = TT(1)
      DIST(2) = TT(2)
      SIDE(1) = LR(1)
      SIDE(2) = LR(2)
      IF (DIST(1).LT.DIST(2)) GO TO 170
      DIST(1) = TT(2)
      DIST(2) = TT(1)
      SIDE(1) = LR(2)
      SIDE(2) = LR(1)
  170 CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE CUT (XCUT,YCUT,NCUT,MXCUT,VS,R1,Z1,R2,Z2)
      use mod_params
      use mod_cgeom
      IMPLICIT NONE
C
C  *********************************************************************
C  *                                                                   *
C  *  CUT:     RETURN THE PROFILE FORMED WHEN THE TWO-DIMENSIONAL      *
C  *           DISTRIBUTION VS IS CUT BY A RAY STARTING AT THE POINT   *
C  *           (R1,Z1) AND PROCEEDING TO THE POINT (R2,Z2).  THE       *
C  *           X VECTOR IS THE DISTANCE FROM THE STARTING POINT WITH   *
C  *           THE FIRST POINT AT (R1,Z1) AND THE LAST AT (R2,Z2).     *
C  *           IF THE RAY DOES NOT INTERSECT THE PLASMA AT ALL, THEN   *
C  *           THESE ARE THE ONLY TWO POINTS RETURNED AND NCUT=2.      *
C  *           SINCE THE TWO END POINTS COULD BE IN THE SAME PLASMA    *
C  *           CELL, THE TEST FOR NO INTERSECTIONS SHOULD BE NCUT=2    *
C  *           AND YCUT(1)=0.0.  IF THE SUPPLIED VECTORS ARE NOT       *
C  *           BIG ENOUGH THE PROFILE IS TRUNCATED AND                 *
C  *           XCUT(NCUT=MXCUT) IS NOT AT (R2,Z2).                     *
C  *                                                                   *
C  *            LORNE HORTON   (JET)         JULY  1993                *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE "PARAMS"
c     include 'params'
C     INCLUDE "CGEOM"
c     include 'cgeom'
      INTEGER NCUT,MXCUT
      REAL    XCUT(MXCUT),YCUT(MXCUT),VS(MAXNKS,MAXNRS),R1,Z1,R2,Z2
C
      INTEGER IR,IK,K,NINT,SIDE(2),IRNEW,IKNEW,SMIN(2)
      REAL    XB(2),WB(2),THETA,DIST(2),DMIN,DMAX,EPS,ATAN2C
      DATA    EPS/1.0E-5/
      EXTERNAL ATAN2C
C
C  MXCUT MUST AT LEAST BE LONG ENOUGH TO HOLD FIRST AND LAST POINTS!
C
      IF (MXCUT.LT.2) THEN
        WRITE(6,*) ' ERROR IN CUT, MXCUT = ',MXCUT
        STOP
      ENDIF
C
      XB(1) = R1
      XB(2) = Z1
      DMAX = SQRT((R2-R1)*(R2-R1) + (Z2-Z1)*(Z2-Z1))
      IF (DMAX.LE. EPS) THEN
        WRITE(6,*) ' ERROR IN CUT, DMAX = ',DMAX
        STOP
      ENDIF
      THETA = ATAN2C(Z2-Z1,R2-R1)
      WB(1) = COS(THETA)
      WB(2) = SIN(THETA)
C
C  IS THE STARTING POINT IN THE PLASMA?
C
      DO 20 IR = 1, NRS
        DO 10 IK = 1, NKS(IR)
          K = KORPG(IK,IR)
          CALL INTERS(NVERTP(K),RVERTP(1,K),ZVERTP(1,K),
     >                XB,WB,NINT,DIST,SIDE)
          IF (NINT.EQ.2 .AND. DIST(1).LE.0.0 .AND. DIST(2).GE.0.0) THEN
C  YES
            NCUT = 1
            XCUT(NCUT) = 0.0
            YCUT(NCUT) = VS(IK,IR)
CCC         WRITE(6,*) IK,IR,XCUT(NCUT),YCUT(NCUT)
            IF (DIST(2).GT.EPS) THEN
              IF (NCUT+1.GT.MXCUT) GOTO 999
              IF (DIST(2).GT.DMAX) GOTO 900
              NCUT = NCUT + 1
              XCUT(NCUT) = DIST(2)-EPS
              YCUT(NCUT) = VS(IK,IR)
CCC           WRITE(6,*) IK,IR,XCUT(NCUT),YCUT(NCUT)
            ENDIF
            GOTO 30
          ENDIF
   10   CONTINUE
   20 CONTINUE
C  NO
      NCUT = 1
      XCUT(NCUT) = 0.0
      YCUT(NCUT) = 0.0
CCC   WRITE(6,*) -1,-1,XCUT(NCUT),YCUT(NCUT)
      GOTO 50
C
C  FIND NEIGHBOURING CELL, IF ANY
C
   30 IF (SIDE(2).EQ.1) THEN
C  BACKWARD
        IF (IK.GT.1 .AND. KORPG(IK-1,IR).NE.0) THEN
          IKNEW = IK - 1
          IRNEW = IR
        ELSE IF (IK.EQ.1 .AND. IR.LT.IRSEP) THEN
          IKNEW = NKS(IR) - 1
          IRNEW = IR
        ELSE
          GOTO 40
        ENDIF
      ELSE IF (SIDE(2).EQ.2) THEN
C  OUTSIDE
        IF (IROUTS(IK,IR).NE.IR .AND.
     >      KORPG(IKOUTS(IK,IR),IROUTS(IK,IR)).NE.0) THEN
          IKNEW = IKOUTS(IK,IR)
          IRNEW = IROUTS(IK,IR)
        ELSE
          GOTO 40
        ENDIF
      ELSE IF (SIDE(2).EQ.3) THEN
C  FORWARD
        IF (IK.LT.NKS(IR) .AND. KORPG(IK+1,IR).NE.0) THEN
          IKNEW = IK + 1
          IRNEW = IR
        ELSE IF (IK.EQ.NKS(IR) .AND. IR.LT.IRSEP) THEN
          IKNEW = 2
          IRNEW = IR
        ELSE
          GOTO 40
        ENDIF
      ELSE IF (SIDE(2).EQ.4) THEN
C  INSIDE
        IF (IRINS(IK,IR).NE.IR .AND.
     >      KORPG(IKINS(IK,IR),IRINS(IK,IR)).NE.0) THEN
          IKNEW = IKINS(IK,IR)
          IRNEW = IRINS(IK,IR)
        ELSE
          GOTO 40
        ENDIF
      ENDIF
C
C  GET PATH THROUGH NEW CELL
C
      IK = IKNEW
      IR = IRNEW
      K = KORPG(IK,IR)
      CALL INTERS(NVERTP(K),RVERTP(1,K),ZVERTP(1,K),
     >            XB,WB,NINT,DIST,SIDE)
      IF (NINT.NE.2) THEN
        WRITE(6,*) ' ERROR IN CUT, (IK,IR) = ',IK,IR,' NINT = ',NINT
        STOP
      ENDIF
C
C  IF PATH LENGTH IS GREATER THAN 2*EPS, ADD NEW POINTS TO PROFILE
C
      IF (DIST(2)-DIST(1).GT.2.0*EPS) THEN
        IF (NCUT+1.GT.MXCUT) GOTO 999
        IF (DIST(1)+EPS.GT.DMAX) GOTO 900
        NCUT = NCUT + 1
        XCUT(NCUT) = DIST(1) + EPS
        YCUT(NCUT) = VS(IK,IR)
CCC     WRITE(6,*) IK,IR,XCUT(NCUT),YCUT(NCUT)
        IF (NCUT+1.GT.MXCUT) GOTO 999
        IF (DIST(2).GT.DMAX) GOTO 900
        NCUT = NCUT + 1
        XCUT(NCUT) = DIST(2) - EPS
        YCUT(NCUT) = VS(IK,IR)
CCC     WRITE(6,*) IK,IR,XCUT(NCUT),YCUT(NCUT)
      ENDIF
C
C  LOOP BACK TO FIND NEXT CELL
C
      GOTO 30
C
C  ADD EXIT POINT FROM PLASMA
C
   40 IF (NCUT+1.GT.MXCUT) GOTO 999
      IF (DIST(2)+EPS.GT.DMAX) GOTO 800
      NCUT = NCUT + 1
      XCUT(NCUT) = DIST(2) + EPS
      YCUT(NCUT) = 0.0
CCC   WRITE(6,*) -IK,-IR,XCUT(NCUT),YCUT(NCUT)
C
C  FIND NEXT ENTRY TO PLASMA
C
   50 DMIN = HI
      IKNEW = 0
      IRNEW = 0
      DO IR = 1,NRS
        DO IK = 1,NKS(IR)
          K = KORPG(IK,IR)
          CALL INTERS(NVERTP(K),RVERTP(1,K),ZVERTP(1,K),
     >                XB,WB,NINT,DIST,SIDE)
          IF (NINT.EQ.2 .AND. DIST(2).GT.XCUT(NCUT)+EPS .AND.
     >        DIST(2).LT.DMIN) THEN
            DMIN = DIST(2)
            IKNEW = IK
            IRNEW = IR
          ENDIF
        ENDDO
      ENDDO
C
      IF (DMIN.LT.HI) THEN
        IK = IKNEW
        IR = IRNEW
        K = KORPG(IK,IR)
        CALL INTERS(NVERTP(K),RVERTP(1,K),ZVERTP(1,K),
     >              XB,WB,NINT,DIST,SIDE)
        IF (NINT.NE.2) THEN
          WRITE(6,*) ' ERROR IN CUT, (IK,IR) = ',IK,IR,' NINT = ',NINT
          STOP
        ENDIF
C
C  ADD ENTRY POINT TO PLASMA
C
        IF (NCUT+1.GT.MXCUT) GOTO 999
        IF (DIST(1).GT.DMAX) GOTO 800
        NCUT = NCUT + 1
        XCUT(NCUT) = DIST(1) - EPS
        YCUT(NCUT) = 0.0
CCC     WRITE(6,*) -IK,-IR,XCUT(NCUT),YCUT(NCUT)
C
C  IF PATH LENGTH IS GREATER THAN 2*EPS, ADD NEW POINTS TO PROFILE
C
        IF (DIST(2)-DIST(1).GT.2.0*EPS) THEN
          IF (NCUT+1.GT.MXCUT) GOTO 999
          IF (DIST(1)+EPS.GT.DMAX) GOTO 900
          NCUT = NCUT + 1
          XCUT(NCUT) = DIST(1) + EPS
          YCUT(NCUT) = VS(IK,IR)
CCC       WRITE(6,*) IK,IR,XCUT(NCUT),YCUT(NCUT)
          IF (NCUT+1.GT.MXCUT) GOTO 999
          IF (DIST(2).GT.DMAX) GOTO 900
          NCUT = NCUT + 1
          XCUT(NCUT) = DIST(2) - EPS
          YCUT(NCUT) = VS(IK,IR)
CCC       WRITE(6,*) IK,IR,XCUT(NCUT),YCUT(NCUT)
        ENDIF
C
C  LOOP BACK TO FIND NEXT CELL
C
        GOTO 30
       endif
c
c      There is no obvious reason for the following block to
c      be attached to the above IF-Block - since flow of control branche
c      elsewhere if the dmin.lt.hi branch is taken. Furthermore - some
c      code branches into label 800 ... which is illegal inside an
c      IF-block.
c
c      ELSE
C
C  NO MORE PLASMA; ADD END POINT AND RETURN
C
  800 IF (NCUT+1.GT.MXCUT) GOTO 999
      NCUT = NCUT + 1
      XCUT(NCUT) = DMAX
      YCUT(NCUT) = 0.0
CCC   WRITE(6,*) -1,-1,XCUT(NCUT),YCUT(NCUT)
      GOTO 999
c
c      ENDIF
C
  900 NCUT = NCUT + 1
      XCUT(NCUT) = DMAX
      YCUT(NCUT) = VS(IK,IR)
CCC   WRITE(6,*) IK,IR,XCUT(NCUT),YCUT(NCUT)
C
  999 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE REFLECT
      use mod_params
      use mod_cgeom
      use mod_grbound
      use mod_comtor
      use mod_pindata
      use mod_dynam4
      use mod_outxy
      IMPLICIT NONE
C
C  *********************************************************************
C  *                                                                   *
C  *  REFLECT: THE VARIOUS JET CODES WILL ONLY WORK WITH UPPER         *
C  *           X-POINTS SO LOWER X-POINT GRIDS ARE REFLECTED IN THE    *
C  *           R=0 PLANE BY GRID2D.  THIS OPERATION IS NOW FLAGGED     *
C  *           IN THE EQUILIBRIUM FILE WITH REFCT=1.  (OLD EQUILIBRIUM *
C  *           FILES WITHOUT THE FLAG SHOULD RETURN REFCT=0 AND WILL   *
C  *           NOT BE REFLECTED HERE - AND THUS MAY GIVE STRANGE       *
C  *           DIAGNOSTIC SIGNAL SIMULATIONS!)  THE VARIABLES WHICH    *
C  *           ARE AFFECTED ARE:                                       *
C  *                  - Z0, ZXP, ZMIN, ZMAX, ZS(IK,IR), ZVERTP,        *
C  *                    IRXYS, IKXYS, IFXYS                            *
c  *                  - ZCW, ZW (and recalulated GA15 wall specs)      *
c  *                  - ZVESM (NIMBUS vessel and pump wall)            *
C  *                  - WALKS(*,2), HWALKS(*,2)                        *
c  *                  - WALLPT array
C  *                                                                   *
C  *            LORNE HORTON   (JET)         JULY  1993                *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE "PARAMS"
c     include 'params'
C     INCLUDE "CGEOM"
c     include 'cgeom'
c
c     include 'grbound'
c     include 'comtor'
c     include 'pindata'
c     include 'dynam4'
c
c     include 'outxy'
C
      INTEGER IR,IK,K,I,IX,IY,IRTMP(MAXGYS),IKTMP(MAXGYS),IFTMP(MAXGYS)
      integer kind
      REAL    ZTMP
c
      real atan2c
      external atan2c 
C
      Z0   = -Z0
      ZXP  = -ZXP
      ZTMP = -ZMIN
      ZMIN = -ZMAX
      ZMAX =  ZTMP
      DO IR = 1,NRS
        DO IK = 1,NKS(IR)
          ZS(IK,IR) = -ZS(IK,IR)
        ENDDO
      ENDDO
      DO K = 1,NPOLYP
        DO I = 1,NVERTP(K)+1
          ZVERTP(I,K) = -ZVERTP(I,K)
        ENDDO
      ENDDO
      DO IX = 1,NXS
        DO IY = 1,NYS
          IRTMP(IY) =  IRXYS(IX,IY)
          IKTMP(IY) =  IKXYS(IX,IY)
          IFTMP(IY) =  IFXYS(IX,IY)
        ENDDO
        DO IY = 1,NYS
          IRXYS(IX,IY) = IRTMP(NYS-IY+1)
          IKXYS(IX,IY) = IKTMP(NYS-IY+1)
          IFXYS(IX,IY) = IFTMP(NYS-IY+1)
        ENDDO
      ENDDO
c
c     Reflect and recalulate both the outer wall and core
c     wall specifications.
c
c     Ion wall - OUTER boundary
c
      do i = 1,ionwpts
         ziw(i) = -ziw(i)
      enddo
c
      KIND = 1
      CALL GA15A(IONWPTS,KIND,iwWORK,4*MAXPTS,iwINDW,MAXPTS,
     >             RIW,ZIW,iwTDUM,iwXDUM,iwYDUM,6)
c
c     Ion wall - CORE boundary
c
      do i = 1,ioncpts
         zcw(i) = -zcw(i)
      end do

      CALL GA15A(IONCPTS,KIND,icWORK,4*MAXPTS,icINDW,MAXPTS,
     >             RCW,ZCW,icTDUM,icXDUM,icYDUM,6)
c
c     Neutral Wall - just for completeness
c     NOTE: ALL of the launch data in the wallpt array
c           dealing with wall orientations have NOT been
c           reflected - any code that uses this will need
c           to add reflection later OR be written to handle
c           the inverted orientation.
c
      do i = 1,pcnt
         zw(i) = -zw(i)
      end do
c
c     Vessel Wall
c
      do i = 1, nves
         zves(i) = - zves(i)
      end do
C
      CALL GA15A(NVES,KIND,nwWORK,4*MAXPTS,nwINDW,MAXPTS,
     >             RVES,ZVES,nwTDUM,nwXDUM,nwYDUM,6)
c
c     NIMBUS vessel and pump walls
c
      do i = 1, nvesm+nvesp
         do k = 1, 2
            zvesm(i,k) = - zvesm(i,k)
         end do
      end do
c
c     Adjust Z values in wallpt array and recalculate angles
c
      do i = 1,wallpts
c
c        Invert Z coordinates for segment
c
         wallpt(i,2)  = -wallpt(i,2)
         wallpt(i,21) = -wallpt(i,21)
         wallpt(i,23) = -wallpt(i,23)
c
c        Recalculate elements 8 and 9 - angles forward and backward
c
         WALLPT(I,8) = ATAN2C(wallpt(i,21)-wallpt(i,2),
     >                          wallpt(i,20)-wallpt(i,1))
         WALLPT(I,9) = ATAN2C(wallpt(i,23)-wallpt(i,2),
     >                          wallpt(i,22)-wallpt(i,1))
c
      end do
c
c     trajectories of ions and neutrals
c
      do i = 1, maxnws
        walks(i,2) = -walks(i,2)
        hwalks(i,2) = -hwalks(i,2)
      enddo
C
      RETURN
      END
c
c
c
      subroutine gridcoords (ix,iy,ik,ir,in)
      use mod_params
      use mod_outxy
      implicit none
      integer ix,iy,ik,ir,in
c     include 'params'
c     include 'outxy'
c
c     Return the ik,ir cooridnates for the ix,iy bin.
c     The separate routine is necessary because div and
c     out now differ in their definition of the following
c     arrays.
c
      ik = ikxys(ix,iy)
      ir = irxys(ix,iy)
      in = ifxys(ix,iy)
      return
      end
c
c
c
      subroutine setkval(kval,alphae,ik,ir)
      use mod_params
      use mod_cgeom
      implicit none
      real kval,alphae
      integer ik,ir
c
c     Required geometry data
c
c     include 'params'
c     include 'cgeom'
c
c     This subroutine assigns a known distribution to the
c     KVALS array in order to facilitate testing of the
c     plotting routines.
c
c
c     The first and only option so far is:
c
c     INT (S1 to S2)  EXP(-alphae * S) which can
c     be integrated to give an analytic expression.
c
c     Local variables
c
      real s1, s2
c
c
c
      if (kss(ik,ir).lt.(0.5*ksmaxs(ir))) then
         if (ik.eq.1) then
            s1 = 0.0
            s2 = kss(ik,ir) + kfords(ik,ir)/2.0
         else
            s1 = kss(ik,ir) - kbacds(ik,ir)/2.0
            s2 = kss(ik,ir) + kfords(ik,ir)/2.0
         endif
         kval = (exp(-alphae*s1)-exp(-alphae*s2))/alphae
      elseif (kss(ik,ir).ge.(0.5*ksmaxs(ir))) then
         if (ik.eq.nks(ir)) then
            s1 = 0.0
            s2 = ksmaxs(ir) - kss(ik,ir) + kbacds(ik,ir)/2.0
         else
            s1 = ksmaxs(ir) - kss(ik,ir) - kfords(ik,ir)/2.0
            s2 = ksmaxs(ir) - kss(ik,ir) + kbacds(ik,ir)/2.0
         endif
         kval = (exp(-alphae*s1)-exp(-alphae*s2))/alphae
      endif
      return
      end


      SUBROUTINE LOS3DA(TVALS,ROBS,POBS,ZOBS,COUTS,WRES,VS,NUMTHE,AVPTS)
      use mod_params
      use mod_cgeom
      IMPLICIT NONE
C
C  *********************************************************************
C  *                                                                   *
C  *  LOS3D :  INTEGRATE THE VARIABLE VS ALONG NUMTHE SIGHT LINES      *
C  *           FROM A COMMON OBSERVATION POINT (ROBS,ZOBS).  SPATIAL   *
C  *           RESOLUTION IS SIMULATED BY 2*(2*AVPTS - 1) CHORDS       *
C  *           WITHIN A CONE CENTERED ON EACH SIGHT LINE, AND ADDED    *
C  *           WITH A WEIGHTING WRES.                                  *
C  *                                                                   *
C  *            CHRIS FARRELL  (HUNTERSKIL)  MARCH 1989                *
C  *            LORNE HORTON   (JET)         JULY  1993                *
C  *            JOHN O'ROURKE  (JET)         AUG   1993                *
C  *                                                                   *
C  *********************************************************************
C
c     include 'params'
C     INCLUDE "CGEOM"
c     include 'cgeom'
c      include 'PPPARA'
c      include 'PPUNIT'
c      include 'PPGEOM'
C
      INTEGER NUMTHE,AVPTS
      INTEGER I,J,K,M,MCP,IFND,IR,IK,IV,IV1,IERR
      integer lerr
      data    lerr/6/
      REAL    TVALS(MAXTHE),
     >        ROBS,POBS,ZOBS,VS(MAXNKS,MAXNRS)
      REAL    COUTS(MAXTHE,MAXCH3,3),WRES(MAXCH3)
C
      REAL DEL,RUN,X,Y,Z,R,RUNMAX
      REAL ANG,ANGT,ANGM,RV1,RV2,ZV1,ZV2
      LOGICAL LIN
C
      M=2*(2*AVPTS - 1)
      MCP=3*AVPTS - 1
      DEL=0.001
      RUNMAX=10000.
C
C  LOOP OVER SIGHT LINES
C
      DO 200 I = 1, NUMTHE
        TVALS(I) = 0.0
C
C  LOOP OVER CHORDS FOR AVERAGING
C
        DO 100 J = 1, M
C
C  BECAUSE THE SIGHT CONE IS SCANNED TWICE (IN PERPENDICULAR DIRECTIONS)
C  THE CENTRAL SIGHT LINE IS COUNTED TWICE -- NO NEED TO RECALCULATE IT
C
          IF (J.EQ.MCP) GOTO 100
C
C  BEGIN MOVING ALONG SIGHT LINE...
C
          RUN=0.
          IK=0
          IR=0
          LIN = .FALSE.
300       RUN=RUN+1.
          X=ROBS*COS(POBS) + (RUN-0.5)*DEL*COUTS(I,J,1)
          Y=ROBS*SIN(POBS) + (RUN-0.5)*DEL*COUTS(I,J,2)
          Z=ZOBS           + (RUN-0.5)*DEL*COUTS(I,J,3)
C
C  ... UNTIL YOU HAVE REACHED A MAXIMUM NUMBER OF STEPS...
C
          IF (RUN.GT.RUNMAX) GOTO 100
C
C  ... OR UNTIL YOU HAVE ENTERED THE VESSEL AND THEN EXITED AGAIN
C
          R = (X*X + Y*Y)**0.5
          IF(Z.LT.-2.1.OR.Z.GT.+2.1.OR.R.LT.+1.7.OR.R.GT.+4.3) THEN
            IF (LIN) THEN
              GOTO 100
            ELSE
              GOTO 300
            ENDIF
          ELSE
            LIN = .TRUE.
          ENDIF
C
C  FIND PLASMA CELL
C
          CALL GETCELLA(R,Z,IK,IR,IERR)
          IF (IERR.NE.0) THEN
            WRITE(LERR,*) ' GETCELL ERROR - POINT NOT LOCATED: (R,Z)= ',
     >                 R,Z
            GOTO 300
          ENDIF
C
C  UPDATE INTEGRAL, COUNTING THE CENTRAL LINE TWICE
C
          TVALS(I) = TVALS(I) + VS(IK,IR)*DEL*WRES(J)
          IF (J.EQ.AVPTS) TVALS(I) = TVALS(I) + VS(IK,IR)*DEL*WRES(J)
C
          GOTO 300
  100   CONTINUE
  200 CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE GETCELLA(R,Z,IK,IR,IERR)
      use mod_params
      use mod_cgeom
      IMPLICIT NONE
C
C  *********************************************************************
C  *                                                                   *
C  *  GETCELL: FIND THE CELL WHICH CONTAINS THE GIVEN POINT (R,Z).     *
C  *           THE INTEGERS IK AND IR ARE BOTH INPUT AND OUTPUT,       *
C  *           AND SHOULD CONTAIN THE LAST KNOWN CELL ALONG A TRACK    *
C  *           IN ORDER TO SPEED THE SEARCH PROCEDURE. IF THE INPUT    *
C  *           POINT IS OUTSIDE THE PLASMA IK AND IR ARE RETURNED      *
C  *           AS -1, IF THE POINT IS IN THE CENTRAL ESCAPE REGION     *
C  *           IK = IR = -2 IS RETURNED.  IF THE ROUTINE FAILS TO      *
C  *           LOCATE THE POINT IERR IS SET TO 1 AND IK = IR = 0.      *
C  *                                                                   *
C  *            LORNE HORTON   (JET)         SEPTEMBER 1993            *
C  *                                                                   *
C  *********************************************************************
C
c     include 'params'
C     INCLUDE "CGEOM"
c     include 'cgeom'
c      include 'PPPARA'
c      include 'PPUNIT'
c      include 'PPGEOM'
C
      INTEGER IK,IR,IERR
      REAL    R,Z
C
      INTEGER I,J,K
      CHARACTER LEQUIL*72
      integer lshot
      LOGICAL INCELLA,SETUP
      DATA    SETUP/.TRUE./,LEQUIL/' '/,lshot/0/
C
      INTEGER   IW,LND
      PARAMETER (LND = 128, IW = 4*LND)
      INTEGER   NP,NE,KIND,INDP(2,LND),INDE(2,LND),LP
      REAL      WORKP(IW),WORKE(IW),
     >          XP(2*(MAXNKS+MAXNRS)+1),YP(2*(MAXNKS+MAXNRS)+1),
     >          XE(MAXNKS),YE(MAXNKS),
     >          T,XD,YD,RESULT
      DATA      KIND/1/,LP/-1/
C
C  WHEN TREATING MULTIPLE FILES FROM THE D.A.M.U. POST PROCESSOR
C  THE EQUILIBRIUM MAY CHANGE.  CHECK FOR THIS AND REINITIALISE, IF
C  NECESSARY.
C
c      IF (lequil.NE.EQUIL) THEN
c        SETUP = .TRUE.
c        LEQUIL = EQUIL
c      ENDIF
      IF (lshot.NE.ishot) THEN
        SETUP = .TRUE.
        Lshot = ishot
      ENDIF
C
C  FIRST TIME IN, SETUP POLYGONS DEFINING PLASMA BOUNDARY AND
C  ESCAPE REGION
C
      IF (SETUP) THEN
C
C  PLASMA BOUNDARY
C
        NP = 0
        J = IRWALL - 1
        DO I = 1,NKS(J)
          K = KORPG(I,J)
C  CHECK TO SEE IF VIRTUAL POINTS EXIST AT ENDS OF RINGS!
          IF (K.NE.0) THEN
            NP = NP + 1
            XP(NP) = RVERTP(2,K)
            YP(NP) = ZVERTP(2,K)
          ENDIF
        ENDDO
        DO J = IRWALL-1,IRSEP,-1
          NP = NP + 1
          I = NKS(J)
          K = KORPG(I,J)
          IF (K.EQ.0) K = KORPG(I-1,J)
          XP(NP) = RVERTP(3,K)
          YP(NP) = ZVERTP(3,K)
        ENDDO
        DO J = NRS,IRTRAP+1,-1
          NP = NP + 1
          I = NKS(J)
          K = KORPG(I,J)
          IF (K.EQ.0) K = KORPG(I-1,J)
          XP(NP) = RVERTP(3,K)
          YP(NP) = ZVERTP(3,K)
        ENDDO
        J = IRTRAP + 1
        DO I = NKS(J),1,-1
          K = KORPG(I,J)
          IF (K.NE.0) THEN
            NP = NP + 1
            XP(NP) = RVERTP(4,K)
            YP(NP) = ZVERTP(4,K)
          ENDIF
        ENDDO
        DO J = IRTRAP+1,NRS
          NP = NP + 1
          I = 1
          K = KORPG(I,J)
          IF (K.EQ.0) K = KORPG(I+1,J)
          XP(NP) = RVERTP(1,K)
          YP(NP) = ZVERTP(1,K)
        ENDDO
        DO J = IRSEP,IRWALL-1
          NP = NP + 1
          I = 1
          K = KORPG(I,J)
          IF (K.EQ.0) K = KORPG(I+1,J)
          XP(NP) = RVERTP(1,K)
          YP(NP) = ZVERTP(1,K)
        ENDDO
        NP = NP + 1
        XP(NP) = RVERTP(2,K)
        YP(NP) = ZVERTP(2,K)
C
C  ESCAPE REGION (NOTE THAT LAST CELL REPEATS FIRST)
C
        NE = 0
        J = 2
        DO I = 1,NKS(J)
          NE = NE + 1
          K = KORPG(I,J)
          XE(NE) = RVERTP(4,K)
          YE(NE) = ZVERTP(4,K)
        ENDDO
C
C  INITIALISE REGIONS
C
        CALL GA15A(NP,KIND,WORKP,IW,INDP,LND,XP,YP,T,XD,YD,LP)
        IF (INDP(1,1).GT.0) THEN
          WRITE(6,*) ' GA15A ERROR ',INDP(1,1),' FOR PLASMA BOUNDARY'
          STOP ' '
        ENDIF
        CALL GA15A(NE,KIND,WORKE,IW,INDE,LND,XE,YE,T,XD,YD,LP)
        IF (INDE(1,1).GT.0) THEN
          WRITE(6,*) ' GA15A ERROR ',INDE(1,1),' FOR ESCAPE REGION'
          STOP ' '
        ENDIF
C
        SETUP = .FALSE.
      ENDIF
C
      IERR = 0
C
C  IF THERE IS A CURRENT CELL, CHECK IT AND ITS NEIGHBOURS
C
      IF (IK.GT.0) THEN
C
C  START WITH LAST KNOWN CELL
C
        IF (INCELLA(R,Z,IK,IR)) RETURN
C
C  OTHERWISE, CHECK NEAREST NEIGHBOURS
C
        IF (INCELLA(R,Z,IKINS(IK,IR),IRINS(IK,IR))) THEN
          I  = IKINS(IK,IR)
          IR = IRINS(IK,IR)
          IK = I
          RETURN
        ENDIF
        IF (INCELLA(R,Z,IKOUTS(IK,IR),IROUTS(IK,IR))) THEN
          I  = IKOUTS(IK,IR)
          IR = IROUTS(IK,IR)
          IK = I
          RETURN
        ENDIF
        IF (IK.GT.1) THEN
          IF (INCELLA(R,Z,IK-1,IR)) THEN
            IK = IK - 1
            RETURN
          ENDIF
        ENDIF
        IF (IK.LT.NKS(IR)) THEN
          IF (INCELLA(R,Z,IK+1,IR)) THEN
            IK = IK + 1
            RETURN
          ENDIF
        ENDIF
      ENDIF
C
C  CHECK IF POINT IS IN PLASMA
C
      CALL GA15B(R,Z,RESULT,NP,KIND,WORKP,IW,INDP,LND,XP,YP,T,XD,YD,LP)
      IF (RESULT.LE.0) THEN
        IK = -1
        IR = -1
        RETURN
      ENDIF
C
C  CHECK IF POINT IS IN ESCAPE REGION
C
      CALL GA15B(R,Z,RESULT,NE,KIND,WORKE,IW,INDE,LND,XE,YE,T,XD,YD,LP)
      IF (RESULT.GE.0) THEN
        IK = -2
        IR = -2
        RETURN
      ENDIF
C
C  FINALLY, LOOP OVER PLASMA CELLS
C
      DO 20 J = 1, NRS
        DO 10 I = 1, NKS(J)
          IF (INCELLA(R,Z,I,J)) THEN
            IK = I
            IR = J
            RETURN
          ENDIF
   10   CONTINUE
   20 CONTINUE
C
C  POINT STILL NOT FOUND
C
      IERR = 1
      IK = 0
      IR = 0
      RETURN
      END
C
C
C
      LOGICAL FUNCTION INCELLA(R,Z,IK,IR)
      use mod_params
      use mod_cgeom
      IMPLICIT NONE
      INTEGER IK,IR
      REAL R,Z
c     include 'params'
C     INCLUDE "CGEOM"
c     include 'cgeom'
c      include 'PPPARA'
c      include 'PPUNIT'
c      include 'PPGEOM'
      INTEGER K,IV,IV1
      REAL*8 R8,Z8,ANG,ANGT,ANGM
      REAL*8 RV1,ZV1,RV2,ZV2
      REAL*8 PI8
      DATA PI8/3.141592654D0/
C
      INCELLA = .FALSE.
C
C  CALCULATE THE ANGLES SUBTENDED AT (R,Z)
C  BY ADJOINING PAIRS OF VERTICES OF A GIVEN POLYGON.
C  IF (R,Z) IS WITHIN THE POLYGON, THE SUM OF THE THREE
C  SMALLEST ANGLES WILL BE GREATER THAN 180 DEGREES.
C
      K = KORPG(IK,IR)
C
C  VIRTUAL POINTS ARE OUTSIDE THE PLASMA AND HAVE NO POLYGON
C  ASSOCIATED WITH THEM
C
      IF (K.EQ.0) RETURN
C
      R8 = R
      Z8 = Z
      ANGT=0.
      ANGM=0.
      DO 10 IV=1,4
        IF(IV.LE.3)THEN
          IV1=IV+1
        ELSE
          IV1=1
        ENDIF
        RV1=RVERTP(IV,K)
        ZV1=ZVERTP(IV,K)
        RV2=RVERTP(IV1,K)
        ZV2=ZVERTP(IV1,K)
        CALL ANGLA(R8,Z8,RV1,ZV1,RV2,ZV2,ANG)
        ANGT=ANGT+ANG
        IF(ANG.GT.ANGM)ANGM=ANG
   10 CONTINUE
      ANGT=ANGT-ANGM
      IF (ANGT.GE.PI8) INCELLA = .TRUE.
C
      RETURN
      END
C
C
C
      SUBROUTINE ANGLA(R0,Z0,R1,Z1,R2,Z2,ANG)
      IMPLICIT NONE
      REAL*8 R0,Z0,R1,Z1,R2,Z2,ANG,V1N,V2N,DOT
      REAL*8 V1(2),V2(2)
      V1(1)=R1-R0
      V1(2)=Z1-Z0
      V2(1)=R2-R0
      V2(2)=Z2-Z0
      V1N=( V1(1)*V1(1) + V1(2)*V1(2) )**0.5
      V2N=( V2(1)*V2(1) + V2(2)*V2(2) )**0.5
      IF(V1N.EQ.0..OR.V2N.EQ.0.)THEN
        ANG=0.
      ELSE
        DOT=( V1(1)*V2(1) + V1(2)*V2(2) ) / (V1N*V2N)
        IF(DOT.GT.1.)DOT=1.
        IF(DOT.LT.-1.)DOT=-1.
        ANG=DACOS(DOT)
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE INTERSA(NV,RV,ZV,XB,WB,NINT,DIST,SIDE)
      IMPLICIT NONE
C
C  *********************************************************************
C  *                                                                   *
C  *  INTERS:  FIND THE NUMBER OF INTERSECTIONS, IF ANY, BETWEEN A     *
C  *           LINE DEFINED BY THE POINT XB AND THE DIRECTION COSINES  *
C  *           WB WITH A POLYGON DEFINED BY ITS VERTICES               *
C  *           (RV(I),ZV(I)),I=1,NV.  THE ALGORITHM IS TAKEN FROM THE  *
C  *           NIMBUS ROUTINE G1.  POLYGONS ARE ASSUMED TO BE CONVEX   *
C  *           SO THAT THERE MUST BE EITHER 0 OR 2 INTERSECTIONS.  IN  *
C  *           ADDITION, IF INTERSECTIONS ARE FOUND, THE DISTANCE FROM *
C  *           XB TO EACH IS RETURNED, AS IS THE SIDE OF THE POLYGON   *
C  *           WHICH IS CROSSED IN EACH CASE.                          *
C  *                                                                   *
C  *            LORNE HORTON   (JET)         JULY  1993                *
C  *                                                                   *
C  *********************************************************************
C
      INTEGER NV,NINT,SIDE(2)
      REAL    RV(NV),ZV(NV),XB(2),WB(2),DIST(2)
C
      INTEGER I,IP1,LR(2)
      REAL XI,YI,XIP1,YIP1,D
      REAL FNUM,DENO,T,XP,YP
      REAL FMIN,FMAX,TT(2),EPS
      DATA EPS/1.0E-5/
C
      NINT = 0
      TT(1) = 0.0
      DO 120 I = 1, NV
        IF (NINT.EQ.2) GOTO 120
        XI = RV(I)
        YI = ZV(I)
        IP1 = I + 1
        IF (IP1.GT.NV) IP1 = 1
        XIP1 = RV(IP1)
        YIP1 = ZV(IP1)
        D = YIP1 - YI
        IF (D) 90,80,90
   80   FNUM = YI - XB(2)
        DENO = WB(2)
        GOTO 100
   90   D = (XIP1-XI)/D
        FNUM = XI - XB(1) - D*(YI-XB(2))
        DENO = WB(1) - D*WB(2)
  100   IF (DENO) 110,120,110
  110   T = FNUM/DENO
        XP = XB(1) + WB(1)*T
        IF (XI.EQ.XIP1) XP = XI
        YP = XB(2) + WB(2)*T
        IF (YI.EQ.YIP1) YP = YI
        FMIN = AMIN1(XI,XIP1)
        FMAX = AMAX1(XI,XIP1)
        IF (XP.LT.FMIN .OR. XP.GT.FMAX) GO TO 120
        FMIN = AMIN1(YI,YIP1)
        FMAX = AMAX1(YI,YIP1)
        IF (YP.LT.FMIN .OR. YP.GT.FMAX) GO TO 120
        IF (NINT.EQ.1 .AND. ABS(T-TT(1)).LT.EPS) GO TO 120
        NINT = NINT + 1
        TT(NINT) = T
        LR(NINT) = I
  120 CONTINUE
      IF (NINT.EQ.0) GO TO 170
      DIST(1) = TT(1)
      DIST(2) = TT(2)
      SIDE(1) = LR(1)
      SIDE(2) = LR(2)
      IF (DIST(1).LT.DIST(2)) GO TO 170
      DIST(1) = TT(2)
      DIST(2) = TT(1)
      SIDE(1) = LR(2)
      SIDE(2) = LR(1)
  170 CONTINUE
C
      RETURN
      END
C
C
      subroutine adjustout(touts,numthe,zadj,robs,zobs)
      use mod_params
      implicit none
c     include 'params'
      integer numthe
      real touts(numthe),zadj,robs,zobs
c
c     This subroutine maps the THETA values in the Touts array from a
c     specific observation position onto a projected R co-ordinate plane
c     at value of Z specified by Z-adjustment.
c
c     David Elder     1995, June 8.
c
      integer i,j
      real    angadj
c
      if (zobs.gt.zadj) then
         angadj = 270.0
      elseif (zobs.lt.zadj) then
         angadj = 90.0
      else
         write (6,*) 'ERROR in mapping THTEA to R plotting coordinates:'
         write (6,*) 'Observation Z =',zobs,
     >               ' is the same as adjustment plane = ',zadj
         return
      endif
c
      do i = 1,numthe
c
c        Adjust Touts
c
         touts(i) = robs + (zobs-zadj) * tan((touts(i)-angadj)*degrad)
c
      end do
c
      return
      end
c
c
c
      subroutine adjustoutz(touts,numthe,radj,robs,zobs)
      use mod_params
      implicit none
c     include 'params'
      integer numthe
      real touts(numthe),radj,robs,zobs
c
c     This subroutine maps the THETA values in the Touts array from a
c     specific observation position onto a projected Z co-ordinate plane
c     at value of R specified by R-adjustment.
c
c     David Elder     1995, June 8.
c
      integer i,j
      real    angadj
c
      if (robs.gt.radj) then
         angadj = 270.0
      elseif (robs.lt.radj) then
         angadj = 90.0
      else
         write (6,*) 'ERROR in mapping THTEA to Z plotting coordinates:'
         write (6,*) 'Observation R =',robs,
     >               ' is the same as adjustment plane = ',radj
         return
      endif
c
      do i = 1,numthe
c
c        Adjust Touts
c
         touts(i) = zobs + (robs-radj) * tan((touts(i)-angadj)*degrad)
c
      end do
c
      return
      end
c
c
c
      logical function checkcell (ik,ir)
      use mod_params
      use mod_cgeom
      implicit none
      integer ik,ir
c     include 'params'
c     include 'cgeom'
c
c     This function checks to see if any of the cell corner points
c     are identical - if it finds degenerate corner points - it returns
c     TRUE.
c
      integer i,j,k,m,n
c
      checkcell = .false.
c
      k = korpg(ik,ir)
c
      n = nvertp(k)
c
c     Cells for which no polygon is defined are considered OK.
c
      if (n.eq.0) return
c
      do i = 1,n-1
         do j = i+1,n
c
            if (rvertp(i,k).eq.rvertp(j,k)) then
c
c                  Check to see if both coordinates of the polygon
c                  vertex are identical.
c
               if (zvertp(i,k).eq.zvertp(j,k)) then
                  checkcell = .true.
                  write (6,*) 'CHECKCELL: cell (',ik,ir,
     >                        ') is degenerate'
                  write(6,*) 'geom1:',nvertp(k)
                  write(6,*) 'geom2:',(rvertp(m,k),m=1,n)
                  write(6,*) 'geom3:',(zvertp(m,k),m=1,n)
                  return
               endif
            endif
c
         enddo
      enddo

      return
      end
C
C
      subroutine region(n,x,w,flag,pltmin,pltmax)
      implicit none
      integer n
      integer flag(n)
      real x(n),w(n),pltmin,pltmax
c
c     This subroutine marks a wall profile plot with the different
c     regions as passed back from NIMBUS.
c
c     Lorne Horton    1995, August 4.
c
      integer   i,j
      real      lold, lnew, labpos
      character lab(8)*2
      data lab/'OT','OC','OD','IT','IC','ID','MS','PV'/

c
c slmod begin
c...  True space:
      CALL PSPACE (0.1, 0.9, 0.1, 0.9)
c
c      call pspace(0.1, 0.9, 0.11, 0.89)
c slmod end
      call map   (x(1)-0.5*w(1),x(n)+0.5*w(n),pltmin,pltmax)
      call broken(5,5,5,5)
      call ctrmag(12)
      lold = x(1) - 0.5*w(1)
      j = 0
      do i = 2, n
        if (flag(i).ne.flag(i-1)) then
          lnew = x(i) - 0.5*w(i)
          call positn(lnew,pltmax)
          call join(lnew,pltmin)
          labpos = pltmin + (pltmax-pltmin)*(1.025+0.025*j)
          if (flag(i-1).gt.0.and.flag(i-1).le.8) then
             call pcscen(0.5*(lold+lnew),labpos,lab(flag(i-1)))
          endif
          lold = lnew
          j = 1-j
        endif
      enddo
      lnew = x(n) + 0.5*w(n)
      call positn(lnew,pltmax)
      call join(lnew,pltmin)
      labpos = pltmin + (pltmax-pltmin)*(1.025+0.025*j)
      if (flag(n).gt.0.and.flag(n).le.8) then
         call pcscen(0.5*(lold+lnew),labpos,lab(flag(n)))
      endif
      call full
c
      return
      end
c
c
c
      subroutine radproc(nizs,job,pradclev)
      use mod_params
      use mod_dynam2
      use mod_dynam3
      use mod_comtor
      use mod_cgeom
      implicit none
c     include    'params'
c     include    'dynam2'
c     include    'dynam3'
c     include    'comtor'
c     include    'cgeom'
c
      integer nizs
      real    pradclev(0:maxizs+1)
      character*(*) job
c
c     RADPROC: This routine processes the impurity ionization array
c              trying to determine specific quantities and
c              characteristics about where the radiation is
c              occuring. These quantities include the volume weighted
c              average of the local densities within the radiating
c              volume. The radiating volume is defined as the volume
c              that radiates the top 2/3's of the energy.
c
      integer ik,ir,in,count(0:maxizs+1),top(0:maxizs+1),num,iz,iz2
      integer bot
c
      real rad(maxnks,maxnrs,0:maxizs+1),radtot(0:maxizs+1)
      real neav(0:maxizs+1),nizav(0:maxizs+1),volav(0:maxizs+1)
      real teav(0:maxizs+1),lzav(0:maxizs+1)
      real pradt(0:maxizs+1),sdtmp
      real sdtot(maxnks,maxnrs)
c
      real radord(maxnks*maxnrs,0:maxizs+1,5)
c
      call rzero(radtot,maxizs+2)
      call rzero(rad,maxnks*maxnrs*(maxizs+2))
      call rzero(sdtot,maxnks*maxnrs)
c
c     Get array containg all Radiation and the TOTAL.
c
      do iz = 0,nizs
         do ir = 1,nrs
            do ik = 1,nks(ir)
               rad(ik,ir,iz) = powls(ik,ir,iz) * kareas(ik,ir)
               radtot(iz) = radtot(iz) + rad(ik,ir,iz)
               rad(ik,ir,nizs+1) = rad(ik,ir,nizs+1) + rad(ik,ir,iz)
               radtot(nizs+1) = radtot(nizs+1) + rad(ik,ir,iz)
               sdtot(ik,ir) = sdtot(ik,ir) + sdlims(ik,ir,iz)
            end do
         end do
      end do
c
c     Sort/order each ionization state and the total
c
c
      do iz = 0,nizs+1
         count(iz) = 0
         do ir = 1,nrs
            do ik = 1,nks(ir)
               if (rad(ik,ir,iz).ne.0.0) then
                  count(iz) = count(iz) + 1
                  radord(count(iz),iz,1) = rad(ik,ir,iz)
                  radord(count(iz),iz,2) = -1
                  radord(count(iz),iz,3) = -1
                  radord(count(iz),iz,4) = ik
                  radord(count(iz),iz,5) = ir
               endif
            end do
         end do
      end do
c
c     Now have unsorted array - need to sort it.
c
      call sortrad(radord,count,top,nizs)
c
c     Now that we have sorted arrays - can take the top 2/3 of the
c     total radiated power and calculated avearged densities.
c
      do iz = 0,nizs+1
c
         num = top(iz)
         bot = 0
c
         pradt(iz) = 0.0
         neav(iz) = 0.0
         teav(iz) = 0.0
         nizav(iz) = 0.0
         volav(iz) = 0.0
c
c         if (count(iz).eq.0) cycle
c
         if (count(iz).ne.0) then
c
c         do while(num.ne.-1
c     >               .and.pradt(iz).lt.0.66*radtot(iz))
c
 30       if (num.eq.-1.or.pradt(iz).ge.0.66*radtot(iz)) goto 20
c
            pradt(iz) = pradt(iz) + radord(num,iz,1)
c
            ik = radord(num,iz,4)
            ir = radord(num,iz,5)
c
            neav(iz) = neav(iz) + kareas(ik,ir) * knbs(ik,ir)
            teav(iz) = teav(iz) + kareas(ik,ir) * ktebs(ik,ir)
c
            if (iz.eq.nizs+1) then
c
               nizav(iz)= nizav(iz)+ kareas(ik,ir) * sdtot(ik,ir)
c
            else
               nizav(iz)= nizav(iz)+ kareas(ik,ir) * sdlims(ik,ir,iz)
            endif
c
            volav(iz)= volav(iz)+ kareas(ik,ir)
c
            bot = num
            num = radord(num,iz,2)
c
         goto 30
c
 20      continue
c
c         end do
c
c        Calculate contour level required for 2/3 of radiated power
c
         pradclev(iz) = radord(bot,iz,1)/radord(top(iz),iz,1)
c
c        Calculate averages
c
         if (volav(iz).ne.0.0) then
            neav(iz) = neav(iz) / volav(iz)
            nizav(iz)= nizav(iz) / volav(iz)
            teav(iz) = teav(iz) / volav(iz)
            lzav(iz) = pradt(iz)/(volav(iz)*neav(iz)*nizav(iz))
         else
            neav(iz) = 0.0
            nizav(iz)= 0.0
            teav(iz) = 0.0
            lzav(iz) = 0.0
         endif
c
c        Endif for count(iz).ne.0
c
         endif
      end do
c
c     Print out summary of data
c
      write (6,*)
      write (6,*) 'SUMMARY of Radiation Source: LOCAL Conditions'
      write (6,'(a)') job
      write (6,*)
      write (6,1020)
c
      do iz = 0,nizs+1
c
         if (iz.eq.nizs+1) then
            write(6,1010) radtot(iz),pradt(iz),neav(iz),
     >         nizav(iz)*absfac,
     >         nizav(iz),teav(iz),volav(iz),lzav(iz)
         else
            write(6,1000) iz,radtot(iz),pradt(iz),neav(iz),
     >         nizav(iz)*absfac,
     >         nizav(iz),teav(iz),volav(iz),lzav(iz)
         endif
c
      end do
c
 1020 format (5x,'Ion State',5x,'Total Rad',4x,'Cutoff Rad',
     >             9x,'Ne AV',8x,'Niz AV',7x,'Niz Raw',9x,'Te AV',
     >             8x,'Volume',8x,'Lz EST')
 1010 format (8x,'ALL',3x,8(1x,g13.5))
 1000 format (7x,i4,3x,8(1x,g13.5))
c
c
      return
      end
c
c
c
      subroutine calc_mfp(lgradti,lgradte,lmfpii,lmfpee)
      use mod_params
      use mod_cgeom
      implicit none
c     include 'params'
c     include 'cgeom'
c
      real lgradte(maxnks,maxnrs),lgradti(maxnks,maxnrs)
      real lmfpii(maxnks,maxnrs),lmfpee(maxnks,maxnrs)
c
c
c     CALC_MFP:
c
c     This routine calculates the electron and ion mean free paths
c     which are then used for plotting along the field lines.
c     It also includes the calaculation of the scale lengths for
c     both the electron and ion temperatures which are also included
c     on the plots.
c
c     David Elder, Nov 19, 1996
c
c
c     This routine includes code to calculate the scale lengths for
c     both density and static pressure - but these have been
c     commented out at this time.
c
      integer ik,ir,id
c
c      real gradnpara(maxnks,maxnrs),gradppara(maxnks,maxnrs)
c
      real gradtepara(maxnks,maxnrs),gradtipara(maxnks,maxnrs)
c
      call calc_grad(gradtepara,ktebs,kteds)
      call calc_grad(gradtipara,ktibs,ktids)
c
c      call calc_grad(gradnpara,knbs,knds)
c
c        Calculate the STATIC pressure in each cell and at the targets.
c
c         do ir = 1,nrs
c            do ik = 1,nks(ir)
c               kpbs(ik,ir) = knbs(ik,ir) *
c     >                          ech * (ktebs(ik,ir)+ktibs(ik,ir))
c            end do
c         end do
c
c         do id = 1,nds
c            kpds(id) = knds(id) * ech * (kteds(id)+ktids(id))
c         end do
c
c         call calc_grad(gradppara,kpbs,kpds)
c
C-----------------------------------------------------------------------
c
c        Calculate the scale lengths
c
      call calc_scale(gradtepara,lgradte,ktebs,nrs,nks)
      call calc_scale(gradtipara,lgradti,ktibs,nrs,nks)
c
c      call calc_scale(gradnpara,lgradn,knbs,nrs,nks)
c      call calc_scale(gradppara,lgradp,kpbs,nrs,nks)
c
c     Calculate the lmfpii and ee values
c
      do ir = 1,nrs
         do ik = 1,nks(ir)
c
c           Calculate the electron and ion mean free paths
c

            if (knbs(ik,ir).le.0.0) then
               lmfpii(ik,ir) = 0.0
               lmfpee(ik,ir) = 0.0
            else
c
               lmfpii(ik,ir)=ktibs(ik,ir)**2 * (1.5e16/knbs(ik,ir))
               lmfpee(ik,ir)=ktebs(ik,ir)**2 * (1.5e16/knbs(ik,ir))
c
            endif
c
         end do
      end do
c
      return
      end
c
c
c
      subroutine calc_grad(valgrad,val,valtarg)
      use mod_params
      use mod_cgeom
      implicit none
c     include 'params'
c     include 'cgeom'
c
      real valgrad(maxnks,maxnrs)
      real val(maxnks,maxnrs)
      real valtarg(maxnds)
      real endval
c
c     CALC_GRAD:
c
c     This subroutine calculates the parallel to the field
c     line gradients of the parameter passed in. It uses the
c     values at the targets passed in the valtarg array to
c     help define the gradients in the first cell.
c
c
      integer ik,ir,id
      real dist1, dist2
c
      do ir = 1, nrs
         do ik = 1,nks(ir)

            if (ik.eq.1) then
c
c              Distinguish between core and SOL
c
               dist2 = (kss(ik+1,ir)-kss(ik,ir))
c
               if (ir.lt.irsep) then
c
                  dist1 = kss(nks(ir),ir) - kss(nks(ir)-1,ir)
                  endval= val(nks(ir)-1,ir)
c
               else
c
c                 IK=1 target (second)
c
                  id = idds(ir,2)

                  dist1 = (kss(ik,ir)-ksb(ik-1,ir))
                  endval = valtarg(id)
c
               endif
c
               if (dist1.ne.0.0.and.dist2.ne.0.0) then
                  valgrad(ik,ir) = ((val(ik+1,ir)-val(ik,ir))
     >                             /dist2
     >                           + (val(ik,ir)-endval)
     >                             /dist1  )/2.0
               elseif (dist1.ne.0.0) then
                  valgrad(ik,ir) = (val(ik,ir)-endval)
     >                                /dist1
               elseif (dist2.ne.0.0) then
                  valgrad(ik,ir) = (val(ik+1,ir)-val(ik,ir))
     >                             /dist2
               else
                  valgrad(ik,ir) = 0.0
               endif
c
            elseif (ik.eq.nks(ir)) then
c
c              Distinguish between core and SOL
c
               dist1 = (kss(ik,ir)-kss(ik-1,ir))

               if (ir.lt.irsep) then
c
                  dist2 = (kss(2,ir) - kss(1,ir))
                  endval= val(2,ir)
c
               else
c
c                 IK=NKS(IR) target (first)
c
                  id = idds(ir,1)
c
                  dist2 = ksb(ik,ir)-kss(ik,ir)
                  endval = valtarg(id)
c
               endif
c
               if (dist1.ne.0.0.and.dist2.ne.0.0) then
                  valgrad(ik,ir) = ((endval-val(ik,ir))
     >                             /dist2
     >                           + (val(ik,ir)-val(ik-1,ir))
     >                             /dist1  )/2.0
               elseif (dist1.ne.0.0) then
                  valgrad(ik,ir) = (val(ik,ir)-val(ik-1,ir))
     >                                /dist1
               elseif (dist2.ne.0.0) then
                  valgrad(ik,ir) = (endval-val(ik,ir))
     >                             /dist2
               else
                  valgrad(ik,ir) = 0.0
               endif
c
            else
c
c              General case - give average of slope to next and last
c              cells.
c
               dist2 = (kss(ik+1,ir)-kss(ik,ir))
               dist1 = (kss(ik,ir)-kss(ik-1,ir))
c
               if (dist1.ne.0.0.and.dist2.ne.0.0) then
                  valgrad(ik,ir) = ((val(ik+1,ir)-val(ik,ir))
     >                             /dist2
     >                           + (val(ik,ir)-val(ik-1,ir))
     >                             /dist1  )/2.0
               elseif (dist1.ne.0.0) then
                  valgrad(ik,ir) = (val(ik,ir)-val(ik-1,ir))
     >                                /dist1
               elseif (dist2.ne.0.0) then
                  valgrad(ik,ir) = (val(ik+1,ir)-val(ik,ir))
     >                             /dist2
               else
                  valgrad(ik,ir) = 0.0
               endif
c
            endif
c
         end do
      end do
c
      return
      end
c
c
c
      subroutine calc_scale(valgrad,valscale,val,nrs,nks)
      use mod_params
      implicit none
c     include 'params'
      real valgrad(maxnks,maxnrs)
      real val(maxnks,maxnrs)
      real valscale(maxnks,maxnrs)
      integer nrs
      integer nks(maxnrs)
c
c     CALC_SCALE:
c
c     This subroutine calculates the cell by cell scale
c     lengths of the quantities passed in.
c
      integer ik,ir
c
      do ir = 1,nrs
         do ik = 1,nks(ir)
c
c           Calculate scale lengths - an error condition will
c           result in a very large scale length being assigned.
c
            if (valgrad(ik,ir).eq.0.0.or.val(ik,ir).eq.0.0) then
               valscale(ik,ir) = hi
            else
               valscale(ik,ir) = abs(
     >                       1.0/((1.0/val(ik,ir))*valgrad(ik,ir)))
            endif
         end do
      end do
c
      return
      end
c
c
C     Krieger IPP 12/96
C     define routine divkill as stub for entry point divkill in
C     DIVIMP. Needed to avoid complains by loader
c
      subroutine divkill
      implicit none
c
c     This routine is the forced EXIT point if OUT is sent a USR
c     KILL signal. It is called divkill for compatibility with
c     the USR1 kill signal code that is included in the sys
c     module and thus used in both DIVIMP and OUT
c
      write (6,*) 'USR1 kill signal received by OUT.'
      write (6,*) 'OUT execution aborted.'
      call prc('USR1 kill signal received by OUT.')
      call prc('OUT execution aborted.')
c
      stop
c
      end
C
C
C
      SUBROUTINE LDADAS(CZ,IZ,ADASID,ADASYR,ADASEX,ISELE,ISELR,ISELX,
     >                  CVALS,WAVE,IRCODE)
      use mod_params
      use mod_cgeom
      use mod_pindata
      use mod_dynam2
      use mod_comtor
      IMPLICIT NONE
C
C  *********************************************************************
C  *                                                                   *
C  *  LDADAS:  CODE TO LOAD THE REQUESTED LINE EMISSION PROFILE        *
C  *           INTO THE MATRIX CVALS.  THE LINE IS SPECIFIED BY THE    *
C  *           NUCLEAR CHARGE OF THE EMITTING ION, THE                 *
C  *           CHARGE OF THE EMITTING IONISATION STATE, AN ID FLAG     *
C  *           WHICH LOCATES THE INPUT FILE, AND THREE BLOCK SELECTOR  *
C  *           NUMBERS, ONE EACH FOR EMISSION BY ELECTRON EXCITATION,  *
C  *           RECOMBINATION FROM THE NEXT HIGHER IONISATION STATE,    *
C  *           AND CHARGE EXCHANGE FROM THE NEXT HIGHER IONISATION     *
C  *           STATE.  IN ADDITION TO THE EMISSION PROFILE, THE        *
C  *           ROUTINE RETURNS THE WAVELENGTH OF THE TRANSITION AND    *
C  *           AN ERROR CODE FROM THE ADAS EXTRACTION ROUTINE, SPEC.   *
C  *           NOTE THAT THERE IS NO CHECKING OF THE SELECTOR NUMBERS! *
C  *                                                                   *
C  *                                                                   *
C  *            LORNE HORTON   (JET)         SEPTEMBER 1993            *
C  *                                                                   *
C  *********************************************************************
C
c     include 'params'
c     INCLUDE 'cgeom'
c     include 'pindata'
c     include 'dynam2'
c     include 'comtor'
C
      INTEGER   CZ,IZ,ADASYR,ISELE,ISELR,ISELX,IRCODE
      REAL      WAVE,CVALS(MAXNKS,MAXNRS)
      CHARACTER ADASID*(*),ADASEX*(*)
C
      INTEGER   IR,IK,IADAS,NPAIRS,IKK
      REAL*8    TADAS(20),DADAS(20)
      REAL*8    WLNGTH,PECAE(20),PECAR(20),PECAX(20)
      LOGICAL*4 LTRNG(20),LDRNG(20)
      CHARACTER ADASGR*8,ADASTY*80,PECTITLE*120
      CHARACTER XFESYM*2
C
      call rzero (cvals,MAXNKS*MAXNRS)
c
      WAVE = 0.0
      IRCODE = 0
      CALL XXUID(ADASID)
      IF (ADASYR.GE.0) THEN
        ADASGR = 'pec??#'//XFESYM(CZ)
        WRITE(ADASGR(4:5),'(I2.2)') ADASYR
      ELSE
        ADASGR = '*'
      ENDIF
      ADASTY = '*'
      CALL XXSPEC(ADASGR,ADASTY,ADASEX)
C
      DO IR = 1,NRS
        DO IK = 1,NKS(IR),20
          NPAIRS = MIN0(20,NKS(IR)-(IK-1))
          DO IADAS = 1,NPAIRS
            TADAS(IADAS) = DBLE(KTEBS(IK+(IADAS-1),IR))
            DADAS(IADAS) = DBLE(1.E-6*RIZB*KNBS(IK+(IADAS-1),IR))
          ENDDO
C
          CALL DZERO(PECAE,NPAIRS)
          IF (ISELE.GT.0) THEN
            CALL SPEC(ISELE,IZ,CZ,NPAIRS,TADAS,DADAS,
     >                WLNGTH,PECAE,LTRNG,LDRNG,PECTITLE,IRCODE)
            IF (IRCODE.NE.0) RETURN
          ELSE IF (ISELE.EQ.-1) THEN
C
C  JUST LOAD EMISSION MEASURE FOR ISEL = -1
C    - SINCE THIS VALUE CAN EXCEED THE UNIX SINGLE PRECISION LIMIT,
C      WORK IN DENSITY UNITS OF 10**18
C
            CALL DINIT(PECAE,NPAIRS,1.D6*1.D-36)
            WLNGTH = 0.0
          ENDIF
C
          CALL DZERO(PECAR,NPAIRS)
          IF (ISELR.GT.0) THEN
            CALL SPEC(ISELR,IZ,CZ,NPAIRS,TADAS,DADAS,
     >                WLNGTH,PECAR,LTRNG,LDRNG,PECTITLE,IRCODE)
            IF (IRCODE.NE.0) RETURN
          ELSE IF (ISELR.EQ.-1) THEN
            CALL DINIT(PECAR,NPAIRS,1.D6*1.D-36)
            WLNGTH = 0.0
          ENDIF
C
C  FOR IMPURITIES USE THE THIRD SWITCH FOR CX, FOR HYDROGEN
C  ADD IN A MOLECULAR CONTRIBUTION INSTEAD
C
          CALL DZERO(PECAX,NPAIRS)
CLDH - USE ION TEMPERATURE FOR CX RATE
          IF (CZ.GT.1.0) THEN
            DO IADAS = 1,NPAIRS
              TADAS(IADAS) = DBLE(KTIBS(IK+(IADAS-1),IR))
            ENDDO
          ENDIF
CLDH
          IF (ISELX.GT.0) THEN
            CALL SPEC(ISELX,IZ,CZ,NPAIRS,TADAS,DADAS,
     >                WLNGTH,PECAX,LTRNG,LDRNG,PECTITLE,IRCODE)
            IF (IRCODE.NE.0) RETURN
          ELSE IF (ISELX.EQ.-1) THEN
            CALL DINIT(PECAX,NPAIRS,1.D6*1.D-36)
            WLNGTH = 0.0
          ENDIF
C
          DO IADAS = 1,NPAIRS
            IKK = IK + (IADAS-1)
c
c            write(6,'(a,3i5,10(1x,g13.5))') 'LDADAS:', ikk,ir,iz,
c     >         rizb,knbs(ikk,ir),ktebs(ikk,ir),sdlims(ikk,ir,iz),
c     >         sdlims(ikk,ir,iz+1),pinatom(ikk,ir),
c     >         pecae(iadas),pecar(iadas),pecax(iadas)
c
            IF (CZ.GT.1.0) THEN
              CVALS(IKK,IR) = 1.D-6*
     >        (PECAE(IADAS) * RIZB * KNBS(IKK,IR) * SDLIMS(IKK,IR,IZ)
     >        +PECAR(IADAS) * RIZB * KNBS(IKK,IR) * SDLIMS(IKK,IR,IZ+1)
     >        +PECAX(IADAS) * PINATOM(IKK,IR) * SDLIMS(IKK,IR,IZ+1))
c
c 
            ELSE
C
C---- HYDROGEN DENSITIES ARE IN DIFFERENT ARRAYS
C
              CVALS(IKK,IR) = 1.D-6*
     >        (PECAE(IADAS) * RIZB * KNBS(IKK,IR) * PINATOM(IKK,IR)
     >        +PECAR(IADAS) * RIZB * KNBS(IKK,IR) * KNBS(IKK,IR)
     >        +PECAX(IADAS) * KNBS(IKK,IR)   * PINMOL(IKK,IR))

c
c              write(6,'(a,3i5,8(1x,g12.5))') 'LDADAS:',ikk,ir,iadas,
c     >                pecae(iadas),pecar(iadas),pecax(iadas),
c     >                knbs(ikk,ir),pinatom(ikk,ir),pinmol(ikk,ir),
c     >                cvals(ikk,ir)
c

            ENDIF
          ENDDO
        ENDDO
      ENDDO
C
      WAVE = WLNGTH
C
      RETURN
      END




c
c
c
      SUBROUTINE LDADAS_RZ(CZ,IZ,ADASID,ADASYR,ADASEX,
     >                     ISELE,ISELR,ISELX,wave,ircode,
     >                     CVALS,exc_den,rec_den,RAXIS,ZAXIS,RPTS,ZPTS)
      use hc_get
      use mod_params
      IMPLICIT NONE
C
C  *********************************************************************
C  *                                                                   *
C  * LDADAS_RZ:CODE TO LOAD THE REQUESTED LINE EMISSION PROFILE        *
C  *           INTO THE MATRIX CVALS.  THE LINE IS SPECIFIED BY THE    *
C  *           NUCLEAR CHARGE OF THE EMITTING ION, THE                 *
C  *           CHARGE OF THE EMITTING IONISATION STATE, AN ID FLAG     *
C  *           WHICH LOCATES THE INPUT FILE, AND THREE BLOCK SELECTOR  *
C  *           NUMBERS, ONE EACH FOR EMISSION BY ELECTRON EXCITATION,  *
C  *           RECOMBINATION FROM THE NEXT HIGHER IONISATION STATE,    *
C  *           AND CHARGE EXCHANGE FROM THE NEXT HIGHER IONISATION     *
C  *           STATE.  IN ADDITION TO THE EMISSION PROFILE, THE        *
C  *           ROUTINE RETURNS THE WAVELENGTH OF THE TRANSITION AND    *
C  *           AN ERROR CODE FROM THE ADAS EXTRACTION ROUTINE, SPEC.   *
C  *           NOTE THAT THERE IS NO CHECKING OF THE SELECTOR NUMBERS! *
C  *                                                                   *
C  *                                                                   *
C  *            LORNE HORTON   (JET)         SEPTEMBER 1993            *
C  *                                                                   *
C  *********************************************************************
C
c     include 'params'
c      INCLUDE 'cgeom'
c      include 'pindata'
c      include 'dynam2'
c      include 'comtor'
C
      INTEGER   CZ,IZ,ADASYR,ISELE,ISELR,ISELX,IRCODE
      REAL      WAVE
      CHARACTER ADASID*(*),ADASEX*(*)
c
c     RZ ARRAY DATA
c
      integer rpts,zpts
      real cvals(rpts,zpts),raxis(rpts),zaxis(zpts)
      real exc_den(rpts,zpts),rec_den(rpts,zpts)
C
      INTEGER   IRIND,IZIND,IADAS,NPAIRS,IKK
      REAL*8    TADAS(20),DADAS(20)
      real*8    nh_data(20),nh_mol_data(20),ni_data(20)
      REAL*8    WLNGTH,PECAE(20),PECAR(20),PECAX(20)
      LOGICAL*4 LTRNG(20),LDRNG(20)
      CHARACTER ADASGR*8,ADASTY*80,PECTITLE*120
      CHARACTER XFESYM*2

      ! Other variables
      real :: rizb

      ! Background plasma data
      real :: ne,te,ti,vb,ef,nh,nh_mol
c
c     Zero data array
c
      cvals = 0.0
c
c     Initialize
c
      rizb = grizb()

      WAVE = 0.0
      IRCODE = 0

      CALL XXUID(ADASID)

      IF (ADASYR.GE.0) THEN
        ADASGR = 'pec??#'//XFESYM(CZ)
        WRITE(ADASGR(4:5),'(I2.2)') ADASYR
      ELSE
        ADASGR = '*'
      ENDIF
      ADASTY = '*'
      CALL XXSPEC(ADASGR,ADASTY,ADASEX)
C
c      write(0,*) 'LDADAS_RZ:',rpts,zpts,cz,iz
c
      DO IRIND = 1,RPTS
        DO IZIND = 1,ZPTS,20
          NPAIRS = MIN0(20,ZPTS-(IZIND-1))

c            write(6,*) 'Loading LDADAS:',irind,izind,npairs

          DO IADAS = 1,NPAIRS

            call get_plasma_rz(raxis(irind),zaxis(izind+iadas-1),
     >                         ne,te,ti,vb,ef,
     >                         nh,nh_mol)

            nh_data(iadas) = nh
            nh_mol_data(iadas) = nh_mol
            ni_data(iadas) = ne

            TADAS(IADAS) = DBLE(te)
            DADAS(IADAS) = DBLE(1.E-6*RIZB*ne)
c            
c            write(6,'(a,2i6,10g12.5)') 'PLASMA:',
c     >                irind,izind+iadas-1,raxis(irind),
c     >                zaxis(izind+iadas-1),ne,te,ti,vb,ef
c
          ENDDO
C
c         Load excitation rates
c
          CALL DZERO(PECAE,NPAIRS)
          IF (ISELE.GT.0) THEN
            CALL SPEC(ISELE,IZ,CZ,NPAIRS,TADAS,DADAS,
     >                WLNGTH,PECAE,LTRNG,LDRNG,PECTITLE,IRCODE)
            IF (IRCODE.NE.0) RETURN
          ELSE IF (ISELE.EQ.-1) THEN
C
C  JUST LOAD EMISSION MEASURE FOR ISEL = -1
C    - SINCE THIS VALUE CAN EXCEED THE UNIX SINGLE PRECISION LIMIT,
C      WORK IN DENSITY UNITS OF 10**18
C
            CALL DINIT(PECAE,NPAIRS,1.D6*1.D-36)
            WLNGTH = 0.0
          ENDIF
C
c         Load recombination rates
c
          CALL DZERO(PECAR,NPAIRS)
          IF (ISELR.GT.0) THEN
            CALL SPEC(ISELR,IZ,CZ,NPAIRS,TADAS,DADAS,
     >                WLNGTH,PECAR,LTRNG,LDRNG,PECTITLE,IRCODE)
            IF (IRCODE.NE.0) RETURN
          ELSE IF (ISELR.EQ.-1) THEN
            CALL DINIT(PECAR,NPAIRS,1.D6*1.D-36)
            WLNGTH = 0.0
          ENDIF
C
C  FOR IMPURITIES USE THE THIRD SWITCH FOR CX, FOR HYDROGEN
C  ADD IN A MOLECULAR CONTRIBUTION INSTEAD
C
          CALL DZERO(PECAX,NPAIRS)
CLDH - USE ION TEMPERATURE FOR CX RATE
          IF (CZ.GT.1.0) THEN
            DO IADAS = 1,NPAIRS
              TADAS(IADAS) = DBLE(ti)
            ENDDO
          ENDIF
CLDH
c
c         Load CX recombination rates
c
          IF (ISELX.GT.0) THEN
            CALL SPEC(ISELX,IZ,CZ,NPAIRS,TADAS,DADAS,
     >                WLNGTH,PECAX,LTRNG,LDRNG,PECTITLE,IRCODE)
            IF (IRCODE.NE.0) RETURN
          ELSE IF (ISELX.EQ.-1) THEN
            CALL DINIT(PECAX,NPAIRS,1.D6*1.D-36)
            WLNGTH = 0.0
          ENDIF
C
c         Calculate total emission for block of cells
c
          DO IADAS = 1,NPAIRS
            IKK = IZIND + (IADAS-1)
c
c            write(6,'(a,3i5,10(1x,g13.5))') 'LDADAS:', ikk,ir,iz,
c     >         rizb,knbs(ikk,ir),ktebs(ikk,ir),sdlims(ikk,ir,iz),
c     >         sdlims(ikk,ir,iz+1),pinatom(ikk,ir),
c     >         pecae(iadas),pecar(iadas),pecax(iadas)
c
            IF (CZ.GT.1.0) THEN

              CVALS(IRind,IKK) = 1.D-6*
     >        (PECAE(IADAS) * RIZB * ni_data(iadas) * exc_den(irind,ikk)
     >        +PECAR(IADAS) * RIZB * ni_data(iadas) * rec_den(irind,ikk) 
     >        +PECAX(IADAS) * nh_data(iadas) * rec_den(irind,ikk))

c
c           write(6,'(a,3i5,9(1x,g12.5))')'LDADAS_RZ:Z:',irind,ikk,iadas,
c     >             pecae(iadas),pecar(iadas),pecax(iadas),
c     >             ni_data(iadas),exc_den(irind,ikk),rec_den(irind,ikk),
c     >             cvals(irind,ikk),nh_data(iadas)
c

c              CVALS(IKK,IR) = 1.D-6*
c     >        (PECAE(IADAS) * RIZB * KNBS(IKK,IR) * SDLIMS(IKK,IR,IZ)
c     >        +PECAR(IADAS) * RIZB * KNBS(IKK,IR) * SDLIMS(IKK,IR,IZ+1)
c     >        +PECAX(IADAS) * PINATOM(IKK,IR) * SDLIMS(IKK,IR,IZ+1))
c
c 
            ELSE
C
C---- HYDROGEN DENSITIES ARE IN DIFFERENT ARRAYS
C
              CVALS(IRind,ikk) = 1.D-6*
     >        (PECAE(IADAS) * RIZB * ni_data(iadas) * nh_data(iadas)
     >        +PECAR(IADAS) * RIZB * ni_data(iadas) * ni_data(iadas)
     >        +PECAX(IADAS) * ni_data(iadas) * nh_mol_data(iadas))

c
c           write(6,'(a,3i5,8(1x,g12.5))')'LDADAS_RZ:H:',irind,ikk,iadas,
c     >                pecae(iadas),pecar(iadas),pecax(iadas),
c     >                ni_data(iadas),nh_data(iadas),nh_mol_data(iadas),
c     >                cvals(irind,ikk)
c

            ENDIF
c
          ENDDO
        ENDDO
      ENDDO
C
      WAVE = WLNGTH
C
      RETURN
      END

c
c
c
      SUBROUTINE LDADAS_RZ_EFF(CZ,IZ,ADASID,ADASYR,ADASEX,
     >                     ISELE,ISELR,ISELX,wave,ircode,
     >                     CVALS,exc_den,rec_den,RAXIS,ZAXIS,RPTS,ZPTS)
      use hc_get
      use mod_params
      IMPLICIT NONE
C
C  *********************************************************************
C  *                                                                   *
C  * LDADAS_RZ:CODE TO LOAD THE REQUESTED LINE EMISSION PROFILE        *
C  *           INTO THE MATRIX CVALS.  THE LINE IS SPECIFIED BY THE    *
C  *           NUCLEAR CHARGE OF THE EMITTING ION, THE                 *
C  *           CHARGE OF THE EMITTING IONISATION STATE, AN ID FLAG     *
C  *           WHICH LOCATES THE INPUT FILE, AND THREE BLOCK SELECTOR  *
C  *           NUMBERS, ONE EACH FOR EMISSION BY ELECTRON EXCITATION,  *
C  *           RECOMBINATION FROM THE NEXT HIGHER IONISATION STATE,    *
C  *           AND CHARGE EXCHANGE FROM THE NEXT HIGHER IONISATION     *
C  *           STATE.  IN ADDITION TO THE EMISSION PROFILE, THE        *
C  *           ROUTINE RETURNS THE WAVELENGTH OF THE TRANSITION AND    *
C  *           AN ERROR CODE FROM THE ADAS EXTRACTION ROUTINE, SPEC.   *
C  *           NOTE THAT THERE IS NO CHECKING OF THE SELECTOR NUMBERS! *
C  *                                                                   *
C  *           Attempt to improve code efficiency for sparse arrays    *
C  *                                                                   *
C  *            LORNE HORTON   (JET)         SEPTEMBER 1993            *
C  *                                                                   *
C  *********************************************************************
C
c     include 'params'
c      INCLUDE 'cgeom'
c      include 'pindata'
c      include 'dynam2'
c      include 'comtor'
C
      INTEGER   CZ,IZ,ADASYR,ISELE,ISELR,ISELX,IRCODE
      REAL      WAVE
      CHARACTER ADASID*(*),ADASEX*(*)
c
c     RZ ARRAY DATA
c
      integer rpts,zpts
      real cvals(rpts,zpts),raxis(rpts),zaxis(zpts)
      real exc_den(rpts,zpts),rec_den(rpts,zpts)
C
      INTEGER   IRIND,IZIND,IADAS,NPAIRS,IKK
      REAL*8    TADAS,DADAS
      REAL*8    WLNGTH,PECAE,PECAR,PECAX
      LOGICAL*4 LTRNG,LDRNG
      CHARACTER ADASGR*8,ADASTY*80,PECTITLE*120
      CHARACTER XFESYM*2

      ! Other variables
      real :: rizb

      ! Background plasma data
      real :: ne,te,ti,vb,ef,nh,nh_mol
c
c     Zero data array
c
      cvals = 0.0
c
c     Initialize
c
      rizb = grizb()

      WAVE = 0.0
      IRCODE = 0

      CALL XXUID(ADASID)

      IF (ADASYR.GE.0) THEN
        ADASGR = 'pec??#'//XFESYM(CZ)
        WRITE(ADASGR(4:5),'(I2.2)') ADASYR
      ELSE
        ADASGR = '*'
      ENDIF
      ADASTY = '*'
      CALL XXSPEC(ADASGR,ADASTY,ADASEX)
C
c      write(0,*) 'LDADAS_RZ:',rpts,zpts,cz,iz
c

      DO IRIND = 1,RPTS
        DO IZIND = 1,ZPTS

           if (cz.gt.1.and.exc_den(irind,izind).le.0.0.and.
     >                     rec_den(irind,izind).le.0.0) cycle
c
c          Get plasma in local cell
c           
           call get_plasma_rz(raxis(irind),zaxis(izind+iadas-1),
     >                         ne,te,ti,vb,ef,
     >                         nh,nh_mol)

           if (cz.eq.1.and.nh.eq.0.0.and.nh_mol.eq.0.0) cycle

           tadas = dble(te)
           dadas = dble(1.0e-6*rizb*ne)

C
c         Load excitation rates
c
           pecae = 0.0

          IF (ISELE.GT.0) THEN
            CALL SPEC(ISELE,IZ,CZ,1,TADAS,DADAS,
     >                WLNGTH,PECAE,LTRNG,LDRNG,PECTITLE,IRCODE)
            IF (IRCODE.NE.0) RETURN
          ELSE IF (ISELE.EQ.-1) THEN
C
C  JUST LOAD EMISSION MEASURE FOR ISEL = -1
C    - SINCE THIS VALUE CAN EXCEED THE UNIX SINGLE PRECISION LIMIT,
C      WORK IN DENSITY UNITS OF 10**18
C
             pecae = 1.D6*1.D-36
             WLNGTH = 0.0
          ENDIF
C
c         Load recombination rates
c
          pecar = 0.0
          IF (ISELR.GT.0) THEN
            CALL SPEC(ISELR,IZ,CZ,1,TADAS,DADAS,
     >                WLNGTH,PECAR,LTRNG,LDRNG,PECTITLE,IRCODE)
            IF (IRCODE.NE.0) RETURN
          ELSE IF (ISELR.EQ.-1) THEN
            pecar = 1.D6*1.D-36
            WLNGTH = 0.0
          ENDIF
C
C  FOR IMPURITIES USE THE THIRD SWITCH FOR CX, FOR HYDROGEN
C  ADD IN A MOLECULAR CONTRIBUTION INSTEAD
C
          pecax = 0.0
CLDH - USE ION TEMPERATURE FOR CX RATE
          tadas = dble(ti)
CLDH
c
c         Load CX recombination rates
c
          IF (ISELX.GT.0) THEN
            CALL SPEC(ISELX,IZ,CZ,1,TADAS,DADAS,
     >                WLNGTH,PECAX,LTRNG,LDRNG,PECTITLE,IRCODE)
            IF (IRCODE.NE.0) RETURN
          ELSE IF (ISELX.EQ.-1) THEN
             pecax = 1.D6*1.D-36
            WLNGTH = 0.0
          ENDIF


c         Calculate total emission for block of cells
c
c
            IF (CZ.GT.1.0) THEN

              CVALS(IRind,izind) = 1.D-6*
     >        (PECAE * RIZB * ne * exc_den(irind,izind)
     >        +PECAR * RIZB * ne * rec_den(irind,izind) 
     >        +PECAX *        nh * rec_den(irind,izind))

c 
            ELSE
C
C---- HYDROGEN DENSITIES ARE IN DIFFERENT ARRAYS
C
              CVALS(IRind,izind) = 1.D-6*
     >        (PECAE * RIZB * ne * nh
     >        +PECAR * RIZB * ne * ne 
     >        +PECAX * ne * nh_mol)

c
c

            ENDIF
c
        ENDDO
      ENDDO
C
      WAVE = WLNGTH
C
      RETURN
      END


C
C
C
      SUBROUTINE LDADAS_TIMEDEP(CZ,IZ,IT,ADASID,ADASYR,ADASEX,
     >                  ISELE,ISELR,ISELX,
     >                  CVALS,WAVE,IRCODE)
      use mod_params
      use mod_cgeom
      use mod_pindata
      use mod_dynam4
      use mod_comtor
      IMPLICIT NONE
C
C  *********************************************************************
C  *                                                                   *
C  *  LDADAS:  CODE TO LOAD THE REQUESTED LINE EMISSION PROFILE        *
C  *           INTO THE MATRIX CVALS.  THE LINE IS SPECIFIED BY THE    *
C  *           NUCLEAR CHARGE OF THE EMITTING ION, THE                 *
C  *           CHARGE OF THE EMITTING IONISATION STATE, AN ID FLAG     *
C  *           WHICH LOCATES THE INPUT FILE, AND THREE BLOCK SELECTOR  *
C  *           NUMBERS, ONE EACH FOR EMISSION BY ELECTRON EXCITATION,  *
C  *           RECOMBINATION FROM THE NEXT HIGHER IONISATION STATE,    *
C  *           AND CHARGE EXCHANGE FROM THE NEXT HIGHER IONISATION     *
C  *           STATE.  IN ADDITION TO THE EMISSION PROFILE, THE        *
C  *           ROUTINE RETURNS THE WAVELENGTH OF THE TRANSITION AND    *
C  *           AN ERROR CODE FROM THE ADAS EXTRACTION ROUTINE, SPEC.   *
C  *           NOTE THAT THERE IS NO CHECKING OF THE SELECTOR NUMBERS! *
c  *                                                                   *
c  *           THIS CODE HAS BEEN MODIFIED TO SUPPORT TIME DEPENDENT   *
C  *           DATA ARRAYS. IN ADDITION IF AN ION WITH ONE CHARGE      *
C  *           IS SPECIFIED - THE CODE STILL USES THE DIVIMP IMPURITY  *
C  *           DATA ASSUMING THAT DIVIMP FOLLOWED HYDROGEN AS AN       *
C  *           SINCE TIME DEPENDENT DATA IS NOT AVAILABLE FROM THE     *
C  *           HYDROGENIC CODES AT THIS TIME.                          *
C  *                                                                   *
C  *                                                                   *
C  *            LORNE HORTON   (JET)         SEPTEMBER 1993            *
C  *            DAVID ELDER    (TORONTO)     MAY 1999                  *
C  *                                                                   *
C  *********************************************************************
C
c     include 'params'
c     INCLUDE 'cgeom'
c     include 'pindata'
c     include 'dynam4'
c     include 'comtor'
C
      INTEGER   CZ,IZ,IT,ADASYR,ISELE,ISELR,ISELX,IRCODE
      REAL      WAVE,CVALS(MAXNKS,MAXNRS)
      CHARACTER ADASID*(*),ADASEX*(*)
C
      INTEGER   IR,IK,IADAS,NPAIRS,IKK
      REAL*8    TADAS(20),DADAS(20)
      REAL*8    WLNGTH,PECAE(20),PECAR(20),PECAX(20)
      LOGICAL*4 LTRNG(20),LDRNG(20)
      CHARACTER ADASGR*8,ADASTY*80,PECTITLE*120
      CHARACTER XFESYM*2
C
      WAVE = 0.0
      IRCODE = 0
      CALL XXUID(ADASID)
      IF (ADASYR.GE.0) THEN
        ADASGR = 'pec??#'//XFESYM(CZ)
        WRITE(ADASGR(4:5),'(I2.2)') ADASYR
      ELSE
        ADASGR = '*'
      ENDIF
      ADASTY = '*'
      CALL XXSPEC(ADASGR,ADASTY,ADASEX)
C
      DO IR = 1,NRS
        DO IK = 1,NKS(IR),20
          NPAIRS = MIN0(20,NKS(IR)-(IK-1))
          DO IADAS = 1,NPAIRS
            TADAS(IADAS) = DBLE(KTEBS(IK+(IADAS-1),IR))
            DADAS(IADAS) = DBLE(1.E-6*RIZB*KNBS(IK+(IADAS-1),IR))
          ENDDO
C
          CALL DZERO(PECAE,NPAIRS)
          IF (ISELE.GT.0) THEN
            CALL SPEC(ISELE,IZ,CZ,NPAIRS,TADAS,DADAS,
     >                WLNGTH,PECAE,LTRNG,LDRNG,PECTITLE,IRCODE)
            IF (IRCODE.NE.0) RETURN
          ELSE IF (ISELE.EQ.-1) THEN
C
C  JUST LOAD EMISSION MEASURE FOR ISEL = -1
C    - SINCE THIS VALUE CAN EXCEED THE UNIX SINGLE PRECISION LIMIT,
C      WORK IN DENSITY UNITS OF 10**18
C
            CALL DINIT(PECAE,NPAIRS,1.D6*1.D-36)
            WLNGTH = 0.0
          ENDIF
C
          CALL DZERO(PECAR,NPAIRS)
          IF (ISELR.GT.0) THEN
            CALL SPEC(ISELR,IZ,CZ,NPAIRS,TADAS,DADAS,
     >                WLNGTH,PECAR,LTRNG,LDRNG,PECTITLE,IRCODE)
            IF (IRCODE.NE.0) RETURN
          ELSE IF (ISELR.EQ.-1) THEN
            CALL DINIT(PECAR,NPAIRS,1.D6*1.D-36)
            WLNGTH = 0.0
          ENDIF
C
C  FOR IMPURITIES USE THE THIRD SWITCH FOR CX, FOR HYDROGEN
C  ADD IN A MOLECULAR CONTRIBUTION INSTEAD
C
          CALL DZERO(PECAX,NPAIRS)
CLDH - USE ION TEMPERATURE FOR CX RATE
          IF (CZ.GT.1.0) THEN
            DO IADAS = 1,NPAIRS
              TADAS(IADAS) = DBLE(KTIBS(IK+(IADAS-1),IR))
            ENDDO
          ENDIF
CLDH
          IF (ISELX.GT.0) THEN
            CALL SPEC(ISELX,IZ,CZ,NPAIRS,TADAS,DADAS,
     >                WLNGTH,PECAX,LTRNG,LDRNG,PECTITLE,IRCODE)
            IF (IRCODE.NE.0) RETURN
          ELSE IF (ISELX.EQ.-1) THEN
            CALL DINIT(PECAX,NPAIRS,1.D6*1.D-36)
            WLNGTH = 0.0
          ENDIF
C
          DO IADAS = 1,NPAIRS
            IKK = IK + (IADAS-1)
C
              CVALS(IKK,IR) = 1.D-6*
     >        (PECAE(IADAS) * RIZB * KNBS(IKK,IR)*LIMS(IKK,IR,IZ,it)
     >        +PECAR(IADAS) * RIZB * KNBS(IKK,IR)*LIMS(IKK,IR,IZ+1,it)
     >        +PECAX(IADAS) * PINATOM(IKK,IR) * LIMS(IKK,IR,IZ+1,it))
C
          ENDDO
        ENDDO
      ENDDO
C
      WAVE = WLNGTH
C
      RETURN
      END
c
C
C     OLD BREM CODE
C
c      SUBROUTINE LDBREM(WAVE,CVALS,IRCODE,NIZS)
c      IMPLICIT NONE
C
C  *********************************************************************
C  *                                                                   *
C  *  LDBREM:  CODE TO LOAD THE BREMSSTRAHLUNG EMISSION PROFILE        *
C  *           INTO THE MATRIX CVALS.  THE WLNGTH OF THE MEASUREMENT   *
C  *           IS SPECIFIED AS INPUT.                                  *
C  *                                                                   *
C  *                                                                   *
C  *            LORNE HORTON   (JET)         SEPTEMBER 1994            *
C  *                                                                   *
C  *********************************************************************
C
c      include 'params'
c      include 'cgeom'
c      include 'dynam2'
c      include 'comtor'
C
c      INTEGER   IRCODE,NIZS
c      REAL      WAVE,CVALS(MAXNKS,MAXNRS)
C
c      INTEGER   IR,IK,IZ
c      REAL      ZEFF, U, GFF, ZNUM, ZDEN
c      REAL*8    GIIIAV, U8, GAM2
C
c      DO IR = 1,NRS
c        DO IK = 1,NKS(IR)
C  ZEFF:
c          ZNUM = CIZB*CIZB*KNBS(IK,IR)
c          ZDEN =      CIZB*KNBS(IK,IR)
c          DO IZ = 1,NIZS
c            ZNUM = ZNUM + IZ*IZ*SDLIMS(IK,IR,IZ)
c            ZDEN = ZDEN +    IZ*SDLIMS(IK,IR,IZ)
c          ENDDO
c          ZEFF = ZNUM/ZDEN
C  REDUCED TRANSITION ENERGY
c          U = 12398./(WAVE*KTEBS(IK,IR))
C  GAUNT FACTOR:
c          U8 = U
c          GAM2 = ZEFF*ZEFF*13.6058/KTEBS(IK,IR)
c          GFF = GIIIAV(U8,GAM2)
C
c          CVALS(IK,IR) = 9.57E-29*RIZB*KNBS(IK,IR)*RIZB*KNBS(IK,IR)*
c     >             ZEFF*GFF*EXP(-U)/(1.E-10*WAVE*SQRT(KTEBS(IK,IR)))
c        ENDDO
c      ENDDO
C
c      RETURN
c      END
c
      SUBROUTINE LDBREM(WAVE,CVALS,IRCODE,nizs)
      use mod_params
      use mod_cgeom
      use mod_dynam2
      use mod_comtor
      IMPLICIT NONE
C
C  *********************************************************************
C  *                                                                   *
C  *  LDBREM:  CODE TO LOAD THE BREMSSTRAHLUNG EMISSION PROFILE        *
C  *           INTO THE MATRIX CVALS.  THE WLNGTH OF THE MEASUREMENT   *
C  *           IS SPECIFIED AS INPUT.                                  *
C  *                                                                   *
C  *                                                                   *
C  *            LORNE HORTON   (JET)         SEPTEMBER 1994            *
C  *                                                                   *
C  *            UPDATED 07/07/97  TO INCLUDE FREE-BOUND EMISSION       *
C  *                              AND A BETTER ROUTINE FOR INCLUDING   *
C  *                              THE IMPURITY CONTRIBUTIONS           *
C  *                                                                   *
C  *********************************************************************
C
c      INCLUDE (PPPARA)
c      INCLUDE (PPUNIT)
c      INCLUDE (PPGEOM)
c      INCLUDE (PPPLAS)
c
c     include 'params'
c     include 'cgeom'
c     include 'dynam2'
c     include 'comtor'
C
      INTEGER   IRCODE,nizs
      REAL      WAVE,CVALS(MAXNKS,MAXNRS)
C
      INTEGER   IR, IK, IZ, J
c
      integer   maxnzs,nzs
      parameter (maxnzs=1,nzs=1)
c
      REAL*8    CONTIN, WAVE8, NE, TEV, NH, NZ(0:MAXIZS,MAXNZS)
C
      call dzero(nz,(maxizs+1)*maxnzs)
c
      WAVE8 = WAVE
c
c      write (0,'(a,2i4,2x,g18.7)') 'BREM: NIZS=',nizs,cion,absfac
c
      DO IR = 1,NRS
        DO IK = 1,NKS(IR)
          NE = KNBS(IK,IR)*1.0E-6
          TEV = KTEBS(IK,IR)
          NH = NE
c
c          DO J = 1, NZS
c
          if (nizs.gt.0) then
c
             J = 1

             DO IZ = 0, cion
c
c              IF (IZ.LE.NIZS(J)) THEN
c
c                NZ(IZ,J) = SDLIMS(IK,IR,IZ,J)*1.0E-6
c
c
                NZ(IZ,J) = SDLIMS(IK,IR,IZ)*1.0E-6 * ABSFAC
c
c
c              ELSE
c                NZ(IZ,J) = 0.0
c              ENDIF
c
                NH = NH - IZ*NZ(IZ,J)
c
             ENDDO
c
          endif
c
c          ENDDO
C
C  REMEMBER TO CONVERT TO M**-3 AND NM**-1
C
          CVALS(IK,IR) = CONTIN(WAVE8, NE, TEV, NZS, MAXIZS, CION,
     >                          NH, NZ)*1.0D7
        ENDDO
      ENDDO
C
      RETURN
      END
c
c
c
      SUBROUTINE LDBREM_SPEC(WAVE,npairs,den,tbrem,
     >                       brempec,IRCODE)
      use mod_params
      use mod_cgeom
      use mod_dynam2
      use mod_comtor
      IMPLICIT NONE
C
C  *********************************************************************
C  *                                                                   *
C  *  LDBREM:  CODE TO LOAD THE BREMSSTRAHLUNG EMISSION PROFILE        *
C  *           INTO THE MATRIX CVALS.  THE WLNGTH OF THE MEASUREMENT   *
C  *           IS SPECIFIED AS INPUT.                                  *
C  *                                                                   *
C  *                                                                   *
C  *            LORNE HORTON   (JET)         SEPTEMBER 1994            *
C  *                                                                   *
C  *            UPDATED 07/07/97  TO INCLUDE FREE-BOUND EMISSION       *
C  *                              AND A BETTER ROUTINE FOR INCLUDING   *
C  *                              THE IMPURITY CONTRIBUTIONS           *
C  *                                                                   *
C  *********************************************************************
C
c      INCLUDE (PPPARA)
c      INCLUDE (PPUNIT)
c      INCLUDE (PPGEOM)
c      INCLUDE (PPPLAS)
c
c     include 'params'
c     include 'cgeom'
c     include 'dynam2'
c     include 'comtor'
C
      INTEGER   IRCODE,npairs
      REAL      WAVE,CVALS(MAXNKS,MAXNRS)
      real      tbrem(npairs),brempec(npairs),den
C
      INTEGER   IR, IK, IZ, J
c
      integer   maxnzs,nzs
      parameter (maxnzs=1,nzs=1)
c
      REAL*8    CONTIN, WAVE8, NE, TEV, NH, NZ(0:MAXIZS,MAXNZS)
C

c
      j=1
      ircode = 0
      WAVE8 = WAVE
c
c      write (0,'(a,2i4,2x,g18.7)') 'BREM: NIZS=',nizs,cion,absfac
c
      DO IR = 1,npairs
         NE = den * 1.0e-6
         TEV = tbrem(ir)
         NH = NE
         DO IZ = 0, cion
            NZ(IZ,J) = 0.0
         END DO
C
C  REMEMBER TO CONVERT TO M**-3 AND NM**-1
C
         brempec(IR) = CONTIN(WAVE8, NE, TEV, NZS, MAXIZS, CION,
     >                          NH, NZ)*1.0D7
c
      ENDDO
C
      RETURN
      END
C
C
C
      SUBROUTINE LDBRFF(WAVE,CVALS,IRCODE,nizs)
      use mod_params
      use mod_cgeom
      use mod_dynam2
      use mod_comtor
      IMPLICIT NONE
C
C  *********************************************************************
C  *                                                                   *
C  *  LDBRFF:  CODE TO LOAD THE FREE-FREE BREMSSTRAHLUNG EMISSION      *
C  *           PROFILE INTO THE MATRIX CVALS.  THE WLNGTH OF THE       *
C  *           MEASUREMENT IS SPECIFIED AS INPUT.  THIS IS A SIMPLE    *
C  *           MOD TO LDBREM WHICH CALLS CONTFF INSTEAD OF CONTIN AND  *
C  *           IS MOSTLY FOR TESTING RATHER THAN FOR DIAGNOSTIC        *
C  *           SIMULATION.                                             *
C  *                                                                   *
C  *            LORNE HORTON   (JET)         JULY 1997                 *
C  *                                                                   *
C  *********************************************************************
C
c      INCLUDE (PPPARA)
c      INCLUDE (PPUNIT)
c      INCLUDE (PPGEOM)
c      INCLUDE (PPPLAS)
c
c     include 'params'
c     include 'cgeom'
c     include 'dynam2'
c     include 'comtor'
C
      integer nzs,maxnzs
      parameter(nzs=1,maxnzs=1)
c
      INTEGER   IRCODE,nizs
      REAL      WAVE,CVALS(MAXNKS,MAXNRS)
C
      INTEGER   IR, IK, IZ, J
      REAL*8    CONTFF, WAVE8, NE, TEV, NH, NZ(0:MAXIZS,MAXNZS)
C
      call dzero(nz,(maxizs+1)*maxnzs)
c
      WAVE8 = WAVE
      DO IR = 1,NRS
        DO IK = 1,NKS(IR)
          NE = KNBS(IK,IR)*1.0E-6
          TEV = KTEBS(IK,IR)
          NH = NE
c
          if (nizs.gt.0) then
c
             J = NZS
c
             DO IZ = 0, CION
c
                NZ(IZ,J) = SDLIMS(IK,IR,IZ)*1.0E-6*ABSFAC
c
                NH = NH - IZ*NZ(IZ,J)
c
             ENDDO
c
          endif
C
C  REMEMBER TO CONVERT TO M**-3 AND NM**-1
C
          CVALS(IK,IR) = CONTFF(WAVE8, NE, TEV, NZS, MAXIZS, CION,
     >                          NH, NZ)*1.0D7
        ENDDO
      ENDDO
C
      RETURN
      END
c
c
c
      subroutine ldpec (den,tadas,dadas,npairs,cz,iz,isel,peca,
     >                  ircode,wave,pectitle)
      implicit none
      real den,wave
      integer npairs,cz,iz,isel,ircode
      real*8 tadas(npairs),dadas(npairs),peca(npairs)
      character*(*) pectitle
c
c     LDPEC: Subroutine to load just one PEC set at a specific
c            density over a range of temperatures.
c
      integer iadas
      real*8  wlngth
      LOGICAL*4 LTRNG(20),LDRNG(20)
c
      wave = 0.0
c
      if (den.ne.0.0) then
         do iadas = 1,npairs
            dadas(iadas) = den * 1.0e-6
         end do
      end if
c
      CALL DZERO(PECA,NPAIRS)
c
      CALL SPEC(ISEL,IZ,CZ,NPAIRS,TADAS,DADAS,
     >          WLNGTH,PECA,LTRNG,LDRNG,PECTITLE,IRCODE)
c
      IF (IRCODE.NE.0) RETURN
c
      wave = wlngth
c
      return
      end
c
c
c
      subroutine calc_expt(iseld,touts,tvals,maxnthe,numthe,
     >                     themin,themax,maxnngs,ngs,datatitle)
      use mod_params
      implicit none
c
c     include 'params'
c
      integer iseld,ngs,numthe,maxnthe,maxnngs
      real themin,themax,touts(maxnthe),tvals(maxnthe,maxnngs)
      character*(*) datatitle
c
c     CALC_EXPT: This subroutine loads the selected
c                experimental data set from disk into
c                the array tvals at the position ngs. It
c                interpolates the experimental data to
c                obtain values at the same points as will
c                be plotted for the DIVIMP data - this
c                was done to enable proper scaling of all
c                the data when the plots are being drawn.
c
c                This data is then loaded into the array
c                TVALS for passing to DRAW. NGS is the
c                plot index to be used by the first set of
c                experimental data.
c
c
c                David Elder, November 24, 1997
c
c
c     Local variables
c
      integer dataunit,maxcols
      parameter (dataunit=13,maxcols=1)
c
c      integer maxdatx,dataunit,maxcols
c      parameter (maxdatx=1000,dataunit=13,maxcols=1)
c
      integer axis_type,in,cur_index,num_expt,ipos,ncols
      integer ind1,ind2
c
c
c     Experimental data
c
      real expt_axis(maxdatx),expt_data(maxdatx),theta
c
c     Write out input -
c
      write(6,*) 'Input:',iseld,maxnthe,numthe,maxnngs,ngs
c
      call load_expt_data(dataunit,iseld,expt_axis,axis_type,expt_data,
     >                    maxcols,maxdatx,num_expt,ncols,
     >                    datatitle)
c
      if (num_expt.le.0) then
c
         write(6,*) 'ERROR IN EXPERIMENTAL DATA: CAN NOT'//
     >                    ' LOAD DATASET # ',
     >                    iseld, ' - NO ELEMENTS FOUND'
         write(0,*) 'ERROR IN EXPERIMENTAL DATA: CAN NOT'//
     >                    ' LOAD DATASET # ',
     >                    iseld, ' - NO ELEMENTS FOUND'
c
         ngs = ngs - 1
c
         return
c
      endif
c
c     Check for single point LOS comparisons
c
      if (numthe.eq.1.and.num_expt.eq.1) then
c
         tvals(1,ngs) = expt_data(1)
         write (6,*) 'Calc_expt: SINGLE:',tvals(1,ngs)
c
      else
c
         do in = 1, numthe
c
            theta = touts(in)
c
            call arrpos(theta,expt_axis,num_expt,ind1,ind2)
c
            if (ind1.eq.-1) then
c
c              Index below position of first experimental data point
c              Set to zero ...
c
               tvals(in,ngs) = 0.0
c
            elseif (ind2.eq.-1) then
c
c              Out of bounds at top of range - again set to zero
c
               tvals(in,ngs) = 0.0
c
            else
c
c              Value in range - perform linear interpolation.
c
               tvals(in,ngs) = expt_data(ind1) +
     >                     (expt_data(ind2)-expt_data(ind1))
     >                    *(theta-expt_axis(ind1))
     >                    /(expt_axis(ind2)-expt_axis(ind1))
c
            endif
c
c           Write out debugging information to start ...
c
c slmod begin - new
            if (ind1.ne.-1.and.ind2.ne.-2) then
              write (6,'(a,3i4,4(1x,g12.5),a,
     >                    3(1x,g12.5))') 'Calc_expt:',in,ind1,ind2,
     >              expt_data(ind1),tvals(in,ngs),
     >              expt_data(ind2),tvals(in,1),'|',
     >              expt_axis(ind1), theta,expt_axis(ind2)
            else
              write (6,'(a)')
     >          'Calc_expt: ind1 or ind2 = -1'
            endif
c
c            write (6,'(a,2i4,4(1x,g12.5),a,
c     >                    3(1x,g12.5))') 'Calc_expt:',in,ind1,ind2,
c     >              expt_data(ind1),tvals(in,ngs),
c     >              expt_data(ind2),tvals(in,1),'|',
c     >              expt_axis(ind1), theta,expt_axis(ind2)
c slmod end
c
c
c
c            if (cur_index.eq.1) then
c
c              Index below position of first experimental data point
c              Set to zero ...
c
c               tvals(in,ngs) = 0.0
c
c            elseif (in.eq.num_expt.and.theta.gt.expt_axis(num_expt))
c     >               then
c
c              Out of bounds at top of range - again set to zero
c
c               tvals(in,ngs) = 0.0
c
c            else
c
c              Value in range - perform linear interpolation.
c
c               tvals(in,ngs) = expt_data(cur_index-1) +
c     >                     (expt_data(cur_index)-expt_data(cur_index-1))
c     >                    *(theta-expt_axis(cur_index-1))
c     >                    /(expt_axis(cur_index)-expt_axis(cur_index-1))
c
c            endif
c
c           Write out debugging information to start ...
c
c slmod begin - new
c            IF (cur_index.GT.1)
c     .        write (6,'(a,2i4,4(1x,g12.5),a,
c
c            write (6,'(a,2i4,4(1x,g12.5),a,
c slmod end
c     >                    3(1x,g12.5))') 'Calc_expt:',in,cur_index,
c     >              expt_data(cur_index-1),tvals(in,ngs),
c     >              expt_data(cur_index),tvals(in,1),'|',
c     >              expt_axis(cur_index-1), theta,expt_axis(cur_index)
c
         end do
c
      endif
c
c     Exit
c
      return
c
      end
c
c
c
      subroutine arrpos(val,vals,nvals,ind1,ind2)
      implicit none
c
      integer nvals,ind1,ind2 
      real val,vals(nvals) 
c
c     ARRPOS: This routine determines the indices of the data in the
c             array vals which are the upper and lower bounds of the
c             input value. The vals array must be ordered - either
c             ascending or descending is fine.  
c
      integer in,ipos,jpos
      external ipos,jpos
c
c
c     Check boundaries as special cases
c
      if (val.eq.vals(1)) then 
         ind1 = 1
         ind2 = 2
         return
      endif
c
      if (val.eq.vals(nvals)) then 
         ind1 = nvals-1
         ind2 = nvals
         return
      endif
c
c     Find specific cell
c
      if (vals(1).lt.vals(nvals)) then 
         in = ipos(val,vals,nvals)   
         if (in.eq.1.and.val.lt.vals(1)) then 
            ind1 = -1
            ind2 =  1 
         elseif (in.eq.nvals.and.val.gt.vals(nvals)) then 
            ind1 = nvals
            ind2 = -1 
         else
            ind1 = in -1
            ind2 = in  
         endif 

      else
         in = jpos(val,vals,nvals)
         if (in.eq.1.and.val.gt.vals(1)) then 
            ind1 = -1
            ind2 =  1 
         elseif (in.eq.nvals.and.val.lt.vals(nvals)) then 
            ind1 = nvals
            ind2 = -1 
         else
            ind1 = in 
            ind2 = in +1
         endif 


      endif
c 
      return
c
      end   
c
c
c

c      subroutine load_expt_data(dataunit,iseld,expt_axis,axis_type,
c     >                    expt_data,maxcols,maxdata,
c     >                    num_expt,ncols,datatitle)
c slmod begin - new
c
c Input:
c ISELD         - data index number
c MAXCOLS       - maximum number of data columns
c MAXDATA       - maximum number of data items
c
c Output:
c EXPT_AXIS     - independent data
c AXIS_TYPE     - read from file 13 (optional listing), default is 1
c EXPT_DATA     - MAXDATA,MAXCOLS dependent data
c NUM_EXPT      - number of data items in the data block, default is 0
c NCOLS         - number of data columns in data block, default is 1
c DATATITLE     - title, default is 'NO TITLE'
c
c slmod end
c      implicit none
c      include 'params'
c      integer dataunit,axis_type,num_expt,iseld,maxdata,colindex
c      integer maxcols,ncols
c      character*(*) datatitle
c      real expt_data(maxdata,maxcols),expt_axis(maxdata)
c
c     LOAD_EXPT_DATA: This routine loads the experimental data
c                     into the local arrays and passes
c                     back the relevant information.
c
c
c     Local variables
c
c      character*200 buffer
c      real xval,yval(maxcols),scalef
c      integer extstr,len,start,in,startn,endn,lenstr,colcnt
c      external extstr,lenstr
c      integer dataset_num,dataset_cnt,total_dataset_cnt
c
c     The following quantities were introduced to allow adjustments
c     to invlaid experimental values without actually changing the
c     contents of the file.
c
c     cutoff_val - this is for the elimination of spikes - any entries
c                  that are greater than this value will be set to zero.
c     offset_val - this is to graphically compensate for non-zero
c                  calibrartion offsets so plot comparison is easier.
c
c      real cutoff_val, offset_val,minexpt
c
c
c     Initialize
c
c      scalef = 1.0
c      startn = 1
c      endn = maxdata
c      offset_val = 0.0
c      cutoff_val = HI
c      minexpt = HI
c      axis_type = 1
c      num_expt = 0
c      dataset_cnt = 0
c      ncols = 1
c      datatitle = 'NO TITLE'
c
c      write (6,*) 'ISELD:',iseld
c
c      rewind (dataunit)
c
c     Scan through data unit looking for INDEX keyword
c     with the appropriate data tag number.
c
c 100  read(dataunit,'(a200)',end=500,err=500) buffer
c
c      len = lenstr(buffer)
c      write(6,*) 'BUFF:',buffer(1:len),':'
c
c     Ignore Empty Lines
c
c      if (buffer.eq.'') goto 100
c
c     Ignore comments
c
c      if (buffer(1:1).eq.'$') goto 100
c
c     Check number of datasets in file
c
c      if (buffer(1:10).eq.'FILECOUNT:') then
c         read (buffer(11:),*) total_dataset_cnt
c
c         write (6,*) 'DATASETS:',total_dataset_cnt
c
c         if (iseld.gt.total_dataset_cnt) then
c            write (6,*) 'REQUESTED EXPERIMENTAL DATA SET # ',
c     >                  iseld,' DOES NOT EXIST'
c            write (6,*) 'ONLY ', total_dataset_cnt,' DATA SETS ARE '//
c     >                       'SPECIFIED BY FILECOUNT:'
c            write (0,*) 'REQUESTED EXPERIMENTAL DATA SET # ',
c     >                  iseld,' DOES NOT EXIST'
c            write (0,*) 'ONLY ', total_dataset_cnt,' DATA SETS ARE '//
c     >                       'SPECIFIED BY FILECOUNT:'
c         endif
c      endif
c
c     Look for INDEX and count for datasets
c
c      if (buffer(1:6).eq.'INDEX:') then
c         dataset_cnt = dataset_cnt+1
c
c         read(buffer(7:),*) dataset_num
c
c         write (6,*) 'INDEX:',dataset_cnt,dataset_num,iseld
c
c        Indexing error - write out error message and continue
c
c         if (dataset_cnt.ne.dataset_num) then
c
c            write (6,*) 'ERROR IN EXPERIMENTAL DATA'//
c     >               ' FILE: DATASET # ',dataset_num, ' IS IN'//
c     >                  ' FILE POSITION ',DATASET_CNT
C
c            write (0,*) 'ERROR IN EXPERIMENTAL DATA'//
c     >               ' FILE: DATASET # ',dataset_num, ' IS IN'//
c     >                  ' FILE POSITION ',DATASET_CNT
C
c         endif
c
c         if (dataset_num.eq.iseld.or.dataset_cnt.eq.iseld) goto 300
c
c      endif
c
c     Loop back and read more data until exit
c
c      goto 100
c
c     Index or count for dataset found - continue processing
c
c 300  continue
c
c     Read in the rest of the HEADER block until data found - then
c     start processing data until EOF or next INDEX.
c
c 350  read(dataunit,'(a200)',end=400,err=400) buffer
c
c      len = lenstr(buffer)
c      write(6,*) 'DATA:',len,':',buffer(1:len),':'
c
c     Ignore Empty lines
c
c      if (buffer.eq.'') goto 350
c
c     Ignore comments
c
c      if (buffer(1:1).eq.'$') goto 350
c
c     Exit if Next INDEX is found
c
c      if (buffer(1:6).eq.'INDEX:') goto 400
c
c     Extract a Title if one is specified
c
c      if (buffer(1:6).eq.'TITLE:') then
c         len = extstr(buffer(7:),start)
c         datatitle = buffer(6+start:len)
c         write(6,*) 'TITLE:',datatitle,':'
c      endif
c
c     Extract Scaling Factor if specified
c
c      if (buffer(1:7).eq.'SCALEF:') read(buffer(8:),*) scalef
c
c     Extract Axis_type if given
c
c      if (buffer(1:5).eq.'AXIS:') read(buffer(6:),*) axis_type
c
c     Extract Number of colums of data if given - one assumed
c
c      if (buffer(1:6).eq.'NCOLS:') then
c
c         read(buffer(7:),*) ncols
c
c         write(6,*) 'NCOLS:',ncols
c
c         if (ncols.gt.maxcols) then
c
c           Issue error message and only load the first column of
c           data
c
c            write (6,*) 'REQUESTED EXPERIMENTAL DATA SET # ',
c     >                  iseld,' HAS MORE COLUMNS THAN MAX =',maxcols
c            write (0,*) 'REQUESTED EXPERIMENTAL DATA SET # ',
c     >                  iseld,' HAS MORE COLUMNS THAN MAX =',maxcols
c            ncols = 1
c
c         endif
c      endif
c
c     Extract Offset Value if specified
c
c      if (buffer(1:7).eq.'OFFSET:') read(buffer(8:),*) offset_val
c
c     Extract Cutoff Value if specified
c
c      if (buffer(1:7).eq.'CUTOFF:') read(buffer(8:),*) cutoff_val
c
c     Extract data limit counters if present
c
c      if (buffer(1:6).eq.'COUNT:') then
c         read(buffer(7:),*) startn,endn
c         write (6,*) 'COUNT:',startn,endn
c      endif
c
c     Other headers will be ignored - check for data and load it.
c
c     Data lines start with a blank
c
c      if (buffer(1:1).eq.' ') then
c
c         if(lenstr(buffer).gt.6) then
c
c
c            write (0,*) 'BUFFER:',buffer,':'
c
c            read(buffer,*) in,xval,
c     >                        (yval(colcnt),colcnt=1,ncols)
c
c            if (in.ge.startn.and.in.le.endn) then
c
c               expt_axis(in-startn+1) = xval
c
c               if (axis_type.eq.5) then
c
c                  expt_axis(in-startn+1) = expt_axis(in-startn+1)
c     >                                     /360.0 * 2.0 * PI
c               endif
c
c               do colcnt = 1,ncols
c                  expt_data(in-startn+1,colcnt) = yval(colcnt)*scalef
c               end do
c
c               num_expt = in-startn+1
c
c              Keep track of minimum value in first experimental data.
c
c               minexpt = min(minexpt,expt_data(in-startn+1,1))
c
c               write (6,*) 'NUM:',num_expt,in,xval,yval(1),scalef
c
c            endif
c
c         endif
c
c      endif
c
c      goto 350
c
c     Continue and wrap up file processing - experimental data
c     has been read.
c
c 400  continue
c
c     Data cleanup processing ...
c
c     Check experimental data for cutoff and adjust by offset.
c
c      if (offset_val.eq.-1.0) offset_val = minexpt
c
c      write (6,*) 'OFFSET_VAL:',offset_val
c
c      do in = 1,num_expt
c
c         do colcnt = 1,ncols
c
c            expt_data(in,colcnt) = expt_data(in,colcnt) - offset_val
c
c            if (expt_data(in,colcnt).gt.cutoff_val)
c     >                                  expt_data(in,colcnt) = 0.0
c
c         end do
c
c      end do
c
c      return
c
c 500  continue
c
c     Error exit condition - DATA SET NOT FOUND
c
c      write (6,*) 'ERROR IN EXPERIMENTAL DATA: EOF OR DATASET # ',
c     >             iseld,' NOT FOUND'
c      write (0,*) 'ERROR IN EXPERIMENTAL DATA: EOF OR DATASET # ',
c     >             iseld,' NOT FOUND'
c
c     Exit
c
c      return
c      end
c
c
c
      real function grid_interpolate(r,z,interp_opt,vs,iis,maxiis,
     >                               vst,vstsw)
      use mod_params
      use mod_cgeom
      implicit none
c     include 'params'
c     include 'cgeom'
      integer interp_opt,iis,maxiis,vstsw
      real r,z,vs(maxnks,maxnrs,maxiis),vst(maxnds)
c
c     GRID_INTERPOLATE: Using the interpolation option specified by
c                       interp_opt - this routine finds an estimate
c                       for the value of vs at the location specified
c                       by R,Z. It only returns usable data for locations
c                       actually on the grid and will return a 0.0
c                       value as an error condition.
c
c     Local variables
c
      integer ik,ir,in
      real p1(3),p2(3),p3(3),pval
      logical griderr
c
c     Set default return value of 0.0 for error conditions
c
      grid_interpolate = 0.0
c
c     First determine what cell (if any) the given R,Z location is in
c
      call gridpos(ik,ir,r,z,.false.,griderr)
c
      write(6,*) 'GRIDPOS:',ik,ir,r,z,griderr,korpg(ik,ir),
     >                      nvertp(korpg(ik,ir))
c slmod begin
      WRITE(0,*) ' STOP: GFORTRAN COMPLAINING, CHECK CODE' 
      STOP
c      write(6,'(8(1x,g12.5))') ((rvertp(in,korpg(ik,ir)),
c     >                         zvertp(in,korpg(ik,ir))),
c     >                         in = 1,4)
c slmod end
c



c
c     Particle was not found on the grid - exit routine
c
      if (griderr) return
c
c     If interpolation option 0 is selcetd return value for cell
c
      if (interp_opt.eq.0) then
c
         grid_interpolate=vs(ik,ir,iis)
c
c     Estimate value at actual point using a planar approximation based
c     on quadrant.
c
      elseif (interp_opt.eq.1) then
c
c        Assign values from cell centre to first vector
c
         p1(1) = rs(ik,ir)
         p1(2) = zs(ik,ir)
         p1(3) = vs(ik,ir,iis)
c
c        Determine which section of the grid cell in which the point lies
c        Use angles to middles of polygon sides.
c
c        Determine if the point lies to inside or outside of flux line.
c
         call cell_section(ik,ir,r,z,p2,p3,vs,iis,maxiis,vst,vstsw)
c
         call calc_plane_value(r,z,pval,p1,p2,p3)
c
         grid_interpolate=pval
c
      endif
c
      return
c
      end
c
c
c
      subroutine cell_section(ik,ir,r,z,p2,p3,
     >                        vs,iis,maxiis,vst,vstsw)
      use mod_params
      use mod_cgeom
      implicit none
c     include 'params'
c     include 'cgeom'
      integer ik,ir,vstsw,iis,maxiis
      real r,z,p2(3),p3(3),vs(maxnks,maxnrs,maxiis),vst(maxnds)
c     
c     This routine estimates the VS values at the middle of the polygon
c     sides and returns the relevant ones in the two vectors.
c
c     Local variables
c
      integer in,ind,nsides,ikn,irn,nvert
      real rm(4),zm(4),vm(4),rvert(4),zvert(4),fact
      external inpoly
      logical inpoly

c
c     Initialization
c
      ind = korpg(ik,ir)
c
      nsides = nvertp(ind)
c
c     Calculate the middle of the polygon sides
c
      do in = 1,nsides
c
         if (in.eq.nsides) then
c
            rm(in) = (rvertp(1,ind) + rvertp(in,ind))/2.0
            zm(in) = (zvertp(1,ind) + zvertp(in,ind))/2.0
c
         else
c
            rm(in) = (rvertp(in,ind) + rvertp(in+1,ind))/2.0
            zm(in) = (zvertp(in,ind) + zvertp(in+1,ind))/2.0
c
         endif
c
      end do
c
c     Calculate the values of the function on each side.
c
c     This code presently assumes a 4-sided polygon with two sides roughly
c     perpendicular and two sides roughly parallel to the field lines.
c
c     Deal with (1,2) side first
c
      if (ik.eq.1) then
c
         if (vstsw.eq.1) then
c
            vm(1) = vst(idds(ir,2))
c
         elseif (vstsw.eq.0) then
c
            vm(1) = vs(ik,ir,iis)
c
         endif
c
      else
c
         vm(1) = vs(ik-1,ir,iis)
     >        + ((ksb(ik-1,ir)-kss(ik-1,ir))
     >          /(kss(ik,ir)-kss(ik-1,ir))
     >          * (vs(ik,ir,iis)-vs(ik-1,ir,iis)))
c
      endif
c
c     Deal with 34 side
c
      if (ik.eq.nks(ir)) then
c
         if (vstsw.eq.1) then
c
            vm(3) = vst(idds(ir,1))
c
         elseif (vstsw.eq.0) then
c
            vm(3) = vs(ik,ir,iis)
c
         endif
c
      else
c
         vm(3) = vs(ik,ir,iis)
     >        + ((ksb(ik,ir)-kss(ik,ir))
     >          /(kss(ik+1,ir)-kss(ik,ir))
     >          * (vs(ik+1,ir,iis)-vs(ik,ir,iis)))
c
      endif
c
c     Deal with parallel sides - inward first
c     This is the 2,3 polygon side.
c
      ikn = ikins(ik,ir)
      irn = irins(ik,ir)
c
c     Check for cases of invalid inward rings - ir.eq.1 may need
c     to need to be changed to ir.eq.ircore when this value is
c     propagated to OUT.
c
      write(6,*) 'VM(2):',ik,ir,ikn,irn,irtrap,irtrap+1,
     >            vs(ik,ir,iis)
c
      if ((ikn.eq.ik.and.irn.eq.ir).or.
     >    (rs(ik,ir).eq.rs(ikn,irn).and.zs(ik,ir).eq.zs(ikn,irn))
     >    .or.ir.eq.irtrap+1.or.ir.eq.2) then
c
c        Assign value of current cell for these conditions.
c
         vm(2) = vs(ik,ir,iis)
c
      else
c
c        When distin/distout tdistin/tdistout are available in OUT -
c        consider replacing the following poloidal distance code
c        with these values.
c
c        Calculate value at cell side midpoint
c
         fact = ((rm(2)-rs(ik,ir))**2+(zm(2)-zs(ik,ir))**2)
     >         / ((rs(ik,ir)-rs(ikn,irn))**2
     >           +(zs(ik,ir)-zs(ikn,irn))**2)
c
         vm(2) = vs(ik,ir,iis) + fact*(vs(ikn,irn,iis)-vs(ik,ir,iis))
c
      endif
c
c     Deal with parallel sides - outward next
c     This is the 1,4 polygon side.
c
      ikn = ikouts(ik,ir)
      irn = irouts(ik,ir)
c
c     Check for cases of invalid inward rings - ir.eq.1 may need
c     to need to be changed to ir.eq.ircore when this value is
c     propagated to OUT.
c
      if ((ikn.eq.ik.and.irn.eq.ir).or.
     >    (rs(ik,ir).eq.rs(ikn,irn).and.zs(ik,ir).eq.zs(ikn,irn))
     >    .or.ir.eq.irwall-1) then
c
c        Assign value of current cell for these conditions.
c
         vm(4) = vs(ik,ir,iis)
c
      else
c
c        When distin/distout tdistin/tdistout are available in OUT -
c        consider replacing the following poloidal distance code
c        with these values.
c
c        Calculate value at cell side midpoint
c
         fact = ((rm(4)-rs(ik,ir))**2+(zm(4)-zs(ik,ir))**2)
     >         / ((rs(ik,ir)-rs(ikn,irn))**2
     >           +(zs(ik,ir)-zs(ikn,irn))**2)
c
         vm(4) = vs(ik,ir,iis) + fact*(vs(ikn,irn,iis)-vs(ik,ir,iis))
c
      endif
c
c     Now find which quadrant of the cell the test point lies in.
c
c
c     Do each quadrant in turn
c
      nvert = 4
      rvert(1) = rs(ik,ir)
      zvert(1) = zs(ik,ir)
c
      do in = 1,nsides
c
         if (in.eq.nsides) then
c
            rvert(2) = rm(in)
            zvert(2) = zm(in)
c
            rvert(3) = rvertp(1,ind)
            zvert(3) = zvertp(1,ind)
c
            rvert(4) = rm(1)
            zvert(4) = zm(1)
c
         else
c
            rvert(2) = rm(in)
            zvert(2) = zm(in)
c
            rvert(3) = rvertp(in+1,ind)
            zvert(3) = zvertp(in+1,ind)
c
            rvert(4) = rm(in+1)
            zvert(4) = zm(in+1)
c
         endif
c
         if (inpoly(r,z,nvert,rvert,zvert)) then
c
c            write (6,*) 'INPOLY:'
c
            if (in.eq.nsides) then

               p2(1) = rm(in)
               p2(2) = zm(in)
               p2(3) = vm(in)

               p3(1) = rm(1)
               p3(2) = zm(1)
               p3(3) = vm(1)

            else

               p2(1) = rm(in)
               p2(2) = zm(in)
               p2(3) = vm(in)

               p3(1) = rm(in+1)
               p3(2) = zm(in+1)
               p3(3) = vm(in+1)

            endif
c
c           Exit loop
c
            exit

         endif
c
      end do
c
c      write (6,*) 'CELL:',r,z,ik,ir,ind
c      do in = 1,4
c         write(6,'(a,6(1x,g12.5))') 'CELL DATA:',
c     >              rm(in),zm(in),vm(in)
c      end do
c
c
      return
c
      end
c
c
c
      subroutine calc_plane_value(r,z,pval,p1,p2,p3)
      implicit none
      real r,z,pval,p1(3),p2(3),p3(3)
c
c     Calculate Plane value: This subroutine uses the three
c     vectors given to define a plane then find the value on
c     this plane at the specified R,Z coordinates. Error
c     checking in this routine is minimal. It assumes that the
c     vectors supplied are not degenerate.
c
c     The first item is to calculate the normal to the plane
c     defined by the three points. Define the vectors from
c     P1 -> P2 and P1-> P3
c
      integer in
      real a(3),b(3),norm(3),d
c
c     Assign vector values
c
c      write (6,*) 'P1:',p1(1),p1(2),p1(3)
c      write (6,*) 'P2:',p2(1),p2(2),p2(3)
c      write (6,*) 'P3:',p3(1),p3(2),p3(3)
c
      do in = 1,3
c
         a(in) = p2(in) - p1(in)
         b(in) = p3(in) - p1(in)
c
      end do
c
c     Calculate normal
c
      norm(1) = a(2) * b(3) - a(3) * b(2)
      norm(2) = a(3) * b(1) - a(1) * b(3)
      norm(3) = a(1) * b(2) - a(2) * b(1)
c
c     Calculate constant
c
      d = norm(1)*p1(1)+norm(2)*p1(2)+norm(3)*p1(3)
c
c     Calculate value at given coordinates
c
      if (norm(3) .eq. 0.0) then
         write(6,*) 'ERROR:calc_plane_value: vertical plane:',
     >           norm(1),norm(2),norm(3),d
         pval = 0.0
      else
c
         pval = (d-norm(1)*r-norm(2)*z)/norm(3)
c
      endif
c
c      write (6,*) 'PVAL:',pval
c
c
c     Exit
c
      return
c
      end
c
c
c
      subroutine load_divdata_array(tmpplot,iselect,istate,itype,
     >     ylab,blab,ref,nizs,ierr)
      use error_handling
      use mod_params
      use mod_cgeom
      use mod_dynam2
      use mod_dynam3
      use mod_comtor
      use mod_pindata
      use mod_reiser_com
      use mod_slcom
      use mod_cedge2d
      use mod_adas_data_spec
      use mod_driftvel
      implicit none
c     include 'params' 
c     include 'cgeom'
c     include 'dynam2'
c     include 'dynam3'
c     include 'comtor'
c     include 'pindata'
c     
c     include 'reiser_com'
c     
c     include 'slcom'
c     include 'cedge2d'
c     include 'adas_data_spec'
c     include 'driftvel'
c     
      real tmpplot(maxnks,maxnrs)
      integer iselect,istate,nizs,ierr,itype      
      character*(*) ylab,blab,ref 
c     
c     LOAD_DIVDATA_ARRAY
c     
c     This routine loads a 2D DIVIMP array of size MAXNKS,MAXNRS with
c     a quantity specified by the values of iselect and istate. 
c     The allowed values of ISELECT are:
c     
c     ITYPE specifies the type of plot - 0 = contour, 1 = integrated
c     ITYPE may also be used to indicate plot specific options 
c     
c     ISELECT = 1 = TOTAL H POWER LOSS  (W)
c     2 = TOTAL IMPURITY POWER LOSS  (W)
c     3 = TOTAL POWER LOSS   (W)
c     4 = SPECIFIED IMPURITY SPECTROSCOPIC LINE 
c     - NEED TO READ ADAS DATA
c     5 = SPECIFIED HYDROGENIC SPECTROSCOPIC LINE 
c     - NEED TO READ ADAS DATA
c     6 = PIN Halpha from PINALPHA array
c     7 = PIN HALPHA - By Component from Eirene - 6 for total
c     - state specifies component
c     1 - H ionisation
c     2 - H+ recombination
c     3 - H2 dissociation
c     4 - H2+ dissociation
c     5 - CX of H and H+
c     6 - TOTAL 
c
c     8 = PIN HGAMMA - By component from Eirene - 6 for total
c     - as above 
c     9 = Hydrogen Neutral Density 
c     10 = Background Plasma Properties
c     1 = density
c     2 = electron temperature
c     3 = ion temperature
c     4 = velocity
c     5 = electric field
c     6 = mach number = velocity/cs (in cell)
c
c     11 = Impurity Species Density - specified by charge state
c     12 = Impurity Species Temperature - specified by charge state
c     13 = Impurity Species Velocity - specified by charge state
c     14 = TOTAL H POWER LOSS (W/m3)
c     15 = TOTAL IMPURITY POWER LOSS (W/m3)
c     16 = TOTAL POWER LOSS (W/m3)
c     17 = Load PLRP (Particular Line Radiation Profile - see PLRP 
c     module for istate values.
c     18 = Fluid code Background Plasma Properties
c     1 = density
c     2 = electron temperature
c     3 = ion temperature
c     4 = velocity
c     5 = electric field
c
c     19 = Fluid code Impurity Species Density - specified by charge state
c     20 = Fluid code Impurity Species Temperature - specified by charge state
c     21 = Fluid code Impurity Species Velocity - specified by charge state
c     22 = SPECIFIED IMPURITY SPECTROSCOPIC LINE AVERAGED TEMPERATURE
c     - MAY NEED TO READ ADAS DATA
c     23 = Impurity Density to Background Ne Ratio
c     Istate = IZ
c     24 = Impurity Temperature to Background Te Ratio
c     Istate = IZ
c     25 = Impurity Velocity to Background Vb Ratio
c     Istate = IZ
c     26 = PIN HBETA - By Component from Eirene - 6 for total
c     - state specifies component
c     1 - H ionisation
c     2 - H+ recombination
c     3 - H2 dissociation
c     4 - H2+ dissociation
c     5 - CX of H and H+
c     6 - TOTAL 
c
c     27 = BRATIO - magnetic field ratios or angles 
c     1 - Ratio of Bpol/Btor 
c     2 - Angle of Btot from "surface" (deg) asin(BRATIO) *180/PI
c
c     28 = HC - Calculation of CD EMISSION (D/XB)
c     istate = specific value 
c     1 - CD Efficiency (D/XB)
c     2 - CD Emissivity (photons/m3)
c
c     29 = HC - HC State density
c     istate = specific HC species 
c     = sum over states for greater than maxstate   
c     1 = C+ (from HC module)
c     2 = C  (from HC module)
c     3 = CH+(from HC module)
c     4 = CH (from HC module)
c
c     30 = HC - HC State Ionization
c     istate = specific HC species (ONLY CH So far)
c     31 = Impurity Ionization - specified by source charge state
c     
c     NOTE: Subgrid Iselect values are loaded by the load_subgrid_array routine found in the
c     subgrid_plots module - these are only listed here for completeness - local code in 
c     the plotting routines has to invoke the appropriate load routine since they require
c     different types of storage.
c     
c     32 = Subgrid impurity density - STATE = IZ
c     33 = Subgrid HC density - STATE = HC STATE INDEX
c     34 = Subgrid impurity ADAS based emissions - additional data read 
c     35 = Subgrid CH emission
c     
c     **** NOTE: When adding new options - increase the value of parameter max_iselect below *****
c     
c     36 = PIN Data 
c     1 = PINION = PIN ionization    
c     2 = PINATOM = PIN Atom density 
c     3 = PINMOL = PIN Molecular density
c     4 = PINIONZ = Impurity ionization
c     5 = PINZ0 = Impurity neutral density  
c     6 = PINQI = Ion heating term
c     7 = PINQE = Electron heating term
c     
c     37 = POWER LOSS (EXCITATION)
c     38 = POWER LOSS (RECOMBINATION/BREM)
c     39 = POWER LOSS (TOTAL) 
c
c     40 = EXB drift related quantities
c        1 = potential
c        2 = Radial E-field
c        3 = Poloidal E-field
c        4 = ExB Radial drift 
c        5 = ExB poloidal Drift 
c        6 = ExB Radial flux (ne x Vexb)
c        7 = ExB Poloidal flux (ne x Vexb)
c     41 = Impurity Emission - Tungsten WI only for now - using defined SXB function
c
c     42 = Power Balance components
c        1 = Ion conduction
c        2 = Ion convection
c        3 = Electron conduction
c        4 = Electron convection
c        5 = Total Conduction
c        6 = Total Convection
c        7 = Total Convection/total conduction
c        8 = Total Convection/electron conduction
c        
c     
c     
      integer max_iselect
      parameter (max_iselect=42)
c     
c     
c     ADAS variables
c     
      CHARACTER ADASID*80,graph3*80
      CHARACTER XFESYM*2
      character adasex*3
      integer   adasyr
      INTEGER ISELE,ISELR,ISELX,iseld,ircode
      integer line 
      REAL WLNGTH
c     
c     Local variables
c     
      real pltmax, pltmin
      integer ik1,ir1 
c     
      real zero_fact
c     
      real mfact,fact
      integer ik,ir,iz,len,lenstr
      external lenstr
c
      real,external :: wi_sxb
      real,external :: power_term
c     
      real cs
c     
c     Calculating radiative power
c     
      real :: ptesa(maxnks),pnesa(maxnks)
      real :: pcoef4(maxnks),pcoef5(maxnks),pnhs(maxnks)
      real :: pnbs(maxnks)
      character*2 :: year
      integer :: iclass
c     
c     Check for subrid ISELECT values which should not be passed
c     to this routine!
c     
      if (iselect.eq.32.or.iselect.eq.33.or.iselect.eq.34.or.
     >     iselect.eq.35) then
         call errmsg('LOAD_DIVDATA_ARRAY:'//
     >        'SUBGRID ISELECT VALUES SPECIFIED AS INPUT',iselect)
         ierr = 1
         return
      endif

c     
c     Check for valid ISELECT as input
c     
      if (iselect.lt.1.or.iselect.gt.max_iselect) then 
         call errmsg('LOAD_DIVDATA:DATA SELECTOR IN'//
     >        ' PLOT FILE IS OUT OF RANGE:',iselect)
         ierr = 1
         return
      endif

c     
c     Echo input
c     
      write(6,'(a,3i5)') 'Loading DIVDATA:',iselect,istate,ierr
      ierr = 0

c     
c     Set the YLAB value
c     
      call set_ylab(iselect,istate,itype,nizs,ylab) 
c     
c     Set the BLAB value
c     
      call set_blab(iselect,istate,itype,nizs,blab) 
c     
c     Initialize scaling factor
c     
      mfact = 1.0
c     
c----------------------------------------------------------
c     
c     Hydrogenic power loss - DIVIMP
c     
c----------------------------------------------------------
c     
      if (iselect.eq.1) then
c     
c     Individual states
c     
c     BLAB = 'BOLO H POWER LOSS (BOLO)'
c     
         if (istate.eq.0.or.istate.eq.1) then 
c     
            do ir = 1,nrs
c     
               do ik = 1, nks(ir)
c     
                  tmpplot(ik,ir) = tmpplot(ik,ir) 
     >                 + hpowls(ik,ir,istate)
     >                 * kareas(ik,ir)
c     
               end do
c     
            end do   
c     
c     Total Hydrogenic 
c     
         else
c     
            do ir = 1,nrs
c     
               do ik = 1, nks(ir)
c     
                  do iz = 0,1
c     
                     tmpplot(ik,ir) = tmpplot(ik,ir) 
     >                    + hpowls(ik,ir,iz)
     >                    * kareas(ik,ir)
c     
                  end do
c     
               end do
c     
            end do   

         endif
c     
c----------------------------------------------------------
c     
c     Impurity power loss - DIVIMP
c     
c----------------------------------------------------------
c     
      elseif (iselect.eq.2) then
c     
c     Individual charge state 
c     
c     
c     Scale by MFACT if required
c     
         IF (ABSFAC.GT.0.0) MFACT = MFACT * ABSFAC
c     
c     BLAB = 'BOLO IMP POW LOSS'
c     
         if (istate.ge.0.and.istate.le.nizs) then 
c     
            do ir = 1,nrs
c     
               do ik = 1, nks(ir)
c     
                  tmpplot(ik,ir) = tmpplot(ik,ir) 
     >                 + powls(ik,ir,istate)*mfact
     >                 * kareas(ik,ir)
c     
               end do
c     
            end do   
c     
c     Total Impurity
c     
         else 
c     
            do ir = 1,nrs
c     
               do ik = 1, nks(ir)
c     
                  do iz = 0,nizs
c     
                     tmpplot(ik,ir) = tmpplot(ik,ir) 
     >                    + powls(ik,ir,iz)*mfact
     >                    * kareas(ik,ir)
c     
                  end do
c     
               end do
c     
            end do   


         endif
c     
c----------------------------------------------------------
c     
c     Total power loss - Hydrogenic + Impurity - DIVIMP 
c     
c----------------------------------------------------------
c     
      elseif (iselect.eq.3) then
c     
c     BLAB = 'BOLO TOTAL POW LOSS'
c     
c     Scale by MFACT if required
c     
         IF (ABSFAC.GT.0.0) MFACT = MFACT * ABSFAC
c     
c     Hydrogenic
c     
c     
         do ir = 1,nrs
c     
            do ik = 1, nks(ir)
c     
               do iz = 0,1
c     
                  tmpplot(ik,ir) = tmpplot(ik,ir) 
     >                 + hpowls(ik,ir,iz)
     >                 * kareas(ik,ir)
c     
               end do
c     
            end do
c     
         end do   
c     
c     Impurity 
c     
         do ir = 1,nrs
c     
            do ik = 1, nks(ir)
c     
               do iz = 0,nizs
c     
                  tmpplot(ik,ir) = tmpplot(ik,ir) 
     >                 + powls(ik,ir,iz)*mfact
     >                 * kareas(ik,ir)
c     
               end do
c     
            end do
c     
         end do   
c     
c----------------------------------------------------------
c     
c     ADAS based - Impurity spectral line
c     
c----------------------------------------------------------
c     
      elseif (iselect.eq.4.or.iselect.eq.22) then
c     
c     BLAB = 'CODE ADAS IMP PLRP'
c     
c     Need to read in ADAS data spec to calculate radiation
c     

         if (cadas_switch.eq.0) then  
            
            CALL RDG1 (GRAPH3,ADASID,adasyr,adasex,
     >           ISELE,ISELR,ISELX,ISELD,IERR)
c     
c     Save the ADAS data read in into the common block for 
c     possible re-use. Do not set the cadas_switch.
c     
            cadasid = adasid
            cadasyr = adasyr
            cadasex = adasex
            cisele  = isele
            ciselr  = iselr
            ciselx  = iselx
            ciseld  = iseld
c     
c     Use ADAS data in common instead of reading from input 
c     
         elseif (cadas_switch.eq.1) then 
c     
            adasid = cadasid
            adasyr = cadasyr
            adasex = cadasex
            isele  = cisele
            iselr  = ciselr
            iselx  = ciselx
            iseld  = ciseld
c     
         endif
c     
         write(6,'(a,i5,a,i5,a,5i5)') 'LOAD_DIVDATA: ADAS:',
     >        cadas_switch,trim(adasid),adasyr,
     >        trim(adasex),isele,iselr,iselx,iseld

c     
         if (ierr.ne.0) return 
c     
         IF (ISTATE.GE.0.AND.ISTATE.LE.NIZS.AND.ISTATE.LT.CION)THEN
c     
            call LDADAS(CION,ISTATE,ADASID,ADASYR,ADASEX,
     >           ISELE,ISELR,ISELX,
     >           tmpplot,Wlngth,IRCODE)
c     
            IF (IRCODE.NE.0) THEN
               WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
               return   
            ENDIF
c     
            if (iselect.eq.4) then   
               REF = 'ADAS PLRP XX XXXXX ('
            elseif (iselect.eq.22) then 
               REF = 'ADAS TEMP XX XXXXX ('
            endif             
c     
            WRITE(REF(11:12),'(I2)') ISTATE
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
c     
         endif 
c     
c     Scale by MFACT
c     
c     write(6,*) 'LOAD:',mfact,absfac
c     
         IF (ABSFAC.GT.0.0) MFACT = MFACT * ABSFAC
c     
         do ir = 1,nrs
            do ik = 1,nks(ir)
               tmpplot(ik,ir) = tmpplot(ik,ir) * mfact
            end do 
         end do
c     
c     Scale by the impurity temperature for 
c     iselect option 22
c     
         if (iselect.eq.22) then 
            do ir = 1,nrs
               do ik = 1,nks(ir)
                  tmpplot(ik,ir) = tmpplot(ik,ir) 
     >                 * sdts(ik,ir,istate)
c     
               end do 
            end do
         endif
c     
c----------------------------------------------------------
c     
c     ADAS based - Hydrogenic spectral lines
c     
c----------------------------------------------------------
c     
      elseif (iselect.eq.5) then
c     
c     BLAB = 'CODE ADAS H PLRP'
c     
c     Need to read in ADAS data spec to calculate radiation
c     
c     
         if (cadas_switch.eq.0) then  
            
            CALL RDG1 (GRAPH3,ADASID,adasyr,adasex,
     >           ISELE,ISELR,ISELX,ISELD,IERR)
c     
c     Use ADAS data in common instead of reading from input 
c     
         elseif (cadas_switch.eq.1) then 
c     
            adasid = cadasid
            adasyr = cadasyr
            adasex = cadasex
            isele  = cisele
            iselr  = ciselr
            iselx  = ciselx
            iseld  = ciseld
c     
         endif
c     
         if (ierr.ne.0) return 
c     
         IF (ISTATE.GE.0 .AND. ISTATE.LE.1) THEN
c     
            call LDADAS(1,ISTATE,ADASID,ADASYR,ADASEX,
     >           ISELE,ISELR,ISELX,
     >           tmpplot,Wlngth,IRCODE)
c     
            IF (IRCODE.NE.0) THEN
               WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
               return   
            ENDIF
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
c     
         endif 
c     
c     Scale by MFACT
c     
c     IF (ABSFAC.GT.0.0) MFACT = MFACT * ABSFAC
c     
         do ir = 1,nrs
            do ik = 1,nks(ir)
               tmpplot(ik,ir) = tmpplot(ik,ir) * mfact

               write(6,'(a,2i6,15(1x,g12.5))') 'HADAS:',ik,ir,mfact,
     >           tmpplot(ik,ir),ktebs(ik,ir),knbs(ik,ir),pinatom(ik,ir)
            end do 
         end do




c     
c----------------------------------------------------------
c     
c     PIN Halpha and Hgamma
c     
c----------------------------------------------------------
c     
      elseif (iselect.eq.6) then
c     
c     PIN Halpha - Total only 
c     
c     BLAB = 'CODE CODE HALPHA'
c     
c     Loop through array
c     
         do ir = 1,nrs
c     
            do ik = 1, nks(ir)
c     
               tmpplot(ik,ir) = pinalpha(ik,ir)
c     
            end do
c     
         end do   
c     
c     
c----------------------------------------------------------
c     
c     PIN Halpha and Hgamma - By component from Eirene 
c     
c----------------------------------------------------------
c     
      elseif (iselect.eq.7.or.iselect.eq.8.or.iselect.eq.26) then
c     
c     PIN Halpha  
c     
         if (iselect.eq.7) then 
c     
c     BLAB = 'CODE CODE HALPHA'
            line = H_BALPHA
c     
         elseif (iselect.eq.8) then 
c     
c     BLAB = 'CODE CODE HGAMMA'
            line = H_BGAMMA
c     
         elseif (iselect.eq.26) then 
c     
c     BLAB = 'CODE CODE HBETA'
            line = H_BBETA
c     
         endif
c     
c     Loop through array
c     
         do ir = 1,nrs
c     
            do ik = 1, nks(ir)
c     
               tmpplot(ik,ir) = pinline(ik,ir,istate,line)
c     
            end do
c     
         end do   
c     
c----------------------------------------------------------
c     
c     PIN Hneutral Density - from Eirene
c     
c----------------------------------------------------------
c     
      elseif (iselect.eq.9) then  
c     
         do ir = 1,nrs
c     
            do ik = 1, nks(ir)
c     
               tmpplot(ik,ir) = pinatom(ik,ir)
c     
            end do
c     
         end do   
c     
c----------------------------------------------------------
c     
c     DIVIMP Background Plasma Properties - Ne, Te, Ti, Vb, E
c     
c----------------------------------------------------------
c     
      elseif (iselect.eq.10) then  
c     
         do ir = 1,nrs
c     
            do ik = 1, nks(ir)
c     
               if (istate.eq.1) then 
                  tmpplot(ik,ir) = knbs(ik,ir)
               elseif (istate.eq.2) then 
                  tmpplot(ik,ir) = ktebs(ik,ir)
               elseif (istate.eq.3) then 
                  tmpplot(ik,ir) = ktibs(ik,ir)
               elseif (istate.eq.4) then 
                  tmpplot(ik,ir) = kvhs(ik,ir) / qtim
               elseif (istate.eq.5) then 
                  tmpplot(ik,ir) = kes(ik,ir)
               elseif (istate.eq.6) then 
                  ! Mach number (signed? yes for now)
                  CS = 9.79E3 * SQRT (0.5*(KTEBS(Ik,IR)+KTIBS(ik,IR))*
     >                (1.0+RIZB)/CRMB)
                  if (cs.ne.0.0) then 
                     tmpplot(ik,ir) = kvhs(ik,ir)/qtim/cs
                  endif
c                  write(0,'(a,2i6,10(1x,g12.5))') 'Mach:',ik,ir,qtim,
c     >                                  cs,kvhs(ik,ir)/qtim,
                  
               endif
c     
            end do
c     
         end do   

c     
c----------------------------------------------------------
c     
c     DIVIMP Impurity Species Densities
c     AND Density Ratio to Background Plasma
c     
c----------------------------------------------------------
c     
      elseif (iselect.eq.11.or.iselect.eq.23) then  
c     
c     Scaling factor 
c     
         IF (ABSFAC.GT.0.0) MFACT = MFACT * ABSFAC
c     
         do ir = 1,nrs
c     
            do ik = 1, nks(ir)
c     
               if (istate.eq.nizs+1) then 

                  do iz = 0,nizs
                     tmpplot(ik,ir) = tmpplot(ik,ir) + 
     >                    sdlims(ik,ir,iz) * mfact
                  end do
               else
                  tmpplot(ik,ir) = sdlims(ik,ir,istate)*mfact
               endif
c     
               if (iselect.eq.23) then 
c     
                  if (knbs(ik,ir).ne.0.0) then 

                     tmpplot(ik,ir) = tmpplot(ik,ir) 
     >                    / knbs(ik,ir) 

                  else

                     tmpplot(ik,ir) = 0.0

                  endif

               endif  

c     
            end do
c     
         end do   

c     
c----------------------------------------------------------
c     
c     DIVIMP Impurity Species Temperatures
c     AND Temperature Ratio to Background Plasma
c     
c----------------------------------------------------------
c     

      elseif (iselect.eq.12.or.iselect.eq.24) then  
c     
         do ir = 1,nrs
c     
            do ik = 1, nks(ir)
c     
               tmpplot(ik,ir) = sdts(ik,ir,istate)
c     
               if (iselect.eq.24) then 

                  if (ktebs(ik,ir).ne.0.0) then 

                     tmpplot(ik,ir) = tmpplot(ik,ir)
     >                    / ktebs(ik,ir) 

                  else

                     tmpplot(ik,ir) = 0.0

                  endif 

               endif 
c     
            end do
c     
         end do   
c     
c     
c----------------------------------------------------------
c     
c     DIVIMP Impurity Species Velocities
c     AND Velocity Ratio to Background Plasma
c     
c----------------------------------------------------------
c     
      elseif (iselect.eq.13.or.iselect.eq.25) then  
c     
         do ir = 1,nrs
c     
            do ik = 1, nks(ir)
c     
               tmpplot(ik,ir) = velavg(ik,ir,istate)
c     
               if (iselect.eq.25) then 
c     
                  if (kvhs(ik,ir).ne.0.0) then 
                     
                     tmpplot(ik,ir) = tmpplot(ik,ir)
     >                    /(kvhs(ik,ir)/qtim)

                  else
                     
                     tmpplot(ik,ir) = 0.0

                  endif   

               endif 
c     
            end do
c     
         end do   


c     
c----------------------------------------------------------
c     
c     Hydrogenic power loss - DIVIMP (W/m3)
c     
c----------------------------------------------------------
c     
      elseif (iselect.eq.14) then
c     
c     Individual states
c     
         if (istate.eq.0.or.istate.eq.1) then 
c     
            do ir = 1,nrs
c     
               do ik = 1, nks(ir)
c     
                  tmpplot(ik,ir) = tmpplot(ik,ir) 
     >                 + hpowls(ik,ir,istate)
c     
               end do
c     
            end do   
c     
c     Total Hydrogenic 
c     
         else
c     
            do ir = 1,nrs
c     
               do ik = 1, nks(ir)
c     
                  do iz = 0,1
c     
                     tmpplot(ik,ir) = tmpplot(ik,ir) 
     >                    + hpowls(ik,ir,iz)
c     
                  end do
c     
               end do
c     
            end do   

         endif
c     
c----------------------------------------------------------
c     
c     Impurity power loss - DIVIMP (W/m3)
c     
c----------------------------------------------------------
c     
      elseif (iselect.eq.15) then
c     
c     Individual charge state 
c     
c     
c     Scale by MFACT if required
c     
         IF (ABSFAC.GT.0.0) MFACT = MFACT * ABSFAC
c
c     
c     BLAB = 'BOLO IMP POW LOSS'
c     
         if (istate.ge.0.and.istate.le.nizs) then 
c     
            do ir = 1,nrs
c     
               do ik = 1, nks(ir)
c     
                  tmpplot(ik,ir) = tmpplot(ik,ir) 
     >                 + powls(ik,ir,istate)*mfact
c     
               end do
c     
            end do   
c     
c     Total Impurity
c     
         else 
c     
            do ir = 1,nrs
c     
               do ik = 1, nks(ir)
c     
                  do iz = 0,nizs
c     
                     tmpplot(ik,ir) = tmpplot(ik,ir) 
     >                    + powls(ik,ir,iz)*mfact
c     
                  end do
c     

               end do
c     
            end do   


         endif

c     
c----------------------------------------------------------
c     
c     Total power loss - Hydrogenic + Impurity - DIVIMP 
c     (W/m3) 
c     
c----------------------------------------------------------
c     
      elseif (iselect.eq.16) then
c     
c     BLAB = 'BOLO TOTAL POW LOSS'
c     
c     Scale by MFACT if required
c     
         IF (ABSFAC.GT.0.0) MFACT = MFACT * ABSFAC
c     
c     Hydrogenic
c     
c     
         do ir = 1,nrs
c     
            do ik = 1, nks(ir)
c     
               do iz = 0,1
c     
                  tmpplot(ik,ir) = tmpplot(ik,ir) 
     >                 + hpowls(ik,ir,iz)
c     
               end do
c     
            end do
c     
         end do   
c     
c     Impurity 
c     
         do ir = 1,nrs
c     
            do ik = 1, nks(ir)
c     
               do iz = 0,nizs
c     
                  tmpplot(ik,ir) = tmpplot(ik,ir) 
     >                 + powls(ik,ir,iz)*mfact
c     
               end do
c     
            end do
c     
         end do   

c     
c----------------------------------------------------------
c     
c     PLRP Calculated by the PLRP.o6a module
c     
c----------------------------------------------------------
c     
      elseif (iselect.eq.17) then  
c     
         do ir = 1,nrs
c     
            do ik = 1, nks(ir)
c     
               tmpplot(ik,ir) = plrps(ik,ir,istate)
c     
            end do
c     
         end do   
c     
c     
c----------------------------------------------------------
c     
c     FLUID CODE Solution Background Plasma Properties - Ne, Te, Ti, Vb, E
c     
c----------------------------------------------------------
c     
      elseif (iselect.eq.18) then  
c     
         do ir = 1,nrs
c     
            do ik = 1, nks(ir)
c     
               if (istate.eq.1) then 
                  tmpplot(ik,ir) = e2dnbs(ik,ir)
               elseif (istate.eq.2) then 
                  tmpplot(ik,ir) = e2dtebs(ik,ir)
               elseif (istate.eq.3) then 
                  tmpplot(ik,ir) = e2dtibs(ik,ir)
               elseif (istate.eq.4) then 
                  tmpplot(ik,ir) = e2dvhs(ik,ir)
               elseif (istate.eq.5) then 
                  tmpplot(ik,ir) = e2des(ik,ir)
               endif
c     
            end do
c     
         end do   

c     
c----------------------------------------------------------
c     
c     FLUID CODE Impurity Species Densities
c     
c----------------------------------------------------------
c     


      elseif (iselect.eq.19) then  
c     
         do ir = 1,nrs
c     
            do ik = 1, nks(ir)
c     
               tmpplot(ik,ir) = e2dnzs(ik,ir,istate)
c     
            end do
c     
         end do   

c     
c----------------------------------------------------------
c     
c     FLUID CODE Impurity Species Temperatures
c     
c----------------------------------------------------------
c     

      elseif (iselect.eq.20) then  
c     
c     Note: This assumes that impurity temperature is 
c     equal to the ion temperature for a fluid 
c     code. 
c     
         do ir = 1,nrs
c     
            do ik = 1, nks(ir)
c     
               tmpplot(ik,ir) = e2dtebs(ik,ir)
c     
            end do
c     
         end do   
c     
c     
c----------------------------------------------------------
c     
c     FLUID CODE Impurity Species Velocities
c     
c----------------------------------------------------------
c     


      elseif (iselect.eq.21) then  
c     
         do ir = 1,nrs
c     
            do ik = 1, nks(ir)
c     
               tmpplot(ik,ir) = e2dvzs(ik,ir,istate)
c     
            end do
c     
         end do   
c     
c     
c     
c----------------------------------------------------------
c     
c     Load HC data related quantities
c     
c----------------------------------------------------------
c     
      elseif (iselect.eq.28.or.iselect.eq.29.or.iselect.eq.30) then 

c     
         call load_hc_data_array(tmpplot,iselect,istate,itype,
     >        ref,nizs,mfact,absfac,ierr)
c     
c     Ionization data TIZS
c     
      elseif (iselect.eq.31) then 
c     
c     Scaling factor 
c     
c
         IF (ABSFAC.GT.0.0) MFACT = MFACT * ABSFAC
c     
         do ir = 1,nrs
c     
            do ik = 1, nks(ir)
c     
               if (istate.eq.nizs+1) then 

                  do iz = 0,nizs
                     tmpplot(ik,ir) = tmpplot(ik,ir) + 
     >                    tizs(ik,ir,iz) * mfact
                  end do
               else
                  tmpplot(ik,ir) = tizs(ik,ir,istate)*mfact
               endif
c     
            end do
c     
         end do   

      elseif (iselect.eq.36) then 
c     
c     Quantities returned from the PIN run and loaded into DIVIMP arrays
c     1 = PINION = PIN ionization    
c     2 = PINATOM = PIN Atom density 
c     3 = PINMOL = PIN Molecular density
c     4 = PINIONZ = Impurity ionization
c     5 = PINZ0 = Impurity neutral density  
c     6 = PINQI = Ion heating term
c     7 = PINQE = Electron heating term
c     
         if (istate.eq.1) then 
c     
c     PINION
c     
            do ir = 1,nrs
c     
               do ik = 1, nks(ir)
c     
                  tmpplot(ik,ir) = pinion(ik,ir)
c     
               end do
c     
            end do   
c     

         elseif (istate.eq.2) then 
c     
c     PINATOM
c     

            do ir = 1,nrs
c     
               do ik = 1, nks(ir)
c     
                  tmpplot(ik,ir) = pinatom(ik,ir)
c     
               end do
c     
            end do   

         elseif (istate.eq.3) then 
c     
c     PINMOL
c     
c     
            do ir = 1,nrs
c     
               do ik = 1, nks(ir)
c     
                  tmpplot(ik,ir) = pinmol(ik,ir)
c     
               end do
c     
            end do   
         elseif (istate.eq.4) then 
c     
c     PINIONZ
c     
            do ir = 1,nrs
c     
               do ik = 1, nks(ir)
c     
                  tmpplot(ik,ir) = pinionz(ik,ir)
c     
               end do
c     
            end do   
         elseif (istate.eq.5) then 
c     
c     PINZ0
c     
            do ir = 1,nrs
c     
               do ik = 1, nks(ir)
c     
                  tmpplot(ik,ir) = pinz0(ik,ir)
c     
               end do
c     
            end do   
         elseif (istate.eq.6) then 
c     
c     PINQI
c     
            do ir = 1,nrs
c     
               do ik = 1, nks(ir)
c     
                  tmpplot(ik,ir) = pinqi(ik,ir)
c     
               end do
c     
            end do   
         elseif (istate.eq.7) then 
c     
c     PINQE
c     
            do ir = 1,nrs
c     
               do ik = 1, nks(ir)
c     
                  tmpplot(ik,ir) = pinqe(ik,ir)
c     
               end do
c     
            end do   


         endif


      elseif (iselect.eq.37) then 
c     
c     Radiated power from hydrogen ... excitation ... recombination/brem or total
c     


         if (cadas_switch.eq.0) then  
            
            CALL RDG1 (GRAPH3,ADASID,adasyr,adasex,
     >           ISELE,ISELR,ISELX,ISELD,IERR)
c     
c     Save the ADAS data read in into the common block for 
c     possible re-use. Do not set the cadas_switch.
c     
            cadasid = adasid
            cadasyr = adasyr
            cadasex = adasex
            cisele  = isele
            ciselr  = iselr
            ciselx  = iselx
            ciseld  = iseld
c     
c     Use ADAS data in common instead of reading from input 
c     
         elseif (cadas_switch.eq.1) then 
c     
            adasid = cadasid
            adasyr = cadasyr
            adasex = cadasex
            isele  = cisele
            iselr  = ciselr
            iselx  = ciselx
            iseld  = ciseld
c     
         endif
c     
         write(6,'(a,i5,a,i5,a,5i5)') 'LOAD_DIVDATA: ADAS:',
     >        cadas_switch,trim(adasid),adasyr,
     >        trim(adasex),isele,iselr,iselx,iseld

         write(year,'(i2.2)') adasyr

         call xxuid(adasid)

         tmpplot = 0.0

         do ir = 1,nrs

            pcoef5 = 0.0
            pcoef4 = 0.0
            
            do ik = 1,nks(ir)
               PTESA(IK) = KTEBS(IK,IR)
               PNESA(IK) = KNBS(IK,IR) * RIZB
               PNBS(IK)  = KNBS(IK,IR)
               PNHS(IK)  = PINATOM(IK,IR)
            end do

!     iclass = 5 is plt
            ICLASS = 5
            CALL ADASRD(YEAR,1,1,ICLASS,NKS(IR),PTESA,PNESA,PCOEF5)


!     iclass = 4 is prb 
            ICLASS = 4
            CALL ADASRD(YEAR,1,1,ICLASS,NKS(IR),PTESA,PNESA,PCOEF4)


            if (istate.eq.0.or.istate.lt.0) then 

               DO IK = 1, NKS(IR)

                  tmpplot(ik,ir) = tmpplot(ik,ir) + 
     >                 PCOEF5(IK)*PNESA(IK)*PNHS(IK)
                  
               end do 

            endif


            if (istate.eq.1.or.istate.lt.0) then 

               DO IK = 1, NKS(IR)
                  tmpplot(ik,ir) = tmpplot(ik,ir) + 
     >                 PCOEF4(IK)*PNESA(IK)*PNBS(IK)

               end do   

            endif


            do ik = 1,nks(ir)
            write(6,'(a,2i6,15(1x,g12.5))') 'BOLO:',ik,ir,
     >           tmpplot(ik,ir),pcoef4(ik),pnesa(ik),pnbs(ik),
     >           pcoef4(ik)*pnesa(ik)*pnbs(ik),
     >           pcoef5(ik),pnesa(ik),pnhs(ik),
     >           pcoef5(ik)*pnesa(ik)*pnhs(ik)
            end do

         end do 

      elseif (iselect.eq.40) then 
c
c         ExB drift related quantities
c         1 - Potential (phi) (V) 
c         2 - Radial Efield (V/m)
c         3 - Poloidal Efield (V/m)
c         4 - Radial ExB drift (m/s)
c         5 - Poloidal ExB drift (m/s)
c         6 - Radial ExB flux   ne x Vexb_rad (/m2/s)
c         7 - Poloidal ExB flux ne x Vexb_pol (/m2/s)
c     
c     
         if (istate.eq.1) then 
c           Potential
            do ir = 1,nrs
               do ik = 1,nks(ir)
                  tmpplot(ik,ir) = osmpot2(ik,ir)
               end do
            end do
         elseif (istate.eq.2) then 
c           Radial electric field
            do ir = 1,nrs
               do ik = 1,nks(ir)
                  tmpplot(ik,ir) = e_rad(ik,ir)
               end do
            end do
         elseif (istate.eq.3) then 
c           Poloidal electric field
            do ir = 1,nrs
               do ik = 1,nks(ir)
                  tmpplot(ik,ir) = e_pol(ik,ir)
               end do
            end do
         elseif (istate.eq.4) then 
c           Exb Radial drift converted back to m/s by dividing by qtim
            do ir = 1,nrs
               do ik = 1,nks(ir)
                  tmpplot(ik,ir) = exb_rad_drft(ik,ir) / qtim
               end do
            end do
         elseif (istate.eq.5) then 
c           Exb Poloidal drift converted back to m/s by dividing by qtim
c           and taking out the geometric factor used to map to 2D
            do ir = 1,nrs
               do ik = 1,nks(ir)
                  fact = sqrt(kbfs(ik,ir)**2-1.0)
                  if (fact.ne.0.0) then 
                     tmpplot(ik,ir) = exb_pol_drft(ik,ir) / qtim /fact
                  else
                     tmpplot(ik,ir) = 0.0
                  endif
               end do
            end do
         elseif (istate.eq.6) then 
c           radial flux
c           Exb Radial drift converted back to m/s by dividing by qtim
            do ir = 1,nrs
               do ik = 1,nks(ir)
                  tmpplot(ik,ir) = knbs(ik,ir)*exb_rad_drft(ik,ir)/qtim
               end do
            end do
         elseif (istate.eq.7) then 
c           poloidal flux
c           Exb Poloidal drift converted back to m/s by dividing by qtim
c           and taking out the geometric factor used to map to 2D
            do ir = 1,nrs
               do ik = 1,nks(ir)
                  fact = sqrt(kbfs(ik,ir)**2-1.0)
                  if (fact.ne.0.0) then 
                     tmpplot(ik,ir) = knbs(ik,ir) * exb_pol_drft(ik,ir) 
     >                                     / qtim /fact
                  else
                     tmpplot(ik,ir) = 0.0
                  endif
               end do
            end do
         endif

c
c     End of ISELECT IF
c
c     
c     Tungsten emission data based on TIZS/SXB
c     
      elseif (iselect.eq.41) then 
c     
c     Tungsten emission using fixed sxb function
c     
c
         IF (ABSFAC.GT.0.0) MFACT = MFACT * ABSFAC
c     
         do ir = 1,nrs
c     
            do ik = 1, nks(ir)
c     
               if (istate.eq.nizs+1) then 

                  do iz = 0,nizs
                     tmpplot(ik,ir) = tmpplot(ik,ir) + 
     >                 tizs(ik,ir,iz) * mfact * 1.0/wi_sxb(ktebs(ik,ir))
                  end do
               else
                  tmpplot(ik,ir) = tizs(ik,ir,istate)*mfact
     >                             * 1.0/wi_sxb(ktebs(ik,ir))
c                  if (tizs(ik,ir,istate).gt.0.0) then 
c                     write(0,'(a,3i6,20(1x,g12.5))')
c     >                  'SXB:',ik,ir,istate,wi_sxb(ktebs(ik,ir)),
c     >                   ktebs(ik,ir),tizs(ik,ir,istate),
c     >                      tmpplot(ik,ir)
c                  endif

               endif
c     
            end do
c     
         end do   

      elseif (iselect.eq.42) then  
c     
c        Power balance terms
c         
         do ir = 1,nrs
c     
            do ik = 1, nks(ir)
c     
               tmpplot(ik,ir) = power_term(ik,ir,istate)
c
            end do
c
         end do
c         
      endif

c     
c     For cells which are exceptionally small compared
c     to the rest of the grid - zero out their contents. 
c     
c     Define exceptionally small as less than 1.0e-5 of 
c     the adjacent cells on the grid. 
c     
      zero_fact = 1.0e-5 
c     
      do ir = 1,nrs
         do ik = 1,nks(ir)

            if (ik.eq.1) then 
               if (kareas(ik,ir).lt.
     >              (zero_fact*kareas(ik+1,ir))) then 
c     
                  tmpplot(ik,ir) = 0.0
                  write(6,'(a,2i6,1p,4(1x,g15.8))')
     >                 'WARNING: DIVDATA:'//
     >                 ' VALUE ZEROED DUE TO GRID SIZE:',ik,ir,
     >                 kareas(ik,ir),kareas(ik+1,ir) 
                  write(0,'(a,2i6,1p,4(1x,g15.8))')
     >                 'WARNING: DIVDATA:'//
     >                 ' VALUE ZEROED DUE TO GRID SIZE:',ik,ir,
     >                 kareas(ik,ir),kareas(ik+1,ir) 
c     
               endif
            elseif (ik.eq.nks(ir)) then 
               if (kareas(ik,ir).lt.
     >              (zero_fact*kareas(ik-1,ir))) then 
c     
                  tmpplot(ik,ir) = 0.0
                  write(6,'(a,2i6,1p,4(1x,g15.8))')
     >                 'WARNING: DIVDATA:'//
     >                 ' VALUE ZEROED DUE TO GRID SIZE:',ik,ir,
     >                 kareas(ik,ir),kareas(ik-1,ir) 
                  write(0,'(a,2i6,1p,4(1x,g15.8))')
     >                 'WARNING: DIVDATA:'//
     >                 ' VALUE ZEROED DUE TO GRID SIZE:',ik,ir,
     >                 kareas(ik,ir),kareas(ik-1,ir) 
c     
               endif
            else
               if (kareas(ik,ir).lt.
     >              (zero_fact*kareas(ik+1,ir)).and.
     >              kareas(ik,ir).lt.
     >              (zero_fact*kareas(ik-1,ir))) then
c     
                  tmpplot(ik,ir) = 0.0
                  write(6,'(a,2i6,1p,4(1x,g15.8))')
     >                 'WARNING: DIVDATA:'//
     >                 ' VALUE ZEROED DUE TO GRID SIZE:',ik,ir,
     >                 kareas(ik+1,ir), kareas(ik,ir),kareas(ik-1,ir) 
                  write(0,'(a,2i6,1p,4(1x,g15.8))')
     >                 'WARNING: DIVDATA:'//
     >                 ' VALUE ZEROED DUE TO GRID SIZE:',ik,ir,
     >                 kareas(ik+1,ir), kareas(ik,ir),kareas(ik-1,ir) 
c     
               endif            

            endif

         end do

      end do

c     
c     DEBUG: 
c     
c     write(6,'(a,i6,a,i6)') 'DIVDATA: ISELECT = ', iselect, 
c     >                       ' ISTATE = ', istate
c     pltmax = -HI
c     pltmin =  HI 
c     
c     do ir = 1,nrs
c     
c     do ik = 1, nks(ir)
c     
c     write(6,'(a,2i6,1p,g15.8)') 'DATA:', ik, ir,
c     >            tmpplot(ik,ir)
c     pltmax = max(tmpplot(ik,ir),pltmax)          
c     pltmin = min(tmpplot(ik,ir),pltmin)          
c     
c     
c     end do
c     
c     end do   
c     write(6,'(a,1p,2(1x,g15.8))') 'DIVDATA: MAX/MIN = ',
c     >           pltmax,pltmin
c     
c     
c     As a post-processing function - copy the values from adjacent rings
c     into the virtual/boundary rings of the grid
c     
c     The boundary rings are IR=1, IR=IRWALL and IR=IRTRAP
c     
c     Look outward for the first ring
c     
      ir =1 
      do ik = 1,nks(ir)
         ik1 = ikouts(ik,ir)
         ir1 = irouts(ik,ir)
         tmpplot(ik,ir) = tmpplot(ik1,ir1)
      end do
c     
c     Look inward from irwall
c     
      ir = irwall
      do ik = 1,nks(ir)
         ik1 = ikins(ik,ir)
         ir1 = irins(ik,ir)
         tmpplot(ik,ir) = tmpplot(ik1,ir1)
      end do
c     
c     Look outward from irtrap if irtrap is not equal to or greater than nrs
c     
      if (irtrap.lt.nrs) then 
         ir = irtrap
         do ik = 1,nks(ir)
            ik1 = ikouts(ik,ir)
            ir1 = irouts(ik,ir)
            tmpplot(ik,ir) = tmpplot(ik1,ir1)
         end do
      endif
c     
      return
      end
c
c
c
      real function wi_sxb(te)
      implicit none
      real te
c
c     jde - SXB formula for WI obtained from spreadsheet from Tyler Abraams
c      
      wi_sxb = 53.1 * (1.0 - 1.04 * exp (-te/22.1))
      if (wi_sxb.lt.0.0) wi_sxb=0.0
      return 
      end
c
c
      subroutine load_divdata_targ(iselect,istate,ir,
     >                  start_targ_val,end_targ_val,ierr)
      use mod_params
      use mod_cgeom
      use mod_dynam2
      use mod_dynam3
      use mod_comtor
      use mod_pindata
      use mod_reiser_com
      use mod_slcom
      use mod_cedge2d
      use mod_adas_data_spec
      implicit none
c     include 'params' 
c     include 'cgeom'
c     include 'dynam2'
c     include 'dynam3'
c     include 'comtor'
c     include 'pindata'
c
c     include 'reiser_com'
c 
c     include 'slcom'
c     include 'cedge2d'
c     include 'adas_data_spec'
c
      integer :: iselect,istate,ierr,ir
      real :: start_targ_val,end_targ_val
      real :: cs
      real, external :: power_term
c     
c     LOAD_DIVTARG_DATA
c
c     This routine loads the two specified variables with the 
c     the values at the two ends of the specified ring. The start
c     target value if for S=0 while the end target value is for 
c     S=SMAX. For quantities without target values - IERR is set to 1. 
c
c     This routine has the same ISELECT and STATE options as 
c     the load_divdata_array routine
c
c     The allowed values of ISELECT are:
c
c     ISELECT = 1 = TOTAL H POWER LOSS  (W)
c               2 = TOTAL IMPURITY POWER LOSS  (W)
c               3 = TOTAL POWER LOSS   (W)
c               4 = SPECIFIED IMPURITY SPECTROSCOPIC LINE 
c                   - NEED TO READ ADAS DATA
c               5 = SPECIFIED HYDROGENIC SPECTROSCOPIC LINE 
c                   - NEED TO READ ADAS DATA
c               6 = PIN Halpha from PINALPHA array
c               7 = PIN HALPHA - By Component from Eirene - 6 for total
c                   - state specifies component
c                     1 - H ionisation
c                     2 - H+ recombination
c                     3 - H2 dissociation
c                     4 - H2+ dissociation
c                     5 - CX of H and H+
c                     6 - TOTAL 
c               8 = PIN HGAMMA - By component from Eirene - 6 for total
c                   - as above 
c               9 = Hydrogen Neutral Density 
c              10 = Background Plasma Properties
c                   1 = density
c                   2 = electron temperature
c                   3 = ion temperature
c                   4 = velocity
c                   5 = electric field
c              11 = Impurity Species Density - specified by charge state
c              12 = Impurity Species Temperature - specified by charge state
c              13 = Impurity Species Velocity - specified by charge state
c              14 = TOTAL H POWER LOSS (W/m3)
c              15 = TOTAL IMPURITY POWER LOSS (W/m3)
c              16 = TOTAL POWER LOSS (W/m3)
c              17 = Load PLRP (Particular Line Radiation Profile - see PLRP 
c                   module for istate values.
c              18 = Fluid code Background Plasma Properties
c                   1 = density
c                   2 = electron temperature
c                   3 = ion temperature
c                   4 = velocity
c                   5 = electric field
c              19 = Fluid code Impurity Species Density - specified by charge state
c              20 = Fluid code Impurity Species Temperature - specified by charge state
c              21 = Fluid code Impurity Species Velocity - specified by charge state
c              22 = SPECIFIED IMPURITY SPECTROSCOPIC LINE AVERAGED TEMPERATURE
c                   - MAY NEED TO READ ADAS DATA
c              23 = Impurity Density to Background Ne Ratio
c                   Istate = IZ
c              24 = Impurity Temperature to Background Te Ratio
c                   Istate = IZ
c              25 = Impurity Velocity to Background Vb Ratio
c                   Istate = IZ
c
c
c     Initialization
c
      start_targ_val = 0.0
      end_targ_val =0.0
      ierr=0

c
c     Check if core ring which doesn't have targets 
c
      if (ir.lt.irsep) then
         ierr = 1
         return
      endif

c
c     Check for valid ISELECT as input
c
      if (.not.(iselect.eq.10.or.iselect.eq.18.or.iselect.eq.42)) then 
c      if (iselect.lt.1.or.iselect.gt.25) then 
         write(6,'(a,i5)') 'LOAD_DIVDATA_TARG:INVALID SELECTOR:',iselect
         ierr = 1
         return
      endif
c
c     Options without target values 
c
      if (iselect.eq.1.or.iselect.eq.2.or.
     >    iselect.eq.3.or.iselect.eq.4.or.
     >    iselect.eq.5.or.iselect.eq.6.or.
     >    iselect.eq.7.or.iselect.eq.8.or.
     >    iselect.eq.9.or.
     >    iselect.eq.11.or.iselect.eq.12.or.
     >    iselect.eq.13.or.iselect.eq.14.or.
     >    iselect.eq.15.or.iselect.eq.16.or.
     >    iselect.eq.17.or.
     >    iselect.eq.19.or.iselect.eq.20.or.
     >    iselect.eq.21.or.iselect.eq.22.or.
     >    iselect.eq.23.or.iselect.eq.24.or.
     >    iselect.eq.25.or.iselect.eq.26
     >    ) then 
c
c         Set ierr =1 for no data
c
          ierr =1 
c
c     DIVIMP plasma background properties - target conditions
c
      elseif (iselect.eq.10) then 
c
c        Ne 
c
         if (istate.eq.1) then 
            start_targ_val = knds(idds(ir,2))
            end_targ_val   = knds(idds(ir,1))
c
c        Te 
c
         elseif (istate.eq.2) then 
            start_targ_val = kteds(idds(ir,2))
            end_targ_val   = kteds(idds(ir,1))
c
c        Ti 
c
         elseif (istate.eq.3) then 
            start_targ_val = ktids(idds(ir,2))
            end_targ_val   = ktids(idds(ir,1))
c
c        Vb
c
         elseif (istate.eq.4) then 
            start_targ_val = kvds(idds(ir,2))
            end_targ_val   = kvds(idds(ir,1))
c
c        E-field
c
         elseif (istate.eq.5) then 
            start_targ_val = keds(idds(ir,2))
            end_targ_val   = keds(idds(ir,1))
c
c        Mach number
c
         elseif (istate.eq.6) then 
            ! in most cases the target mach number should be 1.0
            CS = 9.79E3 * SQRT (0.5*(KTEDS(idds(IR,2))
     >                              +KTIDS(idds(IR,2)))*
     >                (1.0+RIZB)/CRMB)
            if (cs.ne.0.0) then 
               start_targ_val = kvds(idds(ir,2))/cs
            else
               start_targ_val = 0.0
            endif   

            CS = 9.79E3 * SQRT (0.5*(KTEDS(idds(IR,1))
     >                              +KTIDS(idds(IR,1)))*
     >                (1.0+RIZB)/CRMB)
            if (cs.ne.0.0) then 
               end_targ_val = kvds(idds(ir,1))/cs
            else
               end_targ_val = 0.0
            endif   

         endif
c
c     Fluid code background properties - target conditions
c
      elseif (iselect.eq.18) then 
c
c        Ne 
c
         if (istate.eq.1) then 
            start_targ_val = e2dtarg(ir,1,2)
            end_targ_val   = e2dtarg(ir,1,1)
c
c        Te 
c
         elseif (istate.eq.2) then 
            start_targ_val = e2dtarg(ir,2,2)
            end_targ_val   = e2dtarg(ir,2,1)
c
c        Ti 
c
         elseif (istate.eq.3) then 
            start_targ_val = e2dtarg(ir,3,2)
            end_targ_val   = e2dtarg(ir,3,1)
c
c        Vb
c
         elseif (istate.eq.4) then 
            start_targ_val = e2dtarg(ir,4,2)
            end_targ_val   = e2dtarg(ir,4,1)
c
c        E-field
c
         elseif (istate.eq.5) then 
            ierr =1 
         endif
c

      elseif (iselect.eq.42) then
         ! target power terms
         start_targ_val = power_term(0,ir,istate)
         end_targ_val   = power_term(nks(ir)+1,ir,istate)

      endif
c
c
c
      return
      end
c
c
c
      subroutine set_ylab(iselect,istate,itype,nizs,ylab)
      implicit none
      integer iselect,istate,itype,nizs
      character*(*) ylab
c
c     SET_YLAB:
c
c      
c     This is a support routine to the 2D DIVIMP data loading and 
c     integration code. Depending on the values of iselect,istate and 
c     itype - this routine sets the y axis label to something 
c     reasonable.
c
c     Itype specifies the type of plot - 0 = contour, 1 = integrated
c
c  
      integer len,lenstr
      external lenstr
c
c----------------------------------------------------------
c     Hydrogen power loss
c----------------------------------------------------------
c       
      if (iselect.eq.1) then
c
         if (istate.eq.0) then
            YLAB = 'H-NEUTRAL POW LOSS (BOLO)'
         elseif (istate.eq.1) then 
            YLAB = 'H+ POW LOSS (BOLO)'
         else 
            YLAB = 'TOT H POW LOSS (BOLO)'
         endif 
c
         len = lenstr(ylab)
c
         if (itype.eq.0) then 
            ylab = ylab(1:len) // ' (W)'
         elseif (itype.eq.1) then 
            ylab = ylab(1:len) // ' (W*M)'
         endif
c
c----------------------------------------------------------
c     Impurity power loss 
c----------------------------------------------------------
c
      elseif (iselect.eq.2) then
c
c        Individual charge state 
c
         if (istate.ge.0.and.istate.le.nizs) then 
c
            write(ylab,'(''IMP POW LOSS IZ ='',i4)')
     >                                          istate
c
c        Total Impurity
c
         else 
c
            YLAB = 'TOT IMP POW LOSS (BOLO) '
c
         endif 
c
         len = lenstr(ylab)
c
         if (itype.eq.0) then 
            ylab = ylab(1:len) // ' (W)'
         elseif (itype.eq.1) then 
            ylab = ylab(1:len) // ' (W*M)'
         endif
c
c----------------------------------------------------------
c     Total power loss 
c----------------------------------------------------------
c
      elseif (iselect.eq.3) then
c
         YLAB = 'TOT H+IMP POW LOSS (BOLO)'
c
         len = lenstr(ylab)
c
         if (itype.eq.0) then 
            ylab = ylab(1:len) // ' (W)'
         elseif (itype.eq.1) then 
            ylab = ylab(1:len) // ' (W*M)'
         endif
c
c----------------------------------------------------------
c     Line Radiation:
c          iselect = 4 = ADAS impurity line
c          iselect = 5 = ADAS hydrogen line
c          iselect = 6 = PIN total Halpha
c          iselect = 7 = EIRENE Halpha  istate =6 = total
c          iselect = 8 = EIRENE Hgamma  istate =6 = total
c          iselect = 17= PLRP radiation profile
c          iselect = 26= EIRENE Hbeta  istate =6 = total
c          iselect = 34= Subgrid based ADAS impurity line
c----------------------------------------------------------
c
      elseif (iselect.eq.4.or.iselect.eq.5.or.
     >        iselect.eq.6.or.iselect.eq.7.or.
     >        iselect.eq.8.or.iselect.eq.17.or.
     >        iselect.eq.26.or.iselect.eq.34) then
c
         YLAB = 'LINE RADIATION (PHOTONS'
c
         len = lenstr(ylab)
c
         if (itype.eq.0) then 
            ylab = ylab(1:len) // '/M^3/S)'
         elseif (itype.eq.1) then 
            ylab = ylab(1:len) // '/M^2/S)'
         elseif (itype.eq.2) then 
            ylab = ylab(1:len) // '/M^2/S/SR)'
         elseif (itype.eq.3) then 
            ylab = ylab(1:len) // '/M^3/S/SR)'
         elseif (itype.eq.4) then 
            ylab = ylab(1:len) // '/CM^3/S/SR)'
         endif
c
c----------------------------------------------------------
c
c     PIN - Hydrogen Neutral Density 
c
c----------------------------------------------------------
c
      elseif (iselect.eq.9) then   

         YLAB = 'DENSITY ('
c
         len = lenstr(ylab)
c
         if (itype.eq.0) then 
            ylab = ylab(1:len) // '/M^3)'
         elseif (itype.eq.1) then 
            ylab = ylab(1:len) // '/M^2)'
         endif
c
c----------------------------------------------------------
c
c     DIVIMP - Background Plasma Properties
c
c----------------------------------------------------------
c
      elseif (iselect.eq.10) then   

         if (istate.eq.1) then 
            YLAB = 'DENSITY (M^-3)'
         elseif(istate.eq.2) then 
            YLAB = 'ELEC TEMPERATURE (eV)'
         elseif(istate.eq.3) then 
            YLAB = 'ION TEMPERATURE (eV)'
         elseif(istate.eq.4) then 
            YLAB = 'VELOCITY (M/S)'
         elseif(istate.eq.5) then 
            YLAB = 'ELECTRIC FIELD (V/M(?))'
         elseif(istate.eq.6) then 
            YLAB = 'MACH NUMBER'
         endif
c
c----------------------------------------------------------
c
c     DIVIMP - 11 = DIVIMP Impurity Density
c              32 = SUBGRID Impurity Density
c
c----------------------------------------------------------
c
      elseif (iselect.eq.11.or.iselect.eq.32) then   

         write(YLAB,'(''IMP DENSITY: STATE='',i4,
     >                ''(M^-3)'')') istate
c
c
c----------------------------------------------------------
c
c     DIVIMP - Impurity Temperature
c
c----------------------------------------------------------
c
      elseif (iselect.eq.12) then   

         write(YLAB,'(''IMP TEMPERATURE: STATE='',i4,
     >                ''(eV)'')') istate
c
c
c----------------------------------------------------------
c
c     DIVIMP - Impurity Velocity
c
c----------------------------------------------------------
c
      elseif (iselect.eq.13) then   

         write(YLAB,'(''IMP VELOCITY: STATE='',i4,
     >                ''(M/S)'')') istate

c
c----------------------------------------------------------
c     Hydrogen power loss  (W/m3)
c----------------------------------------------------------
c       
      elseif (iselect.eq.14.or.iselect.eq.37) then
c
         if (istate.eq.0) then
            YLAB = 'H-NEUTRAL POW LOSS (BOLO)'
         elseif (istate.eq.1) then 
            YLAB = 'H+ POW LOSS (BOLO)'
         else 
            YLAB = 'TOTAL H POW LOSS (BOLO)'
         endif 
c
         len = lenstr(ylab)
c
         if (itype.eq.0) then 
            ylab = ylab(1:len) // ' (W/M3)'
         elseif (itype.eq.1) then 
            ylab = ylab(1:len) // ' (W/M2)'
         elseif (itype.eq.2) then 
            ylab = ylab(1:len) // ' (W)'
         endif
c
c----------------------------------------------------------
c     Impurity power loss (W/m3)
c----------------------------------------------------------
c
      elseif (iselect.eq.15) then
c
c        Individual charge state 
c
         if (istate.ge.0.and.istate.le.nizs) then 
c
            write(ylab,'(''IMP POW LOSS IZ ='',i4)')
     >                                          istate
c
c        Total Impurity
c
         else 
c
            YLAB = 'TOTAL IMP POW LOSS (BOLO) '
c
         endif 
c
         len = lenstr(ylab)
c
         if (itype.eq.0) then 
            ylab = ylab(1:len) // ' (W/M3)'
         elseif (itype.eq.1) then 
            ylab = ylab(1:len) // ' (W/M2)'
         elseif (itype.eq.2) then 
            ylab = ylab(1:len) // ' (W)'
         endif
c
c----------------------------------------------------------
c     Total power loss  (W/m3)
c----------------------------------------------------------
c
      elseif (iselect.eq.16) then
c
         YLAB = 'TOT H+IMP POW LOSS (BOLO)'
c
         len = lenstr(ylab)
c
         if (itype.eq.0) then 
            ylab = ylab(1:len) // ' (W/M3)'
         elseif (itype.eq.1) then 
            ylab = ylab(1:len) // ' (W/M2)'
         elseif (itype.eq.2) then 
            ylab = ylab(1:len) // ' (W)'
         endif


c
c----------------------------------------------------------
c
c     FLUID CODE - Background Plasma Properties
c
c----------------------------------------------------------
c
      elseif (iselect.eq.18) then   

         if (istate.eq.1) then 
            YLAB = 'FC DENSITY (M^-3)'
         elseif(istate.eq.2) then 
            YLAB = 'FC ELEC TEMPERATURE (eV)'
         elseif(istate.eq.3) then 
            YLAB = 'FC ION TEMPERATURE (eV)'
         elseif(istate.eq.4) then 
            YLAB = 'FC VELOCITY (M/S)'
         elseif(istate.eq.5) then 
            YLAB = 'FC ELECTRIC FIELD (V/M(?))'
         endif
c
c----------------------------------------------------------
c
c     FLUID CODE - Impurity Density
c
c----------------------------------------------------------
c
      elseif (iselect.eq.19) then   

         write(YLAB,'(''FC IMP DENSITY: STATE='',i4,
     >                ''(M^-3)'')') istate
c
c
c----------------------------------------------------------
c
c     FLUID CODE - Impurity Temperature
c
c----------------------------------------------------------
c
      elseif (iselect.eq.20) then   

         write(YLAB,'(''FC IMP TEMPERATURE: STATE='',i4,
     >                ''(eV)'')') istate
c
c
c----------------------------------------------------------
c
c     FLUID CODE - Impurity Velocity
c
c----------------------------------------------------------
c
      elseif (iselect.eq.21) then   

         write(YLAB,'(''FC IMP VELOCITY: STATE='',i4,
     >                ''(M/S)'')') istate

c
c----------------------------------------------------------
c
c     DIVIMP - Emission Weighted impurity temperature
c
c----------------------------------------------------------

      elseif (iselect.eq.22) then 
c
         YLAB = 'EMISSION WEIGHTED AV. ION TEMP (eV)'
c
c
c----------------------------------------------------------
c
c     DIVIMP - Impurity Density Ratio to Background Plasma
c
c----------------------------------------------------------
c
      elseif (iselect.eq.23) then   

         write(YLAB,'(''IMP DENSITY RATIO: STATE='',i4
     >                )') istate
c
c
c----------------------------------------------------------
c
c     DIVIMP - Impurity Temperature Ratio to BG Plasma
c
c----------------------------------------------------------
c
      elseif (iselect.eq.24) then   

         write(YLAB,'(''IMP TEMPERATURE RATIO: STATE='',i4
     >                )') istate
c
c
c----------------------------------------------------------
c
c     DIVIMP - Impurity Velocity Ratio to Background Plasma
c
c----------------------------------------------------------
c
      elseif (iselect.eq.25) then   

         write(YLAB,'(''IMP VELOCITY RATIO: STATE='',i4
     >               )') istate
c
c----------------------------------------------------------
c
c     HC Related quantities
c
c     HC data - 28 = HC CD Emission
c               29 = HC state density
c               30 = HC state ionization
c               33 = SUBGRID HC state density
c               35 = SUBGRID HC emission
c
c----------------------------------------------------------
c
      elseif (iselect.eq.28.or.iselect.eq.29.or.iselect.eq.30.or.
     >        iselect.eq.33.or.iselect.eq.35) then 
c
         call hc_set_ylab(iselect,istate,itype,nizs,ylab)
c
c
c----------------------------------------------------------
c
c     DIVIMP - Impurity Ionization
c
c----------------------------------------------------------
c
      elseif (iselect.eq.31) then   

         write(YLAB,'(''IMP IONIZATION: STATE='',i4,
     >                ''(M^-3)'')') istate
 
      elseif (iselect.eq.36) then   
c
c                   1 = PINION = PIN ionization    
c                   2 = PINATOM = PIN Atom density 
c                   3 = PINMOL = PIN Molecular density
c                   4 = PINIONZ = Impurity ionization
c                   5 = PINZ0 = Impurity neutral density  
c                   6 = PINQI = Ion heating term
c                   7 = PINQE = Electron heating term

         if (istate.eq.1) then 
            YLAB = 'PIN IZ   (/M^3/S)'
         elseif (istate.eq.2) then 
            YLAB = 'PIN ATOM (/M^3)'
         elseif (istate.eq.3) then 
            YLAB = 'PIN MOL  (/M^3)'
         elseif (istate.eq.4) then 
            YLAB = 'PIN ZIZ  (/M^3/S)'
         elseif (istate.eq.5) then 
            YLAB = 'PIN ZDEN (/M^3)'
         elseif (istate.eq.6) then 
            YLAB = 'PIN QI   (W/M^3)'
         elseif (istate.eq.7) then 
            YLAB = 'PIN QE   (W/M^3)'
         endif
c
         len = lenstr(ylab)
c
c         if (itype.eq.0) then 
c            ylab = ylab(1:len) // '/M^3)'
c         elseif (itype.eq.1) then 
c            ylab = ylab(1:len) // '/M^2)'
c         endif

      elseif (iselect.eq.40) then 
c
c         ExB drift related quantities
c         1 - Potential (phi) (V) 
c         2 - Radial Efield (V/m)
c         3 - Poloidal Efield (V/m)
c         4 - Radial ExB drift (m/s)
c         5 - Poloidal ExB drift (m/s)
c         6 - Radial ExB flux   ne x Vexb_rad (/m2/s)
c         7 - Poloidal ExB flux ne x Vexb_pol (/m2/s)
c     
         if (istate.eq.1) then 
            YLAB = 'E-POTENTIAL (V)'
         elseif (istate.eq.2) then 
            YLAB = 'E-RADIAL (V/m)'
         elseif (istate.eq.3) then 
            YLAB = 'E-POLOIDAL (V/m)'
         elseif (istate.eq.4) then 
            YLAB = 'ExB RADIAL DRIFT (m/s)'
         elseif (istate.eq.5) then 
            YLAB = 'ExB POLOIDAL DRIFT (m/s)'
         elseif (istate.eq.6) then 
            YLAB = 'ExB RADIAL FLUX (/m2/s)'
         elseif (istate.eq.7) then 
            YLAB = 'ExB POLOIDAL FLUX (/m2/s)'
         endif

      elseif (iselect.eq.41) then   

         write(YLAB,'(''W0 400.6 EMIS.: STATE='',i4,
     >                ''(PH/M2/S)'')') istate

      elseif (iselect.eq.42) then   

         if (istate.eq.1) then 
            YLAB = 'I-CONDUCTION'
         elseif(istate.eq.2) then 
            YLAB = 'I-CONVECTION'
         elseif(istate.eq.3) then 
            YLAB = 'E-CONDUCTION'
         elseif(istate.eq.4) then 
            YLAB = 'E-CONVECTION'
         elseif(istate.eq.5) then 
            YLAB = 'TOTAL CONDUCTION'
         elseif(istate.eq.6) then 
            YLAB = 'TOTAL CONVECTION'
         elseif(istate.eq.7) then 
            YLAB = 'CONV/COND'
         elseif(istate.eq.8) then 
            YLAB = 'CONV/E-COND'
         endif
 
      endif

c
      return
      end
c
c
c
      subroutine set_blab(iselect,istate,itype,nizs,blab)
      implicit none
      integer iselect,istate,itype,nizs
      character*(*) blab
c
c     SET_BLAB:
c
c      
c     This is a support routine to the 2D DIVIMP data loading and 
c     integration code. Depending on the values of iselect,istate and 
c     itype - this routine sets the plot label to something 
c     reasonable.
c
c     Itype specifies the type of plot - 0 = contour, 1 = integrated
c
c  
      integer len,lenstr
      external lenstr
c
c----------------------------------------------------------
c     Hydrogen power loss
c----------------------------------------------------------
c       
      if (iselect.eq.1) then
c
         if (itype.eq.0) then           
            BLAB = 'H POW LOSS (BOLO)'
         elseif (itype.eq.1) then 
            BLAB = 'CODE H POW LOSS (BOLO)'
         endif
c
c----------------------------------------------------------
c     Impurity power loss 
c----------------------------------------------------------
c
      elseif (iselect.eq.2) then 
c
         if (itype.eq.0) then           
            BLAB = 'IMP POW LOSS (BOLO)'
         elseif (itype.eq.1) then 
            BLAB = 'BOLO IMP POW LOSS'
         endif
c
c----------------------------------------------------------
c     Total power loss 
c----------------------------------------------------------
c
      elseif (iselect.eq.3) then 
c
         if (itype.eq.0) then           
            BLAB = 'TOTAL POW LOSS (BOLO)'
         elseif (itype.eq.1) then 
            BLAB = 'BOLO TOTAL POW LOSS'
         endif
c  
c----------------------------------------------------------
c     ADAS IMPURITY PLRP 
c----------------------------------------------------------
c
c     4 = ADAS Impurity Emission
c    34 = SUBGRID ADAS Impurity Emission 
c
      elseif (iselect.eq.4.or.iselect.eq.34) then  
c
         if (itype.eq.0) then           
            BLAB = 'ADAS IMP PLRP'
         elseif (itype.eq.1.or.itype.eq.2.or.itype.eq.3.or.
     >           itype.eq.4) then 
            BLAB = 'CODE ADAS IMP PLRP'
         endif
c  
c----------------------------------------------------------
c     ADAS HYDROGENIC PLRP 
c----------------------------------------------------------
c
      elseif (iselect.eq.5) then 
c   
         if (itype.eq.0) then           
            BLAB = 'ADAS HYDROGENIC PLRP'
         elseif (itype.eq.1.or.itype.eq.2.or.itype.eq.3.or.
     >           itype.eq.4) then 
            BLAB = 'CODE ADAS H PLRP'
         endif
c  
c----------------------------------------------------------
c     PIN TOTAL HALPHA
c----------------------------------------------------------
c
      elseif (iselect.eq.6) then 
c   
         if (itype.eq.0) then           
            BLAB = 'CODE HALPHA'
         elseif (itype.eq.1.or.itype.eq.2.or.itype.eq.3.or.
     >           itype.eq.4) then 
            BLAB = 'CODE CODE HALPHA'
         endif
c  
c----------------------------------------------------------
c     EIRENE HALPHA
c----------------------------------------------------------
c
      elseif (iselect.eq.7) then 
c   
         if (istate.eq.6) then 
            if (itype.eq.0) then           
               BLAB = 'EIRENE TOTAL HALPHA'
            elseif (itype.eq.1.or.itype.eq.2.or.itype.eq.3.or.
     >           itype.eq.4) then 
               BLAB = 'CODE EIRENE TOT HALPHA'
            endif
         else
            if (itype.eq.0) then           
               write(blab,'(a,i4)') 'EIRENE HALPHA COMP=',istate
            elseif (itype.eq.1.or.itype.eq.2.or.itype.eq.3.or.
     >           itype.eq.4) then 
               write(blab,'(a,i4)') 'CODE EIRENE HALPHA COMP=',istate
            endif
         endif
c  
c----------------------------------------------------------
c     EIRENE HGAMMA 
c----------------------------------------------------------
c
      elseif (iselect.eq.8) then 
c   
         if (istate.eq.6) then 
            if (itype.eq.0) then           
               BLAB = 'EIRENE TOTAL HGAMMA'
            elseif (itype.eq.1.or.itype.eq.2.or.itype.eq.3.or.
     >           itype.eq.4) then 
               BLAB = 'CODE EIRENE TOT HGAMMA'
            endif
         else
            if (itype.eq.0) then           
               write(blab,'(a,i4)') 'EIRENE HGAMMA COMP=',istate
            elseif (itype.eq.1.or.itype.eq.2.or.itype.eq.3.or.
     >           itype.eq.4) then 
               write(blab,'(a,i4)') 'CODE EIRENE HGAMMA COMP=',istate
            endif
         endif
c
c
c----------------------------------------------------------
c     PIN - Hydrogen Neutral Density 
c----------------------------------------------------------
c
      elseif (iselect.eq.9) then   

         BLAB = 'PIN NEUTRAL H DENSITY'
c
c
c----------------------------------------------------------
c
c     DIVIMP - Background Plasma Properties
c
c----------------------------------------------------------
c
      elseif (iselect.eq.10) then   

         if (istate.eq.1) then 
            BLAB = 'BG ION DENSITY'
         elseif(istate.eq.2) then 
            BLAB = 'BG ELECTRON TEMPERATURE'
         elseif(istate.eq.3) then 
            BLAB = 'BG ION TEMPERATURE'
         elseif(istate.eq.4) then 
            BLAB = 'BG VELOCITY'
         elseif(istate.eq.5) then 
            BLAB = 'BG ELECTRIC FIELD'
         elseif(istate.eq.6) then 
            BLAB = 'BG MACH NUMBER'
         endif
c
c----------------------------------------------------------
c
c     DIVIMP - Impurity Density
c
c     11 = DIVIMP Impurity Density
c     32 = SUBGRID Impurity Density
c
c----------------------------------------------------------
c
      elseif (iselect.eq.11.or.iselect.eq.32) then   

         write(BLAB,'(''IMP DENSITY: STATE='',i4,
     >                ''(M^-3)'')') istate
c
c
c----------------------------------------------------------
c
c     DIVIMP - Impurity Temperature
c
c----------------------------------------------------------
c
      elseif (iselect.eq.12) then   

         write(BLAB,'(''IMP TEMPERATURE: STATE='',i4,
     >                ''(eV)'')') istate
c
c
c----------------------------------------------------------
c
c     DIVIMP - Impurity Velocity
c
c----------------------------------------------------------
c
      elseif (iselect.eq.13) then   

         write(BLAB,'(''IMP VELOCITY: STATE='',i4,
     >                ''(M/S)'')') istate

c
c----------------------------------------------------------
c     Hydrogen power loss (W/m3)
c----------------------------------------------------------
c       
      elseif (iselect.eq.14.or.iselect.eq.37) then
c
         if (itype.eq.0) then           
            BLAB = 'H POW LOSS (BOLO)'
         elseif (itype.eq.1.or.itype.eq.2) then 
            BLAB = 'CODE H POW LOSS (BOLO)'
         endif
c
c----------------------------------------------------------
c     Impurity power loss 
c----------------------------------------------------------
c
      elseif (iselect.eq.15) then 
c
         if (itype.eq.0) then           
            BLAB = 'IMP POW LOSS (BOLO)'
         elseif (itype.eq.1.or.itype.eq.2) then 
            BLAB = 'CODE IMP POW LOSS'
         endif
c
c----------------------------------------------------------
c     Total power loss 
c----------------------------------------------------------
c
      elseif (iselect.eq.16) then 
c
         if (itype.eq.0) then           
            BLAB = 'TOTAL POW LOSS (BOLO)'
         elseif (itype.eq.1.or.itype.eq.2) then 
            BLAB = 'CODE TOTAL POW LOSS'
         endif
c
c
      elseif (iselect.eq.17) then  
c
         if (itype.eq.0) then           
            BLAB = 'CUSTOM IMP PLRP'
         elseif (itype.eq.1.or.itype.eq.2) then 
            BLAB = 'CODE CUSTOM IMP PLRP'
         endif
c
c----------------------------------------------------------
c
c     FLUID CODE - Background Plasma Properties
c
c----------------------------------------------------------
c
      elseif (iselect.eq.18) then   

         if (istate.eq.1) then 
            BLAB = 'FC ION DENSITY'
         elseif(istate.eq.2) then 
            BLAB = 'FC ELECTRON TEMPERATURE'
         elseif(istate.eq.3) then 
            BLAB = 'FC ION TEMPERATURE'
         elseif(istate.eq.4) then 
            BLAB = 'FC VELOCITY'
         elseif(istate.eq.5) then 
            BLAB = 'FC ELECTRIC FIELD'
         endif
c
c----------------------------------------------------------
c
c     FLUID CODE - Impurity Density
c
c----------------------------------------------------------
c
      elseif (iselect.eq.19) then   

         write(BLAB,'(''FC IMP DENSITY: STATE='',i4,
     >                ''(M^-3)'')') istate
c
c
c----------------------------------------------------------
c
c     FLUID CODE - Impurity Temperature
c
c----------------------------------------------------------
c
      elseif (iselect.eq.20) then   

         write(BLAB,'(''FC IMP TEMPERATURE: STATE='',i4,
     >                ''(eV)'')') istate
c
c
c----------------------------------------------------------
c
c     FLUID CODE - Impurity Velocity
c
c----------------------------------------------------------
c
      elseif (iselect.eq.21) then   

         write(BLAB,'(''FC IMP VELOCITY: STATE='',i4,
     >                ''(M/S)'')') istate


c  
c----------------------------------------------------------
c     DIVIMP - EMISSION WEIGHTED AVERAGE ION TEMPERATURE 
c----------------------------------------------------------
c
      elseif (iselect.eq.22) then  
c
         BLAB = 'ADAS-BASED AVERAGE ION TEMP'
c
c
c----------------------------------------------------------
c
c     DIVIMP - Impurity Density Ratio to BG Plasma
c
c----------------------------------------------------------
c
      elseif (iselect.eq.23) then   

         write(BLAB,'(''IMP DENSITY RATIO: STATE='',i4
     >                 )') istate
c
c
c----------------------------------------------------------
c
c     DIVIMP - Impurity Temperature Ratio to BG Plasma 
c
c----------------------------------------------------------
c
      elseif (iselect.eq.24) then   

         write(BLAB,'(''IMP TEMPERATURE RATIO: STATE='',i4
     >               )') istate
c
c
c----------------------------------------------------------
c
c     DIVIMP - Impurity Velocity Ratio to BG Plasma 
c
c----------------------------------------------------------
c
      elseif (iselect.eq.25) then   

         write(BLAB,'(''IMP VELOCITY RATIO: STATE='',i4
     >                )') istate

c  
c----------------------------------------------------------
c     EIRENE HBETA
c----------------------------------------------------------
c
      elseif (iselect.eq.26) then 
c   
         if (istate.eq.6) then 
            if (itype.eq.0) then           
               BLAB = 'EIRENE TOTAL HBETA'
            elseif (itype.eq.1.or.itype.eq.2) then 
               BLAB = 'CODE EIRENE TOT HBETA'
            endif
         else
            if (itype.eq.0) then           
               write(blab,'(a,i4)') 'EIRENE HBETA COMP=',istate
            elseif (itype.eq.1.or.itype.eq.2) then 
               write(blab,'(a,i4)') 'CODE EIRENE HBETA COMP=',istate
            endif
         endif
c
c----------------------------------------------------------
c
c     HC Related quantities
c
c----------------------------------------------------------
c
      elseif (iselect.eq.28.or.iselect.eq.29.or.iselect.eq.30.or.
     >        iselect.eq.33.or.iselect.eq.35) then 
c
         call hc_set_blab(iselect,istate,itype,nizs,blab)
c
c
c----------------------------------------------------------
c
c     DIVIMP - Impurity Ionization
c
c----------------------------------------------------------
c
      elseif (iselect.eq.31) then   

         write(BLAB,'(''IMP IONIZATION: STATE='',i4,
     >                ''(M^-3)'')') istate

      elseif (iselect.eq.36) then   
c
c                   1 = PINION = PIN ionization    
c                   2 = PINATOM = PIN Atom density 
c                   3 = PINMOL = PIN Molecular density
c                   4 = PINIONZ = Impurity ionization
c                   5 = PINZ0 = Impurity neutral density  
c                   6 = PINQI = Ion heating term
c                   7 = PINQE = Electron heating term

         if (istate.eq.1) then 
            BLAB = 'PIN IZ   (/M^3/S)'
         elseif (istate.eq.2) then 
            BLAB = 'PIN ATOM (/M^3)'
         elseif (istate.eq.3) then 
            BLAB = 'PIN MOL  (/M^3)'
         elseif (istate.eq.4) then 
            BLAB = 'PIN ZIZ  (/M^3/S)'
         elseif (istate.eq.5) then 
            BLAB = 'PIN ZDEN (/M^3)'
         elseif (istate.eq.6) then 
            BLAB = 'PIN QI   (W/M^3)'
         elseif (istate.eq.7) then 
            BLAB = 'PIN QE   (W/M^3)'
         endif
c
         len = lenstr(blab)
c
c         if (itype.eq.0) then 
c            blab = blab(1:len) // '/M^3)'
c         elseif (itype.eq.1) then 
c            blab = blab(1:len) // '/M^2)'
c         endif

      elseif (iselect.eq.40) then 
c
c         ExB drift related quantities
c         1 - Potential (phi) (V) 
c         2 - Radial Efield (V/m)
c         3 - Poloidal Efield (V/m)
c         4 - Radial ExB drift (m/s)
c         5 - Poloidal ExB drift (m/s)
c         6 - Radial ExB flux   ne x Vexb_rad (/m2/s)
c         7 - Poloidal ExB flux ne x Vexb_pol (/m2/s)
c     
         if (istate.eq.1) then 
            BLAB = 'E-POTENTIAL (V)'
         elseif (istate.eq.2) then 
            BLAB = 'E-RADIAL (V/m)'
         elseif (istate.eq.3) then 
            BLAB = 'E-POLOIDAL (V/m)'
         elseif (istate.eq.4) then 
            BLAB = 'ExB RADIAL DRIFT (m/s)'
         elseif (istate.eq.5) then 
            BLAB = 'ExB POLOIDAL DRIFT (m/s)'
         elseif (istate.eq.6) then 
            BLAB = 'ExB RADIAL FLUX (/m2/s)'
         elseif (istate.eq.7) then 
            BLAB = 'ExB POLOIDAL FLUX (/m2/s)'
         endif

      elseif (iselect.eq.41) then   

         write(BLAB,'(''W0  W0 400.6:ST='',i4,
     >                ''(PH/M2/S)'')') istate

      elseif (iselect.eq.42) then   

         if (istate.eq.1) then 
            BLAB = 'I-CONDUCTION'
         elseif(istate.eq.2) then 
            BLAB = 'I-CONVECTION'
         elseif(istate.eq.3) then 
            BLAB = 'E-CONDUCTION'
         elseif(istate.eq.4) then 
            BLAB = 'E-CONVECTION'
         elseif(istate.eq.5) then 
            BLAB = 'TOTAL CONDUCTION'
         elseif(istate.eq.6) then 
            BLAB = 'TOTAL CONVECTION'
         elseif(istate.eq.7) then 
            BLAB = 'CONV/COND'
         elseif(istate.eq.8) then 
            BLAB = 'CONV/E-COND'
         endif
         
      endif
c
c
c
      return
      end

      subroutine set_elab(iselect,istate,elab)
      implicit none
      integer iselect,istate,iz
      character*(*) elab
c
c     SET_BLAB:
c
c      
c     This is a support routine to the 2D DIVIMP data loading and 
c     integration code. Depending on the values of iselect,istate and 
c     itype - this routine sets the plot label to something 
c     reasonable.
c
c     Itype specifies the type of plot - 0 = contour, 1 = integrated
c
c  
      integer len,lenstr
      external lenstr
c
c----------------------------------------------------------
c     Hydrogen power loss
c----------------------------------------------------------
c       
      if (iselect.eq.1) then
c
c         if (itype.eq.0) then           
            ELAB = 'HpowHpow (bolo)'
c         elseif (itype.eq.1) then 
c            ELAB = 'CODE H POW LOSS (BOLO)'
c         endif
c
c----------------------------------------------------------
c     Impurity power loss 
c----------------------------------------------------------
c
      elseif (iselect.eq.2) then 
c
         ELAB = 'IpowIpow (bolo)'
c         if (itype.eq.0) then           
c            ELAB = 'IMP POW LOSS (BOLO)'
c         elseif (itype.eq.1) then 
c            ELAB = 'BOLO IMP POW LOSS'
c         endif
c
c----------------------------------------------------------
c     Total power loss 
c----------------------------------------------------------
c
      elseif (iselect.eq.3) then 
c
         ELAB = 'TpowTpow (bolo)'
c         if (itype.eq.0) then           
c            ELAB = 'TOTAL POW LOSS (BOLO)'
c         elseif (itype.eq.1) then 
c            ELAB = 'BOLO TOTAL POW LOSS'
c         endif
c  
c----------------------------------------------------------
c     ADAS IMPURITY PLRP 
c----------------------------------------------------------
c
c     4 = ADAS Impurity Emission
c    34 = SUBGRID ADAS Impurity Emission 
c
      elseif (iselect.eq.4.or.iselect.eq.34) then  
c
         ELAB = 'IradIrad (ADAS)'
c         if (itype.eq.0) then           
c            ELAB = 'ADAS IMP PLRP'
c         elseif (itype.eq.1.or.itype.eq.2.or.itype.eq.3.or.
c     >           itype.eq.4) then 
c            ELAB = 'CODE ADAS IMP PLRP'
c         endif
c  
c----------------------------------------------------------
c     ADAS HYDROGENIC PLRP 
c----------------------------------------------------------
c
      elseif (iselect.eq.5) then 
c   
         ELAB = 'HradHrad (ADAS)'
c         if (itype.eq.0) then           
c            ELAB = 'ADAS HYDROGENIC PLRP'
c         elseif (itype.eq.1.or.itype.eq.2.or.itype.eq.3.or.
c     >           itype.eq.4) then 
c            ELAB = 'CODE ADAS H PLRP'
c         endif
c  
c----------------------------------------------------------
c     PIN TOTAL HALPHA
c----------------------------------------------------------
c
      elseif (iselect.eq.6) then 
c   
         ELAB = 'Ha THa T  (PIN)'
c         if (itype.eq.0) then           
c            ELAB = 'CODE HALPHA'
c         elseif (itype.eq.1.or.itype.eq.2.or.itype.eq.3.or.
c     >           itype.eq.4) then 
c            ELAB = 'CODE CODE HALPHA'
c         endif
c  
c----------------------------------------------------------
c     EIRENE HALPHA
c----------------------------------------------------------
c
      elseif (iselect.eq.7) then 
c   
         if (istate.eq.6) then 
            ELAB = 'HA THA T (EIR)'
c            if (itype.eq.0) then           
c               ELAB = 'EIRENE TOTAL HALPHA'
c            elseif (itype.eq.1.or.itype.eq.2.or.itype.eq.3.or.
c     >           itype.eq.4) then 
c               ELAB = 'CODE EIRENE TOT HALPHA'
c            endif
         else
            
            write(elab,'(a,i2,a,i2)') 'HA',istate,'HA',istate

c            if (itype.eq.0) then           
c               write(elab,'(a,i4)') 'EIRENE HALPHA COMP=',istate
c            elseif (itype.eq.1.or.itype.eq.2.or.itype.eq.3.or.
c     >           itype.eq.4) then 
c               write(elab,'(a,i4)') 'CODE EIRENE HALPHA COMP=',istate
c            endif
         endif
c  
c----------------------------------------------------------
c     EIRENE HGAMMA 
c----------------------------------------------------------
c
      elseif (iselect.eq.8) then 
c   
         if (istate.eq.6) then 
            ELAB = 'HG THG T (EIR)'
c            if (itype.eq.0) then           
c               ELAB = 'EIRENE TOTAL HGAMMA'
c            elseif (itype.eq.1.or.itype.eq.2.or.itype.eq.3.or.
c     >           itype.eq.4) then 
c            elseif (itype.eq.1.or.itype.eq.2) then 
c               ELAB = 'CODE EIRENE TOT HGAMMA'
c            endif
         else
            write(elab,'(a,i2,a,i2)') 'HG',istate,'HG',istate
c            if (itype.eq.0) then           
c               write(elab,'(a,i4)') 'EIRENE HGAMMA COMP=',istate
c            elseif (itype.eq.1.or.itype.eq.2.or.itype.eq.3.or.
c     >           itype.eq.4) then 
c               write(elab,'(a,i4)') 'CODE EIRENE HGAMMA COMP=',istate
c            endif
         endif
c
c
c----------------------------------------------------------
c     PIN - Hydrogen Neutral Density 
c----------------------------------------------------------
c
      elseif (iselect.eq.9) then   

         ELAB = 'H0  H0'
c
c
c----------------------------------------------------------
c
c     DIVIMP - Background Plasma Properties
c
c----------------------------------------------------------
c
      elseif (iselect.eq.10) then   

         if (istate.eq.1) then 
            ELAB = 'BGNeBGNe   '
         elseif(istate.eq.2) then 
            ELAB = 'BGTeBGTe   '
         elseif(istate.eq.3) then 
            ELAB = 'BGTiBGTi   '
         elseif(istate.eq.4) then 
            ELAB = 'BGVbBGVb   '
         elseif(istate.eq.5) then 
            ELAB = 'BGEfBGEf   '
         elseif(istate.eq.6) then 
            ELAB = 'BGMaBGMa   '
         endif
c
c----------------------------------------------------------
c
c     DIVIMP - Impurity Density
c
c     11 = DIVIMP Impurity Density
c     32 = SUBGRID Impurity Density
c
c----------------------------------------------------------
c
      elseif (iselect.eq.11.or.iselect.eq.32) then   

         write(ELAB,'(''N'',i3,''N'',i3)')
     >                  istate,istate
c
c
c----------------------------------------------------------
c
c     DIVIMP - Impurity Temperature
c
c----------------------------------------------------------
c
      elseif (iselect.eq.12) then   

         write(ELAB,'(''T'',i3,''T'',i3)')
     >                  istate,istate
c
c
c----------------------------------------------------------
c
c     DIVIMP - Impurity Velocity
c
c----------------------------------------------------------
c
      elseif (iselect.eq.13) then   

         write(ELAB,'(''V'',i3,''V'',i3)')
     >                  istate,istate

c
c----------------------------------------------------------
c     Hydrogen power loss (W/m3)
c----------------------------------------------------------
c       
      elseif (iselect.eq.14.or.iselect.eq.37) then
c
          ELAB = 'HpowHpow (bolo)'
c         if (itype.eq.0) then           
c            ELAB = 'H POW LOSS (BOLO)'
c         elseif (itype.eq.1.or.itype.eq.2) then 
c            ELAB = 'CODE H POW LOSS (BOLO)'
c         endif
c
c----------------------------------------------------------
c     Impurity power loss 
c----------------------------------------------------------
c
      elseif (iselect.eq.15) then 
c
         ELAB = 'IpowIpow (bolo)'
c         if (itype.eq.0) then           
c            ELAB = 'IMP POW LOSS (BOLO)'
c         elseif (itype.eq.1.or.itype.eq.2) then 
c            ELAB = 'CODE IMP POW LOSS'
c         endif
c
c----------------------------------------------------------
c     Total power loss 
c----------------------------------------------------------
c
      elseif (iselect.eq.16) then 
c
         ELAB = 'TpowTpow (bolo)'
c         if (itype.eq.0) then           
c            ELAB = 'TOTAL POW LOSS (BOLO)'
c         elseif (itype.eq.1.or.itype.eq.2) then 
c            ELAB = 'CODE TOTAL POW LOSS'
c         endif
c
c
      elseif (iselect.eq.17) then  
c
         ELAB = 'IradIrad '
c         if (itype.eq.0) then           
c            ELAB = 'CUSTOM IMP PLRP'
c         elseif (itype.eq.1.or.itype.eq.2) then 
c            ELAB = 'CODE CUSTOM IMP PLRP'
c         endif
c
c----------------------------------------------------------
c
c     FLUID CODE - Background Plasma Properties
c
c----------------------------------------------------------
c
      elseif (iselect.eq.18) then   

         if (istate.eq.1) then 
            ELAB = 'FCNeFCNe   '
         elseif(istate.eq.2) then 
            ELAB = 'FCTeFCTe   '
         elseif(istate.eq.3) then 
            ELAB = 'FCTiFCTi   '
         elseif(istate.eq.4) then 
            ELAB = 'FCVbFCVb   '
         elseif(istate.eq.5) then 
            ELAB = 'FCEfFCEf   '
c         elseif(istate.eq.6) then 
c            ELAB = 'FCMaFCMa   '
         endif
c
c----------------------------------------------------------
c
c     FLUID CODE - Impurity Density
c
c----------------------------------------------------------
c
      elseif (iselect.eq.19) then   

         write(ELAB,'(''N'',i3,''N'',i3,'' (FC)'')')
     >                  istate,istate
c
c
c----------------------------------------------------------
c
c     FLUID CODE - Impurity Temperature
c
c----------------------------------------------------------
c
      elseif (iselect.eq.20) then   

         write(ELAB,'(''T'',i3,''T'',i3,'' (FC)'')')
     >                  istate,istate
c
c
c----------------------------------------------------------
c
c     FLUID CODE - Impurity Velocity
c
c----------------------------------------------------------
c
      elseif (iselect.eq.21) then   

         write(ELAB,'(''V'',i3,''V'',i3,'' (FC)'')')
     >                  istate,istate


c  
c----------------------------------------------------------
c     DIVIMP - EMISSION WEIGHTED AVERAGE ION TEMPERATURE 
c----------------------------------------------------------
c
      elseif (iselect.eq.22) then  
c
         ELAB = 'AVTiAVTi  '
c
c
c----------------------------------------------------------
c
c     DIVIMP - Impurity Density Ratio to BG Plasma
c
c----------------------------------------------------------
c
      elseif (iselect.eq.23) then   

         write(ELAB,'(''R'',i3,''R'',i3,'' (N)'')')
     >                 istate,istate
c
c
c----------------------------------------------------------
c
c     DIVIMP - Impurity Temperature Ratio to BG Plasma 
c
c----------------------------------------------------------
c
      elseif (iselect.eq.24) then   

         write(ELAB,'(''R'',i3,''R'',i3,'' (T)'')')
     >                 istate,istate
c
c
c----------------------------------------------------------
c
c     DIVIMP - Impurity Velocity Ratio to BG Plasma 
c
c----------------------------------------------------------
c
      elseif (iselect.eq.25) then   

         write(ELAB,'(''R'',i3,''R'',i3,'' (V)'')')
     >                 istate,istate

c  
c----------------------------------------------------------
c     EIRENE HBETA
c----------------------------------------------------------
c
      elseif (iselect.eq.26) then 
c   
         if (istate.eq.6) then 
            ELAB = 'HB THB T (EIR)'
c            if (itype.eq.0) then           
c               ELAB = 'EIRENE TOTAL HBETA'
c            elseif (itype.eq.1.or.itype.eq.2) then 
c               ELAB = 'CODE EIRENE TOT HBETA'
c            endif
         else
            write(elab,'(a,i2,a,i2)') 'HB',istate,'HB',istate
c            if (itype.eq.0) then           
c               write(elab,'(a,i4)') 'EIRENE HBETA COMP=',istate
c            elseif (itype.eq.1.or.itype.eq.2) then 
c               write(elab,'(a,i4)') 'CODE EIRENE HBETA COMP=',istate
c            endif
         endif
c
c----------------------------------------------------------
c
c     HC Related quantities
c
c----------------------------------------------------------
c
      elseif (iselect.eq.28.or.iselect.eq.29.or.iselect.eq.30.or.
     >        iselect.eq.33.or.iselect.eq.35) then 
c
c         call hc_set_elab(iselect,istate,itype,nizs,elab)
c
c
c----------------------------------------------------------
c
c     DIVIMP - Impurity Ionization
c
c----------------------------------------------------------
c
      elseif (iselect.eq.31) then   

         write(ELAB,'(''I'',i3,''I'',i3)')
     >                 istate,istate

      elseif (iselect.eq.36) then   
c
c                   1 = PINION = PIN ionization    
c                   2 = PINATOM = PIN Atom density 
c                   3 = PINMOL = PIN Molecular density
c                   4 = PINIONZ = Impurity ionization
c                   5 = PINZ0 = Impurity neutral density  
c                   6 = PINQI = Ion heating term
c                   7 = PINQE = Electron heating term

         if (istate.eq.1) then 
            ELAB = 'PIZ PIZ'
         elseif (istate.eq.2) then 
            ELAB = 'PAT PAT'
         elseif (istate.eq.3) then 
            ELAB = 'PMOL PMOL'
         elseif (istate.eq.4) then 
            ELAB = 'PZizPZiz '
         elseif (istate.eq.5) then 
            ELAB = 'PZniPZni'
         elseif (istate.eq.6) then 
            ELAB = 'PQI PQi'
         elseif (istate.eq.7) then 
            ELAB = 'PQe PQe'
         endif
c
c         len = lenstr(elab)
c
c         if (itype.eq.0) then 
c            elab = elab(1:len) // '/M^3)'
c         elseif (itype.eq.1) then 
c            elab = elab(1:len) // '/M^2)'
c         endif

      elseif (iselect.eq.40) then 
c
c         ExB drift related quantities
c         1 - Potential (phi) (V) 
c         2 - Radial Efield (V/m)
c         3 - Poloidal Efield (V/m)
c         4 - Radial ExB drift (m/s)
c         5 - Poloidal ExB drift (m/s)
c         6 - Radial ExB flux   ne x Vexb_rad (/m2/s)
c         7 - Poloidal ExB flux ne x Vexb_pol (/m2/s)
c     
         if (istate.eq.1) then 
            ELAB = 'EpotEpot'
         elseif (istate.eq.2) then 
            ELAB = 'EradErad'
         elseif (istate.eq.3) then 
            ELAB = 'EpolEpol'
         elseif (istate.eq.4) then 
            ELAB = 'EVR EVR'
         elseif (istate.eq.5) then 
            ELAB = 'EVP EVP'
         elseif (istate.eq.6) then 
            ELAB = 'EFR EFR'
         elseif (istate.eq.7) then 
            ELAB = 'EFP EFP'
         endif

      elseif (iselect.eq.41) then   

         write(ELAB,'(''W'',i3,''W'',i3)')
     >                 istate,istate

      elseif (iselect.eq.42) then   

         if (istate.eq.1) then 
            ELAB = 'IcndIcnd'
         elseif(istate.eq.2) then 
            ELAB = 'IcnvIcnv'
         elseif(istate.eq.3) then 
            ELAB = 'EcndEcnd'
         elseif(istate.eq.4) then 
            ELAB = 'EcnvEcnv'
         elseif(istate.eq.5) then 
            ELAB = 'CondCond'
         elseif(istate.eq.6) then 
            ELAB = 'ConvConv'
         elseif(istate.eq.7) then 
            ELAB = 'RVC RVC '
         elseif(istate.eq.8) then 
            ELAB = 'RVCeRVCe'
         endif
         
      endif
c
c
c
      return
      end

c
c     
c
      real function power_term(ik,ir,in)
      use mod_params
      use mod_cgeom
      use mod_dynam2
      use mod_dynam3
      use mod_comtor
      implicit none
c     include 'params' 
c     include 'cgeom'
c     include 'dynam2'
c     include 'dynam3'
c     include 'comtor'
c
      integer ik,ir,in
c
c     Calculate the requested power term for the specific cell
c       
c        1 = Ion conduction
c        2 = Ion convection
c        3 = Electron conduction
c        4 = Electron convection
c        5 = Total Conduction
c        6 = Total Convection
c        7 = Total Convection/total conduction
c        8 = Total Convection/electron conduction
c
c        NOTE: KFEGS and KFIGS are the gradient forces which are stored and loaded in OUT
c        FACT = QTIM * QTIM * EMI/CRMI
c        dTe/ds = KFEGS/FACT      
c
      real :: fact, conde, condi, conve, convi
c     
      fact = qtim * qtim * emi /crmi

      if (ik.eq.0) then
      ! values at first target idds(ir,2)

         condi = -CK0i*KTIDS(IDDS(IR,2))**2.5* KFIDS(idds(ir,2))/fact  
         convi =  2.5*KNDS(IDDS(IR,2))*KVDS(IDDS(IR,2))
     >                  *ECH*KTIDS(IDDS(IR,2)) +
     >         0.5*CRMB*AMU*(KVDS(IDDS(IR,2)))**3*knds(idds(ir,2))
      
         conde = -CK0*KTEDS(IDDS(IR,2))**2.5* KFEDS(idds(ir,2))/fact  
         conve =  2.5*KNDS(IDDS(IR,2))*KVDS(IDDS(IR,2))
     >                  *ECH*KTEDS(IDDS(IR,2)) 


      elseif (ik.eq.nks(ir)+1) then   
      ! values at second target idds(ir,1)

         condi = -CK0i*KTIDS(IDDS(IR,1))**2.5* KFIDS(idds(ir,1))/fact  
         convi =  2.5*KNDS(IDDS(IR,1))*KVDS(IDDS(IR,1))
     >                  *ECH*KTIDS(IDDS(IR,1)) +
     >         0.5*CRMB*AMU*(KVDS(IDDS(IR,1)))**3*knds(idds(ir,1))
      
         conde = -CK0*KTEDS(IDDS(IR,1))**2.5* KFEDS(idds(ir,1))/fact  
         conve =  2.5*KNDS(IDDS(IR,1))*KVDS(IDDS(IR,1))
     >                  *ECH*KTEDS(IDDS(IR,1)) 

      else   

         condi = -CK0i*KTIBS(IK,IR)**2.5* KFIGS(ik,ir)/fact  
         convi =  2.5*KNBS(IK,IR)*KVHS(IK,IR)/QTIM
     >                  *ECH*KTIBS(IK,IR) +
     >            0.5*CRMB*AMU*(KVHS(IK,IR)/QTIM)**3*knbs(ik,ir)
      
         conde = -CK0*KTEBS(IK,IR)**2.5* KFEGS(ik,ir)/fact  
         conve =  2.5*KNBS(IK,IR)*KVHS(IK,IR)/QTIM
     >                  *ECH*KTEBS(IK,IR) 

         write(6,'(a,2i6,20(1x,g12.5))') 'POW:',ik,ir,ck0,ck0i,
     >        fact,condi,conde,conde/condi,ktebs(ik,ir),ktibs(ik,ir),
     >        kfegs(ik,ir)/fact,kfigs(ik,ir)/fact,
     >        kfegs(ik,ir)/kfigs(ik,ir),ck0/ck0i,
     >        ktebs(ik,ir)/ktibs(ik,ir),
     >        kfegs(ik,ir)/kfigs(ik,ir)*ck0/ck0i*
     >        (ktebs(ik,ir)/ktibs(ik,ir))**2.5
         
      endif   
c      
      if (in.eq.1) then 
c        1 = Ion conduction
c            -k0i 5/2 Ti dTi/ds
         power_term = condi
       
      elseif (in.eq.2) then 
c        2 = Ion convection
c             5/2 nv kTi + 1/2 mi v^2 nv
          power_term = convi

      elseif (in.eq.3) then 
c        3 = Electron conduction
c         -k0e 5/2 Te dTe/ds
         power_term = conde

      elseif (in.eq.4) then 
c        4 = Electron convection
c            5/2 nv kTe
          power_term = conve

      elseif (in.eq.5) then 
c        5 = Total Conduction
          power_term = conde+condi

      elseif (in.eq.6) then 
c        6 = Total Convection
          power_term = conve+convi

       elseif (in.eq.7) then 
c        7 = Total Convection/total conduction
          if ((conde+condi).ne.0.0) then
             power_term = (conve+convi)/(conde+condi)
          else
             power_term = 0.0
          endif
       elseif (in.eq.8) then 
c        8 = Total Convection/electron conduction
          if (conde.ne.0.0) then
             power_term = (conve+convi)/conde
          else
             power_term = 0.0
          endif
      endif

      return
      end

      
