C EIRENE07 COMPILATION
C ===== SOURCE: chctrc.f
c------------------------------------------------------------------------
      SUBROUTINE CHCTRC(XPLO,YPLO,ZPLO,IFLAG,ISYM)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CLOGAU
      USE CTRCEI
      USE COMPRT
      USE COMSOU
      USE CCONA
      USE CUPD

      IMPLICIT NONE
C
      INTEGER,PARAMETER :: NTXHST=18

      REAL(DP), INTENT(IN) :: XPLO, YPLO, ZPLO
      INTEGER, INTENT(IN) :: IFLAG, ISYM
      INTEGER :: IGENCS, IGENCV, NLLI
      CHARACTER(20) :: TXTHST(NTXHST)

      DATA TXTHST
     .           /'LOCATE(1)           ',
     .            'ELECTR. IMPACT(2)   ',
     .            'HEAVY PAR. IMPACT(3)',
     .            'PHOTON IMPACT(4)    ',
     .            'ELASTIC COLL.(5)    ',
     .            'CHARGE EXCHANGE(6)  ',
     .            'FOKKER PLANCK(7)    ',
     .            'SURFACE(8)          ',
     .            'SPLITTING(9)        ',
     .            'RUSSIAN ROULETTE(10)',
     .            'PERIODICITY(11)     ',
     .            'RESTART:A. SPLT.(12)',
     .            'SAVE:COND. EXP.(13) ',
     .            'RESTART:COND EXP(14)',
     .            'TIME LIMIT(15)      ',
     .            'GENERATION LIMIT(16)',
     .            'FLUID LIMIT(17)     ',
     .            'ERROR DETECTED      '/
C
C  WRITE TRACK DATA
C
c slmod begin
      IF (.TRUE.) THEN
        CALL USRTRC(XPLO,YPLO,ZPLO,IFLAG,ISYM)
      ENDIF 
c slmod end 
      IF (TRCHST) THEN
        CALL LEER(1)
        WRITE (iunout,*) TXTHST(ISYM)
        IF (ISYM.EQ.1)  CALL MASJ1('NPANU   ',NPANU)
        WRITE (iunout,'(1X,A8)') TEXTS(ISPZ)
        CALL MASJ4 ('ITIME,IFPATH,IUPDTE,ICOL        ',
     .               ITIME,IFPATH,IUPDTE,ICOL)
        CALL MASR3 ('X0,Y0,Z0                ',XPLO,YPLO,ZPLO)
        CALL MASR6 ('VELX,VELY,VELZ,VEL,E0,WEIGHT                    ',
     .               VELX,VELY,VELZ,VEL,E0,WEIGHT)
        CALL MASR1 ('TIME    ',TIME)
        CALL MASJ2 ('IGENCV,IGENCS   ',IGENCV,IGENCS)
        CALL MASJ4 ('NRCELL,IPOLG,NACELL,NBLOCK      ',
     .               NRCELL,IPOLG,NACELL,NBLOCK)
        IF (NLTOR) THEN
          CALL MASJ1 ('NTCELL  ',NTCELL)
        ENDIF
        IF (NLTRA) THEN
          CALL MASJ1R ('NNTCLL,PHI      ',NNTCLL,PHI/DEGRAD)
        ENDIF
        IF (NLPOL) THEN
          CALL MASJ1 ('NPCELL  ',NPCELL)
        ENDIF
        IF ((ISYM.GE.6.AND.ISYM.LE.10).OR.
     .      (ISYM.EQ.1.AND.NLSRF(ISTRA))) THEN
          CALL MASJ5 ('MRSURF,MPSURF,MTSURF,MASURF,NLLI        ',
     .                 MRSURF,MPSURF,MTSURF,MASURF,NLLI)
          CALL MASR1 ('SCOS    ',SCOS)
        ENDIF
      ENDIF

      RETURN
      END
C ===== SOURCE: ellipsoid.f



      SUBROUTINE ELLIPSOID (X0,Y0,Z0,CX,CY,CZ,XLIMS1,YLIMS1,ZLIMS1,
     .                      XLIMS2,YLIMS2,ZLIMS2,RLB,ILCOL,NX,NY,NZ)

      USE PRECISION
      USE CCONA

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: X0, Y0, Z0, CX, CY, CZ, XLIMS1, YLIMS1,
     .                      ZLIMS1, XLIMS2, YLIMS2, ZLIMS2, RLB
      INTEGER, INTENT(IN) :: ILCOL, NX, NY, NZ

      RETURN
      END
C ===== SOURCE: gr3axs.f
c------------------------------------------------------------------------
      subroutine gr3axs(ar,ier,ext,valu,chaxs,ll,i1,i2)
      USE PRECISION
      integer, intent(in) :: ier, i1, i2
      real(sp) :: ar(*), valu(3,2), ext(3,3)
      character(20) :: chaxs(3)
      logical ll
      return
      end
C ===== SOURCE: gr3dim.f
c------------------------------------------------------------------------
      subroutine gr3dim(lar,ier)
      return
      end
C ===== SOURCE: gr3ext.f
c------------------------------------------------------------------------
      subroutine gr3ext(ar,ier,ext)
      return
      end
C ===== SOURCE: gr3net.f
c------------------------------------------------------------------------
      subroutine gr3net(ar,ier,i1,xyz,i2,i3,i4,i5,i6,i7)
      return
      end
C ===== SOURCE: gr3nt1.f
c------------------------------------------------------------------------
      subroutine gr3nt1(ar,ier,i1,x,y,z,i2,i3,i4,i5,i6,i7)
      return
      end
C ===== SOURCE: gr3plo.f
c------------------------------------------------------------------------
      subroutine gr3plo(ar,ier,cha)
      character*(*) cha
      return
      end
C ===== SOURCE: gr3rot.f
c------------------------------------------------------------------------
      subroutine gr3rot(ar,ier,ch1,x1,ch2,x2,ch3,x3)
      character*(*) ch1,ch2,ch3
      return
      end
C ===== SOURCE: grarrw.f
c------------------------------------------------------------------------
      SUBROUTINE GRARRW(XP,YP,XTIP,YTIP,ALEN,AWID,ICODE)
      RETURN
      END
C ===== SOURCE: graxs.f
c-------------------------------------------------------------------------
      SUBROUTINE GRAXS(LOPT,OPTION,LTXTX,TEXTX,LTXTXY,TEXTXY)
      character*(*) option,textx,textxy
      RETURN
      END
C ===== SOURCE: grbld.f
c-------------------------------------------------------------------------
      SUBROUTINE GRBLD(XCM,YCM,ISK,JSK,XMIN,XMAX,YMIN,YMAX,NKURV)
      RETURN
      END
C ===== SOURCE: grchn.f
c-------------------------------------------------------------------------
      SUBROUTINE GRCHN(XX,YY,M,NR)
      RETURN
      END
C ===== SOURCE: grchrc.f
c-------------------------------------------------------------------------
      SUBROUTINE GRCHRC(HEIGHT,ANGLE,INTS)
      RETURN
      END
C ===== SOURCE: grdrdm.f
c-------------------------------------------------------------------------
      SUBROUTINE GRDRDM(PA,NROW,TAB,X,Y)
      RETURN
      END
C ===== SOURCE: grdrhs.f
c-------------------------------------------------------------------------
      SUBROUTINE GRDRHS(PA,NPKNT,PNKT,X,Y)
      RETURN
      END
C ===== SOURCE: grdrlg.f
c-------------------------------------------------------------------------
      SUBROUTINE GRDRLG(PA,TEXTX,TEXTY,TEXTZ,IOPT)
      RETURN
      END
C ===== SOURCE: grdrne.f
c-------------------------------------------------------------------------
      SUBROUTINE GRDRNE(PA,NROW,XYZ)
      RETURN
      END
C ===== SOURCE: grdrw.f
c-------------------------------------------------------------------------
      SUBROUTINE GRDRW(X,Y)
      RETURN
      END
C ===== SOURCE: grdsh.f
c-------------------------------------------------------------------------
      SUBROUTINE GRDSH(A1,A2,A3)
      RETURN
      END
C ===== SOURCE: grend.f
c-------------------------------------------------------------------------
      SUBROUTINE GREND
      RETURN
      END
C ===== SOURCE: grfill.f
c-------------------------------------------------------------------------
      SUBROUTINE GRFILL(N,XX,YY,ISTYLE,ITYPE)
      RETURN
      END
C ===== SOURCE: grfrbn.f
c-------------------------------------------------------------------------
      SUBROUTINE GRFRBN(IFU,IKA,IKS,IRA,IRI)
      RETURN
      END
C ===== SOURCE: grftoc.f
c-------------------------------------------------------------------------
      SUBROUTINE GRFTOC(F,C,L)
      RETURN
      END
C ===== SOURCE: grhhnl.f
c------------------------------------------------------------------------
      SUBROUTINE GRHHNL
      RETURN
      END
C ===== SOURCE: grjmp.f
c-------------------------------------------------------------------------
      SUBROUTINE GRJMP(X,Y)
      RETURN
      END
C ===== SOURCE: grjmps.f
c-------------------------------------------------------------------------
      SUBROUTINE GRJMPS(X,Y,NR)
      RETURN
      END
C ===== SOURCE: grln.f
c-------------------------------------------------------------------------
      SUBROUTINE GRLN(XX,YY,M)
      RETURN
      END
C ===== SOURCE: grmrks.f
c------------------------------------------------------------------------
      SUBROUTINE GRMRKS(X)
      RETURN
      END
C ===== SOURCE: grnwpn.f
c------------------------------------------------------------------------
      SUBROUTINE GRNWPN(I)
      RETURN
      END
C ===== SOURCE: grnxtb.f
c------------------------------------------------------------------------
      SUBROUTINE GRNXTB(K)
      RETURN
      END
C ===== SOURCE: grnxtf.f
c------------------------------------------------------------------------
      SUBROUTINE GRNXTF
      RETURN
      END
C ===== SOURCE: grsclc.f
c------------------------------------------------------------------------
      SUBROUTINE GRSCLC(XA,YA,XB,YB)
      RETURN
      END
C ===== SOURCE: grsclv.f
c------------------------------------------------------------------------
      SUBROUTINE GRSCLV(XA,YA,XB,YB)
      RETURN
      END
C ===== SOURCE: grspts.f
c------------------------------------------------------------------------
      SUBROUTINE GRSPTS(INT)
      RETURN
      END
C ===== SOURCE: grstrt.f
c------------------------------------------------------------------------
      SUBROUTINE GRSTRT(ICAMERA,IDDNUMB)
      RETURN
      END
C ===== SOURCE: grtxt.f
c------------------------------------------------------------------------
      SUBROUTINE GRTXT(X,Y,LTEXT,TEXT)
      character*(*) text
      RETURN
      END
C ===== SOURCE: grtxtc.f
c------------------------------------------------------------------------
      SUBROUTINE GRTXTC(LTEXT,TEXT)
      character*(*) text
      RETURN
      END
C ===== SOURCE: gstxal.f
c------------------------------------------------------------------------
      SUBROUTINE GSTXAL(IALH,IALV)
      RETURN
      END
C ===== SOURCE: kurvef.f
c------------------------------------------------------------------------
      SUBROUTINE KURVEF(X,Y,IST,ISY)
      RETURN
      END
C ===== SOURCE: outpat.f
c------------------------------------------------------------------------
      SUBROUTINE OUTPAT
      RETURN
      END
C ===== SOURCE: pl3d.f


C
C
      SUBROUTINE PL3D(PX,PY,PZ,PP1,PP2)

      USE PRECISION

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: PX, PY, PZ
      REAL(DP), INTENT(OUT) :: PP1, PP2

      PP1=0.d0
      PP2=0.d0
      return
      end
C ===== SOURCE: plt2d.f
      subroutine plt2d
      return
      end
C ===== SOURCE: plteir.f
c------------------------------------------------------------------------
      subroutine plteir(i)
      return
      end
C ===== SOURCE: rpsout.f
c------------------------------------------------------------------------
      SUBROUTINE RPSOUT
      RETURN
      END
C ===== SOURCE: setgks.f
c------------------------------------------------------------------------
      SUBROUTINE SETGKS(IVECT,IERR)
      RETURN
      END
c slmod begin
C ===== SOURCE: exit_own.f
C
C
      SUBROUTINE EXIT_OWN (ICC)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ICC
      CALL GREND
      STOP
      END SUBROUTINE EXIT_OWN
C ===== SOURCE: plttly.f
C  10.6.05:  L_SAME:  USE SAME FRAME AS IN PREVIOUS CALL
C  8.8.06 :  GRPP taken out
C
      SUBROUTINE PLTTLY (X,Y,VBAR,YMN,YMX,IR1,IR2,IRS,NKURV,TXTTAL,
     .                   TXTSPC,TXTUNT,TXTRUN,TXHEAD,
     .                   LBAR,XMI,XMA,YMNLG,YMXLG,LPLOT,LHIST,IERR,
     .                   N1BAR,N1DIM,L_SAME)

      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CPLMSK

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: N1BAR, N1DIM
      REAL(DP), INTENT(IN) :: X(*), VBAR(N1BAR,*),
     .                      YMN(*), YMX(*), YMNLG(*), YMXLG(*)
      REAL(DP), INTENT(INOUT) :: Y(N1DIM,*)
      INTEGER, INTENT(IN) :: IR1(*), IR2(*), IRS(*), NKURV
      LOGICAL, INTENT(IN) :: LBAR(*), LPLOT(*), L_SAME
      LOGICAL, INTENT(IN) :: LHIST
      CHARACTER(LEN=*), INTENT(IN) :: TXTTAL(*),TXTSPC(*),TXTUNT(*),
     .                                TXTRUN, TXHEAD

      REAL(DP) :: YA, YMINY, YMY, AA, FM, ST1, ST2, FP, DMINY, DMAXY,
     .          XMI, XMA
      REAL(DP) :: XMIN, XMAX, YMIN, YMAX
      REAL(SP), SAVE :: PRMSAVE(8)
      INTEGER :: I1, I2, IS, ICURV, IT, ISY, NP, J, NPS, IKURV,
     .           I, IERR, IPEN1, IPEN2
      CHARACTER(10) :: CHR
      CHARACTER(12) :: CHR12
      SAVE YA,IPEN1,IPEN2

      END
C ===== SOURCE: schnitp.f


      SUBROUTINE SCHNITP(X1,Y1,X2,Y2,X3,Y3,X4,Y4,EX,EY)
C   INTERSECTION POINT E OF 2 STRAIGHT LINES G1 AND G2
C   G1 IS DEFINED BY POINTS 1 AND 2
C   G2 IS DEFINED BY POINTS 3 AND 4

      USE PRECISION

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: X1, Y1, X2, Y2, X3, Y3, X4, Y4
      REAL(DP), INTENT(OUT) :: EX, EY
      REAL(DP) :: MUE

      MUE = ((Y2-Y4)*(X3-X4)+(X4-X2)*(Y3-Y4))/
     .      ((X1-X2)*(Y3-Y4)-(Y1-Y2)*(X3-X4)+1.D-20)
      EX = X2 + MUE * (X1-X2)
      EY = Y2 + MUE * (Y1-Y2)
      END
c
c
c
      SUBROUTINE PLTEIR_REINIT
      END
      SUBROUTINE PL3D_REINIT
      END
      SUBROUTINE PLT2D_REINIT
      END
      SUBROUTINE STCOOR_REINIT
      END
C
c slmod end
