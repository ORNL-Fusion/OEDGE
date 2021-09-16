c------------------------------------------------------------------------
      SUBROUTINE CHCTRC(XPLO,YPLO,ZPLO,IFLAG,ISYM)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CLOGAU
      USE CTRCEI
      USE COMPRT
      USE COMSOU

      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: XPLO, YPLO, ZPLO 
      INTEGER, INTENT(IN) :: IFLAG, ISYM
      INTEGER :: IGENCS, IGENCV, NLLI, NNTCLL
      CHARACTER(20) :: TXTHST(15)

      DATA TXTHST
     .           /'LOCATE(1)           ',
     .            'ELECTR. IMPACT(2)   ',
     .            'ION IMPACT(3)       ',
     .            'CHARGE EXCHANGE(4)  ',
     .            'FOKKER PLANCK(5)    ',                                 
     .            'SURFACE(6)          ',                                 
     .            'SPLITTING(7)        ',
     .            'RUSSIAN ROULETTE(8) ',
     .            'PERIODICITY(9)      ',
     .            'RESTART:A. SPLT.(10)',
     .            'SAVE:COND. EXP.(11) ',
     .            'RESTART:COND EXP(12)',
     .            'TIME LIMIT(13)      ',
     .            'GENERATION LIMIT(14)',
     .            'ERROR DETECTED      '/
c slmod begin (sl)
      INTEGER COUNT
      SAVE

      IF (iflag.EQ.0) count = count + 1

      WRITE(80,'(I8,3F9.3,5I4,I7,1P,1E10.2,0P,F8.4,I6,7I5,2X,5I5,
     .  1P,E10.2,0P)') 
     .  count,XPLO,YPLO,ZPLO,IFLAG,ISYM,istra,ntrseg,iperid,
     .  npanu,e0,weight,nacell,ifpath,iupdte,
     .  ityp,nblock,masurf,msurf,mtsurf,
     .  nrcell,npcell,nblock,nntcll,ntcell,
     .  time

      RETURN
c slmod end
C
C  WRITE TRACK DATA
C
      IF (TRCHST) THEN
        CALL LEER(1)
        WRITE (6,*) TXTHST(ISYM)
        IF (ISYM.EQ.1)  CALL MASJ1('NPANU   ',NPANU)
        WRITE (6,'(1X,A8)') TEXTS(ISPZ)
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
          CALL MASJ1R ('NNTCLL,PHI      ',NNTCLL,PHI)
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
c------------------------------------------------------------------------
      subroutine gr3axs(ar,ier,ext,valu,chaxs,ll,i1,i2)
      USE PRECISION
      integer, intent(in) :: ier, i1, i2
      real(sp) :: ar(*), valu(3,2), ext(3,3)
      character(20) :: chaxs(3)      
      logical ll
      return
      end
c------------------------------------------------------------------------
      subroutine gr3dim(lar,ier)
      return
      end
c------------------------------------------------------------------------
      subroutine gr3ext(ar,ier,ext)
      return
      end
c------------------------------------------------------------------------
      subroutine gr3net(ar,ier,i1,xyz,i2,i3,i4,i5,i6,i7)
      return
      end
c------------------------------------------------------------------------
      subroutine gr3nt1(ar,ier,i1,x,y,z,i2,i3,i4,i5,i6,i7)
      return
      end
c------------------------------------------------------------------------
      subroutine gr3plo(ar,ier,cha)
      character*(*) cha
      return
      end
c------------------------------------------------------------------------
      subroutine gr3rot(ar,ier,ch1,x1,ch2,x2,ch3,x3)
      character*(*) ch1,ch2,ch3
      return
      end
c------------------------------------------------------------------------
      SUBROUTINE GRARRW(XP,YP,XTIP,YTIP,ALEN,AWID,ICODE)
      RETURN
      END
c-------------------------------------------------------------------------
      SUBROUTINE GRAXS(LOPT,OPTION,LTXTX,TEXTX,LTXTXY,TEXTXY)
      character*(*) option,textx,textxy
      RETURN
      END
c-------------------------------------------------------------------------
      SUBROUTINE GRBLD(XCM,YCM,ISK,JSK,XMIN,XMAX,YMIN,YMAX,NKURV)
      RETURN
      END
c-------------------------------------------------------------------------
      SUBROUTINE GRCHN(XX,YY,M,NR)
      RETURN
      END
c-------------------------------------------------------------------------
      SUBROUTINE GRCHRC(HEIGHT,ANGLE,INTS)
      RETURN
      END
c-------------------------------------------------------------------------
      SUBROUTINE GRDRDM(PA,NROW,TAB,X,Y)
      RETURN
      END
c-------------------------------------------------------------------------
      SUBROUTINE GRDRHS(PA,NPKNT,PNKT,X,Y)
      RETURN
      END
c-------------------------------------------------------------------------
      SUBROUTINE GRDRLG(PA,TEXTX,TEXTY,TEXTZ,IOPT)
      RETURN
      END
c-------------------------------------------------------------------------
      SUBROUTINE GRDRNE(PA,NROW,XYZ)
      RETURN
      END
c-------------------------------------------------------------------------
      SUBROUTINE GRDRW(X,Y)
      RETURN
      END
c-------------------------------------------------------------------------
      SUBROUTINE GRDSH(A1,A2,A3)
      RETURN
      END
c-------------------------------------------------------------------------
      SUBROUTINE GREND
      RETURN
      END
c-------------------------------------------------------------------------
      SUBROUTINE GRFILL(N,XX,YY,ISTYLE,ITYPE)
      RETURN
      END
c-------------------------------------------------------------------------
      SUBROUTINE GRFRBN(IFU,IKA,IKS,IRA,IRI)
      RETURN
      END
c-------------------------------------------------------------------------
      SUBROUTINE GRFTOC(F,C,L)
      RETURN
      END
c------------------------------------------------------------------------
      SUBROUTINE GRHHNL
      RETURN
      END
c-------------------------------------------------------------------------
      SUBROUTINE GRJMP(X,Y)
      RETURN
      END
c-------------------------------------------------------------------------
      SUBROUTINE GRJMPS(X,Y,NR)
      RETURN
      END
c-------------------------------------------------------------------------
      SUBROUTINE GRLN(XX,YY,M)
      RETURN
      END
c------------------------------------------------------------------------
      SUBROUTINE GRMRKS(X)
      RETURN
      END
c------------------------------------------------------------------------
      SUBROUTINE GRNWPN(I)
      RETURN
      END
c------------------------------------------------------------------------
      SUBROUTINE GRNXTB(K)
      RETURN
      END
c------------------------------------------------------------------------
      SUBROUTINE GRNXTF
      RETURN
      END
c------------------------------------------------------------------------
      SUBROUTINE GRSCLC(XA,YA,XB,YB)
      RETURN
      END
c------------------------------------------------------------------------
      SUBROUTINE GRSCLV(XA,YA,XB,YB)
      RETURN
      END
c------------------------------------------------------------------------
      SUBROUTINE GRSPTS(INT)
      RETURN
      END
c------------------------------------------------------------------------
      SUBROUTINE GRSTRT(ICAMERA,IDDNUMB)
      RETURN
      END
c------------------------------------------------------------------------
      SUBROUTINE GRTXTC(LTEXT,TEXT)
      character*(*) text
      RETURN
      END
c------------------------------------------------------------------------
      SUBROUTINE GRTXT(X,Y,LTEXT,TEXT)
      character*(*) text
      RETURN
      END
c------------------------------------------------------------------------
      SUBROUTINE GSTXAL(IALH,IALV)
      RETURN
      END
c------------------------------------------------------------------------
      SUBROUTINE KURVEF(X,Y,IST,ISY)
      RETURN
      END
c------------------------------------------------------------------------
      SUBROUTINE OUTPAT
      RETURN
      END


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
      subroutine plt2d
      return
      end
c------------------------------------------------------------------------
      subroutine plteir(i)
      return
      end
c------------------------------------------------------------------------
      SUBROUTINE RPSOUT
      RETURN
      END
c------------------------------------------------------------------------
      SUBROUTINE SETGKS(IVECT,IERR)
      RETURN
      END
