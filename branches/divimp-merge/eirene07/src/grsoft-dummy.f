C EIRENE07 COMPILATION
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
C ===== SOURCE: grclp.f
c------------------------------------------------------------------------
      SUBROUTINE GRCLP(K)
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
      SUBROUTINE GRDRLG(PA,TEXTX,TEXTY,TEXTZ,IOPT,a,b,c,i1,i2,i3)
      Real, intent(in) :: pa(1)
      character(1), intent(in) :: textx,texty,textz
      integer, intent(in) :: iopt
      real,intent(in) :: a,b,c
      integer,intent(in) :: i1,i2,i3
      RETURN
      END
C ===== SOURCE: grdrne.f
c-------------------------------------------------------------------------
      SUBROUTINE GRDRNE(PA,NROW,XYZ)
      REAL, INTENT(IN) :: PA(1), XYZ(1,1,1)
      INTEGER, INTENT(IN) :: NROW
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
      REAL XX(N), YY(N)
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
C ===== SOURCE: setgks.f
c------------------------------------------------------------------------
      SUBROUTINE SETGKS(IVECT,IERR)
      RETURN
      END
