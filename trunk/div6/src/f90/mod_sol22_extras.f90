module mod_sol22_extras

! This module contains code that doesn't appear to be currently called from SOL22
  

  implicit none





contains


    SUBROUTINE NoName01(imflag, i,  loopstart,spts,npts,serr,exitcond,  note,te,ti,ne,ne2,ga,vb,vb2,pir,pii,vgradn,&
           act_press,exp_press,vsound,vsubs,vsupers,peiv,peiv2,pradv,pcxv,phelpiv,scxv,convi,conve,condi,conde,errcode,*)
! ======================================================================

! subroutine: NoName01

      use mod_solparams
      use mod_solcommon
      use mod_params
      use mod_slcom
!     INCLUDE 'solparams'
!     INCLUDE 'solcommon'
!     INCLUDE 'params'
!     INCLUDE 'slcom'
      IMPLICIT none
      integer i,i1,k,flag,imflag,lastflag,vcount
      integer loopstart
      integer errcode,npts,ir,irsep
      integer exitcond,negerrflag
      real*8 serr
      REAL*8  spts(mxspts)
      real*8  te(mxspts),ti(mxspts),ne(mxspts),vb(mxspts),exp_press(mxspts),act_press(mxspts),prad(mxspts)
      real*8    ga(mxspts),vb2(mxspts)
      real*8    ne2(mxspts),vsupers(mxspts),vsubs(mxspts)
      real*8    vsound(mxspts),scxv(mxspts)
      real*8    peiv(mxspts),pcxv(mxspts),phelpiv(mxspts)
      real*8    peiv2(mxspts),pradv(mxspts)
      real*8    pir(mxspts),pii(mxspts),vgradn(mxspts)
      real*8    condi(mxspts),conde(mxspts)
      real*8    convi(mxspts),conve(mxspts)
      CHARACTER*2 note(MXSPTS)
      IF (imflag.EQ.1) THEN
        IF (simag1.EQ.HI) THEN
          simag1 = simag - soffset
          simag2 = simag - soffset
        ELSE
          simag2 = simag - soffset
        ENDIF
      ENDIF
!        WRITE(PINOUT,*)
!        WRITE(PINOUT,'(I4,3F10.4,2I4)')
!          i,spts(i),osm_range,osm_range*spts(npts),imflag,exitcond
      IF ((i.GT.loopstart.AND.((imflag.EQ.1.AND.spts(i).GT.osm_range*spts(npts)).OR.exitcond.GE.3)).OR.exitcond.EQ.99) THEN
        ierror = i
        DO i1 = i, npts
          IF (osm_mode.EQ.2) note(i1) = ' e'
          te (i1) = te (i-1)
          ti (i1) = ti (i-1)
          ne (i1) = ne (i-1)
          ne2(i1) = ne2(i-1)
          ga (i1) = ga (i-1)
          vb (i1) = vb (i-1)
          vb2(i1) = vb2(i-1)
          pir   (i1) = pir   (i-1)
          pii   (i1) = pii   (i-1)
          vgradn(i1) = vgradn(i-1)
          act_press(i1) = act_press(i-1)
          exp_press(i1) = exp_press(i-1)
          vsound (i1) = vsound (i-1)
          vsubs  (i1) = vsubs  (i-1)
          vsupers(i1) = vsupers(i-1)
          peiv   (i1) = peiv   (i-1)
          peiv2  (i1) = peiv2  (i-1)
          prad   (i1) = prad   (i-1)
          pradv  (i1) = pradv  (i-1)
          pcxv   (i1) = pcxv   (i-1)
          phelpiv(i1) = phelpiv(i-1)
          scxv   (i1) = scxv   (i-1)
          convi(i1) = convi(i-1)
          conve(i1) = conve(i-1)
          condi(i1) = condi(i-1)
          conde(i1) = conde(i-1)
        ENDDO
        IF (spts(i).GT.osm_range*spts(npts)) THEN
          IF ((exitcond.GE.3.AND.exitcond.LT.99).AND.snegerr.LT.simag) THEN
            errcode = exitcond
            serr    = snegerr
          ELSEIF (imflag.GT.0) THEN
            errcode = imflag
            serr    = simag
          ENDIF
          errcode = -errcode
!          GOTO 2000
          RETURN 1
        ENDIF
      ENDIF
      RETURN
99    STOP



    END SUBROUTINE NoName01


end module mod_sol22_extras
