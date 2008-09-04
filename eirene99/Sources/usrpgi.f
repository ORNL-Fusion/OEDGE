c
c     ==================================================================
c     SYSTEM CALLS
c     ==================================================================
c
      SUBROUTINE Clock(s)
      IMPLICIT none
      DOUBLE PRECISION s
      REAL etime
      REAL rdum1,val(2)
      rdum1 = etime(val)
      s = DBLE(val(1))
      RETURN
      END

      SUBROUTINE Clock_1970(s)
      IMPLICIT none
      DOUBLE PRECISION s
      INTEGER time
      s = DBLE(time())
      RETURN
      END
c
c     ==================================================================
c     SUBSITUTE SUBROUTINES
c     ==================================================================
c
      SUBROUTINE TABPRC(i1,i2,i3)
      IMPLICIT none
      INTEGER i1,i2,i3
      WRITE(0,*) 'WARNING: Calling TABPRC routine substitute'
      RETURN
      END

      SUBROUTINE GRNXTB(i1)
      IMPLICIT none
      INTEGER i1
      WRITE(0,*) 'WARNING: Calling GRNXTB routine substitute'
      RETURN
      END

      SUBROUTINE KURVEF(d1,d2,i1,i2)
      IMPLICIT none
      INTEGER          i1,i2
      DOUBLE PRECISION d1,d2
      WRITE(0,*) 'WARNING: Calling KURVEF routine substitute'
      RETURN
      END

      SUBROUTINE GRBLD(r1,r2,i1,i2,r3,r4,r5,r6)
      IMPLICIT none
      INTEGER i1,i2
      REAL    r1,r2,r3,r4,r5,r6
      WRITE(0,*) 'WARNING: Calling GRBLD  routine substitute'
      RETURN
      END

      SUBROUTINE SM0USR(i1,i2,r1,r2,r3,r4,r5,r6)
      IMPLICIT none
      INTEGER          i1,i2
      DOUBLE PRECISION r1,r2,r3,r4,r5,r6
      WRITE(0,*) 'WARNING: Calling SM0USR routine substitute'
      RETURN
      END

      SUBROUTINE SP1USR
      IMPLICIT none
      WRITE(0,*) 'WARNING: Calling SP1USR routine substitute'
      RETURN
      END

c      SUBROUTINE SM1USR(i1,d1,d2,d3 ,d4 ,d5 ,d6 ,d7 ,d8 ,d9 ,i2,i3,
c     .                  i4,i5,i6,d10,d11,d12,d13,d14,d15,d16)
c      IMPLICIT none
c      INTEGER          i1 ,i2 ,i3, i4,i5,i6
c      DOUBLE PRECISION d1 ,d2 ,d3 ,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,
c     .                 d14,d15,d16
c      WRITE(0,*) 'WARNING: Calling SM1USR routine substitute'
c      RETURN
c      END

      SUBROUTINE GRSCLV(r1,r2,r3,r4)
      IMPLICIT none
      REAL       r1,r2,r3,r4
      WRITE(0,*) 'WARNING: Calling GRSCLV routine substitute'
      RETURN
      END

      SUBROUTINE GRSCLC(r1,r2,r3,r4)
      IMPLICIT none
      REAL       r1,r2,r3,r4
      WRITE(0,*) 'WARNING: Calling GRSCLC routine substitute'
      RETURN
      END

      SUBROUTINE GRJMP(r1,r2)
      IMPLICIT none
      REAL       r1,r2
      WRITE(0,*) 'WARNING: Calling GRJMP  routine substitute'
      RETURN
      END

      SUBROUTINE GRNXTF
      IMPLICIT none
      WRITE(0,*) 'WARNING: Calling GRNXTF routine substitute'
      RETURN
      END

      SUBROUTINE GRNWPN(i1)
      IMPLICIT none
      INTEGER    i1
      WRITE(0,*) 'WARNING: Calling GRNWPN routine substitute'
      RETURN
      END

      SUBROUTINE GRDRW(r1,r2)
      IMPLICIT none
      REAL r1,r2
      WRITE(0,*) 'WARNING: Calling GRDRW  routine substitute'
      RETURN
      END

      SUBROUTINE EXIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CALL GREND
      STOP
      END

      SUBROUTINE CHCTRC(d1,d2,d3,i1,i2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DOUBLE PRECISION d1,d2,d3
      INTEGER          i1,i2,count
      
      DATA count /0/
     
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'CUPD'
      INCLUDE 'CGRID'

      IF (i1.EQ.0) count = count + 1

c...  Old output:
c      WRITE(80,'(I8,3F14.6,4I8,1P,1E10.2,0P,F8.4,7I5)') 
c     .  count,d1,d2,d3,i1,i2,istra,npanu,e0,weight,nacell,ifpath,iupdte,
c     .  ityp,nblock,masurf,msurf

      WRITE(80,'(I8,3F9.3,3I4,I4,I7,1P,1E10.2,0P,F8.4,I6,7I5,2X,5I5,
     .  1P,E10.2,0P)') 
     .  count,d1,d2,d3,i1,i2,istra,ntrseg,
c     .  MAX(0.0,REAL(ntrseg)-1.0)*360.0/REAL(NTTRA-1),
c     .  atan(REAL(d3/d1))*180.0/3.1415+
c     .       MAX(0.0,REAL(ntrseg)-1.0)*360.0/REAL(NTTRA-1),
     .  npanu,e0,weight,nacell,ifpath,iupdte,
     .  ityp,nblock,masurf,msurf,mtsurf,
     .  nrcell,npcell,nblock,nntcll,ntcell,
     .  time

      
c      WRITE(0,*) 'WARNING: Calling CHCTRC routine substitute'
      RETURN
      END

      SUBROUTINE GRSTRT(i1,i2)
      IMPLICIT   none
      INTEGER i1,i2
c      WRITE(0,*) 'WARNING: Calling GRSTRT routine substitute'
      RETURN
      END

      SUBROUTINE GREND
c      STOP 'GREND'
      STOP
      RETURN
      END

      SUBROUTINE OUTPAT
      WRITE(0,*) 'WARNING: Calling OUTPAT routine substitute'
      RETURN
      END

      SUBROUTINE PLT2D
c      WRITE(0,*) 'WARNING: Calling PLT2D  routine substitute'
      RETURN
      END

      SUBROUTINE RPSOUT
      WRITE(0,*) 'WARNING: Calling RPSOUT routine substitute'
      RETURN
      END

      SUBROUTINE PLTEIR(i1)
      IMPLICIT   none
      INTEGER i1
      WRITE(0,*) 'WARNING: Calling PLTEIR routine substitute'
      RETURN
      END

      DOUBLE PRECISION FUNCTION X05BAF()
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMUSR'
      DOUBLE PRECISION s
      CALL Clock(s)
c      CALL Clock_1970(s)
c      CALL SECONDS_SINCE_1980@(s)
      x05baf = s
      RETURN
      END

      SUBROUTINE F04ARF(d1,i1,d2,i2,d3,d4,i3)
      IMPLICIT none
      INTEGER i1,i2,i3
      DOUBLE PRECISION d1(i1,*),d2(*),d3(*),d4(*)
      WRITE(0,*) 'WARNING: Calling F04ARF routine substitute'
      RETURN
      END

      DOUBLE PRECISION FUNCTION S15AEF(d1,i1)
      IMPLICIT none
      INTEGER          i1
      DOUBLE PRECISION d1
      WRITE(0,*) 'WARNING: Calling S15AEF routine substitute'
      d1 = 0.0D0
      i1 = 0
      S15AEF = 0.0
      RETURN
      END

      SUBROUTINE F04ATF(d1,i1,d2,i2,d3,d4,i3,d5,d6,i5)
      IMPLICIT none
      INTEGER          IAA,i1,i2,i3,i4,i5
      PARAMETER       (IAA=50)
      DOUBLE PRECISION d1(i1,*),d2(*),d3(*),d4(IAA,IAA),d5(IAA),d6(IAA)
      WRITE(0,*) 'WARNING: Calling F04ATF routine substitute'
      RETURN
      END

      DOUBLE PRECISION FUNCTION S13AAF(d1,i1)
      IMPLICIT none
      INTEGER          i1
      DOUBLE PRECISION d1
      WRITE(0,*) 'WARNING: Calling S13AAF routine substitute'
      d1 = 1.0
      i1 = 0
      S13AAF = 1.0D+14
      RETURN
      END

      DOUBLE PRECISION FUNCTION X02AJF()
      IMPLICIT NONE
      WRITE(0,*) 'WARNING: Calling X02AJF routine substitute'
      X02AJF = 1.0D+14
      RETURN
      END

      DOUBLE PRECISION FUNCTION X02AKF()
      IMPLICIT NONE
      WRITE(0,*) 'WARNING: Calling X02AKF routine substitute'
      X02AKF = 1.0D+14
      RETURN
      END

      DOUBLE PRECISION FUNCTION X02ALF()
      IMPLICIT NONE
      WRITE(0,*) 'WARNING: Calling F04ALF routine substitute'
      X02ALF = 1.0D+14
      RETURN
      END



