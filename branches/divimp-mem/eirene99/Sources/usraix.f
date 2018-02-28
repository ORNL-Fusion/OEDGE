c
c     ==================================================================
c     SYSTEM CALLS
c     ==================================================================
c
c      SUBROUTINE Clock(s)
c      IMPLICIT none
c      DOUBLE PRECISION s
c      REAL etime
c      REAL rdum1,val(2)
c      rdum1 = etime(val)
c      s = DBLE(val(1))
c      RETURN
c      END
c
      SUBROUTINE Clock_1970(s)
      IMPLICIT none
      DOUBLE PRECISION s
      INTEGER time
c      s = DBLE(time())
      s = 0.0
      write(0,*) 'MOD: NO CLOCK_1970'
      RETURN
      END
c
c
c     ==================================================================
c     SUBSITUTE SUBROUTINES
c     ==================================================================
c

      SUBROUTINE SM0USR(i1,i2,r1,r2,r3,r4,r5,r6)
      IMPLICIT none
      INTEGER          i1,i2
      DOUBLE PRECISION r1,r2,r3,r4,r5,r6
      WRITE(0,*) 'WARNING: Calling SM0USR routine substitute'
      RETURN
      END
c
c
c
      SUBROUTINE SP1USR
      IMPLICIT none
      WRITE(0,*) 'WARNING: Calling SP1USR routine substitute'
      RETURN
      END
c
c
c
      SUBROUTINE SM1USR(i1,d1,d2,d3 ,d4 ,d5 ,d6 ,d7 ,d8 ,d9 ,i2,i3,
     .                  i4,i5,i6,d10,d11,d12,d13,d14,d15,d16)
      IMPLICIT none
      INTEGER          i1 ,i2 ,i3, i4,i5,i6
      DOUBLE PRECISION d1 ,d2 ,d3 ,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,
     .                 d14,d15,d16
      WRITE(0,*) 'WARNING: Calling SM1USR routine substitute'
      RETURN
      END
c
c
c
      REAL*8 FUNCTION TABPRC(i1,i2,i3)
      IMPLICIT none
      INTEGER i1,i2,i3
      WRITE(0,*) 'WARNING: Calling TABPRC routine substitute'
      TABPRC = 0.0D0
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
      INTEGER i1,i2
      REAL*4  d1(*),d2(*)
      WRITE(0,*) 'WARNING: Calling KURVEF routine substitute'
      RETURN
      END

      SUBROUTINE GRBLD(r1,r2,i1,i2,r3,r4,r5,r6,i3)
      IMPLICIT none
      INTEGER i1,i2,i3
      REAL    r1,r2,r3,r4,r5,r6
      WRITE(0,*) 'WARNING: Calling GRBLD  routine substitute'
      RETURN
      END


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
      
      IF (i1.EQ.0) count = count + 1

      WRITE(80,'(I8,3F14.6,4I8,1P,1E10.2,0P,F8.4,6I5)') 
     .  count,d1,d2,d3,i1,i2,istra,npanu,e0,weight,nacell,ifpath,iupdte,
     .  ityp,nblock,masurf

c      WRITE(80,'(I8,3F14.6,3I8,1P,1E10.2,0P,F8.4,6I5)') 
c     .  count,d1,d2,d3,i1,i2,npanu,e0,weight,nacell,ifpath,iupdte,ityp,
c     .  nblock,masurf

c      WRITE(80,'(I8,3F14.6,3I8,1P,1E10.2,0P,F8.4,4I5)') 
c     .  count,d1,d2,d3,i1,i2,npanu,e0,weight,nacell,ifpath,iupdte,ityp

c      WRITE( 0,'(I8,3F14.6,3I8,1P,1E10.2,0P,F8.4,4I5)') 
c     .  count,d1,d2,d3,i1,i2,npanu,e0,weight,nacell,ifpath,iupdte,ityp

c      WRITE(0,*) 'WARNING: Calling CHCTRC routine substitute'
      RETURN
      END

      SUBROUTINE GRSTRT(i1,i2)
      IMPLICIT   none
      INTEGER i1,i2
      WRITE(0,*) 'WARNING: Calling GRSTRT routine substitute'
      RETURN
      END

      SUBROUTINE GREND
      STOP 'GREND'
      RETURN
      END

      SUBROUTINE OUTPAT
      WRITE(0,*) 'WARNING: Calling OUTPAT routine substitute'
      RETURN
      END

      SUBROUTINE PLT2D
      WRITE(0,*) 'WARNING: Calling PLT2D  routine substitute'
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

