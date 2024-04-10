C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
c     Update: 21. 2.1991 Busch
c     Update: 12.12.1991 J.Heinen
C     GR3 Panel GR-Software stehne fuer UNIX Version nicht zur Verfueg.
C     ersatzweise wird gr3plo aufgerufen
C
      SUBROUTINE GR3PAN(ar)
      integer          asf(13)
      real ar(*)
      data asf/0,0,0,0,0,0,0,0,0,0,0,0,0/
      call gsasf(asf)
      call gr3plo(ar,ierr,'HID')
      call grpan
      end
c
      subroutine gr3fpc(wx, wy, wz, ko)
      real wx, wy, wz
      integer ko
      write(*,*) 'dummy call gr3fpc(',wx,',',wy,',',wz,',',ko,')'
      end
