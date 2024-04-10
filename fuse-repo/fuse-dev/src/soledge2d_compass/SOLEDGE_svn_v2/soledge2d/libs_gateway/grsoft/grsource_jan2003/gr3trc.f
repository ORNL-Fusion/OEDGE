C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
CC@PROCESS OPT(3) NOSDUMP NOGOSTMT

C UPDATE: 22. 3.1991 GROTEN
C
C Dieses Programm muss getrennt uebersetzt und in GRNEU TXTLIB gebracht
C werden, weil es unter Umstaenden von einem Benutzerprogramm ersetzt
C wird.
C
      subroutine gr3trc(u,v,w,x,y,z)
      real u,v,w,x,y,z

      x=u
      y=v
      z=w

      end
