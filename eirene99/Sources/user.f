*DK USER
C
C   USER SUPPLIED SUBROUTINES
C
C           ************
C           *  ASDEX-U *
C           ************
C
      SUBROUTINE PROUSR (PRO,INDEX,P0,P1,P2,P3,P4,P5,PROVAC,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PRO(*)
      RETURN
      END
C
      SUBROUTINE REFUSR
C
C   USER SUPPLIED REFLECTION MODEL
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ENTRY RF0USR
C
      ENTRY RF1USR (XMW,XCW,XMP,XCP,IGASF,IGAST,ZCOS,ZSIN,EXPI,
     .              RPROB,E0TERM,*,*,*,*)
C
      ENTRY SPTUSR
      RETURN
      END
C
C
C
      SUBROUTINE PLTUSR(PLABLE,J)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL PLABLE
      RETURN
      END
C
C
      SUBROUTINE GEOUSR
C
C   PREPARE DATA FOR LIMITER-SURFACES
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include 'PARMMOD'
      include 'CADGEO'
      include 'CTRCEI'
      include 'CCONA'
      include 'CGEOM'
      include 'CGRID'
      include 'CLGIN'
      include 'CPOLYG'
      include 'CINIT'
      logical first
      integer onetwo(4),limpos(4)
      character*80 geometry_comment
      data first /.true./
      save first,onetwo,limpos,geometry_comment
c slmod begin - f90 - tr - new
c...bug: LGJUM3 appears to be obsolete.
c
C   LGJUM3(J,I)=.TRUE. :
C   ABSCHALTEN DER FLAECHE I, FALLS TEILCHEN IN ZELLE J SITZT
C
c      DO 10 I=1,NLIMI
c        DO 10 J=1,NSURF
c          LGJUM3(J,I)=.TRUE.
c10    CONTINUE
c
      IF (neutopt.EQ.1) RETURN
c
c It is unnecessary to alter geometry data in this
c routine because...
c
c
c slmod end
      if(first) then
        read(50,'(a80)') geometry_comment
        do i=1,4
          read(50,*) onetwo(i),limpos(i)
        enddo
        first=.false.
        write(*,*) 'GEOMETRY FOR'
        write(*,'(a80)') geometry_comment
        write(*,'(''P'',i1,'' FOR SEGMENT '',i3,
     1       '' LINKED TO INNER LEFT TARGET'')') onetwo(1),limpos(1)
        write(*,'(''P'',i1,'' FOR SEGMENT '',i3,
     1       '' LINKED TO OUTER LEFT TARGET'')') onetwo(2),limpos(2)
        write(*,'(''P'',i1,'' FOR SEGMENT '',i3,
     1      '' LINKED TO INNER RIGHT TARGET'')') onetwo(3),limpos(3)
        write(*,'(''P'',i1,'' FOR SEGMENT '',i3,
     1      '' LINKED TO OUTER RIGHT TARGET'')') onetwo(4),limpos(4)
      endif
c sltmp
      WRITE(0,*) 'GEOUSR: TEMP NON-GENERALIZED GEOMETRY SUPPORTED'
      RETURN
      STOP 'GEOUSR: NON-GENERALIZED GEOMETRY NOT SUPPORTED'
C
C  ANFANG: MODIFY GEOMETRY
C
C LEFT LIMITER
      if(onetwo(1).eq.1) then
        p1(1,limpos(1))=xpol(1,1)
        p1(2,limpos(1))=ypol(1,1)
      else
        p2(1,limpos(1))=xpol(1,1)
        p2(2,limpos(1))=ypol(1,1)
      endif
      if(onetwo(2).eq.1) then
        p1(1,limpos(2))=xpol(nr1st,1)
        p1(2,limpos(2))=ypol(nr1st,1)
      else
        p2(1,limpos(2))=xpol(nr1st,1)
        p2(2,limpos(2))=ypol(nr1st,1)
      endif
C RIGHT LIMITER
      if(onetwo(3).eq.1) then
        p1(1,limpos(3))=xpol(1,np2nd)
        p1(2,limpos(3))=ypol(1,np2nd)
      else
        p2(1,limpos(3))=xpol(1,np2nd)
        p2(2,limpos(3))=ypol(1,np2nd)
      endif
      if(onetwo(4).eq.1) then
        p1(1,limpos(4))=xpol(nr1st,np2nd)
        p1(2,limpos(4))=ypol(nr1st,np2nd)
      else
        p2(1,limpos(4))=xpol(nr1st,np2nd)
        p2(2,limpos(4))=ypol(nr1st,np2nd)
      endif
C
C  ENDE: MODIFY GEOMETRY
C
C
C  ABSCHALTEN NICHT ERREICHBARER ODER DOPPELT VORHANDENER FLAECHEN
C
C   HIERHER: LGJUM1, LGJUM2 SETZEN ZUR BESCHLEUNIGUNG (NICHT UNBEDINGT
C   NOETIG)
C
C   LGJUM1(J,I)=.TRUE. :
C   ABSCHALTEN DER FLAECHE I, FALLS TEILCHEN AUF J SITZT
C
C   LGJUM2(J,I)=.TRUE. :
C   ABSCHALTEN DES ERSTEN SCHNITTPUNKTES MIT FLAECHE I, FALLS
C   TEILCHEN AUF J SITZT (FALLS I EINE FLAECHE ZWEITER ORDNUNG IST)
C
C   DEFAULTS: LGJUM1(J,J)=.TRUE. FUER EBENE FLAECHEN,
C             LGJUM2(J,J)=.TRUE. FUER FLAECHEN ZWEITER ORDNUNG
C
C
c slmod begin
c
c This code has been move...
c
cC   LGJUM3(J,I)=.TRUE. :
cC   ABSCHALTEN DER FLAECHE I, FALLS TEILCHEN IN ZELLE J SITZT
cC
c      DO 10, I=1,NLIMI
c        DO 10,J=1,NSURF
c          LGJUM3(J,I)=.TRUE.
c10    CONTINUE
c slmod end
c
C
C
C
C  SET SOME VOLUMES EXPLIZIT
C
C
C  MODIFY REFLECTION MODEL AT TARGET PLATES
C
C
      RETURN
      END
C
      SUBROUTINE SAMUSR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'CADGEO'


      INTEGER count,numtry
      LOGICAL found,CHKPNT

      ENTRY SM1USR(NLSF,X0,Y0,Z0,INSOR,
     .             SORAD1,SORAD2,
     .             SORAD3,SORAD4,
     .             SORAD5,SORAD6,
     .             IRUSR,IPUSR,ITUSR,IAUSR,IBUSR,
     .             TIWL,TEWL,DIWL,VXWL,VYWL,VZWL,WEISPZ)


c...  To qualify for this lauch routine, the surface must be
c     perpendicular to the z-axis:
      IF (p1(3,insor).NE.p2(3,insor).OR.p2(3,insor).NE.p3(3,insor).OR.
     .    p3(3,insor).NE.p4(3,insor)) THEN
        WRITE(0,*) 'SAMUSR: SURFACE IS INVALID'
        WRITE(6,*) 'SAMUSR: SURFACE IS INVALID'
        CALL EXIT
      ENDIF

      found = .FALSE.

      numtry = 0

      DO WHILE (.NOT.found)

        numtry = numtry + 1

c...    Sample SORADx to get lanuch point (X0,Y0,Z0):
        x0 = sorad1 + RANF() * (sorad2 - sorad1)
        y0 = sorad3 + RANF() * (sorad4 - sorad3)
        z0 = sorad5

c...    Decide if the launch point is inside the polygon or not:
        a1 = x0
        a2 = y0
        b1 = x0 
        b2 = y0 + 100.0D0

        count = 0

        DO i1 = 1, 4
          IF     (i1.EQ.1) THEN
            c1 = p1(1,insor)
            c2 = p1(2,insor)
            d1 = p2(1,insor)
            d2 = p2(2,insor)
          ELSEIF (i1.EQ.2) THEN
            c1 = p2(1,insor)
            c2 = p2(2,insor)
            d1 = p4(1,insor)
            d2 = p4(2,insor)
          ELSEIF (i1.EQ.3) THEN
            c1 = p4(1,insor)
            c2 = p4(2,insor)
            d1 = p3(1,insor)
            d2 = p3(2,insor)
          ELSEIF (i1.EQ.4) THEN
            c1 = p3(1,insor)
            c2 = p3(2,insor)
            d1 = p1(1,insor)
            d2 = p1(2,insor)
          ENDIF

c          WRITE(0,'(A,8F12.4)') 'ABCD:',a1,a2,b1,b2,c1,c2,d1,d2

          IF (CHKPNT(a1,a2,b1,b2,c1,c2,d1,d2)) count = count + 1

c          WRITE(0,*) 'COUNT:',count
        ENDDO      

c        WRITE(0,'(A,3F12.4,2I6)') '-->',x0,y0,z0,count,numtry

        IF     (count.EQ.1) THEN
c...      Launch point is inside polygon:
          found = .TRUE.
        ELSEIF (count.NE.2.AND.count.NE.0) THEN
c...      Problem:
          WRITE(0,*) 'SAMUSR: INVALID INTERECTION COUNT'
          WRITE(0,*) 'SAMUSR: INVALID INTERECTION COUNT'
          CALL EXIT
        ELSEIF (numtry.GE.1000) THEN
c...      Problem:
          WRITE(0,*) 'SAMUSR: CANNOT FIND VALID LAUNCH POINT'
          WRITE(0,*) 'SAMUSR: CANNOT FIND VALID LAUNCH POINT'
          CALL EXIT
        ENDIF        

      ENDDO

      RETURN

c slmod end
      WRITE(0,*) 'ERROR: Calling SAMUSR'
      WRITE(6,*) 'ERROR: Calling SAMUSR'
      STOP 'SAMUSR'
      RETURN
      END
C
      SUBROUTINE MODUSR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
c slmod begin
      INCLUDE 'PARMMOD'
      INCLUDE 'COMUSR'
      INCLUDE 'CPOLYG'
      INCLUDE 'CCOUPL'
      INCLUDE 'CGEOM'
      INCLUDE 'CGRID'

      DIMENSION DUMMY(0:NDXP,0:NDYP),SNID(0:NDXP,0:NDYP,NFL)

      IF (NITER.GT.0) CALL MODBGK

c      CALL ASDUSR

c      WRITE(0,*) 'NOTHING IN MODUSR'

c       ===============================================================
c
c       Output data for BGK iterations.  Some assumptions about the
c       number of ion species is being made, namely that there are 6:
c       D+, C+, D, D2, DD2, D2D.  Only the last 4 are to be saved here
c       and subsequently read by DIVIMP (they will not be processed to
c       the degree that the other PIN quantities are, in order to keep
c       things simple).
c
c
c...use this as the flag to identify when BGK is required -- good?
        IF (NITER.GT.0.AND.IITER.EQ.NITER.AND.ITIMV.GE.NTIME) THEN

          IF (output) WRITE(0,*) 'WRITING BGK BACKGROUND PLASMA DATA'

          DO IT=1,NBMLT

c
c           ION DENSITY
c
            IF (output) WRITE(0 ,*) '[BGK: ION DENSITY]'

            CALL DSET(SNID(0,0,1),(NDXP+1)*(NDYP+1)*NFL,0.0D0)
            DO IFL=1,NFLA
              DO IPLS=1,NPLSI
                IF (IFLB(IPLS).NE.IFL.AND.IFLB(IPLS)+30.NE.IFL) CYCLE
                IF (output) WRITE(0,*) 'IFL,IPLS= ',IFL,IPLS
                DO IX=1,NDXA
                  DO IY=1,NDYA
                    IN=IY+(IX-1)*NR1ST+(IT-1)*NSTRD
                    SNID(IX,IY,IFL)=DIIN(IPLS,IN)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
            IF (NCUTL.NE.NCUTB)
     .        CALL INDMPI (SNID,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .                     NCUTB,NCUTL,NPOINT,NPPLG,1,1)

            NRED=(NPPLG-1)*(NCUTL-NCUTB)
            NDX2=NDXA-NRED
            WRITE(32,'(A,I3,A1)') '[BGK: ION DENSITY, REGION ',IT,']'
            CALL NEUTR(32,NDX2,NDYA,NFLA,SNID,NDX,NDY,NFL,1,1)
c
c           ION TEMPERATURE
c
            IF (output) WRITE(0 ,*) '[BGK: ION TEMPERATURE]'

            CALL DSET(SNID(0,0,1),(NDXP+1)*(NDYP+1)*NFL,0.0D0)
            DO IFL=1,NFLA
              DO IPLS=1,NPLSI
                IF (IFLB(IPLS).NE.IFL.AND.IFLB(IPLS)+30.NE.IFL) CYCLE
                DO IX=1,NDXA
                  DO IY=1,NDYA
                    IN=IY+(IX-1)*NR1ST+(IT-1)*NSTRD
                    SNID(IX,IY,IFL)=TIIN(IPLS,IN)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
            IF (NCUTL.NE.NCUTB)
     .        CALL INDMPI (SNID,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .                     NCUTB,NCUTL,NPOINT,NPPLG,1,1)
            NRED=(NPPLG-1)*(NCUTL-NCUTB)
            NDX2=NDXA-NRED
            WRITE(32,'(A,I3,A)') '[BGK: ION TEMPERATURE, REGION ',IT,']'
            CALL NEUTR(32,NDX2,NDYA,NFLA,SNID,NDX,NDY,NFL,1,1)
c
c           X-VELOCITY
c
            IF (output) WRITE(0 ,*) '[BGK: X-VELOCITY]'

            CALL DSET(SNID(0,0,1),(NDXP+1)*(NDYP+1)*NFL,0.0D0)
            DO IFL=1,NFLA
              DO IPLS=1,NPLSI
                IF (IFLB(IPLS).NE.IFL.AND.IFLB(IPLS)+30.NE.IFL) CYCLE
                DO IX=1,NDXA
                  DO IY=1,NDYA
                    IN=IY+(IX-1)*NR1ST+(IT-1)*NSTRD
                    SNID(IX,IY,IFL)=VXIN(IPLS,IN)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
            IF (NCUTL.NE.NCUTB)
     .        CALL INDMPI (SNID,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .                     NCUTB,NCUTL,NPOINT,NPPLG,1,1)
            NRED=(NPPLG-1)*(NCUTL-NCUTB)
            NDX2=NDXA-NRED
            WRITE(32,'(A,I3,A)') '[BGK: X-VELOCITY, REGION ',IT,']'
            CALL NEUTR(32,NDX2,NDYA,NFLA,SNID,NDX,NDY,NFL,1,1)
c
c           Y-VELOCITY
c
            IF (output) WRITE(0 ,*) '[BGK: Y-VELOCITY]'

            CALL DSET(SNID(0,0,1),(NDXP+1)*(NDYP+1)*NFL,0.0D0)
            DO IFL=1,NFLA
              DO IPLS=1,NPLSI
                IF (IFLB(IPLS).NE.IFL.AND.IFLB(IPLS)+30.NE.IFL) CYCLE
                DO IX=1,NDXA
                  DO IY=1,NDYA
                    IN=IY+(IX-1)*NR1ST+(IT-1)*NSTRD
                    SNID(IX,IY,IFL)=VYIN(IPLS,IN)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
           IF (NCUTL.NE.NCUTB)
     .        CALL INDMPI (SNID,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .                     NCUTB,NCUTL,NPOINT,NPPLG,1,1)
            NRED=(NPPLG-1)*(NCUTL-NCUTB)
            NDX2=NDXA-NRED
            WRITE(32,'(A,I3,A)') '[BGK: Y-VELOCITY, REGION ',IT,']'
            CALL NEUTR(32,NDX2,NDYA,NFLA,SNID,NDX,NDY,NFL,1,1)
c
c           Z-VELOCITY
c
            IF (output) WRITE(0 ,*) '[BGK: Z-VELOCITY]'

            CALL DSET(SNID(0,0,1),(NDXP+1)*(NDYP+1)*NFL,0.0D0)
            DO IFL=1,NFLA
              DO IPLS=1,NPLSI
                IF (IFLB(IPLS).NE.IFL.AND.IFLB(IPLS)+30.NE.IFL) CYCLE
                DO IX=1,NDXA
                  DO IY=1,NDYA
                    IN=IY+(IX-1)*NR1ST+(IT-1)*NSTRD
                    SNID(IX,IY,IFL)=VZIN(IPLS,IN)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
            IF (NCUTL.NE.NCUTB)
     .        CALL INDMPI (SNID,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .                     NCUTB,NCUTL,NPOINT,NPPLG,1,1)
            NRED=(NPPLG-1)*(NCUTL-NCUTB)
            NDX2=NDXA-NRED
            WRITE(32,'(A,I3,A)') '[BGK: Z-VELOCITY, REGION ',IT,']'
            CALL NEUTR(32,NDX2,NDYA,NFLA,SNID,NDX,NDY,NFL,1,1)

          ENDDO

          IF (output) WRITE(0,*) 'DONE'

        ELSE

          IF (output) WRITE(0,*) 'NOT WRITING BGK DATA TO DIVIMP'

        ENDIF
c
c      IF (NFILEL.EQ.1) THEN
c        WRITE(0,*) 'MODIFYING NFILEL'
c        WRITE(6,*) 'MODIFYING NFILEL'
c
c        NFILEL = 2
c      ENDIF
c slmod end
      RETURN
      END
C
C
      SUBROUTINE UPTUSR(XSTOR2,WV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  USER SUPPLIED TRACKLENGTH ESTIMATOR, VOLUME AVERAGED
C
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'CUPD'
      INCLUDE 'COMXS'
      INCLUDE 'CSPEZ'
      INCLUDE 'COMUSR'
      INCLUDE 'CESTIM'
      INCLUDE 'CGRID'
      INCLUDE 'CLOGAU'
      INCLUDE 'CCONA'
      INCLUDE 'CPOLYG'
      INCLUDE 'CZT1'
      DIMENSION CFLAG(6,3)
      DIMENSION XSTOR(NSTOR),XSTOR2(NSTOR,N2ND+N3RD)
      EQUIVALENCE (XSTOR(1),SIGVCX(1))
      DIMENSION CNDYNA(NATM),CNDYNP(NPLS)
C
      RETURN
      END
C
      SUBROUTINE UPSUSR(WT,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'COMUSR'
      INCLUDE 'CLGIN'
      INCLUDE 'COMPRT'
      INCLUDE 'CESTIM'
      INCLUDE 'CGRID'
      INCLUDE 'CADGEO'


      RETURN


c      WRITE(0,*) 'MSURF',msurf

c      IF (y0.GT.-7.0D0) THEN
c        WRITE(0,*) 'MRSURF,NCELL=',msurf,ncell
c      ENDIF


c...  Record time to absorption for all absorbing 
c     additional surfaces:
      IF (MSURF.LE.NLIMI.AND.ILIIN(MSURF).EQ.2) THEN
        slADDS(0,MSURF)=slADDS(0,MSURF)+WT
        slADDS(1,MSURF)=slADDS(1,MSURF)+WT*TIME

c        WRITE(0,*) 'DATA:',MSURF,slADDS(0,MSURF),slADDS(1,MSURF)

      ENDIF





      RETURN

C  SPATIALLY RESOLVED IATM=1 FLUXES AT TARGET SURFACE NO. 1
      IF (IND.NE.1) RETURN
c      WRITE(0,*) 'NRCELL=',nrcell,wt,e0,ityp,iatm,mpsurf,np2nd
      MS=MSURF-NLIM
      IF (MS.EQ.5) THEN
        IF (ITYP.EQ.1.AND.IATM.EQ.1) THEN
          ADDV(1,NRCELL)=ADDV(1,NRCELL)+WT
          ADDV(3,NRCELL)=ADDV(3,NRCELL)+WT*E0
        ENDIF
      ENDIF
C  SPATIALLY RESOLVED IATM=1 FLUXES AT TARGET SURFACE NO. 2
      MS=MSURF-NLIM
      IF (MS.EQ.6.OR.MS.EQ.7.OR.MS.EQ.8.OR.MS.EQ.9) THEN
        IF (ITYP.EQ.1.AND.IATM.EQ.1) THEN
          ADDV(2,NRCELL)=ADDV(2,NRCELL)+WT
          ADDV(4,NRCELL)=ADDV(4,NRCELL)+WT*E0
        ENDIF
      ENDIF


      RETURN
      END
C
      SUBROUTINE UPCUSR(WS,IND)
C
C  USER SUPPLIED COLLISION ESTIMATOR, VOLUME AVERAGED
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      STOP 'UPCUSR called in USER.F'

      RETURN
      END
C
C
c
      subroutine sigusr(i1,i2,d1,d2,d3,d4,d5)
      real*8 d1,d2,d3(*),d4,d5(*)
      return
      end
c
c
      subroutine w7xint
      return
      end

c
c
      subroutine upnusr
      return
      end
c
      SUBROUTINE RETUSR(SIG)
      IMPLICIT REAL*8 (A-H,O-Z)
      RETURN
      END
*//PLASM//
C=======================================================================
C          S U B R O U T I N E   P L A S M
C=======================================================================
      SUBROUTINE PLASM(KARD,NDIMX,NDIMY,NDIMF,N,M,NF,DUMMY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION DUMMY(0:N+1,0:M+1,NF)
c slmod begin
      INCLUDE 'PARMMOD'

      CHARACTER COMMENT2*80,BUFFER2*256

      READ (31,'(A)') COMMENT2
      WRITE( 6,'(A)') ' PLASM: '//comment2
      IF (output)
     .WRITE( 0,'(A)') ' PLASM: '//comment2
c      WRITE( 0,'(A)') ' PLASM: '//comment2
c slmod end
      ND1 = NDIMX + 2
      LIM = (ND1/5)*5 - 4
      DO    110  IF = 1,NDIMF
      DO    110  IY = 0,NDIMY+1
      DO    100  IX = 1,LIM,5
100     READ(KARD,910) (DUMMY(-1+IX-1+III,IY,IF),III = 1,5)
        IF( (LIM+4).EQ.ND1 )     GOTO 110
        READ(KARD,910) (DUMMY(-1+IX,IY,IF),IX = LIM+5,ND1)
110   CONTINUE
      RETURN
910   FORMAT(5(E16.8))
*//END PLASM//
      END
C
C
*//NEUTR//
C=======================================================================
C          S U B R O U T I N E   N E U T R
C=======================================================================
      SUBROUTINE NEUTR(KARD,NDIMX,NDIMY,NDIMF,DUMMY,LDMX,LDMY,LDMF,
     .                 LDNS,IS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER KARD,NDIMX,NDIMY,NDIMF,LDMX,LDMY,LDMF,LDNS,
     .        ND1,LIM,IX,IY,III,IF,IS
      REAL*8 DUMMY(0:LDMX+1,0:LDMY+1,LDMF,LDNS)
C
      ND1 = NDIMX
      LIM = (ND1/5)*5 - 4
      DO  500  IF = 1,NDIMF
        DO  110  IY = 1,NDIMY
          DO  100  IX = 1,LIM,5
  100     WRITE(KARD,910) (DUMMY(IX-1+III,IY,IF,IS),III = 1,5)
          IF( (LIM+4).EQ.ND1 )   GOTO 110
          WRITE(KARD,910) (DUMMY(IX,IY,IF,IS),IX = LIM+5,ND1)
  110   CONTINUE
  500 CONTINUE
      RETURN
  910 FORMAT(5(E16.8))
*//END NEUTR//
      END
C
      SUBROUTINE MSHPROJ(X1,Y1,X2,Y2,X3,Y3,X4,Y4,PUX,PUY,PVX,PVY,
     .                   NDXA,NR1ST,IY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X1(*),Y1(*),X2(*),Y2(*),X3(*),Y3(*),X4(*),Y4(*),PUX(*),
     .          PUY(*),PVX(*),PVY(*)
      EPS60 = 1.D-60
C
C
      DO 1 IX=1,NDXA
C
C  CALCULATE THE NORM OF THE VECTORS (POINT2-POINT1),....
C
        D12 = SQRT((X2(IX)-X1(IX))*(X2(IX)-X1(IX))+(Y2(IX)-Y1(IX))*
     .        (Y2(IX)-Y1(IX)))+EPS60
        D34 = SQRT((X4(IX)-X3(IX))*(X4(IX)-X3(IX))+(Y4(IX)-Y3(IX))*
     .        (Y4(IX)-Y3(IX)))+EPS60
        D13 = SQRT((X3(IX)-X1(IX))*(X3(IX)-X1(IX))+(Y3(IX)-Y1(IX))*
     .        (Y3(IX)-Y1(IX)))+EPS60
        D24 = SQRT((X4(IX)-X2(IX))*(X4(IX)-X2(IX))+(Y4(IX)-Y2(IX))*
     .        (Y4(IX)-Y2(IX)))+EPS60
C
C  CALCULATE THE BISSECTING VECTORS, BUT NOT NORMALISED YET
C
        DUX = (X2(IX)-X1(IX))/D12 + (X4(IX)-X3(IX))/D34
        DUY = (Y2(IX)-Y1(IX))/D12 + (Y4(IX)-Y3(IX))/D34
        DVX = (X3(IX)-X1(IX))/D13 + (X4(IX)-X2(IX))/D24
        DVY = (Y3(IX)-Y1(IX))/D13 + (Y4(IX)-Y2(IX))/D24
C
C  CALCULATE THE COMPONENTS OF THE TWO UNIT VECTOR (= PROJECTION RATE)
C
        IN=IY+(IX-1)*NR1ST
        PUX(IN) = DUX/(SQRT(DUX*DUX+DUY*DUY)+EPS60)
        PUY(IN) = DUY/(SQRT(DUX*DUX+DUY*DUY)+EPS60)
        PVX(IN) = DVX/(SQRT(DVX*DVX+DVY*DVY)+EPS60)
        PVY(IN) = DVY/(SQRT(DVX*DVX+DVY*DVY)+EPS60)
C
C  ORTHOGONORMALIZE, CONSERVE ORIENTATION (E.SCHMIDT)
C
        PUPV=PUX(IN)*PVX(IN)+PUY(IN)*PVY(IN)
        PVX(IN)=PVX(IN)-PUPV*PUX(IN)
        PVY(IN)=PVY(IN)-PUPV*PUY(IN)
        PVPV=SQRT(PVX(IN)*PVX(IN)+PVY(IN)*PVY(IN))+EPS60
        PVX(IN)=PVX(IN)/PVPV
        PVY(IN)=PVY(IN)/PVPV
C
1     CONTINUE
      RETURN
      END
C
C
      function leausr(a,b,c)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      leausr=1
      RETURN
      END
      subroutine timusr(n,x,y,z,vx,vy,vz,n1,n2,t,ic,ie)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      entry norusr(i1,d1,d2,d3,d4,d5,d6,d7,i2)
      RETURN
      END
      subroutine volusr(n,a)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension a(*)
      RETURN
      END
c slmod begin
      SUBROUTINE PLAUSR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      RETURN
      END
c slmod end
      subroutine TMSUSR(TIME0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      RETURN
      END
c slmod begin
c
*//GEOMD//
C=======================================================================
C          S U B R O U T I N E   G E O DIV
C=======================================================================
c
c I am going to try and 'pad' or 'stuff' the grid with extra
c rings to store the geometry data needed for non structured grids.
c There is some suggestion that this is already possible with
c EIRENE with the reference below to NYCUT, but I have not seen
c any implimentation of it in the rest of the code so far.
c
c The trick is to make these padding cells invisible to EIRENE.  If
c there was a flag that I could set to make them invalid, or to make
c EIRENE think that they are in the wall or something.
c
c I will modify the cell finding routine so that it does not look at
c those cells, and set their volume to zero so that they do not
c contribute to the statistics in any way.  I will also have to make
c sure that no particles are launched from target segments belonging
c to those cells.
c
c I doubt that I will be able to make them completely invisible without
c a great deal of work but it is worth a shot.
c
c NDXA and NDYA are read in from the stream IUNIN=50, which is the
c .dat EIRENE input file.
c
c Things to do:
c
c 1) Make sure that this routine can read in the new geometry file
c    generated by DIVIMP
c
      SUBROUTINE GEODIV(NDXA,NDYA,XPLG,YPLG,NPLP,NPNT,NR1ST,
     .                   PUX,PUY,PVX,PVY,NBMLT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include 'PARMMOD'
      include 'CGEOM'
      INCLUDE 'CTRCEI'
      INCLUDE 'CLOGAU'
      DIMENSION XPLG(N1ST,*),YPLG(N1ST,*),NPNT(2,*),
     .          PUX(*),PUY(*),PVX(*),PVY(*)
C
C  GEOMETRY DATA: CELL VERTICES (LINDA ---> EIRENE)
      COMMON /LINEIR/
     R  X1(NDX),Y1(NDX),X2(NDX),Y2(NDX),X3(NDX),Y3(NDX),
     R  X4(NDX),Y4(NDX)
C
      CHARACTER*80 LINE,LINE2
C   DIMENSIONIERUNG FUER GITTER
      INTEGER DIMXH,DIMYH,NNCUT,NXCUT1,NXCUT2,d1,d2,d3,dummy
      CHARACTER c1*12
      DIMENSION DUMMI(3)
      REAL*8 MERK(NDY)
      INTEGER k


      DOUBLE PRECISION xplg2(N1ST,N2ND),yplg2(N1ST,N2ND)


      INTEGER fp,i1,i2,idum1

C  ACTUAL MESH USED IN THIS RUN
C
       REWIND 30

3366  FORMAT(/)

c Temporary...
c       TRCHST = .FALSE.


c
c Read file header:
      IF (debugopt.NE.0)
     .  WRITE(0,*) 'GEODIV: Reading geometry file'

      READ(30,'(A60)'   ) comment
      READ(30,*)
      READ(30,'(A30,I6)') comment,dimxh
      READ(30,'(A30,I6)') comment,dimyh
      READ(30,'(A30,I6)') comment,nncut
      READ(30,'(A30,I6)') comment,nxcut1
      READ(30,'(A30,I6)') comment,nxcut2
      READ(30,'(A30,I6)') comment,dummy
      READ(30,'(A30,I6)') comment,dummy
      READ(30,'(A30,I6)') comment,dummy
      READ(30,*)


      NPLP = 3


      NPNT(1,1) = 1
      NPNT(2,1) = nxcut1
      NPNT(1,2) = nxcut1 + 1
      NPNT(2,2) = nxcut2 + 1
      NPNT(1,3) = nxcut2 + 2
      NPNT(2,3) = dimxh + 3
c
c     Read in geometry data:

      DO ix= 1,dimxh+3
        DO iy = 1, dimyh+1

          READ(30,*,END=25,ERR=20)
     +      xvert(iy,ix,1),yvert(iy,ix,1),
     +      xvert(iy,ix,2),yvert(iy,ix,2),
     +      d1,d2

          xplg(iy,ix) = xvert(iy,ix,1)
          yplg(iy,ix) = yvert(iy,ix,1)
        ENDDO
      ENDDO

c...  Load connection map data - new:
      WRITE(6,*) 'CONNECTION MAP:'
      READ(30,*)
      READ(30,*)
      READ(30,*)
      READ(30,*)
      READ(30,*)
      DO ix= 1,dimxh+3
        READ(30,*,END=25,ERR=20) idum1,(nghpol2(1,iy,ix),iy=1,dimyh)
        WRITE(6,'(100(I6:))') idum1,(nghpol2(1,iy,ix),iy=1,dimyh)
      ENDDO
      READ(30,*)
      READ(30,*)
      READ(30,*)
      DO ix= 1,dimxh+3
        READ(30,*,END=25,ERR=20) idum1,(nghpol2(2,iy,ix),iy=1,dimyh)
        WRITE(6,'(100(I6:))') idum1,(nghpol2(2,iy,ix),iy=1,dimyh)
      ENDDO

3334  FORMAT(4D15.7,5X,2I6,A2)
      goto 30

20    CONTINUE
      write(0,*) 'GEODIV: Error reading geometry file',
     .  ' (IX,IY) ',ix,iy
      GOTO 30

25    CONTINUE
      write(0,*) 'GEODIV: Error EOF while reading geometry file',
     .  ' (IX,IY) ',ix,iy
30    CONTINUE


c...  Also load up the non-default radial surface map which tells
c     EIRENE which side of ring to use when checking for radial
c     surface collisions:
      READ(30,*)
      READ(30,*)
      READ(30,*)
      DO iy= 1,dimyh+1
        READ(30,*,END=35,ERR=40) idum1,radmap(iy)
c        WRITE(0,*) 'RADMAP:',iy,radmap(iy)
      ENDDO
      GOTO 40
35    CONTINUE
      WRITE(0,*) 'GEODIV: Error EOF while radial map'
40    CONTINUE


      IF (GRIDOPT.EQ.1) THEN
c
c Calculate PVRTAG:
c
        pvrtag = 0
        rvrtag = 0


        DO iy = 1, dimyh+1
          DO ii = 1, nplp
            DO ix = npnt(1,ii),npnt(2,ii)-1
c...          Skip the first poloidal surface on each ring:
              IF (ix.EQ.1) CYCLE

              ix0 = ix - 1   
              IF (ix0.LT.npnt(1,ii).AND.ii.GT.1) ix0 = npnt(2,ii-1)-1
              ix1 = ix + 1   

              IF ((xvert(iy,ix,1).EQ.xvert(iy,ix0,1).AND.
     .             yvert(iy,ix,1).EQ.yvert(iy,ix0,1).AND.
     .             xvert(iy,ix,2).EQ.xvert(iy,ix0,2).AND.
     .             yvert(iy,ix,2).EQ.yvert(iy,ix0,2)).OR.
     .            (xvert(iy,ix,1).EQ.xvert(iy,ix1,1).AND.
     .             yvert(iy,ix,1).EQ.yvert(iy,ix1,1).AND.
     .             xvert(iy,ix,2).EQ.xvert(iy,ix1,2).AND.
     .             yvert(iy,ix,2).EQ.yvert(iy,ix1,2))) 
     .        pvrtag(iy,ix) = 1

              WRITE(6,*) 'PVRTAG:',iy,ix,pvrtag(iy,ix)

            ENDDO
          ENDDO
        ENDDO
 
              
c         DO ix = npnt(1,1),npnt(2,1)-1
c           IF (xvert(iy,ix,1).EQ.xvert(iy,ix+1,1).AND.
c    .          yvert(iy,ix,1).EQ.yvert(iy,ix+1,1).AND.
c    .          xvert(iy,ix,2).EQ.xvert(iy,ix+1,2).AND.
c    .          yvert(iy,ix,2).EQ.yvert(iy,ix+1,2)) THEN
c             pvrtag(iy,ix+1) = 1
c           ENDIF
c         ENDDO
c         DO ix = npnt(1,2),npnt(2,2)-1
c           IF (xvert(iy,ix,1).EQ.xvert(iy,ix+1,1).AND.
c    .          yvert(iy,ix,1).EQ.yvert(iy,ix+1,1).AND.
c    .          xvert(iy,ix,2).EQ.xvert(iy,ix+1,2).AND.
c    .          yvert(iy,ix,2).EQ.yvert(iy,ix+1,2)) THEN
c             pvrtag(iy,ix) = 1
c           ENDIF
c         ENDDO
c         DO ix = npnt(1,3),npnt(2,3)-1
c           IF (xvert(iy,ix,1).EQ.xvert(iy,ix+1,1).AND.
c    .          yvert(iy,ix,1).EQ.yvert(iy,ix+1,1).AND.
c    .          xvert(iy,ix,2).EQ.xvert(iy,ix+1,2).AND.
c    .          yvert(iy,ix,2).EQ.yvert(iy,ix+1,2)) THEN
c             pvrtag(iy,ix) = 1
c           ENDIF
c         ENDDO
c       ENDDO
c
c Calcule RVRTAG:
c
        DO iy = 1, dimyh+1
          DO ix = npnt(1,1),npnt(2,1)-1
            IF (xvert(iy,ix,1).EQ.xvert(iy,ix+1,1).AND.
     .          yvert(iy,ix,1).EQ.yvert(iy,ix+1,1).AND.
     .          xvert(iy,ix,2).EQ.xvert(iy,ix+1,2).AND.
     .          yvert(iy,ix,2).EQ.yvert(iy,ix+1,2)) THEN
              rvrtag(iy,ix) = 1
            ENDIF
          ENDDO

c - new
          DO ix = npnt(1,2),npnt(2,2)-1
            IF (xvert(iy,ix,1).EQ.xvert(iy,ix+1,1).AND.
     .          yvert(iy,ix,1).EQ.yvert(iy,ix+1,1).AND.
     .          xvert(iy,ix,2).EQ.xvert(iy,ix+1,2).AND.
     .          yvert(iy,ix,2).EQ.yvert(iy,ix+1,2)) THEN
              rvrtag(iy,ix) = 1
            ENDIF
          ENDDO

          DO ix = npnt(1,3),npnt(2,3)-1
            IF (xvert(iy,ix,1).EQ.xvert(iy,ix+1,1).AND.
     .          yvert(iy,ix,1).EQ.yvert(iy,ix+1,1).AND.
     .          xvert(iy,ix,2).EQ.xvert(iy,ix+1,2).AND.
     .          yvert(iy,ix,2).EQ.yvert(iy,ix+1,2)) THEN
              rvrtag(iy,ix) = 1
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      WRITE(60,*)
      WRITE(60,*) 'Mesh projection:'
      WRITE(60,*)

      DO 1015 IY=1,NDYA
        WRITE(60,*)
        DO 1014 IX=1,NDXA
          IF (gridopt.EQ.1) THEN
            X1(IX)=XVERT(IY,IX  ,1)
            Y1(IX)=YVERT(IY,IX  ,1)
            X2(IX)=XVERT(IY,IX+1,1)
            Y2(IX)=YVERT(IY,IX+1,1)
            X3(IX)=XVERT(IY,IX  ,2)
            Y3(IX)=YVERT(IY,IX  ,2)
            X4(IX)=XVERT(IY,IX+1,2)
            Y4(IX)=YVERT(IY,IX+1,2)
          ELSE

          ENDIF

          WRITE(60,'(2I4,8E11.4)')
     +     iy,ix,
     +     x1(ix),y1(ix),x2(ix),y2(ix),x3(ix),y3(ix),x4(ix),y4(ix)

1014    CONTINUE
        CALL MSHPROJ (X1,Y1,X2,Y2,X3,Y3,X4,Y4,PUX,PUY,PVX,PVY,NDXA,
     .                NR1ST,IY)
1015  CONTINUE
C
c
c       Output geometry data for debugging:
c
        WRITE(60,*)
        WRITE(60,*) 'Geometry data:'

        DO ix = 1, dimyh+1
          WRITE(60,*)
          DO iy = 1, ndya+1
            WRITE(60,'(4F13.8,1X,2I3,5X,2I4,2X,2I4)')
     +      xvert(iy,ix,1),yvert(iy,ix,1),
     +      xvert(iy,ix,2),yvert(iy,ix,2),
     +      pvrtag(iy,ix),rvrtag(iy,ix),iy,ix,iy+1,ix
          ENDDO
        ENDDO


        WRITE(60,*)
        WRITE(60,*) 'Geometry data (DIVIMP style):'

        DO iy = 1, ndya+1
         WRITE(60,*)
         DO ix = 1, dimxh+3
            WRITE(60,'(4F13.8,1X,2I3,5X,2I4,2X,2I4)')
     +      xvert(iy,ix,1),yvert(iy,ix,1),
     +      xvert(iy,ix,2),yvert(iy,ix,2),
     +      pvrtag(iy,ix),rvrtag(iy,ix),ix,iy,ix,iy+1
          ENDDO
        ENDDO




      NP=NPNT(2,NPLP)
      DO 1020 J=1,NDYA+1
        DO 1020 I=1,NP
          XPLG(J,I)=XPLG(J,I)*100.
          YPLG(J,I)=YPLG(J,I)*100.
c          IF (ABS(XPLG(J,I)).LT.5.D-5) XPLG(J,I)=0.
c          IF (ABS(YPLG(J,I)).LT.5.D-5) YPLG(J,I)=0.

          IF (gridopt.EQ.1) THEN
            DO k = 1, 2
              xvert(j,i,k) = xvert(j,i,k) * 100.0
              yvert(j,i,k) = yvert(j,i,k) * 100.0

c              IF (ABS(xvert(j,i,k)).LT.5.D-5) xvert(j,i,k) = 0.0
c              IF (ABS(yvert(j,i,k)).LT.5.D-5) yvert(j,i,k) = 0.0
            ENDDO
          ENDIF

1020  CONTINUE


      close(60)

      WRITE(61,*) 'GEOMETRY DATA - GEODIV'
      WRITE(61,*)
      WRITE(61,*) 'ndxa   ',ndxa
      WRITE(61,*) 'ndya   ',ndya
      WRITE(61,*) 'nplp   ',nplp
      WRITE(61,*)
      WRITE(61,*) '1: ',npnt(1,1),npnt(2,1)
      WRITE(61,*) '2: ',npnt(1,2),npnt(2,2)
      WRITE(61,*) '3: ',npnt(1,3),npnt(2,3)

      DO ix = 1, npnt(2,nplp)
         WRITE(61,*)
         DO iy = 1, ndya+1
          WRITE(61,'(2F10.4,2I6)') xplg(iy,ix),yplg(iy,ix),iy,ix
        ENDDO
      ENDDO
      CLOSE(61)

c...  READ VACUUM GRID POLYGONS:

      optvac = 0

      fp = 98
      OPEN(UNIT=fp,FILE='vacuum_grid.dat',ACCESS='SEQUENTIAL',
     .     STATUS='OLD',ERR=99)

      optvac = 1

      READ(fp,*) 
      READ(fp,*) ascncell,idum1,asc3dmode,asccode,ascncut
c      WRITE(0,*) 'MARK: ASC3DMODE=',asc3dmode,ascncell,idum1
      DO i1 = 1, ascncell       
        READ(fp,*) idum1,asccell(i1),ascnvertex(i1),ascregion(i1)
        READ(fp,*) asczmin3D(i1),asczmax3D(i1)
        DO i2 = 1, ascnvertex(i1)
          READ(fp,*) ascvertex(2*i2-1,i1),ascvertex(2*i2,i1)
        ENDDO
c        IF (i1.EQ.ascncell) THEN
c          DO i2 = 1, ascnvertex(i1)
c            WRITE(6,*) '-->',ascvertex(2*i2-1,i1),ascvertex(2*i2,i1)
c          ENDDO
c        ENDIF
      ENDDO

c...  Attempt to read in the location of the toroidal surfaces
c     that do NBLOCK switching:

      IF (NLMLT) THEN
        READ(fp,*) 
        READ(fp,*) EIRNSDTOR
        DO I1=1,EIRNSDTOR
          READ(fp,*) EIRSDTOR(I1)
        ENDDO

        WRITE(6,*) 'TOROIDAL NBLOCK SWITHCING SURFACES'
c        WRITE(0,*) 'TOROIDAL NBLOCK SWITHCING SURFACES'
        DO I1=1,EIRNSDTOR
          WRITE(6,*) I1,EIRSDTOR(I1)
c          WRITE(0,*) I1,EIRSDTOR(I1)
        ENDDO

        IF (EIRNSDTOR.NE.NBMLT+1) THEN
          WRITE(0,*) 'WARNING: EIRNSDTOR.NE.NBMLT+1'
          WRITE(6,*) 'WARNING: EIRNSDTOR.NE.NBMLT+1'
        ENDIF
      ENDIF

      CLOSE(FP)

      IF (EIRNSDTOR.GT.0.AND.EIRNTRANS.EQ.0) THEN
        WRITE(0,*) 'ERROR: OLD TRANSPAREND ND-SURFACE SYSTEM'
        WRITE(6,*) 'ERROR: OLD TRANSPAREND ND-SURFACE SYSTEM'
        CALL EXIT
      ENDIF

c      WRITE(0,*) 'MARK: ASC3DMODE=',asc3dmode

c      WRITE(6,*)
c      WRITE(6,*) 'VACUUM GRID POLYGON DATA'
c      WRITE(6,*) '  NUMBER OF CELLS= ',ascncell
c      DO i1 = 1, ascncell       
c        WRITE(6,*) i1,asccell(i1),ascnvertex(i1)
c        DO i2 = 1, ascnvertex(i1)
c          WRITE(6,'(2F14.7)') ascvertex(2*i2-1,i1),ascvertex(2*i2,i1)
c        ENDDO
c      ENDDO



99    CONTINUE

      IF (optvac.GT.0) THEN
        WRITE(6,*) 
        WRITE(6,*) '********************'
        WRITE(6,*) 
        WRITE(6,*) 'OPTVAC= ',optvac
        WRITE(6,*) 
        WRITE(6,*) '********************'
        WRITE(6,*) 
      ENDIF


      IF (debugopt.NE.0)
     .  WRITE(0,*) 'GEODIV: Done'

      RETURN



*//END GEOMD//
      END
c
c
c

      SUBROUTINE FindRadialCell(pt,icase)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'CUPD'
      INCLUDE 'CGEOM'
      INCLUDE 'CPOLYG'
      INCLUDE 'CCONA'

c Replace NX2 declaration with appropriate parameter...
      INTEGER          CheckCell,PolPos

      INTEGER          I1,SI,NX,NX2(100),CHK1
      INTEGER          CalcInter
      DOUBLE PRECISION CP(4)


      STOP 'DEFUNCT'


      IF (gridopt.EQ.0) RETURN
c
c Check that IPOLGN is okay:
c
      IF (GRIDOPT.NE.0.AND.ICASE.EQ.3.AND.MRSURF.NE.0) THEN
        NR1 = NRCELL

c this could still mess up for icase = 1...

        IF (printopt.GE.1.AND.printopt.LE.10)
     .    WRITE(6,'(5X,A,I3,A,2F10.5,A)')
     .    'TIMER: IPOLGN check (NR1 = ',nr1,'  X,Y = ',
     .    X0+PT*VELX,Y0+PT*VELY,')'

        chk1 = CheckCell(NR1,IPOLGN,X0+PT*VELX,Y0+PT*VELY,CP)

        IF (chk1.EQ.1) THEN
          IF (printopt.GE.1.AND.printopt.LE.10)
     +      WRITE(6,'(15X,A)') 'OK'

          IF (cp(2).LT.cp(4)) THEN
            NINCX = +1
          ELSE
            NINCX = -1
          ENDIF

        ELSEIF (chk1.EQ.0) THEN
          chk1 = PolPos(nr1,X0+PT*VELX,Y0+PT*VELY,cp,'TIMER 1')

          IF (chk1.EQ.0) THEN
            IF (printopt.GE.1.AND.printopt.LE.10)
     .        WRITE(6,'(15X,A,I3,A,I4,A,I4,A)')
     .          'Outside grid (NR1 = ',
     .           NR1,'  IPOLGN =',ipolgn,'  STATUS = ',chk1,')'
          ELSEIF (chk1.LT.0) THEN
            IF (printopt.GE.1.AND.printopt.LE.10)
     .        WRITE(6,'(15X,A,I3,A,I4,A,I4,A)')
     .          'ERROR - PolPos cannot find cell (NR1 = ',
     .          NR1,'  IPOLGN =',ipolgn,'  STATUS = ',chk1,')'
            STOP 'Cannot find cell in TIMER'
          ELSE
            IF (printopt.GE.1.AND.printopt.LE.10)
     .        WRITE(6,'(15X,A,I4,A,I4,A)')
     .          'Updating (old IPOLGN = ',
     .          ipolgn,'  new IPOLGN = ',chk1,')'

            IF (cp(2).LT.cp(4)) THEN
              NINCX = +1
            ELSE
              NINCX = -1
            ENDIF

            IPOLGN = CHK1
          ENDIF
        ELSE
          IF (printopt.GE.1.AND.printopt.LE.10)
     .      WRITE(6,'(5X,A,I3,A,I4,A,I4,A)')
     .        'TIMER: ERROR from CheckCell (NR1 = ',
     .        NR1,'  IPOLGN =',ipolgn,'  ERR = ',chk1,')'

          STOP 'Error from CHECKCELL in TIMER'
        ENDIF
      ENDIF

c      IF (GRIDOPT.EQ.1.AND.MRSURF.EQ.MRLAST.AND.NLSFRX.EQ.1) THEN
c Change this so that it only resets rings that are outside the
c last unbroken ring...
c        IF (printopt.GE.1.AND.printopt.LE.10)
c     +  WRITE(6,'(5X,A,2I4)')
c     +    'TIMER: Resetting TIMINT (MRSURF,MRLAST) ',mrsurf,mrlast
c        DO J= MRSURF+1,NR1ST
c          TIMINT(J)=0.0
c        ENDDO
c      ENDIF

      IF (GRIDOPT.NE.0.AND.PT.EQ.1.D30.AND.ICASE.EQ.3) THEN


        I00 = PolPos(NRCELL,X00,Y00,cp,'TIMER 2')
c
c      Only look in grid cut region:
c
        DO icut = 1, npplg
          IF (i00.GE.npoint(1,icut).AND.i00.LE.npoint(2,icut)-1)
     .      GOTO 6005
        ENDDO

        WRITE(0,*) 'ERROR (TimeR): Cannot find cut region ',
     .    '( NRCELL NPANU ',nrcell,npanu,')'

 6005   CONTINUE

c
c      Can't find where to go:
c
        TMIN   = 1.D30
        IPMIN  = 0
        DELTAX = 1000.0 * VELX
        DELTAY = 1000.0 * VELY


        cpmin = 1.0D+30
        DO ii = 1, 4
          IF (cp(ii).LT.cpmin) THEN
            cpmin  = cp(ii)
            icpmin = ii
          ENDIF
        ENDDO

        IF (printopt.GE.1.AND.printopt.LE.10) THEN
          WRITE(6,'(5X,A,I4,3X,A,I4,2F10.3,2X,I4)')
     .    'TIMER: Special PolPos = ',i00,
     .    ' (N X00,Y00 ICPMIN) ',nrcell,x00,y00,icpmin
          WRITE(6,'(10X,4F14.8)') cp(1),cp(2),cp(3),cp(4)
        ENDIF

c      Particle is not on bad side... replace with separatrix...
c This isn't fool proff becasue the particle could have entered the
c cell through a poloidal surface... fix...

          IF (printopt.GE.1.AND.printopt.LE.10) THEN
            WRITE(6,'(15X,A,3E10.3,I4)')
     .        'Looking for a surface to cross...'
           ENDIF


            DO K = NPOINT(1,icut), NPOINT(2,icut)-1

              IF (rvrtag(nrcell,k).EQ.1) GOTO 6200

              ISTAT = CalcInter(X00,Y00,X00+DELTAX,Y00+DELTAY,
     .            XVERT(NRCELL,K  ,2),YVERT(NRCELL,K  ,2),
     .            XVERT(NRCELL,K+1,2),YVERT(NRCELL,K+1,2),T1,T2)

              IF ((ISTAT.EQ.1.AND.T1.LT.TMIN).AND.
     .            (ICPMIN.NE.2.OR.I00.NE.K)) THEN
                TMIN  = T1
                IPMIN = K
                ISIDE = 2
                IF (printopt.GE.1.AND.printopt.LE.10) THEN
                  WRITE(6,'(15X,A,3E10.3,I4)')
     .              'Side 2 (T1,T2,TMIN,IPMIN) ',T1,T2,TMIN,IPMIN
                ENDIF
              ELSE
                IF (printopt.EQ.3)
     .            WRITE(6,'(15X,A,3E10.3,I4)')
     .              'Side 2 (T1,T2,TMIN,K) ',T1,T2,TMIN,K
              ENDIF

              ISTAT = CalcInter(X00,Y00,X00+DELTAX,Y00+DELTAY,
     .            XVERT(NRCELL,K  ,1),YVERT(NRCELL,K  ,1),
     .            XVERT(NRCELL,K+1,1),YVERT(NRCELL,K+1,1),T1,T2)

              IF ((ISTAT.EQ.1.AND.T1.LT.TMIN).AND.
     .            (ICPMIN.NE.4.OR.I00.NE.K)) THEN
                TMIN  = T1
                IPMIN = K
                ISIDE = 4
                IF (printopt.GE.1.AND.printopt.LE.10) THEN
                  WRITE(6,'(15X,A,3E10.3,I4)')
     .              'Side 4 (T1,T2,TMIN,IPMIN) ',T1,T2,TMIN,IPMIN
                ENDIF
              ELSE
                IF (printopt.EQ.3)
     .             WRITE(6,'(15X,A,3E10.3,I4)')
     .              'Side 4 (T1,T2,TMIN,K) ',T1,T2,TMIN,K
              ENDIF




 6200      ENDDO



        IF (TMIN.LT.1.D30) THEN
          XNEW = X00 + TMIN * DELTAX
          YNEW = Y00 + TMIN * DELTAY

          IPOLGN = PolPos(NRCELL,XNEW,YNEW,cp,'TIMER 3')

          IF (ABS(VELX).GT.ABS(VELY)) THEN
            PT = (XNEW - X0) / (VELX + EPS60)
          ELSE
            PT = (YNEW - Y0) / (VELY + EPS60)
          ENDIF

          IF (printopt.GE.1.AND.printopt.LE.10) THEN
            WRITE(6,'(5X,A,I4,3X,A,2I4,2F10.3)')
     .        'TIMER: Super special PolPos = ',ipolgn,
     .        ' (N ISIDE XNEW,YNEW) ',nrcell,iside,xnew,ynew
            WRITE(6,'(10X,4F14.8)') cp(1),cp(2),cp(3),cp(4)
          ENDIF

          IF (ISIDE.EQ.2) THEN
            mrsurf = nrcell + 1
            nincx  = +1
          ELSEIF (ISIDE.EQ.4) THEN
            WRITE(6,'(5X,A,I4,3X,A,2I4,2F10.3)')
     .        'TIMER: Invalid?  PolPos = ',ipolgn,
     .        ' (N I00 ISIDE XNEW,YNEW) ',nrcell,I00,iside,xnew,ynew
            WRITE(6,'(10X,4F14.8)') cp(1),cp(2),cp(3),cp(4)

            mrsurf = nrcell
            nincx  = -1
          ELSE
            WRITE(0,*) 'NPANU = ',npanu
            STOP 'CODE: Alpha'
          ENDIF
        ELSE
          IF (printopt.GE.1.AND.printopt.LE.10)
     .      WRITE(6,'(5X,A)')
     .        'TIMER: Leaving though end of ring...?'
        ENDIF
      ENDIF

      RETURN
      END




c
c =======
c
c subroutine: CalcInter
c
c
      INTEGER FUNCTION CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,t1,t2)

      IMPLICIT none

c     Input:
      DOUBLE PRECISION a1,a2,b1,b2,c1,c2,d1,d2

c     Output:
      DOUBLE PRECISION t1,t2

      DOUBLE PRECISION TOL
      PARAMETER (TOL = 1.0E-10)

      DOUBLE PRECISION t0,e1,e2,f1,f2,fact
c
c
c
      IF (ABS(c1-d1).LT.TOL.AND.ABS(c2-d2).LT.TOL) THEN
        WRITE(0,*) 'ERROR (CalcInter): Invalid line segment'
        GOTO 9000
      ENDIF
c
c     Find projection of AB onto CD:
c
      t0 = ((a1 - c1) * (d1 - c1) + (a2 - c2) * (d2 - c2)) /
     .     ((c1 - d1)**2 + (c2 - d2)**2)

      e1 = c1 + t0 * (d1 - c1)
      e2 = c2 + t0 * (d2 - c2)
c
c     Calculate F, the intersection point between AB and CD:
c
c
c
c
      fact = (e1 - a1) * (b1 - a1) + (e2 - a2) * (b2 - a2)

      IF (fact.EQ.0.0) THEN
        t1 = 0.0
      ELSE
        t1 = ((a1 - e1)**2 + (a2 - e2)**2) / fact
      ENDIF
c
c       Determine the parametric location of F on CD:
c
      f1 = a1 + t1 * (b1 - a1)
      f2 = a2 + t1 * (b2 - a2)

      IF (abs(d1-c1).GT.abs(d2-c2)) THEN
        t2 = (f1 - c1) / (d1 - c1)
      ELSE
        t2 = (f2 - c2) / (d2 - c2)
      ENDIF

      IF (t1+TOL.GE.0.0.AND.t1-TOL.LE.1.0.AND.
     .    t2+TOL.GE.0.0.AND.t2-TOL.LE.1.0) THEN
        CalcInter = 1
      ELSE
        CalcInter = 0
      ENDIF

      RETURN
c
c     Error output:
c
9000  CONTINUE
      WRITE(6,'(5X,A,2F13.10,A,2F13.10)')
     .  ' A1,A2 = ',a1,a2,' B1,B2 = ',b1,b2
      WRITE(6,'(5X,A,2F13.10,A,2F13.10)')
     .  ' C1,C2 = ',c1,c2,' D1,D2 = ',d1,d2
      WRITE(6,'(5X,A,2F13.10,A,2F13.10)')
     .  ' E1,E2 = ',e1,e2,' F1,F2 = ',f1,f2
      WRITE(6,'(5X,A,E10.3,A,3E10.3)')
     .  ' FACT  = ',fact,' T0,T1,T2 = ',t0,t1,t2
      STOP
      END
c
c ===
c
      SUBROUTINE FindPoloidalCell(nr,x,y)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'CUPD'
      INCLUDE 'CGEOM'
      INCLUDE 'CPOLYG'

      INTEGER PolPos,CheckCell,ip1,ip2,ii
      DIMENSION cp(4)

        ip1 = PolPos(NR,X,Y,cp,'FOLNEUT 2')

        IF (printopt.GE.1.AND.printopt.LE.10) THEN
          WRITE(6,'(3X,A,I4,3X,A,I4,2F10.3)')
     .      'POLCELL: PolPos = ',ip1,
     .      ' (NR X,Y) ',nr,x,y
          WRITE(6,'(10X,4F14.8)') cp(1),cp(2),cp(3),cp(4)
        ENDIF



c        ip2 = PolPos(nr-nincx,x,Y,cp,'FOLNEUT 3')
c
c        IF (printopt.GE.1.AND.printopt.LE.10) THEN
c          WRITE(6,'(3X,A,I4,3X,A,I4,2F10.3)')
c     .      'FOLNEUT: PolPos = ',ip2,
c     .      ' (NR-NINCX x,Y) ',nr-nincx,x,y
c            WRITE(6,'(10X,4F14.8)') cp(1),cp(2),cp(3),cp(4)
c        ENDIF

        IF (ip1.NE.npcell) THEN
          IF (ip1.GT.0) THEN
            npcell = ip1
            ipolgn = ip1
          ELSEIF (ip1.LT.0) THEN
c            WRITE(0,*) 'ERROR (FindPolCell): Unable to find cell',
c     .        ' ( NR = ',nr,')'
            WRITE(6,*) 'ERROR (FindPolCell): Unable to find cell',
     .        ' ( NR = ',nr,')'

            ip1 = PolPos(nr,x,Y,cp,'FOLNEUT 4')

            itemp1 = printopt
            printopt = 99

c            DO ii1 = 1, npplg
c              DO ii2 = npoint(1,ii1),npoint(2,ii1)-1
c                ip1 = CheckCell(nr-nincx,ii2,x,y,cp)
c              ENDDO
c            ENDDO

            DO ii1 = 1, npplg
              DO ii2 = npoint(1,ii1),npoint(2,ii1)-1
                ip1 = CheckCell(nr,ii2,x,y,cp)
              ENDDO
            ENDDO

            printopt = itemp1

c            WRITE(0,*) 'NPANU = ',npanu
            WRITE(6,*) 'NPANU = ',npanu

            WRITE(6,*) 'PolPos routine - error'
            npcell = -1

c
c
c Can't use LEARC1 at the moment, because there is no differention
c of surface geometry based on which side the particle is on (or coming
c from or going to).  It is okay for particles outside the grid,
c becasue there is crude decision making in LEARC1 based on the ring index.
c
c
c            isl1 = -1
c            isl2 = -1
c            isl1 = LEARC1(x,y,0.0D0,isl2,nr,nr,
c     .                    .TRUE.,.FALSE.,
c     .                    npanu,'TEST 02     ')
c            WRITE(6,'(3X,A,I5)') 'FINDPOLOIDALCELL: LEARC1 = ',isl2
c            npcell = isl2

            RETURN
          ELSE
            WRITE(6,'(3X,A)') 'POLCELL:'
            WRITE(6,'(3X,A)') 'POLCELL: ERROR'
            WRITE(6,'(3X,A)') 'POLCELL:'
          ENDIF
        ENDIF

        IF (printopt.GE.1.AND.printopt.LE.10) THEN
          WRITE(6,'(3X,A)')
     +      'POLCELL: Output (NR,NP IR,IP MRSURF X0,Y0,Z0 X,Y,Z00)'
          WRITE(6,'(5X,A,2I4,A,2I4,A,I4,6F10.4)')
     +      '(',NR,NPCELL,') (',IRCELL,IPCELL,')',MRSURF,
     +      X0,Y0,Z0,x,y,Z00
        ENDIF

c        WRITE(6,'(10X,A,I4,2X,2I4,2X,2I4,2X,2F12.6)')
c     .    'NR NPCELL,IPCELL IP1,IP2 x,Y',
c     .    nr,npcell,ipcell,ip1,ip2,x,y

c      DO ii = 1, npplg
c        DO ip = npoint(1,ii), npoint(2,ii)
c          WRITE(6,'(3X,A,I4,3X,A,2I4,2F10.3)')
c     .      'FOLNEUT: CheckCell = ',CheckCell(nr,ipx00,y,cp),
c     .      ' (NR,IP X,Y) ',nr,ip,x,y
c          WRITE(6,'(10X,4F12.6)') cp(1),cp(2),cp(3),cp(4)
c         ENDDO
c      ENDDO


c        ip1=LEARC2(X00,Y00,NR,NPANU,'Check   1    ')
c        Ip2=LEARC2(X00,Y00,NR-nincx,NPANU,'Check   2    ')

c        WRITE(6,'(10X,A,I4,2X,2I4,2X,2I4,2X,2F12.6)')
c     .    'NRCELL NPCELL,IPCELL IP1,IP2 X00,Y00',
c     .    nrcell,npcell,ipcell,ip1,ip2,x00,y00




      RETURN
      END
c
c
c
c
c
      INTEGER FUNCTION CheckCell(nr,np,xval,yval,cp)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
c Input:
c
      INTEGER          nr,np
      DOUBLE PRECISION xval,yval,cp(4)

      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'CGRID'
      INCLUDE 'CGEOM'
      INCLUDE 'CLOGAU'
      INCLUDE 'CPOLYG'
      INCLUDE 'CCONA'
      INCLUDE 'CLGIN'
      INCLUDE 'CUPD'
      INCLUDE 'CTRIG'

      REAL*8 TOL
      PARAMETER (TOL = -1.0D-03)

c      DOUBLE PRECISION vx(4),vy(4),ax,ay,bx,by,x,y
c sl*16
      REAL*8 vx(4),vy(4),ax,ay,bx,by,x,y,tx(4),ty(4)
      INTEGER          ii,i1,i2,icnt

      CheckCell = 0

      DO ii = 1,4
        cp(ii) = -2.0
      ENDDO
c
c
c
      IF (nr.LT.1.OR.nr.GE.nr1st) THEN
        CheckCell = -1
        RETURN
      ENDIF
c
c
c
      IF (np.LT.NPOINT(1,1).OR.np.EQ.NPOINT(2,1).OR.
     .    np.EQ.NPOINT(2,2).OR.np.GT.NPOINT(2,3)-1) THEN
        CheckCell = -1
        RETURN
      ENDIF
c
c
c
      IF (gridopt.EQ.1.AND.rvrtag(nr,np).EQ.1) RETURN

      IF (printopt.EQ.99)
     .   WRITE(6,'(A,2I4)') 'CHECKCELL: (NR,NP)',
     .    nr,np

      x = xval
      y = yval

      IF (gridopt.EQ.1) THEN
        vx(1) = xvert(nr,np  ,1)
        vx(2) = xvert(nr,np  ,2)
        vx(3) = xvert(nr,np+1,2)
        vx(4) = xvert(nr,np+1,1)

        vy(1) = yvert(nr,np  ,1)
        vy(2) = yvert(nr,np  ,2)
        vy(3) = yvert(nr,np+1,2)
        vy(4) = yvert(nr,np+1,1)
      ELSE
        STOP 'ERROR (CheckCell): Unsupported geometry option'
      ENDIF

      DO ii = 1,4
        cp(ii) = -1.0
      ENDDO

      DO ii = 1, 4
        IF (ii.EQ.4) THEN
          i1 = 4
          i2 = 1
        ELSE
          i1 = ii
          i2 = ii + 1
        ENDIF

        ax = x - vx(i1)
        ay = y - vy(i1)

        bx = vx(i2) - vx(i1)
        by = vy(i2) - vy(i1)

        cp(ii) = ax * by - ay * bx

        IF (bx.NE.0.0) THEN
          tx(ii) = (x - vx(i1)) / bx
        ELSE
          tx(ii) = 0.0
        ENDIF

        IF (by.NE.0.0) THEN
          ty(ii) = (y - vy(i1)) / by
        ELSE
          ty(ii) = 0.0
        ENDIF

        IF (printopt.EQ.99)
     .    WRITE(6,'(3X,3I2,2E10.3,6E10.3,2X,E11.4)')
     .      ii,i1,i2,vx(ii),vy(ii),x,y,ax,ay,bx,by,cp(ii)
      ENDDO

      CheckCell = 1

      icnt = 0

      DO ii = 1,4
        IF (cp(ii).LT.TOL) CheckCell = 0
        IF (cp(ii).LT.0.0) icnt = icnt + 1
      ENDDO

      IF (icnt.GT.1) CheckCell = 0

      IF (CheckCell.EQ.0.AND.icnt.EQ.1) THEN
        icnt = 0
        DO ii = 1, 4
          IF (tx(ii).GE.0.0.AND.tx(ii).LE.1.0.AND.
c...BUG FIX: I DUNNO WHAT HAPPENED HERE.  THIS LINE SHOULD NEVER HAVE BEEN
c   .TRUE. BEFORE SINCE TOL WAS LESS THAN 0 AND ABS( ) COULD NOT BE.  
c   PERHAPS THE SINGLE PRECISION CALL WAS GENERATING SOME STRANGE LARGE
c   NEGATIVE NUMER ALL THE TIME? - SL, NOV 19, 2003
     .        DABS(tx(ii)-ty(ii)).LT.DABS(TOL))
c     .        ABS(tx(ii)-ty(ii)).LT.TOL)
     .      icnt = icnt + 1

          IF (printopt.EQ.99)
     .      WRITE(6,'(A,2F15.7,I6)')
     .        'CHECKCELL: ',tx(ii),ty(ii),icnt
        ENDDO

c Not sure if this is a good plan or not... this should not be necessary...
c I think the problem is in DIVIMP and in the grid refinement code...
        IF (icnt.EQ.1) THEN
          CheckCell = 1

          IF (printopt.EQ.99)
     .      WRITE(6,'(A,2I5)')
     .        'CHECKCELL: Found particle on surface ',np,nr
        ENDIF
      ENDIF

      RETURN
      END
c
c
c
      INTEGER FUNCTION PolPos(nr,x,y,cp,callid)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
c Input:
c
      INTEGER          nr
      DOUBLE PRECISION x,y,cp(4)
      CHARACTER callid*(*)

      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'CGRID'
      INCLUDE 'CGEOM'
      INCLUDE 'CLOGAU'
      INCLUDE 'CPOLYG'
      INCLUDE 'CCONA'
      INCLUDE 'CLGIN'
      INCLUDE 'CUPD'
      INCLUDE 'CTRIG'

      DOUBLE PRECISION TOL
      PARAMETER (TOL = -1.0D-05)

      INTEGER CheckCell

      INTEGER          ii,i2,ip,status,cpminip(200),chkcell
      DOUBLE PRECISION tx,ty,x1,y1,cp1(4),cpmin(200),cpmax
c
c Initialization:
c
      status = 0
      PolPos = -1
c
c Check if NR is outside standard grid:
c
      IF (nr.LT.1.OR.nr.GE.nr1st) THEN
        PolPos = 0
        RETURN
      ENDIF

      cp(1) = -1.0
      cp(2) = -1.0
      cp(3) = -1.0
      cp(4) = -1.0

      DO ii = 1, npplg
        DO ip = npoint(1,ii), npoint(2,ii)-1

           chkcell = CheckCell(nr,ip,x,y,cp1)

           IF (chkcell.EQ.1) THEN
             status = status + 1

             cp(1)  = cp1(1)
             cp(2)  = cp1(2)
             cp(3)  = cp1(3)
             cp(4)  = cp1(4)

             cpmin(status) = 1.0
             DO i2 = 1, 4
               cpmin(status) = MIN(cpmin(status),cp(i2))
             ENDDO
             cpminip(status) = ip

             PolPos = ip
           ELSEIF (chkcell.EQ.-1) THEN
             WRITE(0,*) 'ERROR (CheckCell): Violation ( NRCELL = ',
     .         nr,'  IP = ',ip,'  id = ',callid,')'
             STOP
           ENDIF
         ENDDO
      ENDDO

      IF (status.EQ.1) THEN
        IF (cpmin(1).LT.-1.0D-8) THEN
          IF (printopt.GE.1.AND.printopt.LE.10)
     .      WRITE(6,'(3X,A,1P,E14.7)')
     .        ' POLPOS: *** Questionable *** (CPMAX) ',cpmin(1)
        ENDIF
      ELSEIF (status.LT.10) THEN
        cpmax = -1.0
        DO ii = 1, status
          IF (cpmin(ii).GT.cpmax) THEN
            cpmax  = cpmin  (ii)
            PolPos = cpminip(ii)
          ENDIF

          IF (printopt.GE.1.AND.printopt.LE.10)
     .      WRITE(6,'(3X,A,I2,A,I3,A,E14.7,3A)')
     .        ' POLPOS: Ambiguity (NRCELL = ',
     .        nr,'  NPCELL = ',cpminip(ii),'  CPMIN = ',cpmin(ii),
     .        '  id = ',callid,')'
        ENDDO

        IF (cpmax.LT.TOL)
c        IF (printopt.GE.1.AND.printopt.LE.10.AND.cpmax.LT.-1.0D-8)
     .    WRITE(6,'(3X,A,2I5,1P,E15.7)')
     .      ' POLPOS: Ambiquity POLPOS,NR CPMAX =  ',PolPos,nr,cpmax


      ELSEIF (status.GE.10) THEN

         WRITE(6,'(3X,A)')
     +    ' POLPOS: Output (NR,NP IR,IP MRSURF X0,Y0,Z0 X00,Y00,Z00)'
         WRITE(6,'(5X,A,2I4,A,2I4,A,I4,6F10.4,1X,A)')
     +    '(',NRCELL,NPCELL,') (',IRCELL,IPCELL,')',MRSURF,
     +    X0,Y0,Z0,X00,Y00,Z00,callid

c
c
c
        status = 0
        printopt = 99

        DO ii = 1, npplg
          DO ip = npoint(1,ii), npoint(2,ii)-1
             chkcell = CheckCell(nr,ip,x,y,cp1)
             IF (chkcell.EQ.1) WRITE(6,*) '   MATCH'
           ENDDO
        ENDDO

c
c
c



        WRITE(0,*) 'NPANU = ',npanu
        STOP 'ERROR (PolPos): More than 3 cells identified'
      ENDIF

      RETURN
      END
c
c
c
c
c Jul 22 - There is going to be an additional probelm here, since the
c the assumption that the outer ring of the grid is continuous is being
c violated.  I may not have to worry about this is I can get the
c input file correct.  Still - be careful and add a line here to wanr about wall
c quantities.
c
c
c
      SUBROUTINE CHECKNUM(TAG,I1,I2,X1,X2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      CHARACTER TAG*(*)

      INCLUDE 'PARMMOD'

c      IF (printopt.GE.1.AND.printopt.LE.10) THEN
c        WRITE(6,'(10X,2A,2I4,2G15.7,$)')
c     +    'CHECKNUM: ',TAG,I1,I2,X1,X2
c      ENDIF

      CHKCNT = CHKCNT + 1

      IF (X1.EQ.X2) THEN

c        IF (printopt.GE.1.AND.printopt.LE.10) WRITE(6,*) ' '
c        IF (i1.EQ.25.AND.printopt.GE.1.AND.printopt.LE.10) THEN
c          WRITE(6,'(10X,A,A6,2I4,2G15.7)')
c     +      'CHECKNUM: ',TAG,I1,I2,X1,X2
c        ENDIF



      ELSE
c        IF (printopt.EQ.2.OR.debugopt.EQ.-2) THEN
c          WRITE(6,'(10X,2A,2I4,2G15.7,A)')
c     +      'CHECKNUM: ',TAG,I1,I2,X1,X2,' ERROR'
c        ENDIF
        CHKERR = CHKERR + 1
      ENDIF

      RETURN
      END
c
c
c
      SUBROUTINE SLOUT(slmess,slf1,slf2,slf3)

      CHARACTER slmess*100
      REAL      slf1,slf2,slf3

      RETURN
      END
c
c
c
      SUBROUTINE ReadData(zeile)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c put in subroutine... add a check to see if *** 0 exists...
c
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'

      CHARACTER*1024 line


      CHARACTER zeile*(*)


      geomopt  = 0
      gridopt  = 0
      addopt   = 0
      neutopt  = 0
      debugopt = 0
      gchkopt  = 0
      printopt = 0
      optuser  = 0
      optzmotion = 0

      eirntrans = 0

      optconmap = 0
c      WRITE(0,*) '**************************'
c      WRITE(0,*) 'CONNECTION MAP TURNED ON!'
c      WRITE(0,*) '**************************'
c      WRITE(0,*) 'CONNECTION MAP TURNED OFF!'
c      WRITE(0,*) '**************************'
       WRITE(6,'(1X,A72)') ZEILE


        READ (IUNIN,*) comment,geomopt
        READ (IUNIN,*) comment,gridopt
c        READ (IUNIN,*) comment,addopt
        READ (IUNIN,*) comment,neutopt
        READ (IUNIN,*) comment,debugopt
        READ (IUNIN,*) comment,cxd2opt      
        READ (IUNIN,*) comment,optuser      
        READ (IUNIN,*) comment,eirntrans
        IF (eirntrans.GT.0) THEN
          DO i1 = 1, eirntrans
            READ(IUNIN,*) (eirtrans(i1,i2),i2=1,3)
          ENDDO
        ENDIF
        READ (IUNIN,*) comment,nstrdat
        READ (IUNIN,'(6I6)') (strdat(i1),i1=1,nstrdat)        
        READ (IUNIN,*) comment,volcor2

        READ (IUNIN,'(A1024)') line
        IF (line(2:10).EQ.'Radial co') THEN
          READ (line,*) comment,optconmap
        ELSE
          WRITE(0,*) 'DIVIMP DATA LINE MISSING, STOPPING EIRENE'
          STOP
        ENDIF

       IF (optconmap.EQ.1) THEN
         WRITE(6,*) '*****************************'
         WRITE(6,*) '* CONNECTION MAP TURNED ON! *'
         WRITE(6,*) '*****************************'
       ENDIF

       IF (debugopt.GE.1.AND.debugopt.LE.20) THEN
         printopt = debugopt
       ELSEIF (debugopt.EQ.-41) THEN
         optzmotion = 1
         debugopt = 0
       ELSEIF (debugopt.EQ.-42) THEN
         optzmotion = 2
         debugopt = 0
       ELSEIF (debugopt.EQ.-43) THEN
         optzmotion = 3
         debugopt = 0
       ELSEIF (debugopt.EQ.-44) THEN
         optzmotion = 4
         debugopt = 0
       ELSEIF (debugopt.EQ.-20) THEN
         opttest=-(debugopt+19)
         debugopt=0
         WRITE(0,*) 'OPTTEST= ',opttest
       ELSE
         printopt = 0
       ENDIF

      IF (debugopt.NE.0.OR.optzmotion.NE.0) THEN
        WRITE(0,*)
        WRITE(0,*) '=================================='
        WRITE(0,*)
        WRITE(0,'(1X,A,I6)') 'grid option            ',gridopt
        WRITE(0,'(1X,A,I6)') 'geometry source option ',geomopt
        WRITE(0,'(1X,A,I6)') 'geometry check option  ',gchkopt
c        WRITE(0,'(1X,A,I6)') 'addusr option          ',addopt
        WRITE(0,'(1X,A,I6)') 'neutral wall option    ',neutopt
        WRITE(0,'(1X,A,I6)') 'debug option           ',debugopt
        WRITE(0,'(1X,A,I6)') 'print option           ',printopt
        WRITE(0,'(1X,A,I6)') 'CX D2+ production      ',cxd2opt      
        WRITE(0,'(1X,A,I6)') 'Z motion               ',optzmotion
        WRITE(0,'(1X,A,I6)')
        WRITE(0,*) '=================================='
        WRITE(0,*)
      ENDIF

c Read next 2 lines:
c       READ (IUNIN,'(A72)') ZEILE
c       WRITE(6,'(1X,A72)') ZEILE
c       READ (IUNIN,'(A72)') ZEILE


       RETURN
       END

c
c
c
      SUBROUTINE INIUSR
c      WRITE(0,*) 'WARNING: Calling substitute INIUSR routine'
      RETURN
      END
c
c
c
      SUBROUTINE SP0USR
c      WRITE(0,*) 'WARNING: Calling substitute SP0USR routine'
      RETURN
      END
c
c
c
c ======================================================================
c
c subroutine: ASDUSR
c
c Outputs additional surface data.
c
c
c
c     DIIN
c     TIIN
c     VXIN
c     VYIN
c     VZIN
c
c     PDENA Neutral atom density
c     EDENA    "     "   energy
c     PDENM Neutral molecule density
c     EDENM    "       "     energy
c     PDENI Ionised molecule density
c     EDENI    "       "     energy
c
c
c
      SUBROUTINE ASDUSR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'PARMMOD'
      INCLUDE 'COMSOU'
      INCLUDE 'COMUSR'
      INCLUDE 'COUTAU'
      INCLUDE 'CESTIM'
      INCLUDE 'CGRID'

      COMMON /RESCOM/ EIRRES
      REAL*8          EIRRES(7,NPLS)
c      COMMON /RESCOM/ RRN,RRE,RRM,RATN,RATE,RATM,RESN,RESE,RESM

      INTEGER FP,IN,COUNT

      DATA COUNT /0/
      SAVE

      COUNT = COUNT + 1

      IF (output)
     .WRITE(0,*) 'WRITING ASD'

      FP=98
      OPEN(UNIT=FP,FILE='addsur.dat',ACCESS='SEQUENTIAL',
     .     STATUS='UNKNOWN',POSITION='APPEND')

      WRITE(FP,90) '*'
      WRITE(FP,94) '* DATA FOR ADDITIONAL EIRENE SURFACES ',COUNT
      WRITE(FP,90) '*                                     '
      WRITE(FP,90) '* NRADD (NUMBER OF SURFACES)=         '
      WRITE(FP,92) NRADD

C SURFACE DATA
      WRITE(FP,90) '* SURFACE DATA '
      WRITE(FP,91) 0,'IN','VOLUME (CM-3)'
      DO IN=NSURF+1,NSURF+NRADD
        WRITE(FP,92) IN-NSURF,IN,VOL(IN)
      ENDDO
C NEUTRAL ATOM DATA

      IF (NATMI.GT.2) THEN
        WRITE(0,*) 'ONLY PASSING 1ST AND 2ND ATOMIC INDEX TO DIVIMP'
        WRITE(6,*) 'ONLY PASSING 1ST AND 2ND ATOMIC INDEX TO DIVIMP'
      ENDIF

      WRITE(FP,90) '* NEUTRAL ATOM SPECIES       '
      WRITE(FP,90) '* NATMI (NUMBER OF SPECIES)= '
      WRITE(FP,92) MIN(2,NATMI)
      DO I1=1,MIN(2,NATMI)
        WRITE(FP,91) I1,'IN','PDENA (    )','EDENA (   )'
        DO IN=NSURF+1,NSURF+NRADD
          WRITE(FP,92) IN-NSURF,IN,PDENA(I1,IN),EDENA(I1,IN)
        ENDDO
      ENDDO
C NEUTRAL MOLECULE DATA
      WRITE(FP,90) '* NEUTRAL MOLECULE SPECIES   '
      WRITE(FP,90) '* NMOLI (NUMBER OF SPECIES)= '
      WRITE(FP,92) NMOLI
      DO I1=1,NMOLI
        WRITE(FP,91) I1,'IN','PDENM (    )','EDENM (   )'
        DO IN=NSURF+1,NSURF+NRADD
c...bug - Nov 22, 1999
          WRITE(FP,92) IN-NSURF,IN,PDENM(I1,IN),EDENM(I1,IN)
c          WRITE(FP,92) IN-NSURF,IN,PDENA(I1,IN),EDENA(I1,IN)
        ENDDO
      ENDDO
C TEST ION DATA
      WRITE(FP,90) '* TEST ION SPECIES           '
      WRITE(FP,90) '* NIONI (NUMBER OF SPECIES)= '
      WRITE(FP,92) NIONI
      DO I1=1,NIONI
        WRITE(FP,91) I1,'IN','PDENI (    )','EDENI (   )'
        DO IN=NSURF+1,NSURF+NRADD
          WRITE(FP,92) IN-NSURF,IN,PDENI(I1,IN),EDENI(I1,IN)
        ENDDO
      ENDDO
C BULK ION DATA (FOR BGK ROUTINES)
      WRITE(FP,90) '* BULK ION SPECIES           '
      WRITE(FP,90) '* NPLSI (NUMBER OF SPECIES)= '
      WRITE(FP,92) NPLSI
      DO I1=1,NPLSI
        WRITE(FP,91) I1,'IN','DIIN','TIIN','VXIN','VYIN','VZIN','EDRIFT'
        DO IN=NSURF+1,NSURF+NRADD
          WRITE(FP,92) IN-NSURF,IN,DIIN(I1,IN),TIIN(I1,IN),
     .                 VXIN(I1,IN),VYIN(I1,IN),VZIN(I1,IN),EDRIFT(I1,IN)
        ENDDO
      ENDDO

c RESIDUALS CALCULATED IN MODBGK ROUTINE
      IF (NITER.GE.1) THEN
        WRITE(FP,90) 'BGK RESIDUALS'
        WRITE(FP,*) NPLSI
        DO I1 = 1, NPLSI
          WRITE(FP,93) I1,(EIRRES(I2,I1),I2=1,7)
        ENDDO
      ENDIF

C STRATA DATA
      WRITE(FP,90) 'STRATUM DATA FOR SELECTED ADDITIONAL CELLS'
      WRITE(FP,*) NSTRDAT,NSTRAI,NATMI,NMOLI
      DO ISTR=1,NSTRDAT
        I1=STRDAT(ISTR)
        WRITE(FP,*) STRDAT(ISTR)
        DO ISTRA=1,NSTRAI
          DO IATM=1,NATMI
            WRITE(FP,'(2I6,1P,2E12.4)') ISTRA,IATM, 
     .        PSTRDATA(IATM,ISTR,ISTRA),ESTRDATA(IATM,ISTR,ISTRA)
          ENDDO
          DO IMOL=1,NMOLI
            WRITE(FP,'(2I6,1P,2E12.4)') ISTRA,IMOL, 
     .        PSTRDATM(IMOL,ISTR,ISTRA),ESTRDATM(IMOL,ISTR,ISTRA)
          ENDDO
        ENDDO
c...    Totals:
        DO IATM=1,NATMI
          WRITE(FP,'(2I6,1P,2E12.4)') 0,IATM, 
     .      PDENA(IATM,NSURF+I1),EDENA(IATM,NSURF+I1)
        ENDDO
        DO IMOL=1,NMOLI
          WRITE(FP,'(2I6,1P,2E12.4)') 0,IMOL, 
     .      PDENM(IMOL,NSURF+I1),EDENM(IMOL,NSURF+I1)
        ENDDO
      ENDDO

c...  SNEAK IN THE FLUX FOR EACH STRATUM AS WELL
      WRITE(FP,'(1P,100(E12.4:))') (FLUXT(ISTRA),ISTRA=1,NSTRAI)



c...TEMP
      WRITE(6,*) 'TESTING STRATA DATA'
      DO ISTR=1,NSTRDAT
        I1=STRDAT(ISTR)
        WRITE(6,*)
        WRITE(6,*)  ISTR,I1,':'

        DO ISTRA=1,NSTRAI
          DO IATM=1,NATMI
            WRITE(6,'(2I6,1P,2E12.4)') ISTRA,IATM, 
     .        PSTRDATA(IATM,ISTR,ISTRA),
     .        ESTRDATA(IATM,ISTR,ISTRA)
          ENDDO
          DO IMOL=1,NMOLI
            WRITE(6,'(2I6,1P,2E12.4)') ISTRA,IMOL, 
     .        PSTRDATM(IMOL,ISTR,ISTRA),
     .        ESTRDATM(IMOL,ISTR,ISTRA)
          ENDDO
        ENDDO
 
        WRITE(6,*)
        DO IATM=1,NATMI
          WRITE(6,'(I6,1P,2E12.4)') IATM, 
     .      PDENA(IATM,NSURF+I1),
     .      EDENA(IATM,NSURF+I1)
        ENDDO
        DO IMOL=1,NMOLI
          WRITE(6,'(I6,1P,2E12.4)') IMOL, 
     .      PDENM(IMOL,NSURF+I1),
     .      EDENM(IMOL,NSURF+I1)
        ENDDO
      ENDDO

      CLOSE (FP)

      RETURN
90    FORMAT(A)
91    FORMAT(I3,A5,7X,20(A15:))
92    FORMAT(3X,I5,I7,1P,20(E15.7:),0P,I6)
93    FORMAT(1P,I4,7E12.4,0P)
94    FORMAT(A,I6)
99    STOP
      END
c
c ======================================================================
c
c subroutine: DSet
c
c
      SUBROUTINE DSet(array,ndim,val)

      INTEGER ndim
      DOUBLE PRECISION array(ndim),val


      INTEGER ii

      DO ii = 1, ndim
        array(ii) = val
      ENDDO

      RETURN
      END


c
c ======================================================================
c
c subroutine: LSet
c
c
      SUBROUTINE LSet(array,ndim,val)

      LOGICAL array(*),val
      INTEGER ndim

      INTEGER ii

      DO ii = 1, ndim
        array(ii) = val
      ENDDO

      RETURN
      END


      INTEGER FUNCTION CHKVAC(NEWCELL,OLDCELL,X0,Y0,Z0,ILEARN,NPANU)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'PARMMOD'
c      INCLUDE 'CLOGAU'
c      INCLUDE 'COMPRT'
c      INCLUDE 'CGRID'
c      INCLUDE 'CCONA'

      INTEGER NEWCELL,OLDCELL,NPANU
      REAL*8  X0,Y0,Z0

      INTEGER I1,I2,I3,FP,ISMART(0:NLIM,5),INIT,INDEX,ILEARN,
     .        IMATCH,CUT1,OLDCELL2,ICOUNT,LASTCL,CELLREGION
      LOGICAL NEW,OLD,OUTPUT1,SWITCH,STATUS
      REAL*8  X1,X2,TX,Y1,Y2,TY

      REAL*8     TOL1,TOL2

      DATA INIT,ICOUNT /0,0/
      SAVE

      OLDCELL2 = OLDCELL

      ICOUNT = ICOUNT + 1

      LASTCL = -1

      IF (INIT.EQ.0) THEN
        INIT = 1
        DO I1 = 1, 5
          ISMART(0,I1) = 0
        ENDDO
      ENDIF

      IF (ILEARN.LT.1.OR.ILEARN.GT.4) THEN      
        WRITE(0,*) 'LEARNING INDEX OUT OF RANGE',ILEARN
        STOP 'STOP IN CHKVAC'
      ENDIF

      IF (ISMART(0,ILEARN  ).EQ.NLIM.OR.
     .    ISMART(0,ILEARN+1).EQ.NLIM) THEN
        WRITE(0,*) 'ISMART INDEX OUT OF BOUNDS'
        STOP 'STOP IN CHKVAC'
      ENDIF

c      TOL1 = 2.0D-06
      TOL1 = 2.0D-05
c      TOL2 = 1.0D-05
      TOL2 = 1.0D-04

      NEW = .FALSE.
      OLD = .FALSE.

c  SOME OPTIMIZATION REQUIRED!

      output1 = .FALSE.
c...OUTPUT ON ERRORS IS TURNED OFF!
c      IF (optvac.EQ.2) output1 = .TRUE.

      CELLREGION=-1

      CHKVAC = -1
      IMATCH = -1

      STATUS = .TRUE.

      fp = 6

      TOR=Z0

c...  Find which toroidal block the particle is in, and adjust OLDCELL2
c     accordingly:
       IF (ASC3DMODE.EQ.2.AND.OLDCELL.GT.-1) THEN
        DO CUT1=1,ASCNCUT
          IF (TOR.GE.ASCZMIN3D(CUT1).AND.
     .        TOR.LT.ASCZMAX3D(CUT1)) THEN
            OLDCELL2 = OLDCELL2 - (CUT1 - 1) * ASCNCELL
            EXIT
          ENDIF
        ENDDO
      ELSEIF (OLDCELL.LT.-1) THEN
        OLDCELL2=-1
C RESTRICT ADDITIONAL CELL SEARCH TO
C VACUUM GRID REGION -OLDCELL+1
        CELLREGION=-OLDCELL-1
      ENDIF

c      IF (npanu.EQ.821) THEN
c        DO i1 = 1, ISMART(0,ILEARN)
c          WRITE(0,*) 'SMART:',ilearn,i1,ISMART(i1,ILEARN)
c        ENDDO
c      ENDIF


      DO WHILE(STATUS) 

        IF (output1) WRITE(fp,*) 'CHKVAC: ',tol2

        DO INDEX = -ISMART(0,ILEARN)+1, ASCNCELL

C...      Learning opportunity:
          IF (INDEX.LE.0) THEN
            I1 = ISMART(-INDEX+1,ILEARN)
c            WRITE(0,*) 'LEARNING: ',ilearn,I1,ISMART(0,ilearn)
          ELSE
            I1 = INDEX
          ENDIF

          IF (output1) WRITE(fp,'(A,4I6,2F10.5)') 
     .      'CHKVAC: ',I1,newcell,oldcell,chkvac,asczmin3d(i1),
     .      asczmax3d(i1)

          IF (CELLREGION.NE.-1.AND.ASCREGION(I1).NE.CELLREGION) CYCLE

          DO I2 = 1, ASCNVERTEX(I1)
            IF (I2.EQ.ASCNVERTEX(I1)) THEN
              I3 = 1
            ELSE
              I3 = I2 + 1
            ENDIF
            X1 = ASCVERTEX(2*I2-1,I1)
            Y1 = ASCVERTEX(2*I2  ,I1)
            X2 = ASCVERTEX(2*I3-1,I1)
            Y2 = ASCVERTEX(2*I3  ,I1)

c            IF (i1.EQ.ASCNCELL.and.ilearn.eq.4) THEN
c              WRITE(0,'(A,4F12.4)') '-SEARCH->',x1,y1,x2,y2
c            ENDIF

            SWITCH = .FALSE.
            IF (DABS(X1-X2).LT.1.0D-08) THEN
              IF (DABS(X0-X1).LT.TOL1) THEN
                TX = +99.0D0
              ELSE
                TX = -99.0D0
              ENDIF
            ELSE
              TX = (X0 - X1) / (X2 - X1)
              IF (TX.LT.0.25) THEN
                TX = (X0 - X2) / (X1 - X2)
                SWITCH = .TRUE.
              ENDIF
            ENDIF
            IF (DABS(Y1-Y2).LT.1.0D-08) THEN
c            IF (Y1.EQ.Y2) THEN
              IF (DABS(Y0-Y1).LT.TOL1) THEN
                TY = +99.0D0
              ELSE
                TY = -99.0D0
              ENDIF
            ELSE
              IF (SWITCH) THEN
                TY = (Y0 - Y2) / (Y1 - Y2)
              ELSE
                TY = (Y0 - Y1) / (Y2 - Y1)
              ENDIF
            ENDIF

            IF (((DABS(TX-TY).LT.TOL2.AND.
     .           TX.GE.0.0D0.AND.TX.LE.1.0D0.AND.
     .           TY.GE.0.0D0.AND.TY.LE.1.0D0).OR.
     .          (TX.EQ.+99.0D0.AND.TY.GE.0.0D0.AND.TY.LE.1.0D0).OR.
     .          (TY.EQ.+99.0D0.AND.TX.GE.0.0D0.AND.TX.LE.1.0D0)).AND.
     .          (ASC3DMODE.EQ.0.OR.ASC3DMODE.EQ.2.OR.
     .           (TOR.GT.ASCZMIN3D(I1).AND.TOR.LT.ASCZMAX3D(I1)))) THEN



              IF     (NEWCELL.EQ. 0) THEN
c...            Find additional cell that particle is in:
                CHKVAC = ASCCELL(I1)   
                IMATCH = I1
c...bug?
              ELSEIF (NEWCELL.EQ.-1.AND.OLDCELL2.NE.ASCCELL(I1)) THEN
c              ELSEIF (NEWCELL.EQ.-1.AND.OLDCELL.NE.ASCCELL(I1)) THEN
c...            Find additional cell that particle is in:
                IF     (CHKVAC.EQ.-1) THEN
                  CHKVAC = ASCCELL(I1)   
                  IMATCH = I1
                  LASTCL = CHKVAC
                ELSEIF (ASCCELL(I1).NE.LASTCL) THEN
c                ELSEIF (CHKVAC.NE.ASCCELL(I1)) THEN
                  CHKVAC = -2
                  IMATCH = -2
                ENDIF
              ENDIF


c              IF (npanu.EQ.821) THEN
c                WRITE(0,'(A,5I6,4F10.6,5I6,F10.7)') 
c     .            '-->',npanu,oldcell-1,oldcell2-1,asccell(i1)-1,
c     .                  i2,x0,x1,x2,TOR,icount,chkvac,lastcl,index,i1,
c     .                  tol2
c              ENDIF

c...          Find which toroidal block the particle is in, and increment
c             cell index accordingly:
c...bug?
              IF (ASC3DMODE.EQ.2.AND.CHKVAC.GE.1.AND.
     .            OLDCELL2.NE.ASCCELL(I1)) THEN
c              IF (ASC3DMODE.EQ.2.AND.CHKVAC.GE.1) THEN
                DO CUT1=1,ASCNCUT
                  IF (TOR.GE.ASCZMIN3D(CUT1).AND.
     .                TOR.LT.ASCZMAX3D(CUT1)) THEN
                    CHKVAC = CHKVAC + ASCNCELL * (CUT1 - 1)
                    EXIT
                  ENDIF
                ENDDO
              ENDIF

            ENDIF

            IF (output1) 
     .        WRITE(fp,'(A,2I6,2F12.6,2X,3F10.5,I6,4F12.6)') 
     .          '   ->',I1,I2,TX,TY,X0,Y0,TOR,CHKVAC,X1,X2,Y1,Y2

          ENDDO

c...      Forced exit if using learning curve:
c          IF (INDEX.LE.0) WRITE(0,*) 'TRYING:',i1,index,imatch

          IF (INDEX.LE.0.AND.CHKVAC.GT.0) THEN
c            IF (npanu.EQ.821) WRITE(0,*) 'TRYING TO BOOT',INDEX
            EXIT
          ENDIF

        ENDDO

c       IF (npanu.EQ.821) WRITE(0,*) 'BOOTED',index,i1

        IF (CHKVAC.LT.0.AND.TOL2.LE.1.0D-02) THEN
          TOL2 = TOL2 * 10.0D0
        ELSE

c...      Add cell to learning curve:
          IF (INDEX.GT.0.AND.imatch.GT.0) THEN
c...        About 10% performance improvement...or so...
            ISMART(0               ,ILEARN) = ISMART(0,ILEARN) + 1
            ISMART(ISMART(0,ILEARN),ILEARN) = imatch
c            WRITE(0,*) 'NEW:',INDEX,imatch
          ENDIF

          STATUS = .FALSE.
        ENDIF

      ENDDO



c...  Trying to save the day by looking to see if the particle was launched
c     inside an additional cell, rather than on its boundary.  This can 
c     happen if a polygon surface is not aligned perfectly with the
c     the vacuum grid:


      IF (ILEARN.EQ.4.AND.CHKVAC.LT.0) THEN      



        CHKVAC = -1
        DO INDEX = -ISMART(0,ILEARN+1)+1, ASCNCELL
          IF (INDEX.LE.0) THEN
            I1 = ISMART(-INDEX+1,ILEARN+1)
          ELSE
            I1 = INDEX
          ENDIF
          IF (CELLREGION.NE.-1.AND.ASCREGION(I1).NE.CELLREGION) CYCLE
          XMIN =  1.0D+10
          XMAX = -1.0D+10
          YMIN =  1.0D+10
          YMAX = -1.0D+10
          DO I2 = 1, ASCNVERTEX(I1)
            XMIN = MIN(XMIN,ASCVERTEX(2*I2-1,I1))
            XMAX = MAX(XMAX,ASCVERTEX(2*I2-1,I1))
            YMIN = MIN(YMIN,ASCVERTEX(2*I2  ,I1))
            YMAX = MAX(YMAX,ASCVERTEX(2*I2  ,I1))
          ENDDO

          IF (X0.GT.XMIN.AND.X0.LT.XMAX.AND.
     .        Y0.GT.YMIN.AND.Y0.LT.YMAX) THEN
            IF (CHKVAC.EQ.-1) THEN
              CHKVAC = ASCCELL(I1)

              IF (INDEX.GT.0) THEN
                ISMART(0                 ,ILEARN+1)=ISMART(0,ILEARN+1)+1
                ISMART(ISMART(0,ILEARN+1),ILEARN+1)=I1
              ENDIF
            ELSE
              CHKVAC = -2
            ENDIF
          ENDIF
c...logic is flawed here:
          IF (INDEX.LE.0.AND.CHKVAC.GE.1) EXIT
        ENDDO
        IF (ASC3DMODE.EQ.2.AND.CHKVAC.GE.1) THEN
          DO CUT1=1,ASCNCUT
            IF (TOR.GE.ASCZMIN3D(CUT1).AND.
     .          TOR.LT.ASCZMAX3D(CUT1)) THEN
              CHKVAC = CHKVAC + ASCNCELL * (CUT1 - 1)
              EXIT
            ENDIF
          ENDDO
        ENDIF

      ENDIF



      RETURN
99    STOP
      END



      LOGICAL FUNCTION CHKTRA(Z02,MSURF2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'CLOGAU'
      INCLUDE 'CGRID'
      INCLUDE 'CCONA'

      IF (NLTRA) THEN
        IF (NLTOR) THEN
          STOP 'NLTOR=.TRUE. NOT SUPPORTED IN CHKVAC'
        ELSE
c...      Set the toroidal direction comparison variable to PHI instead
c         of Z0:
          PHISEG=2.0D0*PIA/DBLE(NTTRA-1)
          Z03=(PHI+DBLE(NTRSEG-1)*PHISEG)/DEGRAD
c...      Special case for (almost) full toroidal grid:
          IF (NTRSEG.EQ.1.AND.Z03.LT.0.0) Z03 = Z03 + 360.0
        ENDIF
      ELSE
        Z03 = Z02
      ENDIF

      CHKTRA=.FALSE.

      I1=0
      DO WHILE(.NOT.CHKTRA.AND.I1.LT.EIRNTRANS)
        I1=I1+1
        IF (MSURF2.EQ.NLIM+IDNINT(EIRTRANS(I1,1)).OR.
     .      MSURF2.EQ.-1) THEN
          IF (Z03.GT.EIRTRANS(I1,2).AND.
     .        Z03.LT.EIRTRANS(I1,3)) THEN
            CHKTRA=.TRUE.
          ENDIF
        ENDIF

c          WRITE(0,'(A,3F12.4,I6,L6,2I6)') 
c     .      '-->',z03,EIRTRANS(I1,2),EIRTRANS(I1,3),
c     .      MSURF2,CHKTRA,I1,EIRNTRANS

      ENDDO
      RETURN
      END

      SUBROUTINE CHKSTD(Z02,CELL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'CLOGAU'
      INCLUDE 'CGRID'
      INCLUDE 'CCONA'

c...  Remove NPANU2
      INTEGER CELL,NPANU2

      IF (NLTRA) THEN
        IF (NLTOR) THEN
          STOP 'NLTRA=.TRUE. AND NLTOR=.TRUE. NOT SUPPORTED IN CHKVAC'
        ELSE
c...      Set the toroidal direction comparison variable to PHI instead
c         of Z0:
          PHISEG=2.0D0*PIA/DBLE(NTTRA-1)
          Z03=(PHI+DBLE(NTRSEG-1)*PHISEG)/DEGRAD
c...      Special case for (almost) full toroidal grid:
          IF (NTRSEG.EQ.1.AND.Z03.LT.0.0) Z03 = Z03 + 360.0D0
c          WRITE(0,*) 'SUPERPHI=',Z03,ntrseg
c          WRITE(0,*) 'SUPERPHI=',phi,degrad,phiseg
        ENDIF
      ELSE
        Z03=Z02
      ENDIF

      DO CELL=1,EIRNSDTOR-1
        IF (Z03.GE.EIRSDTOR(CELL).AND.Z03.LT.EIRSDTOR(CELL+1)) EXIT
      ENDDO

      IF (CELL.EQ.EIRNSDTOR) THEN
        IF (.FALSE..AND.
     .      NLTRA.AND.NTRSEG.EQ.1.AND.Z03.GT.-0.5*PHISEG/DEGRAD) THEN
c...      Special case that occurs for a full toroidal grid.  Should really have
c         some better checks here, to make sure that the full grid is in use:
          CELL = 1
        ELSE
          WRITE(6,*) 'CANNOT FIND TOROIDAL REGION ON STANDARD GRID',NPANU
          WRITE(0,*) 'CANNOT FIND TOROIDAL REGION ON STANDARD GRID',NPANU
          WRITE(0,*) 'Z0=',Z03,NTRSEG,-0.5*PHISEG/DEGRAD
          DO CELL=1,EIRNSDTOR-1
            WRITE(0,*) EIRSDTOR(CELL),EIRSDTOR(CELL+1)
          ENDDO
        ENDIF

        CELL=0
      ENDIF

     
      RETURN
99    STOP
      END


c slmod begin


c
c
c
c
c
c
      SUBROUTINE BGKUSR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'COMUSR'
      INCLUDE 'COMPRT'
      INCLUDE 'COMNNL'
      INCLUDE 'COMSOU'
      INCLUDE 'CGRID'
      INCLUDE 'CGEOM'
      INCLUDE 'CINIT'
      INCLUDE 'COMXS'
      INCLUDE 'CSPEI'
      INCLUDE 'CTRCEI'
      INCLUDE 'CCONA'
      INCLUDE 'CESTIM'
      INCLUDE 'CSDVI'
      INCLUDE 'CSDVI_BGK'
      INCLUDE 'CSDVI_COP'
      INCLUDE 'COUTAU'
      INCLUDE 'CLOGAU'
      INCLUDE 'CZT1'
C
      DOUBLE PRECISION, ALLOCATABLE :: PDEN(:), EDEN(:),
     .          PDEN2(:), EDEN2(:), ENERGY(:,:)
      DIMENSION ITYP1(NPLS),ITYP2(NPLS),ISPZ1(NPLS),ISPZ2(NPLS),
     .          IREL1(NPLS),INRC1(NPLS)
      LOGICAL LMARK(NPLS)
      LOGICAL TRCSAV

      NBLCKA=0
      ALLOCATE (PDEN(NRAD))
      ALLOCATE (EDEN(NRAD))
      ALLOCATE (PDEN2(NRAD))
      ALLOCATE (EDEN2(NRAD))
      ALLOCATE (ENERGY(NPLS,NRAD))

      DO 1000 IPLS=1,NPLSI
C .........................................................................
C
C  IS IPLS AN ARTIFICIAL "BGK BACKGROUND SPECIES"?
C
        ITYP1(IPLS)=0
        ISPZ1(IPLS)=0
C
        IF (NPBGKP(IPLS,1).EQ.0) GOTO 1000
c slmod begin
        WRITE (97,*) 'IPLS= ',IPLS
c slmod end
C
C  YES. FIND INCIDENT TEST PARTICLE COLLISION PARTNER: ITYP, ISPZ, IREL

        IBGK1=NPBGKP(IPLS,1)
        IUP1=(IBGK1-1)*3+1
        IUP2=(IBGK1-1)*3+2
        IUP3=(IBGK1-1)*3+3
C
C  TRY ATOMS
        DO IATM=1,NATMI
          IF (NPBGKA(IATM).EQ.IBGK1) THEN
            ITYP1(IPLS)=1
            ISPZ1(IPLS)=IATM
            FACT1=CVRSSA(IATM)
            RMAS1=RMASSA(IATM)
            DO IRAD=1,NRAD
              PDEN(IRAD)=PDENA(IATM,IRAD)
              EDEN(IRAD)=EDENA(IATM,IRAD)
            ENDDO
C  FIND INDEX  NRC
            DO NRC=1,NRCA(IATM)
              IP=IDEZ(IBULKA(IATM,NRC),3,3)
              IF (IP.EQ.IPLS) THEN
                INRC1(IPLS)=NRC
                GOTO 10
              ENDIF
            ENDDO
            GOTO 995
C  FIND INDEX IREL
10          DO IAEL=1,NAELI(IATM)
              IF  (LGAEL(IATM,IAEL,1).EQ.IPLS) THEN
                IREL1(IPLS)=LGAEL(IATM,IAEL,0)
                GOTO 100
              ENDIF
            ENDDO
            GOTO 995
          ENDIF
        ENDDO
C  TRY MOLECULES
        DO IMOL=1,NMOLI
          IF (NPBGKM(IMOL).EQ.IBGK1) THEN
            ITYP1(IPLS)=2
            ISPZ1(IPLS)=IMOL
            FACT1=CVRSSM(IMOL)
            RMAS1=RMASSM(IMOL)
            DO IRAD=1,NRAD
              PDEN(IRAD)=PDENM(IMOL,IRAD)
              EDEN(IRAD)=EDENM(IMOL,IRAD)
            ENDDO
C  FIND INDEX  NRC
            DO NRC=1,NRCM(IMOL)
              IP=IDEZ(IBULKM(IMOL,NRC),3,3)
              IF (IP.EQ.IPLS) THEN
                INRC1(IPLS)=NRC
                GOTO 20
              ENDIF
            ENDDO
            GOTO 995
C  FIND INDEX IREL
20          DO IMEL=1,NMELI(IMOL)
              IF  (LGMEL(IMOL,IMEL,1).EQ.IPLS) THEN
                IREL1(IPLS)=LGMEL(IMOL,IMEL,0)
                GOTO 100
              ENDIF
            ENDDO
            GOTO 995
          ENDIF
        ENDDO
C  TRY TEST IONS
        DO IION=1,NIONI
          IF (NPBGKI(IION).EQ.IBGK1) THEN
            ITYP1(IPLS)=3
            ISPZ1(IPLS)=IION
            FACT1=CVRSSI(IION)
            RMAS1=RMASSI(IION)
            DO IRAD=1,NRAD
              PDEN(IRAD)=PDENI(IION,IRAD)
              EDEN(IRAD)=EDENI(IION,IRAD)
            ENDDO
C  FIND INDEX NRC
            DO NRC=1,NRCI(IION)
              IP=IDEZ(IBULKI(IION,NRC),3,3)
              IF (IP.EQ.IPLS) THEN
                INRC1(IPLS)=NRC
                GOTO 30
              ENDIF
            ENDDO
C  FIND INDEX IREL
30          DO IIEL=1,NIELI(IION)
              IF  (LGIEL(IION,IIEL,1).EQ.IPLS) THEN
                IREL1(IPLS)=LGIEL(IION,IIEL,0)
                GOTO 100
              ENDIF
            ENDDO
            GOTO 995
          ENDIF
        ENDDO
        GOTO 995
C
C  SELF COLLISION OR CROSS COLLISION
C
100     CONTINUE
C
        IF (NPBGKP(IPLS,2).EQ.0) THEN
          ITYP2(IPLS)=0
          ISPZ2(IPLS)=0
C
        ELSEIF (NPBGKP(IPLS,2).NE.0) THEN
C
C  CROSS COLLISION, FIND SECOND COLLISION PARTNER
C  THIS IS NOT THE INGOING COLLIDING TESTPARTICLE,
C  (DETERMINING, E.G., MASS AND DENSITY OF BACKGROUND PARTICLE)
C  BUT THE SECOND, 'ARTIFICIAL', PARTNER
C
          ITYP2(IPLS)=IDEZ(NPBGKP(IPLS,2),1,3)
          ISPZ2(IPLS)=IDEZ(NPBGKP(IPLS,2),3,3)
C
          IF (ITYP2(IPLS).EQ.1) THEN
            IATM2=ISPZ2(IPLS)
            FACT2=CVRSSA(IATM2)
            RMAS2=RMASSA(IATM2)
            DO IRAD=1,NRAD
              PDEN2(IRAD)=PDENA(IATM2,IRAD)
              EDEN2(IRAD)=EDENA(IATM2,IRAD)
            ENDDO
            IBGK2=NPBGKA(IATM2)
          ELSEIF (ITYP2(IPLS).EQ.2) THEN
            IMOL2=ISPZ2(IPLS)
            FACT2=CVRSSM(IMOL2)
            RMAS2=RMASSM(IMOL2)
            DO IRAD=1,NRAD
              PDEN2(IRAD)=PDENM(IMOL2,IRAD)
              EDEN2(IRAD)=EDENM(IMOL2,IRAD)
            ENDDO
            IBGK2=NPBGKM(IMOL2)
          ELSEIF (ITYP2(IPLS).EQ.3) THEN
            IION2=ISPZ2(IPLS)
            FACT2=CVRSSI(IION2)
            RMAS2=RMASSI(IION2)
            DO IRAD=1,NRAD
              PDEN2(IRAD)=PDENI(IION2,IRAD)
              EDEN2(IRAD)=EDENI(IION2,IRAD)
            ENDDO
            IBGK2=NPBGKI(IION2)
          ENDIF
C
        ENDIF
C
        IF (RMASSP(IPLS).NE.RMAS1) THEN
          WRITE (6,*) 'INCONSISTENT MASS FOR IPLS= ',IPLS
          WRITE (6,*) 'RMASSP(IPLS),RMAS1 ',RMASSP(IPLS),RMAS1
          CALL EXIT
        ENDIF
1000  CONTINUE

C.................................................................
C
C  SET INDPRO=7, AND
C  WRITE PLASMA DATA ONTO RWK FOR CALL TO SUBR. PLASMA BELOW
C  TI,NI AND (VX,VY,VZ) FOR IPLS=1,NPLSI
C  PLAY SAVE: WRITE WHOLE RWK ARRAY.
C

      IF (NLMLT) THEN



        NXM=MAX(1,NR1STM)
        NYM=MAX(1,NP2NDM)
c        NZM=MAX(1,NT3RDM)
        NZM=NBMLT


      NRWK1=(6+5*NPLS+NAIN)*NRAD
      IF (NID1 < NRWK1) THEN
        WRITE (6,*) ' RWK-ARRAY IS TOO SMALL TO HOLD PLASMA-DATA '
        WRITE (6,*) ' CHECK PARAMETER NSMSTRA '
        CALL EXIT
      END IF
      DO 510 I=1,NRWK1
        RWK(I)=0.
510   CONTINUE
      DO 550 IR=1,NXM
        DO 550 IP=1,NYM
          DO 550 IT=1,NZM

            IRAD=IR + (IP-1)*NR1P2 + (IT-1)*NSTRD
c            IRAD=IR + ((IP-1)+(IT-1)*NP2T3)*NR1P2 + NBLCKA

            RWK  ((0+0*NPLS)*NRAD+1   +(IRAD-1)*1   )= TEIN(IRAD)
            DO 520 IPLS=1,NPLSI
              RWK((1+0*NPLS)*NRAD+IPLS+(IRAD-1)*NPLS)= TIIN(IPLS,IRAD)
              RWK((1+1*NPLS)*NRAD+IPLS+(IRAD-1)*NPLS)= DIIN(IPLS,IRAD)
              RWK((1+2*NPLS)*NRAD+IPLS+(IRAD-1)*NPLS)= VXIN(IPLS,IRAD)
              RWK((1+3*NPLS)*NRAD+IPLS+(IRAD-1)*NPLS)= VYIN(IPLS,IRAD)
              RWK((1+4*NPLS)*NRAD+IPLS+(IRAD-1)*NPLS)= VZIN(IPLS,IRAD)
520         CONTINUE
            RWK  ((1+5*NPLS)*NRAD+1   +(IRAD-1)*1   )= BXIN(IRAD)
            RWK  ((2+5*NPLS)*NRAD+1   +(IRAD-1)*1   )= BYIN(IRAD)
            RWK  ((3+5*NPLS)*NRAD+1   +(IRAD-1)*1   )= BZIN(IRAD)
            RWK  ((5+5*NPLS)*NRAD+1   +(IRAD-1)*1   )= VOL(IRAD)
            DO 530 IAIN=1,NAINI
              RWK((6+5*NPLS)*NRAD+IAIN+(IRAD-1)*NAIN)= ADIN(IAIN,IRAD)
530         CONTINUE
550   CONTINUE
C
c  same as do 550 loop , for additional cell region
c
      DO 570 IRAD=NSURF+1,NSURF+NRADD
        RWK  ((0+0*NPLS)*NRAD+1   +(IRAD-1)*1   )= TEIN(IRAD)
        DO 560 IPLS=1,NPLSI
          RWK((1+0*NPLS)*NRAD+IPLS+(IRAD-1)*NPLS)= TIIN(IPLS,IRAD)
          RWK((1+1*NPLS)*NRAD+IPLS+(IRAD-1)*NPLS)= DIIN(IPLS,IRAD)
          RWK((1+2*NPLS)*NRAD+IPLS+(IRAD-1)*NPLS)= VXIN(IPLS,IRAD)
          RWK((1+3*NPLS)*NRAD+IPLS+(IRAD-1)*NPLS)= VYIN(IPLS,IRAD)
          RWK((1+4*NPLS)*NRAD+IPLS+(IRAD-1)*NPLS)= VZIN(IPLS,IRAD)
560     CONTINUE
        RWK  ((1+5*NPLS)*NRAD+1   +(IRAD-1)*1   )= BXIN(IRAD)
        RWK  ((2+5*NPLS)*NRAD+1   +(IRAD-1)*1   )= BYIN(IRAD)
        RWK  ((3+5*NPLS)*NRAD+1   +(IRAD-1)*1   )= BZIN(IRAD)
        RWK  ((5+5*NPLS)*NRAD+1   +(IRAD-1)*1   )= VOL(IRAD)
        DO 565 IAIN=1,NAINI
          RWK((6+5*NPLS)*NRAD+IAIN+(IRAD-1)*NAIN)= ADIN(IAIN,IRAD)
565     CONTINUE
570   CONTINUE



      ENDIF


C .........................................................................
C  NOW: NEW COLLISION RATES MUST BE SET FOR THE NEXT ITERATION
C .........................................................................
C
C  IN CASE OF CROSS COLLISION, SOME MODIFICATIONS ON THE
C  BACKGROUND PARAMETERS ARE REQUIRED TEMPORARYLY TO ENFORCE
C  A SPECIFIC RELATION BETWEEN TAU_1,2 AND TAU_2,1
C
C  COMPUTE SOME 'DERIVED' PLASMA DATA PROFILES FROM THE PROFILES
C  E.G.: EDRIFT, NEEDED FOR EPLEL3
C
      CALL PLASMA_DERIV
C
      TRCSAV=TRCAMD
      TRCAMD=.FALSE.
C
      CALL LEER(1)
      DO IPLS=1,NPLSI
        LMARK(IPLS)=.FALSE.
      ENDDO
      DO IPLS1=1,NPLSI
        IF (NPBGKP(IPLS1,2).NE.0) THEN
C  IPLS1 IS A CROSS COLLISION TALLY
C  FIND CORRESPONDING 2ND CROSS COLLISION TALLY
          IPLS2=0
          DO IPLS=1,NPLSI
            IF (ITYP1(IPLS).EQ.ITYP2(IPLS1).AND.
     .          ISPZ1(IPLS).EQ.ISPZ2(IPLS1).AND.
     .          ITYP2(IPLS).EQ.ITYP1(IPLS1).AND.
     .          ISPZ2(IPLS).EQ.ISPZ1(IPLS1)) IPLS2=IPLS
          ENDDO
          IF (IPLS2.EQ.0) GOTO 300
          CALL LEER(1)
          WRITE (6,*) 'CORRESPONDING CROSS COLLISION SPECIES '
          WRITE (6,*) 'IPLS1,IPLS2 ',IPLS1,IPLS2
          IF (LMARK(IPLS1).OR.LMARK(IPLS2)) GOTO 300
C  IPLS2 IS THE SECOND CROSS COLLISION TALLY
          WRITE (6,*) 'MODIFY PARAMETERS FOR CROSS COLLISIONALITIES '
          WRITE (6,*) 'IPLS1,IPLS2 ',IPLS1,IPLS2
          CALL LEER(1)
          LMARK(IPLS1)=.TRUE.
          LMARK(IPLS2)=.TRUE.
          DO IRAD=1,NSBOX
            DS1=DIIN(IPLS1,IRAD)
            DIIN(IPLS1,IRAD)=DIIN(IPLS2,IRAD)
            DIIN(IPLS2,IRAD)=DS1
C
            ENERGY(IPLS1,IRAD)=1.5*TIIN(IPLS1,IRAD)+EDRIFT(IPLS1,IRAD)
            ENERGY(IPLS2,IRAD)=1.5*TIIN(IPLS2,IRAD)+EDRIFT(IPLS2,IRAD)
C
            TS1=0.5*(TIIN(IPLS1,IRAD)+TIIN(IPLS2,IRAD))
            TIIN(IPLS1,IRAD)=TS1
            TIIN(IPLS2,IRAD)=TS1
          ENDDO
300       CONTINUE
        ENDIF
      ENDDO
      CALL LEER(2)
C
C
C  COMPUTE SOME 'DERIVED' PLASMA DATA PROFILES FROM THE MODIFIED PROFILES
C
      CALL PLASMA_DERIV
C
C  RESET BGK-ATOMIC AND MOLECULAR DATA ARRAYS
C
      DO IPLS=1,NPLSI
        IF (NPBGKP(IPLS,1).NE.0) THEN
          ITYP=ITYP1(IPLS)
          ISPZ=ISPZ1(IPLS)
          IREL=IREL1(IPLS)
          NRC=INRC1(IPLS)
          IF (ITYP.EQ.1) THEN
            ISP=ISPZ
            KK=IREACA(ISPZ,NRC)
            EBULK=EBULKA(ISPZ,NRC)
            ISCDE=ISCDEA(ISPZ,NRC)
            IESTM=IESTMA(ISPZ,NRC)
            FACTKK=FREACA(ISPZ,NRC)
          ELSEIF (ITYP.EQ.2) THEN
            ISP=NATMI+ISPZ
            KK=IREACM(ISPZ,NRC)
            EBULK=EBULKM(ISPZ,NRC)
            ISCDE=ISCDEM(ISPZ,NRC)
            IESTM=IESTMM(ISPZ,NRC)
            FACTKK=FREACM(ISPZ,NRC)
          ELSEIF (ITYP.EQ.3) THEN
            ISP=NATMI+NMOLI+ISPZ
            KK=IREACI(ISPZ,NRC)
            EBULK=EBULKI(ISPZ,NRC)
            ISCDE=ISCDEI(ISPZ,NRC)
            IESTM=IESTMI(ISPZ,NRC)
            FACTKK=FREACI(ISPZ,NRC)
          ENDIF
          IF (FACTKK.EQ.0.D0) FACTKK=1.D0
C  BGK COLLISION, RESET TABEL3, EPLEL3
          CALL XSTEL(IREL,ISP,IPLS,EBULK,ISCDE,IESTM,KK,FACTKK)
          IF (NPBGKP(IPLS,2).NE.0) THEN
C  CROSS COLLISION, RESET EPLEL3 FOR TRACKLENGTH ESTIMATOR
            IF (NSTORDR >= NRAD) THEN
              DO J=1,NSBOX
                EPLEL3(IREL,J,1)=ENERGY(IPLS,J)
              ENDDO
            ELSE
              NELREL(IREL)=-3
            END IF
          ENDIF
        ENDIF
      ENDDO
C
      TRCAMD=TRCSAV

      DEALLOCATE (PDEN)
      DEALLOCATE (EDEN)
      DEALLOCATE (PDEN2)
      DEALLOCATE (EDEN2)
      DEALLOCATE (ENERGY)

      CALL PLASMA
      CALL PLASMA_DERIV

      RETURN
C
995   CONTINUE
      WRITE (6,*) 'SPECIES ERROR IN MODBGK'
      CALL EXIT
C
      END


      SUBROUTINE JUMUSR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'PARMMOD'
      INCLUDE 'CADGEO'
      INCLUDE 'CLGIN'
      INCLUDE 'COMUSR'
      INCLUDE 'CGRID'


      REAL*8 V1(3),V2(3),V3(3),V4(3)

      INTEGER ICOUNT,NMATCH,ASCCODE1,NLIMI1,NSBOX1,NSBOX2,FP,
     .        JSCRTCH(NLIM),I,NC

c     ASSIGN IGJUM3

      PARAMETER (TOL=1.0D-03)


c...  Don't proceed if there isn't a vacuum grid:
      IF (SBGKI.GT.EBGKI) RETURN

c...  Attempt to read stored IGJUM3:
c      FP=98
c      OPEN(UNIT=FP,FILE='jummap.dat',FORM='UNFORMATTED',
c     .     STATUS='OLD',ERR=10)
c      READ(FP) ASCCODE1,NLIMI1,NSBOX1,NSBOX2
c      IF (ASCCODE1.EQ.ASCCODE.AND.NLIMI1.EQ.(NLIMI/NOPTM1+1).AND.
c     .    NSBOX1  .EQ.MIN(NOPTIM,NSBOX).AND.NSBOX2.EQ.NSBOX) THEN
c        WRITE(0,*) 'LOADING IGJUM3 FROM DATA FILE'
c        READ(FP) ((IGJUM3(J,I),I=1,NLIMI/NOPTM1+1),
c     .                         J=1,MIN(NOPTIM,NSBOX))
c        READ(FP) ((IGJUM4(J,I),I=0,8),J=1,NSBOX)
c        CLOSE(FP)
c        WRITE(0,*) 'DONE'
c        RETURN
c      ENDIF
c10    CONTINUE
c      CLOSE(FP)
c      WRITE(0,*) 'GENERATING IGJUM3'


      DO J=1,NLIM
        JSCRTCH(J) = 0
      ENDDO
 

c ASSUME ADDITIONAL SURFACES DO NOT OVERLAP THE STANDARD GRID
c SO DO NOT CHECK ADDITIONAL SURFACES WHEN A PARTICLE IS ON THE
c STANDARD GRID.  ONLY TURN OFF POLYGON SURFACES THAT ARE 
c INFINITE IN THE Z DIRECTION, OR THAT ARE PART OF THE VACUUM GRID.

c      WRITE(0,*) 'MARK: JUM PARAMS= ',NSURF,NSURF+ASCNCELL

c      WRITE(0,*) 'WARNING: IGJUM4 SET TO CHECK ALL ADDITIONAL '//
c     .           'SURFACES FROM STANDARD GRID'

      DO NCELL=1,NSURF
        IGJUM4(NCELL,0) =  0
        IF (HSTDI.NE.HADDI) THEN
          IGJUM4(NCELL,0)               =  IGJUM4(NCELL,0) + 1
          IGJUM4(NCELL,IGJUM4(NCELL,0)) = -3
        ENDIF
        IGJUM4(NCELL,0)               =  IGJUM4(NCELL,0) + 1
        IGJUM4(NCELL,IGJUM4(NCELL,0)) = -4

c...    Check the standard additional surfaces as well:
c        IGJUM4(NCELL,0)               =  IGJUM4(NCELL,0) + 1
c        IGJUM4(NCELL,IGJUM4(NCELL,0)) = -1

        DO J=1,NLIMI
          IF (P1(3,J).EQ.-1.0D+20.AND.P2(3,J).EQ.1.0D+20) THEN 
            IF (NLIMPB >= NLIMPS) THEN
              IGJUM3(NCELL,J)=1
            ELSE
              CALL BITSET (IGJUM3,0,NOPTIM,NCELL,J,1,NBITS)
            ENDIF
          ENDIF
        ENDDO
c VACUUM GRID CELLS
        DO J=SBGKI,EBGKI
          IF (NLIMPB >= NLIMPS) THEN
            IGJUM3(NCELL,J)=1
          ELSE
            CALL BITSET (IGJUM3,0,NOPTIM,NCELL,J,1,NBITS)
          ENDIF
        ENDDO
      ENDDO

c...  LGJUM1:
      IF (NLIMPB >= NLIMPS) THEN

      ELSE
        WRITE(0,*) 'IGJUM1 MEMORY SAVING OPTION NOT COMPLETE'
      ENDIF







c LOOP OVER ADDITIONAL SURFACES, AND DECIDE IF THEY ARE 
c COINCIDENT WITH THE POLYGON BOUNDARY OF AN ADDITIONAL CELL

c      DO I1=1,1
      I1 = 0
      NC = 1
      DO I=1,ASCNCELL*ASCNCUT

        I1 = I1 + 1
        IF (I1.EQ.ASCNCELL+1) THEN
          I1 = 1
          NC = NC + 1
        ENDIF

        NCELL=ASCCELL(I1)+(NC-1)*ASCNCELL+NSURF

        IGJUM4(NCELL,0) = 0

        ICOUNT=0

        IF (output)  WRITE(0,*) 'MARK PROGRESS... ',I1,ASCNCELL

c LOOP OVER ALL ADDITIONAL CELLS (VACUUM SURFACES)

        DO J=SBGKI,CBGKI

C
          IF (JSCRTCH(J).EQ.2) CYCLE

c IF AN ADDITIONAL SURFACE IS NOT A SIMPLE LINE SEGMENT WITH
c INFINITE EXTENT, THEN ALWAYS CHECK IT
c          IF (P1(3,J).NE.-1.0D+20.AND.P2(3,J).NE.1.0D+20) CYCLE

c TURN OFF ADDITIONAL SURFACE BY DEFAULT

          IF (NLIMPB >= NLIMPS) THEN
            IGJUM3(NCELL,J)=1
          ELSE
            CALL BITSET (IGJUM3,0,NOPTIM,NCELL,J,1,NBITS)
          ENDIF

c LOOP OVER ALL ADDITIONAL SURFACES AND COMPARE WITH
c ADDITIONAL CELL POLYGON SIDES

          IF (ASC3DMODE.EQ.0.OR.ASC3DMODE.EQ.2) THEN

            DO I3=1,ASCNVERTEX(I1)

              IF (I3.EQ.ASCNVERTEX(I1)) THEN
                I4=1
              ELSE
                I4=I3+1
              ENDIF
              IF ((DABS(P1(1,J)-ASCVERTEX(I3*2-1,I1)).LT.TOL.AND.
     .             DABS(P1(2,J)-ASCVERTEX(I3*2  ,I1)).LT.TOL.AND.
     .             DABS(P2(1,J)-ASCVERTEX(I4*2-1,I1)).LT.TOL.AND.
     .             DABS(P2(2,J)-ASCVERTEX(I4*2  ,I1)).LT.TOL).OR.
     .            (DABS(P2(1,J)-ASCVERTEX(I3*2-1,I1)).LT.TOL.AND.
     .             DABS(P2(2,J)-ASCVERTEX(I3*2  ,I1)).LT.TOL.AND.
     .             DABS(P1(1,J)-ASCVERTEX(I4*2-1,I1)).LT.TOL.AND.
     .             DABS(P1(2,J)-ASCVERTEX(I4*2  ,I1)).LT.TOL)) THEN

                IF (NLIMPB >= NLIMPS) THEN
                  IGJUM3(NCELL,J)=0
                ELSE
                  CALL BITSET (IGJUM3,0,NOPTIM,NCELL,J,0,NBITS)
                ENDIF
                ICOUNT=ICOUNT+1

c...            Add surface to short list of surfaces to check when in 
c               this cell:
                IGJUM4(NCELL,0)               =  IGJUM4(NCELL,0) + 1
                IGJUM4(NCELL,IGJUM4(NCELL,0)) =  J

              ENDIF

c            WRITE(6,'(A,3I6,2F14.8,2X,2F14.8,4I6)') '---->',
c     .        I1,j,I3,
c     .        P1(1,J),P1(2,J),
c     .        ASCVERTEX(I3*2-1,I1),ASCVERTEX(I3*2,I1),
c     .        ncell,j,nsurf,igjum3(ncell,j)
c            WRITE(6,'(A,18X,2F14.8,2X,2F14.8)') '---->',
c     .        P2(1,J),P2(2,J),
c     .        ASCVERTEX(I4*2-1,I1),ASCVERTEX(I4*2,I1)
         
            ENDDO

          ELSEIF (ASC3DMODE.EQ.1) THEN

            STOP 'NO LONGER IN USE'

c SEARCH EACH SIDE OF THE CUBE.  ASSUMPTIONS ARE MADE ABOUT THE ORDER
c OF THE VERTICIES BASED ON HOW THE 3D VACUUM MESH IS WRITTEN
c BY DIVIMP.
            DO I3=1,ASCNVERTEX(I1)
              IF (I3.EQ.ASCNVERTEX(I1)) THEN
                I4=1
              ELSE
                I4=I3+1
              ENDIF

              V1(1) = ASCVERTEX(2*I4-1,I1)
              V1(2) = ASCVERTEX(2*I4  ,I1)
              V1(3) = ASCZMIN3D(       I1)
              V2(1) = ASCVERTEX(2*I4-1,I1)
              V2(2) = ASCVERTEX(2*I4  ,I1)
              V2(3) = ASCZMAX3D(       I1)
              V4(1) = ASCVERTEX(2*I3-1,I1)
              V4(2) = ASCVERTEX(2*I3  ,I1)
              V4(3) = ASCZMAX3D(       I1)
              V3(1) = ASCVERTEX(2*I3-1,I1)
              V3(2) = ASCVERTEX(2*I3  ,I1)
              V3(3) = ASCZMIN3D(       I1)

              NMATCH = 0
              DO I5 = 1, 3
c                WRITE(0,'(A,1P,4E12.5)') 
c     .            'V:',V1(I5),V2(I5),V3(I5),V4(I5)
                IF ((DABS(P1(I5,J)-V1(I5)).LT.TOL.AND.
     .               DABS(P2(I5,J)-V2(I5)).LT.TOL.AND.
     .               DABS(P3(I5,J)-V3(I5)).LT.TOL.AND.
     .               DABS(P4(I5,J)-V4(I5)).LT.TOL).OR.
     .              (DABS(P1(I5,J)-V4(I5)).LT.TOL.AND.
     .               DABS(P2(I5,J)-V3(I5)).LT.TOL.AND.
     .               DABS(P3(I5,J)-V2(I5)).LT.TOL.AND.
     .               DABS(P4(I5,J)-V1(I5)).LT.TOL))  NMATCH = NMATCH + 1
              ENDDO

c MATCH FOUND FOR CELL SIDE
              IF (NMATCH.EQ.3) THEN
c                WRITE(0,*) 'MATCH 0:',i1,j

c...            Add surface to short list of surfaces to check when in 
c               this cell:
                IGJUM4(NCELL,0)               =  IGJUM4(NCELL,0) + 1
                IGJUM4(NCELL,IGJUM4(NCELL,0)) =  J

                JSCRTCH(J) = JSCRTCH(J) + 1

                IF (NLIMPB >= NLIMPS) THEN
                  IGJUM3(NCELL,J)=0
                ELSE
                  CALL BITSET (IGJUM3,0,NOPTIM,NCELL,J,0,NBITS)
                ENDIF
c                IGJUM3(NCELL,J)=0

                ICOUNT=ICOUNT+1
              ENDIF

            ENDDO

c CELL END SURFACE 5:
            XMIN =  1.0D+20
            XMAX = -1.0D+20
            YMIN =  1.0D+20
            YMAX = -1.0D+20
            DO I3=1,ASCNVERTEX(I1)
              XMIN=MIN(XMIN,ASCVERTEX(2*I3-1,I1))
              XMAX=MAX(XMAX,ASCVERTEX(2*I3-1,I1))
              YMIN=MIN(YMIN,ASCVERTEX(2*I3  ,I1))
              YMAX=MAX(YMAX,ASCVERTEX(2*I3  ,I1))
            ENDDO

            V1(1) = XMAX
            V1(2) = YMIN
            V1(3) = ASCZMIN3D(I1)
            V2(1) = XMIN      
            V2(2) = YMIN      
            V2(3) = ASCZMIN3D(I1)
            V3(1) = XMAX      
            V3(2) = YMAX      
            V3(3) = ASCZMIN3D(I1)
            V4(1) = XMIN      
            V4(2) = YMAX      
            V4(3) = ASCZMIN3D(I1)
	    		      
            NMATCH = 0
            DO I5 = 1, 3
c              WRITE(6,'(A,1P,4E12.5,I6)') 
c     .          'P:',V1(I5,J),V2(I5,J),V3(I5,J),V4(I5,J),J
c              WRITE(6,'(A,1P,4E12.5)') 
c     .          'V:',V1(I5),V2(I5),V3(I5),V4(I5)
              IF (DABS(P1(I5,J)-V1(I5)).LT.TOL.AND.
     .            DABS(P2(I5,J)-V2(I5)).LT.TOL.AND.
     .            DABS(P3(I5,J)-V3(I5)).LT.TOL.AND.
     .            DABS(P4(I5,J)-V4(I5)).LT.TOL) NMATCH = NMATCH + 1
            ENDDO
            IF (NMATCH.EQ.3) THEN

c...          Add surface to short list of surfaces to check when in 
c             this cell:
              IGJUM4(NCELL,0)               =  IGJUM4(NCELL,0) + 1
              IGJUM4(NCELL,IGJUM4(NCELL,0)) =  J

              JSCRTCH(J) = JSCRTCH(J) + 1

c               WRITE(0,*) 'MATCH 5:',i1,j
              IF (NLIMPB >= NLIMPS) THEN
                IGJUM3(NCELL,J)=0
              ELSE
                CALL BITSET (IGJUM3,0,NOPTIM,NCELL,J,0,NBITS)
              ENDIF
c              IGJUM3(NCELL,J)=0
c              ICOUNT=ICOUNT+1
            ENDIF

            V1(3) = ASCZMAX3D(I1)
            V2(3) = ASCZMAX3D(I1)
            V4(3) = ASCZMAX3D(I1)
            V3(3) = ASCZMAX3D(I1)
	    		      
            NMATCH = 0
            DO I5 = 1, 3
c              WRITE(0,'(A,1P,4E12.5)') 
c     .          'V:',V1(I5),V2(I5),V3(I5),V4(I5)
              IF (DABS(P1(I5,J)-V1(I5)).LT.TOL.AND.
     .            DABS(P2(I5,J)-V2(I5)).LT.TOL.AND.
     .            DABS(P3(I5,J)-V3(I5)).LT.TOL.AND.
     .            DABS(P4(I5,J)-V4(I5)).LT.TOL) NMATCH = NMATCH + 1
            ENDDO
            IF (NMATCH.EQ.3) THEN

c...          Add surface to short list of surfaces to check when in 
c             this cell:
              IGJUM4(NCELL,0)               =  IGJUM4(NCELL,0) + 1
              IGJUM4(NCELL,IGJUM4(NCELL,0)) =  J

              JSCRTCH(J) = JSCRTCH(J) + 1

c               WRITE(0,*) 'MATCH 6:',i1,j
              IF (NLIMPB >= NLIMPS) THEN
                IGJUM3(NCELL,J)=0
              ELSE
                CALL BITSET (IGJUM3,0,NOPTIM,NCELL,J,0,NBITS)
              ENDIF
c              IGJUM3(NCELL,J)=0
c              ICOUNT=ICOUNT+1
            ENDIF

          ENDIF

        ENDDO

c        WRITE(6,*) 'ICOUNT=',ncell,icount,ascnvertex(i1)
      
        IF (ICOUNT.EQ.ASCNVERTEX(I1)) THEN
c        IF (ASC3DMODE.EQ.0.AND.ICOUNT.EQ.ASCNVERTEX(I1).OR.
c     .      ASC3DMODE.EQ.1.AND.ICOUNT.EQ.ASCNVERTEX(I1)) THEN
c     .      ASC3DMODE.EQ.1.AND.ICOUNT.EQ.ASCNVERTEX(I1)+2)) THEN
          DO J=1,SBGKI-1
            IF (P1(3,J).EQ.-1.0D+20.AND.P2(3,J).EQ.1.0D+20) THEN
              IF (NLIMPB >= NLIMPS) THEN
                IGJUM3(NCELL,J)=1
              ELSE
                CALL BITSET (IGJUM3,0,NOPTIM,NCELL,J,1,NBITS)
              ENDIF
c             IGJUM3(NCELL,J)=1
            ENDIF
          ENDDO
        ELSE
c        ELSEIF (ASC3DMODE.EQ.1.OR.ASC3DMODE.EQ.2) THEN
c...      Make sure that all the wall surfaces are checked:
          IGJUM4(NCELL,0)               =  IGJUM4(NCELL,0) + 1
          IGJUM4(NCELL,IGJUM4(NCELL,0)) =  -1
        ENDIF

c...    Add appropriate toroidal region boundary surfaces:
        IF (ASC3DMODE.EQ.2) THEN
          IF (NC.GT.1      ) THEN
            IGJUM4(NCELL,0)               =  IGJUM4(NCELL,0) + 1
            IGJUM4(NCELL,IGJUM4(NCELL,0)) =  CBGKI+NC-1
          ENDIF
          IF (NC.LT.ASCNCUT) THEN
            IGJUM4(NCELL,0)               =  IGJUM4(NCELL,0) + 1
            IGJUM4(NCELL,IGJUM4(NCELL,0)) =  CBGKI+NC
          ENDIF
        ENDIF                         

c...    Always check the manditory surfaces:      
c        IF (ASC3DMODE.EQ.1.OR.ASC3DMODE.EQ.2) THEN
          IF (HADDI.NE.EBGKI) THEN
            IGJUM4(NCELL,0)               =  IGJUM4(NCELL,0) + 1
            IGJUM4(NCELL,IGJUM4(NCELL,0)) =  -2
          ENDIF
          IGJUM4(NCELL,0)               =  IGJUM4(NCELL,0) + 1
          IGJUM4(NCELL,IGJUM4(NCELL,0)) =  -4
c        ENDIF
      ENDDO

c...  For additional cell 1:
      IGJUM4(NSURF+1,0) =  3
      IGJUM4(NSURF+1,1) = -1
      IGJUM4(NSURF+1,2) = -2
      IGJUM4(NSURF+1,3) = -4


c      DO NCELL=1,NSBOX
c        WRITE(6,'(A,3I6,3X,10I6)') 'MARK: IGJUM4=',
c     .    NCELL,NSURF,MAX(0,NCELL-NSURF),
c     .    (IGJUM4(NCELL,I1),I1=0,9)
c      ENDDO


c...  Dump IJGUM3 to a file:
c      WRITE(0,*) 'WRITING IGJUM3 TO FILE',ASCCODE
c      FP=98
c      OPEN(UNIT=FP,FILE='jummap.dat',FORM='UNFORMATTED',
c     .     STATUS='REPLACE',ERR=98)
c      WRITE(FP) ASCCODE,NLIMI/NOPTM1+1,MIN(NOPTIM,NSBOX),NSBOX
c      WRITE(FP) ((IGJUM3(J,I),I=1,NLIMI/NOPTM1+1),
c     .                        J=1,MIN(NOPTIM,NSBOX))
c      WRITE(FP) ((IGJUM4(J,I),I=0,8),J=1,NSBOX)
c      CLOSE(FP)
c      WRITE(0,*) 'DONE'

      RETURN
 98   WRITE(0,*) 'UNABLE TO OPEN IGJUM3 STORAGE FILE'
 99   STOP
      END


      SUBROUTINE USRTIM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)     
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'CCONA'

      INTEGER i1,i2,fp
      LOGICAL dataloaded
      CHARACTER*1024 cdum1

      DATA dataloaded /.FALSE./
      SAVE

c...  Set up parameters:
      ENTRY USRTI0

      IF (dataloaded) RETURN

      fp = 31

      READ(fp,'(A)',ERR=97,END=98) cdum1

      IF (cdum1(1:10).NE.'[TIME-TO-I') THEN
        WRITE(0,*) 'ERROR: TIME-TO-IONISATION PARAMETERS EXPECTED '//
     .             'BUT NOT FOUND'
        WRITE(0,*) 'CDUM1:',cdum1
        STOP
      ENDIF

      READ(fp,'(A)',ERR=97,END=98) cdum1
      DO WHILE(cdum1(1:1).EQ.'*')
        READ(fp,'(A)',ERR=97,END=98) cdum1     
      ENDDO  
      BACKSPACE fp

c...  Load the number of regions for recording time-to-ionisation
c     information:
      READ(fp,*) timnum

      DO i1 = 1, timnum

c...    Read bounding rectangle(s):
        READ(fp,*) idum1,timx1(i1),timx2(i1),timy1(i1),timy2(i1),
     .             timnbin(i1)
        IF (timnbin(i1).NE.0) THEN
c...      Load bin information:
          timnbin(i1) = timnbin(i1) + 2
          DO i2 = 1, timnbin(i1)+1
            READ(fp,*) timbin(i1,i2)             
          ENDDO      
        ENDIF
c...    Scale units from [m] to [cm]:
        timx1(i1) = timx1(i1) * 100.0D0
        timx2(i1) = timx2(i1) * 100.0D0
        timy1(i1) = timy1(i1) * 100.0D0
        timy2(i1) = timy2(i1) * 100.0D0
      ENDDO

      dataloaded = .TRUE.

      RETURN

c...  Process atom ionisation events: 
      ENTRY USRTI1

c...  Check time recording regions to see if particle
c     qualifies:
      DO i1 = 1, timnum
        IF (x0.GT.timx1(i1).AND.x0.LT.timx2(i1).AND.
     .      y0.GT.timy1(i1).AND.y0.LT.timy2(i1)) THEN

          WRITE(6,*) '--TIME-->',npanu,i1


c...      Update the particle count and log time:
          timcnt(i1) = timcnt(i1) + weight
          timavg(i1) = timavg(i1) + weight * time

c...      Bin the particle time:
          DO i2 = 1, timnbin(i1)
            IF (time.GE.timbin(i1,i2).AND.time.LT.timbin(i1,i2+1)) 
     .        timbinc(i1,i2) = timbinc(i1,i2) + weight               
          ENDDO

        ENDIF

      ENDDO

      RETURN


c...  Dump results to the transfer file:
      ENTRY USRTI2

      fp = 32

      WRITE(fp,'(A)') '[TIME-TO-IONISATION STATISTICS]'

      WRITE(fp,'(A8,4A8,2X,2A14)') 'REGION','X1','X2','Y1','Y2','COUNT',
     .                            'AVERAGE'
      DO i1 = 1, timnum
        timavg(i1) = timavg(i1) / (timcnt(i1) + EPS10)

        WRITE(fp,'(I8,4F8.3,2X,F14.4,1P,E14.4,0P)') 
     .    i1,timx1(i1)*0.01D0,timx2(i1)*0.01D0,
     .       timy1(i1)*0.01D0,timy2(i1)*0.01D0,
     .       timcnt(i1),timavg(i1)

      ENDDO

      DO i1 = 1, timnum
        IF (timnbin(i1).EQ.0) CYCLE

        totbinc = 0.0D0
        DO i2 = 1, timnbin(i1)
          totbinc = totbinc + timbinc(i1,i2)
        ENDDO

        WRITE(fp,'(2A8,4A14)') 'REGION','BIN','T1','T2',
     .                         'FRACTION','COUNT'
        DO i2 = 1, timnbin(i1)
          WRITE(fp,'(2I8,1P,2E12.4,0P,2F14.4)')
     .      i1,i2,timbin(i1,i2),timbin(i1,i2+1),
     .      timbinc(i1,i2)/(totbinc+EPS10)*100.0D0,timbinc(i1,i2)
        ENDDO

      ENDDO

      RETURN
97    WRITE(0,*) 'ERROR: ERROR READING TRANSFER FILE'
      GOTO 99
98    WRITE(0,*) 'ERROR: UNEXPECTED EOF ON TRANSFER FILE'
99    STOP
      END




      LOGICAL FUNCTION CHKPNT(a1,a2,b1,b2,c1,c2,d1,d2)
      IMPLICIT none

      REAL*8 a1,a2,b1,b2,c1,c2,d1,d2

      REAL*8 TOL
      PARAMETER (TOL = 1.0D-07)

      REAL*8 t0,e1,e2,f1,f2,fact,tab,tcd

      CHKPNT = .FALSE.

      IF (DABS(c1-d1).LT.TOL.AND.DABS(c2-d2).LT.TOL) THEN
        WRITE(0,*) 'CHKPNT: ZERO LENGTH LINE SEGMENT' 
        WRITE(6,*) 'CHKPNT: ZERO LENGTH LINE SEGMENT' 
        CALL EXIT
      ENDIF

c...  Calculate cross-product to see if the lines are parallel:
      IF (DABS((a1 - b1) * (c2 - d2)).LT.1.0D-8.AND.
     .    DABS((a2 - b2) * (c1 - d1)).LT.1.0D-8) RETURN

c...  Find projection of AB onto CD:
      t0 = ((a1 - c1) * (d1 - c1) + (a2 - c2) * (d2 - c2)) /
     .     ((c1 - d1)**2 + (c2 - d2)**2)
      e1 = c1 + t0 * (d1 - c1)
      e2 = c2 + t0 * (d2 - c2)

c...  Calculate F, the intersection point between AB and CD:
      fact = (e1 - a1) * (b1 - a1) + (e2 - a2) * (b2 - a2)

c...  Determine the parametric location of F on AB:
      IF (DABS(fact).LT.1.0D-10) THEN
        tab = 0.0
      ELSE
        tab = ((a1 - e1)**2 + (a2 - e2)**2) / fact
      ENDIF

c...  Determine the parametric location of F on CD:
      f1 = a1 + tab * (b1 - a1)
      f2 = a2 + tab * (b2 - a2)

      IF (DABS(d1-c1).GT.DABS(d2-c2)) THEN
        tcd = (f1 - c1) / (d1 - c1)
      ELSE
        tcd = (f2 - c2) / (d2 - c2)
      ENDIF

c      WRITE(0,*) 'CHKPNT:',tab,tcd,
c     .  tab.GT.0.0D0+1.0D-8.AND.tab.LT.1.0D0-1.0D-8,
c     .  tcd.GT.0.0D0+1.0D-8.AND.tcd.LT.1.0D0-1.0D-8

c...  Decide if the lines interect between their respective end points:
      IF (tab.GT.0.0D0+1.0D-8.AND.tab.LT.1.0D0-1.0D-8.AND.
     .    tcd.GT.0.0D0+1.0D-8.AND.tcd.LT.1.0D0-1.0D-8)
     .  CHKPNT = .TRUE.

      RETURN
      END


c slmod end
