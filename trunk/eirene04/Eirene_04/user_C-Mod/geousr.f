C
C
C
      SUBROUTINE GEOUSR
C
C   PREPARE DATA FOR LIMITER-SURFACES
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CADGEO
      USE COMPRT
      USE CTRCEI
      USE CCONA
      USE CGEOM
      USE CGRID
      USE CLGIN
      USE CINIT
      USE CPOLYG
      USE CTRIG
      IMPLICIT NONE
      CHARACTER(80) :: ZEILE
      REAL(DP) :: XCOOR, YCOOR, ZCOOR, xan, xen, yan, yen, zan, zen 
      INTEGER :: NADMOD, NASMOD, NORMOD, NRS, IPUNKT, I,NSSIR, NSSIP,
     .           IDIR, IR, IP, IT, IC, IN, NAS
      logical :: lx1,lx2,lx3,lx4,ly1,ly2,ly3,ly4
C
C MODIFY GEOMETRY
C
C
      READ (IUNIN,'(A80)') ZEILE
      READ (IUNIN,'(3I6)') NADMOD,NASMOD,NORMOD

      DO I=1,NADMOD
        READ (IUNIN,'(2I6,3E12.4)') NRS,IPUNKT,XCOOR,YCOOR,ZCOOR

        SELECT CASE(IPUNKT)

        CASE DEFAULT
           WRITE (6,*) 'WRONG POINTNUMBER IN ADDUSR '
           WRITE (6,*) 'INPUT LINE READING'
           WRITE (6,'(2I6,1P,3E12.4)') NRS,IPUNKT,XCOOR,YCOOR,ZCOOR
           WRITE (6,*) ' IS IGNORED '

        CASE (1)
           P1(1,NRS)=XCOOR
           P1(2,NRS)=YCOOR
           P1(3,NRS)=ZCOOR

        CASE (2)
           P2(1,NRS)=XCOOR
           P2(2,NRS)=YCOOR
           P2(3,NRS)=ZCOOR

        CASE (3)
           P3(1,NRS)=XCOOR
           P3(2,NRS)=YCOOR
           P3(3,NRS)=ZCOOR

        CASE (4)
           P4(1,NRS)=XCOOR
           P4(2,NRS)=YCOOR
           P4(3,NRS)=ZCOOR

        CASE (5)
           P5(1,NRS)=XCOOR
           P5(2,NRS)=YCOOR
           P5(3,NRS)=ZCOOR

        CASE (6)
           P6(1,NRS)=XCOOR
           P6(2,NRS)=YCOOR
           P6(3,NRS)=ZCOOR

        END SELECT
      ENDDO

      DO I=1,NASMOD
        READ (IUNIN,'(5I6)') NAS,IPUNKT,NSSIR,NSSIP
        IF (IPUNKT.EQ.1) THEN
          P1(1,NAS)=XPOL(NSSIR,NSSIP)
          P1(2,NAS)=YPOL(NSSIR,NSSIP)
        ELSEIF (IPUNKT.EQ.2) THEN
          P2(1,NAS)=XPOL(NSSIR,NSSIP)
          P2(2,NAS)=YPOL(NSSIR,NSSIP)
        ELSE
          WRITE (6,*) 'WRONG POINTNUMBER IN ADDUSR '
          WRITE (6,*) 'INPUT LINE READING'
          WRITE (6,'(5I6)') NAS,IPUNKT,NSSIR,NSSIP
          WRITE (6,*) ' IS IGNORED '
        ENDIF
      ENDDO

      DO I = 1, NORMOD
        READ (IUNIN,'(5I6)') IDIR,IR,IP
        IF (IDIR == 1) THEN
          PLNX(IR,IP) = -PLNX(IR,IP)
          PLNY(IR,IP) = -PLNY(IR,IP)
        ELSE IF (IDIR == 2) THEN
          PPLNX(IR,IP) = -PPLNX(IR,IP)
          PPLNY(IR,IP) = -PPLNY(IR,IP)
        ELSE
          WRITE (6,*) ' IDIR =',IDIR,' NOT FORESEEN IN GEOUSR '
          WRITE (6,*) IDIR,IR,IP
          WRITE (6,*) ' IS IGNORED '
        END IF
      END DO
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
C
C  SET SOME VOLUMES EXPLIZIT
C
C
C  MODIFY REFLECTION MODEL AT TARGET PLATES
C
C     do i=1,nlimps
C       do isp=1,natmi+nmoli+nioni
C         recyct(isp,i)=1.
C       enddo
C     enddo

      if (nlimi >= 50) then
      p1(1,1) = xtrian(necke(2,3254))
      p1(2,1) = ytrian(necke(2,3254))

      p2(1,6) = xtrian(necke(3,3395))
      p2(2,6) = ytrian(necke(3,3395))
      

      p2(1,50) = xtrian(necke(3,3083))
      p2(2,50) = ytrian(necke(3,3083))
      
      p1(1,24) = xtrian(necke(2,3251))
      p1(2,24) = ytrian(necke(2,3251))

      p2(1,23) = xtrian(necke(3,3856))
      p2(2,23) = ytrian(necke(3,3856))
      
      p1(1,7) = xtrian(necke(2,3884))
      p1(2,7) = ytrian(necke(2,3884))
      end if

      RETURN
      END
