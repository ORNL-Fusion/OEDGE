C EIRENE06 COMPILATION
C ===== SOURCE: broad_usr.f


      SUBROUTINE BROAD_USR
      IMPLICIT NONE
      RETURN
      END
C ===== SOURCE: geomd.f
C
*//GEOMD//
C=======================================================================
C          S U B R O U T I N E   G E O M D
C=======================================================================
      SUBROUTINE GEOMD(NDXA,NDYA,NPLP,NR1ST,
     .                 PUX,PUY,PVX,PVY)
C
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CGEOM
      IMPLICIT NONE
C
      REAL(DP), INTENT(OUT) :: PUX(*),PUY(*),PVX(*),PVY(*)
      INTEGER, INTENT(INOUT) :: NDXA,NDYA,NPLP,NR1ST

      character(110) :: zeile
      REAL(DP) :: br(0:ndxp,0:ndyp,4),bz(0:ndxp,0:ndyp,4)
C
C  GEOMETRY DATA: CELL VERTICES (LINDA ---> EIRENE)
      REAL(DP) ::
     R  X1(NDX),Y1(NDX),X2(NDX),Y2(NDX),X3(NDX),Y3(NDX),
     R  X4(NDX),Y4(NDX)
      REAL(DP) :: DX, DY
      INTEGER :: I0, I0E, IX, IY, I1, I2, I3, I4, IPART, I, J

      ndxa=0
      ndya=0

      read (30,*)
      read (30,*)
      read (30,*)
      read (30,*)

1     continue
      read (30,'(a110)',end=99) zeile
      i0=index(zeile,'(')
      i0e=index(zeile,')')
      read (zeile(i0+1:i0e-1),*) ix,iy
      ndxa=max(ndxa,ix)
      ndya=max(ndya,iy)
      i1=index(zeile,': (')
      i2=index(zeile(i1+3:),')')+i1+2
      read (zeile(i1+3:i2-1),*) br(ix,iy,4),bz(ix,iy,4)
      i3=index(zeile(i2+1:),'(')+i2
      i4=i3+index(zeile(i3+1:),')')
      read (zeile(i3+1:i4-1),*) br(ix,iy,3),bz(ix,iy,3)

      read (30,'(a110)') zeile

      read (30,'(a110)') zeile
      i1=index(zeile,'(')
      i2=index(zeile,')')
      read (zeile(i1+1:i2-1),*) br(ix,iy,1),bz(ix,iy,1)
      i3=i2+index(zeile(i2+1:),'(')
      i4=i2+index(zeile(i2+1:),')')
      read (zeile(i3+1:i4-1),*) br(ix,iy,2),bz(ix,iy,2)

      read (30,*)
      goto 1


99    continue
      ndxa=ndxa-1
      ndya=ndya-1
C
      DO 1015 IY=1,NDYA
        DO 1014 IX=1,NDXA
          X1(IX)=br(ix,iy,1)
          Y1(IX)=bz(ix,iy,1)
          X2(IX)=br(ix,iy,2)
          Y2(IX)=bz(ix,iy,2)
          X3(IX)=br(ix,iy,4)
          Y3(IX)=bz(ix,iy,4)
          X4(IX)=br(ix,iy,3)
          Y4(IX)=bz(ix,iy,3)
1014    CONTINUE
        CALL MSHPROJ (X1,Y1,X2,Y2,X3,Y3,X4,Y4,PUX,PUY,PVX,PVY,NDXA,
     .                NR1ST,IY)
1015  CONTINUE
C
C SEARCH FOR THE CUTS
C
      IPART=1
      NPOINT(1,IPART)=1
      IY=1
      DO IX=1,NDXA
        DX=BR(IX+1,IY,1)-BR(IX,IY,2)
        DY=BZ(IX+1,IY,1)-BZ(IX,IY,2)
        IF (DX*DX+DY*DY.GT.EPS10) THEN
C CUT GEFUNDEN
          NPOINT(2,IPART)=IX+1+IPART-1
          IPART=IPART+1
          NPOINT(1,IPART)=IX+1+IPART-1
        ENDIF
      ENDDO
      NPOINT(2,IPART)=NDXA+IPART
C
      NPLP=IPART
C
      DO IY=1,NDYA
        DO IPART=1,NPLP
          DO IX=NPOINT(1,IPART),NPOINT(2,IPART)-1
            XPOL(IY,IX)=BR(IX-(IPART-1),IY,1)
            YPOL(IY,IX)=BZ(IX-(IPART-1),IY,1)
          ENDDO
          XPOL(IY,NPOINT(2,IPART))=BR(NPOINT(2,IPART)-IPART,IY,2)
          YPOL(IY,NPOINT(2,IPART))=BZ(NPOINT(2,IPART)-IPART,IY,2)
        ENDDO
      ENDDO
C INTRODUCE OUTERMOST RADIAL POLYGON
      DO IPART=1,NPLP
        DO IX=NPOINT(1,IPART),NPOINT(2,IPART)-1
          XPOL(NDYA+1,IX)=BR(IX-(IPART-1),NDYA,4)
          YPOL(NDYA+1,IX)=BZ(IX-(IPART-1),NDYA,4)
        ENDDO
        XPOL(NDYA+1,NPOINT(2,IPART))=BR(NPOINT(2,IPART)-IPART,NDYA,3)
        YPOL(NDYA+1,NPOINT(2,IPART))=BZ(NPOINT(2,IPART)-IPART,NDYA,3)
      ENDDO
C
      DO J=1,NDYA+1
        DO I=1,NPOINT(2,NPLP)
          XPOL(J,I)=XPOL(J,I)*100.
          YPOL(J,I)=YPOL(J,I)*100.
        END DO
      END DO
C
      ndxa=npoint(2,nplp)-1
c
C     do j=1,ndya+1
C       write (iunout,*)
C       write (iunout,*) 'in geomd polygon ',j
C       write (iunout,'(1p,6e12.4)') 
C    .        (xpol(j,i),ypol(j,i),i=1,npoint(2,nplp))
C     enddo
C
      RETURN
*//END GEOMD//
      END
C ===== SOURCE: geousr.f
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
           WRITE (iunout,*) 'WRONG POINTNUMBER IN ADDUSR '
           WRITE (iunout,*) 'INPUT LINE READING'
           WRITE (iunout,'(2I6,1P,3E12.4)') NRS,IPUNKT,XCOOR,YCOOR,ZCOOR
           WRITE (iunout,*) ' IS IGNORED '

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
          WRITE (iunout,*) 'WRONG POINTNUMBER IN ADDUSR '
          WRITE (iunout,*) 'INPUT LINE READING'
          WRITE (iunout,'(5I6)') NAS,IPUNKT,NSSIR,NSSIP
          WRITE (iunout,*) ' IS IGNORED '
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
          WRITE (iunout,*) ' IDIR =',IDIR,' NOT FORESEEN IN GEOUSR '
          WRITE (iunout,*) IDIR,IR,IP
          WRITE (iunout,*) ' IS IGNORED '
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

      RETURN
      END
C ===== SOURCE: iniusr.f


      SUBROUTINE iniUSR
      IMPLICIT NONE
      RETURN
      END
C ===== SOURCE: leausr.f


      FUNCTION LEAUSR(A,B,C)
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: A, B, C
      INTEGER :: LEAUSR
      LEAUSR=1
      RETURN
      END
C ===== SOURCE: modusr.f
c
c
      subroutine modusr
      return
      end
C ===== SOURCE: mshadj.f
C
C
      SUBROUTINE MSHADJ (X1,Y1,X2,Y2,XPLG,YPLG,NPLP,NPOINT,M1,NDX,IR)
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: X1(*),Y1(*),X2(*),Y2(*)
      INTEGER, INTENT(IN) ::M1, NDX, IR
      REAL(DP), INTENT(INOUT) :: XPLG(M1,*),YPLG(M1,*)
      INTEGER, INTENT(INOUT) :: NPOINT(2,*), NPLP
      REAL(DP) :: EPS, D, D1
      INTEGER :: NPRT, NP, I
      LOGICAL :: LDUMCL,LPRT
C  LPRT=.TRUE.  : VALID PART
C  LDUMCL=.TRUE.: DUMMY CELL
C
      EPS=1.E-20
      NPRT=0
      NP=0
      LDUMCL=.FALSE.
      LPRT=.FALSE.
C
      DO 1 I=1,NDX
        D=ABS((X1(I)-X2(I))**2+(Y1(I)-Y2(I))**2)
        IF (D.LE.EPS) THEN
          IF (LPRT) THEN
C  ENDE DER VALID ZELLEN
            NPOINT(2,NPRT)=NP
            LPRT=.FALSE.
            LDUMCL=.TRUE.
          ENDIF
        ELSEIF (D.GT.EPS) THEN
C  STARTPUNKT DIESES POLYGONS = ERSTER VALID PUNKT?
          IF (.NOT.LPRT.AND..NOT.LDUMCL) THEN
            LPRT=.TRUE.
            NPRT=NPRT+1
            NP=NP+1
            XPLG(IR,NP)=X1(I)
            YPLG(IR,NP)=Y1(I)
            NPOINT(1,NPRT)=NP
          ELSEIF (.NOT.LPRT.AND.LDUMCL) THEN
            D1=ABS((X1(I+1)-X2(I+1))**2+(Y1(I+1)-Y2(I+1))**2)
            IF (D1.GT.EPS) THEN
C  ENDE DER DUMMYZELLEN: SUCHE NAECHSTE ZELLE MIT D1.GT.0.
              LDUMCL=.FALSE.
              LPRT=.TRUE.
              NPRT=NPRT+1
              NPOINT(1,NPRT)=NP
            ENDIF
          ENDIF
        ENDIF
        NP=NP+1
        XPLG(IR,NP)=X2(I)
        YPLG(IR,NP)=Y2(I)
1     CONTINUE
      NPOINT(2,NPRT)=NP
      NPLP=NPRT
      RETURN
      END
C ===== SOURCE: plausr.f


      SUBROUTINE PLAUSR
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CSTEP
      USE CGRID
      USE CINIT
      USE COMSOU
      IMPLICIT NONE
      REAL(DP) :: FACTOR, STEP
      INTEGER :: IG, IR, IT, IN, IP, IPLS, NRWL
      REAL(DP) :: XMACH


      RETURN
      END
C ===== SOURCE: pltusr.f
C
C
      SUBROUTINE PLTUSR(PLABLE,J)
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CADGEO
      USE CCONA
      USE CLGIN
      IMPLICIT NONE
      LOGICAL, INTENT(INOUT) :: PLABLE
      INTEGER, INTENT(IN) :: J
      INTEGER :: IB, MERK, IER
      REAL(DP) :: XS, YS, ZX0, ZY0, ZZ0, CX, CY, CZ, RZYLB, B0B, B1B,
     .          B2B, B3B, F0B, F1B, F2B, F3B, X0, Y0, B0, B1, B2, Z1, Z2
      INTEGER :: IN
      RETURN
      END
C ===== SOURCE: prousr.f
CDK USER
C
C   USER SUPPLIED SUBROUTINES
C
C           ************
C           *  TEXTOR  *  (HELIUM INCLUDED)
C           ************
C
      SUBROUTINE PROUSR (PRO,INDX,P0,P1,P2,P3,P4,P5,PROVAC,N)
C
C     P0 : CENTRAL VALUE
C     P1 : STARTING RADIUS FOR POLYNOMIAL
C     P2 : SWITCH FOR PHASE 1: OH-PHASE
C                           2: NI-PHASE
C     P3 : FACTOR FOR TI: TI=K*TE
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CGRID
      USE CGEOM
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: P0, P1, P2, P3, P4, P5, PROVAC
      REAL(DP), INTENT(OUT) :: PRO(*)
      INTEGER, INTENT(IN) :: INDX, N

      RETURN
      END
C ===== SOURCE: refusr.f


      SUBROUTINE REFUSR
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: XMW,XCW,XMP,XCP,ZCOS,ZSIN,EXPI,RPROB,
     .                        E0TERM
      INTEGER, INTENT(IN) :: IGASF,IGAST
      ENTRY RF0USR
      ENTRY SPTUSR
      ENTRY SP0USR
      ENTRY SP1USR
      ENTRY RF1USR (XMW,XCW,XMP,XCP,IGASF,IGAST,ZCOS,ZSIN,EXPI,
     .              RPROB,E0TERM,*,*,*,*)
      RETURN
      END
C ===== SOURCE: retusr.f
c
c
      subroutine retusr(sig)
      USE PRECISION
      implicit none
      real(dp), intent(in) :: sig
      return
      end
C ===== SOURCE: samusr.f
C
C
      SUBROUTINE SAMUSR (NLSF,X0,Y0,Z0,
     .              SORAD1,SORAD2,SORAD3,SORAD4,SORAD5,SORAD6,
     .              IRUSR,IPUSR,ITUSR,IAUSR,IBUSR,
     .              TIWL,TEWL,DIWL,VXWL,VYWL,VZWL,EFWL,SHWL,WEISPZ)
C
C  SAMPLE INITAL COORDIANTES X,Y,Z ON ADDITIONAL SURFACE NLLI
C
      USE PRECISION
      USE PARMMOD
      USE CADGEO
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: SORAD1,SORAD2,SORAD3,SORAD4,SORAD5,SORAD6
      REAL(DP), INTENT(OUT) :: X0,Y0,Z0,TEWL,TIWL(*),DIWL(*),
     .                         VXWL(*),VYWL(*),VZWL(*),
     .                         EFWL(*), SHWL, WEISPZ(*)
      INTEGER, INTENT(IN) :: NLSF,is1, is2
      INTEGER, INTENT(OUT) :: IRUSR, IPUSR, ITUSR, IAUSR, IBUSR
      REAL(DP) :: X, Y, T, B0, B1, B2, Z1, Z2
      REAL(DP), EXTERNAL :: RANF_EIRENE
      INTEGER :: IER

      entry sm0usr (is1,is2,sorad1,sorad2,sorad3,sorad4,sorad5,sorad6)
      return

      entry SM1USR (NLSF,X0,Y0,Z0,
     .              SORAD1,SORAD2,SORAD3,SORAD4,SORAD5,SORAD6,
     .              IRUSR,IPUSR,ITUSR,IAUSR,IBUSR,
     .              TIWL,TEWL,DIWL,VXWL,VYWL,VZWL,EFWL,SHWL,WEISPZ)

      RETURN
      END
C ===== SOURCE: sigusr.f


      SUBROUTINE SIGUSR(IFIRST,JJJ,ZDS,DUMMY1,PSIG,DUMMY2,ARGST,
     .                  XD0,YD0,ZD0,XD1,YD1,ZD1)
C
C  INPUT:
C          IFIRST: FLAG FOR INITIALISATION
C          NCELL:   INDEX IN TALLY ARRAYS FOR CURRENT ZONE
C          JJJ:    INDEX OF SEGMENT ALONG CHORD
C          ZDS:    LENGTH OF SEGMENT NO. JJJ
C  OUTPUT: CONTRIB. FROM CELL NCELL AND CHORD SEGMENT JJJ TO:
C          THE H ALPHA FLUX PSIG(I),I=0,4 CONTRIBUTIONS
C          FROM ATOMS, MOLECULES, TEST IONS AND BULK IONS
C          THE INTEGRANT ARGST SUCH THAT INTEGR.(ARGST*DL) = PSIG
C
      USE PRECISION
      USE PARMMOD
      USE CESTIM
      USE COMUSR
      USE COMPRT
      USE CSPEI
      USE COMSOU
      USE CLOGAU
      USE COMXS
      USE COMSIG
      USE CUPD
      USE CZT1
      USE CTRCEI
      USE CGRID
      USE CCONA
      USE CADGEO
      USE CLGIN
      USE CSDVI
      USE COUTAU
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: JJJ
      INTEGER, INTENT(INOUT) :: IFIRST
      REAL(DP), INTENT(INOUT) :: PSIG(0:NSPZ+10),ARGST(0:NSPZ+10,NRAD)
      REAL(DP), INTENT(IN) :: ZDS,DUMMY1,DUMMY2,XD0,YD0,ZD0,XD1,YD1,ZD1
      RETURN
      END
C ===== SOURCE: talusr.f
c
c
      subroutine talusr (ICOUNT,VECTOR,TALTOT,TALAV,
     .              TXTTL,TXTSP,TXTUN,ILAST,*)
      USE PRECISION
      USE PARMMOD
      USE CGRID
      USE CGEOM
      implicit NONE
      integer, intent(in) :: icount
      integer, intent(out) :: ilast
      real(dp), intent(in) :: vector(*), TALTOT, TALAV
      character(len=*) :: txttl,txtsp,txtun
      integer :: i

      ilast=1
      return 1
      end
C ===== SOURCE: timusr.f


      SUBROUTINE TIMUSR(N,X,Y,Z,VX,VY,VZ,N1,N2,T,IC,IE,NP,NL)
      USE PRECISION
      USE PARMMOD
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: X,Y,Z,VX,VY,VZ,T,cx,cy,cz,sc
      INTEGER, INTENT(IN) :: N, N1, N2, IC, IE, NP, IS, NRCELL
      LOGICAL :: NL

      ENTRY NORUSR(is,x,y,z,cx,cy,cz,sc,VX,VY,VZ,NRCELL)

      RETURN
      END
C ===== SOURCE: tmsusr.f


      SUBROUTINE TMSUSR (T0)
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: T0
      RETURN
      END
C ===== SOURCE: upcusr.f
C
C
      SUBROUTINE UPCUSR(WS,IND)
C
C  USER SUPPLIED COLLISION ESTIMATOR, VOLUME AVERAGED
C
      USE PRECISION
      USE PARMMOD
      USE CESTIM
      USE COMUSR
      USE COMPRT
      USE COMXS
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: WS
      INTEGER, INTENT(IN) :: IND
C
C     WS=WEIGHT/SIGTOT=WEIGHT/(VEL*ZMFPI)=WEIGHT/(VEL*SIGMA,MACR.)
C
C  FOR PARTICLE DENSITY IN CELL NO. NCELL
C     COLV(1,NCELL)=COLV(1,NCELL)+WS
C
      RETURN
      END
C ===== SOURCE: upnusr.f
c
c
      subroutine upnusr
      return
      end
C ===== SOURCE: upsusr.f
C
C
      SUBROUTINE UPSUSR(WT,IND)
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: WT
      INTEGER, INTENT(IN) :: IND
      RETURN
      END
C ===== SOURCE: uptusr.f
! 23.08.06: VPX, VPY, VRX, VRY changed to ALLOCATABLE, SAVE to speed up
!           subroutine call (save time in storage allocation)
C
C
      SUBROUTINE UPTUSR(XSTOR2,XSTORV2,WV,IFLAG)
C
C  USER SUPPLIED TRACKLENGTH ESTIMATOR, VOLUME AVERAGED
C
      USE PRECISION
      USE PARMMOD
      USE CESTIM
      USE COMUSR
      USE COMPRT
      USE CUPD
      USE COMXS
      USE CSPEZ
      USE CGRID
      USE CLOGAU
      USE CCONA
      USE CPOLYG
      USE CZT1
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: XSTOR2(MSTOR1,MSTOR2,N2ND+N3RD),
     .                         XSTORV2(NSTORV,N2ND+N3RD), WV
      INTEGER, INTENT(IN) :: IFLAG
      REAL(DP), ALLOCATABLE, SAVE :: CNDYNA(:),CNDYNP(:)
CDR
      REAL(DP), ALLOCATABLE, SAVE :: VPX(:),VPY(:),VRX(:),VRY(:)
CDR
      INTEGER :: IAT, IPL, I, IR, IP, IRD
      INTEGER, SAVE :: IFIRST, IA1, IA2, IA3, NA4, INDEXM, INDEXF
      DATA IFIRST/0/

      IF (IFIRST.EQ.0) THEN
        IFIRST=1
        ALLOCATE (CNDYNA(NATM))
        ALLOCATE (CNDYNP(NPLS))
        DO IAT=1,NATMI
          CNDYNA(IAT)=1.D3*AMUA*RMASSA(IAT)
        END DO
        DO IPL=1,NPLSI
          CNDYNP(IPL)=1.D3*AMUA*RMASSP(IPL)
        END DO
C
CDR
CDR  PROVIDE A RADIAL UNIT VECTOR PER CELL
CDR  VPX,VPY,  NEEDED FOR PROJECTING PARTICLE VELOCITIES
CDR  SAME FOR POLOIDAL UNIT VECTOR VRX,VRY
C
        ALLOCATE (VPX(NRAD))
        ALLOCATE (VPY(NRAD))
        ALLOCATE (VRX(NRAD))
        ALLOCATE (VRY(NRAD))
        DO I=1,NRAD
          VPX(I)=0.
          VPY(I)=0.
          VRX(I)=0.
          VRY(I)=0.
        END DO
        DO IR=1,NR1STM
          DO IP=1,NP2NDM
            IRD=IR+(IP-1)*NR1P2
            VPX(IRD)=PLNX(IR,IP)
            VPY(IRD)=PLNY(IR,IP)
            VRX(IRD)=PPLNX(IR,IP)
            VRY(IRD)=PPLNY(IR,IP)
          END DO
        END DO
        IA1=NATMI+NMOLI
        IA2=2*IA1
        IA3=3*IA1
        NA4=4*IA1
        INDEXM=NPLSI
        INDEXF=2*NPLSI
      ENDIF
      RETURN
      END



C ===== SOURCE: vdion.f


      FUNCTION VDION (I)
      USE PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I
      REAL(DP) :: VDION
      VDION=0.
      RETURN
      END
C ===== SOURCE: vecusr.f


      SUBROUTINE VECUSR (I,VX,VY,VZ,IPLS)
      USE PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I, IPLS
      REAL(DP), INTENT(IN) :: VX,VY,VZ
      RETURN
      END
C ===== SOURCE: volusr.f


      SUBROUTINE VOLUSR(N,A)
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: A(*)
      INTEGER, INTENT(IN) :: N
      RETURN
      END
