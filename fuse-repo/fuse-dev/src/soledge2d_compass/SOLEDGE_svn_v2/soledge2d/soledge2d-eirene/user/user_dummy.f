
*DK USER
C
C   USER SUPPLIED SUBROUTINES
C
C           ***********
C           *  TORE   *
C           ***********
C

C
C


      SUBROUTINE EIRENE_GEOUSR

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CADGEO
      USE EIRMOD_CLGIN
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_CGRPTL
      USE EIRMOD_COMUSR

      IMPLICIT NONE


      RETURN
      END
C
      SUBROUTINE Eirene_PROUSR (PRO,INDX,P0,P1,P2,P3,P4,P5,PROVAC,N)
C
      USE Eirmod_PRECISION
      USE Eirmod_PARMMOD
      USE Eirmod_COMUSR
      USE Eirmod_CGRID
      USE Eirmod_CGEOM
      USE Eirmod_CINIT
      USE Eirmod_CCONA
      
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: P0, P1, P2, P3, P4, P5, PROVAC
      REAL(DP), INTENT(OUT) :: PRO(*)
      INTEGER, INTENT(IN) :: INDX, N

      RETURN
      END


      SUBROUTINE EIRENE_PLAUSR
      IMPLICIT NONE
      RETURN
      END

      SUBROUTINE eirene_PLTUSR(PLABLE,J)
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: PLABLE
      INTEGER, INTENT(IN) :: J
      RETURN
      END
C
C
      SUBROUTINE EIRENE_SAMUSR (NLSF,XX0,YY0,ZZ0,
     .              SORD1,SORD2,SORD3,SORD4,SORD5,SORD6,
     .              IRUSR,IPUSR,ITUSR,IAUSR,IBUSR,
     .              TIWL,TEWL,DIWL,VXWL,VYWL,VZWL,WEISPZ)
C
C  SAMPLE INITAL COORDIANTES X,Y,Z ON ADDITIONAL SURFACE NLLI
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CADGEO
      USE EIRMOD_CGRID
      USE EIRMOD_COMUSR
      USE EIRMOD_COMSOU
      USE EIRMOD_COMPRT
      USE EIRMOD_CCONA
      USE EIRMOD_CESTIM
      
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: SORD1,SORD2,SORD3,SORD4,SORD5,SORD6
      REAL(DP) :: XX0,YY0,ZZ0,TEWL,TIWL(*),DIWL(*),
     .                       VXWL(*),VYWL(*),VZWL(*),WEISPZ(*),
     .                       EFWL(*),SHWL
      INTEGER, INTENT(IN) :: NLSF,is1, is2
      INTEGER :: IRUSR, IPUSR, ITUSR, IAUSR, IBUSR
      REAL(DP) :: X, Y, T, B0, B1, B2, Z1, Z2, zep1, rad, pphi, wink,
     .            aout, ain, aslice, frac, tot_area
      real(dp), allocatable :: dummy(:)
      REAL(DP), EXTERNAL :: RANF_EIRENE
      integer, external :: eirene_learca, eirene_learc1
      INTEGER :: IER, i, irc, itc, icell, ic, il, im, iu, nnt

      entry eirene_sm0usr (is1,is2,sord1,sord2,sord3,sord4,sord5,sord6)

      
      return

      entry eirene_SM1USR (NLSF,XX0,YY0,ZZ0,
     .              SORD1,SORD2,SORD3,SORD4,SORD5,SORD6,
     .              IRUSR,IPUSR,ITUSR,IAUSR,IBUSR,
     .              TIWL,TEWL,DIWL,VXWL,VYWL,VZWL,EFWL,SHWL,WEISPZ)
      

      RETURN
      END
C
C
C
      SUBROUTINE EIRENE_UPTUSR(XSTOR2,XSTORV2,WV,IFLAG)
C
C  USER SUPPLIED TRACKLENGTH ESTIMATOR, VOLUME AVERAGED
C
C
CC
C  WV=WEIGHT/VEL
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMXS
      use eirmod_comprt
      use eirmod_cgeom
      use eirmod_cestim
      use eirmod_cupd
      use eirmod_cgrid
      use eirmod_comusr

      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: XSTOR2(MSTOR1,MSTOR2,N2ND+N3RD),
     .                         XSTORV2(NSTORV,N2ND+N3RD), WV
      INTEGER, INTENT(IN) :: IFLAG
      integer :: i,ird,irdo,iacx,ipl1,ipl2,ircx
      integer :: imcx,iicx,iplst
      integer :: iaei,irei,ip,imei,iiei
      real(dp):: dist,wtr,wtrsig

      select case (ityp)

        case (1)  
 
          do I=1,NCOU
            DIST=CLPD(I)
            WTR=WV*DIST
            IRDO=NRCELL+NUPC(I)*NR1P2+NBLCKA
            IRD=NCLTAL(IRDO)

c EI contribution
            do iaei=1,naeii(iatm)
              irei=lgaei(iatm,iaei)
              wtrsig=wtr*sigvei(irei)
              do ip=1,ipplds(irei,0)
                ipls=ipplds(irei,ip)
                IF (LEAPL) addv(ipls,IRD)=addv(ipls,IRD)+
     .                          WTRSIG*ESIGEI(IREI,4)
              enddo
            enddo



c CX contribution

            do iacx=1,nacxi(iatm)
c    need to know which ion ipls reacts and which ion appears
              IRCX=LGACX(IATM,IACX,0)
c pre collision ion
              IPLS=LGACX(IATM,IACX,1)
              WTRSIG=WTR*SIGVCX(IRCX)

              if (LEAPL) addv(ipls,IRD)   = addv(ipls,IRD)
     .                    - WTRSIG*ESIGCX(IRCX,1) 

c post collision (2 secondaries)
              if (N1STX(IRCX,1).EQ.4) THEN
                IPL1=N1STX(IRCX,2)
                if (LEAPL) then
                  addv(IPL1,IRD)= addv(IPL1,IRD)+WTRSIG*ESIGCX(IRCX,1)
                endif
              endif

              if (N2NDX(IRCX,1).EQ.4) THEN
                IPL2=N2NDX(IRCX,2)
                IF (LEAPL) THEN
                  addv(IPL2,IRD)= addv(IPL2,IRD)+WTRSIG*E0
                END IF
              endif 
            enddo
          enddo

        case(2)

          do I=1,NCOU
            DIST=CLPD(I)
            WTR=WV*DIST
            IRDO=NRCELL+NUPC(I)*NR1P2+NBLCKA
            IRD=NCLTAL(IRDO)

c        EI contribution
            do imei=1,nmdsi(imol)
              irei=lgmei(imol,imei)
              wtrsig=wtr*sigvei(irei)
              do ip=1,ipplds(irei,0)
                iplst=ipplds(irei,ip)+nplsi                              
                IF (LEMPL) addv(iplst,IRD)=addv(iplst,IRD)+
     .                          WTRSIG*ESIGEI(IREI,4)
              enddo
            enddo

c         CX contribution

            do imcx=1,nmcxi(imol)
c    need to know which ion ipls reacts and which ion appears
              IRCX=LGMCX(IMOL,IMCX,0)
c pre collision ion
              IPLSt=LGMCX(IMOL,IMCX,1)+nplsi
              WTRSIG=WTR*SIGVCX(IRCX)
              if (LEMPL) addv(iplst,IRD)   = addv(iplst,IRD)
     .                    - WTRSIG*ESIGCX(IRCX,1) 

c post collision (2 secondaries)
              if (N1STX(IRCX,1).EQ.4) THEN
                IPL1=N1STX(IRCX,2)+nplsi
                if (LEMPL) then
                  addv(IPL1,IRD)= addv(IPL1,IRD)+WTRSIG*ESIGCX(IRCX,1)
                endif
              endif

              if (N2NDX(IRCX,1).EQ.4) THEN
                IPL2=N2NDX(IRCX,2)+nplsi
                IF (LEMPL) THEN
                  addv(IPL2,IRD)= addv(IPL2,IRD)+WTRSIG*E0
                END IF
              endif 
            enddo
          enddo

        case(3)

          do I=1,NCOU
            DIST=CLPD(I)
            WTR=WV*DIST
            IRDO=NRCELL+NUPC(I)*NR1P2+NBLCKA
            IRD=NCLTAL(IRDO)

c        EI contribution

            do iiei=1,nidsi(iion)
              irei=lgiei(iion,iiei)
              wtrsig=wtr*sigvei(irei)
              do ip=1,ipplds(irei,0)
                iplst=ipplds(irei,ip)+2*nplsi
                IF (LEIPL) addv(iplst,IRD)=addv(iplst,IRD)+
     .                      WTRSIG*ESIGEI(IREI,4)
              enddo
            enddo

c         CX contribution

            do iicx=1,nicxi(imol)
c    need to know which ion ipls reacts and which ion appears
              IRCX=LGICX(IION,IMCX,0)
c pre collision ion
              IPLSt=LGICX(IION,IMCX,1)+2*nplsi
              WTRSIG=WTR*SIGVCX(IRCX)
              if (LEIPL) addv(iplst,IRD)   = addv(iplst,IRD)
     .                    - WTRSIG*ESIGCX(IRCX,1) 

c post collision (2 secondaries)
              if (N1STX(IRCX,1).EQ.4) THEN
                IPL1=N1STX(IRCX,2)+2*nplsi
                if (LEIPL) then
                  addv(IPL1,IRD)= addv(IPL1,IRD)+WTRSIG*ESIGCX(IRCX,1)
                endif
              endif

              if (N2NDX(IRCX,1).EQ.4) THEN
                IPL2=N2NDX(IRCX,2)+2*nplsi
                IF (LEIPL) THEN
                  addv(IPL2,IRD)= addv(IPL2,IRD)+WTRSIG*E0
                END IF
              endif 
            enddo
          enddo

      end select

      

      RETURN
      END
C
C
      SUBROUTINE EIRENE_UPSUSR(WPR,IND)
C  USER SUPPLIED ESTIMATOR, SURFACE AVERAGED
C  (COLLISION - AND TRACKLENGTH ESTIMATORS ARE IDENTICAL, IF SURFACE
C   AVERAGED)
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CSPEZ
      USE EIRMOD_CUPD
      USE EIRMOD_CGRID
      USE EIRMOD_CZT1
      USE EIRMOD_CLOGAU
      USE EIRMOD_COMXS
      USE EIRMOD_COMPRT
      USE EIRMOD_CCONA
      USE EIRMOD_CGEOM
      USE EIRMOD_CGRPTL
      USE EIRMOD_CESTIM
      USE EIRMOD_COMUSR
      USE EIRMOD_CLGIN

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: WPR
      INTEGER, INTENT(IN) :: IND
 
      RETURN
      END
C
C
      SUBROUTINE EIRENE_UPCUSR(WS,IND)
C
C  USER SUPPLIED COLLISION ESTIMATOR, VOLUME AVERAGED
C
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: WS
      INTEGER, INTENT(IN) :: IND
C
C     WS=WEIGHT/SIGTOT=WEIGHT/(VEL*ZMFPI)=WEIGHT/(VEL*SIGMA,MACR.)
C
      RETURN
      END
C
C
      SUBROUTINE EIRENE_CRVUSR (ILIMI,TA,TB,IND,P,Q)
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: P(*), Q(*)
      REAL(DP), INTENT(IN) :: T, TA, TB
      REAL(DP), INTENT(OUT) :: X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ
      INTEGER, INTENT(IN) :: ILIMI, IND
      ENTRY PNTUSR(T,X,Y,Z,IND,P,Q)
        X=5.+T
        Y=20.
        Z=2.5
      RETURN
      ENTRY DPTUSR(T,DX,DY,DZ,IND,P,Q)
        DX=1.
        DY=0.
        DZ=0.
      RETURN
      ENTRY DDPUSR(T,DDX,DDY,DDZ,IND,P,Q)
        DDX=0.
        DDY=0.
        DDZ=0.
      RETURN
      END

      SUBROUTINE EIRENE_sigusr(ifirst,JJJ,ZDS,PEN,PSIG,TIMAX,ARGST,
     .   XD0,YD0,ZD0,
     .   XD1,YD1,ZD1)
      use eirmod_precision
      use eirmod_parmmod
      implicit none
      integer, intent(in) :: ifirst, jjj
      real(dp), intent(in) :: pen,zds
      REAL(DP), INTENT(IN OUT) :: PSIG(0:NSPZ+10), TIMAX
      REAL(DP), INTENT(IN OUT) :: ARGST(0:NSPZ+10,NRAD)
      real(dp), intent(in) :: XD0,YD0,ZD0,XD1,YD1,ZD1

      return
      end

      SUBROUTINE EIRENE_upnusr
      use eirmod_precision
      use eirmod_parmmod
      use eirmod_cestim
      use eirmod_comprt

c snapshot tally for density
      snapv(1,ncell) = snapv(1,ncell)+WEIGHT


      return
      end
c
c
      subroutine EIRENE_talusr (ICOUNT,VECTOR,TALTOT,TALAV,
     .              TXTTL,TXTSP,TXTUN,ILAST,*)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      implicit NONE
      integer, intent(in) :: icount
      integer, intent(out) :: ilast
      real(dp), intent(in) :: vector(*), TALTOT, TALAV
      character(len=*) :: txttl,txtsp,txtun
      integer :: i
 
      ilast=1
      return 1
      end

c
c     SUBROUTINE EIRENE_diagno
c     return
c     end
c
c
      SUBROUTINE EIRENE_modusr
      return
      end
c
c
      SUBROUTINE EIRENE_retusr(sig)
      USE EIRMOD_PRECISION
      implicit none
      real(dp), intent(in) :: sig
      return
      end
C
C
      SUBROUTINE EIRENE_REFUSR
C
C   USER SUPPLIED REFLECTION MODEL
C
      USE EIRMOD_PRECISION
      use eirmod_clgin
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: XMW,XCW,XMP,XCP,ZCOS,ZSIN,EXPI,RPROB,
     .                        E0TERM
      INTEGER, INTENT(IN) :: IGASF,IGAST

      real :: stick_C
      logical :: particle_mod

      ENTRY eirene_RF0USR
      ENTRY eirene_SPTUSR
      ENTRY eirene_SP0USR
 



      ENTRY eirene_SP1USR
      ENTRY eirene_RF1USR (XMW,XCW,XMP,XCP,IGASF,IGAST,ZCOS,ZSIN,EXPI,
     .              RPROB,E0TERM,*,*,*,*)
      RETURN
      END
C
C
      function eirene_leausr(a,b,c)
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: A, B, C
      INTEGER :: eirene_LEAUSR
      eirene_leausr=1
      END


      SUBROUTINE EIRENE_TIMUSR(N,X,Y,Z,VX,VY,VZ,N1,N2,T,IC,IE,NP,NL)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: X,Y,Z,VX,VY,VZ,T,cx,cy,cz,sc
      INTEGER, INTENT(IN) :: N, N1, N2, IC, IE, NP, IS, NRCELL
      LOGICAL :: NL

      ENTRY eirene_NORUSR(is,x,y,z,cx,cy,cz,sc,VX,VY,VZ,NRCELL)

      RETURN
      END


      SUBROUTINE EIRENE_VOLUSR(N,A)
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: A(*)
      INTEGER, INTENT(IN) :: N
      RETURN
      END


      SUBROUTINE EIRENE_VECUSR (I,VX,VY,VZ,IPLS)
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I, IPLS
      REAL(DP), INTENT(IN) :: VX,VY,VZ
      RETURN
      END


c      SUBROUTINE EIRENE_INIusr
c      IMPLICIT NONE
c      RETURN
c      END


      SUBROUTINE eirene_TMSUSR (T0)
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: T0
      RETURN
      END


      function eirene_vdion(i)
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I
      REAL(DP) :: eirene_VDION
      eirene_vdion=0._dp


      end function
