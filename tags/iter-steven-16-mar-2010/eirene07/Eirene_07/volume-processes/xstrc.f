!pb  30.08.06: data structure for reaction data redefined
!pb  12.10.06: modcol revised

      SUBROUTINE XSTRC(ipls,nrc,idsc,irrc)
cdr
cdr  to replace ph_xsectp, as called from XSECTP.F,
cdr  prepare volume recombination processes: bulk (+ bulk)--> test (+ bulk)
cdr  i.e.                                    H+    +  e   --> H    (+ rad.)
cdr  i.e.                                    H(n=2)       --> Ly-alpha (+H(n=1))
cdr
c    ipls: incident bulk
c    nrc : index of reaction in list of all reactions for IPLS
c    idsc: index of RC reaction for species ipls
c    irrc: index of RC reaction in NRRC arrays
c
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT
      USE COMXS
      USE PHOTON
      IMPLICIT NONE
      integer, intent(in) :: ipls,nrc,idsc,irrc
      integer :: kk,ipl0,ipl1,ipl2,ityp0,ityp1,ityp2,
     .           j
      integer, external :: idez
      real(dp) :: factkk,aik

c  fetch data for process nrc

      kk = IREACP(ipls,nrc)
      if (kk /= idreac) call get_reaction(kk)

      FACTKK=FREACP(IPLS,NRC)
      IF (FACTKK.EQ.0.D0) FACTKK=1.
      aik=reaction%aik


      IPL0 =IDEZ(IBULKP(ipls,nrc),3,3)
      IPL1 =IDEZ(ISCD1P(ipls,nrc),3,3)
      IPL2 =IDEZ(ISCD2P(ipls,nrc),3,3)
      ITYP0=IDEZ(IBULKP(ipls,nrc),1,3)
      ITYP1=IDEZ(ISCD1P(ipls,nrc),1,3)
      ITYP2=IDEZ(ISCD2P(ipls,nrc),1,3)

      LGPRC(IPLS,IDSC)=IRRC

      facrea(kk) = factkk
      NREARC(irrc) = kk
      do j=1,nrad
         tabrc1(irrc,j)=aik*factkk
      enddo
      modcol(6,2,irrc)=1

      select case(ityp1)
      case(0)
         NPHPRC(IRRC)=IPL1
      case(1)
         NATPRC(IRRC)=IPL1
      case(2)
         NMLPRC(IRRC)=IPL1
      case(3)
         NIOPRC(IRRC)=IPL1
      case(4)
         NPLPRC(IRRC)=IPL1
      case default
         write(iunout,*) 'volume-processes: xstrc.f:'
         write(iunout,*) '   ityp1=',ityp1,' not allowed'
         call exit(1)
      end select
      select case(ityp2)
      case(0)
         NPHPRC_2(IRRC)=IPL2
      case(1)
         NATPRC_2(IRRC)=IPL2
      case(2)
         NMLPRC_2(IRRC)=IPL2
      case(3)
         NIOPRC_2(IRRC)=IPL2
      case(4)
         NPLPRC_2(IRRC)=IPL2
      case default
         write(iunout,*) 'volume-processes: xstrc.f:'
         write(iunout,*) '   ityp2=',ityp2,' not allowed'
         call exit(1)
      end select
      return
      END SUBROUTINE XSTRC
