 

      SUBROUTINE SUMSTRA

      USE PRECISION
      USE PARMMOD
      use CESTIM
      use COMUSR
      USE CGRID
      USE CSPEZ
      USE CPES
      USE COMSOU
      USE CSDVI
      USE CSPEI
      USE COUTAU
      IMPLICIT NONE

      include 'mpif.h'
      integer :: ier, ier1, isdv, igrp, my_grp_rnk, icomgrp, mpicw
C

      call mpi_barrier(mpi_comm_world,ier)

      call mpi_comm_group (mpi_comm_world,mpicw,ier)
      call mpi_group_incl (mpicw,nstrai,npesta(1:nstrai),igrp,ier)
      call mpi_comm_create (mpi_comm_world,igrp,icomgrp,ier)
      CALL MPI_GROUP_RANK(IGRP,MY_GRP_RNK,ier)

      IF (MY_GRP_RNK .NE. MPI_UNDEFINED) THEN
        call mpi_reduce(PPMLI(0,0),PPMLI(0,0),NMOLI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(WTOTM(0,0),WTOTM(0,0),NMOLI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PDENMI(0,0),PDENMI(0,0),NMOLI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EDENMI(0,0),EDENMI(0,0),NMOLI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PAMLI(0,0),PAMLI(0,0),NMOLI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PMMLI(0,0),PMMLI(0,0),NMOLI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PIMLI(0,0),PIMLI(0,0),NMOLI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PPHMLI(0,0),PPHMLI(0,0),NMOLI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(POTMLI(0,0),POTMLI(0,0),NMOLI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PRFAMI(0,0),PRFAMI(0,0),NMOLI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PRFMMI(0,0),PRFMMI(0,0),NMOLI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PRFIMI(0,0),PRFIMI(0,0),NMOLI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PRFPHMI(0,0),PRFPHMI(0,0),NMOLI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PRFPMI(0,0),PRFPMI(0,0),NMOLI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EOTMLI(0,0),EOTMLI(0,0),NMOLI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(ERFAMI(0,0),ERFAMI(0,0),NMOLI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(ERFMMI(0,0),ERFMMI(0,0),NMOLI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(ERFIMI(0,0),ERFIMI(0,0),NMOLI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(ERFPHMI(0,0),ERFPHMI(0,0),NMOLI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(ERFPMI(0,0),ERFPMI(0,0),NMOLI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(SPTMLI(0,0),SPTMLI(0,0),NMOLI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PGENMI(0,0),PGENMI(0,0),NMOLI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EGENMI(0,0),EGENMI(0,0),NMOLI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(VGENMI(0,0),VGENMI(0,0),NMOLI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)


        CALL MPI_REDUCE(PPATI(0,0),PPATI(0,0),NATMI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(WTOTA(0,0),WTOTA(0,0),NATMI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PDENAI(0,0),PDENAI(0,0),NATMI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EDENAI(0,0),EDENAI(0,0),NATMI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PAATI(0,0),PAATI(0,0),NATMI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PMATI(0,0),PMATI(0,0),NATMI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PIATI(0,0),PIATI(0,0),NATMI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PPHATI(0,0),PPHATI(0,0),NATMI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(POTATI(0,0),POTATI(0,0),NATMI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PRFAAI(0,0),PRFAAI(0,0),NATMI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PRFMAI(0,0),PRFMAI(0,0),NATMI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PRFIAI(0,0),PRFIAI(0,0),NATMI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PRFPHAI(0,0),PRFPHAI(0,0),NATMI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PRFPAI(0,0),PRFPAI(0,0),NATMI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EOTATI(0,0),EOTATI(0,0),NATMI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(ERFAAI(0,0),ERFAAI(0,0),NATMI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(ERFMAI(0,0),ERFMAI(0,0),NATMI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(ERFIAI(0,0),ERFIAI(0,0),NATMI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(ERFPHAI(0,0),ERFPHAI(0,0),NATMI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(ERFPAI(0,0),ERFPAI(0,0),NATMI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(SPTATI(0,0),SPTATI(0,0),NATMI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PGENAI(0,0),PGENAI(0,0),NATMI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EGENAI(0,0),EGENAI(0,0),NATMI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(VGENAI(0,0),VGENAI(0,0),NATMI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)


        CALL MPI_REDUCE(PPIOI(0,0),PPIOI(0,0),NIONI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EELFI(0,0),EELFI(0,0),NIONI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(WTOTI(0,0),WTOTI(0,0),NIONI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PDENII(0,0),PDENII(0,0),NIONI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EDENII(0,0),EDENII(0,0),NIONI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PAIOI(0,0),PAIOI(0,0),NIONI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PMIOI(0,0),PMIOI(0,0),NIONI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PIIOI(0,0),PIIOI(0,0),NIONI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PPHIOI(0,0),PPHIOI(0,0),NIONI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(POTIOI(0,0),POTIOI(0,0),NIONI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PRFAII(0,0),PRFAII(0,0),NIONI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PRFMII(0,0),PRFMII(0,0),NIONI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PRFIII(0,0),PRFIII(0,0),NIONI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PRFPHII(0,0),PRFPHII(0,0),NIONI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PRFPII(0,0),PRFPII(0,0),NIONI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EOTIOI(0,0),EOTIOI(0,0),NIONI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(ERFAII(0,0),ERFAII(0,0),NIONI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(ERFMII(0,0),ERFMII(0,0),NIONI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(ERFIII(0,0),ERFIII(0,0),NIONI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(ERFPHII(0,0),ERFPHII(0,0),NIONI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(ERFPII(0,0),ERFPII(0,0),NIONI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(SPTIOI(0,0),SPTIOI(0,0),NIONI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PGENII(0,0),PGENII(0,0),NIONI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EGENII(0,0),EGENII(0,0),NIONI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(VGENII(0,0),VGENII(0,0),NIONI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
 

        CALL MPI_REDUCE(PPPHTI(0,0),PPPHTI(0,0),NPHOTI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(WTOTPH(0,0),WTOTPH(0,0),NPHOTI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PDENPHI(0,0),PDENPHI(0,0),NPHOTI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EDENPHI(0,0),EDENPHI(0,0),NPHOTI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PAPHTI(0,0),PAPHTI(0,0),NPHOTI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PMPHTI(0,0),PMPHTI(0,0),NPHOTI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PIPHTI(0,0),PIPHTI(0,0),NPHOTI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PPHPHTI(0,0),PPHPHTI(0,0),NPHOTI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(POTPHTI(0,0),POTPHTI(0,0),NPHOTI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PRFAPHTI(0,0),PRFAPHTI(0,0),NPHOTI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PRFMPHTI(0,0),PRFMPHTI(0,0),NPHOTI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
       CALL MPI_REDUCE(PRFIPHTI(0,0),PRFIPHTI(0,0),NPHOTI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
       CALL MPI_REDUCE(PRFPHPHTI(0,0),PRFPHPHTI(0,0),NPHOTI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PRFPPHTI(0,0),PRFPPHTI(0,0),NPHOTI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
       CALL MPI_REDUCE(EOTPHTI(0,0),EOTPHTI(0,0),NPHOTI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(ERFAPHTI(0,0),ERFAPHTI(0,0),NPHOTI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(ERFMPHTI(0,0),ERFMPHTI(0,0),NPHOTI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(ERFIPHTI(0,0),ERFIPHTI(0,0),NPHOTI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(ERFPHPHTI(0,0),ERFPHPHTI(0,0),NPHOTI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(ERFPPHTI(0,0),ERFPPHTI(0,0),NPHOTI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(SPTPHTI(0,0),SPTPHTI(0,0),NPHOTI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PGENPHI(0,0),PGENPHI(0,0),NPHOTI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EGENPHI(0,0),EGENPHI(0,0),NPHOTI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(VGENPHI(0,0),VGENPHI(0,0),NPHOTI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)


        CALL MPI_REDUCE(PPPLI(0,0),PPPLI(0,0),NPLSI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(WTOTP(0,0),WTOTP(0,0),NPLSI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PAPLI(0,0),PAPLI(0,0),NPLSI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PMPLI(0,0),PMPLI(0,0),NPLSI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PIPLI(0,0),PIPLI(0,0),NPLSI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PPHPLI(0,0),PPHPLI(0,0),NPLSI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(POTPLI(0,0),POTPLI(0,0),NPLSI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EOTPLI(0,0),EOTPLI(0,0),NPLSI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(SPTPLI(0,0),SPTPLI(0,0),NPLSI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)


        CALL MPI_REDUCE(ADDVI(0,0),ADDVI(0,0),NADVI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)

        CALL MPI_REDUCE(COLVI(0,0),COLVI(0,0),NCLVI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)

        CALL MPI_REDUCE(SNAPVI(0,0),SNAPVI(0,0),NSNVI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)

        CALL MPI_REDUCE(COPVI(0,0),COPVI(0,0),NCPVI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)

        CALL MPI_REDUCE(BGKVI(0,0),BGKVI(0,0),NBGVI+1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)

        CALL MPI_REDUCE(PAELI(0),PAELI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PMELI(0),PMELI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PIELI(0),PIELI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PPHELI(0),PPHELI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
C
        CALL MPI_REDUCE(EAELI(0),EAELI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EAATI(0),EAATI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EAMLI(0),EAMLI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EAIOI(0),EAIOI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EAPHTI(0),EAPHTI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EAPLI(0),EAPLI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
C
        CALL MPI_REDUCE(EMELI(0),EMELI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EMATI(0),EMATI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EMMLI(0),EMMLI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EMIOI(0),EMIOI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EMPHTI(0),EMPHTI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EMPLI(0),EMPLI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
C
        CALL MPI_REDUCE(EIELI(0),EIELI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EIATI(0),EIATI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EIMLI(0),EIMLI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EIIOI(0),EIIOI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EIPHTI(0),EIPHTI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EIPLI(0),EIPLI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
C
        CALL MPI_REDUCE(EPHELI(0),EPHELI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EPHATI(0),EPHATI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EPHMLI(0),EPHMLI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EPHIOI(0),EPHIOI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EPHPHTI(0),EPHPHTI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EPHPLI(0),EPHPLI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
C
        CALL MPI_REDUCE(EPATI(0),EPATI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EPMLI(0),EPMLI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EPIOI(0),EPIOI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EPPHTI(0),EPPHTI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(EPPLI(0),EPPLI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
C
C
        CALL MPI_REDUCE(FLUXT(0),FLUXT(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(XMCP(0),XMCP(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(PTRASH(0),PTRASH(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(ETRASH(0),ETRASH(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(ETOTA(0),ETOTA(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(ETOTM(0),ETOTM(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(ETOTI(0),ETOTI(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(ETOTPH(0),ETOTPH(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(ETOTP(0),ETOTP(0),1,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)

C
C
        CALL MPI_REDUCE(SMESTV,SMESTV,NIDV*NRTAL,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(SMESTS,SMESTS,NIDS*NLMPGS,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
 
        IF (NALVI.GT.0) THEN
          CALL MPI_REDUCE(ALGVI(0,0),ALGVI(0,0),NALVI+1,
     .                    mpi_real8,mpi_sum,0,icomgrp,ier1)
        ENDIF
C
        IF (NALSI.GT.0) THEN
          CALL MPI_REDUCE(ALGSI(0,0),ALGSI(0,0),NALSI+1,
     .                    mpi_real8,mpi_sum,0,icomgrp,ier1)
        ENDIF
 
        DO ISDV=1,NSIGCI
          CALL MPI_REDUCE(STVC(0,ISDV,1:NSBOX_TAL),
     .                    STVC(0,ISDV,1:NSBOX_TAL),NSBOX_TAL,mpi_real8,
     .                    mpi_sum,0,icomgrp,ier1)
          CALL MPI_REDUCE(STVC(1,ISDV,1:NSBOX_TAL),
     .                    STVC(1,ISDV,1:NSBOX_TAL),NSBOX_TAL,mpi_real8,
     .                    mpi_sum,0,icomgrp,ier1)
          CALL MPI_REDUCE(STVC(2,ISDV,1:NSBOX_TAL),
     .                    STVC(2,ISDV,1:NSBOX_TAL),NSBOX_TAL,mpi_real8,
     .                    mpi_sum,0,icomgrp,ier1)
        ENDDO   
        IF (NSIGCI > 0) THEN
          CALL MPI_REDUCE(STVCS(0,1:NSIGCI),STVCS(0,1:NSIGCI),NSIGCI,
     .                    mpi_real8,mpi_sum,0,icomgrp,ier1)
          CALL MPI_REDUCE(STVCS(1,1:NSIGCI),STVCS(1,1:NSIGCI),NSIGCI,
     .                    mpi_real8,mpi_sum,0,icomgrp,ier1)
          CALL MPI_REDUCE(STVCS(2,1:NSIGCI),STVCS(2,1:NSIGCI),NSIGCI,
     .                    mpi_real8,mpi_sum,0,icomgrp,ier1)
        END IF
     
        DO ISDV=1,NSIGVI
          CALL MPI_REDUCE(STV(ISDV,1:NSBOX_TAL),STV(ISDV,1:NSBOX_TAL),
     .                    NSBOX_TAL,MPI_REAL8,MPI_SUM,0,ICOMGRP,IER1)
          CALL MPI_REDUCE(EE(ISDV,1:NSBOX_TAL),EE(ISDV,1:NSBOX_TAL),
     .                    NSBOX_TAL,MPI_REAL8,MPI_SUM,0,ICOMGRP,IER1)
        ENDDO   
        IF (NSIGVI > 0) THEN
          CALL MPI_REDUCE(STVS,STVS,NSIGVI,
     .                    MPI_REAL8,MPI_SUM,0,ICOMGRP,IER1)
          CALL MPI_REDUCE(EE,EE,NSIGVI,
     .                    MPI_REAL8,MPI_SUM,0,ICOMGRP,IER1)
        END IF

     
        IF (NSIGSI > 0) THEN
          CALL MPI_REDUCE(STVW(1:NSIGSI,1:NLIMPS),
     .                    STVW(1:NSIGSI,1:NLIMPS),NSIGSI*NLIMPS,
     .                    MPI_REAL8,MPI_SUM,0,ICOMGRP,IER1)
          CALL MPI_REDUCE(FF(1:NSIGSI,1:NLIMPS),
     .                    FF(1:NSIGSI,1:NLIMPS),NSIGSI*NLIMPS,
     .                    MPI_REAL8,MPI_SUM,0,ICOMGRP,IER1)
          CALL MPI_REDUCE(STVWS,STVWS,NSIGCI,
     .                    MPI_REAL8,MPI_SUM,0,ICOMGRP,IER1)
          CALL MPI_REDUCE(FFS,FFS,NSIGCI,
     .                    MPI_REAL8,MPI_SUM,0,ICOMGRP,IER1)
        END IF
    
        call mpi_reduce(LOGMOL,LOGMOL,NMOLI+1,
     .                  mpi_logical,mpi_LOR,0,icomgrp,ier1)
        call mpi_reduce(LOGATM,LOGATM,NATMI+1,
     .                  mpi_logical,mpi_LOR,0,icomgrp,ier1)
        call mpi_reduce(LOGION,LOGION,NIONI+1,
     .                  mpi_logical,mpi_LOR,0,icomgrp,ier1)
        call mpi_reduce(LOGPHOT,LOGPHOT,NPHOTI+1,
     .                  mpi_logical,mpi_LOR,0,icomgrp,ier1)
        call mpi_reduce(LOGPLS,LOGPLS,NPLSI+1,
     .                  mpi_logical,mpi_LOR,0,icomgrp,ier1)

        call mpi_barrier(icomgrp,ier)
        call mpi_group_free (igrp,ier)
        call mpi_comm_free (icomgrp,ier)
      END IF

      RETURN
      END
