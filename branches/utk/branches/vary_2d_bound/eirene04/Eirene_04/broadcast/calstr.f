

      SUBROUTINE CALSTR

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE CSPEZ
      USE COMPRT
      USE CPES
      USE CSDVI
      USE COUTAU
      IMPLICIT NONE

C
      INCLUDE 'mpif.h'
      integer :: igrp(0:nstra), icomgrp(0:nstra)
      integer :: ier1, ier, ir, npean, npeen, i, mpicw, ispc

      npean = npesta(istra)
      npeen = npesta(istra)+npestr(istra)-1
      call mpi_comm_group (mpi_comm_world,mpicw,ier)
      call mpi_group_incl (mpicw,npestr(istra),
     .   (/ (i,i=npean,npeen)/),igrp(istra),ier)
      call mpi_comm_create (mpi_comm_world,igrp(istra),
     .   icomgrp(istra),ier)
c
c   collect data from pe's belonging to one stratum
c
      if (istra.eq.nstrpe(my_pe)) then
        call mpi_barrier(icomgrp(istra),ier)
        call mpi_reduce(WTOTM(0,istra),WTOTM(0,istra),nmoli+1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
        call mpi_reduce(WTOTA(0,istra),WTOTA(0,istra),natmi+1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
        call mpi_reduce(WTOTI(0,istra),WTOTI(0,istra),nioni+1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
        call mpi_reduce(WTOTP(0,istra),WTOTP(0,istra),nplsi+1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
C
C
        call mpi_reduce(XMCP(istra),XMCP(istra),1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
        call mpi_reduce(PTRASH(istra),PTRASH(istra),1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
        call mpi_reduce(ETRASH(istra),ETRASH(istra),1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
        call mpi_reduce(ETOTA(istra),ETOTA(istra),1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
        call mpi_reduce(ETOTM(istra),ETOTM(istra),1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
        call mpi_reduce(ETOTI(istra),ETOTI(istra),1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
        call mpi_reduce(ETOTP(istra),ETOTP(istra),1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
 
        call mpi_reduce(PPMLI(0,istra),PPMLI(0,istra),nmoli+1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
        call mpi_reduce(PPATI(0,istra),PPATI(0,istra),natmi+1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
        call mpi_reduce(PPIOI(0,istra),PPIOI(0,istra),nioni+1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
        call mpi_reduce(PPPHTI(0,istra),PPPHTI(0,istra),nphoti+1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
        call mpi_reduce(PPPLI(0,istra),PPPLI(0,istra),nplsi+1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
 
        call mpi_reduce(EPATI(istra),EPATI(istra),1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
        call mpi_reduce(EPMLI(istra),EPMLI(istra),1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
        call mpi_reduce(EPIOI(istra),EPIOI(istra),1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
        call mpi_reduce(EPPHTI(istra),EPPHTI(istra),1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
        call mpi_reduce(EPPLI(istra),EPPLI(istra),1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
C
C
C        call mpi_reduce(estimv,estimv,nestm1,
C     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
C        do iv=1,nvoltl
C          call mpi_reduce(estimv(iv,1:nrad),estimv(iv,1:nrad),nrad,
C     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
C        end do
        do ir=1,nrtal
          call mpi_reduce(estimv(1,ir),estimv(1,ir),nvoltl,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
        end do
        call mpi_reduce(estims,estims,nestm2,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
        do ispc=1,nadspc
          call mpi_reduce(estiml(ispc)%pspc%spc,estiml(ispc)%pspc%spc,
     .                    estiml(ispc)%pspc%nspc,
     .                    mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
        end do
C
        call mpi_reduce(sdvi1,sdvi1,nsdvi1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
        call mpi_reduce(sdvi2,sdvi2,nsdvi2,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
        call mpi_reduce(sigmac,sigmac,nsdvc1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
        call mpi_reduce(sgmcs,sgmcs,nsdvc2,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)

        call mpi_reduce(LOGMOL(0,ISTRA),LOGMOL(0,ISTRA),NMOLI+1,
     .                  mpi_logical,mpi_LOR,0,icomgrp(istra),ier1)
        call mpi_reduce(LOGATM(0,ISTRA),LOGATM(0,ISTRA),NATMI+1,
     .                  mpi_logical,mpi_LOR,0,icomgrp(istra),ier1)
        call mpi_reduce(LOGION(0,ISTRA),LOGION(0,ISTRA),NIONI+1,
     .                  mpi_logical,mpi_LOR,0,icomgrp(istra),ier1)
        call mpi_reduce(LOGPHOT(0,ISTRA),LOGPHOT(0,ISTRA),NPHOTI+1,
     .                  mpi_logical,mpi_LOR,0,icomgrp(istra),ier1)
        call mpi_reduce(LOGPLS(0,ISTRA),LOGPLS(0,ISTRA),NPLSI+1,
     .                  mpi_logical,mpi_LOR,0,icomgrp(istra),ier1)

        call mpi_barrier(icomgrp(istra),ier)
        call mpi_group_free (igrp(istra),ier)
        call mpi_comm_free (icomgrp(istra),ier)

      endif

      call mpi_barrier(mpi_comm_world,ier)

      RETURN
      END
