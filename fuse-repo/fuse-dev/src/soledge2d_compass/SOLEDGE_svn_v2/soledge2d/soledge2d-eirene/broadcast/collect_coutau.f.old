      subroutine EIRENE_collect_coutau
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CPES
      USE EIRMOD_COUTAU
      USE EIRMOD_CSPEZ
      use EIRMOD_CESTIM
      use EIRMOD_COMUSR
      use EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CGRID
      USE EIRMOD_COMSOU
      USE EIRMOD_CSDVI
      USE EIRMOD_CSPEI
      IMPLICIT NONE
 
      include 'mpif.h'
      REAL(DP), ALLOCATABLE :: OUTAU(:)
      integer :: ier, icolor, istr, icomgrp, ier1, isdv, i
 
      if (nsteff < nprs) then
! collect from group leader pe
        icolor=MPI_UNDEFINED
        do istr=1,nstrai
          if (npts(istr) == 0) cycle
          if (my_pe == npesta(istr)) then
            icolor = 1
            write (0,*) ' my_pe, istr ',my_pe, istr
            exit
          end if
        end do
      else
! collect from all pes
        icolor = 1
      end if
 
      call mpi_barrier(mpi_comm_world,ier)
 
      call mpi_comm_split (mpi_comm_world,icolor,my_pe,icomgrp,ier)
 
      if (icolor == 1) then
        ALLOCATE (OUTAU(NOUTAU))
        CALL EIRENE_WRITE_COUTAU (OUTAU, IUNOUT)
 
        CALL MPI_REDUCE(OUTAU,OUTAU,NOUTAU,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier)
 
        if (my_pe == 0) CALL EIRENE_READ_COUTAU (OUTAU, IUNOUT)
        DEALLOCATE (OUTAU)
 
        call mpi_reduce(LOGMOL,LOGMOL,(NMOLI+1)*(NSTRA+1),
     .                  mpi_logical,mpi_LOR,0,icomgrp,ier)
        call mpi_reduce(LOGATM,LOGATM,(NATMI+1)*(NSTRA+1),
     .                  mpi_logical,mpi_LOR,0,icomgrp,ier)
        call mpi_reduce(LOGION,LOGION,(NIONI+1)*(NSTRA+1),
     .                  mpi_logical,mpi_LOR,0,icomgrp,ier)
        call mpi_reduce(LOGPHOT,LOGPHOT,(NPHOTI+1)*(NSTRA+1),
     .                  mpi_logical,mpi_LOR,0,icomgrp,ier)
        call mpi_reduce(LOGPLS,LOGPLS,(NPLSI+1)*(NSTRA+1),
     .                  mpi_logical,mpi_LOR,0,icomgrp,ier)
 
        if (nsmstra > 0) then
 
        CALL MPI_REDUCE(SMESTV,SMESTV,NIDV*NRTAL,
     .                mpi_real8,mpi_sum,0,icomgrp,ier1)
        CALL MPI_REDUCE(SMESTS,SMESTS,NIDS*NLMPGS,
     .                  mpi_real8,mpi_sum,0,icomgrp,ier1)
 
        DO I=1,NADSPC
          CALL MPI_REDUCE(SMESTL(I)%PSPC%SPC,SMESTL(I)%PSPC%SPC,
     .                    SMESTL(I)%PSPC%NSPC+2,
     .                    MPI_REAL8,MPI_SUM,0,icomgrp,IER1)
          CALL MPI_REDUCE(SMESTL(I)%PSPC%SPCINT,SMESTL(I)%PSPC%SPCINT,
     .                    1,MPI_REAL8,MPI_SUM,0,icomgrp,IER1)
          if (nsigi_spc > 0) then
            call mpi_reduce(smestl(i)%pspc%sdv,smestl(i)%pspc%sdv,
     .                      smestl(i)%pspc%nspc+2,
     .                      mpi_real8,mpi_sum,0,icomgrp,ier1)
            call mpi_reduce(smestl(i)%pspc%sgm,smestl(i)%pspc%sgm,
     .                      smestl(i)%pspc%nspc+2,
     .                      mpi_real8,mpi_sum,0,icomgrp,ier1)
            call mpi_reduce(smestl(i)%pspc%stvs,
     .                      smestl(i)%pspc%stvs,1,
     .                      mpi_real8,mpi_sum,0,icomgrp,ier1)
            call mpi_reduce(smestl(i)%pspc%ees,
     .                      smestl(i)%pspc%ees,1,
     .                      mpi_real8,mpi_sum,0,icomgrp,ier1)
          end if
        END DO
 
        end if
 
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
          CALL MPI_REDUCE(EES,EES,NSIGVI,
     .                    MPI_REAL8,MPI_SUM,0,ICOMGRP,IER1)
        END IF
 
 
        IF (NSIGSI > 0) THEN
          CALL MPI_REDUCE(STVW(1:NSIGSI,1:NLIMPS),
     .                    STVW(1:NSIGSI,1:NLIMPS),NSIGSI*NLIMPS,
     .                    MPI_REAL8,MPI_SUM,0,ICOMGRP,IER1)
          CALL MPI_REDUCE(FF(1:NSIGSI,1:NLIMPS),
     .                    FF(1:NSIGSI,1:NLIMPS),NSIGSI*NLIMPS,
     .                    MPI_REAL8,MPI_SUM,0,ICOMGRP,IER1)
          CALL MPI_REDUCE(STVWS,STVWS,NSIGSI,
     .                    MPI_REAL8,MPI_SUM,0,ICOMGRP,IER1)
          CALL MPI_REDUCE(FFS,FFS,NSIGSI,
     .                    MPI_REAL8,MPI_SUM,0,ICOMGRP,IER1)
        END IF
 
      end if ! icolor=1
 
      call mpi_barrier(mpi_comm_world,ier)
 
      return
      end subroutine EIRENE_collect_coutau
