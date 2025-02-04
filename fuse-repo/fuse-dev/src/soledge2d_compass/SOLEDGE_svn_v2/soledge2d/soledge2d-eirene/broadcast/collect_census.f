      subroutine EIRENE_collect_census
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMNNL
      USE EIRMOD_COUTAU
      USE EIRMOD_COMSOU
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CPES
 
      IMPLICIT NONE
 
      INCLUDE 'mpif.h'
      real(dp), allocatable :: rpselect(:), rand(:), rdistrib(:),
     .                         rscat(:)
      real(dp) :: ra, weight, peflux, totflux, sumrpw, sclfac, add,
     .            totrpw
      real(dp), external :: ranf_eirene
      integer, allocatable :: iranpro(:), ichose(:)
      integer :: ier, i, istr, ncoreal, itotal, il, im, iu, ipe,
     .           ityp, iphot, iatm, imol, iion, is
      integer :: icopro(0:nprs), idistrib(0:nprs), icosend(0:nprs)

      INTEGER, ALLOCATABLE, SAVE :: IPART_TMP(:,:)
      real(dp), allocatable, save ::  RPART_TMP(:,:)
      INTEGER :: RANK_TEST, IERROR_TEST
 
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
!pb   CALL MPI_ALLREDUCE (FLXFAC,FLXFAC,NSTRAI+1,MPI_REAL8,
!pb  .                    MPI_SUM,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (FLXFAC,NSTRAI+1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
 
      RPARTW(0)=0.0
 
! each processor prepares his census for transfer to processor 0
 
!pb      write (iunout,*) ' my_pe, iprnli ',my_pe, iprnli
!pb      write (iunout,*) ' i,x0,y0,z0,weight '
      PEFLUX=0._DP
      DO I=1,IPRNLI
        ISTR=IPART(8,I)
        ITYP=ISPEZI(IPART(9,I),-1)
        WEIGHT=RPART(9,I)
        IF (ITYP.EQ.0) THEN
          IPHOT=ISPEZI(IPART(I,9),0)
          ADD=WEIGHT*FLXFAC(ISTR)*NPRT(IPHOT)
        ELSEIF (ITYP.EQ.1) THEN
          IATM=ISPEZI(IPART(9,I),1)
          ADD=WEIGHT*FLXFAC(ISTR)*NPRT(NSPH+IATM)
        ELSEIF (ITYP.EQ.2) THEN
          IMOL=ISPEZI(IPART(9,I),2)
          ADD=WEIGHT*FLXFAC(ISTR)*NPRT(NSPA+IMOL)
        ELSEIF (ITYP.EQ.3) THEN
          IION=ISPEZI(IPART(9,I),3)
          ADD=WEIGHT*FLXFAC(ISTR)*NPRT(NSPAM+IION)
        ENDIF
        RPARTW(I)=RPARTW(I-1)+WEIGHT*FLXFAC(ISTR)
        PEFLUX = PEFLUX + ADD
!        write (iunout,'(i6,5es12.4)') i,rpart(1:3,i),weight,
!     .                                  WEIGHT*FLXFAC(ISTR)
      END DO
 
      write (iunout,*) ' my_pe ',my_pe
      write (iunout,*) ' flux on census ',peflux
 
!      write (iunout,*) 'rpartw'
!      write (iunout,'(6es12.4)') rpartw(1:iprnli)
! transfer maximum rpartw to processor 0
 
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
 
      call mpi_allreduce(iprnli,itotal,1,MPI_INTEGER,
     .                   MPI_SUM,MPI_COMM_WORLD,ier)
 
      call mpi_allreduce(peflux,totflux,1,MPI_REAL8,
     .                   MPI_SUM,MPI_COMM_WORLD,ier)
 
      write (iunout,*) ' itotal, totflux', itotal, totflux
      if (itotal <= nprnl) then
 
        call mpi_gather(iprnli,1,MPI_INTEGER,
     .                  icopro,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
 
! send all particles to processor 0
 
        if (my_pe == 0) then
          icosend = icopro*npartt
          idistrib(0) = 0
          do ipe = 1, nprs
            idistrib(ipe) = idistrib(ipe-1) + icosend(ipe-1)
          end do
        end if

        call MPI_COMM_RANK(MPI_COMM_WORLD, RANK_TEST, IERROR_TEST)
!        if (RANK_TEST /= 0 ) then
!        write(6,*) 'SIZE(rpart,rpartc): ', SIZE(rpart),SIZE(rpartc)
!        write(6,*) 'MPI process: ',RANK_TEST,
!     .             ' rpart SEND: ',iprnli*npartt
!        end if
!        if (RANK_TEST == 0 ) then
!        write(6,*) '++SIZE(rpart,rpartc): ', SIZE(rpart),SIZE(rpartc)
!        write(6,*) 'rpart RECEIVE: ', icosend
!        write(6,*) 'rpart DISPL: ', idistrib
!        end if

        IF (.NOT.ALLOCATED(RPART_TMP)) THEN
            ALLOCATE (RPART_TMP(NPARTT,NPRNL))
        END IF
        RPART_TMP = RPART

        CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

        call mpi_gatherv(RPART_TMP,iprnli*npartt,MPI_REAL8,
     .                   rpart,icosend,idistrib,MPI_REAL8,
     .                   0,MPI_COMM_WORLD,ier)

        IF (ALLOCATED(RPART_TMP)) THEN
            DEALLOCATE(RPART_TMP)
        END IF

!        write(6,*) 'MPI process: ', RANK_TEST, 'rpart -- SENT'

        CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

        if (my_pe == 0) then
          icosend = icopro*mpartt
          idistrib(0) = 0
          do ipe = 1, nprs
            idistrib(ipe) = idistrib(ipe-1) + icosend(ipe-1)
          end do
        end if

        call MPI_COMM_RANK(MPI_COMM_WORLD, RANK_TEST, IERROR_TEST)
!        if (RANK_TEST /= 0 ) then
!        write(6,*) 'SIZE(ipart,ipartc): ', SIZE(ipart),SIZE(ipartc)
!        write(6,*) 'MPI process: ',RANK_TEST,
!     .             ' irpart SEND: ',iprnli*mpartt
!        end if
!        if (RANK_TEST == 0 ) then
!        write(6,*) '++SIZE(ipart,ipartc): ', SIZE(ipart),SIZE(ipartc)
!        write(6,*) 'ipart RECEIVE: ', icosend
!        write(6,*) 'ipart DISPL: ', idistrib
!        end if

        IF (.NOT.ALLOCATED(IPART_TMP)) THEN
            ALLOCATE (IPART_TMP(MPARTT,NPRNL))
        END IF
        IPART_TMP = IPART

        CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

        call mpi_gatherv(IPART_TMP,iprnli*mpartt,MPI_INTEGER,
     .                   ipart,icosend,idistrib,MPI_INTEGER,
     .                   0,MPI_COMM_WORLD,ier)

        IF (ALLOCATED(IPART_TMP)) THEN
            DEALLOCATE(IPART_TMP)
        END IF

!        write(6,*) 'MPI process: ', RANK_TEST, 'ipart -- SENT'

        CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
 
        iprnli = itotal
        write (iunout,*) 'total flux on census ', totflux
 
      else
 
        itotal = nprnl
 
        allocate (rpselect(-1:nprs))

        allocate (ichose(nprnl))
        ichose = 0

        if (my_pe == 0) then
          rpselect(-1) = 0._dp
          rpselect(0) = RPARTW(iprnli)
        end if
 
        CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
        call mpi_gather(rpartw(iprnli),1,MPI_REAL8,
     .                rpselect(0:),1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
 
!pb        write (iunout,*) ' rpselect before summation '
!pb        write (iunout,'(i6,es12.4)') (ipe,rpselect(ipe),ipe=-1,nprs)
 
        if (my_pe == 0) then
          do ipe=1, nprs-1
            rpselect(ipe) = rpselect(ipe-1) + rpselect(ipe)
          end do
 
!pb          write (iunout,*) 'total flux on census ', totflux
!pb          write (iunout,*) 'rpselect '
!pb          write (iunout,'(i6,es12.4)') (ipe,rpselect(ipe),ipe=-1,nprs-1)
 
! now find random numbers for NPRNL particles
 
          allocate (rand(nprnl))
          allocate (iranpro(nprnl))
          icopro = 0
 
!pb          write (iunout,*) ' randomly chosen particles '
          do i = 1, nprnl
            ra = ranf_eirene() * rpselect(nprs-1)
 
            IL=0
            IU=NPRS
 
            if ( ra <= rpselect(0) ) then
              iu = 0
            else
c  binary search
              DO WHILE (IU-IL.gt.1)
                IM=(IU+IL)*0.5
                IF (RA.GE.rpselect(IM)) THEN
                  IL=IM
                ELSE
                  IU=IM
                ENDIF
              END DO
            end if
 
            icopro(iu) = icopro(iu) + 1
            iranpro(i) = iu
            rand(i) = ra - rpselect(iu-1)
 
!pb            write (iunout,*) i, ra, iu, rand(i)
          end do
 
! setup displacements for distribution of random numbers
          idistrib(0) = 0
          do ipe = 1, nprs-1
            idistrib(ipe) = idistrib(ipe-1) + icopro(ipe-1)
          end do
 
! sort random numbers
          allocate (rdistrib(nprnl))
          do i = 1, nprnl
            ipe = iranpro(i)
            idistrib(ipe) = idistrib(ipe) + 1
            rdistrib(idistrib(ipe)) = rand(i)
          end do
 
! reset displacements for distribution of random numbers
          idistrib(0) = 0
          do ipe = 1, nprs
            idistrib(ipe) = idistrib(ipe-1) + icopro(ipe-1)
          end do
 
        write (iunout,*) 'no of particles per processor '
        write (iunout,'(10i7)') (icopro(ipe),ipe=0,nprs)
 
        end if
 
! broadcast numbers of needed particles per processor
 
! yannick
       	if (.not.allocated(rdistrib)) allocate(rdistrib(nprnl))

        CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
        call mpi_scatter(icopro ,1,MPI_INTEGER,
     .                   ncoreal,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
 
! broadcast random numbers for each processor
        allocate(rscat(ncoreal))
        call mpi_scatterv(rdistrib,icopro,idistrib,MPI_REAL8,
     .                    rscat,ncoreal,MPI_REAL8,0,MPI_COMM_WORLD,ier)
 
 
! on each processor look for the indices of the particles to be
! put into the global census array
 
!pb        write (iunout,*) ' ncoreal ',ncoreal
!pb        write (iunout,*) ' particles selected '
        sumrpw = 0._dp
        is = 0
        do i = 1, ncoreal
 
          RA = RSCAT(I)
 
          IL=0
          IU=IPRNLI
 
c  binary search
          DO WHILE (IU-IL.gt.1)
            IM=(IU+IL)*0.5
            IF (RA.GE.RPARTW(IM)) THEN
              IL=IM
            ELSE
              IU=IM
            ENDIF
          end do

          if (ichose(iu) == 0) then
! particle IU is hit for the first time
            is = is + 1
            rpartc(:,is) = rpart(:,iu)
            ipartc(:,is) = ipart(:,iu)
            ichose(iu) = is
            write (iunout,*) ' is = ',is,' weight = ',rpartc(9,is)
          else
! particle IU is hit repeatedly, increase weight=rpartc(9,i)
            rpartc(9,ichose(iu)) = rpartc(9,ichose(iu)) + rpart(9,iu)
            write (iunout,*) ' is = ',ichose(iu),
     .                       ' weight = ',rpartc(9,ichose(iu))
          end if
        end do

        ncoreal = is

        do i=1, ncoreal
          ISTR=IPARTC(8,I)
          ITYP=ISPEZI(IPARTC(9,I),-1)
          WEIGHT=RPARTC(9,I)
          IF (ITYP.EQ.0) THEN
             IPHOT=ISPEZI(IPARTC(I,9),0)
             ADD=WEIGHT*FLXFAC(ISTR)*NPRT(IPHOT)
          ELSEIF (ITYP.EQ.1) THEN
             IATM=ISPEZI(IPARTC(9,I),1)
             ADD=WEIGHT*FLXFAC(ISTR)*NPRT(NSPH+IATM)
          ELSEIF (ITYP.EQ.2) THEN
             IMOL=ISPEZI(IPARTC(9,I),2)
             ADD=WEIGHT*FLXFAC(ISTR)*NPRT(NSPA+IMOL)
          ELSEIF (ITYP.EQ.3) THEN
             IION=ISPEZI(IPARTC(9,I),3)
             ADD=WEIGHT*FLXFAC(ISTR)*NPRT(NSPAM+IION)
          ENDIF
          sumrpw = sumrpw + add
!pb          write (iunout,'(i6,4es12.4)') i,rpart(1:3,i),weight
 
!          write (iunout,*) i, ra, iu
!          write (iunout,'(15i7)') ipartc(:,i)
!          write (iunout,'(6es12.4,/,6es12.4)') rpartc(:,i)
 
 
        end do


! send chosen particles to processor 0

        call mpi_gather(ncoreal,1,MPI_INTEGER,
     .                  icopro,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
 
        if (my_pe == 0) then
          icosend = icopro*npartt
          idistrib(0) = 0
          do ipe = 1, nprs
            idistrib(ipe) = idistrib(ipe-1) + icosend(ipe-1)
          end do
        end if

        call MPI_COMM_RANK(MPI_COMM_WORLD, RANK_TEST, IERROR_TEST)
!        if (RANK_TEST /= 0 ) then
!        write(6,*) 'SIZE(rpart,rpartc): ', SIZE(rpart),SIZE(rpartc)
!        write(6,*) 'MPI process: ',RANK_TEST,
!     .             ' rpartc SEND: ',ncoreal*npartt
!        end if
!        if (RANK_TEST == 0 ) then
!        write(6,*) '++SIZE(rpart,rpartc): ', SIZE(rpart),SIZE(rpartc)
!        write(6,*) 'rpartc RECEIVE: ', icosend
!        write(6,*) 'rpartc DISPL: ', idistrib
!        end if

        CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

        call mpi_gatherv(rpartc,ncoreal*npartt,MPI_REAL8,
     .                   rpart,icosend,idistrib,MPI_REAL8,
     .                   0,MPI_COMM_WORLD,ier)

!        write(6,*) 'MPI process: ', RANK_TEST, 'rpartc -- SENT'

        CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
 
        if (my_pe == 0) then
          icosend = icopro*mpartt
          idistrib(0) = 0
          do ipe = 1, nprs
            idistrib(ipe) = idistrib(ipe-1) + icosend(ipe-1)
          end do
        end if

        call MPI_COMM_RANK(MPI_COMM_WORLD, RANK_TEST, IERROR_TEST)
!        if (RANK_TEST /= 0 ) then
!        write(6,*) 'SIZE(ipart,ipartc): ', SIZE(ipart),SIZE(ipartc)
!        write(6,*) 'MPI process: ',RANK_TEST,
!     .             ' irpartc SEND: ',ncoreal*mpartt
!        end if
!        if (RANK_TEST == 0 ) then
!        write(6,*) '++SIZE(ipart,ipartc): ', SIZE(ipart),SIZE(ipartc)
!        write(6,*) 'ipartc RECEIVE: ', icosend
!        write(6,*) 'ipartc DISPL: ', idistrib
!        end if

        CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

        call mpi_gatherv(ipartc,ncoreal*mpartt,MPI_INTEGER,
     .                   ipart,icosend,idistrib,MPI_INTEGER,
     .                   0,MPI_COMM_WORLD,ier)
 
!        write(6,*) 'MPI process: ', RANK_TEST, 'ipartc -- SENT'

        CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

        call mpi_reduce (sumrpw,totrpw,1,MPI_REAL8,MPI_SUM,
     .                   0,MPI_COMM_WORLD,ier)
 
!pb        if (my_pe == 0) then
!pb           write (iunout,*) ' rpart collected from all '
!pb           do i=1, nprnl
!pb             write (iunout,'(i6,4es12.4)') i,rpart(1:3,i),weight
!              write (iunout,*) ' i, my_pe ',i,my_pe
!              write (iunout,'(15i6)') ipart(:,i)
!              write (iunout,'(6es12.4,/,6es12.4)') rpart(:,i)
!pb           end do
!pb        end if
 
        iprnli = sum(icopro)
 
! rescale
 
        if (my_pe == 0) then
          sclfac = totflux / totrpw
          write (iunout,*) ' totrpw ',totrpw
          write (iunout,*) ' sclfac ',sclfac
          do i=1,iprnli
            rpart(9,i) = rpart(9,i) * sclfac
          end do
        end if
 
        if (allocated(rscat)) deallocate(rscat)
        if (allocated(rdistrib)) deallocate(rdistrib)
        if (allocated(rand)) deallocate (rand)
        if (allocated(iranpro)) deallocate (iranpro)
        if (allocated(rpselect)) deallocate (rpselect)
        if (allocated(ichose)) deallocate (ichose)
 
      end if
 
      return
 
      end subroutine EIRENE_collect_census
