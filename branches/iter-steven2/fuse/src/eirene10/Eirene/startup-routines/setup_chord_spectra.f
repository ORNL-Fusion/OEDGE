c  14.5.06:  bug fix: 1 line added: if nchtal.ne.1 and. nchtal.ne.3:  cycle

      subroutine setup_chord_spectra

      use precision
      use parmmod
      use cestim
      use comsig
      use comprt
      use cupd
      use ccona

      implicit none
      
      real(dp) :: c1(3), c2(3), PSIG(0:NSPZ+10)
      real(dp) :: ze, timax
      integer :: ichori, ifirst, ichrd, ipvot, nbc2, nac2, iplots, 
     .           ispc, ntot_cell, ntotsp, iprtyp
      
      type(spect_array), allocatable :: svestiml(:), svsmestl(:)
      TYPE(EIRENE_SPECTRUM), POINTER :: ESPEC, SSPEC
      TYPE(CELL_INFO), POINTER :: FIRST, CUR

!  FIND CELLS INTERSECTED BY CHORDS

      ze = 1._dp

      ifirst = -1
      ntot_cell = 0

      do ichori = 1, nchori

        IF (.NOT.NLSTCHR(ICHORI)) CYCLE
        IF ((NCHTAL(ICHORI) /= 1) .AND. (NCHTAL(ICHORI) /= 3) .AND.
     .      (NCHTAL(ICHORI) /= 4) ) CYCLE

        IPVOT=IPIVOT(ICHORI)
        C1(1)=XPIVOT(ICHORI)
        C1(2)=YPIVOT(ICHORI)
        C1(3)=ZPIVOT(ICHORI)
C     
        ICHRD=ICHORD(ICHORI)
        C2(1)=XCHORD(ICHORI)
        C2(2)=YCHORD(ICHORI)
        C2(3)=ZCHORD(ICHORI)
C     
        NBC2=NSPBLC(ICHORI)
        NAC2=NSPADD(ICHORI)

        ALLOCATE(TRAJ(ICHORI)%TRJ)
        TRAJ(ICHORI)%TRJ%P1 = C1
        TRAJ(ICHORI)%TRJ%P2 = C2
        TRAJ(ICHORI)%TRJ%NCOU_CELL = 0
        NULLIFY(TRAJ(ICHORI)%TRJ%CELLS)
        
        CALL LININT (IFIRST,ICHORI,C1,C2,ICHRD,IPVOT,NBC2,NAC2,ZE,
     .               PSIG,TIMAX,1,1,1,IABS(NCHENI))

        ntot_cell = ntot_cell + traj(ichori)%trj%ncou_cell
         
      end do

      IF (NTOT_CELL == 0) RETURN

!  SAVE SPECTRA SPECIFIED VIA INPUT

      IF (NADSPC > 0) THEN

        ALLOCATE(SVESTIML(NADSPC))

        DO ISPC = 1, NADSPC
          SVESTIML(ISPC)%PSPC => ESTIML(ISPC)%PSPC
        END DO
        
        DEALLOCATE(ESTIML)

        IF (ALLOCATED(SMESTL)) THEN
          ALLOCATE(SVSMESTL(NADSPC))
          DO ISPC = 1, NADSPC
            SVSMESTL(ISPC)%PSPC => SMESTL(ISPC)%PSPC
          END DO
          DEALLOCATE(SMESTL)
        END IF        

      END IF

!  set up new arrays for spectra

      NTOTSP = NADSPC + NTOT_CELL

      ALLOCATE(ESTIML(NTOTSP))
      DO ISPC = 1, NADSPC
        ESTIML(ISPC)%PSPC => SVESTIML(ISPC)%PSPC
      END DO

      IF (ALLOCATED(SVSMESTL).or.NSMSTRA.GT.0) THEN
        ALLOCATE(SMESTL(NTOTSP))
        DO ISPC = 1, NADSPC
          SMESTL(ISPC)%PSPC => SVSMESTL(ISPC)%PSPC
        END DO
      END IF
      
!  add spectra for cells along chords

      ispc = nadspc
      do ichori = 1,nchori

        IF ((NCHTAL(ICHORI) /= 1) .AND. (NCHTAL(ICHORI) /= 3) .AND.
     .      (NCHTAL(ICHORI) /= 4) ) CYCLE

         if (.not.associated(traj(ichori)%trj%cells)) cycle
         first => traj(ichori)%trj%cells
         cur => first

         if (nchtal(ichori) == 1) then
           iprtyp = 1
         else if (nchtal(ichori) == 3) then
           iprtyp = 0
         else if (nchtal(ichori) == 4) then
           iprtyp = 1
         else
           write (iunout,*) ' wrong chord type for spectrum '
           write (iunout,*) ' no spectrum set up for cell',cur%no_cell
           cycle
         end if

!  loop over all cells along trajectory
         do 
           allocate(espec)
           espec%isrfcll = 2
           espec%ispcsrf = cur%no_cell
           espec%iprtyp = iprtyp
           espec%iprsp = nspspz(ichori)
           espec%ispctyp = 1
           espec%nspc = abs(ncheni)
           espec%imetsp = 0
           espec%idirec = 1
           if (ncheni > 0) then
             espec%spcmin = emin1(ichori)
             espec%spcmax = emax1(ichori)
             espec%log = .false.
           else
             espec%spcmin = log10(emin1(ichori))
             espec%spcmax = log10(emax1(ichori))
             espec%log = .true.
           end if
           espec%esp_00 = 0._dp
           espec%spc_xplt = 0._dp
           espec%spc_yplt = 0._dp
           espec%spc_same = 0._dp
           espec%spcvx = traj(ichori)%trj%vx
           espec%spcvy = traj(ichori)%trj%vy
           espec%spcvz = traj(ichori)%trj%vz
           espec%esp_min =1.e30_dp
           espec%esp_max = -1.e30_dp
           espec%spcdel=(espec%spcmax-espec%spcmin)/real(espec%nspc,dp)
           espec%spcdeli = 1._dp / (espec%spcdel+eps60)
           allocate(espec%spc(0:espec%nspc+1))
           allocate(espec%sdv(0:espec%nspc+1))
           allocate(espec%sgm(0:espec%nspc+1))
           espec%spc(0:espec%nspc+1) = 0
           
           ispc = ispc + 1
           estiml(ispc)%pspc => espec

           if (allocated(smestl)) then
             allocate(sspec)
             allocate(sspec%spc(0:espec%nspc+1))
             allocate(sspec%sdv(0:espec%nspc+1))
             allocate(sspec%sgm(0:espec%nspc+1))
             sspec = espec
             smestl(ispc)%pspc => sspec
           end if

           cur%no_spect = ispc
           cur => cur%nextc
!pb associated with two arguments tests if both arguments point to the same target
           if (associated(cur,first)) exit
         end do
      end do

      NADSPC = ISPC

      end subroutine setup_chord_spectra
