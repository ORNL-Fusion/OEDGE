*DK INFCOP
* data exchange for styx2D-Eirene
* currently assumes pure plasma (one charged ion species) 
* started April 15, 2011 Y. Marandet
* yannick.marandet@univ-amu.fr
*
      SUBROUTINE EIRENE_INFCOP
      use eirmod_precision
      use eirmod_parmmod
      use eirmod_cspei
      use eirmod_comprt, only : iunout,nltrc
      use styx2eirene
      use eirmod_comsou
      use eirmod_comusr
      use eirmod_cgrid
      use eirmod_ccona
      use eirmod_coutau
      use eirmod_cestim
      use eirmod_comnnl
      use eirmod_ctrig
      use eirmod_cplot
      use eirmod_ctrcei
      use eirmod_clgin
      use eirmod_comxs
      use eirmod_csdvi
      use eirmod_cpes
      use eirmod_ctext
      use all_variables, only : global_parameters

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I1, I2, I3
      integer :: i,ITRI,ISTRA,ITAL,ipls,spz_found
      integer :: iatm,imol
      
      real(dp) :: RP1,THMAX
      real(dp) :: xG,yG,xm,ym,t
      
      real(dp) :: vpar_itri,Temin
      logical :: check_exchg,dir_e,exit_on_nan

      real(dp) :: tst1,tend1,omp_get_wtime
      character(8) :: txtbuff1
      character(2) :: txtbuff2

      

      ENTRY EIRENE_IF0COP
      return

      entry eirene_if1cop

         ! magnetic field
        BXIN(1:Neir_cells)  = BX_tri
        BYIN(1:Neir_cells)  = BY_tri
        BZIN(1:Neir_cells)  = BZ_tri
        BFIN(1:Neir_cells)  = BF_tri

! just to avoid exiting on B=0 condition
        !BXIN(Ntri_styx+1)=1._dp
        !BFIN(Ntri_styx+1)=1._dp


! THIS SHOULD NOT BE NEEDED ANYMORE : CHECK !!
! set up wall material according to eirene_coupling.txt
! atomic mass - nuclear charge (e.g. 12 - 06 for C)

!        select case (wall_mat)
!          case("Be")
!            ZNML(0:NLIMPS)=09._dp
!            ZNCL(0:NLIMPS)=04._dp
!          case("B")
!            ZNML(0:NLIMPS)=11._dp
!            ZNCL(0:NLIMPS)=05._dp
!          case("C")
!            ZNML(0:NLIMPS)=12._dp
!            ZNCL(0:NLIMPS)=06._dp
!          case("Fe")
!            ZNML(0:NLIMPS)=56._dp
!            ZNCL(0:NLIMPS)=26._dp
!          case("Mo")
!            ZNML(0:NLIMPS)=96._dp
!            ZNCL(0:NLIMPS)=42._dp
!          case("W")
!            ZNML(0:NLIMPS)=184._dp
!            ZNCL(0:NLIMPS)=74._dp
!          case default
!            write(*,*) ' Problem in wall material specification 
!     .     (from if2cop), wall_mat = ',wall_mat
!            write(*,*) ' Exit ... '
!            call eirene_exit_own(1)
!        end select      


! highjack the value of NSTEP
! as many step functions as recycling strata * toroidal segments

        NSTEP=Nrecyc*Ntor_cells


! switch off plotting
                
        NLPL2D=.false.
        NLPL3D=.false.
        PLHST=.false.
        TRCHST=.false.

! initialize timer (my_pe = 0)

        if (.not.allocated(time_para)) allocate(time_para(nprs))
        time_para=0

      return

      ENTRY EIRENE_IF1COPls

c make sure there are no nans in plasma background (otherwise strange crashes in EIRENE)

        exit_on_nan=.false.

        if (sum(eiv_e%T) .ne. sum(eiv_e%T)) then
          write(*,*) ' NaNs in Te at Eirene entrance,stop'
          exit_on_nan=.true.
        endif
  
        ! multi-temperature case to be handled properly ...
        if (sum(eiv_i(1)%T) .ne. sum(eiv_i(1)%T)) then
          write(*,*) ' NaNs in Ti at Eirene entrance,stop'
          exit_on_nan=.true.
        endif      

        do ipls=1,nplsi
          if (sum(eiv_i(ipls)%dens) .ne. sum(eiv_i(ipls)%dens)) then
            write(*,*) ' NaNs in Te at Eirene entrance,stop'
            exit_on_nan=.true.
          endif

          if (sum(eiv_i(ipls)%vx) .ne. sum(eiv_i(ipls)%vx)) then
            write(*,*) ' NaNs in vx at Eirene entrance,stop'
            exit_on_nan=.true.
          endif

          if (sum(eiv_i(ipls)%vy) .ne. sum(eiv_i(ipls)%vy)) then
            write(*,*) ' NaNs in vy at Eirene entrance,stop'
            exit_on_nan=.true.
          endif

          if (sum(eiv_i(ipls)%vz) .ne. sum(eiv_i(ipls)%vz)) then
            write(*,*) ' NaNs in vz at Eirene entrance,stop'
            exit_on_nan=.true.
          endif
        enddo

        if (exit_on_nan) call eirene_exit_own(1)

c make sure calculation of rate coefficients does not blow up

c        if (.not.short_cycle) then
c          write(*,*) ' min Te, eirene entrance (eV) = ',
c     .                minval(eiv_e%T)
c          do ipls=1,nplsi
c            write(*,*) ' min Ti, species #',ipls,
c     .               ' eirene entrance (eV) = ',
c     .                minval(eiv_i(ipls)%T)
c          enddo
c        endif

        do itri=1,Neir_cells
           eiv_e%T(itri)=max(eiv_e%T(itri),Temin_eirene)
        enddo 

        do ipls=1,npls
          do itri=1,Neir_cells
            eiv_i(ipls)%T(itri)=max(eiv_i(ipls)%T(itri),Temin_eirene)
          enddo
        enddo


        TEIN(1:Neir_cells)  = eiv_e%T
 
        do ipls=1,nplsi
          DIIN(ipls,1:Neir_cells) = eiv_i(ipls)%dens*1e-6_dp
          VXIN(ipls,1:Neir_cells)= eiv_i(ipls)%vx*1e2_dp
          VYIN(ipls,1:Neir_cells)= eiv_i(ipls)%vy*1e2_dp
          VZIN(ipls,1:Neir_cells)= eiv_i(ipls)%vz*1e2_dp
          TIIN(ipls,1:Neir_cells)= eiv_i(ipls)%T
        enddo
        
      return

      ENTRY EIRENE_IF2COP(ISTRA)

      if (ISTRA <= Nrecyc) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now data for recycling strata
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        NLAVRP(ISTRA)=.false. 
        NLAVRT(ISTRA)=.false.
        NLSYMP(ISTRA)=.false.
        NLSYMT(ISTRA)=.false.
       
        NLRAY(ISTRA)=.false.

        NPTS(ISTRA)=Npart_eirene(1)
        NINITL(ISTRA)=seed_eirene(1)+n_call_eir*seed_gap

      ! source angular distribution
        NAMODS(ISTRA)=1

        if (source_species /= 4) then
      !  monokinetic particles
          NEMODS(ISTRA)=1
        else
          if (sheath_model == 0) then
      ! 7 = with sheath, 6 without sheath
            NEMODS(ISTRA)=7
          else
            NEMODS(ISTRA)=0
          endif
         endif

         FLUX(ISTRA)=1._dp
         SCALV(ISTRA)=0._dp
         IVLSF(ISTRA)=0._dp
         ISCLS(ISTRA)=0._dp
         ISCLT(ISTRA)=0._dp
         ISCL1(ISTRA)=0._dp
         ISCL2(ISTRA)=0._dp
         ISCL3(ISTRA)=0._dp
         ISCLB(ISTRA)=0._dp
         ISCLA(ISTRA)=0._dp


         if (source_species == 1 .and. natmi > 0) then
           NLATM(ISTRA)=.true.
           NLMOL(ISTRA)=.false.
           NLION(ISTRA)=.false.
           NLPLS(ISTRA)=.false.
         elseif (source_species == 2 .and. nmoli > 0) then
           NLATM(ISTRA)=.false.
           NLMOL(ISTRA)=.true.
           NLION(ISTRA)=.false.
           NLPLS(ISTRA)=.false.
         elseif (source_species == 3 .and. nioni > 0) then
           NLATM(ISTRA)=.false.
           NLMOL(ISTRA)=.false.
           NLION(ISTRA)=.true.
           NLPLS(ISTRA)=.false.
        elseif (source_species == 4) then
       ! in this case make sure that the ion is directed towards the wall !!!
          NLATM(ISTRA)=.false.
          NLMOL(ISTRA)=.false.
          NLION(ISTRA)=.false.
          NLPLS(ISTRA)=.true.
        else
          write(iunout,*) 'Error is source species specification, 
     .       infcop'
          call eirene_exit_own(1)
        endif 
 
!!! species index sampled according to fluxes , see WEISPZ setup in samsrf
        NSPEZ(ISTRA)=0

        NLPNT(ISTRA)=.false.
        NLLNE(ISTRA)=.false.
        NLSRF(ISTRA)=.true.
        NLVOL(ISTRA)=.false.
        NLCNS(ISTRA)=.false.

        !!!! test to see whether was not optimal !!!
        !!NSRFSI(ISTRA) = Nsou
        NSRFSI(ISTRA) = Ntor_cells
      
        do i=1,NSRFSI(ISTRA)

          INDIM(i,ISTRA)=4
          INSOR(i,ISTRA)=-1
          INGRDA(i,ISTRA,1)=0
          INGRDE(i,ISTRA,1)=0
          INGRDA(i,ISTRA,2)=0
          INGRDE(i,ISTRA,2)=0
          INGRDE(i,ISTRA,3)=0
          INGRDE(i,ISTRA,3)=0

          SORWGT(i,ISTRA)=1._dp
          ! uniform distribution in Z and step functions for "radial" grid
          SORLIM(i,ISTRA)=204._dp
          ! index istep of step function for this stratum
          SORIND(i,ISTRA)=i+(istra-1)*Ntor_cells
          SOREXP(i,ISTRA)=0._dp
    
          SORIFL(i,ISTRA)=0._dp
        
          NRSOR(i,ISTRA)=-1
          NPSOR(i,ISTRA)=0
          NTSOR(i,ISTRA)=0
          NBSOR(i,ISTRA)=0
          NASOR(i,ISTRA)=0
          NISOR(i,ISTRA)=0
          ISTOR(i,ISTRA)=0

          SORAD1(i,ISTRA)=0._dp
          SORAD2(i,ISTRA)=0._dp
          SORAD3(i,ISTRA)=0._dp
          SORAD4(i,ISTRA)=0._dp
! data for uniform distribution in Z direction
          if (is_3D) then
            SORAD5(i,ISTRA)=Zsurf(i)*RADDEG
            SORAD6(i,ISTRA)=Zsurf(i+1)*RADDEG
          else
            SORAD5(i,ISTRA)=0._dp
            SORAD6(i,ISTRA)=100._dp
          endif

        enddo

        SORENI(ISTRA)=0._dp
        SORENE(ISTRA)=0._dp
        SORVDX(ISTRA)=0._dp
        SORVDY(ISTRA)=0._dp
        SORVDZ(ISTRA)=0._dp

        SORCOS(ISTRA)=1
        SORMAX(ISTRA)=0._dp
        SORCTX(ISTRA)=0._dp
        SORCTY(ISTRA)=0._dp
        SORCTZ(ISTRA)=0._dp
        RAYFRAC(ISTRA)=0._dp

        IF (NLSRF(ISTRA))
     .  THMAX=MAX(0._DP,MIN(PIHA,SORMAX(ISTRA)*DEGRAD))
        IF (NAMODS(ISTRA).EQ.1) THEN
          RP1=SORCOS(ISTRA)+1.
          SORCOS(ISTRA)=1./RP1
!pb          IF (ABS(COS(THMAX)).LE.EPS10) THEN
          IF (ABS(COS(THMAX)).LE.EPS5) THEN
            SORMAX(ISTRA)=1.
          ELSE
            SORMAX(ISTRA)=1.-COS(THMAX)**RP1
          ENDIF
       ELSEIF (NAMODS(ISTRA).EQ.2) THEN
          SORCOS(ISTRA)=SORCOS(ISTRA)*DEGRAD
          SORMAX(ISTRA)=THMAX
       ELSE
          WRITE (iunout,*) 'INPUT ERROR: ISTRA,NAMODS(ISTRA)='
          WRITE (iunout,*) '             ',ISTRA,NAMODS(ISTRA)
          CALL EIRENE_EXIT_OWN(1)
       ENDIF


      elseif(ISTRA > Nrecyc .and. ISTRA<=Nrecyc+Nrecomb) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! recombination strata
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        NLAVRP(ISTRA)=.false.
        NLAVRT(ISTRA)=.false.
        NLSYMP(ISTRA)=.false.
        NLSYMT(ISTRA)=.false.
        NLRAY(ISTRA)=.false.
        SCALV(ISTRA)=0._dp
        IVLSF(ISTRA)=0._dp
        ISCLS(ISTRA)=0._dp
        ISCLT(ISTRA)=0._dp
        ISCL1(ISTRA)=0._dp
        ISCL2(ISTRA)=0._dp
        ISCL3(ISTRA)=0._dp
        ISCLB(ISTRA)=0._dp
        ISCLA(ISTRA)=0._dp

        NLPNT(ISTRA)=.false.
        NLLNE(ISTRA)=.false.
        NLSRF(ISTRA)=.false.
        NLVOL(ISTRA)=.true.
        NLCNS(ISTRA)=.false.

! calculated later, but must be different from zero otherwise strata deactivated
        FLUX(ISTRA)=1._dp

        NPTS(ISTRA)=Npart_eirene(1)
        NINITL(ISTRA)=seed_eirene(1)+n_call_eir*seed_gap

        NLATM(ISTRA)=.false.
        NLMOL(ISTRA)=.false.
        NLION(ISTRA)=.false.
        NLPLS(ISTRA)=.true.

!species index: all ions of charge 1+, thus one per species
        NSPEZ(ISTRA)=singly_charged_ions(ISTRA-Nrecyc)
        
        NSRFSI(ISTRA)=1
        SORWGT(1,ISTRA)=1._dp
  
        INGRDA(1,ISTRA,1)=0
        INGRDE(1,ISTRA,1)=0
! SORLIM > 0 : no user defined model here
        SORLIM(1,ISTRA)=1
! this flag identifies the number IRRC of the recombination channel for species IPLS=1
! 0= sum of processes (e.g. radiative + three body)
! if no recombination reaction defined for bulk ion ipls, use Gordev formula ...
        SORIND(1,ISTRA)=1
! velocity space: use the "old" default case in EIRENE : NEMOD1.neq.1
! first sample bulk ion with shifted maxwelian velocity at Ti
! then use this energy for the ion ...
        NEMODS(ISTRA)=0
        

      elseif ((ISTRA > Nrecyc+Nrecomb).or.
     .         ((ISTRA > Nrecyc+Nrecomb).and.(ISTRA < NSTRATA+1))) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now data for gas puffs
! -number, position(s), flux strength(s), gas temperature(s), divergence(s) 
! The puff rate is for D/s, but molecules are puffed (D2 influx = total_puff/2)
! implemented as point source in the middle of a triangle side
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

        NLAVRP(ISTRA)=.false. 
        NLAVRT(ISTRA)=.false.
        NLSYMP(ISTRA)=.false.
        NLSYMT(ISTRA)=.false.
        NLRAY(ISTRA)=.false.
        SCALV(ISTRA)=0._dp
        IVLSF(ISTRA)=0._dp
        ISCLS(ISTRA)=0._dp
        ISCLT(ISTRA)=0._dp
        ISCL1(ISTRA)=0._dp
        ISCL2(ISTRA)=0._dp
        ISCL3(ISTRA)=0._dp
        ISCLB(ISTRA)=0._dp
        ISCLA(ISTRA)=0._dp

        NLPNT(ISTRA)=.true.
        NLLNE(ISTRA)=.false.
        NLSRF(ISTRA)=.false.
        NLVOL(ISTRA)=.false.
        NLCNS(ISTRA)=.false.

      ! source angular distribution = Gaussian
        NAMODS(ISTRA)=2
      ! gaussian source with T given by SORENI and mean velocity SORVDX, ... (=0)
        NEMODS(ISTRA)=116

        NSRFSI(ISTRA)=1
        SORWGT(1,ISTRA)=1._dp

        SORVDX(ISTRA)=0._dp
        SORVDY(ISTRA)=0._dp
        SORVDZ(ISTRA)=0._dp

        SORENE(ISTRA)=0._dp

! scale the number of particles with flux ?
        NPTS(ISTRA)=Npart_eirene(1)
        NINITL(ISTRA)=seed_eirene(1)+n_call_eir*seed_gap
! puff rate in Amp 
        FLUX(ISTRA)=puffs(ISTRA)%rate*elcha

! now determine species puff
        spz_found =0
        txtbuff1=adjustl(puffs(ISTRA)%species)
        call eirene_lowercase(txtbuff1)       
 
        NLATM(ISTRA)=.false.
        NLMOL(ISTRA)=.false.
        NLION(ISTRA)=.false.
        NLPLS(ISTRA)=.false.
 
        ! atoms ?
        do iatm=1,natmi        
          txtbuff2=adjustl(TEXTS(nsph+iatm))
          call eirene_lowercase(txtbuff2)
          if (txtbuff2(1:2)==txtbuff1(1:2)) then
            spz_found = 1
            NLATM(ISTRA)=.true.
            NSPEZ(ISTRA)=iatm
          endif
        enddo
       
        ! molecules ?
        if (spz_found == 0) then
          do imol=1,nmoli
            txtbuff2=adjustl(TEXTS(nsph+natmi+imol))
            call eirene_lowercase(txtbuff2)
            if (txtbuff2(1:2)==txtbuff1(1:2)) then
              spz_found = 1
              NLMOL(ISTRA)=.true.
              NSPEZ(ISTRA)=imol
            endif
          enddo
        endif

        if (spz_found == 0) then
          write(*,*) ' Check species name for puff # ',
     .                    ISTRA-(Nrecyc+Nrecomb)
          write(*,*) ' Species = ',puffs(istra)%species
          stop
        endif


       SORENI(ISTRA)=puffs(ISTRA)%T0

! determine location of point source; make sure it is inside the triangle !
! source located on the segment MG, where M is the middle of the triangle side (defines the direction)
! and G the center of mass of the triangle. Exact location specified by varaiable t (t=0 at M, t=1 at G)

        t=0.1_dp

        itri=puffs(ISTRA)%itri
        xG=1._dp/3._dp*(xtrian(NECKE(1,itri))
     .            +xtrian(NECKE(2,itri))+xtrian(NECKE(3,itri)))
        yG=1._dp/3._dp*(ytrian(NECKE(1,itri))
     .            +ytrian(NECKE(2,itri))+ytrian(NECKE(3,itri)))

! get triangles vertices needed to calculate source position (SORAD1,...,SORAD3)
! side 1 : vertex 1,2
! side 2 : vertex 2,3
! side 3 : vertex 3,1

      if (puffs(ISTRA)%iside == 1) then  
          xm=0.5_dp*(xtrian(NECKE(1,itri))
     .                  +xtrian(NECKE(2,itri)))
          ym=0.5_dp*(ytrian(NECKE(1,itri))
     .                  +ytrian(NECKE(2,itri)))   
      elseif (puffs(ISTRA)%iside == 2) then
          xm=0.5_dp*(xtrian(NECKE(2,itri))
     .                  +xtrian(NECKE(3,itri)))
          ym=0.5_dp*(ytrian(NECKE(2,itri))
     .                  +ytrian(NECKE(3,itri)))
      else
          xm=0.5_dp*(xtrian(NECKE(3,itri))
     .                  +xtrian(NECKE(1,itri)))
          ym=0.5_dp*(ytrian(NECKE(3,itri))
     .                  +ytrian(NECKE(1,itri)))
      endif

      SORAD1(1,ISTRA)=xm*(1._dp-t)+xG*t
      SORAD2(1,ISTRA)=ym*(1._dp-t)+yG*t

! toroidal position (one puff per strat)
      if (is_3D) then
       ! NTSOR(1,istra) = puffs(ISTRA)%itor
        SORAD3(1,istra) = zzone(puffs(ISTRA)%itor)*RADDEG ! toroidal angle in degrees
      else
        SORAD3(1,ISTRA)=0.5_dp*(ZAA+ZZA)
      endif
! then need to set the normal (outer normal calculated in grid.f, here inner normal, but correct one used for sampling ... )
! normal in plane X,Y, no toroidal component

  	SORAD4(1,ISTRA)= PTRIX(puffs(ISTRA)%iside,puffs(ISTRA)%itri)
  	SORAD5(1,ISTRA)= PTRIY(puffs(ISTRA)%iside,puffs(ISTRA)%iside)
  	SORAD6(1,ISTRA)=0._dp

! beam divergence, beam_div is the standard deviation in degrees

  	SORCOS(ISTRA)=puffs(ISTRA)%divergence*DEGRAD
  	SORMAX(ISTRA)=PIA

  	RAYFRAC(ISTRA)=0._dp

      endif


      return

!999   write(iunout,*) 'Error when reading flux file 
!     .   fort.31, from if2cop'
!      call eirene_exit_own(1)

      ENTRY EIRENE_IF3COP(I1,I2,I3)
! this routine is called on a per srata basis
! mpi reduction for each strata already carried out

 ! this if should be uncommented
      if (sc_level > 1) then

! nsetff is the number of strata with non zero computational time
! the case nsteff>nprs corresponds to one process/stratum
! otherwise the group leader has the whole group information

!       if ((nsteff.ge.nprs).or.(npesta(i1).eq.my_pe)) then

! conversion from cm-3 -> m-3
        eov(i1)%pdena = pdena(:,1:Neir_cells)*1e6_dp
        eov(i1)%pdenm = pdenm(:,1:Neir_cells)*1e6_dp
        eov(i1)%pdeni = pdeni(:,1:Neir_cells)*1e6_dp

! conversion from eV/cm3 to J/m3
        eov(i1)%edena = edena(:,1:Neir_cells)*1e6_dp*elcha
        eov(i1)%edenm = edenm(:,1:Neir_cells)*1e6_dp*elcha
        eov(i1)%edeni = edeni(:,1:Neir_cells)*1e6_dp*elcha

! conversion from g.cm.s-1/cm3 to kg.m.s-1/m3
        eov(i1)%vxdena = vxdena(:,1:Neir_cells)*10._dp
        eov(i1)%vxdenm = vxdenm(:,1:Neir_cells)*10._dp
        eov(i1)%vxdeni = vxdeni(:,1:Neir_cells)*10._dp

        eov(i1)%vydena = vydena(:,1:Neir_cells)*10._dp
        eov(i1)%vydenm = vydenm(:,1:Neir_cells)*10._dp
        eov(i1)%vydeni = vydeni(:,1:Neir_cells)*10._dp

        eov(i1)%vzdena = vzdena(:,1:Neir_cells)*10._dp
        eov(i1)%vzdenm = vzdenm(:,1:Neir_cells)*10._dp
        eov(i1)%vzdeni = vzdeni(:,1:Neir_cells)*10._dp

!! these are the direct coupling data

! conversion from Amp/cm3 to part/m3/s (and Amp to part/s)
        eov(i1)%pael = pael(1:Neir_cells)*1e6_dp/elcha
        eov(i1)%pmel = pmel(1:Neir_cells)*1e6_dp/elcha
        eov(i1)%piel = piel(1:Neir_cells)*1e6_dp/elcha

! bulk ion sources, conversion from Amp/cm3 to part/m3/s (and Amp to part/s)
        
        eov(i1)%papl = papl(:,1:Neir_cells)*1e6_dp/elcha
        eov(i1)%pmpl = pmpl(:,1:Neir_cells)*1e6_dp/elcha
        eov(i1)%pipl = pipl(:,1:Neir_cells)*1e6_dp/elcha

! conversion from Watt/cm3 to Watt/m3
        eov(i1)%eael = eael(1:Neir_cells)*1e6_dp
        eov(i1)%emel = emel(1:Neir_cells)*1e6_dp
        eov(i1)%eiel = eiel(1:Neir_cells)*1e6_dp

        eov(i1)%eapl = eapl(1:Neir_cells)*1e6_dp
        eov(i1)%empl = empl(1:Neir_cells)*1e6_dp
        eov(i1)%eipl = eipl(1:Neir_cells)*1e6_dp

! conversion from Amp.g.cm.s-1/cm3 to kg/m.s-1/m3/s
        eov(i1)%mapl = mapl(:,1:Neir_cells)*10._dp/elcha
        eov(i1)%mmpl = mmpl(:,1:Neir_cells)*10._dp/elcha
        eov(i1)%mipl = mipl(:,1:Neir_cells)*10._dp/elcha

      endif

      return
      ENTRY EIRENE_IF4COP
! this is now the sum over strata

! first volume tallies

      eov(0)%pdena = pdena(:,1:Neir_cells)*1e6_dp
      eov(0)%pdenm = pdenm(:,1:Neir_cells)*1e6_dp
      eov(0)%pdeni = pdeni(:,1:Neir_cells)*1e6_dp

      eov(0)%edena = edena(:,1:Neir_cells)*1e6_dp*elcha
      eov(0)%edenm = edenm(:,1:Neir_cells)*1e6_dp*elcha
      eov(0)%edeni = edeni(:,1:Neir_cells)*1e6_dp*elcha

      eov(0)%vxdena = vxdena(:,1:Neir_cells)*10._dp
      eov(0)%vxdenm = vxdenm(:,1:Neir_cells)*10._dp
      eov(0)%vxdeni = vxdeni(:,1:Neir_cells)*10._dp

      eov(0)%vydena = vydena(:,1:Neir_cells)*10._dp
      eov(0)%vydenm = vydenm(:,1:Neir_cells)*10._dp
      eov(0)%vydeni = vydeni(:,1:Neir_cells)*10._dp
  
      eov(0)%vzdena = vzdena(:,1:Neir_cells)*10._dp
      eov(0)%vzdenm = vzdenm(:,1:Neir_cells)*10._dp
      eov(0)%vzdeni = vzdeni(:,1:Neir_cells)*10._dp

! these are the direct coupling data

! conversion from Amp/cm3 to part/m3/s (and Amp to part/s)
      eov(0)%pael = pael(1:Neir_cells)*1e6_dp/elcha
      eov(0)%pmel = pmel(1:Neir_cells)*1e6_dp/elcha
      eov(0)%piel = piel(1:Neir_cells)*1e6_dp/elcha

      eov(0)%papl = papl(:,1:Neir_cells)*1e6_dp/elcha
      eov(0)%pmpl = pmpl(:,1:Neir_cells)*1e6_dp/elcha
      eov(0)%pipl = pipl(:,1:Neir_cells)*1e6_dp/elcha

! conversion from Watt/cm3 to Watt/m3
      eov(0)%eael = eael(1:Neir_cells)*1e6_dp
      eov(0)%emel = emel(1:Neir_cells)*1e6_dp
      eov(0)%eiel = eiel(1:Neir_cells)*1e6_dp

! total energy sources for ions, for checking purposes
      eov(0)%eapl = eapl(1:Neir_cells)*1e6_dp
      eov(0)%empl = empl(1:Neir_cells)*1e6_dp
      eov(0)%eipl = eipl(1:Neir_cells)*1e6_dp

! energy sources resolved by ions, stored as additionnal tallies

      eov(0)%eaplr = addv(1:nplsi,1:Neir_cells)*1e6_dp
      eov(0)%emplr = addv(nplsi+1:2*nplsi,1:Neir_cells)*1e6_dp
      eov(0)%eiplr = addv(2*nplsi+1:3*nplsi,1:Neir_cells)*1e6_dp


! conversion from Amp.g.cm.s-1/cm3 to kg/m.s-1/m3/s
      eov(0)%mapl = mapl(:,1:Neir_cells)*10._dp/elcha
      eov(0)%mmpl = mmpl(:,1:Neir_cells)*10._dp/elcha
      eov(0)%mipl = mipl(:,1:Neir_cells)*10._dp/elcha


! then surface tallies

! Amp to part/s
      eos(0)%potat  = potat(:,1:Nsurf_tal)/elcha
      eos(0)%prfaat = prfaat(:,1:Nsurf_tal)/elcha
      eos(0)%prfmat = prfmat(:,1:Nsurf_tal)/elcha
      eos(0)%prfpat = prfpat(:,1:Nsurf_tal)/elcha
      eos(0)%potml  = potml(:,1:Nsurf_tal)/elcha
      eos(0)%prfaml = prfaml(:,1:Nsurf_tal)/elcha
      eos(0)%prfmml = prfmml(:,1:Nsurf_tal)/elcha
      eos(0)%prfpml = prfpml(:,1:Nsurf_tal)/elcha
      eos(0)%potio  = potio(:,1:Nsurf_tal)/elcha
      eos(0)%potpl  = potpl(:,1:Nsurf_tal)/elcha

! no conversion needed here (Watt)
      eos(0)%eotat  = eotat(:,1:Nsurf_tal)
      eos(0)%erfaat = erfaat(:,1:Nsurf_tal)
      eos(0)%erfmat = erfmat(:,1:Nsurf_tal)
      eos(0)%erfiat = erfiat(:,1:Nsurf_tal)
      eos(0)%erfpat = erfpat(:,1:Nsurf_tal)
      eos(0)%eotml  = eotml(:,1:Nsurf_tal)
      eos(0)%erfaml = erfaml(:,1:Nsurf_tal)
      eos(0)%erfmml = erfmml(:,1:Nsurf_tal)
      eos(0)%erfiml = erfiml(:,1:Nsurf_tal)
      eos(0)%erfpml = erfpml(:,1:Nsurf_tal)
      eos(0)%eotio  = eotio(:,1:Nsurf_tal)
      eos(0)%erfaio = erfaio(:,1:Nsurf_tal)
      eos(0)%erfmio = erfmio(:,1:Nsurf_tal)
      eos(0)%erfiio = erfiio(:,1:Nsurf_tal)
      eos(0)%erfpio = erfpio(:,1:Nsurf_tal)
      eos(0)%eotpl  = eotpl(:,1:Nsurf_tal)

! Amp to part/s
      eos(0)%sptat  = sptat(:,1:Nsurf_tal)/elcha
      eos(0)%sptml  = sptml(:,1:Nsurf_tal)/elcha
      eos(0)%sptio  = sptio(:,1:Nsurf_tal)/elcha
      eos(0)%sptpl  = sptpl(:,1:Nsurf_tal)/elcha
      eos(0)%spttot = spttot(1:Nsurf_tal)/elcha

      eos(0)%spump  = spump(:,1:Nsurf_tal)/elcha


      if (interface_test_mode) then
! checks on recycling fluxes
        call styx_check_ion_fluxes
! checks on ion energy sources
        call styx_check_ion_energy_sources
      endif

! this assumes no mesh condensation
      ! removed without too much thinking TK3X
   !   if (NTRI .ne. NSBOX_TAL) then
   !   	write(iunout,*) 'NTRI not equal to NSBOX_TAL in if4cop'
   !   	call eirene_exit_own(1)
   !   endif

!     get total flux strength for sum of strata (Amp) = total number of D entering the system/s

      if (NSTRAI>1) then
	fluxt_eirene=FLUXT(0)
      else
      	fluxt_eirene=FLUXT(1)
      endif


!     get actual number of calculated particles

      if (.not.allocated(XMCP_styx)) allocate(XMCP_styx(0:nstrai))
      XMCP_styx = XMCP
      nprs_styx = nprs



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!  statistical noise output at last iteration !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      if (curiter == numiter) then
      
! first statistical noise on integral quantities   
! covariances 
! (format : i, covariance12, sigma1, sigma2, sigmas in physical units in constrast with sgms,sgmws)

	open(unit=555,file='covariance.integrals',status='replace')

	do i=1,NSIGCI
		write(555,'(i3,3es12.4)') i,sgmcs(0,i),sgmcs(1,i),sgmcs(2,i)
        enddo

	close(555)

        inquire(File='statistical_noise',exist=dir_e)
    	if(.not.dir_e) then
       		call system("mkdir statistical_noise")
    	end if

!       volume tallies

      	open(unit=660,file='statistical_noise/sigma_Nat.txt',
     .   	status='replace')
      	open(unit=661,file='statistical_noise/sigma_Nmol.txt',
     .		status='replace')
      	open(unit=662,file='statistical_noise/sigma_Ntion.txt',
     .		status='replace')

        open(unit=663,file='statistical_noise/sigma_Eat.txt',
     . 		status='replace')
        open(unit=664,file='statistical_noise/sigma_Em.txt',
     .		status='replace')
        open(unit=665,file='statistical_noise/sigma_Etion.txt',
     .		status='replace')

        open(unit=666,file='statistical_noise/sigma_Snat.txt',
     .		status='replace')
      	open(unit=667,file='statistical_noise/sigma_Snmol.txt',
     .		status='replace')
      	open(unit=668,file='statistical_noise/sigma_Sntion.txt',
     .		status='replace')

      	open(unit=669,file='statistical_noise/sigma_Smat.txt',
     .		status='replace')
        open(unit=670,file='statistical_noise/sigma_Smmol.txt',
     .		status='replace')
        open(unit=671,file='statistical_noise/sigma_Smtion.txt',
     .		status='replace')

        open(unit=672,file='statistical_noise/sigma_SUeat.txt',
     .		status='replace')
        open(unit=673,file='statistical_noise/sigma_SUemol.txt',
     .		status='replace')
      	open(unit=674,file='statistical_noise/sigma_SUetion.txt',
     .		status='replace')

      	open(unit=675,file='statistical_noise/sigma_SUiat.txt',
     .		status='replace')
        open(unit=676,file='statistical_noise/sigma_SUimol.txt',
     .		status='replace')
        open(unit=677,file='statistical_noise/sigma_SUition.txt',
     .		status='replace')
        

	do ital=1,NSIGVI

		select case (iih(ital))

			case(1)
				do itri=1,NTRI-1
      					write(660,'(3es15.7)') sigma(ital,itri),
     .					sigma(ital,itri),sigma(ital,itri)
				enddo
				close(660)
			case(2)
				do itri=1,NTRI-1
					write(661,'(3es15.7)') sigma(ital,itri),
     .					sigma(ital,itri),sigma(ital,itri)
				enddo
				close(661)
			case(3)
				do itri=1,NTRI-1
					write(662,'(3es15.7)') sigma(ital,itri),
     .					sigma(ital,itri),sigma(ital,itri)
				enddo
				close(662)
			case(5)
				do itri=1,NTRI-1
					write(663,'(3es15.7)') sigma(ital,itri),
     .					sigma(ital,itri),sigma(ital,itri)
				enddo
				close(663)
			case(6)
				do itri=1,NTRI-1
					write(664,'(3es15.7)') sigma(ital,itri),
     .					sigma(ital,itri),sigma(ital,itri)
				enddo
				close(664)
			case(7)
				do itri=1,NTRI-1
					write(665,'(3es15.7)') sigma(ital,itri),
     .					sigma(ital,itri),sigma(ital,itri)
				enddo
				close(665)
			case(9)
				do itri=1,NTRI-1
					write(666,'(3es15.7)') sigma(ital,itri),
     .					sigma(ital,itri),sigma(ital,itri)
				enddo
				close(666)
			case(15)
				do itri=1,NTRI-1
					write(667,'(3es15.7)') sigma(ital,itri),
     .					sigma(ital,itri),sigma(ital,itri)
				enddo
				close(667)
			case(21)
				do itri=1,NTRI-1
					write(668,'(3es15.7)') sigma(ital,itri),
     .					sigma(ital,itri),sigma(ital,itri)
				enddo
				close(668)
			case(97)
				do itri=1,NTRI-1
       					write(669,'(3es15.7)') sigma(ital,itri),
     .					sigma(ital,itri),sigma(ital,itri)
				enddo
				close(669)
			case(98)
				do itri=1,NTRI-1
					write(670,'(3es15.7)') sigma(ital,itri),
     .					sigma(ital,itri),sigma(ital,itri)
				enddo
				close(670)
			case(99)
				do itri=1,NTRI-1
					write(671,'(3es15.7)') sigma(ital,itri),
     .					sigma(ital,itri),sigma(ital,itri)
				enddo
				close(671)
			case(33)
				do itri=1,NTRI-1
					write(672,'(3es15.7)') sigma(ital,itri),
     .					sigma(ital,itri),sigma(ital,itri)
				enddo
				close(672)
			case(39)
				do itri=1,NTRI-1
					write(673,'(3es15.7)') sigma(ital,itri),
     .					sigma(ital,itri),sigma(ital,itri)
				enddo
				close(673)
			case(45)
				do itri=1,NTRI-1
					write(674,'(3es15.7)') sigma(ital,itri),
     .					sigma(ital,itri),sigma(ital,itri)
				enddo
				close(674)
			case(38)
				do itri=1,NTRI-1
					write(675,'(3es15.7)') sigma(ital,itri),
     .					sigma(ital,itri),sigma(ital,itri)
				enddo
				close(675)
			case(44)
				do itri=1,NTRI-1
					write(676,'(3es15.7)') sigma(ital,itri),
     .					sigma(ital,itri),sigma(ital,itri)
				enddo
				close(676)
			case(50)
				do itri=1,NTRI-1
					write(677,'(3es15.7)') sigma(ital,itri),
     .					sigma(ital,itri),sigma(ital,itri)
				enddo
				close(677)

			end select
	enddo

      	
	close(661)
	close(662)
	close(663)
	close(664)
	close(665)
	close(666)
	close(667)
	close(668)
	close(669)
	close(670)
	close(671)
	close(672)
	close(673)
	close(674)
	close(675)
	close(676)
	close(677)


!     report on computationnal time for each nodes in case of a parallel run

     	open(unit=555,file='parallel_timing.log',status='replace')

	do i=1,nprs
       		write(555,'(2i6,es12.4)') i,nstrpe(i-1),time_para(i)
      	enddo

     	close(555)


      endif

!     data on computational times for printout in styx_process_fluxes

      ptimax = maxval(time_para)
      ptimin = minval(time_para)

      rtimax = maxloc(time_para,1)
      rtimin = minloc(time_para,1)

      ptmean = sum(time_para)/nprs


!     further output to check data exchanges (styx -> eirene)

      check_exchg=.true.

      if (check_exchg) then

      	open(unit=660,file='DEIN_tri.txt',status='replace')
      	open(unit=661,file='TEIN_tri.txt',status='replace')
      	open(unit=662,file='TIIN_tri.txt',status='replace')
        open(unit=663,file='TABDS1_tri.txt',status='replace')
        open(unit=664,file='TABRC1_tri.txt',status='replace')
        open(unit=665,file='TABCX3_tri.txt',status='replace')
        open(unit=667,file='vpar_tri.txt',status='replace')

      	do itri=1,NTRI-1
      	    write(660,'(3es15.7)') DEIN(itri),DEIN(itri),DEIN(itri)
  	    write(661,'(3es15.7)') TEIN(itri),TEIN(itri),TEIN(itri)
  	    write(662,'(3es15.7)') TIIN(1,itri),TIIN(1,itri),TIIN(1,itri)
  	   ! ionisation reaction rate #1 (D->D+), s-1
	    write(663,'(3es15.7)') TABDS1(1,itri),TABDS1(1,itri),
     1          TABDS1(1,itri)
           ! recombination reaction rate #1 (D+-D), s-1
  	    write(664,'(3es15.7)') TABRC1(1,itri),TABRC1(1,itri),
     1          TABRC1(1,itri)
   	   ! charge exchange reaction rate #1 (D+ <-> D)
  	    write(665,'(3es15.7)') TABCX3(1,itri,1),TABCX3(1,itri,1),
     1          TABCX3(1,itri,1)
      	    vpar_itri=VXIN(1,itri)*BXIN(itri)+VYIN(1,itri)*BYIN(itri)+
     1          VZIN(1,itri)*BZIN(itri)
  	    write(667,'(3es15.7)') vpar_itri,vpar_itri,vpar_itri
      	enddo

      	close(660)
      	close(661)
      	close(662)
      	close(663)
      	close(664)
      	close(665)
        close(667)
      endif

      END
