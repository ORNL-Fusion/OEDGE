subroutine eirene_get_fluxes()
  use all_variables, only : global_parameters, interp_data2, reference_parameters
  use Mphysics
  USE EIRMOD_PRECISION
  USE EIRMOD_PARMMOD
  USE EIRMOD_COMUSR
  USE EIRMOD_CSTEP
  USE EIRMOD_CGRID
  USE EIRMOD_CINIT
  USE EIRMOD_COMSOU
  USE EIRMOD_CTRIG
  use eirmod_comprt, only : iunout
  use eirmod_ccona
  use styx2eirene

  IMPLICIT NONE
  REAL(DP) :: EIRENE_STEP
  INTEGER :: ITRI, ISIDE, NBIN, ISOR 
  integer :: ITEC1, ITEC2, ITEC3, ISTEP, INDSRF, IS1, IERROR, IPLS
  integer :: IP1,IP2,istra,isp
  INTEGER :: EIRENE_IDEZ
  INTEGER, ALLOCATABLE :: KSTEP(:), INOSRC(:), IPLAN(:), IPLEN(:)
  real(dp) :: DELR, FL,dens,Te,Ti,csx,csy,csz,pflux,fl2
  integer :: icount,isurf,itor
  integer :: icount_neg(3)
  character(50) :: filename

  ! define number of step functions  = global_parameters%n_species
  ! NSTEP defined by number of stratums for which indsrc=6 at the begining of block 9
  ! ISTEP has to be specified in infcop for each recycling strata, variable SORIND
  ! define iZ range for each species: for the moment Z
  ! that is IPLAN(ISTEP) = isp/ Z=1 (singly charged, already available from infcop)
  !         IPLEN(ISTEP) = isp/ Z=Z
  ! then recalculate fluxes for each ion, see Hugo (see how ti, etc ... transfer and what is used to
  ! compute the sheath potential 
  ! check consistency of array sizes with mpi_bcasts !!

  ! the loop below was originally intended as a loop on strata, which makes sense

  icount=0
  icountTe=0
  icountTi=0
  icountDe=0 


  ALLOCATE (KSTEP(NSTEP))
  ALLOCATE (INOSRC(NSTEP))
  ALLOCATE (IPLAN(NSTEP))
  ALLOCATE (IPLEN(NSTEP))

  INOSRC = 0


  do istra=1,Nrecyc
     isp=istra
     fl2=0._dp
     ! loop on set of surfaces per strata (1 in 2D)
     do itor=1,NSRFSI(istra)
       ! check that a step function is specified for this recycling strata (use isurf =1 to get ISOR)
       ISOR=SORLIM(itor,istra)
       ITEC1=EIRENE_IDEZ(ISOR,1,4)
       ITEC2=EIRENE_IDEZ(ISOR,2,4)
       ITEC3=EIRENE_IDEZ(ISOR,3,4)

       ! issue if no step function required for this strata
       IF ((ITEC1 /= 4).AND.(ITEC2 /= 4).AND.(ITEC3 /= 4)) then
         write(*,*) 'No step function required for recycling strata #',istra,' and surface set # =',itor
         stop
       endif

       ISTEP=SORIND(itor,istra)

       IF (ISTEP.EQ.0) THEN
         WRITE (6,*) 'ERROR IN PRIMARY SOURCE DATA '
         WRITE (6,*) 'STEPFUNCTION REQUESTED FOR SOURCE SURFACE '
         WRITE (6,*) 'NO. ', INSOR(isurf,1),' BUT SORIND.EQ.0.'
         CALL EXIT(1)
       ELSEIF (ISTEP.GT.NSTEP) THEN
         CALL EIRENE_MASPRM('NSTEP',5,NSTEP,'ISTEP',5,ISTEP,IERROR)
         CALL EXIT(1)
!     elseif (istep /= istra) then
!        write(*,*) ' Error in strata numbering, istra different from istep, eirene_get_fluxes'
!        call exit(1)
       endif

       INDSRF=INSOR(itor,1)

       IF (INDSRF < 0) INDSRF=NLIM+ABS(INDSRF)
     ! indexes of ion species for this element (hence stratum)
       if (nspez(istra) <= 0) then
          iplan((istra-1)*Ntor_cells+itor)=singly_charged_ions(isp)
          iplen((istra-1)*Ntor_cells+itor)=iplan((istra-1)*Ntor_cells+itor)+global_parameters%element_list(isp)%Z-1
       else
          write(*,*) 'warning, species sampling turned off for stratum #',istra
          iplan((istra-1)*Ntor_cells+itor)=isp
          iplen((istra-1)*Ntor_cells+itor)=isp
       endif

       do ipls = iplan((istra-1)*Ntor_cells+itor),iplen((istra-1)*Ntor_cells+itor)

         kstep(istep) = 0
         write(filename,"(A5,I0)") "wall_",ipls
         open(unit=104,file=trim(adjustl(filename)),status='unknown')

         ! now calculate fluxes
         do isurf=1,Nsou
           ! get vertex number for the segment making the wall
           IP1   = recsurf(isurf)%v1
           IP2   = recsurf(isurf)%v2
           ! corresponding triangle number and side
           ITRI  = recsurf(isurf)%ITRI
           ISIDE = recsurf(isurf)%ISIDE

           call styx_plasma_parameters_on_wall(&
                ip1,ip2,itri,iside,isurf,itor,ipls,isp,dens,Te,Ti,csx,csy,csz,pflux)
           ! store for diagnostics/checks
           pflux_in(itri,iside,itor)=pflux 


           write(104,150) itri, iside, Te, Ti, dens
150        format(I6,I6,es15.7,es15.7,es15.7)

           if (Te<0._dp) then
              write(*,*) 'catched ...'
           endif

           if (INMTI(ISIDE,ITRI) /=0) then
             IF (KSTEP(ISTEP) == 0) RRSTEP(ISTEP,1) = 0._DP
             KSTEP(ISTEP) = KSTEP(ISTEP) + 1
             INOSRC(ISTEP) = istra
             IS1 = ISIDE + 1
             IF (IS1.GT.3) IS1=1             
             IRSTEP(ISTEP,KSTEP(ISTEP))=ITRI
             IPSTEP(ISTEP,KSTEP(ISTEP))=ISIDE
             ITSTEP(ISTEP,KSTEP(ISTEP))=itor
             IASTEP(ISTEP,KSTEP(ISTEP))=0
             IBSTEP(ISTEP,KSTEP(ISTEP))=1
             DELR =  SQRT((XTRIAN(NECKE(ISIDE,ITRI))-XTRIAN(NECKE(IS1,ITRI)))**2+(YTRIAN(NECKE(ISIDE,ITRI))-YTRIAN(NECKE(IS1,ITRI)))**2)
             RRSTEP(ISTEP,KSTEP(ISTEP)+1)=RRSTEP(ISTEP,KSTEP(ISTEP)) + DELR         

             TESTEP(ISTEP,KSTEP(ISTEP)) = Te   

              ! plasma background data need to be stored in all step functions             
             TISTEP(ipls,:,KSTEP(ISTEP)) = Ti
             DISTEP(ipls,:,KSTEP(ISTEP)) = dens
             VXSTEP(ipls,:,KSTEP(ISTEP)) = csx
             VYSTEP(ipls,:,KSTEP(ISTEP)) = csy
             VZSTEP(ipls,:,KSTEP(ISTEP)) = csz

             FLSTEP(ipls,ISTEP,KSTEP(ISTEP)) = ABS(pflux)/DELR

!!! sheath1D data
             if (Te > 0._dp) then
               sheath1D(kstep(istep))%tau = Ti/Te
             else
               sheath1D(kstep(istep))%tau = 0._dp
               write(*,*) 'Warning, Te not strictly positive, from eirene_get_fluxes'
               write(*,*) 'itri = ',itri,' ,iside = ',iside, 'Te = ',Te
             endif

             ! density in cm-3 consistent with ksi0
             if (dens > 0._dp) then
               sheath1D(kstep(istep))%ksi = sheath1D(kstep(istep))%ksi0 * sqrt(dens)        
             else
               sheath1D(kstep(istep))%ksi = 0._dp
               write(*,*) 'Warning, negative density, from eirene_get_fluxes'
               write(*,*) 'itri = ',itri,' ,iside = ',iside, 'Ne = ',dens
             endif

           else
             if (icount==0) then
               write(*,*) 'Inconsistency in wall properties, fluxes given on transparent surfaces'
               write(*,*) 'Check orientation of triangles ... '
             endif
             write(*,'(A7,i6,A10,i2,A10,i2,A14,es12.4)') 'itri = ',itri,' ,iside = ',iside,&
                 ' ,itor = ',itor,' ,iprop = ',inmti(iside,itri),' ,flx (Amp) = ',pflux_in(itri,iside,itor)
             icount=icount+1

           endif
         enddo
       close(104)
     enddo
     ! has to be done after all species in the current strata have been treated
     IF (KSTEP(ISTEP) > 0) THEN
       NBIN=KSTEP(ISTEP)+1
       FL=EIRENE_STEP(IPLAN(ISTEP),IPLEN(ISTEP),NBIN,ISTEP)
       ! accumulate fluxes on the different toroidal rings 
       fl2=fl2+FL
       !FLUX(INOSRC(ISTEP))=fl2
     END IF
    enddo ! itor
    ! is this generally valid ??
    FLUX(istra)=fl2
  end do ! istra



  if (icount > 0) call eirene_exit_own(1)


  ! debugging : output all step functions

  !      open(unit=666,file='step_functions_fast',status='replace')
  !
  !      do igitt=1,ngitt
  !      		write(666,'(2i6,6es12.4)') IRSTEP(1,igitt),IPSTEP(1,igitt), &
  !      			TISTEP(1,1,igitt),DISTEP(1,1,igitt),VXSTEP(1,1,igitt), &
  !      			VYSTEP(1,1,igitt),VZSTEP(1,1,igitt),FLSTEP(1,1,igitt)
  !      enddo

  !      close(666)


  DEALLOCATE (KSTEP)
  DEALLOCATE (INOSRC)

end subroutine eirene_get_fluxes
