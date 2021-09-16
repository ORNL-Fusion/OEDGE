      MODULE CSTEP
c  step functions for specification of incident fluxes on
c                 surfaces, with spatial resolution along surface
C    FLSTEP: total ion particle flux at sheath entrance (s.e.)
C    ELSTEP: total ion energy flux   at sheath entrance (s.e.)
C    RRSTEP: abszissa, i.e. arclength along target (cm)
c
c  JET-2005: patch 1
c  9.9.05:  new step function parameters: 
c    festep,fistep,shstep,vpstep,mcstep
c    festep: electron energy flux at sheath entrance (s.e.) is festep*testep
c            e.g.: flstep=4.5= 2 + 2.5 (2: maxwellian energy flux, 2.5: sheath)
c    fistep: thermal part of ion energy flux at s.e. is fistep*tistep
c            e.g.: fistep= 2.5 (drifting maxwellian thermal energy flux,
c                               with Mach=1 (or: ion thermal veloc=1 ?).
c    shstep: sheath multiplier:  sheath potential is shstep*testep
c            e.g.: shstep= 2.5 (hydrogen, M=1, Te=Ti, single fluid)
c            e.g.: shstep= 2.8 (deuteron, M=1, Te=Ti, single fluid)
c    vpstep: parallel to B-field drift velocity (cm/s) at s.e.
c    mcstep: mach number of parallel flow at s.e.
c
c   there is now a certain redundancy of information on plasma
c   conditions along target surfaces. This may facilitate consistency
c   checks for boundary conditions in case of coupling to edge plasma codes.
C  11.11.05:  ve and eltot introduced, to provide energy flux step function
c
c
      USE PRECISION
      USE PARMMOD

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: ALLOC_CSTEP, DEALLOC_CSTEP, INIT_CSTEP

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R FLSTEP(:,:,:), ELSTEP(:,:,:), FLTOT(:,:), ELTOT(:,:),
     R VF(:,:,:),  VE(:,:,:),   QUOT(:,:,:),   ADD(:,:,:),
     R QUOTI(:,:,:),  ADDIV(:,:,:),
     R TESTEP(:,:),   TISTEP(:,:,:), RRSTEP(:,:),
     R VXSTEP(:,:,:), VYSTEP(:,:,:), VZSTEP(:,:,:),
     R DISTEP(:,:,:), FESTEP(:,:),   FISTEP(:,:,:),
     R SHSTEP(:,:),   VPSTEP(:,:,:), MCSTEP(:,:,:)

      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I IRSTEP(:,:), IPSTEP(:,:), ITSTEP(:,:),
     I IASTEP(:,:), IBSTEP(:,:), IGSTEP(:,:),
     I ISTUF(:),
     I NSMAX(:),    NSPSTI(:),   NSPSTE(:)

      INTEGER, PUBLIC, SAVE ::
     .   NSTPP1, NSTPP2, NSTPP3 , NSTPP4

      CONTAINS


      SUBROUTINE ALLOC_CSTEP

      IF (ALLOCATED(FLSTEP)) RETURN

      NSTPP1=(NSPZ+1)*NSTEP*NGITT
      NSTPP2=NSTEP*NGITT
      NSTPP3=NPLS*NSTEP*NGITT
      NSTPP4=(NSPZ+1)*NSTEP

      ALLOCATE (FLSTEP(0:NSPZ,NSTEP,NGITT))
      ALLOCATE (ELSTEP(0:NSPZ,NSTEP,NGITT))
      ALLOCATE (FLTOT (0:NSPZ,NSTEP))
      ALLOCATE (VF    (0:NSPZ,NSTEP,NGITT))
c   next 2 tallies added nov. 05      !dr
      ALLOCATE (ELTOT (0:NSPZ,NSTEP))
      ALLOCATE (VE    (0:NSPZ,NSTEP,NGITT))
      ALLOCATE (QUOT  (0:NSPZ,NSTEP,NGITT))
      ALLOCATE (ADD   (0:NSPZ,NSTEP,NGITT))
      ALLOCATE (QUOTI (0:NSPZ,NSTEP,NGITT))
      ALLOCATE (ADDIV (0:NSPZ,NSTEP,NGITT))
      ALLOCATE (TESTEP(NSTEP,NGITT))
      ALLOCATE (TISTEP(NPLSTI,NSTEP,NGITT))
      ALLOCATE (RRSTEP(NSTEP,NGITT))
      ALLOCATE (VXSTEP(NPLSV,NSTEP,NGITT))
      ALLOCATE (VYSTEP(NPLSV,NSTEP,NGITT))
      ALLOCATE (VZSTEP(NPLSV,NSTEP,NGITT))
      ALLOCATE (DISTEP(NPLS,NSTEP,NGITT))
c  next 5 tallies added sept. 05     !dr
      ALLOCATE (FESTEP(NSTEP,NGITT))
      ALLOCATE (SHSTEP(NSTEP,NGITT))
      ALLOCATE (FISTEP(NPLS,NSTEP,NGITT))
      ALLOCATE (MCSTEP(NPLS,NSTEP,NGITT))
      ALLOCATE (VPSTEP(NPLSV,NSTEP,NGITT))

      ALLOCATE (IRSTEP(NSTEP,NGITT))
      ALLOCATE (IPSTEP(NSTEP,NGITT))
      ALLOCATE (ITSTEP(NSTEP,NGITT))
      ALLOCATE (IASTEP(NSTEP,NGITT))
      ALLOCATE (IBSTEP(NSTEP,NGITT))
      ALLOCATE (IGSTEP(NSTEP,NGITT))
      ALLOCATE (ISTUF(NSTEP))
      ALLOCATE (NSMAX(NSTEP))
      ALLOCATE (NSPSTI(NSTEP))
      ALLOCATE (NSPSTE(NSTEP))

      WRITE (55,'(A,T25,I15)')
     .       ' CSTEP ',(8*NSTPP1+2*NSTPP4+4*NSTPP2+3*NSTPP3 +
     .                  (NPLSTI+4*NPLSV)*NSTPP2)*8 +
     .                 (6*NSTPP2+4*NSTEP)*4

      CALL INIT_CSTEP

      RETURN
      END SUBROUTINE ALLOC_CSTEP


      SUBROUTINE DEALLOC_CSTEP

      IF (.NOT.ALLOCATED(FLSTEP)) RETURN

      DEALLOCATE (FLSTEP)
      DEALLOCATE (ELSTEP)
      DEALLOCATE (FLTOT)
      DEALLOCATE (VF)
      DEALLOCATE (ELTOT)
      DEALLOCATE (VE)
      DEALLOCATE (QUOT)
      DEALLOCATE (ADD)
      DEALLOCATE (QUOTI)
      DEALLOCATE (ADDIV)
      DEALLOCATE (TESTEP)
      DEALLOCATE (TISTEP)
      DEALLOCATE (RRSTEP)
      DEALLOCATE (VXSTEP)
      DEALLOCATE (VYSTEP)
      DEALLOCATE (VZSTEP)
      DEALLOCATE (DISTEP)
      DEALLOCATE (FESTEP)
      DEALLOCATE (SHSTEP)
      DEALLOCATE (FISTEP)
      DEALLOCATE (MCSTEP)
      DEALLOCATE (VPSTEP)

      DEALLOCATE (IRSTEP)
      DEALLOCATE (IPSTEP)
      DEALLOCATE (ITSTEP)
      DEALLOCATE (IASTEP)
      DEALLOCATE (IBSTEP)
      DEALLOCATE (IGSTEP)
      DEALLOCATE (ISTUF)
      DEALLOCATE (NSMAX)
      DEALLOCATE (NSPSTI)
      DEALLOCATE (NSPSTE)

      RETURN
      END SUBROUTINE DEALLOC_CSTEP


      SUBROUTINE INIT_CSTEP

      FLSTEP = 0._DP
      ELSTEP = 0._DP
      FLTOT  = 0._DP
      VF     = 0._DP
      ELTOT  = 0._DP
      VE     = 0._DP
      QUOT   = 0._DP
      ADD    = 0._DP
      QUOTI  = 0._DP
      ADDIV  = 0._DP
      TESTEP = 0._DP
      TISTEP = 0._DP
      RRSTEP = 0._DP
      VXSTEP = 0._DP
      VYSTEP = 0._DP
      VZSTEP = 0._DP
      DISTEP = 0._DP
      FESTEP = 0._DP
      FISTEP = 0._DP
      SHSTEP = 0._DP
      MCSTEP = 0._DP
      VPSTEP = 0._DP

      IRSTEP = 0
      IPSTEP = 0
      ITSTEP = 0
      IASTEP = 0
      IBSTEP = 0
      IGSTEP = 0
      ISTUF  = 0
      NSMAX  = 0
      NSPSTI = 0
      NSPSTE = 0

      RETURN
      END SUBROUTINE INIT_CSTEP

      END MODULE CSTEP
