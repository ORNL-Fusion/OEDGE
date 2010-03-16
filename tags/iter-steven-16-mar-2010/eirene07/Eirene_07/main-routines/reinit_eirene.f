      SUBROUTINE REINITIALIZATION_OF_EIRENE

      use PRECISION
      use PARMMOD
      use COUTAU
      use CRECH
      
      implicit none

      REAL(DP) :: dummy, H1RN_REINIT,  ranf_eirene_reinit,
     &     ranset_eirene_reinit, sheath_reinit 

C     reinitialization start
      call EIRENE_REINIT
      call SIGHA_REINIT
      
!  OUT Of USE
!      dummy = H1RN_REINIT
!      call H1RNV_REINIT

      dummy = ranf_eirene_reinit()
      dummy = ranset_eirene_reinit()

      call INIT_COUTAU_REINIT
      call CRECH_REINIT

      call STCOOR_REINIT
      call SAMVOL_REINIT
      call VELOCX_REINIT
      call VELOEL_REINIT

      call PL3D_REINIT
      call PLT2D_REINIT
      call PLTEIR_REINIT

      call REFLEC_REINIT

!pb      call UPTCOP_REINIT

C     reinitialization end


      end

