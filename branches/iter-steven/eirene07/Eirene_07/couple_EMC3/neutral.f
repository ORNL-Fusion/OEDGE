      SUBROUTINE RUN_NEUTRAL()
      USE PARALLEL       
      USE PHYSICAL_CELL
      USE BOUNDARY_COND
      USE SOURCE_V_PL
      USE PLASMA_PARAM
      USE ION_SPECIES
      USE NEUTRAL_TRANSPORT

      IMPLICIT NONE
      REAL*8, PARAMETER :: RLF=0.5
      INTEGER           :: IPLS,ITYP,ICA,ICB,I
      REAL*8            :: FLUX,FLUX_M,SP
      LOGICAL           :: LRUN_FIRST

C  Call Eirene
      CALL EIRENE_EMC3_ENTRY(0.,.FALSE.)
      IF(MYPE == 0 ) THEN

C  Post_processing
      FLUX  = PFLUX_TOTAL(1)
      FLUX_M= FLUX/(1.6E-19*1.67E-24*ATMASS(1))

      OUT_N0(1:NCELL_N2,4  ) = OUT_N0(1:NCELL_N2,4  )*FLUX_M
      OUT_N0(1:NCELL_N2,5:6) = OUT_N0(1:NCELL_N2,5:6)*FLUX
      DO I=1,NC_PL
      IF( ZMACH0(I) > 0.) THEN 
        OUT_N0(I,4) = OUT_N0(I,4)/SOUND_S(I)
      ELSE 
        OUT_N0(I,4) =-OUT_N0(I,4)/SOUND_S(I)
      ENDIF 
      ENDDO 
C First run
      DO I=1,MIN(100,NC_PL)
        LRUN_FIRST = VSOUP0(I,1).EQ.0.
        IF(.NOT.LRUN_FIRST) EXIT
      ENDDO

      IF( LRUN_FIRST ) THEN
         VSOUP0(1:NC_PL,1  ) = OUT_N0(1:NC_PL,1)
         VSOUE0(1:NC_PL,0,0) = OUT_N0(1:NC_PL,2)      
         VSOUE0(1:NC_PL,0,1) = OUT_N0(1:NC_PL,3)      
         VSOUM0(1:NC_PL,1  ) = OUT_N0(1:NC_PL,4)
      ELSE
         VSOUP0(1:NC_PL,1  ) =    RLF *VSOUP0(1:NC_PL,1  )  
     &+                       (1.-RLF)*OUT_N0(1:NC_PL,1  )      
         VSOUE0(1:NC_PL,0,0) =    RLF *VSOUE0(1:NC_PL,0,0)  
     &+                       (1.-RLF)*OUT_N0(1:NC_PL,2  )      
         VSOUE0(1:NC_PL,0,1) =    RLF *VSOUE0(1:NC_PL,0,1)  
     &+                       (1.-RLF)*OUT_N0(1:NC_PL,3  )      
         VSOUM0(1:NC_PL,1  ) =    RLF *VSOUM0(1:NC_PL,1  )  
     &+                       (1.-RLF)*OUT_N0(1:NC_PL,4  )      
      ENDIF 

      SP_MAIN_NEW = 0.
      ICA= NCELL_AB(0)+1
      ICB= NCELL_AB(1)
      DO I=ICA,ICB
        SP_MAIN_NEW    = SP_MAIN_NEW+OUT_N0(I,1)*VOLCEL_N0(I)
      ENDDO
      IF(SP_MAIN_OLD.NE.0.) THEN
        SP_MAIN_NEW = RLF*SP_MAIN_OLD+(1.-RLF)*SP_MAIN_NEW
      ENDIF
      SP_MAIN_OLD   = SP_MAIN_NEW

       OPEN (40)
       WRITE(40,'(6E12.4)') VSOUP0(1:NC_PL,1)
       CLOSE(40)

       OPEN (42)
       WRITE(42,'(6E12.4)') VSOUE0(1:NC_PL,0,0)
       CLOSE(42)

       OPEN (43)
       WRITE(43,'(6E12.4)') VSOUE0(1:NC_PL,0,1)
       CLOSE(43)

       OPEN (47)
       WRITE(47,'(6E12.4)') VSOUM0(1:NC_PL,1)
       CLOSE(47)

       OPEN (99,FILE='DENSITY_A')
       WRITE(99,'(6E12.4)') OUT_N0(1:NCELL_N2,5)
       CLOSE(99)  

       OPEN (99,FILE='DENSITY_M')
       WRITE(99,'(6E12.4)') OUT_N0(1:NCELL_N2,6)
       CLOSE(99)  
          
       WRITE(6,*)'    Cell    S_P(Amp.)  S_E_E(Watt) S_E_I(Watt)   
     .    SM_P      SM_N'
          SP          =0.
          SE_SUM_N0_E =0.
          SE_SUM_N0_I =0.
          SM_SUMP_ALL =0.
          SM_SUMN_ALL =0.
          DO I=1,NC_PL
           SP          = SP          + OUT_N0(I,1)*VOLCEL(I)
           SE_SUM_N0_E = SE_SUM_N0_E + OUT_N0(I,2)*VOLCEL(I)
           SE_SUM_N0_I = SE_SUM_N0_I + OUT_N0(I,3)*VOLCEL(I)
           IF( ZMACH0(I) > 0.) THEN 
           SM_SUMP_ALL = SM_SUMP_ALL + OUT_N0(I,4)*VOLCEL(I)
           ELSE
           SM_SUMN_ALL = SM_SUMN_ALL + OUT_N0(I,4)*VOLCEL(I)
           ENDIF
          ENDDO
          WRITE(6,'(2I6,5E12.4)') 1,NC_PL,SP*FLUX,SE_SUM_N0_E*FLUX
     .,   SE_SUM_N0_I*FLUX,SM_SUMP_ALL,SM_SUMN_ALL

          DO ITYP=1,NCTYPE_DEF
          SP          =0.
          SE_SUM_N0_E =0.
          SE_SUM_N0_I =0.
          ICA= NCELL_AB(ITYP-1)+1
          ICB= NCELL_AB(ITYP)
          DO I=ICA,ICB
           SP          = SP          + OUT_N0(I,1)*VOLCEL_N0(I)
           SE_SUM_N0_E = SE_SUM_N0_E + OUT_N0(I,2)*VOLCEL_N0(I)
           SE_SUM_N0_I = SE_SUM_N0_I + OUT_N0(I,3)*VOLCEL_N0(I)
          ENDDO
          WRITE(6,'(2I6,5E12.4)') ICA,ICB,SP*FLUX,SE_SUM_N0_E*FLUX
     .,   SE_SUM_N0_I*FLUX,0.,0.                   
          ENDDO 
      ENDIF 
      IF(NPRS>1) CALL BROADCAST_SOURCE_V_PL()

      RETURN
      END SUBROUTINE RUN_NEUTRAL

