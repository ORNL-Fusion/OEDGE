c -*-Fortran-*-
c
! ammod begin
c
      subroutine global_hc_read_raw_data(version_code,maxrev)
      ! jdemod - comHc is part of HC_init_out_data through HC_out_storage_setup
      !Use ComHC ! HC constants.
      Use HC_Init_Out_Data ! Included to re-set hc_state, hc_density, hc_output, Number_HC_Species.
      implicit none  
      integer :: version_code,maxrev
c
c     jdemod - Reminder - need to bump the DIVIMP version number and add appropriate checks in here
c              when making changes to what gets passed in the raw file. Some changes (removals) have
c              been made here as a part of revising the HC module - note that this will break plots for 
c              old cases
c
c
c     read in HC related data
c
         Call iinout ('R HC_OPT',hc_follow_option,1)
         Call iinout ('R HC_HIG',hc_higher_hcs_option,1)
         Call iinout ('R HC_PRI',hc_evolution_model_primary,1)
         Call iinout ('R HC_SEC',hc_evolution_model_secondary,1)

         Call rinout ('R HC_STA', HC_State_List,
     >     Number_HC_Species)
         Call rinout ('R HC_DEN', HC_Density,
     >     maxnks*maxnrs*(Number_HC_Species))

         !
         ! jdemod - H_density was only writing numer_hc_species and not num_h_states
         !          until 6.47
         if (version_code.ge.(6*maxrev+47)) then 
            Call rinout ('R H_DEN', H_Density,
     >         maxnks*maxnrs*(Num_H_States))
         else
            Call rinout ('R H_DEN', H_Density,
     >          maxnks*maxnrs*(Number_H_Species))
         endif


         Call iinout ('R HC_OUT', HC_Output_List,
     >     Number_HC_Species+1)
         Call rinout ('R HC_WLK', HC_Walks,
     >     Max_Number_Walks*2)
c
c     jdemod - removed indices 0 and -1 from hc_fact
c
         Call rinout ('R HC_FactA', HC_Factor_A,
     >     Number_HC_Species)
         Call rinout ('R HC_FactB', HC_Factor_B,
     >     Number_HC_Species)
c         Call rinout ('R HC_FactA', HC_Factor_A,
c     >     Number_HC_Species+2)
c         Call rinout ('R HC_FactB', HC_Factor_B,
c     >     Number_HC_Species+2)


         Call rinout ('R HC_TIZS_CH', HC_TIZS_CH,
     >     maxnks*maxnrs*2)
c
c         Call rinout ('R HC_TIZS_C2', HC_TIZS_C2,
c     >     maxnks*maxnrs*2)
c         CALL RINOUT ('R FYTOT ', FYTOT,1)      

      return
      end
! ammod end
c
c
c
! ammod begin.
      Subroutine RDG_HC_Trans_Prob (Graph, Start_HC_Species,
     >  End_HC_Species, IErr)

      use mod_reader
      Implicit None
      Integer, Intent (Out) :: Start_HC_Species, End_HC_Species
      Character, Intent (Out) :: Graph*(*)

! *********************************************************************
! *                                                                   *
! * RDG_HC_Prob_Trans: Read in additional line for HC transition plots*
! *                                                                   *
! *********************************************************************

c     Include 'reader'
      
      ! Declare local variables.
      Integer :: IErr
      Character :: Message*72

      IErr = 0
      Message = 'End of file on unit 5'
      Do
         If (IBuf .eq. 0) Then
            Read (5,'(a512)',Err=9998,End=9998) Buffer
	 End If
         Write (9,'(1X,A72,1X,A6)') BUFFER,'RDG1'
         If (Buffer (1:1) .ne. '$') Then
	    Exit
	 End If
      End Do

      Message = 'Expecting 2 integers'
      Read (Buffer, *, Err=9999, End=9999) Graph, Start_HC_Species,
     >  End_HC_Species
      Return
C
 9998 IErr = 1
      Write (6,'(1X,A,4(/1X,A))')
     >  'RDG_HC_Trans_Prob: ERROR READING ',Graph,Message,'LAST LINE 
     >  READ :-',Buffer
      Write (7,'(1X,A,4(/1X,A))')
     >  'RDG_HC_Trans_Prob: ERROR READING ',Graph,Message,'LAST LINE 
     >  READ :-',Buffer
      RETURN
C
 9999 IErr = 1
      Write (6,'(1X,A,4(/1X,A))')
     >  'RDG_HC_Trans_Prob: ERROR READING ',Graph,Message,'LAST LINE 
     >  READ :-',Buffer
      Write (7,'(1X,A,4(/1X,A))')
     >  'RDG_HC_Trans_Prob: ERROR READING ',Graph,Message,'LAST LINE 
     >  READ :-',Buffer
      Return
      End
! ammod end.
