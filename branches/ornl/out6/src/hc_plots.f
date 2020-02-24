c     -*-Fortran-*-
c
c
      subroutine hc_plot_particle_track(iref,iopt)
c
c     Plots HC particle tracks 
c
!      Use ComHC ! Included to use hc_density, Number_HC_Species.
!      Use HC_Init_DIV_Data ! External data.
      Use HC_Init_Out_Data ! Include global data.
      Use Out_Utils ! Includes HC_Ident_Species function.
      use mod_outcom
      use mod_cgeom
      implicit none
      integer iref,iopt
c     include 'outcom'
c     include 'cgeom'
      integer lw,uw,iw
      integer ig
c
C
C-----------------------------------------------------------------------
C     TRACKS OF SELECTED Hydrocarbons (taken from option 181)
C-----------------------------------------------------------------------
C
       If (IOPT .lt. 0 .or. IOPT .gt. Max_Number_Walks) Then
          Write (Output_Unit_HC_Alert,*) 
     >      "Check number of particles to trace specified 
     >      in output file:",IOPT
          Write (6,*) 
     >      "Check number of particles to trace specified 
     >      in output file:",IOPT
          Write (0,*) 
     >      "Check number of particles to trace specified 
     >      in output file:",IOPT
          Return
       End If

c   
c     HC Particle track plots
c
      if (iref.eq.183) then 

        XLAB   = '   X/A  (M)'                                              
        YLAB   = '   Y/A  (M)'                                              
        NPLOTS = NPLOTS - 1                                                 
        LW     = 1                                                          
        UW     = 0                                                          
        IG     = 1                                                          
        NAME = '     TRACK'                                                 

        Do                                                                  

          Do IW = LW, Max_Number_Walks                                             
            If ( HC_Walks (IW,1) .ge. 10.0               
     >        *  RMAX) Then                         
	       Exit                                                         
	    End If                                                          
	    UW = IW                                                         
          End Do                                                            
            
          Write (REF,'(''HC TRACK'',I3,''  ('',I6,'' PTS)'')')              
     >      IG,UW-LW+1                                                      
          If (LW .lt. Max_Number_Walks .and.                                
     >      UW.GT.LW.AND.IG.LE.IOPT) THEN                                   
            NPLOTS = NPLOTS + 1                                             
            WRITE (IPLOT,9012) NPLOTS,REF                                   
            CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,             
     >        YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)              
            CALL SUPIMP ('SELECT')                                          
	    CALL GRTRAC ( HC_Walks (LW,1),               
     >         HC_Walks (LW,2),UW-LW+1,                  
     >        NAME,'LINE',1)                                                
            CALL FRAME                                                      
            IG = IG + 1                                                     
          ENDIF                                                             
	  LW = UW + 2                                                       

          If (LW .ge. Max_Number_Walks) Exit                                

        End Do                                                              
                                                                            
c                                                                           
c ----------------------------------------------------------------------    
c     Magnified hydrocarbon track (taken from option 182)                   
c ----------------------------------------------------------------------    
c                                                                           
                                                                            
      elseif (iref.eq.184) then                                             

        XLAB   = '   X/A  (M)'                                              
        YLAB   = '   Y/A  (M)'                                              
        NPLOTS = NPLOTS - 1                                                 
        LW     = 1                                                          
        UW     = 0                                                          
        IG     = 1                                                          
        NAME = '     TRACK'                                                 
                                                                            
                                                                            
        Do                                                                  
          Do IW = LW, Max_Number_Walks                                      
            IF ( HC_Walks (IW,1).GE.10.0                 
     >        * RMAX) Exit                          
            UW = IW                                                         
          End Do                                                            

          WRITE (REF,'(''HC TRACK'',I3,''  ('',I6,'' PTS)'')')              
     >      IG,UW-LW+1                                                      

          IF (LW .LT. Max_Number_Walks.AND.                                 
     >      UW.GT.LW.AND.IG.LE.IOPT) THEN                                   
            NPLOTS = NPLOTS + 1                                             
            WRITE (6,9012) NPLOTS,REF                                       
            CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,                         
     >        XXMIN-(RXP- HC_Walks (LW+1,1)),            
     >        XXMAX-(RXP- HC_Walks (LW+1,1)),            
     >        YYMIN-(ZXP- HC_Walks (LW+1,2)),            
     >        YYMAX-(ZXP- HC_Walks (LW+1,2)),            
     >        TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)                          
                                                                            
            CALL SUPIMP ('FULL')                                            
            CALL THICK (2)                                                  
            CALL GRTRAC ( HC_Walks (LW,1),               
     >         HC_Walks (LW,2),UW-LW+1,NAME,             
     >        'LINE',1)                                                     
            CALL THICK (1)                                                  
            CALL FRAME                                                      
            IG = IG + 1                                                     
          ENDIF                                                             
          LW = UW + 2                                                       

          If (LW .ge. Max_Number_Walks) Exit                                

	End Do                                                              

c
c
C-----------------------------------------------------------------------
C     PIN - TRACKS OF ALL RECORDED hydrocarbons (taken from option 651)
C-----------------------------------------------------------------------
C

      elseif (iref.eq.185) then 

        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NPLOTS = NPLOTS - 1
        LW     = 1
        UW     = 0
        IG     = 1
        NAME = '    TRACKS'
c
        WRITE (REF,'(''Hydrocarbon Tracks: '',I4)') IOpt 
	WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >      YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL SUPIMP ('SELECT')
C

        Do
          Do IW = LW, Max_Number_Walks
            If ( HC_Walks (IW,1).GE.10.0*
     >         RMAX) Then
	    	Exit
	    End If
            UW = IW
          End Do
          If (LW .lt. Max_Number_Walks
     >      .and. UW .gt. LW .and. IG .le. IOpt) Then
            CALL GRTRAC ( HC_Walks (LW,1),
     >         HC_Walks(LW,2),UW-LW+1,NAME,
     >        'LINE',1)
            IG = IG + 1
          ENDIF
          LW = UW + 2
          IF (LW .ge. Max_Number_Walks) Exit
	End Do

        CALL FRAME


      endif

c
c     Fortmat statements
c

 9012 FORMAT(1X,'PLOT',I3,4X,A)

      return
      end
c
c
c
      subroutine global_hc_plot_setup
! ammod begin.
      Use HC_Init_Out_Data ! Include global data.
! ammod end.
      use mod_outcom
      implicit none
c      include 'params'
c     include 'outcom'
! ammod begin.
      Real HC_NP, HC_NS, HC_NT
! ammod end.
      
! ammod begin.
c
c     jdemod - these quantities that are calculated based
c     on the hc_factors aren't useful since the HC code
c     doesn't even use hc_factor(0) and (-1) since the 
c     HC state array starts at 1. In DIVIMP these were for 
c     primary and secondary neutrals - a concept with 
c     limited applicability in the context of the HC code
c     especially since there is no separate recording of
c     "primary" and "secondary" neutrals.
c
c      IF (FACTA(-1).GT.0.0) Then
c	 HC_NP = 1.0 /  HC_Factor_A(-1)
c      End If
c      IF (FACTA(0) .GT.0.0) Then
c	 HC_NT = 1.0 /  HC_Factor_A(0)
c      End If
c      HC_NS = HC_NT - HC_NP
c      If (HC_NS .gt. 0.001) Then
c         HC_FP = HC_NP / HC_NS
c	 HC_FT = HC_NT / HC_NS
c      End If
	
c      WRITE (6,'(/,'' HCOUT: HC_NP,HC_NS,HC_NT,
c     >   HC_FP, HC_FT='',5F10.3,/)
c     >  ') HC_NP,HC_NS,HC_NT,
c     >   HC_FP, HC_FT
! ammod end.

c 
      return
      end
c
c
c
      subroutine global_hc_read_additional_data(iref,graph,iopt,ierr)
      use comhc
      use mod_outcom
      implicit none
      integer iref,iopt,ierr
      character*(*) graph
c     include 'outcom'


! ammod begin.
      ! Real in additional Hydrocarbon data.
      If (IRef .eq. 513 .or. IRef .eq. 514) Then
         Call RDG_HC_Trans_Prob (Graph,Start_HC_Species,
     >     End_HC_Species, IErr)
	 
	 If (IErr .ne. 0) Then
	    Write (6,*) 'Error reading HC transition probability
     >        data, IRef = ',IRef
            Write (Output_Unit_HC_Alert,*) 'Error reading HC transition probability
     >        data, IRef = ',IRef
            IErr = 1
	    Return
	 End If
      End If
! ammod end.

      return
      end
c
c
c
      subroutine load_hc_data_array(tmpplot,iselect,istate,itype,
     >                           ref,nizs,mfact,absfac,ierr)
      Use ComHC ! HC constants.
!      Use HC_Init_DIV_Data ! External data.
      Use HC_Init_Out_Data ! Included to re-set hc_state, hc_density, hc_output, Number_HC_Species.
      Use Out_Utils ! Includes HC_Ident_Species function.
c
c     LOAD_HC_DATA_ARRAY: 
c
c     This routine loads a standard ik,ir grid with the appropriate
c     HC related quantity specified by iselect, istate and itype
c
c     The current options are:
c 
c              28 = HC - Calculation of CD EMISSION (D/XB)
c                   1 - CD Efficiency (D/XB)
c                   2 - CD Emission (photons/m3) (jdemod)
c              29 = HC - HC State density
c                   istate = specific HC species 
c                          = sum over states for greater than maxstate   
c              30 = HC - HC State Ionization
c                   istate = specific HC species (ONLY CH So far)
c
      use mod_cgeom
      implicit none
c
c     include 'cgeom' 
c
      real tmpplot(maxnks,maxnrs)
      integer iselect,istate,nizs,ierr,itype      
      character*(*) ref 
c
      character lstring*10
      integer len,lenstr
      external lenstr
c
      real mfact,absfac
c
      integer istate_cd

c
c     CD Emission - local variables
c

! ammod begin.
      ! Declare variables for HC LOS calculation.
      Real :: Te
      real :: Xsec
      Real :: ne
      Real :: nhc
      Real :: Gaunt
      Real :: fCD
      Real :: deltaZ
      Real :: deltaE
      Real :: Imethane
      Real :: XCD
      Real :: Aplasma
      Real :: DXBHCmax
      Real :: DXBHCavg
! ammod end.      
c
      integer ik,ir


c
c     Apply scaling factor
c

      IF (ABSFAC.GT.0.0) MFACT = MFACT * ABSFAC

c
c     CD D/XB 
c
      if (iselect.eq.28) then 
c
c jdemod - NOTE: 
c
c          The calculation of D/XB here is not useful in any practical situation - 
c          in the paper the analysis was confined to a predictable near target
c          volume with constant plasma conditions and a known methane influx - this 
c          allowed the D/XB to be estimated experimentally and calculated by WBC.
c          
c          The equivalent quantity calculated on a cell by cell basis on the grid 
c          would require that the CD4 influx into each cell be recorded and then used to 
c          calculate a D/XB for each cell based on the local emission. It might be possible
c          to use the CD4 density in this context - either directly or by assuming a reasonable
c          velocity to obtain a flux estimate. In any case - the numerical values below are just
c          incorrect for anything except trying to compare the near target volume emission.   
c
c        ISTATE=2 should be used for CD emission comparisons and calculations. 
c
         if (istate.eq.1) then 
            ! Calculates (D/XB)CD-430nm using "Investigation of
            ! Carbon Chemical Erosion with Increasing Plasma
	    ! Flux & Density", D.G.Whyte, G.Tynan, R.Doerner,
	    ! J.N.Brooks (NucFus, 2001), eqs. 10/11.
            ! Setup values for equation.
            !
            !
            ! NOTE: CD data is currently state 4 in the stored 
            !       HC data. 
            !
            istate_cd = 4

            write(0,*) 
     >       'WARNING: LOAD_HC_DATA_ARRAY: OPT 28:'//
     >       ' SOME VALUES HARD CODED FOR HC EMISSION'//
     >       ' CALCULATIONS'
            write(6,*) 
     >       'WARNING: LOAD_HC_DATA_ARRAY: OPT 28:'//
     >       ' SOME VALUES HARD CODED FOR HC EMISSION'//
     >       ' CALCULATIONS'

            fCD = 5.3E-3 ! Optical f-value.
	    Gaunt = 1.0 ! Gaunt factor.
	    deltaE = 2.9 ! Photon energy [eV].
	    deltaZ = 0.1 ! Region of interest [m].
	 
            Imethane = 1E11 ! CH4 injection rate [1/s].
            Aplasma = 2.1 ! Plasma cross-sectional area for DIII-D [m].
	 
            DXBHCavg = 0.0
            DXBHCmax = 0.0        

            Do ir = 1,nrs

               Do ik = 1, nks(ir)
                  ne = knbs (ik,ir)
                  Te = ktebs (ik,ir)

                  nhc =  HC_Density (ik,ir,istate_cd)

c
c                 jdemod - possible bug - in the paper the exponential 
c                          factor is exp(-deltaE/Te) not +deltaE
c
c                 XCD=1.6E-11*fCD*Gaunt/(deltaE*SQRT(Te))*EXP(deltaE/Te)
c
                 XCD=1.6E-11*fCD*Gaunt/(deltaE*SQRT(Te))*EXP(-deltaE/Te)


                  tmpplot(ik,ir) = Imethane/(Aplasma*deltaZ*ne*nhc*XCD)
	       
                  DXBHCavg = DXBHCavg + tmpplot(ik,ir)
                  If (tmpplot(ik,ir) .gt. DXBHCmax) Then
                     DXBHCmax = tmpplot(ik,ir)
                  End If

               End Do

            End Do
	 
            ! Print diagnostics.  
            Write (Output_Unit_HC_Data,*) "Maximum D/XB:",DXBHCmax
            Write (Output_Unit_HC_Data,*) "Average D/XB:",
     >                                           DXBHCavg/(nks*nrs)
c
            write (Output_Unit_HC_Data,*) "Number of grid cells?",
     >                                              nks*nrs 
c
c        jdemod  
c
c        CD emission calculation
c
         elseif (istate.eq.2) then 

            ! Calculates XCD-430nm using "Investigation of
            ! Carbon Chemical Erosion with Increasing Plasma
	    ! Flux & Density", D.G.Whyte, G.Tynan, R.Doerner,
	    ! J.N.Brooks (NucFus, 2001), eqs. 10/11.
            ! Setup values for equation.

            !
            ! NOTE: CD data is currently state 4 in the stored 
            !       HC data. 
            !

            fCD = 5.3e-3 ! Optical f-value.
	    Gaunt = 1.0 ! Gaunt factor.
	    deltaE = 2.9 ! Photon energy [eV].
            istate_cd = 4 ! location of CD data - replace with global later
	 
            Do ir = 1,nrs

               Do ik = 1, nks(ir)

                  ne = knbs (ik,ir)
                  Te = ktebs (ik,ir)

                  nhc =  HC_Density (ik,ir,istate_cd)

                  XCD = 1.6e-11*fCD*Gaunt/(deltaE*SQRT(Te))
     >                     * EXP(-deltaE/Te)

                  tmpplot(ik,ir) = ne * nhc * xcd
	       
               End Do

            End Do
	 
            ! Print diagnostics.  
            Write (Output_Unit_HC_Data,*) "Maximum D/XB:",DXBHCmax
            Write (Output_Unit_HC_Data,*) "Average D/XB:",DXBHCavg/
     >                                              (nks*nrs)
c
            write (Output_Unit_HC_Data,*) "Number of grid cells?",
     >                                              nks*nrs 
c
c        CD emission calculation - calculated using the profiles from 
c 
c        http://www.eirene.de/eigen/Sample_ex/ax/ax.html
c
c        U. Fantz J. Nucl. Mat. 337-339 (2005) p1087 - 1091, and: U. Fantz priv. comm. 2005
c
c        Fit <ex> = 1.51e-14 * exp(-3.19/Te) / Te**0.29 * 1.e6   [cm3/s]
c                   remove 1e6 factor for m3/s
c
         elseif (istate.eq.3) then 

            !
            ! NOTE: CD data is currently state 4 in the stored 
            !       HC data. 
            !
	 
            Do ir = 1,nrs

               Do ik = 1, nks(ir)

                  ne = knbs (ik,ir)
                  Te = ktebs (ik,ir)
                  xsec = 1.51e-14 * exp(-3.19/Te) / Te**0.29

                  nhc =  HC_Density (ik,ir,istate_cd)

                  tmpplot(ik,ir) = ne * nhc * xsec
	       
               End Do

            End Do
	 
         endif
c
c        jdemod
c
      elseif (iselect.eq.29) then 
c
c     HC state density
c     
c 
         If (Istate .lt. 0 .or. Istate .gt. 
     >      Number_HC_Species) Then
            write(0,*) 'ERROR: NON-Existant HC Species'//
     >                    ' specified for plot',istate
            write(6,*) 'ERROR: NON-Existant HC Species'//
     >                    ' specified for plot',istate
            Return 
       	 End If
c
         lstring = HC_Ident_Species (istate)
         LEN    = LENSTR (lstring)
         REF    = 'HC Density of ' // lString (1:LEN)

c
c        Assign data array
c
         tmpplot =  HC_Density(:,:,istate) * mfact
c

      elseif (iselect.eq.30) then 

c
c        Only CH Ionization is recorded - presumably entry 0 is the
c        sum for all CH neutrals. 
c
         REF = 'TOTAL CH IONISATION '
c
         tmpplot =  HC_TIZS_CH(:,:,0) * mfact
c
      endif
c
      return
      end
c
c
c
      subroutine hc_set_ylab(iselect,istate,itype,nizs,ylab)
      Use ComHC ! HC constants.
!      Use HC_Init_DIV_Data ! External data.
      Use HC_Init_Out_Data ! Included to re-set hc_state, hc_density, hc_output, Number_HC_Species.
      Use Out_Utils ! Includes HC_Ident_Species function.
c
c     SET Y Axis labels for HC plots
c
c     HC Related quantities ISELECT values
c
c     HC data - 28 = HC CD Emission
c               29 = HC state density
c               30 = HC state ionization
c               33 = SUBGRID HC state density
c               35 = SUBGRID HC emission
c

      implicit none
      integer iselect,istate,itype,nizs
      character*(*) ylab
      Character lstring * 10 ! Character string of HC.
      integer len,lenstr
      external lenstr
c
c
c     CD D/XB and emission
c
      if (iselect.eq.28) then 
c
         if (istate.eq.1) then 
            write(YLAB,'(''CD  CD D/XB'')')
         elseif (istate.eq.2) then 
            write(YLAB,'(''CD  CD-430 EMISSION'')')
            len = lenstr(ylab)
            if (itype.eq.0) then 
               ylab=ylab//' (/M3/S)'
            else
               ylab=ylab//' (/M2/S)'
            endif

         endif
c
c     CD Emission
c
      elseif (iselect.eq.35) then 
c
         if (istate.eq.4) then 
            write(YLAB,'(''CD  CD-430 EMISSION'')')
            len = lenstr(ylab)
            if (itype.eq.0) then 
               ylab=ylab//' (/M3/S)'
            else
               ylab=ylab//' (/M2/S)'
            endif

         endif
c
c     HC Species densities
c
c     29 - HC species density
c     33 - SUBGRID HC species density
c
      elseif (iselect.eq.29.or.iselect.eq.33) then 

         If (Istate .lt. 0 .or. Istate .gt. 
     >      Number_HC_Species) Then
            Return 
       	 End If

        lstring = HC_Ident_Species (istate)
        LEN    = LENSTR (lstring)

         write(YLAB,'(''HC DENSITY OF '',a,
     >                '' (M^-3)'')') lstring(1:len)

c
c     HC Species Ionization rates
c 
      elseif (iselect.eq.30) then 
c
         write(YLAB,'(''TOTAL CH IONIZATION (M^-3S^-1'')')
c
      endif


      return
      end
c
c
c
      subroutine hc_set_blab(iselect,istate,itype,nizs,blab)
      Use ComHC ! HC constants.
!      Use HC_Init_DIV_Data ! External data.
      Use HC_Init_Out_Data ! Included to re-set hc_state, hc_density, hc_output, Number_HC_Species.
      Use Out_Utils ! Includes HC_Ident_Species function.

      implicit none
      integer iselect,istate,itype,nizs
      character*(*) blab
      Character lstring * 10 ! Character string of HC.
      integer len,lenstr
      external lenstr
c
c     SET TABLE Labels for HC quantities
c
c
c     HC Related quantities ISELECT values
c
c     HC data - 28 = HC CD Emission
c               29 = HC state density
c               30 = HC state ionization
c               33 = SUBGRID HC state density
c               35 = SUBGRID HC emission
c
c
c     CD D/XB (possibly emission at a later stage)
c
      if (iselect.eq.28) then 
c
         if (istate.eq.1) then 
            write(BLAB,'(''CD D/XB'')')
         elseif (istate.eq.2) then 
            write(BLAB,'(''CD-430 EMISSION'')')
         endif
c
      elseif (iselect.eq.35) then 
c
         if (istate.eq.4) then 
            write(BLAB,'(''CD-430 EMISSION'')')
         endif
c
c
c     HC Species densities
c
      elseif (iselect.eq.29.or.iselect.eq.33) then 

         If (Istate .lt. 0 .or. Istate .gt. 
     >      Number_HC_Species) Then
            Return 
       	 End If
c
        lstring = HC_Ident_Species (istate)
        LEN    = LENSTR (lstring)

         write(BLAB,'(''HC DENSITY OF '',a,
     >                '' (M^-3)'')') lstring(1:len)

c
c     HC Species Ionization rates
c 
      elseif (iselect.eq.30) then 
c
         write(BLAB,'(''TOTAL CH IONIZATION (M^-3S^-1'')')
c
      endif


      return
      end
c
c
c
      subroutine global_hc_plot_init(crmb,crmi)
      

! ammod begin.
      ! Required to load hydrocarbon data for plots.
      Use ComHC ! Data for loading HC libraries.
      Use HC_Init_Lib_Data ! Function for getting HC charge.
      Use HC_LdDta ! Access hydrocarbon storage data.
      !Use HC_Init_DIV_Data ! External data.
      Use HC_Init_Out_Data ! Include global data.

      implicit none
c
      real :: crmb,crmi

c
! ammod begin.
      ! Call setup routines for hydrocarbon change of state.
      ! jdemod - pass masses to HC initialization routines

      !
      ! jdemod :: 
      ! This needs to be fixed by saving input_HC_H_mass in the
      ! raw file 
      !
      input_HC_H_mass = crmb
      Call Initialize_HC_Data (crmi)
!      Call Initialize_HC_Data (crmb,crmi)
      Call Load_HC_Data ()
! ammod end.


      ! Hydrocarbon data initialization.
      !Call Initialize_Global_Prop_Data
      !Call Initialize_Global_Geom_Data



      return
      end
