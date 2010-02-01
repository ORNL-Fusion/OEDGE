! -*-Fortran-*-
! Hydrocarbon Put/Get.f
! Fixed Format Put/Get Routines 
!
! Adam McLean under Prof. Peter Stangeby
! Research Assistant: David Elder
! July, 2002
!
! Routines written to interact with DIVIMP data arrays in order
! to facilitate use of free format in the main code, but use 
! common block data formated in fixed format transparently, since
! using Include declarations makes the continuation character
! construct a difficulty for mixed source files.

      Module HC_Put

        ! jdemod - removed use of comhc module - duplicates information in the 
        ! included params files
	!Use ComHC ! Hydrocarbon-related variables set in DIVIMP input file.

	! Every good Fortran program has...
        Implicit None

        integer,private:: Output_unit_hc_alert


      Contains

      subroutine hc_put_init_output(outunit)
      implicit none
      integer :: outunit
      !
      ! This routine initializes the output unit for warnings from these 
      ! routines
      !
         output_unit_hc_alert = outunit

      end subroutine hc_put_init_output



     
!--------------------------------
! PUT Routines
!--------------------------------

	Subroutine pcflrex (Value)
		! Routine to store (or put) value into cflrex array in COMMV.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'commv'
		
		! Declare input variables.
		Logical, Intent (In) :: Value
		
		! Store value.
		cflrex = Value
			
	End Subroutine pcflrex

	Subroutine pcflrin (Value)
		! Routine to store (or put) value into cflrin array in COMMV.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'commv'
		
		! Declare input variables.
		Logical, Intent (In) :: Value
		
		! Store value.
		cflrin = Value
	
	End Subroutine pcflrin

	Subroutine pcflrxa (Value)
		! Routine to store (or put) value into cflrxa array in COMMV.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'commv'
		
		! Declare input variables.
		Logical, Intent (In) :: Value
		
		! Store value.
		cflrxa = Value
	
	End Subroutine pcflrxa

	Subroutine pcistizs (Index, Value)
		! Routine to store (or put) value into cistizs array in CNEUT.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value

		! Store value.
		cistizs (Index) = Value
	
	End Subroutine pcistizs

	Subroutine pclll (Index, Value)
		! Routine to store (or put) value into clll array in COMMV.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'commv'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value

		! Store value.
		clll (Index) = Value
	
	End Subroutine pclll

	Subroutine pcmmm (Index, Value)
		! Routine to store (or put) value into cmmm array in COMMV.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'commv'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value

		! Store value.
		cmmm (Index) = Value
	
	End Subroutine pcmmm

	Subroutine pcnnn (Index, Value)
		! Routine to store (or put) value into cnnn array in COMMV.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'commv'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value

		! Store value.
		cnnn (Index) = Value
	
	End Subroutine pcnnn

	Subroutine pcstepl (Value)
		! Routine to store (or put) value into cstepl variable in COMTOR.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'comtor'
		
		! Declare input variables.
		Real, Intent (In) :: Value

		! Store value.
		cstepl = Value
	
	End Subroutine pcstepl

	Subroutine pdebugl (Value)
		! Routine to store (or put) value into debugl variable in COMTOR.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'comtor'
		
		! Declare input variables.
		Logical, Intent (In) :: Value

		! Store value.
		debugl = Value
	
	End Subroutine pdebugl

	Subroutine pddts (Index_1, Index_2, Index_3, Value)
		! Routine to store (or put) value into ddts array in DYNAM1.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'dynam1'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
		Real, Intent (In) :: Value
		
		! Store value.
		ddts (Index_1, Index_2, Index_3) = Value
	
	End Subroutine pddts
	
	Subroutine peprods (Index, Value)
		! Routine to store (or put) value into eprods array in CNEUT.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value
		
		! Store value.
		eprods (Index) = Value
	
	End Subroutine peprods

	Subroutine peranv (Seed, Random_Numbers_Used)
		! Put values into external memory block RANV.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'crand'

		! Declare input variables.
		Double Precision, Intent (In) :: Seed
		Integer, Intent (In) :: Random_Numbers_Used
		
		Call SURAND (Seed,Random_Numbers_Used,Ranv)

	End Subroutine peranv
	
	Subroutine peranva (Seed, Random_Numbers_Used)
		! Put values into external memory block RANV.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'cneut'
		
		! Declare input variables.
		Double Precision, Intent (In) :: Seed
		Integer, Intent (In) :: Random_Numbers_Used

		Call SURAND (Seed,Random_Numbers_Used,Ranva)

	End Subroutine peranva
	
	Subroutine peranvb (Seed, Random_Numbers_Used)
		! Put values into external memory block RANVB.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'cneut'

		! Declare input variables.
		Double Precision, Intent (In) :: Seed
		Integer, Intent (In) :: Random_Numbers_Used
		
		Call SURAND (Seed,Random_Numbers_Used,Ranvb)

	End Subroutine peranvb

	Subroutine peranvc (Seed, Random_Numbers_Used)
		! Put values into external memory block RANVC.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'cneut'
		
		! Declare input variables.
		Double Precision, Intent (In) :: Seed
		Integer, Intent (In) :: Random_Numbers_Used

		Call SURAND (Seed,Random_Numbers_Used,Ranvc)

	End Subroutine peranvc

	Subroutine pidatizs (Index_1, Index_2, Value)
		! Routine to store (or put) value into idatizs array in CNEUT.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Value
		
		! Store value.
		idatizs (Index_1, Index_2) = Value
	
	End Subroutine pidatizs
	
	Subroutine pidprods (Index, Value)
		! Routine to store (or put) value into idprods array in CNEUT.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Integer, Intent (In) :: Value
		
		! Store value.
		idprods (Index) = Value
	
	End Subroutine pidprods

	Subroutine plaunchdat (Index_1, Index_2, Value)
		! Routine to store (or put) value into launchdat array in COMTOR.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'comtor'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Real, Intent (In) :: Value
		
		! Store value.
		launchdat (Index_1, Index_2) = Value
	
	End Subroutine plaunchdat

	Subroutine pkoutds (Index_1, Index_2, Value)
		! Routine to store (or put) value into koutds array in CGEOM.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Value
		
		! Store value.
		koutds (Index_1, Index_2) = Value
	
	End Subroutine pkoutds
	
	Subroutine pkatizs (Index, Value)
		! Routine to store (or put) value into katizs array in CNEUT.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value
		
		! Index should be within 1 - MAXIMP
		If (Index .lt. 1 .or. Index .gt. MAXIMP) Then
			Write (Output_Unit_HC_Alert,*) 
     >                    "Error in subroutine pkatizs:
     >                    Return value outside range of KATIZS."
		End If

		! Store value.
		katizs (Index) = Value
	
	End Subroutine pkatizs

	Subroutine plfps (Index_1, Index_2, Index_3, Value)
		! Routine to store (or put) value into lfps array in CLOCAL.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'clocal'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
		Real, Intent (In) :: Value
		
		! Store value.
		lfps (Index_1, Index_2, Index_3) = Value
	
	End Subroutine plfps

	Subroutine plfss (Index_1, Index_2, Index_3, Value)
		! Routine to store (or put) value into lfss array in CLOCAL.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'clocal'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
		Real, Intent (In) :: Value
		
		! Store value.
		lfss (Index_1, Index_2, Index_3) = Value
	
	End Subroutine plfss

	Subroutine plfts (Index_1, Index_2, Index_3, Value)
		! Routine to store (or put) value into lfts array in CLOCAL.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'clocal'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
		Real, Intent (In) :: Value
		
		! Store value.
		lfts (Index_1, Index_2, Index_3) = Value
	
	End Subroutine plfts

	Subroutine plllfps (Index_1, Index_2, Index_3, Value)
		! Routine to store (or put) value into lllfps array in CLOCAL.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'clocal'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
		Real, Intent (In) :: Value
		
		! Store value.
		lllfps (Index_1, Index_2, Index_3) = Value
	
	End Subroutine plllfps

	Subroutine pltolds (Index_1, Index_2, Index_3, Value)
		! Routine to store (or put) value into ltolds array in CLOCAL.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'clocal'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
		Real, Intent (In) :: Value
		
		! Store value.
		ltolds (Index_1, Index_2, Index_3) = Value
	
	End Subroutine pltolds

	Subroutine pranv (Index, Value)
		! Routine to store (or put) value into ranv array in CRAND.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'crand'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value

		! Store value.
		ranv (Index) = Value
	
	End Subroutine pranv

	Subroutine pranva (Index, Value)
		! Routine to store (or put) value into ranva array in CNEUT.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value

		! Store value.
		ranva (Index) = Value
	
	End Subroutine pranva

	Subroutine pranvb (Index, Value)
		! Routine to store (or put) value into ranvb array in CNEUT.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value

		! Store value.
		ranvb (Index) = Value
	
	End Subroutine pranvb

	Subroutine pranvc (Index, Value)
		! Routine to store (or put) value into ranvc array in CNEUT.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value

		! Store value.
		ranvc (Index) = Value
	
	End Subroutine pranvc
	
	Subroutine psatizs (Index, Value)
		! Routine to store (or put) value into satizs array in CNEUT.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value
		
		! Store value.
		satizs (Index) = Value
	
	End Subroutine psatizs

	Subroutine psputys (Index, Value)
		! Routine to store (or put) value into sputys array in CNEUT.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value
		
		! Index should be within 1 - MAXIMP
		If (Index .lt. 1 .or. Index .gt. MAXIMP) Then
			Write (Output_Unit_HC_Alert,*) 
     >                    "Error in subroutine psputys:
     >                    Return value outside range of SPUTYS."	
		End If

		! Store value.
		sputys (Index) = Value
	
	End Subroutine psputys

	Subroutine psnews (Index, Value)
		! Routine to store (or put) value into snews array in CNEUT.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value

		! Store value.
		snews (Index) = Value
	
	End Subroutine psnews

	Subroutine ptemtizs (Index, Value)
		! Routine to store (or put) value into temtizs array from CNEUT.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value		
				
		! Returned value.
		temtizs (Index) = Value
		
	End Subroutine ptemtizs

	Subroutine ptravel_locations (Index1, Index2, Value)
		! Routine to store (or put) value into travel_locations array in CNEUT.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index1
		Integer, Intent (In) :: Index2
		Logical, Intent (In) :: Value
		
		! Store value.
		travel_locations (Index1,Index2) = Value
	
	End Subroutine ptravel_locations

	Subroutine pvins (Index, Value)
		! Routine to store (or put) value into vins array in CNEUT.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value
		
		! Index should be from 1 to MAXIMP.
		If (Index .lt. 1 .or. Index .gt. MAXIMP) Then
			Write (Output_Unit_HC_Alert,*) 
     >                    "Error in subroutine pvins:
     >                    Return value outside range of VINS: ",Index
		End If

		! Store value.
		vins (Index) = Value
	
	End Subroutine pvins

	Subroutine pwalks (Index_1, Index_2, Value)
		! Routine to store (or put) value into walks array in DYNAM4.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'dynam4'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Real, Intent (In) :: Value
		
		! Store value.
		walks (Index_1, Index_2) = Value
	
	End Subroutine pwalks

	Subroutine pxatizs (Index, Value)
		! Routine to store (or put) value into xatizs array in CNEUT.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value
		
		! Index should be within 1 - MAXIMP
		If (Index .lt. 1 .or. Index .gt. MAXIMP) Then
			Write (Output_Unit_HC_Alert,*) 
     >                    "Error in subroutine pxatizs: "//
     >               " Return value outside range of XPATIZS:"//
     >               " Program STOPPING"
			Write (0,*) 
     >                    "Error in subroutine pxatizs: "//
     >               " Return value outside range of XATIZS:"//
     >               " Program STOPPING"
                     stop
		End If

		! Store value.
		xatizs (Index) = Value
	
	End Subroutine pxatizs
	
	Subroutine pxprods (Index, Value)
		! Routine to store (or put) value into xprods array in CNEUT.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value

                ! jdemod 
                ! Check the array bounds to make sure it is in range
                ! 
		! Index should be within 1 - MAXIMP
		If (Index .lt. 1 .or. Index .gt. MAXIMP) Then
			Write (Output_Unit_HC_Alert,*) 
     >                    "Error in subroutine pxprods: "//
     >               " Return value outside range of XPRODS:"//
     >               " Program STOPPING"
			Write (0,*) 
     >                    "Error in subroutine pxprods: "//
     >               " Return value outside range of XPRODS:"//
     >               " Program STOPPING"
                     stop
		End If

		
		! Store value.
		xprods (Index) = Value
	
	End Subroutine pxprods
	
	Subroutine pyatizs (Index, Value)
		! Routine to store (or put) value into yatizs array in CNEUT.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value
		
		! Index should be within 1 - MAXIMP
		If (Index .lt. 1 .or. Index .gt. MAXIMP) Then
			Write (Output_Unit_HC_Alert,*) 
     >                    "Error in subroutine pyatizs:
     >                    Return value outside range of YATIZS."	
		End If

		! Store value.
		yatizs (Index) = Value
	
	End Subroutine pyatizs

	Subroutine pyprods (Index, Value)
		! Routine to store (or put) value into yprods array in CNEUT.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value
		
		! Store value.
		yprods (Index) = Value
	
	End Subroutine pyprods

! jdemod - remove these subroutines - they don't do anything - I assume they were a first attempt at a 
!          more general put/get set of routines. 
!
!	Subroutine put_divimp_real (Option, Index_1, Index_2, Index_3,
!     > Value)
!		! Every good Fortran routine has...
!		Implicit None
!		
!		! Declare input variables
!		Character (Len=*) :: Option
!		Integer, Intent (In) :: Index_1
!		Integer, Intent (In) :: Index_2
!		Integer, Intent (In) :: Index_3
!		Integer, Intent (Out) :: Value
!		
!		! Required include files
!		Include 'params'
!		
!		
!		! Assignment of array value.
!	
!	End Subroutine put_divimp_real
!	
!	Subroutine put_divimp_integer (Option, Index_1, Index_2,
!     > Index_3, Value)
!		! Every good Fortran routine has...
!		Implicit None
!		
!		! Declare input variables
!		Character (Len=*) :: Option
!		Integer, Intent (In) :: Index_1
!		Integer, Intent (In) :: Index_2
!		Integer, Intent (In) :: Index_3
!		Integer, Intent (Out) :: Value
!		
!		! Required include files
!		Include 'params'
!		
!		
!		! Assignment of array value.
!	
!	End Subroutine put_divimp_integer
	
!--------------------------------
! PUT ADD Routines
!--------------------------------

	Subroutine pachemden (Index_1, Index_2, Value)
		! Routine to store (or put) value by adding it to chemden array in DYNAM3.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'dynam3'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Real, Intent (In) :: Value
		
		! Store value.
		chemden (Index_1, Index_2) = chemden (Index_1, Index_2) + Value
	
	End Subroutine pachemden

	Subroutine pachemizs (Index_1, Index_2, Value)
		! Routine to store (or put) value by adding it to chemizs array in DYNAM3.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'dynam3'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Real, Intent (In) :: Value
		
		! Store value.
		chemizs (Index_1, Index_2) = chemizs (Index_1, Index_2) + Value
	
	End Subroutine pachemizs

	Subroutine pacicuts (Index, Value)
		! Routine to store (or put) value by adding it to cicuts array in COMMV.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'commv'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value
		
		! Store value.
		cicuts (Index) = cicuts (Index) + Value
	
	End Subroutine pacicuts

	Subroutine pacieizs (Index, Value)
		! Routine to store (or put) value by adding it to cieizs array in COMMV.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'commv'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value
		
		! Store value.
		cieizs (Index) = cieizs (Index) + Value
	
	End Subroutine pacieizs

	Subroutine pacitizs (Index, Value)
		! Routine to store (or put) value by adding it to citizs array in COMMV.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'commv'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value
		
		! Store value.
		citizs (Index) = citizs (Index) + Value
	
	End Subroutine pacitizs

	Subroutine pacleakn (Index1, Index2, Value)
		! Routine to store (or put) value by adding it to cleakn array in COMTOR.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'comtor'
		
		! Declare input variables.
		Integer, Intent (In) :: Index1
		Integer, Intent (In) :: Index2
		Real, Intent (In) :: Value
		
		! Store value.
		cleakn (Index1,Index2) = cleakn (Index1,Index2) + 
     >            Value

	End Subroutine pacleakn
	
	Subroutine pacleakpos (Index1, Index2, Value)
		! Routine to store (or put) value by adding it to cleakpos array in COMTOR.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'comtor'
		
		! Declare input variables.
		Integer, Intent (In) :: Index1
		Integer, Intent (In) :: Index2
		Real, Intent (In) :: Value

!
!               jdemod - cleakpos isn't really used at the moment - however, it is 
!                        supposed to contain the starting R and Z ionization locations 
!                        for each particle entering the confined plasma. As
!                        such the values are only assigned not summed. 
!     
                cleakpos(index1,index2) = Value
                
		! Store value.
!		cleakpos (Index1,Index2) = cleakpos (Index1,Index2) + 
!     >            Value
                


	End Subroutine pacleakpos
	
	Subroutine paclll (Index, Value)
		! Routine to store (or put) value by adding it to clll array in COMMV.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'commv'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value
		
		! Store value.
		clll (Index) = clll (Index) + Value
	
	End Subroutine paclll

	Subroutine pacmmm (Index, Value)
		! Routine to store (or put) value by adding it to cmmm array in COMMV.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'commv'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value
		
		! Store value.
		cmmm (Index) = cmmm (Index) + Value
	
	End Subroutine pacmmm

	Subroutine pacnnn (Index, Value)
		! Routine to store (or put) value by adding it to cnnn array in COMMV.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'commv'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value
		
		! Store value.
		cnnn (Index) = cnnn (Index) + Value
	
	End Subroutine pacnnn

	Subroutine pacrtrcs (Index, Value)
		! Routine to store (or put) value by adding it to crtrcs array in COMMV.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'commv'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value
		
		! Store value.
		crtrcs (Index) = crtrcs (Index) + Value
	
	End Subroutine pacrtrcs

	Subroutine paddlims (Index_1, Index_2, Index_3, Value)
		! Routine to store (or put) value by adding it into elims array in DYNAM1.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'dynam1'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
		Double Precision, Intent (In) :: Value
		
		! Store value.
		ddlims (Index_1, Index_2, Index_3) = ddlims 
     > (Index_1, Index_2, Index_3) + Value
	
	End Subroutine paddlims

	Subroutine paddts (Index_1, Index_2, Index_3, Value)
		! Routine to store (or put) value by adding it to ddts array in DYNAM1.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'dynam1'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
		Double Precision, Intent (In) :: Value
		
		! Store value.
		ddts (Index_1,Index_2,Index_3) = ddts 
     > (Index_1,Index_2,Index_3) + Value
	
	End Subroutine paddts

	Subroutine paddvoid (Index, Value)
		! Routine to store (or put) value by adding it to ddvoid array in DYNAM1.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'dynam1'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Double Precision, Intent (In) :: Value
		
		! Store value.
		ddvoid (Index) = ddvoid (Index) + Value
	
	End Subroutine paddvoid

	Subroutine padeps (Index_1, Index_2, Value)
		! Routine to store (or put) value by adding it into deps array in DYNAM3.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'dynam3'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Real, Intent (In) :: Value
		
		! Store value.
		deps (Index_1, Index_2) = deps (Index_1,
     > Index_2) + Value
	
	End Subroutine padeps

	Subroutine padiff (Index_1, Index_2, Index_3, Value)
		! Routine to store (or put) value by adding it into diff array in REISER.
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'reiser_com'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
		Real, Intent (In) :: Value
		
		! Store value.
		diff (Index_1, Index_2, Index_3) = diff 
     > (Index_1, Index_2, Index_3) + Value
	
	End Subroutine padiff

	Subroutine paelims (Index_1, Index_2, Index_3, Value)
		! Routine to store (or put) value by adding it into elims array in DYNAM3.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'dynam3'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
		Real, Intent (In) :: Value
		
		! Store value.
		elims (Index_1, Index_2, Index_3) = elims 
     > (Index_1, Index_2, Index_3) + Value
	
	End Subroutine paelims

	Subroutine pafcell (Index_1, Index_2, Index_3, Value)
		! Routine to store (or put) value by adding it into fcell array in REISER.
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'reiser_com'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
		Real, Intent (In) :: Value
		
		! Store value.
		fcell (Index_1, Index_2, Index_3) = fcell 
     > (Index_1, Index_2, Index_3) + Value
	
	End Subroutine pafcell

	Subroutine paffi (Index_1, Index_2, Index_3, Value)
		! Routine to store (or put) value by adding it into ffi array in REISER.
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'reiser_com'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
		Real, Intent (In) :: Value
		
		! Store value.
		ffi (Index_1, Index_2, Index_3) = ffi 
     > (Index_1, Index_2, Index_3) + Value
	
	End Subroutine paffi

	Subroutine pafthi (Index_1, Index_2, Index_3, Value)
		! Routine to store (or put) value by adding it into fthi array in REISER.
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'reiser_com'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
		Real, Intent (In) :: Value
		
		! Store value.
		fthi (Index_1, Index_2, Index_3) = fthi 
     > (Index_1, Index_2, Index_3) + Value
	
	End Subroutine pafthi

	Subroutine pafvbg (Index_1, Index_2, Index_3, Value)
		! Routine to store (or put) value by adding it into fvbg array in REISER.
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'reiser_com'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
		Real, Intent (In) :: Value
		
		! Store value.
		fvbg (Index_1, Index_2, Index_3) = fvbg 
     > (Index_1, Index_2, Index_3) + Value
	
	End Subroutine pafvbg

	Subroutine paionizdat (Index_1, Index_2, Index_3, Index_4, 
     > Index_5, Value)
		! Routine to store (or put) value by adding it into ionizdat array in COMTOR2.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'comtor'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
		Integer, Intent (In) :: Index_4
		Integer, Intent (In) :: Index_5
		Real, Intent (In) :: Value
		
		! Store value.
		ionizdat (Index_1, Index_2, Index_3, Index_4, Index_5)
     > = ionizdat (Index_1, Index_2, Index_3, Index_4, Index_5) + Value
	
	End Subroutine paionizdat

	Subroutine paionvelavg (Index_1, Index_2, Index_3, Value)
		! Routine to store (or put) value by adding it into velavg array in REISER.
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'reiser_com'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
		Real, Intent (In) :: Value
		
		! Store value.
		velavg (Index_1, Index_2, Index_3) = velavg
     > (Index_1, Index_2, Index_3) + Value
	
	End Subroutine paionvelavg

	Subroutine palims (Index_1, Index_2, Index_3, Index_4, Value)
		! Routine to store (or put) value by adding it into lims array in DYNAM4.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'dynam4'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
		Integer, Intent (In) :: Index_4
		Real, Intent (In) :: Value
		
		! Store value.
		lims (Index_1, Index_2, Index_3, Index_4) = lims
     > (Index_1, Index_2, Index_3, Index_4) + Value
	
	End Subroutine palims

	Subroutine pamtcinf (Index_1, Index_2, Value)
		! Routine to store (or put) value by adding it to mtcinf array in CNEUT.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Real, Intent (In) :: Value
		
		! Store value.
		mtcinf (Index_1, Index_2) = mtcinf (Index_1,
     > Index_2) + Value
	
	End Subroutine pamtcinf

	Subroutine pamtctotcnt (Index_1, Index_2, Value)
		! Routine to store (or put) value by adding it to mtctotcnt array in CNEUT.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Real, Intent (In) :: Value
		
		! Store value.
		mtctotcnt (Index_1, Index_2) = mtctotcnt (Index_1,
     > Index_2) + Value
	
	End Subroutine pamtctotcnt

	Subroutine pancore (Index_1, Index_2, Value)
		! Routine to store (or put) value by adding it to ncore array in DYNAM3.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'dynam3'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Real, Intent (In) :: Value
		
		! Store value.
		ncore (Index_1, Index_2) = ncore (Index_1,
     > Index_2) + Value
	
	End Subroutine pancore

	Subroutine pandivert (Index_1, Index_2, Value)
		! Routine to store (or put) value by adding it to ndivert array in DYNAM3.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'dynam3'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Real, Intent (In) :: Value
		
		! Store value.
		ndivert (Index_1, Index_2) = ndivert (Index_1,
     > Index_2) + Value
	
	End Subroutine pandivert

	Subroutine panedge (Index_1, Index_2, Value)
		! Routine to store (or put) value by adding it to nedge array in DYNAM3.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'dynam3'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Real, Intent (In) :: Value
		
		! Store value.
		nedge (Index_1, Index_2) = nedge (Index_1,
     > Index_2) + Value
	
	End Subroutine panedge

	Subroutine paneros (Index_1, Index_2, Value)
		! Routine to store (or put) value by adding it to neros array in DYNAM3.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'dynam3'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Real, Intent (In) :: Value
		
		! Store value.
		neros (Index_1, Index_2) = neros (Index_1,
     > Index_2) + Value
	
	End Subroutine paneros

	Subroutine panmsol (Index_1, Index_2, Value)
		! Routine to store (or put) value by adding it to nmsol array in DYNAM3.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'dynam3'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Real, Intent (In) :: Value
		
		! Store value.
		nmsol (Index_1, Index_2) = nmsol (Index_1,
     > Index_2) + Value
	
	End Subroutine panmsol

	Subroutine pantrap (Index_1, Index_2, Value)
		! Routine to store (or put) value by adding it to ntrap array in DYNAM3.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'dynam3'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Real, Intent (In) :: Value
		
		! Store value.
		ntrap (Index_1, Index_2) = ntrap (Index_1,
     > Index_2) + Value
	
	End Subroutine pantrap

	Subroutine papromptdeps (Index_1, Index_2, Value)
		! Routine to store (or put) value by adding it to promptdeps array in PROMPTDEP.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'promptdep'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Real, Intent (In) :: Value
		
		! Store value.
		promptdeps (Index_1, Index_2) = promptdeps (Index_1,
     > Index_2) + Value
	
	End Subroutine papromptdeps

	Subroutine patizs (Index_1, Index_2, Index_3, Value)
		! Routine to store (or put) value by adding it into tizs array in DYNAM3.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'dynam3'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
		Real, Intent (In) :: Value
		
		! Store value.
		tizs (Index_1, Index_2, Index_3) = tizs (Index_1,
     > Index_2, Index_3) + Value
	
	End Subroutine patizs

	Subroutine pawalls (Index_1, Index_2, Index_3, Value)
		! Routine to store (or put) value by adding it into walls array in DYNAM3.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'dynam3'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
		Real, Intent (In) :: Value
		
		! Store value.
		walls (Index_1, Index_2, Index_3) = walls (Index_1,
     > Index_2, Index_3) + Value
	
	End Subroutine pawalls

	Subroutine pawallse (Index, Value)
		! Routine to store (or put) value by adding it to wallse array in DYNAM3.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'dynam3'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value
		
		! Store value.
		wallse (Index) = wallse (Index) + Value
	
	End Subroutine pawallse

	Subroutine pawallse_i (Index, Value)
		! Routine to store (or put) value by adding it to wallse_i array in DYNAM3.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'dynam3'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value
		
		! Store value.
		wallse_i (Index) = wallse_i (Index) + Value
	
	End Subroutine pawallse_i

	Subroutine pawallsi (Index, Value)
		! Routine to store (or put) value by adding it to wallsi array in DYNAM3.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'dynam3'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value
		
		! Store value.
		wallsi (Index) = wallsi (Index) + Value
	
	End Subroutine pawallsi

	Subroutine pawallsn (Index, Value)
		! Routine to store (or put) value by adding it to wallsn array in DYNAM3.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'dynam3'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		Real, Intent (In) :: Value
		
		! Store value.
		wallsn (Index) = wallsn (Index) + Value
	
	End Subroutine pawallsn

	Subroutine pawtdep (Index_1, Index_2, Index_3, Value)
		! Routine to store (or put) value by adding it into wtdep array in COMTOR.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'comtor'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
		Real, Intent (In) :: Value
		
		! Store value.
		wtdep (Index_1, Index_2, Index_3) = wtdep (Index_1,
     > Index_2, Index_3) + Value
	
	End Subroutine pawtdep

	Subroutine pawtsource (Index_1, Index_2, Index_3, Index_4,
     > Value)
		! Routine to store (or put) value by adding it into wtsource array in COMTOR.
		
		! Every good Fortran program has...
		Implicit None
		
		! Included common blocks.
		Include 'params'
		Include 'comtor'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
		Integer, Intent (In) :: Index_4
		Real, Intent (In) :: Value
		
		! Store value.
		wtsource (Index_1, Index_2, Index_3, Index_4) =
     > wtsource (Index_1, Index_2, Index_3, Index_4) + Value
	
	End Subroutine pawtsource

      End Module HC_Put
