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

      Module HC_Get

        ! jdemod - this routine should not need comhc since each function 
        ! explicitly include the params common block - using comhc duplicates
        ! information. 

	!Use ComHC ! Hydrocarbon-related variables set in DIVIMP input file.

	! Every good Fortran program has...
        Implicit None

        integer,private :: output_unit_hc_alert

      Contains
     
      subroutine hc_get_init_output(outunit)
      implicit none
      integer :: outunit
      !
      ! This routine initializes the output unit for warnings from these 
      ! routines
      !
         output_unit_hc_alert = outunit

      end subroutine hc_get_init_output

     
!--------------------------------
! GET Routines
!--------------------------------

	Real Function galphai (Index_1,Index_2)
		! Routine to read (or get) value from alphai array from REISER.
		
		! Every good Fortran program has...
      use mod_params
      use mod_reiser_com
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'reiser_com'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		galphai = alphai (Index_1,Index_2)
		
	End Function galphai

	Real Function gbratio (Index_1,Index_2)
		! Routine to read (or get) value from bratio array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gbratio = bratio (Index_1,Index_2)
		
	End Function gbratio

	Real Function gbts (Index_1,Index_2)
		! Routine to read (or get) value from bts array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gbts = bts (Index_1,Index_2)
		
	End Function gbts

	Real Function gcetf (Index_1,Index_2)
		! Routine to read (or get) value from cetf array from CYIELD.
		
		! Every good Fortran program has...
                use mod_cyield
      use mod_params
                Implicit None
		
		! Included common blocks.
c	Include 'params'
		!Include 'cyield'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gcetf = cetf (Index_1,Index_2)
		
	End Function gcetf
	
	Real Function gceth (Index_1,Index_2)
		! Routine to read (or get) value from ceth array from CYIELD.
		
		! Every good Fortran program has...
                use mod_cyield
      use mod_params
                Implicit None
		
		! Included common blocks.
c	Include 'params'
		!Include 'cyield'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gceth = ceth (Index_1,Index_2)
		
	End Function gceth

	Logical Function gcflrex ()
		! Routine to read (or get) value from cflrex array from COMMV.
		
		! Every good Fortran program has...
      use mod_params
      use mod_commv
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'commv'
		
		! Declare input variables.
				
		! Stored value.
		gcflrex = cflrex
		
	End Function gcflrex

	Logical Function gcflrin ()
		! Routine to read (or get) value from cflrin array from COMMV.
		
		! Every good Fortran program has...
      use mod_params
      use mod_commv
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'commv'
		
		! Declare input variables.
				
		! Stored value.
		gcflrin = cflrin
		
	End Function gcflrin

	Logical Function gcflrxa ()
		! Routine to read (or get) value from cflrxa array from COMMV.
		
		! Every good Fortran program has...
      use mod_params
      use mod_commv
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'commv'
		
		! Declare input variables.
				
		! Stored value.
		gcflrxa = cflrxa
		
	End Function gcflrxa

	Logical Function gcheckleak ()
		! Routine to read (or get) value from checkleak variable from COMTOR.
		
		! Every good Fortran program has...
      use mod_params
      use mod_comtor
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'comtor'
		
		! Stored value.
		gcheckleak = checkleak
		
	End Function gcheckleak

	Real Function gcleaks (Index_1)
		! Routine to read (or get) value from cleaks array from COMTOR.
		
		! Every good Fortran program has...
      use mod_params
      use mod_comtor
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'comtor'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
				
		! Stored value.
		gcleaks = cleaks (Index_1)
		
	End Function gcleaks

	Real Function gcleaksn ()
		! Routine to read (or get) value from cleaksn variable from COMTOR.
		
		! Every good Fortran program has...
      use mod_params
      use mod_comtor
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'comtor'
		
		! Declare input variables.
				
		! Stored value.
		gcleaksn = cleaksn
		
	End Function gcleaksn

	Real Function gcq (Index_1,Index_2)
		! Routine to read (or get) value from cq array from CYIELD.
		use mod_cyield
		! Every good Fortran program has...
      use mod_params
		Implicit None
		
		! Included common blocks.
c	Include 'params'
		!Include 'cyield'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gcq = cq (Index_1,Index_2)
		
	End Function gcq

	Real Function gctemav ()
		! Routine to read (or get) value from ctemav array from COMTOR.
		
		! Every good Fortran program has...
      use mod_params
      use mod_comtor
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'comtor'
		
		! Declare input variables.
				
		! Stored value.
		gctemav = ctemav
		
	End Function gctemav

	Real Function gctimes (Index_1,Index_2)
		! Routine to read (or get) value from ctimes array from DYNAM4.
		
		! Every good Fortran program has...
      use mod_params
      use mod_dynam4
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'dynam4'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gctimes = ctimes (Index_1,Index_2)
		
	End Function gctimes

	Real Function gcxsc ()
		! Routine to read (or get) value from cxsc variable from COMTOR.
		
		! Every good Fortran program has...
      use mod_params
      use mod_comtor
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'comtor'
		
		! Declare input variables.
				
		! Stored value.
		gcxsc = cxsc
		
	End Function gcxsc

	Real Function gcysc ()
		! Routine to read (or get) value from cysc variable from COMTOR.
		
		! Every good Fortran program has...
      use mod_params
      use mod_comtor
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'comtor'
		
		! Declare input variables.
				
		! Stored value.
		gcysc = cysc
		
	End Function gcysc

	Double Precision Function gddlims (Index_1, Index_2, Index_3)
		! Routine to read (or get) value from ddlims array from DYNAM1.
		
		! Every good Fortran program has...
      use mod_params
      use mod_dynam1
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'dynam1'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
		
		! Read current value.
		gddlims = ddlims (Index_1, Index_2, Index_3)
	
	End Function gddlims
	
	Real Function gdds (Index)
		! Routine to read (or get) value from dds array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		gdds = dds (Index)
		
	End Function gdds

	Double Precision Function gddts (Index_1, Index_2, Index_3)
		! Routine to read (or get) value from ddts array from DYNAM1.
		
		! Every good Fortran program has...
      use mod_params
      use mod_dynam1
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'dynam1'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
				
		! Read current value.
		gddts = ddts (Index_1, Index_2, Index_3)
	
	End Function gddts

	Real Function gdistin (Index_1,Index_2)
		! Routine to read (or get) value from distin array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gdistin = distin (Index_1,Index_2)
		
	End Function gdistin

	Real Function gdistout (Index_1,Index_2)
		! Routine to read (or get) value from distout array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gdistout = distout (Index_1,Index_2)
		
	End Function gdistout

	Real Function geprods (Index)
		! Routine to read (or get) value from eprods array from CNEUT.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cneut
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		geprods = eprods (Index)
		
	End Function geprods

	Integer Function gidds (Index_1,Index_2)
		! Routine to read (or get) value from idds array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gidds = idds (Index_1,Index_2)
		
	End Function gidds

	Integer Function gidprods (Index)
		! Routine to read (or get) value from idprods array from CNEUT.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cneut
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		gidprods = idprods (Index)
		
	End Function gidprods

	Integer Function gikds (Index)
		! Routine to read (or get) value from ikds array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		gikds = ikds (Index)
		
	End Function gikds

	Integer Function giking (Index_1,Index_2)
		! Routine to read (or get) value from iking array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
						
		! Stored value.
		giking = iking (Index_1,Index_2)
		
	End Function giking

	Integer Function gikins (Index_1,Index_2)
		! Routine to read (or get) value from ikins array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
						
		! Stored value.
		gikins = ikins (Index_1,Index_2)
		
	End Function gikins

	Integer Function gikoutg (Index_1,Index_2)
		! Routine to read (or get) value from ikouts array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
						
		! Stored value.
		gikoutg = ikoutg (Index_1,Index_2)
		
	End Function gikoutg

	Integer Function gikouts (Index_1,Index_2)
		! Routine to read (or get) value from ikouts array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
						
		! Stored value.
		gikouts = ikouts (Index_1,Index_2)
		
	End Function gikouts

	Integer Function girds (Index)
		! Routine to read (or get) value from irds array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		girds = irds (Index)
		
	End Function girds

	Integer Function girins (Index_1,Index_2)
		! Routine to read (or get) value from irins array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
						
		! Stored value.
		girins = irins (Index_1,Index_2)
		
	End Function girins

	Integer Function girouts (Index_1,Index_2)
		! Routine to read (or get) value from irouts array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
						
		! Stored value.
		girouts = irouts (Index_1,Index_2)
		
	End Function girouts

	Integer Function gisprods (Index)
		! Routine to read (or get) value from isprods array from CNEUT.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cneut
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		gisprods = isprods (Index)
		
	End Function gisprods

	Real Function gkalphs (Index)
		! Routine to read (or get) value from kalphs array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		gkalphs = kalphs (Index)
		
	End Function gkalphs

	Real Function gkatizs (Index)
		! Routine to read (or get) value from kaatizs array from CNEUT.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cneut
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		gkatizs = katizs (Index)
		
	End Function gkatizs

	Real Function gkbacds (Index_1,Index_2)
		! Routine to read (or get) value from kbacds array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gkbacds = kbacds (Index_1,Index_2)
		
	End Function gkbacds

	Real Function gkbetas (Index)
		! Routine to read (or get) value from kbetas array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		gkbetas = kbetas (Index)
		
	End Function gkbetas
	
	Real Function gkbfs (Index_1,Index_2)
		! Routine to read (or get) value from kbfs array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gkbfs = kbfs (Index_1,Index_2)
		
	End Function gkbfs

	Real Function gkes (Index_1,Index_2)
		! Routine to read (or get) value from kes array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gkes = kes (Index_1,Index_2)
		
	End Function gkes

	Real Function gkfegs (Index_1,Index_2)
		! Routine to read (or get) value from kfegs array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gkfegs = kfegs (Index_1,Index_2)
		
	End Function gkfegs

	Real Function gkfigs (Index_1,Index_2)
		! Routine to read (or get) value from kfigs array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gkfigs = kfigs (Index_1,Index_2)
		
	End Function gkfigs

	Real Function gkfizs (Index_1,Index_2,Index_3)
		! Routine to read (or get) value from kfizs array from CIONIZ.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cioniz
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cioniz'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
				
		! Stored value.
		gkfizs = kfizs (Index_1,Index_2,Index_3)
		
	End Function gkfizs

	Real Function gkfords (Index_1,Index_2)
		! Routine to read (or get) value from kfords array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gkfords = kfords (Index_1,Index_2)
		
	End Function gkfords

	Real Function gkfps (Index_1,Index_2,Index_3)
		! Routine to read (or get) value from kfps array from CIONIZ.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cioniz
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cioniz'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
				
		! Stored value.
		gkfps = kfps (Index_1,Index_2,Index_3)
		
	End Function gkfps

	Real Function gkfss (Index_1,Index_2,Index_3)
		! Routine to read (or get) value from kfss array from CIONIZ.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cioniz
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cioniz'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
				
		! Stored value.
		gkfss = kfss (Index_1,Index_2,Index_3)
		
	End Function gkfss

	Real Function gkfssmod (Index_1,Index_2)
		! Routine to read (or get) value from kfssmod array from CIONIZ.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cioniz
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cioniz'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gkfssmod = kfssmod (Index_1,Index_2)
		
	End Function gkfssmod

	Real Function gkfts (Index_1,Index_2,Index_3)
		! Routine to read (or get) value from kfts array from CIONIZ.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cioniz
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cioniz'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
				
		! Stored value.
		gkfts = kfts (Index_1,Index_2,Index_3)
		
	End Function gkfts
	
	Real Function gkkkfps (Index_1,Index_2,Index_3)
		! Routine to read (or get) value from kkkfps array from CIONIZ.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cioniz
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cioniz'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
				
		! Stored value.
		gkkkfps = kkkfps (Index_1,Index_2,Index_3)
		
	End Function gkkkfps

	Real Function gkks (Index)
		! Routine to read (or get) value from kks array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		
		! Index should be from 1 to MAXNRS.
		If (Index .lt. 1 .or. Index .gt. MAXNRS) Then
			Write (Output_Unit_HC_Alert,*) 
     >                    "Error in function gkks:
     >                    Return value outside range of KKS:",
     >                    Index,MAXNRS
		End If

		! Read current value.
		gkks = kks (Index)
	
	End Function gkks

	Real Function gkmfps (Index)
		! Routine to read (or get) value from kmfps array from CNEUT2.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cneut2
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cneut2'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		gkmfps = kmfps (Index)
		
	End Function gkmfps

	Real Function gkmfss (Index)
		! Routine to read (or get) value from kmfss array from CNEUT2.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cneut2
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cneut2'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		gkmfss = kmfss (Index)
		
	End Function gkmfss

	Real Function gknbs (Index_1,Index_2)
		! Routine to read (or get) value from ktnbs array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gknbs = knbs (Index_1,Index_2)
		
	End Function gknbs

	Integer Function gkorpg (Index_1,Index_2)
		! Routine to read (or get) value from korpg array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gkorpg = korpg (Index_1,Index_2)
		
	End Function gkorpg

	Real Function gkpchs (Index_1,Index_2,Index_3)
		! Routine to read (or get) value from kpchs array from CIONIZ.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cioniz
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cioniz'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
				
		! Stored value.
		gkpchs = kpchs (Index_1,Index_2,Index_3)
		
	End Function gkpchs

	Real Function gkperps (Index_1,Index_2)
		! Routine to read (or get) value from kperps array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gkperps = kperps (Index_1,Index_2)
		
	End Function gkperps

	Real Function gkplos (Index_1,Index_2,Index_3)
		! Routine to read (or get) value from kplos array from CIONIZ.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cioniz
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cioniz'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
				
		! Stored value.
		gkplos = kplos (Index_1,Index_2,Index_3)
		
	End Function gkplos

	Real Function gkoutds (Index_1,Index_2)
		! Routine to read (or get) value from koutds array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gkoutds = koutds (Index_1,Index_2)
		
	End Function gkoutds

	Real Function gkpizs (Index_1,Index_2)
		! Routine to read (or get) value from kpizs array from CIONIZ.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cioniz
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cioniz'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gkpizs = kpizs (Index_1,Index_2)
		
	End Function gkpizs

	Real Function gkrmax (Index)
		! Routine to read (or get) value from krmax array from CNEUT2.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cneut2
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cneut2'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		gkrmax = krmax (Index)
		
	End Function gkrmax

	Real Function gkrmaxw (Index)
		! Routine to read (or get) value from krmaxw array from CNEUT2.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cneut2
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cneut2'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		gkrmaxw = krmaxw (Index)
		
	End Function gkrmaxw

	Real Function gksmaxs (Index)
		! Routine to read (or get) value from ksmaxs array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		
		! Index should be from 1 to MAXNRS.
		If (Index .lt. 1 .or. Index .gt. MAXNRS) Then
			Write (Output_Unit_HC_Alert,*)
     >                    "Error in function gksmaxs:
     >                    Return value outside range of KSMAXS:",
     >                    Index,MAXNRS
		End If

		! Read current value.
		gksmaxs = ksmaxs (Index)
	
	End Function gksmaxs

	Real Function gksb (Index_1,Index_2)
		! Routine to read (or get) value from ksb array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gksb = ksb (Index_1,Index_2)
		
	End Function gksb

	Real Function gkss (Index_1, Index_2)
		! Routine to read (or get) value from kss array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		
		! Index should be from 1 to MAXNKS, Index_2 should be from 1 to MAXNRS.
		If (Index_1 .lt. 1 .or. Index_1 .gt. MAXNKS
     > .or. Index_2 .lt. 1 .or. Index_2 .gt. MAXNRS) Then
			Write (Output_Unit_HC_Alert,*) 
     >                    "Error in function gkss:
     >                    Return value outside range of KSS:",Index_1,
     >                    MAXNRS,MAXNKS,Index_2
		End If

		! Read current value.
		gkss = kss (Index_1, Index_2)
		
	End Function gkss

	Real Function gktebs (Index_1,Index_2)
		! Routine to read (or get) value from ktebs array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gktebs = ktebs (Index_1,Index_2)
		
	End Function gktebs

	Real Function gktibs (Index_1,Index_2)
		! Routine to read (or get) value from ktibs array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gktibs = ktibs (Index_1,Index_2)
		
	End Function gktibs

	Real Function gktids (Index)
		! Routine to read (or get) value from ktids array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		gktids = ktids (Index)
		
	End Function gktids

	Real Function gkvhs (Index_1,Index_2)
		! Routine to read (or get) value from kvhs array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gkvhs = kvhs (Index_1,Index_2)
		
	End Function gkvhs

	Real Function ginjkind (Index)
		! Routine to read (or get) value from injkind array from COMTOR.
		
		! Every good Fortran program has...
      use mod_params
      use mod_comtor
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'comtor'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		ginjkind = injkind (Index)
		
	End Function ginjkind

	Real Function ginjprob (Index)
		! Routine to read (or get) value from injprob array from COMTOR.
		
		! Every good Fortran program has...
      use mod_params
      use mod_comtor
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'comtor'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		ginjprob = injprob (Index)
		
	End Function ginjprob

	Real Function ginjrind (Index)
		! Routine to read (or get) value from injrind array from COMTOR.
		
		! Every good Fortran program has...
      use mod_params
      use mod_comtor
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'comtor'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		ginjrind = injrind (Index)
		
	End Function ginjrind

	Real Function gkareas (Index_1, Index_2)
		! Routine to return (or get) any value in kareas data array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		
		! Return value.
		gkareas = kareas (Index_1, Index_2)
		
	End Function gkareas

	Real Function glambda1 (Index)
		! Routine to read (or get) value from lambda1 array from REISER.
		
		! Every good Fortran program has...
      use mod_params
      use mod_reiser_com
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'reiser_com'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		glambda1 = lambda1 (Index)
		
	End Function glambda1

	Real Function glaunchdat (Index_1, Index_2)
		! Routine to return (or get) any value in launchdat data array from COMTOR2.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cneut
      !use mod_comtor
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'comtor'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		
		! Return value.
		glaunchdat = launchdat (Index_1, Index_2)
		
	End Function glaunchdat

	Real Function glfps (Index_1,Index_2,Index_3)
		! Routine to read (or get) value from lfps array from CLOCAL.
		
		! Every good Fortran program has...
      use mod_params
      use mod_clocal
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'clocal'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
				
		! Stored value.
		glfps = lfps (Index_1,Index_2,Index_3)
		
	End Function glfps

	Real Function glfss (Index_1,Index_2,Index_3)
		! Routine to read (or get) value from lfss array from CLOCAL.
		
		! Every good Fortran program has...
      use mod_params
      use mod_clocal
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'clocal'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
				
		! Stored value.
		glfss = lfss (Index_1,Index_2,Index_3)
		
	End Function glfss

	Real Function glfts (Index_1,Index_2,Index_3)
		! Routine to read (or get) value from lfts array from CLOCAL.
		
		! Every good Fortran program has...
      use mod_params
      use mod_clocal
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'clocal'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
				
		! Stored value.
		glfts = lfts (Index_1,Index_2,Index_3)
		
	End Function glfts
	
	Real Function glllfps (Index_1,Index_2,Index_3)
		! Routine to read (or get) value from lllfps array from CLOCAL.
		
		! Every good Fortran program has...
      use mod_params
      use mod_clocal
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'clocal'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
				
		! Stored value.
		glllfps = lllfps (Index_1,Index_2,Index_3)
		
	End Function glllfps

	Real Function gltolds (Index_1,Index_2,Index_3)
		! Routine to read (or get) value from ltolds array from CLOCAL.
		
		! Every good Fortran program has...
      use mod_params
      use mod_clocal
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'clocal'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
				
		! Stored value.
		gltolds = ltolds (Index_1,Index_2,Index_3)
		
	End Function gltolds
	
	Integer Function gnks (Index)
		! Routine to read (or get) value from nks array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		gnks = nks (Index)
		
	End Function gnks

	Real Function gpinenz (Index_1,Index_2)
		! Routine to read (or get) value from pinenz array from PINDATA.
		
		! Every good Fortran program has...
      use mod_params
      use mod_pindata
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'pindata'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gpinenz = pinenz (Index_1,Index_2)
		
	End Function gpinenz

	Real Function gnh (Index_1,Index_2)
		! Routine to read (or get) value from pinenz array from PINDATA.
		
		! Every good Fortran program has...
      use mod_params
      use mod_pindata
      use mod_comtor
      use mod_cedge2d
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'pindata'
c               include 'comtor'
c               include 'cedge2d'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
                if (cpinopt.eq.1.or.cpinopt.eq.4) then
                   gnh = pinatom(Index_1,Index_2)
                else
                   gnh = e2datom(index_1,index_2)
                endif

		
	End Function gnh

	Real Function gqtim2 ()
		! Routine to read (or get) value from qtim2 array from REISER.
		
		! Every good Fortran program has...
      use mod_params
      use mod_reiser_com
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'reiser_com'
		
		! Stored value.
		gqtim2 = qtim2
		
	End Function gqtim2

	Real Function gqtim ()
		! Routine to read (or get) value from qtim2 array from REISER.
		
		! Every good Fortran program has...
      use mod_params
      use mod_comtor
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'comtor'
		
		! Stored value.
		gqtim = qtim
		
	End Function gqtim

	Real Function gfsrate ()
		! Routine to read (or get) value from qtim2 array from REISER.
		
		! Every good Fortran program has...
      use mod_params
      use mod_comtor
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'comtor'
		
		! Stored value.
		gfsrate = fsrate
		
	End Function gfsrate

	Real Function gplasma_mass ()
		! Routine to read (or get) value from qtim2 array from REISER.
		
		! Every good Fortran program has...
      use mod_params
      use mod_comtor
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'comtor'
		
		! Stored value.
		gplasma_mass = crmb
		
	End Function gplasma_mass

	Real Function granv (Index)
		! Routine to read (or get) value from ranv array from CRAND.
		
		! Every good Fortran program has...
      use mod_params
      use mod_crand
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'crand'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		granv = ranv(Index)
		
	End Function granv

	Real Function granva (Index)
		! Routine to read (or get) value from ranva array from CNEUT.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cneut
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		granva = ranva (Index)
		
	End Function granva

	Real Function granvb (Index)
		! Routine to read (or get) value from ranvb array from CNEUT.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cneut
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		granvb = ranvb (Index)
		
	End Function granvb

	Real Function granvc (Index)
		! Routine to read (or get) value from ranvc array from CNEUT.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cneut
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		granvc = ranvc (Index)
		
	End Function granvc

	Real Function grp (Index)
		! Routine to read (or get) value from rp array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		grp = rp (Index)
		
	End Function grp

	Real Function grs (Index_1,Index_2)
		! Routine to read (or get) value from rs array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		grs = rs (Index_1,Index_2)
		
	End Function grs

	Real Function grvertp (Index_1,Index_2)
		! Routine to return (or get) value from rvertp array in CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Returned value.
		grvertp = rvertp (Index_1,Index_2)
		
	End Function grvertp

	Real Function grw (Index)
		! Routine to return (or get) value to rw array from COMTOR2.
		
		! Every good Fortran program has...
      use mod_params
      use mod_comtor
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'comtor'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Returned value.
		grw = rw (Index)
		
	End Function grw

	Real Function gsatizs (Index)
		! Routine to read (or get) value from satizs array from CNEUT.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cneut
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		gsatizs = satizs (Index)
		
	End Function gsatizs

	Real Function gsnews (Index)
		! Routine to read (or get) value from snews array from CNEUT.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cneut
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		gsnews = snews (Index)
		
	End Function gsnews

	Real Function gsputys (Index)
		! Routine to read (or get) value from sputys array from CNEUT.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cneut
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		gsputys = sputys (Index)
		
	End Function gsputys

	Integer Function gtagdv (Index_1,Index_2)
		! Routine to read (or get) value from tagdv array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
						
		! Stored value.
		gtagdv = tagdv (Index_1,Index_2)
		
	End Function gtagdv

	Real Function gtempds (Index)
		! Routine to read (or get) value from tempds array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		gtempds = tempds (Index)
		
	End Function gtempds

	Real Function gthetag (Index_1,Index_2)
		! Routine to read (or get) value from thetag array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gthetag = thetag (Index_1,Index_2)
		
	End Function gthetag

	Real Function gthetas (Index)
		! Routine to read (or get) value from thetas array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		gthetas = thetas (Index)
		
	End Function gthetas

	Real Function gthetat (Index)
		! Routine to read (or get) value from thetat array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		gthetat = thetat (Index)
		
	End Function gthetat

	Real Function gvins (Index)
		! Routine to read (or get) value from vins array from CNEUT.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cneut
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		gvins = vins (Index)
		
	End Function gvins

	Integer Function gwallindex (Index)
		! Routine to read (or get) value from wallindex array from COMTOR2.
		
		! Every good Fortran program has...
      use mod_params
      use mod_comtor
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'comtor'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		gwallindex = wallindex (Index)
		
	End Function gwallindex

	Real Function gwallpt (Index_1, Index_2)
		! Routine to return (or get) any value in wallpt data array.
		
		! Every good Fortran program has...
      use mod_params
      use mod_comtor
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'comtor'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		
		! Return value.
		gwallpt = wallpt (Index_1, Index_2)
		
	End Function gwallpt
	
	Real Function gwallse (Index)
		! Routine to return (or get) value to wallpt array for purposes of using the stored preset information.
		
		! Every good Fortran program has...
      use mod_params
      use mod_dynam3
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'dynam3' ! wallse delcared as wallse(maxpts+1).
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Returned value.
		gwallse = wallse (Index)
		
	End Function gwallse
	
	Real Function gwtdep (Index_1, Index_2, Index_3)
		! Routine to load (or get) from the wtdep array in COMTOR.
		
		! Every good Fortran program has...
      use mod_params
      use mod_comtor
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'comtor'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
		Integer, Intent (In) :: Index_3
				
		! Returned value.
		gwtdep = wtdep (Index_1, Index_2, Index_3)
		
	End Function gwtdep

	Real Function gxatizs (Index)
		! Routine to read (or get) value from xatizs array from CNEUT.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cneut
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		
		! Read current value.
		gxatizs = xatizs (Index)
	
	End Function gxatizs
	
	Real Function gxprods (Index)
		! Routine to read (or get) value from xprods array from CNEUT.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cneut
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		
		! Index should be from 1 to MAXIMP.
		If (Index .lt. 1 .or. Index .gt. MAXIMP) Then
			Write (Output_Unit_HC_Alert,*) 
     >                    "Error in function gxprods:
     >                    Return value outside range of XPRODS."
		End If

		! Read current value.
		gxprods = xprods (Index)
	
	End Function gxprods

	Real Function gyatizs (Index)
		! Routine to read (or get) value from yatizs array from CNEUT.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cneut
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		
		! Read current value.
		gyatizs = yatizs (Index)
	
	End Function gyatizs
		
	Real Function gyprods (Index)
		! Routine to read (or get) value from yprods array from CNEUT.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cneut
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cneut'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
		
		! Read current value.
		gyprods = yprods (Index)
	
	End Function gyprods

	Real Function gzp (Index)
		! Routine to read (or get) value from zp array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Stored value.
		gzp = zp (Index)
		
	End Function gzp

	Real Function gzs (Index_1,Index_2)
		! Routine to read (or get) value from zs array from CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Stored value.
		gzs = zs (Index_1,Index_2)
		
	End Function gzs

	Real Function gzvertp (Index_1,Index_2)
		! Routine to return (or get) value from zvertp array in CGEOM.
		
		! Every good Fortran program has...
      use mod_params
      use mod_cgeom
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cgeom'
		
		! Declare input variables.
		Integer, Intent (In) :: Index_1
		Integer, Intent (In) :: Index_2
				
		! Returned value.
		gzvertp = zvertp (Index_1,Index_2)
		
	End Function gzvertp

	Real Function gzw (Index)
		! Routine to return (or get) value to zw array from COMTOR2.
		
		! Every good Fortran program has...
      use mod_params
      use mod_comtor
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'comtor'
		
		! Declare input variables.
		Integer, Intent (In) :: Index
				
		! Returned value.
		gzw = zw (Index)
		
	End Function gzw

	Real Function geran (Seed)
		! Put single random value in RAN
		
		! Every good Fortran program has...
      use mod_params
      use mod_cneut
		Implicit None
		
		! Included common blocks.
c	Include 'params'
c	Include 'cneut'

		! Declare input variables.
		Double Precision, Intent (In) :: Seed

		Call SURAND2 (Seed,1,geran)

	End Function geran

        integer function gcion()
      use mod_params
      use mod_comtor
            implicit none
            !     retrieve the value of CION from the DIVIMP common blocks
c           include 'params'
c           include 'comtor'

            gcion = cion

        end function gcion

        real function gabsfac()
      use mod_params
      use mod_comtor
            implicit none
! retrieve the value of ABSFAC from the DIVIMP common blocks
c           include 'params'
c           include 'comtor'

            gabsfac = absfac

        end function gabsfac

        real function grizb()
      use mod_params
      use mod_comtor
            implicit none
! retrieve the value of ABSFAC from the DIVIMP common blocks
c           include 'params'
c           include 'comtor'

            grizb = rizb

        end function grizb


      End Module HC_Get
