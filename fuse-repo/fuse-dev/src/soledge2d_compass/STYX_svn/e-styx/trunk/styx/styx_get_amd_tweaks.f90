subroutine styx_get_amd_tweaks()
  use styx2eirene
  use eirmod_precision
  use eirmod_parmmod
  use eirmod_comprt
  use eirmod_comusr
  use eirmod_comxs
  use eirmod_comprt
  use eirmod_cinit
  use eirmod_clgin 
  use eirmod_ccona

  implicit none

  integer :: i
  character(100) :: AA

! not needed anymore ?
! request checks on first call to styx_xstcx  
  cx_checks=.true.

! in case of a sensitivity study, get multiplicative factors

  if (tweak_chemistry) then
  	
  	write(*,*) '----------------------------------------------------------------'
        write(*,*) '       chemistry tweaks                                         '
        write(*,*) '----------------------------------------------------------------'

  	select case (am_database)

   		case(1)
  		! AMJUEL level 1
  			open(unit=666,file='tweak_AMJUEL_level_1.txt',status='old')

                	allocate(reac_switch(13),fakt(13))
  			allocate(reac_switchE(13),faktE(13))
 	
 	       		do i=1,7
  				read(666,*) AA
  			enddo
  			! atoms
  			read(666,'(A4,e12.4)') reac_switch(1),fakt(1)
  			read(666,*)
  			read(666,'(A4,e12.4)') reac_switchE(1),faktE(1)
  			read(666,*)
  			read(666,*)
  			read(666,'(A4,e12.4)') reac_switch(2),fakt(2)
  			do i=1,6
  				read(666,*)
  			enddo
  			! molecules
  			read(666,'(A4,e12.4)') reac_switch(3),fakt(3)
  			read(666,*)
  			read(666,*)
  			read(666,'(A4,e12.4)') reac_switch(4),fakt(4)
  			read(666,*)
  			read(666,*)
  			read(666,'(A4,e12.4)') reac_switch(5),fakt(5)
  			
  			do i=1,6
  				read(666,*)
  			enddo
  			! test ions
  			read(666,'(A4,e12.4)') reac_switch(6),fakt(6)
  			read(666,*)
  			read(666,*)
  			read(666,'(A4,e12.4)') reac_switch(7),fakt(7)
  			read(666,*)
  			read(666,*)
  			read(666,'(A4,e12.4)') reac_switch(8),fakt(8)

			do i=1,8
  				if (reac_switch(i) == 'off ') then
					write(*,1) i 
1	format('Reaction #',I3,' has been turned off')   				
					fakt(i)=1e-30_dp
  				else
  					write(*,2) i,fakt(i)
2       format('Multiplicative factor for rate coefficient of reaction #',I3,' = ',ES12.4)
  				endif  			
  			enddo

  			if (reac_switchE(1) == 'off') then
  					write(*,*) 'Hydrogenic radiation losses have been turned off'
  		        else
  					write(*,3) faktE(1)
3       format('Multiplicative factor for radiation losses = ',ES12.4)
  			endif 

        		if (.not.direct_coupling) then
   			! in pre-averaged coupling mode, needed to recalculate sources 
  				! atoms FREACA(iatm,ireaca)
  				FREACA(1,1)=fakt(1)
  				FREACA(1,2)=fakt(2)
  				! molecules
   				FREACM(1,1)=fakt(3)
  				FREACM(1,2)=fakt(4)
   				FREACM(1,3)=fakt(5)
  				! test ions
  				FREACI(1,1)=fakt(6)
  				FREACI(1,2)=fakt(7)
  				FREACI(1,3)=fakt(8)

  			endif
  		
    		case(2)
   		! ADAS level 1
  			open(unit=666,file='tweak_ADAS_level_1.txt',status='old')
   			allocate(reac_switch(13),fakt(13))
  			allocate(reac_switchE(13),faktE(13))
 	
 	       		do i=1,7
  				read(666,*) AA
  			enddo
  			! atoms
  			read(666,'(A4,e12.4)') reac_switch(1),fakt(1)
  			read(666,*)
  			read(666,'(A4,e12.4)') reac_switchE(1),faktE(1)
  			read(666,*)
  			read(666,*)
  			read(666,'(A4,e12.4)') reac_switch(2),fakt(2)
  			do i=1,6
  				read(666,*)
  			enddo
  			! molecules
  			read(666,'(A4,e12.4)') reac_switch(3),fakt(3)
  			read(666,*)
  			read(666,*)
  			read(666,'(A4,e12.4)') reac_switch(4),fakt(4)
  			read(666,*)
  			read(666,*)
  			read(666,'(A4,e12.4)') reac_switch(5),fakt(5)
  			
  			do i=1,6
  				read(666,*)
  			enddo
  			! test ions
  			read(666,'(A4,e12.4)') reac_switch(6),fakt(6)
  			read(666,*)
  			read(666,*)
  			read(666,'(A4,e12.4)') reac_switch(7),fakt(7)
  			read(666,*)
  			read(666,*)
  			read(666,'(A4,e12.4)') reac_switch(8),fakt(8)

			do i=1,8
  				if (reac_switch(i) == 'off ') then
					write(*,1) i   				
					fakt(i)=1e-30_dp
  				else
  					write(*,2) i,fakt(i)
  				endif  			
  			enddo

  			if (reac_switchE(1) == 'off') then
  					write(*,*) 'Hydrogenic radiation losses have been turned off'
  		        else
  					write(*,3) faktE(1)
  			endif 

        		if (.not.direct_coupling) then
   			! in pre-averaged coupling mode, needed to recalculate sources 
  				! atoms FREACA(iatm,ireaca)
  				FREACA(1,1)=fakt(1)
  				FREACA(1,2)=fakt(2)
  				! molecules
   				FREACM(1,1)=fakt(3)
  				FREACM(1,2)=fakt(4)
   				FREACM(1,3)=fakt(5)
  				! test ions
  				FREACI(1,1)=fakt(6)
  				FREACI(1,2)=fakt(7)
  				FREACI(1,3)=fakt(8)

  			endif
  		

        	case(4)
        	! AMJUEL level 2
  			open(unit=666,file='tweak_AMJUEL_level_2.txt',status='old')
  			allocate(reac_switch(13),fakt(13))
  			allocate(reac_switchE(13),faktE(13))
 	
 	       		do i=1,7
  				read(666,*) AA
  			enddo
  			! atoms
  			read(666,'(A4,e12.4)') reac_switch(1),fakt(1)
                        read(666,*)
      			read(666,'(A4,e12.4)') reac_switchE(1),faktE(1)    
  			read(666,*)
  			read(666,*)
  			read(666,'(A4,e12.4)') reac_switch(2),fakt(2)
  			do i=1,6
  				read(666,*)
  			enddo
  			! molecules
  			read(666,'(A4,e12.4)') reac_switch(3),fakt(3)
  			read(666,*)
  			read(666,*)
  			read(666,'(A4,e12.4)') reac_switch(4),fakt(4)
  			read(666,*)
  			read(666,*)
  			read(666,'(A4,e12.4)') reac_switch(5),fakt(5)
  			read(666,*)
  			read(666,*)
  			read(666,'(A4,e12.4)') reac_switch(6),fakt(6)
  			do i=1,6
  				read(666,*)
  			enddo
  			! test ions
  			read(666,'(A4,e12.4)') reac_switch(7),fakt(7)
  			read(666,*)
  			read(666,*)
  			read(666,'(A4,e12.4)') reac_switch(8),fakt(8)
  			read(666,*)
  			read(666,*)
  			read(666,'(A4,e12.4)') reac_switch(9),fakt(9)

			do i=1,9
  				if (reac_switch(i) == 'off ') then
					write(*,1) i 				
					fakt(i)=1e-30_dp
  				else
  					write(*,2) i,fakt(i)
  				endif  			
  			enddo

  			if (reac_switchE(1) == 'off') then
  					write(*,*) 'Hydrogenic radiation losses have been turned off'
  		        else
  					write(*,3) faktE(1)
  			endif 

        		if (.not.direct_coupling) then
   			! in pre-averaged coupling mode, needed to recalculate sources 
  				! atoms FREACA(iatm,ireaca)
  				FREACA(1,1)=fakt(1)
  				FREACA(1,2)=fakt(2)
  				! molecules
   				FREACM(1,1)=fakt(3)
  				FREACM(1,2)=fakt(4)
   				FREACM(1,3)=fakt(5)
  				FREACM(1,4)=fakt(6)
  				! test ions
  				FREACI(1,1)=fakt(7)
  				FREACI(1,2)=fakt(8)
  				FREACI(1,3)=fakt(9)

  			endif


  	end select

  	write(*,*) '----------------------------------------------------------------'

  	close(666)
  endif

  ! fast atomic physics (hardwired, level 1)

  addtls=log(MASSP(2)*PMASSA/RMASSP(1))
  fa1=FREACA(1,1)
  fa2=FREACA(1,2)
  fm1=FREACM(1,1)
  fm2=FREACM(1,2)
  fm3=FREACM(1,3)
  fi1=FREACI(1,1)
  fi2=FREACI(1,2)
  fi3=FREACI(1,3)
  fp1=FREACP(1,1)
  if (fa1 == 0._dp) fa1=1._dp
  if (fa2 == 0._dp) fa2=1._dp
  if (fm2 == 0._dp) fm2=1._dp
  if (fm3 == 0._dp) fm3=1._dp
  if (fm1 == 0._dp) fm1=1._dp
  if (fm2 == 0._dp) fm2=1._dp
  if (fm3 == 0._dp) fm3=1._dp
  if (fi1 == 0._dp) fi1=1._dp
  if (fi2 == 0._dp) fi2=1._dp
  if (fi3 == 0._dp) fi3=1._dp
  if (fp1 == 0._dp) fp1=1._dp


end subroutine
