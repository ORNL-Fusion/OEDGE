 
 
      SUBROUTINE EIRENE_iniUSR
      use eirmod_precision
      use eirmod_parmmod
      use eirmod_cgeom
      use eirmod_cgrid
      use eirmod_comprt

      IMPLICIT NONE

      integer :: NQUAD,IQUAD

! condensate 2 triangles into the soledge2D quadrangle for scoring

      NQUAD=(NTRI-1)/2

      if (NQUAD+1 == NRTAL) then

! index of 'virtual' cell are set to 0 to get correct volume
! this is a problem when ploting tallies (see l. 365, eirene_plteir.f)      
      NCLTAL=0

      do IQUAD=1,NQUAD
      	NCLTAL(2*IQUAD-1)=IQUAD
        NCLTAL(2*IQUAD)=IQUAD
      enddo      

! these numbers must be set here to ensure correct value for NSBOX_TAL      
      NR1TAL=NQUAD+1
      NP2TAL=1
      NT3TAL=1

      elseif (NRTAL == NTRI) then
      	write(iunout,*) 'NRTAL = NTRI = NTRI_soledge+1 no cell
     .    condensation in iniusr ...'
      else
      	write(iunout,*) 'Error in iniusr: NQUAD = ', NQUAD, 
     .   'NTRI =  ',NTRI,'NRTAL = ', NRTAL

        ! to be checked
      	!call eirene_exit_own(1)
      endif

      RETURN
      END
