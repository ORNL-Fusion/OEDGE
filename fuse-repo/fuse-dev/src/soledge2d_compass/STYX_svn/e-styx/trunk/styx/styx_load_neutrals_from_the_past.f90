subroutine styx_load_neutrals_from_the_past
  use eirmod_precision
  use eirmod_parmmod
  use eirmod_comnnl
  use eirmod_comusr
  use eirmod_comsou
  use eirmod_comprt
  implicit none
  real(dp) :: dtimvo
  integer :: i,j

! piece of code taken from input.f to read fort.15
!  READ INITIAL POPULATION FROM PREVIOUS RUN, OVERWRITE DEFAULTS

  IPRNL=0
  IF (NFILEJ.EQ.2.OR.NFILEJ.EQ.3) THEN
    CALL EIRENE_RSNAP
    DTIMVO=DTIMV

    WRITE (iunout,*) 'INITIAL POPULATION FOR FIRST TIMESTEP'
    WRITE (iunout,*) 'READ FROM FILE FT 15 '
    WRITE (iunout,*) 'PARTICLES AND FLUX STORED FOR INITIAL '
    WRITE (iunout,*) 'DISTRIBUTION IN PREVIOUS RUN '
    CALL EIRENE_MASJ1('IPRNL   ',IPRNL)
    CALL EIRENE_MASR1('FLUX    ',FLUX(NSTRA))

    IF (DTIMVN.NE.DTIMVO) THEN
      FLUX(NSTRA)=FLUX(NSTRA)*DTIMVO/DTIMVN
      WRITE (iunout,*) 'FLUX IS RESCALED BY DTIMV_OLD/DTIMV_NEW '
      CALL EIRENE_MASR1('FLUX    ',FLUX(NSTRA))
      CALL EIRENE_LEER(1)
      ENDIF

      CALL EIRENE_LEER(2)
      IF (TIME0.GT.0.) THEN
        DO I=1,IPRNL
          RPARTC(10,I)=TIME0
        ENDDO
        WRITE (iunout,*) 'PARTICLE CLOCK RESET TO TIME0'
        WRITE (iunout,*) 'FIRST TIMESTEP RUNS FROM TIM1 TO TIM2:  '
        CALL EIRENE_MASR2('TIM1, TIM2      ',TIME0,TIME0+DTIMV)
        CALL EIRENE_LEER(2)
      ENDIF

    ENDIF

    DTIMV=DTIMVN

    IF (NPTS(NSTRA).GT.0.AND.FLUX(NSTRA).GT.0) THEN
      NSRFSI(NSTRA)=1
      SORWGT(1,NSTRA)=1.D0
    ENDIF

end subroutine styx_load_neutrals_from_the_past
