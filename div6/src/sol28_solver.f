c     -*-Fortran-*-
c
c ======================================================================
c
      SUBROUTINE SolveFluidEquations
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none

      INTEGER ic,ion
      REAL*8  cs

      ion = 1

      cs = -99.0D0

      IF (nion.NE.1)
     .  CALL ER('SolveFluidEquations','One fluid ion only please',*99)

      IF (logop.GT.0) THEN
        WRITE(logfp,*) 'ABOUT TO CALL PLASMA SOLVER'
        WRITE(logfp,*) 
        WRITE(logfp,*) 'ION DATA :',mi(ion),ai(ion)
        WRITE(logfp,*) '         :',zi(ion),ci(ion)
        WRITE(logfp,*) 
        ic = 0
        WRITE(logfp,*) 'TARGET LO:',ic
        WRITE(logfp,*) '  cs,M      :',cs,machno(ic,ion)
        WRITE(logfp,*) '  jsat,area :',isat(ictarg(LO),ion),area(ic)
        WRITE(logfp,*) '  ne,ni     :',ne(ic),ni(ic,ion)
        WRITE(logfp,*) '  vi        :',vi(ic,ion)
        WRITE(logfp,*) '  pe,pi     :',pe(ic),pi(ic,ion)
        WRITE(logfp,*) '  te,ti     :',te(ic),ti(ic,ion)
        ic = icmax + 1
        WRITE(logfp,*) 'TARGET HI:',ic
        WRITE(logfp,*) '  cs,M      :',cs,machno(ic,ion)
        WRITE(logfp,*) '  jsat,area :',isat(ictarg(HI),ion),area(ic)
        WRITE(logfp,*) '  ne,ni     :',ne(ic),ni(ic,ion)
        WRITE(logfp,*) '  vi        :',vi(ic,ion)
        WRITE(logfp,*) '  pe,pi     :',pe(ic),pi(ic,ion)
        WRITE(logfp,*) '  te,ti     :',te(ic),ti(ic,ion)
c        WRITE(logfp,*) 'MOMINT:',momint(1:icmax,ion)
      ENDIF

c...  Assign target conditions to standard fluid arrays:
      

      DO ic = 1, icmax
        par(ic,ion) = parsrc(0,ion) + parint(ic,ion)
        mom(ic,ion) = momsrc(0,ion) + momint(ic,ion)
      ENDDO


      SELECTCASE (sol_option)
        CASE (28)
          CALL SimpleAsPie(ion)
c        CASE (1)
c...      Portal to SOL23:
c         CALL BranchToSOL23
        CASEDEFAULT
          STOP 'NO CUSTOM SOLVERS YET'
      ENDSELECT


c...  For now:
      ne(1:icmax) = ni(1:icmax,ion)


      RETURN
 99   WRITE(0,*) '  NION=',nion
      STOP
      END
c
c ======================================================================
c
      SUBROUTINE SimpleAsPie(ion)
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none

      INTEGER ion

      INTEGER ic
      LOGICAL cont
      REAL*8  par1x,mom1x,a1,b1,c1,r1,cs


      cont = .TRUE.

      IF (logop.GT.0) WRITE(logfp,*) 'ICMAX:',icmax


      anl_imaginary = .FALSE.

      IF (logop.GT.0) THEN
        WRITE(logfp,*) 'SOL28'
        WRITE(logfp,'(3A4,2A10,A6,3A10,5X,2(2A10,2X),4A10)')
     .    'i','iSL','iSH','ni','vi','M','pe+pi','te','ti',
     .    'par','mom','qcond','qconv','a','b','c','r'
      ENDIF

c...  Density and velocity:
      DO WHILE (cont)
        cont = .FALSE.
        IF (logop.GT.0) THEN
          ic = 0
          WRITE(logfp,'(I4,8X,1P,D10.2,D10.2,0P,F6.2,1P,D10.2,
     .                  0P,2F10.4,5X,1P,2D10.2,2X,2D10.2)')
     .      ic,
     .      ni(ic,ion),vi(ic,ion),machno(ic,ion),pe(ic)+pi(ic,ion),
     .      te(ic),ti(ic,ion),
     .      parsrc(ic,ion),momsrc(ic,ion),
     .      qcond(ic),qconv(ic)
c          WRITE(logfp,'(A,I4,8X,42X,1P,2D10.2,2X,5D10.2,0P,3F8.2)')
c     .      '  SOL28:',ic,
c     .      parsrc(ic,ion),momsrc(ic,ion),
c     .      ni(ic,ion),vi(ic,ion),pe(ic)+pi(ic,ion),
c     .      qcond(ic),qconv(ic),te(ic),ti(ic,ion),
c     .      machno(ic,ion)
        ENDIF

        DO ic = 1, icmax

          par1x = par(ic,ion)
          mom1x = mom(ic,ion)

          a1 = (te(ic) + ti(ic,ion)) * ECH
          b1 = -mom1x
          c1 = AMU * 2.0D0 * par1x**2
          r1 = b1**2 - 4.0D0 * a1 * c1
c          b1 = -(p0 + mom1x) 
c          c1 = AMU * crmb * par1x**2

          IF (r1.GE.0.0D0) THEN
c...        Real solution:
            IF ((opt%super(LO).NE.0.AND.ic.LT.anl_ic_super(LO)).OR.
     .          (opt%super(HI).NE.0.AND.ic.GT.anl_ic_super(HI))) THEN
              ni(ic,ion) = (-b1 - DSQRT(r1)) / (2.0D0 * a1)  ! Supersonic
            ELSE
              ni(ic,ion) = (-b1 + DSQRT(r1)) / (2.0D0 * a1)
            ENDIF
          ELSE
c...        Imaginary, approximate density:
            ni(ic,ion) = -b1 / (2.0D0 * a1)
c...        Record extent of super-sonic regions:
            IF     (ic.LE.icmid) THEN
              anl_ic_super (LO) = ic
              anl_imaginary(LO) = .TRUE.
            ELSEIF (ic.GT.icmid.AND..NOT.anl_imaginary(HI)) THEN
              anl_ic_super (HI) = ic
              anl_imaginary(HI) = .TRUE.
            ENDIF
          ENDIF

c...      Set fluid velocity:
          vi(ic,ion) = par1x / ni(ic,ion)

c...      Set fluid pressure:
          ne(ic)     = ni(ic,ion)
          pe(ic)     = ne(ic) * te(ic) * ECH
          pi(ic,ion) = ni(ic,ion) * 
     .                 (ti(ic,ion) * ECH + mi(ion) * vi(ic,ion)**2)

c...      
          cs = DSQRT( (te(ic) + ti(ic,ion)) * ECH / mi(ion) )  ! Need improved calculation...
          machno(ic,ion) = DABS(vi(ic,ion)) / cs

          IF (logop.GT.0) 
     .       WRITE(logfp,'(3I4,1P,2D10.2,0P,F6.2,1P,D10.2,
     .                     0P,2F10.4,5X,1P,2D10.2,2X,2D10.2,2X,4D10.2)')
     .        ic,anl_ic_super(LO),anl_ic_super(HI),
     .        ni(ic,ion),
     .        vi(ic,ion),machno(ic,ion)*DSIGN(1.0D0,vi(ic,ion)),
     .        pe(ic)+pi(ic,ion),
     .        te(ic),ti(ic,ion),
     .        par(ic,ion),mom(ic,ion),
     .        qcond(ic),qconv(ic),
     .        a1,b1,c1,r1
c     .       WRITE(logfp,'(A,3I4,1P,4D10.2,2X,2D10.2,2X,5D10.2,
c     .                     0P,3F8.2)')
c     .        '  SOL28:',ic,anl_ic_super(LO),anl_ic_super(HI),
c     .        a1,b1,c1,r1,
c     .        par(ic,ion),mom(ic,ion),
c     .        ni(ic,ion),vi(ic,ion),pe(ic)+pi(ic,ion),
c     .        qcond(ic),qconv(ic),te(ic),ti(ic,ion),
c     .        machno(ic,ion)

        ENDDO  ! End of IC loop

        IF (logop.GT.0) THEN
          ic = icmax + 1
          WRITE(logfp,'(I4,8X,1P,D10.2,D10.2,0P,F6.2,1P,D10.2,
     .                  0P,2F10.4,5X,1P,2D10.2,2X,2D10.2)')
     .      ic,
     .      ni(ic,ion),vi(ic,ion),machno(ic,ion),pe(ic)+pi(ic,ion),
     .      te(ic),ti(ic,ion),
     .      parsrc(ic,ion),momsrc(ic,ion),
     .      qcond(ic),qconv(ic)
c          WRITE(logfp,'(A,I4,8X,42X,1P,2D10.2,2X,5D10.2,0P,3F8.2)')
c     .      '  SOL28:',ic,
c     .      parsrc(ic,ion),momsrc(ic,ion),
c     .      ni(ic,ion),vi(ic,ion),pe(ic)+pi(ic,ion),
c     .      qcond(ic),qconv(ic),te(ic),ti(ic,ion),
c     .      machno(ic,ion)
        ENDIF

      ENDDO  ! End of main loop



 
      RETURN
 99   STOP
      END
c
c ======================================================================
c
