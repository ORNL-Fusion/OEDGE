c     -*-Fortran-*-
c
c ======================================================================
c
      SUBROUTINE InterpolateProfile(mode)  ! Pass quantity to be interpolated... so that this routine is non-Te specific...
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none
 
      INTEGER, INTENT(IN) :: mode

      INTEGER i1,ic,inode1,inode2,in1,in2,target,ion
      REAL*8  s1,s2,v1,v2,deltas,deltav,frac,x,A,B,C,t(0:S28_MAXNKS+1)
      REAL*8, POINTER :: s(:)

      ion = 1

      t = 0.0D0

c...  Find first node in series with Te specified:
      DO inode1 = 1, nnode
        IF ((mode.EQ.1.AND.node(inode1)%te     .NE.0.0).OR.
     .      (mode.EQ.2.AND.node(inode1)%ti(ion).NE.0.0).OR.
     .      (opt%bc(LO).EQ.3.AND.inode1.EQ.1)) EXIT
c        IF (node(inode1)%te.NE.0.0.OR.
c     .      (opt%bc(LO).EQ.3.AND.inode1.EQ.1)) EXIT
      ENDDO
      IF (inode1.GT.nnode) CALL ER('InterpolatePro','NODE1 bad',*99)

c...  Scan through subsequent nodes and interpolate between nodes
c     where Te is defined:
      inode2 = 0
      DO WHILE (inode2.LT.nnode)

c...    Find next node in series with Te data:        
        DO inode2 = inode1+1, nnode
          IF ((mode.EQ.1.AND.node(inode2)%te     .NE.0.0).OR.
     .        (mode.EQ.2.AND.node(inode2)%ti(ion).NE.0.0).OR.
     .        (opt%bc(HI).EQ.3.AND.inode2.EQ.nnode)) EXIT
c          IF (node(inode2)%te.NE.0.0.OR.
c     .        (opt%bc(HI).EQ.3.AND.inode2.EQ.nnode)) EXIT
        ENDDO
        IF (inode1.GE.nnode) CALL ER('InterpolatePro','NODE2 bad',*99)

c...    Check which side of the symmetry point the nodes are on
c       and assign 's' (distance along magnetic field line) accordingly:
c        WRITE(logfp,*) 'PRECRIP:',inode1,mnode,inode2,nnode
c        WRITE(0    ,*) 'PRECRIP:',inode1,mnode,inode2,nnode

        IF     (inode1.LT.mnode.AND.inode2.LE.mnode) THEN
          target = LO
          s => sfor
          in1 = inode1
          in2 = inode2
        ELSEIF (inode1.GE.mnode.AND.inode2.GT.mnode) THEN
          target = HI
          s => sbak
          in1 = inode2
          in2 = inode1
        ELSE
          WRITE(0,*) 'ERROR: NODES STRADDLE SYMMETRY POINT'
          STOP
        ENDIF

c...    Cycle if...:  *** NOT GOOD ENOUGH SINCE YOU COULD HAVE Evolve... CALLED UPSTREAM
c       OF THE TARGET IN A DEFINED REGION... what a mess...
c        IF ((mode.EQ.1.AND.node(in2)%par_mode.EQ.6.AND.
c     .       opt%bc(target).GE.2).OR.
c     .      (mode.EQ.2.AND.node(in2)%par_mode.NE.6)) THEN
c          inode1 = inode2
c          CYCLE
c        ENDIF

c...
        s1 = s(node(in1)%icell)
        s2 = s(node(in2)%icell)
        IF (in1.EQ.1.OR.in1.EQ.nnode) s1 = 0.0

        SELECTCASE (mode) 
          CASE (1)
            v1 = DBLE(node(in1)%te)
            v2 = DBLE(node(in2)%te)
          CASE (2)
            v1 = DBLE(node(in1)%ti(ion))
            v2 = DBLE(node(in2)%ti(ion))
          CASE DEFAULT
            CALL ER('InterpolateProfile','Unknown MODE',*99)
        ENDSELECT

c        WRITE(0,*) 'V1,2:',v1,v2

        deltav = v2 - v1
        deltas = s2 - s1

c        SELECTCASE (node_par_mode(in2)) 
        SELECTCASE (node(in2)%par_mode) 
          CASE (1)
c...        Power law:
            x = DBLE(node(in2)%par_exp)
            DO ic = node(inode1)%icell, node(inode2)%icell
              frac = (s(ic) - s1) / deltas
              t(ic) = deltav * frac**x + v1
            ENDDO
          CASE (2)
c...        Exponential:
            x = DBLE(node(in2)%par_exp)
            DO ic = node(inode1)%icell, node(inode2)%icell
              frac = (s(ic) - s1) / deltas
              t(ic) = deltav * frac**x + v1
            ENDDO
          CASE (3)
c...        Pure conduction (no dependence on power distribution):
            x = 2.0D0 / 7.0D0
            A = v1**(1/x)
            B = v2
            C = (B**(1/x) - A) / deltas
            DO ic = node(inode1)%icell, node(inode2)%icell
              t(ic) = (A + C * (s(ic) - s1))**x 
c              frac = (s(ic) - s1) / deltas   ...old, flawed method...
c              te(ic) = deltav * frac**x + v1
            ENDDO
c          CASE (4)
          CASE (5)
c...        Pleasant:
            x = 2.0D0 / 7.0D0
            DO ic = node(inode1)%icell, node(inode2)%icell
              frac = (s(ic) - s1) / deltas
              t(ic) = (1.0D0 - frac) * (deltav * frac**x      + v1) + 
     .                         frac  * (deltav * frac**0.05D0 + v1)         
            ENDDO
          CASE (6)
c...        Te evolution:
            CALL CalculateTeProfile(in1,in2,s,target)
            t = te
            IF (mode.EQ.2) CALL ER('InterpolateProfile','MODE=2 '// 
     .                             'not ready',*99)
c            CALL EvolveTeProfile(in1,in2,s,target)
          CASE DEFAULT
c            WRITE(0,*) 'DAT:',in2,node_par_mode(in2)
            WRITE(0,*) 'DAT:',in2,node(in2)%par_mode
            STOP 'NO DEFAULT AS YET'
        ENDSELECT 
c...    
        inode1 = inode2
      ENDDO

      SELECTCASE (mode) 
        CASE (1)
          te = t
        CASE (2)
          ti(:,ion) = t(:)
        CASE DEFAULT
          CALL ER('InterpolateProfile','Unknown MODE',*99)
      ENDSELECT

      RETURN
 99   WRITE(0,*) 'MODE      =',mode
      WRITE(0,*) 'NNODE     =',nnode
      WRITE(0,*) 'NODE 1 Te =',node(1:nnode)%te
      WRITE(0,*) 'NODE 1 IC =',node(1:nnode)%icell
      WRITE(0,*) 'NODE 1,2  =',inode1,inode2
      STOP
      END
c
c ======================================================================
c
      SUBROUTINE SpecifyDistribution(target,inpic1,inpic2,mode,exponent,
     .                               val)
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none
 
      INTEGER target,inpic1,inpic2,mode
      REAL*8 :: exponent,val(*)

      INTEGER ic,ic1,ic2,region
      REAL*8  deltas,starts,frac,sumval,A,B,C
      REAL*8, POINTER :: s(:)


      IF (inpic1.LE.0) THEN
        region = inpic1
      ELSE
        region = 1
      ENDIF

c...  Identify region:
      SELECTCASE (target)
        CASE (LO)
          s => sfor
          SELECTCASE (region)        
            CASE( 1)  
              ic1 = inpic1
              ic2 = inpic2
            CASE( 0)  
              ic1 = 1
              ic2 = icmid
            CASE(-1)  
              ic1 = 1
              ic2 = icmid / 2
            CASE(-2)  
              ic1 = 1
              DO ic2 = icmid, 1, -1
                IF (s(ic2).LT.0.10D0*smax) EXIT
              ENDDO
              IF (ic2.LE.ic1) ic2 = ic1 + 1
            CASEDEFAULT
              STOP 'BAD BOY'
          ENDSELECT

        CASE (HI)
          s => sbak
          SELECTCASE (region)        
            CASE( 1)  
              ic1 = inpic1
              ic2 = inpic2
            CASE( 0)  
              ic1 = icmid + 1
              ic2 = icmax
            CASE(-1)  
              ic1 = icmid + 1 + (icmax - icmid) / 2
              ic2 = icmax
            CASE(-2)  
              ic2 = icmax
              DO ic1 = icmid+1, icmax
                IF (s(ic1).LT.0.10D0*smax) EXIT
              ENDDO
              IF (ic1.GE.ic2) ic1 = ic2 - 1
          CASEDEFAULT
              STOP 'REALLY BAD BOY'
          ENDSELECT

        CASE (FULL) 
          s => sfor
          SELECTCASE (region)        
            CASE( 1)  
              ic1 = inpic1
              ic2 = inpic2
            CASE(0)  
              ic1 = 1
              ic2 = icmax
            CASEDEFAULT
              STOP 'REALLY, REALLY, REALLY BAD BOY'
          ENDSELECT

        CASEDEFAULT
          STOP 'NO DEFAULTS'
      ENDSELECT

      deltas = sfor(ic2) - sfor(ic1)
      starts = MIN(s(ic1),s(ic2))

c...  Apply distribution:
      SELECTCASE (mode)

        CASE (0)
c...      None:
          val(ic1:ic2) = 0.0D0

        CASE (1)
c...      Power:
          DO ic = ic1, ic2
            frac = 1.0D0 - (s(ic) - starts) / deltas
            val(ic) = frac**exponent
          ENDDO

        CASE (2)
c...      Exponential:
          B = exponent
          A = 1.0D0 / (1.0D0 - EXP(-1.0D0 / B))
          C = -A * EXP(-1.0D0 / B) 
          DO ic = ic1, ic2
            frac =(s(ic) - starts) / deltas
            val(ic) = A * EXP(-frac / B) + C
c            WRITE(0,*) 'FDSDF:',ic,frac,val(ic)
          ENDDO

        CASE (3)
c...      Power centered at the middle of range:
          IF (ic2-ic1.LE.2) THEN
            val(ic1:ic2) = 1.0D0
          ELSE
            DO ic = ic1, ic2
              frac = DABS(s(ic)-starts - 0.5D0*deltas) / (0.5D0*deltas)
              val(ic) = 1.0D0 - MIN(1.0D0,frac)**exponent
c              WRITE(0,*) 'FRAC:',target,frac,ic1,ic2,val(ic)
            ENDDO
          ENDIF

        CASEDEFAULT
          STOP 'YOU NAUGHTY BOY'
      ENDSELECT         


c...  Cell volume normalization:    ! Do I need an area based normalization as well..?
      sumval = 0.0D0
      DO ic = ic1, ic2
        sumval = sumval + val(ic) * vol(ic)
      ENDDO
      IF (sumval.NE.0.0D0) val(ic1:ic2) = val(ic1:ic2) / sumval


c      IF (log.GT.0) THEN
c        WRITE(logfp,*) 'INTEGR:',sumval,ic1,ic2
c        sumval = 0.0D0
c        DO ic = ic1, ic2
c          sumval = sumval + val(ic) * vol(ic)
c        ENDDO    
c        WRITE(logfp,*) '      :',sumval
c      ENDIF

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c
c *** TO BE DELETED? ***
c
c
      SUBROUTINE ScaleArray(val,ic1,ic2,scale,mode)
      USE mod_sol28_solver
      IMPLICIT none
 
      INTEGER ic1,ic2,mode,exponent
      REAL*8 :: val(*),scale

      INTEGER ic
      REAL*8  sumval


c...  Normalize:
      sumval = 0.0D0
      DO ic = ic1, ic2
        sumval = sumval + val(ic) * vol(ic)
      ENDDO
      IF (sumval.NE.0.0D0) val(ic1:ic2) = val(ic1:ic2) / sumval

c...  Check:
      IF (log.GT.0) THEN  
        WRITE(0,*) 'NORM CHECK:',sumval
        sumval = 0.0D0
        DO ic = ic1, ic2
          sumval = sumval + val(ic) * vol(ic)
        ENDDO    
        WRITE(0,*) '      :',sumval
      ENDIF


      SELECTCASE (mode)

        CASE (0)
c...      Full cell scaling over range:
          val(ic1:ic2) = val(ic1:ic2) * scale

        CASEDEFAULT
          STOP 'YOU NAUGHTY BOY, NO SCALING OPTION FOR THIS'
      ENDSELECT         


      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE IntegrateArray(target,array,mode,integ)
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none
 
      INTEGER target,mode
      REAL*8 :: array(*),integ(0:*)

      INTEGER ic,ic1,ic2
      REAL*8  val
      REAL*8, ALLOCATABLE :: tmp_integ(:)

      IF     (target.EQ.LO.OR.target.EQ.HI) THEN
        ic1 = icbnd1(target)
        ic2 = icbnd2(target)
      ELSEIF (target.EQ.FULL) THEN
        ic1 = 1
        ic2 = icmax
      ELSE
        STOP 'REALLY BAD, BAD BOY'
      ENDIF

      SELECTCASE (mode)
        CASE (0)
c...      Return only the total volume integral:
          val = 0.0D0
          DO ic = ic1, ic2
            val = val + vol(ic) * array(ic)
          ENDDO
          integ(TOTAL) = val

        CASE (1:2)
          IF (mode.EQ.2) THEN
            ALLOCATE(tmp_integ(0:ic2))
            tmp_integ(0:ic2) = integ(0:ic2)
          ENDIF
c...      Standard volume integration over range:
          integ(0:icmax) = 0.0D0
          DO ic = ic1, ic2-1
            val = 0.5D0 * vol(ic) * array(ic)
            integ(ic  ) = integ(ic) + val
            integ(ic+1) = integ(ic) + val
          ENDDO
          val = 0.5D0 * vol(ic2) * array(ic2)
          integ(ic2  ) = integ(ic2) + val
          integ(TOTAL) = integ(ic2) + val
          IF (mode.EQ.2) THEN
            integ(ic1:ic2) = integ(ic1:ic2) + tmp_integ(ic1:ic2)
            integ(TOTAL  ) = integ(TOTAL  ) + tmp_integ(TOTAL  )
            DEALLOCATE(tmp_integ)
          ENDIF

        CASEDEFAULT
          STOP 'YOU NAUGHTY BOY, NO SCALING OPTION FOR THIS'
      ENDSELECT         


      RETURN
 99   STOP
      END

