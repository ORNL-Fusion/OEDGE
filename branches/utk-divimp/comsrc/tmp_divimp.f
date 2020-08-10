c
c
c ======================================================================
c
c subroutine: divGetTdepIndex
c
      INTEGER FUNCTION divGetTdepIndex(random_number)
      USE mod_divimp_tdep

      REAL, INTENT(IN) :: random_number

      INTEGER              :: list_n,i
      INTEGER, ALLOCATABLE :: list_i(:)
      SAVE

      IF (random_number.EQ.-1.0) THEN
        IF (ALLOCATED(list_i)) DEALLOCATE(list_i)
        divGetTdepIndex = -1
      ELSE
        IF (.NOT.ALLOCATED(list_i)) THEN
          ALLOCATE(list_i(tdep_load_n))
          list_n = 0
        ENDIF
        IF (list_n.EQ.0) THEN
          list_n = tdep_load_n
          DO i = 1, list_n
            list_i(i) = i
          ENDDO
        ENDIF
        i = MIN(MAX(1,INT(REAL(list_n)*random_number)),list_n)
        divGetTdepIndex = list_i(i)
c        write(50,*) 'index getter',i,list_i(i),list_n
        list_i(i) = list_i(list_n)
        list_n = list_n - 1
      ENDIF

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: divUpdateIterationCounter
c
      SUBROUTINE divUpdateIterationCounter
      USE mod_divimp
      INCLUDE 'params'

      INTEGER   fp
      LOGICAL   exist
      CHARACTER file*256

      file = 'divimp_counter'
      fp = 99
      INQUIRE(FILE=file,EXIST=exist)
      IF (exist) THEN
        OPEN(fp,FILE=file,ACCESS='SEQUENTIAL',STATUS='OLD'    ,ERR=99)
        READ(fp,*) div_iter
        div_iter = div_iter + 1
        CLOSE(fp)
        OPEN(fp,FILE=file,ACCESS='SEQUENTIAL',STATUS='REPLACE',ERR=99)
      ELSE
        OPEN(fp,FILE=file,ACCESS='SEQUENTIAL',STATUS='NEW'    ,ERR=99)
        div_iter = 1
      ENDIF

      WRITE(fp,*) div_iter

      IF (div_iter.LT.opt_div%niter) THEN       
        CLOSE(fp)
      ELSE
        CLOSE(fp,STATUS='delete')
      ENDIF

      IF (sloutput) THEN 
        WRITE(0,*) ' ------ GOING AGAIN ------',div_iter,opt_div%niter
        WRITE(0,*) 
      ENDIF

      RETURN
99    STOP
      END
c
c
c ======================================================================
c ======================================================================
c ======================================================================
c
c subroutine: LoadRibbonData
c
      SUBROUTINE divLoadRibbonData
      USE mod_sol28_io
      USE mod_divimp
      IMPLICIT none
      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'


      LOGICAL osmGetLine

      INTEGER       fp,i,j,n,ntrace,idum1,i1,i2,iring,itrace,count,cnt,
     .              next_ring,bump1,bump2,store_i,store_i1,store_i2
      LOGICAL       status,cont,debug
      REAL          s1,s2,width,limit,location,diff,mindiff,a,b,c,d,
     .              depth,r1,r2,z1,z2,factor
      CHARACTER     file*256,buffer*1024,tag*1024
      CHARACTER*256 buffer_array(100)

      INTEGER ik,ir,id,k
      REAL    brat,dist(0:MAXNKS),s3,s4,frac,length

      INTEGER, PARAMETER :: MAX_NTRACE = 1000, MAX_NPOINTS = 100

      TYPE :: type_ribbon
         REAL*4    :: version
         INTEGER*4 :: ntrace
         REAL*4    :: rho     (MAX_NTRACE)
         INTEGER*4 :: ntangent(MAX_NTRACE)
      ENDTYPE type_ribbon

      TYPE(type_ribbon) :: ribbon

      TYPE :: type_trace
         INTEGER*4 :: n
         INTEGER*4 :: ntangent
         INTEGER*4 :: code  (MAX_NPOINTS)
         INTEGER*4 :: bump  (MAX_NPOINTS)
         REAL*4    :: rho   
         REAL*4    :: s     (MAX_NPOINTS)
         REAL*4    :: impact(MAX_NPOINTS)
      ENDTYPE type_trace

      TYPE(type_trace), ALLOCATABLE :: trace(:)


      TYPE :: type_ring
         INTEGER*4 :: itrace(2)
         INTEGER*4 :: bump  (2)
         INTEGER*4 :: link    
         INTEGER*4 :: bottom
         INTEGER*4 :: scan
         REAL*4    :: srange(2)
         REAL*4    :: target_r(4)
         REAL*4    :: target_z(4)
         REAL*4    :: pol     (100)
      ENDTYPE type_ring

      INTEGER nring
      TYPE(type_ring) :: ring(MAXNRS)


      debug = .FALSE.


c...  Load the ribbon data:
c     ==================================================================

      fp   = 99
      file = TRIM(opt_div%rib_file) ! 'i-rib-9000a.trace.test_000'
      OPEN(UNIT=fp,FILE=file,ACCESS='SEQUENTIAL',STATUS='OLD',ERR=98)

      SELECTCASE (opt_div%rib_format)
c       ----------------------------------------------------------------       
        CASE (1)
          DO WHILE(osmGetLine(fp,buffer,WITH_TAG))
            DO i = 2, LEN_TRIM(buffer)
              IF (buffer(i:i).EQ.'}') EXIT
            ENDDO
            tag = buffer(2:i-1)
            IF (debug) WRITE(0,*) 'TAG:',TRIM(tag)
            SELECTCASE (TRIM(tag))
c             ----------------------------------------------------------
              CASE ('VERSION')
                 status = osmGetLine(fp,buffer,NO_TAG)
                 READ(buffer,*) ribbon%version
c             ----------------------------------------------------------
              CASE ('TRACE SUMMARY')
                 status = osmGetLine(fp,buffer,NO_TAG)
                 READ(buffer,*) ribbon%ntrace
                 DO i = 1, ribbon%ntrace
                   status = osmGetLine(fp,buffer,NO_TAG)
                   CALL SplitBuffer(buffer,buffer_array) 
                   READ(buffer_array(2),*) ribbon%rho(i)
                   READ(buffer_array(9),*) ribbon%ntangent(i)
                 ENDDO
c             ----------------------------------------------------------
              CASE ('TRACE DATA')
                 ALLOCATE(trace(ribbon%ntrace))
                 DO i = 1, ribbon%ntrace
                   trace(i)%rho      = ribbon%rho     (i)
                   trace(i)%ntangent = ribbon%ntangent(i)
                   status = osmGetLine(fp,buffer,NO_TAG)
                   READ(buffer,*) idum1,trace(i)%n
                   DO j = 1, trace(i)%n
                     status = osmGetLine(fp,buffer,NO_TAG)
                     CALL SplitBuffer(buffer,buffer_array) 
                     READ(buffer_array(2),*) trace(i)%code  (j)
                     READ(buffer_array(3),*) trace(i)%bump  (j)
                     READ(buffer_array(4),*) trace(i)%s     (j)
                     READ(buffer_array(8),*) trace(i)%impact(j)
                   ENDDO
                 ENDDO
c             ----------------------------------------------------------
              CASE DEFAULT
                CALL ER('divLoadRibbonGrid','Unknown tag',*99)
c             ----------------------------------------------------------
            ENDSELECT
          ENDDO
c       ----------------------------------------------------------------       
        CASE DEFAULT
          CALL ER('divLoadRibbonData','Unrecognized format',*99)
      ENDSELECT
      CLOSE(fp)


c...  Build a grid:
c     ==================================================================

c...        
      ring(:)%itrace(1) = -1
      ring(:)%itrace(2) = -1
      ring(:)%link      = -1
      ring(:)%bottom    = -1
      ring(:)%scan      = -1

      nring = 1
      ring(nring)%itrace(1) = 1
      ring(nring)%bump  (1) = trace(1)%bump(1)
      ring(nring)%bump  (2) = trace(1)%bump(trace(1)%n)
      ring(nring)%target_r(4) = trace(1)%rho
      ring(nring)%target_z(4) = trace(1)%s(1)
      ring(nring)%target_r(1) = trace(1)%rho
      ring(nring)%target_z(1) = trace(1)%s(trace(1)%n)

      ntrace = ribbon%ntrace
 
      iring = 1

c...  Build up (potentially partial) rungs 
      cnt = 0

      cont = .TRUE.
      DO WHILE (cont)
        cont = .FALSE.

        cnt = cnt + 1
c        IF (cnt.EQ.10) STOP 'cnt stop'

        ring(iring)%scan = 1

c...    Scan inward to find the next tangency point:
        itrace = ring(iring)%itrace(1)

        s1 = ring(iring)%target_z(4)
        s2 = ring(iring)%target_z(1)

        IF (debug) WRITE(0,*) '====================================='
        IF (debug) WRITE(0,*) 'iring,itrace,s1,s2=',iring,itrace,s1,s2

        DO i = itrace+1, ntrace
          n = trace(i)%n
c...      Find a tangency point in the S-range of interest:
          DO j = 1, n
c            WRITE(0,*) 'scan: ',j,trace(i)%code(j),trace(i)%s(j)
c            IF (trace(i)%code(j).EQ.0.AND.
c     .          trace(i)%s(j).GE.s1.AND.trace(i)%s(j).LE.s2) EXIT
            IF (trace(i)%code(j).EQ.0.AND.
     .          trace(i)%s(j).GE.s1.AND.trace(i)%s(j).LE.s2) THEN
c             Scan to either side of the tangency point to see if the
c             bump value is the same as for the origin trace:
              status = .TRUE.
              DO k = j-1, 1
                IF (trace(i)%code(k).EQ.0.OR.
     .              trace(i)%code(k).EQ.1.OR.
     .              trace(i)%code(k).EQ.2.OR.
     .              trace(i)%code(k).EQ.5) THEN
                  IF (trace(i)%bump(k).NE.ring(iring)%bump(1))
     .              status = .FALSE.
                  EXIT
                ENDIF
              ENDDO               
              DO k = j+1, n
                IF (trace(i)%code(k).EQ.0.OR.
     .              trace(i)%code(k).EQ.1.OR.
     .              trace(i)%code(k).EQ.2.OR.
     .              trace(i)%code(k).EQ.5) THEN
                  IF (trace(i)%bump(k).NE.ring(iring)%bump(2))
     .              status = .FALSE.
                  EXIT
                ENDIF
              ENDDO               
              IF (debug) WRITE(0,*) 'tangency points found, status=',
     .                              status
              IF (status) EXIT
            ENDIF
          ENDDO

          IF (j.LE.n) THEN
c...        Suitable tangency point(s) found:
            DO i1 = j, 1, -1
              IF (trace(i)%code(i1).EQ.1.OR.trace(i)%code(i1).EQ.2.OR.
     .            trace(i)%code(i1).EQ.5) EXIT
            ENDDO
            IF (i1.EQ.0  ) CALL ER('divRibbon...','No link 1',*99)
            DO i2 = j, n
              IF (trace(i)%code(i2).EQ.1.OR.trace(i)%code(i2).EQ.2.OR.
     .            trace(i)%code(i2).EQ.5) EXIT
            ENDDO
            IF (i2.EQ.n+1) CALL ER('divRibbon...','No link 2',*99)

            IF (debug) WRITE(0,*) '  tangency point detected, trace=',i
            IF (debug) WRITE(0,*) '  i1,2=',i1,j,i2
c...        Complete the ring:
            ring(iring)%itrace  (2) = i
            ring(iring)%target_r(3) = trace(i)%rho
            ring(iring)%target_z(3) = trace(i)%s(i1)
            ring(iring)%target_r(2) = trace(i)%rho
            ring(iring)%target_z(2) = trace(i)%s(i2)
            IF (debug) WRITE(0,*) '          == ring completed =='
            EXIT
          ENDIF
        ENDDO  ! i loop

c...    Find a ring that has not been completed yet:
        DO iring = 1, nring
          IF (ring(iring)%scan.EQ.-1) THEN
            cont = .TRUE.
            EXIT
          ENDIF
        ENDDO
c...    No rings found, so try to make a new one:
        IF (iring.EQ.nring+1) THEN
          IF (debug) WRITE(0,*) 'trying to make a new ring.......'
          DO iring = 1, nring
            itrace = ring(iring)%itrace(2)
            IF (ring(iring)%link.EQ.1.OR.itrace.EQ.-1) CYCLE
            s1 = ring(iring)%target_z(3)
            s2 = ring(iring)%target_z(2)
            IF (debug) WRITE(0,*) '  iring,itrace      =',iring,itrace
            IF (debug) WRITE(0,*) '  s1,s2             =',s1,s2
            count = 0
            n = trace(itrace)%n
            DO i = 1, n
              IF (trace(itrace)%code(i).NE.0 .OR.   ! looking for tangency points only
     .            trace(itrace)%s   (i).LT.s1.OR.
     .            trace(itrace)%s   (i).GT.s2) CYCLE
              count = count + 1
              ring(iring)%link = 1
              IF (count.EQ.1) THEN
                i2 = i
                DO i1 = i-1, 1, -1
                  IF (trace(itrace)%code(i1).EQ.0.OR.
     .                trace(itrace)%code(i1).EQ.1.OR.
     .                trace(itrace)%code(i1).EQ.2.OR.
     .                trace(itrace)%code(i1).EQ.5) EXIT
                ENDDO
                IF (i1.EQ.0) CALL ER('divRibbon...','No link 3',*99)
                nring = nring + 1
                ring(nring)%itrace(1) = itrace
                ring(nring)%bump  (1) = trace(itrace)%bump(i1)
                ring(nring)%bump  (2) = trace(itrace)%bump(i2)
                ring(nring)%target_r(4) = trace(itrace)%rho
                ring(nring)%target_z(4) = trace(itrace)%s(i1)
                ring(nring)%target_r(1) = trace(itrace)%rho
                ring(nring)%target_z(1) = trace(itrace)%s(i2)
                next_ring = nring
              ENDIF
              i1 = i
              DO i2 = i+1, n
                IF (trace(itrace)%code(i2).EQ.0.OR.
     .              trace(itrace)%code(i2).EQ.1.OR.
     .              trace(itrace)%code(i2).EQ.2.OR.
     .              trace(itrace)%code(i2).EQ.5) EXIT
              ENDDO
              IF (i2.EQ.n+1) CALL ER('divRibbon...','No link 2',*99)
              nring = nring + 1
              ring(nring)%itrace(1) = itrace
              ring(nring)%bump  (1) = trace(itrace)%bump(i1)
              ring(nring)%bump  (2) = trace(itrace)%bump(i2)
              ring(nring)%target_r(4) = trace(itrace)%rho
              ring(nring)%target_z(4) = trace(itrace)%s(i1)
              ring(nring)%target_r(1) = trace(itrace)%rho
              ring(nring)%target_z(1) = trace(itrace)%s(i2)
            ENDDO  ! i loop
            IF (debug) WRITE(0,*) '  done  nring       =',nring
            EXIT
          ENDDO  ! iring loop
          IF (iring.LE.nring) THEN
            iring = next_ring
            cont = .TRUE.
          ENDIF
        ENDIF
      ENDDO

c      stop 'dsfdfdd'

c...  Now close all rings that weren't closed on tangency points:
      IF (debug) WRITE(0,*)
      DO iring = 1, nring
        IF (ring(iring)%itrace(2).NE.-1) CYCLE
        itrace = ring(iring)%itrace(1)
        IF (debug) WRITE(0,*) 'here we go',iring,itrace
        IF (debug) WRITE(0,*) '     bumps',ring(iring)%bump(1:2)
        width = trace(itrace)%s(trace(itrace)%n) -
     .          trace(itrace)%s(1)
        bump1 = ring(iring)%bump(1)
        bump2 = ring(iring)%bump(2)
        count = 0       
        store_i = -1
        DO i = itrace+1, ntrace
          status = .FALSE.
          i1 = -1
          i2 = -1
          DO j = 1, trace(i)%n
            IF     (i1.EQ.-1.AND.
     .              (trace(i)%code(j).EQ.0.OR.
     .               trace(i)%code(j).EQ.2.OR.
     .               trace(i)%code(j).EQ.5)) THEN
              i1 = j
            ELSEIF (i1.NE.-1.AND.
     .              (trace(i)%code(j).EQ.0.OR.
     .               trace(i)%code(j).EQ.1.OR.
     .               trace(i)%code(j).EQ.5)) THEN
              i2 = j
c              IF (debug) WRITE(0,*) '      ....checking',i,i1,i2,count
c              IF (debug) WRITE(0,*) '                  ',
c     .                trace(i)%bump(i1),
c     .                trace(i)%bump(i2)
              IF (trace(i)%bump(i1).NE.bump1.OR.
     .            trace(i)%bump(i2).NE.bump2) THEN
                i1 = -1
              ELSE
c...            Candidate identified:
                count = count + 1
                depth  = trace(i)%rho - trace(itrace)%rho
                length = trace(i)%s(i2) - trace(i)%s(i1)
c                length = trace(itrace)%s(trace(itrace)%n) - 
c     .                   trace(itrace)%s(1)
                factor = (depth  /  1000.0) /   ! vperp = 1000 m/s, for margin
     .                   (length / 30000.0)     ! vparallel for ~10 eV

c                IF (debug) WRITE(0,*) '          good!',depth,factor

                IF (i.EQ.ntrace.OR.factor.GT.1.0) THEN
c...              Complete the ring:
                  ring(iring)%bottom      = 1  ! At the bottom of a PFR
                  ring(iring)%itrace  (2) = i
                  ring(iring)%target_r(3) = trace(i)%rho
                  ring(iring)%target_z(3) = trace(i)%s(i1)
                  ring(iring)%target_r(2) = trace(i)%rho
                  ring(iring)%target_z(2) = trace(i)%s(i2)
                  IF (debug) WRITE(0,*) '          == ring completed =='
                  IF (debug.AND.i.EQ.ntrace) 
     .              WRITE(0,*) '              end of line     '
                  status = .TRUE.
                ENDIF
                store_i  = i
                store_i1 = i1
                store_i2 = i2
                i1 = -1
                EXIT
              ENDIF
            ENDIF
          ENDDO  ! J loop (along trace)
          IF (j.EQ.trace(i)%n+1) THEN
c...        No valid segment was found on current trace because the 
c           search has gone too far (past the end of the points that
c           define the local PFR), so just use the last valid trace:
            IF (store_i.EQ.-1) THEN
              CALL ER('divLoadRibbonData','No valid outer ring '//
     .                'boundary found',*99)
            ELSE
              ring(iring)%bottom      = 1  ! At the bottom of a PFR
              ring(iring)%itrace  (2) = store_i
              ring(iring)%target_r(3) = trace(store_i)%rho
              ring(iring)%target_z(3) = trace(store_i)%s(store_i1)
              ring(iring)%target_r(2) = trace(store_i)%rho
              ring(iring)%target_z(2) = trace(store_i)%s(store_i2)
              status = .TRUE.
              IF (debug) WRITE(0,*) '          == ring saved =='
            ENDIF
          ENDIF
          IF (status) EXIT
        ENDDO  ! I loop (over traces)
      ENDDO  ! IRING loop

c      stop

c...  Region of interest:
      r1 = opt_div%rib_r1
      z1 = opt_div%rib_z1
      r2 = opt_div%rib_r2
      z2 = opt_div%rib_z2

c      debug = .TRUE.

c...  Decide if rings are too wide, and if yes, split them:
c     ------------------------------------------------------------------
      IF (opt_div%rib_rad_opt.NE.0) THEN
c        r1 =  0.20
c        z1 = -5.00
c        r2 =  0.35
c        z2 =  5.00
c        b = 0.30      ! location of minimum
c        c = 0.10      ! spatial scale
c        d = 0.01      ! minimum width
c        a = 0.05 - d  ! maximum width

c...    Width distribution follows a shifted Gaussian:
c         width_limit = a * exp(-((x-b)^2)/(2.0*c^2)) + d
        b  = opt_div%rib_rad_b      ! location of minimum
        c  = opt_div%rib_rad_c      ! spatial scale
        d  = opt_div%rib_rad_d      ! minimum width
        a  = opt_div%rib_rad_a - d  ! maximum width

        cont = .TRUE.
        DO WHILE (cont)
          cont = .FALSE.
          DO iring = 1, nring
c            if (iring.EQ.1) debug = .TRUE.
c            if (iring.GT.1) debug = .false.
            i1 = ring(iring)%itrace(1)
            i2 = ring(iring)%itrace(2)
            IF (debug) write(0,*) '-->',iring,i1,i2
            IF (i2.EQ.i1+1) CYCLE
            width    = trace(i2)%rho - trace(i1)%rho
            location = 0.5 * (trace(i1)%rho + trace(i2)%rho)
            s1 = ring(iring)%target_z(1)
            s2 = ring(iring)%target_z(4)
            IF (location.GE.r1.AND.location.LE.r2.AND.
     .          ((s1.GE.z1.AND.s1.LE.z2).OR.
     .           (s2.GE.z1.AND.s2.LE.z2))) THEN
              limit = EXP(-((location - b)**2) / (2.0 * c**2))
              limit = a * (1.0 - limit) + d
              IF (debug) WRITE(0,*) '---- rio ---',width,limit
            ELSE
              limit = a + d
            ENDIF
            IF (debug) WRITE(0,*) 'refining....',
     .                             iring,width,limit,location
            IF (width.GT.limit) THEN
c...          Make space:
              DO i = nring, iring+1, -1
                ring(i+1) = ring(i)
              ENDDO
              nring = nring + 1
c...          Find trace that's in the middle of the current ring:
              mindiff = 1.0E+6
              DO j = i1+1, i2-1 
                diff = ABS(trace(j)%rho - location)
                IF (diff.LT.mindiff) THEN
                  mindiff = diff
                  i = j
                ENDIF
              ENDDO
              IF (debug) WRITE(0,*) '    split attempt',i1,i,i2
c...          Find sub-region of the trace:
              bump1 = ring(iring)%bump(1)
              bump2 = ring(iring)%bump(2)
              i1 = -1
              i2 = -1
              DO j = 1, trace(i)%n
                IF     (i1.EQ.-1.AND.
     .                  (trace(i)%code(j).EQ.0.OR.
     .                   trace(i)%code(j).EQ.2.OR.
     .                   trace(i)%code(j).EQ.5)) THEN
                  i1 = j
                ELSEIF (i1.NE.-1.AND.
     .                  (trace(i)%code(j).EQ.0.OR.
     .                   trace(i)%code(j).EQ.1.OR.
     .                   trace(i)%code(j).EQ.5)) THEN
                  i2 = j
                  IF (trace(i)%bump(i1).NE.bump1.OR.
     .                trace(i)%bump(i2).NE.bump2) THEN
                    i1 = -1
                  ELSE
                    IF (debug) WRITE(0,*) '          applied',i,i1,i2
c...                Make new ring:
                    ring(iring+1) = ring(iring)
                    ring(iring+1)%itrace  (1) = i
                    ring(iring+1)%target_r(4) = trace(i)%rho
                    ring(iring+1)%target_z(4) = trace(i)%s(i1)
                    ring(iring+1)%target_r(1) = trace(i)%rho
                    ring(iring+1)%target_z(1) = trace(i)%s(i2)
c...                Adjust old ring:
                    ring(iring  )%bottom      = -2
                    ring(iring  )%itrace  (2) = i
                    ring(iring  )%target_r(3) = trace(i)%rho
                    ring(iring  )%target_z(3) = trace(i)%s(i1)
                    ring(iring  )%target_r(2) = trace(i)%rho
                    ring(iring  )%target_z(2) = trace(i)%s(i2)
                    cont = .TRUE.
                    EXIT
                  ENDIF
                ENDIF
              ENDDO  ! J loop (along trace)
            ENDIF
            IF (cont) EXIT
          ENDDO  ! IRING loop
        ENDDO  ! CONT loop
      ENDIF

c      write(0,*) 's,z=',r1,r2,z1,z2
c      debug = .FALSE.

c
c...  Make the DIVIMP grid:
c     ==================================================================
      IF (.TRUE.) THEN
        brat      = 1.0
        maxrings  = nring

        id = 0

        DO ir = 1, maxrings
          iring = ir

          i1 = ring(iring)%itrace(1)
          i2 = ring(iring)%itrace(2)
          location = 0.5 * (trace(i1)%rho + trace(i2)%rho)
          s1 = ring(iring)%target_z(1)
          s2 = ring(iring)%target_z(4)


c...      Setup the poloidal cell boundaries:
c         --------------------------------------------------------------
          IF (opt_div%rib_pol_opt.NE.0.AND.ring(ir)%bottom.NE.1.AND.
     .        location.GE.r1.AND.location.LE.r2.AND.
     .        ((s1.GE.z1.AND.s1.LE.z2).OR.
     .         (s2.GE.z1.AND.s2.LE.z2))) THEN
c need to process each end of the ring separately!
            n = opt_div%rib_pol_n  ! minimum number of cells on each ring
            a = opt_div%rib_pol_a  ! exponent for log distribution from target to midpoint
            b = opt_div%rib_pol_b  ! required spatial scale at the target, i.e. poloidal size of the cell
            c = opt_div%rib_pol_c  ! maximum cell size - absolute
            d = opt_div%rib_pol_d  ! maximum cell size - relative
          ELSE
            n = opt_div%rib_pol_n_def
            a = opt_div%rib_pol_a_def
            b = opt_div%rib_pol_b_def
            c = opt_div%rib_pol_c_def
            d = opt_div%rib_pol_d_def
          ENDIF
 
c...      Iterate until the near-target spatial requirement is met:
          length = ABS(ring(iring)%target_z(1) - 
     .                 ring(iring)%target_z(4))

          cont = .TRUE.
          DO WHILE (cont)
            cont = .FALSE.
            DO i = 0, n
              dist(i) = REAL(i) / REAL(n)
            ENDDO
            dist(0:n) = dist(0:n) - 0.5
            dist(0:n) = (ABS(dist(0:n))**a) * SIGN(1.0,dist(0:n)) 
            dist(0:n) = (dist(0:n) / dist(n)) * 0.5 + 0.5
            j = NINT( (dist(1) * length) / b)
            IF (j.GT.1) THEN
              IF (n.LT.NINT(0.5*REAL(MAXNKS))) THEN
                n = MIN(j * n, NINT(0.5*REAL(MAXNKS)))
                cont = .TRUE.
              ELSEIF (a.GT.0.01) THEN
                a = MAX(0.01,a*0.3)
                cont = .TRUE.
c              ELSEIF (n.LT.MAXNKS-10) THEN
c                n = MIN(j * n, MAXNKS - 10)
c                cont = .TRUE.
              ENDIF
            ENDIF        
          ENDDO
c...      Check that the maximum cell size constraint is satisified:
c          write(0,*) 'max cell size=',c,d
          cont = .TRUE.
          DO WHILE (cont)
            cont = .FALSE.
            DO i = 1, n
c              IF (((dist(i)-dist(i-1))*length.GT.c)) THEN
              IF (((dist(i)-dist(i-1))*length.GT.c).OR.
     .            ((dist(i)-dist(i-1))       .GT.d)) THEN
c                IF ((dist(i)-dist(i-1)).GT.d) THEN
c                  write(0,*) 'relativeness'
c                ENDIF
                DO j = n+1, i+1, -1
                  dist(j) = dist(j-1)
                ENDDO
                dist(i) = 0.5 * (dist(i-1) + dist(i+1))
                n = n + 1
                cont = .TRUE.
c                IF (debug) WRITE(0,*) 'ref',dist(0:n) 
              ENDIF
              IF (n.GT.MAXNKS-5) THEN
                write(0,*) 'stopping early'
                cont = .FALSE.
                EXIT
              ENDIF
            ENDDO
          ENDDO
c          write(0,*) 'dist',dist(0:n) * length
          nks(ir) = n
c
c...      Add poloidal cuts at tangency points:
c         --------------------------------------------------------------
          IF (ring(iring)%link.EQ.1) THEN
            s1 = ring(iring)%target_z(2)
            s2 = ring(iring)%target_z(3)
            itrace = ring(iring)%itrace(2)
c            IF (debug) WRITE(0,'(A,15F7.3)') '  dist',s2,s1
            DO i = 1, trace(itrace)%n
              IF (trace(itrace)%code(i).NE.0.OR.
     .            trace(itrace)%s   (i).LT.s2.OR.
     .            trace(itrace)%s   (i).GT.s1) CYCLE   
               frac = (trace(itrace)%s(i) - s1) / (s2 - s1)
c              IF (debug) WRITE(0,'(A,15F7.3)') '  frac',frac
c...          Insert the point:
              DO j = 0, nks(ir)-1
                IF (frac.GT.dist(j).AND.frac.LT.dist(j+1)) THEN
                  DO k = nks(ir), j+1, -1
                    dist(k+1) = dist(k)
                  ENDDO
                  dist(j+1) = frac
                  nks(ir) = nks(ir) + 1
                  EXIT
                ENDIF
              ENDDO
            ENDDO
          ENDIF
c          IF (debug) WRITE(0,'(A,15F7.3)') '  dist',dist(1:nks(ir)) 

          i1 = ring(iring)%itrace(1)
          i2 = ring(iring)%itrace(2)

c          IF (debug) WRITE(0,*) 'go...',ir,i1,i2

          psitarg(ir,1) = 0.5 * (trace(i1)%rho + trace(i2)%rho)
          psitarg(ir,2) = psitarg(ir,1)

          idring(ir) = TARTOTAR

          DO ik = 1, nks(ir)
            id = id + 1

            bratio(ik,ir) = brat
            kbfs  (ik,ir) = 1.0 / brat
            bts   (ik,ir) = cbphi 

            korpg (ik,ir) = id
            nvertp(id   ) = 4

            s1 = ring(iring)%target_z(1)
            s2 = ring(iring)%target_z(2)
            s3 = ring(iring)%target_z(3)
            s4 = ring(iring)%target_z(4)
 
            IF (ik.EQ.1) THEN
              rvertp(1,id) = trace(i1)%rho
              rvertp(2,id) = trace(i2)%rho
              zvertp(1,id) = s1
              zvertp(2,id) = s2
            ELSE
              rvertp(1,id) = rvertp(4,id-1)
              rvertp(2,id) = rvertp(3,id-1)
              zvertp(1,id) = zvertp(4,id-1)
              zvertp(2,id) = zvertp(3,id-1)
            ENDIF 

            rvertp(3,id) = rvertp(2,id)
            rvertp(4,id) = rvertp(1,id)
            zvertp(3,id) = (1.0 - dist(ik)) * s2 + dist(ik) * s3
            zvertp(4,id) = (1.0 - dist(ik)) * s1 + dist(ik) * s4

            rs(ik,ir) = 0.0
            zs(ik,ir) = 0.0
            DO i1 = 1, nvertp(id)
              rs(ik,ir) = rs(ik,ir) + rvertp(1,id)
              zs(ik,ir) = zs(ik,ir) + zvertp(1,id)
            ENDDO
            rs(ik,ir) = rs(ik,ir) / REAL(nvertp(id))
            zs(ik,ir) = zs(ik,ir) / REAL(nvertp(id))

          ENDDO  ! IK
        ENDDO  ! IR

        npolyp  = id
        vpolmin = (MAXNKS*MAXNRS - npolyp) / 2 + npolyp
        vpolyp  = vpolmin
       
        ikto = 2
        ikti = 3
       
        rves  = 0
        rvesm = 0
       
        irsep  = 1
        irwall = maxrings 
        irtrap = irwall
        nrs    = irwall
        nbr    = 0
       
c...    Necessary..? 
        cutring = 1
        cutpt1 = ikto
        cutpt2 = ikti
       

        CALL OutputData(87,'Ribbon factory')


      ENDIF







      IF (ALLOCATED(trace)) DEALLOCATE(trace)


      WRITE(0,*) 'VERSION!',ribbon%version
      WRITE(0,*) ' FILE= ',TRIM(file)

      RETURN
 98   CALL ER('divLoadRibbonGrid','Data file not found',*99)
 99   WRITE(0,*) ' FILE= ',TRIM(file)
      WRITE(0,*) ' TAG = ',TRIM(tag)
      STOP
      END
c
c ======================================================================
c
c subroutine: CalcTubeDimentions
c
      SUBROUTINE CalcTubeDimensions(tube_3D_data,dangle)
      IMPLICIT none

c      SUBROUTINE CalcTubeDimensions(xin,yin,zin,
c     .             MAXSURFACE,MAXPOINTS,nsur,npts,hsur,vsur)
c      IMPLICIT none
c
c      INTEGER, INTENT(IN) :: MAXSURFACE,MAXPOINTS
c      REAL   , INTENT(IN) :: xin,yin,zin
c      INTEGER, INTENT(OUT) :: nsur,npts(MAXSURFACE),hsur(MAXSURFACE)
c      REAL*8 , INTENT(OUT) :: vsur(3,MAXPOINTS,0:MAXSURFACE)

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      REAL tube_3D_data(5,MAXNKS,MAXNRS),dangle

      INTEGER ir,isegment,nsegment,mode,ncell,icell
      REAL    rho_start,rho_end,rho1,rho2,length_total

      TYPE    type_tube_3D
         REAL      :: dangle
         INTEGER   :: ncell
         REAL*8    :: area  (1000)
         REAL*8    :: length(1000)
         REAL*8    :: volume(1000)
         INTEGER   :: n
         INTEGER   :: index(  200000)
         REAL*8    :: v    (3,200000)
      ENDTYPE type_tube_3D

      TYPE(type_tube_3D) :: tube_3D(0:5)

      WRITE(0,*) '---------------------------------'
      WRITE(0,*) 'NOT CALCULATING TUBE DIMENSIONS!'
      WRITE(0,*) '---------------------------------'
      RETURN


      nsegment = 25

      DO ir = 5, 9 !  irsep, irsep+1

        rho_start = rho(ir,IN14 ) * 1.0001
        rho_end   = rho(ir,OUT23) * 0.9999

        mode = 1

        DO isegment = 1, nsegment

          rho1 = rho_start + REAL(isegment-1) / REAL(nsegment) * 
     .                       (rho_end - rho_start)
          rho2 = rho_start + REAL(isegment  ) / REAL(nsegment) * 
     .                       (rho_end - rho_start)

          CALL ProcessFluxTube(rho1,rho2,REAL(ir),mode,tube_3D(1:5))
c          CALL ProcessFluxTube(rho1,rho2,REAL(ir),
c     .           mode,tube_3D(1:5),
c     .           MAXSURFACE,MAXPOINTS,nsur,npts,hsur,vsur)

          IF (mode.EQ.1) THEN
            tube_3D(0) = tube_3D(1)
            ncell = tube_3D(0)%ncell
          ELSE
            tube_3D(0)%area   = tube_3D(0)%area   + tube_3D(1)%area   
            tube_3D(0)%length = tube_3D(0)%length + tube_3D(1)%length 
            tube_3D(0)%volume = tube_3D(0)%volume + tube_3D(1)%volume
          ENDIF

          mode = 2
        ENDDO

        tube_3D(0)%length = tube_3D(0)%length / REAL(nsegment) 

        dangle = tube_3D(0)%dangle

        length_total = 0.0
        DO icell = 1, ncell

          tube_3D_data(1,icell,ir) = length_total

          length_total = length_total + SNGL(tube_3D(0)%length(icell))

          WRITE(6,'(A,2I6,2F10.5,2X,2F10.5)') 
     .      'TOTAL VOLUME:',icell,ir,
     .      SNGL(tube_3D(0)%volume(icell) *
     .      (360.0D0/tube_3D(0)%dangle)),
     .      kvols(icell,ir),
     .      length_total,
     .      ksb(icell,ir)

          tube_3D_data(2,icell,ir) = length_total
          tube_3D_data(3,icell,ir) = length_total - 
     .                    0.5 * SNGL(tube_3D(0)%length(icell))       
          tube_3D_data(4,icell,ir) = tube_3D(0)%area  (icell)
          tube_3D_data(5,icell,ir) = tube_3D(0)%volume(icell)

        ENDDO
      ENDDO

      RETURN
 99   STOP
      END
c     
c ======================================================================
c
c subroutine: SelectGridRegoin_DIVIMP
c
      SUBROUTINE ProcessFluxTube(xin,yin,zin,mode,tube_3D)
c      SUBROUTINE ProcessFluxTube(xin,yin,zin,
c     .             mode,tube_3D,
c     .             MAXSURFACE,MAXPOINTS,nsur,npts,hsur,vsur)

      IMPLICIT none

c      INTEGER, INTENT(IN) :: MAXSURFACE,MAXPOINTS
      INTEGER, INTENT(IN) :: mode
      REAL   , INTENT(IN) :: xin,yin,zin
c      INTEGER, INTENT(OUT) :: nsur,npts(MAXSURFACE),hsur(MAXSURFACE)
c      REAL*8 , INTENT(OUT) :: vsur(3,MAXPOINTS,0:MAXSURFACE)

      TYPE    type_tube_3D
         REAL      :: dangle
         INTEGER   :: ncell
         REAL*8    :: area  (1000)
         REAL*8    :: length(1000)
         REAL*8    :: volume(1000)
         INTEGER   :: n
         INTEGER   :: index(  200000)
         REAL*8    :: v    (3,200000)
      ENDTYPE type_tube_3D

      TYPE(type_tube_3D) :: tube_3D(5)

      REAL    FindSeparatrixRadius

      INTEGER i1,i2,itube,ik,ike,nvtx,isegment,nsegments,ring,val_ring
      LOGICAL last_point,polygon_good
      REAL    x1,y1,z1,dradius,rseparatrix,rstart,rend,rvalue(3)   ! *** fix the 1000 and the 10000
      REAL*8  len1,len2,x,z,angle,dangle,vtx(1:3,1000),length(1000),
     .        center(3),vector(3),p1(3),p2(3),p3(3),A,B,C,denominator,u,
     .        dist,mindist,vertex(1:3,5),frac,dot_product,side_length,
     .        cross_product(3),diagonal_42(3),diagonal_53(3),area,
     .        volume,volume_total

      INTEGER n,index(10000)
      REAL    fraction(10000)
      REAL*8  v(3,10000)



      INTEGER ntube_3D
      TYPE(type_tube_3D) :: save_tube_3D(2)


      len1 = 1000.0D0
      len2 = 1000.0D0

      IF (mode.EQ.2) THEN
        save_tube_3D(1) = tube_3D(3)
        save_tube_3D(2) = tube_3D(4)
      ENDIF


      ntube_3D = 0

      dradius = 25.994 * 0.001  ! 0.00054 ! 0.0002
      dangle = 1.0D0

c
c     Looking down the field line:
c
c       4     5
c
c          1
c
c       3     2
c

      rstart = xin !  50.96412 * 0.001
      rend   = yin ! 102.95100 * 0.001
      val_ring = NINT(zin)

      rseparatrix = FindSeparatrixRadius(1)

c      x1 = x1 - dradius
 
      nsegments = 1

      volume_total = 0.0D0

      DO isegment = 1, nsegments

        rvalue(1) = rstart + REAL(isegment-1) / REAL(nsegments) * 
     .                       (rend - rstart)
        rvalue(3) = rstart + REAL(isegment  ) / REAL(nsegments) * 
     .                       (rend - rstart)
        rvalue(2) = 0.5 * (rvalue(1) + rvalue(3))


        WRITE(0,*) 'RVALUES:',rvalue(1:3)
        WRITE(0,*) 'RVALUES:',xin-dradius,xin,xin+dradius
        WRITE(0,*) 'RVALUES:',rvalue(1:3)+rseparatrix
        WRITE(0,*) 'RVALUES:',
     .     rseparatrix+xin-dradius,xin+rseparatrix,
     .     rseparatrix+xin+dradius

        x1 = rseparatrix + rvalue(2)
        y1 = 0.0
        z1 = 0.0

        CALL TraceFieldLine_DIVIMP(x1,y1,z1,2,1,len1,len2,0.0D0,n,v,
     .                             index,fraction,ring,10000)
        IF (ring.NE.val_ring) 
     .    CALL ER('ProcessFluxTube','Invalid ring A',*99)
        ntube_3D = ntube_3D + 1
        tube_3D(ntube_3D)%n = n
        tube_3D(ntube_3D)%index(1:n) = index(1:n)
        DO i1 = 1, n
          tube_3D(ntube_3D)%v(1:3,i1) = v(1:3,i1)
c          IF (index(i1).GT.21) 
c     .      WRITE(0,'(A,2I6,F10.2,2X,3F12.5)')
c     .        'index:',i1,index(i1),fraction(i1),
c     .        SNGL(tube_3D(ntube_3D)%v(1:3,i1))
        ENDDO

c        x1 = x1 + dradius
c        x1 = x1 - dradius
        IF (mode.EQ.1) THEN
          x1 = rseparatrix + rvalue(1)
          CALL TraceFieldLine_DIVIMP(x1,y1,z1,2,1,len1,len2,0.0D0,n,v,
     .                               index,fraction,ring,10000)
          IF (ring.NE.val_ring) 
     .      CALL ER('ProcessFluxTube','Invalid ring B',*99)
          ntube_3D = ntube_3D + 1
          tube_3D(ntube_3D)%n = n
          tube_3D(ntube_3D)%index(1:n) = index(1:n)
          DO i1 = 1, n
          tube_3D(ntube_3D)%v(1:3,i1) = v(1:3,i1)
c          IF (index(i1).GT.21) 
c     .      WRITE(0,'(A,2I6,F10.2,2X,3F12.5)')
c     .        'index:',i1,index(i1),fraction(i1),
c     .        SNGL(tube_3D(ntube_3D)%v(1:3,i1))
          ENDDO
c...      Rotate these points toroidally:
          angle = -0.5D0 * dangle * 3.1415927D0 / 180.0D0
          DO i1 = 1, n
            x = tube_3D(ntube_3D)%v(1,i1)
            z = tube_3D(ntube_3D)%v(3,i1)
            tube_3D(ntube_3D)%v(1,i1) =DCOS(angle) * x - DSIN(angle) * z
            tube_3D(ntube_3D)%v(3,i1) =DSIN(angle) * x + DCOS(angle) * z
          ENDDO
        ELSE
          ntube_3D = ntube_3D + 1
          tube_3D(ntube_3D) = save_tube_3D(1)
        ENDIF

c        x1 = x1 + 2.0 * dradius
        x1 = rseparatrix + rvalue(3)
        CALL TraceFieldLine_DIVIMP(x1,y1,z1,2,1,len1,len2,0.0D0,n,v,
     .                             index,fraction,ring,10000)
        IF (ring.NE.val_ring) 
     .    CALL ER('ProcessFluxTube','Invalid ring C',*99)
        ntube_3D = ntube_3D + 1
        tube_3D(ntube_3D)%n = n
        tube_3D(ntube_3D)%index(1:n) = index(1:n)
        DO i1 = 1, n
          tube_3D(ntube_3D)%v(1:3,i1) = v(1:3,i1)
        ENDDO
c...    Rotate these points toroidally:
        angle = -0.5D0 * dangle * 3.1415927D0 / 180.0D0
        DO i1 = 1, n
          x = tube_3D(ntube_3D)%v(1,i1)
          z = tube_3D(ntube_3D)%v(3,i1)
          tube_3D(ntube_3D)%v(1,i1) = DCOS(angle) * x - DSIN(angle) * z
          tube_3D(ntube_3D)%v(3,i1) = DSIN(angle) * x + DCOS(angle) * z
        ENDDO

        ntube_3D = ntube_3D + 1
        tube_3D(ntube_3D) = tube_3D(ntube_3D-1)      
        angle = dangle * 3.1415927D0 / 180.0D0
        DO i1 = 1, n
          x = tube_3D(ntube_3D)%v(1,i1)
          z = tube_3D(ntube_3D)%v(3,i1)
          tube_3D(ntube_3D)%v(1,i1) = DCOS(angle) * x - DSIN(angle) * z
          tube_3D(ntube_3D)%v(3,i1) = DSIN(angle) * x + DCOS(angle) * z
        ENDDO

        IF (mode.EQ.1) THEN
          ntube_3D = ntube_3D + 1
          tube_3D(ntube_3D) = tube_3D(2)      
          angle = dangle * 3.1415927D0 / 180.0D0
          DO i1 = 1, n
            x = tube_3D(ntube_3D)%v(1,i1)
            z = tube_3D(ntube_3D)%v(3,i1)
            tube_3D(ntube_3D)%v(1,i1) =DCOS(angle) * x - DSIN(angle) * z
            tube_3D(ntube_3D)%v(3,i1) =DSIN(angle) * x + DCOS(angle) * z
          ENDDO
        ELSE
          ntube_3D = ntube_3D + 1
          tube_3D(ntube_3D) = save_tube_3D(2)
        ENDIF


c        IF (mode.EQ.1) THEN
c          DO i2 = 1, ntube_3D
c            IF (mode.EQ.2.and.(i1.EQ.2.OR.i1.EQ.5)) CYCLE
c            DO i1 = 1, tube_3D(i2)%n-1
c              nsur = nsur + 1
c              hsur(nsur) = -1
c              npts(nsur) =  2
c              vsur(1:3,1,nsur) = tube_3D(i2)%v(1:3,i1  )
c              vsur(1:3,2,nsur) = tube_3D(i2)%v(1:3,i1+1)
c            ENDDO
c          ENDDO
c        ENDIF

c       CYCLE

c      RETURN




       itube = 1
       ike = 0
       DO i1 = 1, tube_3D(itube)%n 
         ike = MAX(ike,tube_3D(itube)%index(i1))
       ENDDO
       WRITE(0,*) 'IKE:',ike

       tube_3D(1)%ncell = ike
       tube_3D(1)%dangle = dangle
       tube_3D(1)%area  (1:ike) = 0.0
       tube_3D(1)%length(1:ike) = 0.0
       tube_3D(1)%volume(1:ike) = 0.0

c...  START OF LOOP:
c        DO ik = 20, 36
        DO ik = 1, ike

c...      Collect line segments associated with this cell:
          last_point = .FALSE.
          nvtx = 0

        
c...      Select the line segments that are associated with the
c         cell IK:
          itube = 1
          DO i1 = 1, tube_3D(itube)%n 
            IF (tube_3D(itube)%index(i1).EQ.ik) THEN
              last_point = .TRUE.
              nvtx = nvtx + 1
              vtx(1:3,nvtx) =  tube_3D(itube)%v(1:3,i1)
            ELSEIF (last_point.AND.ik.LT.ike) THEN
              nvtx = nvtx + 1
              vtx(1:3,nvtx) = tube_3D(itube)%v(1:3,i1)          
              EXIT
            ENDIF
          ENDDO
        
          DO i1 = 1, nvtx  
c            WRITE(0,'(A,I6,3F12.5)') 'VTX:',i1,SNGL(vtx(1:3,i1))
          ENDDO

c...      Calulate the total length:
          length = 0.0D0
          DO i1 = 1, nvtx-1
            length(i1+1) = length(i1) +
     .               DSQRT((vtx(1,i1+1)-vtx(1,i1))**2 +
     .                     (vtx(2,i1+1)-vtx(2,i1))**2 +
     .                     (vtx(3,i1+1)-vtx(3,i1))**2)
          ENDDO
c...      Find half way point:
          DO i2 = 1, nvtx-1
            IF (length(i2  ).LE.0.5D0*length(nvtx).AND.
     .          length(i2+1).GT.0.5D0*length(nvtx)) THEN
              frac = (0.5D0*length(nvtx) - length(i2)) / 
     .               (      length(i2+1) - length(i2))
              center(1:3) = (1.0D0 - frac) * vtx(1:3,i2  ) +
     .                               frac  * vtx(1:3,i2+1)
              vector(1:3) = vtx(1:3,i2+1) - vtx(1:3,i2)
              EXIT
            ENDIF
          ENDDO
c          DO i1 = 1, nvtx
c            WRITE(0,*) 'LENGTH:',i1,length(i1),i2,frac
c          ENDDO
          WRITE(0,*) 'CENTER:',center(1:3)
          WRITE(0,*) 'VECTOR:',vector(1:3)

          vector(1:3) = vector(1:3) /
     .                 DSQRT(vector(1)**2 + vector(2)**2 + vector(3)**2)
        
c...      Find intersections between the flux-tube boundary field line tracings
c         and the plane defined by the center point of the cell:
          itube = 2
          A = vector(1)
          B = vector(2)
          C = vector(3)
          vertex = 0.0D0
          vertex(1:3,1) = center(1:3)
          DO itube = 2, 5
            mindist = 1.0D+20
            DO i1 = 1, tube_3D(itube)%n-1
              p1(1:3) = tube_3D(itube)%v(1:3,i1  ) - center(1:3)
              p2(1:3) = tube_3D(itube)%v(1:3,i1+1) - center(1:3)
              denominator = A * (p1(1) - p2(1)) + B * (p1(2) - p2(2)) + 
     .                      C * (p1(3) - p2(3))
              IF (DABS(denominator).LT.1.0D-10) THEN
                u = -1.0D0
              ELSE
                u = (A * p1(1) + B * p1(2) + C * p1(3)) / denominator
              ENDIF
              IF (u-1.0D-10.GT.0.0D0.AND.u+1.0D-10.LT.1.0D0) THEN
                p3(1:3) = p1(1:3) + u * (p2(1:3) - p1(1:3))
                dist = DSQRT(p3(1)**2 + p3(2)**2 + p3(3)**2)
                IF (dist.LT.mindist) THEN
                  vertex(1:3,itube) = p3(1:3) + center(1:3)
                  mindist = dist
c                  WRITE(0,*) 'U:',i1,u,dist
c                  WRITE(0,*) 'P3:',p3(1:3) + center(1:3)
c                  WRITE(0,*) 'BINGO!',itube
                ENDIF
              ENDIF
            ENDDO
          ENDDO

c...      Check that the plane made up of the identified points is indeed
c         perpendicular to the central vector of the cell:
          DO i1 = 1, 4
            i2 = i1 + 1
            IF (i2.EQ.5) i2 = 1        
            dot_product = 
     .        vector(1) * (vertex(1,i1)-vertex(1,i2)) + 
     .        vector(2) * (vertex(2,i1)-vertex(2,i2)) + 
     .        vector(3) * (vertex(3,i1)-vertex(3,i2))
            IF (dot_product.GT.1.0D-10) THEN
              WRITE(0,*) 'CROSS_PRODUCT PAIN'
              STOP
            ENDIF
          ENDDO

c...      Check the length of each side:
          polygon_good = .TRUE.
          DO i1 = 2, 5
            i2 = i1 + 1
            IF (i2.EQ.6) i2 = 1               
            side_length = DSQRT((vertex(1,i1)-vertex(1,i2))**2 + 
     .                          (vertex(2,i1)-vertex(2,i2))**2 + 
     .                          (vertex(3,i1)-vertex(3,i2))**2)
            WRITE(0,*) 'SIDE LENGTH:',i1-1,side_length
            IF (side_length.GT.0.5) polygon_good = .FALSE.
          ENDDO

          WRITE(0,*) 'POLYGON_GOOD:',polygon_good

          IF (.FALSE.) THEN
            length(nvtx) = 1.0

            vertex(1,2) = 1.0
            vertex(2,2) = 0.0
            vertex(3,2) = 0.0

            vertex(1,3) = 1.0
            vertex(2,3) = 1.0
            vertex(3,3) = 0.0

            vertex(1,4) = 2.0
            vertex(2,4) = 1.0
            vertex(3,4) = 0.0

            vertex(1,5) = 2.0
            vertex(2,5) = 0.0
            vertex(3,5) = 0.0
 
            vector(1) = 0.0
            vector(2) = 0.0
            vector(3) = 1.0
          ENDIF

          IF (polygon_good) THEN
! AREA of a 3D quadrelateral, n is unit normal vector of the plane...
      ! 2 A = n DOT ( v2 - v0 cross v3 - v1 ) 
          
c...       Calculate the area of a quadrelateral:
           diagonal_42(1:3) = vertex(1:3,4) - vertex(1:3,2)
           diagonal_53(1:3) = vertex(1:3,5) - vertex(1:3,3)           
        
           cross_product(1) = diagonal_42(2) * diagonal_53(3) -   ! I had a function for this somewhere...
     .                        diagonal_42(3) * diagonal_53(2)
           cross_product(2) = diagonal_42(3) * diagonal_53(1) -
     .                        diagonal_42(1) * diagonal_53(3)
           cross_product(3) = diagonal_42(1) * diagonal_53(2) -
     .                        diagonal_42(2) * diagonal_53(1)
        
           area = DABS(0.5D0 * (vector(1) * cross_product(1) + 
     .                          vector(2) * cross_product(2) + 
     .                          vector(3) * cross_product(3)))

          ELSE

            area = 0.0D0

          ENDIF

          volume = area * length(nvtx)

          volume_total = volume_total + volume

          WRITE(0,*) 'AREA,VOLUME:',ik,area,
     .               volume*(360.0D0/dangle)  ! *** LEFT OFF ***

          tube_3D(1)%area  (ik) = area
          tube_3D(1)%length(ik) = length(nvtx)
          tube_3D(1)%volume(ik) = volume

c          DO i1 = 2, 5
c            i2 = i1 + 1
c            IF (i2.EQ.6) i2 = 2
c            nsur = nsur + 1
c            hsur(nsur) = -2
c            npts(nsur) =  2
c            vsur(1:3,1,nsur) = vertex(1:3,i1)
c            vsur(1:3,2,nsur) = vertex(1:3,i2)
c          ENDDO
    
        ENDDO

      ENDDO

c      WRITE(0,*) '*** TOTAL VOLUME ***',volume_total*(360.0D0/dangle)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: SelectGridRegoin_DIVIMP
c
      SUBROUTINE SelectGridRegion_DIVIMP(rhoval,nrings,rings,MAX_IR)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER, INTENT(IN)  :: MAX_IR
      REAL   , INTENT(IN)  :: rhoval
      INTEGER, INTENT(OUT) :: nrings,rings(MAX_IR)

      INTEGER ik,ir,ir1,i1,i2
      LOGICAL addtolist

c...  Find which ring the particular value of RHO corresponds to:
      DO ir = 2, irwall
        IF (rho(ir,CELL1).EQ.0.0) CYCLE
c...    Take the first case since RHO in some secondary PFR could catch things up:
        IF (rhoval.GE.rho(ir,IN14).AND.rhoval.LT.rho(ir,OUT23)) EXIT  
      ENDDO
      IF (ir.EQ.irwall+1) 
     .  CALL ER('SelectGridRegion_DIVIMP','Ring not identified',*99)

      nrings = 1
      rings(nrings) = ir

c...  Assemble a list of immediate neighbours:
      DO ik = 1, nks(ir)
        DO i1 = 1, 2
          IF (i1.EQ.1) ir1 = irins (ik,ir)
          IF (i1.EQ.2) ir1 = irouts(ik,ir)
          addtolist = .TRUE.
          DO i2 = 1, nrings
            IF (rings(i2).EQ.ir1) addtolist = .FALSE.
          ENDDO
          IF (addtolist) THEN
            nrings = nrings + 1
            IF (nrings.GT.MAX_IR) 
     .        CALL ER('SelectGridRegion_DIVIMP','IR bound exceeded',*99)
            rings(nrings) = ir1
          ENDIF
        ENDDO
      ENDDO

      RETURN
99    WRITE(0,*) '  RHOVAL= ',rhoval
      WRITE(0,*) '  RHO   = ',rho(2:irsep,CELL1)
      CALL OutputData(85,'FRUSTRATION...')
      STOP
      END
c
c ======================================================================
c
      REAL FUNCTION FindSeparatrixRadius(mode)   
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER, INTENT(IN) :: mode

      INTEGER FindMidplaneCell

      INTEGER ikm,ik,ir,id
      REAL    rmid,frac

c...  Check if the grid is connected double-null -- lame, should
c     really pass the value of CONNECTED from DIVIMP:
      IF (rs(ikto,irsep).LT.rxp.AND.rs(ikti,irsep).LT.rxp) 
     .  connected = .TRUE.

      IF (connected) THEN
        ir = irsep2
      ELSE
        ir = irsep
      ENDIF

      ikm = FindMidplaneCell(ir)
      id = korpg(ikm,ir)

      frac = (z0 - zvertp(1,id)) / (zvertp(4,id) - zvertp(1,id))

      rmid = (1.0 - frac) * rvertp(1,id) + frac * rvertp(4,id)
c      rmid = MAX(rvertp(1,id),rvertp(4,id))

      write(0,*) 'fucking heell',ikm,id,rmid,ir,frac

      IF (.NOT.connected.AND.rmid.LT.rxp) THEN
c...    Really want the outer midplane radius, but this result
c       suggests the grid is connected, so try again:
        ikm = FindMidplaneCell(irsep2)
        id = korpg(ikm,irsep2)
        rmid = MAX(rvertp(1,id),rvertp(4,id))
        WRITE(0,*) 'WHOA! Looks like a connected grid but CONNECTED '//
     .             'not set'
        STOP 'HALTING CODE'
      ENDIF

      FindSeparatrixRadius = rmid

      RETURN
99    STOP
      END
c
c ======================================================================
c
      INTEGER FUNCTION FindMidplaneCell(ir)   
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'

      INTEGER, INTENT(IN) :: ir

      INTEGER ikm,ik,id
      REAL    rmid,z(2)

c...  Find midplane cell:   ! *** NEED TO DO THIS SO THAT PHI IS REFERENCED PROPERLY TO PHI=0 AT THE OUTER MIDPLANE ***
      ikm = -1              !     NEED A BETTER WAY...
      rmid = -100.0
      DO ik = 1, nks(ir)
        id = korpg(ik,ir)
        z(1) = 0.5 * (zvertp(1,id) + zvertp(2,id))
        z(2) = 0.5 * (zvertp(3,id) + zvertp(4,id))
        IF (((z(1).LT.z0.AND.z(2).GE.z0).OR.
     .       (z(2).LT.z0.AND.z(1).GE.z0)).AND.
c        IF (((z(1).LT.0.0.AND.z(2).GE.0.0).OR.
c     .       (z(2).LT.0.0.AND.z(1).GE.0.0)).AND.
     .      rs(ik,ir).GT.rmid) THEN
          ikm  = ik
          rmid = rs(ik,ir)
        ENDIF
      ENDDO

      FindMidplaneCell = ikm

      RETURN
99    STOP
      END
c
c
c
c ======================================================================
c SHOULD NOT BE HERE!  TEMPORARY...
c ======================================================================
c
c *** AM I USING BRRATIO CORRECTLY? DO I NEED TO MAP THE ANGLE SOMEHOW, 
c SO THAT I'M GETTING THE PITCH ANGLE RIGHT IN THE FIELD LINE
c COORDINATE SYSTEM? IS THIS ALL BOGUS? ***
c
c subroutine: TraceFieldLine_DIVIMP
c
      SUBROUTINE TraceFieldLine_DIVIMP(xin,yin,zin,mode,chop,
     .                                 length1,length2,rlimit,
     .                                 n,v,index,fraction,ring,MAXN)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'colours'
      INCLUDE 'slout'
      INCLUDE 'slcom'

      INTEGER, INTENT(IN) :: mode,MAXN,chop
      REAL   , INTENT(IN) :: xin,yin,zin
      REAL*8              :: length1,length2,rlimit
      INTEGER, INTENT(OUT) :: n,index(MAXN),ring
      REAL   , INTENT(OUT) :: fraction(MAXN)
      REAL*8 , INTENT(OUT) :: v(3,MAXN)

      INTEGER FindMidplaneCell,CalcPoint
      LOGICAL InCell
      REAL    FindSeparatrixRadius, ATAN2C

      INTEGER fp,ik,ir,ikm,id,iobj,isur,ipts,i1,ik1,ir1,ike,status
      LOGICAL finished,debug
      REAL    rsep,rin
      REAL*8  r(2),z(2),deltar,deltaz,deltac,deltap,phi,dphi,rval,pval,
     .        frac1,frac2,angle1,angle2,brat,rfrac,r12frac,z12frac,
     .        rhoin,phiin,frac,rpos,rvalmin,zvalmin,zvalmax,
     .        len1,len2,lenmax1,lenmax2,p1(2),p2(2),t,ar,az,br,bz,
     .        dpol,dtor,alpha,origin(3),d1,d2

      REAL       TOL
      PARAMETER (TOL=1.0D-06)

      debug = .TRUE.

      fp = 0

      dphi = 0.5 !  1.0 ! 5.0 ! 10.0  ! Make this an adjustable parameter...

c      IF (yin.NE.0.0) 
c     .  CALL ER('TraceFieldLine_DIVIMP','Sorry, midplane only',*99)

c...  Convert from Cartesean coordinates to RHO and PHI:
      rsep = FindSeparatrixRadius(1)
      rhoin = DBLE(SQRT(xin**2 + zin**2) - rsep)
      phiin = DBLE(ATAN2C(zin,xin) * 180.0 / PI)
      IF (phiin.LT.0.0D0) phiin = phiin + 360.0D0

c      IF (ABS(xin).LT.1.0E-06) THEN
c        IF (zin.GT.0.0) THEN
c          phiin = 90.0
c        ELSE
c          phiin = 270.0
c        ENDIF
c      ELSE
c        phiin = DBLE(ATAN(ABS(zin / xin)) * 180.0 / PI)
c        IF (xin.LT.0.0.and.zin.GT.0.0) phiin = 180.0 - phiin
c        IF (xin.LT.0.0.and.zin.LT.0.0) phiin = 180.0 + phiin
c        IF (xin.GT.0.0.and.zin.LT.0.0) phiin = 360.0 - phiin
c      ENDIF

      IF (debug) THEN
        WRITE(0,*) 'TRACE  MODE,CHOP : ',mode,chop
        WRITE(0,*) 'TRACE  X,Y,ZIN   : ',xin,yin,zin
        WRITE(0,*) 'TRACE  RHO,PHIIN : ',rhoin,phiin
        WRITE(0,*) 'TRACE  RSEP      : ',rsep
        WRITE(0,*) 'TRACE  LENGTH1,2 : ',length1,length2
      ENDIF

      n = 0

      len1 = 0.0D0
      len2 = 0.0D0

      SELECTCASE (chop)
        CASE(1)    ! No restrictions
          rvalmin =  0.0D0
          zvalmin = -1.0D+20
          zvalmax =  1.0D+20
          lenmax1 =  1.0D+20
          lenmax2 =  1.0D+20
        CASE(2:3,6)  ! No field lines above and below or inside the low field side, main plasma separatrix 
          rvalmin =  0.0D0
          zvalmin = -1.0D+20
          zvalmax =  1.0D+20
          ir = irsep-1
          IF (idring(ir).NE.BOUNDARY) THEN
            rvalmin = DBLE(rxp)
            zvalmin =  1.0D+20
            zvalmax = -1.0D+20
            DO ik = 1, nks(ir)-1
              DO i1 = 1, 4
                id = korpg(ik,ir)
                zvalmin = MIN(zvalmin,DBLE(zvertp(i1,id)))
                zvalmax = MAX(zvalmax,DBLE(zvertp(i1,id)))
              ENDDO
            ENDDO
          ENDIF
          IF (chop.EQ.6) THEN
            lenmax1 = length1
            lenmax2 = length2
          ELSE
            lenmax1 = 1.0D+20
            lenmax2 = 1.0D+20
          ENDIF
c          WRITE(0,*) 'ZVALs:',zvalmin,zvalmax
        CASE(4:5)  ! Restricted length...
          rvalmin =  0.0D0
          zvalmin = -1.0D+20
          zvalmax =  1.0D+20
          lenmax1 =  length1
          lenmax2 =  length2
        CASE(7)  ! No trace below or above vertical boundaries, or inside RLIMIT
          rvalmin =  rlimit
          zvalmin =  length1
          zvalmax =  length2
          lenmax1 =  1.0D+20
          lenmax2 =  1.0D+20
        CASE DEFAULT
          CALL ER('TraceFieldLine_DIVIMP','Unrecognised CHOP',*99)
      ENDSELECT
c      WRITE(0,*) 'CHOP:',chop,rvalmin,zvalmin,zvalmax

c...  Find where we are on the outer midplane:
      IF (.TRUE.) THEN
        rin = SQRT(xin**2 + zin**2)
        DO ir = 2, irwall
          IF (idring(ir).EQ.BOUNDARY) CYCLE
          DO ik = 1, nks(ir)
            IF (InCell(ik,ir,rin,yin)) EXIT
          ENDDO
          IF (ik.NE.nks(ir)+1) EXIT
        ENDDO
        IF (ir.EQ.irwall+1) 
     .    CALL ER('TraceFieldLine_DIVIMP','Cell indices not found',*99)

c       Find the distance from the point in question to each of the cell sides:
        id = korpg(ik,ir)
        ar = DBLE(rvertp(1,id))
        az = DBLE(zvertp(1,id))
        br = DBLE(rvertp(4,id))
        bz = DBLE(zvertp(4,id))
        status = CalcPoint(ar,az,br,bz,DBLE(rin),DBLE(yin),t)
        p1(1) = ar + t * (br - ar)
        p1(2) = az + t * (bz - az)

c        WRITE(0,*) ':',t,status

        IF (debug) THEN
          WRITE(0,*) 'IK,IR  :',ik,ir
          WRITE(0,*) 'AR,AZ  :',ar,az
          WRITE(0,*) 'BR,BZ  :',br,bz
          WRITE(0,*) 'P1     :',p1(1:2)
        ENDIF

        ar = DBLE(rvertp(2,id))
        az = DBLE(zvertp(2,id))
        br = DBLE(rvertp(3,id))
        bz = DBLE(zvertp(3,id))
        status = CalcPoint(ar,az,br,bz,DBLE(rin),DBLE(yin),t)
        p2(1) = ar + t * (br - ar)
        p2(2) = az + t * (bz - az)

c        WRITE(0,*) ':',t,status

        IF (debug) THEN
          WRITE(0,*) 'AR,AZ  :',ar,az
          WRITE(0,*) 'BR,BZ  :',br,bz
          WRITE(0,*) 'P2     :',p2(1:2)
          WRITE(0,*) 'RIN,YIN:',rin,yin
        ENDIF

        d1 = DSQRT( (DBLE(rin)-p1(1))**2 + (DBLE(yin)-p1(2))**2 ) 
        d2 = DSQRT( (DBLE(rin)-p2(1))**2 + (DBLE(yin)-p2(2))**2 ) 

        frac = d1 / (d1 + d2)

c        WRITE(0,*) 'D1,2:',d1,d2
c        WRITE(0,*) 'FRAC:',frac

        rhoin = (1.0D0 - frac) * DBLE(rho(ir,SIDE14)) + 
     .                   frac  * DBLE(rho(ir,SIDE23))

c       Assign the symmetry cell for the field-line trace:
        ikm = ik
        IF (t.GT.0.5D0) ikm = ikm + 1  ! Advance a bit along the ring if the point is closer to the far boundary
 
c...    Decide how to assign the magnetic field line pitch angle information:
        SELECTCASE (mode)
          CASE (1)
c...        Just take the pitch angle at the center of the cell all the time:
            ir1 = 0
            r12frac = frac
            z12frac = frac ! 0.5D0
            rfrac = 0.0D0
          CASE (2)
c...        Interpolate the field line pitch angle, gives a continuous b-field:
            r12frac = frac
            z12frac = frac
            IF (rhoin.LT.rho(ir,CELL1)) THEN
              ir1 = irins(ik,ir)
              rfrac =     -(rhoin          - DBLE(rho(ir,CELL1))) / 
     .                 DBLE(rho(ir1,CELL1) -      rho(ir,CELL1) )
            ELSE
              ir1 = irouts(ik,ir)
              rfrac =      (rhoin          - DBLE(rho(ir,CELL1))) / 
     .                 DBLE(rho(ir1,CELL1) -      rho(ir,CELL1) )
            ENDIF
          CASE DEFAULT
            CALL ER('TraceFieldLine_DIVIMP','Unrecognised MODE',*99)
        ENDSELECT

      ELSE

        STOP 'OLD METHOD'

        DO ir = 2, irwall
c...      For now, ignore the inner SOL for a CDN grid (CONNECTED
c         currently assigned/hacked in FindSeparatrixRadius, which
c         has to be called before this check is made, nasty):
          IF (connected.AND.ir.GE.irsep.AND.ir.LT.irsep2) CYCLE
          
          IF (rho(ir,CELL1).EQ.0.0) CYCLE
          
          IF (debug) WRITE(0,*) 'TRACE IR SCAN : ',
     .      ir,SNGL(rhoin),rho(ir,SIDE14),rho(ir,SIDE23)
          
          IF (rhoin.GE.DBLE(rho(ir,SIDE14)).AND.
     .        rhoin.LT.DBLE(rho(ir,SIDE23))) THEN
c...        Find midplane cell:   ! *** NEED TO DO THIS SO THAT PHI IS REFERENCED PROPERLY TO PHI=0 AT THE OUTER MIDPLANE ***
            ikm = -1              !     NEED A BETTER WAY...          
            ikm = FindMidplaneCell(ir)
            IF (ikm.EQ.-1) CALL ER('TraceFieldLine_DIVIMP','No '//
     .                             'midplane cell found',*99)
          
            IF (debug) WRITE(0,*) 'TRACE RING,IKM : ',ir,ikm
            
c...        Decide how to assign the magnetic field line pitch angle information:
            SELECTCASE (mode)
              CASE (1)
c...            Just take the pitch angle at the center of the cell all the time:
                ir1 = 0
                frac = DBLE(rhoin          - rho(ir,SIDE14)) / 
     .                 DBLE(rho(ir,SIDE23) - rho(ir,SIDE14))
                r12frac = frac
                z12frac = frac ! 0.5D0
                rfrac = 0.0D0
              CASE (2)
c...            Interpolate the field line pitch angle, gives a continuous b-field:
                frac = DBLE(rhoin          - rho(ir,SIDE14)) / 
     .                 DBLE(rho(ir,SIDE23) - rho(ir,SIDE14))
                r12frac = frac
                z12frac = frac
                IF (rhoin.LT.rho(ir,CELL1)) THEN
                  ir1 = irins(ikm,ir)
                  rfrac =     -(rhoin          - DBLE(rho(ir,CELL1))) / 
     .                     DBLE(rho(ir1,CELL1) -      rho(ir,CELL1) )
                ELSE
                  ir1 = irouts(ikm,ir)
                  rfrac =      (rhoin          - DBLE(rho(ir,CELL1))) / 
     .                     DBLE(rho(ir1,CELL1) -      rho(ir,CELL1) )
                ENDIF
              CASE DEFAULT
                CALL ER('TraceFieldLine_DIVIMP','Unrecognised MODE',*99)
            ENDSELECT
            EXIT
          ENDIF
        ENDDO
  
      ENDIF

      IF (ir.EQ.irwall+1) 
     .  CALL ER('TraceFieldLine_DIVIMP','Ring not identified',*99)

      ring = ir

      IF (debug) THEN
        WRITE(0,*) '  RING  =',ring,irsep
        WRITE(0,*) '  RFRAC =',rfrac
      ENDIF 

c...  Work from midplane to low IK target:
      phi = phiin
      DO ik = ikm, 1, -1  ! 1, -1
        id = korpg(ik,ir)
        r(1) =          r12frac  * DBLE(rvertp(3,id)) + 
     .         (1.0D0 - r12frac) * DBLE(rvertp(4,id))
        z(1) =          z12frac  * DBLE(zvertp(3,id)) + 
     .         (1.0D0 - z12frac) * DBLE(zvertp(4,id))
        r(2) =          r12frac  * DBLE(rvertp(2,id)) +
     .         (1.0D0 - r12frac) * DBLE(rvertp(1,id))
        z(2) =          z12frac  * DBLE(zvertp(2,id)) + 
     .         (1.0D0 - z12frac) * DBLE(zvertp(1,id))
        rpos = 0.5D0 * (r(1) + r(2))
        IF (rfrac.LE.0.0D0) THEN
          ik1 = ikins(ik,ir)
          ir1 = irins(ik,ir)
c          WRITE(0,*) 'IK,IR,IK1,IR1=',ik,ir,ik1,ir1
          brat = (1.0D0 + rfrac) * DBLE(bratio(ik,ir)  ) - 
     .                    rfrac  * DBLE(bratio(ik1,ir1))
        ELSE
          ik1 = ikouts(ik,ir)
          ir1 = irouts(ik,ir)
          brat = (1.0D0 - rfrac) * DBLE(bratio(ik,ir)  ) + 
     .                    rfrac *  DBLE(bratio(ik1,ir1))

        ENDIF
c        WRITE(0,'(A,5F10.5)') 
c     .    'BRAT:',REAL(brat),REAL(r(1:2)),REAL(z(1:2))
        deltar = r(2) - r(1)
        deltaz = z(2) - z(1)
        dpol   = DSQRT(deltar**2 + deltaz**2)
        alpha  = DASIN(brat)
        dtor   = dpol / DTAN(alpha)
        deltac = dtor
c        deltac = ABS(deltaz) / brat
        deltap = -1.0D0 * deltac / rpos * 180.0D0 / DBLE(PI)
        angle1 = 0.0D0
        finished = .FALSE.
        DO WHILE (angle1.GT.deltap)
c          IF (ik.EQ.ikm) THEN                         ! *** MAKE THIS AN OPTION *** 
c            angle2 = MAX(angle1-0.1D0*dphi,deltap) 
c          ELSE
            angle2 = MAX(angle1-      dphi,deltap) 
c          ENDIF
c          angle2 = MAX(angle1-dphi,deltap) 
c          frac1 = angle1 / deltap
          frac2 = angle2 / deltap
          IF (angle1.EQ.0.0D0.AND.ik.EQ.ikm) THEN
            n = n + 1
            IF (n.GT.MAXN) CALL ER('TraceFieldLine_DIVIMP','N bust',*99)
            v(1,n) = r(1)
            v(2,n) = z(1)
            v(3,n) = phi
            index(n) = ik + 1
            fraction(n) = frac2
c...        Convert from r,z,phi to x,y,z (y okay already):
            rval = v(1,n)
            pval = v(3,n) * DBLE(PI) / 180.0D0
            v(1,n) = rval * DCOS(pval)
            v(3,n) = rval * DSIN(pval)
c            WRITE(0,*) '==START OF FILAMENT:',SNGL(v(3,n)),SNGL(phi)
          ENDIF
c...      Don't follow the field line at all if length restriction set to zero:
          IF ((chop.EQ.4.OR.chop.EQ.5.OR.chop.EQ.6).AND.
     .        lenmax1.EQ.0.0D0) THEN
            finished = .TRUE.
            EXIT
          ENDIF
          n = n + 1
          IF (n.GT.MAXN) CALL ER('TraceFieldLine_DIVIMP','N bust',*99)
          v(1,n) = r(1) + frac2 * deltar
          v(2,n) = z(1) + frac2 * deltaz
          v(3,n) = phi + angle2
          index(n) = ik
          fraction(n) = frac2
c...      Convert from r,z,phi to x,y,z (y okay already):
          rval = v(1,n)
          pval = v(3,n) * DBLE(PI) / 180.0D0
          v(1,n) = rval * DCOS(pval)
          v(3,n) = rval * DSIN(pval)
          len1 = len1 + DSQRT((v(1,n)-v(1,n-1))**2 +
     .                        (v(2,n)-v(2,n-1))**2 +
     .                        (v(3,n)-v(3,n-1))**2)

c          WRITE(0,*) 'ZVAL-:',v(2,n),zvalmax,zxp
c          WRITE(0,*) 'ZVAL1-:',ik,frac2

          IF ((rval  .LE.rvalmin).OR.
     .        (v(2,n).LE.zvalmin.OR.v(2,n).GE.zvalmax).OR.
c     .        (zxp.LT.0.0.AND.v(2,n).LE.zvalmin).OR.
c     .        (zxp.GT.0.0.AND.v(2,n).GE.zvalmax).OR.
     .        (len1  .GE.lenmax1)) THEN
            finished = .TRUE.
            EXIT
          ENDIF
          angle1 = angle2
        ENDDO 
        IF (finished) EXIT
        phi = phi + deltap
      ENDDO

      IF (debug) THEN
        WRITE(0,*) '  IK    =',ik
      ENDIF 

      origin(1:3) = v(1:3,1)

c...  Swap order of these points, so that they start at the low IK target
c     and proceed to the midplane:
      DO i1 = 1, n/2
        v(1:3,n+1   ) = v(1:3,i1    )
        v(1:3,i1    ) = v(1:3,n-i1+1)
        v(1:3,n-i1+1) = v(1:3,n+1   )
        index(n+1   ) = index(i1    )
        index(i1    ) = index(n-i1+1)
        index(n-i1+1) = index(n+1   )
        fraction(n+1   ) = fraction(i1    )
        fraction(i1    ) = fraction(n-i1+1)
        fraction(n-i1+1) = fraction(n+1   )
      ENDDO

c...  Work from midplane to high IK target:
      phi = phiin
      ike = nks(ir)
      IF (ir.LT.irsep) ike = ike - 1
      DO ik = ikm+1, ike
c...    Don't follow the field line at all if length restriction set to zero:
        IF ((chop.EQ.4.OR.chop.EQ.5.OR.chop.EQ.6).AND.
     .      lenmax2.EQ.0.0D0) EXIT
        id = korpg(ik,ir)
        r(1) =          r12frac  * DBLE(rvertp(2,id)) + 
     .         (1.0D0 - r12frac) * DBLE(rvertp(1,id))
        z(1) =          z12frac  * DBLE(zvertp(2,id)) + 
     .         (1.0D0 - z12frac) * DBLE(zvertp(1,id))
        r(2) =          r12frac  * DBLE(rvertp(3,id)) +
     .         (1.0D0 - r12frac) * DBLE(rvertp(4,id))
        z(2) =          z12frac  * DBLE(zvertp(3,id)) + 
     .         (1.0D0 - z12frac) * DBLE(zvertp(4,id))
        rpos = 0.5D0 * (r(1) + r(2))
        IF (rfrac.LE.0.0D0) THEN
          ik1 = ikins(ik,ir)
          ir1 = irins(ik,ir)
          brat = (1.0D0 + rfrac) * DBLE(bratio(ik,ir  )) -
     .                    rfrac  * DBLE(bratio(ik1,ir1))
        ELSE
          ik1 = ikouts(ik,ir)
          ir1 = irouts(ik,ir)
          brat = (1.0D0 - rfrac) * DBLE(bratio(ik,ir  )) + 
     .                    rfrac  * DBLE(bratio(ik1,ir1))
        ENDIF
        deltar = r(2) - r(1)
        deltaz = z(2) - z(1)
        dpol = DSQRT(deltar**2 + deltaz**2)
        alpha = DASIN(brat)
        dtor = dpol / DTAN(alpha)
        deltac = dtor
        deltap = deltac / rpos * 180.0D0 / DBLE(PI)
        angle1 = 0.0D0
        finished = .FALSE.
        DO WHILE (angle1.LT.deltap)
c          IF (ik.EQ.ikm+1) THEN
c            angle2 = MIN(angle1+0.1D0*dphi,deltap) 
c          ELSE
            angle2 = MIN(angle1+      dphi,deltap) 
c          ENDIF
c          angle2 = MIN(angle1+dphi,deltap) 
c          frac1 = angle1 / deltap
          frac2 = angle2 / deltap
          n = n + 1
          IF (n.GT.MAXN) CALL ER('TraceFieldLine_DIVIMP','N bust',*99)
          v(1,n) = r(1) + frac2 * deltar
          v(2,n) = z(1) + frac2 * deltaz
          v(3,n) = phi + angle2
          index(n) = ik
          fraction(n) = frac2
c...      Convert from r,z,phi to x,y,z (y okay already):
          rval = v(1,n)
          pval = v(3,n) * DBLE(PI) / 180.0D0
          v(1,n) = rval * DCOS(pval)
          v(3,n) = rval * DSIN(pval)
          len2 = len2 + DSQRT((v(1,n)-v(1,n-1))**2 +
     .                        (v(2,n)-v(2,n-1))**2 +
     .                        (v(3,n)-v(3,n-1))**2)
          IF ((rval  .LE.rvalmin).OR.
     .        (v(2,n).LE.zvalmin.OR.v(2,n).GE.zvalmax).OR.
c     .        (zxp.LT.0.0.AND.v(2,n).GE.zvalmax).OR.
c     .        (zxp.GT.0.0.AND.v(2,n).LE.zvalmin).OR.
     .        (len2  .GE.lenmax2 )) THEN
            finished = .TRUE.
            EXIT
          ENDIF

c          WRITE(0,*) 'ZVAL2-:',ik,frac2,ike

          angle1 = angle2
        ENDDO 
        IF (ik.LT.ike) index(n) = ik + 1
        IF (finished) EXIT
        phi = phi + deltap
      ENDDO

      IF (debug) THEN
        WRITE(0,*) '  IK    =',ik,nks(ring)
      ENDIF 

c      DO i1 = 1, n
c        WRITE(0,*) '-->',i1,ike,index(i1)
c      ENDDO

c...  Convert from r,z,phi to x,y,z (y okay already):
c      DO i1 = 1, n
c        rval = v(1,i1)
c        pval = v(3,i1) * DBLE(PI) / 180.0D0
c        v(1,i1) = rval * DSIN(pval)
c        v(3,i1) = rval * DCOS(pval)
cc        WRITE(fp,*) ' V:',v(1:3,i1)
c      ENDDO

c      WRITE(0,*) 'n:',n

c      WRITE(0,*) 'LENGTH:',len1,len2,len1+len2

      SELECTCASE (chop)
        CASE(1)
          length1 = len1
          length2 = len2
        CASE(2)
        CASE(3)
          length1 = len1
          length2 = len2
        CASE(4)
        CASE(5)
          length1 = len1
          length2 = len2
        CASE(6)
        CASE(7)
        CASE DEFAULT
          CALL ER('TraceFieldLine_DIVIMP','Unrecognised CHOP',*99)
      ENDSELECT

c...  Add the origin:
      v(1:3,n+1) = origin(1:3)

c      WRITE(0,* ) '============ORIGIN1',v(1:3,n+1)


      RETURN
 99   WRITE(0,*) ' CHOP= ',chop
      STOP
      END
