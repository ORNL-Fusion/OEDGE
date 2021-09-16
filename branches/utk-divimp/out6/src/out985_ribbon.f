c      
c     -*-Fortran-*-
c
c ======================================================================
c
c subroutine: outCrossSeparatrix
c
      LOGICAL FUNCTION outCrossSeparatrix(ir,p1,p2)
      use mod_params
      use mod_cgeom
      use mod_comtor
      IMPLICIT none
c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'comtor'

      INTEGER, INTENT(IN) :: ir
      REAL*8 , INTENT(IN) :: p1(2),p2(2)

      INTEGER ik,id
      REAL*8  p3(2),p4(2),t1,t2

      outCrossSeparatrix = .FALSE.

      DO ik = 1, nks(ir)
        id = korpg(ik,ir)

        p3(1) = DBLE(rvertp(2,id))
        p3(2) = DBLE(zvertp(2,id))
        p4(1) = DBLE(rvertp(3,id))
        p4(2) = DBLE(zvertp(3,id))

        CALL CalcInter(p1(1),p1(2),p2(1),p2(2),
     .                 p3(1),p3(2),p4(1),p4(2),t1,t2)

        IF (t1.GT.0.0D0.AND.t1.LT.1.0D0.AND.
     .      t2.GT.0.0D0.AND.t2.LT.1.0D0) THEN
          outCrossSeparatrix = .TRUE.
          EXIT
        ENDIF
      ENDDO


      RETURN
 99   STOP
      END
c
c ====================================================================== 
c
      SUBROUTINE rayDeleteTrace(tdat,itrace)
      USE mod_out985
      USE mod_out985_variables
      USE mod_out985_ribbon
      IMPLICIT none

      INTEGER, INTENT(IN   ) :: itrace
      REAL*8 , INTENT(INOUT) :: tdat(6,ntrace)

      INTEGER i,n

      n = trace_i(3,itrace) - trace_i(2,itrace) + 1
    
      CALL rayDeleteTracePoints(tdat,trace_i(2,itrace),n,itrace)

      DO i = itrace, trace_n-1
        trace_i(:,i) = trace_i(:,i+1)
      ENDDO
      trace_n = trace_n - 1

      RETURN
 99   STOP
      END
c
c ====================================================================== 
c
      SUBROUTINE rayDeleteTracePoints(tdat,i,n,itrace)
      USE mod_out985
      USE mod_out985_variables
      USE mod_out985_ribbon
      IMPLICIT none

      INTEGER, INTENT(IN   ) :: i,n,itrace
      REAL*8 , INTENT(INOUT) :: tdat(6,ntrace)

      INTEGER j     

      DO j = i, ntrace-n
        tdat (1,j) = tdat (1,j+n)
        tdat (3,j) = tdat (3,j+n)
        tdat (4,j) = tdat (4,j+n)
        tdat (5,j) = tdat (5,j+n)
        tdat (6,j) = tdat (6,j+n)
        trace(:,j) = trace(:,j+n)
      ENDDO
      ntrace = ntrace - n

      DO j = itrace, trace_n
        IF (trace_i(2,j).GT.i) trace_i(2,j) = trace_i(2,j) - n
        IF (trace_i(3,j).GT.i) trace_i(3,j) = trace_i(3,j) - n
        IF (trace_i(4,j).GT.i) trace_i(4,j) = trace_i(4,j) - n
      ENDDO

      RETURN
 99   STOP
      END
c
c ====================================================================== 
c
      SUBROUTINE rayProcessRibbon(tdat,scale,clean)
      USE mod_out985
      USE mod_out985_variables
      USE mod_out985_ribbon
      IMPLICIT none

      REAL*8  DATAN2C
      LOGICAL outCrossSeparatrix

      REAL, PARAMETER :: PI = 3.141592

      LOGICAL, INTENT(IN) :: clean
      REAL*8, INTENT(INOUT) :: tdat(6,ntrace)
      REAL*8, INTENT(IN   ) :: scale

      INTEGER itrace,i1,i2,io,i,j,k,n,ir,irib
      LOGICAL debug,cont
      REAL*8  diff,s_ref,min_diff,p1(2),p2(2),save_r1

      REAL*8, ALLOCATABLE :: rho_local(:)


      debug = .TRUE.


      DO itrace = 1, trace_n
        i1 = trace_i(2,itrace)
        i2 = trace_i(3,itrace)
        io = trace_i(4,itrace) 
c        io = trace_i(4,itrace) + i1 - 1

c...    Calculate 's' and set CODE:
c
c       0 - tangency point
c       1 - transition from the vacuum in the structure
c       2 - transition from the structure into the vacuum vessel
c       3 - inside the vacuum
c       4 - inside the structure

        IF (trace(4,io).EQ.1.0D0.OR.trace(4,io).EQ.2.0D0) THEN
          tdat(5,io) = trace(4,io) + 2.0
        ELSE
          CALL ER('rayGenerateRibbonGrid','Origin point location'//
     .            'not known',*99)
        ENDIF

c...    First half of the ring:
        DO i = io-1, i1, -1
          tdat(1,i) = tdat(1,i+1) - 
     .                DSQRT( (trace(1,i)-trace(1,i+1))**2 + 
     .                       (trace(2,i)-trace(2,i+1))**2 + 
     .                       (trace(3,i)-trace(3,i+1))**2 ) 
c         Assign CODE variable:
          IF (trace(4,i).EQ.3.0D0) THEN
            IF (tdat(5,i+1).EQ.2.0D0.OR.tdat(5,i+1).EQ.3.0D0) THEN
              tdat(5,i) = 1.0D0
            ELSE    
              tdat(5,i) = 2.0D0
            ENDIF                
          ELSE
            IF     (tdat(5,i+1).EQ.1.0D0) THEN
              tdat(5,i) = 4.0D0
            ELSEIF (tdat(5,i+1).EQ.2.0D0) THEN
              tdat(5,i) = 3.0D0
            ELSE
              tdat(5,i) = tdat(5,i+1)
            ENDIF
          ENDIF
        ENDDO

c...    Second half:
        DO i = io+1, i2
          tdat(1,i) = tdat(1,i-1) + 
     .                DSQRT( (trace(1,i)-trace(1,i-1))**2 + 
     .                       (trace(2,i)-trace(2,i-1))**2 + 
     .                       (trace(3,i)-trace(3,i-1))**2 ) 
c         Assign CODE variable:
          IF (trace(4,i).EQ.3.0D0) THEN
            IF (tdat(5,i-1).EQ.2.0D0.OR.tdat(5,i-1).EQ.3.0D0) THEN
              tdat(5,i) = 1.0D0
            ELSE    
              tdat(5,i) = 2.0D0
            ENDIF                
          ELSE
            IF     (tdat(5,i-1).EQ.1.0D0) THEN
              tdat(5,i) = 4.0D0
            ELSEIF (tdat(5,i-1).EQ.2.0D0) THEN
              tdat(5,i) = 3.0D0
            ELSE
              tdat(5,i) = tdat(5,i-1)
            ENDIF
          ENDIF
        ENDDO

c...    For debugging:
        DO i = i1, i2
          tdat(3,i) = DATAN2C(trace(3,i),trace(1,i)) * 180.0D0/DBLE(PI)
        ENDDO
        DO i = i1+1, i2
          tdat(4,i) = tdat(3,i) - tdat(3,i-1)
        ENDDO

c...    Check for duplicate points and remove any that are found.  I don't
c       really know, at the moment, why these might be here, but need to 
c       get rid of them:
        cont = .TRUE.
        DO WHILE (cont)
          cont = .FALSE.
          i1 = trace_i(2,itrace)
          i2 = trace_i(3,itrace)
          DO i = i1+1, i2
            IF (DABS(tdat(1,i)-tdat(1,i-1)).LT.1.0D-6) THEN
              cont = .TRUE.
              n = 1
              DO j = i+1, i2
                IF (DABS(tdat(1,j)-tdat(1,i)).LT.1.0D-6) n = n + 1
              ENDDO        
c              IF (debug) WRITE(0,'(A,4I8,2X,F10.6,2X,3F10.6)') 
c     .              'deleting!',itrace,i,n,ntrace,tdat(1,i),trace(1:3,i)
              CALL rayDeleteTracePoints(tdat,i,n,itrace)
              EXIT
            ENDIF
          ENDDO
        ENDDO

      ENDDO  ! itrace


c.... Need to calcualte and store the connection lengths here:



c...  Crop to boundary:

c      r1 = 
c      r2 = 
c      s1 = 
c      s2 = 

      write(0,*) 'crop:'
      write(0,*) crop_r1,crop_r2,crop_s1,crop_s2

      IF (DABS(crop_r1+1.0D0).GT.1.0D-5) THEN

        ALLOCATE(rho_local(trace_n))
        rho_local =  0.0D0
        save_r1   = -1.0D0
        DO itrace = 1, trace_n                 
          DO i = trace_i(2,itrace), trace_i(3,itrace)
            IF (tdat(1,i).EQ.0.0D0) THEN
              i1 = i
              EXIT
            ENDIF
          ENDDO
          IF (itrace.GT.1) THEN
            rho_local(itrace) = rho_local(itrace-1) + 
     .                  SNGL(DSQRT((trace(1,i1)-trace(1,i2))**2 + 
     .                             (trace(2,i1)-trace(2,i2))**2 + 
     .                             (trace(3,i1)-trace(3,i2))**2))
c...        Check if the line crosses the separatrix:
            IF (crop_r1.LT.0.0D0) THEN
              p1(1) = DSQRT(trace(1,i2)**2 + trace(3,i2)**2)
              p1(2) = trace(2,i2)
              p2(1) = DSQRT(trace(1,i1)**2 + trace(3,i1)**2)
              p2(2) = trace(2,i1)
              ir = NINT(REAL(DABS(crop_r1)))
              IF (outCrossSeparatrix(ir,p1,p2)) THEN
                write(0,*) 'ir=',ir
                write(0,*) 'p1=',p1(1:2)
                write(0,*) 'p2=',p2(1:2)
                write(0,*) 'rho_local',itrace,rho_local(itrace)
                save_r1 = crop_r1
                crop_r1 = rho_local(itrace-1)
                write(0,*) 'setting separatrix crop',crop_r1
                write(6,*) 'setting separatrix crop',crop_r1
              ENDIF
            ENDIF
          ENDIF
          i2 = i1
        ENDDO

        DO itrace = trace_n, 1, -1
          IF (rho_local(itrace).LT.crop_r1.OR.
     .        rho_local(itrace).GT.crop_r2) THEN
c            write(0,*) 'deleting itrace rho',
c     .        rho_local(itrace),crop_r1,crop_r2
            DO i = itrace, trace_n-1
              rho_local(i) = rho_local(i+1)
            ENDDO
            CALL rayDeleteTrace(tdat,itrace)
          ENDIF
        ENDDO

c...    Restore CROP_R1 if necessary:
        IF (save_r1.NE.-1.0D0) crop_r1 = save_r1

        DO itrace = 1, trace_n
          i1 = -1
          i2 = -1
          DO i = trace_i(2,itrace), trace_i(3,itrace)
            IF (i1.EQ.-1.AND.-tdat(1,i).LT.crop_s1) i1 = i
            IF (             -tdat(1,i).LT.crop_s1) i2 = i
          ENDDO
c          IF (i1.NE.-1) write(0,*) 'h crop',itrace,i1,i2
          IF (i1.NE.-1) 
     .      CALL rayDeleteTracePoints(tdat,i1,i2-i1+1,itrace)
        ENDDO

        DO itrace = 1, trace_n
          i1 = -1
          i2 = -1
          DO i = trace_i(2,itrace), trace_i(3,itrace)
            IF (i1.EQ.-1.AND.-tdat(1,i).GT.crop_s2) i1 = i
            IF (             -tdat(1,i).GT.crop_s2) i2 = i
          ENDDO
c          IF (i1.NE.-1) write(0,*) 'h crop',itrace,i1,i2
          IF (i1.NE.-1) 
     .      CALL rayDeleteTracePoints(tdat,i1,i2-i1+1,itrace)
        ENDDO

        DEALLOCATE(rho_local)

      ENDIF


c      stop


      IF (clean) THEN

c...    Crop the end(s) first trace to nearest intersection(s):

c...    First half of the trace (-ve S values):
        itrace = 1
        i1 = trace_i(2,itrace)
        io = trace_i(4,itrace) 
        i2 = trace_i(3,itrace)
        DO i = io-1, i1, -1
          IF (tdat(5,i).EQ.1.0D0.OR.tdat(5,i).EQ.2.0D0) EXIT
        ENDDO
        IF (i.GT.i1) THEN
          WRITE(0,*) 'crop1',itrace
          WRITE(0,*) '     ',i1,io,i2       
          WRITE(0,*) '     ',i
          WRITE(0,*) '     ',i1,i-i1,itrace
          CALL rayDeleteTracePoints(tdat,i1,i-i1,itrace)
          WRITE(0,*) '     ',trace_i(2,itrace),trace_i(4,itrace),
     .                       trace_i(3,itrace)
          tdat(5,i1) = -1.0D0
c          tdat(5,i1) = -tdat(5,i1) 
          s_ref = tdat(1,i1)
c...      Scan down the rest of the traces and crop them as well, hopefully
c         following the surface contour:
          DO itrace = 2, trace_n
c            write(0,*) 'itrace=',itrace
            i1 = trace_i(2,itrace)
            io = trace_i(4,itrace) 
            DO j = i1, io-2  ! Not allowed to mess with the origin point
              IF (tdat(5,j).NE.1.0D0.AND.tdat(5,j).NE.2.0D0) CYCLE
              diff = tdat(1,j) - s_ref
              IF (diff.GE.0.0D0.AND.diff.LE.scale) THEN
                s_ref = tdat(1,j)
                CALL rayDeleteTracePoints(tdat,i1,j-i1,itrace)
                tdat(5,i1) = -1.0D0
c                tdat(5,i1) = -tdat(5,i1) 
c                write(0,*) '     cut 1',j-i1+1,trace_i(4,itrace)-i1+1
c                write(0,*) '          ',tdat(5,i1)
                EXIT
              ENDIF
            ENDDO
            IF (j.EQ.io-1) THEN
c...          No intersection point found, so just clear out the trace
c             based on S only: 
              k = -1
              min_diff = scale * 10.0D0
              DO j = i1, io-2
                diff = DABS(tdat(1,j) - s_ref)
                IF (diff.LT.min_diff) THEN
                  min_diff = diff
                  k = j
                ENDIF
              ENDDO
              IF (k.EQ.-1) THEN
                CALL ER('rayProcessRibbon','No valid link found 1',*99)
              ELSE
c                write(0,*) '     cut 2',k-i1+1,trace_i(4,itrace)-i1+1
                s_ref = tdat(1,k)
                CALL rayDeleteTracePoints(tdat,i1,k-i1,itrace)
                tdat(5,i1) = -1.0D0
c                write(0,*) '          ',tdat(5,i1)
              ENDIF
            ENDIF
          ENDDO 
        ENDIF
c...    Last half of the trace (+ve S values):
        itrace = 1
        i1 = trace_i(2,itrace)
        io = trace_i(4,itrace) 
        i2 = trace_i(3,itrace)
        DO i = io+1, i2
          IF (tdat(5,i).EQ.1.0D0.OR.tdat(5,i).EQ.2.0D0) EXIT
        ENDDO
        IF (i.LT.i2) THEN
          WRITE(0,*) 'crop2',itrace
          WRITE(0,*) '     ',i1,io,i2       
          WRITE(0,*) '     ',i
          WRITE(0,*) '     ',i2,i2-i,itrace
          CALL rayDeleteTracePoints(tdat,i+1,i2-i,itrace)
          tdat(5,i) = -2.0D0
c          tdat(5,i) = -tdat(5,i) 
          s_ref = tdat(1,i)
          DO itrace = 2, trace_n
            io = trace_i(4,itrace) 
            i2 = trace_i(3,itrace)
            DO j = io+2, i2
              IF (tdat(5,j).NE.1.0D0.AND.tdat(5,j).NE.2.0D0) CYCLE
              diff = s_ref - tdat(1,j)
              IF (diff.GE.0.0D0.AND.diff.LE.scale) THEN
                s_ref = tdat(1,j)
                CALL rayDeleteTracePoints(tdat,j+1,i2-j,itrace)
c                tdat(5,j) = -tdat(5,j) 
                tdat(5,j) = -2.0D0
                EXIT
              ENDIF
            ENDDO
            IF (j.EQ.i2+1) THEN
              k = -1
              min_diff = scale * 10.0D0
              DO j = io+2, i2
                diff = DABS(tdat(1,j) - s_ref)
                IF (diff.LT.min_diff) THEN
                  min_diff = diff
                  k = j
                ENDIF
              ENDDO
              IF (k.EQ.-1) THEN
                CALL ER('rayProcessRibbon','No valid link found 2',*99)
              ELSE
                s_ref = tdat(1,k)
                CALL rayDeleteTracePoints(tdat,k+1,i2-k,itrace)
                tdat(5,k) = -2.0D0
              ENDIF
            ENDIF
          ENDDO 
        ENDIF

c        STOP 'fsdfds'

c...    Wipe portions of the intersection data
        IF (wipe_n.GT.0) THEN
          WRITE(0,*) 'list',wipe_list(1:wipe_n)

          ALLOCATE(rho_local(trace_n))
          rho_local =  0.0D0
          DO itrace = 1, trace_n                 
            DO i = trace_i(2,itrace), trace_i(3,itrace)
              IF (tdat(1,i).EQ.0.0D0) THEN
                i1 = i
                EXIT
              ENDIF
            ENDDO
            IF (itrace.GT.1) 
     .        rho_local(itrace) = rho_local(itrace-1) + 
     .                    SNGL(DSQRT((trace(1,i1)-trace(1,i2))**2 + 
     .                               (trace(2,i1)-trace(2,i2))**2 + 
     .                               (trace(3,i1)-trace(3,i2))**2))
            i2 = i1
          ENDDO

          DO j = 1, wipe_n
            irib = wipe_list(j)
            WRITE(0,*) 'params',opt%rib_r1(irib),opt%rib_r2(irib),
     .                          opt%rib_s1(irib),opt%rib_s2(irib)
 
            DO itrace = 2, trace_n
              i1 = trace_i(2,itrace) 
              i2 = trace_i(3,itrace)
              write(0,*) 'rho',itrace,rho_local(itrace)
              write(0,*) 'i12',i1,i2,ntrace
              write(0,*) 'srg',tdat(1,i1),tdat(1,i2)
              DO i = i1, i2
                IF (rho_local(itrace).GE.opt%rib_r1(irib).AND.
     .              rho_local(itrace).LE.opt%rib_r2(irib).AND.
     .              tdat(1,i)        .GE.opt%rib_s1(irib).AND.
     .              tdat(1,i)        .LE.opt%rib_s2(irib)) THEN
                  write(0,*) 'found',tdat(5,i)
                  tdat(5,i) = 3.0D0
                ENDIF
              ENDDO
            ENDDO
          ENDDO
c         stop 'the promised land'
          DEALLOCATE(rho_local)
        ENDIF


      ENDIF  ! clean.EQ..TRUE.



      RETURN
 99   STOP
      END
c
c ====================================================================== 
c
      SUBROUTINE rayBuildBump(tdat,ntangent,itangent,scale,status)
      USE mod_out985
      USE mod_out985_variables
      USE mod_out985_ribbon
      IMPLICIT none
 
      INTEGER, INTENT(IN)    :: ntangent,itangent(2,1000)
      REAL*8 , INTENT(IN   ) :: scale
      REAL*8 , INTENT(INOUT) :: tdat(6,ntrace)
      LOGICAL, INTENT(OUT  ) :: status

      INTEGER itrace,i,it,i1,i2,j,k,dist1_i,dist2_i,count,count1,count2,
     .        fp,dcount
      LOGICAL debug
      REAL*8  dist,dist1,dist2,slope1,slope2,slope,diff,mindiff,
     .        last_dist,bump,frac,check

      status = .TRUE.

      debug = .TRUE.
      fp    = 6

      dist1_i = itangent(2,ntangent)
      dist2_i = dist1_i
      slope1  = 1.0E+6
      slope2  = 1.0E+6

      bump = DBLE(REAL(ntangent))

      tdat(6,dist1_i) = bump

      count = 0


c        IF (tdat(1,dist1_i).GT.29.49D0.AND.
c     .      tdat(1,dist1_i).LT.29.50D0) 
c     .    WRITE(fp,*) 'rockin!'

      IF (debug) 
     .  WRITE(fp,*) '========tangent s',bump,itangent(1,ntangent),
     .              tdat(1,dist1_i)


c      IF (itangent(1,ntangent).EQ.3) debug = .TRUE.

      dcount    = 0
      last_dist = -1.0E+6

      DO itrace = itangent(1,ntangent)+1, trace_n
        count = count + 1

        IF (debug) WRITE(fp,*) 'itrace',itrace

        dist1 = 1.0D+6
        dist2 = 1.0D+6
        i1 = dist1_i
        i2 = dist2_i
        dist1_i = -1
        dist2_i = -1

        slope1  = 1.0E+6
        slope2  = 1.0E+6

        IF (debug) WRITE(fp,*) 'i1,2',i1,i2

        DO i = trace_i(2,itrace), trace_i(3,itrace)

          IF (i.EQ.it) CYCLE      

          IF (tdat(5,i).EQ.1.0D0.OR.tdat(5,i).EQ.2.0D0) THEN

            IF (debug) WRITE(fp,*) 'tdat(5',tdat(5,i), '-------------'
            IF (debug) WRITE(fp,*) 'tdat(1',tdat(1,i),
     .                 tdat(1,i1),tdat(1,i2)
            IF (debug) WRITE(fp,*) 'count ',count

            IF     (tdat(1,i).LT.0.5D0*(tdat(1,i1)+tdat(1,i2))) THEN

              dist = DSQRT( (trace(1,i)-trace(1,i1))**2 + 
     .                      (trace(2,i)-trace(2,i1))**2 + 
     .                      (trace(3,i)-trace(3,i1))**2 ) 
              slope = DABS(tdat(1,i) - tdat(1,i1))
c              slope = (tdat(1,i) - tdat(1,i1))


              IF (debug) WRITE(fp,*) 'dist1 ',dist ,dist1
              IF (debug) WRITE(fp,*) 'slope1',slope,slope1,
     .                   DABS((slope-slope1)/slope1)
              IF (((count.LE.10.AND.slope.LT.      scale).OR.
     .             (count.GT.10.AND.slope.LT.0.5D0*scale)).AND.
c              IF (((count.LE.10.AND.slope.LT.0.10D0).OR.
c     .             (count.GT.10.AND.slope.LT.0.05D0)).AND.
     .            slope.LT.slope1) THEN
c              IF (dist.LT.0.10D0.AND.dist.LT.dist1.AND.
c     .            (DABS((slope-slope1)/slope1).LT.10.0D0.OR.
c     .             count.LE.2.OR.count1.GT.0)) THEN
                IF (debug) WRITE(fp,*) 'ok1',i
                dist1   = dist   
                slope1  = slope
                dist1_i = i
              ENDIF
            ELSEIF (tdat(1,i).GT.0.5D0*(tdat(1,i1)+tdat(1,i2))) THEN
              dist = DSQRT( (trace(1,i)-trace(1,i2))**2 + 
     .                      (trace(2,i)-trace(2,i2))**2 + 
     .                      (trace(3,i)-trace(3,i2))**2 ) 
              slope = DABS(tdat(1,i) - tdat(1,i2))
c              slope = (tdat(1,i) - tdat(1,i2))

              IF (debug) WRITE(fp,*) 'dist2 ',dist ,dist2
              IF (debug) WRITE(fp,*) 'slope2',slope,slope2,
     .                   DABS((slope-slope2)/slope2)
              IF (((count.LE.10.AND.slope.LT.      scale).OR.
     .             (count.GT.10.AND.slope.LT.0.5D0*scale)).AND.
c              IF (((count.LE.10.AND.slope.LT.0.10D0).OR.
c     .             (count.GT.10.AND.slope.LT.0.05D0)).AND.
     .            slope.LT.slope2) THEN
c              IF (dist.LT.0.10D0.AND.dist.LT.dist2.AND.
c     .            (DABS((slope-slope2)/slope2).LT.10.0D0.OR.
c     .             count.LE.2.OR.count2.GT.0)) THEN
                IF (debug) WRITE(fp,*) 'ok2',i
                dist2   = dist
                slope2  = slope
                dist2_i = i
              ENDIF
            ENDIF

          ENDIF

        ENDDO

c...    Valid intersection point not found, so pick one:
        IF (dist1_i.EQ.-1.AND.dist2_i.NE.-1) THEN 
          count1 = count1 + 1
          mindiff = 10.0D0 * scale
c          mindiff = 1.0D0
          DO i = trace_i(2,itrace), dist2_i-1
            IF (tdat(5,i).LT.-0.000001D0) CYCLE
c          DO i = trace_i(2,itrace), trace_i(3,itrace)
c            IF (tdat(5,i).LT.-0.000001D0.OR.i.EQ.dist2_i) CYCLE
            slope = DABS(tdat(1,i) - tdat(1,i1))
c            slope = (tdat(1,i) - tdat(1,i1))
            IF (.FALSE..AND.count1.EQ.1) THEN
              diff = DABS(slope - slope1)
            ELSE
              diff = DABS(slope)
            ENDIF
            IF (diff.LT.mindiff) THEN
              IF (debug) THEN
                WRITE(fp,*) 'ok1, just picking one',bump,itrace
                WRITE(fp,*) '                     ',i,i2,count1
                WRITE(fp,*) '                     ',tdat(5,i)
                WRITE(fp,*) '                     ',tdat(1,i),tdat(1,i2)
                WRITE(fp,*) '                     ',diff,mindiff
              ENDIF
              mindiff = diff
              dist1_i = i
            ENDIF 
          ENDDO
        ELSE
          count1 = 0
        ENDIF
 
        IF (dist2_i.EQ.-1.AND.dist1_i.NE.-1) THEN 
          count2 = count2 + 1
          mindiff = 10.0D0 * scale
c          mindiff = 1.0D0
          DO i = dist1_i+1, trace_i(3,itrace)
            IF (tdat(5,i).LT.-0.000001D0.OR.i.EQ.dist1_i) CYCLE
c          DO i = trace_i(2,itrace), trace_i(3,itrace)
c            IF (tdat(5,i).LT.-0.000001D0.OR.i.EQ.dist1_i) CYCLE
            slope = DABS(tdat(1,i) - tdat(1,i2))
c            slope = (tdat(1,i) - tdat(1,i2))
            IF (.FALSE..AND.count2.EQ.1) THEN
              diff = DABS(slope - slope2)
            ELSE
              diff = DABS(slope)
            ENDIF
            IF (diff.LT.mindiff) THEN
              IF (debug) THEN
                WRITE(fp,*) 'ok2, just picking one',bump,itrace
                WRITE(fp,*) '                     ',i,i2,count2
                WRITE(fp,*) '                     ',tdat(5,i)
                WRITE(fp,*) '                     ',tdat(1,i),tdat(1,i2)
                WRITE(fp,*) '                     ',diff,mindiff
              ENDIF
              mindiff = diff
              dist2_i = i
            ENDIF 
          ENDDO
        ELSE
          count2 = 0
        ENDIF

c...    Stop the bump if no valid points found:
        IF (dist1_i.EQ.-1.OR.dist2_i.EQ.-1) EXIT
c        IF (dist1_i.EQ.-1.AND.dist2_i.EQ.-1) EXIT

c...    Stop the bump if it is getting narrower when moving out radially:
c        dist = DABS(tdat(1,dist1_i) - tdat(1,dist2_i))
c        write(6,*) 'checking dist',dist,last_dist
c        IF (last_dist.GT.0.0D0.AND.dist-last_dist.LT.-0.001D0) THEN
c          write(6,*) 'dcount!',dcount
c          dcount = dcount + 1
c          IF (dcount.EQ.10) EXIT 
c        ELSE
c          dcount = 0
c        ENDIF
c        last_dist = dist

c...    For the selected points, set them negative and set all points
c       between them negative as well:
c        tdat(5,dist1_i:dist2_i) = -4.0D0

        tdat(5,dist1_i) = -2.0D0
        tdat(5,dist2_i) = -1.0D0
        IF (count1.NE.0) tdat(5,dist1_i) = -20.0D0
        IF (count2.NE.0) tdat(5,dist2_i) = -10.0D0

        tdat(6,dist1_i:dist2_i) = bump

      ENDDO  ! itrace


c...  Fill in any holes in the bump:      
      DO itrace = itangent(1,ntangent)+1, trace_n

        DO i = trace_i(2,itrace), trace_i(3,itrace)
          check = 0.0D0
          IF (tdat(5,i).EQ.-20.0D0) check = -20.0D0
          IF (tdat(5,i).EQ.-10.0D0) check = -10.0D0
          IF (check.EQ.0.0D0) CYCLE
          IF (debug) THEN
            write(fp,*) 'replace',bump,check
            write(fp,*) '       ',itangent(1,ntangent),itrace,
     .                            trace_n
          ENDIF
          i1 = -1
          i2 = -1
c...
          DO k = itrace-1, itangent(1,ntangent), -1
            IF (debug) write(fp,*) '  trying',k
            DO j = trace_i(2,k), trace_i(3,k)          
              IF (tdat(6,j).EQ.bump.AND.(tdat(5,j).EQ.check/10.0D0.OR.
     .                                   tdat(5,j).EQ.      0.0D0)) THEN
                i1 = j
                EXIT
              ENDIF
            ENDDO  
            IF (i1.NE.-1) EXIT
          ENDDO
c... 
          DO k = itrace+1, trace_n
            IF (debug) write(fp,*) '  trying',k
            DO j = trace_i(2,k), trace_i(3,k)          
              IF (tdat(6,j).EQ.bump.AND.(tdat(5,j).EQ.check/10.0D0.OR.
     .                                   tdat(5,j).EQ.      0.0D0)) THEN
                i2 = j
                EXIT
              ENDIF
            ENDDO  
            IF (i2.NE.-1) EXIT
          ENDDO
          IF (i1.EQ.-1.AND.i2.EQ.-1) THEN
            WRITE(0,*) 'WARNING rayBuildBump: Reference points '//
     .                 'not found'
            WRITE(0,*) '       ',i1,i2
            WRITE(0,*) '       ',itrace,bump
            trace(5,i) = -0.999999D0
            IF (debug) THEN
              WRITE(fp,*) 'WARNING rayBuildBump: Reference points '//
     .                   'not found'
              WRITE(fp,*) '       ',i1,i2
              WRITE(fp,*) '       ',itrace,bump
            ENDIF
          ELSEIF (i1.NE.-1.AND.i2.EQ.-1) THEN
            WRITE(fp,*) 'WARNING rayBuildBump: Second reference '//
     .                  'point not found'
            WRITE(fp,*) '       ',i1,i2
            WRITE(fp,*) '       ',itrace,bump
            trace(5,i) = trace(5,i1)
          ELSEIF (i1.EQ.-1.AND.i2.NE.-1) THEN
            WRITE(fp,*) 'WARNING rayBuildBump: First reference '//
     .                  'point not found'
            WRITE(fp,*) '       ',i1,i2
            WRITE(fp,*) '       ',itrace,bump
            trace(5,i) = trace(5,i2)
          ELSE
            frac = 1.0D0 / DBLE(REAL(k - (itrace - 1)))
            tdat (1  ,i) = (1.0D0 - frac) * tdat (1  ,i1) +  ! s, distance along the trace
     .                              frac  * tdat (1  ,i2)  
            trace(1:3,i) = (1.0D0 - frac) * trace(1:3,i1) +  ! x,y,z
     .                              frac  * trace(1:3,i2)
            trace(5  ,i) = (1.0D0 - frac) * trace(5  ,i1) +  ! impact angle on the surface
     .                              frac  * trace(5  ,i2)
            IF (debug) THEN
              WRITE(fp,*) '       ',i1,i2
              WRITE(fp,*) '       ',itrace,itrace-1,k,frac
              IF (trace(5,i).EQ.0.0) WRITE(fp,*) '  --> watch it!'
            ENDIF
          ENDIF
        ENDDO
      ENDDO

c...  Reset code:
      DO i = 1, ntrace
        IF (tdat(5,i).LE.-10.0D0) tdat(5,i) = tdat(5,i) / 10.0D0  
      ENDDO

c...  Reassign bump markers:
c      write(0,*) 'count',count,tdat(1,itangent(2,ntangent))
 
      DO itrace = itangent(1,ntangent)+1, trace_n
        i1 = -1
        i2 = -1
        DO i = trace_i(2,itrace), trace_i(3,itrace)
          IF (tdat(5,i).EQ.-2.0) i1 = i
          IF (tdat(5,i).EQ.-1.0) i2 = i
          IF (i1.NE.-1.AND.i2.NE.-1) THEN
            IF (count.LE.3) THEN
              status = .FALSE.
              IF ((tdat(6,i1).EQ.0.0D0.OR.tdat(6,i1).EQ.bump).OR.
     .            (tdat(6,i2).EQ.0.0D0.OR.tdat(6,i1).EQ.bump)) THEN
                IF (debug) THEN
                  write(fp,*) 'removing false tangency point',itrace
                  write(fp,*) '   itrace=',itrace
                  write(fp,*) '    i1,i2=',i1,i2
                ENDIF
                tdat(5,i1:i2) = 3.0D0
                tdat(6,i1:i2) = 0.0D0
              ENDIF
            ELSE
              IF (debug) THEN
                write(fp,*) 'setting -ve markers'
                write(fp,*) '   itrace=',itrace
                write(fp,*) '    i1,i2=',i1,i2
              ENDIF
              tdat(5,i1:i2) = -4.0D0
              tdat(5,i1   ) = -2.0D0
              tdat(5,i2   ) = -1.0D0
            ENDIF
            EXIT
c            tdat(6,i1:i2) = DBLE(REAL(ntangent))
c            WRITE(0,*) 'i1,2',i1,i2
          ENDIF
        ENDDO
      ENDDO

      IF (debug) WRITE(fp,*) 'count at end',count
c      IF (count.EQ.1) STOP 'got ya!' 

      RETURN
 99   STOP
      END
c
c ====================================================================== 
c
      LOGICAL FUNCTION rayCheckClearance(itrace,ipoint,tdat,scale)
      USE mod_out985
      USE mod_out985_variables
      USE mod_out985_ribbon
      IMPLICIT none
 
      INTEGER, INTENT(IN) :: itrace,ipoint
      REAL*8 , INTENT(IN) :: tdat(6,ntrace)
      REAL*8 , INTENT(IN) :: scale

      INTEGER fp,i1,i2,i3,i
      REAL*8  dist1,dist2

      rayCheckClearance = .FALSE.

      fp = 6

      write(fp,*) 'checkclearance',itrace,ipoint,tdat(1,ipoint)

      DO i1 = itrace-1, 1, -1
        dist1 = 1.0D+6
        dist2 = 1.0D+6
        i2 = trace_i(2,i1)
        i3 = trace_i(3,i1)
        DO i = i2, i3
c          write(fp,*) 'ipopint',ipoint
c          write(fp,*) 'i      ',i,i2,i3
c          write(fp,*) '       ',tdat(1,i)
c          write(fp,*) '       ',tdat(1,ipoint)
c          write(fp,*) '       ',tdat(5,i)
c          write(fp,*) 'what:',i,ipoint,tdat(1,i),tdat(1,ipoint),
c     .                        tdat(5,i)
          IF (tdat(1,i).LE.tdat(1,ipoint).AND.tdat(5,i).LT.0.0D0)
     .      dist1 = MIN(dist1,tdat(1,ipoint)-tdat(1,i))
          IF (tdat(1,i).GE.tdat(1,ipoint).AND.tdat(5,i).LT.0.0D0)
     .      dist2 = MIN(dist2,tdat(1,i)-tdat(1,ipoint))
        ENDDO
        WRITE(fp,*) 'tan check',i1,SNGL(dist1),SNGL(dist2),
     .     SNGL(dist1+dist2)
        IF (dist1+dist2.LT.5.0D0*scale) THEN
          write(fp,*) 'KILLER!',dist1,dist2,tdat(1,ipoint)
          rayCheckClearance = .TRUE.
          EXIT
        ENDIF
      ENDDO

      RETURN
 99   STOP
      END
c
c ====================================================================== 
c
      SUBROUTINE rayCleanRibbon(tdat,scale)
      USE mod_out985
      USE mod_out985_variables
      USE mod_out985_ribbon
      IMPLICIT none
 
      REAL*8 , INTENT(INOUT) :: tdat(6,ntrace)
      REAL*8 , INTENT(IN   ) :: scale

      LOGICAL rayCheckClearance

      INTEGER fp,ntangent,i,j,i1,i2,i3,itrace,n,
     .        itangent(2,1000)  ! 1=itrace,2=index
      LOGICAL debug,status

      fp = 6

      ntangent = 0

      DO itrace = 1, trace_n-1
        i1 = trace_i(2,itrace)
        i2 = trace_i(3,itrace)

        IF (fp.NE.-1) write(fp,*) 'checking',itrace,i1,i2,ntrace
        
c...    Look for new tangency points:
        i3 = -1
        DO i = i1, i2
          IF (i.GT.trace_i(3,itrace)) EXIT

          IF     ((tdat(5,i).EQ.1.0D0.OR.
     .             tdat(5,i).EQ.2.0D0).AND.i3.EQ.-1) THEN  
            i3 = i                                       
          ELSEIF ((tdat(5,i).EQ.1.0D0.OR.
     .             tdat(5,i).EQ.2.0D0).AND.i3.NE.-1) THEN  ! reversed (compared to the output file).
            IF (fp.NE.-1) write(fp,*) 'candidate',i3,i

            IF (DABS(tdat(1,i)-tdat(1,i3)).GT.scale) THEN
c            IF (DABS(tdat(1,i)-tdat(1,i3)).GT.0.10D0) THEN
c            IF (DABS(tdat(1,i)-tdat(1,i3)).GT.0.05D0) THEN  
c             The points are too far apart to be a tangency pair:
              IF (fp.NE.-1) write(fp,*) 'too far apart'
              i3 = i
            ELSEIF (rayCheckClearance(itrace,i3,tdat,scale)) THEN
c...          Make sure the space above the tangency point is legit:
            ELSE
              IF (fp.NE.-1) write(fp,*) 'new point!'
c             Generate a new tangency point:
              ntangent = ntangent + 1
              itangent(1,ntangent) = itrace
              itangent(2,ntangent) = i3
              tdat(1,i3) = 0.5D0 * (tdat(1,i3) + tdat(1,i))
              tdat(5,i3) = 0.0D0
              DO j = 1, 3
                trace(j,i3) = 0.5D0 * (trace(j,i3) + trace(j,i))
              ENDDO
              trace(5,i3) = 0.0D0
c             Delete points spanning tangency pair (inclusive):
              n = i - i3
              IF (fp.NE.-1) write(fp,*) 'new point - n',n
              CALL rayDeleteTracePoints(tdat,i3+1,n,itrace)
              CALL rayBuildBump(tdat,ntangent,itangent,scale,status)
              IF (.NOT.status) THEN
c...            Something was invalid about the tangency point, so delete it:
                tdat(5,i3) = 3.0D0
                tdat(6,i3) = 0.0D0
                ntangent = ntangent - 1
                IF (fp.NE.-1) write(fp,*) '   *** DELETING POINT ***'
              ENDIF
              i3 = -1          
            ENDIF
          ENDIF
        ENDDO
      ENDDO

c      stop 'sdfsd'

c        DO i = 1, ntrace
c          IF (tdat(5,i).EQ.0.0D0) STOP 'AH HA 2!'
c        ENDDO

c...  Assign all points with a +ve CODE value to 3:
      DO i = 1, ntrace
        IF (tdat(5,i).GT.0.0D0) tdat(5,i) = 3.0D0
        IF (tdat(5,i).LT.0.0D0) tdat(5,i) = -tdat(5,i)
      ENDDO

c        DO i = 1, ntrace
c          IF (tdat(5,i).EQ.0.0D0) STOP 'AH HA 3!'
c        ENDDO

      RETURN
 99   STOP
      END
c
c ====================================================================== 
c
      SUBROUTINE rayStoreRibbon(irib,tdat,clean)
      USE mod_out985
      USE mod_out985_variables
      USE mod_out985_ribbon
      IMPLICIT none
 
      INTEGER, INTENT(IN   ) :: irib
      LOGICAL, INTENT(IN   ) :: clean
      REAL*8 , INTENT(INOUT) :: tdat(6,ntrace)

      INTEGER i1,i2,i3,i4,i5,i6,iloop,fp,itrace,i,npts,code,bump
      LOGICAL inside
      REAL    rho_local,angle

      CHARACTER fname*512,sufx*64

      fp = 99

      DO iloop = 1, 2  !2
c...    Store data:

        SELECTCASE (iloop)
          CASE(1)
            sufx = '_000_full'
          CASE(2)
            sufx = '_000'
        ENDSELECT


        WRITE(fname,'(A)') 'trace.'//TRIM(opt%rib_tag(irib))//TRIM(sufx)
        OPEN(UNIT=fp,FILE=TRIM(fname),STATUS='REPLACE',
     .       ACCESS='SEQUENTIAL')

        WRITE(fp,'(A)') '* Ribbon grid data file from RAY, which, '//
     .                  'like, totally rocks'
        WRITE(fp,'(A)') '*'
        WRITE(fp,'(A)') '{VERSION}'
        WRITE(fp,'(F5.1)') 2.1

        WRITE(fp,'(A)') '*'
        WRITE(fp,'(A)') '*  CODE:'
        WRITE(fp,'(A)') '*'
        WRITE(fp,'(A)') '*    0 - tangency point'
        WRITE(fp,'(A)') '*    1 - wall intersection point, where '//
     .                           'the field line goes from inside '//
     .                           'the vacuum vessel volume to outside'
        WRITE(fp,'(A)') '*    2 - wall intersection point, where '//
     .                           'the field line goes from outside '//
     .                           'the vacuum vessel volume to inside'
        WRITE(fp,'(A)') '*    3 - point is inside the vacuum vessel '//
     .                  'volume'
        WRITE(fp,'(A)') '*    4 - point is outside'
        WRITE(fp,'(A)') '*    5 - end point for the field line trace'
        WRITE(fp,'(A)') '*'
        WRITE(fp,'(A)') '* Field lines where the origin is on the '//
     .                  'LFS run from the bottom of the machine to '//
     .                  'the top.'
        WRITE(fp,'(A)') '*'

        WRITE(fp,'(A  )') '{TRACE SUMMARY}'
        WRITE(fp,'(A  )') '*'
        WRITE(fp,'(A  )') '* number of traces'
        WRITE(fp,'(I10)') trace_n
        WRITE(fp,'(A)') '*'
        WRITE(fp,'(A)') '*  INDEX     - field line trace number'
        WRITE(fp,'(A)') '*  RHO_REL   - cross-field coordinate for '//
     .    'the ribbon grid, with the origin on the inner most surface'
        WRITE(fp,'(A)') '*  RHO_ABS   - distance from the separatri'//
     .    'x mapped to the outer midplane, from the true magnetic '//
     .    'equilibrium'
        WRITE(fp,'(A)') '*  ORIGIN    - index of origin point on '//
     .    'each field line (only useful for the _full data file)'
        WRITE(fp,'(A)') '*  X,Y,Z     - location of the origin '//
     .    'point in Cartesian coordinates'
        WRITE(fp,'(A)') '*  #TANGENT  - number of tangency points'
        WRITE(fp,'(A)') '*  #INTER_LO - number of intersection '//
     .    'points for S < 0'
        WRITE(fp,'(A)') '*  #INTER_HI - number of intersection '//
     .    'points for S > 0'
        WRITE(fp,'(A)') '*'

        WRITE(fp,'(A,A8,2A12,2A8,3A12,2X,3A12)')
     .    '*','index','rho_rel','rho_abs','origin','code','x','y','z',
     .    '#tangent','#inter_lo','#inter_hi'
        WRITE(fp,'(A,8X,2A12,16X,3A12,2X)')
     .    '*','(m)','(m)','(m)','(m)','(m)'
        rho_local = 0.0
        DO itrace = 1, trace_n                 
          i1 =  0  ! index of origin
          i3 =  0  ! # tangency
          i5 =  0  ! # s < 0 intersection points
          i6 =  0  ! # s > 0 intersection points
          DO i = trace_i(2,itrace), trace_i(3,itrace)
            IF (tdat(1,i).EQ.0.0D0) THEN
              i1 = trace_i(3,itrace) - i + 1
              i2 = i
            ENDIF
            IF (tdat(5,i).EQ.0.0D0) i3 = i3 + 1
            IF (tdat(1,i).GT.0.0D0) THEN
              IF (tdat(5,i).EQ.1.0D0.OR.tdat(5,i).EQ.2.0D0) i5 = i5 + 1
            ELSE
              IF (tdat(5,i).EQ.1.0D0.OR.tdat(5,i).EQ.2.0D0) i6 = i6 + 1
            ENDIF
          ENDDO
          IF (itrace.GT.1) 
     .      rho_local = rho_local + 
     .                  SNGL(DSQRT( (trace(1,i2)-trace(1,i4))**2 + 
     .                              (trace(2,i2)-trace(2,i4))**2 + 
     .                              (trace(3,i2)-trace(3,i4))**2 ))
          WRITE(fp,'(1X,I8,2F12.7,I8,I8,3F12.7,2X,3I12)')
     .      itrace,rho_local,-1.0,i1,
     .      NINT(REAL(trace(4,i2)))+2,
     .      trace(1:3,i2),i3,i5,i6

          i4 = i2
        ENDDO

        WRITE(fp,'(A)') '{TRACE DATA}'
        WRITE(fp,'(A)') '*'
        WRITE(fp,'(A)') '*  BUMP   - index of bump-like feature '//
     .    'in RHO,S space (not a great description, I know)'
        WRITE(fp,'(A)') '*  S      - distance along the field '//
     .    'the trace, with S=0 at the origin point'
        WRITE(fp,'(A)') '*  IMPACT - impact angle at the '//
     .    'surface for each intersection point'
        WRITE(fp,'(A)') '*  WIDTH  - for radial grid morphing'
        WRITE(fp,'(A)') '*  METRIC - for radial grid morphing'
        WRITE(fp,'(A)') '*'

        DO itrace = 1, trace_n                 
          i1 = trace_i(2,itrace)
          i2 = trace_i(3,itrace)

          IF (tdat(5,i1).NE.1.D0.AND.tdat(5,i1).NE.2.D0) tdat(5,i1)=5.D0
          IF (tdat(5,i2).NE.1.D0.AND.tdat(5,i2).NE.2.D0) tdat(5,i2)=5.D0

          IF (iloop.EQ.1) THEN
            npts = i2 - i1 + 1
          ELSE
            npts = 0
            DO i = i1, i2
              IF (NINT(REAL(tdat(5,i))).LE.2.OR.
     .            NINT(REAL(tdat(5,i))).EQ.5) npts = npts + 1
            ENDDO
          ENDIF

c          IF (npts.EQ.0) CYCLE

          WRITE(fp,'(A ,2A8)') '*','trace','npts'
          WRITE(fp,'(1X,2I8,70X,A,2I8,A)') itrace,npts,
     .       '(',i1,i2,')'
          WRITE(fp,'(A ,A8,2(2X,A6),2X,A12,2X,3A12,2X,A12,2X,2A12,
     .               2X,A12)') 
     .       '*','index','code','bump','s','x','y','z','impact',
     .       'width','metric','dphi'
          WRITE(fp,'(A ,26X,A12,2X,3A12,2X,A12,2X,24X,2X,A12)')
     .       '*','(m)','(m)','(m)','(m)','(degrees)','(degrees)'

          inside = .FALSE.  ! Inside the wall

          DO i = i2, i1, -1  ! In reverse order, so that field lines on the outboard side go from the bottom to top
            code  = NINT(REAL(tdat(5,i)))
            bump  = NINT(REAL(tdat(6,i)))
            angle = SNGL(trace(5,i))

            IF (iloop.EQ.2.AND.NINT(REAL(tdat(5,i))).GT.2.AND.
     .                         NINT(REAL(tdat(5,i))).NE.5) CYCLE
            
            WRITE(fp,'(1X,I8,2(2X,I6),2X,F12.7,2X,3F12.7,2X,
     .                 F12.7,2X,2F12.7,2X,2F8.3,2X,I8,I3)')
     .        i2-i+1,
     .        code,bump,
     .        -tdat(1,i),
     .        trace(1:3,i),
     .        angle,
     .        -1.0,-1.0,
     .        SNGL(tdat(3:4,i)),
     .        i,NINT(REAL(trace(4,i)))

c...        Check if the CODE sequencing makes sense:
            IF (clean.AND.iloop.EQ.1.AND.i.NE.i1.AND.i.NE.i2.AND.
     .          (.NOT.inside.AND.code.EQ.2.OR.
     .                inside.AND.code.EQ.1)) THEN
              WRITE(0,*) 'WARNING rayStoreRibbon: Invalid sequence '//
     .                   'of intersection points'
              WRITE(0,*) '  ITRACE=',itrace
              WRITE(0,*) '  INDEX =',i2-i+1
              WRITE(0,*) '  BUMP  =',bump
            ENDIF
            IF (.NOT.inside.AND.code.EQ.1) inside = .TRUE.
            IF (     inside.AND.code.EQ.2) inside = .FALSE.
c...        Check if the BUMP sequencing makes sense:

c...        Check if the intersection angles make sense:
            IF (clean.AND.iloop.EQ.1.AND.
     .          (angle.LT.0.0.OR.angle.GT.90.0.OR.
     .           (angle.EQ.0.0.AND.(code.EQ.1.OR.code.EQ.2)))
     .         ) THEN
              WRITE(0,*) 'WARNING rayStoreRibbon: Invalid angle '//
     .                   'of intersection'
              WRITE(0,*) '  ???',angle.LT.0.0,angle.GT.90.0,
     .           angle.EQ.0.0,(code.EQ.1.OR.code.EQ.2)
              WRITE(0,*) '  ITRACE=',itrace
              WRITE(0,*) '  INDEX =',i2-i+1
              WRITE(0,*) '  BUMP  =',bump
              WRITE(0,*) '  CODE  =',code
              WRITE(0,*) '  ANGLE =',angle
            ENDIF
          ENDDO

        ENDDO  ! itrace

        CLOSE(fp)

      ENDDO

      RETURN
 99   STOP
      END
c
c ====================================================================== 
c
      SUBROUTINE rayLoadCastemTrace(irib,irad,MAX_NTRACE,file,n,v)
      IMPLICIT none
 
      INTEGER  , INTENT(IN ) :: irib,irad,MAX_NTRACE
      CHARACTER, INTENT(IN ) :: file*(*)
      INTEGER  , INTENT(OUT) :: n
      REAL*8   , INTENT(OUT) :: v(3,MAX_NTRACE)

      INTEGER   fp,ctrace,i,num1,num2
      LOGICAL   skip_origin
      REAL      rdum1
      REAL*8    origin(3)
      CHARACTER dummy*1024

      skip_origin = .TRUE.

      fp   = 99
      IF (irad.EQ.1) WRITE(0,*) 'FILE:',TRIM(file)

      OPEN(fp,FILE=TRIM(file),FORM='FORMATTED',STATUS='OLD',ERR=98)       

      ctrace = 0
      n      = 0
      DO WHILE (.TRUE.)
        READ(fp,'(A1024)',ERR=98) dummy
c        WRITE(0,*) 'dum ',TRIM(dummy),' ',LEN_TRIM(dummy)
        IF (dummy(1:1).EQ.'*'.OR.LEN_TRIM(dummy).LE.1) CYCLE
        READ(dummy,*) num1
        READ(fp   ,*) num2
        EXIT
      ENDDO

      IF (irad.EQ.1) WRITE(0,*) ' num1,2:',num1,num2

      DO WHILE (.TRUE.)
        READ(fp,'(A)',END=10,ERR=98) dummy
        IF (dummy(1:1).EQ.'*'.OR.LEN_TRIM(dummy).LT.5) CYCLE
c        IF (n.EQ.0) WRITE(0,*) 'dummy>',TRIM(dummy)//'<'
        n = n + 1
        READ(dummy,*) rdum1,rdum1,v(1,n),v(3,n),v(2,n)  ! Some reordering to account for the different coordinate systems in CASTEM and RAY
        IF (n.EQ.num2+1.AND.skip_origin) THEN
c        IF (n.EQ.341.AND.skip_origin) THEN
          n = n - 1
          skip_origin = .FALSE.
          CYCLE
        ENDIF
        IF (n.EQ.num2) THEN                              
c        IF (n.EQ.340) THEN                              
c...      Swap the order of the points:
          origin(1:3) = v(1:3,1)
          DO i = 1, n/2
            v(1:3,n  +1) = v(1:3,  i  )
            v(1:3,  i  ) = v(1:3,n-i+1)
            v(1:3,n-i+1) = v(1:3,n  +1)
          ENDDO
        ENDIF
        IF (n.EQ.num1+num2-1) THEN
c        IF (n.EQ.340+590-1) THEN
          ctrace = ctrace + 1
c          write(0,*) 'ctrace',ctrace,irib,irad
          IF (ctrace.EQ.irad) THEN
c          IF (ctrace.EQ.irib) THEN
            EXIT
          ELSE
c            WRITE(0,*) 'n=',n
            n = 0
            skip_origin = .TRUE.
          ENDIF
        ENDIF
      ENDDO          

 10   CONTINUE

      IF (ctrace.NE.irad)
c      IF (ctrace.NE.irib)
     .  CALL ER('rayLoadCastemTrace','Trace not found',*99)
 
c...  And swap again:
      DO i = 1, n/2
        v(1:3,n  +1) = v(1:3,  i  )
        v(1:3,  i  ) = v(1:3,n-i+1)
        v(1:3,n-i+1) = v(1:3,n  +1)
      ENDDO

      v(1:3,n+1) = origin(1:3)

      CLOSE(fp)

      RETURN
 98   CALL ER('rayLoadCastemTrace','Data file not found',*99)
 99   WRITE(0,*) '  FILE= ',TRIM(file)
      STOP
      END
c
c ====================================================================== 
c
      SUBROUTINE rayGenerateRibbonGrid
      USE mod_out985
      USE mod_out985_variables
      USE mod_out985_ribbon
      IMPLICIT none 


      REAL   FindSeparatrixRadius

      INTEGER, PARAMETER :: MAX_NLIST = 100
      REAL   , PARAMETER :: PI        = 3.141592

      INTEGER irib,iphi,irad,i1,i2,i3,i4,fp,i,j,o,n,itrace
      LOGICAL clean
      REAL    xin,yin,zin,x1,y1,z1,x2,y2,z2,phi1,phi2,phi,
     .        frac,frac_point,length,frac_step,frac_mark
      REAL*8  mat(3,3)


      INTEGER n_list
      REAL*8  v_list(3,MAX_NLIST),d_list(MAX_NLIST)

      INTEGER ninter,iobj,iside,isrf
      REAL*8  vinter(3,MAXINTER),dinter(MAXINTER),DTOL,v1(3),v2(3)

      REAL*8, ALLOCATABLE :: tdat(:,:)

      DTOL = 1.0D-10

      WRITE(0,*) 'n     =',opt%rib_n
      WRITE(0,*) 'option=',opt%rib_option(1)
      WRITE(0,*) 'nrad  =',opt%rib_nrad(1)
      WRITE(0,*) 'nphi  =',opt%rib_nphi(1)
      WRITE(0,*) 'limit =',opt%rib_limit(1)
      WRITE(0,*) 'tag   =',TRIM(opt%rib_tag(1))

      crop_r1 = -1.0D0
      crop_r2 = -1.0D0
      crop_s1 = -1.0D0
      crop_s2 = -1.0D0

      wipe_n = 0

      DO irib = 1, opt%rib_n

c...    Set active region of interset for cropping:
        IF (opt%rib_option(irib).EQ.2) THEN
          crop_r1 = DBLE(opt%rib_r1(irib))
          crop_r2 = DBLE(opt%rib_r2(irib))
          crop_s1 = DBLE(opt%rib_s1(irib))
          crop_s2 = DBLE(opt%rib_s2(irib))
          CYCLE
        ENDIF

c...    Set regions to wipe clean of all intersection points:
        IF (opt%rib_option(irib).EQ.3) THEN
          SELECTCASE (opt%rib_wipe(irib))
            CASE (-1) 
              wipe_n = 0
            CASE ( 1,2) 
              wipe_n = wipe_n + 1
              wipe_list(wipe_n) = irib
            CASE DEFAULT
              CALL ER('rayGenerateRibbonGrid','Unrecognized '//
     .              'wipe option',*99)
          ENDSELECT
        ENDIF

        x1 = opt%rib_r(1,irib)
        y1 = opt%rib_z(1,irib)
        x2 = opt%rib_r(2,irib)
        y2 = opt%rib_z(2,irib) 
        z1 = 0.0
        z2 = 0.0

        phi1 = opt%rib_phi(1,irib)
        phi2 = opt%rib_phi(2,irib)

        DO iphi = 1, opt%rib_nphi(irib)
c...      Setup storage for the field line traces:
          IF (ALLOCATED(trace)) DEALLOCATE(trace)
          CALL AllocTrace(-1,MP_INITIALIZE)
c...      Count up the total number of traces there will be:      
          IF (ALLOCATED(trace_i)) DEALLOCATE(trace_i)
          ALLOCATE(trace_i(4,opt%rib_nrad(irib)))
          trace_n = 0

c         Need to rotate the points according to phi:
          frac = REAL(iphi - 1) / REAL(MAX(1,opt%rib_nphi(iphi)-1))
          phi = phi1 + frac * (phi2 - phi1)

          SELECTCASE (opt%rib_trace(irib))
c           --------------------------------------------------------------
            CASE(1)
c...          Collect the field line data (including wall intersection):          
              v1 = ( / x1, y1, z1 / )
              v2 = ( / x2, y2, z2 / )

c             Rotate around the y-axis:
              CALL Calc_Transform2(mat,0.0D0,1,0)
              CALL Calc_Transform2(mat,DBLE(phi*PI/180.0),2,1)
              CALL Transform_Vect (mat,v1)
              CALL Transform_Vect (mat,v2)

              n_list = 0
              d_list = 1.0D+6
              
              WRITE(0,*) 'phi:',phi
              WRITE(0,*) 'v1 :',v1
              WRITE(0,*) 'v2 :',v2
              
              DO iobj = 1, nobj
                DO iside = 1, obj(iobj)%nside
                  DO isrf = obj(iobj)%iside(iside,1), 
     .                      obj(iobj)%iside(iside,2) 
              
                    ninter = 0
                    CALL LineThroughSurface(v1,v2,iobj,iside,isrf,
     .                                      ninter,vinter,dinter,0,DTOL)
              
                    IF     (ninter.GT.1) THEN
                      CALL ER('rayGenerateRibbonGrid','More than one '//
     .                        'intersection detected',*99)
                    ELSEIF (ninter.EQ.1) THEN 
                      n_list = n_list + 1
                      IF (n_list.GT.MAX_NLIST) 
     .                  CALL ER('rayGenerateRibbonGrid','N_LIST '//
     .                          'too large, increase array size',*99)
                      d_list(    n_list) = dinter(    1)
                      v_list(1:3,n_list) = vinter(1:3,1)
                    ENDIF
              
                  ENDDO
                ENDDO
              ENDDO

              WRITE(0,*) 'N_LIST:',n_list
              WRITE(0,*) '      :',d_list(1:n_list)

              i1 = MINLOC(d_list(1:n_list),1)           
c              WRITE(0,*) '  chosen one',i1

              length = SQRT( (v1(1)-v2(1))**2 + (v1(2)-v2(2))**2 +
     .                       (v1(3)-v2(3))**2 )   

              frac_point = d_list(i1) / length
c           --------------------------------------------------------------
            CASE DEFAULT
              frac_point = 1.0
c           --------------------------------------------------------------
          ENDSELECT

c         STOP

          frac_step = 0.1
          frac_mark = frac_step

          DO irad = 1, opt%rib_nrad(irib)

            frac = REAL(irad - 1) / REAL(MAX(1,opt%rib_nrad(irib)-1))

            xin = v1(1) + frac * (v2(1) - v1(1))
            yin = v1(2) + frac * (v2(2) - v1(2))
            zin = v1(3) + frac * (v2(3) - v1(3))
          
            IF (frac.GE.frac_mark) THEN
              WRITE(0,*) 'frac:',frac,frac_point
              frac_mark = frac_mark + frac_step
            ENDIF

            CALL rayProcessTrace(irib,irad,xin,yin,zin,phi,
     .                           frac.LT.frac_point)

          ENDDO  ! irad

          ALLOCATE(tdat(6,ntrace))
          tdat = 0.0D0

          clean = .TRUE.

          CALL rayProcessRibbon(tdat,DBLE(opt%rib_scale(irib)),clean)

          IF (clean) THEN
            CALL rayCleanRibbon(tdat,DBLE(opt%rib_scale(irib)))
          ELSE
            DO i = 1, ntrace
              IF (tdat(5,i).LT.0.0D0) tdat(5,i) = -tdat(5,i)
            ENDDO
          ENDIF

          CALL rayStoreRibbon(irib,tdat,clean)

          DEALLOCATE(tdat)

        ENDDO  ! iphi

      ENDDO  ! irib



      RETURN
      STOP 'WHA-WHO! RIBBINNESS!'

      RETURN
 99   STOP
      END
c
c ====================================================================== 
c

      SUBROUTINE rayProcessTrace(irib,irad,xin,yin,zin,phi,inside)
      USE mod_out985
      USE mod_out985_variables
      USE mod_out985_ribbon
      IMPLICIT none 

      INTEGER, INTENT(IN) :: irib,irad
      LOGICAL, INTENT(IN) :: inside
      REAL   , INTENT(IN) :: xin,yin,zin,phi

      REAL   ATAN2C
      REAL*8 DATAN2C

      INTEGER, PARAMETER :: MAX_NTRACE = 200000 , MAX_NLIST = 100
      REAL   , PARAMETER :: PI         = 3.141592

      INTEGER iv,i1,i2,i3,i4,i5,i6,nv,ring,index(MAX_NTRACE),ishift,
     .        iobj,iside,isrf,ivtx,status,nr,nz,nphi,nplace,i,j,
     .        r1i,z1i,phi1i,
     .        r2i,z2i,phi2i,hold_r1i,hold_z1i,hold_phi1i,
     .        v_dat(3,MAX_NTRACE/100),n_dat,
     .        n_list,s_list(MAX_NLIST)
      LOGICAL cont,check_origin
      REAL    fraction(MAX_NTRACE),r1,z1,phi1,rhoin,
     .        r_min,z_min,r_max,z_max,phi3
      CHARACTER file*512

      REAL*8  mat(3,3)
      REAL*8  len1,len2,rlimit,v(3,MAX_NTRACE),norm(3),vect(3),
     .        length,dprod,angle,
     .        v_list(3,MAX_NLIST),d_list(MAX_NLIST),origin(3)
      REAL*8, ALLOCATABLE :: vdat(:,:)

      INTEGER, ALLOCATABLE :: srf_object(:),srf_check(:,:,:,:)
      LOGICAL, ALLOCATABLE :: srf_mask(:)

      REAL*8  r,z

      INTEGER ninter
      REAL*8  vinter(3,MAXINTER),dinter(MAXINTER)
      REAL*8  DTOL,v1(3),v2(3)

      SAVE

      DTOL = 1.0D-10

      ninter = 0
      n_dat  = 0
      v_dat  = 0

c...  Do the field line trace:
      SELECTCASE (opt%rib_trace(irib))
c       ----------------------------------------------------------------
        CASE(1)  ! Lame DIVIMP field line tracing
          len1   = -3.5D0
          len2   =  6.0D0  ! 5.0D0
          rlimit =  4.6D0
          v_dat = 0
          CALL TraceFieldLine_DIVIMP(xin,yin,zin,2,7,len1,len2,rlimit,
     .                            nv,v,index,fraction,ring,MAX_NTRACE)
          origin(1:3) = v(1:3,nv+1)

          IF (.TRUE.) THEN

c...        Find the point on the trace that's closest to the original origin
c           point in R,Z, and then rotate the field line so that this matches 
c           up with the origin point.  This is necessary because the DIVIMP
c           field line tracing routine is a bit crude:
            
            ALLOCATE(vdat(2,nv))
            
c            WRITE(0,*) 'dat:',nv
            
            vdat(1,1:nv) = DSQRT( v(1,1:nv)**2 + v(3,1:nv)**2 )  ! RHO
            rhoin = SQRT(xin**2 + zin**2)
            vdat(2,1:nv) = DSQRT( (vdat(1,1:nv) - DBLE(rhoin))**2 +  ! distance to origin point in R,Z
     .                            (v   (2,1:nv) - DBLE(yin  ))**2)
            
            i = MINLOC(vdat(2,1:nv),1)           
            
            phi3 = SNGL(DATAN2C(v(3,i),v(1,i))) * 180.0 / PI
            IF (phi-phi3.LT. 360.0) phi3 = phi3 + 360.0
            IF (phi-phi3.GT.-360.0) phi3 = phi3 - 360.0
            
c            WRITE(0,*) 'dat:',xin,yin,zin
c            WRITE(0,*) 'dat:',tdat(2,i-1:i+1)
c            WRITE(0,*) 'dat:',phi,phi3
            
c            WRITE(0,*) 'dat------------:',SNGL(v(1:3,i))
            
c           Rotate around the y-axis:
            CALL Calc_Transform2(mat,0.0D0,1,0)
            CALL Calc_Transform2(mat,DBLE((phi-phi3)*PI/180.0),2,1)
            DO j = 1, nv
              CALL Transform_Vect(mat,v(1:3,j))
            ENDDO
            
c            WRITE(0,*) 'dat------------:',xin,yin,zin
c            WRITE(0,*) 'dat------------:',SNGL(v(1:3,i))
            
            origin(1:3) = v(1:3,i)
            
            DEALLOCATE(vdat)

          ENDIF
c       ----------------------------------------------------------------
        CASE(2)  ! Load data from CASTEM
          file = opt%rib_tfile(irib)
          CALL rayLoadCastemTrace(irib,irad,MAX_NTRACE,file,nv,v)
          origin(1:3) = v(1:3,nv+1)
c       ----------------------------------------------------------------

        CASE DEFAULT
          CALL ER('rayGenerateRibbonGrid','Unrecognized '//
     .            'field line tracing option',*99)
      ENDSELECT            


      IF (.NOT.ALLOCATED(srf_check)) THEN
        nr   = 20
        nz   = 40
        nphi = 60
        nplace = MAX(10000,nsrf / ((nr * nz * nphi) / 1000 + 1))
        ALLOCATE(srf_check (nr,nz,nphi,0:nplace))
        ALLOCATE(srf_object(nsrf))
        ALLOCATE(srf_mask  (nsrf))

        r_min =   3.0
        r_max =   9.0
        z_min =  -4.0
        z_max =   6.5

c        WRITE(0,*) 'rmin,max=',r_min,r_max
c        WRITE(0,*) 'zmin,max=',z_min,z_max

        cont = .TRUE.
        DO WHILE (cont)
          cont = .FALSE.

          DO isrf = 1, nsrf
            DO ivtx = 1, srf(isrf)%nvtx
              v1 = vtx(1:3,srf(isrf)%ivtx(ivtx))           
 
              r1 = SQRT(SNGL(v1(1))**2 + SNGL(v1(3))**2)
              z1 = SNGL(v1(2))
              phi1 = ATAN2C(SNGL(v1(3)),SNGL(v1(1))) * 180. / PI + 180.0

              r1i   = INT(REAL(nr  ) * (r1-r_min) / (r_max-r_min)) + 1
              z1i   = INT(REAL(nz  ) * (z1-z_min) / (z_max-z_min)) + 1
              phi1i = INT(REAL(nphi) * ((phi1+0.00001) / 360.001)) + 1 

              WRITE(6,*) 'r1i,z1i,phi1i=',r1i,z1i,phi1i,phi1

              IF (ivtx.EQ.1) THEN
                hold_r1i   = r1i  
                hold_z1i   = z1i  
                hold_phi1i = phi1i
              ELSEIF (r1i.NE.hold_r1i.OR.phi1i.NE.hold_phi1i.OR.
     .                z1i.NE.hold_z1i) THEN
c...            Exit early from the loop if a point is found that's not
c               located in the same zone as the first point:
                EXIT
              ENDIF

            ENDDO  ! IVTX
            
            IF (ivtx.EQ.srf(isrf)%nvtx+1) THEN
c...          The surface is completely inside the zone:
              i1 = srf_check(r1i,z1i,phi1i,0) + 1
              IF (i1.GT.nplace)
     .          CALL ER('rayProcessTrace','NPLACE is too '//
     .                  'small, be less aggressive',*99)
              srf_check(r1i,z1i,phi1i,0 ) = i1
              srf_check(r1i,z1i,phi1i,i1) = isrf
            ELSE
c...          Location is ill-defined, so add the surface to all of the 
c             surrounding zones:
              DO r1i = hold_r1i-1, hold_r1i+1
                DO z1i = hold_z1i-1, hold_z1i+1
                  DO phi1i = hold_phi1i-1, hold_phi1i+1
                    i1 = r1i
                    IF (i1.EQ.0   ) i1 = nr
                    IF (i1.EQ.nr+1) i1 = 1
                    i2 = z1i
                    IF (i2.EQ.0   ) i2 = nz
                    IF (i2.EQ.nz+1) i2 = 1
                    i3 = phi1i
                    IF (i3.EQ.0     ) i3 = nphi
                    IF (i3.EQ.nphi+1) i3 = 1
                    i4 = srf_check(i1,i2,i3,0) + 1
                    IF (i4.GT.nplace)
     .                CALL ER('rayProcessTrace','NPLACE is too '//
     .                        'small, be less aggressive',*99)
                    srf_check(i1,i2,i3,0 ) = i4
                    srf_check(i1,i2,i3,i4) = isrf          
                  ENDDO
                ENDDO
              ENDDO
            ENDIF

          ENDDO  ! ISRF
        ENDDO  ! CONT loop

c        DO i1 = 1, nr
c          DO i2 = 1, nz
c            DO i3 = 1, nphi
c              WRITE(6,*) 'check:',i1,i2,i3,srf_check(i1,i2,i3,0)
c            ENDDO
c          ENDDO
c        ENDDO

c...    Map all surfaces to objects (assumes that all surfaces are associated
c       with ISIDE=1, which probably won't be true in general):
        DO isrf = 1, nsrf
          DO iobj = 1, nobj
            IF (isrf.GE.obj(iobj)%iside(1,1).AND.
     .          isrf.LE.obj(iobj)%iside(1,2)) THEN
              srf_object(isrf) = iobj
              EXIT
            ENDIF
          ENDDO
        ENDDO

      ENDIF


c...  Find wall intersections:
      DO iv = nv, 2, -1
        v1 = v(1:3,iv  )
        v2 = v(1:3,iv-1)

c      DO iv = 1, nv-1 
c       v1 = v(1:3,iv  )
c       v2 = v(1:3,iv+1)

        n_list = 0

        srf_mask = .FALSE.

        r1   = SQRT(SNGL(v1(1))**2 + SNGL(v1(3))**2)
        z1   = SNGL(v1(2))
        phi1 = ATAN2C(SNGL(v1(3)),SNGL(v1(1))) * 180.0 / PI + 180.0

c             WRITE(0,*) 'segment',r1,z1,phi1

        r1i   = INT(REAL(nr  ) * (r1-r_min) / (r_max-r_min)) + 1
        z1i   = INT(REAL(nz  ) * (z1-z_min) / (z_max-z_min)) + 1
        phi1i = INT(REAL(nphi) * ((phi1+0.00001) / 360.001)) + 1 

        r1   = SQRT(SNGL(v2(1))**2 + SNGL(v2(3))**2)
        z1   = SNGL(v2(2))
        phi1 = ATAN2C(SNGL(v2(3)),SNGL(v2(1))) * 180.0 / PI + 180.0

c              WRITE(0,*) '       ',r1,z1,phi1

        r2i   = INT(REAL(nr  ) * (r1-r_min) / (r_max-r_min)) + 1
        z2i   = INT(REAL(nz  ) * (z1-z_min) / (z_max-z_min)) + 1
        phi2i = INT(REAL(nphi) * ((phi1+0.00001) / 360.001)) + 1 

        IF (r1i.EQ.r2i.AND.z1i.EQ.z2i.AND.phi1i.EQ.phi2i) THEN
          ishift = 0 
        ELSE
          ishift = 1
        ENDIF

c        WRITE(0,*) 'segment',i1,r1,z1,phi1,INT(phi1/10)

        DO r2i = r1i-ishift, r1i+ishift
          DO z2i = z1i-ishift, z1i+ishift
            DO phi2i = phi1i-ishift, phi1i+ishift
              i1 = r2i
              i2 = z2i
              i3 = phi2i

c              WRITE(0,*) '       ',r1i,z1i,phi1i
c              WRITE(0,*) '       ',i1,i2,i3,ishift

              IF (i1.EQ.0     ) i1 = nr
              IF (i1.EQ.nr+1  ) i1 = 1
              IF (i2.EQ.0     ) i2 = nz
              IF (i2.EQ.nz+1  ) i2 = 1
              IF (i3.EQ.0     ) i3 = nphi
              IF (i3.EQ.nphi+1) i3 = 1

c              WRITE(0,*) 'points ',v1
c              WRITE(0,*) '       ',v2
c              WRITE(0,*) '       ',nr,nz,nphi
              

              IF (i1.LT.1.OR.i1.GT.nr.OR.
     .            i2.LT.1.OR.i2.GT.nz.OR.
     .            i3.LT.1.OR.i3.GT.nphi)
     .          CALL ER('rayProcessTrace','Vertex out-of-bounds, '//
     .                  'increas size of domain',*99)
              
              DO i4 = 1, srf_check(i1,i2,i3,0)

                isrf  = srf_check(i1,i2,i3,i4)

                IF (srf_mask(isrf)) CYCLE
             
                iobj  = srf_object(isrf)
                iside = 1

c                WRITE(0,*) '       ',isrf,iobj

                ninter = 0
                CALL LineThroughSurface(v1,v2,iobj,iside,isrf,
     .                                  ninter,vinter,dinter,0,DTOL)

c                WRITE(0,'(A,4I6,2(2X,3F8.3),I6)') 
c     .                '   :',i1,iobj,iside,isrf,v1,v2,n

                IF     (ninter.GT.1) THEN
                  CALL ER('rayProcessTrace','More than one '//
     .                    'intersection detected',*99)
                ELSEIF (ninter.EQ.1) THEN 
                  n_list = n_list + 1
                  IF (n_list.GT.MAX_NLIST) 
     .              CALL ER('rayProcessTrace','N_LIST '//
     .                      'too large, increase array size',*99)
                  d_list(    n_list) = dinter(    1)
                  v_list(1:3,n_list) = vinter(1:3,1)
                  s_list(    n_list) = isrf         
                ENDIF

c                IF (ninter.GE.1) WRITE(0,*) '  cut!'
c                IF (ninter.GE.1) EXIT

                srf_mask(isrf) = .TRUE.

              ENDDO  ! surface
c              IF (ninter.GE.1) EXIT
            ENDDO  ! zone
c            IF (ninter.GE.1) EXIT
          ENDDO  ! zone
c          IF (ninter.GE.1) EXIT
        ENDDO  ! zone
c        IF (ninter.GE.1) EXIT

        IF (n_list.GT.0) THEN
c...      Add the intersection point and register the intersection
c         and surface in a list:
          DO i5 = nv, iv, -1      
            v(1:3,i5+n_list) = v(1:3,i5) 
          ENDDO
          nv = nv + 1
          DO i5 = MAX(1,n_dat), 1, -1
            v_dat(1  ,i5+n_list) = v_dat(1  ,i5) + n_list  ! Shift the pointer since new vertices are being added
            v_dat(2:3,i5+n_list) = v_dat(2:3,i5) 
          ENDDO
          n_dat = n_dat + n_list
          
c          WRITE(0,*) 'register'
c          WRITE(0,*) 'd_list',d_list(1:n_list)
c          WRITE(0,*) 's_list',s_list(1:n_list)

          DO i6 = 1, n_list
            i5 = MAXLOC(d_list(1:n_list),1)           
c            WRITE(0,*) '  chosen one',i5,iv+i6-1,i6   
            v(1:3,iv+i6-1) = v_list(1:3,i5)
            v_dat(1,i6) = iv + i6 - 1        ! index of the new point in the V array
            v_dat(2,i6) = 3                  ! 3=intersection point (=1 origin point inside vessel, =2 origin outside)
            v_dat(3,i6) = s_list(i5)         ! surface index where the intersection took place 
            d_list(i5) = -1.0E+6
          ENDDO
        ENDIF
      ENDDO  ! field line trace line segments

c...  Store the data:
      trace_n = trace_n + 1
      trace_i(1,trace_n) = irib
      trace_i(2,trace_n) = ntrace + 1
      trace_i(3,trace_n) = ntrace + nv

      IF (ntrace+nv.GT.maxntrace) CALL AllocTrace(-1,MP_INCREASE_SIZE)

c      ntrace = 0
c      trace = 0.0D0

      DO i1 = 1, 3
        trace(i1,ntrace+1:ntrace+nv) = v(i1,1:nv)
      ENDDO

c...  Mark the intersection points and calculate the intersection 
c     angles:
      DO i1 = 1, n_dat
        i2 = v_dat(1,i1)
        trace(4,ntrace+i2) = DBLE(REAL(v_dat(2,i1)))

        isrf = v_dat(3,i1)
        CALL raySurfaceNormal(isrf,norm)

        vect(1:3) = trace(1:3,ntrace+i2+1) - trace(1:3,ntrace+i2-1)

        write(6,*) ' dprod',ntrace+i2
        write(6,*) '      ',norm
        write(6,*) '      ',vect
        write(6,*) '      ',trace(1:3,ntrace+i2+1)
        write(6,*) '      ',trace(1:3,ntrace+i2  )
        write(6,*) '      ',trace(1:3,ntrace+i2-1)

        length = DSQRT(vect(1)**2 + vect(2)**2 + vect(3)**2)
        vect = vect / length 

        CALL gmCalcDotProduct(norm,vect,dprod)

        angle = DABS(DACOS(dprod) * 180.0D0 / PI)     ! Get the angle between the normal and the trace
        IF (angle.GT.90.0D0) angle = 180.0D0 - angle  ! Want the smallest angle to the surface normal
        angle = 90.0D0 - angle                        ! Want the angle to the surface, not the normal

        trace(5,ntrace+i2) = angle

c        IF (ntrace+i2.GT.95200.AND.ntrace+i2.LT.9500) THEN

          write(6,*) '      ',length,dprod,angle
c        ENDIF
 

c        write(0,*) 'angle:',angle
      ENDDO

c...  Mark the origin:
      check_origin = .FALSE.
      DO i1 = 1, nv
        IF (trace(1,ntrace+i1).EQ.origin(1).AND.
     .      trace(2,ntrace+i1).EQ.origin(2).AND.
     .      trace(3,ntrace+i1).EQ.origin(3)) THEN

          trace_i(4,trace_n) = i1 + ntrace
          IF (inside) THEN
c            WRITE(0,*) 'FOUND ORIGIN! -- inside!'
            trace(4,ntrace+i1) = 1.0D0
          ELSE
c            WRITE(0,*) 'FOUND ORIGIN! -- outside!'
            trace(4,ntrace+i1) = 2.0D0            
          ENDIF
          check_origin = .TRUE.
          EXIT
        ENDIF
      ENDDO
      IF (.NOT.check_origin) 
     .  CALL ER('rayProcessTrace','Origin not found',*99)

      ntrace = ntrace + nv

c      WRITE(0,*) 'N_DAT:',n_dat
c      WRITE(0,*) 'V_DAT:',v_dat(1,1:n_dat)
c      WRITE(0,*) '     :',v_dat(2,1:n_dat)
c      WRITE(0,*) '     :',v_dat(3,1:n_dat)


c      DEALLOCATE(srf_check)
c      DEALLOCATE(srf_object)



      RETURN
 99   STOP
      END
c
c ====================================================================== 
c
