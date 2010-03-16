c     -*-Fortran-*-
c
c ====================================================================
c
      SUBROUTINE LineCutTube(p1,p2,itube,p)
      USE mod_geometry
      USE mod_sol28_global
      IMPLICIT none

      INTEGER, INTENT(IN)  :: itube      
      REAL*8 , INTENT(IN)  :: p1(2),p2(2)
      REAL*8 , INTENT(OUT) :: p(2)

      INTEGER GetObject

      INTEGER ic1,ic2,icell,iobj
      REAL*8  p3(2),p4(2),t12,t34

      ic1 = tube(itube)%cell_index(1)
      ic2 = tube(itube)%cell_index(2)
      DO icell = ic1, ic2
        iobj = GetObject(icell,IND_CELL)                   
        CALL GetVertex(iobj,2,p3(1),p3(2))          
        CALL GetVertex(iobj,3,p4(1),p4(2))          
        CALL CalcInter(p1(1),p1(2),p2(1),p2(2),
     .                 p3(1),p3(2),p4(1),p4(2),t12,t34)
        IF (t12.GE.0.0D0.AND.t34.GE.0.0D0.AND.t34.LT.1.0D0) THEN
          p(1) = p3(1) + t34 * (p4(1) - p3(1))
          p(2) = p3(2) + t34 * (p4(2) - p3(2))
          EXIT
        ENDIF
      ENDDO
      IF (icell.EQ.ic2+1) 
     .  CALL ER('LineCutTube','Intersection not found',*99)

      RETURN
 99   STOP
      END
c
c ====================================================================
c
      SUBROUTINE SimplePlasmaProfile(type,index,val,coord,result)
      USE mod_geometry
      USE mod_sol28_global
      USE mod_interface
      IMPLICIT none
      
      INTEGER, INTENT(IN)  :: type,index,coord
      REAL   , INTENT(IN)  :: val
      REAL   , INTENT(OUT) :: result(3)
 
      INTEGER GetObject
      REAL*8  CalcPolygonToroidalVolume

      INTEGER, PARAMETER :: NSTEP  = 200, IND_R  = 1, IND_VOL = 2, 
     .                      IND_NE = 1  , IND_TE = 2, IND_TI  = 3
      REAL*8 , PARAMETER :: PI  = 3.14159265358979323846D0,
     .                      ECH = 1.602D-19 

      INTEGER fp,i,j,k,itube,icell,iobj,n,ic1,ic2,count,npro,cross(2),
     .        i_end
      LOGICAL cont
      REAL*8  r0,z0,p1(2),p2(2),p(2),
     .        a_sep,a_end,xdata(NSTEP,2),ydata(0:NSTEP,3),volume,step,
     .        total,x_end,param1,param2,a_ped1,a_ped2,
     .        r,mtanh,a_knee,a_slope,a_etb,a_delta,a_SOL,t_SOL,b_SOL,
     .        adjust,target,metric,diff,param,frac,pro1,val1,
     .        slope_mtanh,slope_exp,diff_min,xval,A,B,C,D
      REAL*8, ALLOCATABLE :: x(:),y(:),pro_a(:),pro_val(:)

      LOGICAL firstcall, storedata
      DATA    firstcall, storedata / .TRUE. , .TRUE./
      SAVE

      fp = 88

      result = -1.0


c...  Setup up the radial coordinate:

c     coord = 1 - linear on line segment
c           = 2 - linear on line segment, but from first to last ring intersection
c           = 3 - PSIn over range of applicability (like coord=2) 
c           = 4 - RHO
c           = 5 - PSIn (raw)
c           = 6 - linear on line segment, but from tube link to infinity

c...  Initialisation:      
      IF (firstcall) THEN
        xdata = 0.0D0
        ydata = 0.0D0

        r0 = grid%r0
        z0 = 0.0D0
c...    Find the outer midplane radius and the outer radial extent of the c
c       interpolation region:

        p1(1) = r0
        p1(2) = z0
        p2(1) = r0 + 100.0D0
        p2(2) = z0
        itube = grid%isep-1
        WRITE(0,*) 'ITUBE 1=',itube
        CALL LineCutTube(p1,p2,itube,p)
        a_sep = p(1) - p1(1)  

        itube = osmnode(index)%tube_range(2) 
        WRITE(0,*) 'ITUBE 2=',itube,index
        CALL LineCutTube(p1,p2,itube,p)
        a_end = p(1) - p1(1)
        
c...    Calculate core volume:
        itube = grid%isep-1
        ic1 = tube(itube)%cell_index(1)
        ic2 = tube(itube)%cell_index(2)
        n = ic2 - ic1 + 1
        ALLOCATE(x(n))
        ALLOCATE(y(n))
        n = 0
        DO icell = ic1, ic2
          n = n + 1
          iobj = GetObject(icell,IND_CELL)                   
          CALL GetVertex(iobj,2,x(n),y(n))          
        ENDDO
c Debug:
c        DEALLOCATE(x)
c        DEALLOCATE(y)
c        n = 4
c        ALLOCATE(x(n))
c        ALLOCATE(y(n))
c        x(1) = 0.0D0
c        y(1) = 0.0D0
c        x(2) = 0.0D0
c        y(2) = 1.0D0
c        x(3) = 1.0D0
c        y(3) = 1.0D0
c        x(4) = 1.0D0
c        y(4) = 0.0D0

        volume = CalcPolygonToroidalVolume(x,y,n)

        WRITE(fp,*) 'volume:',volume
        DO i = 1, n
          WRITE(fp,*) 'x,y:',i,x(i),y(i)
        ENDDO

        DEALLOCATE(x)
        DEALLOCATE(y)

        step = (2.0D0 * a_sep) / DBLE(NSTEP-1)

c        a_etb = DBLE(osmnode(index)%fit_p(3)) * a_sep

        WRITE(fp,*) 'STEP: ',a_sep,a_end,step

c       Distribute the independent variable so that the points are concentrated
c       near the separatrix:
        DO i = 2, NSTEP
          step = MAX(400.0,DBLE(ABS((i-1)-NSTEP/2))**2)
          xdata(i,1) = xdata(i-1,1) + step
c          WRITE(0,*) 'STEP:',i,step
        ENDDO
c       Recale so the total length covered is 2a:
        xdata(:,1) = xdata(:,1) * (2.0D0 * a_sep / xdata(NSTEP,1))

        total = 0.0D0
        DO i = 1, NSTEP/2
          xdata(i,IND_VOL) = 0.5D0*(xdata(i  ,IND_R)+xdata(i+1,IND_R)) * 
     .                             (xdata(i+1,IND_R)-xdata(i  ,IND_R))  
          total = total + xdata(i,IND_VOL)                       
        ENDDO
        xdata(:,IND_VOL) = xdata(:,IND_VOL) / total

        firstcall = .FALSE.
      ENDIF


      SELECTCASE (NINT(osmnode(index)%fit_quantity))
        CASE (1)
          j = IND_NE
        CASE (4)
          j = IND_TE
        CASE DEFAULT
          CALL ER('SamplePlasmaProfile','Unknown quantity',*99)
      ENDSELECT

      IF (ydata(0,j).EQ.0.0D0) THEN
        ydata(0,j) = -1.0     

        SELECTCASE (type)
          CASE (1)
            target  = DBLE(osmnode(index)%fit_p(1))
            a_knee  = DBLE(osmnode(index)%fit_p(2))
            a_etb   = DBLE(osmnode(index)%fit_p(3)) + a_sep
            a_delta = DBLE(osmnode(index)%fit_p(4))      
            a_SOL   = DBLE(osmnode(index)%fit_p(5))
            t_SOL   = DBLE(osmnode(index)%fit_p(6))
            b_SOL   = DBLE(osmnode(index)%fit_p(7))
            param1  = DBLE(osmnode(index)%fit_p(8))
            param2  = DBLE(osmnode(index)%fit_p(9))
            a_ped1 = a_etb - 2.0D0 * a_delta
            a_ped2 = a_etb + 2.0D0 * a_delta
          CASE (2)
            target  = DBLE(osmnode(index)%fit_p(1))
            a_etb   = DBLE(osmnode(index)%fit_p(2)) + a_sep
            a_SOL   = DBLE(osmnode(index)%fit_p(3))
            t_SOL   = DBLE(osmnode(index)%fit_p(4))
            b_SOL   = DBLE(osmnode(index)%fit_p(5))
            param1  = DBLE(osmnode(index)%fit_p(6))
            param2  = DBLE(osmnode(index)%fit_p(7))
          CASE DEFAULT
            CALL ER('SimplePlasmaProfile','Unknown TYPE',*99) 
        ENDSELECT          

        IF (j.EQ.IND_TE.AND.ydata(0,1).EQ.0.0D0) 
     .    CALL ER('SamplePlasmaProfile','NE data must be assigned '//
     .            'before T profile can be calculated',*99)

c...    Convert pedestal shape parameters into r/a coordinates:
c        SELECTCASE (coord)
c          (4) ! RHO : or (m), so need to update 

        a_slope = 0.0D0

        adjust = 0.0D0

        cont = .TRUE.
        DO WHILE(cont)
          cont = .FALSE.
c...      Calculate the core plasma profile:
          SELECTCASE (type)
            CASE (1)
              DO i = 1, NSTEP
                r         =(a_etb - xdata(i,IND_R)) / (2.0D0 * a_delta)
                mtanh     =((1.0D0 + a_slope * r) * EXP(r) - EXP(-r)) / 
     .                     (EXP(r) + EXP(-r))
                ydata(i,j)=(a_knee - a_SOL)/2.0D0 * (mtanh+1.D0) + a_SOL
              ENDDO
            CASE (2)
              DO i = 1, NSTEP 
                r = (a_etb - xdata(i,IND_R)) / a_etb
                ydata(i,j) = a_slope * a_SOL * r + a_SOL
c                WRITE(88,*) 'PRO:',i,ydata(i,j),diff
              ENDDO
            CASE DEFAULT
              CALL ER('SimplePlasmaProfile','Unknown TYPE',*99) 
          ENDSELECT          
c...      Calculate the metric for adjusting the slope of the core profile
c         to match the specified line averaged density or stored energy:
          SELECTCASE (j)
            CASE (IND_NE)  ! Line averaged density:   **** IS THIS RIGHT? ***
              metric = 0.0D0
              DO i = 1, NSTEP
                metric = metric + ydata(i,IND_NE) * xdata(i,IND_VOL)
              ENDDO
            CASE (IND_TE)  ! Stored energy 
              DO i = 1, NSTEP
                IF     (xdata(i,IND_R).LT.a_ped1) THEN
                  param = param1
                ELSEIF (xdata(i,IND_R).LT.a_ped2) THEN
                  frac = (xdata(i,IND_R) - a_ped1) / (a_ped2 - a_ped1)
                  frac = frac**3.0
                  param = (1.0D0 - frac) * param1 + frac * param2
                ELSE
                  param = param2
                ENDIF
                ydata(i,IND_TI) = ydata(i,IND_TE) * param
              ENDDO
              metric = 0.0D0
              DO i = 1, NSTEP
                metric = metric + 1.0D-06 * 1.5D0 *
     .                   xdata(i,IND_VOL) * volume * ydata(i,IND_NE) * 
     .                   (ydata(i,IND_TE) + ydata(i,IND_TI)) * ECH
              ENDDO
          ENDSELECT
c...      Evaluate the proximity of the metric to the requested value (TARGET) and
c         adjust the slope of the core profile accordingly:
          diff = (metric - target) / target
          IF (DABS(diff).GT.0.001) THEN
            cont = .TRUE.
            IF (adjust.EQ.0.0D0) THEN
              count = 0
              IF (diff.LT.0.0D0) THEN
                adjust = -1.0D0
              ELSE
                adjust =  1.0D0
              ENDIF
            ELSEIF (diff.GT.0.0D0.AND.adjust.LT.0.0D0.OR.
     .              diff.LT.0.0D0.AND.adjust.GT.0.0D0) THEN
              count = count + 1
              IF (count.EQ.10) THEN
                count = 0
                adjust = adjust * 10.0D0
              ENDIF
            ELSE
              count = 0
              adjust = -0.3D0 * adjust
            ENDIF
            a_slope = a_slope + adjust
          ENDIF
        ENDDO

c...    Add the exponential tail that extends into the SOL:
        SELECTCASE (type)
          CASE (1)
            DO i_end = 1, NSTEP
              IF (xdata(i_end,IND_R).GT.a_end) EXIT
            ENDDO
            IF (i_end.EQ.NSTEP+1)
     .        CALL ER('SamplePlasmaProfile','A_END outside range',*99)    
            i_end = i_end - 1
            x_end = xdata(i_end,IND_R)
            
            diff_min = 1.0D+20
            cross(j) = -1
            DO i = 1, NSTEP-1 
              IF (xdata(i,IND_R).LT.a_etb+2.0D0*a_delta.OR.
     .            xdata(i,IND_R).GT.a_etb+8.0D0*a_delta) CYCLE
            
              slope_mtanh = (ydata(i+1,j) - ydata(i,j)) / 
     .                      (xdata(i+1,1) - xdata(i,1) + 1.0D-10)
            
c             For y(x) = A e**(-1/t) + B ; C = y(0), D = y(1)
              xval = x_end - xdata(i,IND_R) 
            
              C = ydata(i,j)
              D = b_SOL
              A = (D - C) / (DEXP(-xval / t_SOL) - 1.0D0)
              B = C - A
              slope_exp  = -1.0D0 * A / t_SOL
            
              diff = DABS((slope_mtanh - slope_exp) / slope_exp)
              IF (diff.LT.diff_min) THEN
                diff_min = diff
                cross(j) = i + 1
              ENDIF
            ENDDO
          CASE (2)
            DO i_end = 1, NSTEP
              IF (xdata(i_end,IND_R).GT.a_end) EXIT
            ENDDO
            IF (i_end.EQ.NSTEP+1)
     .        CALL ER('SamplePlasmaProfile','A_END outside range',*99)    
            i_end = i_end - 1
            x_end = xdata(i_end,IND_R)
            DO i = 1, NSTEP
              IF (xdata(i,IND_R).GT.a_etb) EXIT
            ENDDO
            IF (i.EQ.NSTEP+1)
     .        CALL ER('SamplePlasmaProfile','A_ETB outside range',*99)    
            cross(j) = i
          CASE DEFAULT
            CALL ER('SimplePlasmaProfile','Unknown TYPE',*99) 
        ENDSELECT

        IF (cross(j).EQ.-1)
     .    CALL ER('SamplePlasmaProfile','No core/SOL cross-over',*99)
c...    Calculate exponential fall off into the SOL:
        ydata(cross(j):NSTEP,j) = 0.0D0

        C    = ydata(cross(j)-1,j)
        D    = b_SOL
        xval = x_end - xdata(cross(j)-1,IND_R)
        A    = (D - C) / (DEXP(-xval / t_SOL)- 1.0D0)
        B    = C - A

        DO i = cross(j), i_end
          xval = xdata(i,IND_R) - xdata(cross(j)-1,IND_R)
          ydata(i,j) = A * DEXP(-xval / t_SOL) + B
          WRITE(88,*) 'XVAL:',i,xval,ydata(i,j)
        ENDDO

c...    Set YDATA beyone A_END to a constant:
        ydata(i_end+1:NSTEP,j) = ydata(i_end,j)

c...    Set Ti from Te:
        IF (j.EQ.IND_TE) THEN 
          DO i = 1, NSTEP
            IF     (xdata(i,IND_R).LT.a_ped1) THEN
              param = param1
            ELSEIF (xdata(i,IND_R).LT.a_ped2) THEN
              frac = (xdata(i,IND_R) - a_ped1) / (a_ped2 - a_ped1)
              frac = frac**3.0
              param = (1.0D0 - frac) * param1 + frac * param2
            ELSE
              param = param2
            ENDIF
            ydata(i,IND_TI) = ydata(i,IND_TE) * param
c            WRITE(0,*) 'PARAM B:',i,param
          ENDDO
        ENDIF




        WRITE(fp,*) 'A     :',a_sep
        WRITE(fp,*) 'X_END :',x_end
        WRITE(fp,*) 'SUM   :',SUM(xdata(:,j))
        WRITE(fp,*) 'CROSS :',cross(j),a_SOL
        WRITE(fp,*) 'METRIC:',j,metric,target
        DO i = 1, NSTEP
          WRITE(fp,'(A,2I6,3F12.4,1P,E12.4,0P,2F10.1)') 
     .      'DATA:',i,j,xdata(i,1),xdata(i,1)/a_sep,xdata(i,2),
     .      ydata(i,1:3)
        ENDDO
        WRITE(fp,*) 
        DO i = 1, NSTEP-1
          WRITE(fp,'(A,2I6,F12.4,1P,3E12.4,0P)') 
     .      'SLOPE:',i,j,xdata(i,1),
     .       (ydata(i+1,1)-ydata(i,1))/(xdata(i+1,1)-xdata(i,1)),
     .       -a_SOL/t_SOL*EXP(-(xdata(i,1)-a_sep)/t_SOL),
     .       -(ydata(i,1)-b_SOL)/t_SOL
        ENDDO

      ENDIF


c...  Map the given coordinate to distance in meters along the minor radius:
      ALLOCATE(pro_a  (ntube))
      ALLOCATE(pro_val(ntube))
      SELECTCASE (coord)
        CASE (4) ! rho
          npro = 0
          DO itube = 1, grid%n
            IF (tube(itube)%rho.EQ.0.0) CYCLE
            npro = npro + 1
            pro_a  (npro) = DBLE(tube(itube)%rho) + a_sep
            pro_val(npro) = DBLE(tube(itube)%rho)
          ENDDO
        CASE (5) ! PSIn
          npro = 0
          DO itube = 1, grid%n
            IF (tube(itube)%rho.EQ.0.0) CYCLE
            npro = npro + 1
            pro_a  (npro) = DBLE(tube(itube)%rho ) + a_sep
            pro_val(npro) = DBLE(tube(itube)%psin)
c            WRITE(0,*) 'pro:',npro,pro_a(npro),pro_val(npro)
          ENDDO
        CASE DEFAULT
          CALL ER('SamplePlasmaProfile','Unknown COORD value',*99)
      ENDSELECT

      DO i = 1, npro-1
c        WRITE(0,*) '?:',val,pro_val(i  )*0.999D0,
c     .                      pro_val(i+1)*1.001D0
        IF (val.GT.pro_val(i  )-DABS(pro_val(i  ))*0.001D0.AND.
     .      val.LT.pro_val(i+1)+DABS(pro_val(i+1))*0.001D0) THEN
          frac = (val - pro_val(i)) / (pro_val(i+1) - pro_val(i))
          val1 = (1.0D0 - frac) * pro_a(i) + frac * pro_a(i+1)
          EXIT
        ENDIF
      ENDDO
      IF (i.EQ.npro) 
     .  CALL ER('SamplePlasmaProfile','Independent coordinate '//
     .          'not adjusted',*99)

      DEALLOCATE(pro_a  )
      DEALLOCATE(pro_val)

c...  Sample the appropriate profile:
      DO i = 1, NSTEP
        IF (val1.GE.xdata(i,IND_R).AND.val1.LE.xdata(i+1,IND_R)) THEN
          frac = (val1             - xdata(i,IND_R)) / 
     .           (xdata(i+1,IND_R) - xdata(i,IND_R))
          SELECTCASE (NINT(osmnode(index)%fit_quantity))
            CASE (1) ! ne
              result(1)=SNGL((1.0D0-frac)*ydata(i,1)+frac*ydata(i+1,1))
            CASE (4) ! Te and Ti
              result(2)=SNGL((1.0D0-frac)*ydata(i,2)+frac*ydata(i+1,2))
              result(3)=SNGL((1.0D0-frac)*ydata(i,3)+frac*ydata(i+1,3))
            CASEDEFAULT
              CALL ER('SamplePlasmaProfile','Unknown quantity',*99)
          ENDSELECT

          WRITE(fp,10) 'val1:',i,a_sep,val1,frac,result
 10       FORMAT(A,I6,3F10.2,1P,3E10.2,0P)

          EXIT
        ENDIF
      ENDDO
      IF (i.EQ.NSTEP+1) 
     .  CALL ER('SamplePlasmaProfile','Plasma data not found',*99)


      IF (ydata(1,3).NE.0.0D0.AND.storedata) THEN 
        storedata = .FALSE.
        CALL inOpenInterface('osm.idl.pedestal')
        CALL inPutData(a_sep                  ,'PED_A'       ,'m')   
        CALL inPutData(a_end                  ,'PED_A_END'   ,'m')   
        CALL inPutData(xdata(cross(1),IND_R  ),'PED_CROSS_NE','m')   
        CALL inPutData(xdata(cross(2),IND_R  ),'PED_CROSS_TE','m')   
        CALL inPutData(xdata(cross(2),IND_R  ),'PED_CROSS_TI','m')   
        CALL inPutData(xdata(1:NSTEP ,IND_R  ),'PED_R'       ,'m')   
        CALL inPutData(xdata(1:NSTEP ,IND_VOL),'PED_VOLFR'   ,'none')   
        CALL inPutData(ydata(1:NSTEP ,IND_NE ),'PED_NE'      ,'m-3')   
        CALL inPutData(ydata(1:NSTEP ,IND_TE ),'PED_TE'      ,'eV')   
        CALL inPutData(ydata(1:NSTEP ,IND_TI ),'PED_TI'      ,'eV')   
        CALL inCloseInterface
      ENDIF

      RETURN
 99   WRITE(0,*) '  VAL = ',val
      IF (ALLOCATED(pro_val)) THEN
        WRITE(0,*) '  PRO_VAL(1   ) = ',pro_val(1)
        WRITE(0,*) '  PRO_VAL(npro) = ',pro_val(npro)
      ENDIF
      STOP
      END



c
c ====================================================================
c
      SUBROUTINE RoutineInDevelopment
      IMPLICIT none
      RETURN
 99   STOP
      END
c
c ====================================================================
c
      SUBROUTINE LocalGridRefinement
      USE mod_geometry
      IMPLICIT none

      TYPE(type_srf   ) newsrf
      TYPE(type_object) newobj




      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE CalculateTeProfile(inode1,inode2,s,target)
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none
 
      INTEGER, INTENT(IN) :: inode1,inode2,target
      REAL*8 , INTENT(IN) :: s(0:icmax+1)

      INTEGER ion,ic,ic1,ic2,i1,i2,count
      INTEGER icstep
      LOGICAL cont
      REAL*8  te1,te2,Psol,L,k,stotal,vsign,adjust

      ion = 1

      ic1 = node(inode1)%icell
      ic2 = node(inode2)%icell
      icstep = SIGN(1,ic1-ic2)
c      te1 = DBLE(node(inode1)%te)
      IF (opt%bc(target).EQ.3) THEN
        IF (target.EQ.LO) te1 = te(ictarg(LO))
        IF (target.EQ.HI) te1 = te(ictarg(HI))
      ELSE
        te1 = DBLE(node(inode1)%te)
      ENDIF
      te2 = DBLE(node(inode2)%te)
c      IF (opt%bc(target).EQ.3) THEN
c        IF (target.EQ.LO) te2 = te(ictarg(LO))
c        IF (target.EQ.HI) te2 = te(ictarg(HI))
c      ELSE
c        te2 = DBLE(node(inode2)%te)
c      ENDIF

      k = DBLE(opt%te_kappa(target))
      L = 0.5D0 * DBLE(tube%smax)

c...  Set total power into flux tube:
      SELECTCASE (opt%te_ano_psol(target)) 
        CASE (1000)           ! ie ENEANO array from ref_plasma...? 
          Psol = 0.0D0  ! ???
c          STOP 'OKAY, READY'
        CASE (0)
c...      Conduction model estimate, all in upstream (initial guess only):
          Psol = 2.0D0/7.0D0 * (te2**3.5D0 - te1**3.5D0) * k / L   ! Stangeby, p 190
c          WRITE(0,*) 'PSOL:',psol,ic1,ic2
        CASE (1)
c...      From reference solution:
          Psol = ref_tube%eneano(target)
        CASE DEFAULT
          STOP 'USER OPTION NOT READY, YOU NINNY'
      ENDSELECT


      adjust = 0.0D0
      count = 0

              WRITE(logfp,'(A,I6)') 
     .          'TARGET:',target

      cont = .TRUE.
      DO WHILE (cont)
        cont = .FALSE.

        count = count + 1

        CALL EvolveTeProfile(inode1,inode2,s,target,Psol,te1,te2)


c...    Analyse electron profiles and modify parameters, if necessary:

        IF (DABS(te(ic1)-te1).GT.MIN(0.001D0,0.05D0*te1)) THEN
c          vsign = SIGN(1.0D0,te(ic1)-te1) 
c          IF (target.EQ.LO) vsign = vsign * SIGN(1.0D0,te2-te1)
          vsign = SIGN(1.0D0,te(ic1)-te1) * SIGN(1.0D0,te2-te1)
          IF (adjust.EQ.0.0D0) THEN
            adjust = 0.3D0 * vsign
          ELSEIF ((adjust.GT.0.0D0.AND.vsign.EQ.-1.0D0).OR.
     .            (adjust.LT.0.0D0.AND.vsign.EQ. 1.0D0)) THEN
            adjust = -0.3D0 * adjust
          ENDIF
          cont = .TRUE.
        ENDIF

c...    Analyse electron profiles and modify parameters, if necessary:

        SELECTCASE (opt%bc(target))  ! (opt%te_ano_psol(target))  ! *** CHANGE VARIABLE *** 
          CASE (1)
c...        Adjust the anomalous power term on the half ring:
c...        Scale Psol to improve match to specified target temperature:
            SELECTCASE (1)
              CASE (1)  ! Power into SOL:
                Psol = Psol * (1.0D0 + adjust)
              CASE DEFAULT
                STOP 'NOT READY, SORRY'
            ENDSELECT
            IF (log.GE.2) THEN
              WRITE(logfp,'(A,I3,3F11.4,1P,2E14.6,0P)') 
     .          'ADJUST:',count,te(ic1),te1,vsign,adjust,Psol
            ENDIF

          CASE (3)

            node(inode2)%te = node(inode2)%te - adjust
            te2 = DBLE(node(inode2)%te)

            IF (log.GE.2) THEN
              WRITE(logfp,'(A,3F11.4,1P,E14.6,0P,F11.4)') 
     .          'ADJUST:',te(ic1),te1,vsign,adjust,te2
            ENDIF

          CASE DEFAULT
            CALL ER('EvolveTeProfile','Bad BC',*99)
        ENDSELECT

c...    Analyse ion profiles:

c...    If convection being used, need to run for extra iterations
c       to make sure the solution is converged:

        IF (count.EQ.50) THEN
          te(ic1:ic2) = te1
          WRITE(0,*) '*** ENERGY MODEL FAILED ***'
          WRITE(logfp,*) '*** ENERGY MODEL FAILED ***'
          EXIT
c          STOP 'EXCESS ITERATIONS'
        ENDIF

      ENDDO





c...  Store Psol information:
      tube%eneano(target) = Psol   ! Keep this up...?

c...  Linear distribution of anomalous power for now...
      SELECTCASE (opt%te_ano(target))
        CASE(1000)
        CASE(1)
        CASE(2)
c...      Everything in at the midplane:
          eneano(ic2+icstep) = Psol / sdelta(ic2+icstep)  ! Approximate..!
c          WRITE(0,*) 'ENEA:',psol,ic2+icstep
        CASE(3)
c...      Distributed uniformly along the half ring:
          IF (ic1.GT.ic2) THEN
            ic  = ic1
            ic1 = ic2
            ic2 = ic
          ENDIF
          ic1 = MAX(ic1,1)
          ic2 = MIN(ic2,icmax)
          IF (target.EQ.LO) ic2 = ic2 - 1
          IF (target.EQ.HI) ic1 = ic1 + 1
          stotal = SUM(sdelta(ic1:ic2))
          eneano(ic1:ic2) = Psol / stotal
c          WRITE(0,*) 'ENEA:',psol,ic1,ic2,icstep,stotal
        CASE DEFAULT
          CALL ER('EvolveTeProfile','Bad TE_ANO',*99)
      ENDSELECT


      IF (log.GE.2) THEN
        WRITE(logfp,'(A)') 'DONE TE ITERATIONS'
      ENDIF

      IF (target.EQ.HI.AND.opt%bc(target).EQ.1) THEN
        WRITE(logfp,*) 'INTEGRATING ENERGY SOURCES AGAIN:'
        CALL IntegrateSources(3)

        ic1 = ictarg(LO)
        ic2 = ictarg(HI)

c        qe(ic1) = -eneint(icmid  ,1)
c        qe(ic2) = -(eneint(TOTAL,1)-eneint(icmid,1))

        WRITE(logfp,*) 'QE   :',qe(0),qe(icmax+1)
        WRITE(logfp,*) 'ISAT :',isat(ictarg(LO:HI),ion)
        WRITE(logfp,*) 'TE   :',te(ictarg(LO:HI))

        gamma(LO) = qe(0)       / (isat(ic1,ion)*ECH) / te(ic1)
        gamma(HI) = qe(icmax+1) / (isat(ic2,ion)*ECH) / te(ic2)

c        gamma = 3.0D0
c        qe(ic1) = 3.0D0 * (ECH*te(ic1)) * isat(ic1,ion)
c        qe(ic2) = 3.0D0 * (ECH*te(ic2)) * isat(ic1,ion)

        WRITE(logfp,*) 'GAMMA:',gamma(LO:HI)
        WRITE(logfp,*) '...and again...'

        CALL IntegrateSources(3)
c        gamma(LO:HI) = 3.0D0
c        STOP 'sdgsdg'
      ENDIF


      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE EvolveTeProfile(inode1,inode2,s,target,Psol,te1,te2)
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none
 
      INTEGER, INTENT(IN) :: inode1,inode2,target
      REAL*8 , INTENT(IN) :: s(0:icmax+1),Psol,te1,te2

      INTEGER :: ic,ic1,ic2,icstep,count,ion
      REAL*8  :: tgrad,x,k,frac,tarsign,
     .           flux,qano(0:S28_MAXNKS+1),
     .           qe_src(0:S28_MAXNKS+1),
     .           qconv_e(0:S28_MAXNKS+1),qconv_i(0:S28_MAXNKS+1)
      LOGICAL :: cont

c      te1 = DBLE(node(inode1)%te)
c      te2 = DBLE(node(inode2)%te)
      ic1 = node(inode1)%icell
      ic2 = node(inode2)%icell
      icstep = SIGN(1,ic1-ic2)
      tarsign = DBLE(icstep)

c      IF (opt%bc(target).GE.2.AND.
c     .    (ic1.EQ.ictarg(target).AND.te(ic1).EQ.0.0D0)) THEN
c        WRITE(0,*) 'CYCLING -------------------------'
c        RETURN
c      ENDIF


c      WRITE(logfp,*) 'TE: ',ic1,ic2,icstep,target
c      WRITE(logfp,*) '    ',te1,te2

      ion = 1


c...  Integrate e/i cf/vol sources/sinks:
      qe_src = 0.0D0
      IF (icstep.EQ.-1) THEN
        qe_src(ic1:ic2) = eneint(ic2,1) - eneint(ic1:ic2,1)  ! Should I set the INT arrays from 0:icmax+1, and have 
        IF (ic1.EQ.0) qe_src(ic1) = eneint(ic2,1)            ! TOTAL = icmax+1 to make things like this simpler?
c        WRITE(0,*) 'QE:',qe(0:icmax+1)
c        WRITE(0,*) 'QE:',enesrc(0:icmax+1,1)
      ELSE
        qe_src(ic2:ic1) = eneint(ic2:ic1,1) - eneint(ic2,1)  
        IF (ic1.EQ.icmax+1) qe_src(ic1) = eneint(TOTAL,1)-eneint(ic2,1)            
c        WRITE(0,*) 'QE:',qe(0:icmax+1)
c        WRITE(0,*) 'QE:',eneint(ic2,1),ic1,icmax+1
c        WRITE(0,*) 'QE:',enesrc(0:icmax+1,1)
      ENDIF

c...  Initializations:
      k = DBLE(opt%te_kappa(target))
      x = 5.0D0 / 2.0D0
      qano  = 0.0D0
      te(ic2) = te2
 
c...  Setup anomalous energy input:    NEED TO ADD EFFICIENCY, ie DON'T ALWAYS HAVE TO CALCULATE EVERYTHING...
      SELECTCASE (opt%te_ano(target))
        CASE(1000)
        CASE(0)
c...      None:
        CASE(1)
c...      From reference solution - ANO source incorporated into QE_SRC via ENEINT:
        CASE(2)
c...      All power in at the symmetry point:
          qano = Psol 
        CASE(3)
c...      Power distributed uniformly about the symmetry point(!):
          DO ic = ic2+icstep, ic1, icstep  ! No power assigned to the symmetry point cell...
            frac = (s(ic ) - (s(ic2) - 0.5D0 * sdelta(ic2))) / 
     .             (s(ic1) - (s(ic2) - 0.5D0 * sdelta(ic2)))
            qano(ic) = Psol * frac                      
          ENDDO
c        CASE(3)
c...      Between x-points:
c        CASE(4)
c...      B-field depenence... use CalcProfile?   What if heat is not centred at symmetry point? 
        CASE DEFAULT
          STOP 'USER TE OPTION NOT READY'
      ENDSELECT

c...  Electron convected heat flux:
      SELECTCASE (opt%te_conv(target))
        CASE(0)
c...      None:
          DO ic = ic2, ic1, icstep
            qconv_e(ic) = 0.0D0
          ENDDO
        CASE(1)
          te(ic1) = te1
          DO ic = ic2, ic1, icstep
c            WRITE(logfp,*) 'CONV:',ic,te(ic),ne(ic),vi(ic,ion)
            IF (te(ic).NE.0.0D0) THEN
              flux = ne(ic) * vi(ic,ion) * tarsign
              qconv_e(ic) = 5.0D0 / 2.0D0 * ECH * te(ic) * flux
            ENDIF
          ENDDO
        CASE DEFAULT
          STOP 'USER TE OPTION NOT READY'
      ENDSELECT

c...  Ion convected heat flux:
      SELECTCASE (0)  ! (opt%ti_conv(target))  ! Don't need this when solving the electron channel...
        CASE(0)
c...      None:
          DO ic = ic2, ic1, icstep
            qconv_i(ic) = 0.0D0
          ENDDO
        CASE(1)
          te(ic1) = te1                        ! *** ASSUMES Ti = Te, on this line and 4 below ***
          DO ic = ic2, ic1, icstep
            IF (te(ic).NE.0.0D0) THEN
              flux = ne(ic) * vi(ic,ion) * tarsign
              qconv_i(ic) = (5.0D0 / 2.0D0 * ECH * te(ic) +  
     .                       0.5D0 * mi(ion) * vi(ic,ion)**2) * flux
            ENDIF
          ENDDO
        CASE DEFAULT
          STOP 'USER TE OPTION NOT READY'
      ENDSELECT

      qconv = qconv_e + qconv_i

c...  Electron-ion equilibration:
      IF (.FALSE.) THEN
      ENDIF

c...  Calculate Te profile:
      te(ic1) = 0.0D0
      DO ic = ic2+icstep, ic1, icstep
        qcond(ic) = qano(ic) + qe_src(ic) - qconv(ic) 
        tgrad = -qcond(ic) / (k * te(ic-icstep)**x)                ! Reference...
c...    Note: This is not strictly correct since the convection is a cell centered quantity
c       but the temperature evolution occurs between cell centres, i.e. the convecion
c       term should be a wieghted average, or some such, of the convection
c       across the 2 cells -- similar inconcistencies can likely be found elsewhere
c       and will likely need to be resolved before a full cross-field transport model
c       can work... might even run into trouble before then.  For a similar slop, the
c       treatment of the last half-cell on each ring also needs careful review.
        te(ic) = te(ic-icstep) + (s(ic-icstep) - s(ic)) * tgrad

        qcond(ic) = qcond(ic) !* tarsign
        qconv(ic) = qconv(ic) !* tarsign

c        qe(ic) = qcond(ic) + qconv(ic)


        IF (log.GE.2) THEN
          IF (ic.EQ.ic2+icstep) THEN
            WRITE(logfp,'(A4,A8,2X,3A12,2X,2A12)')
     .        'IND','Te','qcond','qconv_e','qconv_i','qano','qe_src'   
          ENDIF
          WRITE(logfp,'(I4,F8.2,2X,1P,3D12.4,2X,2D12.4,1P,
     .                  2E10.2)')
     .      ic,
     .      te(ic),
     .      qcond(ic),qconv_e(ic),qconv_i(ic),
     .      qano(ic),qe_src(ic),
     .      tgrad,(s(ic-icstep) - s(ic))
        ENDIF
c          IF (target.EQ.LO)
c     .      WRITE(0,'(A,2I6,2F10.2,1P,6E10.2,2X,2E10.2,0P)') 
c     .        '  Te->',ic,ic-icstep,s(ic),te(ic),vi(ic,ion),flux,
c     .        qano(ic),qe(ic),qconv(ic),qcond(ic),vi(ic,ion),ne(ic)
c        IF (te(ic).LE.0.1D0*te1) EXIT  
c        IF (te(ic).LE.0.5D0*te1.OR.te(ic).GT.1.5D0*te2) THEN
c          te(ic) = MIN(MAX(te(ic),te1),te2)
c          EXIT  ! There is no way to avoid the near target dip in 
c        ENDIF
        IF (te(ic).LE.0.5D0*MIN(te1,te2)) EXIT  ! There is no way to avoid the near target dip in 
c        IF (te(ic).LE.0.5D0*te1) EXIT  ! There is no way to avoid the near target dip in 
c        IF (te(ic).LT.te1) EXIT        ! Te if Qe is ill-posed?
      ENDDO


c      WRITE(0,*) 'QCOND:',qcond(ic1),qcond(ic2)
c      WRITE(0,*) 'TE   :',te(ic1),te(ic2)
c      WRITE(0,*) 'PSOL :',psol,ic1,ic2,opt%te_ano_psol(target)



      RETURN
 99   STOP
      END





