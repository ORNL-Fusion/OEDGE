c     -*-Fortran-*-
c
c ====================================================================
c
      SUBROUTINE osm_UpstreamProfile(type,index,val,coord,result)
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
c...    Find the outer midplane radius and the outer radial extent of the 
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
          xdata(i,IND_R) = xdata(i-1,IND_R) + step
c          WRITE(0,*) 'STEP:',i,step
        ENDDO
c       Recale so the total length covered is 2a:
        xdata(:,IND_R) = xdata(:,IND_R)*(2.0D0 * a_sep / xdata(NSTEP,1))

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
          CALL ER('osm_UpstreamProfile','Unknown quantity',*99)
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
            CALL ER('osm_UpstreamProfile','Unknown TYPE',*99) 
        ENDSELECT          

        IF (j.EQ.IND_TE.AND.ydata(0,1).EQ.0.0D0) 
     .    CALL ER('osm_UpstreamProfile','NE data must be assigned '//
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
              CALL ER('osm_UpstreamProfile','Unknown TYPE',*99) 
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
     .        CALL ER('osm_UpstreamProfile','A_END outside range',*99)    
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
     .        CALL ER('osm_UpstreamProfile','A_END outside range',*99)    
            i_end = i_end - 1
            x_end = xdata(i_end,IND_R)
            DO i = 1, NSTEP
              IF (xdata(i,IND_R).GT.a_etb) EXIT
            ENDDO
            IF (i.EQ.NSTEP+1)
     .        CALL ER('osm_UpstreamProfile','A_ETB outside range',*99)    
            cross(j) = i
          CASE DEFAULT
            CALL ER('osm_UpstreamProfile','Unknown TYPE',*99) 
        ENDSELECT

        IF (cross(j).EQ.-1)
     .    CALL ER('osm_UpstreamProfile','No core/SOL cross-over',*99)
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
          CALL ER('osm_UpstreamProfile','Unknown COORD value',*99)
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
     .  CALL ER('osm_UpstreamProfile','Independent coordinate '//
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
              CALL ER('osm_UpstreamProfile','Unknown quantity',*99)
          ENDSELECT

          WRITE(fp,10) 'val1:',i,a_sep,val1,frac,result
 10       FORMAT(A,I6,3F10.2,1P,3E10.2,0P)

          EXIT
        ENDIF
      ENDDO
      IF (i.EQ.NSTEP+1) 
     .  CALL ER('osm_UpstreamProfile','Plasma data not found',*99)


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

