c     -*-Fortran-*-c
c ======================================================================
c
      SUBROUTINE MarkEdge(p1,p2,pdat,ndat,min_dist)
      IMPLICIT none

      INTEGER, INTENT(IN)    :: ndat
      REAL*8 , INTENT(IN)    :: p1(3),p2(3),min_dist
      REAL*8 , INTENT(INOUT) :: pdat(5,0:*)

      REAL*8  CalcPerp2

      INTEGER i
      REAL*8  t

      DO i = 1, ndat
        pdat(4,i) = CalcPerp2(p1,p2,pdat(1:3,i),t)
        IF (pdat(4,i).LT.min_dist) pdat(5,i) = 1.0D0
      ENDDO       

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c
      SUBROUTINE rayLoadITERFWP(ielement)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

c...  Input:
      INTEGER ielement

      INTEGER AddVertex,AddSurface
      REAL*8  CalcPerp2

      TYPE(type_surface) newsrf
      INTEGER   idum1,istart,fp,ndat,i,j,k,ny,nz,iy,iz,psum,
     .          check_n,dist_list(10),dist_n,tcnt,cnt,
     .          dtri(3,1000),dsid(3,1000),ntri,i1,i2,i3,i4,ishift,
     .          ibegin,
     .          stri(1000),match,list_n,list_i(0:1000),col_i(1000),
     .          col_max,
     .          col_n,col_min,col_cnt,icol,new_s,new_i(1000),new_n
      LOGICAL   cont,status
      CHARACTER dummy*1024,file*1024
      REAL*8    newvtx(3,25),mat(3,3),angle,frac ,tmpvtx(3,25),
     .          pdat(5,0:20000),rad,rdum,check,check_list(0:100),
     .          avg(3),p1(3),p2(3),p3(3),angle1,angle2,angle3,len1,
     .          len2,len3,focus(3),v1(3),v2(3),v3(3),t,min_dist,dist,
     .          step,v_panel(3),dotprod,v_norm(3),v_left(3),v_right(3),
     .          corner(3,4),d1,d2,d3

      newvtx = 0.0D0
      istart = nsrf + 1


      SELECTCASE (opt%obj_option(ielement))
c       ----------------------------------------------------------------
        CASE (1)
c...      ITER wall (no small feat):
          IF (.FALSE.) THEN
            newvtx(1,1) =  4.000D0
            newvtx(2,1) =  1.000D0
            newvtx(1,2) =  4.000D0
            newvtx(2,2) = -1.000D0
            DO i1 = 2, 2, -1
              newsrf%type = SP_LINE_SEGMENT
              newsrf%nvtx = 2
              newsrf%ivtx(1) = AddVertex(newvtx(1,i1  ))
              newsrf%ivtx(2) = AddVertex(newvtx(1,i1-1))
              idum1 = AddSurface(newsrf)
            ENDDO
            IF (nobj+1.GT.MAX3D) 
     .        CALL ER('LoadVesselStructures','Insufficient array '//
     .                'bounds for all objects',*99)    
            IF (istart.GT.nsrf) THEN
              WRITE(0,*) 'LoadVesselStructures: Strange, no objects'
              RETURN
            ENDIF
            nobj = nobj + 1
            WRITE(0,*) 'VESSEL STRUCTURE IOBJ:',nobj
            obj(nobj)%index       = ielement  ! nobj
            obj(nobj)%type        = OP_EMPTY
            obj(nobj)%mode        = 0      
            obj(nobj)%surface     = 1      ! SOLID
            obj(nobj)%wedge1      = 0
            obj(nobj)%wedge2      = 0
            obj(nobj)%colour      = 1
            obj(nobj)%orientation = 1      ! CW
            obj(nobj)%ik          = 0
            obj(nobj)%ir          = 0
            obj(nobj)%in          = -1  ! What should this be?
            obj(nobj)%ivolume     = 0
            obj(nobj)%nside       = 1
            obj(nobj)%iside(1,1)  = istart ! Start index of range of surfaces in surface array, from loading code above
            obj(nobj)%iside(1,2)  = nsrf   ! End index of range of surfaces in surface array
            obj(nobj)%gsur(1)     = GT_TC
            obj(nobj)%tsur(1)     = SP_VESSEL_WALL
            obj(nobj)%reflec(1)   = opt%obj_reflec(ielement)
c..         Defunct:
            obj(nobj)%nsur        = 0
            obj(nobj)%ipts(2,1)   = 0
            obj(nobj)%nmap(1)     = 0
          ENDIF



          fp = 99
          OPEN(fp,FILE='3d_wall.dat',FORM='FORMATTED',STATUS='OLD',
     .         ERR=98)       
          ndat = 0
          DO WHILE (.TRUE.)
            READ(fp,'(A)',END=10,ERR=98) dummy
            IF (dummy(1:1).EQ.'*'.OR.LEN_TRIM(dummy).LT.5) CYCLE
            WRITE(0,*) 'dummy>',TRIM(dummy)//'<'
            ndat = ndat + 1
            READ(dummy,*) pdat(1:2,ndat)
          ENDDO
 10       CONTINUE
          CLOSE(fp)
          WRITE(0,*) 'Ndat:',ndat

c          pdat(1:2,1) =  (/ 8.09730, -1.48350 /)
c          pdat(1:2,2) =  (/ 7.87680, -1.36560 /)
c          pdat(1:2,3) =  (/ 7.28420, -2.24640 /)
c          pdat(1:2,4) =  (/ 7.46680, -2.41720 /)

          DO i = 1, ndat, 4
            j = i - 1
            DO angle = 0.0D0, 359.0D0, 20.0D0
              newvtx(1:3,1) = (/ pdat(1,1+j), pdat(2,1+j), -0.50D0 /)
              newvtx(1:3,2) = (/ pdat(1,2+j), pdat(2,2+j),  0.00D0 /)
              newvtx(1:3,3) = (/ pdat(1,3+j), pdat(2,3+j),  0.00D0 /)
              newvtx(1:3,4) = (/ pdat(1,4+j), pdat(2,4+j), -0.50D0 /)
c             Rotate about y-axis (swing):
              CALL RotateVertices(angle,newvtx(1,1),4)
              CALL AddPolygon(SP_PLANAR_POLYGON,newvtx(1,1),4)
              newvtx(1:3,1) = (/ pdat(1,2+j), pdat(2,2+j),  0.00D0 /)
              newvtx(1:3,2) = (/ pdat(1,1+j), pdat(2,1+j),  0.50D0 /)
              newvtx(1:3,3) = (/ pdat(1,4+j), pdat(2,4+j),  0.50D0 /)
              newvtx(1:3,4) = (/ pdat(1,3+j), pdat(2,3+j),  0.00D0 /)
c             Rotate about y-axis (swing):
              CALL RotateVertices(angle,newvtx(1,1),4)
              CALL AddPolygon(SP_PLANAR_POLYGON,newvtx(1,1),4)
c             Top cap:
              newvtx(1:3,1) = (/ pdat(1,1+j), pdat(2,1+j), +0.50D0 /)
              newvtx(1:3,2) = (/ pdat(1,2+j), pdat(2,2+j),  0.00D0 /)
              newvtx(1:3,3) = (/ pdat(1,1+j), pdat(2,1+j), -0.50D0 /)
c             Rotate about y-axis (swing):
              CALL RotateVertices(angle,newvtx(1,1),3)
              CALL AddPolygon(SP_PLANAR_POLYGON,newvtx(1,1),3)
c             Bottom:
              newvtx(1:3,1) = (/ pdat(1,4+j), pdat(2,4+j), -0.50D0 /)
              newvtx(1:3,2) = (/ pdat(1,3+j), pdat(2,3+j),  0.00D0 /)
              newvtx(1:3,3) = (/ pdat(1,4+j), pdat(2,4+j), +0.50D0 /)
c             Rotate about y-axis (swing):
              CALL RotateVertices(angle,newvtx(1,1),3)
              CALL AddPolygon(SP_PLANAR_POLYGON,newvtx(1,1),3)
            ENDDO
          ENDDO

          IF (nobj+1.GT.MAX3D) 
     .      CALL ER('LoadVesselStructures','Insufficient array '//
     .              'bounds for all objects',*99)    
          IF (istart.GT.nsrf) THEN
            WRITE(0,*) 'LoadVesselStructures: Strange, no objects'
            RETURN
          ENDIF

          nobj = nobj + 1
          WRITE(0,*) 'VESSEL STRUCTURE IOBJ:',nobj
          obj(nobj)%index       = ielement  ! nobj
          obj(nobj)%type        = OP_EMPTY
          obj(nobj)%mode        = 0      
          obj(nobj)%surface     = 1      ! SOLID
          obj(nobj)%wedge1      = 0
          obj(nobj)%wedge2      = 0
          obj(nobj)%colour      = 1
          obj(nobj)%orientation = 1      ! CW
          obj(nobj)%ik          = 0
          obj(nobj)%ir          = 0
          obj(nobj)%in          = -1  ! What should this be?
          obj(nobj)%ivolume     = 0
          obj(nobj)%nside       = 1
          obj(nobj)%iside(1,1)  = istart ! Start index of range of surfaces in surface array, from loading code above
          obj(nobj)%iside(1,2)  = nsrf   ! End index of range of surfaces in surface array
          obj(nobj)%gsur(1)     = GT_TD
          obj(nobj)%tsur(1)     = SP_VESSEL_WALL
          obj(nobj)%reflec(1)   = opt%obj_reflec(ielement)
c..       Defunct:
          obj(nobj)%nsur        = 0
          obj(nobj)%ipts(2,1)   = 0
          obj(nobj)%nmap(1)     = 0
c       ----------------------------------------------------------------
        CASE (2)
          pdat = -999.0D0

          fp = 99
          file = opt%obj_fname(ielement)
          OPEN(fp,FILE=TRIM(file),FORM='FORMATTED',STATUS='OLD',
     .         ERR=98)       
          ndat = 0
          DO WHILE (.TRUE.)
            READ(fp,'(A)',END=20,ERR=98) dummy
            IF (dummy(1:1).EQ.'*'.OR.LEN_TRIM(dummy).LT.5) CYCLE
            WRITE(0,*) 'dummy>',TRIM(dummy)//'<'
            ndat = ndat + 1
            READ(dummy,*) rdum,rdum,pdat(1:3,ndat)
          ENDDO
 20       CONTINUE
          CLOSE(fp)
          WRITE(0,*) 'Ndat:',ndat

c...      Swap Y and Z to match RAY coordinate system:
          pdat(5,:) = pdat(3,:)
          pdat(3,:) = pdat(2,:)
          pdat(2,:) = pdat(5,:)

          DO angle = 0.0D0, 359.0D0, 360.0D0 
            DO i = 1, ndat, 4
              newvtx(1:3,1) = pdat(1:3,i+0)
              newvtx(1:3,2) = pdat(1:3,i+1)
              newvtx(1:3,3) = pdat(1:3,i+2)
              newvtx(1:3,4) = pdat(1:3,i+3)
c             Rotate about y-axis (swing):
              CALL RotateVertices(angle,newvtx(1,1),4)
              CALL AddPolygon(SP_PLANAR_POLYGON,newvtx(1,1),4)
            ENDDO
          ENDDO

          IF (nobj+1.GT.MAX3D) 
     .      CALL ER('LoadVesselStructures','Insufficient array '//
     .              'bounds for all objects',*99)    
          IF (istart.GT.nsrf) THEN
            WRITE(0,*) 'LoadVesselStructures: Strange, no objects'
            RETURN
          ENDIF

          nobj = nobj + 1
          WRITE(0,*) 'VESSEL STRUCTURE IOBJ:',nobj,istart,nsrf
          obj(nobj)%index       = ielement  ! nobj
          obj(nobj)%type        = OP_EMPTY
          obj(nobj)%mode        = 0      
          obj(nobj)%surface     = 1      ! SOLID
          obj(nobj)%wedge1      = 0
          obj(nobj)%wedge2      = 0
          obj(nobj)%colour      = 1
          obj(nobj)%orientation = 1      ! CW
          obj(nobj)%ik          = 0
          obj(nobj)%ir          = 0
          obj(nobj)%in          = -1  ! What should this be?
          obj(nobj)%ivolume     = 0
          obj(nobj)%nside       = 1
          obj(nobj)%iside(1,1)  = istart ! Start index of range of surfaces in surface array, from loading code above
          obj(nobj)%iside(1,2)  = nsrf   ! End index of range of surfaces in surface array
          obj(nobj)%gsur(1)     = GT_TD
          obj(nobj)%tsur(1)     = SP_VESSEL_WALL
          obj(nobj)%reflec(1)   = opt%obj_reflec(ielement)
c..       Defunct:
          obj(nobj)%nsur        = 0
          obj(nobj)%ipts(2,1)   = 0
          obj(nobj)%nmap(1)     = 0

c       ----------------------------------------------------------------
        CASE (3)
          pdat = -999.0D0

          fp = 99
          file = opt%obj_fname(ielement)
          OPEN(fp,FILE=TRIM(file),FORM='FORMATTED',STATUS='OLD',
     .         ERR=98)       
          ndat = 0
          DO WHILE (.TRUE.)
            READ(fp,'(A)',END=30,ERR=98) dummy
            IF (dummy(1:1).EQ.'*'.OR.LEN_TRIM(dummy).LT.5) CYCLE
            WRITE(0,*) 'dummy>',TRIM(dummy)//'<'
            ndat = ndat + 1
            READ(dummy,*) rdum,pdat(1:3,ndat)
          ENDDO
 30       CONTINUE
          CLOSE(fp)
          WRITE(0,*) 'Ndat:',ndat

c...      Sort the data:
          WRITE(0,*) 'SORTING 1...'
          cont = .TRUE.
          DO WHILE(cont)
            cont = .FALSE.
            DO i = 1, ndat-1
              IF (pdat(2,i).GT.pdat(2,i+1)) THEN 
                pdat(:,0  ) = pdat(:,i  )
                pdat(:,i  ) = pdat(:,i+1)
                pdat(:,i+1) = pdat(:,0  )
                cont = .TRUE. 
              ENDIF
            ENDDO
          ENDDO
          WRITE(0,*) 'SORTING 2...'
          cont = .TRUE.
          DO WHILE(cont)
            cont = .FALSE.
            DO i = 1, ndat-1
              IF (pdat(3,i).GT.pdat(3,i+1).AND.
     .            pdat(2,i).EQ.pdat(2,i+1)) THEN 
                pdat(:,0  ) = pdat(:,i  )
                pdat(:,i  ) = pdat(:,i+1)
                pdat(:,i+1) = pdat(:,0  )
                cont = .TRUE.
              ENDIF
            ENDDO
          ENDDO

c         Define the back-plane of the panel:          

c         Find average point for the tile:
          avg(1) = SUM(pdat(1,1:ndat)) / DBLE(ndat)
          avg(2) = SUM(pdat(2,1:ndat)) / DBLE(ndat)
          avg(3) = SUM(pdat(3,1:ndat)) / DBLE(ndat)
          pdat(4,1:ndat) = DSQRT( (avg(1) - pdat(1,1:ndat))**2 + 
     .                            (avg(2) - pdat(2,1:ndat))**2 + 
     .                            (avg(3) - pdat(3,1:ndat))**2 )

          pdat(5,:) = 0.0D0  ! Used to mark boundary points

c          DO i = 1, 100
c            WRITE(0,'(A,I6,3F10.3,F12.6)') 'pdat:',i,pdat(1:4,i)
c          ENDDO

          list_n = 0
          DO i = 1, 4
            j = MAXLOC(pdat(4,1:ndat),1) 
            pdat(4,j) = -1.0D6
            pdat(5,j) =  1.0D0
            list_n         = list_n + 1
            list_i(list_n) = j
          ENDDO       
c         Sort:          
          DO i = 1, list_n-1
            DO j = i+1, list_n
              IF (list_i(i).GT.list_i(j)) THEN
                list_i(0) = list_i(i)
                list_i(i) = list_i(j)
                list_i(j) = list_i(0)
              ENDIF
            ENDDO
          ENDDO

c         Divide all the points into column groups:
          check = -999.0 
          col_n = 0
          DO i = 1, ndat
            IF (DABS(pdat(2,i)-check).GT.1.0D-6) THEN 
              col_n = col_n + 1
              col_i(col_n) = i
              check = pdat(2,i)
            ENDIF
          ENDDO
          col_n = col_n + 1
          col_i(col_n) = ndat + 1


          min_dist = 0.001D0

          p1 = pdat(1:3,list_i(1))
          p2 = pdat(1:3,list_i(2))
          CALL MarkEdge(p1,p2,pdat,ndat,min_dist)
          v1 = pdat(1:3,list_i(2)) - pdat(1:3,list_i(1))
          v2 = pdat(1:3,list_i(3)) - pdat(1:3,list_i(1))
          CALL gmCalcCrossProduct(v1,v2,v3)
          p1 = pdat(1:3,list_i(1))
          p2 = p1 + v3 
          CALL MarkEdge(p1,p2,pdat,ndat,min_dist)
          p1 = pdat(1:3,list_i(2))
          p2 = p1 + v3 
          CALL MarkEdge(p1,p2,pdat,ndat,min_dist)

          p1 = pdat(1:3,list_i(3))
          p2 = pdat(1:3,list_i(4))
          CALL MarkEdge(p1,p2,pdat,ndat,min_dist)
          v1 = pdat(1:3,list_i(4)) - pdat(1:3,list_i(3))
          v2 = pdat(1:3,list_i(1)) - pdat(1:3,list_i(3))
          CALL gmCalcCrossProduct(v1,v2,v3)
          p1 = pdat(1:3,list_i(3))
          p2 = p1 + v3 
          CALL MarkEdge(p1,p2,pdat,ndat,min_dist)
          p1 = pdat(1:3,list_i(4))
          p2 = p1 + v3 
          CALL MarkEdge(p1,p2,pdat,ndat,min_dist)

c         Also mark the first and last points of the columns with fewer
c         points:
          DO i = 1, col_n-1
            p1 = pdat(1:3,list_i(1))
            p2 = pdat(1:3,list_i(3))
            DO j = col_i(i), col_i(i+1)-1
              pdat(4,j) = CalcPerp2(p1,p2,pdat(1:3,j),t)
            ENDDO
            j = MINLOC(pdat(4,col_i(i):col_i(i+1)-1),1)
            pdat(5,col_i(i)-1+j) = 1.0D0
            p1 = pdat(1:3,list_i(2))
            p2 = pdat(1:3,list_i(4))
            DO j = col_i(i), col_i(i+1)-1
              pdat(4,j) = CalcPerp2(p1,p2,pdat(1:3,j),t)
            ENDDO
            j = MINLOC(pdat(4,col_i(i):col_i(i+1)-1),1)
            pdat(5,col_i(i)-1+j) = 1.0D0
          ENDDO       

          DO i = 1, 100 ! ndat-100, ndat
            WRITE(0,'(A,I6,3F10.3,F12.6,F5.1)') 
     .        'pdat:',i,pdat(1:3,i),pdat(4:5,i)
          ENDDO

          WRITE(0,*) 'avg :',SNGL(avg(1:3))

          WRITE(0,*) 'dist:',list_i(1:4)

          col_min = 1E6
          DO i = 1, col_n-1
            col_min = MIN(col_i(i+1)-col_i(i),col_min)
          ENDDO
          

          DO i = 1, col_n-1

            WRITE(0,'(A,3I6)') 
     .        'col:',i,col_i(i),col_i(i+1)-col_i(i)
          ENDDO




          ntri = 0
          dtri = 0
          stri = 0
          dsid = 0


c          focus(1:3) = pdat(1:3,2)

          cnt = 0

          i1 = 200
          i2 = 201
          i3 = 0
c          focus(1:3) = pdat(1:3,i1)
          focus(1:3) = (pdat(1:3,i1) + pdat(1:3,i2)) / 2.0D0

          min_dist = 0.10D0

c         Set the vector 'normal' for the panel:
          v1 = pdat(1:3,list_i(2)) - pdat(1:3,list_i(1))
          v2 = pdat(1:3,list_i(3)) - pdat(1:3,list_i(1))
          CALL gmCalcCrossProduct(v1,v2,v_panel)
          len1 = DSQRT( v_panel(1)**2 + v_panel(2)**2 + v_panel(3)**2 )
          v_panel= v_panel / len1  ! Normalize

          v1 = pdat(1:3,list_i(1)) - pdat(1:3,list_i(3))
          v1 = v1 / DSQRT( v1(1)**2 + v1(2)**2 + v1(3)**2 )
          v1 = 0.5D0 * (v1 + v_panel)
          v1 = v1 / DSQRT( v1(1)**2 + v1(2)**2 + v1(3)**2 )
          v_left = v1

          v1 = pdat(1:3,list_i(3)) - pdat(1:3,list_i(1))
          v1 = v1 / DSQRT( v1(1)**2 + v1(2)**2 + v1(3)**2 )
          v1 = 0.5D0 * (v1 + v_panel)
          v1 = v1 / DSQRT( v1(1)**2 + v1(2)**2 + v1(3)**2 )
          v_right = v1

          WRITE(0,*) v_panel
          WRITE(0,*) v_left
          WRITE(0,*) v_right

          DO WHILE (.TRUE.)

            cnt = cnt + 1

            pdat(4,1:ndat) = SQRT( (focus(1) - pdat(1,1:ndat))**2 + 
     .                             (focus(2) - pdat(2,1:ndat))**2 + 
     .                             (focus(3) - pdat(3,1:ndat))**2 )

c           Set correct vector:
            v_norm = v_panel
            col_cnt = 0
            DO i = 1, col_n-1
              IF ((col_i(i+1)-col_i(i)).GT.col_min) THEN
                col_cnt = col_cnt+1
                IF (i1.GE.col_i(i).AND.i1.LT.col_i(i+1).AND.
     .              i2.GE.col_i(i).AND.i2.LT.col_i(i+1)) THEN
                  IF (col_cnt.EQ.1.OR.col_cnt.EQ.3) THEN
                    WRITE(0,*) 'activating left!'
                    v_norm = v_left
                  ELSE
                    WRITE(0,*) 'activating right!'
                    v_norm = v_right
                  ENDIF
                  EXIT
                ENDIF
              ENDIF
            ENDDO              


c           Collect nearby points to use when trying to make a new triangle:
            dist_n = 0
            j = 0
            DO WHILE (.TRUE.) 
              j = j + 1
              IF (j.EQ.ndat+1) EXIT
 
              i = MINLOC(pdat(4,1:ndat),1) 

c              IF (j.LT.10) WRITE(0,*) '->',i,j,pdat(4,i)

              IF (pdat(4,i).LT.min_dist.AND.
     .            .NOT.(i.EQ.i1.OR.i.EQ.i2.OR.i.EQ.i3)) THEN  ! Can't duplicate any points on the reference) 
                dist_n = dist_n + 1
                dist_list(dist_n) = i
                IF (dist_n.EQ.10) EXIT
              ENDIF

              pdat(4,i) = 1.0D+6
            ENDDO

            IF (dist_n.EQ.0) 
     .        CALL ER('rayLoadITERFWP','Unable to find nearby '//
     .                'points',*99)

            WRITE(0,*)
            WRITE(0,'(A,I6,2X,10I6)') '>',ntri,dist_list(1:dist_n)

c           Try to make a triangle:
            DO tcnt = 1, dist_n

              status = .FALSE.

              dtri(1,ntri+1) = i1
              dtri(2,ntri+1) = i2
              dtri(3,ntri+1) = dist_list(tcnt)            

              p1 = pdat(1:3,dtri(1,ntri+1))
              p2 = pdat(1:3,dtri(2,ntri+1))
              p3 = pdat(1:3,dtri(3,ntri+1))
              
              len1 = DSQRT ( (p3(1) - p2(1))**2 + (p3(2) - p2(2))**2 + 
     .                       (p3(3) - p2(3))**2 )
              len2 = DSQRT ( (p1(1) - p3(1))**2 + (p1(2) - p3(2))**2 + 
     .                       (p1(3) - p3(3))**2 )
              len3 = DSQRT ( (p2(1) - p1(1))**2 + (p2(2) - p1(2))**2 + 
     .                       (p2(3) - p1(3))**2 )
              
              angle1 = (len2**2 + len3**2 - len1**2) / (2.0D0*len2*len3)
              angle2 = (len1**2 + len3**2 - len2**2) / (2.0D0*len1*len3)
              angle3 = (len1**2 + len2**2 - len3**2) / (2.0D0*len1*len2)
              angle1 = DACOS(angle1) * 180.0D0 / 3.151492D0
              angle2 = DACOS(angle2) * 180.0D0 / 3.151492D0
              angle3 = DACOS(angle3) * 180.0D0 / 3.151492D0
          
              WRITE(0,'(A,I6,3F10.2,2X,3I6,2X,3F10.4)') 
     .            ' ANG:',tcnt,angle1,angle2,angle3,
     .                    dtri(1:3,ntri+1),
     .                    len1,len2,len3          


c             
              IF (angle1.GT.170.0D0.OR.angle2.GT.170.0D0.OR.
     .            angle3.GT.170.0D0) THEN
                WRITE(0,*) '    ---> angle to large'
                status = .TRUE.
              ENDIF

c...          Check if the triangle already exists:
              IF (.NOT.status) THEN
                stri(ntri+1) = SUM(dtri(1:3,ntri+1))
                DO i = 1, ntri
                  IF (stri(ntri+1).EQ.stri(i)) THEN
                    match = 0
                    DO j = 1, 3
                      DO k = 1, 3
                        IF (dtri(j,ntri+1).EQ.dtri(k,i)) match = match+1
                      ENDDO
                    ENDDO
                    IF (match.EQ.3) status = .TRUE.
                  ENDIF
                ENDDO
              ENDIF

              IF (.NOT.status) THEN  ! Check that centre point is on the correct side
                v1 = pdat(1:3,dtri(1,ntri+1)) - 
     .               pdat(1:3,dtri(3,ntri+1))
                v2 = pdat(1:3,dtri(2,ntri+1)) - 
     .               pdat(1:3,dtri(3,ntri+1))
                CALL gmCalcCrossProduct(v1,v2,v3)
c...            Get the angle between V3 and V_NORM:
                dotprod = v3(1) * v_norm(1) + 
     .                    v3(2) * v_norm(2) + 
     .                    v3(3) * v_norm(3)
                len3 = DSQRT(v3(1)**2 + v3(2)**2 + v3(3)**2)
                angle1 = dotprod / (len3) ! V_NORM is normalized
 
c                WRITE(0,*) 'angle:',angle1
c                WRITE(0,*) '     :',dtri(1:3,ntri+1)
c                WRITE(0,*) '     :',pdat(1,dtri(1:3,ntri+1))
c                WRITE(0,*) '     :',pdat(2,dtri(1:3,ntri+1))
c                WRITE(0,*) '     :',pdat(3,dtri(1:3,ntri+1))
c                WRITE(0,*) '  dot:',dotprod
c                WRITE(0,*) '     :',v1
c                WRITE(0,*) '     :',v2
c                WRITE(0,*) '     :',v3
c                WRITE(0,*) '     :',v_norm

                WRITE(0,*) 'angle:',DACOS(angle1)*180.0D0/3.141592D0

                angle1 = DACOS(angle1) * 180.0D0 / 3.141592D0

                IF (angle1.GT.89.9D0) THEN
                  WRITE(0,*) '    ---> facing wrong way'
                  status = .TRUE.
                ENDIF
              ENDIF

              IF (.NOT.status) THEN  ! Triangle accepted
                ntri = ntri + 1
                IF (ntri.GT.1) THEN
                  dsid(1,i) = 1  
                  status = .FALSE.
                  DO i = 1, ntri-1
                    DO j = 1, 3
                      k = j + 1
                      IF (k.EQ.4) k = 1
                      IF ((dtri(2,ntri).EQ.dtri(j,i).AND.
     .                     dtri(3,ntri).EQ.dtri(k,i)).OR.
     .                    (dtri(2,ntri).EQ.dtri(k,i).AND.
     .                     dtri(3,ntri).EQ.dtri(j,i))) dsid(2,ntri) = 1
                      IF ((dtri(3,ntri).EQ.dtri(j,i).AND.
     .                     dtri(1,ntri).EQ.dtri(k,i)).OR.
     .                    (dtri(3,ntri).EQ.dtri(k,i).AND.
     .                     dtri(1,ntri).EQ.dtri(j,i))) dsid(3,ntri) = 1
                    ENDDO
                  ENDDO
                ENDIF

                DO i = 1, 3
                  j = i + 1
                  IF (j.EQ.4) j = 1                
                  IF (pdat(5,dtri(i,ntri)).EQ.1.0D0.AND.
     .                pdat(5,dtri(j,ntri)).EQ.1.0D0) dsid(i,ntri) = 1
                ENDDO


                WRITE(0,'(A,3I6,2X,3I6)') 
     .            ' GOD:',dtri(1:3,ntri),dsid(1:3,ntri)
                WRITE(0,'(A,3F6.1)') 
     .            '    :',pdat(5,dtri(1:3,ntri))

                EXIT


              ELSE
c                WRITE(0,'(A,3I6)') ' BAD:',stri(ntri+1)
              ENDIF 

            ENDDO

            IF (cnt.EQ.200) EXIT

            IF (tcnt.EQ.dist_n+1) THEN
              WRITE(0,*) 'STOPPING BEFCASJG OF'
              STOP
            ENDIF

c           Exit condition:

            status = .FALSE.
            DO i = 1, ntri
              DO j = 1, 3
                IF (dsid(j,i).EQ.0) THEN
                  status = .TRUE.
                  EXIT
                ENDIF
              ENDDO
              IF (status) EXIT
            ENDDO

            IF (j.EQ.4.AND.i.EQ.ntri+1) THEN
              WRITE(0,*) 'NO MORE TRIANGLES TO MAKE',ntri
              EXIT  ! Exit condition...
            ENDIF

            dsid(j,i) = 1  ! Mark this side as processed

            WRITE(0,'(A,2I6)') ' SET:',i,j

            j = j + 1
            IF (j.EQ.4) j = 1
            i1 = dtri(j,i)
            j = j - 1
            IF (j.EQ.0) j = 3
            i2 = dtri(j,i)
            j = j - 1
            IF (j.EQ.0) j = 3
            i3 = dtri(j,i)


            focus(1:3) = 0.5D0 * (pdat(1:3,i1) + pdat(1:3,i2))

            WRITE(0,'(A,3I6,2X,3F10.4)') '    :',i1,i2,i3,focus


          ENDDO

c          stop 'sdgfsdfsd'



          DO angle = 0.0D0, 359.0D0, 360.0D0 

            DO i = 1, ntri
              newvtx(1:3,1) = pdat(1:3,dtri(1,i))
              newvtx(1:3,2) = pdat(1:3,dtri(2,i))
              newvtx(1:3,3) = pdat(1:3,dtri(3,i))

c             Rotate about y-axis (swing):
              CALL RotateVertices(angle,newvtx(1,1),3)
              WRITE(0,'(I6,3F12.4,2X,3I6)') 
     .                i,pdat(1:3,dtri(1,i)),
     .                    dtri(1:3,i)
              CALL AddPolygon(SP_PLANAR_POLYGON,newvtx(1,1),3)
            ENDDO

          ENDDO

          IF (nobj+1.GT.MAX3D) 
     .      CALL ER('LoadVesselStructures','Insufficient array '//
     .              'bounds for all objects',*99)    
          IF (istart.GT.nsrf) THEN
            WRITE(0,*) 'LoadVesselStructures: Strange, no objects'
            RETURN
          ENDIF

          nobj = nobj + 1
          WRITE(0,*) 'VESSEL STRUCTURE IOBJ:',nobj,istart,nsrf
          obj(nobj)%index       = ielement  ! nobj
          obj(nobj)%type        = OP_EMPTY
          obj(nobj)%mode        = 0      
          obj(nobj)%surface     = 1      ! SOLID
          obj(nobj)%wedge1      = 0
          obj(nobj)%wedge2      = 0
          obj(nobj)%colour      = 1
          obj(nobj)%orientation = 1      ! CW
          obj(nobj)%ik          = 0
          obj(nobj)%ir          = 0
          obj(nobj)%in          = -1  ! What should this be?
          obj(nobj)%ivolume     = 0
          obj(nobj)%nside       = 1
          obj(nobj)%iside(1,1)  = istart ! Start index of range of surfaces in surface array, from loading code above
          obj(nobj)%iside(1,2)  = nsrf   ! End index of range of surfaces in surface array
          obj(nobj)%gsur(1)     = GT_TD
          obj(nobj)%tsur(1)     = SP_VESSEL_WALL
          obj(nobj)%reflec(1)   = opt%obj_reflec(ielement)
c..       Defunct:
          obj(nobj)%nsur        = 0
          obj(nobj)%ipts(2,1)   = 0
          obj(nobj)%nmap(1)     = 0
c       ----------------------------------------------------------------
        CASE (4)
          pdat = -999.0D0

          fp = 99
          file = opt%obj_fname(ielement)
          OPEN(fp,FILE=TRIM(file),FORM='FORMATTED',STATUS='OLD',
     .         ERR=98)       
          ndat = 0
          DO WHILE (.TRUE.)
            READ(fp,'(A)',END=40,ERR=98) dummy
            IF (dummy(1:1).EQ.'*'.OR.LEN_TRIM(dummy).LT.5) CYCLE
c            WRITE(0,*) 'dummy>',TRIM(dummy)//'<'
            ndat = ndat + 1
            READ(dummy,*) rdum,pdat(1:3,ndat)
          ENDDO
 40       CONTINUE
          CLOSE(fp)
          WRITE(0,*) 'Ndat:',ndat

c...      Sort the data:
          WRITE(0,*) 'SORTING 1...'
          cont = .TRUE.
          DO WHILE(cont)
            cont = .FALSE.
            DO i = 1, ndat-1
c              IF (pdat(1,i).GT.pdat(1,i+1)) THEN 
              IF (pdat(2,i).GT.pdat(2,i+1)) THEN 
                pdat(:,0  ) = pdat(:,i  )
                pdat(:,i  ) = pdat(:,i+1)
                pdat(:,i+1) = pdat(:,0  )
                cont = .TRUE. 
              ENDIF
            ENDDO
          ENDDO

c LEFT OFF
c well seems I've just wasted time with all this since I can't load BM7-10 since they
c don't follow the same pattern, i.e. they don't seem to have any pattern at all
c will scrap all this and try creating the panels from the equations...

c         Define the back-plane of the panel:          

c         Find average point for the tile:
          avg(1) = SUM(pdat(1,1:ndat)) / DBLE(ndat)
          avg(2) = SUM(pdat(2,1:ndat)) / DBLE(ndat)
          avg(3) = SUM(pdat(3,1:ndat)) / DBLE(ndat)
          pdat(4,1:ndat) = DSQRT( (avg(1) - pdat(1,1:ndat))**2 + 
     .                            (avg(2) - pdat(2,1:ndat))**2 + 
     .                            (avg(3) - pdat(3,1:ndat))**2 )

          pdat(5,:) = 0.0D0  ! Used to mark boundary points

          DO i = 1, 100
            WRITE(0,'(A,I6,3F10.3,F12.6)') 'pdat:',i,pdat(1:4,i)
          ENDDO


          list_n = 0
          DO i = 1, 4
            j = MAXLOC(pdat(4,1:ndat),1) 
            pdat(4,j) = -1.0D6
            pdat(5,j) =  1.0D0
            list_n         = list_n + 1
            list_i(list_n) = j
          ENDDO       
c         Sort:          
          DO i = 1, list_n-1
            DO j = i+1, list_n
              IF (list_i(i).GT.list_i(j)) THEN
                list_i(0) = list_i(i)
                list_i(i) = list_i(j)
                list_i(j) = list_i(0)
              ENDIF
            ENDDO
          ENDDO
c         Store the corner points since LIST_I will be messed up after this:
          DO i = 1, 4
            corner(1:3,i) = pdat(1:3,list_i(i))
          ENDDO

c         The corners seem to be ordered properly, but there's no guarantee that this will always be the case...
c          DO j = 1, 4
c            i = list_i(j)
c            WRITE(0,'(A,I6,3F10.3,F12.6)') 'pdat:',i,pdat(1:4,i)
c          ENDDO
c          STOP

c         Divide all the points into column groups:
          check = -999.0 
          col_n = 0
          DO i = 1, ndat
            IF (DABS(pdat(2,i)-check).GT.1.0D-6) THEN 
              col_n = col_n + 1
              col_i(col_n) = i
              check = pdat(2,i)
            ENDIF
          ENDDO
          col_n = col_n + 1
          col_i(col_n) = ndat + 1
       
          DO i = 1, col_n-1
            WRITE(0,'(A,4I6)') 
     .        'col:',i,col_i(i),col_i(i+1),col_i(i+1)-col_i(i)
          ENDDO
          STOP 'sdfsdfsd'
 
          WRITE(0,*) 'SORTING 2...'
c          cont = .TRUE.
c          DO WHILE(cont)
c            cont = .FALSE.
c            DO i = 1, ndat-1
c              IF (pdat(3,i).GT.pdat(3,i+1).AND.
c     .            pdat(2,i).EQ.pdat(2,i+1)) THEN 
c                pdat(:,0  ) = pdat(:,i  )
c                pdat(:,i  ) = pdat(:,i+1)
c                pdat(:,i+1) = pdat(:,0  )
c                cont = .TRUE.
c              ENDIF
c            ENDDO
c          ENDDO


c         Bottom edge first:

          min_dist = 0.001D0
          p1(1:3) = corner(1:3,1)  ! Bottom left corner (facing the panel)
          p2(1:3) = corner(1:3,3)  ! Bottom right
          CALL MarkEdge(p1,p2,pdat,ndat,min_dist)  
          DO icol = 1, col_n-1
            i1 = col_i(icol )
            i2 = col_i(icol+1) - 1
            DO i = i1, i2-1
              DO j = i+1, i2
                IF (pdat(4,i).GT.pdat(4,j)) THEN
                  pdat(:,0) = pdat(:,i)
                  pdat(:,i) = pdat(:,j)
                  pdat(:,j) = pdat(:,0)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
c         Top edge:
          min_dist = 0.001D0
          p1(1:3) = corner(1:3,2)  ! Top left corner
          p2(1:3) = corner(1:3,4)  ! Top left
          CALL MarkEdge(p1,p2,pdat,ndat,min_dist)  
          DO icol = 1, col_n-1
            i1 = col_i(icol )
            i2 = col_i(icol+1) - 1
            DO i = (i1+i2)/2, i2-1  ! Only sort the top half of the column
              DO j = i+1, i2
                IF (pdat(4,i).LT.pdat(4,j)) THEN
                  pdat(:,0) = pdat(:,i)
                  pdat(:,i) = pdat(:,j)
                  pdat(:,j) = pdat(:,0)
                ENDIF
              ENDDO
            ENDDO
          ENDDO

c         Sort the columns that need sorting:
          col_max = 0
          DO i = 1, col_n-1
            col_max = MAX(col_max,col_i(i+1)-col_i(i))
          ENDDO
          WRITE(0,*) 'col_max:',col_max

          cnt = 0
          DO icol = col_n-1, 1, -1
            i1 = col_i(icol )
            i2 = col_i(icol+1) - 1
            IF (i2-i1+1.LT.col_max) CYCLE

            p1 = pdat(1:3,i1)
            p2 = pdat(1:3,i2)
            CALL MarkEdge(p1,p2,pdat,ndat,min_dist)  ! just do all points for now...

            step  = 0.005
            new_s = i1

            cnt = cnt + 1
            IF (MOD(cnt,2).EQ.1) THEN 
              WRITE(0,*) 'backward'
              d1 = 10.5D0*step
              d2 = -0.5D0*step
              d3 = -step
            ELSE
              WRITE(0,*) 'forward'
              d1 = -0.5D0*step
              d2 = 10.5D0*step
              d3 = step
            ENDIF

            j = 0
            DO dist = d1, d2, d3
              new_n = 0
              DO i = i1, i2
                IF (pdat(4,i).GE.dist.AND.pdat(4,i).LT.dist+step) THEN
                  new_n = new_n + 1
                  new_i(new_n) = i
                ENDIF
              ENDDO
              IF (new_n.GT.0) THEN
                j = j + 1

                DO i = 1, new_n
                  pdat(:,0        ) = pdat(:,new_s+i-1)
                  pdat(:,new_s+i-1) = pdat(:,new_i(i) )
                  pdat(:,new_i(i) ) = pdat(:,0        )
                ENDDO


                WRITE(0,*) 'pass...',new_n+new_s,i2+1
                IF (new_n+new_s.NE.i2+1) THEN  ! Don't add a point if this is the last series of points
                  col_n = col_n + 1
                  DO i = col_n, icol+2+(j-1), -1
                    col_i(i) = col_i(i-1)
                  ENDDO
                  col_i(i) = col_i(i-1) + new_n  ! -j
                ENDIF

                WRITE(0,*) 'pass...',j,cnt,new_s,new_n
                DO i = i1, i2
                  WRITE(0,'(A,2I6,3F12.4)') 
     .              'DIST:',i,new_n,pdat(4,i),dist,dist+step
                ENDDO

                new_s = new_s + new_n
              ENDIF              


            ENDDO

c            IF (cnt.EQ.2) EXIT

          ENDDO          

          DO k = 1, col_n-1
            WRITE(0,'(A,4I6)') 
     .        'col:',k,col_i(k),col_i(k+1),col_i(k+1)-col_i(k)
          ENDDO




c         Sort again because the above loop can mess things up:

c         Bottom edge first:
          min_dist = 0.001D0
          p1(1:3) = corner(1:3,1)  ! Bottom left corner (facing the panel)
          p2(1:3) = corner(1:3,3)  ! Bottom right
          CALL MarkEdge(p1,p2,pdat,ndat,min_dist)  
          DO icol = 1, col_n-1
            i1 = col_i(icol )
            i2 = col_i(icol+1) - 1
            DO i = i1, i2-1
              DO j = i+1, i2
                IF (pdat(4,i).GT.pdat(4,j)) THEN
                  pdat(:,0) = pdat(:,i)
                  pdat(:,i) = pdat(:,j)
                  pdat(:,j) = pdat(:,0)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
c         Top edge:
          min_dist = 0.001D0
          p1(1:3) = corner(1:3,2)  ! Top left corner
          p2(1:3) = corner(1:3,4)  ! Top left
          CALL MarkEdge(p1,p2,pdat,ndat,min_dist)  
          DO icol = 1, col_n-1
            i1 = col_i(icol )
            i2 = col_i(icol+1) - 1
            DO i = (i1+i2)/2, i2-1  ! Only sort the top half of the column
              DO j = i+1, i2
                IF (pdat(4,i).LT.pdat(4,j)) THEN
                  pdat(:,0) = pdat(:,i)
                  pdat(:,i) = pdat(:,j)
                  pdat(:,j) = pdat(:,0)
                ENDIF
              ENDDO
            ENDDO
          ENDDO


          min_dist = 0.001D0
          p1(1:3) = corner(1:3,1)  ! Bottom left corner (facing the panel)
          p2(1:3) = corner(1:3,3)  ! Bottom right
          CALL MarkEdge(p1,p2,pdat,ndat,min_dist)  
          DO i = 1, 50
            WRITE(0,'(A,I6,3F10.3,F12.6)') 'pdat:',i,pdat(1:4,i)
          ENDDO


c...      Swap Y and Z to match RAY coordinate system:
          pdat(5,:) = pdat(3,:)
          pdat(3,:) = pdat(2,:)
          pdat(2,:) = pdat(5,:)

          DO angle = 0.0D0, 359.0D0, 360.0D0 

            DO j = 1, col_n-2
              i1 = col_i(j )
              i2 = col_i(j+1) - 1
              i3 = col_i(j+1)
              i4 = col_i(j+2) - 1

                write(0,*) j,i1,i2,i3,i4

              ibegin = 0
              ishift = 0
              IF (ABS(MOD((i2-i1)-(i4-i3),2)).NE.0) 
     .          CALL ER('rayLoadITERFWP','Odd point differential',*99)
              IF     (i2-i1.LT.i4-i3) THEN
                ishift = ((i4-i3) - (i2-i1)) / 2
              ELSEIF (i2-i1.GT.i4-i3) THEN
                ibegin = ((i2-i1) - (i4-i3)) / 2
              ENDIF

              WRITE(0,*) 'ibegin,ishift',ibegin,ishift

              DO i = ibegin, i2-i1-1-ibegin
                IF (ibegin.NE.0.OR.ishift.NE.0) 
     .            write(0,*) j,i,i1+i,i3+i+ishift
                newvtx(1:3,1) = pdat(1:3,i1+i                )
                newvtx(1:3,2) = pdat(1:3,i1+i+1              )
                newvtx(1:3,3) = pdat(1:3,i3+i+1+ishift-ibegin)
                newvtx(1:3,4) = pdat(1:3,i3+i  +ishift-ibegin)
c               Rotate about y-axis (swing):
                CALL RotateVertices(angle,newvtx(1,1),4)
                CALL AddPolygon(SP_PLANAR_POLYGON,newvtx(1,1),4)
              ENDDO

            ENDDO

          ENDDO

          IF (nobj+1.GT.MAX3D) 
     .      CALL ER('LoadVesselStructures','Insufficient array '//
     .              'bounds for all objects',*99)    
          IF (istart.GT.nsrf) THEN
            WRITE(0,*) 'LoadVesselStructures: Strange, no objects'
            RETURN
          ENDIF

          nobj = nobj + 1
          WRITE(0,*) 'VESSEL STRUCTURE IOBJ:',nobj,istart,nsrf
          obj(nobj)%index       = ielement  ! nobj
          obj(nobj)%type        = OP_EMPTY
          obj(nobj)%mode        = 0      
          obj(nobj)%surface     = 1      ! SOLID
          obj(nobj)%wedge1      = 0
          obj(nobj)%wedge2      = 0
          obj(nobj)%colour      = 1
          obj(nobj)%orientation = 1      ! CW
          obj(nobj)%ik          = 0
          obj(nobj)%ir          = 0
          obj(nobj)%in          = -1  ! What should this be?
          obj(nobj)%ivolume     = 0
          obj(nobj)%nside       = 1
          obj(nobj)%iside(1,1)  = istart ! Start index of range of surfaces in surface array, from loading code above
          obj(nobj)%iside(1,2)  = nsrf   ! End index of range of surfaces in surface array
          obj(nobj)%gsur(1)     = GT_TD
          obj(nobj)%tsur(1)     = SP_VESSEL_WALL
          obj(nobj)%reflec(1)   = opt%obj_reflec(ielement)
c..       Defunct:
          obj(nobj)%nsur        = 0
          obj(nobj)%ipts(2,1)   = 0
          obj(nobj)%nmap(1)     = 0
c       ----------------------------------------------------------------
        CASE DEFAULT
          CALL ER('LoadITERFWP','Unknown geometry sub-option',*99)
      ENDSELECT


      RETURN
 98   CALL ER('LoadITERFWP','File not found',*99)
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE LoadImageReconstruction(ielement)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

c...  Input:
      INTEGER, INTENT(IN) :: ielement

      INTEGER AddVertex,AddSurface

      INTEGER   fp,count,istart,idum1,i1,i2,ivol,fobj,
     .          iobj,iobj2,iside,iside2,isrf,isrf2
      LOGICAL   map_found
      REAL*8    val,newvtx(3,4),v1(3),v2(3),v3(3),v4(3)
      CHARACTER file*1024,buffer*2048
      TYPE(type_surface) newsrf

      REAL*8, PARAMETER :: DTOL = 1.0D-10

      fp = 99


      WRITE(0,*) 'LOADING IMAGE RECONSTRUCTIONS'//
     .  opt%obj_fname(ielement)


c...  Load up surface and vertex arrays:

      istart = nsrf + 1
      fobj   = nobj + 1
      ivol   = 0

      file = opt%obj_fname(ielement)

      OPEN(fp,FILE=TRIM(file),FORM='FORMATTED',STATUS='OLD',ERR=98)     
 
      count = 0

      READ(fp,*)  ! Clear index line

      DO WHILE (.TRUE.)       

        READ(fp,'(A2048)',END=10) buffer         

        newvtx = 0.0D0
        READ(buffer,*,ERR=10) 
     .    idum1,val,idum1,(newvtx(1:2,i1),i1=1,idum1)

c        newvtx(1:2,1) = newvtx(1:2,1) + 0.0000005
c        newvtx(1:2,2) = newvtx(1:2,2) + 0.0000005
c        newvtx(1:2,3) = newvtx(1:2,3) + 0.0000005
c        newvtx(1:2,4) = newvtx(1:2,4) + 0.0000005

        IF (nobj+1.GT.MAX3D) 
     .    CALL ER('LoadImageRecon','Insufficient array bounds '//
     .            'for all objects',*99)    

        nobj = nobj + 1
        ivol = ivol + 1

c...    Extract vertices:
c        WRITE(0,*) 'IMAGE RECON. NOBJ:',nobj
c        WRITE(0,*) '   :',idum1,val

        obj(nobj)%index       = ielement  ! nobj
        obj(nobj)%type        = OP_INTEGRATION_VOLUME
        obj(nobj)%subtype     = OP_INVERSION_GRID  ! OK for now...
        obj(nobj)%mode        = 0      
        obj(nobj)%surface     = 1      ! SOLID ???
        obj(nobj)%wedge1      = 0
        obj(nobj)%wedge2      = 0
        obj(nobj)%colour      = 3
        obj(nobj)%orientation = 1      ! CW
        obj(nobj)%ik          = 0
        obj(nobj)%ir          = 0
        obj(nobj)%in          = 0
        obj(nobj)%ivolume     = ivol
        obj(nobj)%nside       = idum1
        obj(nobj)%reflec      = 0
        obj(nobj)%quantity(1) = 1.0D0 ! val ! / 1.0E+18 ! 1.0 ! val


c        WRITE(0,*) 'quantit:',nobj,obj(nobj)%quantity(1)
c..     Defunct:
        obj(nobj)%nsur        = 0
        obj(nobj)%ipts(2,1)   = 0
        obj(nobj)%nmap(1)     = 0

        DO i1 = 1, obj(nobj)%nside
          i2 = i1 + 1
          IF (i2.GT.obj(nobj)%nside) i2 = 1  
          newsrf%type = -1
          newsrf%nvtx =  2
          newsrf%ivtx(1) = AddVertex(newvtx(1,i1))
          newsrf%ivtx(2) = AddVertex(newvtx(1,i2))
          obj(nobj)%iside(i1,1) = AddSurface(newsrf)
          obj(nobj)%iside(i1,2) = obj(nobj)%iside(i1,1)
          obj(nobj)%gsur(i1)    = GT_TC
          obj(nobj)%tsur(i1)    = -1
c          WRITE(0,*) '   -',i1,obj(nobj)%iside(i1,1)
        ENDDO

      ENDDO
 10   CONTINUE
      CLOSE(fp)


      IF (.TRUE.) THEN
c...    Force rebuild of connection map:

        DO iobj = fobj, nobj
          IF (MOD(iobj-fobj+1,1000).EQ.0) 
     .      WRITE(0,*) '  BUILDING CONNECTION MAP',iobj,nobj-fobj+1
c          IF (iobj.EQ.5) STOP ' debugging...'
          DO iside = 1, obj(iobj)%nside
            isrf = obj(iobj)%iside(iside,1)     
            v1(1:3) = vtx(1:3,srf(isrf)%ivtx(1))
            v2(1:3) = vtx(1:3,srf(isrf)%ivtx(2))

            obj(iobj)%tsur(iside) = SP_GRID_BOUNDARY
            obj(iobj)%imap(1,iside) = iobj
            obj(iobj)%isur(1,iside) = iside

c            WRITE(0,*) '  IOBJ :',iobj,iside,isrf
c            WRITE(0,*) '       :',v1(1:2)
c            WRITE(0,*) '       :',v2(1:2)

            map_found = .FALSE. 
            DO iobj2 = MAX(fobj,iobj-1050),MIN(nobj,iobj+1050) ! Crude...
              IF (iobj.EQ.iobj2) CYCLE
              DO iside2 = 1, obj(iobj2)%nside
                isrf2 = obj(iobj2)%iside(iside2,1)     
                v3(1:3) = vtx(1:3,srf(isrf2)%ivtx(1))
                v4(1:3) = vtx(1:3,srf(isrf2)%ivtx(2))
c                WRITE(0,*) '  IOBJ2:',iobj2,iside2,isrf2
c                WRITE(0,*) '       :',v1(1:2)
c                WRITE(0,*) '       :',v2(1:2)
c                WRITE(0,*) '       :',v3(1:2)
c                WRITE(0,*) '       :',v4(1:2)
                IF ((DABS(v1(1)-v3(1)).LT.DTOL.AND.
     .               DABS(v1(2)-v3(2)).LT.DTOL.AND.
     .               DABS(v2(1)-v4(1)).LT.DTOL.AND.
     .               DABS(v2(2)-v4(2)).LT.DTOL).OR.
     .              (DABS(v1(1)-v4(1)).LT.DTOL.AND.
     .               DABS(v1(2)-v4(2)).LT.DTOL.AND.
     .               DABS(v2(1)-v3(1)).LT.DTOL.AND.
     .               DABS(v2(2)-v3(2)).LT.DTOL)) THEN
                  obj(iobj)%tsur(iside) = SP_GRID_SURFACE
                  obj(iobj)%imap(1,iside) = iobj2
                  obj(iobj)%isur(1,iside) = iside2
                  map_found = .TRUE.
c                  WRITE(0,*) '  MAP:',iobj,iside,iobj2,iside2
                  EXIT
                ENDIF
              ENDDO
              IF (map_found) EXIT
            ENDDO 
c            WRITE(0,*) '  :',iobj,iside,
c     .        obj(iobj)%imap(1,iside),obj(iobj)%isur(1,iside)
c            IF (iobj.EQ.5) STOP 'sdfsdf'
          ENDDO
        ENDDO
c          WRITE(0,*) 'CEL:',iobj,obj(iobj)%tsur(1:4)
c          WRITE(0,*) '   :',iobj,obj(iobj)%imap(1,1:4)
c          WRITE(0,*) '   :',iobj,obj(iobj)%isur(1,1:4)
      ENDIF

      WRITE(0,*) 'DONE'

      WRITE(6,*) 'DESPERATE'
      DO iobj = 1, nobj
        DO iside = 1, 4
          isrf = obj(iobj)%iside(iside,1)     
          WRITE(6,*) iobj,iside,vtx(1:2,srf(isrf)%ivtx(1))
          WRITE(6,*) iobj,iside,vtx(1:2,srf(isrf)%ivtx(2))
c          obj(iobj)%v(1,iside) = vtx(1,srf(isrf)%ivtx(1))
c          obj(iobj)%v(2,iside) = vtx(2,srf(isrf)%ivtx(1))
c          obj(iobj)%v(3,iside) = 0.0D0
c          obj(nobj)%npts(iside) = 2
c          obj(nobj)%ipts(1,iside) = 1
c          obj(nobj)%ipts(2,iside) = 2
c          obj(nobj)%nmap(iside)   = 1
        ENDDO
c        obj(nobj)%nside = 0
c        obj(nobj)%nsur  = 4
      ENDDO
c      nvtx = 0
c      nsrf = 0


      RETURN
 98   WRITE(0,*) 'ERROR LoadImageReconstruction: File not found'
      WRITE(0,*) '   '//TRIM(file)
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE LoadLineSegmentFile(ielement)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

c...  Input:
      INTEGER ielement

      INTEGER AddVertex,AddSurface

      INTEGER   fp,count,istart,n1,n2,idum1,option
      REAL*8    newvtx(3,2)
      CHARACTER file*1024,buffer*2048
      TYPE(type_surface) newsrf


      fp = 99

      option = opt%obj_option(ielement)

      WRITE(0,*) 'LOADING LINE SEGMENTS'//
     .  opt%obj_fname(ielement)


c...  Load up surface and vertex arrays:

      istart = nsrf + 1

      file = opt%obj_fname(ielement)

      OPEN(fp,FILE=file(1:LEN_TRIM(file)),
     .     FORM='FORMATTED',STATUS='OLD',ERR=98)     
 
      n1 = opt%obj_n(ielement,1)
      n2 = opt%obj_n(ielement,2)

      count = 0

      DO WHILE (.TRUE.)       

        READ(fp,'(A2048)',END=10) buffer         

        IF (buffer(1:1).EQ.'*'.OR.LEN_TRIM(buffer).EQ.0.OR.
     .      buffer(1:1).EQ.'#') THEN
          count = 0
          CYCLE
        ENDIF

c...    Extract vertices:
        newvtx(1:3,1) = 0.0D0

        SELECTCASE (option)
          CASE(1)
            READ(buffer,*,ERR=10,END=10) newvtx(1,1),newvtx(2,1)
          CASE(2)
            READ(buffer,*,ERR=10,END=10) newvtx(2,1),newvtx(1,1)
          CASE(3)
            READ(buffer,*,ERR=10,END=10) newvtx(1:2,1),newvtx(1:2,2)
          CASE DEFAULT
            CALL ER('LoadLineSegmentFile','Unknown orientation '//
     .              'option',*99)
        ENDSELECT

        IF (option.EQ.3.OR.
     .      ((option.EQ.1.OR.option.EQ.2).AND.count.GT.0.AND.
     .       ((count.GE.n1).OR.(-1.EQ.n1)).AND.
     .       ((count.LE.n2).OR.(-1.EQ.n2)))) THEN
          newsrf%type = SP_LINE_SEGMENT  
          newsrf%nvtx = 2
          newsrf%ivtx(1) = AddVertex(newvtx(1,1))
          newsrf%ivtx(2) = AddVertex(newvtx(1,2))
          idum1 = AddSurface(newsrf)
        ENDIF
        newvtx(1:3,2) = newvtx(1:3,1)
        count = count + 1
      ENDDO
 10   CONTINUE
      CLOSE(fp)

c      WRITE(0,*) 'VERTICES? SURFACES?',nvtx,nsrf
c      WRITE(0,*) '                   ',istart

c...  Assign object(s):    

c      WRITE(0,*) 'NOBJ:',nobj,MAX3D

      IF (nobj+1.GT.MAX3D) 
     .  CALL ER('LoadVesselStructures','Insufficient array bounds '//
     .          'for all objects',*99)    

      IF (istart.GT.nsrf) THEN
        WRITE(0,*) 'LoadVesselStructures: Strange, no objects loaded'
        RETURN
      ENDIF

      nobj = nobj + 1
      WRITE(0,*) 'VESSEL STRUCTURE IOBJ:',nobj

      obj(nobj)%index       = ielement  ! nobj
      obj(nobj)%type        = OP_EMPTY
      obj(nobj)%mode        = 0      
      obj(nobj)%surface     = 1      ! SOLID
      obj(nobj)%wedge1      = 0
      obj(nobj)%wedge2      = 0
      obj(nobj)%colour      = MAX(1,opt%obj_colour(ielement)/100)
      obj(nobj)%orientation = 1      ! CW
      obj(nobj)%ik          = 0
      obj(nobj)%ir          = 0
      obj(nobj)%in          = -1  ! What should this be?
      obj(nobj)%ivolume     = 0
      obj(nobj)%nside       = 1
      obj(nobj)%iside(1,1)  = istart ! Start index of range of surfaces in surface array, from loading code above
      obj(nobj)%iside(1,2)  = nsrf   ! End index of range of surfaces in surface array
      obj(nobj)%gsur(1)     = GT_TC
      obj(nobj)%tsur(1)     = SP_VESSEL_WALL
      obj(nobj)%reflec(1)   = opt%obj_reflec(ielement)
c..   Defunct:
      obj(nobj)%nsur        = 0
      obj(nobj)%ipts(2,1)   = 0
      obj(nobj)%nmap(1)     = 0

      WRITE(0,*) 'DONE'

      RETURN
 98   WRITE(0,*) 'ERROR LoadLineSegmentFile: File not found'
      WRITE(0,*) '   '//file(1:LEN_TRIM(file))
 99   STOP
      END
c
c ====================================================================== 
c
c
c             side 3
c         3------------4            1  2  3  4  5 ....   nxbin
c         |            |            nxbin+1 ......     2*nxbin
c         |            |            .... 
c         |            |            ....           nybin*nxbin
c  side 2 |            | side 4
c         |            |
c         |            |
c         |            |
c         2------------1
c             side 1
c
c
c 
c
c
      SUBROUTINE BuildInversionMesh(ielement)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

      INTEGER ielement



      LOGICAL CheckInversionCell

      INTEGER i1,i2,i3,ix,iy,ncells,isector,nsector,nxbin,nybin,ivol,
     .        istart,iend,fobj,iobj,isid,iobj2,ik,ir,iside
      LOGICAL outofgrid
      REAL    angle,dangle,ang
      REAL*8  p1(3,8),p2(3,8),xcen,ycen,xwidth,ywidth,xdelta,ydelta,
     .        x1,z1,xorigin,yorigin,distsep,distwal,
     .        frac,maxdist,maxdiststd,maxdistxpt,
     .        x(4),y(4),x2(4),y2(4)

      REAL, PARAMETER :: RADDEG = 57.29577952

      REAL*8, PARAMETER :: DTOL = 1.0D-10

c...    Number of objects per toroidal segment:

c      WRITE(0,*) 'BUILDING RADIAL MAP'


      nsector = opt%obj_nsector
      IF (nsector.EQ.-1) nsector = grd_ntorseg ! eirntorseg

      IF (nsector.EQ.0) THEN
        WRITE(0,*) 'OVER-RIDING nsector IN BUILDINVERSIONMESH'
        nsector = 48
      ENDIF

      dangle = 360.0 / REAL(nsector) / RADDEG

      fobj = nobj + 1

      IF (ielement.NE.0) THEN
        nxbin = opt%obj_n(ielement,1)
        nybin = opt%obj_n(ielement,2)
        ncells = nxbin * nybin

        xcen = 0.5D0*DBLE(opt%obj_r(ielement,1) + opt%obj_r(ielement,2))
        ycen = 0.5D0*DBLE(opt%obj_z(ielement,1) + opt%obj_z(ielement,2))
c        xcen = 0.5D0 * (opt%obj_r(ielement,1) + opt%obj_r(ielement,2))
c        ycen = 0.5D0 * (opt%obj_z(ielement,1) + opt%obj_z(ielement,2))

        xwidth = DBLE(opt%obj_r(ielement,2) - opt%obj_r(ielement,1))
        ywidth = DBLE(opt%obj_z(ielement,2) - opt%obj_z(ielement,1))
c        xwidth = opt%obj_r(ielement,2) - opt%obj_r(ielement,1)
c        ywidth = opt%obj_z(ielement,2) - opt%obj_z(ielement,1)

        xdelta = xwidth / DBLE(nxbin)
        ydelta = ywidth / DBLE(nybin) 

        xorigin = xcen - 0.5D0 * DBLE(nxbin) * xdelta
        yorigin = ycen - 0.5D0 * DBLE(nybin) * ydelta

        WRITE(6,*) 'XCEN,YCEN',xcen,ycen
        WRITE(6,*) 'delta',xdelta,ydelta
        WRITE(6,*) 'wid',xwidth,ywidth
        WRITE(6,*) 'origin',xorigin,yorigin
         

        IF (opt%obj_option(ielement).EQ.6) nybin = 2 * nybin - 1
      ELSE
c...    Original code:
        STOP 'CODE OBSOLETE'
      ENDIF


      IF ((ielement.NE.0.AND.
     .     (opt%obj_option(ielement).EQ.2.OR.
     .      opt%obj_option(ielement).EQ.4.OR.
     .      opt%obj_option(ielement).EQ.6))) THEN
c...    (Perfect) toroidal symmetry:

        ivol = 0  ! Integration volume index 

        maxdiststd = 0.15D0
        maxdistxpt = 0.25D0

        DO iy = nybin, 1, -1
          DO ix = 1, nxbin

c...        Check if the inversion cell is inside the OSM fluid grid:
            IF (opt%obj_option(ielement).EQ.4) THEN
              xcen = xorigin + xdelta * (DBLE(ix) - 0.5D0)
              ycen = yorigin + ydelta * (DBLE(iy) - 0.5D0)           
              IF (.NOT.CheckInversionCell(1,xcen,ycen)) CYCLE
            ENDIF

            IF (nobj+1.GT.MAX3D) 
     .        CALL ER('BuildObjects','Insufficient array bounds '//
     .                'for all objects',*98)     

            ivol = ivol + 1

            nobj = nobj + 1
            obj(nobj)%index       = ielement  ! nobj
            obj(nobj)%type        = OP_INTEGRATION_VOLUME
            obj(nobj)%subtype     = OP_INVERSION_GRID
            obj(nobj)%mode        = 0      
            obj(nobj)%surface     = 1      ! SOLID
            obj(nobj)%wedge1      = 0
            obj(nobj)%wedge2      = 0
            obj(nobj)%colour      = 3
            obj(nobj)%orientation = 1      ! CW
            obj(nobj)%ik          = ix
            obj(nobj)%ir          = iy
            obj(nobj)%in          = 0
            obj(nobj)%ivolume     = ivol
            obj(nobj)%nsur        = 4
            obj(nobj)%gsur(1:4)   = GT_TC
            obj(nobj)%nver        = 4
            obj(nobj)%reflec(1:4) = 0
 
            obj(nobj)%quantity = 1.0
c            obj(nobj)%quantity(1) = 1.0

c...        Checkerboard:
            IF (.FALSE.) THEN
              IF (MOD(ix,2).EQ.0) THEN
                IF (MOD(iy,2).EQ.0) THEN 
                  obj(nobj)%quantity(1) = 1.0
                ELSE
                  obj(nobj)%quantity(1) = 0.0
                ENDIF
              ELSE
                IF (MOD(iy,2).EQ.0) THEN 
                  obj(nobj)%quantity(1) = 0.0
                ELSE
                  obj(nobj)%quantity(1) = 1.0
                ENDIF
              ENDIF
            ENDIF
c            IF (ix.GT.nxbin/2) obj(nobj)%quantity(1) = 0.0
c            obj(3)%quantity(1) = 10.0
c            obj(nobj)%quantity(1) = REAL(ivol)

            SELECTCASE (opt%obj_option(ielement)) 
              CASE(2,4)  ! Rectangles
c...            Vertices:
                obj(nobj)%v(1,1) = xorigin + xdelta * DBLE(ix)   
                obj(nobj)%v(2,1) = yorigin + ydelta * DBLE(iy-1)
                obj(nobj)%v(3,1) = 0.0D0
                obj(nobj)%v(1,2) = xorigin + xdelta * DBLE(ix-1)
                obj(nobj)%v(2,2) = yorigin + ydelta * DBLE(iy-1)
                obj(nobj)%v(3,2) = 0.0D0
                obj(nobj)%v(1,3) = xorigin + xdelta * DBLE(ix-1)
                obj(nobj)%v(2,3) = yorigin + ydelta * DBLE(iy)
                obj(nobj)%v(3,3) = 0.0D0
                obj(nobj)%v(1,4) = xorigin + xdelta * DBLE(ix)
                obj(nobj)%v(2,4) = yorigin + ydelta * DBLE(iy)
                obj(nobj)%v(3,4) = 0.0D0

c...            Surface 1:
                IF (iy.EQ.1) THEN    
                  obj(nobj)%tsur(1) = SP_GRID_BOUNDARY
                  obj(nobj)%nmap(1) = 1
                  obj(nobj)%imap(1,1) = nobj
                  obj(nobj)%isur(1,1) = 1
                ELSE
                  obj(nobj)%tsur(1) = SP_GRID_SURFACE
                  obj(nobj)%nmap(1) = 1
                  obj(nobj)%imap(1,1) = nobj + nxbin
                  obj(nobj)%isur(1,1) = 3
                ENDIF
          
                obj(nobj)%npts(1) = 2
                obj(nobj)%ipts(1,1) = 1
                obj(nobj)%ipts(2,1) = 2
c...            Surface 2:
                IF (ix.EQ.1) THEN
                  obj(nobj)%tsur(2) = SP_GRID_BOUNDARY
                  obj(nobj)%nmap(2) = 1
                  obj(nobj)%imap(1,2) = nobj
                  obj(nobj)%isur(1,2) = 2
                ELSE
                  obj(nobj)%tsur(2) = SP_GRID_SURFACE
                  obj(nobj)%nmap(2) = 1
                  obj(nobj)%imap(1,2) = nobj - 1
                  obj(nobj)%isur(1,2) = 4
                ENDIF
                obj(nobj)%npts(2) = 2
                obj(nobj)%ipts(1,2) = 2
                obj(nobj)%ipts(2,2) = 3
c...            Surface 4:
                IF (iy.EQ.nybin) THEN   
                  obj(nobj)%tsur(3) = SP_GRID_BOUNDARY
                  obj(nobj)%nmap(3) = 1
                  obj(nobj)%imap(1,3) = nobj
                  obj(nobj)%isur(1,3) = 3
                ELSE
                  obj(nobj)%tsur(3) = SP_GRID_SURFACE
                  obj(nobj)%nmap(3) = 1
                  obj(nobj)%imap(1,3) = nobj - nxbin
                  obj(nobj)%isur(1,3) = 1
                ENDIF
                obj(nobj)%npts(3) = 2
                obj(nobj)%ipts(1,3) = 3
                obj(nobj)%ipts(2,3) = 4
c...            Surface 4:
                IF (ix.EQ.nxbin) THEN
                  obj(nobj)%tsur(4) = SP_GRID_BOUNDARY
                  obj(nobj)%nmap(4) = 1
                  obj(nobj)%imap(1,4) = nobj
                  obj(nobj)%isur(1,4) = 4
                ELSE
                  obj(nobj)%tsur(4) = SP_GRID_SURFACE
                  obj(nobj)%nmap(4) = 1
                  obj(nobj)%imap(1,4) = nobj + 1
                  obj(nobj)%isur(1,4) = 2
                ENDIF
                obj(nobj)%npts(4) = 2
                obj(nobj)%ipts(1,4) = 4
                obj(nobj)%ipts(2,4) = 1

              CASE(6) ! Diamonds
c...            Vertices:
                obj(nobj)%v(1,1) = xorigin + xdelta * (DBLE(ix  )-0.5D0)
                obj(nobj)%v(2,1) = yorigin + ydelta *  DBLE(iy-1)
                obj(nobj)%v(3,1) = 0.0D0
                obj(nobj)%v(1,2) = xorigin + xdelta *  DBLE(ix-1)
                obj(nobj)%v(2,2) = yorigin + ydelta * (DBLE(iy-1)+0.5D0)
                obj(nobj)%v(3,2) = 0.0D0
                obj(nobj)%v(1,3) = xorigin + xdelta * (DBLE(ix-1)+0.5D0)
                obj(nobj)%v(2,3) = yorigin + ydelta *  DBLE(iy  )
                obj(nobj)%v(3,3) = 0.0D0
                obj(nobj)%v(1,4) = xorigin + xdelta *  DBLE(ix  )
                obj(nobj)%v(2,4) = yorigin + ydelta * (DBLE(iy  )-0.5D0)
                obj(nobj)%v(3,4) = 0.0D0

                IF (MOD(iy,2).EQ.0) THEN 
                  obj(nobj)%v(1,1:4) = obj(nobj)%v(1,1:4) + 0.5D0*xdelta   
                ELSE
c                  obj(nobj)%v(2,1:4) = obj(nobj)%v(2,1:4) + 0.5D0*ydelta   
                ENDIF

                obj(nobj)%v(2,1:4) = obj(nobj)%v(2,1:4) - 
     .                               0.5D0*DBLE(iy-1)*ydelta   

                obj(nobj)%npts(1:4) = 2
                obj(nobj)%nmap(1:4) = 1

                obj(nobj)%ipts(1,1) = 1
                obj(nobj)%ipts(2,1) = 2
                obj(nobj)%ipts(1,2) = 2
                obj(nobj)%ipts(2,2) = 3
                obj(nobj)%ipts(1,3) = 3
                obj(nobj)%ipts(2,3) = 4
                obj(nobj)%ipts(1,4) = 4
                obj(nobj)%ipts(2,4) = 1

              CASE DEFAULT
                CALL ER('BuildInversionMesh','Unknown option',*99)
            ENDSELECT  
        
          ENDDO
        ENDDO
        

        IF (opt%obj_option(ielement).EQ.4.OR.
     .      opt%obj_option(ielement).EQ.6) THEN
c...      Force rebuild of connection map:

          DO iobj = fobj, nobj
            DO isid = 1, obj(iobj)%nsur

              x(1:4) = obj(iobj)%v(1,1:4)
              y(1:4) = obj(iobj)%v(2,1:4)

              obj(iobj)%tsur(isid) = SP_GRID_BOUNDARY
              obj(iobj)%imap(1,isid) = iobj
              obj(iobj)%isur(1,isid) = isid
 
              SELECTCASE (isid)
                CASE(1)
                  DO iobj2 = MIN(nobj,iobj+1), MIN(nobj,iobj+nxbin+1) 
                    x2(1:4) = obj(iobj2)%v(1,1:4)
                    y2(1:4) = obj(iobj2)%v(2,1:4)
                    IF (DABS(x(1)-x2(4)).LT.DTOL.AND.
     .                  DABS(y(1)-y2(4)).LT.DTOL.AND.
     .                  DABS(x(2)-x2(3)).LT.DTOL.AND.
     .                  DABS(y(2)-y2(3)).LT.DTOL) THEN
                      obj(iobj)%tsur(isid) = SP_GRID_SURFACE
                      obj(iobj)%imap(1,isid) = iobj2
                      obj(iobj)%isur(1,isid) = 3                   
                      EXIT
                    ENDIF
                  ENDDO

                CASE(2)
                  DO iobj2 = MAX(fobj,iobj-1),MAX(fobj,iobj-nxbin-1),-1
                    x2(1:4) = obj(iobj2)%v(1,1:4)
                    y2(1:4) = obj(iobj2)%v(2,1:4)
                    IF (DABS(x(2)-x2(1)).LT.DTOL.AND.
     .                  DABS(y(2)-y2(1)).LT.DTOL.AND.
     .                  DABS(x(3)-x2(4)).LT.DTOL.AND.
     .                  DABS(y(3)-y2(4)).LT.DTOL) THEN
                      obj(iobj)%tsur(isid) = SP_GRID_SURFACE
                      obj(iobj)%imap(1,isid) = iobj2
                      obj(iobj)%isur(1,isid) = 4                    
                      EXIT
                    ENDIF
                  ENDDO

                CASE(3)
                  DO iobj2 = MAX(fobj,iobj-1),MAX(fobj,iobj-nxbin-1),-1
                    x2(1:4) = obj(iobj2)%v(1,1:4)
                    y2(1:4) = obj(iobj2)%v(2,1:4)
                    IF (DABS(x(3)-x2(2)).LT.DTOL.AND.
     .                  DABS(y(3)-y2(2)).LT.DTOL.AND.
     .                  DABS(x(4)-x2(1)).LT.DTOL.AND.
     .                  DABS(y(4)-y2(1)).LT.DTOL) THEN
                      obj(iobj)%tsur(isid) = SP_GRID_SURFACE
                      obj(iobj)%imap(1,isid) = iobj2
                      obj(iobj)%isur(1,isid) = 1                  
                      EXIT
                    ENDIF
                  ENDDO               

                CASE(4)
                  DO iobj2 = MIN(nobj,iobj+1),MIN(nobj,iobj+nxbin+1)
                    x2(1:4) = obj(iobj2)%v(1,1:4)
                    y2(1:4) = obj(iobj2)%v(2,1:4)
                    IF (DABS(x(1)-x2(2)).LT.DTOL.AND.
     .                  DABS(y(1)-y2(2)).LT.DTOL.AND.
     .                  DABS(x(4)-x2(3)).LT.DTOL.AND.
     .                  DABS(y(4)-y2(3)).LT.DTOL) THEN
                      obj(iobj)%tsur(isid) = SP_GRID_SURFACE
                      obj(iobj)%imap(1,isid) = iobj2
                      obj(iobj)%isur(1,isid) = 2                   
                      EXIT
                    ENDIF
                  ENDDO
                CASEDEFAULT
                  CALL ER('BuildInversionMesh','Too many sides',*99)
              ENDSELECT       

            ENDDO
 
c            WRITE(0,*) 'CEL:',iobj,obj(iobj)%tsur(1:4)
c            WRITE(0,*) '   :',iobj,obj(iobj)%imap(1,1:4)
c            WRITE(0,*) '   :',iobj,obj(iobj)%isur(1,1:4)

          ENDDO

        ENDIF


      ELSEIF ((ielement.NE.0.AND.
     .         opt%obj_option(ielement).EQ.1)) THEN
c     .        opt%ob_invgrd.EQ.1) THEN

c...    Toroidal approximation via discretization:

        istart = 1
        iend   = nsector 

c...    *HACK* (more hacks below)
        WRITE(0,*)
        WRITE(0,*) '--------------------------------------------------'
        WRITE(0,*) ' TOROIDAL INVERSION MESH HACK!'
        WRITE(0,*) '--------------------------------------------------'
        istart = 0
        iend   = 2
        ivol   = 0 

        IF (iend-istart+1.GT.nsector) 
     .    CALL ER('BuildInversionMesh','An excess of toroidal '//
     .            'sectors requested',*99)

        DO isector = istart, iend

          ang = REAL(isector - 1) * dangle

c...      *HACK* 
c          ivol = 0  ! Integration volume index 

          DO iy = nybin, 1, -1
            DO ix = 1, nxbin

              IF (nobj+1.GT.MAX3D) 
     .          CALL ER('BuildObjects','Insufficient array bounds '//
     .                  'for all objects',*98)     

              ivol = ivol + 1

              nobj = nobj + 1
              obj(nobj)%index       = ielement  ! nobj
              obj(nobj)%type        = OP_INTEGRATION_VOLUME
              obj(nobj)%mode        = 0      
              obj(nobj)%surface     = 1      ! SOLID
              obj(nobj)%phi         = ang
              obj(nobj)%wedge1      = 0
              obj(nobj)%wedge2      = 0
              obj(nobj)%colour      = 3
              obj(nobj)%orientation = 1      ! CW
              obj(nobj)%ik          = ix
              obj(nobj)%ir          = iy
              obj(nobj)%in          = 0
              obj(nobj)%ivolume     = ivol
              obj(nobj)%nsur        = 6
              obj(nobj)%gsur(1:6)   = GT_TD
              obj(nobj)%nver        = 8

c             obj(nobj)%quantity = 1.0

c...          *HACK* 
              IF (iy.EQ.nybin/2.AND.
     .            ix.GT.nxbin/3.AND.ix.LT.nxbin-nxbin/3+1.AND.
     .            isector.EQ.1) THEN
                obj(nobj)%quantity(1) = 1.0
              ELSE
                obj(nobj)%quantity(1) = 0.0
              ENDIF

c              IF (MOD(ix,2).EQ.0) THEN
c                IF (MOD(iy,2).EQ.0) THEN 
c                  obj(nobj)%quantity(1) = 1.0
c                ELSE
c                  obj(nobj)%quantity(1) = 0.0
c                ENDIF
c              ELSE
c                IF (MOD(iy,2).EQ.0) THEN 
c                  obj(nobj)%quantity(1) = 0.0
c                ELSE
c                  obj(nobj)%quantity(1) = 1.0
c                ENDIF
c              ENDIF
c              obj(3)%quantity(1) = 10.0
c            obj(nobj)%quantity(1) = REAL(ivol)
c            obj(nobj)%quantity(2) = 1.0

              obj(nobj)%nmap = 0

c...          ...
              p1(1,1)   = xorigin + xdelta * DBLE(ix)   ! Clean this up...
              p1(2,1)   = yorigin + ydelta * DBLE(iy-1)
              p1(3,1)   = DBLE(p1(1,1))*DTAN(DBLE(-0.5*dangle))
              p1(1,1+4) = p1(1,1)
              p1(2,1+4) = p1(2,1)
              p1(3,1+4) = DBLE(p1(1,1))*DTAN(DBLE(+0.5*dangle))

              p1(1,2)   = xorigin + xdelta * DBLE(ix-1)
              p1(2,2)   = yorigin + ydelta * DBLE(iy-1)
              p1(3,2)   = DBLE(p1(1,2))*DTAN(DBLE(-0.5*dangle))
              p1(1,2+4) = p1(1,2)
              p1(2,2+4) = p1(2,2)
              p1(3,2+4) = DBLE(p1(1,2))*DTAN(DBLE(+0.5*dangle))

              p1(1,3)   = xorigin + xdelta * DBLE(ix-1)
              p1(2,3)   = yorigin + ydelta * DBLE(iy)
              p1(3,3)   = DBLE(p1(1,3))*DTAN(DBLE(-0.5*dangle))
              p1(1,3+4) = p1(1,3)
              p1(2,3+4) = p1(2,3)
              p1(3,3+4) = DBLE(p1(1,3))*DTAN(DBLE(+0.5*dangle))

              p1(1,4)   = xorigin + xdelta * DBLE(ix)
              p1(2,4)   = yorigin + ydelta * DBLE(iy)
              p1(3,4)   = DBLE(p1(1,4))*DTAN(DBLE(-0.5*dangle))
              p1(1,4+4) = p1(1,4)
              p1(2,4+4) = p1(2,4)
              p1(3,4+4) = DBLE(p1(1,4))*DTAN(DBLE(+0.5*dangle))
c... 
              DO i1 = 1, 8
                obj(nobj)%v(1:3,i1) = p1(1:3,i1)
              ENDDO
c...          Rotate vertices:
              DO i1 = 1, 8
                x1 = obj(nobj)%v(1,i1)
                z1 = obj(nobj)%v(3,i1)
                obj(nobj)%v(1,i1) = DCOS(DBLE(ang)) * x1 -
     .                              DSIN(DBLE(ang)) * z1
                obj(nobj)%v(3,i1) = DSIN(DBLE(ang)) * x1 +
     .                              DCOS(DBLE(ang)) * z1
              ENDDO
c...          Surface 1:
              IF (isector.EQ.istart) THEN
                IF (iend-istart+1.EQ.nsector) THEN   
c...              Full torus:
                  obj(nobj)%tsur(1) = SP_GRID_SURFACE
                  obj(nobj)%nmap(1) = 1
                  obj(nobj)%imap(1,1) = nobj + ncells * (nsector - 1)
                  obj(nobj)%isur(1,1) = 6
                ELSE
                  obj(nobj)%tsur(1) = SP_GRID_BOUNDARY
                  obj(nobj)%nmap(1) = 1
                  obj(nobj)%imap(1,1) = nobj
                  obj(nobj)%isur(1,1) = 1
                ENDIF
              ELSE
                obj(nobj)%tsur(1) = SP_GRID_SURFACE
                obj(nobj)%nmap(1) = 1
                obj(nobj)%imap(1,1) = nobj - ncells
                obj(nobj)%isur(1,1) = 6
              ENDIF
              obj(nobj)%npts(1) = 4
              obj(nobj)%ipts(1,1) = 1
              obj(nobj)%ipts(2,1) = 2
              obj(nobj)%ipts(3,1) = 3
              obj(nobj)%ipts(4,1) = 4
c...          Surface 2:
              IF (iy.EQ.1) THEN     ! *** THIS IS NOT TRUE IF TARGET TRANSPARENT!
                obj(nobj)%tsur(2) = SP_GRID_BOUNDARY
                obj(nobj)%nmap(2) = 1
                obj(nobj)%imap(1,2) = nobj
                obj(nobj)%isur(1,2) = 2
              ELSE
                obj(nobj)%tsur(2) = SP_GRID_SURFACE
                obj(nobj)%nmap(2) = 1
                obj(nobj)%imap(1,2) = nobj + nxbin
                obj(nobj)%isur(1,2) = 4
              ENDIF
              obj(nobj)%npts(2) = 4
              obj(nobj)%ipts(1,2) = 2   ! If this is right, then the rotation is backwards! 
              obj(nobj)%ipts(2,2) = 1
              obj(nobj)%ipts(3,2) = 5
              obj(nobj)%ipts(4,2) = 6
c...          Surface 3:
              IF (ix.EQ.1) THEN
                obj(nobj)%tsur(3) = SP_GRID_BOUNDARY
                obj(nobj)%nmap(3) = 1
                obj(nobj)%imap(1,3) = nobj
                obj(nobj)%isur(1,3) = 3
              ELSE
                obj(nobj)%tsur(3) = SP_GRID_SURFACE
                obj(nobj)%nmap(3) = 1
                obj(nobj)%imap(1,3) = nobj - 1
                obj(nobj)%isur(1,3) = 5
              ENDIF
              obj(nobj)%npts(3) = 4
              obj(nobj)%ipts(1,3) = 3
              obj(nobj)%ipts(2,3) = 2
              obj(nobj)%ipts(3,3) = 6
              obj(nobj)%ipts(4,3) = 7
c...          Surface 4:
              IF (iy.EQ.nybin) THEN   
                obj(nobj)%tsur(4) = SP_GRID_BOUNDARY
                obj(nobj)%nmap(4) = 1
                obj(nobj)%imap(1,4) = nobj
                obj(nobj)%isur(1,4) = 4
              ELSE
                obj(nobj)%tsur(4) = SP_GRID_SURFACE
                obj(nobj)%nmap(4) = 1
                obj(nobj)%imap(1,4) = nobj - nxbin
                obj(nobj)%isur(1,4) = 2
              ENDIF
              obj(nobj)%npts(4) = 4
              obj(nobj)%ipts(1,4) = 4    ! If this is right, then the rotation is backwards! 
              obj(nobj)%ipts(2,4) = 3
              obj(nobj)%ipts(3,4) = 7
              obj(nobj)%ipts(4,4) = 8
c...          Surface 5:
              IF (ix.EQ.nxbin) THEN
                obj(nobj)%tsur(5) = SP_GRID_BOUNDARY
                obj(nobj)%nmap(5) = 1
                obj(nobj)%imap(1,5) = nobj
                obj(nobj)%isur(1,5) = 5
              ELSE
                obj(nobj)%tsur(5) = SP_GRID_SURFACE
                obj(nobj)%nmap(5) = 1
                obj(nobj)%imap(1,5) = nobj + 1
                obj(nobj)%isur(1,5) = 3
              ENDIF
              obj(nobj)%npts(5) = 4
              obj(nobj)%ipts(1,5) = 1
              obj(nobj)%ipts(2,5) = 4
              obj(nobj)%ipts(3,5) = 8
              obj(nobj)%ipts(4,5) = 5
c...          Surface 6:
              IF (isector.EQ.iend) THEN
                IF (iend-istart+1.EQ.nsector) THEN   
c...              Full torus:
                  obj(nobj)%tsur(6) = SP_GRID_SURFACE
                  obj(nobj)%nmap(6) = 1
                  obj(nobj)%imap(1,6) = nobj - ncells * (nsector - 1)
                  obj(nobj)%isur(1,6) = 1
                ELSE
                  obj(nobj)%tsur(6) = SP_GRID_BOUNDARY
                  obj(nobj)%nmap(6) = 1
                  obj(nobj)%imap(1,6) = nobj
                  obj(nobj)%isur(1,6) = 6
                ENDIF
              ELSE
                obj(nobj)%tsur(6) = SP_GRID_SURFACE
                obj(nobj)%nmap(6) = 1
                obj(nobj)%imap(1,6) = nobj + ncells
                obj(nobj)%isur(1,6) = 1
              ENDIF
              obj(nobj)%npts(6) = 4
              obj(nobj)%ipts(1,6) = 8
              obj(nobj)%ipts(2,6) = 7
              obj(nobj)%ipts(3,6) = 6
              obj(nobj)%ipts(4,6) = 5
            ENDDO
          ENDDO
 
        ENDDO

      ELSE
        CALL ER('BuildInversionMesh','Unknown option',*99)
      ENDIF

c...  Local mesh refinement adjustment based on mask, from a previous 
c     iteration:


c...  Generalized cell mapping method:

c      WRITE(0,*) '???',obj(8)%gsur(1)
      WRITE(6,*) 'DESPERATE'
      DO iobj = 1, nobj
        DO i1 = 1, 4
          i2 = i1 + 1
          IF (i2.GT.4) i2 = 1
          WRITE(6,*) iobj,i1,obj(iobj)%v(1:2,i1)
          WRITE(6,*) iobj,i1,obj(iobj)%v(1:2,i2)
        ENDDO
      ENDDO

 98   RETURN
 99   STOP
      END

