c     -*-Fortran-*-
c
c ======================================================================
c
c
c
c
c
c
c
      SUBROUTINE User_SurfaceReflectivity(imodel,scale)
      USE mod_out985
      IMPLICIT none

      TYPE(type_options985) :: opt
      INTEGER imodel
      REAL*8  scale

      SELECTCASE (opt%ref_wlgth(imodel))
        CASE (-1)

        CASE DEFAULT
          CALL ER('User_SurfaceReflectivity','Unknown option',*99)
      ENDSELECT

      RETURN
 99   STOP
      END
c
c ======================================================================
c



c
c ======================================================================
c
c
c
c
c
c
c
      SUBROUTINE User_LoadVesselGeometry(opt,idum2,dummy)
      USE mod_out985
      IMPLICIT none

      TYPE(type_options985) :: opt
      INTEGER idum2
      CHARACTER dummy*(*)

      CHARACTER cdum1*1024

      SELECTCASE (idum2)
        CASE (-1)
c...      Hard coded vessel geometry:
          READ(dummy,*) cdum1,
     .      opt%obj_type  (opt%obj_num),
     .      opt%obj_option(opt%obj_num),
     .      opt%obj_colour(opt%obj_num),
     .      opt%obj_reflec(opt%obj_num),
     .      opt%obj_fudge (opt%obj_num),
     .      opt%obj_factor(opt%obj_num),
     .      opt%obj_n     (opt%obj_num,1),
     .      opt%obj_n     (opt%obj_num,2)

        CASE DEFAULT
          WRITE(0,*) 'OFFENDING OPTION:',idum2
          CALL ER('User_LoadVesselGeometry','Unknown option',*99)
      ENDSELECT

      RETURN
 99   STOP
      END
c
c ======================================================================
c


c
c ======================================================================
c
c
c
c
c
c
c
      SUBROUTINE User_CustomObjects(ielement)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

c...  Input:
      INTEGER ielement

      INTEGER AddVertex,AddSurface

      TYPE(type_surface) newsrf
      INTEGER i1,idum1,istart
      REAL*8  newvtx(3,25),mat(3,3),angle,frac ,tmpvtx(3,25)

      newvtx = 0.0D0
      istart = nsrf + 1


      SELECTCASE (opt%obj_option(ielement))
        CASE (-1)
c...      Wall and floor:

          newvtx(1,1) =  0.33D0
          newvtx(2,1) =  2.19D0
          newvtx(1,2) =  2.00D0
          newvtx(2,2) =  2.19D0
          newvtx(1,3) =  2.00D0
          newvtx(2,3) = -2.19D0
          newvtx(1,4) =  0.33D0
          newvtx(2,4) = -2.19D0

          DO i1 = 4, 2, -1   ! ** REVERSED 'CAUSE PLOT IS BACKWARDS, BELOW ALSO... ***
            newsrf%type = SP_LINE_SEGMENT
            newsrf%nvtx = 2
            newsrf%ivtx(1) = AddVertex(newvtx(1,i1  ))
            newsrf%ivtx(2) = AddVertex(newvtx(1,i1-1))
            idum1 = AddSurface(newsrf)
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
          obj(nobj)%gsur(1)     = GT_TC
          obj(nobj)%tsur(1)     = SP_VESSEL_WALL
          obj(nobj)%reflec(1)   = opt%obj_reflec(ielement)
c..       Defunct:
          obj(nobj)%nsur        = 0
          obj(nobj)%ipts(2,1)   = 0
          obj(nobj)%nmap(1)     = 0

        CASE (-2)
c...      Center column:

          newvtx(1,1 ) =  0.330D0
          newvtx(2,1 ) = -2.190D0
          newvtx(1,2 ) =  0.330D0
          newvtx(2,2 ) = -2.150D0
          newvtx(1,3 ) =  0.190D0
          newvtx(2,3 ) = -1.905D0
          newvtx(1,4 ) =  0.190D0
          newvtx(2,4 ) = -1.820D0
          newvtx(1,5 ) =  0.175D0
          newvtx(2,5 ) = -1.720D0
          newvtx(1,6 ) =  0.175D0
          newvtx(2,6 ) = -1.670D0

          newvtx(1,7 ) =  0.280D0  ! 0.280D0
          newvtx(2,7 ) = -1.670D0
          newvtx(1,8 ) =  0.280D0  ! 0.280D0
          newvtx(2,8 ) = -1.450D0
          newvtx(1,9 ) =  0.280D0  ! 0.280D0
          newvtx(2,9 ) = -1.220D0

          newvtx(1,10) =  0.195D0
          newvtx(2,10) = -1.080D0
c...      Mirror:
          DO i1 = 1, 10    
            newvtx(1,i1+10) =  newvtx(1,11-i1)
            newvtx(2,i1+10) = -newvtx(2,11-i1)        
          ENDDO

          DO i1 = 20, 2, -1
            newsrf%type = SP_LINE_SEGMENT
            newsrf%nvtx = 2
            newsrf%ivtx(1) = AddVertex(newvtx(1,i1  ))
            newsrf%ivtx(2) = AddVertex(newvtx(1,i1-1))
            idum1 = AddSurface(newsrf)
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
          obj(nobj)%gsur(1)     = GT_TC
          obj(nobj)%tsur(1)     = SP_VESSEL_WALL
          obj(nobj)%reflec(1)   = opt%obj_reflec(ielement)
c..       Defunct:
          obj(nobj)%nsur        = 0
          obj(nobj)%ipts(2,1)   = 0
          obj(nobj)%nmap(1)     = 0

        CASE (-3)
c...      Divertor targets:


c   0.00400000     0.757000     -1.83102
c    0.0170000      1.69800     -1.83775
 
c    0.0900000     0.753000     -1.82500
c     0.199000      1.68800     -1.82500
c 
c   0.00400000     0.757000     -1.86500
c    0.0170000      1.69800     -1.86500
c 
c    0.0900000     0.753000     -1.86500
c     0.199000      1.68800     -1.86500


c...      Rotate 45 times:
c          DO angle = 0.0D0, 359D0, 7500D0
          DO angle = 0.0D0, 359D0, 7.5D0

c            newvtx(1,1) =  0.75700  ! LEON
c            newvtx(2,1) = -1.86500 
c            newvtx(3,1) =  0.00400
c            newvtx(1,2) =  1.69800
c            newvtx(2,2) = -1.86500
c            newvtx(3,2) =  0.01700 
c            newvtx(1,3) =  1.68800
c            newvtx(2,3) = -1.83102
c            newvtx(3,3) =  0.19900
c            newvtx(1,4) =  0.75300
c            newvtx(2,4) = -1.83775
c            newvtx(3,4) =  0.09000


            newvtx(1,4) =  0.73350   ! DO
            newvtx(2,4) = -1.83102
            newvtx(3,4) = -0.19100

            newvtx(1,3) =  0.75090
            newvtx(2,3) = -1.82500
            newvtx(3,3) = -0.10770

            newvtx(1,2) =  1.64220  
            newvtx(2,2) = -1.82500 
            newvtx(3,2) = -0.23300

            newvtx(1,1) =  1.60480
            newvtx(2,1) = -1.83780
            newvtx(3,1) = -0.41210 

c...        Rotate about y-axis (swing):
            CALL Calc_Transform2(mat,0.0D0,1,0)
            CALL Calc_Transform2(mat,angle*3.141592/180.0,2,1)
            DO i1 = 1, 4
              CALL Transform_Vect(mat,newvtx(1,i1))
            ENDDO

c            IF (angle.LT.10.0D0) THEN
c              WRITE(0,*)
c              WRITE(0,*) 'DATA:',newvtx(1:3,4)
c              WRITE(0,*) '    :',newvtx(1:3,3)
c              WRITE(0,*) '    :',newvtx(1:3,2)
c              WRITE(0,*) '    :',newvtx(1:3,1)
c            ENDIF

            newsrf%type = SP_PLANAR_POLYGON
            newsrf%nvtx = 4
            DO i1 = 4, 1, -1   ! ** REVERSED 'CAUSE PLOT IS BACKWARDS, BELOW ALSO... ***
              newsrf%ivtx(i1) = AddVertex(newvtx(1,i1))
            ENDDO
            idum1 = AddSurface(newsrf)


c...        Mirror:
            tmpvtx = newvtx
            newvtx(1:3,4) =  tmpvtx(1:3,3)
            newvtx(1:3,3) =  tmpvtx(1:3,4)
            newvtx(1:3,2) =  tmpvtx(1:3,1)
            newvtx(1:3,1) =  tmpvtx(1:3,2)
            newvtx(2,1:4) = -tmpvtx(2,1:4)
            newsrf%type = SP_PLANAR_POLYGON
            newsrf%nvtx = 4
            DO i1 = 4, 1, -1   ! ** REVERSED 'CAUSE PLOT IS BACKWARDS, BELOW ALSO... ***
              newsrf%ivtx(i1) = AddVertex(newvtx(1,i1))
            ENDDO
            idum1 = AddSurface(newsrf)

            IF (.TRUE.) THEN

              newvtx(1,1) =  1.64220  
              newvtx(2,1) = -1.82500 
              newvtx(3,1) = -0.23300
              newvtx(1,2) =  1.60480
              newvtx(2,2) = -1.83780
              newvtx(3,2) = -0.41210 
              newvtx(1,3) =  1.60480
              newvtx(2,3) = -1.90000
              newvtx(3,3) = -0.41210 
              newvtx(1,4) =  1.64220  
              newvtx(2,4) = -1.90000 
              newvtx(3,4) = -0.23300
c...          Rotate about y-axis (swing):
              DO i1 = 1, 4
                CALL Transform_Vect(mat,newvtx(1,i1))
              ENDDO
              newsrf%type = SP_PLANAR_POLYGON
              newsrf%nvtx = 4
              DO i1 = 4, 1, -1   ! ** REVERSED 'CAUSE PLOT IS BACKWARDS, BELOW ALSO... ***
                newsrf%ivtx(i1) = AddVertex(newvtx(1,i1))
              ENDDO
              idum1 = AddSurface(newsrf)
         
c...          Mirror:
              tmpvtx = newvtx
              newvtx(1:3,4) =  tmpvtx(1:3,3)
              newvtx(1:3,3) =  tmpvtx(1:3,4)
              newvtx(1:3,2) =  tmpvtx(1:3,1)
              newvtx(1:3,1) =  tmpvtx(1:3,2)
              newvtx(2,1:4) = -tmpvtx(2,1:4)
              newsrf%type = SP_PLANAR_POLYGON
              newsrf%nvtx = 4
              DO i1 = 4, 1, -1   ! ** REVERSED 'CAUSE PLOT IS BACKWARDS, BELOW ALSO... ***
                newsrf%ivtx(i1) = AddVertex(newvtx(1,i1))
              ENDDO
              idum1 = AddSurface(newsrf)
         
         
              newvtx(1,1) =  1.60480
              newvtx(2,1) = -1.83780
              newvtx(3,1) = -0.41210 
              newvtx(1,2) =  0.73350   ! DO
              newvtx(2,2) = -1.83102
              newvtx(3,2) = -0.19100
              newvtx(1,3) =  0.73350   
              newvtx(2,3) = -1.90000
              newvtx(3,3) = -0.19100
              newvtx(1,4) =  1.60480
              newvtx(2,4) = -1.90000
              newvtx(3,4) = -0.41210 
c...          Rotate about y-axis (swing):
              DO i1 = 1, 4
                CALL Transform_Vect(mat,newvtx(1,i1))
              ENDDO
              newsrf%type = SP_PLANAR_POLYGON
              newsrf%nvtx = 4
              DO i1 = 4, 1, -1   ! ** REVERSED 'CAUSE PLOT IS BACKWARDS, BELOW ALSO... ***
                newsrf%ivtx(i1) = AddVertex(newvtx(1,i1))
              ENDDO
              idum1 = AddSurface(newsrf)
c...          Mirror:
              newvtx(2,1) = -1.82500
              newvtx(2,2) = -1.82500 
              newvtx(2,1:4) = -newvtx(2,1:4)
              newsrf%type = SP_PLANAR_POLYGON
              newsrf%nvtx = 4
              DO i1 = 4, 1, -1   
                newsrf%ivtx(5-i1) = AddVertex(newvtx(1,i1))
              ENDDO
              idum1 = AddSurface(newsrf)


              newvtx(1,1) =  0.75090
              newvtx(2,1) = -1.82500
              newvtx(3,1) = -0.10770
              newvtx(1,2) =  1.64220  
              newvtx(2,2) = -1.82500 
              newvtx(3,2) = -0.23300
              newvtx(1,3) =  1.64220  
              newvtx(2,3) = -1.90000
              newvtx(3,3) = -0.23300
              newvtx(1,4) =  0.75090
              newvtx(2,4) = -1.90000
              newvtx(3,4) = -0.10770         
c...          Rotate about y-axis (swing):
              DO i1 = 1, 4
                CALL Transform_Vect(mat,newvtx(1,i1))
              ENDDO
              newsrf%type = SP_PLANAR_POLYGON
              newsrf%nvtx = 4
              DO i1 = 4, 1, -1   ! ** REVERSED 'CAUSE PLOT IS BACKWARDS, BELOW ALSO... ***
                newsrf%ivtx(i1) = AddVertex(newvtx(1,i1))
              ENDDO
              idum1 = AddSurface(newsrf)
c...          Mirror:
              newvtx(2,1) = -1.83102
              newvtx(2,2) = -1.83780 
              newvtx(2,1:4) = -newvtx(2,1:4)
              newsrf%type = SP_PLANAR_POLYGON
              newsrf%nvtx = 4
              DO i1 = 4, 1, -1 
                newsrf%ivtx(5-i1) = AddVertex(newvtx(1,i1))
              ENDDO
              idum1 = AddSurface(newsrf)


              newvtx(1,1) =  0.73350   ! DO
              newvtx(2,1) = -1.83102
              newvtx(3,1) = -0.19100
              newvtx(1,2) =  0.75090
              newvtx(2,2) = -1.82500
              newvtx(3,2) = -0.10770
              newvtx(1,3) =  0.75090
              newvtx(2,3) = -1.90000
              newvtx(3,3) = -0.10770
              newvtx(1,4) =  0.73350
              newvtx(2,4) = -1.90000
              newvtx(3,4) = -0.19100
c...          Rotate about y-axis (swing):
              DO i1 = 1, 4
                CALL Transform_Vect(mat,newvtx(1,i1))
              ENDDO
              newsrf%type = SP_PLANAR_POLYGON
              newsrf%nvtx = 4
              DO i1 = 4, 1, -1 
                newsrf%ivtx(i1) = AddVertex(newvtx(1,i1))
              ENDDO
              idum1 = AddSurface(newsrf)
c...          Mirror:
              tmpvtx = newvtx
              newvtx(1:3,4) =  tmpvtx(1:3,3)
              newvtx(1:3,3) =  tmpvtx(1:3,4)
              newvtx(1:3,2) =  tmpvtx(1:3,1)
              newvtx(1:3,1) =  tmpvtx(1:3,2)
              newvtx(2,1:4) = -tmpvtx(2,1:4)
              newsrf%type = SP_PLANAR_POLYGON
              newsrf%nvtx = 4
              DO i1 = 4, 1, -1 
                newsrf%ivtx(i1) = AddVertex(newvtx(1,i1))
              ENDDO
              idum1 = AddSurface(newsrf)
         
            ENDIF
         
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

        CASE (-4)
c...      MID flat plate only

          newvtx(1,1) =  0.742D0
          newvtx(2,1) = -1.950D0
          newvtx(1,2) =  0.742D0
          newvtx(2,2) = -1.838D0  ! -1.830D0
          newvtx(1,3) =  1.623D0
          newvtx(2,3) = -1.838D0  ! -1.830D0
          newvtx(1,4) =  1.623D0
          newvtx(2,4) = -1.950D0
c...      Mirror:
          DO i1 = 1, 4
            newvtx(1,9-i1) =  newvtx(1,i1)
            newvtx(2,9-i1) = -newvtx(2,i1)        
          ENDDO
          DO i1 = 1, 3
            newsrf%type = SP_LINE_SEGMENT
            newsrf%nvtx = 2
            newsrf%ivtx(1) = AddVertex(newvtx(1,i1  ))
            newsrf%ivtx(2) = AddVertex(newvtx(1,i1+1))
            idum1 = AddSurface(newsrf)
          ENDDO
          DO i1 = 5, 7
            newsrf%type = SP_LINE_SEGMENT
            newsrf%nvtx = 2
            newsrf%ivtx(1) = AddVertex(newvtx(1,i1  ))
            newsrf%ivtx(2) = AddVertex(newvtx(1,i1+1))
            idum1 = AddSurface(newsrf)
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
          obj(nobj)%gsur(1)     = GT_TC
          obj(nobj)%tsur(1)     = SP_VESSEL_WALL
          obj(nobj)%reflec(1)   = opt%obj_reflec(ielement)
c..       Defunct:
          obj(nobj)%nsur        = 0
          obj(nobj)%ipts(2,1)   = 0
          obj(nobj)%nmap(1)     = 0

        CASE (-5)
c...      P5:

          newvtx(1,1 ) =  1.560D0
          newvtx(2,1 ) =  0.580D0
          newvtx(1,2 ) =  1.570D0
          newvtx(2,2 ) =  0.590D0
          newvtx(1,3 ) =  1.735D0
          newvtx(2,3 ) =  0.590D0
          newvtx(1,4 ) =  1.745D0
          newvtx(2,4 ) =  0.580D0
          newvtx(1,5 ) =  1.745D0
          newvtx(2,5 ) =  0.415D0
          newvtx(1,6 ) =  1.735D0
          newvtx(2,6 ) =  0.405D0
          newvtx(1,7 ) =  1.570D0
          newvtx(2,7 ) =  0.405D0
          newvtx(1,8 ) =  1.560D0
          newvtx(2,8 ) =  0.415D0
          newvtx(1,9 ) =  1.560D0
          newvtx(2,9 ) =  0.580D0
c...      Mirror:
c          DO i1 = 1, 9
c            newvtx(1,i1+9) =  newvtx(1,i1)
c            newvtx(2,i1+9) = -newvtx(2,i1)        
c          ENDDO

          DO i1 = 1, 8
c          DO i1 = 9, 2, -1
            newsrf%type = SP_LINE_SEGMENT
            newsrf%nvtx = 2
            newsrf%ivtx(1) = AddVertex(newvtx(1,i1  ))
            newsrf%ivtx(2) = AddVertex(newvtx(1,i1+1))
            idum1 = AddSurface(newsrf)
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
          obj(nobj)%gsur(1)     = GT_TC
          obj(nobj)%tsur(1)     = SP_VESSEL_WALL
          obj(nobj)%reflec(1)   = opt%obj_reflec(ielement)
c..       Defunct:
          obj(nobj)%nsur        = 0
          obj(nobj)%ipts(2,1)   = 0
          obj(nobj)%nmap(1)     = 0

        CASE (-6)
c...      Gaps between divertor targets:


c...      Rotate 45 times:
          frac = 0.15

          DO angle = 0.0D0, 359D0, 7.5D0

            newvtx(1,4) =  0.75090
            newvtx(2,4) = -1.82500
            newvtx(3,4) = -0.10770

            newvtx(1,3) = (1.0 - frac) *   0.752155309817899 + 
     .                           frac  *   0.758533587564741
            newvtx(2,3) = (1.0 - frac) * (-1.83101999759674) + 
     .                           frac  * (-1.82500004768372) 
            newvtx(3,3) = (1.0 - frac) * (-9.362502535266594E-002) + 
     .                           frac  * (-8.766515263936292E-003)

            newvtx(1,2) = (1.0 - frac) *   1.64486053648259 + 
     .                           frac  *   1.65856334628506
            newvtx(2,2) = (1.0 - frac) * (-1.83780002593994) + 
     .                           frac  * (-1.82500004768372)  
            newvtx(3,2) = (1.0 - frac) * (-0.199106026346142) + 
     .                           frac  * (-1.665657951072375E-002)

            newvtx(1,1) =  1.64220  
            newvtx(2,1) = -1.82500 
            newvtx(3,1) = -0.23300

cvtx:  0.752155309817899       -1.83101999759674      -9.362502535266594E-002
cvtx:  0.758533587564741       -1.82500004768372      -8.766515263936292E-003

cvtx:   1.64486053648259       -1.83780002593994      -0.199106026346142
cvtx:   1.65856334628506       -1.82500004768372      -1.665657951072375E-002

c            newvtx(1,3) =  0.75215531
c            newvtx(2,3) = -1.83102
c            newvtx(3,3) = -0.0936250

c            newvtx(1,2) =  1.64486054
c            newvtx(2,2) = -1.83780
c            newvtx(3,2) = -0.19910603

c            newvtx(1,1) =  1.64220  
c            newvtx(2,1) = -1.82500 
c            newvtx(3,1) = -0.23300

c...        Rotate about y-axis (swing):
            CALL Calc_Transform2(mat,0.0D0,1,0)
            CALL Calc_Transform2(mat,angle*3.141592/180.0,2,1)
            DO i1 = 1, 4
              CALL Transform_Vect(mat,newvtx(1,i1))
            ENDDO

            IF (angle.EQ.7.5D0) THEN
              WRITE(0,*) 'vtx:',newvtx(1:3,3)
              WRITE(0,*) 'vtx:',newvtx(1:3,2)
            ENDIF

            newsrf%type = SP_PLANAR_POLYGON
            newsrf%nvtx = 4
            DO i1 = 4, 1, -1   ! ** REVERSED 'CAUSE PLOT IS BACKWARDS, BELOW ALSO... ***
              newsrf%ivtx(i1) = AddVertex(newvtx(1,i1))
            ENDDO
            idum1 = AddSurface(newsrf)

c...        Mirror:
            newvtx(1,4) = (1.0 - frac) *   0.75090  + frac *   0.73350
            newvtx(2,4) = (1.0 - frac) *   1.83102  + frac *   1.82500
            newvtx(3,4) = (1.0 - frac) * (-0.10770) + frac * (-0.19100)

            newvtx(1,3) =  0.75215530981790
            newvtx(2,3) =  1.82500004768372 
            newvtx(3,3) = -9.36250253526659E-002

            newvtx(1,2) =  1.64486053648259 
            newvtx(2,2) =  1.82500004768372  
            newvtx(3,2) = -0.19910602634614 

            newvtx(1,1) = (1.0 - frac) *   1.64220  + frac *   1.60480
            newvtx(2,1) = (1.0 - frac) *   1.83780  + frac *   1.82500
            newvtx(3,1) = (1.0 - frac) * (-0.23300) + frac * (-0.41210)

            WRITE(0,*) 'NEWVTX:',newvtx(1:3,4)

c...        Rotate about y-axis (swing):
            CALL Calc_Transform2(mat,0.0D0,1,0)
            CALL Calc_Transform2(mat,angle*3.141592/180.0,2,1)
            DO i1 = 1, 4
              CALL Transform_Vect(mat,newvtx(1,i1))
            ENDDO
            newsrf%type = SP_PLANAR_POLYGON
            newsrf%nvtx = 4
            DO i1 = 4, 1, -1   ! ** REVERSED 'CAUSE PLOT IS BACKWARDS, BELOW ALSO... ***
              newsrf%ivtx(5-i1) = AddVertex(newvtx(1,i1))
            ENDDO
            idum1 = AddSurface(newsrf)

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

        CASE (-7)
c...      Big screen, to catch rays from incomplete 3D grids:

          newvtx(1,1) =  0.000D0
          newvtx(2,1) =  3.000D0
          newvtx(1,2) =  3.000D0
          newvtx(2,2) =  3.000D0
          newvtx(1,3) =  3.000D0
          newvtx(2,3) = -3.000D0
          newvtx(1,4) =  0.000D0
          newvtx(2,4) = -3.000D0

          DO i1 = 4, 2, -1
            newsrf%type = SP_LINE_SEGMENT
            newsrf%nvtx = 2
            newsrf%ivtx(1) = AddVertex(newvtx(1,i1  ))
            newsrf%ivtx(2) = AddVertex(newvtx(1,i1-1))
            idum1 = AddSurface(newsrf)
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
          obj(nobj)%gsur(1)     = GT_TC
          obj(nobj)%tsur(1)     = SP_VESSEL_WALL
          obj(nobj)%reflec(1)   = opt%obj_reflec(ielement)
c..       Defunct:
          obj(nobj)%nsur        = 0
          obj(nobj)%ipts(2,1)   = 0
          obj(nobj)%nmap(1)     = 0

        CASE DEFAULT
          WRITE(0,*) 'OFFENDING OPTION:',opt%obj_type(ielement)
          CALL ER('User_CustomObjects','Unknown option',*99)
      ENDSELECT

      RETURN
 99   STOP
      END
c
c ======================================================================
c