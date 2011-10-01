c     -*-Fortran-*-
c
c ======================================================================
c
c function: AddSurface
c
c
      INTEGER FUNCTION AddSurface(newsrf)
      USE MOD_OUT985
      IMPLICIT none

      TYPE(type_surface) newsrf

      IF (.NOT.ALLOCATED(srf)) THEN
c...    Check that array allocated.  This must be allocated/initialized
c       before the call to this array to be sure that NSRF=0:
        CALL ER('AddSurface','SRF array not allocated',*99)
      ELSEIF (nsrf+1.GE.maxsrf) THEN
c...    The size of SRF is not sufficient, so make it bigger:
        CALL ALLOC_SURFACE(-1,MP_INCREASE_SIZE)
      ENDIF

c...  Scan over existing vertices to see if it already exists:

c...  Add the vertex to the list:
      nsrf = nsrf + 1
      srf(nsrf) = newsrf
      AddSurface = nsrf

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c function: AddVertex
c
c
      INTEGER FUNCTION AddVertex(newvtx)
      USE MOD_OUT985
      IMPLICIT none

      REAL*8 newvtx(3)

      INTEGER, PARAMETER :: VTXSTEP = 10000
      REAL*8 , PARAMETER :: DTOL    = 1.0D-07

      INTEGER ivtx
      LOGICAL message

      DATA message /.TRUE./
      SAVE

      IF (.NOT.ALLOCATED(vtx)) THEN
c...    Initial allocation of memory for VTX (vertex list) array:
        CALL Alloc_Vertex(-1,MP_INITIALIZE)
      ELSEIF (nvtx+1.GE.maxvtx) THEN
c...    The size of VTX is not sufficient, so make it bigger:
        CALL Alloc_Vertex(-1,MP_INCREASE_SIZE)
      ENDIF

c...  Scan over existing vertices to see if it already exists:
      IF (.FALSE.) THEN
        DO ivtx = 1, nvtx   
          IF (DABS(vtx(1,ivtx)-newvtx(1)).LT.DTOL.AND.
     .        DABS(vtx(2,ivtx)-newvtx(2)).LT.DTOL.AND.
     .        DABS(vtx(3,ivtx)-newvtx(3)).LT.DTOL) THEN
c...        Found one:
c            WRITE(0,*) 'ADDVERTEX: FOUND DUPLICATE!'
            AddVertex = ivtx
            RETURN
         ENDIF
        ENDDO
      ELSEIF (message) THEN
        WRITE(0,*) 
        WRITE(0,*) '--------------------------------------------------'
        WRITE(0,*) '     *** NOT CHECKING FOR DUPLICATE VTXS ***      '
        WRITE(0,*) '--------------------------------------------------'
        WRITE(0,*) 
        message = .FALSE.  
      ENDIF

c...  Add the vertex to the list:
      nvtx = nvtx + 1
      vtx(1:3,nvtx) = newvtx(1:3)
      AddVertex = nvtx

      RETURN
 99   STOP
      END
c
c ====================================================================== 
c
      SUBROUTINE BuildSeparatrix(time,nsep,xsep,ysep,mode)
      IMPLICIT none

      INTEGER nsep,mode
      REAL    time,xsep(nsep),ysep(nsep)  ! Can I send pointers for these?

      INTEGER   fp,i1,npts
      CHARACTER dummy*1024

      REAL, ALLOCATABLE :: r(:),z(:),psin(:)

c...  Load 2D magnetic equilibrium data (psin):      
      fp = 99
      OPEN(fp,FILE='equilibrium.dat',FORM='FORMATTED',
     .     STATUS='OLD',ERR=98) 
      DO WHILE (.TRUE.)
        READ(fp,'(A)',END=10,ERR=96) dummy
        IF     (dummy(1:6).EQ.'[HEAD]') THEN
        ELSEIF (dummy(1:6).EQ.'[SHOT]') THEN
        ELSEIF (dummy(1:6).EQ.'[NPTS]') THEN
          READ(dummy(7:),*) npts
          IF (.NOT.ALLOCATED(r))    ALLOCATE(r(npts))
          IF (.NOT.ALLOCATED(z))    ALLOCATE(z(npts))
          IF (.NOT.ALLOCATED(psin)) ALLOCATE(psin(npts))
        ELSEIF (dummy(1:6).EQ.'[TIME]') THEN
c...      Need option for interpolating between time slices:
        ELSEIF (dummy(1:6).EQ.'[DATA]') THEN
c...      Need to check if this is the correct time:
          DO i1 = 1, npts
            READ(fp,*,END=97) r(i1),z(i1),psin(i1)
          ENDDO
        ELSE
          CALL ER('BuildSeparatrix','Invalid data file format',*99)
        ENDIF
      ENDDO
 10   CONTINUE

c...  Divide the area where there is equilibrium data into a regular
c     grid and identify the ... 
      



c...  Clear memory:
      IF (ALLOCATED(r))    DEALLOCATE(r)
      IF (ALLOCATED(z))    DEALLOCATE(z)
      IF (ALLOCATED(psin)) DEALLOCATE(psin)

      RETURN
 96   CALL ER('BuildSeparatrix','Problem reading equilibrium file',*99)
 97   CALL ER('BuildSeparatrix','Unexpected end of file',*99)
 98   CALL ER('BuildSeparatrix','Unable to open equilibrium file',*99)
 99   STOP
      END
c
c ======================================================================
c
      subroutine calc_transform2(m,angle,axis,init)
      implicit none
c
      real*8 m(3,3)
      real*8 angle
      integer axis,init
c
c     CALC_TRANSFORM: calculate the transformation matrix resulting from a 
c                     rotation of magnitude angle about the given axis. 
c
c      
      real*8 trans(3,3)
      real*8 tempm(3,3)
      integer i,j
c
c     Initialize to identity matrix if initialization option on
c
      if (init.eq.0) then 
         m = 0.0D0
c         call dzero(m,9)
         m(1,1) = 1.0D0 
         m(2,2) = 1.0D0 
         m(3,3) = 1.0D0 
      endif
c
      trans = 0.0D0
      tempm = 0.0D0
c      call dzero(trans,9)
c      call dzero(tempm,9)
c
c     Calculate transformation matrix
c
c     Rotation about first-axis 
c
      if (axis.eq.1) then 
         trans(1,1) =  1.0
         trans(2,2) =  cos(angle)
         trans(3,3) =  cos(angle)
         trans(2,3) = -sin(angle) 
         trans(3,2) =  sin(angle) 
c
c     Rotation about second-axis 
c
      elseif (axis.eq.2) then 
         trans(1,1) =  cos(angle)
         trans(2,2) =  1.0
         trans(3,3) =  cos(angle)
         trans(1,3) = -sin(angle) 
         trans(3,1) =  sin(angle) 
c
c     Rotation about third-axis
c
      elseif (axis.eq.3) then 
         trans(1,1) =  cos(angle)
         trans(2,2) =  cos(angle)
         trans(3,3) =  1.0
         trans(1,2) = -sin(angle) 
         trans(2,1) =  sin(angle) 
      endif 
c
c      
c     Multiply transformation and matrix - pre-multiply 
c
      do i=1,3
         do j = 1,3
c                         
            tempm(i,j) =   trans(i,1) * m(1,j) 
     >                   + trans(i,2) * m(2,j) 
     >                   + trans(i,3) * m(3,j)
c
         end do 
      end do
c  
c     Assign resulting transformation to m
c      
c     m = tempm
c
      do i = 1,3
         do j = 1,3
            m(i,j) = tempm(i,j)
         end do
      end do 
c
c      write(0,'(a,i4,12(1x,g12.5))') 'M-CALC:',axis,
c     >                     ((m(i,j),j=1,3),i=1,3)
c
      return
      end
c
c
c ======================================================================
c
      SUBROUTINE PointInPoly(p3,iobj,iside,isur,plane,result,status)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

      INTEGER iobj,iside,isur,plane,status
      LOGICAL result
      REAL*8  p3(3)

      integer nextv,nv,v
      LOGICAL inpoly,output 
      real*8  cp,lastcp,x0,x1,x2,y0,y1,y2

c...  Project the point and polygon onto the closes plane: 

c      output = .TRUE.

c
      lastcp = 0.0D0
c
      inpoly = .false.

      IF (obj(iobj)%nside.NE.0) THEN
        nv = srf(isur)%nvtx
      ELSE
        nv = obj(iobj)%npts(iside)
      ENDIF

      DO v = 1, nv
         if (v.eq.nv) then
            nextv = 1
         else
            nextv = v+1
         endif
         IF     (plane.EQ.1) THEN
           x0 = p3(2) 
           y0 = p3(3) 
           IF (obj(iobj)%nside.NE.0) THEN
c...         Need to know which surface for the side: (rename isur to iside)
             x1 = vtx(2,srf(isur)%ivtx(v    ))
             y1 = vtx(3,srf(isur)%ivtx(v    ))
             x2 = vtx(2,srf(isur)%ivtx(nextv))
             y2 = vtx(3,srf(isur)%ivtx(nextv))
           ELSE
             x1 = obj(iobj)%v(2,obj(iobj)%ipts(v    ,iside))
             y1 = obj(iobj)%v(3,obj(iobj)%ipts(v    ,iside))
             x2 = obj(iobj)%v(2,obj(iobj)%ipts(nextv,iside))
             y2 = obj(iobj)%v(3,obj(iobj)%ipts(nextv,iside))
           ENDIF
         ELSEIF (plane.EQ.2) THEN
           x0 = p3(1) 
           y0 = p3(3) 
           IF (obj(iobj)%nside.NE.0) THEN
c...         Need to know which surface for the side:
             x1 = vtx(1,srf(isur)%ivtx(v    ))
             y1 = vtx(3,srf(isur)%ivtx(v    ))
             x2 = vtx(1,srf(isur)%ivtx(nextv))
             y2 = vtx(3,srf(isur)%ivtx(nextv))
           ELSE
             x1 = obj(iobj)%v(1,obj(iobj)%ipts(v    ,iside))
             y1 = obj(iobj)%v(3,obj(iobj)%ipts(v    ,iside))
             x2 = obj(iobj)%v(1,obj(iobj)%ipts(nextv,iside))
             y2 = obj(iobj)%v(3,obj(iobj)%ipts(nextv,iside)) 
           ENDIF
         ELSEIF (plane.EQ.3) THEN
           x0 = p3(1) 
           y0 = p3(2) 
           IF (obj(iobj)%nside.NE.0) THEN
             x1 = vtx(1,srf(isur)%ivtx(v    ))
             y1 = vtx(2,srf(isur)%ivtx(v    ))
             x2 = vtx(1,srf(isur)%ivtx(nextv))
             y2 = vtx(2,srf(isur)%ivtx(nextv))
           ELSE
             x1 = obj(iobj)%v(1,obj(iobj)%ipts(v    ,iside))
             y1 = obj(iobj)%v(2,obj(iobj)%ipts(v    ,iside))
             x2 = obj(iobj)%v(1,obj(iobj)%ipts(nextv,iside))
             y2 = obj(iobj)%v(2,obj(iobj)%ipts(nextv,iside)) 
           ENDIF
         ENDIF

         IF (status.EQ.1.AND..FALSE.) THEN
           WRITE(0,*) 'PLAN:',plane
           WRITE(0,*) 'VER0:',x0,y0
           WRITE(0,*) 'VER1:',x1,y1
           WRITE(0,*) 'VER2:',x2,y2

           WRITE(0,*) 'VECS1:',x0-x1,y2-y1
           WRITE(0,*) 'VECS2:',y0-y1,x2-x1
         ENDIF

         cp = ( (x0-x1) * (y2-y1) ) - ( (y0-y1) * (x2-x1) )
c
c
         if (v.eq.1.and.cp.ne.0.0D0) lastcp = cp

         IF (status.EQ.1) THEN
            WRITE(0,'(A,2I6,2F12.6)') 'CP:',v,nextv,cp,lastcp
          ENDIF
c         WRITE(6,'(A,3D18.6)') 'CP:',cp,lastcp,cp*lastcp
c
          if ((lastcp * cp).lt.0.0D0) GOTO 10
c
          if (cp.ne.0.0) lastcp = cp
c
      ENDDO

      inpoly = .true.

      IF (status.EQ.1) THEN
         WRITE(0,'(A)') 'GOOD!'
      ENDIF
c      WRITE(6,'(A)') 'GOOD!'

 10   CONTINUE


c      IF (obj(iobj)%nside.NE.0) inpoly = .FALSE.
 
      IF (inpoly) THEN
        result = .TRUE.
      ELSE
        result = .FALSE.
      ENDIF


      RETURN
 99   STOP
      END
c
c ====================================================================== 
c
      SUBROUTINE LineThroughSurface(v1,v2,iobj,iside,isur,n,v,d,status)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

      REAL*8 v1(3),v2(3)
      INTEGER n
      INTEGER iobj,iside,isur,mode,status
      REAL*8  v(3,MAXINTER),d(MAXINTER)

      REAL*8     DTOL 
      PARAMETER (DTOL=1.0D-10)

      INTEGER i1,is,plane,fp
      LOGICAL result,output
      REAL*8  av(3),bv(3),nv(3),p31(3),p21(3),p3(3),u,mind(5),nmind,
     .        denom,DTOL2,
     .        x1,x2,y1,y2,z1,z2,dx12,dy12,dz12,r3,r4,y3,y4,dr34,dy34,
     .        r,y,a,b,c,s(2),t,beta,gamma,b24ac

      fp = 6

      IF (nchord.EQ.-1) THEN
        output = .TRUE.
      ELSE
        output = .FALSE.
      ENDIF


c...  Provide some scaling that is sensitive to the length of the viewing chord:  Maybe the viewing chord shouldn't shrink?
c      DTOL2 = 1.0D-09 /    ! *** Hopefully this is strong enough...
c     .       ((v1(1)-v2(1))**2 + (v1(2)-v2(2))**2 + (v1(3)-v2(3))**2)
c      IF (nchord.EQ.dchord) WRITE(fp,*) '     DTOL2:',DTOL2
      DTOL2 = 0.0D0

c...  Find intersection of the chord with the plane


      IF     (obj(iobj)%gsur(iside).EQ.GT_TC) THEN

c
c       line in x,y,z space is converted to r,y space
c       surface is in r,y space already, with r=x for toroidally coninuous surfaces
c 
c       view is p1(x1,y1,z1),p2(x2,y2,z3)
c          dx12 is x2 - x1, etc.
c          equation of line is  x(s) = x1 + s * dx12
c
c       surface is p3(r3,y3),p4(r4,y4)
c          equation of line is  x(t) = x3 + t * dx34
c
c       line equation for r(s) in r,y system is
c          r(s) = SQRT ( (x1 + s * dx12)**2 + (z1 + s * dz12)**2 )
c
c          set r(s) = r3 + t * dr34 
c
c          and get rif of t from y = (
c
c          algebra gets
c
c          r(s) = beta * s + gamma
c
c


        IF (output) THEN
          WRITE(0,*) '   TC:'
          WRITE(0,*) '   TC is:',iobj,iside
          WRITE(0,*) '   TC  t:',obj(iobj)%type.EQ.OP_INTEGRATION_VOLUME
          WRITE(0,'(1X,A,3F10.4)') '   TC v2:',v1
          WRITE(0,'(1X,A,3F10.4)') '   TC v1:',v2
        ENDIF

        s = -999.0
c        v =  0.0D0
c        d = -1.0D0         
 
        x1 = v1(1)
        x2 = v2(1)
        y1 = v1(2)
        y2 = v2(2)
        z1 = v1(3)
        z2 = v2(3)
        dx12 = x2 - x1
        dy12 = y2 - y1
        dz12 = z2 - z1
        IF (obj(iobj)%nside.NE.0) THEN
          r3 = vtx(1,srf(isur)%ivtx(1))
          r4 = vtx(1,srf(isur)%ivtx(2))
          y3 = vtx(2,srf(isur)%ivtx(1))
          y4 = vtx(2,srf(isur)%ivtx(2))
        ELSE
          r3 = obj(iobj)%v(1,obj(iobj)%ipts(1,iside))
          r4 = obj(iobj)%v(1,obj(iobj)%ipts(2,iside))
          y3 = obj(iobj)%v(2,obj(iobj)%ipts(1,iside))
          y4 = obj(iobj)%v(2,obj(iobj)%ipts(2,iside))
        ENDIF
        dr34 = r4 - r3
        dy34 = y4 - y3

        IF (ABS(dy34).LT.DTOL) THEN
c...      Surface is horizontal:          

          IF (ABS(dy12).LT.DTOL) THEN
c...        Viewing chord is also horizontal:          
            IF (ABS(y1-y3).LT.DTOL) THEN
              WRITE(0,*) 'NCHORD:',nchord
              STOP 'THIS IS A BAD SITUATION' 
            ELSE
c...          No intersection because...
            ENDIF
          ELSE
            s(1) = (y3 - y1) / dy12
 
            IF (nchord.EQ.dchord) WRITE(fp,*) '   s1:',s

c                                           *** GONNA WORK? ***
            IF (s(1).GE.DTOL.AND.s(1).LT.1.0D0+DTOL2) THEN 

              r = SQRT((x1+s(1)*dx12)**2 + (z1+s(1)*dz12)**2)

              t = (r - r3) / dr34

              IF (nchord.EQ.dchord) WRITE(fp,*) '   t1:',t

              IF (t.GT.0.0D0.AND.t.LT.1.0D0) THEN
                n = n + 1
                v(1,n) = x1 + s(1) * dx12
                v(2,n) = y1 + s(1) * dy12
                v(3,n) = z1 + s(1) * dz12
                d(n) =  DSQRT((v1(1) - v(1,n))**2 +
     .                        (v1(2) - v(2,n))**2 +
     .                        (v1(3) - v(3,n))**2)
c                IF (d(n).LT.RTOL) THEN
c                  WRITE(0,*) 'JIGGY A',nchord
c                  n = n - 1
c                ENDIF
              ELSE
c...            No intersection because ...
              ENDIF            
            ELSE
c...          No intersection because ...
            ENDIF
          ENDIF

        ELSE
          beta  = dr34 * dy12 / dy34
          gamma = ((y1 - y3) * dr34 / dy34) + r3
        
          a = dx12**2 + dz12**2 - beta**2
          b = 2 * (x1 * dx12 + z1 * dz12 - beta * gamma)
          c = x1**2 + z1**2 - gamma**2

          b24ac = b**2 - 4 * a * c        

          IF (b24ac.LE.0.0D0) THEN
c...        No intersection because...
          ELSE
            s(2) = SQRT(b24ac)
            s(1) = (-b - s(2)) / (2.0D0 * a)
            s(2) = (-b + s(2)) / (2.0D0 * a)

            DO is = 1, 2

              IF (nchord.EQ.dchord) WRITE(fp,*) '   s2:',s(is)
c                                         *** GONNA SCREW THINGS UP? ***
              IF (s(is).GE.DTOL.AND.s(is).LE.1.0D0+DTOL2) THEN

                y = y1 + s(is) * dy12

                t = (y - y3) / dy34
 
                IF (nchord.EQ.dchord) WRITE(fp,*) '   t2:',t
                IF (t.GT.0.0D0.AND.t.LT.1.0D0) THEN
                  n = n + 1
                  v(1,n) = x1 + s(is) * dx12
                  v(2,n) = y
                  v(3,n) = z1 + s(is) * dz12
                  d(n) =  DSQRT((v1(1) - v(1,n))**2 +
     .                          (v1(2) - v(2,n))**2 +
     .                          (v1(3) - v(3,n))**2)
c                  IF (d(n).LT.RTOL) THEN
c                    WRITE(0,*) 'JIGGY B',nchord,d(n)
c                    n = n - 1
c                  ENDIF

                ELSE
c...              No intersection because ...
                ENDIF
              ELSE
c...            No intersection because ...
              ENDIF
            ENDDO

          ENDIF
        ENDIF

        IF (output) THEN
          WRITE(0,*) '   TC  s:',s(1),s(2)
          WRITE(0,*) '   TC  t:',t

          WRITE(0,'(1X,A,3F10.4)') '   TC  v:',v(1:3,1)
          WRITE(0,'(1X,A,3F10.4)') '   TC  v:',v(1:3,2)
          WRITE(0,*) '   TC  d:',d(1)
          WRITE(0,*) '   TC  d:',d(2)
        ENDIF

      ELSEIF (obj(iobj)%gsur(iside).EQ.GT_TD) THEN
c...    Floating 3D cartesian surface (most general toroidally 
c       discretized representation):

        DO i1 = 1, 3
c...      Store this?  Certainly...
          IF (obj(iobj)%nside.NE.0) THEN
            IF (srf(isur)%ivtx(1).EQ.0.OR.
     .          srf(isur)%ivtx(1).EQ.0) THEN
              WRITE(0,*) 'DATA:',iobj,iside,isur
              WRITE(0,*) 'DATA:',srf(isur)%ivtx(1)
              WRITE(0,*) 'DATA:',srf(isur)%ivtx(2)
              WRITE(0,*) 'DATA:',srf(isur)%ivtx(3)
              STOP 'PROBLEM BBB'
            ENDIF
            av(i1) = vtx(i1,srf(isur)%ivtx(1)) - 
     .               vtx(i1,srf(isur)%ivtx(2))
            bv(i1) = vtx(i1,srf(isur)%ivtx(3)) - 
     .               vtx(i1,srf(isur)%ivtx(2))
          ELSE
            IF (obj(iobj)%ipts(1,iside).EQ.0.OR.
     .          obj(iobj)%ipts(1,iside).EQ.0) THEN
              WRITE(0,*) 'DATA:',iobj,iside
              WRITE(0,*) 'DATA:',obj(iobj)%ipts(1,iside)
              WRITE(0,*) 'DATA:',obj(iobj)%ipts(2,iside)
              WRITE(0,*) 'DATA:',obj(iobj)%ipts(3,iside)
              STOP 'PROBLEM AAA'
            ENDIF
            av(i1) = obj(iobj)%v(i1,obj(iobj)%ipts(1,iside)) - 
     .               obj(iobj)%v(i1,obj(iobj)%ipts(2,iside))
            bv(i1) = obj(iobj)%v(i1,obj(iobj)%ipts(3,iside)) - 
     .               obj(iobj)%v(i1,obj(iobj)%ipts(2,iside))
          ENDIF
        ENDDO
c...    Surface normal:
        nv(1) =  av(2) * bv(3) - av(3) * bv(2)
        nv(2) = -av(1) * bv(3) + av(3) * bv(1)
        nv(3) =  av(1) * bv(2) - av(2) * bv(1) 
        IF (output.AND.iobj.EQ.38257) THEN
          WRITE(0,'(A,3F10.5)') 'MAGIC N:',nv(1),nv(2),nv(3)
        ENDIF
c...    Vectors:
        DO i1 = 1, 3
          IF (obj(iobj)%nside.NE.0) THEN
            p31(i1) = vtx(i1,srf(isur)%ivtx(2)) - v1(i1)
            p21(i1) = v2(i1)                    - v1(i1)
          ELSE
            p31(i1) = obj(iobj)%v(i1,obj(iobj)%ipts(2,iside)) - v1(i1)
            p21(i1) = v2(i1)                                  - v1(i1)
          ENDIF
        ENDDO

c...    Need checks in case chord is parallel with surface...
        denom = nv(1) * p21(1) + nv(2) * p21(2) + nv(3) * p21(3)
        IF (DABS(denom).GT.DTOL) THEN
          u = (nv(1) * p31(1) + nv(2) * p31(2) + nv(3) * p31(3)) / denom
        ELSE
          u = -1.0D0
        ENDIF

c...    Should also pass minu so that if u.GT.minu it is not even checked if the
c       intersection point is on the polygon:

c...    This will always be satisfied, unless the chord and surface normal are orthogonal?

        IF (status.EQ.1) WRITE(0,*) 'U:',u
c        WRITE(6,*) '  u:',u

        IF (u.GT.0.0D0.AND.u.LT.1.0D0+DTOL) THEN

c...      This point is on the plane:
          IF (output.AND.iobj.EQ.38257) THEN
            WRITE(0,*) 'MAGIC U:',u
          ENDIF
          DO i1 = 1, 3
            p3(i1) = v1(i1) + u * (v2(i1) - v1(i1))
          ENDDO
c...      Project onto the appropriate plane:
          IF     (ABS(nv(1)).GT.ABS(nv(2)).AND.
     .            ABS(nv(1)).GT.ABS(nv(3))) THEN
            plane = 1
          ELSEIF (ABS(nv(2)).GT.ABS(nv(1)).AND.
     .            ABS(nv(2)).GT.ABS(nv(3))) THEN
            plane = 2
          ELSEIF (ABS(nv(3)).GT.ABS(nv(1)).AND.
     .            ABS(nv(3)).GT.ABS(nv(2))) THEN
            plane = 3
          ELSEIF (ABS(nv(1)).EQ.ABS(nv(3))) THEN
            plane = 1
          ELSEIF (ABS(nv(1)).EQ.ABS(nv(2))) THEN
            plane = 2
          ELSEIF (ABS(nv(2)).EQ.ABS(nv(3))) THEN
            plane = 3
          ELSE
            CALL ER('CheckSurfaceIntesection','Bad projection',*99)
          ENDIF

          IF (output.AND.iobj.EQ.38257) THEN
            WRITE(0,'(A,3F10.5)') 'MAGIC P:',p3(1),p3(2),p3(3)
            WRITE(0,*) 'MAGIC PLANE:',plane
          ENDIF

c...      Check if the point is inside the surface polygon:
          CALL PointInPoly(p3,iobj,iside,isur,plane,result,status)

          IF (result) THEN
            n = n + 1
            v(1:3,n) = p3(1:3)
            d(n) =  DSQRT((v1(1) - v(1,n))**2 +
     .                    (v1(2) - v(2,n))**2 +
     .                    (v1(3) - v(3,n))**2)
          ELSE
c            v =  0.0D0
c            d = -1.0D0
          ENDIF

        ELSE
c          v =  0.0D0
c          d = -1.0D0
c         CALL ER('LineThroughSurface','This no happen',*99)
        ENDIF     
     
      ELSE
        CALL ER('LineThroughSurface','Unknown surface geometry '//
     .          'type',*99)
      ENDIF



      RETURN
 99   WRITE(0,*) 'NORMAL:',nv(1),nv(2),nv(3)
      WRITE(0,*) '   TC  s:',s(1),s(2)
      STOP
      END
c
c ====================================================================== 
c
      SUBROUTINE FindSurfaceIntersections(v1,v2,mode,status)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none
c...  Input:
      REAL*8  v1(3),v2(3)
      INTEGER mode,status

      INTEGER iobj,iside,isrf,i1,i2,i3,i4,isrf1,isrf2,fp,
     .        ninter,ointer(0:MAXINTER),sinter(0:MAXINTER),
     .        srfinter(0:MAXINTER),
     .        nsurlist,count,n
      LOGICAL cont     
      REAL*8  v(3,MAXINTER),d(MAXINTER),
     .        vinter(3,0:MAXINTER),dinter(0:MAXINTER)

      INTEGER, POINTER :: surlist(:,:)

      fp = 6

      ninter = 0
      dinter = 0.0D0


      IF     (mode.EQ.IT_VWINTER) THEN
        nsurlist = nvwlist
        surlist => vwlist
        IF (nchord.EQ.dchord) WRITE(fp,*) 'SEARCH: vessel wall'
      ELSEIF (mode.EQ.IT_GBINTER) THEN
        nsurlist = ngblist
        surlist => gblist
        IF (nchord.EQ.dchord) WRITE(fp,*) 'SEARCH: grid boundary'
      ELSEIF (mode.EQ.IT_OBINTER) THEN
        nsurlist = noblist
        surlist => oblist
        IF (nchord.EQ.dchord) WRITE(fp,*) 'SEARCH: object map'
      ELSE
        CALL ER('FindSurfaceIntersections','MODE problem',*99)
      ENDIF
      IF (nchord.EQ.dchord) THEN
        WRITE(fp,*) '  NSURLIST:',nsurlist
        IF (nsurlist.GT.0) THEN
          WRITE(fp,*) '  SURLIST1:',surlist(1:nsurlist,1)
          WRITE(fp,*) '  SURLIST2:',surlist(1:nsurlist,2)
        ENDIF
      ENDIF


      DO i3 = 1, nsurlist
        iobj  = surlist(i3,1)
        iside = surlist(i3,2)

        IF (nchord.EQ.dchord) THEN
          WRITE(fp,*) ' INTER?:',mode,i3,iobj,iside
        ENDIF

        IF (obj(iobj)%nside.NE.0) THEN
          IF (obj(iobj)%nside.GT.1) THEN
            isrf1 = obj(iobj)%iside(iside,1)    
            isrf2 = obj(iobj)%iside(iside,2)            
c            STOP 'CANNOT DO NSIDE.GT.1 IT SEEMS...'
          ELSE
            isrf1 = obj(iobj)%iside(1,1)    
            isrf2 = obj(iobj)%iside(1,2)
          ENDIF
        ELSE
          isrf1 = 1
          isrf2 = 1
        ENDIF

        DO isrf = isrf1, isrf2
          n = 0
          v =  0.0D0
          d = -1.0D0         
          CALL LineThroughSurface(v1,v2,iobj,iside,isrf,n,v,d,status) 
c          WRITE(fp,*) 'D:',d(1:n)
          DO i4 = 1, n
            IF (ninter.LT.MAXINTER) THEN  
              IF (obj(iobj)%gsur(iside).EQ.GT_TC.AND.  ! Don't include intersection where 
     .            d(i4).GT.1.0D-10.OR.                 ! the point is *on* a surface...
     .            obj(iobj)%gsur(iside).EQ.GT_TD.AND.  
     .            d(i4).GT.1.0D-20) THEN
c              IF (d(i4).GT.1.0D-10) THEN  
                ninter = ninter + 1        
                dinter(ninter) = d(i4)
                vinter(1:3,ninter) = v(1:3,i4)

                ointer(ninter) = iobj    ! Reform naming here...
                sinter(ninter) = iside
                srfinter(ninter) = isrf

                IF (nchord.EQ.dchord) THEN
                  WRITE(fp,*) '    HIT:',iobj,iside,isrf
                ENDIF

              ENDIF
            ELSE
              CALL ER('FindSurfaceIntersections','Zounds!',*99)
            ENDIF

          ENDDO
        ENDDO
      ENDDO

          IF (nchord.EQ.dchord) THEN
            WRITE(fp,*) 'BEFORE SORT:',ninter,nchord
            WRITE(fp,*) dinter(1:ninter)
            WRITE(fp,*) ointer(1:ninter)
            WRITE(fp,*) sinter(1:ninter)
          ENDIF

      IF (ninter.GT.1) THEN
c...    Sort intersections with respect to distance from the viewing location:
        cont = .TRUE.
        DO WHILE (cont)  ! *** only if niter.gt.1...
          cont = .FALSE.
          DO i1 = 1, ninter-1
            DO i2 = i1+1, ninter
              IF (dinter(i1).GT.dinter(i2)) THEN
                dinter(0)  = dinter(i2)
                dinter(i2) = dinter(i1)
                dinter(i1) = dinter(0)
                vinter(1:3,0)  = vinter(1:3,i2)  
                vinter(1:3,i2) = vinter(1:3,i1)  
                vinter(1:3,i1) = vinter(1:3,0) 
                cont = .TRUE.

                ointer(0)  = ointer(i2)
                ointer(i2) = ointer(i1)
                ointer(i1) = ointer(0)
                sinter(0)  = sinter(i2)
                sinter(i2) = sinter(i1)
                sinter(i1) = sinter(0)
                srfinter(0)  = srfinter(i2)   ! Put all these into a structure...!
                srfinter(i2) = srfinter(i1)
                srfinter(i1) = srfinter(0)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
c...    Check if 2 are the same, which can happen if the chord cuts at a seam:
        DO i1 = ninter-1, 1, -1
          IF (DABS(dinter(i1)-dinter(i1+1)).LT.1.0D-05) THEN
            DO i2 = i1+1, ninter-1
              dinter(i2)  = dinter(i2+1)
              ointer(i2)  = ointer(i2+1)
              sinter(i2)  = sinter(i2+1)
              srfinter(i2)  = srfinter(i2+1)
c              sinter(i1)  = sinter(i2+1)               ! BUG! Well, about 1 day to find this one... 19.04.06 -SL
c              srfinter(i1)  = srfinter(i2+1)
              vinter(1:3,i2) = vinter(1:3,i2+1)  
            ENDDO
            ninter = ninter - 1
c            WRITE(0,*) 'WARNING! CUTTING INTERSECTION...'
          ENDIF
        ENDDO          

c        IF (mode.EQ.2) THEN 
c          WRITE(0,'(A,I3,20(F8.5,I7,I3:))') 
c     .      'CUTS2:',ninter,
c     .      (dinter(i1),ointer(i1),sinter(i1),i1=1,ninter)
c        ENDIF
      ENDIF
 

      IF (nchord.EQ.dchord) THEN
        WRITE(fp,*) 'AFTER SORT:',ninter,nchord
        WRITE(fp,*) dinter(1:ninter)
        WRITE(fp,*) ointer(1:ninter)
        WRITE(fp,*) sinter(1:ninter)
      ENDIF

      IF     (mode.EQ.IT_VWINTER) THEN
c...    Wall intersections:
        IF (ninter.GE.1) THEN
c        IF (ninter.GE.2) THEN  ! *MARK*
          nvwinter = ninter
          vwinter(1:ninter)%dist = dinter(1:ninter)
          vwinter(1:ninter)%obj  = ointer(1:ninter)
          vwinter(1:ninter)%sur  = sinter(1:ninter)
          vwinter(1:ninter)%srf  = srfinter(1:ninter)
          DO i1 = 1, ninter
            vwinter(i1)%v(1:3) = vinter(1:3,i1)     
          ENDDO
          IF (nchord.EQ.dchord) THEN
            WRITE(fp,*) 'NUMBER OF WALL INTESECTIONS:',ninter,nchord
            WRITE(fp,*) dinter(1:ninter)
            WRITE(fp,*) ointer(1:ninter)
            WRITE(fp,*) sinter(1:ninter)
          ENDIF
        ELSE
          WRITE(0,*) 'NUMBER OF WALL INTESECTIONS:',ninter,nchord
          WRITE(0,*) 'EXPECTING AT LEAST 1 INTERSECTION, TROUBLE'
          status = -1
          RETURN
c          CALL ER('FindSurfaceIntersections','Expecting 1 wall '//
c     .            'intersection',*99)
        ENDIF 

      ELSEIF (mode.EQ.IT_GBINTER) THEN
c...    Magnetic grid boundary surfaces:
        IF (.TRUE.) THEN
          IF (nchord.EQ.dchord) THEN
            WRITE(fp,*) 'NUMBER OF GRID BOUNDARY INTERSECTIONS:',
     .      ninter,nchord
          ENDIF
          ngbinter = ninter
          gbinter(1:ninter)%dist = dinter(1:ninter)
          gbinter(1:ninter)%obj  = ointer(1:ninter)
          gbinter(1:ninter)%sur  = sinter(1:ninter)
          gbinter(1:ninter)%srf  = srfinter(1:ninter)
          DO i1 = 1, ninter
            gbinter(i1)%v(1:3) = vinter(1:3,i1)     
          ENDDO
c          WRITE(0,*) 'NUMBER OF GRID INTESECTIONS:',ninter,nchord
c          WRITE(0,*) '                            ',dinter(1:ninter)
        ELSE
        ENDIF

      ELSEIF (mode.EQ.IT_OBINTER) THEN
c...    Surfaces bounding the current object:
        IF (.TRUE.) THEN
c        IF (ninter.EQ.1) THEN
          IF (nchord.EQ.dchord) THEN
            WRITE(fp,*) 'NUMBER OF GRID INTERSECTIONS:',
     .      ninter,nchord
          ENDIF
          nobinter = ninter
          obinter(1:ninter)%dist = dinter(1:ninter)
          obinter(1:ninter)%obj  = ointer(1:ninter)
          obinter(1:ninter)%sur  = sinter(1:ninter)
          obinter(1:ninter)%srf  = srfinter(1:ninter)
          DO i1 = 1, ninter
            obinter(i1)%v(1:3) = vinter(1:3,i1)     
          ENDDO
c...      Check that only 1 surface is identified, unless the surface
c         in question is of geometry type GT_TC ('toroidally conintuous')
c         which can be curved in x,y,z space:
          count = 0
          DO i1 = 1, nobinter
            IF (obj(obinter(i1)%obj)%gsur(obinter(i1)%sur).NE.GT_TC)
     .        count = count + 1
          ENDDO
          IF (nobinter.GT.0.AND.count.LE.2) THEN
c          IF (nobinter.GT.0.AND.count.LE.1) THEN
c...        All is good:
          ELSE
c            WRITE(0,*) 'NOBINTER PROBLEM:',ninter
            status = -2
c            RETURN
c            CALL ER('FindSurfaceIntersections','Expecting 1 intersect'//
c     .              'ion for each object',*99)
          ENDIF
        ENDIF

      ELSE
        CALL ER('FindSurfaceIntersections','Unknown mode',*99)
      ENDIF



      RETURN
 99   WRITE(0,*) '   N     =',n
      WRITE(0,*) '   NINTER=',ninter,MAXINTER,mode
      STOP
      END
