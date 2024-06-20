C@PROCESS OPT(3) NOSDUMP NOGOSTMT
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
************************************************************************
       SUBROUTINE gzclp1(n,xyz,nrm,tedg,  abcd)
************************************************************************
*
*  dummy arguments:
*      n - number of input and output polygon vertices;
*      xyz - coordinates of input output polygon vertices;
*      nrm - normals of input output polygon vertices;
*      tedg - presence of input output polygon;
*
*
*      abcd - clipping plane parameters;
*
*  function: this subroutine clipps the input polygon against clipping
*            plane produceing the output polygon;
*
*  CHANGED. M. BUSCH    12.4.91
*                       INTEGER*4  --> INTEGER
*                       REAL   *4  --> REAL
*  CHANGED. M. BUSCH    6.1.93 Common Bloecke nicht mehr ueber include
*                       sondern als Source eingefuegt
*  CHANGED. M. BUSCH    10.3.93 MINREAL, MAXREAL +- 1.e.38 statt
*                               +-1.e50 (Anpassung an IEEE Format)
************************************************************************
       IMPLICIT NONE
       INTEGER    M
       INTEGER    INSIDE,OUT
       PARAMETER (M=12)
       PARAMETER (OUT=1, INSIDE=2)
*
*                          declare dummy arguments
       INTEGER    n
       REAL       xyz(1:3,1:M), nrm(1:3,1:M)
       LOGICAL    tedg(1:M)
*
       REAL       abcd(1:4)
*-----------------------------------------------------------------------
*                          declare local variables
       INTEGER    gzinsd
       INTEGER    np
       REAL       xyzp(1:3,1:M), nrmp(1:3,1:M)
       LOGICAL    pedg(1:M)
*
*
*
       REAL       point(1:3), norm(1:3)
       INTEGER    first, next
       LOGICAL    p1, in
       INTEGER    i,j
*----------------------------------------------------------------------
*
       np=0
       first=1
*
*
*
       IF(gzinsd(xyz(1,first), abcd).EQ.INSIDE) THEN
*
          in=.TRUE.
          p1=.TRUE.
          np=np+1
          CALL gzoutp(xyzp(1,np), xyz(1, first))
          CALL gzoutp(nrmp(1,np), nrm(1, first))
       ELSE
          in=.FALSE.
          p1=.FALSE.
       ENDIF
*
       DO 100 next=2,n
         IF(in) THEN
           IF(gzinsd(xyz(1, next), abcd).EQ.INSIDE) THEN
*
             np=np+1
             CALL gzoutp(xyzp(1, np), xyz(1, next))
             CALL gzoutp(nrmp(1, np), nrm(1, next))
*
             pedg(np-1)=tedg(next-1)
*
           ELSE
             in=.FALSE.
             CALL gzints(point, xyz(1,first),xyz(1,next),
     +                   norm,  nrm(1,first),nrm(1,next), abcd)
*
             np=np+1
             CALL gzoutp(xyzp(1, np), point)
             CALL gzoutp(nrmp(1, np), norm)
*
             pedg(np-1)=tedg(next-1)
           ENDIF
*
       ELSE
           IF(gzinsd(xyz(1,next), abcd).EQ.INSIDE) THEN
             in=.TRUE.
             CALL gzints(point, xyz(1,first),xyz(1,next),
     +                   norm,  nrm(1,first),nrm(1,next), abcd)
*
             np=np+1
             CALL gzoutp(xyzp(1, np), point)
             CALL gzoutp(nrmp(1, np), norm)
*
             IF(np-1.NE.0) pedg(np-1)=.FALSE.
*
*
             np=np+1
             CALL gzoutp(xyzp(1, np), xyz(1, next))
             CALL gzoutp(nrmp(1, np), nrm(1, next))
*
             pedg(np-1)=tedg(next-1)
*
           ENDIF
       ENDIF
       first=next
 100   CONTINUE
*
*
*
       IF(n.EQ.0) RETURN
*
       IF(in.EQV.p1) THEN
          IF(in) pedg(np)=tedg(n)
       ELSE
             CALL gzints(point, xyz(1,n),xyz(1,1),
     +                   norm,  nrm(1,n),nrm(1,1), abcd)
*
             np=np+1
             CALL gzoutp(xyzp(1, np), point)
             CALL gzoutp(nrmp(1, np), norm)
*
             IF(p1) THEN
                pedg(np-1)=.FALSE.
                pedg(np)=tedg(n)
             ELSE
                pedg(np-1)=tedg(n)
                pedg(np)=.FALSE.
             ENDIF
       ENDIF
*
*
*
       n=np
*
       IF(n.EQ.0) RETURN
       DO 120 i=1,n
         DO 130 j=1, 3
           xyz(j,i)=xyzp(j,i)
           nrm(j,i)=nrmp(j,i)
 130     CONTINUE
         tedg(i)=pedg(i)
 120   CONTINUE
*
*
*---------------------------------------------------------  end gzclp1
       END
*
*
*
       SUBROUTINE gzoutp(a,b)
       REAL       a(1:3),b(1:3)
       a(1)=b(1)
       a(2)=b(2)
       a(3)=b(3)
       END
*
************************************************************************
*
*
*
*
C@PROCESS OPT(3) NOSDUMP NOGOSTMT
       INTEGER FUNCTION gzinsd(a,abcd)
*
       INTEGER   IN,OUT
       PARAMETER (OUT=1,IN=2)
       REAL        a(1:3), abcd(1:4), aa
*
*
       aa=(a(1)*abcd(1)+a(2)*abcd(2)+a(3)*abcd(3)+abcd(4))
       IF(aa.LT.0.0) THEN
        gzinsd=OUT
       ELSE
        gzinsd=IN
       ENDIF
*
       RETURN
       END
*
*
************************************************************************
*
*
C@PROCESS OPT(3) NOSDUMP NOGOSTMT
       SUBROUTINE gzints(point, p1,p2,  norm,n1,n2,  abcd)
************************************************************************
*
*  dummy arguments:
*      point - intesection point of a line segment with a plane;
*      p1,p2 - ednpoits of the line segment;
*      norm - normal of intersection point;
*      n1,n2 - normals of points p1,p2;
*      abcd - plane parameters;
*
*  function:
*      this subroutine computes an intersection point of specifed
*      line segment and specifed plane;
*
************************************************************************
       IMPLICIT NONE
*                          declare dummy arguments
       REAL       point(1:3), p1(1:3), p2(1:3)
       REAL       norm(1:3),  n1(1:3), n2(1:3)
       REAL       abcd(1:4)
*
*-----------------------------------------------------------------------
       REAL       gzunity
       REAL       mi, nmi
*
*
       mi=(abcd(1)*p1(1) + abcd(2)*p1(2) + abcd(3)*p1(3) + abcd(4))/
     +    (abcd(1)*(p1(1)-p2(1)) + abcd(2)*(p1(2)-p2(2)) +
     +     abcd(3)*(p1(3)-p2(3)) )
*
       point(1)=p1(1) + mi*(p2(1)-p1(1))
       point(2)=p1(2) + mi*(p2(2)-p1(2))
       point(3)=p1(3) + mi*(p2(3)-p1(3))
*
*
       nmi=1.0-mi
*
       norm(1)=n1(1)     + mi*(n2(1)-n1(1))
       norm(2)=n1(2)     + mi*(n2(2)-n1(2))
       norm(3)=n1(3)     + mi*(n2(3)-n1(3))
       nmi=gzunity(norm)
*
*------------------------------------------------------   end gzints
       END
*
*
************************************************************************
*
C@PROCESS OPT(3) NOSDUMP NOGOSTMT
       SUBROUTINE gzclpp(trg, nrm, edg,   n,trgs,nrms,edgs,  view)
************************************************************************
*
*  dummy arguments:
*      trg - vertices coordinates of the input triangle;
*      nrm - normals of triangle vertex;
*      edg - presence of triangle edges;
*
*      n - triangle number (output);
*      trgs - vertices coordinates of triangles (output);
*      nrms - normals of vertices of triangles (output);
*      edgs - presence of edges of triangles (output);
*
*      view - view index;
*
*  function:
*      this subroutine clips the specifed triangle against clipping
*      volume of specifed view producing a set of output triangles
*      which represent a visible pice of the input triangle;
*
************************************************************************
       IMPLICIT NONE
*  include global constants and variables
c      include(gzcst11)
***********************************************************************
*
*    autor: Zenon Zowierucha, zdv069, 20.Nov.1989
*
*    file name, type :  GZCST1 COPY
*
*    function :  this file defines global constants and is destinated
*                for use with INCLUDE directive;
*
*   update: Marlene Busch, zdv131, 10. 4.91
*           DEVIND = 9 --> DEVIND =81
***********************************************************************
*----------------------------------------------------------------------
*  declare constants
*----------------------------------------------------------------------
*
*                         maximal number of vertices of all faces
       INTEGER      VINFACE
*
*                    maximal number of edges that belong to one vertex
       INTEGER      EOFVERT
*
*                    maximal number of faces that belong to one vertex
       INTEGER      FOFVERT
*
*                         maximal number of views
       INTEGER      MAXVIEW
*
*                         maximal number of face colors (hues)
       INTEGER      MAXFCLR
*
*
*                         maximal number of pixels on the x-axis
       INTEGER      MAXPKX
*
*                         maximal number of pixels on the y-axis
       INTEGER      MAXPKY
*
*                          maximal gloss value
       INTEGER      MAXGLS
*
*                          maximal number of objects
       INTEGER      MAXOBJ
*
*                          maximal number of colors in look up table
       INTEGER      MAXCLR
*
*                          error file identification number
       INTEGER      DEVIND
*
*----------------------------------------------------------------------
*  seting values
*----------------------------------------------------------------------
*
*
       PARAMETER (VINFACE=3   )
       PARAMETER (EOFVERT=8   )
       PARAMETER (FOFVERT=8   )
*
       PARAMETER (MAXVIEW=32  )
*
       PARAMETER (MAXFCLR=60  )
*
       PARAMETER (MAXPKX =1024)
       PARAMETER (MAXPKY =1024)
       PARAMETER (MAXGLS =32  )
*
       PARAMETER (MAXOBJ =1200 )
       PARAMETER (MAXCLR =128 )
       PARAMETER (DEVIND= 81  )
c      include(gzcst21)
***********************************************************************
*
*  autor: Zenon Zowierucha, zdv069
*
*  file name, type :   GZCST2 COPY
*
*  function :  this file defines global constants;
*              this file is destined for use with INCLUDE directive;
*
*
***********************************************************************
*----------------------------------------------------------------------
*  declare constants
*----------------------------------------------------------------------
*                          projection types
       INTEGER    PARALLEL
       INTEGER    PERSPECTIVE
*
*                          viewing styles
       INTEGER    WIRE
       INTEGER    CONST
       INTEGER    WCONST
       INTEGER    GWSHADED
       INTEGER    GSHADED
       INTEGER    PWSHADED
       INTEGER    PSHADED
*
*                          light source types
       INTEGER    SUN
       INTEGER    POINT
*
       INTEGER    OFF
       INTEGER    ON
*
       INTEGER    PRECONCATENATE
       INTEGER    POSTCONCATENATE
       INTEGER    REPLACE
*
       INTEGER    FRONT
       INTEGER    BACK
*
*
       REAL       MINREAL
       REAL       MAXREAL
*
       INTEGER    ONLINE
       INTEGER    OFFLINE
*
*
*
*----------------------------------------------------------------------
*  seting values
*----------------------------------------------------------------------
*
       PARAMETER (PARALLEL    =1)
       PARAMETER (PERSPECTIVE=2)
*
       PARAMETER (WIRE   =1)
       PARAMETER (CONST   =2)
       PARAMETER (WCONST  =3)
       PARAMETER (GWSHADED=4)
       PARAMETER (GSHADED =5)
       PARAMETER (PWSHADED=6)
       PARAMETER (PSHADED =7)
*
       PARAMETER (SUN  =1)
       PARAMETER (POINT=2)
*
       PARAMETER (OFF=1)
       PARAMETER (ON =2)
*
       PARAMETER (PRECONCATENATE= 1)
       PARAMETER (POSTCONCATENATE=2)
       PARAMETER (REPLACE=        3)
*
       PARAMETER (FRONT=1)
       PARAMETER (BACK =2)
*
*
       PARAMETER (MINREAL =-1.0E38)
       PARAMETER (MAXREAL = 1.0E38)
*
       PARAMETER (ONLINE = 1)
       PARAMETER (OFFLINE= 2)
c      include(gzvchr)
***********************************************************************
*
*  autor:  Zenon Zowierucha, zdv069
*
*  file name type : GZVCHR COPY
*
*  function: this file declares global variables, which describe a view
*            characteristic; this file is destinated for use with
*            INCLUDE directive;
*
*
***********************************************************************
*
*                          activity of view (viewactiv);
*                          vactiv(i) is equal TRUE if view nr i is
*                          activ or equal FALSE in other case;
       LOGICAL   vactiv(1:MAXVIEW)
*
*                          view style (WIRE,HIDDEN,WSHADED or SHADED);
       INTEGER    vstyle(1:MAXVIEW, 1:MAXOBJ)
*
*                          near clipping flag (OFF or ON);
       INTEGER    nearfl(1:MAXVIEW)
*
*                          far clipping flag (OFF or ON);
       INTEGER    farfl(1:MAXVIEW)
*
*                          of background color of view nr. i);
       REAL       bcgclr(1:3, 1:MAXVIEW)
*
*                          background color index;
       INTEGER    bcgind(1:MAXVIEW)
*
*                          viewport border color
       REAL       brdclr(1:3, 1:MAXVIEW)
*
*                          viewport border color index;
       INTEGER    brdind(1:MAXVIEW)
*
*                          viewport border activity (brdact(i) is equal
*                          ON if border for viewport nr. i is
*                          visible or is equal OFF in other case);
       INTEGER    brdact(1:MAXVIEW)
*
*
*
       COMMON/gzvchr/ vstyle, nearfl, farfl,
     +                bcgclr, brdclr, bcgind, brdind, brdact, vactiv
       save /gzvchr/
*
c      include(gzvmpg)
***********************************************************************
*
*   autor:  Zenon Zowierucha, zdv069
*
*   file name, type :  GZVMPG COPY
*
*   function :  this file declares global variables which define
*               view mapping characteristic;
*               this file is destinated for use with INCLUDE directive;
*
*
***********************************************************************
*
*                          window in view coordinate system (VCS),
*                          window(1,i), window(2,i), window(3,i),
*                          window(4,i) are xmin, ymin, xmax, ymax
*                          of window of view number i
       REAL       window(1:4, 1:MAXVIEW)
*
*                          viewport in device coordinate system (DCS),
*                          vieprt(1,i), vieprt(2,i), vieprt(3,i),
*                          vieprt(4,i) are respectivly xmin, ymin,
*                          xmax, ymax of viewport of view number i
       REAL       vieprt(1:4, 1:MAXVIEW)
*
*                          projection type (PARALEL or PERSPECTIVE)
       INTEGER    prtype(1:MAXVIEW)
*
*                          projection reference point (VCS)
       REAL       prp(1:3, 1:MAXVIEW)
*
*                          projection plane distance (VCS)
*                          (from the origin of the coordinate system)
       REAL       dist(1:MAXVIEW)
*
*                          near clipping plane distance (VCS)
*                          (from the orogin of the coordinate system)
       REAL       near(1:MAXVIEW)
*
*                          far clipping plane distance (VCS)
*                          (from the origin of the coordinate system)
       REAL       far(1:MAXVIEW)
*
*                          parameters of up, left, down and right
*                          clipping planes equations; abcd(1,i,j)*x +
*                          abcd(2,i,j)*y+abcd(3,i,j)*z + abcd(4,i,j)
*                          (i= 1-up, 2-left, 3-down, 4-right clipping
*                          plane; j is view index);
       REAL       abcd(1:4, 1:4, 1:MAXVIEW)
*
*
       COMMON /gzvmpg/window, vieprt,
     +              prp, dist, near, far,
     +              prtype, abcd
*
       SAVE  /gzvmpg/
       INTEGER    M
       PARAMETER (M=12)
*
*  declare dummy arguments
       REAL       trg(1:3,1:3), nrm(1:3,1:3)
       LOGICAL    edg(1:3)
       INTEGER    n
       REAL       trgs(1:3,1:M), nrms(1:3,1:M)
       LOGICAL    edgs(1:M)
       INTEGER    view
*
*  declare local variables
       INTEGER    nk
       REAL       xyz(1:3,1:M), nrml(1:3,1:M), abcd1(1:4,1:4)
       INTEGER    i,j,k
       LOGICAL    tedg(1:M)
       REAL       xleft
       INTEGER    ii,il,ip
*
*-----------------------------------------------------------------------
*
       DO 100 i=1,3
         DO 110 j=1,3
           xyz(j,i)=trg(j,i)
           nrml(j,i)=nrm(j,i)
 110     CONTINUE
         tedg(i)=edg(i)
 100   CONTINUE
*
*
       nk=3
*
       IF(prtype(view).EQ.PARALLEL) THEN
*
         abcd1(1,1)=0.0
         abcd1(2,1)=-1.0
         abcd1(4,1)= window(4,view)
*
         abcd1(1,2)=1.0
         abcd1(2,2)= 0.0
         abcd1(4,2)=-window(1,view)
*
*
         abcd1(1,3)=0.0
         abcd1(2,3)= 1.0
         abcd1(4,3)=-window(2,view)
*
         abcd1(1,4)=-1.0
         abcd1(2,4)= 0.0
         abcd1(4,4)= window(3,view)
*
         DO 120 k=1, 4
           abcd1(3,k)=0.0
           CALL gzclp1(nk, xyz, nrml, tedg,  abcd1(1,k))
*
           IF(nk.EQ.0) THEN
                  n=0
                  RETURN
           ENDIF
 120     CONTINUE
*
*
       ELSE
*      ---------------------------------    perspective projection
           DO 150 k=1,4
             CALL gzclp1(nk, xyz, nrml, tedg, abcd(1,k,view))
*
             IF(nk.EQ.0) THEN
                  n=0
                  RETURN
             ENDIF
*
 150       CONTINUE
*
       ENDIF
*
*
       IF(nearfl(view).EQ.ON) THEN
*                          front clipping
         abcd1(1,1)=0.0
         abcd1(2,1)=0.0
         abcd1(3,1)=-1.0
         abcd1(4,1)= near(view)
*
         CALL gzclp1(nk, xyz, nrml, tedg,  abcd1(1,1))
*
         IF(nk.EQ.0) THEN
                n=0
                RETURN
         ENDIF
       ENDIF
*
*
*
       IF(farfl(view).EQ.ON) THEN
*                          back clipping
         abcd1(1,1)=0.0
         abcd1(2,1)=0.0
         abcd1(3,1)=1.0
         abcd1(4,1)=-far(view)
*
         CALL gzclp1(nk, xyz, nrml, tedg,  abcd1(1,1))
*
         IF(nk.EQ.0) THEN
                n=0
                RETURN
         ENDIF
       ENDIF
*
*--------------------------------    split polygon xyz into triangles
*
       IF(nk.EQ.3) THEN
          n=1
          DO 2 i=1,3
            trgs(i,1)=xyz(i,1)
            trgs(i,2)=xyz(i,2)
            trgs(i,3)=xyz(i,3)
            nrms(i,1)=nrml(i,1)
            nrms(i,2)=nrml(i,2)
            nrms(i,3)=nrml(i,3)
            edgs(i)=tedg(i)
 2        CONTINUE
       ELSE
       n=0
 1000  xleft=xyz(1,1)
       ii=1
*
       DO 190 j=2, nk
         IF(xleft.GT.xyz(1,j)) THEN
            xleft=xyz(1,j)
            ii=j
         ENDIF
 190   CONTINUE
*
       il=ii-1
       IF(il.EQ.0) il=nk
       ip=ii+1
       IF(ip.GT.nk) ip=1
*
       k=n*3 + 1
       DO 200 i=1,3
          trgs(i, k) = xyz(i, il)
          trgs(i, k+1) = xyz(i, ii)
          trgs(i, k+2) = xyz(i, ip)
*
          nrms(i, k)=nrml(i, il)
          nrms(i, k+1)=nrml(i, ii)
          nrms(i, k+2)=nrml(i, ip)
 200   CONTINUE
*
       edgs(k)=tedg(il)
       edgs(k+1)=tedg(ii)
       IF(nk.GT.3) THEN
          edgs(k+2)=.FALSE.
       ELSE
          edgs(k+2)=tedg(ip)
       ENDIF
       tedg(il)=.FALSE.
*
       n=n+1
       nk=nk-1
*
       DO 210 k=ii, nk
          DO 220 j=1,3
*
             xyz(j, k)=xyz(j, k+1)
             nrml(j, k)=nrml(j, k+1)
 220      CONTINUE
*
          tedg(k)=tedg(k+1)
 210   CONTINUE
*
       IF(nk.GT.2) GOTO 1000
       ENDIF
*----------------------------------------------------   end gzclpp
       END
