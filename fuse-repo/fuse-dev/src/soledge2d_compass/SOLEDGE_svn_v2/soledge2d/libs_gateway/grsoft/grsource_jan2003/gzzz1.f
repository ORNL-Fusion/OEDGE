C@PROCESS OPT(3) NOSDUMP NOGOSTMT
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
       SUBROUTINE gztrg(n, mem, tria, visi, objind, noroom)
*
************************************************************************
*
*  autor: Zenon Zowierucha, zdv069, 23.Nov.1989
*  file name, type:    GZZZ FORTRAN
*
*  dummy arguments:
*              n - number of user supplyed memory cells;
*
*              mem - user supplyed memory area;
*
*              tria   - a triangle, tria(1,i),tria(2,i),tria(3,i) are
*                    x,y,z coordinates of vertex nr i of the triangle;
*              visi   - visibility value of triangle edges;
*                     the edge nr. i  is visible, if visi(i) is equal
*                     1  or is invisible if  visi(i)  is equal  0;
*
*              objind - object index;
*
*  function: this subroutine changes the triangle informations into
*            more suitable form (see file GZPTCH COPY);
*
*  called by:  main program
*
*  CHANGED. M. BUSCH    12.4.91
*                       INTEGER*4  --> INTEGER
*                       REAL   *4  --> REAL
*  CHANGED. M. BUSCH    6.1.93 Common Bloecke nicht mehr ueber include
*                       sondern als Source eingefuegt
*  CHANGED. M. BUSCH    10.3.93 MINREAL, MAXREAL +- 1.e.38 statt
*                               +-1.e50 (Anpassung an IEEE Format)
*
************************************************************************
       IMPLICIT NONE
*-----------------------------------------------------------------------
*  include global constants and variables declarations
*-----------------------------------------------------------------------
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
c      include(gzptch1)
***********************************************************************
*
*   autor: Zenon Zowierucha, zdv069, 20.Nov.1989
*
*   file name, type :  GZPTCH COPY
*
*   function :  this file declares global variables which describe
*               geometral and topological dates of patches;
*               this file is destinated for use with INCLUDE directive;
*
*
***********************************************************************
*
*                          index of the last object in database;
       INTEGER    lastobj
*
*                          number of vertices of patches
*                          (vertno(i) is number of vertices of patch
*                           nr. i);
       INTEGER    vertno(1:MAXOBJ)
*
*                          patch indices (vrtind(i) is the index of the
*                          first vertex of patch (object) nr. i);
       INTEGER    vrtind(1:MAXOBJ)
*
*
*                          index of the last vertex in table vertex;
       INTEGER    lastvrt
*
*
*                          number of edges of patches;
       INTEGER    edgeno(1:MAXOBJ)
*
*                          patch indices (edgind(i) is the index of the
*                          first edge of patch (object) nr. i)
       INTEGER    edgind(1:MAXOBJ)
*
*
*                          index of the last edge in table edges;
       INTEGER    lastedg
*
*
*                          number of faces of patches;
       INTEGER    faceno(1:MAXOBJ)
*
*                          face indexes (fceind(i) is index of the
*                          first face of object nr. i);
       INTEGER    fceind(1:MAXOBJ)
*
*
*                          index of the last face in table faces;
       INTEGER    lastfce
*
*                          list of faces which belong to a vertex;
*                          fcevrt(1,j),fcevrt(1,k1),...,fcevrt(1,ki)
*                          are indexes of faces which belong to vertex
*                          nr j; k1=fcevrt(2,j),...,ki=fcevrt(2,k(i-1))
*                          until ki=0 ;
       INTEGER    lastfv
*
*
*
*
*                          object activity (objact(i,j) is .TRUE. if
*                          object nr j is visible (activ) in view nr i
*                          or is equal .FALSE. in other case);
       LOGICAL    objact(1:MAXVIEW, 1:MAXOBJ)
*
*
*                          object presence (objprs(i) is equal .TRUE.
*                          if object nr. i is in memory (is defined) or
*                          is equal .FALSE. in other case);
       LOGICAL    objprs(1:MAXOBJ)
*
*                          begin of vertex table in storage defined
*                          by the user;
       INTEGER    vrtbeg
*
*                          begin of the edge table;
       INTEGER    edgbeg
*
*                          begin of edge presence table;
       INTEGER    prsbeg
*
*                          begin of the face table;
       INTEGER    fcebeg
*
*                          begin of face normal table;
       INTEGER    fvbeg
*
*                          begin of vertex normal table;
       INTEGER    vnbeg
       INTEGER    fnbeg
*                          maximal number of vertices;
       INTEGER    maxvrt
*
*                          maximal number of edges;
       INTEGER    maxedge
*
*                          maximal number of faces;
       INTEGER    maxface
*
*
*
       COMMON /gzptch/vertno, edgeno, faceno, lastobj,
     +                lastvrt, lastedg, lastfce,
     +                lastfv, fceind, vrtind, edgind,
     +                vrtbeg, edgbeg, prsbeg, fcebeg, fvbeg, vnbeg,
     +                fnbeg,
     +                maxvrt, maxedge, maxface,
     +                objact, objprs
      SAVE /gzptch/
*
*
*
*
*
*
*                          faces of patches; faces(1,i), faces(2,i),
*                          ...,faces(VINFACE,i) are indexes of vertexes
*                          which build a face number i;
*                          faces(VINFACE+1,i),...,faces(2*VINFACE,i)
*                          are indexes of edges which build face nr i;
*      INTEGER    faces(1:2*VINFACE, 1:MAXFACE)
*
*                          edges of patches; edges(1,i), edges(2,i)
*                          are indices of vertices which are joined by
*                          edge number i; edges(3,i), edges(4,i) are
*                          indices of faces which are divided by edge
*                          number i;
*      INTEGER    edges(1:4, 1:MAXEDGE)
*
*                          coordinates of patch vertices;
*                          vertex(1,i),vertex(2,i),vertex(3,i) = x,y,z
*                          coordinates of vertex number i;
*      REAL       vertex(1:3, 1:MAXVERT)
*
*                          list of faces which belong to a vertex;
*                          fcevrt(1,j),fcevrt(1,k1),...,fcevrt(1,ki)
*                          are indexes of faces which belong to vertex
*                          nr j; k1=fcevrt(2,j),...,ki=fcevrt(2,k(i-1))
*                          until ki=0 ;
*      INTEGER    fcevrt(1:2, 1:FOFVERT*MAXVERT)
*
*
*                          normal vectors for vertices;,
*                          vnrm(1,i),vnrm(2,i),vnrm(3,i) are x,y,z
*                          components of normal vector for vertex nr i;
*      REAL       vnrm(1:3, 1:MAXVERT)
*                          normal vectors for faces, fnrm(1,i),
*                          fnrm(2,i),fnrm(3,i) are x,y,z components
*                          of normal vector for face number i;
*      REAL       fnrm(1:3, 1:MAXFACE)
*                          presence of edges; edgprs(i) is equal.TRUE.
*                          if edge number i is visible or equal .FALSE.
*                          in other case;
*      LOGICAL    edgprs(1:MAXEDGE)
*
*-----------------------------------------------------------------------
*  declare dummy arguments
*-----------------------------------------------------------------------
*
       INTEGER    n
       INTEGER    mem(n)
C Busch 25.6.99
       REAL rmem(n)
       LOGICAL lmem(n)
C Busch
       REAL       tria(1:3, 1:3)
       INTEGER    visi(1:3)
       INTEGER    objind
       INTEGER    noroom
*
*
*      write(*,*)vrtbeg, edgbeg, prsbeg, fcebeg, fvbeg
C Busch 25.6.99
       CALL gz1trg(tria, visi, objind,
     +             rmem(vrtbeg), mem(edgbeg), lmem(prsbeg), mem(fcebeg),
     +             mem(fvbeg), noroom)
*
       RETURN
       END
*
*
*----------------------------------------------------------------------
*
C@PROCESS OPT(3) NOSDUMP NOGOSTMT
       SUBROUTINE gz1trg(tria, visi, objind,
     +                   vertex, edges, edgprs, faces, fcevrt, noroom)
*
*
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
c      include(gzptch1)
***********************************************************************
*
*   autor: Zenon Zowierucha, zdv069, 20.Nov.1989
*
*   file name, type :  GZPTCH COPY
*
*   function :  this file declares global variables which describe
*               geometral and topological dates of patches;
*               this file is destinated for use with INCLUDE directive;
*
*
***********************************************************************
*
*                          index of the last object in database;
       INTEGER    lastobj
*
*                          number of vertices of patches
*                          (vertno(i) is number of vertices of patch
*                           nr. i);
       INTEGER    vertno(1:MAXOBJ)
*
*                          patch indices (vrtind(i) is the index of the
*                          first vertex of patch (object) nr. i);
       INTEGER    vrtind(1:MAXOBJ)
*
*
*                          index of the last vertex in table vertex;
       INTEGER    lastvrt
*
*
*                          number of edges of patches;
       INTEGER    edgeno(1:MAXOBJ)
*
*                          patch indices (edgind(i) is the index of the
*                          first edge of patch (object) nr. i)
       INTEGER    edgind(1:MAXOBJ)
*
*
*                          index of the last edge in table edges;
       INTEGER    lastedg
*
*
*                          number of faces of patches;
       INTEGER    faceno(1:MAXOBJ)
*
*                          face indexes (fceind(i) is index of the
*                          first face of object nr. i);
       INTEGER    fceind(1:MAXOBJ)
*
*
*                          index of the last face in table faces;
       INTEGER    lastfce
*
*                          list of faces which belong to a vertex;
*                          fcevrt(1,j),fcevrt(1,k1),...,fcevrt(1,ki)
*                          are indexes of faces which belong to vertex
*                          nr j; k1=fcevrt(2,j),...,ki=fcevrt(2,k(i-1))
*                          until ki=0 ;
       INTEGER    lastfv
*
*
*
*
*                          object activity (objact(i,j) is .TRUE. if
*                          object nr j is visible (activ) in view nr i
*                          or is equal .FALSE. in other case);
       LOGICAL    objact(1:MAXVIEW, 1:MAXOBJ)
*
*
*                          object presence (objprs(i) is equal .TRUE.
*                          if object nr. i is in memory (is defined) or
*                          is equal .FALSE. in other case);
       LOGICAL    objprs(1:MAXOBJ)
*
*                          begin of vertex table in storage defined
*                          by the user;
       INTEGER    vrtbeg
*
*                          begin of the edge table;
       INTEGER    edgbeg
*
*                          begin of edge presence table;
       INTEGER    prsbeg
*
*                          begin of the face table;
       INTEGER    fcebeg
*
*                          begin of face normal table;
       INTEGER    fvbeg
*
*                          begin of vertex normal table;
       INTEGER    vnbeg
       INTEGER    fnbeg
*                          maximal number of vertices;
       INTEGER    maxvrt
*
*                          maximal number of edges;
       INTEGER    maxedge
*
*                          maximal number of faces;
       INTEGER    maxface
*
*
*
       COMMON /gzptch/vertno, edgeno, faceno, lastobj,
     +                lastvrt, lastedg, lastfce,
     +                lastfv, fceind, vrtind, edgind,
     +                vrtbeg, edgbeg, prsbeg, fcebeg, fvbeg, vnbeg,
     +                fnbeg,
     +                maxvrt, maxedge, maxface,
     +                objact, objprs
      SAVE /gzptch/
*
*
*
*
*
*
*                          faces of patches; faces(1,i), faces(2,i),
*                          ...,faces(VINFACE,i) are indexes of vertexes
*                          which build a face number i;
*                          faces(VINFACE+1,i),...,faces(2*VINFACE,i)
*                          are indexes of edges which build face nr i;
*      INTEGER    faces(1:2*VINFACE, 1:MAXFACE)
*
*                          edges of patches; edges(1,i), edges(2,i)
*                          are indices of vertices which are joined by
*                          edge number i; edges(3,i), edges(4,i) are
*                          indices of faces which are divided by edge
*                          number i;
*      INTEGER    edges(1:4, 1:MAXEDGE)
*
*                          coordinates of patch vertices;
*                          vertex(1,i),vertex(2,i),vertex(3,i) = x,y,z
*                          coordinates of vertex number i;
*      REAL       vertex(1:3, 1:MAXVERT)
*
*                          list of faces which belong to a vertex;
*                          fcevrt(1,j),fcevrt(1,k1),...,fcevrt(1,ki)
*                          are indexes of faces which belong to vertex
*                          nr j; k1=fcevrt(2,j),...,ki=fcevrt(2,k(i-1))
*                          until ki=0 ;
*      INTEGER    fcevrt(1:2, 1:FOFVERT*MAXVERT)
*
*
*                          normal vectors for vertices;,
*                          vnrm(1,i),vnrm(2,i),vnrm(3,i) are x,y,z
*                          components of normal vector for vertex nr i;
*      REAL       vnrm(1:3, 1:MAXVERT)
*                          normal vectors for faces, fnrm(1,i),
*                          fnrm(2,i),fnrm(3,i) are x,y,z components
*                          of normal vector for face number i;
*      REAL       fnrm(1:3, 1:MAXFACE)
*                          presence of edges; edgprs(i) is equal.TRUE.
*                          if edge number i is visible or equal .FALSE.
*                          in other case;
*      LOGICAL    edgprs(1:MAXEDGE)
*
*  declare dummy arguments
*
       REAL       tria(1:3, 1:VINFACE)
       INTEGER    visi(1:3)
       INTEGER    objind
       REAL   vertex(1:3, *)
       INTEGER   edges(1:4, *)
       LOGICAL   edgprs(*)
       INTEGER   faces(1:6, *)
       INTEGER   fcevrt(1:2, *)
       INTEGER   noroom
*
*-----------------------------------------------------------------------
*  declare local variables
*-----------------------------------------------------------------------
       INTEGER     pointer(1:VINFACE)
*
       INTEGER     a, b, c, j, k, m, ii
*
       REAL        x, y, z
*
       LOGICAL     new(1:VINFACE)
       LOGICAL     present(1:3), ret,ok
*
       LOGICAL     donenrm(1:MAXOBJ)
       COMMON/gzdone/donenrm
       SAVE /gzdone/
*
       LOGICAL     OVERFLOW
       REAL        XMIN, YMIN, ZMIN, XMAX, YMAX,ZMAX

       DATA        OVERFLOW/.FALSE./
       DATA        XMIN/0.0    /
       DATA        ymin/0.0    /
       DATA        zmin/0.0    /
       DATA        xmax/0.0    /
       DATA        ymax/0.0   /
       DATA        zmax/0.0    /
*
*----------------------------------------------------------------------
       IF (OVERFLOW)  THEN
       noroom=1
       RETURN
       ENDIF
*
       ret=.FALSE.
*
       IF(objind.LT.1.0 .OR. objind.GT.MAXOBJ) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'  GZTRG INVALID OBJECT INDEX'
       ELSE
          IF(objind.NE.lastobj .AND. objprs(objind)) THEN
             ret=.TRUE.
       WRITE(DEVIND,*)'  GZTRG YOU CAN APPEND A TRIANGLE ONLY TO THE '
       WRITE(DEVIND,*)'LAST OBJECT OR CREATE NOT EXISTING OBJECT '
          ENDIF
       ENDIF
*
       IF(lastvrt+3.GE.maxvrt .OR. lastedg+3.GE.maxedge .OR.
     +    lastfce+1.GE.maxface) THEN
          overflow=.TRUE.
          noroom=1
          ret=.TRUE.
          WRITE(DEVIND,*)'  GZTRG TO MANY TRIANGLES, DATABASE IS FULL'
       ENDIF
*
*
       IF(ret) THEN
          WRITE(DEVIND,*)'  GZTRG NO ACTION PERFORMED '
          RETURN
       ENDIF
*
       noroom=0
       lastfce=lastfce+1
       IF(.NOT.objprs(objind)) THEN
          objprs(objind)=.TRUE.
          lastobj=objind
          vrtind(objind)=lastvrt+1
          edgind(objind)=lastedg+1
          fceind(objind)=lastfce
       ENDIF
*
       donenrm(objind)=.FALSE.
       faceno(objind)=faceno(objind)+1
         DO 2  j=1, VINFACE
*
          x = tria(1, j)
          y = tria(2, j)
          z = tria(3, j)
*
          ok=.FALSE.

          IF(x.GT.xmax) THEN
             xmax=x
             ok=.TRUE.
          ENDIF
          IF(x.LT.xmin) THEN
             xmin=x
             ok=.TRUE.
          ENDIF
          IF(y.GT.ymax) THEN
             ymax=y
             ok=.TRUE.
          ENDIF
          IF(y.LT.ymin) THEN
             ymin=y
             ok=.TRUE.
          ENDIF
          IF(z.GT.zmax) THEN
             zmax=z
             ok=.TRUE.
          ENDIF
          IF(z.LT.zmin) THEN
             zmin=z
             ok=.TRUE.
          ENDIF
*
          IF(ok) THEN
           new(j) = .TRUE.
           lastvrt= lastvrt+ 1
           vertno(objind)=vertno(objind)+1
           pointer(j) = lastvrt
           vertex(1, lastvrt)= x
           vertex(2, lastvrt)= y
           vertex(3, lastvrt)= z

          ELSE
*------------------------------------------------------  insert vertex
          k = lastvrt
 3        IF (k .EQ. vrtind(objind)-1) THEN
           new(j) = .TRUE.
           lastvrt= lastvrt+ 1
           vertno(objind)=vertno(objind)+1
           pointer(j) = lastvrt
           vertex(1, lastvrt)= x
           vertex(2, lastvrt)= y
           vertex(3, lastvrt)= z
          ELSE IF (x. NE .vertex(1, k)) THEN
              k = k-1
              GOTO 3
          ELSE IF (y .NE. vertex(2, k)) THEN
              k = k-1
              GOTO 3
          ELSE IF (z .NE. vertex(3, k)) THEN
              k = k-1
              GOTO 3
          ELSE
              pointer(j) = k
              new(j) = .FALSE.
          ENDIF
        ENDIF
*-----------------------------------  insert a new vertex of face nr i
         faces(j, lastfce)= pointer(j)
*
 2      CONTINUE
*--------------------   insert new face into fcevrt for three vertexes
       DO 80 j=1, VINFACE
         k=fcevrt(1, pointer(j))
         IF(k.EQ.0) THEN
            fcevrt(1, pointer(j))=lastfce
            fcevrt(2, pointer(j))=0
         ELSE
            k=pointer(j)
 1000       m=k
            k=fcevrt(2, k)
            IF(k.NE.0) GOTO 1000
            lastfv=lastfv+1

            IF (lastfv .GE. (vnbeg-fvbeg)/2 ) THEN
             overflow=.TRUE.
             noroom=1
             ret=.TRUE.
             WRITE(DEVIND,*)
     >                     '  GZTRG TO MANY TRIANGLES, DATABASE IS FULL'
             WRITE(DEVIND,*)'  GZTRG NO ACTION PERFORMED '
             RETURN

            ENDIF
            fcevrt(2, m)=lastfv
            fcevrt(1, lastfv)=lastfce
            fcevrt(2, lastfv)=0
         ENDIF
 80    CONTINUE
*-----------------------------------------------------------------------
        DO 1200 ii=1,3
         IF(visi(ii).NE.0) THEN
           present(ii)= .TRUE.
         ELSE
           present(ii)= .FALSE.
         ENDIF
 1200    CONTINUE
*
        k = lastedg
*
*----------------------------------------------------  insert new edges
        DO 4  j=1,VINFACE
*
         c = j+1
         IF(c .EQ. 4) c=1
*
*-------------------   consider an edge from vertex nr a to vertex nr b
         a = pointer(j)
         b = pointer(c)
*
         IF (new(j) .OR. new(c)) THEN
*         ------------------------------   this edge is not in edges yet
*                                          append this edge to edges
          edgeno(objind)=edgeno(objind)+1
          lastedg=lastedg+1
          edges(1, lastedg)= a
          edges(2, lastedg)= b
*         --------------------------   append current face for this edge
          edges(3, lastedg)=lastfce
          edges(4, lastedg)=0
*         ---------------------------   append this edge to current face
          faces(VINFACE+j,lastfce)=lastedg
*
          edgprs(lastedg) = present(j)
*
*         --------------------------------------------------------------
         ELSE
*
*         ----------------   establish if the edge was in table of edges
          k = lastedg
 5        IF (k .EQ. 0) THEN
*             ----------------------   the edge is not faund; append it
              lastedg= lastedg+1
              edgeno(objind)=edgeno(objind)+1
              edges(1,lastedg) = a
              edges(2,lastedg) = b
*             ------------------   append current face to table of edges
              edges(3,lastedg)=lastfce
              edges(4,lastedg)=0
*             -------------------   append the edge to table of faces
              faces(VINFACE+j,lastfce)=lastedg
*
              edgprs(lastedg) = present(j)
*             ----------------------------------------------------------
*
          ELSE
            IF (.NOT.((edges(1,k).EQ.a .AND. edges(2,k).EQ.b) .OR.
     +                (edges(1,k).EQ.b .AND. edges(2,k).EQ.a))) THEN
                           k = k-1
                           GOTO 5
            ELSE
*              ------------------   the edge is found in table of edges
               edges(4, k)=lastfce
               faces(VINFACE+j,lastfce)=k
               IF (present(j)) THEN
                           edgprs(k) = .TRUE.
               ENDIF
            ENDIF
          ENDIF
        ENDIF
*

 4     CONTINUE
       RETURN
       END
*--------------------------------------------------------   end gz1trg
*
*
*
*-----------------------------------------------------------------------
*
*
C@PROCESS OPT(3) NOSDUMP NOGOSTMT
       SUBROUTINE gz1nrml(vertex, faces, fcevrt, vnrm, fnrm, objind)
       IMPLICIT NONE
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
c      include(gzptch1)
***********************************************************************
*
*   autor: Zenon Zowierucha, zdv069, 20.Nov.1989
*
*   file name, type :  GZPTCH COPY
*
*   function :  this file declares global variables which describe
*               geometral and topological dates of patches;
*               this file is destinated for use with INCLUDE directive;
*
*
***********************************************************************
*
*                          index of the last object in database;
       INTEGER    lastobj
*
*                          number of vertices of patches
*                          (vertno(i) is number of vertices of patch
*                           nr. i);
       INTEGER    vertno(1:MAXOBJ)
*
*                          patch indices (vrtind(i) is the index of the
*                          first vertex of patch (object) nr. i);
       INTEGER    vrtind(1:MAXOBJ)
*
*
*                          index of the last vertex in table vertex;
       INTEGER    lastvrt
*
*
*                          number of edges of patches;
       INTEGER    edgeno(1:MAXOBJ)
*
*                          patch indices (edgind(i) is the index of the
*                          first edge of patch (object) nr. i)
       INTEGER    edgind(1:MAXOBJ)
*
*
*                          index of the last edge in table edges;
       INTEGER    lastedg
*
*
*                          number of faces of patches;
       INTEGER    faceno(1:MAXOBJ)
*
*                          face indexes (fceind(i) is index of the
*                          first face of object nr. i);
       INTEGER    fceind(1:MAXOBJ)
*
*
*                          index of the last face in table faces;
       INTEGER    lastfce
*
*                          list of faces which belong to a vertex;
*                          fcevrt(1,j),fcevrt(1,k1),...,fcevrt(1,ki)
*                          are indexes of faces which belong to vertex
*                          nr j; k1=fcevrt(2,j),...,ki=fcevrt(2,k(i-1))
*                          until ki=0 ;
       INTEGER    lastfv
*
*
*
*
*                          object activity (objact(i,j) is .TRUE. if
*                          object nr j is visible (activ) in view nr i
*                          or is equal .FALSE. in other case);
       LOGICAL    objact(1:MAXVIEW, 1:MAXOBJ)
*
*
*                          object presence (objprs(i) is equal .TRUE.
*                          if object nr. i is in memory (is defined) or
*                          is equal .FALSE. in other case);
       LOGICAL    objprs(1:MAXOBJ)
*
*                          begin of vertex table in storage defined
*                          by the user;
       INTEGER    vrtbeg
*
*                          begin of the edge table;
       INTEGER    edgbeg
*
*                          begin of edge presence table;
       INTEGER    prsbeg
*
*                          begin of the face table;
       INTEGER    fcebeg
*
*                          begin of face normal table;
       INTEGER    fvbeg
*
*                          begin of vertex normal table;
       INTEGER    vnbeg
       INTEGER    fnbeg
*                          maximal number of vertices;
       INTEGER    maxvrt
*
*                          maximal number of edges;
       INTEGER    maxedge
*
*                          maximal number of faces;
       INTEGER    maxface
*
*
*
       COMMON /gzptch/vertno, edgeno, faceno, lastobj,
     +                lastvrt, lastedg, lastfce,
     +                lastfv, fceind, vrtind, edgind,
     +                vrtbeg, edgbeg, prsbeg, fcebeg, fvbeg, vnbeg,
     +                fnbeg,
     +                maxvrt, maxedge, maxface,
     +                objact, objprs
      SAVE /gzptch/
*
*
*
*
*
*
*                          faces of patches; faces(1,i), faces(2,i),
*                          ...,faces(VINFACE,i) are indexes of vertexes
*                          which build a face number i;
*                          faces(VINFACE+1,i),...,faces(2*VINFACE,i)
*                          are indexes of edges which build face nr i;
*      INTEGER    faces(1:2*VINFACE, 1:MAXFACE)
*
*                          edges of patches; edges(1,i), edges(2,i)
*                          are indices of vertices which are joined by
*                          edge number i; edges(3,i), edges(4,i) are
*                          indices of faces which are divided by edge
*                          number i;
*      INTEGER    edges(1:4, 1:MAXEDGE)
*
*                          coordinates of patch vertices;
*                          vertex(1,i),vertex(2,i),vertex(3,i) = x,y,z
*                          coordinates of vertex number i;
*      REAL       vertex(1:3, 1:MAXVERT)
*
*                          list of faces which belong to a vertex;
*                          fcevrt(1,j),fcevrt(1,k1),...,fcevrt(1,ki)
*                          are indexes of faces which belong to vertex
*                          nr j; k1=fcevrt(2,j),...,ki=fcevrt(2,k(i-1))
*                          until ki=0 ;
*      INTEGER    fcevrt(1:2, 1:FOFVERT*MAXVERT)
*
*
*                          normal vectors for vertices;,
*                          vnrm(1,i),vnrm(2,i),vnrm(3,i) are x,y,z
*                          components of normal vector for vertex nr i;
*      REAL       vnrm(1:3, 1:MAXVERT)
*                          normal vectors for faces, fnrm(1,i),
*                          fnrm(2,i),fnrm(3,i) are x,y,z components
*                          of normal vector for face number i;
*      REAL       fnrm(1:3, 1:MAXFACE)
*                          presence of edges; edgprs(i) is equal.TRUE.
*                          if edge number i is visible or equal .FALSE.
*                          in other case;
*      LOGICAL    edgprs(1:MAXEDGE)
*
       REAL   vertex(1:3, *)
       INTEGER   faces(1:6, *)
       INTEGER   fcevrt(1:2, *)
       REAL      vnrm(1:3, *)
       REAL      fnrm(1:3, *)
       INTEGER   objind
*
*-----------------------------------------------------------------------
*  declare local variables
*-----------------------------------------------------------------------
       REAL       dummy, xyz(1:3, 1:VINFACE), gzunity
       INTEGER    i,j,k,m,  fi,fno,  ei,eno
*-----------------------------------------------------------------------
       fi=fceind(objind)
       fno= fi + faceno(objind) - 1
*
       DO 100 i=fi, fno
          DO 110 j=1, VINFACE
            DO 120 k=1, 3
               xyz(k, j)=vertex(k, faces(j, i))
 120        CONTINUE
 110      CONTINUE
          CALL gznormal(3 , xyz, fnrm(1, i))
          dummy=gzunity(fnrm(1, i))
 100   CONTINUE
*
*----------------------------------------   compute normals for vertexes
       ei=vrtind(objind)
       eno= ei + vertno(objind) - 1
*
       DO 130 i=ei, eno
            DO 140 m=1, 3
              vnrm(m, i)=fnrm(m, fcevrt(1, i))
 140        CONTINUE
*
           j=1
           k=fcevrt(2, i)
 160       IF(k.NE.0) THEN
              j=j+1
              DO 150 m=1, 3
                vnrm(m, i)=vnrm(m, i) + fnrm(m, fcevrt(1, k))
 150          CONTINUE
              k=fcevrt(2, k)
              GOTO 160
           ENDIF
*
           dummy=gzunity(vnrm(1, i))
 130   CONTINUE
*
       RETURN
       END
*
*----------------------------------------------------------   end gznrm
*
************************************************************************
*
C@PROCESS OPT(3) NOSDUMP NOGOSTMT
       SUBROUTINE gznormal(vertno, xyz,
     +                   norm)
*
************************************************************************
*  dummy arguments:
*               vertno - number of vertexes of a polygon;
*               xyz    - table of coordinates of polygon vertexes;
*               norm   - normal vektor for the polygon;
*
*  function:    this subroutine compute a normal vektor for the
*               polygon;
*
*
************************************************************************
       IMPLICIT NONE
*-----------------------------------------------------------------------
*  declare dummy arguments
*-----------------------------------------------------------------------
       INTEGER        vertno
       REAL           xyz(1:3, 1:3     )
       REAL           norm(1:3)
*
*-----------------------------------------------------------------------
*  declare local variables
*-----------------------------------------------------------------------
       INTEGER        m, ma, mb, i
       REAL           rightx, xv1, yv1, zv1, xv2, yv2, zv2
*
       rightx = xyz(1, 1)
       m = 1
*--------------------------   find the most right vertex of the polygon
       DO 10  i=2, vertno
          IF (xyz(1, i) .GT. rightx) THEN
              rightx = xyz(1, i)
              m = i
          ENDIF
 10    CONTINUE
*
*-----------------------   find precede vertex for the most right vertex
       ma = m-1
       IF (ma .EQ. 0) THEN
           ma = vertno
       ENDIF
*----------------------   find the next vertex for the most right vertex
       mb = m+1
       IF (mb .GT. vertno) THEN
           mb = 1
       ENDIF
*
       xv1 = xyz(1, m) - xyz(1, ma)
       yv1 = xyz(2, m) - xyz(2, ma)
       zv1 = xyz(3, m) - xyz(3, ma)
*
       xv2 = xyz(1, mb) - xyz(1, m)
       yv2 = xyz(2, mb) - xyz(2, m)
       zv2 = xyz(3, mb) - xyz(3, m)
*
       norm(1) = yv1*zv2 - zv1*yv2
       norm(2) = zv1*xv2 - xv1*zv2
       norm(3) = xv1*yv2 - yv1*xv2
*
       RETURN
       END
*---------------------------------------------------------   end normal
*
************************************************************************
*
*
       REAL      FUNCTION gzunity(norm)
*
************************************************************************
*  dummy arguments:
*               norm - normal vektor;
*
*  function:  this function returns a lenght of a vektor norm and
*             transformates the vektor to a gzunity vektor;
*
*
*
************************************************************************
       IMPLICIT NONE
*-----------------------------------------------------------------------
*  declare dummy arguments
*-----------------------------------------------------------------------
       REAL      norm(1:3)
*
*-----------------------------------------------------------------------
*  define local constants
*-----------------------------------------------------------------------
*
       REAL      ROUNDOFF
       PARAMETER (ROUNDOFF=0.999E-30)
*
*-----------------------------------------------------------------------
*  declare local variables
*-----------------------------------------------------------------------
       REAL      length
*
*-----------------------------------------------------------------------
       length = SQRT(norm(1)*norm(1) + norm(2)*norm(2) +
     +        norm(3)*norm(3))
       IF (length .LT. ROUNDOFF) THEN
           gzunity = 0.0
       ELSE
           norm(1) = norm(1)/length
           norm(2) = norm(2)/length
           norm(3) = norm(3)/length
           gzunity = length
       ENDIF
*--------------------------------------------------------   end gzunity
       END
