************************************************************************
C@PROCESS OPT(3) NOSDUMP NOGOSTMT
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
       SUBROUTINE gzinit(n, mem)
*
************************************************************************
*
*   autor:  Zenon Zowierucha, zdv069, 20.Nov.1989
*
*   file name, type:   GZINIT FORTRAN
*
*   dummy arguments:  n - number of four-byte words;
*                     mem - storage;
*
*
*   function:   this subroutine assigns default values to all global
*               variables;
*
*   called by:  main program;
*
*  CHANGED. M. BUSCH    12.4.91
*                       INTEGER*4  --> INTEGER
*                       REAL   *4  --> REAL
*  CHANGED. M. BUSCH    6.1.93 Common Bloecke nicht mehr ueber include
*                       sondern als Source eingefuegt
*
*
************************************************************************
       IMPLICIT NONE
*-----------------------------------------------------------------------
*  include global constants definitions
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
       PARAMETER (MINREAL =-1.0E30)
       PARAMETER (MAXREAL = 1.0E30)
*
       PARAMETER (ONLINE = 1)
       PARAMETER (OFFLINE= 2)
*-----------------------------------------------------------------------
*  include global variables declarations
*-----------------------------------------------------------------------
c      include(gzclrs)
***********************************************************************
*
*   autor:  Zenon Zowierucha, zdv069
*
*   file name, type :  GZCLRS COPY
*
*   function :  this file declares global variables which describe
*               patch attributes;
*               this file is destinated for use with INCLUDE
*               directive;
*
*
***********************************************************************
*
*                          front face colors; fouclr(1,i),
*                          fouclr(2,i),fouclr(3,i) are recpectivily
*                          RGB components of front face color of
*                          object number i;
       REAL       fouclr(1:3, 1:MAXOBJ)
*
*                          back face colors;
       REAL       finclr(1:3, 1:MAXOBJ)
*
*                          front edge colors;
       REAL       eouclr(1:3, 1:MAXOBJ)
*
*                          back edge color ;
       REAL       einclr(1:3, 1:MAXOBJ)
*
*
*                          front face colors low indices;
*                          in device look up table, high index is equal
*                          fouind(i)+foucno(i)-1);
       INTEGER    fouind(1:MAXOBJ)
*
*                          number of front face colors (hues);
       INTEGER    foucno(1:MAXOBJ)
*
*
*                          back face colors low indices;
*                          (high index is equal finind(i)+fincno(i)-1);
       INTEGER    finind(1:MAXOBJ)
*
*                          number of back face colors (hues);
       INTEGER    fincno(1:MAXOBJ)
*
*
*                          front edge color indices;
       INTEGER    eouind(1:MAXOBJ)
*
*                          back edge color indices;
       INTEGER    einind(1:MAXOBJ)
*
*                          front face glosses; (in range <1:MAXGLS>,
*                          1 - lowest, MAXGLS - highest gloss);
       INTEGER    fougls(1:MAXOBJ)
*
*                          back face glosses;
       INTEGER    fingls(1:MAXOBJ)
*
       INTEGER    dclro(1:MAXOBJ)
       INTEGER    dclri(1:MAXOBJ)
*
*                          look up table;
       REAL       rgb(1:3, 0:MAXCLR-1)
*
*
       COMMON/gzclrs/fouind, finind,
     +               eouind, einind, foucno, fincno,
     +               fouclr, finclr,
     +               eouclr, einclr,
     +               fougls, fingls, dclro, dclri, rgb
       SAVE /gzclrs/
*
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
*-----------------------------------------------------------------------
*  declare dummy arguments;
*-----------------------------------------------------------------------

       INTEGER   n
       INTEGER   mem(n)
*

*-----------------------------------------------------------------------
*  declare local variables; setting values
*-----------------------------------------------------------------------
*
*      ---------------------------------------   default view attributes
*                          default window
       REAL       DWINDOW(1:4)
*                          default viewport
       REAL       DVIEPRT(1:4)
*                          default projection reference point
       REAL       DPRP(1:3)
*                          default view plane distance
       REAL       DDIST
*                          default near clipping plane
       REAL       DNEAR
*                          default far clipping plane
       REAL       DFAR

*      --------------------------------------   default light attributes
*
*                          default light source position
       REAL       DSRCPOS(1:3)
*
*                          default source power
       REAL       DSRCPWR
*                          default ambient illumination
       REAL       DAMBNT
*                          default transformation matrix
       REAL       dmatrix(1:4, 1:4)
*      -----------------------------------------------------------------
*                          work variables
       INTEGER    i, j
*
       DATA       DWINDOW     /-660.0, -660.0, 660.0, 660.0/
*
       DATA       DVIEPRT     /0.0, 0.0, 1.0, 1.0 /
*
       DATA       DPRP     /0.0, 0.0,  11000.0/
*
       DATA       DDIST /512.0/
*
       DATA       DNEAR /350.0/
*
       DATA       DFAR /-580.0/
*
       DATA       DSRCPOS      /1024.0,  600.0, 1024.0/
*
       DATA       DSRCPWR /0.87/
*
       DATA       DAMBNT /0.2/
*
       DATA       DMATRIX          /1.0, 0.0, 0.0, 0.0,
     +                              0.0, 1.0, 0.0, 0.0,
     +                              0.0, 0.0, 1.0, 0.0,
     +                              0.0, 0.0, 0.0, 1.0/
*-----------------------------------------------------------------------

*                          create error file
C     CALL OBEY('FILEDEF 82 DISK GZERR GZ A ( RECFM F LRECL 80')
      OPEN(82, FILE='gzerr.gz',status='unknown')
*
       DO 100 i=1, MAXVIEW
         CALL gzdctv(i)
         CALL gzghvmpg(i,dwindow,dvieprt,PERSPECTIVE,dprp,ddist,
     +                 dnear, dfar)
         CALL gzvier(i,ON,ON, 0, 1, ON)
         CALL gzsrcv(i,dsrcpos,dsrcpwr,POINT,dambnt)
         CALL gzvtrn(i, dmatrix, REPLACE)
 100   CONTINUE
*
       lastobj=0
       lastvrt=0
       lastedg=0
       lastfce=0
*
*
       DO 110 i=1, MAXOBJ
       CALL gzsfgl(i, FRONT, 14)
       CALL gzsfgl(i, BACK, 14)
*
*
*
       vertno(i)=0
       edgeno(i)=0
       faceno(i)=0
*
       vrtind(i)=0
       edgind(i)=0
       fceind(i)=0
*
       objprs(i)=.FALSE.
*
       CALL gzmtrn(i, dmatrix, REPLACE)
       CALL gzfclr(1, 1)
*
*
       CALL gzeclr(i, FRONT, 8)
       CALL gzeclr(i, BACK,  7)
*
       DO 120 j=1,MAXVIEW
         CALL gzdptv(j,i)
 120     CONTINUE
*
       CALL gzmtrn(i, dmatrix, REPLACE)
*
 110   CONTINUE
       CALL gzactv(1)
       CALL gzaptv(1,1,PSHADED)
*
*
       maxface=n/28 - 4
       maxvrt=maxface
       maxedge= 1.4 * maxvrt
*
       vrtbeg=1
       edgbeg= vrtbeg + 3*maxvrt
       prsbeg= edgbeg + 4*maxedge
       fcebeg= prsbeg + maxedge
       fvbeg = fcebeg + 6*maxface
       vnbeg = fvbeg  + 6*maxvrt
       fnbeg = vnbeg  + 3*maxvrt
*
       lastfv=maxvrt
       CALL gzint(mem(fvbeg), maxvrt)
*
*---------------------------------------------------------   end gzinit
       END
*
*
C@PROCESS OPT(3) NOSDUMP NOGOSTMT
       SUBROUTINE gzint(fcevrt, maxvrt)
       IMPLICIT NONE
       INTEGER   fcevrt(1:2, *)
       INTEGER   maxvrt, ii
*
       DO 11 ii=1, maxvrt
          fcevrt(1,ii)=0
          fcevrt(2,ii)=0
 11    CONTINUE
       END
*
*
*
************************************************************************
C@PROCESS OPT(3) NOSDUMP NOGOSTMT
       SUBROUTINE gzfclr(mirror,icol)
*
************************************************************************
*
*  Set Face Color Palette
*
*  parameters:  objind - object index;
*               ftype  - face type (1=FRONT, 2=BACK);
*               frgb   - face color;
*               clrind - index of the first color in the palette;
*               clrno  - number of colors in the palette;
*
*  function:    this subroutine defines a color palette in workstation
*               look up table;
*
*
*  called by:   main program
*               gzinit (file GZINIT FORTRAN)
*
*
*
************************************************************************
       IMPLICIT NONE
*-----------------------------------------------------------------------
*  include gobal constants definitions and global variables deklarations
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
       PARAMETER (MINREAL =-1.0E30)
       PARAMETER (MAXREAL = 1.0E30)
*
       PARAMETER (ONLINE = 1)
       PARAMETER (OFFLINE= 2)
c      include(gzclrs)
***********************************************************************
*
*   autor:  Zenon Zowierucha, zdv069
*
*   file name, type :  GZCLRS COPY
*
*   function :  this file declares global variables which describe
*               patch attributes;
*               this file is destinated for use with INCLUDE
*               directive;
*
*
***********************************************************************
*
*                          front face colors; fouclr(1,i),
*                          fouclr(2,i),fouclr(3,i) are recpectivily
*                          RGB components of front face color of
*                          object number i;
       REAL       fouclr(1:3, 1:MAXOBJ)
*
*                          back face colors;
       REAL       finclr(1:3, 1:MAXOBJ)
*
*                          front edge colors;
       REAL       eouclr(1:3, 1:MAXOBJ)
*
*                          back edge color ;
       REAL       einclr(1:3, 1:MAXOBJ)
*
*
*                          front face colors low indices;
*                          in device look up table, high index is equal
*                          fouind(i)+foucno(i)-1);
       INTEGER    fouind(1:MAXOBJ)
*
*                          number of front face colors (hues);
       INTEGER    foucno(1:MAXOBJ)
*
*
*                          back face colors low indices;
*                          (high index is equal finind(i)+fincno(i)-1);
       INTEGER    finind(1:MAXOBJ)
*
*                          number of back face colors (hues);
       INTEGER    fincno(1:MAXOBJ)
*
*
*                          front edge color indices;
       INTEGER    eouind(1:MAXOBJ)
*
*                          back edge color indices;
       INTEGER    einind(1:MAXOBJ)
*
*                          front face glosses; (in range <1:MAXGLS>,
*                          1 - lowest, MAXGLS - highest gloss);
       INTEGER    fougls(1:MAXOBJ)
*
*                          back face glosses;
       INTEGER    fingls(1:MAXOBJ)
*
       INTEGER    dclro(1:MAXOBJ)
       INTEGER    dclri(1:MAXOBJ)
*
*                          look up table;
       REAL       rgb(1:3, 0:MAXCLR-1)
*
*
       COMMON/gzclrs/fouind, finind,
     +               eouind, einind, foucno, fincno,
     +               fouclr, finclr,
     +               eouclr, einclr,
     +               fougls, fingls, dclro, dclri, rgb
      SAVE /gzclrs/
*
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
*-----------------------------------------------------------------------
*  declare parameters
*-----------------------------------------------------------------------
*
*
       INTEGER    mirror, icol
       INTEGER    objind
*
*      INTEGER    clrind
*      INTEGER    clrno
*-----------------------------------------------------------------------
*  declare local variables
*-----------------------------------------------------------------------
*      REAL       r, pice(1:3), fou(1:3,0:MAXFCLR-1), foc(1:3)
*      INTEGER    i, j, i1
       INTEGER    i1
*
       LOGICAL    ret
*
*
*
*-----------------------------------------------------------------------
       ret=.FALSE.
*
*      IF(objind.LE.0 .OR. objind.GT.MAXOBJ) THEN
*         ret=.TRUE.
*         WRITE(DEVIND,*)'02    GZFCLR   INVALID OBJECT INDEX '
*      ENDIF
*
*      IF(ftype.NE.FRONT .AND. ftype.NE.BACK) THEN
*         ret=.TRUE.
*         WRITE(DEVIND,*)'03    GZFCLR   INVALID FACE TYPE '
*      ENDIF
*
*
*        IF(frgb(1).GT.1.0 .OR. frgb(1).LT.0.0) THEN
*            ret=.TRUE.
*            WRITE(DEVIND,*)
*    +       '04    GZFCLR   RED COLOR COMPONENT OUT OF RANGE'
*        ENDIF
*
*        IF(frgb(2).GT.1.0 .OR. frgb(2).LT.0.0) THEN
*            ret=.TRUE.
*         WRITE(DEVIND,*)
*    +    '05    GZFCLR   GREEN COLOR COMPONENT OUT OF RANGE'
*        ENDIF
*
*        IF(frgb(3).GT.1.0 .OR. frgb(3).LT.0.0) THEN
*            ret=.TRUE.
*          WRITE(DEVIND,*)
*    +     '06    GZFCLR   BLUE COLOR COMPONENT OUT OF RANGE'
*        ENDIF
*
*        IF(clrind.LT.0 .OR. clrind.GE.MAXCLR) THEN
*            ret=.TRUE.
*            WRITE(DEVIND,*)'07    GZFCLR   COLOR INDEX OUT OF RANGE'
*        ENDIF
*
*        IF(clrno.LE.3 .OR. clrno.GT.MAXFCLR) THEN
*            ret=.TRUE.
*            WRITE(DEVIND,*)
*    +       '08    GZFCLR   NUMBER OF COLORS OUT OF RANGE'
*        ENDIF
*
*        IF(clrind+clrno-1.GE.MAXCLR) THEN
*            ret=.TRUE.
*      WRITE(DEVIND,*)
*    + '09    GZFCLR   TO MANY COLORS TO DEFINE THE PALETTE'
*        ENDIF
*
*
*
*      IF(ret) THEN
*            WRITE(DEVIND,*)' GZFCLR   NO ACTION IS PERFORMED'
*            RETURN
*      ENDIF
*
*      IF(ftype.EQ.FRONT) THEN
*          dclro(objind)=1
*          r=frgb(1)
*          DO 12 i=2, 3
*            IF(frgb(i).GT.r) THEN
*              r=frgb(i)
*              dclro(objind)=i
*            ENDIF
*12        CONTINUE
*
*          DO 14 i=1, 3
*14         foc(i)=frgb(i) * 1.0/r
*      ELSE
*          dclri(objind)=1
*          r=frgb(1)
*          DO 16 i=2, 3
*            IF(frgb(i).GT.r) THEN
*               r=frgb(i)
*               dclri(objind)=i
*            ENDIF
*16        CONTINUE
*
*          DO 18 i=1, 3
*18         foc(i)=frgb(i)* 1.0/r
*
*      ENDIF
*      DO 100 i=1, 3
*100     fou(i,clrno-2)=foc(i)
*
*
*      DO 110 i=1, 3
*        ----------------------------------    set gloss color to white
*        fou(i, clrno-1) = 1.0
*        ------------------------------   compute steps of color change
*110     pice(i) = fou(i, clrno-2)/(clrno-1)
*----------------------------   set color palette for outside of a face
*      DO 120 i=clrno-3, 0, -1
*        DO 130 j=1, 3
*130       fou(j, i) = fou(j, i+1) - pice(j)
*120   CONTINUE
*
*--------------------------------------------------   set look up table
*
*
*      DO 170 i=0, clrno-1
*        DO 180 j=1, 3
*180        rgb(j, i+clrind) = fou(j, i)
*170   CONTINUE
*
*      IF(ftype.EQ.FRONT) THEN
*        DO 200 j=1, 3
* 200      fouclr(j,objind)=frgb(j)
*          fouind(objind)=clrind
*          foucno(objind)=clrno
*      ELSE
*        DO 210 j=1, 3
* 210      finclr(j,objind)=frgb(j)
*          finind(objind)=clrind
*          fincno(objind)=clrno
*      ENDIF
*
       CALL gzmirr(mirror,icol)
       DO  6789 objind=1,MAXOBJ
       
       select case(icol)

       case(1)

       dclro(objind)=1
       dclri(objind)=2
*
       case(2)

       dclro(objind)=2
       dclri(objind)=1
*
       case(3)

       dclro(objind)=3
       dclri(objind)=1
*
       case(4)


       dclro(objind)=1
       dclri(objind)=3
*
       case(5)

       dclro(objind)=2
       dclri(objind)=3
*
       case(6)

*
       dclro(objind)=3
       dclri(objind)=2
       
       end select
*
       fouind(objind)=8
       foucno(objind)=60
*
       finind(objind)=68
       fincno(objind)=60
*

       DO  121  i1=1,3
        fouclr(i1,objind)=rgb(i1,63)
        finclr(i1,objind)=rgb(i1,125)
 121   CONTINUE
 6789  CONTINUE
       END
*----------------------------------------------------------   end gzfclr
************************************************************************
C@PROCESS OPT(3) NOSDUMP NOGOSTMT
       SUBROUTINE gzeclr(objind,etype,clrind)
*
************************************************************************
*
*  Set Edge Color
*
*  parameters:  objind - object index;
*               ftype  - edge type (1=FRONT, 2=BACK);
*               frgb   - edge color;
*               clrind - index of edge color in look up table;
*
*  function:    this subroutine defines a edge color for specifed
*               object;
*
*
*  called by:   main program
*               gzinit (file GZINIT FORTRAN)
*
*
*
************************************************************************
       IMPLICIT NONE
*-----------------------------------------------------------------------
*  include gobal constants definitions and global variables deklarations
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
       PARAMETER (MINREAL =-1.0E30)
       PARAMETER (MAXREAL = 1.0E30)
*
       PARAMETER (ONLINE = 1)
       PARAMETER (OFFLINE= 2)
c      include(gzclrs)
***********************************************************************
*
*   autor:  Zenon Zowierucha, zdv069
*
*   file name, type :  GZCLRS COPY
*
*   function :  this file declares global variables which describe
*               patch attributes;
*               this file is destinated for use with INCLUDE
*               directive;
*
*
***********************************************************************
*
*                          front face colors; fouclr(1,i),
*                          fouclr(2,i),fouclr(3,i) are recpectivily
*                          RGB components of front face color of
*                          object number i;
       REAL       fouclr(1:3, 1:MAXOBJ)
*
*                          back face colors;
       REAL       finclr(1:3, 1:MAXOBJ)
*
*                          front edge colors;
       REAL       eouclr(1:3, 1:MAXOBJ)
*
*                          back edge color ;
       REAL       einclr(1:3, 1:MAXOBJ)
*
*
*                          front face colors low indices;
*                          in device look up table, high index is equal
*                          fouind(i)+foucno(i)-1);
       INTEGER    fouind(1:MAXOBJ)
*
*                          number of front face colors (hues);
       INTEGER    foucno(1:MAXOBJ)
*
*
*                          back face colors low indices;
*                          (high index is equal finind(i)+fincno(i)-1);
       INTEGER    finind(1:MAXOBJ)
*
*                          number of back face colors (hues);
       INTEGER    fincno(1:MAXOBJ)
*
*
*                          front edge color indices;
       INTEGER    eouind(1:MAXOBJ)
*
*                          back edge color indices;
       INTEGER    einind(1:MAXOBJ)
*
*                          front face glosses; (in range <1:MAXGLS>,
*                          1 - lowest, MAXGLS - highest gloss);
       INTEGER    fougls(1:MAXOBJ)
*
*                          back face glosses;
       INTEGER    fingls(1:MAXOBJ)
*
       INTEGER    dclro(1:MAXOBJ)
       INTEGER    dclri(1:MAXOBJ)
*
*                          look up table;
       REAL       rgb(1:3, 0:MAXCLR-1)
*
*
       COMMON/gzclrs/fouind, finind,
     +               eouind, einind, foucno, fincno,
     +               fouclr, finclr,
     +               eouclr, einclr,
     +               fougls, fingls, dclro, dclri, rgb
       SAVE /gzclrs/
*
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
*-----------------------------------------------------------------------
*  declare parameters
*-----------------------------------------------------------------------
*
       INTEGER    objind
       INTEGER    etype
*
       INTEGER    clrind
*-----------------------------------------------------------------------
*  declare local variables
*-----------------------------------------------------------------------
*
       LOGICAL    ret
*
*
*
*-----------------------------------------------------------------------
       ret=.FALSE.
*
       IF(objind.LE.0 .OR. objind.GT.MAXOBJ) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'02    GZECLR   INVALID OBJECT INDEX '
       ENDIF
*
       IF(etype.NE.FRONT .AND. etype.NE.BACK) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'10    GZECLR   INVALID EDGE TYPE '
       ENDIF
*
*
*        IF(ergb(1).GT.1.0 .OR. ergb(1).LT.0.0) THEN
*            ret=.TRUE.
*            WRITE(DEVIND,*)
*    +       '04    GZECLR   RED COLOR COMPONENT OUT OF RANGE'
*        ENDIF
*
*        IF(ergb(2).GT.1.0 .OR. ergb(2).LT.0.0) THEN
*            ret=.TRUE.
*         WRITE(DEVIND,*)
*    +    '05    GZECLR   GREEN COLOR COMPONENT OUT OF RANGE'
*        ENDIF
*
*        IF(ergb(3).GT.1.0 .OR. ergb(3).LT.0.0) THEN
*            ret=.TRUE.
*          WRITE(DEVIND,*)
*    +     '06    GZECLR   BLUE COLOR COMPONENT OUT OF RANGE'
*        ENDIF
*
         IF(clrind.LT.0 .OR. clrind.GE.MAXCLR) THEN
             ret=.TRUE.
             WRITE(DEVIND,*)'07    GZECLR   COLOR INDEX OUT OF RANGE'
         ENDIF
*
*
*
*
*
       IF(ret) THEN
             WRITE(DEVIND,*)' GZECLR   NO ACTION IS PERFORMED'
             RETURN
       ENDIF
*
*
       IF(etype.EQ.FRONT) THEN
*        DO 200 j=1, 3
*          rgb(j,clrind)=ergb(j)
* 200      eouclr(j,objind)=ergb(j)
           eouind(objind)=clrind
       ELSE
*        DO 210 j=1, 3
*          rgb(j,clrind)=ergb(j)
*          einclr(j,objind)=ergb(j)
           einind(objind)=clrind
       ENDIF
*

       END
*----------------------------------------------------------   end gzeclr
************************************************************************
C@PROCESS OPT(3) NOSDUMP NOGOSTMT
       SUBROUTINE gzvier(view, anearfl, afarfl, bcind,
     +                   brind, brdflg)
*
************************************************************************
*
*  Set View Characteristic
*
*  dummy arguments:  view - index of view;
*                    anearfl - near clipping flag (OFF or ON);
*                    afarfl  - far clipping flag (OFF or ON);
*                    bcgrgb  - background color;
*                    bcind   - background color index;
*                    brdrgb  - border color;
*                    brind   - border color index;
*                    brdflg  - border flag (OFF or ON);
*
*  function:  this subroutine set a view characteristic;
*
*  called by:  main program
*              gzinit (file GZINIT FORTRAN)
*
************************************************************************
       IMPLICIT NONE
*-----------------------------------------------------------------------
*  include global constants and variables
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
       PARAMETER (MINREAL =-1.0E30)
       PARAMETER (MAXREAL = 1.0E30)
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
       SAVE /gzvchr/
*
c      include(gzclrs)
***********************************************************************
*
*   autor:  Zenon Zowierucha, zdv069
*
*   file name, type :  GZCLRS COPY
*
*   function :  this file declares global variables which describe
*               patch attributes;
*               this file is destinated for use with INCLUDE
*               directive;
*
*
***********************************************************************
*
*                          front face colors; fouclr(1,i),
*                          fouclr(2,i),fouclr(3,i) are recpectivily
*                          RGB components of front face color of
*                          object number i;
       REAL       fouclr(1:3, 1:MAXOBJ)
*
*                          back face colors;
       REAL       finclr(1:3, 1:MAXOBJ)
*
*                          front edge colors;
       REAL       eouclr(1:3, 1:MAXOBJ)
*
*                          back edge color ;
       REAL       einclr(1:3, 1:MAXOBJ)
*
*
*                          front face colors low indices;
*                          in device look up table, high index is equal
*                          fouind(i)+foucno(i)-1);
       INTEGER    fouind(1:MAXOBJ)
*
*                          number of front face colors (hues);
       INTEGER    foucno(1:MAXOBJ)
*
*
*                          back face colors low indices;
*                          (high index is equal finind(i)+fincno(i)-1);
       INTEGER    finind(1:MAXOBJ)
*
*                          number of back face colors (hues);
       INTEGER    fincno(1:MAXOBJ)
*
*
*                          front edge color indices;
       INTEGER    eouind(1:MAXOBJ)
*
*                          back edge color indices;
       INTEGER    einind(1:MAXOBJ)
*
*                          front face glosses; (in range <1:MAXGLS>,
*                          1 - lowest, MAXGLS - highest gloss);
       INTEGER    fougls(1:MAXOBJ)
*
*                          back face glosses;
       INTEGER    fingls(1:MAXOBJ)
*
       INTEGER    dclro(1:MAXOBJ)
       INTEGER    dclri(1:MAXOBJ)
*
*                          look up table;
       REAL       rgb(1:3, 0:MAXCLR-1)
*
*
       COMMON/gzclrs/fouind, finind,
     +               eouind, einind, foucno, fincno,
     +               fouclr, finclr,
     +               eouclr, einclr,
     +               fougls, fingls, dclro, dclri, rgb
       SAVE /gzclrs/
*
*-----------------------------------------------------------------------
*  declare dummy arguments
*-----------------------------------------------------------------------
       INTEGER    view, anearfl, afarfl
*
       INTEGER    bcind
*
       INTEGER    brind
       INTEGER    brdflg
*-----------------------------------------------------------------------
       LOGICAL    ret
*-----------------------------------------------------------------------
*
       ret=.FALSE.
*
       IF(view.LE.0 .OR. view.GT.MAXVIEW) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'01    GZVIER   INVALID VIEW INDEX '
       ENDIF
*
       IF(anearfl.NE.OFF .AND. anearfl.NE.ON) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'11    GZVIER   INVALID NEAR FLAG VALUE'
       ENDIF
       IF(afarfl.NE.OFF .AND. afarfl.NE.ON) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'12    GZVIER   INVALID FAR FLAG VALUE'
       ENDIF
*
*
*
         IF(bcind.LT.0 .OR. bcind.GE.MAXCLR) THEN
             ret=.TRUE.
             WRITE(DEVIND,*)
     +       '16    GZVIER   BACKGROUND COLOR INDEX OUT RANGE'
         ENDIF
*
*
         IF(brind.LT.0 .OR. brind.GE.MAXCLR) THEN
             ret=.TRUE.
         WRITE(DEVIND,*)
     +   '20    GZVIER   BORDER COLOR INDEX OUT OF RANGE'
         ENDIF
*
       IF(brdflg.NE.OFF .AND. brdflg.NE.ON) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'21    GZVIER   INVALID BORDER FLAG VALUE'
       ENDIF
*
*
*
*
       IF(ret) THEN
             WRITE(DEVIND,*)' GZVIER   NO ACTION IS PERFORMED'
             RETURN
       ENDIF
*
*
       nearfl(view) = anearfl
       farfl(view) = afarfl
*
       bcgind(view)=bcind
       brdind(view)=brind
       brdact(view)=brdflg
*
       RETURN
       END
*----------------------------------------------------------   end ghvier
*
************************************************************************
C@PROCESS OPT(3) NOSDUMP NOGOSTMT DC(feld)
       SUBROUTINE gzghvmpg(view, awindow, avieprt, aprtype, aprp, adist,
     +                     anear, afar)
************************************************************************
*
*  Set View Mapping
*
*  dummy arguments: view - view index;
*                   awindow - view window;
*                   avieprt - viewport;
*                   aprtype - projection type;
*                   aprp    - projection reference point;
*                   adist   - view plane distance;
*                   anear   - near clipping plane distance;
*                   afar    - far clipping plane distance;
*
*  function:  this subroutine set global view mapping values;
*
*  called by: main program
*             gzinit (file GZINIT FORTRAN)
*
************************************************************************
       IMPLICIT NONE
*-----------------------------------------------------------------------
*  include global constants and variables
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
       PARAMETER (MINREAL =-1.0E30)
       PARAMETER (MAXREAL = 1.0E30)
*
       PARAMETER (ONLINE = 1)
       PARAMETER (OFFLINE= 2)
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
       SAVE /gzvmpg/
c      include(gzclrs)
***********************************************************************
*
*   autor:  Zenon Zowierucha, zdv069
*
*   file name, type :  GZCLRS COPY
*
*   function :  this file declares global variables which describe
*               patch attributes;
*               this file is destinated for use with INCLUDE
*               directive;
*
*
***********************************************************************
*
*                          front face colors; fouclr(1,i),
*                          fouclr(2,i),fouclr(3,i) are recpectivily
*                          RGB components of front face color of
*                          object number i;
       REAL       fouclr(1:3, 1:MAXOBJ)
*
*                          back face colors;
       REAL       finclr(1:3, 1:MAXOBJ)
*
*                          front edge colors;
       REAL       eouclr(1:3, 1:MAXOBJ)
*
*                          back edge color ;
       REAL       einclr(1:3, 1:MAXOBJ)
*
*
*                          front face colors low indices;
*                          in device look up table, high index is equal
*                          fouind(i)+foucno(i)-1);
       INTEGER    fouind(1:MAXOBJ)
*
*                          number of front face colors (hues);
       INTEGER    foucno(1:MAXOBJ)
*
*
*                          back face colors low indices;
*                          (high index is equal finind(i)+fincno(i)-1);
       INTEGER    finind(1:MAXOBJ)
*
*                          number of back face colors (hues);
       INTEGER    fincno(1:MAXOBJ)
*
*
*                          front edge color indices;
       INTEGER    eouind(1:MAXOBJ)
*
*                          back edge color indices;
       INTEGER    einind(1:MAXOBJ)
*
*                          front face glosses; (in range <1:MAXGLS>,
*                          1 - lowest, MAXGLS - highest gloss);
       INTEGER    fougls(1:MAXOBJ)
*
*                          back face glosses;
       INTEGER    fingls(1:MAXOBJ)
*
       INTEGER    dclro(1:MAXOBJ)
       INTEGER    dclri(1:MAXOBJ)
*
*                          look up table;
       REAL       rgb(1:3, 0:MAXCLR-1)
*
*
       COMMON/gzclrs/fouind, finind,
     +               eouind, einind, foucno, fincno,
     +               fouclr, finclr,
     +               eouclr, einclr,
     +               fougls, fingls, dclro, dclri, rgb
       SAVE /gzclrs/
*
c      include(gzzbuf)
***********************************************************************
*
*  autor:  Zenon Zowierucha, zdv069
*
*  file name, type:  GZZBUF COPY
*
*  function:  this file declares global variables which describe
*             z-buffer table and picture table;
*             this file is destinated for use with INCLUDE directive;
*
***********************************************************************
*                          z-buffer; (first index is x coordinate,
*                          second -  y coordinate);
       REAL       zbuff(0:MAXPKX-1, 0:MAXPKY-1)
*
*                          frame buffer (picture);
       CHARACTER (len=1) picture(0:MAXPKX-1, 0:MAXPKY-1)
*
       COMMON/gzzbuff/zbuff
       COMMON/feld/picture
       SAVE /gzzbuff/,/feld/
*
*
c      include(rdoff)
       REAL     ROUNDOFF
       PARAMETER (ROUNDOFF=1.0E-10)
*-----------------------------------------------------------------------
*  declare dummy arguments
*-----------------------------------------------------------------------
       INTEGER    view
       REAL       awindow(1:4), avieprt(1:4)
       INTEGER    aprtype
       REAL       aprp(1:3), adist, anear, afar
*-----------------------------------------------------------------------
*  declare local variables
*-----------------------------------------------------------------------
       INTEGER    i, j
       REAL       p(1:3, 1:4)
       REAL       nrm(1:3)
       REAL       dummy, gzunity
       LOGICAL    ret
*
*-----------------------------------------------------------------------
       ret=.FALSE.
*
       IF(view.LT.1 .OR. view.GT.MAXVIEW) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'01    GZGHVMPG  INVALID VIEW INDEX'
       ENDIF
*
       IF(awindow(1).GE.awindow(3) .OR. awindow(2).GE.awindow(4)) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'22    GZGHVMPG  INVALID WINDOW'
       ENDIF
*
*
       IF(avieprt(1).GE.avieprt(3) .OR. avieprt(2).GE.avieprt(4) .OR.
     +    avieprt(1).LT.0.0 .OR. avieprt(2).LT.0.0 .OR.
     +    avieprt(3).LT.0.0 .OR. avieprt(4).LT.0.0 .OR.
     +    avieprt(1).GT.1.0 .OR. avieprt(2).GT.1.0 .OR.
     +    avieprt(3).GT.1.0 .OR. avieprt(4).GT.1.0) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'23    GZGHVMPG  INVALID VIEWPORT'
       ENDIF
*
       IF(aprtype.NE.PARALLEL .AND. aprtype.NE.PERSPECTIVE) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)
     +    '24    GZGHVMPG  INVALID PROJECTION TYPE VALUE'
       ENDIF
*
       IF(ABS(aprp(3)-adist).LE.ROUNDOFF) THEN
          ret=.TRUE.
       WRITE(DEVIND,*)
     + '25    GZGHVMPG  PRP IS POSITIONED ON THE VIEW PLANE'
       ENDIF
*
*
       IF(afar.GE.anear) THEN
          ret=.TRUE.
       WRITE(DEVIND,*)
     +'26 GZGHVMPG  FAR CLIPPING PLANE IS FRONT OF NEAR CLIPPING PLANE'
       ENDIF
*
       IF(ret) THEN
          WRITE(DEVIND,*)' GZGHVMPG  NO ACTION IS PERFORMED'
          RETURN
       ENDIF
*
           vieprt(1, view) = avieprt(1)*(MAXPKX-1.0)
           vieprt(2, view) = avieprt(2)*(MAXPKY-1.0)
           vieprt(3, view) = avieprt(3)*(MAXPKX-1.0)
           vieprt(4, view) = avieprt(4)*(MAXPKY-1.0)
*
           window(1, view) = awindow(1)
           window(2, view) = awindow(2)
           window(3, view) = awindow(3)
           window(4, view) = awindow(4)
*
*
       prtype(view) = aprtype
*
       DO 120 i=1, 3
          prp(i, view) = aprp(i)
 120   CONTINUE
*
       dist(view) = adist
       near(view) = anear
       far(view) = afar
*
*---------------------   compute parameters of clipping planes equations
*                        up right vertex of the window has index 1;
*                        up left- 2; down left- 3; down right- 4;
       p(1,1)= window(3,view)
       p(2,1)= window(4,view)
*                          up left window vertex
       p(1,2)= window(1,view)
       p(2,2)= window(4,view)
*                          down left
       p(1,3)= window(1,view)
       p(2,3)= window(2,view)
*                          down right
       p(1,4)= window(3,view)
       p(2,4)= window(2,view)
*                          all window vertices have the same z values
       DO 130 i=1, 4
         p(3,i)=adist
 130   CONTINUE
*
       DO 140 i=1, 4
         DO 150 j=1, 3
           p(j,i)=p(j,i)-aprp(j)
 150     CONTINUE
 140   CONTINUE
*
       DO 160 i=1, 4
       j=i+1
       IF(i.EQ.4) j=1
       nrm(1)= p(3,i)*p(2,j) - p(2,i)*p(3,j)
       nrm(2)= p(1,i)*p(3,j) - p(3,i)*p(1,j)
       nrm(3)= p(2,i)*p(1,j) - p(1,i)*p(2,j)
*
       dummy=gzunity(nrm)
*
       abcd(1,i,view)= nrm(1)
       abcd(2,i,view)= nrm(2)
       abcd(3,i,view)= nrm(3)
       abcd(4,i,view)= -(nrm(1)*aprp(1) + nrm(2)*aprp(2) +
     +                   nrm(3)*aprp(3))
 160   CONTINUE
*
       END
*---------------------------------------------------------   end gzghvmp
*
************************************************************************
C@PROCESS OPT(3) NOSDUMP NOGOSTMT DC(feld)
       SUBROUTINE gzvtrn(view, matrix, type)
*
************************************************************************
*
*  Set View Transformation
*
*  dummy arguments: view   - view index;
*                   matrix - view transformation matrix;
*                   type   - composition type ( 1 - PRECONCATENATE,
*                            2 - POSTCONCATENATE, 3 - REPLACE);
*
*  function:  this subroutine sets global view tranformation;
*
*  called by: main program
*             gzinit (file GZINIT FORTRAN)
*
*  calls:     gzmult (file GZINIT FORTRAN)
*
************************************************************************
       IMPLICIT NONE
*-----------------------------------------------------------------------
*  include global constants and variables
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
       PARAMETER (MINREAL =-1.0E30)
       PARAMETER (MAXREAL = 1.0E30)
*
       PARAMETER (ONLINE = 1)
       PARAMETER (OFFLINE= 2)
c      include(gztrns)
***********************************************************************
*
*  autor: Zenon Zowierucha, zdv069
*
*  file name, type:   GZTRNS COPY
*
*  function:  this file declares global variables that define three
*             dimmensional transformations;
*             this file is destinated for a use with INCLUDE directive;
*
*
***********************************************************************
*
*                          view transformation matrix
       REAL       vtrans(1:4, 1:4, 1:MAXVIEW)
*
*                          modeling transformation matrix
*                          (mtrans(1,1,i) is modeling transformation
*                           matrix for patch (object) nr. i)
       REAL       mtrans(1:4, 1:4, 1:MAXOBJ)
*
*
       COMMON/gztrns/vtrans, mtrans
       SAVE /gztrns/
*
*-----------------------------------------------------------------------
*  declare dummy arguments
*-----------------------------------------------------------------------
       INTEGER    view
       REAL       matrix(1:4, 1:4)
       INTEGER    type
*-----------------------------------------------------------------------
*  declare local variables
*-----------------------------------------------------------------------
       INTEGER    i, j
       LOGICAL    ret
*-----------------------------------------------------------------------
       ret=.FALSE.
*
       IF(view.LT.1 .OR. view.GT.MAXVIEW) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'01    GZVTRN  INVALID VIEW INDEX'
       ENDIF
*
*------------------------------------   verify if matrix is zero matrix
       i=1
       j=1
 50    IF(matrix(j, i).NE.0.0) GOTO 80
       j=j+1
       IF(j.EQ.4) i=i+1
       IF(i.EQ.5) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'27    GZVTRN  ZERO TRANFORMATION MATRIX'
          GOTO 80
       ENDIF
*
       GOTO 50
*
 80    IF(type.NE.PRECONCATENATE .AND. type.NE.POSTCONCATENATE .AND.
     +    type.NE.REPLACE) THEN
              ret=.TRUE.
           WRITE(DEVIND,*)'28    GZVTRN  INVALID COMPOSITION TYPE'
       ENDIF
*
       IF(ret) THEN
         WRITE(DEVIND,*)' GZVTRN  NO ACTION IS PERFORMED'
         RETURN
       ENDIF
*
*-----------------------------------------------------------------------
      
       select case(type)

      
*      -----------------------------------------------   preconcatenate
       case(1)
       CALL gzmult(vtrans(1,1,view), matrix, vtrans(1,1,view))
      
*      ----------------------------------------------   postconcatenate
       case(2)
       CALL gzmult(matrix, vtrans(1,1,view), vtrans(1,1,view))
      
*      ------------------------------------------------------   replace
       case(3)
       DO 310 i=1, 4
         DO 320 j=1, 4
           vtrans(j, i, view) = matrix(j, i)
 320     CONTINUE
 310   CONTINUE

       case default
      
       end select 

       RETURN
       END
*---------------------------------------------------------   end gzvtrn
*
************************************************************************
C@PROCESS OPT(3) NOSDUMP NOGOSTMT DC(feld)
       SUBROUTINE gzmtrn(objind, matrix, type)
*
************************************************************************
*
*  Set Modeling Transformation
*
*  dummy arguments: objind - object index;
*                   matrix - view transformation matrix;
*                   type   - composition type ( 1 - PRECONCATENATE,
*                            2 - POSTCONCATENATE, 3 - REPLACE);
*
*  function:  this subroutine sets global modeling tranformation
*             for specifed object;
*
*  called by: main program
*             gzinit (file GZINIT FORTRAN)
*
*  calls:     gzmult (file GZINIT FORTRAN)
*
************************************************************************
*-----------------------------------------------------------------------
*  include global constants and variables
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
       PARAMETER (MINREAL =-1.0E30)
       PARAMETER (MAXREAL = 1.0E30)
*
       PARAMETER (ONLINE = 1)
       PARAMETER (OFFLINE= 2)
c      include(gztrns)
***********************************************************************
*
*  autor: Zenon Zowierucha, zdv069
*
*  file name, type:   GZTRNS COPY
*
*  function:  this file declares global variables that define three
*             dimmensional transformations;
*             this file is destinated for a use with INCLUDE directive;
*
*
***********************************************************************
*
*                          view transformation matrix
       REAL       vtrans(1:4, 1:4, 1:MAXVIEW)
*
*                          modeling transformation matrix
*                          (mtrans(1,1,i) is modeling transformation
*                           matrix for patch (object) nr. i)
       REAL       mtrans(1:4, 1:4, 1:MAXOBJ)
*
*
       COMMON/gztrns/vtrans, mtrans
       SAVE /gztrns/
*
*-----------------------------------------------------------------------
*  declare dummy arguments
*-----------------------------------------------------------------------
       INTEGER    objind
       REAL       matrix(1:4, 1:4)
       INTEGER    type
*-----------------------------------------------------------------------
*  declare local variables
*-----------------------------------------------------------------------
       INTEGER    i, j
       LOGICAL    ret
*-----------------------------------------------------------------------
       ret=.FALSE.
*
*
       IF(objind.LE.0 .OR. objind.GT.MAXOBJ) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'02    GZMTRN   INVALID OBJECT INDEX '
       ENDIF
*
*------------------------------------   verify if matrix is zero matrix
       i=1
       j=1
 50    IF(matrix(j, i).NE.0.0) GOTO 80
       j=j+1
       IF(j.EQ.4) i=i+1
       IF(i.EQ.5) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'27    GZMTRN  ZERO TRANFORMATION MATRIX'
          GOTO 80
       ENDIF
*
       GOTO 50
*
 80    IF(type.NE.PRECONCATENATE .AND. type.NE.POSTCONCATENATE .AND.
     +    type.NE.REPLACE) THEN
              ret=.TRUE.
           WRITE(DEVIND,*)'28    GZMTRN  INVALID COMPOSITION TYPE'
       ENDIF
*
       IF(ret) THEN
         WRITE(DEVIND,*)' GZMTRN  NO ACTION IS PERFORMED'
         RETURN
       ENDIF
*
*-----------------------------------------------------------------------
     
       select case(type)
       
*      -----------------------------------------------   preconcatenate
       case(1)

      CALL gzmult(mtrans(1,1,objind), matrix, mtrans(1,1,objind))
     
*      ----------------------------------------------   postconcatenate
       case(2)

      CALL gzmult(matrix, mtrans(1,1,objind), mtrans(1,1,objind))
    
*      ------------------------------------------------------   replace
       case(3)

       DO 310 i=1, 4
         DO 320 j=1, 4
           mtrans(j, i, objind) = matrix(j, i)
 320     CONTINUE      
 310   CONTINUE
       case default
       end select
      
       RETURN
       END
*---------------------------------------------------------   end gzmtrn
*
************************************************************************
C@PROCESS OPT(3) NOSDUMP NOGOSTMT DC(feld)
       SUBROUTINE gzmult(mtr1, mtr2, mtr3)
*
************************************************************************
*
*  dummy arguments: mtr1, mtr2, mtr3 - tranformation matrices;
*
*  function:  this subroutine gzmultiplys matrixes mtr1 and mtr2 and set
*             a result in matrix mtr3;
*
*  called by:  gzvtrn (file GZINIT FORTRAN)
*
*              gzmtrn (file GZINIT FORTRAN)
*
************************************************************************
       IMPLICIT NONE
*-----------------------------------------------------------------------
*  declare dummy arguments
*-----------------------------------------------------------------------
       REAL       mtr1(1:4,1:4), mtr2(1:4,1:4), mtr3(1:4,1:4)
*-----------------------------------------------------------------------
*  declare local variables
*-----------------------------------------------------------------------
       REAL       mtr4(1:4,1:4), temp
       INTEGER    i, j, k
*-----------------------------------------------------------------------
       DO 100 i=1, 4
         DO 110 j=1, 4
            temp = 0.0
              DO 120 k=1, 4
                 temp = temp + mtr1(k, i)*mtr2(j, k)
 120          CONTINUE
            mtr4(j, i) = temp
 110     CONTINUE
 100   CONTINUE
*
       DO 130 i=1, 4
         DO 140 j=1, 4
           mtr3(j, i) = mtr4(j, i)
 140     CONTINUE
 130   CONTINUE
*
       RETURN
       END
*-----------------------------------------------------------   end gzmul
*
************************************************************************
C@PROCESS OPT(3) NOSDUMP NOGOSTMT DC(feld)
       SUBROUTINE gzsrcv(view, asrcpos, asrcpwr, asrctyp, aambnt)
*
************************************************************************
*
*  Associate Source With View
*
*  dummy arguments: view - view index:
*                   asrcpos - position of a light source;
*                   asrcpwr - power of a light source;
*                   asrctyp - type of a light source ( SUN or POINT);
*                   aambnt  - ambient illumination;
*
*  function:  this subroutine sets light characteristics for a view;
*
*  called by: main program
*             gzinit (file GZINIT FORTRAN)
*
************************************************************************
       IMPLICIT NONE
*-----------------------------------------------------------------------
*  include global constants and variables
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
       PARAMETER (MINREAL =-1.0E30)
       PARAMETER (MAXREAL = 1.0E30)
*
       PARAMETER (ONLINE = 1)
       PARAMETER (OFFLINE= 2)
c      include(gzsrc)
***********************************************************************
*
*   autor: Zenon Zowierucha, zdv069
*
*   file name, type :   GZSRC COPY
*
*   function :   this file declares global variables which describe
*                light sources and illumination characteristics;
*                this file is destinated for use with INCLUDE
*                directive;
*
*
***********************************************************************
*
*                          light source position (VCS)
       REAL       srcpos(1:3, 1:MAXVIEW)
*
*                          light source power (in range <0:1>,
*                          1 - highest, 0 - lowest power)
       REAL       srcpwr(1:MAXVIEW)
*
*                          light source type (SUN or POINT)
       INTEGER    srctyp(1:MAXVIEW)
*
*                          ambient light intensity (in range <0:1>,
*                          1 - highest, 0 - lowest intensity)
       REAL       ambnt(1:MAXVIEW)
*
*
*
       COMMON/gzsrc/srctyp, srcpos, srcpwr, ambnt
       SAVE /gzsrc/
*
*-----------------------------------------------------------------------
*  declare dummy arguments
*-----------------------------------------------------------------------
       INTEGER    view
       REAL       asrcpos(1:3), asrcpwr
       INTEGER    asrctyp
       REAL       aambnt
*-----------------------------------------------------------------------
*  declare local variables
*-----------------------------------------------------------------------
       INTEGER    i
       REAL       dummy, gzunity
       LOGICAL    ret
*-----------------------------------------------------------------------
       ret=.FALSE.
*
*
       IF(view.LT.1 .OR. view.GT.MAXVIEW) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'01    GZSRCV  INVALID VIEW INDEX'
       ENDIF
*
*
       IF(asrcpwr.LT.0.0 .OR. asrcpwr.GT.1.0) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'29    GZSRCV  INVALID SOURCE POWER VALUE'
       ENDIF
*
*
       IF(asrctyp.NE.SUN .AND. asrctyp.NE.POINT) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'30    GZSRCV  INVALID SOURCE TYPE'
       ENDIF
*
*
       IF(aambnt.LT.0.0 .OR. aambnt.GT.1.0) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)
     +    '31    GZSRCV  INVALID AMBIENT ILLUMINATION VALUE'
       ENDIF
*
*
       IF(ret) THEN
         WRITE(DEVIND,*)' GZSRCV  NO ACTION IS PERFORMED'
         RETURN
       ENDIF
*
*
       DO 100 i=1, 3
         srcpos(i, view) = asrcpos(i)
 100   CONTINUE
       IF(srctyp(view).EQ.SUN) dummy=gzunity(srcpos(1,view))
*
       srctyp(view) = asrctyp
       srcpwr(view) = asrcpwr
       ambnt(view) = aambnt
*
       RETURN
       END
*---------------------------------------------------------  end gzsrcv
*
************************************************************************
C@PROCESS OPT(3) NOSDUMP NOGOSTMT DC(feld)
       SUBROUTINE gzaptv(view, objind, alg)
*
************************************************************************
*
*  Associate Patch With View
*
*  dummy arguments:
*                   view  - view index;
*                   objind- object index;
*                   alg   - view algorithm;
*
*  function:  this subroutine acctivates specifed view and object;
*
*  called by: main program
*             gzinit (file GZINIT FORTRAN)
*
************************************************************************
       IMPLICIT NONE
*-----------------------------------------------------------------------
*  iclude global constants and variables
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
       PARAMETER (MINREAL =-1.0E30)
       PARAMETER (MAXREAL = 1.0E30)
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
       SAVE /gzvchr/
*
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
      SAVE  /gzptch/
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
*-----------------------------------------------------------------------
*  declare dummy arguments
*-----------------------------------------------------------------------
       INTEGER    view
       INTEGER    objind
       INTEGER    alg
*-----------------------------------------------------------------------
       LOGICAL    ret
*-----------------------------------------------------------------------
       ret=.FALSE.
*
*
       IF(view.LT.1 .OR. view.GT.MAXVIEW) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'01    GZAPTV  INVALID VIEW INDEX'
       ENDIF
*
       IF(objind.LE.0 .OR. objind.GT.MAXOBJ) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'02    GZAPTV   INVALID OBJECT INDEX '
       ENDIF
*
       IF(alg.LT.WIRE .OR. alg.GT.PSHADED) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'32    GZAPTV   INVALID ALGORITHM TYPE'
       ENDIF
*
       IF(ret) THEN
         WRITE(DEVIND,*)' GZAPTV  NO ACTION IS PERFORMED'
         RETURN
       ENDIF
*
           vstyle(view, objind) = alg
*
           objact(view, objind)=.TRUE.
*
*
       RETURN
       END
*---------------------------------------------------------   end gzaptv
*
************************************************************************
C@PROCESS OPT(3) NOSDUMP NOGOSTMT DC(feld)
       SUBROUTINE gzdptv(view, objind)
*
************************************************************************
*
*  Dissassociate Patch With View
*
*  dummy arguments: view  - view index;
*                   objind- object index
*
*  function:  this subroutine dissactivates specifed object in view;
*
*  called by: main program
*
************************************************************************
       IMPLICIT NONE
*-----------------------------------------------------------------------
*  iclude global constants and variables
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
       PARAMETER (MINREAL =-1.0E30)
       PARAMETER (MAXREAL = 1.0E30)
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
       SAVE /gzvchr/
*
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
*-----------------------------------------------------------------------
*  declare dummy arguments
*-----------------------------------------------------------------------
       INTEGER    view
       INTEGER    objind
*-----------------------------------------------------------------------
       LOGICAL    ret
*-----------------------------------------------------------------------
       ret=.FALSE.
*
       IF(view.LT.1 .OR. view.GT.MAXVIEW) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'01    GZDPTV  INVALID VIEW INDEX'
       ENDIF
*
       IF(objind.LE.0 .OR. objind.GT.MAXOBJ) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'02    GZDPTV   INVALID OBJECT INDEX '
       ENDIF
*
*
       IF(ret) THEN
         WRITE(DEVIND,*)' GZDPTV  NO ACTION IS PERFORMED'
         RETURN
       ENDIF
*
*
*
       objact(view,objind)=.FALSE.
       RETURN
*---------------------------------------------------------   end gzdptv
       END
*
*
************************************************************************
C@PROCESS OPT(3) NOSDUMP NOGOSTMT
       SUBROUTINE gzactv(view)
*
************************************************************************
*
*  Activate a View
*
*  dummy arguments: view  - view index;
*
*  function:  this subroutine activates specifed view;
*
*  called by: main program
*
************************************************************************
       IMPLICIT NONE
*-----------------------------------------------------------------------
*  iclude global constants and variables
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
       PARAMETER (MINREAL =-1.0E30)
       PARAMETER (MAXREAL = 1.0E30)
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
       SAVE /gzvchr/
*
*-----------------------------------------------------------------------
*  declare dummy arguments
*-----------------------------------------------------------------------
       INTEGER    view
*-----------------------------------------------------------------------
       LOGICAL    ret
*-----------------------------------------------------------------------
       ret=.FALSE.
*
       IF(view.LT.1 .OR. view.GT.MAXVIEW) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'01    GZACTV  INVALID VIEW INDEX'
       ENDIF
*
*
*
       IF(ret) THEN
         WRITE(DEVIND,*)' GZACTV  NO ACTION IS PERFORMED'
         RETURN
       ENDIF
*
*
*
       vactiv(view)=.TRUE.
       RETURN
*---------------------------------------------------------   end gzactv
       END
*
*
************************************************************************
C@PROCESS OPT(3) NOSDUMP NOGOSTMT
       SUBROUTINE gzdctv(view)
*
************************************************************************
*
*  Dissactivate a View
*
*  dummy arguments: view  - view index;
*
*  function:  this subroutine activates specifed view;
*
*  called by: main program
*
************************************************************************
       IMPLICIT NONE
*-----------------------------------------------------------------------
*  iclude global constants and variables
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
       PARAMETER (MINREAL =-1.0E30)
       PARAMETER (MAXREAL = 1.0E30)
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
       SAVE /gzvchr/
*
*-----------------------------------------------------------------------
*  declare dummy arguments
*-----------------------------------------------------------------------
       INTEGER    view
*-----------------------------------------------------------------------
       LOGICAL    ret
*-----------------------------------------------------------------------
       ret=.FALSE.
*
       IF(view.LT.1 .OR. view.GT.MAXVIEW) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'01    GZDCTV  INVALID VIEW INDEX'
       ENDIF
*
*
*
       IF(ret) THEN
         WRITE(DEVIND,*)' GZDCTV  NO ACTION IS PERFORMED'
         RETURN
       ENDIF
*
*
*
       vactiv(view)=.FALSE.
       RETURN
*---------------------------------------------------------   end gzdctv
       END
*
*
************************************************************************
C@PROCESS OPT(3) NOSDUMP NOGOSTMT
       SUBROUTINE gzsfgl(objind, fcetype, gloss)
*
************************************************************************
*
*  Set Face Gloss
*
*  dummy arguments: objind - objec index ;
*                   fcetype - face type (FRONT or BACK);
*                   gloss - surface gloss;
*
*  function:  this subroutine sets surface gloss for specifed object;
*
*  called by: main program
*
************************************************************************
       IMPLICIT NONE
*-----------------------------------------------------------------------
*  iclude global constants and variables
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
       PARAMETER (MINREAL =-1.0E30)
       PARAMETER (MAXREAL = 1.0E30)
*
       PARAMETER (ONLINE = 1)
       PARAMETER (OFFLINE= 2)
c      include(gzclrs)
***********************************************************************
*
*   autor:  Zenon Zowierucha, zdv069
*
*   file name, type :  GZCLRS COPY
*
*   function :  this file declares global variables which describe
*               patch attributes;
*               this file is destinated for use with INCLUDE
*               directive;
*
*
***********************************************************************
*
*                          front face colors; fouclr(1,i),
*                          fouclr(2,i),fouclr(3,i) are recpectivily
*                          RGB components of front face color of
*                          object number i;
       REAL       fouclr(1:3, 1:MAXOBJ)
*
*                          back face colors;
       REAL       finclr(1:3, 1:MAXOBJ)
*
*                          front edge colors;
       REAL       eouclr(1:3, 1:MAXOBJ)
*
*                          back edge color ;
       REAL       einclr(1:3, 1:MAXOBJ)
*
*
*                          front face colors low indices;
*                          in device look up table, high index is equal
*                          fouind(i)+foucno(i)-1);
       INTEGER    fouind(1:MAXOBJ)
*
*                          number of front face colors (hues);
       INTEGER    foucno(1:MAXOBJ)
*
*
*                          back face colors low indices;
*                          (high index is equal finind(i)+fincno(i)-1);
       INTEGER    finind(1:MAXOBJ)
*
*                          number of back face colors (hues);
       INTEGER    fincno(1:MAXOBJ)
*
*
*                          front edge color indices;
       INTEGER    eouind(1:MAXOBJ)
*
*                          back edge color indices;
       INTEGER    einind(1:MAXOBJ)
*
*                          front face glosses; (in range <1:MAXGLS>,
*                          1 - lowest, MAXGLS - highest gloss);
       INTEGER    fougls(1:MAXOBJ)
*
*                          back face glosses;
       INTEGER    fingls(1:MAXOBJ)
*
       INTEGER    dclro(1:MAXOBJ)
       INTEGER    dclri(1:MAXOBJ)
*
*                          look up table;
       REAL       rgb(1:3, 0:MAXCLR-1)
*
*
       COMMON/gzclrs/fouind, finind,
     +               eouind, einind, foucno, fincno,
     +               fouclr, finclr,
     +               eouclr, einclr,
     +               fougls, fingls, dclro, dclri, rgb
      SAVE /gzclrs/
*
*-----------------------------------------------------------------------
*  declare dummy arguments
*-----------------------------------------------------------------------
       INTEGER    objind, fcetype, gloss
*-----------------------------------------------------------------------
       LOGICAL    ret
*-----------------------------------------------------------------------
       ret=.FALSE.
*
       IF(objind.LT.1   .OR. objind.GT.MAXOBJ) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'02    GZSFGL  INVALID OBJECT INDEX'
       ENDIF
*
       IF(fcetype.NE.FRONT .AND. fcetype.NE.BACK) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'03    GZSFGL  INVALID FACE TYPE'
       ENDIF
*
       IF(gloss.LT.1 .OR. gloss.GT.MAXGLS) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'33    GZSFGL  GLOSS VALUE OUT OF RANGE'
       ENDIF
*
*
       IF(ret) THEN
         WRITE(DEVIND,*)' GZSFGL  NO ACTION IS PERFORMED'
         RETURN
       ENDIF
*
       IF(fcetype.EQ.BACK) THEN
          fougls(objind)=gloss
       ELSE
          fingls(objind)=gloss
       ENDIF
*
*

       RETURN
*---------------------------------------------------------   end gzsfgl
       END
C@PROCESS OPT(3) NOSDUMP NOGOSTMT
      SUBROUTINE GZMIRR(MIRROR,ICOL)
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
       PARAMETER (MINREAL =-1.0E30)
       PARAMETER (MAXREAL = 1.0E30)
*
       PARAMETER (ONLINE = 1)
       PARAMETER (OFFLINE= 2)
c      include(gzclrs)
***********************************************************************
*
*   autor:  Zenon Zowierucha, zdv069
*
*   file name, type :  GZCLRS COPY
*
*   function :  this file declares global variables which describe
*               patch attributes;
*               this file is destinated for use with INCLUDE
*               directive;
*
*
***********************************************************************
*
*                          front face colors; fouclr(1,i),
*                          fouclr(2,i),fouclr(3,i) are recpectivily
*                          RGB components of front face color of
*                          object number i;
       REAL       fouclr(1:3, 1:MAXOBJ)
*
*                          back face colors;
       REAL       finclr(1:3, 1:MAXOBJ)
*
*                          front edge colors;
       REAL       eouclr(1:3, 1:MAXOBJ)
*
*                          back edge color ;
       REAL       einclr(1:3, 1:MAXOBJ)
*
*
*                          front face colors low indices;
*                          in device look up table, high index is equal
*                          fouind(i)+foucno(i)-1);
       INTEGER    fouind(1:MAXOBJ)
*
*                          number of front face colors (hues);
       INTEGER    foucno(1:MAXOBJ)
*
*
*                          back face colors low indices;
*                          (high index is equal finind(i)+fincno(i)-1);
       INTEGER    finind(1:MAXOBJ)
*
*                          number of back face colors (hues);
       INTEGER    fincno(1:MAXOBJ)
*
*
*                          front edge color indices;
       INTEGER    eouind(1:MAXOBJ)
*
*                          back edge color indices;
       INTEGER    einind(1:MAXOBJ)
*
*                          front face glosses; (in range <1:MAXGLS>,
*                          1 - lowest, MAXGLS - highest gloss);
       INTEGER    fougls(1:MAXOBJ)
*
*                          back face glosses;
       INTEGER    fingls(1:MAXOBJ)
*
       INTEGER    dclro(1:MAXOBJ)
       INTEGER    dclri(1:MAXOBJ)
*
*                          look up table;
       REAL       rgb(1:3, 0:MAXCLR-1)
*
*
       COMMON/gzclrs/fouind, finind,
     +               eouind, einind, foucno, fincno,
     +               fouclr, finclr,
     +               eouclr, einclr,
     +               fougls, fingls, dclro, dclri, rgb
      SAVE /gzclrs/
*
C---- MIRROR zwischen 1 und 5: 1 Starke Spiegelreflexion, 5 geringe
C---- COLOR TABLE FOR SURFACES
C---- WENN HIER AENDERUNG NHELL, DANN AUCH IN DRAWSURFACE.
C---- Q : BETWEEN 0 AND 1
      PARAMETER(NHELL=60,NCOL=128,Q=.25,MIRMAX=5,
     >          MIFA=NHELL/(5*(MIRMAX-1)) )
      REAL    RGB1(3,NCOL)
      EQUIVALENCE (rgb, rgb1)
COTA: RED-GREEN, GREEN-RED, BLUE-RED, RED-BLUE, GREEN-BLUE, BLUE-GREEN
      REAL    COTA(6,6)
      DATA    COTA      /.3,.04,.03,  .04,.3,.03,
     >                   .03,.3,.03,  .3,.04,.03,
     >                   .03,.04,.3,  .3,.03,.04,
     >                   .3,.03,.04,  .03,.04,.3,
     >                   .03,.3,.04,  .03,.04,.3,
     >                   .03,.04,.3,  .03,.3,.04/
      K=NHELL-(MIRROR-1)*MIFA
      H1=(NHELL-K)/FLOAT(NHELL-1)
C---- Alle 3 Farbanteile sind beteiligt
      DO 4 J=1,3
         IF (COTA(J,ICOL) .GT. 0.1 ) THEN
            A = COTA(J,ICOL)
         ELSE
            A = COTA(J,ICOL)*(MIRMAX+1-MIRROR)/MIRMAX
         ENDIF
         IF (COTA(J+3,ICOL) .GT. 0.1 ) THEN
            B = COTA(J+3,ICOL)
         ELSE
            B = COTA(J+3,ICOL)*(MIRMAX+1-MIRROR)/MIRMAX
         ENDIF
         IF (K.NE.NHELL) THEN
            ALN=LOG(A)
            AH1=EXP(ALN*H1)
            A0=AH1*(1-Q)+Q
            A1=2*(AH1*Q-Q)/H1-AH1*ALN
            A2=(AH1*ALN-(AH1*Q-Q)/H1)/H1
            BLN=LOG(B)
            BH1=EXP(BLN*H1)
            B0=BH1*(1-Q)+Q
            B1=2*(BH1*Q-Q)/H1-BH1*BLN
            B2=(BH1*BLN-(BH1*Q-Q)/H1)/H1
         ENDIF
         DO  5 I = 1,NHELL
C----------- HELL von 1 nach 0
             HELL=(NHELL-I)/FLOAT(NHELL-1)
             IF (I .LE. K) THEN
C----------     A**HELL zwischen A und 1
                RGB1(J,8+ I)          = A**HELL
                RGB1(J,NHELL + 8 + I) = B**HELL
             ELSE
                RGB1(J,8+ I) =          A0+HELL*(A1+HELL*A2)
                RGB1(J,NHELL + 8 + I) = B0+HELL*(B1+HELL*B2)
             ENDIF
 5       CONTINUE
 4    CONTINUE
      END
*
*
*
*
       SUBROUTINE gzdefclr(clr, index)
       REAL    clr(1:3)
       INTEGER   index
*
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
       PARAMETER (MINREAL =-1.0E30)
       PARAMETER (MAXREAL = 1.0E30)
*
       PARAMETER (ONLINE = 1)
       PARAMETER (OFFLINE= 2)
c      include(gzclrs)
***********************************************************************
*
*   autor:  Zenon Zowierucha, zdv069
*
*   file name, type :  GZCLRS COPY
*
*   function :  this file declares global variables which describe
*               patch attributes;
*               this file is destinated for use with INCLUDE
*               directive;
*
*
***********************************************************************
*
*                          front face colors; fouclr(1,i),
*                          fouclr(2,i),fouclr(3,i) are recpectivily
*                          RGB components of front face color of
*                          object number i;
       REAL       fouclr(1:3, 1:MAXOBJ)
*
*                          back face colors;
       REAL       finclr(1:3, 1:MAXOBJ)
*
*                          front edge colors;
       REAL       eouclr(1:3, 1:MAXOBJ)
*
*                          back edge color ;
       REAL       einclr(1:3, 1:MAXOBJ)
*
*
*                          front face colors low indices;
*                          in device look up table, high index is equal
*                          fouind(i)+foucno(i)-1);
       INTEGER    fouind(1:MAXOBJ)
*
*                          number of front face colors (hues);
       INTEGER    foucno(1:MAXOBJ)
*
*
*                          back face colors low indices;
*                          (high index is equal finind(i)+fincno(i)-1);
       INTEGER    finind(1:MAXOBJ)
*
*                          number of back face colors (hues);
       INTEGER    fincno(1:MAXOBJ)
*
*
*                          front edge color indices;
       INTEGER    eouind(1:MAXOBJ)
*
*                          back edge color indices;
       INTEGER    einind(1:MAXOBJ)
*
*                          front face glosses; (in range <1:MAXGLS>,
*                          1 - lowest, MAXGLS - highest gloss);
       INTEGER    fougls(1:MAXOBJ)
*
*                          back face glosses;
       INTEGER    fingls(1:MAXOBJ)
*
       INTEGER    dclro(1:MAXOBJ)
       INTEGER    dclri(1:MAXOBJ)
*
*                          look up table;
       REAL       rgb(1:3, 0:MAXCLR-1)
*
*
       COMMON/gzclrs/fouind, finind,
     +               eouind, einind, foucno, fincno,
     +               fouclr, finclr,
     +               eouclr, einclr,
     +               fougls, fingls, dclro, dclri, rgb
      SAVE /gzclrs/
*
*
*
       rgb(1,index)=clr(1)
       rgb(2,index)=clr(2)
       rgb(3,index)=clr(3)
*
       RETURN
       END
