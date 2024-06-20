************************************************************************
C@PROCESS OPT(3) IL(DIM) NOSDUMP NOGOSTMT
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
       SUBROUTINE gzdrf1(xyz1,par,shd1)
************************************************************************
*
*  file name, type:  GZPLOT FORTRAN
*
*  dummy arguments: xyz - a triangle vertices coordinates;
*                   (DC after perspective projection);
*                   xyz(1,j),xyz(2,j),xyz(3,j) are x,y,z coordinates of
*                   vertex nr j of the triangle;
*
*                   par - triangle plane equation parameters
*                         (par(1)*x + par(2)*y + par(3)*z + par(4)=0);
*
*                   shd - illumination intensity of triangle vertices;
*                         shd(1,j),...,shd(3,j) are RGB components of
*                         illumination intensity of a vertex nr j;
*
*  function:  this subroutine fills a shaded triangle with Gouraud's
*             shading algorithm (the triangle was perspectivily
*             projected);
*
*  called by:
*
*  calls: gzlin1 (file GZPLOT FORTRAN)
*
*  CHANGED. M. BUSCH    12.4.91
*                       INTEGER*4  --> INTEGER
*                       REAL   *4  --> REAL
*  CHANGED. M. BUSCH    6.1.93 Common Bloecke nicht mehr ueber include
*                       sondern als Source eingefuegt
*  CHANGED. G. Groten   5.1.93 Zufallszahlengenerator GZNGEN/GZSEED
*  CHANGED. M. BUSCH    10.3.93 MINREAL, MAXREAL +- 1.e.38 statt
*                               +-1.e50 (Anpassung an IEEE Format)
************************************************************************
       IMPLICIT NONE
*-----------------------------------------------------------------------
*  declare dummy arguments
*-----------------------------------------------------------------------
       REAL       xyz1(1:3, 1:3), par(1:4), shd1(1:3, 1:3)
*
*-----------------------------------------------------------------------
*  define common blocs
*-----------------------------------------------------------------------
*                          common with gzlin1 subroutine
       REAL    a, b, czp, mia, c,d
       COMMON/gzl/a,b,czp,mia,c,d
       SAVE /gzl/
*
*                          projection reference point (DC)
       REAL       dcprp(1:3)
*
*                          view plane distance (DC)
       REAL       dcdist
*
       COMMON/gzdc/dcprp,dcdist
       SAVE /gzdc/
*
       REAL    wdw(1:4),vpt(1:4),prp1(1:3),sc, dst
       COMMON/gzw/wdw,vpt,prp1,sc,dst
        SAVE /gzw/
c      include(rdoff)
       REAL     ROUNDOFF
       PARAMETER (ROUNDOFF=1.0E-10)
*
*-----------------------------------------------------------------------
*  declare local variables
*-----------------------------------------------------------------------
       REAL    xyz(1:3,1:3),shd(1:3,1:3)
       REAL       temp(1:3)
*
       REAL        dy(1:3)
       REAL        x(1:3)
       REAL        y(1:3)
       REAL        s(1:3,1:3)
       REAL        sx(1:3)
       REAL        ss(1:3,1:3)
       REAL        yi
       INTEGER     edg, no,i,j,k
*-----------------------------------------------------------------------
*
       a=par(1)
       b=par(2)
       c=par(3)
       d=par(4)
       czp=par(3)*dcdist+ par(4)
       mia=a*dcprp(1) + b*dcprp(2) + c*dcprp(3) + par(4)
       IF(ABS(mia).LE.ROUNDOFF) RETURN
*
*---------------------------------------------   sort triangle vertices
       DO 10 i=1,3
        DO 11 j=1,3
         xyz(j,i)=xyz1(j,i)
         shd(j,i)=shd1(j,i)
 11     CONTINUE
 10    CONTINUE
       DO 100 k=2, 3
          DO 110 i=3, k, -1
            IF(xyz(2, i).GT.xyz(2, i-1)) THEN
               DO 120 j=1, 3
                   temp(j)=xyz(j, i-1)
                   xyz(j, i-1)=xyz(j, i)
                   xyz(j, i)=temp(j)
                   temp(j)=shd(j, i-1)
                   shd(j, i-1)=shd(j, i)
                   shd(j, i)=temp(j)
 120           CONTINUE
            ENDIF
 110      CONTINUE
 100   CONTINUE
*
*-----------------------------------------------------------------------
       y(1) =AINT(xyz(2, 1)) + 0.5
       y(2) =AINT(xyz(2, 2)) + 0.5
       y(3) =AINT(xyz(2, 3)) + 0.5
*
       dy(1) = y(1)-y(3)
       dy(2) = y(1)-y(2)
       dy(3) = y(2)-y(3)
*
       edg=0
       no=0
*
       IF(dy(1).GT.0.0) THEN
          sx(1)= (xyz(1,3)-xyz(1,1))/dy(1)
          ss(1,1)=(shd(1,3)-shd(1,1))/dy(1)
          ss(2,1)=(shd(2,3)-shd(2,1))/dy(1)
          ss(3,1)=(shd(3,3)-shd(3,1))/dy(1)
*
          x(1)=xyz(1,1) + 0.5*sx(1)
*
          s(1,1)=shd(1,1) + 0.5*ss(1,1)
          s(2,1)=shd(2,1) + 0.5*ss(2,1)
          s(3,1)=shd(3,1) + 0.5*ss(3,1)
*
          no=no+1
          edg=edg+2
       ENDIF
*
       IF(dy(2).GT.0.0) THEN
          sx(2)= (xyz(1,2)-xyz(1,1))/dy(2)
          ss(1,2)=(shd(1,2)-shd(1,1))/dy(2)
          ss(2,2)=(shd(2,2)-shd(2,1))/dy(2)
          ss(3,2)=(shd(3,2)-shd(3,1))/dy(2)
*
          x(2)=xyz(1,1) + 0.5*sx(2)
*
          s(1,2)=shd(1,1) + 0.5*ss(1,2)
          s(2,2)=shd(2,1) + 0.5*ss(2,2)
          s(3,2)=shd(3,1) + 0.5*ss(3,2)
*
          no=no+1
          edg=edg+4
       ENDIF
*
       IF(dy(3).GT.0.0) THEN
          sx(3)= (xyz(1,3)-xyz(1,2))/dy(3)
          ss(1,3)=(shd(1,3)-shd(1,2))/dy(3)
          ss(2,3)=(shd(2,3)-shd(2,2))/dy(3)
          ss(3,3)=(shd(3,3)-shd(3,2))/dy(3)
*
          x(3)=xyz(1,2) + 0.5*sx(3)
*
          s(1,3)=shd(1,2) + 0.5*ss(1,3)
          s(2,3)=shd(2,2) + 0.5*ss(2,3)
          s(3,3)=shd(3,2) + 0.5*ss(3,3)
*
          no=no+1
          edg=edg+8
       ENDIF
       IF(no.EQ.0) THEN
                   CALL gzplt1(NINT(xyz(1,1)), NINT(xyz(2,1)),
     +             (xyz(3,1)+xyz(3,2)+xyz(3,3))/3.0, shd(1,1))
       ELSE
           IF(no.EQ.2) THEN
              i=(edg-2)/4 + 1
           ELSE
              i=2
           ENDIF
           yi=AINT(y(1))
 50        IF(yi.GE.y(i)) THEN
                  IF(x(1).LT.x(i)) THEN
                  CALL gzlin1(yi,ANINT(x(1)),s(1,1),
     +                           ANINT(x(i)),s(1,i))
                  ELSE
                  CALL gzlin1(yi,ANINT(x(i)),s(1,i),
     +                         ANINT(x(1)),s(1,1))
                  ENDIF
             x(1)=x(1)+sx(1)
             s(1,1)=s(1,1)+ss(1,1)
             s(2,1)=s(2,1)+ss(2,1)
             s(3,1)=s(3,1)+ss(3,1)
*
             x(i)=x(i)+sx(i)
             s(1,i)=s(1,i)+ss(1,i)
             s(2,i)=s(2,i)+ss(2,i)
             s(3,i)=s(3,i)+ss(3,i)
*
             yi=yi-1.0
             GOTO 50
         ENDIF
*
*
             IF(no.EQ.3) THEN
*
 60            IF(yi.GE.y(3)) THEN
                  IF(x(1).LT.x(3)) THEN
                  CALL gzlin1(yi,ANINT(x(1)),s(1,1),
     +                           ANINT(x(3)),s(1,3))
                  ELSE
                  CALL gzlin1(yi,ANINT(x(3)),s(1,3),
     +                           ANINT(x(1)),s(1,1))
                  ENDIF
                 x(1)=x(1)+sx(1)
                 s(1,1)=s(1,1)+ss(1,1)
                 s(2,1)=s(2,1)+ss(2,1)
                 s(3,1)=s(3,1)+ss(3,1)
*
                 x(3)=x(3)+sx(3)
                 s(1,3)=s(1,3)+ss(1,3)
                 s(2,3)=s(2,3)+ss(2,3)
                 s(3,3)=s(3,3)+ss(3,3)
*
                 yi=yi-1.0
                 GOTO 60
              ENDIF
           ENDIF
*         ENDIF
         ENDIF
 643   RETURN
       END
*----------------------------------------------------------   end gzdrf1
************************************************************************
C@PROCESS OPT(3) IL(DIM) NOSDUMP NOGOSTMT DC(FELD)
       SUBROUTINE gzplt1(x,y,z,si)
************************************************************************
*
*  dummy arguments:
*       x  -  x-coordinate of the pixel (DCS);
*       y  -  y-coordinate of the pixel (DCS);
*       z  -  z-coordinate of a point (DCS);
*       si -  color of the pixel; si(1),si(2),si(3) are RGB components;
*
*  function:  this subroutine set one pixel of a face in frame buffer;
*
*  called by: gzlin1 (file GZDRF1 FORTRAN)
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
        SAVE /gzzbuff/, /feld/
*
*
*-----------------------------------------------------------------------
*  declare dummy arguments
       INTEGER    x, y
       REAL       z
       REAL       si(1:3)
*-----------------------------------------------------------------------
*                          current face color low index
       INTEGER    cfcli
*
*                          current face color number
       INTEGER    cfcno
*
*                          current activ color component (1-red,
*                          2-green, 3-blue);
       INTEGER    cacc
*
       COMMON/gzcurr/cfcli, cfcno, cacc
       REAL      dzz
       COMMON/gzpt/dzz
       SAVE /gzcurr/,/gzpt/
*-----------------------------------------------------------------------
*  declare local constants
       REAL    MINCLS
       PARAMETER(MINCLS=0.90)
*-----------------------------------------------------------------------
       INTEGER   colind
       REAL      rand, delta, col
*      REAL      sii
       DOUBLE PRECISION r
       IF(z+dzz.GE.zbuff(x, y)) THEN
          zbuff(x, y)=z
*
*            IF(si(1).GE.MINCLS.AND.si(2).GE.MINCLS.AND.si(3).GE.MINCLS)
*    +       THEN
*               sii=(si(1)+si(2)+si(3))/3.0
*               sii=0.79
*               CALL GZNGEN(R,1)
*               rand=r
*               IF(sii.GT.rand) THEN
*                     colind= cfcno
*               ELSE
*                     colind= cfcno - 1
*               ENDIF
*               picture(x, y)= CHAR(colind + cfcli)
*            ELSE
                col= (cfcno-1) * si(cacc) - 0.45
                delta= col - AINT(col)
                CALL GZNGEN(r,1)
                rand=r
                IF(delta.GT.rand) THEN
                      colind= INT(col +1.0)
                ELSE
                      colind= INT(col)
                ENDIF
                picture(x, y)= CHAR(colind + cfcli)
             ENDIF
*      ENDIF
       RETURN
       END
*---------------------------------------------------------   end gzplt1
************************************************************************
C@PROCESS OPT(3) IL(DIM) NOSDUMP NOGOSTMT DC(FELD)
       SUBROUTINE gzplt2(x,y,z,si)
************************************************************************
*
*  dummy arguments:
*       x  -  x-coordinate of the pixel (DCS);
*       y  -  y-coordinate of the pixel (DCS);
*       z  -  z-coordinate of a point (DCS);
*       si -  color of the pixel; si(1),si(2),si(3) are RGB components;
*
*  function:  this subroutine set one pixel of a face in frame buffer;
*
*  called by: gzlin2 (file GZDRF1 FORTRAN)
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
*-----------------------------------------------------------------------
*  declare dummy arguments
       INTEGER    x, y
       REAL       z
       REAL       si(1:3)
*-----------------------------------------------------------------------
*                          current face color low index
       INTEGER    cfcli
*
*                          current face color number
       INTEGER    cfcno
*
*                          current activ color component (1-red,
*                          2-green, 3-blue);
       INTEGER    cacc
*
       COMMON/gzcurr/cfcli, cfcno, cacc
       REAL      dzz
       COMMON/gzpt/dzz
        SAVE /gzcurr/,/gzpt/
*-----------------------------------------------------------------------
*  declare local constants
       REAL    MINCLS
       PARAMETER(MINCLS=0.99)
*-----------------------------------------------------------------------
       INTEGER   colind
       REAL      rand, delta, col
*      REAL      sii, dz
       DOUBLE PRECISION r
*     if(z+dzz.lt.zbuff(x,y).and.picture(x,y).gt.
*    +char(64).and.dzz.gt.0.0)
*    + write(*,*)' front z', z+dzz, 'back z',zbuff(x,y)
*
       IF(z+dzz.GE.zbuff(x, y)) THEN
          zbuff(x, y)=z
*
*            IF(si(1).GE.MINCLS.AND.si(2).GE.MINCLS.AND.si(3).GE.MINCLS)
*    +       THEN
*               sii=(si(1)+si(2)+si(3))/16.0
*               sii=0.79
*               CALL GZNGEN(r,1)
*               rand=r
*               IF(sii.GT.rand) THEN
                      colind= cfcno-1
*               ELSE
*                     colind= cfcno-1
*               ENDIF
*               picture(x, y)= CHAR(colind + cfcli)
*            ELSE
                col= (cfcno-1) * si(cacc) - 0.45
                delta= col - AINT(col)
                CALL gzngen(r,1)
                rand=r
                IF(delta.GT.rand) THEN
                      colind= INT(col +1.0)
                ELSE
                      colind= INT(col)
                ENDIF
                picture(x, y)= CHAR(colind + cfcli)
*            ENDIF
       ENDIF
       RETURN
       END
*---------------------------------------------------------   end gzplt1
************************************************************************
C@PROCESS OPT(3) IL(DIM) NOSDUMP NOGOSTMT DC(FELD)
       SUBROUTINE gzlin1(y, x1,s1, x2,s2)
************************************************************************
*
*  dummy arguments:
*      y - y coordinate of line segment;
*      x1,x2 - x coordinates of ends of line segment (x2>x1)
*      s1,s2 - illumination intensites if edpoints of line segment;
*
*  function: this subroutine draws a shaded horizontal line segment
*            from punkt (x1,y) to (x2,y); this subroutine is suitable
*            only for perspective projected lines;
*
************************************************************************
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
c      include(rdoff)
       REAL     ROUNDOFF
       PARAMETER (ROUNDOFF=1.0E-10)
*
       REAL       y
       REAL       x1
       REAL       s1(1:3)
       REAL       x2
       REAL       s2(1:3)
*-----------------------------------------------   define common blocs
       REAL       a, b, c, d
       REAL       czp
       REAL       mia
       COMMON/gzl/a, b, czp, mia, c,d
*
*                          projection reference point (DCS)
       REAL       dcprp(1:3)
*
*                          view plane distance (DCS)
       REAL       dcdist
       COMMON/gzdc/dcprp, dcdist
       REAL    wdw(1:4),vpt(1:4),prp1(1:3),sc, dst
       COMMON/gzw/wdw,vpt,prp1,sc,dst
       SAVE /gzl/,/gzdc/,/gzw/
*
*-------------------------------------------   declare local variables
       REAL       z
       REAL       si(1:3)
       REAL       ds(1:3)
       REAL       mi
       REAL       dmi
       REAL       x
       REAL       dx
       INTEGER    yi
*
*-----------------------------------------------------------------------
       mi=(a*x1 + b*y + czp)/mia
       IF(ABS(1.0-mi).LE.ROUNDOFF) RETURN
       dx=x2-x1
       IF(dx.EQ.0.0) THEN
         CALL gzplt1(NINT(x1),MAXPKY - NINT(y) - 1,
     +   (dcdist-mi*dcprp(3))/(1.0-mi),s1)
         RETURN
       ENDIF
       dmi=a/mia
       ds(1)=(s2(1)-s1(1))/dx
       ds(2)=(s2(2)-s1(2))/dx
       ds(3)=(s2(3)-s1(3))/dx
       si(1)=s1(1)
       si(2)=s1(2)
       si(3)=s1(3)
       yi= MAXPKY - NINT(y) - 1
       x=x1

 310   IF(x.LE.x2) THEN
          z = (dcdist-mi*dcprp(3))/(1.0-mi)
*         z = (-a*x - b*y - d)/c
          CALL gzplt1(NINT(x), yi, z, si)
*

          si(1)=si(1)+ds(1)
          si(2)=si(2)+ds(2)
          si(3)=si(3)+ds(3)
*
*      mi=(a*x  + b*y + czp)/mia
          mi=mi+dmi
          x=x+1.0
          GOTO 310
       ENDIF
*
       RETURN
*-------------------------------------------------------   end gzlin1
       END
*
************************************************************************
C@PROCESS OPT(3) IL(DIM) NOSDUMP NOGOSTMT DC(FELD)
       SUBROUTINE gzdda1(x1,y1, x2,y2, par, color)
*
************************************************************************
*
*  dummy arguments:
*      x1,y1  - coordinates of begin point of line segment (DC after
*               perspective projection);
*      x2,y2  - coordinates of end point of line segment;
*      par    - line segment plane parameters
*               (par(1)*xi + par(2)*yi + par(3)*zi + par(4)=0, where
*               xi,yi,zi belong to the line segment;
*
*  function:  this subroutine draws a perspective projected line segment
*             in z-buffer with color 'color';
*
************************************************************************
       IMPLICIT NONE
*
c      include (gzcst11)
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
c      include (gzcst21)
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
c      include(rdoff)
       REAL     ROUNDOFF
       PARAMETER (ROUNDOFF=1.0E-10)
*-------------------------------------------   declare dummy arguments
       REAL       x1,y1, x2,y2, par(1:4)
       INTEGER    color
*
*-------------------------------------------   define common blocs
*
       REAL       dcprp(1:3), dcdist
       COMMON/gzdc/ dcprp, dcdist
       SAVE /gzdc/
*
*-------------------------------------------   declare local variables
*
       REAL       a,b,czp,mia, mi
       REAL       yy1,yy2, xx1,xx2, dx,dy, absdx,absdy, sx,sy, x,y
       REAL       rob, temp
       INTEGER    i, steps
*
*-----------------------------------------------------------------------
*
       a=par(1)
       b=par(2)
       czp=par(3)*dcdist + par(4)
       mia=a*dcprp(1) + b*dcprp(2) + par(3)*dcprp(3) + par(4)
       IF(ABS(mia).LT.ROUNDOFF) RETURN
*
       yy1=AINT(y1) + 0.5
       yy2=AINT(y2) + 0.5
*
       xx1=x1
       xx2=x2
*
       IF(yy1.LT.yy2) THEN
         rob=yy1
         yy1=yy2
         yy2=rob
         rob=xx1
         xx1=xx2
         xx2=rob
       ENDIF
*
       dx=xx2-xx1
       dy=yy2-yy1
*
       absdx=ABS(dx)
       absdy=ABS(dy)
*
       IF(absdx.LT.1.0 .AND. absdy.LT.1.0) THEN
         mi=(a*ANINT(x1) + b*ANINT(y1) + czp)/mia
        CALL gzplt5(NINT(x1),-NINT(y1)+MAXPKY-1,
     +  (dcdist-mi*dcprp(3))/(1.0-mi), color)
         RETURN
       ENDIF
*
       temp=ANINT(AMAX1(absdx,absdy))
       sx=dx/temp
       sy=dy/temp
       x=xx1 + 0.5*sx
       y=AINT(yy1)
*
       steps=NINT(temp)
*
       DO 100 i=1,steps
         mi=(a*ANINT(x) + b*ANINT(y) + czp)/mia
*
         CALL gzplt5(NINT(x),-NINT(y)+MAXPKY-1,
     +    (dcdist-mi*dcprp(3))/(1.0-mi),color)
*
         x=x+sx
         y=y+sy
 100   CONTINUE
*
*-----------------------------------------------------   end gzdda1
       END
*
*
************************************************************************
C@PROCESS OPT(3) IL(DIM) NOSDUMP NOGOSTMT DC(FELD)
       SUBROUTINE gzplt5(x,y,z,color)
*
************************************************************************
*
*  dummy arguments:
*      x,y - coordinates of a pixel on the screen (DC);
*      z   - z coordinate of the pixel;
*      color - color index;
*
*  function:
*      this subroutine draws a pixel in z-buffer;
*
************************************************************************
       IMPLICIT NONE
       INTEGER    x,y
       REAL       z
       INTEGER    color
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
*
*
*
       REAL      dzz
       COMMON/gzpt/dzz
        SAVE /gzzbuff/,/feld/,/gzpt/
*-----------------------------------------------------------------------
*
       IF(z+dzz+0.9.GE.zbuff(x,y)) THEN
              zbuff(x,y)=z
              picture(x,y)= CHAR(color)
       ENDIF
*------------------------------------------------------   end gzplt2
       END
*
*
*
*
*
*
*
************************************************************************
C@PROCESS OPT(3) IL(DIM) NOSDUMP NOGOSTMT DC(FELD)
       SUBROUTINE gzshade(xyz, vnrm, rgb, rgb1, gloss, prp, prtype,
     +            srcpos, srctyp, srcpwr, ambnt)
************************************************************************
*
*  dummy parameters:
*      xyz - point coordinates (WCS);
*      vnrm - normal vektor of the plane of the point;
*      rgb - color components of the point(input);
*      rgb1- color components of the point(output);
*      gloss - gloss of the point (in range <1:MAXGLS>);
*      prp - projection reference point (WCS);
*      prptyp - type of projection (PARALEL or PERSPECTIVE);
*      srcpos - position of light source (or gzunity vektor if srctyp=SU
*      srctyp - type of light source (SUN or POINT);
*      srcpwr - power of light source (in range <0:1>);
*      ambnt - ambient illumination intensity(in range <0:1>);
*
*  function:  this subroutine computes a illumination intensity
*             (RGB values) in dependence on light source and observer
*             (prp) position;
*
************************************************************************
       IMPLICIT NONE
*  include global constants
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
*
*  declare dummy arguments
       REAL       xyz(1:3)
       REAL       vnrm(1:3)
       REAL       rgb(1:3)
       REAL       rgb1(1:3)
       INTEGER    gloss
       REAL       prp(1:3)
       INTEGER    prtype
       REAL       srcpos(1:3)
       INTEGER    srctyp
       REAL       srcpwr
       REAL       ambnt
*-----------------------------------------------------------------------
*  declare local variables
*
*                          light source vektor
       REAL       li(1:3)
*
*                          observer vektor
       REAL       e(1:3)
*
       REAL       h(1:30)
       REAL       cosi, cosn, dummy,gzunity
*
*-----------------------------------------------------------------------
       IF(prtype.EQ.PARALLEL) THEN
          e(1)=0.0
          e(2)=0.0
          e(3)=1.0
       ELSE
          e(1) = prp(1) - xyz(1)
          e(2) = prp(2) - xyz(2)
          e(3) = prp(3) - xyz(3)
          dummy = gzunity(e)
       ENDIF
       IF(srctyp.EQ.POINT) THEN
*        --------------------------------------    a point light source
          li(1) = srcpos(1) - xyz(1)
          li(2) = srcpos(2) - xyz(2)
          li(3) = srcpos(3) - xyz(3)
          dummy = gzunity(li)
       ELSE
*        -------------------------------------   a paralel light source
          li(1)=srcpos(1)
          li(2)=srcpos(2)
          li(3)=srcpos(3)
       ENDIF
*
       h(1)=e(1) + li(1)
       h(2)=e(2) + li(2)
       h(3)=e(3) + li(3)
       dummy = gzunity(h)
       dummy=gzunity(vnrm)
       cosi=vnrm(1)*li(1) + vnrm(2)*li(2) + vnrm(3)*li(3)
       IF(cosi.GT.0.0) THEN
            cosn = vnrm(1)*h(1) + vnrm(2)*h(2) + vnrm(3)*h(3)
       ELSE
            cosi = 0.0
            cosn = 0.0
       ENDIF
       IF(cosn .GT. 0.111) THEN
        dummy = srcpwr * cosn**gloss
*       dummy=0.0
       ELSE
        dummy=0.0
       ENDIF
       rgb1(1) = ambnt*rgb(1) + srcpwr*rgb(1)*cosi + dummy
       IF(rgb1(1).GT.1.0) rgb1(1)=1.0
       rgb1(2) = ambnt*rgb(2) + srcpwr*rgb(2)*cosi + dummy
       IF(rgb1(2).GT.1.0) rgb1(2)=1.0
       rgb1(3) = ambnt*rgb(3) + srcpwr*rgb(3)*cosi + dummy
       IF(rgb1(3).GT.1.0) rgb1(3)=1.0
*
       RETURN
       END
*-------------------------------------------------------    end gzshade
*
*
*
*
*
*
*
*
*
*
*
*
*
*
*
*
*
*
*
*
************************************************************************
C@PROCESS OPT(3) IL(DIM) NOSDUMP NOGOSTMT DC(FELD)
       SUBROUTINE gzdrf2(xyz1,par,nrm1)
************************************************************************
*
*  file name, type:  GZPLOT FORTRAN
*
*  dummy arguments: xyz - a triangle vertices coordinates;
*                   (DC after perspective projection);
*                   xyz(1,j),xyz(2,j),xyz(3,j) are x,y,z coordinates of
*                   vertex nr j of the triangle;
*
*                   par - triangle plane equation parameters
*                         (par(1)*x + par(2)*y + par(3)*z + par(4)=0);
*
*                   nrm - vetex normals;
*
*  function:  this subroutine fills a shaded triangle with Phong's
*             shading algorithm (the triangle was perspectivily
*             projected);
*
*  called by:
*
*  calls: gzlin2 (file GZPLOT FORTRAN)
*
************************************************************************
       IMPLICIT NONE
c      include(rdoff)
       REAL     ROUNDOFF
       PARAMETER (ROUNDOFF=1.0E-10)
*-----------------------------------------------------------------------
*  declare dummy arguments
*-----------------------------------------------------------------------
       REAL       xyz1(1:3, 1:3), par(1:4), nrm1(1:3, 1:3)
*
*-----------------------------------------------------------------------
*  define common blocs
*-----------------------------------------------------------------------
*                          common with gzlin2 subroutine
       REAL    a, b, czp, mia, c,d
       COMMON/gzl/a,b,czp,mia,c,d
*
*                          projection reference point (DC)
       REAL       dcprp(1:3)
*
*                          view plane distance (DC)
       REAL       dcdist
*
       COMMON/gzdc/dcprp,dcdist
*
       REAL    wdw(1:4),vpt(1:4),prp1(1:3),sc, dst
       COMMON/gzw/wdw,vpt,prp1,sc,dst
       SAVE /gzl/,/gzdc/,/gzw/
*
*
*-----------------------------------------------------------------------
*  declare local variables
*-----------------------------------------------------------------------
       REAL    xyz(1:3,1:3),nrm(1:3,1:3)
       REAL       temp(1:3)
*
       REAL        dy(1:3)
       REAL        x(1:3)
       REAL        y(1:3)
       REAL        n(1:3,1:3)
       REAL        sx(1:3)
       REAL        sn(1:3,1:3)
       REAL        yi
*
       INTEGER     edg, no,i,j,k
*-----------------------------------------------------------------------
*
       a=par(1)
       b=par(2)
       c=par(3)
       d=par(4)
       czp=par(3)*dcdist+ par(4)
       mia=a*dcprp(1) + b*dcprp(2) + c*dcprp(3) + par(4)
       IF(ABS(mia).LE.ROUNDOFF) RETURN
*
*---------------------------------------------   sort triangle vertices
       DO 10 i=1,3
*       DO 11 j=1,3
         xyz(1,i)=     xyz1(1,i)
         xyz(2,i)=xyz1(2,i)
         xyz(3,i)=xyz1(3,i)
         nrm(1,i)=nrm1(1,i)
         nrm(2,i)=nrm1(2,i)
         nrm(3,i)=nrm1(3,i)
 10    CONTINUE
       DO 100 k=2, 3
          DO 110 i=3, k, -1
            IF(xyz(2, i).GT.xyz(2, i-1)) THEN
               DO 120 j=1, 3
                   temp(j)=xyz(j, i-1)
                   xyz(j, i-1)=xyz(j, i)
                   xyz(j, i)=temp(j)
                   temp(j)=nrm(j, i-1)
                   nrm(j, i-1)=nrm(j, i)
                   nrm(j, i)=temp(j)
 120           CONTINUE
            ENDIF
 110      CONTINUE
 100   CONTINUE
*
*-----------------------------------------------------------------------
       y(1) =AINT(xyz(2, 1)) + 0.5
       y(2) =AINT(xyz(2, 2)) + 0.5
       y(3) =AINT(xyz(2, 3)) + 0.5
*
       dy(1) = y(1)-y(3)
       dy(2) = y(1)-y(2)
       dy(3) = y(2)-y(3)
*
       edg=0
       no=0
*
       IF(dy(1).GT.0.0) THEN
          sx(1)= (xyz(1,3)-xyz(1,1))/dy(1)
          sn(1,1)=(nrm(1,3)-nrm(1,1))/dy(1)
          sn(2,1)=(nrm(2,3)-nrm(2,1))/dy(1)
          sn(3,1)=(nrm(3,3)-nrm(3,1))/dy(1)
*
          x(1)=xyz(1,1) +0.5* sx(1)
*
          n(1,1)=nrm(1,1) + 0.5*sn(1,1)
          n(2,1)=nrm(2,1) + 0.5*sn(2,1)
          n(3,1)=nrm(3,1) + 0.5*sn(3,1)
*
          no=no+1
          edg=edg+2
       ENDIF
*
       IF(dy(2).GT.0.0) THEN
          sx(2)= (xyz(1,2)-xyz(1,1))/dy(2)
          sn(1,2)=(nrm(1,2)-nrm(1,1))/dy(2)
          sn(2,2)=(nrm(2,2)-nrm(2,1))/dy(2)
          sn(3,2)=(nrm(3,2)-nrm(3,1))/dy(2)
*
          x(2)=xyz(1,1) +0.5* sx(2)
*
          n(1,2)=nrm(1,1) + 0.5*sn(1,2)
          n(2,2)=nrm(2,1) + 0.5*sn(2,2)
          n(3,2)=nrm(3,1) + 0.5*sn(3,2)
*
          no=no+1
          edg=edg+4
       ENDIF
*
       IF(dy(3).GT.0.0) THEN
          sx(3)= (xyz(1,3)-xyz(1,2))/dy(3)
          sn(1,3)=(nrm(1,3)-nrm(1,2))/dy(3)
          sn(2,3)=(nrm(2,3)-nrm(2,2))/dy(3)
          sn(3,3)=(nrm(3,3)-nrm(3,2))/dy(3)
*
          x(3)=xyz(1,2) +0.5* sx(3)
*
          n(1,3)=nrm(1,2) + 0.5*sn(1,3)
          n(2,3)=nrm(2,2) + 0.5*sn(2,3)
          n(3,3)=nrm(3,2) + 0.5*sn(3,3)
*
          no=no+1
          edg=edg+8
       ENDIF
       IF(no.EQ.0) THEN
                CALL gzlin2(ANINT(xyz(2,1)),ANINT(xyz(1,1)),nrm(1,1),
     +                                      ANINT(xyz(1,2)),nrm(1,1))
       ELSE
           IF(no.EQ.2) THEN
              i=(edg-2)/4 + 1
           ELSE
              i=2
           ENDIF
           yi=AINT(y(1))
*
 50        IF(yi.GE.y(i)) THEN
                  IF(x(1).LT.x(i)) THEN
                  CALL gzlin2(yi,     (x(1)),n(1,1),
     +                                (x(i)),n(1,i))
                  ELSE
                  CALL gzlin2(yi,     (x(i)),n(1,i),
     +                              (x(1)),n(1,1))
                  ENDIF
             x(1)=x(1)+sx(1)
             n(1,1)=n(1,1)+sn(1,1)
             n(2,1)=n(2,1)+sn(2,1)
             n(3,1)=n(3,1)+sn(3,1)
*
             x(i)=x(i)+sx(i)
             n(1,i)=n(1,i)+sn(1,i)
             n(2,i)=n(2,i)+sn(2,i)
             n(3,i)=n(3,i)+sn(3,i)
*
             yi=yi-1.0
             GOTO 50
         ENDIF
*
*
             IF(no.EQ.3) THEN
*
 60            IF(yi.GE.y(3)) THEN
                  IF(x(1).LT.x(3)) THEN
                  CALL gzlin2(yi,     (x(1)),n(1,1),
     +                                (x(3)),n(1,3))
                  ELSE
                  CALL gzlin2(yi,     (x(3)),n(1,3),
     +                                (x(1)),n(1,1))
                  ENDIF
                 x(1)=x(1)+sx(1)
                 n(1,1)=n(1,1)+sn(1,1)
                 n(2,1)=n(2,1)+sn(2,1)
                 n(3,1)=n(3,1)+sn(3,1)
*
                 x(3)=x(3)+sx(3)
                 n(1,3)=n(1,3)+sn(1,3)
                 n(2,3)=n(2,3)+sn(2,3)
                 n(3,3)=n(3,3)+sn(3,3)
*
                 yi=yi-1.0
                 GOTO 60
              ENDIF
           ENDIF
*         ENDIF
         ENDIF
 643   RETURN
       END
*----------------------------------------------------------   end gzdrf1
************************************************************************
C@PROCESS OPT(3) IL(DIM) NOSDUMP NOGOSTMT DC(FELD)
       SUBROUTINE gzlin2(y, x1,n1, x2,n2)
************************************************************************
*
*  dummy arguments:
*      y - y coordinate of line segment;
*      x1,x2 - x coordinates of ends of line segment (x2>x1)
*      n1,s2 - illumination intensites if edpoints of line segment;
*
*  function: this subroutine draws a shaded horizontal line segment
*            from punkt (x1,y) to (x2,y); this subroutine is suitable
*            only for perspective projected lines;
*
************************************************************************
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
c      include(rdoff)
       REAL     ROUNDOFF
       PARAMETER (ROUNDOFF=1.0E-10)
*
       REAL       y
       REAL       x1
       REAL       n1(1:3)
       REAL       x2
       REAL       n2(1:3)
*-----------------------------------------------   define common blocs
       REAL       a, b, c, d
       REAL       czp
       REAL       mia
       COMMON/gzl/a, b, czp, mia, c,d
*
*                          projection reference point (DCS)
       REAL       dcprp(1:3)
*
*                          view plane distance (DCS)
       REAL       dcdist
       COMMON/gzdc/dcprp, dcdist
       REAL    wdw(1:4),vpt(1:4),prp1(1:3),sc, dst
       COMMON/gzw/wdw,vpt,prp1,sc,dst
       SAVE /gzl/,/gzdc/,/gzw/
*
*
       REAL        rgb1(1:3),dcsrcpos(1:3),dcsrcpwr,dcambnt
       INTEGER     stype,ptype,gloss
       COMMON/gzdcsrc/rgb1,dcsrcpos,dcsrcpwr,dcambnt,stype,ptype,gloss
       SAVE /gzdcsrc/
*
*
*-------------------------------------------   declare local variables
       REAL     rgb2(1:3)
       REAL     xyz(1:3)
       REAL       z
       REAL       sn(1:3)
       REAL       dn(1:3)
       REAL       mi
       REAL       dmi,ami
       REAL       x
       REAL       dx
       INTEGER    yi
*
*-----------------------------------------------------------------------
       mi=(a*x1 + b*y + czp)/mia
       ami=1.0-mi
       IF(ABS(ami).LE.ROUNDOFF) RETURN
       dx=x2-x1
       IF(dx.EQ.0.0) THEN
          xyz(1)=(x1-mi*dcprp(1))/ami
          xyz(2)=(y-mi*dcprp(2))/ami
          z=     (dcdist-mi*dcprp(3))/ami
          xyz(3)=z
          CALL gzshade(xyz,n1,rgb1,rgb2,gloss,dcprp,ptype,
     +                 dcsrcpos,stype,dcsrcpwr,dcambnt)
*
         CALL gzplt2(NINT(x1),MAXPKY - NINT(y) - 1, z, rgb2)
         RETURN
       ENDIF
       dmi=a/mia
       dn(1)=(n2(1)-n1(1))/dx
       dn(2)=(n2(2)-n1(2))/dx
       dn(3)=(n2(3)-n1(3))/dx
       sn(1)=n1(1)
       sn(2)=n1(2)
       sn(3)=n1(3)
       yi= MAXPKY - NINT(y) - 1
       x=x1
*
 310   IF(x.LE.x2) THEN
          ami=1.0-mi
          xyz(1)=(x-mi*dcprp(1))/ami
          xyz(2)=(y-mi*dcprp(2))/ami
          z=     (dcdist-mi*dcprp(3))/ami
          xyz(3)=z
          CALL gzshade(xyz,sn,rgb1,rgb2,gloss,dcprp,ptype,
     +                 dcsrcpos,stype,dcsrcpwr,dcambnt)
*
          CALL gzplt2(NINT(x),yi,z,rgb2)
*
          sn(1)=sn(1)+dn(1)
          sn(2)=sn(2)+dn(2)
          sn(3)=sn(3)+dn(3)
*
          mi=mi+dmi
          x=x+1.0
          GOTO 310
       ENDIF
*
       RETURN
*-------------------------------------------------------   end gzlin1
       END
*
*
*
*
*
*
*
*
*
*
*
*
*
*
*
*
************************************************************************
C@PROCESS OPT(3) IL(DIM) NOSDUMP NOGOSTMT DC(FELD)
       SUBROUTINE gzdrf3(xyz1,par,shd1)
************************************************************************
*
*  file name, type:  GZPLOT FORTRAN
*
*  dummy arguments: xyz - a triangle vertices coordinates;
*                   (DC after parallel projection);
*                   xyz(1,j),xyz(2,j),xyz(3,j) are x,y,z coordinates of
*                   vertex nr j of the triangle;
*
*                   par - triangle plane equation parameters
*                         (par(1)*x + par(2)*y + par(3)*z + par(4)=0);
*
*                   shd - illumination intensity of triangle vertices;
*                         shd(1,j),...,shd(3,j) are RGB components of
*                         illumination intensity of a vertex nr j;
*
*  function:  this subroutine fills a shaded triangle with Gouraud's
*             shading algorithm (the triangle was perspectivily
*             projected);
*
*  called by:
*
*  calls: gzlin1 (file GZPLOT FORTRAN)
*
************************************************************************
       IMPLICIT NONE
c      include(rdoff)
       REAL     ROUNDOFF
       PARAMETER (ROUNDOFF=1.0E-10)
*-----------------------------------------------------------------------
*  declare dummy arguments
*-----------------------------------------------------------------------
       REAL       xyz1(1:3, 1:3), par(1:4), shd1(1:3, 1:3)
*
*-----------------------------------------------------------------------
*  define common blocs
*-----------------------------------------------------------------------
*                          common with gzlin1 subroutine
       REAL    a, b, czp, mia, c,d
       COMMON/gzl/a,b,czp,mia,c,d
      save /gzl/
*
*
*
*
*-----------------------------------------------------------------------
*  declare local variables
*-----------------------------------------------------------------------
       REAL    xyz(1:3,1:3),shd(1:3,1:3)
       REAL       temp(1:3)
*
       REAL        dy(1:3)
       REAL        x(1:3)
       REAL        y(1:3)
       REAL        s(1:3,1:3)
       REAL        sx(1:3)
       REAL        ss(1:3,1:3)
       REAL        yi
*
       INTEGER     edg, no,i,j,k
*-----------------------------------------------------------------------
*
       a=par(1)
       b=par(2)
       c=par(3)
       d=par(4)
       IF(ABS(c).LE.ROUNDOFF) RETURN
*
*---------------------------------------------   sort triangle vertices
       DO 10 i=1,3
        DO 11 j=1,3
         xyz(j,i)=xyz1(j,i)
         shd(j,i)=shd1(j,i)
 11     CONTINUE
 10    CONTINUE
       DO 100 k=2, 3
          DO 110 i=3, k, -1
            IF(xyz(2, i).GT.xyz(2, i-1)) THEN
               DO 120 j=1, 3
                   temp(j)=xyz(j, i-1)
                   xyz(j, i-1)=xyz(j, i)
                   xyz(j, i)=temp(j)
                   temp(j)=shd(j, i-1)
                   shd(j, i-1)=shd(j, i)
                   shd(j, i)=temp(j)
 120           CONTINUE
            ENDIF
 110      CONTINUE
 100   CONTINUE
*
*-----------------------------------------------------------------------
       y(1) =AINT(xyz(2, 1)) + 0.5
       y(2) =AINT(xyz(2, 2)) + 0.5
       y(3) =AINT(xyz(2, 3)) + 0.5
*
       dy(1) = y(1)-y(3)
       dy(2) = y(1)-y(2)
       dy(3) = y(2)-y(3)
*
       edg=0
       no=0
*
       IF(dy(1).GT.0.0) THEN
          sx(1)= (xyz(1,3)-xyz(1,1))/dy(1)
          ss(1,1)=(shd(1,3)-shd(1,1))/dy(1)
          ss(2,1)=(shd(2,3)-shd(2,1))/dy(1)
          ss(3,1)=(shd(3,3)-shd(3,1))/dy(1)
*
          x(1)=xyz(1,1) + 0.5*sx(1)
*
          s(1,1)=shd(1,1) + 0.5*ss(1,1)
          s(2,1)=shd(2,1) + 0.5*ss(2,1)
          s(3,1)=shd(3,1) + 0.5*ss(3,1)
*
          no=no+1
          edg=edg+2
       ENDIF
*
       IF(dy(2).GT.0.0) THEN
          sx(2)= (xyz(1,2)-xyz(1,1))/dy(2)
          ss(1,2)=(shd(1,2)-shd(1,1))/dy(2)
          ss(2,2)=(shd(2,2)-shd(2,1))/dy(2)
          ss(3,2)=(shd(3,2)-shd(3,1))/dy(2)
*
          x(2)=xyz(1,1) + 0.5*sx(2)
*
          s(1,2)=shd(1,1) + 0.5*ss(1,2)
          s(2,2)=shd(2,1) + 0.5*ss(2,2)
          s(3,2)=shd(3,1) + 0.5*ss(3,2)
*
          no=no+1
          edg=edg+4
       ENDIF
*
       IF(dy(3).GT.0.0) THEN
          sx(3)= (xyz(1,3)-xyz(1,2))/dy(3)
          ss(1,3)=(shd(1,3)-shd(1,2))/dy(3)
          ss(2,3)=(shd(2,3)-shd(2,2))/dy(3)
          ss(3,3)=(shd(3,3)-shd(3,2))/dy(3)
*
          x(3)=xyz(1,2) + 0.5*sx(3)
*
          s(1,3)=shd(1,2) + 0.5*ss(1,3)
          s(2,3)=shd(2,2) + 0.5*ss(2,3)
          s(3,3)=shd(3,2) + 0.5*ss(3,3)
*
          no=no+1
          edg=edg+8
       ENDIF
       IF(no.EQ.0) THEN
                   CALL gzplt1(NINT(xyz(1,1)), NINT(xyz(2,1)),
     +             (xyz(3,1)+xyz(3,2)+xyz(3,3))/3.0, shd(1,1))
       ELSE
           IF(no.EQ.2) THEN
              i=(edg-2)/4 + 1
           ELSE
              i=2
           ENDIF
           yi=AINT(y(1))
 50        IF(yi.GE.y(i)) THEN
                  IF(x(1).LT.x(i)) THEN
                  CALL gzlin3(yi,ANINT(x(1)),s(1,1),
     +                           ANINT(x(i)),s(1,i))
                  ELSE
                  CALL gzlin3(yi,ANINT(x(i)),s(1,i),
     +                         ANINT(x(1)),s(1,1))
                  ENDIF
             x(1)=x(1)+sx(1)
             s(1,1)=s(1,1)+ss(1,1)
             s(2,1)=s(2,1)+ss(2,1)
             s(3,1)=s(3,1)+ss(3,1)
*
             x(i)=x(i)+sx(i)
             s(1,i)=s(1,i)+ss(1,i)
             s(2,i)=s(2,i)+ss(2,i)
             s(3,i)=s(3,i)+ss(3,i)
*
             yi=yi-1.0
             GOTO 50
         ENDIF
*
*
             IF(no.EQ.3) THEN
*
 60            IF(yi.GE.y(3)) THEN
                  IF(x(1).LT.x(3)) THEN
                  CALL gzlin3(yi,ANINT(x(1)),s(1,1),
     +                           ANINT(x(3)),s(1,3))
                  ELSE
                  CALL gzlin3(yi,ANINT(x(3)),s(1,3),
     +                           ANINT(x(1)),s(1,1))
                  ENDIF
                 x(1)=x(1)+sx(1)
                 s(1,1)=s(1,1)+ss(1,1)
                 s(2,1)=s(2,1)+ss(2,1)
                 s(3,1)=s(3,1)+ss(3,1)
*
                 x(3)=x(3)+sx(3)
                 s(1,3)=s(1,3)+ss(1,3)
                 s(2,3)=s(2,3)+ss(2,3)
                 s(3,3)=s(3,3)+ss(3,3)
*
                 yi=yi-1.0
                 GOTO 60
              ENDIF
           ENDIF
*         ENDIF
         ENDIF
 643   RETURN
       END
*----------------------------------------------------------   end gzdrf3
************************************************************************
************************************************************************
C@PROCESS OPT(3) IL(DIM) NOSDUMP NOGOSTMT DC(FELD)
       SUBROUTINE gzlin3(y, x1,s1, x2,s2)
************************************************************************
*
*  dummy arguments:
*      y - y coordinate of line segment;
*      x1,x2 - x coordinates of ends of line segment (x2>x1)
*      s1,s2 - illumination intensites if edpoints of line segment;
*
*  function: this subroutine draws a shaded horizontal line segment
*            from punkt (x1,y) to (x2,y); this subroutine is suitable
*            only for perspective projected lines;
*
************************************************************************
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
*
       REAL       y
       REAL       x1
       REAL       s1(1:3)
       REAL       x2
       REAL       s2(1:3)
*-----------------------------------------------   define common blocs
       REAL       a, b, c, d
       REAL       czp
       REAL       mia
       COMMON/gzl/a, b, czp, mia, c,d
       SAVE /gzl/
*
*
*-------------------------------------------   declare local variables
       REAL       z
       REAL       si(1:3)
       REAL       ds(1:3)
       REAL       dz
       REAL       x
       REAL       dx
       INTEGER    yi
*
*-----------------------------------------------------------------------
       dx=x2-x1
       IF(dx.EQ.0.0) THEN
         CALL gzplt1(NINT(x1),MAXPKY - NINT(y) - 1,
     +   (-a*x1 - b*y - d)/c, s1)
         RETURN
       ENDIF
*
       z=(-a*x1 - b*y - d)/c
       dz=-a/c
       ds(1)=(s2(1)-s1(1))/dx
       ds(2)=(s2(2)-s1(2))/dx
       ds(3)=(s2(3)-s1(3))/dx
       si(1)=s1(1)
       si(2)=s1(2)
       si(3)=s1(3)
       yi= MAXPKY - NINT(y) - 1
       x=x1
*
 310   IF(x.LE.x2) THEN
          CALL gzplt1(NINT(x), yi, z, si)
*
          si(1)=si(1)+ds(1)
          si(2)=si(2)+ds(2)
          si(3)=si(3)+ds(3)
*
          z=z+dz
          x=x+1.0
          GOTO 310
       ENDIF
*
       RETURN
*-------------------------------------------------------   end gzlin3
       END
*
************************************************************************
C@PROCESS OPT(3) IL(DIM) NOSDUMP NOGOSTMT DC(FELD)
       SUBROUTINE gzdda2(x1,y1, x2,y2, par, color)
*
************************************************************************
*
*  dummy arguments:
*      x1,y1  - coordinates of begin point of line segment (DC after
*               parallel projection);
*      x2,y2  - coordinates of end point of line segment;
*      par    - line segment plane parameters
*               (par(1)*xi + par(2)*yi + par(3)*zi + par(4)=0, where
*               xi,yi,zi belong to the line segment;
*
*  function:  this subroutine draws a perspective projected line segment
*             in z-buffer with color 'color';
*
************************************************************************
       IMPLICIT NONE
*
c      include (gzcst11)
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
c      include (gzcst21)
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
c      include(rdoff)
       REAL     ROUNDOFF
       PARAMETER (ROUNDOFF=1.0E-10)
*-------------------------------------------   declare dummy arguments
       REAL       x1,y1, x2,y2, par(1:4)
       INTEGER    color
*
*-------------------------------------------   declare local variables
*
       REAL       a,b,c,d
       REAL       yy1,yy2, xx1,xx2, dx,dy, absdx,absdy, sx,sy, x,y
       REAL       rob, temp, z
       INTEGER    i, steps
*
*-----------------------------------------------------------------------
*
       a=par(1)
       b=par(2)
       c=par(3)
       d=par(4)
       IF(ABS(c).LT.ROUNDOFF) RETURN
*
       yy1=AINT(y1) + 0.5
       yy2=AINT(y2) + 0.5
*
       xx1=x1
       xx2=x2
*
       IF(yy1.LT.yy2) THEN
         rob=yy1
         yy1=yy2
         yy2=rob
         rob=xx1
         xx1=xx2
         xx2=rob
       ENDIF
*
       dx=xx2-xx1
       dy=yy2-yy1
*
       absdx=ABS(dx)
       absdy=ABS(dy)
*
       IF(absdx.LT.1.0 .AND. absdy.LT.1.0) THEN

        CALL gzplt5(NINT(x1),-NINT(y1)+MAXPKY-1,
     +  (-a*x1- b*y1- d)/c , color)
         RETURN
       ENDIF
*
       temp=ANINT(AMAX1(absdx,absdy))
       sx=dx/temp
       sy=dy/temp
       x=xx1 + 0.5*sx
       y=AINT(yy1)
*
       z=(-a*x - b*y -d)/c
*       dz=-a/c
       steps=NINT(temp)
*
       DO 100 i=1,steps
*
         CALL gzplt5(NINT(x),-NINT(y)+MAXPKY-1,
     +    z+.2,color)
*
         x=x+sx
         y=y+sy
         z=(-a*x - b*y -d)/c
*
 100   CONTINUE
*
*-----------------------------------------------------   end gzdda1
       END
************************************************************************
C@PROCESS OPT(3) IL(DIM) NOSDUMP NOGOSTMT DC(FELD)
       SUBROUTINE gzdrf4(xyz1,par,nrm1)
************************************************************************
*
*  file name, type:  GZPLOT FORTRAN
*
*  dummy arguments: xyz - a triangle vertices coordinates;
*                   (DC after parallel projection);
*                   xyz(1,j),xyz(2,j),xyz(3,j) are x,y,z coordinates of
*                   vertex nr j of the triangle;
*
*                   par - triangle plane equation parameters
*                         (par(1)*x + par(2)*y + par(3)*z + par(4)=0);
*
*                   nrm - vetex normals;
*
*  function:  this subroutine fills a shaded triangle with Phong's
*             shading algorithm (the triangle was perspectivily
*             projected);
*
*  called by:
*
*  calls: gzlin2 (file GZPLOT FORTRAN)
*
************************************************************************
       IMPLICIT NONE
c      include(rdoff)
       REAL     ROUNDOFF
       PARAMETER (ROUNDOFF=1.0E-10)
*-----------------------------------------------------------------------
*  declare dummy arguments
*-----------------------------------------------------------------------
       REAL       xyz1(1:3, 1:3), par(1:4), nrm1(1:3, 1:3)
*
*-----------------------------------------------------------------------
*  define common blocs
*-----------------------------------------------------------------------
*                          common with gzlin2 subroutine
       REAL    a, b, czp, mia, c,d
       COMMON/gzl/a,b,czp,mia,c,d
*
*
       REAL    wdw(1:4),vpt(1:4),prp1(1:3),sc, dst
       COMMON/gzw/wdw,vpt,prp1,sc,dst
       SAVE /gzl/,/gzw/
*
*
*-----------------------------------------------------------------------
*  declare local variables
*-----------------------------------------------------------------------
       REAL    xyz(1:3,1:3),nrm(1:3,1:3)
       REAL       temp(1:3)
*
       REAL        dy(1:3)
       REAL        x(1:3)
       REAL        y(1:3)
       REAL        n(1:3,1:3)
       REAL        sx(1:3)
       REAL        sn(1:3,1:3)
       REAL        yi
*
       INTEGER     edg, no,i,j,k
*-----------------------------------------------------------------------
*
       a=par(1)
       b=par(2)
       c=par(3)
       d=par(4)
       IF(ABS(c).LE.ROUNDOFF) RETURN
*
*---------------------------------------------   sort triangle vertices
       DO 10 i=1,3
*       DO 11 j=1,3
         xyz(1,i)=     xyz1(1,i)
         xyz(2,i)=xyz1(2,i)
         xyz(3,i)=xyz1(3,i)
         nrm(1,i)=nrm1(1,i)
         nrm(2,i)=nrm1(2,i)
         nrm(3,i)=nrm1(3,i)
 10    CONTINUE
       DO 100 k=2, 3
          DO 110 i=3, k, -1
            IF(xyz(2, i).GT.xyz(2, i-1)) THEN
               DO 120 j=1, 3
                   temp(j)=xyz(j, i-1)
                   xyz(j, i-1)=xyz(j, i)
                   xyz(j, i)=temp(j)
                   temp(j)=nrm(j, i-1)
                   nrm(j, i-1)=nrm(j, i)
                   nrm(j, i)=temp(j)
 120           CONTINUE
            ENDIF
 110      CONTINUE
 100   CONTINUE
*
*-----------------------------------------------------------------------
       y(1) =AINT(xyz(2, 1)) + 0.5
       y(2) =AINT(xyz(2, 2)) + 0.5
       y(3) =AINT(xyz(2, 3)) + 0.5
*
       dy(1) = y(1)-y(3)
       dy(2) = y(1)-y(2)
       dy(3) = y(2)-y(3)
*
       edg=0
       no=0
*
       IF(dy(1).GT.0.0) THEN
          sx(1)= (xyz(1,3)-xyz(1,1))/dy(1)
          sn(1,1)=(nrm(1,3)-nrm(1,1))/dy(1)
          sn(2,1)=(nrm(2,3)-nrm(2,1))/dy(1)
          sn(3,1)=(nrm(3,3)-nrm(3,1))/dy(1)
*
          x(1)=xyz(1,1) +0.5* sx(1)
*
          n(1,1)=nrm(1,1) + 0.5*sn(1,1)
          n(2,1)=nrm(2,1) + 0.5*sn(2,1)
          n(3,1)=nrm(3,1) + 0.5*sn(3,1)
*
          no=no+1
          edg=edg+2
       ENDIF
*
       IF(dy(2).GT.0.0) THEN
          sx(2)= (xyz(1,2)-xyz(1,1))/dy(2)
          sn(1,2)=(nrm(1,2)-nrm(1,1))/dy(2)
          sn(2,2)=(nrm(2,2)-nrm(2,1))/dy(2)
          sn(3,2)=(nrm(3,2)-nrm(3,1))/dy(2)
*
          x(2)=xyz(1,1) +0.5* sx(2)
*
          n(1,2)=nrm(1,1) + 0.5*sn(1,2)
          n(2,2)=nrm(2,1) + 0.5*sn(2,2)
          n(3,2)=nrm(3,1) + 0.5*sn(3,2)
*
          no=no+1
          edg=edg+4
       ENDIF
*
       IF(dy(3).GT.0.0) THEN
          sx(3)= (xyz(1,3)-xyz(1,2))/dy(3)
          sn(1,3)=(nrm(1,3)-nrm(1,2))/dy(3)
          sn(2,3)=(nrm(2,3)-nrm(2,2))/dy(3)
          sn(3,3)=(nrm(3,3)-nrm(3,2))/dy(3)
*
          x(3)=xyz(1,2) +0.5* sx(3)
*
          n(1,3)=nrm(1,2) + 0.5*sn(1,3)
          n(2,3)=nrm(2,2) + 0.5*sn(2,3)
          n(3,3)=nrm(3,2) + 0.5*sn(3,3)
*
          no=no+1
          edg=edg+8
       ENDIF
       IF(no.EQ.0) THEN
                CALL gzlin4(ANINT(xyz(2,1)),ANINT(xyz(1,1)),nrm(1,1),
     +                                      ANINT(xyz(1,2)),nrm(1,1))
       ELSE
           IF(no.EQ.2) THEN
              i=(edg-2)/4 + 1
           ELSE
              i=2
           ENDIF
           yi=AINT(y(1))
*
 50        IF(yi.GE.y(i)) THEN
                  IF(x(1).LT.x(i)) THEN
                  CALL gzlin4(yi,ANINT(x(1)),n(1,1),
     +                           ANINT(x(i)),n(1,i))
                  ELSE
                  CALL gzlin4(yi,ANINT(x(i)),n(1,i),
     +                          ANINT(x(1)),n(1,1))
                  ENDIF
             x(1)=x(1)+sx(1)
             n(1,1)=n(1,1)+sn(1,1)
             n(2,1)=n(2,1)+sn(2,1)
             n(3,1)=n(3,1)+sn(3,1)
*
             x(i)=x(i)+sx(i)
             n(1,i)=n(1,i)+sn(1,i)
             n(2,i)=n(2,i)+sn(2,i)
             n(3,i)=n(3,i)+sn(3,i)
*
             yi=yi-1.0
             GOTO 50
         ENDIF
*
*
             IF(no.EQ.3) THEN
*
 60            IF(yi.GE.y(3)) THEN
                  IF(x(1).LT.x(3)) THEN
                  CALL gzlin4(yi,ANINT(x(1)),n(1,1),
     +                           ANINT(x(3)),n(1,3))
                  ELSE
                  CALL gzlin4(yi,ANINT(x(3)),n(1,3),
     +                           ANINT(x(1)),n(1,1))
                  ENDIF
                 x(1)=x(1)+sx(1)
                 n(1,1)=n(1,1)+sn(1,1)
                 n(2,1)=n(2,1)+sn(2,1)
                 n(3,1)=n(3,1)+sn(3,1)
*
                 x(3)=x(3)+sx(3)
                 n(1,3)=n(1,3)+sn(1,3)
                 n(2,3)=n(2,3)+sn(2,3)
                 n(3,3)=n(3,3)+sn(3,3)
*
                 yi=yi-1.0
                 GOTO 60
              ENDIF
           ENDIF
*         ENDIF
         ENDIF
 643   RETURN
       END
*----------------------------------------------------------   end gzdrf4
************************************************************************
C@PROCESS OPT(3) IL(DIM) NOSDUMP NOGOSTMT DC(FELD)
       SUBROUTINE gzlin4(y, x1,n1, x2,n2)
************************************************************************
*
*  dummy arguments:
*      y - y coordinate of line segment;
*      x1,x2 - x coordinates of ends of line segment (x2>x1)
*      n1,s2 - illumination intensites if edpoints of line segment;
*
*  function: this subroutine draws a shaded horizontal line segment
*            from punkt (x1,y) to (x2,y); this subroutine is suitable
*            only for parallel projected lines;
*
************************************************************************
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
*
       REAL       y
       REAL       x1
       REAL       n1(1:3)
       REAL       x2
       REAL       n2(1:3)
*-----------------------------------------------   define common blocs
       REAL       a, b, c, d
       REAL       czp
       REAL       mia
       COMMON/gzl/a, b, czp, mia, c,d
*
*                          projection reference point (DC)
       REAL       dcprp(1:3)
*
*                          view plane distance (DC)
       REAL       dcdist
*
       COMMON/gzdc/dcprp,dcdist
       SAVE /gzl/,/gzdc/
*
*
*
       REAL        rgb1(1:3),dcsrcpos(1:3),dcsrcpwr,dcambnt
       INTEGER     stype,ptype,gloss
       COMMON/gzdcsrc/rgb1,dcsrcpos,dcsrcpwr,dcambnt,stype,ptype,gloss
       SAVE /gzdcsrc/
*
*
*-------------------------------------------   declare local variables
       REAL     rgb2(1:3)
       REAL     xyz(1:3)
       REAL       z, dz
       REAL       sn(1:3)
       REAL       dn(1:3)
       REAL       x
       REAL       dx
       INTEGER    yi
*
*-----------------------------------------------------------------------
       dx=x2-x1
       IF(dx.EQ.0.0) THEN
          xyz(1)=x1
          xyz(2)=y
          z=     (-a*x1 - b*y -d)/c
          xyz(3)=z
          CALL gzshade(xyz,n1,rgb1,rgb2,gloss,dcprp,ptype,
     +                 dcsrcpos,stype,dcsrcpwr,dcambnt)
*
         CALL gzplt2(NINT(x1),MAXPKY - NINT(y) - 1, z, rgb2)
         RETURN
       ENDIF
*
       z=(-a*x1 - b*y - d)/c
       dz=-a/c
       dn(1)=(n2(1)-n1(1))/dx
       dn(2)=(n2(2)-n1(2))/dx
       dn(3)=(n2(3)-n1(3))/dx
       sn(1)=n1(1)
       sn(2)=n1(2)
       sn(3)=n1(3)
       yi= MAXPKY - NINT(y) - 1
       x=x1
       xyz(2)=y
 310   IF(x.LE.x2) THEN
          xyz(1)=x
          xyz(3)=z
          CALL gzshade(xyz,sn,rgb1,rgb2,gloss,dcprp,ptype,
     +                 dcsrcpos,stype,dcsrcpwr,dcambnt)
*
          CALL gzplt2(NINT(x),yi,z,rgb2)
*
          sn(1)=sn(1)+dn(1)
          sn(2)=sn(2)+dn(2)
          sn(3)=sn(3)+dn(3)
*
          z=z+dz
          x=x+1.0
          GOTO 310
       ENDIF
*
       RETURN
*-------------------------------------------------------   end gzlin4
       END
*
*
*
*
*
*
*
*
*
************************************************************************
C@PROCESS OPT(3) IL(DIM) NOSDUMP NOGOSTMT DC(FELD)
       SUBROUTINE gzdrf5(xyz1,par,shd)
************************************************************************
*
*  file name, type:  GZPLOT FORTRAN
*
*  dummy arguments: xyz - a triangle vertices coordinates;
*                   (DC after perspective projection);
*                   xyz(1,j),xyz(2,j),xyz(3,j) are x,y,z coordinates of
*                   vertex nr j of the triangle;
*
*                   par - triangle plane equation parameters
*                         (par(1)*x + par(2)*y + par(3)*z + par(4)=0);
*
*                   shd - illumination intensity of a triangle;
*                         shd(1),...,shd(3) are RGB components of
*                         illumination intensity of a triangle;
*
*  function:  this subroutine fills a shaded triangle with constant
*             shading algorithm (the triangle was perspectivily
*             projected);
*
*  called by:
*
*  calls: gzlin5 (file GZPLOT FORTRAN)
*
************************************************************************
       IMPLICIT NONE
c      include(rdoff)
       REAL     ROUNDOFF
       PARAMETER (ROUNDOFF=1.0E-10)
*-----------------------------------------------------------------------
*  declare dummy arguments
*-----------------------------------------------------------------------
       REAL       xyz1(1:3, 1:3), par(1:4), shd(1:3)
*
*-----------------------------------------------------------------------
*  define common blocs
*-----------------------------------------------------------------------
*                          common with gzlin1 subroutine
       REAL    a, b, czp, mia, c,d
       COMMON/gzl/a,b,czp,mia,c,d
*
*                          projection reference point (DC)
       REAL       dcprp(1:3)
*
*                          view plane distance (DC)
       REAL       dcdist
*
       COMMON/gzdc/dcprp,dcdist
*
       REAL    wdw(1:4),vpt(1:4),prp1(1:3),sc, dst
       COMMON/gzw/wdw,vpt,prp1,sc,dst
       SAVE /gzl/,/gzdc/,/gzw/
*
*
*-----------------------------------------------------------------------
*  declare local variables
*-----------------------------------------------------------------------
       REAL    xyz(1:3,1:3)
       REAL       temp(1:3)
*
       REAL        dy(1:3)
       REAL        x(1:3)
       REAL        y(1:3)
       REAL        sx(1:3)
       REAL        yi
*
       INTEGER     edg, no,i,j,k
*-----------------------------------------------------------------------
*
       a=par(1)
       b=par(2)
       c=par(3)
       d=par(4)
       czp=par(3)*dcdist+ par(4)
       mia=a*dcprp(1) + b*dcprp(2) + c*dcprp(3) + par(4)
       IF(ABS(mia).LE.ROUNDOFF) RETURN
*
*---------------------------------------------   sort triangle vertices
       DO 10 i=1,3
        DO 11 j=1,3
         xyz(j,i)=xyz1(j,i)
 11     CONTINUE
 10    CONTINUE
       DO 100 k=2, 3
          DO 110 i=3, k, -1
            IF(xyz(2, i).GT.xyz(2, i-1)) THEN
               DO 120 j=1, 3
                   temp(j)=xyz(j, i-1)
                   xyz(j, i-1)=xyz(j, i)
                   xyz(j, i)=temp(j)
 120           CONTINUE
            ENDIF
 110      CONTINUE
 100   CONTINUE
*
*-----------------------------------------------------------------------
       y(1) =AINT(xyz(2, 1)) + 0.5
       y(2) =AINT(xyz(2, 2)) + 0.5
       y(3) =AINT(xyz(2, 3)) + 0.5
*
       dy(1) = y(1)-y(3)
       dy(2) = y(1)-y(2)
       dy(3) = y(2)-y(3)
*
       edg=0
       no=0
*
       IF(dy(1).GT.0.0) THEN
          sx(1)= (xyz(1,3)-xyz(1,1))/dy(1)
*
          x(1)=xyz(1,1) + 0.5*sx(1)
*
*
          no=no+1
          edg=edg+2
       ENDIF
*
       IF(dy(2).GT.0.0) THEN
          sx(2)= (xyz(1,2)-xyz(1,1))/dy(2)
*
          x(2)=xyz(1,1) + 0.5*sx(2)
*
          no=no+1
          edg=edg+4
       ENDIF
*
       IF(dy(3).GT.0.0) THEN
          sx(3)= (xyz(1,3)-xyz(1,2))/dy(3)
*
          x(3)=xyz(1,2) + 0.5*sx(3)
*
*
          no=no+1
          edg=edg+8
       ENDIF
       IF(no.EQ.0) THEN
                   CALL gzplt1(NINT(xyz(1,1)), NINT(xyz(2,1)),
     +             (xyz(3,1)+xyz(3,2)+xyz(3,3))/3.0, shd)
       ELSE
           IF(no.EQ.2) THEN
              i=(edg-2)/4 + 1
           ELSE
              i=2
           ENDIF
           yi=AINT(y(1))
 50        IF(yi.GE.y(i)) THEN
                  IF(x(1).LT.x(i)) THEN
                  CALL gzlin5(yi,ANINT(x(1)), ANINT(x(i)),shd)
                  ELSE
                  CALL gzlin5(yi,ANINT(x(i)),ANINT(x(1)),shd)
                  ENDIF
             x(1)=x(1)+sx(1)
*
             x(i)=x(i)+sx(i)
*
             yi=yi-1.0
             GOTO 50
         ENDIF
*
*
             IF(no.EQ.3) THEN
*
 60            IF(yi.GE.y(3)) THEN
                  IF(x(1).LT.x(3)) THEN
                  CALL gzlin5(yi,ANINT(x(1)), ANINT(x(3)),shd)
                  ELSE
                  CALL gzlin5(yi,ANINT(x(3)), ANINT(x(1)),shd)
                  ENDIF
                 x(1)=x(1)+sx(1)
*
                 x(3)=x(3)+sx(3)
*
                 yi=yi-1.0
                 GOTO 60
              ENDIF
           ENDIF
*         ENDIF
         ENDIF
 643   RETURN
       END
*----------------------------------------------------------   end gzdrf5
************************************************************************
C@PROCESS OPT(3) IL(DIM) NOSDUMP NOGOSTMT DC(FELD)
       SUBROUTINE gzlin5(y, x1, x2, s)
************************************************************************
*
*  dummy arguments:
*      y - y coordinate of line segment;
*      x1,x2 - x coordinates of ends of line segment (x2>x1)
*      s     - illumination intensity of a line segment;
*
*  function: this subroutine draws a shaded horizontal line segment
*            from punkt (x1,y) to (x2,y); this subroutine is suitable
*            only for perspective projected lines;
*
************************************************************************
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
c      include(rdoff)
       REAL     ROUNDOFF
       PARAMETER (ROUNDOFF=1.0E-10)
*
       REAL       y
       REAL       x1
       REAL       x2
       REAL       s(1:3)
*-----------------------------------------------   define common blocs
       REAL       a, b, c, d
       REAL       czp
       REAL       mia
       COMMON/gzl/a, b, czp, mia, c,d
*
*                          projection reference point (DCS)
       REAL       dcprp(1:3)
*
*                          view plane distance (DCS)
       REAL       dcdist
       COMMON/gzdc/dcprp, dcdist
       REAL    wdw(1:4),vpt(1:4),prp1(1:3),sc, dst
       COMMON/gzw/wdw,vpt,prp1,sc,dst
        SAVE /gzl/,/gzdc/,/gzw/
*
*-------------------------------------------   declare local variables
       REAL       z
       REAL       mi
       REAL       dmi
       REAL       x
       REAL       dx
       INTEGER    yi
*
*-----------------------------------------------------------------------
       mi=(a*x1 + b*y + czp)/mia
       IF(ABS(1.0-mi).LE.ROUNDOFF) RETURN
       dx=x2-x1
       IF(dx.EQ.0.0) THEN
         CALL gzplt1(NINT(x1),MAXPKY - NINT(y) - 1,
     +   (dcdist-mi*dcprp(3))/(1.0-mi),s)
         RETURN
       ENDIF
       dmi=a/mia
       yi= MAXPKY - NINT(y) - 1
       x=x1

 310   IF(x.LE.x2) THEN
          z = (dcdist-mi*dcprp(3))/(1.0-mi)
*         z = (-a*x - b*y - d)/c
          CALL gzplt1(NINT(x), yi, z, s)
*
          mi=mi+dmi
          x=x+1.0
          GOTO 310
       ENDIF
*
       RETURN
*-------------------------------------------------------   end gzlin5
       END
************************************************************************
C@PROCESS OPT(3) IL(DIM) NOSDUMP NOGOSTMT DC(FELD)
       SUBROUTINE gzdrf6(xyz1,par,shd)
************************************************************************
*
*  file name, type:  GZPLOT FORTRAN
*
*  dummy arguments: xyz - a triangle vertices coordinates;
*                   (DC after parallel projection);
*                   xyz(1,j),xyz(2,j),xyz(3,j) are x,y,z coordinates of
*                   vertex nr j of the triangle;
*
*                   par - triangle plane equation parameters
*                         (par(1)*x + par(2)*y + par(3)*z + par(4)=0);
*
*                   shd - illumination intensity of a triangle;
*                         shd(1),...,shd(3) are RGB components of
*                         illumination intensity of a triangle;
*
*  function:  this subroutine fills a shaded triangle with constant
*             shading algorithm (the triangle was parallel
*             projected);
*
*  called by:
*
*  calls: gzlin6 (file GZPLOT FORTRAN)
*
************************************************************************
       IMPLICIT NONE
c      include(rdoff)
       REAL     ROUNDOFF
       PARAMETER (ROUNDOFF=1.0E-10)
*-----------------------------------------------------------------------
*  declare dummy arguments
*-----------------------------------------------------------------------
       REAL       xyz1(1:3, 1:3), par(1:4), shd(1:3)
*
*-----------------------------------------------------------------------
*  define common blocs
*-----------------------------------------------------------------------
*                          common with gzlin1 subroutine
       REAL    a, b, czp, mia, c,d
       COMMON/gzl/a,b,czp,mia,c,d
       save /gzl/
*
*
*
*
*-----------------------------------------------------------------------
*  declare local variables
*-----------------------------------------------------------------------
       REAL    xyz(1:3,1:3)
       REAL       temp(1:3)
*
       REAL        dy(1:3)
       REAL        x(1:3)
       REAL        y(1:3)
       REAL        sx(1:3)
       REAL        yi
*
       INTEGER     edg, no,i,j,k
*-----------------------------------------------------------------------
*
       a=par(1)
       b=par(2)
       c=par(3)
       d=par(4)
       IF(ABS(c).LE.ROUNDOFF) RETURN
*
*---------------------------------------------   sort triangle vertices
       DO 10 i=1,3
        DO 11 j=1,3
         xyz(j,i)=xyz1(j,i)
 11     CONTINUE
 10    CONTINUE
       DO 100 k=2, 3
          DO 110 i=3, k, -1
            IF(xyz(2, i).GT.xyz(2, i-1)) THEN
               DO 120 j=1, 3
                   temp(j)=xyz(j, i-1)
                   xyz(j, i-1)=xyz(j, i)
                   xyz(j, i)=temp(j)
 120           CONTINUE
            ENDIF
 110      CONTINUE
 100   CONTINUE
*
*-----------------------------------------------------------------------
       y(1) =AINT(xyz(2, 1)) + 0.5
       y(2) =AINT(xyz(2, 2)) + 0.5
       y(3) =AINT(xyz(2, 3)) + 0.5
*
       dy(1) = y(1)-y(3)
       dy(2) = y(1)-y(2)
       dy(3) = y(2)-y(3)
*
       edg=0
       no=0
*
       IF(dy(1).GT.0.0) THEN
          sx(1)= (xyz(1,3)-xyz(1,1))/dy(1)
*
          x(1)=xyz(1,1) + 0.5*sx(1)
*
          no=no+1
          edg=edg+2
       ENDIF
*
       IF(dy(2).GT.0.0) THEN
          sx(2)= (xyz(1,2)-xyz(1,1))/dy(2)
*
          x(2)=xyz(1,1) + 0.5*sx(2)
*
*
          no=no+1
          edg=edg+4
       ENDIF
*
       IF(dy(3).GT.0.0) THEN
          sx(3)= (xyz(1,3)-xyz(1,2))/dy(3)
*
          x(3)=xyz(1,2) + 0.5*sx(3)
*
*
          no=no+1
          edg=edg+8
       ENDIF
       IF(no.EQ.0) THEN
                   CALL gzplt1(NINT(xyz(1,1)), NINT(xyz(2,1)),
     +             (xyz(3,1)+xyz(3,2)+xyz(3,3))/3.0, shd)
       ELSE
           IF(no.EQ.2) THEN
              i=(edg-2)/4 + 1
           ELSE
              i=2
           ENDIF
           yi=AINT(y(1))
 50        IF(yi.GE.y(i)) THEN
                  IF(x(1).LT.x(i)) THEN
                  CALL gzlin6(yi,ANINT(x(1)), ANINT(x(i)),shd)
                  ELSE
                  CALL gzlin6(yi,ANINT(x(i)), ANINT(x(1)),shd)
                  ENDIF
             x(1)=x(1)+sx(1)
*
             x(i)=x(i)+sx(i)
*
             yi=yi-1.0
             GOTO 50
         ENDIF
*
*
             IF(no.EQ.3) THEN
*
 60            IF(yi.GE.y(3)) THEN
                  IF(x(1).LT.x(3)) THEN
                  CALL gzlin6(yi,ANINT(x(1)), ANINT(x(3)),shd)
                  ELSE
                  CALL gzlin6(yi,ANINT(x(3)), ANINT(x(1)),shd)
                  ENDIF
                 x(1)=x(1)+sx(1)
*
                 x(3)=x(3)+sx(3)
*
                 yi=yi-1.0
                 GOTO 60
              ENDIF
           ENDIF
*         ENDIF
         ENDIF
 643   RETURN
       END
*----------------------------------------------------------   end gzdrf3
************************************************************************
************************************************************************
C@PROCESS OPT(3) IL(DIM) NOSDUMP NOGOSTMT DC(FELD)
       SUBROUTINE gzlin6(y, x1, x2,s)
************************************************************************
*
*  dummy arguments:
*      y - y coordinate of line segment;
*      x1,x2 - x coordinates of ends of line segment (x2>x1)
*      s     - illumination intensity if a line segment;
*
*  function: this subroutine draws a shaded horizontal line segment
*            from punkt (x1,y) to (x2,y); this subroutine is suitable
*            only for parallel projected lines;
*
************************************************************************
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
*
       REAL       y
       REAL       x1
       REAL       x2
       REAL       s(1:3)
*-----------------------------------------------   define common blocs
       REAL       a, b, c, d
       REAL       czp
       REAL       mia
       COMMON/gzl/a, b, czp, mia, c,d
       SAVE /gzl/
*
*
*-------------------------------------------   declare local variables
       REAL       z
       REAL       dz
       REAL       x
       REAL       dx
       INTEGER    yi
*
*-----------------------------------------------------------------------
       dx=x2-x1
       IF(dx.EQ.0.0) THEN
         CALL gzplt1(NINT(x1),MAXPKY - NINT(y) - 1,
     +   (-a*x1 - b*y - d)/c, s)
         RETURN
       ENDIF
*
       z=(-a*x1 - b*y - d)/c
       dz=-a/c
       yi= MAXPKY - NINT(y) - 1
       x=x1
*
 310   IF(x.LE.x2) THEN
          CALL gzplt1(NINT(x), yi, z, s)
*
          z=z+dz
          x=x+1.0
          GOTO 310
       ENDIF
*
       RETURN
*-------------------------------------------------------   end gzlin6
       END
