C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
       SUBROUTINE gzcol(rgb1, ncol)
***********************************************************************
*
* dummy arguments
*  input:  none
*  output:
*          rgb1 - color look-up table;
*          ncol - number of entries of the color look-up
*
*  function: this subrotine defines the color look-up table and its
*            size;
*
*  called by: show
*
*  calls: none
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
       REAL   rgb1(1:3, 0:127)
       INTEGER    ncol
       INTEGER   i,j
*
       DO 100 i=0, 127
         DO 110 j=1, 3
           rgb1(j,i)=rgb(j,i)
 110     CONTINUE
 100   CONTINUE
*
       ncol=128
       END
