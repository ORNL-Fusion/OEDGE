CC@process opt(3) nosdump nogostmt
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C     UPDATE M. BUSCH 10. 4. 91
C     Unitnummer 79 --> 82 geaendert ---> 85 Groten 12.6.91
C     UPDATE 19. 7.91 GROTEN
C     UPDATE  8. 1.93 Busch SECOND nicht unter UNIX - stattdessen
C                           AIX: mclock
C                           MCLOCK() - verbrauchte Zeit in 100tel Sek
C                           Namen fuer Datei in Kleinbuchstaben
Cjh                         Open fuer Input: gr3zbu.input
C                           Open fuer Input: gr3zbu.inp
C     Update: 28.9.93       AIX: mclock
C                           MCLOCK() - verbrauchte Zeit in 100tel Sek
C     nur EIN Unix Source - daher Zeitmessroutinen in Kommentar
C***********************************************************************
      subroutine gr3zbu(a,c)
      real a(*)
      character(len=*) c
      COMMON /GR3SAV/ WIN,MNP,MNL,MNF,IL0,IF0,IP0,ISY,ISZ,IF1,IF2,IF3,
     >        IF4,IL,L,IW,LW,ZENPRO,MAFL,RT(3,3),GR(3,3),NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR
      SAVE /GR3SAV/
      SAVE
C---- Wie in GR3DIM :
      ifree=15*mnp+1
      CALL gr3zbf(A,A(ISY),A(ISZ),A(IF1),A(IF2),A(IF3),A(IF4),A(IL),
     >            A(L),a(ifree),c)
      end
CC@process opt(3) nosdump nogostmt
      subroutine gr3zbf(x,y,z,lv1,lv2,lv3,lv4,ll,llv,lin,c)
      parameter(mfl=256,maobj=4,in=85,pi=3.141593)

      REAL X(*),Y(*),Z(*),tria(3,3),srcpos(3),prp(3),ex33(3,3),window(4)
!      REAL viepr1(1:4),clr(3,8),windo(4)
      REAL clr(3,8),windo(4)
      INTEGER LV1(*),LV2(*),LV3(*),LV4(*),LL(2,*),LLV(2,*),lin(4,*)
      integer ivisi(3),mersur(mfl),ifcgl(2,maobj),ialg(maobj)
      character(len=2) cin
      character(len=255) cbild
      character(len=255) ccolor
      character(len=*) c
!      logical first,noroom
      logical noroom

      COMMON /gr3bac/ icolba(8),centr,wind(4)
CDEC$ PSECT /gr3bac/ NOSHR
      COMMON /GR4COM/ IFL,JGR(7,MFL),XC(500),YC(500),COMA(3,2)
CDEC$ PSECT /GR4COM/ NOSHR
      COMMON /GR3SAV/ WIN,MNP,MNL,MNF,IL0,IF0,IP0,ISY,ISZ,IF1,IF2,IF3,
     >        IF4,IL,L,IW,LW,ZENPRO,MAFL,RT(3,3),GR(3,3),NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR

      SAVE /gr3bac/, /GR4COM/, /GR3SAV/,FIRST
!      SAVE
!      data first /.true./
      integer :: first=.true.
!      data viepr1 / 0.001, 0.001,0.999,0.999 /
      real :: viepr1(1:4)=(/ 0.001, 0.001,0.999,0.999 /)

      if (first) then
         first=.false.
         write(cin,'(i2.2)') in
CC       open (in, file='gr3zbu.input')
      endif

      lmem=31*mnp-4*if0
      call gzinit(lmem,lin(1,if0+1))
Cjh   open (in, file='gr3zbu.input', status='unknown')
      open (in, file='gr3zbu.inp', status='unknown')
      nsurf=min(maobj,ifl)
*                                   set mirror and face color
      read(in,*)
      read(in,*)
      read(in,*) mirror,icofac
*                                   set background color
      read(in,*)
      read(in,*) ibcgr
*                                   set lightsource position
      read(in,*)
      read(in,*) phi,theta,r
*                                   clipping planes
      read(in,*)
      read(in,*) nearf,fnear,ifarf,far
*                                   border, color
      read(in,*)
      read(in,*) ibofla,ibord
*                                   first 8 colors
      read(in,*)
      read(in,*) ((clr(k,i),k=1,3),i=1,8)
*                                   set rendering algorithms
      read(in,*)
      read(in,*)
      read(in,*) (ialg(i),i=1,nsurf)
*                                        set face gloss
      read(in,*)
      read(in,*) (ifcgl(1,i),ifcgl(2,i),i=1,nsurf)
      if (c.NE.'NOFILE') write(*,*)' transforming ...'
C AIX
c     ITIME=MCLOCK()
C CRAY
c     ITIME=SECOND()
      call gr3ext(x,ier,ex33)
      dx=(ex33(1,3)-ex33(1,1))*.05
      dy=(ex33(2,3)-ex33(2,1))*.05
      if (dx.gt.dy) then
         window(1)=ex33(1,1)-dx
         window(3)=ex33(1,3)+dx
         window(2)=(ex33(2,3)+ex33(2,1)-ex33(1,3)+ex33(1,1))/2-dx
         window(4)=(ex33(2,3)+ex33(2,1)+ex33(1,3)-ex33(1,1))/2+dx
      else
         window(2)=ex33(2,1)-dy
         window(4)=ex33(2,3)+dy
         window(1)=(ex33(1,3)+ex33(1,1)-ex33(2,3)+ex33(2,1))/2-dy
         window(3)=(ex33(1,3)+ex33(1,1)+ex33(2,3)-ex33(2,1))/2+dy
      endif
      dxx = (wind(3)-wind(1))*.05
      dyy = (wind(4)-wind(2))*.05
      windo(1)=wind(1)-dxx
      windo(2)=wind(2)-dyy
      windo(3)=wind(3)+dxx
      windo(4)=wind(4)+dyy
      del=sqrt((window(3)-window(1))**2+(window(4)-window(2))**2)*.0005
      xmit=(ex33(1,3)+ex33(1,1))/2
      ymit=(ex33(2,3)+ex33(2,1))/2
      zmit=(ex33(3,3)+ex33(3,1))/2
*
*                          set edge colors
      do 20 i=1,8
         call gzdefclr(clr(1,i),i-1)
   20 continue
      do 16 j=1,if0
         do 15 i=1,4
            lin(i,j)=0
  15     continue
  16  continue

      ivisi(1)=1
      ivisi(2)=1
      ivisi(3)=1
      isurf=1
*
*                                   set rendering algorithms
      CALL gzaptv(1,isurf,ialg(isurf))
*
*                                       set edge colors
      icol=jgr(3,isurf)
      icolo=icolba(icol)
      CALL gzeclr(isurf,1,icol-1)
      CALL gzeclr(isurf,2,icolo-1)
*
*                                        set face gloss
      CALL gzsfgl(isurf,1,ifcgl(1,isurf))
      CALL gzsfgl(isurf,2,ifcgl(2,isurf))
      iober=jgr(4,isurf)
      ivo=-2 000 000 000
      do 18 i=1,il0
         if (i.gt.iober) then
            mersur(isurf)=ivo
            isurf=isurf+1
            if (isurf.gt.nsurf) then
               write(*,*) 'Zahl vorhandener Flaechen (Objekte):',ifl
               write(*,*) 'Erlaubte Maximalzahl               :',maobj
               write(*,*) 'Die ueberzaehligen werden weggelassen.'
               goto 29
            endif
            iober=jgr(4,isurf)
*
*                                         set rendering algorithm
            CALL gzaptv(1,isurf,ialg(isurf))
*
*                                       set edge colors
            icol=jgr(3,isurf)
            icolo=icolba(icol)
            CALL gzeclr(isurf,1,icol-1)
            CALL gzeclr(isurf,2,icolo-1)
*
*                                        set face gloss
            CALL gzsfgl(isurf,1,ifcgl(1,isurf))
            CALL gzsfgl(isurf,2,ifcgl(2,isurf))
         endif
         linsol=0
         do 17 k=1,2
            iv=llv(k,i)
            if (iv.ne.0) then
               ivo=max(iv,ivo)
               if (abs(lv1(iv)).eq.ll(1,i).and.abs(lv2(iv)).eq.ll(2,i)
     >         .or.abs(lv1(iv)).eq.ll(2,i).and.abs(lv2(iv)).eq.ll(1,i))
     >         then
                   ili=1
               else if
     >         (abs(lv2(iv)).eq.ll(1,i).and.lv3(iv).eq.ll(2,i) .or.
     >          abs(lv2(iv)).eq.ll(2,i).and.lv3(iv).eq.ll(1,i))then
                   ili=2
               else if (lv3(iv).eq.ll(1,i) .and. lv4(iv).eq.ll(2,i) .or.
     >                  lv3(iv).eq.ll(2,i) .and. lv4(iv).eq.ll(1,i))then
                   ili=3
               else if (lv4(iv).eq.ll(1,i).and.abs(lv1(iv)).eq.ll(2,i)
     >         .or. lv4(iv).eq.ll(2,i).and.abs(lv1(iv)).eq.ll(1,i))then
                   ili=4
               else
                   ili=0
               endif
               if (ili.ne.0) lin(ili,iv)=1
            else
               linsol=linsol+1
            endif
  17     continue
         if (linsol.eq.2) then
            ivo=-1
            l1=ll(1,i)
            l2=ll(2,i)
            a=y(l2)-y(l1)
            b=x(l1)-x(l2)
c-----------Naeherung fuer SQRT(A**2+B**2) :       (Hypothenuse)
c-----------maximaler relativer Fehler: 1%
            A1=ABS(a)+1E-37
            B1=ABS(B)+1E-37
            HYPO=A1+B1-1.14265*A1*B1/(A1+B1)
            if (hypo.ne.0) then
               a=a/hypo*del
               b=b/hypo*del
               tria(1,1)=x(l1)-a
               tria(2,1)=y(l1)-b
               tria(3,1)=z(l1)
               tria(1,2)=x(l1)+a
               tria(2,2)=y(l1)+b
               tria(3,2)=z(l1)
               tria(1,3)=x(l2)-a
               tria(2,3)=y(l2)-b
               tria(3,3)=z(l2)
               call gztrg(lmem,lin(1,if0+1),tria,ivisi,isurf,noroom)
               if (noroom) then
                  write(*,*) 'GR3ZBU: Arbeitsspeicher zu klein'
                  write(*,*) 1,'lmem=',lmem,'mnp=',mnp,'if0',if0
                  stop
               endif
               tria(1,1)=x(l1)+a
               tria(2,1)=y(l1)+b
               tria(3,1)=z(l1)
               tria(1,2)=x(l2)+a
               tria(2,2)=y(l2)+b
               tria(3,2)=z(l2)
               tria(1,3)=x(l2)-a
               tria(2,3)=y(l2)-b
               tria(3,3)=z(l2)
               call gztrg(lmem,lin(1,if0+1),tria,ivisi,isurf,noroom)
               if (noroom) then
                  write(*,*) 'GR3ZBU: Arbeitsspeicher zu klein'
                  write(*,*) 2,'lmem=',lmem,'mnp=',mnp,'if0',if0
                  stop
               endif
            endif
         endif
  18  continue
      mersur(isurf)=ivo

  29  isurf=0
 777  isurf=isurf+1
      if (isurf.gt.nsurf) goto 22
      iober=mersur(isurf)
      if (iober.le.0) goto 777
      do 19 i=1,if0
         if (c.NE.'NOFILE' .and. mod(i,500).eq.0)
     >      write(*,*) 'face ',i,' of ',if0
         if (i.gt.iober) then
 555        isurf=isurf+1
            if (isurf.gt.nsurf) goto 22
            iober=mersur(isurf)
            if (iober.le.0) goto 555
         endif
         i1=abs(lv1(i))
         i2=abs(lv2(i))
         i3=lv3(i)
         i4=lv4(i)
         tria(1,1)=x(i1)
         tria(2,1)=y(i1)
         tria(3,1)=z(i1)
         tria(1,2)=x(i2)
         tria(2,2)=y(i2)
         tria(3,2)=z(i2)
         tria(1,3)=x(i3)
         tria(2,3)=y(i3)
         tria(3,3)=z(i3)
         ivisi(1)=lin(1,i)
         ivisi(2)=lin(2,i)
         if (i1.ne.i4 .and. i3.ne.i4) then
            ivisi(3)=0
            call gztrg(lmem,lin(1,if0+1),tria,ivisi,isurf,noroom)
            if (noroom) then
               write(*,*) 'GR3ZBU: Arbeitsspeicher zu klein'
               write(*,*) 3,'lmem=',lmem,'mnp=',mnp,'if0',if0
               stop
            endif
            tria(1,2)=x(i3)
            tria(2,2)=y(i3)
            tria(3,2)=z(i3)
            tria(1,3)=x(i4)
            tria(2,3)=y(i4)
            tria(3,3)=z(i4)
            ivisi(1)=0
            ivisi(2)=lin(3,i)
            ivisi(3)=lin(4,i)
         else
            ivisi(3)=lin(3,i)
         endif
         call gztrg(lmem,lin(1,if0+1),tria,ivisi,isurf,noroom)
         if (noroom) then
            write(*,*) 'GR3ZBU: Arbeitsspeicher zu klein'
            write(*,*) 4,'lmem=',lmem,'mnp=',mnp,'if0',if0
            stop
         endif
   19 continue
   22 continue
*
*                          set face colors
      CALL gzfclr(mirror,icofac)
*
*
*                          activate views
      CALL gzactv(1)
*
*                          set view mapping
      prp(1)=0.
      prp(2)=0.
c     dist=ex33(3,1)
      dist=0
      if (centr.eq.0.) then
         prp(3)=10000.*abs(ex33(3,1)-ex33(3,3))
         CALL gzghvmpg(1,window, viepr1,2, prp, dist, fnear, far)
      else
         prp(3)=centr
         CALL gzghvmpg(1,windo , viepr1,2, prp, dist, fnear, far)
      endif
*
*                          set light source
      srcpos(1)=r*cos(theta/180*pi)*cos(phi/180*pi)+xmit
      srcpos(2)=r*cos(theta/180*pi)*sin(phi/180*pi)+ymit
      srcpos(3)=r*sin(theta/180*pi)                +zmit
      CALL gzsrcv(1,srcpos,0.40,2,.17)
*
*                          set mapping attributes
      CALL gzvier(1,nearf,ifarf,ibcgr-1,ibord-1,ibofla)
*
c     if (c.NE.'NOFILE') then
C AIX
c        TIME =(MCLOCK()-ITIME)/100
C CRAY
c        TIME = second()-ITIME
c        WRITE(*,*)' CPU SECONDS USED FOR TRANSFORMING:',TIME
c        write(*,*)' rendering ...'
c     endif
C AIX
c     ITIME=MCLOCK()
C CRAY
c     ITIME=second()
      if (C.ne.'NOFILE') then
C Busch 8.1.92 KNamen ohne Vorsetzen der Zeichen p_ und c_
C        cbild=      c(:min(6,len(c)))
C        ccolor=      c(:min(6,len(c)))
         cbild=      c
         ccolor=     c
*                          draw
         call gzappv(lmem,lin(1,if0+1),2,cbild,ccolor)
C AIX
c        TIME =(MCLOCK()-ITIME)/100
C CRAY
c        TIME = second()-ITIME
c        WRITE(*,*)' CPU SECONDS USED FOR RENDERING:',TIME
      else
         call gzappv(lmem,lin(1,if0+1),1,cbild,ccolor)
      endif
      close (in)
      end
************************************************************************
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
       SUBROUTINE gzappv(n, mem, picmod, pname, cname)
************************************************************************
*
*  Associate Patch With Workstation
*
*  function:
*    this subroutine draws a patch on the workstation(picmod=ONLINE)
*    or archives the picture in the file name picfile(picmod=OFFLINE);
*
*  CHANGED : M. BUSCH ZDV131
*            10. 4. 91 Unitnummern zum Wegschreiben der Farbmatrix
*                      und der Farbtabelle von 10,11 auf 83,84 geaendert
*                      REAL   *4 --> REAL
*                      INTEGER*4 --> INTEGER
*    Update: 16.8.91 Busch
*    UPDATE: 21.2.92 BUSCH OWNFRB ALS LETZTES ARUGUMENT BEI PHFRB
*                          UND EXTERNAL OWNFRB
*  CHANGED. M. BUSCH a) 6.1.93 Common Bloecke nicht mehr ueber include
*                       sondern als Source eingefuegt
*                    b) 6.1.93 alternatives Abspeichern zur Weiterver-
*                              arbeitung unter VM (graPHIGS) oder
*                              UNIX mit vx,pbmplus
*  CHANGED. G. Groten   7.1.93 Zufallszahlen mit RASEED/RANGEN
*                       AIX
*                       Programme wurde aus /KFA/tools/forf77/source
*                       und umbenannt: GZSEED/GZNGEN
*  CHANGED. M. BUSCH    10.3.93 MINREAL, MAXREAL +- 1.e.38 statt
*                               +-1.e50 (Anpassung an IEEE Format)
*  CHANGED. G. Groten   10.3.98 unused argument edges bei gz1appv entfernt
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
*                          edge number i; edges(3,i), edges(4,i) aregztran1.f
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
       INTEGER   n
       INTEGER   mem(n)
C Busch  25.6.99
        REAL rmem(n)
        LOGICAL lmem(n)
C Busch
       INTEGER    picmod
       CHARACTER(len=*) pname, cname
*
*

C Busch 25.6.99*
       CALL gz1appv(picmod, pname, cname,
     +              rmem(vrtbeg), lmem(prsbeg), mem(fcebeg),
     +              mem(fvbeg), rmem(vnbeg),rmem(fnbeg))
*
*
       END
*
*---------------------------------------------------------------------
*
       SUBROUTINE gz1appv(picmod, pname, cname,
     +                    vertex, edgprs, faces, fcevrt, vnrm, fnrm)
       IMPLICIT NONE
C      EXTERNAL OWNFRB
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
       CHARACTER(len=1) picture(0:MAXPKX-1, 0:MAXPKY-1)
*
       COMMON/gzzbuff/zbuff
       COMMON/feld/picture
       SAVE /gzzbuff/ ,/feld/
*
*
*
       INTEGER    picmod
       CHARACTER(len=*)  pname, cname
       REAL    vertex(1:3,*)
       LOGICAL    edgprs(*)
       INTEGER    faces(1:6, *)
       INTEGER    fcevrt(1:2, *)
       REAL       vnrm(1:3, *)
       REAL       fnrm(1:3, *)
*
*
*
*
*  declare common blocs
       INTEGER    cfcli,cfcno,cacc
       COMMON/gzcurr/cfcli,cfcno,cacc
       SAVE /gzcurr/
*
       REAL    dzz,ro
       COMMON/gzpt/dzz
       SAVE /gzpt/
*
       LOGICAL donenrm(1:MAXOBJ)
       COMMON/gzdone/donenrm
       SAVE /gzdone/
*
       REAL       dcprp(1:3), dcdist
       COMMON/gzdc/dcprp,dcdist
       SAVE /gzdc/
*
*
*  declare local variables
       LOGICAL  ret, w, gor, cnst
       CHARACTER(len=8)  nullchr
       CHARACTER(len=1)  chr
       INTEGER      xmin,xmax, ymin,ymax, x,y
*                          current view index
       INTEGER      view
*                          current object index
       INTEGER      obj
*
*                          current transformation matrix
       REAL         matrix(1:4,1:4)
*
*                   work vertex,vertex normal and face normal tables;
       REAL     fn(1:3)
*
*                          current face index
       INTEGER      fce
       INTEGER      fcel
*
*                current triangle(vertices,normals and edges presence);
       REAL     trg(1:3,1:3),nrm(1:3,1:3),v1(1:3),v2(1:3),fnorm(1:3)
       LOGICAL  edg(1:3)
*
*                triangles after clipping;
       REAL     trgs(1:3,1:12), nrms(1:3,1:12), ctrg(1:3)
       LOGICAL  edgs(1:12)
       INTEGER   ntrgs
*
*                current window, viewport and prp;
       REAL     wdw(1:4), vpt(1:4), prp1(1:3)
*
       REAL     scy,scx, sc
       REAL     mi
*
*                current triangle number;
       INTEGER    trgn
*
*                triangle in DC;
       REAL      dctria(1:3,1:3)
*
*                triangle after projection;
       REAL      ptria(1:3,1:3)
*
*                triangle vertices color;
       REAL      rgb1(1:3,1:3), rgb2(1:3,1:3)
       REAL      rgb3(1:3)
*                current face gloss
       INTEGER   gls
*
*                current source attributes;
       REAL      spos(1:3), spwr, sambnt
*
*                plane parameters
       REAL      par(1:4)
*
*
       INTEGER    sprtype, stype, vst
*
       INTEGER    i,j,k, linind, kk, trgind
*
       REAL       dst
       COMMON/gzw/wdw,vpt,prp1,sc,dst
       SAVE /gzw/
*
       REAL    dcsrcpos(1:3), d1src(1:3)
       COMMON/gzdcsrc/rgb3,dcsrcpos,spwr,sambnt,stype,sprtype,gls
       SAVE /gzdcsrc/
*
*----------------------------------------------------------------------
       CALL GZSEED(0)
       chr=CHAR(0)
       DO 100 i=1,8
         nullchr(i:i)=chr
 100   CONTINUE
*
*  check dummy arguments
       ret=.FALSE.
       IF(picmod.NE.ONLINE .AND. picmod.NE.OFFLINE) THEN
          ret=.TRUE.
          WRITE(DEVIND,*)'  GZAPPV  INVALID PICTURE MODE'
       ENDIF
*
*      IF(picmod.EQ.ONLINE) THEN
*         IF(LLE(picfile,nullchr) .OR. LEN(picfile).NE. 8) THEN
*         ret=.TRUE.
*         WRITE(DEVIND,*)'  GZAPPV  INVALID PICTURE NAME'
*         ENDIF
*      ENDIF
*
*
       IF(ret) THEN
          WRITE(DEVIND,*)'  GZAPPV  NO ACTION IS PERFORMED'
          RETURN
       ENDIF
*
*-----------------------------------------------------------------------
*
*
       DO 110 view=MAXVIEW, 1, -1
*
       IF(.NOT.vactiv(view)) GOTO 110
*
       DO 400 i=1,4
         vpt(i)=vieprt(i,view)
         wdw(i)=window(i,view)
 400   CONTINUE
*
       scy=(vpt(4)-vpt(2))/(wdw(4)-wdw(2))
       scx=(vpt(3)-vpt(1))/(wdw(3)-wdw(1))
       IF(scy.GT.scx) THEN
          sc=scx
       ELSE
          sc=scy
       ENDIF
*
       ro =(vpt(3)-vpt(1))/150.0
       DO 410 i=1,3
          prp1(i)=prp(i,view)
          d1src(i)=srcpos(i,view)
          spos(i)=srcpos(i,view)
 410   CONTINUE
       dst=dist(view)
       sprtype=prtype(view)
       spwr=srcpwr(view)
       sambnt=ambnt(view)
       stype=srctyp(view)
*
       IF(stype.EQ.POINT) THEN
          dcsrcpos(1)= sc*(d1src(1)-wdw(1)) + vpt(1)
          dcsrcpos(2)= sc*(d1src(2)-wdw(2)) + vpt(2)
          dcsrcpos(3)= sc*d1src(3)
       ENDIF
*
       dcprp(1)= sc*(prp1(1) - wdw(1))+ vpt(1)
       dcprp(2)= sc*(prp1(2) - wdw(2))+vpt(2)
       dcprp(3)= sc*prp1(3)
       dcdist = sc*dist(view)
*
*
*
*                          clear the viewport
         xmin= NINT(vpt(1))
         ymax=-NINT(vpt(2))+MAXPKY-1
         xmax= NINT(vpt(3))
         ymin=-NINT(vpt(4))+MAXPKY-1
*
         chr= CHAR(bcgind(view))
         DO 500 y=ymin, ymax
          DO 510 x=xmin, xmax
            zbuff(x,y)=MINREAL
            picture(x,y)=chr
 510      CONTINUE
 500     CONTINUE
*
*
*
       DO 120 obj=1, MAXOBJ
        IF(.NOT.objprs(obj)) GOTO 120
        IF(.NOT.objact(view,obj)) GOTO 120
*
       vst = vstyle(view,obj)
       IF(vst.EQ.GSHADED .OR. vst.EQ.PSHADED .OR. vst.EQ.CONST)
     + THEN
          w=.FALSE.
       ELSE
          w=.TRUE.
       ENDIF
       IF(vst.EQ.GWSHADED .OR. vst.EQ.GSHADED) THEN
          gor = .TRUE.
       ELSE
          gor= .FALSE.
       ENDIF
*
       IF(vst.EQ.WCONST .OR. vst.EQ.CONST) THEN
          cnst= .TRUE.
       ELSE
          cnst=.FALSE.
       ENDIF
*
*                          compute normals
        IF(.NOT.donenrm(obj)) THEN
            CALL gz1nrml(vertex, faces, fcevrt, vnrm, fnrm, obj)
            donenrm(obj)=.TRUE.
        ENDIF
*
*                          compute current transformation matrix;
       CALL gzmult(vtrans(1,1,view), mtrans(1,1,obj), matrix)
*
*
       fcel=fceind(obj) + faceno(obj) - 1
       DO 130  fce=fceind(obj), fcel
*
        CALL gzdotran(fnrm(1,fce), fn, matrix, 0)

        DO 600 i=1,3
          k=faces(i,fce)
          CALL gzdotran(vertex(1,k), trg(1,i), matrix, 1)
          CALL gzdotran(vnrm(1,k), nrm(1,i), matrix, 0)
*
*
          k=faces(i+3,fce)
          edg(i)=edgprs(k)
 600    CONTINUE
*
        CALL gzclpp(trg,nrm,edg,  ntrgs,trgs,nrms,edgs, view)
        IF(ntrgs.EQ.0) GOTO 130
        kk=3*ntrgs
*                          compute plane parameters
        par(1)=fn(1)
        par(2)=fn(2)
        par(3)=fn(3)
        par(4)=-( par(1)*(sc*(trg(1,1)-wdw(1)) + vpt(1)) +
     +         par(2)*( sc*(trg(2,1)-wdw(2)) + vpt(2))     +
     +            par(3)*sc*trg(3,1))
*
       IF(sprtype.EQ.PERSPECTIVE) THEN
         v1(1)= fn(1)
         v1(2)= fn(2)
         v1(3)= fn(3)
         v2(1)= prp1(1) - trgs(1,1)
         v2(2)= prp1(2) - trgs(2,1)
         v2(3)= prp1(3) - trgs(3,1)
*
             IF(v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3) .GT.0.0) THEN
              DO 620 i=1, 3
*
                 rgb1(1,i)=fouclr(1,obj)
                 rgb1(2,i)=fouclr(2,obj)
                 rgb1(3,i)=fouclr(3,obj)
 620          CONTINUE
*
                 fnorm(1)= fn(1)
                 fnorm(2)= fn(2)
                 fnorm(3)= fn(3)
*
                 rgb3(1)=rgb1(1,1)
                 rgb3(2)=rgb1(2,1)
                 rgb3(3)=rgb1(3,1)
                 dzz=ro
*
                 gls=fougls(obj)
                 cfcli=fouind(obj)
                 cfcno=foucno(obj)-1
                 cacc =dclro(obj)
                 linind=eouind(obj)
             ELSE
                DO 630 i=1, 3
*
                 rgb1(1,i)=finclr(1,obj)
                 rgb1(2,i)=finclr(2,obj)
                 rgb1(3,i)=finclr(3,obj)
 630            CONTINUE
*
                 fnorm(1)=-fn(1)
                 fnorm(2)=-fn(2)
                 fnorm(3)=-fn(3)
*
                 rgb3(1)=rgb1(1,1)
                 rgb3(2)=rgb1(2,1)
                 rgb3(3)=rgb1(3,1)
                 dzz=0
                 gls=fingls(obj)
                 cfcli=finind(obj)
                 cfcno=fincno(obj)-1
                 cacc =dclri(obj)
                 linind=einind(obj)
                 DO 640 k=1,3*ntrgs
*                   k=3*(j-1)+1
                    nrms(1,k)= -nrms(1,k)
                    nrms(2,k)= -nrms(2,k)
                    nrms(3,k)= -nrms(3,k)
 640             CONTINUE
*
*
*
*
              ENDIF
*
*
          DO 140 trgn=1, ntrgs
           trgind = 3*(trgn-1)
*
           DO 610 i=1,3
            k=trgind+i
           mi = (dist(view) - trgs(3, trgind+i))/(prp1(3) -
     +            trgs(3, trgind+i))
*
            dctria(1,i)=trgs(1,k) + mi*(prp1(1)-trgs(1,k))
            dctria(2,i)=trgs(2,k) + mi*(prp1(2)-trgs(2,k))
*
*
*
*
*
            ptria(1,i)= sc*(dctria(1,i) - wdw(1))+ vpt(1)
            ptria(2,i)=  (sc*(dctria(2,i)-wdw(2))+ vpt(2))
            ptria(3,i)= sc*trgs(3,k)
 610      CONTINUE
*
*
          k=trgind
*
       IF(vst.NE.WIRE) THEN
         IF(gor) THEN

*                          shading
          CALL gzshade(trgs(1,trgind+1), nrms(1,trgind+1), rgb1(1,1),
     +             rgb2(1,1), gls,prp1,sprtype,spos,stype,spwr,sambnt)
*
*
          CALL gzshade(trgs(1,trgind+2), nrms(1,trgind+2), rgb1(1,2),
     +             rgb2(1,2), gls,prp1,sprtype,spos,stype,spwr,sambnt)
*
          CALL gzshade(trgs(1,trgind+3), nrms(1,trgind+3), rgb1(1,3),
     +             rgb2(1,3), gls,prp1,sprtype,spos,stype,spwr,sambnt)
*
*                          draw a triangle
          CALL gzdrf1(ptria,par,rgb2)
        ELSEIF(cnst) THEN
          ctrg(1)= (trgs(1,trgind+1) + trgs(1,trgind+2) +
     +              trgs(1,trgind+3))/3.0
*
          ctrg(2)= (trgs(2,trgind+1) + trgs(2,trgind+2) +
     +              trgs(2,trgind+3))/3.0
          ctrg(3)= (trgs(3,trgind+1) + trgs(3,trgind+2) +
     +              trgs(3,trgind+3))/3.0
*
          CALL gzshade(ctrg, fnorm, rgb1(1,1), rgb2, gls, prp1,
     +                 sprtype, spos, stype, spwr, sambnt)
          CALL gzdrf5(ptria, par, rgb2)
         ELSE
              CALL gzdrf2(ptria,par,nrms(1,trgind+1))
         ENDIF
       ENDIF
          k = trgind+1
*
          IF(w) THEN
            IF(edgs(k))
     +      CALL gzdda1(ptria(1,1),ptria(2,1),ptria(1,2),ptria(2,2),
     +                    par, linind)
*
            IF(edgs(k+1))
     +     CALL gzdda1(ptria(1,2),ptria(2,2),ptria(1,3),ptria(2,3),
     +                    par, linind)
*
            IF(edgs(k+2))
     +        CALL gzdda1(ptria(1,3),ptria(2,3),ptria(1,1),ptria(2,1),
     +                    par, linind)
          ENDIF
*
 140     CONTINUE
         ELSE
*
             IF(fn(3).GE.0.0) THEN
              DO 622 i=1, 3
*
                 rgb1(1,i)=fouclr(1,obj)
                 rgb1(2,i)=fouclr(2,obj)
                 rgb1(3,i)=fouclr(3,obj)
 622          CONTINUE
*
                 fnorm(1)= fn(1)
                 fnorm(2)= fn(2)
                 fnorm(3)= fn(3)
*
                 rgb3(1)=rgb1(1,1)
                 rgb3(2)=rgb1(2,1)
                 rgb3(3)=rgb1(3,1)
                 dzz=ro
*
                 gls=fougls(obj)
                 cfcli=fouind(obj)
                 cfcno=foucno(obj)-1
                 cacc =dclro(obj)
                 linind=eouind(obj)
             ELSE
                DO 632 i=1, 3
*
                 rgb1(1,i)=finclr(1,obj)
                 rgb1(2,i)=finclr(2,obj)
                 rgb1(3,i)=finclr(3,obj)
 632            CONTINUE
*
                 fnorm(1)=-fn(1)
                 fnorm(2)=-fn(2)
                 fnorm(3)=-fn(3)
*
                 rgb3(1)=rgb1(1,1)
                 rgb3(2)=rgb1(2,1)
                 rgb3(3)=rgb1(3,1)
                 dzz=0
                 gls=fingls(obj)
                 cfcli=finind(obj)
                 cfcno=fincno(obj)-1
                 cacc =dclri(obj)
                 linind=einind(obj)
                 DO 642 k=1,3*ntrgs
*                   k=3*(j-1)+1
                    nrms(1,k)= -nrms(1,k)
                    nrms(2,k)= -nrms(2,k)
                    nrms(3,k)= -nrms(3,k)
 642             CONTINUE
*
*
*
*
              ENDIF
*
*
          DO 142 trgn=1, ntrgs
           trgind = 3*(trgn-1)
*
           DO 612 i=1,3
            k=trgind+i
*
*
            ptria(1,i)= sc*(trgs(1,k) - wdw(1))+ vpt(1)
            ptria(2,i)=  (sc*(trgs(2,k)-wdw(2))+ vpt(2))
            ptria(3,i)= sc*trgs(3,k)
 612      CONTINUE
*
*
          k=trgind
*
       IF(vst.NE.WIRE) THEN
         IF(gor) THEN

*                          shading
          CALL gzshade(trgs(1,trgind+1), nrms(1,trgind+1), rgb1(1,1),
     +             rgb2(1,1), gls,prp1,sprtype,spos,stype,spwr,sambnt)
*
*
          CALL gzshade(trgs(1,trgind+2), nrms(1,trgind+2), rgb1(1,2),
     +             rgb2(1,2), gls,prp1,sprtype,spos,stype,spwr,sambnt)
*
          CALL gzshade(trgs(1,trgind+3), nrms(1,trgind+3), rgb1(1,3),
     +             rgb2(1,3), gls,prp1,sprtype,spos,stype,spwr,sambnt)
*
*                          draw a triangle
          CALL gzdrf3(ptria,par,rgb2)
        ELSEIF(cnst) THEN
          ctrg(1)= (trgs(1,trgind+1) + trgs(1,trgind+2) +
     +              trgs(1,trgind+3))/3.0
*
          ctrg(2)= (trgs(2,trgind+1) + trgs(2,trgind+2) +
     +              trgs(2,trgind+3))/3.0
          ctrg(3)= (trgs(3,trgind+1) + trgs(3,trgind+2) +
     +              trgs(3,trgind+3))/3.0
*
          CALL gzshade(ctrg, fnorm, rgb1(1,1), rgb2, gls, prp1,
     +                 sprtype, spos, stype, spwr, sambnt)
          CALL gzdrf5(ptria, par, rgb2)
         ELSE
              CALL gzdrf4(ptria,par,nrms(1,trgind+1))
         ENDIF
       ENDIF
          k = trgind+1
*
          IF(w) THEN
            IF(edgs(k))
     +      CALL gzdda2(ptria(1,1),ptria(2,1),ptria(1,2),ptria(2,2),
     +                    par, linind)
*
            IF(edgs(k+1))
     +     CALL gzdda2(ptria(1,2),ptria(2,2),ptria(1,3),ptria(2,3),
     +                    par, linind)
*
            IF(edgs(k+2))
     +        CALL gzdda2(ptria(1,3),ptria(2,3),ptria(1,1),ptria(2,1),
     +                    par, linind)
          ENDIF
*
 142     CONTINUE
        ENDIF
 130    CONTINUE
 120    CONTINUE
*                          draw a viewport boundary
       IF(brdact(view).EQ.ON) THEN
         chr=CHAR(brdind(view))
         DO 530 y=ymin, ymax
            picture(xmin,y)=chr
            picture(xmax,y)=chr
 530     CONTINUE
         DO 550 x=xmin, xmax
            picture(x,ymax)=chr
            picture(x,ymin)=chr
 550     CONTINUE
       ENDIF
*
 110    CONTINUE
*
       x=MAXPKX
       y=MAXPKY
       i=1
       j=1
       k=1
       IF(picmod.EQ.ONLINE) THEN
C Busch unter UNIX nicht aufrufen
C       CALL PHFRBMA(X,Y,I,J,K,OWNFRB)
       ELSE
        CALL gzb(pname, cname)
       ENDIF
       END
       SUBROUTINE gzb(pname, cname)
       CHARACTER(len=*)  pname, cname
       CHARACTER(len=255)  p, c

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
*                          look up table;
*
       COMMON/gzclrs/fouind, finind,
     +               eouind, einind, foucno, fincno,
     +               fouclr, finclr,
     +               eouclr, einclr,
     +               fougls, fingls, dclro, dclri, rgb
       SAVE /gzclrs/
*
       character(len=1024) a(1:1024)
C Busch 6.1.93
       COMMON/feld/a
       SAVE /feld/
       integer i
       p=pname
       c=cname
       l=len(pname)

       do 201 i=1,l
          if(pname(i:i).eq.' ') goto 202
 201   CONTINUE
 202   ll=i-1

       call grgif(p(1:ll)//'.gif', MAXCLR, rgb, 1024, 1024, a)
       END
