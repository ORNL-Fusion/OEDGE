      SUBROUTINE GR3HID(X,Y,Z,IFACE1,IFACE2,IFACE3,IFACE4,MNF,INDF,
     $                  ILINE,LNFC,MNL,INDL,EX4,WORK,LENW,IER,MNP)

c  x,y,z           is an array-group of three dimensional points.
c  iface 1 2 3 4   is an array-group of integers which indicate
c                  the points in x,y,z which belong to the indf
c                  faces under consideration.
c  indf            is the current number of faces in iface. if faces are
c                  actually triangles, one point can be entered twice.
c  iline(2,mnl)    is an array which indicates the endpoints of line
c                  segments and mnl is the maximum number of lines
c                  which will be considered.
c  indl            indicates the current number of lines under
c                  consideration. in practice indl should be somewhat
c                  less than mnl since some lines may be added to
c                  improve the picture when faces are non-planar.
c                  the points for iline are in x,y,z.
c  lnfc(2,mnl)     is an array which indicates the indicies of the faces
c                  to which a line belongs.  zero is entered for
c                  nonexistant faces. if lnfc(1,i)=j and lnfc(2,i)=k,
c                  the i-th line will not be compared to the j-th or
c                  k-th face and the points making up line i are assumed
c                  to be among points numbered 1,2, and 4 in face j and
c                  2,3, and 4 in face k. the last is irrelevant unless
c                  the normals of the trianges formed by 1,2,4 and 2,3,4
c                  have z-coordinate of different sign.
c  work(lenw)      working space with length at least 46*mnp



c  hdden3 compares faces to lines and sends visible segments to printer.
c  if a line is compared to a face on which it lies, it may dissappear
c  due to roundoff error.

c     it would be fairly easy to eliminate one half of the
c     faces and the edges in case the surface to be represented
c     is closed.  that would cause a speedup by a factor of 4.


c---- integer iface1(mnf),iface2(mnf),iface3(mnf),iface4(mnf)
c---- integer iline(2,mnl),lnfc(2,mnl)


c     .. scalar arguments ..
      INTEGER           IER,INDF,INDL,LENW,MNF,MNL,MNP
c     ..
c     .. array arguments ..
      REAL              EX4(4),WORK(LENW),X(MNP),Y(MNP),Z(MNP)
      INTEGER           IFACE1(*),IFACE2(*),IFACE3(*),IFACE4(*),
     $                  ILINE(2,*),LNFC(2,*)
c     ..
c     .. local scalars ..
      REAL              XIDONE,XMNF
      INTEGER           IDONE,IW,IW1,IW2,IW3,JW,LW,LW1,LW2,NR
      LOGICAL           FLAG
c     ..
c     .. external subroutines ..
      EXTERNAL          GR3HDD,GR3MMF,GR3SOR
c     ..
c     .. intrinsic functions ..
      INTRINSIC         MAX,SQRT
c     ..
c     .. save statement ..
c     ..
c     .. executable statements ..

      IF (LENW.LT.23*MNF+4*MNL) THEN
         IF (IER.NE.2) WRITE (6,FMT=*)
     $       'GR3HID RC=10: WORK IS NOT BIG ENOUGH'
         IF (IER.NE.2 .AND. IER.NE.1) STOP
         IER    = 10
         GO TO 20
*
      ELSE
         CALL GR3MMF(X,Y,Z,IFACE1,IFACE2,IFACE3,IFACE4,WORK,INDF)

         XMNF   = INDF
         NR     = MAX(SQRT(XMNF*0.5),1.0)
         IW     = 13*MNF + 1
         IW1    = IW + NR**2 + 2
         IW2    = IW1 + 4*MNF
c---- do until( .not. flag ) ------------------------------------------
   10    CONTINUE
         FLAG   = .FALSE.

         CALL GR3SOR(NR,WORK,INDF,WORK(IW2),WORK(IW1),WORK(IW),IDONE,
     $               FLAG,MNF,EX4)

         IF (FLAG) THEN
            XIDONE = IDONE
            XMNF   = MNF
            NR     = .75*NR*SQRT(XIDONE/XMNF)
            IF (NR.EQ.0) THEN
               WRITE (*,FMT=*) 'NR .EQ. 0 IN GR3HID'
               STOP
*
            END IF
*
            GO TO 10
*
         END IF
c---- end do until -----------------------------------------------------
         IW3    = IW2 + MNF
         LW     = 1 + 22*MNF
         LW1    = LW + 2*MNL
         LW2    = LW1 + 2*MNL
         JW     = 1 + 5*MNF
         CALL GR3HDD(NR,X,Y,Z,IFACE1,IFACE2,IFACE3,IFACE4,ILINE,INDF,
     $               INDL,LNFC,MNP,MNL,MNF,EX4,WORK,WORK(JW),WORK(LW),
     $               WORK(LW1),WORK(LW2),WORK(IW),WORK(IW1),WORK(IW2),
     $               WORK(IW3))
         IER    = 0
      END IF

   20 CONTINUE
      END
C@PROCESS nosdump nogostmt OPT(3) IL(DIM) fips(f)
C     DEBUG SUBCHK
C     END DEBUG
      SUBROUTINE GR3MMF(X,Y,Z,IFACE1,IFACE2,IFACE3,IFACE4,FMM,INDF)

c finds maximum and minimum coordinates for faces.
c the entries in fmm(j,i) are the min x for face i if j = 1,
c the max x for the face if j = 2, the min y if j = 3, the max y
c if j = 4, and the max z if j = 5.



c     .. scalar arguments ..
      INTEGER           INDF
c     ..
c     .. array arguments ..
      INTEGER           IF0,IF1,IF2,IF3,IF4,IL,IL0,IP0,ISY,ISZ,IW,LQ,LW,
     $                  MAFL,MNF,MNL,MNP,NOTDIM
      REAL              FMM(5,MNF),X(MNP),Y(MNP),Z(MNP)
      INTEGER           IFACE1(MNF),IFACE2(MNF),IFACE3(MNF),IFACE4(MNF)
c     ..
c     .. local scalars ..
      REAL              QX1,QX2,QX3,QX4,QY1,QY2,QY3,QY4,QZ1,QZ2,QZ3,QZ4,
     $                  TMAX1,TMAX2,TMIN1,TMIN2
      INTEGER           I,IND1,IND2,IND3,IND4
c     ..
c     .. intrinsic functions ..
      INTRINSIC         ABS,MAX,MIN
c     ..
c     .. common blocks ..
      COMMON            /GR3SAV/WIN,MNP,MNL,MNF,IL0,IF0,IP0,ISY,ISZ,IF1,
     $                  IF2,IF3,IF4,IL,LQ,IW,LW,ZENPRO,MAFL,RT,GR,NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR
      REAL              WIN,ZENPRO
      REAL              GR(3,3),RT(3,3)
c     ..
c     .. save statement ..
      SAVE /GR3SAV/
      SAVE
c     ..
c     .. executable statements ..

      DO 10 I = 1,INDF
         IND1   = ABS(IFACE1(I))
         IND2   = ABS(IFACE2(I))
         IND3   = IFACE3(I)
         IND4   = IFACE4(I)

         QZ1    = Z(IND1)
         QZ2    = Z(IND2)
         QZ3    = Z(IND3)
         QZ4    = Z(IND4)

         FMM(5,I) = MAX(MAX(QZ1,QZ2),MAX(QZ3,QZ4))

         QX1    = X(IND1)
         QX2    = X(IND2)
         QX3    = X(IND3)
         QX4    = X(IND4)

         TMAX1  = QX2
         IF (QX1.GE.QX2) TMAX1  = QX1
         TMIN1  = QX1
         IF (QX1.GE.QX2) TMIN1  = QX2

         TMAX2  = QX4
         IF (QX3.GE.QX4) TMAX2  = QX3
         TMIN2  = QX3
         IF (QX3.GE.QX4) TMIN2  = QX4

         FMM(1,I) = MIN(TMIN1,TMIN2)
         FMM(2,I) = MAX(TMAX1,TMAX2)

         QY1    = Y(IND1)
         QY2    = Y(IND2)
         QY3    = Y(IND3)
         QY4    = Y(IND4)

         TMAX1  = QY2
         IF (QY1.GE.QY2) TMAX1  = QY1
         TMIN1  = QY1
         IF (QY1.GE.QY2) TMIN1  = QY2

         TMAX2  = QY4
         IF (QY3.GE.QY4) TMAX2  = QY3
         TMIN2  = QY3
         IF (QY3.GE.QY4) TMIN2  = QY4

         FMM(4,I) = MAX(TMAX1,TMAX2)
         FMM(3,I) = MIN(TMIN1,TMIN2)
   10 CONTINUE

      END
C@PROCESS nosdump nogostmt OPT(3) IL(DIM) fips(f)
C     DEBUG SUBCHK
C     END DEBUG
      SUBROUTINE GR3NRM(X,Y,Z,IFACE1,IFACE2,IFACE3,IFACE4,ENRM,INDF,
     $                  DIFF)

c      finds the normal equations for the two trianges composing
c      each of the faces.

***********************************************



c     .. scalar arguments ..
      INTEGER           INDF
c     ..
c     .. array arguments ..
      INTEGER           IF0,IF1,IF2,IF3,IF4,IL,IL0,IP0,ISY,ISZ,IW,LQ,LW,
     $                  MAFL,MNF,MNL,MNP,NOTDIM
      REAL              ENRM(2,4,INDF),X(MNP),Y(MNP),Z(MNP)
      INTEGER           IFACE1(INDF),IFACE2(INDF),IFACE3(INDF),
     $                  IFACE4(INDF)
      LOGICAL           DIFF(0:INDF)
c     ..
c     .. local scalars ..
      REAL              TEMP1,TEMP2,TEMP3,TEMP4,TEMP5,TEMP6
      INTEGER           I,IND1,IND2,IND3,IND4
c     ..
c     .. intrinsic functions ..
      INTRINSIC         ABS
c     ..
c     .. common blocks ..
      COMMON            /GR3SAV/WIN,MNP,MNL,MNF,IL0,IF0,IP0,ISY,ISZ,IF1,
     $                  IF2,IF3,IF4,IL,LQ,IW,LW,ZENPRO,MAFL,RT,GR,NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR
      REAL              WIN,ZENPRO
      REAL              GR(3,3),RT(3,3)
c     ..
c     .. save statement ..
      SAVE /GR3SAV/
      SAVE
c     ..
c     .. executable statements ..
***********************************************
      DO 10 I = 1,INDF
         IFACE1(I) = ABS(IFACE1(I))
         IFACE2(I) = ABS(IFACE2(I))
         IND1   = IFACE1(I)
         IND2   = IFACE2(I)
         IND3   = IFACE3(I)
         IND4   = IFACE4(I)
         TEMP1  = (Y(IND2)-Y(IND1))* (Z(IND4)-Z(IND2)) -
     $            (Y(IND4)-Y(IND2))* (Z(IND2)-Z(IND1))

         TEMP2  = (X(IND2)-X(IND1))* (Z(IND4)-Z(IND2)) -
     $            (X(IND4)-X(IND2))* (Z(IND2)-Z(IND1))

         TEMP3  = (X(IND2)-X(IND1))* (Y(IND4)-Y(IND2)) -
     $            (X(IND4)-X(IND2))* (Y(IND2)-Y(IND1))
         TEMP4  = (Y(IND3)-Y(IND2))* (Z(IND4)-Z(IND2)) -
     $            (Y(IND4)-Y(IND2))* (Z(IND3)-Z(IND2))

         TEMP5  = (X(IND3)-X(IND2))* (Z(IND4)-Z(IND2)) -
     $            (X(IND4)-X(IND2))* (Z(IND3)-Z(IND2))

         TEMP6  = (X(IND3)-X(IND2))* (Y(IND4)-Y(IND2)) -
     $            (X(IND4)-X(IND2))* (Y(IND3)-Y(IND2))

         ENRM(1,1,I) = TEMP1
         IF (TEMP3.LT.0.0) THEN
            ENRM(1,1,I) = -TEMP1
            IFACE1(I) = -IND1
            IF (TEMP6.EQ.0.) IFACE2(I) = -IND2
         END IF
*
         ENRM(1,2,I) = -TEMP2
         IF (TEMP3.LT.0.0) ENRM(1,2,I) = TEMP2
         ENRM(1,3,I) = TEMP3
         IF (TEMP3.LT.0.0) ENRM(1,3,I) = -TEMP3
         ENRM(2,1,I) = TEMP4
         IF (TEMP6.LT.0.0) THEN
            ENRM(2,1,I) = -TEMP4
            IFACE2(I) = -IND2
            IF (TEMP3.EQ.0.) IFACE1(I) = -IND1
         END IF
*
         ENRM(2,2,I) = -TEMP5
         IF (TEMP6.LT.0.0) ENRM(2,2,I) = TEMP5
         ENRM(2,3,I) = TEMP6
         IF (TEMP6.LT.0.0) ENRM(2,3,I) = -TEMP6

         DIFF(I) = TEMP3*TEMP6 .LT. 0.0

         ENRM(1,4,I) = ENRM(1,1,I)*X(IND2) + ENRM(1,2,I)*Y(IND2) +
     $                 ENRM(1,3,I)*Z(IND2)

         ENRM(2,4,I) = ENRM(2,1,I)*X(IND2) + ENRM(2,2,I)*Y(IND2) +
     $                 ENRM(2,3,I)*Z(IND2)

   10 CONTINUE

      END
C@PROCESS nosdump nogostmt OPT(3) IL(DIM) fips(f)
C     DEBUG SUBCHK
C     END DEBUG
      SUBROUTINE GR3BUC(EX4,X0,Y0,X1,Y1,NUM,LIND,NX,NY)

c     gets the indexes of the rectangles (buckets) met by the line
c     segment between (x0,y0) and (x1,y1).




c     .. scalar arguments ..
      REAL              X0,X1,Y0,Y1
      INTEGER           NUM,NX,NY
c     ..
c     .. array arguments ..
      REAL              EX4(4)
      INTEGER           LIND(43,2)
c     ..
c     .. local scalars ..
      REAL              TEMP,TI,XDIV,XINEW,XIOLD,XNULL,XNY,XT,XX0,XX1,
     $                  YDIV,YINEW,YIOLD,YNULL,YT,YY0,YY1
      INTEGER           ITEMP,IX,IXBOT,IXT,IXTOP,IY,IYBOT,IYTOP
      LOGICAL           LGCL
c     ..
c     .. intrinsic functions ..
      INTRINSIC         MIN
c     ..
c     .. save statement ..
      SAVE
c     ..
c     .. executable statements ..

      XDIV   = NX/ (EX4(3)-EX4(1))
      YDIV   = NY/ (EX4(4)-EX4(2))
      XNY    = NY

      IF (X0.LT.X1) THEN
         XX0    = X0
         XX1    = X1
         YY0    = Y0
         YY1    = Y1
*
      ELSE
         XX0    = X1
         XX1    = X0
         YY0    = Y1
         YY1    = Y0
      END IF

      XNULL  = (XX0-EX4(1))*XDIV
      XT     = (XX1-EX4(1))*XDIV

      YNULL  = (YY0-EX4(2))*YDIV
      YT     = (YY1-EX4(2))*YDIV

      LGCL   = YNULL .LE. YT

      XIOLD  = XNULL
      YIOLD  = YNULL

      ITEMP  = XNULL
      TEMP   = ITEMP

      IXBOT  = ITEMP + 1

      ITEMP  = XT
      TEMP   = ITEMP

      IF (TEMP.LT.XT) THEN
         IXTOP  = ITEMP
*
      ELSE
         IXTOP  = ITEMP - 1
      END IF

      IXTOP  = MIN(IXTOP,NX-1)

      NUM    = 0

      DO 20 IX = IXBOT,IXTOP
         XINEW  = IX
         TI     = (XINEW-XNULL)/ (XT-XNULL)
         YINEW  = (1.-TI)*YNULL + TI*YT

         IF (LGCL) THEN
            IYBOT  = MIN(YIOLD+1.,XNY)
            IYTOP  = MIN(YINEW+1.,XNY)
*
         ELSE
            IYTOP  = MIN(YIOLD+1.,XNY)
            IYBOT  = MIN(YINEW+1.,XNY)
         END IF

         XIOLD  = XINEW
         YIOLD  = YINEW

         DO 10 IY = IYBOT,IYTOP
            NUM    = NUM + 1
            LIND(NUM,1) = IX
            LIND(NUM,2) = IY
   10    CONTINUE

   20 CONTINUE

      IXT    = XT
      TEMP   = IXT

      IF (TEMP.LT.XT) THEN
         IX     = IXT + 1
*
      ELSE
         IX     = IXT
      END IF

      YINEW  = YT

      IF (LGCL) THEN
         IYBOT  = MIN(YIOLD+1.,XNY)
         IYTOP  = MIN(YINEW+1.,XNY)
*
      ELSE
         IYTOP  = MIN(YIOLD+1.,XNY)
         IYBOT  = MIN(YINEW+1.,XNY)
      END IF

      DO 30 IY = IYBOT,IYTOP
         NUM    = NUM + 1
         LIND(NUM,1) = IX
         LIND(NUM,2) = IY
   30 CONTINUE

      END
C@PROCESS nosdump nogostmt OPT(3) IL(DIM) fips(f)
C     DEBUG SUBCHK
C     END DEBUG
c----------------------------------------------------------------------
c ibm-version of gr3cbc (must be replaced on cray: gr3cbc cft77 )
      SUBROUTINE GR3CBC(NUM,NR,LIND,INDEXF,IFPT,IFLIST,IFREE)

c cbc-collect buckets
c     the entries in the multibucket sort corresponding to buckets
c     having indicies stored in lind are collected and returned in
c     iflist.  the last face found has index num.  the last entry in
c     lind corresponding to a pair of indicies of a bucket is also num.

c     .. scalar arguments ..
      INTEGER           NR,NUM
c     ..
c     .. array arguments ..
      INTEGER           IF0,IF1,IF2,IF3,IF4,IL,IL0,IP0,ISY,ISZ,IW,LQ,LW,
     $                  MAFL,MNF,MNL,MNP,NOTDIM
      INTEGER         IFLIST(MNF),IFPT(NR*NR+1),INDEXF(MNF*4),LIND(43,2)
      LOGICAL           IFREE(MNF)
c     ..
c     .. local scalars ..
      INTEGER           II,IND,IX,IY,JJ,NUMF
c     ..
c     .. common blocks ..
      COMMON            /GR3SAV/WIN,MNP,MNL,MNF,IL0,IF0,IP0,ISY,ISZ,IF1,
     $                  IF2,IF3,IF4,IL,LQ,IW,LW,ZENPRO,MAFL,RT,GR,NOTDIM
CDEC$ PSECT /GR3SAV/ NOSHR
      REAL              WIN,ZENPRO
      REAL              GR(3,3),RT(3,3)
c     ..
c     .. save statement ..
      SAVE /GR3SAV/
      SAVE
c     ..
c     .. executable statements ..

      NUMF   = 0

      DO 20 II = 1,NUM
         IX     = LIND(II,1)
         IY     = LIND(II,2)

         IF( NR*(IY-1)+IX .EQ. 0 ) GOTO 20
         DO 10 JJ = IFPT(NR*(IY-1)+IX),IFPT(NR*(IY-1)+IX+1) - 1
            IND    = INDEXF(JJ)

            IF (IFREE(IND)) THEN
               IFREE(IND) = .FALSE.
               NUMF   = NUMF + 1
               IFLIST(NUMF) = IND
            END IF

   10    CONTINUE
   20 CONTINUE

      NUM    = NUMF

      DO 30 II = 1,NUMF
         IFREE(IFLIST(II)) = .TRUE.
   30 CONTINUE

      END
C@PROCESS nosdump nogostmt OPT(3) IL(DIM) fips(f)
C     DEBUG SUBCHK
C     END DEBUG
      SUBROUTINE GR3FXL(RMAX,RMIN,PAVAIL,AVAIL,PSTART,PRM,RM,FLAG)
*********************************************************************
* the upper and lower bounds of visible intervals are stored in rm. *
* rm(1,*) is the left endpoint and rm(2,*) is the right endpoint.   *
* pstart points at the visible interval with the smallest left end  *
* point.  prm(*) is the index of the visible interval to the right  *
* of the interval stored in position rm(.,*) in rm.  the array      *
* avail indicates the indexes of available locations in rm.  upon   *
* initialization each index of rm should be in avail and pavail     *
* should point at the greatest index of avail.  the indexes of rm   *
* in avail need not be ordered.  pavail is the                      *
* index of the cell in avail with greatest index which contains an  *
* index of an available cell in rm.  upon entry                     *
* rmin and rmax are the endpoints of an invisible interval on the   *
* line segment under consideration.  if pstart is zero on exit,     *
* the segment under consideration is invisible. if flag is true,    *
* there is no more free space in rm.                                *
*********************************************************************




c     .. parameters ..
      INTEGER           NDI
      PARAMETER         (NDI=2*65)
c     ..
c     .. scalar arguments ..
      REAL              RMAX,RMIN
      INTEGER           PAVAIL,PSTART
      LOGICAL           FLAG
c     ..
c     .. array arguments ..
      REAL              RM(2,NDI)
      INTEGER           AVAIL(NDI),PRM(NDI)
c     ..
c     .. local scalars ..
      REAL              R1,R2
      INTEGER           INDEX,LAST,NEXT
      LOGICAL           LGCL1,LGCL2,LGCL3
c     ..
c     .. save statement ..
      SAVE
c     ..
c     .. executable statements ..

      LAST   = 0

   10 CONTINUE

      NEXT   = PSTART

   20 CONTINUE

      R1     = RM(1,NEXT)
      R2     = RM(2,NEXT)

      LGCL3  = R1 .GE. RMAX

      IF (.NOT. (LGCL3.OR.R2.LE.RMIN)) THEN

         LGCL1  = R1 .GE. RMIN
         LGCL2  = R2 .LE. RMAX

         IF (LGCL1 .AND. LGCL2) THEN
            PAVAIL = PAVAIL + 1
            AVAIL(PAVAIL) = NEXT
            IF (LAST.NE.0) THEN
               PRM(LAST) = PRM(NEXT)
*
            ELSE IF (PRM(NEXT).EQ.0) THEN
               PSTART = 0
               GO TO 30
*
            ELSE
               PSTART = PRM(NEXT)
               GO TO 10
*
            END IF
*
            NEXT   = LAST
*
         ELSE IF (.NOT. (LGCL1.OR.LGCL2)) THEN
            INDEX  = AVAIL(PAVAIL)
            PRM(INDEX) = PRM(NEXT)
            PRM(NEXT) = INDEX
            LAST   = NEXT
            PAVAIL = PAVAIL - 1

            IF (PAVAIL.EQ.0) THEN
               FLAG   = .TRUE.
               GO TO 30
*
            END IF

            RM(1,INDEX) = RMAX
            RM(2,INDEX) = R2
            RM(2,NEXT) = RMIN
            NEXT   = INDEX
            LGCL3  = .TRUE.
*
         ELSE IF (.NOT.LGCL1) THEN
            RM(2,NEXT) = RMIN
*
         ELSE
            RM(1,NEXT) = RMAX
         END IF
*
      END IF

      IF (.NOT.LGCL3) THEN
         LAST   = NEXT
         NEXT   = PRM(NEXT)
         IF (NEXT.NE.0) GO TO 20
      END IF

   30 CONTINUE
      END
C@PROCESS nosdump nogostmt OPT(3) IL(DIM) fips(f)
C     DEBUG SUBCHK
C     END DEBUG
c----- fuer ibm cms ist dies ein leeres unterprogramm
c----- bei mvs (ibm) und bei der cray, also bei batch, ist es anderswo.
cjh      SUBROUTINE GR3OUT
cjh      END
C@PROCESS nosdump nogostmt OPT(3) IL(DIM) fips(f)
C     DEBUG SUBCHK
C     END DEBUG
c-----------------------------------------------------------------------
c draw an arrow (pfeil) with x,y- direction from point xu,yu
c-----------------------------------------------------------------------
      SUBROUTINE GR3PFL(XU,YU,X,Y,ART)

c     .. scalar arguments ..
      REAL              X,XU,Y,YU
      CHARACTER         ART
c     ..
c     .. local scalars ..
      REAL              AM,XN1,XN2,YN1,YN2
      INTEGER           I,J
c     ..
c     .. local arrays ..
      REAL              XX(2,7),XY(2,7),XZ(2,7),YX(2,7),YY(2,7),YZ(2,7)
c     ..
c     .. external subroutines ..
      EXTERNAL          GRDRW,GRJMP
c     ..
c     .. intrinsic functions ..
      INTRINSIC         ABS
c     ..
c     .. save statement ..
      SAVE
c     ..
c     .. data statements ..

      DATA              XX/.0,1.2,1.275,1.5,1.725,1.275,1.725,1.5,1.8,
     $                  3.,3.,2.85,3.,2.85/
      DATA              YX/.0,0.0,0.225,0.0,0.225,-.225,-.225,0.0,0.0,
     $                  0.,0.,0.15,0.,-.15/
      DATA              XY/.0,1.2,1.275,1.5,1.500,1.500,1.725,1.5,1.8,
     $                  3.,3.,2.85,3.,2.85/
      DATA              YY/.0,0.0,0.225,0.0,-.225,0.000,0.225,0.0,0.0,
     $                  0.,0.,0.15,0.,-.15/
      DATA              XZ/.0,1.2,1.275,1.725,1.725,1.275,1.275,1.725,
     $                  1.8,3.,3.,2.85,3.,2.85/
      DATA              YZ/.0,0.0,0.225,0.225,0.225,-.225,-.225,-.225,
     $                  0.0,0.,0.,0.15,0.,-.15/
c     ..
c     .. executable statements ..
      IF (ART.EQ.'X') THEN
         DO 10 I = 1,7
            XN1    = XX(1,I)*X - YX(1,I)*Y
            YN1    = XX(1,I)*Y + YX(1,I)*X
            XN2    = XX(2,I)*X - YX(2,I)*Y
            YN2    = XX(2,I)*Y + YX(2,I)*X
            CALL GRJMP(XU+XN1,YU+YN1)
            CALL GRDRW(XU+XN2,YU+YN2)
   10    CONTINUE
*
      ELSE IF (ART.EQ.'Y') THEN
         DO 30 I = 1,7
            IF (I.EQ.2 .AND. XN2.LT.0.) THEN
               DO 20 J = 2,4
                  YY(1,J) = -YY(1,J)
   20          CONTINUE
            END IF
*
            XN1    = XY(1,I)*X - YY(1,I)*Y
            YN1    = XY(1,I)*Y + YY(1,I)*X
            XN2    = XY(2,I)*X - YY(2,I)*Y
            YN2    = XY(2,I)*Y + YY(2,I)*X
            CALL GRJMP(XU+XN1,YU+YN1)
            CALL GRDRW(XU+XN2,YU+YN2)
   30    CONTINUE
         YY(1,2) = ABS(YY(1,2))
         YY(1,3) = -ABS(YY(1,3))
         YY(1,4) = ABS(YY(1,4))
*
      ELSE IF (ART.EQ.'Z') THEN
         AM     = XZ(2,2)
         DO 40 I = 1,7
            IF (I.EQ.2 .AND. ABS(XN2).LT.ABS(YN2)) THEN
               XZ(2,2) = XZ(1,2)
               YZ(2,2) = -YZ(2,2)
               YZ(1,3) = -YZ(1,3)
               YZ(2,3) = -YZ(2,3)
               XZ(1,4) = XZ(2,4)
               YZ(1,4) = -YZ(1,4)
            END IF
*
            XN1    = XZ(1,I)*X - YZ(1,I)*Y
            YN1    = XZ(1,I)*Y + YZ(1,I)*X
            XN2    = XZ(2,I)*X - YZ(2,I)*Y
            YN2    = XZ(2,I)*Y + YZ(2,I)*X
            CALL GRJMP(XU+XN1,YU+YN1)
            CALL GRDRW(XU+XN2,YU+YN2)
   40    CONTINUE
         XZ(2,2) = AM
         YZ(2,2) = ABS(YZ(2,2))
         YZ(1,3) = ABS(YZ(1,3))
         YZ(2,3) = -ABS(YZ(2,3))
         XZ(1,4) = XZ(2,1)
         YZ(1,4) = -ABS(YZ(1,4))
*
      ELSE
         WRITE (*,FMT=*) 'GR3PFL: ACHSE NICHT X,Y ODER Z'
      END IF
*
      END
