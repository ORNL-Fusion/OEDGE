      SUBROUTINE WriteVacuumGrid
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      INTEGER fp,i1,i2,cut1
      REAL    zlen,lencyl

      COMMON /WRVACCOM/ convertz
      LOGICAL           convertz

      IF (eirntorseg.NE.0.AND..NOT.convertz) THEN
c...    Convert cylindrical approximation toroidal cut surfaces to 
c       EIRENE toroidal approximation:

        IF (eirzaa.GE.0.0.OR.eirzaa.EQ.-1.0) 
     .    CALL ER('WriteVacuumGrid','Invalid EIRZAA value',*99)

        lencyl = -eirzaa
        DO i1 = 1, MAXASC3D
          asc_zmin3D(i1) = asc_zmin3D(i1) / lencyl * 360.0
          asc_zmax3D(i1) = asc_zmax3D(i1) / lencyl * 360.0           
        ENDDO
        convertz = .TRUE.
      ENDIF


      fp = 98

      OPEN(UNIT=fp,FILE='vacuum_grid.dat',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE')

      WRITE(fp,'(A)') 'DIVIMP ADDITIONAL CELL VACUUM GRID POLYGONS'
      WRITE(fp,*) asc_ncell,MAXASC,asc_3dmode,asccode,ascncut
      DO i1 = 1, asc_ncell
        WRITE(fp,*) asc_cell(i1),asc_cell(i1)+1+eirnpgdat,
     .              ascnvertex(i1),asc_region(i1)
        IF (eirntorseg.NE.0) THEN 
          WRITE(fp,'(2E16.7)') asc_zmin3D(MIN(i1,MAXASC3D)),
     .                         asc_zmax3D(MIN(i1,MAXASC3D))
        ELSE
c...      This is a bizarre method of indexing this array in the transfer file...?  Fix?
          WRITE(fp,'(2E16.7)') asc_zmin3D(MIN(i1,MAXASC3D))*100.0,
     .                         asc_zmax3D(MIN(i1,MAXASC3D))*100.0
        ENDIF
        DO i2 = 1, ascnvertex(i1)
          WRITE(fp,'(2E16.7)') ascvertex(2*i2-1,i1)*100.0,
     .                         ascvertex(2*i2  ,i1)*100.0
        ENDDO
      ENDDO

c...  Sticking this here -- better think about a better place:
      IF (eirnsdtor.GT.1) THEN
        WRITE(fp,'(A)') 'LOCATION OF STANDARD GRID NBLOCK SWITCHING '//
     .                  'SURFACES'

        WRITE(fp,*) eirnsdtor+1

        DO i1 = 1, eirnsdtor
          IF (eirntorseg.NE.0) THEN
            WRITE(fp,*) eirsdtor(i1)
          ELSE
            WRITE(fp,*) eirsdtor(i1)*100.0
          ENDIF
        ENDDO

        IF     (eirzaa.EQ.-1.0) THEN
c...      Cylindrical approximation with the "toroidal circumference" based
c         on the x-point r-coordinate:
          zlen = 2.0 * PI * rxp * eirtorfrac * 100.0
        ELSEIF (eirzaa.LT.0.0) THEN
c          CALL ER('SetupToroidalSurfaces','Length in z-direction not '//
c     .            'compatible with 3D standard grid C',*99)
c         SHOULD DO A BETTER JOB HERE
          zlen = eirtorfrac * 360.0
        ELSE
          zlen = eirzaa * 100.0
        ENDIF

        WRITE(fp,*) zlen

      ENDIF

      CLOSE(fp)

      RETURN
99    STOP
      END

c
c ======================================================================
c ======================================================================
c
c Generate vacuum grid:
c
c ======================================================================
c ======================================================================
c
c subroutine: LocalGridRefinement
c
      SUBROUTINE LocalGridRefinement_Old(xmin,xmax,ymin,ymax,region)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      INTEGER region
      REAL    xmin,xmax,ymin,ymax

      REAL       TOL1       ,TOL2
      PARAMETER (TOL1=1.0D-6,TOL2=1.0D-4)

      INTEGER cell,i1,addcount,hcell,hgrid1,hgrid2,hregion,hnvp,iend,
     .        cell2,i2,i3,i4,hascnvertex,icorner,icount,istatus(4),
     .        iloop,ilist(80),nlist,ncpoly,icut1,icut2,hlink(4)
      LOGICAL easycell,centerinside,firstinside,secondinside,
     .        inserthere,output
      REAL    temp,xcen,ycen,xsidemin,xsidemax,ysidemin,ysidemax,
     .        hrvp(8),hzvp(8),hascvertex(80),maxsidelen,vcoef(8,4),
     .        vx(8,4),vy(8,4),cpoly(160)
      REAL*8  t1,t2

      DATA  vcoef / +0.5, -0.5,  0.0, -0.5,  0.0,  0.0, +0.5,  0.0,
     .               0.0, -0.5, -0.5, -0.5, -0.5,  0.0,  0.0,  0.0,
     .               0.0,  0.0, -0.5,  0.0, -0.5, +0.5,  0.0, +0.5,
     .              +0.5,  0.0,  0.0,  0.0,  0.0, +0.5, +0.5, +0.5 /


      output = .FALSE.

c...  Sort bounds if required:
      IF (xmin.GT.xmax) THEN
        temp = xmin
        xmin = xmax
        xmax = temp
      ENDIF
      IF (ymin.GT.ymax) THEN
        temp = ymin
        ymin = ymax
        ymax = temp
      ENDIF

      addcount = 0
      cell    = 0

      iend = asc_ncell
      
      DO iloop = 1, iend

        cell = cell + 1

c...    Find approximate cell center and bounds of the cell:
        xsidemin =  HI
        xsidemax = -HI
        ysidemin =  HI
        ysidemax = -HI
        DO i1 = 1, ascnvertex(cell)
          xsidemin = MIN(xsidemin,ascvertex(2*i1-1,cell))
          xsidemax = MAX(xsidemax,ascvertex(2*i1-1,cell))
          ysidemin = MIN(ysidemin,ascvertex(2*i1  ,cell))
          ysidemax = MAX(ysidemax,ascvertex(2*i1  ,cell))
        ENDDO
        xcen = 0.5 * (xsidemin + xsidemax)
        ycen = 0.5 * (ysidemin + ysidemax)

        IF (xcen.GT.xmin.AND.xcen.LT.xmax.AND.
     .      ycen.GT.ymin.AND.ycen.LT.ymax.AND.
     .      (region.EQ.asc_region(cell).OR.
     .       region.EQ.-1)) THEN

           IF (output) WRITE(0,*) 'TARGET AQUIRED:',cell

           easycell = .TRUE.
           DO i1 = 1, asc_nvp(cell)
             easycell = easycell.AND.asc_link(i1,cell).NE.-1
           ENDDO


c...       Record cell information:
           hcell       = asc_cell  (  cell)     
           hgrid1      = asc_grid  (1,cell)
           hgrid2      = asc_grid  (2,cell)
           hregion     = asc_region(  cell)     
           hnvp        = asc_nvp   (  cell)
           hascnvertex = ascnvertex(cell) 
           DO i1 = 1, 4
             hlink(i1) = asc_link(i1,cell)
           ENDDO
           DO i1 = 1, 2*asc_nvp(  cell)
             hrvp(i1) = asc_rvp(i1,cell)
             hzvp(i1) = asc_zvp(i1,cell)
           ENDDO
           DO i1 = 1, 2*hascnvertex
             hascvertex(i1) = ascvertex(i1,cell)
           ENDDO
c...       Add another point to close HASCVERTEX:
           hascvertex(2*(hascnvertex+1)-1) = hascvertex(1)
           hascvertex(2*(hascnvertex+1)  ) = hascvertex(2)


           IF (.FALSE..AND.easycell) THEN
c...         Easy cell. Divide by four:

             IF (output) WRITE(0,*) 'YEAH! EASY CELL'

c...         Replace focus cell:
             asc_rvp(1+hnvp,cell) = xcen
             asc_rvp(2     ,cell) = xcen
             asc_rvp(2+hnvp,cell) = xcen
             asc_zvp(2+hnvp,cell) = ycen
             asc_rvp(3     ,cell) = xcen
             asc_zvp(3     ,cell) = ycen
             asc_zvp(3+hnvp,cell) = ycen
             asc_zvp(4     ,cell) = ycen
             DO i1 = 1, asc_nvp(cell)
               asc_link(i1,cell) = 999999 
             ENDDO
             ascnvertex(cell) = hnvp
             DO i1 = 1, hnvp
               ascvertex(2*i1-1,cell) = asc_rvp(i1,cell)
               ascvertex(2*i1  ,cell) = asc_zvp(i1,cell)
             ENDDO




c...         Make space for new cell:
             DO cell2 = asc_ncell, cell+1, -1
               asc_cell  (  cell2+3) = cell2 + 3
               asc_grid  (1,cell2+3) = asc_grid  (1,cell2)
               asc_grid  (2,cell2+3) = asc_grid  (2,cell2)
               asc_region(  cell2+3) = asc_region(  cell2)     
               asc_nvp   (  cell2+3) = asc_nvp   (  cell2)
               ascnvertex(  cell2+3) = ascnvertex(  cell2)
               DO i1 = 1, asc_nvp(cell2)
                 asc_link(i1,cell2+3) = asc_link(i1,cell2)
               ENDDO
               DO i1 = 1, 2*asc_nvp(cell2)
                 asc_rvp(i1,cell2+3) = asc_rvp(i1,cell2)
                 asc_zvp(i1,cell2+3) = asc_zvp(i1,cell2)
               ENDDO          
               DO i1 = 1, 2*ascnvertex(cell2) 
                 ascvertex(i1,cell2+3) = ascvertex(i1,cell2)
               ENDDO
             ENDDO

             asc_ncell = asc_ncell + 3

c...         Generate new cells:
             asc_rvp(1,cell+1) = xcen
             asc_zvp(1,cell+1) = hzvp(1)
             asc_rvp(2,cell+1) = hrvp(2)
             asc_zvp(2,cell+1) = hzvp(2)
             asc_rvp(3,cell+1) = hrvp(3)
             asc_zvp(3,cell+1) = ycen
             asc_rvp(4,cell+1) = xcen
             asc_zvp(4,cell+1) = ycen

             asc_rvp(1,cell+2) = xcen
             asc_zvp(1,cell+2) = ycen
             asc_rvp(2,cell+2) = hrvp(2)
             asc_zvp(2,cell+2) = ycen
             asc_rvp(3,cell+2) = hrvp(3)
             asc_zvp(3,cell+2) = hzvp(3)
             asc_rvp(4,cell+2) = xcen
             asc_zvp(4,cell+2) = hzvp(4)

             asc_rvp(1,cell+3) = hrvp(1)
             asc_zvp(1,cell+3) = ycen
             asc_rvp(2,cell+3) = xcen
             asc_zvp(2,cell+3) = ycen
             asc_rvp(3,cell+3) = xcen
             asc_zvp(3,cell+3) = hzvp(3)
             asc_rvp(4,cell+3) = hrvp(4)
             asc_zvp(4,cell+3) = hzvp(4)

c...         Assign polygon sides:
             DO i1 = 1, 3
               asc_cell  (  cell+i1) = cell + i1
               asc_grid  (1,cell+i1) = asc_grid  (1,cell   )
               asc_grid  (2,cell+i1) = asc_grid  (2,cell   )
               asc_region(  cell+i1) = asc_region(  cell   )     
               asc_nvp   (  cell+i1) = asc_nvp   (  cell   )
               ascnvertex(  cell+i1) = asc_nvp   (  cell+i1)
               DO i2 = 1, asc_nvp(cell+i1)
                 asc_link(i2,cell+i1) = 999999
                 i4 = i2 + 1
                 IF (i2.EQ.asc_nvp(cell+i1)) i4 = 1
                 asc_rvp(i2+4,cell+i1) = asc_rvp(i4,cell+i1)
                 asc_zvp(i2+4,cell+i1) = asc_zvp(i4,cell+i1)
               ENDDO
               DO i2 = 1, ascnvertex(cell+i1)
                 ascvertex(2*i2-1,cell+i1) = asc_rvp(i2,cell+i1)
                 ascvertex(2*i2  ,cell+i1) = asc_zvp(i2,cell+i1)
               ENDDO
             ENDDO

c...         Move index along:
             cell = cell + 3

           ELSE
c...         Outch! Not an easy cell:


c...         Attempt to determine the maximum dimension of this cell.  If 
c            there are not 2 measurements that match, then label the 
c            cell as poorly defined and crash:

             maxsidelen = 0.0
             icorner    = 0
             DO i1 = 1, hnvp
               i2 = i1 + 1
               i4 = i1 - 1
               IF (i1.EQ.hnvp) i2 = 1
               IF (i1.EQ.1   ) i4 = hnvp

               maxsidelen = MAX(maxsidelen,ABS(hrvp(i1)-hrvp(i2)))
               maxsidelen = MAX(maxsidelen,ABS(hzvp(i1)-hzvp(i2)))

c...           Find a good corner to use to establish the center of the cell, but
c              prefer to use a cell side rather than a boundary side, since
c              a boundary side can be poorly defined (at present) if the region has
c              already been refined:
               IF (hrvp(i1).EQ.hrvp(i4+hnvp).AND.
     .             hzvp(i1).EQ.hzvp(i4+hnvp))        icorner = i1
               IF (icorner.EQ.0.AND.hlink(i1).EQ.-1) icorner = i1
             ENDDO

c...         Stop the run if a reference corner could not be found.  This should
c            be a very rare event:
             IF (icorner.EQ.0) 
     .         CALL ER('LocalMeshRefinemnt','Unable to identify '//
     .                 'reference vertex for cell',*99)

c...         Establish center of cell along which the cell will be divided:
             IF (icorner.EQ.1.OR.icorner.EQ.4) THEN
               xcen = hrvp(icorner) - 0.5 * maxsidelen
             ELSE
               xcen = hrvp(icorner) + 0.5 * maxsidelen
             ENDIF
             IF (icorner.EQ.1.OR.icorner.EQ.2) THEN
               ycen = hzvp(icorner) + 0.5 * maxsidelen
             ELSE
               ycen = hzvp(icorner) - 0.5 * maxsidelen
             ENDIF

c...         Define corner verticies of new cells:
             CALL RZero(vx,8*4)
             CALL RZero(vy,8*4)
             DO i1 = 1, 4
               DO i2 = 1, 4
                 i3 = i2 + 1
                 IF (i3.GT.hnvp) i3 = 1
                 vx(i2     ,i1) = xcen + vcoef(2*i2-1,i1) * maxsidelen
                 vy(i2     ,i1) = ycen + vcoef(2*i2  ,i1) * maxsidelen
                 vx(i2+hnvp,i1) = xcen + vcoef(2*i3-1,i1) * maxsidelen
                 vy(i2+hnvp,i1) = ycen + vcoef(2*i3  ,i1) * maxsidelen

                 IF (output) WRITE(0,*) ':::',i1,i2,vx(i2,i1),vy(i2,i1)
               ENDDO
             ENDDO

c...         Check if at least one vertex from the parent polygon is inside
c            the new cell.  If there is not, then mark the child cell to be
c            discarded:
             DO i1 = 1, 4 
               istatus(i1) = 0
               DO i2 = 1, hascnvertex   

                 IF (output.AND.i1.EQ.1)
     .             WRITE(0,*) 'HASCV:',hascvertex(2*i2-1),
     .                                 hascvertex(2*i2)

                 IF (hascvertex(2*i2-1).GE.vx(2,i1)-TOL1.AND.
     .               hascvertex(2*i2-1).LE.vx(1,i1)+TOL1.AND.
     .               hascvertex(2*i2  ).GE.vy(2,i1)-TOL1.AND.
     .               hascvertex(2*i2  ).LE.vy(3,i1)+TOL1) istatus(i1)=1
               ENDDO
               IF (output) WRITE(0,*) 'ISTATUS:',istatus(i1)
             ENDDO

c...         Start the rather complicated process of generating the daughter cells from the 
c            parent cell polygon:

c...         Decide if cell center is inside or outside the parent cell polygon:
             icount = 0
             DO i1 = 1, hascnvertex
               i2 = i1 + 1
               IF (i1.EQ.hascnvertex) i2 = 1

               CALL CalcInter
     .           (DBLE(hascvertex(2*i1-1)),DBLE(hascvertex(2*i1  )),
     .            DBLE(hascvertex(2*i2-1)),DBLE(hascvertex(2*i2  )),
     .            DBLE(xcen),DBLE(ycen),DBLE(xcen+100.0),DBLE(ycen),
     .            t1,t2)

               IF (t1.GT.0.0D0.AND.t1.LT.1.0D0.AND.
     .             t2.GT.0.0D0.AND.t2.LT.1.0D0) icount = icount + 1
             ENDDO

             IF (icount.EQ.0.OR.MOD(icount,2).EQ.0) THEN
               centerinside = .FALSE.
               IF (output) WRITE(0,*) 'outside',icount,centerinside
             ELSE
               centerinside = .TRUE.
c...           If the cell center is inside the parent polygon, then
c              there will be at least some part of each child cell
c              inside the polygon, so turn all the child cells on:
               DO i1 = 1,4
                 istatus(i1) = 1
               ENDDO
               IF (output) WRITE(0,*) 'inside ',icount,centerinside
             ENDIF
 
c...         Loop over daughted cell, avoiding those that are outside the parent
c            polygon:
             icount = 0
             DO i1 = 1, 4
               IF (istatus(i1).EQ.0) CYCLE

c...           Make a list of parent polygon sides that are at least partially
c              inside the child focus cell:
               nlist = 0
               DO i2 = 1, hascnvertex
                 IF ((hascvertex(2* i2   -1).GE.vx(2,i1)-TOL1.AND.
     .                hascvertex(2* i2   -1).LE.vx(1,i1)+TOL1.AND.
     .                hascvertex(2* i2     ).GE.vy(2,i1)-TOL1.AND.
     .                hascvertex(2* i2     ).LE.vy(3,i1)+TOL1).OR.
     .               (hascvertex(2*(i2+1)-1).GE.vx(2,i1)-TOL1.AND.
     .                hascvertex(2*(i2+1)-1).LE.vx(1,i1)+TOL1.AND.
     .                hascvertex(2*(i2+1)  ).GE.vy(2,i1)-TOL1.AND.
     .                hascvertex(2*(i2+1)  ).LE.vy(3,i1)+TOL1)) THEN
                   nlist        = nlist + 1
                   ilist(nlist) = i2
                   IF (output) WRITE(0,*) '---->',i1,i2,nlist
                 ELSE

c...               If is possible for a segment to pass through the child cell
c                  without having end points in the cell, so just for fun,
c                  also look for intersections with the spoke sides of the child
c                  cell:

                   icut1 = i1 + 1
                   IF (icut1.GT.4) icut1 = icut1 - 4
                   icut2 = icut1 + 1
                   IF (icut2.GT.4) icut2 = icut2 - 4
		 
                   CALL CalcInter
     .               (DBLE(hascvertex(2* i2   -1)),
     .                DBLE(hascvertex(2* i2     )),
     .                DBLE(hascvertex(2*(i2+1)-1)),
     .                DBLE(hascvertex(2*(i2+1)  )),
     .                DBLE(vx(icut1,i1)),DBLE(vy(icut1,i1)),
     .                DBLE(vx(icut2,i1)),DBLE(vy(icut2,i1)),
     .                t1,t2)
                 
                   IF (t1.GT.0.0D0     .AND.t1.LT.1.0D0     .AND.
     .                 t2.GT.0.0D0-TOL2.AND.t2.LT.1.0D0+TOL2) THEN
                     nlist        = nlist + 1
                     ilist(nlist) = i2
                     IF (output) WRITE(0,*) '-!-->',i1,i2,nlist
                   ENDIF
		 
                   icut1 = i1 + 2
                   IF (icut1.GT.4) icut1 = icut1 - 4
                   icut2 = icut1 + 1
                   IF (icut2.GT.4) icut2 = icut2 - 4
		 
                   CALL CalcInter
     .               (DBLE(hascvertex(2* i2   -1)),
     .                DBLE(hascvertex(2* i2     )),
     .                DBLE(hascvertex(2*(i2+1)-1)),
     .                DBLE(hascvertex(2*(i2+1)  )),
     .                DBLE(vx(icut1,i1)),DBLE(vy(icut1,i1)),
     .                DBLE(vx(icut2,i1)),DBLE(vy(icut2,i1)),
     .                t1,t2)
                 
                   IF (t1.GT.0.0D0     .AND.t1.LT.1.0D0     .AND.
     .                 t2.GT.0.0D0-TOL2.AND.t2.LT.1.0D0+TOL2) THEN
                     nlist        = nlist + 1
                     ilist(nlist) = i2
                     IF (output) WRITE(0,*) '--!->',i1,i2,nlist
                   ENDIF

                 ENDIF

               ENDDO

c...           Build a list of verticies and crop flagged polygon segments to
c              child cell boundaries:

               ncpoly = 0

               DO i2 = 1, nlist

                 i3 = ilist(i2)

c...             Decide if end points are inside the child cell:
                 firstinside  = .FALSE.
                 secondinside = .FALSE.
                 IF (hascvertex(2* i3   -1).GE.vx(2,i1)-TOL1.AND.
     .               hascvertex(2* i3   -1).LE.vx(1,i1)+TOL1.AND.
     .               hascvertex(2* i3     ).GE.vy(2,i1)-TOL1.AND.
     .               hascvertex(2* i3     ).LE.vy(3,i1)+TOL1) 
     .             firstinside  = .TRUE.
                 IF (hascvertex(2*(i3+1)-1).GE.vx(2,i1)-TOL1.AND.
     .               hascvertex(2*(i3+1)-1).LE.vx(1,i1)+TOL1.AND.
     .               hascvertex(2*(i3+1)  ).GE.vy(2,i1)-TOL1.AND.
     .               hascvertex(2*(i3+1)  ).LE.vy(3,i1)+TOL1)
     .             secondinside = .TRUE.

                 IF (output) WRITE(0,*) 'END POINT:',i3,firstinside,
     .                                                  secondinside

                 IF (firstinside) THEN
c...               Add the first point and then splice the segment (note that
c                  there can be more than one vertex from a given segment):
                   ncpoly = ncpoly + 1
                   cpoly(2*ncpoly-1) = hascvertex(2* i3   -1)
                   cpoly(2*ncpoly  ) = hascvertex(2* i3     )
                   IF (output) WRITE(0,*) 'ADD FIRS: ',i3
                 ENDIF
		 
c...             Try to cut the line with the first center point segment:
		 icut1 = i1 + 1
                 IF (icut1.GT.4) icut1 = icut1 - 4
                 icut2 = icut1 + 1
                 IF (icut2.GT.4) icut2 = icut2 - 4
		 
                 CALL CalcInter
     .             (DBLE(hascvertex(2* i3   -1)),
     .              DBLE(hascvertex(2* i3     )),
     .              DBLE(hascvertex(2*(i3+1)-1)),
     .              DBLE(hascvertex(2*(i3+1)  )),
     .              DBLE(vx(icut1,i1)),DBLE(vy(icut1,i1)),
     .              DBLE(vx(icut2,i1)),DBLE(vy(icut2,i1)),
     .              t1,t2)

	         IF (output) WRITE(0,'(A,I6,2E23.8)') 'T1,2:',i3,t1,t2
                 
                 IF (t1.GT.0.0D0     .AND.t1.LT.1.0D0     .AND.
     .               t2.GT.0.0D0-TOL2.AND.t2.LT.1.0D0+TOL2) THEN
                   ncpoly = ncpoly + 1
                   cpoly(2*ncpoly-1) = SNGL(DBLE(vx(icut1,i1)) + 
     .               t2 * (DBLE(vx(icut2,i1)) - DBLE(vx(icut1,i1))))
                   cpoly(2*ncpoly  ) = SNGL(DBLE(vy(icut1,i1)) + 
     .               t2 * (DBLE(vy(icut2,i1)) - DBLE(vy(icut1,i1))))

                   IF (output) WRITE(0,*) 'CUTTING1: ',i3,t1,t2
                 ENDIF
		 
                 icut1 = i1 + 2
                 IF (icut1.GT.4) icut1 = icut1 - 4
                 icut2 = icut1 + 1
                 IF (icut2.GT.4) icut2 = icut2 - 4
		 
                 CALL CalcInter
     .             (DBLE(hascvertex(2* i3   -1)),
     .              DBLE(hascvertex(2* i3     )),
     .              DBLE(hascvertex(2*(i3+1)-1)),
     .              DBLE(hascvertex(2*(i3+1)  )),
     .              DBLE(vx(icut1,i1)),DBLE(vy(icut1,i1)),
     .              DBLE(vx(icut2,i1)),DBLE(vy(icut2,i1)),
     .              t1,t2)
                 
                 IF (t1.GT.0.0D0     .AND.t1.LT.1.0D0     .AND.
     .               t2.GT.0.0D0-TOL2.AND.t2.LT.1.0D0+TOL2) THEN
                   ncpoly = ncpoly + 1
                   cpoly(2*ncpoly-1) = SNGL(DBLE(vx(icut1,i1)) + 
     .               t2 * (DBLE(vx(icut2,i1)) - DBLE(vx(icut1,i1))))
                   cpoly(2*ncpoly  ) = SNGL(DBLE(vy(icut1,i1)) + 
     .               t2 * (DBLE(vy(icut2,i1)) - DBLE(vy(icut1,i1))))
                   IF (output) WRITE(0,*) 'CUTTING2: ',i3,t1,t2
                 ENDIF
		 
                 IF (secondinside) THEN
c...               Add the second vertex:
                   ncpoly = ncpoly + 1
                   cpoly(2*ncpoly-1) = hascvertex(2*(i3+1)-1)
                   cpoly(2*ncpoly  ) = hascvertex(2*(i3+1)  )
                   IF (output) WRITE(0,*) 'ADD LAST: ',i3
                 ENDIF

               ENDDO

c...           Cull the duplicates:
               DO i2 = 1, ncpoly-1
                 DO i3 = i2+1, ncpoly
                   IF (cpoly(2*i2-1).EQ.cpoly(2*i3-1).AND.
     .                 cpoly(2*i2  ).EQ.cpoly(2*i3  )) THEN
                     DO i4 = i3, ncpoly-1
                       cpoly(2*i4-1) = cpoly(2*(i4+1)-1)
                       cpoly(2*i4  ) = cpoly(2*(i4+1)  )
                     ENDDO
                     ncpoly = ncpoly - 1                     
                   ENDIF
                 ENDDO
               ENDDO

               IF (centerinside) THEN
c...             Have to insert the center point of the parent cell:

                 i2 = i1 + 1
                 IF (i2.GT.4) i2 = i2 - 4

                 inserthere = .FALSE.
                 i3 = 0                     
                 DO WHILE (.NOT.inserthere.AND.i3.LT.ncpoly)
                   i3 = i3 + 1
                   IF (i1.EQ.1.OR.i1.EQ.3) THEN
                     IF (cpoly(2*i3-1).EQ.vx(i2,i1)) inserthere = .TRUE.
                   ELSE
                     IF (cpoly(2*i3  ).EQ.vy(i2,i1)) inserthere = .TRUE.
                   ENDIF 
                   IF (inserthere) THEN
c...                 Make room for new point:
                     DO i4 = ncpoly+1, i3+1, -1
                       cpoly(2*i4-1) = cpoly(2*(i4-1)-1)
                       cpoly(2*i4  ) = cpoly(2*(i4-1)  )                       
                     ENDDO
                     ncpoly = ncpoly + 1
                     cpoly(2*(i3+1)-1) = xcen
                     cpoly(2*(i3+1)  ) = ycen                                 

                     IF (output) WRITE(0,*) 'INSERTING AFTER',i3
                   ENDIF
                 ENDDO

                 IF (.NOT.inserthere) 
     .             CALL ER('LocalMeshRefinement','Unable to position '//
     .                     'parent center point in child cell',*99)
               ENDIF

c...           Close cpoly for convenience:
               cpoly(2*(ncpoly+1)-1) = cpoly(1) 
               cpoly(2*(ncpoly+1)  ) = cpoly(2) 


               IF (output) WRITE(0,*) 'XCEN.YCEN:',xcen,ycen,icount

               DO i2 = 1, ncpoly+1
                 IF (output) WRITE(0,*) 'POLY:',i1,cpoly(2*i2-1),
     .                                             cpoly(2*i2)
               ENDDO


               IF (icount.GT.0) THEN
c...             Make space for new cell:
                  DO cell2 = asc_ncell, cell+icount, -1
                   asc_cell  (  cell2+1) = cell2 + 1
                   asc_grid  (1,cell2+1) = asc_grid  (1,cell2)
                   asc_grid  (2,cell2+1) = asc_grid  (2,cell2)
                   asc_region(  cell2+1) = asc_region(  cell2)     
                   asc_nvp   (  cell2+1) = asc_nvp   (  cell2)
                   ascnvertex(  cell2+1) = ascnvertex(  cell2)
                   DO i2 = 1, asc_nvp(cell2)
                     asc_link(i2,cell2+1) = asc_link(i2,cell2)
                   ENDDO
                   DO i2 = 1, 2*asc_nvp(cell2)
                     asc_rvp(i2,cell2+1) = asc_rvp(i2,cell2)
                     asc_zvp(i2,cell2+1) = asc_zvp(i2,cell2)
                   ENDDO          
                   DO i2 = 1, 2*ascnvertex(cell2) 
                     ascvertex(i2,cell2+1) = ascvertex(i2,cell2)
                   ENDDO
                 ENDDO
                 asc_ncell = asc_ncell + 1
               ENDIF

c...           Okay, so build the cell already:
               asc_cell  (  cell+icount) = cell + icount
               asc_grid  (1,cell+icount) = hgrid1      
               asc_grid  (2,cell+icount) = hgrid2      
               asc_region(  cell+icount) = hregion          
               asc_nvp   (  cell+icount) = hnvp        
               DO i2 = 1, hnvp
                 asc_link(i2,cell+icount) = 999999
               ENDDO

               ascnvertex(  cell+icount) = ncpoly
               DO i2 = 1, ncpoly
                 ascvertex(2*i2-1,cell+icount) = cpoly(2*i2-1)
                 ascvertex(2*i2  ,cell+icount) = cpoly(2*i2  )
               ENDDO

c...           Assign default standard sides, to be overwritten in the 
c              next block of code if side data is found from the polygon
c              data:
               DO i2 = 1, hnvp
                 asc_rvp(i2     ,cell+icount) = vx(i2     ,i1)
                 asc_zvp(i2     ,cell+icount) = vy(i2     ,i1)
                 asc_rvp(i2+hnvp,cell+icount) = vx(i2+hnvp,i1)
                 asc_zvp(i2+hnvp,cell+icount) = vy(i2+hnvp,i1)
               ENDDO


c...           Last big task is to build the standard sides from the
c              polygon:
               asc_link(1,cell+icount) = -1
               asc_link(2,cell+icount) = -1
               asc_link(3,cell+icount) = -1
               asc_link(4,cell+icount) = -1
               DO i2 = 1, ncpoly
                 IF     (ABS(cpoly(2* i2     )-vy(1,i1)).LT.TOL1.AND.
     .                   ABS(cpoly(2*(i2+1)  )-vy(1,i1)).LT.TOL1) THEN
                   asc_rvp(1     ,cell+icount) = cpoly(2* i2   -1)
                   asc_zvp(1     ,cell+icount) = cpoly(2* i2     )
                   asc_rvp(1+hnvp,cell+icount) = cpoly(2*(i2+1)-1)
                   asc_zvp(1+hnvp,cell+icount) = cpoly(2*(i2+1)  )
                   IF (ABS(asc_zvp(1,cell+icount)-ymin).GT.TOL1)
     .               asc_link(1,cell+icount) = 999999
                   IF (output) WRITE(0,*) 'FOUND 1'
                 ELSEIF (ABS(cpoly(2* i2   -1)-vx(2,i1)).LT.TOL1.AND.
     .                   ABS(cpoly(2*(i2+1)-1)-vx(2,i1)).LT.TOL1) THEN
                   asc_rvp(2     ,cell+icount) = cpoly(2* i2   -1)
                   asc_zvp(2     ,cell+icount) = cpoly(2* i2     )
                   asc_rvp(2+hnvp,cell+icount) = cpoly(2*(i2+1)-1)
                   asc_zvp(2+hnvp,cell+icount) = cpoly(2*(i2+1)  )
                   IF (ABS(asc_rvp(2,cell+icount)-xmin).GT.TOL1) 
     .               asc_link(2,cell+icount) = 999999
                   IF (output) WRITE(0,*) 'FOUND 2'
                 ELSEIF (ABS(cpoly(2* i2     )-vy(3,i1)).LT.TOL1.AND.
     .                   ABS(cpoly(2*(i2+1)  )-vy(3,i1)).LT.TOL1) THEN
                   asc_rvp(3     ,cell+icount) = cpoly(2* i2   -1)
                   asc_zvp(3     ,cell+icount) = cpoly(2* i2     )
                   asc_rvp(3+hnvp,cell+icount) = cpoly(2*(i2+1)-1)
                   asc_zvp(3+hnvp,cell+icount) = cpoly(2*(i2+1)  )
                   IF (ABS(asc_zvp(3,cell+icount)-ymax).GT.TOL1) 
     .               asc_link(3,cell+icount) = 999999
                   IF (output) WRITE(0,*) 'FOUND 3'
                 ELSEIF (ABS(cpoly(2* i2   -1)-vx(4,i1)).LT.TOL1.AND.
     .                   ABS(cpoly(2*(i2+1)-1)-vx(4,i1)).LT.TOL1) THEN
                   asc_rvp(4     ,cell+icount) = cpoly(2* i2   -1)
                   asc_zvp(4     ,cell+icount) = cpoly(2* i2     )
                   asc_rvp(4+hnvp,cell+icount) = cpoly(2*(i2+1)-1)
                   asc_zvp(4+hnvp,cell+icount) = cpoly(2*(i2+1)  )
                   IF (ABS(asc_rvp(4,cell+icount)-xmax).GT.TOL1) 
     .               asc_link(4,cell+icount) = 999999
                   IF (output) WRITE(0,*) 'FOUND 4'
                 ENDIF
               ENDDO

c...           Increase the number of new cells added:
               icount = icount + 1
           
             ENDDO 

             IF (icount.EQ.0) THEN
               CALL ER('LocalGridRefinement','No daughter cells '//
     .                 'added',*99)
             ELSE
c...           Move index along:
               cell = cell + icount - 1
             ENDIF

           ENDIF

         ENDIF
        
      ENDDO


c      STOP 'sljkdfdsf'



      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: Build3DVacuumGrid
c
      SUBROUTINE Build3DVacuumGrid(xcutmin,xcutmax,ycutmin,ycutmax,zcut,
     .                             mode)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      INTEGER cell,cell2,v1,v2,nvp,s1,ncell,fp3,i1,i2,cut1,ncut,mode,
     .        cut2
      REAL    zval,zcut,xmin,xmax,ymin,ymax,
     .        xcutmin,xcutmax,ycutmin,ycutmax

      REAL TOL
      PARAMETER (TOL=1.0E-6)

      WRITE(SLOUT,'(A,1P,4E15.7,0P)') 
     .  'BOUNDS:',xcutmin,xcutmax,ycutmin,ycutmax


      IF     (mode.EQ.2) THEN


        IF (asc_3dmode.EQ.0) THEN

          asc_3Dmode = 2

          IF     (eirzaa.EQ.-1.0) THEN
c...        Cylindrical approximation with the "toroidal circumference" based
c           on the x-point r-coordinate:
            asc_zmin3D(1) = 0.0
            asc_zmax3D(1) = 2.0 * PI * rxp * eirtorfrac
          ELSEIF (eirzaa.LT.0.0) THEN
c            CALL ER('Build3DVacuumGrid','Length in z-direction not '//
c     .              'compatible with 3D vacuum grid B',*99)
            asc_zmin3D(1) = 0.0
            asc_zmax3D(1) = -eirzaa * eirtorfrac
          ELSE
            asc_zmin3D(1) = 0.0
            asc_zmax3D(1) = eirzaa
          ENDIF
 
          ascncut = 1

        ENDIF

c...    Set vacuum grid quantities:
        ncut = ascncut
        DO cut1 = 1, ncut
          IF (.NOT.(asc_zmin3D(cut1).LT.zcut.AND.
     .              asc_zmax3D(cut1).GT.zcut)) CYCLE 

c...      Make room for new cut:
          ascncut = ascncut + 1
          DO cut2 = ascncut, cut1+2, -1
            asc_zmin3D(cut2) = asc_zmin3D(cut2-1)
            asc_zmax3D(cut2) = asc_zmax3D(cut2-1)             
            DO i1 = 1, 4
              asc_xvp3D(1,i1,cut2-1) = asc_xvp3D(1,i1,cut2-2)
              asc_yvp3D(1,i1,cut2-1) = asc_yvp3D(1,i1,cut2-2)
              asc_zvp3D(1,i1,cut2-1) = asc_zvp3D(1,i1,cut2-2)
            ENDDO
          ENDDO

c...      Assign the toroidal boundaries of this segment:
          asc_zmin3D(cut1+1) = zcut
          asc_zmax3D(cut1+1) = asc_zmax3D(cut1)

c...      Define the surface the cuts the grid at ZCUT:
          asc_xvp3D(1,1,cut1) = xcutmax
          asc_yvp3D(1,1,cut1) = ycutmin
          asc_zvp3D(1,1,cut1) = zcut
          asc_xvp3D(1,2,cut1) = xcutmin
          asc_yvp3D(1,2,cut1) = ycutmin
          asc_zvp3D(1,2,cut1) = zcut
          asc_xvp3D(1,3,cut1) = xcutmin
          asc_yvp3D(1,3,cut1) = ycutmax
          asc_zvp3D(1,3,cut1) = zcut
          asc_xvp3D(1,4,cut1) = xcutmax
          asc_yvp3D(1,4,cut1) = ycutmax
          asc_zvp3D(1,4,cut1) = zcut

c...      Adjust toroidal location of the segment that has just been 
c         divided:
          asc_zmax3D(cut1) = zcut
        ENDDO

      ELSEIF (mode.EQ.1) THEN

        IF (asc_3Dmode.EQ.0) THEN

          asc_3Dmode = 1

c...      Set default 3D quantities for all existing 2D vacuum cells:

          DO cell = 1, asc_ncell
            nvp = asc_nvp(cell)

c...        Add 2 sides to the cell:
            asc_nvp3D(cell) = nvp + 2

            DO s1 = 1, nvp
              asc_link3D(s1,cell) = asc_link(s1,cell)
            ENDDO          
            asc_link3D(nvp+1,cell) = -1
            asc_link3D(nvp+2,cell) = -1

            IF (eirzaa.EQ.-1.0) THEN
c...          Cylindrical approximation with the "toroidal circumference" based
c             on the x-point r-coordinate:
              asc_zmin3D(cell) = 0.0
              asc_zmax3D(cell) = 2.0 * PI * rxp * eirtorfrac
            ELSEIF (eirzaa.LT.0.0) THEN
c              CALL ER('Build3DVacuumGrid','Length in z-direction not '//
c     .                'compatible with 3D vacuum grid A',*99)
              asc_zmin3D(cell) = 0.0
              asc_zmax3D(cell) = -eirzaa * eirtorfrac
            ELSE
              asc_zmin3D(cell) = 0.0
              asc_zmax3D(cell) = eirzaa
            ENDIF

c...        Surfaces 1 to nvp:
            DO v1 = 1, nvp
              IF (asc_link(v1,cell).LT.0) CYCLE

              asc_xvp3D(v1,1    ,cell) = asc_rvp(v1+nvp,cell) 
              asc_yvp3D(v1,1    ,cell) = asc_zvp(v1+nvp,cell) 
              asc_zvp3D(v1,1    ,cell) = 0.0
              asc_xvp3D(v1,1+nvp,cell) = asc_rvp(v1+nvp,cell) 
              asc_yvp3D(v1,1+nvp,cell) = asc_zvp(v1+nvp,cell) 
              asc_zvp3D(v1,1+nvp,cell) = asc_zmax3D(cell)
                                                  
              asc_xvp3D(v1,2    ,cell) = asc_rvp(v1+nvp,cell) 
              asc_yvp3D(v1,2    ,cell) = asc_zvp(v1+nvp,cell) 
              asc_zvp3D(v1,2    ,cell) = asc_zmax3D(cell)
              asc_xvp3D(v1,2+nvp,cell) = asc_rvp(v1    ,cell) 
              asc_yvp3D(v1,2+nvp,cell) = asc_zvp(v1    ,cell) 
              asc_zvp3D(v1,2+nvp,cell) = asc_zmax3D(cell)
                                                  
              asc_xvp3D(v1,3    ,cell) = asc_rvp(v1    ,cell) 
              asc_yvp3D(v1,3    ,cell) = asc_zvp(v1    ,cell) 
              asc_zvp3D(v1,3    ,cell) = asc_zmax3D(cell)
              asc_xvp3D(v1,3+nvp,cell) = asc_rvp(v1    ,cell) 
              asc_yvp3D(v1,3+nvp,cell) = asc_zvp(v1    ,cell) 
              asc_zvp3D(v1,3+nvp,cell) = 0.0
                                                  
              asc_xvp3D(v1,4    ,cell) = asc_rvp(v1    ,cell) 
              asc_yvp3D(v1,4    ,cell) = asc_zvp(v1    ,cell) 
              asc_zvp3D(v1,4    ,cell) = 0.0
              asc_xvp3D(v1,4+nvp,cell) = asc_rvp(v1+nvp,cell) 
              asc_yvp3D(v1,4+nvp,cell) = asc_zvp(v1+nvp,cell) 
              asc_zvp3D(v1,4+nvp,cell) = 0.0
            ENDDO

c...        Surfaces with normals in the -z and  +z direction::
            DO s1 = nvp+1, nvp+2
              IF (s1.EQ.nvp+1) zval = 0.0
              IF (s1.EQ.nvp+2) zval = asc_zmax3D(cell)
c...          Find boundary of rectangle that can bound the cell toroidally:
              xmin =  HI
              ymin =  HI
              xmax = -HI
              ymax = -HI
              DO v1 = 1, ascnvertex(cell)
                xmin = MIN(ascvertex(2*v1-1,cell),xmin)
                xmax = MAX(ascvertex(2*v1-1,cell),xmax)
                ymin = MIN(ascvertex(2*v1  ,cell),ymin)
                ymax = MAX(ascvertex(2*v1  ,cell),ymax)
              ENDDO
              asc_xvp3D(s1,1,cell) = xmax
              asc_yvp3D(s1,1,cell) = ymin
              asc_zvp3D(s1,1,cell) = zval 	  
              asc_xvp3D(s1,2,cell) = xmin
              asc_yvp3D(s1,2,cell) = ymin
              asc_zvp3D(s1,2,cell) = zval 	  
              asc_xvp3D(s1,3,cell) = xmin
              asc_yvp3D(s1,3,cell) = ymax
              asc_zvp3D(s1,3,cell) = zval 	  
              asc_xvp3D(s1,4,cell) = xmax
              asc_yvp3D(s1,4,cell) = ymax
              asc_zvp3D(s1,4,cell) = zval 	  
              DO v1 = nvp+1, nvp+4
                IF (v2.GT.nvp) THEN
                  v2 = 1
                ELSE
                  v2 = v1 - nvp + 1
                ENDIF
                asc_xvp3D(s1,v1,cell) = asc_xvp3D(s1,v2,cell)
                asc_yvp3D(s1,v1,cell) = asc_yvp3D(s1,v2,cell)
                asc_zvp3D(s1,v1,cell) = asc_zvp3D(s1,v2,cell)	  
              ENDDO
            ENDDO						       	  
  
          ENDDO
      
        ENDIF

c...    Loop over all cells and subdivide any that are within the 
c       prescribed bounds:
        ncell = asc_ncell

        DO cell = 1, ncell
 
          nvp = asc_nvp(cell)
 
c...      Do not split cell if its span in the z-direction, as determined
c         by side 1 on surface 1 of the cell, does not include ZCUT:

          IF (.NOT.(asc_zmin3D(cell).LT.zcut.AND.
     .              asc_zmax3D(cell).GT.zcut)) CYCLE

          IF (xcutmin.NE.99.0) THEN
c...        Find boundary of rectangle that can bound the cell toroidally:
            xmin =  HI
            ymin =  HI
            xmax = -HI
            ymax = -HI
            DO v1 = 1, ascnvertex(cell)
              xmin = MIN(ascvertex(2*v1-1,cell),xmin)
              xmax = MAX(ascvertex(2*v1-1,cell),xmax)
              ymin = MIN(ascvertex(2*v1  ,cell),ymin)
              ymax = MAX(ascvertex(2*v1  ,cell),ymax)
            ENDDO
          ENDIF

          IF (xcutmin.NE.99.0.AND.
     .        (xmin.LT.xcutmin-TOL.OR.xmax.GT.xcutmax+TOL.OR.
     .         ymin.LT.ycutmin-TOL.OR.ymax.GT.ycutmax+TOL)) CYCLE

c...      Create a new vacuum cell:
          cell2 = asc_ncell + 1

          WRITE(SLOUT,*) 'ADDING CELL:',cell,cell2

c...      Assign 2D quantities:
          asc_cell  (cell2) = cell2
          asc_region(cell2) = asc_region(cell)
          asc_nvp   (cell2) = asc_nvp   (cell)
          DO v1 = 1, 2*asc_nvp(cell)
            asc_rvp(v1,cell2) = asc_rvp(v1,cell)
            asc_zvp(v1,cell2) = asc_zvp(v1,cell)
          ENDDO

c...      Assign 3D quantites:
          asc_nvp3D(cell2) = nvp + 2

c...      Extend connection map:
          DO s1 = 1, nvp
            IF (asc_link(s1,cell).LT.0) THEN
              asc_link  (s1,cell2) = -1
              asc_link3D(s1,cell2) = -1
            ELSE
              asc_link  (s1,cell2) = cell2 + (asc_link(s1,cell)-cell)
              asc_link3D(s1,cell2) = cell2 + (asc_link(s1,cell)-cell)
            ENDIF
          ENDDO
          asc_link3D(nvp+1,cell2) = cell
          asc_link3D(nvp+2,cell2) = asc_link3D(nvp+2,cell)
          asc_link3D(nvp+2,cell) = cell2


c...      Assign verticies to CELL2 and adjust Z coordinates:
          DO s1 = 1, nvp+2
            IF (asc_link3D(s1,cell2).LT.0.AND.s1.LE.nvp) CYCLE
            DO v1 = 1, 2*nvp
              asc_xvp3D(s1,v1,cell2) = asc_xvp3D(s1,v1,cell)
              asc_yvp3D(s1,v1,cell2) = asc_yvp3D(s1,v1,cell)
              asc_zvp3D(s1,v1,cell2) = asc_zvp3D(s1,v1,cell)
              IF (asc_zvp3D(s1,v1,cell2).LT.zcut) 
     .          asc_zvp3D(s1,v1,cell2) = zcut
            ENDDO
          ENDDO
          DO s1 = 1, nvp+2
            IF (asc_link3D(s1,cell).LT.0) CYCLE
            DO v1 = 1, 2*nvp
              IF (asc_zvp3D(s1,v1,cell).GT.zcut) 
     .          asc_zvp3D(s1,v1,cell) = zcut
            ENDDO
          ENDDO

c...      Set z-direction length of cell:
          asc_zmin3D(cell2) = zcut
          asc_zmax3D(cell2) = asc_zmax3D(cell)
          asc_zmax3D(cell) = zcut

c...      Recalculate cell volumes:
          CALL CalcPolygonVolume(asc_rvp (1,cell),asc_zvp(1,cell),4,
     .                           asc_link(1,cell),asc_vol(cell),cell)
          CALL CalcPolygonVolume(asc_rvp (1,cell2),asc_zvp(1,cell2),4,
     .                           asc_link(1,cell2),asc_vol(cell2),cell2)

          asc_ncell = cell2
 
        ENDDO

        fp3=98
        OPEN(UNIT=fp3,FILE='dump.dat',ACCESS='SEQUENTIAL',
     .       STATUS='REPLACE',ERR=96)
        DO cell = 1, asc_ncell
          DO s1 = 1, 6
            IF (asc_link3D(s1,cell).LT.0) CYCLE
            DO v1 = 1, 4        
              WRITE(fp3,'(I6,2(3E14.6,2X),2I6)') 
     .          cell,
     .          asc_xvp3D(s1,v1  ,cell),asc_yvp3D(s1,v1  ,cell),
     .          asc_zvp3D(s1,v1  ,cell),
     .          asc_xvp3D(s1,v1+4,cell),asc_yvp3D(s1,v1+4,cell),
     .          asc_zvp3D(s1,v1+4,cell),s1,v1
            ENDDO
          ENDDO
        ENDDO
        CLOSE(fp3)

      ENDIF

      RETURN
96    WRITE(0,*) 'TROUBLE MAKING BUMP FILE'
99    STOP
      END
c
c ======================================================================
c
c subroutine: AssignAddCellLink
c
      SUBROUTINE AssignAddCellLink
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slcom'

      REAL       TOL
      PARAMETER (TOL=1.0E-06)

      INTEGER
     .  link_old(6,MAXASC)

      INTEGER cell,cell2,v1,v2,v3,v4,n1,s1,s2,i1,i2
 

      IF (asc_3Dmode.EQ.1) THEN


c...    Build new connection map for vacuum grid:

        DO cell = 1, asc_ncell
          DO s1 = 1, asc_nvp3D(cell)
            link_old(s1,cell) = asc_link3D(s1,cell)
          ENDDO
        ENDDO

c...    Reinitialize cell link array:
        DO cell = 1, asc_ncell
          DO s1 = 1, asc_nvp3D(cell)
            IF (asc_link3D(s1,cell).GT.0) asc_link3D(s1,cell) = 999999
          ENDDO
        ENDDO

c...  For all additional cells, search every other cell and find 
c     coindicent sides.  Link these via asc_link:

      DO cell = 1, asc_ncell

c        WRITE(0,*) 'LINKING: ',cell,asc_ncell

        DO s1 = 1, asc_nvp3D(cell)

          IF (asc_link3D(s1,cell).LT.0.OR.
     .        asc_link3D(s1,cell).NE.999999) CYCLE

          DO cell2 = 1, asc_ncell

	    IF (cell.EQ.cell2) CYCLE

            IF (s1.LE.asc_nvp(cell)) THEN

c...          Side surfaces:

              DO s2 = 1, asc_nvp(cell2)
	  
                IF (asc_link3D(s2,cell2).LT.0.OR.
     .              asc_link3D(s2,cell2).NE.999999) CYCLE

                IF (asc_xvp3D(s1,1,cell).EQ.asc_xvp3D(s2,4,cell2).AND.
     .              asc_yvp3D(s1,1,cell).EQ.asc_yvp3D(s2,4,cell2).AND.
     .              asc_zvp3D(s1,1,cell).EQ.asc_zvp3D(s2,4,cell2).AND.
     .              asc_xvp3D(s1,2,cell).EQ.asc_xvp3D(s2,3,cell2).AND.
     .              asc_yvp3D(s1,2,cell).EQ.asc_yvp3D(s2,3,cell2).AND.
     .              asc_zvp3D(s1,2,cell).EQ.asc_zvp3D(s2,3,cell2).AND.
     .              asc_xvp3D(s1,3,cell).EQ.asc_xvp3D(s2,2,cell2).AND.
     .              asc_yvp3D(s1,3,cell).EQ.asc_yvp3D(s2,2,cell2).AND.
     .              asc_zvp3D(s1,3,cell).EQ.asc_zvp3D(s2,2,cell2).AND.
     .              asc_xvp3D(s1,4,cell).EQ.asc_xvp3D(s2,1,cell2).AND.
     .              asc_yvp3D(s1,4,cell).EQ.asc_yvp3D(s2,1,cell2).AND.
     .              asc_zvp3D(s1,4,cell).EQ.asc_zvp3D(s2,1,cell2)) THEN
 	  
c                   WRITE(0,'(A,4I6)') 'MATCH    : ',cell,s1,cell2,s2
                   asc_link3D(s1,cell) = asc_cell(cell2)
                   asc_link3D(s2,cell2) = asc_cell(cell)
  	  
                ENDIF

              ENDDO

            ELSE

c...          End surfaces:

              IF (s1.EQ.5) s2 = 6
              IF (s1.EQ.6) s2 = 5

              IF (asc_xvp3D(s1,1,cell).EQ.asc_xvp3D(s2,1,cell2).AND.
     .            asc_yvp3D(s1,1,cell).EQ.asc_yvp3D(s2,1,cell2).AND.
     .            asc_zvp3D(s1,1,cell).EQ.asc_zvp3D(s2,1,cell2).AND.
     .            asc_xvp3D(s1,2,cell).EQ.asc_xvp3D(s2,2,cell2).AND.
     .            asc_yvp3D(s1,2,cell).EQ.asc_yvp3D(s2,2,cell2).AND.
     .            asc_zvp3D(s1,2,cell).EQ.asc_zvp3D(s2,2,cell2).AND.
     .            asc_xvp3D(s1,3,cell).EQ.asc_xvp3D(s2,3,cell2).AND.
     .            asc_yvp3D(s1,3,cell).EQ.asc_yvp3D(s2,3,cell2).AND.
     .            asc_zvp3D(s1,3,cell).EQ.asc_zvp3D(s2,3,cell2).AND.
     .            asc_xvp3D(s1,4,cell).EQ.asc_xvp3D(s2,4,cell2).AND.
     .            asc_yvp3D(s1,4,cell).EQ.asc_yvp3D(s2,4,cell2).AND.
     .            asc_zvp3D(s1,4,cell).EQ.asc_zvp3D(s2,4,cell2)) THEN

c                 WRITE(0,'(A,4I6)') 'MATCH END: ',cell,s1,cell2,s2
                 asc_link3D(s1,cell) = asc_cell(cell2)

              ENDIF

            ENDIF

          ENDDO

        ENDDO

      ENDDO


c      WRITE(0,*) 'FINDING HANGIN BOUNDARIES'

      DO cell = 1, asc_ncell
        DO s1 = 1, asc_nvp3D(cell)
          IF (asc_link3D(s1,cell).EQ.999999) asc_link3D(s1,cell) = -1
        ENDDO
      ENDDO

      DO cell = 1, asc_ncell
        DO v1 = 1, asc_nvp3D(cell)

          WRITE(SLOUT,'(A,4I7)') 'CHECKING:',
     .      cell,v1,asc_link3D(v1,cell),link_old(v1,cell)

          IF (asc_link3D(v1,cell).NE.link_old(v1,cell))
     .      WRITE(SLOUT,*) 'NOPE'

        ENDDO

      ENDDO



c      WRITE(0,*) 'WRITING TO SLOUT'
    





c      WRITE(0,*) 'DONE LINKING'


c      DO cell = 1, asc_ncell
c        WRITE(SLOUT,'(A,4I10,3F10.5)') '3DCELL:',
c     .    cell,asc_nvp(cell),ascnvertex(cell),asc_cell(cell),
c     .    asc_zmin3D(cell),asc_zmax3D(cell),asc_vol(cell)
c        DO s1= 1, asc_nvp3D(cell)
c          WRITE(SLOUT,'(6I10)') asc_link3D(s1,cell)
c          DO n1 = 1, 4
c            WRITE(SLOUT,'(5X,A,I6,2(3F10.5,2X))')
c     .     'S: ',n1,asc_xvp3D(s1,n1  ,cell),asc_yvp3D(s1,n1  ,cell),
c     .        asc_zvp3D(s1,n1  ,cell),
c     .        asc_xvp3D(s1,n1+4,cell),asc_yvp3D(s1,n1+4,cell),
c     .        asc_zvp3D(s1,n1+4,cell)
c          ENDDO
c          DO n1 = 1, ascnvertex(cell)
c            WRITE(SLOUT,'(5X,A,I6,2F10.5))') 'V: ',
c     .         n1,ascvertex(2*n1-1,cell),ascvertex(2*n1,cell)
c          ENDDO
c        ENDDO
c      ENDDO



        RETURN


      ENDIF








c      RETURN
c      IF (asc_3Dmode.NE.2) RETURN


      WRITE(SLOUT,*) 'CHECKING ASSIGNCELLLINK 1'


      DO cell = 1, asc_ncell
        DO v1 = 1, asc_nvp(cell)

          link_old(v1,cell) = asc_link(v1,cell)

        ENDDO
      ENDDO

      DO cell = 1, asc_ncell
        DO v1 = 1, asc_nvp(cell)

          WRITE(SLOUT,'(5I6)')
     .     cell,v1,asc_link(v1,cell),link_old(v1,cell),asc_nvp(cell)

          IF (asc_link(v1,cell).NE.link_old(v1,cell))
     .      WRITE(SLOUT,*) 'NOPE'

        ENDDO
      ENDDO

c...  Reinitialize cell link array:
      DO cell = 1, asc_ncell
        DO v1 = 1, asc_nvp(cell)
          IF (asc_link(v1,cell).GT.0) asc_link(v1,cell) = 999999
        ENDDO
      ENDDO


c...  For all additional cells, search every other cell and find 
c     coindicent sides.  Link these via asc_link:
      DO cell = 1, asc_ncell
c        WRITE(0,*) 'CELL:',cell
        DO cell2 = 1, asc_ncell

	  IF (cell.EQ.cell2) CYCLE

          DO v1 = 1, asc_nvp(cell)

            IF (asc_link(v1,cell).LT.0.OR.
     .          asc_link(v1,cell).NE.999999) CYCLE

            v3 = v1 + asc_nvp(cell)
	    
            DO v2 = 1, asc_nvp(cell2)
	    
              IF (asc_link(v2,cell2).LT.0.OR.
     .            asc_link(v2,cell2).NE.999999) CYCLE
	    
              v4 = v2 + asc_nvp(cell2)
             
	      IF (cell.EQ.326.AND.cell2.EQ.327) 
     .          WRITE(SLOUT,'(A,3X,4I6,4(2F10.6,2X))') 'LOOK:', 
     .            cell,cell2,v1,v2,
     .            asc_rvp(v1,cell),asc_rvp(v4,cell2),
     .            asc_zvp(v1,cell),asc_zvp(v4,cell2),
     .            asc_rvp(v3,cell),asc_rvp(v2,cell2),
     .            asc_zvp(v3,cell),asc_zvp(v2,cell2)
	    
              IF (ABS(asc_rvp(v1,cell)-asc_rvp(v4,cell2)).LT.TOL.AND.
     .            ABS(asc_zvp(v1,cell)-asc_zvp(v4,cell2)).LT.TOL.AND.
     .            ABS(asc_rvp(v3,cell)-asc_rvp(v2,cell2)).LT.TOL.AND.
     .            ABS(asc_zvp(v3,cell)-asc_zvp(v2,cell2)).LT.TOL) THEN

c              IF (asc_rvp(v1,cell).EQ.asc_rvp(v4,cell2).AND.
c     .            asc_zvp(v1,cell).EQ.asc_zvp(v4,cell2).AND.
c     .            asc_rvp(v3,cell).EQ.asc_rvp(v2,cell2).AND.
c     .            asc_zvp(v3,cell).EQ.asc_zvp(v2,cell2)) THEN
	    
                  asc_link(v1,cell) = cell2
                  asc_link(v2,cell2) = cell
	    
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO


      DO cell = 1, asc_ncell
        DO v1 = 1, asc_nvp(cell)
          IF (asc_link(v1,cell).EQ.999999) asc_link(v1,cell) = -1
        ENDDO
      ENDDO


      WRITE(SLOUT,*) 'CHECKING ASSIGNCELLLINK 2'


      DO cell = 1, asc_ncell
        DO v1 = 1, asc_nvp(cell)


          IF (asc_link(v1,cell).NE.link_old(v1,cell).AND.
     .        link_old(v1,cell).NE.999999) THEN
            WRITE(SLOUT,'(4I8,A)')
     .        cell,v1,asc_link(v1,cell),link_old(v1,cell),' NOPE'
          ELSE
            WRITE(SLOUT,'(4I8)')
     .        cell,v1,asc_link(v1,cell),link_old(v1,cell)
          ENDIF

        ENDDO
      ENDDO

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: CalcPolygonVolume
c
c a better solution would be to scan each vertex to be the common vertex,
c instead of the "center" point, for intersections with cell sides when going
c from that vertex to all of the other verticies.  If there is an intersection
c then you don't use that vertex, otherwise simple convex polygons are okay.
c
c
      SUBROUTINE CalcPolygonVolume(rp1,zp1,np,tp,v,cellnum)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      COMMON /EASMESH/ nseg,nwal,rseg,zseg,rwal,zwal
      INTEGER          nseg,nwal
      DOUBLE PRECISION rseg(MAXNKS),zseg(MAXNKS),
     .                 rwal(MAXPTS),zwal(MAXPTS)


c...  TEMP:
      INTEGER, SAVE :: count
      DATA             count /0/
c      SAVE

c...  Arguments:
      INTEGER np,tp(np),cellnum,cut1
      REAL    rp1(2*np),zp1(2*np),v
      REAL*8  rp2(2*np),zp2(2*np)

      DOUBLE PRECISION TOL
      PARAMETER       (TOL=1.0E-7)

      INTEGER CalcPoint
      REAL*8  CalcTriangleArea

      INTEGER i1,i2,i3,i4,nv,code,indbnd,tagbnd,tagtmp,cell
      LOGICAL status
      REAL    sum
      REAL*8  x(10),y(10),xc,yc,vold,rv(MAXPTS),zv(MAXPTS),t,rend,zend,
     .        area,r1,z1,rcen

c...  Move to double precision:
      DO i1 = 1, 2*np
        rp2(i1) = rp1(i1)
        zp2(i1) = zp1(i1)
      ENDDO


c      IF (.FALSE.) THEN

        i4 = 0
        DO i1 = 1, np
          IF (tp(i1).LT.0) CYCLE
          DO i2 = i1, i1 + np, np
            status = .TRUE.
            DO i3 = 1, i4
              IF (DABS(rp2(i2)-x(i3)).LT.TOL.AND.
     .            DABS(zp2(i2)-y(i3)).LT.TOL) status = .FALSE.
            ENDDO
            IF (status) THEN
              i4 = i4 + 1
              IF (i4.GT.10) CALL ER('CalcPolygonVolume',
     .                              'Array bounds exceeded',*99)
              x(i4) = rp2(i2)
              y(i4) = zp2(i2)
            ENDIF
          ENDDO
        ENDDO

c...    Find center point:
        xc = 0.0D0
        yc = 0.0D0
        DO i1 = 1, i4
          xc = xc + x(i1)
          yc = yc + y(i1)
        ENDDO
        xc = xc / DBLE(i4)
        yc = yc / DBLE(i4)

c        DO i1 = 1, i4
c          WRITE(50,*) 'VOLUME: ',x(i1),y(i1)
c        ENDDO
c        WRITE(50,*) 'VOLUME: ',xc,yc

c...Calculate area:
        v = 0
        DO i1 = 1, i4
          IF (i1.EQ.i4) THEN
            i2 = 1
          ELSE
            i2 = i1 + 1
          ENDIF

          v = v + CalcTriangleArea(x(i1),y(i1),x(i2),y(i2),xc,yc)
        ENDDO

        v = v * 2.0 * PI * rxp * eirtorfrac
c        v = v * 2.0 * 3.1415 * rxp

c        WRITE(50,*) 'VOLUME: ',v
c        WRITE( 0,*) 'VOLUME: ',v

c      ELSE
c...    New way:

c...    Been here before?
        cell = MOD(cellnum,asc_ncell)            
        IF (cell.EQ.0) cell = asc_ncell
        IF (ascnvertex(cell).GT.0) THEN
          nv = ascnvertex(cell)
          DO i1 = 1, nv
            rv(i1) = DBLE(ascvertex(2*i1-1,cell))
            zv(i1) = DBLE(ascvertex(2*i1  ,cell))
          ENDDO                  
c...      Spagetti:
          GOTO 50
        ENDIF

c...    Build list of ploygon verticies.  Cells that border the standard 
c       grid or neutral wall may have additional verticies taken from the
c       grid and wall coordinate arrays:
 
        count = count + 1

        vold = v / (2.0 * PI * rxp * eirtorfrac)

        IF (rp2(8).EQ.rp2(1).AND.rp2(5).EQ.rp2(2).AND.
     .      rp2(6).EQ.rp2(3).AND.rp2(7).EQ.rp2(4).AND.
     .      zp2(8).EQ.zp2(1).AND.zp2(5).EQ.zp2(2).AND.
     .      zp2(6).EQ.zp2(3).AND.zp2(7).EQ.zp2(4)) THEN
c...      Complete rectangular cell:

c          xc = 0.5 * (rp2(2) + rp2(1))
c          v = (rp2(1) - rp2(2)) * (zp2(3) - zp2(2)) * 2.0 * PI * xc

          nv = np
          DO i1 = 1, np
            rv(i1) = rp2(i1)
            zv(i1) = zp2(i1)
          ENDDO 
        ELSE
c...      More complicated cell that intersects the standard grid or 
c         neutral wall:

          nv = 0
          DO i1 = 1, np

            IF (tp(i1).NE.-1) THEN
c...          Start with the first qualified side:

              nv = nv + 1
              rv(nv) = rp2(i1)
              zv(nv) = zp2(i1)

              IF (i1.EQ.np) THEN
                i2 = 1
              ELSE
                i2 = i1 + 1
              ENDIF  
             
              IF (rp2(i1+4).NE.rp2(i2).OR.tp(i2).EQ.-1.OR.
     .            zp2(i1+4).NE.zp2(i2)) THEN
c...            End vertex of this side is not coincident with the 
c               starting vertex of the next side, so search for 
c               intersection with the grid or wall:                 
                nv = nv + 1
                rv(nv) = rp2(i1+4)
                zv(nv) = zp2(i1+4)
c...            Search standard grid boundary:
                indbnd = 0
                tagbnd = 0
                DO i3 = 1, nseg
                  code = CalcPoint(rseg(i3  ),zseg(i3  ),
     .                             rseg(i3+1),zseg(i3+1),
     .                             rv(nv),zv(nv),t)
                  IF (code.EQ.1.AND.t.GT.-TOL.AND.t.LT.1.0-TOL) THEN
                    IF (indbnd.EQ.0) THEN
                      indbnd = i3
                      tagbnd = 1
                    ELSE
                      CALL WN('CalcPolygonVolume','More than one '//
     .                        'intersection point found on standard '//
     .                        'grid boundary')
                    ENDIF
                  ENDIF
                ENDDO

c...            Search neutral wall bounary:
                DO i3 = 1, nwal
                  code = CalcPoint(rwal(i3  ),zwal(i3  ),
     .                             rwal(i3+1),zwal(i3+1),
     .                             rv(nv),zv(nv),t)
                  IF (code.EQ.1.AND.t.GT.-TOL.AND.t.LT.1.0-TOL) THEN
                    IF (indbnd.EQ.0) THEN
                      indbnd = i3
                      tagbnd = 2
                    ELSE
                      CALL WN('CalcPolygonVolume','More than one '//
     .                        'intersection point found on neutral '//
     .                        'grid boundary')
                    ENDIF
                  ENDIF
                ENDDO

                IF (tagbnd.NE.0) THEN
c...              Find next valid vertex in the cell:
                  DO WHILE(tp(i2).EQ.-1) 
                    IF (i2.EQ.4) THEN
                      i2 = 1
                    ELSE
                      i2 = i2 + 1
                    ENDIF
                  ENDDO
                  rend = rp2(i2)
                  zend = zp2(i2)                  
c...              Add boundary segment verticies to the polygon:
                  status = .TRUE.            
                  i4     = indbnd
                  tagtmp = tagbnd
                  DO WHILE (status)
c                    WRITE(0,*) '     SCANNING: ',i2,i4
                    IF     (tagbnd.EQ.1) THEN
c...                  Standard grid boundary:
                      code = CalcPoint(rseg(i4  ),zseg(i4  ),
     .                                 rseg(i4+1),zseg(i4+1),
     .                                 rend,zend,t)
                      IF (code.EQ.1.AND.t.GT.-TOL.AND.t.LT.1.0-TOL) THEN
                        status = .FALSE.
                      ELSE
                        nv = nv + 1
                        rv(nv) = rseg(i4+1)
                        zv(nv) = zseg(i4+1)                       
c...                    Advance the index, and switch to the neutral wall
c                       boundary segments if necessary:
                        IF (i4.EQ.nseg) THEN
                          tagbnd = 2
                          i4     = 1                         
                        ELSE
                          i4 = i4 + 1
                        ENDIF
                      ENDIF
                    ELSEIF (tagbnd.EQ.2) THEN
c...                  Neutral wall boundary:
                      code = CalcPoint(rwal(i4  ),zwal(i4  ),
     .                                 rwal(i4+1),zwal(i4+1),
     .                                 rend,zend,t)
                      IF (code.EQ.1.AND.t.GT.-TOL.AND.t.LT.1.0-TOL) THEN
                        status = .FALSE.
                      ELSE
                        nv = nv + 1
                        rv(nv) = rwal(i4+1)
                        zv(nv) = zwal(i4+1)                       
c...                    Advance the index, and switch to the standard grid
c                       boundary segments if necessary:
                        IF (i4.EQ.nwal) THEN
                          tagbnd = 1
                          i4     = 1                         
                        ELSE
                          i4 = i4 + 1
                        ENDIF
                      ENDIF
                    ELSE
                      CALL ER('CalcPolygonVolume','Lost',*99)
                    ENDIF
                  ENDDO

                ELSE
                  CALL ER('CalcPolygonVolume','No intersection with '//
     .                                        'boundary found',*99)
                ENDIF

c                WRITE(0,*) 'ODD    CELL: ',count,i1,i2,i4
              ELSE
c...            Next side is fine, so don't do anything:
c                WRITE(0,*) 'ODD    CELL: ',count,i1,i2,' CLEAR'
              ENDIF
            ENDIF          
          ENDDO
        ENDIF

        DO i1 = 1, nv
c          WRITE(0,*) '            V1:',i1,rv(i1),zv(i1)
        ENDDO

c...    Eliminate degenerate verticies:
        i3 = 1
        DO WHILE(i3.LT.nv)
          i2 = 0
          DO i1 = i3+1, nv
            IF (rv(i1).EQ.rv(i3).AND.zv(i1).EQ.zv(i3)) i2 = i1
          ENDDO
          IF (i2.GT.0) THEN
c...        Delete the degenerate point:
            DO i1 = i2, nv-1
              rv(i1) = rv(i1+1)
              zv(i1) = zv(i1+1)
            ENDDO
            nv = nv -1
          ELSE
c...        No degenerate point found so advance the index:
            i3 = i3 + 1
          ENDIF
        ENDDO

        DO i1 = 1, nv
          WRITE(SLOUT,'(A,2I4,2F12.4)')
     .       '   V2:',cellnum,i1,rv(i1),zv(i1)
        ENDDO

        IF (nv.GE.40) 
     .    CALL ER('CalcPolygonVolume','VERTEX array out of bounds',*99)

c...    Assign verticies to additional cell array:
        IF     (asc_3Dmode.EQ.2.AND.cellnum.LE.asc_ncell) THEN
          ascnvertex(cellnum) = nv
          DO i1 = 1, nv
            ascvertex(2*i1-1,cellnum) = SNGL(rv(i1))
            ascvertex(2*i1  ,cellnum) = SNGL(zv(i1))
          ENDDO                  
        ELSEIF (asc_3dmode.NE.2) THEN
          ascnvertex(cellnum) = nv
          DO i1 = 1, nv
            ascvertex(2*i1-1,cellnum) = SNGL(rv(i1))
            ascvertex(2*i1  ,cellnum) = SNGL(zv(i1))
          ENDDO
        ENDIF

50      CONTINUE

c...    Calculate volume:
        area = 0.0
        DO i1 = 1, nv
          IF (i1.EQ.nv) THEN
            i2 = 1
          ELSE
            i2 = i1 + 1
          ENDIF
          area = area + (rv(i2)*(zv(i1)+100.0D0)) -
     .                  (rv(i1)*(zv(i2)+100.0D0))  
c          area = area + (DABS(rv(i1)*zv(i2)) - DABS(rv(i2)*zv(i1)))
        ENDDO
        area = 0.5 * area


        IF     (asc_3dmode.EQ.1) THEN
          STOP 'NO LONGER SUPPORTED 0003'
          xc = (asc_zmax3D(cellnum) - asc_zmin3D(cellnum))/(2.0*PI)
        ELSEIF (asc_3dmode.EQ.2) THEN
          cut1 = -1
          DO i1 = 1, ascncut
            IF (cellnum.GE.asc_ncell*(i1-1)+1.AND.
     .          cellnum.LE.asc_ncell* i1     ) cut1 = i1
          ENDDO
          xc = (asc_zmax3D(cut1) - asc_zmin3D(cut1))/(2.0*PI)

          IF (eirntorseg.GT.0) THEN
c NEED

c...        Calculate R coordinate of the center of the cell:
            cell = MOD(cellnum,asc_ncell)            
            IF (cell.EQ.0) cell = asc_ncell
            rcen = 0.0D0
            DO i1 = 1, ascnvertex(cell)
              rcen = rcen + DBLE(ascvertex(2*i2-1,cell))
            ENDDO
            rcen = rcen / DBLE(ascnvertex(cell))

            r1 = -eirzaa / (2.0D0 * DBLE(PI))
            z1 = xc * (2.0D0 * DBLE(PI))
            
            xc = rcen * z1 / r1 

c            IF (cellnum.EQ.1) sum=0.0
c            IF (cell.EQ.2) THEN
c              sum = sum + xc   
c              WRITE(0,'(A,1P,4E12.4)') 'SUM:',sum,xc,area,xc*area
c            ENDIF

            xc = xc / (2.0D0 * DBLE(PI))

c            WRITE(0,*) 'r1,z1',r1,z1
c            WRITE(0,*) 'rcen,xc',rcen,xc
c            STOP
          ENDIF

        ELSEIF (eirzaa.EQ.-1.0) THEN
c...      Cylindrical approximation with the "toroidal circumference" based
c         on the x-point r-coordinate:
          xc = rxp * eirtorfrac
        ELSEIF (iflexopt(7).EQ.21.OR.eirzaa.LT.0.0) THEN
          IF (cellnum.EQ.1) THEN
            WRITE(0     ,*) 'USING TOROIDAL ADD CELL VOLUME'
            WRITE(PINOUT,*) 'USING TOROIDAL ADD CELL VOLUME'
          ENDIF
          xc = 0.0
          DO i1 = 1, nv
            xc = xc + rv(i1)
          ENDDO
          xc = xc / DBLE(nv) * DBLE(eirtorfrac)

c            IF (cell.EQ.2) THEN
c              WRITE(0,*) 'SUM:',xc*2.0D0*DBLE(PI),area
c            ENDIF
        ELSE
c...      "Toroidal" circumference specified in the DIVIMP input file:
          IF (cellnum.EQ.1) THEN
            WRITE(0     ,*) 'USING CUSTOM CYLINDRICAL ADD CELL VOLUME'
            WRITE(PINOUT,*) 'USING CUSTOM CYLINDRICAL ADD CELL VOLUME'
          ENDIF
          xc = eirzaa / (2.0 * PI)
        ENDIF
        
        v = area * 2.0D0 * DBLE(PI) * xc

        vold = vold * 2.0 * PI * xc

        asc_area(cellnum) = area
      
        if (asc_3dmode.eq.2) then
          WRITE(SLOUT,'(A,2I6,1P,3E12.4)') 
     .      'VOLUME CELL: ',cellnum,cut1,v,SNGL(vold),area
        else 
          WRITE(SLOUT,'(A,I6,1P,4E12.4)') 
     .      'VOLUME CELL: ',cellnum,v,SNGL(vold),area,xc*2.0*PI
        endif

c      ENDIF

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: ProcessCell
c
      SUBROUTINE ProcessCell(scel,tcel,acel,lcel,
     .                       vcel,rcel,zcel)
      IMPLICIT none

      INCLUDE 'params'

      COMMON /EASMESH/ nseg,nwal,rseg,zseg,rwal,zwal
      INTEGER          nseg,nwal
      DOUBLE PRECISION rseg(MAXNKS),zseg(MAXNKS),
     .                 rwal(MAXPTS),zwal(MAXPTS)

      INTEGER          tcel(4),acel(4),scel
      LOGICAL          lcel(4)
      DOUBLE PRECISION rcel(8),zcel(8),vcel

      LOGICAL Intersect

      INTEGER          i1,i2,n,s,nv
      LOGICAL          ldum,output
      DOUBLE PRECISION t,rv(MAXNKS),zv(MAXNKS)

      output = .FALSE.
c      output = .TRUE.

c...Prep
      DO i1 = 1, 4
        IF (i1.LT.4) THEN
          rcel(i1+4) = rcel(i1+1)
          zcel(i1+4) = zcel(i1+1)
        ELSE
          rcel(i1+4) = rcel(1)
          zcel(i1+4) = zcel(1)
        ENDIF
      ENDDO

c...  Assess intersections between the cell sides and the region boundary:
      DO i1 = 1, 4
        IF (Intersect(2,rcel(i1  ),zcel(i1  ),
     .                  rcel(i1+4),zcel(i1+4),t,n,s)) tcel(i1) = n
        
        IF (output) 
     .     WRITE(50,'(A,I6,G10.2,2I4,4F12.4)') '   -> i1,t= ',
     .       i1,t,n,s,rcel(i1),rcel(i1+4),zcel(i1),zcel(i1+4)
      ENDDO

c...  Decide if cells is outside region of interest
      IF     (tcel(1).EQ.0.AND.tcel(2).EQ.0.AND.
     .        tcel(3).EQ.0.AND.tcel(4).EQ.0.AND.
     .        .NOT.(lcel(1).OR.lcel(2).OR.lcel(3).OR.lcel(4))) THEN
        scel = -1

      ELSEIF ( tcel(1).EQ.0.AND.tcel(2).EQ.0.AND.
     .         tcel(3).EQ.0.AND.tcel(4).EQ.0.AND.
     .        (lcel(1).AND.lcel(2).AND.lcel(3).AND.lcel(4))) THEN
        scel = 0

      ELSE
c...    Well, here we have a real mess:
        scel = 0

        DO i1 = 1, 4
          IF (i1.LT.4) THEN
            i2 = i1 + 1
          ELSE
            i2 = 1
          ENDIF

          IF (output) THEN
            WRITE( 0,*) 'MES: I1,I2 = ',i1,i2,tcel(i1),lcel(i1),lcel(i2)
            WRITE(50,*) 'MES: I1,I2 = ',i1,i2,tcel(i1),lcel(i1),lcel(i2)
          ENDIF

          IF     (tcel(i1).EQ.0) THEN
            IF     (      lcel(i1).AND.lcel(i2) ) THEN
              tcel(i1) =  0
            ELSEIF (.NOT.(lcel(i1).OR. lcel(i2))) THEN
              tcel(i1) = -1
            ELSE
              CALL ER('ProcessCell','Unknown TCEL=0 event',*99)
            ENDIF

          ELSEIF (tcel(i1).EQ.1) THEN
            IF     (      lcel(i1).AND..NOT.lcel(i2))  THEN
              ldum =Intersect(2,rcel(i1  ),zcel(i1  ),
     .                          rcel(i1+4),zcel(i1+4),t,n,s)
              rcel(i1+4) = rcel(i1) + t * (rcel(i1+4) - rcel(i1))
              zcel(i1+4) = zcel(i1) + t * (zcel(i1+4) - zcel(i1))
              tcel(i1) = 1
            ELSEIF (.NOT.lcel(i1).AND.     lcel(i2)) THEN
              ldum =Intersect(2,rcel(i1  ),zcel(i1  ),
     .                          rcel(i1+4),zcel(i1+4),t,n,s)
              rcel(i1  ) = rcel(i1) + t * (rcel(i1+4) - rcel(i1))
              zcel(i1  ) = zcel(i1) + t * (zcel(i1+4) - zcel(i1))
              tcel(i1) = 2
            ELSEIF (     lcel(i1).AND.     lcel(i2)) THEN
c...          Surface passes through corner of cell, nothing to do:
              tcel(i1) = 0
            ELSEIF (.NOT.lcel(i1).AND..NOT.lcel(i2)) THEN
c...          Intersection at an acute vertex, ignore side:
              tcel(i1) = -1
            ELSE
              CALL ER('ProcessCell','Unknown TCEL=1 event',*99)
            ENDIF

          ELSEIF (tcel(i1).EQ.2) THEN
            IF     (.NOT.(lcel(i1).OR.lcel(i2)))  THEN
              ldum =Intersect(2,rcel(i1  ),zcel(i1  ),
     .                          rcel(i1+4),zcel(i1+4),t,n,s)
              rcel(i1  ) = rcel(i1) + t * (rcel(i1+4) - rcel(i1))
              zcel(i1  ) = zcel(i1) + t * (zcel(i1+4) - zcel(i1))

              ldum =Intersect(2,rcel(i1+4),zcel(i1+4),
     .                          rcel(i1  ),zcel(i1  ),t,n,s)
              rcel(i1+4) = rcel(i1+4) + t * (rcel(i1) - rcel(i1+4))
              zcel(i1+4) = zcel(i1+4) + t * (zcel(i1) - zcel(i1+4))

              tcel(i1)  = 3
            ELSE
              CALL ER('ProcessCell','Unknown TCEL=3 event',*99)
            ENDIF

          ELSE
c...Could be the wall breaks only one side of the cell.  Fixing this requires
c   that more than 2 verticies are stored for each side.  Hopefully, this will
c   not happen often, and there is always the option of resetting the additional
c   mesh parameters so that the anomalous intersection does not occur.
            CALL ER('ProcessCell','Unknown intersection event',*99)
          ENDIF
        ENDDO
c...    An intersection is registered, but all cells are outside the grid
c       so the intersection occured at a vertex.  Ignore cell:
        IF (tcel(1).EQ.-1.AND.tcel(2).EQ.-1.AND.
     .      tcel(3).EQ.-1.AND.tcel(4).EQ.-1.AND.
     .      .NOT.(lcel(1).OR.lcel(2).OR.lcel(3).OR.lcel(4)))
     .    scel = -1
      ENDIF

      RETURN
99    CONTINUE
      WRITE(0,*) 'R,Z=',rcel(1),zcel(1)

      STOP
      END
c
c ======================================================================
c
c function: Interesct
c
      LOGICAL FUNCTION Intersect(mode,r1,z1,r2,z2,t1,num,sur)
      IMPLICIT none

      INCLUDE 'params'

      COMMON /EASMESH/ nseg,nwal,rseg,zseg,rwal,zwal
      INTEGER          nseg,nwal
      DOUBLE PRECISION rseg(MAXNKS),zseg(MAXNKS),
     .                 rwal(MAXPTS),zwal(MAXPTS)

      DOUBLE PRECISION r1,z1,r2,z2,t1
      INTEGER          mode,num,sur

      DOUBLE PRECISION TOL
      PARAMETER       (TOL=1.0D-12)

      INTEGER          i1
      DOUBLE PRECISION t2,tmin

      Intersect = .FALSE.

      t1   = -1.0D0
      num  = 0
      sur  = 0

      tmin = HI

      DO i1 = 1, nseg
        CALL CalcInter(r1,z1,r2,z2,rseg(i1  ),zseg(i1  ),
     .                             rseg(i1+1),zseg(i1+1),t1,t2)

        IF ((mode.EQ.1.AND.
     .       t1.GT.0.0D0+TOL.AND.t1.LE.1.0D0+TOL.AND.
     .       t2.GT.0.0D0+TOL.AND.t2.LE.1.0D0+TOL).OR.
     .      (mode.EQ.2.AND.
     .       t1.GE.0.0D0-TOL.AND.t1.LE.1.0D0+TOL.AND.
     .       t2.GT.0.0D0+TOL.AND.t2.LE.1.0D0+TOL)) THEN
          num = num + 1
          IF (t1.LT.tmin) tmin = t1
          IF     (sur.EQ.2) THEN
            sur = 3
          ELSEIF (sur.EQ.0) THEN
            sur = 1
          ENDIF
        ENDIF
c        WRITE(50,90) '         S-->> ',i1,t1,t2,sur,num,tmin,
c     .               rseg(i1),zseg(i1),rseg(i1+1),zseg(i1+1)
90      FORMAT(A,I4,1P,2E16.8,0P,2I4,1P,E10.2,0P,4F10.4)
      ENDDO

      DO i1 = 1, nwal
        CALL CalcInter(r1,z1,r2,z2,rwal(i1  ),zwal(i1  ),
     .                             rwal(i1+1),zwal(i1+1),t1,t2)

        IF ((mode.EQ.1.AND.
     .       t1.GT.0.0D0+TOL.AND.t1.LE.1.0D0+TOL.AND.
     .       t2.GT.0.0D0+TOL.AND.t2.LE.1.0D0+TOL).OR.
     .      (mode.EQ.2.AND.
     .       t1.GE.0.0D0-TOL.AND.t1.LE.1.0D0+TOL.AND.
     .       t2.GT.0.0D0+TOL.AND.t2.LE.1.0D0+TOL)) THEN
          num = num + 1
          IF (t1.LT.tmin) tmin = t1
          IF     (sur.EQ.1) THEN
            sur = 3
          ELSEIF (sur.EQ.0) THEN
            sur = 2
          ENDIF
        ENDIF
c        WRITE(50,90) '         W-->> ',i1,t1,t2,sur,num,tmin,
c     .               rwal(i1),zwal(i1),rwal(i1+1),zwal(i1+1)
      ENDDO


      IF (num.GT.0) THEN
        t1        = tmin
        Intersect = .TRUE.
      ENDIF

      RETURN
99    STOP
      END
cc
c ======================================================================
c
c subroutine: AssembleGrid
c
c
c
c     +--------3--------+
c     |3               4|
c     |                 |
c     |                 |
c     2                 4
c     |                 |
c     |                 |
c     |2               1|
c     +--------1--------+
c
c
      SUBROUTINE AssembleGrid(nver,rver,nhor,zhor,code,region)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slcom'

      COMMON /EASMESH/ nseg,nwal,rseg,zseg,rwal,zwal
      INTEGER          nseg,nwal
      DOUBLE PRECISION rseg(MAXNKS),zseg(MAXNKS),
     .                 rwal(MAXPTS),zwal(MAXPTS)


      INTEGER          nver,nhor,region
      DOUBLE PRECISION rver(nver),zhor(nhor)

      INTEGER CalcPoint
      LOGICAL Intersect

      REAL       TOL
      PARAMETER (TOL=1.0E-07)


      INTEGER    MAXVER    ,MAXCEL
      PARAMETER (MAXVER=100,MAXCEL=5000)

      INTEGER code,i1,i2,i3,i4,i5,i6,i7,imax,imin,nadd,n,s,reg
      LOGICAL status,tagseg1,tagseg2,output
      REAL*8  rmin,zmax,dmin,t,r1,z1

      INTEGER          npoi
      DOUBLE PRECISION rpoi(MAXVER,MAXVER),zpoi(MAXVER,MAXVER)
      LOGICAL          lpoi(MAXVER,MAXVER)

      COMMON /CELCOM/  rcel,zcel,vcel,lcel,acel,tcel,scel,ncel
      DOUBLE PRECISION rcel(8,MAXCEL),zcel(8,MAXCEL),vcel(MAXCEL)
      LOGICAL          lcel(4,MAXCEL)
      INTEGER          acel(4,MAXCEL),tcel(4,MAXCEL),scel(MAXCEL),ncel

c      output = .TRUE.
      output = .FALSE.

c _cell - additional surface cell identifier
c _side -
c _link - neighbouring additional cell
c _grid - is the surface associated with a standard grid boundary surface, and
c         if so, what are the index boundaries

      CALL IZero(scel,MAXCEL)
      CALL DZero(vcel,MAXCEL)
      CALL IZero(acel,MAXCEL*4)
      CALL IZero(tcel,MAXCEL*4)
      CALL DZero(rcel,MAXCEL*8)
      CALL DZero(zcel,MAXCEL*8)

      DO i1 = 1, MAXCEL
        DO i2 = 1, 4
          lcel(i2,i1) = .FALSE.
        ENDDO
      ENDDO

      nadd = 0


      IF     (code.EQ.1) THEN
c...    PFZ:
c
c
c
c
c
c


c...    Find peak vertex (assuming there are no local maxima, which is
c       unlikely) and set center line:
        IF (output) WRITE(0,*) 'FIND CENTER LINE'
        zmax = -HI
        imax = 1
        DO i1 = 1, nseg+1
          IF (zseg(i1).GT.zmax) THEN
            zmax = MAX(zmax,zseg(i1))
            imax = i1
          ENDIF
        ENDDO
        IF (output) WRITE(SLOUT,*) 'ZMAX: ',zmax,rseg(imax)
c...Find closest vertical gridline to RSEG(IMAX) and shift grid to match:
        IF (output) WRITE(0,*) 'SHIFT GRID TO CENTER LINE'
        dmin = HI
        imin = 1
        DO i1 = 1, nver
          IF (DABS(rver(i1)-rseg(imax)).LT.DABS(dmin)) THEN
            dmin = rver(i1)-rseg(imax)
            imin = i1
          ENDIF
        ENDDO

        IF (eirbgk.EQ.2) THEN
          DO i1 = 1, nver
            rver(i1) = rver(i1) - dmin
          ENDDO
          IF (output) WRITE(SLOUT,*) 'IMIN: ',imin,rseg(imin)
c...Shift the rest of the vertical mesh lines that intersect the grid:
          DO i1 = 1, nver
            IF (i1.EQ.imin) CYCLE
            IF (Intersect(1,rver(i1),+10.0D0,
     .                      rver(i1),-10.0D0,t,n,s)) THEN
              IF (s.EQ.1.OR.s.EQ.3) THEN
                dmin = HI
                i3   = 0
                DO i2 = 1, nseg+1
                  IF (DABS(rver(i1)-rseg(i2)).LT.DABS(dmin)) THEN
                    dmin = rver(i1)-rseg(i2)
                    i3 = i2
                  ENDIF
                ENDDO
c              WRITE(0,*) 'GRID SEG = ',i3
                rver(i1) = rseg(i3)
              ENDIF

            ENDIF
c            WRITE(0,*) 'VERTICAL SEGMENT: ',i1,t,n,s
          ENDDO

c...      Shift the rest of HORIZONTAL -- DO IT __ vertical mesh lines that
c         intersect the grid:
          DO i1 = 1, nhor
            IF (Intersect(1,+10.0D0,zhor(i1),
     .                      -10.0D0,zhor(i1),t,n,s)) THEN
              IF (s.EQ.1.OR.s.EQ.3) THEN
                dmin = HI
                i3   = 0
                DO i2 = 1, nseg+1
                  IF (DABS(zhor(i1)-zseg(i2)).LT.DABS(dmin)) THEN
                    dmin = zhor(i1) - zseg(i2)
                    i3 = i2
                  ENDIF
                ENDDO
c                WRITE(0,*) 'GRID SEG = ',i3
                zhor(i1) = zseg(i3)
              ENDIF
            ENDIF
c            WRITE(0,*) 'VERTICAL SEGMENT: ',i1,t,n,s
          ENDDO
        ENDIF

c...Delete any degenerate grid lines....
        IF (output) THEN
          WRITE(0,*) '***WORK TO DO*** DELETE DEGENERATE GRID LINES'
c          WRITE(0,*) '***WORK TO DO*** IMPROVE CELL VOLUME CALCULATION'
        ENDIF





c debug
        IF (output) THEN
          WRITE(0,*) 'OUTPUT VERTICAL LINE DATA'
          DO i1 = 1, nver
            WRITE(SLOUT,*) 'VER: ',i1,rver(i1),imax,rseg(imax)
          ENDDO
          DO i1 = 1, nhor
            WRITE(SLOUT,*) 'HOR: ',i1,zhor(i1)
          ENDDO
        ENDIF
c...Construct mesh:
        DO i2 = 1, nhor
          DO i1 = 1, nver
            rpoi(i1,i2) = rver(i1)
            zpoi(i1,i2) = zhor(i2)
          ENDDO
        ENDDO
c...Determine which cells are inside the region of interest:
c...Find starting point (segment containing max. mesh z value):
        i3 = 0
        t  = 0.0D0
        i1 = imin
        DO i2 = 1, nhor-1
          IF (CalcPoint(rpoi(i1,i2  ),zpoi(i1,i2  ),
     .                  rpoi(i1,i2+1),zpoi(i1,i2+1),
     .                  rpoi(i1,i2  ),zmax         ,t).EQ.1) i3 = i2
        ENDDO
        IF (output)
     .    WRITE(SLOUT,*) 'START SEGMENT: = ',
     .      imin,i3,rpoi(i1,i3),rpoi(i1,i3+1)
        IF (i3.EQ.0) CALL ER('AssembleGrid','Cannot find starting '//
     .                                      'point',*99)

c...    Scan along starting segment (knowing that point i3, and all points above
c       that are outside the region, so set lpoi=.FALSE.)
        DO i2 = i3, 1, -1
          lpoi(i1,i2) = .FALSE.
        ENDDO

        status = .FALSE.
        DO i2 = i3, nhor-1
          IF (Intersect(1,rpoi(i1,i2  ),zpoi(i1,i2  ),
     .                    rpoi(i1,i2+1),zpoi(i1,i2+1),t,n,s)) THEN
            IF (n.GT.0.AND.MOD(n,2).EQ.1) status = .NOT.status
          ENDIF
          lpoi(i1,i2+1) = status
        ENDDO
c...Scan perpendicular to start segment
        DO i2 = 1, nhor
          status = lpoi(imin,i2)
          DO i1 = imin, 2, -1
            IF (Intersect(1,rpoi(i1  ,i2),zpoi(i1  ,i2),
     .                      rpoi(i1-1,i2),zpoi(i1-1,i2),t,n,s)) THEN
              IF (n.GT.0.AND.MOD(n,2).EQ.1) status = .NOT.status
            ENDIF
            lpoi(i1-1,i2) = status
          ENDDO

          status = lpoi(imin,i2)
          DO i1 = imin, nver-1
            IF (output) WRITE(50,*) 'SCANNING RIGHT i1 =',i1,i2,status
            IF (Intersect(1,rpoi(i1  ,i2),zpoi(i1  ,i2),
     .                      rpoi(i1+1,i2),zpoi(i1+1,i2),t,n,s)) THEN
              IF (n.GT.0.AND.MOD(n,2).EQ.1) status = .NOT.status
            ENDIF
            lpoi(i1+1,i2) = status
          ENDDO
        ENDDO
c...debug
        IF (output) THEN
          DO i1 = 1, nver
            WRITE(SLOUT,*)
            DO i2 = 1, nhor
              WRITE(SLOUT,*) 'GRID: ',i1,i2,lpoi(i1,i2),rpoi(i1,i2),
     .                       zpoi(i1,i2),lpoi(i1,i2)
            ENDDO
          ENDDO
          WRITE(0,*) 'BUILDING CELLS'
        ENDIF



c...build cells:
        i3 = 0
        DO i2 = 1, nhor-1
          DO i1 = 1, nver-1
            i3 = i3 + 1
            IF (i2.LT.nhor-1) acel(1,i3) = (i2    ) * (nver - 1) + i1
            IF (i1.GT.1     ) acel(2,i3) =  i3 - 1
            IF (i2.GT.1     ) acel(3,i3) = (i2 - 2) * (nver - 1) + i1
            IF (i1.LT.nver-1) acel(4,i3) =  i3 + 1
            rcel(1,i3) = rpoi(i1+1,i2+1)
            rcel(2,i3) = rpoi(i1  ,i2+1)
            rcel(3,i3) = rpoi(i1  ,i2  )
            rcel(4,i3) = rpoi(i1+1,i2  )
            zcel(1,i3) = zpoi(i1+1,i2+1)
            zcel(2,i3) = zpoi(i1  ,i2+1)
            zcel(3,i3) = zpoi(i1  ,i2  )
            zcel(4,i3) = zpoi(i1+1,i2  )
            lcel(1,i3) = lpoi(i1+1,i2+1)
            lcel(2,i3) = lpoi(i1  ,i2+1)
            lcel(3,i3) = lpoi(i1  ,i2  )
            lcel(4,i3) = lpoi(i1+1,i2  )
          ENDDO
        ENDDO
        ncel = i3
c...process each cell

c...debug
        IF (output) THEN
          DO i1 = 1, ncel
            WRITE(SLOUT,91) 'CELL: ',
     .        i1,scel(i1),(tcel(i2,i1),i2=1,4),(acel(i2,i1),i2=1,4),
     .                    (lcel(i2,i1),i2=1,4),
     .                    (rcel(i2,i1),zcel(i2,i1),i2=1,4),
     .                    vcel(i1)
            WRITE(SLOUT,92) 'CELL: ',
     .                    (rcel(i2,i1),zcel(i2,i1),i2=5,8)
          ENDDO
        ENDIF

        DO i1 = 1, ncel
          IF (output) THEN
c            WRITE(0    ,*) 'PROCESS: ',i1
            WRITE(SLOUT,*) 'PROCESS: ',i1
          ENDIF
          CALL ProcessCell(scel(i1),tcel(1,i1),acel(1,i1),lcel(1,i1),
     .                     vcel(i1),rcel(1,i1),zcel(1,i1))
        ENDDO


c...Delete any degenerate grid lines....

c        goto 35

        IF (eirbgk.EQ.2) THEN
c...      Move grid points that are on the grid boundary, but are not co-
c         incident with grid vertices:
          DO i1 = 1, ncel
            DO i2 = 1, 8
              DO i3 = 1, nseg
                i6 = CalcPoint(rseg(i3),zseg(i3),rseg(i3+1),zseg(i3+1),
     .                         rcel(i2,i1),zcel(i2,i1),t)

                IF (i6.EQ.1.AND.t.GT.TOL.AND.t.LT.1.0D0-TOL) THEN

                  IF (output) THEN
                    WRITE(SLOUT,*) 'FOUND ONE:',i1,i2,i3,t
                    WRITE(SLOUT,*) '   ',rcel(i2,i1),zcel(i2,i1)
                  ENDIF
c...              Scan the verticies on cell I1 to see if moving vertex I2
c                 will create a zero volume cell:
                  tagseg1 = .TRUE.
                  tagseg2 = .TRUE.
                  DO i4 = 1, 8
c                  DO i4 = i2, 8
                    IF (i2.EQ.i4) CYCLE

                    IF (DABS(rcel(i4,i1)-rseg(i3  )).LT.TOL.AND.
     .                  DABS(zcel(i4,i1)-zseg(i3  )).LT.TOL) THEN
                      IF (output) WRITE(SLOUT,*) '   SEG 1 DETECT ',i4
                      tagseg1 = .FALSE.
                    ENDIF
                    IF (DABS(rcel(i4,i1)-rseg(i3+1)).LT.TOL.AND.
     .                  DABS(zcel(i4,i1)-zseg(i3+1)).LT.TOL) THEN
                      tagseg2 = .FALSE.
                      IF (output) WRITE(SLOUT,*) '   SEG 2 DETECT ',i4
                    ENDIF
                  ENDDO
                  IF (output) WRITE(SLOUT,*) '   ',tagseg1,tagseg2
c...              Move vacuum grid vertex to match grid boundary vertex:
                  IF ((tagseg1.AND.t.LT.0.5).OR..NOT.tagseg2) THEN
                    IF (output)
     .                WRITE(SLOUT,*) '  ASSIGNING 1',rseg(i3),zseg(i3)
                    i5 = i3
                  ELSE
                    IF (output) 
     .              WRITE(SLOUT,*) '  ASSIGNING 2',rseg(i3+1),zseg(i3+1)
                    i5 = i3 + 1
                  ENDIF
c...              Check that the grid boundary vertex is available:
                  tagseg1 = .TRUE.
                  DO i6 = 1, ncel
                    DO i7 = 1, 8
                      IF (DABS(rcel(i7,i6)-rseg(i5)).LT.TOL.AND.
     .                    DABS(zcel(i7,i6)-zseg(i5)).LT.TOL) THEN
                        IF (output)
     .                    WRITE(SLOUT,*) 
     .                      '   BUMMER - VIOLATION ',i7,i6,i5
                        tagseg2 = .FALSE.
                      ENDIF
                    ENDDO
                  ENDDO
c...              Store location of vertex I2,I1:
                  r1 = rcel(i2,i1)
                  z1 = zcel(i2,i1)
c...              Move vertex I2,I1:
                  rcel(i2,i1) = rseg(i5)
                  zcel(i2,i1) = zseg(i5)
c...              Move all verticies that are coincident with I2,I1 as well:
                  DO i6 = 1, ncel
                    DO i7 = 1, 8
                      IF (DABS(rcel(i7,i6)-r1).LT.TOL.AND.
     .                    DABS(zcel(i7,i6)-z1).LT.TOL) THEN
                        IF (output) 
     .                    WRITE(SLOUT,*) '   ALSO MOVING ',i6,i7
                        rcel(i7,i6) = rseg(i5)
                        zcel(i7,i6) = zseg(i5)
                      ENDIF
                    ENDDO
                  ENDDO
                  IF (output) 
     .              WRITE(SLOUT,*) '   ',rcel(i2,i1),zcel(i2,i1)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF



c...Watch for 2+ grid surface intersections...
        IF (output)
     .    WRITE(0,*) '***WORK TO DO*** NOT CHECKING FOR MULTIPLE GRID'//
     .             ' INTERSECTIONS'

c...Check for any sides that have zero length (can happen for verticies on
c   surfaces):
        DO i1 = 1, ncel
          DO i2 = 1, 4
            IF (SQRT((rcel(i2,i1)-rcel(i2+4,i1))**2.0 +
     .               (zcel(i2,i1)-zcel(i2+4,i1))**2.0).LT.1.0D-10) THEN
              tcel(i2,i1) = -2
c              WRITE( 0,*) 'ZERO LENGTH SIDE FOUND  I1,I2= ',i1,i2
              IF (output) WRITE(50,*) 'ZERO LENGTH SIDE FOUND  I1,I2= ',
     .                                i1,i2
            ENDIF
          ENDDO
        ENDDO

c...debug
        IF (output) THEN
          DO i1 = 1, ncel
            WRITE(SLOUT,91) 'CELL2: ',
     .        i1,scel(i1),(tcel(i2,i1),i2=1,4),(acel(i2,i1),i2=1,4),
     .                    (lcel(i2,i1),i2=1,4),
     .                    (rcel(i2,i1),zcel(i2,i1),i2=1,4),
     .                    vcel(i1)
            WRITE(SLOUT,92) 'CELL2: ',
     .                    (rcel(i2,i1),zcel(i2,i1),i2=5,8)
          ENDDO
        ENDIF
91      FORMAT(A,I4,I4,1X,4I4,1X,4I4,1X,4L2,4(2F7.3,1X),1P,E10.2,0P)
92      FORMAT(A,4X,4X,1X,16X,1X,16X,1X,8X ,4(2F7.3,1X))





      ELSEIF (code.EQ.2) THEN
c...    DIII-D plenum grid:
c
c
c
c
c
c

        output = .TRUE.


c...    Shift the rest of the vertical mesh lines that intersect the grid:
        IF (eirbgk.EQ.2) THEN
          DO i1 = 1, nver
            IF (Intersect(1,rver(i1),+10.0D0,
     .                      rver(i1),-10.0D0,t,n,s)) THEN
              IF (s.EQ.1.OR.s.EQ.3) THEN
                dmin = HI
                i3   = 0
                DO i2 = 1, nseg+1
                  IF (DABS(rver(i1)-rseg(i2)).LT.DABS(dmin)) THEN
                   dmin = rver(i1)-rseg(i2)
                   i3 = i2
                  ENDIF
                ENDDO
c                WRITE(0,*) 'GRID SEG = ',i3
                rver(i1) = rseg(i3)
              ENDIF

            ENDIF
c            WRITE(0,*) 'VERTICAL SEGMENT: ',i1,t,n,s
          ENDDO

c...      Shift the rest of HORIZONTAL -- DO IT __ vertical mesh lines that
c         intersect the grid:
          DO i1 = 1, nhor
            IF (Intersect(1, 0.0D0,zhor(i1),
     .                      10.0D0,zhor(i1),t,n,s)) THEN
              IF (s.EQ.1.OR.s.EQ.3) THEN
                dmin = HI
                i3   = 0
                DO i2 = 1, nseg+1
                  IF (DABS(zhor(i1)-zseg(i2)).LT.DABS(dmin)) THEN
                    dmin = zhor(i1) - zseg(i2)
                    i3 = i2
                  ENDIF
                ENDDO
c                WRITE(0,*) 'GRID SEG = ',i3
                zhor(i1) = zseg(i3)
              ENDIF
            ENDIF
c            WRITE(0,*) 'VERTICAL SEGMENT: ',i1,t,n,s
          ENDDO
        ENDIF



c...    Find leftmost intersection between the horizontal gridlines 
c       and the standard mesh:
        rmin = HI
        imin = 0
        DO i1 = 1, nhor
          IF (Intersect(1, 0.0D0,zhor(i1),
     .                    10.0D0,zhor(i1),t,n,s)) THEN
            IF (s.EQ.1.OR.s.EQ.3) THEN
              IF (t*10.0D0.LT.rmin) THEN
c                WRITE(0,*) '-->',t * 10.0D0
	
                rmin = t * 10.0
                imin = i1
              ENDIF
            ENDIF
          ENDIF
c          WRITE(0,*) 'VERTICAL SEGMENT: ',i1,t,n,s
        ENDDO


c...Delete any degenerate grid lines....
        IF (output) THEN
          WRITE(0,*) '***WORK TO DO*** DELETE DEGENERATE GRID LINES'
c          WRITE(0,*) '***WORK TO DO*** IMPROVE CELL VOLUME CALCULATION'
        ENDIF

        IF (output) WRITE(SLOUT,*) 'RMIN: ',imin,rmin,rseg(imin)

c debug
        IF (output) THEN
          WRITE(0,*) 'OUTPUT VERTICAL LINE DATA'
          DO i1 = 1, nver
            WRITE(SLOUT,*) 'VER: ',i1,rver(i1)
          ENDDO
          DO i1 = 1, nhor
            WRITE(SLOUT,*) 'HOR: ',i1,zhor(i1),imin
          ENDDO
        ENDIF


c...    Construct mesh:

        DO i2 = 1, nhor
          DO i1 = 1, nver
            rpoi(i1,i2) = rver(i1)
            zpoi(i1,i2) = zhor(i2)
          ENDDO
        ENDDO

c...    Determine which cells are inside the region of interest:

c...    Find starting point (segment containing min. mesh r value).  Assume
c       that the middle line only cuts the standard grid segments once:
        i3 = 0
        t  = 0.0D0
        i2 = imin
        DO i1 = 1, nver-1
          IF (CalcPoint(rpoi(i1  ,i2),zpoi(i1  ,i2),
     .                  rpoi(i1+1,i2),zpoi(i1+1,i2),
     .                  rmin         ,zpoi(i1  ,i2),t).EQ.1) i3 = i2

          WRITE(SLOUT,'(A,I6,3(2F10.4,2X))') '-->',
     .      i1,rpoi(i1,i2),zpoi(i1,i2),rpoi(i1+1,i2),zpoi(i1+1,i2),
     .         rmin,zpoi(i1,i2)
        ENDDO

        IF (i3.EQ.0) CALL ER('AssembleGrid','Cannot find starting '//
     .                                      'point',*99)

        IF (output)
     .    WRITE(SLOUT,*) 'START SEGMENT: = ',
     .      imin,i3,rpoi(i3,i1),rpoi(i3+1,i1)


c...    Scan along starting segment (it is assumed that the left-most point 
c       is outside the grid):
        status = .FALSE.
        DO i1 = 1, nver-1
          IF (Intersect(1,rpoi(i1  ,i2),zpoi(i1  ,i2),
     .                    rpoi(i1+1,i2),zpoi(i1+1,i2),t,n,s)) THEN
            IF (n.GT.0.AND.MOD(n,2).EQ.1) status = .NOT.status
          ENDIF
          lpoi(i1+1,i2) = status
        ENDDO

c...    Scan perpendicular to start segment:
        DO i1 = 1, nver
          status = lpoi(i1,imin)
          DO i2 = imin, 2, -1
            IF (Intersect(1,rpoi(i1,i2  ),zpoi(i1,i2  ),
     .                      rpoi(i1,i2-1),zpoi(i1,i2-1),t,n,s)) THEN
              IF (n.GT.0.AND.MOD(n,2).EQ.1) status = .NOT.status
            ENDIF
            lpoi(i1,i2-1) = status
          ENDDO
          status = lpoi(i1,imin)
          DO i2 = imin, nhor-1
            IF (Intersect(1,rpoi(i1,i2  ),zpoi(i1,i2  ),
     .                      rpoi(i1,i2+1),zpoi(i1,i2+1),t,n,s)) THEN
              IF (n.GT.0.AND.MOD(n,2).EQ.1) status = .NOT.status
            ENDIF
            lpoi(i1,i2+1) = status
          ENDDO
        ENDDO
c...debug
        IF (output) THEN
          DO i1 = 1, nver
            WRITE(SLOUT,*)
            DO i2 = 1, nhor
              WRITE(SLOUT,'(A,2I4,L3,2F10.4)') 
     .          'GRID: ',i1,i2,lpoi(i1,i2),rpoi(i1,i2),
     .                   zpoi(i1,i2)
            ENDDO
          ENDDO
          WRITE(0,*) 'BUILDING CELLS'
        ENDIF



c  OKAY BELOW?


c...    Build cells:
        i3 = 0
        DO i2 = 1, nhor-1
          DO i1 = 1, nver-1
            i3 = i3 + 1
            IF (i2.LT.nhor-1) acel(1,i3) = (i2    ) * (nver - 1) + i1
            IF (i1.GT.1     ) acel(2,i3) =  i3 - 1
            IF (i2.GT.1     ) acel(3,i3) = (i2 - 2) * (nver - 1) + i1
            IF (i1.LT.nver-1) acel(4,i3) =  i3 + 1
            rcel(1,i3) = rpoi(i1+1,i2+1)
            rcel(2,i3) = rpoi(i1  ,i2+1)
            rcel(3,i3) = rpoi(i1  ,i2  )
            rcel(4,i3) = rpoi(i1+1,i2  )
            zcel(1,i3) = zpoi(i1+1,i2+1)
            zcel(2,i3) = zpoi(i1  ,i2+1)
            zcel(3,i3) = zpoi(i1  ,i2  )
            zcel(4,i3) = zpoi(i1+1,i2  )
            lcel(1,i3) = lpoi(i1+1,i2+1)
            lcel(2,i3) = lpoi(i1  ,i2+1)
            lcel(3,i3) = lpoi(i1  ,i2  )
            lcel(4,i3) = lpoi(i1+1,i2  )
          ENDDO
        ENDDO
        ncel = i3
c...process each cell

c...debug
        IF (output) THEN
          DO i1 = 1, ncel
            WRITE(SLOUT,91) 'CELL: ',
     .        i1,scel(i1),(tcel(i2,i1),i2=1,4),(acel(i2,i1),i2=1,4),
     .                    (lcel(i2,i1),i2=1,4),
     .                    (rcel(i2,i1),zcel(i2,i1),i2=1,4),
     .                    vcel(i1)
            WRITE(SLOUT,92) 'CELL: ',
     .                    (rcel(i2,i1),zcel(i2,i1),i2=5,8)
          ENDDO
        ENDIF


        DO i1 = 1, ncel
          IF (output) THEN
c            WRITE(0    ,*) 'PROCESS: ',i1
            WRITE(SLOUT,*) 'PROCESS: ',i1
          ENDIF
          CALL ProcessCell(scel(i1),tcel(1,i1),acel(1,i1),lcel(1,i1),
     .                     vcel(i1),rcel(1,i1),zcel(1,i1))
        ENDDO


c...Delete any degenerate grid lines....

c        goto 35

        IF (eirbgk.EQ.2) THEN
c...      Move grid points that are on the grid boundary, but are not co-
c         incident with grid vertices (surfaces?):
          DO i1 = 1, ncel
            DO i2 = 1, 8
              DO i3 = 1, nseg
                i6 = CalcPoint(rseg(i3),zseg(i3),rseg(i3+1),zseg(i3+1),
     .                         rcel(i2,i1),zcel(i2,i1),t)

                IF (i6.EQ.1.AND.t.GT.TOL.AND.t.LT.1.0D0-TOL) THEN

                  IF (output) THEN
                    WRITE(SLOUT,*) 'FOUND ONE:',i1,i2,i3,t
                    WRITE(SLOUT,*) '   ',rcel(i2,i1),zcel(i2,i1)
                  ENDIF
c...              Scan the verticies on cell I1 to see if moving vertex I2
c                 will create a zero volume cell:
                  tagseg1 = .TRUE.
                  tagseg2 = .TRUE.
                  DO i4 = 1, 8
c                  DO i4 = i2, 8
                    IF (i2.EQ.i4) CYCLE

                    IF (DABS(rcel(i4,i1)-rseg(i3  )).LT.TOL.AND.
     .                  DABS(zcel(i4,i1)-zseg(i3  )).LT.TOL) THEN
                      IF (output) WRITE(SLOUT,*) '   SEG 1 DETECT ',i4
                      tagseg1 = .FALSE.
                    ENDIF
                    IF (DABS(rcel(i4,i1)-rseg(i3+1)).LT.TOL.AND.
     .                  DABS(zcel(i4,i1)-zseg(i3+1)).LT.TOL) THEN
                      tagseg2 = .FALSE.
                      IF (output) WRITE(SLOUT,*) '   SEG 2 DETECT ',i4
                    ENDIF
                  ENDDO
                  IF (output) WRITE(SLOUT,*) '   ',tagseg1,tagseg2
c...              Move vacuum grid vertex to match grid boundary vertex:
                  IF ((tagseg1.AND.t.LT.0.5).OR..NOT.tagseg2) THEN
                    IF (output)
     .                WRITE(SLOUT,*) '  ASSIGN 1',rseg(i3),zseg(i3)
                    i5 = i3
                  ELSE
                    IF (output)
     .                WRITE(SLOUT,*) '  ASSIGN 2',rseg(i3+1),zseg(i3+1)
                    i5 = i3 + 1
                  ENDIF
c...              Check that the grid boundary vertex is available:
                  tagseg1 = .TRUE.
                  DO i6 = 1, ncel
                    DO i7 = 1, 8
                      IF (DABS(rcel(i7,i6)-rseg(i5)).LT.TOL.AND.
     .                    DABS(zcel(i7,i6)-zseg(i5)).LT.TOL) THEN
                        IF (output)
     .                    WRITE(SLOUT,*) '   VIOLATION ',i7,i6,i5
                        tagseg2 = .FALSE.
                      ENDIF
                    ENDDO
                  ENDDO
c...              Store location of vertex I2,I1:
                  r1 = rcel(i2,i1)
                  z1 = zcel(i2,i1)
c...              Move vertex I2,I1:
                  rcel(i2,i1) = rseg(i5)
                  zcel(i2,i1) = zseg(i5)
c...              Move all verticies that are coincident with I2,I1 as well:
                  DO i6 = 1, ncel
                    DO i7 = 1, 8
                      IF (DABS(rcel(i7,i6)-r1).LT.TOL.AND.
     .                    DABS(zcel(i7,i6)-z1).LT.TOL) THEN
                        IF (output) WRITE(SLOUT,*) '   MOVING ',i6,i7
                        rcel(i7,i6) = rseg(i5)
                        zcel(i7,i6) = zseg(i5)
                      ENDIF
                    ENDDO
                  ENDDO
                  IF (output) 
     .              WRITE(SLOUT,*) '   ',rcel(i2,i1),zcel(i2,i1)
                ENDIF
              ENDDO
            ENDDO
          ENDDO

        ENDIF

c 35     continue



c...Watch for 2+ grid surface intersections...
       IF (output)
     .   WRITE(0,*) '***WORK TO DO*** NOT CHECKING FOR MULTIPLE GRID'//
     .             ' INTERSECTIONS'

c...Check for any sides that have zero length (can happen for verticies on
c   surfaces):
        DO i1 = 1, ncel
          DO i2 = 1, 4
            IF (SQRT((rcel(i2,i1)-rcel(i2+4,i1))**2.0 +
     .               (zcel(i2,i1)-zcel(i2+4,i1))**2.0).LT.1.0D-10) THEN
              tcel(i2,i1) = -2
c              WRITE( 0,*) 'ZERO LENGTH SIDE FOUND  I1,I2= ',i1,i2
              IF (output) WRITE(50,*) 'ZERO LENGTH SIDE FOUND  I1,I2= ',
     .                                i1,i2
            ENDIF
          ENDDO
        ENDDO



c...debug
        IF (output) THEN
          DO i1 = 1, ncel
            WRITE(SLOUT,91) 'CELL2: ',
     .        i1,scel(i1),(tcel(i2,i1),i2=1,4),(acel(i2,i1),i2=1,4),
     .                    (lcel(i2,i1),i2=1,4),
     .                    (rcel(i2,i1),zcel(i2,i1),i2=1,4),
     .                    vcel(i1)
            WRITE(SLOUT,92) 'CELL2: ',
     .                    (rcel(i2,i1),zcel(i2,i1),i2=5,8)
          ENDDO
        ENDIF

      ENDIF




c     ******** COMMON CODE! ********



c...Assign additional surfaces:

      DO i1 = 1, ncel
        IF (scel(i1).EQ.-1) CYCLE

        asc_ncell = asc_ncell + 1
        asc_cell(asc_ncell) = asc_ncell
        asc_nvp (asc_ncell) = 4
        asc_vol (asc_ncell) = vcel(i1)
        asc_region(asc_ncell) = region
        DO i2 = 1, 4
          asc_link(i2,asc_ncell) = MIN(tcel(i2,i1),0)
        ENDDO
c...okay? - NOT!
c        DO i2 = 1, 4
c          asc_rvp(i2,asc_ncell) = rcel(i2,i1)
c          asc_zvp(i2,asc_ncell) = zcel(i2,i1)
c        ENDDO
c        DO i2 = 1, 4
c          i3 = i2 + 1
c          IF (i3.EQ.5) i3 = 1
c          asc_rvp(i2+4,asc_ncell) = asc_rvp(i3,asc_ncell)
c          asc_zvp(i2+4,asc_ncell) = asc_zvp(i3,asc_ncell)
c        ENDDO
        DO i2 = 1, 8
          asc_rvp(i2,asc_ncell) = rcel(i2,i1)
          asc_zvp(i2,asc_ncell) = zcel(i2,i1)
        ENDDO
      ENDDO

      IF (output) THEN
        DO i1 = asc_rstart(asc_nregion), asc_ncell
          WRITE(SLOUT,93) 'ASC1: ',
     .      i1,asc_cell(i1),asc_region(i1),
     .      (asc_link(i2,i1),i2=1,4),
     .      (asc_grid(i2,i1),i2=1,2),
     .      (asc_rvp (i2,i1),
     .       asc_zvp (i2,i1),i2=1,asc_nvp(i1))
          WRITE(SLOUT,94) 'ASC1: ',
     .      (asc_rvp (i2,i1),
     .       asc_zvp (i2,i1),i2=asc_nvp(i1)+1,2*asc_nvp(i1))
        ENDDO
      ENDIF

c...Find neighbours:

c...  This will have to be revamped, and moved likely, once neighbouring
c     vacuum cell regions are connected to each other:

      DO i2 = asc_rstart(asc_nregion), asc_ncell
        DO i3 = asc_rstart(asc_nregion), asc_ncell
          IF (i2.EQ.i3) CYCLE
          DO i4 = 1, asc_nvp(i2)
            IF (asc_link(i4,i2).LE.-1) CYCLE
            DO i5 = 1, asc_nvp(i3)
             i6 = i4 + asc_nvp(i2)
             i7 = i5 + asc_nvp(i3)
c             IF (i2.EQ.1) THEN
c               WRITE(0,*) i4,i6,i2,'  ',i5,i7,i3,i1,'  ',
c   .            ABS(asc_rvp(i4,i2)-asc_rvp(i7,i3)),
c   .            ABS(asc_zvp(i4,i2)-asc_zvp(i7,i3)),
c   .            ABS(asc_rvp(i6,i2)-asc_rvp(i5,i3)),
c   .            ABS(asc_zvp(i6,i2)-asc_zvp(i5,i3))
c             ENDIF
             IF (ABS(asc_rvp(i4,i2)-asc_rvp(i7,i3)).LT.TOL.AND.
     .           ABS(asc_zvp(i4,i2)-asc_zvp(i7,i3)).LT.TOL.AND.
     .           ABS(asc_rvp(i6,i2)-asc_rvp(i5,i3)).LT.TOL.AND.
     .           ABS(asc_zvp(i6,i2)-asc_zvp(i5,i3)).LT.TOL
     .         ) THEN
               IF (asc_link(i4,i2).EQ.0) THEN
                 asc_link(i4,i2) = asc_cell(i3)
               ELSE
c                 CALL WN('AssembleGrid','Crap 01')
                 WRITE(SLOUT,*) 'I2,I4,I3,I5= ',i2,i4,i3,i5,
     .                 asc_link(i4,i2)
                 CALL ER('AssembleGrid','Crap 01',*99)
               ENDIF
             ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO

c...  Calculate cell volume (assume a reasonably concave surface):

c...  TEMP!
c      asc_ncell = 18

      DO i1 = asc_rstart(asc_nregion), asc_ncell
        CALL CalcPolygonVolume(asc_rvp (1,i1),asc_zvp(1,i1),4,
     .                         asc_link(1,i1),asc_vol(i1),i1)
      ENDDO





      IF (output) THEN
        DO i1 = asc_rstart(asc_nregion), asc_ncell
          WRITE(SLOUT,93) 'ASC2: ',
     .      i1,asc_cell(i1),asc_region(i1),
     .      (asc_link(i2,i1),i2=1,4),
     .      (asc_grid(i2,i1),i2=1,2),
     .      (asc_rvp (i2,i1),
     .       asc_zvp (i2,i1),i2=1,asc_nvp(i1))
          WRITE(SLOUT,94) 'ASC2: ',
     .      (asc_rvp (i2,i1),
     .       asc_zvp (i2,i1),i2=asc_nvp(i1)+1,2*asc_nvp(i1))
        ENDDO
        WRITE(0,*) 'REG,ASC_NCELL= ',asc_ncell
      ENDIF






      RETURN
90    FORMAT(A,I6,2F12.8)

93    FORMAT(A,3I4,1X,4I4,1X,2I4,1X,10(2F10.6,1X:))
94    FORMAT(A,12X,1X,16X,1X, 8X,1X,10(2F10.6,1X:))


99    STOP
      END
c
c ======================================================================
c
c subroutine: BuildVacuumGrid
c
      SUBROUTINE BuildVacuumGrid
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER i1,ir,id,code,seg1
      LOGICAL mode,output
      REAL    r1,r2,z1,z2,zcut

      COMMON /WRVACCOM/ convertz
      LOGICAL           convertz

      COMMON /EASMESH/ nseg,nwal,rseg,zseg,rwal,zwal
      INTEGER          nseg,nwal
      DOUBLE PRECISION rseg(MAXNKS),zseg(MAXNKS),
     .                 rwal(MAXPTS),zwal(MAXPTS)

      INTEGER nver,nhor,ikmark,inmark,ik,i2,i3,cell,cell2
      LOGICAL status
      REAL    sver,iver,ever,rtmp,rval,
     .        shor,ihor,ehor,ztmp,zval,
     .        xcutmin,xcutmax,ycutmin,ycutmax,sum
c...MAX something here...
      REAL*8  rver(100),zhor(100)
      REAL*8  a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd,tmin

      output = .FALSE.

     

      IF (output) THEN
        WRITE(0,*)
        WRITE(0,*) 'CALCULATING ADDITIONAL SURFACE GEOMETRY'
        WRITE(0,*)
      ENDIF

      convertz = .FALSE.

      asc_nregion = 0
      asc_ncell   = 0
      asc_3Dmode  = 0
      ascncut     = 1
      CALL IZero(asc_rstart,10)
      CALL IZero(asc_rend  ,10)
      CALL IZero(asc_cell  ,MAXASC)
      CALL IZero(asc_region,MAXASC)
      CALL IZero(asc_link  ,MAXASC*4)
      CALL IZero(asc_grid  ,MAXASC*2)
      CALL IZero(asc_nvp   ,MAXASC)
      CALL RZero(asc_rvp   ,MAXASC*8)
      CALL RZero(asc_zvp   ,MAXASC*8)
      CALL RZero(asc_vol   ,MAXASCDAT)

      CALL RZero(asc_zmin3d,MAXASC)
      CALL RZero(asc_zmax3d,MAXASC)


c...  Checks:
      IF (nks(irtrap+1)+1.GT.MAXNKS) CALL ER('BVG','Bad 1',*99)
      IF (wallpts      +1.GT.MAXPTS) CALL ER('BVG','Bad 2',*99)


      CALL IZero(vacregion,MAXVACREGION)



      IF (stopopt2.EQ.21) THEN
        WRITE(PINOUT,*) 'LOADING VACUUM GRID'
 
        CALL ProcessVacuumGrid('LOAD')
   
        RETURN
      ENDIF



      IF     (stopopt2.NE.9.AND.stopopt2.NE.20) THEN
c...    PFZ:

        WRITE(0,*) 'THIS METHOD GENERATING THE VACUUM '//
     .             'REGION IS DEFUNCT'
        STOP 'STOP: HALTING CODE IN BUILDVACUUMGRID'

        code = 1

        vacregion = 1

        asc_nregion             = asc_nregion + 1
        asc_rstart(asc_nregion) = asc_ncell + 1

        ir = irtrap + 1
        DO i1 = 1, nks(ir)
          id = korpg(i1,ir)
          rseg(i1) = rvertp(1,id)
          zseg(i1) = zvertp(1,id)
        ENDDO
        rseg(i1) = rvertp(4,id)
        zseg(i1) = zvertp(4,id)
        nseg = nks(ir)
c...Extract relevant wall segments.
c...note: Expecting a continuous wall.
        mode = .FALSE.
        nwal = 0
        i1   = 0
        DO WHILE (.NOT.(.NOT.mode.AND.nwal.GT.0))
          i1 = i1 + 1
          IF (i1.GT.wallpts) i1 = 1

          IF (thesis) THEN
            IF     (wallpt2(i1,1).EQ.rvertp(4,korpg(nks(ir),ir)).AND.
     .              wallpt2(i1,2).EQ.zvertp(4,korpg(nks(ir),ir))) THEN
              mode = .TRUE.
            ELSEIF (wallpt2(i1,1).EQ.rvertp(1,korpg(1      ,ir)).AND.
     .              wallpt2(i1,2).EQ.zvertp(1,korpg(1      ,ir))) THEN
              mode = .FALSE.
              rwal(nwal+1) = wallpt2(i1,1)
              zwal(nwal+1) = wallpt2(i1,2)
            ENDIF

            IF (mode) THEN
              nwal = nwal + 1
              rwal(nwal) = wallpt2(i1,1)
              zwal(nwal) = wallpt2(i1,2)
            ENDIF

          ELSE
c...        DIII-D:
            IF     (wallpt(i1,20).EQ.rvertp(4,korpg(nks(ir),ir)).AND.
     .              wallpt(i1,21).EQ.zvertp(4,korpg(nks(ir),ir))) THEN
              mode = .TRUE.
            ELSEIF (wallpt(i1,20).EQ.rvertp(1,korpg(1      ,ir)).AND.
     .              wallpt(i1,21).EQ.zvertp(1,korpg(1      ,ir))) THEN
              mode = .FALSE.
              rwal(nwal+1) = wallpt(i1,20)
              zwal(nwal+1) = wallpt(i1,21)
            ENDIF

            IF (mode) THEN
              nwal = nwal + 1
              rwal(nwal) = wallpt(i1,20)
              zwal(nwal) = wallpt(i1,21)
            ENDIF
            IF (output) WRITE(0,*) 'WALL: ',i1,mode,nwal,wallpts
          ENDIF
        ENDDO

        IF (thesis) THEN
          IF (nrs.EQ.6) THEN
            nver = 7
            rver( 1) = 0.450
            rver( 2) = 0.528
            rver( 3) = 0.540
            rver( 4) = 0.550
            rver( 5) = 0.560
            rver( 6) = 0.572
            rver( 7) = 0.650

            nhor = 11
            zhor( 1) = -0.150
            zhor( 2) = -0.370
            zhor( 3) = -0.400
            zhor( 4) = -0.430
            zhor( 5) = -0.460
            zhor( 6) = -0.500
            zhor( 7) = -0.530
            zhor( 8) = -0.550
            zhor( 9) = -0.590
            zhor(10) = -0.620
            zhor(11) = -0.750
c            nver = 5
c            rver( 1) = 0.600
c            rver( 2) = 0.636
c            rver( 3) = 0.650
c            rver( 4) = 0.664
c            rver( 5) = 0.700
cc            nver = 3
cc            rver( 1) = 0.600
cc            rver( 2) = 0.650
cc            rver( 3) = 0.700
c
c            nhor = 10
c            zhor( 1) =  0.100
c            zhor( 2) = -0.100
c            zhor( 3) = -0.125
c            zhor( 4) = -0.175
c            zhor( 5) = -0.200
c            zhor( 6) = -0.225
c            zhor( 7) = -0.260
c            zhor( 8) = -0.300
c            zhor( 9) = -0.350
c            zhor(10) = -0.500
          ELSE
            nver = 16
            rver( 1) = 0.440
            rver( 2) = 0.460
            rver( 3) = 0.480
            rver( 4) = 0.500
            rver( 5) = 0.520
            rver( 6) = 0.540
            rver( 7) = 0.560
            rver( 8) = 0.580
            rver( 9) = 0.600
            rver(10) = 0.620
            rver(11) = 0.640
            rver(12) = 0.660
            rver(13) = 0.680
            rver(14) = 0.700
            rver(15) = 0.720
            rver(16) = 0.740

            nhor = 18
            zhor( 1) = -0.400
            zhor( 2) = -0.420
            zhor( 3) = -0.440
            zhor( 4) = -0.460
            zhor( 5) = -0.480
            zhor( 6) = -0.500
            zhor( 7) = -0.520
            zhor( 8) = -0.540
            zhor( 9) = -0.560
            zhor(10) = -0.580
            zhor(11) = -0.600
            zhor(12) = -0.620
            zhor(13) = -0.640
            zhor(14) = -0.660
            zhor(15) = -0.680
            zhor(16) = -0.700
            zhor(17) = -0.720
            zhor(18) = -0.720
          ENDIF

        ELSE
c...      DIII-D (sort of specific to sonnet.oskn_b03 at the moment):
          IF (.TRUE.) THEN
c...        For "realistic" pressure gauge geometry:

            nver = 12
            rver( 1) = 1.050
            rver( 2) = 1.130
            rver( 3) = 1.170
            rver( 4) = 1.200
            rver( 5) = 1.250
            rver( 6) = 1.270
            rver( 7) = 1.330
            rver( 8) = 1.360
            rver( 9) = 1.390
            rver(10) = 1.420
            rver(11) = 1.440
            rver(12) = 1.500

            nhor = 15
            zhor( 1) = -1.200
            zhor( 2) = -1.265
            zhor( 3) = -1.285
            zhor( 4) = -1.305
            zhor( 5) = -1.325
            zhor( 6) = -1.345
            zhor( 7) = -1.3662  
            zhor( 8) = -1.375
            zhor( 9) = -1.380
            zhor(10) = -1.385
            zhor(11) = -1.390
            zhor(12) = -1.395
            zhor(13) = -1.410
            zhor(14) = -1.415
            zhor(15) = -1.800
          ELSE          
c...        For straight PFZ pressure gauge tube wall geometry:
            nver = 14
            rver( 1) = 1.050
            rver( 2) = 1.110
            rver( 3) = 1.115
            rver( 4) = 1.150
            rver( 5) = 1.175
            rver( 6) = 1.225
            rver( 7) = 1.275
            rver( 8) = 1.300
            rver( 9) = 1.350
            rver(10) = 1.375
            rver(11) = 1.400
            rver(12) = 1.425
            rver(13) = 1.465
            rver(14) = 1.500

            nhor = 14
            zhor( 1) = -1.200
            zhor( 2) = -1.275
            zhor( 3) = -1.300
            zhor( 4) = -1.330
            zhor( 5) = -1.360
            zhor( 6) = -1.400
            zhor( 7) = -1.450
            zhor( 8) = -1.500
            zhor( 9) = -1.550
            zhor(10) = -1.600
            zhor(11) = -1.650
            zhor(12) = -1.700
            zhor(13) = -1.750
            zhor(14) = -1.800
          ENDIF
        ENDIF

c...debug
        IF (output) THEN
          DO i1 = 1, nseg+1
            WRITE(SLOUT,90) 'SEG: ',i1,rseg(i1),zseg(i1)
          ENDDO
          DO i1 = 1, nwal+1
            WRITE(SLOUT,90) 'WAL: ',i1,rwal(i1),zwal(i1)
          ENDDO
        ENDIF

        CALL AssembleGrid(nver,rver,nhor,zhor,code)

        asc_rend(asc_nregion) = asc_ncell 

      ENDIF


c-- chunk 1

      IF (stopopt2.EQ.9.OR.stopopt2.EQ.10) THEN
c...    DIII-D plenum:
c
c
c
c
c
c
c
c
c

        WRITE(0,*) 'THIS METHOD GENERATING THE VACUUM '//
     .             'REGION IS DEFUNCT'
        STOP 'STOP: HALTING CODE IN BUILDVACUUMGRID'



        asc_nregion             = asc_nregion + 1
        asc_rstart(asc_nregion) = asc_ncell + 1


        code = 2

        vacregion = vacregion + 2

c        output = .TRUE.

c...    Find plenum:
        STOP 'PFZ BREAK F'
        ir   = irbreak - 1
        nseg = 0
        DO i1 = nks(ir), ikbreak(ir), -1
          id = korpg(i1,ir)

          nseg = nseg + 1
          rseg(nseg) = rvertp(3,id)
          zseg(nseg) = zvertp(3,id)
        ENDDO
        nseg = nseg - 1

c...    Extract relevant wall segments:
        mode = .FALSE.
        nwal = 0
        i1   = 0
        r1   = rseg(1)
        z1   = zseg(1)
        r2   = rseg(nseg+1)
        z2   = zseg(nseg+1)

c
c     jdemod - this code is only executed for .not.thesis - again trying to obtain the same
c              functionality without causing a compiler bug - just commenting out the 
c              internal block executes the code if thesis is true
c             
c
        if (.not.thesis) then 

        DO WHILE (mode.OR.nwal.EQ.0)
          i1 = i1 + 1
          IF (i1.GT.nvesm+nvesp) i1 = 1

c
c jdemod - removed blank IF statement since it seemed to be exercising some pgf compiler bug
c
c          IF (thesis) THEN
c
c          ELSE
c
c          if (.not.thesis) then 
c...        DIII-D:
            IF     ((rvesm(i1,1).EQ.r2).AND.
     .              (zvesm(i1,1).EQ.z2)) THEN
              mode = .TRUE.
            ELSEIF ((rvesm(i1,2).EQ.r1).AND.
     .              (zvesm(i1,2).EQ.z1)) THEN
              mode = .FALSE.
              rwal(nwal+1) = rvesm(i1,2)
              zwal(nwal+1) = zvesm(i1,2)
            ENDIF
c
            IF (mode) THEN
              nwal = nwal + 1
              rwal(nwal) = rvesm(i1,1)
              zwal(nwal) = zvesm(i1,1)
            ENDIF

            IF (output) WRITE(0,*) 'WALL: ',i1,mode,nwal,nvesm+nvesp

c          ENDIF

        ENDDO

        endif


        nver = 15
        rver( 1) = 1.600
        rver( 2) = 1.675
        rver( 3) = 1.700
        rver( 4) = 1.725
        rver( 5) = 1.750
        rver( 6) = 1.775
        rver( 7) = 1.800
        rver( 8) = 1.850
        rver( 9) = 1.900
        rver(10) = 1.950
        rver(11) = 2.000
        rver(12) = 2.050
        rver(13) = 2.100
        rver(14) = 2.127
        rver(15) = 2.200        
	
        nhor = 14
        zhor( 1) = -0.950
        zhor( 2) = -1.050
        zhor( 3) = -1.100
        zhor( 4) = -1.150
        zhor( 5) = -1.200
        zhor( 6) = -1.250
        zhor( 7) = -1.300
        zhor( 8) = -1.325
        zhor( 9) = -1.337
        zhor(10) = -1.350
        zhor(11) = -1.360
        zhor(12) = -1.375
        zhor(13) = -1.400
        zhor(14) = -1.500

c...debug
        IF (output) THEN
          DO i1 = 1, nseg+1
            WRITE(SLOUT,90) 'SEG: ',i1,rseg(i1),zseg(i1)
          ENDDO
          DO i1 = 1, nwal+1
            WRITE(SLOUT,90) 'WAL: ',i1,rwal(i1),zwal(i1)
          ENDDO
        ENDIF

        CALL AssembleGrid(nver,rver,nhor,zhor,code)

        asc_rend(asc_nregion) = asc_ncell 

      ENDIF


c--- chunk 1 end



      
      IF (stopopt2.EQ.20) THEN



        WRITE(0,*) 'RIGHT PLACE'

        output = .TRUE.
c
c
c
c
c
c
c
c
c...    PFZ vacuum grid with the grid data taken from the DIVIMP 
c       input file:
        nhor = 0
        nver = 0
        shor = 0.0
        ihor = 0.0
        ehor = 0.0
        sver = 0.0
        iver = 0.0
        ever = 0.0
        DO i1 = 1, vacnseg
          IF (vacseg(i1,1).EQ.1.0) THEN
c...        Assign grid segments for the PFZ region:
            IF     (vacseg(i1,4).EQ.0.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                nver       = nver + 1
                rver(nver) = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                nhor       = nhor + 1
                zhor(nhor) = vacseg(i1,3)
              ENDIF
            ELSEIF (vacseg(i1,4).EQ.1.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                sver = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                shor = vacseg(i1,3)
              ENDIF
            ELSEIF (vacseg(i1,4).EQ.2.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                iver = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                ihor = vacseg(i1,3)
              ENDIF
            ELSEIF (vacseg(i1,4).EQ.3.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                ever = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                ehor = vacseg(i1,3)
              ENDIF
            ENDIF
          ENDIF
        ENDDO

        IF (sver.NE.0.0.AND.iver.NE.0.0.AND.ever.NE.0.0) THEN
c...      Build list of vertical segments:       
          rval = sver
          DO WHILE (rval.LT.ever) 
            nver       = nver + 1
            rver(nver) = rval
            rval       = rval + iver
          ENDDO
          nver       = nver + 1
          rver(nver) = ever

c...      Sort the segments from smallest R-value to largest:
          status = .TRUE.
          DO WHILE (status)
            status = .FALSE.
            DO i1 = 1, nver-1
              IF (rver(i1).GT.rver(i1+1)) THEN
                rtmp       = rver(i1)
                rver(i1)   = rver(i1+1)
                rver(i1+1) = rtmp
                status     = .TRUE.
              ENDIF
            ENDDO
          ENDDO
        ENDIF

        IF (shor.NE.0.0.AND.ihor.NE.0.0.AND.shor.NE.0.0) THEN
c...      Build list of horizontal segments:
          zval = shor
          DO WHILE (zval.GT.ehor) 
            nhor       = nhor + 1
            zhor(nhor) = zval
            zval       = zval + ihor
          ENDDO
          nhor       = nhor + 1
          zhor(nhor) = ehor

c...      Sort the segments from smallest R-value to largest:
          status = .TRUE.
          DO WHILE (status)
            status = .FALSE.
            DO i1 = 1, nhor-1
              IF (zhor(i1).LT.zhor(i1+1)) THEN
                ztmp       = zhor(i1)
                zhor(i1)   = zhor(i1+1)
                zhor(i1+1) = ztmp
                status     = .TRUE.
              ENDIF
            ENDDO
          ENDDO
        ENDIF

        IF (nver.GT.0.AND.nhor.GT.0) THEN    
c...      PFZ region vacuum grid lines exist in the input file,
c         so continue:

          code = 1

          vacregion(1) = 1

          IF (output) WRITE(0,*) 'BUILDING PFZ GRID'

c...      Extract segments from the standard grid:
          ir = irtrap + 1
          DO i1 = 1, nks(ir)
            id = korpg(i1,ir)
            rseg(i1) = rvertp(1,id)
            zseg(i1) = zvertp(1,id)
          ENDDO
          rseg(i1) = rvertp(4,id)
          zseg(i1) = zvertp(4,id)
          nseg = nks(ir)

c...      Extract relevant wall segments (expecting a continuous wall):
          mode = .FALSE.
          nwal = 0
          i1   = 0
          DO WHILE (mode.OR.nwal.EQ.0)
            i1 = i1 + 1
            IF (i1.GT.wallpts) i1 = 1

            IF     (wallpt(i1,20).EQ.rvertp(4,korpg(nks(ir),ir)).AND.
     .              wallpt(i1,21).EQ.zvertp(4,korpg(nks(ir),ir))) THEN
              mode = .TRUE.
            ELSEIF (wallpt(i1,20).EQ.rvertp(1,korpg(1      ,ir)).AND.
     .              wallpt(i1,21).EQ.zvertp(1,korpg(1      ,ir))) THEN
              mode = .FALSE.
              rwal(nwal+1) = wallpt(i1,20)
              zwal(nwal+1) = wallpt(i1,21)
            ENDIF
	  
            IF (mode) THEN
              nwal = nwal + 1
              rwal(nwal) = wallpt(i1,20)
              zwal(nwal) = wallpt(i1,21)
            ENDIF
          ENDDO

          IF (output) THEN
            DO i1 = 1, nver
              WRITE(SLOUT,90) 'VER-PFZ: ',i1,rver(i1)
            ENDDO           
            DO i1 = 1, nhor
              WRITE(SLOUT,90) 'HOR-PFZ: ',i1,zhor(i1)
            ENDDO           

            DO i1 = 1, nseg+1
              WRITE(SLOUT,90) 'SEG-PFZ: ',i1,rseg(i1),zseg(i1)
            ENDDO
            DO i1 = 1, nwal+1
              WRITE(SLOUT,90) 'WAL-PFZ: ',i1,rwal(i1),zwal(i1)
            ENDDO
          ENDIF

          asc_nregion             = asc_nregion + 1
          asc_rstart(asc_nregion) = asc_ncell + 1

          CALL AssembleGrid(nver,rver,nhor,zhor,code,1)

          asc_rend(asc_nregion) = asc_ncell 
        ENDIF



c
c
c
c
c
c
c
c
c...    DIII-D plenum vacuum grid with the grid data taken from the DIVIMP 
c       input file:
        nhor = 0
        nver = 0
        shor = 0.0
        ihor = 0.0
        ehor = 0.0
        sver = 0.0
        iver = 0.0
        ever = 0.0
        DO i1 = 1, vacnseg
          IF (vacseg(i1,1).EQ.2.0) THEN
c...        Assign grid segments for the PFZ region:
            IF     (vacseg(i1,4).EQ.0.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                nver       = nver + 1
                rver(nver) = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                nhor       = nhor + 1
                zhor(nhor) = vacseg(i1,3)
              ENDIF
            ELSEIF (vacseg(i1,4).EQ.1.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                sver = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                shor = vacseg(i1,3)
              ENDIF
            ELSEIF (vacseg(i1,4).EQ.2.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                iver = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                ihor = vacseg(i1,3)
              ENDIF
            ELSEIF (vacseg(i1,4).EQ.3.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                ever = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                ehor = vacseg(i1,3)
              ENDIF
            ENDIF
          ENDIF
        ENDDO

        IF (sver.NE.0.0.AND.iver.NE.0.0.AND.ever.NE.0.0) THEN
c...      Build list of vertical segments:       
          rval = sver
          DO WHILE (rval.LT.ever) 
            nver       = nver + 1
            rver(nver) = rval
            rval       = rval + iver
          ENDDO
          nver       = nver + 1
          rver(nver) = ever

c...      Sort the segments from smallest R-value to largest:
          status = .TRUE.
          DO WHILE (status)
            status = .FALSE.
            DO i1 = 1, nver-1
              IF (rver(i1).GT.rver(i1+1)) THEN
                rtmp       = rver(i1)
                rver(i1)   = rver(i1+1)
                rver(i1+1) = rtmp
                status     = .TRUE.
              ENDIF
            ENDDO
          ENDDO
        ENDIF

        IF (shor.NE.0.0.AND.ihor.NE.0.0.AND.shor.NE.0.0) THEN
c...      Build list of horizontal segments:
          zval = shor
          DO WHILE (zval.GT.ehor) 
            nhor       = nhor + 1
            zhor(nhor) = zval
            zval       = zval + ihor
          ENDDO
          nhor       = nhor + 1
          zhor(nhor) = ehor

c...      Sort the segments from smallest R-value to largest:
          status = .TRUE.
          DO WHILE (status)
            status = .FALSE.
            DO i1 = 1, nhor-1
              IF (zhor(i1).LT.zhor(i1+1)) THEN
                ztmp       = zhor(i1)
                zhor(i1)   = zhor(i1+1)
                zhor(i1+1) = ztmp
                status     = .TRUE.
              ENDIF
            ENDDO
          ENDDO
        ENDIF

        IF (nver.GT.0.AND.nhor.GT.0) THEN    
c...      Plenum region vacuum grid lines exist in the input file,
c         so continue:

          code = 2

          vacregion(2) = 1

          IF (output) WRITE(0,*) 'BUILDING PLENUM GRID',irbreak

c...      Extract segments from the standard grid:
          ir   = irbreak - 1
          nseg = 0
c...TROUBLE:
          DO i1 = nks(ir), ikbreak(ir), -1
            id = korpg(i1,ir)
            nseg = nseg + 1
            rseg(nseg) = rvertp(3,id)
            zseg(nseg) = zvertp(3,id)
          ENDDO
          nseg = nseg - 1

c          IF (output) WRITE(0,*) 'BUILDING PLENUM GRID A'

c...      Extract relevant wall segments:
          mode = .FALSE.
          nwal = 0
          i1   = 0
          r1   = rseg(1)
          z1   = zseg(1)
          r2   = rseg(nseg+1)
          z2   = zseg(nseg+1)
          DO WHILE (mode.OR.nwal.EQ.0)
            i1 = i1 + 1

            IF (i1.GT.nvesm+nvesp.AND.nwal.EQ.0) 
     .        CALL ER('BuildVacuumGrid','Unable to locate wall '//
     .                                  'connection',*99)

            IF (i1.GT.nvesm+nvesp) i1 = 1

            IF     (rvesm(i1,1).EQ.r2.AND.
     .              zvesm(i1,1).EQ.z2) THEN
              mode = .TRUE.
            ELSEIF (rvesm(i1,1).EQ.r1.AND.
     .              zvesm(i1,1).EQ.z1) THEN
              mode = .FALSE.
              rwal(nwal+1) = rvesm(i1,1)
              zwal(nwal+1) = zvesm(i1,1)
            ENDIF
	    IF (mode) THEN
              nwal = nwal + 1
              rwal(nwal) = rvesm(i1,1)
              zwal(nwal) = zvesm(i1,1)
            ENDIF

c            WRITE(0,'(A,2I6,L3,4F16.8)') 
c     .        '-->',i1,nwal,mode,rvesm(i1,1),r2,zvesm(i1,1),z2

          ENDDO

c          IF (output) WRITE(0,*) 'BUILDING PLENUM GRID B'

          IF (output) THEN
            DO i1 = 1, nver
              WRITE(SLOUT,90) 'VER-PLE: ',i1,rver(i1)
            ENDDO           
            DO i1 = 1, nhor
              WRITE(SLOUT,90) 'HOR-PLE: ',i1,zhor(i1)
            ENDDO           

            DO i1 = 1, nseg+1
              WRITE(SLOUT,90) 'SEG-PLE: ',i1,rseg(i1),zseg(i1)
            ENDDO
            DO i1 = 1, nwal+1
              WRITE(SLOUT,90) 'WAL-PLE: ',i1,rwal(i1),zwal(i1)
            ENDDO
          ENDIF

          asc_nregion             = asc_nregion + 1
          asc_rstart(asc_nregion) = asc_ncell + 1

c          IF (output) WRITE(0,*) 'BUILDING PLENUM GRID C'

          CALL AssembleGrid(nver,rver,nhor,zhor,code,2)

c          IF (output) WRITE(0,*) 'BUILDING PLENUM GRID D'

          asc_rend(asc_nregion) = asc_ncell 
        ENDIF
c
c
c
c
c
c
c
c
c
c
c...    Outer SOL (region not completely filled)
        nhor = 0
        nver = 0
        shor = 0.0
        ihor = 0.0
        ehor = 0.0
        sver = 0.0
        iver = 0.0
        ever = 0.0
        DO i1 = 1, vacnseg
          IF (vacseg(i1,1).EQ.3.0) THEN
c...        Assign grid segments for the PFZ region:
            IF     (vacseg(i1,4).EQ.0.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                nver       = nver + 1
                rver(nver) = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                nhor       = nhor + 1
                zhor(nhor) = vacseg(i1,3)
              ENDIF
            ELSEIF (vacseg(i1,4).EQ.1.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                sver = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                shor = vacseg(i1,3)
              ENDIF
            ELSEIF (vacseg(i1,4).EQ.2.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                iver = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                ihor = vacseg(i1,3)
              ENDIF
            ELSEIF (vacseg(i1,4).EQ.3.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                ever = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                ehor = vacseg(i1,3)
              ENDIF
            ENDIF
          ENDIF
        ENDDO

        IF (sver.NE.0.0.AND.iver.NE.0.0.AND.ever.NE.0.0) THEN
c...      Build list of vertical segments:       
          rval = sver
          DO WHILE (rval.LT.ever) 
            nver       = nver + 1
            rver(nver) = rval
            rval       = rval + iver
          ENDDO
          nver       = nver + 1
          rver(nver) = ever

c...      Sort the segments from smallest R-value to largest:
          status = .TRUE.
          DO WHILE (status)
            status = .FALSE.
            DO i1 = 1, nver-1
              IF (rver(i1).GT.rver(i1+1)) THEN
                rtmp       = rver(i1)
                rver(i1)   = rver(i1+1)
                rver(i1+1) = rtmp
                status     = .TRUE.
              ENDIF
            ENDDO
          ENDDO
        ENDIF



        IF (shor.NE.0.0.AND.ihor.NE.0.0.AND.shor.NE.0.0) THEN
c...      Build list of horizontal segments:
          zval = shor
          DO WHILE (zval.GT.ehor) 
            nhor       = nhor + 1
            zhor(nhor) = zval
            zval       = zval + ihor
          ENDDO
          nhor       = nhor + 1
          zhor(nhor) = ehor

c...      Sort the segments from smallest Z-value to largest:
          status = .TRUE.
          DO WHILE (status)
            status = .FALSE.
            DO i1 = 1, nhor-1
              IF (zhor(i1).LT.zhor(i1+1)) THEN
                ztmp       = zhor(i1)
                zhor(i1)   = zhor(i1+1)
                zhor(i1+1) = ztmp
                status     = .TRUE.
              ENDIF
            ENDDO
          ENDDO
        ENDIF


        IF (nver.GT.0.AND.nhor.GT.0) THEN    
c...      Outer SOL region vacuum grid lines exist in the input file,
c         so continue:

c...      JET grids are inverted, so the following will not work (see
c         'JET' tage below):
          IF (cgridopt.EQ.0) CALL ER('BuildVacuumGrid',
     .                               'Broken for JET grids',*99)

          code = 2

          vacregion(4) = 1

          IF (output) WRITE(0,*) 'BUILDING OUTER SOL GRID'

c...      Extract segments from the standard grid:
          ir     = irwall-1
          nseg   = 0
          ikmark = 0

          IF (.TRUE.) THEN

            DO i1 = nks(ir), nks(ir)/2, -1
              id = korpg(i1,ir)
              IF (zvertp(3,id).LT.SNGL(zhor(1))) THEN
c...            Use all the segments that are below the outer midplane:
                nseg = nseg + 1
                rseg(nseg) = rvertp(3,id)
                zseg(nseg) = zvertp(3,id)
                ikmark = i1
              ENDIF
            ENDDO

c...        Add segment to bridge the standard grid and the neutral
c           wall:

            a1 = DBLE(rxp)
            a2 = zhor(1)
            b1 = 100.0D0
            b2 = zhor(1)

            id = korpg(ikmark,ir) 
            c1 = rvertp(3,id)
            c2 = zvertp(3,id)
            d1 = rvertp(2,id)
            d2 = zvertp(2,id)

            CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)

            rseg(nseg+1) = c1 + tcd * (d1 - c1)
            zseg(nseg+1) = c2 + tcd * (d2 - c2)

            rwal(1) = rseg(nseg+1)
            zwal(1) = zseg(nseg+1)

            tmin = HI
            DO i1 = 1, nvesm
              c1 = rvesm(i1,1)
              c2 = zvesm(i1,1)
              d1 = rvesm(i1,2)
              d2 = zvesm(i1,2)

              CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)

              IF (tcd.GT.0.0D0.AND.tcd.LT.1.0D0.AND.
     .            tab.GT.0.0D0.AND.tab.LT.tmin) THEN
                tmin    = tab
                inmark  = i1
                rwal(2) = c1 + tcd * (d1 - c1)
                zwal(2) = c2 + tcd * (d2 - c2)
              ENDIF

            ENDDO
          ELSE
c...JET 
          ENDIF

c...      Extract relevant wall segments:
          mode = .TRUE.
          i1   = inmark-1
          nwal = 2
          DO WHILE (mode.OR.nwal.EQ.0)
            i1 = i1 + 1
            IF (i1.GT.nvesm+nvesp) i1 = 1
	  
            IF (rvesm(i1,2).EQ.rseg(1).AND.
     .          zvesm(i1,2).EQ.zseg(1)) THEN
              mode = .FALSE.
              rwal(nwal+1) = rvesm(i1,2)
              zwal(nwal+1) = zvesm(i1,2)
            ELSE
              nwal = nwal + 1
              rwal(nwal) = rvesm(i1,2)
              zwal(nwal) = zvesm(i1,2)
            ENDIF
          ENDDO

c...      An additional surface must be added to the EIRENE
c         simulation at the main chamber boundary of this 
c         vacuum grid region, so that the additional cell 
c         index search can be triggered in EIRENE.  This is the
c         most flexible solution (at the moment anyway) that
c         that is able to account for the presence, or lack there
c         of, of other vacuum grid regions:

c...TAG-A 
c          eirnasdat = eirnasdat + 1
c          eirasdat(eirnasdat,1) = 2.0
c          eirasdat(eirnasdat,2) = 1.0
c          eirasdat(eirnasdat,3) = rwal(1)
c          eirasdat(eirnasdat,4) = zwal(1)
c          eirasdat(eirnasdat,5) = 1.0
c          eirasdat(eirnasdat,6) = rwal(2)
c          eirasdat(eirnasdat,7) = zwal(2)

c          CALL AssignNIMBUSWall


c...      Shift zhor(1) up slightly so that it is outside the 
c         vacuum cell region:
          zhor(1) = zhor(1) + 0.00001D0


          IF (output) THEN
            DO i1 = 1, nver
              WRITE(SLOUT,90) 'VER-PLE: ',i1,rver(i1)
            ENDDO           
            DO i1 = 1, nhor
              WRITE(SLOUT,90) 'HOR-PLE: ',i1,zhor(i1)
            ENDDO           

            DO i1 = 1, nseg+1
              WRITE(SLOUT,90) 'SEG-PLE: ',i1,rseg(i1),zseg(i1)
            ENDDO
            DO i1 = 1, nwal+1
              WRITE(SLOUT,90) 'WAL-PLE: ',i1,rwal(i1),zwal(i1)
            ENDDO
          ENDIF

          asc_nregion             = asc_nregion + 1
          asc_rstart(asc_nregion) = asc_ncell + 1



          CALL AssembleGrid(nver,rver,nhor,zhor,code,3)




          asc_rend(asc_nregion) = asc_ncell 
        ENDIF




c
c
c
c
c
c
c
c
c
c...    Lame hard-coded vacuum grid option:
        nhor = 0
        nver = 0
        shor = 0.0
        ihor = 0.0
        ehor = 0.0
        sver = 0.0
        iver = 0.0
        ever = 0.0
        DO i1 = 1, vacnseg
          IF (vacseg(i1,1).EQ.99.0) THEN
            IF     (vacseg(i1,4).EQ.0.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                nver       = nver + 1
                rver(nver) = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                nhor       = nhor + 1
                zhor(nhor) = vacseg(i1,3)
              ENDIF
            ELSEIF (vacseg(i1,4).EQ.1.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                sver = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                shor = vacseg(i1,3)
              ENDIF
            ELSEIF (vacseg(i1,4).EQ.2.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                iver = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                ihor = vacseg(i1,3)
              ENDIF
            ELSEIF (vacseg(i1,4).EQ.3.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                ever = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                ehor = vacseg(i1,3)
              ENDIF
            ENDIF
          ENDIF
        ENDDO

        IF (sver.NE.0.0.AND.iver.NE.0.0.AND.ever.NE.0.0) THEN
c...      Build list of vertical segments:       
          rval = sver
          DO WHILE (rval.LT.ever) 
            nver       = nver + 1
            rver(nver) = rval
            rval       = rval + iver
          ENDDO
          nver       = nver + 1
          rver(nver) = ever

c...      Sort the segments from smallest R-value to largest:
          status = .TRUE.
          DO WHILE (status)
            status = .FALSE.
            DO i1 = 1, nver-1
              IF (rver(i1).GT.rver(i1+1)) THEN
                rtmp       = rver(i1)
                rver(i1)   = rver(i1+1)
                rver(i1+1) = rtmp
                status     = .TRUE.
              ENDIF
            ENDDO
          ENDDO
        ENDIF

        IF (shor.NE.0.0.AND.ihor.NE.0.0.AND.shor.NE.0.0) THEN
c...      Build list of horizontal segments:
          zval = shor
          DO WHILE (zval.GT.ehor) 
            nhor       = nhor + 1
            zhor(nhor) = zval
            zval       = zval + ihor
          ENDDO
          nhor       = nhor + 1
          zhor(nhor) = ehor

c...      Sort the segments from smallest R-value to largest:
          status = .TRUE.
          DO WHILE (status)
            status = .FALSE.
            DO i1 = 1, nhor-1
              IF (zhor(i1).LT.zhor(i1+1)) THEN
                ztmp       = zhor(i1)
                zhor(i1)   = zhor(i1+1)
                zhor(i1+1) = ztmp
                status     = .TRUE.
              ENDIF
            ENDDO
          ENDDO
        ENDIF

        IF (nver.GT.0.AND.nhor.GT.0) THEN    
c...      Plenum region vacuum grid lines exist in the input file,
c         so continue:

          code = 1

          vacregion(8) = 1

          IF (output) WRITE(0,*) 'LAME CUSTOM GRID'

c...      Extract relevant wall segments.  Assume all the additional
c         surfaces are to be included:
          nwal = 0
          DO i1 = nvesm+1, nvesm+nvesp
            nwal = nwal + 1
            rwal(nwal) = rvesm(i1,1)
            zwal(nwal) = zvesm(i1,1)
          ENDDO
          rwal(nwal+1) = rvesm(nvesm+nvesp,2)
          zwal(nwal+1) = zvesm(nvesm+nvesp,2)

c          IF (output) WRITE(0,*) 'WALL ASSIGNED'

c...      Extract segments from the standard grid (outer 
c         target segments):
          mode = .FALSE.
          nseg = 0
          i1   = irsep-1
          DO WHILE (mode.OR.nseg.EQ.0)
            i1 = i1 + 1
            IF (i1.EQ.irwall) 
     .        CALL ER('BuildVacuumGrid','Problem with lame grid',*99)
            id = korpg(nks(i1),i1)
            IF (rvertp(4,id).EQ.rwal(nwal+1).AND.
     .          zvertp(4,id).EQ.zwal(nwal+1)) THEN
              mode = .TRUE.
            ENDIF
	    IF (mode) THEN
              nseg = nseg + 1
              rseg(nseg) = rvertp(4,id)
              zseg(nseg) = zvertp(4,id)
            ENDIF
            IF (rvertp(3,id).EQ.rwal(1).AND.
     .          zvertp(3,id).EQ.zwal(1)) THEN
              mode = .FALSE.
              rseg(nseg+1) = rvertp(3,id)
              zseg(nseg+1) = zvertp(3,id)
            ENDIF
          ENDDO

c          IF (output) WRITE(0,*) 'SEGMENTS ASSIGNED'

c          nseg = 0
c          DO ir = irsep, irwall-1
c            nseg = nseg + 1
c            id = korpg(nks(ir),ir)
c            rseg(1) = rvertp(4,id)
c            zseg(1) = zvertp(4,id)
c            rseg(2) = rvertp(3,id)
c            zseg(2) = zvertp(3,id)
c          ENDDO



          IF (output) THEN
            DO i1 = 1, nver
              WRITE(SLOUT,90) 'VER-PLE: ',i1,rver(i1)
            ENDDO           
            DO i1 = 1, nhor
              WRITE(SLOUT,90) 'HOR-PLE: ',i1,zhor(i1)
            ENDDO           

            DO i1 = 1, nseg+1
              WRITE(SLOUT,90) 'SEG-PLE: ',i1,rseg(i1),zseg(i1)
            ENDDO
            DO i1 = 1, nwal+1
              WRITE(SLOUT,90) 'WAL-PLE: ',i1,rwal(i1),zwal(i1)
            ENDDO
          ENDIF

c          IF (output) WRITE(0,*) 'BUILDING VACUUM GRID'

          asc_nregion             = asc_nregion + 1
          asc_rstart(asc_nregion) = asc_ncell + 1

          CALL AssembleGrid(nver,rver,nhor,zhor,code,99)

           


          asc_rend(asc_nregion) = asc_ncell 


        ENDIF






c
c
c
c
c
c
c
c...    Lame hard-coded vacuum grid option #3 - the box:
        nhor = 0
        nver = 0
        shor = 0.0
        ihor = 0.0
        ehor = 0.0
        sver = 0.0
        iver = 0.0
        ever = 0.0
        DO i1 = 1, vacnseg
          IF (vacseg(i1,1).EQ.96.0) THEN
            IF     (vacseg(i1,4).EQ.0.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                nver       = nver + 1
                rver(nver) = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                nhor       = nhor + 1
                zhor(nhor) = vacseg(i1,3)
              ENDIF
            ELSEIF (vacseg(i1,4).EQ.1.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                sver = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                shor = vacseg(i1,3)
              ENDIF
            ELSEIF (vacseg(i1,4).EQ.2.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                iver = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                ihor = vacseg(i1,3)
              ENDIF
            ELSEIF (vacseg(i1,4).EQ.3.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                ever = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                ehor = vacseg(i1,3)
              ENDIF
            ENDIF
          ENDIF
        ENDDO

        IF (sver.NE.0.0.AND.iver.NE.0.0.AND.ever.NE.0.0) THEN
c...      Build list of vertical segments:       
          rval = sver
          DO WHILE (rval.LT.ever) 
            nver       = nver + 1
            rver(nver) = rval
            rval       = rval + iver
          ENDDO
          nver       = nver + 1
          rver(nver) = ever

c...      Sort the segments from smallest R-value to largest:
          status = .TRUE.
          DO WHILE (status)
            status = .FALSE.
            DO i1 = 1, nver-1
              IF (rver(i1).GT.rver(i1+1)) THEN
                rtmp       = rver(i1)
                rver(i1)   = rver(i1+1)
                rver(i1+1) = rtmp
                status     = .TRUE.
              ENDIF
            ENDDO
          ENDDO
        ENDIF

        IF (shor.NE.0.0.AND.ihor.NE.0.0.AND.shor.NE.0.0) THEN
c...      Build list of horizontal segments:
          zval = shor
          DO WHILE (zval.GT.ehor) 
            nhor       = nhor + 1
            zhor(nhor) = zval
            zval       = zval + ihor
          ENDDO
          nhor       = nhor + 1
          zhor(nhor) = ehor

c...      Sort the segments from smallest R-value to largest:
          status = .TRUE.
          DO WHILE (status)
            status = .FALSE.
            DO i1 = 1, nhor-1
              IF (zhor(i1).LT.zhor(i1+1)) THEN
                ztmp       = zhor(i1)
                zhor(i1)   = zhor(i1+1)
                zhor(i1+1) = ztmp
                status     = .TRUE.
              ENDIF
            ENDDO
          ENDDO
        ENDIF

        IF (nver.GT.0.AND.nhor.GT.0) THEN    
c...      In-the-gap region exists so continue:

          code = 1

          vacregion(10) = 1

          IF (output) WRITE(0,*) 'LAME CUSTOM GRID #3'

c...      Extract relevant wall segments.  Assume all the additional
c         surfaces are to be included:
          seg1 = 41
          nwal = 0
          DO i1 = 1, 3
            nwal = nwal + 1
            rwal(nwal) = rvesm(seg1,1)
            zwal(nwal) = zvesm(seg1,1)
            seg1 = seg1 + 1
          ENDDO
          rwal(nwal+1) = rvesm(seg1-1,2)
          zwal(nwal+1) = zvesm(seg1-1,2)           

c...      Extract segments from the standard grid (outer 
c         target segments):
          nseg = 0
          seg1 = 44
          DO i1 = 1, 3
            nseg = nseg + 1
            rseg(nseg) = rvesm(seg1,1)
            zseg(nseg) = zvesm(seg1,1)
            seg1 = seg1 + 1
          ENDDO
          rseg(nseg+1) = rvesm(seg1-1,2)
          zseg(nseg+1) = zvesm(seg1-1,2)

          IF (output) THEN
            DO i1 = 1, nver
              WRITE(SLOUT,90) 'VER-PLE 96: ',i1,rver(i1)
            ENDDO           
            DO i1 = 1, nhor
              WRITE(SLOUT,90) 'HOR-PLE 96: ',i1,zhor(i1)
            ENDDO           

            DO i1 = 1, nseg+1
              WRITE(SLOUT,90) 'SEG-PLE 96: ',i1,rseg(i1),zseg(i1)
            ENDDO
            DO i1 = 1, nwal+1
              WRITE(SLOUT,90) 'WAL-PLE 96: ',i1,rwal(i1),zwal(i1)
            ENDDO
          ENDIF

          asc_nregion             = asc_nregion + 1
          asc_rstart(asc_nregion) = asc_ncell + 1
          CALL AssembleGrid(nver,rver,nhor,zhor,code,96)
          asc_rend(asc_nregion) = asc_ncell 
        ENDIF


c
c...    Lame hard-coded vacuum grid option #2:
        nhor = 0
        nver = 0
        shor = 0.0
        ihor = 0.0
        ehor = 0.0
        sver = 0.0
        iver = 0.0
        ever = 0.0
        DO i1 = 1, vacnseg
          IF (vacseg(i1,1).EQ.97.0) THEN
            IF     (vacseg(i1,4).EQ.0.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                nver       = nver + 1
                rver(nver) = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                nhor       = nhor + 1
                zhor(nhor) = vacseg(i1,3)
              ENDIF
            ELSEIF (vacseg(i1,4).EQ.1.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                sver = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                shor = vacseg(i1,3)
              ENDIF
            ELSEIF (vacseg(i1,4).EQ.2.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                iver = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                ihor = vacseg(i1,3)
              ENDIF
            ELSEIF (vacseg(i1,4).EQ.3.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                ever = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                ehor = vacseg(i1,3)
              ENDIF
            ENDIF
          ENDIF
        ENDDO

        IF (sver.NE.0.0.AND.iver.NE.0.0.AND.ever.NE.0.0) THEN
c...      Build list of vertical segments:       
          rval = sver
          DO WHILE (rval.LT.ever) 
            nver       = nver + 1
            rver(nver) = rval
            rval       = rval + iver
          ENDDO
          nver       = nver + 1
          rver(nver) = ever

c...      Sort the segments from smallest R-value to largest:
          status = .TRUE.
          DO WHILE (status)
            status = .FALSE.
            DO i1 = 1, nver-1
              IF (rver(i1).GT.rver(i1+1)) THEN
                rtmp       = rver(i1)
                rver(i1)   = rver(i1+1)
                rver(i1+1) = rtmp
                status     = .TRUE.
              ENDIF
            ENDDO
          ENDDO
        ENDIF

        IF (shor.NE.0.0.AND.ihor.NE.0.0.AND.shor.NE.0.0) THEN
c...      Build list of horizontal segments:
          zval = shor
          DO WHILE (zval.GT.ehor) 
            nhor       = nhor + 1
            zhor(nhor) = zval
            zval       = zval + ihor
          ENDDO
          nhor       = nhor + 1
          zhor(nhor) = ehor

c...      Sort the segments from smallest R-value to largest:
          status = .TRUE.
          DO WHILE (status)
            status = .FALSE.
            DO i1 = 1, nhor-1
              IF (zhor(i1).LT.zhor(i1+1)) THEN
                ztmp       = zhor(i1)
                zhor(i1)   = zhor(i1+1)
                zhor(i1+1) = ztmp
                status     = .TRUE.
              ENDIF
            ENDDO
          ENDDO
        ENDIF

        IF (nver.GT.0.AND.nhor.GT.0) THEN    
c...      In-the-gap region exists so continue:

          code = 1

          vacregion(9) = 1

          IF (output) WRITE(0,*) 'LAME CUSTOM GRID #2'

c...      Extract relevant wall segments.  Assume all the additional
c         surfaces are to be included:
          seg1 = 34
          nwal = 0
          DO i1 = 1, 7
            nwal = nwal + 1
            rwal(nwal) = rvesm(seg1,2)
            zwal(nwal) = zvesm(seg1,2)
            seg1 = seg1 - 1
            IF (seg1.EQ.31) seg1 = 57
          ENDDO
          rwal(nwal+1) = rvesm(seg1+1,1)
          zwal(nwal+1) = zvesm(seg1+1,1)           

c...      Extract segments from the standard grid (outer 
c         target segments):
          nseg = 0
          i1   = irbreak
          DO WHILE (i1.NE.irwall)
            id = korpg(nks(i1),i1)
            nseg = nseg + 1
            rseg(nseg) = rvertp(4,id)
            zseg(nseg) = zvertp(4,id)
            i1 = i1 + 1
            IF (i1.EQ.nrs+1) i1 = irsep
          ENDDO
          id = korpg(nks(i1-1),i1-1)
          rseg(nseg+1) = rvertp(3,id)
          zseg(nseg+1) = zvertp(3,id)

          IF (output) THEN
            DO i1 = 1, nver
              WRITE(SLOUT,90) 'VER-PLE 97: ',i1,rver(i1)
            ENDDO           
            DO i1 = 1, nhor
              WRITE(SLOUT,90) 'HOR-PLE 97: ',i1,zhor(i1)
            ENDDO           

            DO i1 = 1, nseg+1
              WRITE(SLOUT,90) 'SEG-PLE 97: ',i1,rseg(i1),zseg(i1)
            ENDDO
            DO i1 = 1, nwal+1
              WRITE(SLOUT,90) 'WAL-PLE 97: ',i1,rwal(i1),zwal(i1)
            ENDDO
          ENDIF

          asc_nregion             = asc_nregion + 1
          asc_rstart(asc_nregion) = asc_ncell + 1
          CALL AssembleGrid(nver,rver,nhor,zhor,code,97)
          asc_rend(asc_nregion) = asc_ncell 
        ENDIF
c
c
c
c
c...    Lame hard-coded vacuum grid option #4:
        nhor = 0
        nver = 0
        shor = 0.0
        ihor = 0.0
        ehor = 0.0
        sver = 0.0
        iver = 0.0
        ever = 0.0
        DO i1 = 1, vacnseg
          IF (vacseg(i1,1).EQ.95.0) THEN
            IF     (vacseg(i1,4).EQ.0.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                nver       = nver + 1
                rver(nver) = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                nhor       = nhor + 1
                zhor(nhor) = vacseg(i1,3)
              ENDIF
            ELSEIF (vacseg(i1,4).EQ.1.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                sver = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                shor = vacseg(i1,3)
              ENDIF
            ELSEIF (vacseg(i1,4).EQ.2.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                iver = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                ihor = vacseg(i1,3)
              ENDIF
            ELSEIF (vacseg(i1,4).EQ.3.0) THEN
              IF     (vacseg(i1,2).EQ.1.0) THEN
                ever = vacseg(i1,3)
              ELSEIF (vacseg(i1,2).EQ.2.0) THEN
                ehor = vacseg(i1,3)
              ENDIF
            ENDIF
          ENDIF
        ENDDO

        IF (sver.NE.0.0.AND.iver.NE.0.0.AND.ever.NE.0.0) THEN
c...      Build list of vertical segments:       
          rval = sver
          DO WHILE (rval.LT.ever) 
            nver       = nver + 1
            rver(nver) = rval
            rval       = rval + iver
          ENDDO
          nver       = nver + 1
          rver(nver) = ever

c...      Sort the segments from smallest R-value to largest:
          status = .TRUE.
          DO WHILE (status)
            status = .FALSE.
            DO i1 = 1, nver-1
              IF (rver(i1).GT.rver(i1+1)) THEN
                rtmp       = rver(i1)
                rver(i1)   = rver(i1+1)
                rver(i1+1) = rtmp
                status     = .TRUE.
              ENDIF
            ENDDO
          ENDDO
        ENDIF

        IF (shor.NE.0.0.AND.ihor.NE.0.0.AND.shor.NE.0.0) THEN
c...      Build list of horizontal segments:
          zval = shor
          DO WHILE (zval.GT.ehor) 
            nhor       = nhor + 1
            zhor(nhor) = zval
            zval       = zval + ihor
          ENDDO
          nhor       = nhor + 1
          zhor(nhor) = ehor

c...      Sort the segments from smallest R-value to largest:
          status = .TRUE.
          DO WHILE (status)
            status = .FALSE.
            DO i1 = 1, nhor-1
              IF (zhor(i1).LT.zhor(i1+1)) THEN
                ztmp       = zhor(i1)
                zhor(i1)   = zhor(i1+1)
                zhor(i1+1) = ztmp
                status     = .TRUE.
              ENDIF
            ENDDO
          ENDDO
        ENDIF

        IF (nver.GT.0.AND.nhor.GT.0) THEN    
c...      In-the-gap region exists so continue:

          code = 1

          vacregion(9) = 1

          IF (output) WRITE(0,*) 'LAME CUSTOM GRID #4'

c...      Extract relevant wall segments.  Assume all the additional
c         surfaces are to be included:
          seg1 = 30
          nwal = 0

c          CALL DumpGrid('CUSTOM CRAP')

          IF (nvesm.EQ.99) THEN
            DO i1 = 1, 6
              nwal = nwal + 1
              rwal(nwal) = rvesm(seg1,2)
              zwal(nwal) = zvesm(seg1,2)
              seg1 = seg1 - 1
              IF (seg1.EQ.28) seg1 = 53
            ENDDO
          ELSEIF (nvesm.EQ.107) THEN
            DO i1 = 1, 6
              nwal = nwal + 1
              rwal(nwal) = rvesm(seg1,2)
              zwal(nwal) = zvesm(seg1,2)
              seg1 = seg1 - 1
              IF (seg1.EQ.28) seg1 = 57
            ENDDO
          ELSEIF (nvesm.EQ.111) THEN
            DO i1 = 1, 6
              nwal = nwal + 1
              rwal(nwal) = rvesm(seg1,2)
              zwal(nwal) = zvesm(seg1,2)
              seg1 = seg1 - 1
              IF (seg1.EQ.28) seg1 = 53
            ENDDO
          ELSEIF (nvesm.EQ.102) THEN
            DO i1 = 1, 9
              nwal = nwal + 1
              rwal(nwal) = rvesm(seg1,2)
              zwal(nwal) = zvesm(seg1,2)
              seg1 = seg1 - 1
              IF (seg1.EQ.28) seg1 = 56
            ENDDO
          ELSE
            CALL ER('BuildVaccumGrid','Confused',*99)
          ENDIF
          rwal(nwal+1) = rvesm(seg1+1,1)
          zwal(nwal+1) = zvesm(seg1+1,1)           

c...      Extract segments from the standard grid (outer 
c         target segments):
          nseg = 0
          i1   = irbreak
          DO WHILE (i1.NE.irwall)
            id = korpg(nks(i1),i1)
            nseg = nseg + 1
            rseg(nseg) = rvertp(4,id)
            zseg(nseg) = zvertp(4,id)
            i1 = i1 + 1
            IF (i1.EQ.nrs+1) i1 = irsep
          ENDDO
          id = korpg(nks(i1-1),i1-1)
          rseg(nseg+1) = rvertp(3,id)
          zseg(nseg+1) = zvertp(3,id)

          IF (output) THEN
            DO i1 = 1, nver
              WRITE(SLOUT,90) 'VER-PLE 95: ',i1,rver(i1)
            ENDDO           
            DO i1 = 1, nhor
              WRITE(SLOUT,90) 'HOR-PLE 95: ',i1,zhor(i1)
            ENDDO           

            DO i1 = 1, nseg+1
              WRITE(SLOUT,90) 'SEG-PLE 95: ',i1,rseg(i1),zseg(i1)
            ENDDO
            DO i1 = 1, nwal+1
              WRITE(SLOUT,90) 'WAL-PLE 95: ',i1,rwal(i1),zwal(i1)
            ENDDO
          ENDIF

          asc_nregion             = asc_nregion + 1
          asc_rstart(asc_nregion) = asc_ncell + 1
          CALL AssembleGrid(nver,rver,nhor,zhor,code,95)
          asc_rend(asc_nregion) = asc_ncell 
        ENDIF




c...    Search for local refinement of vacuum grid:
        DO i1 = 1, vacnseg

          IF (vacseg(i1,1).EQ.98.0.AND.vacseg(i1,2).EQ.1.0) THEN

            DO i2 = 1, NINT(vacseg(i1,4))

C...   IF A REGION THAT HAS ALREADY BEEN REFINED IS REFINED AGAIN, THEN
C      IT WILL FAIL, SINCE THE REFERENCE DIMENTION FOR THE REFINED 
C      CELL MAY BE INCORRECT, SINCE THE REFERENCE SIDE CAN BE CHOSEN
C      FROM ASC_LINK=-1, BUT THAT IS NO GOOD FOR REFINED CELLS, SINCE ASC_LINK=-1
C      CAN BE FOR A CELL ON A SWTICHING SURFACE, AND NOT A WALL BOUNDARY:

c... SORT OF WORKS NOW!

c              IF (i2.GT.1) STOP 'CANNOT DO THIS YET'

              CALL LocalGridRefinement_Old
     .               (vacseg(i1  ,3),vacseg(i1+1,3),
     .                vacseg(i1+2,3),vacseg(i1+3,3),
     .                NINT(vacseg(i1+1,4)))

            ENDDO

            CALL AssignAddCellLink
         
c...        Calculate vacuum cell volumes:
            DO cell = 1, asc_ncell
              CALL CalcPolygonVolume
     .               (asc_rvp (1,cell),asc_zvp(1,cell),4,
     .                asc_link(1,cell),asc_vol(cell),cell)
            ENDDO

            DO i3 = asc_rstart(asc_nregion), asc_ncell
              WRITE(SLOUT,93) 'ASC?: ',
     .          i3,asc_cell(i3),asc_region(i3),
     .          (asc_link(i2,i3),i2=1,4),
     .          (asc_grid(i2,i3),i2=1,2),
     .          (asc_rvp (i2,i3),
     .           asc_zvp (i2,i3),i2=1,asc_nvp(i3))
              WRITE(SLOUT,94) 'ASC?: ',
     .          (asc_rvp (i2,i3),
     .           asc_zvp (i2,i3),i2=asc_nvp(i3)+1,2*asc_nvp(i3))
            ENDDO

	
          ENDIF
	
        ENDDO








c...    Search for toroidal cuts:
        DO i1 = 1, vacnseg
          IF     (vacseg(i1,1).EQ.99.0.AND.vacseg(i1,2).EQ.3.0) THEN
c...        High definition of toroidal resolution:

            WRITE(0,*) 'WARNING: NO LONGER SUPPORTED!  WONT WORK IN'//   
     .                 ' ALL SITUATIONS BECAUSE NO LONGER EXCLUSIVELY '
            WRITE(0,*) '         EXECUTED FROM THE LAME GRID IF BLOCK'

            IF (vacseg(i1,4).EQ.-1.0) THEN
              xcutmin = vacseg(i1+2,3)
              xcutmax = vacseg(i1+2,4)
              ycutmin = vacseg(i1+3,3)
              ycutmax = vacseg(i1+3,4)
              DO i2 = 1, 1+NINT(vacseg(i1+1,4))
                zcut = vacseg(i1,3) + REAL(i2-1) * vacseg(i1+1,3)
	
c                WRITE(0     ,*) 'COOL CUTTING LAME GRID AT Z=',zcut 
                WRITE(PINOUT,*) 'COOL CUTTING LAME GRID AT Z=',zcut 
	
                CALL Build3DVacuumGrid(xcutmin,xcutmax,
     .                                 ycutmin,ycutmax,
     .                                 zcut,1)
              ENDDO
              CALL AssignAddCellLink
            ELSE
              xcutmin = 99.0
              xcutmax = 99.0
              ycutmin = 99.0
              ycutmax = 99.0
              zcut = vacseg(i1,3)
	
              WRITE(0,*) 'CUTTING LAME GRID AT Z=',zcut
	
              CALL Build3DVacuumGrid(xcutmin,xcutmax,
     .                               ycutmin,ycutmax,
     .                               zcut,1)
            ENDIF
          ELSEIF (vacseg(i1,1).EQ.99.0.AND.vacseg(i1,2).EQ.4.0) THEN
c...        More memory efficient method of toroidal resolution:
	
            xcutmin = vacseg(i1+2,3)
            xcutmax = vacseg(i1+2,4)
            ycutmin = vacseg(i1+3,3)
            ycutmax = vacseg(i1+3,4)
	    
            DO i2 = 1, 1+NINT(vacseg(i1+1,4))
	    
              zcut = vacseg(i1,3) + REAL(i2-1) * vacseg(i1+1,3)
	    
c              WRITE(0     ,*) 'NEW CUTTING LAME GRID AT Z=',zcut 
              WRITE(PINOUT,*) 'NEW CUTTING LAME GRID AT Z=',zcut 	      
	
              CALL Build3DVacuumGrid(xcutmin,xcutmax,
     .                               ycutmin,ycutmax,
     .                               zcut,2)
	    ENDDO
	
c...        Calculate vacuum cell volumes:
            cell2 = 0
            DO cell = 1, asc_ncell*ascncut
              cell2 = cell2 + 1
              IF (cell2.GT.asc_ncell) cell2 = 1
                CALL CalcPolygonVolume
     .                 (asc_rvp (1,cell2),asc_zvp(1,cell2),4,
     .                  asc_link(1,cell2),asc_vol(cell),cell)
            ENDDO
	
          ENDIF
	
	
	
        ENDDO








      ENDIF

c      CALL Build3DVacuumGrid 
c...  Has to be called from within the vacuum grid region if block because
c     calculating the volumes requires that the SEG's and WAL's be defined properly:
c      IF (asc_3dmode.EQ.1) CALL AssignAddCellLink

c      WRITE(0,*) 'DONE MAKING VACUUM GRID'


        DO i1 = 1, asc_ncell
          WRITE(SLOUT,93) 'ASC3: ',
     .      i1,asc_cell(i1),asc_region(i1),
     .      (asc_link(i2,i1),i2=1,4),
     .      (asc_grid(i2,i1),i2=1,2),
     .      (asc_rvp (i2,i1),
     .       asc_zvp (i2,i1),i2=1,asc_nvp(i1))
          WRITE(SLOUT,94) 'ASC3: ',
     .      (asc_rvp (i2,i1),
     .       asc_zvp (i2,i1),i2=asc_nvp(i1)+1,2*asc_nvp(i1))
        ENDDO
        WRITE(0,*) 'REG,ASC_NCELL= ',asc_ncell

      CALL ProcessVacuumGrid('SAVE')




c...  OUTPUT SOME VOLUMES:

c      i2 = 2
c      sum=0.0
c      DO i1 = 1, ascncut
c        sum=sum+ asc_vol(i2+(i1-1)*asc_ncell)
c        WRITE(0,*) 'VBOLCHECK:',i2,i2+(i1-1)*asc_ncell,
c     .                          asc_vol(i2+(i1-1)*asc_ncell)
c      ENDDO
c      WRITE(0,*) 'sum:',sum 
c      STOP 'sdfd'

      RETURN

90    FORMAT(A,I6,2F12.8)

93    FORMAT(A,3I4,1X,4I4,1X,2I4,1X,10(2F12.8,1X:))
94    FORMAT(A,12X,1X,16X,1X, 8X,1X,10(2F12.8,1X:))


99    CALL OutputGrid(87,'Failure in BuildVacuumGrid')
      DO i1 = 1, nver
        WRITE(SLOUT,90) 'VER-PLE: ',i1,rver(i1)
      ENDDO           
      DO i1 = 1, nhor
        WRITE(SLOUT,90) 'HOR-PLE: ',i1,zhor(i1)
      ENDDO           
      DO i1 = 1, nseg+1
        WRITE(SLOUT,90) 'SEG-PLE: ',i1,rseg(i1),zseg(i1)
      ENDDO
      DO i1 = 1, nwal+1
        WRITE(SLOUT,90) 'WAL-PLE: ',i1,rwal(i1),zwal(i1)
      ENDDO
      STOP
      END


      
