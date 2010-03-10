c     -*-Fortran-*-
c
c
c ======================================================================
c
c subroutine: DrawGrid
c
c
      SUBROUTINE DrawGrid(iopt1)
      USE mod_eirene06_parameters ! 04
      USE mod_eirene06
      IMPLICIT none
      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'comgra'
      INCLUDE 'colours'
      INCLUDE 'slcom'
      INCLUDE 'slout'

      INTEGER iopt1

      REAL       TOL
      PARAMETER (TOL=1.0E-06)

      REAL       DTOL
      PARAMETER (DTOL=1.0D-06)

      INTEGER iopt,i1,i2,nline,lastcolour,v1,v2,ik,ir,id,idum1
      LOGICAL drawseparatrix
      INTEGER, ALLOCATABLE :: lines(:,:),lcolour(:)

      drawseparatrix = .FALSE.

      iopt = iopt1                      ! FIX REQUIRED FOR COMPILER HANGUP ON CHANGING SIGN OF IOPT, SEE BELOW...

      IF (iopt.LT.0) THEN
        iopt = -iopt
        drawseparatrix = .TRUE.
      ENDIF

      WRITE(0,*) 'IOPT:',iopt

c      CALL THICK2(10)

      IF (iopt.LT.500) THEN

        ALLOCATE(lines(3*ntri,2))
        ALLOCATE(lcolour(3*ntri))

        nline = 0
        DO i1 = 1, ntri
          DO v1 = 1, 3
            v2 = v1 + 1
            IF (v1.EQ.3) v2 = 1
c *TEMP*
            IF     (iopt.EQ.95) THEN

              IF (tri(i1)%type.EQ.VACUUM_GRID.AND.     ! Walls
     .            tri(i1)%sur(v1).NE.0) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 2
                lines(nline,1)=tri(i1)%ver(v1)
                lines(nline,2)=tri(i1)%ver(v2)
              ENDIF

c              drawseparatrix = .TRUE.

              IF (tri(i1)%type.EQ.MAGNETIC_GRID.AND.   ! Targets
     .            tri(i1)%sur(v1).NE.0.AND.
     .            tri(i1)%sideindex(2,v1).NE.0) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 3
                lines(nline,1)=tri(i1)%ver(v1)
                lines(nline,2)=tri(i1)%ver(v2)
              ENDIF

            ELSEIF (iopt.EQ.96) THEN

c             All triangles:
              nline = nline + 1
              lcolour(nline) = ncols + 21
              lines(nline,1)=tri(i1)%ver(v1)
              lines(nline,2)=tri(i1)%ver(v2)

            ELSEIF (iopt.EQ.97) THEN
              IF (i1.EQ.6509.OR.i1.EQ.6513) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 3
                lines(nline,1)=tri(i1)%ver(v1)
                lines(nline,2)=tri(i1)%ver(v2)
              ENDIF

            ELSEIF (iopt.EQ.98) THEN

              IF (tri(i1)%type.EQ.VACUUM_GRID) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 3
                lines(nline,1)=tri(i1)%ver(v1)
                lines(nline,2)=tri(i1)%ver(v2)
              ENDIF

              IF (tri(i1)%type.EQ.MAGNETIC_GRID) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 3
                lines(nline,1)=tri(i1)%ver(v1)
                lines(nline,2)=tri(i1)%ver(v2)
              ENDIF

            ELSEIF (iopt.EQ.99) THEN
c...          Magnetic grid triangles + wall surfaces:
              IF (tri(i1)%type.EQ.VACUUM_GRID.AND.
c     .            (tri(i1)%sur(v1).EQ.4.OR.
c     .             tri(i1)%sur(v1).GT.10)) THEN
     .             tri(i1)%sur(v1).NE.0) THEN


                nline = nline + 1
                lcolour(nline) = 1
                lines(nline,1)=tri(i1)%ver(v1)
                lines(nline,2)=tri(i1)%ver(v2)
              ENDIF

              IF (tri(i1)%type.EQ.MAGNETIC_GRID.AND.
     .            v1.NE.1) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 1
                lines(nline,1)=tri(i1)%ver(v1)
                lines(nline,2)=tri(i1)%ver(v2)
              ENDIF

            ELSE

              IF (tri(i1)%type.EQ.VACUUM_GRID.AND.
c     .            (tri(i1)%sur(v1).EQ.4.OR.
c     .             tri(i1)%sur(v1).GT.10)) THEN
     .             tri(i1)%sur(v1).NE.0) THEN
                nline = nline + 1
                lcolour(nline) = 1
                lines(nline,1)=tri(i1)%ver(v1)
                lines(nline,2)=tri(i1)%ver(v2)
              ENDIF

              IF (tri(i1)%type.EQ.MAGNETIC_GRID.AND.
     .            tri(i1)%sur(v1).NE.0) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 1
                lines(nline,1)=tri(i1)%ver(v1)
                lines(nline,2)=tri(i1)%ver(v2)
              ENDIF

              IF (tri(i1)%type.EQ.MAGNETIC_GRID.AND.
     .            tri(i1)%index(2).EQ.irsep.AND.
     .            tri(i1)%sideindex(1,v1).EQ.14) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 3
                lines(nline,1)=tri(i1)%ver(v1)
                lines(nline,2)=tri(i1)%ver(v2)
              ENDIF

           ENDIF

c            nline = nline + 1
c            lines(nline,1)=tri(i1)%ver(v1)
c            lines(nline,2)=tri(i1)%ver(v2)

          ENDDO
        ENDDO
c...    Remove duplicates:
        DO i1 = 1, nline-1
          DO i2 = i1+1, nline
            IF (lines(i1,1).NE.-999.0.AND.lines(i2,1).NE.-999.0) THEN
c * IMPROVE THIS CHECK! *
              IF ((DABS(ver(lines(i1,1),1)-
     .                  ver(lines(i2,2),1)).LT.DTOL.AND.
     .             DABS(ver(lines(i1,1),2)-
     .                  ver(lines(i2,2),2)).LT.DTOL.AND.
     .             DABS(ver(lines(i1,2),1)-
     .                  ver(lines(i2,1),1)).LT.DTOL.AND.
     .             DABS(ver(lines(i1,2),2)-
     .                  ver(lines(i2,1),2)).LT.DTOL).OR.
     .            (DABS(ver(lines(i1,1),1)-
     .                  ver(lines(i2,1),1)).LT.DTOL.AND.
     .             DABS(ver(lines(i1,1),2)-
     .                  ver(lines(i2,1),2)).LT.DTOL.AND.
     .             DABS(ver(lines(i1,2),1)-
     .                  ver(lines(i2,2),1)).LT.DTOL.AND.
     .             DABS(ver(lines(i1,2),2)-
     .                  ver(lines(i2,2),2)).LT.DTOL)) THEN
                IF (lcolour(i1).EQ.1) THEN
                  lines(i1,1) = -999.0
                ELSE
                  lines(i2,1) = -999.0
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        DO i1 = nline, 1, -1
          IF (lines(i1,1).EQ.-999.0) THEN
c            WRITE(0,*) 'DELETING:',i1
            DO i2 = i1, nline-1
              lines(i2,1) = lines(i2+1,1)
              lines(i2,2) = lines(i2+1,2)
              lcolour(i2) = lcolour(i2+1)
            ENDDO
            nline = nline - 1
          ENDIF
        ENDDO
      ENDIF




      IF (.TRUE.) THEN
c...    Plot polygons:
        CALL PSPACE (map1x,map2x,map1y,map2y)      
        CALL MAP    (cxmin,cxmax,cymin,cymax)
        lastcolour = -1
        DO i1 = 1, nline
          IF (lastcolour.NE.lcolour(i1)) THEN
            CALL LINCOL(lcolour(i1)) 
            lastcolour = lcolour(i1)
          ENDIF
          CALL POSITN(SNGL(ver(lines(i1,1),1)),SNGL(ver(lines(i1,1),2)))
          CALL JOIN  (SNGL(ver(lines(i1,2),1)),SNGL(ver(lines(i1,2),2)))
        ENDDO
        DEALLOCATE(lines)
        DEALLOCATE(lcolour)
      ENDIF



      IF (drawseparatrix) THEN

        CALL THICK2(1)
        CALL THICK (1)
        CALL BROKEN(5,5,5,5)

c        CALL LINCOL(ncols+2) 
        CALL LINCOL(ncols+4) 
  
        ir = irsep
        DO ik = 1, nks(ir)
          id = korpg(ik,ir)
          CALL POSITN(rvertp(1,id),zvertp(1,id))
          CALL JOIN  (rvertp(4,id),zvertp(4,id))
        ENDDO

        IF (nrs.EQ.65) THEN
          ir = 38
          DO ik = 1, nks(ir)
            id = korpg(ik,ir)
            CALL POSITN(rvertp(2,id),zvertp(2,id))
            CALL JOIN  (rvertp(3,id),zvertp(3,id))
          ENDDO
        ENDIF

        IF (nrs.EQ.72) THEN
          ir = 47
          DO ik = 1, nks(ir)
            id = korpg(ik,ir)
            CALL POSITN(rvertp(2,id),zvertp(2,id))
            CALL JOIN  (rvertp(3,id),zvertp(3,id))
          ENDDO
        ENDIF

c        IF (irwall.GT.23) THEN
c          ir = 31
c          DO ik = 1, nks(ir)
c            id = korpg(ik,ir)
c            CALL POSITN(rvertp(1,id),zvertp(1,id))
c            CALL JOIN  (rvertp(4,id),zvertp(4,id))
c          ENDDO
c          ir = 40
c          DO ik = 1, nks(ir)
c            id = korpg(ik,ir)
c            CALL POSITN(rvertp(1,id),zvertp(1,id))
c            CALL JOIN  (rvertp(4,id),zvertp(4,id))
c          ENDDO
c        ENDIF

        CALL FULL
      ENDIF




c...  Frame:
      CALL DrawFrame


      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: DrawColourScale
c
c
      SUBROUTINE DrawColourScale(mode,colmode,qmin,qmax,label)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slout'

      INTEGER   mode,colmode 
      REAL      qmin,qmax
      CHARACTER label*(*)

      INTEGER CH1

      REAL      qval,dspot,spot,dist,minx,maxx,miny,maxy,stepsize,i1r,
     .          rscale,dscale,xmove,ymove,xpos,ypos,nchar
      LOGICAL   yshift
      CHARACTER nums*256      
 
      CALL PSPACE (0.0, 1.35, 0.0,1.0)
      CALL CSPACE (0.0, 1.35, 0.0,1.0)
      CALL MAP    (0.0, 1.35, 0.0,1.0)
      CALL FULL
      CALL HSV  ! This should be moved...

      IF     (mode.EQ.0) THEN  ! Don't draw a scale
      ELSEIF (mode.EQ.1) THEN  ! Vertical scale drawn to the left of the plot
        dspot = 0.016         
        minx = map2x + 0.02
        maxx = map2x + 0.04
        miny =  HI
        maxy = map2y - dspot/2.0
        dscale = 2.0
        DO rscale = 100.0-dscale, 0.0, -dscale
          qval = (rscale + 0.5 * dscale) / 100.0 * (qmax - qmin) + qmin
          CALL SetCol255_04(colmode,qval,qmin,qmax)
          CALL FILCOL(255)
          CALL LINCOL(255)
          SPOT  = maxy - (100.0 - dscale - rscale) / 4.0 * DSPOT
          CALL BOX (minx,maxx,SPOT-DSPOT/2.0,SPOT+DSPOT/2.0)
          miny = MIN(miny,SPOT-DSPOT/2.0)
          maxy = MAX(maxy,SPOT+DSPOT/2.0)
        ENDDO
c...    Box:
        CALL LINCOL(1)
        CALL POSITN (minx,miny)
        CALL JOIN   (minx,maxy)
        CALL POSITN (minx,maxy)
        CALL JOIN   (maxx,maxy)
        CALL POSITN (maxx,maxy)
        CALL JOIN   (maxx,miny)
        CALL POSITN (maxx,miny)
        CALL JOIN   (minx,miny)
c...    Text:
        CALL CTRMAG(12)
        dscale = 20.0
        DO rscale = 100.0, 0.0, -dscale
          qval = qmin + (qmax - qmin) * rscale / 100.0
          IF     (qval.GT.-1.0.AND.qval.LT.10.0) THEN
            WRITE(nums,'(F3.1)') qval
          ELSEIF (ABS(qval).LT.100.0) THEN
            WRITE(nums,'(I3)') NINT(qval)
          ELSEIF (ABS(qval).LT.1000.0) THEN
            WRITE(nums,'(I4)') NINT(qval)
          ELSE
            WRITE(nums,'(1P,E9.1)') qval
          ENDIF
c          IF (qmax.GT.1.0.AND.qmax.LT.100.0) THEN
c            WRITE(nums,'(F4.1)') 
c     .        qval
cc            WRITE(nums,'(F4.1,A,I3,A)') 
cc     .        qval,' (',NINT(qval/qmax*100.0),'%)'
c          ELSE
c            WRITE(nums,'(1P,E10.2,0P,A,I3,A)') 
c     .        qval,' (',NINT(qval/qmax*100.0),'%)'
c          ENDIF
          spot = rscale / 100.0 * (maxy - miny) + miny
          IF (label.NE.'none') THEN
c...        Tick:
            CALL POSITN(maxx      ,spot)
            CALL JOIN  (maxx+0.005,spot)
c...        Label:
            CALL PLOTST(maxx+0.015,spot,nums(1:LEN_TRIM(nums)))
c            CALL PLOTST(maxx+0.005,spot,nums(1:LEN_TRIM(nums)))
          ENDIF
        ENDDO
        IF (label.NE.'none') THEN
          CALL CTRORI(90.0)
          CALL PLOTST(maxx+0.070,0.5*(miny+maxy),
     .                label(1:LEN_TRIM(label)))
          CALL CTRORI(0.0)
        ENDIF
c        CALL PLOTST(maxx+0.010,0.5*(miny+maxy),'T')


      ELSEIF (mode.EQ.2) THEN  ! Horizontal scale below plot
        dspot = 0.010 ! for dist = 0.70 987 plots                
c        dspot = 0.012 ! for dist = 0.80 987 plots                
c        dspot = 0.0134 
c        dspot = 0.016
        IF (label.EQ.'none') THEN
          minx = HI                ! map2x + 0.02
          maxx = map2x - dspot/2.0 ! map2x + 0.04
          miny = map1y - 0.030     ! map1y - 0.04      ! HI
          maxy = map1y - 0.010     ! map1y - 0.02      ! map2y - dspot/2.0
        ELSE
          minx = HI                ! map2x + 0.02
          maxx = map2x - dspot/2.0 ! map2x + 0.04
          miny = map1y - 0.075     ! map1y - 0.04      ! HI
          maxy = map1y - 0.055     ! map1y - 0.02      ! map2y - dspot/2.0
        ENDIF
        dscale = 2.0
        DO rscale = 100.0-dscale, 0.0, -dscale
          qval = (rscale + 0.5 * dscale) / 100.0 * (qmax - qmin) + qmin
          IF (qmin.NE.qmax) CALL SetCol255_04(colmode,qval,qmin,qmax)
          CALL FILCOL(255)
          CALL LINCOL(255)
          SPOT  = maxx - (100.0 - dscale - rscale) / 4.0 * DSPOT
          CALL BOX (spot-dspot/2.0,spot+dspot/2.0,miny,maxy)
          minx = MIN(minx,SPOT-DSPOT/2.0)
          maxx = MAX(maxx,SPOT+DSPOT/2.0)
        ENDDO
c...    Box:
        CALL LINCOL(1)
        CALL POSITN (minx,miny)
        CALL JOIN   (minx,maxy)
        CALL POSITN (minx,maxy)
        CALL JOIN   (maxx,maxy)
        CALL POSITN (maxx,maxy)
        CALL JOIN   (maxx,miny)
        CALL POSITN (maxx,miny)
        CALL JOIN   (minx,miny)
c...    Text:
        CALL CTRMAG(12)
        dscale = 20.0
        yshift = .FALSE.
        ymove = 0.0
        DO rscale = 100.0, 0.0, -dscale
          qval = qmin + (qmax - qmin) * rscale / 100.0
c          qval = qval / qmax
          IF     (qval.GT.-0.999.AND.qval.LT.10.0) THEN
c          IF     (qmax.GT.0.1.AND.qmax.LE.1.0) THEN
            WRITE(nums,'(F3.1)') qval
          ELSEIF (ABS(qval).LT.999.0) THEN
c          ELSEIF ((qmax.GT. 1.0.AND.qmax.LE. 999.9).OR.
c     .            (qmax.LE.-1.0.AND.qmax.GT.-999.9)) THEN
            WRITE(nums,'(I4)') NINT(qval)
          ELSE
            WRITE(nums,'(1P,E9.1)') qval
            yshift = .TRUE.
          ENDIF
          spot = rscale / 100.0 * (maxx - minx) + minx

          IF (label.NE.'none') THEN
c...        Tick:
            CALL POSITN(spot,miny      )
            CALL JOIN  (spot,miny-0.005)
c...        Label:
            ypos = miny - 0.020 + ymove
            CALL PCSCEN(spot,ypos,nums(CH1(nums):LEN_TRIM(nums)))
            IF (yshift.AND.ymove.EQ.0.0) THEN
              ymove = -0.02
            ELSE
              ymove = 0.0
            ENDIF       
          ENDIF
        ENDDO

        IF (label.NE.'none') THEN
          IF (yshift) ymove = -0.02
          IF (.TRUE.) THEN
            nchar = REAL(LEN_TRIM(label) - CH1(label) + 1) 
            xpos = 0.5 * (minx + maxx) - 0.5 * nchar * 0.0090
          ELSE
            xpos = 0.5 * (minx + maxx)  !need to not count formatting characters...
          ENDIF
          ypos = miny - 0.045 + ymove 
          CALL PLOTST(xpos,ypos,label(CH1(label):LEN_TRIM(label)))
        ENDIF

      ELSE
        CALL ER('DrawColourScale','Invalid mode',*99)
      ENDIF

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: SetCol255
c
c
      SUBROUTINE SetCol255_04(mode,qval,qmin,qmax)
      IMPLICIT none

      INTEGER mode 
      REAL    qval,qmin,qmax
    

c * OLD *
      REAL frac1,hue,pastel,midfrac,scale
      COMMON /PLOT982COM/ hardscale,colscale
      LOGICAL             hardscale,colscale
      LOGICAL grayscale

      INTEGER last_mode
      REAL bright
      REAL frac,frac5,fmod5,last_frac,last_bright   ! mode.eq.1

      DATA last_mode,last_frac,last_bright /-1, -1.0, -1.0/    

      SAVE


      IF     (mode.EQ.1) THEN

        frac = (qval - qmin) / (qmax - qmin + 1.0E-10)
        frac5 = 100.0*frac
        fmod5 = AMOD(frac5,2.0)
        frac = MIN(0.98,(frac5-fmod5)/100.0)

        bright = 1.0-(0.98-frac)**20

        IF (mode.NE.last_mode.OR.frac.NE.last_frac.OR.
     .      bright.NE.last_bright) 
     .    CALL ColSet(1.0-0.75*frac,1.0,bright,255)
c     .    CALL ColSet(0.75*frac+0.25,1.0,bright,255)

c      ELSEIF (mode.EQ.2) THEN

      ELSEIF (mode.EQ.2) THEN
c...    Grayscale:

        frac = (qval - qmin) / (qmax - qmin)

c        frac = frac * 0.9 + 0.1

        IF (mode.NE.last_mode.OR.frac.NE.last_frac) 
c     .    CALL ColSet(0.0,0.0,frac,255)  ! Default
c     .    CALL ColSet(0.50,1.0-frac,frac, ! Original DIVCAM
c     .                255)
     .    CALL ColSet(0.30+0.230*frac,
     .                MIN(1.0,4.0*(1.0-frac)    ),
     .                MIN(1.0,4.0*     frac     ),
     .                255)
c     .    CALL ColSet(0.25+1.000*frac,
c     .                MIN(1.0,4.0*(1.0-frac)    ),
c     .                MIN(1.0,4.0*     frac     ),
c     .                255)
c     .    CALL ColSet(0.50,1.0-2.0*(MAX(0.0,frac**0.5-0.5)),frac**0.5,
c     .                255)

c        CALL ColSet(0.0,0.0,frac,255)
c        CALL ColSet(0.0,0.0,1.0-frac,255)


      ELSE
        CALL ER('SetCol255_04','Invalid mode',*99)
      ENDIF


      last_mode = mode
      last_frac = frac
      last_bright = bright

      RETURN
c
c * OLD *
c

      grayscale = .FALSE.

      frac = frac1
      IF (qmin.LT.-1.0E-6.AND.qmax.GT.1.0E-06) THEN
c...
        midfrac = -qmin / (qmax - qmin) 
        scale = MAX(midfrac,1.0-midfrac)

        IF (frac.LT.midfrac) THEN
          hue    = 0.34
          pastel = (midfrac - frac) / scale
          bright = 1.0
        ELSE
          hue    = 0.0
          pastel = (frac - midfrac) / scale
          bright = 1.0
        ENDIF
        CALL ColSet(hue,pastel,bright,255)

      ELSEIF (grayscale) THEN
        IF (hardscale.AND.frac.LT.0.07) THEN
          CALL ColSet(0.0,0.0,1.0,255)
        ELSE
          IF (frac.GT.1.0) STOP 'sdfsdsgsd'
          IF (frac.LT.0.0) STOP 'sdfsdsgsd adsfasd'
          frac = frac * 0.9 + 0.1
          CALL ColSet(0.0,0.0,1.0-frac,255)
        ENDIF

      ELSEIF (hardscale) THEN

        IF (.NOT..TRUE.) THEN
          IF (frac.LT.0.01) THEN
            bright = 1.0
            frac   = 1.0
            pastel = 0.0
          ELSE
            frac = frac * 0.93 + 0.07
            IF (frac.LE.0.27) THEN
              bright = 1.0-((0.27-frac)/(0.27-0.07))**2
              bright = MAX(0.0,bright)
            ELSE
              bright = 1.0
            ENDIF
            frac = (1.0 - frac) * 0.90
            frac = frac + 0.34
            IF (frac.GT.1.0) frac = frac - 1.0
            pastel = 1.0
          ENDIF
        ELSE
          IF (frac.LT.0.07) THEN
            bright = 1.0
            frac   = 1.0
            pastel = 0.0
          ELSE
            IF (frac.LE.0.27) THEN
              bright = 1.0-((0.27-frac)/(0.27-0.07))**2
              bright = MAX(0.0,bright)
            ELSE
              bright = 1.0
            ENDIF
            frac = (1.0 - frac) * 0.90
            frac = frac + 0.34
            IF (frac.GT.1.0) frac = frac - 1.0
            pastel = 1.0
          ENDIF
        ENDIF
        CALL ColSet(frac,pastel,bright,255)

      ELSE
        bright = 1.0-(0.98-frac)**20
        CALL ColSet(frac,1.0,bright,255)

      ENDIF

      RETURN
99    STOP
      END
c
c ======================================================================
c
c
      SUBROUTINE NextLine(fp1,ntally,icount,rdum)
      IMPLICIT none

      INTEGER   fp,fp1,ntally,icount,i1
      REAL      rdum(*)   
      CHARACTER buffer*512
      LOGICAL output

c      output = .FALSE.
c      IF (fp1.EQ.44) THEN
c        fp = 99
c        output = .TRUE.
c      ELSE
        fp = fp1
c      ENDIF

      DO WHILE (.TRUE.) 
        READ(fp,'(A512)',END=98) buffer               
c        IF (output) WRITE(0,*) 'BUFFER:',icount,buffer(1:50)
        IF (buffer(1:1).EQ.'*') CYCLE 
        READ(buffer,*,ERR=97) icount,(rdum(i1),i1=1,ntally)          
c        IF (output) WRITE(0,*) 'BUFFER:',icount,buffer(1:50)
        RETURN
      ENDDO
      
 97   CALL ER('NextLine','Data format error',*99)
 98   CALL ER('NextLine','Unexpected end-of-file',*99)
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE Plot987(job,graph,ref,title,iopt,
     .                   xxmin,xxmax,yymin,yymax,ft,fp,zadj,
     .                   ismoth,ignors,itec,avs,navs,nizs)
      USE mod_eirene06_parameters
      USE mod_eirene06 
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'cedge2d'
      INCLUDE 'pindata'
      INCLUDE 'dynam2'
      INCLUDE 'dynam3'
      INCLUDE 'comgra'
      INCLUDE 'colours'
      INCLUDE 'slcom'
      INCLUDE 'slout'

      REAL       TOL
      PARAMETER (TOL=1.0E-06)

      REAL*8     DTOL
      PARAMETER (DTOL=1.0D-06)

      INTEGER CH1
      REAL    GetMach

      INTEGER   ismoth,IGNORS(MAXNGS),ITEC,NAVS,iopt,nizs
      REAL      XXMIN,XXMAX,YYMIN,YYMAX,ft,fp,zadj,AVS(0:100)
      CHARACTER TITLE*(*),JOB*(*),GRAPH*(*),REF*(*)

      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost

      COMMON /PLOT982COM/ hardscale,colscale
      LOGICAL             hardscale,colscale

      CHARACTER table*36,nview*36,anly*36,plane*36,
     .          YLABEL*256,XLABEL*36,smooth*64,glabel*512

      INTEGER i1,i2,v1,v2,nc,ik,ir,id,iz,
     .        nline,lastcolour,scaleopt,colouropt,posopt
      LOGICAL setqmin,setqmax,inside,scale_set
      REAL    qmin,qmax,frac,frac5,fmod5,scalefact,fact,
     .        posx,posy,poswidth,posheight,rdum(4),taus
      CHARACTER label*512,cdum1*512,cdum2*512

      REAL, POINTER :: gdata(:,:)

      INTEGER, ALLOCATABLE :: lines2(:,:),nv(:),lcolour(:)
      REAL   , ALLOCATABLE :: tdata(:),tdata1(:)
      REAL, TARGET, ALLOCATABLE :: gdata1(:,:)
      REAL   , ALLOCATABLE :: rv(:,:),zv(:,:),cq(:)

c...  ADAS:
      CHARACTER ADASID*80,adasex*3,graph3*80
      integer   adasyr,ISELE,ISELR,ISELX,iseld,ierr,ircode
      REAL      wlngth

      INTEGER numplots
      REAL    dx,dy,dist,save_map2x
      DATA    numplots /0/
 
      LOGICAL :: reset_origin = .TRUE.

      SAVE

c      WRITE(0,*) 'DATA:',job
c      WRITE(0,*) 'DATA:',graph
c      WRITE(0,*) 'DATA:',title
c      WRITE(0,*) 'DATA:',ref

      WRITE(glabel,'(512(A:))') (' ',i1=1,LEN(glabel)) 

      IF (numplots.EQ.0) 
     .  glabel = job  (CH1(job  ):LEN_TRIM(job  ))//'     '//
     .           graph(CH1(graph):LEN_TRIM(graph))

      WRITE(nview ,'(512(A:))') (' ',i1=1,LEN(nview )) 
      WRITE(plane ,'(512(A:))') (' ',i1=1,LEN(plane )) 
      WRITE(anly  ,'(512(A:))') (' ',i1=1,LEN(anly  )) 
      WRITE(table ,'(512(A:))') (' ',i1=1,LEN(table )) 
      WRITE(xlabel,'(512(A:))') (' ',i1=1,LEN(xlabel)) 
      WRITE(ylabel,'(512(A:))') (' ',i1=1,LEN(ylabel)) 
      WRITE(smooth,'(512(A:))') (' ',i1=1,LEN(smooth)) 
c      WRITE(graph ,'(512(A:))') (' ',i1=1,LEN(graph )) 
c      WRITE(job   ,'(512(A:))') (' ',i1=1,LEN(job   )) 


      xlabel = 'z (m)'
      ylabel = 'r (m)'
c      xlabel = 'none'
c      ylabel = 'none'

c...  Use GRTSET_TRIM:
      slopt4 = 1
c...  Stopping resizing of scale font in ghost1.o6a:
      iopt_ghost = 1

c      cxmin=xxmin     
c      cxmax=xxmax 
c      cymin=yymin 
c      cymax=yymax 

c      CALL SLSET (0.05,0.95,0.05,0.95,
c     .            xxmin,xxmax,yymin,yymax)
c      CALL PSPACE (map1x,map2x,map1y,map2y)

      scale_set = .FALSE.
      scaleopt   = 2
      colouropt  = 1
      scalefact  = 1.0
      label = 'default'

      qmin =  HI
      qmax = -HI
c...  Read scale information:
      READ(5,'(A512)') cdum1
      IF   (cdum1(8:12).EQ.'Scale'.OR.cdum1(8:12).EQ.'scale'.OR.
     .      cdum1(8:12).EQ.'SCALE') THEN
        scale_set = .TRUE.
        READ(cdum1,*) cdum2,scaleopt,colouropt,scalefact,qmin,qmax,
     .                label
        IF (qmin.EQ.-99.0) qmin =  HI
        IF (qmax.EQ.-99.0) qmax = -HI
c        WRITE(0,*) 'SCALE:',scaleopt,colouropt,scalefact,
c     .              label(1:LEN_TRIM(label))
        IF (label.EQ.'default') 
     .    label = graph(CH1(graph):LEN_TRIM(graph))
      ELSE
        BACKSPACE 5
      ENDIF

c...  Plot location:      
      READ(5,'(A512)') cdum1
      IF   (cdum1(8:15).EQ.'Position'.OR.cdum1(8:15).EQ.'position'.OR.
     .      cdum1(8:15).EQ.'POSITION') THEN
        READ(cdum1,*) cdum2,posopt,posx,posy,poswidth,posheight

        WRITE(0,*) 'POSITION',posopt,posx,posy,poswidth,posheight

        IF     (posopt.EQ.1) THEN
          map1x = 0.05
          map2x = map1x + 0.80 * (xxmax-xxmin) / (yymax-yymin)
c          map2x = map1x + 0.75 * (xxmax-xxmin) / (yymax-yymin)
          map1y = 0.15
          map2y = map1y + 0.80
c          map2y = map1y + 0.75
c          map1x = 0.05
c          map2x = map1x + 0.40 * (xxmax-xxmin) / (yymax-yymin)
c          map1y = 0.15
c          map2y = map1y + 0.40
        ELSEIF (posopt.EQ.2) THEN
          map1x = posx
          map2x = map1x + poswidth * (xxmax-xxmin) / (yymax-yymin)
          map1y = posy
          map2y = map1y + poswidth
        ENDIF

      ELSE
        BACKSPACE 5
        dist = 0.70 ! 0.80
        dx = dist * (xxmax-xxmin) / (yymax-yymin)
        IF (reset_origin) THEN
          map1x = 0.05 
        ELSE
          map1x = save_map2x
        ENDIF
c        map1x = 0.05 + REAL(numplots) * dx
        map2x = map1x + dx
        map1y = 0.20
        map2y = map1y + dist
        save_map2x = map2x
      ENDIF

      CALL FULL
c      CALL THICK2(6)
      CALL THICK2(1)
      CALL THICK(1)

      IF (numplots.NE.0) ylabel = 'none'

      CALL GRTSET_TRIM (TITLE,' ',' ',' ',glabel,
     >                  xXMIN,xXMAX,yYMIN,yYMAX,
     .                  ' ',xlabel,ylabel,
     .                  0,' ',0,' ',1)

c      CALL GRTSET_TRIM (TITLE,REF,nVIEW,PLANE,glabel,
c     >                  xXMIN,xXMAX,yYMIN,yYMAX,
c     .                  TABLE,XLABEL,YLABEL,
c     .                  0,smooth,0,ANLY,1)


      IF (iopt.LT.500) CALL LoadTriangles_06  ! just LoadTriangles before...

c...problem is some small triangles driving up the neutral density? ... check against magnetic 
c grid neutral density plot? or check agains ionisation plot to see if scales are the same... 
c...check .eirdat file...
c      WRITE(0,*) 'DEBUG: IOPT',iopt

      IF (iopt.GE.500) THEN

        ALLOCATE(gdata1(MAXNKS,MAXNRS))
        gdata1 = 0.0

c        WRITE(0,*) 'DEBUG: IOPT',iopt

        SELECTCASE (iopt)
          CASE (501)
            gdata => e2dion
          CASE (502)
            gdata => e2drec
          CASE (520)  
            gdata => e2dnbs
          CASE (521) 
            gdata => e2dtebs
          CASE (522) 
            gdata => e2dtibs
          CASE (530)  ! Cross-field metric THETAG 
            gdata1 = 0.0
            DO ir = 2, nrs ! irsep-1
              gdata1(1:nks(ir),ir) = thetag(1:nks(ir),ir)
            ENDDO
            gdata => gdata1
          CASE (531)  ! Plot of BRATIO
            gdata1 = 0.0
            DO ir = 2, nrs ! irsep-1
              DO ik = 1, nks(ir)
                gdata1(ik,ir) = 1.0 / bratio(ik,ir) ! ksb(ik,ir) - ksb(ik-1,ir)
              ENDDO
            ENDDO
            gdata => gdata1

          CASE (540) 
            gdata1 = 0.0
            DO ir = 2, nrs ! irsep-1
              gdata1(1:nks(ir),ir) = pinion(1:nks(ir),ir)
            ENDDO
            gdata => gdata1
          CASE (542) 
            gdata1 = 0.0
            DO ir = 2, nrs
              IF (idring(ir).EQ.BOUNDARY) CYCLE
              gdata1(1:nks(ir),ir) = pinalpha(1:nks(ir),ir)
            ENDDO
            gdata => gdata1
          CASE (599) 
            gdata => e2dcxrec  ! Testing
          CASE (600)  ! Atom density
            gdata => pinatom
          CASE (601)  ! Molecule density
            gdata => pinmol
          CASE (700)  ! Mach no. (absolute)
c            DO ir = irsep, nrs
            DO ir = 2, nrs
              IF (idring(ir).EQ.BOUNDARY) CYCLE
              DO ik = 1, nks(ir)
                gdata1(ik,ir) = GetMach
     .                    (kvhs(ik,ir)/qtim,ktebs(ik,ir),ktibs(ik,ir)) * 
     .                    SIGN(1.0,kvhs(ik,ir))
              ENDDO
            ENDDO
            gdata => gdata1
          CASE (701)  ! Velocity (absolute)
c            DO ir = irsep, nrs
            DO ir = 2, nrs
              IF (idring(ir).EQ.BOUNDARY) CYCLE
              gdata1(:,ir) = ABS(kvhs(:,ir)) / qtim
            ENDDO
            gdata => gdata1
          CASE (708) 
            gdata1 = 0.0
c            DO ir = irsep, nrs
            DO ir = 2, nrs
              IF (idring(ir).EQ.BOUNDARY) CYCLE
              gdata1(1:nks(ir),ir) = ktebs(1:nks(ir),ir)
            ENDDO
            gdata => gdata1
          CASE (710) 
            gdata1 = 0.0
c            DO ir = irsep, nrs
            DO ir = 2, nrs
              IF (idring(ir).EQ.BOUNDARY) CYCLE
              gdata1(1:nks(ir),ir) = knbs(1:nks(ir),ir)
            ENDDO
            gdata => gdata1
          CASE (720) 
            gdata1 = 0.0
c            DO ir = irsep, nrs
            DO ir = 2, nrs
              IF (idring(ir).EQ.BOUNDARY) CYCLE
              gdata1(1:nks(ir),ir) = ktibs(1:nks(ir),ir)
            ENDDO
            gdata => gdata1
          CASE (799:875) 
c              WRITE(0,*) 'DEBUG: here 1'
              iz = iopt - 800
              READ(5,'(A512)') cdum1
c             ----------------------------------------------------------
              IF (cdum1(8:11).EQ.'Adas'.OR.cdum1(8:11).EQ.'ADAS'.OR.
     .            cdum1(8:11).EQ.'adas') THEN
c...            Load PLRP data from ADAS:
                BACKSPACE 5
                CALL RDG1 (GRAPH3,ADASID,adasyr,adasex,
     .                     ISELE,ISELR,ISELX,ISELD,IERR)
                WRITE(0,*) 'ADAS SETTINGS:',adasid,adasyr,adasex
                WRITE(0,*) '             :',isele,iselr,iselx,iseld
                IF (IERR.NE.0) THEN
                  WRITE(6,*) '987: ERROR READING ADAS DETAILS, '//
     .                       'IERR = ',IERR
                  IERR = 0
                  GOTO 99
                ENDIF
                WRITE(0,*) '             :',cion,iz
                CALL LDADAS(cion,IZ,ADASID,ADASYR,ADASEX,ISELE,ISELR,
     .                      ISELX,gdata1,Wlngth,IRCODE)
                WRITE(0,*) 'ADAS DATA:',iz,wlngth,ircode
c             ----------------------------------------------------------
              ELSEIF (cdum1(8:10).EQ.'Ion'.OR.cdum1(8:10).EQ.'ION'.OR.
     .                cdum1(8:10).EQ.'ion') THEN
c                WRITE(0,*) 'Loading ionisation data'
                DO ir = 2, nrs
                  IF (idring(ir).EQ.BOUNDARY) CYCLE
                  gdata1(1:nks(ir),ir) = tizs(1:nks(ir),ir,iz)
                ENDDO
c             ----------------------------------------------------------
              ELSEIF (cdum1(8:12).EQ.'Power'.OR.cdum1(8:12).EQ.'POWER'
     .                .OR.cdum1(8:12).EQ.'power') THEN
                WRITE(0,*) 'Loading total radiated power data',iz,
     .                     MIN(cion,nizs)
                IF (iz.EQ.-1) THEN
                  DO iz = 0, MIN(cion,nizs)
                    DO ir = 2, nrs
                      IF (idring(ir).EQ.BOUNDARY) CYCLE
                      gdata1(1:nks(ir),ir) = gdata1(1:nks(ir),ir   ) + 
     .                                       powls (1:nks(ir),ir,iz)
                    ENDDO
                  ENDDO
                ELSE
                  DO ir = 2, nrs
                    IF (idring(ir).EQ.BOUNDARY) CYCLE
                    gdata1(1:nks(ir),ir) = powls(1:nks(ir),ir,iz) 
                  ENDDO
                ENDIF
                gdata1 = gdata1 * absfac
        WRITE(6,*) 'powls ioout :',powls(1,irsep,:)
        WRITE(6,*) 'ddlims ioout:',sdlims(1,irsep,:)
        WRITE(6,*) 'gdata1      :',gdata1(1,irsep)
c             ----------------------------------------------------------
              ELSEIF (cdum1(8:15).EQ.'Legrange'.OR.    ! Net force on impurities, from OUT 
     .                cdum1(8:15).EQ.'LEGRANGE'.OR.    ! plot 669/670
     .                cdum1(8:15).EQ.'legrange') THEN
                FACT = QTIM**2 * EMI / CRMI
                DO ir = 2, nrs ! irsep-1
                  DO ik = 1, nks(ir)
                    TAUS = CRMI * KTIBS(IK,IR)**1.5 * SQRT(1.0/CRMB) /
     +                     (6.8E-14 * (1 + CRMB / CRMI) * KNBS(IK,IR) *
     +                     REAL(IZ)**2.0 * RIZB**2 * 15.0)
                    RDUM(1) = AMU * CRMI * KVHS(IK,IR) / QTIM / TAUS
                    RDUM(2) = KFIGS(IK,IR) * KBETAS(IZ) * ECH / FACT
                    RDUM(3) = KFEGS(IK,IR) * KALPHS(IZ) * ECH / FACT
                    RDUM(4) = REAL(IZ) * KES(IK,IR) * ECH / FACT
                    WRITE(6,'(A,2I6,5E10.2)') 
     .                'FORCES:',ik,ir,rdum(1:4),
     .                ABS(sum(rdum(1:4)))*scalefact
                    gdata1(ik,ir) = ABS(SUM(rdum(1:4)))
                  ENDDO
                ENDDO
                gdata => gdata1
c             ----------------------------------------------------------
              ELSE
c...            Load impurity density data:
c                WRITE(0,*) 'DEBUG: here 2'
c                   IF (iz.EQ.0) 
c     .                WRITE(0,*) sdlims(1:nks(109),109,iz)
                BACKSPACE 5
                DO ir = 2, nrs
                  IF (idring(ir).EQ.BOUNDARY) CYCLE
                  gdata1(1:nks(ir),ir) = sdlims(1:nks(ir),ir,iz)
                ENDDO
              ENDIF
              gdata => gdata1
c             ----------------------------------------------------------
          CASE DEFAULT 
            CALL ER('Plot987','Unrecognized option',*99)
        ENDSELECT

        ALLOCATE(nv(MAXNKS*MAXNRS))
        ALLOCATE(rv(4,MAXNKS*MAXNKS))
        ALLOCATE(zv(4,MAXNKS*MAXNKS))
        ALLOCATE(cq(MAXNKS*MAXNRS))

c...    Load up data:
        nc = 0
        nv = 0
c       DO ir = irsep, nrs
        DO ir = 1, nrs
          IF (idring(ir).EQ.BOUNDARY) CYCLE
c          IF (ir.LT.irsep.OR.ir.GE.irtrap) CYCLE
          IF ((iopt.EQ.520.OR.iopt.EQ.521.OR.iopt.EQ.522.OR.
     .         iopt.EQ.599).AND.
     .        ir.LT.irsep) CYCLE

          DO ik = 1, nks(ir)
            IF (ir.LT.irsep.AND.ik.EQ.nks(ir)) CYCLE
            id = korpg(ik,ir)
            nc = nc + 1
            nv(nc) = nvertp(id)
            DO i1 = 1, nvertp(id)
              rv(i1,nc) = rvertp(i1,id)
              zv(i1,nc) = zvertp(i1,id)
            ENDDO
            cq(nc) = gdata(ik,ir)
          ENDDO
        ENDDO

      ELSEIF (iopt.LT.90) THEN

        ALLOCATE(tdata(ntri))
        tdata = 0.0
        IF (iopt.EQ.1 ) CALL LoadTriangleData(2,1,1 ,1,tdata,'default')  ! D  density
        IF (.FALSE..AND.iopt.EQ.1) THEN
          DO i1 = 1, ntri 
            IF (ver(tri(i1)%ver(1),1).GT.1.80D0) THEN
c     .          ver(tri(i1)%ver(1),1).LT.1.98D0) THEN
              IF (ver(tri(i1)%ver(1),2).LT.-1.80D0) 
     .          WRITE(0,*) 'LOWER DIVERTOR n_D :',tdata(i1)
              IF (ver(tri(i1)%ver(1),2).GT.-0.3D0.AND.
     .            ver(tri(i1)%ver(1),2).LT. 0.3D0) 
     .          WRITE(0,*) 'MIDPLANE       n_D :',tdata(i1)
              IF (ver(tri(i1)%ver(1),2).GT. 1.80D0) 
     .          WRITE(0,*) 'UPPER DIVERTOR n_D :',tdata(i1)
            ENDIF
          ENDDO
        ENDIF
        IF (iopt.EQ.2 ) CALL LoadTriangleData(3,1,1 ,1,tdata,'default')  ! D2 density
        IF (.FALSE..AND.iopt.EQ.2) THEN
          DO i1 = 1, ntri 
            IF (ver(tri(i1)%ver(1),1).GT. 1.80D0) THEN
              IF (ver(tri(i1)%ver(1),2).LT.-1.80D0) 
     .          WRITE(0,*) 'LOWER DIVERTOR n_D2:',tdata(i1)
              IF (ver(tri(i1)%ver(1),2).GT.-0.3D0.AND.
     .            ver(tri(i1)%ver(1),2).LT. 0.3D0) 
     .          WRITE(0,*) 'MIDPLANE       n_D2:',tdata(i1)
              IF (ver(tri(i1)%ver(1),2).GT. 1.80D0) 
     .          WRITE(0,*) 'UPPER DIVERTOR n_D2:',tdata(i1)
            ENDIF
          ENDDO
        ENDIF

        IF (iopt.EQ.3 ) CALL LoadTriangleData(6,1,7 ,1,tdata,'default')  ! Dalpha
        IF (iopt.EQ.4 ) CALL LoadTriangleData(1,1,5 ,1,tdata,'default')  ! Ionisation
        IF (.FALSE..AND.iopt.EQ.4) THEN
          DO i1 = 1, ntri
            IF (tri(i1)%type.EQ.MAGNETIC_GRID.AND.
     .          tri(i1)%index(2).GE.irsep)
     .        tdata(i1) = 0.0
          ENDDO
        ENDIF

        IF (iopt.EQ.5 ) CALL LoadTriangleData(1,1,10,0,tdata,'default')  ! D+ density
        IF (iopt.EQ.6 ) CALL LoadTriangleData(1,1,15,0,tdata,'default')  ! D+ temperature
        IF (iopt.EQ.8 ) THEN                                             ! D2 pressure from D2 energy density
          CALL LoadTriangleData(3,1,6,1,tdata,'default') 
          fact = ECH * 0.667 * 7.502
c          fact = ECH * 1.0E+06 * 7.502
          DO i1 = 1, ntri 
            tdata(i1) = tdata(i1) * fact
            IF (.FALSE..AND.
     .          ver(tri(i1)%ver(1),1).LT. 0.65D0.AND.
     .          ver(tri(i1)%ver(1),2).LT.-0.62D0) THEN
              WRITE(0,*) 'PRESSURE:',
     .          ver(tri(i1)%ver(1),1),ver(tri(i1)%ver(1),2),tdata(i1)
            ENDIF
c            IF (tri(i1)%type.EQ.MAGNETIC_GRID) tdata(i1) = 0.0
          ENDDO
        ENDIF
        IF (iopt.EQ.9 ) CALL LoadTriangleData(3,1,7,0,tdata,'default')  ! D2 average energy
        IF (iopt.EQ.9) THEN
          DO i1 = 1, ntri 
            tdata(i1) = tdata(i1) * 0.667
            IF (.FALSE..AND.
     .          ver(tri(i1)%ver(1),1).LT. 0.65D0.AND.
     .          ver(tri(i1)%ver(1),2).LT.-0.62D0) THEN
              WRITE(0,*) 'AVG ENG:',
     .          ver(tri(i1)%ver(1),1),ver(tri(i1)%ver(1),2),tdata(i1)
            ENDIF
          ENDDO
        ENDIF

        IF (iopt.EQ.10) THEN 
          ALLOCATE(tdata1(ntri))
          CALL LoadTriangleData(2,1,1,1,tdata ,'default')  ! D  density
          CALL LoadTriangleData(3,1,1,1,tdata1,'default')  ! D2 density
          tdata(1:ntri) = tdata(1:ntri) + 2.0 * tdata1(1:ntri)
          DEALLOCATE(tdata1)
        ENDIF

        IF (iopt.EQ.21) CALL LoadTriangleData(2,1,1,0,tdata,'default')  ! D  density, no volume scaling
        IF (iopt.EQ.22) CALL LoadTriangleData(3,1,1,0,tdata,'default')  ! D2 density, no volume scaling
        IF (iopt.EQ.23) CALL LoadTriangleData(6,2,6,1,tdata,'default')  ! Dgamma (total)

        IF (iopt.EQ.60) CALL LoadTriangleData(5,1,1,1,tdata,'default')  ! Balmer alpha
        IF (iopt.EQ.61) CALL LoadTriangleData(5,2,1,1,tdata,'default')  ! Lyman alpha
        IF (iopt.EQ.62) CALL LoadTriangleData(5,3,1,1,tdata,'default')  ! Lyman beta
        IF (iopt.EQ.63) CALL LoadTriangleData(5,4,1,1,tdata,'default')  ! Lyman gamma
        IF (iopt.EQ.64) CALL LoadTriangleData(5,5,1,1,tdata,'default')  ! Lyman delta
        IF (iopt.EQ.65) CALL LoadTriangleData(5,6,1,1,tdata,'default')  ! Lyman epsilon

        IF (iopt.EQ.80) CALL LoadTriangleData(2,2,1,1,tdata,'default')  ! He(1|1) density
        IF (iopt.EQ.81) CALL LoadTriangleData(2,3,1,1,tdata,'default')  ! He(2|1)
        IF (iopt.EQ.82) CALL LoadTriangleData(2,4,1,1,tdata,'default')  ! He(2|3)

        ALLOCATE(nv(ntri))
        ALLOCATE(rv(3,ntri))
        ALLOCATE(zv(3,ntri))
        ALLOCATE(cq(ntri))

c...    Load up data:
        nc = 0
        nv = 0
        DO i1 = 1, ntri
          IF ((iopt.EQ.5.OR.iopt.EQ.6).AND.
     .        (tri(i1)%type.NE.MAGNETIC_GRID.OR.
     .         tri(i1)%type.EQ.MAGNETIC_GRID.AND.
     .         tri(i1)%index(2).LT.irsep)) CYCLE
          nc = nc + 1
          DO v1 = 1, 3
            nv(nc) = nv(nc) + 1
            rv(v1,nc) = SNGL(ver(tri(i1)%ver(v1),1))
            zv(v1,nc) = SNGL(ver(tri(i1)%ver(v1),2))
          ENDDO
          cq(nc) = tdata(i1)
        ENDDO
      ENDIF

c...  Scale the data:
      IF (ALLOCATED(cq)) THEN
        IF (scalefact.EQ.-1.0) THEN
         cq = LOG10(MAX(1.0E-10,cq))
        ELSE
         cq = cq * scalefact
        ENDIF
      ENDIF

      IF (iopt.LT.90.OR.iopt.GE.500) THEN
c...    Find plotting quantity minimum and maximum:
        setqmin = .TRUE.
        setqmax = .TRUE.
        IF (qmin.NE. HI) setqmin = .FALSE.
        IF (qmax.NE.-HI) setqmax = .FALSE.
        IF (setqmin.OR.setqmax) THEN
          DO i1 = 1, nc
c...        Decide if the cell is within the viewing range:
            inside = .FALSE.
            DO i2 = 1, nv(i1)
              IF (rv(i2,i1).GE.xxmin.AND.rv(i2,i1).LE.xxmax.AND.
     .            zv(i2,i1).GE.yymin.AND.zv(i2,i1).LE.yymax) 
     .          inside=.TRUE.
            ENDDO
            IF (inside) THEN
              IF (setqmin.AND.
     .            cq(i1).NE.0.0.AND.qmin.GT.cq(i1)) qmin=cq(i1)
              IF (setqmax.AND.
     .            cq(i1).NE.0.0.AND.qmax.LT.cq(i1)) qmax=cq(i1)
            ENDIF
          ENDDO
          IF (setqmin.AND.qmin.GT.-1.0E-10) qmin = MAX(qmin,0.01*qmax)
        ENDIF
c        WRITE(0,*) 'QMIN,QMAX:',qmin,qmax

c...    Draw polygons:
        CALL PSPACE(MAP1X,MAP2X,MAP1Y,MAP2Y)
        CALL MAP   (CXMIN,CXMAX,CYMIN,CYMAX)
        CALL HSV
c         hardscale = .TRUE.
        IF (qmin.NE.qmax) THEN  ! This check is required to avoid a strange seg fault - SL, 22/01/2010
          DO i1 = 1, nc
            IF (cq(i1).LT.qmin) cq(i1) = qmin
c            IF (cq(i1).GE.qmin) THEN
c              WRITE(0,*) colouropt,i1,cq(i1),qmin,qmax
              CALL SetCol255_04(colouropt,cq(i1),qmin,qmax)
c              CALL SetCol255_04(colouropt,cq(i1),qmin,qmax)
              CALL FILCOL(255)
              CALL LINCOL(255) 
c              WRITE(6,*) 'PLOT:',rv(1,i1),zv(1,i1),cq(i1),nv(i1)
c              WRITE(6,*) '    :',rv(2,i1),zv(2,i1)
c              WRITE(6,*) '    :',rv(3,i1),zv(3,i1)
c              WRITE(6,*) '    :',rv(4,i1),zv(4,i1)
              CALL PTPLOT(rv(1,i1),zv(1,i1),1,nv(i1),1)
c            ENDIF
          ENDDO
        ENDIF

        IF (scale_set) THEN
c...      Process tags:
          IF (label(1:14).EQ.'<charge state>') 
     .     WRITE(label,'(A,I2,A)') '+',iz,' '//label(15:LEN_TRIM(label))
        ENDIF
        CALL DrawColourScale(scaleopt,colouropt,qmin,qmax,label)

        IF (ALLOCATED(tdata)) DEALLOCATE(tdata)
        DEALLOCATE(nv)
        DEALLOCATE(rv)
        DEALLOCATE(zv)
        DEALLOCATE(cq)
      ENDIF



c...  Draw Vessel and grid outline:

      IF (iopt.GE.500) THEN
c...    Magnetic grid:
        nline = 0
        ALLOCATE(lines2(4*MAXNKS*MAXNRS,2))
        ALLOCATE(lcolour(4*MAXNKS*MAXNRS))

        nver = 0
        CALL ALLOC_VERTEX(10*MAXNKS) ! Borrowed from mod_triangle

        DO ir = 1, nrs
          IF (idring(ir).EQ.BOUNDARY) CYCLE
          DO ik = 1, nks(ir)

            id = korpg(ik,ir)

            IF     (idring(irouts(ik,ir)).EQ.BOUNDARY.OR.
     .              irouts(ik,ir).EQ.ir) THEN  ! The connection map for CDN is not working...
              nver = nver + 1
              ver(nver,1) = DBLE(rvertp(2,id))
              ver(nver,2) = DBLE(zvertp(2,id))
              nver = nver + 1
              ver(nver,1) = DBLE(rvertp(3,id))
              ver(nver,2) = DBLE(zvertp(3,id))
              nline = nline + 1
              lines2(nline,1) = nver - 1
              lines2(nline,2) = nver 
              lcolour(nline) = ncols + 1
            ELSEIF (idring(irins(ik,ir)).EQ.BOUNDARY.OR.ir.EQ.irsep.OR.
     .            (irsep.NE.irsep2.AND.irsep2.GT.0.AND.
     .             ir.EQ.irouts(1                 ,MAX(1,irsep2)).OR.
     .             ir.EQ.irouts(nks(MAX(1,irsep2)),MAX(1,irsep2)))) THEN
              nver = nver + 1
              ver(nver,1) = DBLE(rvertp(1,id))
              ver(nver,2) = DBLE(zvertp(1,id))
              nver = nver + 1
              ver(nver,1) = DBLE(rvertp(4,id))
              ver(nver,2) = DBLE(zvertp(4,id))
              nline = nline + 1
              lines2(nline,1) = nver - 1
              lines2(nline,2) = nver 
              IF (idring(irins(ik,ir)).EQ.BOUNDARY) THEN
                lcolour(nline) = ncols + 1
              ELSE
                lcolour(nline) = ncols + 3
              ENDIF
            ENDIF

          ENDDO

c...      Targets:
          IF (ir.GT.irsep) THEN
            id = korpg(1,ir)
            nver = nver + 1
            ver(nver,1) = DBLE(rvertp(1,id))
            ver(nver,2) = DBLE(zvertp(1,id))
            nver = nver + 1
            ver(nver,1) = DBLE(rvertp(2,id))
            ver(nver,2) = DBLE(zvertp(2,id))
            nline = nline + 1
            lines2(nline,1) = nver - 1
            lines2(nline,2) = nver 
            lcolour(nline) = ncols + 1

            id = korpg(nks(ir),ir)
            nver = nver + 1
            ver(nver,1) = DBLE(rvertp(3,id))
            ver(nver,2) = DBLE(zvertp(3,id))
            nver = nver + 1
            ver(nver,1) = DBLE(rvertp(4,id))
            ver(nver,2) = DBLE(zvertp(4,id))
            nline = nline + 1
            lines2(nline,1) = nver - 1
            lines2(nline,2) = nver 
            lcolour(nline) = ncols + 1
          ENDIF

        ENDDO

c...    Wall:
        DO i1 = 1, wallpts
          IF (wallpt(i1,18).NE.0.0) CYCLE

          nver = nver + 1
          ver(nver,1) = DBLE(wallpt(i1,20))
          ver(nver,2) = DBLE(wallpt(i1,21))
          nver = nver + 1
          ver(nver,1) = DBLE(wallpt(i1,22))
          ver(nver,2) = DBLE(wallpt(i1,23))
          nline = nline + 1
          lines2(nline,1) = nver - 1
          lines2(nline,2) = nver 
          lcolour(nline) = 1
        ENDDO

      ELSEIF (iopt.LT.500) THEN
c      IF (iopt.EQ.99) THEN

        ALLOCATE(lines2(3*ntri,2))
        ALLOCATE(lcolour(3*ntri))
        nline = 0
        DO i1 = 1, ntri
          DO v1 = 1, 3
            v2 = v1 + 1
            IF (v1.EQ.3) v2 = 1
c *TEMP*
            IF     (iopt.EQ.95) THEN

              IF (tri(i1)%type.EQ.VACUUM_GRID.AND.     ! Walls
     .             tri(i1)%sur(v1).NE.0) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 1
                lines2(nline,1)=tri(i1)%ver(v1)
                lines2(nline,2)=tri(i1)%ver(v2)
              ENDIF

              IF (tri(i1)%type.EQ.MAGNETIC_GRID.AND.   ! Targets
     .            tri(i1)%sur(v1).NE.0) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 3
                lines2(nline,1)=tri(i1)%ver(v1)
                lines2(nline,2)=tri(i1)%ver(v2)
              ENDIF

            ELSEIF (iopt.EQ.96) THEN
c             All triangles:
              nline = nline + 1
              lcolour(nline) = ncols + 1
              lines2(nline,1)=tri(i1)%ver(v1)
              lines2(nline,2)=tri(i1)%ver(v2)

c              WRITE(0,*) 'VER:',i1,v1,tri(i1)%ver(v1),v2,tri(i1)%ver(v2)

            ELSEIF (iopt.EQ.97) THEN
              IF (i1.EQ.6268.OR.i1.EQ.6416) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 3
                lines2(nline,1)=tri(i1)%ver(v1)
                lines2(nline,2)=tri(i1)%ver(v2)
              ENDIF

            ELSEIF (iopt.EQ.98) THEN

c              IF (i1.NE.1.AND.i1.NE.380.AND.
c     .            i1.NE.423.AND.i1.NE.448.AND.i1.NE.450) CYCLE
c              IF (tri(i1)%type.EQ.VACUUM_GRID) WRITE(0,*) ' *** VAC!'

              IF (tri(i1)%type.EQ.VACUUM_GRID) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 3
                lines2(nline,1)=tri(i1)%ver(v1)
                lines2(nline,2)=tri(i1)%ver(v2)
              ENDIF

              IF (tri(i1)%type.EQ.MAGNETIC_GRID) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 3
                lines2(nline,1)=tri(i1)%ver(v1)
                lines2(nline,2)=tri(i1)%ver(v2)
              ENDIF

            ELSEIF (iopt.EQ.99) THEN
c...          Magnetic grid triangles + wall surfaces:
              IF (tri(i1)%type.EQ.VACUUM_GRID.AND.
c     .            (tri(i1)%sur(v1).EQ.4.OR.
c     .             tri(i1)%sur(v1).GT.10)) THEN
     .             tri(i1)%sur(v1).NE.0) THEN


                nline = nline + 1
                lcolour(nline) = 1
                lines2(nline,1)=tri(i1)%ver(v1)
                lines2(nline,2)=tri(i1)%ver(v2)
              ENDIF

              IF (tri(i1)%type.EQ.MAGNETIC_GRID.AND.
     .            v1.NE.1) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 1
                lines2(nline,1)=tri(i1)%ver(v1)
                lines2(nline,2)=tri(i1)%ver(v2)
              ENDIF

            ELSE

              IF (tri(i1)%type.EQ.VACUUM_GRID.AND.
c     .            (tri(i1)%sur(v1).EQ.4.OR.
c     .             tri(i1)%sur(v1).GT.10)) THEN
     .             tri(i1)%sur(v1).NE.0) THEN

                nline = nline + 1
                lcolour(nline) = 1  ! ncols + 1
                lines2(nline,1)=tri(i1)%ver(v1)
                lines2(nline,2)=tri(i1)%ver(v2)
              ENDIF

              IF (tri(i1)%type.EQ.MAGNETIC_GRID.AND.
     .            tri(i1)%sur(v1).NE.0) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 3
                lines2(nline,1)=tri(i1)%ver(v1)
                lines2(nline,2)=tri(i1)%ver(v2)
              ENDIF

              IF (tri(i1)%type.EQ.MAGNETIC_GRID.AND.
     .            tri(i1)%index(2).EQ.irsep.AND.
     .            tri(i1)%sideindex(1,v1).EQ.14) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 1
                lines2(nline,1)=tri(i1)%ver(v1)
                lines2(nline,2)=tri(i1)%ver(v2)
              ENDIF

           ENDIF
c            nline = nline + 1
c            lines2(nline,1)=tri(i1)%ver(v1)
c            lines2(nline,2)=tri(i1)%ver(v2)
          ENDDO
        ENDDO
c...    Remove duplicates:
        DO i1 = 1, nline-1
          DO i2 = i1+1, nline
            IF (lines2(i1,1).NE.-999.0.AND.lines2(i2,1).NE.-999.0) THEN
c              IF (lines2(i1,1).LE.0.OR.lines2(i1,2).LE.0.OR.
c     .            lines2(i2,1).LE.0.OR.lines2(i2,2).LE.0) CYCLE
c              IF (lines2(i1,1).GT.288.OR.lines2(i1,2).GT.288.OR.  
c     .            lines2(i2,1).GT.288.OR.lines2(i2,2).GT.288) CYCLE

              IF
c * IMPROVE THIS CHECK! *
     .          ((DABS(ver(lines2(i1,1),1)-
     .                 ver(lines2(i2,2),1)).LT.DTOL.AND.
     .            DABS(ver(lines2(i1,1),2)-
     .                 ver(lines2(i2,2),2)).LT.DTOL.AND.
     .            DABS(ver(lines2(i1,2),1)-
     .                 ver(lines2(i2,1),1)).LT.DTOL.AND.
     .            DABS(ver(lines2(i1,2),2)-
     .                 ver(lines2(i2,1),2)).LT.DTOL).OR.
     .           (DABS(ver(lines2(i1,1),1)-
     .                 ver(lines2(i2,1),1)).LT.DTOL.AND.
     .            DABS(ver(lines2(i1,1),2)-
     .                 ver(lines2(i2,1),2)).LT.DTOL.AND.
     .            DABS(ver(lines2(i1,2),1)-
     .                 ver(lines2(i2,2),1)).LT.DTOL.AND.
     .            DABS(ver(lines2(i1,2),2)-
     .                 ver(lines2(i2,2),2)).LT.DTOL)) THEN
                IF (lcolour(i1).EQ.1) THEN
                  lines2(i1,1) = -999.0
                ELSE
                  lines2(i2,1) = -999.0
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        DO i1 = nline, 1, -1
          IF (lines2(i1,1).EQ.-999.0) THEN
c            WRITE(0,*) 'DELETING:',i1
            DO i2 = i1, nline-1
              lines2(i2,1) = lines2(i2+1,1)
              lines2(i2,2) = lines2(i2+1,2)
              lcolour(i2) = lcolour(i2+1)
            ENDDO
            nline = nline - 1
          ENDIF
        ENDDO
      ENDIF




      IF (iopt.EQ.95) THEN
        CALL DrawGrid(iopt)        
      ELSE
c...    Plot polygons:
        CALL PSPACE (map1x,map2x,map1y,map2y)      
        CALL MAP    (cxmin,cxmax,cymin,cymax)
        lastcolour = -1
        DO i1 = 1, nline
          IF (lastcolour.NE.lcolour(i1)) THEN
            CALL LINCOL(lcolour(i1)) 
            lastcolour = lcolour(i1)
          ENDIF

c          IF (lines2(i1,1).LE.0.OR.lines2(i1,2).LE.0.OR.
c     .        lines2(i2,1).LE.0.OR.lines2(i2,2).LE.0) CYCLE
c          IF (lines2(i1,1).GT.288.OR.lines2(i1,2).GT.288.OR. 
c     .        lines2(i2,1).GT.288.OR.lines2(i2,2).GT.288) CYCLE

          CALL POSITN(SNGL(ver(lines2(i1,1),1)),
     .                SNGL(ver(lines2(i1,1),2)))
          CALL JOIN  (SNGL(ver(lines2(i1,2),1)),
     .                SNGL(ver(lines2(i1,2),2)))
        ENDDO

c...    Frame:
        CALL DrawFrame
      ENDIF





c...     Draw the views on the LOS plot:
c         CALL PSPACE (map1x,map2x,map1y,map2y)
c         CALL MAP    (cxmin,cxmax,0.0,1.0)
c         CALL BROKEN(6,6,6,6)
c         CALL LinCol(1)
c         DO i2 = 1, nline
c           CALL POSITN (lines2(i2),0.0)
c           CALL JOIN   (lines2(i2),1.0)        
c         ENDDO




c...  Add a caption to the plot:
      READ(5,'(A256)') cdum1
      IF   (cdum1(8:14).EQ.'Noframe'.OR.cdum1(8:12).EQ.'noframe'.OR.
     .      cdum1(8:14).EQ.'NOFRAME') THEN
        numplots = numplots + 1
        reset_origin = .FALSE.
      ELSE
        numplots = 0
        BACKSPACE 5
        CALL FRAME
        reset_origin = .TRUE.
      ENDIF

c...  Clear arrays:
      CALL DEALLOC_VERTEX  
      CALL DEALLOC_SURFACE   
      CALL DEALLOC_TRIANGLE

      IF (ALLOCATED(gdata1)) DEALLOCATE(gdata1)

      IF (ALLOCATED(lines2))   DEALLOCATE(lines2)
      IF (ALLOCATED(lcolour)) DEALLOCATE(lcolour)

      RETURN
99    WRITE(0,*) 'IOPT =',iopt
      STOP
      END
