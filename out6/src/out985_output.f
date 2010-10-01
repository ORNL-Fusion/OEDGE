c
c ======================================================================
c
c subroutine: PlotLineShapes
c
      SUBROUTINE PlotLineShapes(istart,iend,MAXPIXEL,npixel,pixel)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

      INTEGER istart,iend,MAXPIXEL,npixel
      TYPE(type_view) :: pixel(MAXPIXEL)

      INCLUDE 'params'
      INCLUDE 'comgra'
      INCLUDE 'colours'
      INCLUDE 'slout'

      INTEGER   iint,ipixel,ispectrum,count,idat,ndat,i1,maxcnt
      REAL      xdat(MAXSPECBIN),ydat(MAXSPECBIN)
      CHARACTER xlabel*36,ylabel*64,label*5,caption*256
      CHARACTER*32 avg_label(10)
c      CHARACTER avg_label*(32)(10) ! Works for ifort, but not pgi


c...  Line shapes:
      count = 0
      DO iint = 1, opt%int_num
        IF (opt%int_type(iint).NE.2) CYCLE
        count = count + 1
        DO ipixel = istart, iend

          ispectrum = pixel(ipixel)%index_spe(iint)

          map1x = 0.10 + REAL(count-1) * 0.45
          map2x = map1x + 0.35 
          map2y = 0.95 - REAL(ipixel-istart) * 0.08 
          map1y = map2y - 0.08

          ndat = MAXSPECBIN
          xdat(1:ndat) = opt%int_wlngthbin(1:ndat,iint) 
          ydat(1:ndat) = spectrum(1:ndat,ispectrum)

          cxmin =  xdat(1)
          cxmax =  xdat(ndat)
          cymin =  HI
          cymax = -HI
          DO idat = 1, ndat
            cymin = MIN(cymin,ydat(idat))
            cymax = MAX(cymax,ydat(idat))
          ENDDO
          cymax = cymax * 1.05

          xlabel = 'lambda (nm)'
          IF (ipixel.NE.iend) xlabel = 'none      ' 
          ylabel = '... (...)'
          CALL GRTSET_TRIM(' ',' ',' ',' ',' ',      ! TITLE,REF,nVIEW,PLANE,glabel,
     .                     cxmin,cxmax,cymin,cymax,
     .                     ' ',xlabel,ylabel,        ! TABLE,XLAB,YLAB,
     .                     0,' ',0,' ',1)            ! 0,smooth,0,ANLY,1)

          WRITE(caption,'(1X,I2,64X)') ipixel
          IF (ipixel.EQ.istart)
     .      WRITE(caption(4: ),'(F10.2)') opt%int_wlngth(iint)
          CALL PLOTST(map1x,0.2*map1y+0.8*map2y,
     .                caption(1:LEN_TRIM(caption)))

          WRITE(label,'(I5,1X)') ipixel

          CALL GRTRAC(xdat,ydat,ndat,label,'LINE',1)

          CALL LINCOL(1)         ! What a mess, but couldn't get things to 
          CALL BROKEN(2,6,2,6)
          CALL PSPACE(MAP1X,MAP2X,MAP1Y,MAP2Y)
          CALL MAP   (CXMIN,CXMAX,CYMIN,CYMAX) 
          CALL LINCOL(1)
          CALL POSITN(0.0,cymin)
          CALL JOIN  (0.0,cymax)

          CALL FULL
          CALL DrawFrame




        ENDDO


      ENDDO

c...  Line-of-sight averages:

      CALL FRAME     

      avg_label(1) = 'avg ne'
      avg_label(2) = 'avg Te'
      avg_label(3) = 'avg nb'
      avg_label(4) = 'avg vb'
      avg_label(5) = 'avg Tb'
      avg_label(6) = 'avg ni'
      avg_label(7) = 'avg vi'
      avg_label(8) = 'avg Ti'

c...  How many plots?
      maxcnt = 0
      DO iint = 1, opt%int_num
        IF (opt%int_type(iint).EQ.3) maxcnt = maxcnt + 1
      ENDDO

      count = 0
      DO iint = 1, opt%int_num
        IF (opt%int_type(iint).NE.3) CYCLE
        count = count + 1

        ispectrum = pixel(ipixel)%index_spe(iint)

        map1x = 0.10                         !  + REAL(count-1) * 0.45
        map2x = map1x + 0.55 
        map2y = 0.95 - REAL(count-1) * 0.12 
        map1y = map2y - 0.12

        ndat = iend - istart + 1
        DO i1 = istart, iend
          xdat(i1-istart+1) = REAL(i1)
        ENDDO
        ydat(1:ndat) = pixel(istart:iend)%average(iint)

        WRITE(0,*) ndat,xdat(1:ndat),ydat(1:ndat)

        cxmin =  xdat(1)
        cxmax =  xdat(ndat)
        cymin =  HI
        cymax = -HI
        DO idat = 1, ndat
          cymin = MIN(cymin,ydat(idat))
          cymax = MAX(cymax,ydat(idat))
        ENDDO
        cymax = cymax * 1.05

        xlabel = 'pixel'
        IF (count.NE.maxcnt) xlabel = 'none      ' 
        ylabel = '... (...)'
        CALL GRTSET_TRIM(' ',' ',' ',' ',' ',      ! TITLE,REF,nVIEW,PLANE,glabel,
     .                   cxmin,cxmax,cymin,cymax,
     .                   ' ',xlabel,ylabel,        ! TABLE,XLAB,YLAB,
     .                   0,' ',0,' ',1)            ! 0,smooth,0,ANLY,1)

        WRITE(caption,'(256X)') 
        WRITE(caption,'(I2,1X,A6,1X,F10.2)') iint,
     .    avg_label(opt%int_average(iint)),opt%int_wlngth(iint)
        CALL PLOTST(map1x,0.2*map1y+0.8*map2y,
     .              caption(1:LEN_TRIM(caption)))

        WRITE(label,'(I5,1X)') ipixel

        CALL GRTRAC(xdat,ydat,ndat,label,'LINE',1)

        CALL LINCOL(1)
        CALL BROKEN(2,6,2,6)
        CALL PSPACE(MAP1X,MAP2X,MAP1Y,MAP2Y)
        CALL MAP   (CXMIN,CXMAX,CYMIN,CYMAX) 
        CALL LINCOL(1)
        CALL POSITN(0.0,cymin)
        CALL JOIN  (0.0,cymax)

        CALL FULL
        CALL DrawFrame

      ENDDO



      RETURN
 9    STOP
      END
c
c ======================================================================
c
c subroutine: GetSchematics
c
      SUBROUTINE GetSchematics(xin,yin,zin,
     .                         mode,MAXSURFACE,MAXPOINTS,
     .                         opt,nobj,obj,
     .                         nsur,npts,hsur,vsur,len1,len2)
      USE MOD_OUT985
      USE mod_eirene04
      IMPLICIT none

c...  Input:
      INTEGER mode,MAXSURFACE,MAXPOINTS,nobj,
     .        nsur,npts(MAXSURFACE),hsur(MAXSURFACE)
      REAL    xin,yin,zin
      REAL*8  vsur(3,MAXPOINTS,0:MAXSURFACE),len1,len2
      TYPE(type_options985) :: opt
      TYPE(type_3D_object)  :: obj(nobj)

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'colours'
      INCLUDE 'pindata'
      INCLUDE 'slout'
      INCLUDE 'slcom'


      INTEGER fp,ik,ir,ikm,id,iobj,isur,ipts,i1,i2,nstart,count,
     .        ik1,ir1,ike,v1,v2,ring
      REAL    r(2),z(2),deltar,deltaz,deltac,deltap,phi,dphi,rval,pval,
     .        frac1,frac2,angle1,angle2,rmid,
     .        brat,rfilament,rfrac

      INTEGER n,index(200000)
      REAL    fraction(200000)
      REAL*8  v(3,200000)

      fp = 0

      SELECTCASE (mode)
        CASE (6)            ! Triangle grid
          CALL LoadTriangles
          DO i1 = 1, ntri
            DO v1 = 1, 3
              v2 = v1 + 1
              IF (v1.EQ.3) v2 = 1
              IF (tri(i1)%type.EQ.MAGNETIC_GRID) THEN
                nsur = nsur + 1 
                hsur(nsur) = -3
                npts(nsur) = 2
                vsur(1,1,nsur) = DBLE(ver(tri(i1)%ver(v1),1))
                vsur(2,1,nsur) = DBLE(ver(tri(i1)%ver(v1),2))
                vsur(3,1,nsur) = 0.0D0
                vsur(1,2,nsur) = DBLE(ver(tri(i1)%ver(v2),1))
                vsur(2,2,nsur) = DBLE(ver(tri(i1)%ver(v2),2))
                vsur(3,2,nsur) = 0.0D0
              ENDIF
            ENDDO
          ENDDO

        CASE (5)            ! Fluid grid
          DO ir = 1, nrs
            DO ik = 1, nks(ir)
              id = korpg(ik,ir)
              DO i1 = 1, nvertp(id)
                i2 = i1 + 1
                IF (i2.GT.nvertp(id)) i2 = 1
                nsur = nsur + 1 
                hsur(nsur) = -3
                npts(nsur) = 2
                vsur(1,1,nsur) = DBLE(rvertp(i1,id))
                vsur(2,1,nsur) = DBLE(zvertp(i1,id))
                vsur(3,1,nsur) = 0.0D0
                vsur(1,2,nsur) = DBLE(rvertp(i2,id))
                vsur(2,2,nsur) = DBLE(zvertp(i2,id))
                vsur(3,2,nsur) = 0.0D0
              ENDDO
            ENDDO
          ENDDO

        CASE (1)            ! Separatrix and vessel wall
          ir = irsep
          DO ik = 1, nks(ir)
            id = korpg(ik,ir)
            nsur = nsur + 1 
            hsur(nsur) = -3
            npts(nsur) = 2
            vsur(1,1,nsur) = DBLE(rvertp(1,id))
            vsur(2,1,nsur) = DBLE(zvertp(1,id))
            vsur(3,1,nsur) = 0.0D0
            vsur(1,2,nsur) = DBLE(rvertp(4,id))
            vsur(2,2,nsur) = DBLE(zvertp(4,id))
            vsur(3,2,nsur) = 0.0D0
          ENDDO

          IF (.FALSE.) THEN
            ir = irsep2
            DO ik = 1, nks(ir)
              id = korpg(ik,ir)
              nsur = nsur + 1 
              hsur(nsur) = -3
              npts(nsur) = 2
              vsur(1,1,nsur) = DBLE(rvertp(1,id))
              vsur(2,1,nsur) = DBLE(zvertp(1,id))
              vsur(3,1,nsur) = 0.0D0
              vsur(1,2,nsur) = DBLE(rvertp(4,id))
              vsur(2,2,nsur) = DBLE(zvertp(4,id))
              vsur(3,2,nsur) = 0.0D0
            ENDDO
          ENDIF

          DO i1 = 1, nvesm
            nsur = nsur + 1 
            hsur(nsur) = 4
            npts(nsur) = 2
            vsur(1,1,nsur) = DBLE(rvesm(i1,1))
            vsur(2,1,nsur) = DBLE(zvesm(i1,1))
            vsur(3,1,nsur) = 0.0D0
            vsur(1,2,nsur) = DBLE(rvesm(i1,2))
            vsur(2,2,nsur) = DBLE(zvesm(i1,2))
            vsur(3,2,nsur) = 0.0D0
          ENDDO

        CASE (2)            ! Field line
          nstart = nsur + 1
          count = 1  ! -1
c          DO ir = irsep, irsep+15, 15 ! DIII-D
c          DO ir = irsep, irsep+19, 20  ! MAST
c          DO ir = irsep, irsep+29, 29  ! MAST
          DO ir = irsep, irsep
            count = count + 1
            ikm = -1
            rmid = -100.0
c...        Find midplane of lines:
            DO ik = 1, nks(ir)
              id = korpg(ik,ir)
              z(1) = 0.5 * (zvertp(1,id) + zvertp(2,id))
              z(2) = 0.5 * (zvertp(3,id) + zvertp(4,id))
              IF (((z(1).LT.z0.AND.z(2).GE.z0).OR.
     .             (z(2).LT.z0.AND.z(1).GE.z0)).AND.
c              IF (((z(1).LT.0.0D0.AND.z(2).GE.0.0).OR.
c     .             (z(2).LT.0.0D0.AND.z(1).GE.0.0)).AND.
     .            rs(ik,ir).GT.rmid) THEN
                ikm  = ik
                rmid = rs(ik,ir)
              ENDIF
            ENDDO
            IF (ikm.EQ.-1) THEN
              WRITE(0,*) 'WHAO! NO MIDPLANE CELL FOUND',ir
              CYCLE 
            ELSE
              WRITE(0,*) 'MIDPLANE:',ikm,ir
            ENDIF
c...        Work from midplane to low IK target:
            phi = 0.0   ! z-axis
            DO ik = ikm, 1, -1
              id = korpg(ik,ir)
              r(1) = 0.5 * (rvertp(3,id) + rvertp(4,id))
              z(1) = 0.5 * (zvertp(3,id) + zvertp(4,id))
              r(2) = 0.5 * (rvertp(1,id) + rvertp(2,id))
              z(2) = 0.5 * (zvertp(1,id) + zvertp(2,id))
              deltar = r(2) - r(1)
              deltaz = z(2) - z(1)
              deltac = ABS(deltaz) / bratio(ik,ir)
              deltap = deltac / rs(ik,ir) * 180.0 / PI
              dphi   = 5.0
              angle1 = 0.0
              DO WHILE (angle1.LT.deltap)
                angle2 = MIN(angle1+dphi,deltap) 
                frac1 = angle1 / deltap
                frac2 = angle2 / deltap
                nsur = nsur + 1
                hsur(nsur) = -count  ! -3 + count 
                npts(nsur) =  2
                vsur(1,1,nsur) = r(1) + frac1 * deltar
                vsur(2,1,nsur) = z(1) + frac1 * deltaz
                vsur(3,1,nsur) = phi + angle1
                vsur(1,2,nsur) = r(1) + frac2 * deltar
                vsur(2,2,nsur) = z(1) + frac2 * deltaz
                vsur(3,2,nsur) = phi + angle2

c                WRITE(0,'(A,4I6,4F12.4)')
c     .            'DATA 1:',nsur,npts(nsur),
c     .            ik,ir,angle1,angle2,phi,deltap
c                WRITE(0,'(3F12.4)')
c     .            vsur(1:3,1,nsur)
c                WRITE(0,'(3F12.4)')
c     .            vsur(1:3,2,nsur)

                angle1 = angle2
              ENDDO 
              phi = phi + deltap
            ENDDO

            phi = 0.0   ! z-axis
            DO ik = ikm+1, nks(ir)
              id = korpg(ik,ir)
              r(1) = 0.5 * (rvertp(1,id) + rvertp(2,id))
              z(1) = 0.5 * (zvertp(1,id) + zvertp(2,id))
              r(2) = 0.5 * (rvertp(3,id) + rvertp(4,id))
              z(2) = 0.5 * (zvertp(3,id) + zvertp(4,id))
              deltar = r(2) - r(1)
              deltaz = z(2) - z(1)
              deltac = ABS(deltaz) / bratio(ik,ir)
              deltap = -deltac / rs(ik,ir) * 180.0 / PI

              angle1 = 0.0
              DO WHILE (angle1.GT.deltap)
                angle2 = MAX(angle1-dphi,deltap) 
                frac1 = angle1 / deltap
                frac2 = angle2 / deltap
                nsur = nsur + 1
                hsur(nsur) = -count  ! -3 + count
                npts(nsur) =  2
                vsur(1,1,nsur) = r(1) + frac1 * deltar
                vsur(2,1,nsur) = z(1) + frac1 * deltaz
                vsur(3,1,nsur) = phi + angle1
                vsur(1,2,nsur) = r(1) + frac2 * deltar
                vsur(2,2,nsur) = z(1) + frac2 * deltaz
                vsur(3,2,nsur) = phi + angle2

c                WRITE(0,'(A,4I6,4F12.4)')
c     .            'DATA 1:',nsur,npts(nsur),
c     .            ik,ir,angle1,angle2,phi,deltap
c                WRITE(0,'(3F12.4)')
c     .            vsur(1:3,1,nsur)
c                WRITE(0,'(3F12.4)')
c     .            vsur(1:3,2,nsur)

                angle1 = angle2
              ENDDO 
              phi = phi + deltap
            ENDDO

          ENDDO   
c...      Convert from r,z,phi to x,y,z (y okay already):
          DO isur = nstart, nsur
            DO i1 = 1, npts(isur)
              rval = vsur(1,i1,isur)
              pval = vsur(3,i1,isur) * PI / 180.0
              vsur(1,i1,isur) = rval * SIN(pval)
              vsur(3,i1,isur) = rval * COS(pval)
c              WRITE(0,'(A,I6,3F12.4)') 
c     .          'DATA 2:',isur,vsur(1:3,i1,isur)
            ENDDO
          ENDDO

        CASE (3)

          nstart = nsur + 1
          count = 2

          rfilament = rho(irsep,SIDE14) 
          WRITE(fp,*) '  RFILAMENT =',rfilament

          DO ir = 2, irwall
            IF (rho(ir,CELL1).EQ.0.0) CYCLE

            IF (rfilament.GE.rho(ir,SIDE14).AND.
     .          rfilament.LT.rho(ir,SIDE23)) THEN
c...          Find midplane cell:
              ikm = -1
              rmid = -100.0 
              DO ik = 1, nks(ir)
                id = korpg(ik,ir)
                z(1) = 0.5 * (zvertp(1,id) + zvertp(2,id))
                z(2) = 0.5 * (zvertp(3,id) + zvertp(4,id))
                IF (((z(1).LT.z0.AND.z(2).GE.z0).OR.
     .               (z(2).LT.z0.AND.z(1).GE.z0)).AND.
c                IF (((z(1).LT.0.0D0.AND.z(2).GE.0.0).OR.
c     .               (z(2).LT.0.0D0.AND.z(1).GE.0.0)).AND.
     .              rs(ik,ir).GT.rmid) THEN
                  ikm  = ik
                  rmid = rs(ik,ir)
                ENDIF
              ENDDO
              IF (ikm.EQ.-1) 
     .          CALL ER('BLAH','No midplane cell found',*99)

              IF (rfilament.LT.rho(ir,CELL1)) THEN
                ir1 = irins(ikm,ir)
                rfrac = -(rfilament      - rho(ir,CELL1)) / 
     .                   (rho(ir1,CELL1) - rho(ir,CELL1))
              ELSE
                ir1 = irouts(ikm,ir)
                rfrac =  (rfilament      - rho(ir,CELL1)) / 
     .                   (rho(ir1,CELL1) - rho(ir,CELL1))
               WRITE(fp,*) rfilament,rho(ir,CELL1),rho(ir1,CELL1)
              ENDIF

              WRITE(fp,*) 'MIDPLANE:',ikm,ir,rfrac,ir1,rho(ir,CELL1)
              EXIT
            ENDIF

          ENDDO


c...      Work from midplane to low IK target:
          phi = 0.0   ! z-axis
          DO ik = ikm, 1, -1
            id = korpg(ik,ir)
            r(1) = 0.5 * (rvertp(3,id) + rvertp(4,id))
            z(1) = 0.5 * (zvertp(3,id) + zvertp(4,id))
            r(2) = 0.5 * (rvertp(1,id) + rvertp(2,id))
            z(2) = 0.5 * (zvertp(1,id) + zvertp(2,id))
            IF (rfrac.LE.0.0) THEN
              ik1 = ikins(ik,ir)
              ir1 = irins(ik,ir)
              brat = (1.0+rfrac)*bratio(ik,ir) - rfrac*bratio(ik1,ir1)
            ELSE
              ik1 = ikouts(ik,ir)
              ir1 = irouts(ik,ir)
              brat = (1.0-rfrac)*bratio(ik,ir) + rfrac*bratio(ik1,ir1)
            ENDIF
            deltar = r(2) - r(1)
            deltaz = z(2) - z(1)
            deltac = ABS(deltaz) / brat
            deltap = deltac / rs(ik,ir) * 180.0 / PI
            dphi   = 5.0
            angle1 = 0.0
            DO WHILE (angle1.LT.deltap)
              angle2 = MIN(angle1+dphi,deltap) 
              frac1 = angle1 / deltap
              frac2 = angle2 / deltap
              nsur = nsur + 1
              hsur(nsur) = -count  
              npts(nsur) =  2
              vsur(1,1,nsur) = r(1) + frac1 * deltar
              vsur(2,1,nsur) = z(1) + frac1 * deltaz
              vsur(3,1,nsur) = phi + angle1
              vsur(1,2,nsur) = r(1) + frac2 * deltar
              vsur(2,2,nsur) = z(1) + frac2 * deltaz
              vsur(3,2,nsur) = phi + angle2
              angle1 = angle2
            ENDDO 
            phi = phi + deltap
          ENDDO
          phi = 0.0   ! z-axis
          ike = nks(ir)
          IF (ir.LT.irsep) ike = ike - 1
          DO ik = ikm+1, ike
            id = korpg(ik,ir)
            r(1) = 0.5 * (rvertp(1,id) + rvertp(2,id))
            z(1) = 0.5 * (zvertp(1,id) + zvertp(2,id))
            r(2) = 0.5 * (rvertp(3,id) + rvertp(4,id))
            z(2) = 0.5 * (zvertp(3,id) + zvertp(4,id))
            IF (rfrac.LE.0.0) THEN
              ik1 = ikins(ik,ir)
              ir1 = irins(ik,ir)
              brat = (1.0+rfrac)*bratio(ik,ir) - rfrac*bratio(ik1,ir1)
            ELSE
              ik1 = ikouts(ik,ir)
              ir1 = irouts(ik,ir)
              brat = (1.0-rfrac)*bratio(ik,ir) + rfrac*bratio(ik1,ir1)
            ENDIF
            deltar = r(2) - r(1)
            deltaz = z(2) - z(1)
            deltac = ABS(deltaz) / brat
            deltap = -deltac / rs(ik,ir) * 180.0 / PI
            angle1 = 0.0
            DO WHILE (angle1.GT.deltap)
              angle2 = MAX(angle1-dphi,deltap) 
              frac1 = angle1 / deltap
              frac2 = angle2 / deltap
              nsur = nsur + 1
              hsur(nsur) = -count 
              npts(nsur) =  2
              vsur(1,1,nsur) = r(1) + frac1 * deltar
              vsur(2,1,nsur) = z(1) + frac1 * deltaz
              vsur(3,1,nsur) = phi + angle1
              vsur(1,2,nsur) = r(1) + frac2 * deltar
              vsur(2,2,nsur) = z(1) + frac2 * deltaz
              vsur(3,2,nsur) = phi + angle2
              angle1 = angle2
            ENDDO 
            phi = phi + deltap
          ENDDO
c...      Convert from r,z,phi to x,y,z (y okay already):
          DO isur = nstart, nsur
            DO i1 = 1, npts(isur)
              rval = vsur(1,i1,isur)
              pval = vsur(3,i1,isur) * PI / 180.0
              vsur(1,i1,isur) = rval * SIN(pval)
              vsur(3,i1,isur) = rval * COS(pval)
            ENDDO
          ENDDO

        CASE (4)
c          CALL TraceFieldLine_DIVIMP(xin,yin,zin,2,3,len1,len2,
         CALL TraceFieldLine_DIVIMP(xin,yin,zin,2,5,len1,len2,
c          CALL TraceFieldLine_DIVIMP(xin,yin,zin,2,1,len1,len2,
     .                               n,v,index,fraction,ring,200000)
          WRITE(0,*) 'FIELD LINE TRACE:',len1,len2
          DO i1 = 1, n-1
            nsur = nsur + 1
            hsur(nsur) = -1
            npts(nsur) =  2
            vsur(1:3,1,nsur) = v(1:3,i1  )
            vsur(1:3,2,nsur) = v(1:3,i1+1)
          ENDDO

        CASE (7)            ! Separatrix only
          ir = irsep
          DO ik = 1, nks(ir)
            id = korpg(ik,ir)
            nsur = nsur + 1 
            hsur(nsur) = -3
            npts(nsur) = 2
            vsur(1,1,nsur) = DBLE(rvertp(1,id))
            vsur(2,1,nsur) = DBLE(zvertp(1,id))
            vsur(3,1,nsur) = 0.0D0
            vsur(1,2,nsur) = DBLE(rvertp(4,id))
            vsur(2,2,nsur) = DBLE(zvertp(4,id))
            vsur(3,2,nsur) = 0.0D0
          ENDDO



        CASE DEFAULT
          IF (mode.LT.0) THEN
            DO iobj = 1, nobj
              IF (obj(iobj)%index.EQ.ABS(mode)) THEN
                DO isur = 1, MAX(obj(iobj)%nsur,obj(iobj)%nside)
                  IF (obj(iobj)%tsur(isur).NE.SP_VESSEL_WALL) CYCLE  ! *TEMP*
                  IF (obj(iobj)%gsur(isur).EQ.GT_TC) THEN
                    nsur = nsur + 1 
                    npts(nsur) = 4
                    hsur(nsur) = -1
                    IF (obj(iobj)%nside.NE.0) THEN
                    ELSE
                      ipts = obj(iobj)%ipts(1,isur)
                      vsur(1,1,nsur) = obj(iobj)%v(1,ipts) 
                      vsur(2,1,nsur) = obj(iobj)%v(2,ipts) 
                      vsur(3,1,nsur) = 0.0D0
                      ipts = obj(iobj)%ipts(2,isur)
                      vsur(1,2,nsur) = obj(iobj)%v(1,ipts) 
                      vsur(2,2,nsur) = obj(iobj)%v(2,ipts) 
                      vsur(3,2,nsur) = 0.0D0
                      vsur(1,3,nsur) = vsur(1,2,nsur)
                      vsur(2,3,nsur) = vsur(2,2,nsur)
                      vsur(3,3,nsur) = 0.00001D0
                      vsur(1,4,nsur) = vsur(1,1,nsur)
                      vsur(2,4,nsur) = vsur(2,1,nsur) 
                      vsur(3,4,nsur) = 0.00001D0
                    ENDIF
                  ENDIF
                ENDDO
              ENDIF
            ENDDO       
          ENDIF

      ENDSELECT

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: SelectTetrahedrons
c
c
      SUBROUTINE SelectTetrahedrons(nsur,npts2,vsur,
     .                              MAXSURFACE,MAXPOINTS,status)    
      USE mod_eirene06_locals
      IMPLICIT none

      INTEGER, INTENT(IN) :: nsur,MAXSURFACE,MAXPOINTS,status,
     .                       npts2(0:MAXSURFACE)
      REAL*8 , INTENT(IN) :: vsur(3,MAXPOINTS,0:MAXSURFACE)


      REAL*8 CalcPerp

      INTEGER iobj,npts,isur,old_nobj
      REAL*8  r(2),phi(2),v(3,1000),p(3),t,pdist


      SELECTCASE(1) 
        CASE (0)

c          IF (.TRUE.) CALL SetupSimpleTetrahedrons
          obj(1:nobj)%segment(1) = 0

          CALL CheckTetrahedronStructure

c          CALL CleanObjectArray

          CALL BuildConnectionMap_New
c          WRITE(0,*) 'MAP:',obj(1)%omap(1:4)
c          WRITE(0,*) '   :',obj(1)%smap(1:4)
c          WRITE(0,*) 'MAP:',obj(2)%omap(1:4)
c          WRITE(0,*) '   :',obj(2)%smap(1:4)

c          obj(1)%segment(1) = 1    
c          obj(2)%segment(1) = 1    
c          CALL DivideTetrahedron(1) 
c          CALL DivideTetrahedron(2) 

          obj(10000:11000)%segment(1) = 1    
          DO iobj = 10000, 11000
c            CALL DivideTetrahedron(iobj) 
          ENDDO

          CALL CheckTetrahedronStructure

          CALL CleanObjectArray

          CALL BuildConnectionMap_New

           
c          obj(1:nobj)%segment(1) = 0
c          obj(54999)%segment(1) = 1    
c          obj(55008)%segment(1) = 1    

          WRITE(0,*) 'NOBJ:',nobj

cSURFACE:       54999      125595      125589      125596      125597
cSURFACE:       55008     -125553     -125589     -125607     -125608

c          DO iobj = 1, nobj
c            IF (obj(iobj)%segment(1).EQ.-1) CYCLE
c            WRITE(0,*) 'MAP:',iobj
c            WRITE(0,*) '   :',obj(iobj)%omap(1:4)
c            WRITE(0,*) '   :',INT(obj(iobj)%smap(1:4))
c          ENDDO

c          obj(1:nobj)%segment(1) = 1

        CASE (1)

c          CALL CheckTetrahedronStructure
c          CALL CleanObjectArray
c          CALL BuildConnectionMap_New


          npts = 2
          r  (1:2) = 0.5D0
          phi(1:2) = 10.0D0
          v(1,1:npts) =  r(1:npts) * DCOS(phi(1:npts)*D_DEGRAD)
          v(3,1:npts) =  r(1:npts) * DSIN(phi(1:npts)*D_DEGRAD)
          v(2,1) =  1.0D0
          v(2,2) = -1.0D0
c          v(1,1) =  0.2D0
c          v(2,1) =  1.0D0
c          v(3,1) =  0.0D0
c          v(1,2) =  0.8D0
c          v(2,2) = -1.0D0
c          v(3,2) =  0.6D0



          WRITE(0,*) 'x:',v(1,1:npts)
          WRITE(0,*) 'y:',v(2,1:npts)
          WRITE(0,*) 'z:',v(3,1:npts)
c          p(1) =  1.0D0
c          p(2) =  0.9999999D0
c          p(3) =  0.0D0
c          pdist = CalcPerp(v(1,1),v(1,2),p,t)


          obj(1:nobj)%segment(1) = 0


          WRITE(0,*) 'NSUR:',nsur


          DO iobj = 1, 0 ! nobj
            IF (MOD(iobj,100000).EQ.0) WRITE(0,*) '  PROCESSING:',iobj
 
            IF (obj(iobj)%index(IND_IR).NE.5) CYCLE
c            IF (obj(iobj)%index(IND_IR).NE.17) CYCLE
            IF (obj(iobj)%segment(1).NE.1) CYCLE

            p(1) = DBLE(obj(iobj)%x)
            p(2) = DBLE(obj(iobj)%y)
            p(3) = DBLE(obj(iobj)%z)

            DO isur = 1, nsur
c              pdist = CalcPerp(v(1,1),v(1,2),p,t)
              pdist = CalcPerp(vsur(1,1,isur),vsur(1,2,isur),p,t)

c              WRITE(0,*) 'PDIST=',isur,iobj,pdist
c              WRITE(0,*) '     =',vsur(1:3,1,isur)
c              WRITE(0,*) '     =',vsur(1:3,2,isur)

              IF (pdist.GE.0.0D0.AND.pdist.LT.0.1D0) THEN
                obj(iobj)%segment(1) = 1  ! *** HACK ***
                EXIT
              ELSE
                obj(iobj)%segment(1) = 0
              ENDIF
            ENDDO
          ENDDO


          old_nobj = nobj
          DO iobj = 1, old_nobj
            IF (obj(iobj)%segment(1).EQ.1) CALL DivideTetrahedron(iobj) 
          ENDDO

c          CALL CheckTetrahedronStructure
c          CALL CleanObjectArray
c          CALL BuildConnectionMap_New

      ENDSELECT

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: TestTetrahedrons
c
c
      SUBROUTINE TestTetrahedrons(option,fname,nsur,npts,vsur,hsur,
     .                            MAXSURFACE,MAXPOINTS,status)
      USE mod_eirene06_locals
      IMPLICIT none

      INTEGER option,nsur,MAXSURFACE,MAXPOINTS,status,
     .        npts(0:MAXSURFACE),hsur(0:MAXSURFACE)
      REAL*8  vsur(3,MAXPOINTS,0:MAXSURFACE)
      CHARACTER fname*(*)

      INTEGER iobj,iside,isrf,i1,iobj1

      CALL LoadObjects(fname(1:LEN_TRIM(fname)),status)
c      CALL LoadObjects('tetrahedrons.001.raw',status)

      WRITE(0,*) 'NOBJ:',nobj

         WRITE(0,*) 'NPTS 490:',nsur,npts(nsur)

      IF (status.EQ.-1) RETURN

c      nsur = 0

c      WRITE(0,*) 'NOBJ:',nobj

      WRITE(0,*) 'NSUR START TETRAHEDRONS:',nsur

      SELECTCASE (option) 
c       ----------------------------------------------------------------   
        CASE (1)
          DO isrf = 1, nsrf
            IF (nsur.EQ.MAXSURFACE) THEN
              WRITE(0,*) 'ERROR TestTetrahedrons: MAXSURFACE exceeded'
              RETURN
            ENDIF
            nsur = nsur + 1
            hsur(nsur) = 301
            npts(nsur) = srf(isrf)%nvtx
            DO i1 = 1, npts(nsur)
              vsur(1:3,i1,nsur) = vtx(1:3,srf(isrf)%ivtx(i1))
            ENDDO
          ENDDO
c       ----------------------------------------------------------------   
        CASE (2)  ! Vessel wall
          WRITE(0,*) 'MESSAGE TestTetrahedrons: OPTION=2'
          DO iobj = 1, nobj
c           IF (grp(obj(iobj)%group)%origin.EQ.GRP_VACUUM_GRID) CYCLE     
            DO iside = 1, obj(iobj)%nside
              isrf = obj(iobj)%iside(iside)
              isrf = ABS(isrf)
              IF (obj(iobj)%omap(iside).NE.0) CYCLE
c               IF (srf(isrf)%index(IND_SURFACE).EQ.0) CYCLE
              DO WHILE(isrf.GT.0)
                IF (nsur.GE.MAXSURFACE-10000) THEN
                  WRITE(0,*) 'ERROR TestTetra...: MAXSURFACE exceeded'
                  RETURN
                ENDIF
                nsur = nsur + 1
                IF (grp(obj(iobj)%group)%origin.EQ.GRP_VACUUM_GRID) 
     .            hsur(nsur) = -2 ! 301 ! -2 ! 301
                IF (grp(obj(iobj)%group)%origin.EQ.GRP_MAGNETIC_GRID) 
     .            hsur(nsur) = -2 ! 301 ! -3 ! 301
                npts(nsur) = srf(isrf)%nvtx
                IF (npts(nsur).NE.3) STOP 'sdgfsdgsdsd'
                DO i1 = 1, npts(nsur)
                  vsur(1:3,i1,nsur) = vtx(1:3,srf(isrf)%ivtx(i1))
                ENDDO
                isrf = srf(isrf)%link
                IF (isrf.NE.0) STOP 'USING LINK CODE - BRAVO'
              ENDDO
            ENDDO
          ENDDO
c       ----------------------------------------------------------------   
        CASE (3)
c...      Search for tetrahedrons that are close to a particular trajectory:
          CALL SelectTetrahedrons(nsur,npts,vsur,
     .                            MAXSURFACE,MAXPOINTS,status)

          DO iobj = 1, nobj
            DO iside = 1, obj(iobj)%nside
              isrf = ABS(obj(iobj)%iside(iside))

c              IF (iobj.EQ.55008) WRITE(0,*) ' ORIGIN:',
c     .           grp(obj(iobj)%group)%origin.EQ.GRP_VACUUM_GRID,
c     .           grp(obj(iobj)%group)%origin.EQ.GRP_MAGNETIC_GRID

              IF (grp(obj(iobj)%group)%origin.EQ.GRP_VACUUM_GRID.AND.
     .            (srf(isrf)%index(IND_SURFACE).NE.8.OR..TRUE.)) CYCLE

              IF (grp(obj(iobj)%group)%origin.EQ.GRP_MAGNETIC_GRID.AND.
     .            obj(iobj)%segment(1).EQ.0) CYCLE  ! *** HACK ***
     
              isrf = ABS(obj(iobj)%iside(iside))

c              DO WHILE(isrf.GT.0)
              IF (nsur.GE.MAXSURFACE) THEN
                WRITE(0,*) 'ERROR TestTetra...: MAXSURFACE exceeded'
                RETURN
              ENDIF
              nsur = nsur + 1
              IF (grp(obj(iobj)%group)%origin.EQ.GRP_VACUUM_GRID) 
     .          hsur(nsur) = -2 ! 301
              IF (grp(obj(iobj)%group)%origin.EQ.GRP_MAGNETIC_GRID) 
     .          hsur(nsur) = -3 ! 301
              npts(nsur) = srf(isrf)%nvtx
              DO i1 = 1, npts(nsur)
                vsur(1:3,i1,nsur) = vtx(1:3,srf(isrf)%ivtx(i1))
              ENDDO
c                isrf = srf(isrf)%link
c                IF (isrf.NE.0) STOP 'USING LINK CODE - BRAVO'
c              ENDDO
            ENDDO
          ENDDO
c       ----------------------------------------------------------------   
        CASE (4)
          DO iobj = 1, nobj
c              IF (grp(obj(iobj)%group)%origin.EQ.GRP_MAGNETIC_GRID.AND.
c     .            (obj(iobj)%index(IND_IR).NE.17)) CYCLE
c     .            (obj(iobj)%index(IND_IR).NE.5)) CYCLE
c            IF (grp(obj(iobj)%group)%origin.NE.GRP_MAGNETIC_GRID) CYCLE     
c            IF (grp(obj(iobj)%group)%origin.NE.GRP_VACUUM_GRID) CYCLE     
c            IF (obj(iobj)%index(IND_IR).NE.5) CYCLE
c            IF (obj(iobj)%segment(1).EQ.0) CYCLE

            DO iside = 1, obj(iobj)%nside
              isrf = obj(iobj)%iside(iside)
c              IF (isrf.LT.0) CYCLE  ! *** NOTE! ***

c              iobj1 = obj(iobj)%omap(iside)
c              IF (iobj1.NE.0.AND.obj(iobj)%segment(1).EQ.0) THEN
c                IF (grp(obj(iobj1)%group)%origin.EQ.GRP_MAGNETIC_GRID) 
c     .            CYCLE          
c              ENDIF
 
              isrf = ABS(isrf)
              IF (srf(isrf)%index(IND_SURFACE).NE.4.AND.     ! Midplane port
     .            srf(isrf)%index(IND_SURFACE).NE.6.AND. 
     .            srf(isrf)%index(IND_SURFACE).NE.8.AND. 
     .            srf(isrf)%index(IND_SURFACE).NE.26) CYCLE
c              IF (srf(isrf)%index(IND_SURFACE).NE.4.AND.     ! Upper port
c     .            srf(isrf)%index(IND_SURFACE).NE.22) CYCLE

c              IF (grp(obj(iobj)%group)%origin.EQ.GRP_VACUUM_GRID.AND.
c     .            srf(isrf)%index(IND_SURFACE).NE.8) CYCLE

c              isrf = ABS(obj(iobj)%iside(iside))
              DO WHILE(isrf.GT.0)
                IF (nsur.GE.MAXSURFACE-10000) THEN
                  WRITE(0,*) 'ERROR TestTetra...: MAXSURFACE exceeded'
                  RETURN
                ENDIF
                nsur = nsur + 1
                IF (grp(obj(iobj)%group)%origin.EQ.GRP_VACUUM_GRID) 
     .            hsur(nsur) = 301 ! -2 ! 301
                IF (grp(obj(iobj)%group)%origin.EQ.GRP_MAGNETIC_GRID) 
     .            hsur(nsur) = 301 ! -3 ! 301
                npts(nsur) = srf(isrf)%nvtx
                IF (npts(nsur).NE.3) STOP 'sdgfsdgsdsd'
                DO i1 = 1, npts(nsur)
                  vsur(1:3,i1,nsur) = vtx(1:3,srf(isrf)%ivtx(i1))
                ENDDO
                isrf = srf(isrf)%link
                IF (isrf.NE.0) STOP 'USING LINK CODE - BRAVO'
              ENDDO
            ENDDO
          ENDDO          
c       ----------------------------------------------------------------   
        CASE DEFAULT
          CALL WN('TestTetrahedrons','Unknown option')
      ENDSELECT

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: SolidPlot985
c
c
      SUBROUTINE SolidPlot985(opt,nobj,obj,iplot)
      USE MOD_OUT985
      USE mod_filament
      USE mod_options
      IMPLICIT none

c...  Input:
      TYPE(type_options985) :: opt
      INTEGER, INTENT(IN) :: nobj, iplot
      TYPE(type_3D_object) :: obj(nobj)

      INCLUDE 'params'
      INCLUDE 'comgra'
      INCLUDE 'cgeom'
      INCLUDE 'colours'
      INCLUDE 'slout'
      INCLUDE 'slcom'

      INTEGER   MAXSURFACE       ,MAXPOINTS
c      PARAMETER(MAXSURFACE=5000000,MAXPOINTS=10)
      PARAMETER(MAXSURFACE=310000,MAXPOINTS=10)

      INTEGER FindMidplaneCell
      REAL    FindSeparatrixRadius


      INTEGER nsur,iobj,isur,isrf,i1,i2,i3,ipts,ipts1,ipts2,pass,nx,ny,
     .        ix,iy,iver,count,fp,ntmp,nlight,isrf1,isrf2,status,
     .        icolour,last_icolour,ir,ifilament,ikm,id,ibin,nbin,nlist,
     .        option,sub_option,iplot1,idum1,nsur_solid
      LOGICAL cont,solid
      REAL    XXMIN,XXMAX,YYMIN,YYMAX,theta,theta1,theta2,
     .        xsur(MAXPOINTS),ysur(MAXPOINTS),frac,
     .        xin,yin,zin,rsep,rhoval,
     .        xcen,ycen,xnear,ynear,angle(3),delr,dely,delphi
      REAL*8  a(3),b(3),n(3),p1(3,6),p2(3,6),light(3,10),view(3),p3(3),
     .        mat(3,3),deltax,deltay,res,xs1,xs2,ys1,ys2,
     .        dangle,ang,x1,z1,r,x,y,z,
     .        sv(3),lv(3,10),nv(3),vv(3),ndotv,ldots,
     .        mindsur,maxdsur,frac1,frac2,limit1,limit2,len1,len2
      CHARACTER fname*512,buffer*2048,cdum1*512


      INTEGER, ALLOCATABLE :: npts(:),hsur(:),dlist(:),ilist(:)
      REAL*8,  ALLOCATABLE :: csur(:,:),vsur(:,:,:),bsur(:,:),dsur(:)


      CALL setup_col(16,5)
   
      CALL THICK2(5)

      solid = .FALSE.
      nsur_solid = MAXSURFACE + 1

      ALLOCATE(npts(0:MAXSURFACE))
      ALLOCATE(hsur(0:MAXSURFACE))
      ALLOCATE(dsur(0:MAXSURFACE))
      ALLOCATE(csur(3,0:MAXSURFACE))
      ALLOCATE(bsur(3,0:MAXSURFACE))
      ALLOCATE(vsur(3,MAXPOINTS,0:MAXSURFACE))

      nsur = 0
      csur = 0.0D0
      vsur = 0.0D0

c...  Distance from observer to surface center:
      view(1) = 0.0D0  
      view(2) = 0.0D0
      view(3) = 5.0D0  

      light = 1.0D0

      nlight = 1
      light(1,1) =  10.0D0
      light(2,1) =   5.0D0 
      light(3,1) =   2.0D0 
      light(1,2) =   2.0D0
      light(2,2) =  10.0D0 
      light(3,2) =   2.0D0 
      light(1,3) =  -5.0D0
      light(2,3) =  10.0D0 
      light(3,3) =  -2.0D0 

      xxmin =  HI
      xxmax = -HI
      yymin =  HI
      yymax = -HI

      angle = 0.0

      DO iplot1 = iplot+1, opt%nplots
        READ(opt%plots(iplot1),*) cdum1
        IF (cdum1(1:4).EQ.'item') THEN
              xin = 0.0
              yin = 0.0
              zin = 0.0
              len1 = 0.0D0
              len2 = 0.0D0

c...      Add geometry item to list of objects to be plotted:
          READ(opt%plots(iplot1),*) cdum1,option
          SELECTCASE (option)
c           ------------------------------------------------------------
            CASE (001) ! LOS integration geometry objects
              IF (nsur_solid.EQ.MAXSURFACE+1) nsur_solid = nsur + 1
              READ(opt%plots(iplot1),*) cdum1,option,sub_option,icolour

              DO iobj = 1, nobj
                DO isur = 1, MAX(obj(iobj)%nsur,obj(iobj)%nside)
c                  IF (obj(iobj)%tsur(isur).NE.SP_VESSEL_WALL) CYCLE  ! *TEMP*

                  IF     (obj(iobj)%gsur(isur).EQ.GT_TC) THEN
c...                Create polygons for toroidally continuous surfaces:
                    IF (obj(iobj)%nside.NE.0) THEN
                      isrf1 = obj(iobj)%iside(isur,1)
                      isrf2 = obj(iobj)%iside(isur,2)
                      r = 0.0
                    ELSE
                      isrf1 = 1
                      isrf2 = 1
                      x = obj(iobj)%v(1,obj(iobj)%ipts(1,isur)) 
                      z = 0.0D0
                      r = DSQRT( x**2 + z**2 )
                    ENDIF
                    IF (r.LT.0.5D0) THEN
c                      dangle =  2.0D0 * DBLE(PI) / 180.0D0
                      dangle =  5.0D0 * DBLE(PI) / 180.0D0
                    ELSE
                      dangle = 30.0D0 * DBLE(PI) / 180.0D0
                    ENDIF
                    DO ang = 0.0D0, 359.0D0 * DBLE(PI) / 180.0D0, dangle
                      DO isrf = isrf1, isrf2
                        IF (obj(iobj)%nside.NE.0) THEN
                          p1(1,1) = vtx(1,srf(isrf)%ivtx(1))
                          p1(2,1) = vtx(2,srf(isrf)%ivtx(1))
                          p1(3,1) = p1(1,1) * DTAN(-0.5D0*dangle)
                          p2(1,1) = p1(1,1)
                          p2(2,1) = p1(2,1)
                          p2(3,1) = p2(1,1) * DTAN(+0.5D0*dangle)
                          p1(1,2) = vtx(1,srf(isrf)%ivtx(2))
                          p1(2,2) = vtx(2,srf(isrf)%ivtx(2))
                          p1(3,2) = p1(1,2) * DTAN(-0.5D0*dangle)
                          p2(1,2) = p1(1,2)
                          p2(2,2) = p1(2,2)
                          p2(3,2) = p2(1,2) * DTAN(+0.5D0*dangle)
                        ELSE
                          p1(1,1) =obj(iobj)%v(1,obj(iobj)%ipts(1,isur)) 
                          p1(2,1) =obj(iobj)%v(2,obj(iobj)%ipts(1,isur)) 
                          p1(3,1) = p1(1,1) * DTAN(-0.5D0*dangle)
                          p2(1,1) = p1(1,1)
                          p2(2,1) = p1(2,1)
                          p2(3,1) = p2(1,1) * DTAN(+0.5D0*dangle)
                          p1(1,2) =obj(iobj)%v(1,obj(iobj)%ipts(2,isur)) 
                          p1(2,2) =obj(iobj)%v(2,obj(iobj)%ipts(2,isur)) 
                          p1(3,2) = p1(1,2) * DTAN(-0.5D0*dangle)
                          p2(1,2) = p1(1,2)
                          p2(2,2) = p1(2,2)
                          p2(3,2) = p2(1,2) * DTAN(+0.5D0*dangle)
                        ENDIF
c...                    Rotate vertices toroidally:
                        DO i3 = 1, 2
                          x1 = p1(1,i3)
                          z1 = p1(3,i3)
                          p1(1,i3) = DCOS(ang) * x1 - DSIN(ang) * z1
                          p1(3,i3) = DSIN(ang) * x1 + DCOS(ang) * z1
                          x1 = p2(1,i3)
                          z1 = p2(3,i3)
                          p2(1,i3) = DCOS(ang) * x1 - DSIN(ang) * z1
                          p2(3,i3) = DSIN(ang) * x1 + DCOS(ang) * z1
                        ENDDO
                        nsur = nsur + 1
                        npts(nsur) = 4
                        hsur(nsur) = opt%obj_colour(obj(iobj)%index)
c                        IF (MOD(nsur,1000).EQ.0)
c     .                    WRITE(0,*) '  - ',iobj,obj(iobj)%index,
c     .                                      hsur(nsur),r
                        csur(1:3,nsur) = 0.0D0
                        vsur(1:3,4,nsur) = p1(1:3,1)
                        vsur(1:3,3,nsur) = p2(1:3,1)
                        vsur(1:3,2,nsur) = p2(1:3,2)
                        vsur(1:3,1,nsur) = p1(1:3,2)
                      ENDDO
                    ENDDO
                  ELSEIF (obj(iobj)%gsur(isur).EQ.GT_TD) THEN
c...                Load polygons from toroidally descrete surfaces:
c                    IF (obj(iobj)%tsur(isur).NE.SP_GRID_BOUNDARY) CYCLE
c     . ..  .            obj(iobj)%ik.NE.1) CYCLE
c     . ..  .            obj(iobj)%ik.NE.nks(obj(iobj)%ir)) CYCLE
c     . ..  .            obj(iobj)%ir.NE.2) CYCLE
c                    IF (obj(iobj)%ivolume.NE.1) CYCLE   ! *** HACK *** 
c                    IF (obj(iobj)%ir     .NE.9) CYCLE 
c                    IF (isur.EQ.1) 

c                      WRITE(0,*) 'OBJ:',iobj,obj(iobj)%ik, 
c     .                  obj(iobj)%tsur(1).EQ.SP_GRID_BOUNDARY,
c     .                  obj(iobj)%nside,
c     .                  obj(iobj)%iside(isur,1:2)

                    IF (obj(iobj)%nside.NE.0) THEN
                       DO isrf = obj(iobj)%iside(isur,1),
     .                          obj(iobj)%iside(isur,2)
                        IF (nsur+1.GT.MAXSURFACE) 
     .                    CALL ER('DrawSoildPlot','Array bounds',*99)
                        nsur = nsur + 1
                        npts(nsur) = srf(isrf)%nvtx
                        hsur(nsur) = opt%obj_colour(obj(iobj)%index)
                        csur(1:3,nsur) = 0.0D0
                        DO i1 = 1, npts(nsur)
                          vsur(1:3,i1,nsur) =vtx(1:3,srf(isrf)%ivtx(i1))
c                          WRITE(0,*) 'WORKING',isrf,i1,vsur(1:3,i1,nsur)         
                        ENDDO
                      ENDDO
                    ELSE
                      nsur = nsur + 1
                      npts(nsur) = obj(iobj)%npts(isur)
                      hsur(nsur) = opt%obj_colour(obj(iobj)%index)
                      csur(1:3,nsur) = 0.0D0
                      DO i1 = 1, npts(nsur)
                        iver = obj(iobj)%ipts(i1,isur)
                        vsur(1:3,i1,nsur) = obj(iobj)%v(1:3,iver)
                      ENDDO
                    ENDIF
                  ELSE
                    CALL ER('SolidPlot985','Unknown surface type',*99)
                  ENDIF
            
                ENDDO
              ENDDO
c           ------------------------------------------------------------
            CASE (002) ! 2D annotation
              READ(opt%plots(iplot1),*) cdum1,option,sub_option,icolour
              IF     (sub_option.EQ.4) THEN
                READ(opt%plots(iplot1),*) cdum1,(idum1,i1=1,3),
     .                                    delr,dely,delphi
                xin = FindSeparatrixRadius(1) + delr
                yin = 0.0
                zin = 0.0
                len1 = 1.0E+20
                len2 = 1.0E+20
              ELSEIF (sub_option.EQ.6) THEN
                READ(opt%plots(iplot1),*) cdum1,(idum1,i1=1,3),
     .                                    delr,dely,delphi
                xin = FindSeparatrixRadius(1) + delr
                len1 = DBLE(dely)
                len2 = DBLE(delphi) 
                sub_option = 4
              ENDIF
              CALL GetSchematics(xin,yin,zin,
     .                           sub_option,MAXSURFACE,MAXPOINTS,
     .                           opt,nobj,obj,
     .                           nsur,npts,hsur,vsur,len1,len2)
              npts(nsur) = 2 ! This appears necessary -- compiler bug? or something naught somewhere...
              IF (sub_option.EQ.4) THEN
              ENDIF
c           ------------------------------------------------------------
            CASE (003) ! 3D flux-tube
              READ(opt%plots(iplot1),*) cdum1,option,sub_option,icolour,
     .                                  xin
              yin = 0.0
              zin = 0.0
c              CALL CalcTubeDimensions(xin,yin,zin,
c     .               MAXSURFACE,MAXPOINTS,nsur,npts,hsur,vsur)
c              CALL ProcessFluxTube(xin,yin,zin,
c     .               MAXSURFACE,MAXPOINTS,nsur,npts,hsur,vsur)
c              npts(nsur) = 2 ! This appears necessary -- compiler bug? or something naught somewhere...
c           ------------------------------------------------------------
            CASE (020) ! Filaments
              READ(opt%plots(iplot1),*) cdum1,option,sub_option,icolour
              IF (sub_option.EQ.1) THEN
                READ(opt%plots(iplot1),*) cdum1,idum1,idum1,idum1,fname
                CALL LoadFilamentData(fname(1:LEN_TRIM(fname)),status)
                IF (status.EQ.-1) 
     .            CALL ER('SolidPlot985','Unable to open filament '//
     .                    'data file',*98)
              ELSEIF (sub_option.EQ.2) THEN
                opt_fil%opt  = 2
                opt_fil%start_time = 0.0
                CALL DefineFilaments
                CALL SetupFilaments(0.0D0)
              ELSE
                CALL ER('SolidPlot985','Unrecognised filament '//
     .                  'option',*99)
              ENDIF
              WRITE(6,*) '==FIELD LINE DATA:'
              DO ifilament = 1, nfilament
                IF (filament(ifilament)%status.EQ.0) CYCLE
                DO i1 = 1, filament(ifilament)%nvtx
                  IF (sub_option.EQ.1) THEN
                    len1 = filament(ifilament)%lcell(1,i1)   ! *** HACK ***
                    len2 = filament(ifilament)%lcell(2,i1)
                  ELSE
                    len1 = 1000.0D0
                    len2 = 1000.0D0
                  ENDIF
                  xin = SNGL(filament(ifilament)%vtx(1,i1))
                  yin = 0.0
                  zin = SNGL(filament(ifilament)%vtx(3,i1))
c                  WRITE(0,*) 'FILAMENT',ifilament,i1,xin,zin
                  rsep = FindSeparatrixRadius(1)
                  rhoval = SQRT(xin**2 + zin**2) - rsep
c                  WRITE(0,'(A,I6,3F10.6,2(2X,2F12.6))') ' ',
c     .              ifilament,rsep,rhoval,
c     .              len1,len2,xin,zin
                  CALL GetSchematics(xin,yin,zin,
     .                            4,MAXSURFACE,MAXPOINTS,opt,nobj,obj,
     .                            nsur,npts,hsur,vsur,len1,len2)
                  npts(nsur) = 2 ! This appears necessary -- compiler bug? or something naught somewhere...
c...              Output:
                  WRITE(6,'(A,2I6,3F10.6,2X,2F12.6)') ' ',
     .              ifilament,i1,rsep,rhoval,rsep+rhoval,
     .              len1,len2
                ENDDO
              ENDDO
c           ------------------------------------------------------------
            CASE (021) ! Test tetrahedrons
              READ(opt%plots(iplot1),*) cdum1,option,sub_option,icolour,
     .                                  fname
              IF (sub_option.NE.2.AND.
     .            nsur_solid.EQ.MAXSURFACE+1) nsur_solid = nsur + 1
              CALL TestTetrahedrons(sub_option,fname,nsur,npts,vsur,
     .                              hsur,MAXSURFACE,MAXPOINTS,status)
              CALL ClearTetArrays
              WRITE(0,*) 'TETRAHEDRON FILE STATUS:',status
              IF (status.EQ.-1) RETURN
            CASE DEFAULT
          ENDSELECT

        ELSEIF (cdum1(1:4).EQ.'zoom') THEN
          READ(opt%plots(iplot1),*) cdum1,xcen,ycen,xnear,ynear
          xxmin = xcen - xnear
          xxmax = xcen + xnear
          yymin = ycen - ynear
          yymax = ycen + ynear

        ELSEIF (cdum1(1:4).EQ.'axis') THEN
          READ(opt%plots(iplot1),*) cdum1,(angle(i1),i1=1,3)
          angle = angle * PI / 180.0
        ELSE
          EXIT 
        ENDIF   
      ENDDO      

      IF (.TRUE.) THEN

      ELSEIF (.FALSE.) THEN

        fp = 99
        fname = 'test-p3.raw'
        OPEN(fp,FILE=fname(1:LEN_TRIM(fname)),
     .       FORM='FORMATTED',STATUS='OLD',ERR=98)     
        
        count = 0
        DO WHILE (count.LE.4396)
c        DO WHILE (count.LE.341)

          READ(fp,'(A2048)') buffer         
          WRITE(0,*) 'buffer: '//buffer(1:10)
          IF (buffer(1:3).EQ.'obj') THEN
            DO WHILE (.TRUE.) 
              READ(fp,'(A2048)',END=10) buffer 
c              WRITE(0,*) 'buffer: '//buffer(1:100)
c              WRITE(0,*) 'buffer: '//buffer(101:200)
              IF (buffer(1:3).NE.'obj') THEN
                nsur = nsur + 1
                npts(nsur) = 3
                csur(1:3,nsur) = 0.0D0
                READ(buffer,*) 
     .            (vsur(1,i1,nsur),
     .             vsur(3,i1,nsur),
     .             vsur(2,i1,nsur),i1=1,3)
c                READ(buffer,*) (vsur(1:3,i1,nsur),i1=1,3)
c...            Convert units from mm to m:
                DO i1 = 1, npts(nsur)
                  vsur(1:3,i1,nsur) = vsur(1:3,i1,nsur) * 0.001
                ENDDO
c                WRITE(0,'(1X,A,3F10.4)') '  DATA:',vsur(1:3,1,nsur)
c                WRITE(0,'(1X,A,3F10.4)') '  DATA:',vsur(1:3,2,nsur)
c                WRITE(0,'(1X,A,3F10.4)') '  DATA:',vsur(1:3,3,nsur)
              ELSE
                WRITE(0,*) '  GETTING OUT',count,nsur
                BACKSPACE fp
                EXIT
              ENDIF
            ENDDO
            count = count + 1
          ENDIF
          
        ENDDO
 10     CONTINUE
      ENDIF


c.... Estimate center of surface (before rotation/translation):
      csur = 0.0D0
      DO isur = 1, nsur
c        WRITE(0,*) 'NPTS:',isur,npts(isur)
        DO ipts = 1, npts(isur)
          csur(1:3,isur) = csur(1:3,isur) + vsur(1:3,ipts,isur)
        ENDDO
        csur(1:3,isur) = csur(1:3,isur) / DBLE(REAL(npts(isur)))
      ENDDO


      IF (.TRUE.) THEN
        WRITE(0,*) '    SELECTING'

c...    Decide which surfaces are to be deleted:
        DO isur = 1, nsur
          IF (hsur(isur).LT.-1.OR.npts(isur).EQ.2) CYCLE

          IF (MOD(ABS(hsur(isur))/100,10).EQ.1.OR.
     .        MOD(ABS(hsur(isur))/100,10).EQ.2) CYCLE

          r = DSQRT( csur(1,isur)**2 + csur(3,isur)**2 )

c          IF (csur(3,isur).GT.0.01D0.AND.                          ! MAST
c     .        r.GT.0.75.AND.DABS(csur(2,isur)).LT.1.5) THEN        
c          IF (r.GT.1.20.AND.                                       ! DIII-D
c     .        csur(2,isur).GT.-1.0.AND.csur(2,isur).LT.1.2) THEN
c...        Tag for deletion:
c            npts(isur) = -999
c          ENDIF
        ENDDO

c...    Delete:
        WRITE(0,*) '    DELETING  NSUR:',nsur
        DO i1 = nsur, 1, -1
          IF (npts(i1).EQ.-999) THEN
            DO i2 = i1, nsur-1
              npts(i2) = npts(i2+1)
              hsur(i2) = hsur(i2+1)
              DO i3 = 1, npts(i2)
                vsur(1:3,i3,i2) = vsur(1:3,i3,i2+1)
              ENDDO
            ENDDO
            nsur = nsur - 1
          ENDIF
        ENDDO
        WRITE(0,*) '              NSUR:',nsur
      ENDIF


c.... Estimate center of surface:
      csur = 0.0D0
      DO isur = 1, nsur
        DO ipts = 1, npts(isur)
          csur(1:3,isur) = csur(1:3,isur) + vsur(1:3,ipts,isur)
        ENDDO
        csur(1:3,isur) = csur(1:3,isur) / DBLE(REAL(npts(isur)))
      ENDDO


       IF (angle(1).NE.0.0.OR.angle(2).NE.0.0.OR.angle(3).NE.0.0) THEN
         WRITE(0,*) '    TRANSFORMING'

c...    Warp:

c...    Rotate/tilt/swing/displace:
c            ...or...
c       Displace/swing/tilt/rotate:

        DO isur = 1, nsur
          DO ipts = 1, npts(isur)
c...        Rotate about z-axis (roll):
            CALL Calc_Transform2(mat,0.0D0,1,0)  ! move outside loop...
cc            angle = DBLE(-10.0*PI/180.0)
c            angle = DBLE(  0.0*PI/180.0)
            CALL Calc_Transform2(mat,DBLE(angle(1)),3,1)
            CALL Transform_Vect(mat,vsur(1,ipts,isur))

            IF (isur.EQ.1.AND.ipts.EQ.1) THEN    ! Lame
              DO i1 = 1, nlight
                CALL Transform_Vect(mat,light(1,i1))
              ENDDO
            ENDIF

c...        Rotate about x-axis (tilt):               
            CALL Calc_Transform2(mat,0.0D0,1,0)
cc            angle = DBLE(+00.0*PI/180.0)
c            angle = DBLE( 20.0*PI/180.0)
cc            angle = DBLE(+90.0*PI/180.0)
            CALL Calc_Transform2(mat,DBLE(angle(2)),1,1)
            CALL Transform_Vect(mat,vsur(1,ipts,isur))

            IF (isur.EQ.1.AND.ipts.EQ.1) THEN
              DO i1 = 1, nlight
                CALL Transform_Vect(mat,light(1,i1))
              ENDDO
            ENDIF

c...        Rotate about y-axis (swing):
            CALL Calc_Transform2(mat,0.0D0,1,0)
cc            angle = DBLE(-90.0*PI/180.0)
c            angle = DBLE( 00.0*PI/180.0)
            CALL Calc_Transform2(mat,DBLE(angle(3)),2,1)
            CALL Transform_Vect(mat,vsur(1,ipts,isur))

            IF (isur.EQ.1.AND.ipts.EQ.1) THEN
              DO i1 = 1, nlight
                CALL Transform_Vect(mat,light(1,i1))
              ENDDO
            ENDIF

c...        Translate:
c            vsur(1,ipts,isur) = vsur(1,ipts,isur) + 0.0D0
c            vsur(2,ipts,isur) = vsur(2,ipts,isur) + 0.0D0
c            vsur(3,ipts,isur) = vsur(3,ipts,isur) - 10.0D0

c...        Calculate distance of each point from the viewing plane:
 
c            WRITE(0,*) 'FRAC:',frac

c            frac = 1.0 / (1.0 + DABS(vsur(3,ipts,isur)))**1
c            vsur(1:2,ipts,isur) = vsur(1:2,ipts,isur) * frac
          ENDDO
        ENDDO    
      ENDIF



c.... Estimate center of surface:
      csur = 0.0D0
      DO isur = 1, nsur
        DO ipts = 1, npts(isur)
          csur(1:3,isur) = csur(1:3,isur) + vsur(1:3,ipts,isur)
        ENDDO
        csur(1:3,isur) = csur(1:3,isur) / DBLE(REAL(npts(isur)))
      ENDDO


      IF (.TRUE..OR.solid) THEN
        WRITE(0,*) '    CHECKING VISIBILITY'

c...    Decide which surfaces are visible:
        bsur = 0.0D0
        DO isur = 1, nsur
          IF (isur.LT.nsur_solid) CYCLE
          IF (hsur(isur).LT.0.OR.npts(isur).EQ.2) CYCLE
          a(1:3) = vsur(1:3,1,isur) - vsur(1:3,2,isur) 
          b(1:3) = vsur(1:3,3,isur) - vsur(1:3,2,isur) 
c...      Surface normal:
          n(1) =  a(2) * b(3) - a(3) * b(2)
          n(2) =  a(3) * b(1) - a(1) * b(3)
          n(3) =  a(1) * b(2) - a(2) * b(1) 
          IF (.FALSE.) THEN
c...        Check if surface is visible from both sides, and if yes, flip the normal if surface
c           invisible:
          ELSEIF (MOD(ABS(hsur(isur))/100,10).NE.1.AND.
     .            MOD(ABS(hsur(isur))/100,10).NE.3.AND.
     .            n(3).LE.0.0D0) THEN
c...        Tag for deletion:
            npts(isur) = -999
          ELSE
c...        Store normalized surface normal for lighting calculation later on (or do it now?):
            bsur(1:3,isur) = n(1:3) / DSQRT(n(1)**2 + n(2)**2 + n(3)**2)
          ENDIF
        ENDDO


c ***ALSO CLIP BASED ON WHETHER OR NOT THE SURFACE IS BEHIND THE OBSERVER***

c...    Delete surfaces:
        WRITE(0,*) 'NSUR:',nsur
        DO i1 = nsur, 1, -1
          IF (npts(i1).EQ.-999) THEN
            DO i2 = i1, nsur-1
              npts(i2) = npts(i2+1)
              hsur(i2) = hsur(i2+1)
              bsur(1:3,i2) = bsur(1:3,i2+1)
              DO i3 = 1, npts(i2)
                vsur(1:3,i3,i2) = vsur(1:3,i3,i2+1)
              ENDDO
            ENDDO
            nsur = nsur - 1
          ENDIF
        ENDDO
        WRITE(0,*) 'NSUR:',nsur

      ENDIF


c.... Estimate center of surface:
      csur = 0.0D0
      DO isur = 1, nsur
        DO ipts = 1, npts(isur)
          csur(1:3,isur) = csur(1:3,isur) + vsur(1:3,ipts,isur)
        ENDDO
        csur(1:3,isur) = csur(1:3,isur) / DBLE(REAL(npts(isur)))
      ENDDO


c...  Subdivide surface to give a minimum surface size, so that the chances of
c     incorrect distance ordering of surfaces is reduced:
      IF (.TRUE..OR.solid) THEN
        res = 0.20D0
        ntmp = nsur
        DO i1 = 1, ntmp
          IF (i1.LT.nsur_solid) CYCLE
          IF (npts(i1).EQ.4) THEN

            deltax = MAX(DSQRT((vsur(1,1,i1) - vsur(1,2,i1))**2 +
     .                         (vsur(2,1,i1) - vsur(2,2,i1))**2 +
     .                         (vsur(3,1,i1) - vsur(3,2,i1))**2),
     .                   DSQRT((vsur(1,3,i1) - vsur(1,4,i1))**2 +
     .                         (vsur(2,3,i1) - vsur(2,4,i1))**2 +
     .                         (vsur(3,3,i1) - vsur(3,4,i1))**2))
            deltay = MAX(DSQRT((vsur(1,2,i1) - vsur(1,3,i1))**2 +
     .                         (vsur(2,2,i1) - vsur(2,3,i1))**2 +
     .                         (vsur(3,2,i1) - vsur(3,3,i1))**2),
     .                   DSQRT((vsur(1,4,i1) - vsur(1,1,i1))**2 +
     .                         (vsur(2,4,i1) - vsur(2,1,i1))**2 +
     .                         (vsur(3,4,i1) - vsur(3,1,i1))**2))

            IF (deltax.GT.res.OR.deltay.GT.res) THEN

c...          Store vertex data:
              DO i2 = 1, npts(i1)
                p1(1:3,i2) = vsur(1:3,i2,i1)
              ENDDO

              nx = INT(deltax/res) + 1
              ny = INT(deltay/res) + 1

              DO ix = 1, nx
                DO iy = 1, ny
                  xs1 = DBLE(ix - 1) / DBLE(nx)
                  xs2 = DBLE(ix    ) / DBLE(nx)
                  ys1 = DBLE(iy - 1) / DBLE(ny)
                  ys2 = DBLE(iy    ) / DBLE(ny)

                  IF (ix.EQ.1.AND.iy.EQ.1) THEN
                    isur = i1
                  ELSE
                    nsur = nsur + 1
                    isur = nsur
                    npts(isur)     = npts(i1) 
                    hsur(isur)     = hsur(i1) 
                    bsur(1:3,isur) = bsur(1:3,i1) 
                  ENDIF

                  vsur(1:3,1,isur) = (1.0D0-xs1)*p1(1:3,1)+xs1*p1(1:3,2)
                  vsur(1:3,2,isur) = (1.0D0-xs2)*p1(1:3,1)+xs2*p1(1:3,2)
                  vsur(1:3,3,isur) = (1.0D0-xs2)*p1(1:3,4)+xs2*p1(1:3,3)
                  vsur(1:3,4,isur) = (1.0D0-xs1)*p1(1:3,4)+xs1*p1(1:3,3)

                  DO i2 = 1, npts(i1)
                    p2(1:3,i2) = vsur(1:3,i2,isur)
                  ENDDO          

                  vsur(1:3,1,isur) = (1.0D0-ys1)*p2(1:3,1)+ys1*p2(1:3,4)
                  vsur(1:3,2,isur) = (1.0D0-ys1)*p2(1:3,2)+ys1*p2(1:3,3)
                  vsur(1:3,3,isur) = (1.0D0-ys2)*p2(1:3,2)+ys2*p2(1:3,3)
                  vsur(1:3,4,isur) = (1.0D0-ys2)*p2(1:3,1)+ys2*p2(1:3,4)

                ENDDO
              ENDDO

            ENDIF
 
          ELSE
c            CALL ER('DrawSolidPlot','Can only subdivide quadrangles',*99)
          ENDIF
        ENDDO

      ENDIF

      WRITE(0,*) 'NSUR:',nsur

c...  Decide how close each surface is to the viewing location and sort accordingly:

c.... Estimate center of surface:
      csur = 0.0D0
      DO isur = 1, nsur
        DO ipts = 1, npts(isur)
          csur(1:3,isur) = csur(1:3,isur) + vsur(1:3,ipts,isur)
        ENDDO
        csur(1:3,isur) = csur(1:3,isur) / DBLE(npts(isur))
      ENDDO

c...  Want to find out which surfaces are closest to the observer, but in this scheme, 
c     where the objects are rotated and translated with the observer fixed, it is
c     enough just to take the z-coordinate, i.e. the eye is assumed to be along the
c     z-axis (at infinity?):
      DO i1 = 1, nsur
        dsur(i1) = csur(3,i1)
c        dsur(i1) = DSQRT((csur(1,i1) - view(1))**2 + 
c     .                   (csur(2,i1) - view(2))**2 +
c     .                   (csur(3,i1) - view(3))**2)
      ENDDO


c...  Sort surfaces, listing (and therefore drawing) the farthest surface from the
c     observer first:
      IF (.TRUE.) THEN
c...    Find range of distances:
        mindsur =  1.0D+20
        maxdsur = -1.0D+20
        DO isur = 1, nsur        
          mindsur = MIN(mindsur,dsur(isur))
          maxdsur = MAX(maxdsur,dsur(isur))
        ENDDO
        mindsur = mindsur - 0.00001D0
        maxdsur = maxdsur + 0.00001D0
c        WRITE(0,*) 'MIN,MAX:',mindsur,maxdsur
c...    Bin objects by distance from viewer:
        nbin = 1000
        nlist = 0
        ALLOCATE(ilist(nbin))
        ALLOCATE(dlist(0:nsur))
        ilist(1) = 1
        DO ibin = 1, nbin-1
          frac1 = DBLE(ibin - 1) / DBLE(nbin - 1)
          frac2 = DBLE(ibin    ) / DBLE(nbin - 1)
          limit1 = mindsur + frac1 * (maxdsur - mindsur) 
          limit2 = mindsur + frac2 * (maxdsur - mindsur) 
          DO isur = 1, nsur
            IF (dsur(isur).GE.limit1.AND.dsur(isur).LT.limit2) THEN
              nlist = nlist + 1
              dlist(nlist) = isur
            ENDIF 
          ENDDO
          ilist(ibin+1) = nlist + 1
c          WRITE(0,*) 'LIST:',ilist(ibin+1)-ilist(ibin),nlist,nsur
        ENDDO
c...    Sort the list of objects:
        DO ibin = 1, nbin-1
          IF (MOD(ibin,nbin/10).EQ.0) WRITE(0,*) 'PROCESSING BIN ',ibin
          DO i1 = ilist(ibin), ilist(ibin+1)-2
            DO i2 = i1+1, ilist(ibin+1)-1
              IF (dsur(dlist(i1)).GT.dsur(dlist(i2))) THEN
c...            Swap:
                dlist(0)  = dlist(i1)
                dlist(i1) = dlist(i2)
                dlist(i2) = dlist(0)
              ENDIF
            ENDDO
          ENDDO
        ENDDO

        DEALLOCATE(ilist)

      ENDIF

c.... Estimate center of surface:
c      csur = 0.0D0
c      DO isur = 1, nsur
c        DO ipts = 1, npts(isur)
c          csur(1:3,isur) = csur(1:3,isur) + vsur(1:3,ipts,isur)
c        ENDDO
c        csur(1:3,isur) = csur(1:3,isur) / DBLE(REAL(npts(isur)))
c      ENDDO

c...  *TEMP* -- how do I do this properly?
      
      IF (xxmin.EQ.HI) THEN
        DO isur = 1, nsur
          DO ipts = 1, npts(isur)
            xxmin = MIN(xxmin,SNGL(vsur(1,ipts,isur)))
            xxmax = MAX(xxmax,SNGL(vsur(1,ipts,isur)))
            yymin = MIN(yymin,SNGL(vsur(2,ipts,isur)))
            yymax = MAX(yymax,SNGL(vsur(2,ipts,isur)))
          ENDDO
        ENDDO
        IF (xxmin.LT.0.0) xxmin = xxmin * 1.10
        IF (xxmin.GT.0.0) xxmin = xxmin * 0.90
        IF (xxmax.GT.0.0) xxmax = xxmax * 1.10
        IF (yymin.LT.0.0) yymin = yymin * 1.10
        IF (yymin.GT.0.0) yymin = yymin * 0.90
        IF (yymax.GT.0.0) yymax = yymax * 1.10
        IF (yymax.LT.0.0) yymax = yymax * 0.90
      ENDIF

      map1x = 0.05
      map1y = 0.05

      IF (xxmax-xxmin.GT.yymax-yymin) THEN
c      IF (.FALSE..AND.xxmax-xxmin.GT.yymax-yymin) THEN
        map2x = map1x + 0.90
        map2y = map1y + 0.90 * (yymax - yymin) / (xxmax - xxmin)
      ELSE
        map2x = map1x + 0.90 * (xxmax - xxmin) / (yymax - yymin)
        map2y = map1y + 0.90 
      ENDIF
c      map2x = map1x + 0.90
c      map2y = map1y + 0.90

      CALL GRTSET_TRIM(' ',' ',' ',' ',' ',    ! TITLE,REF,nVIEW,PLANE,glabel,
     .                 xXMIN,xXMAX,yYMIN,yYMAX,
     .                 ' ',' ',' ',            ! TABLE,XLAB,YLAB,
     .                 0,' ',0,' ',1)          ! 0,smooth,0,ANLY,1)

c      WRITE(0,*) 'CX,CY:',cxmin,cxmax,cymin,cymax

c...  Draw polygons:
      CALL PSPACE(MAP1X,MAP2X,MAP1Y,MAP2Y)
      CALL MAP   (CXMIN,CXMAX,CYMIN,CYMAX)  ! CX set in GRTSET_TRIM
      CALL LINCOL(1) 

      vv(1:3) = view (1:3) / DSQRT( view(1)**2+ view(2)**2+ view(3)**2)
      DO i1 = 1, nlight
        lv(1:3,i1) = light(1:3,i1) / 
     .               DSQRT(light(1,i1)**2+light(2,i1)**2+light(3,i1)**2)
      ENDDO

      CALL RGB

      last_icolour = -999

c      IF (solid) THEN
c...    Solid:

      WRITE(0,*)
c      WRITE(0,*) 'NSUR,NSUR_SOLID:',nsur,nsur_solid

      DO i2 = 1, nlist
        isur = dlist(i2)

c        IF (vsur(3,1,isur).GT.0.0D0) CYCLE
         
        IF (isur.GE.nsur_solid) THEN
          last_icolour = -999
c          WRITE(0,*) 'ISUR:',i2,dlist(i2),nsur,nlist

c...      Determine light shading:   *** NEEDS WORK! ***

c...      Find angle between surface and light vector:
          IF (.TRUE.) THEN
            IF (hsur(isur).LT.0.OR.npts(isur).EQ.2) THEN
              ldots = 1.0D0
              r = DSQRT( csur(1,isur)**2 + csur(3,isur)**2 )
              IF (hsur(isur).EQ.-3.OR.hsur(isur).EQ.-2)   
c     .          ldots = MAX(0.8D0,MIN(1.0D0,(r - 1.0D0) / 0.5D0))   ! DIII-D
     .          ldots = MAX(0.8D0,MIN(1.0D0,(r - 0.2D0) / 0.7D0))   ! MAST
            ELSE
              nv(1:3) = bsur(1:3,isur)
              ndotv = nv(1) * vv(1) + nv(2) * vv(2) + nv(3) * vv(3)
              sv(1:3) = 2.0D0 * ndotv * nv(1:3) - vv(1:3)
              ldots = 0.0D0
              DO i1 = 1, nlight
                ldots = ldots + 
     .                  lv(1,i1)*sv(1) + lv(2,i1)*sv(2) + lv(3,i1)*sv(3)
              ENDDO
              IF (ldots.LT.0.0) ldots = 0.0
              IF (ldots.GT.1.0) ldots = 1.0
    
c              ldots = 0.6 * ldots + 0.4
              ldots = 0.7 * ldots + 0.3
c              ldots = 0.90 * ldots + 0.10

              IF (.FALSE.) THEN           
                IF (MOD(isur,1000).EQ.0.OR.isur.GT.15400) 
     .            WRITE(0,*) '  - ',isur,dsur(isur),
     .                       hsur(isur)
              ELSEIF (.FALSE..AND.isur.LT.100) THEN
                WRITE(0,*) '  - ',isur,dsur(isur),ldots
                WRITE(0,*) '    ',vv
                WRITE(0,*) '    ',nv
                WRITE(0,*) '    ',sv
                WRITE(0,*) '  l1',lv(1:3,1)
                WRITE(0,*) '  l2',lv(1:3,2)
                WRITE(0,*) '  l3',lv(1:3,3)
              ENDIF

            ENDIF
c...        Set colour:
            SELECTCASE (MOD(ABS(hsur(isur)),100))  ! *** HACK *** No idea why the -1 is necessary...
c            SELECTCASE (MOD(ABS(hsur(isur-1)),100))  ! *** HACK *** No idea why the -1 is necessary...
              CASE (1)
                CALL ColSet(SNGL(ldots),0.0,0.0,255)
              CASE (2)
                CALL ColSet(0.0,SNGL(ldots),0.0,255)
              CASE (3)
                CALL ColSet(0.0,0.0,SNGL(ldots),255)
              CASE (4)  ! Black
                CALL ColSet(0.0        ,0.0,0.0,255)
              CASE DEFAULT
                CALL ColSet(SNGL(ldots),SNGL(ldots),SNGL(ldots),255)
            ENDSELECT
            CALL FILCOL(255)
            CALL LINCOL(255) 
          ELSE
          ENDIF
c...      Plot polygon:
          DO ipts = 1, npts(isur)
            xsur(ipts) = SNGL(vsur(1,ipts,isur))  
            ysur(ipts) = SNGL(vsur(2,ipts,isur))
          ENDDO
          IF (npts(isur).EQ.2) THEN
            CALL PTPLOT(xsur,ysur,1,1,1)             ! NEED THIS FOR COLOUR TO TAKE, FOR SOME REASON
            CALL PTJOIN(xsur,ysur,1,npts(isur),1)
          ELSE
            CALL PTPLOT(xsur,ysur,1,npts(isur),1)
          ENDIF
c        ENDDO

        ELSE
c...    Wireframe:
c        DO isur = 1, nsur

          IF (hsur(isur).LT.0) THEN
            icolour = ncols + ABS(hsur(isur))
          ELSE
            icolour = hsur(isur) 
          ENDIF
          IF (icolour.NE.last_icolour) THEN
            CALL LINCOL(icolour) 
            last_icolour = icolour
          ENDIF
c          CALL LINCOL(ncols+1) 

          DO ipts1 = 1, npts(isur)  
c          DO ipts1 = 1, npts(nsur)  
            ipts2 = ipts1 + 1 
            IF (ipts2.GT.npts(isur)) ipts2 = 1
c            IF (ipts2.GT.npts(nsur)) ipts2 = 1
            IF (ipts1.EQ.2.AND.npts(isur).EQ.2) CYCLE

            p1(1:3,1) = vsur(1:3,ipts1,isur)
            p2(1:3,1) = vsur(1:3,ipts2,isur)

            CALL POSITN (SNGL(p1(1,1)),SNGL(p1(2,1)))
            CALL JOIN   (SNGL(p2(1,1)),SNGL(p2(2,1))) 
          ENDDO

        ENDIF

      ENDDO




      CALL DrawFrame

      WRITE(0,*) 'DONE DRAWING'

      WRITE(0,*) 'DEALLOCATING MEMORY'

c...  Clear memory:
      IF (ALLOCATED(npts))  DEALLOCATE(npts)
      IF (ALLOCATED(hsur))  DEALLOCATE(hsur)
      IF (ALLOCATED(dsur))  DEALLOCATE(dsur)
      IF (ALLOCATED(csur))  DEALLOCATE(csur)
      IF (ALLOCATED(bsur))  DEALLOCATE(bsur)
      IF (ALLOCATED(vsur))  DEALLOCATE(vsur)
      IF (ALLOCATED(dlist)) DEALLOCATE(dlist)

      WRITE(0,*) 'LEAVING DRAWING ROUTINE'

      RETURN
 98   WRITE(0,*) 'ERROR DrawSolidPlot: File not found'
      WRITE(0,*) '   FNAME= "'//fname(1:LEN_TRIM(fname))//'"'
 99   STOP
      END
c
c ======================================================================
c
c subroutine: Output985
c
      SUBROUTINE WireframePlot985(iopt,MAXPIXEL,npixel,pixel,image,
     .                            iplot)
      USE mod_out985
      USE mod_out985_variables
      USE mod_interface
      IMPLICIT none

      INTEGER iopt,MAXPIXEL,npixel         ! Put this into _variables... 
      INTEGER, INTENT(IN) :: iplot
      TYPE(type_view) :: pixel(MAXPIXEL)
      REAL*8 image(1100,1100)

      INCLUDE 'params'      
      INCLUDE 'comgra'
      INCLUDE 'colours'
      INCLUDE 'slout'

      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost

      CHARACTER glabel*512,caption*1024
      CHARACTER TITLE*80,TITLE2*80,JOB*72,GRAPH*80,GRAPH1*256
      CHARACTER*36 REF,NVIEW,PLANE,ANLY,TABLE,ZLABS(-2:MAXIZS+1)
      CHARACTER*36 XLAB
      CHARACTER*72 YLAB
      CHARACTER*72 SMOOTH
      CHARACTER*256 cdum1

      INTEGER axis1(3)
      REAL    angle1(3)

      INTEGER CH1

      INTEGER   nxbin,nybin,ipixel,nbin,isrf,isrf1,isrf2,isid,
     .          idet,nplot,iobj,i1,i2,i3,i4,ix,iy,iplot1
      REAL      deltar,deltaz,qmin,qmax,qval,ang1,ang2,ang,
     .          frac1,frac2,xcen,ycen,xnear,ynear,count
      REAL      XXMIN,XXMAX,YYMIN,YYMAX,dangle
      REAL*8    angle,p1(3,6),p2(3,6),x1,x2,z1,mat(3,3)

      INTEGER dat1,dat2
      REAL, ALLOCATABLE :: xdat(:),ydat(:)
      CHARACTER xlabel*256,ylabel*256,tag_x*2,tag_y*2,file*512


c *TEMP*
      INTEGER, ALLOCATABLE :: nv(:)
      REAL   , ALLOCATABLE :: rv(:,:),zv(:,:),cq(:)



c...  PLOT A: 3D surface geometry plot for a given X,Y,Z along a particular view
c     -Translate and rotate as necessary:
c     -Project onto the R,Z plane:

      glabel = job  (CH1(job  ):LEN_TRIM(job  ))//'     '//
     .         graph(CH1(graph):LEN_TRIM(graph))

      ! jdemod - changed the : to ) to balance the parentheses in the format specifiers
      !          NOTE: these write statements will produce lines of 512 characters each
      !                right-padded with blanks 
      WRITE(nview ,'(512(A))') (' ',i1=1,LEN(nview )) 
      WRITE(plane ,'(512(A))') (' ',i1=1,LEN(plane )) 
      WRITE(anly  ,'(512(A))') (' ',i1=1,LEN(anly  )) 
      WRITE(table ,'(512(A))') (' ',i1=1,LEN(table )) 
      WRITE(xlab  ,'(512(A))') (' ',i1=1,LEN(xlab  )) 
      WRITE(ylab  ,'(512(A))') (' ',i1=1,LEN(ylab  )) 
      WRITE(smooth,'(512(A))') (' ',i1=1,LEN(smooth)) 
      WRITE(graph ,'(512(A))') (' ',i1=1,LEN(graph )) 
      WRITE(job   ,'(512(A))') (' ',i1=1,LEN(job   )) 

      XLAB = 'none'
      YLAB = 'none'
c      XLAB = '   R  (M)'
c      YLAB = '   Z  (M)'

      CALL THICK (1)
      CALL THICK2(1)

      slopt = 1

c...  Set the "zoom" for the plot:
      xxmin =  0.0  ! Not sure if these are any good...
      xxmax =  6.0
      yymin = -3.0
      yymax =  3.0
      DO iplot1 = iplot+1, opt%nplots
        READ(opt%plots(iplot1),*) cdum1
        IF (cdum1(1:4).EQ.'zoom') THEN
          READ(opt%plots(iplot1),*) cdum1,xcen,ycen,xnear,ynear
          xxmin = xcen - xnear
          xxmax = xcen + xnear
          yymin = ycen - ynear
          yymax = ycen + ynear
        ELSE
          EXIT
        ENDIF   
      ENDDO      

      axis1 (1) = 1
      axis1 (2) = 2
      axis1 (3) = 3
      angle1(1) = 0.0
      angle1(2) = 0.0
      angle1(3) = 0.0

      nplot = 0
      iplot1 = iplot1 - 1
      DO WHILE (.TRUE.) 
c...    Wireframe:
        iplot1 = iplot1 + 1
        READ(opt%plots(iplot1),*) cdum1
c        WRITE(0,*) 'AXIS:',cdum1(1:LEN_TRIM(cdum1))
        IF (cdum1(1:4).EQ.'axis') THEN
          READ(opt%plots(iplot1),*) cdum1,(axis1(i1),angle1(i1),i1=1,3)
          nplot = nplot + 1
          IF (nplot.EQ.5) THEN
           nplot = 1
           CALL FRAME
          ENDIF
        ELSE
          EXIT
        ENDIF   

c...    Setup transformation matrix:
        CALL Calc_Transform2(mat,0.0D0,1,0)
        DO i1 = 1, 3
          angle = DBLE(angle1(i1)*PI/180.0)
          CALL Calc_Transform2(mat,angle,axis1(i1),1)
        ENDDO

c...    Use GRTSET_TRIM:
        slopt4 = 1
c...    Stopping resizing of scale font in ghost1.o6a:
        iopt_ghost = 1

        SELECTCASE (nplot)
          CASE(1)
            map1x = 0.07
            map1y = 0.55
          CASE(2)
            map1x = 0.57
            map1y = 0.55
          CASE(3)
            map1x = 0.07
            map1y = 0.08
          CASE(4)
            map1x = 0.57
            map1y = 0.08
          CASEDEFAULT
        ENDSELECT

        map2x = map1x + 0.40
        map2y = map1y + 0.40

        CALL GRTSET_TRIM(TITLE,REF,nVIEW,PLANE,glabel,
     >                   xXMIN,xXMAX,yYMIN,yYMAX,
     .                   TABLE,XLAB,YLAB,
     .                   0,smooth,0,ANLY,1)
 
c...    Draw polygons:
        CALL PSPACE(MAP1X,MAP2X,MAP1Y,MAP2Y)
        CALL MAP   (CXMIN,CXMAX,CYMIN,CYMAX)

        count = 0.0
        DO iobj = 1, nobj
c          IF (iobj.NE.1035) CYCLE  ! *TEMP*
c          IF (iobj.NE.471) CYCLE  ! *TEMP*
c          IF (iobj.NE.471.AND.iobj.NE.472.AND.iobj.NE.458) CYCLE  ! *TEMP*

          CALL LINCOL(ncols+obj(iobj)%colour) 

          DO isid = 1, MAX(obj(iobj)%nsur,obj(iobj)%nside)

c             IF (obj(iobj)%tsur(isid).NE.SP_GRID_BOUNDARY) CYCLE

c            IF (isid.NE.2) CYCLE

c            IF (obj(iobj)%flag(isid).NE.-1.AND.
c     .          obj(iobj)%tsur(isid).NE.SP_VESSEL_WALL) CYCLE  ! *TEMP*

            IF (obj(iobj)%tsur(isid).NE.SP_VESSEL_WALL.AND.
     .          obj(iobj)%tsur(isid).NE.SP_GRID_BOUNDARY) CYCLE

c              WRITE(6,*) 'BOUNDARY?',obj(iobj)%ik,obj(iobj)%ir

c            IF (obj(iobj)%tsur(isid).NE.SP_VESSEL_WALL) CYCLE   ! Recent filter...

c            IF (obj(iobj)%imap(1,isid).NE.0) CYCLE
c            IF (obj(iobj)%ik.NE.16) CYCLE


            IF (obj(iobj)%type.EQ.OP_INTEGRATION_VOLUME) THEN
              count = count + 
     .                1.0 / REAL(MAX(obj(iobj)%nsur,obj(iobj)%nside))
c              IF (count.GT.100.0) CYCLE
c              IF (count.GT.15000.0) CYCLE
            ENDIF

c *** NEED TO STORE LINE SEGMENTS AND DELETE DUPLICATES... 

            IF (obj(iobj)%nside.NE.0) THEN
              isrf1 = obj(iobj)%iside(isid,1)
              isrf2 = obj(iobj)%iside(isid,2)
            ELSE
              isrf1 = 1
              isrf2 = 1
            ENDIF

            IF (obj(iobj)%gsur(isid).EQ.GT_TC) THEN
              dangle =  15.0 * PI / 180.0
              IF (.TRUE..OR.nplot.EQ.1) THEN
                ang1 = 0.0
                ang2 = 1.0 * PI / 180.0
              ELSE
                ang1 = 0.0
                ang2 = 359.0 * PI / 180.0
              ENDIF
c              DO ang = 0.0, 1.0 * PI / 180.0, dangle
              DO ang = ang1, ang2, dangle
c                IF (ang.GT.0.0.AND.
c     .              obj(iobj)%type.NE.OP_INTEGRATION_VOLUME) CYCLE
c                DO isur = 1, 1  ! THIS IS HERE FOR WHEN SIDES ARE USED!
                DO isrf = isrf1, isrf2
                  IF (obj(iobj)%nside.NE.0) THEN
                    p1(1,1) = vtx(1,srf(isrf)%ivtx(1))
                    p1(2,1) = vtx(2,srf(isrf)%ivtx(1))
                    p1(3,1) = p1(1,1) * DTAN(-0.5D0*dangle)
                    p2(1,1) = p1(1,1)
                    p2(2,1) = p1(2,1)
                    p2(3,1) = p2(1,1) * DTAN(+0.5D0*dangle)

                    p1(1,2) = vtx(1,srf(isrf)%ivtx(2))
                    p1(2,2) = vtx(2,srf(isrf)%ivtx(2))
                    p1(3,2) = p1(1,2) * DTAN(-0.5D0*dangle)
                    p2(1,2) = p1(1,2)
                    p2(2,2) = p1(2,2)
                    p2(3,2) = p2(1,2) * DTAN(+0.5D0*dangle)
                  ELSE
                    p1(1,1) = obj(iobj)%v(1,obj(iobj)%ipts(1,isid)) 
                    p1(2,1) = obj(iobj)%v(2,obj(iobj)%ipts(1,isid)) 
                    p1(3,1) = DBLE(p1(1,1))*DTAN(DBLE(-0.5*dangle))
                    p2(1,1) = p1(1,1)
                    p2(2,1) = p1(2,1)
                    p2(3,1) = DBLE(p2(1,1))*DTAN(DBLE(+0.5*dangle))

                    p1(1,2) = obj(iobj)%v(1,obj(iobj)%ipts(2,isid)) 
                    p1(2,2) = obj(iobj)%v(2,obj(iobj)%ipts(2,isid)) 
                    p1(3,2) = DBLE(p1(1,2))*DTAN(DBLE(-0.5*dangle))
                    p2(1,2) = p1(1,2)
                    p2(2,2) = p1(2,2)
                    p2(3,2) = DBLE(p2(1,2))*DTAN(DBLE(+0.5*dangle))
                  ENDIF
c...              Rotate vertices:
                  DO i3 = 1, 2
                    x1 = p1(1,i3)
                    z1 = p1(3,i3)
                    p1(1,i3) = DCOS(DBLE(ang)) * x1 - DSIN(DBLE(ang))*z1
                    p1(3,i3) = DSIN(DBLE(ang)) * x1 + DCOS(DBLE(ang))*z1
                    x1 = p2(1,i3)
                    z1 = p2(3,i3)
                    p2(1,i3) = DCOS(DBLE(ang)) * x1 - DSIN(DBLE(ang))*z1
                    p2(3,i3) = DSIN(DBLE(ang)) * x1 + DCOS(DBLE(ang))*z1
                    call transform_vect(mat,p1(1,i3))
                    call transform_vect(mat,p2(1,i3))
                  ENDDO
                  CALL POSITN (SNGL(p1(1,1)),SNGL(p1(2,1)))
                  CALL JOIN   (SNGL(p1(1,2)),SNGL(p1(2,2))) 
                  CALL POSITN (SNGL(p1(1,1)),SNGL(p1(2,1)))
                  CALL JOIN   (SNGL(p2(1,1)),SNGL(p2(2,1))) 
                  CALL POSITN (SNGL(p1(1,2)),SNGL(p1(2,2)))
                  CALL JOIN   (SNGL(p2(1,2)),SNGL(p2(2,2))) 
                  CALL POSITN (SNGL(p2(1,1)),SNGL(p2(2,1)))
                  CALL JOIN   (SNGL(p2(1,2)),SNGL(p2(2,2))) 
                ENDDO
              ENDDO
            ELSEIF (obj(iobj)%gsur(isid).EQ.GT_TD) THEN
              IF (obj(iobj)%nside.NE.0) THEN

c...            Filter:
                count = 0.0
                IF (obj(iobj)%tsur(isid).NE.SP_GRID_BOUNDARY) CYCLE
c                WRITE(0,*) 'GO MAN',iobj,isid
                DO isrf = isrf1, isrf2
                  DO i3 = 1, srf(isrf)%nvtx
                    i4 = i3 + 1
                    IF (i4.GT.srf(isrf)%nvtx) i4 = 1
                    p1(1:3,1) = vtx(1:3,srf(isrf)%ivtx(i3))
                    p2(1:3,1) = vtx(1:3,srf(isrf)%ivtx(i4))
                    call transform_vect(mat,p1(1,1))
                    call transform_vect(mat,p2(1,1))
                    CALL POSITN (SNGL(p1(1,1)),SNGL(p1(2,1)))
                    CALL JOIN   (SNGL(p2(1,1)),SNGL(p2(2,1))) 
                  ENDDO
                ENDDO
              ELSE
                DO i3 = 1, obj(iobj)%npts(isid)
                  i4 = i3 + 1
                  IF (i4.EQ.obj(iobj)%npts(isid)+1) i4 = 1
                  p1(1,1) = obj(iobj)%v(1,obj(iobj)%ipts(i3,isid))
                  p1(2,1) = obj(iobj)%v(2,obj(iobj)%ipts(i3,isid))
                  p1(3,1) = obj(iobj)%v(3,obj(iobj)%ipts(i3,isid))
                  p2(1,1) = obj(iobj)%v(1,obj(iobj)%ipts(i4,isid))
                  p2(2,1) = obj(iobj)%v(2,obj(iobj)%ipts(i4,isid))
                  p2(3,1) = obj(iobj)%v(3,obj(iobj)%ipts(i4,isid))
                  call transform_vect(mat,p1(1,1))
                  call transform_vect(mat,p2(1,1))
                  CALL POSITN (SNGL(p1(1,1)),SNGL(p1(2,1)))
                  CALL JOIN   (SNGL(p2(1,1)),SNGL(p2(2,1))) 
                ENDDO
              ENDIF
            ELSE
              CALL ER('Plot985','Unknown surface geometry type',*98)
            ENDIF
          ENDDO
        ENDDO 
c...    Surface highlighting:
        DO iobj = 1, 0 ! nobj
          CALL LINCOL(ncols+obj(iobj)%colour+2) 
          DO isid = 1, MAX(obj(iobj)%nsur,obj(iobj)%nside)

            IF ((obj(iobj)%tsur(isid).NE.SP_VESSEL_WALL.AND.
     .           obj(iobj)%tsur(isid).NE.SP_GRID_BOUNDARY).OR.
     .          obj(iobj)%esurf(isid).NE.21) CYCLE

            IF (obj(iobj)%nside.NE.0) THEN
              isrf1 = obj(iobj)%iside(isid,1)
              isrf2 = obj(iobj)%iside(isid,2)
            ELSE
              isrf1 = 1
              isrf2 = 1
            ENDIF

            IF     (obj(iobj)%gsur(isid).EQ.GT_TC) THEN
            ELSEIF (obj(iobj)%gsur(isid).EQ.GT_TD) THEN
              IF (obj(iobj)%nside.NE.0) THEN
                DO isrf = isrf1, isrf2
                  DO i3 = 1, srf(isrf)%nvtx
                    i4 = i3 + 1
                    IF (i4.GT.srf(isrf)%nvtx) i4 = 1
                    p1(1:3,1) = vtx(1:3,srf(isrf)%ivtx(i3))
                    p2(1:3,1) = vtx(1:3,srf(isrf)%ivtx(i4))
                    call transform_vect(mat,p1(1,1))
                    call transform_vect(mat,p2(1,1))
                    CALL POSITN (SNGL(p1(1,1)),SNGL(p1(2,1)))
                    CALL JOIN   (SNGL(p2(1,1)),SNGL(p2(2,1))) 
                  ENDDO
                ENDDO
              ELSE
                STOP 'OUT OF DATE'
              ENDIF
            ELSE
              CALL ER('Plot985','Unknown surface geometry type',*98)
            ENDIF
          ENDDO
        ENDDO 
c...    Draw pixel views:
        IF (.TRUE.) THEN
          CALL LINCOL(ncols+55) 
          DO i1 = 1, MIN(500,npixel)
            p1(1,1) = pixel(i1)%global_v1(1)
            p1(2,1) = pixel(i1)%global_v1(2)
            p1(3,1) = pixel(i1)%global_v1(3)
            p2(1,1) = pixel(i1)%global_v2(1)
            p2(2,1) = pixel(i1)%global_v2(2)
            p2(3,1) = pixel(i1)%global_v2(3)
            call transform_vect(mat,p1(1,1))
            call transform_vect(mat,p2(1,1))
            CALL POSITN (SNGL(p1(1,1)),SNGL(p1(2,1)))
            CALL JOIN   (SNGL(p2(1,1)),SNGL(p2(2,1))) 
          ENDDO
        ELSEIF (.TRUE.) THEN
c        ELSEIF (.TRUE..AND.nplot.GT.1) THEN
          CALL LINCOL(ncols+55) 
c          DO i1 = 1, MIN(1,nchord)
c          DO i1 = 1, MIN(500,nchord)
          DO i1 = 1, MIN(500,npixel)
            p1(1,1) = s_chord(i1)%v1(1)
            p1(2,1) = s_chord(i1)%v1(2)
            p1(3,1) = s_chord(i1)%v1(3)
            p2(1,1) = s_chord(i1)%v2(1)
            p2(2,1) = s_chord(i1)%v2(2)
            p2(3,1) = s_chord(i1)%v2(3)
            call transform_vect(mat,p1(1,1))
            call transform_vect(mat,p2(1,1))
            CALL POSITN (SNGL(p1(1,1)),SNGL(p1(2,1)))
            CALL JOIN   (SNGL(p2(1,1)),SNGL(p2(2,1))) 
          ENDDO
        ENDIF
c...    Frame:
        CALL LINCOL(1)
        CALL DrawFrame
      ENDDO


      RETURN
98    WRITE(0,*) 'OBJECT,SIDE,TYPE=',i1,i2,obj(i1)%gsur(i2)
99    STOP
      END

c
c ======================================================================
c
c subroutine: Output985
c
      SUBROUTINE ImagePlot985(iopt,MAXPIXEL,npixel,pixel,image,
     .                        iplot)
      USE mod_out985
      USE mod_out985_variables
      USE mod_interface
      IMPLICIT none

      INTEGER iopt,MAXPIXEL,npixel         ! Put this into _variables... 
      INTEGER, INTENT(IN) :: iplot
      TYPE(type_view) :: pixel(MAXPIXEL)
      REAL*8 image(1100,1100)

      INCLUDE 'params'      
      INCLUDE 'comgra'
      INCLUDE 'colours'
      INCLUDE 'slout'

      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost

      CHARACTER glabel*512,caption*1024
      CHARACTER TITLE*80,TITLE2*80,JOB*72,GRAPH*80,GRAPH1*256
      CHARACTER*36 REF,NVIEW,PLANE,ANLY,TABLE,ZLABS(-2:MAXIZS+1)
      CHARACTER*36 XLAB
      CHARACTER*72 YLAB
      CHARACTER*72 SMOOTH
      CHARACTER*256 cdum1

      INTEGER axis1(3)
      REAL    angle1(3)

      INTEGER CH1

      INTEGER   nxbin,nybin,ipixel,nbin,isrf,isrf1,isrf2,isid,
     .          idet,nplot,iobj,i1,i2,i3,i4,ix,iy,iplot1
      REAL      deltar,deltaz,qmin,qmax,qval,ang1,ang2,ang,
     .          frac1,frac2,xcen,ycen,xnear,ynear,count
      REAL      XXMIN,XXMAX,YYMIN,YYMAX,dangle
      REAL*8    angle,p1(3,6),p2(3,6),x1,x2,z1,mat(3,3)

      INTEGER dat1,dat2
      REAL, ALLOCATABLE :: xdat(:),ydat(:)
      CHARACTER xlabel*256,ylabel*256,tag_x*2,tag_y*2,file*512


c *TEMP*
      INTEGER, ALLOCATABLE :: nv(:)
      REAL   , ALLOCATABLE :: rv(:,:),zv(:,:),cq(:)



c...  PLOT B: 3D LOS integration through the 3D objects
c     -build list of primary chords
c     -build list of secondary chords (reflections)
c     -fast routine for determining the intersection between a line and a surface
c     -ray trace (using connection map and wedge index to speed things up)



      IF (npixel.GT.1) THEN
c...    Image:

        DO idet = 1, opt%ndet

c          CALL FRAME

          nxbin = 0
          nybin = 0
          qmin =  0.0
c          qmin =  HI
          qmax = -HI
          DO ipixel = opt%det_istart(idet), opt%det_iend(idet)
c          DO ipixel = 1, npixel
            ix = pixel(ipixel)%xindex
            iy = pixel(ipixel)%yindex
            nxbin = MAX(nxbin,ix)
            nybin = MAX(nybin,iy)
c            qmin = MIN(qmin,pixel(ipixel)%integral)
            qmax = MAX(qmax,SNGL(pixel(ipixel)%integral(1)))
          ENDDO
c          qmin = MIN(qmin,0.0)

          WRITE(0,*) 'QMAD:',qmin,qmax

          IF (qmin.EQ.qmax) THEN
            WRITE(0,*) 'PROBLEM... CYCLING'
            CYCLE
          ENDIF

          IF (nxbin.GE.nybin) THEN
            map1x = 0.05            
            map2x = map1x + 0.55
            map2y = 0.95
            map1y = map2y - 0.55 * REAL(nybin) / REAL(nxbin)
          ELSE
            map1x = 0.05            
            map2x = map1x + 0.55 * REAL(nxbin) / REAL(nybin)
            map1y = 0.40 
            map2y = map1y + 0.55 
          ENDIF 

          CALL PSPACE(map1x,map2x,map1y,map2y)
          CALL MAP   (0.0,1.0,1.0,0.0)

          ALLOCATE(nv(1))
          ALLOCATE(rv(4,1))  
          ALLOCATE(zv(4,1))
          ALLOCATE(cq(1))

          IF (.NOT..TRUE.) THEN
            DO frac1 = 0.0, 1.0*0.99999, 0.01                  ! Needs more work, grouping regons of common colour. 
              CALL SetCol255_04(2,frac1+0.005,0.0,1.0)         ! Otherwise, not much savings... 
              CALL FILCOL(255)                         
              CALL LINCOL(255)                         
              DO ipixel = opt%det_istart(idet), opt%det_iend(idet)
c              DO ipixel = 1, npixel
                cq(1) = SNGL(pixel(ipixel)%integral(1))
                frac2 = (cq(1) - qmin) / (qmax - qmin)
                IF (frac2.LT.frac1.OR.frac2.GE.frac1+0.01) CYCLE

                ix = pixel(ipixel)%xindex
                iy = pixel(ipixel)%yindex
                deltar = 1.0 / REAL(nxbin)
                deltaz = 1.0 / REAL(nybin)
                i1 = 1
                rv(1,i1) = (ix - 1) * deltar 
                zv(1,i1) = (iy - 1) * deltaz
                rv(2,i1) = (ix - 1) * deltar 
                zv(2,i1) = (iy    ) * deltaz
                rv(3,i1) = (ix    ) * deltar 
                zv(3,i1) = (iy    ) * deltaz
                rv(4,i1) = (ix    ) * deltar 
                zv(4,i1) = (iy - 1) * deltaz
                CALL PTPLOT(rv(1,i1),zv(1,i1),1,4,1)
              ENDDO
            ENDDO
          ELSE
c            DO ipixel = 1, npixel
            WRITE(file,'(1024X)')          
            WRITE(file,'(A,I0.2)') 
     .        'idl.'//TRIM(opt%fmap)//'_los_',idet  
            CALL inOpenInterface(TRIM(file),ITF_WRITE)
            DO ipixel = opt%det_istart(idet), opt%det_iend(idet)
              ix = pixel(ipixel)%xindex
              iy = pixel(ipixel)%yindex
              deltar = 1.0 / REAL(nxbin)
              deltaz = 1.0 / REAL(nybin)
              i1 = 1
              rv(1,i1) = (ix - 1) * deltar 
              zv(1,i1) = (iy - 1) * deltaz
              rv(2,i1) = (ix - 1) * deltar 
              zv(2,i1) = (iy    ) * deltaz
              rv(3,i1) = (ix    ) * deltar 
              zv(3,i1) = (iy    ) * deltaz
              rv(4,i1) = (ix    ) * deltar 
              zv(4,i1) = (iy - 1) * deltaz
              cq(1) = SNGL(pixel(ipixel)%integral(1))
              CALL SetCol255_04(2,cq(i1),qmin,qmax)    ! Should really reorganize this loop, to draw all pixels of a 
              CALL FILCOL(255)                         ! particular shade, from a limited set of contours, to 
              CALL LINCOL(255)                         ! try and shrink the size of the .ps file... 
              CALL PTPLOT(rv(1,i1),zv(1,i1),1,4,1)
              CALL inPutData(ix,'i','none')
              CALL inPutData(iy,'j','none')
              CALL inPutData(cq(i1),'data','ph m-2 s-1')
              CALL inPutData(pixel(ipixel)%global_v1(1),'x1','m')
              CALL inPutData(pixel(ipixel)%global_v1(2),'y1','m')
              CALL inPutData(pixel(ipixel)%global_v1(3),'z1','m')
              CALL inPutData(pixel(ipixel)%global_v2(1),'x2','m')
              CALL inPutData(pixel(ipixel)%global_v2(2),'y2','m')
              CALL inPutData(pixel(ipixel)%global_v2(3),'z2','m')
            ENDDO
            CALL inCloseInterface
          ENDIF
c...      Clear arrays:
          DEALLOCATE(nv)
          DEALLOCATE(rv)
          DEALLOCATE(zv)
          DEALLOCATE(cq)
          CALL DrawColourScale(1,2,qmin,qmax,'none')
c...      Frame the plot:
          CALL DrawFrame

c...      Inversion mesh coverage:
          IF (.TRUE.) THEN
            cxmin =  HI
            cxmax = -HI
            cymin =  HI
            cymax = -HI
            qmin = 0.0
            qmax = 0.0
            file = 'output.trc'
            WRITE(0,*) 'DUMP INVERSION COVERAGE:',file(1:LEN_TRIM(file))
            CALL inOpenInterface(file,ITF_WRITE)   ! TRIM(file) would not work, compiler bug...
            DO iobj = 1, nobj
              IF (obj(iobj)%type.NE.OP_INTEGRATION_VOLUME) CYCLE
              IF (obj(iobj)%type.NE.GT_TC) CYCLE
              CALL inPutData(obj(iobj)%nver,'npts','none')            
              DO i2 = 1, obj(iobj)%nver
                WRITE(tag_x,'(A,I1)') 'x',i2
                WRITE(tag_y,'(A,I1)') 'y',i2
                CALL inPutData(obj(iobj)%v(1,i2),tag_x,'m')            
                CALL inPutData(obj(iobj)%v(2,i2),tag_y,'m')            
                cxmin = MIN(cxmin,obj(iobj)%v(1,i2))
                cxmax = MAX(cxmax,obj(iobj)%v(1,i2))
                cymin = MIN(cymin,obj(iobj)%v(2,i2))
                cymax = MAX(cymax,obj(iobj)%v(2,i2))
              ENDDO
              CALL inPutData(obj(iobj)%path,'path','m')            
              qmax = MAX(qmax,obj(iobj)%path)
            ENDDO     
            CALL inCloseInterface
          ENDIF
          WRITE(0,*) 'PROFILE QMIN,QMAX :',qmin,qmax
          IF (qmax.GT.0.0) THEN
            map1x = 0.05
            map1y = 0.02
            IF (cxmax-cxmin.GT.cymax-cymin) THEN
              map2x = map1x + 0.35
              map2y = map1y + 0.35 * (cymax - cymin) / (cxmax - cxmin)
            ELSE
              map2x = map1x + 0.35 * (cxmax - cxmin) / (cymax - cymin)
              map2y = map1y + 0.35 
            ENDIF
            CALL PSPACE(map1x,map2x,map1y,map2y)
            CALL MAP   (cxmin,cxmax,cymin,cymax)
            DO iobj = 1, nobj
              IF (obj(iobj)%type.NE.OP_INTEGRATION_VOLUME) CYCLE
              IF     (obj(iobj)%path.GT.0.01*qmax) THEN
                CALL SetCol255_04(2,obj(iobj)%path,qmin,qmax)
              ELSEIF (obj(iobj)%path.GT.0.0) THEN
                CALL SetCol255_04(1,qmax*0.9      ,qmin,qmax)
              ELSE
                CALL SetCol255_04(1,qmax*0.6      ,qmin,qmax)
              ENDIF
              CALL FILCOL(255)
              CALL LINCOL(255) 
              CALL PTPLOT(REAL(obj(iobj)%v(1,1:4)),
     .                    REAL(obj(iobj)%v(2,1:4)),
     .                    1,obj(iobj)%nver,1)
c             WRITE(0,*) 'plot:',obj(iobj)%v(1,1:4),
c    .                           obj(iobj)%v(2,1:4)
            ENDDO
c...      Annotate:
c          CALL LinCol(ncols+1)
c          CALL CTRMAG(12)
c          WRITE(caption,'(A,3I4)') 'NX,NY,BIN:',opt%nxbin,opt%nybin,nbin
c          CALL PLOTST(0.02,0.02,caption(1:LEN_TRIM(caption)))
c...      Frame:
c          CALL DrawGrid(-95)
c          CALL Supimp('PARTIAL')
            Call DrawFrame
          ENDIF


c...      Line plot:
          IF (.FALSE.) THEN
            slopt2 = 1
            iopt_ghost = 1  ! 2
            plottype(1) = 2
            plottype(2) = 3
            map1x = 0.65             
            map2x = map1x + 0.55
            map1y = 0.77 
            map2y = map1y + 0.15
c...        Assign data:
            dat1 = 35600+1
            dat2 = 35600+200
            ALLOCATE(xdat(dat2-dat1+1))
            ALLOCATE(ydat(dat2-dat1+1))
            xlabel = 'pixel    '
            ylabel = 'signal   '
            DO i1 = dat1, dat2
              xdat(i1-dat1+1) = REAL(i1)
            ENDDO
            ydat(1:dat2-dat1+1) = SNGL(pixel(dat1:dat2)%integral(1))
            IF (idet.EQ.1) THEN
              WRITE(6,*) '*PIXEL DATA'
              WRITE(6,*) '*'
              WRITE(6,*) dat2-dat1+1
              DO i1 = 1, dat2-dat1+1
                WRITE(6,*) REAL(i1),xdat(i1),ydat(i1)
              ENDDO
            ENDIF
            cxmin = 1.0
            cxmax = REAL(dat2-dat1+1)
            cymin =  HI
            cymax = -HI
            DO i1 = 1, dat2-dat1+1
              cymin = MIN(cymin,ydat(i1))
              cymax = MAX(cymax,ydat(i1))
            ENDDO
            CALL GRTSET_TRIM(' ',' ',' ',' ',' ',      ! TITLE,REF,nVIEW,PLANE,glabel,
     .                       cxmin,cxmax,cymin,cymax,
     .                       ' ',xlabel,ylabel,        ! TABLE,XLAB,YLAB,
     .                       0,' ',0,' ',1)            ! 0,smooth,0,ANLY,1)
            CALL GRTRAC(xdat,ydat,dat2-dat1+1,'ref ','LINE',1)        
            DEALLOCATE(xdat)
            DEALLOCATE(ydat)
            CALL DrawFrame

          ENDIF

        ENDDO
      ENDIF


      IF (.TRUE..AND.opt%img_opt.NE.0) THEN
c...    Plot loaded camera image (if there is one):

c        IF (opt%img_nxbin.GE.100) THEN
c          nbin = opt%img_nxbin / 50          ! Should also check opt%nybin
c          nbin = 1
c        ELSE
          nbin = 1
c        ENDIF

c...    Plot image:
c        nbin = 1 !6  ! Binning (NBIN=1 is no binning)
        IF (nbin.GT.1.AND.
     .      (nbin.GE.opt%img_nxbin.OR.
     .       nbin.GE.opt%img_nybin))
     .    CALL ER('985','Bin request greater than image resolution',
     .            *99)

        nxbin = opt%img_nxbin / nbin
        nybin = opt%img_nybin / nbin

        WRITE(0,*) ' *** NXBIN',nxbin,nybin

        qmin =  HI
        qmax = -HI
        DO ix = 1, nxbin
          DO iy = 1, nybin
            qval = opt%img_image(ix,iy)
c            qval = 0.0
c            DO i2 = 0, nbin-1
c              DO i3 = 0, nbin-1
c                qval = qval + SNGL(opt%img_image(nbin*(ix-1)+1+i2,
c     .                                           nbin*(iy-1)+1+i3))
c              ENDDO
c            ENDDO           
            qmin = MIN(qmin,qval)  ! Need to average here, and below?  
            qmax = MAX(qmax,qval)  ! Need to average here, and below?  
          ENDDO
        ENDDO

        WRITE(0,*) 'QVAL MASK IMAGE:',qmin,qmax,nxbin,nybin
c        STOP 'sdfsd'


        IF (nxbin.GE.nybin) THEN
          map1x = 0.65            
          map2x = map1x + 0.55
          map2y = 0.95
          map1y = map2y - 0.55 * REAL(nybin) / REAL(nxbin)
        ELSE
          map1x = 0.65            
          map2x = map1x + 0.55 * REAL(nxbin) / REAL(nybin)
          map1y = 0.40
          map2y = map1y + 0.55 
        ENDIF 

        CALL PSPACE(map1x,map2x,map1y,map2y)
        CALL MAP   (0.0,1.0,1.0,0.0)
        ALLOCATE(nv(1))
        ALLOCATE(rv(4,1))  
        ALLOCATE(zv(4,1))
        ALLOCATE(cq(1))
        deltar = 1.0 / REAL(nxbin)
        deltaz = 1.0 / REAL(nybin)
        DO ix = 1, nxbin
          DO iy = 1, nybin
            i1 = 1
            rv(1,i1) = REAL(ix - 1) * deltar 
            zv(1,i1) = REAL(iy - 1) * deltaz
            rv(2,i1) = REAL(ix - 1) * deltar 
            zv(2,i1) = REAL(iy    ) * deltaz
            rv(3,i1) = REAL(ix    ) * deltar 
            zv(3,i1) = REAL(iy    ) * deltaz
            rv(4,i1) = REAL(ix    ) * deltar 
            zv(4,i1) = REAL(iy - 1) * deltaz
            cq(i1) = opt%img_image(ix,iy)
c            cq(i1) = 0.0
c            DO i2 = 0, nbin-1
c              DO i3 = 0, nbin-1
c                cq(i1) = cq(i1) + SNGL(opt%img_image(nbin*(ix-1)+1+i2,
c     .                                               nbin*(iy-1)+1+i3))
c              ENDDO
c            ENDDO
            CALL SetCol255_04(2,cq(i1),qmin,qmax)   ! See above note on reducing size of .ps file...
            CALL FILCOL(255)
            CALL LINCOL(255) 
            CALL PTPLOT(rv(1,i1),zv(1,i1),1,4,1)
          ENDDO
        ENDDO
c...    Clear arrays:
        IF (ALLOCATED(nv)) DEALLOCATE(nv)
        IF (ALLOCATED(rv)) DEALLOCATE(rv)
        IF (ALLOCATED(zv)) DEALLOCATE(zv)
        IF (ALLOCATED(cq)) DEALLOCATE(cq)
c...    Annotate:
        CALL LinCol(ncols+1)
        CALL CTRMAG(12)
        WRITE(caption,'(A,2I5)') 'XBIN,YBIN:',nxbin,nybin
        CALL PLOTST(0.02,0.03,caption(1:LEN_TRIM(caption)))
        CALL DrawFrame
c...    Clear memory:
c        IF (ALLOCATED(opt%img_image)) DEALLOCATE(opt%img_image)     
      ENDIF


      IF (.TRUE.) THEN
c...    Spectra:
        slopt2 = 1
        iopt_ghost = 1
        DO idet = 1, opt%ndet
          WRITE(0,*) 'PLOTTING LINE SHAPES',
     .        opt%det_istart(idet),
     .        opt%det_iend  (idet)
          CALL FRAME
          CALL PlotLineShapes(opt%det_istart(idet),
     .                        opt%det_iend  (idet),
     .                        MAXPIXEL,npixel,pixel)
        ENDDO
        WRITE(0,*) 'DONE'
      ENDIF


      RETURN
98    WRITE(0,*) 'OBJECT,SIDE,TYPE=',i1,i2,obj(i1)%gsur(i2)
99    STOP
      END


c
c ======================================================================
c
c subroutine: Output985
c
      SUBROUTINE Output985(iopt,MAXPIXEL,npixel,pixel,image)
      USE mod_out985
      USE mod_out985_variables
      USE mod_interface
      IMPLICIT none

      INTEGER iopt,MAXPIXEL,npixel         ! Put this into _variables... 
      TYPE(type_view) :: pixel(MAXPIXEL)
      REAL*8 image(1100,1100)

      INTEGER   iplot,option
      LOGICAL   debug
      CHARACTER buffer*1024

      debug = .FALSE.

      IF (debug) THEN
        WRITE(0,*) 'PLOT LIST:',opt%nplots
        DO iplot = 1, opt%nplots
          WRITE(0,*) TRIM(opt%plots(iplot))
        ENDDO
      ENDIF

      DO iplot = 1, opt%nplots
        WRITE(buffer,'(1024X)')
c        WRITE(0,*) 'PLOTS:',opt%plots(iplot)
        READ(opt%plots(iplot),*) buffer
        IF (buffer(1:4).NE.'plot') CYCLE
        READ(opt%plots(iplot),*) opt%plots(iplot),option
        SELECTCASE (option)
          CASE (000)
          CASE (001)
c....       Wireframe (all plots currently...):
            CALL WireframePlot985
     .             (iopt,MAXPIXEL,npixel,pixel,image,iplot)
            CALL FRAME
          CASE (002)
c....       Solid:
c            WRITE(0,*) 'DRAWING SOLID 3D PLOT'
            CALL SolidPlot985(opt,nobj,obj,iplot)
            CALL FRAME
          CASE (003)
c....       Image plot, and line shapes for the moment:
            CALL ImagePlot985
     .             (iopt,MAXPIXEL,npixel,pixel,image,iplot)
            CALL FRAME
          CASE DEFAULT
            CALL WN('Output985','Unrecognised plot option')
            WRITE(0,*) '  BUFFER= "',buffer(1:LEN_TRIM(buffer)),'"'
        ENDSELECT
      ENDDO      

c      CALL FRAME

      RETURN
 99   STOP
      END
