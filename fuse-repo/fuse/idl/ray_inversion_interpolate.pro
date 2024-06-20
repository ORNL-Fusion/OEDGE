;
; New
; ===
;
; helium
; ------
; a=ray(shot=17107,channel=1,frame=12,ext='_diamonds.cgm',region=2,/full)
; a=ray(shot=18843,channel=2,frame=13,ext='_500x500_5mm_diamonds.cgm',region=3,/full) 
;
;
;
; OLD...
; 
; RGB
; ---
; a=ray(shot=18360,channel=4,frame=56,ext='cgm')
; save,a,filename='RGB_18360_56_4_cgm_a.sav',/compress
;
; Helium
; ------
;
; a=ray(shot=17107,channel=1,frame=12,ext='cgm')
; a=ray(shot=[17107,17108,17107],channel=[2,2,1],frame=12,ext='_diamonds.cgm')
;
; a=ray(shot=[18845,18845],channel=[1,2],frame=[18,13],ext='_500x500_5mm_diamonds.cgm')
;
; a=ray(shot=15021,channel=2,frame=14,ext='cgm')
; a=ray(shot=15021,channel=2,frame=14,ext='_test.osm')
; a=ray(shot=[15021,15021],channel=2,frame=14,ext=['_500x500_5mm_diamonds.cgm','_500x500_5mm_diamonds.osm'])
;
; Pellets
; -------
; a=ray(shot=[17985,17985],channel=[1,2],frame=15,ext='sav',/test,row=600,yrange=[0.0,0.05],xrange=[500,1000])
;
;
; a=ray(shot=[17107,17107],channel=[2,1],frame=12,col=45,xrange=[50,170],yrange=[4,10],ps='668_728_in')
; a=ray(shot=[17107,17107],channel=[2,1],frame=12,row=220,xrange=[10,220],yrange=[2,8],ps='668_728_out')
;
; a=ray(shot=[17107,17108],channel=[1,2],frame=12,col=170,xrange=[50,170],yrange=[0.1,0.5],fudge=[1.0,0.8],ps='728_707_in')
; a=ray(shot=[17107,17108],channel=[1,2],frame=12,row=230,xrange=[10,230],yrange=[0.1,0.3],fudge=[1.0,0.8],ps='728_707_out')
; 
;
;

; 
; a=ray(shot=[17107,17108],channel=[1,2],frame=12,row=100,yr=[0.,0.5],fudge=[1.0,0.8])
;
;


;ray,shot=15169,lines=-1,row=250,scale=1,frame=13
;  some hardcoding: fname1=fname2,image2[*,*]=1.0, see "Midplane"
;
;
; ray,shot=15018,lines=2,row=17,scale=4,frame=10


; need to line up the peaks?  or at least to a careful/automated check that the peaks line up properly 
;
; ======================================================================
;
FUNCTION MakeImageBig,image,scale

  dim1 = N_ELEMENTS(image[*,0])
  dim2 = N_ELEMENTS(image[0,*])

  print,dim1,dim2

  image_big = MAKE_ARRAY(scale*dim1,scale*dim2,/FLOAT,value=0.0)
  FOR idim2 = 0, dim2-1 DO BEGIN
    FOR idim1 = 0, dim1-1 DO BEGIN
      image_big[scale*idim1:scale*idim1+scale-1,scale*idim2:scale*idim2+scale-1] = image[idim1,idim2]
    ENDFOR
  ENDFOR

  RETURN,image_big

END
;
; ======================================================================
; ======================================================================
;
FUNCTION ProcessImage, inv, chisq_limit, fit_sample=fit_sample, fit_cutoff=fit_cutoff, progress=progress, plots=plots

  colors = ['Black','Red','Blue','Orange','Purple', 'Hotpink', 'Green']

  IF (NOT KEYWORD_SET(fit_sample)) THEN fit_sample = 20

  nseg  = inv.contour.n
  nxpts = 250
  nypts = MAX([1,fit_sample/nseg])
  IF (KEYWORD_SET(fit_cutoff)) THEN threshold = fit_cutoff ELSE threshold = 0.3

  print,'nseg      = ',nseg
  print,'fit_sample= ',fit_sample
  print,'nypts     = ',nypts

  res       = MAKE_ARRAY(nxpts,(nseg-1)*nypts,/FLOAT)
  res_fit   = MAKE_ARRAY(nxpts,(nseg-1)*nypts,/FLOAT)
  res_chisq = MAKE_ARRAY((nseg-1)*nypts,/FLOAT)
  res_shift = MAKE_ARRAY((nseg-1)*nypts,/LONG)
  res_i     = MAKE_ARRAY((nseg-1)*nypts,/LONG)
  res_frac  = MAKE_ARRAY((nseg-1)*nypts,/FLOAT)
  res_dist  = MAKE_ARRAY((nseg-1)*nypts,/FLOAT)
  res_ix    = MAKE_ARRAY((nseg-1)*nypts,/FLOAT)
  res_iy    = MAKE_ARRAY((nseg-1)*nypts,/FLOAT)

  ix = inv.contour.ix
  iy = inv.contour.iy

  icolour = 0
  IF (KEYWORD_SET(plots)) THEN BEGIN
    window,1
    !P.BACKGROUND = TrueColor('White')
    plot,[0.0],[0.0],/nodata,xrange=[0.25,0.75],yrange=[0.0,5.0E+21],color=Truecolor('Black')
  ENDIF

  dist  = 0.0
  xpos1 = inv.contour.x[0]
  ypos1 = inv.contour.y[0]

  FOR i = 1, N_ELEMENTS(ix)-1 DO BEGIN

    a = FLOAT([ix[i-1],iy[i-1]])
    b = FLOAT([ix[i  ],iy[i  ]])
    c = [0.0,0.0]
    d = [0.0,0.0]

    IF (KEYWORD_SET(plots)) THEN BEGIN
      wset,0
      plots,LONG([a[0],b[0]]),LONG([a[1],b[1]]),color=Truecolor('White'),thick=1,/device
    ENDIF

    CASE ABS(inv.region) OF
      1: e = b - a
      2: e = a - b
      3: e = b - a
      4: e = a - b
    ENDCASE

;    IF (a[1] LT 0.0) THEN e = a - b  ; Crap.  Quick fix...
;    IF (a[1] GE 0.0) THEN e = b - a
;    e = a - b

    len_e = SQRT(e[0]^2 + e[1]^2)
    theta = ACOS(e[0] / len_e) 

;    print,theta*180/3.1415

    length = 100.0
    d[0] = length * SIN(theta)
    d[1] = length * COS(theta)

    FOR t = 0.0, 0.999, 1.0/FLOAT(nypts) DO BEGIN
      c = t * (b - a) + a
      v1 = -d + c  ; d
      v2 =  d + c  
      e     = FINDGEN(nxpts) / FLOAT(nxpts-1)
      edata = FLTARR (nxpts) 

      xpos2 = inv.contour.x[i-1] + t * (inv.contour.x[i] - inv.contour.x[i-1])
      ypos2 = inv.contour.y[i-1] + t * (inv.contour.y[i] - inv.contour.y[i-1])
      dist = dist + SQRT((xpos2 - xpos1)^2 + (ypos2 - ypos1)^2)
      xpos1 = xpos2
      ypos1 = ypos2

      ixpos = FLOAT(ix[i-1]) + t * FLOAT(ix[i] - ix[i-1])
      iypos = FLOAT(ix[i-1]) + t * FLOAT(ix[i] - ix[i-1])

;      print,'THETA:',theta,c,v1,v2

      n  = N_ELEMENTS(inv.x) ; inv.grid.nx
      nx = N_ELEMENTS(inv.x) 
      ny = N_ELEMENTS(inv.y) 
      m = FINDGEN(n)+1.0
      FOR j = 0, N_ELEMENTS(e)-1 DO BEGIN
        k = (v2-v1) * e[j] + v1
        IF (LONG(k[0]) GE 0 AND LONG(k[0]) LE nx-1 AND    $
            LONG(k[1]) GE 0 AND LONG(k[1]) LE ny-1) THEN  $
          edata[j] = INTERPOL(inv.data[*,LONG(k[1])],m,k[0])
      ENDFOR

      IF (KEYWORD_SET(plots)) THEN BEGIN
        wset,0
        plots,LONG([c[0],v1[0]]),LONG([c[1],v1[1]]),color=TrueColor('White'),thick=1,/device
        plots,LONG([c[0],v2[0]]),LONG([c[1],v2[1]]),color=TrueColor('White'),thick=1,/device
      ENDIF

      icolour = icolour + 1

;      IF (icolour LT 4) THEN CONTINUE
 
      CASE 2 OF
        1: result = MAX(edata,imax)
        2: BEGIN
          i1 = WHERE(edata GE threshold*MAX(edata)) ;  AND edata LE 0.8*MAX(edata))
          xdata = e    [i1]
          ydata = edata[i1]
          ymax = MAX(ydata,imax)

          IF (icolour EQ 1) THEN BEGIN
            A_param = [1.0,10000.0,xdata[imax]-0.10,xdata[imax],0.01]  ; Resonable guesses...
          ENDIF ELSE BEGIN
            A_param[0] = 1.0
            A_param[2] = xdata[imax]-0.05
            A_param[3] = xdata[imax]
          ENDELSE

          weights = MAKE_ARRAY(N_ELEMENTS(xdata),/FLOAT,value=1.0)                
;          weights = (1.0 / (ydata/ymax))^2

          chisq = 1.0E+06
          FOR u = 0.0, 1.00001, 0.2 DO BEGIN
            A_param1 = A_param
            A_param1[2] = A_param[2] + u * 0.1
 
            yfit1 = MPCURVEFIT(xdata, ydata/ymax, weights, A_param1, sigma,   $
                               CHISQ=chisq1,FUNCTION_NAME='gfunct_new',/NODERIVATIVE,/QUIET)          

;             print,'t:',t,chisq,chisq1,u
            IF (chisq1 LT chisq) THEN BEGIN
              if (KEYWORD_SET(progress) OR   $
                  (u GT 0.0 AND t GT 0.01 AND chisq - chisq1 GT 0.0001) ) THEN begin
                IF ( KEYWORD_SET(progress) AND u EQ 0.0) THEN yfit = yfit1
                IF (KEYWORD_SET(plots)) THEN BEGIN
                  wset,1
                  plot ,xdata,ydata,xstyle=1,ystyle=1,yrange=[0.0,1.1*MAX(ydata)],psym=6,color=Truecolor('Black')
                  oplot,xdata,yfit *ymax,color=Truecolor(colors[2])
                  oplot,xdata,yfit1*ymax,color=Truecolor(colors[3])
                  XYOUTS, 0.9 * (MIN(xdata)) + 0.1 * (MAX(xdata)),  $
                          0.95* ymax, STRTRIM(STRING(chisq1),2),  $
                          color=Truecolor('Black'),CHARSIZE=2.0
                ENDIF
;               stop
              endif
              A_fit = A_param1
              yfit  = yfit1
              chisq = chisq1
            ENDIF
          ENDFOR
    

          IF (chisq GT chisq_limit) THEN BEGIN
;           The fit is poor so use peak from the raw data:
            result = MAX(ydata,imax)

;            print,'chisq:',chisq

            IF (KEYWORD_SET(plots)) THEN BEGIN
              wset,1
              plot,xdata,ydata,ystyle=1,yrange=[0.0,1.1*MAX(ydata)],psym=6
              oplot,xdata,yfit*ymax,color=Truecolor(colors[4])
            ENDIF
;            stop

          ENDIF ELSE BEGIN 
            result = MAX(yfit,imax)
          ENDELSE
          result = MIN(ABS(e-xdata[imax]),imax)

          IF (icolour EQ 8) THEN BEGIN
            IF (KEYWORD_SET(plots)) THEN BEGIN
              wset,1
              oplot,e,edata,color=Truecolor(colors[2])
              oplot,xdata,yfit*ymax,color=Truecolor(colors[4])
              oplot,xdata,ydata,psym=6
;              oplot,e[imax],edata[imax],psym=5,color=2
;              stop
            ENDIF
          ENDIF

          END
      ENDCASE

      eshift = 0.5 - e[imax]

      shift = nxpts/2 - imax
      FOR i1 = 0, nxpts-1 DO BEGIN
        IF (i1-shift GE 0 AND i1-shift LE nxpts-1) THEN BEGIN
          res    [i1,icolour-1] = edata[i1-shift]

          gfunct_new, e, A_fit, result
          res_fit[i1,icolour-1] = result[i1-shift] * ymax

          res_shift[icolour-1] = shift
          res_chisq[icolour-1] = chisq
          res_i    [icolour-1] = i
          res_frac [icolour-1] = t

          res_dist[icolour-1] = dist
          res_ix  [icolour-1] = ixpos
          res_iy  [icolour-1] = iypos
        ENDIF
      ENDFOR

      IF (KEYWORD_SET(plots)) THEN BEGIN
        wset,1
;        oplot,e+eshift,edata,color=icolour
;        oplot,e,res[*,icolour-1],color=icolour
      ENDIF

;      result = MAX(res[*,icolour-1],imax)

    ENDFOR

  ENDFOR

  fit = {  $
    nx          : nxpts                   ,  $  ;
    ny          : nypts * (nseg-1)        ,  $  ; 
    segment     : res_i                   ,  $  ;
    dist        : res_dist                ,  $  ;
    ix          : res_ix                  ,  $  ;
    iy          : res_iy                  ,  $  ;
    peak        : res_fit[nxpts/2,*]      ,  $  ;
    distance    : res_frac                ,  $  ;
    chisq_limit : chisq_limit             ,  $  ; 
    chisq       : res_chisq               ,  $  ;
    data        : res_fit                 ,  $  ;
    raw         : res                     ,  $  ;
    threshold   : threshold               ,  $  ;
    shift       : res_shift               ,  $  ;
    quality     : res_chisq LT chisq_limit }

  result = CREATE_STRUCT(inv, 'fit', fit)

  RETURN, result

END  
;
; ======================================================================
; ======================================================================
;
FUNCTION ContourImage, inv, sp_peak=sp_peak, spt=spt, xpt=xpt, plots=plots

  nseg = 10

; If overplotting 2 inversions then need to make sure the plots are 
; on the same spatial map, although this really needs to be done when setting
; up the inversion mesh:
  nx = N_ELEMENTS(inv.x)
  ny = N_ELEMENTS(inv.y)
  x = [inv.x[0],inv.x[nx-1]]
  y = [inv.y[0],inv.y[ny-1]]

;  print,x,y

  IF (KEYWORD_SET(plots)) THEN BEGIN
    WINDOW,0,XSIZE=nx,YSIZE=ny
    DEVICE, DECOMPOSED=0
    LOADCT,5
    TVSCL,inv.data

    !P.BACKGROUND = TrueColor('White')
    PLOT,[0],[0],xrange =x    ,yrange =y,      $
                 xmargin=[0,0],ymargin=[0,0],  $
	         xstyle =5    ,ystyle =5    ,/NOERASE,/NODATA
  ENDIF

; Hardcoding hell...
  CASE inv.machine OF
    'MAST': BEGIN
      CASE inv.shot OF
        24860: BEGIN
          CASE ABS(inv.region) OF
              2: BEGIN  ; Lower inner (24860...)
                x_xpt =  0.673424  ; from m-det-0000c
                y_xpt = -1.175227
                x_spt =  0.280668 
                y_spt = -1.429251 
                END
            ENDCASE
            END 
        ELSE: BEGIN 
          PRINT,''
          PRINT,'WARNING ray_ContourImage: Using default strike-point position'
          PRINT,''
          CASE ABS(inv.region) OF
            1: BEGIN  ; Upper inner (18843, 259 ms)
              x_xpt =  0.587          ; EFK_XPOINT in xpad
              y_xpt =  1.171 - 0.03 
              x_spt =  0.280
              y_spt =  1.376 - 0.04 
              END
            2: BEGIN  ; Lower inner (17107, 239 ms)
              x_xpt =  0.658      
              y_xpt = -1.193 
              x_spt =  0.280
              y_spt = -1.380 + 0.02
              END
            3: BEGIN  ; Upper outer (18843, 259ms)
              x_xpt =  0.587      
              y_xpt =  1.171 - 0.060
              x_spt =  1.084 + 0.040
              y_spt =  1.825
              END
            4: BEGIN  ; Lower outer (17107, 239 ms)
              x_xpt =  0.658         ;  0.58
              y_xpt = -1.193         ; -1.10
              x_spt =  1.145 - 0.02  ;  1.00
              y_spt = -1.825
              END
            ELSE: BEGIN
              PRINT,'ERROR ContourImage: Unrecognized region'
              PRINT,'  REGION=',inv.region
              STOP
              END
          ENDCASE
          END
      ENDCASE
      END
  ENDCASE

  ; Search for peak emission and set that as the strike-point location:
  IF (KEYWORD_SET(sp_peak)) THEN BEGIN
    print,'Peak over-ride!'
    SWITCH inv.region OF
      1:
      2: BEGIN
        a = MIN(ABS(inv.x-x_spt) ,ixmax)
        a = MAX(inv.data[ixmax,*],iymax)
        x_spt = inv.x[ixmax]
        y_spt = inv.y[iymax]
        BREAK
        END
      3:
      4: BEGIN
        a = MIN(ABS(inv.y-y_spt) ,iymax)
        a = MAX(inv.data[*,iymax],ixmax)
        x_spt = inv.x[ixmax]
        y_spt = inv.y[iymax]
        BREAK
        END
      ELSE: BREAK
    ENDSWITCH
;    ixmax = -1
;    iymax = -1
;    maxval = -1.0
;    FOR iy = 0, ny-1 DO BEGIN
;      val = MAX(inv.data[*,iy],ix)
;      IF (val GT maxval) THEN BEGIN
;        ixmax = ix
;        iymax = iy
;        maxval = val
;      ENDIF
;    ENDFOR
;    x_spt = inv.x[ixmax]
;    y_spt = inv.y[iymax]
  ENDIF

  IF (KEYWORD_SET(xpt)) THEN BEGIN
    x_xpt = xpt[0]
    y_xpt = xpt[1]
  ENDIF
  IF (KEYWORD_SET(spt)) THEN BEGIN
    x_spt = spt[0]
    y_spt = spt[1]
  ENDIF

  print,'STRIKE POINT= ',x_spt,y_spt
  print,'X-POINT     = ',x_xpt,y_xpt

  CASE 0 OF
    0: BEGIN
      FOR y1 = y_spt, 0.95*(y_xpt-y_spt)+y_spt, (y_xpt-y_spt)/FLOAT(nseg) DO BEGIN
        t  = (y1 - y_spt) / (y_xpt - y_spt)
        x1 = (1.0 - t) * x_spt + t * x_xpt
        a  = MIN(ABS(inv.x-x1),ix1)
        a  = MIN(ABS(inv.y-y1),iy1)
        IF (y1 EQ y_spt) THEN BEGIN
          ix = ix1
          iy = iy1
          x  = x1
          y  = y1
        ENDIF ELSE BEGIN
          ix = [ix,ix1]
          iy = [iy,iy1]
          x  = [x  ,x1]
          y  = [y  ,y1]
        ENDELSE
        IF (KEYWORD_SET(plots)) THEN BEGIN 
          wset,0
          plots,ix1,iy1,thick=2,psym=6,/device,COLOR=TrueColor('White')
        ENDIF
      ENDFOR
      END
 
;    1: BEGIN
;      i = 0
;      a = MIN(ABS(inv.grid.x-x_xpt),ix1)
;      ix2 = inv.grid.nx-1
;      print,'ix1,2=',ix1,ix2
;      FOR y = y_spt, 0.95*(y_xpt-y_spt)+y_spt, (y_xpt-y_spt)/FLOAT(nseg) DO BEGIN
;        a = MIN(ABS(inv.grid.y-y),imin)
;        a = MAX(inv.grid.data[ix1:ix2,imin],imax)
;        imax = imax + ix1 - 1
;        IF (y EQ y_spt) THEN BEGIN
;          ix = imax
;          iy = imin
;        ENDIF ELSE BEGIN
;          ix = [ix,imax]
;          iy = [iy,imin]
;          i = N_ELEMENTS(ix)
;          plots,ix[i-2:i-1],iy[i-2:i-1],color=3,thick=2,psym=6,/device
;        ENDELSE
;      ENDFOR
;      END
;    2: BEGIN
;      i = 0
;      a = MIN(ABS(inv.grid.x-x_xpt),ix1)
;      ix2 = inv.grid.nx-1
;      print,'ix1,2=',ix1,ix2
;
;      FOR y = y_spt, 0.90*(y_xpt-y_spt)+y_spt, (y_xpt-y_spt)/FLOAT(nseg) DO BEGIN
;
;        t = (y - y_spt) / (y_xpt - y_spt)
;        x = (1.0 - t) * x_spt + t * x_xpt
;
;;        IF (t LT 0.02 ) THEN CONTINUE
;
;        a = MIN(ABS(inv.grid.x-(x-0.10)),ix0)
;        a = MIN(ABS(inv.grid.x- x      ),ix1)
;        a = MIN(ABS(inv.grid.x-(x+0.10)),ix2)
;        a = MIN(ABS(inv.grid.y- y      ),iy1)
;
;        print,ix0,ix1,ix2,x,y,t
;
;        xdata = inv.grid.x   [ix0:ix2]
;        ydata = inv.grid.data[ix0:ix2,iy1]
;
;;        wset,0
;;        loadct,5 ; 3
;;        inv.grid.data[ix0:ix2,iy1] = MAX(inv.grid.data)
;;        tvscl,inv.grid.data
;
;        i = WHERE(ydata GE 0.5*MAX(ydata))
;        ymax = MAX(ydata,imax)
;
;        a = [1.0,10000.0,xdata[imax],xdata[imax],0.01]  ; Resonable guesses...
;
;        weights = MAKE_ARRAY(N_ELEMENTS(xdata[i]),/FLOAT,value=1.0)                
;           print,'STRANGE PROBLEM HERE, AGAIN!'        
;           stop
;;        yfit = MPCURVEFIT(xdata[i], ydata[i]/ymax, weights, a, sigma,   $
;;                          CHISQ=chisq,FUNCTION_NAME='gfunct_new',/NODERIVATIVE,/QUIET)          
;
;        IF (chisq GT 0.1) THEN BEGIN
;
;          a = MAX(ydata[i],imax)
;
;;          wset,1
;;          plot,xdata,ydata,ystyle=1,yrange=[0.0,1.1*MAX(ydata)]
;;          oplot,xdata[i],ydata[i],psym=6
;
;        ENDIF ELSE BEGIN
;;          print,a,xdata[imax]
;
;;          safe_colors,/first
;;          wset,1
;;          plot,xdata,ydata,ystyle=1,yrange=[0.0,1.1*MAX(ydata)]
;;          oplot,xdata[i],ydata[i],psym=6
;;          oplot,xdata[i],yfit*ymax,color=2
;;          stop
;
;          a = MAX(yfit,imax)
;        ENDELSE
;
;        xdata1 = xdata[i]                           ; Cleaner way to code this..?
;        a = MIN(ABS(inv.grid.x-xdata1[imax]),xind)
;        a = MIN(ABS(inv.grid.y-y           ),yind)
;
;        IF (y EQ y_spt) THEN BEGIN
;          ix = xind
;          iy = yind
;        ENDIF ELSE BEGIN
;          ix = [ix,xind]
;          iy = [iy,yind]
;          i = N_ELEMENTS(ix)
;          wset,0
;          plots,ix[i-2:i-1],iy[i-2:i-1],color=3,thick=2,psym=6,/device
;        ENDELSE
;      ENDFOR
;      END
  ENDCASE

; Do a quick fit along the divertor leg:
  IF (1) THEN BEGIN
;    nx = N_ELEMENTS(inv.x)
;    ny = N_ELEMENTS(inv.y)
;     xtri = MAKE_ARRAY(nx*ny,/FLOAT)
;    ytri = MAKE_ARRAY(nx*ny,/FLOAT)
;    dtri = MAKE_ARRAY(nx*ny,/FLOAT)
;    FOR j = 0, ny-1 DO BEGIN
;      FOR i = 0, nx-1 DO BEGIN
;        xtri[i+j*nx] = inv.x[i]
;        ytri[i+j*nx] = inv.y[j]
;        dtri[i+j*nx] = inv.data[i,j]
;      ENDFOR
;    ENDFOR
;    TRIANGULATE, xtri, ytri, tr, b
    xout = x[0] + (FINDGEN(101)) / 100.0 * (x[nseg-1] - x[0])
    yout = y[0] + (FINDGEN(101)) / 100.0 * (y[nseg-1] - y[0]) 
    dist = SQRT((xout - x[0])^2 + (yout - y[0])^2)

    ix_out = FLOAT(ix[0]) + (FINDGEN(101)) / 100.0 * FLOAT(ix[nseg-1] - x[0])
    iy_out = FLOAT(iy[0]) + (FINDGEN(101)) / 100.0 * FLOAT(iy[nseg-1] - y[0])

    n = N_ELEMENTS(xout)
;   Interpolation:
;    fit = TRIGRID(xtri, ytri, dtri, tr, xout=xout, yout=yout)  
;    line_fit = MAKE_ARRAY(n,/FLOAT,VALUE=0.0)
;    FOR i = 0, n-1 DO line_fit[i] = fit[i,i]
;   Simple fit:  - just as good as the interpolation if the spatial resolution is high, which it is...
    line_fit2 = MAKE_ARRAY(n,/FLOAT,VALUE=0.0)
    FOR i = 0, n-1 DO BEGIN
      a = MIN(ABS(xout[i]-inv.x),ixmin)
      a = MIN(ABS(yout[i]-inv.y),iymin)
      line_fit2[i] = inv.data[ixmin,iymin] 
    ENDFOR    
  ENDIF

  contour = {  $
    n     : nseg      ,  $  ;
    ix    : ix        ,  $  ;
    iy    : iy        ,  $  ; 
    x     : x         ,  $  ; 
    y     : y         ,  $  ;
    x_spt : x_spt     ,  $  ; 
    y_spt : y_spt     ,  $  ;
    x_xpt : x_xpt     ,  $  ; 
    y_xpt : y_xpt     ,  $  ;
    fit_x    : xout   ,  $  ;
    fit_y    : yout   ,  $  ;
    dist     : dist   ,  $  ;
    data_ix  : ix_out ,  $  ;
    data_iy  : iy_out ,  $  ;
    data     : line_fit2 }     ;
;    fit2  : line_fit2 }     ;
 
  result = CREATE_STRUCT(inv, 'contour', contour)

  RETURN, result

END


;
; ======================================================================
;
PRO gfunct_new, x, a, f, pder

  alpha = A[0] / (1.0 + A[1] * (x - A[2])^2)

  beta = EXP(-(x - A[3])^2 / (2.0 * A[4]^2))

  f = alpha * beta

;  print,a
;  print,x
;  print,alpha


  IF N_PARAMS() GE 4 THEN BEGIN

    pder = [ [(a[1]*(x-a[2]))^2],  [2.0*a[0]*a[1]*(x-a[2])^2],  [-2.0*a[0]*a[1]^2*(x-a[1])], [replicate(1.0, N_ELEMENTS(X))] ]




;    help,pder,/struct
;    PRINT,pder
;    stop
  ENDIF
  
END  
;
; ======================================================================
;
FUNCTION ray_old,                 $
    shot=shot              ,  $  ;
    frame=frame            ,  $  ;
    channel=channel        ,  $  ;
    ext=ext                ,  $  ;
    casename=casename      ,  $  ;
    datatag=datatag        ,  $  ;
    region=region          ,  $  ;  Region of interest
    sp_peak=sp_peak        ,  $  ;  Region of interest
    full=full              ,  $  ;
    save=save              ,  $  ;
    chisq_limit=chisq_limit      ;

;  PRINT,shot,frame,channel,region

  IF (NOT KEYWORD_SET(chisq_limit)) THEN chisq_limit = 0.2

;  inv = GetInversion(shot=shot,channel=channel,frame=frame,   $
;                     ext=ext,casename=casename,region=region, $
;                     datatag=datatag,/full)

  inv = ContourImage(inv,sp_peak=sp_peak)

  inv = ProcessImage(inv, chisq_limit)

  IF (KEYWORD_SET(save)) THEN BEGIN
    fname = './PSI08_INV_'+                            $
            STRTRIM(STRING(inv.shot   ),1)+'_'+  $
            STRTRIM(STRING(inv.frame  ),1)+'_'+  $
            STRTRIM(STRING(inv.channel),1)+'_'+  $
            STRTRIM(STRING(region     ),1)+'.sav'
    SAVE,filename=fname,inv
  ENDIF

  return, inv

END


;
;
;

