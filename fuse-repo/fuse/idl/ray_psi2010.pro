;
; ======================================================================
;
; To start:
;
; > ln -s ~slisgo/fuse_data/mast/idl_data ~/fuse/idl/ray_data
; > cd ~/fuse/idl      
; > idl
; > @cortex_make  ; *** NEW ***
; > @ray_make
; > .r grid_separatrix  ; *** NEW ***
; > .r ray_psi2010
;
; Examples:
;
; result=ray_psi2010_process(1,/a_only,/plots)  align the strike- and x-points
; result=ray_psi2010_process(1)                 process both lines and store the data in ./ray_data/ray_<ID>.sav
; ray_psi2010_plots,A,B,param1=param1		A - reconstruction data index (see below in file), 
;					        B - plot number : 
;     1 - Basic comparison plot.  Aspect ration of the poloidal reconstructions not correct.
;     	  param1 - change the vertical line where the reconstructions are sampled (horizontal pixel number)
;
;         ray_psi2010_plots,1,1,cutoff=0.1,/equ
; 
;     2 - Plot of poloidal reconstruction to a screen window with proper aspect ratio.
;         YOU NEED TO SPECIFY /show_a or /show_b
;         you can turn off the dotted line with /no_line
;
;         ray_psi2010_plots,12,2,/equ,/show_a,/no_line,cutoff=0.1
;         ray_psi2010_plots,1,1,/equ,/show_b,/no_line,cutoff=0.5   
;         ray_psi2010_plots,51,2,/equ,/show_a,/no_line               x-point almost looks right, shifted in a bit perhaps
;         ray_psi2010_plots,51,2,/equ,/show_b
;         ray_psi2010_plots,51,2,equ='24861/24861_240.equ',/show_b
;         ray_psi2010_plots,51,2,/equ,/show_a
;
;     3 - Shaded surface plot, need to say /show_a or /show_b.
;     4 - Plot of reconstruction with separatrix overlayed -- work in progress.
;     5 - Plot of the two reconstructions showing where they intersect.
;         line 'a' - dark gray, line 'b' - light gray, both added - white
;         Can specify /show_a or /show_b.
;         Need to set the cutoff for each line, i.e. the fraction of the peak value above which
;         profile is set to a given value and below which it is set to zero, to isolate the basic 
;         contour of the emission -- this will be different for each line -- see example below.         
;
;         ray_psi2010_plots,71,5,cutoff=[0.05,0.20]
;         ray_psi2010_plots,71,5,cutoff=[0.05,0.20],/show_a
;         ray_psi2010_plots,1,5,cutoff=[0.02,0.015]           trying to see if Da and Dg overlap for attached case, for relative spatial calibration
;         ray_psi2010_plots,51,5,cutoff=[0.10,0.05]           overlap of CII and CIII for detached reference
;         ray_psi2010_plots,22,5,cutoff=[0.015,0.0125]        D_a and D_g overlap
;
;
;         OSM inner divertor contours:
;         ----------------------------
;         ray_psi2010_plots,51,5,cutoff=[0.30,0.20],/equ      
;         ray_psi2010_plots,51,2,/show_a,cutoff=0.60,/equ,/no_line   CII drops to background levels (approx.) at R=0.43, vertical alignment with target (approx.)
;         ray_psi2010_plots,51,2,/show_b,cutoff=0.40,/equ,/no_line   CIII
;         ray_psi2010_plots,22,2,/show_a,cutoff=0.70,/equ,/no_line   Dalpha
;         ray_psi2010_plots,22,2,/show_b,cutoff=0.40,/equ,/no_line   Dgamma
;
;
;     	OSM:
;       ---
;	ray_psi2010_plots,52,2,/equ,/show_a,cutoff=[0.7],/no_line
; 	ray_psi2010_plots,52,2,/equ,/show_b,cutoff=[0.5],/no_line
;	ray_psi2010_plots,23,2,/equ,/show_b,cutoff=[0.5],/no_line
;	ray_psi2010_plots,23,2,/equ,/show_a,cutoff=[0.9],/no_line
;
;	ray_psi2010_plots,1 ,1,/equ,/no_line
;	ray_psi2010_plots,1 ,1,/equ,/no_line,/ratio
;	ray_psi2010_plots,23,1,/equ,/no_line
;	ray_psi2010_plots,23,1,/equ,/no_line,/ratio
;
;
; ray_psi2010_pass				reprocess and save all reconstructions
; ray_psi2010_output,'filename'			put all B=1 plots into a postscript file in ./ray_data
;
;
; ======================================================================
;
PRO ray_psi2010_plot_surfaces, surfaces, color=color

  IF (NOT KEYWORD_SET(color)) THEN color = 'White'

  line_color = color

; Use discontinuities in the contours to isolate flux surfaces:
  ndata = N_ELEMENTS(TAG_NAMES(surfaces))      
print,'ndata= ',ndata
  help,surfaces,/struct

  FOR l = 1, ndata DO BEGIN
    data  = cortex_ExtractStructure(surfaces,l)
    ncont = N_ELEMENTS(TAG_NAMES(data))      
    FOR i = 1, ncont DO BEGIN
      cont = cortex_ExtractStructure(data,i)
      j = WHERE(cont.dist GT 0.1, count_break)
      j = [-1,j,cont.n-1]
      FOR k = 0, count_break-2 DO  $
        OPLOT,cont.x[j[k]+1:j[k+1]],cont.y[j[k]+1:j[k+1]],  $
              COLOR=TrueColor(line_color),LINESTYLE=2        
    ENDFOR
    line_color = 'Lightgreen'
  ENDFOR
END
;
; ======================================================================
;
FUNCTION ray_psi2010_flux_surfaces, equ

  shifts = [+0.004, +0.002, 0.000, -0.002, -0.004]
;  shifts = 0.000

  print,equ

  FOR j = 0, N_ELEMENTS(equ)-1 DO BEGIN
    xpoint_zone = [0.5,0.8,-1.4,1.4]
    b = grid_ReadEQUFile('~/fuse_data/mast/shots/'+ equ[j]) ; 25029/25029_312.equ')
    b = grid_FindNullPoints (b,xpoint_zone,1)
    b = grid_AnalyseBoundary(b,xpoint_zone,1)

    FOR i = 0, N_ELEMENTS(shifts)-1 DO BEGIN
      cont = grid_ExtractContour(b.psi, b.x, b.y, b.psi_b+shifts[i])
      name = 'data' + STRING(i+1,FORMAT='(I0)')
      IF (i EQ 0) THEN conts = CREATE_STRUCT(      name,cont) ELSE  $
                       conts = CREATE_STRUCT(conts,name,cont)
    ENDFOR

    name = 'data' + STRING(j+1,FORMAT='(I0)')
    IF (j EQ 0) THEN result = CREATE_STRUCT(       name,conts) ELSE  $
                     result = CREATE_STRUCT(result,name,conts)
  ENDFOR

  RETURN, result
END
;
;
; ======================================================================
;
PRO ray_psi2010_contour, data_2, colorct, color, nlevels, c_colors, xpos, ypos,  $
                         title=title, fill=fill, no_line=no_line, equ=equ, a=a, b=b, large=large,  $
                         text_color=text_color

  IF (KEYWORD_SET(a)) THEN BEGIN
    data = data_2.a
    file = data_2.afile
  ENDIF ELSE BEGIN
    data = data_2.b
    file = data_2.bfile
  ENDELSE

  IF (KEYWORD_SET(large)) THEN charsize = 2.0 ELSE charsize = 1.0

  xmin = MIN(data.x)
  xmax = MAX(data.x)
  ymin = MIN(data.y)
  ymax = MAX(data.y)
  xdelta = xmax - xmin
  ydelta = ymax - ymin
  xmin = xmin - 0.03 * xdelta
  xmax = xmax + 0.03 * xdelta
  ymin = ymin - 0.03 * ydelta
  ymax = ymax + 0.03 * ydelta

  PLOT,[xmax,xmin],[ymin,ymax], /NODATA, CHARSIZE=(0.5*(charsize-1.0)+1.0) , $ 
       XSTYLE=1,YSTYLE=1,color=Truecolor(color),TITLE=title,XTITLE='R (m)',YTITLE='Z (m)'
  LOADCT,colorct
  levels = (FINDGEN(nlevels) / FLOAT(nlevels)) * (MAX(data.data) - MIN(data.data)) + MIN(data.data)
  CONTOUR,data.data,data.x,data.y,  $
          levels=levels,c_colors=c_colors,/FILL,/OVERPLOT
  IF (NOT KEYWORD_SET(no_line)) THEN  $
    OPLOT,data.contour.x,data.contour.y,LINESTYLE=1,color=Truecolor('White')

  str = STRSPLIT(file,'/',/EXTRACT)
  n = N_ELEMENTS(str)
  IF (NOT KEYWORD_SET(text_color)) THEN text_color = 'White'
  FOR i = 0, n-1 DO BEGIN
    XYOUTS,xpos,ypos-0.015*charsize*FLOAT(i),STRTRIM(str[i],2),/NORMAL,  $
           color=Truecolor(text_color),CHARSIZE=charsize
  ENDFOR
  IF (KEYWORD_SET(equ)) THEN BEGIN
    FOR j = 0, N_ELEMENTS(equ)-1 DO BEGIN
      XYOUTS,xpos,ypos-0.015*charsize*FLOAT(i+j),equ[j],/NORMAL,  $
       color=Truecolor(text_color),CHARSIZE=charsize
    ENDFOR
  ENDIF
END
;
; ======================================================================
;
PRO ray_psi2010_plots, data, option, ascale=ascale, bscale=bscale, param1=param1, fill=fill, ps=ps,  $
               show_a=show_a,show_b=show_b,  $
               nlevels=nlevels,cutoff=cutoff,no_line=no_line, equ=equ, ratio=ratio

  CASE (SIZE(data,/TNAME)) OF
    'INT': BEGIN
       file = 'ray_data/ray_'+STRTRIM(STRING(data),2)+'.sav'
       RESTORE,file,/VERBOSE
       data = result
       END
    'LONG': BEGIN
       file = 'ray_data/ray_'+STRTRIM(STRING(data),2)+'.sav'
       RESTORE,file,/VERBOSE
       data = result
       END
     ELSE: 
  ENDCASE

  IF (NOT KEYWORD_SET(title  )) THEN title   = data.title
  IF (NOT KEYWORD_SET(ascale )) THEN ascale  = data.ascale ELSE ascale = 1.0
  IF (NOT KEYWORD_SET(bscale )) THEN bscale  = data.bscale ELSE ascale = 1.0
  IF (NOT KEYWORD_SET(color  )) THEN  $
    IF (option EQ 5) THEN color = 0 ELSE color = 5
  IF (NOT KEYWORD_SET(nlevels)) THEN nlevels = 50
  IF (KEYWORD_SET(equ)) THEN BEGIN
    IF (N_ELEMENTS(equ) EQ 1) THEN IF (equ EQ 1) THEN equ = data.equ 
    surfaces = ray_psi2010_flux_surfaces(equ)
  ENDIF
  IF (KEYWORD_SET(cutoff)) THEN BEGIN
    IF (N_ELEMENTS(cutoff) EQ 1) THEN cutoff = [cutoff,cutoff]
    vala = cutoff[0] * MAX(data.a.data)
    valb = cutoff[1] * MAX(data.b.data)
    i = WHERE(data.a.data GE vala, count)
    IF (count GT 0) THEN data.a.data[i] = vala
    i = WHERE(data.b.data GE valb, count)
    IF (count GT 0) THEN data.b.data[i] = valb
  ENDIF


  c_colors = LONG(FINDGEN(nlevels) * (255.0 / FLOAT(nlevels-1)) )

  !P.BACKGROUND = Truecolor('White')

  CASE option OF
;   --------------------------------------------------------------------
    1: BEGIN  ; comparison plot
      IF (NOT KEYWORD_SET(ps)) THEN BEGIN
        WINDOW,2,RETAIN=2,XSIZE=800,YSIZE=800
        DEVICE, DECOMPOSED=0
      ENDIF

      IF (NOT KEYWORD_SET(param1)) THEN ix = data.a.contour.ix[1] ELSE  $
                                        ix = LONG(param1)
      PRINT,'IX=',ix

      !P.MULTI = [0,2,2]

      ray_psi2010_contour, data, /a, color, 'Black', nlevels, c_colors, 0.11, 0.93, fill=fill, equ=equ, no_line=no_line
      IF (KEYWORD_SET(equ)) THEN ray_psi2010_plot_surfaces, surfaces
      ray_psi2010_contour, data, /b, color, 'Red'  , nlevels, c_colors, 0.61, 0.93, fill=fill, equ=equ, no_line=no_line
      IF (KEYWORD_SET(equ)) THEN ray_psi2010_plot_surfaces, surfaces

      XYOUTS,0.5,0.980,title,/NORMAL,color=Truecolor('Black'), ALIGNMENT=0.5, CHARSIZE=1.2

      ymax = MAX([data.a.contour.data*ascale, data.b.contour.data*bscale])
      CASE 1 OF
        0: BEGIN
          PLOT ,data.a.contour.data_ix,data.a.contour.data*ascale,YRANGE=[0.0,ymax],  $
                XTITLE='x index',color=Truecolor('Black')
          OPLOT,data.b.contour.data_ix,data.b.contour.data*bscale,color=Truecolor('Red')
          OPLOT,[ix,ix],[0,10000],color=Truecolor('Blue'),LINESTYLE=1
          END
        1: BEGIN
          x0 = data.b.contour.x[0]
          xdata1 = data.a.contour.fit_x - x0
          xdata2 = data.b.contour.fit_x - x0
          ydata1 = data.a.contour.data * ascale
          ydata2 = data.b.contour.data * bscale
          PLOT ,xdata1,ydata1,YRANGE=[0.0,ymax],color=Truecolor('Black'),XSTYLE=1,YSTYLE=1, $
                XTITLE='radial distance from inner target (m)',YTITLE='(arb)', $
                TITLE='Profile along inner divertor leg (dotted line)'
          OPLOT,xdata2,ydata2,color=Truecolor('Red')
          IF (NOT KEYWORD_SET(ratio)) THEN BEGIN
            xdata = [data.b.x[ix]-x0, data.b.x[ix]-x0]
            ydata = [0.0, 10000.0]
            OPLOT,xdata,ydata,color=Truecolor('Blue'),LINESTYLE=1
          ENDIF
          END
      ENDCASE

      ymax = MAX([data.a.data[ix,*]*ascale,data.b.data[ix,*]*bscale])
      option = 1
      IF (KEYWORD_SET(ratio)) THEN option = 2
      CASE option OF
        0: BEGIN
          PLOT ,data.a.data[ix,*]*ascale,YRANGE=[0.0,ymax],color=Truecolor('Black')
          OPLOT,data.b.data[ix,*]*bscale,                  color=Truecolor('Red')
          END
        1: BEGIN
          xdata1 = data.a.y
          xdata2 = data.b.y
          ydata1 = data.a.data[ix,*] * ascale
          ydata2 = data.b.data[ix,*] * bscale
          PLOT ,xdata1,ydata1,YRANGE=[0.0,ymax],color=Truecolor('Black'), XSTYLE=1,YSTYLE=1, $
                XTITLE='vertical position (m)',YTITLE='(arb)',  $
                TITLE='Vertical profile along blue dotted line'
          OPLOT,xdata2,ydata2,color=Truecolor('Red')
          END
        2: BEGIN
          x0 = data.b.contour.x[0]
          xdata1 = data.a.contour.fit_x - x0
          xdata2 = data.b.contour.fit_x - x0
          ydata1 = data.a.contour.data * ascale
          ydata2 = data.b.contour.data * bscale
          ydata = SMOOTH(ydata1 / ydata2 * 80.0,11)
          PLOT ,xdata1,ydata,YRANGE=[0.0,140.0],color=Truecolor('Black'),XSTYLE=1,YSTYLE=1, $
                XTITLE='radial distance from inner target (m)',YTITLE='(arb)', $
                TITLE='D_alpha / D_gamma ratio (scaled, smoothed)'
          END
      ENDCASE
      !P.MULTI = 0
      END
;   --------------------------------------------------------------------
    2: BEGIN  ; Big contour plot

      IF (NOT KEYWORD_SET(ps)) THEN BEGIN
        IF (KEYWORD_SET(show_a)) THEN dim = SIZE(data.a.data,/DIMENSIONS) 
        IF (KEYWORD_SET(show_b)) THEN dim = SIZE(data.b.data,/DIMENSIONS) 
        WINDOW,2,RETAIN=2,XSIZE=dim[0]*1.8,YSIZE=dim[1]*1.8
        DEVICE, DECOMPOSED=0
      ENDIF
      IF (KEYWORD_SET(show_a)) THEN ray_psi2010_contour, data, /a, color, 'Black', nlevels,c_colors, 0.16, 0.90, fill=fill, title=title, no_line=no_line, equ=equ, /large
      IF (KEYWORD_SET(show_b)) THEN ray_psi2010_contour, data, /b, color, 'Black', nlevels,c_colors, 0.16, 0.90, fill=fill, title=title, no_line=no_line, equ=equ, /large

      IF (KEYWORD_SET(equ)) THEN ray_psi2010_plot_surfaces, surfaces

      END
;   --------------------------------------------------------------------
    3: BEGIN  ; Surface plot
      IF (NOT KEYWORD_SET(ps)) THEN BEGIN
        WINDOW,2,RETAIN=2,XSIZE=800,YSIZE=800
        DEVICE, DECOMPOSED=0
      ENDIF

      IF (NOT (KEYWORD_SET(show_a) OR KEYWORD_SET(show_b))) THEN BEGIN
        PRINT, 'NEED TO SAY /show_a OR /show_b'
        STOP
      ENDIF
      IF (KEYWORD_SET(show_a)) THEN file = data.afile
      IF (KEYWORD_SET(show_b)) THEN file = data.bfile
      IF (KEYWORD_SET(show_a)) THEN data = data.a
      IF (KEYWORD_SET(show_b)) THEN data = data.b

      LOADCT,3
      !P.BACKGROUND = 255
      SHADE_SURF, data.data, data.x, data.y, AZ=-45.0, AX=30.0, color=0,  $
                  XTITLE='R (m)', YTITLE='Z (m)', ZTITLE = 'arb',  $
                  CHARSIZE = 3.0, ZSTYLE=1,YSTYLE=1,XSTYLE=1

      XYOUTS,0.50,0.98,title,/NORMAL,color=Truecolor('Black'),CHARSIZE=1.5, ALIGNMENT=0.5
      XYOUTS,0.50,0.96,file ,/NORMAL,color=Truecolor('Black') , ALIGNMENT=0.5
      END
;   --------------------------------------------------------------------
    4: BEGIN  ; Big contour plot with separatrix over plotted
      IF (NOT KEYWORD_SET(ps)) THEN BEGIN
        IF (KEYWORD_SET(show_a)) THEN dim = SIZE(data.a.data,/DIMENSIONS) 
        IF (KEYWORD_SET(show_b)) THEN dim = SIZE(data.b.data,/DIMENSIONS) 
        WINDOW,2,RETAIN=2,XSIZE=dim[0]*1.5,YSIZE=dim[1]*1.5
        DEVICE, DECOMPOSED=0
      ENDIF

      IF (NOT (KEYWORD_SET(show_a) OR KEYWORD_SET(show_b))) THEN BEGIN
        PRINT, 'NEED TO SAY /show_a OR /show_b'
        STOP
      ENDIF
      IF (KEYWORD_SET(show_a)) THEN ray_psi2010_contour, data.a, data.afile, color, 'Black', nlevels,c_colors, 0.16, 0.90, fill=fill, title=title
      IF (KEYWORD_SET(show_b)) THEN ray_psi2010_contour, data.b, data.bfile, color, 'Black', nlevels,c_colors, 0.16, 0.90, fill=fill, title=title

      file = '~/fuse_data/mast/shots/25028/25028_312.equ'
      b = grid_readequfile(file)
      CONTOUR, b.psi, b.x, b.y, levels=[b.psi_boundary], color=TrueColor('Red'), /OVERPLOT

      file = '~/fuse_data/mast/shots/25029/25029_312.equ'
      b = grid_readequfile(file)
      CONTOUR, b.psi, b.x, b.y, levels=[b.psi_boundary], color=TrueColor('Green'), /OVERPLOT

      END
;   --------------------------------------------------------------------
    5: BEGIN  ; show intersection between images

      IF (NOT KEYWORD_SET(ps)) THEN BEGIN
        dim = SIZE(data.a.data,/DIMENSIONS)
        WINDOW,2,RETAIN=2,XSIZE=dim[0]*1.8,YSIZE=dim[1]*1.8
        DEVICE, DECOMPOSED=0
      ENDIF

      IF (NOT KEYWORD_SET(cutoff)) THEN cutoff = [0.1,0.1]

      vala = cutoff[0] * MAX(data.a.data)
      valb = cutoff[1] * MAX(data.b.data)

      data.a.data[WHERE(data.a.data LT vala)] = 0.0
      data.a.data[WHERE(data.a.data GE vala)] = 1.0
      data.b.data[WHERE(data.b.data LT valb)] = 0.0
      data.b.data[WHERE(data.b.data GE valb)] = 2.0

      IF (NOT KEYWORD_SET(show_a) AND NOT KEYWORD_SET(show_b)) THEN BEGIN
        data.a.data = data.a.data + data.b.data
        data.afile = data.afile + ' / ' + data.bfile
      ENDIF ELSE BEGIN
        IF (KEYWORD_SET(show_b)) THEN BEGIN
          data.a     = data.b
          data.afile = data.bfile
        ENDIF
      ENDELSE

      ray_psi2010_contour, data, /a, color, 'Black', nlevels,c_colors,  $
                           0.16, 0.90, fill=fill, title=title, /no_line, equ=equ, /large,  $
                           text_color = 'Red'

      IF (KEYWORD_SET(equ)) THEN ray_psi2010_plot_surfaces, surfaces, color='Red' 

      END
;   --------------------------------------------------------------------
   ENDCASE
END
;
; ======================================================================
;
FUNCTION ray_psi2010_process,option,plots=plots,a_only=a_only,b_only=b_only,uberplots=uberplots

  result = -1

  !P.MULTI = [0,1,1]

  PRINT,'PROCESSING ID ',option

  CASE option OF
   1: BEGIN  
     title = 'D_a/D_g CALIBRATION, ATTACHED CONDITIONS: 24867 at 340 ms'
     equ = '24867/24867_335.equ'
     plot_option = 1
     fit_sample=10
     ascale = 1.0 
     bscale = 2.0
     aspt = [0.280,-1.435]
     axpt = [0.73 ,-1.20 ]
     bspt = [0.280,-1.425]
     bxpt = [0.73 ,-1.20 ]
     afile='FFC_24867_1225_3_HL01_rbc_Dalpha.cgm'
     bfile='FFC_24867_1225_1_HL07_rba_Dgamma.cgm'
     END
   10: BEGIN 
     title = 'REFERENCE D_alpha / D_gamma : 24861 at 201 ms'
     equ = '24861/24861_200.equ'
     plot_option = 1
     fit_sample=10
     ascale = 1.0 * (100.0 / 67.0) ; to match 24867, where the relative calibration was done
     bscale = 2.0 
     aspt = [0.280,-1.470]
     axpt = [0.70 ,-1.14 ]
     bspt = [0.280,-1.460]
     bxpt = [0.70 ,-1.14 ]
     afile='FFC_24861_554_3_HL01_rbc_Dalpha.cgm'
     bfile='FFC_24861_553_1_HL07_rba_Dgamma.cgm'
     END
   11: BEGIN 
     title = 'REFERENCE D_alpha / D_gamma : 24861 at 242 ms'
     equ = '24861/24861_240.equ'
     plot_option = 1
     fit_sample=10
     ascale = 1.0 * (100.0 / 67.0)
     bscale = 2.0 
     aspt = [0.280,-1.425]
     axpt = [0.67 ,-1.15 ]
     bspt = [0.280,-1.415]
     bxpt = [0.67 ,-1.15 ]
     afile='FFC_24861_762_3_HL01_rbc_Dalpha.cgm'
     bfile='FFC_24861_762_1_HL07_rba_Dgamma.cgm'
     END
   12: BEGIN 
     title = 'REFERENCE D_alpha / D_gamma : 24861 at 312 ms'
     equ = '24861/24861_315.equ'
     plot_option = 1
     fit_sample=10
     ascale = 1.0 * (100.0 / 67.0)
     bscale = 2.0 
     aspt = [0.280,-1.415]
     axpt = [0.71 ,-1.18 ]
     bspt = [0.280,-1.405]
     bxpt = [0.71 ,-1.18 ]
     afile='FFC_24861_1109_3_HL01_rbc_Dalpha.cgm'
     bfile='FFC_24861_1108_1_HL07_rba_Dgamma.cgm'
     END
   20: BEGIN
     title = 'REFERENCE REPEAT D_alpha : 24861 and 25028 at 201 ms'
     equ = '25028/25028_200equ'
     plot_option = 1
     fit_sample=10
     ascale = 1.0 
     bscale = 1.0 * (67.0 / 50.0)
     aspt = [0.280,-1.470]
     axpt = [0.70 ,-1.14 ]
     bspt = [0.280,-1.475]
     bxpt = [0.70 ,-1.14 ]
     afile='FFC_24861_554_3_HL01_rbc_Dalpha.cgm'
     bfile='FFC_25028_603_3_HL01_rbc_Dalpha.cgm'
     END    
;
;         24861    24862   24866   24867   25028
;    rbc   67 us    67 us   50 us  100 us   50 us
;    rdb    2 ms     3 ms    1 ms            3 ms
;    rba  200 us   200 us  100 us  200 us  167 ms
;    rdd  600 us   200 us  200 us          200 ms
;
   21: BEGIN
     title = 'REFERENCE REPEAT D_alpha / D_gamma : 25028 at 201 ms'             ; 21
     equ = '25028/25028_200.equ'
     plot_option = 1
     fit_sample=10
     ascale = 1.0 * (67.0  / 50.0 )
     bscale = 2.0 * (200.0 / 167.0) 
     aspt = [0.280,-1.4775]
     axpt = [0.70 ,-1.14 ]
     bspt = [0.280,-1.465]
     bxpt = [0.70 ,-1.14 ]
     afile='FFC_25028_603_3_HL01_rbc_Dalpha.cgm'
     bfile='FFC_25028_603_1_HL07_rba_Dgamma.cgm'
     END    
   22: BEGIN
     title = 'REFERENCE REPEAT D_alpha / D_gamma : 25028 at 242 ms'             ; 22
     equ = '25028/25028_250.equ'
     plot_option = 1
     fit_sample=10
     ascale = 1.0 * (67.0  / 50.0 )
     bscale = 2.0 * (200.0 / 167.0) 
     aspt = [0.280,-1.430]
     axpt = [0.705,-1.155]
     bspt = [0.280,-1.420]
     bxpt = [0.710,-1.140]
     afile='FFC_25028_810_3_HL01_rbc_Dalpha.cgm'
     bfile='FFC_25028_810_1_HL07_rba_Dgamma.cgm'
     END    
   23: BEGIN
     title = 'REFERENCE REPEAT D_alpha / D_gamma : 25028 at 312 ms'             ; 23
     equ = '25028/25028_310.equ'
     plot_option = 1
     fit_sample=10
     ascale = 1.0 * (67.0  / 50.0 )
     bscale = 2.0 * (200.0 / 167.0) 
     aspt = [0.280,-1.420]
     axpt = [0.725,-1.175]
     bspt = [0.280,-1.4075]
     bxpt = [0.730,-1.160]
     afile='FFC_25028_1158_3_HL01_rbc_Dalpha.cgm'
     bfile='FFC_25028_1158_1_HL07_rba_Dgamma.cgm'
     END    
   30: BEGIN
     title = 'REFERENCE and NO i/b GAS D_alpha : 24861 and 24862 at 201 ms'     ; 30
     equ = '24861/24861_200.equ'
     plot_option = 1
     fit_sample=10
     ascale = 1.0 
     bscale = 1.0 
     aspt = [0.280,-1.470]
     axpt = [0.70 ,-1.14 ]
     bspt = [0.280,-1.475]
     bxpt = [0.70 ,-1.16 ]
     afile='FFC_24861_554_3_HL01_rbc_Dalpha.cgm'
     bfile='FFC_24862_554_3_HL01_rbc_Dalpha.cgm'
     END
   31: BEGIN
     title = 'NO i/b GAS D_alpha / D_gamma : 24862 at 201 ms'
     equ = '24862/24862_200.equ'
     plot_option = 1
     fit_sample=10
     ascale = 1.0 * (100.0 / 67.0)
     bscale = 2.0 
     aspt = [0.280,-1.475]
     axpt = [0.70 ,-1.16 ]
     bspt = [0.280,-1.460]
     bxpt = [0.70 ,-1.16 ]
     afile='FFC_24862_554_3_HL01_rbc_Dalpha.cgm'
     bfile='FFC_24862_554_1_HL07_rba_Dgamma.cgm'
     END    
   40: BEGIN
     title = 'REFERENCE and DENSITY RAMP D_alpha : 24861 and 24866 at 201 ms'   ; 40
     equ = '24861/24861_200.equ'
     plot_option = 1
     fit_sample=10
     ascale = 1.0 
     bscale = 1.0 * (67.0 / 50.0)
     aspt = [0.280,-1.475]
     axpt = [0.70 ,-1.14 ]
     bspt = [0.280,-1.4725]
     bxpt = [0.70 ,-1.14 ]
     afile='FFC_24861_554_3_HL01_rbc_Dalpha.cgm'
     bfile='FFC_24866_554_3_HL01_rbc_Dalpha.cgm'
     END    
   41: BEGIN
     title = 'REFERENCE and DENSITY RAMP D_alpha : 24861 and 24866 at 242 ms'
     equ = '24861/24861_240.equ'
     plot_option = 1
     fit_sample=10
     ascale = 1.0 
     bscale = 1.0 * (67.0 / 50.0)
     aspt = [0.280,-1.425]
     axpt = [0.67 ,-1.15 ]
     bspt = [0.280,-1.415]
     bxpt = [0.685,-1.135]
     afile='FFC_24861_762_3_HL01_rbc_Dalpha.cgm'
     bfile='FFC_24866_766_3_HL01_rbc_Dalpha.cgm'
     END    
   42: BEGIN
     title = 'DENSITY RAMP D_alpha / D_gamma : 24866 at 201 ms'
     equ = '24866/24866_200.equ'
     plot_option = 1
     fit_sample=10
     ascale = 1.0 * (67.0  /  50.0)
     bscale = 2.0 * (200.0 / 100.0)
     aspt = [0.280,-1.4725]
     axpt = [0.70 ,-1.14 ]
     bspt = [0.280,-1.4575]
     bxpt = [0.70 ,-1.14 ]
     afile='FFC_24866_554_3_HL01_rbc_Dalpha.cgm'
     bfile='FFC_24866_553_1_HL07_rba_Dgamma.cgm'
     END    
   43: BEGIN
     title = 'DENSITY RAMP D_alpha / D_gamma : 24866 at 242 ms'
     equ = '24866/24866_240.equ'
     plot_option = 1
     fit_sample=10
     ascale = 1.0 * (67.0  /  50.0)
     bscale = 2.0 * (200.0 / 100.0)
     aspt = [0.280,-1.420]
     axpt = [0.685,-1.135]
     bspt = [0.280,-1.405]
     bxpt = [0.685,-1.135]
     afile='FFC_24866_766_3_HL01_rbc_Dalpha.cgm'
     bfile='FFC_24866_766_1_HL07_rba_Dgamma.cgm'
     END    
   44: BEGIN
     title = 'DENSIT RAMP D_alpha / D_gamma : 24866 at 271 ms'
     equ = '24866/24866_265.equ'
     plot_option = 1
     fit_sample=10
     ascale = 1.0 * (67.0  /  50.0)
     bscale = 2.0 * (200.0 / 100.0)
     aspt = [0.280,-1.410]
     axpt = [0.685,-1.135]
     bspt = [0.280,-1.390]
     bxpt = [0.685,-1.135]
     afile='FFC_24866_903_3_HL01_rbc_Dalpha.cgm'
     bfile='FFC_24866_903_1_HL07_rba_Dgamma.cgm'
     END    
   45: BEGIN
     title = 'DENSITY RAMP D_alpha / D_gamma : 24866 at 312 ms'
     equ = '24866/24866_315.equ'
     plot_option = 1
     fit_sample=10
     ascale = 1.0 * (67.0  /  50.0)
     bscale = 2.0 * (200.0 / 100.0)
     aspt = [0.280,-1.395]
     axpt = [0.685,-1.135]
     bspt = [0.280,-1.380]
     bxpt = [0.685,-1.135]
     afile='FFC_24866_1109_3_HL01_rbc_Dalpha.cgm'
     bfile='FFC_24866_1109_1_HL07_rba_Dgamma.cgm'
     END    
;  =====================================================================
;  =====================================================================
;  Carbon:
   50: BEGIN 
     title = 'REFERENCE CII and CIII : 25029 at 201 ms'
     equ = '25029/25029_200.equ'
     plot_option = 1
     fit_sample=10
     ascale = 5.0
     bscale = 1.0
     aspt = [0.280,-1.439]
     axpt = [0.705,-1.180]
     bspt = [0.280,-1.410]
     bxpt = [0.715,-1.165]
     afile='FFC_25029_603_3_HL01_rbc_CII.cgm'
     bfile='FFC_25029_603_1_HL07_rba_CIII.cgm'
     END
   51: BEGIN 
     title = 'REFERENCE CII and CIII : 25029 at 242 ms'
     equ = '25029/25029_255.equ'
     plot_option = 1
     fit_sample=10
     ascale = 5.0
     bscale = 1.0
     aspt = [0.280,-1.402]
     axpt = [0.705,-1.155]
     bspt = [0.280,-1.382]
     bxpt = [0.710,-1.140]
     afile='FFC_25029_810_3_HL01_rbc_CII.cgm'
     bfile='FFC_25029_810_1_HL07_rba_CIII.cgm'
     END
   52: BEGIN 
     title = 'REFERENCE CII and CIII : 25029 at 312 ms'
     equ = '25029/25029_310.equ'
     plot_option = 1
     fit_sample=10
     ascale = 5.0
     bscale = 1.0
     aspt = [0.280,-1.395]
     axpt = [0.725,-1.175]
     bspt = [0.280,-1.370]
     bxpt = [0.730,-1.160]
     afile='FFC_25029_1158_3_HL01_rbc_CII.cgm'
     bfile='FFC_25029_1158_1_HL07_rba_CIII.cgm'
     END
   60: BEGIN
     title = 'NO LOWER i/b PUFF CII and CIII : 24869 at 201 ms'
     equ = '24869/24869_205.equ'
     plot_option = 1
     fit_sample=10
     ascale = 5.0
     bscale = 2.0 
     aspt = [0.280,-1.470]
     axpt = [0.68 ,-1.18 ]
     bspt = [0.280,-1.467]
     bxpt = [0.70 ,-1.165]
     afile='FFC_24869_552_3_HL01_rbc_CII.cgm'
     bfile='FFC_24869_552_1_HL07_rba_CIII.cgm'
     END    
   70: BEGIN 
     title = 'REFERENCE Dgamma and CIII : 24028 and 25029 at 312 ms'
     equ = '24028/24028_310.equ'
     plot_option = 1
     fit_sample=10
     ascale = 1.0
     bscale = 2.0
     aspt = [0.280,-1.4075]
     axpt = [0.730,-1.160]
     bspt = [0.280,-1.370]
     bxpt = [0.730,-1.160]
     afile='FFC_25028_1158_1_HL07_rba_Dgamma.cgm'
     bfile='FFC_25029_1158_1_HL07_rba_CIII.cgm'
     END
   71: BEGIN 
     title = 'REFERENCE Dalpha and CII : 24028 and 25029 at 312 ms'
     equ = '25028/25028_310.equ'
     plot_option = 1
     fit_sample=10
     ascale = 1.0
     bscale = 5.0
     aspt = [0.280,-1.420]
     axpt = [0.725,-1.175]
     bspt = [0.280,-1.395]
     bxpt = [0.725,-1.175]
     afile='FFC_25028_1158_3_HL01_rbc_Dalpha.cgm'
     bfile='FFC_25029_1158_3_HL01_rbc_CII.cgm'
     END
   72: BEGIN 
     title = 'REFERENCE Dgamma and CII : 24028 and 25029 at 312 ms'
     equ = '25028/25028_310.equ'
     plot_option = 1
     fit_sample=10
     ascale = 1.0
     bscale = 2.0
     aspt = [0.280,-1.4075]
     axpt = [0.730,-1.160]
     bspt = [0.280,-1.395]
     bxpt = [0.725,-1.175]
     afile='FFC_25028_1158_1_HL07_rba_Dgamma.cgm'
     bfile='FFC_25029_1158_3_HL01_rbc_CII.cgm'
     END
   73: BEGIN 
     title = 'REFERENCE Dalpha and CIII : 24028 and 25029 at 312 ms'
     equ = '25028/25028_310.equ'
     plot_option = 1
     fit_sample=10
     ascale = 1.0
     bscale = 5.0
     aspt = [0.280,-1.420]
     axpt = [0.725,-1.175]
     bspt = [0.280,-1.370]
     bxpt = [0.730,-1.160]
     afile='FFC_25028_1158_3_HL01_rbc_Dalpha.cgm'
     bfile='FFC_25029_1158_1_HL07_rba_CIII.cgm'
     END
   80: BEGIN 
     title = 'REFERENCE Ddelta: 25028 and 25029 at 312 ms'
     equ = '25028/25028_310.equ'
     plot_option = 1
     fit_sample=10
     ascale = 1.0
     bscale = 1.0
     aspt = [0.280,-1.420]
     axpt = [0.725,-1.175]
     bspt = [0.280,-1.420]
     bxpt = [0.725,-1.175]
     afile='DIVCAM_25028_17_2_HL01_rdb_Ddelta.cgm'
     bfile='DIVCAM_25029_17_2_HL01_rdb_Ddelta.cgm'
     END
   81: BEGIN 
     title = 'REFERENCE Dbeta: 25028 and 25029 at 312 ms'
     equ = '25028/25028_210.equ'
     plot_option = 1
     fit_sample=10
     ascale = 1.0
     bscale = 1.0
     aspt = [0.280,-1.4075]
     axpt = [0.730,-1.160]
     bspt = [0.280,-1.4075]
     bxpt = [0.730,-1.160]
     afile='DIVCAM_25028_17_4_HL07_rdd_Dbeta.cgm'
     bfile='DIVCAM_25029_17_4_HL07_rdd_Dbeta.cgm'
     END

;  =====================================================================
;  =====================================================================
;  Comparing new and old inversions of 20100413:
   100: BEGIN
     title = 'REFERENCE REPASS D_alpha : 25028 at 312 ms'
     equ = '25028/25028_310.equ'
     plot_option = 1
     fit_sample=10
     ascale = 1.0 * (67.0  / 50.0 )
     bscale = 1.0 * (67.0  / 50.0 )
     aspt = [0.280,-1.420]
     axpt = [0.725,-1.175]
     bspt = [0.280,-1.420]
     bxpt = [0.725,-1.175]
     afile='FFC_25028_1158_3_HL01_rbc_Dalpha.cgm'
     bfile='old_FFC_25028_1158_3_HL01_rbc_Dalpha.cgm'
     END    
   101: BEGIN
     title = 'REFERENCE REPASS D_gamma : 25028 at 312 ms'
     equ = '25028/25028_310.equ'
     plot_option = 1
     fit_sample=10
     ascale = 2.0 * (200.0 / 167.0) 
     bscale = 2.0 * (200.0 / 167.0) 
     aspt = [0.280,-1.4075]
     axpt = [0.730,-1.160]
     bspt = [0.280,-1.4075]
     bxpt = [0.730,-1.160]
     afile='FFC_25028_1158_1_HL07_rba_Dgamma.cgm'
     bfile='old_FFC_25028_1158_1_HL07_rba_Dgamma.cgm'
     END    
   110: BEGIN 
     title = 'REFERENCE REPASS CII : 25029 at 312 ms'
     equ = '25029/25029_310.equ'
     plot_option = 1
     fit_sample=10
     ascale = 5.0
     bscale = 5.0
     aspt = [0.280,-1.395]
     axpt = [0.725,-1.175]
     bspt = [0.280,-1.395]
     bxpt = [0.725,-1.175]
     afile='FFC_25029_1158_3_HL01_rbc_CII.cgm'
     bfile='old_FFC_25029_1158_3_HL01_rbc_CII.cgm'
     END
   111: BEGIN 
     title = 'REFERENCE Repass CIII : 25029 at 312 ms'
     equ = '25029/25029_310.equ'
     plot_option = 1
     fit_sample=10
     ascale = 1.0
     bscale = 1.0
     aspt = [0.280,-1.4075]
     axpt = [0.730,-1.160]
     bspt = [0.280,-1.4075]
     bxpt = [0.730,-1.160]
     afile='FFC_25029_1158_1_HL07_rba_CIII.cgm'
     bfile='old_FFC_25029_1158_1_HL07_rba_CIII.cgm'
     END

  ENDCASE

  title = STRTRIM(STRING(option),2) + ': ' + title

  PRINT, '----------------------------------------------------------------------'
  PRINT, 'PROCESSING: ',title
  PRINT, '----------------------------------------------------------------------'

  a = -1
  b = -1

  IF (NOT KEYWORD_SET(b_only)) THEN  $
    a = ray(file=afile,shot=24860,region=2,fit_cutoff=0.5,fit_sample=fit_sample,spt=aspt,xpt=axpt,plots=plots)
  IF (NOT KEYWORD_SET(a_only)) THEN  $
    b = ray(file=bfile,shot=24860,region=2,fit_cutoff=0.5,fit_sample=fit_sample,spt=bspt,xpt=bxpt,plots=plots)

  result = { title : title, afile : afile, bfile : bfile , equ : equ, $
             a : a, b : b ,     $ 
             ascale : ascale ,  $
             bscale : bscale }

  SAVE,result,filename='ray_data/ray_'+STRTRIM(STRING(option),2)+'.sav'

  IF (KEYWORD_SET(uberplots)) THEN ray_psi2010_plots,result,plot_option

  RETURN, result

END
;
; ======================================================================
;
PRO ray_psi2010_pass, option=option

  IF (NOT KEYWORD_SET(option)) THEN  $
    option = [1, 10,11,12, 20,21,22,23, 30,31, 40,41,42,43,44,45, 50,51,52, 60, 70,71, 80,81, 100,101, 110,111]

  FOR i = 0, N_ELEMENTS(option)-1 DO result = ray_psi2010_process(option[i])
END
;
; ======================================================================
;
PRO ray_psi2010_output, ps, option=option, nlevels=nlevels, pos=pos

  PSOPEN, filename = 'ray_data/' + ps + '.ps'

  IF (NOT KEYWORD_SET(option)) THEN  $
    option = [1, 10,11,12, 20,21,22,23, 30,31, 40,41,42,43,44,45, 50,51,52, 60, 70,71, 80,81, 100,101, 110,111]

  FOR i = 0, N_ELEMENTS(option)-1 DO ray_psi2010_plots,option[i],1,param1=pos, ps='on'

  PSCLOSE
END
;
; ======================================================================
;
