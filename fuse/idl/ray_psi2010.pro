;
; ======================================================================
;
; To start:
;
; > ln -s ~slisgo/fuse_data/mast/idl_data ~/fuse/idl/data_ray
; > cd ~/fuse/idl      
; > idl
; > @ray_make
; > .r psi_2010
;
; Examples:
;
; result=ray_psi2010_process(1,/a_only,/plots)  align the strike- and x-points
; result=ray_psi2010_process(1)                 process both lines and store the data in ./data_ray/ray_<ID>.sav
; ray_psi2010_plots,A,B,param1=param1		A - reconstruction data index (see below in file), 
;					        B - plot number : 
;     1 - basic comparison plot
;     	  param1 - change the vertical line where the reconstructions are sampled (horizontal pixel number)
;     2 - plot of reconstruction to window with proper aspect ratio
;         YOU NEED TO SPECIFY /show_a or /show_b
;         you can turn off the dotten line with /no_line
;     3 - shaded surface plot, need to say /show_a or /show_b
;     4 - plot of reconstruction with separatrix overlayed -- work in progress
;     5 - plot of the two reconstructions showing where they intersect
;         line 'a' - dark gray, line 'b' - light gray, both added - white
;         Can specify /show_a or /show_b.
;         Need to set the cutoff for each line, i.e. the fraction of the peak value below which the
;         profile is set to zero, to isolate the basic contour of the emission -- this will be
;         different for each line -- see example below.
;
;        ray_psi2010_plots,71,5,cutoff=[0.05,0.20]
;        ray_psi2010_plots,71,5,cutoff=[0.05,0.20],/show_a
;
; ray_psi2010_pass				reprocess and save all reconstructions
; ray_psi2010_output,'filename'			put all B=1 plots into a postscript file in ./data_ray
;
;
;
; ======================================================================
;
PRO ray_psi2010_contour, data, file, colorct, color, nlevels, c_colors, xpos, ypos,  $
                         title=title, fill=fill, no_line=no_line

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

  PLOT,[xmax,xmin],[ymin,ymax], /NODATA, $ 
       XSTYLE=1,YSTYLE=1,color=Truecolor(color),TITLE=title,XTITLE='R (m)',YTITLE='Z (m)'
  LOADCT,colorct
  levels = (FINDGEN(nlevels) / FLOAT(nlevels)) * (MAX(data.data) - MIN(data.data)) + MIN(data.data)
  CONTOUR,data.data,data.x,data.y,  $
          levels=levels,c_colors=c_colors,/FILL,/OVERPLOT
  IF (NOT KEYWORD_SET(no_line)) THEN  $
    OPLOT,data.contour.x,data.contour.y,LINESTYLE=1,color=Truecolor('White')

;  IF (KEYWORD_SET(fill)) THEN file_colour = 'White' ELSE file_colour = 'Black'
  file_colour = 'White' 
  XYOUTS,xpos,ypos,file,/NORMAL,color=Truecolor(file_colour)

END
;
; ======================================================================
;
PRO ray_psi2010_plots, data, option, ascale=ascale, bscale=bscale, param1=param1, fill=fill, ps=ps,  $
               show_a=show_a,show_b=show_b,  $
               nlevels=nlevels,cutoff=cutoff,no_line=no_line

  CASE (SIZE(data,/TNAME)) OF
    'INT': BEGIN
       file = 'data_ray/ray_'+STRTRIM(STRING(data),2)+'.sav'
       RESTORE,file,/VERBOSE
       data = result
       END
    'LONG': BEGIN
       file = 'data_ray/ray_'+STRTRIM(STRING(data),2)+'.sav'
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

  c_colors = LONG(FINDGEN(nlevels) * (255.0 / FLOAT(nlevels-1)) )

  !P.BACKGROUND = Truecolor('White')

  CASE option OF
;   --------------------------------------------------------------------
    1: BEGIN  ; comparison plot
      IF (NOT KEYWORD_SET(ps)) THEN BEGIN
        WINDOW,2,RETAIN=2
        DEVICE, DECOMPOSED=0
      ENDIF

      IF (NOT KEYWORD_SET(param1)) THEN ix = data.a.contour.ix[1] ELSE  $
                                        ix = LONG(param1)
      PRINT,'IX=',ix

      !P.MULTI = [0,2,2]

      ray_psi2010_contour, data.a, data.afile, color, 'Black', nlevels, c_colors, 0.11, 0.93, fill=fill
      ray_psi2010_contour, data.b, data.bfile, color, 'Red'  , nlevels, c_colors, 0.61, 0.93, fill=fill

      XYOUTS,0.5,0.975,title,/NORMAL,color=Truecolor('Black'), ALIGNMENT=0.5, CHARSIZE=1.5

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
          xdata = [data.b.x[ix]-x0, data.b.x[ix]-x0]
          ydata = [0.0, 10000.0]
          OPLOT,xdata,ydata,color=Truecolor('Blue'),LINESTYLE=1
          END
      ENDCASE

      ymax = MAX([data.a.data[ix,*]*ascale,data.b.data[ix,*]*bscale])
      CASE 1 OF
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
      ENDCASE
     !P.MULTI = 0
      END
;   --------------------------------------------------------------------
    2: BEGIN  ; Big contour plot
      IF (NOT KEYWORD_SET(ps)) THEN BEGIN
        IF (KEYWORD_SET(show_a)) THEN dim = SIZE(data.a.data,/DIMENSIONS) 
        IF (KEYWORD_SET(show_b)) THEN dim = SIZE(data.b.data,/DIMENSIONS) 
        WINDOW,2,RETAIN=2,XSIZE=dim[0]*1.5,YSIZE=dim[1]*1.5
        DEVICE, DECOMPOSED=0
      ENDIF
      IF (KEYWORD_SET(show_a)) THEN ray_psi2010_contour, data.a, data.afile, color, 'Black', nlevels,c_colors, 0.16, 0.90, fill=fill, title=title, no_line=no_line
      IF (KEYWORD_SET(show_b)) THEN ray_psi2010_contour, data.b, data.bfile, color, 'Black', nlevels,c_colors, 0.16, 0.90, fill=fill, title=title, no_line=no_line
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
        WINDOW,2,RETAIN=2,XSIZE=dim[0]*1.7,YSIZE=dim[1]*1.7
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
        data_c = data.a
        data_c.data = data.a.data + data.b.data
      ENDIF ELSE BEGIN
        IF (KEYWORD_SET(show_a)) THEN data_c = data.a
        IF (KEYWORD_SET(show_b)) THEN data_c = data.b
      ENDELSE

      ray_psi2010_contour, data_c, file, color, 'Black', nlevels,c_colors,  $
                           0.16, 0.90, fill=fill, title=title, /no_line

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
     title = 'REFERENCE REPEAT D_alpha / D_gamma : 25028 at 201 ms'
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
     title = 'REFERENCE REPEAT D_alpha / D_gamma : 25028 at 242 ms'
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
     title = 'REFERENCE REPEAT D_alpha / D_gamma : 25028 at 312 ms'
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
     title = 'REFERENCE and NO i/b GAS D_alpha : 24861 and 24862 at 201 ms'
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
     title = 'REFERENCE and DENSITY RAMP D_alpha : 24861 and 24866 at 201 ms'
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

  result = { title : title, afile : afile, bfile : bfile ,  $
             a : a, b : b ,     $ 
             ascale : ascale ,  $
             bscale : bscale }

  SAVE,result,filename='data_ray/ray_'+STRTRIM(STRING(option),2)+'.sav'

  IF (KEYWORD_SET(uberplots)) THEN ray_psi2010_plots,result,plot_option

  RETURN, result

END
;
; ======================================================================
;
PRO ray_psi2010_pass, option=option

  IF (NOT KEYWORD_SET(option)) THEN  $
    option = [1, 10,11,12, 20,21,22,23, 30,31, 40,41,42,43,44,45, 50,51,52, 60]

  FOR i = 0, N_ELEMENTS(option)-1 DO result = ray_psi2010_process(option[i])
END
;
; ======================================================================
;
PRO ray_psi2010_output, ps, option=option, nlevels=nlevels, pos=pos

  PSOPEN, filename = 'data_ray/' + ps + '.ps'

  IF (NOT KEYWORD_SET(option)) THEN  $
    option = [1, 10,11,12, 20,21,22,23, 30,31, 40,41,42,43,44,45, 50,51,52, 60]

  FOR i = 0, N_ELEMENTS(option)-1 DO ray_psi2010_plots,option[i],1,param1=pos, ps='on'

  PSCLOSE
END
;
; ======================================================================
;
