; Example:
; ts=get_ts(shots=[15165,15169],times=[0.180,0.250],fp=0.9,xrange=[.8,1.3],yrange=[3.E+19,200.])
; ts=get_ts(shots=[15024],times=[0.210,0.240],fp=[0.90,0.80,1.10,0.985],xrange=[.7,1.4],yrange=[3.E+19,100.])
;
;
;
; 
; ts=get_ts(shots=[15021,17466,17469],times=[0.220,0.280],xrange=[.0,1.3],yrange=[6.E+19,1200.],/nofit)
;
;
; From RS code
;
PRO GFUNC_EDGE,x,a,f,dfda

        f = EDGEFUNCTIONATS(x,a)

        dfda = EDGEFUNCTIONATS(x,a,/DA)

END

;
;
;

FUNCTION get_ts, shots=shots,         $ ; List of shots
                 times=times,         $ ; Time interval over which to collect data
                 xrange=xrange,       $ ; PSIn x-axis range for plots
                 yrange=yrange,       $ ; Element 1: max n_e, 2: max Te
                 fp=fp,               $ ; Transition point between core poly fit and edge/SOL exp fit (default = 0.9)
                 edge=edge,           $ ; /EDGE to plot only data from edgeTS system, and not the standard YAG
                 noplot=noplot,       $ ; /NOPLOT to turn off plots
                 output = output,     $ ; /OUTPUT to turn on some screen text output
                 nofit=nofit,         $ ; /NOFIT to stop the fits
                 nocolour=nocolour,   $ ; /NOCOLOUR to turn off standard colour coding on plots
                 ps=ps                  ; /PS send output to postscipt file 




  MAXPTS = 10000
  struct ={                               $
           nfp         :  0            ,  $  ; Number of PSIn fit transition points (2 for L-mode, 3 for H-mode)
           fp          : fltarr(4)     ,  $  ; Transition points between poly (core), tanh (edge), exp (SOL) fits
           n           : -1            ,  $  ; Number of data points (not yet assigned)
           shots       : intarr(50)    ,  $  ; List of shots where data was requested
           times       : fltarr(2)     ,  $  ; Time interval over which data was requested
           shot        : intarr(MAXPTS),  $  ; Shot number
           pname       : intarr(MAXPTS),  $  ; "Polyname", used to tag data from the EDGE system
           t           : fltarr(MAXPTS),  $  ; Time (s)
           r           : fltarr(MAXPTS),  $  ; R (m)
           z           : fltarr(MAXPTS),  $  ; Z (m) = 0.0 
           psin        : fltarr(MAXPTS),  $  ; Normalized magnetic flux coordinate
           n_e         : fltarr(MAXPTS),  $  ; Electron density
           n_e_e       : fltarr(MAXPTS),  $  ; Electron density error
           n_e_f_psin  : fltarr(MAXPTS),  $  ; PSIn values for fit to n_e data
           n_e_f       : fltarr(MAXPTS),  $  ; Fit to n_e data
           n_e_f0      : fltarr(MAXPTS),  $  ; Debug - temporary
           n_e_f1      : fltarr(MAXPTS),  $  ; Debug
           n_e_A       : fltarr(3)     ,  $  ; Exponential n_e fit parameters (for PSIn > fp[nfp-1])
           n_e_TH      : fltarr(5)     ,  $  ; TANH fit to H-mode pedestal (for fp[0] < PSIn < fp[1])
           n_e_C       : fltarr(6)     ,  $  ; Polynomial n_e fit parameters (for PSIn < fp[0])
           te          : fltarr(MAXPTS),  $  ; Te
           te_e        : fltarr(MAXPTS),  $  ; Te error
           te_f_psin   : fltarr(MAXPTS),  $  ; PSIn values for fit to Te data
           te_f        : fltarr(MAXPTS),  $  ; Fit to Te data
           te_f0       : fltarr(MAXPTS),  $  ; Debug - temporary
           te_f1       : fltarr(MAXPTS),  $  ; Debug
           te_A        : fltarr(3)     ,  $  ; Exponential n_e fit parameters (for PSIn > fp[nfp-1])
           te_TH       : fltarr(5)     ,  $  ; TANH fit to H-mode pedestal (for fp[0] < PSIn < fp[1])
           te_C        : fltarr(6)      }    ; Polynomial n_e fit parameters (for PSIn < fp[0])


  ts  = struct

  IF (NOT KEYWORD_SET(fp)) THEN fp = 0.9

  nfp = N_ELEMENTS(fp)

  ts.shots = shots
  ts.times = times
  ts.nfp   = nfp
  ts.fp    = fp


; Collect data
  FOR is = 0, N_ELEMENTS(shots)-1 DO BEGIN

    shot = shots[is]

    atm = atm_ida_read(shot)
    time1 = times[0]
    time2 = times[1]

    i = WHERE(atm.t.n_e GE time1 AND atm.t.n_e LE time2)  ; Correct for Te as well?


;  left off , need to cycle  if no data

  
    IF (KEYWORD_SET(edge)) THEN BEGIN
      j = WHERE(atm.polyname GT 200)
    ENDIF ELSE BEGIN
      j = WHERE(atm.polyname GT -1)
    ENDELSE

    n = N_ELEMENTS(j)

    nstart = N_ELEMENTS(WHERE(ts.t NE 0.0))  ; A bit arbitrary

    FOR iburst = 0, N_ELEMENTS(i)-1 DO BEGIN
      i1 = nstart +  iburst    * n
      i2 = nstart + (iburst+1) * n - 1

      ts.shot [i1:i2] = shot
      ts.t    [i1:i2] = atm.t.n_e   [  i[iburst]] 
      ts.r    [i1:i2] = atm.r2      [j,i[iburst]]
      ts.z    [i1:i2] = 0.0
      ts.n_e  [i1:i2] = atm.n_e     [j,i[iburst]] 
      ts.n_e_e[i1:i2] = atm.d.n_e   [j,i[iburst]] 
      ts.te   [i1:i2] = atm.t_e     [j,i[iburst]]
      ts.te_e [i1:i2] = atm.d.t_e   [j,i[iburst]]
      ts.pname[i1:i2] = atm.polyname[j]
    ENDFOR

    nfinish = N_ELEMENTS(WHERE(ts.t NE 0.0)) - 1

;   Get PSIn:
    ts.psin[nstart:nfinish] =  $
           get_psin(shot,ts.r[nstart:nfinish],  $
                         ts.z[nstart:nfinish],  $
                         ts.t[nstart:nfinish],/MAST)
  ENDFOR

; Fitting:
  IF (NOT KEYWORD_SET(nofit)) THEN BEGIN
;	edge_poly_goodness = where( atm.polyname gt 200 AND atm.n_e gt 5e17)
    
    FitProfiles, ts, 1
    FitProfiles, ts, 2

    print,'Density fit coeff:'
    print,'C   :',ts.n_e_C
    PRINT,'TANH:',[ts.n_e_TH,0.0]
    print,'A   :',[ts.n_e_A,0.0,0.0,0.0]
    print,'Te fit coeff:'
    print,'C   :',ts.te_C
    PRINT,'TANH:',[ts.te_TH,0.0]
    print,'A   :',[ts.te_A,0.0,0.0,0.0]


  ENDIF

; Plots:
  IF (NOT KEYWORD_SET(noplot)) THEN BEGIN
    PlotProfiles, ts, shots, xrange, yrange, nofit, nocolour, ps
  ENDIF

; Output:
  IF (KEYWORD_SET(output)) THEN BEGIN
    psin = MAKE_ARRAY(20,/FLOAT,value=0.0)
    n_e  = MAKE_ARRAY(20,/FLOAT,value=0.0)
    te   = MAKE_ARRAY(20,/FLOAT,value=0.0)
    i = WHERE(ts.te_f_psin NE 0.0) 
    max_psin = ts.te_f_psin[N_ELEMENTS(i)-1]
    print,max_psin
    FOR i = 0, 19 DO BEGIN
      psin1 = FLOAT(i) / FLOAT(19) * max_psin
      IF (psin1 LT fp[0]) THEN BEGIN
        n_e[i] = POLY(psin1,ts.n_e_C)        
        te [i] = POLY(psin1,ts.te_C )        
      ENDIF ELSE BEGIN
        IF (psin1 GT fp[nfp-1]) THEN BEGIN
          A = ts.n_e_A
          n_e[i] = A[0] * EXP( (psin1-fp[nfp-1]) * A[1]) + A[2]        
          A = ts.te_A
          te [i] = A[0] * EXP( (psin1-fp[nfp-1]) * A[1]) + A[2]        
        ENDIF ELSE BEGIN
        ENDELSE
      ENDELSE
      print,psin1,n_e[i],te[i]
    ENDFOR
  ENDIF  

return, ts

END
;
; ======================================================================
;
PRO PlotProfiles, ts, shots, xrange, yrange, nofit, nocolour, ps


  IF (KEYWORD_SET(ps)) THEN BEGIN
    !P.MULTI = 0
;    !P.MULTI = [0, 1, 2]
    SET_PLOT, 'PS'
    DEVICE, FILE=ps+'.ps',/encapsulated,/color,bits=8
    !P.FONT=0
    loadct,4 ; 5
  ENDIF ELSE BEGIN
    !P.MULTI = [0, 1, 2]
    safe_colors,/first             
  ENDELSE

  i = WHERE(ts.psin NE 0.0) 

  IF (KEYWORD_SET(xrange)) THEN BEGIN
    min_x = xrange[0]
    max_x = xrange[1]
  ENDIF ELSE BEGIN
    min_x = MIN(ts.psin[i])
    max_x = MAX(ts.psin[i])
  ENDELSE
;
; Density:
;  GOTO, JUMP_TEMPERATURE

  IF (KEYWORD_SET(yrange)) THEN BEGIN
    min_y = 0.0
    max_y = yrange[0]
  ENDIF ELSE BEGIN
    i = WHERE(ts.psin GE min_x AND ts.psin LE max_x)
    min_y = MIN(ts.n_e[i])
    max_y = MAX(ts.n_e[i])
  ENDELSE
  plot, [min_x,max_x], [min_y,max_y], /NODATA, title='ne', xstyle=1, ystyle=1
  FOR is = 0, N_ELEMENTS(shots)-1 DO BEGIN
    psym_colour = is + 2
    xyouts,xco(0.8),yco(0.5+0.05*is),STRING(shots(is)),col=psym_colour  ; From Rory, bless
    IF (KEYWORD_SET(nocolour)) THEN psym_colour = 1
    i = WHERE(ts.shot EQ shots[is] AND FINITE(ts.n_e) EQ 1 AND ts.pname LE 200) 
    print,'WTF:',N_ELEMENTS(i)
    IF (N_ELEMENTS(i) GT 1) THEN oplot, ts.psin[i], ts.n_e[i], psym=5, color=psym_colour
    i = WHERE(ts.shot EQ shots[is] AND FINITE(ts.n_e) EQ 1 AND ts.pname GT 200) 
    oplot, ts.psin[i], ts.n_e[i], psym=6, color=psym_colour
  ENDFOR

  IF (NOT KEYWORD_SET(nofit)) THEN BEGIN

    data_x = FINDGEN(200) / 199.0 * (max_x - min_x) + min_x

    ts.n_e_f_psin = data_x
;
;   Core polynomial fit:
    i = WHERE(data_x LE ts.fp[0]) 
    IF (i[0] NE -1) THEN BEGIN
      C = ts.n_e_C
      data_y = POLY(data_x[i],C)
      oplot, data_x[i], data_y
      ts.n_e_f[i] = data_y
    ENDIF
;
;   Edge TANH fit:
    IF (ts.nfp EQ 4) THEN BEGIN
      a0 = ts.n_e_TH[0]  ; position
      a1 = ts.n_e_TH[1]  ; width
      a2 = ts.n_e_TH[2]  ; height
      a3 = ts.n_e_TH[3]  ; slope
      a4 = ts.n_e_TH[4]  ; offset

      i = WHERE(data_x GE ts.fp[1] AND data_x LE ts.fp[2]) 
      x = data_x[i] 
      ;- function (from /home/mastts/lib/edgefunctionats.pro)
      z = (a0 - x)/(2.*a1)
      data_y = (a2-a4)/2.*(MTANH(z,a3) + 1) + a4
      oplot, data_x[i], data_y*1.0E+19, color=2

      i = WHERE(data_x GE ts.fp[0] AND data_x LE ts.fp[3]) 
      x = data_x[i] 
      z = (a0 - x)/(2.*a1)
      data_y = (a2-a4)/2.*(MTANH(z,a3) + 1) + a4
      oplot, data_x[i], data_y*1.0E+19

      ts.n_e_f[i] = data_y*1.0E+19

    ENDIF
;
;   SOL Exponential fit:
    i = WHERE(data_x GE ts.fp[ts.nfp-1]) 
    A = ts.n_e_A
    data_y = A[0] * EXP(A[1] * (data_x[i]-ts.fp[ts.nfp-1])) + A[2]  
    oplot, data_x[i], data_y

    ts.n_e_f[i] = data_y

    IF (NOT KEYWORD_SET(ps)) THEN BEGIN
      i = WHERE(ts.n_e_f0 NE 0.0)
      oplot, ts.n_e_f0[i], ts.n_e_f1[i], psym = 3
    ENDIF
  ENDIF
;
; Temperature:
;  GOTO, JUMP_END

  JUMP_TEMPERATURE:

  IF (KEYWORD_SET(yrange)) THEN BEGIN
    min_y = 0.0
    max_y = yrange[1]
  ENDIF ELSE BEGIN
    i = WHERE(ts.psin GE min_x AND ts.psin LE max_x)
    min_y = MIN(ts.te[i])
    max_y = MAX(ts.te[i])
  ENDELSE
  plot, [min_x,max_x], [min_y,max_y], /NODATA, title='Te', xstyle=1, ystyle=1
  FOR is = 0, N_ELEMENTS(shots)-1 DO BEGIN
    psym_colour = is + 2
    IF (KEYWORD_SET(nocolour)) THEN psym_colour = 1
    i = WHERE(ts.shot EQ shots[is] AND FINITE(ts.te) EQ 1 AND ts.pname LE 200)
    IF (N_ELEMENTS(i) GT 1) THEN oplot, ts.psin[i], ts.te[i], psym=5, color=psym_colour
    i = WHERE(ts.shot EQ shots[is] AND FINITE(ts.te) EQ 1 AND ts.pname GT 200)
    oplot, ts.psin[i], ts.te[i], psym=6, color=psym_colour
  ENDFOR

  IF (NOT KEYWORD_SET(nofit)) THEN BEGIN

    data_x = FINDGEN(200) / 199.0 * (max_x - min_x) + min_x

    ts.te_f_psin = data_x

    i = WHERE(data_x LE ts.fp[0]) 
    IF (i[0] NE -1) THEN BEGIN
      C = ts.te_C
      data_y = POLY(data_x[i],C)
      oplot, data_x[i], data_y
      ts.te_f[i] = data_y
    ENDIF

    IF (ts.nfp EQ 4) THEN BEGIN
      a0 = ts.te_TH[0]  ; position
      a1 = ts.te_TH[1]  ; width
      a2 = ts.te_TH[2]  ; height
      a3 = ts.te_TH[3]  ; slope
      a4 = ts.te_TH[4]  ; offset

      i = WHERE(data_x GE ts.fp[1] AND data_x LE ts.fp[2]) 
      x = data_x[i] 
      z = (a0 - x)/(2.*a1)
      data_y = (a2-a4)/2.*(MTANH(z,a3) + 1) + a4
      oplot, data_x[i], data_y, color=2

      
      i = WHERE(data_x GE ts.fp[0] AND data_x LE ts.fp[3]) 
      x = data_x[i] 
      z = (a0 - x)/(2.*a1)
      data_y = (a2-a4)/2.*(MTANH(z,a3) + 1) + a4
      oplot, data_x[i], data_y

      ts.te_f[i] = data_y

    ENDIF

    i = WHERE(data_x GE ts.fp[ts.nfp-1]) 
    A = ts.te_A
    data_y = A[0] * EXP(A[1] * (data_x[i]-ts.fp[ts.nfp-1])) + A[2]  
    oplot, data_x[i], data_y

    ts.te_f[i] = data_y

    IF (NOT KEYWORD_SET(ps)) THEN BEGIN
      i = WHERE(ts.te_f0 NE 0.0)
      oplot, ts.te_f0[i], ts.te_f1[i], psym = 3
    ENDIF
  ENDIF

  JUMP_END:
  IF (KEYWORD_SET(ps)) THEN BEGIN
    DEVICE, /CLOSE  
    SET_PLOT, 'X'  
  ENDIF 

  !P.MULTI = 0  
END
;
; ======================================================================
;
PRO FitProfiles, ts, mode 
; ts temporary hopefully, for debugging at the moment...

  TH = 0

  fp = ts.fp
  nfp = ts.nfp
  xdata = ts.psin
  pname = ts.pname

  CASE mode OF
    1: BEGIN  ; ne
       ydata = ts.n_e
       END
    2: BEGIN  ; Te
       ydata = ts.te
       END
  ENDCASE


  CASE 2 OF  
    0:
    1: BEGIN
       END
    2: BEGIN

       FOR loop = 0, 2 DO BEGIN
;
;        ---------------------------------------------------------------
;        TANH fit to edge:
         IF (loop EQ 0 AND nfp EQ 4) THEN BEGIN

           i = WHERE(xdata GE fp[1] AND xdata LE fp[2] AND  $
                     FINITE(ydata) EQ 1 AND pname GT 200)

           psin = xdata[i] 

           center = 0.5 * (fp[1] + fp[2])

           CASE (1) OF
             1: BEGIN
                CASE mode OF
                  1: BEGIN  ; Ne
     	             ;      center, width, pedestal value, core_slope, sol value
                     fact = 1.0E+19
                     data = ydata[i] / fact
                     weights = 1.0 / (ts.n_e_e[i] / fact)^2
                     guess = [center, 0.05, 3.3  , -0.01, 0.1]
                     END
                  2: BEGIN  ; Te
                     fact = 1.0
                     data = ydata[i] 
                     weights = 1.0 / ts.te_e[i]^2
                     guess = [center, 0.05, 100.0, 0.3  , 5.0]
                     END
                ENDCASE

                print,"STARTING GUESS:",guess

                yfit = LMCURVEFIT(psin,data,weights,guess,da,  $
                                  NCHISQ=chisq,                $
                                  FUNCTION_NAME='gfunc_edge',  $
                                  TOL=1.0E-5,                  $
                                  COVARIANCE=covariance,       $
                                  ERROR=error)

                TH = guess

                print,"MODE :",mode
                print,"PSIN :",psin
                print,"DATA :",data
                print,"WEIGHTS :",weights
                print,"YFIT :",yfit
                print,"GUESS:",guess

                i = SORT(psin)
                psin2 = psin[i]
                data2 = yfit[i]

                END
           ENDCASE

         ENDIF
;
;        ---------------------------------------------------------------
;        Exponential fit to SOL:
         IF (loop EQ 0) THEN BEGIN 
           i = WHERE(xdata GE fp[nfp-1] AND FINITE(ydata) EQ 1)
           n = N_ELEMENTS(i)

           IF (nfp EQ 4) THEN BEGIN
             dummy= MIN(ABS(psin2-fp[3]),imin)
             GFUNC_EDGE,fp[3],TH,val
             psin = MAKE_ARRAY(3*n,/FLOAT,value=psin2[imin])
             data = MAKE_ARRAY(3*n,/FLOAT,value=val*fact)          
             psin[n:2*n-1] = 1.2 ; 1.35
             CASE mode OF 
               1: data[n:2*n+1] = 1.0E+17
               2: data[n:2*n+1] = 1.0
             ENDCASE
             weights = MAKE_ARRAY(3*n,/DOUBLE,value=1.0) 
             psin[0:n-1] = xdata[i]
             data[0:n-1] = ydata[i]
           ENDIF ELSE BEGIN
             psin = MAKE_ARRAY(2*n,/FLOAT,value=1.35)
             CASE mode OF 
               1: data = MAKE_ARRAY(2*n,/FLOAT,value=1.0E+17)
               2: data = MAKE_ARRAY(2*n,/FLOAT,value=1.0)
             ENDCASE
             psin[0:n-1] = xdata[i]
             data[0:n-1] = ydata[i]
             weights = MAKE_ARRAY(2*n,/DOUBLE,value=1.0) 
           ENDELSE

           CASE mode OF
             1: BEGIN  ; fitting ne
;                weights[0:n-1] = (0.2-(xdata[i]-DOUBLE(fp[nfp-1])))^5
;                weights[0:n-1] = 2.0
;                weights1 = 1.0 / ts.n_e_e[i]^2
;                weights[0:n-1  ] = weights1[0:n-1]
;                weights[n:2*n-1] = MAX(0.1*weights1)
                END
             2: BEGIN  ; fitting Te
                weights1 = 1.0 / ts.te_e[i]^2
                weights[0:n-1] = weights1[0:n-1]
                IF (nfp EQ 4) THEN weights[n:3*n-1] = MAX(weights1)
;                IF (nfp EQ 4) THEN weights[n:2*n-1] = MAX(weights1)
                i = WHERE(psin LE 1.0 OR (psin GT 1.0 AND data LT 75.0))
                psin = psin[i]
                data = data[i]
                weights = weights[i]
                END
           ENDCASE
         ENDIF ELSE BEGIN
           i = WHERE(ABS(data4-data3) LT 2.*data4)
           psin = psin3[i]
           data = data3[i]
           weights = weights3[i]
         ENDELSE
     
         j = SORT(psin)
         CASE mode OF
           1: A = DOUBLE([1.0E+19, -10., 1.E+18])  ; ne
           2: A = DOUBLE([50.0   , -10., 5.0   ])  ; Te
         ENDCASE
         psin3    = psin[j] 
         data3    = data[j]  ; Store ydata used in fit for descrimination on subsequent loop pass
         weights3 = weights[j]
         data4 = MPCURVEFIT(psin[j]-DOUBLE(fp[nfp-1]), data[j], weights[j], A,  $
                            sigma, FUNCTION_NAME='gfunct_ts',status=status) 

;        For debugging:
         CASE mode OF
           1: BEGIN
              ts.n_e_f0 = 0.0
              ts.n_e_f1 = 0.0
              ts.n_e_f0 = psin3
              ts.n_e_f1 = data3
              END
           2: BEGIN
              ts.te_f0 = 0.0
              ts.te_f1 = 0.0
              ts.te_f0 = psin3
              ts.te_f1 = data3
              END
         ENDCASE
;
;        ---------------------------------------------------------------
;        Polynomial fit to core:
         i = WHERE(xdata GT 0.0 AND xdata LE fp[0] AND FINITE(ydata) EQ 1)
         n = N_ELEMENTS(i)
         IF (nfp EQ 4) THEN BEGIN
           dummy= MIN(ABS(psin2-fp[0]),imin)
           psin = MAKE_ARRAY(2*n,/FLOAT,value=psin2[imin])
           data = MAKE_ARRAY(2*n,/FLOAT,value=data2[imin]*fact)          
         ENDIF ELSE BEGIN
           psin = MAKE_ARRAY(2*n,/FLOAT,value=psin3[0])
           data = MAKE_ARRAY(2*n,/FLOAT,value=data4[0])
         ENDELSE
         psin[0:n-1] = xdata[i]
         data[0:n-1] = ydata[i]
         IF (nfp EQ 4) THEN BEGIN
           CASE mode OF
             1: degree = 2
             2: degree = 3
           ENDCASE
         ENDIF ELSE BEGIN
           CASE mode OF
             1: degree = 3 ; 2
             2: degree = 4 ; 3
           ENDCASE
         ENDELSE
         C = POLY_FIT(psin, data, degree, status=status)

         i = WHERE(psin LT fp[0])
         psin = psin[i]
         data = data[i]
         j = SORT(psin)
         psin1 = psin[j] 
         data1 = POLY(psin[j],C)

         print,status
       ENDFOR 

       END
  ENDCASE  


; Output:
  CASE mode OF
    1: BEGIN
       ts.n_e_A  = A
       ts.n_e_C  = C
       ts.n_e_TH = TH
       END
    2: BEGIN
       ts.te_A = A
       ts.te_C = C
       ts.te_TH = TH
       END
  ENDCASE

END
;
; ----------------------------------------------------------------------
;
PRO gfunct_ts, X, A, F, pder  

  bx = EXP(A[1] * X)  

  F = A[0] * bx + A[2]  
;  F = A[0] * bx   
  
;If the procedure is called with four parameters, calculate the  
;partial derivatives.  
  IF N_PARAMS() GE 4 THEN $  
    pder = [[bx], [A[0] * X * bx], [replicate(1.0, N_ELEMENTS(X))]]  
END  
