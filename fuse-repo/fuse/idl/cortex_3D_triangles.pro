
;
; ======================================================================
; from https://www.idlcoyote.com/tips/make_circle.html
;
FUNCTION cortex_PlotCircle, zcenter, ycenter, radius

  points = (2 * !PI / 49.0) * FINDGEN(50)
  x = MAKE_ARRAY(50,/FLOAT,VALUE=-0.0001)
  y = ycenter + radius * SIN( points )
  z = zcenter + radius * COS( points )

;print,y
;print,z

  PLOTS, x, y, z, /DATA, COLOR=TrueColor('Red'), /T3D

  RETURN, 0 ; TRANSPOSE([[x],[y]])
END

;
; ======================================================================
;
FUNCTION cortex_MapToWall, data_tor_R, data_tor_Z, wall_n, wall_x, wall_y, wall_s

  n = N_ELEMENTS(data_tor_R)
  data_pol_d = MAKE_ARRAY(n,/FLOAT,VALUE=1.0E+20)
  data_pol_s = data_pol_d
  
  FOR i = 0L, n-1 DO BEGIN
  
   if (i MOD 100000 EQ 0) then print,'  i',i,' of',n-1
  
   p1 = [ data_tor_R[i], data_tor_Z[i] ]

   FOR j = 0L, wall_n-2 DO BEGIN

     k = j + 1
     IF (k EQ wall_n) THEN k = 0

     l0 = [ wall_x[j], wall_y[j] ]
     l1 = [ wall_x[k], wall_y[k] ] 

     d = PNT_LINE(p1,l0,l1,pl,/INTERVAL)

     IF (d LT data_pol_d[i]) THEN BEGIN
       data_pol_d[i] = d

       IF (ABS(l0[0]-l1[0]) GT 1.0E-06) THEN s = (pl[0] - l0[0]) / (l1[0] - l0[0]) ELSE  $
                                             s = (pl[1] - l0[1]) / (l1[1] - l0[1])
       s = MAX([0.0,MIN([1.0,s])])
       s0 = wall_s[j  ]
       s1 = wall_s[j+1]
       data_pol_s[i] = s0 + s * (s1 - s0)
;        IF (i GE 26000 AND i LT 26100) THEN BEGIN
;          print,'  >  ', p1,l0,l1,j,data_pol_d[i],s,data_pol_s[i],FORMAT='(A,3(2F12.6,2X),I6,3F12.6)'
;        ENDIF
      ENDIF

    ENDFOR

;    IF (i GE 26000 AND i LT 26100) THEN  $
;      print,'data',i,p1,data_pol_d[i],data_pol_s[i],FORMAT='(A,I6,2F12.6,2X,2F12.6)'



  ENDFOR

  result = { data_pol_d : data_pol_d, data_pol_s : data_pol_s }

  RETURN, result

END
;
; ======================================================================
;
FUNCTION cortex_DOTP, v1, v2

  result = 0.0

  FOR i = 0, N_ELEMENTS(v1)-1 DO result = result + v1[i] * v2[i]

  RETURN, result
END
;
; ======================================================================
;
FUNCTION cortex_VectorLength, vector

  result = SQRT( vector[0]^2 + vector[1]^2 + vector[2]^2 )

  RETURN, result
END
;
; ======================================================================
;
FUNCTION cortex_RotateVector, vector, angles

  result_vector = vector

  IF (angles[0] NE 0.0) THEN BEGIN
    angle_rad = angles[0] * !PI / 180.0D
    matrix = [ [ 1.0D , 0.0D           ,  0.0D           ] ,  $
               [ 0.0D , COS(angle_rad) , -SIN(angle_rad) ] ,  $ 
               [ 0.0D , SIN(angle_rad) ,  COS(angle_rad) ] ]  
    result_vector = matrix ## result_vector
  ENDIF

  IF (angles[1] NE 0.0) THEN BEGIN
    angle_rad = angles[1] * !PI / 180.0
    matrix = [ [ COS(angle_rad) , 0.0 ,  SIN(angle_rad) ] ,  $
               [ 0.0            , 1.0 ,  0.0            ] ,  $
               [-SIN(angle_rad) , 0.0 ,  COS(angle_rad) ] ]
    result_vector = matrix ## result_vector
  ENDIF

  IF (angles[2] NE 0.0) THEN BEGIN
    angle_rad = angles[2] * !PI / 180.0
    matrix = [ [ COS(angle_rad) , -SIN(angle_rad) , 0.0 ] ,  $
               [ SIN(angle_rad) ,  COS(angle_rad) , 0.0 ] ,  $ 
               [ 0.0            ,  0.0            , 1.0 ] ]  
    result_vector = matrix ## result_vector
  ENDIF


  RETURN, result_vector

END
;
; ======================================================================
;
FUNCTION cortex_Plot3DTriangles, plot, plot_data, annotate=annotate, ps=ps, time_slice=time_slice

  window_id = 0
  window_xsize = 700
  window_ysize = 700

  dev_xsize = !D.X_SIZE
  dev_ysize = !D.Y_SIZE 
  !P.BACKGROUND = 255 ; TrueColor('White')

  title    = plot.title
  notes    = plot.notes
  charsize = plot.charsize

  scale_label = 'test'

  !P.CHARSIZE  = charsize
  !P.CHARTHICK = plot.thick
  !P.THICK     = plot.thick
  !X.THICK     = plot.thick
  !Y.THICK     = plot.thick
  !Z.THICK     = plot.thick

  color_table = plot.color_table

; Build a master list of polygons:

  ndata = N_ELEMENTS(TAG_NAMES(plot_data))
  IF (ndata LE 0) THEN BEGIN
    PRINT, 'ERROR cortex_Plot3DTriangles: No data found'
    RETURN, -1
  ENDIF

  print, '------------------------'
  print, 'ndata',ndata
  print, '------------------------'

  FOR idata = 1, ndata DO BEGIN
    data = cortex_ExtractStructure(plot_data,idata)
    CASE idata OF
;     ------------------------------------------------------------------
;
;     Triangle geometry.
;
      1: BEGIN

        CASE plot.option OF      
          3:
          ELSE: BEGIN
            data_x     = data.x
            data_y     = data.y
            data_z     = data.z
            data_v     = data.v
            data_i     = data.i
            data_n     = data.n
            data_val   = MAKE_ARRAY(data_n,/FLOAT,VALUE=-999.0)
            data_area  = MAKE_ARRAY(data_n,/FLOAT,VALUE=-999.0)
            data_type  = MAKE_ARRAY(data_n,/LONG ,VALUE=1     )
            data_iliin = MAKE_ARRAY(data_n,/LONG ,VALUE=1     )
            END
        ENDCASE



print,'index check ------------------ ',N_ELEMENTS(data_x),N_ELEMENTS(data_i)
help,data,/struct

        END
;     ------------------------------------------------------------------
;
;     Triangle data.
;
      2: BEGIN
        shift = N_ELEMENTS(data_x)
;       Selection of data:
        CASE plot.option OF      
          1: BEGIN
            data_select = -data.nt_par_mol_1_0
            i = WHERE(data.iliin EQ -2, n)
            END
          2: BEGIN
            data_select = data.nt_par_imp_1_0
            i = WHERE(data_select GT 0.0, n)
            END
          3: BEGIN

            scale_label = 'D2 particle flux (...)'
            data_select = data.nt_par_mol_1_0 ; 

;            scale_label = 'D2 heat flux (kW m-2)'
;            data_select = data.in_ene_mol_1_0 * 1.0E-3 ; kW

;            scale_label = 'D heat flux (kW m-2)'
;            data_select = data.in_ene_atm_1_0 * 1.0E-3 ; kW

            data_tor_R   = SQRT( data.x^2 + data.z^2 )
            data_tor_Z   = data.y
            data_tor_phi = ASIN( data.z / data_tor_R ) * 180.0 / !PI

;            i = WHERE( ABS(data_tor_phi[data.v[0,*]] - data_tor_phi[data.v[1,*]]) LT 1.0E-3 AND  $
;                       ABS(data_tor_phi[data.v[2,*]] - data_tor_phi[data.v[1,*]]) LT 1.0E-3 AND  $
;                       data_tor_phi[data.v[0,*]] GT 0.0 AND  $
;                       data_tor_phi[data.v[0,*]] LT 1.0 , n )


            i = WHERE(data.in_add EQ 82 AND data.iliin EQ 2, n)

;            i = WHERE(data.in_add EQ 79 OR data.in_add EQ 81, n)

;            i = WHERE(data.in_add EQ 14, n) ! linear
            print, 'size of i',n
print,MIN(data_tor_phi),MAX(data_tor_phi)
print,MIN(data_tor_phi[data.v[0,i]]),MAX(data_tor_phi[data.v[0,i]])
            END
          ELSE:
        ENDCASE
;       Selection of surface:
        IF (n LT 1) THEN BEGIN
          PRINT,'ERROR cortex_Plot3DTriangles: No flux accounting surfaces found'
        ENDIF

        IF (plot.option EQ 3) THEN BEGIN
          n_index = 3 * LINDGEN(n)
          data_x = REFORM(data.x[data.v[0:2,i]],3*n)
          data_y = REFORM(data.y[data.v[0:2,i]],3*n)
          data_z = REFORM(data.z[data.v[0:2,i]],3*n)
          data_v = TRANSPOSE([ [ n_index   ],  $
                               [ n_index+1 ],  $
                               [ n_index+2 ] ])
          data_i     = data.i     [i]
          data_val   = data_select[i]
          data_area  = data.area  [i]
          data_type  = REPLICATE(2,n)  
          data_iliin = data.iliin [i]  
          data_n = N_ELEMENTS(data_i)
          data_index = i
        ENDIF ELSE BEGIN
          n_index = 3 * LINDGEN(n) + N_ELEMENTS(data_x)
          data_x = [ data_x, REFORM(data.x[data.v[0:2,i]],3*n) ]
          data_y = [ data_y, REFORM(data.y[data.v[0:2,i]],3*n) ]
          data_z = [ data_z, REFORM(data.z[data.v[0:2,i]],3*n) ]
          data_v = TRANSPOSE([ [ REFORM(data_v[0,*]), n_index   ],  $
                               [ REFORM(data_v[1,*]), n_index+1 ],  $
                               [ REFORM(data_v[2,*]), n_index+2 ] ])
          data_i     = [ data_i    , data.i     [i] ]
          data_val   = [ data_val  , data_select[i] ]
          data_area  = [ data_area , data.area  [i] ]
          data_type  = [ data_type , REPLICATE(2,n) ]  
          data_iliin = [ data_iliin, data.iliin [i] ]  
          data_n = N_ELEMENTS(data_i)
        ENDELSE

        i = WHERE(data_val NE -999.0, count)
        net1 = 0.0
        net2 = 0.0
        net3 = 0.0
        IF (count GE 1) THEN BEGIN
          net1 = net1 + TOTAL(data_val[i]*data_area[i])
        ENDIF
        i = WHERE(data.iliin EQ -2, count)
        IF (count GT 0) THEN BEGIN
          net2 = TOTAL(data.nt_par_atm_1_0[i]*data.area[i])
          print, 'net',2.0*net1,net2,2.0*net1+net2,2.0*net1-net2
        ENDIF

print,'index check ------------------ ',N_ELEMENTS(data_x),N_ELEMENTS(data_i)
help,data,/struct
        END
;     ------------------------------------------------------------------
;
;     3D ionisation data.
;
      3: BEGIN

        CASE plot.option OF      
          3:
          ELSE: BEGIN
            shift = N_ELEMENTS(data_x)
;           Selection of data:
            CASE 2 OF      
              1: data_select = data.dens
              2: data_select = data.s_ion
              ELSE:
            ENDCASE

            i = WHERE(data.region EQ 1, n)
            IF (n LT 1) THEN BEGIN
              PRINT,'ERROR cortex_Plot3DTriangles: No magnetic grid volumes found'
            ENDIF

            n_index = LINDGEN(n) + N_ELEMENTS(data_x)

            data_x = [ data_x, data.x[i] ]
            data_y = [ data_y, data.y[i] ]
            data_z = [ data_z, data.z[i] ]

            data_v = TRANSPOSE([ [ REFORM(data_v[0,*]), n_index ],  $
                                 [ REFORM(data_v[1,*]), n_index ],  $
                                 [ REFORM(data_v[2,*]), n_index ] ])

            data_i    = [ data_i   , data.iobj  [i] ]
            data_val  = [ data_val , data_select[i] ]
            data_type = [ data_type, REPLICATE(3,n) ]

            data_n = N_ELEMENTS(data_i)

            print,'index check ------------------ ',N_ELEMENTS(data_x),N_ELEMENTS(data_i)
          END
        ENDCASE
        END
;     ------------------------------------------------------------------
;
;     Field line information.
;
      4: BEGIN

        CASE plot.option OF      
          3:
          ELSE: BEGIN

            n = N_ELEMENTS(data.v[*,0])

            line_index = data.index

help,data,/struct
print,'>>>>> n',n,N_ELEMENTS(line_index)
        
            data_type = [ data_type, REPLICATE(4,n) ]
            data_i    = [ data_i, N_ELEMENTS(data_x) + INDGEN(n) ]
            data_x    = [ data_x, data.v[*,0] ]
            data_y    = [ data_y, data.v[*,1] ]
            data_z    = [ data_z, data.v[*,2] ]

print,'n',n

print, 'whtt',N_ELEMENTS(data_type),N_ELEMENTS(data_x),N_ELEMENTS(data_i)

;            print,data.v[*,0]
            print,'index check ------------------ ',N_ELEMENTS(data_x)
            print,'data_i',data_i[N_ELEMENTS(data_i)-100:N_ELEMENTS(data_i)-1]

;            print, i
;            print,data_type[i]

          END
        ENDCASE
        END
;     ------------------------------------------------------------------
      ELSE: BEGIN
        PRINT,'ERROR cortex_Plot3DTriangles: Unrecognized data set'
        RETURN,-1
        END
    ENDCASE
  ENDFOR

  print,'calculating data_tor'
  IF (1 EQ 1) THEN BEGIN
    data_tor_R   = SQRT( data_x^2 + data_z^2 )
    data_tor_Z   = data_y
    data_tor_phi = ASIN( data_z / data_tor_R ) * 180.0 / !PI


;    i = WHERE(data_tor_phi LT 0.0, count)
;    IF (count GT 0) THEN data_tor_phi[i] += 360.0
  ENDIF

;print,data_tor_phi[data_i[i]]
;stop

  print,'done'
;
; ----------------------------------------------------------------------
;
  IF (KEYWORD_SET(annotate)) THEN BEGIN
    help,annotate.data1,/struct

    print,'calculating data_pol'

    wall_x = annotate.data1.x
    wall_y = annotate.data1.y
    wall_n = N_ELEMENTS(wall_x)
    wall_s = MAKE_ARRAY(wall_n,/FLOAT,VALUE=0.0)

    FOR j = wall_n-2, 0, -1 DO BEGIN
      k = j + 1
      wall_s[j] = SQRT( (wall_x[j]-wall_x[k])^2 + (wall_y[j]-wall_y[k])^2 ) + wall_s[j+1]
    ENDFOR

print,'wall_s length',wall_s

    IF (1 eq 1 and plot.quick EQ 1) THEN BEGIN
      RESTORE, FILENAME='output/'+STRTRIM(plot.case_name[0],2)+'_2.sav'
    ENDIF ELSE BEGIN
      result = cortex_MapToWall(data_tor_R, data_tor_Z, wall_n, wall_x, wall_y, wall_s)
      data_pol_d = result.data_pol_d
      data_pol_s = result.data_pol_s
      cortex_Undefine, result
      SAVE, FILENAME='output/'+STRTRIM(plot.case_name[0],2)+'_2.sav', data_pol_d, data_pol_s
    ENDELSE
    
    print,'  data_pol_c'
    
;    data_pol_c = (data_tor_phi / 360.0) * (2.0 * !PI * data_tor_R) ; 0.1)


    data_pol_c = (data_tor_phi / 360.0) * (2.0 * !PI * 6.0 ) ; 0.02 ) ; 7.0) ; 0.1)  ! *** gas mod hack *** 


    print,'data_pol_d',MIN(data_pol_d),MAX(data_pol_d)
    print,'data_pol_s',MIN(data_pol_s),MAX(data_pol_s)
    print,'data_pol_c',MIN(data_pol_c),MAX(data_pol_c)
    
    print,'done'
  ENDIF

;
; Stretch out the field lines in circumference, rather than having them circle back on their nasty selves due to their toroidal origins:
; ----------------------------------------------------------------------
;
  i = WHERE(data_type eq 4, n)

  IF (1 EQ 1 AND n GT 0) THEN BEGIN

   line_z       = data_z      [data_i[i]]
   line_tor_phi = data_tor_phi[data_i[i]]
   line_pol_c   = data_pol_c  [data_i[i]]

;print,i

print,N_ELEMENTS(data_z),N_ELEMENTS(line_index), N_ELEMENTS(line_z)
print,data_z[N_ELEMENTS(data_z)-10:N_ELEMENTS(data_z)-1]
print,line_z[N_ELEMENTS(line_z)-10:N_ELEMENTS(line_z)-1]
print,line_index[N_ELEMENTS(line_index)-10:N_ELEMENTS(line_index)-1]

;print,data_z[i]

    FOR line = 1, 2 DO BEGIN

      j = WHERE(line_index EQ line, count)

      IF (count EQ 0) THEN CONTINUE

      line_z_j       = line_z      [j] 
      line_tor_phi_j = line_tor_phi[j]
      line_pol_c_j   = line_pol_c  [j]

      k = WHERE(line_z_j EQ 0.0, count)
      
      print, '<><><><><><><><><><> count',count,line,N_ELEMENTS(i)

 print,line_z_j[N_ELEMENTS(line_z_j)-10:N_ELEMENTS(line_z_j)-1]

;rint,line_z_j

      n = N_ELEMENTS(line_z_j)
      
      FOR l = k[0]+2, n-1 DO BEGIN
        IF (line_pol_c_j[l] GT line_pol_c_j[l-1]) THEN BEGIN
          delta = line_pol_c_j[l-1] - line_pol_c_j[l]
          line_pol_c_j[l:n-1] = line_pol_c_j[l:n-1] + 2 * delta
        ENDIF
      ENDFOR

      FOR l = k[0]-2, 0, -1 DO BEGIN
        IF (line_pol_c_j[l] LT line_pol_c_j[l+1]) THEN BEGIN
          delta = line_pol_c_j[l] - line_pol_c_j[l+1]
          line_pol_c_j[0:l] = line_pol_c_j[0:l] - 2 * delta
        ENDIF
      ENDFOR

      line_z      [j] = line_z_j      
      line_tor_phi[j] = line_tor_phi_j
      line_pol_c  [j] = line_pol_c_j  

    ENDFOR

    data_pol_c[data_i[i]] = line_pol_c

    print,'---------',data_pol_c[N_ELEMENTS(data_pol_c)-10:N_ELEMENTS(data_pol_c)-1]

  ENDIF

;
; Ionisation over plot:
; ----------------------------------------------------------------------
;
  i = WHERE(data_type EQ 3, n)
 
  IF (1 EQ 1 AND n GT 0) THEN BEGIN

    print,'>>>>>>>>>>>>> ionisation plot <<<<<<<<<<<<<<<<'

    j = REFORM(data_v[0,i])

help,j

;    print,'n',n
;    print,'data_n',data_n
;    print,'i',i[0:10],i[N_ELEMENTS(i)-10:N_ELEMENTS(i)-1]

;    print,'data.x',data.x[  0:10 ],data.x[N_ELEMENTS(data.x)-10:N_ELEMENTS(data.x)-1]
;    print,'data.y',data.y[  0:10 ],data.y[N_ELEMENTS(data.x)-10:N_ELEMENTS(data.x)-1]
;    print,'data.z',data.z[  0:10 ],data.z[N_ELEMENTS(data.x)-10:N_ELEMENTS(data.x)-1]

;    print,'data_x',data_x[j[0:10]],data_x[j[N_ELEMENTS(j)-10:N_ELEMENTS(j)-1]]
;    print,'data_y',data_y[j[0:10]],data_y[j[N_ELEMENTS(j)-10:N_ELEMENTS(j)-1]]
;    print,'data_z',data_z[j[0:10]],data_z[j[N_ELEMENTS(j)-10:N_ELEMENTS(j)-1]]

;    print,'data_pol_d',data_pol_d[j[0:10]],data_pol_d[j[N_ELEMENTS(j)-10:N_ELEMENTS(j)-1]]
;    print,'data_pol_s',data_pol_s[j[0:10]],data_pol_s[j[N_ELEMENTS(j)-10:N_ELEMENTS(j)-1]]
;    print,'data_pol_c',data_pol_c[j[0:10]],data_pol_c[j[N_ELEMENTS(j)-10:N_ELEMENTS(j)-1]]

    n_points = 100

    x = data_x  [j] / MAX(data_x  [j])
    y = data_y  [j] / MAX(data_y  [j])
    z = data_z  [j] / MAX(data_z  [j])
    f = data_val[i] / MAX(data_val[i])

print,'min/max_f',MIN(f),MAX(f)
help,x
help,y
help,z
help,result
print,MIN(x),MAX(x)
print,MIN(y),MAX(y)
print,MIN(z),MAX(z)

    IF (1 EQ 1 AND plot.quick EQ 1) THEN BEGIN

      RESTORE, FILENAME='output/'+STRTRIM(plot.case_name[0],2)+'_3.sav'

    ENDIF ELSE BEGIN

      print,'qhull'
      QHULL, x, y, z, tet, /DELAUNAY 
      print,'qgrid3'
      volume = QGRID3(x,y,z,f,tet,DIMENSION=n_points)

      SAVE, FILENAME='output/'+STRTRIM(plot.case_name[0],2)+'_3.sav', volume

    ENDELSE


    help,volume,/struct


    print,'min/max d',MIN(x),MAX(x)
    print,'min/max s',MIN(y),MAX(y)
    print,'min/max c',MIN(z),MAX(z)
  
    print,'xvolume'
loadct,0
;    XVOLUME, BYTSCL(volume)

    print,'min/max d',MIN(x),MAX(x)
    print,'min/max s',MIN(y),MAX(y)
    print,'min/max c',MIN(z),MAX(z)

    print,'done'

  ENDIF
;
; 
; ----------------------------------------------------------------------
;
  geometry_mode = 2

  CASE geometry_mode OF
    1: BEGIN  ; "realistic"
      x = data_x
      y = data_y
      z = data_z
      END
    2: BEGIN  ; flat
      x = data_pol_d ; + 0.1
      y = data_pol_s
      z = data_pol_c
      END
    ELSE: BEGIN
      PRINT,'ERROR cortex_Plot3DTriangles: Unrecognized geometry_mode'
      STOP
      END    
  ENDCASE

; *** hack ***
  angle1 = [  0.0, 80.0, 0.0]
  angle2 = [  0.0,  0.0,  0.0]
  angle3 = [  0.0,   0.0,   0.0]
;  angle1 = [  0.0,   0.0,   0.0]
;  angle2 = [  0.0,   0.0,   0.0]
;  angle3 = [  0.0,   0.0,   0.0]

;
; 
; ----------------------------------------------------------------------
;
  print,'calculating visibility'

  view_vector = [ 0.0 , 0.0 , 10.0D ]
  view_vector = cortex_RotateVector(view_vector,-angle3)
  view_vector = cortex_RotateVector(view_vector,-angle2)
  view_vector = cortex_RotateVector(view_vector,-angle1)

  norm = MAKE_ARRAY(3,data_n,/FLOAT              )
  dist = MAKE_ARRAY(  data_n,/FLOAT,VALUE=-1.0E+6)
  list = MAKE_ARRAY(  data_n,/LONG ,VALUE=-1     )
  n = -1L
  FOR i = 0L, data_n-1 DO BEGIN

    IF (data_type[i] EQ 3 OR data_type[i] EQ 4) THEN CONTINUE

    ; Decide if the surface is visible by calculating the surface normal and checking whether or not the
    ; angle between it and the view vector, which is pointing into your eye, is less than 180 degrees:
    ivtx1 = data_v[0,i]
    ivtx2 = data_v[1,i]
    ivtx3 = data_v[2,i]

    v1 = [ x[ivtx1] - x[ivtx2], y[ivtx1] - y[ivtx2], z[ivtx1] - z[ivtx2] ]
    v2 = [ x[ivtx3] - x[ivtx2], y[ivtx3] - y[ivtx2], z[ivtx3] - z[ivtx2] ]

    norm[*,i] = CROSSP(v1,v2) 
    dot_product = cortex_DOTP(view_vector,norm[*,i])

; *** hack ***     IF (dot_product LT 0.0) THEN CONTINUE 

    ; Calculate the distance from the centre of the surface to the location of your eye:
    x1 = TOTAL(x[data_v[0:2,i]]) / 3.0
    y1 = TOTAL(y[data_v[0:2,i]]) / 3.0
    z1 = TOTAL(z[data_v[0:2,i]]) / 3.0      
    n++
    dist[n] = cortex_VectorLength([x1,y1,z1]-view_vector)
    list[n] = i
  ENDFOR
  n++
  ; Sort the order in which the surfaces will be drawn, based on how far they are from your eye:
  order = REVERSE(SORT(dist))

help,dist
print,n,data_n

  dist_min = MIN(dist[0:n-1])
  dist_max = MAX(dist[0:n-1])

  print,'done'

  length_view_vector = cortex_VectorLength(view_vector)

  light_vector = [1.0,0.0,1.0]
  length_light_vector = cortex_VectorLength(light_vector)
;
; ----------------------------------------------------------------------
;
  xdelta = 1.0
  ydelta = 1.0

  IF (KEYWORD_SET(ps)) THEN BEGIN
    xsize = 297.0  ; Landscape A4
    ysize = 210.0  
    display_ratio = 4.2/4.4 ; 4.45 / 5.2  ; Small correction based on measurements from the GV display
    aspect_ratio = (xdelta / xsize) / (ydelta / ysize) * display_ratio
  ENDIF ELSE BEGIN
    PRINT, 'NEED TO FIX WINDOW ASPECT RATIO'
    RETURN, -1
    display_ratio = 4.4
    aspect_ratio = xdelta / ydelta / display_ratio *            $ 
                   (FLOAT(window_ysize) / FLOAT(window_xsize))
    WINDOW, window_id, XSIZE=window_xsize, YSIZE=window_ysize
  ENDELSE
  IF (plot.aspect_ratio GT 0.0) THEN aspect_ratio = plot.aspect_ratio

  T3D, /RESET  
;  T3D, ROTATE=angles
  T3D, ROTATE=angle1
  T3D, ROTATE=angle2
  T3D, ROTATE=angle3
; *** hack ***
  T3D, SCALE    =[0.9,0.9,0.9]
  T3D, TRANSLATE=[0.05,0.075,0.0]
;  T3D, SCALE    =[0.6,0.6,0.6]  ; 3D gas puff from upper port area
;  T3D, TRANSLATE=[0.2,0.4,0.0]

;  T3D, PERSPECTIVE=5.0

  CASE geometry_mode OF
    1: BEGIN

      delta_x = MAX(x) - MIN(x)
      delta_y = MAX(y) - MIN(y)
      delta_z = MAX(z) - MIN(z)

      mid_x = 0.5 * (MAX(x) + MIN(x))
      mid_y = 0.5 * (MAX(y) + MIN(y))
      mid_z = 0.5 * (MAX(z) + MIN(z))

      delta = 0.5 * MAX([delta_x,delta_y,delta_z])

      xrange = [ mid_x - delta / aspect_ratio, mid_x + delta / aspect_ratio ]
      yrange = [ mid_y - delta, mid_y + delta ]
      zrange = [ mid_z - delta, mid_z + delta ]

; *** hack ***      
;      xrange = [      MIN(y),     MAX(y)]
;      yrange = [      MIN(y),     MAX(y)]
;      zrange = [-0.75*MAX(y),0.75*MAX(y)]
      END
    2: BEGIN
      xrange = [     MIN(y),    MAX(y)]
      yrange = [     MIN(y),    MAX(y)]
      zrange = [  -0.5*(MAX(y)-MIN(y))/aspect_ratio,0.5*(MAX(y)-MIN(y))/aspect_ratio]
;      zrange = [-0.71*MAX(y),0.71*MAX(y)]

print,'xrange',xrange
print,'yrange',yrange
print,'zrange',zrange

      END
    ELSE: BEGIN
      PRINT,'ERROR cortex_Plot3DTriangles: Unrecognized geometry_mode'
      STOP
      END    
  ENDCASE


;  SURFACE, [[0.0,0.0],[0.0,0.0]], /NODATA,  $
;           XRANGE=xrange, YRANGE=yrange, ZRANGE=zrange,  $
;           XSTYLE=1     , YSTYLE=1     , ZSTYLE=1     ,  $
;           /T3D ; , BACKGROUND=Truecolor('Green')
;;           POSITION=[0.1, 0.1, 1.9, 1.9, -1.0, 1.0]
;;           AX=30.0, AZ=30.0 , /SAVE ; /T3D ; /SAVE
; *** hack ***
  SURFACE, [[0.0,0.0],[0.0,0.0]], /NODATA,  $
           XRANGE=xrange, YRANGE=yrange, ZRANGE=zrange,  $
;           XSTYLE=1     , YSTYLE=1     , ZSTYLE=1     ,  $
   XTITLE='x', $
   YTITLE='y', $
   ZTITLE='z', $
           XSTYLE=5     , YSTYLE=5     , ZSTYLE=5     ,  $
           /T3D ; , BACKGROUND=Truecolor('Green')
;           POSITION=[0.1, 0.1, 1.9, 1.9, -1.0, 1.0]
;           AX=30.0, AZ=30.0 , /SAVE ; /T3D ; /SAVE



; Establish range of values:
  minval =  1.0E+20
  maxval = -1.0E+20
  FOR j = 0L, n-1 DO BEGIN
    i = list[order[j]]
    IF (data_type[i] NE 2) THEN CONTINUE
    minval = MIN([minval,data_val[i]])
    maxval = MAX([maxval,data_val[i]])
  ENDFOR
;print,data_val

; minval = -1.0E+16
; maxval = -1.0E+15


;  minval = 0.0
;  maxval = maxval * 0.01

minval = minval * 0.0
maxval = maxval * 0.00001

; *** gas mod hack ***  maxval = maxval * 0.01
; maxval = maxval * 1.0

; maxval = maxval * 10.0
;  maxval = maxval * 0.1

  print,'min/max ----------------------------',minval,maxval

  print,'plotting surfaces'
   
;
; Setup shading.
;

  lastct = -1
  frac = MAKE_ARRAY(n,/FLOAT,VALUE=0.0)
  FOR l = 0L, n-1 DO BEGIN
    i = list[order[l]]

    IF (data_type[i] EQ 3 OR data_type[i] EQ 4) THEN CONTINUE

    IF (data_val[i] NE -999.0) THEN BEGIN

      frac1 = (data_val[i] - minval) / (maxval - minval)
      frac1 = MAX([frac1,0.0])
      frac1 = MIN([frac1,1.0])
      IF (plot.log AND data_val[i] EQ 0.0) THEN frac1 = 0.0

      frac[l] = frac1

    ENDIF ELSE BEGIN
      ; Get the angle between the surface normal and the light vector: 

      dot_product = ABS(cortex_DOTP(light_vector,norm[*,i]))
      light_angle = ACOS(dot_product / length_light_vector / cortex_VectorLength(norm[*,i])) * 180.0 / !PI
      frac1 = (light_angle / 180.0)^1.0
      frac2 = ((dist[order[l]] - dist_min) / (dist_max - dist_min) ) *0.3 + 0.7

      frac[l] = 1.0 - frac1 * frac2

    ENDELSE
  ENDFOR

  max_x = -1E+20

;
; ----------------------------------------------------------------------
; Draw the structural surfaces.
;
  IF (1 EQ 1 AND plot.option NE 3) THEN BEGIN

    i = where(data_type EQ 1) 

    print,'x',MIN(x[i]),MAX(x[i])
    print,'y',MIN(y[i]),MAX(y[i])
    print,'z',MIN(z[i]),MAX(z[i])

    max_y = 0.0

    LOADCT, 0, /SILENT
    FOR l = 0L, n-1 DO BEGIN
      i = list[order[l]]

      IF (data_type[i] NE      1) THEN CONTINUE
      IF (data_val [i] NE -999.0) THEN CONTINUE

      IF (geometry_mode EQ 2) THEN BEGIN
        x_check = MAX(x[data_v[0:2,i]])
        IF (x_check GT 0.01) THEN CONTINUE
      ENDIF

      max_y = MAX([max_y,MAX(y[data_v[0:2,i]])])

      max_x = MAX([max_x,MAX(x[data_v[0:2,i]])])

      frac1 = (frac[l] - MIN(frac)) / (MAX(frac) - MIN(frac))

      color = FIX( (1.0 - frac1) * 250.0) 

    ; print,'light',light_angle,frac1,frac2,color,dist_min,dist_max

       POLYFILL, x[data_v[0:2,i]], y[data_v[0:2,i]], z[data_v[0:2,i]], /DATA, NOCLIP=1, /T3D, COLOR=color ;  COLOR=TrueColor('Red')
    ENDFOR

;oplot,[-1.0,1.0],[12.0,12.0],color=TrueColor('Red')
;       PLOTS, [0.0,0.0],[12.8,12.8],[-1.0,1.0], /DATA, COLOR=TrueColor('Red'), /T3D
;       PLOTS, [0.0,0.0],[ 4.8, 4.8],[-1.0,1.0], /DATA, COLOR=TrueColor('Red'), /T3D

;rint,'max_y',max_y

  ENDIF

;
; ----------------------------------------------------------------------
; Draw the transparent surfaces.
;
  dump_data = 1

  IF (dump_data EQ 1) THEN BEGIN

    fp = 3
    FREE_LUN,fp        
    file_name = 'cortex_data/triangle_dump.txt'
    OPENW,fp,file_name, ERROR=err
    IF (err NE 0) THEN BEGIN
      PRINT,'ERROR cortex_Plot3DTriangles: Problem opening data stream (to file)'
      PRINT,'FILE_NAME= ',file_name
      PRINTF,-2,!err.msg
      FREE_LUN,fp
      STOP
    ENDIF    

    count = 0
    FOR l = 0L, n-1 DO BEGIN
      i = list[order[l]]
      IF (data_type[i] NE      2) THEN CONTINUE
      IF (data_val [i] EQ -999.0) THEN CONTINUE
      count++
    ENDFOR

    PRINTF,fp,'*'
    PRINTF,fp,count
  ENDIF


  FOR l = 0L, n-1 DO BEGIN

    i = list[order[l]]

    IF (data_type[i] NE      2) THEN CONTINUE
    IF (data_val [i] EQ -999.0) THEN CONTINUE

    IF (geometry_mode EQ 2) THEN BEGIN
      x_check = MAX(x[data_v[0:2,i]])
      IF (x_check GT 0.01) THEN CONTINUE
    ENDIF

    max_x = MAX([max_x,MAX(x[data_v[0:2,i]])])

    frac1 = frac[l]
    IF (plot.log AND data_val[i] EQ 0.0) THEN frac1 = 0.0
    color = cortex_SetColour(xrange,frac1,lastct)
    LOADCT, color_table, /SILENT

    POLYFILL, x[data_v[0:2,i]], y[data_v[0:2,i]], z[data_v[0:2,i]], /DATA, NOCLIP=1, /T3D, COLOR=color ;  COLOR=TrueColor('Red')

    IF (dump_data EQ 1) THEN BEGIN
 
      ivtx1 = data_v[0,i]
      ivtx2 = data_v[1,i]
      ivtx3 = data_v[2,i]

      v1 = [ data_x[ivtx1] - data_x[ivtx2], data_y[ivtx1] - data_y[ivtx2], data_z[ivtx1] - data_z[ivtx2] ]
      v2 = [ data_x[ivtx3] - data_x[ivtx2], data_y[ivtx3] - data_y[ivtx2], data_z[ivtx3] - data_z[ivtx2] ]

      norm = CROSSP(v1,v2) 
      norm = norm / cortex_VectorLength(norm)

      PRINTF,fp,i,  $
        MEAN(data_x[data_v[0:2,i]]),  $
        MEAN(data_y[data_v[0:2,i]]),  $
        MEAN(data_z[data_v[0:2,i]]),  $
        norm,  $
        data_area[i],  $
        data_val [i],  $
        data_x[data_v[0:2,i]],data_y[data_v[0:2,i]],data_z[data_v[0:2,i]],FORMAT='(I9,2(3F10.6,5X),2E10.2,5X,3(3F10.6))'

    ENDIF


;    Draw surface normals:
;    PLOTS, [data_x[ivtx2],data_x[ivtx2]+normal[0]*100.0],  $
;           [data_y[ivtx2],data_y[ivtx2]+normal[1]*100.0],  $
;           [data_z[ivtx2],data_z[ivtx2]+normal[2]*100.0],  $
;           /DATA, COLOR=TrueColor('Red'), /T3D
  ENDFOR

;  PLOTS, [0.0,view_vector[0]],  $
;         [0.0,view_vector[1]],  $
;         [0.0,view_vector[2]],  $
;         /DATA, COLOR=TrueColor('Blue'), /T3D, THICK=3

  print,'done',max_x

  IF (dump_data EQ 1) THEN CLOSE, fp

;
; ----------------------------------------------------------------------
; Draw the magnetic field lines.
;

  i = WHERE(data_type EQ 4, n)

  IF (0 EQ 1 AND n GT 0) THEN BEGIN

    line_x = x[ data_i[i] ]
    line_y = y[ data_i[i] ]
    line_z = z[ data_i[i] ]

    FOR line = 1, 2 DO BEGIN

      j = WHERE(line_index EQ line)  

      print, '> FIELD LINE >>>>>>>>>>>>>>',line,N_ELEMENTS(j)

      yshift = 0.0

      FOR k = 0, 0 DO BEGIN

;print,line_x[j]
print,line_y[j]
;print,line_z[j]

        PLOTS, line_x[j], line_y[j] + yshift, line_z[j], /DATA, COLOR=TrueColor('Red'), /T3D
        yshift = yshift + 0.50 

      ENDFOR

    ENDFOR

  ENDIF

;
; ----------------------------------------------------------------------
; Calculate radius for injection of gas.
;
  i = WHERE(data_type EQ 2, n)

  IF (0 EQ 1 AND n GE 1) THEN BEGIN

    x_cen = MAKE_ARRAY(n,/FLOAT,VALUE=0.0) & y_cen = x_cen & z_cen = x_cen
    FOR j = 0, n-1 DO BEGIN
      x_cen[j] = MEAN( x[data_v[0:2,i[j]]] )
      y_cen[j] = MEAN( y[data_v[0:2,i[j]]] )
      z_cen[j] = MEAN( z[data_v[0:2,i[j]]] )
    ENDFOR

    max_val = MAX(data_val[i],i_max)

    flux_val = data_val[i] * data_area[i]
    flux_tot = TOTAL(flux_val)

    print,'total',2.0*flux_tot

    radius = SQRT( (x_cen-x_cen[i_max])^2 + (y_cen-y_cen[i_max])^2 + (z_cen-z_cen[i_max])^2 )

    dist_end  = 5.1
    dist_step = 0.01
    k = -1 & dist = 0.0
    WHILE (1) DO BEGIN
      k++
      dist += dist_step
      IF (dist GT dist_end) THEN BREAK
      j = WHERE( radius LT dist, count)
      flux_total = TOTAL(flux_val[j])  
      IF (count GE 1) THEN BEGIN
         IF (k EQ 0) THEN BEGIN
           flux_sum  = flux_total
           flux_dist = dist
         ENDIF ELSE BEGIN
           IF (flux_total NE flux_sum[N_ELEMENTS(flux_sum-1)-1]) THEN BEGIN
             flux_sum = [flux_sum , flux_total]
             flux_dist =[flux_dist, dist      ]
           ENDIF
         ENDELSE
      ENDIF
    ENDWHILE

    cont = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]

    cont_dist = INTERPOL(flux_dist,flux_sum/flux_tot,cont)

    fp1 = 3
    FREE_LUN,fp1        
    file_name = 'output/' + STRTRIM(plot.case_name[0],2) + '.idl.tet_puff_radius'
    OPENW,fp1,file_name, ERROR=err
    IF (err NE 0) THEN BEGIN
      PRINT,'ERROR cortex_Plot3DTrianlges: Problem opening data stream (to file)'
      PRINT,'  FILE_NAME= ',file_name
      PRINTF,-2,!err.msg
      STOP
    ENDIF
    FOR i = 0, N_ELEMENTS(cont)-1 DO  $
      PRINTF,fp1,cont[i],cont_dist[i],FORMAT='(2F12.4)'
    CLOSE,fp1
    FREE_LUN,fp1

    FOR j = 0, N_ELEMENTS(cont)-1 DO BEGIN

      k = WHERE( radius LT cont_dist[j], count)
      IF (count GE 1) THEN BEGIN
         flux_sum = TOTAL(flux_val[k])
         print, 'within',flux_sum,flux_sum/flux_tot
      ENDIF
    
      dummy = cortex_PlotCircle(z_cen[i_max], y_cen[i_max], cont_dist[j])

    ENDFOR

    print,'maxes',z_cen[i_max], y_cen[i_max]

  ENDIF






;
; Draw axes:
; ----------------------------------------------------------------------
;
;  AXIS, /XAxis, 0, 0, 0, /T3D, CharSize=1.5, Color=TrueColor('Black'), XSTYLE=1, XRANGE=xrange
;  AXIS, /YAxis, 0, 0, 0, /T3D, CharSize=1.5, Color=TrueColor('Black'), YSTYLE=1, YRANGE=yrange
;  AXIS, /ZAxis, 0, 0, 0, /T3D, CharSize=1.5, Color=TrueColor('Black'), ZSTYLE=1, ZRANGE=zrange
;
; Setup plot area:
; ----------------------------------------------------------------------
;
  size = 0.7 ; 0.8
  IF (plot.size NE 0.0) THEN size = plot.size

;  T3D, SCALE    =[1.0,1.0,1.0]
;  T3D, TRANSLATE=[0.0,0.0,0.0]

  IF (N_ELEMENTS(WHERE(plot.center EQ 0.0)) EQ 2) THEN BEGIN 
    xcen = 0.5 * TOTAL(plot.frame_bnds) 
    ycen = 0.5
  ENDIF ELSE BEGIN
    IF (plot.center[0] NE 0.0) THEN xcen = plot.center[0]  $
                               ELSE xcen = 0.5 * TOTAL(plot.frame_bnds)
    IF (plot.center[1] NE 0.0) THEN ycen = plot.center[1]  $
                               ELSE ycen = 0.5
  ENDELSE  
;
; Setup plot zoom:
; ----------------------------------------------------------------------
;
  IF (N_ELEMENTS(WHERE(plot.zoom EQ 0.0)) EQ 4) THEN BEGIN 
    xmin = 0.0
    xmax = 1.0
    ymin = 0.0
    ymax = 1.0
  ENDIF ELSE BEGIN
    xmin = plot.zoom[0]
    xmax = plot.zoom[3]
    ymin = plot.zoom[1]
    ymax = plot.zoom[4]
    zmin = plot.zoom[2]
    zmax = plot.zoom[5]
  ENDELSE
  xdelta = xmax - xmin
  ydelta = ymax - ymin
;
; ----------------------------------------------------------------------
;
  IF (KEYWORD_SET(ps)) THEN BEGIN
    xsize = 297.0  ; Landscape A4
    ysize = 210.0  
    display_ratio = 3.25 / 3.35  ; Small correction based on measurements from the GV display

    display_ratio = 0.5 ; *** HACK ***

    aspect_ratio = (xdelta / xsize) / (ydelta / ysize) * display_ratio
  ENDIF ELSE BEGIN
    PRINT, 'NEED TO FIX WINDOW ASPECT RATIO'
    RETURN, -1
    display_ratio = 4.4
    aspect_ratio = xdelta / ydelta / display_ratio *            $ 
                   (FLOAT(window_ysize) / FLOAT(window_xsize))
    WINDOW, window_id, XSIZE=window_xsize, YSIZE=window_ysize
  ENDELSE
  IF (plot.aspect_ratio GT 0.0) THEN aspect_ratio = plot.aspect_ratio
;
;
  IF ((xdelta / xsize) LT (ydelta / ysize)) THEN BEGIN
    xpos = [xcen - 0.7 * size * aspect_ratio ,  $
            xcen + 0.7 * size * aspect_ratio ]
    ypos = [ycen - 0.5 * size, ycen + 0.5 * size]
  ENDIF ELSE BEGIN
    xpos = [xcen - 0.5 * size, xcen + 0.5 * size]
    ypos = [ycen - 0.5 * size / aspect_ratio,  $
            ycen + 0.5 * size / aspect_ratio]
  ENDELSE
;
; Check if the plot fits in the allocated portion of the page, and scale the plot
; as necessary:
  IF (NOT plot.rigid AND (xpos[0] LT plot.frame_bnds[0] OR xpos[1] GT plot.frame_bnds[1])) THEN BEGIN
    ratio = 0.8 * (plot.frame_bnds[1] - plot.frame_bnds[0]) / (xpos[1] - xpos[0]) 
    xpos = [xcen - ratio * 0.5 * (xpos[1] - xpos[0]),  $
            xcen + ratio * 0.5 * (xpos[1] - xpos[0])]
    ypos = [ycen - ratio * 0.5 * (ypos[1] - ypos[0]),  $
            ycen + ratio * 0.5 * (ypos[1] - ypos[0])]
  ENDIF
;
; Plot the scale: 
; ----------------------------------------------------------------------
  nsteps = 100.0

  IF (aspect_ratio LE 1.0) THEN BEGIN
    x0 = (0.95 * xpos[0] + 0.05 * xpos[1]) * dev_xsize
    x1 = (0.05 * xpos[0] + 0.95 * xpos[1]) * dev_xsize
    y0 = (ypos[0] - 0.11 * (ypos[1] - ypos[0])) * dev_ysize
    y1 = (ypos[0] - 0.09 * (ypos[1] - ypos[0])) * dev_ysize
;    y0 = (ypos[0] - 0.15 * (ypos[1] - ypos[0])) * dev_ysize
;    y1 = (ypos[0] - 0.13 * (ypos[1] - ypos[0])) * dev_ysize
;    y0 = (0.10 * ypos[0]                 ) * dev_ysize
;    y1 = (0.30 * ypos[0]                 ) * dev_ysize
    LOADCT, color_table
    lastct = -1
    FOR i = 0, nsteps-1 DO BEGIN
      frac1 = FLOAT(i  ) / nsteps
      frac2 = FLOAT(i+1) / nsteps
      x = [x0 + frac1 * (x1 - x0), x0 + frac2 * (x1 - x0),  $
           x0 + frac2 * (x1 - x0), x0 + frac1 * (x1 - x0)]
      y = [y0, y0, y1 ,y1]
      color = cortex_SetColour(xrange,0.5 * (frac1 + frac2),lastct)
      POLYFILL, x, y, /DEVICE, COLOR=color
    ENDFOR
;   This is nasty...
    log_prefix = ''
    IF (plot.log NE 0) THEN log_prefix = 'LOG10 '
    PLOT, [minval,maxval], [0.0,0.0], POSITION=[x0,y0,x1,y1], YTICKFORMAT='(A1)',  $
          TICKLEN=0, /NODATA, /DEVICE, /NOERASE,  $
          XTITLE=log_prefix+scale_label,COLOR=TrueColor('Black'), XSTYLE=1
  ENDIF ELSE BEGIN
    x0 = (xpos[1] + 0.09 * (xpos[1] - xpos[0])) * dev_xsize
    x1 = (xpos[1] + 0.11 * (xpos[1] - xpos[0])) * dev_xsize
    y0 = (0.95 * ypos[0] + 0.05 * ypos[1]) * dev_ysize
    y1 = (0.05 * ypos[0] + 0.95 * ypos[1]) * dev_ysize
    LOADCT, color_table
    lastct = -1
    FOR i = 0, nsteps-1 DO BEGIN
      frac1 = FLOAT(i  ) / nsteps
      frac2 = FLOAT(i+1) / nsteps
      x = [x0, x0, x1 ,x1]
      y = [y0 + frac1 * (y1 - y0), y0 + frac2 * (y1 - y0),  $
           y0 + frac2 * (y1 - y0), y0 + frac1 * (y1 - y0)]
      color = cortex_SetColour(xrange,0.5 * (frac1 + frac2),lastct)
      POLYFILL, x, y, /DEVICE, COLOR=color
    ENDFOR
    log_prefix = ''
    IF (plot.log NE 0) THEN log_prefix = 'LOG10 '
    PLOT, [0.0,0.0], xrange, POSITION=[x0,y0,x1,y1], XTICKFORMAT='(A1)',  $
          TICKLEN=0, /NODATA, /DEVICE, /NOERASE,  $
          YTITLE=log_prefix+scale_label,COLOR=TrueColor('Black'), YSTYLE=1
;    OPLOT, [x0,x0,x1,x1], [y0,y1,y1,y0], COLOR=TrueColor('Black')
  ENDELSE
;


;
;
;
;
  IF (1 EQ 1 AND plot.option EQ 3) THEN BEGIN

    i = WHERE(data_type EQ 2, n)

    print, 'size of i',n,N_ELEMENTS(data_type)

    i = WHERE( data_pol_c[data_v[0,*]] GT -0.05 AND  $
               data_pol_c[data_v[0,*]] LT  0.05 AND  $
               data_pol_c[data_v[1,*]] GT -0.05 AND  $
               data_pol_c[data_v[1,*]] LT  0.05 AND  $
               data_pol_c[data_v[2,*]] GT -0.05 AND  $
               data_pol_c[data_v[2,*]] LT  0.05 , n )

    print, 'size of i',n

    s = REFORM((data_pol_s[data_v[0,i]] + data_pol_s[data_v[1,i]] + data_pol_s[data_v[2,i]]) / 3.0,n)
    v = data_val[i]

    
    print,'v' ,MIN(v),MAX(v)
    print,'s' ,MIN(data_pol_s[data_v[0,i]]),MAX(data_pol_s[data_v[0,i]])
    print,'c1',MIN(data_pol_c[data_v[0,i]]),MAX(data_pol_c[data_v[0,i]])
    print,'c2',MIN(data_pol_c[data_v[1,i]]),MAX(data_pol_c[data_v[1,i]])
    print,'c2',MIN(data_pol_c[data_v[2,i]]),MAX(data_pol_c[data_v[2,i]])


    j = SORT(s)
    xtitle = 'poloidal distance along vacuum vessel surface from cassette-blanket gap (m)'
    ytitle = 'poloidal profile along central cassette gap (kW)'
    PLOT,s[j],v[j],xstyle=1,ystyle=1,xtitle=xtitle,ytitle=ytitle,title=scale_label    


    val = 5.812
    del = 0.2
    bnd1 = val-del
    bnd2 = val+del
    i = WHERE( data_pol_s[data_v[0,*]] GT bnd1 AND  $
               data_pol_s[data_v[0,*]] LT bnd2 AND  $
               data_pol_s[data_v[1,*]] GT bnd1 AND  $
               data_pol_s[data_v[1,*]] LT bnd2 AND  $
               data_pol_s[data_v[2,*]] GT bnd1 AND  $
               data_pol_s[data_v[2,*]] LT bnd2 , n )

    print, 'size of i',n

    s = REFORM((data_pol_c[data_v[0,i]] +  $
                data_pol_c[data_v[1,i]] +  $
                data_pol_c[data_v[2,i]]) / 3.0, n)
    v = data_val[i]
    j = SORT(s)
    xtitle = '(approximate) toroidal distance from central cassette gap (m)'
    ytitle = 'toroidal profile behind inner strike-point (kW)'
    PLOT,s[j],v[j],xstyle=1,ystyle=1,xtitle=xtitle,ytitle=ytitle,title=scale_label    


    val = 2.53
    del = 0.2
    bnd1 = val-del
    bnd2 = val+del
    i = WHERE( data_pol_s[data_v[0,*]] GT bnd1 AND  $
               data_pol_s[data_v[0,*]] LT bnd2 AND  $
               data_pol_s[data_v[1,*]] GT bnd1 AND  $
               data_pol_s[data_v[1,*]] LT bnd2 AND  $
               data_pol_s[data_v[2,*]] GT bnd1 AND  $
               data_pol_s[data_v[2,*]] LT bnd2 , n )

    print, 'size of i',n

    s = REFORM((data_pol_c[data_v[0,i]] +  $
                data_pol_c[data_v[1,i]] +  $
                data_pol_c[data_v[2,i]]) / 3.0, n)
    v = data_val[i]
    j = SORT(s)
    xtitle = '(approximate) toroidal distance from central cassette gap (m)'
    ytitle = 'toroidal profile behind outer strike-point (kW)'
    PLOT,s[j],v[j],xstyle=1,ystyle=1,xtitle=xtitle,ytitle=ytitle,title=scale_label    



;    val = 0.383
;    del = 0.20
;    bnd1 = val-del
;    bnd2 = val+del
;    i = WHERE( data_pol_s[data_v[0,*]] GT bnd1 AND  $
;               data_pol_s[data_v[0,*]] LT bnd2 AND  $
;               data_pol_s[data_v[1,*]] GT bnd1 AND  $
;               data_pol_s[data_v[1,*]] LT bnd2 AND  $
;               data_pol_s[data_v[2,*]] GT bnd1 AND  $
;               data_pol_s[data_v[2,*]] LT bnd2 , n )
;    print, 'size of i',n


;                data_pol_c[data_v[1,i]] +  $
;                data_pol_c[data_v[2,i]]) / 3.0, n)
;    v = data_val[i]
;    j = SORT(s)
;    xtitle = '(approximate) toroidal distance from central cassette gap (m)'
;    ytitle = 'toroidal profile behind outer cassette-blanket gap (kW)'
;    PLOT,s[j],v[j],xstyle=1,ystyle=1,xtitle=xtitle,ytitle=ytitle,title=scale_label    



  ENDIF







  RETURN, 0

END
;
; ======================================================================
;
