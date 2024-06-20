;
; ======================================================================
;
; ======================================================================
;
FUNCTION grid_AddCoreRings, c_array, wall, b, debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON commons, param, option, state
  ; ------------------------------------------------------------------

  tags  = STRUPCASE(TAG_NAMES(c_array))
  nctr  = N_ELEMENTS(tags)

  PRINT,'---------- ADDING param.CORE RINGS ------------'

  dist = grid_RadialDistribution(c_array, wall, b, param.CORE, debug=debug, xrange=xrange, yrange=yrange)


help,dist,/struct
print,dist.psi
print,b.psi_1st_xpoint
;stop

  core_n = N_ELEMENTS(dist.dist)      ; number of core rings to add
;  core_n = 3      ; number of core rings to add

  psi   = b.psi ; psi_raw
  psi_x = b.x
  psi_y = b.y

;  help,c_array,/struct   

;  psi_val = b.psi_1st_xpoint * 0.99999D

  core_i = 0

  FOR ictr = nctr+1, nctr+core_n DO BEGIN

    cont = 1

    psi_val = dist.psi[core_i] ; core_n - 1 - core_i]

      print,'================================================'
      print,'psi_val =', psi_val,core_i

;    core_step_dist = 0.01   ; spatial separatrion (m) between rings, at the outer midplane
;    core_step_psi  = 0.010  ; step in PSI to start with  
;    psi_val = psi_val + core_step_psi 
 
    count = 0
;
;   --------------------------------------------------------------------
;   GET CONTOUR FROM PSI MAP AND SEARCH FOR THE SEGMENT OF INTEREST 
;  
    WHILE cont DO BEGIN

      count = count + 1

      ctr = grid_ExtractContour(psi, psi_x, psi_y, psi_val)

      ibrk = WHERE(ctr.dist GT 0.10, cbrk)  ; PARAMETER
      nseg = N_ELEMENTS(ibrk) 
      ibrk = [-1,ibrk,ctr.n-1]

;        PLOT,ctr.x,ctr.y,color=Truecolor('Orange'),  $
;             XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE


      print,'nseg =',nseg

;
;     ------------------------------------------------------------------
;     LOOP OVER INDIVIDUAL CONTOUR SEGMENTS AND EXAMINE
;  
      FOR iseg = 0, nseg-1 DO BEGIN

        print,'----------------------------------------'
        print,'contour seg:',iseg,ibrk[iseg]+1,ibrk[iseg+1]

        ; Not a real segment so skip it, may need to strengthen this check:
        IF (cbrk GT 0) THEN BEGIN
          IF (ibrk[iseg]+2 EQ ibrk[iseg+1] OR  $
              ibrk[iseg]+1 EQ ibrk[iseg+1] OR  $
              ibrk[iseg]   EQ ibrk[iseg+1]) THEN CONTINUE  

          x = ctr.x[ibrk[iseg]+1:ibrk[iseg+1]]
          y = ctr.y[ibrk[iseg]+1:ibrk[iseg+1]]
        ENDIF ELSE BEGIN
          x = ctr.x
          y = ctr.y
        ENDELSE

        PLOT,[x,x[0]],[y,y[0]],color=Truecolor('Orange'),  $
             XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

;        print, 'psi_val = ',psi_val
;
;       ------------------------------------------------------------------
;       TEST IF THE CONTOUR IS IN THE RIGHT PLACE
;  
        status = 1  ; No real check performed yet, obviously -- should just make sure it's inside the sep core bit I guess... just looking at points in polygon...

        IF (status) THEN BEGIN
;
;       ------------------------------------------------------------------
;       IF YES, ADD THE CONTOUR TO THE LIST AND EXIT THE LOOPS
;  

          print,'adding...'

          core_i++
          IF (core_i EQ 1     ) THEN map_in  = [-1,-1] ELSE map_in  = [ictr - 1, -1]
          IF (core_i EQ core_n) THEN map_out = [ 1,-1] ELSE map_out = [ictr + 1, -1] 
          ctr = {  $
            state      : 0          ,  $
            origin     : 0          ,  $
            separatrix : 0          ,  $
            region     : param.CORE ,  $
            psi        : psi_val    ,  $
            map_in     : map_in     ,  $
            map_out    : map_out    ,  $
            x          : x          ,  $
            y          : y       }

          name = 'contour' + STRING(ictr,FORMAT='(I0)')

          c_array = CREATE_STRUCT(c_array,name,ctr)

          cont = 0
          BREAK
        ENDIF ELSE BEGIN
;
;       ------------------------------------------------------------------
;       IF NO, MODIFY THE SEARCH
;  

        ENDELSE

      ENDFOR

      IF (count EQ 5) THEN stop

      cont = 0

    ENDWHILE

  ENDFOR
;
; ----------------------------------------------------------------------
; ----------------------------------------------------------------------
; SETUP THE GRID POINTS
;  
; ----------------------------------------------------------------------
; COLLECT REQUIRED INFORMATION FROM THE PRIMARY SEPARATRIX
;  
  ctr = grid_ExtractStructure(c_array,'contour1')      

;  help,ctr,/struct

  IF (state.geometry EQ param.LIMITED) THEN BEGIN
    i = [0,N_ELEMENTS(ctr.grid_x)-1]
;    print,ctr.grid_x[i[0]],ctr.grid_y[i[0]]
;    print,ctr.grid_x[i[1]],ctr.grid_y[i[1]]
    xpt_x = ctr.grid_x[0]
    xpt_y = ctr.grid_y[0]
  ENDIF ELSE BEGIN
    i = WHERE(ctr.grid_x EQ ctr.null_x AND ctr.grid_y EQ ctr.null_y, count)
    print,i
    print,ctr.grid_s[i]
    IF (count NE 2) THEN BEGIN
      PRINT,'ERROR grid_AddCoreRings: x-point not properly identified'
      STOP
    ENDIF
    xpt_x = ctr.null_x
    xpt_y = ctr.null_y
  ENDELSE
;
; ----------------------------------------------------------------------
; DETERMINE THE DISTRIBUTION OF POINTS ON THE CORE PORTION OF SEPARATRIX
;
  sep_x = ctr.grid_x[i[0]:i[1]]  ; a bit crude at low resolution
  sep_y = ctr.grid_y[i[0]:i[1]]
;
; ----------------------------------------------------------------------
; LOOP OVER THE CONTOURS
;
  FOR ictr = nctr+1, nctr+core_n DO BEGIN

    name = 'contour' + STRING(ictr,FORMAT='(I0)')

    ctr = grid_ExtractStructure(c_array,name)      
;
;   --------------------------------------------------------------------
;   IDENTIFY THE START OF THE CORE RING (JUST THE CLOSEST POINT FOR NOW)
;
    x = ctr.x
    y = ctr.y

    ; Set the origin of the contour to the midplane to help refinement
    ; code (next bit):
    dummy = MIN(ABS(y), i)    
    x = [x[i:N_ELEMENTS(x)-1],x[0:i-1]]
    y = [y[i:N_ELEMENTS(y)-1],y[0:i-1]]

    ; Identify the point on the contour that's closest to the primary
    ; separatrix, then refine the contour there and search again:
    FOR j = 0, 1 DO BEGIN
      dummy = MIN( SQRT( (x - xpt_x)^2 + (y - xpt_y)^2) ,i)    
      IF (j EQ 0) THEN grid_RefineContour, x, y, i
    ENDFOR

    PLOT,[x[i]],[y[i]],color=Truecolor('Darkgrey'), PSYM=3,  $ 
         XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

    print,'------->',x[i],y[i]
    print,i
;   Start the contour with the point that's closest to the x-point:
    new_x = x
    new_y = y
    IF (i GT 0) THEN BEGIN
      new_x = [new_x[i:N_ELEMENTS(new_x)-1],new_x[0:i-1]]
      new_y = [new_y[i:N_ELEMENTS(new_y)-1],new_y[0:i-1]]
    ENDIF

    PLOT,new_x,new_y,color=Truecolor('Darkgrey'),  $ 
         XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
    PLOT,[new_x[0]],[new_y[0]],color=Truecolor('Darkgrey'), PSYM=7,  $ 
         XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

    x = [new_x,new_x[0]] ; close the loop
    y = [new_y,new_y[0]]
;
;   --------------------------------------------------------------------
;   SPLICE THE GRID AND ASSIGN THE GRID POINTS TO THE CONTOUR
;
    pos_x = [x[0]]
    pos_y = [y[0]]

    FOR i = 1, N_ELEMENTS(sep_x)-2 DO BEGIN  

      OPLOT,[sep_x[i]],[sep_y[i]],PSYM=6,COLOR=Truecolor('Lightblue')

      ; Find the points that are closest to the boundary points along the 
      ; separatrix:
      search_x = x
      search_y = y
      niterations = 3  ; was 5, 10/01/2012
      FOR j = 0, niterations DO BEGIN  
        proximity = SQRT( ((sep_x[i])[0] - search_x)^2 + ((sep_y[i])[0] - search_y)^2 )
        dummy = MIN(proximity,imin)
        IF (j LT niterations) THEN grid_RefineContour, search_x, search_y, imin
      ENDFOR
      
;      ; Check that the identified grid point is not degenerate (or close to it) with 
;      ; any other grid points:
;      proximity = SQRT( (pos_x - search_x[imin])^2 + (pos_y - search_y[imin])^2 )
;      dummy = MIN(proximity,jmin)
;
;      IF (proximity[jmin] LT 1.0D-4) THEN BEGIN                                         ;  parameter (see a few lines below as well)
;
;        OPLOT,[search_x[imin]],[search_y[imin]],PSYM=6,COLOR=Truecolor('Black')
;        print, 'shit',proximity[jmin]
;
;        ; Needed in some cases:
;        grid_RefineContour, search_x, search_y, imin
;
;        n = N_ELEMENTS(pos_x)
;
;        IF (n EQ 1) THEN BEGIN
;          print, 'do some work -- core ring problems at start'
;          stop
;        ENDIF ELSE BEGIN
;
;          ; Collect the points that are close to the failed point, using the distance
;          ; to the previous grid point as the spatial scale:
;          distance = SQRT( (pos_x[jmin] - pos_x[n-1])^2 +  $ 
;                           (pos_y[jmin] - pos_y[n-1])^2 )
;          proximity = SQRT( (pos_x[jmin] - search_x)^2 + (pos_y[jmin] - search_y)^2 )
;
;          j = WHERE(proximity LT 0.3D*distance AND proximity GT 0.0D, count)
;          IF (count LE 1) THEN BEGIN
;            PRINT,'ERROR grid_AddCoreRings: Points are all too far away'
;            STOP
;          ENDIF
;          search_x = search_x[j]
;          search_y = search_y[j]
;          ; Now select the one that is closest to the previous grid point:
;          proximity = SQRT( (pos_x[n-1] - search_x)^2 + (pos_y[n-1] - search_y)^2 )
;          dummy = MIN(proximity,kmin)
;          print, 'proximity range',proximity
;
;          OPLOT,[search_x[kmin]],[search_y[kmin]],PSYM=6,COLOR=Truecolor('Pink')          
;          print,'proximity',proximity
;
;          proximity = SQRT( (pos_x - search_x[kmin])^2 + (pos_y - search_y[kmin])^2 )
;          dummy = MIN(proximity,lmin)
;          print, 'proximity',proximity[lmin]
;          IF (proximity[lmin] LT 1.0D-5) THEN BEGIN                                         ;  parameter 
;            PRINT,'ERROR grid_AddCoreRings: Correction algorithm failed'
;            STOP
;          ENDIF
;
;          imin = kmin
;        ENDELSE
;      ENDIF

      pos_x = [pos_x,search_x[imin]]
      pos_y = [pos_y,search_y[imin]]

;      PLOT,[search_x[imin]],[search_y[imin]],color=Truecolor('Lightgreen'), PSYM=7,  $ 
;           XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
    ENDFOR

    pos_x = [pos_x,x[0]]
    pos_y = [pos_y,y[0]]
;    pos_x = [x[0],pos_x,x[0]]
;    pos_y = [y[0],pos_y,y[0]]



    ; Check for degenerate points:
    
    n = N_ELEMENTS(pos_x)

    FOR i = 1, n-1 DO BEGIN
      proximity = SQRT( (pos_x[i]-pos_x[i-1])^2 + (pos_y[i]-pos_y[i-1])^2 )

      IF (proximity LT 1.0D-4) THEN BEGIN

        IF (i EQ n-1) THEN j = i - 1 ELSE j = i

        PLOT,[pos_x[j]],[pos_y[j]],color=Truecolor('Lightgreen'), PSYM=6,  $ 
             XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

        pos_x[j] = 0.5D * (pos_x[j-1] + pos_x[j+1])
        pos_y[j] = 0.5D * (pos_y[j-1] + pos_y[j+1])

;        stop
      ENDIF 

    ENDFOR

    last_pos_x = pos_x
    last_pos_y = pos_y


    grid_x = pos_x
    grid_y = pos_y


    PLOT,grid_x,grid_y,color=Truecolor('Lightgreen'), PSYM=7,  $ 
         XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE


;    help,result,/struct
;    print,result.x
;    print,result.y
   
    ; Add the grid point data to the contour data structure:
    ctr = CREATE_STRUCT(ctr,'grid_x',grid_x,'grid_y',grid_y)
    c_array = grid_UpdateStructure(c_array,name,ctr)

  ENDFOR

;stop
  result = c_array

  RETURN, result

END
;
; ======================================================================
;
; ======================================================================
;


