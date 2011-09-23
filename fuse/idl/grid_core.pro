;
; ======================================================================
;
; ======================================================================
;
FUNCTION grid_AddCoreRings, c_array, wall, b, debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON params, CORE, SOL, PFZ, FORWARD, BACKWARD, LOWER_NULL, UPPER_NULL, HFS, LFS 
  ; ------------------------------------------------------------------

  tags  = STRUPCASE(TAG_NAMES(c_array))
  nctr  = N_ELEMENTS(tags)

  PRINT,'---------- ADDING CORE RINGS ------------'

  core_n = 3      ; number of core rings to add

  psi   = b.psi ; psi_raw
  psi_x = b.x
  psi_y = b.y

;  help,c_array,/struct   

  psi_val = b.psi_1st_xpoint * 0.99999D

  core_i = 0

  FOR ictr = nctr+1, nctr+core_n DO BEGIN

    cont = 1

    core_step_dist = 0.01   ; spatial separatrion (m) between rings, at the outer midplane
    core_step_psi  = 0.010  ; step in PSI to start with  

    psi_val = psi_val + core_step_psi 
 
    count = 0
;
;   --------------------------------------------------------------------
;   GET CONTOUR FROM PSI MAP AND SEARCH FOR THE SEGMENT OF INTEREST 
;  
    WHILE cont DO BEGIN

      count = count + 1

      ctr = grid_ExtractContour(psi, psi_x, psi_y, psi_val)

      ibrk = WHERE(ctr.dist GT 0.10)  ; PARAMETER
      nseg = N_ELEMENTS(ibrk) 
      ibrk = [-1,ibrk,ctr.n-1]
;
;     ------------------------------------------------------------------
;     LOOP OVER INDIVIDUAL CONTOUR SEGMENTS AND EXAMINE
;  
      FOR iseg = 0, 0 DO BEGIN

        print,'----------------------------------------'
        print,'contour seg:',iseg,ibrk[iseg]+1,ibrk[iseg+1]

        ; Not a real segment so skip it, may need to strengthen this check:
        IF (ibrk[iseg]+2 EQ ibrk[iseg+1] OR  $
            ibrk[iseg]+1 EQ ibrk[iseg+1] OR  $
            ibrk[iseg]   EQ ibrk[iseg+1]) THEN CONTINUE  

        x = ctr.x[ibrk[iseg]+1:ibrk[iseg+1]]
        y = ctr.y[ibrk[iseg]+1:ibrk[iseg+1]]

        PLOT,[x,x[0]],[y,y[0]],color=Truecolor('Orange'),  $
             XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

        print, 'psi_val = ',psi_val
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
          core_i++
          IF (core_i EQ 1     ) THEN map_in  = [-1,-1] ELSE map_in  = [ictr - 1, -1]
          IF (core_i EQ core_n) THEN map_out = [ 1,-1] ELSE map_out = [ictr + 1, -1] 
          ctr = {  $
            state      : 0       ,  $
            origin     : 0       ,  $
            separatrix : 0       ,  $
            region     : CORE    ,  $
            psi        : psi_val ,  $
            map_in     : map_in  ,  $
            map_out    : map_out ,  $
            x          : x       ,  $
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

  i = WHERE(ctr.grid_x EQ ctr.null_x AND ctr.grid_y EQ ctr.null_y, count)
  print,i
  print,ctr.grid_s[i]
  IF (count NE 2) THEN BEGIN
    PRINT,'ERROR grid_AddCoreRings: x-point not properly identified'
    STOP
  ENDIF

  xpt_x  = ctr.null_x
  xpt_y  = ctr.null_y
;
; ----------------------------------------------------------------------
; DETERMINE THE DISTRIBUTION OF POINTS ON THE CORE PORTION OF SEPARATRIX
;
  x = ctr.grid_x[i[0]:i[1]]  ; a bit crude at low resolution
  y = ctr.grid_y[i[0]:i[1]]
  n = N_ELEMENTS(x)

  length = MAKE_ARRAY(n,/DOUBLE,VALUE=0.0D)
  FOR j = 1, N_ELEMENTS(x)-1 DO  $
    length[j] = length[j-1] + SQRT( (x[j] - x[j-1])^2 + (y[j] - y[j-1])^2 )     

  position = (length / length[n-1])[1:n-2]

  print,position
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
    dummy = MIN( ABS(y), i)    
    x = [x[i:N_ELEMENTS(x)-1],x[0:i-1]]
    y = [y[i:N_ELEMENTS(y)-1],y[0:i-1]]

    ; Identify the point on the contour that's closest to the primary
    ; separatrix, then refine the contour there and search again:
    FOR j = 0, 1 DO BEGIN
      dummy = MIN( SQRT( (x - xpt_x)^2 + (y - xpt_y)^2) ,i)    
      IF (j EQ 0) THEN grid_RefineContour, x, y, i
    ENDFOR

    PLOT,[x[i]],[y[i]],color=Truecolor('White'), PSYM=6,  $ 
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

    PLOT,new_x,new_y,color=Truecolor('White'),  $ 
         XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
    PLOT,[new_x[0]],[new_y[0]],color=Truecolor('White'), PSYM=7,  $ 
         XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
;
;   --------------------------------------------------------------------
;   SPLICE THE GRID AND ASSIGN THE GRID POINTS TO THE CONTOUR
;
    x = [new_x,new_x[0]] ; close the loop
    y = [new_y,new_y[0]]
    result = grid_SpliceContour(x,y,position=position)  

    grid_x = result.x
    grid_y = result.y
;    grid_x = result.x[1:N_ELEMENTS(result.x)-2]
;    grid_y = result.y[1:N_ELEMENTS(result.y)-2]

    PLOT,grid_x,grid_y,color=Truecolor('Blue'), PSYM=6,  $ 
         XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

;    help,result,/struct
    print,result.x
    print,result.y
   
    ; Add the grid point data to the contour data structure:
    ctr = CREATE_STRUCT(ctr,'grid_x',grid_x,'grid_y',grid_y)
    c_array = grid_UpdateStructure(c_array,name,ctr)

  ENDFOR

  result = c_array

  RETURN, result

END
;
; ======================================================================
;
; ======================================================================
;


