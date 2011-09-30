;
; ======================================================================
;
; ======================================================================
;
FUNCTION grid_CheckIntegrity, c_array, wall, b, debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON params, CORE, SOL, PFZ, FORWARD, BACKWARD, LOWER_NULL, UPPER_NULL, HFS, LFS 
  ; ------------------------------------------------------------------


  ; not sure what to do here yet, other than check integrity of the initial contours and the assembled grid

END
;
; ======================================================================
;
; ======================================================================
;
FUNCTION grid_BuildGrid, c_array, wall, b, debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON params, CORE, SOL, PFZ, FORWARD, BACKWARD, LOWER_NULL, UPPER_NULL, HFS, LFS 
  ; ------------------------------------------------------------------

  tags  = STRUPCASE(TAG_NAMES(c_array))
  nctr  = N_ELEMENTS(tags)

  PRINT,'---------- BUILDING THE FUCKING GRID ------------'


  IF (KEYWORD_SET(debug)) THEN BEGIN
    PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
    OPLOT,wall.x,wall.y,color=Truecolor('Yellow')
    FOR i = 0, nctr-1 DO BEGIN
      ctr = grid_ExtractStructure(c_array,tags[i])      
      OPLOT,ctr.x,ctr.y,color=Truecolor('Red')    
      OPLOT,ctr.grid_x,ctr.grid_y,color=Truecolor('Pink'),PSYM=3  
      npts = N_ELEMENTS(ctr.x)
      XYOUTS,ctr.x[npts/2],ctr.y[npts/2],STRTRIM(STRING(i+1),2),color=Truecolor('White')
    ENDFOR
  ENDIF


;
; ----------------------------------------------------------------------
; COLLECTION SOME INFORMATION ABOUT SEPARATRICES
;
  ; isep[1] = primary separatrix
  ; isep[2] = secondary separatrix, with x-point located on the vessel wall
  ; isep[3] = secondary separatrix, x-point inside the vessel, first leg
  ; isep[4] =  secondary separatrix, x-point inside the vessel, second leg

  isep = MAKE_ARRAY(4,  /INT,VALUE=-1)
  ixpt = MAKE_ARRAY(4,2,/INT,VALUE=-1)

  isearch = [1, 2, 2, 2]
  iresult = [2, 1, 1, 1]

  FOR ictr = 1, nctr DO BEGIN

    ctr = grid_ExtractStructure(c_array,tags[ictr-1])      

    IF (ctr.separatrix EQ 0) THEN CONTINUE

    i = ctr.separatrix - 1

    isep[i] = ictr

    ixpt[i,*] = WHERE(ctr.grid_x EQ b.x[b.null_i[isearch[i]]] AND  $
                      ctr.grid_y EQ b.y[b.null_j[isearch[i]]], count)
    IF (count NE iresult[i]) THEN BEGIN
      PRINT,'ERROR grid_BuildGrid: Mangled primary x-point search'
      PRINT,'  ICTR      =',ictr
      PRINT,'  TAGS      =',tags[ictr-1]
      PRINT,'  SEPATATRIX=',i+1
      PRINT,'  IXPT      =',ixpt[i]
      STOP 
    ENDIF

  ENDFOR

  print,'separatrix data..............'
  print,isep
  print,ixpt[*,0]
  print,ixpt[*,1]

;
; ----------------------------------------------------------------------
; IDENTIFY THE CORE CONTOURS
;
  ; Identify the core contours:
  ilist = -1
  FOR ictr = nctr, 1, -1 DO BEGIN
    ctr = grid_ExtractStructure(c_array,tags[ictr-1])          
    IF (ctr.region EQ CORE) THEN ilist = [ilist,ictr]
  ENDFOR
  ilist = ilist[1:N_ELEMENTS(ilist)-1] ; Get rid of first entry

  print,'core contours', ilist

;
; ----------------------------------------------------------------------
; FIND WHERE TO START THE GRID WHEN FOLLOWING THE WALL
;
  ; Start with the contour that's just outside the low-index primary
  ; strike-point (the inner target strike-point for lower single null):
  ctr  = grid_ExtractStructure(c_array,'contour1')  ; separatrix
  ictr = (ctr.map_out[0])[0]  
  iwall1 = WHERE(wall.ptc EQ ictr AND wall.ptt EQ 1, count)
  IF (count NE 1) THEN BEGIN
    PRINT,'ERROR grid_BuildGrid: Something is wrong with the wall'
    PRINT,'  IWALL1=',iwall1
    STOP 
  ENDIF
  istart = iwall1

  grid_n = 0

  core_active = 1

  i1 = 0
  i2 = 0

  icore = 0

  tube_n  = 0
;
; ----------------------------------------------------------------------
; FOLLOW THE WALL AND BUILD THE GRID RINGS
;
  WHILE (1) DO BEGIN  

    IF (core_active) THEN BEGIN
;  
;     --------------------------------------------------------------------
;     PROCESS THE CORE RINGS FIRST
;  
      IF (icore EQ N_ELEMENTS(ilist)-1) THEN BEGIN
        ictr1 = ilist[icore]
        ictr2 = 1
        i2 = ixpt[0,0]  ; Set the outer starting point to the x-point
      ENDIF ELSE BEGIN
        ictr1 = ilist[icore  ]
        ictr2 = ilist[icore+1]
      ENDELSE
      ctr1 = grid_ExtractStructure(c_array,tags[ictr1-1])      
      ctr2 = grid_ExtractStructure(c_array,tags[ictr2-1])      

      print,'new contour ---core---- ',ictr1,'-----------',FORMAT='(A,I6,A)'
      print,'            ----------- ',ictr2,'-----------',FORMAT='(A,I6,A)'

    ENDIF ELSE BEGIN
;  
;     --------------------------------------------------------------------
;     SCAN BACKWARD ALONG THE WALL TO FIND THE NEIGHBOURING CONTOUR
;  
      print,'new contour ----------- ',iwall1,ictr1,'-----------',FORMAT='(A,2I6,A)'

      iwall2 = iwall1
      WHILE (1) DO BEGIN
        iwall2--
        IF (iwall2 LT 0) THEN iwall2 = N_ELEMENTS(wall.ptc) - 1
        IF (wall.ptc[iwall2] NE -1 AND (wall.ptt[iwall2] EQ 1 OR  $
                                        wall.ptt[iwall2] EQ 3)) THEN BREAK
      ENDWHILE
      ictr2 = (wall.ptc[iwall2])[0]
      ctr2  = grid_ExtractStructure(c_array,tags[ictr2-1])      
   
      print,'            ----------- ',iwall2,ictr2,'-----------',FORMAT='(A,2I6,A)'

      ; Special case when the neighbouring starting location is a tangency point:
      IF (wall.ptt[iwall2] EQ 3) THEN BEGIN
        print,'   scramble!'
        dummy = MIN( SQRT( ((wall.pt1[0,iwall2])[0]-ctr2.grid_x)^2 +  $
                           ((wall.pt1[1,iwall2])[0]-ctr2.grid_y)^2), i2)
        print,wall.pt1[0:1,iwall2]
        print,ctr2.grid_x[i2],ctr2.grid_y[i2],ctr2.grid_s[i2]
      ENDIF

    ENDELSE  ; CORE_ACTIVE block
;  
;   ----------------------------------------------------------------------
;   ASSEMBLE RING
;  
    tube_x1 = -1.0D
    tube_y1 = -1.0D
    tube_x2 = -1.0D
    tube_y2 = -1.0D

    WHILE (1) DO BEGIN
 
      IF (ctr1.separatrix NE 0) THEN BEGIN
        i = ctr1.separatrix - 1

;        print, 'data check',i,ixpt[i,*]

        IF (i EQ 0 AND i1 EQ ixpt[i,0]) THEN BEGIN
          print,'----1 blasterizing a new contour dude----'
          i1 = ixpt[i,1]
        ENDIF

        IF (i EQ 2 AND i1 EQ ixpt[i,0]) THEN BEGIN
          print,'----2 blasterizing a new contour dude----'
          ctr1 = grid_ExtractStructure(c_array,tags[isep[i+1]-1])                
          i1 = ixpt[i+1,0]
        ENDIF

        IF (i EQ 3 AND i1 EQ ixpt[i,0]) THEN BEGIN
          print,'----3 blasterizing a new contour dude----'
          ctr1 = grid_ExtractStructure(c_array,tags[isep[i-1]-1])                
          i1 = ixpt[i-1,0]
        ENDIF

      ENDIF   

;      IF (grid_n EQ 0) THEN BEGIN
        grid_n++
;        print,'i1,i2',i1,i2,N_ELEMENTS(ctr1.grid_x),N_ELEMENTS(ctr2.grid_x)
        IF (core_active) THEN BEGIN
          cell_x = [ctr1.grid_x[i1  ], ctr2.grid_x[i2  ], ctr2.grid_x[i2+1], ctr1.grid_x[i1+1]]
          cell_y = [ctr1.grid_y[i1  ], ctr2.grid_y[i2  ], ctr2.grid_y[i2+1], ctr1.grid_y[i1+1]]
        ENDIF ELSE BEGIN
          cell_x = [ctr2.grid_x[i2  ], ctr1.grid_x[i1  ], ctr1.grid_x[i1+1], ctr2.grid_x[i2+1]]
          cell_y = [ctr2.grid_y[i2  ], ctr1.grid_y[i1  ], ctr1.grid_y[i1+1], ctr2.grid_y[i2+1]]
        ENDELSE
;      ENDIF ELSE BEGIN
;      ENDELSE

      IF (KEYWORD_SET(debug)) THEN  $
        OPLOT,[cell_x,cell_x[0]],  $
              [cell_y,cell_y[0]],color=Truecolor('Lightblue')

      ; Collect the cell data:      
      tube_x1 = [tube_x1,cell_x[0]]
      tube_y1 = [tube_y1,cell_y[0]]
      tube_x2 = [tube_x2,cell_x[1]]
      tube_y2 = [tube_y2,cell_y[1]]

      i1++
      i2++

      IF ((i1 EQ N_ELEMENTS(ctr1.grid_x)-1                       ) OR  $
          (i2 EQ N_ELEMENTS(ctr2.grid_x)-1 AND ctr2.region EQ PFZ)) THEN BEGIN  ; Special exit for broken PFR contours:

        tube_x1 = [tube_x1[1:N_ELEMENTS(tube_x1)-1],cell_x[3]]
        tube_y1 = [tube_y1[1:N_ELEMENTS(tube_y1)-1],cell_y[3]]
        tube_x2 = [tube_x2[1:N_ELEMENTS(tube_x2)-1],cell_x[2]]
        tube_y2 = [tube_y2[1:N_ELEMENTS(tube_y2)-1],cell_y[2]]

        IF (core_active) THEN BEGIN
          ictr0 = ictr1
          ctr0  = ctr1
          ictr1 = ictr2
          ctr1  = ctr2
          ictr2 = ictr0
          ctr2  = ctr0
        ENDIF

        tube = {  $
          region      : ctr2.region           ,  $ 
          ictr1       : ictr2                 ,  $
          ictr2       : ictr1                 ,  $
          psi1        : ctr2.psi              ,  $
          psi2        : ctr1.psi              ,  $
          separatrix1 : ctr2.separatrix       ,  $
          separatrix2 : ctr1.separatrix       ,  $
          map_in      : ctr2.map_in           ,  $
          map_out     : ctr1.map_out          ,  $
          n           : N_ELEMENTS(tube_x1)-1 ,  $
          x1          : tube_x1               ,  $
          y1          : tube_y1               ,  $
          x2          : tube_x2               ,  $
          y2          : tube_y2    }
        tube_n++
        tag = 'tube' + STRING(tube_n,FORMAT='(I0)')
        IF (KEYWORD_SET(tubes)) THEN tubes = CREATE_STRUCT(tubes,tag,tube) ELSE  $
                                     tubes = CREATE_STRUCT(      tag,tube)

        BREAK
      ENDIF

    ENDWHILE
;
;   --------------------------------------------------------------------
;   ADVANCE TO THE NEXT CONTOUR
;
    i1 = 0
    i2 = 0

    IF (core_active) THEN BEGIN
      icore++

;      print, 'checking core completiong',icore,N_ELEMENTS(ilist)

      ; Switch over to following the wall if the core is complete:
      IF (icore EQ N_ELEMENTS(ilist)) THEN BEGIN
        iwall1 = istart
        ictr1  = (wall.ptc[iwall1])[0]
        ctr1   = grid_ExtractStructure(c_array,tags[ictr1-1])      

        core_active = 0
      ENDIF

    ENDIF ELSE BEGIN

      WHILE (1) DO BEGIN
        status = 0
        iwall1++
        IF (iwall1 GT N_ELEMENTS(wall.ptc)-1) THEN iwall1 = 0

;        print,'>>> ------- checking',iwall1
;        print,'>>> pti',wall.pti[iwall1]
;        print,'>>> ptc',wall.ptc[iwall1]
;        print,'>>> ptt',wall.ptt[iwall1]
;help,wall,/struct
;stop

        IF (wall.ptc[iwall1] NE -1 AND wall.ptt[iwall1] EQ 1) THEN BEGIN
          ictr1 = (wall.ptc[iwall1])[0]
          ctr1  = grid_ExtractStructure(c_array,tags[ictr1-1])              
;          print,'ctr region   ',ictr1,ctr1.region
;          print,'ctr in ',ctr1.map_in ,FORMAT='(A,2I6)'
;          print,'ctr out',ctr1.map_out,FORMAT='(A,2I6)'
          CASE ctr1.region OF
            SOL: status = 1
            PFZ: IF (ctr1.map_in[0] NE -1) THEN status = 1 
          ENDCASE
        ENDIF 

;       Special case for tangency points in private flux regions:
        IF (wall.ptc[iwall1] NE -1 AND wall.ptt[iwall1] EQ 3) THEN BEGIN
          ictr1 = (wall.ptc[iwall1])[0]
          ctr1  = grid_ExtractStructure(c_array,tags[ictr1-1])              
;          print,'ctr region ah!',ictr1,ctr1.region
;          print,'ctr in ',ctr1.map_in ,FORMAT='(A,2I6)'
;          print,'ctr out',ctr1.map_out,FORMAT='(A,2I6)'
          CASE ctr1.region OF
            SOL: 
            PFZ: BEGIN
              IF (ctr1.map_in[1] NE -1) THEN BEGIN
                ; Set I1 so that the the ring starts at the tangency point that's
                ; part way down the contour:
                dummy = MIN( SQRT( ((wall.pt1[0,iwall1])[0]-ctr1.grid_x)^2 +  $
                                   ((wall.pt1[1,iwall1])[0]-ctr1.grid_y)^2), i1)
                status = 1
              ENDIF ELSE BEGIN
                IF (ctr1.map_in[1] EQ -1) THEN BEGIN
                  ; Do nothing, because there's no ring outside the section of
                  ; the contour past the tangency point:
                ENDIF ELSE BEGIN
                  print, 'case not handled yet'
                  stop
                ENDELSE
              ENDELSE
             END
          ENDCASE
        ENDIF 
      
        IF (status) THEN BREAK
      ENDWHILE

      ; Exit the loop if back around to the start of the primary separatrix:
      IF (iwall1 EQ istart) THEN BREAK

    ENDELSE

  ENDWHILE

help,tubes,/struct

  result = tubes

  RETURN, result

END
;
; ======================================================================
;
; ======================================================================
;


