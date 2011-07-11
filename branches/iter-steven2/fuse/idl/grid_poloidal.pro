;
; ======================================================================
;
FUNCTION grid_GetRegion, c_array, tag

  tags  = STRUPCASE(TAG_NAMES(c_array))

;    ctr = grid_ExtractStructure(c_array,'CONTOUR'+STRTRIM(STRING(ictr),2))      
    print,'tag:'+tag
;    help,c_array,/struct
    j = -1 
    FOR i = 0, N_ELEMENTS(tags)-1 DO BEGIN
;      print,tags[i],tag eq tags[i]
      IF (tag EQ tags[i]) THEN BEGIN
        j = i
        BREAK
      ENDIF
    ENDFOR
    IF (j EQ -1) THEN STOP
;    i = WHERE(TAG_NAMES(c_array) EQ STRUPCASE(tag), count)
;    IF (count EQ 0) THEN STOP
;    ctr = c_array.(i[0]) 
    ctr = c_array.(j) 
;help,ctr,/struct
;stop
  result = ctr.region

  RETURN, result
END
;
; ======================================================================
;
FUNCTION grid_BuildRadialMap, c_array, wall, debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON params, CORE, SOL, PFZ, FORWARD, BACKWARD, LOWER_NULL, UPPER_NULL, HFS, LFS
  ; ------------------------------------------------------------------

  tags  = STRUPCASE(TAG_NAMES(c_array))
  nctr  = N_ELEMENTS(tags)
  nwall = N_ELEMENTS(wall.ptc)

  istart = WHERE(wall.ptc EQ 1 AND wall.ptt EQ 1, count)
  IF (count EQ 0) THEN BEGIN
    PRINT, 'ERROR grid_BuildRadialMap: Starting contour not found'
    STOP
  ENDIF ELSE IF (KEYWORD_SET(debug)) THEN BEGIN
    PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
    OPLOT,wall.x,wall.y,color=Truecolor('Yellow')
    OPLOT,[wall.pt1[0,istart]],[wall.pt1[1,istart]],color=Truecolor('Hotpink'),PSYM=6
    FOR i = 0, nctr-1 DO BEGIN
      ctr = grid_ExtractStructure(c_array,tags[i])      
      OPLOT,ctr.x,ctr.y,color=Truecolor('Red')    
      npts = N_ELEMENTS(ctr.x)
      XYOUTS,ctr.x[npts/2],ctr.y[npts/2],STRTRIM(STRING(i+1),2),color=Truecolor('White')
    ENDFOR
  ENDIF



  ; There should be up to 2 contours mapped to any other contour, which
  ; occurs when there was a tagnecy point:
  map_in  = MAKE_ARRAY(nctr,3,/LONG,VALUE=-1L)
  map_out = map_in

  iwall = istart
  cwall = 1
  ictr_last = 1
  itar_last = 1

  PRINT,'---------- START RADIAL MAP ------------'

  WHILE (cwall LE nwall) DO BEGIN
    cwall++
    iwall++
    IF (iwall GT nwall-1) THEN iwall = 0

    print,nwall,iwall,cwall,FORMAT='(3I8)'    

    ictr = wall.ptc[iwall]  ; contour index of wall point
    itar = wall.ptt[iwall]  ; target 

    ; Cycle to the next wall point if this one isn't associated with a contour:
    IF (ictr EQ -1) THEN CONTINUE

    IF (KEYWORD_SET(debug)) THEN  $
      OPLOT,[wall.pt1[0,iwall]],[wall.pt1[1,iwall]],color=Truecolor('Hotpink'),PSYM=1

    print,'OK---',ictr,itar,ictr_last,itar_last,FORMAT='(A8,4I8)'    

    IF ((itar EQ 1 OR itar EQ 3) AND (itar_last EQ 1 OR itar_last EQ 3)) THEN BEGIN  ; NEW

      tag = 'CONTOUR'+STRTRIM(STRING(ictr),2)
      region = grid_GetRegion(c_array,tag)

      IF (itar NE 3 AND itar_last NE 3) THEN BEGIN
        IF (region EQ SOL OR itar_last EQ 1) THEN i = 0 ELSE i = 1
        IF (map_in[ictr-1,i] NE -1) THEN BEGIN
          PRINT, 'ERROR grid_BuildRadialMap: Inward overload (1)'
          FOR i = 0, nctr-1 DO print,i+1,map_in[i,0:1],format='(3I6)' 
          STOP
        ENDIF
        map_in[ictr-1,i] = ictr_last 
        map_in[ictr-1,2] = region
      ENDIF

      IF ((region EQ SOL AND itar EQ 1 AND itar_last EQ 3) OR  $
          (region EQ PFZ AND itar EQ 3 AND itar_last EQ 1)) THEN BEGIN
        IF (region EQ PFZ) THEN i = 1 ELSE i = 0
        IF (map_in[ictr-1,i] NE -1) THEN BEGIN
          PRINT, 'ERROR grid_BuildRadialMap: Inward overload (2)'
          FOR i = 0, nctr-1 DO print,i+1,map_in[i,0:1],format='(3I6)' 
          STOP
        ENDIF
        map_in[ictr-1,i] = ictr_last 
        map_in[ictr-1,2] = region
      ENDIF

      tag = 'CONTOUR'+STRTRIM(STRING(ictr_last),2)
      region = grid_GetRegion(c_array,tag)

      IF (itar NE 3 AND itar_last NE 3) THEN BEGIN
        IF (region EQ PFZ OR itar_last EQ 1) THEN i = 0 ELSE i = 1
        IF (map_out[ictr_last-1,i] NE -1) THEN BEGIN
          PRINT, 'ERROR grid_BuildRadialMap: Outward overload'
          FOR i = 0, nctr-1 DO print,i+1,map_out[i,0:1],format='(3I6)' 
          STOP
        ENDIF
        map_out[ictr_last-1,i] = ictr 
        map_out[ictr_last-1,2] = region
      ENDIF

      IF ((region EQ SOL AND itar EQ 1 AND itar_last EQ 3) OR  $
          (region EQ PFZ AND itar EQ 3 AND itar_last EQ 1)) THEN BEGIN
        IF (region EQ PFZ) THEN i = 0 ELSE i = 1
        IF (map_out[ictr_last-1,i] NE -1) THEN BEGIN
          PRINT, 'ERROR grid_BuildRadialMap: Outward overload (2)'
          FOR i = 0, nctr-1 DO print,i+1,map_out[i,0:1],format='(3I6)' 
          STOP
        ENDIF
        map_out[ictr_last-1,i] = ictr 
        map_out[ictr_last-1,2] = region
      ENDIF

    ENDIF

    print,'  MAP',map_in [ictr-1,0:1],FORMAT='(A8,3I8)'    
    print,'     ',map_out[ictr-1,0:1],FORMAT='(A8,3I8)'    

    ictr_last = ictr
    itar_last = itar

  ENDWHILE


  ; Small update to map for secondary separatrix:
  FOR ictr = 1, nctr DO BEGIN 
    ctr = grid_ExtractStructure(c_array,tags[ictr-1])      
    IF (ctr.separatrix LT 3) THEN CONTINUE
    ; Scan the wall to find the next point that references a ring and
    ; assign that ring to the connection map:
    iwall = WHERE(wall.ptc EQ ictr AND wall.ptt EQ 2, count)
    IF (count LE 0) THEN BEGIN
      PRINT,'ERROR grid_BuildRadialMap: Unable to find separatrix segment'
      PRINT,'  ICTR=',ictr
      PRINT,'  SEP =',ctr.separatrix 
      STOP      
    ENDIF
    cnt = 0
    WHILE (1 EQ 1) DO BEGIN
      cnt++
      IF (cnt GT nwall) THEN BEGIN
        PRINT,'ERROR grid_BuildRadialMap: Unable to find wall reference'
        STOP      
      ENDIF
      iwall = grid_FindWallNeighbour(wall.pt1,wall.pt2,wall.pti,iwall,2)
      IF (wall.ptc[iwall] GT 0 AND wall.ptt[iwall] EQ 2) THEN BEGIN
        map_in[ictr-1,1] = wall.ptc[iwall]
        IF (ctr.separatrix EQ 4) THEN  map_out[wall.ptc[iwall]-1,1] = ictr
        BREAK
      ENDIF
    ENDWHILE
  ENDFOR


  ; Update the connection map information for each contour:
  FOR i = 0, nctr-1 DO BEGIN 
    ctr = grid_ExtractStructure(c_array,tags[i])      

    dummy = WHERE(STRUPCASE(TAG_NAMES(ctr)) EQ 'MAP_IN',count)
    IF (count EQ 0) THEN BEGIN
      ctr = CREATE_STRUCT(ctr,'map_in',map_in[i,0:1],'map_out',map_out[i,0:1]) 
    ENDIF ELSE BEGIN
      ctr = grid_UpdateStructure(ctr,'map_in' ,map_in [i,0:1])
      ctr = grid_UpdateStructure(ctr,'map_out',map_out[i,0:1])
    ENDELSE
    c_array = grid_UpdateStructure(c_array,tags[i],ctr)
  ENDFOR


  print,'in'  
  FOR i = 0, nctr-1 DO BEGIN 
    ctr = grid_ExtractStructure(c_array,tags[i])      
    print,i+1,map_in[i,0:2],'  ',ctr.map_in[0:1],format='(4I6,A,2I6)'
  ENDFOR

  print,'out'  
  FOR i = 0, nctr-1 DO BEGIN 
    ctr = grid_ExtractStructure(c_array,tags[i])      
    print,i+1,map_out[i,0:2],'  ',ctr.map_out[0:1],format='(4I6,A,2I6)'
  ENDFOR

  result = c_array

  RETURN, result

END
;
; ======================================================================
;
FUNCTION grid_SetPoloidalDomains, c_array, wall, b, debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON params, CORE, SOL, PFZ, FORWARD, BACKWARD, LOWER_NULL, UPPER_NULL, HFS, LFS
  ; ------------------------------------------------------------------

  tags  = STRUPCASE(TAG_NAMES(c_array))
  nctr  = N_ELEMENTS(tags)


  IF (KEYWORD_SET(debug)) THEN BEGIN
    PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
    OPLOT,wall.x,wall.y,color=Truecolor('Yellow')
    FOR i = 0, nctr-1 DO BEGIN
      ctr = grid_ExtractStructure(c_array,tags[i])      
      OPLOT,ctr.x,ctr.y,color=Truecolor('Red')    
      npts = N_ELEMENTS(ctr.x)
      XYOUTS,ctr.x[npts/2],ctr.y[npts/2],STRTRIM(STRING(i+1),2),color=Truecolor('White')
    ENDFOR
  ENDIF


  PRINT,'---------- START DOMAIN IDENTIFICATION ------------'


  ; Process tangency points:

  vector_n = 0  

  section = 0

  FOR ictr = 1, nctr DO BEGIN 
    ctr = grid_ExtractStructure(c_array,tags[ictr-1])      

    IF (ctr.tangent_i EQ -1 AND  $
        (ctr.separatrix LT 3 OR ctr.tangent_p1[0] EQ 0.0)) THEN CONTINUE

    print, ' =================== ictr ================= ',ictr,' ',tags[ictr-1]

    OPLOT,ctr.x,ctr.y,color=Truecolor('Pink')    
    OPLOT,[ctr.tangent_p1[0]],[ctr.tangent_p1[1]],color=Truecolor('Lightgreen'),PSYM=6    

    ; Follow the trajectory into the separatrix:
    p1 =        ctr.tangent_p1
    p2 = 2.0D * ctr.tangent_p1 - ctr.tangent_p2

    OPLOT,[p1[0],p2[0]],[p1[1],p2[1]],color=Truecolor('Lightgreen')    

    section++
    dummy = WHERE(STRUPCASE(TAG_NAMES(ctr)) EQ 'SECTION',count)
    IF (count EQ 0) THEN BEGIN
      ctr = CREATE_STRUCT(ctr,'section',section,'section_x',p1[0],  $
                                                'section_y',p1[1])
    ENDIF ELSE BEGIN
      ctr = grid_UpdateStructure(ctr,'section'  ,[ctr.section  ,section])        
      ctr = grid_UpdateStructure(ctr,'section_x',[ctr.section_x,p1[0]  ])
      ctr = grid_UpdateStructure(ctr,'section_y',[ctr.section_y,p1[1]  ])      
    ENDELSE
    c_array = grid_UpdateStructure(c_array,tags[ictr-1],ctr)

    ; Special case where processing the secondary separatrix:
    IF (ctr.separatrix GE 3) THEN BEGIN
      FOR ictr2 = 1, nctr DO BEGIN
        IF (ictr2 EQ ictr) THEN CONTINUE
        ctr2 = grid_ExtractStructure(c_array,tags[ictr2-1])      
        IF (ctr2.separatrix GE 3) THEN BEGIN
          dummy = WHERE(STRUPCASE(TAG_NAMES(ctr2)) EQ 'SECTION', count)
          IF (count EQ 0) THEN BEGIN
            ctr2 = CREATE_STRUCT(ctr2,'section',section,'section_x',p1[0],  $
                                                        'section_y',p2[1])
          ENDIF ELSE BEGIN
            ctr2 = grid_UpdateStructure(ctr2,'section'  ,[ctr2.section  ,section])        
            ctr2 = grid_UpdateStructure(ctr2,'section_x',[ctr2.section_x,p1[0]  ])
            ctr2 = grid_UpdateStructure(ctr2,'section_y',[ctr2.section_y,p1[1]  ])
          ENDELSE     
          c_array = grid_UpdateStructure(c_array,tags[ictr2-1],ctr2)
          BREAK
        ENDIF
      ENDFOR
    ENDIF

    cnt=0
    ictr3 = ictr
   
    ; For each contour with a wall tangency point, scan inward through the
    ; other contours looking for a corresponding intersection with the inner-
    ; facing normal of the tangency point:

    WHILE (1 EQ 1) DO BEGIN
      
      print, 'round =========>',cnt,ictr3

      ; Not trying to use the connection map here -- why not do you think?
      min_ictr2 = -1
      min_ctr2  = -1
      min_inter = -1
      min_s     = 1.0D+06
      FOR ictr2 = 1, nctr DO BEGIN
        IF (ictr2 EQ ictr3) THEN CONTINUE
        ctr2 = grid_ExtractStructure(c_array,tags[ictr2-1])      
        ; Special case where processing the secondary separatrix
        ; intersections:
        IF (ctr.separatrix GE 3 AND ctr2.separatrix GE 3) THEN CONTINUE  
        n = N_ELEMENTS(ctr2.x)
        seg1 = MAKE_ARRAY(2,n-1,/DOUBLE,VALUE=0.0D)
        seg2 = seg1
        seg1[0,*] = ctr2.x[0:n-2]
        seg1[1,*] = ctr2.y[0:n-2]
        seg2[0,*] = ctr2.x[1:n-1]
        seg2[1,*] = ctr2.y[1:n-1]
        inter = grid_Intersection(p1, p2, seg1, seg2, 0, status=status)
        IF (status) THEN BEGIN
          print, 'intersection',ictr2,N_ELEMENTS(inter.s)
          IF (MIN(inter.s) LT min_s) THEN BEGIN
            min_inter = inter
            min_s     = MIN(inter.s)
            min_ictr2 = ictr2
            min_ctr2  = ctr2
          ENDIF
        ENDIF
      ENDFOR
      print,'taken -->',min_ictr2
      ictr2 = min_ictr2
      ctr2  = min_ctr2
      inter = min_inter
      dummy = MIN(inter.s,min_i)
      inter_i = inter.i[min_i]
      inter_x = inter.x[min_i]
      inter_y = inter.y[min_i]

      ; Check if there's an intersection with a previous cross-field vector:
      IF (vector_n EQ 0) THEN BEGIN
        vector_seg1 = MAKE_ARRAY(2,1000,/DOUBLE,VALUE=0.0D)
        vector_seg2 = MAKE_ARRAY(2,1000,/DOUBLE,VALUE=0.0D)
        vector_seg1[0,0] = p1[0] 
        vector_seg1[1,0] = p1[1]
        vector_seg2[0,0] = inter_x
        vector_seg2[1,0] = inter_y
        vector_c  = [ictr]
        vector_n  = 1
      ENDIF ELSE BEGIN
        FOR ivec = 0, vector_n-1 DO BEGIN
          IF (ictr EQ vector_c[ivec]) THEN CONTINUE 
          vector_inter = grid_Intersection(p1, [inter_x,inter_y], vector_seg1[*,ivec],  $
                                                                  vector_seg2[*,ivec], 0, status=status)          
          IF (status) THEN BREAK
        ENDFOR
        IF (ivec NE vector_n) THEN BEGIN
          print, 'collision detected!'
          OPLOT,[p1[0],p2[0]],[p1[1],p2[1]],COLOR=Truecolor('Gold')
          OPLOT,[vector_seg1[0,ivec:ivec+1], vector_seg2[0,ivec:ivec+1]],  $
                [vector_seg1[1,ivec:ivec+1], vector_seg2[1,ivec:ivec+1]],  $
                COLOR=Truecolor('Orange')

          v1 = p2                         - p1 
          v2 = vector_seg2[*,ivec] - vector_seg1[*,ivec]
          len1 = SQRT( v1[0]^2 + v1[1]^2 )
          len2 = SQRT( v2[0]^2 + v2[1]^2 )

          v2 = (len1 / len2) * v2

          frac = (vector_inter.t[0]) * 0.7D   ; parameter

          v3 = frac * v1 + (1.0D - frac) * v2

          p2 = p1 + v3

print,vector_inter.s
print,vector_inter.t
help,frac
          OPLOT,[p1[0],p2[0]],[p1[1],p2[1]],COLOR=Truecolor('Gold')
          CONTINUE
        ENDIF ELSE BEGIN
          vector_n++
          vector_seg1[0,vector_n-1] = p1[0] 
          vector_seg1[1,vector_n-1] = p1[1]
          vector_seg2[0,vector_n-1] = inter_x
          vector_seg2[1,vector_n-1] = inter_y
          vector_c                  = [vector_c,ictr]
        ENDELSE
      ENDELSE

      print, 'separatrix ====>',ictr,ctr2.separatrix

      OPLOT,ctr2.x,ctr2.y,color=Truecolor('White')    

      OPLOT,[inter_x],[inter_y],color=Truecolor('White'),PSYM=6
      OPLOT,[ctr2.x[inter_i  ]],[ctr2.y[inter_i  ]],color=Truecolor('Red'),PSYM=6
      OPLOT,[ctr2.x[inter_i+1]],[ctr2.y[inter_i+1]],color=Truecolor('Red'),PSYM=6

      ; Add point to the contour and register the designation of the 
      ; tangency point:

      ;   Note that this intersection point is the results of a linear
      ;   interpolation, so the curvature of the contour will be a bit
      ;   off.  Might be worth refining the contour after the intersection 
      ;   is identified and then finding a new intersection, in order
      ;   to minimize this effect:

      n = N_ELEMENTS(ctr2.x)
      x = [ctr2.x[0:inter_i],inter_x,ctr2.x[inter_i+1:n-1]]
      y = [ctr2.y[0:inter_i],inter_y,ctr2.y[inter_i+1:n-1]]
      ctr2 = grid_UpdateStructure(ctr2,'x',x)
      ctr2 = grid_UpdateStructure(ctr2,'y',y)
      dummy = WHERE(STRUPCASE(TAG_NAMES(ctr2)) EQ 'SECTION',count)
      IF (count EQ 0) THEN BEGIN
        ctr2 = CREATE_STRUCT(ctr2,'section',section,'section_x',inter_x,  $
                                                    'section_y',inter_y)
      ENDIF ELSE BEGIN
        ctr2 = grid_UpdateStructure(ctr2,'section'  ,[ctr2.section  ,section])        
        ctr2 = grid_UpdateStructure(ctr2,'section_x',[ctr2.section_x,inter_x])
        ctr2 = grid_UpdateStructure(ctr2,'section_y',[ctr2.section_y,inter_y])      
      ENDELSE
      c_array = grid_UpdateStructure(c_array,tags[ictr2-1],ctr2)

      ; Quit when the separatrix has been reached:
      IF (ictr2 EQ 1) THEN BREAK

      ; Calculate the new cross-field trajectory:
      p1_new = [inter_x,inter_y]
      IF (ctr.region EQ PFZ) THEN direction = 2 ELSE direction = 1
      p2_new = grid_GetOrthogonal(ctr2.x,ctr2.y,p1_new,direction,1.5D,ictr=ictr2,ctrs=c_array)

      OPLOT,[p1_new[0],p2_new[0]],[p1_new[1],p2_new[1]],color=Truecolor('Blue')

      cnt++

      p1 = p1_new
      p2 = p2_new              
 
      ictr3 = ictr2
    ENDWHILE

  ENDFOR

  ; Need to add two more slice points at the primary x-point:
  ictr = grid_GetIndex(c_array,'separatrix',1)
  ctr  = grid_ExtractStructure(c_array,tags[ictr-1])        
  focus_x = ctr.null_x
  focus_y = ctr.null_y
  i = WHERE(ctr.x EQ focus_x AND ctr.y EQ focus_y)
  section++
  dummy = WHERE(STRUPCASE(TAG_NAMES(ctr2)) EQ 'SECTION',count)
  IF (count EQ 0) THEN BEGIN
    ctr = CREATE_STRUCT(ctr,'section'  ,[section,section+1],  $
                            'section_x',[focus_x,focus_x  ],  $
                            'section_y',[focus_y,focus_y  ])
  ENDIF ELSE BEGIN
    ctr = grid_UpdateStructure(ctr,'section'  ,[ctr.section  ,section,section+1])        
    ctr = grid_UpdateStructure(ctr,'section_x',[ctr.section_x,focus_x,focus_x  ])
    ctr = grid_UpdateStructure(ctr,'section_y',[ctr.section_y,focus_y,focus_y  ])      
  ENDELSE
  section++
  c_array = grid_UpdateStructure(c_array,tags[ictr-1],ctr)

  ; Sort the sections registered, starting from the inner target:
  special = -1
  FOR ictr = 1, nctr DO BEGIN 
    ctr = grid_ExtractStructure(c_array,tags[ictr-1])      

    dummy = WHERE(STRUPCASE(TAG_NAMES(ctr)) EQ 'SECTION', count)
    IF (count EQ 0) THEN CONTINUE

    n = N_ELEMENTS(ctr.section)
    length = MAKE_ARRAY(n,/FLOAT,VALUE=-1.0)
    FOR i = 0, n-1 DO BEGIN 
      j = WHERE(ctr.x EQ ctr.section_x[i] AND ctr.y EQ ctr.section_y[i], count)
      IF (count NE 1 AND NOT (count EQ 2 AND ctr.separatrix EQ 1)) THEN BEGIN
        PRINT,'ERROR grid_CreatePoloidalPoints: Problem finding section marker(s) on contour'
        PRINT,'  I       =',i
        PRINT,'  SECTION =',section
        STOP
      ENDIF
      ; Special case for the primary separatrix:
      IF (count EQ 2) THEN BEGIN
        special++
        length[i] = grid_Length(ctr.x[0:j[0]],ctr.y[0:j[0]])  
        IF (special EQ 1) THEN BEGIN
          ; Need to add the distance around the core for the second 
          ; section marker:  
          j = WHERE(ctr.x EQ ctr.null_x AND ctr.y EQ ctr.null_y, count)
          IF (count NE 2) THEN BEGIN
            PRINT,'ERROR grid_CreatePoloidalPoints: Malformed separatrix'
            STOP
          ENDIF
          length[i] = length[i] + grid_Length(ctr.x[j[0]:j[1]],ctr.y[j[0]:j[1]])            
        ENDIF
      ENDIF ELSE  $
        length[i] = grid_Length(ctr.x[0:j],ctr.y[0:j]) 
    ENDFOR
    i = SORT(length)
    ctr = grid_UpdateStructure(ctr,'section'  ,[-999 ,ctr.section  [i],-999 ])        
    ctr = grid_UpdateStructure(ctr,'section_x',[-1.0D,ctr.section_x[i],-1.0D])
    ctr = grid_UpdateStructure(ctr,'section_y',[-1.0D,ctr.section_y[i],-1.0D])
    c_array = grid_UpdateStructure(c_array,tags[ictr-1],ctr)
  ENDFOR

  ; Set end points codes for separatricies:
  ctr = grid_ExtractStructure(c_array,tags[0])      
  ctr.section[0                        ] = -1
  ctr.section[N_ELEMENTS(ctr.section)-1] = -2
  c_array = grid_UpdateStructure(c_array,tags[0],ctr)
  FOR iseparatrix = 3, 4 DO BEGIN
    ictr = grid_GetIndex(c_array,'separatrix',iseparatrix)
    IF (ictr EQ -1) THEN BREAK
    ctr = grid_ExtractStructure(c_array,tags[ictr-1])      
    ctr.section[0                        ] = -2*(iseparatrix-3)-3  ; so....
    ctr.section[N_ELEMENTS(ctr.section)-1] = -2*(iseparatrix-3)-4
    c_array = grid_UpdateStructure(c_array,tags[ictr-1],ctr)
  ENDFOR
;  FOR ictr = 1, nctr DO BEGIN
;    ctr = grid_ExtractStructure(c_array,tags[ictr-1])      
;    IF (ctr.separatrix LT 3) THEN CONTINUE
;    val = ctr.separatrix
;    ctr.section[0                        ] = -2*(val-3)-3
;    ctr.section[N_ELEMENTS(ctr.section)-1] = -2*(val-3)-4
;    c_array = grid_UpdateStructure(c_array,tags[ictr-1],ctr)
;  ENDFOR

  ; Setup section structure elements for all rings that don't have them yet:
  FOR ictr = 1, nctr DO BEGIN
    ctr = grid_ExtractStructure(c_array,tags[ictr-1])
    dummy = WHERE(STRUPCASE(TAG_NAMES(ctr)) EQ 'SECTION', count)
    IF (count GT 0) THEN CONTINUE
    ctr = CREATE_STRUCT(ctr,'section'  ,[-999 ,-999 ],  $
                            'section_x',[-1.0D,-1.0D],  $
                            'section_y',[-1.0D,-1.0D])
    c_array = grid_UpdateStructure(c_array,tags[ictr-1],ctr)
  ENDFOR


  WHILE (1 EQ 1) DO BEGIN
    print, ' '
    print, ' ============= pass ============='
    print, ' '

    status = 0

    FOR ictr = 2, nctr DO BEGIN
      ctr = grid_ExtractStructure(c_array,tags[ictr-1])      

      dummy = WHERE(STRUPCASE(TAG_NAMES(ctr)) EQ 'SECTION', count)  ; *** TEMP ***
      IF (count EQ 0) THEN CONTINUE
     
      IF (ctr.section[0] NE -999) THEN CONTINUE

      status = 1

      print, 'processing end points ---- ',ictr,' -------'
    
      IF (ctr.region EQ PFZ) THEN ictr2 = ctr.map_out[0] ELSE  $
                                  ictr2 = ctr.map_in [0]

      print,'  region,map=',ictr ,ctr.region,ctr.map_in,ctr.map_out,format='(a,7I6)'

      ctr2  = grid_ExtractStructure(c_array,tags[ictr2-1])

      print,'  region,map=',ictr2,ctr.region,ctr2.map_in,ctr2.map_out,format='(a,7I6)'

      IF (ctr2.section[0] EQ -999) THEN CONTINUE

      PRINT, '     >>>>> go >>>> '

      status2 = 0

      ; Primary separatrix is a special case:
      IF (ctr2.separatrix EQ 1) THEN BEGIN
        ctr.section[0                        ] = ctr2.section[0                         ] 
        ctr.section[N_ELEMENTS(ctr.section)-1] = ctr2.section[N_ELEMENTS(ctr2.section)-1]
        status2 = 1
      ENDIF

      ; The reference ring doesn't have a tangency point so just copy the end points:
      IF (NOT status2 AND ctr2.tangent_i EQ -1) THEN BEGIN
        ctr.section[0                        ] = ctr2.section[0                         ] 
        ctr.section[N_ELEMENTS(ctr.section)-1] = ctr2.section[N_ELEMENTS(ctr2.section)-1]
        status2 = 1
      ENDIF

      ; The reference ring does have a tangency point, so need to figure out how
      ; the end points of the focus ring map to this point:
      IF (NOT status2 AND ctr2.tangent_i NE -1) THEN BEGIN
        location = -1
        IF (ctr.region EQ PFZ) THEN  $
          map = ctr2.map_in  ELSE  $
          map = ctr2.map_out
        IF (map[0] EQ ictr) THEN location = 1
        IF (map[1] EQ ictr) THEN location = 2
        IF (location EQ -1) THEN BEGIN
          PRINT,'ERROR grid_SetPoloidalDomains: Something wrong with radial mapping'
          PRINT,'  ICTR =',ictr
          PRINT,'  ICTR2=',ictr2
          STOP
        ENDIF
        print,'  location  =',location
        ; Need to identify the section marker associated with the tangency point:
        i = WHERE(ctr2.section_x EQ ctr2.tangent_p1[0] AND  $
                  ctr2.section_y EQ ctr2.tangent_p1[1], count)
        IF (count EQ 0) THEN BEGIN
          PRINT,'ERROR grid_SetPoloidalDomains: Unable to find tangency point'
          PRINT,'  ICTR2=',ictr2
          STOP
        ENDIF
        CASE location OF
          1: BEGIN
            ctr.section[0                        ] = ctr2.section[0] 
            ctr.section[N_ELEMENTS(ctr.section)-1] = ctr2.section[i] 
            END
          2: BEGIN
            ctr.section[0                        ] = ctr2.section[i                         ] 
            ctr.section[N_ELEMENTS(ctr.section)-1] = ctr2.section[N_ELEMENTS(ctr2.section)-1] 
            END
        ENDCASE
        status2 = 1
      ENDIF

      IF (status2 EQ 0) THEN BEGIN
        PRINT,'ERROR grid_SetPoloidalDomains: Unrecognized section assignment case'
        PRINT,'  ICTR =',ictr
        PRINT,'  ICTR2=',ictr2
        STOP
      ENDIF

      print,'  section=',ctr.section

      c_array = grid_UpdateStructure(c_array,tags[ictr-1],ctr)
    ENDFOR

    IF (NOT status) THEN BREAK

  ENDWHILE


  FOR ictr = 1, nctr DO BEGIN 
    ctr = grid_ExtractStructure(c_array,tags[ictr-1])      

    print, 'sections',ictr,'  ',ctr.section
  ENDFOR


  ; Update for secondary separatrix, since the end points are not set 
  ; properly at this point:

  ; Correct the secondary PFZ:
  ictr1 = grid_GetIndex(c_array,'separatrix',3)
  IF (ictr1 NE -1) THEN BEGIN
    ictr2 = grid_GetIndex(c_array,'separatrix',4)
    ctr1  = grid_ExtractStructure(c_array,tags[ictr1-1])     
    ctr2  = grid_ExtractStructure(c_array,tags[ictr2-1])     
    FOR ictr = 1, nctr DO BEGIN
      ctr = grid_ExtractStructure(c_array,tags[ictr-1])     
      IF (ctr.region NE PFZ) THEN CONTINUE
      IF ( ctr.section[N_ELEMENTS( ctr.section)-1] EQ  $
          ctr2.section[N_ELEMENTS(ctr2.section)-1]) THEN BEGIN
        ctr.section[N_ELEMENTS(ctr.section)-1] = ctr1.section[N_ELEMENTS(ctr1.section)-1]
        c_array = grid_UpdateStructure(c_array,tags[ictr-1],ctr)              
      ENDIF
    ENDFOR
  ENDIF

  print,' '
  FOR ictr = 1, nctr DO BEGIN 
    ctr = grid_ExtractStructure(c_array,tags[ictr-1])      
    print, 'sections',ictr,'  ',ctr.section
  ENDFOR

  ; Need to adjust the end values in the SOL when there's a secondary
  ; x-point
  FOR ictr = 1, nctr DO BEGIN 
    ctr = grid_ExtractStructure(c_array,tags[ictr-1])      

    IF (ctr.separatrix LT 3) THEN CONTINUE

    IF (ctr.separatrix EQ 3) THEN end_point = 0 ELSE end_point = 1
 
    ictr2 = ctr.map_in[end_point] 
    ctr2  = grid_ExtractStructure(c_array,tags[ictr2-1])      

    location = -1
    IF (ctr2.map_out[0] EQ ictr) THEN location = 1
    IF (ctr2.map_out[1] EQ ictr) THEN location = 2
    IF (location EQ -1) THEN BEGIN
      PRINT,'ERROR grid_SetPoloidalDomains: Something wrong with radial mapping (2)'
      PRINT,'  ICTR =',ictr
      PRINT,'  ICTR2=',ictr2
      STOP
    ENDIF
    IF ((end_point EQ 0 AND location EQ 2) OR  $
        (end_point EQ 1 AND location EQ 1)) THEN BEGIN
      ; Need to identify the section marker associated with the tangency point:
      i = WHERE(ctr2.section_x EQ ctr2.tangent_p1[0] AND  $
                ctr2.section_y EQ ctr2.tangent_p1[1], count)
      IF (count EQ 0) THEN BEGIN
        PRINT,'ERROR grid_SetPoloidalDomains: Unable to find tangency point'
        PRINT,'  ICTR2=',ictr2
        STOP
      ENDIF
    ENDIF ELSE  $
      IF (end_point EQ 0) THEN i = 0 ELSE i = N_ELEMENTS(ctr2.section)-1

    IF (end_point EQ 0) THEN old_value = ctr.section[0                        ] ELSE  $
                             old_value = ctr.section[N_ELEMENTS(ctr.section)-1]
    new_value = ctr2.section[i]

    ; Now need to scan all other rings and update any references to the old
    ; target index:
    FOR ictr2 = 1, nctr DO BEGIN 
      ctr2 = grid_ExtractStructure(c_array,tags[ictr2-1])      
      i = WHERE(ctr2.section EQ old_value, count)
      IF (count EQ 0) THEN CONTINUE
      ctr2.section[i] = new_value
      c_array = grid_UpdateStructure(c_array,tags[ictr2-1],ctr2)
    ENDFOR

  ENDFOR

  print,' '
  FOR ictr = 1, nctr DO BEGIN 
    ctr = grid_ExtractStructure(c_array,tags[ictr-1])      
    print, 'sections',ictr,'  ',ctr.section
  ENDFOR


; Check that the ordering of section points on each ring is the same, i.e. that
; the radial intersections computed for one point don't intersect any others:

; I don't know if there's an easy way to do this actually... perhaps better to just make the intersection
; finding code more robust to start with...

;  section_n    = 0
;  FOR ictr = 1, nctr DO BEGIN 
;    ctr = grid_ExtractStructure(c_array,tags[ictr-1])      
;    IF (section_n EQ 0) THEN BEGIN
;      section_list = ctr.section
;      section_n    = N_ELEMENTS(section_list)
;    ENDIF ELSE BEGIN
;      FOR i = 0, N_ELEMENTS(ctr.section)-1 DO BEGIN
;        j = WHERE(ctr.section[i] EQ section_list, count)
;        IF (count EQ 0) THEN BEGIn
;          section_n 
;          section_list 
;        ENDIF ELSE BEGIN
;
;        ENDIF
;      ENDFOR
;    ENDELSE 
;  ENDFOR
  

  result = c_array

  RETURN, result
END
;
; ======================================================================
;
FUNCTION grid_SpliceContour, val_x, val_y,  $
                             mode     = mode    ,  $
                             spacing  = spacing ,  $
                             position = position

;  IF (1 EQ 1) THEN BEGIN
  IF (KEYWORD_SET(position)) THEN BEGIN

    print,'position',position

  ENDIF ELSE BEGIN

    IF (NOT KEYWORD_SET(mode)) THEN BEGIN 
      PRINT,'ERROR grid_SpliceContour: MODE not set'
      STOP
    ENDIF

;   Setup the distribution of points along the contour:
    CASE mode OF
      1: BEGIN

        IF (NOT KEYWORD_SET(spacing)) THEN BEGIN 
          PRINT,'ERROR grid_SpliceContour: SPACING not set'
          STOP
        ENDIF

        frac = grid_GetFrac(val_x,val_y,dist=dist)
        n    = MAX([3,LONG(dist[N_ELEMENTS(dist)-1] / spacing)+ 1])  ; parameter
        step = 1.0D / (DOUBLE(n-1))
;        print,'n',n
;        print,'step',step,dist[N_ELEMENTS(dist)-1]
        position = [step]
        FOR i = 1, n-3 DO position = [position,position[i-1]+step]
;        print,'poition',position
        END
      ELSE: BEGIN
        PRINT,'ERROR grid_SpliceContour: Unrecognised mode'
        PRINT,'  MODE=',mode
        STOP
        END
    ENDCASE

  ENDELSE

  oplot,val_x,val_y,color=Truecolor('Lightblue')

;print,'  position',position

  ; Start with end point:
  result_x = val_x[0]
  result_y = val_y[0]

  FOR i = 0, N_ELEMENTS(position)-1 DO BEGIN
    x = val_x
    y = val_y
    ; Pass through twice to find where to extract the point, refining
    ; the contour on the first iteration:
    FOR j = 0, 1 DO BEGIN
      frac = grid_GetFrac(x,y,dist=dist)
      FOR k = 0, N_ELEMENTS(frac)-2 DO  $
        IF (frac[k] LT position[i] AND frac[k+1] GT position[i]) THEN BREAK
      IF (k EQ N_ELEMENTS(frac)-1) THEN BEGIN
        PRINT,'ERROR grid_SpliceContour: Unable to location position'
        STOP
      ENDIF
      IF (j EQ 0) THEN grid_RefineContour, x, y, k
    ENDFOR

    ; Linearly interpolate to extract the point:  (*** should really use a SPLINE ***)
    frac = (position[i] - frac[k]) / (frac[k+1] - frac[k])    

;    print,frac
    new_x = x[k] + frac * (x[k+1] - x[k])
    new_y = y[k] + frac * (y[k+1] - y[k])

    result_x = [result_x,new_x]
    result_y = [result_y,new_y]
  ENDFOR

  ; And, finally, add the end point:
  result_x = [result_x,val_x[N_ELEMENTS(val_x)-1]]
  result_y = [result_y,val_y[N_ELEMENTS(val_y)-1]]

;  oplot,result_x,result_y,color=Truecolor('Lightblue'),PSYM=6

  result = { x      : result_x                ,  $
             y      : result_y                ,  $
             t      : [0.0, position, 1.0]    ,  $
             length : dist[N_ELEMENTS(dist)-1]}

  RETURN, result

END
;
; ======================================================================
;
FUNCTION grid_CreatePoloidalPoints, c_array, wall, debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON params, CORE, SOL, PFZ, FORWARD, BACKWARD, LOWER_NULL, UPPER_NULL, HFS, LFS
  ; ------------------------------------------------------------------

  tags  = STRUPCASE(TAG_NAMES(c_array))
  nctr  = N_ELEMENTS(tags)

  IF (KEYWORD_SET(debug)) THEN BEGIN
    PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
    OPLOT,wall.x,wall.y,color=Truecolor('Yellow')
;    OPLOT,wall.pt1[0,*],wall.pt1[1,*],color=Truecolor('Yellow'),PSYM=7
    FOR i = 0, nctr-1 DO BEGIN
      ctr = grid_ExtractStructure(c_array,tags[i])      
      OPLOT,ctr.x,ctr.y,color=Truecolor('Red')    
      npts = N_ELEMENTS(ctr.x)
      XYOUTS,ctr.x[npts/2],ctr.y[npts/2],STRTRIM(STRING(i+1),2),color=Truecolor('White')
    ENDFOR
  ENDIF

  PRINT,'---------- SETTING UP BASE POLOIDAL POINTS ------------'

  separatrix_flag = [0,0]

  nsection = 0

        stopper = 0

  FOR ipass = 1, 2 DO BEGIN

    FOR ictr = 1, nctr DO BEGIN 
      ctr = grid_ExtractStructure(c_array,tags[ictr-1])      

      IF ((ipass EQ 1 AND ctr.separatrix EQ 0) OR  $
          (ipass EQ 2 AND ctr.separatrix NE 0)) THEN CONTINUE

      print, '========================================== '
      print, 'ictr ----- ',ictr,' ',tags[ictr-1]
      print, '========================================== '
   

;      print,'map',ctr.map_in , FORMAT='(A,2I6)'
;      print,'   ',ctr.map_out, FORMAT='(A,2I6)'
   
      dummy = WHERE(STRUPCASE(TAG_NAMES(ctr)) EQ 'SECTION',count)
      IF (count EQ 0) THEN BEGIN
        CONTINUE
      ENDIF
   
      print,'section',ctr.section
   
      x = ctr.x
      y = ctr.y
   
      ; Add in the end points from the wall array:
      FOR j = 1, 2 DO BEGIN
        i = WHERE(wall.ptc EQ ictr AND wall.ptt EQ j)
        p1 = wall.pt1[0:1,i]
        i = grid_LocatePoint(x,y,p1)
        x = [x[0:i],p1[0],x[i+1:N_ELEMENTS(x)-1]]
        y = [y[0:i],p1[1],y[i+1:N_ELEMENTS(y)-1]]
        IF (j EQ 1) THEN i = 0 ELSE i = N_ELEMENTS(ctr.section)-1
        ctr.section_x[i] = p1[0]
        ctr.section_y[i] = p1[1]
      ENDFOR
   
      OPLOT,x,y,color=Truecolor('Pink')    
   
;      IF (ctr.tangent_i NE -1) THEN  $
;        OPLOT,[ctr.tangent_p1[0]],[ctr.tangent_p1[1]],color=Truecolor('Lightgreen'),PSYM=6    
   
      n = N_ELEMENTS(ctr.section)
   
      FOR isec = 0, n-2 DO BEGIN
   
        i1 = ctr.section[isec  ]  
        i2 = ctr.section[isec+1]
   
        x1 = ctr.section_x[isec]
        y1 = ctr.section_y[isec]
        j = WHERE(x EQ x1 AND y EQ y1, count)
        IF (count GT 1) THEN BEGIN
          IF (ctr.separatrix EQ 1) THEN BEGIN
            IF (separatrix_flag[0] EQ 0) THEN BEGIN
              j = j[0] & separatrix_flag[0] = 1
            ENDIF ELSE  $
              j = j[1]
          ENDIF ELSE BEGIN
            PRINT,'ERROR grid_CreatePoloidalPoints: Multiple points identified (1)'
            STOP
          ENDELSE
        ENDIF ELSE BEGIN
          IF (count EQ 0) THEN BEGIN
            PRINT,'ERROR grid_CreatePoloidalPoints: Point not found (1)'
            PRINT,'x',x
            PRINT,'y',y
            PRINT,'x1,y1',x1,y1
            STOP
          ENDIF
        ENDELSE
   
        x2 = ctr.section_x[isec+1]
        y2 = ctr.section_y[isec+1]
        k = WHERE(x EQ x2 AND y EQ y2, count)
        IF (count GT 1) THEN BEGIN
          IF (ctr.separatrix EQ 1) THEN BEGIN
            IF (separatrix_flag[1] EQ 0) THEN BEGIN
              k = k[0] & separatrix_flag[1] = 1
            ENDIF ELSE  $
              k = k[1]
          ENDIF ELSE BEGIN
            PRINT,'ERROR grid_CreatePoloidalPoints: Multiple points identified (2)'
            STOP
          ENDELSE
        ENDIF ELSE BEGIN
          IF (count EQ 0) THEN BEGIN
            PRINT,'ERROR grid_CreatePoloidalPoints: Point not found (2)'
            STOP
          ENDIF
        ENDELSE
       
        print,'isec,j,k',isec,j,k,i1,i2,N_ELEMENTS(ctr.x),x1,y1,x2,y2,FORMAT='(A,3I6,6X,3I6,6X,4F10.4)'
   
        IF (k EQ j+1) THEN BEGIN
          ; The identified points are consecutive along the ring, so 
          ; need to interpolate before trying to splice:
          frac = DINDGEN(11) / 10.0D
          new_x = x[j]
          new_y = y[j]
          FOR l = 1, N_ELEMENTS(frac)-2 DO BEGIN
            new_x = [new_x, x[j] + frac[l] * (x[k] - x[j])]  ; Shorthand way of coding this wouldn't work...
            new_y = [new_y, y[j] + frac[l] * (y[k] - y[j])]
          ENDFOR
          new_x = [new_x,x[k]]
          new_y = [new_y,y[k]]
        ENDIF ELSE BEGIN
          new_x = x[j:k]
          new_y = y[j:k]
        ENDELSE
   
   
;          help,section_t,/struct
   
        ; Scan to look for the current interval in the list:
        IF (ictr EQ 1) THEN BEGIN
          count1 = 0 & i3 = -1 
          count2 = 0 & i4 = -1 
        ENDIF ELSE BEGIN
          i3 = WHERE(section_i1 EQ i1, count1)
          i4 = WHERE(section_i2 EQ i2, count2)

          print,'section_i1',section_i1
          print,'section_i2',section_i2
        ENDELSE

        print,'i1,3',i1,i3
        print,'i2,4',i2,i4
   


        IF (count1 EQ 0 OR count2 EQ 0) THEN BEGIN

          print,'         ============ path 1 ============'

          ; Process the separatrix, which provides the point distrubution reference to the
          ; rest of the grid (true?):
          result = grid_SpliceContour(new_x,new_y,mode=1,spacing=0.20D)   ; parameter
   
          ; Collect point distributions for each section:
          nsection++
          IF (nsection EQ 1) THEN BEGIN
            section_i1 = [i1]
            section_i2 = [i2]
            section_l  = result.length
            section_t  = CREATE_STRUCT('data0', result.t)
          ENDIF ELSE BEGIN
            tag = 'data' + STRING(nsection-1,FORMAT='(I0)')
            section_i1 = [section_i1,i1]
            section_i2 = [section_i2,i2]
            section_l  = [section_l ,result.length]
            section_t  = CREATE_STRUCT(section_t, tag, result.t)
          ENDELSE
         
;          print,'section_i1',section_i1
;          print,'section_i2',section_i2
;          help,section_t,/struct
;          stop   
        ENDIF ELSE BEGIN

          ; Processing separatrices, but looking at the second leg of the 
          ; second separatrix for disconnected double null and a decision
          ; needs to be made at the second x-point about which direction to
          ; take -- toward the far target in the secondary PFR, or to the
          ; primary divertor.  Choosing the first occurance of the 2nd
          ; x-point in the list should give the second option, which is the
          ; desired one:
          IF (count1 EQ 2 AND ipass EQ 1) THEN BEGIN
      
            i3 = (i3[0])[0]
           ; help,i3
           ; stop
   
          ENDIF


; LEFT OFF
; something's not right because a point on the farther-out rings that corresponds to the primary x-point is not included...
; either need to enforce orthogonaity near the x-points by doing intersections with the nearest rings, or need to post-process
; the points, sliding them along the rings so that they are more orthogonal to each other..?

          IF (i3 EQ i4) THEN BEGIN

            print,'         ============ path 2 ============'

;            help,section_t,/struct
print, 'i3=', i3
            tag = 'data' + STRING(i3,FORMAT='(I0)')
print, 'tag=',tag
            section = grid_ExtractStructure(section_t,tag)      
;            help,section
;            print,section
            position = section[1:N_ELEMENTS(section)-2]  ; leave off the end points, which are added in grid_SpliceContour

          ENDIF ELSE BEGIN 

            print,'         ============ path 3 ============'

            ; Build a path from the first point to the second:
            j = i3
            k = 0
            length = 0.0
            WHILE (1) DO BEGIN
              print,'j :',j,k,nsection,format='(A,3I6)'

              IF (k EQ 0) THEN list_j = j ELSE list_j = [list_j,j]

              length = length + section_l[j]

              IF (section_i2[j] EQ i2) THEN BREAK

              print,'section_i1   ',section_i1
              print,'section_i2   ',section_i2
              print,'section_i2[j]',section_i2[j]
              help,section_i1
              help,LONG(section_i2[j])

              j = WHERE( section_i1 EQ (section_i2[j])[0] )

              k++
              IF (j EQ -1 OR k EQ nsection) THEN BREAK
            ENDWHILE
    
            print,'j-:',j,k,nsection,format='(A,3I6)'

            IF (section_i2[j] NE i2) THEN BEGIN
              PRINT,'ERROR grid_CreatePoloidalPoints: Bad situation'
              STOP
            ENDIF
            
            print,'list_j',list_j
            print,'length',length
;            help,section_t,/struct
            shift = 0.0D
            FOR j = 0, N_ELEMENTS(list_j)-1 DO BEGIN

              frac = section_l[list_j[j]] / length

              tag = 'data' + STRING(list_j[j],FORMAT='(I0)')
              print, 'tag ',tag
              section = grid_ExtractStructure(section_t,tag)      
              help,section
;              print,section
              print,'  ---',section_l[list_j[j]],frac,tag
              n = N_ELEMENTS(section)
              print,'section',n,section
              print,'section',n,section[1:n-2]
              print,'section',n,section[1:n-2]*2.0

              IF (j EQ 0) THEN position =            frac[0] * section[1:n-2] ELSE  $
                               position = [position, frac[0] * section[1:n-2] + shift]    

              shift = shift + frac[0] 
              print,'position',position

              stopper = stopper + 1
            ENDFOR

          ENDELSE

          result = grid_SpliceContour(new_x,new_y,position=position)  
        ENDELSE
   
;        help,grid,/struct
;        print,grid.t
;        oplot,grid.x,grid.y,color=Truecolor('Lightblue'),PSYM=6
;        stop
   
        oplot,result.x,result.y,color=Truecolor('Lightblue'),PSYM=6

        if (stopper EQ 3) THEN stop

        IF (isec EQ 0) THEN BEGIN
          grid_x = result.x
          grid_y = result.y
        ENDIF ELSE BEGIn
          nx = N_ELEMENTS(result.x) - 1
          grid_x = [grid_x,result.x[1:nx]] ; drop the first point for all segments except the first
          grid_y = [grid_y,result.y[1:nx]]
        ENDELSE
   
   
      ENDFOR  ; isec

;      oplot,grid_x,grid_y,color=Truecolor('Lightblue'),PSYM=6

    ENDFOR  ; ictr
print, 'stopping after pass 1' 
stop

  ENDFOR  ; ipass


  result = c_array

  RETURN, result

END


