;
; ======================================================================
;
; ======================================================================
;
FUNCTION grid_SpliceContour, val_x, val_y,  $
                             mode     = mode    ,  $
                             spacing  = spacing ,  $
                             position = position

;  IF (1 EQ 1) THEN BEGIN
  IF (KEYWORD_SET(position)) THEN BEGIN

    print,'position > ',position

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
        n    = MAX([3,LONG(dist[N_ELEMENTS(dist)-1] / spacing) + 1])  ; parameter
        step = 1.0D / (DOUBLE(n-1))
        print,'n',n
        print,'step',step,dist[N_ELEMENTS(dist)-1]
        position = [step]
        FOR i = 1, n-3 DO position = [position,position[i-1]+step]
        print,'position',position
stop
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
;print,'val_x',val_x
;print,'val_y',val_y

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
;print,'frac',frac
;print,'dist',dist
      FOR k = 0, N_ELEMENTS(frac)-2 DO  $
        IF (frac[k] LE position[i] AND frac[k+1] GT position[i]) THEN BREAK  
      IF (k EQ N_ELEMENTS(frac)-1) THEN BEGIN
        PRINT,'ERROR grid_SpliceContour: Unable to locate position'
        PRINT,'  position ',position
        PRINT,'  frac     ',frac
        STOP
      ENDIF
      IF (j LT 1) THEN grid_RefineContour, x, y, k
    ENDFOR



    ; Linearly interpolate to extract the point:  (*** should really use a SPLINE ***)
    frac = (position[i] - frac[k]) / (frac[k+1] - frac[k])    

;print, 'i,k',i,k,frac
;print,'x',x
;print,'y',y

;    print,frac
    new_x = x[k] + frac * (x[k+1] - x[k])
    new_y = y[k] + frac * (y[k+1] - y[k])

    result_x = [result_x,new_x]
    result_y = [result_y,new_y]
  ENDFOR

  ; And, finally, add the end point:
  result_x = [result_x,val_x[N_ELEMENTS(val_x)-1]]
  result_y = [result_y,val_y[N_ELEMENTS(val_y)-1]]

;print, 'x',result_x
;print, 'y',result_y

;  oplot,result_x,result_y,color=Truecolor('Lightblue'),PSYM=6,SYMSIZE=0.5

  result = { x      : result_x                ,  $
             y      : result_y                ,  $
             t      : [0.0, position, 1.0]    ,  $
             length : dist[N_ELEMENTS(dist)-1]}

  RETURN, result

END
;
; ======================================================================
;
; ======================================================================
;
FUNCTION grid_EnforceMidplaneAlignment, contours, wall, b, tubes
  ; ------------------------------------------------------------------
  COMMON commons, param, option, state
  ; ------------------------------------------------------------------

  debug  = option.debug
  xrange = option.xrange
  yrange = option.yrange

  geometry = state.geometry

  tags = STRUPCASE(TAG_NAMES(contours))
  nctr = N_ELEMENTS(tags)

  ; Get separatrix contour:
  ctr = grid_ExtractStructure(contours,'contour1')  

  i = WHERE( ctr.grid_x EQ b.x[b.null_i[1]] AND  $
             ctr.grid_y EQ b.y[b.null_j[1]], count)
  IF (count NE 2) THEN BEGIN
    PRINT,'ERROR grid_EnforceMidplaneAlignment: Cannot find primary x-point points'
    STOP
  ENDIF

  grid_s = ctr.grid_s
  grid_s[i] = 0  ; get rid of the primary x-point markers, which don't appear on other contours

  x = ctr.grid_x[i[0]:i[1]]
  y = ctr.grid_y[i[0]:i[1]]

  i = WHERE(x LT b.x[b.null_i[1]])
  inner_x = x[i]
  inner_y = y[i]
  dummy = MIN(ABS(inner_x),inner_i)  ; 'midplane' is the maximum inner point
  i = WHERE(x GT b.x[b.null_i[1]])
  outer_x = x[i]
  outer_y = y[i]
  dummy = MAX(ABS(outer_x),outer_i)  ; 'midplane' is the maxium outer point

  ; Set the y-values for the midplane, based on the separatrix points:
  mid_inner_y = inner_y(inner_i)
  mid_outer_y = outer_y(outer_i)

  mid_sol =          WHERE(ctr.grid_x EQ inner_x[inner_i] AND ctr.grid_y EQ inner_y[inner_i])
  mid_sol = [mid_sol,WHERE(ctr.grid_x EQ outer_x[outer_i] AND ctr.grid_y EQ outer_y[outer_i])]

  print,inner_x[inner_i],inner_y[inner_i]
  print,outer_x[outer_i],outer_y[outer_i]

  print,ctr.grid_x[mid_sol[0]],ctr.grid_y[mid_sol[0]]
  print,ctr.grid_x[mid_sol[1]],ctr.grid_y[mid_sol[1]]

  mid_lfs_range = [-999,-999,-999,-999]
  mid_hfs_range = [-999,-999,-999,-999]

  FOR ipass = 1, 2 DO BEGIN ; do separatrices first, as per usual

    FOR iside = 0, 1 DO BEGIN  ; inner (HFS) then outer (LFS)

      i = mid_sol[iside]

      ; Scan forward and backward to find the marker point span:
      FOR i1 = i-1,                    0, -1 DO IF (grid_s[i1] NE 0) THEN BREAK  
      FOR i2 = i  , N_ELEMENTS(grid_s)-1     DO IF (grid_s[i2] NE 0) THEN BREAK  
      IF (i1 EQ -1 OR i2 EQ N_ELEMENTS(grid_s)) THEN BEGIN
        PRINT,'ERROR grid_EnforceMidplaneAlignment: Unexpected this is...'
        STOP
      ENDIF

      mid_sol_range = [grid_s[i1],grid_s[i2],i-i1,-999]

      IF (ipass EQ 2 AND iside EQ 0 AND mid_hfs_range[0] EQ -999) THEN  $
        mid_hfs_range = mid_sol_range
      IF (ipass EQ 2 AND iside EQ 1 AND mid_lfs_range[0] EQ -999) THEN  $
        mid_lfs_range = mid_sol_range


      print,'--------side-------',iside
      print,'mid_sol_range',mid_sol_range
      print,'grid_s',grid_s[WHERE(grid_s NE 0)]
      print,'searching',i,i1,i2,iside
      print,'pair',grid_s[i1],grid_s[i2]


      ; Search over contours to find the same sequence of cut points:
      FOR ictr = 2, nctr DO BEGIN
        ctr = grid_ExtractStructure(contours,tags[ictr-1])

        print, '-------------------------------ictr--------------------------------',ictr,ipass,iside,ctr.separatrix,ctr.region

        IF ((ipass EQ 1 AND ctr.separatrix EQ 0) OR  $
            (ipass EQ 2 AND ctr.separatrix NE 0)) THEN CONTINUE

        IF (ctr.region GE param.PFZ) THEN BEGIN
           PRINT,'************************* SKIPPING PFZ ',ctr.region,ictr,' ***************************'
           CONTINUE
        ENDIF

        IF ((iside EQ 0 AND ctr.region EQ param.SOL_LFS) OR  $
            (iside EQ 1 AND ctr.region EQ param.SOL_HFS)) THEN BEGIN
           PRINT,'************************* SKIPPING WRONG SIDE ',ctr.region,ictr,' ***************************'

           CONTINUE
        ENDIF

        ; Check if the contour region actually crosses the midplane:
        IF ((iside EQ 0 AND (MIN(ctr.grid_y) GT mid_inner_y OR MAX(ctr.grid_y) LT mid_inner_y)) OR  $
            (iside EQ 1 AND (MIN(ctr.grid_y) GT mid_outer_y OR MAX(ctr.grid_y) LT mid_outer_y)))  THEN BEGIN
          PRINT,'************************* SKIPPING ABOVE/BELOW ',ctr.region,ictr,' ***************************'
          CONTINUE
        ENDIF

        print,'damn',iside EQ 0,MIN(ctr.grid_y) GT mid_inner_y,MAX(ctr.grid_y) LT mid_inner_y

        IF (ipass EQ 1) THEN BEGIN
          oplot,ctr.grid_x,ctr.grid_y,color=Truecolor('Orange')
        ENDIF

        IF (ipass EQ 2) THEN BEGIN
          CASE ctr.region OF
             param.SOL    : mid_range = mid_sol_range
             param.SOL_LFS: mid_range = mid_lfs_range
             param.SOL_HFS: mid_range = mid_hfs_range
          ENDCASE          
          IF (mid_range[0] EQ -999) THEN BEGIN
            PRINT,'************************* DATA NOT SET, SKIPPING ',ctr.region,ictr,' ***************************'
            CONTINUE
          ENDIF           
        ENDIF ELSE mid_range = mid_sol_range

        j1 = WHERE(ctr.grid_s EQ mid_range[0], count1)
        j2 = WHERE(ctr.grid_s EQ mid_range[1], count2)

        print,'count1,2 A=',count1,count2,ipass

        ; Skip to next contour if the marker points are not found:
        IF (count1 EQ 0 AND count2 EQ 0) THEN BEGIN
          PRINT,'************************* SKIPPING NO MARKERS FOUND ',ctr.region,ictr,' ***************************'
          print,'grid_s',ctr.grid_s
          print,'range ',mid_range[0:1]
          print,'mid_lfs_range',mid_lfs_range
          print,'mid_hfs_range',mid_hfs_range
          CONTINUE
        ENDIF

        IF ( ipass EQ 2 AND count1 EQ 0 ) THEN BEGIN
          PRINT,'************************* BHAM! J1 BASED ON DELTA! ',ctr.region,ictr,ipass,' ***************************'
          print, 'count1,2 B= ',count1,count2,mid_range[0:1]
          IF (mid_range[3] EQ -999) THEN BEGIN
            PRINT,'************************* j1 : not fix defined : skipping  ***************************'
            CONTINUE
          ENDIF
          j1 = (j2)[0] - mid_range[3]
          print, 'oh boy, j1,2=',j1,j2,mid_range[3]

          IF (j1 LT 0) THEN BEGIN
            oplot,ctr.grid_x,ctr.grid_y,color=Truecolor('Darkgrey'),PSYM=6
            print,'_s=',ctr.grid_s
;            stop
            PRINT,'************************* j1 < 0 skipping  ***************************'
            CONTINUE
          ENDIF 
        ENDIF

        IF ( ipass EQ 2 AND count2 EQ 0 ) THEN BEGIN
          PRINT,'************************* BHAM! J2 BASED ON DELTA! ',ctr.region,ictr,ipass,' ***************************'
          print, 'count1,2 C= ',count1,count2,mid_range[0:1]
          print, 'j1,2,n    = ',j1,j2,N_ELEMENTS(ctr.grid_x)-1
          IF (j1 EQ N_ELEMENTS(ctr.grid_x)-1) THEN BEGIN
            PRINT,'************************* j2 : end point detected : skipping  ***************************'
            CONTINUE
          ENDIF 
          IF (mid_range[3] EQ -999) THEN BEGIN
            PRINT,'************************* j2 : no fix defined : skipping  ***************************'
            oplot,ctr.grid_x,ctr.grid_y,color=Truecolor('Green'),PSYM=6
            PRINT,'AW, DAMN (2)'
            STOP
            CONTINUE
          ENDIF
          j2 = (j1)[0] + mid_range[3]
          print, 'oh boy, j1,2=',j1,j2,mid_range[3]
          IF (j2 GT N_ELEMENTS(ctr.grid_x)-1) THEN BEGIN
            oplot,ctr.grid_x,ctr.grid_y,color=Truecolor('Darkgrey'),PSYM=6
            print,'_s=',ctr.grid_s
;            stop
            PRINT,'************************* j2 < N_MAX skipping  ***************************'
            CONTINUE
          ENDIF 
        ENDIF

        j1 = (j1)[0]
        j2 = (j2)[0]

        ; Check that there are no marker points between these ones, which 
        ; shouldn't happen, but need to be sure:
        FOR j = j1+1, j2-1 DO   $
          IF (ctr.grid_s[j] NE 0) THEN BEGIN
            PRINT,'ERROR grid_EnforceMidplaneAlignment: Inconsistent sequence'
            STOP
          ENDIF

        j = mid_range[2]

        print,'pair',mid_range[0:2]
        print,'j1',j1
        print,'j2',j2 
        print,'j ',j

        print,ctr.grid_y[j]
        print,'ctr.grid_s',ctr.grid_s[WHERE(ctr.grid_s NE 0)]

        ; Store the separatrix values:
        IF (ipass EQ 1) THEN BEGIN

          IF ((b.y[b.null_j[1]] LT 0.0D AND ctr.separatrix EQ 4) OR  $
              (b.y[b.null_j[1]] GT 0.0D AND ctr.separatrix EQ 3)) THEN  $   ; is this correct 3 always comes first?
            region = param.SOL_LFS ELSE region = param.SOL_HFS

          print,'region',region,ctr.separatrix

          CASE ctr.separatrix OF
;           ------------------------------------------------------------
            2:        ; secondary x-point is on the wall
            3: BEGIN  ; secondary x-point is inside the vessel, the first contour segment   

              xpt_i = WHERE(ctr.grid_x EQ b.x[b.null_i[2]] AND  $
                            ctr.grid_y EQ b.y[b.null_j[2]], count)
              IF (count EQ 0) THEN BEGIN 
                PRINT,'ERROR grid_EnforceMidplaneAlignment: Secondary x-point not found (1)'
                STOP               
              ENDIF
              ; 
              ; --------------------------------------------------------
              ; CHECK IF THE POINT j2 CORRESPONDS TO THE SECONDARY X-
              ; POINT, AND IF YES, REMAP TO THE NEXT TANGENCY POINT 
              ; ALONG THE CONTOUR.  THIS IS A PROBLEM SINCE THE 
              ; SECONDARY X-POINT DOES NOT APPEAR AS A TANGENCY
              IF (0 EQ 1 AND j2 EQ xpt_i) THEN BEGIN
                print, '!!! do some damn work !!! (1)'
                arop
                FOR j3 = j2+1, N_ELEMENTS(ctr.grid_s)-1 DO IF (ctr.grid_s[j3] NE 0) THEN BREAK  
                IF (j3 EQ N_ELEMENTS(ctr.grid_s)) THEN BEGIN
                  PRINT,'ERROR grid_EnforceMidplaneAlignment: Search failed (1)'
                  STOP
                ENDIF
              ENDIF ELSE j3 = j2

;              new_range = mid_sol_range
;              new_range = [ctr.grid_s[j1],ctr.grid_s[j3],j,-999]
              new_range = [ctr.grid_s[j1],ctr.grid_s[j3],j,j3-j1]

              IF (region EQ param.SOL_HFS) THEN  $
                 mid_hfs_range = new_range ELSE  $
                 mid_lfs_range = new_range

              END
;           ------------------------------------------------------------
            4: BEGIN  ; secondary x-point is inside the vessel, the second contour segment   

              xpt_i = WHERE(ctr.grid_x EQ b.x[b.null_i[2]] AND  $
                            ctr.grid_y EQ b.y[b.null_j[2]], count)
              IF (count EQ 0) THEN BEGIN 
                PRINT,'ERROR grid_EnforceMidplaneAlignment: Secondary x-point not found (2)'
                STOP               
              ENDIF
              ; 
              ; --------------------------------------------------------
              ; CHECK IF THE POINT j1 CORRESPONDS TO THE SECONDARY X-
              ; POINT, AND IF YES, REMAP TO THE NEXT TANGENCY POINT 
              ; EARLIER ALONG THE CONTOUR.
              ;
              IF (0 EQ 1 AND j1 EQ xpt_i) THEN BEGIN
                stop
              ENDIF ELSE j3 = j1
              ; new_range = mid_sol_range
              new_range = [ctr.grid_s[j1],ctr.grid_s[j2],j,j2-j1]
;              new_range = [ctr.grid_s[j3],ctr.grid_s[j2],j,j2-j3]
              IF (region EQ param.SOL_HFS) THEN  $
                 mid_hfs_range = new_range ELSE  $
                 mid_lfs_range = new_range
              END
;           ------------------------------------------------------------
          ENDCASE
        ENDIF


        ; Check if the contour region actually crosses the midplane:
        IF ((iside EQ 0 AND ctr.grid_y[j1] LT mid_inner_y AND ctr.grid_y[j2] LT mid_inner_y) OR  $
            (iside EQ 0 AND ctr.grid_y[j1] GT mid_inner_y AND ctr.grid_y[j2] GT mid_inner_y) OR  $
            (iside EQ 1 AND ctr.grid_y[j1] LT mid_outer_y AND ctr.grid_y[j2] LT mid_outer_y) OR  $
            (iside EQ 1 AND ctr.grid_y[j1] GT mid_outer_y AND ctr.grid_y[j2] GT mid_outer_y)) THEN BEGIN
          PRINT,'************************* SKIPPING ABOVE/BELOW ',ctr.region,ictr,ipass,iside,' ***************************'


          PRINT,'STILL NEEDED -- DAMN!'
;          STOP
;          CONTINUE  ; ***TEMP***
        ENDIF

;
;       ----------------------------------------------------------------
;       MOVE THE "MIDPLANE POINT" TO THE MIDPLANE
;
        IF (j1 EQ 0) THEN BEGIN
;        IF (ctr.grid_s[j1] LT 0) THEN BEGIN
          iwall1 = WHERE(wall.ptc EQ ictr AND wall.ptt EQ 1, count)
          IF (count EQ 0) THEN BEGIN
            PRINT,'ERROR grid_EnforceMidplaneAlignment: Target point not found (1)'
            STOP
          ENDIF
          i1     =  1
          count1 = -1
        ENDIF ELSE  $
          i1 = WHERE(ctr.grid_x[j1] EQ ctr.x AND ctr.grid_y[j1] EQ ctr.y, count1)

        IF (count1 EQ 0) THEN BEGIN
          ; Need to search for the point since it is not a tangency or end point:
          dummy = grid_PerpDistance(ctr.grid_x[j1], ctr.grid_y[j1],  $
                                    ctr.x,ctr.y,s=s,index=index)
          print,'search index 2',index,s,dummy
          i1 = index + 1
          count1 = -99
        ENDIF

        IF (j2 EQ N_ELEMENTS(ctr.grid_s)-1) THEN BEGIN
;        IF (ctr.grid_s[j2] LT 0) THEN BEGIN
          iwall2 = WHERE(wall.ptc EQ ictr AND wall.ptt EQ 2, count)
          IF (count EQ 0) THEN BEGIN
            PRINT,'ERROR grid_EnforceMidplaneAlignment: Target point not found (2)'
            STOP
          ENDIF
          i2     = N_ELEMENTS(ctr.x) - 2
          count2 = -1
        ENDIF ELSE  $
          i2 = WHERE(ctr.grid_x[j2] EQ ctr.x AND ctr.grid_y[j2] EQ ctr.y, count2)

        IF (count2 EQ 0) THEN BEGIN
;          print,'shit -- need to check that this is working properly'
;          stop
          ; Need to search for the point since it is not a tangency or end point:
          dummy = grid_PerpDistance(ctr.grid_x[j2], ctr.grid_y[j2],  $
                                    ctr.x,ctr.y,s=s,index=index)
          print,'search index 2',index,s,dummy
          i2 = index
          count2 = -99
        ENDIF

        IF (count1 EQ -1) THEN BEGIN
          res_x = [wall.pt1[0,iwall1],ctr.x[i1]]  
          res_y = [wall.pt1[1,iwall1],ctr.y[i1]]  
        ENDIF                                     
        IF (count1 EQ -99) THEN BEGIN
          res_x = [ctr.grid_x[j1],ctr.x[i1]]
          res_y = [ctr.grid_y[j1],ctr.y[i1]]
        ENDIF
        IF (count1 GT 0) THEN BEGIN
          res_x = [ctr.x[i1]]
          res_y = [ctr.y[i1]]
        ENDIF

        res_x = [res_x,ctr.x[i1+1:i2]]
        res_y = [res_y,ctr.y[i1+1:i2]]

        IF (count2 EQ -1) THEN BEGIN
          res_x = [res_x,wall.pt1[0,iwall2]]
          res_y = [res_y,wall.pt1[1,iwall2]]
        ENDIF
        IF (count2 EQ -99) THEN BEGIN
          res_x = [res_x,ctr.grid_x[j2]]
          res_y = [res_y,ctr.grid_y[j2]]
        ENDIF


        res_frac = grid_GetFrac(res_x,res_y,dist=res_dist)

print,'coun',count1,count2
print,'numb',i1,i2
;print,'frac',res_frac
print,'dist',res_dist[N_ELEMENTS(res_dist)-1]
oplot,res_x,res_y,psym=3,color=Truecolor('Orange')


        x = ctr.grid_x[j1:j2]
        y = ctr.grid_y[j1:j2]
        frac = grid_GetFrac(x,y,dist=dist)

print,'j',j,ctr.grid_s[j1]
print,'frac',frac
print,'dist',dist[N_ELEMENTS(dist)-1]
print,'midpoint',x[j],y[j]
;print,'below',x[j-1],y[j-1]
;print,'above',x[j+1],y[j+1]
oplot,x,y,PSYM=3
oplot,[x[j]],[y[j]],PSYM=6,COLOR=Truecolor('Orange')



print, 'inner_mid',mid_inner_y,y[j],ctr.region
print, 'outer_mid',mid_outer_y,y[j]


        IF (iside EQ 0) THEN mid_y =  mid_inner_y ELSE  $
                             mid_y =  mid_outer_y

        IF (ABS(y[j] - mid_y) LT 0.001D) THEN BEGIN
          PRINT,'************************* SKIPPING TOO CLOSE ',ctr.region,ictr,' ***************************'
          CONTINUE
        ENDIF

        IF (iside EQ 0) THEN i = WHERE(res_x LT (b.x[b.null_i[1]])[0], count) ELSE  $
                             i = WHERE(res_x GT (b.x[b.null_i[1]])[0], count) 

        IF (count EQ 0) THEN BEGIN
          PRINT,'************************* SKIPPING WRONG SIDE ',ctr.region,ictr,' ***************************'
          print,'this should not go off, yes?'
          CONTINUE
print,'iside',iside
print,'x    ',x
print,'res_x',res_x
print,'xpt_x',(b.x[b.null_i[1]])[0]
          STOP
        ENDIF

;print,res_x
;print,b.x[b.null_i[1]],iside
;print,'i',i
        dummy = MIN(ABS(res_y[i] - mid_y),k)
        i = WHERE(res_x[i[k]] EQ res_x AND res_y[i[k]] EQ res_y)
oplot,[res_x[i]],[res_y[i]],psym=6,color=Truecolor('Cyan') 
print,res_x[i],res_y[i],mid_y
print,'fracs',frac[j],res_frac[i]
print,'j',j
print,'N_ELEMENTS(frac)',N_ELEMENTS(frac)
print,'N_ELEMENTS(frac2)',N_ELEMENTS(frac2)  

        frac1 = frac / frac[j]
        frac2 = (frac    - frac[N_ELEMENTS(frac)-1]) /  $
                (frac[j] - frac[N_ELEMENTS(frac)-1])
        frac_adjust = ABS([frac1[0:j], frac2[j+1:N_ELEMENTS(frac)-1]])

        position = frac + frac_adjust * ((res_frac[i] - frac[j]))[0]

help, frac_adjust
help, frac_adjust * (res_frac[i] - frac[j])
print,'frac    ',frac
print,'position',position

;IF (ipass EQ 2 and iside eq 1) THEN stop

        result = grid_SpliceContour(res_x,res_y,position=position[1:N_ELEMENTS(position)-2])

oplot,result.x,result.y,psym=6,color=Truecolor('Magenta') 

        ; Update the points in the contour:
        ctr.grid_x[j1:j2] = result.x    
        ctr.grid_y[j1:j2] = result.y    
        contours = grid_UpdateStructure(contours,tags[ictr-1],ctr)

help,result,/struct

print,'frac1',frac1
print,'frac2',frac2[j:N_ELEMENTS(frac)-1]
;print,'frac_new',frac_new
print,N_ELEMENTS(frac),N_ELEMENTS(frac_new)

;IF (ipass EQ 2 AND iside EQ 1) THEN stop

      ENDFOR  ; ICTR

    ENDFOR  ; ISIDE

  ENDFOR  ; IPASS

;  stop

  RETURN,contours

  END
;
; ======================================================================
;
; ======================================================================
;
FUNCTION grid_PoloidalDistribution, length
  ; ------------------------------------------------------------------
  COMMON commons, param, option, state
  ; ------------------------------------------------------------------

;
; ----------------------------------------------------------------------
;
;

  res_n      = 10
  res_adjust = 10

;  res_min = 0.0010D ; 0.001D ; 0.001D  ; iter
;  res_max = 0.500D ; 0.020D ; 0.050D
;  res_min = 0.0010D ; 0.001D ; 0.001D  ; d3d
;  res_max = 0.050D ; 0.020D ; 0.050D
;  res_exp = 1.0D

  res_min = option.pol_res_min
  res_max = option.pol_res_max
  res_exp = option.pol_res_exp

  status = 1

  WHILE (status) DO BEGIN


    CASE 1 OF
;   --------------------------------------------------------------------
      1: BEGIN

;        print,'res_n',res_n,res_adjust

        res = DINDGEN(res_n+1) / DOUBLE(res_n)

        ; Apply the distribution function of choice:

;        res = [res[0],(res[1:N_ELEMENTS(res)-1])^0.0]
;        scale = (ALOG(res_min + 1.0D) / res_exp) / res[1]
    
;print, 'res1',res

         res = res(1:res_n-1)

         CASE option.pol_res_opt OF
           param.EXPONENTIAL: res = [EXP(res/res_exp) - 1.0D]
           ELSE: BEGIN 
             PRINT,'ERROR grid_PoloidalDistribution: Unrecognised distribution'
             PRINT,'  OPT = ', option.pol_res_opt
             STOP
             END
         ENDCASE
;        res = [EXP(res) - 1.0D]
;        res = [EXP(res - res_exp) - EXP(-res_exp)] ; no change
;        res = [res^res_exp]
;       plot,res,psym=6

;        x = res * 10.0D 
;        shift = -5.0D
;        res = [res[0],( EXP(0.5D*res_exp*(x+shift)) - 1.0D ) /  $
;                      ( EXP(0.5D*res_exp*(x+shift)) + 1.0D ) ]
;        res = [res[0], (res[1:res_n-1] - res[1]) * (res_max - res_min) / MAX(res) + res_min]
;print,'res3',res

;plot,[0.0,x],res,psym=6
;stop

        ; Rescale the point nearest the separatrix (RES = 0.0) so that it
        ; is equal to RES_MIN:
;print,'before',res     
        scale = res_min / res[0]
        res = res * scale
;print,'after ',res     
        ; Limit the maxium spatial resolution:
        FOR i = 0, N_ELEMENTS(res)-1 DO res[i] = MIN([res[i],res_max])
;print,'after limited',res     

;plot,res,psym=6
;stop

        ; Mirror:
        res = [res,REVERSE(res)]
;print,'mirror',res
;print,res[izero]
        ; Generate distribution of points in space:
        dist = [0.0D0,res]
        n = N_ELEMENTS(dist)
        FOR i = 1, n-1 DO dist[i] = dist[i-1] + res[i-1]
;print,'dist  ',dist
;i=where(dist GE 0.0D)
;plot,dist,res,psym=6,xrange=[0.0,0.05]

        END
;    ------------------------------------------------------------------
    ENDCASE

; print,'comparison',dist[n-1],length

    IF            (res_adjust GT 0 AND dist[n-1] GT length) THEN BEGIN
      res_adjust = -res_adjust / 2
    ENDIF ELSE IF (res_adjust LT 0 AND dist[n-1] LT length) THEN BEGIN
      res_adjust = -res_adjust / 2
    ENDIF

    IF (ABS(res_adjust) EQ 0) THEN BEGIN
      status = 0 
      dist = dist * (length / (dist[n-1])[0])
    ENDIF ELSE  $
      res_n = res_n + res_adjust

  ENDWHILE

; print,'comparison final',dist[n-1],length

; print,dist[1]-dist[0]
; print,dist[n-1]-dist[n-2]


  result = dist / dist[n-1]

;  print,'dist',dist

  RETURN, result

END
;
; ======================================================================
;
; Analyse the distribution of "tangency" points and 
;
;
;
FUNCTION grid_SetPoloidalDomains, c_array, wall, b, debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON commons, param, option, state
  ; ------------------------------------------------------------------

  IF (debug) THEN BEGIN
    PRINT, ' '
    PRINT, '----------------------------------------------------------------------'
    PRINT, ''
    PRINT, 'SETTING POLOIDAL DOMAINS'
    PRINT, ''
    PRINT, '----------------------------------------------------------------------'
    PRINT, ' '
  ENDIF

  tags  = STRUPCASE(TAG_NAMES(c_array))
  nctr  = N_ELEMENTS(tags)

  IF (KEYWORD_SET(debug)) THEN BEGIN
    PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1, CHARSIZE=1.0
    OPLOT,wall.x,wall.y,color=Truecolor('Black')
    FOR i = 0, nctr-1 DO BEGIN
      ctr = grid_ExtractStructure(c_array,tags[i])      
      print,'contour debug',i,ctr.separatrix
      OPLOT,ctr.x,ctr.y,color=Truecolor('Red')    
      npts = N_ELEMENTS(ctr.x)
      XYOUTS,ctr.x[npts/2],ctr.y[npts/2],STRTRIM(STRING(i+1),2),color=Truecolor('Pink'),CHARSIZE=2.0
    ENDFOR
  ENDIF

  PRINT,'---------- START DOMAIN IDENTIFICATION ------------'

  ; Process tangency points:

  vector_n = 0  

  section = 0

  FOR ipass = 1, 1 DO BEGIN  ; no double pass for now... which I wanted so that the 2nd xpoint othogonal was correct and not in competition with other tangency points

    FOR ictr = 1, nctr DO BEGIN 
      ctr = grid_ExtractStructure(c_array,tags[ictr-1])      

      IF (ictr EQ 1 AND state.geometry EQ param.LIMITED) THEN CONTINUE
    
;      IF ((ipass EQ 1 AND ctr.separatrix EQ 0) OR  $
;          (ipass EQ 2 AND ctr.separatrix NE 0)) THEN CONTINUE

      IF (ctr.tangent_i EQ -1 AND  $
          (ctr.separatrix LT 3 OR ctr.tangent_p1[0] EQ 0.0)) THEN CONTINUE
    
      print, ' =================== ictr ================= ',ictr,' ',tags[ictr-1]

    
      OPLOT,ctr.x,ctr.y,color=Truecolor('Pink')    
      OPLOT,[ctr.tangent_p1[0]],[ctr.tangent_p1[1]],color=Truecolor('Lightgreen'),PSYM=6    
    
      ; Follow the trajectory into the separatrix:
      p1 =        ctr.tangent_p1
      p2 = 3.0D * ctr.tangent_p1 - 2.0D * ctr.tangent_p2  ; dev
  ;    p2 = 2.0D * ctr.tangent_p1 - ctr.tangent_p2  
    
      OPLOT,[p1[0],p2[0]],[p1[1],p2[1]],color=Truecolor('Cyan')    
    
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
                                                          'section_y',p1[1])  ; 'section_y',p1[1]) -bug, SL, 12/10/2011
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
        
        print, 'round =========>',cnt,ictr3,ctr.separatrix
    
        ; Not trying to use the connection map here -- why not do you think?
        min_ictr2 = -1
        min_ctr2  = -1
        min_inter = -1
        min_s     = 1.0D+06
        FOR ictr2 = nctr, 1, -1 DO BEGIN

          IF (ictr2 EQ ictr3) THEN CONTINUE

          ctr2 = grid_ExtractStructure(c_array,tags[ictr2-1])      

          ; Special case where processing the secondary separatrix
          ; intersections:
          IF (ctr.separatrix GE 3 AND ctr2.separatrix GE 3) THEN CONTINUE  

          ; Only look at rings with larger PSI values:
          IF (ctr.region LT param.PFZ AND ctr.psi GT ctr2.psi) THEN CONTINUE

          n = N_ELEMENTS(ctr2.x)

          IF (n LE 5) THEN BEGIN
            print,'forcing enhanced resolution on a contour -- super lame',n
            stop
            x = ctr2.x
            y = ctr2.y
            grid_RefineContour, x, y, n/2
            ctr2    = grid_UpdateStructure(ctr2,'x',x)
            ctr2    = grid_UpdateStructure(ctr2,'y',y)            
            c_array = grid_UpdateStructure(c_array,tags[ictr2-1],ctr2)
            grid_UpdateWall, b, c_array, wall, $
                             debug=debug, xrange=xrange, yrange=yrange
            c_array = grid_BuildRadialMap(c_array, wall,  $
                                          debug=debug, xrange=xrange, yrange=yrange)
            n = N_ELEMENTS(ctr2.x)
            print,'n is ',n
          ENDIF 

          ishift = 1
          seg1 = MAKE_ARRAY(2,n-2*ishift-1,/DOUBLE,VALUE=0.0D)
          seg2 = seg1
          seg1[0,*] = ctr2.x[ishift  :n-ishift-2]
          seg1[1,*] = ctr2.y[ishift  :n-ishift-2]
          seg2[0,*] = ctr2.x[ishift+1:n-ishift-1]
          seg2[1,*] = ctr2.y[ishift+1:n-ishift-1]

          speed_res = 5 ; 50  ! parameter

          IF (n GT 2*speed_res) THEN BEGIN

            status1 = 0

            segn = N_ELEMENTS(seg1[0,*])
            segi = (INDGEN((segn-1)/speed_res) + 1)* speed_res - 1
;            segi = [0,segi]
            segi = [0,segi,segn-1]

            j = [-1]
            SWITCH ctr2.separatrix OF 
;             ----------------------------------------------------------
              0: BREAK
;             ----------------------------------------------------------
              1: BEGIN  ; primary separatrix
                IF (state.geometry EQ param.LIMITED) THEN BREAK

                j = WHERE( seg1[0,*] EQ b.x[b.null_i[1]] AND seg1[1,*] EQ b.y[b.null_j[1]], count)
                IF (count NE 2) THEN BEGIN
                  PRINT,'ERROR grid_SetPoloidalDomains: Primary x-point not found'
                  STOP
                ENDIF

                BREAK
                END
;             ----------------------------------------------------------
              2: BEGIN  ; secondary, x-point on vessel wall
                proximity = SQRT ( (seg1[0,*]-b.x[b.null_i[2]])^2 +  $
                                   (seg1[1,*]-b.y[b.null_j[2]])^2 )
                dummy = MIN(proximity,j)
                j = [j]

                BREAK
                END
;             ----------------------------------------------------------
              3:        ; secondary, HFS and LFS legs
              4: BEGIN
                j = WHERE( seg1[0,*] EQ b.x[b.null_i[2]] AND seg1[1,*] EQ b.y[b.null_j[2]], count)
                IF (count NE 1) THEN BEGIN
                  PRINT,'ERROR grid_SetPoloidalDomains: Primary x-point not found'
                  STOP
                ENDIF
                j = [j]
                BREAK
                END
;             ----------------------------------------------------------
              ELSE: BEGIN
                PRINT,'ERROR grid_SetPoloidalDomains: Unknown separatrix designation'
                STOP
                END
;             ----------------------------------------------------------
            ENDSWITCH

            ; Place the x-point indeces in the list, if necessary:
            IF (j[0] NE -1) THEN BEGIN
              FOR k = 0, N_ELEMENTS(j)-1 DO BEGIN
                FOR l = 0, N_ELEMENTS(segi)-2 DO BEGIN
                  IF (segi[l] LT j[k] AND segi[l+1] GT j[k]) THEN BEGIN 
                    segi = [segi[0:l],j[k],segi[l+1:N_ELEMENTS(segi)-1]]
                    BREAK
                  ENDIF
                ENDFOR
              ENDFOR  
            ENDIF

            seg3 =   seg1[*,segi[0:N_ELEMENTS(segi)-2]]
            seg4 = [[seg1[*,segi[1:N_ELEMENTS(segi)-2]]],[seg2[*,segn-1]]]

            inter = grid_Intersection(p1, p2, seg3, seg4, 0, status=status)
            IF (0 EQ 1 AND NOT status AND status1) THEN begin
              for k = 0, MIN([10,N_ELEMENTS(seg1[0,*])])-1 DO Begin
                print, 'seg1,2: ,',k,seg1[0:1,k],seg2[0:1,k],format='(A,I6,2(2F10.4,2X))'
              endfor
              for k = MAX([0,N_ELEMENTS(seg1[0,*])-10]), N_ELEMENTS(seg1[0,*])-1 DO Begin
                print, 'seg1,2: ,',k,seg1[0:1,k],seg2[0:1,k],format='(A,I6,2(2F10.4,2X))'
              endfor
              for k = 0, N_ELEMENTS(seg3[0,*])-1 DO Begin
                print, 'seg3,4: ,',k,seg3[0:1,k],seg4[0:1,k],format='(A,I6,2(2F10.4,2X))'
              endfor
              inter = grid_Intersection(p1, p2, seg3, seg4, 0, status=status, /show)
              STOP
            endif

          ENDIF ELSE  $
            inter = grid_Intersection(p1, p2, seg1, seg2, 0, status=status)

          IF (status) THEN BEGIN

            IF (n GT 2*speed_res) THEN BEGIN
              ; Find the actual intersection:

              IF (N_ELEMENTS(inter.i) GT 1) THEN BEGIN
                dummy = MIN(inter.s,imin)
                inter.i = inter.i[imin]
                inter.s = inter.s[imin]
                inter.t = inter.t[imin]
                inter.x = inter.x[imin]
                inter.y = inter.y[imin]
              ENDIF

              segj = (inter.i[0])[0]

              seg3_save = seg3
              seg4_save = seg4

              seg3 = seg1[*,segi[segj]:segi[segj+1]]
              seg4 = seg2[*,segi[segj]:segi[segj+1]]

              inter = grid_Intersection(p1, p2, seg3, seg4, 0, status=status)
              IF (NOT status) THEN BEGIN
                ; A false positive, due to the approximate nature of the initial check:
                inter = { s : min_s + 1.0D }
              ENDIF ELSE BEGIN
                inter.i = inter.i + segi[segj]
              ENDELSE
            ENDIF

            IF (MIN(inter.s) LT min_s) THEN BEGIN
              inter.i = inter.i + ishift
              min_inter = inter 
              min_s     = MIN(inter.s)
              min_ictr2 = ictr2
              min_ctr2  = ctr2
            ENDIF
          ENDIF
        ENDFOR
        IF (min_s EQ 1.0D+6) THEN BEGIN
;   
;         ----------------------------------------------------------------
;         NO INTERSECTIONS FOUND AT ALL, SO FORCE ONE (RISKY)
;   
          IF (b.y[b.null_j[1]] LT 0.0 AND p1[1] LT 0.0) THEN BEGIN
            x = ctr2.x
            y = ctr2.y
            dist = 1.0D+6
            FOR i = 0, N_ELEMENTS(x)-1 DO  $
              dist = [dist,grid_PerpDistance(x[i],y[i],[p1[0],p2[0]],[p1[1],p2[1]],/nostop)]
            dist = dist[1:N_ELEMENTS(dist)-1]
            dummy = MIN(dist,i)
            IF (i LT N_ELEMENTS(x)/2) THEN i=MAX([i,3]) ELSE i=MIN([i-1,N_ELEMENTS(x)-4])
            new_x = MEAN([x[i],x[i+1]])
            new_y = MEAN([y[i],y[i+1]])
            p2 = [new_x,new_y]
            inter  = {          $
              i    : [i  ]   ,  $
              s    : [0.0]   ,  $
              t    : [0.0]   ,  $ 
              x    : [new_x] ,  $
              y    : [new_y] }
            print,'having to fake it!',i,ictr2
            min_ictr2 = 1
            min_ctr2  = ctr2
            min_inter = inter
          ENDIF ELSE BEGIN
            PRINT,'ERROR grid_SetPoloidalDomains: Strange animal'
            STOP
          ENDELSE
        ENDIF

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
    
            v1 = p2                  - p1 
            v2 = vector_seg2[*,ivec] - vector_seg1[*,ivec]
            len1 = SQRT( v1[0]^2 + v1[1]^2 )
            len2 = SQRT( v2[0]^2 + v2[1]^2 )
    
            v2 = (len1 / len2) * v2
    
            frac = (vector_inter.t[0]) * 0.7D   ; parameter
    
            v3 = frac * v1 + (1.0D - frac) * v2
    
            p2 = p1 + v3
    
            ;print,vector_inter.s
            ;print,vector_inter.t
            ;help,frac
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
    
        OPLOT,[inter_x],[inter_y],color=Truecolor('Darkgrey'),PSYM=6
        OPLOT,[ctr2.x[inter_i  ]],[ctr2.y[inter_i  ]],color=Truecolor('Red'),PSYM=6
        OPLOT,[ctr2.x[inter_i+1]],[ctr2.y[inter_i+1]],color=Truecolor('Red'),PSYM=6
    
        ; Add point to the contour and register the designation of the 
        ; tangency point:
    
        ;   Note that this intersection point is the result of a linear
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
        IF (ctr.region GE param.PFZ) THEN direction = 2 ELSE direction = 1

        p2_new = grid_GetOrthogonal(ctr2.x,ctr2.y,p1_new,direction,1.5D,ictr=ictr2,ctrs=c_array)

        ; Also quit if a leg of the secondary separatrix is encountered in the secondary PFZ,
        ; so that an intersection with the primary separatrix will not occur:
        IF (p2_new[0] EQ -999.0D) THEN BREAK
    
        OPLOT,[p1_new[0],p2_new[0]],[p1_new[1],p2_new[1]],color=Truecolor('Blue')
    
        cnt++
    
        p1 = p1_new
        p2 = p2_new              
    
        ictr3 = ictr2
      ENDWHILE
    
    ENDFOR  ; ICTR loop

  ENDFOR  ; IPASS loop



  IF (state.geometry NE param.LIMITED) THEN BEGIN

    print,'MODDING THE PRIMARY X-POINT--------------------'

    ; Need to add two more slice points at the primary x-point:
    ictr = grid_GetIndex(c_array,'separatrix',1)
    ctr  = grid_ExtractStructure(c_array,tags[ictr-1])        
    focus_x = ctr.null_x
    focus_y = ctr.null_y
    i = WHERE(ctr.x EQ focus_x AND ctr.y EQ focus_y)
    section++
    IF (N_ELEMENTS(ctr2) EQ 0) THEN count = 0 ELSE  $
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
  ENDIF

  IF (debug) THEN BEGIN
    PRINT,' '
    PRINT,'------------------------------------------------------------------------'
    PRINT,'SORTING THE SECTIONS'
    PRINT,'------------------------------------------------------------------------'
    PRINT,' '
  ENDIF

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
        PRINT,'  ICTR     =',ictr
        PRINT,'  I        =',i
        PRINT,'  N        =',n
        PRINT,'  SECTION LIST  =',section
        PRINT,'  SECTION       =',ctr.section
        PRINT,'  SECTION_X     =',ctr.section_x[i]
        PRINT,'  SECTION_Y     =',ctr.section_y[i]
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





  ; Set end points codes for secondary separatrices:
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
;
; ----------------------------------------------------------------------
; ASSIGN THE END POINTS TO EACH CONTOUR SECTION LIST
;
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

      IF (ctr.region GE param.PFZ) THEN ictr2 = ctr.map_out[0] ELSE  $
                                        ictr2 = ctr.map_in [0]

      ctr2  = grid_ExtractStructure(c_array,tags[ictr2-1])

      print,'  region,map=',ictr2,ctr.region,ctr2.map_in,ctr2.map_out,ctr2.section[0],format='(a,8I6)'

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
        IF (ctr.region GE param.PFZ) THEN  $
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
;
; ----------------------------------------------------------------------
; SET THE ...
;
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
      IF (ctr.region LT param.PFZ) THEN CONTINUE

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

;
; ----------------------------------------------------------------------
;
;
  ; Need to adjust the end values in the SOL when there's a secondary
  ; x-point

  FOR ictr = 1, nctr DO BEGIN 
    ctr = grid_ExtractStructure(c_array,tags[ictr-1])      

    print, 'ictr',ictr,ctr.separatrix

    IF (ctr.separatrix LT 3) THEN CONTINUE

    IF (ctr.separatrix EQ 3) THEN end_point = 0 ELSE end_point = 1
 
    CASE end_point OF
;     -------------------------------------------------------------------
      0: BEGIN

        ictr2 = ctr.map_in[end_point] 
        ctr2  = grid_ExtractStructure(c_array,tags[ictr2-1])      

        location = -1
        IF (ctr2.map_out[0] EQ ictr) THEN location = 1
        IF (ctr2.map_out[1] EQ ictr) THEN location = 2
        IF (location EQ -1) THEN BEGIN
          ; This could either be a problem with the connection map or the unusual
          ; case where a tangency point on CTR2 is adjacent to CTR, so that CTR2 
          ; doesn't map to CTR.  This can happen when the inner wall contact is just
          ; inside the secondary separatrix, which is nasty.  So, check for a tangency
          ; point on CTR2 which is next to CTR on the wall, otherwise register an error:
          IF (1 EQ 1) THEN BEGIN
            i = WHERE(ctr2.section_x EQ ctr2.tangent_p1[0] AND  $
                      ctr2.section_y EQ ctr2.tangent_p1[1], count)
            IF (count EQ 0) THEN BEGIN
              PRINT,'ERROR grid_SetPoloidalDomains: Unable to find tangency point (1)'
              PRINT,'  ICTR2=',ictr2
              STOP
            ENDIF
          ENDIF ELSE BEGIN
            PRINT,'ERROR grid_SetPoloidalDomains: Something wrong with radial mapping (2)'
            PRINT,' ICTR =',ictr
            PRINT,' ICTR2=',ictr2
            oplot,ctr.x ,ctr.y ,COLOR=Truecolor('Green')
            oplot,ctr2.x,ctr2.y,COLOR=Truecolor('Orange')
            STOP
          ENDELSE
        ENDIF ELSE BEGIN
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

        ENDELSE    

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

       END
;     -------------------------------------------------------------------
      1: BEGIN

        ; Scan along the wall until either a tangency point or the primary separatrix is reached:
 
        iwall = WHERE(wall.ptc EQ ictr AND wall.ptt EQ 2, count)
        IF (count EQ 0) THEN BEGIN
          PRINT,'ERROR grid_SetPoloidalDomains: Cannot locate secondary separatrix on the wall'
          STOP
        ENDIF
        
        new_value = -999 
        FOR idummy = 0, N_ELEMENTS(wall.ptc)-2 DO BEGIN
          iwall++
          IF (iwall GT N_ELEMENTS(wall.ptc)-1) THEN iwall = 0
;          print,iwall,wall.ptc[iwall],wall.ptt[iwall],FORMAT='(3I6)'

          IF (wall.ptt[iwall] EQ 3) THEN BEGIN
            ; Found a tangency point, and now need to find out which section index
            ; this tangency point corresponds to:
            ictr2 = (wall.ptc[iwall])[0]
            ctr2 = grid_ExtractStructure(c_array,tags[ictr2-1])      
            i = WHERE(ABS(ctr2.section_x-(wall.pt1[0,iwall])[0]) LT 1.0D-10 AND  $
                      ABS(ctr2.section_y-(wall.pt1[1,iwall])[0]) LT 1.0D-10, count)
            IF (count EQ 0) THEN BEGIN
              PRINT,'ERROR grid_SetPoloidalDomains: Cannot find tangency point'
              PRINT,'  SECTION_X = ',ctr2.section_x
              PRINT,'  SECTION_Y = ',ctr2.section_y
              PRINT,'  WALL_X    = ',wall.pt1[0,iwall]
              PRINT,'  WALL_Y    = ',wall.pt1[1,iwall]
              STOP
            ENDIF
            new_value = ctr2.section[i]
          ENDIF ELSE  $
          IF (wall.ptt[iwall] EQ 2 AND wall.ptc[iwall] EQ 1) THEN BEGIN
            ; Found the primary separatrix:
            new_value = -2
          ENDIF
          IF (new_value NE -999) THEN BREAK
        ENDFOR
        IF (new_value EQ -999) THEN BEGIN
          PRINT,'ERROR grid_SetPoloidalDomains: Unable to identify end point mapping'
          STOP
        ENDIF

        IF (end_point EQ 0) THEN old_value = ctr.section[0                        ] ELSE  $
                                 old_value = ctr.section[N_ELEMENTS(ctr.section)-1]
   
        ; Scan all other rings and update any references to the old
        ; target index:
        FOR ictr2 = 1, nctr DO BEGIN 
          ctr2 = grid_ExtractStructure(c_array,tags[ictr2-1])      
          i = WHERE(ctr2.section EQ old_value, count)
          IF (count EQ 0) THEN CONTINUE
          print,'updating',ictr2,i,new_value
          ctr2.section[i] = new_value
          c_array = grid_UpdateStructure(c_array,tags[ictr2-1],ctr2)
        ENDFOR

        END
;     -------------------------------------------------------------------
    ENDCASE


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
; ======================================================================
;
PRO grid_ChooseIndex, geometry, b, ctr, section_l, j, i2, xpoint1_i1
  ; ------------------------------------------------------------------
  COMMON commons, param, option, state
  ; ------------------------------------------------------------------

  ; Set stuff:  *** PUT IN A SUBROUTINE, STORE IN b STRUCTURE ***
  xpoint1_psi = b.psi_1st_xpoint
  xpoint1_x   = b.x[b.null_i[1]]
  xpoint1_y   = b.y[b.null_j[1]]
  IF (N_ELEMENTS(b.null_i) EQ 2) THEN BEGIN
    xpoint2_psi = -1.0D+6
    xpoint2_x   =  1.0D+6
    xpoint2_y   =  1.0D+6
  ENDIF ELSE BEGIN
    xpoint2_psi = b.psi_2nd_xpoint
    xpoint2_x   = b.x[b.null_i[2]]
    xpoint2_y   = b.y[b.null_j[2]]
  ENDELSE
;  LSN = 1
;  USN = 2
;  CDN = 3
;  LDN = 4
;  UDN = 5
;  print, 'i2,xpoint_i1 = ',i2,xpoint1_i1
  print, 'geometry     = ',geometry
;  print, 'psi = ',xpoint1_psi,ctr.psi,xpoint2_psi
  print, 'sep = ',ctr.separatrix

  IF (geometry    EQ param.LDN  AND                     $  ; lower double null
      xpoint2_psi GT ctr.psi    AND                     $  ; outside second separatrix
      xpoint2_x   GT MEAN(ctr.x)) THEN BEGIN            $  ; to the left of the secondary x-point
    IF (section_l[j[1]] EQ 0.0) THEN j = (j[0])[0] ELSE j = (j[1])[0]
    print, ' *** funny business A *** '
  ENDIF ELSE BEGIN IF $
     (geometry    EQ param.LDN  AND                     $  ; lower double null
      xpoint2_psi GT ctr.psi    AND                     $  ; outside second separatrix
      xpoint2_x   LT MEAN(ctr.x)) THEN BEGIN            $  ; to the right of the secondary x-point
    j = (j[0])[0]
    print, ' *** funny business B *** '
  ENDIF ELSE BEGIN IF $
      (geometry    EQ param.UDN  AND                    $  ; upper double null
       xpoint2_psi GT ctr.psi    AND                    $  ; outside second separatrix
       xpoint2_x   LT MEAN(ctr.x)) THEN BEGIN           $  ; to the right of the secondary x-point
    IF (i2 EQ xpoint1_i1[0]) THEN j = (j[0])[0] ELSE j = (j[1])[0]
    print, ' *** funny business C *** '
  ENDIF ELSE BEGIN IF $
      (geometry    EQ param.UDN  AND                    $  ; upper double null
       xpoint2_psi GT ctr.psi    AND                    $  ; outside second separatrix
       xpoint2_x   GT MEAN(ctr.x)) THEN BEGIN           $  ; to the left of the secondary x-point
    j = (j[0])[0]
    print, ' *** funny business D *** '

  ENDIF ELSE BEGIN IF                                       $  ; PRIMARY PFR
    (((geometry EQ param.LSN OR geometry EQ param.LDN) AND  $  ; lower primary null
      xpoint1_psi LT ctr.psi AND                            $  ; private flux region
      xpoint1_y   GT MEAN(ctr.y))                           $  ; below primary x-point
    OR  $
     ((geometry EQ param.USN OR geometry EQ param.UDN) AND  $  ; upper primary null
      xpoint1_psi LT ctr.psi AND                            $  ; private flux region
      xpoint1_y   LT MEAN(ctr.y))) THEN BEGIN               $  ; above primary x-point
    j = (j[1])[0]
    print, ' *** funny business E *** '

  ENDIF ELSE BEGIN IF                                   $  ; SECONDARY PFR
    ((geometry    EQ param.LDN  AND                     $  ; lower primary null
      xpoint2_psi LT ctr.psi    AND                     $  ; private flux region
      xpoint2_y   LT MEAN(ctr.y))                       $  ; above secondary x-point
    OR  $
     (geometry    EQ param.UDN  AND                     $  ; upper primary null
      xpoint2_psi LT ctr.psi    AND                     $  ; private flux region
      xpoint2_y   GT MEAN(ctr.y))) THEN BEGIN           $  ; above primary x-point
    j = (j[1])[0]
    print, ' *** funny business F *** '
    print, ' data = ',xpoint1_psi,ctr.psi,xpoint2_psi
  ENDIF ELSE BEGIN IF                                      $  ; MAIN SOL
   (((geometry EQ param.LSN OR geometry EQ param.USN) AND  $  ; lower/upper single null
     xpoint1_psi GT ctr.psi)                               $  ; 
;     xpoint1_psi LT ctr.psi)                               $  ;  changed 15/11/12 -- not sure why this ever worked...
   OR  $
    ((geometry EQ param.LDN OR geometry EQ param.UDN) AND  $  ; lower/upper double null
      xpoint1_psi GT ctr.psi AND                           $  ; 
      xpoint2_psi LT ctr.psi)) THEN BEGIN                  $  ; 
    j = (j[0])[0]
    print, ' *** funny business MAIN param.SOL *** '
  ENDIF ELSE BEGIN
    print, ' *** funny business, bad news *** '
    print, 'geometry',geometry
    print, 'xpoint1_psi',xpoint1_psi
    print, 'ctr.psi    ',ctr.psi
    print, 'xpoint1_y  ',xpoint1_y   
    print, 'ctr.y      ',MEAN(ctr.y)
    print, 'ctr.sep    ',ctr.separatrix
    print, geometry EQ param.USN
    print, xpoint1_psi LT ctr.psi
    stop
  ENDELSE
  ENDELSE
  ENDELSE
  ENDELSE
  ENDELSE
  ENDELSE
  ENDELSE

END
;
; ======================================================================
;
; ======================================================================
;
FUNCTION grid_CreatePoloidalPoints, c_array, wall, b, debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON commons, param, option, state
  ; ------------------------------------------------------------------

  IF (debug) THEN BEGIN
    PRINT, ' '
    PRINT, '----------------------------------------------------------------------'
    PRINT, ''
    PRINT, 'CREATING POLOIDAL POINTS'
    PRINT, ''
    PRINT, '----------------------------------------------------------------------'
    PRINT, ' '
  ENDIF


  tags  = STRUPCASE(TAG_NAMES(c_array))
  nctr  = N_ELEMENTS(tags)

  IF (KEYWORD_SET(debug)) THEN BEGIN
    PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
    OPLOT,wall.x,wall.y,color=Truecolor('Black')
;    OPLOT,wall.pt1[0,*],wall.pt1[1,*],color=Truecolor('Black'),PSYM=7
    FOR i = 0, nctr-1 DO BEGIN
      ctr = grid_ExtractStructure(c_array,tags[i])      
      OPLOT,ctr.x,ctr.y,color=Truecolor('Red')    
      npts = N_ELEMENTS(ctr.x)
      XYOUTS,ctr.x[npts/2],ctr.y[npts/2],STRTRIM(STRING(i+1),2),color=Truecolor('Darkgrey')
    ENDFOR
  ENDIF

  PRINT,'---------- SETTING UP BASE POLOIDAL POINTS ------------'

  ; Set stuff:  *** PUT IN A SUBROUTINE, STORE IN b STRUCTURE ***

  geometry = state.geometry

  IF (geometry EQ param.LIMITED OR N_ELEMENTS(b.null_i) EQ 2) THEN BEGIN
    xpoint2_psi = -1.0D+6
    xpoint2_x   =  1.0D+6
    xpoint2_y   =  1.0D+6
  ENDIF ELSE BEGIN
    xpoint2_psi = b.psi_2nd_xpoint
    xpoint2_x   = b.x[b.null_i[2]]
    xpoint2_y   = b.y[b.null_j[2]]
  ENDELSE

;
; ----------------------------------------------------------------------
; BUILD A LIST OF SEGMENT POINTS THAT ARE ASSOCIATED WITH X-POINTS
;
  IF (geometry NE param.LIMITED) THEN BEGIN

    xpt_list = [0]
    FOR ictr = 1, nctr DO BEGIN 
      ctr = grid_ExtractStructure(c_array,tags[ictr-1])      
      IF (ctr.separatrix EQ 0) THEN CONTINUE
    
      IF (ctr.separatrix EQ 1) THEN i = 1 ELSE i = 2
    
      focus_x = b.x[b.null_i[i]]
      focus_y = b.y[b.null_j[i]]
    
      j = WHERE( ctr.section_x EQ focus_x AND ctr.section_y EQ focus_y, count)
      IF (count EQ 0) THEN BEGIN
        PRINT,'ERROR grid_CreatePoloidalPoints: x-point not located in section list'
        STOP
      ENDIF
    
      xpt_list = [xpt_list,ctr.section[j]]
    
    ENDFOR
    xpt_list = xpt_list[1:N_ELEMENTS(xpt_list)-1]
    
    print,'xpt_list',xpt_list

  ENDIF
;
; ----------------------------------------------------------------------
; CALCULATE THE LENGTH OF THE CORE POLOIDAL CONTOUR ON THE PRIMARY SEP.
;
  ctr = grid_ExtractStructure(c_array,tags[0])      

  IF (geometry EQ param.LIMITED) THEN BEGIN
    i = [0,N_ELEMENTS(ctr.x)-1]
  ENDIF ELSE BEGIN
    i = WHERE( ctr.x EQ b.x[b.null_i[1]] AND ctr.y EQ b.y[b.null_j[1]], count)
    IF (count NE 2) THEN BEGIN
      PRINT,'ERROR grid_CreatePoloidalPoints: Cannot find primary x-point points (1)'
      STOP
    ENDIF
  ENDELSE
  dummy = grid_GetFrac(ctr.x[i[0]:i[1]],ctr.y[i[0]:i[1]],dist=dist)
  core_boundary = dist[N_ELEMENTS(dist)-1]

  print,'core_boundary',core_boundary
  oplot,ctr.x[i[0]:i[1]],ctr.y[i[0]:i[1]],COLOR=Truecolor('Lightgreen')
;
; ----------------------------------------------------------------------
; 
;

  separatrix_flag = [0,0]

  nsection = 0

  stopper = 0

  FOR ipass = 1, 2 DO BEGIN

;   if (ipass eq 2) THEN stop

    FOR ictr = 1, nctr DO BEGIN 
      ctr = grid_ExtractStructure(c_array,tags[ictr-1])      

      IF ((ipass EQ 1 AND ctr.separatrix EQ 0) OR  $
          (ipass EQ 2 AND ctr.separatrix NE 0)) THEN CONTINUE

      print, '========================================== '
      print, 'ictr ----- ',ictr,' ',tags[ictr-1]
      print, '========================================== '
      print, 'ctr section = ',ctr.section

      FOR i = 0, N_ELEMENTS(ctr.section)-1 DO BEGIN
        OPLOT,ctr.x,ctr.y,color=Truecolor('Cyan')
        XYOUTS,[ctr.section_x[i]], [ctr.section_y[i]], CHARSIZE=2.0,  $
               STRTRIM(STRING(ctr.section[i]),2),color=Truecolor('Black')
      ENDFOR

      dummy = WHERE(STRUPCASE(TAG_NAMES(ctr)) EQ 'SECTION',count)
      IF (count EQ 0) THEN BEGIN
        CONTINUE
      ENDIF
   
      x = ctr.x
      y = ctr.y

      IF (ictr EQ 1 AND geometry EQ param.LIMITED) THEN BEGIN
        ; For the limited grid, the start and end points of the LCFS (and the tangency points) 
        ; are all the same point:
        ctr.section_x[0                        ] = ctr.x[0]
        ctr.section_y[0                        ] = ctr.y[0]
        ctr.section_x[N_ELEMENTS(ctr.section)-1] = ctr.x[0]
        ctr.section_y[N_ELEMENTS(ctr.section)-1] = ctr.y[0]
      ENDIF ELSE  BEGIN
        ; Add in the end points from the wall array:
        FOR j = 1, 2 DO BEGIN
          i = WHERE(wall.ptc EQ ictr AND wall.ptt EQ j)
          p1 = wall.pt1[0:1,i]
          i = grid_LocatePoint(x,y,p1)
          ; If the grid has been processed properly up to this point in the code, 
          ; the following checks should be unnecessary, but just making sure:
          IF ((j EQ 1 AND i NE 0              ) OR  $
              (j EQ 2 AND i NE N_ELEMENTS(x)-2)) THEN BEGIN
            n = N_ELEMENTS(x)
            PRINT,'ERROR grid_CreatePoloidalPoints: Malformed contour (1)'
            PRINT,'  ictr     = ',ictr
            PRINT,'  j        = ',j
            PRINT,'  i        = ',i
            PRINT,'  n        = ',n-1
            PRINT,'  x,y[0]   = ',x[0],y[0]
            PRINT,'  x,y[n-2] = ',x[n-2],y[n-2]
            PRINT,'  p1       = ',p1
            OPLOT,ctr.x,ctr.y,COLOR=Truecolor('Green')
            STOP
          ENDIF
          x = [x[0:i],p1[0],x[i+1:N_ELEMENTS(x)-1]]
          y = [y[0:i],p1[1],y[i+1:N_ELEMENTS(y)-1]]
          IF (j EQ 1) THEN i = 0 ELSE i = N_ELEMENTS(ctr.section)-1
          ctr.section_x[i] = p1[0]
          ctr.section_y[i] = p1[1]
        ENDFOR
        ; Strip the end points of the contour, which are outside the wall (by convention):
        x = x[1:N_ELEMENTS(x)-2]
        y = y[1:N_ELEMENTS(y)-2]
      ENDELSE

      OPLOT,x,y,color=Truecolor('Pink')    
   
      n = N_ELEMENTS(ctr.section)
   
      FOR isec = 0, n-2 DO BEGIN
   
        i1 = ctr.section[isec  ]  
        i2 = ctr.section[isec+1]
   
        x1 = ctr.section_x[isec]
        y1 = ctr.section_y[isec]

print,'x1,y1',isec  ,x1,y1,separatrix_flag
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

print,'x2,y2',isec+1,x2,y2,separatrix_flag

        IF (count GT 1) THEN BEGIN
print,'count >>>>>>', count,k
          IF (ctr.separatrix EQ 1) THEN BEGIN
            IF (state.geometry NE param.LIMITED AND separatrix_flag[1] EQ 0) THEN BEGIN
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

print, '******** MANUAL REFINEMENT! *******'

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
   
;print, 'merdi-licious x',j,new_x
;print, 'merdi-licious y',k,new_y
;print, 'chit i',ctr.section
;print, 'chit x',ctr.section_x
;print, 'chit y',ctr.section_y

;          help,section_t,/struct
   
        ; Scan to look for the current interval in the list:
        IF (ictr EQ 1) THEN BEGIN
          count1 = 0 & i3 = -1 
          count2 = 0 & i4 = -1 
        ENDIF ELSE BEGIN
          i3 = WHERE(section_i1 EQ i1, count1)
          i4 = WHERE(section_i2 EQ i2, count2)

          print,'section_i1 = ',section_i1
          print,'section_i2 = ',section_i2
        ENDELSE

        print,'         ------------ new try -----------'
        print,'i1,3',i1,i3
        print,'i2,4',i2,i4
   
        IF (count1 EQ 0 OR count2 EQ 0) THEN BEGIN

          print,'         ============ path 1 ============'

          ; Process separatrix data, which provides the point distrubution reference to the
          ; rest of the grid:

          ; Setup the poloidal distrubution:

          pol_frac = grid_GetFrac(x,y,dist=dist)
          pol_pos  = grid_PoloidalDistribution(dist[N_ELEMENTS(dist)-1])

;print, 'pol_pos',pol_pos
;print, 'dist',dist
;print,'length', dist[N_ELEMENTS(dist)-1]
;print, 'ictr',ictr


          pol_bnd1 = (pol_frac[j])[0]
          pol_bnd2 = (pol_frac[k])[0]
;
;         --------------------------------------------------------------
;         INCREASE POLOIDAL RESOLUTION NEAR THE SECONDARY X-POINT
;
          IF (1 EQ 1 AND ctr.separatrix EQ 1 AND N_ELEMENTS(b.null_i) GT 2) THEN BEGIN
            xpt_x = b.x[b.null_i[2]]
            xpt_y = b.y[b.null_j[2]]
   print,'2nd x-point',xpt_x,xpt_y

            proximity = SQRT ( (x - xpt_x)^2 + (y - xpt_y)^2 )
            dummy = MIN(proximity,imin)
   print,'x,y',x[imin],y[imin],pol_frac[imin]
            IF (KEYWORD_SET(debug)) THEN BEGIN
              OPLOT,[xpt_x  ],[xpt_y  ],COLOR=Truecolor('Blue'),PSYM=6
              OPLOT,[x[imin]],[y[imin]],COLOR=Truecolor('Blue'),PSYM=6
            ENDIF

            interest  = option.pol_2nd_roi  ; 0.03D  ; parameter
            iteration = 0
            WHILE (1) DO BEGIN
              irange = WHERE( ABS(pol_pos - pol_frac[imin]) LT interest, count)
              IF (count EQ 0) THEN BEGIN
                interest = interest * 2.0D
                IF (interest GE 0.25D) THEN BEGIN
                  ; Something is funny, i.e. no points found to refine, which shouldn't
                  ; happen if the spatial resolution is reasonable:
                  PRINT,'ERROR grid_CreatePoloidalPoints: Issues around 2nd x-point resolution'
                  STOP
                ENDIF
              ENDIF ELSE BEGIN
;print,'interest',interest,N_ELEMENTS(pol_frac)
              print,irange
              print,pol_pos[irange]
; print,'pol1',pol_pos
;plot,pol_frac,psym=3
                FOR l = N_ELEMENTS(irange)-2, 0, -1 DO BEGIN
                  m = irange[l]
                  new_pos = 0.5D * (pol_pos[m] + pol_pos[m+1])
                  pol_pos = [pol_pos[0:m],new_pos,pol_pos[m+1:N_ELEMENTS(pol_pos)-1]]
                ENDFOR
; print,'pol2',pol_pos
;oplot,pol_frac,psym=6
                iteration++
                interest = 0.5D * interest 
                IF (iteration EQ option.pol_2nd_niter) THEN BREAK 
              ENDELSE

            ENDWHILE
;stop
          ENDIF
;
;         --------------------------------------------------------------
;         FORCE SEPARATRIX POINTS AT THE INNER AND OUTER MIDPLANES
;
          IF (ctr.separatrix EQ 1 AND state.geometry NE param.LIMITED) THEN BEGIN

            mid_i = WHERE( x EQ b.x[b.null_i[1]] AND y EQ b.y[b.null_j[1]], count)
            IF (count NE 2) THEN BEGIN
              PRINT,'ERROR grid_CreatePoloidalPoints: Cannot find primary x-point points (2)'
              STOP
            ENDIF

            core_x = x[mid_i[0]:mid_i[1]]
            core_y = y[mid_i[0]:mid_i[1]]
            oplot,core_x,core_y,COLOR=Truecolor('Black')

            dummy = MIN(core_x,imin)
            dummy = MAX(core_x,imax)
            oplot,[core_x[imin]],[core_y[imin]],PSYM=6
            oplot,[core_x[imax]],[core_y[imax]],PSYM=6

;            oplot,[x[imin+i[0]]],[y[imin+i[0]]],PSYM=2
;            oplot,[x[imax+i[0]]],[y[imax+i[0]]],PSYM=2

            ;         inner          outer
            mid_i = [ imin+mid_i[0], imax+mid_i[0] ]

            oplot,[x[mid_i[0]]],[y[mid_i[0]]],PSYM=2
            oplot,[x[mid_i[1]]],[y[mid_i[1]]],PSYM=2

            print,pol_frac[mid_i]

            FOR mid_k = 0, 1 DO BEGIN
              FOR mid_j = 0, N_ELEMENTS(pol_pos)-2 DO  $
                IF (pol_pos[mid_j  ] LT pol_frac[mid_i[mid_k]] AND  $
                    pol_pos[mid_j+1] GT pol_frac[mid_i[mid_k]]) THEN BREAK
              IF (mid_j EQ N_ELEMENTS(pol_pos)-1) THEN BEGIN
                PRINT,'ERROR grid_CreatePoloidalPoints: Cannot find midplane',mid_k
                STOP
              ENDIF
              print,mid_j
              pol_pos = [pol_pos[0:mid_j],pol_frac[mid_i[mid_k]],pol_pos[mid_j+1:N_ELEMENTS(pol_pos)-1]]
            ENDFOR
          ENDIF
;
;         --------------------------------------------------------------
;         SELECT POINTS THAT ARE IN THE CURRENT SEGMENT
;
;          pol_padding = (pol_bnd2 - pol_bnd1) * 0.05D ; PARAMETER
;          pol_padding = 0.0D

        ; in effect briefly for d3d_10 I think (or _9), where the "stop" commands were removed
        ;  print, 'woooooooooooooooooooooooooooooooooooooooo!',j,k
        ;  print,'num',N_ELEMENTS(x)-1
        ;  pol_padding_1 = MIN([0.01D,10.0D*option.pol_res_min])  ; parameter   
        ;  pol_padding_2 = pol_padding_1
        ;
        ;  IF (i1 LT 0) THEN pol_padding_1 = 0.0D ELSE BEGIN $
        ;    dummy = WHERE( i1 EQ xpt_list, count)
        ;    IF (count NE 0) THEN pol_padding_1 = 0.007D * core_boundary  ; parameter
        ;    IF (count NE 0) THEN begin
        ;      print, 'holy shit 1',pol_padding_1
        ;      stop
        ;    endif
        ;  ENDELSE
	;
        ;  IF (i2 LT 0) THEN pol_padding_2 = 0.0D ELSE BEGIN $
        ;    dummy = WHERE( i2 EQ xpt_list, count)
        ;    IF (count NE 0) THEN pol_padding_2 = 0.007D * core_boundary  ; parameter
        ;    IF (count NE 0) THEN begin
        ;      print, 'holy shit 2',pol_padding_2
        ;      stop
        ;    endif
        ;  ENDELSE
        ;
        ; pol_padding_1 = pol_padding_1 / dist[N_ELEMENTS(dist)-1]  ; parameter          
        ; pol_padding_2 = pol_padding_2 / dist[N_ELEMENTS(dist)-1]  ; parameter          

 ;        in effect until 10/01/2012
 ;        and brought back 25/01/2012
 ;        and out again on 02/03/2012 -- kills points near the target, and near x-point problem now solved in _core?
          ; pol_padding_1 = MIN([0.01D,10.0D*option.pol_res_min]) / dist[N_ELEMENTS(dist)-1]  ; parameter
          ; pol_padding_2 = pol_padding_1
 
          IF (i1 LT 0) THEN pol_padding_1 = 0.0D ELSE  $
            pol_padding_1 = MIN([0.01D,10.0D*option.pol_res_min]) / dist[N_ELEMENTS(dist)-1]  ; parameter

          IF (i2 LT 0) THEN pol_padding_2 = 0.0D ELSE  $
            pol_padding_2 = MIN([0.01D,10.0D*option.pol_res_min]) / dist[N_ELEMENTS(dist)-1]  ; parameter

   ;       pol_padding = MIN([0.01D,10.0D*option.pol_res_min]) / dist[N_ELEMENTS(dist)-1]  ; parameter
   ;       pol_padding = 0.1D*MIN([pol_pos[0],1.0D-pol_pos[N_ELEMENTS(pol_pos)-1]]) ; parameter

          l = WHERE(pol_pos GT (pol_bnd1 + pol_padding_1) AND pol_pos LT (pol_bnd2 - pol_padding_2), count)
;          l = WHERE(pol_pos GT (pol_bnd1 + pol_padding) AND pol_pos LT (pol_bnd2 - pol_padding), count)
          IF (count EQ 0) THEN BEGIN
             ; No points available, so add a couple of token points so that the
             ; code below doesn't break, and because you could end up with a tube
             ; that has only 2 cells for extended grids (rare, but possible):
             pol_pos_local = [0.333D,0.667D]
          ENDIF ELSE BEGIN

;          pol_delta = MAKE_ARRAY(N_ELEMENTS(pol_pos)-2,/DOUBLE)
;          FOR l = 1, N_ELEMENTS(pol_pos)-2 DO BEGIN pol_delta[l-1] = pol_pos[l+1] - pol_pos[l]

            pol_pos_local = (pol_pos[l] - pol_bnd1) / (pol_bnd2 - pol_bnd1)

            ; Check that the first and last points are not right next to the 
            ; boundary for this segment:
            IF (N_ELEMENTS(pol_pos_local) GE 2) THEN BEGIN                            
;print,pol_pos_local
;n_save = N_ELEMENTS(pol_pos_local)
              IF (        pol_pos_local[0  ]  LT 0.2D*(pol_pos_local[1  ] - pol_pos_local[0  ])) THEN  $  ; parameter new
                pol_pos_local = pol_pos_local[1:N_ELEMENTS(pol_pos_local)-1]
              n = N_ELEMENTS(pol_pos_local)
              IF (n GE 2) THEN  $
                IF ((1.0D - pol_pos_local[n-1]) LT 0.2D*(pol_pos_local[n-1] - pol_pos_local[n-2])) THEN  $  ; parameter new
                  pol_pos_local = pol_pos_local[0:n-2]
;print,pol_pos_local
;IF (n_save NE N_ELEMENTS(pol_pos_local)) THEN stop
            ENDIF
            ; If there are less than 2 points, assign the minimum required number with 
            ; a nominal (uniform) distribution:
            IF (N_ELEMENTS(pol_pos_local) LT 2) THEN BEGIN
              pol_pos_local = [0.333D,0.667D]
            ENDIF ELSE BEGIN
            ENDELSE
          ENDELSE
;print, pol_pos_local
;if (i1 eq 13) then stop
;
;         --------------------------------------------------------------
;         CHECK FOR ITERATIVE REFINEMENT DUE TO A MALFORMED GRID
;
          IF (state.refinement_pol GT 0) THEN BEGIN
            print, 'the promised land',i1,i2
;            help,state,/struct

;            print,'state pol_i1',state.ref_pol_i1
;            print,'state pol_i2',state.ref_pol_i2


            FOR ref_i = 0, N_ELEMENTS(state.ref_pol_i1)-1 DO BEGIN


;              print,'>>>>>>>>>>>>>>>. wtf?',ref_i,state.ref_pol_i1[ref_i],state.ref_pol_i2[ref_i]

              IF (state.ref_pol_i1[ref_i] EQ i1 AND state.ref_pol_i2[ref_i] EQ i2) THEN BEGIN
;              IF (state.ref_pol_i1[ref_i] EQ i1 AND state.ref_pol_i2[ref_i]) THEN BEGIN  ; BUG! was missing the "EQ i2", so sad - SL, 29/02/2012
                print, 'holy crapoly'


                pol_pos_local = [0.0D,pol_pos_local,1.0D]  ; Add the end points back on for now

print,'pol_pos',pol_pos_local
print,'sta_pos',state.ref_pol_pos[ref_i]

                n = N_ELEMENTS(pol_pos_local)
  
                FOR seg_i = 0, n-2 DO  $
                  IF (pol_pos_local[seg_i  ] LE state.ref_pol_pos[ref_i] AND  $
                      pol_pos_local[seg_i+1] GT state.ref_pol_pos[ref_i]) THEN BREAK

                IF (ABS(pol_pos_local[seg_i  ] - state.ref_pol_pos[ref_i]) LT  $
                    ABS(pol_pos_local[seg_i+1] - state.ref_pol_pos[ref_i])) THEN  $
                  seg_i = [MAX([0,seg_i-1]),     seg_i        ] ELSE  $
                  seg_i = [       seg_i    ,MIN([seg_i+1,n-2])]

print,'seg_i',seg_i

                FOR seg_j = seg_i[1], seg_i[0], -1 DO  $
                  pol_pos_local = [        pol_pos_local[0:seg_j]                          ,  $
                                   0.5D * (pol_pos_local[  seg_j] + pol_pos_local[seg_j+1]),  $
                                           pol_pos_local[seg_j+1:N_ELEMENTS(pol_pos_local)-1]]

print,'pol_pos',pol_pos_local               


                pol_pos_local = pol_pos_local[1:N_ELEMENTS(pol_pos_local)-2]


print,'pol_pos',pol_pos_local               

; stop

              ENDIF



            ENDFOR
;            if (die and state.refinement eq 2) then stop

          ENDIF
;
;         --------------------------------------------------------------
;         COLLECT THE GRID POINTS BASED ON THE POLOIDAL DESCRETISATION
;
          result = grid_SpliceContour(new_x,new_y,position=pol_pos_local)

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
         
          print,'section_i1 = ',section_i1
          print,'section_i2 = ',section_i2
;          help,section_t,/struct
;IF (isec EQ 2) THEN          stop   

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

          ; Seems this check is also required for main param.SOL rings that are inside
          ; the 2nd x-point (all seems a bit arbitrary):

          IF (count1 EQ 2 AND ipass EQ 2 AND ctr.psi GT xpoint2_psi) THEN BEGIN
            i3 = (i3[0])[0]
          ENDIF


          IF (i3 EQ i4) THEN BEGIN

            ; Build a path from the first point to the second, where the reference
            ; points are consecutive in the list:

            print,'         ============ path 2 ============'

            print, 'i3=',i3
            print, 'i4=',i4

            tag = 'data' + STRING(i3,FORMAT='(I0)')

            section = grid_ExtractStructure(section_t,tag)      
;            help,section
;            print,section
            position = section[1:N_ELEMENTS(section)-2]  ; leave off the end points, which are added in grid_SpliceContour

          ENDIF ELSE BEGIN 

            print,'         ============ path 3 ============'

            ; Build a path from the first point to the second, where there are reference
            ; points that are not consecutive in the list:

            j = i3
            k = 0
            length = 0.0

            print,'i3 = ',i3

            WHILE (1) DO BEGIN
              print,'j :',j,k,nsection,format='(A,3I6)'

              IF (k EQ 0) THEN list_j = j ELSE list_j = [list_j,j]

              length = length + section_l[j]

              IF (section_i2[j] EQ i2) THEN BREAK

              print,'section_i1   ',section_i1
              print,'section_i2   ',section_i2
              print,'section_i2[j]',section_i2[j]
;              help,section_i1
;              help,LONG(section_i2[j])

              i2_save = (section_i2[j])[0]

              j = WHERE( section_i1 EQ (section_i2[j])[0], count )

              IF (count EQ 2) THEN BEGIN

;                print, b.psi_1st_xpoint
;                print, b.psi_2nd_xpoint
;                print, ctr.psi
;                print, b.x[b.null_i[0:2]]
;                print, b.y[b.null_j[0:2]]

                print,'j mess',j
                print,'section_i1[j]',section_i1[j]


                grid_ChooseIndex, geometry, b, ctr, section_l, j, i2_save, xpoint1_i1

              ENDIF

              k++
              IF (j EQ -1 OR k EQ nsection) THEN BREAK
            ENDWHILE
    
            print,'j-:',j,k,nsection,format='(A,3I6)'

            IF (count EQ 0) THEN BEGIN
              PRINT,'ERROR grid_CreatePoloidalPoints: What is happening'
              STOP
            ENDIF

            IF (section_i2[j] NE i2) THEN BEGIN
              PRINT,'ERROR grid_CreatePoloidalPoints: Bad situation'
              STOP
            ENDIF
            
            print,'list_j',list_j
            print,'length',length
;            help,section_t,/struct
            shift = 0.0D
            FOR j = 0, N_ELEMENTS(list_j)-1 DO BEGIN

              ; Skip this section if the length is zero, i.e. as is the case when linking
              ; the x-point points for the private flux regions:
              IF (section_l[list_j[j]] EQ 0.0) THEN CONTINUE

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

              IF (j LT N_ELEMENTS(list_j)-1) THEN take = 1 ELSE take = 2

              IF (j EQ 0) THEN position =            frac[0] * section[1:n-take] ELSE  $
                               position = [position, frac[0] * section[1:n-take] + shift]    

              shift = shift + frac[0] 
              print,'position = ',position

              stopper = stopper + 1
            ENDFOR

          ENDELSE

          result = grid_SpliceContour(new_x,new_y,position=position)  
        ENDELSE
   
;        help,grid,/struct
;        print,grid.t
;        oplot,grid.x,grid.y,color=Truecolor('Lightblue'),PSYM=6
;        stop
   
;        oplot,result.x,result.y,color=Truecolor('Lightblue'),PSYM=3


;        if (stopper EQ 3) THEN stop

        IF (isec EQ 0) THEN BEGIN
          grid_x = result.x
          grid_y = result.y
          grid_s = MAKE_ARRAY(N_ELEMENTS(result.x-1),/LONG,VALUE=0L)
          grid_s[0                   ] = i1
          grid_s[N_ELEMENTS(grid_s)-1] = i2
        ENDIF ELSE BEGIN
          nx = N_ELEMENTS(result.x) - 1
          grid_x = [grid_x,result.x[1:nx]] ; drop the first point for all segments except the first
          grid_y = [grid_y,result.y[1:nx]]
          grid_s = [grid_s,MAKE_ARRAY(nx,/LONG,VALUE=0L)]
          grid_s[N_ELEMENTS(grid_s)-1] = i2
        ENDELSE
      
      ENDFOR  ; isec

      oplot,grid_x,grid_y,color=Truecolor('Lightblue'),PSYM=6,SYMSIZE=0.3

;      print,'grid_x',grid_x
;      print,'grid_y',grid_y
;      print,'grid_s',grid_s

      ; Add the grid point data to the contour data structure:
      ctr = CREATE_STRUCT(ctr,'grid_s',grid_s,'grid_x',grid_x,'grid_y',grid_y)
      c_array = grid_UpdateStructure(c_array,tags[ictr-1],ctr)

;      IF (ipass EQ 2 AND ictr EQ 99) THEN BEGIN
;        print, 'stopping as requested',ictr
;        stop
;      ENDIF

;     Need to link the inner and outer PFZ regions by cobbling together
;     the x-point points on the primary separatrix:
      IF (ipass EQ 1 AND ictr EQ 1 AND state.geometry NE param.LIMITED) THEN BEGIN

        print, 'adding primary pfz connector'

        i1 = WHERE(ctr.section_x EQ b.x[b.null_i[1]] AND  $
                   ctr.section_y EQ b.y[b.null_j[1]], count)

        IF (count NE 2) THEN BEGIN
          PRINT, 'ERROR grid_CreatePoloidalPoints: x-points not found'
          PRINT, '  I1 = ',i1
          STOP
        ENDIF

        print, 'b.x[b.null_i[1]]',b.x[b.null_i[1]]
        print, 'b.y[b.null_j[1]]',b.y[b.null_j[1]]
        print, 'ctr.section_x',ctr.section_x
        print, 'ctr.section_y',ctr.section_y
        print, 'ctr.section',ctr.section[i1]
        PRINT, '  I1 = ',i1
  
        nsection++
        tag = 'data' + STRING(nsection-1,FORMAT='(I0)')
        section_i1 = [section_i1,ctr.section[i1[0]]]
        section_i2 = [section_i2,ctr.section[i1[1]]]
        section_l  = [section_l ,0.0]
        section_t  = CREATE_STRUCT(section_t, tag, 0.0)

        xpoint1_i1 = ctr.section[i1[0:1]]
      ENDIF 

;IF (ictr EQ 28) THEN BEGIN
;  print, 'stopping after contour' 
;  stop
;ENDIF
;          IF (state.refinement_pol GT 0 and ictr EQ 1) THEN stop

    ENDFOR  ; ictr
;print, 'stopping after pass 1' 
;stop

  ENDFOR  ; ipass



  state = grid_UpdateStructure(state,'section_i1',section_i1)
  state = grid_UpdateStructure(state,'section_i2',section_i2)

;  print, 'done for now'
;  stop

  result = c_array



  RETURN, result

END


