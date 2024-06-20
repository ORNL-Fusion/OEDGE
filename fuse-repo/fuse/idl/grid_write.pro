;
; ======================================================================
;
; ======================================================================
;
FUNCTION grid_NormalisedFlux, psi, b, mode=mode


;  psi_norm = 1.0D - b.psi_1st_xpoint
;  psi_min = MIN(b.psi) + psi_norm
;  psin = 1.0D / ( (psi + psi_norm - psi_min) / (psi_norm - psi_min) )  ; DIII-D definition (I think)

;  IF (psib EQ 0.0) THEN BEGIN
;    PRINT,'Fising...'
;    psib = 1.0
;    psi = psi + psib
;  ENDIF
;  psi_min = MIN(psi) + psib
;  psin = 1.0 / ((psi+psib-psi_min) / (psib-psi_min))  ; DIII-D definition (I think)
;  print,min(psi+psib-psi_min),max(psi+psib-psi_min)

;  min_psi = MIN(b.psi)
;  psin = 1.0D + psi / min_psi  

  min_psi = MIN(-b.psi)
  psin = 1.0D + psi / min_psi  

  result = psin

  RETURN, result

END
;
; ======================================================================
;
; ======================================================================
;
FUNCTION grid_WriteGridFile, fname, tubes, c_array, wall, b, debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON commons, param, option, state
  ; ------------------------------------------------------------------

; help,b,/struct

  print,' '
  print,'FNAME: ','grid_data/'+fname

  fp = 3
  FREE_LUN, fp
  OPENW, fp, 'grid_data/'+fname, error=err
  IF (err NE 0) THEN BEGIN
    PRINT,'ERROR grid_WriteGridFile: Unable to open file'
    RETURN, -1
  ENDIF

  tags   = STRUPCASE(TAG_NAMES(tubes))
  ntubes = N_ELEMENTS(tags)

  PRINTF,fp,'GRID_OSM'
  PRINTf,fp,'1.00'

  PRINTF,fp,'* ----------------------------------------------------------------------'
  PRINTF,fp,'*'
  PRINTF,fp,'* Grid generation parameters:'
  PRINTF,fp,'*'
  PRINTF,fp,'* boundary[0]     ',option.ctr_boundary[0] ,FORMAT='(A,F10.5)'  
  PRINTF,fp,'* boundary[1]     ',option.ctr_boundary[1] ,FORMAT='(A,F10.5)'  
  PRINTF,fp,'* boundary[2]     ',option.ctr_boundary[2] ,FORMAT='(A,F10.5)'  
  PRINTF,fp,'* boundary[3]     ',option.ctr_boundary[3] ,FORMAT='(A,F10.5)'  
  PRINTF,fp,'* boundary[4]     ',option.ctr_boundary[4] ,FORMAT='(A,F10.5)'  
  PRINTF,fp,'* boundary[5]     ',option.ctr_boundary[5] ,FORMAT='(A,F10.5)'  
  PRINTF,fp,'*'
  PRINTF,fp,'* pol_res_min     ',option.pol_res_min     ,FORMAT='(A,F10.5)'
  PRINTF,fp,'* pol_res_max     ',option.pol_res_max     ,FORMAT='(A,F10.5)'
  PRINTF,fp,'* rad_res_core_min',option.rad_res_core_min,FORMAT='(A,F10.5)'
  PRINTF,fp,'* rad_res_core_max',option.rad_res_core_max,FORMAT='(A,F10.5)'
  PRINTF,fp,'* rad_res_sol_min ',option.rad_res_sol_min ,FORMAT='(A,F10.5)'
  PRINTF,fp,'* rad_res_sol_max ',option.rad_res_sol_max ,FORMAT='(A,F10.5)'
  PRINTF,fp,'*'
  PRINTF,fp,'* ----------------------------------------------------------------------'
  PRINTF,fp,'{NUMBER OF TUBES} ',ntubes
  PRINTF,fp,'* ----------------------------------------------------------------------'
  PRINTF,fp,'{TUBE DATA}'
  PRINTF,fp,'* headers'
  FOR itube = 1, ntubes DO BEGIN
    tube = grid_ExtractStructure(tubes,tags[itube-1])      
    psin1 = grid_NormalisedFlux(tube.psi1,b)
    psin2 = grid_NormalisedFlux(tube.psi2,b)
    PRINTF,fp,itube,tube.n,tube.region,tube.separatrix1,tube.separatrix2,  $
              tube.psi1,tube.psi2, psin1, psin2,  $
              tube.ictr1,tube.ictr2,tube.map_in,tube.map_out,  $
              FORMAT='(5I6,4F15.10,6I6)'
  ENDFOR

  note_reg = ['CORE','SOL','SOL','SOL','PFZ','PFZ']
  note_sep = [' ',', PRIMARY SEPARATRIX',', SECONDARY SEPARATRIX',', SECONDARY SEPARATRIX (BRANCH 1)',', SECONDARY SEPARATRIX (BRANCH 2)']

  nx = N_ELEMENTS(b.x)
  ny = N_ELEMENTS(b.y)

  PRINTF,fp,'* ----------------------------------------------------------------------'
  PRINTF,fp,'{CELL DATA}'
  FOR itube = 1, ntubes DO BEGIN

    tube = grid_ExtractStructure(tubes,tags[itube-1])      
    PRINTF,fp,'* ',note_reg[tube.region-1],note_sep[tube.separatrix1]
    PRINTF,fp,'* headers'

    FOR icell = 0, tube.n-1 DO BEGIN

      x1 = 0.5D * (tube.x2[icell  ] + tube.x1[icell  ])
      y1 = 0.5D * (tube.y2[icell  ] + tube.y1[icell  ])
      x2 = 0.5D * (tube.x2[icell+1] + tube.x1[icell+1])
      y2 = 0.5D * (tube.y2[icell+1] + tube.y1[icell+1])

      x1 = (x1 - MIN(b.x)) / (MAX(b.x)-MIN(b.x)) * (nx-1)
      y1 = (y1 - MIN(b.y)) / (MAX(b.y)-MIN(b.y)) * (ny-1)
      x2 = (x2 - MIN(b.x)) / (MAX(b.x)-MIN(b.x)) * (nx-1)
      y2 = (y2 - MIN(b.y)) / (MAX(b.y)-MIN(b.y)) * (ny-1)

      bratio = INTERPOLATE(b.b_ratio,[x1,x2],[y1,y2]) ; ,CUBIC=-0.5D)

      bratio = MEAN(bratio)

      PRINTF,fp,icell+1,itube,bratio,  $
                tube.x2[icell+1],tube.y2[icell+1],tube.x1[icell+1],tube.y1[icell+1],  $
                FORMAT='(2I6,F12.7,2X,4F15.10)'
      PRINTF,fp,tube.x2[icell  ],tube.y2[icell  ],tube.x1[icell  ],tube.y1[icell  ],  $
                FORMAT='(24X      ,2X,4F15.10)'
    ENDFOR
  ENDFOR

  PRINTF,fp,'* ----------------------------------------------------------------------'
  PRINTF,fp,'{WALL DATA}'

  PRINTF,fp,N_ELEMENTS(wall.pti),FORMAT='(I6)'

  FOR iwall = 0, N_ELEMENTS(wall.pti)-1 DO BEGIN
    PRINTF,fp,iwall,  $
              wall.ptc[iwall],wall.ptt[iwall],  $
              wall.pt1[0:1,iwall],wall.pt2[0:1,iwall],  $
              FORMAT='(I6,2I8,2X,2F15.10,2X,2F15.10)'

  ENDFOR

  PRINTF,fp,'* ----------------------------------------------------------------------'
  PRINTF,fp,'{END}'
  PRINTF,fp,' '

  CLOSE,fp

END

