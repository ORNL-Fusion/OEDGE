;
; ======================================================================
;
; ======================================================================
;
FUNCTION grid_WriteGridFile, fname, tubes, c_array, wall, b, debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON params, CORE, SOL, PFZ, FORWARD, BACKWARD, LOWER_NULL, UPPER_NULL, HFS, LFS 
  ; ------------------------------------------------------------------

  print,' '
  print,'FNAME: ',fname

  fp = 3
  FREE_LUN, fp
  OPENW, fp, fname, error=err
  IF (err NE 0) THEN BEGIN
    PRINT,'ERROR grid_WriteGridFile: Unable to open file'
    RETURN, -1
  ENDIF

  tags   = STRUPCASE(TAG_NAMES(tubes))
  ntubes = N_ELEMENTS(tags)

  PRINTF,fp,'GRID_OSM'
  PRINTf,fp,'1.00'

  PRINTF,fp,'* ----------------------------------------------------------------------'
  PRINTF,fp,'{NUMBER OF TUBES} ',ntubes
  PRINTF,fp,'* ----------------------------------------------------------------------'
  PRINTF,fp,'{TUBE DATA}'
  PRINTF,fp,'* headers'
  FOR itube = 1, ntubes DO BEGIN
    tube = grid_ExtractStructure(tubes,tags[itube-1])      
    PRINTF,fp,itube,tube.n,tube.region,tube.separatrix1,tube.separatrix2,  $
              tube.psi1,tube.psi2,  $
              tube.ictr1,tube.ictr2,tube.map_in,tube.map_out,  $
              FORMAT='(5I6,2F15.10,6I6)'
  ENDFOR

  note_reg = ['none','CORE','SOL','none','none','PFZ']
  note_sep = [' ',', PRIMARY SEPARATRIX',', SECONDARY SEPARATRIX',', SECONDARY SEPARATRIX (BRANCH 1)',', SECONDARY SEPARATRIX (BRANCH 2)']

  PRINTF,fp,'* ----------------------------------------------------------------------'
  PRINTF,fp,'{CELL DATA}'
  FOR itube = 1, ntubes DO BEGIN
    tube = grid_ExtractStructure(tubes,tags[itube-1])      
    PRINTF,fp,'* ',note_reg[tube.region],note_sep[tube.separatrix1]
    PRINTF,fp,'* headers'
    FOR icell = 0, tube.n-1 DO BEGIN
      PRINTF,fp,icell+1,itube,-1.0D,  $
                tube.x2[icell+1],tube.y2[icell+1],tube.x1[icell+1],tube.y1[icell+1],  $
                FORMAT='(2I6,F12.7,2X,4F15.10)'
      PRINTF,fp,tube.x2[icell  ],tube.y2[icell  ],tube.x1[icell  ],tube.y1[icell  ],  $
                FORMAT='(24X      ,2X,4F15.10)'
    ENDFOR
  ENDFOR

  PRINTF,fp,'* ----------------------------------------------------------------------'
  PRINTF,fp,'{END}'
  PRINTF,fp,' '


  CLOSE,fp

END

