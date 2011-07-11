;
; ======================================================================
;
PRO grid_SetupParameters
  ; ------------------------------------------------------------------
  COMMON params, CORE, SOL, PFZ, FORWARD, BACKWARD, LOWER_NULL, UPPER_NULL, HFS, LFS
                 
  FORWARD       =  1.0D
  BACKWARD      = -1.0D
  CORE          = 1
  SOL           = 2
  SOL_LFS       = 3  
  SOL_HFS       = 4  
  PFZ           = 5  
  PFZ_SECONDARY = 6  

  HFS = 1
  LFS = 2

  LOWER_NULL    = 1
  UPPER_NULL    = 2
  ; ------------------------------------------------------------------
END
