;
; ======================================================================
;
FUNCTION cortex_LoadWalldyn, plot

  PRINT,'Loading WALLDYN data'

  CASE 0 OF
    ; ----------------------------------------------------------------------
    0: BEGIN
      ; --------------------------------------------------    
      buffer = ' '
      idum = MAKE_ARRAY(10000,1000,/LONG  ,VALUE=-1L )
      rdum = MAKE_ARRAY(10000,1000,/FLOAT ,VALUE=-1.0)
      sdum = MAKE_ARRAY(10000     ,/STRING,VALUE='null')
      ; --------------------------------------------------    
      fp1 = 3
      FREE_LUN,fp1        
      file_name = 'cortex_walldyn/0903_wall.dat'
      OPENR,fp1,file_name, ERROR=err
      IF (err NE 0) THEN BEGIN
        PRINT,'ERROR cortex_LoadWalldyn: Unable to open ',file_name
        PRINT,'  FILE_NAME= ',file_name
        PRINTF,-2,!err.msg
        FREE_LUN,fp1
        STOP
      ENDIF   
      ; --------------------------------------------------    
      READF,fp1,buffer
      i = -1
      WHILE NOT EOF(fp1) DO BEGIN  
        i++
        READF, fp1, buffer
        str = STRSPLIT(buffer,/EXTRACT)
        idum[i,0:1] = LONG (str[0:1])  ; I can see to assign the variables directly in the READF statement... strange...
        rdum[i,0:4] = FLOAT(str[2:6])
      ENDWHILE        
      dist = MAKE_ARRAY(i+1,/FLOAT,VALUE=0.0)
      dist[0] = 0.5 * rdum[0,4]
      FOR j = 1, i DO dist[j] = dist[j-1] + 0.5 * (rdum[j-1,4] + rdum[j,4])
      wall_0903 = {  $
        version   : 1.0         ,  $
        file      : file_name   ,  $
        n         : i+1         ,  $ 
        bin_start : idum[0:i,0] ,  $
        bin_end   : idum[0:i,1] ,  $
        R_start   : rdum[0:i,0] ,  $
        Z_start   : rdum[0:i,1] ,  $
        R_end     : rdum[0:i,2] ,  $
        Z_end     : rdum[0:i,3] ,  $
        length    : rdum[0:i,4] ,  $
        dist      : dist        }
      CLOSE,fp1
      ; --------------------------------------------------
      file_name = 'cortex_walldyn/0904_wall.dat'
      OPENR,fp1,file_name, ERROR=err
      READF,fp1,buffer
      i = -1
      WHILE NOT EOF(fp1) DO BEGIN  
        i++
        READF, fp1, buffer
        str = STRSPLIT(buffer,/EXTRACT)
        idum[i,0:1] = LONG (str[0:1])
        rdum[i,0:4] = FLOAT(str[2:6])
      ENDWHILE  
      dist = MAKE_ARRAY(i+1,/FLOAT,VALUE=0.0)
      dist[0] = 0.5 * rdum[0,4]
      FOR j = 1, i DO dist[j] = dist[j-1] + 0.5 * (rdum[j-1,4] + rdum[j,4])
      wall_0904 = {  $
        version   : 1.0         ,  $
        file      : file_name   ,  $
        n         : i+1         ,  $ 
        bin_start : idum(0:i,0) ,  $
        bin_end   : idum(0:i,1) ,  $
        R_start   : rdum(0:i,0) ,  $
        Z_start   : rdum(0:i,1) ,  $
        R_end     : rdum(0:i,2) ,  $
        Z_end     : rdum(0:i,3) ,  $
        length    : rdum(0:i,4) ,  $
        dist      : dist        }
      CLOSE,fp1
      ; --------------------------------------------------
      file_name = 'cortex_walldyn/NetBeDep.dat'
      OPENR,fp1,file_name, ERROR=err
      READF,fp1,buffer
      i = -1
      WHILE NOT EOF(fp1) DO BEGIN  
        i++
        READF, fp1, buffer
        str = STRSPLIT(buffer,/EXTRACT)
        j = N_ELEMENTS(str)
        sdum(i      ) =       str(0    )
        rdum(i,0:j-2) = FLOAT(str(1:j-1)) 
      ENDWHILE  
      netbedep  = {  $
        version : 1.0             ,  $
        file    : file_name       ,  $
        n       : j-1             ,  $ 
        run     : sdum(0:i)       ,  $
        data    : rdum(0:i,0:j-2) }
      CLOSE,fp1
      ; --------------------------------------------------
      file_name = 'cortex_walldyn/codepdata.dat'
      OPENR,fp1,file_name, ERROR=err
      READF,fp1,buffer
      i = -1
      WHILE NOT EOF(fp1) DO BEGIN  
        i++
        READF, fp1, buffer
        str = STRSPLIT(buffer,/EXTRACT)
        j = N_ELEMENTS(str)
        sdum(i      ) =       str(0    )
        rdum(i,0:j-2) = FLOAT(str(1:j-1)) 
      ENDWHILE  
      codep  = {  $
        version : 1.0             ,  $
        file    : file_name       ,  $
        n       : j-1             ,  $ 
        run     : sdum(0:i)       ,  $
        data    : rdum(0:i,0:j-2) }
      CLOSE,fp1
      ; --------------------------------------------------
      file_name = 'cortex_walldyn/walltempdata.dat'
      OPENR,fp1,file_name, ERROR=err
      READF,fp1,buffer
      i = -1
      WHILE NOT EOF(fp1) DO BEGIN  
        i++
        READF, fp1, buffer
        str = STRSPLIT(buffer,/EXTRACT)
        j = N_ELEMENTS(str)
        sdum(i      ) =       str(0    )
        rdum(i,0:j-2) = FLOAT(str(1:j-1)) 
      ENDWHILE  
      walltemp  = {  $
        version : 1.0             ,  $
        file    : file_name       ,  $
        n       : j-1             ,  $ 
        run     : sdum(0:i)       ,  $
        data    : rdum(0:i,0:j-2) }
      CLOSE,fp1
      ; --------------------------------------------------
      file_name = 'cortex_walldyn/i-dib-0903-1514-01g_ipp_Be_fluxsum.dat'
      OPENR,fp1,file_name, ERROR=err
      READF,fp1,buffer
      str = STRSPLIT(buffer,/EXTRACT)
      n = N_ELEMENTS(str)
      FOR i = 0, n-2 DO BEGIN
        READF, fp1, buffer
        str = STRSPLIT(buffer,/EXTRACT)
        rdum[0:n-2,i] = FLOAT(str[1:n-1]) 
      ENDFOR
      data1 = {  $
        version : 1.0               ,  $
        file    : file_name         ,  $
        n       : n-1               ,  $ 
        data    : rdum[0:n-2,0:n-2] }
      CLOSE,fp1
      file_name = 'cortex_walldyn/i-dib-0903-1514-00m_ipp_Be_fluxsum.dat'
      OPENR,fp1,file_name, ERROR=err
      READF,fp1,buffer
      str = STRSPLIT(buffer,/EXTRACT)
      n = N_ELEMENTS(str)
      FOR i = 0, n-2 DO BEGIN
        READF, fp1, buffer
        str = STRSPLIT(buffer,/EXTRACT)
        rdum[0:n-2,i] = FLOAT(str[1:n-1]) 
      ENDFOR
      data2 = {  $
        version : 1.0               ,  $
        file    : file_name         ,  $
        n       : n-1               ,  $ 
        data    : rdum[0:n-2,0:n-2] }
      CLOSE,fp1
      redep_matrix = { data1 : data1, data2 : data2 }
      ; --------------------------------------------------
      END
    ; ----------------------------------------------------------------------
    ELSE: BEGIN  
      PRINT,'big bad wolf'
      STOP
      END
    ; ----------------------------------------------------------------------
  ENDCASE

  result = {  $
    wall_0903    : wall_0903    ,  $  
    wall_0904    : wall_0904    ,  $  
    netbedep     : netbedep     ,  $
    codep        : codep        ,  $
    walltemp     : walltemp     ,  $
    redep_matrix : redep_matrix }

  RETURN,result
END
