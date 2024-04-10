;
; ======================================================================
;
FUNCTION inExtractStructure, struct, tag

  index = WHERE(TAG_NAMES(struct) EQ STRUPCASE(tag), count)

  IF (count NE 0) THEN BEGIN
    result = struct.(index[0]) 
  ENDIF ELSE BEGIN
    PRINT, 'ERROR grid_ExtractStructure: TAG not found'
    PRINT, '  TAG = >'+STRTRIM(tag,2)+'<'
    PRINT, '  NAMES = ',TAG_NAMES(struct)
;    HELP, struct, /struct
;    HELP, struct.contour12,/struct
    result = -1
  ENDELSE

  RETURN, result
END

;
; ======================================================================
;
; inGetDataStructure
; 
; A function that returns an empty data structure, to be filled when
; reading in the data file.  The type and size of the array is set from
; information in the data file header.
;
; DATA_TYPE	INTEGER	  A code that specifies whether or not the data is
;                         INTEGER*4 / LONG, REAL*4 / FLOAT, or REAL*8 /
;                         DOUBLE.
; N		INTEGER   Specifies the number of elements in the data
;                         array.
;
FUNCTION inGetDataStructure, data_type, n

  header = {                                        $
    header_version : 0.1                         ,  $  ;
    file           : STRING(' ',FORMAT='(A256)') ,  $  ;  
    tag            : STRING(' ',FORMAT='(A256)') ,  $  ;  
    units          : STRING(' ',FORMAT='(A256)') ,  $  ;
    type           : data_type                   ,  $  ;
    n              : n}                                ;

  CASE data_type OF  
;    1: 
    2: struct = CREATE_STRUCT(header,'data',LONARR(n))  ; INTEGER*4
    3: struct = CREATE_STRUCT(header,'data',FLTARR(n))  ; REAL*4
    4: struct = CREATE_STRUCT(header,'data',DBLARR(n))  ; REAL*8
    5: struct = CREATE_STRUCT(header)                   ; STRING
    ELSE: BEGIN
      PRINT,'ERROR inGetDataStructure: Unknown data type'
      STOP
      END
  ENDCASE

  RETURN, struct
END
;
; ========================================================================
;
; FUNCTION inReadLine
;
; Reads in a single line from the data file.  Comment characters are 
; processed, i.e. everything on a line that is after a "$" or "*" is
; ignored.  If an entire line is ignored (the comment character is at the
; very beginning or there are only spaces) then another line is read, 
; until a line with something on it is found or the end of the file is
; reached.
;
; FP		INTEGER  The file pointer, which was set when the file 
;                        was opened
; BUFFER 	STRING   The contents of the data line from the file.
;
FUNCTION inReadLine, fp, buffer

  WHILE (1) DO BEGIN

    READF,fp,buffer

    i = STRPOS(buffer,'$')
    IF (i EQ  0) THEN CONTINUE 
    IF (i NE -1) THEN BEGIN
      str    = STRSPLIT(buffer,'$',/EXTRACT)
      buffer = str[0]
    ENDIF

    i = STRPOS(buffer,'*')   ; Combine with above check somehow...
    IF (i EQ  0) THEN CONTINUE 
    IF (i NE -1) THEN BEGIN
      str    = STRSPLIT(buffer,'*',/EXTRACT)
      buffer = str[0]
    ENDIF

    IF (STRMATCH(buffer,'*{FILE END}*') EQ 1) THEN RETURN, 0
    IF (STRMATCH(buffer, '*{*'        ) EQ 1) THEN RETURN, 1
    RETURN, 2

  ENDWHILE
END 

;
; ========================================================================
;
; FUNCTION inLoadData
;
; Reads in all the data in an INTERFACE data file.
;
PRO inLoadData, file_name, debug=debug, no_check=no_check

  COMMON in_commons, param, master


;  print,'loading data for ', file_name

; Initializations:
  fp          = 2
  buffer      = ' '


; Check if the data is already loaded:
  IF (N_ELEMENTS(master) NE 0 AND NOT KEYWORD_SET(no_check)) THEN BEGIN
    n = N_ELEMENTS(TAG_NAMES(master))
    FOR i = 0, n-1 DO BEGIN
      tag = 'data' + STRING(i+1,FORMAT='(I0)')
      struct = inExtractStructure(master,tag)
      IF (struct.file EQ file_name) THEN BEGIN
;        PRINT,'MESSAGE inLoadData: Data file already loaded!'
        POINT_LUN, fp, 0
        RETURN
      ENDIF  
    ENDFOR
  ENDIF


  cont = 1  
  WHILE (cont GE 1) DO BEGIN
;   Read a line from the data file:
    IF (cont EQ 1) THEN status = inReadLine(fp, buffer)
    cont = 1

;   Extract tag and isolate data:
    i1   = STRPOS(buffer,'{')           
    i2   = STRPOS(buffer,'}')           
    tag  = STRMID(buffer,i1+1,i2-i1-1)  
    data = STRMID(buffer,i2+1)          

    IF (KEYWORD_SET(debug)) THEN PRINT,'TAG:',tag,'>'+data+'<'

    CASE tag OF
;     Collect information about how the data is arranged in the file: 
      'FILE INDENT'      : file_indent  = LONG(data)
      'FILE INDEX'       : file_index   = LONG(data)
      'FILE COLUMNS'     : file_columns = LONG(data)
      'FILE COLUMN WIDTH': column_width = LONG(STRSPLIT(data,/EXTRACT))
      ; a (bit of a) hack for loading string data
      'STRINGS'          : BEGIN
         status = inReadLine(fp,buffer)
         n = LONG(buffer)
         FOR i = 0, n-1 DO BEGIN
           status = inReadLine(fp,buffer)            
           IF (i EQ 0) THEN string_array = buffer ELSE string_array = [string_array,buffer]
         ENDFOR

         struct = inGetDataStructure(5,n)
         struct.file = file_name
         struct.tag  = 'STRINGS'
         struct = CREATE_STRUCT(struct,'data',string_array)
         name = 'data1'
         struct_array = CREATE_STRUCT('data1',struct)

         END
;     Data related tags:
      'DATA TAG': BEGIN
;        Search for a match to the requested data tag in the current 
;        set of data columns:
         FOR i = 0, file_columns-1 DO BEGIN
           i1   = STRPOS(data,'{')
           i2   = STRPOS(data,'}')
           tag  = STRMID(data,i1+1,i2-i1-1)
           data = STRMID(data,i2+1)

           IF (i EQ 0) THEN data_tag = tag ELSE data_tag = [data_tag,tag] 

         ENDFOR
         END
      'DATA TYPE'  : data_type = LONG(STRSPLIT(data,/EXTRACT))
      'DATA N'     : data_n    = LONG(STRSPLIT(data,/EXTRACT))
      'DATA VALUES': BEGIN
;        Read and assign the numerical data stored in the file:

         ; Build the structure array that will store the data:
         FOR i = 0, file_columns-1 DO BEGIN
           struct = inGetDataStructure(data_type[i],data_n[i])
           struct.file = file_name
           struct.tag  = data_tag[i]
           name = 'data' + STRING(i+1,FORMAT='(I0)')
           IF (i EQ 0) THEN struct_array = CREATE_STRUCT(             name,struct) ELSE  $
                            struct_array = CREATE_STRUCT(struct_array,name,struct)
         ENDFOR


         n = -1L
         WHILE (inReadLine(fp,buffer) EQ 2) DO BEGIN

           n = n + 1L

           i1 = file_indent

           FOR i = 0, file_columns-1 DO BEGIN

             data = STRMID(buffer,i1,column_width[i])

             IF (KEYWORD_SET(debug)) THEN PRINT,i1,i2,'>'+data+'<',n

             IF (n LT struct_array.(i).n) THEN BEGIN

               CASE struct_array.(i).type OF
                 ; 1:
                 2: struct_array.(i).data[n] = LONG  (data)
                 3: struct_array.(i).data[n] = FLOAT (data)
                 4: struct_array.(i).data[n] = DOUBLE(data)
                 ELSE: BEGIN
                   PRINT,'ERROR inLoadData: Unknown data type'
                   STOP
                   END
               ENDCASE

             ENDIF

             i1 = i1 + column_width[i]

           ENDFOR

         ENDWHILE

         cont =  2
         END
      ELSE:  
;       Keep looking through the file:
    ENDCASE

    IF (N_ELEMENTS(struct_array) NE 0) THEN BEGIN

      IF (N_ELEMENTS(master) EQ 0) THEN  $
        master = struct_array  $
      ELSE BEGIN
;        IF (KEYWORD_SET(test)) THEN BEGIN
;          i = 10
;          name = 'data' + STRING(i,FORMAT='(I0)')
;          master = struct_array  
;        ENDIF ELSE BEGIN
          n1 = N_ELEMENTS(TAG_NAMES(master))
          n2 = N_ELEMENTS(TAG_NAMES(struct_array))
          FOR i = n1+1, n1+n2 DO BEGIN
            name = 'data' + STRING(i,FORMAT='(I0)')
            master = CREATE_STRUCT(master,name,struct_array.(i-(n1+1)))
          ENDFOR
;        ENDELSE
      ENDELSE

      UNDEFINE, struct_array

    ENDIF

;   End of file marker:
    IF (STRMATCH(buffer,'*{FILE END}*') EQ 1) THEN BREAK
  ENDWHILE

; Rewind the data file for the next read-through:
  POINT_LUN, fp, 0

  RETURN
END
;
; ========================================================================
;
; FUNCTION inOpenInterface
;
; Opens the data file.
;
FUNCTION inOpenInterface, file_name, debug=debug, no_check=no_check

  COMMON in_commons, param, master

  param = { file_name : file_name }

;  PRINT, '***************************'
;  PRINT, '******** OPENING! *********',file_name
;  PRINT, '***************************'

; Check that FILE_NAME is valid:
;   STILL TO DO

  IF (N_ELEMENTS(master) GT 0) THEN  $   ; limit how much data is stored (without this, large numbers of datasets can make things very slow) -- SL, 13/06/2017
    IF (N_ELEMENTS(TAG_NAMES(master)) GT 100) THEN UNDEFINE, master

; Open the data stream:
  fp = 2
  FREE_LUN,fp

  IF (KEYWORD_SET(debug)) THEN PRINT, 'Opening interface for ',file_name

  OPENR,fp,file_name,error=error
  IF (error NE 0) THEN BEGIN
    PRINT,'ERROR inOpenInterface: Unable to access data file'
    PRINT,' FILE= ',file_name
    RETURN, -1
  ENDIF  

  READF,fp,version

  inLoadData, file_name, debug=debug, no_check=no_check

; Perhaps need to check the file integrity, ie are there duplicate tags...
;   STILL TO DO

;   print,'master', N_ELEMENTS(TAG_NAMES(master)), N_ELEMENTS(master)

  RETURN, 0

END
;
; ========================================================================
;
PRO inCloseInterface
  fp = 2
  CLOSE, fp
END
;
; ========================================================================
; 
; FUNCTION inGetData
;
; Main routine for accessing the formatted CORTEX data file.  A tag is 
; supplied and the file that's currently open for access (from a previous
; call to inOpenInterface) is scanned to see if the tag exists. If yes,
; then great.
;
; DATA_TAG	STRING  Uh, the tag used to identify the data.
; FULL		        Send back a data structure that includes the meta-
;                       data, like the data type and the file where the
;                       data was loaded from.  If FULL is not specified
;                       then only the data itself is returned, i.e. a 
;                       one dimensional array of the appropriate type.
; FLUSH			Clear the data from previous calls to inGetData
;                       that's currently begin held in memory.  NOT 
;                       IMPLIMENTED YET.
;
FUNCTION inGetData,data_tag,full=full,flush=flush,debug=debug

  COMMON in_commons, param, master

; Check that the input parameters are valid:
;   STILL TO DO

; Initializations:
  fp          = 2
  buffer      = ' '
  data_column = -1
  data_found  = 0


; Search to see if the requested data is already in memory:
  IF (N_ELEMENTS(master) NE 0) THEN BEGIN

;help,master,/struct

;      print,'master ',data_tag

    n = N_ELEMENTS(TAG_NAMES(master))
    FOR i = 0, n-1 DO BEGIN

; print,'--------------------',i,n
;help,master.(i),/struct
;      print,'master ',data_tag,'   ',master.(i).tag

      IF (master.(i).file EQ param.file_name AND master.(i).tag EQ data_tag) THEN BEGIN
        data_found = 1
        struct = master.(i)
        BREAK        
      ENDIF  

    ENDFOR

  ENDIF

; Nothing found so access the data file:
  IF (NOT data_found) THEN BEGIN


    PRINT,'SHOULD NOT BE HERE!'
    print,param.file_name
    print,data_tag
    STOP

    cont = 1  
    WHILE (cont GE 1) DO BEGIN
;     Read a line from the data file:
      IF (cont EQ 1) THEN status = inReadLine(fp, buffer)
      cont = 1

;     Extract tag and isolate data:
      i1   = STRPOS(buffer,'{')           
      i2   = STRPOS(buffer,'}')           
      tag  = STRMID(buffer,i1+1,i2-i1-1)  
      data = STRMID(buffer,i2+1)          

      IF (KEYWORD_SET(debug)) THEN PRINT,'TAG:',tag,'>'+data+'<'

      CASE tag OF
;       Collect information about how the data is arranged in the file: 
        'FILE INDENT'      : file_indent  = LONG(data)
        'FILE INDEX'       : file_index   = LONG(data)
        'FILE COLUMNS'     : file_columns = LONG(data)
        'FILE COLUMN WIDTH': column_width = LONG(STRSPLIT(data,/EXTRACT))
;       Data related tags:
        'DATA TAG': BEGIN
;          Search for a match to the requested data tag in the current 
;          set of data columns:
           FOR i = 0, file_columns-1 DO BEGIN
             i1   = STRPOS(data,'{')
             i2   = STRPOS(data,'}')
             tag  = STRMID(data,i1+1,i2-i1-1)
             data = STRMID(data,i2+1)
             IF (STRMATCH(tag,data_tag)) THEN data_column = i
           ENDFOR
           END
        'DATA TYPE'  : data_type = LONG(STRSPLIT(data,/EXTRACT))
        'DATA N'     : data_n    = LONG(STRSPLIT(data,/EXTRACT))
        'DATA VALUES': BEGIN
;          Read and assign the numerical data stored in the file:
           IF (data_column NE -1) THEN BEGIN
             struct = inGetDataStructure(data_type[data_column],data_n[data_column])
             struct.tag = data_tag
           ENDIF
           n = -1L
           WHILE (inReadLine(fp,buffer) EQ 2) DO BEGIN
             IF (data_column EQ -1) THEN CONTINUE
             n = n + 1
             IF (n LT struct.n) THEN BEGIN
               i1 = file_indent
               IF (data_column GT 0) THEN i1 = i1 + LONG(TOTAL(column_width[0:data_column-1]))
               data = STRMID(buffer,i1,column_width[data_column])
               IF (KEYWORD_SET(debug)) THEN PRINT,i1,i2,'>'+data+'<',n
               CASE struct.type OF
                 ; 1:
                 2: struct.data[n] = LONG  (data)
                 3: struct.data[n] = FLOAT (data)
                 4: struct.data[n] = DOUBLE(data)
                 ELSE: BEGIN
                   PRINT,'ERROR inGetData: Unknown data type'
                   STOP
                   END
               ENDCASE
             ENDIF
           ENDWHILE
           data_column = -1
           cont        =  2
           END
        ELSE:  
;         Keep looking through the file:
      ENDCASE
;     End of file marker:
      IF (STRMATCH(buffer,'*{FILE END}*') EQ 1) THEN BREAK
    ENDWHILE
;   Rewind the data file for the next read-through:
    POINT_LUN, fp, 0
  ENDIF

;  IF (NOT KEYWORD_SET(struct)) THEN RETURN, -1
  IF (NOT KEYWORD_SET(struct)) THEN BEGIN
    PRINT,'ERROR inGetData: Data tag ',data_tag,' not found'
    PRINT,'  FILE_NAME = ',param.file_name
    STOP
  ENDIF

  ; Pass something back:
  IF (KEYWORD_SET(full)           ) THEN RETURN, struct            ; the entire structure, with meta-data
  IF (N_ELEMENTS(struct.data) EQ 1) THEN RETURN, (struct.data)[0]  ; just a single value, if that's all there is
                                         RETURN, struct.data       ; the data array
END
;
; ======================================================================
;
PRO inTest,tube=tube

  path = '/home/slisgo/divimp/results/'

  file_name = path + 'u-lin-0003d.test.idl'
;  file_name = path + 'm-hmr-0012b.test.idl'

  inOpenInterface, file_name
  b_tube = inGetData('TUBE',/full)
  b_s    = inGetData('s')
  b_Te   = inGetData('Te')
  b_ne   = inGetData('ne')
  inCloseInterface

;  file_name = path + 'u-lin-0002i.test.idl'
;  inOpenInterface, file_name
;  c_tube = inGetData('TUBE',/full)
;  c_s    = inGetData('s')
;  c_Te   = inGetData('Te')
;  c_ne   = inGetData('ne')
;  inCloseInterface

  
;  b_s = REVERSE(b_s)
;  c_s = REVERSE(c_s)

  safe_colors,/first
  !P.MULTI = [0, 1, 2] 
  IF (NOT KEYWORD_SET(tube)) THEN tube = 1
  b_i = WHERE(b_tube.data EQ tube)  
;  c_i = WHERE(c_tube.data EQ tube)  

  print,b_s[b_i]
  print,b_Te[b_i]

  plot ,b_s[b_i],b_Te[b_i],yrange=[MIN(b_Te[b_i]),MAX(b_Te[b_i])*1.02]
;  oplot,c_s[c_i],c_Te[c_i],color=2

  plot ,b_s[b_i],b_ne[b_i]
;  oplot,c_s[c_i],c_ne[c_i],color=2

  !P.MULTI = 0

END
;
; ======================================================================
;
; This is just here to make IDL happy, i.e. to have a subroutine in the
; file that has the same name as the file itself.
;  
PRO interface
END
;
; ======================================================================
;
