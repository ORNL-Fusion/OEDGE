FUNCTION GETCOLOR, thisColor, TRUE=truecolor
     
      ; Set up the color vectors.
      
   names   = ['CHARCOAL', 'RED', 'GREEN', 'BLUE', 'YELLOW']
   rvalue  = [    70,      255,      0,      0,     255   ]
   gvalue  = [    70,        0,    255,      0,     255   ]
   bvalue  = [    70,        0,      0,    255,       0   ]

      ; Did the user ask for a specific color? If not, return
      ; all the colors. If the user asked for a specific color,
      ; find out if a 24-bit value is required. Return to main
      ; IDL level if an error occurs.
      
   ON_ERROR, 1
   np = N_PARAMS()
   IF np EQ 1 THEN BEGIN

         ; Make sure the parameter is an uppercase string.
      
      varInfo = SIZE(thisColor)
      IF varInfo(varInfo(0) + 1) NE 7 THEN $
         MESSAGE, 'The color name must be a string.'
      thisColor = STRUPCASE(thisColor)
      
         ; Get the color triple for this color.
         
      colorIndex = WHERE(names EQ thisColor)

         ; If you can't find it. Issue an infomational message,
         ; set the index to a YELLOW color, and continue.
         
      IF colorIndex(0) LT 0 THEN BEGIN
         MESSAGE, "Can't find color. Returning yellow.", /INFORMATIONAL
         colorIndex = 4
      ENDIF
      
         ; Get the color triple.
      
      r = rvalue(colorIndex)
      g = gvalue(colorIndex)
      b = bvalue(colorIndex)
      returnColor = REFORM([r, g, b], 1, 3)
      
         ; Did the user want a 24-bit value? If so, call COLOR24.
         
      IF KEYWORD_SET(trueColor) THEN returnColor = COLOR24(returnColor)
      RETURN, returnColor
   ENDIF

      ; If you got here. Return all the colors.
      
   RETURN, REFORM([rvalue, gvalue, bvalue], 5, 3)
   END