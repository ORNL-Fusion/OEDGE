;-----------------------------------------------------------------------------------
; Routine: TRUECOLOR
; Version: 1.02
; Author : R.Martin
; Date   : 07.09.05
;-----------------------------------------------------------------------------------
; TRUECOLOR
;
; Function returns the true color RGB value for a named color
;
; Calling sequence:
;
; 1) Using a string value
;
;   col=TRUECOLOR(colorname)
;
; col       : long integer cotaining RGB color index.
; colorname : string containing any one of the set color name, names may be 
;             abrieviated.
;
; 
; 2) Using a keyword. Alternatively the code can be called using a keyword
;
;   col=TRUECOLOR(/colorname)
;
;
; 3) Calculate rgb value
;
;   col=TRUECOLOR(r, g, b)
;   col=TRUECOLOR([r,g,b])
;
;
; Example Code
;
;   x=findgen(100)*0.1
;   plot, x, x*sin(x), col=TRUECOLOR('Green')
;   oplot, x, x*cos(x), col=TRUCOLOR(/red)
;
;-----------------------------------------------------------------------------------
;

function truecolor, arg, arg2, arg3,  _extra=extra,       $
                    help=showlist,               $
                    loadct=setcolortable

  common loc_truecolor, list
  common colors, r1, g1, b1, r2, g2, b2

  if (n_elements(r1) gt 0) then begin
    r2=r1
    g2=g1
    b2=b1
  endif

  if (n_elements(list) eq 0) then begin
    ncol=141
    list=replicate({name:'', hex:'', r:0B, g:0B, b:0B, rgb:0L}, ncol)
    tmp=strarr(ncol)
    tmp(0)='Aliceblue F0F8FF'
    tmp(1)='Antiquewhite FAEBD7'
    tmp(2)='Aqua 00FFFF'
    tmp(3)='Aquamarine 7FFFD4'
    tmp(4)='Azure F0FFFF'
    tmp(5)='Beige F5F5DC'
    tmp(6)='Bisque FFE4C4'
    tmp(7)='Black 000000'
    tmp(8)='Blanchedalmond FFEBCD'
    tmp(9)='Blue 0000FF'
    tmp(10)='Blueviolet 8A2BE2'
    tmp(11)='Brown A52A2A'
    tmp(12)='Burlywood DEB887'
    tmp(13)='Cadetblue 5F9EA0'
    tmp(14)='Chartreuse 7FFF00'
    tmp(15)='Chocolate D2691E'
    tmp(16)='Coral FF7F50'
    tmp(17)='Cornflowerblue 6495ED'
    tmp(18)='Cornsilk FFF8DC'
    tmp(19)='Crimson DC143C'
    tmp(20)='Cyan 00FFFF'
    tmp(21)='Darkblue 00008B'
    tmp(22)='Darkcyan 008B8B'
    tmp(23)='Darkgoldenrod B8860B'
    tmp(24)='Darkgray A9A9A9'
    tmp(25)='Darkgreen 006400'
    tmp(26)='Darkkhaki BDB76B'
    tmp(27)='Darkmagenta 8B008B'
    tmp(28)='Darkolivegreen 556B2F'
    tmp(29)='Darkorange FF8C00'
    tmp(30)='Darkorchid 9932CC'
    tmp(31)='Darkred 8B0000'
    tmp(32)='Darksalmon E9967A'
    tmp(33)='Darkseagreen 8FBC8F'
    tmp(34)='Darkslateblue 483D8B'
    tmp(35)='Darkslategray 2F4F4F'
    tmp(36)='Darkturquoise 00CED1'
    tmp(37)='Darkviolet 9400D3'
    tmp(38)='Deeppink FF1493'
    tmp(39)='Deepskyblue 00BFFF'
    tmp(40)='Dimgray 696969'
    tmp(41)='Dodgerblue 1E90FF'
    tmp(42)='Firebrick B22222'
    tmp(43)='Floralwhite FFFAF0'
    tmp(44)='Forestgreen 228B22'
    tmp(45)='Fuchsia FF00FF'
    tmp(46)='Gainsboro DCDCDC'
    tmp(47)='Ghostwhite F8F8FF'
    tmp(48)='Gold FFD700'
    tmp(49)='Goldenrod DAA520'
    tmp(50)='Gray 808080'
    tmp(51)='Green 008000'
    tmp(52)='Greenyellow ADFF2F'
    tmp(53)='Honeydew F0FFF0'
    tmp(54)='Hotpink FF69B4'
    tmp(55)='Indianred CD5C5C'
    tmp(56)='Indigo 4B0082'
    tmp(57)='Ivory FFFFF0'
    tmp(58)='Khaki F0E68C'
    tmp(59)='Lavender E6E6FA'
    tmp(60)='Lavenderblush FFF0F5'
    tmp(61)='Lawngreen 7CFC00'
    tmp(62)='Lemonchiffon FFFACD'
    tmp(63)='Lightblue ADD8E6'
    tmp(64)='Lightcoral F08080'
    tmp(65)='Lightcyan E0FFFF'
    tmp(66)='Lightgoldenrodyellow FAFAD2'
    tmp(67)='Lightgreen 90EE90'
    tmp(68)='Lightgrey D3D3D3'
    tmp(69)='Lightpink FFB6C1'
    tmp(70)='Lightsalmon FFA07A'
    tmp(71)='Lightseagreen 20B2AA'
    tmp(72)='Lightskyblue 87CEFA'
    tmp(73)='Lightslategray 778899'
    tmp(74)='Lightsteelblue B0C4DE'
    tmp(75)='Lightyellow FFFFE0'
    tmp(76)='Lime 00FF00'
    tmp(77)='Limegreen 32CD32'
    tmp(78)='Linen FAF0E6'
    tmp(79)='Magenta FF00FF'
    tmp(80)='Maroon 800000'
    tmp(81)='Mediumauqamarine 66CDAA'
    tmp(82)='Mediumblue 0000CD'
    tmp(83)='Mediumorchid BA55D3'
    tmp(84)='Mediumpurple 9370D8'
    tmp(85)='Mediumseagreen 3CB371'
    tmp(86)='Mediumslateblue 7B68EE'
    tmp(87)='Mediumspringgreen 00FA9A'
    tmp(88)='Mediumturquoise 48D1CC'
    tmp(89)='Mediumvioletred C71585'
    tmp(90)='Midnightblue 191970'
    tmp(91)='Mintcream F5FFFA'
    tmp(92)='Mistyrose FFE4E1'
    tmp(93)='Moccasin FFE4B5'
    tmp(94)='Navajowhite FFDEAD'
    tmp(95)='Navy 000080'
    tmp(96)='Oldlace FDF5E6'
    tmp(97)='Olive 808000'
    tmp(98)='Olivedrab 688E23'
    tmp(99)='Orange FFA500'
    tmp(100)='Orangered FF4500'
    tmp(101)='Orchid DA70D6'
    tmp(102)='Palegoldenrod EEE8AA'
    tmp(103)='Palegreen 98FB98'
    tmp(104)='Paleturquoise AFEEEE'
    tmp(105)='Palevioletred D87093'
    tmp(106)='Papayawhip FFEFD5'
    tmp(107)='Peachpuff FFDAB9'
    tmp(108)='Peru CD853F'
    tmp(109)='Pink FFC0CB'
    tmp(110)='Plum DDA0DD'
    tmp(111)='Powderblue B0E0E6'
    tmp(112)='Purple 800080'
    tmp(113)='Red FF0000'
    tmp(114)='Rosybrown BC8F8F'
    tmp(115)='Royalblue 4169E1'
    tmp(116)='Saddlebrown 8B4513'
    tmp(117)='Salmon FA8072'
    tmp(118)='Sandybrown F4A460'
    tmp(119)='Seagreen 2E8B57'
    tmp(120)='Seashell FFF5EE'
    tmp(121)='Sienna A0522D'
    tmp(122)='Silver C0C0C0'
    tmp(123)='Skyblue 87CEEB'
    tmp(124)='Slateblue 6A5ACD'
    tmp(125)='Slategray 708090'
    tmp(126)='Snow FFFAFA'
    tmp(127)='Springgreen 00FF7F'
    tmp(128)='Steelblue 4682B4'
    tmp(129)='Tan D2B48C'
    tmp(130)='Teal 008080'
    tmp(131)='Thistle D8BFD8'
    tmp(132)='Tomato FF6347'
    tmp(133)='Turquoise 40E0D0'
    tmp(134)='Violet EE82EE'
    tmp(135)='Wheat F5DEB3'
    tmp(136)='White FFFFFF'
    tmp(137)='Whitesmoke F5F5F5'
    tmp(138)='Yellow FFFF00'
    tmp(139)='YellowGreen 9ACD32'
    tmp(140)='#BackGround# 000000'

    tmp=strupcase(tmp)
    
    tmp=strsplit(strjoin(strupcase(tmp), ' '), /extract)
    tmp=transpose(reform(tmp, 2, ncol))
    list.name=tmp(*,0)
    list.hex=tmp(*,1)

    rgb=list.hex
    list.hex=strmid(rgb, 4, 2)+strmid(rgb, 2, 2)+strmid(rgb, 0, 2)
    rgb=list.rgb
    reads, list.hex, rgb, format='(z)'
    list.rgb=rgb
    
    list.b=rgb/65536
    list.r=rgb mod 256
    list.g=(rgb/256) mod 256

    list=temporary(shift(list, 1))
  endif

  if keyword_set(setcolortable) or (!d.name eq 'PS') then begin
    tvlct, r, g, b, /get
    if ~array_equal([r,g,b], [list.r,list.g,list.b]) then tvlct, list.r, list.g, list.b
  endif

  if keyword_set(showlist) then return, list[1:*].name

  if is_integer(arg) then begin
    if undefined(arg2) and undefined(arg3) then begin
      if (n_elements(arg) ge 3) then begin
        rgb=byte(arg)
        return, rgb[0]+256L*rgb[1]+65536L*rgb[2]
      endif else return, arg
    endif

    if not_integer(arg2) or not_integer(arg3) then return, !p.color

    return, byte(arg)+256L*byte(arg2)+65536L*byte(arg3)
  endif

  if not_string(arg) then begin
    if (n_elements(extra) eq 0) then return, !p.color

    name=tag_names(extra)
    if (n_elements(name) gt 1) then $
      print, 'Warning TRUECOLOR will return result in alphabetic order'

    name=name[0]
      
  endif else name=strupcase(arg)

  if (n_elements(name) eq 1) then rgb=0 else rgb=lonarr(n_elements(name))
  for i=0, n_elements(name)-1 do begin
    pos=where(strpos(list.name, name[i]) eq 0, flag)
    if (flag ne 0) then rgb[i]=pos[0]
  endfor

  if (!d.name eq 'X') or (!d.name eq 'WIN') then begin
    device, get_decompose=decompsetting
    if (decompsetting eq 0) then begin
      tvlct, r, g, b, /get
      if ~array_equal([r,g,b], [list.r,list.g,list.b]) then tvlct, list.r, list.g, list.b
    endif
    if (decompsetting) then return, list[rgb].rgb
  endif

  return, rgb

end

;-----------------------------------------------------------------------------------
; Modification history
;
; Version 1.02 08.11.05
;  - Returns !p.color if no other answer can be found
;  - Allow RGB values to be input - arrays possible
;  - Return multiple colour indicies
;
;
;
;
; Vesrion 1.01 - R.Martin 
;   07.09.05 - Automatically load color table if DECOMP=0
;            - Make color 0 black to return default IDL behaviour
;
; Version 1.00 - R.Martin
;


    
