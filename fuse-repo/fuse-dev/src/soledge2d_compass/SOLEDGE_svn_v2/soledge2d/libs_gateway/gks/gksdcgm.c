/*
 *
 * Copyright @ 1984 - 1999   Jochen Werner, Josef Heinen
 *
 * Permission to use, copy, and distribute this software and its
 * documentation for any purpose with or without fee is hereby granted,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.
 *
 * Permission to modify the software is granted, but not the right to
 * distribute the modified code.  Modifications are to be distributed
 * as patches to released version.
 *
 * This software is provided "as is" without express or implied warranty.
 *
 * Send your comments or suggestions to
 *  J.Heinen@kfa-juelich.de.
 *
 *
 * FACILITY:
 *
 *	GLI GKS V4.5
 *
 * ABSTRACT:
 *
 *	This module contains a logical device driver for CGM (Computer
 *	Graphics Metafile).
 *
 * AUTHOR:
 *
 *	Jochen Werner, Josef Heinen
 *
 * VERSION:
 *
 *	V1.0
 *
 */



#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#ifndef MSDOS
#include <sys/types.h>
#endif

#ifdef VMS
#include <descrip.h>
#endif

#ifdef cray
#include <fortran.h>
#endif

#ifdef VMS
#define CHARARG(a) struct dsc$descriptor *a
#else
#ifdef cray
#define CHARARG(a) _fcd a
#else
#if defined(_WIN32) && !defined(__GNUC__)
#define CHARARG(a) char *(a), unsigned short a##_len
#else
#define CHARARG(a) char *a
#endif /* _WIN32 */
#endif /* cray */
#endif /* VMS */

#if defined(_WIN32) && !defined(__GNUC__)
#define STDCALL __stdcall
#else
#define STDCALL
#endif

#include "cgm.h"
#include "gksdefs.h"


#define TRUE 1
#define FALSE 0
#define NIL 0


#define odd(number)   ((number) & 01)
#define nint(a)       ((int)((a) + 0.5))

#define Color8Bit(c) \
    (c) >= 588 ? 80 + (c-588)/56 * 12 + nint((c-588)%56 * 11.0/56.0) : \
    (c) >= 257 ? 8 + nint((c-257)/330.0 * (72-1)) : (c)

#if defined(__cplusplus) || defined(c_plusplus)
#define CGM_FUNC (void (*)(...))
#else
#define CGM_FUNC (void (*)())
#endif


typedef struct {
    float a, b, c, d;
    } norm_xform;

typedef struct {
    int type;
    float width;
    int color;} line_attributes;

typedef struct {
    int type;
    float width;
    int color;} marker_attributes;

typedef struct {
    int font;
    int prec;
    float expfac;
    float spacing;
    int color;
    float height;
    int upx;
    int upy;
    int path;
    int halign;
    int valign;} text_attributes;

typedef struct {
    int intstyle;
    int hatch_index;
    int pattern_index;
    int color;} fill_attributes;

typedef struct {
    norm_xform xform;			    /* internal transformation */
    line_attributes pline;		    /* current polyline attributes */
    marker_attributes pmark;		    /* current marker attributes */
    text_attributes text;		    /* current text attributes */
    fill_attributes fill;		    /* current fill area attributes */
    int buffer_ind;			    /* output buffer index */
    char buffer[max_buffer+2];		    /* output buffer */
    float color_t[max_colorxs*3];	    /* color table */
    int conid;				    /* GKS connection id */
    unsigned active;			    /* indicates active workstation */
    unsigned begin_page;                    /* indicates begin page */
    float vp[4];			    /* current GKS viewport */
    float wn[4];			    /* current GKS window */
    int xext, yext;                         /* VDC extent */
    float mm;				    /* metric size in mm */
    char cmd_buffer[hdr_long + max_long];   /* where we buffer output */
    char *cmd_hdr;                          /* the command header */
    char *cmd_data;                         /* the command data */
    int cmd_index;                          /* index into the command data */
    int bfr_index;                          /* index into the buffer */
    int partition;                          /* which partition in the output */
    enum encode_enum encode;		    /* type of encoding */
#if defined(__cplusplus) || defined(c_plusplus)
    void (*cgm[n_melements])(...);	    /* cgm functions and procedures */
#else
    void (*cgm[n_melements])();		    /* cgm functions and procedures */
#endif
    } cgm_context;

static cgm_context *p;

static char digits[10] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};

static char *fonts[max_std_textfont] = {
    "AvantGarde-Book", "Courier", "Helvetica", "LubalinGraph-Book",
    "NewCenturySchlbk-Roman", "Souvenir-Light", "Symbol", "Times-Roman",
    "AvantGarde-Demi", "Courier-Bold", "Helvetica-Bold", "LubalinGraph-Demi",
    "NewCenturySchlbk-Bold", "Souvenir-Demi", "Symbol", "Times-Bold",
    "AvantGarde-BookOblique", "Courier-Oblique", "Helvetica-Oblique",
    "LubalinGraph-BookOblique", "NewCenturySchlbk-Italic",
    "Souvenir-LightItalic", "Symbol", "Times-Italic", "AvantGarde-DemiOblique",
    "Courier-BoldOblique", "Helvetica-BoldOblique", "LubalinGraph-DemiOblique",
    "NewCenturySchlbk-BoldItalic", "Souvenir-DemiItalic", "Symbol",
    "Times-BoldItalic"
    };

/* Transform world coordinates to virtual device coordinates */

static void WC_to_VDC (float xin, float yin, int *xout, int *yout)
{
    float x, y;

/* Normalization transformation */

    x = p->xform.a * xin + p->xform.b;
    y = p->xform.c * yin + p->xform.d;

/* Segment transformation */

    GST (&x, &y);

/* Virtual device transformation */

    *xout = (int)(x * max_coord);
    *yout = (int)(y * max_coord);
}



/* Flush output buffer */

static void cgmt_fb (void)
{
#ifdef VMS
    struct dsc$descriptor_s text;
#else
#ifdef cray
    _fcd text;
#endif /* cray */
#endif /* VMS */

    if (p->buffer_ind != 0)
	{
	p->buffer[p->buffer_ind++] = '\n';
	p->buffer[p->buffer_ind] = '\0';
#ifdef VMS
	text.dsc$b_dtype = DSC$K_DTYPE_T;
	text.dsc$b_class = DSC$K_CLASS_S;
	text.dsc$w_length = p->buffer_ind;
	text.dsc$a_pointer = p->buffer;

	BUFOUT (&p->conid, &p->buffer_ind, &text);
#else
#ifdef cray
        text = _cptofcd(p->buffer, p->buffer_ind);

	BUFOUT (&p->conid, &p->buffer_ind, text);
#else
	BUFOUT (&p->conid, &p->buffer_ind, p->buffer,
	    (unsigned short)p->buffer_ind);
#endif /* cray */
#endif /* VMS */

	p->buffer_ind = 0;
	p->buffer[0] = '\0';
	}
}



/* Write a character to CGM clear text */

static void cgmt_outc (char chr)
{
    if (p->buffer_ind >= cgmt_recl)
	cgmt_fb ();

    p->buffer[p->buffer_ind++] = chr;
    p->buffer[p->buffer_ind] = '\0';
}



/* Write string to CGM clear text */

static void cgmt_out_string (char *string)
{
    if ((int) (p->buffer_ind + strlen (string)) >= cgmt_recl)
	{
	cgmt_fb ();
	strcpy (p->buffer, "   ");
	p->buffer_ind = 3;
	}

    strcat (p->buffer, string);
    p->buffer_ind = p->buffer_ind + strlen (string);
}



/* Start output command */

static void cgmt_start_cmd (int cl, int el)
{
    cgmt_out_string (cgmt_cptr[cl][el]);
}



/* Flush output command */

static void cgmt_flush_cmd (int this_flush)
{
    cgmt_outc (term_char);

    cgmt_fb ();
}



/* Write a CGM clear text string */

static void cgmt_string (char *cptr, int slen)
{
    int i;

    cgmt_outc (' ');
    cgmt_outc (quote_char);

    for (i  = 0; i < slen; ++i) 
	{
	if (cptr[i] == quote_char) 
	    {
	    cgmt_outc (quote_char);
	    }

	cgmt_outc (cptr[i]);
	}

    cgmt_outc (quote_char);
}



/* Write a signed integer variable */

static void cgmt_int (int xin)
{
    static char buf[max_pwrs + 2];
    register char *cptr;
    register int is_neg, j;
      
    cptr = buf + max_pwrs + 1;
    *cptr = '\0';

    if (xin < 0) 
	{
	is_neg = 1;
	xin = -xin;
	} 
    else 
	is_neg = 0;

    if (xin == 0) 
	{
	*--cptr = digits[0];

	if ((int) (p->buffer_ind + strlen (cptr)) < cgmt_recl)
	    cgmt_outc (' ');

	cgmt_out_string (cptr); /* all done */
	return;
	}

    while (xin) 
	{
	j = xin % 10;
	*--cptr = digits[j];
	xin /= 10;
	}

    if (is_neg) 
	*--cptr = '-';

    if ((int) (p->buffer_ind + strlen (cptr)) < cgmt_recl)
	cgmt_outc (' ');

    cgmt_out_string (cptr);
}



/* Write a real variable */

static void cgmt_real (float xin)
{
    char buffer[max_str];

    sprintf (buffer, " %.6f", xin);
    cgmt_out_string (buffer);
}



/* Write an integer point */

static void cgmt_ipoint (int x, int y)
{
    char buffer[max_str];

    sprintf (buffer, " %d,%d", x, y);
    cgmt_out_string (buffer);
}



/* Begin metafile */

static void cgmt_begin (char *comment)
{
    cgmt_start_cmd (0, (int) B_Mf);

    if (*comment) 
	cgmt_string (comment, strlen (comment));
    else 
	cgmt_string (NULL, 0);

    cgmt_flush_cmd (final_flush);
}



/* End metafile */

static void cgmt_end (void)
{
    cgmt_start_cmd (0, (int) E_Mf);

    cgmt_flush_cmd (final_flush);

    cgmt_fb ();
}



/* Begin picture */

static void cgmt_bp (char *pic_name)
{
    cgmt_start_cmd (0, (int) B_Pic);

    if (*pic_name) 
	cgmt_string (pic_name, strlen (pic_name));
    else 
	cgmt_string (NULL, 0);

    cgmt_flush_cmd (final_flush);
}



/* Begin picture body */

static void cgmt_bpage (void)
{
    cgmt_start_cmd (0, (int) B_Pic_Body);

    cgmt_flush_cmd (final_flush);
}



/* End picture */

static void cgmt_epage (void)
{
    cgmt_start_cmd (0, (int) E_Pic);

    cgmt_flush_cmd (final_flush);
}



/* Metafile version */

static void cgmt_mfversion (void)
{
    cgmt_start_cmd (1, (int) MfVersion);

    cgmt_int (1);

    cgmt_flush_cmd (final_flush);
}



/* Metafile description */

static void cgmt_mfdescrip (void)
{
    char *descrip;

    cgmt_start_cmd (1, MfDescrip);

    descrip = "GLI GKS V4.5 CGM Clear Text";
    cgmt_string (descrip, strlen (descrip));

    cgmt_flush_cmd (final_flush);
}



/* VDC type */

static void cgmt_vdctype (void)
{
    cgmt_start_cmd (1, (int) vdcType);

    cgmt_out_string (" Integer");

    cgmt_flush_cmd (final_flush);
}



/* Integer precision */

static void cgmt_intprec (void)
{
    cgmt_start_cmd (1, (int) IntPrec);

    cgmt_int (-32768);
    cgmt_int (32767);

    cgmt_flush_cmd (final_flush);
}



/* Real precision */

static void cgmt_realprec (void)
{

    cgmt_start_cmd (1, (int) RealPrec);

    cgmt_real (-32768.);
    cgmt_real (32768.);
    cgmt_int (4);

    cgmt_flush_cmd (final_flush);
}



/* Index precision */

static void cgmt_indexprec (void)
{
    cgmt_start_cmd (1, (int) IndexPrec);

    cgmt_int (-32768);
    cgmt_int (32767);

    cgmt_flush_cmd (final_flush);
}



/* Colour precision */

static void cgmt_colprec (void)
{
    cgmt_start_cmd (1, (int) ColPrec);

    cgmt_int ((1 << cprec) -1);

    cgmt_flush_cmd (final_flush);
}



/* Colour index precision */

static void cgmt_cindprec (void)
{
    cgmt_start_cmd (1, (int) CIndPrec);

    cgmt_int ((1 << cxprec) -1);

    cgmt_flush_cmd (final_flush);
}



/* Colour value extent */

static void cgmt_cvextent (void)
{
    int i;

    cgmt_start_cmd (1, (int) CVExtent);

    for (i = 0; i < 3; ++i) 
	{
	cgmt_int (0);
	}

    for (i = 0; i < 3; ++i) 
	{
	cgmt_int (max_colors-1);
	}

    cgmt_flush_cmd (final_flush);
}



/* Maximum colour index */

static void cgmt_maxcind (void)
{
    cgmt_start_cmd (1, (int) MaxCInd);

    cgmt_int (max_colorxs-1);

    cgmt_flush_cmd (final_flush);
}



/* Metafile element list */

static void cgmt_mfellist (void)
{
    int i;

    cgmt_start_cmd(1, (int) MfElList);

    cgmt_outc (' ');
    cgmt_outc (quote_char);

    for (i = 2; i < 2 * n_melements; i += 2) 
	{
	cgmt_out_string (cgmt_cptr[element_list[i]][element_list[i+1]]);
	
	if (i < 2 * (n_melements - 1)) 
	    cgmt_outc (' ');
	}

    cgmt_outc (quote_char);
    cgmt_flush_cmd (final_flush);
}



/* Font List */

static void cgmt_fontlist (void)
{
    int i;
    char s[max_str];

    cgmt_start_cmd(1, (int) FontList);

    cgmt_outc (' ');

    for (i = 0; i < max_std_textfont; i++) 
	{
	sprintf (s, "'%s'%s", fonts[i], (i < max_std_textfont-1) ? ", " : "");
	cgmt_out_string (s);
	}

    cgmt_flush_cmd (final_flush);
}



/* Character announcer */

static void cgmt_cannounce (void)
{
    cgmt_start_cmd (1, (int) CharAnnounce);

    cgmt_out_string (" Extd8Bit"); 

    cgmt_flush_cmd (final_flush);
}



/* Scaling mode */

static void cgmt_scalmode (void)
{
    cgmt_start_cmd (2, (int) ScalMode);

    if (p->mm > 0)
	{
	cgmt_out_string (" Metric");
	cgmt_real (p->mm);
	}
    else
	cgmt_out_string (" Abstract");

    cgmt_flush_cmd (final_flush);
}



/* Colour selection mode */

static void cgmt_colselmode (void)
{
    cgmt_start_cmd (2, (int) ColSelMode);

    cgmt_out_string (" Indexed"); 

    cgmt_flush_cmd (final_flush);
}



/* Line width specification mode */

static void cgmt_lwsmode (void)
{
    cgmt_start_cmd (2, (int) LWidSpecMode);

    cgmt_out_string (" Scaled"); 

    cgmt_flush_cmd (final_flush);
}



/* Marker size specification mode */

static void cgmt_msmode (void)
{
    cgmt_start_cmd (2, (int) MarkSizSpecMode);

    cgmt_out_string (" Scaled"); 

    cgmt_flush_cmd (final_flush);
}



/* VDC extent */

static void cgmt_vdcextent (void)
{
    cgmt_start_cmd (2, (int) vdcExtent);

    cgmt_ipoint (0, 0);
    cgmt_ipoint (p->xext, p->yext);

    cgmt_flush_cmd (final_flush);
}



/* Background colour */

static void cgmt_backcol (void)
{
    cgmt_start_cmd (2, (int) BackCol);

    cgmt_int (255);
    cgmt_int (255);
    cgmt_int (255);

    cgmt_flush_cmd (final_flush);
}



/* VDC integer precision */

static void cgmt_vdcintprec (void)
{
    cgmt_start_cmd (3, (int) vdcIntPrec);

    cgmt_int (-32768);
    cgmt_int (32767);

    cgmt_flush_cmd (final_flush);
}



/* Clip rectangle */

static void cgmt_cliprect (int *int_coords)
{
    cgmt_start_cmd (3, (int) ClipRect);

    cgmt_ipoint (int_coords[0], int_coords[1]);
    cgmt_ipoint (int_coords[2], int_coords[3]); 

    cgmt_flush_cmd (final_flush);
}



/* Clip indicator */

static void cgmt_clipindic (int clip_ind)
{
    cgmt_start_cmd (3, (int) ClipIndic);

    if (clip_ind)
	{
	cgmt_out_string (" On"); 
	}
    else
	{
	cgmt_out_string (" Off"); 
	}

    cgmt_flush_cmd (final_flush);
}



/* Polyline */

static void cgmt_pline (int no_pairs, int *x1_ptr, int *y1_ptr)
{
    int i;

    cgmt_start_cmd (4, (int) PolyLine);

    for (i = 0; i < no_pairs; ++i) 
	{
	cgmt_ipoint (x1_ptr[i], y1_ptr[i]);
	}

    cgmt_flush_cmd (final_flush);
}



/* Polymarker */

static void cgmt_pmarker (int no_pairs, int *x1_ptr, int *y1_ptr)
{
    int i;

    cgmt_start_cmd (4, (int) PolyMarker);

    for (i = 0; i < no_pairs; ++i) 
	{
	cgmt_ipoint (x1_ptr[i], y1_ptr[i]);
	}   

    cgmt_flush_cmd (final_flush);
}



/* Text */

static void cgmt_text (int x, int y, int final, char *buffer)
{
    cgmt_start_cmd(4, (int) Text);

    cgmt_ipoint (x, y);

    if (final)
	{ 
	cgmt_out_string (" Final");
	}
    else 
	{
	cgmt_out_string (" NotFinal");
	}

    cgmt_string (buffer, strlen (buffer));

    cgmt_flush_cmd (final_flush);
}



/* Polygon */

static void cgmt_pgon (int no_pairs, int *x1_ptr, int *y1_ptr)
{
    int i;

    cgmt_start_cmd (4, (int) Polygon);

    for (i = 0; i < no_pairs; ++i) 
	{
	cgmt_ipoint(x1_ptr[i], y1_ptr[i]);
	}

    cgmt_flush_cmd (final_flush);
}



/* Line type */

static void cgmt_ltype (int line_type)
{
    cgmt_start_cmd (5, (int) LType);

    cgmt_int ((int) line_type);

    cgmt_flush_cmd (final_flush);
}



/* Line width */

#ifdef linux
static void cgmt_lwidth(rmul)
    float rmul;
#else
static void cgmt_lwidth(float rmul)
#endif
{
    cgmt_start_cmd (5, (int) LWidth);

    cgmt_real (rmul); 

    cgmt_flush_cmd (final_flush);
}



/* Line colour */

static void cgmt_lcolour (int index)
{
    cgmt_start_cmd(5, (int) LColour);

    cgmt_int (index); 

    cgmt_flush_cmd (final_flush);
}



/* Marker type */

static void cgmt_mtype (int marker)
{
    cgmt_start_cmd (5, (int) MType);

    cgmt_int (marker);

    cgmt_flush_cmd (final_flush);
}



/* Marker size */

#ifdef linux
static void cgmt_msize(rmul)
    float rmul;
#else
static void cgmt_msize(float rmul)
#endif
{
    cgmt_start_cmd (5, (int) MSize);

    cgmt_real (rmul); 

    cgmt_flush_cmd (final_flush);
}



/* Marker colour */

static void cgmt_mcolour (int index)
{
    cgmt_start_cmd (5, (int) MColour);

    cgmt_int (index); 

    cgmt_flush_cmd (final_flush);
}



/* Text font index */

static void cgmt_tfindex (int index)
{
    cgmt_start_cmd (5, (int) TFIndex);

    cgmt_int (index);

    cgmt_flush_cmd (final_flush);
}



/* Text precision */

static void cgmt_tprec (int precision)
{
    cgmt_start_cmd (5, (int) TPrec);

    switch (precision) 
	{

	case string:
	    {   
	    cgmt_out_string (" String"); 
	    break;
	    }

	case character: 
	    {
	    cgmt_out_string (" Character"); 
	    break;
	    }

	case stroke:    
	    {
	    cgmt_out_string (" Stroke"); 
	    break;
	    }
	}

    cgmt_flush_cmd (final_flush);
}



/* Character expansion factor */

#ifdef linux
static void cgmt_cexpfac(factor)
    float factor;
#else
static void cgmt_cexpfac(float factor)
#endif
{
    cgmt_start_cmd (5, (int) CExpFac);

    cgmt_real (factor);

    cgmt_flush_cmd (final_flush);
}



/* Character space */

#ifdef linux
static void cgmt_cspace(space)
    float space;
#else
static void cgmt_cspace(float space)
#endif
{
    cgmt_start_cmd (5, (int) CSpace);

    cgmt_real (space);

    cgmt_flush_cmd (final_flush);
}



/* Text colour */

static void cgmt_tcolour (int index)
{
    cgmt_start_cmd (5, (int) TColour);

    cgmt_int (index); 
	
    cgmt_flush_cmd (final_flush);
}



/* Character height */

static void cgmt_cheight (int height)
{
    cgmt_start_cmd (5, (int) CHeight);

    cgmt_int (height);

    cgmt_flush_cmd (final_flush);
}



/* Character orientation */

static void cgmt_corient (int x_up, int y_up, int x_base, int y_base)
{
    cgmt_start_cmd (5, (int) COrient);

    cgmt_int (x_up);
    cgmt_int (y_up);
    cgmt_int (x_base);
    cgmt_int (y_base);

    cgmt_flush_cmd (final_flush);
}



/* Text path */

static void cgmt_tpath (int new_path)
{
    cgmt_start_cmd (5, (int) TPath);

    switch (new_path) 
	{
	case right:     
	    {
	    cgmt_out_string (" Right"); 
	    break;
	    }

	case left:
	    {   
	    cgmt_out_string (" Left"); 
	    break;
	    }

	case up:
	    {   
	    cgmt_out_string (" Up"); 
	    break;
	    }

	case down:
	    {   
	    cgmt_out_string (" Down"); 
	    break;
	    }
    }

    cgmt_flush_cmd(final_flush);
}




/* Text alignment */

static void cgmt_talign (int hor, int ver)
{
    cgmt_start_cmd (5, (int) TAlign);

    switch (hor) 
	{
	case normal_h:
	    {
	    cgmt_out_string (" NormHoriz"); 
	    break;
	    }

	case left_h:
	    {
	    cgmt_out_string (" Left"); 
	    break;
	    }

	case center_h:
	    {
	    cgmt_out_string (" Ctr"); 
	    break;
	    }

	case right_h:
	    {
	    cgmt_out_string (" Right"); 
	    break;
	    }

	case cont_h:    
	    {
	    cgmt_out_string (" ContHoriz"); 
	    break;
	    }
	}

    switch (ver) 
	{

	case normal_v:  
	    {
	    cgmt_out_string (" NormVert"); 
	    break;
	    }

	case top_v:     
	    {
	    cgmt_out_string (" Top");    
	    break;
	    }

	case cap_v:     
	    {
	    cgmt_out_string (" Cap"); 
	    break;
	    }

	case half_v:    
	    {
	    cgmt_out_string (" Half"); 
	    break;
	    }

	case base_v:    
	    {
	    cgmt_out_string (" Base"); 
	    break;
	    }

	case bottom_v:  
	    {
	    cgmt_out_string (" Bottom"); 
	    break;
	    }

	case cont_v:    
	    {
	    cgmt_out_string (" ContVert"); 
	    break;
	    }
	}

    cgmt_real (0.0);
    cgmt_real (0.0);

    cgmt_flush_cmd (final_flush);
}



/* Interior style */

static void cgmt_intstyle (int style)
{
    cgmt_start_cmd (5, (int) IntStyle);

    switch (style) 
	{

	case hollow:
	    {
	    cgmt_out_string (" Hollow"); 
	    break;
	    }

	case solid_i:   
	    {
	    cgmt_out_string (" Solid"); 
	    break;
	    }

	case pattern:   
	    {
	    cgmt_out_string (" Pat"); 
	    break;
	    }

	case hatch:
	    {
	    cgmt_out_string (" Hatch"); 
	    break;
	    }

	case empty:     
	    {
	    cgmt_out_string (" Empty"); 
	    break;
	    }
	}

    cgmt_flush_cmd (final_flush);
}



/* Fill colour */

static void cgmt_fillcolour (int index)
{
    cgmt_start_cmd (5, (int) FillColour);

    cgmt_int (index); 

    cgmt_flush_cmd (final_flush);
}



/* Hatch index */

static void cgmt_hindex (int new_index)
{
    cgmt_start_cmd (5, (int) HatchIndex);
    cgmt_int (new_index);
    cgmt_flush_cmd (final_flush);
}



/* Pattern index */

static void cgmt_pindex (int new_index)
{
    cgmt_start_cmd (5, (int) PatIndex);
    cgmt_int (new_index);
    cgmt_flush_cmd (final_flush);
}



/* Colour table */

static void cgmt_coltab (int beg_index, int no_entries, float *ctab)
{
    int i, j;

    cgmt_start_cmd (5, (int) ColTab);
    cgmt_int (beg_index);

    for (i = beg_index; i < (beg_index + no_entries); ++i) 
	{
	for (j = 0; j < 3; ++j) 
	    {
	    cgmt_int ((int)(ctab[(i - beg_index) * 3 + j] * (max_colors - 1)));
	    }
	}

    cgmt_flush_cmd (final_flush);
}



/* Cell array */

static void cgmt_carray (int xmin, int xmax, int ymin, int ymax, int dx, int dy,
    int dimx, int *array)
{
    int ix, iy, c;

    cgmt_start_cmd (4, (int) Cell_Array);

    cgmt_ipoint (xmin, ymin);
    cgmt_ipoint (xmax, ymax);
    cgmt_ipoint (xmax, ymin);
    cgmt_int (dx);
    cgmt_int (dy);
    cgmt_int (max_colors-1);

    for (iy = 0; iy < dy; iy++)
	{
	cgmt_fb ();

	for (ix = 0; ix < dx; ix++)
	    {
            c = array[dimx * iy + ix];
            c = Color8Bit(c);
	    cgmt_int (c);

	    if (ix < dx-1)
		cgmt_outc (',');
	    }
	}
    
    cgmt_flush_cmd (final_flush);
}



/* Flush output buffer */

static void cgmb_fb (void)

{
#ifdef VMS
    struct dsc$descriptor_s text;
#else
#ifdef cray
    _fcd text;
#endif /* cray */
#endif /* VMS */

    if (p->buffer_ind != 0)
	{
	p->buffer[p->buffer_ind] = '\0';
#ifdef VMS
	text.dsc$b_dtype = DSC$K_DTYPE_T;
	text.dsc$b_class = DSC$K_CLASS_S;
	text.dsc$w_length = p->buffer_ind;
	text.dsc$a_pointer = p->buffer;

	BUFOUT (&p->conid, &p->buffer_ind, &text);
#else
#ifdef cray
        text = _cptofcd(p->buffer, p->buffer_ind);

	BUFOUT (&p->conid, &p->buffer_ind, text);
#else
	BUFOUT (&p->conid, &p->buffer_ind, p->buffer,
	    (unsigned short)p->buffer_ind);
#endif /* cray */
#endif /* VMS */

	p->buffer_ind = 0;
	p->buffer[0] = '\0';
	}
}


/* Write one byte to buffer */

static void cgmb_outc (char chr)
{
    if (p->buffer_ind >= max_buffer)
	cgmb_fb ();

    p->buffer[p->buffer_ind++] = chr;
}



/* Start output command */

static void cgmb_start_cmd (int cl, int el)
{
#define cl_max 15
#define el_max 127

    p->cmd_hdr = p->cmd_buffer + p->bfr_index;
    p->cmd_data = p->cmd_hdr + hdr_long;
    p->bfr_index += hdr_long;

    p->cmd_hdr[0] = (cl << 4) | (el >>3);
    p->cmd_hdr[1] = el << 5;
    p->cmd_index = 0;
    p->partition = 1;

#undef cl_max
#undef el_max
}



/* Flush output command */

static void cgmb_flush_cmd (int this_flush)
{
    int i;

    if ((this_flush == final_flush) && (p->partition == 1) && 
	(p->cmd_index <= max_short)) 
	{
	p->cmd_hdr[1] |= p->cmd_index;
	
	/* flush out the header */

	for (i = 0; i < hdr_short; ++i) 
	    {
	    cgmb_outc (p->cmd_hdr[i]);
	    }

	}
    else 
	{
	/* need a long form */

	if (p->partition == 1) 
	    {   
	    /* first one */

	    p->cmd_hdr[1] |= 31;

	    for (i = 0; i < hdr_short; ++i) 
		{
		cgmb_outc (p->cmd_hdr[i]);
		}
	    }

	p->cmd_hdr[2] = p->cmd_index >> 8;
	p->cmd_hdr[3] = p->cmd_index & 255;

	if (this_flush == int_flush)
	    {
	    p->cmd_hdr[2] |= 1 << 7;    /* more come */
	    }

	/* flush out the header */

	for (i = hdr_short; i < hdr_long; ++i)
	    {   
	    cgmb_outc (p->cmd_hdr[i]);
	    }
	}


    /* now flush out the data */

    for (i = 0; i < p->cmd_index; ++i) 
	{
	cgmb_outc (p->cmd_data[i]);
	}

    if (p->cmd_index % 2) 
	{
	cgmb_outc ('\0');
	}

    p->cmd_index = 0;
    p->bfr_index = 0;
    ++p->partition;
}




/* Write one byte */

static void cgmb_out_bc (int c)
{
    if (p->cmd_index >= max_long) 
	{
	cgmb_flush_cmd (int_flush);
	}

    p->cmd_data[p->cmd_index++] = c;
}




/* Write multiple bytes */

static void cgmb_out_bs (char *cptr, int n)
{
    int to_do, space_left, i;

    to_do = n;
    space_left = max_long - p->cmd_index;

    while (to_do > space_left) 
	{
	for (i = 0; i < space_left; ++i) 
	    {
	    p->cmd_data[p->cmd_index++] = *cptr++;
	    }

	cgmb_flush_cmd (int_flush);
	to_do -= space_left;
	space_left = max_long;
	}

    for (i = 0; i < to_do; ++i) 
	{
	p->cmd_data[p->cmd_index++] = *cptr++;
	}
}




/* Write a CGM binary string */

static void cgmb_string (char *cptr, int slen)
{
    int to_do;
    unsigned char byte1, byte2;

    /* put out null strings, however */

    if (slen == 0) 
	{
	cgmb_out_bc (0);
	return;
	}

    /* now non-trivial stuff */

    if (slen < 255) 
	{
	/* simple case */

	cgmb_out_bc (slen);
	cgmb_out_bs (cptr, slen);
	}
    else 
	{
	cgmb_out_bc (255);
	to_do = slen;

	while (to_do > 0) 
	    {
	    if (to_do < max_long) 
		{
		/* last one */

		byte1 = to_do >> 8;
		byte2 = to_do & 255;

		cgmb_out_bc (byte1);
		cgmb_out_bc (byte2);
		cgmb_out_bs (cptr, to_do);

		to_do = 0;
		} 
	    else 
		{
		byte1 = (max_long >> 8) | (1 << 7);
		byte2 = max_long & 255;

		cgmb_out_bc (byte1);
		cgmb_out_bc (byte2);
		cgmb_out_bs (cptr, max_long);

		to_do -= max_long;
		}
	    }
	}
}




/* Write a signed integer variable */

static void cgmb_gint (int xin, int precision)
{

    int i, no_out, xshifted;
    char buffer[4];

    no_out = precision / byte_size;

    xshifted = xin;
    for (i = no_out - 1; i >= 0; --i)
	{
	buffer[i] = xshifted & byte_mask;
	xshifted >>= byte_size;
	}

    if ((xin < 0) && (buffer[0] > '\0'))        /* maybe after truncation */
	{
	buffer[0] |= 1 << 7;            /* assuming two's complement */
	}

    cgmb_out_bs (buffer, no_out);
}




/* Write an unsigned integer variable */

static void cgmb_uint (unsigned int xin, int precision)
{
    int i, no_out;
    unsigned char buffer[4];

    no_out = precision / byte_size;

    for (i = no_out - 1; i >= 0; --i) 
	{
	buffer[i] = xin & byte_mask;
	xin >>= byte_size;
	}

    cgmb_out_bs ((char *) buffer, no_out);
}



/* Write fixed point variable */

static void cgmb_fixed (float xin)
{
    int exp_part, fract_part;
    float fract_real;

    exp_part = (int)xin;
    if (exp_part > xin) 
	{
	exp_part -= 1;  /* take it below xin */
	}

    fract_real = xin - exp_part;
    fract_part = (int)(fract_real * (01 << real_prec_fract));

    cgmb_gint (exp_part, real_prec_exp);
    cgmb_uint (fract_part, real_prec_fract);
}



/* Write IEEE floating point variable */

static void cgmb_float (float xin)
{
    unsigned char arry[8];
    int sign_bit, i;
    unsigned int exponent;
    unsigned long fract;
    double dfract;

    if (xin < 0.0) 
	{
	sign_bit = 1;
	xin = -xin;
	}
    else 
	{
	sign_bit = 0;
	}

    /* first calculate the exponent and fraction */

    if (xin == 0.0) 
	{
	exponent = 0;
	fract = 0;
	} 
    else 
	{
	switch (real_prec_exp + real_prec_fract) 
	    {

	    /* first 32 bit precision */

	    case 32:
		{
		if (xin < 1.0) 
		    {
		    for (i = 0; (xin < 1.0) && (i < 128); ++i) 
			{
			xin *= 2.0;
			}

		    exponent = 127 - i;
		    } 
		else
		    {
		    if (xin >= 2.0) 
			{
			for (i=0; (xin >= 2.0) && (i < 127); ++i) 
			    {
			    xin /= 2.0;
			    }
			exponent = 127 + i;
			}
		    else 
			{
			exponent = 127;
			}
		    }
		dfract = xin - 1.0;
		
		for (i = 0; i < 23; ++i) 
		    {
		    dfract *= 2.0;
		    }

		fract = (unsigned long) dfract;
		break;
		}

		/* now 64 bit precision */

	    case 64:
		{
		break;
		}
	    }
	}

    switch (real_prec_exp + real_prec_fract) 
	{
    
	/* first 32 bit precision */

	case 32:
	    {
	    arry[0] = ((sign_bit & 1) << 7) | ((exponent >> 1) & 127);
	    arry[1] = ((exponent & 1) << 7) | ((fract >> 16) & 127);
	    arry[2] = (fract >> 8) & 255;
	    arry[3] = fract & 255;
	    cgmb_out_bs ((char *) arry, 4);
	    break;
	    }

	/* now 64 bit precision */

	case 64:
	    {
	    break;
	    }
	}
}



/* Write direct colour value */

static void cgmb_dcint (int xin)
{
    cgmb_uint (xin, cprec);
}



/* Write a signed int at VDC integer precision */

static void cgmb_vint (int xin)
{
    cgmb_gint (xin, 16);
}



/* Write a standard CGM signed int */

static void cgmb_sint (int xin)
{
    cgmb_gint (xin, 16);
}



/* Write a signed int at index precision */

static void cgmb_xint (int xin)
{
    cgmb_gint (xin, 16);
}



/* Write an unsigned integer at colour index precision */

static void cgmb_cxint (int xin)
{
    cgmb_uint((unsigned) xin, cxprec);
}



/* Write an integer at fixed (16 bit) precision */

static void cgmb_eint (int xin)
{
    char byte1;
    unsigned char byte2;

    byte1 = xin / 256;
    byte2 = xin & 255;
    cgmb_out_bc (byte1);
    cgmb_out_bc (byte2);
}



/* Begin metafile */

static void cgmb_begin (char *comment)
{
    cgmb_start_cmd (0, (int) B_Mf);

    if (*comment) 
	{
	cgmb_string (comment, strlen (comment));
	}
    else 
	{
	cgmb_string (NULL, 0);
	}

    cgmb_flush_cmd (final_flush);

    cgmb_fb ();
}



/* End metafile */

static void cgmb_end (void)
{
    /* put out the end metafile command */

    cgmb_start_cmd (0, (int) E_Mf);
    cgmb_flush_cmd (final_flush);

    /* flush out the buffer */

    cgmb_fb ();
}



/* Start picture */

static void cgmb_bp (char *pic_name)
{
    cgmb_start_cmd (0, (int) B_Pic);

    if (*pic_name) 
	{
	cgmb_string (pic_name, strlen (pic_name));
	}
    else 
	{
	cgmb_string (NULL, 0);
	}

    cgmb_flush_cmd (final_flush);
}



/* Start picture body */

static void cgmb_bpage (void)
{
    cgmb_start_cmd (0, (int) B_Pic_Body);

    cgmb_flush_cmd (final_flush);
}



/* End picture */

static void cgmb_epage (void)
{
    cgmb_start_cmd (0, (int) E_Pic);

    cgmb_flush_cmd (final_flush);
}



/* Metafile version */

static void cgmb_mfversion (void)
{
    cgmb_start_cmd (1, (int) MfVersion);

    cgmb_sint (1);

    cgmb_flush_cmd (final_flush);
}



/* Metafile description */

static void cgmb_mfdescrip (void)
{
    char *descrip;

    cgmb_start_cmd (1, MfDescrip);

    descrip = "GLI GKS V4.5 CGM Binary";
    cgmb_string (descrip, strlen (descrip));

    cgmb_flush_cmd (final_flush);
}



/* VDC type */

static void cgmb_vdctype (void)
{
    cgmb_start_cmd (1, (int) vdcType);

    cgmb_eint ((int) vdc_int);

    cgmb_flush_cmd (final_flush);
}



/* Integer precision */

static void cgmb_intprec (void)
{
    cgmb_start_cmd (1, (int) IntPrec);

    cgmb_sint (16);

    cgmb_flush_cmd (final_flush);
}



/* Real precision */

static void cgmb_realprec (void)
{
    cgmb_start_cmd (1, (int) RealPrec);

    cgmb_sint (1);
    cgmb_sint (16);
    cgmb_sint (16);

    cgmb_flush_cmd (final_flush);
}



/* Index precision */

static void cgmb_indexprec (void)
{
    cgmb_start_cmd (1, (int) IndexPrec);

    cgmb_sint (16);

    cgmb_flush_cmd (final_flush);
}



/* Colour precision */

static void cgmb_colprec (void)
{
    cgmb_start_cmd (1, (int) ColPrec);

    cgmb_sint (cprec);

    cgmb_flush_cmd (final_flush);
}



/* Colour index precision */

static void cgmb_cindprec (void)
{
    cgmb_start_cmd (1, (int) CIndPrec);

    cgmb_sint (cxprec);

    cgmb_flush_cmd (final_flush);
}



/* Colour value extent */

static void cgmb_cvextent (void)
{
    cgmb_start_cmd (1, (int) CVExtent);

    cgmb_dcint (0);
    cgmb_dcint (0);
    cgmb_dcint (0);
    cgmb_dcint (max_colors-1);
    cgmb_dcint (max_colors-1);
    cgmb_dcint (max_colors-1);

    cgmb_flush_cmd (final_flush);
}



/* Maximum colour index */

static void cgmb_maxcind (void)
{
    cgmb_start_cmd (1, (int) MaxCInd);

    cgmb_cxint (max_colorxs-1);

    cgmb_flush_cmd (final_flush);
}



/* Metafile element list */

static void cgmb_mfellist (void)
{
    int i;

    cgmb_start_cmd (1, (int) MfElList);
    cgmb_sint (n_melements);

    for (i = 2; i < 2 * n_melements; ++i) 
	{
	cgmb_xint (element_list[i]);
	}

    cgmb_flush_cmd (final_flush);
}



/* Font List */

static void cgmb_fontlist (void)
{
    register int i, slen;
    register char *s;

    for (i = 0, slen = 0; i < max_std_textfont; i++) 
	{
	slen += strlen (fonts[i]) + 1;
	}
    s = (char *) malloc (slen);

    *s = '\0';
    for (i = 0; i < max_std_textfont; i++) 
	{
	strcat (s, fonts[i]);
	if (i < max_std_textfont - 1)
	    strcat (s, " ");
	}

    cgmb_start_cmd (1, (int) FontList);

    cgmb_string (s, strlen (s));

    cgmb_flush_cmd (final_flush);

    free (s);
}



/* Character announcer */

static void cgmb_cannounce (void)
{
    cgmb_start_cmd (1, (int) CharAnnounce);

    cgmb_eint (3);

    cgmb_flush_cmd (final_flush);
}



/* Scaling mode */

static void cgmb_scalmode (void)
{
    cgmb_start_cmd (2, (int) ScalMode);

    if (p->mm > 0)
	{
	cgmb_eint (1);
	cgmb_float (p->mm);
	}
    else
	{
	cgmb_eint (0);
	cgmb_float (0.);
	}

    cgmb_flush_cmd (final_flush);
}



/* Colour selection mode */

static void cgmb_colselmode (void)
{
    cgmb_start_cmd (2, (int) ColSelMode);

    cgmb_eint ((int) i_c_mode);

    cgmb_flush_cmd (final_flush);
}




/* line width specification mode */

static void cgmb_lwsmode (void)
{
    cgmb_start_cmd (2, (int) LWidSpecMode);

    cgmb_eint ((int) scaled);

    cgmb_flush_cmd (final_flush);
}



/* marker size specification mode */

static void cgmb_msmode (void)
{
    cgmb_start_cmd (2, (int) MarkSizSpecMode);

    cgmb_eint ((int) scaled);

    cgmb_flush_cmd (final_flush);
}



/* VDC extent */

static void cgmb_vdcextent (void)
{
    cgmb_start_cmd (2, (int) vdcExtent);

    cgmb_vint (0); cgmb_vint (0); 
    cgmb_vint (p->xext); cgmb_vint (p->yext);

    cgmb_flush_cmd (final_flush);
}



/* Background colour */

static void cgmb_backcol (void)
{
    cgmb_start_cmd (2, (int) BackCol);

    cgmb_dcint (255);
    cgmb_dcint (255);
    cgmb_dcint (255);

    cgmb_flush_cmd (final_flush);
}



/* VDC integer precision */

static void cgmb_vdcintprec (void)
{
    cgmb_start_cmd (3, (int) vdcIntPrec);

    cgmb_sint (16);

    cgmb_flush_cmd (final_flush);
}



/* Clip rectangle */

static void cgmb_cliprect (int *int_coords)
{
    int i;

    cgmb_start_cmd (3, (int) ClipRect);

    for (i = 0; i < 4; ++i) 
	cgmb_vint (int_coords[i]); 

    cgmb_flush_cmd (final_flush);
}



/* Clip indicator */

static void cgmb_clipindic (int clip_ind)
{
    cgmb_start_cmd (3, (int) ClipIndic);

    cgmb_eint (clip_ind);

    cgmb_flush_cmd (final_flush);
}



/* Polyline */

static void cgmb_pline (int no_pairs, int *x1_ptr, int *y1_ptr)
{
    int i;

    cgmb_start_cmd (4, (int) PolyLine);

    for (i = 0; i < no_pairs; ++i) 
	{
	cgmb_vint (x1_ptr[i]);
	cgmb_vint (y1_ptr[i]);
	}

    cgmb_flush_cmd (final_flush);
}



/* Polymarker */

static void cgmb_pmarker (int no_pairs, int *x1_ptr, int *y1_ptr)
{
    int i;

    cgmb_start_cmd (4, (int) PolyMarker);

    for (i = 0; i < no_pairs; ++i) 
	{
	cgmb_vint (x1_ptr[i]);
	cgmb_vint (y1_ptr[i]);
	}

    cgmb_flush_cmd (final_flush);
}



/* Text */

static void cgmb_text (int x, int y, int final, char *buffer)
{
    cgmb_start_cmd (4, (int) Text);

    cgmb_vint (x);
    cgmb_vint (y);

    cgmb_eint (final);
    cgmb_string (buffer, strlen (buffer));

    cgmb_flush_cmd (final_flush);
}



/* Polygon */

static void cgmb_pgon (int no_pairs, int *x1_ptr, int *y1_ptr)
{
    int i;

    cgmb_start_cmd (4, (int) Polygon);

    for (i = 0; i < no_pairs; ++i) 
	{
	cgmb_vint (x1_ptr[i]);
	cgmb_vint (y1_ptr[i]);
	}

    cgmb_flush_cmd (final_flush);
}



/* Cell array */

static void cgmb_carray (int xmin, int xmax, int ymin, int ymax, int dx, int dy,
    int dimx, int *array)
{
    int ix, iy, c;

    cgmb_start_cmd (4, (int) Cell_Array);

    cgmb_vint (xmin);
    cgmb_vint (ymin);
    cgmb_vint (xmax);
    cgmb_vint (ymax);
    cgmb_vint (xmax);
    cgmb_vint (ymin);

    cgmb_sint (dx);
    cgmb_sint (dy);
    cgmb_sint (cprec);
    cgmb_eint (1);

    for (iy = 0; iy < dy; iy++)
	{
	for (ix = 0; ix < dx; ix++)
	    {
            c = array[dimx * iy + ix];
            c = Color8Bit(c);
	    cgmb_out_bc (c);
	    }

	if (odd (dx))
	    cgmb_out_bc (0);
	}

    cgmb_flush_cmd (final_flush);
}



/* Line type */

static void cgmb_ltype (int line_type)
{
    cgmb_start_cmd (5, (int) LType);

    cgmb_xint (line_type);

    cgmb_flush_cmd(final_flush);
}




/* Line width */

#ifdef linux
static void cgmb_lwidth(rmul)
    float rmul;
#else
static void cgmb_lwidth(float rmul)
#endif
{
    cgmb_start_cmd (5, (int) LWidth);

    cgmb_fixed (rmul); 

    cgmb_flush_cmd (final_flush);
}



/* Line colour */

static void cgmb_lcolour (int index)
{
    cgmb_start_cmd (5, (int) LColour);

    cgmb_cxint (index); 

    cgmb_flush_cmd (final_flush);
}



/* Marker type */

static void cgmb_mtype (int marker)
{
    cgmb_start_cmd (5, (int) MType);

    cgmb_xint (marker); 

    cgmb_flush_cmd (final_flush);
}



/* Marker size */

#ifdef linux
static void cgmb_msize(rmul)
    float rmul;
#else
static void cgmb_msize(float rmul)
#endif
{
    cgmb_start_cmd (5, (int) MSize);

    cgmb_fixed (rmul); 

    cgmb_flush_cmd (final_flush);
}



/* Marker colour */

static void cgmb_mcolour (int index)
{

    cgmb_start_cmd (5, (int) MColour);

    cgmb_cxint (index); 

    cgmb_flush_cmd (final_flush);
}



/* Text font index */

static void cgmb_tfindex (int index)
{
    cgmb_start_cmd (5, (int) TFIndex);

    cgmb_xint (index);

    cgmb_flush_cmd (final_flush);
}



/* Text precision */

static void cgmb_tprec (int precision)
{
    cgmb_start_cmd (5, (int) TPrec);

    cgmb_eint (precision);

    cgmb_flush_cmd (final_flush);
}



/* Character expansion factor */

#ifdef linux
static void cgmb_cexpfac(factor)
    float factor;
#else
static void cgmb_cexpfac(float factor)
#endif
{
    cgmb_start_cmd (5, (int) CExpFac);

    cgmb_fixed (factor);

    cgmb_flush_cmd (final_flush);
}



/* Character space */

#ifdef linux
static void cgmb_cspace(space)
    float space;
#else
static void cgmb_cspace(float space)
#endif
{
    cgmb_start_cmd (5, (int) CSpace);

    cgmb_fixed (space);

    cgmb_flush_cmd (final_flush);
}



/* Text colour */

static void cgmb_tcolour (int index)
{
    cgmb_start_cmd(5, (int) TColour);

    cgmb_cxint(index); 

    cgmb_flush_cmd(final_flush);
}




/* Character height */

static void cgmb_cheight (int height)
{
    cgmb_start_cmd (5, (int) CHeight);

    cgmb_vint (height);

    cgmb_flush_cmd (final_flush);
}



/* Character orientation */

static void cgmb_corient (int x_up, int y_up, int x_base, int y_base)
{
    cgmb_start_cmd (5, (int) COrient);

    cgmb_vint (x_up);
    cgmb_vint (y_up);
    cgmb_vint (x_base);
    cgmb_vint (y_base);

    cgmb_flush_cmd (final_flush);
}




/* Text path */

static void cgmb_tpath (int new_path)
{
    cgmb_start_cmd (5, (int) TPath);

    cgmb_eint (new_path);

    cgmb_flush_cmd (final_flush);
}



/* Text alignment */

static void cgmb_talign (int hor, int ver) 
{
    cgmb_start_cmd (5, (int) TAlign);

    cgmb_eint (hor);
    cgmb_eint (ver);
    cgmb_fixed (0.);
    cgmb_fixed (0.);

    cgmb_flush_cmd (final_flush);
}



/* Interior style */

static void cgmb_intstyle (int style)
{
    cgmb_start_cmd (5, (int) IntStyle);

    cgmb_eint (style);

    cgmb_flush_cmd (final_flush);
}



/* Fill colour */

static void cgmb_fillcolour (int index)
{
    cgmb_start_cmd (5, (int) FillColour);

    cgmb_cxint (index); 

    cgmb_flush_cmd(final_flush);
}




/* Hatch index */

static void cgmb_hindex (int new_index)
{
    cgmb_start_cmd (5, (int) HatchIndex);

    cgmb_xint (new_index);

    cgmb_flush_cmd (final_flush);
}


/* Pattern index */

static void cgmb_pindex (int new_index)
{
    cgmb_start_cmd (5, (int) PatIndex);

    cgmb_xint (new_index);

    cgmb_flush_cmd (final_flush);
}



/* Colour table */

static void cgmb_coltab (int beg_index, int no_entries, float *ctab)
{
    int i, j;

    cgmb_start_cmd (5, (int) ColTab);
    cgmb_cxint (beg_index);

    for (i = beg_index; i < (beg_index + no_entries); ++i)
	{
	for (j = 0; j < 3; ++j) 
	    {
	    cgmb_dcint ((int)(ctab[(i - beg_index) * 3 + j] * (max_colors-1)));
	    }
	}

    cgmb_flush_cmd (final_flush);
}



static void init_color_table (void)
{
    int i, j;

    for (i = 0; i < max_colorxs; i++)
        {
        j = i;
        GQRGB (&j, &p->color_t [i*3], &p->color_t [i*3+1], &p->color_t [i*3+2]);
        }
}



static void setup_colors (void)
{
    int i;

    for (i = 0; i < max_colorxs; i++)
	p->cgm[coltab] (i, 1, &p->color_t[3*i]);
}



static void set_xform (unsigned init)
{
    int errind, tnr;
    float vp_new[4], wn_new[4];
    float clprt[4];
    static int clip_old;
    int clip_new, clip_rect[4];
    int i;
    unsigned update = FALSE;

    if (init)
	{
	GQCNTN (&errind, &tnr);
	GQNT (&tnr, &errind, p->wn, p->vp);
	GQCLIP (&errind, &clip_old, clprt);
	}

    GQCNTN (&errind, &tnr);
    GQNT (&tnr, &errind, wn_new, vp_new);
    GQCLIP (&errind, &clip_new, clprt);

    for (i = 0; i < 4; i++)
	{
	if (vp_new[i] != p->vp[i])
	    {
	    p->vp[i] = vp_new[i];
	    update = TRUE;
	    }
	if (wn_new[i] != p->wn[i])
	    {
	    p->wn[i] = wn_new[i];
	    update = TRUE;
	    }
	}

    if (init || update || (clip_old != clip_new))
	{
	p->xform.a = (vp_new[1] - vp_new[0]) / (wn_new[1] - wn_new[0]);
	p->xform.b = vp_new[0] - wn_new[0] * p->xform.a;
	p->xform.c = (vp_new[3] - vp_new[2]) / (wn_new[3] - wn_new[2]);
	p->xform.d = vp_new[2] - wn_new[2] * p->xform.c;

	if (init)
	    {
	    if (clip_new)
		{
		clip_rect[0] = (int)(vp_new[0] * max_coord);
		clip_rect[1] = (int)(vp_new[2] * max_coord);
		clip_rect[2] = (int)(vp_new[1] * max_coord);
		clip_rect[3] = (int)(vp_new[3] * max_coord);

		p->cgm[cliprect] (clip_rect);
		p->cgm[clipindic] (TRUE);
		}
	    else
		{
		p->cgm[clipindic] (FALSE);
		}                   
	    clip_old = clip_new;
	    }
	else
	    {
	    if ((clip_old != clip_new) || update)
		{
		if (clip_new)
		    {
		    clip_rect[0] = (int)(vp_new[0] * max_coord);
		    clip_rect[1] = (int)(vp_new[2] * max_coord);
		    clip_rect[2] = (int)(vp_new[1] * max_coord);
		    clip_rect[3] = (int)(vp_new[3] * max_coord);

		    p->cgm[cliprect] (clip_rect);
		    p->cgm[clipindic] (TRUE);
		    }
		else
		    {
		    p->cgm[clipindic] (FALSE);
		    }               
		clip_old = clip_new;
		}
	    }
	}
}



static void output_points (void (*output_func)(int, int *, int *), int n_points,
    float *x, float *y)
{
    int i;
    static int x_buffer[max_pbuffer], y_buffer[max_pbuffer];
    int *d_x_buffer, *d_y_buffer;

    set_xform (FALSE);

    if (n_points > max_pbuffer)
	{
	d_x_buffer = (int *)malloc (sizeof(float) * n_points);
	d_y_buffer = (int *)malloc (sizeof(float) * n_points);

	for (i = 0; i < n_points; i++)
	    {
            WC_to_VDC(x[i], y[i], &d_x_buffer[i], &d_y_buffer[i]);
	    }

	output_func (n_points, d_x_buffer, d_y_buffer);

	free (d_y_buffer);
	free (d_x_buffer);
	}
    else
	{
	for (i = 0; i < n_points; i++)
	    {
            WC_to_VDC(x[i], y[i], &x_buffer[i], &y_buffer[i]);
	    }

	output_func (n_points, x_buffer, y_buffer);
	}
}



static void setup_polyline_attributes (unsigned init)
{
    line_attributes newpline;
    int errind;

    if (init)
	{
	p->pline.type = 1;
	p->pline.width = 1.0;
	p->pline.color = 1;
	}
    else
	{
	GQLN (&errind, &newpline.type);
	GQLWSC (&errind, &newpline.width);
	GQPLCI (&errind, &newpline.color);

	if (p->encode == cgm_grafkit)
	    {
	    if (newpline.type < 0)
		newpline.type = max_std_linetype - newpline.type;
	    }

	if (newpline.type != p->pline.type)
	    {
	    p->cgm[ltype] (newpline.type);
	    p->pline.type = newpline.type;
	    }
	if (newpline.width != p->pline.width)
	    {
	    p->cgm[lwidth] (newpline.width);
	    p->pline.width = newpline.width;
	    }
	if (newpline.color != p->pline.color)
	    {
	    p->cgm[lcolour] (newpline.color);
	    p->pline.color = newpline.color;
	    }
	}
}


    
static void setup_polymarker_attributes (unsigned init)
{
    marker_attributes newpmark;
    int errind;

    if (init)
	{
	p->pmark.type = 3;
	p->pmark.width = 1.0;
	p->pmark.color = 1;
	}
    else
	{
	GQMK (&errind, &newpmark.type);
	GQMKSC (&errind, &newpmark.width);
	GQPMCI (&errind, &newpmark.color);

	if (p->encode == cgm_grafkit)
	    {
	    if (newpmark.type < 0)
		newpmark.type = max_std_markertype - newpmark.type;
	    if (newpmark.type > 5)
		newpmark.type = 3;
	    }

	if (newpmark.type != p->pmark.type)
	    {
	    p->cgm[mtype] (newpmark.type);
	    p->pmark.type = newpmark.type;
	    }
	if (newpmark.width != p->pmark.width)
	    {
	    p->cgm[msize] (newpmark.width);
	    p->pmark.width = newpmark.width;
	    }
	if (newpmark.color != p->pmark.color)
	    {
	    p->cgm[mcolour] (newpmark.color);
	    p->pmark.color = newpmark.color;
	    }
	}
}



static void setup_text_attributes (unsigned init)
{
    text_attributes newtext;
    int errind;
    float upx, upy, norm;

    if (init)
	{
	p->text.font = 1;
	p->text.prec = 0;
	p->text.expfac = 1.0;
	p->text.spacing = 0.0;
	p->text.color = 1;
	p->text.height = 0.01;
	p->text.upx = 0; 
	p->text.upy = max_coord;
	p->text.path = 0;
	p->text.halign = 0;
	p->text.valign = 0;
	}
    else
	{
	GQTXFP (&errind, &newtext.font, &newtext.prec);
	GQCHXP (&errind, &newtext.expfac);
	GQCHSP (&errind, &newtext.spacing);
	GQTXCI (&errind, &newtext.color);
        GSCT ();
        GCHH (&newtext.height);
	GQCHUP (&errind, &upx, &upy);
        upx *= p->xform.a;
        upy *= p->xform.c;
	GCST (&upx, &upy);
        norm = fabs(upx) > fabs(upy) ? fabs(upx) : fabs(upy);
        newtext.upx = (int)(upx / norm * max_coord);
        newtext.upy = (int)(upy / norm * max_coord);
	GQTXP (&errind, &newtext.path);
	GQTXAL (&errind, &newtext.halign, &newtext.valign);

	if (p->encode == cgm_grafkit)
	    {
	    if (newtext.font < 0)
		newtext.font = max_std_textfont - newtext.font;
	    newtext.prec = 2;
	    }

	if (newtext.font != p->text.font)
	    {
	    p->cgm[tfindex] (newtext.font);
	    p->text.font = newtext.font;
	    }
	if (newtext.prec != p->text.prec)
	    {
	    p->cgm[tprec] (newtext.prec);
	    p->text.prec = newtext.prec;
	    }
	if (newtext.expfac != p->text.expfac)
	    {
	    p->cgm[cexpfac] (newtext.expfac);
	    p->text.expfac = newtext.expfac;
	    }
	if (newtext.spacing != p->text.spacing)
	    {
	    p->cgm[cspace] (newtext.spacing);
	    p->text.spacing = newtext.spacing;
	    }
	if (newtext.color != p->text.color)
	    {
	    p->cgm[tcolour] (newtext.color);
	    p->text.color = newtext.color;
	    }
	if (newtext.height != p->text.height)
	    {
	    p->cgm[cheight] ((int)(newtext.height * max_coord));
	    p->text.height = newtext.height;
	    }
	if ((newtext.upx != p->text.upx) || (newtext.upy != p->text.upy))
	    {
	    p->cgm[corient] (newtext.upx, newtext.upy, newtext.upy,
                -newtext.upx);
	    p->text.upx = newtext.upx;
	    p->text.upy = newtext.upy;
	    }
	if (newtext.path != p->text.path)
	    {
	    p->cgm[tpath] (newtext.path);
	    p->text.path = newtext.path;
	    }
	if ((newtext.halign != p->text.halign) || (newtext.valign !=
	    p->text.valign))
	    {
	    p->cgm[talign] (newtext.halign, newtext.valign);
	    p->text.halign = newtext.halign;
	    p->text.valign = newtext.valign;
	    }
	}
}

    

static void setup_fill_attributes (unsigned init)
{
    fill_attributes newfill;
    int errind;

    if (init)
	{
	p->fill.intstyle = 0;
	p->fill.color = 1;
	p->fill.pattern_index = 1;
	p->fill.hatch_index = 1;
	}
    else
	{
	GQFAIS (&errind, &newfill.intstyle);
	GQFACI (&errind, &newfill.color);
	GQFASI (&errind, &newfill.pattern_index);
	GQFASI (&errind, &newfill.hatch_index);
	
	if (newfill.intstyle != p->fill.intstyle)
	    {
	    p->cgm[intstyle] (newfill.intstyle);
	    p->fill.intstyle = newfill.intstyle;
	    }
	if (newfill.color != p->fill.color)
	    {
	    p->cgm[fillcolour] (newfill.color);
	    p->fill.color = newfill.color;
	    }
	if (newfill.pattern_index != p->fill.pattern_index)
	    {
	    p->cgm[pindex] (newfill.pattern_index);
	    p->fill.pattern_index = newfill.pattern_index;
	    }
	if (newfill.hatch_index != p->fill.hatch_index)
	    {
	    p->cgm[hindex] (newfill.hatch_index);
	    p->fill.hatch_index = newfill.hatch_index;
	    }
	}
}



static char *local_time (void)
{
    struct tm *time_structure;
    time_t time_val;
    static char *weekday[7] = {
	"Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday",
	"Saturday"};
    static char *month[12] = {
	"January", "February", "March", "April", "May", "June", "July", 
	"August", "September", "October", "November", "Dezember"};
    static char time_string[81];

    time (&time_val);
    time_structure = localtime (&time_val);

    sprintf (time_string, "%s, %s %d, 19%d %d:%02d:%02d", 
	weekday[time_structure->tm_wday], month[time_structure->tm_mon],
	time_structure->tm_mday, time_structure->tm_year,
	time_structure->tm_hour, time_structure->tm_min,
	time_structure->tm_sec);

    return (time_string);
}
    


static void setup_clear_text_context (void)
{
    p->cgm[begin] = CGM_FUNC cgmt_begin;
    p->cgm[end] = CGM_FUNC cgmt_end;
    p->cgm[bp] = CGM_FUNC cgmt_bp;
    p->cgm[bpage] = CGM_FUNC cgmt_bpage;
    p->cgm[epage] = CGM_FUNC cgmt_epage;
    p->cgm[mfversion] = CGM_FUNC cgmt_mfversion;
    p->cgm[mfdescrip] = CGM_FUNC cgmt_mfdescrip;
    p->cgm[vdctype] = CGM_FUNC cgmt_vdctype;
    p->cgm[intprec] = CGM_FUNC cgmt_intprec;
    p->cgm[realprec] = CGM_FUNC cgmt_realprec;
    p->cgm[indexprec] = CGM_FUNC cgmt_indexprec;
    p->cgm[colprec] = CGM_FUNC cgmt_colprec;
    p->cgm[cindprec] = CGM_FUNC cgmt_cindprec;
    p->cgm[cvextent] = CGM_FUNC cgmt_cvextent;
    p->cgm[maxcind] = CGM_FUNC cgmt_maxcind;
    p->cgm[mfellist] = CGM_FUNC cgmt_mfellist;
    p->cgm[fontlist] = CGM_FUNC cgmt_fontlist;
    p->cgm[cannounce] = CGM_FUNC cgmt_cannounce;
    p->cgm[scalmode] = CGM_FUNC cgmt_scalmode;
    p->cgm[colselmode] = CGM_FUNC cgmt_colselmode;
    p->cgm[lwsmode] = CGM_FUNC cgmt_lwsmode;
    p->cgm[msmode] = CGM_FUNC cgmt_msmode;
    p->cgm[vdcextent] = CGM_FUNC cgmt_vdcextent;
    p->cgm[backcol] = CGM_FUNC cgmt_backcol;
    p->cgm[vdcintprec] = CGM_FUNC cgmt_vdcintprec;
    p->cgm[cliprect] = CGM_FUNC cgmt_cliprect;
    p->cgm[clipindic] = CGM_FUNC cgmt_clipindic;
    p->cgm[pline] = CGM_FUNC cgmt_pline;
    p->cgm[pmarker] = CGM_FUNC cgmt_pmarker;
    p->cgm[text] = CGM_FUNC cgmt_text;
    p->cgm[pgon] = CGM_FUNC cgmt_pgon;
    p->cgm[ltype] = CGM_FUNC cgmt_ltype;
    p->cgm[lwidth] = CGM_FUNC cgmt_lwidth;
    p->cgm[lcolour] = CGM_FUNC cgmt_lcolour;
    p->cgm[mtype] = CGM_FUNC cgmt_mtype;
    p->cgm[msize] = CGM_FUNC cgmt_msize;
    p->cgm[mcolour] = CGM_FUNC cgmt_mcolour;
    p->cgm[tfindex] = CGM_FUNC cgmt_tfindex;
    p->cgm[tprec] = CGM_FUNC cgmt_tprec;
    p->cgm[cexpfac] = CGM_FUNC cgmt_cexpfac;
    p->cgm[cspace] = CGM_FUNC cgmt_cspace;
    p->cgm[tcolour] = CGM_FUNC cgmt_tcolour;
    p->cgm[cheight] = CGM_FUNC cgmt_cheight;
    p->cgm[corient] = CGM_FUNC cgmt_corient;
    p->cgm[tpath] = CGM_FUNC cgmt_tpath;
    p->cgm[talign] = CGM_FUNC cgmt_talign;
    p->cgm[intstyle] = CGM_FUNC cgmt_intstyle;
    p->cgm[fillcolour] = CGM_FUNC cgmt_fillcolour;
    p->cgm[hindex] = CGM_FUNC cgmt_hindex;
    p->cgm[pindex] = CGM_FUNC cgmt_pindex;
    p->cgm[coltab] = CGM_FUNC cgmt_coltab;
    p->cgm[carray] = CGM_FUNC cgmt_carray;

    p->buffer_ind = 0;
    p->buffer[0] = '\0';
}



static void setup_binary_context (void)
{
    p->cgm[begin] = CGM_FUNC cgmb_begin;
    p->cgm[end] = CGM_FUNC cgmb_end;
    p->cgm[bp] = CGM_FUNC cgmb_bp;
    p->cgm[bpage] = CGM_FUNC cgmb_bpage;
    p->cgm[epage] = CGM_FUNC cgmb_epage;
    p->cgm[mfversion] = CGM_FUNC cgmb_mfversion;
    p->cgm[mfdescrip] = CGM_FUNC cgmb_mfdescrip;
    p->cgm[vdctype] = CGM_FUNC cgmb_vdctype;
    p->cgm[intprec] = CGM_FUNC cgmb_intprec;
    p->cgm[realprec] = CGM_FUNC cgmb_realprec;
    p->cgm[indexprec] = CGM_FUNC cgmb_indexprec;
    p->cgm[colprec] = CGM_FUNC cgmb_colprec;
    p->cgm[cindprec] = CGM_FUNC cgmb_cindprec;
    p->cgm[cvextent] = CGM_FUNC cgmb_cvextent;
    p->cgm[maxcind] = CGM_FUNC cgmb_maxcind;
    p->cgm[mfellist] = CGM_FUNC cgmb_mfellist;
    p->cgm[fontlist] = CGM_FUNC cgmb_fontlist;
    p->cgm[cannounce] = CGM_FUNC cgmb_cannounce;
    p->cgm[scalmode] = CGM_FUNC cgmb_scalmode;
    p->cgm[colselmode] = CGM_FUNC cgmb_colselmode;
    p->cgm[lwsmode] = CGM_FUNC cgmb_lwsmode;
    p->cgm[msmode] = CGM_FUNC cgmb_msmode;
    p->cgm[vdcextent] = CGM_FUNC cgmb_vdcextent;
    p->cgm[backcol] = CGM_FUNC cgmb_backcol;
    p->cgm[vdcintprec] = CGM_FUNC cgmb_vdcintprec;
    p->cgm[cliprect] = CGM_FUNC cgmb_cliprect;
    p->cgm[clipindic] = CGM_FUNC cgmb_clipindic;
    p->cgm[pline] = CGM_FUNC cgmb_pline;
    p->cgm[pmarker] = CGM_FUNC cgmb_pmarker;
    p->cgm[text] = CGM_FUNC cgmb_text;
    p->cgm[pgon] = CGM_FUNC cgmb_pgon;
    p->cgm[ltype] = CGM_FUNC cgmb_ltype;
    p->cgm[lwidth] = CGM_FUNC cgmb_lwidth;
    p->cgm[lcolour] = CGM_FUNC cgmb_lcolour;
    p->cgm[mtype] = CGM_FUNC cgmb_mtype;
    p->cgm[msize] = CGM_FUNC cgmb_msize;
    p->cgm[mcolour] = CGM_FUNC cgmb_mcolour;
    p->cgm[tfindex] = CGM_FUNC cgmb_tfindex;
    p->cgm[tprec] = CGM_FUNC cgmb_tprec;
    p->cgm[cexpfac] = CGM_FUNC cgmb_cexpfac;
    p->cgm[cspace] = CGM_FUNC cgmb_cspace;
    p->cgm[tcolour] = CGM_FUNC cgmb_tcolour;
    p->cgm[cheight] = CGM_FUNC cgmb_cheight;
    p->cgm[corient] = CGM_FUNC cgmb_corient;
    p->cgm[tpath] = CGM_FUNC cgmb_tpath;
    p->cgm[talign] = CGM_FUNC cgmb_talign;
    p->cgm[intstyle] = CGM_FUNC cgmb_intstyle;
    p->cgm[fillcolour] = CGM_FUNC cgmb_fillcolour;
    p->cgm[hindex] = CGM_FUNC cgmb_hindex;
    p->cgm[pindex] = CGM_FUNC cgmb_pindex;
    p->cgm[coltab] = CGM_FUNC cgmb_coltab;
    p->cgm[carray] = CGM_FUNC cgmb_carray;

    p->buffer_ind = 0;
    p->buffer[0] = '\0';

    p->bfr_index = 0;
}



static void cgm_begin_page (void)
{
    p->cgm[bp] (local_time ());

    if (p->encode != cgm_grafkit)
        p->cgm[scalmode] ();

    p->cgm[colselmode] ();

    if (p->encode != cgm_grafkit)
        {
        p->cgm[lwsmode] ();
        p->cgm[msmode] ();
        }

    p->cgm[vdcextent] ();
    p->cgm[backcol] ();

    p->cgm[bpage] ();
    p->cgm[vdcintprec] ();

    setup_colors ();

    set_xform (TRUE);

    setup_polyline_attributes (TRUE);
    setup_polymarker_attributes (TRUE);
    setup_text_attributes (TRUE);
    setup_fill_attributes (TRUE);

    p->begin_page = FALSE;
}



void STDCALL GKDCGM (
    int *fctid, int *dx, int *dy, int *dimx, int *ia, int *lr1, float *r1,
    int *lr2, float *r2, int *lc, CHARARG(chars), cgm_context **context)
{
    char *env, *buffer;

    p = *context;

    switch (*fctid) {

      case 2:
/*
 *  Open workstation
 */
	p = (cgm_context *) malloc (sizeof(cgm_context));

	p->conid = ia[1];

	if ((ia[2] == GCGMO) || (ia[2] == (GCGMO | GCGMB)))
	    {
	    p->encode = cgm_binary;
            setup_binary_context ();
	    }
	else if (ia[2] == 8 || ia[2] == (GCGMO | GCGMT))
	    {
	    p->encode = cgm_clear_text;
            setup_clear_text_context ();
	    }
	else if (ia[2] == (GCGMO | GCGMC))
	    {
	    gks_fprintf (stderr, "GKS: unsupported CGM encoding (character)\n");
	    ia[0] = ia[1] = 0;
            return;
	    }
	else if (ia[2] == (GCGMO | GCGMGK))
	    {
	    p->encode = cgm_grafkit;
	    setup_binary_context ();
	    }
	else
	    {
	    gks_fprintf (stderr, "GKS: invalid bit mask (%x)\n", ia[2]);
	    ia[0] = ia[1] = 0;
            return;
	    }

	buffer = "GLI GKS, Copyright @ 1986-1995, Josef Heinen, Jochen Werner";
        if (env = (char *) getenv("GLI_GKS"))
            if (strncmp(env, "GKSGRAL", 7) == 0)
                buffer = "GKSGRAL CGM produced by GLI GKS";
        
	if (((char *) getenv("GLI_GKS_SCALE_MODE_METRIC")) != NULL)
	    p->mm = 0.19685/max_coord*1000;
	else
	    p->mm = 0;

	p->cgm[begin] (buffer);
	p->cgm[mfversion] ();
	p->cgm[mfdescrip] ();

	if (p->encode != cgm_grafkit)
	    {
	    p->cgm[vdctype] ();
	    p->cgm[intprec] ();
#if 0
	    p->cgm[realprec] (); /* causes problems with RALCGM */
#endif
	    p->cgm[indexprec] ();
	    p->cgm[colprec] ();
	    p->cgm[cindprec] ();
	    p->cgm[maxcind] ();
	    p->cgm[cvextent] ();
	    }

	p->cgm[mfellist] ();
	p->cgm[fontlist] ();

	if (p->encode != cgm_grafkit)
	    p->cgm[cannounce] ();

        init_color_table ();

        p->xext = p->yext = max_coord;

	p->begin_page = TRUE;
	p->active = FALSE;

	*context = p;
	if (sizeof(char *) > sizeof(int))
	    {
            long *la = (long *)ia;
            *la = (long)p;
            }
	else
	    *ia = (long)p;
	break;

      case 3:
/*
 *  Close workstation
 */
	p->cgm[epage] ();
	p->cgm[end] ();

	free (p);
	break;

      case 4:
/*
 *  Activate workstation
 */
	p->active = TRUE;
	break;

      case 5:
/*
 *  Deactivate workstation
 */
	p->active = FALSE;
	break;

      case 6:
/*
 *  Clear workstation
 */
	if (!p->begin_page)
	    {
	    p->cgm[epage] ();
	    p->begin_page = TRUE;
	    }
	break;

      case 12:
/*
 *  Polyline
 */
	if (p->active)
	    {
	    if (p->begin_page)
                cgm_begin_page ();

	    setup_polyline_attributes (FALSE);
	    output_points ((void (*)(int, int *, int *))p->cgm[pline],
		ia[0], r1, r2);
	    }
	break;

      case 13:
/*
 *  Polymarker
 */
	if (p->active)
	    {
	    if (p->begin_page)
                cgm_begin_page ();

	    setup_polymarker_attributes (FALSE);
	    output_points ((void (*)(int, int *, int *))p->cgm[pmarker],
		ia[0], r1, r2);
	    }
	break;

      case 14:
/*
 *  Text
 */
	if (p->active)
	    {
            int x, y;

	    if (p->begin_page)
                cgm_begin_page ();

	    set_xform (FALSE);
	    setup_text_attributes (FALSE);

            WC_to_VDC(r1[0], r2[0], &x, &y);
#ifdef VMS
	    p->cgm[text] (x, y, TRUE, chars->dsc$a_pointer);
#else
	    p->cgm[text] (x, y, TRUE, chars);
#endif
	    }
	break;

      case 15:
/*
 *  Fill Area
 */
	if (p->active)
	    {
	    if (p->begin_page)
                cgm_begin_page ();

	    setup_fill_attributes (FALSE);
	    output_points ((void (*)(int, int *, int *))p->cgm[pgon],
		ia[0], r1, r2);
	    }
	break;

      case 16:
/*
 *  Cell Array
 */
	if (p->active)
	    {
            int xmin, xmax, ymin, ymax;

	    if (p->begin_page)
                cgm_begin_page ();

	    set_xform (FALSE);

            WC_to_VDC(r1[0], r2[0], &xmin, &ymin);
            WC_to_VDC(r1[1], r2[1], &xmax, &ymax);

	    p->cgm[carray] (xmin, xmax, ymin, ymax, *dx, *dy, *dimx, ia);
	    }
	break;

      case 48:
/*
 *  Set color representation
 */
	if (p->begin_page)
            {
            p->color_t[ia[1]*3] = r1[0];
            p->color_t[ia[1]*3+1] = r1[1];
            p->color_t[ia[1]*3+2] = r1[2];
            }
	break;

      case 54:
/*
 *  Set workstation window
 */
	if (p->begin_page)
            {
            p->xext = (int)(max_coord * (r1[1]-r1[0]));
            p->yext = (int)(max_coord * (r2[1]-r2[0]));
            }
        break;

      case 55:
/*
 *  Set workstation viewport
 */
	if (p->begin_page)
	    {
	    if (p->mm > 0)
		p->mm = (r1[1]-r1[0])/max_coord*1000;
	    }
        break;

      }
}
