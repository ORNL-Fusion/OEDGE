/*
 * Copyright @ 1999  Josef Heinen
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
 *      GLI GKS V4.5
 *
 * ABSTRACT:
 *
 *      This module contains a logical device driver for Adobe's
 *      Portable Document Format (PDF)
 *
 * AUTHOR:
 *
 *      Josef Heinen
 *
 * VERSION:
 *
 *      V1.0
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>

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

#include "gksdefs.h"

#ifdef ZLIB
#include <zlib.h>
#else
typedef unsigned char Byte;
typedef unsigned long uLong;
#endif

#define MAX_TNR 9
#define MAX_COLOR 980
#define MAX_FONT 43
#define MAX_BITMAP 120
#define HATCH_STYLE 108

#define MEMORY_INCREMENT 32768
#define MAX_OBJECTS 10000
#define MAX_PAGES 100
#define NO_OF_BUFS 10

#define Pi 3.141592654

#ifndef min
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#define WC_to_NDC(xw, yw, tnr, xn, yn) \
  xn = a[tnr] * (xw) + b[tnr]; \
  yn = c[tnr] * (yw) + d[tnr]

#define WC_to_NDC_rel(xw, yw, tnr, xn, yn) \
  xn = a[tnr] * (xw); \
  yn = c[tnr] * (yw)

#define NDC_to_DC(xn, yn, xd, yd) \
  xd = p->a * (xn) + p->b; \
  yd = p->c * (yn) + p->d

#define DC_to_NDC(xd, yd, xn, yn) \
  xn = ((xd) - p->b) / p->a; \
  yn = ((yd) - p->d) / p->c

#define CharXform(phi, xrel, yrel, x, y) \
  x = cos(phi) * (xrel) - sin(phi) * (yrel); \
  y = sin(phi) * (xrel) + cos(phi) * (yrel);

#define nint(a) ((int)(a + 0.5))

#define Color8Bit(c) \
  c >= 588 ? 80 + (c - 588) / 56 * 12 + nint((c - 588) % 56 * 11.0 / 56.0) : \
  c >= 257 ? 8 + nint((c - 257) / 330.0 * (72 - 1)) : c

#define pdf_obj(p, id) \
  p->byte_offset[id] = p->stream->length; \
  pdf_printf(p->stream, "%ld 0 obj\n", id);

#define pdf_endobj(p)		pdf_printf(p->stream, "endobj\n")
#define pdf_dict(p)		pdf_printf(p->stream, "<<\n")
#define pdf_enddict(p)		pdf_printf(p->stream, ">>\n")
#define pdf_stream(p)		pdf_printf(p->stream, "stream\n")
#define pdf_endstream(p)	pdf_printf(p->stream, "endstream\n")

#define pdf_save(p)		pdf_printf(p->content, "q\n")
#define pdf_restore(p)		pdf_printf(p->content, "Q\n")
#define pdf_clip(p)		pdf_printf(p->content, "W\n")
#define pdf_moveto(p, x, y)	pdf_printf(p->content, "%.2f %.2f m\n", x, y)
#define pdf_lineto(p, x, y)	pdf_printf(p->content, "%.2f %.2f l\n", x, y)
#define pdf_stroke(p)		pdf_printf(p->content, "S\n")
#define pdf_eofill(p)		pdf_printf(p->content, "f*\n")
#define pdf_point(p, x, y)	pdf_printf(p->content, "%.2f %.2f ", x, y)
#define pdf_curveto(p)		pdf_printf(p->content, "c\n")
#define pdf_setdash(p, dash)	pdf_printf(p->content, "%s 0 d\n", dash)

#define pdf_setlinewidth(p, width) \
  pdf_printf(p->content, "%s w\n", pdf_float(width))

#define pdf_text(p, xorg, yorg, text) \
  pdf_printf(p->content, "BT\n/F%d %d Tf\n%.2f %.2f Td\n(%s) Tj\nET\n", \
  p->font, p->size, xorg, yorg, text)

#define pdf_setrgbcolor(p, red, green, blue) \
  pdf_printf(p->content, "%s %s %s RG\n", \
    pdf_float(red), pdf_float(green), pdf_float(blue))

#define pdf_setfillcolor(p, red, green, blue) \
  pdf_printf(p->content, "%s %s %s rg\n", \
    pdf_float(red), pdf_float(green), pdf_float(blue))

#define PDF ws_state_list

typedef struct PDF_stream_t
  {
    Byte *buffer;
    uLong size, length;
  }
PDF_stream;

typedef struct PDF_page_t
  {
    long object, contents, fonts[MAX_FONT];
    float width, height;
    PDF_stream *stream;
  }
PDF_page;

typedef struct ws_state_list_t
  {
    int state;
    int fd;
    float window[4], viewport[4];
    int empty;
    float width, height;
    float a, b, c, d;
    int stroke;
    float lastx, lasty;
    float red[MAX_COLOR], green[MAX_COLOR], blue[MAX_COLOR];
    int color, fillcolor, ltype, font, size;
    float lwidth, angle;
    PDF_stream *stream;
    long object_number;
    long info, root, outlines, pages;
    long byte_offset[MAX_OBJECTS];
    PDF_page *page[MAX_PAGES];    
    int current_page;
    PDF_stream *content;
    int compress;
  }
ws_state_list;

static
ws_state_list *p;

static
gks_state_list *gksl;

static
float a[MAX_TNR], b[MAX_TNR], c[MAX_TNR], d[MAX_TNR];

static
char *dash_pattern[35] =
  {
    "[4 3 4 3 4 3 4 6]", "[4 3 4 3 4 6]",
    "[4 3 4 6]", "[3 2 3 2 3 2 3 6]",
    "[3 2 3 2 3 6]", "[3 2 3 6]",
    "[3 2 3 2 3 2 3 4]", "[3 2 3 2 3 4]",
    "[3 2 3 4]", "[0 2]", "[0 4]", "[0 8]", "[0 10]",
    "[0 3 0 3 0 6]", "[0 3 0 6]", "[6 3 0 3 0 3 0 3]",
    "[6 3 0 3 0 3]", "[6 3 0 3]", "[12 4 6 4]",
    "[12 4]", "[6 6]", "[6 4]",
    "[0 6 0 6 0 12]", "[0 6 0 12]", "[0 12]", "[12 12]",
    "[14 8 12 8]", "[24 8]", "[12 6 0 6 0 6 0 6]",
    "[12 6 0 6 0 6]", "[]",
    "[]", "[12 8]", "[0 6]", "[12 6 0 6]"
  };

static
char *fonts[MAX_FONT] = 
  {
    "Times-Roman", "Times-Italic", "Times-Bold", "Times-BoldItalic",
    "Helvetica", "Helvetica-Oblique", "Helvetica-Bold",
    "Helvetica-BoldOblique", "Courier", "Courier-Oblique",
    "Courier-Bold", "Courier-BoldOblique", "Symbol",
    "LubalinGraph-Book", "LubalinGraph-BookOblique",
    "LubalinGraph-Demi", "LubalinGraph-DemiOblique",
    "NewCenturySchlbk-Roman", "NewCenturySchlbk-Italic",
    "NewCenturySchlbk-Bold", "NewCenturySchlbk-BoldItalic",
    "AvantGarde-Book", "AvantGarde-BookOblique", "AvantGarde-Demi",
    "AvantGarde-DemiOblique", "Souvenir-Light",
    "Souvenir-LightItalic", "Souvenir-Demi", "Souvenir-DemiItalic",
    "Helvetica-Narrow", "Helvetica-Narrow-Oblique",
    "Helvetica-Narrow-Bold", "Helvetica-Narrow-BoldOblique",
    "Bookman-Light", "Bookman-LightItalic", "Bookman-Demi",
    "Bookman-DemiItalic", "Palatino-Roman", "Palatino-Italic",
    "Palatino-Bold", "Palatino-BoldItalic",
    "ZapfChancery-MediumItalic", "ZapfDingbats"
  };

static
int map[32] =
  {
    22,  9,  5, 14, 18, 26, 13,  1,
    24, 11,  7, 16, 20, 28, 13,  3,
    23, 10,  6, 15, 19, 27, 13,  2,
    25, 12,  8, 17, 21, 29, 13,  4
  };

static
char bitmap[MAX_BITMAP][65];

static
int predef_font[] = { 1, 1, 1, -2, -3, -4 };

static
int predef_prec[] = { 0, 1, 2, 2, 2, 2 };

static
float xfac[4] = { 0, 0, -0.5, -1 };

static
float yfac[6] = { 0, -1.2, -1, -0.5, 0, 0.2 };

#if defined(hpux) && !defined(NAGware)
static int (*line_routine_a)();
static int (*fill_routine_a)();
static int (*text_routine_a)();
#endif

static
void fill_routine(int *n, float *px, float *py, int *tnr);

static char buf_array[NO_OF_BUFS][20];
static int current_buf = 0;

static
char *pdf_float(float f)
{
  char *buf = buf_array[(current_buf++) % NO_OF_BUFS];

  if (fabs(f) < 0.00001)
    return "0";

  sprintf(buf, "%.4g", f);
  if (strchr(buf, 'e'))
    {
      if (fabs(f) < 1)
	sprintf(buf, "%1.5f", f);
      else if (fabs(f) < 1000)
	sprintf(buf, "%1.2f", f);
      else
	sprintf(buf, "%1.0f", f);
    }

  return buf;
}

static
void pdf_memcpy(PDF_stream *p, char *s, size_t n)
{
  if (p->length + n >= p->size)
    {
      while (p->length + n >= p->size)
	p->size += MEMORY_INCREMENT;
      p->buffer = (Byte *)realloc(p->buffer, p->size);
    }

  memcpy(p->buffer + p->length, s, n);
  p->length += n;
}

static
void pdf_printf(PDF_stream *p, char *args,...)
{
  va_list ap;
  char fmt[BUFSIZ], s[BUFSIZ];

  strcpy(fmt, args);

  va_start(ap, args);
  vsprintf(s, fmt, ap);
  va_end(ap);

  pdf_memcpy(p, s, strlen(s));
}

static
void pdf_error(char *args,...)
{
  va_list ap;
  char fmt[BUFSIZ], s[BUFSIZ];

  strcpy(fmt, args);

  va_start(ap, args);
  vsprintf(s, fmt, ap);
  va_end(ap);

  gks_fprintf(stderr, s);
  exit(-1);
}

static
PDF_stream *pdf_alloc_stream(void)
{
  PDF_stream *p;

  p = (PDF_stream *)calloc(1, sizeof(PDF_stream));
  p->buffer = NULL;
  p->size = p->length = 0;

  return p;
}

static
long pdf_alloc_id(PDF *p)
{
  if (p->object_number >= MAX_OBJECTS)
    pdf_error("Too many objects (%ld)", p->object_number);

  return ++(p->object_number);
}

static
void pdf_open(int fd)
{
  p->fd = fd;

  p->stream = pdf_alloc_stream();

  p->object_number = p->current_page = 0;

  p->info = pdf_alloc_id(p);
  p->root = pdf_alloc_id(p);
  p->outlines = pdf_alloc_id(p);
  p->pages = pdf_alloc_id(p);
}

static
void pdf_page(PDF *p, float height, float width)
{
  PDF_page *page;
  int font;

  if (++(p->current_page) >= MAX_PAGES)
    pdf_error("Too many pages in document (%d)\n", p->current_page);

  page = (PDF_page *)calloc(1, sizeof(PDF_page));

  page->object = pdf_alloc_id(p);
  page->contents = pdf_alloc_id(p);
  page->width = width;
  page->height = height;
  page->stream = pdf_alloc_stream();

  p->page[p->current_page] = page;
  p->content = page->stream;

  for (font = 0; font < MAX_FONT; font++)
    page->fonts[font] = 0;
}

static
void pdf_close(PDF *p)
{
  time_t timer;
  struct tm ltime;
  long start_xref;
  int count, object, font;

  pdf_printf(p->stream, "%%PDF-1.%d\n", p->compress ? 2 : 0);
  pdf_printf(p->stream, "%%âãÏÓ\n");

  time(&timer);
  ltime = *localtime(&timer);

  pdf_obj(p, p->info);
  pdf_dict(p);
  pdf_printf(p->stream, "/Creator (GLI GKS)\n");
  pdf_printf(p->stream, "/CreationDate (D:%04d%02d%02d%02d%02d%02d)\n",
    ltime.tm_year + 1900, ltime.tm_mon + 1, ltime.tm_mday, ltime.tm_hour,
    ltime.tm_min, ltime.tm_sec);
  pdf_printf(p->stream, "/Producer (%s)\n", "GLI GKS V4.5 PDF driver");
  pdf_enddict(p);
  pdf_endobj(p);

  pdf_obj(p, p->root);
  pdf_dict(p);
  pdf_printf(p->stream, "/Type /Catalog\n");
  pdf_printf(p->stream, "/Pages %ld 0 R\n", p->pages);
  pdf_printf(p->stream, "/Outlines %ld 0 R\n", p->outlines);
  pdf_enddict(p);
  pdf_endobj(p);

  pdf_obj(p, p->outlines);
  pdf_dict(p);
  pdf_printf(p->stream, "/Type /Outlines\n");
  pdf_printf(p->stream, "/Count 0\n");
  pdf_enddict(p);
  pdf_endobj(p);

  pdf_obj(p, p->pages);
  pdf_dict(p);
  pdf_printf(p->stream, "/Type /Pages\n");
  pdf_printf(p->stream, "/Count %d\n", p->current_page);
  pdf_printf(p->stream, "/Kids [");

  for (count = 1; count <= p->current_page; count++)
    {
      pdf_printf(p->stream, "%ld 0 R", p->page[count]->object);
      if (count < p->current_page)
	pdf_printf(p->stream, count % 6 ? " " : "\n");
    }

  pdf_printf(p->stream, "]\n");
  pdf_enddict(p);
  pdf_endobj(p);

  for (count = 1; count <= p->current_page; count++)
    {
      PDF_page *page = p->page[count];

      pdf_obj(p, page->object);
      pdf_dict(p);
      pdf_printf(p->stream, "/Type /Page\n");
      pdf_printf(p->stream, "/Parent %ld 0 R\n", p->pages);
      pdf_printf(p->stream, "/Resources << /Font <<");
      for (font = 0; font < MAX_FONT; font++)
	{
	  if (page->fonts[font])
	    pdf_printf(p->stream, " /F%d %ld 0 R", font, page->fonts[font]);
	}
      pdf_printf(p->stream, " >> /ProcSet [/PDF /Text /ImageC] >>\n");
      pdf_printf(p->stream, "/MediaBox [0 0 %g %g]\n",
	page->height, page->width);
      pdf_printf(p->stream, "/Contents %ld 0 R\n", page->contents);
      pdf_enddict(p);
      pdf_endobj(p);

      p->content = page->stream;
      pdf_obj(p, page->contents);
      pdf_dict(p);

#ifdef ZLIB
      if (p->compress)
	{
	  Byte *buffer;
	  uLong length;
	  int err;

	  length = p->content->length + 1024;
	  buffer = (Byte *)calloc((int)length, 1);
	  if ((err = compress(buffer, &length, p->content->buffer,
	    p->content->length)) != Z_OK)
	    {
	      pdf_error("Compression failed (err=%d)\n", err);
	    }
	  free(p->content->buffer);

	  p->content->buffer = buffer;
	  p->content->size = p->content->length = length;
	  pdf_printf(p->stream, "/Length %ld\n", p->content->length);
	  pdf_printf(p->stream, "/Filter [/FlateDecode]\n");
	  buffer[p->content->length++] = '\n';
	}
      else
	{
	  pdf_printf(p->stream, "/Length %ld\n", p->content->length);
	}
#else
      pdf_printf(p->stream, "/Length %ld\n", p->content->length);
#endif

      pdf_enddict(p);
      pdf_stream(p);
      pdf_memcpy(p->stream, (char *)p->content->buffer, p->content->length);
      pdf_endstream(p);
      pdf_endobj(p);

      for (font = 0; font < MAX_FONT; font++)
	{
	  if (page->fonts[font])
	    {
	      pdf_obj(p, page->fonts[font]);
	      pdf_dict(p);
	      pdf_printf(p->stream, "/Type /Font\n");
	      pdf_printf(p->stream, "/Subtype /Type1\n");
	      pdf_printf(p->stream, "/Name /F%d\n", font);
	      pdf_printf(p->stream, "/BaseFont /%s\n", fonts[font]);
	      pdf_printf(p->stream, "/Encoding /WinAnsiEncoding\n");
	      pdf_enddict(p);
	      pdf_endobj(p);
	    }
	}
      free(p->content->buffer);
    }

  start_xref = p->stream->length;
  pdf_printf(p->stream, "xref\n");
  pdf_printf(p->stream, "0 %ld\n", p->object_number + 1);
  pdf_printf(p->stream, "0000000000 65535 f \n");
  for (object = 1; object <= p->object_number; object++)
    pdf_printf(p->stream, "%010ld 00000 n \n", p->byte_offset[object]);

  pdf_printf(p->stream, "trailer\n");
  pdf_dict(p);
  pdf_printf(p->stream, "/Size %ld\n", p->object_number + 1);
  pdf_printf(p->stream, "/Root %ld 0 R\n", p->root);
  pdf_printf(p->stream, "/Info %ld 0 R\n", p->info);
  pdf_enddict(p);
  pdf_printf(p->stream, "startxref\n");
  pdf_printf(p->stream, "%ld\n", start_xref);

  pdf_printf(p->stream, "%%%%EOF\n");

#if defined(_WIN32) && !defined(__GNUC__)
  BUFOUT(&p->fd, &p->stream->length, p->stream->buffer,
    (unsigned short)p->stream->length);
#else
  BUFOUT(&p->fd, (int *)&p->stream->length, (char *)p->stream->buffer,
    (unsigned short)p->stream->length);
#endif

  free(p->stream->buffer);
}

static
void pdf_text_ex(PDF *p, float xorg, float yorg, char *text)
{
  float rad, c, s;

  rad = p->angle * Pi / 180;
  c = cos(rad);
  s = sin(rad);

  pdf_printf(p->content,
    "BT\n/F%d %d Tf\n%s %s %s %s %.2f %.2f Tm\n(%s) Tj\nET\n",
    p->font, p->size, pdf_float(c), pdf_float(s), pdf_float(-s), pdf_float(c),
    xorg, yorg, text);
}

static
void set_norm_xform(int tnr, float *wn, float *vp)
{
  a[tnr] = (vp[1] - vp[0]) / (wn[1] - wn[0]);
  b[tnr] = vp[0] - wn[0] * a[tnr];
  c[tnr] = (vp[3] - vp[2]) / (wn[3] - wn[2]);
  d[tnr] = vp[2] - wn[2] * c[tnr];
}

static
void init_norm_xform(void)
{
  int tnr;

  for (tnr = 0; tnr < MAX_TNR; tnr++)
    set_norm_xform(tnr, gksl->window[tnr], gksl->viewport[tnr]);
}

static
void set_xform(void)
{
  float a, b, c, d;

  a = (p->viewport[1] - p->viewport[0]) / (p->window[1] - p->window[0]);
  b = 810 / 0.288;
  p->a = a * b;
  p->b = b * (p->viewport[0] - p->window[0] * a);
  c = (p->viewport[3] - p->viewport[2]) / (p->window[3] - p->window[2]);
  d = 558 / 0.1984;
  p->c = c * d;
  p->d = d * (p->viewport[2] - p->window[2] * c);
}

static
void seg_xform(float *x, float *y)
{
  float xx;

  xx = *x * gksl->mat[0][0] + *y * gksl->mat[0][1] + gksl->mat[2][0];
  *y = *x * gksl->mat[1][0] + *y * gksl->mat[1][1] + gksl->mat[2][1];
  *x = xx;
}

static
void seg_xform_rel(float *x, float *y)
{
  float xx;

  xx = *x * gksl->mat[0][0] + *y * gksl->mat[0][1];
  *y = *x * gksl->mat[1][0] + *y * gksl->mat[1][1];
  *x = xx;
}

static
void set_color_rep(int color, float red, float green, float blue)
{
  if (color >= 0 && color < MAX_COLOR)
    {
      p->red[color] = red;
      p->green[color] = green;
      p->blue[color] = blue;
    }
}

static
void init_colors(void)
{
  int color;
  float red, green, blue;

  for (color = 0; color < MAX_COLOR; color++)
    {
      GQRGB (&color, &red, &green, &blue);
      p->red[color] = red;
      p->green[color] = green;
      p->blue[color] = blue;
    }
}

static
void create_patterns(void)
{
  register int i, j, k;
  int pattern, parray[33];

  for (i = 0; i < MAX_BITMAP; i++)
    {
      pattern = i;
      GKQPA(&pattern, parray);    
      for (j = 0, k = 1; j < 32; j += 2)
	{
	  sprintf(bitmap[i] + j, "%02x", parray[k]);
	  if (++k > *parray)
	    k = 1;
	}
      bitmap[i][64] = '\0';
    }
}

static
void open_ws(int fd, int wstype)
{
  p = (ws_state_list *)calloc(1, sizeof(struct ws_state_list_t));

  p->compress = wstype == 102;

  p->window[0] = p->window[2] = 0.0;
  p->window[1] = p->window[3] = 1.0;
  p->viewport[0] = p->viewport[2] = 0;
  p->viewport[1] = p->viewport[3] = 0.1984;
  p->width = p->height = 558;

  p->empty = 1;

  p->stroke = 0;
  p->lastx = p->lasty = -1;

  p->color = p->fillcolor = -1;
  p->ltype = -999; p->lwidth = -1.0;
  p->font = 1; p->size = 24; p->angle = 0;

  set_xform();

  pdf_open(fd);
}

static
void set_clip(float *clrt)
{
#if 0
  float x0, x1, y0, y1;

  NDC_to_DC(clrt[0], clrt[2], x0, y0);
  NDC_to_DC(clrt[1], clrt[3], x1, y1);

  pdf_moveto(p, x0, y0);
  pdf_lineto(p, x1, y0);
  pdf_lineto(p, x1, y1);
  pdf_lineto(p, x0, y1);
  pdf_clip(p);
#endif
}

static
void set_color(int color)
{
  if (color < MAX_COLOR)
    {
      if (p->color != color)
	{
	  pdf_setrgbcolor(p, p->red[color], p->green[color], p->blue[color]);
	  p->color = color;
	}
    }
}

static
void set_fillcolor(int color)
{
  if (color < MAX_COLOR)
    {
      if (p->fillcolor != color)
	{
	  pdf_setfillcolor(p, p->red[color], p->green[color], p->blue[color]);
	  p->fillcolor = color;
	}
    }
}

static
void begin_page(void)
{
  pdf_page(p, p->width, p->height);
  set_clip(p->window);
  p->empty = 0;

  p->color = p->fillcolor = -1;
  set_color(1);
  set_fillcolor(1);
}

static
void close_ws(void)
{
  pdf_close(p);

  free(p);
}

static
void clear_ws(void)
{
  p->empty = 1;
}

static
void stroke(void)
{
  if (p->stroke)
    {
      pdf_stroke(p);
      p->stroke = 0;
    }
}

static
void eofill(void)
{
  if (p->stroke)
    {
      pdf_eofill(p);
      p->stroke = 0;
    }
}

static
void move_routine(float *x, float *y)
{
  float xdev, ydev;

  stroke();

  NDC_to_DC(*x, *y, xdev, ydev);
  pdf_moveto(p, xdev, ydev);

  p->lastx = xdev;
  p->lasty = ydev;
}

static
void draw_routine(float *x, float *y)
{
  float xdev, ydev;

  NDC_to_DC(*x, *y, xdev, ydev);
  if (xdev != p->lastx || ydev != p->lasty)
    {
      pdf_lineto(p, xdev, ydev);
      p->lastx = xdev;
      p->lasty = ydev;
      p->stroke = 1;
    }
}

static
void line_routine(int *n, float *px, float *py, int *ltype, int *tnr)
{
  register int i;
  float x, y, xdev, ydev;

  for (i = 0; i < *n; i++)
    {
      WC_to_NDC(px[i], py[i], *tnr, x, y);
      seg_xform(&x, &y);
      NDC_to_DC(x, y, xdev, ydev);

      if (i == 0)
	pdf_moveto(p, xdev, ydev);
      else
	pdf_lineto(p, xdev, ydev);
    }

  p->stroke = 1;
  stroke();
}

static
void set_linetype(int ltype)
{
  if (p->ltype != ltype)
    {
      pdf_setdash(p, dash_pattern[ltype + 30]);
      p->ltype = ltype;
    }
}

static
void set_linewidth(float lwidth)
{
  if (p->lwidth != lwidth)
    {
      pdf_setlinewidth(p, lwidth);
      p->lwidth = lwidth;
    }
}

static
void polyline(int *n, float *px, float *py)
{
  int ln_type, ln_color;
  float ln_width;

  ln_type  = gksl->asf[0] ? gksl->ltype  : gksl->lindex;
  ln_width = gksl->asf[1] ? gksl->lwidth : 1;
  ln_color = gksl->asf[2] ? Color8Bit(gksl->plcoli) : 1;

  set_linetype(ln_type);
  set_linewidth(ln_width);
  set_color(ln_color);

  GSDT(p->window, p->viewport);
  GPOLIN(n, px, py, &ln_type, &gksl->cntnr, move_routine, draw_routine);
  stroke();
}

static
void draw_marker(float xn, float yn, int mtype, float mscale, int mcolor)
{
  int r, curve, i;
  float scale, x, y, xr, yr;
  int pc, op;

  static int marker[26][57] =
  {
    { 5, 9, -4, 7, 4, 7, 7, 4, 7, -4, 	/* omark */
      4, -7, -4, -7, -7, -4, -7, 4, 
      -4, 7,  3, 9, -4, 7, 4, 7, 7, 4, 
      7, -4, 4, -7, -4, -7, -7, -4, 
      -7, 4, -4, 7,  0 }, 
    { 5, 13, -2, 8, 2, 8, 2, 2, 8, 2, 	/* hollow plus */
      8, -2, 2, -2, 2, -8, -2, -8, 
      -2, -2, -8, -2, -8, 2, -2, 2, 
      -2, 8,  3, 13, -2, 8, 2, 8, 
      2, 2, 8, 2, 8, -2, 2, -2, 2, -8, 
      -2, -8, -2, -2, -8, -2, -8, 2, 
      -2, 2, -2, 8,  0 }, 
    { 4, 4, -8, 0, 4, 7, 4, -7, 	/* solid triangle right */
      -8, 0,  0 }, 
    { 4, 4, 8, 0, -4, -7, -4, 7,        /* solid triangle left */
      8, 0,  0 }, 
    { 5, 4, 0, 8, 7, -4, -7, -4, 0, 8,  /* triangle up down */
      5, 4, 0, -8, -7, 4, 7, 4, 0, -8, 
      3, 4, 0, 8, 7, -4, -7, -4, 0, 8, 
      3, 4, 0, -8, -7, 4, 7, 4, 0, -8, 
      0 }, 
    { 4, 11, 0, 9, 2, 2, 9, 3, 3, -1, 	/* solid star */
      6, -8, 0, -3, -6, -8, -3, -1, 
      -9, 3, -2, 2, 0, 9,  0 }, 
    { 5, 11, 0, 9, 2, 2, 9, 3, 3, -1,   /* hollow star */
      6, -8, 0, -3, -6, -8, -3, -1, 
      -9, 3, -2, 2, 0, 9, 
      3, 11, 0, 9, 2, 2, 9, 3, 3, -1, 
      6, -8, 0, -3, -6, -8, -3, -1, 
      -9, 3, -2, 2, 0, 9,  0 }, 
    { 4, 5, 0, 9, 9, 0, 0, -9, -9, 0,   /* solid diamond */
      0, 9,  0 }, 
    { 5, 5, 0, 9, 9, 0, 0, -9, -9, 0,   /* hollow diamond */
      0, 9,  3, 5, 0, 9, 9, 0, 0, -9, 
      -9, 0, 0, 9,  0 }, 
    { 4, 5, 9, 9, -9, -9, 9, -9, -9, 9, /* solid hourglass */
      9, 9,  0 }, 
    { 5, 5, 9, 9, -9, -9, 9, -9, -9, 9, /* hollow hourglass */
      9, 9,  3, 5, 9, 9, -9, -9, 9, -9, 
      -9, 9, 9, 9,  0 }, 
    { 4, 5, 9, 9, 9, -9, -9, 9, -9, -9, /* solid bowtie */
      9, 9,  0 }, 
    { 5, 5, 9, 9, 9, -9, -9, 9, -9, -9, /* hollow bowtie */
      9, 9,  3, 5, 9, 9, 9, -9, -9, 9, 
      -9, -9, 9, 9,  0 }, 
    { 4, 5, 9, 9, 9, -9, -9, -9, -9, 9, /* solid square */
      9, 9,  0 }, 
    { 5, 5, 9, 9, 9, -9, -9, -9, -9, 9, /* hollow square */
      9, 9,  3, 5, 9, 9, 9, -9, -9, -9, 
      -9, 9, 9, 9,  0 }, 
    { 4, 4, -9, 9, 9, 9, 0, -9, -9, 9,  /* solid triangle down */
      0 }, 
    { 5, 4, -9, 9, 9, 9, 0, -9, -9, 9,  /* hollow triangle down */
      3, 4, -9, 9, 9, 9, 0, -9, -9, 9, 
      0 }, 
    { 4, 4, 0, 9, 9, -9, -9, -9, 0, 9,  /* solid triangle up */
      0 }, 
    { 5, 4, 0, 9, 9, -9, -9, -9, 0, 9,  /* hollow triangle up */
      3, 4, 0, 9, 9, -9, -9, -9, 0, 9,
      0 }, 
    { 7,  0 }, 				/* solid circle */
    { 0 },     	        		/* not used */
    { 1,  0 },     	        	/* dot */
    { 2,  0, 0, 0, 9,  2, 0, 0, 9, 0,   /* plus */
      2, 0, 0, 0, -9,  2, 0, 0, -9, 0,  
      0 }, 
    { 2, 0, 0, 0, 9,  2, 0, 0, 9, 3,    /* asterisk */
      2, 0, 0, 6, -9,  2, 0, 0, -6, -9,  
      2, 0, 0, -9, 3,
      0 }, 
    { 8,  6,  0 }, 			/* circle */
    { 2, 0, 0, 9, 9,  2, 0, 0, 9, -9,   /* diagonal cross */
      2, 0, 0, -9, -9,  2, 0, 0, -9, 9, 
      0 }
  };

  static float cx[4][3] = {
    { 0.5523, 1, 1 },
    { 1, 0.5523, 0 },
    { -0.5523, -1, -1 },
    { -1, -0.5523, 0 }
  };

  static float cy[4][3] = {
    { -1, -0.5523, 0 },
    { 0.5523, 1, 1 },
    { 1, 0.5523, 0 },
    { -0.5523, -1, -1 }
  };

  r = (int)(3 * mscale);
  scale = mscale / 3.0;

  NDC_to_DC(xn, yn, x, y);

  pc = 0;    
  mtype = (r > 0.5) ? mtype + 20 : 21;

  do 
  {
    op = marker[mtype][pc];
    switch (op)
    {
      case 1: /* point */
	pdf_moveto(p, x, y);
	pdf_lineto(p, x, y);
	pdf_stroke(p);
	break;

      case 2: /* line */
	for (i = 0; i < 2; i++)
        {
          xr =  scale * marker[mtype][pc + 2 * i + 1];
          yr = -scale * marker[mtype][pc + 2 * i + 2];
          seg_xform_rel(&xr, &yr);
          if (i == 0)
	    pdf_moveto(p, x - xr, y - yr);
	  else
	    pdf_lineto(p, x - xr, y - yr);
        }
	pdf_stroke(p);
	pc += 4;
	break;

      case 3: /* polyline */
	for (i = 0; i < marker[mtype][pc + 1]; i++)
	{
	  xr =  scale * marker[mtype][pc + 2 + 2 * i];
	  yr = -scale * marker[mtype][pc + 3 + 2 * i];
          seg_xform_rel(&xr, &yr);
          if (i == 0)
	    pdf_moveto(p, x - xr, y - yr);
	  else
	    pdf_lineto(p, x - xr, y - yr);
	}
	pdf_stroke(p);
	pc += 1 + 2 * marker[mtype][pc + 1];
	break;

      case 4: /* filled polygon */
      case 5: /* hollow polygon */
	if (op == 5)
	  set_fillcolor(0);
	for (i = 0; i < marker[mtype][pc + 1]; i++)
	{
	  xr =  scale * marker[mtype][pc + 2 + 2 * i];
	  yr = -scale * marker[mtype][pc + 3 + 2 * i];
          seg_xform_rel(&xr, &yr);
          if (i == 0)
	    pdf_moveto(p, x - xr, y - yr);
	  else
	    pdf_lineto(p, x - xr, y - yr);
	}
	pdf_eofill(p);
        pc += 1 + 2 * marker[mtype][pc + 1];
	if (op == 5)
	  set_fillcolor(mcolor);
	break;

      case 6: /* arc */
	xr =  0;
	yr =  -r;
	seg_xform_rel(&xr, &yr);
	pdf_moveto(p, x - xr, y - yr);
	for (curve = 0; curve < 4; curve++)
	{
	  for (i = 0; i < 3; i++)
	  {
	    xr = r * cx[curve][i];
	    yr = r * cy[curve][i];
            seg_xform_rel(&xr, &yr);
	    pdf_point(p, x - xr, y - yr);
	  }
	  pdf_curveto(p);        
	}
	pdf_stroke(p);
	break;

      case 7: /* filled arc */
      case 8: /* hollow arc */
	if (op == 8)
	  set_fillcolor(0);
	xr =  0;
	yr =  -r;
	seg_xform_rel(&xr, &yr);
	pdf_moveto(p, x - xr, y - yr);
	for (curve = 0; curve < 4; curve++)
	{
	  for (i = 0; i < 3; i++)
	  {
	    xr = r * cx[curve][i];
	    yr = r * cy[curve][i];
            seg_xform_rel(&xr, &yr);
	    pdf_point(p, x - xr, y - yr);
	  }
	  pdf_curveto(p);        
	}
	pdf_eofill(p);
	if (op == 8)
	  set_fillcolor(mcolor);
	break;
    }
    pc++;
  }
  while (op != 0);
}


static
void marker_routine(
  int n, float *px, float *py, int mtype, float mscale, int mcolor)
{
  float x, y;
  float *clrt = gksl->viewport[gksl->cntnr];
  register int i, draw;

  for (i = 0; i < n; i++)
  {
    WC_to_NDC(px[i], py[i], gksl->cntnr, x, y);
    seg_xform(&x, &y);

    if (gksl->clip == GCLIP)
      draw = (x >= clrt[0] && x <= clrt[1] && y >= clrt[2] && y <= clrt[3]);
    else
      draw = 1;

    if (draw)
      draw_marker(x, y, mtype, mscale, mcolor);
  }
}

static
void polymarker(int *n, float *px, float *py)
{
  int mk_type, mk_color;
  float mk_size;

  mk_type  = gksl->asf[3] ? gksl->mtype  : gksl->mindex;
  mk_size  = gksl->asf[4] ? gksl->mszsc  : 1;
  mk_color = gksl->asf[5] ? Color8Bit(gksl->pmcoli) : 1;

  set_linetype(GLSOLI);
  set_linewidth(1.0);
  set_color(mk_color);

  marker_routine(*n, px, py, mk_type, mk_size, mk_color);
}

static
void set_font(int font)
{
  float ux, uy, scale;
  float width, height;
  PDF_page *page = p->page[p->current_page];

  font = abs(font);
  if (font >= 101 && font <= 143)
    font -= 100;
  else if (font >= 1 && font <= 32)
    font = map[font - 1];
  else
    font = 9;

  font--;
  if (!page->fonts[font])
    page->fonts[font] = pdf_alloc_id(p);
  p->font = font;
  
  WC_to_NDC_rel(gksl->chup[0], gksl->chup[1], gksl->cntnr, ux, uy);
  seg_xform_rel(&ux, &uy);

  p->angle = -atan2(ux, uy) * 180 / Pi;
  if (p->angle < 0) p->angle += 360;
  p->angle = (int)(p->angle + 0.5);

  scale = sqrt(gksl->chup[0]*gksl->chup[0] + gksl->chup[1]*gksl->chup[1]);
  ux = gksl->chup[0] / scale * gksl->chh;
  uy = gksl->chup[1] / scale * gksl->chh;
  WC_to_NDC_rel(ux, uy, gksl->cntnr, ux, uy);

  width = 0;
  height = sqrt(ux*ux + uy*uy);
  seg_xform_rel(&width, &height);

  height = sqrt(width*width + height*height);
  p->size = (int)(height * fabs(p->c) / 0.72 + 0.5);
}

static
void text_routine(float *x, float *y, int *nchars,
#ifdef VMS
    struct dsc$descriptor *text)
#else
#ifdef cray
    _fcd text)
#else
    char *text)
#endif /* cray */
#endif /* VMS */
{
  char *chars, s[BUFSIZ], *cp;
  float xrel, yrel, xorg, yorg;
  int tx_font, tx_prec;
  float phi, ax, ay;
  int width, ch, buff[256];
  register int i;

#ifdef VMS
  chars = text->dsc$a_pointer;
  *nchars = text->dsc$w_length;
#else
#ifdef cray
  chars = _fcdtocp(text);
  *nchars = _fcdlen(text);
#else
  chars = text;
#endif /* cray */
#endif /* VMS */

  NDC_to_DC(*x, *y, xorg, yorg);

  tx_font  = gksl->asf[6] ? gksl->txfont : predef_font[gksl->tindex - 1];
  tx_prec  = gksl->asf[6] ? gksl->txprec : predef_prec[gksl->tindex - 1];

  if (tx_prec == GSTRP)
    {
      width = 0;
      for (i = 0; i < *nchars; i++)
	{
	  ch = chars[i];
	  GKAFM(&tx_font, &ch, buff);
	  width += buff[1] - buff[0];
	}
      width = (int)(width * p->size * 0.72 / buff[2]);

      phi = p->angle * Pi / 180;
      xrel = width * xfac[gksl->txal[0]];
      yrel = p->size * yfac[gksl->txal[1]];

      CharXform(phi, xrel, yrel, ax, ay);

      xorg += ax;
      yorg += ay;
    }

  cp = s;
  for (i = 0; i < *nchars; i++)
    {
      ch = chars[i];
      if (ch == '(' || ch == ')' || ch == '\\')
	*cp++ = '\\';
      *cp++ = ch;
    }
  *cp = '\0';

  if (p->angle != 0)
    pdf_text_ex(p, xorg, yorg, s);
  else
    pdf_text(p, xorg, yorg, s);
}

static
void text(float *px, float *py, int nchars,
#ifdef VMS
    struct dsc$descriptor *chars)
#else
#ifdef cray
    _fcd chars)
#else
    char *chars)
#endif /* cray */
#endif /* VMS */
{
  int tx_font, tx_prec, tx_color;
  float x, y;
  int avail = 1;

  tx_font  = gksl->asf[6] ? gksl->txfont : predef_font[gksl->tindex - 1];
  tx_prec  = gksl->asf[6] ? gksl->txprec : predef_prec[gksl->tindex - 1];
  tx_color = gksl->asf[9] ? Color8Bit(gksl->txcoli) : 1;

  set_linetype(GLSOLI);
  set_linewidth(1.0);
  set_fillcolor(tx_color);

  if (tx_prec != GSTRKP)
    set_font(tx_font);

  if (tx_prec == GSTRP)
    {
      WC_to_NDC(*px, *py, gksl->cntnr, x, y);
      seg_xform(&x, &y);

      text_routine(&x, &y, &nchars, chars);
    }
  else
    {
#if defined(hpux) && !defined(NAGware)
      GTEXTS(px, py, &nchars, chars, &avail,
#if defined(__hp9000s700) || defined(__hp9000s300)
	line_routine_a, fill_routine_a, text_routine);
#else
	&line_routine_a, &fill_routine_a, &text_routine_a);
#endif /* __hp9000s700 || __hp9000s300 */
#else
      GTEXTS(px, py, &nchars, chars, &avail,
	line_routine, fill_routine, text_routine);
#endif
    }
}

static
void fill_routine(int *n, float *px, float *py, int *tnr)
{
  int ln_type = 0;

  GSDT(p->window, p->viewport);
  GPOLIN(n, px, py, &ln_type, tnr, move_routine, draw_routine);
  eofill();
}

static
void fillarea(int *n, float *px, float *py)
{
  int fl_color;

  fl_color = gksl->asf[12] ? Color8Bit(gksl->facoli) : 1;

  set_fillcolor(fl_color);

  fill_routine(n, px, py, &gksl->cntnr);
}

static
void cellarray(
  float xmin, float xmax, float ymin, float ymax,
  int dx, int dy, int dimx, int *colia)
{
  float x1, y1, x2, y2;
  int x, y, width, height;
  float rx1, rx2, ry1, ry2;
  register int i, j, ix, iy, ind, color;
  int swapx, swapy, count, chars_per_line;
  unsigned char data[3];

  WC_to_NDC(xmin, ymax, gksl->cntnr, x1, y1);
  seg_xform(&x1, &y1);
  NDC_to_DC(x1, y1, rx1, ry1);

  WC_to_NDC(xmax, ymin, gksl->cntnr, x2, y2);
  seg_xform(&x2, &y2);
  NDC_to_DC(x2, y2, rx2, ry2);

  width = (int)(fabs(rx2 - rx1) + 1);
  height = (int)(fabs(ry2 - ry1) + 1);
  x = (int)min(rx1, rx2);
  y = (int)min(ry1, ry2);

  swapx = rx1 > rx2;
  swapy = ry1 > ry2;

  pdf_save(p);
  pdf_printf(p->content, "%d 0 0 %d %d %d cm\n", width, height, x, y);
  pdf_printf(p->content, "BI\n");
  pdf_printf(p->content, "/W %d\n", dx);
  pdf_printf(p->content, "/H %d\n", dy);
  pdf_printf(p->content, "/BPC %d\n", 8);
  pdf_printf(p->content, "/CS /RGB\n");
  pdf_printf(p->content, "/F /AHx\n");
  pdf_printf(p->content, "ID ");

  chars_per_line = 0;
  for (j = 0; j < dy; j++)
    {
      iy = swapy ? dy - 1 - j : j;
      for (i = 0; i < dx; i++)
	{
	  ix = swapx ? dx - 1 - i : i;
	  ind = colia[iy * dimx + ix];
	  color = Color8Bit(ind);

	  data[0] = (Byte)(p->red[color] * 255);
	  data[1] = (Byte)(p->green[color] * 255);
	  data[2] = (Byte)(p->blue[color] * 255);
	  for (count = 0; count < 3; count++)
	    {
	      pdf_printf(p->content, "%02x", data[count]);
	      if ((chars_per_line += 2) >= 64)
		{
		  pdf_printf(p->content, "\n");
		  chars_per_line = 0;
		}
	    }
	}
    }

  pdf_printf(p->content, ">EI\n");
  pdf_restore(p);
}

void STDCALL GKDPDF(
  int *fctid, int *dx, int *dy, int *dimx, int *ia, int *lr1, float *r1,
  int *lr2, float *r2, int *lc, CHARARG(chars), ws_state_list **ptr)
{
  p = *ptr;

  switch (*fctid)
    {
      case  2:
/* open workstation */
	open_ws(ia[1], ia[2]);

	gksl = (gks_state_list *)(ia + 4);

	init_norm_xform();
	init_colors();
	create_patterns();

        if (sizeof(char *) > sizeof(int))
	  {
	    long *la = (long *)ia;
	    *la = (long)p;
	  }
        else
	  *ia = (long)p;

#if defined(hpux) && !defined(NAGware)
        line_routine_a = (int (*)())line_routine;
        fill_routine_a = (int (*)())fill_routine;
        text_routine_a = (int (*)())text_routine;
#endif
	break;

      case  3:
/* close workstation */
	close_ws();
	break;

      case  4:
/* activate workstation */
	p->state = GACTIV;
	break;

      case  5:
/* deactivate workstation */
	p->state = GINACT;
	break;

      case  6:
/* clear workstation */
	clear_ws();
	break;

      case  8:
/* update workstation */
	break;

      case 12:
/* polyline */
	if (p->state == GACTIV)
	  {
	    if (p->empty)
	      begin_page();
	    polyline(ia, r1, r2);
	  }
	break;

      case 13:
/* polymarker */
	if (p->state == GACTIV)
	  {
	    if (p->empty)
	      begin_page();
	    polymarker(ia, r1, r2);
	  }
	break;

      case 14:
/* text */
	if (p->state == GACTIV)
	  {
	    if (p->empty)
	      begin_page();
#ifdef VMS
            text(r1, r2, chars->dsc$w_length, chars);
#else
#ifdef cray
            text(r1, r2, _fcdlen(chars), chars);
#else
            text(r1, r2, strlen(chars), chars);
#endif /* cray */
#endif /* VMS */
	  }
	break;

      case 15:
/* fill area */
	if (p->state == GACTIV)
	  {
	    if (p->empty)
	      begin_page();
	    fillarea(ia, r1, r2);
	  }
	break;

    case 16:
/* cell array */
	if (p->state == GACTIV)
	  {
	    if (p->empty)
	      begin_page();
	    cellarray(r1[0], r1[1], r2[0], r2[1], *dx, *dy, *dimx, ia);
	  }
	break;

      case 48:
/* set color representation */
	set_color_rep(ia[1], r1[0], r1[1], r1[2]);
	break;

      case 49:
/* set window */
	set_norm_xform(*ia, gksl->window[*ia], gksl->viewport[*ia]);
	break;

      case 50:
/* set viewport */
	set_norm_xform(*ia, gksl->window[*ia], gksl->viewport[*ia]);
	break;

      case 54:
/* set workstation window */
	p->window[0] = r1[0];
	p->window[1] = r1[1];
	p->window[2] = r2[0];
	p->window[3] = r2[1];

	set_xform();
	init_norm_xform();
	if (!p->empty)
	  set_clip(p->window);
	break;

      case 55:
/* set workstation viewport */
	p->viewport[0] = r1[0];
	p->viewport[1] = r1[1];
	p->viewport[2] = r2[0];
	p->viewport[3] = r2[1];

	set_xform();
	init_norm_xform();
	break;

      default:
	;
    }
}
