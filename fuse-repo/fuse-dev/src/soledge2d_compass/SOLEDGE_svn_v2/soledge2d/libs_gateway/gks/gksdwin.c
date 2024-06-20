/*
 * Copyright @ 1997 - 1998   Josef Heinen
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
 *      This module contains a logical device driver for Microsoft Windows
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

#include "gksdefs.h"

#ifdef _WIN32

#include <math.h>

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#undef WIN32_LEAN_AND_MEAN

#if !defined(__GNUC__)
#define CHARARG(a) char *(a), unsigned short a##_len
#define STDCALL __stdcall
#else
#define CHARARG(a) char *a
#endif

#define MAX_TNR 9
#define MAX_COLOR 256
#define MAX_POINTS 2048
#define MAX_MESSAGES 80
#define MAX_BITMAP 120
#define HATCH_STYLE 108

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
  xd = (int)(p->a * (xn) + p->b + 0.5); \
  yd = (int)(p->c * (yn) + p->d + 0.5)

#define DC_to_NDC(xd, yd, xn, yn) \
    xn = ((xd) - p->b) / p->a; \
    yn = ((yd) - p->d) / p->c

#define CharXform(xrel, yrel, x, y) \
  x = cos_f[p->path] * (xrel) - sin_f[p->path] * (yrel); \
  y = sin_f[p->path] * (xrel) + cos_f[p->path] * (yrel);

#define nint(a) ((int)(a + 0.5))

#define Color8Bit(c) \
  c >= 588 ? 80 + (c - 588) / 56 * 12 + nint((c - 588) % 56 * 11.0 / 56.0) : \
  c >= 257 ? 8 + nint((c - 257) / 330.0 * (72 - 1)) : c

typedef struct ws_state_list_t
{
  int state;
  float window[4], viewport[4];
  int show;
  float width, height, mwidth, mheight;
  int swidth, sheight;
  float a, b, c, d;
  int path, capheight;
  DWORD thread;
  HWND win;
  HDC dc, memdc;
  HBITMAP bm, bm_old;
  int double_buffering;
  COLORREF palette[MAX_COLOR];
  DWORD pixel[MAX_COLOR];
  HBRUSH bg, brush;
  RECT rc[MAX_TNR];
  HRGN rgn;
  HFONT font;
  MSG msg[MAX_MESSAGES];
  int read_index, write_index;
} ws_state_list;

typedef struct win_t
{
  HANDLE instance, previnstance;
} win;

win gksw;

static
HANDLE threadevent;

static
ws_state_list *p;

static
gks_state_list *gksl;

static
POINT *points = NULL;

static
int max_points = MAX_POINTS;

static
float a[MAX_TNR], b[MAX_TNR], c[MAX_TNR], d[MAX_TNR];

static
float xfac[4] = { 0, 0, -0.5, -1 };

static
float yfac[6] = { 0, -1.2, -1, -0.5, 0, 0.2 };

static
float sin_f[] = { 0, 1, 0, -1 };

static
float cos_f[] = { 1, 0, -1, 0 };

static
int map[32] =
{
  22,  9,  5, 14, 18, 26, 13,  1,
  24, 11,  7, 16, 20, 28, 13,  3,
  23, 10,  6, 15, 19, 27, 13,  2,
  25, 12,  8, 17, 21, 29, 13,  4
};

static
HBITMAP bitmap[MAX_BITMAP];

static
int predef_font[] = { 1, 1, 1, -2, -3, -4 };

static
int predef_prec[] = { 0, 1, 2, 2, 2, 2 };

static
int predef_ints[] = { 0, 1, 3, 3, 3 };

static
int predef_styli[] = { 1, 1, 1, 2, 3 };

static
char *fonts[] = 
{
  "Times", "Helvetica", "Courier", "Symbol",
  "ITC Bookman", "NewCenturySchlbk", "AvantGarde", "Palatino"
};

static
float capheights[8] =
{
  0.72, 0.76, 0.6, 0.76, 0.76, 0.76, 0.76, 0.76
};

static
void fill_routine(int *n, float *px, float *py, int *tnr);

static
void set_clip_region(int tnr)
{
  RECT *rc;

  if (gksl->clip == GCLIP)
    rc = p->rc + tnr;
  else
    rc = p->rc;

  SetRectRgn(p->rgn, rc->left, rc->top, rc->right, rc->bottom);
}

static
void set_norm_xform(int tnr, float *wn, float *vp)
{
  RECT *rc = &p->rc[tnr];

  a[tnr] = (vp[1] - vp[0]) / (wn[1] - wn[0]);
  b[tnr] = vp[0] - wn[0] * a[tnr];
  c[tnr] = (vp[3] - vp[2]) / (wn[3] - wn[2]);
  d[tnr] = vp[2] - wn[2] * c[tnr];

  NDC_to_DC(vp[0], vp[3], rc->left, rc->top);
  NDC_to_DC(vp[1], vp[2], rc->right, rc->bottom);
  rc->right  += 1;
  rc->bottom += 1;
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
  p->a = (p->width - 1) / (p->window[1] - p->window[0]);
  p->b = -p->window[0] * p->a;
  p->c = (p->height - 1) / (p->window[2] - p->window[3]);
  p->d = p->height - 1 - p->window[2] * p->c;
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
void set_color(int color, float red, float green, float blue)
{
  BYTE r, g, b, tmp, c[4];

  if (color >= 0 && color < MAX_COLOR)
  {
    r = (BYTE)(red * 255);
    g = (BYTE)(green * 255);
    b = (BYTE)(blue * 255);
    p->palette[color] = RGB(r, g, b);
    memcpy(c, p->palette + color, 4);
    tmp = c[0]; c[0] = c[2]; c[2] = tmp;
    memcpy(p->pixel + color, c, 4);
  }
}

static
void init_colors(void)
{
  int color;
  float red, green, blue;
  BYTE r, g, b, tmp, c[4];

  for (color = 0; color < MAX_COLOR; color++)
  {
    GQRGB (&color, &red, &green, &blue);
    r = (BYTE)(red * 255);
    g = (BYTE)(green * 255);
    b = (BYTE)(blue * 255);
    p->palette[color] = RGB(r, g, b);
    memcpy(c, p->palette + color, 4);
    tmp = c[0]; c[0] = c[2]; c[2] = tmp;
    memcpy(p->pixel + color, c, 4);
  }
  p->bg = GetStockObject(WHITE_BRUSH);
}

static
void create_patterns(void)
{
  register int i, j;
  int pattern, parray[33], width, height;
  short int bits[32];

  for (i = 0; i < MAX_BITMAP; i++)
  {
    pattern = i;
    GKQPA(&pattern, parray);    
    width = height = (*parray == 32) ? 16 : (*parray == 4) ? 8 : *parray;
    for (j = 0; j < width; j++)
      bits[j] = parray[(j % *parray) + 1];
    bitmap[i] = CreateBitmap(width, height, 1, 1, bits);
  }
}

static
void create_bitmap(void)
{
  RECT rc;

  p->dc = GetDC(p->win);
  GetClientRect(p->win, &rc);
  p->bm = CreateCompatibleBitmap(p->dc, rc.right, rc.bottom);
  p->memdc = CreateCompatibleDC(p->dc);
  p->bm_old = SelectObject(p->memdc, p->bm);
  FillRect(p->memdc, &rc, p->bg);
  ReleaseDC(p->win, p->dc);
}

static
LONG FAR PASCAL windowproc(HWND win, UINT message, UINT wParam, LONG lParam)
{ 
  RECT rc;
  HDC dc;
   
  switch (message)
  {
    case WM_SIZE:
      if (wParam)
      {
        MoveWindow(p->win, 50, 50, LOWORD(lParam) + 7, HIWORD(lParam) + 26,
          TRUE);
        UpdateWindow(p->win);
        create_bitmap();
      }
      break;

    case WM_PAINT:
      if (p->bm)
      {
	dc = GetDC(p->win);
	GetClientRect(p->win, &rc);
	BitBlt(dc, 0, 0, rc.right, rc.bottom, p->memdc, 0, 0, SRCCOPY);
	ReleaseDC(p->win, dc);
      }
      p->msg[p->write_index].message = message;
      p->write_index = (p->write_index + 1) % MAX_MESSAGES;
      break;

    case WM_MOUSEMOVE:
    case WM_LBUTTONDOWN:
    case WM_MBUTTONDOWN:
    case WM_RBUTTONDOWN:
      p->msg[p->write_index].message = message;
      p->msg[p->write_index].wParam = wParam;
      p->msg[p->write_index].lParam = lParam;
      p->write_index = (p->write_index + 1) % MAX_MESSAGES;
      break;

    case WM_DESTROY:
      PostQuitMessage(0);
      return 0;
  }

  return DefWindowProc(win, message, wParam, lParam);
}

static
void get_window_extent(void)
{
  RECT rc;

  p->dc = GetDC(p->win);

  GetClientRect(p->win, &rc);
  p->width   = rc.right - rc.left;
  p->height  = rc.bottom - rc.top;
  p->mwidth  = 0.001 * GetDeviceCaps(p->dc, HORZSIZE);
  p->mheight = 0.001 * GetDeviceCaps(p->dc, VERTSIZE);
  p->swidth  = GetDeviceCaps(p->dc, HORZRES);
  p->sheight = GetDeviceCaps(p->dc, VERTRES);

  ReleaseDC(p->win, p->dc);
}

static
DWORD WINAPI threadfunc(LPVOID parm)
{
  WNDCLASS wndclass; 
  MSG msg;  
  RECT rc;
 
  if (!gksw.previnstance)
  {
    wndclass.style         = CS_HREDRAW | CS_VREDRAW | CS_CLASSDC;
    wndclass.lpfnWndProc   = windowproc;
    wndclass.cbClsExtra    = 0;
    wndclass.cbWndExtra    = 0;
    wndclass.hInstance     = gksw.instance;
    wndclass.hIcon         = LoadIcon(NULL, IDI_APPLICATION);
    wndclass.hCursor       = LoadCursor(NULL, IDC_ARROW); 
    wndclass.hbrBackground = GetStockObject(WHITE_BRUSH);
    wndclass.lpszMenuName  = NULL;
    wndclass.lpszClassName = "Gli";

    RegisterClass(&wndclass);    
  }
                                                      
  if (parm == NULL)
    p->win = CreateWindow("Gli", "gli", WS_OVERLAPPEDWINDOW,
      50, 50, p->width + 7, p->height + 26, NULL, NULL, gksw.instance, NULL);
  else
    p->win = (HWND)parm;
  
  GetClientRect(p->win, &rc);
  p->rgn = CreateRectRgn(rc.left, rc.top, rc.right, rc.bottom);

  UpdateWindow(p->win);

  create_bitmap();
  get_window_extent();

  SetEvent(threadevent);

  while (GetMessage(&msg, NULL, 0, 0))
  {
    TranslateMessage(&msg);
    DispatchMessage(&msg);
  }
  
  return msg.wParam;
}

static
void show_window(void)
{
  ShowWindow(p->win, SW_SHOWNORMAL);
  p->show = 1;
}

static
void open_ws(void)
{
  char *env;

  p = (ws_state_list *)malloc(sizeof(struct ws_state_list_t));

  max_points = MAX_POINTS;
  points = (POINT *)malloc(max_points * sizeof(POINT));

  p->window[0] = p->window[2] = 0.0;
  p->window[1] = p->window[3] = 1.0;
  p->viewport[0] = p->viewport[2] = 0;
  p->width = p->height = 500;

  p->bm_old = NULL;
  p->double_buffering = (char *)getenv("GLI_GKS_DOUBLE_BUF") != NULL ? 1 : 0;

  p->read_index = p->write_index = 0;

  env = (char *)getenv("GKSconid");
  if (env != NULL)
    sscanf(env, "%lx", &p->win);

  threadevent = CreateEvent(NULL, FALSE, FALSE, NULL);

  if (CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)threadfunc,	env != NULL ?
	p->win : NULL, 0, &p->thread) == NULL)
  {
    gks_fprintf(stderr, "can't create thread\n");
    exit(-1);
  }
  WaitForSingleObject(threadevent, INFINITE);
  CloseHandle(threadevent);

  p->show = 0;

  set_xform();

  p->viewport[1] = (float)p->width * p->mwidth / p->swidth;
  p->viewport[3] = (float)p->height * p->mheight / p->sheight;
}

static
void resize_window(void)
{
  float max_width, max_height;
  int width, height;
  RECT rc;

  max_width = p->mwidth * 0.925;
  max_height = p->mheight * 0.925;

  GKFVP(p->viewport, &max_width, &max_height);
  width = nint((p->viewport[1] - p->viewport[0]) / p->mwidth * p->swidth);
  height = nint((p->viewport[3] - p->viewport[2]) / p->mheight * p->sheight);

  if (p->width != width || p->height != height)
  {
    p->width = width;
    p->height = height;
    if (p->bm)
    {
      SelectObject(p->memdc, p->bm_old);
      DeleteObject(p->bm);
      p->bm = p->bm_old = NULL;
      DeleteDC(p->memdc);
      SendMessage(p->win, WM_SIZE, 1, MAKELONG(p->width, p->height));
    }
  }
}

static
void close_ws(void)
{
  int i, tnr;

  for (i = 0; i < MAX_BITMAP; i++)
    DeleteObject(bitmap[i]);

  DeleteObject(p->rgn);
  if (p->thread != 0)
  {
    SendMessage(p->win, WM_DESTROY, 0, 0L);
    CloseHandle(&p->thread);
  }

  free(points);
  free(p);
}

static
void clear_ws(void)
{
  RECT rc;

  set_clip_region(0);

  p->dc = GetDC(p->win);
  GetClientRect(p->win, &rc);

  if (!p->double_buffering)
  {
    SelectClipRgn(p->dc, p->rgn);
    FillRect(p->dc, &rc, p->bg);
  }
  if (p->bm)
  {
    SelectClipRgn(p->memdc, p->rgn);
    FillRect(p->memdc, &rc, p->bg);
  }
  ReleaseDC(p->win, p->dc);

  set_clip_region(gksl->cntnr);
}

static
void update_ws(void)
{
  RECT rc;

  if (p->bm)
  {
    p->dc = GetDC(p->win);
    GetClientRect(p->win, &rc);
    BitBlt(p->dc, 0, 0, rc.right, rc.bottom, p->memdc, 0, 0, SRCCOPY);
    ReleaseDC(p->win, p->dc);
  }
}

static
void move_to(float *x, float *y)
{
  int ix, iy;

  NDC_to_DC(*x, *y, ix, iy);
  if (!p->double_buffering)
    MoveToEx(p->dc, ix, iy, NULL);
  if (p->bm)
    MoveToEx(p->memdc, ix, iy, NULL);
}

static
void line_to(float *x, float *y)
{
  int ix, iy;

  NDC_to_DC(*x, *y, ix, iy);
  if (!p->double_buffering)
    LineTo(p->dc, ix, iy);
  if (p->bm)
    LineTo(p->memdc, ix, iy);
}

static
void move_routine(float *x, float *y)
{
  GMOVE(x, y, move_to);
}

static
void draw_routine(float *x, float *y)
{
  GDASH(x, y, move_to, line_to);
}

static
void line_routine(int *n, float *px, float *py, int *linetype, int *tnr)
{
  float x, y;
  int i, ix, iy, ix0, iy0;

  WC_to_NDC(px[0], py[0], *tnr, x, y);
  seg_xform(&x, &y);
  NDC_to_DC(x, y, ix, iy);
  ix0 = ix;
  iy0 = iy;

  if (!p->double_buffering)
    MoveToEx(p->dc, ix, iy, NULL);
  if (p->bm)
    MoveToEx(p->memdc, ix, iy, NULL);
  for (i = 1;  i < *n;  i++)
  {
    WC_to_NDC(px[i], py[i], *tnr, x, y);
    seg_xform(&x, &y);
    NDC_to_DC(x, y, ix, iy);

    if (!p->double_buffering)
      LineTo(p->dc, ix, iy);
    if (p->bm)
      LineTo(p->memdc, ix, iy);
  }
  if (*linetype == 0)
  {
    if (!p->double_buffering)
      LineTo(p->dc, ix0, iy0);
    if (p->bm)
      LineTo(p->memdc, ix0, iy0);
  }
}

static
void polyline(int *n, float *px, float *py)
{
  int ln_type, ln_color;
  float ln_width;
  HPEN pen, dc_pen, memdc_pen;

  ln_type  = gksl->asf[0] ? gksl->ltype  : gksl->lindex;
  ln_width = gksl->asf[1] ? gksl->lwidth : 1;
  ln_color = gksl->asf[1] ? Color8Bit(gksl->plcoli) : 1;

  p->dc = GetDC(p->win);
  if (!p->double_buffering)
    SelectClipRgn(p->dc, p->rgn);
  if (p->bm)
    SelectClipRgn(p->memdc, p->rgn);
  if (ln_width > 1)
  {
    DWORD style = PS_GEOMETRIC | PS_COSMETIC | PS_ENDCAP_FLAT | PS_JOIN_ROUND;
    LOGBRUSH lb;

    lb.lbStyle = BS_SOLID;
    lb.lbColor = p->palette[ln_color];
    lb.lbHatch = 0;

    pen = ExtCreatePen(style, ln_width, &lb, 0, NULL);
  }
  else
    pen = CreatePen(PS_SOLID, ln_width, p->palette[ln_color]);
  if (!p->double_buffering)
  {
    dc_pen = SelectObject(p->dc, pen);
    SetBkMode(p->dc, TRANSPARENT);
  }
  if (p->bm)
  {
    memdc_pen = SelectObject(p->memdc, pen);
    SetBkMode(p->memdc, TRANSPARENT);
  }
  if (ln_type < 0 || ln_type > 1)
  {
    GSDT(p->window, p->viewport);
    GPOLIN(n, px, py, &ln_type, &gksl->cntnr, move_routine, draw_routine);
  }
  else
    line_routine(n, px, py, &ln_type, &gksl->cntnr);

  if (!p->double_buffering)
    SelectObject(p->dc, dc_pen);
  if (p->bm)
    SelectObject(p->memdc, memdc_pen);
  DeleteObject(pen);
  ReleaseDC(p->win, p->dc);
}

static
void draw_marker(float xn, float yn, int mtype, int mscale, int mcolor)
{
  int r, x, y, i;
  int pc, op;
  float scale, xr, yr;
  POINT points[13];

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
	if (!p->double_buffering)
	  SetPixel(p->dc, x, y, p->palette[mcolor]);
        if (p->bm)
          SetPixel(p->memdc, x, y, p->palette[mcolor]);
	break;

      case 2: /* line */
	for (i = 0; i < 2; i++)
        {
          xr =  scale * marker[mtype][pc + 2 * i + 1];
          yr = -scale * marker[mtype][pc + 2 * i + 2];
          seg_xform_rel(&xr, &yr);
          points[i].x = nint(x - xr);
          points[i].y = nint(y + yr);
        }
	if (!p->double_buffering)
	  Polyline(p->dc, points, 2);
	if (p->bm)
	  Polyline(p->memdc, points, 2);
	pc += 4;
	break;

      case 3: /* polyline */
	for (i = 0; i < marker[mtype][pc + 1]; i++)
	{
	  xr =  scale * marker[mtype][pc + 2 + 2 * i];
	  yr = -scale * marker[mtype][pc + 3 + 2 * i];
          seg_xform_rel(&xr, &yr);
          points[i].x = nint(x - xr);
          points[i].y = nint(y + yr);
	}
	if (!p->double_buffering)
	  Polyline(p->dc, points, marker[mtype][pc + 1]);
	if (p->bm)
	  Polyline(p->memdc, points, marker[mtype][pc + 1]);
	pc += 1 + 2 * marker[mtype][pc + 1];
	break;

      case 4: /* filled polygon */
      case 5: /* hollow polygon */
	if (op == 5)
	{
	  if (!p->double_buffering)
	    SelectObject(p->dc, p->bg);
	  if (p->bm)
	    SelectObject(p->memdc, p->bg);
	}
	for (i = 0; i < marker[mtype][pc + 1]; i++)
	{
	  xr =  scale * marker[mtype][pc + 2 + 2 * i];
	  yr = -scale * marker[mtype][pc + 3 + 2 * i];
          seg_xform_rel(&xr, &yr);
          points[i].x = nint(x - xr);
          points[i].y = nint(y + yr);
	}
	if (!p->double_buffering)
	  Polygon(p->dc, points, marker[mtype][pc + 1]);
	if (p->bm)
	  Polygon(p->memdc, points, marker[mtype][pc + 1]);
        pc += 1 + 2 * marker[mtype][pc + 1];
	if (op == 5)
	{
	  if (!p->double_buffering)
	    SelectObject(p->dc, p->brush);
	  if (p->bm)
	    SelectObject(p->memdc, p->brush);
	}
	break;

      case 6: /* arc */
	if (!p->double_buffering)
	  Arc(p->dc, x - r, y - r, x + r, y + r, x - r, y - r, x - r, y - r);
	if (p->bm)
	  Arc(p->memdc, x - r, y - r, x + r, y + r, x - r, y - r, x - r, y - r);
	break;

      case 7: /* filled arc */
      case 8: /* hollow arc */
	if (op == 8)
	{
	  if (!p->double_buffering)
	    SelectObject(p->dc, p->bg);
	  if (p->bm)
	    SelectObject(p->memdc, p->bg);
	}
	if (!p->double_buffering)
	  Chord(p->dc, x - r, y - r, x + r, y + r, x - r, y - r, x - r, y - r);
	if (p->bm)
	  Chord(p->memdc, x - r, y - r, x + r, y + r, x - r, y - r, x - r,
	    y - r);
	if (op == 8)
	{
	  if (!p->double_buffering)
	    SelectObject(p->dc, p->brush);
	  if (p->bm)
	    SelectObject(p->memdc, p->brush);
	}
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
  HPEN pen, dc_pen, memdc_pen;
  HBRUSH dc_brush, memdc_brush;

  mk_type  = gksl->asf[3] ? gksl->mtype  : gksl->mindex;
  mk_size  = gksl->asf[4] ? gksl->mszsc  : 1;
  mk_color = gksl->asf[5] ? Color8Bit(gksl->pmcoli) : 1;

  p->dc = GetDC(p->win);
  pen = CreatePen(PS_SOLID, 1, p->palette[mk_color]);
  p->brush = CreateSolidBrush(p->palette[mk_color]);
  if (!p->double_buffering)
  {
    dc_pen = SelectObject(p->dc, pen);
    SetBkMode(p->dc, OPAQUE);
    dc_brush = SelectObject(p->dc, p->brush);
  }
  if (p->bm)
  {
    memdc_pen = SelectObject(p->memdc, pen);
    SetBkMode(p->memdc, OPAQUE);
    memdc_brush = SelectObject(p->memdc, p->brush);
  }
  marker_routine(*n, px, py, mk_type, mk_size, mk_color);

  if (!p->double_buffering)
    SelectObject(p->dc, dc_brush);
  if (p->bm)
    SelectObject(p->memdc, memdc_brush);
  DeleteObject(p->brush);
  if (!p->double_buffering)
    SelectObject(p->dc, dc_pen);
  if (p->bm)
    SelectObject(p->memdc, memdc_pen);
  DeleteObject(pen);
  ReleaseDC(p->win, p->dc);
}

static
void draw_text(int x, int y, int width, char *chars, int nchars)
{
  TEXTMETRIC tm;
  HDC from, to;
  HBITMAP src, dest;
  RECT rc;
  COLORREF pixel, bg;
  register int height, descent, w, h;
  register int i, j, ii, jj;

  if (p->path != 0)
  {
    GetTextMetrics(p->dc, &tm);
    height = tm.tmHeight;
    descent = tm.tmDescent;

    switch (p->path)
    {
      case 1: x -= height - descent; y -= width; w = height; h = width; break;
      case 2: x -= width; y -= descent; w = width; h = height; break;
      case 3: x -= descent; w = height; h = width; break;
    }

    from = CreateCompatibleDC(p->dc);
    src = CreateCompatibleBitmap(p->dc, width, height);
    SelectObject(from, src);

    SelectObject(from, p->font);
    SetTextColor(from, GetTextColor(p->dc));
    SetTextAlign(from, TA_LEFT | TA_BASELINE);
    SetBkMode(from, TRANSPARENT);

    rc.left = 0;
    rc.top = 0;
    rc.right = width;
    rc.bottom = height;
    FillRect(from, &rc, p->bg);
    bg = GetPixel(from, 0, 0);
    ExtTextOut(from, 0, height - descent, 0, NULL, chars, nchars, NULL);

    to = CreateCompatibleDC(p->dc);
    dest = CreateCompatibleBitmap(p->dc, w, h);
    SelectObject(to, dest);

    BitBlt(to, 0, 0, w, h, p->memdc != NULL ? p->memdc : p->dc, x, y, SRCCOPY);

    for (i = 0; i < width; i++)
    {
      for (j = 0; j < height; j++)
      {
        switch (p->path)
        {
          case 1: ii = j; jj = h - i - 1; break;
          case 2: ii = w - i - 1; jj = h - j - 1; break;
          case 3: ii = w - j - 1; jj = i; break;
        }
        pixel = GetPixel(from, i, j);
        if (pixel != bg)
          SetPixel(to, ii, jj, pixel);
      }
    }

    if (!p->double_buffering)
      BitBlt(p->dc, x, y, w, h, to, 0, 0, SRCCOPY);
    if (p->bm)
      BitBlt(p->memdc, x, y, w, h, to, 0, 0, SRCCOPY);

    DeleteDC(to);
    DeleteObject(dest);
    DeleteDC(from);
    DeleteObject(src);
  }
  else
  {
    if (!p->double_buffering)
      ExtTextOut(p->dc, x, y, 0, NULL, chars, nchars, NULL);
    if (p->bm)
      ExtTextOut(p->memdc, x, y, 0, NULL, chars, nchars, NULL);
  }
}

static
void text_routine(float x, float y, int nchars, char *chars)
{
  int xstart, ystart;
  SIZE size;
  float xrel, yrel, ax, ay;

  NDC_to_DC(x, y, xstart, ystart);

  GetTextExtentPoint32(p->dc, chars, nchars, &size);
  xrel = size.cx * xfac[gksl->txal[0]];
  yrel = p->capheight * yfac[gksl->txal[1]];
  CharXform(xrel, yrel, ax, ay);
  xstart += (int)ax;
  ystart -= (int)ay;

  draw_text(xstart, ystart, size.cx, chars, nchars);
}

static
void create_font(int font)
{
  float rad, scale, ux, uy;
  int family, size, angle;
  float width, height, capheight;
  LOGFONT lf;

  font = abs(font);
  if (font >= 101 && font <= 129)
    font -= 100;
  else if (font >= 1 && font <= 32)
    font = map[font - 1];
  else
    font = 9;

  if (font > 13) font += 3;
  family = (font - 1) / 4;

  WC_to_NDC_rel(gksl->chup[0], gksl->chup[1], gksl->cntnr, ux, uy);
  seg_xform_rel(&ux, &uy);

  rad = -atan2(ux, uy);
  angle = (int)(rad * 180 / Pi + 0.5);
  if (angle < 0) angle += 360;
  p->path = ((angle + 45) / 90) % 4;

  scale = sqrt(gksl->chup[0] * gksl->chup[0] + gksl->chup[1] * gksl->chup[1]);
  ux = gksl->chup[0] / scale * gksl->chh;
  uy = gksl->chup[1] / scale * gksl->chh;
  WC_to_NDC_rel(ux, uy, gksl->cntnr, ux, uy);

  width = 0;
  height = sqrt(ux * ux + uy * uy);
  seg_xform_rel(&width, &height);

  height = sqrt(width * width + height * height);
  capheight = nint(height * (fabs(p->c) + 1));
  p->capheight = nint(capheight);

  memset((void *)&lf, 0, sizeof(LOGFONT));
  lf.lfHeight = - nint(capheight / capheights[family]);
  lf.lfWeight = (font % 4 == 1 || font % 4 == 2) ? FW_NORMAL : FW_BOLD;
  lf.lfItalic = (font % 4 == 2 || font % 4 == 0);
  strcpy(lf.lfFaceName, fonts[family]);
  p->font = CreateFontIndirect(&lf);
}

static
void text(float *px, float *py, int nchars, char *chars)
{
  int tx_font, tx_prec, tx_color;
  float x, y;
  unsigned short chars_len = nchars;
  HPEN pen, dc_pen, memdc_pen;
  HFONT dc_font, memdc_font;

  tx_font  = gksl->asf[6] ? gksl->txfont : predef_font[gksl->tindex - 1];
  tx_prec  = gksl->asf[6] ? gksl->txprec : predef_prec[gksl->tindex - 1];
  tx_color = gksl->asf[9] ? Color8Bit(gksl->txcoli) : 1;

  p->dc = GetDC(p->win);
  if (!p->double_buffering)
    SelectClipRgn(p->dc, p->rgn);
  if (p->bm)
    SelectClipRgn(p->memdc, p->rgn);
  if (tx_prec == GSTRP)
  {
    create_font(tx_font);

    WC_to_NDC(*px, *py, gksl->cntnr, x, y);
    seg_xform(&x, &y);

    if (!p->double_buffering)
    {
      dc_font = SelectObject(p->dc, p->font);
      SetTextColor(p->dc, p->palette[tx_color]);
      SetTextAlign(p->dc, TA_LEFT | TA_BASELINE);
      SetBkMode(p->dc, TRANSPARENT);
    }
    if (p->bm)
    {
      memdc_font = SelectObject(p->memdc, p->font);
      SetTextColor(p->memdc, p->palette[tx_color]);
      SetTextAlign(p->memdc, TA_LEFT | TA_BASELINE);
      SetBkMode(p->memdc, TRANSPARENT);
    }
    text_routine(x, y, nchars, chars);

    if (!p->double_buffering)
      SelectObject(p->dc, dc_font);
    if (p->bm)
      SelectObject(p->memdc, memdc_font);
    DeleteObject(p->font);
  }
  else
  {
    pen = CreatePen(PS_SOLID, 1, p->palette[tx_color]);
    if (!p->double_buffering)
      dc_pen = SelectObject(p->dc, pen);
    if (p->bm)
      memdc_pen = SelectObject(p->memdc, pen);

#ifndef __GNUC__
    GSIMTX(px, py, &nchars, chars, chars_len, line_routine, fill_routine);
#else
    GSIMTX(px, py, &nchars, chars, line_routine, fill_routine, chars_len);
#endif

    if (!p->double_buffering)
      SelectObject(p->memdc, dc_pen);
    if (p->bm)
      SelectObject(p->memdc, memdc_pen);
    DeleteObject(pen);
  }

  ReleaseDC(p->win, p->dc);
}

static
void fill_routine(int *n, float *px, float *py, int *tnr)
{
  register int i;
  float x, y;

  if (*n > max_points)
  {
    points = (POINT *)realloc(points, *n * sizeof(POINT));
    max_points = *n;
  }

  for (i = 0; i < *n; i++)
  {
    WC_to_NDC(px[i], py[i], *tnr, x, y);
    seg_xform(&x, &y);
    NDC_to_DC(x, y, points[i].x, points[i].y);
  }

  if (!p->double_buffering)
    Polygon(p->dc, points, *n);
  if (p->bm)
    Polygon(p->memdc, points, *n);
}

static
void fillarea(int *n, float *px, float *py)
{
  int fl_inter, fl_style, fl_color;
  HPEN pen, dc_pen, memdc_pen;
  HBRUSH dc_brush, memdc_brush;

  fl_inter = gksl->asf[10] ? gksl->ints   : predef_ints[gksl->findex - 1];
  fl_style = gksl->asf[11] ? gksl->styli  : predef_styli[gksl->findex - 1];
  fl_color = gksl->asf[12] ? Color8Bit(gksl->facoli) : 1;

  p->dc = GetDC(p->win);
  if (!p->double_buffering)
    SelectClipRgn(p->dc, p->rgn);
  if (p->bm)
    SelectClipRgn(p->memdc, p->rgn);
  pen = CreatePen(PS_NULL, 1, p->palette[fl_color]);
  if (fl_inter == GHOLLO)
    p->brush = p->bg;
  else if (fl_inter == GSOLID)
    p->brush = CreateSolidBrush(p->palette[fl_color]);
  else if (fl_inter == GPATTR || fl_inter == GHATCH)
  {
    if (fl_inter == GHATCH)
      fl_style += HATCH_STYLE;
    if (fl_style >= MAX_BITMAP)
      fl_style = 1;
    p->brush = CreatePatternBrush(bitmap[fl_style]);
  }
  if (!p->double_buffering)
  {
    dc_pen = SelectObject(p->dc, pen);
    SetBkMode(p->dc, OPAQUE);
    dc_brush = SelectObject(p->dc, p->brush);
  }
  if (p->bm)
  {
    memdc_pen = SelectObject(p->memdc, pen);
    SetBkMode(p->memdc, OPAQUE);
    memdc_brush = SelectObject(p->memdc, p->brush);
  }
  fill_routine(n, px, py, &gksl->cntnr);

  if (!p->double_buffering)
    SelectObject(p->dc, dc_brush);
  if (p->bm)
    SelectObject(p->memdc, memdc_brush);
  DeleteObject(p->brush);
  if (!p->double_buffering)
    SelectObject(p->dc, dc_pen);
  if (p->bm)
    SelectObject(p->memdc, memdc_pen);
  DeleteObject(pen);
  ReleaseDC(p->win, p->dc);
}

static
void cellarray(
  float xmin, float xmax, float ymin, float ymax,
  int dx, int dy, int dimx, int *colia)
{
  float x1, y1, x2, y2;
  int ix1, ix2, iy1, iy2;
  int x, y, width, height;
  register int i, j, ix, iy, ind;
  int swapx, swapy;
  HDC hdc;
  HBITMAP hbm_old, hbm;
  DWORD *pix = NULL;
  LPBITMAPINFO bmi;

  WC_to_NDC(xmin, ymax, gksl->cntnr, x1, y1);
  seg_xform(&x1, &y1);
  NDC_to_DC(x1, y1, ix1, iy1);

  WC_to_NDC(xmax, ymin, gksl->cntnr, x2, y2);
  seg_xform(&x2, &y2);
  NDC_to_DC(x2, y2, ix2, iy2);

  width = abs(ix2 - ix1) + 1;
  height = abs(iy2 - iy1) + 1;
  x = min(ix1, ix2);
  y = min(iy1, iy2);

  swapx = ix1 > ix2;
  swapy = iy1 > iy2;

  p->dc = GetDC(p->win);
  if (!p->double_buffering)
    SelectClipRgn(p->dc, p->rgn);
  if (p->bm)
    SelectClipRgn(p->memdc, p->rgn);

  bmi = (BITMAPINFO *)malloc(sizeof(BITMAPINFO) + (256 * sizeof(RGBQUAD)));
  bmi->bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
  bmi->bmiHeader.biWidth = width;
  bmi->bmiHeader.biHeight = height;
  bmi->bmiHeader.biPlanes = 1;
  bmi->bmiHeader.biBitCount = 32;
  bmi->bmiHeader.biCompression = BI_RGB;
  bmi->bmiHeader.biSizeImage = 0;
  bmi->bmiHeader.biXPelsPerMeter = 0;
  bmi->bmiHeader.biYPelsPerMeter = 0;
  bmi->bmiHeader.biClrUsed = 0;
  bmi->bmiHeader.biClrImportant = 0;

  hdc = CreateCompatibleDC(p->dc);
  hbm = CreateDIBSection(hdc, bmi, DIB_RGB_COLORS, (void *)&pix, NULL, 0);
  hbm_old = SelectObject(hdc, hbm);

  for (j = 0; j < height; j++)
  {
    iy = dy * j / height;
    if (swapy)
      iy = dy - 1 - iy;
    for (i = 0; i < width; i++)
    {
      ix = dx * i / width;
      if (swapx)
        ix = dx - 1 - ix;
      ind = colia[iy * dimx + ix];
      pix[j * width + i] = p->pixel[Color8Bit(ind)];
    }
  }

  if (!p->double_buffering)
    BitBlt(p->dc, x, y, width, height, hdc, 0, 0, SRCCOPY);  
  if (p->bm)
    BitBlt(p->memdc, x, y, width, height, hdc, 0, 0, SRCCOPY);  

  SelectObject(hdc, hbm_old);
  DeleteObject(hbm);
  DeleteDC(hdc);

  free(bmi);

  ReleaseDC(p->win, p->dc);
}

static
void crosshair(int mousex, int mousey)
{
  int drawmode;

  p->dc = GetDC(p->win);
  drawmode = SetROP2(p->dc, R2_NOT);

  MoveToEx(p->dc, mousex, 0, NULL); LineTo(p->dc, mousex, p->height);
  MoveToEx(p->dc, 0, mousey, NULL); LineTo(p->dc, p->width, mousey);

  SetROP2(p->dc, drawmode);
  ReleaseDC(p->win, p->dc);
}

static
void get_pointer(float *x, float *y, int *state)
{
  MSG *msg;
  int mousex = -1, mousey = -1;

  p->read_index = p->write_index;
  for (;;)
  {
    while (p->read_index == p->write_index)
      Sleep(10);

    msg = p->msg + p->read_index;
    if (msg->message != WM_PAINT)
    {  
      if (mousex >= 0 && mousey >= 0)
        crosshair(mousex, mousey);

      mousex = LOWORD(msg->lParam);
      mousey = HIWORD(msg->lParam);
    }
    if (msg->message == WM_MOUSEMOVE || msg->message == WM_PAINT)
    {
      if (mousex >= 0 && mousey >= 0)
        crosshair(mousex, mousey);
    }
    else
      break;

    p->read_index = (p->read_index + 1) % MAX_MESSAGES;
  }

  DC_to_NDC(mousex, mousey, *x, *y);
  if (msg->message == WM_LBUTTONDOWN)
    *state = GOK;
  else
    *state = GNONE;
}

static
void get_string(int *n, char *chars, int *state)
{
  *n = 0;
  *chars = '\0';
  *state = GNONE;
}

#if !defined(__GNUC__)
void STDCALL GKDGA(
#else
void GKDGA(
#endif
    int *fctid, int *dx, int *dy, int *dimx, int *ia, int *lr1, float *r1,
    int *lr2, float *r2, int *lc, CHARARG(chars))
{
  switch (*fctid)
  {

    case  2:
// open workstation
      open_ws();

      gksl = (gks_state_list *)(ia + 4);

      init_norm_xform();
      init_colors();
      create_patterns();
      break;

    case  3:
// close workstation
      close_ws();
      break;

    case  4:
// activate workstation
      p->state = GACTIV;
      break;

    case  5:
// deactivate workstation
      p->state = GINACT;
      break;

    case  6:
// clear workstation
      clear_ws();
      break;

    case  8:
// update workstation
      if (p->double_buffering && ia[1] == GPERFO)
	update_ws();
      break;

    case 12:
// polyline
      if (!p->show)
        show_window();
      if (p->state == GACTIV)
        polyline(ia, r1, r2);
      break;

    case 13:
// polymarker
      if (!p->show)
        show_window();
      if (p->state == GACTIV)
        polymarker(ia, r1, r2);
      break;

    case 14:
// text
      if (!p->show)
        show_window();
      if (p->state == GACTIV)
        text(r1, r2, strlen(chars), chars);
      break;

    case 15:
// fill area
      if (!p->show)
        show_window();
      if (p->state == GACTIV)
        fillarea(ia, r1, r2);
      break;

    case 16:
// cell array
      if (!p->show)
        show_window();
      if (p->state == GACTIV)
        cellarray(r1[0], r1[1], r2[0], r2[1], *dx, *dy, *dimx, ia);
      break;

    case 48:
// set color representation
      set_color(ia[1], r1[0], r1[1], r1[2]);
      break;

    case 49:
// set window
      set_norm_xform(*ia, gksl->window[*ia], gksl->viewport[*ia]);
      break;

    case 50:
// set viewport
      set_norm_xform(*ia, gksl->window[*ia], gksl->viewport[*ia]);
      if (*ia == gksl->cntnr)
        set_clip_region(*ia);
      break;

    case 52:
// select normalization transformation
    case 53:
// set clipping inidicator
      set_clip_region(gksl->cntnr);
      break;

    case 54:
// set workstation window
      p->window[0] = r1[0];
      p->window[1] = r1[1];
      p->window[2] = r2[0];
      p->window[3] = r2[1];

      set_xform();
      init_norm_xform();
      break;

    case 55:
// set workstation viewport
      p->viewport[0] = r1[0];
      p->viewport[1] = r1[1];
      p->viewport[2] = r2[0];
      p->viewport[3] = r2[1];

      resize_window();
      set_xform();
      init_norm_xform();
      break;

    case 81:
// request locator
      if (!p->show)
        show_window();
      get_pointer(r1, r2, &ia[0]);
      break;

    case 83:
// request string
      if (!p->show)
        show_window();
      get_string(&ia[1], chars, &ia[0]);
      break;

    default:
      ;
  }
}

#else

void GKDGA(
  int *fctid, int *dx, int *dy, int *dimx, int *ia,
  int *lr1, float *r1, int *lr2, float *r2, int *lc, char *chars)
{
  if (*fctid == 2)
  {
    gks_fprintf(stderr, "GKS: GDI32 not supported for this type of machine\n");
    ia[0] = 0;
  }
}

#endif /* _WIN32 */
