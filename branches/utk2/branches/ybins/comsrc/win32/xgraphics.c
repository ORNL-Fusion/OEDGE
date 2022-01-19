/*    xGraphics.c
 *
 *    This file contains all the general purpose graphics output routines.
 *    g1tran calls these routines. For most graphical operations the data
 *    is stored in the Replay Buffer via the xReplay.c routines and then
 *    the appropriate Xlib routine is called to do the graphics "live".
 *    In the case of an expose (etc.) then the eventhandler in xEvent is
 *    called and the Xlib routines are called again there.
 *
 *    In all cases the routine calls Xlib routines only if the 
 *    server connection is active.
 *
 *    See VERSION and graphicversion (below) for history.
 *
 */

#include <stdio.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>

#include "opcodes.h"
#include "defines.h"

extern Display *GHdisplay;
extern Window GHwin;
extern GC GHgc;
extern XFontStruct *font_info;

/* current window addressability (size) DECLARED: xInit.c */
extern int Dresx, Dresy;
extern int origDresx, origDresy;
extern int RealDresy;

/* user default items DECLARED xInit.c */
extern Bool keep_ratio;

/* has X failed for some reason? DECLARED: xInit.c */
extern Bool xActive;

/* what events have happened    DECLARED: xEvent.c */
extern Bool xExpose, xKeyPress, xButtonPress, xMap;

/* colour information    DECLARED: xInit.c */
extern long xPixels[MAXCOLS];
extern Bool allocated[MAXCOLS];
extern Colormap colourmap;

/* cursor declarations DECLARED xInit.c */
extern Cursor normalcursor, inputcursor;

/* global pointer position and keybuffer DECLARED xEvent.c */
extern unsigned int xx, xy;
extern char keybuffer[];
extern int xButton;

XPoint points[MAXPOINTS];
XPoint tpoints[MAXPOINTS];

/* tracing    DECLARED: xTrace.c */
extern FILE *graphictraceFILE;
extern Bool graphictrace;
char graphicversion[] = "GHOST X11 postprocessor Version 2.3 Graphics \
Interface";
/*
*
*
*
*        CALL XPOLYL(IPOLYX,IPOLYY,INPOLY)
*
*        Called by: G1TRAN 
*
*        Draw a polyline (GRID function)
*
*
*
*/

int XPOLYL(xpoints, ypoints, numpoints)
int *xpoints, *ypoints, *numpoints;
{
    XPoint *p, *tp;
    int i;
    int *xp, *yp;
    unsigned int replaydata;
    register int ytemp;
    
    if (!xActive) return;    /* return if X is not active */

    xp = xpoints;
    yp = ypoints;
    p = points;
    tp = tpoints;

    for (i=0; i<*numpoints; i++)
    {

       tp->x = (short) (Dresx * (*xp) / origDresx);
       ytemp = Dresy * (*yp) / origDresy;
       tp->y = (short) (RealDresy - ytemp);
                                              /* have to DRAW vectors with */
                                              /* transform to new resized  */
                                              /* window   --               */
       p->x = (short) *xp;                    /* BUT RECORD UNtransformed  */
       p->y = (short) *yp;                    /* version                   */

       ++p;
       ++tp;
       ++xp;
       ++yp;
    }

    sigio_off();
    replaydata = XDRAWLINES;
    RecordValues(&replaydata,1);
    RecordValues(points, *numpoints);  /* XPoint is 4 bytes */

    XDrawLines(GHdisplay, GHwin, GHgc, tpoints, *numpoints, CoordModeOrigin);
#ifdef DEBUG_SYNC
    XFlush(GHdisplay);
#endif
    sigio_on();
}


/*
*
*
*
*        CALL XPOLYM(IPOLYX,IPOLYY,1)
*
*        Called by: G1TRAN
*
*        Draw a polymarker (GRID function)
*
*
*/

int XPOLYM(xpoints, ypoints, numpoints)
int *xpoints, *ypoints, *numpoints;
{
    XPoint *p, *tp;
    int i;
    int *xp, *yp;
    unsigned int replaydata;
    register int ytemp;
    
    if (!xActive) return;    /* return if X is not active */

    xp = xpoints;
    yp = ypoints;
    p = points;
    tp = tpoints;

    for (i=0; i<*numpoints; i++)
    {

       tp->x = (short) (Dresx * (*xp) / origDresx);
       ytemp = Dresy * (*yp) / origDresy;
       tp->y = (short) (RealDresy - ytemp);
                                              /* have to DRAW markers with */
                                              /* transform to new resized  */
                                              /* window   --               */
       p->x = (short) *xp;                    /* BUT RECORD UNtransformed  */
       p->y = (short) *yp;                    /* version                   */

       ++p;
       ++tp;
       ++xp;
       ++yp;
    }

    sigio_off();
    replaydata = XDRAWPOINTS;
    RecordValues(&replaydata,1);
    RecordValues(points, *numpoints);

    XDrawPoints(GHdisplay, GHwin, GHgc, tpoints, *numpoints, CoordModeOrigin);
#ifdef DEBUG_SYNC
    XFlush(GHdisplay);
#endif
    sigio_on();
}


/*
*
*
*
*        CALL XTEXT(IPOLYX(1),IPOLYY(1),CHRSTR)
*
*        Called by: G1TRAN
*
*        Draw a single character a the point specified (GRID function)
*
*
*/

int XTEXT(xpoint, ypoint, character)
int *xpoint, *ypoint;
unsigned int *character;
{
    unsigned int replaydata;
    int x,y;
    register int ytemp;
    unsigned int shifted_char;   /* shifted because character stored in low byte */

    if (!xActive) return;        /* return if X is not active */
    shifted_char = *character << 24;

    sigio_off();

    replaydata = XDRAWSTRING;
    RecordValues(&replaydata,1);

    x = Dresx * (*xpoint)/ origDresx;
    ytemp = Dresy * (*ypoint) / origDresy;
    y = RealDresy - ytemp;

    RecordValues(xpoint,1);
    RecordValues(ypoint,1);
    RecordValues((unsigned int *) &shifted_char,1);

    XDrawString(GHdisplay, GHwin, GHgc, x, y, &shifted_char, 1);

#ifdef DEBUG_SYNC
    XFlush(GHdisplay);
#endif
    sigio_on();
}


/*
*
*
*
*        CALL XCHARH(IYVAL)
*
*        ** NOT IMPLEMENTED **
*
*
*
*
*/

void XCHARH(height)
int *height;
{
    if (!xActive) return;    /* return if X is not active */

    return;
}

/*
*
*
*
*        CALL XCHARO(CHUPY,-CHUPX,CHUPX,CHUPY)
*
*        ** NOT IMPLEMENTED **
*
*
*
*
*/

void XCHARO(a,b,c,d)
int *a, *b, *c, *d;
{
    if (!xActive) return;    /* return if X is not active */

    return;
}

/*
*
*
*
*        CALL XCHARB(CHEXPF)
*
*        ** NOT IMPLEMENTED **
*
*
*
*
*/

void XCHARB(oblate)
int *oblate;
{
    if (!xActive) return;    /* return if X is not active */

    return;
}

/*
*
*
*
*        CALL XLINET(IYVAL)
*
*        ** NOT IMPLEMENTED **
*
*
*
*
*
*/

void XLINET(type)
int *type;
{
    if (!xActive) return;    /* return if X is not active */

    return;
}

/*
*
*
*
*        CALL XFILLP(IPOLYX,IPOLYY,INPOLY)
*
*
*        Called by: G1TRAN
*
*        Fill polygon defined by closed list of points (GRID function)
*
*/

int XFILLP(xpoints, ypoints, numpoints)
int *xpoints, *ypoints, *numpoints;
{
    XPoint *p, *tp;
    int i;
    int *xp, *yp;
    unsigned int replaydata;
    register int ytemp;
    
    if (!xActive) return;    /* return if X is not active */

    xp = xpoints;
    yp = ypoints;
    p = points;
    tp = tpoints;

    for (i=0; i<*numpoints; i++)
    {

       tp->x = (short) (Dresx * (*xp) / origDresx);
       ytemp = Dresy * (*yp) / origDresy;
       tp->y = (short) (RealDresy - ytemp);
                                              /* have to DRAW vectors with */
                                              /* transform to new resized  */
                                              /* window   --               */
       p->x = (short) *xp;                    /* BUT RECORD UNtransformed  */
       p->y = (short) *yp;                    /* version                   */

       ++p;
       ++tp;
       ++xp;
       ++yp;
    }


    sigio_off();
    replaydata = XFILLPOLY;
    RecordValues(&replaydata,1);
    RecordValues(points, *numpoints);  /* XPoint is 4 bytes */

    XFillPolygon(GHdisplay, GHwin, GHgc, tpoints, *numpoints,
                    Complex,CoordModeOrigin);
#ifdef DEBUG_SYNC
    XFlush(GHdisplay);
#endif
    sigio_on();
}


/*
*
*
*
*        CALL XCOLTB(KOLNUM,IREDIN,IGRNIN,IBLUIN)
*
*        Called by: G1TRAN G1HRDW
*
*        Is not called if monochrome.
*
*        Set a colour table entry. If GHOST colour number has already
*        been set then the values are simply changed. If not a new
*        X Colour Cell is allocated from the PseudoColour colourmap
*        and the values set.
*
*
*/

void XCOLTB(colournumber, red, green, blue)
int *colournumber;
int *red, *green, *blue;
{
    XColor xColour;
    unsigned long plane_masks; /* UNUSED */

    if (!xActive) return;    /* return if X is not active */
#ifdef DEBUG
    if (graphictrace)
    {
       fprintf(graphictraceFILE, "Set colour table entry\n"); 
       fflush(graphictraceFILE);
    }
#endif
    sigio_off();
    if (allocated[*colournumber] == False)
    {

#ifdef DEBUG
       if (graphictrace)
       {
          fprintf(graphictraceFILE, "Colour: alloc new cell ");
          fprintf(graphictraceFILE, "colournumber=%d\n",*colournumber);
          fflush(graphictraceFILE);
       }
#endif
       XAllocColorCells(GHdisplay, colourmap, False , &plane_masks, 0,
                           &xColour, 1);
       xPixels[*colournumber] = xColour.pixel;
       allocated[*colournumber] = True;
    }
    else
    {
       xColour.pixel = xPixels[*colournumber];
    }

#ifdef DEBUG
   if (graphictrace)
   {
      fprintf(graphictraceFILE, "Colour: store colour values ");
      fprintf(graphictraceFILE, "colournumber=%d\n",*colournumber);
      fprintf(graphictraceFILE, "colour pixel number=%d\n",xColour.pixel);
      fprintf(graphictraceFILE, "red=%d ",*red);
      fprintf(graphictraceFILE, "green=%d ",*green);
      fprintf(graphictraceFILE, "blue=%d\n",*blue);
      fflush(graphictraceFILE);
   }
#endif

    xColour.red   = (unsigned short) (65535 * (*red) / 255);
    xColour.green = (unsigned short) (65535 * (*green) / 255);
    xColour.blue  = (unsigned short) (65535 * (*blue) / 255);
    xColour.flags = DoRed | DoGreen | DoBlue;

#ifdef DEBUG
    if (graphictrace)
    {
       fprintf(graphictraceFILE, "Store colour table entry\n"); 
       fflush(graphictraceFILE);
    }
#endif

    XStoreColors(GHdisplay, colourmap, &xColour, 1);
    XFlush(GHdisplay);

    sigio_on();
    return;
}

/*
*
*
*
*        CALL XSETFC(IYVAL)
*
*        Called by: G1TRAN G1HRDW
*
*        Sets the foreground colour in the GC.
*
*
*/

void XSETFC(colour)
int *colour;
{
    unsigned int replaydata;

    if (!xActive) return;    /* return if X is not active */
    
    sigio_off();
    replaydata = XSETCOLOUR;
    RecordValues(&replaydata,1);

    RecordValues(&xPixels[*colour], 1);
    XSetForeground(GHdisplay, GHgc, xPixels[*colour]);


#ifdef DEBUG
    if (graphictrace)
    {
       fprintf(graphictraceFILE, "Set to colour number=%d\n",*colour);
       fprintf(graphictraceFILE, "colour pixel number=%d\n",xPixels[*colour]);
       fflush(graphictraceFILE);
    }
#endif

    sigio_on();

    return;
}
/*
*
*
*
*        CALL XGETIN(IX,IY,KCHAR)
*
*        Called by: G1TRAN
*
*        Change the cursor to the input prompt and wait for
*        either a mouse keyclick or keyboard event.
*
*
*/

void XGETIN(xp, yp, character)
unsigned int *character;
int *xp, *yp;

{
    XEvent event;
    if (!xActive) return;    /* return if X is not active */

    sigio_off();

    XDefineCursor(GHdisplay, GHwin, inputcursor);
    XFlush(GHdisplay);
    xButtonPress = False;
    xKeyPress = False;

    XBell(GHdisplay,100);
    while (!xButtonPress && !xKeyPress)
    {
       XPeekEvent(GHdisplay, &event);      /* see if theres an event */
       xeventhandler();
    }

    /* get back the coords of the pointer in the space of the original
       window (ghost doesnt know we've been resized) */

    *xp = (int)(origDresx * (int)xx / Dresx);
    *yp = (int)(origDresy- (origDresy*( (int)xy - RealDresy + Dresy)/Dresy));

    if (xButtonPress)
    {
       switch (xButton)
       {
          case 1:
               *character = 'l';
               break;

          case 2:
               *character = 'm';
               break;

          case 3:
               *character = 'r';
               break;

          default:
               *character = 'U';
               break;
       }
    }
    else
    {
       /* its a keypress */
       *character = keybuffer[0];
    }

#ifdef DEBUG
    if (graphictrace)
    {
       fprintf(graphictraceFILE, "xgetin: ");
       fprintf(graphictraceFILE, "cursor x=%d ",xx);
       fprintf(graphictraceFILE, "cursor y=%d\n",xy);
       fprintf(graphictraceFILE, "character= %c (0x%x)\n",
                                  *character,*character);
       fflush(graphictraceFILE);
    }
#endif

    XDefineCursor(GHdisplay, GHwin, normalcursor);
    XFlush(GHdisplay);
    sigio_on();
    return;
}

/*
*
*
*
*        CALL XCLRWN()
*
*        Called by: G1TRAN
*
*        Clears the GHOST window.
*
*
*/

void XCLRWN()
{
    if (!xActive) return;    /* return if X is not active */

    XClearWindow(GHdisplay, GHwin);
    EraseReplay();
    return;
}

/*
*
*
*
*        CALL XWNDIM(x,y)
*
*        Used by GHOST Rel. 8 for change in DRESX/Y after FRAME.
*
*
*
*
*/

void XWNDIM(winx,winy)
int *winx, *winy;

{
    XWindowAttributes windowattrs;

    *winx = 0;
    *winy = 0;
    if (!xActive) return;    /* return if X is not active */

    sigio_off();

    XGetWindowAttributes(GHdisplay, GHwin, &windowattrs);
    *winx = windowattrs.width - 1;
    *winy = windowattrs.height - 1;
    Dresx = *winx;
    Dresy = *winy;
    RealDresy = Dresy;

    if (keep_ratio)
    {
       Dresx = Min(*winx, *winy);
       Dresy = Dresx;
    }

    origDresx = Dresx;
    origDresy = Dresy;

    sigio_on();
    return;
}

/*
*
*
*
*        CALL XCLOSE
*
*        Called by: G1TRAN
*
*        Closedown, free up Server Resources and free allocated
*        replay memory.
*
*
*/

void XCLOSE()
{
    if (!xActive) return;    /* return if X is not active */

    sigio_off();
    XUnloadFont(GHdisplay, font_info->fid);
    XFreeGC(GHdisplay, GHgc);
    EraseReplay();
    XCloseDisplay(GHdisplay);
    return;
}

/*
*
*
*
*        CALL XFLSBF
*
*
*        Called by: G1TRAN
*
*        Flushes *output* only.
*
*
*/

void XFLSBF()
{
    if (!xActive) return;    /* return if X is not active */

    XFlush(GHdisplay);

    return;
}
