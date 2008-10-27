/*
*
*   Initialise, switch on and off, signal handler with
*   XEvent loop.
*
*   Look at eventversion below for version number and
*   in the file VERSION for history.
*
*
*/

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>

#include <stdio.h>
#include <signal.h>
#include <fcntl.h>

#include "opcodes.h"
#include "defines.h"

#ifdef streams
#include <stropts.h>
#endif

extern Display *GHdisplay;
extern Window GHwin;
extern GC GHgc;

/* which events occured during event loop. */
Bool xExpose, xKeyPress, xButtonPress, xMap;

/* global pointer position and keybuffer */
unsigned int xx, xy;
char keybuffer[1];
int xButton;

/* current, original and real y window addressability (size) DECLARED: xInit.c */
extern int Dresx, Dresy;
extern int origDresx, origDresy;
extern int RealDresy;
static int newDresx,newDresy;
/* user default items DECLARED xInit.c*/
extern Bool keep_ratio;

XEvent report;
int fflags;
Bool IO_on = False;

/* replay workspace  MAYBE malloc later */
XPoint replay[MAXPOINTS];

extern FILE *eventtraceFILE;
extern Bool eventtrace;
char eventversion[] = "GHOST X11 postprocessor Version 2.3 Event Handler";
xeventhandler()
{
    XComposeStatus status; /* UNUSED */
    KeySym keysym;         /* UNUSED */

    xMap = False;
    xExpose = False;
    xKeyPress = False;
    xButtonPress = False;
#ifdef DEBUG_LOW
    if (eventtrace)
    {
       fprintf(eventtraceFILE,"event handler\n");
       fflush(eventtraceFILE);
    }
#endif

#ifdef SysV
    signal(SIGIO, SIG_IGN);
#endif

    while (XPending(GHdisplay))
    {
#ifdef DEBUG_LOW
       if (eventtrace)
       {
          fprintf(eventtraceFILE,"Events waiting...\n");
          fflush(eventtraceFILE);
       }
#endif
       XNextEvent(GHdisplay, &report);
#ifdef DEBUG_LOW
       if (eventtrace)
       {
          fprintf(eventtraceFILE,"...got next event\n");
          fflush(eventtraceFILE);
       }
#endif
       switch (report.type)
       {
       case MapNotify:     /* -------------------------------*/
            /* get rid of all Expose events on the queue */
            while (XCheckTypedEvent(GHdisplay, Expose, &report));
            xMap = True;
#ifdef DEBUG
            if (eventtrace)
            {
               fprintf(eventtraceFILE,"MapNotifyEvent\n");
               fflush(eventtraceFILE);
            }
#endif
	    
	    break;
       case Expose:     /* -------------------------------*/
            /* get rid of all other Expose events on the queue */
            while (XCheckTypedEvent(GHdisplay, Expose, &report));
            xExpose = True;
#ifdef DEBUG
            if (eventtrace)
            {
               fprintf(eventtraceFILE,"ExposeEvent\n");
               fflush(eventtraceFILE);
            }
#endif
	    
	    ReplayGraphics();
	    break;
       case ConfigureNotify: /* --------------------------*/
            /* window has been resized, change width and height to
             * send to draw_text and draw_graphics in next Expose */

            newDresx = report.xconfigure.width - 1;
            newDresy = report.xconfigure.height - 1;
            RealDresy = newDresy;

            if (keep_ratio)
            {
               newDresx = Min(newDresx, newDresy);
               newDresy = Min(newDresx, newDresy);
            }
	    if(newDresx != Dresx || newDresy != Dresy){
	      Dresx = newDresx;
	      Dresy = newDresy;
	      /* Force an expose event */
	      XClearArea(GHdisplay,GHwin,0,0,0,0,True);     
#ifdef DEBUG
	      if (eventtrace)
		{
		  fprintf(eventtraceFILE,"ConfigureEvent\n");
		  fprintf(eventtraceFILE,"Dresx=%d Dresy=%d\n",Dresx,Dresy);
		  fflush(eventtraceFILE);
		}
#endif
	    }

            break;
       case ButtonPress: /* -------------------------------*/

            /* get button number and position */

            xButtonPress = True;
            xButton = report.xbutton.button;
            xx = report.xbutton.x;
            xy = report.xbutton.y;

#ifdef DEBUG
            if (eventtrace)
            {
               fprintf(eventtraceFILE,"ButtonPress\n");
               fflush(eventtraceFILE);
            }
#endif
            break;
       case KeyPress: /* ----------------------------------*/

            /* set the input flag for GHOST. etc.. */
            xKeyPress = True;
            XLookupString(&report, keybuffer, 1, &keysym, &status);
            xx = report.xbutton.x;
            xy = report.xbutton.y;

#ifdef DEBUG
           if (eventtrace)
           {
              fprintf(eventtraceFILE,"KeyPress\n");
              fflush(eventtraceFILE);
           }
#endif

           break;
       default:      /* ----------------------------------*/

            /* all events selected by StructureNotifyMask
             * except ConfigureNotify are thrown away here,
             * since nothing is done with them. */
#ifdef DEBUG
            if (eventtrace)
            {
               fprintf(eventtraceFILE,"Untrapped XEvent\n");
               fflush(eventtraceFILE);
            }
#endif
            break;
       } /* end switch */

    } /* end while */

#ifdef DEBUG_LOW
    if (eventtrace)
    {
       fprintf(eventtraceFILE,"return from event loop.\n");
       fflush(eventtraceFILE);
    }
#endif

    return;
}

sigio_init()
{
    int fflags;

#ifdef SysV
    fflags = fcntl( GHdisplay->fd, O_NONBLOCK, 0);
#else
    fflags = fcntl( GHdisplay->fd, F_GETFL, 0);
    fflags = fcntl( GHdisplay->fd, F_SETFL, fflags | FASYNC);
    fflags = fcntl( GHdisplay->fd, F_SETOWN, getpid());
#endif

#ifdef streams
    (void) ioctl(GHdisplay->fd, I_SETSIG, S_INPUT | S_OUTPUT);
#endif

#ifdef DEBUG
    if (eventtrace)
    {
       fprintf(eventtraceFILE,"sigio_init\n");
       fflush(eventtraceFILE);
    }
#endif

}

sigio_on()
{
    if (IO_on) return;
    signal(SIGIO, xeventhandler);
    IO_on = True;

#ifdef DEBUG_LOW
    if (eventtrace)
    {
       fprintf(eventtraceFILE,"sigio_on\n");
       fflush(eventtraceFILE);
    }
#endif
}

sigio_off()
{
    if (!IO_on) return;
    signal(SIGIO, SIG_IGN);
    IO_on = False;

#ifdef DEBUG_LOW
    if (eventtrace)
    {
       fprintf(eventtraceFILE,"sigio_off\n");
       fflush(eventtraceFILE);
    }
#endif
}


ReplayGraphics()
{
    int len;               /* length of data fron ReplayValues          */
    int i;                 /* temporary used in loop                    */
    unsigned int *ptr;     /* pick up address of replay to get opcode.  */
                           /* also used for replay ptr                  */
    XPoint *rp;            /* pointer to replay for scaling loop        */
    unsigned long pixel;   /* the pixel used to call X                  */
    int x, y;              /* position of DrawString                    */
    char *character;       /* DrawString                                */
    register int ytemp;    /* temporary for transforming endpoints      */

#ifdef DEBUG
    if (eventtrace)
    {
       fprintf(eventtraceFILE,"ReplayGraphics\n");
       fflush(eventtraceFILE);
    }
#endif

/*    XClearWindow(GHdisplay,GHwin); */

    while ((len = ReplayValues(replay)) != 0)
    {
       ptr = (unsigned int *) replay;
       switch (*ptr)
       {
          case XDRAWLINES:   /*---------------------------------------*/

#ifdef DEBUG_LOW
               if (eventtrace)
               {
                  fprintf(eventtraceFILE,"Drawlines\n");
                  fflush(eventtraceFILE);
               }
#endif

               if ((len = ReplayValues(replay)) == 0) break;
#ifdef DEBUG_LOW
               if (eventtrace)
               {
                  fprintf(eventtraceFILE,"Len= %d\n", len);
                  fflush(eventtraceFILE);
               }
#endif

               rp = replay;
               for (i=0; i<len; i++)
               {
                  rp->x = (short) (Dresx * (int)(rp->x) / origDresx);
                  ytemp = Dresy * (int)(rp->y) / origDresy;
                  rp->y = (short) (RealDresy - ytemp);
                  ++rp;
               }
               
               XDrawLines(GHdisplay, GHwin, GHgc, replay, len, CoordModeOrigin);

#ifdef DEBUG_SYNC
               XFlush(GHdisplay);
#endif

               break;
          case XDRAWPOINTS:  /*---------------------------------------*/

#ifdef DEBUG_LOW
               if (eventtrace)
               {
                  fprintf(eventtraceFILE,"Drawpoints\n");
                  fflush(eventtraceFILE);
               }
#endif
               if ((len = ReplayValues(replay)) == 0) break;

               rp = replay;
               for (i=0; i<len; i++)
               {
                  rp->x = (short) (Dresx * (int)(rp->x) / origDresx);
                  ytemp = Dresy * (int)(rp->y) / origDresy;
                  rp->y = (short) (RealDresy - ytemp);
                  ++rp;
               }

               XDrawPoints(GHdisplay, GHwin, GHgc, replay, len, CoordModeOrigin);
#ifdef DEBUG_SYNC
               XFlush(GHdisplay);
#endif
               break;
          case XSETCOLOUR:      /*---------------------------------------*/

#ifdef DEBUG_LOW
               if (eventtrace)
               {
                  fprintf(eventtraceFILE,"Setcolour\n");
                  fflush(eventtraceFILE);
               }
#endif
               if ((len = ReplayValues(&pixel)) == 0) break;
#ifdef DEBUG
               if (eventtrace)
               {
                  fprintf(eventtraceFILE,"Setcolour pixel=%d\n",pixel);
                  fflush(eventtraceFILE);
               }
#endif
               XSetForeground(GHdisplay, GHgc, pixel);

#ifdef DEBUG_SYNC
               XFlush(GHdisplay);
#endif
               break;

          case XFILLPOLY:      /*---------------------------------------*/

#ifdef DEBUG_LOW
               if (eventtrace)
               {
                  fprintf(eventtraceFILE,"FillPoly\n");
                  fflush(eventtraceFILE);
               }
#endif
               if ((len = ReplayValues(replay)) == 0) break;
#ifdef DEBUG_LOW
               if (eventtrace)
               {
                  fprintf(eventtraceFILE,"Len= %d\n", len);
                  fflush(eventtraceFILE);
               }
#endif
               rp = replay;
               for (i=0; i<len; i++)
               {
                  rp->x = (short) (Dresx * (int)(rp->x) / origDresx);
                  ytemp = Dresy * (int)(rp->y) / origDresy;
                  rp->y = (short) (RealDresy - ytemp);
                  ++rp;
               }
               
               XFillPolygon(GHdisplay, GHwin, GHgc, replay, len, Complex,
                               CoordModeOrigin);
#ifdef DEBUG_SYNC
               XFlush(GHdisplay);
#endif
               break;

          case XDRAWSTRING:      /*---------------------------------------*/

#ifdef DEBUG_LOW
               if (eventtrace)
               {
                  fprintf(eventtraceFILE,"Drawtext\n");
                  fflush(eventtraceFILE);
               }
#endif
               if ((len = ReplayValues(replay)) == 0) break;
               ptr = (unsigned int *) replay;
               x = Dresx * (int) *ptr / origDresx;

               if ((len = ReplayValues(replay)) == 0) break;
               ptr = (unsigned int *) replay;
               y = Dresy * (int) *ptr / origDresy;
               y = RealDresy - y;

               if ((len = ReplayValues(replay)) == 0) break;
               ptr = (unsigned int *) replay;
               character = (char *) ptr;

               XDrawString(GHdisplay, GHwin, GHgc, x, y, character, 1);

#ifdef DEBUG_SYNC
               XFlush(GHdisplay);
#endif
               break;


          default:
#ifdef DEBUG
               if (eventtrace)
               {
                  fprintf(eventtraceFILE,"Unknown Graphics Type in Replay!\n");
                  fflush(eventtraceFILE);
               }
#endif
               ;

       } /* end switch */
    } /* end while */
    XFlush(GHdisplay);    /* make sure redraw complete */
}
