/*   xInit.c
 *
 *   Royd Whittington, SERC Daresbury Laboratory 1989.
 *
 *   X11 Initialisation is done here for GHOST80 postprocessor.
 *
 *   Look at initversion below for version number and
 *   in the file VERSION for history.
 *
 *
 */
/**IGDS********************************************************
 * Modded IGDS 26-4-91 to incorporate command line args
 * All IGDS mods are flagged as this comment
 **************************************************************/
#include <stdio.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/cursorfont.h>

#include "defines.h"

#include "bitmaps/icon_bitmap"
#define BITMAPDEPTH 1

/**IGDS******************************************************************
 * Define function codes as per GHOST manual
 ************************************************************************/
#define GHT  0    /* Window height                                      */
#define GWID 1    /* Window width                                       */
#define GXOR 2    /* X-coordinate of window origin on screen            */
#define GYOR 3    /* Y-coordinate of window origin on screen            */
#define GTIT 4    /* Window title (appears in Window manager title bar) */


Display *GHdisplay;
int GHscreen;
Window GHwin;
GC GHgc;
XFontStruct *font_info;
Bool xActive = False;

/* what events have happened    DECLARED: xEvent.c */
extern Bool xExpose, xKeyPress, xButtonPress, xMap;

/* current, old and real y (for keep aspect ratio) window addressability (size) */
int Dresx, Dresy;
int origDresx, origDresy;
int RealDresy;

/* colour information */
long xPixels[MAXCOLS];
Bool allocated[MAXCOLS];
Colormap colourmap;

/* cursor declarations */
Cursor normalcursor, inputcursor;

/* user default items */
Bool keep_ratio = True;

/*tracing    DECLARED: xTrace.c */
extern FILE *inittraceFILE;
extern Bool inittrace;
char ProgName[256]; /**IGDS program name passed from FORTRAN **/
char initversion[] = "GHOST X11 postprocessor Version 2.3 Initialisation";
/**IGDS Modded argument declaration *******/
XINIT(flags,comdata,depth, dresx, dresy,win_name,prog_name,win_namlen,prog_namlen)
int *flags,*comdata,*depth, *dresx, *dresy;
char *win_name,*prog_name;  /* NB FORTRAN strings do not end in \0 use next arg for length */
int win_namlen,prog_namlen;      /* NB this is passed as a 'hidden argument' by FORTRAN */

{
    char *myGetResource();
    unsigned int width, height;           /* window size               */
    int x = 0, y = 0;                     /* window position           */
    unsigned int border_width = 4;        /* border four pixels wide   */
    int PositionHint =0;              /* Does WM intercept ?       */
    unsigned int display_width, display_height;
    char window_name[256];         /**IGDS set default in code  */
    char *icon_name = "GHOST80";
    Pixmap icon_pixmap;
    XSizeHints size_hints;

    /* General Colour Setup */
    XVisualInfo vTemplate;                /* template of the >8 visual  */
    XVisualInfo *visualList;              /* list of visuals that match */
    int visualsMatched;                   /* and the number that match  */
    Visual *defaultVisual;                /* used to contain default    */
    XColor BorderColour, ColourDummy;
    XColor xColour;
    unsigned long plane_masks;

    XSetWindowAttributes attributes;      /* for XCreateWindow */
    unsigned long mask;

    char *defaultstring;
    char *display_name = NULL;
    XWindowAttributes win_attributes;
    XEvent event;      /* dummy return for wait for initial expose loop */
    int i;                                /* loop variable              */

    char **myargv=NULL;     /* dummy for XSetProperties */
    char myargc=0;


#ifdef DEBUG
    DEBUG_init();  /* if debugging initialise. */
#endif

    /* connect to X server */

    xActive = False;
    /**IGDS copy program name to the global array ProgName ***/
    for(i=prog_namlen-1;i>=0;i--)if(prog_name[i] != ' ')break;
    ProgName[i+1]=0;
    for(;i>=0;i--)ProgName[i] = prog_name[i];
    if (( GHdisplay = XOpenDisplay(display_name)) == NULL )
    {
       (void) fprintf(stderr,"GHOST80 X: ");
       (void) fprintf(stderr,"cannot connect to X server %s\n",
                        XDisplayName(display_name));
       (void) fprintf(stderr,"Check \"DISPLAY\" environment variable.\n");
       return 1;
    }

    xActive = True;

    /* get screen size from display structure macro  */
    GHscreen = DefaultScreen(GHdisplay);

#ifdef DEBUG
    if (inittrace)
    {
       fprintf(inittraceFILE, "Server name: %s\n", XDisplayName(display_name));
       fprintf(inittraceFILE, "Vendor: %s         Server Release: %d\n",
               ServerVendor(GHdisplay), VendorRelease(GHdisplay));
       fprintf(inittraceFILE, "Using X Version %d Release %d\n",
               ProtocolVersion(GHdisplay),ProtocolRevision(GHdisplay));
       fflush(inittraceFILE);
    }
#endif

    display_width  = DisplayWidth(GHdisplay, GHscreen);
    display_height = DisplayHeight(GHdisplay, GHscreen);
    *depth = (int) DefaultDepth(GHdisplay,GHscreen);

    /* default window size for GHOST */
    height = display_width/2.5;
    width  = (int) ((float)height * 1.4);

    /* change that if user has his own defaults */
    defaultstring = (char *) myGetResource( "gwid");
    if (defaultstring != NULL)
       if (atoi(defaultstring) >= 0)
       {  
          width = atoi(defaultstring) < display_width ? atoi(defaultstring) :
                                                     display_width;
       }
    defaultstring = (char *) myGetResource( "ght");
    if (defaultstring != NULL)
       if (atoi(defaultstring) >= 0)
       {
          height = atoi(defaultstring) < display_height ? atoi(defaultstring) :
                                                        display_height;
       }
    /**IGDS  Allow user to specify x and y origin as well **/
    defaultstring = (char *)myGetResource( "gxor" );
    if(defaultstring != NULL )
	if(atoi(defaultstring) >= 0)
	{
          PositionHint=1;
          x = atoi(defaultstring) < display_width - width ? atoi(defaultstring) :
							display_width-width;
	}

    defaultstring = (char *)myGetResource( "gyor" );
    if(defaultstring != NULL )
	if(atoi(defaultstring) >= 0)
	{
          PositionHint=1;
          y = atoi(defaultstring) < display_height - height ? atoi(defaultstring) :
							display_height-height;
	}

    /**IGDS now do command line arguments which override Xresource defaults **/
    if(flags[GHT] == 1)height = comdata[GHT] < display_height ? comdata[GHT] : display_height;
    if(flags[GWID] == 1)width = comdata[GWID] < display_width ? comdata[GWID] : display_width;
    if(flags[GXOR] == 1){
          PositionHint=1;
          x = comdata[GXOR] < display_width -width ? comdata[GXOR] : display_width-width;
	}
    if(flags[GYOR] == 1){
         PositionHint=1;
         y = comdata[GYOR] < display_height-height ? comdata[GYOR] : display_height-height;
       }      
    /* set up colour system: default to no allocated cells
       and create screen colourmap */
    for (i=0; i<MAXCOLS; i++) allocated[i]=False;

    switch (*depth) {
       case 1:
#ifdef DEBUG
               if (inittrace)
               {
                  fprintf(inittraceFILE,
                          "1-bit display\n"); 
                  fflush(inittraceFILE); XFlush(GHdisplay);
               }
#endif
               BorderColour.pixel = BlackPixel(GHdisplay,GHscreen);
               for (i=1; i<MAXCOLS; i++) xPixels[i]=BlackPixel(GHdisplay,GHscreen);
               for (i=1; i<MAXCOLS; i++) allocated[i]=True;
               xPixels[0] = WhitePixel(GHdisplay,GHscreen);

               visualList = XGetVisualInfo(GHdisplay, VisualNoMask,
                                            &vTemplate, &visualsMatched);
               visualList->visual = DefaultVisual(GHdisplay,GHscreen);
               vTemplate.depth = 1;
               colourmap = DefaultColormap(GHdisplay, GHscreen);

               break;

       case 2:
       case 3:
       case 4:
       case 5:
       case 6:
       case 7:
#ifdef DEBUG
               if (inittrace)
               {
                  fprintf(inittraceFILE,
                          "<8 bit display\n");
                  fflush(inittraceFILE); XFlush(GHdisplay);
               }
#endif

                  vTemplate.screen = GHscreen;     /* this is what we want */
                  vTemplate.depth  = *depth;
                  vTemplate.class  = PseudoColor;  /* with pseudocolour    */
                  visualList = XGetVisualInfo(GHdisplay,
                     VisualScreenMask | VisualDepthMask | VisualClassMask, 
                     &vTemplate, &visualsMatched);

                  if (visualsMatched == 0)
                  {
#ifdef DEBUG
                     if (inittrace)
                     {
                        fprintf(inittraceFILE,
                                "cannot get pseudocolour visual "); 
                        fprintf(inittraceFILE,
                                "from <8 bit display - using DefaultColormap");
                        fflush(inittraceFILE); XFlush(GHdisplay);
                     }
#endif
                    colourmap = DefaultColormap(GHdisplay, GHscreen);
                    break;
                  }

                  /* create a colourmap with the first matching visual */
                  colourmap = XCreateColormap(GHdisplay,
                                 RootWindow(GHdisplay,GHscreen),
                                 visualList[0].visual,
                                 AllocNone);

                  XInstallColormap(GHdisplay, colourmap);

               defaultstring = (char *) myGetResource( "borderColor");
               XAllocNamedColor(GHdisplay, colourmap, defaultstring,
                                   &BorderColour, &ColourDummy);

               /* allocate colour for background and call it GHOST colour zero*/

               if (XAllocColorCells(GHdisplay, colourmap, False, &plane_masks, 0,
                        &xColour, 1) == False)
               {
#ifdef DEBUG
                  if (inittrace)
                  {
                     fprintf(inittraceFILE, "Can't alloc colour cell\n");
                     fflush(inittraceFILE); XFlush(GHdisplay);
                     return 1;
                  }
#endif
               }
               xPixels[0] = xColour.pixel;
               allocated[0] = True;

               xColour.red   = (unsigned short) 0;
               xColour.green = (unsigned short) 0; /*default black background*/
               xColour.blue  = (unsigned short) 0;
               xColour.flags = DoRed | DoGreen | DoBlue;
#ifdef DEBUG
               if (inittrace)
               {
                  fprintf(inittraceFILE, "XStoreColors for background...\n"); 
                  fflush(inittraceFILE); XFlush(GHdisplay);
               }
#endif

               XStoreColors(GHdisplay, colourmap, &xColour, 1);

               break;

       default:
               /* procedure is to determine whether the default visual is
                * class PseudoColor and use it. If not then obtain a PseudoColor
                * visual from the server and create/install the appropriate
                * colourmap. If a PseudoColor visual does not exist then
                * the default colourmap is used (which would probably not
                * work well later on)
                */

#ifdef DEBUG
               if (inittrace)
               {
                  fprintf(inittraceFILE,
                          ">=8 bit display\n"); 
                  fflush(inittraceFILE); XFlush(GHdisplay);
               }
#endif

               defaultVisual = DefaultVisual(GHdisplay, GHscreen);
               if (defaultVisual->class == PseudoColor)
	       {
#ifdef DEBUG
                  if (inittrace)
                  {
                     XSync(GHdisplay, False);
                     fprintf(inittraceFILE,
                             "Default Visual class = PseudoColor\n"); 
                     fflush(inittraceFILE);
                  }
#endif
                  vTemplate.screen = GHscreen;     /* this is what we got  */
                  vTemplate.depth  = 8;            /* eight bit,           */
                  vTemplate.class  = PseudoColor;  /* with pseudocolour    */
                  visualList = XGetVisualInfo(GHdisplay, VisualNoMask,
                                               &vTemplate, &visualsMatched);
                  visualList->visual = DefaultVisual(GHdisplay,GHscreen);
                  colourmap = DefaultColormap(GHdisplay, GHscreen);
	       }
               else
	       {
#ifdef DEBUG
                     if (inittrace)
                     {
                        XSync(GHdisplay,False);
                        fprintf(inittraceFILE,
                                "Matching Visual...\n"); 
                        fflush(inittraceFILE);
                     }
#endif
                  vTemplate.screen = GHscreen;     /* this is what we want */
                  vTemplate.depth  = 8;            /* eight bit,           */
                  vTemplate.class  = PseudoColor;  /* with pseudocolour    */
                  visualList = XGetVisualInfo(GHdisplay,
                     VisualScreenMask | VisualDepthMask | VisualClassMask, 
                     &vTemplate, &visualsMatched);

                  if (visualsMatched == 0)
                  {
#ifdef DEBUG
                     if (inittrace)
                     {
                        fprintf(inittraceFILE,
                                "cannot get pseudocolour visual "); 
                        fprintf(inittraceFILE,
                                "from >=8 bit display - using DefaultColormap");
                        fflush(inittraceFILE); XFlush(GHdisplay);
                     }
#endif
                    colourmap = DefaultColormap(GHdisplay, GHscreen);
                    break;
                  }

                  /* create a colourmap with the first matching visual */
                  colourmap = XCreateColormap(GHdisplay,
                                 RootWindow(GHdisplay,GHscreen),
                                 visualList[0].visual,
                                 AllocNone);

                  XInstallColormap(GHdisplay, colourmap);
	       }

               defaultstring = myGetResource("borderColor");
               XAllocNamedColor(GHdisplay, colourmap, defaultstring,
                                   &BorderColour, &ColourDummy);

               /* allocate colour for background and call it GHOST colour zero*/

               if (XAllocColorCells(GHdisplay, colourmap, False, &plane_masks, 0,
                        &xColour, 1) == False)
               {
#ifdef DEBUG
                  if (inittrace)
                  {
                     fprintf(inittraceFILE, "Can't alloc colour cell\n");
                     fflush(inittraceFILE); XFlush(GHdisplay);
                     return 1;
                  }
#endif
               }
               xPixels[0] = xColour.pixel;
               allocated[0] = True;

               xColour.red   = (unsigned short) 0;
               xColour.green = (unsigned short) 0; /*default black background*/
               xColour.blue  = (unsigned short) 0;
               xColour.flags = DoRed | DoGreen | DoBlue;
#ifdef DEBUG
               if (inittrace)
               {
                  fprintf(inittraceFILE, "XStoreColors for background...\n"); 
                  fflush(inittraceFILE); XFlush(GHdisplay);
               }
#endif

               XStoreColors(GHdisplay, colourmap, &xColour, 1);

               break;

    } /* end switch  */ 

    /* create opaque window */

    mask = CWBackPixel | CWBorderPixel | CWColormap ;
    attributes.background_pixel = xPixels[0];
    attributes.border_pixel = BorderColour.pixel;
    attributes.colormap = colourmap;
    
#ifdef DEBUG
               if (inittrace)
               {
                  XSync(GHdisplay,False);
                  fprintf(inittraceFILE, "XCreatewindow...\n"); 
                  fflush(inittraceFILE);
               }
#endif

    GHwin = XCreateWindow(GHdisplay, RootWindow(GHdisplay,GHscreen),
                x, y, width, height, border_width, vTemplate.depth,
                InputOutput, visualList[0].visual, mask, &attributes);

#ifdef DEBUG
    if (inittrace)
    {
       XSync(GHdisplay,False);
       fprintf(inittraceFILE, "Free visualList...\n"); 
       fflush(inittraceFILE);
    }
#endif
    XFree(visualList);

#ifdef DEBUG
    if (inittrace)
    {
       XSync(GHdisplay,False);
       fprintf(inittraceFILE, "Create Pixmap for icon...\n"); 
       fflush(inittraceFILE);
    }
#endif
    /* create pixmap of depth 1 (bitmap) for icon */
    icon_pixmap = XCreateBitmapFromData(GHdisplay, GHwin, icon_bitmap_bits,
                      icon_bitmap_width, icon_bitmap_height);

    /* initialise hint properties for window manager */
    /*IGDS Get min and max widths and heights from resource database **/
    defaultstring = myGetResource("gMinHeight");
    if(defaultstring){
           size_hints.min_height = atoi(defaultstring);
	 }
    else {
           size_hints.min_height = 16;
	 }
    defaultstring = myGetResource("gMinWidth");
    if(defaultstring){
           size_hints.min_width = atoi(defaultstring);
	 }
    else {
           size_hints.min_width = 16;
	 }
    defaultstring = myGetResource("gMaxHeight");
    if(defaultstring){
           size_hints.max_height = atoi(defaultstring);
	 }
    else {
           size_hints.max_height = display_height;
	 }
    defaultstring = myGetResource("gMaxWidth");
    if(defaultstring){
           size_hints.max_width = atoi(defaultstring);
	 }
    else {
           size_hints.max_width = display_width;
	 }
    size_hints.flags = PSize | PMinSize | PMaxSize;
    if(PositionHint) size_hints.flags |= (USPosition | PPosition);
    size_hints.x = x;
    size_hints.y = y;
    size_hints.width = width;
    size_hints.height = height;

    /* set properties for window manager (always before mapping) */
    /**IGDS set up window name ***/
    strcpy(window_name,"GHOST-80 X-window");
    defaultstring = myGetResource("gtit");
   
    if(defaultstring != NULL){ /* Name specified in resource database */
       strcpy(window_name,defaultstring);
     }
    /**IGDS ... now see if it's in the command line args **/
    if(flags[GTIT] == 1){
       int i;
       /* Copy FORTRAN string into window_name */
       for(i=0;i<comdata[GTIT];i++)window_name[i] = win_name[i];
       window_name[comdata[GTIT]] = '\0';
     }
    XSetStandardProperties(GHdisplay, GHwin, window_name, icon_name,
                             icon_pixmap, myargv, myargc, &size_hints);
    /* select event types wanted */
    XSelectInput(GHdisplay, GHwin, EVENTSWANTED);

#ifdef DEBUG
    if (inittrace)
    {
       XSync(GHdisplay,False);
       fprintf(inittraceFILE, "Loading font...\n"); 
       fflush(inittraceFILE);
    }
#endif
    load_font(&font_info);

#ifdef DEBUG
    if (inittrace)
    {
       fprintf(inittraceFILE, "Creating GC...\n"); 
       fflush(inittraceFILE); XFlush(GHdisplay);
    }
#endif
     /* create GC for text and drawing */
    get_GC(GHwin, &GHgc, font_info);

#ifdef DEBUG
    if (inittrace)
    {
       fprintf(inittraceFILE, "Setup cursors...\n"); 
       fflush(inittraceFILE); XFlush(GHdisplay);
    }
#endif
    /* setup cursor character */
    normalcursor = XCreateFontCursor(GHdisplay, XC_diamond_cross);
    inputcursor  = XCreateFontCursor(GHdisplay, XC_crosshair);
    XDefineCursor(GHdisplay, GHwin, normalcursor);

#ifdef DEBUG
    if (inittrace)
    {
       fprintf(inittraceFILE,"...and XmapWindow.\n"); 
       fflush(inittraceFILE); XFlush(GHdisplay);
    }
#endif
    /* GHdisplay window */
    XMapWindow(GHdisplay, GHwin);

#ifdef DEBUG
    if (inittrace)
    {
       fprintf(inittraceFILE, "End of X initialisation\n"); 
       fflush(inittraceFILE); XFlush(GHdisplay);
    }
#endif

    xMap = False;
    while (!xMap)
    {
      XPeekEvent(GHdisplay, &event);   /* see if theres an event */
      xeventhandler();
    }

    XGetWindowAttributes(GHdisplay, GHwin, &win_attributes);
    *dresx = win_attributes.width - 1;
    *dresy = win_attributes.height - 1;
    Dresx     = *dresx;                     /* current window x size    */
    Dresy     = *dresy;                     /* current window y size    */
    RealDresy = *dresy;

    if (keep_ratio)
    {
       Dresx = Min(*dresx, *dresy);
       Dresy = Dresx;
    }

    origDresx = Dresx;                      /* original window x size   */
    origDresy = Dresy;                      /* original window y size   */

    sigio_init();
    sigio_on();

#ifdef DEBUG
    if (inittrace)
    {
       fprintf(inittraceFILE, "Initialisation Complete\n");
       fprintf(inittraceFILE, " win_attributes.width=%d ",
                                            win_attributes.width);
       fprintf(inittraceFILE, " win_attributes.height=%d\n",
                                           win_attributes.height);
       fprintf(inittraceFILE, " origDresx=%d origDresy=%d\n",
                                             origDresx,origDresy); 
       fprintf(inittraceFILE, " GHdisplay depth=%d\n",*depth);
       fflush(inittraceFILE); XFlush(GHdisplay);
    }
#endif

    return;
}

get_GC(GHwin, GHgc, font_info)
Window GHwin;
GC *GHgc;
XFontStruct *font_info;
{
    unsigned long valuemask = 0;    /* ignore XGCvalues and use defaults  */
    XGCValues values;
    unsigned int line_width = 0;
    int line_style = LineSolid;
    int cap_style = CapButt;
    int join_style = JoinRound;

    /* create default graphics context */
    *GHgc = XCreateGC(GHdisplay, GHwin, valuemask, &values);

    /* specify font */
    if (font_info != (XFontStruct *) NULL)
       XSetFont(GHdisplay, *GHgc, font_info->fid);

    /* specify white foreground */
    XSetForeground(GHdisplay, *GHgc, WhitePixel(GHdisplay,GHscreen));

    /* set line attributes */
    XSetLineAttributes(GHdisplay, *GHgc, line_width, line_style, cap_style,
       join_style);

}

load_font(font_info_ptr)
XFontStruct **font_info_ptr;
{ 
    char *myGetResource();
    char *fontname = myGetResource("gfont");

    /* access font */
    if(fontname== NULL){  /** No font specified in resource database **/
       if ((*font_info_ptr = XLoadQueryFont(GHdisplay,"9x15")) == NULL)
       {
          (void) fprintf(stderr, "GHOST80-X11: cannot find 9x15 ");
          (void) fprintf(stderr, "GHOST80-X11: using GC default\n");
       }
     }
    else {
      if ((*font_info_ptr = XLoadQueryFont(GHdisplay,fontname)) == NULL)
	{
	  (void) fprintf(stderr, "GHOST80-X11: cannot open font: %s\n", fontname);
	  (void) fprintf(stderr, "GHOST80-X11: using 9x15\n");
	  
	  if ((*font_info_ptr = XLoadQueryFont(GHdisplay,"9x15")) == NULL)
	    {
	      (void) fprintf(stderr, "GHOST80-X11: cannot find 9x15 ");
	      (void) fprintf(stderr, "GHOST80-X11: using GC default\n");
	    }
	}
    }
}
/**IGDS GetResource.
 * This routine is a temporary solution to the Class/Instance differentiation in
 * Resource files, using XGetDefault. First resource specifiers are scanned for
 * entries starting with the string contained in ProgName, which is taken as the
 * "instance name". Next we look for specifiers with the class name (in this case "Ghost"),
 * if a default was not found with the instance name.
 *
 * NB this should really be done by using the Xrm routines, but these require access
 *    to all the command line arguments, and these are not all parsed by the FORTRAN.
 **/

char *myGetResource( Option ) char *Option; {
 char *result;
 result = XGetDefault(GHdisplay,ProgName,Option);
 if(result== NULL)
   result = XGetDefault(GHdisplay,"Ghost",Option);
 return(result);
}
