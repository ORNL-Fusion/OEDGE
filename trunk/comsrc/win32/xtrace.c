/*
*
*  Initialise Debug log files.
*
*  Throughout the source there are cpp ifdef's. If defined at compilation:
*
*  DEBUG      - gives the highest level of debug information.
*  DEBUG_SYNC - performs an XFlush after most graphical output operations.
*  DEBUG_LOW  - gives verbose low level debug information.
*
*  Given that, the existence of the file, named below, will cause debug
*  information to be generated and output into this file. Note that this
*  could be a ln -s /dev/tty!
*
*
*/

#include <stdio.h>
#include <time.h>

#include "defines.h"

#define REPLAYTRACEFILE "GHreplay.log"
#define EVENTTRACEFILE "GHevent.log"
#define INITTRACEFILE "GHinit.log"
#define GRAPHICTRACEFILE "GHgraphics.log"

/* trace replay buffer recording code */
FILE *replaytraceFILE;
Bool replaytrace=False;
extern char replayversion[];

/* trace event handling code */
FILE *eventtraceFILE;
Bool eventtrace=False;
extern char eventversion[];

/* trace initialisation */
FILE *inittraceFILE;
Bool inittrace=False;
extern char initversion[];

/* trace graphics interface */
FILE *graphictraceFILE;
Bool graphictrace=False;
extern char graphicversion[];

DEBUG_init()
{
    char *ctime();
    long timevar;

    replaytrace = False;
    if (fopen(REPLAYTRACEFILE,"r") != NULL)
    {
       if ((replaytraceFILE = fopen(REPLAYTRACEFILE,"a")) != NULL)
       {
          replaytrace = True;
          fprintf(replaytraceFILE,"\n%s\n",replayversion);
          timevar = time(0);
          fprintf(replaytraceFILE, "%s\n", ctime(&timevar) );
       }
    }

    eventtrace = False;
    if (fopen(EVENTTRACEFILE,"r") != NULL)
    {
       if ((eventtraceFILE = fopen(EVENTTRACEFILE,"a")) != NULL)
       {
          eventtrace = True;
          fprintf(eventtraceFILE,"\n%s\n",eventversion);
          timevar = time(0);
          fprintf(eventtraceFILE, "%s\n", ctime(&timevar) );
       } /* end if */
    } /* end if */

    inittrace = False;
    if (fopen(INITTRACEFILE,"r") != NULL)
    {
       if ((inittraceFILE = fopen(INITTRACEFILE,"a")) != NULL)
       {
          inittrace = True;
          fprintf(inittraceFILE,"\n%s\n",initversion);
          timevar = time(0);
          fprintf(inittraceFILE, "%s\n", ctime(&timevar) );
       } /* end if */
    } /* end if */

    graphictrace = False;
    if (fopen(GRAPHICTRACEFILE,"r") != NULL)
    {
       if ((graphictraceFILE = fopen(GRAPHICTRACEFILE,"a")) != NULL)
       {
          graphictrace = True;
          fprintf(graphictraceFILE,"\n%s\n",graphicversion);
          timevar = time(0);
          fprintf(graphictraceFILE, "%s\n", ctime(&timevar) );
       } /* end if */
    } /* end if */
 }

