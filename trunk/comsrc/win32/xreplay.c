/*
*
*  3 routines provide the Record/replay facility required by 
*  the Expose Event interrupt routine. As Xlib calls are made
*  a copy of the arguments are copied into the Replay Buffer.
*  This can be read at any time to retrieve and Replay all Xlib
*  calls made so far.
*
*  RecordValues(unsigned int *values, unsigned int len) => 
*                0 = success, 1=failure
*  ReplayValues(unsigned int *values) => unsigned int len
*                0 = no data, # = length of data in words (ints)...
*                                        ...pointed to by values
*                                 
*  EraseReplay()
*
*
*  Logging is provided by a log file called "GHreplay.log". If this
*  file exists, in the pwd, then Replay log information is appended
*  to it with a time stamp. Watch out! It can get very large.
*
*  See VERSION for history.
*
*/


#include <stdio.h>

#include "defines.h"

#define LINK_OPCODE 10
#define VALUES_OPCODE 20

/* define the size of the replay buffers and the total number.
   The CHUNKSIZE is in units of words (ints).                    */

/* special behaviour of ardent due to broken malloc/free hence only
   ever allocate one *large* replay buffer  */

#ifdef ardent
#define CHUNKSIZE 2048*2048
#define MAXCHUNKS 1
#endif
#ifdef sg
#define CHUNKSIZE 2048*2048
#define MAXCHUNKS 1
#endif
#ifndef ardent && sg
#define CHUNKSIZE 2048
#define MAXCHUNKS 2048
#endif

/*replay control*/
unsigned int *NextFreeElement, *NextStart, *BaseElement, *ReplayCounter;
unsigned int *BaseElements[MAXCHUNKS], CountItem=0, ReplayItem=0;
int ElementCount = CHUNKSIZE;
int ChunkCount = 0;

/*tracing control*/
extern FILE *replaytraceFILE;
extern Bool replaytrace;
char replayversion[] = "GHOST X11 postprocessor Version 2.3 Replay Buffer";


/* 
*  RecordValues(values, len) adds recorded int's into the Replay
*  buffer. The data is stored in the form:
*
*     ---------------------------------------------------
*     |  type   |   address/length  |  values  |  type  ...
*     ---------------------------------------------------
*
*     type is either LINK_OPCODE or VALUES_OPCODE. LINK_OPCODEs
*     have the address of the next chunk of data as the next word.
*     VALUES_OPCODEs have the number of words of data in the length
*     word.
*
*/



int RecordValues(values, len)
unsigned int *values;
unsigned int len;
{
    if (len <= 0) return;

    /* if the values to store exeeds the amount laft in this buffer
       start a new one. There's two words overhead as well.         */

    if (ElementCount + len + 2 >= CHUNKSIZE)
    {
       if ((NextStart=(unsigned int *)malloc(CHUNKSIZE*sizeof(int))) == NULL)
       {
          fprintf(stderr,"GHOST80-X11: Out of memory\n");
          fflush(stderr);
       }

       if (ChunkCount == 0)
       {
          if (replaytrace)
          {
             fprintf(replaytraceFILE, "CHUNKSIZE= %d\n",CHUNKSIZE);
             fprintf(replaytraceFILE, "MAXCHUNKS= %d\n",MAXCHUNKS);
             fflush(replaytraceFILE);
          }
          BaseElement     = NextStart;
          NextFreeElement = NextStart;
          ReplayCounter   = NextStart;
       }

       if (ChunkCount >= MAXCHUNKS)
       {
          fprintf(stderr, "GHOST80-X11: Internal error: picture ");
          fprintf(stderr, "too complex, increase MAXCHUNKS!\n");
          (void) free(NextStart);
          XCLOSE();
          exit(1);
       }
       BaseElements[ChunkCount] = NextStart;

       ChunkCount++;

      *NextFreeElement++ = LINK_OPCODE;
      *NextFreeElement++ = (unsigned int) NextStart;
       NextFreeElement = NextStart;
       ElementCount = 0;

#ifdef DEBUG
       if (replaytrace)
       {
          fprintf(replaytraceFILE, "ChunkCount= %d  ",ChunkCount);
          fprintf(replaytraceFILE, "ReplayCounter= %x ",ReplayCounter);
          fprintf(replaytraceFILE, "BaseElement= %x ",BaseElement);
          fprintf(replaytraceFILE, "NextFreeElement= %x\n",NextFreeElement);
          fprintf(replaytraceFILE, "CountItem= %d\n",CountItem);
          fflush(replaytraceFILE);
       }
#endif

    } /* end create new chunk clause */

   *NextFreeElement++ = VALUES_OPCODE;
   *NextFreeElement++ = len;
    bcopy(values, NextFreeElement, len * sizeof(int));
    NextFreeElement += len;
    ElementCount += len+2;
    CountItem += len;
#ifdef DEBUG_LOW
    if (replaytrace)
    {
       fprintf(replaytraceFILE, "Replay chunk used= %d  ",ElementCount);
       fprintf(replaytraceFILE, "ReplayCounter= %x ",ReplayCounter);
       fprintf(replaytraceFILE, "BaseElement= %x ",BaseElement);
       fprintf(replaytraceFILE, "NextFreeElement= %x\n",NextFreeElement);
       fprintf(replaytraceFILE, "CountItem= %d\n",CountItem);
       fflush(replaytraceFILE);
    }
#endif

}

int ReplayValues(values)
unsigned int *values;
{
    int len;

    if (ReplayItem >= CountItem)
    {
#ifdef DEBUG
       if (replaytrace)
       {
          fprintf(replaytraceFILE, "Replay RESET. ");
          fprintf(replaytraceFILE, "ReplayCounter= %x ",ReplayCounter);
          fprintf(replaytraceFILE, "BaseElement= %x ",BaseElement);
          fprintf(replaytraceFILE, "NextFreeElement= %x\n",NextFreeElement);
          fprintf(replaytraceFILE, "ReplayItem= %d\n",ReplayItem);
          fprintf(replaytraceFILE, "CountItem= %d\n",CountItem);
          fflush(replaytraceFILE);
       }
#endif
       ReplayCounter = BaseElement;  /* reset pointer in buffer */
       ReplayItem = 0;               /* reset item counter */
       return 0; /* end of Replay buffer */
    }

    if (*ReplayCounter++ == LINK_OPCODE)
    {
       /* deal with a link to a new block */
       ReplayCounter = (unsigned int *) *ReplayCounter;
#ifdef DEBUG
       if (replaytrace)
       {
          fprintf(replaytraceFILE, "Replay new block. ");
          fprintf(replaytraceFILE, "ReplayCounter= %x ",ReplayCounter);
          fprintf(replaytraceFILE, "BaseElement= %x ",BaseElement);
          fprintf(replaytraceFILE, "NextFreeElement= %x\n",NextFreeElement);
          fprintf(replaytraceFILE, "ReplayItem= %d\n",ReplayItem);
          fflush(replaytraceFILE);
       }
#endif
       if (*ReplayCounter++ != VALUES_OPCODE)
       {
          /* the first word in a new chunk MUST be a VALUES_OPCODE */
          fprintf(stderr,
          "GHOST80-X11: Replay Buffer contains garbage at");
          fprintf(stderr,         " beginning of chunk\n");
          fflush(stderr);
          return 0;
       }
    }
    len = *ReplayCounter++; /* get values */
    bcopy(ReplayCounter, values, len*sizeof(int));
    ReplayCounter += len;
    ReplayItem += len;
    return len;
}

int EraseReplay()
{
    int i;

#ifdef DEBUG
    if (replaytrace)
    {
       fprintf(replaytraceFILE, "EraseReplay\n");
       fprintf(replaytraceFILE, "Chunks to free= %d\n", ChunkCount);
       fprintf(replaytraceFILE, "Items in Replay Buffer= %d\n",CountItem);
       fflush(replaytraceFILE);
    }
#endif

    for (i=0; i<ChunkCount; i++) {
#ifdef DEBUG
       if (replaytrace)
       {
          fprintf(replaytraceFILE, "freeing %x...\n", BaseElements[i]);
          fflush(replaytraceFILE);
       }
#endif
       (void) free(BaseElements[i]);
    }

    ChunkCount = 0;
    ElementCount = CHUNKSIZE; /* forces new allocation in RecordValues */
    BaseElement     = NULL;
    NextFreeElement = NULL;
    ReplayCounter   = NULL;
    CountItem       = 0;
}
