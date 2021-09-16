/*
 * Self contained interface routine between GHOST-80 and the Real
 * World Graphics PC4000 frame buffer. The hardware routines called
 * are in the libpc.a library. (i.e. -lpc is required when loading
 * programs)
 * Terry Martin 12/04/89
 */

#define PORT    0x300

struct pc4000 {
   unsigned short status;
   unsigned short data;
 };

static volatile struct pc4000 *regs;

g1dvio_ (icomnd, idata, nchars)
int *icomnd; char *idata; int *nchars;
{
int count = 1; char c; unsigned short data_s; unsigned long data_l;
   switch (*icomnd) {
   case 1:              /* ICOMND = 1 : Initialise pc4000 */

      if (openpc4000() < 0) {
      printf("Unable to open pc4000\n");
      exit (1);
      }

   case 2:              /* ICOMND = 2 : Close pc4000 */
      return;

   case 3:              /* Write data to pc4000 */
      while (count <= *nchars) {
         c = *idata++;
         count++;
         switch (c) {
	    case 1:        /* Write one byte to interface */
               out_byte (*idata++);
               count++;
               break;

            case 2:        /* Writes two bytes to interface */
               data_s = *idata << 8 | *(idata+1);
               idata += 2;
               count += 2;
               out_word (data_s);
               break;

	    case 4:        /* Write four bytes to interface */
               data_l = *idata     << 24 |
                        *(idata+1) << 16 |
                        *(idata+2) <<  8 |
                        *(idata+3);
               idata += 4;
               count += 4;
               out_long (data_l);
               break;

            default:       /* Invalid number of bytes */
               printf("Illegal entry in array (g1dvio)\n");
               exit (1);
         }
      }
      return;
   }
} 

openpc4000()
{
    if ((regs = (struct pc4000 *)mappcio(PORT,
                sizeof(struct pc4000))) == 0)
      return(-1);
    else
      return(0);
}

unsigned char inp(addr)
unsigned char *addr;
{
    return(*addr);
}

int outp(addr, data)
unsigned char *addr;
unsigned char data;
{
    *addr = data; wbflush();
    return(0);
}

int outp16(addr, data)
unsigned short *addr;
unsigned short data;
{
    *addr = data; wbflush();
    return(0);
}

out_byte (x)
unsigned char x ;
{
        while (!(inp(&regs->status) & 0x40))
                continue;
        outp(&regs->status , x);
}

out_word (x)
unsigned short x;
{
    while (!(inp(&regs->status) & 0x40))
      continue;
    outp16(&regs->data, x);
}

out_long (x)
unsigned long x;
{
   unsigned short *y;
   y = (unsigned short *) &x;
   out_word (*y++);
   out_word (*y);
}


