/* grsoftsys.h   */

#ifdef __cplusplus
extern "C" {
#endif

extern void grarrw(float *, float *, float *, float *, float *, float *, int *);
extern void graxlin(float *, float *, float *, float *, float *, float *, int *, int *);
extern void graxlog(float *, float *, float *, float *, float *, float *, int *, int *);
extern void graxsl(int *, int *, int *);
extern void grbld(float *, float *, int *, int *, float *, float *, float *, float *, int *);
extern void grchn(float *, float *, int *, int *);
extern void grchnc(float *, float *, int *, int *);
extern void grchrc(float *, float *, int *);
extern void grclp(int *);
extern void grcrcl(float *, float *, float *, float *, float *);
extern void grdel(void);
extern void grdcur(int *, float *, float *);
extern void grdn(int *, float *, float *);
extern void grdrdm(float *, int *, float *, float *, float *);
extern void grdrku(float *, int *, float *);
extern void grdrdu(float *, int *, float *, float *, float *);
extern void grdrhs(float *, int *, float *, float *, float *);
extern void grdrne(float *, int *, float *);
extern void grdrw(float *, float * );
extern void grdrws(float *, float *, int *);
extern void grdsh(float *, float *, float *);
extern void grend(void);
extern void grende(void);
extern void grfrbn(int *, int *, int *, int *, int *);
extern void grfill(int *, float *, float *, int *, int *);
extern void grfont(int *);
extern void grgfld(int *, float *, int *, int *, int *, int *, int *, int *, int *, float *, float *, float *, float *, int *);
extern void grhhfl(int *, int *, float *, int *, int *, float *,int * , float *, int *, int *, float *, int *);
extern void grjmp(float *, float *);
extern void grjmps(float *, float *, int *);
extern void grlgnd(char *,int *);
extern void grln(float *, float *, int *);
extern void grlncn(float *, float *, int *);
extern void grlnlg(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, int *);
extern void grmrks(float *);
extern void grmskf(void);
extern void grmskn(void);
extern void grnwpn(int *);
extern void grnxtf(void);
extern void grpctr(int *, int *, char);
extern void grpts(float *, float *, int *, int *);
extern void grscla(float *,float *, int *, float *);
extern void grsclc(float *, float *, float *, float *);
extern void grsclp(float *, float *, int *);
extern void grsclv(float *, float *, float *, float *);
extern void grshd(float *, float *, float *, float *, float *, float *, int *,int *);
extern void grshow(void);
extern void grsphr(int *, int *, float *, int *, float *, float *, float *, float * , int *, int *, int *);
extern void grspts(int *);
extern void grstrt(int *, int *);
extern void grvar( float *, float *,float *,int *,int *);
extern void grvfld(int *, float *, float *, int *, int *, int *, int *, int *, int *, int *, int *, float *, float *, float *, float *, int *);
extern void grwin(float *, float *, float *, float *, int *);
extern void grwinc(float *, float *, float *, float *, int *);
extern void grwinv(float *, float *, float *, float *, int *);
extern void gr3arr(float *, int *, float *, float *, int *, int *, float *, int *, int *);
extern void gr3axs(float *, int *, float *, float *, char *, int *, int *, int *);
extern void gr3bks(int *, int *, int *, int *, int *, int *, int *, int *);
extern void gr3blk(float *, int *, float *, float *, int *, int *, int *);
extern void gr3cen(float *, int *, float *);
extern void gr3dim(int *, int *);
extern void gr3ext(float *, int *, float *);
extern void gr3fpc(float *, float *, float *, int *);
extern void gr3kil(int *);
extern void gr3net(float *, int *, int *, float *, int *, int *, int *, int *, int *, int *);
extern void gr3nt1(float *, int *, int *, float *, float *, float *,
       int *, int *, int *, int *, int *, int *);
extern void gr3fun( float *, int *, int *, float *, int *, int *, float *, int *, int *, float *, int *, int *);
extern void gr3hhl( float *,int * , int *, int *, float *, int *, int *, float *, int *, float *, int *, float *, int *, int *, int *);
extern void gr3imp( float *, int *, int *, int *, float *, int *, int *, float *, int *, int *, float *, int *, int *, float *, float *, int *, int *);
extern void gr3mrk(float *, int *, float *, float *, int *, int *, int *, int *);
extern void gr3obb(void);
extern void gr3obe(float *, int *, int *, int *);
extern void gr3obp(float *, int *, float *, int *, int *, int *);
extern void gr3trc(float *, float *, float *, float *, float *, float *);
extern void gr3tri(float *, int *, float *, int *, int *, int *, int *, int *,int *);
extern void gr3tub( float *, int *, float *, int *, int *, float *, int *, int *, int *);
extern void gr3ori(float *, int *);
extern void gr3rot(float *, int *, char *, float *, char *, float *, char *, float *,int ,int , int );
extern void gr3pan(float *);
extern void gr3rtx(float *, float *, char *,int );
extern void gr3trc(float *, float *, float *, float *, float *, float *);
extern void gr3trl(float *, int *, float *);
extern void gr3trf(float *, float *, float *, float *, int *, int *);
extern void gr3win(float *, int *);
extern void gr00dg(void);
extern void gr90dg(void);
extern void gstxal(int *, int *);
extern void kurvef(float *, float *, int *, int *);
extern void pskurf(float *, float *, int *, float *, float *);
extern void pskurl(float *, float *, int *, float *, float *, int *);
extern void skurf(float *, float *, int *, float *, float *);
extern void skurl(float *, float *, int *, float *, float *, int *);

 

extern void graxs(int *, char *, int *, char *, int * ,char *,                                    int , int , int);
extern void grdrax(int * ,char *, int * ,char *, int * ,char *, float *,                           int ,int ,int );
extern void grdrlg(float *, char *, char *, char *, int * ,float *,                         float * ,float *,int * ,int *,int *,int *,int *,int *,                          int ,int ,int );
extern void grhhnl(int, float *, int * ,float *, int * ,float *, int *,                     float *,int *, char *, int * ,int * ,int * ,int * ,int);
extern void grhpan(int * ,float *, int * ,float *, int * ,float *,                                 int *, float *,int *, char *, int *,int *,int *,                                int *,int);
extern void grscax(char *, char *, char *, char *,int, int , int, int);
extern void grscdl(int *, float *, float *,int *, char *, int *,                                   float *, float *, int );
extern void grtxt(float *, float *, int *, char *,int);
extern void grtxtc(int *, char *,int );
extern void gr3ant(float *, int *, float *, char *, float *, float *,                              float *, int *, int *, int *, int *,int);
extern void gr3plo(float *, int *, char *,int );
extern void gr3ste(float *, int *, float *,char *, int *,int );
extern void gr3zbu(float *, char *,int );
 
#ifdef __cplusplus
}
#endif
