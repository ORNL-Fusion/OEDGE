#ifndef GRSOFT_H_INCLUDED
#define GRSOFT_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

#ifndef  cray
#if  defined  (VMS) || defined (hpux) || defined (aix)
#define GR00DG gr00dg 
#define GR90DG gr90dg
#define GRARRW grarrw
#define GRAXLIN graxlin
#define GRAXLOG graxlog
#define GRAXS graxs
#define GRAXSL graxsl
#define GRBLD 	grbld
#define GRCHN 	grchn
#define GRCHNC 	grchnc
#define GRCHRC 	grchrc
#define GRCLP 	grclp
#define GRCRCL 	grcrcl
#define GRDCUR 	grdcur
#define GRDEL 	grdel
#define GRDN 	grdn
#define GRDRAX 	grdrax
#define GRDRDM 	grdrdm
#define GRDRDU 	grdrdu
#define GRDRHS 	grdrhs
#define GRDRKU 	grdrku
#define GRDRLG 	grdrlg
#define GRDRNE 	grdrne
#define GRDRW 	grdrw
#define GRDRWS 	grdrws
#define GRDSH 	grdsh
#define GREND 	grend
#define GRENDE 	grende
#define GRFILL 	grfill
#define GRFONT 	grfont
#define GRFRBN 	grfrbn
#define GRGFLD 	grgfld
#define GRJMP 	grjmp
#define GRJMPS 	grjmps
#define GRHHFL 	grhhfl
#define GRHHNL 	grhhnl
#define GRHPAN 	grhpan
#define GRLGND 	grlgnd
#define GRLN 	grln
#define GRLNCN	grlncn
#define GRMRKS 	grmrks
#define GRMSKN 	grmskn
#define GRMSKF 	grmskf
#define GRNWPN 	grnwpn
#define GRNXTF 	grnxtf
#define GRPTS 	grpts
#define GRSCAX 	grscax
#define GRSCDL 	grscdl
#define GRSCLC 	grsclc
#define GRSCLP 	grsclp
#define GRSCLV 	grsclv
#define GRSHD 	grshd
#define GRSHOW 	grshow
#define GRSPHR 	grsphr
#define GRSPTS 	grspts
#define GRSTRT 	grstrt
#define GRTXT 	grtxt
#define GRTXTC 	grtxtc
#define GRVAR 	grvar
#define GRVFLD 	grvfld
#define GRWIN 	grwin
#define GRWINC 	grwinc
#define GRWINV 	grwinv
#define GSTXAL 	gstxal
#define KURVEF 	kurvef
#define PSKURF 	pskurf
#define PSKURL 	pskurl
#define SKURF 	skurf
#define SKURL 	skurl
#define GR3ANT 	gr3ant
#define GR3ARR 	gr3arr
#define GR3AXS 	gr3axs
#define GR3BKS 	gr3bks
#define GR3BLK 	gr3blk
#define GR3CEN 	gr3cen
#define GR3DIM 	gr3dim
#define GR3EXT 	gr3ext
#define GR3FPC 	gr3fpc
#define GR3FUN 	gr3fun
#define GR3HHL 	gr3hhl
#define GR3IMP 	gr3imp
#define GR3KIL 	gr3kil
#define GR3MRK 	gr3mrk
#define GR3NET 	gr3net
#define GR3NT1 	gr3nt1
#define GR3OBB 	gr3obb
#define GR3OBE 	gr3obe
#define GR3OBP 	gr3obp
#define GR3ORI 	gr3ori
#define GR3PAN 	gr3pan
#define GR3PLO 	gr3plo
#define GR3ROT 	gr3rot
#define GR3RTX 	gr3rtx
#define GR3STE 	gr3ste
#define GR3TRC 	gr3trc
#define GR3TRF 	gr3trf
#define GR3TRI 	gr3tri
#define GR3TRL 	gr3trl
#define GR3TUB 	gr3tub
#define GR3WIN 	gr3win
#define GR3ZBU 	gr3zbu
	
#else

#define GR00DG gr00dg_
#define GR90DG gr90dg_
#define GRARRW grarrw_
#define GRAXLIN graxlin_
#define GRAXLOG graxlog_
#define GRAXS graxs_
#define GRAXSL graxsl_
#define GRBLD 	grbld_
#define GRCHN 	grchn_
#define GRCHNC 	grchnc_
#define GRCHRC 	grchrc_
#define GRCLP 	grclp_
#define GRCRCL 	grcrcl_
#define GRDCUR 	grdcur_
#define GRDEL 	grdel_
#define GRDN 	grdn_
#define GRDRAX 	grdrax_
#define GRDRDM 	grdrdm_
#define GRDRDU 	grdrdu_
#define GRDRHS 	grdrhs_
#define GRDRKU 	grdrku_
#define GRDRLG 	grdrlg_
#define GRDRNE 	grdrne_
#define GRDRW 	grdrw_
#define GRDRWS 	grdrws_
#define GRDSH 	grdsh_
#define GREND 	grend_
#define GRENDE 	grende_
#define GRFILL 	grfill_
#define GRFONT 	grfont_
#define GRFRBN 	grfrbn_
#define GRGFLD 	grgfld_
#define GRJMP 	grjmp_
#define GRJMPS 	grjmps_
#define GRHHFL 	grhhfl_
#define GRHHNL 	grhhnl_
#define GRHPAN 	grhpan_
#define GRLGND 	grlgnd_
#define GRLN 	grln_
#define GRLNCN	grlncn_
#define GRMRKS 	grmrks_
#define GRMSKN 	grmskn_
#define GRMSKF 	grmskf_
#define GRNWPN 	grnwpn_
#define GRNXTF 	grnxtf_
#define GRPTS 	grpts_
#define GRSCAX 	grscax_
#define GRSCDL 	grscdl_
#define GRSCLC 	grsclc_
#define GRSCLP 	grsclp_
#define GRSCLV 	grsclv_
#define GRSHD 	grshd_
#define GRSHOW 	grshow_
#define GRSPHR 	grsphr_
#define GRSPTS 	grspts_
#define GRSTRT 	grstrt_
#define GRTXT 	grtxt_
#define GRTXTC 	grtxtc_
#define GRVAR 	grvar_
#define GRVFLD 	grvfld_
#define GRWIN 	grwin_
#define GRWINC 	grwinc_
#define GRWINV 	grwinv_
#define GSTXAL 	gstxal_
#define KURVEF 	kurvef_
#define PSKURF 	pskurf_
#define PSKURL 	pskurl_
#define SKURF 	skurf_
#define SKURL 	skurl_
#define GR3ANT 	gr3ant_
#define GR3ARR 	gr3arr_
#define GR3AXS 	gr3axs_
#define GR3BKS 	gr3bks_
#define GR3BLK 	gr3blk_
#define GR3CEN 	gr3cen_
#define GR3DIM 	gr3dim_
#define GR3EXT 	gr3ext_
#define GR3FPC 	gr3fpc_
#define GR3FUN 	gr3fun_
#define GR3HHL 	gr3hhl_
#define GR3IMP 	gr3imp_
#define GR3KIL 	gr3kil_
#define GR3MRK 	gr3mrk_
#define GR3NET 	gr3net_
#define GR3NT1 	gr3nt1_
#define GR3OBB 	gr3obb_
#define GR3OBE 	gr3obe_
#define GR3OBP 	gr3obp_
#define GR3ORI 	gr3ori_
#define GR3PAN 	gr3pan_
#define GR3PLO 	gr3plo_
#define GR3ROT 	gr3rot_
#define GR3RTX 	gr3rtx_
#define GR3STE 	gr3ste_
#define GR3TRC 	gr3trc_
#define GR3TRF 	gr3trf_
#define GR3TRI 	gr3tri_
#define GR3TRL 	gr3trl_
#define GR3TUB 	gr3tub_
#define GR3WIN 	gr3win_
#define GR3ZBU 	gr3zbu_




#endif
#endif /*cray */

 
#include "grsoftsys.h"
 
/*******************************/
/*                             */
/* C-Interface zur GR-Software */
/* grsoft.h                    */
/*                             */
/*******************************/

#define FLOAT_MATRIX void *
#define INT_MATRIX void *
 

extern void grstrt_c(int camera,int ddnumb);
extern void grend_c(void);
extern void grende_c(void);
extern void grnxtf_c(void);
extern void gr3obb_c(void);
extern void grmskn_c(void);
extern void grmskf_c(void);
extern void grshow_c(void);
extern void grsclp_c(float xcm, float ycm, int rahmen);
extern void grsclc_c(float xa,float ya,float xb,float yb);
extern void grsclv_c(float xa,float ya,float xb,float yb);
extern void grdn_c(int idin,float *xcm, float *ycm);
extern void grclp_c(int iclip);
extern void grwin_c(float x1,float y1,float x2,float y2, int idummy);
extern void grwinc_c(float x1,float y1,float x2,float y2, int idummy);
extern void grwinv_c(float x1,float y1,float x2,float y2, int idummy);
extern void grchrc_c(float height,float angle,int idummy);
extern void grdsh_c(float a1,float a2,float a3);
extern void grfont_c(int ifont);
extern void grmrks_c(float height);
extern void grnwpn_c(int ipen);
extern void grspts_c(int int1);
extern void grln_c(float *xx,float *yy,int m);
extern void grlncn_c(float *xx,float *yy,int m);
extern void grchn_c(float *xx,float *yy,int m,int nr);
extern void grchnc_c(float *xx,float *yy,int m,int nr);
extern void grpts_c(float *xx,float *yy,int m,int nr);
extern void grjmp_c(float x,float y);
extern void grjmps_c(float x,float y,int nr);
extern void grdrw_c(float x, float y);
extern void grdrws_c(float x, float y,int nr);
extern void grvar_c(float *xx,float *yy,float *var,int n,int koax);
extern void grcrcl_c(float xm, float ym, float r,float phi1,float phi2);
extern void grarrw_c(float xp, float yp, float xtip,float ytip,                                      float alen,float awid,int icode);
extern void gstxal_c(int ialh,int ialv);
extern void grfill_c(int n,float *xx,float *yy,int istyle,int itype);
extern void grshd_c(float *xx1,float *yy1,float *xx2,float *yy2,float abst,                             float winkel,int n1,int n2);
extern void graxsl_c(int kindx,int kindy,int kinddec);
extern void graxlin_c(float x1,float y1,float x2,float y2,float s1,                                   float s2,int links,int kindax);
extern void graxlog_c(float x1,float y1,float x2,float y2,float s1,                                   float s2,int links,int kindax);
extern void grdcur_c(int np,float *xxp, float *yyp);
extern void grhhfl_c(int nx, int ix, float *vx,int ny,int iy,float *yy,                              int nw,float *vw,int *ivfa,int ndx,FLOAT_MATRIX tab,                            int intact);
extern void grgfld_c(int imax,FLOAT_MATRIX f,int istax, int incx,int nx,                             int istay,int incy,int ny,int kind,float hoehe,                                 float breite,float winkel,float dicke,int *ier);
extern void grvfld_c(int ndim1,FLOAT_MATRIX va,FLOAT_MATRIX vb,                                      INT_MATRIX iifa,int istax,int incx,int nx,int istay,                             int incy,int ny,int kind,float hoehe,float breite,                              float winkel,float dicke,int *ier);
extern void grdrdm_c(float *drdmpa,int nrow,FLOAT_MATRIX tab,float *xx,                                    float *yy);
extern void grdrdu_c(float *drdmpa,int nrow,FLOAT_MATRIX tab,float *xx,                                    float *yy);
extern void grdrku_c(float *drdmpa,int n,float *xyz);
extern void grdrne_c(float *parm,int nrow,float *xyz);
extern void grdrhs_c(float *drdmpa,int npn,FLOAT_MATRIX pnkt,float *xx,                                    float *yy);
extern void grfrbn_c(int ifa,int ika,int iks,int ira,int iri);
extern void grsphr_c(int ns,int nr,FLOAT_MATRIX vs,int *infa,float chi,                                     float psi,float zp,float *hs,int *ipo,int mdm,                                  int *ier);
extern void gr3dim_c(int lar,int *ier);
extern void gr3net_c(float *ar,int *ier, int n1dm,FLOAT_MATRIX xyz,int n1,                                  int i1,int n2,int i2,int ntrnsp,int icol);
extern void gr3nt1_c(float *ar,int *ier, int n1dm,FLOAT_MATRIX x,                                    FLOAT_MATRIX y,FLOAT_MATRIX z,int n1,int i1,int n2,                             int i2,int ntrnsp, int icol);
extern void gr3fun_c(float *ar,int *ier, int nxd,FLOAT_MATRIX fumat,int nx,                                 int ix,float *xx,int ny,int iy,float *yy,int ntrnsp,                              int icol);
extern void gr3hhl_c(float *ar,int *ier,int nx,int ix,float *xx,int ny,                                int iy,float *yy,int nv,float *vv,int nxd,                                        FLOAT_MATRIX fumat,int kind,int ntrnsp,int icol);
extern void gr3imp_c(float *ar,int *ier,int nxd,int nyd,FLOAT_MATRIX f,                                     int nx,int ix,float *xx,int ny,int iy,float *yy,                                  int nz,int iz,float *zz,float c,int ntrnsp,int icol);
extern void gr3tri_c(float *ar,int *ier,FLOAT_MATRIX xyz,int nparts, int *npoint,int *nbound,int ivord,int ntrnsp,int icol);
extern void gr3mrk_c(float *ar,int *ier,FLOAT_MATRIX zent,float *rr,int *kkind,                                int np,int ntrnsp,int icol);
extern void gr3arr_c(float *ar,int *ier,FLOAT_MATRIX xyz,FLOAT_MATRIX vec,int nv,                                 int kind,float rmax,int ntrnsp,int icol);
 extern void gr3blk_c(float *ar,int *ier,FLOAT_MATRIX cenbas,FLOAT_MATRIX dims,                                    int nq,int ntrnsp,int icol);
extern void gr3tub_c(float *ar,int *ier,FLOAT_MATRIX spacur,int ncur,int ncirc,                             float rcirc,int closed,int ntrnsp,int icol);
extern void gr3obp_c(float *ar,int *ier,FLOAT_MATRIX elem,int m34,int nelem,INT_MATRIX lin);
extern void gr3obe_c(float *ar,int *ier,int ntrnsp,int icol);
extern void gr3axs_c(float *ar,int *ier,FLOAT_MATRIX ex33,FLOAT_MATRIX valu,char chaxs[3][20],int chaori,int kindax,int icol);
extern void gr3pan_c(float *ar);
extern void gr3ext_c(float *ar,int *ier,FLOAT_MATRIX ext);
extern void gr3trc_c(float u,float v,float w, float x,float y,float z);
extern void gr3bks_c(int i1,int i2,int i3,int i4,int i5,int i6,int i7,                               int i8);
extern void gr3fpc_c(float wx,float wy,float wz,int k0);
extern void gr3kil_c(int *ier);
extern void gr3rot_c(float *ar,int *ier,char *ax1,float wi1,                                          char *ax2, float wi2,char *ax3,float wi3);
extern void gr3ori_c(float *ar, int *ier);
extern void gr3rtx_c(FLOAT_MATRIX rot,float wi,char *ax);
extern void gr3trf_c(FLOAT_MATRIX rot,float *xx,float *yy,float *zz,int n,                                    int inc);
extern void gr3trl_c(float *rot,int *ier,float *tv);
extern void gr3win_c(FLOAT_MATRIX cor,int *ier);
extern void gr3cen_c(float *ar,int *ier,float za);
extern void kurvef_c(float *xx,float *yy,int ist,int isy);
extern void skurf_c(float *xx,float *yy,int ist,float form,float *zz);
extern void skurl_c(float *xx,float *yy,int ist,float form,float *zz,                                  int ksk);
extern void pskurf_c(float *xx,float *yy,int ist,float form,float *zz);
extern void pskurl_c(float *xx,float *yy,int ist,float form,float *zz,                                  int ksk);
extern void grbld_c(float xcm,float ycm,int isk,int jsk, float xmin,                                float xmax,float ymin, float ymax,int nkurv);
extern void grlgnd_c(char tx[14][48]);



extern void graxs_c(int lopt, char *option, int ltxtx, char *textx,                                 int ltxtxy,char *textxy);
extern void grdrax_c(int lx,char *txtx,int ly, char *txty, int lz,                                   char *txtz,float brtxt);
extern void grdrlg_c(float *drdmpa,char *txtx,char *txty,char *txtz,                                 int iopt,float x,float y, float z, int ix,                                   int iy, int iz,int lx,int ly, int lz);
extern void grhhnl_c(int ndimx,FLOAT_MATRIX tab,int nx, float *xx,int ny,                            float *yy, int l1,float *wert,int *integ,                                       char *option,int istax,int incx, int istay,                                     int incy);
extern void grhpan_c(int n1,float *tab,int nx, float *xx,int ny,                                      float *yy,int l1,float *wert, int *int1,                                        char *option, int istax, int incx,                                              int istay,int incy);
extern void grscax_c(char *cx,char *cy,char *cz1, char *cz2);
extern void grscdl_c(int n,float *xx, float *yy,int idash,char *cha,int nh,                           float *xh, float *yh);
extern void grtxt_c(float x,float y, int len,char *str);
extern void grtxtc_c( int len,char *str);
extern void gr3ant_c(float *ar,int *ier,float *ptxt,char *antxt,                                      float siztxt,float ptrlen, float ptrang, int idireg,                            int ibound,int ifont, int icol);
extern void gr3plo_c(float *ar, int *ier, char *how);
extern void gr3ste_c(float *ar, int *ier, float centr, char *art,                                      int ioculi);
extern void gr3zbu_c(float *ar, char *cc);

#ifdef __cplusplus
}
#endif

#endif /* grsoft_h_included */
