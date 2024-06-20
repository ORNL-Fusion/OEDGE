/* C Interface zur GR-Software fuer AIX                     */
/* Name: g_grinterface.c                                    */
/*                                                          */
/* Autor: M. Busch   email:   ma.busch@kfa-juelich.de       */
/* Datum: 02.01.95                                          */
/* Update: 3.4.95                                           */

#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef cray
#include <fortran.h>
#endif

#include "grsoft.h"


 void spiegeln( float ( * mat), int dim1,int dim2,int dim3)

{
int ll=0;
int ii,jj,kk;

if ( dim3 > 0 )
{
auto float *vv3 = (float *) malloc(dim1*dim2*dim3 *sizeof(float) );

   for (kk=0;kk<dim3;kk++)
   {
       for (jj=0;jj<dim2;jj++)
         {
         for (ii=0;ii<dim1;ii++)
           {
            vv3[ll]=mat[ dim3*dim2*ii + dim3*jj + kk ];
            ll++;
           }
        }
      }
for (ii=0;ii<ll;ii++)
mat[ii]=vv3[ii];
}

else

{
auto float *vv2 = (float *) malloc(dim1*dim2 *sizeof(float) );

   for (jj=0;jj<dim2;jj++)
      {
      for (ii=0;ii<dim1;ii++)
         {
         vv2[ll]=mat[ dim2*ii+jj ];
         ll++;
        }
      }
for (ii=0;ii<ll;ii++)
mat[ii]=vv2[ii];
}
return;
}

void grstrt_c(int camera,int ddnumb)
{
GRSTRT(&camera,&ddnumb);
return;
}

void grend_c(void)
{
GREND();
return;
}

void grende_c(void)
{
GRENDE();
return;
}

void grnxtf_c(void)
{
GRNXTF();
return;
}

void grmskn_c(void)
{
GRMSKN();
return;
}

void grmskf_c(void)
{
GRMSKF();
return;
}

void gr3obb_c(void)
{
GR3OBB();
return;
}

void grshow_c(void)
{
GRSHOW();
return;
}


void grsclp_c(float xcm, float ycm, int rahmen)
{
GRSCLP(&xcm,&ycm,&rahmen);
return;
}

void grsclc_c(float xa,float ya,float xb,float yb)
{
GRSCLC(&xa,&ya,&xb,&yb);
return;
}

void grsclv_c(float xa,float ya,float xb,float yb)
{
GRSCLV(&xa,&ya,&xb,&yb);
return;
}

void grdn_c(int idin,float *xcm, float *ycm)
{
GRDN(&idin,xcm,ycm);
return;
}

void grclp_c(int iclip)
{
GRCLP(&iclip);
return;
}

void grwin_c(float x1,float y1,float x2,float y2, int idummy)
{
GRWIN(&x1,&y1,&x2,&y2,&idummy);
return;
}

void grwinc_c(float x1,float y1,float x2,float y2, int idummy)
{
GRWINC(&x1,&y1,&x2,&y2,&idummy);
return;
}

void grwinv_c(float x1,float y1,float x2,float y2, int idummy)
{
GRWINV(&x1,&y1,&x2,&y2,&idummy);
return;
}

void grchrc_c(float height,float angle,int idummy)
{
GRCHRC(&height,&angle,&idummy);
return;
}

void grdsh_c(float a1,float a2,float a3)
{
GRDSH(&a1,&a2,&a3);
return;
}

void grfont_c(int ifont)
{
GRFONT(&ifont);
return;
}

void grmrks_c(float height)
{
GRMRKS(&height);
return;
}

void grnwpn_c(int ipen)
{
GRNWPN(&ipen);
return;
}

void grspts_c(int int1)
{
GRSPTS(&int1);
return;
}

void grln_c(float *xx,float *yy,int m)
{
GRLN(xx,yy,&m);
return;
}

void grlncn_c(float *xx,float *yy,int m)
{
GRLNCN(xx,yy,&m);
return;
}

void grchn_c(float *xx,float *yy,int m,int nr)
{
GRCHN(xx,yy,&m,&nr);
return;
}

void grchnc_c(float *xx,float *yy,int m,int nr)
{
GRCHNC(xx,yy,&m,&nr);
return;
}

void grpts_c(float *xx,float *yy,int m,int nr)
{
GRPTS(xx,yy,&m,&nr);
return;
}

void grjmp_c(float x,float y)
{
GRJMP(&x,&y);
return;
}

void grjmps_c(float x,float y,int nr)
{
GRJMPS(&x,&y,&nr);
return;
}

void grdrw_c(float x,float y)
{
GRDRW(&x,&y);
return;
}

void grdrws_c(float x,float y,int nr)
{
GRDRWS(&x,&y,&nr);
return;
}

void grvar_c(float *xx,float *yy,float *var,int n,int koax)
{
GRVAR(xx,yy,var,&n,&koax);
return;
}

void grcrcl_c(float xm, float ym, float r,float phi1,float phi2)
{
GRCRCL(&xm,&ym,&r,&phi1,&phi2);
return;
}

void grarrw_c(float xp, float yp, float xtip,float ytip,float alen,                           float awid,int icode)
{
GRARRW(&xp,&yp,&xtip,&ytip,&alen,&awid,&icode);
return;
}

void gstxal_c(int ialh,int ialv)
{
GSTXAL(&ialh,&ialv);
return;
}

void grfill_c(int n,float *xx,float *yy,int istyle,int itype)
{
GRFILL(&n,xx,yy,&istyle,&itype);
return;
}

void grshd_c(float *xx1,float *yy1,float *xx2,float *yy2,float abst,                             float winkel,int n1,int n2)
{
GRSHD(xx1,yy1,xx2,yy2,&abst,&winkel,&n1,&n2);
return;
}

void graxsl_c(int kindx,int kindy,int kinddec)
{
GRAXSL(&kindx,&kindy,&kinddec);
return;
}

void graxlin_c(float x1,float y1,float x2,float y2,float s1,                                   float s2,int links,int kindax)
{
GRAXLIN(&x1,&y1,&x2,&y2,&s1,&s2,&links,&kindax);
return;
}

void graxlog_c(float x1,float y1,float x2,float y2,float s1,                                   float s2,int links,int kindax)
{
GRAXLOG(&x1,&y1,&x2,&y2,&s1,&s2,&links,&kindax);
return;
}

void grdcur_c(int np,float *xxp, float *yyp)
{
GRDCUR(&np,xxp,yyp);
return;
}

void grhhfl_c(int nx, int ix, float *vx,int ny,int iy,float *yy,                              int nw,float *vw,int *ivfa,int ndx,FLOAT_MATRIX tab,                            int intact)
{
GRHHFL(&nx,&ix,vx,&ny,&iy,yy,&nw,vw,ivfa,&ndx,tab,&intact);
return;
}

void  grgfld_c(int imax,FLOAT_MATRIX f,int istax, int incx,int nx,                                  int istay,int incy,int ny,int kind,float hoehe,                                 float breite,float winkel,float dicke,int *ier)
{
GRGFLD(&imax,f,&istax,&incx,&nx,&istay,&incy,&ny,&kind,&hoehe,                         &breite,&winkel,&dicke,ier);
return;
}

void  grvfld_c(int ndim1,FLOAT_MATRIX va,FLOAT_MATRIX vb,INT_MATRIX iifa,                       int istax,int incx,int nx,int istay,int incy,int ny,                            int kind,float hoehe,float breite,float winkel,                                 float dicke,int *ier)
{
GRVFLD(&ndim1,va,vb,iifa,&istax,&incx,&nx,&istay,&incy,&ny,&kind,                       &hoehe,&breite,&winkel,&dicke,ier);
return;
}

void  grdrdm_c(float *drdmpa,int nrow,FLOAT_MATRIX tab,float *xx,                      float *yy)
{
GRDRDM(drdmpa,&nrow,tab,xx,yy);
return;
}

void  grdrdu_c(float *drdmpa,int nrow,FLOAT_MATRIX tab,float *xx,                              float *yy)
{
GRDRDU(drdmpa,&nrow,tab,xx,yy);
return;
}

void  grdrku_c(float *drdmpa,int n,float *xyz)
{
GRDRKU(drdmpa,&n,xyz);
return;
}

void  grdrne_c(float *parm,int nrow,float *xyz)
{
GRDRNE(parm,&nrow,xyz);
return;
}

void  grdrhs_c(float *drdmpa,int npn,FLOAT_MATRIX pnkt,float *xx,                              float *yy)
{
GRDRHS(drdmpa,&npn,pnkt,xx,yy);
return;
}

void  grfrbn_c(int ifa,int ika,int iks,int ira,int iri)
{
GRFRBN(&ifa,&ika,&iks,&ira,&iri);
return;
}

void  grsphr_c(int ns,int nr,FLOAT_MATRIX vs,int *infa,float chi,                              float psi,float zp,float *hs,int *ipo,int mdm,int *ier)
{
GRSPHR(&ns,&nr,vs,infa,&chi,&psi,&zp,hs,ipo,&mdm,ier);
return;
}

void  gr3dim_c(int lar,int *ier)
{
GR3DIM(&lar,ier);
return;
}

void  gr3net_c(float *ar,int *ier, int n1dm,FLOAT_MATRIX xyz,int n1,                            int i1,int n2,int i2,int ntrnsp,int icol)
{
GR3NET(ar,ier, &n1dm,xyz,&n1,&i1,&n2,&i2,&ntrnsp,&icol);
return;
}

void  gr3nt1_c(float *ar,int *ier, int n1dm,FLOAT_MATRIX x,
              FLOAT_MATRIX y,FLOAT_MATRIX z,int n1,int i1,int n2,int i2,int ntrnsp,int icol)
{
GR3NT1(ar,ier, &n1dm,x,y,z,&n1,&i1,&n2,&i2,&ntrnsp,&icol);
return;
}

void  gr3fun_c(float *ar,int *ier, int nxd,FLOAT_MATRIX fumat,int nx,                           int ix,float *xx,int ny,int iy,float *yy,int ntrnsp,                              int icol)
{
GR3FUN(ar,ier,&nxd,fumat,&nx,&ix,xx,&ny,&iy,yy,&ntrnsp,&icol);
return;
}

void  gr3hhl_c(float *ar,int *ier,int nx,int ix,float *xx,int ny,                                int iy,float *yy,int nv,float *vv,int nxd,                                        FLOAT_MATRIX fumat,int kind,int ntrnsp,int icol)
{
GR3HHL(ar,ier,&nx,&ix,xx,&ny,&iy,yy,&nv,vv,&nxd,fumat,&kind,&ntrnsp,                   &icol);
return;
}

void  gr3imp_c(float *ar,int *ier,int nxd,int nyd,FLOAT_MATRIX f,                               int nx,int ix,float *xx,int ny,int iy,float *yy,                                  int nz,int iz,float *zz,float c,int ntrnsp,int icol)
{
GR3IMP(ar,ier,&nxd,&nyd,f,&nx,&ix,xx,&ny,&iy,yy,&nz,&iz,zz,&c,&ntrnsp,                   &icol);
return;
}

void  gr3tri_c(float *ar,int *ier,FLOAT_MATRIX xyz,int nparts, int *npoint,int *nbound,int ivord,int ntrnsp,int icol)
{
GR3TRI(ar,ier,xyz,&nparts,npoint,nbound,&ivord,&ntrnsp,&icol);
return;
}

void  gr3mrk_c(float *ar,int *ier,FLOAT_MATRIX zent,float *rr,int *kkind,                                int np,int ntrnsp,int icol)
{
GR3MRK(ar,ier,zent,rr,kkind,&np,&ntrnsp,&icol);
return;
}

void  gr3arr_c(float *ar,int *ier,FLOAT_MATRIX xyz,FLOAT_MATRIX vec,                            int nv,int kind,float rmax,int ntrnsp,int icol)
{
GR3ARR(ar,ier,xyz,vec,&nv,&kind,&rmax,&ntrnsp,&icol);
return;
}

void  gr3blk_c(float *ar,int *ier,FLOAT_MATRIX cenbas,FLOAT_MATRIX dims,                        int nq,int ntrnsp,int icol)
{
GR3BLK(ar,ier,cenbas,dims,&nq,&ntrnsp,&icol);
return;
}

void  gr3tub_c(float *ar,int *ier,FLOAT_MATRIX spacur,int ncur,int ncirc,                       float rcirc,int closed,int ntrnsp,int icol)
{
GR3TUB(ar,ier,spacur,&ncur,&ncirc,&rcirc,&closed,&ntrnsp,&icol);
return;
}

void  gr3obp_c(float *ar,int *ier,FLOAT_MATRIX elem,int m34,int nelem,                           INT_MATRIX lin)
{
GR3OBP(ar,ier,elem,&m34,&nelem,lin);
return;
}

void  gr3obe_c(float *ar,int *ier,int ntrnsp,int icol)
{
GR3OBE(ar,ier,&ntrnsp,&icol);
return;
}

void  gr3pan_c(float *ar)
{
GR3PAN(ar);
return;
}

void  gr3ext_c(float *ar,int *ier,FLOAT_MATRIX ext)
{
GR3EXT(ar,ier,ext);
return;
}

void  gr3bks_c(int i1,int i2,int i3,int i4,int i5,int i6,int i7,int i8)
{
GR3BKS(&i1,&i2,&i3,&i4,&i5,&i6,&i7,&i8);
return;
}

void  gr3fpc_c(float wx,float wy,float wz,int k0)
{
GR3FPC(&wx,&wy,&wx,&k0);
return;
}

void  gr3kil_c(int *ier)
{
GR3KIL(ier);
return;
}

void  gr3rot_c(float *ar,int *ier,char *ax1,float wi1,                                          char *ax2, float wi2,char *ax3,float wi3)
{
#ifdef cray
    _fcd text1;
    _fcd text2;
    _fcd text3;
     text1 = _cptofcd(ax1, strlen(ax1) );
     text2 = _cptofcd(ax2, strlen(ax2) );
     text3 = _cptofcd(ax3, strlen(ax3) );
     GR3ROT(ar,ier,text1,&wi1,text2,&wi2,text3,&wi3,strlen(ax1),                                  strlen(ax2),strlen(ax3) );
#else
     GR3ROT(ar,ier,ax1,&wi1,ax2,&wi2,ax3,&wi3,strlen(ax1),                                  strlen(ax2),strlen(ax3) );
#endif
return;
}

void  gr3ori_c(float *ar,int *ier)
{
GR3ORI(ar,ier);
return;
}

void  gr3rtx_c(FLOAT_MATRIX rot,float wi,char *ax)
{
#ifdef cray
    _fcd text;
     text = _cptofcd(ax, strlen(ax) );
     GR3RTX(rot,&wi,text,strlen(ax));
#else
     GR3RTX(rot,&wi,ax,strlen(ax));
#endif
return;
}

void  gr3trf_c(FLOAT_MATRIX rot,float *xx,float *yy,float *zz,int n,int inc)
{
GR3TRF(rot,xx,yy,zz,&n,&inc);
return;
}

void  gr3trl_c(float *rot,int *ier,float *tv)
{
GR3TRL(rot,ier,tv);
return;
}

void  gr3win_c(FLOAT_MATRIX cor,int *ier)
{
GR3WIN(cor,ier);
return;
}

void  gr3cen_c(float *ar,int *ier,float za)
{
GR3CEN(ar,ier,&za);
return;
}

void  kurvef_c(float *xx,float *yy,int ist,int isy)
{
KURVEF(xx,yy,&ist,&isy);
return;
}

void  skurf_c(float *xx,float *yy,int ist,float form,float *zz)
{
SKURF(xx,yy,&ist,&form,zz);
return;
}

void  skurl_c(float *xx,float *yy,int ist,float form,float *zz,int ksk)
{
SKURL(xx,yy,&ist,&form,zz,&ksk);
return;
}

void  pskurf_c(float *xx,float *yy,int ist,float form,float *zz)
{
PSKURF(xx,yy,&ist,&form,zz);
return;
}

void  pskurl_c(float *xx,float *yy,int ist,float form,float *zz,int ksk)
{
PSKURL(xx,yy,&ist,&form,zz,&ksk);
return;
}

void  grbld_c(float xcm,float ycm,int isk,int jsk, float xmin,                                float xmax,float ymin, float ymax,int nkurv)
{
GRBLD(&xcm,&ycm,&isk,&jsk,&xmin,&xmax,&ymin,&ymax,&nkurv);
return;
}

void  grlgnd_c(char tx[14][48])
{
GRLGND((char *)tx,strlen((char *)tx));
return;
}

void  gr3axs_c(float *ar,int *ier,FLOAT_MATRIX ex33,FLOAT_MATRIX valu,                          char chaxs[3][20],int chaori,int kindax,int icol)
{
GR3AXS(ar,ier,ex33,valu,chaxs,&chaori,&kindax,&icol);
return;
}


void graxs_c(int lopt, char *option, int ltxtx, char *textx,                                 int ltxtxy,char *textxy)
{
#ifdef cray
     _fcd text1;
     _fcd text2;
     _fcd text3;
     text1 = _cptofcd(option, strlen(option) );
     text2 = _cptofcd(textx, strlen(textx) );
     text3 = _cptofcd(textxy, strlen(textxy) );
     GRAXS(&lopt,text1,&ltxtx,text2,&ltxtxy,text3,                                      (lopt == -1) ? strlen(option) : lopt,                                           (ltxtx == -1) ? strlen(textx) : ltxtx,                                          (ltxtxy == -1) ? strlen(textxy) : ltxtxy);
#else
     GRAXS(&lopt,option,&ltxtx,textx,&ltxtxy,textxy,                                      (lopt == -1) ? strlen(option) : lopt,                                           (ltxtx == -1) ? strlen(textx) : ltxtx,                                          (ltxtxy == -1) ? strlen(textxy) : ltxtxy);
#endif
return;
}

void grdrax_c(int lx,char *txtx,int ly, char *txty, int lz,char *txtz,                        float brtxt)
{
#ifdef cray
    _fcd text1;
    _fcd text2;
    _fcd text3;
    text1 = _cptofcd(txtx, strlen(txtx) );
    text2 = _cptofcd(txty, strlen(txty) );
    text3 = _cptofcd(txtz, strlen(txtz) );

   GRDRAX(&lx,text1,&ly,text2,&lz,text3,&brtxt,                                               (lx == -1) ? strlen(txtx) : lx,                                                 (ly == -1) ? strlen(txty) : ly,                                                 (lz == -1) ? strlen(txtz) : lz);
#else
   GRDRAX(&lx,txtx,&ly,txty,&lz,txtz,&brtxt,                                               (lx == -1) ? strlen(txtx) : lx,                                                 (ly == -1) ? strlen(txty) : ly,                                                 (lz == -1) ? strlen(txtz) : lz);
#endif
      
return;
}

void grdrlg_c(float *drdmpa,char *txtx,char *txty,char *txtz,int iopt,                        float x,float y, float z, int ix, int iy, int iz,                            int lx,int ly, int lz)
{
#ifdef cray
    _fcd text1;
    _fcd text2;
    _fcd text3;
    text1 = _cptofcd(txtx, strlen(txtx) );
    text2 = _cptofcd(txty, strlen(txty) );
    text3 = _cptofcd(txtz, strlen(txtz) );
    GRDRLG(drdmpa,text1,text2,text3,&iopt,&x,&y,&z,&ix,&iy,&iz,&lx,&ly,&lz,                (lx == -1) ? strlen(txtx) : lx,                                                 (ly == -1) ? strlen(txty) : ly,                                                 (lz == -1) ? strlen(txtz) : lz);
#else
    GRDRLG(drdmpa,txtx,txty,txtz,&iopt,&x,&y,&z,&ix,&iy,&iz,&lx,&ly,&lz,                (lx == -1) ? strlen(txtx) : lx,                                                 (ly == -1) ? strlen(txty) : ly,                                                 (lz == -1) ? strlen(txtz) : lz);
#endif
return;
}

void grhhnl_c(int ndimx,FLOAT_MATRIX tab,int nx, float *xx,int ny,                            float *yy,int l1,float *wert,int *integ,char *option,                           int istax,int incx, int istay, int incy)
{
#ifdef cray
    _fcd text;
    text = _cptofcd(option, strlen(option) );
    GRHHNL(&ndimx,tab , &nx, xx, &ny, yy,&l1,wert, integ,text,&istax,                    &incx,&istay,&incy, strlen(option) );
#else
    GRHHNL(&ndimx,tab , &nx, xx, &ny, yy,&l1,wert, integ,option,&istax,                    &incx,&istay,&incy, strlen(option) );
#endif
return;
}

void grhpan_c(int n1,float *tab,int nx,float *xx,int ny,float *yy,int l1,                    float *wert, int *int1, char *option, int istax, int incx,                      int istay,int incy)
{
#ifdef cray
    _fcd text;
    text = _cptofcd(option, strlen(option) );
    GRHPAN(&n1,tab,&nx,xx,&ny,yy,&l1,wert,&int1,text,&istax,&incx,                       &istay,&incy,strlen(option));
#else
    GRHPAN(&n1,tab,&nx,xx,&ny,yy,&l1,wert,&int1,option,&istax,&incx,                       &istay,&incy,strlen(option));
#endif
return;
}

void grscax_c(char *cx,char *cy,char *cz1, char *cz2)
{ 
#ifdef cray
    _fcd text1;
    _fcd text2;
    _fcd text3;
    _fcd text4;
    text1 = _cptofcd(cx, strlen(cx) );
    text2 = _cptofcd(cy, strlen(cy) );
    text3 = _cptofcd(cz1, strlen(cz1) );
    text4 = _cptofcd(cz2, strlen(cz2) );
    GRSCAX(text1,text2,text3,text4,strlen(cx),strlen(cy),strlen(cz1),strlen(cz2));
#else
    GRSCAX(cx,cy,cz1,cz2,strlen(cx),strlen(cy),strlen(cz1),strlen(cz2));
#endif
return;
}

void grscdl_c(int n, float *xx, float *yy,int idash,char *cha,int nh,                           float *xh, float *yh)
{
#ifdef cray
    _fcd text;
    text = _cptofcd(cha, strlen(cha) );
    GRSCDL(&n,xx,yy,&idash,text,&nh,xh,yh,strlen(cha));
#else
    GRSCDL(&n,xx,yy,&idash,cha,&nh,xh,yh,strlen(cha));
#endif
return;
}


void grtxt_c(float x,float y, int len,char *str)
{
#ifdef cray
    _fcd text;
    text = _cptofcd(str, strlen(str ));
    GRTXT(&x,&y, &len ,text, (len == -1) ? strlen(str) : len);
#else
    GRTXT(&x,&y, &len ,str, (len == -1) ? strlen(str) : len);
#endif
return;
}

void grtxtc_c( int len,char *str)
{
#ifdef cray
    _fcd text;
    text = _cptofcd(str, strlen(str ));
    GRTXTC(&len ,text, (len == -1) ? strlen(str) : len);
#else
    GRTXTC(&len ,str, (len == -1) ? strlen(str) : len);
#endif
return;
}

void gr3ant_c(float *ar,int *ier,float *ptxt,char *antxt,float siztxt,                         float ptrlen, float ptrang, int idireg, int ibound,                             int ifont, int icol)
{
#ifdef cray
    _fcd text;
    text = _cptofcd(antxt, strlen(antxt ));
    GR3ANT(ar,ier,ptxt, text, &siztxt,&ptrlen,&ptrang,&idireg,&ibound,                   &ifont,&icol,strlen(antxt));
#else
   GR3ANT(ar,ier,ptxt, antxt, &siztxt,&ptrlen,&ptrang,&idireg,&ibound,                   &ifont,&icol,strlen(antxt));
#endif
return;
}

void gr3plo_c(float *ar, int *ier, char *how)
{
#ifdef cray
    _fcd text;
    text = _cptofcd(how, strlen(how ));
    GR3PLO( ar, ier, text,strlen(how));
#else
    GR3PLO( ar, ier, how,strlen(how));
#endif
return;
}

void gr3ste_c(float *ar, int *ier, float centr, char *art, int ioculi)
{
#ifdef cray
    _fcd text;
    text = _cptofcd(art, strlen(art ));
    GR3STE( ar, ier, &centr,text,&ioculi,strlen(art));
#else
    GR3STE( ar, ier, &centr,art,&ioculi,strlen(art));
#endif
return;
}

void gr3zbu_c(float *ar, char *cc)
{
#ifdef cray
    _fcd text;
    text = _cptofcd(cc, strlen(cc) );
    GR3ZBU( ar ,text,strlen(cc));
#else
    GR3ZBU( ar ,cc,strlen(cc));
#endif
return;
}
