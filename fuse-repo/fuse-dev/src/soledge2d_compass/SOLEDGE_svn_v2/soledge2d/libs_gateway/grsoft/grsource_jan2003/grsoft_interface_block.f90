! Diese Version ist F-vertraeglich

module grsoft_interface_block
  implicit none
  public :: skurf,pskurf,kurvef,grwin,grvfld,grtxt,grtxtc,grsphr,grstrt,     &
            grsclv,grsclc,grln,grjmp,grjmps,grend,grende,grvar,grspts,grshd, &
            grshow,grsclp,grscdl,grscax,grrahm,grpts,grnxtf,grnwpn,grmskn,   &
            grmskf,grmrks,grlncn,grlgnd,grhhnl,grhhfl,grfont,grgfld,grfill,  &
            grdsh,grdrws,grdrw,grdrlg,grdrhs,grdrdm,grdrax,grdn,grcrcl,grclp,&
            grchrc,grnwsg,grchnc,grbld,grchn,graxs,graxsl,graxlin,grarrw,    &
            graxlog,gr3trc,gr90dg,gr3ste,gr3zbu,gr3pan,gr3plo,gr3cen,gr3dim, &
            gr3tri,gr3obb,gr3rtx,gr3trf,gr3arr,gr3blk,gr3imp,gr3hhl,gr3mrk,  &
            gr3ant,gr3bks,gr3tub

  interface

    subroutine skurf (x, y, ist, form, z)                              
      integer,intent(in) :: ist
      real,intent(in)    :: form
      real,intent(in),dimension(ist)   :: x,y
      real,intent(in),dimension(7,ist) :: z
    end subroutine skurf

    subroutine pskurf (x, y, ist, form, z)                             
      integer,intent(in) :: ist
      real,intent(in)    :: form
      real,intent(in),dimension(ist)   :: x,y
      real,intent(in),dimension(7,ist) :: z 
    end subroutine pskurf

    subroutine kurvef (x, y, ist, isy)                                 
      integer,intent(in)             :: ist, isy
      real,intent(in),dimension(ist) :: x,y
    end subroutine kurvef

    subroutine grwin (wx1, wy1, wx2, wy2, idummy)                      
      real,intent(in) :: wx1,wy1,wx2,wy2
      integer,intent(in) :: idummy
    end subroutine grwin

    subroutine grvfld (ndim1, va, vb, ifa, istax, incx, nx, istay,    &
         incy, ny, knd, hoehe, breite, winkel, dicke, ier)
      integer,intent(in) :: ndim1, istax, incx, nx, istay, incy, ny,    &
                                 knd
      integer,intent(inout) :: ier
      real,intent(in) :: hoehe, breite, winkel, dicke                                  
      real,intent(in),dimension(ndim1,*) :: va , vb 
      integer,intent(in),dimension(ndim1, * ) :: ifa
    end subroutine grvfld

    subroutine grtxt (x, y, ltext, text)                               
      real,intent(in) :: x, y                                                          
      integer,intent(in) :: ltext                                                      
      character(len=*),intent(in) :: text                                               
    end subroutine grtxt

    subroutine grtxtc (ltext, text)                                    
      integer,intent(in) :: ltext
      character(len=*),intent(in) :: text
    end subroutine grtxtc

    subroutine grsphr (ns, nr, vs, infa, chi, psi, zp, hs, ipo, mdm,  &
         ier)
      integer,intent(in) :: ns, nr, mdm
      integer,intent(in),dimension(ns) :: infa
      integer,intent(in),dimension(*)  :: ipo
      integer,intent(inout)            :: ier
      real,intent(in)                  :: zp, psi, chi
      real,intent(in),dimension(4, * ) :: vs 
      real,intent(in),dimension( * )   :: hs 
    end subroutine grsphr

    subroutine grstrt (camera, ddnumb)                                 
      integer,intent(in) :: camera,ddnumb
    end subroutine grstrt

    subroutine grsclv (xmin, ymin, xmax, ymax)                         
      real,intent(in) :: xmin, ymin, xmax, ymax                                        
    end subroutine grsclv

    subroutine grsclc (xmin, ymin, xmax, ymax)                         
      real,intent(in) :: xmin, ymin, xmax, ymax                                        
    end subroutine grsclc

    subroutine grln (xx, yy, m)                                        
      integer,intent(in) :: m
      real,intent(in),dimension(m) :: xx, yy
    end subroutine grln

    subroutine grjmp (x, y)                                            
      real,intent(in) :: x, y                                                          
    end subroutine grjmp

    subroutine grjmps (x, y, nr)                                       
      integer,intent(in) :: nr
      real,intent(in)    :: x, y
    end subroutine grjmps

    subroutine grend()
    end subroutine grend

    subroutine grende()
    end subroutine grende

    subroutine grvar (x, y, var, n, klip)                              
      integer,intent(in) :: n, klip                                                    
      real,intent(in),dimension(n) :: x,y,var
    end subroutine grvar

    subroutine grspts (ispots)                                         
      integer,intent(in) :: ispots
    end subroutine grspts

    subroutine grshd (x1, y1, x2, y2, d, th, n1, n2)                   
      real,intent(in),dimension(*) :: x1,y1,x2,y2
      real,intent(in)              :: d, th
      integer,intent(in)           :: n1, n2
    end subroutine grshd

    subroutine grshow()
    end subroutine grshow

    subroutine grsclp (xcm, ycm, ibox)                                 
      integer,intent(in) :: ibox
      real,intent(in) :: xcm,ycm
    end subroutine grsclp

    subroutine grscdl (n, x, y, idash, cha, nh, xh, yh)                
      integer,intent(in) :: n,nh,idash
      real,intent(in),dimension(nh) :: xh, yh
      real,intent(in),dimension(n)  :: x, y
      character(len=*),intent(in) :: cha
    end subroutine grscdl

    subroutine grscax (cx, cy, cz1, cz2)                               
      character(len=*),intent(in) :: cx, cy, cz1, cz2
    end subroutine grscax

    subroutine grrahm()
    end subroutine grrahm

    subroutine grpts (x, y, n, nr)                                     
      integer,intent(in) :: n, nr
      real,intent(in),dimension(*) :: x , y
    end subroutine grpts

    subroutine grnxtf()
    end subroutine grnxtf

    subroutine grnwpn (ipen)                                           
      integer,intent(in) :: ipen
    end subroutine grnwpn

    subroutine grmskn()
    end subroutine grmskn

    subroutine grmskf()
    end subroutine grmskf

    subroutine grmrks (hoehe)                                          
      real,intent(in) :: hoehe
    end subroutine grmrks

    subroutine grlncn (xx, yy, m)                                      
      integer,intent(in) :: m
      real,intent(in),dimension(m) ::  xx, yy
    end subroutine grlncn

    subroutine grlgnd (rt)                                             
      character(len=*),intent(in),dimension(6,14) :: rt
    end subroutine grlgnd

    subroutine grhhnl (ndimx, tab, nx, xx, ny, yy, l1, wert, integ,   &
           option, istax, incx, istay, incy)
      integer,intent(in) ::  incx, incy, istax, istay, l1, ndimx, nx, ny
      integer,intent(in),dimension(*) :: integ
      character(len=*),intent(in) :: option
      real,intent(in),dimension(ndimx, * ) :: tab
      real,intent(in),dimension(*) :: xx,yy,wert
    end subroutine grhhnl

    subroutine grhhfl (nx, ix, x, ny, iy, y, nw, w, jco, n1, f, intact)
      integer,intent(in) :: ix, iy, nw, n1, nx, ny, intact
      integer,intent(in),dimension(nw + 1) :: jco
      real,intent(in),dimension(n1, (ny - 1) * iy + 1) :: f 
      real,intent(in),dimension(nw) :: w
      real,intent(in),dimension( (nx - 1) * ix + 1) ::  x
      real,intent(in),dimension( (ny - 1) * iy + 1) ::  y 
    end subroutine grhhfl

    subroutine grfont (ifont)                                          
      integer,intent(in) :: ifont
    end subroutine grfont

    subroutine grgfld (imax, f, istax, incx, nx, istay, incy, ny,     &
         knd, hoehe, breite, winkel, dicke, ier)
      integer,intent(in) :: incy, istay, ny, knd, istax, imax,        &
                                 nx, incx
      integer,intent(inout) :: ier
      real,intent(in) :: hoehe, breite, winkel, dicke
      real,intent(in),dimension(imax, * ) :: f
    end subroutine grgfld

    subroutine grfill (n, xx, yy, istyle, itype)                       
      integer,intent(in) :: n, istyle, itype                                           
      real,intent(in),dimension(n) :: xx, yy
    end subroutine grfill

    subroutine grdsh (a1, a2, a3)                                      
      real,intent(in) :: a1, a2, a3                                                    
    end subroutine grdsh

    subroutine grdrws (x, y, nr)                                       
      integer,intent(in) :: nr
      real,intent(in) :: x, y                                                          
    end subroutine grdrws

    subroutine grdrw (x, y)                                            
      real,intent(in) :: x, y                                                          
    end subroutine grdrw

    subroutine grdrlg (drdmpa, textx, texty, textz, iopt, xx, yy, zz, &
         ix, iy, iz)
      integer,intent(in)              :: iopt, ix, iy, iz
      character(len=*),intent(in)     :: textx, texty, textz
      real,intent(in)                 :: xx, yy, zz
      real,intent(inout),dimension(*) :: drdmpa 
    end subroutine grdrlg

    subroutine grdrhs (pa, np, pt, xx, yy)                             
      integer,intent(in)               :: np
      real,intent(in),dimension(3, np) :: pt
      real,intent(in),dimension(*)     :: xx,yy
      real,intent(inout),dimension(*)  :: pa
    end subroutine grdrhs

    subroutine grdrdm (pa, nrow, tab, xx, yy)                          
      integer,intent(in) :: nrow
      real,intent(in),dimension(nrow, * ) :: tab
      real,intent(in),dimension(*) :: xx,yy
      real,intent(inout),dimension(*) :: pa 
    end subroutine grdrdm

    subroutine grdrax (lx, txtx, ly, txty, lz, txtz, brtxt)            
      integer,intent(in) :: lx, ly, lz
      character(len=*),intent(in) :: txtx, txty, txtz
      real,intent(in) :: brtxt
    end subroutine grdrax

    subroutine grdn (idin, x, y)                                       
      integer,intent(in) :: idin                                                       
      real,intent(out)   :: x,y
    end subroutine grdn

    subroutine grcrcl (xm, ym, r, ph1, ph2)                            
      real,intent(in) :: xm,ym,r,ph1,ph2
    end subroutine grcrcl

    subroutine grclp (iclip)                                           
      integer,intent(in) :: iclip
    end subroutine grclp

    subroutine grchrc (hoehe, winkel, idummy)                          
      integer,intent(in) :: idummy
      real,intent(in) :: hoehe,winkel
    end subroutine grchrc

    subroutine grnwsg()
    end subroutine grnwsg

    ! --- subroutine grhpan
    ! --- end subroutine grhpan

    subroutine grchnc (xx, yy, m, nr)                                  
      integer,intent(in) :: m,nr
      real,intent(in),dimension(*) :: xx, yy
    end subroutine grchnc

    subroutine grbld (xxcm, yycm, iisk, jjsk, xmi, xma, ymi, yma,     &
         nbil)
      integer,intent(in) :: iisk, jjsk, nbil
      real,intent(in) :: xxcm, yycm,xmi,xma,ymi,yma
    end subroutine grbld

    subroutine grchn (xx, yy, m, nr)                                   
      integer,intent(in) :: m, nr
      real,intent(in),dimension(*) :: xx, yy
    end subroutine grchn

    subroutine graxs (lopt, option, ltxtx, textx, ltxtxy, textxy)      
      integer,intent(in) :: lopt, ltxtx, ltxtxy                                        
      character(len=*),intent(in) :: option,textx,textxy
    end subroutine graxs

    subroutine graxsl (ix0, iy0, i10)                                  
      integer,intent(in) :: ix0,iy0,i10
    end subroutine graxsl

    subroutine graxlin (von, wo, bis, hier, unt, ob, links, achse)     
      real,intent(in) :: von, wo, bis, hier, unt, ob                                   
      logical,intent(in) :: links
      integer,intent(in) :: achse
    end subroutine graxlin

    subroutine grarrw (xanf, yanf, xend, yend, alen, awid, icode)      
      real,intent(in) :: xanf, yanf, xend, yend, alen, awid                            
      integer,intent(in) :: icode
    end subroutine grarrw

    subroutine graxlog (x10, y10, x20, y20, s10, s20, links0, achse)   
      real,intent(in) :: x10, y10, x20, y20, s10, s20
      logical,intent(in) :: links0
      integer,intent(in) :: achse
    end subroutine graxlog

    subroutine gr90dg()
    end subroutine gr90dg

    subroutine gr3trc (u, v, w, x, y, z)
      real,intent(in)  :: u, v, w
      real,intent(out) :: x, y, z
    end subroutine gr3trc

    subroutine gr3ste (ar, ier, centr, art, ioculi)                    
      integer,intent(in)      :: ioculi
      integer,intent(inout)   :: ier
      character(len=*),intent(in) :: art
      real,intent(inout),dimension(*) :: ar
      real,intent(in)         :: centr
    end subroutine gr3ste
                       
    subroutine gr3zbu (a, c)
      real,intent(in),dimension(*) :: a
      character(len=*),intent(in)  :: c
    end subroutine gr3zbu

    subroutine gr3pan (ar)
      real,intent(inout),dimension(*) :: ar
    end subroutine gr3pan

    subroutine gr3plo (a, ier, how)                                    
      integer,intent(inout)   :: ier
      real,intent(in),dimension(*) :: a
      character(len=*),intent(in) :: how
    end subroutine gr3plo

    subroutine gr3cen (a, ier, ab)                                     
      integer,intent(inout) :: ier
      real,intent(in)       :: ab
      real,intent(inout),dimension(*) :: a
    end subroutine gr3cen

    subroutine gr3dim (lwork, ier)                                     
      integer,intent(in)    :: lwork
      integer,intent(inout) :: ier
    end subroutine gr3dim

    subroutine gr3tri (ar, ier, xyz, ng, npo, nb, lik, ifaces, jco)    
      integer,intent(in)    :: ifaces, jco, lik, ng
      integer,intent(inout) :: ier
      real,intent(in),dimension(3, * ) :: xyz
      real,intent(inout),dimension(*)  :: ar
      integer,intent(in),dimension(ng) :: nb, npo
    end subroutine gr3tri

    subroutine gr3obb()
    end subroutine gr3obb

    subroutine gr3rtx (rot, theta, ax)                                 
      real,intent(in)      :: theta
      real,intent(out),dimension(3,3) :: rot
      character(len=*),intent(in) :: ax
    end subroutine gr3rtx

    subroutine gr3trf (rot, x, y, z, np, inc)                          
      integer,intent(in) :: inc, np                                                    
      real,intent(in),dimension(3,3) :: rot
      real,intent(inout),dimension(np) :: x, y, z
    end subroutine gr3trf

    subroutine gr3arr (ar, ier, xyz, riclan, n, knd, rmax, ifaces,   &
         jco)
      real,intent(in)       :: rmax
      real,intent(inout),dimension(*) :: ar
      integer,intent(in)    :: ifaces, jco, knd, n
      integer,intent(inout) :: ier
      real,intent(in),dimension(3,n) :: riclan, xyz
    end subroutine gr3arr

    subroutine gr3blk (ar, ier, zent, bdh, n, ifaces, jco)
      real,intent(inout),dimension(*) :: ar
      integer,intent(in)    :: ifaces, jco, n
      integer,intent(inout) :: ier
      real,intent(in),dimension(3,n) :: bdh, zent
    end subroutine gr3blk

    subroutine gr3imp (ar, ie, n1, n2, f, nx, ix, x, ny, iy, y, nz,   &
         iz, z, cf, iff, jco)
      real,intent(inout),dimension(*) :: ar
      real,intent(in)       :: cf, f, x, y, z
      integer,intent(in)    :: iff, ix, iy, iz, jco, n1, n2, nx, ny, nz
      integer,intent(inout) :: ie
    end subroutine gr3imp

    subroutine gr3hhl (ar, ier, nx, ix, x, ny, iy, y, nw, we, n1, f,  &
         ig, ifa, jco)
      integer,intent(in)    :: ifa, ig, ix, iy, jco, n1, nw, nx, ny
      real,intent(in),dimension(n1,*) :: f
      real,intent(in),dimension(nw)   :: we
      real,intent(in),dimension(n1)   :: x
      real,intent(in),dimension(*)    :: y
      real,intent(inout),dimension(*) :: ar
      integer,intent(inout) :: ier
    end subroutine gr3hhl

    subroutine gr3mrk (ar, ier, zent, r, knd, n, ifaces, jco)         
      integer,intent(in) :: ier, ifaces, jco, n                                        
      real,intent(inout),dimension(*) :: ar
      real,intent(in),dimension(n)    :: r
      real,intent(in),dimension(3,n)  :: zent
      integer,intent(in),dimension(n) :: knd
    end subroutine gr3mrk

    subroutine gr3ant (ar, ier, ptxt, antxt, siztxt, ptrlen, ptrang,  &
         idirec, ibound, ifont, icol)
      real,intent(inout),dimension(*) :: ar
      real,intent(in),dimension(3)    :: ptxt
      real,intent(in)                 :: siztxt, ptrlen, ptrang
      character(len=*),intent(in)     :: antxt
      integer,intent(in)              :: idirec, ibound, ifont, icol
      integer,intent(inout)           :: ier
    end subroutine gr3ant

    subroutine gr3bks (i1, i2, i3, i4, i5, i6, i7, i8)                 
      integer,intent(in) :: i1, i2, i3, i4, i5, i6, i7, i8                             
    end subroutine gr3bks

    subroutine gr3tub (ar, ier, poi, n, m, r, closed, ifaces, jco)     
      integer,intent(in)    :: ifaces, jco, n, m
      integer,intent(inout) :: ier
      real,intent(inout),dimension(*) :: ar
      real,intent(in)       :: r
      real,intent(in),dimension(3,n) :: poi
      logical,intent(in)    :: closed
    end subroutine gr3tub

  end interface

end module grsoft_interface_block
