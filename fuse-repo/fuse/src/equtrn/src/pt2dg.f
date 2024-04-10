      program pt2dg

c  version : 25.11.2006 00:33
c=======================================================================
c*** Conversion of equilibrium file produced with PROTEUS into DG format
c*** (courtesy Prof. O. de Barbieri)
c=======================================================================
c                                                                       
#include "pt2dg.inc"                                              
c                                                                       
c-----------------------------------------------------------------------
c                                                                       
      dimension lko(6)
      logical ex,exi
      dimension nrnz(2),rmnmx(2),zmnmx(2),btrt(2)
      character*256 in_mesh,in_psi,out_file,hlp_txt*8
      namelist/input/nrnz,rmnmx,zmnmx,btrt,in_mesh,in_psi,out_file
      data nrnz, rmnmx, zmnmx, btrt, in_mesh, in_psi, out_file  /
     /	    2*0,  2*0.,  2*0., 2*0.,  ' '   , 	' ' , 	' '   	/
c
c                        -----------------------------------------------
c      
c     Predefined function: double of the area of a triangle
c     ------------------- 
      artri(x1,x2,x3,y1,y2,y3) = abs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))
c                                                                 
c-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c                                                                       
cank{
c*** Input the parameters
      ex=.true.
      read(5,input,err=10,end=10)
      ex=.false.
 10   if(ex) then !{
      	write(0,*) 'pt2dg:  error in the parameter file format.'
      	write(0,*) '  	    It must conform to the NAMELIST input.'
	stop
      end if !}
c      read(nparm,*) mpr,mpz,rmin,rmax,zmin,zmax,btf,rtf
c      print *,mpr,mpz,rmin,rmax,zmin,zmax,btf,rtf
c*** check the parameters
      write(0,*) 'Checking the data from the parameter file'
      ex=.false.
      if(nrnz(1).le.0 .or. nrnz(2).le.0) then !{
        write(0,*) 'nrnz must be positive : ',nrnz
        ex=.true.
      end if !}
      if(rmnmx(1).ge.rmnmx(2)) then !{
        write(0,*) 'rmnmx must be accending: ',rmnmx
        ex=.true.
      end if !}
      if(rmnmx(1).le.0.) then !{
        write(0,*) 'rmnmx must be positive: ',rmnmx
        ex=.true.
      end if !}
      if(zmnmx(1).ge.zmnmx(2)) then !{
        write(0,*) 'zmnmx must be accending: ',zmnmx
        ex=.true.
      end if !}
      if(btrt(1).le.0. .or. btrt(2).le.0.) then !{
        write(0,*) 'btrt must be positive : ',btrt
        ex=.true.
      end if !}
      if(in_mesh.eq.' ') then !{
        write(0,*) 'in_mesh file must be specified'
        ex=.true.
      end if !}
      if(in_psi.eq.' ') then !{
        write(0,*) 'in_psi file must be specified'
        ex=.true.
      end if !}
      if(out_file.eq.' ') then !{
        write(0,*) 'out_file file must be specified'
        ex=.true.
      end if !}
      if(ex) then !{
      	write(0,*) 'Please correct the parameter file!'
	stop
      end if  !}
c*** Assign the parameters and check further
      mpr=nrnz(1)
      mpz=nrnz(2)
      rmin=rmnmx(1)
      rmax=rmnmx(2)
      zmin=zmnmx(1)
      zmax=zmnmx(2)
      btf=btrt(1)
      rtf=btrt(2)
      if(mpr.gt.npr) then !{
        ex=.true.
        write(0,*) 'mpr > npr : ',mpr,npr
      end if !}
      if(mpz.gt.npz) then !{
        ex=.true.
        write(0,*) 'mpz > npz : ',mpz,npz
      end if !}
      inquire(file=in_mesh, exist=exi)
      if(.not.exi) then !{
        ex=.true.
        write(0,*) 'in_mesh file not found:',in_mesh
      end if !}
      inquire(file=in_psi, exist=exi)
      if(.not.exi) then !{
        ex=.true.
        write(0,*) 'in_psi file not found:',in_psi
      end if !}
      if(ex) then !{
      	write(0,*) 'Please correct the parameter file!'
	stop
      end if  !}
      write(0,*) '  - look OK'
c*** Open the files
      hlp_txt='in_mesh'
      open(nread,file=in_mesh,err=15)
      hlp_txt='in_psi'
      open(ndisk,file=in_psi,err=15)
      hlp_txt='in_psi'
      open(nwreqdg,file=out_file,err=15)
      ex=.true.
 15   if(.not.ex) then !{
      	write(0,*) 'Failed opening ',hlp_txt
      end if  !}
cank}
c      
c      
c                        The following is the flag for B.pol calculation
c                        If it is less than zero ----> no B.pol
c                        -----------------------------------------------
            kbpol      = 1
c
c
c                        The following is the flag for standard output
c                        (without B.pol even if kbpol=1)
c                        -----------------------------------------------
            kwrite     =-1
c           ++++++
c
                         if(kwrite.gt.0)                  then
                         write(nwrite,2200) npr,npz,mpr,mpz
                                                          endif
c
c                                                                       
c                        +---------------------------+
c--------------------->  | Read mesh from unit NREAD |  <---------------
c                        +---------------------------+
c                                                                       
      read(nread,1003)   (cmesht(lm:lm),lm = 1,80)                      
      read(nread,1000)   mnodes,melems,mhmesh                           
                         if(mnodes.gt.nnodes)               then       
                         write(nwrite,1060) mnodes,nnodes               
                         lkill = 1                                      
                                                            endif       
                         if(melems.gt.nelems)               then       
                         write(nwrite,1070) melems,nelems               
                         lkill = 1                                      
                                                            endif       
      do 100      j100 = 1,mnodes                                       
      read(nread,1001)   lnode, lbound, zr, zz                          
                         rznods(1,lnode)  = zr                          
                         rznods(2,lnode)  = zz                          
                         kbound(lnode)    = lbound                      
  100 continue                                                          
      do 120      j120 = 1,melems                                       
      read(nread,1002)   modele, lelem, keyreg(lelem)                   
     +,                                (lko(lke),lke=1,modele)        
      do 110      j110 = 1,6                                            
                         keynod(j110,lelem) = lko(j110)               
  110 continue                                                          
  120 continue                                                          
c                                                                       
c                                                                       
c                                                                       
c                        +--------------------------+
c--------------------->  | Read PSI from unit NDISK |  <----------------
c                        +--------------------------+
c                                                                       
                         rewind ndisk                                   
                         read(ndisk,*)    ( psi(j1), j1=1,mnodes )      
c                                                                       
c-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c                                                                       
c                                                                       
c                        +--------------------+
c--------------------->  | Quadrilateral mesh |  <----------------------
c                        +--------------------+
c
            lcount     = 0  
            zdr        = (rmax-rmin)/float(mpr-1)
            zdz        = (zmax-zmin)/float(mpz-1)
      do 500       j50 = 1,mpr
            xrp(j50)   = rmin + zdr*float(j50-1)
      do 490       j49 = 1,mpz
            xzp(j49)   = zmin + zdz*float(j49-1)
      do 290       j29 = 1,melems 
      do 200       j20 = 1,6
            lko(j20)   = keynod(j20,j29)
            l20        = lko(j20)
            pnodr(j20) = rznods(1,l20)
            pnodz(j20) = rznods(2,l20)
            psiele(j20)= psi(l20)
  200 continue
            lnwt       = -1
            lreg       = keyreg(j29)
                         if(lreg.le.-100)   then
            lnwt       = 1
                         go to 25
                                            endif
      do 230       j23 = 1,melems
                         if(j23.eq.j29) go to 23
            lr23       = keyreg(j23)
                         if(lr23.le.-100.and.
     +                      lreg.eq.0)         then
      do 220       j22 = 4,6
            lk22       = keynod(j22,j23)
      do 210       j21 = 4,6
            lk21       = keynod(j21,j29)
                         if(lk21.eq.lk22) then
            lnwt       = 1
                         go to 25
                                          endif
  210 continue
  220 continue
                                               endif
   23                    continue
  230 continue
                         go to 27
   25                    continue
                         if(lnwt.gt.0)      then
            kel        = j29
            yrp        = xrp(j50)
            yzp        = xzp(j49)
                         CALL TRINWT
c                        -----=====
                         if(kerror.lt.0) then
      psi4(j50,j49)    = ypsi
            lcount     = lcount + 1
c
                         if(kwrite.gt.0)                  then
                         write(nwrite,2300) j50,j49,j29
     +,                        xrp(j50),xzp(j49),psi4(j50,j49)
                         write(nwrite,2400)
     +                         pnodr(1),pnodr(2),pnodr(3)
     +,                        pnodz(1),pnodz(2),pnodz(3)
     +,                        psiele(1),psiele(2),psiele(3)
                         write(nwrite,2600) csi,eta
                                                          endif
c
                         go to 49
                                         endif
                         go to 29
                                            endif
   27                    continue
            kel        = j29
            yrp        = xrp(j50)
            yzp        = xzp(j49)
            zarea      = artri(pnodr(1),pnodr(2),pnodr(3),
     +                         pnodz(1),pnodz(2),pnodz(3))
            zar12      = artri(pnodr(1),pnodr(2),xrp(j50),
     +                         pnodz(1),pnodz(2),xzp(j49))
            zar13      = artri(pnodr(1),xrp(j50),pnodr(3),
     +                         pnodz(1),xzp(j49),pnodz(3))
            zar23      = artri(xrp(j50),pnodr(2),pnodr(3),
     +                         xzp(j49),pnodz(2),pnodz(3))
            zarsum     = zar12 + zar13 + zar23
            zardif     = abs(zarea - zarsum)
                         if(zardif.lt.1.e-03)                     then
            zr0        = pnodr(3)                                       
            zz0        = pnodz(3)                                       
            za11       = pnodr(1) - pnodr(3)                            
            za12       = pnodr(2) - pnodr(3)                            
            za21       = pnodz(1) - pnodz(3)                            
            za22       = pnodz(2) - pnodz(3)                            
            zdet       = za11*za22 - za12*za21                          
            zb11       = za22/zdet                                      
            zb12       = -za12/zdet                                     
            zb21       = -za21/zdet                                     
            zb22       = za11/zdet                                      
            zx0        = (zz0*za12-zr0*za22)/zdet                       
            ze0        = (zr0*za21-zz0*za11)/zdet                       
            z00e       = psiele(3)                                      
            zx1e       = 4.e00*psiele(5)-3.e00*psiele(3)-psiele(1)      
            zy1e       = 4.e00*psiele(4)-3.e00*psiele(3)-psiele(2)      
            zx2e       = 2.e00*(psiele(1)+psiele(3)-2.e00*psiele(5))    
            zxye       = 4.e00*(psiele(3)+psiele(6)-                    
     +                          psiele(4)-psiele(5))                    
            zy2e       = 2.e00*(psiele(2)+psiele(3)-2.e00*psiele(4))    
            zp00       = z00e + zx1e*zx0 + zy1e*ze0 + zx2e*zx0*zx0 + 
     +                   zxye*zx0*ze0 + zy2e*ze0*ze0                   
            zpr1       = zx1e*zb11 + zy1e*zb21 +
     +                   zxye*(zb11*ze0+zb21*zx0) +
     +                   2.e00*(zx2e*zb11*zx0+zy2e*zb21*ze0)            
            zpz1       = zx1e*zb12 + zy1e*zb22 + 
     +                   zxye*(zb12*ze0+zb22*zx0) +
     +                   2.e00*(zx2e*zb12*zx0+zy2e*zb22*ze0)           
            zpr2       = zx2e*zb11*zb11 +zxye*zb11*zb21 +zy2e*zb21*zb21 
            zprz       = 2.e00*(zx2e*zb11*zb12+zy2e*zb21*zb22) +      
     +                   zxye*(zb11*zb22+zb12*zb21)                    
            zpz2       = zx2e*zb12*zb12 +zxye*zb12*zb22 +zy2e*zb22*zb22
      psi4(j50,j49)    = zp00 + zpr1*xrp(j50) + zpz1*xzp(j49) +
     +                   zpr2*xrp(j50)*xrp(j50)               + 
     +                   zprz*xrp(j50)*xzp(j49)               + 
     +                   zpz2*xzp(j49)*xzp(j49) 
            lcount     = lcount + 1
c
                         if(kwrite.gt.0)                  then
                         write(nwrite,2300) j50,j49,j29
     +,                        xrp(j50),xzp(j49),psi4(j50,j49)
                         write(nwrite,2400) 
     +                         pnodr(1),pnodr(2),pnodr(3)
     +,                        pnodz(1),pnodz(2),pnodz(3)
     +,                        psiele(1),psiele(2),psiele(3)
                         write(nwrite,2500) zarea,zarsum,zardif
                                                          endif
c
            csi        = zx0 + zb11*xrp(j50) + zb12*xzp(j49)
            eta        = ze0 + zb21*xrp(j50) + zb22*xzp(j49)
                                                 go to 49
                                                                  endif
   29                    continue
  290 continue
   49                    continue
c
                         if(kbpol.gt.0)  then
            yrp        = xrp(j50)
            yzp        = xzp(j49)
                         CALL BPOLRZ
c                        -----======
      bpr(j50,j49)     = ybpolr    
      bpz(j50,j49)     = ybpolz   
c--b
cK-cK-cK-cK-cK-cK-cK-cK- if(j49.eq.1) write(nwrite,5000)
cK-cK-cK-cK-cK-cK-cK-cK- write(nwrite,5001) j50,j49,ybpolr,ybpolz 
cK-  +,cK-cK-cK-cK-cK-cK-cK-cK-cK-cK-cK-cK-cK-cK-cK-yrp,cK-yzp
cK-cK-cK-cK-cK-cK-cK-cK- write(nwrite,5002) lnnk,csi,eta
c--e
                                         endif
c
  490 continue
  500 continue
c
c
c
                         if(kwrite.gt.0)                  then
                         write(nwrite,3000) (xrp(jr),jr=1,6)
      do 600       jz = 1,mpz
                         write(nwrite,3100) xzp(jz)
     +,                  psi4(1,jz),psi4(2,jz),psi4(3,jz)
     +,                  psi4(4,jz),psi4(5,jz),psi4(6,jz)
  600 continue
                         write(6,2000) mpr,mpz,lcount
                                                          endif
c
c                                
            lwreqdg    = nwreqdg      
c
                         CALL WREQDG(lwreqdg,npr,npz,kret,mpr,mpz,psib
     +,                              btf,rtf,xrp,xzp,psi4)
c
                         write(nwrite,2100) kret
c                                                                       
c
                         if(kbpol.gt.0)  then
                         CALL WBPOLRZ
c                        -----=======
                                         endif
c
c                                                                       
c-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
c                                                                       
c                        +--------------------------------------------+ 
c                        | FORMATs to read the mesh data (unit NREAD) | 
c                        +--------------------------------------------+ 
c                                                                       
 1000 format(3i5)                                                       
 1001 format(2i5,1p,2e20.8)                                             
 1002 format(9i5)                                                       
 1003 format(80a1)                                                      
 1060 format(10x,'Stop in AAMAIN :',/       
     +,      10x,'mnodes = ',i5,' ###   nnodes = ',i5,/                 
     +,      10x,'Please set nnodes .gt. mnodes ')                      
 1070 format(10x,'Stop in AAMAIN :',/  
     +,      10x,'melems = ',i5,' ###   nelems = ',i5,/                 
     +,      10x,'Please set nelems .gt. melems ')                      
 1100 format(10x,'Stop in AAMAIN :',/
     +,      10x,'label of nodes should be 1',/                         
     +,      10x,'label = ',i2)                                         
c
 2000 format(1h1,5('   ',/),20x,'Message from AAMAIN:'
     +,3('   ',/),t10,'mpr =',i4,t35,'mpz =',i4,t50,'lcount =',i4)
 2100 format(4('   ',/),t15,'After call of WREQDG: kret =',i3)
 2200 format(1h1,5('   ',/),t20,'npr =',i4,t40,'npz =',i4
     +,/,                   t20,'mpr =',i4,t40,'mpz =',i4)
 2300 format(4('   ',/),t13,'j50 =',i3,t43,'j49 =',i3
     +,          t91,'inside element No. ',i5
     +,/,        t10,'xrp(*) =',1p,e12.4,t40,'xzp(*) =',e12.4
     +,          t67,'psi4(*,*) =',e12.4)
 2400 format('   ',/
     +,          t10,'  r(1) =',1p,e12.4,t40,'  r(2) =',e12.4
     +,          t70,'  r(3) =',e12.4,/
     +,          t10,'  z(1) =',   e12.4,t40,'  z(2) =',e12.4
     +,          t70,'  z(3) =',e12.4,/
     +,          t10,'psi(1) =',   e12.4,t40,'psi(2) =',e12.4
     +,          t70,'psi(3) =',e12.4)
 2500 format(t10,' zarea =',1p,e12.4,t40,'zarsum =',e12.4
     +,          t70,'zardif =',e12.4)
 2600 format(t10,'Newton'
     +,          t40,'   csi =',1p,e12.4
     +,          t70,'   eta =',e12.4)
 3000 format(1h1,5('   ',/),t18,1p,e16.4,5e16.4)
 3100 format(               t2 ,1p,e16.4,6e16.4)
c--b
 5000 format('     ',/,'    ')
 5001 format('  ',/
     +,      5x,'j50 =',i3,5x,'j49 =',i3
     +,     10x,'bpolr =',1p,e12.4,10x,'bpolz =',e12.4
     +,/,   t39,'yrp =',e12.4,12x,'yzp =',e12.4)
 5002 format(4x,'lnnk =',i3,t39,'csi =',e12.4,12x,'eta =',e12.4)
c--e
c                                                                       
c-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
c                                                                       
 9999 format(1H1,'          ',/)                                        
c                                                                       
c-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
c                                                                       
c                                                                       
      END                                                             
      SUBROUTINE TRINWT                          !  RECTANGLE.AAMAIN
c     =================  
c                                                                       
c-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c                                                                       
c     Decides whether the point yr,yz is inside the curved side
c     triangle pnodr(*),pnodz(*).
c
c     Input quantities:  yr,yz,pnodr(*),pnodz(*).
c     Output quantities: kerror, csi,eta.
c
c     If the point is inside:  kerror=-1, otherwise kerror=+1.
c                                                                       
c-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c                                                                       
#include "pt2dg.inc"                                              
c                                                                
c-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c
c
            ztol       = 1.e-05
            kerror     = 1
            zr3        = pnodr(3)
            zr1c       = -pnodr(1) - 3.e00*pnodr(3) + 4.e00*pnodr(5)
            zr1e       = -pnodr(2) - 3.e00*pnodr(3) + 4.e00*pnodr(4)
            zr2c       = 2.e00*(pnodr(1)+pnodr(3))  - 4.e00*pnodr(5)
            zrce       = 4.e00*(pnodr(3)-pnodr(4)-pnodr(5)+pnodr(6))
            zr2e       = 2.e00*(pnodr(2)+pnodr(3))  - 4.e00*pnodr(4)
            zz3        = pnodz(3)
            zz1c       = -pnodz(1) - 3.e00*pnodz(3) + 4.e00*pnodz(5)
            zz1e       = -pnodz(2) - 3.e00*pnodz(3) + 4.e00*pnodz(4)
            zz2c       = 2.e00*(pnodz(1)+pnodz(3))  - 4.e00*pnodz(5)
            zzce       = 4.e00*(pnodz(3)-pnodz(4)-pnodz(5)+pnodz(6))
            zz2e       = 2.e00*(pnodz(2)+pnodz(3))  - 4.e00*pnodz(4)
            csi        = 0.33333e00   
            eta        = 0.33333e00  
      do 100       j10 = 1,100
            zfr        = -yrp + zr3 + zr1c*csi + zr1e*eta           +
     +                   zr2c*csi*csi + zrce*csi*eta + zr2e*eta*eta 
            zfz        = -yzp + zz3 + zz1c*csi + zz1e*eta           +
     +                   zz2c*csi*csi + zzce*csi*eta + zz2e*eta*eta 
            zdfrc      = zr1c + 2.e00*zr2c*csi + zrce*eta
            zdfre      = zr1e + zrce*csi + 2.e00*zr2e*eta   
            zdfzc      = zz1c + 2.e00*zz2c*csi + zzce*eta
            zdfze      = zz1e + zzce*csi + 2.e00*zz2e*eta  
            zdelta     = zdfrc*zdfze - zdfre*zdfzc
            zdcsi      = (zfz*zdfre - zfr*zdfze)/zdelta     
            zdeta      = (zfr*zdfzc - zfz*zdfrc)/zdelta   
            ztest      = sqrt(zdcsi*zdcsi + zdeta*zdeta + 1.e-36)  
                         if(ztest.le.ztol) then
                         go to 99
                                           endif
            csi       = csi + zdcsi
            eta       = eta + zdeta
  100 continue
c     ++++++
      RETURN
c     ++++++
   99                   continue      
            csieta    = csi + eta
                        if(-1.e-03.le.csi.     and.
     +                     csi.le.1.001e00.    and.   
     +                     -1.e-03.le.eta.     and.
     +                     eta.le.1.001e00.    and.
     +                     -1.e-03.le.csieta.  and.
     +                     csieta.le.1.001e00         ) then
            zp3        = psiele(3)
            zp1c       = -psiele(1) - 3.e00*psiele(3) + 4.e00*psiele(5)
            zp1e       = -psiele(2) - 3.e00*psiele(3) + 4.e00*psiele(4)
            zp2c       = 2.e00*(psiele(1)+psiele(3))  - 4.e00*psiele(5)
            zpce       = 4.e00*(psiele(3)-psiele(4)-psiele(5)+psiele(6))
            zp2e       = 2.e00*(psiele(2)+psiele(3))  - 4.e00*psiele(4)
            ypsi       = zp3 + zp1c*csi + zp1e*eta                  +
     +                   zp2c*csi*csi + zpce*csi*eta + zp2e*eta*eta 
            kerror     = -1   
                                                        endif
c                                                                       
c-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c                                                                       
      RETURN                                                            
      END                                                               
      SUBROUTINE BPOLRZ                          !  RECTANGLE.AAMAIN
c     =================  
c                                                                       
c-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c                                                                       
c     Computes the two (-R- & -Z-) components of Bpoloidal.
c
c     Input quantities:  csi,eta,yrp,yzp,pnodr(*),pnodz(*),psiele(*)
c     Output quantities: ybpolr, ybpolz
c                                                                       
c-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c                                                                       
#include "pt2dg.inc"                                              
c
      dimension          zf(6),          zd(2,6),        zjac(2,2)
c                                                                
c-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c
c
      do 100        j1 = 1,6
                zf(j1) = 0.e00                                          
              zd(1,j1) = 0.e00                                          
              zd(2,j1) = 0.e00                                          
  100 continue 
c
c                                                         
            zl1        = csi 
            zl2        = eta  
            zl3        = 1.e00 - zl1 - zl2  
            zf(1)      = zl1*(2.e00*zl1-1.e00) 
            zf(2)      = zl2*(2.e00*zl2-1.e00)  
            zf(3)      = zl3*(2.e00*zl3-1.e00) 
            zf(4)      = 4.e00*zl2*zl3          
            zf(5)      = 4.e00*zl3*zl1         
            zf(6)      = 4.e00*zl1*zl2        
            zd(1,1)    = 4.e00*zl1 - 1.e00       
            zd(2,2)    = 4.e00*zl2 - 1.e00      
            zd(1,3)    = -(4.e00*zl3 - 1.e00)    
            zd(2,3)    = -(4.e00*zl3 - 1.e00)   
            zd(1,4)    = -4.e00*zl2              
            zd(2,4)    = 4.e00*(zl3 - zl2)     
            zd(1,5)    = 4.e00*(zl3 - zl1)    
            zd(2,5)    = -4.e00*zl1           
            zd(1,6)    = 4.e00*zl2            
            zd(2,6)    = 4.e00*zl1           
c
c
            zrx1       = 4.e00*pnodr(5)-3.e00*pnodr(3)-pnodr(1)      
            zrx2       = 2.e00*(pnodr(1)+pnodr(3)-2.e00*pnodr(5))    
            zrxy       = 4.e00*(pnodr(3)+pnodr(6)-                    
     +                          pnodr(4)-pnodr(5))  
            zjac(1,1)  = zrx1 + 2.e00*zrx2*csi + zrxy*eta  
c                
            zry1       = 4.e00*pnodr(4)-3.e00*pnodr(3)-pnodr(2)      
            zry2       = 2.e00*(pnodr(2)+pnodr(3)-2.e00*pnodr(4))    
            zjac(1,2)  = zry1 + zrxy*csi + 2.e00*zry2*eta
c
            zzx1       = 4.e00*pnodz(5)-3.e00*pnodz(3)-pnodz(1)      
            zzx2       = 2.e00*(pnodz(1)+pnodz(3)-2.e00*pnodz(5))    
            zzxy       = 4.e00*(pnodz(3)+pnodz(6)-                    
     +                          pnodz(4)-pnodz(5))  
            zjac(2,1)  = zzx1 + 2.e00*zzx2*csi + zzxy*eta    
c              
            zzy1       = 4.e00*pnodz(4)-3.e00*pnodz(3)-pnodz(2)      
            zzy2       = 2.e00*(pnodz(2)+pnodz(3)-2.e00*pnodz(4))    
            zjac(2,2)  = zzy1 + zzxy*csi + 2.e00*zzy2*eta
c
            zjdet      = zjac(1,1)*zjac(2,2) - zjac(2,1)*zjac(1,2)
c
c
            zpsir      = 0.e00
            zpsiz      = 0.e00
      do 200        j2 = 1,6
            zpsir      = zpsir + psiele(j2)*
     +                         (zjac(2,1)*zd(1,j2)-zjac(2,2)*zd(2,j2))
            zpsiz      = zpsiz - psiele(j2)*
     +                         (zjac(1,1)*zd(1,j2)-zjac(1,2)*zd(2,j2))
  200 continue   
            zpsir      = zpsir/zjdet 
            zpsiz      = zpsiz/zjdet
c
c
            ybpolr     = zpsiz/yrp 
            ybpolz     =-zpsir/yrp 
c                                                                       
c                                                                       
c-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c                                                                       
      RETURN                                                            
      END                                                               
      SUBROUTINE WBPOLRZ
c     ==================
c                                                                       
c-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c                                                                       
c     Writes the two (-R- & -Z-) components of Bpoloidal.
c
c-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c                                                                       
#include "pt2dg.inc"                                              
c                                                                
c-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c
      lun = nwreqdg 
c
      write(lun,*) 
      write(lun,*) '     ((bpr(j,k),j=1,jm),k=1,km)'
      write(lun,8000) ((bpr(j,k),j=1,mpr),k=1,mpz)
c
c                                                                       
      write(lun,*) 
      write(lun,*) '     ((bpz(j,k),j=1,jm),k=1,km)'
      write(lun,8000) ((bpz(j,k),j=1,mpr),k=1,mpz)
c                                                                       
c                                                                       
c-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c                                                                       
      RETURN    
c
 8000 format(5(3x,e14.8))
c                                                        
      END                                                               
c
c
c      subroutine wreqdg(lun,ngpr,ngpz,iret, nr,nz, psib,
c     ,                                             btf,rtf,rgr,zgr,pfm)
cc=====================================================
cc*** Write the equilibrium data in the dg-compatible format.
cc***
cc*** Input:
cc***  lun     the logical unit number for the output
cc***  ngpr    the maximum number of points in R direction
cc***  nr      the actual number of points in R direction
cc***  nz      the actual number of points in Z direction
cc***  psib    the poloidal flux at the separatrix
cc***  btf     the toroidal magnetic field at the R=rtf
cc***  rgr     the R values for the grid points
cc***  zgr     the Z values for the grid points
cc***  pfm     the values of the poloidal flux
cc***
cc*** Output:
cc***  iret    return code (0 means OK)
cc=====================================================
cc
cc  version : 23.06.97 17:23
cc
c      dimension rgr(ngpr), zgr(ngpz), pfm(ngpr,ngpz)
cc... toroidal field in tesla, radius in m
cc         btf, rtf, psib
cc
cc***  Write the toroidal field and the corresponding radius...
cc*** -- now obsolette!
cc
cc      write(3,*,err=99) 'Toroidal field in Tesla, radius in m'
cc      write(3,*,err=99) btf, rtf
cc
cc***  ... then the plasma equilibrium ...
cc
c      iret=0
c      write(lun,*,err=99)
c     /    '   jm   :=  no. of grid points in radial direction;'
c      write(lun,*,err=99)
c     /    '   km   :=  no. of grid points in vertical direction;'
c      write(lun,*,err=99)
c     /    '   r    :=  radial   co-ordinates of grid points  [m];'
c      write(lun,*,err=99)
c     /    '   z    :=  vertical co-ordinates of grid points  [m];'
c      write(lun,*,err=99)
c     /    '   psi  :=  flux per radiant at grid points     [Wb/rad];'
c      write(lun,*,err=99)
c     /    '   psib :=  psi at plasma boundary              [Wb/rad];'
c      write(lun,*,err=99)
c     /    '   btf  :=  toroidal magnetic field                  [t];'
c      write(lun,*,err=99)
c     /    '   rtf  :=  major radius at which btf is specified   [m];'
c      write(lun,*,err=99)
c      write(lun,*,err=99)
c      write(lun,*,err=99) '   jm    = ', nr,';'
c      write(lun,*,err=99) '   km    = ', nz,';'
c      write(lun,*,err=99) '   psib  = ',psib,' Wb/rad;'
c      write(lun,*,err=99) '   btf   = ',btf,' t;'
c      write(lun,*,err=99) '   rtf   = ',rtf,' m;'
c      write(lun,*,err=99)
c      write(lun,*,err=99) '   r(1:jm);'
c      write(lun,8000) (rgr(i),i=1,nr)
c      write(lun,*)
c      write(lun,*,err=99) '   z(1:km);'
c      write(lun,8000) (zgr(i),i=1,nz)
c      write(lun,*)
c      write(lun,*) '     ((psi(j,k)-psib,j=1,jm),k=1,km)'
c      write(lun,8000) ((pfm(i,j)-psib,i=1,nr),j=1,nz)
c 8000 format(5(3x,e14.8))
c      iret=0
c      return
cc-----------------------------------------------------
cc
c 99   print *,'==== wreqdg: error writing the files'
c      iret=8
cc
c      end
