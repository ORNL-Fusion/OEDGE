module mod_sol22_plots

  use mod_sol22_sources
  use mod_sol22_output


  implicit none



contains


  subroutine mkplot(dspts,npts,dte,dti,dne,dvb,xflag)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon

    !     Calls GHOST routines to plot the temperature and other results.

    !     include 'solparams'
    !     include 'solswitch'
    !     include 'solcommon'

    implicit none
    integer npts,xflag


    real*8 dspts(mxspts),dte(mxspts),dti(mxspts),dne(mxspts),dvb(mxspts)
    real  spts(mxspts),te(mxspts),ti(mxspts),ne(mxspts),vb(mxspts)
    integer i,ibrok,l,ispot,ncases
    character*50 xlabel,ylabel
    character*36 table,ref(4)
    real plots(mxspts,maxcases),scales(maxcases)
    real cxmin,cxmax,cymin,cymax


    !     local variables


    !real fmax,fmin
    !external fmax,fmin
    integer lenstr
    external lenstr

    !     Initialization

    do i = 1, mxspts
       spts(i) = dspts(i)
       te(i) = dte(i)
       ti(i) = dti(i)
       ne(i) = dne(i)
       vb(i) = dvb(i)

    end do
    ncases = 4

    ISPOT  = 12
    scales(1) = fmax1d(te,npts) ! slmod fmax(te,npts,1)
    scales(2) = fmax1d(ti,npts) ! slmod fmax(ti,npts,1)
    scales(3) = fmax1d(ne,npts) ! slmod fmax(ne,npts,1)

    scales(4) = abs(fmin1d(vb,npts)) ! slmod abs(fmin(vb,npts,1))
    xlabel = 'Distance along the field line (m)'
    ylabel = 'Normalized Quantities'

    table = 'Table of Values'
    write(ref(1),'(''Scale: Te = '',g13.6)') scales(1)
    write(ref(2),'(''Scale: Ti = '',g13.6)') scales(2)
    write(ref(3),'(''Scale: Ne = '',g13.6)') scales(3)

    !     draw titles

    write(ref(4),'(''Scale: Vb = '',g13.6)') -scales(4)

    !     Draw Frame

    call drawtitles(title,table,xlabel,ylabel)

    !     draw scales

    call drawframe
    cymin = -1.0
    cymax =  1.0
    cxmin = 0.0
    if (xflag.eq.0) then
       cxmax = fmax1d(spts,npts) ! fmax(spts,npts,1)
    elseif (xflag.eq.1)  then
       cxmax = graphran

    endif


    !     draw plots - load data into plot array

    call drawscales (cxmin,cxmax,cymin,cymax)
    do i = startn,npts
       plots(i,1) = te(i)
       plots(i,2) = ti(i)
       plots(i,3) = ne(i)
       plots(i,4) = vb(i)

    end do

    call plotdata(spts,plots,npts,ref,ncases,cxmin,cymin,cxmax,cymax,scales)
    return



  end subroutine mkplot

  subroutine mkauxplot(dspts,npts,dte,dti,dne,dvb,dphelpiv,dpeiv,dpeiv2,dpcxv,dpradv,xflag)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon

    !     This routine calculates the various quantities that go
    !     into calculating Te, Ti, n,v and plots them as a function of
    !     S.

    !     include 'solparams'
    !     include 'solswitch'
    !     include 'solcommon'

    implicit none
    integer npts,xflag
    real*8 dspts(mxspts),dte(mxspts),dti(mxspts),dne(mxspts),dvb(mxspts)
    real*8 dphelpiv(mxspts),dpeiv(mxspts),dpcxv(mxspts)

    !     local variables

    real*8 dpeiv2(mxspts),dpradv(mxspts)
    real spts(mxspts),te(mxspts),ti(mxspts),ne(mxspts),vb(mxspts)
    real phelpiv(mxspts),peiv(mxspts),pcxv(mxspts)
    real peiv2(mxspts)
    integer i,ibrok,l,ispot
    character*50 xlabel,ylabel
    character*36 table,ref1,ref2,ref3,ref4,ref5,ref6
    character*112 comment
    integer lenstr

    external lenstr

    real cxmin,cxmax,cymin,cymax
    real condi(mxspts),conde(mxspts),conv1i(mxspts)
    real conv2i(mxspts),convi(mxspts),conve(mxspts)
    real pradv(mxspts)
    real powi(mxspts),powe(mxspts)
    real sume(mxspts), sumi(mxspts)
    real totci(mxspts),totpi(mxspts)
    real totce(mxspts),totpe(mxspts)
    real*8 s,t1i,t1e,n

    !real*8 fegrad, figrad,pcx,pei,phelpi,conv,cond
    !external fegrad,figrad,pcx,pei,phelpi,conv,cond

    !real*8 paes,pais
    !external paes,pais

    !     Copy everything into local variables

    do i = 1,mxspts
       spts(i)   = dspts(i)
       te(i)     =  dte(i)
       ti(i)     =  dti(i)
       ne(i)     =  dne(i)
       vb(i)     =  dvb(i)
       phelpiv(i)=  dphelpiv(i)
       peiv(i)   =  dpeiv(i)
       peiv2(i)  =  dpeiv2(i)
       pcxv(i)   =  dpcxv(i)
       pradv(i)  =  dpradv(i)

       !     Calculate all the auxiliary values - at each S-point
       !     And the scale for the plot.

    end do
    cymin =  1.0e10
    cymax = -1.0e10


    ISPOT  = 12
    do i = startn,npts
       s   = spts(i)
       t1e = te(i)
       t1i = ti(i)


       n   = ne(i)
       if (actswcond.eq.0.0) then
          conve(i) = 0.0
          conv1i(i) = 0.0
       else
          conve(i) = cond(s,t1e)
          conv1i(i) = cond(s,t1i)

       endif
       if (actswconv.eq.0.0) then
          conv2i(i) = 0.0
       else
          conv2i(i) = conv(s,n,t1i)

          !         if (actswprad.eq.0.0) then
          !            pradv(i)  = 0.0
          !         else
          !            pradv(i)  = prad(s)
          !         endif

       endif
       powe(i) = paes(s)


       powi(i) = pais(s)

       convi(i) = conv1i(i) + conv2i(i)

       conde(i)= -(conve(i)+pradv(i)+phelpiv(i)+peiv(i)+powe(i))

       condi(i)  = - (convi(i)+ pcxv(i)-peiv(i)+powi(i))
       sume(i) = conde(i) + conve(i) + pradv(i)+ phelpiv(i) + peiv(i)

       sumi(i) = condi(i) + convi(i)+ pcxv(i) - peiv(i)
       totci(i) = condi(i) + convi(i)

       totce(i) = conde(i) + conve(i)
       totpi(i) = pcxv(i) - peiv(i) + powi(i)

       totpe(i) = phelpiv(i) + peiv(i) + powe(i) + pradv(i)



    end do
    if (actswpei.eq.3) then
       do i = startn,npts
          peiv(i) = peiv2(i)
       end do

    endif

    !     Print out tables of Auxiliary values

    if (xflag.eq.0) then
       if (actswpei.eq.3.0) then
          call prbs
          call prs('Pei values are estimates NOT included in Totals:')

       endif
       call prbs
       call prs('Tables of calculated SOL Auxiliary values - E')
       call prbs
       write(comment,200)
       call prs(comment)
       do i = startn, npts
          write(comment,100) spts(i),conde(i),conve(i),pradv(i),phelpiv(i),peiv(i),sume(i),powe(i),totce(i),totpe(i)
          call prs (comment)

       end do
       call prbs
       call prs('Tables of calculated SOL Auxiliary values - I')
       call prbs
       write(comment,300)
       call prs(comment)
       do i = startn, npts
          write(comment,100) spts(i),condi(i),conv1i(i),conv2i(i), pcxv(i),-peiv(i),sumi(i),powi(i),totci(i), totpi(i)
          call prs (comment)

       end do



    endif
100 format (g10.3,9g11.4)
200 format (5x,'S',9x,'Cond',7x,'Conv',7x,'Prad',7x,'Phelp',6x,'Pei',8x,'Sum E',6x,'Pow e',5x,'Cond+Conv',4x,'Srcs')


    !     Set up the graph. E-power balance.


300 format (5x,'S',9x,'Cond',7x,'Conv1',5x,'Kinetic',5x,'Pcx',8x,'Pei',8x,'Sum I',6x,'Pow I',5x,'Cond+Conv',4x,'Srcs')
    xlabel = 'Distance along the field line (m)'
    ylabel = 'E-power balance quantities'

    table = 'Table of Values'
    write(ref1,'(''Conduction term '')')
    write(ref2,'(''5/2nvkTe        '')')
    write(ref3,'(''Prad            '')')
    write(ref4,'(''Phelpi          '')')
    write(ref5,'(''Pei             '')')

    !     calculate cymin and cymax ... if they are the same - skip plotting


    !     Find Y range for the E-power plot

    write(ref6,'(''Sum of all terms'')')
    cymin =  1.0e30

    cymax = -1.0e30

    call fminmax(cymin,cymax,conde,npts)

    if (actswcond.ne.0.0)call fminmax(cymin,cymax,conve,npts)

    if (actswprad.ne.0.0)call fminmax(cymin,cymax,pradv,npts)

    if (actswpei.ne.0.0)call fminmax(cymin,cymax,peiv,npts)

    if (actswphelp.ne.0.0)call fminmax(cymin,cymax,phelpiv,npts)


    call fminmax(cymin,cymax,sume,npts)

    !---- DRAW TITLES

    if (cymin.lt.cymax) then
       CALL PSPACE (0.0, 1.35, 0.0, 1.0)
       CALL MAP    (0.0, 1.35, 0.0, 1.0)
       CALL CTRMAG (20)
       CALL LINCOL (1)
       CALL THICK  (2)
       L = LENSTR(TITLE)
       CALL PCSCEN (0.8, 0.95, TITLE(:L))
       L = LENSTR(TABLE)
       CALL PCSCEN (1.18, 0.87, TABLE(:L))
       L = LENSTR (XLABEL)
       CALL PCSCEN (0.5,0.025,XLABEL(:L))
       CALL THICK  (1)
       CALL CTRMAG (12)
       L = LENSTR (REF1)
       CALL plotst (1.06, 0.77, REF1(:L))
       L = LENSTR (REF2)
       CALL plotst (1.06, 0.72, REF2(:L))
       L = LENSTR (REF3)
       CALL plotst (1.06, 0.67, REF3(:L))
       L = LENSTR (REF4)
       CALL plotst (1.06, 0.62, REF4(:L))
       L = LENSTR (REF5)
       CALL plotst (1.06, 0.57, REF5(:L))
       L = LENSTR (REF6)
       CALL plotst (1.06, 0.52, REF6(:L))

       !---- DRAW FRAMES

       CALL CTRMAG (ISPOT)
       CALL LINCOL (3)
       CALL POSITN (0.1, 0.1)
       CALL   JOIN (0.1, 0.9)
       CALL   JOIN (0.9, 0.9)
       CALL   JOIN (0.9, 0.1)
       CALL   JOIN (0.1, 0.1)
       CALL POSITN (0.93, 0.1)
       CALL   JOIN (0.93, 0.9)
       CALL   JOIN (1.35, 0.9)
       CALL   JOIN (1.35, 0.1)
       CALL   JOIN (0.93, 0.1)
       CALL POSITN (0.93, 0.85)
       CALL   JOIN (1.35, 0.85)
       CALL POSITN (0.93, 0.15)
       CALL   JOIN (1.35, 0.15)
       CALL POSITN (0.93, 0.20)
       CALL   JOIN (1.35, 0.20)
       CALL POSITN (0.93, 0.25)

       !     Draw scales - xscale

       CALL   JOIN (1.35, 0.25)
       CXMIN = spts(1)
       if (xflag.eq.0) then
          cxmax = spts(npts)
       elseif (xflag.eq.1)  then
          cxmax = graphran

       endif
       CALL PSPACE (0.1, 0.9, 0.1, 0.9)
       CALL MAP    (CXMIN, CXMAX, 0.1, 0.9)
       CALL CTRMAG (10)

       !     yscale


       !     Set Y-scales

       CALL XSCALE
       CALL PSPACE (0.1, 0.9, 0.11, 0.89)
       CALL MAP    (0.0, 1.0, CYMIN, CYMAX)
       CALL LINCOL (3)
       CALL CTRMAG (10)
       CALL YSCALE
       CALL PSPACE (0.0, 1.35, 0.11, 0.89)
       CALL MAP    (0.0, 1.35, 0.0, 1.0)
       CALL CTRORI (90.0)
       CALL LINCOL (1)
       CALL CTRMAG (14)
       CALL THICK  (2)
       L = LENSTR (YLABEL)
       CALL PCSEND (0.02,0.99,YLABEL(:L))
       CALL CTRORI (0.0)

       CALL THICK  (1)

       !---- DRAW PLOTS - E-power balance

       !     Cond

       CALL FULL

       !     Conv

       call plotln(1,0.77,npts,spts,conde,cxmin,cxmax,cymin,cymax,1.0)

       !     Prad

       call plotln(2,0.72,npts,spts,conve,cxmin,cxmax,cymin,cymax,1.0)

       !     Phelpi

       call plotln(3,0.67,npts,spts,pradv,cxmin,cxmax,cymin,cymax,1.0)

       !     Pei

       call plotln(4,0.62,npts,spts,phelpiv,cxmin,cxmax,cymin,cymax,1.0)

       !     Sum of E-balance terms

       call plotln(5,0.57,npts,spts,peiv,cxmin,cxmax,cymin,cymax,1.0)

       call plotln(6,0.52,npts,spts,sume,cxmin,cxmax,cymin,cymax,1.0)

       !     End of E-power balance plot

       call frame



       !     Set up the graph. I-power balance.


       !     If cymin = cymax then skip the plot as there will be no range
       !     to plot.

       !     Find Y range for the plot

    endif
    cymin =  1.0e30

    cymax = -1.0e30

    call fminmax(cymin,cymax,condi,npts)

    if (actswcond.ne.0.0.or.actswconv.ne.0.0)call fminmax(cymin,cymax,convi,npts)

    if (actswpcx.ne.0.0)call fminmax(cymin,cymax,pcxv,npts)

    if (actswpei.ne.0.0)call fminmax(cymin,cymax,peiv,npts)

    !     Test cymin,cymax

    call fminmax(cymin,cymax,sumi,npts)


    if (cymin.lt.cymax) then
       xlabel = 'Distance along the field line (m)'
       ylabel = 'I-power balance quantities'
       ISPOT  = 12

       table = 'Table of Values'
       write(ref1,'(''Conduction term '')')
       write(ref2,'(''Convection term '')')
       write(ref3,'(''Pcx             '')')
       write(ref4,'(''Pei             '')')

       !---- DRAW TITLES

       write(ref5,'(''Sum of all terms'')')
       CALL PSPACE (0.0, 1.35, 0.0, 1.0)
       CALL MAP    (0.0, 1.35, 0.0, 1.0)
       CALL CTRMAG (20)
       CALL LINCOL (1)
       CALL THICK  (2)
       L = LENSTR(TITLE)
       CALL PCSCEN (0.8, 0.95, TITLE(:L))
       L = LENSTR(TABLE)
       CALL PCSCEN (1.18, 0.87, TABLE(:L))
       L = LENSTR (XLABEL)
       CALL PCSCEN (0.5,0.025,XLABEL(:L))
       CALL THICK  (1)
       CALL CTRMAG (12)
       L = LENSTR (REF1)
       CALL plotst (1.06, 0.77, REF1(:L))
       L = LENSTR (REF2)
       CALL plotst (1.06, 0.72, REF2(:L))
       L = LENSTR (REF3)
       CALL plotst (1.06, 0.67, REF3(:L))
       L = LENSTR (REF4)
       CALL plotst (1.06, 0.62, REF4(:L))
       L = LENSTR (REF5)
       CALL plotst (1.06, 0.57, REF5(:L))

       !---- DRAW FRAMES

       CALL CTRMAG (ISPOT)
       CALL LINCOL (3)
       CALL POSITN (0.1, 0.1)
       CALL   JOIN (0.1, 0.9)
       CALL   JOIN (0.9, 0.9)
       CALL   JOIN (0.9, 0.1)
       CALL   JOIN (0.1, 0.1)
       CALL POSITN (0.93, 0.1)
       CALL   JOIN (0.93, 0.9)
       CALL   JOIN (1.35, 0.9)
       CALL   JOIN (1.35, 0.1)
       CALL   JOIN (0.93, 0.1)
       CALL POSITN (0.93, 0.85)
       CALL   JOIN (1.35, 0.85)
       CALL POSITN (0.93, 0.15)
       CALL   JOIN (1.35, 0.15)
       CALL POSITN (0.93, 0.20)
       CALL   JOIN (1.35, 0.20)
       CALL POSITN (0.93, 0.25)

       !     Draw scales - xscale

       CALL   JOIN (1.35, 0.25)
       CXMIN = spts(1)
       if (xflag.eq.0) then
          cxmax = spts(npts)
       elseif (xflag.eq.1)  then
          cxmax = graphran

       endif
       CALL PSPACE (0.1, 0.9, 0.1, 0.9)
       CALL MAP    (CXMIN, CXMAX, 0.1, 0.9)
       CALL CTRMAG (10)

       !     yscale


       !     Set Y-scales

       CALL XSCALE
       CALL PSPACE (0.1, 0.9, 0.11, 0.89)
       CALL MAP    (0.0, 1.0, CYMIN, CYMAX)
       CALL LINCOL (3)
       CALL CTRMAG (10)
       CALL YSCALE
       CALL PSPACE (0.0, 1.35, 0.11, 0.89)
       CALL MAP    (0.0, 1.35, 0.0, 1.0)
       CALL CTRORI (90.0)
       CALL LINCOL (1)
       CALL CTRMAG (14)
       CALL THICK  (2)
       L = LENSTR (YLABEL)
       CALL PCSEND (0.02,0.99,YLABEL(:L))
       CALL CTRORI (0.0)

       CALL THICK  (1)

       !---- DRAW PLOTS -  I-Power Balance

       !     Cond

       CALL FULL

       !     Conv

       call plotln(1,0.77,npts,spts,condi,cxmin,cxmax,cymin,cymax,1.0)

       !     Pcx

       call plotln(2,0.72,npts,spts,convi,cxmin,cxmax,cymin,cymax,1.0)

       !     Pei

       call plotln(3,0.67,npts,spts,pcxv,cxmin,cxmax,cymin,cymax,1.0)

       !     Sum of I-balance terms

       call plotln(4,0.62,npts,spts,peiv,cxmin,cxmax,cymin,cymax,1.0)

       call plotln(5,0.57,npts,spts,sumi,cxmin,cxmax,cymin,cymax,1.0)

       !     End of I-power plot

       call frame


    endif
    return



  end subroutine mkauxplot
  subroutine mkvelplot(dspts,npts,dte,dti,dne,dvb,dne2,dga,xflag)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon

    !     This routine calculates the various quantities that go
    !     into calculating Te, Ti, n,v and plots them as a function of
    !     S.

    !     include 'solparams'
    !     include 'solswitch'
    !     include 'solcommon'

    implicit none
    integer npts,xflag
    real*8 dspts(mxspts),dte(mxspts),dti(mxspts),dne(mxspts),dvb(mxspts)

    !     local variables

    real*8 dne2(mxspts),dga(mxspts)
    real spts(mxspts),te(mxspts),ti(mxspts),ne(mxspts),vb(mxspts)
    real ne2(mxspts),ga(mxspts)
    character*50 xlabel,ylabel
    character*36 table,ref(3)

    character*78 comment

    !real fmin,fmax
    !external fmin,fmax

    integer lenstr,i
    external lenstr


    real cxmin,cxmax,cymin,cymax
    real vplots(mxspts,maxcases),scales(maxcases)

    !     Initialize


    !     Copy everything into local varaibles

    integer ncases
    do i = 1,mxspts
       spts(i)   =  dspts(i)
       te(i)     =  dte(i)
       ti(i)     =  dti(i)
       ne(i)     =  dne(i)
       vb(i)     =  dvb(i)
       ne2(i)    =  dne2(i)
       ga(i)     =  dga(i)

    end do
    ncases = 3
    do i = 1,ncases
       scales(i) = 0.0

       !     Calculate all the velocities (super and sub sonic at each S-point)
       !     And the scale for the plot.

    end do
    cymin =  1.0e10

    cymax = -1.0e10

    write (6,*) 'vsuper:',lastiters
    !         if ((spts(i).lt.lastiters).or.(spts(i).eq.0.0)) then
    do i = startn,npts
       if ((spts(i).lt.lastiters)) then
          if (ne2(i).gt.0.0) then
             vplots(i,1) = ga(i)/ne2(i)
          else
             vplots(i,1) = 0.0
          endif
          vplots(i,2) = vb(i)
       elseif (spts(i).ge.lastiters) then
          vplots(i,1) = vb(i)
          if (ne2(i).gt.0.0) then
             vplots(i,2) = ga(i)/ne2(i)
          else
             vplots(i,2) = 0.0
          endif
       endif
       vplots(i,3) = -sqrt((te(i)+ti(i))/mb * econv/mconv)

       !     Print

    end do

    !     Print out tables of the velocity values

    if (xflag.eq.0) then
       call prbs
       call pris('Tables of calculated SOL Velocity values ',ringnum)
       call prbs
       write(comment,300)
       call prs(comment)
       do i = startn, npts
          write(comment,100) spts(i),vplots(i,1),vplots(i,2),vplots(i,3),ga(i),ne(i)
          call prs(comment)

       end do



    endif
100 format (6g13.5)

    !     Set up the graph. Velocities


300 format (6x,'S',11x,'Vsub',10x,'Vsuper',9x,'Cs',8x,'Gamma',9x,'Ne')
    xlabel = 'Distance along the field line (m)'
    ylabel = 'Velocity in m/s'

    !     draw titles

    table = 'Table of Values'

    !     Draw Frame

    call drawtitles(title,table,xlabel,ylabel)

    !     Calculate Min and Max values

    call drawframe
    cxmin = 0.0


    !      cxmax = fmax(spts,npts,1)

    cymax = 0.0
    if (xflag.eq.0) then

       !         cxmax = spts(npts)

       cxmax = 10.0 * lastiters
    elseif (xflag.eq.1)  then
       cxmax = graphran

       !      cymin = fmin(vplots,npts,ncases)


       !      Need to establish a more usable value of cymin ... perhaps
       !      four times the minimum value of Vb

    endif


    !     draw scales

    cymin= 4.0 * fmin1d(vb,npts) ! slmod     cymin= 4.0 * fmin(vb,npts,1)

    call drawscales (cxmin,cxmax,cymin,cymax)
    write(ref(1),'(''Sub-sonic branch'')')
    write(ref(2),'(''Sup-sonic branch'')')
    !     draw plots

    write(ref(3),'(''Sound Speed (Cs)'')')


    call plotdata(spts,vplots,npts,ref,ncases,cxmin,cymin,cxmax,cymax,scales)
    return



  end subroutine mkvelplot



  subroutine plotln(ibrok,spot,npts,x,y,xmin,xmax,ymin,ymax,scale)

    !     This subroutine draws a curve on the graph.

    implicit none
    integer ibrok,npts,i

    real spot,x(npts),y(npts),xmin,xmax,ymin,ymax,scale
    CALL PSPACE (0.1, 0.9, 0.11, 0.89)

    CALL MAP    (XMIN,XMAX,YMIN,YMAX)

    CALL BROKEN (3*IBROK,2*IBROK,3*IBROK,2*IBROK)

    CALL POSITN (x(1), y(1)/scale)
    DO  I = 2,npts
       CALL JOIN (x(i), y(i)/scale)

    end do
    CALL LINCOL (2)
    CALL CTRMAG (12)
    CALL PSPACE (0.0, 1.35, 0.0,1.0)
    CALL CSPACE (0.0, 1.35, 0.0,1.0)
    CALL MAP    (0.0, 1.35, 0.0,1.0)
    CALL POSITN (0.95,SPOT)
    CALL JOIN   (1.03,SPOT)
    CALL FULL
    return



  end subroutine plotln



  subroutine plotdata(sdata,pvalues,pts,cases,ncases,cxmin,cymin,cxmax,cymax,scales)
    use mod_solparams
    use mod_solcommon

    !     This routine plots the data in the two dimensional array.

    !     include 'solparams'
    !     include 'solcommon'

    implicit none
    integer ncases,pts
    real cxmin,cxmax,cymin,cymax
    real sdata(mxspts)
    real pvalues(mxspts,maxcases),scales(maxcases)

    character*(*) cases(ncases)
    integer i,j

    real values(mxspts),startpos,diffpos,plotpos
    startpos = 0.77

    diffpos = 0.05
    do i = 1,ncases
       plotpos= startpos - (i-1) * diffpos
       do j = 1,pts
          values(j) = pvalues(j,i)
       enddo
       call plotvals(i,plotpos,pts,sdata,values,cxmin,cxmax,cymin,cymax,cases(i),scales(i))

    enddo

    call frame
    return



  end subroutine plotdata



  subroutine plotvals(ibrok,spot,npts,x,y,xmin,xmax,ymin,ymax,name,scale)

    !     This subroutine draws a curve on the graph.

    implicit none
    integer ibrok,npts,i
    real spot,x(npts),y(npts),xmin,xmax,ymin,ymax,scale
    character*(*) name
    integer lenstr,l

    external lenstr
    CALL PSPACE (0.1, 0.9, 0.11, 0.89)

    CALL MAP    (XMIN,XMAX,YMIN,YMAX)

    CALL BROKEN (3*IBROK,2*IBROK,3*IBROK,2*IBROK)
    if (scale.eq.0.0) then

       CALL POSITN (x(1), y(1))
       DO  I = 2,npts
          CALL JOIN (x(i), y(i))
       end do
    elseif (scale.ne.0.0) then

       CALL POSITN (x(1), y(1)/scale)
       DO  I = 2,npts
          CALL JOIN (x(i), y(i)/scale)
       end do

    endif
    CALL LINCOL (2)
    CALL CTRMAG (12)
    CALL PSPACE (0.0, 1.35, 0.0,1.0)
    CALL CSPACE (0.0, 1.35, 0.0,1.0)
    CALL MAP    (0.0, 1.35, 0.0,1.0)
    CALL POSITN (0.95,SPOT)
    CALL JOIN   (1.03,SPOT)
    CALL FULL
    CALL THICK  (1)
    CALL CTRMAG (12)
    L = LENSTR (name)
    CALL plotst (1.06, spot, name(:L))

    call full
    return



  end subroutine plotvals



  subroutine drawscales(cxmin,cxmax,cymin,cymax)

    !     Draw scales

    implicit none

    !     xscale


    real cxmin,cxmax,cymin,cymax
    CALL PSPACE (0.1, 0.9, 0.1, 0.9)
    CALL MAP    (CXMIN, CXMAX, 0.1, 0.9)
    CALL CTRMAG (10)

    !     yscale

    CALL XSCALE
    CALL PSPACE (0.1, 0.9, 0.11, 0.89)
    CALL MAP    (0.0, 1.0, CYMIN, CYMAX)
    CALL LINCOL (3)
    CALL CTRMAG (10)

    CALL YSCALE
    CALL FULL
    return



  end subroutine drawscales
  subroutine drawframe

    !     Draws plot frame work

    implicit none
    CALL PSPACE (0.0, 1.35, 0.0, 1.0)
    CALL MAP    (0.0, 1.35, 0.0, 1.0)
    call thick(1)
    call ctrmag(12)
    CALL LINCOL (3)
    CALL POSITN (0.1, 0.1)
    CALL   JOIN (0.1, 0.9)
    CALL   JOIN (0.9, 0.9)
    CALL   JOIN (0.9, 0.1)
    CALL   JOIN (0.1, 0.1)
    CALL POSITN (0.93, 0.1)
    CALL   JOIN (0.93, 0.9)
    CALL   JOIN (1.35, 0.9)
    CALL   JOIN (1.35, 0.1)
    CALL   JOIN (0.93, 0.1)
    CALL POSITN (0.93, 0.85)
    CALL   JOIN (1.35, 0.85)
    CALL POSITN (0.93, 0.15)
    CALL   JOIN (1.35, 0.15)
    CALL POSITN (0.93, 0.20)
    CALL   JOIN (1.35, 0.20)
    CALL POSITN (0.93, 0.25)


    CALL   JOIN (1.35, 0.25)
    return



  end subroutine drawframe

  subroutine drawtitles(title,table,xlabel,ylabel)

    !     Draw labels

    implicit none
    integer l,lenstr
    external lenstr

    character*(*) title,xlabel,ylabel,table
    CALL PSPACE (0.0, 1.35, 0.0, 1.0)
    CALL MAP    (0.0, 1.35, 0.0, 1.0)
    CALL CTRMAG (20)
    CALL LINCOL (1)
    CALL THICK  (2)
    L = LENSTR(TITLE)
    CALL PCSCEN (0.8, 0.95, TITLE(:L))
    L = LENSTR(TABLE)
    CALL PCSCEN (1.18, 0.87, TABLE(:L))
    L = LENSTR (XLABEL)
    CALL PCSCEN (0.5,0.025,XLABEL(:L))
    CALL THICK  (1)

    !     Ylabel

    CALL CTRMAG (12)
    CALL PSPACE (0.0, 1.35, 0.11, 0.89)
    CALL MAP    (0.0, 1.35, 0.0, 1.0)
    CALL CTRORI (90.0)
    CALL LINCOL (1)
    CALL CTRMAG (14)
    CALL THICK  (2)
    L = LENSTR (YLABEL)
    CALL PCSEND (0.02,0.99,YLABEL(:L))
    CALL CTRORI (0.0)

    CALL THICK  (1)
    return

  end subroutine drawtitles

  ! slmod begin
  real function fmin1d(values,pts)
    use mod_solparams
    !     Finds minimum in array
    implicit none
    real values(mxspts)
    integer pts
    integer i,j
    fmin1d = 1.0e30
    !do i = 1,ncases
       do j = 1,pts
          fmin1d = amin1(fmin1d,values(j))!values(j,i))
       end do
    !end do
    return
  end function fmin1d
  ! slmod end
  real function fmin(values,pts,ncases)
    use mod_solparams

    !     Finds minimum in array

    !     include 'solparams'

    implicit none
    real values(mxspts,maxcases)

    integer pts,ncases

    integer i,j
    fmin = 1.0e30
    do i = 1,ncases
       do j = 1,pts
          fmin = amin1(fmin,values(j,i))
       end do
    end do
    return



  end function fmin
  ! slmod begin
  real function fmax1d(values,pts)
    use mod_solparams
    !     Finds maximum in array
    implicit none
    real values(mxspts)
    integer pts,ncases
    integer i,j
    fmax1d = -1.0e30
    !do i = 1,ncases
       do j = 1,pts
          fmax1d = amax1(fmax1d,values(j))
       end do
    !end do
    return
  end function fmax1d
  ! slmod end  

  real function fmax(values,pts,ncases)
    use mod_solparams

    !     Finds maximum in array

    !     include 'solparams'

    implicit none
    real values(mxspts,maxcases)

    integer pts,ncases

    integer i,j
    fmax = -1.0e30
    do i = 1,ncases
       do j = 1,pts
          fmax = amax1(fmax,values(j,i))
       end do
    end do
    return



  end function fmax


  subroutine fminmax(ymin,ymax,vals,npts)

    !     Finds the maximum and minimum values in the
    !     array vals ... if they are greater/less than the
    !     values passed in in ymin,ymax.

    implicit none
    integer npts,i

    real ymin,ymax,vals(npts)
    do i = 1,npts
       ymin = amin1(ymin,vals(i))
       ymax = amax1(ymax,vals(i))
    enddo
    return

  end subroutine fminmax


end module mod_sol22_plots
