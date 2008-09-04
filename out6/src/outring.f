      subroutine plot_ring(iselect,istate,
     >                  iexpt,minfrac,maxfrac,axis_type,
     >                  plot_type,
     >                  iopt,job,title,table,nplots,
     >                  iplot,nizs,
     >                  ierr) 
      implicit none
c
      include 'params'
      include 'cgeom'
c
      integer iselect,istate,iexpt,ierr 
      real minfrac,maxfrac
      integer axis_type,plot_type
      character*(*) job,title,table
      integer iopt,nplots,iplot,nizs
c
c     PLOT_RING: This routine plots a generalized ALONG RING plot. 
c
c     Local variables
c
      real tmpplot(maxnks,maxnrs)
c
c     Along ring variables
c
      real,allocatable::tvals(:,:),touts(:),twids(:)
c
c     Plot labels   
c
      character*36 blabs(2),xlab,ylab
      character*44 blabd,refd,refe,plane,anly,nview
c
      integer ik,ir,in
      integer flag,targ_flag,nvals,ndata,offset_n
      real plotmin,plotmax,tmpval
      real start_targ_val,end_targ_val
c
      INTEGER IGNORS(1),ITEC,NAVS
      REAL AVS(0:100)
c
c     Initialization
c
      ierr  = 0     
      targ_flag = 0
      navs = 0
      ignors(1) = 1
c
      refd  = ' '
      blabd = ' '  
c
c     Loading DIV DATA:
c
c
c     Do not load DIVIMP data if ISELECT is zero - only plotting experimental 
c     data. 
c
      if (iselect.ne.0) then 
         call load_divdata_array(tmpplot,iselect,istate,0,
     >                         ylab,blabd,refd,nizs,ierr)
      endif 
c
c     If there is an error loading the DIVIMP data to plot then 
c     exit
c
      if (ierr.ne.0) return
c
c     Load target data if any is available 
c     - targ_flag is 0 if data available
c
      call load_divdata_targ(iselect,istate,ir,
     >            start_targ_val,end_targ_val,targ_flag)
c
c     Allocate the storage arrays 
c
      ir = iopt
      nvals = nks(ir) 
      ndata = 1 
      offset_n = 0
c
c     Adjust values if target data is to be included 
c 
      if (targ_flag.eq.0) then 
         offset_n = 1
         nvals = nvals + 2 
      endif
c 
c     Allocate plot arrays
c
      allocate(tvals(nvals,ndata),stat=flag)
      allocate(touts(nvals),stat=flag)
      allocate(twids(nvals),stat=flag)
c
      if (flag.ne.0) then 
         write(6,'(a,i5)') 'ERROR IN PLOT_RING:'
     >               //' CAN NOT ALLOCATE STORAGE: ',ierr
         write(0,'(a,i5)') 'ERROR IN PLOT_RING:'
     >               //' CAN NOT ALLOCATE STORAGE: ',ierr
         return
      endif
c
c     Set XLAB and assign axis ends if targets are included
c
c     VS. S
c 
      if (axis_type.eq.0.or.axis_type.eq.2) then 
c
         XLAB = 'S (M)'
c
         plotmin = minfrac * ksmaxs(ir)
         plotmax = maxfrac * ksmaxs(ir)
c
         if (targ_flag.eq.0) then 
            touts(1) = 0.0
            touts(nvals) = ksmaxs(ir)
            twids(1) = 1.0
            twids(nvals) = 1.0
         endif
c
c     VS. P
c
      elseif (axis_type.eq.1.or.axis_type.eq.3) then 
c
         XLAB = 'P (M)'
c
         plotmin = minfrac * kpmaxs(ir)
         plotmax = maxfrac * kpmaxs(ir)
c
         if (targ_flag.eq.0) then 
            touts(1) = 0.0
            touts(nvals) = kpmaxs(ir)
            twids(1) = 1.0
            twids(nvals) = 1.0
         endif
c
      endif
c
c     Assign ring values to plot variables      
c
      if (targ_flag.eq.0) then 
         tvals(1,ndata) = start_targ_val
         tvals(nvals,ndata) = end_targ_val   
      endif 
c
      do ik = 1,nks(ir)
c
         tvals(ik+offset_n,ndata) = tmpplot(ik,ir) 
c
         if (axis_type.eq.0.or.axis_type.eq.2) then 
            touts(ik+offset_n) = kss(ik,ir)
            twids(ik+offset_n) = ksb(ik,ir)-ksb(ik-1,ir)
         elseif (axis_type.eq.1.or.axis_type.eq.4) then 
            touts(ik+offset_n) = kps(ik,ir)
            twids(ik+offset_n) = kpb(ik,ir)-kpb(ik-1,ir)
         endif
c
      enddo
c
c     Invert axes for axis type 2 or 3 
c     and rescale to zero at low end
c
      if (axis_type.eq.2.or.axis_type.eq.3) then  
c
         do ik = 1,nvals/2
c
c           Swap data for each quantity 
c 
c           Values
c
            do in = 1,ndata
               tmpval = tvals(ik,in)
               tvals(ik,in) = tvals(nvals-ik+1,in) 
               tvals(nvals-ik+1,in) = tmpval 
            end do
c
c           Axis
c
            tmpval = touts(ik)
            touts(ik) = touts(nvals-ik+1) 
            touts(nvals-ik+1) = tmpval 
c
c           Widths
c
            tmpval = twids(ik)
            twids(ik) = twids(nvals-ik+1) 
            twids(nvals-ik+1) = tmpval 
c
         end do 
c
c        Now adjust axis so it scales from zero again
c
         do ik = 1,nvals
            if (axis_type.eq.2) then 
               touts(ik) = ksmaxs(ir) - touts(ik)
            elseif (axis_type.eq.3) then 
               touts(ik) = kpmaxs(ir) - touts(ik)
            endif
         end do
c
c        Adjust plotting range - note that meaning of min and max 
c        changes
c
         if (axis_type.eq.2) then 
            plotmin = (1.0-maxfrac) * ksmaxs(ir)
            plotmax = (1.0-minfrac) * ksmaxs(ir) 
         elseif (axis_type.eq.3) then 
            plotmin = (1.0-maxfrac) * kpmaxs(ir)
            plotmax = (1.0-minfrac) * kpmaxs(ir) 
         endif
c
      endif
c
c     Write plot data
c
c      write(6,*) 'NVALS:',nvals,plotmin,plotmax
c      do ik = 1,nvals
c         write(6,'(i5,3(1x,g12.5))') ik,twids(ik),touts(ik),
c     >             tvals(ik,1)
c      end do
c
c     Set rest of labels
c
c
c     Set NVIEW, PLANE ...
c
      NVIEW  = 'GENERALIZED ALONG RING PLOT'
      write(plane,'(a,i5)')  'ALONG RING: ',ir
c
c     Draw plot
c
      CALL DRAW (tOUTS,tWIDS,tVALS,nvals,nvals,ANLY,
     >    ndata,99,plotmin,plotmax,-HI,HI,IGNORS,ITEC,
     >    AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,BLABd,REFd,NVIEW,PLANE,
     >    TABLE,plot_type,2,1.0,iexpt)
c
c 
c     De-Allocate plot arrays
c
      deallocate(tvals)
      deallocate(touts)
      deallocate(twids)
c
      return
      end 






