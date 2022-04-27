      subroutine setup_signal_output
      use mod_signal_com
      implicit none
c     include 'signal_com'
      integer ierr
      open(unit=signal_unit,file='signal_output.dat',
     >     form='formatted',access='sequential',iostat=ierr)
      

      return 
      end


      subroutine write_signal(as,cs,ia1,ia2,mxxnas,nbs,igs,title,
     >     blabs,aaxlab,baxlab,ref,view,plane,anly,job)
      use mod_signal_com
      implicit  none
c     include 'signal_com'
      INTEGER   MXXNAS,NBS,IGS(*),ia1,ia2
      REAL      AS(*),CS(MXXNAS,*)
      CHARACTER TITLE*80
      character job*72
      character aaxlab*36
      character baxlab*(*)
      CHARACTER BLABS(nbs)*(*)
      character REF*(*),VIEW*(*),PLANE*(*),anly*(*)
c     
c     This routine prints the values that are plotted in the graph in a
c     two column format ... the first column is the "X" axis and the
c     second is the "Y" axis. These columns of data can then be
c     extracted and plotted using a spreadsheet program. The purpose
c     of this is to allow some "easy" collation of results from
c     different cases without incurring a lot of overhead since
c     one would expect to need this type of collation primarily
c     for presentations and reports.
c     
c     Note!: the routine plot_expt now also adds data to the column data formatted plot file - as a 
c     result iout must be kept the same in both routines. 
c     
      integer in,ia,ib,icnt

      character*512 line,tag,label,units
      character*36 val
      character*1 delim

      icnt = 0
      do ib = 1,nbs
         if (igs(ib).gt.0) then 
            icnt = icnt+1
         endif
      enddo 

      delim = ','


      tag = 'SIGNAL:'
      write(signal_unit,'(a,1x,a)') trim(tag),trim(job(1:36))
      

      tag = 'LABELS:'
      label = trim(tag)
      label = trim(label) // trim(aaxlab)
      tag='UNITS:'
      units = trim(tag) // trim(aaxlab)
c      write(0,*) 'Label:',trim(label),':'

      do ib = 1,nbs
c         write(0,*) 'LabelB:',ib,':',trim(blabs(ib)),':'
c         write(0,*) 'LabelB:',ib,':',trim(baxlab),':'
         if (igs(ib).gt.0) then 
            label = trim(label) // delim // trim(blabs(ib))
c            write(0,*) 'Label2:',trim(label),':',trim(blabs(ib)),':'
            units = trim(units) // delim // trim(baxlab)
         endif 
      enddo

c      write(0,*) 'Label3:',trim(label),':'
      write(signal_unit,'(a)') trim(label)
      write(signal_unit,'(a)') trim(units)

      tag = 'DESC:'
      write(signal_unit,'(a,1x,a,1x,50('' / '',a))') 
     >     trim(tag),trim(title),
     >     trim(view),trim(plane),trim(anly),trim(ref)

      tag = 'DATA:'
      write(signal_unit,'(a,1x,2(1x,i10))') 
     >      trim(tag), ia2-ia1 + 1, icnt+1

c      delim = ' '
      do ia = ia1,ia2
         write(val,'(1x,g18.8,1x)') as(ia) 
         line = trim(val)
c         write(0,*) 'L1:',trim(line),':',trim(val),':'
         do ib = 1,nbs
            if (igs(ib).gt.0) then 
               write(val,'(1x,g18.8,1x)') cs(ia,ib) 
               line = trim(line)//delim//trim(val)
c               write(0,*) 'L2:',trim(line),':',trim(val),':'
            endif
         end do

c         write(0,*) 'L3:',trim(line),':'

         write(signal_unit,'(a)') trim(line)
      end do 

      return
      end


