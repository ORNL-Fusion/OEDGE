
      subroutine  Read_AdditionalPlotData(buffer)
      use mod_params
      use mod_outcom
      implicit none
      character*(*) buffer
c     include 'params'
c     include 'outcom'
c
c     Read_AdditionalPlotData:
c
c     This routine reads in addtional optional plot data that
c     is defined in the input file by a "#" followed by a two
c     character designation. Each option is described below - an
c     unexpected option will generate a message but not an error. 
c
      integer tmp_nsets,in
c      integer lenstr
c      external lenstr
      character*80 graph
c
c
c     #01 - read in contour plot scaling data
c
      if (buffer(2:4).eq.'#01') then   

         READ (BUFFER,*,ERR=9999,END=9999) graph,
     >                      minscale,maxscale,localcngs
c
c     #02 - read in contour plot zoom data
c
      elseif (buffer(2:4).eq.'#02') then   

         READ(buffer,*,err=9999,end=9999) graph,
     >                        xcen,ycen,xnear2,ynear2

c
c     #03 - reads in a list of experimental datasets to be included
c           on the plot if possible. 
c

      elseif (buffer(2:4).eq.'#03') then   

         READ(buffer,*,err=9999,end=9999) graph,
     >                        tmp_nsets
c
c        Limit number of datasets to max specified in code
c 
         expt_nsets= min(tmp_nsets,max_expt_datasets)
c
         if (expt_nsets.gt.0) then 

            READ(buffer,*,err=9999,end=9999) graph,
     >          tmp_nsets,(expt_datasets(in),in=1,expt_nsets)
       
         else
            expt_nsets = 0
         endif
c
c     Issue message for unrecognized option 
c
      else 
c

         len = lenstr(buffer) 
         write(6,'(a,a)') 'Unrecognized Optional Plot Input Line:',
     >                    buffer(1:len)
         write(0,'(a,a)') 'Unrecognized Optional Plot Input Line:',
     >                    buffer(1:len)


      endif

      return 
c
c     Trap I/O errors and issue some message
c

 9999 len = lenstr(buffer)

      write(6,'(a,a)') 'ERROR reading additional plot data:',
     >                 buffer(1:len)
      write(0,'(a,a)') 'ERROR reading additional plot data:',
     >                 buffer(1:len)
     
      return 
      end 

