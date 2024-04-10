      program prinequ
c
c  version : 02.09.99 14:21
c
c=====================================================
c*** Printout of dg-compatible equilibrium data in a table
c***
c*** The file is read from standard input and the result sent to the
c*** standard output
c=====================================================
      parameter (ngpr=129, ngpz=257)
      real*8 fg(ngpr),pg(ngpr),ffg(ngpr),ppg(ngpr)
      real*8 pfm(ngpr,ngpz),rgr(ngpr),zgr(ngpz)
      real*8 rdim,zdim,rcntc,redge,zmsmid,rma,zma,psimin,psilim,btorc
      character title*40, date*9, page
c=====================================================
c
      page=char(12)
c      call date2(date)
      call date_and_time(date)
      call rdeqdg(1,ngpr,ngpz,iret, nr,nz,btorc,rcntc,rgr,zgr,pfm)
      if(iret.ne.0) then
          print *,'==== dg2ef: error in rdeqdg. iret =',iret
          stop
      end if
c
      nw=30
      ne=0
      nb=1
      do k=1,nr
          ne=min0(ne+nw,nr)
          print '(6x,30i6)',(i,i=nb,ne)
          print '(8x,30f6.3)',(rgr(i),i=nb,ne)
          print *
          do j=1,nz
              print '(2x,31f6.3)',zgr(j),(pfm(i,j),i=nb,ne)
          end do
          print *,page
          nb=ne+1
          if(nb.gt.nr) stop
      end do
c
      end
