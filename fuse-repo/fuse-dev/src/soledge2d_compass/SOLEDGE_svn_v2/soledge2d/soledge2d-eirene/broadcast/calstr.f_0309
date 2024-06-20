!pb 300806  reduce commands for WTOTE and EELFI added
!           reduction of spectrum data corrected
!pb 181206  group management changed
!pb 110707  calls to mpi_reduce corrected 

      SUBROUTINE EIRENE_CALSTR

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CSPEZ
      USE EIRMOD_COMPRT
      USE EIRMOD_CPES
      USE EIRMOD_CSDVI
      USE EIRMOD_COUTAU
      IMPLICIT NONE

C
      INCLUDE 'mpif.h'
      real(dp), allocatable :: help(:)
      integer :: igrp(0:nstra), icomgrp(0:nstra)
      integer :: ier1, ier, ir, npean, npeen, i, mpicw, ispc, my_pe_gr,
     .           mxdim, ns
      logical, allocatable :: lhelp(:)

      npean = npesta(istra)
      npeen = npesta(istra)+npestr(istra)-1
      call mpi_comm_group (mpi_comm_world,mpicw,ier)
      call mpi_comm_split (mpi_comm_world,istra,my_pe-npesta(istra),
     .                     icomgrp(istra),ier)
c
c   collect data from pe's belonging to one stratum
c
      if (istra.eq.nstrpe(my_pe)) then
        call mpi_barrier(icomgrp(istra),ier)
        my_pe_gr = my_pe-npesta(istra)
        
        mxdim = max(nvoltl,nsrftl,nsd,nsdw,
     .              nmoli+1,natmi+1,nioni+1,nphoti+1,nplsi+1)
     	allocate (help(mxdim))
     	
        call mpi_reduce(WTOTM(0,istra),help,nmoli+1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) WTOTM(0:nmoli,istra) = help(1:nmoli+1)

        call mpi_reduce(WTOTA(0,istra),help,natmi+1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) WTOTA(0:natmi,istra) = help(1:natmi+1)

        call mpi_reduce(WTOTI(0,istra),help,nioni+1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) WTOTI(0:nioni,istra) = help(1:nioni+1)

        call mpi_reduce(WTOTP(0,istra),help,nplsi+1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) WTOTP(0:nplsi,istra) = help(1:nplsi+1)

        call mpi_reduce(WTOTE(istra),help,1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) WTOTE(istra) = help(1)
C
C
        call mpi_reduce(XMCP(istra),help,1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) XMCP(istra) = help(1)

        call mpi_reduce(PTRASH(istra),help,1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) PTRASH(istra) = help(1)

        call mpi_reduce(ETRASH(istra),help,1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) ETRASH(istra) = help(1)

        call mpi_reduce(ETOTA(istra),help,1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) ETOTA(istra) = help(1)

        call mpi_reduce(ETOTM(istra),help,1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) ETOTM(istra) = help(1)

        call mpi_reduce(ETOTI(istra),help,1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) ETOTI(istra) = help(1)

        call mpi_reduce(ETOTP(istra),help,1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) ETOTP(istra) = help(1)

        call mpi_reduce(EELFI(0,istra),help,nioni+1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) EELFI(0:nioni,istra) = help(1:nioni+1)


        call mpi_reduce(PPMLI(0,istra),help,nmoli+1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) PPMLI(0:nmoli,istra) = help(1:nmoli+1)

        call mpi_reduce(PPATI(0,istra),help,natmi+1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) PPATI(0:natmi,istra) = help(1:natmi+1)

        call mpi_reduce(PPIOI(0,istra),help,nioni+1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) PPIOI(0:nioni,istra) = help(1:nioni+1)

        call mpi_reduce(PPPHTI(0,istra),help,nphoti+1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) PPPHTI(0:nphoti,istra) = help(1:nphoti+1)

        call mpi_reduce(PPPLI(0,istra),help,nplsi+1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) PPPLI(0:nplsi,istra) = help(1:nplsi+1)


        call mpi_reduce(EPATI(istra),help,1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) EPATI(istra) = help(1)

        call mpi_reduce(EPMLI(istra),help,1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) EPMLI(istra) = help(1)

        call mpi_reduce(EPIOI(istra),help,1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) EPIOI(istra) = help(1)

        call mpi_reduce(EPPHTI(istra),help,1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) EPPHTI(istra) = help(1)

        call mpi_reduce(EPPLI(istra),help,1,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) EPPLI(istra) = help(1)

C
C
        do ir=1,nrtal
          call mpi_reduce(estimv(1,ir),help,nvoltl,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	  if (my_pe_gr==0) estimv(1:nvoltl,ir) = help(1:nvoltl)
        end do

	do ir=1,nlmpgs
          call mpi_reduce(estims(1,ir),help,nsrftl,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	  if (my_pe_gr==0) estims(1:nsrftl,ir) = help(1:nsrftl)
        end do

        do ispc=1,nadspc
          ns = estiml(ispc)%pspc%nspc
          call mpi_reduce(estiml(ispc)%pspc%spc,help,
     .                    estiml(ispc)%pspc%nspc+2,
     .                    mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
          if (my_pe_gr==0) estiml(ispc)%pspc%spc(0:ns+1) = help(1:ns+2)
          
          if (nsigi_spc > 0) then
            call mpi_reduce(estiml(ispc)%pspc%sdv,help,
     .                      estiml(ispc)%pspc%nspc+2,
     .                      mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
            if (my_pe_gr==0) 
     .        estiml(ispc)%pspc%sdv(0:ns+1) = help(1:ns+2)
            
            call mpi_reduce(estiml(ispc)%pspc%sgm,help,
     .                      estiml(ispc)%pspc%nspc+2,
     .                      mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
            if (my_pe_gr==0) 
     .        estiml(ispc)%pspc%sgm(0:ns+1) = help(1:ns+2)
            
            call mpi_reduce(estiml(ispc)%pspc%sgms,help,
     .                      1,mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
            if (my_pe_gr==0) estiml(ispc)%pspc%sgms = help(1)
          end if
        end do
C
	if (nsd > 0) then
          do ir=1,nrtal+1
            call mpi_reduce(sdvi1(1,ir),help,nsd,
     .                      mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) sdvi1(1:nsd,ir) = help(1:nsd)
          end do
	end if     

	if (nsdw > 0) then
          do ir=1,nlimps+1
            call mpi_reduce(sdvi2(1,ir),help,nsdw,
     .                      mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) sdvi2(1:nsdw,ir) = help(1:nsdw)
          end do
	end if

	if (ncv > 0) then
          do ir=1,nrtal
     	    call mpi_reduce(sigmac(0,1,ir),help,3*ncv,
     .                  mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) sigmac(0:2,1:ncv,ir) = 
     .                       reshape(help(1:3*ncv), (/ 3, ncv /))
          end do

          call mpi_reduce(sgmcs,help,nsdvc2,
     .                    mpi_real8,mpi_sum,0,icomgrp(istra),ier1)
	  if (my_pe_gr==0) sgmcs(0:2,1:ncv) = 
     .                     reshape(help(1:3*ncv), (/ 3, ncv /))
	end if

	deallocate(help)
	mxdim = max (nmoli+1,natmi+1,nioni+1,nphoti+1,nplsi+1)
	allocate (lhelp(mxdim))

        call mpi_reduce(LOGMOL(0,ISTRA),lhelp,NMOLI+1,
     .                  mpi_logical,mpi_LOR,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) LOGMOL(0:nmoli,ISTRA) = lhelp(1:nmoli+1)

        call mpi_reduce(LOGATM(0,ISTRA),lhelp,NATMI+1,
     .                  mpi_logical,mpi_LOR,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) LOGATM(0:natmi,ISTRA) = lhelp(1:natmi+1)

        call mpi_reduce(LOGION(0,ISTRA),lhelp,NIONI+1,
     .                  mpi_logical,mpi_LOR,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) LOGION(0:nioni,ISTRA) = lhelp(1:nioni+1)

        call mpi_reduce(LOGPHOT(0,ISTRA),lhelp,NPHOTI+1,
     .                  mpi_logical,mpi_LOR,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) LOGPHOT(0:nphoti,ISTRA) = lhelp(1:nphoti+1)

        call mpi_reduce(LOGPLS(0,ISTRA),lhelp,NPLSI+1,
     .                  mpi_logical,mpi_LOR,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) LOGPLS(0:nplsi,ISTRA) = lhelp(1:nplsi+1)
	
	deallocate(lhelp)

        call mpi_barrier(icomgrp(istra),ier)
        call mpi_comm_free (icomgrp(istra),ier)

      endif

      call mpi_barrier(mpi_comm_world,ier)
      
      call eirene_calstr_usr

      call mpi_barrier(mpi_comm_world,ier)

      RETURN
      END
