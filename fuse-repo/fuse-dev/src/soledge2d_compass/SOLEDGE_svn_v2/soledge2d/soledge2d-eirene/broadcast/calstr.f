!pb 300806  reduce commands for WTOTE and EELFI added
!           reduction of spectrum data corrected
!pb 181206  group management changed
!pb 110707  calls to mpi_reduce corrected 
!pb 060309  mpi_real8 --> mpi_double_precision
!pb 090309  rewritten to use automatic arrays as output buffer in mpi_reduce
!pb 090309  loops reorganized          
!pb 270309  typos corrected             

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
      real(dp), allocatable :: help(:), helpest(:)
      real(dp) :: helpa(0:natm), helpm(0:nmol), helpi(0:nion),
     .            helpp(0:npls), helpph(0:nphot), helpv(nrtal+1),
     .            helps(nlmpgs+1), helpc
      real(dp) :: dummyv(nrtal+1), dummys(nlmpgs+1) 
      integer :: igrp(0:nstra), icomgrp(0:nstra)
      integer :: ier1, ier, ir, npean, npeen, i, mpicw, ispc, my_pe_gr,
     .           mxdim, ns, j
      logical, allocatable :: lhelp(:)
      logical :: lhelpa(0:natm),lhelpm(0:nmol), lhelpi(0:nion),
     .           lhelpp(0:npls), lhelpph(0:nphot)

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
     	
        call mpi_reduce(WTOTM(0,istra),helpm,nmoli+1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) WTOTM(0:nmoli,istra) = helpm(0:nmoli)

        call mpi_reduce(WTOTA(0,istra),helpa,natmi+1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) WTOTA(0:natmi,istra) = helpa(0:natmi)

        call mpi_reduce(WTOTI(0,istra),helpi,nioni+1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) WTOTI(0:nioni,istra) = helpi(00:nioni)

        call mpi_reduce(WTOTP(0,istra),helpp,nplsi+1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) WTOTP(0:nplsi,istra) = helpp(0:nplsi)

        call mpi_reduce(WTOTE(istra),helpc,1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) WTOTE(istra) = helpc
C
C
        call mpi_reduce(XMCP(istra),helpc,1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) XMCP(istra) = helpc

        call mpi_reduce(PTRASH(istra),helpc,1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) PTRASH(istra) = helpc

        call mpi_reduce(ETRASH(istra),helpc,1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) ETRASH(istra) = helpc

        call mpi_reduce(ETOTA(istra),helpc,1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) ETOTA(istra) = helpc

        call mpi_reduce(ETOTM(istra),helpc,1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) ETOTM(istra) = helpc

        call mpi_reduce(ETOTI(istra),helpc,1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) ETOTI(istra) = helpc

        call mpi_reduce(ETOTP(istra),helpc,1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) ETOTP(istra) = helpc

        call mpi_reduce(EELFI(0,istra),helpi,nioni+1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) EELFI(0:nioni,istra) = helpi(0:nioni)


        call mpi_reduce(PPMLI(0,istra),helpm,nmoli+1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) PPMLI(0:nmoli,istra) = helpm(0:nmoli)

        call mpi_reduce(PPATI(0,istra),helpa,natmi+1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) PPATI(0:natmi,istra) = helpa(0:natmi)

        call mpi_reduce(PPIOI(0,istra),helpi,nioni+1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) PPIOI(0:nioni,istra) = helpi(0:nioni)

        if (nphoti > 0) then
        call mpi_reduce(PPPHTI(0,istra),helpph,nphoti+1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) PPPHTI(0:nphoti,istra) = helpph(0:nphoti)
        end if

        call mpi_reduce(PPPLI(0,istra),helpp,nplsi+1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) PPPLI(0:nplsi,istra) = helpp(0:nplsi)


        call mpi_reduce(EPATI(istra),helpc,1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) EPATI(istra) = helpc

        call mpi_reduce(EPMLI(istra),helpc,1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) EPMLI(istra) = helpc

        call mpi_reduce(EPIOI(istra),helpc,1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) EPIOI(istra) = helpc

        call mpi_reduce(EPPHTI(istra),helpc,1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) EPPHTI(istra) = helpc

        call mpi_reduce(EPPLI(istra),helpc,1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) EPPLI(istra) = helpc

C
C
        do ir=1,nvoltl
          dummyv(1:nrtal) = estimv(ir,1:nrtal)
          call mpi_reduce(dummyv,helpv,nrtal,
     .         mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	  if (my_pe_gr==0) estimv(ir,1:nrtal) = helpv(1:nrtal)
        end do

	do ir=1,nsrftl
          dummys(1:nlmpgs) = estims(ir,1:nlmpgs)
          call mpi_reduce(dummys,helps,nlmpgs,
     .         mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	  if (my_pe_gr==0) estims(ir,1:nlmpgs) = helps(1:nlmpgs)
        end do

        do ispc=1,nadspc
          ns = estiml(ispc)%pspc%nspc
          allocate (helpest(ns+2))
          call mpi_reduce(estiml(ispc)%pspc%spc,helpest,
     .                    estiml(ispc)%pspc%nspc+2,
     .         mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
          if (my_pe_gr==0) estiml(ispc)%pspc%spc(0:ns+1)=helpest(1:ns+2)
          
          if (nsigi_spc > 0) then
            call mpi_reduce(estiml(ispc)%pspc%sdv,helpest,
     .                      estiml(ispc)%pspc%nspc+2,
     .           mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
            if (my_pe_gr==0) 
     .        estiml(ispc)%pspc%sdv(0:ns+1) = helpest(1:ns+2)
            
            call mpi_reduce(estiml(ispc)%pspc%sgm,helpest,
     .                      estiml(ispc)%pspc%nspc+2,
     .           mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
            if (my_pe_gr==0) 
     .        estiml(ispc)%pspc%sgm(0:ns+1) = helpest(1:ns+2)
            
            call mpi_reduce(estiml(ispc)%pspc%sgms,helpest,1,
     .           mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
            if (my_pe_gr==0) estiml(ispc)%pspc%sgms = helpest(1)
          end if
          deallocate (helpest)
        end do
C
	if (nsd > 0) then
          do ir=1,nsd
            dummyv(1:nrtal+1) = sdvi1(ir,1:nrtal+1)
            call mpi_reduce(dummyv,helpv,nrtal+1,
     .           mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) sdvi1(ir,1:nrtal+1) = helpv(1:nrtal+1)
          end do
	end if     

	if (nsdw > 0) then
          do ir=1,nsdw
            dummys(1:nlimps+1) = sdvi2(ir,1:nlimps+1)
            call mpi_reduce(dummys,helps,nlimps+1,
     .           mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) sdvi2(ir,1:nlimps+1) = helps(1:nlimps+1)
          end do
	end if

	if (ncv > 0) then
          do i=0,2
            do j=1,ncv
              dummyv(1:nrtal) = sigmac(i,j,1:nrtal)
     	      call mpi_reduce(dummyv,helpv,nrtal,
     .           mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	      if (my_pe_gr==0) sigmac(i,j,1:nrtal) = helpv(1:nrtal)
            end do
          end do

          call mpi_reduce(sgmcs(0,1:ncv),helpv,ncv,
     .         mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	  if (my_pe_gr==0) sgmcs(0,1:ncv) = helpv(1:ncv)

          call mpi_reduce(sgmcs(1,1:ncv),helpv,ncv,
     .         mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	  if (my_pe_gr==0) sgmcs(1,1:ncv) = helpv(1:ncv)

          call mpi_reduce(sgmcs(2,1:ncv),helpv,ncv,
     .         mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	  if (my_pe_gr==0) sgmcs(2,1:ncv) = helpv(1:ncv)
	end if

	deallocate(help)
	mxdim = max (nmoli+1,natmi+1,nioni+1,nphoti+1,nplsi+1)
	allocate (lhelp(mxdim))

        call mpi_reduce(LOGMOL(0,ISTRA),lhelpm,NMOLI+1,
     .       mpi_logical,mpi_LOR,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) LOGMOL(0:nmoli,ISTRA) = lhelpm(0:nmoli)

        call mpi_reduce(LOGATM(0,ISTRA),lhelpa,NATMI+1,
     .       mpi_logical,mpi_LOR,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) LOGATM(0:natmi,ISTRA) = lhelpa(0:natmi)

        call mpi_reduce(LOGION(0,ISTRA),lhelpi,NIONI+1,
     .       mpi_logical,mpi_LOR,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) LOGION(0:nioni,ISTRA) = lhelpi(0:nioni)

        if (nphoti > 0) then
        call mpi_reduce(LOGPHOT(0,ISTRA),lhelpph,NPHOTI+1,
     .       mpi_logical,mpi_LOR,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) LOGPHOT(0:nphoti,ISTRA) = lhelpph(0:nphoti)
        end if

        call mpi_reduce(LOGPLS(0,ISTRA),lhelpp,NPLSI+1,
     .       mpi_logical,mpi_LOR,0,icomgrp(istra),ier1)
	if (my_pe_gr==0) LOGPLS(0:nplsi,ISTRA) = lhelpp(0:nplsi)
	
	deallocate(lhelp)

        call mpi_barrier(icomgrp(istra),ier)
        call mpi_comm_free (icomgrp(istra),ier)

      endif

      call mpi_barrier(mpi_comm_world,ier)
      
      call eirene_calstr_usr

      call mpi_barrier(mpi_comm_world,ier)

      RETURN
      END
