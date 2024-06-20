subroutine styx_mpi_reduce_eirene_variables
  use eirmod_precision
  use eirmod_parmmod
  use eirmod_comusr
  use eirmod_cpes
  use eirmod_comsou
  use styx2eirene
  implicit none

  integer :: icolor,icomgrp,my_pe_gr,mxdim,istr,j,ier
  real(dp), allocatable :: help(:),alpha_V(:),beta_V(:)
  integer, allocatable :: helpi(:),hit_V(:)

  include 'mpif.h'

  if (sc_level > 1) then

     icolor=0

     if (nsteff < nprs) then
        ! collect from group leader pe
        icolor=MPI_UNDEFINED
        do istr=1,nstrata
           if (npts(istr) == 0) cycle
           if (my_pe == npesta(istr)) then
              icolor = 1
              exit
           end if
        end do
     else
        ! collect from all pes
        icolor = 1
     end if

     call mpi_barrier(mpi_comm_world,ier)

     call mpi_comm_split (mpi_comm_world,icolor,my_pe,icomgrp,ier)

     if (icolor == 1) then

        mxdim = max(NMOLI,NATMI,NIONI,NPHOTI,NPLSI)*Neir_cells

        allocate (help(mxdim))

        call mpi_comm_rank(mpi_comm_world,my_pe_gr,ier)

        do istr=1,NSTRATA
           if (npts(istr) == 0) cycle

           call mpi_barrier(icomgrp,ier)
           call mpi_reduce(eov(istr)%pdena,help,natmi*Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
           if (my_pe == 0) eov(istr)%pdena=reshape(help(1:natmi*Neir_cells),(/natmi,Neir_cells/))

           call mpi_barrier(icomgrp,ier)
           call mpi_reduce(eov(istr)%pdenm,help,nmoli*Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
           if (my_pe == 0) eov(istr)%pdenm=reshape(help(1:nmoli*Neir_cells),(/nmoli,Neir_cells/))

           call mpi_barrier(icomgrp,ier)
           call mpi_reduce(eov(istr)%pdeni,help,nioni*Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
           if (my_pe == 0) eov(istr)%pdeni=reshape(help(1:nioni*Neir_cells),(/nioni,Neir_cells/))

           call mpi_barrier(icomgrp,ier)
           call mpi_reduce(eov(istr)%edena,help,natmi*Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
           if (my_pe == 0) eov(istr)%edena=reshape(help(1:natmi*Neir_cells),(/natmi,Neir_cells/))

           call mpi_barrier(icomgrp,ier)
           call mpi_reduce(eov(istr)%edenm,help,nmoli*Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
           if (my_pe == 0) eov(istr)%edenm=reshape(help(1:nmoli*Neir_cells),(/nmoli,Neir_cells/))

           call mpi_barrier(icomgrp,ier)
           call mpi_reduce(eov(istr)%edeni,help,nioni*Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
           if (my_pe == 0) eov(istr)%edeni=reshape(help(1:nioni*Neir_cells),(/nioni,Neir_cells/))

           call mpi_barrier(icomgrp,ier)
           call mpi_reduce(eov(istr)%vxdena,help,natmi*Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
           if (my_pe == 0) eov(istr)%vxdena=reshape(help(1:natmi*Neir_cells),(/natmi,Neir_cells/))

           call mpi_barrier(icomgrp,ier)
           call mpi_reduce(eov(istr)%vxdenm,help,nmoli*Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
           if (my_pe == 0) eov(istr)%vxdenm=reshape(help(1:nmoli*Neir_cells),(/nmoli,Neir_cells/))

           call mpi_barrier(icomgrp,ier)
           call mpi_reduce(eov(istr)%vxdeni,help,nioni*Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
           if (my_pe == 0) eov(istr)%vxdeni=reshape(help(1:nioni*Neir_cells),(/nioni,Neir_cells/))

           call mpi_barrier(icomgrp,ier)
           call mpi_reduce(eov(istr)%vydena,help,natmi*Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
           if (my_pe == 0) eov(istr)%vydena=reshape(help(1:natmi*Neir_cells),(/natmi,Neir_cells/))

           call mpi_barrier(icomgrp,ier)
           call mpi_reduce(eov(istr)%vydenm,help,nmoli*Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
           if (my_pe == 0) eov(istr)%vydenm=reshape(help(1:nmoli*Neir_cells),(/nmoli,Neir_cells/))

           call mpi_barrier(icomgrp,ier)
           call mpi_reduce(eov(istr)%vydeni,help,nioni*Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
           if (my_pe == 0) eov(istr)%vydeni=reshape(help(1:nioni*Neir_cells),(/nioni,Neir_cells/))

           call mpi_barrier(icomgrp,ier)
           call mpi_reduce(eov(istr)%vzdena,help,natmi*Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
           if (my_pe == 0) eov(istr)%vzdena=reshape(help(1:natmi*Neir_cells),(/natmi,Neir_cells/))

           call mpi_barrier(icomgrp,ier)
           call mpi_reduce(eov(istr)%vzdenm,help,nmoli*Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
           if (my_pe == 0) eov(istr)%vzdenm=reshape(help(1:nmoli*Neir_cells),(/nmoli,Neir_cells/))

           call mpi_barrier(icomgrp,ier)
           call mpi_reduce(eov(istr)%vzdeni,help,nioni*Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
           if (my_pe == 0) eov(istr)%vzdeni=reshape(help(1:nioni*Neir_cells),(/nioni,Neir_cells/))

           if (direct_coupling) then

              call mpi_barrier(icomgrp,ier)
              call mpi_reduce(eov(istr)%pael,help,Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
              if (my_pe == 0) eov(istr)%pael=help(1:Neir_cells)

              call mpi_barrier(icomgrp,ier)
              call mpi_reduce(eov(istr)%pmel,help,Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
              if (my_pe == 0) eov(istr)%pmel=help(1:Neir_cells)

              call mpi_barrier(icomgrp,ier)
              call mpi_reduce(eov(istr)%piel,help,Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
              if (my_pe == 0) eov(istr)%piel=help(1:Neir_cells)

              call mpi_barrier(icomgrp,ier)
              call mpi_reduce(eov(istr)%papl,help,nplsi*Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
              if (my_pe == 0) eov(istr)%papl=reshape(help(1:nplsi*Neir_cells),(/nplsi,Neir_cells/))

              call mpi_barrier(icomgrp,ier)
              call mpi_reduce(eov(istr)%pmpl,help,nplsi*Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
              if (my_pe == 0) eov(istr)%pmpl=reshape(help(1:nplsi*Neir_cells),(/nplsi,Neir_cells/))

              call mpi_barrier(icomgrp,ier)
              call mpi_reduce(eov(istr)%pipl,help,nplsi*Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
              if (my_pe == 0) eov(istr)%pipl=reshape(help(1:nplsi*Neir_cells),(/nplsi,Neir_cells/))

              call mpi_barrier(icomgrp,ier)
              call mpi_reduce(eov(istr)%eael,help,Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
              if (my_pe == 0) eov(istr)%eael=help(1:Neir_cells)

              call mpi_barrier(icomgrp,ier)
              call mpi_reduce(eov(istr)%emel,help,Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
              if (my_pe == 0) eov(istr)%emel=help(1:Neir_cells)

              call mpi_barrier(icomgrp,ier)
              call mpi_reduce(eov(istr)%eiel,help,Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
              if (my_pe == 0) eov(istr)%eiel=help(1:Neir_cells)

              call mpi_barrier(icomgrp,ier)
              call mpi_reduce(eov(istr)%eapl,help,Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
              if (my_pe == 0) eov(istr)%eapl=help(1:Neir_cells)

              call mpi_barrier(icomgrp,ier)
              call mpi_reduce(eov(istr)%empl,help,Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
              if (my_pe == 0) eov(istr)%empl=help(1:Neir_cells)

              call mpi_barrier(icomgrp,ier)
              call mpi_reduce(eov(istr)%eipl,help,Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
              if (my_pe == 0) eov(istr)%eipl=help(1:Neir_cells)

              call mpi_barrier(icomgrp,ier)
              call mpi_reduce(eov(istr)%eaplr,help,nplsi*Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
              if (my_pe == 0) eov(istr)%eaplr=reshape(help(1:nplsi*Neir_cells),(/nplsi,Neir_cells/))

              call mpi_barrier(icomgrp,ier)
              call mpi_reduce(eov(istr)%emplr,help,nplsi*Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
              if (my_pe == 0) eov(istr)%emplr=reshape(help(1:nplsi*Neir_cells),(/nplsi,Neir_cells/))

              call mpi_barrier(icomgrp,ier)
              call mpi_reduce(eov(istr)%eiplr,help,nplsi*Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
              if (my_pe == 0) eov(istr)%eiplr=reshape(help(1:nplsi*Neir_cells),(/nplsi,Neir_cells/))

              call mpi_barrier(icomgrp,ier)
              call mpi_reduce(eov(istr)%mapl,help,natmi*Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
              if (my_pe == 0) eov(istr)%mapl=reshape(help(1:natmi*Neir_cells),(/natmi,Neir_cells/))

              call mpi_barrier(icomgrp,ier)
              call mpi_reduce(eov(istr)%mmpl,help,nmoli*Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
              if (my_pe == 0) eov(istr)%mmpl=reshape(help(1:nmoli*Neir_cells),(/nmoli,Neir_cells/))

              call mpi_barrier(icomgrp,ier)
              call mpi_reduce(eov(istr)%mipl,help,nioni*Neir_cells,mpi_double_precision,mpi_sum,0,icomgrp,ier)
              if (my_pe == 0) eov(istr)%mipl=reshape(help(1:nioni*Neir_cells),(/nioni,Neir_cells/))

           endif

        enddo

        deallocate(help)

        call mpi_barrier(icomgrp,ier)
        call mpi_comm_free(icomgrp,ier)

     endif

  endif

  
  call mpi_barrier(mpi_comm_world,ier)

!!!! reducing incidence angle data !!!!!!
!!!! this should only involve processes working on recycling strata !!!
!!!! but ok since arrays are initialized to 0 for other strata      !!!

  allocate(alpha_V(Nsou),beta_V(Nsou),hit_V(Nsou))
  do j=1,Nsou
    alpha_V(j) = sheath1D(j)%alpha_V
    beta_V(j)  = sheath1D(j)%beta_V
    hit_V(j)   = sheath1D(j)%hit_V
  enddo

  allocate(helpi(nsou))

  call mpi_reduce(hit_V,helpi,nsou,mpi_integer,mpi_sum,0,mpi_comm_world,ier)
  if (my_pe == 0) hit_V=helpi

  deallocate(helpi)
  allocate(help(nsou))

  call mpi_reduce(alpha_V,help,nsou,mpi_double_precision,mpi_sum,0,mpi_comm_world,ier)
  if (my_pe == 0) alpha_V=help
  
  call mpi_reduce(beta_V,help,nsou,mpi_double_precision,mpi_sum,0,mpi_comm_world,ier)
  if (my_pe == 0) beta_V=help

  do j=1,Nsou
    sheath1D(j)%alpha_V= alpha_V(j)
    sheath1D(j)%beta_V = beta_V(j)
    sheath1D(j)%hit_V  = hit_V(j)
  enddo

  deallocate(help,alpha_V,beta_V,hit_V)

  




end subroutine styx_mpi_reduce_eirene_variables

