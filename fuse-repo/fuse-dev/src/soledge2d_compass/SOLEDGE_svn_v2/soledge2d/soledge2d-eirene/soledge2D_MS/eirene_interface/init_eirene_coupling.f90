subroutine init_eirene_coupling()

  use all_variables
  use styx2eirene, only : Sn_intg,dt_eirene
  use eirmod_cpes
  implicit none
  integer  :: op,ot,omt,omp_get_num_procs,omp_get_num_threads,omp_get_max_threads  
  integer*4 ier

  include 'mpif.h'

  omt=omp_get_max_threads()
  if(my_pe==0) then
     call omp_set_num_threads(omt)
  else
     call omp_set_num_threads(1)
  end if
  op=omp_get_num_procs()
  ot=omp_get_num_threads()
  omt=omp_get_max_threads()
  write(*,*) 'I am process ',my_pe,'   :  I use ', op, ' cpus'
  write(*,*) 'I am process ',my_pe,'   :  I use now ', ot, ' threads'
  write(*,*) 'I am process ',my_pe,'   :  I may use max ', omt, ' threads'
  
  if (my_pe == 0)  then
     call styx_load_param_eirene()
     call styx_init_eirene(1)
  endif
  call eirene_init_phase(dt_eirene,.false.,.true.,1,.true.)
  call styx_init_eirene(2)
  if (my_pe == 0) then              
     call check_volume_consistency()
     write(*,*) 'Eirene : OK'
  endif


end subroutine init_eirene_coupling
