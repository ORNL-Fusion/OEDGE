      subroutine init_modules(nizs)
      use mtc
      use velocity_dist
      use subgrid
      implicit none
      integer :: nizs
      include 'params'
      include 'comtor'
      include 'hc_global_opts'

c
c     This subroutine is used to initialize data in modules which is used 
c     locally but where it is not easy or possible or desirable to include
c     common blocks with a lot of extraneous data that is not needed in 
c     the specific module. 
c
c     local variables
c
      integer :: max_hc_state
c
c     There were several examples where common blocks with a large number
c     of variables were being included in modules in order to access one
c     or two quantities. Until the code it reorganized it seems better
c     to just duplicate the limited number of variables that fall into this
c     category and initialize them at a global level. 
c
      call load_mtc_options(mtcopt,kelighi,kelighg)

c
c     Initialize the options to the debug neutral velocity routines
c
c     This code must be called here
c
c      call init_velocity_dist(debug_neutv,debug_neutv_einmax,
c     >                        debug_neutv_nbins)
c      


c
c     Initialize the subgrid module
c
c     global_hc_option - turns on and off hc subgrid
c     nizs = max ionization state on the regular grid
c     4 = max HC state number to be recorded on the hc subgrid
c
      max_hc_state = 4
      write(0,*) 'Init modules:',global_hc_follow_option,
     >              nizs,max_hc_state

      call init_subgrid(global_hc_follow_option,nizs,max_hc_state)


      return
      end
