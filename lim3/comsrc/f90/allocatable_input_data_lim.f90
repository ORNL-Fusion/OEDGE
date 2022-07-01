module allocatable_input_data
      implicit none

      public
      
      ! This module contains the declaration for all arrays that could be allocated prior to mass allocation
      ! This includes most unstructured input arrays
      

      real,public,allocatable:: extfluxdata(:,:)
      real,public,allocatable:: ss_cymfs(:,:),mnbg(:,:)
        
      real, allocatable,public:: pbin_bnds(:),surf_bnds(:,:)

      real,public,allocatable :: absorb_surf_data(:,:)
      real,public,allocatable :: absorb_plasma(:,:)
 
      
      contains


        subroutine deallocate_allocatable_input
          implicit none

          ! from mod_comxyt
          if (allocated(pbin_bnds)) deallocate(pbin_bnds)
          if (allocated(surf_bnds)) deallocate(surf_bnds)

          
          ! from mod_comtor
          if (allocated(ss_cymfs)) deallocate(ss_cymfs)
          if (allocated(mnbg)) deallocate(mnbg)
          if (allocated(extfluxdata)) deallocate(extfluxdata)

          ! from yreflection
          if (allocated(absorb_surf_data)) deallocate(absorb_surf_data)
          if (allocated(absorb_plasma)) deallocate(absorb_plasma)

        end subroutine deallocate_allocatable_input


        
  end module allocatable_input_data
