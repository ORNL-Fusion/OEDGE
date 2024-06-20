module Mlist

  type cell
     integer*4 row
     integer*4 col
     real*8 val
     type(cell), pointer :: p => null()
  end type cell

contains

  recursive subroutine list_all(ptr)
    implicit none
    type(cell), pointer :: ptr
    if(associated(ptr%p)) call list_all(ptr%p)
    print *, ptr%row, ptr%col, ptr%val
  end subroutine list_all

  recursive subroutine free_all(ptr)
    implicit none
    type(cell), pointer :: ptr
    if(associated(ptr%p)) call free_all(ptr%p)
    deallocate(ptr)
  end subroutine free_all

  recursive subroutine list2array(rows,cols,vals,nnz,ind,ptr)
    implicit none
    type(cell), pointer :: ptr
    integer*4 nnz
    integer*4 rows(nnz)
    integer*4 cols(nnz)
    real*8 vals(nnz)
    integer*4 ind
    if(associated(ptr%p)) then
       call list2array(rows,cols,vals,nnz,ind+1,ptr%p)
    end if
    rows(ind)=ptr%row
    cols(ind)=ptr%col
    vals(ind)=ptr%val
  end subroutine list2array

end module Mlist
