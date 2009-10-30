module divimp_utilities
  use error_handling
  implicit none

  private

  public :: sort_arrays


contains


  subroutine sort_arrays(order_opt,n_elements,arr1,arr2,arr3,arr4,arr5,arr6)
    implicit none
    integer :: order_opt,n_elements
    real :: arr1(n_elements)
    real,optional :: arr2(n_elements),arr3(n_elements),arr4(n_elements),arr5(n_elements),arr6(n_elements)
    real :: temp_val
    integer :: in1,in2

    ! Sort the data in arr1 in either ascending (order_opt=0) or descending (order_opt=1) order
    ! Data in arrays 2 to 6 will be sorted to match arr1



    do in1 = 1,n_elements
       do in2 = in1+1, n_elements
          if (order_opt.eq.0) then 
             if (arr1(in2).lt.arr1(in1)) then 
                ! swap
                call swap(arr1,n_elements,in1,in2)
                if (present(arr2)) call swap(arr2,n_elements,in1,in2)
                if (present(arr3)) call swap(arr3,n_elements,in1,in2)
                if (present(arr4)) call swap(arr4,n_elements,in1,in2)
                if (present(arr5)) call swap(arr5,n_elements,in1,in2)
                if (present(arr6)) call swap(arr6,n_elements,in1,in2)

             endif
          elseif (order_opt.eq.1) then 
             if (arr1(in2).gt.arr1(in1)) then 
                ! swap
                call swap(arr1,n_elements,in1,in2)
                if (present(arr2)) call swap(arr2,n_elements,in1,in2)
                if (present(arr3)) call swap(arr3,n_elements,in1,in2)
                if (present(arr4)) call swap(arr4,n_elements,in1,in2)
                if (present(arr5)) call swap(arr5,n_elements,in1,in2)
                if (present(arr6)) call swap(arr6,n_elements,in1,in2)

             endif
          endif

       end do
    end do

  end subroutine sort_arrays

  subroutine swap(arr,nele,in1,in2)
    implicit none
    integer :: nele, in1, in2
    real :: arr(nele)
    real :: temp

    temp = arr(in2)
    arr(in2) = arr(in1)
    arr(in1) = temp

  end subroutine swap

end module divimp_utilities
