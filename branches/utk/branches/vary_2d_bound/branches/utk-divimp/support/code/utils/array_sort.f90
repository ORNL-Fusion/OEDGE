module array_sort
  use error_handling
  implicit none

  interface sort_arrays

     module procedure r_sort_arrays,r8_sort_arrays,i_sort_arrays

  end interface


  private

  public :: sort_arrays


contains


  subroutine r_sort_arrays(order_opt,n_elements,arr1,arr2,arr3,arr4,arr5,arr6)
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
                ! r_swap
                call r_swap(arr1,n_elements,in1,in2)
                if (present(arr2)) call r_swap(arr2,n_elements,in1,in2)
                if (present(arr3)) call r_swap(arr3,n_elements,in1,in2)
                if (present(arr4)) call r_swap(arr4,n_elements,in1,in2)
                if (present(arr5)) call r_swap(arr5,n_elements,in1,in2)
                if (present(arr6)) call r_swap(arr6,n_elements,in1,in2)

             endif
          elseif (order_opt.eq.1) then 
             if (arr1(in2).gt.arr1(in1)) then 
                ! r_swap
                call r_swap(arr1,n_elements,in1,in2)
                if (present(arr2)) call r_swap(arr2,n_elements,in1,in2)
                if (present(arr3)) call r_swap(arr3,n_elements,in1,in2)
                if (present(arr4)) call r_swap(arr4,n_elements,in1,in2)
                if (present(arr5)) call r_swap(arr5,n_elements,in1,in2)
                if (present(arr6)) call r_swap(arr6,n_elements,in1,in2)

             endif
          endif

       end do
    end do

  end subroutine r_sort_arrays


  subroutine r8_sort_arrays(order_opt,n_elements,arr1,arr2,arr3,arr4,arr5,arr6)
    implicit none
    integer :: order_opt,n_elements
    real*8 :: arr1(n_elements)
    real*8,optional :: arr2(n_elements),arr3(n_elements),arr4(n_elements),arr5(n_elements),arr6(n_elements)
    real*8 :: temp_val
    integer :: in1,in2

    ! Sort the data in arr1 in either ascending (order_opt=0) or descending (order_opt=1) order
    ! Data in arrays 2 to 6 will be sorted to match arr1



    do in1 = 1,n_elements
       do in2 = in1+1, n_elements
          if (order_opt.eq.0) then 
             if (arr1(in2).lt.arr1(in1)) then 
                ! r8_swap
                call r8_swap(arr1,n_elements,in1,in2)
                if (present(arr2)) call r8_swap(arr2,n_elements,in1,in2)
                if (present(arr3)) call r8_swap(arr3,n_elements,in1,in2)
                if (present(arr4)) call r8_swap(arr4,n_elements,in1,in2)
                if (present(arr5)) call r8_swap(arr5,n_elements,in1,in2)
                if (present(arr6)) call r8_swap(arr6,n_elements,in1,in2)

             endif
          elseif (order_opt.eq.1) then 
             if (arr1(in2).gt.arr1(in1)) then 
                ! r8_swap
                call r8_swap(arr1,n_elements,in1,in2)
                if (present(arr2)) call r8_swap(arr2,n_elements,in1,in2)
                if (present(arr3)) call r8_swap(arr3,n_elements,in1,in2)
                if (present(arr4)) call r8_swap(arr4,n_elements,in1,in2)
                if (present(arr5)) call r8_swap(arr5,n_elements,in1,in2)
                if (present(arr6)) call r8_swap(arr6,n_elements,in1,in2)

             endif
          endif

       end do
    end do

  end subroutine r8_sort_arrays


  subroutine i_sort_arrays(order_opt,n_elements,arr1,arr2,arr3,arr4,arr5,arr6)
    implicit none
    integer :: order_opt,n_elements
    integer :: arr1(n_elements)
    integer,optional :: arr2(n_elements),arr3(n_elements),arr4(n_elements),arr5(n_elements),arr6(n_elements)
    integer :: temp_val
    integer :: in1,in2

    ! Sort the data in arr1 in either ascending (order_opt=0) or descending (order_opt=1) order
    ! Data in arrays 2 to 6 will be sorted to match arr1



    do in1 = 1,n_elements
       do in2 = in1+1, n_elements
          if (order_opt.eq.0) then 
             if (arr1(in2).lt.arr1(in1)) then 
                ! i_swap
                call i_swap(arr1,n_elements,in1,in2)
                if (present(arr2)) call i_swap(arr2,n_elements,in1,in2)
                if (present(arr3)) call i_swap(arr3,n_elements,in1,in2)
                if (present(arr4)) call i_swap(arr4,n_elements,in1,in2)
                if (present(arr5)) call i_swap(arr5,n_elements,in1,in2)
                if (present(arr6)) call i_swap(arr6,n_elements,in1,in2)

             endif
          elseif (order_opt.eq.1) then 
             if (arr1(in2).gt.arr1(in1)) then 
                ! i_swap
                call i_swap(arr1,n_elements,in1,in2)
                if (present(arr2)) call i_swap(arr2,n_elements,in1,in2)
                if (present(arr3)) call i_swap(arr3,n_elements,in1,in2)
                if (present(arr4)) call i_swap(arr4,n_elements,in1,in2)
                if (present(arr5)) call i_swap(arr5,n_elements,in1,in2)
                if (present(arr6)) call i_swap(arr6,n_elements,in1,in2)

             endif
          endif

       end do
    end do

  end subroutine i_sort_arrays


  subroutine r_swap(arr,nele,in1,in2)
    implicit none
    integer :: nele, in1, in2
    real :: arr(nele)
    real :: temp

    temp = arr(in2)
    arr(in2) = arr(in1)
    arr(in1) = temp

  end subroutine r_swap

  subroutine r8_swap(arr,nele,in1,in2)
    implicit none
    integer :: nele, in1, in2
    real*8 :: arr(nele)
    real*8 :: temp

    temp = arr(in2)
    arr(in2) = arr(in1)
    arr(in1) = temp

  end subroutine r8_swap

  subroutine i_swap(arr,nele,in1,in2)
    implicit none
    integer :: nele, in1, in2
    integer :: arr(nele)
    integer :: temp

    temp = arr(in2)
    arr(in2) = arr(in1)
    arr(in1) = temp

  end subroutine i_swap

end module array_sort
