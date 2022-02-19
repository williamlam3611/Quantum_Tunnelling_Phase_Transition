module sort
implicit none

private

public :: sort_bubble, sort_mode, sort_normalise

interface sort_normalise
    module procedure sort_normalise_, sort_normalise_array
end interface sort_normalise

contains
    function sort_normalise_(array) result(ret)
        real*8, intent(in) :: array(:)
        real*8             :: ret(size(array))
        ret = array / sum(array)
    
    end function sort_normalise_
    
    function sort_normalise_array(array) result(ret)
        real*8, intent(in) :: array(:, :)
        real*8             :: ret(size(array(:, 1)), size(array(1, :)))
        integer            :: i
        do i = 1, size(array(:, 1))
            ret(i, :) = sort_normalise(array(i, :))
        end do
        
    end function sort_normalise_array

    subroutine sort_bubble(array)
        integer, intent(inout) :: array(:)
        integer                :: temp
        integer                :: i, j
        logical                :: swapped
        do j = size(array) - 1, 1, -1
            swapped = .FALSE.
            do i = 1, j
                if (array(i) > array(i + 1)) then
                    temp       = array(i)
                    array(i)   = array(i + 1)
                    array(i+1) = temp
                    swapped    = .true.
                end if
            end do
            if (.not. swapped) then
                return
            end if
        end do
      
    end subroutine sort_bubble
    
    function sort_mode(array) result(val)
        integer, intent(in) :: array(:)
        integer             :: val
        integer             :: temp_array(size(array))
        integer             :: max_count, count, i
        temp_array = array
        call sort_bubble(temp_array) 
        val        = temp_array(1)
        max_count  = 1
        count      = 1
        do i = 2, size(temp_array)
            if (temp_array(i) == temp_array(i - 1)) then
                count = count + 1
            else
                if (count > max_count) then
                    max_count = count
                    val       = temp_array(i - 1)
                end if
                count = 1
            end if
        end do
        if (count > max_count) then
            val = temp_array(size(temp_array))
        end if
    
    end function sort_mode


end module sort
