module potential
    implicit none
    
    public :: potential_add_well
    
contains
    subroutine potential_add_well(array, start_ly, end_ly, initial_depth, final_depth)
    real*8,  intent(inout) :: array(:)
    integer, intent(in)    :: start_ly, end_ly
    real*8,  intent(in)    :: initial_depth, final_depth
    real*8                 :: potential(abs(end_ly - start_ly) + 1)
    integer                :: i, start, finish
    start = min(start_ly, end_ly)
    finish = max(start_ly, end_ly)
    potential = initial_depth
    if (finish /= start) then
        do i = 1, finish - start + 1
            potential(i) = potential(i) + (final_depth - initial_depth) * (i - 1) / (finish - start)
        end do
    end if
    do i = 1, size(potential)
        if (1 <= start + i - 1 .and. start + i - 1 <= size(array)) then
            array(start + i - 1) = array(start + i -1) + potential(i)
        else 
            print*, " Warning! The layer no. inputted is out of bound for potential. The well has been cropped. "
        end if
    end do
    
    end subroutine potential_add_well    
    
    subroutine potential_add_well_curve(array, start_ly,initial_depth, a, c)
    real*8,  intent(inout) :: array(:)
    integer, intent(in)    :: start_ly !end_ly,       
    real*8,  intent(in)    :: initial_depth !, final_depth
    real*8, allocatable    :: potential(:)
    integer                :: i, start, finish,  approx_ly_invl, end_ly
    real*8                 :: b
    real*8, intent(in)     :: a, c
      

    b = - LOG(c-initial_depth)/LOG(a)

    
    !approx_ly_invl = int(EXP(LOG(1/(abs(initial_depth)/1000))/b) -a )
    approx_ly_invl = int( EXP( -LOG(c) /b  ) -a )

   
    end_ly = start_ly + approx_ly_invl
    
    allocate(potential(abs(end_ly - start_ly) + 1))
    
    start = min(start_ly, end_ly)
    finish = max(start_ly, end_ly)
    potential = 0d0!initial_depth
    
    
    if (finish /= start) then
        do i = 1, finish - start + 1
            !potential(i) = potential(i) + (final_depth - initial_depth) * (i - 1) / (finish - start)
            
            potential(i) =  (-1/(( (i-1) + a)**(b))) + c
        end do
    end if
    do i = 1, size(potential)
        if (1 <= start + i - 1 .and. start + i - 1 <= size(array)) then
            array(start + i - 1) = array(start + i -1) + potential(i)
        else 
            print*, " Warning! The layer no. inputted is out of bound for potential. The well has been cropped. "
        end if
    end do
    
    end subroutine potential_add_well_curve       
    
    
end module potential
