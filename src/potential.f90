module hqt_potential
    use hqt_constants, only: hqt_dp
    use hqt_tool,      only: hqt_range
    use hqt_type,      only: hqt_heterostructure
    implicit none
    
    public :: hqt_add_well, hqt_add_well_curve_width, hqt_add_well_curve_depth
    
    interface hqt_add_well
        module procedure hqt_add_well_array, hqt_add_well_heterostructure 
    end interface hqt_add_well  
    
contains
    subroutine hqt_add_well_array(array, initial_depth, final_depth, start_layer, stop_layer)
        real(hqt_dp),           intent(inout) :: array(:)
        real(hqt_dp),           intent(in)    :: initial_depth
        real(hqt_dp), optional, intent(in)    :: final_depth
        integer,      optional, intent(in)    :: start_layer
        integer,      optional, intent(in)    :: stop_layer
        real(hqt_dp)                          :: depth
        integer                               :: start
        integer                               :: stop
        depth = initial_depth
        start = 1
        stop = size(array)
        if (present(final_depth)) depth = final_depth
        if (present(start_layer)) start = start_layer
        if (present(stop_layer)) stop = stop_layer
        array(start:stop) = array(start:stop) + hqt_range(initial_depth, depth, start - stop)
    
    end subroutine hqt_add_well_array    
    
    subroutine hqt_add_well_heterostructure(heterostructure, initial_depth, final_depth, start_layer, stop_layer)
        type(hqt_heterostructure), intent(inout) :: heterostructure
        real(hqt_dp),              intent(in)    :: initial_depth
        real(hqt_dp), optional,    intent(in)    :: final_depth
        integer,      optional,    intent(in)    :: start_layer
        integer,      optional,    intent(in)    :: stop_layer
        real(hqt_dp)                             :: depth
        integer                                  :: start
        integer                                  :: stop
        depth = initial_depth
        start = 1
        stop = size(heterostructure%potential)
        if (present(final_depth)) depth = final_depth
        if (present(start_layer)) start = start_layer
        if (present(stop_layer)) stop = stop_layer
        heterostructure%potential(start:stop) = heterostructure%potential(start:stop) &
                                                    + hqt_range(initial_depth, depth, start - stop)
    
    end subroutine hqt_add_well_heterostructure   
    
    subroutine hqt_add_well_curve_width(array, start_ly,initial_depth, a, c, width)
        ! for varying width 
        real(hqt_dp),  intent(inout) :: array(:)
        integer, intent(in)    :: start_ly !end_ly,       
        real(hqt_dp),  intent(in)    :: initial_depth !, final_depth
        real(hqt_dp), allocatable    :: potential(:)
        integer                :: i, start, finish,   end_ly, approx_ly_invl
        real(hqt_dp)                 :: b
        real(hqt_dp), intent(in)     :: a, c
        real(hqt_dp), intent(out)    :: width
          

        b = - LOG(c-initial_depth)/LOG(a)

        
        !approx_ly_invl = int(EXP(LOG(1/(abs(initial_depth)/1000))/b) -a )
        approx_ly_invl = int( EXP( -LOG(c) /b  ) -a )

       
        end_ly = start_ly + approx_ly_invl
        
        width = EXP( -LOG(c) /b  ) -a +1 
        
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
    
    end subroutine hqt_add_well_curve_width        
    
   
   
    subroutine hqt_add_well_curve_depth(array, start_ly,end_ly, b, c, depth)
        ! for varying depth
        real(hqt_dp),  intent(inout) :: array(:)
        integer, intent(in)    :: start_ly, end_ly       
        real(hqt_dp), allocatable    :: potential(:)
        integer                :: i, start, finish
        real(hqt_dp)                 :: a
        real(hqt_dp), intent(in)     :: b, c
        real(hqt_dp)                 :: width
        real(hqt_dp), intent(out)    :: depth
        
        width = end_ly - start_ly + 1
        
        
        a = EXP(LOG(1/c)/b ) - width
        
        depth = - 1/((a) ** b) + c

     
        allocate(potential(abs(end_ly - start_ly) + 1))
        
        start = min(start_ly, end_ly)
        finish = max(start_ly, end_ly)
        potential = 0d0
        
        
        if (finish /= start) then
            do i = 1, finish - start + 1
                
                potential(i) = (-1/(( (i-1) + a)**(b))) + c
            end do
        end if
        do i = 1, size(potential)
            if (1 <= start + i - 1 .and. start + i - 1 <= size(array)) then
                array(start + i - 1) = array(start + i -1) + potential(i)
            else 
                print*, " Warning! The layer no. inputted is out of bound for potential. The well has been cropped. "
            end if
        end do
        
    end subroutine hqt_add_well_curve_depth    

    
end module hqt_potential
