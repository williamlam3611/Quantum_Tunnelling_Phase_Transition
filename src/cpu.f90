module hqt_cpu
    use mpi
    use hqt_constants, only: hqt_dp
    use hqt_tool,      only: hqt_read
    implicit none
    
    private
    
    public  :: hqt_start, hqt_stop, &
               hqt_get_id, hqt_get_num, hqt_get_master_id, hqt_is_master, &
               hqt_broadcast, hqt_sum, hqt_split_work, &
               hqt_send_double, hqt_recv_double, hqt_broadcast_file_content
    
    logical, private, save      :: hqt_started   = .false.
    integer, private, parameter :: hqt_master_id = 0
    integer, private, save      :: hqt_id, hqt_num
    integer, private            :: hqt_error
    
    interface hqt_sum
        module procedure hqt_sum_double, hqt_sum_integer
    end interface hqt_sum
    
    interface hqt_broadcast
        module procedure hqt_broadcast_integer, hqt_broadcast_double, hqt_broadcast_double_complex, hqt_broadcast_string
    end interface hqt_broadcast
    
contains
    subroutine hqt_start()
        call MPI_INIT(hqt_error)
        call MPI_COMM_RANK(MPI_COMM_WORLD, hqt_id, hqt_error)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, hqt_num, hqt_error)
        hqt_started = .true.
    
    end subroutine hqt_start
    
    subroutine hqt_stop()
        if (.not. hqt_started) then
            call hqt_start()
        end if
        call MPI_FINALIZE(hqt_error)
        
    end subroutine hqt_stop
    
    subroutine hqt_split_work(start, stop, ammount)
        integer, intent(in)  :: ammount
        integer, intent(out) :: start, stop
        if (.not. hqt_started) then
            call hqt_start()
        end if
        if (hqt_num - mod(hqt_id, hqt_num) < hqt_id) then
            start = ammount - (hqt_num - hqt_id) * (floor(dble(ammount) / dble(hqt_num)) + 1) + 1
            stop = start + floor(dble(ammount) / dble(hqt_num))
        else
            start = hqt_id * floor(dble(ammount) / dble(hqt_num)) + 1
            stop = start + floor(dble(ammount) / dble(hqt_num)) - 1
        end if
        
    end subroutine hqt_split_work
        
    subroutine hqt_sum_double(v_in, v_out, target_id)
        real(hqt_dp), intent(in)  :: v_in(..)
        real(hqt_dp), intent(out) :: v_out(..)
        integer, optional   :: target_id
        integer             :: id = hqt_master_id
        if (.not. hqt_started) then
            call hqt_start()
        end if
        if (present(target_id)) then
            id = target_id
        end if
        call MPI_REDUCE(v_in, v_out, size(v_in), MPI_DOUBLE_PRECISION, MPI_SUM, id, MPI_COMM_WORLD, hqt_error)
    
    end subroutine hqt_sum_double
    
    subroutine hqt_sum_integer(v_in, v_out, target_id)
        integer, intent(in)  :: v_in(..)
        integer, intent(out) :: v_out(..)
        integer, optional    :: target_id
        integer              :: id = hqt_master_id
        if (.not. hqt_started) then
            call hqt_start()
        end if
        if (present(target_id)) then
            id = target_id
        end if
        call MPI_REDUCE(v_in, v_out, size(v_in), MPI_INT, MPI_SUM, id, MPI_COMM_WORLD, hqt_error)
    
    end subroutine hqt_sum_integer
    
    subroutine hqt_broadcast_integer(v, length, target_id)
        integer, intent(in)    :: length
        integer, intent(inout) :: v(..)
        integer, optional      :: target_id
        integer                :: id = hqt_master_id
        if (.not. hqt_started) then
            call hqt_start()
        end if
        if (present(target_id)) then
            id = target_id
        end if
        call MPI_BCAST(v, length, MPI_INT, id, MPI_COMM_WORLD, hqt_error)
    
    end subroutine hqt_broadcast_integer
    
    subroutine hqt_broadcast_string(v, length, target_id)
        integer,      intent(in)                 :: length
        character(:), allocatable, intent(inout) :: v
        integer, optional                        :: target_id
        integer                                  :: id = hqt_master_id
        if (.not. hqt_started) then
            call hqt_start()
        end if
        if (present(target_id)) then
            id = target_id
        end if
        call MPI_BCAST(v, length, MPI_CHAR, id, MPI_COMM_WORLD, hqt_error)
    
    end subroutine hqt_broadcast_string
    
    subroutine hqt_broadcast_double(v, length, target_id)
        integer, intent(in)    :: length
        real(hqt_dp),  intent(inout) :: v(..)
        integer, optional      :: target_id
        integer                :: id = hqt_master_id 
        if (.not. hqt_started) then
            call hqt_start()
        end if
        if (present(target_id)) then
            id = target_id
        end if
        call MPI_BCAST(v, length, MPI_DOUBLE_PRECISION, id, MPI_COMM_WORLD, hqt_error)
    
    end subroutine hqt_broadcast_double
    
    subroutine hqt_broadcast_double_complex(v, length, target_id)
        integer,    intent(in)    :: length
        complex(hqt_dp), intent(inout) :: v(..)
        integer,    optional      :: target_id
        integer                   :: id = hqt_master_id 
        if (.not. hqt_started) then
            call hqt_start()
        end if
        if (present(target_id)) then
            id = target_id
        end if
        call MPI_BCAST(v, length, MPI_DOUBLE_COMPLEX, id, MPI_COMM_WORLD, hqt_error)
    
    end subroutine hqt_broadcast_double_complex
    
    function hqt_broadcast_file_content(file_name, target_id) result(content)
        character(*), intent(in)      :: file_name
        integer, optional, intent(in) :: target_id
        integer                       :: id = hqt_master_id 
        integer                       :: content_length
        character(:), allocatable     :: content
        if (hqt_get_id() == id) then
            content = hqt_read(file_name)
            content_length = len(content)
        end if
        call hqt_broadcast(content_length, 1, id)
        if (hqt_get_id() .ne. id) allocate(character(content_length) :: content)
        call hqt_broadcast(content, content_length, id)

    end function hqt_broadcast_file_content
    
    
    subroutine hqt_send_double(data, count, to, tag)
        real(hqt_dp), intent(inout)           :: data(:)
        integer, intent(in)        :: tag      
        integer                         :: ierr
        integer, intent(in)             :: to, count

        if (.not. hqt_started) then
            call hqt_start()
        end if
    
        call MPI_SEND(data, count, MPI_DOUBLE_PRECISION, to , tag, MPI_COMM_WORLD, ierr ) 

    end subroutine hqt_send_double
    
        subroutine hqt_recv_double(data, count, from, tag)
        real(hqt_dp), intent(inout)           :: data(:)
        integer, intent(in)        :: tag      
        integer                         :: ierr, status(MPI_STATUS_SIZE) 
        integer, intent(in)             :: from, count

        if (.not. hqt_started) then
            call hqt_start()
        end if
    
        call MPI_RECV(data, count, MPI_DOUBLE_PRECISION, from , tag, MPI_COMM_WORLD, status, ierr) 

    end subroutine hqt_recv_double
    
    
    
    
    function hqt_get_id()
        integer :: hqt_get_id 
        if (.not. hqt_started) then
            call hqt_start()
        end if
        hqt_get_id = hqt_id
    
    end function hqt_get_id
    
    function hqt_get_num()
        integer :: hqt_get_num
        if (.not. hqt_started) then
            call hqt_start()
        end if
        hqt_get_num = hqt_num
    
    end function hqt_get_num
    
    function hqt_get_master_id()
        integer :: hqt_get_master_id
        if (.not. hqt_started) then
            call hqt_start()
        end if
        hqt_get_master_id = hqt_master_id
        
    end function hqt_get_master_id
    
    function hqt_is_master()
        logical :: hqt_is_master
        if (.not. hqt_started) then
            call hqt_start()
        end if
        hqt_is_master = (hqt_get_id() == hqt_get_master_id())
        
    end function hqt_is_master

end module hqt_cpu
