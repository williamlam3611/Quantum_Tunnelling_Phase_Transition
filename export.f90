module export
implicit none
    
private

public  :: export_hstack, export_vstack, export_create_dir, export_to_string

interface export_hstack
    module procedure export_hstack_character, export_hstack_character_array, export_hstack_character_array2, &
                     export_hstack_integer,   export_hstack_integer_array,   export_hstack_integer_array2, &
                     export_hstack_double,    export_hstack_double_array,    export_hstack_double_array2, &
                     export_hstack_double_array_sixteen, export_hstack_double_sixteen, export_hstack_double_array2_sixteen
end interface export_hstack

interface export_vstack
    module procedure export_vstack_character, export_vstack_character_array, export_vstack_character_array2, &
                     export_vstack_integer,   export_vstack_integer_array,   export_vstack_integer_array2, &
                     export_vstack_double,    export_vstack_double_array,    export_vstack_double_array2, &
                     export_vstack_double_array2_sixteen, export_vstack_double_array_sixteen, export_vstack_double_sixteen
end interface export_vstack

interface export_to_string
    module procedure export_to_string_integer, export_to_string_double, export_to_string_logical
end interface export_to_string


contains
    subroutine export_open_append(path, file_number)
        character(*), intent(in)  :: path
        integer,      intent(out) :: file_number
        logical                   :: exists
        inquire(file = trim(path), exist = exists)
        if (exists) then
            open(newunit = file_number, file = trim(path), status = "old", position = "append", action = "write", &
                 form = "formatted")
        else
            open(newunit = file_number, file = trim(path), status = "new", action = "write", form = "formatted")
        end if
    
    end subroutine export_open_append
    
    subroutine export_open_read(path, file_number)
        character(*), intent(in)  :: path
        integer,      intent(out) :: file_number
        logical                   :: exists
        inquire(file = trim(path), exist = exists)
        if (exists) then
            open(newunit = file_number, file = trim(path), status = "old", access = "sequential", &
                 action = "read", form = "formatted")
        else
            open(newunit = file_number, file = trim(path), status = "new", access = "sequential", &
                 action = "read", form = "formatted")
        end if
    
    end subroutine export_open_read
    
    function export_create_dir(path, prefix) result(dir)
        character(*), intent(in)  :: path
        character(*), intent(in)  :: prefix
        integer                   :: file_label
        logical                   :: directory_exists
        character(:), allocatable :: dir
        inquire(file = path//prefix//"/.", exist = directory_exists)
        if (.not. directory_exists) then
            dir = path//prefix//"/"
            call execute_command_line("mkdir -p """//dir//"""")
        else
            file_label       = 0
            do while (directory_exists)
                file_label = file_label + 1
                inquire(file = path//prefix//export_to_string(file_label)//"/.", exist = directory_exists)
            end do 
            dir = path//prefix//export_to_string(file_label)//"/"
            call execute_command_line("mkdir -p """//dir//"""")
        end if
        
    
    end function export_create_dir
        
    function export_to_string_integer(value, filler) result(string)
        integer,      intent(in)  :: value
        character,    optional    :: filler
        character                 :: filler_character
        character(16)             :: buffer
        character(:), allocatable :: string
        integer                   :: i
        filler_character = " "
        if (present(filler)) filler_character = filler
        write(buffer, "(I16)") value
        i = index(buffer, " ", back = .true.)
        buffer = adjustl(buffer)
        string = trim(adjustl(repeat(filler_character, i)//buffer))
    
    end function export_to_string_integer
    
    function export_to_string_double(value, filler) result(string)
        real*8,       intent(in)  :: value
        character,    optional    :: filler
        character                 :: filler_character
        character(16)             :: buffer
        character(:), allocatable :: string
        integer                   :: i
        filler_character = " "
        if (present(filler)) filler_character = filler
        write(buffer, "(F16.8)") value
        i = index(buffer, " ", back = .true.)
        buffer = adjustl(buffer)
        string = trim(adjustl(repeat(filler_character, i)//buffer))
    
    end function export_to_string_double
    
    function export_to_string_logical(value) result(string)
        logical,      intent(in)  :: value
        character(16)             :: buffer
        character(:), allocatable :: string
        integer                   :: i
        if (value) then
            write(buffer, "(I1)") 1
        else
            write(buffer, "(I1)") 0
        end if
        string = trim(adjustl(buffer))
    
    end function export_to_string_logical
    
    subroutine export_vstack_character(path, data, fmt)
        character(*), intent(in)  :: path
        character(*), intent(in)  :: data
        character(*), optional    :: fmt
        if (present(fmt)) then
            call export_vstack(path, (/ data /), fmt)
        else
            call export_vstack(path, (/ data /))
        end if
        
    end subroutine export_vstack_character
    
    subroutine export_vstack_character_array(path, data, fmt)
        character(*), intent(in)  :: path
        character(*), intent(in)  :: data(:)
        character(*), optional    :: fmt
        if (present(fmt)) then
            call export_vstack(path, reshape(data, (/ 1, size(data)/)), fmt)
        else
            call export_vstack(path, reshape(data, (/ 1, size(data)/)))
        end if
        
    end subroutine export_vstack_character_array
    
    subroutine export_vstack_character_array2(path, data, fmt)
        character(*), intent(in)  :: path
        character(*), intent(in)  :: data(:, :)
        character(*), optional    :: fmt
        character(:), allocatable :: format
        integer                   :: file_number
        integer                   :: i, j
        format = "A"
        if (present(fmt)) then
            format = fmt
        end if
        call export_open_append(path, file_number)
        do i = 1, size(data(:, 1))
            do j = 1, size(data(1, :))
                write(file_number, fmt = "("//format//")", advance = "no") data(i, j)
            end do
            write(file_number, fmt = *)
        end do
        close(file_number)
    
    end subroutine export_vstack_character_array2
    
    subroutine export_vstack_integer(path, data, fmt)
        character(*), intent(in)  :: path
        integer,      intent(in)  :: data
        character(*), optional    :: fmt
        if (present(fmt)) then
            call export_vstack(path, (/ data /), fmt)
        else
            call export_vstack(path, (/ data /))
        end if
        
    end subroutine export_vstack_integer
    
    subroutine export_vstack_integer_array(path, data, fmt)
        character(*), intent(in)  :: path
        integer,      intent(in)  :: data(:)
        character(*), optional    :: fmt
        if (present(fmt)) then
            call export_vstack(path, reshape(data, (/ 1, size(data)/)), fmt)
        else
            call export_vstack(path, reshape(data, (/ 1, size(data)/)))
        end if
        
    end subroutine export_vstack_integer_array
    
    subroutine export_vstack_integer_array2(path, data, fmt)
        character(*), intent(in)  :: path
        integer,      intent(in)  :: data(:, :)
        character(*), optional    :: fmt
        character(:), allocatable :: format
        integer                   :: file_number
        integer                   :: i, j
        format = "I16"
        if (present(fmt)) then
            format = fmt
        end if
        call export_open_append(path, file_number)
        do i = 1, size(data(:, 1))
            do j = 1, size(data(1, :))
                write(file_number, fmt = "("//format//")", advance = "no") data(i, j)
            end do
            write(file_number, fmt = *)
        end do
        close(file_number)
    
    end subroutine export_vstack_integer_array2
    
    subroutine export_vstack_double(path, data, fmt)
        character(*), intent(in)  :: path
        real*8,       intent(in)  :: data
        character(*), optional    :: fmt
        if (present(fmt)) then
            call export_vstack(path, (/ data /), fmt)
        else
            call export_vstack(path, (/ data /))
        end if
        
    end subroutine export_vstack_double
    
    subroutine export_vstack_double_sixteen(path, data, fmt)
        character(*), intent(in)  :: path
        real*16,       intent(in)  :: data
        character(*), optional    :: fmt
        if (present(fmt)) then
            call export_vstack(path, (/ data /), fmt)
        else
            call export_vstack(path, (/ data /))
        end if
        
    end subroutine export_vstack_double_sixteen
    
    subroutine export_vstack_double_array(path, data, fmt)
        character(*), intent(in)  :: path
        real*8,       intent(in)  :: data(:)
        character(*), optional    :: fmt
        if (present(fmt)) then
            call export_vstack(path, reshape(data, (/ 1, size(data)/)), fmt)
        else
            call export_vstack(path, reshape(data, (/ 1, size(data)/)))
        end if
        
    end subroutine export_vstack_double_array
    
    subroutine export_vstack_double_array_sixteen(path, data, fmt)
        character(*), intent(in)  :: path
        real*16,       intent(in)  :: data(:)
        character(*), optional    :: fmt
        if (present(fmt)) then
            call export_vstack(path, reshape(data, (/ 1, size(data)/)), fmt)
        else
            call export_vstack(path, reshape(data, (/ 1, size(data)/)))
        end if
        
    end subroutine export_vstack_double_array_sixteen
    
    subroutine export_vstack_double_array2(path, data, fmt)
        character(*), intent(in)  :: path
        real*8,       intent(in)  :: data(:, :)
        character(*), optional    :: fmt
        character(:), allocatable :: format
        integer                   :: file_number
        integer                   :: i, j
        format = "F64.16"!"F16.8"
        if (present(fmt)) then
            format = fmt
        end if
        call export_open_append(path, file_number)
        do i = 1, size(data(:, 1))
            do j = 1, size(data(1, :))
                write(file_number, fmt = "("//format//")", advance = "no") data(i, j)
            end do
            write(file_number, fmt = *)
        end do
        close(file_number)
    
    end subroutine export_vstack_double_array2
    
        subroutine export_vstack_double_array2_sixteen(path, data, fmt)
        character(*), intent(in)  :: path
        real*16,       intent(in)  :: data(:, :)
        character(*), optional    :: fmt
        character(:), allocatable :: format
        integer                   :: file_number
        integer                   :: i, j
        format = "F64.16"!"F16.8"
        if (present(fmt)) then
            format = fmt
        end if
        call export_open_append(path, file_number)
        do i = 1, size(data(:, 1))
            do j = 1, size(data(1, :))
                write(file_number, fmt = "("//format//")", advance = "no") data(i, j)
            end do
            write(file_number, fmt = *)
        end do
        close(file_number)
    
    end subroutine export_vstack_double_array2_sixteen
    
    subroutine export_hstack_character(path, data, fmt)
        character(*), intent(in)  :: path
        character(*), intent(in)  :: data
        character(*), optional    :: fmt
        if (present(fmt)) then
            call export_hstack(path, (/ data /), fmt)
        else
            call export_hstack(path, (/ data /))
        end if
        
    end subroutine export_hstack_character
    
    subroutine export_hstack_character_array(path, data, fmt)
        character(*), intent(in)  :: path
        character(*), intent(in)  :: data(:)
        character(*), optional    :: fmt
        if (present(fmt)) then
            call export_hstack(path, reshape(data, (/ 1, size(data)/)), fmt)
        else
            call export_hstack(path, reshape(data, (/ 1, size(data)/)))
        end if
        
    end subroutine export_hstack_character_array
    
    subroutine export_hstack_character_array2(path, data, fmt)
        character(*), intent(in)  :: path
        character(*), intent(in)  :: data(:, :)
        character(*), optional    :: fmt
        character(:), allocatable :: format
        character                 :: c
        integer                   :: status
        integer                   :: file_number_read, file_number_write
        integer                   :: i, j
        format = "A"
        if (present(fmt)) then
            format = fmt
        end if
        status = 0
        call export_open_read(path, file_number_read)
        call export_open_append(path//".temp", file_number_write)
        do i = 1, size(data(1, :))
            read(file_number_read, fmt = "(A1)", iostat = status, advance = "no") c
            do while (status == 0)
                write(file_number_write, fmt = "(A1)", advance = "no") c
                read(file_number_read, fmt = "(A1)", advance = "no", iostat = status) c
            end do
            backspace(file_number_read)
            read(file_number_read, fmt = *, iostat = status)
            do j = 1, size(data(:, 1))
                write(file_number_write, fmt = "("//format//")", advance = "no") data(j, i)
            end do
            write(file_number_write, fmt = *)
        end do
        close(file_number_read, status='delete')
        close(file_number_write)
        call rename(path//".temp", path)
    
    end subroutine export_hstack_character_array2
    
    subroutine export_hstack_integer(path, data, fmt)
        character(*), intent(in)  :: path
        integer,      intent(in)  :: data
        character(*), optional    :: fmt
        if (present(fmt)) then
            call export_hstack(path, (/ data /), fmt)
        else
            call export_hstack(path, (/ data /))
        end if
        
    end subroutine export_hstack_integer
    
    subroutine export_hstack_integer_array(path, data, fmt)
        character(*), intent(in)  :: path
        integer,      intent(in)  :: data(:)
        character(*), optional    :: fmt
        if (present(fmt)) then
            call export_hstack(path, reshape(data, (/ 1, size(data)/)), fmt)
        else
            call export_hstack(path, reshape(data, (/ 1, size(data)/)))
        end if
        
    end subroutine export_hstack_integer_array
    
    subroutine export_hstack_integer_array2(path, data, fmt)
        character(*), intent(in)  :: path
        integer,      intent(in)  :: data(:, :)
        character(*), optional    :: fmt
        character(:), allocatable :: format
        character                 :: c
        integer                   :: status
        integer                   :: file_number_read, file_number_write
        integer                   :: i, j
        format = "I16"
        if (present(fmt)) then
            format = fmt
        end if
        status = 0
        call export_open_read(path, file_number_read)
        call export_open_append(path//".temp", file_number_write)
        do i = 1, size(data(1, :))
            read(file_number_read, fmt = "(A1)", iostat = status, advance = "no") c
            do while (status == 0)
                write(file_number_write, fmt = "(A1)", advance = "no") c
                read(file_number_read, fmt = "(A1)", advance = "no", iostat = status) c
            end do
            backspace(file_number_read)
            read(file_number_read, fmt = *, iostat = status)
            do j = 1, size(data(:, 1))
                write(file_number_write, fmt = "("//format//")", advance = "no") data(j, i)
            end do
            write(file_number_write, fmt = *)
        end do
        close(file_number_read, status='delete')
        close(file_number_write)
        call rename(path//".temp", path)
    
    end subroutine export_hstack_integer_array2
    
    subroutine export_hstack_double(path, data, fmt)
        character(*), intent(in)  :: path
        real*8,       intent(in)  :: data
        character(*), optional    :: fmt
        if (present(fmt)) then
            call export_hstack(path, (/ data /), fmt)
        else
            call export_hstack(path, (/ data /))
        end if
        
    end subroutine export_hstack_double
    
    
    subroutine export_hstack_double_sixteen(path, data, fmt)
        character(*), intent(in)  :: path
        real*16,       intent(in)  :: data
        character(*), optional    :: fmt
        if (present(fmt)) then
            call export_hstack(path, (/ data /), fmt)
        else
            call export_hstack(path, (/ data /))
        end if
        
    end subroutine export_hstack_double_sixteen
    
    
    subroutine export_hstack_double_array(path, data, fmt)
        character(*), intent(in)  :: path
        real*8,       intent(in)  :: data(:)
        character(*), optional    :: fmt
        if (present(fmt)) then
            call export_hstack(path, reshape(data, (/ 1, size(data)/)), fmt)
        else
            call export_hstack(path, reshape(data, (/ 1, size(data)/)))
        end if
        
    end subroutine export_hstack_double_array
    
    subroutine export_hstack_double_array_sixteen(path, data, fmt)
        character(*), intent(in)  :: path
        real*16,       intent(in)  :: data(:)
        character(*), optional    :: fmt
        if (present(fmt)) then
            call export_hstack(path, reshape(data, (/ 1, size(data)/)), fmt)
        else
            call export_hstack(path, reshape(data, (/ 1, size(data)/)))
        end if
        
    end subroutine export_hstack_double_array_sixteen
    
    subroutine export_hstack_double_array2(path, data, fmt)
        character(*), intent(in)  :: path
        real*8,       intent(in)  :: data(:, :)
        character(*), optional    :: fmt
        character(:), allocatable :: format
        character                 :: c
        integer                   :: status
        integer                   :: file_number_read, file_number_write
        integer                   :: i, j
        format = "F64.16" !"F16.8"
        if (present(fmt)) then
            format = fmt
        end if
        status = 0
        call export_open_read(path, file_number_read)
        call export_open_append(path//".temp", file_number_write)
        do i = 1, size(data(1, :))
            read(file_number_read, fmt = "(A1)", iostat = status, advance = "no") c
            do while (status == 0)
                write(file_number_write, fmt = "(A1)", advance = "no") c
                read(file_number_read, fmt = "(A1)", advance = "no", iostat = status) c
            end do
            backspace(file_number_read)
            read(file_number_read, fmt = *, iostat = status)
            do j = 1, size(data(:, 1))
                write(file_number_write, fmt = "("//format//")", advance = "no") data(j, i)
            end do
            write(file_number_write, fmt = *)
        end do
        close(file_number_read, status='delete')
        close(file_number_write)
        call rename(path//".temp", path)
    
    end subroutine export_hstack_double_array2
    
    subroutine export_hstack_double_array2_sixteen(path, data, fmt)
        character(*), intent(in)  :: path
        real*16,       intent(in)  :: data(:, :)
        character(*), optional    :: fmt
        character(:), allocatable :: format
        character                 :: c
        integer                   :: status
        integer                   :: file_number_read, file_number_write
        integer                   :: i, j
        format = "F64.16" !"F16.8"
        if (present(fmt)) then
            format = fmt
        end if
        status = 0
        call export_open_read(path, file_number_read)
        call export_open_append(path//".temp", file_number_write)
        do i = 1, size(data(1, :))
            read(file_number_read, fmt = "(A1)", iostat = status, advance = "no") c
            do while (status == 0)
                write(file_number_write, fmt = "(A1)", advance = "no") c
                read(file_number_read, fmt = "(A1)", advance = "no", iostat = status) c
            end do
            backspace(file_number_read)
            read(file_number_read, fmt = *, iostat = status)
            do j = 1, size(data(:, 1))
                write(file_number_write, fmt = "("//format//")", advance = "no") data(j, i)
            end do
            write(file_number_write, fmt = *)
        end do
        close(file_number_read, status='delete')
        close(file_number_write)
        call rename(path//".temp", path)
    
    end subroutine export_hstack_double_array2_sixteen


end module export
