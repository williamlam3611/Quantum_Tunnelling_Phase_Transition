module tool
use iso_fortran_env, only : iostat_end
use constants, only: dp, sp
use error,     only: error_type, error_failed, error_passed, &
                     error_none, error_not_found, &
                     error_end_of_file, error_cannot_access, &
                     error_file_exists
implicit none

private

public :: tool_token_write, &
          tool_token_read, &
          tool_open_append, &
          tool_open_read, &
          tool_create_directory, &
          tool_range, &
          tool_to_string, &
          tool_count_unique, &
          tool_location, &
          tool_delete, &
          tool_to_integer, &
          tool_to_real, &
          tool_to_complex, &
          operator(+), &
          tool_write, &
          tool_read, &
          tool_index


interface tool_range
    module procedure tool_range_real_double, &
                     tool_range_real, &
                     tool_range_complex_double, &
                     tool_range_complex, &
                     tool_range_integer
end interface tool_range

interface operator(+)
    procedure concatinate
end interface

interface tool_token_write
    module procedure :: tool_token_write_string, tool_token_write_file
end interface tool_token_write

interface tool_to_string
    module procedure tool_to_string_integer, tool_to_string_integer_array, &
                     tool_to_string_real, tool_to_string_real_array, &
                     tool_to_string_complex, tool_to_string_complex_array
end interface tool_to_string

interface tool_count_unique
    module procedure tool_count_unique_integer
end interface tool_count_unique

interface tool_count
    module procedure tool_count_integer, tool_count_string
end interface tool_count

interface tool_location
    module procedure tool_location_integer
end interface tool_location

interface tool_index
    module procedure tool_index_integer
end interface tool_index


contains
    function tool_index_integer(array, value) result(ret)
        integer, intent(in) :: array(:)
        integer, intent(in) :: value
        integer             :: i
        integer             :: ret
        do i = 1, size(array)
            if (array(i) == value) then
                ret = i
                return
            end if
        end do
        ret = 0
        
    end function tool_index_integer

    subroutine tool_write(file_name, string, replace)
        character(*), intent(in)      :: file_name
        character(*), intent(in)      :: string
        logical, optional, intent(in) :: replace
        integer                       :: file_number
        if (present(replace) .and. replace) call tool_delete(file_name)
        call tool_open_append(file_name, file_number)
        write(file_number, "(A)") string
        close(file_number)
        
    end subroutine tool_write
    
    function tool_read(file_name) result(string)
        character(*), intent(in)      :: file_name
        type(error_type)              :: file_error
        integer                       :: i
        integer                       :: file_number
        character(:), allocatable     :: string
        call tool_open_read_stream(file_name, file_number, err = file_error)
        if (error_failed(file_error)) then
            string = ""
            return
        end if
        inquire(unit = file_number, size = i)
        allocate(character(i) :: string)
        read(file_number) string
        close(file_number)
        
    end function tool_read

    function tool_token_write_string(token, value) result(string)
        character(*), intent(in)  :: token
        character(*), intent(in)  :: value
        character(:), allocatable :: string
        string = trim(adjustl(token)) + " = " + trim(adjustl(value)) + ";"
        
    end function tool_token_write_string

    function tool_token_write_file(path, token, value) result(string)
        character(*),     intent(in) :: path
        character(*),     intent(in) :: token
        character(*),     intent(in) :: value
        integer                      :: file_number
        type(error_type)             :: file_err
        character(:), allocatable    :: string
        call tool_open_append(path, file_number, file_err)
        if (error_failed(file_err)) then
            return
        end if
        
        string = tool_token_write(token, value)
        write(file_number, fmt = "(A)") string
        close(file_number)
    
    end function tool_token_write_file
    
    function tool_token_read(string_path, token) result (value)
        character(*), intent(in)  :: token
        character(*), intent(in)  :: string_path
        integer                   :: file_number
        type(error_type)          :: file_err
        character(:), allocatable :: value
        call tool_open_read(string_path, file_number, file_err)
        close(file_number)
        if (error_passed(file_err)) then
            value = tool_token_read_file(string_path, token)
        else
            value = tool_token_read_string(string_path, token)
        end if
        
    end function tool_token_read
    
    function tool_token_read_string(string, token) result(value)
        character(*), intent(in)  :: token
        character(*), intent(in)  :: string
        integer                   :: i
        integer                   :: j
        character(:), allocatable :: buffer
        character(:), allocatable :: value
        i = index(string, trim(adjustl(token))//" ") + len(trim(adjustl(token)))
        if (i == len(trim(adjustl(token)))) i = index(string, trim(adjustl(token))//"=") + len(trim(adjustl(token)))
        if (i == len(trim(adjustl(token)))) then
            value = ""
            return
        end if
        j = index(string(i:), "=") + i
        if (j == i) then
            value = ""
            return
        end if
        
        i = j
        do while (.true.)
            if (index(string(i:), ";") == 0 .or. i > len(string)) then
                i = len(string) - 1
                exit
            end if
            i = index(string(i:), ";") + i - 2
            if (tool_count(string(j:i), "{") == tool_count(string(j:i), "}")) then
                exit
            else
                i = i + 2
            end if
        end do
        
        value = trim(adjustl(string(j:i)))
        if (value(1:1) == "{") value = trim(adjustl(value(2:len(value) - 1)))
        
    end function tool_token_read_string
    
    function tool_token_read_file(path, token) result(value)
        character(*),     intent(in) :: path
        character(*),     intent(in) :: token
        character                    :: c
        integer                      :: file_number
        integer                      :: i
        type(error_type)             :: file_err
        character(:), allocatable    :: value
        call tool_open_read(path, file_number, file_err)
        if (error_failed(file_err)) then
            value = ""
            return
        end if
        
        do while (.true.)
            read(file_number, fmt = "(A1)", iostat = i, advance = "no") c
            if (i == iostat_end) then
                value = ""
                exit
            end if
            if (c == "!" .or. c == "#" .or. c == "\n" .or. c == "\r") then
                read(file_number, fmt = *)
                cycle
            end if
            
            if (c == ";") then
                value = tool_token_read_string(value//";", token)
                if (value == "") cycle
                exit
            else
                value = value + c
            end if
        end do
        
        close(file_number)
        
    end function tool_token_read_file
    
    
    function concatinate(str1, str2) result(str)
        character(*), intent(in)         :: str1
        character(*), intent(in)         :: str2
        character(len(str1) + len(str2)) :: str
        str = str1//str2
        
    end function concatinate
    
    
    subroutine tool_delete(path, err)
        character(*), intent(in)                :: path
        integer                                 :: file_number
        integer                                 :: iostat
        logical                                 :: exists
        type(error_type), optional, intent(out) :: err
        if (present(err)) err = error_none
        inquire(file = trim(path), exist = exists)
        if (exists) then
            open(newunit = file_number, iostat = iostat, file = trim(path), status='old')
            if (iostat == 0) then
                close(file_number, status = "delete")
            else
                if (present(err)) err = error_cannot_access
            end if
        else
            if (present(err)) err = error_not_found
        end if
        
    end subroutine tool_delete
    
    
    subroutine tool_open_append(path, file_number, err)
        character(*),     intent(in)            :: path
        logical                                 :: exists        
        integer,          intent(out)           :: file_number
        type(error_type), intent(out), optional :: err
        if (present(err)) err = error_none
        inquire(file = trim(path), exist = exists)
        if (exists) then
            open(newunit = file_number, file = trim(path), status = "old", position = "append", action = "write", &
                 form = "formatted")
        else
            open(newunit = file_number, file = trim(path), status = "new", action = "write", form = "formatted")
        end if
    
    end subroutine tool_open_append
    
    
    subroutine tool_open_read(path, file_number, err)
        character(*),     intent(in)            :: path
        logical                                 :: exists
        integer,          intent(out)           :: file_number
        type(error_type), intent(out), optional :: err
        if (present(err)) err = error_none
        inquire(file = trim(path), exist = exists)
        if (exists) then
            open(newunit = file_number, file = trim(path), status = "old", access = "sequential", &
                 action = "read", form = "formatted")
        else
            if (present(err)) err = error_not_found
        end if
    
    end subroutine tool_open_read
    
    subroutine tool_open_read_stream(path, file_number, rec_length, err)
        character(*),     intent(in)            :: path
        integer, optional, intent(in)           :: rec_length
        logical                                 :: exists
        integer,          intent(out)           :: file_number
        type(error_type), intent(out), optional :: err
        if (present(err)) err = error_none
        inquire(file = trim(path), exist = exists)
        if (exists) then
            if (present(rec_length)) then
                open(newunit = file_number, file = trim(path), status = "old", access = "stream", &
                     action = "read", form = "unformatted")
            else
                open(newunit = file_number, file = trim(path), status = "old", access = "stream", &
                     action = "read", form = "unformatted")
            end if
            
        else
            if (present(err)) err = error_not_found
        end if
    
    end subroutine tool_open_read_stream
    
    
    subroutine tool_create_directory(path, err)
        character(*),     intent(in)            :: path
        character(:), allocatable               :: dir
        logical                                 :: exists
        type(error_type), intent(out), optional :: err
        if (present(err)) err = error_none
        dir = trim(adjustl(path))
        if (dir(len(dir):) .ne. "/") dir = dir//"/"
        inquire(file = dir//".", exist = exists)
        if (.not. exists) then
            call execute_command_line("mkdir -p """//dir//"""")
        else
            if (present(err)) err = error_file_exists
        end if
    
    end subroutine tool_create_directory
    
    
    function tool_count_unique_integer(array) result(number)
        integer, intent(in) :: array(:)
        integer             :: unique(size(array))
        integer             :: i, j
        integer             :: number
        number = 0
        unique = 0
        if (any(array == 0)) then
            number = number + 1
        end if
        do i = 1, size(array)
            if (array(i) /= 0 .and. .not. any(unique == array(i))) then
                number = number + 1
                do j = 1, size(unique)
                    if (unique(j) == 0) then
                        unique(j) = array(i)
                        exit
                    end if
                end do
            end if
        end do
    
    end function tool_count_unique_integer
    
    
    function tool_count_integer(array, value) result(number)
        integer, intent(in) :: array(:)
        integer, intent(in) :: value
        integer             :: i
        integer             :: number
        number = 0
        do i = 1, size(array)
            if (array(i) == value) number = number + 1
        end do
    
    end function tool_count_integer
    
    function tool_count_string(string, value) result(number)
        character(*), intent(in) :: string
        character(*), intent(in) :: value
        integer                  :: i
        integer                  :: j
        integer                  :: number
        number = 0
        j = 1
        do while (.true.)
            if (j > len(string)) exit
            i = index(string(j:), value)
            if (i == 0) exit
            number = number + 1
            j = j + i + len(value) - 1
        end do
    
    end function tool_count_string
    
    
    function tool_location_integer(array, val) result(location)
        integer, intent(in) :: array(:)
        integer, intent(in) :: val
        integer             :: i
        integer             :: location
        location = 0
        do i = 1, size(array)
            if (array(i) == val) then
                location = i
                exit
            end if
        end do
    
    end function tool_location_integer
    

    function tool_range_real_double(start, stop, number) result(array)
        integer, intent(in)  :: number
        real(dp), intent(in) :: start
        real(dp), intent(in) :: stop
        integer              :: i
        real(dp)             :: array(number)
        array = start
        do i = 2, number
            array(i) = array(i) + (i - 1) * (stop - start) / (number - 1)
        end do
    
    end function tool_range_real_double
    
    function tool_range_real(start, stop, number) result(array)
        integer, intent(in)  :: number
        real(sp), intent(in) :: start
        real(sp), intent(in) :: stop
        integer              :: i
        real(sp)             :: array(number)
        array = start
        do i = 2, number
            array(i) = array(i) + (i - 1) * (stop - start) / (number - 1)
        end do
    
    end function tool_range_real
    
    function tool_range_complex_double(start, stop, number) result(array)
        integer, intent(in)     :: number
        complex(dp), intent(in) :: start
        complex(dp), intent(in) :: stop
        integer                 :: i
        complex(dp)             :: array(number)
        array = start
        do i = 2, number
            array(i) = array(i) + (i - 1) * (stop - start) / (number - 1)
        end do
    
    end function tool_range_complex_double
    
    function tool_range_complex(start, stop, number) result(array)
        integer, intent(in)     :: number
        complex(sp), intent(in) :: start
        complex(sp), intent(in) :: stop
        integer                 :: i
        complex(sp)             :: array(number)
        array = start
        do i = 2, number
            array(i) = array(i) + (i - 1) * (stop - start) / (number - 1)
        end do
    
    end function tool_range_complex
    
    function tool_range_integer(start, stop, number) result(array)
        integer, intent(in) :: number
        integer, intent(in) :: start
        integer, intent(in) :: stop
        integer             :: i
        integer             :: array(number)
        array = start
        do i = 2, number
            array(i) = array(i) + (i - 1) * (stop - start) / (number - 1)
        end do
    
    end function tool_range_integer
    

    function tool_to_string_integer(value, fmt) result(string)
        integer,                intent(in) :: value
        character(*), optional, intent(in) :: fmt
        character(:), allocatable          :: string
        if (present(fmt)) then
            string = tool_to_string((/ value /), fmt)
        else
            string = tool_to_string((/ value /))
        end if
        
    end function tool_to_string_integer
    
    function tool_to_string_integer_array(value, fmt) result(string)
        integer,                intent(in) :: value(:)
        character(*), optional, intent(in) :: fmt
        integer                            :: i
        integer                            :: j
        character(:), allocatable          :: format
        character(:), allocatable          :: buffer
        character(:), allocatable          :: string
        if (present(fmt)) then
            format = trim(adjustl(fmt))
            if (format(1:1) == "(") format = format(2:len(format) - 1)
            format = trim(adjustl(format))
        else
            format = "I19"
        end if
        
        i = scan(format, "0123456789")
        j = scan(format, ".") - 1
        if (j == -1) j = len(format)
        read(format(i:j), fmt = *) i
        allocate(character(i) :: buffer)
        string = ""
        
        do i = 1, size(value)
            write(buffer, fmt = "("//format//")") value(i)
            string = string + trim(adjustl(buffer))
            if (i < size(value)) string = string + ", "
        end do
        string = trim(adjustl(string))
    
    end function tool_to_string_integer_array
    
    
    function tool_to_string_real(value, fmt) result(string)
        real(dp),               intent(in) :: value
        character(*), optional, intent(in) :: fmt
        character(:), allocatable          :: string
        if (present(fmt)) then
            string = tool_to_string((/ value /), fmt)
        else
            string = tool_to_string((/ value /))
        end if
        
    end function tool_to_string_real
    
    function tool_to_string_real_array(value, fmt) result(string)
        real(dp),               intent(in) :: value(:)
        character(*), optional, intent(in) :: fmt
        integer                            :: i
        integer                            :: j
        character(:), allocatable          :: format
        character(:), allocatable          :: buffer
        character(:), allocatable          :: string
        if (present(fmt)) then
            format = trim(adjustl(fmt))
            if (format(1:1) == "(") format = format(2:len(format) - 1)
            format = trim(adjustl(format))
        else
            format = "ES16.8"
        end if
        
        i = scan(format, "0123456789")
        j = scan(format, ".") - 1
        if (j == -1) j = len(format)
        read(format(i:j), fmt = *) i
        allocate(character(i) :: buffer)
        string = ""
        
        do i = 1, size(value)
            write(buffer, fmt = "("//format//")") value(i)
            string = string + trim(adjustl(buffer))
            if (i < size(value)) string = string + ", "
        end do
        string = trim(adjustl(string))
    
    end function tool_to_string_real_array
    
    
    function tool_to_string_complex(value, fmt) result(string)
        complex(dp),            intent(in) :: value
        character(*), optional, intent(in) :: fmt
        character(:), allocatable          :: string
        if (present(fmt)) then
            string = tool_to_string((/ value /), fmt)
        else
            string = tool_to_string((/ value /))
        end if
        
    end function tool_to_string_complex
    
    function tool_to_string_complex_array(value, fmt) result(string)
        complex(dp),            intent(in) :: value(:)
        character(*), optional, intent(in) :: fmt
        integer                            :: i
        integer                            :: j
        character(:), allocatable          :: format
        character(:), allocatable          :: buffer
        character(:), allocatable          :: string
        if (present(fmt)) then
            format = trim(adjustl(fmt))
            if (format(1:1) == "(") format = format(2:len(format) - 1)
            format = trim(adjustl(format))
        else
            format = "ES16.8, ES16.8"
        end if
        
        i = scan(format, "0123456789")
        j = scan(format, ".") - 1
        if (j == -1) j = len(format)
        read(format(i:j), fmt = *) i
        allocate(character(i) :: buffer)
        string = ""
        
        do i = 1, size(value)
            write(buffer, fmt = "("//format//")") real(value(i))
            string = string + "(" + trim(adjustl(buffer)) + ", "
            write(buffer, fmt = "("//format//")") aimag(value(i))
            string = string + trim(adjustl(buffer)) + ")"
            if (i < size(value)) string = string + ", "
        end do
        string = trim(adjustl(string))
    
    end function tool_to_string_complex_array
    
    
    function tool_to_integer(value) result(res)
        character(*), intent(in)  :: value
        integer                   :: i
        integer,      allocatable :: res(:)
        i = tool_count(value, ",") + 1
        allocate(res(i))
        read(value, fmt = *, iostat = i) res
        if (i .ne. 0) deallocate(res)
        
    end function tool_to_integer
    
    
    function tool_to_real(value) result(res)
        character(*), intent(in)  :: value
        integer                   :: i
        real(dp),     allocatable :: res(:)
        i = tool_count(value, ",") + 1
        allocate(res(i))
        read(value, fmt = *, iostat = i) res
        if (i .ne. 0) deallocate(res)
        
    end function tool_to_real
    
    
    function tool_to_complex(value) result(res)
        character(*), intent(in)  :: value
        integer                   :: i
        complex(dp),  allocatable :: res(:)
        i = (tool_count(value, ",") + 1) / 2
        allocate(res(i))
        read(value, fmt = *, iostat = i) res
        if (i .ne. 0) deallocate(res)
        
    end function tool_to_complex
    
    
end module tool
