module hqt_tool
use iso_fortran_env, only : iostat_end
use hqt_constants, only: hqt_dp
use hqt_error,     only: hqt_error_type, hqt_error_failed, hqt_error_passed, &
                         hqt_error_none, hqt_error_not_found, &
                         hqt_error_cannot_access, hqt_error_file_exists
implicit none

private

public :: hqt_token_write, &
          hqt_token_read, &
          hqt_open_append, &
          hqt_open_read, &
          hqt_create_directory, &
          hqt_range, &
          hqt_to_string, &
          hqt_count_unique, &
          hqt_location, &
          hqt_delete, &
          hqt_to_integer, &
          hqt_to_real, &
          hqt_to_complex, &
          operator(+), &
          hqt_write, &
          hqt_read, &
          hqt_index


interface hqt_range
    module procedure hqt_range_real, &
                     hqt_range_complex, &
                     hqt_range_integer
end interface hqt_range

interface operator(+)
    procedure hqt_concatinate
end interface

interface hqt_token_write
    module procedure :: hqt_token_write_string, hqt_token_write_file
end interface hqt_token_write

interface hqt_to_string
    module procedure hqt_to_string_integer, hqt_to_string_integer_array, &
                     hqt_to_string_real, hqt_to_string_real_array, &
                     hqt_to_string_complex, hqt_to_string_complex_array
end interface hqt_to_string

interface hqt_count_unique
    module procedure hqt_count_unique_integer
end interface hqt_count_unique

interface hqt_count
    module procedure hqt_count_integer, hqt_count_string
end interface hqt_count

interface hqt_location
    module procedure hqt_location_integer
end interface hqt_location

interface hqt_index
    module procedure hqt_index_integer
end interface hqt_index


contains
    function hqt_index_integer(array, value) result(ret)
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
        
    end function hqt_index_integer

    subroutine hqt_write(file_name, string, replace)
        character(*), intent(in)      :: file_name
        character(*), intent(in)      :: string
        logical, optional, intent(in) :: replace
        integer                       :: file_number
        if (present(replace) .and. replace) call hqt_delete(file_name)
        call hqt_open_append(file_name, file_number)
        write(file_number, "(A)") string
        close(file_number)
        
    end subroutine hqt_write
    
    function hqt_read(file_name) result(string)
        character(*), intent(in)      :: file_name
        type(hqt_error_type)          :: file_error
        integer                       :: i
        integer                       :: file_number
        character(:), allocatable     :: string
        call hqt_open_read_stream(file_name, file_number, err = file_error)
        if (hqt_error_failed(file_error)) then
            string = ""
            return
        end if
        inquire(unit = file_number, size = i)
        allocate(character(i) :: string)
        read(file_number) string
        close(file_number)
        
    end function hqt_read

    function hqt_token_write_string(token, value) result(string)
        character(*), intent(in)  :: token
        character(*), intent(in)  :: value
        character(:), allocatable :: string
        string = trim(adjustl(token)) + " = " + trim(adjustl(value)) + ";"
        
    end function hqt_token_write_string

    function hqt_token_write_file(path, token, value) result(string)
        character(*),     intent(in) :: path
        character(*),     intent(in) :: token
        character(*),     intent(in) :: value
        integer                      :: file_number
        type(hqt_error_type)         :: file_err
        character(:), allocatable    :: string
        call hqt_open_append(path, file_number, file_err)
        if (hqt_error_failed(file_err)) then
            return
        end if
        
        string = hqt_token_write(token, value)
        write(file_number, fmt = "(A)") string
        close(file_number)
    
    end function hqt_token_write_file
    
    function hqt_token_read(string_path, token) result (value)
        character(*), intent(in)  :: token
        character(*), intent(in)  :: string_path
        integer                   :: file_number
        type(hqt_error_type)          :: file_err
        character(:), allocatable :: value
        call hqt_open_read(string_path, file_number, file_err)
        close(file_number)
        if (hqt_error_passed(file_err)) then
            value = hqt_token_read_file(string_path, token)
        else
            value = hqt_token_read_string(string_path, token)
        end if
        
    end function hqt_token_read
    
    function hqt_token_read_string(string, token) result(value)
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
            if (hqt_count(string(j:i), "{") == hqt_count(string(j:i), "}")) then
                exit
            else
                i = i + 2
            end if
        end do
        
        value = trim(adjustl(string(j:i)))
        if (value(1:1) == "{") value = trim(adjustl(value(2:len(value) - 1)))
        
    end function hqt_token_read_string
    
    function hqt_token_read_file(path, token) result(value)
        character(*),     intent(in) :: path
        character(*),     intent(in) :: token
        character                    :: c
        integer                      :: file_number
        integer                      :: i
        type(hqt_error_type)         :: file_err
        character(:), allocatable    :: value
        call hqt_open_read(path, file_number, file_err)
        if (hqt_error_failed(file_err)) then
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
                value = hqt_token_read_string(value//";", token)
                if (value == "") cycle
                exit
            else
                value = value + c
            end if
        end do
        
        close(file_number)
        
    end function hqt_token_read_file
    
    
    function hqt_concatinate(str1, str2) result(str)
        character(*), intent(in)         :: str1
        character(*), intent(in)         :: str2
        character(len(str1) + len(str2)) :: str
        str = str1//str2
        
    end function hqt_concatinate
    
    
    subroutine hqt_delete(path, err)
        character(*), intent(in)                :: path
        integer                                 :: file_number
        integer                                 :: iostat
        logical                                 :: exists
        type(hqt_error_type), optional, intent(out) :: err
        if (present(err)) err = hqt_error_none
        inquire(file = trim(path), exist = exists)
        if (exists) then
            open(newunit = file_number, iostat = iostat, file = trim(path), status='old')
            if (iostat == 0) then
                close(file_number, status = "delete")
            else
                if (present(err)) err = hqt_error_cannot_access
            end if
        else
            if (present(err)) err = hqt_error_not_found
        end if
        
    end subroutine hqt_delete
    
    
    subroutine hqt_open_append(path, file_number, err)
        character(*),     intent(in)            :: path
        logical                                 :: exists        
        integer,          intent(out)           :: file_number
        type(hqt_error_type), intent(out), optional :: err
        if (present(err)) err = hqt_error_none
        inquire(file = trim(path), exist = exists)
        if (exists) then
            open(newunit = file_number, file = trim(path), status = "old", position = "append", action = "write", &
                 form = "formatted")
        else
            open(newunit = file_number, file = trim(path), status = "new", action = "write", form = "formatted")
        end if
    
    end subroutine hqt_open_append
    
    
    subroutine hqt_open_read(path, file_number, err)
        character(*),     intent(in)            :: path
        logical                                 :: exists
        integer,          intent(out)           :: file_number
        type(hqt_error_type), intent(out), optional :: err
        if (present(err)) err = hqt_error_none
        inquire(file = trim(path), exist = exists)
        if (exists) then
            open(newunit = file_number, file = trim(path), status = "old", access = "sequential", &
                 action = "read", form = "formatted")
        else
            if (present(err)) err = hqt_error_not_found
        end if
    
    end subroutine hqt_open_read
    
    subroutine hqt_open_read_stream(path, file_number, rec_length, err)
        character(*),     intent(in)            :: path
        integer, optional, intent(in)           :: rec_length
        logical                                 :: exists
        integer,          intent(out)           :: file_number
        type(hqt_error_type), intent(out), optional :: err
        if (present(err)) err = hqt_error_none
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
            if (present(err)) err = hqt_error_not_found
        end if
    
    end subroutine hqt_open_read_stream
    
    
    subroutine hqt_create_directory(path, err)
        character(*),     intent(in)            :: path
        character(:), allocatable               :: dir
        logical                                 :: exists
        type(hqt_error_type), intent(out), optional :: err
        if (present(err)) err = hqt_error_none
        dir = trim(adjustl(path))
        if (dir(len(dir):) .ne. "/") dir = dir//"/"
        inquire(file = dir//".", exist = exists)
        if (.not. exists) then
            call execute_command_line("mkdir -p """//dir//"""")
        else
            if (present(err)) err = hqt_error_file_exists
        end if
    
    end subroutine hqt_create_directory
    
    
    function hqt_count_unique_integer(array) result(number)
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
    
    end function hqt_count_unique_integer
    
    
    function hqt_count_integer(array, value) result(number)
        integer, intent(in) :: array(:)
        integer, intent(in) :: value
        integer             :: i
        integer             :: number
        number = 0
        do i = 1, size(array)
            if (array(i) == value) number = number + 1
        end do
    
    end function hqt_count_integer
    
    function hqt_count_string(string, value) result(number)
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
    
    end function hqt_count_string
    
    
    function hqt_location_integer(array, val) result(location)
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
    
    end function hqt_location_integer
    
    
    function hqt_range_real(start, stop, number) result(array)
        integer, intent(in)      :: number
        real(hqt_dp), intent(in) :: start
        real(hqt_dp), intent(in) :: stop
        integer                  :: i
        real(hqt_dp)             :: array(number)
        array = start
        do i = 2, number
            array(i) = array(i) + (i - 1) * (stop - start) / (number - 1)
        end do
    
    end function hqt_range_real
    
    function hqt_range_complex(start, stop, number) result(array)
        integer, intent(in)     :: number
        complex(hqt_dp), intent(in) :: start
        complex(hqt_dp), intent(in) :: stop
        integer                 :: i
        complex(hqt_dp)             :: array(number)
        array = start
        do i = 2, number
            array(i) = array(i) + (i - 1) * (stop - start) / (number - 1)
        end do
    
    end function hqt_range_complex
    
    function hqt_range_integer(start, stop, number) result(array)
        integer, intent(in) :: number
        integer, intent(in) :: start
        integer, intent(in) :: stop
        integer             :: i
        integer             :: array(number)
        array = start
        do i = 2, number
            array(i) = array(i) + (i - 1) * (stop - start) / (number - 1)
        end do
    
    end function hqt_range_integer
    

    function hqt_to_string_integer(value, fmt) result(string)
        integer,                intent(in) :: value
        character(*), optional, intent(in) :: fmt
        character(:), allocatable          :: string
        if (present(fmt)) then
            string = hqt_to_string((/ value /), fmt)
        else
            string = hqt_to_string((/ value /))
        end if
        
    end function hqt_to_string_integer
    
    function hqt_to_string_integer_array(value, fmt) result(string)
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
    
    end function hqt_to_string_integer_array
    
    
    function hqt_to_string_real(value, fmt) result(string)
        real(hqt_dp),               intent(in) :: value
        character(*), optional, intent(in) :: fmt
        character(:), allocatable          :: string
        if (present(fmt)) then
            string = hqt_to_string((/ value /), fmt)
        else
            string = hqt_to_string((/ value /))
        end if
        
    end function hqt_to_string_real
    
    function hqt_to_string_real_array(value, fmt) result(string)
        real(hqt_dp),               intent(in) :: value(:)
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
    
    end function hqt_to_string_real_array
    
    
    function hqt_to_string_complex(value, fmt) result(string)
        complex(hqt_dp),            intent(in) :: value
        character(*), optional, intent(in) :: fmt
        character(:), allocatable          :: string
        if (present(fmt)) then
            string = hqt_to_string((/ value /), fmt)
        else
            string = hqt_to_string((/ value /))
        end if
        
    end function hqt_to_string_complex
    
    function hqt_to_string_complex_array(value, fmt) result(string)
        complex(hqt_dp),            intent(in) :: value(:)
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
    
    end function hqt_to_string_complex_array
    
    
    function hqt_to_integer(value) result(res)
        character(*), intent(in)  :: value
        integer                   :: i
        integer,      allocatable :: res(:)
        i = hqt_count(value, ",") + 1
        allocate(res(i))
        read(value, fmt = *, iostat = i) res
        if (i .ne. 0) deallocate(res)
        
    end function hqt_to_integer
    
    
    function hqt_to_real(value) result(res)
        character(*), intent(in)  :: value
        integer                   :: i
        real(hqt_dp),     allocatable :: res(:)
        i = hqt_count(value, ",") + 1
        allocate(res(i))
        read(value, fmt = *, iostat = i) res
        if (i .ne. 0) deallocate(res)
        
    end function hqt_to_real
    
    
    function hqt_to_complex(value) result(res)
        character(*), intent(in)  :: value
        integer                   :: i
        complex(hqt_dp),  allocatable :: res(:)
        i = (hqt_count(value, ",") + 1) / 2
        allocate(res(i))
        read(value, fmt = *, iostat = i) res
        if (i .ne. 0) deallocate(res)
        
    end function hqt_to_complex
    
    
end module hqt_tool
