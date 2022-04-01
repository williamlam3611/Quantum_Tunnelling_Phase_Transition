module hqt_error
implicit none

public

type :: hqt_error_type(length)
    integer           :: code
    integer, len      :: length = 64
    character(length) :: message

end type hqt_error_type


interface operator (==)
    module procedure hqt_error_equals
end interface


type(hqt_error_type), parameter :: hqt_error_none = hqt_error_type(          code =  0, message = "No error")
type(hqt_error_type), parameter :: hqt_error_general = hqt_error_type(       code = -1, message = "General error")

type(hqt_error_type), parameter :: hqt_error_end_of_file = hqt_error_type(   code = -2, message = "End of file")
type(hqt_error_type), parameter :: hqt_error_not_found = hqt_error_type(     code = -3, message = "Not found")
type(hqt_error_type), parameter :: hqt_error_file_exists = hqt_error_type(   code = -4, message = "File exists")
type(hqt_error_type), parameter :: hqt_error_cannot_access = hqt_error_type( code = -5, message = "Cannot access")

contains
    pure function hqt_error_passed(err) result(res)
        type(hqt_error_type), intent(in) :: err
        logical                      :: res
        if (err%code >= 0) then
            res = .true.
        else
            res = .false.
        end if
    
    end function hqt_error_passed
    
    pure function hqt_error_failed(err) result(res)
        type(hqt_error_type), intent(in) :: err
        logical                      :: res
        if (err%code < 0) then
            res = .true.
        else
            res = .false.
        end if
    
    end function hqt_error_failed
    
    pure function hqt_error_equals(err1, err2) result(res)
        type(hqt_error_type), intent(in) :: err1
        type(hqt_error_type), intent(in) :: err2
        logical                      :: res
        if (err1%code == err2%code) then
            res = .true.
        else
            res = .false.
        end if
        
    end function hqt_error_equals


end module hqt_error
