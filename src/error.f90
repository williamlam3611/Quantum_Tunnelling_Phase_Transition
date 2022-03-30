module error
implicit none

public

type :: error_type(length)
    integer           :: code
    integer, len      :: length = 64
    character(length) :: message

end type error_type


interface operator (==)
    module procedure error_equals
end interface


type(error_type), parameter :: error_none = error_type(          code =  0, message = "No error")
type(error_type), parameter :: error_general = error_type(       code = -1, message = "General error")

type(error_type), parameter :: error_end_of_file = error_type(   code = -2, message = "End of file")
type(error_type), parameter :: error_not_found = error_type(     code = -3, message = "Not found")
type(error_type), parameter :: error_file_exists = error_type(   code = -4, message = "File exists")
type(error_type), parameter :: error_cannot_access = error_type( code = -5, message = "Cannot access")

contains
    pure function error_passed(err) result(res)
        type(error_type), intent(in) :: err
        logical                      :: res
        if (err%code >= 0) then
            res = .true.
        else
            res = .false.
        end if
    
    end function error_passed
    
    pure function error_failed(err) result(res)
        type(error_type), intent(in) :: err
        logical                      :: res
        if (err%code < 0) then
            res = .true.
        else
            res = .false.
        end if
    
    end function error_failed
    
    pure function error_equals(err1, err2) result(res)
        type(error_type), intent(in) :: err1
        type(error_type), intent(in) :: err2
        logical                      :: res
        if (err1%code == err2%code) then
            res = .true.
        else
            res = .false.
        end if
        
    end function error_equals


end module
