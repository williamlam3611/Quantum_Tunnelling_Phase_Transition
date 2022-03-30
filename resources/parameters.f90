module parameters
use constants, only: dp
use error,     only: error_type, error_passed
use tool,      only: tool_token_read, tool_token_write, tool_write, tool_to_string, tool_to_real, tool_to_integer, tool_to_complex
implicit none

public

real(dp), save :: broadening   = 0.0025_dp
integer,  save :: num_energy   = 128


contains
    subroutine parameters_read(file_name)
        character(*), intent(in)  :: file_name
        character(:), allocatable :: string
        real(dp),     allocatable :: temp_real(:)
        integer,      allocatable :: temp_integer(:)
        temp_real = tool_to_real(tool_token_read(file_name, "broadening"))
        if (allocated(temp_real)) broadening = temp_real(1)
        
        temp_integer = tool_to_integer(tool_token_read(file_name, "num_energy"))
        if (allocated(temp_integer)) num_energy = temp_integer(1)
        
    end subroutine parameters_read
    
    subroutine parameters_write(file_name)
        character(*), intent(in)  :: file_name
        character(:), allocatable :: string
        call tool_write(file_name, parameters_serialise(), .true.)
    
    end subroutine parameters_write
             
    function parameters_serialise() result(string)
        character(:), allocatable :: string
        string = tool_token_write(string, "broadening",   tool_to_string(broadening))
        string = tool_token_write(string, "num_energy",   tool_to_string(num_energy))
    
    end function parameters_serialise
    
    subroutine parameters_deserialise(string)
        character(*), intent(in) :: string
        broadening   = minval(tool_to_real(tool_token_read(string, "broadening")), dim = 1)
        num_energy   = minval(tool_to_integer(tool_token_read(string, "num_energy")), dim = 1)
    
    end subroutine parameters_deserialise

end module parameters
