module hqt_spectral
use hqt_constants,  only: hqt_dp
use hqt_type,       only: hqt_heterostructure
use hqt_tool,       only: hqt_range
implicit none

private

public :: hqt_spectral_function


interface hqt_spectral_function
    module procedure hqt_spectral_function_none, &
                     hqt_spectral_function_range
end interface


contains
    function hqt_spectral_function_range(energy, weight, structure, energy_range, broadening, &
                                   min_layer, max_layer, min_band, max_band, min_state, max_state) result(spectral_function)
        type(hqt_heterostructure(*)), intent(in) :: structure
        real(hqt_dp), intent(in)    :: energy(:)
        complex(hqt_dp), intent(in) :: weight(:, :)
        real(hqt_dp), intent(in)    :: energy_range(:)
        real(hqt_dp), intent(in)    :: broadening
        integer, intent(in)     :: min_layer
        integer, intent(in)     :: max_layer
        integer, intent(in)     :: min_band
        integer, intent(in)     :: max_band
        integer, intent(in)     :: min_state
        integer, intent(in)     :: max_state
        integer                 :: l, b, n, e
        real(hqt_dp)                :: spectral_function(max_layer - min_layer + 1, &
                                                     max_band - min_band + 1, &
                                                     max_state - min_state + 1, &
                                                     size(energy_range))
        do l = min_layer, max_layer
            do b = min_band, max_band
                do n = min_state, max_state
                    do e = 1, size(energy_range)
                        spectral_function(l, b, n, e) = &
                            -1.0_hqt_dp * aimag(weight((l - 1) * structure%material%num_bands + b, n)**2 &
                             / (energy_range(e) - energy(n) + cmplx(0.0_hqt_dp, broadening, kind = hqt_dp)))
                    end do
                end do
            end do
        end do
    
    end function hqt_spectral_function_range

    function hqt_spectral_function_none(energy, weight, structure, energy_range, broadening) result(spectral_function)
        type(hqt_heterostructure(*)), intent(in) :: structure
        real(hqt_dp), intent(in)    :: energy(:)
        complex(hqt_dp), intent(in) :: weight(:, :)
        real(hqt_dp), intent(in)    :: energy_range(:)
        real(hqt_dp), intent(in)    :: broadening
        integer                 :: l, b, n
        real(hqt_dp)                :: spectral_function(structure%num_layers, &
                                                     structure%material%num_bands, &
                                                     size(energy), &
                                                     size(energy_range))
        spectral_function = hqt_spectral_function(energy, weight, structure, energy_range, broadening, &
                                   1, size(spectral_function(:, 1, 1, 1)), &
                                   1, size(spectral_function(1, :, 1, 1)), &
                                   1, size(spectral_function(1, 1, :, 1)))
                                   
    end function hqt_spectral_function_none
    

end module hqt_spectral
