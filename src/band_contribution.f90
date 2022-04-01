module hqt_contribution
use hqt_constants,  only: hqt_dp
use hqt_type,       only: hqt_heterostructure
implicit none

private

public :: hqt_band_contribution

interface hqt_band_contribution
    module procedure hqt_band_contribution_none, &
                     hqt_band_contribution_range
end interface


contains
    function hqt_band_contribution_range(energy, weight, structure, &
                                   min_layer, max_layer, min_band, max_band, min_state, max_state) result(band_contribution)
        type(hqt_heterostructure(*)), intent(in) :: structure
        real(hqt_dp), intent(in)    :: energy(:)
        complex(hqt_dp), intent(in) :: weight(:, :)
        integer, intent(in)     :: min_layer
        integer, intent(in)     :: max_layer
        integer, intent(in)     :: min_band
        integer, intent(in)     :: max_band
        integer, intent(in)     :: min_state
        integer, intent(in)     :: max_state
        integer                 :: l, b, n
        real(hqt_dp)                :: band_contribution(max_layer - min_layer + 1, &
                                                     max_band - min_band + 1, &
                                                     max_state - min_state + 1)
        do l = min_layer, max_layer
            do b = min_band, max_band
                do n = min_state, max_state
                    band_contribution(l, b, n) = weight((l - 1) * structure%material%num_bands + b, n)**2
                end do
            end do
        end do
    
    end function hqt_band_contribution_range

    function hqt_band_contribution_none(energy, weight, structure) result(band_contribution)
        type(hqt_heterostructure(*)), intent(in) :: structure
        real(hqt_dp), intent(in)    :: energy(:)
        complex(hqt_dp), intent(in) :: weight(:, :)
        real(hqt_dp)                :: band_contribution(structure%num_layers, &
                                                     structure%material%num_bands, &
                                                     size(energy))
        band_contribution = hqt_band_contribution(energy, weight, structure, &
                                   1, size(band_contribution(:, 1, 1)), &
                                   1, size(band_contribution(1, :, 1)), &
                                   1, size(band_contribution(1, 1, :)))
                                   
    end function hqt_band_contribution_none


end module hqt_contribution
