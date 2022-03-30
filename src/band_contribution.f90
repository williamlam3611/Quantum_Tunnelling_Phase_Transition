module band_contribution
use constants,  only: dp
use type,       only: type_heterostructure
implicit none

private

public :: band_contribution_band_contribution

interface band_contribution_band_contribution
    module procedure band_contribution_band_contribution_none, &
                     band_contribution_band_contribution_range
end interface


contains
    function band_contribution_band_contribution_range(energy, weight, structure, &
                                   min_layer, max_layer, min_band, max_band, min_state, max_state) result(band_contribution)
        type(type_heterostructure(*)), intent(in) :: structure
        real(dp), intent(in)    :: energy(:)
        complex(dp), intent(in) :: weight(:, :)
        integer, intent(in)     :: min_layer
        integer, intent(in)     :: max_layer
        integer, intent(in)     :: min_band
        integer, intent(in)     :: max_band
        integer, intent(in)     :: min_state
        integer, intent(in)     :: max_state
        integer                 :: l, b, n
        real(dp)                :: band_contribution(max_layer - min_layer + 1, &
                                                     max_band - min_band + 1, &
                                                     max_state - min_state + 1)
        do l = min_layer, max_layer
            do b = min_band, max_band
                do n = min_state, max_state
                    band_contribution(l, b, n) = weight((l - 1) * structure%material%num_bands + b, n)**2
                end do
            end do
        end do
    
    end function band_contribution_band_contribution_range

    function band_contribution_band_contribution_none(energy, weight, structure) result(band_contribution)
        type(type_heterostructure(*)), intent(in) :: structure
        real(dp), intent(in)    :: energy(:)
        complex(dp), intent(in) :: weight(:, :)
        real(dp)                :: band_contribution(structure%num_layers, &
                                                     structure%material%num_bands, &
                                                     size(energy))
        band_contribution = band_contribution_band_contribution(energy, weight, structure, &
                                   1, size(band_contribution(:, 1, 1)), &
                                   1, size(band_contribution(1, :, 1)), &
                                   1, size(band_contribution(1, 1, :)))
                                   
    end function band_contribution_band_contribution_none


end module band_contribution
