module hqt_density
use hqt_constants,         only: hqt_pi, hqt_dp
use hqt_type,              only: hqt_heterostructure,  hqt_k_mesh
use hqt_spectral,          only: hqt_spectral_function
implicit none

private

public :: hqt_number_density

interface hqt_number_density
    module procedure hqt_number_density_none, &
                     hqt_number_density_range
end interface


contains
    function hqt_number_density_range(energy, weight, structure, mesh, energy_range, broadening, &
                                   min_layer, max_layer, min_band, max_band, min_state, max_state) result (number_density)
        type( hqt_heterostructure(*)), intent(in) :: structure
        type( hqt_k_mesh),             intent(in) :: mesh
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
        real(hqt_dp)                :: number_density(max_layer - min_layer + 1, &
                                                  max_band - min_band + 1, &
                                                  max_state - min_state + 1)
        number_density = sum( hqt_spectral_function(energy, weight, structure, energy_range, broadening, &
                                   min_layer, max_layer, min_band, max_band, min_state, max_state), 4)
        number_density = number_density * (1.0_hqt_dp / hqt_pi) &
                          * (1.0_hqt_dp / mesh%length**2) &
                          * (abs(energy_range(2) - energy_range(1))) &
                          * (1.0_hqt_dp / structure%material%crystal_length**2)
    
    end function hqt_number_density_range

    function hqt_number_density_none(energy, weight, structure, mesh, energy_range, broadening) result(number_density)
        type( hqt_heterostructure(*)), intent(in) :: structure
        type( hqt_k_mesh),             intent(in) :: mesh
        real(hqt_dp), intent(in)    :: energy(:)
        complex(hqt_dp), intent(in) :: weight(:, :)
        real(hqt_dp), intent(in)    :: energy_range(:)
        real(hqt_dp), intent(in)    :: broadening
        real(hqt_dp)                :: number_density(structure%num_layers, &
                                                  structure%material%num_bands, &
                                                  size(energy_range))
        number_density = hqt_number_density(energy, weight, structure, mesh, energy_range, broadening, &
                                   1, size(number_density(:, 1, 1)), &
                                   1, size(number_density(1, :, 1)), &
                                   1, size(number_density(1, 1, :)))
                                   
    end function hqt_number_density_none


end module hqt_density
