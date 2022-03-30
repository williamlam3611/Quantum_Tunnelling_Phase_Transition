module number_density
use constants,         only: pi, dp
use type,              only: type_heterostructure, type_k_mesh
use spectral_function, only: spectral_function_spectral_function
implicit none

private

public :: number_density_number_density

interface number_density_number_density
    module procedure number_density_number_density_none, &
                     number_density_number_density_range
end interface


contains
    function number_density_number_density_range(energy, weight, structure, mesh, energy_range, broadening, &
                                   min_layer, max_layer, min_band, max_band, min_state, max_state) result (number_density)
        type(type_heterostructure(*)), intent(in) :: structure
        type(type_k_mesh),             intent(in) :: mesh
        real(dp), intent(in)    :: energy(:)
        complex(dp), intent(in) :: weight(:, :)
        real(dp), intent(in)    :: energy_range(:)
        real(dp), intent(in)    :: broadening
        integer, intent(in)     :: min_layer
        integer, intent(in)     :: max_layer
        integer, intent(in)     :: min_band
        integer, intent(in)     :: max_band
        integer, intent(in)     :: min_state
        integer, intent(in)     :: max_state
        real(dp)                :: number_density(max_layer - min_layer + 1, &
                                                  max_band - min_band + 1, &
                                                  max_state - min_state + 1)
        number_density = sum(spectral_function_spectral_function(energy, weight, structure, energy_range, broadening, &
                                   min_layer, max_layer, min_band, max_band, min_state, max_state), 4)
        number_density = number_density * (1.0_dp / pi) &
                          * (1.0_dp / mesh%length**2) &
                          * (abs(energy_range(2) - energy_range(1))) &
                          * (1.0_dp / structure%material%crystal_length**2)
    
    end function number_density_number_density_range

    function number_density_number_density_none(energy, weight, structure, mesh, energy_range, broadening) result(number_density)
        type(type_heterostructure(*)), intent(in) :: structure
        type(type_k_mesh),             intent(in) :: mesh
        real(dp), intent(in)    :: energy(:)
        complex(dp), intent(in) :: weight(:, :)
        real(dp), intent(in)    :: energy_range(:)
        real(dp), intent(in)    :: broadening
        real(dp)                :: number_density(structure%num_layers, &
                                                  structure%material%num_bands, &
                                                  size(energy_range))
        number_density = number_density_number_density(energy, weight, structure, mesh, energy_range, broadening, &
                                   1, size(number_density(:, 1, 1)), &
                                   1, size(number_density(1, :, 1)), &
                                   1, size(number_density(1, 1, :)))
                                   
    end function number_density_number_density_none


end module number_density
