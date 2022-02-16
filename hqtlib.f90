module hqtlib
use cpu
use transform
use bulk
use matrix
use spectral
use route
implicit none

private

public :: hqtlib_find_energy_and_weight, hqtlib_find_max_num_states, &
          hqtlib_find_energy_max, hqtlib_find_energy_min, &
          hqtlib_find_density, hqtlib_find_spectral_distribution

interface hqtlib_find_energy_and_weight
    module procedure hqtlib_find_energy_and_weight_integer_pot, hqtlib_find_energy_and_weight_real_pot, &
                     hqtlib_find_energy_and_weight_integer_layer, hqtlib_find_energy_and_weight_real_layer
end interface hqtlib_find_energy_and_weight

interface hqtlib_find_energy_max
    module procedure hqtlib_find_energy_max_layer, hqtlib_find_energy_max_pot
end interface hqtlib_find_energy_max

interface hqtlib_find_energy_min
    module procedure hqtlib_find_energy_min_layer, hqtlib_find_energy_min_pot
end interface hqtlib_find_energy_min

contains
    function hqtlib_find_energy_min_pot(hr, hrw, pot, broadening) result(energy_min)
        complex*16, intent(in)       :: hr(:, :, :, :, :)
        integer, intent(in)          :: hrw(:, :, :)
        real*8, intent(in)           :: pot(:)
        real*8, optional, intent(in) :: broadening
        real*8                       :: energy(size(pot) * size(hr(1, 1, 1, :, 1)))
        complex*16                   :: weight(size(pot) * size(hr(1, 1, 1, :, 1)), &
                                               size(pot) * size(hr(1, 1, 1, :, 1)))
        integer                      :: num_found
        real*8                       :: energy_min
        call hqtlib_find_energy_and_weight(energy, weight, num_found, 0d0, 0d0, hr, hrw, pot, &
                                           1, 1)
        energy_min = energy(1) 
        if (present(broadening)) energy_min = energy_min &
                                              - sqrt((broadening / 10d0) - broadening**2) - sqrt(broadening - broadening**2)
        
    end function hqtlib_find_energy_min_pot
    
    function hqtlib_find_energy_min_layer(hr, hrw, num_layers, broadening) result(energy_min)
        complex*16, intent(in)       :: hr(:, :, :, :, :)
        integer, intent(in)          :: hrw(:, :, :)
        integer, intent(in)          :: num_layers
        real*8, optional, intent(in) :: broadening
        real*8                       :: pot(num_layers)
        real*8                       :: energy_min
        pot = 0d0
        if (present(broadening)) then
            energy_min = hqtlib_find_energy_min(hr, hrw, pot, broadening)
        else
            energy_min = hqtlib_find_energy_min(hr, hrw, pot)
        end if
        
    end function hqtlib_find_energy_min_layer
    
    function hqtlib_find_energy_max_pot(hr, hrw, pot, length_scale, broadening) result(energy_max)
        complex*16, intent(in)       :: hr(:, :, :, :, :)
        integer, intent(in)          :: hrw(:, :, :)
        real*8, intent(in)           :: pot(:)
        real*8, intent(in)           :: length_scale
        real*8, optional, intent(in) :: broadening
        real*8                       :: energy(size(pot) * size(hr(1, 1, 1, :, 1)))
        complex*16                   :: weight(size(pot) * size(hr(1, 1, 1, :, 1)), &
                                               size(pot) * size(hr(1, 1, 1, :, 1)))
        integer                      :: num_found
        real*8                       :: energy_max
        call hqtlib_find_energy_and_weight(energy, weight, num_found, &
                                           3.141592653589793d0 * length_scale, 3.141592653589793d0 * length_scale, hr, hrw, pot, &
                                           size(pot) * size(hr(1, 1, 1, :, 1)), size(pot) * size(hr(1, 1, 1, :, 1)))
        energy_max = energy(1)
        if (present(broadening)) energy_max = energy_max &
                                              + sqrt((broadening / 10d0) - broadening**2) + sqrt(broadening - broadening**2)
        
    end function hqtlib_find_energy_max_pot
    
    function hqtlib_find_energy_max_layer(hr, hrw, num_layers, length_scale, broadening) result(energy_max)
        complex*16, intent(in)       :: hr(:, :, :, :, :)
        integer, intent(in)          :: hrw(:, :, :)
        integer, intent(in)          :: num_layers
        real*8, intent(in)           :: length_scale
        real*8, optional, intent(in) :: broadening
        real*8                       :: pot(num_layers)
        real*8                       :: energy_max
        pot = 0d0
        if (present(broadening)) then
            energy_max = hqtlib_find_energy_max(hr, hrw, pot, length_scale, broadening)
        else
            energy_max = hqtlib_find_energy_max(hr, hrw, pot, length_scale)
        end if
        
    end function hqtlib_find_energy_max_layer
    
    function hqtlib_find_max_num_states(hr, hrw, pot, energy_min, energy_max) result(max_num_states)
        complex*16, intent(in) :: hr(:, :, :, :, :)
        integer, intent(in)    :: hrw(:, :, :)
        real*8, intent(in)     :: pot(:)
        real*8, intent(in)     :: energy_min, energy_max
        real*8                 :: energy(size(pot) * size(hr(1, 1, 1, :, 1)))
        complex*16             :: weight(size(pot) * size(hr(1, 1, 1, :, 1)), &
                                         size(pot) * size(hr(1, 1, 1, :, 1)))
        integer                :: max_num_states
        call hqtlib_find_energy_and_weight(energy, weight, max_num_states, 0d0, 0d0, hr, hrw, pot, &
                                           energy_min, energy_max)
        
    end function hqtlib_find_max_num_states

    subroutine hqtlib_find_energy_and_weight_integer_pot(energy, weight, num_found, kx, ky, hr, hrw, pot, &
                                                   min_energy, max_energy)
        complex*16, intent(in)  :: hr(:, :, :, :, :)
        integer, intent(in)     :: hrw(:, :, :)
        real*8, intent(in)      :: kx, ky
        real*8, intent(in)      :: pot(:)
        integer, intent(in)     :: min_energy, max_energy
        complex*16              :: hb(size(pot) * size(hr(1, 1, 1, :, 1)), size(pot) * size(hr(1, 1, 1, :, 1)))
        real*8, intent(out)     :: energy(size(pot) * size(hr(1, 1, 1, :, 1)))
        complex*16, intent(out) :: weight(size(pot) * size(hr(1, 1, 1, :, 1)), size(pot) * size(hr(1, 1, 1, :, 1)))
        integer, intent(out)    :: num_found
        hb = bulk_build(transform_r_to_kz(hr, hrw, kx, ky), size(pot))
        call bulk_add_potential(hb, pot)
        call matrix_get_eigen(hb, energy, weight, num_found, min_energy, max_energy)
        
    end subroutine hqtlib_find_energy_and_weight_integer_pot
    
    subroutine hqtlib_find_energy_and_weight_integer_layer(energy, weight, num_found, kx, ky, hr, hrw, num_layers, &
                                                   min_energy, max_energy)
        complex*16, intent(in)  :: hr(:, :, :, :, :)
        integer, intent(in)     :: hrw(:, :, :)
        real*8, intent(in)      :: kx, ky
        integer, intent(in)     :: num_layers
        integer, intent(in)     :: min_energy, max_energy
        real*8                  :: pot(num_layers)
        real*8, intent(out)     :: energy(num_layers * size(hr(1, 1, 1, :, 1)))
        complex*16, intent(out) :: weight(num_layers * size(hr(1, 1, 1, :, 1)), num_layers * size(hr(1, 1, 1, :, 1)))
        integer, intent(out)    :: num_found
        pot = 0d0
        call hqtlib_find_energy_and_weight(energy, weight, num_found, kx, ky, hr, hrw, pot, min_energy, max_energy)
        
    end subroutine hqtlib_find_energy_and_weight_integer_layer
    
    subroutine hqtlib_find_energy_and_weight_real_pot(energy, weight, num_found, kx, ky, hr, hrw, pot, &
                                                  min_energy, max_energy)
        complex*16, intent(in)  :: hr(:, :, :, :, :)
        integer, intent(in)     :: hrw(:, :, :)
        real*8, intent(in)      :: kx, ky
        real*8, intent(in)      :: pot(:)
        real*8, intent(in)      :: min_energy, max_energy
        complex*16              :: hb(size(pot) * size(hr(1, 1, 1, :, 1)), size(pot) * size(hr(1, 1, 1, :, 1)))
        real*8, intent(out)     :: energy(size(pot) * size(hr(1, 1, 1, :, 1)))
        complex*16, intent(out) :: weight(size(pot) * size(hr(1, 1, 1, :, 1)), size(pot) * size(hr(1, 1, 1, :, 1)))
        integer, intent(out)    :: num_found
        hb = bulk_build(transform_r_to_kz(hr, hrw, kx, ky), size(pot))
        call bulk_add_potential(hb, pot)
        call matrix_get_eigen(hb, energy, weight, num_found, min_energy, max_energy)
        
    end subroutine hqtlib_find_energy_and_weight_real_pot
    
    subroutine hqtlib_find_energy_and_weight_real_layer(energy, weight, num_found, kx, ky, hr, hrw, num_layers, &
                                                   min_energy, max_energy)
        complex*16, intent(in)  :: hr(:, :, :, :, :)
        integer, intent(in)     :: hrw(:, :, :)
        real*8, intent(in)      :: kx, ky
        integer, intent(in)     :: num_layers
        real*8, intent(in)      :: min_energy, max_energy
        real*8                  :: pot(num_layers)
        real*8, intent(out)     :: energy(num_layers * size(hr(1, 1, 1, :, 1)))
        complex*16, intent(out) :: weight(num_layers * size(hr(1, 1, 1, :, 1)), num_layers * size(hr(1, 1, 1, :, 1)))
        integer, intent(out)    :: num_found
        pot = 0d0
        call hqtlib_find_energy_and_weight(energy, weight, num_found, kx, ky, hr, hrw, pot, min_energy, max_energy)
        
    end subroutine hqtlib_find_energy_and_weight_real_layer

    subroutine hqtlib_find_density(density, energy, weight, energy_range, broadening, total_num_bands, crystal_length, num_k, &
                                   layers, bands, energy_levels)
        real*8, intent(in)               :: energy(:)
        complex*16, intent(in)           :: weight(:, :)
        real*8, intent(in)               :: energy_range(:)
        real*8, intent(in)               :: broadening
        integer, intent(in)              :: total_num_bands
        real*8, intent(in)               :: crystal_length
        integer, intent(in)              :: num_k
        integer, optional, intent(in)    :: layers(:)
        integer, optional, intent(in)    :: bands(:)
        integer, optional, intent(in)    :: energy_levels(:)
        integer                          :: total_num_layers
        integer, allocatable             :: temp_layers(:)
        integer, allocatable             :: temp_bands(:)
        integer, allocatable             :: temp_energy_levels(:)
        real*8, allocatable              :: spect(:, :, :, :)
        real*8, allocatable, intent(inout) :: density(:, :, :)
        total_num_layers = size(weight(:, 1)) / total_num_bands
        if (.not. (present(layers) .or. present(bands) .or. present(energy_levels))) then
            if (allocated(density)) then
                if (any(shape(density) .ne. (/ total_num_layers, total_num_bands, size(energy) /))) then
                    deallocate(density)
                    allocate(density(total_num_layers, total_num_bands, size(energy)))
                end if
            else
                allocate(density(total_num_layers, total_num_bands, size(energy)))
            end if
            density = 0d0
            call hqtlib_find_spectral_distribution(spect, energy, weight, energy_range, broadening, total_num_bands)
        else
            if (present(layers)) then
                temp_layers = layers
            else
                temp_layers = route_range_integer(1, total_num_layers, total_num_layers)
            end if
            if (present(bands)) then
                temp_bands = bands
            else
                temp_bands = route_range_integer(1, total_num_bands, total_num_bands)
            end if
            if (present(energy_levels)) then
                temp_energy_levels = energy_levels
            else
                temp_energy_levels = route_range_integer(1, size(energy), size(energy))
            end if
            if (allocated(density)) then
                if (any(shape(density) .ne. (/ size(temp_layers), size(temp_bands), size(temp_energy_levels) /))) then
                    deallocate(density)
                    allocate(density(size(temp_layers), size(temp_bands), size(temp_energy_levels)))
                end if
            else
                allocate(density(size(temp_layers), size(temp_bands), size(temp_energy_levels)))
            end if
            density = 0d0
            call hqtlib_find_spectral_distribution(spect, energy, weight, energy_range, broadening, total_num_bands, &
                                                   temp_layers, temp_bands, temp_energy_levels)
        end if
        density = sum(spect, 4)
        density = density * (1d0 / 3.141592653589793d0) &
                          * (1d0 / dble(num_k)) &
                          * ((energy_range(2) - energy_range(1)) / size(energy_range)) &
                          * (1d0 / crystal_length**2)
    
    end subroutine hqtlib_find_density
    
    subroutine hqtlib_find_spectral_distribution(spect, energy, weight, energy_range, broadening, total_num_bands, &
                                                 layers, bands, energy_levels)
        real*8, intent(in)               :: energy(:)
        complex*16, intent(in)           :: weight(:, :)
        real*8, intent(in)               :: energy_range(:)
        real*8, intent(in)               :: broadening
        integer, intent(in)              :: total_num_bands
        integer, optional, intent(in)    :: layers(:)
        integer, optional, intent(in)    :: bands(:)
        integer, optional, intent(in)    :: energy_levels(:)
        integer                          :: total_num_layers
        integer                          :: l, b, n, e
        integer, allocatable             :: temp_layers(:)
        integer, allocatable             :: temp_bands(:)
        integer, allocatable             :: temp_energy_levels(:)
        real*8, allocatable, intent(inout) :: spect(:, :, :, :)
        total_num_layers = size(weight(:, 1)) / total_num_bands
        if (.not. (present(layers) .or. present(bands) .or. present(energy_levels))) then
            if (allocated(spect)) then
                if (any(shape(spect) .ne. (/ total_num_layers, total_num_bands, size(energy), size(energy_range) /))) then
                    deallocate(spect)
                    allocate(spect(total_num_layers, total_num_bands, size(energy), size(energy_range)))
                end if
            else
                allocate(spect(total_num_layers, total_num_bands, size(energy), size(energy_range)))
            end if
            do l = 1, total_num_layers
                do b = 1, total_num_bands
                    do n = 1, size(energy)
                        do e = 1, size(energy_range)
                            spect(l, b, n, e) = spectral_function(energy(n), dble(abs(weight((l - 1) * total_num_bands + b, &
                                                                                             n)**2)), energy_range(e), broadening)
                        end do
                    end do
                end do
            end do
        else
            if (present(layers)) then
                temp_layers = layers
            else
                temp_layers = route_range_integer(1, total_num_layers, total_num_layers)
            end if
            if (present(bands)) then
                temp_bands = bands
            else
                temp_bands = route_range_integer(1, total_num_bands, total_num_bands)
            end if
            if (present(energy_levels)) then
                temp_energy_levels = energy_levels
            else
                temp_energy_levels = route_range_integer(1, size(energy), size(energy))
            end if
            if (allocated(spect)) then
                if (any(shape(spect) .ne. (/ size(layers), size(bands), size(energy_levels), size(energy_range) /))) then
                    deallocate(spect)
                    allocate(spect(size(layers), size(bands), size(energy_levels), size(energy_range)))
                end if
            else
                allocate(spect(size(layers), size(bands), size(energy_levels), size(energy_range)))
            end if
            do l = 1, size(temp_layers)
                do b = 1, size(temp_bands)
                    do n = 1, size(temp_energy_levels)
                        do e = 1, size(energy_range)
                            spect(l, b, n, e) = spectral_function(energy(temp_energy_levels(n)), &
                                                         dble(abs(weight((temp_layers(l) - 1) * total_num_bands + temp_bands(b), &
                                                                  temp_energy_levels(n))**2)), &
                                                                  energy_range(e), broadening)
                        end do
                    end do
                end do
            end do
        end if
    
    end subroutine hqtlib_find_spectral_distribution

end module hqtlib
