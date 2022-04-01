module hqt_schrodinger
use hqt_constants, only: hqt_dp
use hqt_type,      only: hqt_heterostructure
use hqt_matrix,    only: hqt_solve_eigenvalue_equation

implicit none

private

public :: hqt_solve_schrodinger


interface hqt_solve_schrodinger
    module procedure hqt_solve_schrodinger_none, &
                     hqt_solve_schrodinger_integer, &
                     hqt_solve_schrodinger_real
end interface


contains
    subroutine hqt_solve_schrodinger_none(energy, weight, structure, kx, ky)
        type(hqt_heterostructure(*)), intent(in) :: structure
        real(hqt_dp), intent(in)              :: kx
        real(hqt_dp), intent(in)              :: ky
        real(hqt_dp), intent(out)             :: energy(structure%num_layers * structure%material%num_bands)
        complex(hqt_dp), intent(out)          :: weight(structure%num_layers * structure%material%num_bands, &
                                                    structure%num_layers * structure%material%num_bands)
        call hqt_solve_eigenvalue_equation(hqt_supercell(kx, ky, structure), energy, weight)
        
    end subroutine hqt_solve_schrodinger_none
    
    subroutine hqt_solve_schrodinger_integer(energy, weight, num_found, structure, kx, ky, minimum, maximum)
        type(hqt_heterostructure(*)), intent(in) :: structure
        real(hqt_dp), intent(in)              :: kx
        real(hqt_dp), intent(in)              :: ky
        integer, intent(in)               :: minimum
        integer, intent(in)               :: maximum
        integer, intent(out)              :: num_found
        real(hqt_dp), intent(out)             :: energy(structure%num_layers * structure%material%num_bands)
        complex(hqt_dp), intent(out)          :: weight(structure%num_layers * structure%material%num_bands, &
                                                    structure%num_layers * structure%material%num_bands)
        call hqt_solve_eigenvalue_equation(hqt_supercell(kx, ky, structure), energy, weight, num_found, minimum, maximum)
        
    end subroutine hqt_solve_schrodinger_integer
    
    subroutine hqt_solve_schrodinger_real(energy, weight, num_found, structure, kx, ky, minimum, maximum)
        type(hqt_heterostructure(*)), intent(in) :: structure
        real(hqt_dp), intent(in)              :: kx
        real(hqt_dp), intent(in)              :: ky
        real(hqt_dp), intent(in)              :: minimum
        real(hqt_dp), intent(in)              :: maximum
        integer, intent(out)              :: num_found
        real(hqt_dp), intent(out)             :: energy(structure%num_layers * structure%material%num_bands)
        complex(hqt_dp), intent(out)          :: weight(structure%num_layers * structure%material%num_bands, &
                                                    structure%num_layers * structure%material%num_bands)
        call hqt_solve_eigenvalue_equation(hqt_supercell(kx, ky, structure), energy, weight, num_found, minimum, maximum)
        
    end subroutine hqt_solve_schrodinger_real
    
    
    function hqt_supercell(kx, ky, structure) result(supercell)
        type(hqt_heterostructure(*)), intent(in) :: structure
        real(hqt_dp), intent(in)                   :: kx
        real(hqt_dp), intent(in)                   :: ky
        complex(hqt_dp)                            :: zcell(2 * structure%material%max_hopping + 1, &
                                                        structure%material%num_bands, &
                                                        structure%material%num_bands)
        real(hqt_dp)                               :: phase
        integer                                :: x, y, z, i, j, k
        complex(hqt_dp)                            :: supercell(structure%num_layers * structure%material%num_bands, &
                                                            structure%num_layers * structure%material%num_bands)
        ! Build zcell
        zcell = cmplx(0.0_hqt_dp, 0.0_hqt_dp, kind = hqt_dp)
        do x = 1, 2 * structure%material%max_hopping + 1
            do y = 1, 2 * structure%material%max_hopping + 1
                do z = 1, 2 * structure%material%max_hopping + 1
                    do i = 1, structure%material%num_bands
                        do j = 1, structure%material%num_bands
                            phase = -1.0_hqt_dp * (real(x - structure%material%max_hopping - 1, kind = hqt_dp) * kx &
                                             + real(y - structure%material%max_hopping - 1, kind = hqt_dp) * ky)
                            zcell(z, i, j) = zcell(z, i, j) &
                                + structure%material%tunneling(x, y, z, i, j) * dcmplx(cos(phase), sin(phase)) & 
                                    / dble(structure%material%tunneling_weight(x, y, z))
                        end do
                    end do
                end do
            end do
        end do
        ! Build supercell
        supercell = cmplx(0.0_hqt_dp, 0.0_hqt_dp, kind = hqt_dp)
        do i = 1, structure%num_layers
            do j = 1, structure%num_layers
                if (abs(i - j) <= structure%material%max_hopping) then
                    ! Add zcell
                    supercell(structure%material%num_bands * (i - 1) + 1:structure%material%num_bands * i, &
                              structure%material%num_bands * (j - 1) + 1:structure%material%num_bands * j) &
                        = zcell(i - j + structure%material%max_hopping + 1, :, :)
                    ! Add Potential
                    do k = 1, structure%material%num_bands
                        supercell(structure%material%num_bands * (i - 1) + k, &
                                  structure%material%num_bands * (j - 1) + k) &
                            = supercell(structure%material%num_bands * (i - 1) + k, &
                                        structure%material%num_bands * (j - 1) + k) + structure%potential(k)
                    end do
                end if
            end do
        end do
        
    end function hqt_supercell


end module hqt_schrodinger
