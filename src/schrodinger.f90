module schrodinger
use constants, only: dp
use type,      only: type_heterostructure
use matrix,    only: matrix_solve_eigenvalue_equation

implicit none

private

public :: schrodinger_solve_schrodinger


interface schrodinger_solve_schrodinger
    module procedure schrodinger_solve_schrodinger_none, &
                     schrodinger_solve_schrodinger_integer, &
                     schrodinger_solve_schrodinger_real
end interface


contains
    subroutine schrodinger_solve_schrodinger_none(energy, weight, structure, kx, ky)
        type(type_heterostructure(*)), intent(in) :: structure
        real(dp), intent(in)              :: kx
        real(dp), intent(in)              :: ky
        real(dp), intent(out)             :: energy(structure%num_layers * structure%material%num_bands)
        complex(dp), intent(out)          :: weight(structure%num_layers * structure%material%num_bands, &
                                                    structure%num_layers * structure%material%num_bands)
        call matrix_solve_eigenvalue_equation(schrodinger_supercell(kx, ky, structure), energy, weight)
        
    end subroutine schrodinger_solve_schrodinger_none
    
    subroutine schrodinger_solve_schrodinger_integer(energy, weight, num_found, structure, kx, ky, minimum, maximum)
        type(type_heterostructure(*)), intent(in) :: structure
        real(dp), intent(in)              :: kx
        real(dp), intent(in)              :: ky
        integer, intent(in)               :: minimum
        integer, intent(in)               :: maximum
        integer, intent(out)              :: num_found
        real(dp), intent(out)             :: energy(structure%num_layers * structure%material%num_bands)
        complex(dp), intent(out)          :: weight(structure%num_layers * structure%material%num_bands, &
                                                    structure%num_layers * structure%material%num_bands)
        call matrix_solve_eigenvalue_equation(schrodinger_supercell(kx, ky, structure), energy, weight, num_found, minimum, maximum)
        
    end subroutine schrodinger_solve_schrodinger_integer
    
    subroutine schrodinger_solve_schrodinger_real(energy, weight, num_found, structure, kx, ky, minimum, maximum)
        type(type_heterostructure(*)), intent(in) :: structure
        real(dp), intent(in)              :: kx
        real(dp), intent(in)              :: ky
        real(dp), intent(in)              :: minimum
        real(dp), intent(in)              :: maximum
        integer, intent(out)              :: num_found
        real(dp), intent(out)             :: energy(structure%num_layers * structure%material%num_bands)
        complex(dp), intent(out)          :: weight(structure%num_layers * structure%material%num_bands, &
                                                    structure%num_layers * structure%material%num_bands)
        call matrix_solve_eigenvalue_equation(schrodinger_supercell(kx, ky, structure), energy, weight, num_found, minimum, maximum)
        
    end subroutine schrodinger_solve_schrodinger_real
    
    
    function schrodinger_supercell(kx, ky, structure) result(supercell)
        type(type_heterostructure(*)), intent(in) :: structure
        real(dp), intent(in)                   :: kx
        real(dp), intent(in)                   :: ky
        complex(dp)                            :: zcell(2 * structure%material%max_hopping + 1, &
                                                        structure%material%num_bands, &
                                                        structure%material%num_bands)
        real(dp)                               :: phase
        integer                                :: x, y, z, i, j, k
        complex(dp)                            :: supercell(structure%num_layers * structure%material%num_bands, &
                                                            structure%num_layers * structure%material%num_bands)
        ! Build zcell
        zcell = cmplx(0.0_dp, 0.0_dp, kind = dp)
        do x = 1, 2 * structure%material%max_hopping + 1
            do y = 1, 2 * structure%material%max_hopping + 1
                do z = 1, 2 * structure%material%max_hopping + 1
                    do i = 1, structure%material%num_bands
                        do j = 1, structure%material%num_bands
                            phase = -1.0_dp * (real(x - structure%material%max_hopping - 1, kind = dp) * kx &
                                             + real(y - structure%material%max_hopping - 1, kind = dp) * ky)
                            zcell(z, i, j) = zcell(z, i, j) &
                                + structure%material%tunneling(x, y, z, i, j) * dcmplx(cos(phase), sin(phase)) & 
                                    / dble(structure%material%tunneling_weight(x, y, z))
                        end do
                    end do
                end do
            end do
        end do
        ! Build supercell
        supercell = cmplx(0.0_dp, 0.0_dp, kind = dp)
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
        
    end function schrodinger_supercell


end module schrodinger
