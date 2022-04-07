program poisson
use hqt
implicit none


type(hqt_crystal)            :: srtio3
type(hqt_heterostructure)    :: system
type(hqt_k_mesh)             :: mesh

integer, parameter           :: aspect    = 256

integer                      :: path(aspect)
integer                      :: path_nodes(3)
real(hqt_dp),    allocatable :: energy(:)
complex(hqt_dp), allocatable :: weight(:, :)
real(hqt_dp)                 :: conducting_band_minimum
integer                      :: num_found

integer                      :: i


complex(hqt_dp)              :: test(2,2)

call hqt_start()

! Parameters
srtio3 = hqt_deserialise_crystal(hqt_broadcast_file_content("./srtio3.dat"))
system = hqt_heterostructure(    srtio3, 2)
mesh   = hqt_k_mesh(             aspect, srtio3, 1.0_hqt_dp)

call hqt_path_map(mesh, hqt_pi * reshape((/ 1.0_hqt_dp, 1.0_hqt_dp, &                 ! M
                                            0.0_hqt_dp, 0.0_hqt_dp, &                 ! Gamma
                                            1.0_hqt_dp, 0.0_hqt_dp /), (/ 2, 3 /)), & ! X
                                            size(path), path, path_nodes)

allocate(energy(srtio3%num_bands * system%num_layers))
allocate(weight(srtio3%num_bands * system%num_layers, srtio3%num_bands * system%num_layers))                                            

call hqt_solve_eigenvalue_equation(hqt_supercell(0.0_hqt_dp, 0.0_hqt_dp, system), energy, weight)
print *, energy
!print *, real(hqt_supercell(0.0_hqt_dp, 0.0_hqt_dp, system))

!call hqt_solve_schrodinger(energy, weight, num_found, system, &
!                           srtio3%conducting_band_minimum_kx, srtio3%conducting_band_minimum_ky, 1, 1)
!conducting_band_minimum = energy(1)
!print *, energy

call hqt_stop()

end program poisson
