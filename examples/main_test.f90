program main_test
use constants
use cpu
use hqtlib
implicit none


type(hqtlib_crystal)         :: srtio3
type(hqtlib_heterostructure) :: system
type(hqtlib_k_mesh)          :: mesh


call cpu_start()
! Parameters
srtio3 = hqtlib_deserialise_crystal(cpu_broadcast_file_content("srtio3.dat"))
system = hqtlib_heterostructure(    srtio3, 150)
mesh   = hqtlib_k_mesh(             256,    srtio3, 1.0_dp)

call cpu_stop()


end program main_test
