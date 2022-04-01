program poisson
use hqt
implicit none


type(hqt_crystal)         :: srtio3
type(hqt_heterostructure) :: system
type(hqt_k_mesh)          :: mesh

call hqt_start()

! Parameters
srtio3 = hqt_deserialise_crystal(hqt_broadcast_file_content("./crystals/srtio3.dat"))
system = hqt_heterostructure(    srtio3, 150)
mesh   = hqt_k_mesh(             256,    srtio3, 1.0_hqt_dp)

if (hqt_is_master()) then
    print *, hqt_serialise(mesh)
end if

call hqt_stop()

end program poisson
